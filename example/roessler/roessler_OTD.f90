module Roessler_OTD
   ! Standard Library.
   use stdlib_stats_distribution_normal, only: normal => rvs_normal
   use stdlib_optval, only: optval
   use stdlib_sorting, only: sort
   use stdlib_linalg, only: eig, svdvals, eye
   ! RKLIB module for time integration.
   use rklib_module
   ! LightKrylov for linear algebra.
   use LightKrylov
   use LightKrylov_Constants
   use LightKrylov_AbstractVectors
   ! Roessler
   use Roessler
   implicit none

   character(len=*), parameter, private :: this_module = 'Roessler_OTD'

   public :: a, b, c

   !------------------------------
   !-----     PARAMETERS     -----
   !------------------------------

   integer, parameter :: r = 2
   real(dp), parameter :: t_GS = 5.0_dp ! In finite-precision arithmetic we need to reorthonormalize sometimes
   character(len=*), parameter :: report_file_OTD = 'example/roessler/output_roessler_OTD.txt'
   character(len=*), parameter :: report_file_OTD_LE = 'example/roessler/output_roessler_OTD_LE.txt'

   ! Reference values (https://chaosbook.org/extras/simon/Roessler.html, orbit 1)
   real(dp), dimension(2), parameter :: EV_ref = (/0.097000856_dp, 0.097000856_dp/)
   real(dp), dimension(2), parameter :: LE_ref = (/0.0_dp, 0.149141556_dp/)

   !-------------------------------------------
   !-----     LIGHTKRYLOV VECTOR TYPE     -----
   !-------------------------------------------

   type, extends(abstract_vector_rdp), public :: pos_vector
      real(dp) :: x = 0.0_dp
      real(dp) :: y = 0.0_dp
      real(dp) :: z = 0.0_dp
   contains
      private
      procedure, pass(self), public :: zero => zero_p
      procedure, pass(self), public :: dot => dot_p
      procedure, pass(self), public :: scal => scal_p
      procedure, pass(self), public :: axpby => axpby_p
      procedure, pass(self), public :: rand => rand_p
      procedure, pass(self), public :: get_size => get_size_p
   end type pos_vector

contains

   !=========================================================
   !=========================================================
   !=====                                               =====
   !=====     LIGHTKRYLOV MANDATORY IMPLEMENTATIONS     =====
   !=====                                               =====
   !=========================================================
   !=========================================================

   !----------------------------------------------------
   !-----     TYPE-BOUND PROCEDURE FOR VECTORS     -----
   !----------------------------------------------------

   subroutine zero_p(self)
      class(pos_vector), intent(inout) :: self
      ! spatial coordinates of initial condition for orbit
      self%x = 0.0_dp
      self%y = 0.0_dp
      self%z = 0.0_dp
   end subroutine zero_p

   real(dp) function dot_p(self, vec) result(alpha)
      class(pos_vector), intent(in) :: self
      class(abstract_vector_rdp), intent(in) :: vec
      select type (vec)
      type is (pos_vector)
         alpha = self%x*vec%x + self%y*vec%y + self%z*vec%z
      class default
         call type_error('vec','pos_vector','IN',this_module,'dot_p')
      end select
   end function dot_p

   subroutine scal_p(self, alpha)
      class(pos_vector), intent(inout) :: self
      real(dp), intent(in)    :: alpha
      self%x = self%x*alpha
      self%y = self%y*alpha
      self%z = self%z*alpha
   end subroutine scal_p

   subroutine axpby_p(alpha, vec, beta, self)
      class(pos_vector), intent(inout) :: self
      class(abstract_vector_rdp), intent(in)    :: vec
      real(dp), intent(in)    :: alpha, beta
      select type (vec)
      type is (pos_vector)
         self%x = beta*self%x + alpha*vec%x
         self%y = beta*self%y + alpha*vec%y
         self%z = beta*self%z + alpha*vec%z
      class default
         call type_error('vec','pos_vector','IN',this_module,'axpby_p')
      end select
   end subroutine axpby_p

   integer function get_size_p(self) result(N)
      class(pos_vector), intent(in) :: self
      N = npts + 1
   end function get_size_p

   subroutine rand_p(self, ifnorm)
      class(pos_vector), intent(inout) :: self
      logical, optional, intent(in)    :: ifnorm
      logical :: normalized
      real(dp) :: mu, var
      real(dp) :: alpha

      mu = 0.0_dp
      var = 1.0_dp
      self%x = normal(mu, var)
      self%y = normal(mu, var)
      self%z = normal(mu, var)

      normalized = optval(ifnorm, .false.)
      if (normalized) then
         alpha = self%norm()
         call self%scal(1.0_dp/alpha)
      end if
   end subroutine rand_p

   subroutine OTD_rhs(me, t, x, f)
      ! Time-integrator.
      class(rk_class), intent(inout)             :: me
      ! Current time.
      real(dp), intent(in)                :: t
      ! State vector.
      real(dp), dimension(:), intent(in)  :: x
      ! Time-derivative.
      real(dp), dimension(:), intent(out) :: f
      ! internal
      real(dp), dimension(npts)   :: bf
      real(dp), dimension(npts, r) :: q, Lq, fp
      real(dp), dimension(r, r)    :: Lr, Phi
      integer :: i, j

      bf = x(:npts)
      q = reshape(x(npts + 1:(r + 1)*npts), shape(q))

      call nonlinear_roessler(bf, f(:npts))
      do i = 1, r
         call linear_roessler(q(:, i), bf, Lq(:, i))
      end do
      ! build reduced operator and rotation matrix
      Lr = 0.0_dp
      Phi = 0.0_dp
      do i = 1, r
         do j = 1, r
            Lr(i, j) = dot_product(q(:, i), Lq(:, j))
         end do
         do j = i + 1, r
            Phi(i, j) = Lr(i, j)
            Phi(j, i) = -Phi(i, j)
         end do
      end do
      ! Construct forcing and add the rhs
      fp = Lq - matmul(q, Lr - Phi)
      f(npts + 1:(r + 1)*npts) = reshape(fp, [size(fp)])
      ! Construct dFTLE(i)/dt = Lr(i,i)
      do i = 1, r
         f((r + 1)*npts + i) = Lr(i, i)
      end do
   end subroutine OTD_rhs

   !------------------------------
   !-----     INTEGRATOR     -----
   !------------------------------

   subroutine OTD_step(integrator, bf, vec_in, FTLE_in, time, Tstep, vec_out, FTLE_out)
      ! Integrator
      class(rk_class), intent(inout)  :: integrator
      ! Basic state
      class(abstract_vector_rdp), intent(inout)  :: bf
      ! Input vector.
      class(abstract_vector_rdp), intent(in)     :: vec_in(r)
      ! Fundamental solution matrix
      real(dp), intent(in)     :: FTLE_in(r)
      ! Current simulation time
      real(dp), intent(inout)  :: time
      ! Integration time for this step.
      real(dp), intent(in)     :: Tstep
      ! Output vector.
      class(abstract_vector_rdp), intent(out)    :: vec_out(r)
      ! Fundamental solution matrix
      real(dp), intent(out)    :: FTLE_out(r)

      ! internals
      real(dp)                          :: dt = 1.0_dp
      real(dp), dimension((r + 1)*npts + r) :: pos_in, pos_out
      integer                           :: i

      select type (integrator)
      type is (rks54_class)
         ! Get the state.
         call get_pos(bf, pos_in(:npts)) ! bf is a state vector
         do i = 1, r
            call get_pos(vec_in(i), pos_in(npts*i + 1:npts*(i + 1)))
         end do
         pos_in(npts*(r + 1) + 1:) = FTLE_in
         ! Integrate.
         call integrator%integrate(time, pos_in, dt, time + Tstep, pos_out)
         ! Pass-back the state.
         call set_pos(pos_out(:npts), bf)
         do i = 1, r
            call set_pos(pos_out(npts*i + 1:npts*(i + 1)), vec_out(i))
         end do
         FTLE_out = pos_out(npts*(r + 1) + 1:)
         time = time + Tstep
      class default
         call type_error('integrator','rks54_class','INOUT',this_module,'OTD_step')
      end select
   end subroutine OTD_step

   subroutine OTD_map(bf, vec_in, Tend, vec_out, t_FTLE, if_rst)
      ! Basic state
      class(abstract_vector_rdp), intent(inout) :: bf
      ! Input vector.
      class(abstract_vector_rdp), intent(inout) :: vec_in(r)
      ! Integration time.
      real(dp), intent(in)    :: Tend
      ! Output vector.
      class(abstract_vector_rdp), intent(out)   :: vec_out(r)
      ! Integration time for FLTEs
      real(dp), intent(in)    :: t_FTLE
      ! restart trajectory at t_FLTE?
      logical, optional, intent(in)    :: if_rst

      ! internals
      type(rks54_class)         :: OTD_roessler
      integer                   :: idx(1)
      real(dp)                  :: time, t_complete, t1, t2, tvec(2)
      ! logical                   :: if_GS
      integer                   :: i, j, p_cnt
      real(dp), dimension(r)    :: FTLE_in, FTLE_out
      real(dp), dimension(npts) :: bf_bkp

      integer, parameter :: ndof = npts*(r + 1) + r

      ! Initialize integrator.
      call OTD_roessler%initialize(n=ndof, f=OTD_rhs, report=OTD_report_file, report_rate=20)

      ! Save IC
      call get_pos(bf, bf_bkp)

      ! Initialization
      time = 0.0_dp
      FTLE_in = 0.0_dp
      p_cnt = 0
      tvec = (/t_GS, t_FTLE/)
      idx = minloc(tvec)
      t1 = minval(tvec)
      t2 = maxval(tvec)
      ! Regular integration with much less frequent reorthonormalization
      do j = 1, floor(Tend/t2)
         do i = 1, floor(t2/t1)
            call OTD_step(OTD_roessler, bf, vec_in, FTLE_in, time, t1, vec_out, FTLE_out)
            if (t1 == t_GS) then
               call orthonormalize_basis(vec_out)                 ! Reorthonormalize
               FTLE_in = FTLE_out
            else
               p_cnt = p_cnt + 1
               call report_LE(FTLE_out, time, t_FTLE, p_cnt)
               FTLE_in = 0.0_dp                                   ! reset FTLE computation
               if (if_rst) call set_pos(bf_bkp, bf)               ! reset orbit to avoid orbit drift
            end if
            call copy(vec_in, vec_out)
         end do
         t_complete = modulo(t2, t1)
         if (t_complete > 1e-4_dp) call OTD_step(OTD_roessler, bf, vec_in, FTLE_in, time, t_complete, vec_out, FTLE_out)
         if (t2 == t_GS) then
            call orthonormalize_basis(vec_out)                    ! Reorthonormalize
            FTLE_in = FTLE_out
         else
            p_cnt = p_cnt + 1
            call report_LE(FTLE_out, time, t_FTLE, p_cnt)
            FTLE_in = 0.0_dp                                      ! reset FTLE computation
            if (if_rst) call set_pos(bf_bkp, bf)                  ! reset orbit to avoid orbit drift
         end if
         call copy(vec_in, vec_out)
      end do
      t_complete = modulo(Tend, t2)
      if (t_complete > 1e-4_dp) call OTD_step(OTD_roessler, bf, vec_in, FTLE_in, time, t_complete, vec_out, FTLE_out)
   end subroutine OTD_map

   !-------------------------------------------
   !-----     MISCELLANEOUS UTILITIES     -----
   !-------------------------------------------

   subroutine get_pos(vec_in, pos)
      class(abstract_vector_rdp), intent(in)  :: vec_in
      real(dp), dimension(npts), intent(out) :: pos
      pos = 0.0_dp
      select type (vec_in)
      type is (pos_vector)
         pos(1) = vec_in%x
         pos(2) = vec_in%y
         pos(3) = vec_in%z
      class default
         call type_error('vec','pos_vector','IN',this_module,'get_pos')
      end select
   end subroutine get_pos

   subroutine set_pos(pos, vec_out)
      real(dp), dimension(npts), intent(in)  :: pos
      class(abstract_vector_rdp), intent(out) :: vec_out
      select type (vec_out)
      type is (pos_vector)
         vec_out%x = pos(1)
         vec_out%y = pos(2)
         vec_out%z = pos(3)
      class default
         call type_error('vec_out','pos_vector','OUT',this_module,'set_pos')
      end select
   end subroutine set_pos

   subroutine OTD_report_file(me, t, x)
      class(rk_class), intent(inout)      :: me
      real(dp), intent(in)                :: t
      real(dp), dimension(:), intent(in)  :: x

      ! internal
      real(dp), dimension(npts)      :: bf
      real(dp), dimension(npts, r)    :: q, Lq
      real(dp), dimension(r, r)       :: Lr, Lsym, qTq
      real(dp), dimension(r)         :: FTLE, s
      complex(dp), dimension(r)      :: l
      complex(dp), dimension(r, r)    :: v
      complex(dp), dimension(npts, r) :: u, su
      integer                        :: i, j, iunit
      logical                        :: is_cc
      integer, allocatable           :: idx(:)

      bf = x(:npts)
      q = reshape(x(npts + 1:(r + 1)*npts), shape(q))
      FTLE = x(npts*(r + 1) + 1:)/t

      do i = 1, r
         call linear_roessler(q(:, i), bf, Lq(:, i))
      end do
      ! build reduced operator
      Lr = 0.0_dp
      qTq = 0.0_dp
      do i = 1, r
         do j = 1, r
            Lr(i, j) = dot_product(Lq(:, i), q(:, j))
            qTq(i, j) = dot_product(q(:, i), q(:, j))
         end do
      end do
      ! spectral analysis
      Lsym = 0.5_dp*(Lr + transpose(Lr))
      call eig(Lsym, l, right=v)
      s = real(l)
      idx = maxloc(s)
      i = idx(1)
      su = matmul(q, v)
      call eig(Lr, l, right=v)
      u = matmul(q, v)
      is_cc = .false.
      if (abs(aimag(l(1))) > 0.0_dp) is_cc = .true.

      open (newunit=iunit, file=report_file_OTD, status='old', action='write', position='append')
      write (iunit, '(*(E16.9,1X))', ADVANCE='NO') t, bf, q
      if (is_cc) then
         write (iunit, '(I2,1X,*(E16.9,1X))', ADVANCE='NO') 1, s(i), real(su(:, i)), &
            real(l(1)), aimag(l(1)), real(u(:, 1)), aimag(u(:, 1))
      else
         write (iunit, '(I2,1X,*(E16.9,1X))', ADVANCE='NO') 0, s(i), real(su(:, i)), &
            real(l), real(u)
      end if
      write (iunit, '(*(E16.9,1X))') FTLE, (qTq(i, i) - 1, i=1, r), &
         ((qTq(i, j), j=1, i - 1), i=2, r)
      close (iunit)
   end subroutine OTD_report_file

   subroutine report_LE(FTLE, time, period, p_cnt)
      ! Integrated FTLE values
      real(dp), dimension(r), intent(in) :: FTLE
      ! simulation
      real(dp), intent(in) :: time
      ! period
      real(dp), intent(in) :: period
      ! period counter
      integer, intent(in) :: p_cnt
      ! internal
      integer :: iunit
      real(dp), dimension(r) :: LE
      LE = FTLE/period
      call sort(LE)
      print '(A10,I3,A3,1X,*(F16.9,1X))', 'OTD:  t=', p_cnt, 'T ', LE
      open (newunit=iunit, file=report_file_OTD_LE, status='old', action='write', position='append')
      write (iunit, '(F16.9,1X,I16,1X,*(F16.9,1X))') time, p_cnt, LE, LE_ref
      close (iunit)
   end subroutine report_LE

   subroutine write_header()
      ! internals
      integer :: i, j, iunit
      open (newunit=iunit, file=report_file_OTD, status='new', action='write')
      ! time, baseflow
      write (iunit, '(*(A16,1X))', ADVANCE='NO') 't', 'BF_x', 'BF_y', 'BF_z'
      ! basis vectors
      do i = 1, r
         write (iunit, '(*(A13,I1,A2,1X))', ADVANCE='NO') 'q', i, '_x', 'q', i, '_y', 'q', i, '_z'
      end do
      ! instantaneous numerical abscissa of the reduced operator
      write (iunit, '(A2,1X,A16,1X)', ADVANCE='NO') 'cc', 'sigma_1'
      ! instantaneous direction of largest possible growth
      write (iunit, '(*(A13,I1,A2,1X))', ADVANCE='NO') 'us', i, '_x', 'us', i, '_y', 'us', i, '_z'
      ! instantaneous eigenvalues of the reduced operator
      do i = 1, r
         write (iunit, '(A15,I1,1X)', ADVANCE='NO') 'l_', i
      end do
      ! instantaneous projected eigenvectors of the reduced operator
      do i = 1, r
         write (iunit, '(*(A13,I1,A2,1X))', ADVANCE='NO') 'u', i, '_x', 'u', i, '_y', 'u', i, '_z'
      end do
      ! instantaneous FLTEs
      do i = 1, r
         write (iunit, '(A15,I1,1X)', ADVANCE='NO') 'FTLE_', i
      end do
      ! orthonormality of the basis vectors
      do i = 1, r
         write (iunit, '(A7,I1,A2,I1,A5,1X)', ADVANCE='NO') '<q', i, ',q', i, '> - 1'
      end do
      do i = 2, r
         do j = 1, i - 1
            write (iunit, '(A11,I1,A2,I1,A1,1X)', ADVANCE='NO') '<q', j, ',q', i, '>'
         end do
      end do
      write (iunit, *) ''; close (iunit)
   end subroutine write_header

   subroutine write_header_LE()
      ! internals
      integer :: i, iunit
      open (newunit=iunit, file=report_file_OTD_LE, status='new', action='write')
      ! time, baseflow
      write (iunit, '(*(A16,1X))', ADVANCE='NO') 't', 'period'
      ! LE
      do i = 1, r
         write (iunit, '(A15,I1,1X)', ADVANCE='NO') 'LE_', i
      end do
      ! LE_ref
      do i = 1, r
         write (iunit, '(A15,I1,1X)', ADVANCE='NO') 'LEref_', i
      end do
      write (iunit, *) ''; close (iunit)
   end subroutine write_header_LE

end module Roessler_OTD
