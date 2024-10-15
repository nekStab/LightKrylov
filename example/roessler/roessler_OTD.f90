module Roessler_OTD
   ! Standard Library.
   use stdlib_stats_distribution_normal, only: normal => rvs_normal
   use stdlib_optval, only: optval
   use stdlib_linalg, only: eigvals, svdvals, eye
   ! RKLIB module for time integration.
   use rklib_module
   ! LightKrylov for linear algebra.
   use LightKrylov
   use LightKrylov, only: wp => dp
   use LightKrylov_AbstractVectors
   ! Roessler
   use Roessler
   implicit none
 
   character*128, parameter, private :: this_module = 'Roessler_OTD'
 
   public :: a, b, c
 
   !------------------------------
   !-----     PARAMETERS     -----
   !------------------------------
 
   integer, parameter :: r = 2
   character(len=128), parameter :: file = 'example/roessler/output_roessler_OTD.txt'

   !-------------------------------------------
   !-----     LIGHTKRYLOV VECTOR TYPE     -----
   !-------------------------------------------
 
   type, extends(abstract_vector_rdp), public :: pos_vector
      real(wp) :: x = 0.0_wp
      real(wp) :: y = 0.0_wp
      real(wp) :: z = 0.0_wp
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
      self%x = 0.0_wp
      self%y = 0.0_wp
      self%z = 0.0_wp
      return
   end subroutine zero_p
   
   real(wp) function dot_p(self, vec) result(alpha)
      class(pos_vector)       , intent(in) :: self
      class(abstract_vector_rdp), intent(in) :: vec
      select type(vec)
      type is(pos_vector)
         alpha = self%x*vec%x + self%y*vec%y + self%z*vec%z
      end select
      return
   end function dot_p
   
   subroutine scal_p(self, alpha)
      class(pos_vector), intent(inout) :: self
      real(wp)           , intent(in)    :: alpha
      self%x = self%x * alpha
      self%y = self%y * alpha
      self%z = self%z * alpha
      return
   end subroutine scal_p
   
   subroutine axpby_p(self, alpha, vec, beta)
      class(pos_vector)       , intent(inout) :: self
      class(abstract_vector_rdp), intent(in)    :: vec
      real(wp)                  , intent(in)    :: alpha, beta
      select type(vec)
      type is(pos_vector)
         self%x = alpha*self%x + beta*vec%x
         self%y = alpha*self%y + beta*vec%y
         self%z = alpha*self%z + beta*vec%z
      end select
      return
   end subroutine axpby_p
   
   integer function get_size_p(self) result(N)
      class(pos_vector), intent(in) :: self
      N = npts+1
      return
   end function get_size_p
   
   subroutine rand_p(self, ifnorm)
      class(pos_vector), intent(inout) :: self
      logical, optional,   intent(in)    :: ifnorm
      logical :: normalized
      real(wp) :: mu, var
      real(wp) :: alpha
   
      mu = 0.0_wp
      var = 1.0_wp
      self%x = normal(mu, var)
      self%y = normal(mu, var)
      self%z = normal(mu, var)
 
      normalized = optval(ifnorm, .false.)
      if (normalized) then
         alpha = self%norm()
         call self%scal(1.0_wp/alpha)
      endif
      return
   end subroutine rand_p
   
   subroutine OTD_rhs(me, t, x, f)
      ! Time-integrator.
      class(rk_class), intent(inout)             :: me
      ! Current time.
      real(kind=wp)  , intent(in)                :: t
      ! State vector.
      real(kind=wp)  , dimension(:), intent(in)  :: x
      ! Time-derivative.
      real(kind=wp)  , dimension(:), intent(out) :: f
      ! internal
      real(kind=wp)  , dimension(npts)   :: bf
      real(kind=wp)  , dimension(npts,r) :: u, Lu, fp
      real(kind=wp)  , dimension(r,r)    :: Lr, Phi
      integer :: i, j
   
      bf = x(:npts)
      u = reshape(x(npts+1:(r+1)*npts), shape(u))
   
      call    nonlinear_roessler(        bf, f(:npts))
      do i = 1, r
         call linear_roessler   (u(:,i), bf, Lu(:,i))
      end do
      ! build reduced operator and rotation matrix
      Lr = 0.0_wp
      Phi = 0.0_wp
      do i = 1, r
         do j = 1, r
            Lr(i,j) = dot_product(Lu(:,i), u(:,j))
         end do
         do j = i+1, r
            Phi(i,j) = dot_product(Lu(:,i), u(:,j))
            Phi(j,i) = -Phi(i,j)
         enddo
      end do
      ! Construct forcing and add the rhs
      fp = Lu - matmul(u, Lr - Phi)
      f(npts+1:(r+1)*npts) = reshape(fp, [size(fp)])
      ! Construct dF/dt = Lr for the fundamental solution matrix
      f((r+1)*npts+1:) = reshape(Lr, [size(Lr)])
      return
   end subroutine OTD_rhs
   
   !------------------------------
   !-----     INTEGRATOR     -----
   !------------------------------

   subroutine OTD_step(integrator, bf, vec_in, F_in, time, Tstep, vec_out, F_out)
      ! Integrator
      class(rk_class)           , intent(inout)  :: integrator
      ! Basic state
      class(abstract_vector_rdp), intent(inout)  :: bf
      ! Input vector.
      class(abstract_vector_rdp), intent(in)     :: vec_in(r)
      ! Fundamental solution matrix
      real(wp),                   intent(in)     :: F_in(r,r)
      ! Current simulation time
      real(wp),                   intent(inout)  :: time
      ! Integration time for this step.
      real(wp),                   intent(in)     :: Tstep
      ! Output vector.
      class(abstract_vector_rdp), intent(out)    :: vec_out(r)
      ! Fundamental solution matrix
      real(wp),                   intent(out)    :: F_out(r,r)
      
      ! internals      
      real(wp)                             :: dt = 1.0_wp
      real(wp), dimension((r+1)*npts+r**2) :: pos_in, pos_out
      integer                              :: i
      
      select type(integrator)
      type is(rks54_class)
         select type(bf)
         type is(state_vector)
            select type(vec_in)
            type is(pos_vector)
               select type(vec_out)
               type is(pos_vector)
                  ! Get the state.
                  call get_position(bf, pos_in(:npts)) ! bf is a state vector
                  do i = 1, r
                     call get_pos(vec_in(i), pos_in(npts*i+1:npts*(i+1)))
                  end do
                  pos_in(npts*(r+1)+1:) = reshape(F_in, [size(F_in)])
                  ! Integrate.
                  call integrator%integrate(time, pos_in, dt, time+Tstep, pos_out)
                  ! Pass-back the state.
                  call set_position(pos_out(:npts), bf)
                  do i = 1, r
                     call set_pos(pos_out(npts*i+1:npts*(i+1)), vec_out(i))
                  end do
                  F_out = reshape(pos_out(npts*(r+1)+1:), shape(F_out))
                  time = time + Tstep     
               end select
            end select
         end select
      end select
      
      return
   end subroutine OTD_step

   subroutine OTD_map(bf, vec_in, Tend, vec_out, t_FTLE, if_report_stdout)
      ! Basic state
      class(abstract_vector_rdp), intent(inout) :: bf
      ! Input vector.
      class(abstract_vector_rdp), intent(inout) :: vec_in(r)
      ! Integration time.
      real(wp),                   intent(in)    :: Tend
      ! Output vector.
      class(abstract_vector_rdp), intent(out)   :: vec_out(r)
      ! Integration time for FLTEs
      real(wp),                   intent(in)    :: t_FTLE
      ! report?
      logical, optional,          intent(in)    :: if_report_stdout

      ! internals
      type(rks54_class)        :: OTD_roessler
      integer                  :: idx(1)
      real(wp)                 :: time, t_complete, t1, t2, tvec(2)
      logical                  :: if_GS
      integer                  :: i, j
      real(wp), dimension(r,r) :: F_in, F_out, G
      real(wp), dimension(r)   :: FTLEv

      integer, parameter :: ndof = npts*(r+1) + r**2

      ! In finite-precision arithmetic we need to reorthonormalize sometimes
      real(wp), parameter :: t_GS = 1.0_wp

      ! Initialize integrator.
      if (present(if_report_stdout)) then
         if (if_report_stdout) then
            call OTD_roessler%initialize(n=ndof, f=OTD_rhs, report=OTD_report_stdout, report_rate=100)
         else
            call OTD_roessler%initialize(n=ndof, f=OTD_rhs, report=OTD_report_file, report_rate=1)
         end if
      else
         call OTD_roessler%initialize(n=ndof, f=OTD_rhs)
      end if

      ! Integrate
      time = 0.0_wp
      tvec = (/ t_GS, t_FTLE /)
      idx = minloc(tvec)
      t1 = minval(tvec)
      t2 = maxval(tvec)
      F_in = eye(r)
      do j = 1, floor(Tend/t2)
         ! Regular integration with much less frequent reorthonormalization
         do i = 1, floor(t2/t1)
            call OTD_step(OTD_roessler, bf, vec_in, F_in, time, t1, vec_out, F_out)
            if (t1 == t_GS) then
               call orthonormalize_basis(vec_out)
               F_in = F_out
            else
               FTLEv = compute_FTLEs(F_out, t_FTLE)
               print *, 'FTLE (', t_FTLE,'): ', FTLEv
               F_in = eye(r)
            endif
            call copy_basis(vec_in, vec_out)
         end do
         ! Complete to the full integration time if not a multiple of TGS
         t_complete = modulo(t2, t1)
         if (t_complete > 1e-4_wp) call OTD_step(OTD_roessler, bf, vec_in, F_in, time, t_complete, vec_out, F_out)
         if (t2 == t_GS) then
            call orthonormalize_basis(vec_out)
            F_in = F_out
         else
            FTLEv = compute_FTLEs(F_out, t_FTLE)
            print *, 'FTLE (', t_FTLE,'): ', FTLEv
            F_in = eye(r)
         endif
         call copy_basis(vec_in, vec_out)
      end do
      t_complete = modulo(Tend, t2)
      if (t_complete > 1e-4_wp) call OTD_step(OTD_roessler, bf, vec_in, F_in, time, t_complete, vec_out, F_out)
      
      return
   end subroutine OTD_map

   function compute_FTLEs(Phi, T) result(FTLE)
      real(wp), dimension(r,r), intent(in) :: Phi
      real(wp),                 intent(in) :: T
      real(wp), dimension(r) :: FTLE
      ! internal
      real(wp), dimension(r,r) :: G
      G = matmul(transpose(Phi), Phi)
      FTLE = log(sqrt(svdvals(G)))/T
      !FTLE = log(sqrt(real(eigvals(G))))/T
   end function

   !-------------------------------------------
   !-----     MISCELLANEOUS UTILITIES     -----
   !-------------------------------------------
 
   subroutine get_pos(vec_in, pos)
      class(abstract_vector_rdp), intent(in)  :: vec_in
      real(wp), dimension(npts),  intent(out) :: pos
      pos = 0.0_wp
      select type (vec_in)
      type is (pos_vector)
          pos(1) = vec_in%x
          pos(2) = vec_in%y
          pos(3) = vec_in%z
      end select
      return
   end subroutine get_pos

   subroutine set_pos(pos, vec_out)
      real(wp), dimension(npts),  intent(in)  :: pos
      class(abstract_vector_rdp), intent(out) :: vec_out
      select type (vec_out)
      type is (pos_vector)
         vec_out%x = pos(1)
         vec_out%y = pos(2)
         vec_out%z = pos(3)
      end select
      return
   end subroutine set_pos

   subroutine OTD_report_stdout(me, t, x)
      class(rk_class), intent(inout)      :: me
      real(wp), intent(in)                :: t
      real(wp), dimension(:), intent(in)  :: x
      
      ! internal
      real(wp), dimension(npts)   :: bf
      real(wp), dimension(npts,r) :: u, Lu
      real(wp), dimension(r,r)    :: Lr, uTu, F
      real(wp), dimension(r)      :: FTLEv
      integer                     :: i, j

      bf = x(:npts)
      u = reshape(x(npts+1:(r+1)*npts), shape(u))
      F = reshape(x(npts*(r+1)+1:), shape(F))

      do i = 1, r
         call linear_roessler   (u(:,i), bf, Lu(:,i))
      end do
      ! build reduced operator
      Lr = 0.0_wp
      uTu = 0.0_wp
      do i = 1, r
         do j = 1, r
            Lr(i,j)  = dot_product(Lu(:,i), u(:,j))
            uTu(i,j) = dot_product( u(:,i), u(:,j))
         end do
      end do
      FTLEv = compute_FTLEs(F, t)

      print '(*(E16.9,1X))', t, bf, u, real(eigvals(Lr)), log10(FTLEv), ( uTu(i,i) - 1.0_wp, i = 1, r )

      return
   end subroutine OTD_report_stdout

   subroutine OTD_report_file(me, t, x)
      class(rk_class), intent(inout)      :: me
      real(wp), intent(in)                :: t
      real(wp), dimension(:), intent(in)  :: x
      
      ! internal
      real(wp), dimension(npts)   :: bf
      real(wp), dimension(npts,r) :: u, Lu
      real(wp), dimension(r,r)    :: Lr, uTu, F
      real(wp), dimension(r)      :: FTLEv
      integer                     :: i, j, iunit

      bf = x(:npts)
      u = reshape(x(npts+1:(r+1)*npts), shape(u))
      F = reshape(x(npts*(r+1)+1:), shape(F))

      do i = 1, r
         call linear_roessler   (u(:,i), bf, Lu(:,i))
      end do
      ! build reduced operator
      Lr = 0.0_wp
      uTu = 0.0_wp
      do i = 1, r
         do j = 1, r
            Lr(i,j)  = dot_product(Lu(:,i), u(:,j))
            uTu(i,j) = dot_product( u(:,i), u(:,j))
         end do
      end do
      FTLEv = compute_FTLEs(F, t)

      open(newunit=iunit, file=file, status='old', action='write', position='append')
      write(iunit, '(*(E16.9,1X))') t, bf, u, real(eigvals(Lr)),  log10(FTLEv), ( uTu(i,i) - 1.0_wp, i = 1, r )
      close(iunit)

      return
   end subroutine OTD_report_file

   subroutine write_header()
      ! internals
      integer :: i, iunit
      open(newunit=iunit, file=file, status='new', action='write')
      ! time, baseflow
      write(iunit,'(*(A16,1X))', ADVANCE='NO') 't', 'x_BF', 'y_BF', 'z_BF'
      ! basis vectors
      do i = 1, r 
         write(iunit,'(*(A15,I1,1X))', ADVANCE='NO') 'xp_', i, 'yp_', i, 'zp_', i
      end do
      ! instantaneous eigenvalues of the reduced operator
      do i = 1, r
         write(iunit,'(A15,I1,1X)', ADVANCE='NO') 'lambda_', i
      end do
      ! orthogonality of the basis vectors
      do i = 1, r
         write(iunit,'(A5,I1,A3,I0,A5,1X)', ADVANCE='NO') '<u_', i, ',u_', i, '> - 1'
      end do
      write(iunit,*) ''; close(iunit)
      return
   end subroutine write_header
 
end module Roessler_OTD