module Roessler
   ! Standard Library.
   use stdlib_stats_distribution_normal, only: normal => rvs_normal
   use stdlib_optval, only: optval
   ! RKLIB module for time integration.
   use rklib_module
   ! LightKrylov for linear algebra.
   use LightKrylov
   use LightKrylov_Logger
   use LightKrylov_Constants
   implicit none

   character(len=*), parameter, private :: this_module = 'Roessler'
   character(len=*), parameter :: report_file = 'example/roessler/roessler_output.txt'

   public :: a, b, c

   !------------------------------
   !-----     PARAMETERS     -----
   !------------------------------

   integer, parameter :: npts = 3
   real(dp), parameter :: a = 0.2_dp
   real(dp), parameter :: b = 0.2_dp
   real(dp), parameter :: c = 5.7_dp

   !-------------------------------------------
   !-----     LIGHTKRYLOV VECTOR TYPE     -----
   !-------------------------------------------

   type, extends(abstract_vector_rdp), public :: state_vector
      real(dp) :: x = 0.0_dp
      real(dp) :: y = 0.0_dp
      real(dp) :: z = 0.0_dp
      real(dp) :: T = 0.0_dp
   contains
      private
      procedure, pass(self), public :: zero
      procedure, pass(self), public :: dot
      procedure, pass(self), public :: scal
      procedure, pass(self), public :: axpby
      procedure, pass(self), public :: rand
      procedure, pass(self), public :: get_size
   end type state_vector

   !-------------------------------------------
   !-----     LIGHTKRYLOV SYSTEM TYPE     -----
   !-------------------------------------------

   type, extends(abstract_system_rdp), public :: roessler_upo
   contains
      private
      procedure, pass(self), public :: response => nonlinear_map
   end type roessler_upo

   type, extends(abstract_jacobian_linop_rdp), public :: jacobian
   contains
      private
      procedure, pass(self), public :: matvec => linear_map
      procedure, pass(self), public :: rmatvec => linear_map ! dummy, we do not need the adjoint of the jacobian
   end type jacobian

   type, extends(abstract_jacobian_linop_rdp), public :: floquet_operator
   contains
      private
      procedure, pass(self), public :: matvec => monodromy_map
      procedure, pass(self), public :: rmatvec => monodromy_map ! dummy, we do not need the adjoint of the monodromy
   end type floquet_operator

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

   subroutine zero(self)
      class(state_vector), intent(inout) :: self
      ! spatial coordinates of initial condition for orbit
      self%x = 0.0_dp
      self%y = 0.0_dp
      self%z = 0.0_dp
      ! period
      self%T = 0.0_dp
   end subroutine zero

   real(dp) function dot(self, vec) result(alpha)
      class(state_vector), intent(in) :: self
      class(abstract_vector_rdp), intent(in) :: vec
      select type (vec)
      type is (state_vector)
         alpha = self%x*vec%x + self%y*vec%y + self%z*vec%z + self%T*vec%T
      class default
         call type_error('vec','state_vector','IN',this_module,'dot')
      end select
   end function dot

   subroutine scal(self, alpha)
      class(state_vector), intent(inout) :: self
      real(dp), intent(in)    :: alpha
      self%x = self%x*alpha
      self%y = self%y*alpha
      self%z = self%z*alpha
      self%T = self%T*alpha
   end subroutine scal

   subroutine axpby(alpha, vec, beta, self)
      class(state_vector), intent(inout) :: self
      class(abstract_vector_rdp), intent(in)    :: vec
      real(dp), intent(in)    :: alpha, beta
      select type (vec)
      type is (state_vector)
         self%x = beta*self%x + alpha*vec%x
         self%y = beta*self%y + alpha*vec%y
         self%z = beta*self%z + alpha*vec%z
         self%T = beta*self%T + alpha*vec%T
      class default
         call type_error('vec','state_vector','IN',this_module,'axpby')
      end select
   end subroutine axpby

   integer function get_size(self) result(N)
      class(state_vector), intent(in) :: self
      N = npts + 1
   end function get_size

   subroutine rand(self, ifnorm)
      class(state_vector), intent(inout) :: self
      logical, optional, intent(in)    :: ifnorm
      logical :: normalized
      real(dp) :: mu, var
      real(dp) :: alpha

      mu = 0.0_dp
      var = 1.0_dp
      self%x = normal(mu, var)
      self%y = normal(mu, var)
      self%z = normal(mu, var)
      self%T = normal(mu, var)

      normalized = optval(ifnorm, .false.)
      if (normalized) then
         alpha = self%norm()
         call self%scal(1.0_dp/alpha)
      end if
   end subroutine rand

   !===================================
   !===================================
   !=====                         =====
   !=====     ROESSLER SYSTEM     =====
   !=====                         =====
   !===================================
   !===================================

   !----------------------------------------
   !-----      NONLINEAR EQUATIONS     -----
   !----------------------------------------

   subroutine nonlinear_roessler(x, f)
      ! State vector.
      real(dp), dimension(npts), intent(in)  :: x
      ! Time-derivative.
      real(dp), dimension(npts), intent(out) :: f

      f = 0.0_dp
      f(1) = -x(2) - x(3)
      f(2) = x(1) + a*x(2)
      f(3) = b + x(3)*(x(1) - c)

   end subroutine nonlinear_roessler

   !-----------------------------
   !-----      JACOBIAN     -----
   !-----------------------------

   subroutine linear_roessler(xp, bf, f)
      ! State vector.
      real(dp), dimension(:), intent(in)  :: xp
      ! Base state.
      real(dp), dimension(:), intent(in)  :: bf
      ! Time-derivative.
      real(dp), dimension(:), intent(out) :: f

      f = 0.0_dp
      f(1) = -xp(2) - xp(3)
      f(2) = xp(1) + a*xp(2)
      f(3) = xp(1)*bf(3) + xp(3)*(bf(1) - c)

   end subroutine linear_roessler

   !-----------------------------------------------------------------
   !-----      WRAPPER FOR NONLINEAR INTEGRATION WITH RKLIB     -----
   !-----------------------------------------------------------------

   subroutine NL_rhs(me, t, x, f)
      ! Time-integrator.
      class(rk_class), intent(inout)             :: me
      ! Current time.
      real(dp), intent(in)                :: t
      ! State vector.
      real(dp), dimension(:), intent(in)  :: x
      ! Time-derivative.
      real(dp), dimension(:), intent(out) :: f

      call nonlinear_roessler(x, f)

   end subroutine NL_rhs

   !-------------------------------------------------------------
   !-----      WRAPPER FOR LINEAR INTEGRATION WITH RKLIB    -----
   !-------------------------------------------------------------

   subroutine combined_rhs(me, t, x, f)
      ! Time-integrator.
      class(rk_class), intent(inout)             :: me
      ! Current time.
      real(dp), intent(in)                :: t
      ! State vector.
      real(dp), dimension(:), intent(in)  :: x
      ! Time-derivative.
      real(dp), dimension(:), intent(out) :: f

      call nonlinear_roessler(x(:npts), f(:npts))
      call linear_roessler(x(npts + 1:), x(:npts), f(npts + 1:))

   end subroutine combined_rhs

   !-------------------------------------------------------------
   !-----     TYPE-BOUND PROCEDURES FOR THE INTEGRATORS     -----
   !-------------------------------------------------------------

   subroutine nonlinear_map(self, vec_in, vec_out, atol)
      ! Dynamical system.
      class(roessler_upo), intent(inout)  :: self
      ! Input vector.
      class(abstract_vector_rdp), intent(in)  :: vec_in
      ! Output vector.
      class(abstract_vector_rdp), intent(out) :: vec_out
      ! Solver tolerances if needed
      real(dp), intent(in)  :: atol

      ! Time-integrator.
      type(rks54_class)         :: nonlinear_integrator
      real(dp)                  :: dt = 1.0_dp
      real(dp), dimension(npts) :: pos_in, pos_out

      ! Evaluate F(X).
      select type (vec_in)
      type is (state_vector)
         select type (vec_out)
         type is (state_vector)
            ! Get state vector.
            call get_position(vec_in, pos_in)
            ! Initialize integrator.
            call nonlinear_integrator%initialize(n=npts, f=NL_rhs)!, report=roessler_report_stdout, report_rate=1)
            ! Integrate forward in time.
            call nonlinear_integrator%integrate(0.0_dp, pos_in, dt, vec_in%T, pos_out)
            ! Pass-back the state vector.
            call set_position(pos_out, vec_out)

            ! Evaluate residual F(X) - X.
            call vec_out%sub(vec_in)

            ! Add period residual
            vec_out%T = 0.0_dp
         class default
            call type_error('vec_out','state_vector','OUT',this_module,'nonlinear_map')
         end select
      class default
         call type_error('vec_in','state_vector','IN',this_module,'nonlinear_map')
      end select
   end subroutine nonlinear_map

   subroutine linear_map(self, vec_in, vec_out)
      ! Linear Operator.
      class(jacobian), intent(inout) :: self
      ! Input vector.
      class(abstract_vector_rdp), intent(in)  :: vec_in
      ! Output vector.
      class(abstract_vector_rdp), intent(out) :: vec_out

      ! Time-integrator.
      type(rks54_class)           :: combined_roessler
      real(dp)                    :: dt = 1.0_dp
      real(dp)                    :: period
      real(dp), dimension(2*npts) :: pos_in, pos_out
      type(state_vector)          :: vec

      select type (vec_in)
      type is (state_vector)
         select type (vec_out)
         type is (state_vector)
            ! Get the state.
            call get_position(self%X, pos_in(:npts))
            call get_period(self%X, period)
            call get_position(vec_in, pos_in(npts + 1:))
            ! Initialize integrator.
            call combined_roessler%initialize(n=2*npts, f=combined_rhs)!, report=roessler_report_stdout, report_rate=1)
            ! Evaluate:
            ! 1. F(X)
            ! 2. exp(tau*J) @ dx
            call combined_roessler%integrate(0.0_dp, pos_in, dt, period, pos_out)
            ! Pass-back the state.
            call set_position(pos_out(npts + 1:), vec_out)

            ! Evaluate [ exp(tau*J) - I ] @ dx.
            call vec_out%sub(vec_in)

            ! Evaluate f'(X(T), T) * dT and add it to the position residual
            call compute_fdot(pos_out(:npts), vec)
            call vec_out%axpby(vec_in%T, vec, 1.0_dp)

            ! Evaluate f'(X(0), 0).T @ dx and add phase condition
            call compute_fdot(pos_in(:npts), vec)
            vec_out%T = vec_in%dot(vec)
         class default
            call type_error('vec_out','state_vector','OUT',this_module,'linear_map')
         end select
      class default
         call type_error('vec_in','state_vector','IN',this_module,'linear_map')
      end select
   end subroutine linear_map

   subroutine monodromy_map(self, vec_in, vec_out)
      ! Linear Operator.
      class(floquet_operator), intent(inout) :: self
      ! Input vector.
      class(abstract_vector_rdp), intent(in)  :: vec_in
      ! Output vector.
      class(abstract_vector_rdp), intent(out) :: vec_out

      ! Time-integrator.
      type(rks54_class)           :: combined_roessler
      real(dp)                    :: dt = 1.0_dp
      real(dp)                    :: period
      real(dp), dimension(2*npts) :: pos_in, pos_out

      select type (vec_in)
      type is (state_vector)
         select type (vec_out)
         type is (state_vector)
            ! Get the state.
            call get_position(self%X, pos_in(:npts))
            call get_period(self%X, period)
            call get_position(vec_in, pos_in(npts + 1:))
            ! Initialize integrator.
            call combined_roessler%initialize(n=2*npts, f=combined_rhs)!, report=roessler_report_stdout, report_rate=1)
            ! Evaluate:
            ! 1. F(X)
            ! 2. exp(tau*J) @ dx
            call combined_roessler%integrate(0.0_dp, pos_in, dt, period, pos_out)
            ! Pass-back the state.
            call set_position(pos_out(npts + 1:), vec_out)
         class default
            call type_error('vec_out','state_vector','OUT',this_module,'monodromy_map')
         end select
      class default
         call type_error('vec_in','state_vector','IN',this_module,'monodromy_map')
      end select
   end subroutine monodromy_map

   !-------------------------------------------
   !-----     MISCELLANEOUS UTILITIES     -----
   !-------------------------------------------

   subroutine get_position(vec_in, pos)
      class(abstract_vector_rdp), intent(in)  :: vec_in
      real(dp), dimension(npts), intent(out) :: pos
      pos = 0.0_dp
      select type (vec_in)
      type is (state_vector)
         pos(1) = vec_in%x
         pos(2) = vec_in%y
         pos(3) = vec_in%z
      class default
         call type_error('vec','state_vector','IN',this_module,'get_position')
      end select
   end subroutine get_position

   subroutine get_period(vec_in, period)
      class(abstract_vector_rdp), intent(in)  :: vec_in
      real(dp), intent(out) :: period
      select type (vec_in)
      type is (state_vector)
         period = vec_in%T
      class default
         call type_error('vec_in','state_vector','IN',this_module,'get_period')
      end select
   end subroutine get_period

   subroutine set_position(pos, vec_out)
      real(dp), dimension(npts), intent(in)  :: pos
      class(abstract_vector_rdp), intent(out) :: vec_out
      select type (vec_out)
      type is (state_vector)
         vec_out%x = pos(1)
         vec_out%y = pos(2)
         vec_out%z = pos(3)
      class default
         call type_error('vec_out','state_vector','OUT',this_module,'set_position')
      end select
   end subroutine set_position

   subroutine compute_fdot(pos, vec_out)
      real(dp), dimension(npts), intent(in)  :: pos
      class(abstract_vector_rdp), intent(out) :: vec_out
      ! internal
      real(dp), dimension(npts) :: wrk
      call vec_out%zero() ! ensure that the period shift vec_out%T is zero
      call nonlinear_roessler(pos, wrk)
      call set_position(wrk, vec_out)
   end subroutine compute_fdot

   subroutine roessler_report_stdout(me, t, x)
      class(rk_class), intent(inout)      :: me
      real(dp), intent(in)                :: t
      real(dp), dimension(:), intent(in)  :: x

      print '(*(F15.6,1X))', t, x

   end subroutine roessler_report_stdout

   subroutine write_report_header()
      ! internals
      integer :: iunit
      open (newunit=iunit, file=report_file, status='new', action='write')
      ! time, baseflow
      write (iunit, '(*(A16,1X))') 't', 'BF_x', 'BF_y', 'BF_z'
      close (iunit)
   end subroutine write_report_header

   subroutine roessler_report_file(me, t, x)
      class(rk_class), intent(inout)      :: me
      real(dp), intent(in)                :: t
      real(dp), dimension(:), intent(in)  :: x
      ! internals
      integer :: iunit
      open (newunit=iunit, file=report_file, status='old', action='write', position='append')
      write (iunit, '(*(F15.6,1X))') t, x
      close (iunit)
   end subroutine roessler_report_file

end module Roessler
