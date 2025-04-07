module Ginzburg_Landau
   ! RKLIB module for time integration.
   use rklib_module
   ! LightKrylov for linear algebra.
   use LightKrylov
   use LightKrylov_Logger
   ! Standard Library.
   use stdlib_math, only: linspace
   use stdlib_optval, only: optval
   implicit none

   character(len=*), parameter, private :: this_module = 'Ginzburg_Landau'

   public :: nx
   public :: initialize_parameters

   !------------------------------
   !-----     PARAMETERS     -----
   !------------------------------

   ! Mesh related parameters.
   real(kind=dp), parameter :: L = 200.0_dp ! Domain length
   integer, parameter :: nx = 512      ! Number of grid points (excluding boundaries).
   real(kind=dp), parameter :: dx = L/(nx + 1) ! Grid size.

   ! Physical parameters.
   complex(kind=dp), parameter :: nu = cmplx(2.0_dp, 0.2_dp, kind=dp)
   complex(kind=dp), parameter :: gamma = cmplx(1.0_dp, -1.0_dp, kind=dp)
   real(kind=dp), parameter :: mu_0 = 0.38_dp
   real(kind=dp), parameter :: c_mu = 0.2_dp
   real(kind=dp), parameter :: mu_2 = -0.01_dp
   real(kind=dp)               :: mu(1:nx)

   !-------------------------------------------
   !-----     LIGHTKRYLOV VECTOR TYPE     -----
   !-------------------------------------------

   type, extends(abstract_vector_cdp), public :: state_vector
      complex(kind=dp) :: state(nx) = 0.0_dp
   contains
      private
      procedure, pass(self), public :: zero
      procedure, pass(self), public :: dot
      procedure, pass(self), public :: scal
      procedure, pass(self), public :: axpby
      procedure, pass(self), public :: rand
      procedure, pass(self), public :: get_size
   end type state_vector

   !------------------------------------------
   !-----     EXPONENTIAL PROPAGATOR     -----
   !------------------------------------------

   type, extends(abstract_exptA_linop_cdp), public :: exponential_prop
   contains
      private
      procedure, pass(self), public :: matvec => direct_solver
      procedure, pass(self), public :: rmatvec => adjoint_solver
   end type exponential_prop

   interface exponential_prop
      module function construct_exptA(tau) result(A)
         real(kind=dp), intent(in) :: tau
         type(exponential_prop) :: A
      end function
   end interface

contains

   module procedure construct_exptA
      A%tau = tau
   end procedure

   !========================================================================
   !========================================================================
   !=====                                                              =====
   !=====     PHYSICAL MODEL : LINEARIZED GINZBURG-LANDAU EQUATION     =====
   !=====                                                              =====
   !========================================================================
   !========================================================================

   !--------------------------------------------------------------
   !-----     CONSTRUCT THE MESH AND PHYSICAL PARAMETERS     -----
   !--------------------------------------------------------------

   subroutine initialize_parameters()
      implicit none
      ! Mesh array.
      real(kind=dp), allocatable :: x(:)

      ! Construct mesh.
      x = linspace(-L/2, L/2, nx + 2)

      ! Construct mu(x)
      mu(:) = (mu_0 - c_mu**2) + (mu_2/2.0_dp)*x(2:nx + 1)**2
   end subroutine initialize_parameters

   !---------------------------------------------------------
   !-----      LINEARIZED GINZBURG-LANDAU EQUATIONS     -----
   !---------------------------------------------------------

   subroutine rhs(me, t, x, f)
      ! Time-integrator.
      class(rk_class), intent(inout)             :: me
      ! Current time.
      real(kind=dp), intent(in)                :: t
      ! State vector.
      real(kind=dp), dimension(:), intent(in)  :: x
      ! Time-derivative.
      real(kind=dp), dimension(:), intent(out) :: f

      ! Internal variables.
      integer :: i
      real(kind=dp), dimension(nx) :: u, du
      real(kind=dp), dimension(nx) :: v, dv
      real(kind=dp)                :: d2u, d2v, cu, cv

      ! Sets the internal variables.
      f = 0.0_dp
      u = x(1:nx); du = f(1:nx)
      v = x(nx + 1:2*nx); dv = f(nx + 1:2*nx)

      !---------------------------------------------------
      !-----     Linear Ginzburg Landau Equation     -----
      !---------------------------------------------------

      ! Left most boundary points.
      cu = u(2)/(2*dx); cv = v(2)/(2*dx)
      du(1) = -(real(nu)*cu - aimag(nu)*cv) ! Convective term.
      dv(1) = -(aimag(nu)*cu + real(nu)*cv) ! Convective term.

      d2u = (u(2) - 2*u(1))/dx**2; d2v = (v(2) - 2*v(1))/dx**2
      du(1) = du(1) + real(gamma)*d2u - aimag(gamma)*d2v ! Diffusion term.
      dv(1) = dv(1) + aimag(gamma)*d2u + real(gamma)*d2v ! Diffusion term.

      du(1) = du(1) + mu(1)*u(1) ! Non-parallel term.
      dv(1) = dv(1) + mu(1)*v(1) ! Non-parallel term.

      ! Interior nodes.
      do i = 2, nx - 1
         ! Convective term.
         cu = (u(i + 1) - u(i - 1))/(2*dx)
         cv = (v(i + 1) - v(i - 1))/(2*dx)
         du(i) = -(real(nu)*cu - aimag(nu)*cv)
         dv(i) = -(aimag(nu)*cu + real(nu)*cv)

         ! Diffusion term.
         d2u = (u(i + 1) - 2*u(i) + u(i - 1))/dx**2
         d2v = (v(i + 1) - 2*v(i) + v(i - 1))/dx**2
         du(i) = du(i) + real(gamma)*d2u - aimag(gamma)*d2v
         dv(i) = dv(i) + aimag(gamma)*d2u + real(gamma)*d2v

         ! Non-parallel term.
         du(i) = du(i) + mu(i)*u(i)
         dv(i) = dv(i) + mu(i)*v(i)
      end do

      ! Right most boundary points.
      cu = -u(nx - 1)/(2*dx); cv = -v(nx - 1)/(2*dx)
      du(nx) = -(real(nu)*cu - aimag(nu)*cv) ! Convective term.
      dv(nx) = -(aimag(nu)*cu + real(nu)*cv) ! Convective term.

      d2u = (-2*u(nx) + u(nx - 1))/dx**2; d2v = (-2*v(nx) + v(nx - 1))/dx**2
      du(nx) = du(nx) + real(gamma)*d2u - aimag(gamma)*d2v ! Diffusion term.
      dv(nx) = dv(nx) + aimag(gamma)*d2u + real(gamma)*d2v ! Diffusion term.

      du(nx) = du(nx) + mu(nx)*u(nx) ! Non-parallel term.
      dv(nx) = dv(nx) + mu(nx)*v(nx) ! Non-parallel term.

      ! Copy results to the output array.
      f(1:nx) = du; f(nx + 1:2*nx) = dv
   end subroutine rhs

   !-----------------------------------------------------------
   !-----     Adjoint linear Ginzburg-Landau equation     -----
   !-----------------------------------------------------------

   subroutine adjoint_rhs(me, t, x, f)
      ! Time-integrator.
      class(rk_class), intent(inout)           :: me
      ! Current time.
      real(kind=dp), intent(in)                :: t
      ! State vector.
      real(kind=dp), dimension(:), intent(in)  :: x
      ! Time-derivative.
      real(kind=dp), dimension(:), intent(out) :: f

      ! Internal variables.
      integer :: i
      real(kind=dp), dimension(nx) :: u, du
      real(kind=dp), dimension(nx) :: v, dv
      real(kind=dp)                :: d2u, d2v, cu, cv

      ! Sets the internal variables.
      f = 0.0_dp
      u = x(1:nx); du = f(1:nx)
      v = x(nx + 1:2*nx); dv = f(nx + 1:2*nx)

      !---------------------------------------------------
      !-----     Linear Ginzburg Landau Equation     -----
      !---------------------------------------------------

      ! Left most boundary points.
      cu = u(2)/(2*dx); cv = v(2)/(2*dx)
      du(1) = (real(nu)*cu + aimag(nu)*cv) ! Convective term.
      dv(1) = (-aimag(nu)*cu + real(nu)*cv) ! Convective term.

      d2u = (u(2) - 2*u(1))/dx**2; d2v = (v(2) - 2*v(1))/dx**2
      du(1) = du(1) + real(gamma)*d2u + aimag(gamma)*d2v ! Diffusion term.
      dv(1) = dv(1) - aimag(gamma)*d2u + real(gamma)*d2v ! Diffusion term.

      du(1) = du(1) + mu(1)*u(1) ! Non-parallel term.
      dv(1) = dv(1) + mu(1)*v(1) ! Non-parallel term.

      ! Interior nodes.
      do i = 2, nx - 1
         ! Convective term.
         cu = (u(i + 1) - u(i - 1))/(2*dx)
         cv = (v(i + 1) - v(i - 1))/(2*dx)
         du(i) = (real(nu)*cu + aimag(nu)*cv)
         dv(i) = (-aimag(nu)*cu + real(nu)*cv)

         ! Diffusion term.
         d2u = (u(i + 1) - 2*u(i) + u(i - 1))/dx**2
         d2v = (v(i + 1) - 2*v(i) + v(i - 1))/dx**2
         du(i) = du(i) + real(gamma)*d2u + aimag(gamma)*d2v
         dv(i) = dv(i) - aimag(gamma)*d2u + real(gamma)*d2v

         ! Non-parallel term.
         du(i) = du(i) + mu(i)*u(i)
         dv(i) = dv(i) + mu(i)*v(i)
      end do

      ! Right most boundary points.
      cu = -u(nx - 1)/(2*dx); cv = -v(nx - 1)/(2*dx)
      du(nx) = (real(nu)*cu + aimag(nu)*cv) ! Convective term.
      dv(nx) = (-aimag(nu)*cu + real(nu)*cv) ! Convective term.

      d2u = (-2*u(nx) + u(nx - 1))/dx**2; d2v = (-2*v(nx) + v(nx - 1))/dx**2
      du(nx) = du(nx) + real(gamma)*d2u + aimag(gamma)*d2v ! Diffusion term.
      dv(nx) = dv(nx) - aimag(gamma)*d2u + real(gamma)*d2v ! Diffusion term.

      du(nx) = du(nx) + mu(nx)*u(nx) ! Non-parallel term.
      dv(nx) = dv(nx) + mu(nx)*v(nx) ! Non-parallel term.

      ! Copy results to the output array.
      f(1:nx) = du; f(nx + 1:2*nx) = dv
   end subroutine adjoint_rhs

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
      self%state = 0.0_dp
   end subroutine zero

   complex(kind=dp) function dot(self, vec) result(alpha)
      class(state_vector), intent(in) :: self
      class(abstract_vector_cdp), intent(in) :: vec
      select type (vec)
      type is (state_vector)
         alpha = dot_product(self%state, vec%state)
      class default
         call type_error('vec','state_vector','IN',this_module,'dot')
      end select
   end function dot

   subroutine scal(self, alpha)
      class(state_vector), intent(inout) :: self
      complex(kind=dp), intent(in)    :: alpha
      self%state = self%state*alpha
   end subroutine scal

   subroutine axpby(alpha, vec, beta, self)
      class(state_vector), intent(inout) :: self
      class(abstract_vector_cdp), intent(in)    :: vec
      complex(kind=dp), intent(in)    :: alpha, beta
      select type (vec)
      type is (state_vector)
         self%state = beta*self%state + alpha*vec%state
      class default
         call type_error('vec','state_vector','IN',this_module,'axpby')
      end select
   end subroutine axpby

   integer function get_size(self) result(N)
      class(state_vector), intent(in) :: self
      N = nx
   end function get_size

   subroutine rand(self, ifnorm)
      class(state_vector), intent(inout) :: self
      logical, optional, intent(in)    :: ifnorm
      real(kind=dp) :: tmp(nx, 2)
      ! internals
      logical :: normalize
      real(kind=dp) :: alpha
      normalize = optval(ifnorm, .true.)
      call random_number(tmp)
      self%state%re = tmp(:, 1); self%state%im = tmp(:, 2)
      if (normalize) then
         alpha = self%norm()
         call self%scal(cmplx(1.0_dp, 0.0_dp, kind=dp)/alpha)
      end if
   end subroutine rand

   !------------------------------------------------------------------------
   !-----     TYPE-BOUND PROCEDURES FOR THE EXPONENTIAL PROPAGATOR     -----
   !------------------------------------------------------------------------

   subroutine direct_solver(self, vec_in, vec_out)
      ! Linear Operator.
      class(exponential_prop), intent(inout)  :: self
      ! Input vector.
      class(abstract_vector_cdp), intent(in)  :: vec_in
      ! Output vector.
      class(abstract_vector_cdp), intent(out) :: vec_out

      ! Time-integrator.
      type(rks54_class) :: prop
      real(kind=dp)     :: dt = 1.0_dp
      real(kind=dp)     :: state_ic(2*nx), state_fc(2*nx)

      select type (vec_in)
      type is (state_vector)
         select type (vec_out)
         type is (state_vector)
            ! Get state vector.
            state_ic(:nx) = vec_in%state%re
            state_ic(nx + 1:) = vec_in%state%im
            ! Initialize propagator.
            call prop%initialize(n=2*nx, f=rhs)
            ! Integrate forward in time.
            call prop%integrate(0.0_dp, state_ic, dt, self%tau, state_fc)
            ! Pass-back the state vector.
            vec_out%state%re = state_fc(:nx)
            vec_out%state%im = state_fc(nx + 1:)
         class default
            call type_error('vec_out','state_vector','OUT',this_module,'direct_solver')
         end select
      class default
         call type_error('vec_in','state_vector','IN',this_module,'direct_solver')
      end select
   end subroutine direct_solver

   subroutine adjoint_solver(self, vec_in, vec_out)
      ! Linear Operator.
      class(exponential_prop), intent(inout)  :: self
      ! Input vector.
      class(abstract_vector_cdp), intent(in)  :: vec_in
      ! Output vector.
      class(abstract_vector_cdp), intent(out) :: vec_out

      ! Time-integrator.
      type(rks54_class) :: prop
      real(kind=dp)     :: dt = 1.0_dp
      real(kind=dp)     :: state_ic(2*nx), state_fc(2*nx)

      select type (vec_in)
      type is (state_vector)
         select type (vec_out)
         type is (state_vector)
            ! Get the state.
            state_ic(:nx) = vec_in%state%re
            state_fc(nx + 1:) = vec_in%state%im
            ! Initialize propagator.
            call prop%initialize(n=2*nx, f=adjoint_rhs)
            ! Integrate forward in time.
            call prop%integrate(0.0_dp, state_ic, dt, self%tau, state_fc)
            ! Pass-back the state.
            vec_out%state%re = state_fc(:nx)
            vec_out%state%im = state_fc(nx + 1:)
         class default
            call type_error('vec_out','state_vector','OUT',this_module,'adjoint_solver')
         end select
      class default
         call type_error('vec_in','state_vector','IN',this_module,'adjoint_solver')
      end select
   end subroutine adjoint_solver

end module Ginzburg_Landau
