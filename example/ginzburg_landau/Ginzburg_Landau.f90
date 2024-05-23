module Ginzburg_Landau
  ! RKLIB module for time integration.
  use rklib_module
  ! LightKrylov for linear algebra.
  use LightKrylov
  use LightKrylov, only: wp => dp
  ! Standard Library.
  use stdlib_math, only : linspace
  use stdlib_optval, only : optval
  implicit none

  private
  public :: nx
  public :: initialize_parameters

  !------------------------------
  !-----     PARAMETERS     -----
  !------------------------------

  ! Mesh related parameters.
  real(kind=wp), parameter :: L  = 200.0_wp ! Domain length
  integer      , parameter :: nx = 512      ! Number of grid points (excluding boundaries).
  real(kind=wp), parameter :: dx = L/(nx+1) ! Grid size.

  ! Physical parameters.
  complex(kind=wp), parameter :: nu    = cmplx(2.0_wp, 0.2_wp, kind=wp)
  complex(kind=wp), parameter :: gamma = cmplx(1.0_wp, -1.0_wp, kind=wp)
  real(kind=wp)   , parameter :: mu_0  = 0.38_wp
  real(kind=wp)   , parameter :: c_mu  = 0.2_wp
  real(kind=wp)   , parameter :: mu_2  = -0.01_wp
  real(kind=wp)               :: mu(1:nx)

  !-------------------------------------------
  !-----     LIGHTKRYLOV VECTOR TYPE     -----
  !-------------------------------------------

  type, extends(abstract_vector_cdp), public :: state_vector
     complex(kind=wp) :: state(nx) = 0.0_wp
   contains
     private
     procedure, pass(self), public :: zero
     procedure, pass(self), public :: dot
     procedure, pass(self), public :: scal
     procedure, pass(self), public :: axpby
     procedure, pass(self), public :: rand
  end type state_vector

  !------------------------------------------
  !-----     EXPONENTIAL PROPAGATOR     -----
  !------------------------------------------

  type, extends(abstract_linop_cdp), public :: exponential_prop
     real(kind=wp) :: tau ! Integration time.
   contains
     private
     procedure, pass(self), public :: matvec => direct_solver
     procedure, pass(self), public :: rmatvec => adjoint_solver
  end type exponential_prop

contains

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
    real(kind=wp), allocatable :: x(:)

    ! Construct mesh.
    x = linspace(-L/2, L/2, nx+2)

    ! Construct mu(x)
    mu(:) = (mu_0 - c_mu**2) + (mu_2 / 2.0_wp) * x(2:nx+1)**2

    return
  end subroutine initialize_parameters


  !---------------------------------------------------------
  !-----      LINEARIZED GINZBURG-LANDAU EQUATIONS     -----
  !---------------------------------------------------------

  subroutine rhs(me, t, x, f)
    ! Time-integrator.
    class(rk_class), intent(inout)             :: me
    ! Current time.
    real(kind=wp)  , intent(in)                :: t
    ! State vector.
    real(kind=wp)  , dimension(:), intent(in)  :: x
    ! Time-derivative.
    real(kind=wp)  , dimension(:), intent(out) :: f

    ! Internal variables.
    integer :: i, j, k
    real(kind=wp), dimension(nx) :: u, du
    real(kind=wp), dimension(nx) :: v, dv
    real(kind=wp)                :: d2u, d2v, cu, cv

    ! Sets the internal variables.
    f = 0.0_wp
    u = x(1:nx)      ; du = f(1:nx)
    v = x(nx+1:2*nx) ; dv = f(nx+1:2*nx)

    !---------------------------------------------------
    !-----     Linear Ginzburg Landau Equation     -----
    !---------------------------------------------------

    ! Left most boundary points.
    cu = u(2) / (2*dx) ; cv = v(2) / (2*dx)
    du(1) = -(real(nu)*cu - aimag(nu)*cv) ! Convective term.
    dv(1) = -(aimag(nu)*cu + real(nu)*cv) ! Convective term.

    d2u = (u(2) - 2*u(1)) / dx**2 ; d2v = (v(2) - 2*v(1)) / dx**2
    du(1) = du(1) + real(gamma)*d2u - aimag(gamma)*d2v ! Diffusion term.
    dv(1) = dv(1) + aimag(gamma)*d2u + real(gamma)*d2v ! Diffusion term.

    du(1) = du(1) + mu(1)*u(1) ! Non-parallel term.
    dv(1) = dv(1) + mu(1)*v(1) ! Non-parallel term.

    ! Interior nodes.
    do i = 2, nx-1
       ! Convective term.
       cu = (u(i+1) - u(i-1)) / (2*dx)
       cv = (v(i+1) - v(i-1)) / (2*dx)
       du(i) = -(real(nu)*cu - aimag(nu)*cv)
       dv(i) = -(aimag(nu)*cu + real(nu)*cv)

       ! Diffusion term.
       d2u = (u(i+1) - 2*u(i) + u(i-1)) / dx**2
       d2v = (v(i+1) - 2*v(i) + v(i-1)) / dx**2
       du(i) = du(i) + real(gamma)*d2u - aimag(gamma)*d2v
       dv(i) = dv(i) + aimag(gamma)*d2u + real(gamma)*d2v

       ! Non-parallel term.
       du(i) = du(i) + mu(i)*u(i)
       dv(i) = dv(i) + mu(i)*v(i)
    enddo

    ! Right most boundary points.
    cu = -u(nx-1) / (2*dx) ; cv = -v(nx-1) / (2*dx)
    du(nx) = -(real(nu)*cu - aimag(nu)*cv) ! Convective term.
    dv(nx) = -(aimag(nu)*cu + real(nu)*cv) ! Convective term.

    d2u = (-2*u(nx) + u(nx-1)) / dx**2 ; d2v = (-2*v(nx) + v(nx-1)) / dx**2
    du(nx) = du(nx) + real(gamma)*d2u - aimag(gamma)*d2v ! Diffusion term.
    dv(nx) = dv(nx) + aimag(gamma)*d2u + real(gamma)*d2v ! Diffusion term.

    du(nx) = du(nx) + mu(nx)*u(nx) ! Non-parallel term.
    dv(nx) = dv(nx) + mu(nx)*v(nx) ! Non-parallel term.

    ! Copy results to the output array.
    f(1:nx) = du ; f(nx+1:2*nx) = dv

    return
  end subroutine rhs

  !-----------------------------------------------------------
  !-----     Adjoint linear Ginzburg-Landau equation     -----
  !-----------------------------------------------------------

  subroutine adjoint_rhs(me, t, x, f)
    ! Time-integrator.
    class(rk_class), intent(inout)             :: me
    ! Current time.
    real(kind=wp)  , intent(in)                :: t
    ! State vector.
    real(kind=wp)  , dimension(:), intent(in)  :: x
    ! Time-derivative.
    real(kind=wp)  , dimension(:), intent(out) :: f

    ! Internal variables.
    integer :: i, j, k
    real(kind=wp), dimension(nx) :: u, du
    real(kind=wp), dimension(nx) :: v, dv
    real(kind=wp)                :: d2u, d2v, cu, cv

    ! Sets the internal variables.
    f = 0.0_wp
    u = x(1:nx)      ; du = f(1:nx)
    v = x(nx+1:2*nx) ; dv = f(nx+1:2*nx)

    !---------------------------------------------------
    !-----     Linear Ginzburg Landau Equation     -----
    !---------------------------------------------------

    ! Left most boundary points.
    cu = u(2) / (2*dx) ; cv = v(2) / (2*dx)
    du(1) = (real(nu)*cu + aimag(nu)*cv) ! Convective term.
    dv(1) = (-aimag(nu)*cu + real(nu)*cv) ! Convective term.

    d2u = (u(2) - 2*u(1)) / dx**2 ; d2v = (v(2) - 2*v(1)) / dx**2
    du(1) = du(1) + real(gamma)*d2u + aimag(gamma)*d2v ! Diffusion term.
    dv(1) = dv(1) - aimag(gamma)*d2u + real(gamma)*d2v ! Diffusion term.

    du(1) = du(1) + mu(1)*u(1) ! Non-parallel term.
    dv(1) = dv(1) + mu(1)*v(1) ! Non-parallel term.

    ! Interior nodes.
    do i = 2, nx-1
       ! Convective term.
       cu = (u(i+1) - u(i-1)) / (2*dx)
       cv = (v(i+1) - v(i-1)) / (2*dx)
       du(i) = (real(nu)*cu + aimag(nu)*cv)
       dv(i) = (-aimag(nu)*cu + real(nu)*cv)

       ! Diffusion term.
       d2u = (u(i+1) - 2*u(i) + u(i-1)) / dx**2
       d2v = (v(i+1) - 2*v(i) + v(i-1)) / dx**2
       du(i) = du(i) + real(gamma)*d2u + aimag(gamma)*d2v
       dv(i) = dv(i) - aimag(gamma)*d2u + real(gamma)*d2v

       ! Non-parallel term.
       du(i) = du(i) + mu(i)*u(i)
       dv(i) = dv(i) + mu(i)*v(i)
    enddo

    ! Right most boundary points.
    cu = -u(nx-1) / (2*dx) ; cv = -v(nx-1) / (2*dx)
    du(nx) = (real(nu)*cu + aimag(nu)*cv) ! Convective term.
    dv(nx) = (-aimag(nu)*cu + real(nu)*cv) ! Convective term.

    d2u = (-2*u(nx) + u(nx-1)) / dx**2 ; d2v = (-2*v(nx) + v(nx-1)) / dx**2
    du(nx) = du(nx) + real(gamma)*d2u + aimag(gamma)*d2v ! Diffusion term.
    dv(nx) = dv(nx) - aimag(gamma)*d2u + real(gamma)*d2v ! Diffusion term.

    du(nx) = du(nx) + mu(nx)*u(nx) ! Non-parallel term.
    dv(nx) = dv(nx) + mu(nx)*v(nx) ! Non-parallel term.

    ! Copy results to the output array.
    f(1:nx) = du ; f(nx+1:2*nx) = dv

    return
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
    self%state = 0.0_wp
    return
  end subroutine zero

  complex(kind=wp) function dot(self, vec) result(alpha)
    class(state_vector)   , intent(in) :: self
    class(abstract_vector_cdp), intent(in) :: vec
    select type(vec)
    type is(state_vector)
       alpha = dot_product(self%state, vec%state)
    end select
    return
  end function dot

  subroutine scal(self, alpha)
    class(state_vector), intent(inout) :: self
    complex(kind=wp)      , intent(in)    :: alpha
    self%state = self%state * alpha
    return
  end subroutine scal

  subroutine axpby(self, alpha, vec, beta)
    class(state_vector)   , intent(inout) :: self
    class(abstract_vector_cdp), intent(in)    :: vec
    complex(kind=wp)         , intent(in)    :: alpha, beta
    select type(vec)
    type is(state_vector)
       self%state = alpha*self%state + beta*vec%state
    end select
    return
  end subroutine axpby

  subroutine rand(self, ifnorm)
    class(state_vector), intent(inout) :: self
    logical, optional,   intent(in)    :: ifnorm
    real(kind=wp) :: tmp(nx, 2)
    ! internals
    logical :: normalize
    real(kind=wp) :: alpha
    normalize = optval(ifnorm,.true.)
    call random_number(tmp)
    self%state%re = tmp(:, 1) ; self%state%im = tmp(:, 2)
    if (normalize) then
      alpha = self%norm()
      call self%scal(cmplx(1.0_wp, 0.0_wp, kind=wp)/alpha)
    endif
    return
  end subroutine rand

  !------------------------------------------------------------------------
  !-----     TYPE-BOUND PROCEDURES FOR THE EXPONENTIAL PROPAGATOR     -----
  !------------------------------------------------------------------------

  subroutine direct_solver(self, vec_in, vec_out)
    ! Linear Operator.
    class(exponential_prop), intent(in)  :: self
    ! Input vector.
    class(abstract_vector_cdp) , intent(in)  :: vec_in
    ! Output vector.
    class(abstract_vector_cdp) , intent(out) :: vec_out

    ! Time-integrator.
    type(rks54_class) :: prop
    real(kind=wp)     :: dt = 1.0_wp
    real(kind=wp)     :: state_ic(2*nx), state_fc(2*nx)

    select type(vec_in)
    type is(state_vector)
       select type(vec_out)
       type is(state_vector)
          ! Get state vector.
          state_ic(:nx) = vec_in%state%re
          state_ic(nx+1:) = vec_in%state%im
          ! Initialize propagator.
          call prop%initialize(n=2*nx, f=rhs)
          ! Integrate forward in time.
          call prop%integrate(0.0_wp, state_ic, dt, self%tau, state_fc)
          ! Pass-back the state vector.
          vec_out%state%re = state_fc(:nx)
          vec_out%state%im = state_fc(nx+1:)
       end select
    end select
    return
  end subroutine direct_solver

  subroutine adjoint_solver(self, vec_in, vec_out)
    ! Linear Operator.
    class(exponential_prop), intent(in)  :: self
    ! Input vector.
    class(abstract_vector_cdp) , intent(in)  :: vec_in
    ! Output vector.
    class(abstract_vector_cdp) , intent(out) :: vec_out

    ! Time-integrator.
    type(rks54_class) :: prop
    real(kind=wp)     :: dt = 1.0_wp
    real(kind=wp)     :: state_ic(2*nx), state_fc(2*nx)

    select type(vec_in)
    type is(state_vector)
       select type(vec_out)
       type is(state_vector)
          ! Get the state.
          state_ic(:nx) = vec_in%state%re
          state_fc(nx+1:) = vec_in%state%im
          ! Initialize propagator.
          call prop%initialize(n=2*nx, f=adjoint_rhs)
          ! Integrate forward in time.
          call prop%integrate(0.0_wp, state_ic, dt, self%tau, state_fc)
          ! Pass-back the state.
          vec_out%state%re = state_fc(:nx)
          vec_out%state%im = state_fc(nx+1:)
       end select
    end select
    return
  end subroutine adjoint_solver

end module Ginzburg_Landau
