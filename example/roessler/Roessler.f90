module Roessler
   ! RKLIB module for time integration.
   use rklib_module
   ! LightKrylov for linear algebra.
   use LightKrylov
   use LightKrylov, only: wp => dp
   ! Standard Library.
   use stdlib_math, only : linspace
   use stdlib_optval, only : optval
   implicit none
 
   character*128, parameter, private :: this_module = 'Roessler'
 
   public :: a, b, c
 
   !------------------------------
   !-----     PARAMETERS     -----
   !------------------------------
 
   real(dp), parameter :: a = 0.2
   real(dp), parameter :: b = 0.2
   real(dp), parameter :: c = 5.7
 
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
      procedure, pass(self), public :: get_size
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
 
   !===================================
   !===================================
   !=====                         =====
   !=====     ROESSLER SYSTEM     =====
   !=====                         =====
   !===================================
   !===================================
 
   !--------------------------------------------------------------
   !-----     CONSTRUCT THE MESH AND PHYSICAL PARAMETERS     -----
   !--------------------------------------------------------------
  
 
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
      self%x = 0.0_dp
      self%y = 0.0_dp
      self%z = 0.0_dp
      return
   end subroutine zero
  
   real(dp) function dot(self, vec) result(alpha)
      class(state_vector)       , intent(in) :: self
      class(abstract_vector_rdp), intent(in) :: vec
      select type(vec)
      type is(state_vector)
         alpha = self%x*vec%x + self%y*vec%y + self%z*vec%z
      end select
      return
   end function dot
  
   subroutine scal(self, alpha)
      class(state_vector), intent(inout) :: self
      real(dp)           , intent(in)    :: alpha
      self%x = self%x * alpha
      self%y = self%y * alpha
      self%z = self%z * alpha
      return
   end subroutine scal
  
   subroutine axpby(self, alpha, vec, beta)
      class(state_vector)       , intent(inout) :: self
      class(abstract_vector_rdp), intent(in)    :: vec
      real(dp)                  , intent(in)    :: alpha, beta
      select type(vec)
      type is(state_vector)
         self%x = alpha*self%x + beta*vec%x
         self%y = alpha*self%y + beta*vec%y
         self%z = alpha*self%z + beta*vec%z
      end select
      return
   end subroutine axpby
  
   integer function get_size(self) result(N)
      class(state_vector), intent(in) :: self
      N = 3
      return
   end function get_size
  
   subroutine rand(self, ifnorm)
      class(state_vector), intent(inout) :: self
      logical, optional,   intent(in)    :: ifnorm
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
 
 end module Roessler