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
 
   integer,  parameter :: npts = 3
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

   !-------------------------------------------
   !-----     LIGHTKRYLOV SYSTEM TYPE     -----
   !-------------------------------------------

   type, extends(abstract_system_rdp), public :: roessler
   contains
      private
      procedure, pass(self), public :: eval => nonlinear_map
   end type roessler

   !type, extends(abstract_jacobian_linop_rdp), public :: jacobian
   !   real(kind=wp) :: tau ! Integration time.
   !contains
   !   private
   !   procedure, pass(self), public :: matvec => linear_map
   !   procedure, pass(self), public :: rmatvec => adj_linear_map
   !end type jacobian
 
 contains
 
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
 
   subroutine NLrhs(me, t, x, f)
      ! Time-integrator.
      class(rk_class), intent(inout)             :: me
      ! Current time.
      real(kind=wp)  , intent(in)                :: t
      ! State vector.
      real(kind=wp)  , dimension(:), intent(in)  :: x
      ! Time-derivative.
      real(kind=wp)  , dimension(:), intent(out) :: f

      f  = 0.0_wp
      f(1) = -x(2) - x(3)
      f(2) = x(1) + a * x(2)
      f(3) = b + x(3) * (x(1) - c)
      
      return
    end subroutine NLrhs
 
   !-----------------------------
   !-----      JACOBIAN     -----
   !-----------------------------
 
   subroutine Lrhs(me, t, x, f)
      ! Time-integrator.
      class(rk_class), intent(inout)             :: me
      ! Current time.
      real(kind=wp)  , intent(in)                :: t
      ! State vector.
      real(kind=wp)  , dimension(:), intent(in)  :: x
      ! Time-derivative.
      real(kind=wp)  , dimension(:), intent(out) :: f

      ! internals
      real(kind=wp), dimension(npts) :: X_bf, xp

      ! extract nonlinear and linear states
      X_bf = x(1:npts)
      xp   = x(npts+1:2*npts)
 
      f = 0.0_wp
      f(1) = -xp(2) - xp(3)
      f(2) =  xp(1) + a*xp(2)
      f(3) =  xp(1)*X_bf(3) + xp(3)*(X_bf(1) - c)
     
     return
   end subroutine Lrhs
  
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
 
   !-------------------------------------------------------------
   !-----     TYPE-BOUND PROCEDURES FOR THE INTEGRATORS     -----
   !-------------------------------------------------------------
 
   subroutine nonlinear_map(self, vec_in, vec_out)
      ! Linear Operator.
      class(roessler), intent(in)  :: self
      ! Input vector.
      class(abstract_vector_cdp) , intent(in)  :: vec_in
      ! Output vector.
      class(abstract_vector_cdp) , intent(out) :: vec_out
      
      ! Time-integrator.
      type(rks54_class) :: nonlinear_roessler
      real(kind=wp)     :: dt = 1.0_wp
      real(kind=wp)     :: state_in(npts), state_out(npts)
      
      select type(vec_in)
      type is(state_vector)
         select type(vec_out)
         type is(state_vector)
            ! Get state vector.
            call get_state_vector(vec_in, state_in)
            ! Initialize propagator.
            call nonlinear_roessler%initialize(n=npts, f=NLrhs)
            ! Integrate forward in time.
            call nonlinear_roessler%integrate(0.0_wp, state_in, dt, self%tau, state_out)
            ! Pass-back the state vector.
            call set_state_vector(state_out, vec_out)
         end select
      end select
      return
   end subroutine nonlinear_map
 
!   subroutine linear_map(self, vec_in, vec_out)
!      ! Linear Operator.
!      class(exponential_prop), intent(in)  :: self
!      ! Input vector.
!      class(abstract_vector_cdp) , intent(in)  :: vec_in
!      ! Output vector.
!      class(abstract_vector_cdp) , intent(out) :: vec_out
!      
!      ! Time-integrator.
!      type(rks54_class) :: prop
!      real(kind=wp)     :: dt = 1.0_wp
!      real(kind=wp)     :: state_ic(2*nx), state_fc(2*nx)
!      
!      select type(vec_in)
!      type is(state_vector)
!         select type(vec_out)
!         type is(state_vector)
!            ! Get the state.
!            state_ic(:nx) = vec_in%state%re
!            state_fc(nx+1:) = vec_in%state%im
!            ! Initialize propagator.
!            call prop%initialize(n=2*nx, f=adjoint_rhs)
!            ! Integrate forward in time.
!            call prop%integrate(0.0_wp, state_ic, dt, self%tau, state_fc)
!            ! Pass-back the state.
!            vec_out%state%re = state_fc(:nx)
!            vec_out%state%im = state_fc(nx+1:)
!         end select
!      end select
!      return
!   end subroutine linear_map

   !-------------------------------------------
   !-----     MISCELLANEOUS UTILITIES     -----
   !-------------------------------------------
 

   subroutine get_state_vector(vec_in, state)
      class(abstract_vector),    intent(in)  :: vec_in
      real(dp), dimension(npts), intent(out) :: state

      state = 0.0_wp
      select type (vec_in)
      type is (state_vector)
          state(1) = vec_in%x
          state(2) = vec_in%y
          state(3) = vec_in%z        
      end select

      return
   end subroutine get_state_vector

   subroutine set_state_vector(state, vec_out)
      real(dp), dimension(npts), intent(in)  :: state
      class(abstract_vector),    intent(out) :: vec_out

      select type (vec_out)
      type is (state_vector)
         vec_out%x = state(1)
         vec_out%y = state(2)
         vec_out%z = state(3)        
      end select

      return
   end subroutine get_state_vector
 
end module Roessler