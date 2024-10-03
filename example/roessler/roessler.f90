module Roessler
   ! Standard Library.
   use stdlib_stats_distribution_normal, only: normal => rvs_normal
   use stdlib_optval, only: optval
   ! RKLIB module for time integration.
   use rklib_module
   ! LightKrylov for linear algebra.
   use LightKrylov
   use LightKrylov, only: wp => dp
   implicit none
 
   character*128, parameter, private :: this_module = 'Roessler'
 
   public :: a, b, c
 
   !------------------------------
   !-----     PARAMETERS     -----
   !------------------------------
 
   integer,  parameter :: npts = 3
   real(wp), parameter :: a = 0.2
   real(wp), parameter :: b = 0.2
   real(wp), parameter :: c = 5.7
 
   !-------------------------------------------
   !-----     LIGHTKRYLOV VECTOR TYPE     -----
   !-------------------------------------------
 
   type, extends(abstract_vector_rdp), public :: state_vector
      real(wp) :: x = 0.0_wp
      real(wp) :: y = 0.0_wp
      real(wp) :: z = 0.0_wp
      real(wp) :: T = 0.0_wp
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
      !real(wp) :: tau ! Integration time
   contains
      private
      procedure, pass(self), public :: eval => nonlinear_map
   end type roessler_upo

   type, extends(abstract_jacobian_linop_rdp), public :: jacobian
      !real(wp) :: tau ! Integration time.
   contains
      private
      procedure, pass(self), public :: matvec => linear_map
      procedure, pass(self), public :: rmatvec => linear_map ! dummy, we do not need the adjoint of the jacobian
   end type jacobian
 
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
 
   subroutine nonlinear_roessler(x, f)
      ! State vector.
      real(kind=wp)  , dimension(npts), intent(in)  :: x
      ! Time-derivative.
      real(kind=wp)  , dimension(npts), intent(out) :: f
      
      f  = 0.0_wp
      f(1) = -x(2) - x(3)
      f(2) = x(1) + a * x(2)
      f(3) = b + x(3) * (x(1) - c)
      
      return
   end subroutine nonlinear_roessler

   !-----------------------------
   !-----      JACOBIAN     -----
   !-----------------------------

   subroutine linear_roessler(xp, bf, f)
      ! State vector.
      real(kind=wp)  , dimension(:), intent(in)  :: xp
      ! Base state.
      real(kind=wp)  , dimension(:), intent(in)  :: bf
      ! Time-derivative.
      real(kind=wp)  , dimension(:), intent(out) :: f

      f = 0.0_wp
      f(1) = -xp(2) - xp(3)
      f(2) =  xp(1) + a*xp(2)
      f(3) =  xp(1)*bf(3) + xp(3)*(bf(1) - c)
   
      return
   end subroutine linear_roessler

   !-----------------------------------------------------------------
   !-----      WRAPPER FOR NONLINEAR INTEGRATION WITH RKLIB     -----
   !-----------------------------------------------------------------
 
   subroutine NL_rhs(me, t, x, f)
      ! Time-integrator.
      class(rk_class), intent(inout)             :: me
      ! Current time.
      real(kind=wp)  , intent(in)                :: t
      ! State vector.
      real(kind=wp)  , dimension(:), intent(in)  :: x
      ! Time-derivative.
      real(kind=wp)  , dimension(:), intent(out) :: f

      call nonlinear_roessler(x, f)
      
      return
   end subroutine NL_rhs
 
   !-------------------------------------------------------------
   !-----      WRAPPER FOR LINEAR INTEGRATION WITH RKLIB    -----
   !-------------------------------------------------------------
 
   subroutine combined_rhs(me, t, x, f)
      ! Time-integrator.
      class(rk_class), intent(inout)             :: me
      ! Current time.
      real(kind=wp)  , intent(in)                :: t
      ! State vector.
      real(kind=wp)  , dimension(:), intent(in)  :: x
      ! Time-derivative.
      real(kind=wp)  , dimension(:), intent(out) :: f
 
      call nonlinear_roessler(            x(:npts), f(:npts))
      call    linear_roessler(x(npts+1:), x(:npts), f(npts+1:))
     
      return
   end subroutine combined_rhs
  
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
      self%x = 0.0_wp
      self%y = 0.0_wp
      self%z = 0.0_wp
      ! period
      self%T = 0.0_wp
      return
   end subroutine zero
  
   real(wp) function dot(self, vec) result(alpha)
      class(state_vector)       , intent(in) :: self
      class(abstract_vector_rdp), intent(in) :: vec
      select type(vec)
      type is(state_vector)
         alpha = self%x*vec%x + self%y*vec%y + self%z*vec%z + self%T*vec%T
      end select
      return
   end function dot
  
   subroutine scal(self, alpha)
      class(state_vector), intent(inout) :: self
      real(wp)           , intent(in)    :: alpha
      self%x = self%x * alpha
      self%y = self%y * alpha
      self%z = self%z * alpha
      self%T = self%T * alpha
      return
   end subroutine scal
  
   subroutine axpby(self, alpha, vec, beta)
      class(state_vector)       , intent(inout) :: self
      class(abstract_vector_rdp), intent(in)    :: vec
      real(wp)                  , intent(in)    :: alpha, beta
      select type(vec)
      type is(state_vector)
         self%x = alpha*self%x + beta*vec%x
         self%y = alpha*self%y + beta*vec%y
         self%z = alpha*self%z + beta*vec%z
         self%T = alpha*self%T + beta*vec%T
      end select
      return
   end subroutine axpby
  
   integer function get_size(self) result(N)
      class(state_vector), intent(in) :: self
      N = npts+1
      return
   end function get_size
  
   subroutine rand(self, ifnorm)
      class(state_vector), intent(inout) :: self
      logical, optional,   intent(in)    :: ifnorm
      logical :: normalized
      real(wp) :: mu, var
      real(wp) :: alpha
  
      mu = 0.0_wp
      var = 1.0_wp
      self%x = normal(mu, var)
      self%y = normal(mu, var)
      self%z = normal(mu, var)
      self%T = normal(mu, var)

      normalized = optval(ifnorm, .false.)
      if (normalized) then
         alpha = self%norm()
         call self%scal(1.0_wp/alpha)
      endif
      return
   end subroutine rand
 
   !-------------------------------------------------------------
   !-----     TYPE-BOUND PROCEDURES FOR THE INTEGRATORS     -----
   !-------------------------------------------------------------
 
   subroutine nonlinear_map(self, vec_in, vec_out)
      ! Linear Operator.
      class(roessler_upo),        intent(in)  :: self
      ! Input vector.
      class(abstract_vector_rdp), intent(in)  :: vec_in
      ! Output vector.
      class(abstract_vector_rdp), intent(out) :: vec_out
      
      ! Time-integrator.
      type(rks54_class)         :: nonlinear_roessler
      real(wp)                  :: dt = 1.0_wp
      real(wp)                  :: period
      real(wp), dimension(npts) :: pos_in, pos_out
      
      ! Evaluate F(X).
      select type(vec_in)
      type is(state_vector)
         select type(vec_out)
         type is(state_vector)
            ! Get state vector.
            call get_position(vec_in, pos_in)
            ! Initialize integrator.
            call nonlinear_roessler%initialize(n=npts, f=NL_rhs)!, report=roessler_report_stdout, report_rate=1)
            ! Integrate forward in time.
            call nonlinear_roessler%integrate(0.0_wp, pos_in, dt, vec_in%T, pos_out)
            ! Pass-back the state vector.
            call set_position(pos_out, vec_out)

            ! Evaluate residual F(X) - X.
            call vec_out%sub(vec_in)

            ! Add period residual
            vec_out%T = 0.0_wp            
         end select
      end select

      return
   end subroutine nonlinear_map
 
   subroutine linear_map(self, vec_in, vec_out)
      ! Linear Operator.
      class(jacobian), intent(in) :: self
      ! Input vector.
      class(abstract_vector_rdp) , intent(in)  :: vec_in
      ! Output vector.
      class(abstract_vector_rdp) , intent(out) :: vec_out
      
      ! Time-integrator.
      type(rks54_class)           :: combined_roessler
      real(wp)                    :: dt = 1.0_wp
      real(wp)                    :: period
      real(wp), dimension(2*npts) :: pos_in, pos_out
      type(state_vector)          :: vec
      
      select type(vec_in)
      type is(state_vector)
         select type(vec_out)
         type is(state_vector)
            ! Get the state.
            call get_position(self%X, pos_in(:npts))
            call get_period(  self%X, period)
            call get_position(vec_in, pos_in(npts+1:))
            ! Initialize integrator.
            call combined_roessler%initialize(n=2*npts, f=combined_rhs)!, report=roessler_report_stdout, report_rate=1)
            ! Evaluate:
            ! 1. F(X)
            ! 2. exp(tau*J) @ dx
            call combined_roessler%integrate(0.0_wp, pos_in, dt, period, pos_out)
            ! Pass-back the state.
            call set_position(pos_out(npts+1:), vec_out)

            ! Evaluate [ exp(tau*J) - I ] @ dx.
            call vec_out%sub(vec_in)

            ! Evaluate f'(X(T), T) * dT and add it to the position residual
            call compute_fdot(pos_out(:npts), vec)
            call vec_out%axpby(1.0_wp, vec, vec_in%T)

            ! Evaluate f'(X(0), 0).T @ dx and add phase condition
            call compute_fdot(pos_in(:npts), vec) 
            vec_out%T = vec%dot(self%X)
            
         end select
      end select
      
      return
   end subroutine linear_map

   !-------------------------------------------
   !-----     MISCELLANEOUS UTILITIES     -----
   !-------------------------------------------
 
   subroutine get_position(vec_in, pos)
      class(abstract_vector_rdp), intent(in)  :: vec_in
      real(wp), dimension(npts),  intent(out) :: pos

      pos = 0.0_wp
      select type (vec_in)
      type is (state_vector)
          pos(1) = vec_in%x
          pos(2) = vec_in%y
          pos(3) = vec_in%z
      end select

      return
   end subroutine get_position

   subroutine get_period(vec_in, period)
      class(abstract_vector_rdp), intent(in)  :: vec_in
      real(wp),                   intent(out) :: period

      select type (vec_in)
      type is (state_vector)
          period = vec_in%T
      end select

      return
   end subroutine get_period

   subroutine set_position(pos, vec_out)
      real(wp), dimension(npts),  intent(in)  :: pos
      class(abstract_vector_rdp), intent(out) :: vec_out

      select type (vec_out)
      type is (state_vector)
         vec_out%x = pos(1)
         vec_out%y = pos(2)
         vec_out%z = pos(3)
      end select

      return
   end subroutine set_position

   subroutine compute_fdot(pos, vec_out)
      real(wp), dimension(npts) , intent(in)  :: pos
      class(abstract_vector_rdp), intent(out) :: vec_out
      ! internal
      real(wp), dimension(npts) :: wrk
      call nonlinear_roessler(pos, wrk)
      call vec_out%zero()
      call set_position(wrk, vec_out)
      return
   end subroutine compute_fdot

   subroutine roessler_report_stdout(me, t, x)
      class(rk_class), intent(inout)      :: me
      real(wp), intent(in)                :: t
      real(wp), dimension(:), intent(in)  :: x
      
      print '(*(F15.6,1X))', t, x

      return
   end subroutine

   subroutine roessler_report_file(me, t, x)
      class(rk_class), intent(inout)      :: me
      real(wp), intent(in)                :: t
      real(wp), dimension(:), intent(in)  :: x
      ! internals
      integer :: iunit
      open(newunit=iunit, file='new_orbit.txt', status='old', action='write', position='append')
      write(iunit, '(*(F15.6,1X))') t, x
      close(iunit)
      return
   end subroutine
 
end module Roessler