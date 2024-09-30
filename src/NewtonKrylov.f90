module LightKrylov_NewtonKrylov
   use stdlib_optval, only: optval
   use LightKrylov_Constants
   use LightKrylov_Logger
   use LightKrylov_AbstractVectors
   use LightKrylov_AbstractLinops
   use LightKrylov_AbstractSystems
   use LightKrylov_IterativeSolvers
   use LightKrylov_Utils

   implicit none
   private

   character*128, parameter :: this_module = 'LightKrylov_NewtonKrylov'

   public :: newton
   public :: constant_atol, dynamic_tol

   interface newton
      module procedure newton_rdp
   end interface

   interface constant_atol
      module procedure constant_atol_rdp
   end interface

   interface dynamic_tol
      module procedure dynamic_tol_rdp
   end interface

   abstract interface 
      subroutine abstract_scheduler_rdp(tol, rnorm, iter, info)
         import dp
         !! Abstract interface to define tolerance scheduler for the Newton iteration
         real(dp), intent(out) :: tol
         !! Tolerance to be used
         real(dp), intent(in)  :: rnorm
         !! Norm of the residual of the current iterate
         integer,  intent(in)  :: iter
         !! Newton iteration count
         integer,  intent(out)  :: info
         !! Information flag
      end subroutine abstract_scheduler_rdp
   end interface

   type, extends(abstract_opts), public :: newton_dp_opts
      !! Options for Newton-Krylov fixed-point iteration.
      integer :: maxiter = 100
      !! Maximum number of Newton iterations (default = 100)
      logical :: ifbisect = .false.
      !! Bisection toggle to enforce residual reduction (default = .false.)
      integer :: maxstep_bisection = 5
      !! Maximum number of bisections (evaluations of F) for step selection (default = 5)
      !! Ignored if ifbisect = .false.
      logical :: verbose = .false.
      !! Verbosity control (default: `.false.`)
      class(gmres_dp_opts), allocatable :: gmres_opts
      !! GMRES solver options
    end type

contains

   subroutine newton_rdp(sys, X, info, options, scheduler)
      !! Classic no-frills implementation of the Newton-Krylov root-finding algorithm
      class(abstract_system_rdp), intent(inout)  :: sys
      !! Dynamical system for which we wish to compute a fixed point
      class(abstract_vector_rdp), intent(inout)  :: X
      !! Initial guess for the fixed point, will be overwritten with solution
      integer,                    intent(out)    :: info
      !! Information flag
      type(newton_dp_opts),  optional, intent(in)   :: options
      type(newton_dp_opts)                          :: opts
      !! Options for the Newton-Krylov iteration

      procedure(abstract_scheduler_rdp), optional :: scheduler

      !--------------------------------------
      !-----     Internal variables     -----
      !--------------------------------------
      
      ! residual vector
      class(abstract_vector_rdp), allocatable :: residual, increment
      real(dp) :: rnorm, tol
      logical :: converged, verb
      integer :: i, maxiter, maxstep_bisection
      type(gmres_dp_opts) :: gmres_opts
      procedure(abstract_scheduler_rdp), pointer :: tolerance_scheduler => null()

      ! Options
      if (present(options)) then
         opts       = options
         gmres_opts = options%gmres_opts
      else ! default
         opts       = newton_dp_opts()
         gmres_opts = gmres_dp_opts()

      end if
      if (present(scheduler)) then
         tolerance_scheduler => scheduler
      else
         tolerance_scheduler => constant_atol_rdp
      endif

      ! Initialisation
      verb = opts%verbose
      
      maxiter = opts%maxiter
      maxstep_bisection = opts%maxstep_bisection
      converged = .false.
      allocate(residual, source=X); call residual%zero()
      allocate(increment,source=X); call increment%zero()

      if (verb) write(*,*) 'Starting Newton ...'
      ! Newton iteration
      newton: do i = 1, maxiter

         call sys%eval(X, residual)
         rnorm = residual%norm()
         if (verb) write(*,*) "Iteration", i, ": Residual norm = ", rnorm

         ! Check for convergence.
         call tolerance_scheduler(tol, rnorm, i, info)
         if (rnorm < tol) then
            if (verb) write(*,*) "Convergence achieved."
            converged = .true.
            exit newton
         end if

         ! Define the Jacobian
         call sys%set_base_state(X)
        
         ! Solve the linear system using GMRES.
         call residual%chsgn()
         call gmres(sys%jacobian, residual, increment, info)
         call check_info(info, 'gmres', module=this_module, procedure='newton_classic_rdp')

         ! Update the solution and overwrite X0
         if (opts%ifbisect) then
            call increment_bisection(X, sys, increment, rnorm, maxstep_bisection, verb)
         else
            call X%add(increment)
         endif

      enddo newton

      if (.not.converged) then
         if (verb) write(*,*) 'Newton iteration did not converge within', maxiter, 'steps.'
         info = -1
      endif

      return
   end subroutine newton_rdp

   subroutine increment_bisection(X, sys, increment, rnorm_old, maxstep, verb)
      class(abstract_vector_rdp), intent(inout) :: X
      class(abstract_system_rdp), intent(in)    :: sys
      class(abstract_vector_rdp), intent(in)    :: increment
      real(dp),                   intent(in)    :: rnorm_old
      integer,                    intent(in)    :: maxstep
      logical,                    intent(in)    :: verb
      ! internals
      integer :: i
      real(dp) :: alpha, step, rnorm_ref, rnorm0, rnorm1
      class(abstract_vector_rdp), allocatable :: Xin, residual
      character(len=128) :: fmt

      write(fmt,*) '(A,I3,A,F8.6,A,E9.3)'

      allocate(Xin, source=X)
      allocate(residual, source=X); call residual%zero()

      if (verb) then
         print *, 'Current residual norm:', rnorm_old
         print *, 'Start Newton step bisection ...'
      end if

      alpha = 1.0_dp
      call X%axpby(one_rdp, increment, alpha)
      ! evaluate residual norm
      call sys%eval(X, residual)
      rnorm1 = residual%norm()    
      if (verb) print fmt, '   step', 1, ': alpha = ', alpha, ' => rnorm = ', rnorm1
      rnorm_ref = rnorm1 ! Save norm of full step as reference

      alpha = alpha/2
      call X%axpby(one_rdp, increment, alpha)
      ! evaluate residual norm
      call sys%eval(X, residual)
      rnorm0 = residual%norm()    
      if (verb) print fmt, '   step ', 2, ': alpha = ', alpha, ' => rnorm = ', rnorm0

      step = 0.5_dp
      do i = 3, maxstep
         ! decide on bisection step
         step  = step/2
         if (rnorm1 < rnorm0) then
            alpha = alpha + step
         else
            alpha = alpha - step
            rnorm0 = rnorm1
         end if
         ! compute new trial solution
         call X%axpby(zero_rdp, Xin, one_rdp)
         call X%axpby(one_rdp, increment, alpha)
         ! evaluate residual norm
         call sys%eval(X, residual)
         rnorm1 = residual%norm()
         if (verb) print fmt, '   step', i, ': alpha = ', alpha, ' => rnorm = ', rnorm1
      enddo
      if (rnorm1 > rnorm_ref) then
         print *, 'Newton step bisection: Reverting to full Newton step.'
         ! compute new trial solution
         call X%axpby(zero_rdp, Xin, one_rdp)
         call X%axpby(one_rdp, increment, one_rdp)
      end if

      return
   end subroutine

   subroutine constant_atol_rdp(tol, rnorm, iter, info)
      !! Abstract interface to define tolerance scheduler for the Newton iteration
      real(dp), intent(out) :: tol
      !! Tolerance to be used
      real(dp), intent(in)  :: rnorm
      !! Norm of the residual of the current iterate
      integer,  intent(in)  :: iter
      !! Newton iteration count
      integer,  intent(out)  :: info
      !! Information flag
      tol = atol_dp
      return
   end subroutine constant_atol_rdp

   subroutine dynamic_tol_rdp(tol, rnorm, iter, info)
      !! Abstract interface to define tolerance scheduler for the Newton iteration
      real(dp), intent(out) :: tol
      !! Tolerance to be used
      real(dp), intent(in)  :: rnorm
      !! Norm of the residual of the current iterate
      integer,  intent(in)  :: iter
      !! Newton iteration count
      integer,  intent(out)  :: info
      !! Information flag
      tol = max(0.1*rnorm, atol_dp)
      return
   end subroutine dynamic_tol_rdp

end module LightKrylov_NewtonKrylov
