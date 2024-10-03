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
      !! Implements the simple Newton-Krylov method for finding roots (fixed points) of a nonlinear vector-valued function
      !! \( F(\mathbf{X}) \), i.e. solutions \( \mathbf{X}^* \) such that \( F(\mathbf{X}^*) - \mathbf{X}^* = \mathbf{0} \) 
      !! starting from an initial guess via successive solution increments based on local linearization \( \mathbf{J}_\mathbf{X} \) 
      !! (the Jacobian) of the nonlinear function in the vicinity of the current solution.
      !!
      !! **Algorthmic Features**
      !! - At each iteration, the next Newton step \( \delta\mathbf{x}\) is computed as the solution of the linear system
      !!   
      !! $$ \mathbf{J}_\mathbf{X} \delta \mathbf{x} = \mathbf{r} $$
      !!
      !!   where \( \mathbf{r} = -F(\mathbf{X}) \) is the residual of the nonlinear function.
      !! - The Jacobian is never assembled and the linear system is solved using one of the available iterative solvers.
      !! - When the residual norm does not decrease during iteration indicating that the linearization is not a very
      !!   accurate model of the function's behaviour, which often happens during the initial iterations, a 1D step bisection
      !!   method based on the golden ratio is implemented to dampen the step and improve convergence of the method.
      !! - The implementation allows for dynamic tolerances (also known as inexact Newton), where the approximation for 
      !!   the residual and the linear system can be solved with relaxed tolerances to reduce overall time to solution.
      !! - The method is suitable to both fixed points and periodic orbits via the choice of residual and corresponding
      !!   Jacobian matrix. In the case of unforced periodic orbits, the period is itself an unknown that must be included
      !!   in the iteration.
      !!
      !! **Advantages**
      !! - The iterative solution of the linear systems has a comparatively low storage footprint.
      !! - If the Newton iteration converges, the convergence is formally asymptotically of second order. Using dynamic
      !!   tolerances and line searches slightly reduce this convergence rate in exchange for a larger convergence region.
      !! 
      !!
      !! **Limitations**
      !! - The method is not guaranteed to converge if the initial guess is too far from the fixed point. 
      !!   If the Newton iteration diverges even with step bisection, the best suggestion is to find a 
      !!   better initial guess. If this is not feasible, some alternatives to improve the convergence 
      !!   of the Newton iteration are possible (but not implemented to date), including various line search
      !!   algorithms and trust region methods (doglog, double dogleg, hookstep, ...).
      !!
      !! **References**
      !! - Viswanath, D. (2007). "Recurrent motions within plane Couette turbulence". Journal of Fluid Mechanics, 580, 339-358.
      !! - Sánchez, J., Net, M., Garcıa-Archilla, B., & Simó, C. (2004). "Newton–Krylov continuation of periodic orbits 
      !!   for Navier–Stokes flows". Journal of Computational Physics, 201(1), 13-33.
      !! - Frantz, R. A., Loiseau, J. C., & Robinet, J. C. (2023). "Krylov methods for large-scale dynamical systems: Application 
      !!   in fluid dynamics". Applied Mechanics Reviews, 75(3), 030802.
      !!
      module procedure newton_rsp
      module procedure newton_rdp
      module procedure newton_csp
      module procedure newton_cdp
   end interface

   interface constant_atol
      module procedure constant_atol_sp
      module procedure constant_atol_dp
   end interface

   interface dynamic_tol
      module procedure dynamic_tol_sp
      module procedure dynamic_tol_dp
   end interface

   abstract interface
      subroutine abstract_scheduler_sp(tol, rnorm, iter, info)
         import sp
         !! Abstract interface to define a tolerance scheduler for the Newton iteration
         real(sp), intent(out) :: tol
         !! Tolerance to be used
         real(sp), intent(in)  :: rnorm
         !! Norm of the residual of the current iterate
         integer,  intent(in)  :: iter
         !! Newton iteration count
         integer,  intent(out)  :: info
         !! Information flag
      end subroutine abstract_scheduler_sp
      subroutine abstract_scheduler_dp(tol, rnorm, iter, info)
         import dp
         !! Abstract interface to define a tolerance scheduler for the Newton iteration
         real(dp), intent(out) :: tol
         !! Tolerance to be used
         real(dp), intent(in)  :: rnorm
         !! Norm of the residual of the current iterate
         integer,  intent(in)  :: iter
         !! Newton iteration count
         integer,  intent(out)  :: info
         !! Information flag
      end subroutine abstract_scheduler_dp
   end interface

   type, extends(abstract_opts), public :: newton_opts
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
    end type

contains

   subroutine newton_rsp(sys, X, info, options, linear_solver, linear_solver_options, preconditioner, scheduler)
      class(abstract_system_rsp),                         intent(inout) :: sys
      !! Dynamical system for which we wish to compute a fixed point
      class(abstract_vector_rsp),                         intent(inout) :: X
      !! Initial guess for the fixed point, will be overwritten with solution
      integer,                                            intent(out)   :: info
      !! Information flag
      type(newton_opts),                        optional, intent(in)    :: options
      type(newton_opts)                                                 :: opts
      !! Options for the Newton-Krylov iteration
      procedure(abstract_linear_solver_rsp),    optional                :: linear_solver
      !! Linear solver to be used to find Newton step
      class(abstract_opts),                     optional, intent(in)    :: linear_solver_options
      class(abstract_opts), allocatable                                 :: solver_opts
      !! Options for the linear solver
      class(abstract_precond_rsp),              optional, intent(in)    :: preconditioner
      class(abstract_precond_rsp), allocatable                          :: precond
      !! Preconditioner for the linear solver
      procedure(abstract_scheduler_sp),        optional                :: scheduler

      !--------------------------------------
      !-----     Internal variables     -----
      !--------------------------------------
      
      ! residual vector
      class(abstract_vector_rsp), allocatable :: residual, increment
      real(sp) :: rnorm, tol
      logical :: converged, has_precond, has_solver_opts, verb
      integer :: i, maxiter, maxstep_bisection
      procedure(abstract_linear_solver_rsp), pointer :: solver => null()
      procedure(abstract_scheduler_sp), pointer :: tolerance_scheduler => null()

      ! Newton-Krylov options
      if (present(options)) then
         opts = options
      else
         opts = newton_opts()
      end if
      ! Linear solver
      if (present(linear_solver)) then
         solver => linear_solver
      else
         solver => gmres_rsp
      end if
      ! Linear solver options ?
      if (present(linear_solver_options)) then
         has_solver_opts = .true.
         allocate(solver_opts, source=linear_solver_options)
      else
         has_solver_opts = .false.
      end if
      ! Preconditioner ?
      if (present(preconditioner)) then
         has_precond = .true.
         allocate(precond, source=preconditioner)
      else
         has_precond = .false.
      end if
      ! Scheduler
      if (present(scheduler)) then
         tolerance_scheduler => scheduler
      else
         tolerance_scheduler => constant_atol_sp
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
         sys%jacobian%X = X
        
         ! Solve the linear system using GMRES.
         call residual%chsgn()
         if (.not. has_precond .and. .not. has_solver_opts) then
            call solver(sys%jacobian, residual, increment, info,                                              transpose=.false.)
         elseif (.not. has_precond .and. has_solver_opts) then
            call solver(sys%jacobian, residual, increment, info,                         options=solver_opts, transpose=.false.)
         elseif (has_precond .and. .not. has_solver_opts) then
            call solver(sys%jacobian, residual, increment, info, preconditioner=precond,                      transpose=.false.)
         else
            call solver(sys%jacobian, residual, increment, info, preconditioner=precond, options=solver_opts, transpose=.false.)
         end if
         print *, 'DEBUG: Newton-Krylov iteration', i, ': solver info:', info
         !call check_info(info, 'linear_solver', module=this_module, procedure='newton_rsp')

         ! Update the solution and overwrite X0
         if (opts%ifbisect) then
            call increment_bisection_rsp(X, sys, increment, rnorm, maxstep_bisection, verb)
         else
            call X%add(increment)
         endif

      enddo newton

      if (.not.converged) then
         if (verb) write(*,*) 'Newton iteration did not converge within', maxiter, 'steps.'
         info = -1
      endif

      return
   end subroutine newton_rsp

   subroutine newton_rdp(sys, X, info, options, linear_solver, linear_solver_options, preconditioner, scheduler)
      class(abstract_system_rdp),                         intent(inout) :: sys
      !! Dynamical system for which we wish to compute a fixed point
      class(abstract_vector_rdp),                         intent(inout) :: X
      !! Initial guess for the fixed point, will be overwritten with solution
      integer,                                            intent(out)   :: info
      !! Information flag
      type(newton_opts),                        optional, intent(in)    :: options
      type(newton_opts)                                                 :: opts
      !! Options for the Newton-Krylov iteration
      procedure(abstract_linear_solver_rdp),    optional                :: linear_solver
      !! Linear solver to be used to find Newton step
      class(abstract_opts),                     optional, intent(in)    :: linear_solver_options
      class(abstract_opts), allocatable                                 :: solver_opts
      !! Options for the linear solver
      class(abstract_precond_rdp),              optional, intent(in)    :: preconditioner
      class(abstract_precond_rdp), allocatable                          :: precond
      !! Preconditioner for the linear solver
      procedure(abstract_scheduler_dp),        optional                :: scheduler

      !--------------------------------------
      !-----     Internal variables     -----
      !--------------------------------------
      
      ! residual vector
      class(abstract_vector_rdp), allocatable :: residual, increment
      real(dp) :: rnorm, tol
      logical :: converged, has_precond, has_solver_opts, verb
      integer :: i, maxiter, maxstep_bisection
      procedure(abstract_linear_solver_rdp), pointer :: solver => null()
      procedure(abstract_scheduler_dp), pointer :: tolerance_scheduler => null()

      ! Newton-Krylov options
      if (present(options)) then
         opts = options
      else
         opts = newton_opts()
      end if
      ! Linear solver
      if (present(linear_solver)) then
         solver => linear_solver
      else
         solver => gmres_rdp
      end if
      ! Linear solver options ?
      if (present(linear_solver_options)) then
         has_solver_opts = .true.
         allocate(solver_opts, source=linear_solver_options)
      else
         has_solver_opts = .false.
      end if
      ! Preconditioner ?
      if (present(preconditioner)) then
         has_precond = .true.
         allocate(precond, source=preconditioner)
      else
         has_precond = .false.
      end if
      ! Scheduler
      if (present(scheduler)) then
         tolerance_scheduler => scheduler
      else
         tolerance_scheduler => constant_atol_dp
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
         sys%jacobian%X = X
        
         ! Solve the linear system using GMRES.
         call residual%chsgn()
         if (.not. has_precond .and. .not. has_solver_opts) then
            call solver(sys%jacobian, residual, increment, info,                                              transpose=.false.)
         elseif (.not. has_precond .and. has_solver_opts) then
            call solver(sys%jacobian, residual, increment, info,                         options=solver_opts, transpose=.false.)
         elseif (has_precond .and. .not. has_solver_opts) then
            call solver(sys%jacobian, residual, increment, info, preconditioner=precond,                      transpose=.false.)
         else
            call solver(sys%jacobian, residual, increment, info, preconditioner=precond, options=solver_opts, transpose=.false.)
         end if
         print *, 'DEBUG: Newton-Krylov iteration', i, ': solver info:', info
         !call check_info(info, 'linear_solver', module=this_module, procedure='newton_rdp')

         ! Update the solution and overwrite X0
         if (opts%ifbisect) then
            call increment_bisection_rdp(X, sys, increment, rnorm, maxstep_bisection, verb)
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

   subroutine newton_csp(sys, X, info, options, linear_solver, linear_solver_options, preconditioner, scheduler)
      class(abstract_system_csp),                         intent(inout) :: sys
      !! Dynamical system for which we wish to compute a fixed point
      class(abstract_vector_csp),                         intent(inout) :: X
      !! Initial guess for the fixed point, will be overwritten with solution
      integer,                                            intent(out)   :: info
      !! Information flag
      type(newton_opts),                        optional, intent(in)    :: options
      type(newton_opts)                                                 :: opts
      !! Options for the Newton-Krylov iteration
      procedure(abstract_linear_solver_csp),    optional                :: linear_solver
      !! Linear solver to be used to find Newton step
      class(abstract_opts),                     optional, intent(in)    :: linear_solver_options
      class(abstract_opts), allocatable                                 :: solver_opts
      !! Options for the linear solver
      class(abstract_precond_csp),              optional, intent(in)    :: preconditioner
      class(abstract_precond_csp), allocatable                          :: precond
      !! Preconditioner for the linear solver
      procedure(abstract_scheduler_sp),        optional                :: scheduler

      !--------------------------------------
      !-----     Internal variables     -----
      !--------------------------------------
      
      ! residual vector
      class(abstract_vector_csp), allocatable :: residual, increment
      real(sp) :: rnorm, tol
      logical :: converged, has_precond, has_solver_opts, verb
      integer :: i, maxiter, maxstep_bisection
      procedure(abstract_linear_solver_csp), pointer :: solver => null()
      procedure(abstract_scheduler_sp), pointer :: tolerance_scheduler => null()

      ! Newton-Krylov options
      if (present(options)) then
         opts = options
      else
         opts = newton_opts()
      end if
      ! Linear solver
      if (present(linear_solver)) then
         solver => linear_solver
      else
         solver => gmres_csp
      end if
      ! Linear solver options ?
      if (present(linear_solver_options)) then
         has_solver_opts = .true.
         allocate(solver_opts, source=linear_solver_options)
      else
         has_solver_opts = .false.
      end if
      ! Preconditioner ?
      if (present(preconditioner)) then
         has_precond = .true.
         allocate(precond, source=preconditioner)
      else
         has_precond = .false.
      end if
      ! Scheduler
      if (present(scheduler)) then
         tolerance_scheduler => scheduler
      else
         tolerance_scheduler => constant_atol_sp
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
         sys%jacobian%X = X
        
         ! Solve the linear system using GMRES.
         call residual%chsgn()
         if (.not. has_precond .and. .not. has_solver_opts) then
            call solver(sys%jacobian, residual, increment, info,                                              transpose=.false.)
         elseif (.not. has_precond .and. has_solver_opts) then
            call solver(sys%jacobian, residual, increment, info,                         options=solver_opts, transpose=.false.)
         elseif (has_precond .and. .not. has_solver_opts) then
            call solver(sys%jacobian, residual, increment, info, preconditioner=precond,                      transpose=.false.)
         else
            call solver(sys%jacobian, residual, increment, info, preconditioner=precond, options=solver_opts, transpose=.false.)
         end if
         print *, 'DEBUG: Newton-Krylov iteration', i, ': solver info:', info
         !call check_info(info, 'linear_solver', module=this_module, procedure='newton_csp')

         ! Update the solution and overwrite X0
         if (opts%ifbisect) then
            call increment_bisection_csp(X, sys, increment, rnorm, maxstep_bisection, verb)
         else
            call X%add(increment)
         endif

      enddo newton

      if (.not.converged) then
         if (verb) write(*,*) 'Newton iteration did not converge within', maxiter, 'steps.'
         info = -1
      endif

      return
   end subroutine newton_csp

   subroutine newton_cdp(sys, X, info, options, linear_solver, linear_solver_options, preconditioner, scheduler)
      class(abstract_system_cdp),                         intent(inout) :: sys
      !! Dynamical system for which we wish to compute a fixed point
      class(abstract_vector_cdp),                         intent(inout) :: X
      !! Initial guess for the fixed point, will be overwritten with solution
      integer,                                            intent(out)   :: info
      !! Information flag
      type(newton_opts),                        optional, intent(in)    :: options
      type(newton_opts)                                                 :: opts
      !! Options for the Newton-Krylov iteration
      procedure(abstract_linear_solver_cdp),    optional                :: linear_solver
      !! Linear solver to be used to find Newton step
      class(abstract_opts),                     optional, intent(in)    :: linear_solver_options
      class(abstract_opts), allocatable                                 :: solver_opts
      !! Options for the linear solver
      class(abstract_precond_cdp),              optional, intent(in)    :: preconditioner
      class(abstract_precond_cdp), allocatable                          :: precond
      !! Preconditioner for the linear solver
      procedure(abstract_scheduler_dp),        optional                :: scheduler

      !--------------------------------------
      !-----     Internal variables     -----
      !--------------------------------------
      
      ! residual vector
      class(abstract_vector_cdp), allocatable :: residual, increment
      real(dp) :: rnorm, tol
      logical :: converged, has_precond, has_solver_opts, verb
      integer :: i, maxiter, maxstep_bisection
      procedure(abstract_linear_solver_cdp), pointer :: solver => null()
      procedure(abstract_scheduler_dp), pointer :: tolerance_scheduler => null()

      ! Newton-Krylov options
      if (present(options)) then
         opts = options
      else
         opts = newton_opts()
      end if
      ! Linear solver
      if (present(linear_solver)) then
         solver => linear_solver
      else
         solver => gmres_cdp
      end if
      ! Linear solver options ?
      if (present(linear_solver_options)) then
         has_solver_opts = .true.
         allocate(solver_opts, source=linear_solver_options)
      else
         has_solver_opts = .false.
      end if
      ! Preconditioner ?
      if (present(preconditioner)) then
         has_precond = .true.
         allocate(precond, source=preconditioner)
      else
         has_precond = .false.
      end if
      ! Scheduler
      if (present(scheduler)) then
         tolerance_scheduler => scheduler
      else
         tolerance_scheduler => constant_atol_dp
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
         sys%jacobian%X = X
        
         ! Solve the linear system using GMRES.
         call residual%chsgn()
         if (.not. has_precond .and. .not. has_solver_opts) then
            call solver(sys%jacobian, residual, increment, info,                                              transpose=.false.)
         elseif (.not. has_precond .and. has_solver_opts) then
            call solver(sys%jacobian, residual, increment, info,                         options=solver_opts, transpose=.false.)
         elseif (has_precond .and. .not. has_solver_opts) then
            call solver(sys%jacobian, residual, increment, info, preconditioner=precond,                      transpose=.false.)
         else
            call solver(sys%jacobian, residual, increment, info, preconditioner=precond, options=solver_opts, transpose=.false.)
         end if
         print *, 'DEBUG: Newton-Krylov iteration', i, ': solver info:', info
         !call check_info(info, 'linear_solver', module=this_module, procedure='newton_cdp')

         ! Update the solution and overwrite X0
         if (opts%ifbisect) then
            call increment_bisection_cdp(X, sys, increment, rnorm, maxstep_bisection, verb)
         else
            call X%add(increment)
         endif

      enddo newton

      if (.not.converged) then
         if (verb) write(*,*) 'Newton iteration did not converge within', maxiter, 'steps.'
         info = -1
      endif

      return
   end subroutine newton_cdp


   subroutine increment_bisection_rsp(X, sys, increment, rold, maxstep, verb)
      !! Classic 1D bisection method based on the golden ratio to damped the Newton step in 
      !! order to maximally reduce the residual at each iteration.
      class(abstract_vector_rsp), intent(inout) :: X
      !! Current system state to be updated
      class(abstract_system_rsp), intent(in)    :: sys
      !! Dynamical system for which the residual is minimized
      class(abstract_vector_rsp), intent(in)    :: increment
      !! Newton step computed from the standard method
      real(sp),                   intent(in)    :: rold
      !! Residual of the current system state to determine improvement
      integer,                    intent(in)    :: maxstep
      !! Maximum number of bisection steps. Each additional bisection step requires an evaluation of the nonlinear function
      logical,                    intent(in)    :: verb
      !! Verbosity toggle
      ! internals
      integer :: i, j, idx(1)
      real(sp) :: invphi, invphi2
      real(sp) :: alpha(4), step
      real(sp) :: res(4)
      class(abstract_vector_rsp), allocatable :: Xin, residual

      allocate(Xin, source=X)
      allocate(residual, source=X); call residual%zero()
      step    = one_rsp
      invphi  = (sqrt(5.0) - 1.0)/2.0  ! 1 / phi
      invphi2 = (3.0 - sqrt(5.0))/2.0  ! 1 / phi**2
      alpha = (/ zero_rsp, invphi2*one_rsp, invphi*one_rsp, one_rsp /)
      res   = (/ rold, zero_rsp, zero_rsp, zero_rsp /)

      call X%add(increment)
      ! evaluate residual norm
      call sys%eval(X, residual)
      res(4) = residual%norm()

      if (res(4) > rold) then
         if (verb) print *, 'Start Newton step bisection ...'
         ! compute new trial solutions
         do j = 2, 3
            call X%axpby(zero_rsp, Xin, one_rsp)
            call X%axpby(one_rsp, increment, alpha(j))
            call sys%eval(X, residual)
            res(j) = residual%norm()
         end do

         do i = 1, maxstep
            step = step * invphi
            if (res(2) < res(3)) then
               ! alphas
               ! a1 is kept
               alpha(3:4) = alpha(2:3)           
               ! residuals
               ! r1 is kept
               res(3:4) = res(2:3)
               ! new point --> a2, r2
               alpha(2) = alpha(1) + step * invphi2
               call X%axpby(zero_rsp, Xin, one_rsp)
               call X%axpby(one_rsp, increment, alpha(2))
               call sys%eval(X, residual)
               res(2) = residual%norm()
            else
               ! alphas
               alpha(1:2) = alpha(2:3)
               ! a4 is kept
               ! residuals
               res(1:2) = res(2:3)
               ! r4 is kept
               ! new point --> a3, r3
               alpha(3) = alpha(1) + step * invphi
               call X%axpby(zero_rsp, Xin, one_rsp)
               call X%axpby(one_rsp, increment, alpha(3))
               call sys%eval(X, residual)
               res(3) = residual%norm()
            end if
         end do
         ! set new vector to optimal step
         idx = minloc(res)
         call X%axpby(zero_rsp, Xin, one_rsp)
         call X%axpby(one_rsp, increment, alpha(idx(1)))
      end if

      return
   end subroutine

   subroutine increment_bisection_rdp(X, sys, increment, rold, maxstep, verb)
      !! Classic 1D bisection method based on the golden ratio to damped the Newton step in 
      !! order to maximally reduce the residual at each iteration.
      class(abstract_vector_rdp), intent(inout) :: X
      !! Current system state to be updated
      class(abstract_system_rdp), intent(in)    :: sys
      !! Dynamical system for which the residual is minimized
      class(abstract_vector_rdp), intent(in)    :: increment
      !! Newton step computed from the standard method
      real(dp),                   intent(in)    :: rold
      !! Residual of the current system state to determine improvement
      integer,                    intent(in)    :: maxstep
      !! Maximum number of bisection steps. Each additional bisection step requires an evaluation of the nonlinear function
      logical,                    intent(in)    :: verb
      !! Verbosity toggle
      ! internals
      integer :: i, j, idx(1)
      real(dp) :: invphi, invphi2
      real(dp) :: alpha(4), step
      real(dp) :: res(4)
      class(abstract_vector_rdp), allocatable :: Xin, residual

      allocate(Xin, source=X)
      allocate(residual, source=X); call residual%zero()
      step    = one_rdp
      invphi  = (sqrt(5.0) - 1.0)/2.0  ! 1 / phi
      invphi2 = (3.0 - sqrt(5.0))/2.0  ! 1 / phi**2
      alpha = (/ zero_rdp, invphi2*one_rdp, invphi*one_rdp, one_rdp /)
      res   = (/ rold, zero_rdp, zero_rdp, zero_rdp /)

      call X%add(increment)
      ! evaluate residual norm
      call sys%eval(X, residual)
      res(4) = residual%norm()

      if (res(4) > rold) then
         if (verb) print *, 'Start Newton step bisection ...'
         ! compute new trial solutions
         do j = 2, 3
            call X%axpby(zero_rdp, Xin, one_rdp)
            call X%axpby(one_rdp, increment, alpha(j))
            call sys%eval(X, residual)
            res(j) = residual%norm()
         end do

         do i = 1, maxstep
            step = step * invphi
            if (res(2) < res(3)) then
               ! alphas
               ! a1 is kept
               alpha(3:4) = alpha(2:3)           
               ! residuals
               ! r1 is kept
               res(3:4) = res(2:3)
               ! new point --> a2, r2
               alpha(2) = alpha(1) + step * invphi2
               call X%axpby(zero_rdp, Xin, one_rdp)
               call X%axpby(one_rdp, increment, alpha(2))
               call sys%eval(X, residual)
               res(2) = residual%norm()
            else
               ! alphas
               alpha(1:2) = alpha(2:3)
               ! a4 is kept
               ! residuals
               res(1:2) = res(2:3)
               ! r4 is kept
               ! new point --> a3, r3
               alpha(3) = alpha(1) + step * invphi
               call X%axpby(zero_rdp, Xin, one_rdp)
               call X%axpby(one_rdp, increment, alpha(3))
               call sys%eval(X, residual)
               res(3) = residual%norm()
            end if
         end do
         ! set new vector to optimal step
         idx = minloc(res)
         call X%axpby(zero_rdp, Xin, one_rdp)
         call X%axpby(one_rdp, increment, alpha(idx(1)))
      end if

      return
   end subroutine

   subroutine increment_bisection_csp(X, sys, increment, rold, maxstep, verb)
      !! Classic 1D bisection method based on the golden ratio to damped the Newton step in 
      !! order to maximally reduce the residual at each iteration.
      class(abstract_vector_csp), intent(inout) :: X
      !! Current system state to be updated
      class(abstract_system_csp), intent(in)    :: sys
      !! Dynamical system for which the residual is minimized
      class(abstract_vector_csp), intent(in)    :: increment
      !! Newton step computed from the standard method
      real(sp),                   intent(in)    :: rold
      !! Residual of the current system state to determine improvement
      integer,                    intent(in)    :: maxstep
      !! Maximum number of bisection steps. Each additional bisection step requires an evaluation of the nonlinear function
      logical,                    intent(in)    :: verb
      !! Verbosity toggle
      ! internals
      integer :: i, j, idx(1)
      real(sp) :: invphi, invphi2
      complex(sp) :: alpha(4), step
      real(sp) :: res(4)
      class(abstract_vector_csp), allocatable :: Xin, residual

      allocate(Xin, source=X)
      allocate(residual, source=X); call residual%zero()
      step    = one_csp
      invphi  = (sqrt(5.0) - 1.0)/2.0  ! 1 / phi
      invphi2 = (3.0 - sqrt(5.0))/2.0  ! 1 / phi**2
      alpha = (/ zero_csp, invphi2*one_csp, invphi*one_csp, one_csp /)
      res   = (/ rold, zero_rsp, zero_rsp, zero_rsp /)

      call X%add(increment)
      ! evaluate residual norm
      call sys%eval(X, residual)
      res(4) = residual%norm()

      if (res(4) > rold) then
         if (verb) print *, 'Start Newton step bisection ...'
         ! compute new trial solutions
         do j = 2, 3
            call X%axpby(zero_csp, Xin, one_csp)
            call X%axpby(one_csp, increment, alpha(j))
            call sys%eval(X, residual)
            res(j) = residual%norm()
         end do

         do i = 1, maxstep
            step = step * invphi
            if (res(2) < res(3)) then
               ! alphas
               ! a1 is kept
               alpha(3:4) = alpha(2:3)           
               ! residuals
               ! r1 is kept
               res(3:4) = res(2:3)
               ! new point --> a2, r2
               alpha(2) = alpha(1) + step * invphi2
               call X%axpby(zero_csp, Xin, one_csp)
               call X%axpby(one_csp, increment, alpha(2))
               call sys%eval(X, residual)
               res(2) = residual%norm()
            else
               ! alphas
               alpha(1:2) = alpha(2:3)
               ! a4 is kept
               ! residuals
               res(1:2) = res(2:3)
               ! r4 is kept
               ! new point --> a3, r3
               alpha(3) = alpha(1) + step * invphi
               call X%axpby(zero_csp, Xin, one_csp)
               call X%axpby(one_csp, increment, alpha(3))
               call sys%eval(X, residual)
               res(3) = residual%norm()
            end if
         end do
         ! set new vector to optimal step
         idx = minloc(res)
         call X%axpby(zero_csp, Xin, one_csp)
         call X%axpby(one_csp, increment, alpha(idx(1)))
      end if

      return
   end subroutine

   subroutine increment_bisection_cdp(X, sys, increment, rold, maxstep, verb)
      !! Classic 1D bisection method based on the golden ratio to damped the Newton step in 
      !! order to maximally reduce the residual at each iteration.
      class(abstract_vector_cdp), intent(inout) :: X
      !! Current system state to be updated
      class(abstract_system_cdp), intent(in)    :: sys
      !! Dynamical system for which the residual is minimized
      class(abstract_vector_cdp), intent(in)    :: increment
      !! Newton step computed from the standard method
      real(dp),                   intent(in)    :: rold
      !! Residual of the current system state to determine improvement
      integer,                    intent(in)    :: maxstep
      !! Maximum number of bisection steps. Each additional bisection step requires an evaluation of the nonlinear function
      logical,                    intent(in)    :: verb
      !! Verbosity toggle
      ! internals
      integer :: i, j, idx(1)
      real(dp) :: invphi, invphi2
      complex(dp) :: alpha(4), step
      real(dp) :: res(4)
      class(abstract_vector_cdp), allocatable :: Xin, residual

      allocate(Xin, source=X)
      allocate(residual, source=X); call residual%zero()
      step    = one_cdp
      invphi  = (sqrt(5.0) - 1.0)/2.0  ! 1 / phi
      invphi2 = (3.0 - sqrt(5.0))/2.0  ! 1 / phi**2
      alpha = (/ zero_cdp, invphi2*one_cdp, invphi*one_cdp, one_cdp /)
      res   = (/ rold, zero_rdp, zero_rdp, zero_rdp /)

      call X%add(increment)
      ! evaluate residual norm
      call sys%eval(X, residual)
      res(4) = residual%norm()

      if (res(4) > rold) then
         if (verb) print *, 'Start Newton step bisection ...'
         ! compute new trial solutions
         do j = 2, 3
            call X%axpby(zero_cdp, Xin, one_cdp)
            call X%axpby(one_cdp, increment, alpha(j))
            call sys%eval(X, residual)
            res(j) = residual%norm()
         end do

         do i = 1, maxstep
            step = step * invphi
            if (res(2) < res(3)) then
               ! alphas
               ! a1 is kept
               alpha(3:4) = alpha(2:3)           
               ! residuals
               ! r1 is kept
               res(3:4) = res(2:3)
               ! new point --> a2, r2
               alpha(2) = alpha(1) + step * invphi2
               call X%axpby(zero_cdp, Xin, one_cdp)
               call X%axpby(one_cdp, increment, alpha(2))
               call sys%eval(X, residual)
               res(2) = residual%norm()
            else
               ! alphas
               alpha(1:2) = alpha(2:3)
               ! a4 is kept
               ! residuals
               res(1:2) = res(2:3)
               ! r4 is kept
               ! new point --> a3, r3
               alpha(3) = alpha(1) + step * invphi
               call X%axpby(zero_cdp, Xin, one_cdp)
               call X%axpby(one_cdp, increment, alpha(3))
               call sys%eval(X, residual)
               res(3) = residual%norm()
            end if
         end do
         ! set new vector to optimal step
         idx = minloc(res)
         call X%axpby(zero_cdp, Xin, one_cdp)
         call X%axpby(one_cdp, increment, alpha(idx(1)))
      end if

      return
   end subroutine


   !--------------------------------------------------------------------
   !-----     Definition of two basic tolerance schedulers (sp)    -----
   !--------------------------------------------------------------------

   subroutine constant_atol_sp(tol, rnorm, iter, info)
      !! Abstract interface to define tolerance scheduler for the Newton iteration
      real(sp), intent(out) :: tol
      !! Tolerance to be used
      real(sp), intent(in)  :: rnorm
      !! Norm of the residual of the current iterate
      integer,  intent(in)  :: iter
      !! Newton iteration count
      integer,  intent(out)  :: info
      !! Information flag
      tol = 10*atol_sp
      return
   end subroutine constant_atol_sp

   subroutine dynamic_tol_sp(tol, rnorm, iter, info)
      !! Abstract interface to define tolerance scheduler for the Newton iteration
      real(sp), intent(out) :: tol
      !! Tolerance to be used
      real(sp), intent(in)  :: rnorm
      !! Norm of the residual of the current iterate
      integer,  intent(in)  :: iter
      !! Newton iteration count
      integer,  intent(out)  :: info
      !! Information flag
      tol = max(0.1*rnorm, atol_sp)
      return
   end subroutine dynamic_tol_sp

   !--------------------------------------------------------------------
   !-----     Definition of two basic tolerance schedulers (dp)    -----
   !--------------------------------------------------------------------

   subroutine constant_atol_dp(tol, rnorm, iter, info)
      !! Abstract interface to define tolerance scheduler for the Newton iteration
      real(dp), intent(out) :: tol
      !! Tolerance to be used
      real(dp), intent(in)  :: rnorm
      !! Norm of the residual of the current iterate
      integer,  intent(in)  :: iter
      !! Newton iteration count
      integer,  intent(out)  :: info
      !! Information flag
      tol = 10*atol_dp
      return
   end subroutine constant_atol_dp

   subroutine dynamic_tol_dp(tol, rnorm, iter, info)
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
   end subroutine dynamic_tol_dp


end module LightKrylov_NewtonKrylov
