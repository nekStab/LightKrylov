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

   character(len=128), parameter :: this_module = 'LightKrylov_NewtonKrylov'

   public :: newton
   public :: constant_atol_sp
   public :: dynamic_tol_sp
   public :: constant_atol_dp
   public :: dynamic_tol_dp

   interface newton
      !! Implements the simple Newton-Krylov method for finding roots (fixed points) of a nonlinear vector-valued function
      !! \( F(\mathbf{X}) \), i.e. solutions \( \mathbf{X}^* \) such that \( F(\mathbf{X}^*) - \mathbf{X}^* = \mathbf{0} \) 
      !! starting from an initial guess via successive solution increments based on local linearization \( \mathbf{J}_\mathbf{X} \) 
      !! (the Jacobian) of the nonlinear function in the vicinity of the current solution.
      !!
      !! **Algorthmic Features**
      !!
      !! - At iteration \(k\), the standard Newton step \( \mathbf{\delta x}_k\) is computed as the solution of the linear system
      !!   
      !! $$ \mathbf{J}_\mathbf{X_k} \mathbf{\delta x}_k = \mathbf{r}_k $$
      !!
      !!   where \( \mathbf{r}_k = -F(\mathbf{X}_k) \) is the residual of the nonlinear function. The new guess for the fixed
      !!   point is then given by:
      !!
      !! $$ \mathbf{X}_{k+1} = \mathbf{X}_k + \alpha \mathbf{\delta x}_k$$
      !!
      !!   where \( \alpha \in \left( 0, 1 \right] \) parametrizes the step length. The standard Newton algorithm sets \( \alpha = 1 \).
      !!
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
      !!
      !! - The iterative solution of the linear systems has a comparatively low storage footprint.
      !! - If the Newton iteration converges, the convergence is formally asymptotically of second order. Using dynamic
      !!   tolerances and line searches slightly reduce this convergence rate in exchange for a larger convergence region.
      !! 
      !! **Limitations**
      !!
      !! - The method is not guaranteed to converge if the initial guess is too far from the fixed point. 
      !!   If the Newton iteration diverges even with step bisection, the best suggestion is to find a 
      !!   better initial guess. If this is not feasible, some alternatives to improve the convergence 
      !!   of the Newton iteration are possible (but not implemented to date), including various line search
      !!   algorithms and trust region methods (doglog, double dogleg, hookstep, ...).
      !!
      !! **References**
      !!
      !! - Sánchez, J., Net, M., Garcıa-Archilla, B., & Simó, C. (2004). "Newton–Krylov continuation of periodic orbits 
      !!   for Navier–Stokes flows". Journal of Computational Physics, 201(1), 13-33.
      !! - Viswanath, D. (2007). "Recurrent motions within plane Couette turbulence". Journal of Fluid Mechanics, 580, 339-358.
      !! - Duguet, Y., Pringle, C. C. T., Kerswell, R. R. (2008). "Relative periodic orbits in transitional pipe flow". Physics
      !!   of Fluids, 20(11), 114102.
      !! - Frantz, R. A., Loiseau, J. C., & Robinet, J. C. (2023). "Krylov methods for large-scale dynamical systems: Application 
      !!   in fluid dynamics". Applied Mechanics Reviews, 75(3), 030802.
      module procedure newton_rsp
      module procedure newton_rdp
      module procedure newton_csp
      module procedure newton_cdp
   end interface

   abstract interface
      subroutine abstract_scheduler_sp(tol, target_tol, rnorm, iter, info)
         import sp
         !! Abstract interface to define a tolerance scheduler for the Newton iteration
         real(sp), intent(out) :: tol
         !! Tolerance to be used
         real(sp), intent(in) :: target_tol
         !! Target tolerance
         real(sp), intent(in)  :: rnorm
         !! Norm of the residual of the current iterate
         integer,  intent(in)  :: iter
         !! Newton iteration count
         integer,  intent(out)  :: info
         !! Information flag
      end subroutine abstract_scheduler_sp

      subroutine abstract_scheduler_dp(tol, target_tol, rnorm, iter, info)
         import dp
         !! Abstract interface to define a tolerance scheduler for the Newton iteration
         real(dp), intent(out) :: tol
         !! Tolerance to be used
         real(dp), intent(in) :: target_tol
         !! Target tolerance
         real(dp), intent(in)  :: rnorm
         !! Norm of the residual of the current iterate
         integer,  intent(in)  :: iter
         !! Newton iteration count
         integer,  intent(out)  :: info
         !! Information flag
      end subroutine abstract_scheduler_dp

   end interface

contains

   subroutine newton_rsp(sys, X, solver, info, tolerance, options, linear_solver_options, preconditioner, scheduler)
      class(abstract_system_rsp),                         intent(inout) :: sys
      !! Dynamical system for which we wish to compute a fixed point
      class(abstract_vector_rsp),                         intent(inout) :: X
      !! Initial guess for the fixed point, will be overwritten with solution
      procedure(abstract_linear_solver_rsp)                :: solver
      !! Linear solver to be used to find Newton step
      integer,                                            intent(out)   :: info
      !! Information flag
      real(sp),                                 optional, intent(in)    :: tolerance
      real(sp)                                                          :: target_tol
      !! Target absolute solver tolerance
      type(newton_sp_opts),                     optional, intent(in)    :: options
      type(newton_sp_opts)                                              :: opts
      !! Options for the Newton-Krylov iteration
      class(abstract_opts),                     optional, intent(in)    :: linear_solver_options
      !! Options for the linear solver
      class(abstract_precond_rsp),              optional, intent(in)    :: preconditioner
      !! Preconditioner for the linear solver
      procedure(abstract_scheduler_sp),         optional                :: scheduler

      !--------------------------------------
      !-----     Internal variables     -----
      !--------------------------------------
      
      procedure(abstract_scheduler_sp),      pointer :: tolerance_scheduler => null()
      class(abstract_vector_rsp), allocatable        :: residual, increment
      real(sp)           :: rnorm, tol
      logical            :: converged
      integer            :: i, maxiter, maxstep_bisection
      character(len=256) :: msg
      
      ! Newton-solver tolerance
      target_tol = optval(tolerance, atol_sp)
      ! Newton-Krylov options
      if (present(options)) then
         opts = options
      else
         opts = newton_sp_opts()
      end if
      ! Scheduler
      if (present(scheduler)) then
         tolerance_scheduler => scheduler
      else
         tolerance_scheduler => constant_atol_sp
      endif

      ! Initialisation      
      maxiter           = opts%maxiter
      maxstep_bisection = opts%maxstep_bisection
      converged         = .false.
      allocate(residual, source=X); call residual%zero()
      allocate(increment,source=X); call increment%zero()

      call sys%eval(X, residual, target_tol)
      rnorm = residual%norm()
      ! Check for lucky convergence.
      if (rnorm < target_tol) then
         write(msg,'(A)') 'Initial guess is a fixed point to tolerance!'
         call logger%log_warning(msg, module=this_module, procedure='newton_rsp')
         converged = .true.
         return
      end if

      write(msg,'(A)') 'Starting Newton iteration ...'
      call logger%log_information(msg, module=this_module, procedure='newton_rsp')
      ! Newton iteration
      newton: do i = 1, maxiter

         ! Set dynamic tolerances for Newton iteration and linear solves.
         call tolerance_scheduler(tol, target_tol, rnorm, i, info)
         write(msg,"(A,I0,3(A,E9.2))") 'Start step ', i, ': rnorm= ', rnorm, ', tol= ', tol, ', target= ', target_tol
         call logger%log_message(msg, module=this_module, procedure='newton_rsp')

         ! Define the Jacobian
         sys%jacobian%X = X
        
         ! Solve the linear system using GMRES.
         call residual%chsgn()
         call solver(sys%jacobian, residual, increment, info, atol=tol, &
            & preconditioner=preconditioner, options=linear_solver_options, transpose=.false.)
         call check_info(info, 'linear_solver', module=this_module, procedure='newton_rsp')

         ! Update the solution and overwrite X0
         if (opts%ifbisect) then
            call increment_bisection_rsp(X, sys, increment, rnorm, tol, maxstep_bisection)
         else
            call X%add(increment)
         endif

         ! Evaluate new residual
         call sys%eval(X, residual, tol)
         rnorm = residual%norm()

         ! Check for convergence.
         if (rnorm < tol) then
            if (tol >= target_tol .and. tol < 100.0_sp*target_tol) then
               ! the tolerances are not at the target, check the accurate residual                  
               call sys%eval(X, residual, target_tol)
               if (rnorm < target_tol) then
                  write(msg,'(A,I0,A)') 'Newton iteration converged after ',i,' iterations.'
                  call logger%log_message(msg, module=this_module, procedure='newton_rsp')
                  converged = .true.
                  exit newton
               else
                  write(msg,'(A)') 'Dynamic tolerance but not target tolerance reached. Continue iteration.'
                  call logger%log_warning(msg, module=this_module, procedure='newton_rsp')
               end if
            end if
         end if

      enddo newton

      if (.not.converged) then
         write(msg,'(A,I0,A)') 'Newton iteration did not converge within', maxiter, 'steps.'
         call logger%log_warning(msg, module=this_module, procedure='newton_rsp')
         info = -1
      endif

      return
   end subroutine newton_rsp

   subroutine newton_rdp(sys, X, solver, info, tolerance, options, linear_solver_options, preconditioner, scheduler)
      class(abstract_system_rdp),                         intent(inout) :: sys
      !! Dynamical system for which we wish to compute a fixed point
      class(abstract_vector_rdp),                         intent(inout) :: X
      !! Initial guess for the fixed point, will be overwritten with solution
      procedure(abstract_linear_solver_rdp)                :: solver
      !! Linear solver to be used to find Newton step
      integer,                                            intent(out)   :: info
      !! Information flag
      real(dp),                                 optional, intent(in)    :: tolerance
      real(dp)                                                          :: target_tol
      !! Target absolute solver tolerance
      type(newton_dp_opts),                     optional, intent(in)    :: options
      type(newton_dp_opts)                                              :: opts
      !! Options for the Newton-Krylov iteration
      class(abstract_opts),                     optional, intent(in)    :: linear_solver_options
      !! Options for the linear solver
      class(abstract_precond_rdp),              optional, intent(in)    :: preconditioner
      !! Preconditioner for the linear solver
      procedure(abstract_scheduler_dp),         optional                :: scheduler

      !--------------------------------------
      !-----     Internal variables     -----
      !--------------------------------------
      
      procedure(abstract_scheduler_dp),      pointer :: tolerance_scheduler => null()
      class(abstract_vector_rdp), allocatable        :: residual, increment
      real(dp)           :: rnorm, tol
      logical            :: converged
      integer            :: i, maxiter, maxstep_bisection
      character(len=256) :: msg
      
      ! Newton-solver tolerance
      target_tol = optval(tolerance, atol_dp)
      ! Newton-Krylov options
      if (present(options)) then
         opts = options
      else
         opts = newton_dp_opts()
      end if
      ! Scheduler
      if (present(scheduler)) then
         tolerance_scheduler => scheduler
      else
         tolerance_scheduler => constant_atol_dp
      endif

      ! Initialisation      
      maxiter           = opts%maxiter
      maxstep_bisection = opts%maxstep_bisection
      converged         = .false.
      allocate(residual, source=X); call residual%zero()
      allocate(increment,source=X); call increment%zero()

      call sys%eval(X, residual, target_tol)
      rnorm = residual%norm()
      ! Check for lucky convergence.
      if (rnorm < target_tol) then
         write(msg,'(A)') 'Initial guess is a fixed point to tolerance!'
         call logger%log_warning(msg, module=this_module, procedure='newton_rdp')
         converged = .true.
         return
      end if

      write(msg,'(A)') 'Starting Newton iteration ...'
      call logger%log_information(msg, module=this_module, procedure='newton_rdp')
      ! Newton iteration
      newton: do i = 1, maxiter

         ! Set dynamic tolerances for Newton iteration and linear solves.
         call tolerance_scheduler(tol, target_tol, rnorm, i, info)
         write(msg,"(A,I0,3(A,E9.2))") 'Start step ', i, ': rnorm= ', rnorm, ', tol= ', tol, ', target= ', target_tol
         call logger%log_message(msg, module=this_module, procedure='newton_rdp')

         ! Define the Jacobian
         sys%jacobian%X = X
        
         ! Solve the linear system using GMRES.
         call residual%chsgn()
         call solver(sys%jacobian, residual, increment, info, atol=tol, &
            & preconditioner=preconditioner, options=linear_solver_options, transpose=.false.)
         call check_info(info, 'linear_solver', module=this_module, procedure='newton_rdp')

         ! Update the solution and overwrite X0
         if (opts%ifbisect) then
            call increment_bisection_rdp(X, sys, increment, rnorm, tol, maxstep_bisection)
         else
            call X%add(increment)
         endif

         ! Evaluate new residual
         call sys%eval(X, residual, tol)
         rnorm = residual%norm()

         ! Check for convergence.
         if (rnorm < tol) then
            if (tol >= target_tol .and. tol < 100.0_dp*target_tol) then
               ! the tolerances are not at the target, check the accurate residual                  
               call sys%eval(X, residual, target_tol)
               if (rnorm < target_tol) then
                  write(msg,'(A,I0,A)') 'Newton iteration converged after ',i,' iterations.'
                  call logger%log_message(msg, module=this_module, procedure='newton_rdp')
                  converged = .true.
                  exit newton
               else
                  write(msg,'(A)') 'Dynamic tolerance but not target tolerance reached. Continue iteration.'
                  call logger%log_warning(msg, module=this_module, procedure='newton_rdp')
               end if
            end if
         end if

      enddo newton

      if (.not.converged) then
         write(msg,'(A,I0,A)') 'Newton iteration did not converge within', maxiter, 'steps.'
         call logger%log_warning(msg, module=this_module, procedure='newton_rdp')
         info = -1
      endif

      return
   end subroutine newton_rdp

   subroutine newton_csp(sys, X, solver, info, tolerance, options, linear_solver_options, preconditioner, scheduler)
      class(abstract_system_csp),                         intent(inout) :: sys
      !! Dynamical system for which we wish to compute a fixed point
      class(abstract_vector_csp),                         intent(inout) :: X
      !! Initial guess for the fixed point, will be overwritten with solution
      procedure(abstract_linear_solver_csp)                :: solver
      !! Linear solver to be used to find Newton step
      integer,                                            intent(out)   :: info
      !! Information flag
      real(sp),                                 optional, intent(in)    :: tolerance
      real(sp)                                                          :: target_tol
      !! Target absolute solver tolerance
      type(newton_sp_opts),                     optional, intent(in)    :: options
      type(newton_sp_opts)                                              :: opts
      !! Options for the Newton-Krylov iteration
      class(abstract_opts),                     optional, intent(in)    :: linear_solver_options
      !! Options for the linear solver
      class(abstract_precond_csp),              optional, intent(in)    :: preconditioner
      !! Preconditioner for the linear solver
      procedure(abstract_scheduler_sp),         optional                :: scheduler

      !--------------------------------------
      !-----     Internal variables     -----
      !--------------------------------------
      
      procedure(abstract_scheduler_sp),      pointer :: tolerance_scheduler => null()
      class(abstract_vector_csp), allocatable        :: residual, increment
      real(sp)           :: rnorm, tol
      logical            :: converged
      integer            :: i, maxiter, maxstep_bisection
      character(len=256) :: msg
      
      ! Newton-solver tolerance
      target_tol = optval(tolerance, atol_sp)
      ! Newton-Krylov options
      if (present(options)) then
         opts = options
      else
         opts = newton_sp_opts()
      end if
      ! Scheduler
      if (present(scheduler)) then
         tolerance_scheduler => scheduler
      else
         tolerance_scheduler => constant_atol_sp
      endif

      ! Initialisation      
      maxiter           = opts%maxiter
      maxstep_bisection = opts%maxstep_bisection
      converged         = .false.
      allocate(residual, source=X); call residual%zero()
      allocate(increment,source=X); call increment%zero()

      call sys%eval(X, residual, target_tol)
      rnorm = residual%norm()
      ! Check for lucky convergence.
      if (rnorm < target_tol) then
         write(msg,'(A)') 'Initial guess is a fixed point to tolerance!'
         call logger%log_warning(msg, module=this_module, procedure='newton_csp')
         converged = .true.
         return
      end if

      write(msg,'(A)') 'Starting Newton iteration ...'
      call logger%log_information(msg, module=this_module, procedure='newton_csp')
      ! Newton iteration
      newton: do i = 1, maxiter

         ! Set dynamic tolerances for Newton iteration and linear solves.
         call tolerance_scheduler(tol, target_tol, rnorm, i, info)
         write(msg,"(A,I0,3(A,E9.2))") 'Start step ', i, ': rnorm= ', rnorm, ', tol= ', tol, ', target= ', target_tol
         call logger%log_message(msg, module=this_module, procedure='newton_csp')

         ! Define the Jacobian
         sys%jacobian%X = X
        
         ! Solve the linear system using GMRES.
         call residual%chsgn()
         call solver(sys%jacobian, residual, increment, info, atol=tol, &
            & preconditioner=preconditioner, options=linear_solver_options, transpose=.false.)
         call check_info(info, 'linear_solver', module=this_module, procedure='newton_csp')

         ! Update the solution and overwrite X0
         if (opts%ifbisect) then
            call increment_bisection_csp(X, sys, increment, rnorm, tol, maxstep_bisection)
         else
            call X%add(increment)
         endif

         ! Evaluate new residual
         call sys%eval(X, residual, tol)
         rnorm = residual%norm()

         ! Check for convergence.
         if (rnorm < tol) then
            if (tol >= target_tol .and. tol < 100.0_sp*target_tol) then
               ! the tolerances are not at the target, check the accurate residual                  
               call sys%eval(X, residual, target_tol)
               if (rnorm < target_tol) then
                  write(msg,'(A,I0,A)') 'Newton iteration converged after ',i,' iterations.'
                  call logger%log_message(msg, module=this_module, procedure='newton_csp')
                  converged = .true.
                  exit newton
               else
                  write(msg,'(A)') 'Dynamic tolerance but not target tolerance reached. Continue iteration.'
                  call logger%log_warning(msg, module=this_module, procedure='newton_csp')
               end if
            end if
         end if

      enddo newton

      if (.not.converged) then
         write(msg,'(A,I0,A)') 'Newton iteration did not converge within', maxiter, 'steps.'
         call logger%log_warning(msg, module=this_module, procedure='newton_csp')
         info = -1
      endif

      return
   end subroutine newton_csp

   subroutine newton_cdp(sys, X, solver, info, tolerance, options, linear_solver_options, preconditioner, scheduler)
      class(abstract_system_cdp),                         intent(inout) :: sys
      !! Dynamical system for which we wish to compute a fixed point
      class(abstract_vector_cdp),                         intent(inout) :: X
      !! Initial guess for the fixed point, will be overwritten with solution
      procedure(abstract_linear_solver_cdp)                :: solver
      !! Linear solver to be used to find Newton step
      integer,                                            intent(out)   :: info
      !! Information flag
      real(dp),                                 optional, intent(in)    :: tolerance
      real(dp)                                                          :: target_tol
      !! Target absolute solver tolerance
      type(newton_dp_opts),                     optional, intent(in)    :: options
      type(newton_dp_opts)                                              :: opts
      !! Options for the Newton-Krylov iteration
      class(abstract_opts),                     optional, intent(in)    :: linear_solver_options
      !! Options for the linear solver
      class(abstract_precond_cdp),              optional, intent(in)    :: preconditioner
      !! Preconditioner for the linear solver
      procedure(abstract_scheduler_dp),         optional                :: scheduler

      !--------------------------------------
      !-----     Internal variables     -----
      !--------------------------------------
      
      procedure(abstract_scheduler_dp),      pointer :: tolerance_scheduler => null()
      class(abstract_vector_cdp), allocatable        :: residual, increment
      real(dp)           :: rnorm, tol
      logical            :: converged
      integer            :: i, maxiter, maxstep_bisection
      character(len=256) :: msg
      
      ! Newton-solver tolerance
      target_tol = optval(tolerance, atol_dp)
      ! Newton-Krylov options
      if (present(options)) then
         opts = options
      else
         opts = newton_dp_opts()
      end if
      ! Scheduler
      if (present(scheduler)) then
         tolerance_scheduler => scheduler
      else
         tolerance_scheduler => constant_atol_dp
      endif

      ! Initialisation      
      maxiter           = opts%maxiter
      maxstep_bisection = opts%maxstep_bisection
      converged         = .false.
      allocate(residual, source=X); call residual%zero()
      allocate(increment,source=X); call increment%zero()

      call sys%eval(X, residual, target_tol)
      rnorm = residual%norm()
      ! Check for lucky convergence.
      if (rnorm < target_tol) then
         write(msg,'(A)') 'Initial guess is a fixed point to tolerance!'
         call logger%log_warning(msg, module=this_module, procedure='newton_cdp')
         converged = .true.
         return
      end if

      write(msg,'(A)') 'Starting Newton iteration ...'
      call logger%log_information(msg, module=this_module, procedure='newton_cdp')
      ! Newton iteration
      newton: do i = 1, maxiter

         ! Set dynamic tolerances for Newton iteration and linear solves.
         call tolerance_scheduler(tol, target_tol, rnorm, i, info)
         write(msg,"(A,I0,3(A,E9.2))") 'Start step ', i, ': rnorm= ', rnorm, ', tol= ', tol, ', target= ', target_tol
         call logger%log_message(msg, module=this_module, procedure='newton_cdp')

         ! Define the Jacobian
         sys%jacobian%X = X
        
         ! Solve the linear system using GMRES.
         call residual%chsgn()
         call solver(sys%jacobian, residual, increment, info, atol=tol, &
            & preconditioner=preconditioner, options=linear_solver_options, transpose=.false.)
         call check_info(info, 'linear_solver', module=this_module, procedure='newton_cdp')

         ! Update the solution and overwrite X0
         if (opts%ifbisect) then
            call increment_bisection_cdp(X, sys, increment, rnorm, tol, maxstep_bisection)
         else
            call X%add(increment)
         endif

         ! Evaluate new residual
         call sys%eval(X, residual, tol)
         rnorm = residual%norm()

         ! Check for convergence.
         if (rnorm < tol) then
            if (tol >= target_tol .and. tol < 100.0_dp*target_tol) then
               ! the tolerances are not at the target, check the accurate residual                  
               call sys%eval(X, residual, target_tol)
               if (rnorm < target_tol) then
                  write(msg,'(A,I0,A)') 'Newton iteration converged after ',i,' iterations.'
                  call logger%log_message(msg, module=this_module, procedure='newton_cdp')
                  converged = .true.
                  exit newton
               else
                  write(msg,'(A)') 'Dynamic tolerance but not target tolerance reached. Continue iteration.'
                  call logger%log_warning(msg, module=this_module, procedure='newton_cdp')
               end if
            end if
         end if

      enddo newton

      if (.not.converged) then
         write(msg,'(A,I0,A)') 'Newton iteration did not converge within', maxiter, 'steps.'
         call logger%log_warning(msg, module=this_module, procedure='newton_cdp')
         info = -1
      endif

      return
   end subroutine newton_cdp


   subroutine increment_bisection_rsp(X, sys, increment, rold, tol, maxstep)
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
      real(sp),                   intent(in)    :: tol
      integer,                    intent(in)    :: maxstep
      !! Maximum number of bisection steps. Each additional bisection step requires an evaluation of the nonlinear function

      ! internals
      integer :: i, j, idx(1)
      real(sp) :: invphi, invphi2
      real(sp) :: alpha(4), step
      real(sp) :: res(4)
      class(abstract_vector_rsp), allocatable :: Xin, residual
      character(len=256) :: msg

      allocate(Xin, source=X)
      allocate(residual, source=X); call residual%zero()
      step    = one_rsp
      invphi  = (sqrt(5.0) - 1.0)/2.0  ! 1 / phi
      invphi2 = (3.0 - sqrt(5.0))/2.0  ! 1 / phi**2
      alpha = (/ zero_rsp, invphi2*one_rsp, invphi*one_rsp, one_rsp /)
      res   = (/ rold, zero_rsp, zero_rsp, zero_rsp /)

      call X%add(increment)
      ! evaluate residual norm
      call sys%eval(X, residual, tol)
      res(4) = residual%norm()

      if (res(4) > rold) then
         write(msg,'(A)') 'Start Newton step bisection ... '
         call logger%log_information(msg, module=this_module, procedure='increment_bisection_rsp')
         ! compute new trial solutions
         do j = 2, 3
            call copy(X, Xin)
            call X%axpby(one_rsp, increment, alpha(j))
            call sys%eval(X, residual, tol)
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
               call copy(X, Xin)
               call X%axpby(one_rsp, increment, alpha(2))
               call sys%eval(X, residual, tol)
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
               call copy(X, Xin)
               call X%axpby(one_rsp, increment, alpha(3))
               call sys%eval(X, residual, tol)
               res(3) = residual%norm()
            end if
            write(msg,'(4X,I0,2(A,F6.4))') i, ': New interval: ', alpha(1), ' <= alpha <= ', alpha(4)
            call logger%log_information(msg, module=this_module, procedure='increment_bisection_rsp')
         end do
         ! set new vector to optimal step
         idx = minloc(res)
         write(msg,'(A,F6.4)') 'Optimal damping: alpha= ', alpha(idx(1))
         call logger%log_information(msg, module=this_module, procedure='increment_bisection_rsp')
         call copy(X, Xin)
         call X%axpby(one_rsp, increment, alpha(idx(1)))
      else
         write(msg,'(A)') 'Full Newton step reduces the residual. Skip bisection.'
         call logger%log_information(msg, module=this_module, procedure='increment_bisection_rsp')
      end if

      return
   end subroutine

   subroutine increment_bisection_rdp(X, sys, increment, rold, tol, maxstep)
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
      real(dp),                   intent(in)    :: tol
      integer,                    intent(in)    :: maxstep
      !! Maximum number of bisection steps. Each additional bisection step requires an evaluation of the nonlinear function

      ! internals
      integer :: i, j, idx(1)
      real(dp) :: invphi, invphi2
      real(dp) :: alpha(4), step
      real(dp) :: res(4)
      class(abstract_vector_rdp), allocatable :: Xin, residual
      character(len=256) :: msg

      allocate(Xin, source=X)
      allocate(residual, source=X); call residual%zero()
      step    = one_rdp
      invphi  = (sqrt(5.0) - 1.0)/2.0  ! 1 / phi
      invphi2 = (3.0 - sqrt(5.0))/2.0  ! 1 / phi**2
      alpha = (/ zero_rdp, invphi2*one_rdp, invphi*one_rdp, one_rdp /)
      res   = (/ rold, zero_rdp, zero_rdp, zero_rdp /)

      call X%add(increment)
      ! evaluate residual norm
      call sys%eval(X, residual, tol)
      res(4) = residual%norm()

      if (res(4) > rold) then
         write(msg,'(A)') 'Start Newton step bisection ... '
         call logger%log_information(msg, module=this_module, procedure='increment_bisection_rdp')
         ! compute new trial solutions
         do j = 2, 3
            call copy(X, Xin)
            call X%axpby(one_rdp, increment, alpha(j))
            call sys%eval(X, residual, tol)
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
               call copy(X, Xin)
               call X%axpby(one_rdp, increment, alpha(2))
               call sys%eval(X, residual, tol)
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
               call copy(X, Xin)
               call X%axpby(one_rdp, increment, alpha(3))
               call sys%eval(X, residual, tol)
               res(3) = residual%norm()
            end if
            write(msg,'(4X,I0,2(A,F6.4))') i, ': New interval: ', alpha(1), ' <= alpha <= ', alpha(4)
            call logger%log_information(msg, module=this_module, procedure='increment_bisection_rdp')
         end do
         ! set new vector to optimal step
         idx = minloc(res)
         write(msg,'(A,F6.4)') 'Optimal damping: alpha= ', alpha(idx(1))
         call logger%log_information(msg, module=this_module, procedure='increment_bisection_rdp')
         call copy(X, Xin)
         call X%axpby(one_rdp, increment, alpha(idx(1)))
      else
         write(msg,'(A)') 'Full Newton step reduces the residual. Skip bisection.'
         call logger%log_information(msg, module=this_module, procedure='increment_bisection_rdp')
      end if

      return
   end subroutine

   subroutine increment_bisection_csp(X, sys, increment, rold, tol, maxstep)
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
      real(sp),                   intent(in)    :: tol
      integer,                    intent(in)    :: maxstep
      !! Maximum number of bisection steps. Each additional bisection step requires an evaluation of the nonlinear function

      ! internals
      integer :: i, j, idx(1)
      real(sp) :: invphi, invphi2
      complex(sp) :: alpha(4), step
      real(sp) :: res(4)
      class(abstract_vector_csp), allocatable :: Xin, residual
      character(len=256) :: msg

      allocate(Xin, source=X)
      allocate(residual, source=X); call residual%zero()
      step    = one_csp
      invphi  = (sqrt(5.0) - 1.0)/2.0  ! 1 / phi
      invphi2 = (3.0 - sqrt(5.0))/2.0  ! 1 / phi**2
      alpha = (/ zero_csp, invphi2*one_csp, invphi*one_csp, one_csp /)
      res   = (/ rold, zero_rsp, zero_rsp, zero_rsp /)

      call X%add(increment)
      ! evaluate residual norm
      call sys%eval(X, residual, tol)
      res(4) = residual%norm()

      if (res(4) > rold) then
         write(msg,'(A)') 'Start Newton step bisection ... '
         call logger%log_information(msg, module=this_module, procedure='increment_bisection_csp')
         ! compute new trial solutions
         do j = 2, 3
            call copy(X, Xin)
            call X%axpby(one_csp, increment, alpha(j))
            call sys%eval(X, residual, tol)
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
               call copy(X, Xin)
               call X%axpby(one_csp, increment, alpha(2))
               call sys%eval(X, residual, tol)
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
               call copy(X, Xin)
               call X%axpby(one_csp, increment, alpha(3))
               call sys%eval(X, residual, tol)
               res(3) = residual%norm()
            end if
            write(msg,'(4X,I0,2(A,F6.4))') i, ': New interval: ', alpha(1), ' <= alpha <= ', alpha(4)
            call logger%log_information(msg, module=this_module, procedure='increment_bisection_csp')
         end do
         ! set new vector to optimal step
         idx = minloc(res)
         write(msg,'(A,F6.4)') 'Optimal damping: alpha= ', alpha(idx(1))
         call logger%log_information(msg, module=this_module, procedure='increment_bisection_csp')
         call copy(X, Xin)
         call X%axpby(one_csp, increment, alpha(idx(1)))
      else
         write(msg,'(A)') 'Full Newton step reduces the residual. Skip bisection.'
         call logger%log_information(msg, module=this_module, procedure='increment_bisection_csp')
      end if

      return
   end subroutine

   subroutine increment_bisection_cdp(X, sys, increment, rold, tol, maxstep)
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
      real(dp),                   intent(in)    :: tol
      integer,                    intent(in)    :: maxstep
      !! Maximum number of bisection steps. Each additional bisection step requires an evaluation of the nonlinear function

      ! internals
      integer :: i, j, idx(1)
      real(dp) :: invphi, invphi2
      complex(dp) :: alpha(4), step
      real(dp) :: res(4)
      class(abstract_vector_cdp), allocatable :: Xin, residual
      character(len=256) :: msg

      allocate(Xin, source=X)
      allocate(residual, source=X); call residual%zero()
      step    = one_cdp
      invphi  = (sqrt(5.0) - 1.0)/2.0  ! 1 / phi
      invphi2 = (3.0 - sqrt(5.0))/2.0  ! 1 / phi**2
      alpha = (/ zero_cdp, invphi2*one_cdp, invphi*one_cdp, one_cdp /)
      res   = (/ rold, zero_rdp, zero_rdp, zero_rdp /)

      call X%add(increment)
      ! evaluate residual norm
      call sys%eval(X, residual, tol)
      res(4) = residual%norm()

      if (res(4) > rold) then
         write(msg,'(A)') 'Start Newton step bisection ... '
         call logger%log_information(msg, module=this_module, procedure='increment_bisection_cdp')
         ! compute new trial solutions
         do j = 2, 3
            call copy(X, Xin)
            call X%axpby(one_cdp, increment, alpha(j))
            call sys%eval(X, residual, tol)
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
               call copy(X, Xin)
               call X%axpby(one_cdp, increment, alpha(2))
               call sys%eval(X, residual, tol)
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
               call copy(X, Xin)
               call X%axpby(one_cdp, increment, alpha(3))
               call sys%eval(X, residual, tol)
               res(3) = residual%norm()
            end if
            write(msg,'(4X,I0,2(A,F6.4))') i, ': New interval: ', alpha(1), ' <= alpha <= ', alpha(4)
            call logger%log_information(msg, module=this_module, procedure='increment_bisection_cdp')
         end do
         ! set new vector to optimal step
         idx = minloc(res)
         write(msg,'(A,F6.4)') 'Optimal damping: alpha= ', alpha(idx(1))
         call logger%log_information(msg, module=this_module, procedure='increment_bisection_cdp')
         call copy(X, Xin)
         call X%axpby(one_cdp, increment, alpha(idx(1)))
      else
         write(msg,'(A)') 'Full Newton step reduces the residual. Skip bisection.'
         call logger%log_information(msg, module=this_module, procedure='increment_bisection_cdp')
      end if

      return
   end subroutine


   !--------------------------------------------------------------------
   !-----     Definition of two basic tolerance schedulers (sp)    -----
   !--------------------------------------------------------------------

   subroutine constant_atol_sp(tol, target_tol, rnorm, iter, info)
      !! Abstract interface to define tolerance scheduler for the Newton iteration
      real(sp), intent(out) :: tol
      !! Tolerance to be used
      real(sp), intent(in) :: target_tol
      !! Target tolerance
      real(sp), intent(in)  :: rnorm
      !! Norm of the residual of the current iterate
      integer,  intent(in)  :: iter
      !! Newton iteration count
      integer,  intent(out)  :: info
      !! Information flag
      character(len=256) :: msg
      tol = target_tol
      write(msg,'(A,E9.2)') 'Solver tolerance set to tol= ', tol
      call logger%log_information(msg, module=this_module, procedure='constant_atol_sp')
      return
   end subroutine constant_atol_sp

   subroutine dynamic_tol_sp(tol, target_tol, rnorm, iter, info)
      !! Abstract interface to define tolerance scheduler for the Newton iteration
      real(sp), intent(out) :: tol
      !! Tolerance to be used
      real(sp), intent(in) :: target_tol
      !! Target tolerance
      real(sp), intent(in)  :: rnorm
      !! Norm of the residual of the current iterate
      integer,  intent(in)  :: iter
      !! Newton iteration count
      integer,  intent(out)  :: info
      !! Information flag
      ! internals
      real(sp) :: tol_old
      character(len=256) :: msg
      
      tol_old = tol
      tol = max(0.1*rnorm, target_tol)

      if (tol /= tol_old) then
         if (tol == target_tol) then
            write(msg,'(A,E9.2)') 'Solver tolerance set to input target. tol= ', tol
         else
            write(msg,'(A,E9.2)') 'Solver tolerance set to tol= ', tol
         end if
         call logger%log_information(msg, module=this_module, procedure='dynamic_tol_sp')
      end if
      return
   end subroutine dynamic_tol_sp

   !--------------------------------------------------------------------
   !-----     Definition of two basic tolerance schedulers (dp)    -----
   !--------------------------------------------------------------------

   subroutine constant_atol_dp(tol, target_tol, rnorm, iter, info)
      !! Abstract interface to define tolerance scheduler for the Newton iteration
      real(dp), intent(out) :: tol
      !! Tolerance to be used
      real(dp), intent(in) :: target_tol
      !! Target tolerance
      real(dp), intent(in)  :: rnorm
      !! Norm of the residual of the current iterate
      integer,  intent(in)  :: iter
      !! Newton iteration count
      integer,  intent(out)  :: info
      !! Information flag
      character(len=256) :: msg
      tol = target_tol
      write(msg,'(A,E9.2)') 'Solver tolerance set to tol= ', tol
      call logger%log_information(msg, module=this_module, procedure='constant_atol_dp')
      return
   end subroutine constant_atol_dp

   subroutine dynamic_tol_dp(tol, target_tol, rnorm, iter, info)
      !! Abstract interface to define tolerance scheduler for the Newton iteration
      real(dp), intent(out) :: tol
      !! Tolerance to be used
      real(dp), intent(in) :: target_tol
      !! Target tolerance
      real(dp), intent(in)  :: rnorm
      !! Norm of the residual of the current iterate
      integer,  intent(in)  :: iter
      !! Newton iteration count
      integer,  intent(out)  :: info
      !! Information flag
      ! internals
      real(dp) :: tol_old
      character(len=256) :: msg
      
      tol_old = tol
      tol = max(0.1*rnorm, target_tol)

      if (tol /= tol_old) then
         if (tol == target_tol) then
            write(msg,'(A,E9.2)') 'Solver tolerance set to input target. tol= ', tol
         else
            write(msg,'(A,E9.2)') 'Solver tolerance set to tol= ', tol
         end if
         call logger%log_information(msg, module=this_module, procedure='dynamic_tol_dp')
      end if
      return
   end subroutine dynamic_tol_dp


end module LightKrylov_NewtonKrylov
