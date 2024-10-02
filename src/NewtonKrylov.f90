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
    end type

contains

   subroutine newton_rdp(sys, X, info, options, linear_solver, linear_solver_options, preconditioner, scheduler)
      !! Classic no-frills implementation of the Newton-Krylov root-finding algorithm
      class(abstract_system_rdp), intent(inout)  :: sys
      !! Dynamical system for which we wish to compute a fixed point
      class(abstract_vector_rdp), intent(inout)  :: X
      !! Initial guess for the fixed point, will be overwritten with solution
      integer,                    intent(out)    :: info
      !! Information flag
      type(newton_dp_opts), optional, intent(in)   :: options
      type(newton_dp_opts)                          :: opts
      !! Options for the Newton-Krylov iteration
      procedure(abstract_linear_solver_rdp), optional :: linear_solver
      !! Linear solver to be used to find Newton step
      class(abstract_opts), optional, intent(in) :: linear_solver_options
      class(abstract_opts), allocatable :: solver_opts
      !! Options for the linear solver
      class(abstract_precond_rdp), optional, intent(in) :: preconditioner
      class(abstract_precond_rdp), allocatable :: precond
      !! Preconditioner for the linear solver
      procedure(abstract_scheduler_rdp), optional :: scheduler

      !--------------------------------------
      !-----     Internal variables     -----
      !--------------------------------------
      
      ! residual vector
      class(abstract_vector_rdp), allocatable :: residual, increment
      real(dp) :: rnorm, tol
      logical :: converged, has_precond, has_solver_opts, verb
      integer :: i, maxiter, maxstep_bisection
      procedure(abstract_linear_solver_rdp), pointer :: solver => null()
      procedure(abstract_scheduler_rdp), pointer :: tolerance_scheduler => null()
         

      ! Newton-Krylov options
      if (present(options)) then
         opts = options
      else
         opts = newton_dp_opts()
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
         !call gmres(sys%jacobian, residual, increment, info)
         !call check_info(info, 'gmres', module=this_module, procedure='newton_rdp')

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

   subroutine increment_bisection(X, sys, increment, rold, maxstep, verb)
      class(abstract_vector_rdp), intent(inout) :: X
      class(abstract_system_rdp), intent(in)    :: sys
      class(abstract_vector_rdp), intent(in)    :: increment
      real(dp),                   intent(in)    :: rold
      integer,                    intent(in)    :: maxstep
      logical,                    intent(in)    :: verb
      ! internals
      integer :: i, j, idx(1)
      real(dp) :: invphi, invphi2, step
      real(dp) :: alpha(4), res(4)
      class(abstract_vector_rdp), allocatable :: Xin, residual

      allocate(Xin, source=X)
      allocate(residual, source=X); call residual%zero()
      step    = 1.0_dp
      invphi  = (sqrt(5.0) - 1.0)/2.0  ! 1 / phi
      invphi2 = (3.0 - sqrt(5.0))/2.0  ! 1 / phi**2
      alpha = (/ 0.0_dp, invphi2, invphi, 1.0_dp /)
      res   = (/ rold, 0.0_dp, 0.0_dp, 0.0_dp /)

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
      tol = 10*atol_dp
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
