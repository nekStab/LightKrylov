module LightKrylov_NewtonKrylov
    use stdlib_optval, only: optval
    use LightKrylov_Constants
    use LightKrylov_Logger
    use LightKrylov_AbstractVectors
    use LightKrylov_AbstractLinops
    use LightKrylov_AbstractSystems
    use LightKrylov_IterativeSolvers

    implicit none
    private

    character*128, parameter :: this_module = 'LightKrylov_NewtonKrylov'

    public :: newton

    interface newton
       module procedure newton_classic_rdp
    end interface

contains

    subroutine newton_classic_rdp(sys, X0, tol, verb, info)
       !! Classic no-frills implementation of the Newton-Krylov root-finding algorithm
       class(abstract_system_rdp), intent(inout)  :: sys
       !! Dynamical system for which we wish to compute a fixed point
       class(abstract_vector_rdp), intent(inout)  :: X0
       !! Initial guess for the fixed point, will be overwritten with solution
       real(dp),                   intent(in)     :: tol
       !! Solver tolerance
       logical, optional,          intent(in)     :: verb
       !! verbosity toggle
       integer,                    intent(out)    :: info
       !! Information flag

       !--------------------------------------
       !-----     Internal variables     -----
       !--------------------------------------
       ! residual vector
       class(abstract_vector_rdp), allocatable :: residual, increment
       real(dp) :: rnorm
       logical :: converged, verb_
       integer :: i, maxiter

       ! Optional parameters
       verb_   = optval(verb, .false.)

       ! Initialisation
       maxiter = 100
       converged = .false.
       allocate(residual, source=X0); call residual%zero()
       allocate(increment,source=X0); call increment%zero()

       if (verb) write(*,*) 'Starting Newton ...'
       ! Newton iteration
       newton: do i = 1, maxiter

          call sys%eval(X0, residual)
          rnorm = residual%norm()
          if (verb) write(*,*) "Iteration", i, ": Residual norm = ", rnorm

          ! Check for convergence.
          if (rnorm < tol) then
            if (verb) write(*,*) "Convergence achieved."
            converged = .true.
            exit newton
          end if

          ! Define the Jacobian
          call sys%set_base_state(X0)
        
          ! Solve the linear system using GMRES.
          call residual%chsgn()
          call gmres(sys%jacobian, residual, increment, info)
          call check_info(info, 'gmres', module=this_module, procedure='newton_classic_rdp')

          ! Update the solution and overwrite X0
          call X0%add(increment)

       enddo newton

       if (.not.converged) then
         if (verb) write(*,*) 'Newton iteration did not converge within', maxiter, 'steps.'
          info = -1
       endif

       return
    end subroutine newton_classic_rdp

end module LightKrylov_NewtonKrylov
