module LightKrylov_AbstractSystems
    !!  This module provides the abstract types necessary to define an algebraic system of
    !!  nonlinear equations to be solved using the Newton method.
    use stdlib_optval, only: optval
    use LightKrylov_Logger
    use LightKrylov_Constants
    use LightKrylov_Timer_Utils, only: lightkrylov_timer
    use LightKrylov_AbstractVectors
    use LightKrylov_AbstractLinops
    implicit none(type, external)
    private

    character(len=*), parameter :: this_module      = 'LK_Systems'
    character(len=*), parameter :: this_module_long = 'LK_AbstractSystems'

    public :: verify_system_linearization

    ! Base type for abstract systems.
    type, abstract, public :: abstract_system
        !!  Base type to define an abstract system. All other system types defined
        !!  in `LightKrylov` derive from this fundamental one.
        !!
        !!  @warning
        !!  Users should not extend this abstract class to define their own types.
        !!  @endwarning
        private
        integer :: eval_counter = 0
        type(lightkrylov_timer) :: eval_timer = lightkrylov_timer('system eval timer')
    contains
        procedure, pass(self), public :: get_eval_counter
        !! Return eval counter value
        procedure, pass(self), public :: reset_eval_counter
        !! Reset eval counter
        procedure, pass(self), public :: print_timer_info
        !! Print current timing data
        procedure, pass(self), public :: reset_timer => reset_eval_timer
        !! Reset current timing data
        procedure, pass(self), public :: finalize_timer => finalize_eval_timer
        !! Finalize timer and print complete history
    end type abstract_system

    !----------------------------------------------------------------------------
    !-----     ABSTRACT GENERAL real(sp) SYSTEM DEFINITION WITH kind=sp     -----
    !----------------------------------------------------------------------------

    ! Abstract Jacobian linop for kind=sp
    type, abstract, extends(abstract_linop_rsp), public :: abstract_jacobian_linop_rsp
        !! Abstract type for the local linearization of the system around the state X.
        !! This type is intended for use within the abstract_system type.
        class(abstract_vector_rsp), allocatable :: X
        !! System state around which the equatons are linearized.
    contains
    end type abstract_jacobian_linop_rsp

    ! Abstract system for kind=sp.
    type, abstract, extends(abstract_system), public :: abstract_system_rsp
        !! Base type to extend in order to define a real(sp)-valued system for
        !! Newton fixed-point iteration via the Jacobian.
        class(abstract_jacobian_linop_rsp), allocatable :: jacobian
        !! System Jacobian \( \left. \frac{\partial \mathbf{F}}{\partial \mathbf{X}} \right|_{X^*} \).
    contains
        private
        procedure(abstract_eval_rsp), pass(self), deferred, public :: response
        !! Procedure to evaluate the system residual \( \mathbf{Y} = \mathbf{F}(\mathbf{X}) \).
        ! Wrapper including counter increment
        procedure, pass(self), public :: eval => eval_rsp
        !! Wrapper for response including the counter increment
    end type abstract_system_rsp

    abstract interface
        subroutine abstract_eval_rsp(self, vec_in, vec_out, atol)
            !! Interface for the evaluation of the system response.
            use LightKrylov_AbstractVectors
            import abstract_system_rsp, sp
            implicit none(type, external)
            class(abstract_system_rsp), intent(inout)  :: self
            !! System
            class(abstract_vector_rsp), intent(in)  :: vec_in
            !! State
            class(abstract_vector_rsp), intent(out) :: vec_out
            !! Response
            real(sp),                   intent(in)  :: atol
            !! Solver tolerance
        end subroutine abstract_eval_rsp
    end interface

    !----------------------------------------------------------------------------
    !-----     ABSTRACT GENERAL real(dp) SYSTEM DEFINITION WITH kind=dp     -----
    !----------------------------------------------------------------------------

    ! Abstract Jacobian linop for kind=dp
    type, abstract, extends(abstract_linop_rdp), public :: abstract_jacobian_linop_rdp
        !! Abstract type for the local linearization of the system around the state X.
        !! This type is intended for use within the abstract_system type.
        class(abstract_vector_rdp), allocatable :: X
        !! System state around which the equatons are linearized.
    contains
    end type abstract_jacobian_linop_rdp

    ! Abstract system for kind=dp.
    type, abstract, extends(abstract_system), public :: abstract_system_rdp
        !! Base type to extend in order to define a real(dp)-valued system for
        !! Newton fixed-point iteration via the Jacobian.
        class(abstract_jacobian_linop_rdp), allocatable :: jacobian
        !! System Jacobian \( \left. \frac{\partial \mathbf{F}}{\partial \mathbf{X}} \right|_{X^*} \).
    contains
        private
        procedure(abstract_eval_rdp), pass(self), deferred, public :: response
        !! Procedure to evaluate the system residual \( \mathbf{Y} = \mathbf{F}(\mathbf{X}) \).
        ! Wrapper including counter increment
        procedure, pass(self), public :: eval => eval_rdp
        !! Wrapper for response including the counter increment
    end type abstract_system_rdp

    abstract interface
        subroutine abstract_eval_rdp(self, vec_in, vec_out, atol)
            !! Interface for the evaluation of the system response.
            use LightKrylov_AbstractVectors
            import abstract_system_rdp, dp
            implicit none(type, external)
            class(abstract_system_rdp), intent(inout)  :: self
            !! System
            class(abstract_vector_rdp), intent(in)  :: vec_in
            !! State
            class(abstract_vector_rdp), intent(out) :: vec_out
            !! Response
            real(dp),                   intent(in)  :: atol
            !! Solver tolerance
        end subroutine abstract_eval_rdp
    end interface

    !----------------------------------------------------------------------------
    !-----     ABSTRACT GENERAL complex(sp) SYSTEM DEFINITION WITH kind=sp     -----
    !----------------------------------------------------------------------------

    ! Abstract Jacobian linop for kind=sp
    type, abstract, extends(abstract_linop_csp), public :: abstract_jacobian_linop_csp
        !! Abstract type for the local linearization of the system around the state X.
        !! This type is intended for use within the abstract_system type.
        class(abstract_vector_csp), allocatable :: X
        !! System state around which the equatons are linearized.
    contains
    end type abstract_jacobian_linop_csp

    ! Abstract system for kind=sp.
    type, abstract, extends(abstract_system), public :: abstract_system_csp
        !! Base type to extend in order to define a complex(sp)-valued system for
        !! Newton fixed-point iteration via the Jacobian.
        class(abstract_jacobian_linop_csp), allocatable :: jacobian
        !! System Jacobian \( \left. \frac{\partial \mathbf{F}}{\partial \mathbf{X}} \right|_{X^*} \).
    contains
        private
        procedure(abstract_eval_csp), pass(self), deferred, public :: response
        !! Procedure to evaluate the system residual \( \mathbf{Y} = \mathbf{F}(\mathbf{X}) \).
        ! Wrapper including counter increment
        procedure, pass(self), public :: eval => eval_csp
        !! Wrapper for response including the counter increment
    end type abstract_system_csp

    abstract interface
        subroutine abstract_eval_csp(self, vec_in, vec_out, atol)
            !! Interface for the evaluation of the system response.
            use LightKrylov_AbstractVectors
            import abstract_system_csp, sp
            implicit none(type, external)
            class(abstract_system_csp), intent(inout)  :: self
            !! System
            class(abstract_vector_csp), intent(in)  :: vec_in
            !! State
            class(abstract_vector_csp), intent(out) :: vec_out
            !! Response
            real(sp),                   intent(in)  :: atol
            !! Solver tolerance
        end subroutine abstract_eval_csp
    end interface

    !----------------------------------------------------------------------------
    !-----     ABSTRACT GENERAL complex(dp) SYSTEM DEFINITION WITH kind=dp     -----
    !----------------------------------------------------------------------------

    ! Abstract Jacobian linop for kind=dp
    type, abstract, extends(abstract_linop_cdp), public :: abstract_jacobian_linop_cdp
        !! Abstract type for the local linearization of the system around the state X.
        !! This type is intended for use within the abstract_system type.
        class(abstract_vector_cdp), allocatable :: X
        !! System state around which the equatons are linearized.
    contains
    end type abstract_jacobian_linop_cdp

    ! Abstract system for kind=dp.
    type, abstract, extends(abstract_system), public :: abstract_system_cdp
        !! Base type to extend in order to define a complex(dp)-valued system for
        !! Newton fixed-point iteration via the Jacobian.
        class(abstract_jacobian_linop_cdp), allocatable :: jacobian
        !! System Jacobian \( \left. \frac{\partial \mathbf{F}}{\partial \mathbf{X}} \right|_{X^*} \).
    contains
        private
        procedure(abstract_eval_cdp), pass(self), deferred, public :: response
        !! Procedure to evaluate the system residual \( \mathbf{Y} = \mathbf{F}(\mathbf{X}) \).
        ! Wrapper including counter increment
        procedure, pass(self), public :: eval => eval_cdp
        !! Wrapper for response including the counter increment
    end type abstract_system_cdp

    abstract interface
        subroutine abstract_eval_cdp(self, vec_in, vec_out, atol)
            !! Interface for the evaluation of the system response.
            use LightKrylov_AbstractVectors
            import abstract_system_cdp, dp
            implicit none(type, external)
            class(abstract_system_cdp), intent(inout)  :: self
            !! System
            class(abstract_vector_cdp), intent(in)  :: vec_in
            !! State
            class(abstract_vector_cdp), intent(out) :: vec_out
            !! Response
            real(dp),                   intent(in)  :: atol
            !! Solver tolerance
        end subroutine abstract_eval_cdp
    end interface


    !! Taylor test of the system linearization around the state x.
    interface verify_system_linearization
        module procedure taylor_test_system_linearization_rsp
        module procedure taylor_test_system_linearization_rdp
        module procedure taylor_test_system_linearization_csp
        module procedure taylor_test_system_linearization_cdp
    end interface verify_system_linearization
contains

    !---------------------------------------------------------------
    !-----     Getter/Setter routines for abstract_systems     -----
    !---------------------------------------------------------------

    pure integer function get_eval_counter(self) result(count)
        !! Getter function for the number of system evaluations.
        implicit none(type, external)
        class(abstract_system), intent(in) :: self
        count = self%eval_counter
    end function get_eval_counter

    subroutine reset_eval_counter(self, procedure, counter, reset_timer, soft_reset, clean_timer)
        !! Setter function to reset the system evaluation counter.
        implicit none(type, external)
        class(abstract_system), intent(inout) :: self
        character(len=*), intent(in) :: procedure
        !! name of the caller routine
        integer, optional, intent(in) :: counter
        !! optional flag to reset to an integer other than zero.
        logical, optional, intent(in) :: reset_timer
        !! optional flag to reset also the timer
        logical, optional, intent(in) :: soft_reset
        !! optional flag to choose whether to save previous timing data (default: .true.)
        logical, optional, intent(in) :: clean_timer
        !! optional flag to choose whether to fully reset the timer (default: .false.)
        ! internals
        integer :: counter_, count_old
        logical :: reset_timer_
        character(len=128) :: msg
        counter_ = optval(counter, 0)
        count_old = self%get_eval_counter()
        reset_timer_ = optval(reset_timer, .true.)
        if (count_old /= 0 .or. counter_ /= 0) then
            write(msg,'(A,I0,A,I0,A)') 'Total number of evals: ', count_old, '. Resetting counter to ', counter_, '.'
            call log_message(msg, this_module, 'reset_eval_counter('//trim(procedure)//')')
            self%eval_counter = counter_
        end if
        if (reset_timer_) call self%reset_timer(soft_reset, clean_timer)
        return
    end subroutine reset_eval_counter

    subroutine print_timer_info(self)
        !! Print the current timing data for the system evaluation.
        !! Note: Wrapper of the corresponding routine from lightkrylov_timer
        implicit none(type, external)
        class(abstract_system), intent(inout) :: self
        call self%eval_timer%print_info()
    end subroutine print_timer_info

    subroutine reset_eval_timer(self, soft, clean)
        !! Setter routine to reset the system evaluation timer.
        !! Note: Wrapper of the corresponding routine from lightkrylov_timer
        implicit none(type, external)
        class(abstract_system), intent(inout) :: self
        logical, optional, intent(in) :: soft
        logical, optional, intent(in) :: clean
        call self%eval_timer%reset(soft, clean, verbose=.true.)
    end subroutine reset_eval_timer

    subroutine finalize_eval_timer(self)
        !! Finalize the system evaluation timer and print summary.
        !! Note: Wrapper of the corresponding routine from lightkrylov_timer
        implicit none(type, external)
        class(abstract_system), intent(inout) :: self
        call self%eval_timer%finalize()
    end subroutine finalize_eval_timer

    !----------------------------------------------------------------------
    !-----     Wrappers for system response to increment counters     -----
    !----------------------------------------------------------------------

    subroutine eval_rsp(self, vec_in, vec_out, atol)
        implicit none(type, external)
        class(abstract_system_rsp), intent(inout) :: self
        class(abstract_vector_rsp), intent(in)    :: vec_in
        class(abstract_vector_rsp), intent(out)   :: vec_out
        real(sp),                             intent(in)    :: atol
        ! internal
        character(len=128) :: msg
        self%eval_counter = self%eval_counter + 1
        write(msg,'(I0,1X,A)') self%eval_counter, 'start'
        call log_debug(msg, this_module, 'response')
        call self%eval_timer%start()
        call self%response(vec_in, vec_out, atol)
        call self%eval_timer%stop()
        write(msg,'(I0,1X,A)') self%eval_counter, 'end'
        call log_debug(msg, this_module, 'response')
        return
    end subroutine eval_rsp
    subroutine eval_rdp(self, vec_in, vec_out, atol)
        implicit none(type, external)
        class(abstract_system_rdp), intent(inout) :: self
        class(abstract_vector_rdp), intent(in)    :: vec_in
        class(abstract_vector_rdp), intent(out)   :: vec_out
        real(dp),                             intent(in)    :: atol
        ! internal
        character(len=128) :: msg
        self%eval_counter = self%eval_counter + 1
        write(msg,'(I0,1X,A)') self%eval_counter, 'start'
        call log_debug(msg, this_module, 'response')
        call self%eval_timer%start()
        call self%response(vec_in, vec_out, atol)
        call self%eval_timer%stop()
        write(msg,'(I0,1X,A)') self%eval_counter, 'end'
        call log_debug(msg, this_module, 'response')
        return
    end subroutine eval_rdp
    subroutine eval_csp(self, vec_in, vec_out, atol)
        implicit none(type, external)
        class(abstract_system_csp), intent(inout) :: self
        class(abstract_vector_csp), intent(in)    :: vec_in
        class(abstract_vector_csp), intent(out)   :: vec_out
        real(sp),                             intent(in)    :: atol
        ! internal
        character(len=128) :: msg
        self%eval_counter = self%eval_counter + 1
        write(msg,'(I0,1X,A)') self%eval_counter, 'start'
        call log_debug(msg, this_module, 'response')
        call self%eval_timer%start()
        call self%response(vec_in, vec_out, atol)
        call self%eval_timer%stop()
        write(msg,'(I0,1X,A)') self%eval_counter, 'end'
        call log_debug(msg, this_module, 'response')
        return
    end subroutine eval_csp
    subroutine eval_cdp(self, vec_in, vec_out, atol)
        implicit none(type, external)
        class(abstract_system_cdp), intent(inout) :: self
        class(abstract_vector_cdp), intent(in)    :: vec_in
        class(abstract_vector_cdp), intent(out)   :: vec_out
        real(dp),                             intent(in)    :: atol
        ! internal
        character(len=128) :: msg
        self%eval_counter = self%eval_counter + 1
        write(msg,'(I0,1X,A)') self%eval_counter, 'start'
        call log_debug(msg, this_module, 'response')
        call self%eval_timer%start()
        call self%response(vec_in, vec_out, atol)
        call self%eval_timer%stop()
        write(msg,'(I0,1X,A)') self%eval_counter, 'end'
        call log_debug(msg, this_module, 'response')
        return
    end subroutine eval_cdp

    !--------------------------------------
    !-----      UTILITY FUNCTIONS     -----
    !--------------------------------------

    subroutine taylor_test_system_linearization_rsp(sys, nsteps, eps0, stepratio, tolerance, verbose, orders)
        !! Taylor test of the system linearization around the state sys%jacobian%x.
        !!
        !! For a perturbation direction d and decreasing steps h, prints
        !!   ||F(x + alpha delta) - F(x)||                    ~ O(|alpha|)    for any F,
        !!   ||F(x + alpha delta) - F(x) - alpha J(x) delta|| ~ O(|alpha|^2)  if J is the correct linearization,
        !! together with the observed convergence order of the remainder.
        !! NOTE: The test assumes that the system is allocated, the jacobian and the baseflow are set and initialized
        class(abstract_system_rsp), intent(inout)  :: sys
        !! System whose linearization is to be tested
        integer, optional,          intent(in)  :: nsteps
        !! number of steps to use in the Taylor test (default: 5)
        real(sp), optional,         intent(in)  :: eps0
        !! initial step size for the Taylor test (default: 1.0e-2)
        real(sp), optional,         intent(in)  :: stepratio
        !! ratio between consecutive step sizes (default: 0.1)
        real(sp), optional,         intent(in)  :: tolerance
        !! tolerance for the system evaluation (default: 1.0e-8)
        logical, optional,            intent(in)  :: verbose
        !! flag to print the results of the Taylor test (default: .true.)
        real(sp), optional, allocatable, intent(out) :: orders(:)
        !! observed convergence orders of the remainder between consecutive steps (size nsteps-1)
        ! internal
        class(abstract_vector_rsp), allocatable :: X0, Fx, dF, alpha_x, delta, Jdelta, wrk
        integer :: nsteps_, k, ipass, npass
        real(sp) :: tol, xnorm, eps, eps0_, q, r0, r1, r1_prev, order
        logical :: verb
        real(sp) :: alpha
        character(len=256) :: msg
        character(len=*), parameter :: this_procedure = 'verify_system_linearization_rsp'

        if (.not. allocated(sys%jacobian)) then
            call stop_error('System Jacobian is not allocated.', this_module, this_procedure)
        end if
        if (.not. allocated(sys%jacobian%X)) then
            call stop_error('Jacobian base state X is not set.', this_module, this_procedure)
        end if

        !> Deals with optional arguments.
        nsteps_  = optval(nsteps, 4)
        q        = optval(stepratio, 0.1_sp)
        eps0_    = optval(eps0, 5.0e-2_sp)
        tol      = optval(tolerance, 10.0_sp**(-(precision(1.0_sp)-1)))
        verb     = optval(verbose, .true.)
        npass    = 1
        if (present(orders)) then
            allocate(orders(max(nsteps_-1, 0)*npass))
            orders = 0.0_sp
        end if

        !> Reference response.
        allocate(X0, Fx, dF, alpha_x, delta, Jdelta, wrk, mold=sys%jacobian%X)
        call X0%zero()
        call X0%add(sys%jacobian%X) 
        xnorm = X0%norm()
        call sys%eval(X0, Fx, tol)

        if (verb) call log_message('Taylor test of the system linearization', this_module, this_procedure)
        write(msg, '(3(A,ES12.5))') '   ||x|| = ', xnorm, ', ||F(x)|| = ', Fx%norm(), ', atol = ', tol
        if (verb) call log_message(msg, this_module, this_procedure)

        !> Perturbation direction: delta = J u / ||J u|| (random vector, apply J to satisfy BCs).
        call wrk%rand(ifnorm=.true.)
        call sys%jacobian%apply_matvec(wrk, delta)
        call delta%scal(one_rsp / delta%norm())

        !> Linearized response J delta.
        call sys%jacobian%apply_matvec(delta, Jdelta)

        write(msg, '(A,ES12.5)') '   ||J delta|| = ', Jdelta%norm()
        if (verb) call log_message(msg, this_module, this_procedure)

        do ipass = 1, npass
            write(msg, '(3X,A12,2X,A16,2X,A16,2X,A6)') 'e', '||dF||', '||dF - e*Jd||', 'order'

            r1_prev = 0.0_sp
            steps: do k = 1, nsteps_
                eps = eps0_* q**(k-1)
                alpha = eps
                
                !> Perturbed state x + alpha * delta.
                call alpha_x%zero()
                call alpha_x%add(X0)
                call alpha_x%axpby(alpha, delta, one_rsp)

                !> dF = F(x + alpha * delta) - F(x).
                call sys%eval(alpha_x, dF, tol)
                call dF%sub(Fx)
                r0 = dF%norm()

                !> Remainder dF - alpha * J * delta.
                call dF%axpby(-alpha, Jdelta, one_rsp)
                r1 = dF%norm()

                !> Print results and convergence order.
                if (k == 1 .and. verb) call log_message(msg, this_module, this_procedure)
                if (k > 1 .and. r1 > 0.0_sp .and. r1_prev > 0.0_sp) then
                    order = log(r1_prev / r1) / log(1.0_sp / q)
                    if (present(orders)) orders((ipass-1)*(nsteps_-1) + k-1) = order
                    write(msg, '(3X,ES12.5,2X,ES16.5,2X,ES16.5,2X,F6.3)') eps, r0, r1, order
                else
                    write(msg, '(3X,ES12.5,2X,ES16.5,2X,ES16.5,2X,A6)') eps, r0, r1, '-'
                end if
                if (verb) call log_message(msg, this_module, this_procedure)

                r1_prev = r1
            end do steps
        end do
        deallocate(X0, Fx, dF, alpha_x, delta, Jdelta, wrk)
    end subroutine taylor_test_system_linearization_rsp

    subroutine taylor_test_system_linearization_rdp(sys, nsteps, eps0, stepratio, tolerance, verbose, orders)
        !! Taylor test of the system linearization around the state sys%jacobian%x.
        !!
        !! For a perturbation direction d and decreasing steps h, prints
        !!   ||F(x + alpha delta) - F(x)||                    ~ O(|alpha|)    for any F,
        !!   ||F(x + alpha delta) - F(x) - alpha J(x) delta|| ~ O(|alpha|^2)  if J is the correct linearization,
        !! together with the observed convergence order of the remainder.
        !! NOTE: The test assumes that the system is allocated, the jacobian and the baseflow are set and initialized
        class(abstract_system_rdp), intent(inout)  :: sys
        !! System whose linearization is to be tested
        integer, optional,          intent(in)  :: nsteps
        !! number of steps to use in the Taylor test (default: 5)
        real(dp), optional,         intent(in)  :: eps0
        !! initial step size for the Taylor test (default: 1.0e-2)
        real(dp), optional,         intent(in)  :: stepratio
        !! ratio between consecutive step sizes (default: 0.1)
        real(dp), optional,         intent(in)  :: tolerance
        !! tolerance for the system evaluation (default: 1.0e-8)
        logical, optional,            intent(in)  :: verbose
        !! flag to print the results of the Taylor test (default: .true.)
        real(dp), optional, allocatable, intent(out) :: orders(:)
        !! observed convergence orders of the remainder between consecutive steps (size nsteps-1)
        ! internal
        class(abstract_vector_rdp), allocatable :: X0, Fx, dF, alpha_x, delta, Jdelta, wrk
        integer :: nsteps_, k, ipass, npass
        real(dp) :: tol, xnorm, eps, eps0_, q, r0, r1, r1_prev, order
        logical :: verb
        real(dp) :: alpha
        character(len=256) :: msg
        character(len=*), parameter :: this_procedure = 'verify_system_linearization_rdp'

        if (.not. allocated(sys%jacobian)) then
            call stop_error('System Jacobian is not allocated.', this_module, this_procedure)
        end if
        if (.not. allocated(sys%jacobian%X)) then
            call stop_error('Jacobian base state X is not set.', this_module, this_procedure)
        end if

        !> Deals with optional arguments.
        nsteps_  = optval(nsteps, 4)
        q        = optval(stepratio, 0.1_dp)
        eps0_    = optval(eps0, 5.0e-2_dp)
        tol      = optval(tolerance, 10.0_dp**(-(precision(1.0_dp)-1)))
        verb     = optval(verbose, .true.)
        npass    = 1
        if (present(orders)) then
            allocate(orders(max(nsteps_-1, 0)*npass))
            orders = 0.0_dp
        end if

        !> Reference response.
        allocate(X0, Fx, dF, alpha_x, delta, Jdelta, wrk, mold=sys%jacobian%X)
        call X0%zero()
        call X0%add(sys%jacobian%X) 
        xnorm = X0%norm()
        call sys%eval(X0, Fx, tol)

        if (verb) call log_message('Taylor test of the system linearization', this_module, this_procedure)
        write(msg, '(3(A,ES12.5))') '   ||x|| = ', xnorm, ', ||F(x)|| = ', Fx%norm(), ', atol = ', tol
        if (verb) call log_message(msg, this_module, this_procedure)

        !> Perturbation direction: delta = J u / ||J u|| (random vector, apply J to satisfy BCs).
        call wrk%rand(ifnorm=.true.)
        call sys%jacobian%apply_matvec(wrk, delta)
        call delta%scal(one_rdp / delta%norm())

        !> Linearized response J delta.
        call sys%jacobian%apply_matvec(delta, Jdelta)

        write(msg, '(A,ES12.5)') '   ||J delta|| = ', Jdelta%norm()
        if (verb) call log_message(msg, this_module, this_procedure)

        do ipass = 1, npass
            write(msg, '(3X,A12,2X,A16,2X,A16,2X,A6)') 'e', '||dF||', '||dF - e*Jd||', 'order'

            r1_prev = 0.0_dp
            steps: do k = 1, nsteps_
                eps = eps0_* q**(k-1)
                alpha = eps
                
                !> Perturbed state x + alpha * delta.
                call alpha_x%zero()
                call alpha_x%add(X0)
                call alpha_x%axpby(alpha, delta, one_rdp)

                !> dF = F(x + alpha * delta) - F(x).
                call sys%eval(alpha_x, dF, tol)
                call dF%sub(Fx)
                r0 = dF%norm()

                !> Remainder dF - alpha * J * delta.
                call dF%axpby(-alpha, Jdelta, one_rdp)
                r1 = dF%norm()

                !> Print results and convergence order.
                if (k == 1 .and. verb) call log_message(msg, this_module, this_procedure)
                if (k > 1 .and. r1 > 0.0_dp .and. r1_prev > 0.0_dp) then
                    order = log(r1_prev / r1) / log(1.0_dp / q)
                    if (present(orders)) orders((ipass-1)*(nsteps_-1) + k-1) = order
                    write(msg, '(3X,ES12.5,2X,ES16.5,2X,ES16.5,2X,F6.3)') eps, r0, r1, order
                else
                    write(msg, '(3X,ES12.5,2X,ES16.5,2X,ES16.5,2X,A6)') eps, r0, r1, '-'
                end if
                if (verb) call log_message(msg, this_module, this_procedure)

                r1_prev = r1
            end do steps
        end do
        deallocate(X0, Fx, dF, alpha_x, delta, Jdelta, wrk)
    end subroutine taylor_test_system_linearization_rdp

    subroutine taylor_test_system_linearization_csp(sys, nsteps, eps0, stepratio, tolerance, verbose, orders)
        !! Taylor test of the system linearization around the state sys%jacobian%x.
        !!
        !! For a perturbation direction d and decreasing steps h, prints
        !!   ||F(x + alpha delta) - F(x)||                    ~ O(|alpha|)    for any F,
        !!   ||F(x + alpha delta) - F(x) - alpha J(x) delta|| ~ O(|alpha|^2)  if J is the correct linearization,
        !! together with the observed convergence order of the remainder.
        !! NOTE: The test assumes that the system is allocated, the jacobian and the baseflow are set and initialized
        class(abstract_system_csp), intent(inout)  :: sys
        !! System whose linearization is to be tested
        integer, optional,          intent(in)  :: nsteps
        !! number of steps to use in the Taylor test (default: 5)
        real(sp), optional,         intent(in)  :: eps0
        !! initial step size for the Taylor test (default: 1.0e-2)
        real(sp), optional,         intent(in)  :: stepratio
        !! ratio between consecutive step sizes (default: 0.1)
        real(sp), optional,         intent(in)  :: tolerance
        !! tolerance for the system evaluation (default: 1.0e-8)
        logical, optional,            intent(in)  :: verbose
        !! flag to print the results of the Taylor test (default: .true.)
        real(sp), optional, allocatable, intent(out) :: orders(:)
        !! observed convergence orders of the remainder between consecutive steps (size nsteps-1)
        ! internal
        class(abstract_vector_csp), allocatable :: X0, Fx, dF, alpha_x, delta, Jdelta, wrk
        integer :: nsteps_, k, ipass, npass
        real(sp) :: tol, xnorm, eps, eps0_, q, r0, r1, r1_prev, order
        logical :: verb
        complex(sp) :: alpha
        character(len=256) :: msg
        character(len=*), parameter :: this_procedure = 'verify_system_linearization_csp'

        if (.not. allocated(sys%jacobian)) then
            call stop_error('System Jacobian is not allocated.', this_module, this_procedure)
        end if
        if (.not. allocated(sys%jacobian%X)) then
            call stop_error('Jacobian base state X is not set.', this_module, this_procedure)
        end if

        !> Deals with optional arguments.
        nsteps_  = optval(nsteps, 4)
        q        = optval(stepratio, 0.1_sp)
        eps0_    = optval(eps0, 5.0e-2_sp)
        tol      = optval(tolerance, 10.0_sp**(-(precision(1.0_sp)-1)))
        verb     = optval(verbose, .true.)
        npass    = 2 ! real and imag perturbations run separately for complex systems.
        if (present(orders)) then
            allocate(orders(max(nsteps_-1, 0)*npass))
            orders = 0.0_sp
        end if

        !> Reference response.
        allocate(X0, Fx, dF, alpha_x, delta, Jdelta, wrk, mold=sys%jacobian%X)
        call X0%zero()
        call X0%add(sys%jacobian%X) 
        xnorm = X0%norm()
        call sys%eval(X0, Fx, tol)

        if (verb) call log_message('Taylor test of the system linearization', this_module, this_procedure)
        write(msg, '(3(A,ES12.5))') '   ||x|| = ', xnorm, ', ||F(x)|| = ', Fx%norm(), ', atol = ', tol
        if (verb) call log_message(msg, this_module, this_procedure)

        !> Perturbation direction: delta = J u / ||J u|| (random vector, apply J to satisfy BCs).
        call wrk%rand(ifnorm=.true.)
        call sys%jacobian%apply_matvec(wrk, delta)
        call delta%scal(one_csp / delta%norm())

        !> Linearized response J delta.
        call sys%jacobian%apply_matvec(delta, Jdelta)

        write(msg, '(A,ES12.5)') '   ||J delta|| = ', Jdelta%norm()
        if (verb) call log_message(msg, this_module, this_procedure)

        do ipass = 1, npass
            msg = merge('Pass 1: real perturbation', 'Pass 2: imag perturbation', ipass == 1)
            if (verb) call log_message(msg, this_module, this_procedure)
            write(msg, '(3X,A12,2X,A16,2X,A16,2X,A6)') 'e', '||dF||', '||dF - e*Jd||', 'order'

            r1_prev = 0.0_sp
            steps: do k = 1, nsteps_
                eps = eps0_* q**(k-1)
                alpha = merge(eps*one_csp, eps*one_im_csp, ipass == 1)
                
                !> Perturbed state x + alpha * delta.
                call alpha_x%zero()
                call alpha_x%add(X0)
                call alpha_x%axpby(alpha, delta, one_csp)

                !> dF = F(x + alpha * delta) - F(x).
                call sys%eval(alpha_x, dF, tol)
                call dF%sub(Fx)
                r0 = dF%norm()

                !> Remainder dF - alpha * J * delta.
                call dF%axpby(-alpha, Jdelta, one_csp)
                r1 = dF%norm()

                !> Print results and convergence order.
                if (k == 1 .and. verb) call log_message(msg, this_module, this_procedure)
                if (k > 1 .and. r1 > 0.0_sp .and. r1_prev > 0.0_sp) then
                    order = log(r1_prev / r1) / log(1.0_sp / q)
                    if (present(orders)) orders((ipass-1)*(nsteps_-1) + k-1) = order
                    write(msg, '(3X,ES12.5,2X,ES16.5,2X,ES16.5,2X,F6.3)') eps, r0, r1, order
                else
                    write(msg, '(3X,ES12.5,2X,ES16.5,2X,ES16.5,2X,A6)') eps, r0, r1, '-'
                end if
                if (verb) call log_message(msg, this_module, this_procedure)

                r1_prev = r1
            end do steps
        end do
        deallocate(X0, Fx, dF, alpha_x, delta, Jdelta, wrk)
    end subroutine taylor_test_system_linearization_csp

    subroutine taylor_test_system_linearization_cdp(sys, nsteps, eps0, stepratio, tolerance, verbose, orders)
        !! Taylor test of the system linearization around the state sys%jacobian%x.
        !!
        !! For a perturbation direction d and decreasing steps h, prints
        !!   ||F(x + alpha delta) - F(x)||                    ~ O(|alpha|)    for any F,
        !!   ||F(x + alpha delta) - F(x) - alpha J(x) delta|| ~ O(|alpha|^2)  if J is the correct linearization,
        !! together with the observed convergence order of the remainder.
        !! NOTE: The test assumes that the system is allocated, the jacobian and the baseflow are set and initialized
        class(abstract_system_cdp), intent(inout)  :: sys
        !! System whose linearization is to be tested
        integer, optional,          intent(in)  :: nsteps
        !! number of steps to use in the Taylor test (default: 5)
        real(dp), optional,         intent(in)  :: eps0
        !! initial step size for the Taylor test (default: 1.0e-2)
        real(dp), optional,         intent(in)  :: stepratio
        !! ratio between consecutive step sizes (default: 0.1)
        real(dp), optional,         intent(in)  :: tolerance
        !! tolerance for the system evaluation (default: 1.0e-8)
        logical, optional,            intent(in)  :: verbose
        !! flag to print the results of the Taylor test (default: .true.)
        real(dp), optional, allocatable, intent(out) :: orders(:)
        !! observed convergence orders of the remainder between consecutive steps (size nsteps-1)
        ! internal
        class(abstract_vector_cdp), allocatable :: X0, Fx, dF, alpha_x, delta, Jdelta, wrk
        integer :: nsteps_, k, ipass, npass
        real(dp) :: tol, xnorm, eps, eps0_, q, r0, r1, r1_prev, order
        logical :: verb
        complex(dp) :: alpha
        character(len=256) :: msg
        character(len=*), parameter :: this_procedure = 'verify_system_linearization_cdp'

        if (.not. allocated(sys%jacobian)) then
            call stop_error('System Jacobian is not allocated.', this_module, this_procedure)
        end if
        if (.not. allocated(sys%jacobian%X)) then
            call stop_error('Jacobian base state X is not set.', this_module, this_procedure)
        end if

        !> Deals with optional arguments.
        nsteps_  = optval(nsteps, 4)
        q        = optval(stepratio, 0.1_dp)
        eps0_    = optval(eps0, 5.0e-2_dp)
        tol      = optval(tolerance, 10.0_dp**(-(precision(1.0_dp)-1)))
        verb     = optval(verbose, .true.)
        npass    = 2 ! real and imag perturbations run separately for complex systems.
        if (present(orders)) then
            allocate(orders(max(nsteps_-1, 0)*npass))
            orders = 0.0_dp
        end if

        !> Reference response.
        allocate(X0, Fx, dF, alpha_x, delta, Jdelta, wrk, mold=sys%jacobian%X)
        call X0%zero()
        call X0%add(sys%jacobian%X) 
        xnorm = X0%norm()
        call sys%eval(X0, Fx, tol)

        if (verb) call log_message('Taylor test of the system linearization', this_module, this_procedure)
        write(msg, '(3(A,ES12.5))') '   ||x|| = ', xnorm, ', ||F(x)|| = ', Fx%norm(), ', atol = ', tol
        if (verb) call log_message(msg, this_module, this_procedure)

        !> Perturbation direction: delta = J u / ||J u|| (random vector, apply J to satisfy BCs).
        call wrk%rand(ifnorm=.true.)
        call sys%jacobian%apply_matvec(wrk, delta)
        call delta%scal(one_cdp / delta%norm())

        !> Linearized response J delta.
        call sys%jacobian%apply_matvec(delta, Jdelta)

        write(msg, '(A,ES12.5)') '   ||J delta|| = ', Jdelta%norm()
        if (verb) call log_message(msg, this_module, this_procedure)

        do ipass = 1, npass
            msg = merge('Pass 1: real perturbation', 'Pass 2: imag perturbation', ipass == 1)
            if (verb) call log_message(msg, this_module, this_procedure)
            write(msg, '(3X,A12,2X,A16,2X,A16,2X,A6)') 'e', '||dF||', '||dF - e*Jd||', 'order'

            r1_prev = 0.0_dp
            steps: do k = 1, nsteps_
                eps = eps0_* q**(k-1)
                alpha = merge(eps*one_cdp, eps*one_im_cdp, ipass == 1)
                
                !> Perturbed state x + alpha * delta.
                call alpha_x%zero()
                call alpha_x%add(X0)
                call alpha_x%axpby(alpha, delta, one_cdp)

                !> dF = F(x + alpha * delta) - F(x).
                call sys%eval(alpha_x, dF, tol)
                call dF%sub(Fx)
                r0 = dF%norm()

                !> Remainder dF - alpha * J * delta.
                call dF%axpby(-alpha, Jdelta, one_cdp)
                r1 = dF%norm()

                !> Print results and convergence order.
                if (k == 1 .and. verb) call log_message(msg, this_module, this_procedure)
                if (k > 1 .and. r1 > 0.0_dp .and. r1_prev > 0.0_dp) then
                    order = log(r1_prev / r1) / log(1.0_dp / q)
                    if (present(orders)) orders((ipass-1)*(nsteps_-1) + k-1) = order
                    write(msg, '(3X,ES12.5,2X,ES16.5,2X,ES16.5,2X,F6.3)') eps, r0, r1, order
                else
                    write(msg, '(3X,ES12.5,2X,ES16.5,2X,ES16.5,2X,A6)') eps, r0, r1, '-'
                end if
                if (verb) call log_message(msg, this_module, this_procedure)

                r1_prev = r1
            end do steps
        end do
        deallocate(X0, Fx, dF, alpha_x, delta, Jdelta, wrk)
    end subroutine taylor_test_system_linearization_cdp


end module LightKrylov_AbstractSystems
