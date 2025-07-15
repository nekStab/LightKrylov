module LightKrylov_AbstractSystems
    !!  This module provides the abstract types necessary to define an algebraic system of
    !!  nonlinear equations to be solved using the Newton method.
    use stdlib_optval, only: optval
    use LightKrylov_Logger
    use LightKrylov_Constants
    use LightKrylov_Timer_Utils, only: lightkrylov_timer
    use LightKrylov_Timing, only: time_lightkrylov
    use LightKrylov_AbstractVectors
    use LightKrylov_AbstractLinops
    implicit none
    private

    character(len=*), parameter :: this_module      = 'LK_Systems'
    character(len=*), parameter :: this_module_long = 'LK_AbstractSystems'

    ! Base type for abstract systems.
    type, abstract, public :: abstract_system
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
        !! Abstract type for the local linearization of the system around the state X
        class(abstract_vector_rsp), allocatable :: X
        !! System state around which the equatons are linearized.
    contains
    end type

    ! Abstract system for kind=sp.
    type, abstract, extends(abstract_system), public :: abstract_system_rsp
        !! System for Newton fixed-point iteration via the Jacobian
        class(abstract_jacobian_linop_rsp), allocatable :: jacobian
        !! System Jacobian \( \left. \frac{\partial \mathbf{F}}{\partial \mathbf{X}} \right|_{X^*} \).
    contains
        private
        procedure(abstract_eval_rsp), pass(self), deferred, public :: response
        !! Procedure to evaluate the system response \( \mathbf{Y} = \mathbf{F}(\mathbf{X}) \).
        ! Wrapper including counter increment
        procedure, pass(self), public :: eval => eval_rsp
        !! Wrapper for response including the counter increment
    end type

    abstract interface
        subroutine abstract_eval_rsp(self, vec_in, vec_out, atol)
            !! Interface for the evaluation of the system response.
            use LightKrylov_AbstractVectors
            import abstract_system_rsp, sp
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
        !! Abstract type for the local linearization of the system around the state X
        class(abstract_vector_rdp), allocatable :: X
        !! System state around which the equatons are linearized.
    contains
    end type

    ! Abstract system for kind=dp.
    type, abstract, extends(abstract_system), public :: abstract_system_rdp
        !! System for Newton fixed-point iteration via the Jacobian
        class(abstract_jacobian_linop_rdp), allocatable :: jacobian
        !! System Jacobian \( \left. \frac{\partial \mathbf{F}}{\partial \mathbf{X}} \right|_{X^*} \).
    contains
        private
        procedure(abstract_eval_rdp), pass(self), deferred, public :: response
        !! Procedure to evaluate the system response \( \mathbf{Y} = \mathbf{F}(\mathbf{X}) \).
        ! Wrapper including counter increment
        procedure, pass(self), public :: eval => eval_rdp
        !! Wrapper for response including the counter increment
    end type

    abstract interface
        subroutine abstract_eval_rdp(self, vec_in, vec_out, atol)
            !! Interface for the evaluation of the system response.
            use LightKrylov_AbstractVectors
            import abstract_system_rdp, dp
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
        !! Abstract type for the local linearization of the system around the state X
        class(abstract_vector_csp), allocatable :: X
        !! System state around which the equatons are linearized.
    contains
    end type

    ! Abstract system for kind=sp.
    type, abstract, extends(abstract_system), public :: abstract_system_csp
        !! System for Newton fixed-point iteration via the Jacobian
        class(abstract_jacobian_linop_csp), allocatable :: jacobian
        !! System Jacobian \( \left. \frac{\partial \mathbf{F}}{\partial \mathbf{X}} \right|_{X^*} \).
    contains
        private
        procedure(abstract_eval_csp), pass(self), deferred, public :: response
        !! Procedure to evaluate the system response \( \mathbf{Y} = \mathbf{F}(\mathbf{X}) \).
        ! Wrapper including counter increment
        procedure, pass(self), public :: eval => eval_csp
        !! Wrapper for response including the counter increment
    end type

    abstract interface
        subroutine abstract_eval_csp(self, vec_in, vec_out, atol)
            !! Interface for the evaluation of the system response.
            use LightKrylov_AbstractVectors
            import abstract_system_csp, sp
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
        !! Abstract type for the local linearization of the system around the state X
        class(abstract_vector_cdp), allocatable :: X
        !! System state around which the equatons are linearized.
    contains
    end type

    ! Abstract system for kind=dp.
    type, abstract, extends(abstract_system), public :: abstract_system_cdp
        !! System for Newton fixed-point iteration via the Jacobian
        class(abstract_jacobian_linop_cdp), allocatable :: jacobian
        !! System Jacobian \( \left. \frac{\partial \mathbf{F}}{\partial \mathbf{X}} \right|_{X^*} \).
    contains
        private
        procedure(abstract_eval_cdp), pass(self), deferred, public :: response
        !! Procedure to evaluate the system response \( \mathbf{Y} = \mathbf{F}(\mathbf{X}) \).
        ! Wrapper including counter increment
        procedure, pass(self), public :: eval => eval_cdp
        !! Wrapper for response including the counter increment
    end type

    abstract interface
        subroutine abstract_eval_cdp(self, vec_in, vec_out, atol)
            !! Interface for the evaluation of the system response.
            use LightKrylov_AbstractVectors
            import abstract_system_cdp, dp
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

contains

    !---------------------------------------------------------------
    !-----     Getter/Setter routines for abstract_systems     -----
    !---------------------------------------------------------------
 
    pure integer function get_eval_counter(self) result(count)
      !! Getter function for the number of eval calls
      class(abstract_system), intent(in) :: self
      count = self%eval_counter
    end function get_eval_counter

    subroutine reset_eval_counter(self, procedure, counter, reset_timer, soft_reset, clean_timer)
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
      !! Print the current timing data for the system evaluation
      !! Note: Wrapper of the corresponding routine from lightkrylov_timer
      class(abstract_system), intent(inout) :: self
      call self%eval_timer%print_info()
    end subroutine print_timer_info

    subroutine reset_eval_timer(self, soft, clean)
      !! Setter routine to reset the system evaluation timer
      !! Note: Wrapper of the corresponding routine from lightkrylov_timer
      class(abstract_system), intent(inout) :: self
      logical, optional, intent(in) :: soft
      logical, optional, intent(in) :: clean
      call self%eval_timer%reset(soft, clean, verbose=.true.)
    end subroutine reset_eval_timer

    subroutine finalize_eval_timer(self)
      !! Finalize the system evaluation timer and print summary
      !! Note: Wrapper of the corresponding routine from lightkrylov_timer
      class(abstract_system), intent(inout) :: self
      call self%eval_timer%finalize()
    end subroutine finalize_eval_timer

    !---------------------------------------------------------------------
    !-----     Wrapper for system response to increment counters     -----
    !---------------------------------------------------------------------

    subroutine eval_rsp(self, vec_in, vec_out, atol)
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

end module LightKrylov_AbstractSystems
