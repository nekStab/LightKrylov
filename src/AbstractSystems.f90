module LightKrylov_AbstractSystems
    !!  This module provides the abstract types necessary to define an algebraic system of
    !!  nonlinear equations to be solved using the Newton method.
    use stdlib_optval, only: optval
    use LightKrylov_Logger
    use LightKrylov_Constants
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
    contains
        procedure, pass(self), public :: get_eval_counter
        !! Return eval counter value
        procedure, pass(self), public :: reset_eval_counter
        !! Reset eval counter
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

    subroutine reset_eval_counter(self, counter, ifprint)
      class(abstract_system), intent(inout) :: self
      integer, optional, intent(in) :: counter
      !! optional flag to reset to an integer other than zero.
      logical, optional, intent(in) :: ifprint
      !! optional flag to print the number of evals to log prior to resetting.
      ! internals
      integer :: counter_
      logical :: ifprint_
      character(len=128) :: msg
      counter_ = optval(counter, 0)
      ifprint_ = optval(ifprint, .false.)
      if (ifprint_) then
         write(msg,'(A,I0,A,I0,A)') 'Total number of evals: ', self%get_eval_counter(), '. Resetting counter to ', counter_, '.'
         call logger%log_message(msg, module=this_module, procedure='reset_eval_counter')
      end if
      self%eval_counter = counter_
      return
    end subroutine reset_eval_counter

    !---------------------------------------------------------------------
    !-----     Wrapper for system response to increment counters     -----
    !---------------------------------------------------------------------

    subroutine eval_rsp(self, vec_in, vec_out, atol)
        class(abstract_system_rsp), intent(inout) :: self
        class(abstract_vector_rsp), intent(in)    :: vec_in
        class(abstract_vector_rsp), intent(out)   :: vec_out
        real(sp),                             intent(in)    :: atol
        self%eval_counter = self%eval_counter + 1
        call self%response(vec_in, vec_out, atol)
        return
    end subroutine eval_rsp
    subroutine eval_rdp(self, vec_in, vec_out, atol)
        class(abstract_system_rdp), intent(inout) :: self
        class(abstract_vector_rdp), intent(in)    :: vec_in
        class(abstract_vector_rdp), intent(out)   :: vec_out
        real(dp),                             intent(in)    :: atol
        self%eval_counter = self%eval_counter + 1
        call self%response(vec_in, vec_out, atol)
        return
    end subroutine eval_rdp
    subroutine eval_csp(self, vec_in, vec_out, atol)
        class(abstract_system_csp), intent(inout) :: self
        class(abstract_vector_csp), intent(in)    :: vec_in
        class(abstract_vector_csp), intent(out)   :: vec_out
        real(sp),                             intent(in)    :: atol
        self%eval_counter = self%eval_counter + 1
        call self%response(vec_in, vec_out, atol)
        return
    end subroutine eval_csp
    subroutine eval_cdp(self, vec_in, vec_out, atol)
        class(abstract_system_cdp), intent(inout) :: self
        class(abstract_vector_cdp), intent(in)    :: vec_in
        class(abstract_vector_cdp), intent(out)   :: vec_out
        real(dp),                             intent(in)    :: atol
        self%eval_counter = self%eval_counter + 1
        call self%response(vec_in, vec_out, atol)
        return
    end subroutine eval_cdp

end module LightKrylov_AbstractSystems
