module LightKrylov_AbstractSystems
    !!  This module provides the abstract types necessary to define an algebraic system of
    !!  nonlinear equations to be solved using the Newton method.
    use LightKrylov_Constants
    use LightKrylov_AbstractVectors
    use LightKrylov_AbstractLinops
    implicit none
    private

    character(len=128), parameter :: this_module = 'LightKrylov_AbstractSystems'

    ! Base type for abstract systems.
    type, abstract, public :: abstract_system
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
        procedure(abstract_eval_rsp), pass(self), deferred, public :: eval
        !! Procedure to evaluate the system response \( \mathbf{Y} = \mathbf{F}(\mathbf{X}) \).
    end type

    abstract interface
        subroutine abstract_eval_rsp(self, vec_in, vec_out, atol)
            !! Interface for the evaluation of the system response.
            use LightKrylov_AbstractVectors
            import abstract_system_rsp, sp
            class(abstract_system_rsp), intent(in)  :: self
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
        procedure(abstract_eval_rdp), pass(self), deferred, public :: eval
        !! Procedure to evaluate the system response \( \mathbf{Y} = \mathbf{F}(\mathbf{X}) \).
    end type

    abstract interface
        subroutine abstract_eval_rdp(self, vec_in, vec_out, atol)
            !! Interface for the evaluation of the system response.
            use LightKrylov_AbstractVectors
            import abstract_system_rdp, dp
            class(abstract_system_rdp), intent(in)  :: self
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
        procedure(abstract_eval_csp), pass(self), deferred, public :: eval
        !! Procedure to evaluate the system response \( \mathbf{Y} = \mathbf{F}(\mathbf{X}) \).
    end type

    abstract interface
        subroutine abstract_eval_csp(self, vec_in, vec_out, atol)
            !! Interface for the evaluation of the system response.
            use LightKrylov_AbstractVectors
            import abstract_system_csp, sp
            class(abstract_system_csp), intent(in)  :: self
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
        procedure(abstract_eval_cdp), pass(self), deferred, public :: eval
        !! Procedure to evaluate the system response \( \mathbf{Y} = \mathbf{F}(\mathbf{X}) \).
    end type

    abstract interface
        subroutine abstract_eval_cdp(self, vec_in, vec_out, atol)
            !! Interface for the evaluation of the system response.
            use LightKrylov_AbstractVectors
            import abstract_system_cdp, dp
            class(abstract_system_cdp), intent(in)  :: self
            !! System
            class(abstract_vector_cdp), intent(in)  :: vec_in
            !! State
            class(abstract_vector_cdp), intent(out) :: vec_out
            !! Response
            real(dp),                   intent(in)  :: atol
            !! Solver tolerance
        end subroutine abstract_eval_cdp
    end interface

end module LightKrylov_AbstractSystems
