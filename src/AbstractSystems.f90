module LightKrylov_AbstractSystems
    use LightKrylov_AbstractVectors
    use LightKrylov_AbstractLinops
    implicit none
    private

    character*128, parameter :: this_module = 'LightKrylov_AbstractSystems'

    !> General abstract type for general systems.
    type, abstract, public :: abstract_dynamical_system
    end type abstract_dynamical_system

    !-----------------------------------------------------------
    !-----     ABSTRACT GENERAL SYSTEM TYPE DEFINITION     -----
    !-----------------------------------------------------------

    !> Abstract Jacobian linop.
    type, abstract, extends(abstract_linop_rdp), public :: abstract_jacobian_linop_rdp
        class(abstract_vector_rdp), allocatable :: X
    contains
        private
        procedure(abstract_jac_matvec_rdp), pass(self), deferred, public :: matvec
        procedure(abstract_jac_matvec_rdp), pass(self), deferred, public :: rmatvec
    end type

    abstract interface
        subroutine abstract_jac_matvec_rdp(self, vec_in, vec_out)
            !! Interface for the matrix-vector product.
            use lightkrylov_AbstractVectors
            import abstract_jacobian_linop_rdp
            class(abstract_jacobian_linop_rdp) , intent(in)  :: self
            !! Linear operator \(\mathbf{A}\).
            class(abstract_vector_rdp), intent(in)  :: vec_in
            !! Vector to be multiplied by \(\mathbf{A}\).
            class(abstract_vector_rdp), intent(out) :: vec_out
            !! Result of the matrix-vector product.
        end subroutine abstract_jac_matvec_rdp
    end interface

    !> Abstract continuous system.
    type, abstract, extends(abstract_dynamical_system), public :: abstract_system_rdp
        class(abstract_jacobian_linop_rdp), allocatable :: jacobian
        !! System Jacobian \( \left. \frac{\partial \mathbf{F}}{\partial \mathbf{X}} \right|_{X^*} \).
    contains
        private
        procedure(abstract_eval_rdp), pass(self), deferred, public :: eval
        !! Procedure to evaluate the system response \( \mathbf{Y} = \mathbf{F}(\mathbf{X}) \).
    end type

    abstract interface
        subroutine abstract_eval_rdp(self, vec_in, vec_out)
            !! Interface for the evaluation of the system response.
            use LightKrylov_AbstractVectors
            import abstract_system_rdp
            class(abstract_system_rdp), intent(in)  :: self
            !! System
            class(abstract_vector_rdp), intent(in)  :: vec_in
            !! Base state
            class(abstract_vector_rdp), intent(out) :: vec_out
            !! System response
        end subroutine abstract_eval_rdp
    end interface

end module LightKrylov_AbstractSystems