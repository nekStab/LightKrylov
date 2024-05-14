module lightkrylov_AbstractLinops
    use lightkrylov_constants
    use lightkrylov_AbstractVectors
    implicit none
    private

    type, abstract, public :: abstract_linop
    contains
    private
        procedure, pass(from), public :: copy
        generic, public :: assignment(=) => copy
    end type abstract_linop





    !------------------------------------------------------------------------------
    !-----     Definition of an abstract real(sp) operator with kind=sp     -----
    !------------------------------------------------------------------------------
    type, abstract, extends(abstract_linop), public :: abstract_linop_rsp
    contains
        private
        procedure(abstract_matvec_rsp), pass(self), deferred, public :: matvec
        procedure(abstract_matvec_rsp), pass(self), deferred, public :: rmatvec
    end type

    abstract interface
        subroutine abstract_matvec_rsp(self, vec_in, vec_out)
            !! Interface for the matrix-vector product.
            use lightkrylov_AbstractVectors
            import abstract_linop_rsp
            class(abstract_linop_rsp) , intent(in)  :: self
            !! Linear operator \(\mathbf{A}\).
            class(abstract_vector_rsp), intent(in)  :: vec_in
            !! Vector to be multiplied by \(\mathbf{A}\).
            class(abstract_vector_rsp), intent(out) :: vec_out
            !! Result of the matrix-vector product.
        end subroutine abstract_matvec_rsp
    end interface

    type, extends(abstract_linop_rsp), public :: adjoint_linop_rsp
        class(abstract_linop_rsp), allocatable :: A
    contains
        private
        procedure, pass(self), public :: matvec => adjoint_matvec_rsp
        procedure, pass(self), public :: rmatvec => adjoint_rmatvec_rsp
    end type



    !------------------------------------------------------------------------------
    !-----     Definition of an abstract real(dp) operator with kind=dp     -----
    !------------------------------------------------------------------------------
    type, abstract, extends(abstract_linop), public :: abstract_linop_rdp
    contains
        private
        procedure(abstract_matvec_rdp), pass(self), deferred, public :: matvec
        procedure(abstract_matvec_rdp), pass(self), deferred, public :: rmatvec
    end type

    abstract interface
        subroutine abstract_matvec_rdp(self, vec_in, vec_out)
            !! Interface for the matrix-vector product.
            use lightkrylov_AbstractVectors
            import abstract_linop_rdp
            class(abstract_linop_rdp) , intent(in)  :: self
            !! Linear operator \(\mathbf{A}\).
            class(abstract_vector_rdp), intent(in)  :: vec_in
            !! Vector to be multiplied by \(\mathbf{A}\).
            class(abstract_vector_rdp), intent(out) :: vec_out
            !! Result of the matrix-vector product.
        end subroutine abstract_matvec_rdp
    end interface

    type, extends(abstract_linop_rdp), public :: adjoint_linop_rdp
        class(abstract_linop_rdp), allocatable :: A
    contains
        private
        procedure, pass(self), public :: matvec => adjoint_matvec_rdp
        procedure, pass(self), public :: rmatvec => adjoint_rmatvec_rdp
    end type



    !------------------------------------------------------------------------------
    !-----     Definition of an abstract complex(sp) operator with kind=sp     -----
    !------------------------------------------------------------------------------
    type, abstract, extends(abstract_linop), public :: abstract_linop_csp
    contains
        private
        procedure(abstract_matvec_csp), pass(self), deferred, public :: matvec
        procedure(abstract_matvec_csp), pass(self), deferred, public :: rmatvec
    end type

    abstract interface
        subroutine abstract_matvec_csp(self, vec_in, vec_out)
            !! Interface for the matrix-vector product.
            use lightkrylov_AbstractVectors
            import abstract_linop_csp
            class(abstract_linop_csp) , intent(in)  :: self
            !! Linear operator \(\mathbf{A}\).
            class(abstract_vector_csp), intent(in)  :: vec_in
            !! Vector to be multiplied by \(\mathbf{A}\).
            class(abstract_vector_csp), intent(out) :: vec_out
            !! Result of the matrix-vector product.
        end subroutine abstract_matvec_csp
    end interface

    type, extends(abstract_linop_csp), public :: adjoint_linop_csp
        class(abstract_linop_csp), allocatable :: A
    contains
        private
        procedure, pass(self), public :: matvec => adjoint_matvec_csp
        procedure, pass(self), public :: rmatvec => adjoint_rmatvec_csp
    end type



    !------------------------------------------------------------------------------
    !-----     Definition of an abstract complex(dp) operator with kind=dp     -----
    !------------------------------------------------------------------------------
    type, abstract, extends(abstract_linop), public :: abstract_linop_cdp
    contains
        private
        procedure(abstract_matvec_cdp), pass(self), deferred, public :: matvec
        procedure(abstract_matvec_cdp), pass(self), deferred, public :: rmatvec
    end type

    abstract interface
        subroutine abstract_matvec_cdp(self, vec_in, vec_out)
            !! Interface for the matrix-vector product.
            use lightkrylov_AbstractVectors
            import abstract_linop_cdp
            class(abstract_linop_cdp) , intent(in)  :: self
            !! Linear operator \(\mathbf{A}\).
            class(abstract_vector_cdp), intent(in)  :: vec_in
            !! Vector to be multiplied by \(\mathbf{A}\).
            class(abstract_vector_cdp), intent(out) :: vec_out
            !! Result of the matrix-vector product.
        end subroutine abstract_matvec_cdp
    end interface

    type, extends(abstract_linop_cdp), public :: adjoint_linop_cdp
        class(abstract_linop_cdp), allocatable :: A
    contains
        private
        procedure, pass(self), public :: matvec => adjoint_matvec_cdp
        procedure, pass(self), public :: rmatvec => adjoint_rmatvec_cdp
    end type




    !----------------------------------------------------------------------------------
    !-----     Definition of an abstract symmetric positive definite operator     -----
    !----------------------------------------------------------------------------------
    type, abstract, extends(abstract_linop_rsp), public :: abstract_spd_linop_rsp
    contains
    end type




    !----------------------------------------------------------------------------------
    !-----     Definition of an abstract symmetric positive definite operator     -----
    !----------------------------------------------------------------------------------
    type, abstract, extends(abstract_linop_rdp), public :: abstract_spd_linop_rdp
    contains
    end type




    !----------------------------------------------------------------------------------
    !-----     Definition of an abstract Hermitian positive definite operator     -----
    !----------------------------------------------------------------------------------
    type, abstract, extends(abstract_linop_csp), public :: abstract_hermitian_linop_csp
    contains
    end type
 




    !----------------------------------------------------------------------------------
    !-----     Definition of an abstract Hermitian positive definite operator     -----
    !----------------------------------------------------------------------------------
    type, abstract, extends(abstract_linop_cdp), public :: abstract_hermitian_linop_cdp
    contains
    end type
 




contains

    subroutine copy(out, from)
        class(abstract_linop), intent(in) :: from
        class(abstract_linop), allocatable, intent(out) :: out
        if (allocated(out)) deallocate(out)
        allocate(out, source=from)
        return
    end subroutine copy

    subroutine adjoint_matvec_rsp(self, vec_in, vec_out)
        class(adjoint_linop_rsp), intent(in) :: self
        class(abstract_vector_rsp), intent(in) :: vec_in
        class(abstract_vector_rsp), intent(out) :: vec_out

        call self%A%rmatvec(vec_in, vec_out)

        return
    end subroutine adjoint_matvec_rsp

    subroutine adjoint_rmatvec_rsp(self, vec_in, vec_out)
        class(adjoint_linop_rsp), intent(in) :: self
        class(abstract_vector_rsp), intent(in) :: vec_in
        class(abstract_vector_rsp), intent(out) :: vec_out

        call self%A%matvec(vec_in, vec_out)

        return
    end subroutine adjoint_rmatvec_rsp

    subroutine adjoint_matvec_rdp(self, vec_in, vec_out)
        class(adjoint_linop_rdp), intent(in) :: self
        class(abstract_vector_rdp), intent(in) :: vec_in
        class(abstract_vector_rdp), intent(out) :: vec_out

        call self%A%rmatvec(vec_in, vec_out)

        return
    end subroutine adjoint_matvec_rdp

    subroutine adjoint_rmatvec_rdp(self, vec_in, vec_out)
        class(adjoint_linop_rdp), intent(in) :: self
        class(abstract_vector_rdp), intent(in) :: vec_in
        class(abstract_vector_rdp), intent(out) :: vec_out

        call self%A%matvec(vec_in, vec_out)

        return
    end subroutine adjoint_rmatvec_rdp

    subroutine adjoint_matvec_csp(self, vec_in, vec_out)
        class(adjoint_linop_csp), intent(in) :: self
        class(abstract_vector_csp), intent(in) :: vec_in
        class(abstract_vector_csp), intent(out) :: vec_out

        call self%A%rmatvec(vec_in, vec_out)

        return
    end subroutine adjoint_matvec_csp

    subroutine adjoint_rmatvec_csp(self, vec_in, vec_out)
        class(adjoint_linop_csp), intent(in) :: self
        class(abstract_vector_csp), intent(in) :: vec_in
        class(abstract_vector_csp), intent(out) :: vec_out

        call self%A%matvec(vec_in, vec_out)

        return
    end subroutine adjoint_rmatvec_csp

    subroutine adjoint_matvec_cdp(self, vec_in, vec_out)
        class(adjoint_linop_cdp), intent(in) :: self
        class(abstract_vector_cdp), intent(in) :: vec_in
        class(abstract_vector_cdp), intent(out) :: vec_out

        call self%A%rmatvec(vec_in, vec_out)

        return
    end subroutine adjoint_matvec_cdp

    subroutine adjoint_rmatvec_cdp(self, vec_in, vec_out)
        class(adjoint_linop_cdp), intent(in) :: self
        class(abstract_vector_cdp), intent(in) :: vec_in
        class(abstract_vector_cdp), intent(out) :: vec_out

        call self%A%matvec(vec_in, vec_out)

        return
    end subroutine adjoint_rmatvec_cdp


end module lightkrylov_AbstractLinops
