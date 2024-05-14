module lightkrylov_AbstractLinops
    use lightkrylov_constants
    use lightkrylov_utils
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




    !--------------------------------------------------
    !-----     Definition of the Identity map     -----
    !--------------------------------------------------

    type, extends(abstract_linop_rsp), public :: Id_rsp
        contains
        private
        procedure, pass(self), public :: matvec => id_matvec_rsp
        procedure, pass(self), public :: rmatvec => id_matvec_rsp
    end type

    type, extends(abstract_linop_rdp), public :: Id_rdp
        contains
        private
        procedure, pass(self), public :: matvec => id_matvec_rdp
        procedure, pass(self), public :: rmatvec => id_matvec_rdp
    end type

    type, extends(abstract_linop_csp), public :: Id_csp
        contains
        private
        procedure, pass(self), public :: matvec => id_matvec_csp
        procedure, pass(self), public :: rmatvec => id_matvec_csp
    end type

    type, extends(abstract_linop_cdp), public :: Id_cdp
        contains
        private
        procedure, pass(self), public :: matvec => id_matvec_cdp
        procedure, pass(self), public :: rmatvec => id_matvec_cdp
    end type


    !----------------------------------------------
    !-----     Definition of scaled linop     -----
    !----------------------------------------------

    type, extends(abstract_linop_rsp), public :: scaled_linop_rsp
        !! Defines a scaled linear operator \( \mathbf{B} = \sigma \mathbf{A} \) with \( \mathbf{A} \) a real-valued operator and \( \sigma \in \mathbb{R} \).
        class(abstract_linop_rsp), allocatable :: A
        !! Base linear operator to be scaled.
        real(sp) :: sigma
        !! Scaling factor.
    contains
        private
        procedure, pass(self), public :: matvec => scaled_matvec_rsp
        procedure, pass(self), public :: rmatvec => scaled_rmatvec_rsp
    end type
    type, extends(abstract_linop_rdp), public :: scaled_linop_rdp
        !! Defines a scaled linear operator \( \mathbf{B} = \sigma \mathbf{A} \) with \( \mathbf{A} \) a real-valued operator and \( \sigma \in \mathbb{R} \).
        class(abstract_linop_rdp), allocatable :: A
        !! Base linear operator to be scaled.
        real(dp) :: sigma
        !! Scaling factor.
    contains
        private
        procedure, pass(self), public :: matvec => scaled_matvec_rdp
        procedure, pass(self), public :: rmatvec => scaled_rmatvec_rdp
    end type
    type, extends(abstract_linop_csp), public :: scaled_linop_csp
        !! Defines a scaled linear operator \( \mathbf{B} = \sigma \mathbf{A} \) with \( \mathbf{A} \) a complex-valued operator and \( \sigma \in \mathbb{C} \).
        class(abstract_linop_csp), allocatable :: A
        !! Base linear operator to be scaled.
        complex(sp) :: sigma
        !! Scaling factor.
    contains
        private
        procedure, pass(self), public :: matvec => scaled_matvec_csp
        procedure, pass(self), public :: rmatvec => scaled_rmatvec_csp
    end type
    type, extends(abstract_linop_cdp), public :: scaled_linop_cdp
        !! Defines a scaled linear operator \( \mathbf{B} = \sigma \mathbf{A} \) with \( \mathbf{A} \) a complex-valued operator and \( \sigma \in \mathbb{C} \).
        class(abstract_linop_cdp), allocatable :: A
        !! Base linear operator to be scaled.
        complex(dp) :: sigma
        !! Scaling factor.
    contains
        private
        procedure, pass(self), public :: matvec => scaled_matvec_cdp
        procedure, pass(self), public :: rmatvec => scaled_rmatvec_cdp
    end type

    !------------------------------------------------
    !-----     Definition of axpby operator     -----
    !------------------------------------------------

    type, extends(abstract_linop_rsp), public :: axpby_linop_rsp
        class(abstract_linop_rsp), allocatable :: A, B
        real(sp) :: alpha, beta
        logical :: transA = .false., transB = .false.
    contains
        private
        procedure, pass(self), public :: matvec => axpby_matvec_rsp
        procedure, pass(self), public :: rmatvec => axpby_rmatvec_rsp
    end type
    type, extends(abstract_linop_rdp), public :: axpby_linop_rdp
        class(abstract_linop_rdp), allocatable :: A, B
        real(dp) :: alpha, beta
        logical :: transA = .false., transB = .false.
    contains
        private
        procedure, pass(self), public :: matvec => axpby_matvec_rdp
        procedure, pass(self), public :: rmatvec => axpby_rmatvec_rdp
    end type
    type, extends(abstract_linop_csp), public :: axpby_linop_csp
        class(abstract_linop_csp), allocatable :: A, B
        complex(sp) :: alpha, beta
        logical :: transA = .false., transB = .false.
    contains
        private
        procedure, pass(self), public :: matvec => axpby_matvec_csp
        procedure, pass(self), public :: rmatvec => axpby_rmatvec_csp
    end type
    type, extends(abstract_linop_cdp), public :: axpby_linop_cdp
        class(abstract_linop_cdp), allocatable :: A, B
        complex(dp) :: alpha, beta
        logical :: transA = .false., transB = .false.
    contains
        private
        procedure, pass(self), public :: matvec => axpby_matvec_cdp
        procedure, pass(self), public :: rmatvec => axpby_rmatvec_cdp
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

    subroutine id_matvec_rsp(self, vec_in, vec_out)
        class(Id_rsp), intent(in) :: self
        class(abstract_vector_rsp), intent(in) :: vec_in
        class(abstract_vector_rsp), intent(out) :: vec_out
        call vec_out%axpby(zero_rsp, vec_in, one_rsp)
        return
    end subroutine id_matvec_rsp
    subroutine id_matvec_rdp(self, vec_in, vec_out)
        class(Id_rdp), intent(in) :: self
        class(abstract_vector_rdp), intent(in) :: vec_in
        class(abstract_vector_rdp), intent(out) :: vec_out
        call vec_out%axpby(zero_rdp, vec_in, one_rdp)
        return
    end subroutine id_matvec_rdp
    subroutine id_matvec_csp(self, vec_in, vec_out)
        class(Id_csp), intent(in) :: self
        class(abstract_vector_csp), intent(in) :: vec_in
        class(abstract_vector_csp), intent(out) :: vec_out
        call vec_out%axpby(zero_csp, vec_in, one_csp)
        return
    end subroutine id_matvec_csp
    subroutine id_matvec_cdp(self, vec_in, vec_out)
        class(Id_cdp), intent(in) :: self
        class(abstract_vector_cdp), intent(in) :: vec_in
        class(abstract_vector_cdp), intent(out) :: vec_out
        call vec_out%axpby(zero_cdp, vec_in, one_cdp)
        return
    end subroutine id_matvec_cdp

    subroutine scaled_matvec_rsp(self, vec_in, vec_out)
        class(scaled_linop_rsp), intent(in) :: self
        class(abstract_vector_rsp), intent(in) :: vec_in
        class(abstract_vector_rsp), intent(out) :: vec_out

        call self%A%matvec(vec_in, vec_out) ; call vec_out%scal(self%sigma)
        return
    end subroutine scaled_matvec_rsp

    subroutine scaled_rmatvec_rsp(self, vec_in, vec_out)
        class(scaled_linop_rsp), intent(in) :: self
        class(abstract_vector_rsp), intent(in) :: vec_in
        class(abstract_vector_rsp), intent(out) :: vec_out

        call self%A%rmatvec(vec_in, vec_out) ; call vec_out%scal(self%sigma)
        return
    end subroutine scaled_rmatvec_rsp
    subroutine scaled_matvec_rdp(self, vec_in, vec_out)
        class(scaled_linop_rdp), intent(in) :: self
        class(abstract_vector_rdp), intent(in) :: vec_in
        class(abstract_vector_rdp), intent(out) :: vec_out

        call self%A%matvec(vec_in, vec_out) ; call vec_out%scal(self%sigma)
        return
    end subroutine scaled_matvec_rdp

    subroutine scaled_rmatvec_rdp(self, vec_in, vec_out)
        class(scaled_linop_rdp), intent(in) :: self
        class(abstract_vector_rdp), intent(in) :: vec_in
        class(abstract_vector_rdp), intent(out) :: vec_out

        call self%A%rmatvec(vec_in, vec_out) ; call vec_out%scal(self%sigma)
        return
    end subroutine scaled_rmatvec_rdp
    subroutine scaled_matvec_csp(self, vec_in, vec_out)
        class(scaled_linop_csp), intent(in) :: self
        class(abstract_vector_csp), intent(in) :: vec_in
        class(abstract_vector_csp), intent(out) :: vec_out

        call self%A%matvec(vec_in, vec_out) ; call vec_out%scal(self%sigma)
        return
    end subroutine scaled_matvec_csp

    subroutine scaled_rmatvec_csp(self, vec_in, vec_out)
        class(scaled_linop_csp), intent(in) :: self
        class(abstract_vector_csp), intent(in) :: vec_in
        class(abstract_vector_csp), intent(out) :: vec_out

        call self%A%rmatvec(vec_in, vec_out) ; call vec_out%scal(self%sigma)
        return
    end subroutine scaled_rmatvec_csp
    subroutine scaled_matvec_cdp(self, vec_in, vec_out)
        class(scaled_linop_cdp), intent(in) :: self
        class(abstract_vector_cdp), intent(in) :: vec_in
        class(abstract_vector_cdp), intent(out) :: vec_out

        call self%A%matvec(vec_in, vec_out) ; call vec_out%scal(self%sigma)
        return
    end subroutine scaled_matvec_cdp

    subroutine scaled_rmatvec_cdp(self, vec_in, vec_out)
        class(scaled_linop_cdp), intent(in) :: self
        class(abstract_vector_cdp), intent(in) :: vec_in
        class(abstract_vector_cdp), intent(out) :: vec_out

        call self%A%rmatvec(vec_in, vec_out) ; call vec_out%scal(self%sigma)
        return
    end subroutine scaled_rmatvec_cdp

    subroutine axpby_matvec_rsp(self, vec_in, vec_out)
        class(axpby_linop_rsp), intent(in) :: self
        class(abstract_vector_rsp), intent(in) :: vec_in
        class(abstract_vector_rsp), intent(out) :: vec_out

        ! Working array.
        class(abstract_vector_rsp), allocatable :: wrk

        ! Allocate working array.
        allocate(wrk, source=vec_in) ; call wrk%zero()

        ! w = A @ x
        if (self%transA) then
            call self%A%rmatvec(vec_in, wrk)
        else
            call self%A%matvec(vec_in, wrk)
        endif

        ! y = B @ x
        if (self%transB) then
            call self%B%rmatvec(vec_in, vec_out)
        else
            call self%B%matvec(vec_in, vec_out)
        endif

        ! y = alpha*w + beta*y
        call vec_out%axpby(self%beta, wrk, self%alpha)

        return
    end subroutine axpby_matvec_rsp

    subroutine axpby_rmatvec_rsp(self, vec_in, vec_out)
        class(axpby_linop_rsp), intent(in) :: self
        class(abstract_vector_rsp), intent(in) :: vec_in
        class(abstract_vector_rsp), intent(out) :: vec_out

        ! Working array.
        class(abstract_vector_rsp), allocatable :: wrk

        ! Allocate working array.
        allocate(wrk, source=vec_in) ; call wrk%zero()

        ! w = A @ x
        if (self%transA) then
            call self%A%matvec(vec_in, wrk)
        else
            call self%A%rmatvec(vec_in, wrk)
        endif

        ! y = B @ x
        if (self%transB) then
            call self%B%matvec(vec_in, vec_out)
        else
            call self%B%rmatvec(vec_in, vec_out)
        endif

        ! y = alpha*w + beta*y
        call vec_out%axpby(self%beta, wrk, self%alpha)

        return
    end subroutine axpby_rmatvec_rsp

    subroutine axpby_matvec_rdp(self, vec_in, vec_out)
        class(axpby_linop_rdp), intent(in) :: self
        class(abstract_vector_rdp), intent(in) :: vec_in
        class(abstract_vector_rdp), intent(out) :: vec_out

        ! Working array.
        class(abstract_vector_rdp), allocatable :: wrk

        ! Allocate working array.
        allocate(wrk, source=vec_in) ; call wrk%zero()

        ! w = A @ x
        if (self%transA) then
            call self%A%rmatvec(vec_in, wrk)
        else
            call self%A%matvec(vec_in, wrk)
        endif

        ! y = B @ x
        if (self%transB) then
            call self%B%rmatvec(vec_in, vec_out)
        else
            call self%B%matvec(vec_in, vec_out)
        endif

        ! y = alpha*w + beta*y
        call vec_out%axpby(self%beta, wrk, self%alpha)

        return
    end subroutine axpby_matvec_rdp

    subroutine axpby_rmatvec_rdp(self, vec_in, vec_out)
        class(axpby_linop_rdp), intent(in) :: self
        class(abstract_vector_rdp), intent(in) :: vec_in
        class(abstract_vector_rdp), intent(out) :: vec_out

        ! Working array.
        class(abstract_vector_rdp), allocatable :: wrk

        ! Allocate working array.
        allocate(wrk, source=vec_in) ; call wrk%zero()

        ! w = A @ x
        if (self%transA) then
            call self%A%matvec(vec_in, wrk)
        else
            call self%A%rmatvec(vec_in, wrk)
        endif

        ! y = B @ x
        if (self%transB) then
            call self%B%matvec(vec_in, vec_out)
        else
            call self%B%rmatvec(vec_in, vec_out)
        endif

        ! y = alpha*w + beta*y
        call vec_out%axpby(self%beta, wrk, self%alpha)

        return
    end subroutine axpby_rmatvec_rdp

    subroutine axpby_matvec_csp(self, vec_in, vec_out)
        class(axpby_linop_csp), intent(in) :: self
        class(abstract_vector_csp), intent(in) :: vec_in
        class(abstract_vector_csp), intent(out) :: vec_out

        ! Working array.
        class(abstract_vector_csp), allocatable :: wrk

        ! Allocate working array.
        allocate(wrk, source=vec_in) ; call wrk%zero()

        ! w = A @ x
        if (self%transA) then
            call self%A%rmatvec(vec_in, wrk)
        else
            call self%A%matvec(vec_in, wrk)
        endif

        ! y = B @ x
        if (self%transB) then
            call self%B%rmatvec(vec_in, vec_out)
        else
            call self%B%matvec(vec_in, vec_out)
        endif

        ! y = alpha*w + beta*y
        call vec_out%axpby(self%beta, wrk, self%alpha)

        return
    end subroutine axpby_matvec_csp

    subroutine axpby_rmatvec_csp(self, vec_in, vec_out)
        class(axpby_linop_csp), intent(in) :: self
        class(abstract_vector_csp), intent(in) :: vec_in
        class(abstract_vector_csp), intent(out) :: vec_out

        ! Working array.
        class(abstract_vector_csp), allocatable :: wrk

        ! Allocate working array.
        allocate(wrk, source=vec_in) ; call wrk%zero()

        ! w = A @ x
        if (self%transA) then
            call self%A%matvec(vec_in, wrk)
        else
            call self%A%rmatvec(vec_in, wrk)
        endif

        ! y = B @ x
        if (self%transB) then
            call self%B%matvec(vec_in, vec_out)
        else
            call self%B%rmatvec(vec_in, vec_out)
        endif

        ! y = alpha*w + beta*y
        call vec_out%axpby(self%beta, wrk, self%alpha)

        return
    end subroutine axpby_rmatvec_csp

    subroutine axpby_matvec_cdp(self, vec_in, vec_out)
        class(axpby_linop_cdp), intent(in) :: self
        class(abstract_vector_cdp), intent(in) :: vec_in
        class(abstract_vector_cdp), intent(out) :: vec_out

        ! Working array.
        class(abstract_vector_cdp), allocatable :: wrk

        ! Allocate working array.
        allocate(wrk, source=vec_in) ; call wrk%zero()

        ! w = A @ x
        if (self%transA) then
            call self%A%rmatvec(vec_in, wrk)
        else
            call self%A%matvec(vec_in, wrk)
        endif

        ! y = B @ x
        if (self%transB) then
            call self%B%rmatvec(vec_in, vec_out)
        else
            call self%B%matvec(vec_in, vec_out)
        endif

        ! y = alpha*w + beta*y
        call vec_out%axpby(self%beta, wrk, self%alpha)

        return
    end subroutine axpby_matvec_cdp

    subroutine axpby_rmatvec_cdp(self, vec_in, vec_out)
        class(axpby_linop_cdp), intent(in) :: self
        class(abstract_vector_cdp), intent(in) :: vec_in
        class(abstract_vector_cdp), intent(out) :: vec_out

        ! Working array.
        class(abstract_vector_cdp), allocatable :: wrk

        ! Allocate working array.
        allocate(wrk, source=vec_in) ; call wrk%zero()

        ! w = A @ x
        if (self%transA) then
            call self%A%matvec(vec_in, wrk)
        else
            call self%A%rmatvec(vec_in, wrk)
        endif

        ! y = B @ x
        if (self%transB) then
            call self%B%matvec(vec_in, vec_out)
        else
            call self%B%rmatvec(vec_in, vec_out)
        endif

        ! y = alpha*w + beta*y
        call vec_out%axpby(self%beta, wrk, self%alpha)

        return
    end subroutine axpby_rmatvec_cdp


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
