module LightKrylov_AbstractLinops
    !!  This module provides the base classes `abtract_linop_rsp`, `abstract_linop_rdp`,
    !!  `abstract_linop_csp` and `abstract_linop_cdp` which can be used to define your own
    !!  linear operators. To do so, you simply need to provide two type-bound procedures:
    !!  
    !!  - `matvec(self, vec_in, vec_out)` : Computes the matrix-vector product.
    !!  - `rmatvec(self, vec_in, vec_out)`: Computes the transpose matrix-vector product.
    !!
    !!  It also provides extended types to define the identity operator, symmetric linear
    !!  operators, scalar-multiplication of a linear multiplication, as well as addition
    !!  of two linear operators.
    use LightKrylov_Constants
    use LightKrylov_Utils
    use LightKrylov_AbstractVectors
    implicit none
    private

    character(len=128), parameter :: this_module = 'Lightkrylov_AbstractLinops'

    type, abstract, public :: abstract_linop
        !!  Base type to define an abstract linear operator. All other types defined in
        !!  `LightKrylov` derive from this fundamental one.
        !!
        !!  @warning
        !!  Users should not extend this abstract class to define their own types.
        !!  @endwarning
    end type abstract_linop

    !------------------------------------------------------------------------------
    !-----     Definition of an abstract real(sp) operator with kind=sp     -----
    !------------------------------------------------------------------------------
    type, abstract, extends(abstract_linop), public :: abstract_linop_rsp
        !! Base type to extend in order to define a real(sp)-valued linear operator.
    contains
        private
        procedure(abstract_matvec_rsp), pass(self), deferred, public :: matvec
        !! Procedure to compute the matrix-vector product \( \mathbf{y} = \mathbf{Ax} \).
        procedure(abstract_matvec_rsp), pass(self), deferred, public :: rmatvec
        !! Procedure to compute the reversed matrix-vector product \( \mathbf{y} = \mathbf{A}^H \mathbf{x} \).
    end type

    abstract interface
        subroutine abstract_matvec_rsp(self, vec_in, vec_out)
            !! Interface for the matrix-vector product.
            use LightKrylov_AbstractVectors
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
        !! Utility type to define an adjoint linear operator. The definition of `matvec` and `rmatvec`
        !! are directly inherited from those used to define `A`. Note that this utility does not
        !! compute the adjoint for you. It simply provides a utility to define a new operator
        !! with `matvec` and `rmatvec` being switched.
        class(abstract_linop_rsp), allocatable :: A
        !! Linear operator whose adjoint needs to be defined.
    contains
        private
        procedure, pass(self), public :: matvec => adjoint_matvec_rsp
        procedure, pass(self), public :: rmatvec => adjoint_rmatvec_rsp
    end type



    !------------------------------------------------------------------------------
    !-----     Definition of an abstract real(dp) operator with kind=dp     -----
    !------------------------------------------------------------------------------
    type, abstract, extends(abstract_linop), public :: abstract_linop_rdp
        !! Base type to extend in order to define a real(dp)-valued linear operator.
    contains
        private
        procedure(abstract_matvec_rdp), pass(self), deferred, public :: matvec
        !! Procedure to compute the matrix-vector product \( \mathbf{y} = \mathbf{Ax} \).
        procedure(abstract_matvec_rdp), pass(self), deferred, public :: rmatvec
        !! Procedure to compute the reversed matrix-vector product \( \mathbf{y} = \mathbf{A}^H \mathbf{x} \).
    end type

    abstract interface
        subroutine abstract_matvec_rdp(self, vec_in, vec_out)
            !! Interface for the matrix-vector product.
            use LightKrylov_AbstractVectors
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
        !! Utility type to define an adjoint linear operator. The definition of `matvec` and `rmatvec`
        !! are directly inherited from those used to define `A`. Note that this utility does not
        !! compute the adjoint for you. It simply provides a utility to define a new operator
        !! with `matvec` and `rmatvec` being switched.
        class(abstract_linop_rdp), allocatable :: A
        !! Linear operator whose adjoint needs to be defined.
    contains
        private
        procedure, pass(self), public :: matvec => adjoint_matvec_rdp
        procedure, pass(self), public :: rmatvec => adjoint_rmatvec_rdp
    end type



    !------------------------------------------------------------------------------
    !-----     Definition of an abstract complex(sp) operator with kind=sp     -----
    !------------------------------------------------------------------------------
    type, abstract, extends(abstract_linop), public :: abstract_linop_csp
        !! Base type to extend in order to define a complex(sp)-valued linear operator.
    contains
        private
        procedure(abstract_matvec_csp), pass(self), deferred, public :: matvec
        !! Procedure to compute the matrix-vector product \( \mathbf{y} = \mathbf{Ax} \).
        procedure(abstract_matvec_csp), pass(self), deferred, public :: rmatvec
        !! Procedure to compute the reversed matrix-vector product \( \mathbf{y} = \mathbf{A}^H \mathbf{x} \).
    end type

    abstract interface
        subroutine abstract_matvec_csp(self, vec_in, vec_out)
            !! Interface for the matrix-vector product.
            use LightKrylov_AbstractVectors
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
        !! Utility type to define an adjoint linear operator. The definition of `matvec` and `rmatvec`
        !! are directly inherited from those used to define `A`. Note that this utility does not
        !! compute the adjoint for you. It simply provides a utility to define a new operator
        !! with `matvec` and `rmatvec` being switched.
        class(abstract_linop_csp), allocatable :: A
        !! Linear operator whose adjoint needs to be defined.
    contains
        private
        procedure, pass(self), public :: matvec => adjoint_matvec_csp
        procedure, pass(self), public :: rmatvec => adjoint_rmatvec_csp
    end type



    !------------------------------------------------------------------------------
    !-----     Definition of an abstract complex(dp) operator with kind=dp     -----
    !------------------------------------------------------------------------------
    type, abstract, extends(abstract_linop), public :: abstract_linop_cdp
        !! Base type to extend in order to define a complex(dp)-valued linear operator.
    contains
        private
        procedure(abstract_matvec_cdp), pass(self), deferred, public :: matvec
        !! Procedure to compute the matrix-vector product \( \mathbf{y} = \mathbf{Ax} \).
        procedure(abstract_matvec_cdp), pass(self), deferred, public :: rmatvec
        !! Procedure to compute the reversed matrix-vector product \( \mathbf{y} = \mathbf{A}^H \mathbf{x} \).
    end type

    abstract interface
        subroutine abstract_matvec_cdp(self, vec_in, vec_out)
            !! Interface for the matrix-vector product.
            use LightKrylov_AbstractVectors
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
        !! Utility type to define an adjoint linear operator. The definition of `matvec` and `rmatvec`
        !! are directly inherited from those used to define `A`. Note that this utility does not
        !! compute the adjoint for you. It simply provides a utility to define a new operator
        !! with `matvec` and `rmatvec` being switched.
        class(abstract_linop_cdp), allocatable :: A
        !! Linear operator whose adjoint needs to be defined.
    contains
        private
        procedure, pass(self), public :: matvec => adjoint_matvec_cdp
        procedure, pass(self), public :: rmatvec => adjoint_rmatvec_cdp
    end type




    !--------------------------------------------------
    !-----     Definition of the Identity map     -----
    !--------------------------------------------------

    type, extends(abstract_linop_rsp), public :: Id_rsp
        !! Utility type to define the Identity operator. Note that the type-bound procedures
        !! for `matvec` and `rmatvec` do not have to be defined by the user.
        contains
        private
        procedure, pass(self), public :: matvec => id_matvec_rsp
        procedure, pass(self), public :: rmatvec => id_matvec_rsp
    end type

    type, extends(abstract_linop_rdp), public :: Id_rdp
        !! Utility type to define the Identity operator. Note that the type-bound procedures
        !! for `matvec` and `rmatvec` do not have to be defined by the user.
        contains
        private
        procedure, pass(self), public :: matvec => id_matvec_rdp
        procedure, pass(self), public :: rmatvec => id_matvec_rdp
    end type

    type, extends(abstract_linop_csp), public :: Id_csp
        !! Utility type to define the Identity operator. Note that the type-bound procedures
        !! for `matvec` and `rmatvec` do not have to be defined by the user.
        contains
        private
        procedure, pass(self), public :: matvec => id_matvec_csp
        procedure, pass(self), public :: rmatvec => id_matvec_csp
    end type

    type, extends(abstract_linop_cdp), public :: Id_cdp
        !! Utility type to define the Identity operator. Note that the type-bound procedures
        !! for `matvec` and `rmatvec` do not have to be defined by the user.
        contains
        private
        procedure, pass(self), public :: matvec => id_matvec_cdp
        procedure, pass(self), public :: rmatvec => id_matvec_cdp
    end type


    !----------------------------------------------
    !-----     Definition of scaled linop     -----
    !----------------------------------------------

    type, extends(abstract_linop_rsp), public :: scaled_linop_rsp
        !! Defines a scaled linear operator \( \mathbf{B} = \sigma \mathbf{A} \) with \( \mathbf{A} \) a real-valued operator and \( \sigma \in \mathbb{R} \). The definitions of `matvec`
        !! and `rmatvec` are directly inherited from those used to define `A` and do not have to
        !! be defined by the user.
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
        !! Defines a scaled linear operator \( \mathbf{B} = \sigma \mathbf{A} \) with \( \mathbf{A} \) a real-valued operator and \( \sigma \in \mathbb{R} \). The definitions of `matvec`
        !! and `rmatvec` are directly inherited from those used to define `A` and do not have to
        !! be defined by the user.
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
        !! Utility type to define a composite linear operator \( \mathbf{C} = \alpha \mathbf{A} + \beta \mathbf{B} \).
        !! The definitions of `matvec` and `rmatvec` are directly inherited from those used to
        !! define `A` and `B`.
        class(abstract_linop_rsp), allocatable :: A, B
        !! Underlying linear operators.
        real(sp) :: alpha, beta
        !! Scaling factors.
        logical :: transA = .false., transB = .false.
        !! Logical flag to control whether \( \mathbf{A} \) and/or \( \mathbf{B} \) need to be transposed.
    contains
        private
        procedure, pass(self), public :: matvec => axpby_matvec_rsp
        procedure, pass(self), public :: rmatvec => axpby_rmatvec_rsp
    end type
    type, extends(abstract_linop_rdp), public :: axpby_linop_rdp
        !! Utility type to define a composite linear operator \( \mathbf{C} = \alpha \mathbf{A} + \beta \mathbf{B} \).
        !! The definitions of `matvec` and `rmatvec` are directly inherited from those used to
        !! define `A` and `B`.
        class(abstract_linop_rdp), allocatable :: A, B
        !! Underlying linear operators.
        real(dp) :: alpha, beta
        !! Scaling factors.
        logical :: transA = .false., transB = .false.
        !! Logical flag to control whether \( \mathbf{A} \) and/or \( \mathbf{B} \) need to be transposed.
    contains
        private
        procedure, pass(self), public :: matvec => axpby_matvec_rdp
        procedure, pass(self), public :: rmatvec => axpby_rmatvec_rdp
    end type
    type, extends(abstract_linop_csp), public :: axpby_linop_csp
        !! Utility type to define a composite linear operator \( \mathbf{C} = \alpha \mathbf{A} + \beta \mathbf{B} \).
        !! The definitions of `matvec` and `rmatvec` are directly inherited from those used to
        !! define `A` and `B`.
        class(abstract_linop_csp), allocatable :: A, B
        !! Underlying linear operators.
        complex(sp) :: alpha, beta
        !! Scaling factors.
        logical :: transA = .false., transB = .false.
        !! Logical flag to control whether \( \mathbf{A} \) and/or \( \mathbf{B} \) need to be transposed.
    contains
        private
        procedure, pass(self), public :: matvec => axpby_matvec_csp
        procedure, pass(self), public :: rmatvec => axpby_rmatvec_csp
    end type
    type, extends(abstract_linop_cdp), public :: axpby_linop_cdp
        !! Utility type to define a composite linear operator \( \mathbf{C} = \alpha \mathbf{A} + \beta \mathbf{B} \).
        !! The definitions of `matvec` and `rmatvec` are directly inherited from those used to
        !! define `A` and `B`.
        class(abstract_linop_cdp), allocatable :: A, B
        !! Underlying linear operators.
        complex(dp) :: alpha, beta
        !! Scaling factors.
        logical :: transA = .false., transB = .false.
        !! Logical flag to control whether \( \mathbf{A} \) and/or \( \mathbf{B} \) need to be transposed.
    contains
        private
        procedure, pass(self), public :: matvec => axpby_matvec_cdp
        procedure, pass(self), public :: rmatvec => axpby_rmatvec_cdp
    end type

    !----------------------------------------------------------------
    !-----     Definition of an abstract symmetric operator     -----
    !----------------------------------------------------------------
    type, abstract, extends(abstract_linop_rsp), public :: abstract_sym_linop_rsp
        !! Abstract representation of an abstract symmetric (real valued) linear operator.
    contains
    end type



    !----------------------------------------------------------------
    !-----     Definition of an abstract symmetric operator     -----
    !----------------------------------------------------------------
    type, abstract, extends(abstract_linop_rdp), public :: abstract_sym_linop_rdp
        !! Abstract representation of an abstract symmetric (real valued) linear operator.
    contains
    end type



    !----------------------------------------------------------------------------------
    !-----     Definition of an abstract Hermitian positive definite operator     -----
    !----------------------------------------------------------------------------------
    type, abstract, extends(abstract_linop_csp), public :: abstract_hermitian_linop_csp
        !! Abstract representation of an abstract hermitian (complex-valued) linear operator.
    contains
    end type
 



    !----------------------------------------------------------------------------------
    !-----     Definition of an abstract Hermitian positive definite operator     -----
    !----------------------------------------------------------------------------------
    type, abstract, extends(abstract_linop_cdp), public :: abstract_hermitian_linop_cdp
        !! Abstract representation of an abstract hermitian (complex-valued) linear operator.
    contains
    end type
 



contains

    subroutine id_matvec_rsp(self, vec_in, vec_out)
        class(Id_rsp), intent(in) :: self
        class(abstract_vector_rsp), intent(in) :: vec_in
        class(abstract_vector_rsp), intent(out) :: vec_out
        call copy(vec_out, vec_in)
        return
    end subroutine id_matvec_rsp
    subroutine id_matvec_rdp(self, vec_in, vec_out)
        class(Id_rdp), intent(in) :: self
        class(abstract_vector_rdp), intent(in) :: vec_in
        class(abstract_vector_rdp), intent(out) :: vec_out
        call copy(vec_out, vec_in)
        return
    end subroutine id_matvec_rdp
    subroutine id_matvec_csp(self, vec_in, vec_out)
        class(Id_csp), intent(in) :: self
        class(abstract_vector_csp), intent(in) :: vec_in
        class(abstract_vector_csp), intent(out) :: vec_out
        call copy(vec_out, vec_in)
        return
    end subroutine id_matvec_csp
    subroutine id_matvec_cdp(self, vec_in, vec_out)
        class(Id_cdp), intent(in) :: self
        class(abstract_vector_cdp), intent(in) :: vec_in
        class(abstract_vector_cdp), intent(out) :: vec_out
        call copy(vec_out, vec_in)
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


end module LightKrylov_AbstractLinops
