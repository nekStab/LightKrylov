module LightKrylov_AbstractLinops
    !!  This module provides the base classes `abtract_linop_rsp`, `abstract_linop_rdp`,
    !!  `abstract_linop_csp` and `abstract_linop_cdp` which can be used to define your own
    !!  linear operators. To do so, you simply need to provide two type-bound procedures:
    !!  
    !!  - `matvec(self, vec_in, vec_out)`   :   Computes the matrix-vector product.
    !!  - `rmatvec(self, vec_in, vec_out)   :   Computes the transpose matrix-vector product.
    !!
    !!  It also provides extended types to define the identity operator, symmetric linear
    !!  operators, scalar-multiplication of a linear multiplication, as well as addition
    !!  of two linear operators.
    use stdlib_optval, only: optval
    use stdlib_linalg_blas, only: gemv
    use LightKrylov_Logger
    use LightKrylov_Constants
    use LightKrylov_Timer_Utils, only: lightkrylov_timer
    use LightKrylov_Timing, only: time_lightkrylov
    use LightKrylov_Utils
    use LightKrylov_AbstractVectors
    implicit none
    private

    character(len=*), parameter :: this_module      = 'LK_Linops'
    character(len=*), parameter :: this_module_long = 'Lightkrylov_AbstractLinops'

    type, abstract, public :: abstract_linop
        !!  Base type to define an abstract linear operator. All other types defined in
        !!  `LightKrylov` derive from this fundamental one.
        !!
        !!  @warning
        !!  Users should not extend this abstract class to define their own types.
        !!  @endwarning
        integer, private :: matvec_counter  = 0
        integer, private :: rmatvec_counter = 0
        type(lightkrylov_timer) :: matvec_timer  = lightkrylov_timer('matvec timer')
        type(lightkrylov_timer) :: rmatvec_timer = lightkrylov_timer('rmatvec timer')
    contains
        procedure, pass(self), public :: get_counter
        !! Return matvec/rmatvec counter value
        procedure, pass(self), public :: reset_counter
        !! Reset matvec/rmatvec counter
        procedure, pass(self), public :: print_timer_info
        !! Print current timing data
        procedure, pass(self), public :: reset_timer => reset_linop_timer
        !! Reset current timing data
        procedure, pass(self), public :: finalize_timer => finalize_linop_timer
        !! Finalize timers and print complete history_info
    end type abstract_linop

    !------------------------------------------------------------------------------
    !-----     Definition of an abstract real(sp) operator with kind=sp     -----
    !------------------------------------------------------------------------------
    type, abstract, extends(abstract_linop), public :: abstract_linop_rsp
        !! Base type to extend in order to define a real(sp)-valued linear operator.
    contains
        private
        ! User defined procedures
        procedure(abstract_matvec_rsp), pass(self), deferred, public :: matvec
        !! Procedure to compute the matrix-vector product \( \mathbf{y} = \mathbf{Ax} \).
        procedure(abstract_matvec_rsp), pass(self), deferred, public :: rmatvec
        !! Procedure to compute the reversed matrix-vector product \( \mathbf{y} = \mathbf{A}^H \mathbf{x} \).
        ! Wrappers including counter increment
        procedure, pass(self), public :: apply_matvec => apply_matvec_rsp
        !! Wrapper for matvec including the counter increment
        procedure, pass(self), public :: apply_rmatvec => apply_rmatvec_rsp
        !! Wrapper for rmatvec including the counter increment
    end type

    abstract interface
        subroutine abstract_matvec_rsp(self, vec_in, vec_out)
            !! Interface for the matrix-vector product.
            use lightkrylov_AbstractVectors
            import abstract_linop_rsp
            class(abstract_linop_rsp) , intent(inout)  :: self
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

    !--------------------------------------------------------------------------------------------
    !-----     Definition of an abstract real(sp) exponential propagator with kind=sp     -----
    !--------------------------------------------------------------------------------------------
    type, abstract, extends(abstract_linop_rsp), public :: abstract_exptA_linop_rsp
        !! Utility type to define the exponential propagator \( \mathbf{\Phi}_\tau \) which is the linear map 
        !! corresponding to the matrix exponential of the (possibly time-dependent) system Jacobian 
        !! \( \mathbf{L}(t) \) over a time horizon \( \tau \) as:
        !! $$ \mathbf{\Phi}_\tau = \int_0^\tau \mathbf{L}(t) \: \text{d}t $$
        !! Note that explicit knowledge or definition of the Jacobian is NOT required. This utility function
        !! is intended for the use in conjuction with a time-stepper algorithm that computes the integral
        !! directly.
        !!
        !!  @warning
        !!  While it is not necessary to use this utility operator, it is strongly recommended for operators
        !!  that correspond to exponential propagators to extend from this abstract type to allow for more
        !!  rigorous type checks in the application.
        !!  @endwarning
        real(sp), public :: tau
        !! Time horizon for the temporal integration. This variable must be set when the operator is instantiated.
    end type
    !------------------------------------------------------------------------------
    !-----     Definition of an abstract real(dp) operator with kind=dp     -----
    !------------------------------------------------------------------------------
    type, abstract, extends(abstract_linop), public :: abstract_linop_rdp
        !! Base type to extend in order to define a real(dp)-valued linear operator.
    contains
        private
        ! User defined procedures
        procedure(abstract_matvec_rdp), pass(self), deferred, public :: matvec
        !! Procedure to compute the matrix-vector product \( \mathbf{y} = \mathbf{Ax} \).
        procedure(abstract_matvec_rdp), pass(self), deferred, public :: rmatvec
        !! Procedure to compute the reversed matrix-vector product \( \mathbf{y} = \mathbf{A}^H \mathbf{x} \).
        ! Wrappers including counter increment
        procedure, pass(self), public :: apply_matvec => apply_matvec_rdp
        !! Wrapper for matvec including the counter increment
        procedure, pass(self), public :: apply_rmatvec => apply_rmatvec_rdp
        !! Wrapper for rmatvec including the counter increment
    end type

    abstract interface
        subroutine abstract_matvec_rdp(self, vec_in, vec_out)
            !! Interface for the matrix-vector product.
            use lightkrylov_AbstractVectors
            import abstract_linop_rdp
            class(abstract_linop_rdp) , intent(inout)  :: self
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

    !--------------------------------------------------------------------------------------------
    !-----     Definition of an abstract real(dp) exponential propagator with kind=dp     -----
    !--------------------------------------------------------------------------------------------
    type, abstract, extends(abstract_linop_rdp), public :: abstract_exptA_linop_rdp
        !! Utility type to define the exponential propagator \( \mathbf{\Phi}_\tau \) which is the linear map 
        !! corresponding to the matrix exponential of the (possibly time-dependent) system Jacobian 
        !! \( \mathbf{L}(t) \) over a time horizon \( \tau \) as:
        !! $$ \mathbf{\Phi}_\tau = \int_0^\tau \mathbf{L}(t) \: \text{d}t $$
        !! Note that explicit knowledge or definition of the Jacobian is NOT required. This utility function
        !! is intended for the use in conjuction with a time-stepper algorithm that computes the integral
        !! directly.
        !!
        !!  @warning
        !!  While it is not necessary to use this utility operator, it is strongly recommended for operators
        !!  that correspond to exponential propagators to extend from this abstract type to allow for more
        !!  rigorous type checks in the application.
        !!  @endwarning
        real(dp), public :: tau
        !! Time horizon for the temporal integration. This variable must be set when the operator is instantiated.
    end type
    !------------------------------------------------------------------------------
    !-----     Definition of an abstract complex(sp) operator with kind=sp     -----
    !------------------------------------------------------------------------------
    type, abstract, extends(abstract_linop), public :: abstract_linop_csp
        !! Base type to extend in order to define a complex(sp)-valued linear operator.
    contains
        private
        ! User defined procedures
        procedure(abstract_matvec_csp), pass(self), deferred, public :: matvec
        !! Procedure to compute the matrix-vector product \( \mathbf{y} = \mathbf{Ax} \).
        procedure(abstract_matvec_csp), pass(self), deferred, public :: rmatvec
        !! Procedure to compute the reversed matrix-vector product \( \mathbf{y} = \mathbf{A}^H \mathbf{x} \).
        ! Wrappers including counter increment
        procedure, pass(self), public :: apply_matvec => apply_matvec_csp
        !! Wrapper for matvec including the counter increment
        procedure, pass(self), public :: apply_rmatvec => apply_rmatvec_csp
        !! Wrapper for rmatvec including the counter increment
    end type

    abstract interface
        subroutine abstract_matvec_csp(self, vec_in, vec_out)
            !! Interface for the matrix-vector product.
            use lightkrylov_AbstractVectors
            import abstract_linop_csp
            class(abstract_linop_csp) , intent(inout)  :: self
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

    !--------------------------------------------------------------------------------------------
    !-----     Definition of an abstract complex(sp) exponential propagator with kind=sp     -----
    !--------------------------------------------------------------------------------------------
    type, abstract, extends(abstract_linop_csp), public :: abstract_exptA_linop_csp
        !! Utility type to define the exponential propagator \( \mathbf{\Phi}_\tau \) which is the linear map 
        !! corresponding to the matrix exponential of the (possibly time-dependent) system Jacobian 
        !! \( \mathbf{L}(t) \) over a time horizon \( \tau \) as:
        !! $$ \mathbf{\Phi}_\tau = \int_0^\tau \mathbf{L}(t) \: \text{d}t $$
        !! Note that explicit knowledge or definition of the Jacobian is NOT required. This utility function
        !! is intended for the use in conjuction with a time-stepper algorithm that computes the integral
        !! directly.
        !!
        !!  @warning
        !!  While it is not necessary to use this utility operator, it is strongly recommended for operators
        !!  that correspond to exponential propagators to extend from this abstract type to allow for more
        !!  rigorous type checks in the application.
        !!  @endwarning
        real(sp), public :: tau
        !! Time horizon for the temporal integration. This variable must be set when the operator is instantiated.
    end type
    !------------------------------------------------------------------------------
    !-----     Definition of an abstract complex(dp) operator with kind=dp     -----
    !------------------------------------------------------------------------------
    type, abstract, extends(abstract_linop), public :: abstract_linop_cdp
        !! Base type to extend in order to define a complex(dp)-valued linear operator.
    contains
        private
        ! User defined procedures
        procedure(abstract_matvec_cdp), pass(self), deferred, public :: matvec
        !! Procedure to compute the matrix-vector product \( \mathbf{y} = \mathbf{Ax} \).
        procedure(abstract_matvec_cdp), pass(self), deferred, public :: rmatvec
        !! Procedure to compute the reversed matrix-vector product \( \mathbf{y} = \mathbf{A}^H \mathbf{x} \).
        ! Wrappers including counter increment
        procedure, pass(self), public :: apply_matvec => apply_matvec_cdp
        !! Wrapper for matvec including the counter increment
        procedure, pass(self), public :: apply_rmatvec => apply_rmatvec_cdp
        !! Wrapper for rmatvec including the counter increment
    end type

    abstract interface
        subroutine abstract_matvec_cdp(self, vec_in, vec_out)
            !! Interface for the matrix-vector product.
            use lightkrylov_AbstractVectors
            import abstract_linop_cdp
            class(abstract_linop_cdp) , intent(inout)  :: self
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

    !--------------------------------------------------------------------------------------------
    !-----     Definition of an abstract complex(dp) exponential propagator with kind=dp     -----
    !--------------------------------------------------------------------------------------------
    type, abstract, extends(abstract_linop_cdp), public :: abstract_exptA_linop_cdp
        !! Utility type to define the exponential propagator \( \mathbf{\Phi}_\tau \) which is the linear map 
        !! corresponding to the matrix exponential of the (possibly time-dependent) system Jacobian 
        !! \( \mathbf{L}(t) \) over a time horizon \( \tau \) as:
        !! $$ \mathbf{\Phi}_\tau = \int_0^\tau \mathbf{L}(t) \: \text{d}t $$
        !! Note that explicit knowledge or definition of the Jacobian is NOT required. This utility function
        !! is intended for the use in conjuction with a time-stepper algorithm that computes the integral
        !! directly.
        !!
        !!  @warning
        !!  While it is not necessary to use this utility operator, it is strongly recommended for operators
        !!  that correspond to exponential propagators to extend from this abstract type to allow for more
        !!  rigorous type checks in the application.
        !!  @endwarning
        real(dp), public :: tau
        !! Time horizon for the temporal integration. This variable must be set when the operator is instantiated.
    end type

    interface adjoint
        module procedure initialize_adjoint_rsp
        module procedure initialize_adjoint_rdp
        module procedure initialize_adjoint_csp
        module procedure initialize_adjoint_cdp
    end interface
    public :: adjoint

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
        !! Defines a scaled linear operator \( \mathbf{B} = \sigma \mathbf{A} \) with \( \mathbf{A} \) a real-valued operator and
        !! \( \sigma \in \mathbb{R} \). The definitions of `matvec` and `rmatvec` are directly inherited from those used to define
        !! `A` and do not have to be defined by the user.
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
        !! Defines a scaled linear operator \( \mathbf{B} = \sigma \mathbf{A} \) with \( \mathbf{A} \) a real-valued operator and
        !! \( \sigma \in \mathbb{R} \). The definitions of `matvec` and `rmatvec` are directly inherited from those used to define
        !! `A` and do not have to be defined by the user.
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
        !! Defines a scaled linear operator \( \mathbf{B} = \sigma \mathbf{A} \) with \( \mathbf{A} \) a complex(sp)-valued operator
        !! and \( \sigma \in \mathbb{R}\ ) or \( \mathbb{C} \).
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
        !! Defines a scaled linear operator \( \mathbf{B} = \sigma \mathbf{A} \) with \( \mathbf{A} \) a complex(dp)-valued operator
        !! and \( \sigma \in \mathbb{R}\ ) or \( \mathbb{C} \).
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
        !! The definitions of `matvec` and `rmatvec` are directly inherited from those used to define `A` and `B`.
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
        !! The definitions of `matvec` and `rmatvec` are directly inherited from those used to define `A` and `B`.
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
        !! The definitions of `matvec` and `rmatvec` are directly inherited from those used to define `A` and `B`.
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
        !! The definitions of `matvec` and `rmatvec` are directly inherited from those used to define `A` and `B`.
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

    !------------------------------------------------
    !-----     Convenience dense linop type     -----
    !------------------------------------------------

    type, extends(abstract_linop_rsp), public :: dense_linop_rsp
        real(sp), allocatable :: data(:, :)
    contains
        procedure, pass(self), public :: matvec => dense_matvec_rsp
        procedure, pass(self), public :: rmatvec => dense_rmatvec_rsp
    end type
    type, extends(abstract_linop_rdp), public :: dense_linop_rdp
        real(dp), allocatable :: data(:, :)
    contains
        procedure, pass(self), public :: matvec => dense_matvec_rdp
        procedure, pass(self), public :: rmatvec => dense_rmatvec_rdp
    end type
    type, extends(abstract_linop_csp), public :: dense_linop_csp
        complex(sp), allocatable :: data(:, :)
    contains
        procedure, pass(self), public :: matvec => dense_matvec_csp
        procedure, pass(self), public :: rmatvec => dense_rmatvec_csp
    end type
    type, extends(abstract_linop_cdp), public :: dense_linop_cdp
        complex(dp), allocatable :: data(:, :)
    contains
        procedure, pass(self), public :: matvec => dense_matvec_cdp
        procedure, pass(self), public :: rmatvec => dense_rmatvec_cdp
    end type

    interface dense_linop
        module procedure initialize_dense_linop_from_array_rsp
        module procedure initialize_dense_linop_from_array_rdp
        module procedure initialize_dense_linop_from_array_csp
        module procedure initialize_dense_linop_from_array_cdp
    end interface
    public :: dense_linop
   
contains

    !--------------------------------------------------------------
    !-----     Getter/Setter routines for abstract_linops     -----
    !--------------------------------------------------------------

    pure integer function get_counter(self, trans) result(count)
      !! Getter function for the number of matvec calls
      class(abstract_linop), intent(in) :: self
      logical, intent(in) :: trans
      !! matvec or rmatvec?
      if (trans) then
         count = self%rmatvec_counter
      else
         count = self%matvec_counter
      end if
    end function get_counter

    subroutine reset_counter(self, trans, procedure, counter, reset_timer, soft_reset, clean_timer)
      !! Setter routine to reset the matvec counter and reset timers
      class(abstract_linop), intent(inout) :: self
      logical, intent(in) :: trans
      !! matvec or rmatvec?
      character(len=*), intent(in) :: procedure
      !! name of the caller routine
      integer, optional, intent(in) :: counter
      !! optional flag to reset to an integer other than zero.
      logical, optional, intent(in) :: reset_timer
      !! optional flag to reset also the timers
      logical, optional, intent(in) :: soft_reset
      !! optional flag to choose whether to save previous timing data (default: .true.)
      logical, optional, intent(in) :: clean_timer
      !! optional flag to choose whether to fully reset the timer (default: .false.)
      ! internals
      integer :: counter_, count_old
      logical :: reset_timer_ 
      character(len=128) :: msg
      counter_ = optval(counter, 0)
      count_old = self%get_counter(trans)
      reset_timer_ = optval(reset_timer, .false.)
      if ( count_old /= 0 .or. counter_ /= 0) then
        if (trans) then
            write(msg,'(A,I0,A,I0,A)') 'Total number of rmatvecs: ', count_old, '. Resetting counter to ', counter_, '.'
            call log_message(msg, this_module, 'reset_counter('//trim(procedure)//')')
            self%rmatvec_counter = counter_
        else
            write(msg,'(A,I0,A,I0,A)') 'Total number of matvecs: ', count_old, '. Resetting counter to ', counter_, '.'
            call log_message(msg, this_module, 'reset_counter('//trim(procedure)//')')
            self%matvec_counter = counter_
        end if
      end if
      if (reset_timer_) call self%reset_timer(trans, soft_reset, clean_timer)
      return
    end subroutine reset_counter

    subroutine print_timer_info(self, trans)
      !! Getter routine to print the current timing information for matvec/rmatvec
      !! Note: Wrapper of the corresponding routine from lightkrylov_timer
      class(abstract_linop), intent(inout) :: self
      logical, optional, intent(in) :: trans
      !! matvec or rmatvec?
      ! internal
      logical :: transpose
      transpose = optval(trans, .false.)
      if (transpose) then
        call self%rmatvec_timer%print_info()
      else
        call self%matvec_timer%print_info()
      end if
    end subroutine print_timer_info

    subroutine reset_linop_timer(self, trans, soft, clean)
      !! Setter routine to reset the matvec/rmatvec timers
      !! Note: Wrapper of the corresponding routine from lightkrylov_timer
      class(abstract_linop), intent(inout) :: self
      logical, optional, intent(in) :: trans
      !! matvec or rmatvec?
      logical, optional, intent(in) :: soft
      logical, optional, intent(in) :: clean
      ! internal
      logical :: transpose
      transpose = optval(trans, .false.)
      if (present(trans)) then
        if (transpose) then
            call self%rmatvec_timer%reset(soft, clean, verbose=.true.)
        else
            call self%matvec_timer%reset(soft, clean, verbose=.true.)
        end if
      else
        call self%rmatvec_timer%reset(soft, clean, verbose=.true.)
        call self%matvec_timer%reset(soft, clean, verbose=.true.)
      end if
    end subroutine reset_linop_timer

    subroutine finalize_linop_timer(self)
      !! Finalize the matvec/rmatvec timers
      !! Note: Wrapper of the corresponding routine from lightkrylov_timer
      class(abstract_linop), intent(inout) :: self
      call self%matvec_timer%finalize()
      call self%rmatvec_timer%finalize()
    end subroutine finalize_linop_timer

    !---------------------------------------------------------------------
    !-----     Wrappers for matvec/rmatvec to increment counters     -----
    !---------------------------------------------------------------------

    subroutine apply_matvec_rsp(self, vec_in, vec_out)
        class(abstract_linop_rsp), intent(inout) :: self
        class(abstract_vector_rsp), intent(in) :: vec_in
        class(abstract_vector_rsp), intent(out) :: vec_out
        ! internal
        character(len=128) :: msg
        self%matvec_counter = self%matvec_counter + 1
        write(msg,'(I0,1X,A)') self%matvec_counter, 'start'
        call log_debug(msg, this_module, 'matvec')
        call self%matvec_timer%start()
        call self%matvec(vec_in, vec_out)
        call self%matvec_timer%stop()
        write(msg,'(I0,1X,A)') self%matvec_counter, 'end'
        call log_debug(msg, this_module, 'matvec')
        return
    end subroutine apply_matvec_rsp

    subroutine apply_rmatvec_rsp(self, vec_in, vec_out)
        class(abstract_linop_rsp), intent(inout) :: self
        class(abstract_vector_rsp), intent(in) :: vec_in
        class(abstract_vector_rsp), intent(out) :: vec_out
        ! internal 
        character(len=128) :: msg
        self%rmatvec_counter = self%rmatvec_counter + 1
        write(msg,'(I0,1X,A)') self%rmatvec_counter, 'start'
        call log_debug(msg, this_module, 'rmatvec')
        call self%rmatvec_timer%start()
        call self%rmatvec(vec_in, vec_out)
        call self%rmatvec_timer%stop()
        write(msg,'(I0,1X,A)') self%rmatvec_counter, 'end'
        call log_debug(msg, this_module, 'rmatvec')
        return
    end subroutine apply_rmatvec_rsp
    subroutine apply_matvec_rdp(self, vec_in, vec_out)
        class(abstract_linop_rdp), intent(inout) :: self
        class(abstract_vector_rdp), intent(in) :: vec_in
        class(abstract_vector_rdp), intent(out) :: vec_out
        ! internal
        character(len=128) :: msg
        self%matvec_counter = self%matvec_counter + 1
        write(msg,'(I0,1X,A)') self%matvec_counter, 'start'
        call log_debug(msg, this_module, 'matvec')
        call self%matvec_timer%start()
        call self%matvec(vec_in, vec_out)
        call self%matvec_timer%stop()
        write(msg,'(I0,1X,A)') self%matvec_counter, 'end'
        call log_debug(msg, this_module, 'matvec')
        return
    end subroutine apply_matvec_rdp

    subroutine apply_rmatvec_rdp(self, vec_in, vec_out)
        class(abstract_linop_rdp), intent(inout) :: self
        class(abstract_vector_rdp), intent(in) :: vec_in
        class(abstract_vector_rdp), intent(out) :: vec_out
        ! internal 
        character(len=128) :: msg
        self%rmatvec_counter = self%rmatvec_counter + 1
        write(msg,'(I0,1X,A)') self%rmatvec_counter, 'start'
        call log_debug(msg, this_module, 'rmatvec')
        call self%rmatvec_timer%start()
        call self%rmatvec(vec_in, vec_out)
        call self%rmatvec_timer%stop()
        write(msg,'(I0,1X,A)') self%rmatvec_counter, 'end'
        call log_debug(msg, this_module, 'rmatvec')
        return
    end subroutine apply_rmatvec_rdp
    subroutine apply_matvec_csp(self, vec_in, vec_out)
        class(abstract_linop_csp), intent(inout) :: self
        class(abstract_vector_csp), intent(in) :: vec_in
        class(abstract_vector_csp), intent(out) :: vec_out
        ! internal
        character(len=128) :: msg
        self%matvec_counter = self%matvec_counter + 1
        write(msg,'(I0,1X,A)') self%matvec_counter, 'start'
        call log_debug(msg, this_module, 'matvec')
        call self%matvec_timer%start()
        call self%matvec(vec_in, vec_out)
        call self%matvec_timer%stop()
        write(msg,'(I0,1X,A)') self%matvec_counter, 'end'
        call log_debug(msg, this_module, 'matvec')
        return
    end subroutine apply_matvec_csp

    subroutine apply_rmatvec_csp(self, vec_in, vec_out)
        class(abstract_linop_csp), intent(inout) :: self
        class(abstract_vector_csp), intent(in) :: vec_in
        class(abstract_vector_csp), intent(out) :: vec_out
        ! internal 
        character(len=128) :: msg
        self%rmatvec_counter = self%rmatvec_counter + 1
        write(msg,'(I0,1X,A)') self%rmatvec_counter, 'start'
        call log_debug(msg, this_module, 'rmatvec')
        call self%rmatvec_timer%start()
        call self%rmatvec(vec_in, vec_out)
        call self%rmatvec_timer%stop()
        write(msg,'(I0,1X,A)') self%rmatvec_counter, 'end'
        call log_debug(msg, this_module, 'rmatvec')
        return
    end subroutine apply_rmatvec_csp
    subroutine apply_matvec_cdp(self, vec_in, vec_out)
        class(abstract_linop_cdp), intent(inout) :: self
        class(abstract_vector_cdp), intent(in) :: vec_in
        class(abstract_vector_cdp), intent(out) :: vec_out
        ! internal
        character(len=128) :: msg
        self%matvec_counter = self%matvec_counter + 1
        write(msg,'(I0,1X,A)') self%matvec_counter, 'start'
        call log_debug(msg, this_module, 'matvec')
        call self%matvec_timer%start()
        call self%matvec(vec_in, vec_out)
        call self%matvec_timer%stop()
        write(msg,'(I0,1X,A)') self%matvec_counter, 'end'
        call log_debug(msg, this_module, 'matvec')
        return
    end subroutine apply_matvec_cdp

    subroutine apply_rmatvec_cdp(self, vec_in, vec_out)
        class(abstract_linop_cdp), intent(inout) :: self
        class(abstract_vector_cdp), intent(in) :: vec_in
        class(abstract_vector_cdp), intent(out) :: vec_out
        ! internal 
        character(len=128) :: msg
        self%rmatvec_counter = self%rmatvec_counter + 1
        write(msg,'(I0,1X,A)') self%rmatvec_counter, 'start'
        call log_debug(msg, this_module, 'rmatvec')
        call self%rmatvec_timer%start()
        call self%rmatvec(vec_in, vec_out)
        call self%rmatvec_timer%stop()
        write(msg,'(I0,1X,A)') self%rmatvec_counter, 'end'
        call log_debug(msg, this_module, 'rmatvec')
        return
    end subroutine apply_rmatvec_cdp

    !------------------------------------------------------------------------------
    !-----     Concrete matvec/rmatvec implementations for special linops     -----
    !------------------------------------------------------------------------------

    subroutine id_matvec_rsp(self, vec_in, vec_out)
        class(Id_rsp), intent(inout) :: self
        class(abstract_vector_rsp), intent(in) :: vec_in
        class(abstract_vector_rsp), intent(out) :: vec_out
        call copy(vec_out, vec_in)
        return
    end subroutine id_matvec_rsp
    subroutine id_matvec_rdp(self, vec_in, vec_out)
        class(Id_rdp), intent(inout) :: self
        class(abstract_vector_rdp), intent(in) :: vec_in
        class(abstract_vector_rdp), intent(out) :: vec_out
        call copy(vec_out, vec_in)
        return
    end subroutine id_matvec_rdp
    subroutine id_matvec_csp(self, vec_in, vec_out)
        class(Id_csp), intent(inout) :: self
        class(abstract_vector_csp), intent(in) :: vec_in
        class(abstract_vector_csp), intent(out) :: vec_out
        call copy(vec_out, vec_in)
        return
    end subroutine id_matvec_csp
    subroutine id_matvec_cdp(self, vec_in, vec_out)
        class(Id_cdp), intent(inout) :: self
        class(abstract_vector_cdp), intent(in) :: vec_in
        class(abstract_vector_cdp), intent(out) :: vec_out
        call copy(vec_out, vec_in)
        return
    end subroutine id_matvec_cdp

    subroutine scaled_matvec_rsp(self, vec_in, vec_out)
        class(scaled_linop_rsp), intent(inout) :: self
        class(abstract_vector_rsp), intent(in) :: vec_in
        class(abstract_vector_rsp), intent(out) :: vec_out
        call self%A%apply_matvec(vec_in, vec_out) ; call vec_out%scal(self%sigma)
        return
    end subroutine scaled_matvec_rsp

    subroutine scaled_rmatvec_rsp(self, vec_in, vec_out)
        class(scaled_linop_rsp), intent(inout) :: self
        class(abstract_vector_rsp), intent(in) :: vec_in
        class(abstract_vector_rsp), intent(out) :: vec_out
        call self%A%apply_rmatvec(vec_in, vec_out) ; call vec_out%scal(self%sigma)
        return
    end subroutine scaled_rmatvec_rsp
    subroutine scaled_matvec_rdp(self, vec_in, vec_out)
        class(scaled_linop_rdp), intent(inout) :: self
        class(abstract_vector_rdp), intent(in) :: vec_in
        class(abstract_vector_rdp), intent(out) :: vec_out
        call self%A%apply_matvec(vec_in, vec_out) ; call vec_out%scal(self%sigma)
        return
    end subroutine scaled_matvec_rdp

    subroutine scaled_rmatvec_rdp(self, vec_in, vec_out)
        class(scaled_linop_rdp), intent(inout) :: self
        class(abstract_vector_rdp), intent(in) :: vec_in
        class(abstract_vector_rdp), intent(out) :: vec_out
        call self%A%apply_rmatvec(vec_in, vec_out) ; call vec_out%scal(self%sigma)
        return
    end subroutine scaled_rmatvec_rdp
    subroutine scaled_matvec_csp(self, vec_in, vec_out)
        class(scaled_linop_csp), intent(inout) :: self
        class(abstract_vector_csp), intent(in) :: vec_in
        class(abstract_vector_csp), intent(out) :: vec_out
        call self%A%apply_matvec(vec_in, vec_out) ; call vec_out%scal(self%sigma)
        return
    end subroutine scaled_matvec_csp

    subroutine scaled_rmatvec_csp(self, vec_in, vec_out)
        class(scaled_linop_csp), intent(inout) :: self
        class(abstract_vector_csp), intent(in) :: vec_in
        class(abstract_vector_csp), intent(out) :: vec_out
        call self%A%apply_rmatvec(vec_in, vec_out) ; call vec_out%scal(self%sigma)
        return
    end subroutine scaled_rmatvec_csp
    subroutine scaled_matvec_cdp(self, vec_in, vec_out)
        class(scaled_linop_cdp), intent(inout) :: self
        class(abstract_vector_cdp), intent(in) :: vec_in
        class(abstract_vector_cdp), intent(out) :: vec_out
        call self%A%apply_matvec(vec_in, vec_out) ; call vec_out%scal(self%sigma)
        return
    end subroutine scaled_matvec_cdp

    subroutine scaled_rmatvec_cdp(self, vec_in, vec_out)
        class(scaled_linop_cdp), intent(inout) :: self
        class(abstract_vector_cdp), intent(in) :: vec_in
        class(abstract_vector_cdp), intent(out) :: vec_out
        call self%A%apply_rmatvec(vec_in, vec_out) ; call vec_out%scal(self%sigma)
        return
    end subroutine scaled_rmatvec_cdp

    subroutine axpby_matvec_rsp(self, vec_in, vec_out)
        class(axpby_linop_rsp), intent(inout) :: self
        class(abstract_vector_rsp), intent(in) :: vec_in
        class(abstract_vector_rsp), intent(out) :: vec_out

        ! Working array.
        class(abstract_vector_rsp), allocatable :: wrk

        ! Allocate working array.
        allocate(wrk, mold=vec_in) ; call wrk%zero()

        ! w = A @ x
        if (self%transA) then
            call self%A%apply_rmatvec(vec_in, wrk)
        else
            call self%A%apply_matvec(vec_in, wrk)
        endif

        ! y = B @ x
        if (self%transB) then
            call self%B%apply_rmatvec(vec_in, vec_out)
        else
            call self%B%apply_matvec(vec_in, vec_out)
        endif

        ! y = alpha*w + beta*y
        call vec_out%axpby(self%alpha, wrk, self%beta)

        return
    end subroutine axpby_matvec_rsp

    subroutine axpby_rmatvec_rsp(self, vec_in, vec_out)
        class(axpby_linop_rsp), intent(inout) :: self
        class(abstract_vector_rsp), intent(in) :: vec_in
        class(abstract_vector_rsp), intent(out) :: vec_out

        ! Working array.
        class(abstract_vector_rsp), allocatable :: wrk

        ! Allocate working array.
        allocate(wrk, mold=vec_in) ; call wrk%zero()

        ! w = A @ x
        if (self%transA) then
            call self%A%apply_matvec(vec_in, wrk)
        else
            call self%A%apply_rmatvec(vec_in, wrk)
        endif

        ! y = B @ x
        if (self%transB) then
            call self%B%apply_matvec(vec_in, vec_out)
        else
            call self%B%apply_rmatvec(vec_in, vec_out)
        endif

        ! y = alpha*w + beta*y
        call vec_out%axpby(self%alpha, wrk, self%beta)

        return
    end subroutine axpby_rmatvec_rsp

    subroutine axpby_matvec_rdp(self, vec_in, vec_out)
        class(axpby_linop_rdp), intent(inout) :: self
        class(abstract_vector_rdp), intent(in) :: vec_in
        class(abstract_vector_rdp), intent(out) :: vec_out

        ! Working array.
        class(abstract_vector_rdp), allocatable :: wrk

        ! Allocate working array.
        allocate(wrk, mold=vec_in) ; call wrk%zero()

        ! w = A @ x
        if (self%transA) then
            call self%A%apply_rmatvec(vec_in, wrk)
        else
            call self%A%apply_matvec(vec_in, wrk)
        endif

        ! y = B @ x
        if (self%transB) then
            call self%B%apply_rmatvec(vec_in, vec_out)
        else
            call self%B%apply_matvec(vec_in, vec_out)
        endif

        ! y = alpha*w + beta*y
        call vec_out%axpby(self%alpha, wrk, self%beta)

        return
    end subroutine axpby_matvec_rdp

    subroutine axpby_rmatvec_rdp(self, vec_in, vec_out)
        class(axpby_linop_rdp), intent(inout) :: self
        class(abstract_vector_rdp), intent(in) :: vec_in
        class(abstract_vector_rdp), intent(out) :: vec_out

        ! Working array.
        class(abstract_vector_rdp), allocatable :: wrk

        ! Allocate working array.
        allocate(wrk, mold=vec_in) ; call wrk%zero()

        ! w = A @ x
        if (self%transA) then
            call self%A%apply_matvec(vec_in, wrk)
        else
            call self%A%apply_rmatvec(vec_in, wrk)
        endif

        ! y = B @ x
        if (self%transB) then
            call self%B%apply_matvec(vec_in, vec_out)
        else
            call self%B%apply_rmatvec(vec_in, vec_out)
        endif

        ! y = alpha*w + beta*y
        call vec_out%axpby(self%alpha, wrk, self%beta)

        return
    end subroutine axpby_rmatvec_rdp

    subroutine axpby_matvec_csp(self, vec_in, vec_out)
        class(axpby_linop_csp), intent(inout) :: self
        class(abstract_vector_csp), intent(in) :: vec_in
        class(abstract_vector_csp), intent(out) :: vec_out

        ! Working array.
        class(abstract_vector_csp), allocatable :: wrk

        ! Allocate working array.
        allocate(wrk, mold=vec_in) ; call wrk%zero()

        ! w = A @ x
        if (self%transA) then
            call self%A%apply_rmatvec(vec_in, wrk)
        else
            call self%A%apply_matvec(vec_in, wrk)
        endif

        ! y = B @ x
        if (self%transB) then
            call self%B%apply_rmatvec(vec_in, vec_out)
        else
            call self%B%apply_matvec(vec_in, vec_out)
        endif

        ! y = alpha*w + beta*y
        call vec_out%axpby(self%alpha, wrk, self%beta)

        return
    end subroutine axpby_matvec_csp

    subroutine axpby_rmatvec_csp(self, vec_in, vec_out)
        class(axpby_linop_csp), intent(inout) :: self
        class(abstract_vector_csp), intent(in) :: vec_in
        class(abstract_vector_csp), intent(out) :: vec_out

        ! Working array.
        class(abstract_vector_csp), allocatable :: wrk

        ! Allocate working array.
        allocate(wrk, mold=vec_in) ; call wrk%zero()

        ! w = A @ x
        if (self%transA) then
            call self%A%apply_matvec(vec_in, wrk)
        else
            call self%A%apply_rmatvec(vec_in, wrk)
        endif

        ! y = B @ x
        if (self%transB) then
            call self%B%apply_matvec(vec_in, vec_out)
        else
            call self%B%apply_rmatvec(vec_in, vec_out)
        endif

        ! y = alpha*w + beta*y
        call vec_out%axpby(self%alpha, wrk, self%beta)

        return
    end subroutine axpby_rmatvec_csp

    subroutine axpby_matvec_cdp(self, vec_in, vec_out)
        class(axpby_linop_cdp), intent(inout) :: self
        class(abstract_vector_cdp), intent(in) :: vec_in
        class(abstract_vector_cdp), intent(out) :: vec_out

        ! Working array.
        class(abstract_vector_cdp), allocatable :: wrk

        ! Allocate working array.
        allocate(wrk, mold=vec_in) ; call wrk%zero()

        ! w = A @ x
        if (self%transA) then
            call self%A%apply_rmatvec(vec_in, wrk)
        else
            call self%A%apply_matvec(vec_in, wrk)
        endif

        ! y = B @ x
        if (self%transB) then
            call self%B%apply_rmatvec(vec_in, vec_out)
        else
            call self%B%apply_matvec(vec_in, vec_out)
        endif

        ! y = alpha*w + beta*y
        call vec_out%axpby(self%alpha, wrk, self%beta)

        return
    end subroutine axpby_matvec_cdp

    subroutine axpby_rmatvec_cdp(self, vec_in, vec_out)
        class(axpby_linop_cdp), intent(inout) :: self
        class(abstract_vector_cdp), intent(in) :: vec_in
        class(abstract_vector_cdp), intent(out) :: vec_out

        ! Working array.
        class(abstract_vector_cdp), allocatable :: wrk

        ! Allocate working array.
        allocate(wrk, mold=vec_in) ; call wrk%zero()

        ! w = A @ x
        if (self%transA) then
            call self%A%apply_matvec(vec_in, wrk)
        else
            call self%A%apply_rmatvec(vec_in, wrk)
        endif

        ! y = B @ x
        if (self%transB) then
            call self%B%apply_matvec(vec_in, vec_out)
        else
            call self%B%apply_rmatvec(vec_in, vec_out)
        endif

        ! y = alpha*w + beta*y
        call vec_out%axpby(self%alpha, wrk, self%beta)

        return
    end subroutine axpby_rmatvec_cdp


    !-------------------------------------------------
    !-----     ADJOINT TYPE-BOUND PROCEDURES     -----
    !-------------------------------------------------

    function initialize_adjoint_rsp(A) result(B)
        class(abstract_linop_rsp), intent(in) :: A
        class(adjoint_linop_rsp), allocatable :: B
        allocate(B) ; B%A = A
        return
    end function

    subroutine adjoint_matvec_rsp(self, vec_in, vec_out)
        class(adjoint_linop_rsp), intent(inout) :: self
        class(abstract_vector_rsp), intent(in) :: vec_in
        class(abstract_vector_rsp), intent(out) :: vec_out

        call self%A%apply_rmatvec(vec_in, vec_out)

        return
    end subroutine adjoint_matvec_rsp

    subroutine adjoint_rmatvec_rsp(self, vec_in, vec_out)
        class(adjoint_linop_rsp), intent(inout) :: self
        class(abstract_vector_rsp), intent(in) :: vec_in
        class(abstract_vector_rsp), intent(out) :: vec_out

        call self%A%apply_matvec(vec_in, vec_out)

        return
    end subroutine adjoint_rmatvec_rsp

    function initialize_adjoint_rdp(A) result(B)
        class(abstract_linop_rdp), intent(in) :: A
        class(adjoint_linop_rdp), allocatable :: B
        allocate(B) ; B%A = A
        return
    end function

    subroutine adjoint_matvec_rdp(self, vec_in, vec_out)
        class(adjoint_linop_rdp), intent(inout) :: self
        class(abstract_vector_rdp), intent(in) :: vec_in
        class(abstract_vector_rdp), intent(out) :: vec_out

        call self%A%apply_rmatvec(vec_in, vec_out)

        return
    end subroutine adjoint_matvec_rdp

    subroutine adjoint_rmatvec_rdp(self, vec_in, vec_out)
        class(adjoint_linop_rdp), intent(inout) :: self
        class(abstract_vector_rdp), intent(in) :: vec_in
        class(abstract_vector_rdp), intent(out) :: vec_out

        call self%A%apply_matvec(vec_in, vec_out)

        return
    end subroutine adjoint_rmatvec_rdp

    function initialize_adjoint_csp(A) result(B)
        class(abstract_linop_csp), intent(in) :: A
        class(adjoint_linop_csp), allocatable :: B
        allocate(B) ; B%A = A
        return
    end function

    subroutine adjoint_matvec_csp(self, vec_in, vec_out)
        class(adjoint_linop_csp), intent(inout) :: self
        class(abstract_vector_csp), intent(in) :: vec_in
        class(abstract_vector_csp), intent(out) :: vec_out

        call self%A%apply_rmatvec(vec_in, vec_out)

        return
    end subroutine adjoint_matvec_csp

    subroutine adjoint_rmatvec_csp(self, vec_in, vec_out)
        class(adjoint_linop_csp), intent(inout) :: self
        class(abstract_vector_csp), intent(in) :: vec_in
        class(abstract_vector_csp), intent(out) :: vec_out

        call self%A%apply_matvec(vec_in, vec_out)

        return
    end subroutine adjoint_rmatvec_csp

    function initialize_adjoint_cdp(A) result(B)
        class(abstract_linop_cdp), intent(in) :: A
        class(adjoint_linop_cdp), allocatable :: B
        allocate(B) ; B%A = A
        return
    end function

    subroutine adjoint_matvec_cdp(self, vec_in, vec_out)
        class(adjoint_linop_cdp), intent(inout) :: self
        class(abstract_vector_cdp), intent(in) :: vec_in
        class(abstract_vector_cdp), intent(out) :: vec_out

        call self%A%apply_rmatvec(vec_in, vec_out)

        return
    end subroutine adjoint_matvec_cdp

    subroutine adjoint_rmatvec_cdp(self, vec_in, vec_out)
        class(adjoint_linop_cdp), intent(inout) :: self
        class(abstract_vector_cdp), intent(in) :: vec_in
        class(abstract_vector_cdp), intent(out) :: vec_out

        call self%A%apply_matvec(vec_in, vec_out)

        return
    end subroutine adjoint_rmatvec_cdp


    !--------------------------------------------------------------------------
    !-----     Type-bound procedures for convenience dense linop type     -----
    !--------------------------------------------------------------------------

    subroutine dense_matvec_rsp(self, vec_in, vec_out)
        class(dense_linop_rsp), intent(inout) :: self
        class(abstract_vector_rsp), intent(in) :: vec_in
        class(abstract_vector_rsp), intent(out) :: vec_out
        select type(vec_in)
        type is(dense_vector_rsp)
            select type(vec_out)
            type is(dense_vector_rsp)
                block
                integer :: m, n
                real(sp) :: alpha, beta
                m = size(self%data, 1) ; n = size(self%data, 2)
                alpha = one_rsp ; beta = zero_rsp
                vec_out = vec_in
                call gemv("N", m, n, alpha, self%data, m, vec_in%data, 1, beta, vec_out%data, 1)
                end block
            class default
                call stop_error("The intent [OUT] argument 'vec_out' must be of type 'dense_vector'", this_module, 'matvec')
            end select
        class default
            call stop_error("The intent [IN] argument 'vec_in' must be of type 'dense_vector'", this_module, 'matvec')
        end select
        return
    end subroutine

    subroutine dense_rmatvec_rsp(self, vec_in, vec_out)
        class(dense_linop_rsp), intent(inout) :: self
        class(abstract_vector_rsp), intent(in) :: vec_in
        class(abstract_vector_rsp), intent(out) :: vec_out
        select type(vec_in)
        type is(dense_vector_rsp)
            select type(vec_out)
            type is(dense_vector_rsp)
                block
                integer :: m, n
                real(sp) :: alpha, beta
                m = size(self%data, 1) ; n = size(self%data, 2)
                alpha = one_rsp ; beta = zero_rsp
                vec_out = vec_in
                call gemv("T", m, n, alpha, self%data, m, vec_in%data, 1, beta, vec_out%data, 1)
                end block
            class default
                call stop_error("The intent [OUT] argument 'vec_out' must be of type 'dense_vector'", this_module, 'matvec')
            end select
        class default
            call stop_error("The intent [IN] argument 'vec_in' must be of type 'dense_vector'", this_module, 'matvec')
        end select
         return
    end subroutine

    subroutine dense_matvec_rdp(self, vec_in, vec_out)
        class(dense_linop_rdp), intent(inout) :: self
        class(abstract_vector_rdp), intent(in) :: vec_in
        class(abstract_vector_rdp), intent(out) :: vec_out
        select type(vec_in)
        type is(dense_vector_rdp)
            select type(vec_out)
            type is(dense_vector_rdp)
                block
                integer :: m, n
                real(dp) :: alpha, beta
                m = size(self%data, 1) ; n = size(self%data, 2)
                alpha = one_rdp ; beta = zero_rdp
                vec_out = vec_in
                call gemv("N", m, n, alpha, self%data, m, vec_in%data, 1, beta, vec_out%data, 1)
                end block
            class default
                call stop_error("The intent [OUT] argument 'vec_out' must be of type 'dense_vector'", this_module, 'matvec')
            end select
        class default
            call stop_error("The intent [IN] argument 'vec_in' must be of type 'dense_vector'", this_module, 'matvec')
        end select
        return
    end subroutine

    subroutine dense_rmatvec_rdp(self, vec_in, vec_out)
        class(dense_linop_rdp), intent(inout) :: self
        class(abstract_vector_rdp), intent(in) :: vec_in
        class(abstract_vector_rdp), intent(out) :: vec_out
        select type(vec_in)
        type is(dense_vector_rdp)
            select type(vec_out)
            type is(dense_vector_rdp)
                block
                integer :: m, n
                real(dp) :: alpha, beta
                m = size(self%data, 1) ; n = size(self%data, 2)
                alpha = one_rdp ; beta = zero_rdp
                vec_out = vec_in
                call gemv("T", m, n, alpha, self%data, m, vec_in%data, 1, beta, vec_out%data, 1)
                end block
            class default
                call stop_error("The intent [OUT] argument 'vec_out' must be of type 'dense_vector'", this_module, 'matvec')
            end select
        class default
            call stop_error("The intent [IN] argument 'vec_in' must be of type 'dense_vector'", this_module, 'matvec')
        end select
         return
    end subroutine

    subroutine dense_matvec_csp(self, vec_in, vec_out)
        class(dense_linop_csp), intent(inout) :: self
        class(abstract_vector_csp), intent(in) :: vec_in
        class(abstract_vector_csp), intent(out) :: vec_out
        select type(vec_in)
        type is(dense_vector_csp)
            select type(vec_out)
            type is(dense_vector_csp)
                block
                integer :: m, n
                complex(sp) :: alpha, beta
                m = size(self%data, 1) ; n = size(self%data, 2)
                alpha = one_csp ; beta = zero_csp
                vec_out = vec_in
                call gemv("N", m, n, alpha, self%data, m, vec_in%data, 1, beta, vec_out%data, 1)
                end block
            class default
                call stop_error("The intent [OUT] argument 'vec_out' must be of type 'dense_vector'", this_module, 'matvec')
            end select
        class default
            call stop_error("The intent [IN] argument 'vec_in' must be of type 'dense_vector'", this_module, 'matvec')
        end select
        return
    end subroutine

    subroutine dense_rmatvec_csp(self, vec_in, vec_out)
        class(dense_linop_csp), intent(inout) :: self
        class(abstract_vector_csp), intent(in) :: vec_in
        class(abstract_vector_csp), intent(out) :: vec_out
        select type(vec_in)
        type is(dense_vector_csp)
            select type(vec_out)
            type is(dense_vector_csp)
                block
                integer :: m, n
                complex(sp) :: alpha, beta
                m = size(self%data, 1) ; n = size(self%data, 2)
                alpha = one_csp ; beta = zero_csp
                vec_out = vec_in
                call gemv("C", m, n, alpha, self%data, m, vec_in%data, 1, beta, vec_out%data, 1)
                end block
            class default
                call stop_error("The intent [OUT] argument 'vec_out' must be of type 'dense_vector'", this_module, 'matvec')
            end select
        class default
            call stop_error("The intent [IN] argument 'vec_in' must be of type 'dense_vector'", this_module, 'matvec')
        end select
         return
    end subroutine

    subroutine dense_matvec_cdp(self, vec_in, vec_out)
        class(dense_linop_cdp), intent(inout) :: self
        class(abstract_vector_cdp), intent(in) :: vec_in
        class(abstract_vector_cdp), intent(out) :: vec_out
        select type(vec_in)
        type is(dense_vector_cdp)
            select type(vec_out)
            type is(dense_vector_cdp)
                block
                integer :: m, n
                complex(dp) :: alpha, beta
                m = size(self%data, 1) ; n = size(self%data, 2)
                alpha = one_cdp ; beta = zero_cdp
                vec_out = vec_in
                call gemv("N", m, n, alpha, self%data, m, vec_in%data, 1, beta, vec_out%data, 1)
                end block
            class default
                call stop_error("The intent [OUT] argument 'vec_out' must be of type 'dense_vector'", this_module, 'matvec')
            end select
        class default
            call stop_error("The intent [IN] argument 'vec_in' must be of type 'dense_vector'", this_module, 'matvec')
        end select
        return
    end subroutine

    subroutine dense_rmatvec_cdp(self, vec_in, vec_out)
        class(dense_linop_cdp), intent(inout) :: self
        class(abstract_vector_cdp), intent(in) :: vec_in
        class(abstract_vector_cdp), intent(out) :: vec_out
        select type(vec_in)
        type is(dense_vector_cdp)
            select type(vec_out)
            type is(dense_vector_cdp)
                block
                integer :: m, n
                complex(dp) :: alpha, beta
                m = size(self%data, 1) ; n = size(self%data, 2)
                alpha = one_cdp ; beta = zero_cdp
                vec_out = vec_in
                call gemv("C", m, n, alpha, self%data, m, vec_in%data, 1, beta, vec_out%data, 1)
                end block
            class default
                call stop_error("The intent [OUT] argument 'vec_out' must be of type 'dense_vector'", this_module, 'matvec')
            end select
        class default
            call stop_error("The intent [IN] argument 'vec_in' must be of type 'dense_vector'", this_module, 'matvec')
        end select
         return
    end subroutine

    
    function initialize_dense_linop_from_array_rsp(A) result(linop)
        real(sp), intent(in) :: A(:, :)
        type(dense_linop_rsp) :: linop
        linop%data = A
        return
    end function
    function initialize_dense_linop_from_array_rdp(A) result(linop)
        real(dp), intent(in) :: A(:, :)
        type(dense_linop_rdp) :: linop
        linop%data = A
        return
    end function
    function initialize_dense_linop_from_array_csp(A) result(linop)
        complex(sp), intent(in) :: A(:, :)
        type(dense_linop_csp) :: linop
        linop%data = A
        return
    end function
    function initialize_dense_linop_from_array_cdp(A) result(linop)
        complex(dp), intent(in) :: A(:, :)
        type(dense_linop_cdp) :: linop
        linop%data = A
        return
    end function

end module LightKrylov_AbstractLinops
