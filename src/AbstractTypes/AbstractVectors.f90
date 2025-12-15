module LightKrylov_AbstractVectors
    !! This module provides the base class `absract_vector` from which all Krylov vectors
    !! needs to be derived. To use `LightKrylov`, you need to extend one of the
    !! followings:
    !!
    !! - `abstract_vector_rsp`  :   Real-valued vector with single precision arithmetic.
    !! - `abstract_vector_rdp`  :   Real-valued vector with double precision arithmetic.
    !! - `abstract_vector_csp`  :   Complex-valued vector with single precision arithmetic.
    !! - `abstract_vector_cdp`  :   Complex-valued vector with double precision arithmetic.
    !!
    !! To extend either of these abstract types, you need to provide an associated implementation
    !! for the following type-bound procedures:
    !!
    !! - `zero(self)`                   :   A subroutine zeroing-out the vector.
    !! - `rand(self, ifnorm)`           :   A subroutine creating a random vector, possibily normalized to have unit-norm (`ifnorm = .true.`).
    !! - `scal(self, alpha)`            :   A subroutine computing *in-place* the scalar multiplication \( \mathbf{x} \leftarrow \alpha \mathbf{x} \).
    !! - `axpby(alpha, vec, beta, self) :   A subroutine computing *in-place* the product \( \mathbf{y} \leftarrow \alpha \mathbf{x} + \beta \mathbf{y} \).
    !! - `dot(self, vec)`               :   A function computing the inner product \( \alpha = \langle \mathbf{x} \vert \mathbf{y} \rangle \).
    !! - `get_size(self)`               :   A function returning the dimension \( n \) of the vector \( \mathbf{x} \).
    !!
    !! Once these type-bound procedures have been implemented by the user, they will automatically 
    !! be used to define:
    !!
    !! - vector addition    :   `add(self, vec) = axpby(1, vec, 1, self)`
    !! - vector subtraction :   `sub(self, vec) = axpby(-1, vec, 1, self)`
    !! - vector norm        :   `norm(self)     = sqrt(self%dot(self))`
    !!
    !! This module also provides the following utility subroutines:
    !!
    !! - `innerprod(X, Y)`                  : Function computing the product \(\mathbf{X}^H \mathbf{y} \) between a Krylov basis `X` and a Krylov vector  (resp. basis) `Y`.
    !! - `linear_combination(Y, X, V)`      : Subroutine computing the linear combination \( \mathbf{y}_j = \sum_{i=1}^n \mathbf{x}_i v_{ij} \).
    !! - `axpby_basis(alpha, X, beta, Y)`   : In-place computation of \( \mathbf{Y} \leftarrow \alpha \mathbf{X} + \beta \mathbf{Y} \) where `X` and `Y` are arrays of `abstract_vector`.
    !! - `zero_basis(X)`                    : Zero-out a collection of `abstract_vectors`.
    !! - `copy_basis(out, from)`            : Copy a collection of `abstract_vectors`.
    !! - `rand_basis(X, ifnorm)`            : Create a collection of random `abstract_vectors`. If `ifnorm = .true.`, the vectors are normalized to have unit-norm.

    use stdlib_optval, only: optval
    use stdlib_linalg_blas, only: scal, axpy, dot, dotc
    use LightKrylov_Constants
    use LightKrylov_Utils
    use LightKrylov_Logger
    implicit none(type, external)
    private

    character(len=*), parameter :: this_module      = 'LK_Vectors'
    character(len=*), parameter :: this_module_long = 'Lightkrylov_AbstractVectors'

    public :: innerprod, Gram
    public :: linear_combination
    public :: axpby_basis
    public :: zero_basis
    public :: copy
    public :: rand_basis

    interface innerprod
        !!  Compute the inner product vector \( \mathbf{v} = \mathbf{X}^H \mathbf{y} \) or matrix
        !!  \( \mathbf{M} = \mathbf{X}^H \mathbf{Y} \).
        !!
        !!  ### Description
        !!
        !!  This interface provides methods for computing the inner products between a basis
        !!  of `real` or `complex` vectors \( \mathbf{X} \) and a single vector 
        !!  \( \mathbf{y} \) or another basis \( \mathbf{Y} \). Depending on the case, it
        !!  returns a one-dimensional array \( \mathbf{v} \) or a two-dimensional array
        !!  \( \mathbf{M} \) with the same type as \( \mathbf{X} \).
        !!
        !!  ### Example
        !!
        !!  The example below assumes that you have already extended the `abstract_vector_rdp`
        !!  class to define your own `my_real_vector` type.
        !!
        !!  ```fortran
        !!      type(my_real_vector), dimension(10) :: X
        !!      type(my_real_vector)                :: y
        !!      real(dp), dimension(:), allocatable :: v
        !!
        !!      ! ... Part of your code where you initialize everything ...
        !!
        !!      v = innerprod(X, y)
        !!
        !!      ! ... Rest of your code ...
        !!  ```
        !!
        !!  Similarly, for computing the matrix of inner products between two bases
        !!
        !!  ```fortran
        !!      type(my_real_vector), dimension(10) :: X
        !!      type(my_real_vector), dimension(10) :: Y
        !!      real(dp), dimension(:, :), allocatable :: M
        !!
        !!      ! ... Part of your code where you initialize everything ...
        !!
        !!      M = innerprod(X, Y)
        !!
        !!      ! ... Rest of your code ...
        !!  ```
        module procedure innerprod_vector_rsp
        module procedure innerprod_matrix_rsp
        module procedure innerprod_vector_rdp
        module procedure innerprod_matrix_rdp
        module procedure innerprod_vector_csp
        module procedure innerprod_matrix_csp
        module procedure innerprod_vector_cdp
        module procedure innerprod_matrix_cdp
    end interface

    interface Gram
        !!  Compute the Gram matrix \( \mathbf{G} = \mathbf{X}^H \mathbf{X} \).
        !!
        !!  ### Description
        !!
        !!  This interface provides methods for computing the Gram matrix associated to a basis of 
        !!  `abstract_vector` \( \mathbf{X} \).
        !!
        !!  ### Example
        !!
        !!  The example below assumes that you have already extended the `abstract_vector_rdp`
        !!  class to define your own `my_real_vector` type.
        !!
        !!  ```fortran
        !!      type(my_real_vector), dimension(10) :: X
        !!      real(dp), dimension(:, :), allocatable :: G
        !!
        !!      ! ... Part of your code where you initialize everything ...
        !!
        !!      G = Gram(X)
        !!
        !!      ! ... Rest of your code ...
        !!  ```
        module procedure gram_matrix_rsp
        module procedure gram_matrix_rdp
        module procedure gram_matrix_csp
        module procedure gram_matrix_cdp
    end interface

    interface linear_combination
        !!  Given a set of extended `abstract_vectors` and coefficients, return the corresponding
        !!  linear combinations.
        !!
        !!  ### Description
        !!
        !!  This interface provides methods for computing linear combinations of a set of
        !!  `abstract_vectors`. Depending on its input, it either computes
        !!
        !!  \[
        !!      \mathbf{y} = \sum_{i=1}^n \alpha_i \mathbf{x}_i,
        !!  \]
        !!
        !!  i.e. a single vector, or
        !!
        !!  \[
        !!      \mathbf{y}_j = \sum_{i=1}^n \alpha_{ij} \mathbf{x}_i,
        !!  \]
        !!
        !!  i.e. a set of vectors of the same type as \( \mathbf{X} \).
        !!
        !!  ### Example
        !!
        !!  ```fortran
        !!      type(my_real_vector), dimension(10) :: X
        !!      real(dp), dimension(m, n)           :: B
        !!      type(my_real_vector)                :: Y
        !!
        !!      ! ... Whatever your code is doing ...
        !!
        !!      call linear_combination(Y, X, B)
        !!
        !!      ! ... Rest of your code ...
        !!  ```
        module procedure linear_combination_vector_rsp
        module procedure linear_combination_matrix_rsp
        module procedure linear_combination_vector_rdp
        module procedure linear_combination_matrix_rdp
        module procedure linear_combination_vector_csp
        module procedure linear_combination_matrix_csp
        module procedure linear_combination_vector_cdp
        module procedure linear_combination_matrix_cdp
    end interface

    interface axpby_basis
        !!  In-place addition of two arrays of extended `abstract_vector`.
        !!  
        !!  ### Description
        !!
        !!  This interface provides methods to add in-place two arrays of
        !!  extended `abstract_vector`, i.e.
        !!
        !!  \[
        !!      \mathbf{Y}_i \leftarrow \alpha \mathbf{X}_i + \beta \mathbf{Y}_i.
        !!  \]
        !!
        !!  No out-of-place alternative is currently available in `LightKrylov`.
        !!  If you do need an out-of-place version, you can combine `axpby_basis`
        !!  with `copy`.
        !!
        !!  ### Example
        !!
        !!  ```fortran
        !!      type(my_real_vector), dimension(10) :: X
        !!      type(my_real_vector), dimension(10) :: Y
        !!      real(dp), dimension(10)             :: alpha, beta
        !!
        !!      ! ... Whatever your code is doing ...
        !!
        !!      call axpby_basis(alpha, X, beta, Y)
        !!
        !!      ! ... Rest of your code ...
        !!  ```
        module procedure axpby_basis_rsp
        module procedure axpby_basis_rdp
        module procedure axpby_basis_csp
        module procedure axpby_basis_cdp
    end interface

    interface zero_basis
        !!  This interface provides methods to zero-out a collection of `abstract_vector` `X`.
        !!  It is a simple wrapper around `X(i)%zero()`.
        !!
        !!  ### Example
        !!
        !!  ```fortran
        !!      type(my_real_vector), dimension(10) :: X
        !!
        !!      ! ... Your code ...
        !!
        !!      call zero_basis(X)
        !!
        !!      ! ... Your code ...
        !!  ```
        module procedure zero_basis_rsp
        module procedure zero_basis_rdp
        module procedure zero_basis_csp
        module procedure zero_basis_cdp
    end interface

    interface copy
        !!  This interface provides methods to copy an array `X` of `abstract_vector` into
        !!  another array `Y`. Note that `Y` needs to be pre-allocated.
        !!
        !!  ### Example
        !!
        !!  ```fortran
        !!      type(my_real_vector), dimension(10) :: X
        !!      type(my_real_vector), dimension(10) :: Y
        !!
        !!      ! ... Your code ...
        !!
        !!      call copy(Y, X)
        !!
        !!      ! ... Your code ...
        !!  ```
        module procedure copy_vector_rsp
        ! module procedure copy_basis_rsp
        module procedure copy_vector_rdp
        ! module procedure copy_basis_rdp
        module procedure copy_vector_csp
        ! module procedure copy_basis_csp
        module procedure copy_vector_cdp
        ! module procedure copy_basis_cdp
    end interface

    interface rand_basis
        !!  This interface provides methods to create an array `X` of random `abstract_vector`.
        !!  It is a simple wrapper around `X(i)%rand(ifnorm)`.
        !!
        !!  ### Example
        !!
        !!  ```fortran
        !!      type(my_real_vector), dimension(10) :: X
        !!      logical                             :: ifnorm = .true.
        !!
        !!      ! ... Your code ...
        !!
        !!      call rand_basis(X, ifnorm)
        !!
        !!      ! ... Your code ...
        !!  ```
        module procedure rand_basis_rsp
        module procedure rand_basis_rdp
        module procedure rand_basis_csp
        module procedure rand_basis_cdp
    end interface

    type, abstract, public :: abstract_vector
        !!  Base abstract type from which all other types of vectors used in `LightKrylov`
        !!  are being derived from.
        !!
        !!  @warning
        !!  Users should not extend this abstract class to define their own types.
        !!  @endwarning
    end type abstract_vector

    !----------------------------------------------------------------------------
    !-----     Definition of an abstract real(sp) vector with kind=sp     -----
    !----------------------------------------------------------------------------

    type, abstract, extends(abstract_vector), public :: abstract_vector_rsp
        !!  Abstract type to define real(sp)-valued vectors.
        !!  Derived-types defined by the user should be extending one such class.
    contains
        private
        procedure(abstract_zero_rsp), pass(self), deferred, public :: zero
        !! Sets an `abstract_vector_rsp` to zero.
        procedure(abstract_rand_rsp), pass(self), deferred, public :: rand
        !! Creates a random `abstract_vector_rsp`.
        procedure(abstract_scal_rsp), pass(self), deferred, public :: scal
        !! Compute the scalar-vector product.
        procedure(abstract_axpby_rsp), pass(self), deferred, public :: axpby
        !! In-place computation of \( \mathbf{y} \leftarrow \alpha \mathbf{x} + \beta \mathbf{y} \).
        procedure(abstract_dot_rsp), pass(self), deferred, public :: dot
        !! Computes the dot product between two `abstract_vector_rsp`.
        procedure(abstract_get_size_rsp), pass(self), deferred, public :: get_size
        !! Return size of specific abstract vector
        procedure, pass(self), public :: norm => norm_rsp
        !! Computes the norm of the `abstract_vector`.
        procedure, pass(self), public :: add => add_rsp
        !! Adds two `abstract_vector`, i.e. \( \mathbf{y} \leftarrow \mathbf{x} + \mathbf{y}\).
        procedure, pass(self), public :: sub => sub_rsp
        !! Subtracts two `abstract_vector`, i.e. \( \mathbf{y} \leftarrow \mathbf{y} - \mathbf{x} \).
        procedure, pass(self), public :: chsgn => chsgn_rsp
        !! Change the sign of a vector, i.e. \( \mathbf{x} \leftarrow -\mathbf{x} \).
    end type abstract_vector_rsp

    abstract interface
        subroutine abstract_zero_rsp(self)
            !! Abstract interface to zero-out a vector in-place.
            import abstract_vector_rsp
            implicit none(type, external)
            class(abstract_vector_rsp), intent(inout) :: self
            !! Vector to be zeroed-out.
        end subroutine abstract_zero_rsp

        subroutine abstract_rand_rsp(self, ifnorm)
            !! Abstract interface to generate a random (normalized) vector.
            import abstract_vector_rsp
            implicit none(type, external)
            class(abstract_vector_rsp), intent(inout) :: self
            !! Vector to be initialized.
            logical, optional, intent(in) :: ifnorm
        end subroutine abstract_rand_rsp

        subroutine abstract_scal_rsp(self, alpha)
            !! Abstract interface to scale a vector.
            import abstract_vector_rsp, sp
            implicit none(type, external)
            class(abstract_vector_rsp), intent(inout) :: self
            !! Vector to be scaled.
            real(sp), intent(in) :: alpha
            !! Scaling factor.
        end subroutine abstract_scal_rsp

        subroutine abstract_axpby_rsp(alpha, vec, beta, self)
            !! Abstract interface to add/scale two vectors in-place.
            import abstract_vector_rsp, sp
            implicit none(type, external)
            class(abstract_vector_rsp), intent(inout) :: self
            !! Input/Output vector.
            class(abstract_vector_rsp), intent(in) :: vec
            !! Vector to be added/subtracted.
            real(sp), intent(in) :: alpha, beta
        end subroutine abstract_axpby_rsp

        function abstract_dot_rsp(self, vec) result(alpha)
            !! Abstract interface to compute the dot product.
            import abstract_vector_rsp, sp
            implicit none(type, external)
            class(abstract_vector_rsp), intent(in) :: self, vec
            !! Vectors whose dot product will be computed.
            real(sp) :: alpha
            !! Result of the dot product.
        end function abstract_dot_rsp

        function abstract_get_size_rsp(self) result(N)
            !! Abstract interface to return the size of the specific abstract vector.
            import abstract_vector_rsp
            implicit none(type, external)
            class(abstract_vector_rsp), intent(in) :: self
            !! Vector for which to return the size.
            integer :: N
            !! Size of the vector
        end function abstract_get_size_rsp

    end interface

    !----------------------------------------------------------------------------
    !-----     Definition of an abstract real(dp) vector with kind=dp     -----
    !----------------------------------------------------------------------------

    type, abstract, extends(abstract_vector), public :: abstract_vector_rdp
        !!  Abstract type to define real(dp)-valued vectors.
        !!  Derived-types defined by the user should be extending one such class.
    contains
        private
        procedure(abstract_zero_rdp), pass(self), deferred, public :: zero
        !! Sets an `abstract_vector_rdp` to zero.
        procedure(abstract_rand_rdp), pass(self), deferred, public :: rand
        !! Creates a random `abstract_vector_rdp`.
        procedure(abstract_scal_rdp), pass(self), deferred, public :: scal
        !! Compute the scalar-vector product.
        procedure(abstract_axpby_rdp), pass(self), deferred, public :: axpby
        !! In-place computation of \( \mathbf{y} \leftarrow \alpha \mathbf{x} + \beta \mathbf{y} \).
        procedure(abstract_dot_rdp), pass(self), deferred, public :: dot
        !! Computes the dot product between two `abstract_vector_rdp`.
        procedure(abstract_get_size_rdp), pass(self), deferred, public :: get_size
        !! Return size of specific abstract vector
        procedure, pass(self), public :: norm => norm_rdp
        !! Computes the norm of the `abstract_vector`.
        procedure, pass(self), public :: add => add_rdp
        !! Adds two `abstract_vector`, i.e. \( \mathbf{y} \leftarrow \mathbf{x} + \mathbf{y}\).
        procedure, pass(self), public :: sub => sub_rdp
        !! Subtracts two `abstract_vector`, i.e. \( \mathbf{y} \leftarrow \mathbf{y} - \mathbf{x} \).
        procedure, pass(self), public :: chsgn => chsgn_rdp
        !! Change the sign of a vector, i.e. \( \mathbf{x} \leftarrow -\mathbf{x} \).
    end type abstract_vector_rdp

    abstract interface
        subroutine abstract_zero_rdp(self)
            !! Abstract interface to zero-out a vector in-place.
            import abstract_vector_rdp
            implicit none(type, external)
            class(abstract_vector_rdp), intent(inout) :: self
            !! Vector to be zeroed-out.
        end subroutine abstract_zero_rdp

        subroutine abstract_rand_rdp(self, ifnorm)
            !! Abstract interface to generate a random (normalized) vector.
            import abstract_vector_rdp
            implicit none(type, external)
            class(abstract_vector_rdp), intent(inout) :: self
            !! Vector to be initialized.
            logical, optional, intent(in) :: ifnorm
        end subroutine abstract_rand_rdp

        subroutine abstract_scal_rdp(self, alpha)
            !! Abstract interface to scale a vector.
            import abstract_vector_rdp, dp
            implicit none(type, external)
            class(abstract_vector_rdp), intent(inout) :: self
            !! Vector to be scaled.
            real(dp), intent(in) :: alpha
            !! Scaling factor.
        end subroutine abstract_scal_rdp

        subroutine abstract_axpby_rdp(alpha, vec, beta, self)
            !! Abstract interface to add/scale two vectors in-place.
            import abstract_vector_rdp, dp
            implicit none(type, external)
            class(abstract_vector_rdp), intent(inout) :: self
            !! Input/Output vector.
            class(abstract_vector_rdp), intent(in) :: vec
            !! Vector to be added/subtracted.
            real(dp), intent(in) :: alpha, beta
        end subroutine abstract_axpby_rdp

        function abstract_dot_rdp(self, vec) result(alpha)
            !! Abstract interface to compute the dot product.
            import abstract_vector_rdp, dp
            implicit none(type, external)
            class(abstract_vector_rdp), intent(in) :: self, vec
            !! Vectors whose dot product will be computed.
            real(dp) :: alpha
            !! Result of the dot product.
        end function abstract_dot_rdp

        function abstract_get_size_rdp(self) result(N)
            !! Abstract interface to return the size of the specific abstract vector.
            import abstract_vector_rdp
            implicit none(type, external)
            class(abstract_vector_rdp), intent(in) :: self
            !! Vector for which to return the size.
            integer :: N
            !! Size of the vector
        end function abstract_get_size_rdp

    end interface

    !----------------------------------------------------------------------------
    !-----     Definition of an abstract complex(sp) vector with kind=sp     -----
    !----------------------------------------------------------------------------

    type, abstract, extends(abstract_vector), public :: abstract_vector_csp
        !!  Abstract type to define complex(sp)-valued vectors.
        !!  Derived-types defined by the user should be extending one such class.
    contains
        private
        procedure(abstract_zero_csp), pass(self), deferred, public :: zero
        !! Sets an `abstract_vector_csp` to zero.
        procedure(abstract_rand_csp), pass(self), deferred, public :: rand
        !! Creates a random `abstract_vector_csp`.
        procedure(abstract_scal_csp), pass(self), deferred, public :: scal
        !! Compute the scalar-vector product.
        procedure(abstract_axpby_csp), pass(self), deferred, public :: axpby
        !! In-place computation of \( \mathbf{y} \leftarrow \alpha \mathbf{x} + \beta \mathbf{y} \).
        procedure(abstract_dot_csp), pass(self), deferred, public :: dot
        !! Computes the dot product between two `abstract_vector_csp`.
        procedure(abstract_get_size_csp), pass(self), deferred, public :: get_size
        !! Return size of specific abstract vector
        procedure, pass(self), public :: norm => norm_csp
        !! Computes the norm of the `abstract_vector`.
        procedure, pass(self), public :: add => add_csp
        !! Adds two `abstract_vector`, i.e. \( \mathbf{y} \leftarrow \mathbf{x} + \mathbf{y}\).
        procedure, pass(self), public :: sub => sub_csp
        !! Subtracts two `abstract_vector`, i.e. \( \mathbf{y} \leftarrow \mathbf{y} - \mathbf{x} \).
        procedure, pass(self), public :: chsgn => chsgn_csp
        !! Change the sign of a vector, i.e. \( \mathbf{x} \leftarrow -\mathbf{x} \).
    end type abstract_vector_csp

    abstract interface
        subroutine abstract_zero_csp(self)
            !! Abstract interface to zero-out a vector in-place.
            import abstract_vector_csp
            implicit none(type, external)
            class(abstract_vector_csp), intent(inout) :: self
            !! Vector to be zeroed-out.
        end subroutine abstract_zero_csp

        subroutine abstract_rand_csp(self, ifnorm)
            !! Abstract interface to generate a random (normalized) vector.
            import abstract_vector_csp
            implicit none(type, external)
            class(abstract_vector_csp), intent(inout) :: self
            !! Vector to be initialized.
            logical, optional, intent(in) :: ifnorm
        end subroutine abstract_rand_csp

        subroutine abstract_scal_csp(self, alpha)
            !! Abstract interface to scale a vector.
            import abstract_vector_csp, sp
            implicit none(type, external)
            class(abstract_vector_csp), intent(inout) :: self
            !! Vector to be scaled.
            complex(sp), intent(in) :: alpha
            !! Scaling factor.
        end subroutine abstract_scal_csp

        subroutine abstract_axpby_csp(alpha, vec, beta, self)
            !! Abstract interface to add/scale two vectors in-place.
            import abstract_vector_csp, sp
            implicit none(type, external)
            class(abstract_vector_csp), intent(inout) :: self
            !! Input/Output vector.
            class(abstract_vector_csp), intent(in) :: vec
            !! Vector to be added/subtracted.
            complex(sp), intent(in) :: alpha, beta
        end subroutine abstract_axpby_csp

        function abstract_dot_csp(self, vec) result(alpha)
            !! Abstract interface to compute the dot product.
            import abstract_vector_csp, sp
            implicit none(type, external)
            class(abstract_vector_csp), intent(in) :: self, vec
            !! Vectors whose dot product will be computed.
            complex(sp) :: alpha
            !! Result of the dot product.
        end function abstract_dot_csp

        function abstract_get_size_csp(self) result(N)
            !! Abstract interface to return the size of the specific abstract vector.
            import abstract_vector_csp
            implicit none(type, external)
            class(abstract_vector_csp), intent(in) :: self
            !! Vector for which to return the size.
            integer :: N
            !! Size of the vector
        end function abstract_get_size_csp

    end interface

    !----------------------------------------------------------------------------
    !-----     Definition of an abstract complex(dp) vector with kind=dp     -----
    !----------------------------------------------------------------------------

    type, abstract, extends(abstract_vector), public :: abstract_vector_cdp
        !!  Abstract type to define complex(dp)-valued vectors.
        !!  Derived-types defined by the user should be extending one such class.
    contains
        private
        procedure(abstract_zero_cdp), pass(self), deferred, public :: zero
        !! Sets an `abstract_vector_cdp` to zero.
        procedure(abstract_rand_cdp), pass(self), deferred, public :: rand
        !! Creates a random `abstract_vector_cdp`.
        procedure(abstract_scal_cdp), pass(self), deferred, public :: scal
        !! Compute the scalar-vector product.
        procedure(abstract_axpby_cdp), pass(self), deferred, public :: axpby
        !! In-place computation of \( \mathbf{y} \leftarrow \alpha \mathbf{x} + \beta \mathbf{y} \).
        procedure(abstract_dot_cdp), pass(self), deferred, public :: dot
        !! Computes the dot product between two `abstract_vector_cdp`.
        procedure(abstract_get_size_cdp), pass(self), deferred, public :: get_size
        !! Return size of specific abstract vector
        procedure, pass(self), public :: norm => norm_cdp
        !! Computes the norm of the `abstract_vector`.
        procedure, pass(self), public :: add => add_cdp
        !! Adds two `abstract_vector`, i.e. \( \mathbf{y} \leftarrow \mathbf{x} + \mathbf{y}\).
        procedure, pass(self), public :: sub => sub_cdp
        !! Subtracts two `abstract_vector`, i.e. \( \mathbf{y} \leftarrow \mathbf{y} - \mathbf{x} \).
        procedure, pass(self), public :: chsgn => chsgn_cdp
        !! Change the sign of a vector, i.e. \( \mathbf{x} \leftarrow -\mathbf{x} \).
    end type abstract_vector_cdp

    abstract interface
        subroutine abstract_zero_cdp(self)
            !! Abstract interface to zero-out a vector in-place.
            import abstract_vector_cdp
            implicit none(type, external)
            class(abstract_vector_cdp), intent(inout) :: self
            !! Vector to be zeroed-out.
        end subroutine abstract_zero_cdp

        subroutine abstract_rand_cdp(self, ifnorm)
            !! Abstract interface to generate a random (normalized) vector.
            import abstract_vector_cdp
            implicit none(type, external)
            class(abstract_vector_cdp), intent(inout) :: self
            !! Vector to be initialized.
            logical, optional, intent(in) :: ifnorm
        end subroutine abstract_rand_cdp

        subroutine abstract_scal_cdp(self, alpha)
            !! Abstract interface to scale a vector.
            import abstract_vector_cdp, dp
            implicit none(type, external)
            class(abstract_vector_cdp), intent(inout) :: self
            !! Vector to be scaled.
            complex(dp), intent(in) :: alpha
            !! Scaling factor.
        end subroutine abstract_scal_cdp

        subroutine abstract_axpby_cdp(alpha, vec, beta, self)
            !! Abstract interface to add/scale two vectors in-place.
            import abstract_vector_cdp, dp
            implicit none(type, external)
            class(abstract_vector_cdp), intent(inout) :: self
            !! Input/Output vector.
            class(abstract_vector_cdp), intent(in) :: vec
            !! Vector to be added/subtracted.
            complex(dp), intent(in) :: alpha, beta
        end subroutine abstract_axpby_cdp

        function abstract_dot_cdp(self, vec) result(alpha)
            !! Abstract interface to compute the dot product.
            import abstract_vector_cdp, dp
            implicit none(type, external)
            class(abstract_vector_cdp), intent(in) :: self, vec
            !! Vectors whose dot product will be computed.
            complex(dp) :: alpha
            !! Result of the dot product.
        end function abstract_dot_cdp

        function abstract_get_size_cdp(self) result(N)
            !! Abstract interface to return the size of the specific abstract vector.
            import abstract_vector_cdp
            implicit none(type, external)
            class(abstract_vector_cdp), intent(in) :: self
            !! Vector for which to return the size.
            integer :: N
            !! Size of the vector
        end function abstract_get_size_cdp

    end interface


    !----------------------------------------------------------------------------------
    !-----     Convenience vector type to wrap standard Fortran rank-1 arrays     -----
    !----------------------------------------------------------------------------------

    type, extends(abstract_vector_rsp), public :: dense_vector_rsp
        integer :: n
        real(sp), allocatable :: data(:)
    contains
        private
        procedure, pass(self), public :: zero => dense_zero_rsp
        !! Sets an `abstract_vector_rsp` to zero.
        procedure, pass(self), public :: rand => dense_rand_rsp
        !! Creates a random `abstract_vector_rsp`.
        procedure, pass(self), public :: scal => dense_scal_rsp
        !! Compute the scalar-vector product.
        procedure, pass(self), public :: axpby => dense_axpby_rsp
        !! In-place computation of \( \mathbf{y} \leftarrow \alpha \mathbf{x} + \beta \mathbf{y} \).
        procedure, pass(self), public :: dot => dense_dot_rsp
        !! Computes the dot product between two `abstract_vector_rsp`.
        procedure, pass(self), public :: get_size => dense_get_size_rsp
        !! Return size of specific abstract vector
    end type dense_vector_rsp
    !----------------------------------------------------------------------------------
    !-----     Convenience vector type to wrap standard Fortran rank-1 arrays     -----
    !----------------------------------------------------------------------------------

    type, extends(abstract_vector_rdp), public :: dense_vector_rdp
        integer :: n
        real(dp), allocatable :: data(:)
    contains
        private
        procedure, pass(self), public :: zero => dense_zero_rdp
        !! Sets an `abstract_vector_rdp` to zero.
        procedure, pass(self), public :: rand => dense_rand_rdp
        !! Creates a random `abstract_vector_rdp`.
        procedure, pass(self), public :: scal => dense_scal_rdp
        !! Compute the scalar-vector product.
        procedure, pass(self), public :: axpby => dense_axpby_rdp
        !! In-place computation of \( \mathbf{y} \leftarrow \alpha \mathbf{x} + \beta \mathbf{y} \).
        procedure, pass(self), public :: dot => dense_dot_rdp
        !! Computes the dot product between two `abstract_vector_rdp`.
        procedure, pass(self), public :: get_size => dense_get_size_rdp
        !! Return size of specific abstract vector
    end type dense_vector_rdp
    !----------------------------------------------------------------------------------
    !-----     Convenience vector type to wrap standard Fortran rank-1 arrays     -----
    !----------------------------------------------------------------------------------

    type, extends(abstract_vector_csp), public :: dense_vector_csp
        integer :: n
        complex(sp), allocatable :: data(:)
    contains
        private
        procedure, pass(self), public :: zero => dense_zero_csp
        !! Sets an `abstract_vector_csp` to zero.
        procedure, pass(self), public :: rand => dense_rand_csp
        !! Creates a random `abstract_vector_csp`.
        procedure, pass(self), public :: scal => dense_scal_csp
        !! Compute the scalar-vector product.
        procedure, pass(self), public :: axpby => dense_axpby_csp
        !! In-place computation of \( \mathbf{y} \leftarrow \alpha \mathbf{x} + \beta \mathbf{y} \).
        procedure, pass(self), public :: dot => dense_dot_csp
        !! Computes the dot product between two `abstract_vector_csp`.
        procedure, pass(self), public :: get_size => dense_get_size_csp
        !! Return size of specific abstract vector
    end type dense_vector_csp
    !----------------------------------------------------------------------------------
    !-----     Convenience vector type to wrap standard Fortran rank-1 arrays     -----
    !----------------------------------------------------------------------------------

    type, extends(abstract_vector_cdp), public :: dense_vector_cdp
        integer :: n
        complex(dp), allocatable :: data(:)
    contains
        private
        procedure, pass(self), public :: zero => dense_zero_cdp
        !! Sets an `abstract_vector_cdp` to zero.
        procedure, pass(self), public :: rand => dense_rand_cdp
        !! Creates a random `abstract_vector_cdp`.
        procedure, pass(self), public :: scal => dense_scal_cdp
        !! Compute the scalar-vector product.
        procedure, pass(self), public :: axpby => dense_axpby_cdp
        !! In-place computation of \( \mathbf{y} \leftarrow \alpha \mathbf{x} + \beta \mathbf{y} \).
        procedure, pass(self), public :: dot => dense_dot_cdp
        !! Computes the dot product between two `abstract_vector_cdp`.
        procedure, pass(self), public :: get_size => dense_get_size_cdp
        !! Return size of specific abstract vector
    end type dense_vector_cdp

    interface dense_vector
        module procedure initialize_dense_vector_from_array_rsp
        module procedure initialize_dense_vector_from_array_rdp
        module procedure initialize_dense_vector_from_array_csp
        module procedure initialize_dense_vector_from_array_cdp
    end interface
    public :: dense_vector

contains

    !-----------------------------------------------------------------------
    !-----     TYPE-BOUND PROCEDURES FOR THE ABSTRACT VECTOR TYPES     -----
    !-----------------------------------------------------------------------

    function norm_rsp(self) result(alpha)
        implicit none(type, external)
        !! Compute the norm of an `abstract_vector`.
        class(abstract_vector_rsp), intent(in) :: self
        !! Vector whose norm needs to be computed.
        real(sp) :: alpha
        !! Norm of the vector.
        alpha = abs(self%dot(self)) ; alpha = sqrt(alpha)
    end function norm_rsp

    subroutine sub_rsp(self, vec)
        implicit none(type, external)
        !! Subtract two `abstract_vector` in-place.
        class(abstract_vector_rsp), intent(inout) :: self
        !! Input/Output vector.
        class(abstract_vector_rsp), intent(in) :: vec
        !! Vector to be subtracted.
        call self%axpby(-one_rsp, vec, one_rsp)
    end subroutine sub_rsp

    subroutine add_rsp(self, vec)
        implicit none(type, external)
        !! Add two `abstract_vector` in-place.
        class(abstract_vector_rsp), intent(inout) :: self
        !! Input/Output vector.
        class(abstract_vector_rsp), intent(in) :: vec
        !! Vector to be added.
        call self%axpby(one_rsp, vec, one_rsp)
    end subroutine add_rsp

    subroutine chsgn_rsp(self)
        implicit none(type, external)
        !! Changes the sign of the `abstract_vector`.
        class(abstract_vector_rsp), intent(inout) :: self
        !! Vector whose entries need to change sign.
        call self%scal(-one_rsp)
    end subroutine chsgn_rsp

    function norm_rdp(self) result(alpha)
        implicit none(type, external)
        !! Compute the norm of an `abstract_vector`.
        class(abstract_vector_rdp), intent(in) :: self
        !! Vector whose norm needs to be computed.
        real(dp) :: alpha
        !! Norm of the vector.
        alpha = abs(self%dot(self)) ; alpha = sqrt(alpha)
    end function norm_rdp

    subroutine sub_rdp(self, vec)
        implicit none(type, external)
        !! Subtract two `abstract_vector` in-place.
        class(abstract_vector_rdp), intent(inout) :: self
        !! Input/Output vector.
        class(abstract_vector_rdp), intent(in) :: vec
        !! Vector to be subtracted.
        call self%axpby(-one_rdp, vec, one_rdp)
    end subroutine sub_rdp

    subroutine add_rdp(self, vec)
        implicit none(type, external)
        !! Add two `abstract_vector` in-place.
        class(abstract_vector_rdp), intent(inout) :: self
        !! Input/Output vector.
        class(abstract_vector_rdp), intent(in) :: vec
        !! Vector to be added.
        call self%axpby(one_rdp, vec, one_rdp)
    end subroutine add_rdp

    subroutine chsgn_rdp(self)
        implicit none(type, external)
        !! Changes the sign of the `abstract_vector`.
        class(abstract_vector_rdp), intent(inout) :: self
        !! Vector whose entries need to change sign.
        call self%scal(-one_rdp)
    end subroutine chsgn_rdp

    function norm_csp(self) result(alpha)
        implicit none(type, external)
        !! Compute the norm of an `abstract_vector`.
        class(abstract_vector_csp), intent(in) :: self
        !! Vector whose norm needs to be computed.
        real(sp) :: alpha
        !! Norm of the vector.
        alpha = abs(self%dot(self)) ; alpha = sqrt(alpha)
    end function norm_csp

    subroutine sub_csp(self, vec)
        implicit none(type, external)
        !! Subtract two `abstract_vector` in-place.
        class(abstract_vector_csp), intent(inout) :: self
        !! Input/Output vector.
        class(abstract_vector_csp), intent(in) :: vec
        !! Vector to be subtracted.
        call self%axpby(-one_csp, vec, one_csp)
    end subroutine sub_csp

    subroutine add_csp(self, vec)
        implicit none(type, external)
        !! Add two `abstract_vector` in-place.
        class(abstract_vector_csp), intent(inout) :: self
        !! Input/Output vector.
        class(abstract_vector_csp), intent(in) :: vec
        !! Vector to be added.
        call self%axpby(one_csp, vec, one_csp)
    end subroutine add_csp

    subroutine chsgn_csp(self)
        implicit none(type, external)
        !! Changes the sign of the `abstract_vector`.
        class(abstract_vector_csp), intent(inout) :: self
        !! Vector whose entries need to change sign.
        call self%scal(-one_csp)
    end subroutine chsgn_csp

    function norm_cdp(self) result(alpha)
        implicit none(type, external)
        !! Compute the norm of an `abstract_vector`.
        class(abstract_vector_cdp), intent(in) :: self
        !! Vector whose norm needs to be computed.
        real(dp) :: alpha
        !! Norm of the vector.
        alpha = abs(self%dot(self)) ; alpha = sqrt(alpha)
    end function norm_cdp

    subroutine sub_cdp(self, vec)
        implicit none(type, external)
        !! Subtract two `abstract_vector` in-place.
        class(abstract_vector_cdp), intent(inout) :: self
        !! Input/Output vector.
        class(abstract_vector_cdp), intent(in) :: vec
        !! Vector to be subtracted.
        call self%axpby(-one_cdp, vec, one_cdp)
    end subroutine sub_cdp

    subroutine add_cdp(self, vec)
        implicit none(type, external)
        !! Add two `abstract_vector` in-place.
        class(abstract_vector_cdp), intent(inout) :: self
        !! Input/Output vector.
        class(abstract_vector_cdp), intent(in) :: vec
        !! Vector to be added.
        call self%axpby(one_cdp, vec, one_cdp)
    end subroutine add_cdp

    subroutine chsgn_cdp(self)
        implicit none(type, external)
        !! Changes the sign of the `abstract_vector`.
        class(abstract_vector_cdp), intent(inout) :: self
        !! Vector whose entries need to change sign.
        call self%scal(-one_cdp)
    end subroutine chsgn_cdp


    !--------------------------------------------------------------------------------
    !-----     TYPE-BOUND PROCEDURES FOR THE CONVENIENCE DENSE VECTOR TYPES     -----
    !--------------------------------------------------------------------------------
    
    function initialize_dense_vector_from_array_rsp(x) result(vec)
        implicit none(type, external)
        real(sp), intent(in) :: x(:)
        type(dense_vector_rsp) :: vec
        vec%n = size(x) ; vec%data = x
    end function initialize_dense_vector_from_array_rsp

    subroutine dense_zero_rsp(self)
        implicit none(type, external)
        class(dense_vector_rsp), intent(inout) :: self
        integer :: iostat
        character(len=100) :: errmsg
        if(.not. allocated(self%data)) then
            allocate(self%data(self%n), stat=iostat, errmsg=errmsg)
            if (iostat /= 0) call stop_error(errmsg, &
                                             module=this_module, &
                                             procedure="dense_zero_rsp")
        endif
        self%data = zero_rsp
    end subroutine dense_zero_rsp

    subroutine dense_rand_rsp(self, ifnorm)
        implicit none(type, external)
        class(dense_vector_rsp), intent(inout) :: self
        logical, optional, intent(in) :: ifnorm
        integer :: iostat
        character(len=100) :: errmsg
        call random_number(self%data)
    end subroutine dense_rand_rsp

    subroutine dense_scal_rsp(self, alpha)
        implicit none(type, external)
        class(dense_vector_rsp), intent(inout) :: self
        real(sp), intent(in) :: alpha
        integer :: n
        n = self%get_size()
        call scal(n, alpha, self%data, 1)
    end subroutine dense_scal_rsp

    subroutine dense_axpby_rsp(alpha, vec, beta, self)
        implicit none(type, external)
        real(sp), intent(in) :: alpha, beta
        class(dense_vector_rsp), intent(inout) :: self
        class(abstract_vector_rsp), intent(in) :: vec
        integer :: n, m, iostat
        character(len=100) :: errmsg
        m = vec%get_size()
        if(.not. allocated(self%data)) then
            allocate(self%data(m), source=zero_rsp, stat=iostat, errmsg=errmsg)
            if (iostat /= 0) call stop_error(errmsg, &
                                             module=this_module, &
                                             procedure="dense_axpby_rsp")
        endif
        n = self%get_size()
        if (m /= n) call stop_error("Inconsistent size between the two vectors.")

        select type (vec)
        type is(dense_vector_rsp)
            if (beta /= zero_rsp) call self%scal(beta)
            call axpy(n, alpha, vec%data, 1, self%data, 1)
        class default
            call type_error('vec','dense_vector_rsp','IN',this_module,'dense_axpby_rsp')
        end select
    end subroutine dense_axpby_rsp

    function dense_dot_rsp(self, vec) result(alpha)
        implicit none(type, external)
        class(dense_vector_rsp), intent(in) :: self
        class(abstract_vector_rsp), intent(in) :: vec
        real(sp) :: alpha
        integer :: n
        n = self%get_size()
        select type (vec)
        type is(dense_vector_rsp)
            alpha = dot(n, self%data, 1, vec%data, 1)
        class default
            call type_error('vec','dense_vector_rsp','IN',this_module,'dense_dot_rsp')
        end select
    end function dense_dot_rsp

    function dense_get_size_rsp(self) result(n)
        implicit none(type, external)
        class(dense_vector_rsp), intent(in) :: self
        integer :: n
        n = size(self%data)
    end function dense_get_size_rsp

    function initialize_dense_vector_from_array_rdp(x) result(vec)
        implicit none(type, external)
        real(dp), intent(in) :: x(:)
        type(dense_vector_rdp) :: vec
        vec%n = size(x) ; vec%data = x
    end function initialize_dense_vector_from_array_rdp

    subroutine dense_zero_rdp(self)
        implicit none(type, external)
        class(dense_vector_rdp), intent(inout) :: self
        integer :: iostat
        character(len=100) :: errmsg
        if(.not. allocated(self%data)) then
            allocate(self%data(self%n), stat=iostat, errmsg=errmsg)
            if (iostat /= 0) call stop_error(errmsg, &
                                             module=this_module, &
                                             procedure="dense_zero_rdp")
        endif
        self%data = zero_rdp
    end subroutine dense_zero_rdp

    subroutine dense_rand_rdp(self, ifnorm)
        implicit none(type, external)
        class(dense_vector_rdp), intent(inout) :: self
        logical, optional, intent(in) :: ifnorm
        integer :: iostat
        character(len=100) :: errmsg
        call random_number(self%data)
    end subroutine dense_rand_rdp

    subroutine dense_scal_rdp(self, alpha)
        implicit none(type, external)
        class(dense_vector_rdp), intent(inout) :: self
        real(dp), intent(in) :: alpha
        integer :: n
        n = self%get_size()
        call scal(n, alpha, self%data, 1)
    end subroutine dense_scal_rdp

    subroutine dense_axpby_rdp(alpha, vec, beta, self)
        implicit none(type, external)
        real(dp), intent(in) :: alpha, beta
        class(dense_vector_rdp), intent(inout) :: self
        class(abstract_vector_rdp), intent(in) :: vec
        integer :: n, m, iostat
        character(len=100) :: errmsg
        m = vec%get_size()
        if(.not. allocated(self%data)) then
            allocate(self%data(m), source=zero_rdp, stat=iostat, errmsg=errmsg)
            if (iostat /= 0) call stop_error(errmsg, &
                                             module=this_module, &
                                             procedure="dense_axpby_rdp")
        endif
        n = self%get_size()
        if (m /= n) call stop_error("Inconsistent size between the two vectors.")

        select type (vec)
        type is(dense_vector_rdp)
            if (beta /= zero_rdp) call self%scal(beta)
            call axpy(n, alpha, vec%data, 1, self%data, 1)
        class default
            call type_error('vec','dense_vector_rdp','IN',this_module,'dense_axpby_rdp')
        end select
    end subroutine dense_axpby_rdp

    function dense_dot_rdp(self, vec) result(alpha)
        implicit none(type, external)
        class(dense_vector_rdp), intent(in) :: self
        class(abstract_vector_rdp), intent(in) :: vec
        real(dp) :: alpha
        integer :: n
        n = self%get_size()
        select type (vec)
        type is(dense_vector_rdp)
            alpha = dot(n, self%data, 1, vec%data, 1)
        class default
            call type_error('vec','dense_vector_rdp','IN',this_module,'dense_dot_rdp')
        end select
    end function dense_dot_rdp

    function dense_get_size_rdp(self) result(n)
        implicit none(type, external)
        class(dense_vector_rdp), intent(in) :: self
        integer :: n
        n = size(self%data)
    end function dense_get_size_rdp

    function initialize_dense_vector_from_array_csp(x) result(vec)
        implicit none(type, external)
        complex(sp), intent(in) :: x(:)
        type(dense_vector_csp) :: vec
        vec%n = size(x) ; vec%data = x
    end function initialize_dense_vector_from_array_csp

    subroutine dense_zero_csp(self)
        implicit none(type, external)
        class(dense_vector_csp), intent(inout) :: self
        integer :: iostat
        character(len=100) :: errmsg
        if(.not. allocated(self%data)) then
            allocate(self%data(self%n), stat=iostat, errmsg=errmsg)
            if (iostat /= 0) call stop_error(errmsg, &
                                             module=this_module, &
                                             procedure="dense_zero_csp")
        endif
        self%data = zero_csp
    end subroutine dense_zero_csp

    subroutine dense_rand_csp(self, ifnorm)
        implicit none(type, external)
        class(dense_vector_csp), intent(inout) :: self
        logical, optional, intent(in) :: ifnorm
        integer :: iostat
        character(len=100) :: errmsg
        real(sp), allocatable :: y(:, :)
        allocate(y(size(self%data), 2), stat=iostat, errmsg=errmsg)
        if (iostat /= 0) call stop_error(errmsg, &
                                         module=this_module, &
                                         procedure="dense_rand_csp")
        call random_number(y)
        self%data%re = y(:, 1) ; self%data%im = y(:, 2)
    end subroutine dense_rand_csp

    subroutine dense_scal_csp(self, alpha)
        implicit none(type, external)
        class(dense_vector_csp), intent(inout) :: self
        complex(sp), intent(in) :: alpha
        integer :: n
        n = self%get_size()
        call scal(n, alpha, self%data, 1)
    end subroutine dense_scal_csp

    subroutine dense_axpby_csp(alpha, vec, beta, self)
        implicit none(type, external)
        complex(sp), intent(in) :: alpha, beta
        class(dense_vector_csp), intent(inout) :: self
        class(abstract_vector_csp), intent(in) :: vec
        integer :: n, m, iostat
        character(len=100) :: errmsg
        m = vec%get_size()
        if(.not. allocated(self%data)) then
            allocate(self%data(m), source=zero_csp, stat=iostat, errmsg=errmsg)
            if (iostat /= 0) call stop_error(errmsg, &
                                             module=this_module, &
                                             procedure="dense_axpby_csp")
        endif
        n = self%get_size()
        if (m /= n) call stop_error("Inconsistent size between the two vectors.")

        select type (vec)
        type is(dense_vector_csp)
            if (beta /= zero_csp) call self%scal(beta)
            call axpy(n, alpha, vec%data, 1, self%data, 1)
        class default
            call type_error('vec','dense_vector_csp','IN',this_module,'dense_axpby_csp')
        end select
    end subroutine dense_axpby_csp

    function dense_dot_csp(self, vec) result(alpha)
        implicit none(type, external)
        class(dense_vector_csp), intent(in) :: self
        class(abstract_vector_csp), intent(in) :: vec
        complex(sp) :: alpha
        integer :: n
        n = self%get_size()
        select type (vec)
        type is(dense_vector_csp)
            alpha = dotc(n, self%data, 1, vec%data, 1)
        class default
            call type_error('vec','dense_vector_csp','IN',this_module,'dense_dot_csp')
        end select
    end function dense_dot_csp

    function dense_get_size_csp(self) result(n)
        implicit none(type, external)
        class(dense_vector_csp), intent(in) :: self
        integer :: n
        n = size(self%data)
    end function dense_get_size_csp

    function initialize_dense_vector_from_array_cdp(x) result(vec)
        implicit none(type, external)
        complex(dp), intent(in) :: x(:)
        type(dense_vector_cdp) :: vec
        vec%n = size(x) ; vec%data = x
    end function initialize_dense_vector_from_array_cdp

    subroutine dense_zero_cdp(self)
        implicit none(type, external)
        class(dense_vector_cdp), intent(inout) :: self
        integer :: iostat
        character(len=100) :: errmsg
        if(.not. allocated(self%data)) then
            allocate(self%data(self%n), stat=iostat, errmsg=errmsg)
            if (iostat /= 0) call stop_error(errmsg, &
                                             module=this_module, &
                                             procedure="dense_zero_cdp")
        endif
        self%data = zero_cdp
    end subroutine dense_zero_cdp

    subroutine dense_rand_cdp(self, ifnorm)
        implicit none(type, external)
        class(dense_vector_cdp), intent(inout) :: self
        logical, optional, intent(in) :: ifnorm
        integer :: iostat
        character(len=100) :: errmsg
        real(dp), allocatable :: y(:, :)
        allocate(y(size(self%data), 2), stat=iostat, errmsg=errmsg)
        if (iostat /= 0) call stop_error(errmsg, &
                                         module=this_module, &
                                         procedure="dense_rand_cdp")
        call random_number(y)
        self%data%re = y(:, 1) ; self%data%im = y(:, 2)
    end subroutine dense_rand_cdp

    subroutine dense_scal_cdp(self, alpha)
        implicit none(type, external)
        class(dense_vector_cdp), intent(inout) :: self
        complex(dp), intent(in) :: alpha
        integer :: n
        n = self%get_size()
        call scal(n, alpha, self%data, 1)
    end subroutine dense_scal_cdp

    subroutine dense_axpby_cdp(alpha, vec, beta, self)
        implicit none(type, external)
        complex(dp), intent(in) :: alpha, beta
        class(dense_vector_cdp), intent(inout) :: self
        class(abstract_vector_cdp), intent(in) :: vec
        integer :: n, m, iostat
        character(len=100) :: errmsg
        m = vec%get_size()
        if(.not. allocated(self%data)) then
            allocate(self%data(m), source=zero_cdp, stat=iostat, errmsg=errmsg)
            if (iostat /= 0) call stop_error(errmsg, &
                                             module=this_module, &
                                             procedure="dense_axpby_cdp")
        endif
        n = self%get_size()
        if (m /= n) call stop_error("Inconsistent size between the two vectors.")

        select type (vec)
        type is(dense_vector_cdp)
            if (beta /= zero_cdp) call self%scal(beta)
            call axpy(n, alpha, vec%data, 1, self%data, 1)
        class default
            call type_error('vec','dense_vector_cdp','IN',this_module,'dense_axpby_cdp')
        end select
    end subroutine dense_axpby_cdp

    function dense_dot_cdp(self, vec) result(alpha)
        implicit none(type, external)
        class(dense_vector_cdp), intent(in) :: self
        class(abstract_vector_cdp), intent(in) :: vec
        complex(dp) :: alpha
        integer :: n
        n = self%get_size()
        select type (vec)
        type is(dense_vector_cdp)
            alpha = dotc(n, self%data, 1, vec%data, 1)
        class default
            call type_error('vec','dense_vector_cdp','IN',this_module,'dense_dot_cdp')
        end select
    end function dense_dot_cdp

    function dense_get_size_cdp(self) result(n)
        implicit none(type, external)
        class(dense_vector_cdp), intent(in) :: self
        integer :: n
        n = size(self%data)
    end function dense_get_size_cdp

    
    !--------------------------------------
    !-----      UTILITY FUNCTIONS     -----
    !--------------------------------------

    subroutine linear_combination_vector_rsp(y, X, v)
        !! Given `X` and `v`, this function return \( \mathbf{y} = \mathbf{Xv} \) where
        !! `y` is an `abstract_vector`, `X` an array of `abstract_vector` and `v` a
        !! Fortran array containing the coefficients of the linear combination.
        implicit none(type, external)
        class(abstract_vector_rsp), allocatable, intent(out) :: y
        !! Ouput vector.
        class(abstract_vector_rsp), intent(in) :: X(:)
        !! Krylov basis.
        real(sp), intent(in) :: v(:)
        !! Coordinates of `y` in the Krylov basis `X`.

        ! Internal variables
        integer :: i, iostat
        character(len=100) :: errmsg

        ! Check sizes.
        if (size(X) /= size(v)) then
            call stop_error("Krylov basis X and low-dimensional vector v have different sizes.", &
                              & this_module, 'linear_combination_vector_rsp')
        endif

        ! Initialize output vector.
        if (.not. allocated(y)) then
            allocate(y, source=X(1), stat=iostat, errmsg=errmsg)
            if (iostat /= 0) call stop_error(errmsg, &
                                             module=this_module, &
                                             procedure="linear_combination_vector_rsp")
        endif
        call y%zero()
        ! Compute linear combination.
        do i = 1, size(X)
            call y%axpby(v(i), X(i), one_rsp) ! y = y + X[i]*v[i]
        enddo
    end subroutine linear_combination_vector_rsp

    subroutine linear_combination_matrix_rsp(Y, X, B)
        !! Given `X` and `B`, this function computes \(\mathbf{Y} = \mathbf{XB} \) where
        !! `X` and `Y` are arrays of `abstract_vector`, and `B` is a 2D Fortran array.
        implicit none(type, external)
        class(abstract_vector_rsp), allocatable, intent(out) :: Y(:)
        !! Output matrix.
        class(abstract_vector_rsp), intent(in) :: X(:)
        !! Krylov basis.
        real(sp), intent(in) :: B(:, :)
        !! Coefficients of the linear combinations.

        ! Internal variables.
        integer :: i, j, iostat
        character(len=100) :: errmsg
    
        ! Check sizes.
        if (size(X) /= size(B, 1)) then
            call stop_error("Krylov basis X and combination matrix B have incompatible sizes.", &
                            this_module, 'linear_combination_matrix_rsp')
        endif

        ! Initialize output basis.
        if (.not. allocated(Y)) then
            allocate(Y(size(B, 2)), source=X(1), stat=iostat, errmsg=errmsg)
            if (iostat /= 0) call stop_error(errmsg, &
                                             module=this_module, &
                                             procedure="linear_combination_matrix_rsp")
        else
            if (size(Y) /= size(B, 2)) then
                call stop_error("Krylov basis Y and combination matrix B have incompatible sizes.", &
                                this_module, 'linear_combination_matrix_rsp')
            endif
        endif

        do j = 1, size(Y)
            call Y(j)%zero()
            do i = 1, size(X)
                call Y(j)%axpby(B(i, j), X(i), one_rsp) ! y(j) = B(i,j)*X(i) + y(j)
            enddo
        enddo
    end subroutine linear_combination_matrix_rsp

    function gram_matrix_rsp(X) result(G)
        !! Computes the inner product/Gram matrix associated with the basis \( \mathbf{X} \).
        implicit none(type, external)
        class(abstract_vector_rsp), intent(in) :: X(:)
        real(sp) :: G(size(X), size(X))
        integer :: i, j
        do i = 1, size(X)
            do j = i, size(X)
                G(i, j) = X(i)%dot(X(j))
                G(j, i) = G(i, j)
            enddo
        enddo
    end function gram_matrix_rsp

    function innerprod_vector_rsp(X, v) result(y)
        !! Computes the inner product vector \( \mathbf{y} = \mathbf{X}^H \mathbf{v} \) between
        !! a basis `X` of `abstract_vector` and `v`, a single `abstract_vector`.
        implicit none(type, external)
        class(abstract_vector_rsp), intent(in) :: X(:), v
        !! Basis and single instance of `abstract_vector` whose inner products need to be computed.
        real(sp) :: y(size(X))
        !! Resulting inner-product vector.

        ! Local variables.
        integer :: i

        y = zero_rsp
        do i = 1, size(X)
            y(i) = X(i)%dot(v)
        enddo
    end function innerprod_vector_rsp

    function innerprod_matrix_rsp(X, Y) result(M)
        !! Computes the inner product matrix \( \mathbf{M} = \mathbf{X}^H \mathbf{Y} \) between
        !! two bases of `abstract_vector`.
        implicit none(type, external)
        class(abstract_vector_rsp), intent(in) :: X(:), Y(:)
        !! Bases of `abstract_vector` whose inner products need to be computed.
        real(sp) :: M(size(X), size(Y))
        !! Resulting inner-product matrix.

        ! Local variables.
        integer :: i, j

        M = zero_rsp
        do j = 1, size(Y)
            do i = 1, size(X)
                M(i, j) = X(i)%dot(Y(j))
            enddo
        enddo
    end function innerprod_matrix_rsp

    impure elemental subroutine axpby_basis_rsp(alpha, X, beta, Y)
        !! Compute in-place \( \mathbf{Y} \leftarrow \alpha \mathbf{X} + \beta \mathbf{Y} \) where
        !! `X` and `Y` are arrays of `abstract_vector` and `alpha` and `beta` are real(sp)
        !! numbers.
        implicit none(type, external)
        class(abstract_vector_rsp), intent(in) :: X
        !! Input/Ouput array of `abstract_vector`.
        class(abstract_vector_rsp), intent(inout) :: Y
        !! Array of `abstract_vector` to be added/subtracted to/from `X`.
        real(sp), intent(in) :: alpha, beta
        !! Scalar multipliers.
        call Y%axpby(alpha, X, beta)
    end subroutine axpby_basis_rsp

    impure elemental subroutine zero_basis_rsp(X)
        implicit none(type, external)
        class(abstract_vector_rsp), intent(inout) :: X
        call X%zero()
    end subroutine zero_basis_rsp

    impure elemental subroutine copy_vector_rsp(out, from)
        implicit none(type, external)
        class(abstract_vector_rsp), intent(in) :: from
        class(abstract_vector_rsp), intent(out) :: out
        ! Copy array.
        call out%axpby(one_rsp, from, zero_rsp)
    end subroutine copy_vector_rsp

    impure elemental subroutine rand_basis_rsp(X, ifnorm)
        implicit none(type, external)
        class(abstract_vector_rsp), intent(inout) :: X
        logical, optional, intent(in) :: ifnorm
        call X%rand(ifnorm=ifnorm)
    end subroutine rand_basis_rsp

    subroutine linear_combination_vector_rdp(y, X, v)
        !! Given `X` and `v`, this function return \( \mathbf{y} = \mathbf{Xv} \) where
        !! `y` is an `abstract_vector`, `X` an array of `abstract_vector` and `v` a
        !! Fortran array containing the coefficients of the linear combination.
        implicit none(type, external)
        class(abstract_vector_rdp), allocatable, intent(out) :: y
        !! Ouput vector.
        class(abstract_vector_rdp), intent(in) :: X(:)
        !! Krylov basis.
        real(dp), intent(in) :: v(:)
        !! Coordinates of `y` in the Krylov basis `X`.

        ! Internal variables
        integer :: i, iostat
        character(len=100) :: errmsg

        ! Check sizes.
        if (size(X) /= size(v)) then
            call stop_error("Krylov basis X and low-dimensional vector v have different sizes.", &
                              & this_module, 'linear_combination_vector_rdp')
        endif

        ! Initialize output vector.
        if (.not. allocated(y)) then
            allocate(y, source=X(1), stat=iostat, errmsg=errmsg)
            if (iostat /= 0) call stop_error(errmsg, &
                                             module=this_module, &
                                             procedure="linear_combination_vector_rdp")
        endif
        call y%zero()
        ! Compute linear combination.
        do i = 1, size(X)
            call y%axpby(v(i), X(i), one_rdp) ! y = y + X[i]*v[i]
        enddo
    end subroutine linear_combination_vector_rdp

    subroutine linear_combination_matrix_rdp(Y, X, B)
        !! Given `X` and `B`, this function computes \(\mathbf{Y} = \mathbf{XB} \) where
        !! `X` and `Y` are arrays of `abstract_vector`, and `B` is a 2D Fortran array.
        implicit none(type, external)
        class(abstract_vector_rdp), allocatable, intent(out) :: Y(:)
        !! Output matrix.
        class(abstract_vector_rdp), intent(in) :: X(:)
        !! Krylov basis.
        real(dp), intent(in) :: B(:, :)
        !! Coefficients of the linear combinations.

        ! Internal variables.
        integer :: i, j, iostat
        character(len=100) :: errmsg
    
        ! Check sizes.
        if (size(X) /= size(B, 1)) then
            call stop_error("Krylov basis X and combination matrix B have incompatible sizes.", &
                            this_module, 'linear_combination_matrix_rdp')
        endif

        ! Initialize output basis.
        if (.not. allocated(Y)) then
            allocate(Y(size(B, 2)), source=X(1), stat=iostat, errmsg=errmsg)
            if (iostat /= 0) call stop_error(errmsg, &
                                             module=this_module, &
                                             procedure="linear_combination_matrix_rdp")
        else
            if (size(Y) /= size(B, 2)) then
                call stop_error("Krylov basis Y and combination matrix B have incompatible sizes.", &
                                this_module, 'linear_combination_matrix_rdp')
            endif
        endif

        do j = 1, size(Y)
            call Y(j)%zero()
            do i = 1, size(X)
                call Y(j)%axpby(B(i, j), X(i), one_rdp) ! y(j) = B(i,j)*X(i) + y(j)
            enddo
        enddo
    end subroutine linear_combination_matrix_rdp

    function gram_matrix_rdp(X) result(G)
        !! Computes the inner product/Gram matrix associated with the basis \( \mathbf{X} \).
        implicit none(type, external)
        class(abstract_vector_rdp), intent(in) :: X(:)
        real(dp) :: G(size(X), size(X))
        integer :: i, j
        do i = 1, size(X)
            do j = i, size(X)
                G(i, j) = X(i)%dot(X(j))
                G(j, i) = G(i, j)
            enddo
        enddo
    end function gram_matrix_rdp

    function innerprod_vector_rdp(X, v) result(y)
        !! Computes the inner product vector \( \mathbf{y} = \mathbf{X}^H \mathbf{v} \) between
        !! a basis `X` of `abstract_vector` and `v`, a single `abstract_vector`.
        implicit none(type, external)
        class(abstract_vector_rdp), intent(in) :: X(:), v
        !! Basis and single instance of `abstract_vector` whose inner products need to be computed.
        real(dp) :: y(size(X))
        !! Resulting inner-product vector.

        ! Local variables.
        integer :: i

        y = zero_rdp
        do i = 1, size(X)
            y(i) = X(i)%dot(v)
        enddo
    end function innerprod_vector_rdp

    function innerprod_matrix_rdp(X, Y) result(M)
        !! Computes the inner product matrix \( \mathbf{M} = \mathbf{X}^H \mathbf{Y} \) between
        !! two bases of `abstract_vector`.
        implicit none(type, external)
        class(abstract_vector_rdp), intent(in) :: X(:), Y(:)
        !! Bases of `abstract_vector` whose inner products need to be computed.
        real(dp) :: M(size(X), size(Y))
        !! Resulting inner-product matrix.

        ! Local variables.
        integer :: i, j

        M = zero_rdp
        do j = 1, size(Y)
            do i = 1, size(X)
                M(i, j) = X(i)%dot(Y(j))
            enddo
        enddo
    end function innerprod_matrix_rdp

    impure elemental subroutine axpby_basis_rdp(alpha, X, beta, Y)
        !! Compute in-place \( \mathbf{Y} \leftarrow \alpha \mathbf{X} + \beta \mathbf{Y} \) where
        !! `X` and `Y` are arrays of `abstract_vector` and `alpha` and `beta` are real(dp)
        !! numbers.
        implicit none(type, external)
        class(abstract_vector_rdp), intent(in) :: X
        !! Input/Ouput array of `abstract_vector`.
        class(abstract_vector_rdp), intent(inout) :: Y
        !! Array of `abstract_vector` to be added/subtracted to/from `X`.
        real(dp), intent(in) :: alpha, beta
        !! Scalar multipliers.
        call Y%axpby(alpha, X, beta)
    end subroutine axpby_basis_rdp

    impure elemental subroutine zero_basis_rdp(X)
        implicit none(type, external)
        class(abstract_vector_rdp), intent(inout) :: X
        call X%zero()
    end subroutine zero_basis_rdp

    impure elemental subroutine copy_vector_rdp(out, from)
        implicit none(type, external)
        class(abstract_vector_rdp), intent(in) :: from
        class(abstract_vector_rdp), intent(out) :: out
        ! Copy array.
        call out%axpby(one_rdp, from, zero_rdp)
    end subroutine copy_vector_rdp

    impure elemental subroutine rand_basis_rdp(X, ifnorm)
        implicit none(type, external)
        class(abstract_vector_rdp), intent(inout) :: X
        logical, optional, intent(in) :: ifnorm
        call X%rand(ifnorm=ifnorm)
    end subroutine rand_basis_rdp

    subroutine linear_combination_vector_csp(y, X, v)
        !! Given `X` and `v`, this function return \( \mathbf{y} = \mathbf{Xv} \) where
        !! `y` is an `abstract_vector`, `X` an array of `abstract_vector` and `v` a
        !! Fortran array containing the coefficients of the linear combination.
        implicit none(type, external)
        class(abstract_vector_csp), allocatable, intent(out) :: y
        !! Ouput vector.
        class(abstract_vector_csp), intent(in) :: X(:)
        !! Krylov basis.
        complex(sp), intent(in) :: v(:)
        !! Coordinates of `y` in the Krylov basis `X`.

        ! Internal variables
        integer :: i, iostat
        character(len=100) :: errmsg

        ! Check sizes.
        if (size(X) /= size(v)) then
            call stop_error("Krylov basis X and low-dimensional vector v have different sizes.", &
                              & this_module, 'linear_combination_vector_csp')
        endif

        ! Initialize output vector.
        if (.not. allocated(y)) then
            allocate(y, source=X(1), stat=iostat, errmsg=errmsg)
            if (iostat /= 0) call stop_error(errmsg, &
                                             module=this_module, &
                                             procedure="linear_combination_vector_csp")
        endif
        call y%zero()
        ! Compute linear combination.
        do i = 1, size(X)
            call y%axpby(v(i), X(i), one_csp) ! y = y + X[i]*v[i]
        enddo
    end subroutine linear_combination_vector_csp

    subroutine linear_combination_matrix_csp(Y, X, B)
        !! Given `X` and `B`, this function computes \(\mathbf{Y} = \mathbf{XB} \) where
        !! `X` and `Y` are arrays of `abstract_vector`, and `B` is a 2D Fortran array.
        implicit none(type, external)
        class(abstract_vector_csp), allocatable, intent(out) :: Y(:)
        !! Output matrix.
        class(abstract_vector_csp), intent(in) :: X(:)
        !! Krylov basis.
        complex(sp), intent(in) :: B(:, :)
        !! Coefficients of the linear combinations.

        ! Internal variables.
        integer :: i, j, iostat
        character(len=100) :: errmsg
    
        ! Check sizes.
        if (size(X) /= size(B, 1)) then
            call stop_error("Krylov basis X and combination matrix B have incompatible sizes.", &
                            this_module, 'linear_combination_matrix_csp')
        endif

        ! Initialize output basis.
        if (.not. allocated(Y)) then
            allocate(Y(size(B, 2)), source=X(1), stat=iostat, errmsg=errmsg)
            if (iostat /= 0) call stop_error(errmsg, &
                                             module=this_module, &
                                             procedure="linear_combination_matrix_csp")
        else
            if (size(Y) /= size(B, 2)) then
                call stop_error("Krylov basis Y and combination matrix B have incompatible sizes.", &
                                this_module, 'linear_combination_matrix_csp')
            endif
        endif

        do j = 1, size(Y)
            call Y(j)%zero()
            do i = 1, size(X)
                call Y(j)%axpby(B(i, j), X(i), one_csp) ! y(j) = B(i,j)*X(i) + y(j)
            enddo
        enddo
    end subroutine linear_combination_matrix_csp

    function gram_matrix_csp(X) result(G)
        !! Computes the inner product/Gram matrix associated with the basis \( \mathbf{X} \).
        implicit none(type, external)
        class(abstract_vector_csp), intent(in) :: X(:)
        complex(sp) :: G(size(X), size(X))
        integer :: i, j
        do i = 1, size(X)
            do j = i, size(X)
                G(i, j) = X(i)%dot(X(j))
                G(j, i) = G(i, j)
            enddo
        enddo
    end function gram_matrix_csp

    function innerprod_vector_csp(X, v) result(y)
        !! Computes the inner product vector \( \mathbf{y} = \mathbf{X}^H \mathbf{v} \) between
        !! a basis `X` of `abstract_vector` and `v`, a single `abstract_vector`.
        implicit none(type, external)
        class(abstract_vector_csp), intent(in) :: X(:), v
        !! Basis and single instance of `abstract_vector` whose inner products need to be computed.
        complex(sp) :: y(size(X))
        !! Resulting inner-product vector.

        ! Local variables.
        integer :: i

        y = zero_csp
        do i = 1, size(X)
            y(i) = X(i)%dot(v)
        enddo
    end function innerprod_vector_csp

    function innerprod_matrix_csp(X, Y) result(M)
        !! Computes the inner product matrix \( \mathbf{M} = \mathbf{X}^H \mathbf{Y} \) between
        !! two bases of `abstract_vector`.
        implicit none(type, external)
        class(abstract_vector_csp), intent(in) :: X(:), Y(:)
        !! Bases of `abstract_vector` whose inner products need to be computed.
        complex(sp) :: M(size(X), size(Y))
        !! Resulting inner-product matrix.

        ! Local variables.
        integer :: i, j

        M = zero_csp
        do j = 1, size(Y)
            do i = 1, size(X)
                M(i, j) = X(i)%dot(Y(j))
            enddo
        enddo
    end function innerprod_matrix_csp

    impure elemental subroutine axpby_basis_csp(alpha, X, beta, Y)
        !! Compute in-place \( \mathbf{Y} \leftarrow \alpha \mathbf{X} + \beta \mathbf{Y} \) where
        !! `X` and `Y` are arrays of `abstract_vector` and `alpha` and `beta` are complex(sp)
        !! numbers.
        implicit none(type, external)
        class(abstract_vector_csp), intent(in) :: X
        !! Input/Ouput array of `abstract_vector`.
        class(abstract_vector_csp), intent(inout) :: Y
        !! Array of `abstract_vector` to be added/subtracted to/from `X`.
        complex(sp), intent(in) :: alpha, beta
        !! Scalar multipliers.
        call Y%axpby(alpha, X, beta)
    end subroutine axpby_basis_csp

    impure elemental subroutine zero_basis_csp(X)
        implicit none(type, external)
        class(abstract_vector_csp), intent(inout) :: X
        call X%zero()
    end subroutine zero_basis_csp

    impure elemental subroutine copy_vector_csp(out, from)
        implicit none(type, external)
        class(abstract_vector_csp), intent(in) :: from
        class(abstract_vector_csp), intent(out) :: out
        ! Copy array.
        call out%axpby(one_csp, from, zero_csp)
    end subroutine copy_vector_csp

    impure elemental subroutine rand_basis_csp(X, ifnorm)
        implicit none(type, external)
        class(abstract_vector_csp), intent(inout) :: X
        logical, optional, intent(in) :: ifnorm
        call X%rand(ifnorm=ifnorm)
    end subroutine rand_basis_csp

    subroutine linear_combination_vector_cdp(y, X, v)
        !! Given `X` and `v`, this function return \( \mathbf{y} = \mathbf{Xv} \) where
        !! `y` is an `abstract_vector`, `X` an array of `abstract_vector` and `v` a
        !! Fortran array containing the coefficients of the linear combination.
        implicit none(type, external)
        class(abstract_vector_cdp), allocatable, intent(out) :: y
        !! Ouput vector.
        class(abstract_vector_cdp), intent(in) :: X(:)
        !! Krylov basis.
        complex(dp), intent(in) :: v(:)
        !! Coordinates of `y` in the Krylov basis `X`.

        ! Internal variables
        integer :: i, iostat
        character(len=100) :: errmsg

        ! Check sizes.
        if (size(X) /= size(v)) then
            call stop_error("Krylov basis X and low-dimensional vector v have different sizes.", &
                              & this_module, 'linear_combination_vector_cdp')
        endif

        ! Initialize output vector.
        if (.not. allocated(y)) then
            allocate(y, source=X(1), stat=iostat, errmsg=errmsg)
            if (iostat /= 0) call stop_error(errmsg, &
                                             module=this_module, &
                                             procedure="linear_combination_vector_cdp")
        endif
        call y%zero()
        ! Compute linear combination.
        do i = 1, size(X)
            call y%axpby(v(i), X(i), one_cdp) ! y = y + X[i]*v[i]
        enddo
    end subroutine linear_combination_vector_cdp

    subroutine linear_combination_matrix_cdp(Y, X, B)
        !! Given `X` and `B`, this function computes \(\mathbf{Y} = \mathbf{XB} \) where
        !! `X` and `Y` are arrays of `abstract_vector`, and `B` is a 2D Fortran array.
        implicit none(type, external)
        class(abstract_vector_cdp), allocatable, intent(out) :: Y(:)
        !! Output matrix.
        class(abstract_vector_cdp), intent(in) :: X(:)
        !! Krylov basis.
        complex(dp), intent(in) :: B(:, :)
        !! Coefficients of the linear combinations.

        ! Internal variables.
        integer :: i, j, iostat
        character(len=100) :: errmsg
    
        ! Check sizes.
        if (size(X) /= size(B, 1)) then
            call stop_error("Krylov basis X and combination matrix B have incompatible sizes.", &
                            this_module, 'linear_combination_matrix_cdp')
        endif

        ! Initialize output basis.
        if (.not. allocated(Y)) then
            allocate(Y(size(B, 2)), source=X(1), stat=iostat, errmsg=errmsg)
            if (iostat /= 0) call stop_error(errmsg, &
                                             module=this_module, &
                                             procedure="linear_combination_matrix_cdp")
        else
            if (size(Y) /= size(B, 2)) then
                call stop_error("Krylov basis Y and combination matrix B have incompatible sizes.", &
                                this_module, 'linear_combination_matrix_cdp')
            endif
        endif

        do j = 1, size(Y)
            call Y(j)%zero()
            do i = 1, size(X)
                call Y(j)%axpby(B(i, j), X(i), one_cdp) ! y(j) = B(i,j)*X(i) + y(j)
            enddo
        enddo
    end subroutine linear_combination_matrix_cdp

    function gram_matrix_cdp(X) result(G)
        !! Computes the inner product/Gram matrix associated with the basis \( \mathbf{X} \).
        implicit none(type, external)
        class(abstract_vector_cdp), intent(in) :: X(:)
        complex(dp) :: G(size(X), size(X))
        integer :: i, j
        do i = 1, size(X)
            do j = i, size(X)
                G(i, j) = X(i)%dot(X(j))
                G(j, i) = G(i, j)
            enddo
        enddo
    end function gram_matrix_cdp

    function innerprod_vector_cdp(X, v) result(y)
        !! Computes the inner product vector \( \mathbf{y} = \mathbf{X}^H \mathbf{v} \) between
        !! a basis `X` of `abstract_vector` and `v`, a single `abstract_vector`.
        implicit none(type, external)
        class(abstract_vector_cdp), intent(in) :: X(:), v
        !! Basis and single instance of `abstract_vector` whose inner products need to be computed.
        complex(dp) :: y(size(X))
        !! Resulting inner-product vector.

        ! Local variables.
        integer :: i

        y = zero_cdp
        do i = 1, size(X)
            y(i) = X(i)%dot(v)
        enddo
    end function innerprod_vector_cdp

    function innerprod_matrix_cdp(X, Y) result(M)
        !! Computes the inner product matrix \( \mathbf{M} = \mathbf{X}^H \mathbf{Y} \) between
        !! two bases of `abstract_vector`.
        implicit none(type, external)
        class(abstract_vector_cdp), intent(in) :: X(:), Y(:)
        !! Bases of `abstract_vector` whose inner products need to be computed.
        complex(dp) :: M(size(X), size(Y))
        !! Resulting inner-product matrix.

        ! Local variables.
        integer :: i, j

        M = zero_cdp
        do j = 1, size(Y)
            do i = 1, size(X)
                M(i, j) = X(i)%dot(Y(j))
            enddo
        enddo
    end function innerprod_matrix_cdp

    impure elemental subroutine axpby_basis_cdp(alpha, X, beta, Y)
        !! Compute in-place \( \mathbf{Y} \leftarrow \alpha \mathbf{X} + \beta \mathbf{Y} \) where
        !! `X` and `Y` are arrays of `abstract_vector` and `alpha` and `beta` are complex(dp)
        !! numbers.
        implicit none(type, external)
        class(abstract_vector_cdp), intent(in) :: X
        !! Input/Ouput array of `abstract_vector`.
        class(abstract_vector_cdp), intent(inout) :: Y
        !! Array of `abstract_vector` to be added/subtracted to/from `X`.
        complex(dp), intent(in) :: alpha, beta
        !! Scalar multipliers.
        call Y%axpby(alpha, X, beta)
    end subroutine axpby_basis_cdp

    impure elemental subroutine zero_basis_cdp(X)
        implicit none(type, external)
        class(abstract_vector_cdp), intent(inout) :: X
        call X%zero()
    end subroutine zero_basis_cdp

    impure elemental subroutine copy_vector_cdp(out, from)
        implicit none(type, external)
        class(abstract_vector_cdp), intent(in) :: from
        class(abstract_vector_cdp), intent(out) :: out
        ! Copy array.
        call out%axpby(one_cdp, from, zero_cdp)
    end subroutine copy_vector_cdp

    impure elemental subroutine rand_basis_cdp(X, ifnorm)
        implicit none(type, external)
        class(abstract_vector_cdp), intent(inout) :: X
        logical, optional, intent(in) :: ifnorm
        call X%rand(ifnorm=ifnorm)
    end subroutine rand_basis_cdp

end module LightKrylov_AbstractVectors
