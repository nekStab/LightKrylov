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
    !! - `init_like(mold)`              :   A subroutine to initialise the vector after allocation to match the size and structure of the mold.
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
    !! - `copy_basis(out, from)`            : Copy a collection of `abstract_vectors`. `out` must be pre-allocated.
    !! - `init_like_basis(out, mold)`       : Initialize size and structure of a collection of `abstract_vectors` based on mold.
    !! - `free_basis()         `            : Free the memory of a collection of `abstract_vectors`.
    !! - `rand_basis(X, ifnorm)`            : Create a collection of random `abstract_vectors`. If `ifnorm = .true.`, the vectors are normalized to have unit-norm.
    !! @warning
    !! The resulting vectors do not form an orthonormal basis. For this purpose use the utility function `initialize_random_orthonormal_basis`.
    !! @endwarning

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
    public :: init_like_basis
    public :: free_basis
    public :: verify_vector_axioms

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
        !!  Internally, copy will run Y(i)%init_like(X(1)) to ensure that Y is conformant with X.
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

    interface init_like_basis
        !!  This interface provides methods to initialize an array `X` of `abstract_vector` 
        !!  based on the shape of another array `mold` of `abstract_vector`.
        !!  It is a simple wrapper around `X(i)%init_like(y)`.
        !!
        !!  ### Example
        !!
        !!  ```fortran
        !!      type(my_real_vector), dimension(10), allocatable :: X
        !!      type(my_real_vector) :: mold
        !!
        !!      ! ... Your code ...
        !!
        !!      allocate(X, mold) ; call init_like_basis(X, mold)
        !!
        !!      ! ... Your code ...
        !!  ```
        module procedure init_like_basis_rsp
        module procedure init_like_basis_rdp
        module procedure init_like_basis_csp
        module procedure init_like_basis_cdp
    end interface

    interface free_basis
        !!  This interface provides methods to free the memory of an array `X` of `abstract_vector`.
        !!  It is a simple wrapper around `X(i)%free()`.
        !!
        !!  ### Example
        !!
        !!  ```fortran
        !!      type(my_real_vector), dimension(10), allocatable :: X
        !!      type(my_real_vector) :: mold
        !!
        !!      allocate(X, mold) ; call init_like_basis(X, mold)
        !!
        !!      ! ... Your code ...
        !!
        !!      call free_basis(X)
        !!
        !!      ! ... Your code ...
        !!  ```
        module procedure free_basis_rsp
        module procedure free_basis_rdp
        module procedure free_basis_csp
        module procedure free_basis_cdp
    end interface

    interface verify_vector_axioms
        module procedure verify_vector_axioms_rsp
        module procedure verify_vector_axioms_rdp
        module procedure verify_vector_axioms_csp
        module procedure verify_vector_axioms_cdp
    end interface


    type, abstract, public :: abstract_vector
        !!  Base abstract type from which all other types of vectors used in `LightKrylov`
        !!  are being derived from.
        !!
        !!  @warning
        !!  Users should not extend this abstract class to define their own types.
        !!  @endwarning
        logical :: is_initialized = .false.
        !! Set by `init_like`, cleared by `free`. Bookkeeping only for the managed CPU types
        !! here; used by unmanaged-resource (e.g. GPU) extensions to make `free` idempotent.
        logical :: owns_data      = .true.
        !! `.false.` marks an alias/view so `free` won't release shared storage. For GPU/
        !! resource-managing extensions; inert and untested for the CPU types shipped here.
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
        procedure(abstract_init_like_rsp), pass(self), deferred, public :: init_like
        !! Shapes `self` like `mold`: allocates if unallocated, no-op if already conformant,
        !! releases and reallocates if the shape differs. Sole guaranteed allocator for
        !! type-agnostic code; for resource-managing (e.g. GPU) types, this is where device
        !! buffers are allocated/reshaped/released.
        procedure, pass(self), public :: free => free_rsp
        !! Releases any unmanaged resources owned by `self`. Terminal: after `free`, the only
        !! valid operation on `self` is another `free`. Idempotent: safe to call on an already-
        !! freed vector. The default is a no-op, correct for types whose storage is a managed
        !! Fortran allocatable (cleaned up automatically). Types holding unmanaged resources
        !! (e.g. device buffers) override this to release them, behind a guard so the explicit
        !! end-of-use call and any `final` backstop never double-release.
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

        subroutine abstract_init_like_rsp(self, mold)
            !! Abstract interface to shape `self` like `mold`.
            !!
            !! Contract — the implementation must leave `self` shaped like `mold`, handling
            !! three cases:
            !!
            !! - `self` unallocated            : allocate internal storage from `mold`'s shape.
            !! - `self` allocated, conformant  : no-op (cheap fast path; do not reallocate).
            !! - `self` allocated, nonconformant: release prior storage, then reallocate.
            !!
            !! Must be idempotent. This is the sole guaranteed allocator used by type-agnostic
            !! algorithms, which call it immediately after `allocate(self, mold=...)`. For types
            !! holding unmanaged resources (e.g. device buffers), this routine owns their
            !! allocation, reshaping, and release.
            import abstract_vector_rsp
            implicit none(type, external)
            class(abstract_vector_rsp), intent(inout) :: self
            !! Vector to be shaped.
            class(abstract_vector_rsp), intent(in) :: mold
            !! Vector whose shape `self` should match.
        end subroutine abstract_init_like_rsp

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
        procedure(abstract_init_like_rdp), pass(self), deferred, public :: init_like
        !! Shapes `self` like `mold`: allocates if unallocated, no-op if already conformant,
        !! releases and reallocates if the shape differs. Sole guaranteed allocator for
        !! type-agnostic code; for resource-managing (e.g. GPU) types, this is where device
        !! buffers are allocated/reshaped/released.
        procedure, pass(self), public :: free => free_rdp
        !! Releases any unmanaged resources owned by `self`. Terminal: after `free`, the only
        !! valid operation on `self` is another `free`. Idempotent: safe to call on an already-
        !! freed vector. The default is a no-op, correct for types whose storage is a managed
        !! Fortran allocatable (cleaned up automatically). Types holding unmanaged resources
        !! (e.g. device buffers) override this to release them, behind a guard so the explicit
        !! end-of-use call and any `final` backstop never double-release.
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

        subroutine abstract_init_like_rdp(self, mold)
            !! Abstract interface to shape `self` like `mold`.
            !!
            !! Contract — the implementation must leave `self` shaped like `mold`, handling
            !! three cases:
            !!
            !! - `self` unallocated            : allocate internal storage from `mold`'s shape.
            !! - `self` allocated, conformant  : no-op (cheap fast path; do not reallocate).
            !! - `self` allocated, nonconformant: release prior storage, then reallocate.
            !!
            !! Must be idempotent. This is the sole guaranteed allocator used by type-agnostic
            !! algorithms, which call it immediately after `allocate(self, mold=...)`. For types
            !! holding unmanaged resources (e.g. device buffers), this routine owns their
            !! allocation, reshaping, and release.
            import abstract_vector_rdp
            implicit none(type, external)
            class(abstract_vector_rdp), intent(inout) :: self
            !! Vector to be shaped.
            class(abstract_vector_rdp), intent(in) :: mold
            !! Vector whose shape `self` should match.
        end subroutine abstract_init_like_rdp

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
        procedure(abstract_init_like_csp), pass(self), deferred, public :: init_like
        !! Shapes `self` like `mold`: allocates if unallocated, no-op if already conformant,
        !! releases and reallocates if the shape differs. Sole guaranteed allocator for
        !! type-agnostic code; for resource-managing (e.g. GPU) types, this is where device
        !! buffers are allocated/reshaped/released.
        procedure, pass(self), public :: free => free_csp
        !! Releases any unmanaged resources owned by `self`. Terminal: after `free`, the only
        !! valid operation on `self` is another `free`. Idempotent: safe to call on an already-
        !! freed vector. The default is a no-op, correct for types whose storage is a managed
        !! Fortran allocatable (cleaned up automatically). Types holding unmanaged resources
        !! (e.g. device buffers) override this to release them, behind a guard so the explicit
        !! end-of-use call and any `final` backstop never double-release.
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

        subroutine abstract_init_like_csp(self, mold)
            !! Abstract interface to shape `self` like `mold`.
            !!
            !! Contract — the implementation must leave `self` shaped like `mold`, handling
            !! three cases:
            !!
            !! - `self` unallocated            : allocate internal storage from `mold`'s shape.
            !! - `self` allocated, conformant  : no-op (cheap fast path; do not reallocate).
            !! - `self` allocated, nonconformant: release prior storage, then reallocate.
            !!
            !! Must be idempotent. This is the sole guaranteed allocator used by type-agnostic
            !! algorithms, which call it immediately after `allocate(self, mold=...)`. For types
            !! holding unmanaged resources (e.g. device buffers), this routine owns their
            !! allocation, reshaping, and release.
            import abstract_vector_csp
            implicit none(type, external)
            class(abstract_vector_csp), intent(inout) :: self
            !! Vector to be shaped.
            class(abstract_vector_csp), intent(in) :: mold
            !! Vector whose shape `self` should match.
        end subroutine abstract_init_like_csp

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
        procedure(abstract_init_like_cdp), pass(self), deferred, public :: init_like
        !! Shapes `self` like `mold`: allocates if unallocated, no-op if already conformant,
        !! releases and reallocates if the shape differs. Sole guaranteed allocator for
        !! type-agnostic code; for resource-managing (e.g. GPU) types, this is where device
        !! buffers are allocated/reshaped/released.
        procedure, pass(self), public :: free => free_cdp
        !! Releases any unmanaged resources owned by `self`. Terminal: after `free`, the only
        !! valid operation on `self` is another `free`. Idempotent: safe to call on an already-
        !! freed vector. The default is a no-op, correct for types whose storage is a managed
        !! Fortran allocatable (cleaned up automatically). Types holding unmanaged resources
        !! (e.g. device buffers) override this to release them, behind a guard so the explicit
        !! end-of-use call and any `final` backstop never double-release.
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

        subroutine abstract_init_like_cdp(self, mold)
            !! Abstract interface to shape `self` like `mold`.
            !!
            !! Contract — the implementation must leave `self` shaped like `mold`, handling
            !! three cases:
            !!
            !! - `self` unallocated            : allocate internal storage from `mold`'s shape.
            !! - `self` allocated, conformant  : no-op (cheap fast path; do not reallocate).
            !! - `self` allocated, nonconformant: release prior storage, then reallocate.
            !!
            !! Must be idempotent. This is the sole guaranteed allocator used by type-agnostic
            !! algorithms, which call it immediately after `allocate(self, mold=...)`. For types
            !! holding unmanaged resources (e.g. device buffers), this routine owns their
            !! allocation, reshaping, and release.
            import abstract_vector_cdp
            implicit none(type, external)
            class(abstract_vector_cdp), intent(inout) :: self
            !! Vector to be shaped.
            class(abstract_vector_cdp), intent(in) :: mold
            !! Vector whose shape `self` should match.
        end subroutine abstract_init_like_cdp

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
        procedure, pass(self), public :: init_like => dense_init_like_rsp
        !! Initialize `self` like `mold`: allocate if unallocated, no-op if already conformant, release and reallocate if the shape differs.
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
        procedure, pass(self), public :: init_like => dense_init_like_rdp
        !! Initialize `self` like `mold`: allocate if unallocated, no-op if already conformant, release and reallocate if the shape differs.
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
        procedure, pass(self), public :: init_like => dense_init_like_csp
        !! Initialize `self` like `mold`: allocate if unallocated, no-op if already conformant, release and reallocate if the shape differs.
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
        procedure, pass(self), public :: init_like => dense_init_like_cdp
        !! Initialize `self` like `mold`: allocate if unallocated, no-op if already conformant, release and reallocate if the shape differs.
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

    subroutine free_rsp(self)
        implicit none(type, external)
        !! Default `free`: no-op. Correct for types whose only storage is a managed Fortran
        !! allocatable, which is released automatically on deallocation / scope exit. Types
        !! holding unmanaged resources override this to release them, guarded for idempotency.
        !! Clears the `is_initialized` flag so that the action is terminal.
        class(abstract_vector_rsp), intent(inout) :: self
        self%is_initialized = .false.
    end subroutine free_rsp

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

    subroutine free_rdp(self)
        implicit none(type, external)
        !! Default `free`: no-op. Correct for types whose only storage is a managed Fortran
        !! allocatable, which is released automatically on deallocation / scope exit. Types
        !! holding unmanaged resources override this to release them, guarded for idempotency.
        !! Clears the `is_initialized` flag so that the action is terminal.
        class(abstract_vector_rdp), intent(inout) :: self
        self%is_initialized = .false.
    end subroutine free_rdp

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

    subroutine free_csp(self)
        implicit none(type, external)
        !! Default `free`: no-op. Correct for types whose only storage is a managed Fortran
        !! allocatable, which is released automatically on deallocation / scope exit. Types
        !! holding unmanaged resources override this to release them, guarded for idempotency.
        !! Clears the `is_initialized` flag so that the action is terminal.
        class(abstract_vector_csp), intent(inout) :: self
        self%is_initialized = .false.
    end subroutine free_csp

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

    subroutine free_cdp(self)
        implicit none(type, external)
        !! Default `free`: no-op. Correct for types whose only storage is a managed Fortran
        !! allocatable, which is released automatically on deallocation / scope exit. Types
        !! holding unmanaged resources override this to release them, guarded for idempotency.
        !! Clears the `is_initialized` flag so that the action is terminal.
        class(abstract_vector_cdp), intent(inout) :: self
        self%is_initialized = .false.
    end subroutine free_cdp


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
            call check_allocation(iostat, errmsg, this_module, "dense_zero_rsp")
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
            call check_allocation(iostat, errmsg, this_module, "dense_axpby_rsp")
        endif
        n = self%get_size()
        if (m /= n) call stop_error("Inconsistent size between the two vectors.")

        select type (vec)
        type is(dense_vector_rsp)
            if (beta == zero_rsp) then
                self%data = zero_rsp
            else
                call self%scal(beta)
            end if
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

    subroutine dense_init_like_rsp(self, mold)
        implicit none(type, external)
        !! Shape `self` like `mold` (three cases: allocate / no-op / reallocate).
        class(dense_vector_rsp), intent(inout) :: self
        class(abstract_vector_rsp), intent(in)  :: mold
        integer :: m, iostat
        character(len=100) :: errmsg

        ! mold must be the same concrete type to read its size.
        select type (mold)
        type is (dense_vector_rsp)
            m = mold%get_size()
        class default
            call type_error('mold','dense_vector_rsp','IN',this_module,'dense_init_like_rsp')
        end select

        if (.not. allocated(self%data)) then
            ! Case 1: unallocated -> allocate.
            allocate(self%data(m), stat=iostat, errmsg=errmsg)
            call check_allocation(iostat, errmsg, this_module, "dense_init_like_rsp")
        else if (size(self%data) /= m) then
            ! Case 3: allocated but nonconformant -> release, reallocate.
            deallocate(self%data)
            allocate(self%data(m), stat=iostat, errmsg=errmsg)
            call check_allocation(iostat, errmsg, this_module, "dense_init_like_rsp")
        end if
        ! Case 2: allocated and conformant -> fall through, no realloc.

        self%n = m
        self%is_initialized = .true.
    end subroutine dense_init_like_rsp

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
            call check_allocation(iostat, errmsg, this_module, "dense_zero_rdp")
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
            call check_allocation(iostat, errmsg, this_module, "dense_axpby_rdp")
        endif
        n = self%get_size()
        if (m /= n) call stop_error("Inconsistent size between the two vectors.")

        select type (vec)
        type is(dense_vector_rdp)
            if (beta == zero_rdp) then
                self%data = zero_rdp
            else
                call self%scal(beta)
            end if
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

    subroutine dense_init_like_rdp(self, mold)
        implicit none(type, external)
        !! Shape `self` like `mold` (three cases: allocate / no-op / reallocate).
        class(dense_vector_rdp), intent(inout) :: self
        class(abstract_vector_rdp), intent(in)  :: mold
        integer :: m, iostat
        character(len=100) :: errmsg

        ! mold must be the same concrete type to read its size.
        select type (mold)
        type is (dense_vector_rdp)
            m = mold%get_size()
        class default
            call type_error('mold','dense_vector_rdp','IN',this_module,'dense_init_like_rdp')
        end select

        if (.not. allocated(self%data)) then
            ! Case 1: unallocated -> allocate.
            allocate(self%data(m), stat=iostat, errmsg=errmsg)
            call check_allocation(iostat, errmsg, this_module, "dense_init_like_rdp")
        else if (size(self%data) /= m) then
            ! Case 3: allocated but nonconformant -> release, reallocate.
            deallocate(self%data)
            allocate(self%data(m), stat=iostat, errmsg=errmsg)
            call check_allocation(iostat, errmsg, this_module, "dense_init_like_rdp")
        end if
        ! Case 2: allocated and conformant -> fall through, no realloc.

        self%n = m
        self%is_initialized = .true.
    end subroutine dense_init_like_rdp

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
            call check_allocation(iostat, errmsg, this_module, "dense_zero_csp")
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
        call check_allocation(iostat, errmsg, this_module, "dense_rand_csp")
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
            call check_allocation(iostat, errmsg, this_module, "dense_axpby_csp")
        endif
        n = self%get_size()
        if (m /= n) call stop_error("Inconsistent size between the two vectors.")

        select type (vec)
        type is(dense_vector_csp)
            if (beta == zero_csp) then
                self%data = zero_csp
            else
                call self%scal(beta)
            end if
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

    subroutine dense_init_like_csp(self, mold)
        implicit none(type, external)
        !! Shape `self` like `mold` (three cases: allocate / no-op / reallocate).
        class(dense_vector_csp), intent(inout) :: self
        class(abstract_vector_csp), intent(in)  :: mold
        integer :: m, iostat
        character(len=100) :: errmsg

        ! mold must be the same concrete type to read its size.
        select type (mold)
        type is (dense_vector_csp)
            m = mold%get_size()
        class default
            call type_error('mold','dense_vector_csp','IN',this_module,'dense_init_like_csp')
        end select

        if (.not. allocated(self%data)) then
            ! Case 1: unallocated -> allocate.
            allocate(self%data(m), stat=iostat, errmsg=errmsg)
            call check_allocation(iostat, errmsg, this_module, "dense_init_like_csp")
        else if (size(self%data) /= m) then
            ! Case 3: allocated but nonconformant -> release, reallocate.
            deallocate(self%data)
            allocate(self%data(m), stat=iostat, errmsg=errmsg)
            call check_allocation(iostat, errmsg, this_module, "dense_init_like_csp")
        end if
        ! Case 2: allocated and conformant -> fall through, no realloc.

        self%n = m
        self%is_initialized = .true.
    end subroutine dense_init_like_csp

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
            call check_allocation(iostat, errmsg, this_module, "dense_zero_cdp")
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
        call check_allocation(iostat, errmsg, this_module, "dense_rand_cdp")
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
            call check_allocation(iostat, errmsg, this_module, "dense_axpby_cdp")
        endif
        n = self%get_size()
        if (m /= n) call stop_error("Inconsistent size between the two vectors.")

        select type (vec)
        type is(dense_vector_cdp)
            if (beta == zero_cdp) then
                self%data = zero_cdp
            else
                call self%scal(beta)
            end if
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

    subroutine dense_init_like_cdp(self, mold)
        implicit none(type, external)
        !! Shape `self` like `mold` (three cases: allocate / no-op / reallocate).
        class(dense_vector_cdp), intent(inout) :: self
        class(abstract_vector_cdp), intent(in)  :: mold
        integer :: m, iostat
        character(len=100) :: errmsg

        ! mold must be the same concrete type to read its size.
        select type (mold)
        type is (dense_vector_cdp)
            m = mold%get_size()
        class default
            call type_error('mold','dense_vector_cdp','IN',this_module,'dense_init_like_cdp')
        end select

        if (.not. allocated(self%data)) then
            ! Case 1: unallocated -> allocate.
            allocate(self%data(m), stat=iostat, errmsg=errmsg)
            call check_allocation(iostat, errmsg, this_module, "dense_init_like_cdp")
        else if (size(self%data) /= m) then
            ! Case 3: allocated but nonconformant -> release, reallocate.
            deallocate(self%data)
            allocate(self%data(m), stat=iostat, errmsg=errmsg)
            call check_allocation(iostat, errmsg, this_module, "dense_init_like_cdp")
        end if
        ! Case 2: allocated and conformant -> fall through, no realloc.

        self%n = m
        self%is_initialized = .true.
    end subroutine dense_init_like_cdp


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
        allocate(y, mold=X(1), stat=iostat, errmsg=errmsg)
        call check_allocation(iostat, errmsg, this_module, "linear_combination_vector_rsp")
        call y%init_like(X(1))
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
        allocate(Y(size(B, 2)), mold=X(1), stat=iostat, errmsg=errmsg)
        call check_allocation(iostat, errmsg, this_module, "linear_combination_matrix_rsp")
        call init_like_basis(Y, X(1))

        call zero_basis(Y)
        do j = 1, size(Y)
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
        class(abstract_vector_rsp), intent(inout) :: out
        ! Reset output based on input.
        call out%init_like(from)
        ! Copy array.
        call out%axpby(one_rsp, from, zero_rsp)
    end subroutine copy_vector_rsp

    impure elemental subroutine init_like_basis_rsp(X, mold)
        implicit none(type, external)
        class(abstract_vector_rsp), intent(inout) :: X
        class(abstract_vector_rsp), intent(in)    :: mold
        call X%init_like(mold)
    end subroutine init_like_basis_rsp

    impure elemental subroutine free_basis_rsp(X)
        implicit none(type, external)
        class(abstract_vector_rsp), intent(inout) :: X
        call X%free()
    end subroutine free_basis_rsp
    
    impure elemental subroutine rand_basis_rsp(X, ifnorm)
        implicit none(type, external)
        class(abstract_vector_rsp), intent(inout) :: X
        logical, optional, intent(in) :: ifnorm
        call X%rand(ifnorm=ifnorm)
    end subroutine rand_basis_rsp

    logical function verify_vector_axioms_rsp(x, ntrials, tolerance) result(success)
        implicit none(type, external)
        class(abstract_vector_rsp), intent(in) :: x
        !! Derived-type whose implementation needs to be tested.
        integer, optional, intent(in) :: ntrials
        !! Number of random samples generated for the tests.
        real(sp), optional, intent(in) :: tolerance

        integer :: ntrials_, i
        real(sp) :: tol
        character(len=128) :: failed_test
        character(len=256) :: msg

        !> Deals with optional argument.
        ntrials_ = optval(ntrials, 100)
        tol = optval(tolerance, 10.0_sp**(-(precision(1.0_sp)-1)))

        !> Run all tests to verify axioms.
        success = .false.
        verification: do i = 1, ntrials_

            !-----------------------------------
            !-----     VECTOR ADDITION     -----
            !-----------------------------------

            addition_distributivity: block
                class(abstract_vector_rsp), allocatable :: u, v, w
                class(abstract_vector_rsp), allocatable :: wrk1, wrk2

                !> Generate random vectors.
                allocate(u, v, w, wrk1, wrk2, mold=x)
                call u%init_like(x) ; call u%rand()
                call v%init_like(x) ; call v%rand()
                call w%init_like(x) ; call w%rand()

                !> Check distributivity.
                call copy(wrk1, v)
                call copy(wrk2, v)
                call wrk1%add(w)    ! v + w
                call wrk2%add(u)    ! u + v

                call u%add(wrk1)    ! u + (v + w)
                call w%add(wrk2)    ! (u + v) + w

                call u%sub(w)

                !> Check correctness.
                success = merge(.true., .false., u%norm() <= tol)

                !> Cleanup.
                call wrk1%free() ; call wrk2%free()
                call u%free() ; call v%free() ; call w%free()

                !> Exit if the test fails.
                if (.not. success) then
                    failed_test = 'addition_distributivity'
                    exit verification
                end if
            end block addition_distributivity

            addition_commutativity: block
                class(abstract_vector_rsp), allocatable :: u, v, w

                !> Generate random vectors.
                allocate(u, v, w, mold=x)
                call u%init_like(x) ; call u%rand()
                call v%init_like(x) ; call v%rand()
                call copy(w, v)

                !> Check commutativity.
                call v%add(u)
                call u%add(w)

                call u%sub(v)

                !> Check correctness.
                success = merge(.true., .false., u%norm() <= tol)

                !> Cleanup.
                call u%free() ; call v%free() ; call w%free()

                !> Exit if the test fails.
                if (.not. success) then
                    failed_test = 'addition_commutativity'
                    exit verification
                end if
            end block addition_commutativity

            addition_zero: block
                class(abstract_vector_rsp), allocatable :: u, v, z

                !> Generate random vector.
                allocate(u, v, z, mold=x)
                call u%init_like(x) ; call u%rand()
                call copy(v, u)
                call z%init_like(x) ; call z%zero()

                !> Check zero element.
                call u%add(z)
                call u%sub(v)

                !> Check correctness
                success = merge(.true., .false., u%norm() <= tol)

                !> Cleanup.
                call u%free() ; call v%free() ; call z%free()

                !> Exit if the test fails.
                if (.not. success) then
                    failed_test = 'addition_zero'
                    exit verification
                end if
            end block addition_zero

            additive_inverse: block
                class(abstract_vector_rsp), allocatable :: u, v
                allocate(u, v, mold=x)
                call u%init_like(x) ; call u%rand()
                call copy(v, u)
                call u%sub(v)

                !> Check correctness.
                success = merge(.true., .false., u%norm() <= tol)

                !> Cleanup.
                call u%free() ; call v%free()

                !> Exit if the test fails.
                if (.not. success) then
                    failed_test = 'additive_inverse'
                    exit verification
                end if
            end block additive_inverse

            !-----------------------------------------
            !-----     SCALAR MULTIPLICATION     -----
            !-----------------------------------------
            scaling_identity: block
                class(abstract_vector_rsp), allocatable :: u, v
                real(sp), parameter :: one = 1.0_sp

                !> Generate random vector.
                allocate(u, v, mold=x)
                call u%init_like(x) ; call u%rand()
                call copy(v, u)
                
                call v%scal(one)
                call u%sub(v)

                !> Check correctness
                success = merge(.true., .false., u%norm() <= tol)
                
                !> Cleanup.
                call u%free() ; call v%free()
                
                !> Exit if the test fails.
                if (.not. success) then
                    failed_test = 'scaling_identity'
                    exit verification
                end if
            end block scaling_identity

            scaling_compatibility: block
                class(abstract_vector_rsp), allocatable :: u, v
                real(sp) :: a, b
                call random_number(a)
                call random_number(b)

                !> Generate random vectors.
                allocate(u, v, mold=x)
                call u%init_like(x) ; call u%rand(ifnorm=.true.)
                call copy(v, u)

                !> Check associativity.
                call v%scal(b)
                call v%scal(a)
                call u%scal(a*b)
                call u%sub(v)

                !> Check correctness.
                success = merge(.true., .false., u%norm() <= tol)

                !> Cleanup.
                call u%free() ; call v%free()
                
                !> Exit if the test fails.
                if (.not. success) then
                    failed_test = 'scaling_compatibility'
                    exit verification
                end if
            end block scaling_compatibility

            scaling_distributivity: block
                class(abstract_vector_rsp), allocatable :: u, v, w
                real(sp) :: a
                call random_number(a)

                !> Generate random vectors.
                allocate(u, v, w, mold=x)
                call u%init_like(x) ; call u%rand()
                call v%init_like(x) ; call v%rand()
                call copy(w, u)

                !> Check distributivity.
                call w%add(v)
                call w%scal(a)

                call u%scal(a)
                call v%scal(a)
                call v%add(u)

                call v%sub(w)

                !> Check correctness.
                success = merge(.true., .false., v%norm() <= tol)

                !> Cleanup.
                call u%free() ; call v%free() ; call w%free()

                !> Exit if the test fails.
                if (.not. success) then
                    failed_test = 'scaling_distributivity'
                    exit verification
                end if
            end block scaling_distributivity

            scaling_distributivity_bis: block
                class(abstract_vector_rsp), allocatable :: u, v
                real(sp) :: a, b
                call random_number(a)
                call random_number(b)

                !> Generate random vector.
                allocate(u, v, mold=x)
                call u%init_like(x) ; call u%rand()
                call copy(v, u)

                !> Check distributivity.
                call v%axpby(a, u, b)
                call u%scal(a+b)

                call v%sub(u)

                !> Check correctness.
                success = merge(.true., .false., v%norm() <= tol)

                !> Cleanup.
                call u%free() ; call v%free()

                !>  Exit if the test fails.
                if (.not. success) then
                    failed_test = 'scaling_distributivity_bis'
                    exit verification
                end if
            end block scaling_distributivity_bis

            mold_independence: block
                !! guard for resource-sharing overrides (GPU), void for CPU applications.
                class(abstract_vector_rsp), allocatable :: u
                real(sp) :: xnorm

                !> Generate random vector.
                xnorm = x%norm()
                allocate(u, mold=x)
                call u%init_like(x) ; call u%rand()

                ! x must be unaffected by operations on u (checking init_like compliance).
                call u%scal(2*one_rsp)

                !> Check correctness.
                success = merge(.true., .false., abs(x%norm() - xnorm) <= tol)

                ! Cleanup.
                call u%free()

                !> Exit if the test fails.
                if (.not. success) then
                    failed_test = 'mold_independence'
                    exit verification
                end if
            end block mold_independence
        enddo verification
        if (success) then
            write(msg, '(A,I0,A)') 'All vector axioms verified (', ntrials_, ' trials).'
            call log_information(msg, this_module, 'verify_vector_axioms_rsp')
        else
            write(msg, '(A,I0,A)') 'Vector axiom check FAILED at trial ', i, ', test: '//trim(failed_test)
            call log_warning(msg, this_module, 'verify_vector_axioms_rsp')
        end if
    end function verify_vector_axioms_rsp

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
        allocate(y, mold=X(1), stat=iostat, errmsg=errmsg)
        call check_allocation(iostat, errmsg, this_module, "linear_combination_vector_rdp")
        call y%init_like(X(1))
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
        allocate(Y(size(B, 2)), mold=X(1), stat=iostat, errmsg=errmsg)
        call check_allocation(iostat, errmsg, this_module, "linear_combination_matrix_rdp")
        call init_like_basis(Y, X(1))

        call zero_basis(Y)
        do j = 1, size(Y)
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
        class(abstract_vector_rdp), intent(inout) :: out
        ! Reset output based on input.
        call out%init_like(from)
        ! Copy array.
        call out%axpby(one_rdp, from, zero_rdp)
    end subroutine copy_vector_rdp

    impure elemental subroutine init_like_basis_rdp(X, mold)
        implicit none(type, external)
        class(abstract_vector_rdp), intent(inout) :: X
        class(abstract_vector_rdp), intent(in)    :: mold
        call X%init_like(mold)
    end subroutine init_like_basis_rdp

    impure elemental subroutine free_basis_rdp(X)
        implicit none(type, external)
        class(abstract_vector_rdp), intent(inout) :: X
        call X%free()
    end subroutine free_basis_rdp
    
    impure elemental subroutine rand_basis_rdp(X, ifnorm)
        implicit none(type, external)
        class(abstract_vector_rdp), intent(inout) :: X
        logical, optional, intent(in) :: ifnorm
        call X%rand(ifnorm=ifnorm)
    end subroutine rand_basis_rdp

    logical function verify_vector_axioms_rdp(x, ntrials, tolerance) result(success)
        implicit none(type, external)
        class(abstract_vector_rdp), intent(in) :: x
        !! Derived-type whose implementation needs to be tested.
        integer, optional, intent(in) :: ntrials
        !! Number of random samples generated for the tests.
        real(dp), optional, intent(in) :: tolerance

        integer :: ntrials_, i
        real(dp) :: tol
        character(len=128) :: failed_test
        character(len=256) :: msg

        !> Deals with optional argument.
        ntrials_ = optval(ntrials, 100)
        tol = optval(tolerance, 10.0_dp**(-(precision(1.0_dp)-1)))

        !> Run all tests to verify axioms.
        success = .false.
        verification: do i = 1, ntrials_

            !-----------------------------------
            !-----     VECTOR ADDITION     -----
            !-----------------------------------

            addition_distributivity: block
                class(abstract_vector_rdp), allocatable :: u, v, w
                class(abstract_vector_rdp), allocatable :: wrk1, wrk2

                !> Generate random vectors.
                allocate(u, v, w, wrk1, wrk2, mold=x)
                call u%init_like(x) ; call u%rand()
                call v%init_like(x) ; call v%rand()
                call w%init_like(x) ; call w%rand()

                !> Check distributivity.
                call copy(wrk1, v)
                call copy(wrk2, v)
                call wrk1%add(w)    ! v + w
                call wrk2%add(u)    ! u + v

                call u%add(wrk1)    ! u + (v + w)
                call w%add(wrk2)    ! (u + v) + w

                call u%sub(w)

                !> Check correctness.
                success = merge(.true., .false., u%norm() <= tol)

                !> Cleanup.
                call wrk1%free() ; call wrk2%free()
                call u%free() ; call v%free() ; call w%free()

                !> Exit if the test fails.
                if (.not. success) then
                    failed_test = 'addition_distributivity'
                    exit verification
                end if
            end block addition_distributivity

            addition_commutativity: block
                class(abstract_vector_rdp), allocatable :: u, v, w

                !> Generate random vectors.
                allocate(u, v, w, mold=x)
                call u%init_like(x) ; call u%rand()
                call v%init_like(x) ; call v%rand()
                call copy(w, v)

                !> Check commutativity.
                call v%add(u)
                call u%add(w)

                call u%sub(v)

                !> Check correctness.
                success = merge(.true., .false., u%norm() <= tol)

                !> Cleanup.
                call u%free() ; call v%free() ; call w%free()

                !> Exit if the test fails.
                if (.not. success) then
                    failed_test = 'addition_commutativity'
                    exit verification
                end if
            end block addition_commutativity

            addition_zero: block
                class(abstract_vector_rdp), allocatable :: u, v, z

                !> Generate random vector.
                allocate(u, v, z, mold=x)
                call u%init_like(x) ; call u%rand()
                call copy(v, u)
                call z%init_like(x) ; call z%zero()

                !> Check zero element.
                call u%add(z)
                call u%sub(v)

                !> Check correctness
                success = merge(.true., .false., u%norm() <= tol)

                !> Cleanup.
                call u%free() ; call v%free() ; call z%free()

                !> Exit if the test fails.
                if (.not. success) then
                    failed_test = 'addition_zero'
                    exit verification
                end if
            end block addition_zero

            additive_inverse: block
                class(abstract_vector_rdp), allocatable :: u, v
                allocate(u, v, mold=x)
                call u%init_like(x) ; call u%rand()
                call copy(v, u)
                call u%sub(v)

                !> Check correctness.
                success = merge(.true., .false., u%norm() <= tol)

                !> Cleanup.
                call u%free() ; call v%free()

                !> Exit if the test fails.
                if (.not. success) then
                    failed_test = 'additive_inverse'
                    exit verification
                end if
            end block additive_inverse

            !-----------------------------------------
            !-----     SCALAR MULTIPLICATION     -----
            !-----------------------------------------
            scaling_identity: block
                class(abstract_vector_rdp), allocatable :: u, v
                real(dp), parameter :: one = 1.0_dp

                !> Generate random vector.
                allocate(u, v, mold=x)
                call u%init_like(x) ; call u%rand()
                call copy(v, u)
                
                call v%scal(one)
                call u%sub(v)

                !> Check correctness
                success = merge(.true., .false., u%norm() <= tol)
                
                !> Cleanup.
                call u%free() ; call v%free()
                
                !> Exit if the test fails.
                if (.not. success) then
                    failed_test = 'scaling_identity'
                    exit verification
                end if
            end block scaling_identity

            scaling_compatibility: block
                class(abstract_vector_rdp), allocatable :: u, v
                real(dp) :: a, b
                call random_number(a)
                call random_number(b)

                !> Generate random vectors.
                allocate(u, v, mold=x)
                call u%init_like(x) ; call u%rand(ifnorm=.true.)
                call copy(v, u)

                !> Check associativity.
                call v%scal(b)
                call v%scal(a)
                call u%scal(a*b)
                call u%sub(v)

                !> Check correctness.
                success = merge(.true., .false., u%norm() <= tol)

                !> Cleanup.
                call u%free() ; call v%free()
                
                !> Exit if the test fails.
                if (.not. success) then
                    failed_test = 'scaling_compatibility'
                    exit verification
                end if
            end block scaling_compatibility

            scaling_distributivity: block
                class(abstract_vector_rdp), allocatable :: u, v, w
                real(dp) :: a
                call random_number(a)

                !> Generate random vectors.
                allocate(u, v, w, mold=x)
                call u%init_like(x) ; call u%rand()
                call v%init_like(x) ; call v%rand()
                call copy(w, u)

                !> Check distributivity.
                call w%add(v)
                call w%scal(a)

                call u%scal(a)
                call v%scal(a)
                call v%add(u)

                call v%sub(w)

                !> Check correctness.
                success = merge(.true., .false., v%norm() <= tol)

                !> Cleanup.
                call u%free() ; call v%free() ; call w%free()

                !> Exit if the test fails.
                if (.not. success) then
                    failed_test = 'scaling_distributivity'
                    exit verification
                end if
            end block scaling_distributivity

            scaling_distributivity_bis: block
                class(abstract_vector_rdp), allocatable :: u, v
                real(dp) :: a, b
                call random_number(a)
                call random_number(b)

                !> Generate random vector.
                allocate(u, v, mold=x)
                call u%init_like(x) ; call u%rand()
                call copy(v, u)

                !> Check distributivity.
                call v%axpby(a, u, b)
                call u%scal(a+b)

                call v%sub(u)

                !> Check correctness.
                success = merge(.true., .false., v%norm() <= tol)

                !> Cleanup.
                call u%free() ; call v%free()

                !>  Exit if the test fails.
                if (.not. success) then
                    failed_test = 'scaling_distributivity_bis'
                    exit verification
                end if
            end block scaling_distributivity_bis

            mold_independence: block
                !! guard for resource-sharing overrides (GPU), void for CPU applications.
                class(abstract_vector_rdp), allocatable :: u
                real(dp) :: xnorm

                !> Generate random vector.
                xnorm = x%norm()
                allocate(u, mold=x)
                call u%init_like(x) ; call u%rand()

                ! x must be unaffected by operations on u (checking init_like compliance).
                call u%scal(2*one_rdp)

                !> Check correctness.
                success = merge(.true., .false., abs(x%norm() - xnorm) <= tol)

                ! Cleanup.
                call u%free()

                !> Exit if the test fails.
                if (.not. success) then
                    failed_test = 'mold_independence'
                    exit verification
                end if
            end block mold_independence
        enddo verification
        if (success) then
            write(msg, '(A,I0,A)') 'All vector axioms verified (', ntrials_, ' trials).'
            call log_information(msg, this_module, 'verify_vector_axioms_rdp')
        else
            write(msg, '(A,I0,A)') 'Vector axiom check FAILED at trial ', i, ', test: '//trim(failed_test)
            call log_warning(msg, this_module, 'verify_vector_axioms_rdp')
        end if
    end function verify_vector_axioms_rdp

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
        allocate(y, mold=X(1), stat=iostat, errmsg=errmsg)
        call check_allocation(iostat, errmsg, this_module, "linear_combination_vector_csp")
        call y%init_like(X(1))
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
        allocate(Y(size(B, 2)), mold=X(1), stat=iostat, errmsg=errmsg)
        call check_allocation(iostat, errmsg, this_module, "linear_combination_matrix_csp")
        call init_like_basis(Y, X(1))

        call zero_basis(Y)
        do j = 1, size(Y)
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
        class(abstract_vector_csp), intent(inout) :: out
        ! Reset output based on input.
        call out%init_like(from)
        ! Copy array.
        call out%axpby(one_csp, from, zero_csp)
    end subroutine copy_vector_csp

    impure elemental subroutine init_like_basis_csp(X, mold)
        implicit none(type, external)
        class(abstract_vector_csp), intent(inout) :: X
        class(abstract_vector_csp), intent(in)    :: mold
        call X%init_like(mold)
    end subroutine init_like_basis_csp

    impure elemental subroutine free_basis_csp(X)
        implicit none(type, external)
        class(abstract_vector_csp), intent(inout) :: X
        call X%free()
    end subroutine free_basis_csp
    
    impure elemental subroutine rand_basis_csp(X, ifnorm)
        implicit none(type, external)
        class(abstract_vector_csp), intent(inout) :: X
        logical, optional, intent(in) :: ifnorm
        call X%rand(ifnorm=ifnorm)
    end subroutine rand_basis_csp

    logical function verify_vector_axioms_csp(x, ntrials, tolerance) result(success)
        implicit none(type, external)
        class(abstract_vector_csp), intent(in) :: x
        !! Derived-type whose implementation needs to be tested.
        integer, optional, intent(in) :: ntrials
        !! Number of random samples generated for the tests.
        real(sp), optional, intent(in) :: tolerance

        integer :: ntrials_, i
        real(sp) :: tol
        character(len=128) :: failed_test
        character(len=256) :: msg

        !> Deals with optional argument.
        ntrials_ = optval(ntrials, 100)
        tol = optval(tolerance, 10.0_sp**(-(precision(1.0_sp)-1)))

        !> Run all tests to verify axioms.
        success = .false.
        verification: do i = 1, ntrials_

            !-----------------------------------
            !-----     VECTOR ADDITION     -----
            !-----------------------------------

            addition_distributivity: block
                class(abstract_vector_csp), allocatable :: u, v, w
                class(abstract_vector_csp), allocatable :: wrk1, wrk2

                !> Generate random vectors.
                allocate(u, v, w, wrk1, wrk2, mold=x)
                call u%init_like(x) ; call u%rand()
                call v%init_like(x) ; call v%rand()
                call w%init_like(x) ; call w%rand()

                !> Check distributivity.
                call copy(wrk1, v)
                call copy(wrk2, v)
                call wrk1%add(w)    ! v + w
                call wrk2%add(u)    ! u + v

                call u%add(wrk1)    ! u + (v + w)
                call w%add(wrk2)    ! (u + v) + w

                call u%sub(w)

                !> Check correctness.
                success = merge(.true., .false., u%norm() <= tol)

                !> Cleanup.
                call wrk1%free() ; call wrk2%free()
                call u%free() ; call v%free() ; call w%free()

                !> Exit if the test fails.
                if (.not. success) then
                    failed_test = 'addition_distributivity'
                    exit verification
                end if
            end block addition_distributivity

            addition_commutativity: block
                class(abstract_vector_csp), allocatable :: u, v, w

                !> Generate random vectors.
                allocate(u, v, w, mold=x)
                call u%init_like(x) ; call u%rand()
                call v%init_like(x) ; call v%rand()
                call copy(w, v)

                !> Check commutativity.
                call v%add(u)
                call u%add(w)

                call u%sub(v)

                !> Check correctness.
                success = merge(.true., .false., u%norm() <= tol)

                !> Cleanup.
                call u%free() ; call v%free() ; call w%free()

                !> Exit if the test fails.
                if (.not. success) then
                    failed_test = 'addition_commutativity'
                    exit verification
                end if
            end block addition_commutativity

            addition_zero: block
                class(abstract_vector_csp), allocatable :: u, v, z

                !> Generate random vector.
                allocate(u, v, z, mold=x)
                call u%init_like(x) ; call u%rand()
                call copy(v, u)
                call z%init_like(x) ; call z%zero()

                !> Check zero element.
                call u%add(z)
                call u%sub(v)

                !> Check correctness
                success = merge(.true., .false., u%norm() <= tol)

                !> Cleanup.
                call u%free() ; call v%free() ; call z%free()

                !> Exit if the test fails.
                if (.not. success) then
                    failed_test = 'addition_zero'
                    exit verification
                end if
            end block addition_zero

            additive_inverse: block
                class(abstract_vector_csp), allocatable :: u, v
                allocate(u, v, mold=x)
                call u%init_like(x) ; call u%rand()
                call copy(v, u)
                call u%sub(v)

                !> Check correctness.
                success = merge(.true., .false., u%norm() <= tol)

                !> Cleanup.
                call u%free() ; call v%free()

                !> Exit if the test fails.
                if (.not. success) then
                    failed_test = 'additive_inverse'
                    exit verification
                end if
            end block additive_inverse

            !-----------------------------------------
            !-----     SCALAR MULTIPLICATION     -----
            !-----------------------------------------
            scaling_identity: block
                class(abstract_vector_csp), allocatable :: u, v
                complex(sp), parameter :: one = 1.0_sp

                !> Generate random vector.
                allocate(u, v, mold=x)
                call u%init_like(x) ; call u%rand()
                call copy(v, u)
                
                call v%scal(one)
                call u%sub(v)

                !> Check correctness
                success = merge(.true., .false., u%norm() <= tol)
                
                !> Cleanup.
                call u%free() ; call v%free()
                
                !> Exit if the test fails.
                if (.not. success) then
                    failed_test = 'scaling_identity'
                    exit verification
                end if
            end block scaling_identity

            scaling_compatibility: block
                class(abstract_vector_csp), allocatable :: u, v
                complex(sp) :: a, b
                real(sp) :: c(2)
                call random_number(c)
                a = cmplx(c(1), c(2), kind=sp)
                call random_number(c)
                b = cmplx(c(1), c(2), kind=sp)

                !> Generate random vectors.
                allocate(u, v, mold=x)
                call u%init_like(x) ; call u%rand(ifnorm=.true.)
                call copy(v, u)

                !> Check associativity.
                call v%scal(b)
                call v%scal(a)
                call u%scal(a*b)
                call u%sub(v)

                !> Check correctness.
                success = merge(.true., .false., u%norm() <= tol)

                !> Cleanup.
                call u%free() ; call v%free()
                
                !> Exit if the test fails.
                if (.not. success) then
                    failed_test = 'scaling_compatibility'
                    exit verification
                end if
            end block scaling_compatibility

            scaling_distributivity: block
                class(abstract_vector_csp), allocatable :: u, v, w
                complex(sp) :: a
                real(sp) :: b(2)
                call random_number(b)
                a = cmplx(b(1), b(2), kind=sp)

                !> Generate random vectors.
                allocate(u, v, w, mold=x)
                call u%init_like(x) ; call u%rand()
                call v%init_like(x) ; call v%rand()
                call copy(w, u)

                !> Check distributivity.
                call w%add(v)
                call w%scal(a)

                call u%scal(a)
                call v%scal(a)
                call v%add(u)

                call v%sub(w)

                !> Check correctness.
                success = merge(.true., .false., v%norm() <= tol)

                !> Cleanup.
                call u%free() ; call v%free() ; call w%free()

                !> Exit if the test fails.
                if (.not. success) then
                    failed_test = 'scaling_distributivity'
                    exit verification
                end if
            end block scaling_distributivity

            scaling_distributivity_bis: block
                class(abstract_vector_csp), allocatable :: u, v
                complex(sp) :: a, b
                real(sp) :: c(2)
                call random_number(c)
                a = cmplx(c(1), c(2), kind=sp)
                call random_number(c)
                b = cmplx(c(1), c(2), kind=sp)

                !> Generate random vector.
                allocate(u, v, mold=x)
                call u%init_like(x) ; call u%rand()
                call copy(v, u)

                !> Check distributivity.
                call v%axpby(a, u, b)
                call u%scal(a+b)

                call v%sub(u)

                !> Check correctness.
                success = merge(.true., .false., v%norm() <= tol)

                !> Cleanup.
                call u%free() ; call v%free()

                !>  Exit if the test fails.
                if (.not. success) then
                    failed_test = 'scaling_distributivity_bis'
                    exit verification
                end if
            end block scaling_distributivity_bis

            mold_independence: block
                !! guard for resource-sharing overrides (GPU), void for CPU applications.
                class(abstract_vector_csp), allocatable :: u
                real(sp) :: xnorm

                !> Generate random vector.
                xnorm = x%norm()
                allocate(u, mold=x)
                call u%init_like(x) ; call u%rand()

                ! x must be unaffected by operations on u (checking init_like compliance).
                call u%scal(2*one_csp)

                !> Check correctness.
                success = merge(.true., .false., abs(x%norm() - xnorm) <= tol)

                ! Cleanup.
                call u%free()

                !> Exit if the test fails.
                if (.not. success) then
                    failed_test = 'mold_independence'
                    exit verification
                end if
            end block mold_independence
        enddo verification
        if (success) then
            write(msg, '(A,I0,A)') 'All vector axioms verified (', ntrials_, ' trials).'
            call log_information(msg, this_module, 'verify_vector_axioms_csp')
        else
            write(msg, '(A,I0,A)') 'Vector axiom check FAILED at trial ', i, ', test: '//trim(failed_test)
            call log_warning(msg, this_module, 'verify_vector_axioms_csp')
        end if
    end function verify_vector_axioms_csp

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
        allocate(y, mold=X(1), stat=iostat, errmsg=errmsg)
        call check_allocation(iostat, errmsg, this_module, "linear_combination_vector_cdp")
        call y%init_like(X(1))
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
        allocate(Y(size(B, 2)), mold=X(1), stat=iostat, errmsg=errmsg)
        call check_allocation(iostat, errmsg, this_module, "linear_combination_matrix_cdp")
        call init_like_basis(Y, X(1))

        call zero_basis(Y)
        do j = 1, size(Y)
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
        class(abstract_vector_cdp), intent(inout) :: out
        ! Reset output based on input.
        call out%init_like(from)
        ! Copy array.
        call out%axpby(one_cdp, from, zero_cdp)
    end subroutine copy_vector_cdp

    impure elemental subroutine init_like_basis_cdp(X, mold)
        implicit none(type, external)
        class(abstract_vector_cdp), intent(inout) :: X
        class(abstract_vector_cdp), intent(in)    :: mold
        call X%init_like(mold)
    end subroutine init_like_basis_cdp

    impure elemental subroutine free_basis_cdp(X)
        implicit none(type, external)
        class(abstract_vector_cdp), intent(inout) :: X
        call X%free()
    end subroutine free_basis_cdp
    
    impure elemental subroutine rand_basis_cdp(X, ifnorm)
        implicit none(type, external)
        class(abstract_vector_cdp), intent(inout) :: X
        logical, optional, intent(in) :: ifnorm
        call X%rand(ifnorm=ifnorm)
    end subroutine rand_basis_cdp

    logical function verify_vector_axioms_cdp(x, ntrials, tolerance) result(success)
        implicit none(type, external)
        class(abstract_vector_cdp), intent(in) :: x
        !! Derived-type whose implementation needs to be tested.
        integer, optional, intent(in) :: ntrials
        !! Number of random samples generated for the tests.
        real(dp), optional, intent(in) :: tolerance

        integer :: ntrials_, i
        real(dp) :: tol
        character(len=128) :: failed_test
        character(len=256) :: msg

        !> Deals with optional argument.
        ntrials_ = optval(ntrials, 100)
        tol = optval(tolerance, 10.0_dp**(-(precision(1.0_dp)-1)))

        !> Run all tests to verify axioms.
        success = .false.
        verification: do i = 1, ntrials_

            !-----------------------------------
            !-----     VECTOR ADDITION     -----
            !-----------------------------------

            addition_distributivity: block
                class(abstract_vector_cdp), allocatable :: u, v, w
                class(abstract_vector_cdp), allocatable :: wrk1, wrk2

                !> Generate random vectors.
                allocate(u, v, w, wrk1, wrk2, mold=x)
                call u%init_like(x) ; call u%rand()
                call v%init_like(x) ; call v%rand()
                call w%init_like(x) ; call w%rand()

                !> Check distributivity.
                call copy(wrk1, v)
                call copy(wrk2, v)
                call wrk1%add(w)    ! v + w
                call wrk2%add(u)    ! u + v

                call u%add(wrk1)    ! u + (v + w)
                call w%add(wrk2)    ! (u + v) + w

                call u%sub(w)

                !> Check correctness.
                success = merge(.true., .false., u%norm() <= tol)

                !> Cleanup.
                call wrk1%free() ; call wrk2%free()
                call u%free() ; call v%free() ; call w%free()

                !> Exit if the test fails.
                if (.not. success) then
                    failed_test = 'addition_distributivity'
                    exit verification
                end if
            end block addition_distributivity

            addition_commutativity: block
                class(abstract_vector_cdp), allocatable :: u, v, w

                !> Generate random vectors.
                allocate(u, v, w, mold=x)
                call u%init_like(x) ; call u%rand()
                call v%init_like(x) ; call v%rand()
                call copy(w, v)

                !> Check commutativity.
                call v%add(u)
                call u%add(w)

                call u%sub(v)

                !> Check correctness.
                success = merge(.true., .false., u%norm() <= tol)

                !> Cleanup.
                call u%free() ; call v%free() ; call w%free()

                !> Exit if the test fails.
                if (.not. success) then
                    failed_test = 'addition_commutativity'
                    exit verification
                end if
            end block addition_commutativity

            addition_zero: block
                class(abstract_vector_cdp), allocatable :: u, v, z

                !> Generate random vector.
                allocate(u, v, z, mold=x)
                call u%init_like(x) ; call u%rand()
                call copy(v, u)
                call z%init_like(x) ; call z%zero()

                !> Check zero element.
                call u%add(z)
                call u%sub(v)

                !> Check correctness
                success = merge(.true., .false., u%norm() <= tol)

                !> Cleanup.
                call u%free() ; call v%free() ; call z%free()

                !> Exit if the test fails.
                if (.not. success) then
                    failed_test = 'addition_zero'
                    exit verification
                end if
            end block addition_zero

            additive_inverse: block
                class(abstract_vector_cdp), allocatable :: u, v
                allocate(u, v, mold=x)
                call u%init_like(x) ; call u%rand()
                call copy(v, u)
                call u%sub(v)

                !> Check correctness.
                success = merge(.true., .false., u%norm() <= tol)

                !> Cleanup.
                call u%free() ; call v%free()

                !> Exit if the test fails.
                if (.not. success) then
                    failed_test = 'additive_inverse'
                    exit verification
                end if
            end block additive_inverse

            !-----------------------------------------
            !-----     SCALAR MULTIPLICATION     -----
            !-----------------------------------------
            scaling_identity: block
                class(abstract_vector_cdp), allocatable :: u, v
                complex(dp), parameter :: one = 1.0_dp

                !> Generate random vector.
                allocate(u, v, mold=x)
                call u%init_like(x) ; call u%rand()
                call copy(v, u)
                
                call v%scal(one)
                call u%sub(v)

                !> Check correctness
                success = merge(.true., .false., u%norm() <= tol)
                
                !> Cleanup.
                call u%free() ; call v%free()
                
                !> Exit if the test fails.
                if (.not. success) then
                    failed_test = 'scaling_identity'
                    exit verification
                end if
            end block scaling_identity

            scaling_compatibility: block
                class(abstract_vector_cdp), allocatable :: u, v
                complex(dp) :: a, b
                real(dp) :: c(2)
                call random_number(c)
                a = cmplx(c(1), c(2), kind=dp)
                call random_number(c)
                b = cmplx(c(1), c(2), kind=dp)

                !> Generate random vectors.
                allocate(u, v, mold=x)
                call u%init_like(x) ; call u%rand(ifnorm=.true.)
                call copy(v, u)

                !> Check associativity.
                call v%scal(b)
                call v%scal(a)
                call u%scal(a*b)
                call u%sub(v)

                !> Check correctness.
                success = merge(.true., .false., u%norm() <= tol)

                !> Cleanup.
                call u%free() ; call v%free()
                
                !> Exit if the test fails.
                if (.not. success) then
                    failed_test = 'scaling_compatibility'
                    exit verification
                end if
            end block scaling_compatibility

            scaling_distributivity: block
                class(abstract_vector_cdp), allocatable :: u, v, w
                complex(dp) :: a
                real(dp) :: b(2)
                call random_number(b)
                a = cmplx(b(1), b(2), kind=dp)

                !> Generate random vectors.
                allocate(u, v, w, mold=x)
                call u%init_like(x) ; call u%rand()
                call v%init_like(x) ; call v%rand()
                call copy(w, u)

                !> Check distributivity.
                call w%add(v)
                call w%scal(a)

                call u%scal(a)
                call v%scal(a)
                call v%add(u)

                call v%sub(w)

                !> Check correctness.
                success = merge(.true., .false., v%norm() <= tol)

                !> Cleanup.
                call u%free() ; call v%free() ; call w%free()

                !> Exit if the test fails.
                if (.not. success) then
                    failed_test = 'scaling_distributivity'
                    exit verification
                end if
            end block scaling_distributivity

            scaling_distributivity_bis: block
                class(abstract_vector_cdp), allocatable :: u, v
                complex(dp) :: a, b
                real(dp) :: c(2)
                call random_number(c)
                a = cmplx(c(1), c(2), kind=dp)
                call random_number(c)
                b = cmplx(c(1), c(2), kind=dp)

                !> Generate random vector.
                allocate(u, v, mold=x)
                call u%init_like(x) ; call u%rand()
                call copy(v, u)

                !> Check distributivity.
                call v%axpby(a, u, b)
                call u%scal(a+b)

                call v%sub(u)

                !> Check correctness.
                success = merge(.true., .false., v%norm() <= tol)

                !> Cleanup.
                call u%free() ; call v%free()

                !>  Exit if the test fails.
                if (.not. success) then
                    failed_test = 'scaling_distributivity_bis'
                    exit verification
                end if
            end block scaling_distributivity_bis

            mold_independence: block
                !! guard for resource-sharing overrides (GPU), void for CPU applications.
                class(abstract_vector_cdp), allocatable :: u
                real(dp) :: xnorm

                !> Generate random vector.
                xnorm = x%norm()
                allocate(u, mold=x)
                call u%init_like(x) ; call u%rand()

                ! x must be unaffected by operations on u (checking init_like compliance).
                call u%scal(2*one_cdp)

                !> Check correctness.
                success = merge(.true., .false., abs(x%norm() - xnorm) <= tol)

                ! Cleanup.
                call u%free()

                !> Exit if the test fails.
                if (.not. success) then
                    failed_test = 'mold_independence'
                    exit verification
                end if
            end block mold_independence
        enddo verification
        if (success) then
            write(msg, '(A,I0,A)') 'All vector axioms verified (', ntrials_, ' trials).'
            call log_information(msg, this_module, 'verify_vector_axioms_cdp')
        else
            write(msg, '(A,I0,A)') 'Vector axiom check FAILED at trial ', i, ', test: '//trim(failed_test)
            call log_warning(msg, this_module, 'verify_vector_axioms_cdp')
        end if
    end function verify_vector_axioms_cdp

end module LightKrylov_AbstractVectors
