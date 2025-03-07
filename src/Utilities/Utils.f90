module LightKrylov_Utils
    !!  This module provides a set of utility functions used throughout `LightKrylov`.
    !!  It includes:
    !!
    !!  - `assert_shape`: Assert that the shape of the argument is the expected shape.
    !!  - `eig`: Compute the eigenvalue decomposition of a general matrix.
    !!  - `sqrtm`: Compute the non-negative square root of a symmetric positive definite matrix using its SVD.
    !!  - `ordschur`: Re-order the Schur factorization to have the selected eigenvalues in the upper left block.
    !!
    !!  Note that as the development of `stdlib` progresses, some of these functions
    !!  will be deprecated in favor of the `stdlib` implementations.

    !--------------------------------------------
    !-----     Standard Fortran Library     -----
    !--------------------------------------------
    use iso_fortran_env, only: output_unit

    !-------------------------------
    !-----     LightKrylov     -----
    !-------------------------------
    use LightKrylov_Constants
    ! use LightKrylov_Logger
    use LightKrylov_Logger, only: log_warning, log_error, log_message, log_information, &
    &                             stop_error, check_info

    implicit none
    private

    character(len=*), parameter :: this_module      = 'LK_Utils'
    character(len=*), parameter :: this_module_long = 'LightKrylov_Utils'

    !----------------------------------
    !-----     Public exports     -----
    !----------------------------------

    public :: assert_shape
    public :: log2
    public :: eig
    public :: ordschur
    public :: sqrtm
    public :: expm
    public :: givens_rotation
    public :: apply_givens_rotation
    public :: solve_triangular

    !-------------------------------------------------
    !-----     Options for iterative solvers     -----
    !-------------------------------------------------

    type, abstract, public :: abstract_opts
        !! Abstract type for options from which all others are extended.
    end type

    type, abstract, public :: abstract_metadata
        !! Abstract type for solver metadata from which all others are extended.
        private
        contains
        procedure(abstract_print_metadata), pass(self), deferred, public :: print
        procedure(abstract_reset_metadata), pass(self), deferred, public :: reset
    end type

    abstract interface
        subroutine abstract_print_metadata(self, reset_counters, verbose)
            import abstract_metadata
            class(abstract_metadata), intent(inout) :: self
            logical, optional, intent(in) :: reset_counters
            logical, optional, intent(in) :: verbose
        end subroutine

        subroutine abstract_reset_metadata(self)
            import abstract_metadata
            class(abstract_metadata), intent(inout) :: self
        end subroutine
    end interface

    !-------------------------------------
    !-----     Utility functions     -----
    !-------------------------------------

    ! NOTE : Most of these functions will gradually disappear as more stable
    !        versions make their ways into the Fortran stdlib library.

    interface assert_shape
        !! This interface provides methods to assert tha thte shape of its input vector or
        !! matrix is as expected. It throws an error if not.
        module subroutine assert_shape_vector_rsp(v, size, vecname, module, procedure)
            real(sp), intent(in) :: v(:)
            integer, intent(in) :: size(:)
            character(len=*), intent(in) :: vecname
            character(len=*), intent(in) :: module
            character(len=*), intent(in) :: procedure
        end subroutine

        module subroutine assert_shape_matrix_rsp(A, size, matname, module, procedure)
            real(sp), intent(in) :: A(:, :)
            integer, intent(in) :: size(:)
            character(len=*), intent(in) :: matname
            character(len=*), intent(in) :: module
            character(len=*), intent(in) :: procedure
        end subroutine
        module subroutine assert_shape_vector_rdp(v, size, vecname, module, procedure)
            real(dp), intent(in) :: v(:)
            integer, intent(in) :: size(:)
            character(len=*), intent(in) :: vecname
            character(len=*), intent(in) :: module
            character(len=*), intent(in) :: procedure
        end subroutine

        module subroutine assert_shape_matrix_rdp(A, size, matname, module, procedure)
            real(dp), intent(in) :: A(:, :)
            integer, intent(in) :: size(:)
            character(len=*), intent(in) :: matname
            character(len=*), intent(in) :: module
            character(len=*), intent(in) :: procedure
        end subroutine
        module subroutine assert_shape_vector_csp(v, size, vecname, module, procedure)
            complex(sp), intent(in) :: v(:)
            integer, intent(in) :: size(:)
            character(len=*), intent(in) :: vecname
            character(len=*), intent(in) :: module
            character(len=*), intent(in) :: procedure
        end subroutine

        module subroutine assert_shape_matrix_csp(A, size, matname, module, procedure)
            complex(sp), intent(in) :: A(:, :)
            integer, intent(in) :: size(:)
            character(len=*), intent(in) :: matname
            character(len=*), intent(in) :: module
            character(len=*), intent(in) :: procedure
        end subroutine
        module subroutine assert_shape_vector_cdp(v, size, vecname, module, procedure)
            complex(dp), intent(in) :: v(:)
            integer, intent(in) :: size(:)
            character(len=*), intent(in) :: vecname
            character(len=*), intent(in) :: module
            character(len=*), intent(in) :: procedure
        end subroutine

        module subroutine assert_shape_matrix_cdp(A, size, matname, module, procedure)
            complex(dp), intent(in) :: A(:, :)
            integer, intent(in) :: size(:)
            character(len=*), intent(in) :: matname
            character(len=*), intent(in) :: module
            character(len=*), intent(in) :: procedure
        end subroutine
    end interface

    interface log2
        !! Utility function to compute the base-2 logarithm of a real number.
        elemental real(sp) module function log2_rsp(x) result(y)
            real(sp), intent(in) :: x
        end function
        elemental real(dp) module function log2_rdp(x) result(y)
            real(dp), intent(in) :: x
        end function
    end interface

    interface eig
        !!  Computes the eigenvalue decomposition of a general square matrix.
        !!
        !!  ### Description
        !!
        !!  This interface provides methods to compute the solution to the eigenproblem
        !!  \( \mathbf{Ax} = \lambda \mathbf{x} \), where $\mathbf{A}$ is a square `real`
        !!  or `complex` matrix.
        !!
        !!  Result array `lambda` returns the eigenvalues of \( \mathbf{A} \), while `vecs`
        !!  returns the corresponding eigenvectors. Note that it follows the LAPACK convention
        !!  when \( \mathbf{A} \) is `real`. The solver is based on LAPACK's `*GEEV` backends.
        !!
        !!  ### Syntax
        !!
        !!  `call eig(A, vecs, lambda)`
        !!
        !!  ### Arguments
        !!
        !!  `A`: `real` or `complex` square array containing the coefficient matrix. It is an `intent(in)` argument.
        !!
        !!  `vecs`: Square array of the same size, type, and kind as `A` containing the eigenvectors
        !!  (following LAPACK's convention for `real` matrices). It is an `intent(out)` argument.
        !!
        !!  `lambda`: `complex` rank-1 array of the same kind as `A` containing the eigenvalues.
        !!  It is an `intent(out)` argument.
        !!
        !!  @note
        !!  Due to the abstrct nature of the vector types defined in `LightKrylov`, it is unlikely
        !!  that this implementation will be superseeded in favor of the `stdlib` one as the latter
        !!  does not follow the LAPACK's convention.
        !!  @endnote
        module subroutine eig_rsp(A, vecs, vals)
            real(sp), intent(in) :: A(:, :)
            real(sp), intent(out) :: vecs(:, :)
            complex(sp), intent(out) :: vals(:)
        end subroutine
        module subroutine eig_rdp(A, vecs, vals)
            real(dp), intent(in) :: A(:, :)
            real(dp), intent(out) :: vecs(:, :)
            complex(dp), intent(out) :: vals(:)
        end subroutine
        module subroutine eig_csp(A, vecs, vals)
            complex(sp), intent(in) :: A(:, :)
            complex(sp), intent(out) :: vecs(:, :)
            complex(sp), intent(out) :: vals(:)
        end subroutine
        module subroutine eig_cdp(A, vecs, vals)
            complex(dp), intent(in) :: A(:, :)
            complex(dp), intent(out) :: vecs(:, :)
            complex(dp), intent(out) :: vals(:)
        end subroutine
    end interface

    interface ordschur
        !!  Given the Schur factorization and basis of a matrix, reorders it to have the selected
        !!  eigenvalues in the upper left block.
        !!
        !!  ### Description
        !!
        !!  This interface provides methods to re-order the Schur factorization of a `real` or
        !!  `complex` square matrix. Note that, if \( \mathbf{A} \) is `real`, it returns the
        !!  real Schur form.
        !!
        !!  ### Syntax
        !!
        !!  `call ordschur(T, Q, selected)`
        !!
        !!  ### Arguments
        !!
        !!  `T`: `real` or `complex` square array containing the Schur factorization of a matrix. 
        !!  On exit, it is overwritten with its re-ordered counterpart. It is an `intent(inout)` argument.
        !!  
        !!  `Q`: Two-dimensional square array of the same size, type and kind as `A`. It contains
        !!  the original Schur basis on entry and the re-ordered one on exit.
        !!  It is an `intent(inout)` argument.
        !!
        !!  `selected`: `logical` rank-1 array selecting which eigenvalues need to be moved in the
        !!  upper left block of the Schur factorization.
        !!  It is an `intent(in)` arguement.
        module subroutine ordschur_rsp(T, Q, selected)
            real(sp), intent(inout) :: T(:, :)
            real(sp), intent(inout) :: Q(:, :)
            logical, intent(in) :: selected(:)
        end subroutine
        module subroutine ordschur_rdp(T, Q, selected)
            real(dp), intent(inout) :: T(:, :)
            real(dp), intent(inout) :: Q(:, :)
            logical, intent(in) :: selected(:)
        end subroutine
        module subroutine ordschur_csp(T, Q, selected)
            complex(sp), intent(inout) :: T(:, :)
            complex(sp), intent(inout) :: Q(:, :)
            logical, intent(in) :: selected(:)
        end subroutine
        module subroutine ordschur_cdp(T, Q, selected)
            complex(dp), intent(inout) :: T(:, :)
            complex(dp), intent(inout) :: Q(:, :)
            logical, intent(in) :: selected(:)
        end subroutine
    end interface

    interface sqrtm
        !!  Computes the non-negative square root of a symmetric positive definite matrix
        !!  using its singular value decomposition.
        !!
        !!  ### Description
        !!
        !!  This interface provides methods to compute the non-negative square root of a symmetric
        !!  (hermitian) positive definite matrix \( \mathbf{A} \).
        !!
        !!  ### Syntax
        !!
        !!  `call sqrtm(A, sqrtmA, info)`
        !!
        !!  ### Arguments
        !!  
        !!  `A`: Symmetric (hermitian) positive definite matrix whose non-negative square root
        !!  needs to be computed. It is an `intent(in)` argument.
        !!
        !!  `sqrtmA`: Non-negative square root of `A`. It has the same size, kind and type as `A`.
        !!  It is an `intent(out)` argument.
        !!
        !!  `info`: Information flag. It is an `intent(out)` argument. 
        module subroutine sqrtm_rsp(A, sqrtA, info)
            real(sp), intent(inout) :: A(:, :)
            real(sp), intent(out) :: sqrtA(:, :)
            integer, intent(out) :: info
        end subroutine
        module subroutine sqrtm_rdp(A, sqrtA, info)
            real(dp), intent(inout) :: A(:, :)
            real(dp), intent(out) :: sqrtA(:, :)
            integer, intent(out) :: info
        end subroutine
        module subroutine sqrtm_csp(A, sqrtA, info)
            complex(sp), intent(inout) :: A(:, :)
            complex(sp), intent(out) :: sqrtA(:, :)
            integer, intent(out) :: info
        end subroutine
        module subroutine sqrtm_cdp(A, sqrtA, info)
            complex(dp), intent(inout) :: A(:, :)
            complex(dp), intent(out) :: sqrtA(:, :)
            integer, intent(out) :: info
        end subroutine
    end interface

    interface expm
        !!  ### Description
        !!
        !!  Evaluate the exponential of a dense matrix using Pade approximations.
        !!
        !!  ### Syntax
        !!
        !!  ```fortran
        !!      E = expm(A, order)
        !!  ```
        !!
        !!  ### Arguments
        !!
        !!  `E` : `real` or `complex` rank-2 array with \( E = \exp(A) \).
        !!
        !!  `A` : `real` or `complex` matrix that needs to be exponentiated.
        !!
        !!  `order` (optional) : Order of the Pade approximation. By default `order = 10`.
        module function expm_rsp(A, order) result(E)
            real(sp), intent(in) :: A(:, :)
            !! Matrix to be exponentiated.
            real(sp) :: E(size(A, 1), size(A, 1))
            !! Output matrix E = exp(tA).
            integer, intent(in), optional :: order
            !! Order of the Pade approximation.
        end function
        module function expm_rdp(A, order) result(E)
            real(dp), intent(in) :: A(:, :)
            !! Matrix to be exponentiated.
            real(dp) :: E(size(A, 1), size(A, 1))
            !! Output matrix E = exp(tA).
            integer, intent(in), optional :: order
            !! Order of the Pade approximation.
        end function
        module function expm_csp(A, order) result(E)
            complex(sp), intent(in) :: A(:, :)
            !! Matrix to be exponentiated.
            complex(sp) :: E(size(A, 1), size(A, 1))
            !! Output matrix E = exp(tA).
            integer, intent(in), optional :: order
            !! Order of the Pade approximation.
        end function
        module function expm_cdp(A, order) result(E)
            complex(dp), intent(in) :: A(:, :)
            !! Matrix to be exponentiated.
            complex(dp) :: E(size(A, 1), size(A, 1))
            !! Output matrix E = exp(tA).
            integer, intent(in), optional :: order
            !! Order of the Pade approximation.
        end function
    end interface

    interface givens_rotation
        pure module function givens_rotation_rsp(x) result(g)
            real(sp), intent(in) :: x(2)
            !! Vector whose second needs to be eliminated.
            real(sp)             :: g(2)
            !! Entries of the Givens rotation matrix.
        end function
        pure module function givens_rotation_rdp(x) result(g)
            real(dp), intent(in) :: x(2)
            !! Vector whose second needs to be eliminated.
            real(dp)             :: g(2)
            !! Entries of the Givens rotation matrix.
        end function
        pure module function givens_rotation_csp(x) result(g)
            complex(sp), intent(in) :: x(2)
            !! Vector whose second needs to be eliminated.
            complex(sp)             :: g(2)
            !! Entries of the Givens rotation matrix.
        end function
        pure module function givens_rotation_cdp(x) result(g)
            complex(dp), intent(in) :: x(2)
            !! Vector whose second needs to be eliminated.
            complex(dp)             :: g(2)
            !! Entries of the Givens rotation matrix.
        end function
    end interface

    interface apply_givens_rotation
        pure module subroutine apply_givens_rotation_rsp(h, c, s)
            real(sp), intent(inout) :: h(:)
            !! k-th column of the Hessenberg matrix.
            real(sp), intent(inout) :: c(:)
            !! Cosine components of the Givens rotations.
            real(sp), intent(inout) :: s(:)
            !! Sine components of the Givens rotations.
        end subroutine
        pure module subroutine apply_givens_rotation_rdp(h, c, s)
            real(dp), intent(inout) :: h(:)
            !! k-th column of the Hessenberg matrix.
            real(dp), intent(inout) :: c(:)
            !! Cosine components of the Givens rotations.
            real(dp), intent(inout) :: s(:)
            !! Sine components of the Givens rotations.
        end subroutine
        pure module subroutine apply_givens_rotation_csp(h, c, s)
            complex(sp), intent(inout) :: h(:)
            !! k-th column of the Hessenberg matrix.
            complex(sp), intent(inout) :: c(:)
            !! Cosine components of the Givens rotations.
            complex(sp), intent(inout) :: s(:)
            !! Sine components of the Givens rotations.
        end subroutine
        pure module subroutine apply_givens_rotation_cdp(h, c, s)
            complex(dp), intent(inout) :: h(:)
            !! k-th column of the Hessenberg matrix.
            complex(dp), intent(inout) :: c(:)
            !! Cosine components of the Givens rotations.
            complex(dp), intent(inout) :: s(:)
            !! Sine components of the Givens rotations.
        end subroutine
    end interface

    interface solve_triangular
        pure module function solve_triangular_rsp(A, b) result(x)
            real(sp), intent(in) :: A(:, :)
            !! Matrix to invert.
            real(sp), intent(in) :: b(:)
            !! Right-hand side vector.
            real(sp), allocatable :: x(:)
            !! Solution vector.
        end function
        pure module function solve_triangular_rdp(A, b) result(x)
            real(dp), intent(in) :: A(:, :)
            !! Matrix to invert.
            real(dp), intent(in) :: b(:)
            !! Right-hand side vector.
            real(dp), allocatable :: x(:)
            !! Solution vector.
        end function
        pure module function solve_triangular_csp(A, b) result(x)
            complex(sp), intent(in) :: A(:, :)
            !! Matrix to invert.
            complex(sp), intent(in) :: b(:)
            !! Right-hand side vector.
            complex(sp), allocatable :: x(:)
            !! Solution vector.
        end function
        pure module function solve_triangular_cdp(A, b) result(x)
            complex(dp), intent(in) :: A(:, :)
            !! Matrix to invert.
            complex(dp), intent(in) :: b(:)
            !! Right-hand side vector.
            complex(dp), allocatable :: x(:)
            !! Solution vector.
        end function
    end interface
end module
