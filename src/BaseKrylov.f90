module LightKrylov_BaseKrylov
    !!  This module provides a collection of Krylov-based factorizations forming the
    !!  computational core of `LightKrylov`. It also provides a set of utility functions
    !!  to operate on arrays of `abstract_vector`. The most important ones are:
    !!
    !!  - `arnoldi(A, X, H, info)`: Arnoldi factorization for general square matrices.
    !!  - `lanczos(A, X, H, info)`: Lanczos factorization for general symmetric/hermitian matrices.
    !!  - `bidiagonalization(A, U, V, B)`: Lanczos bidiagonalization for arbitrary matrices.
    !!  - `qr(X, R, perm, info)`: QR factorization (with and without column pivoting) of an array of `abstract_vector`.


    !--------------------------------------------
    !-----     Standard Fortran Library     -----
    !--------------------------------------------
    use iso_fortran_env
    use stdlib_optval, only: optval
    use stdlib_linalg, only: eye

    !-------------------------------
    !-----     LightKrylov     -----
    !-------------------------------
    use LightKrylov_Constants
    use LightKrylov_Logger
    use LightKrylov_Utils
    use LightKrylov_AbstractVectors
    use LightKrylov_AbstractLinops

    implicit none
    private
    
    character(len=128), parameter :: this_module = 'LightKrylov_BaseKrylov'

    public :: qr
    public :: apply_permutation_matrix
    public :: apply_inverse_permutation_matrix
    public :: arnoldi
    public :: initialize_krylov_subspace
    public :: is_orthonormal
    public :: orthonormalize_basis
    public :: orthogonalize_against_basis
    public :: bidiagonalization
    public :: lanczos
    public :: krylov_schur
    public :: double_gram_schmidt_step
    public :: eigvals_select_sp
    public :: eigvals_select_dp

    interface qr
        !!  ### Description
        !!
        !!  Given an array \( X \) of types derived from `abstract_vector`, it computes the
        !!  *in-place* QR factorization of \( X \), i.e.
        !!
        !!  \[
        !!      X = QR,
        !!  \]
        !!
        !!  where \( Q \) is an orthonormal arrays of vectors such that \( Q^H Q = I \) and
        !!  \( R \) is upper triangular. Note that it can also perform the QR factorization
        !!  with column pivoting
        !!
        !!  \[
        !!      XP = QR
        !!  \]
        !!
        !!  where \( P \) is a permutation matrix ensuring that the diagonal entries of \( R \)
        !!  have non-increasing absolute values. This amounts to using the pivoting QR as a
        !!  rank-revealing factorization.
        !!
        !!  **References**
        !!
        !!  - G. H. Golub & C. F. Van Loan. "Matrix Computations". 4th edition, The John Hopkins
        !!   University Press, 2013.
        !!   See Chapter 5.2.8: Modified Gram-Schmidt algorithm.
        !!
        !!  ### Syntax
        !!
        !!  ```fortran
        !!      call qr(Q [, R] [, perm], info [, tol])
        !!  ```
        !!
        !!  ### Arguments
        !!
        !!  `Q`: Array of types derived from one of the base types provided in the
        !!       `AbstractVectors` module. On entry, it contains the original array.
        !!       On exit, it is overwritten by the orthogonal basis for its span.
        !!       It is an `intent(inout)` argument.
        !!
        !!  `R`: `real` or `complex` rank-2 array. On exit, its contains the upper triangular
        !!        matrix resulting from the QR factorization. It is an `intent(out)` argument.
        !!
        !!  `perm` (*optional*): Rank-1 array of `integer` corresponding to the indices of
        !!                       permuted columns. If `perm` is absent, the naive QR factorization
        !!                       is being computed.
        !!
        !!  `info`: `integer` information flag.
        !!
        !!  `tol` (*optional*): Numerical tolerance to determine whether two vectors are colinear
        !!                      or not. Default `tol = atol_sp` or `tol = atol_dp`.
        module procedure qr_no_pivoting_rsp
        module procedure qr_with_pivoting_rsp
        module procedure qr_no_pivoting_rdp
        module procedure qr_with_pivoting_rdp
        module procedure qr_no_pivoting_csp
        module procedure qr_with_pivoting_csp
        module procedure qr_no_pivoting_cdp
        module procedure qr_with_pivoting_cdp
    end interface

    interface swap_columns
        module procedure swap_columns_rsp
        module procedure swap_columns_rdp
        module procedure swap_columns_csp
        module procedure swap_columns_cdp
    end interface

    interface apply_permutation_matrix
        !!  ### Description
        !!
        !!  Given an array \( X \) and a permutation vector \( p \), this function computes
        !!  *in-place* the column-permuted matrix
        !!
        !!  \[
        !!      X = X P
        !!  \]
        !!
        !!  where \( P \) is the column-permutation matrix constructed from the permutation
        !!  vector \( p \).
        !!
        !!  ### Syntax
        !!
        !!  ```fortran
        !!      call apply_permutation_matrix(X, perm)
        !!  ```
        !!
        !!  ### Arguments
        !!
        !!  `Q` : Array of vectors derived from the base types defined in the `AbstractVectors`
        !!       module. On entry, it is the original array. On exit, it contains the
        !!       column-permuted version computed in-place. It is an `intent(inout)` argument.
        !!
        !!  `perm` : Rank-1 array of `integer` corresponding to the desired permutation vector.
        !!          It is an `intent(in)` argument.
        module procedure apply_permutation_matrix_rsp
        module procedure apply_permutation_matrix_array_rsp
        module procedure apply_permutation_matrix_rdp
        module procedure apply_permutation_matrix_array_rdp
        module procedure apply_permutation_matrix_csp
        module procedure apply_permutation_matrix_array_csp
        module procedure apply_permutation_matrix_cdp
        module procedure apply_permutation_matrix_array_cdp
    end interface

    interface apply_inverse_permutation_matrix
        !!  ### Description
        !!
        !!  Given an array \( X \) and a permutation vector \( p \), this function computes
        !!  *in-place* the column-permuted matrix
        !!
        !!  \[
        !!      X = X P^{-1}
        !!  \]
        !!
        !!  where \( P \) is the column-permutation matrix constructed from the permutation
        !!  vector \( p \) and \( P^{-1} \) its inverse.
        !!
        !!  ### Syntax
        !!
        !!  ```fortran
        !!      call apply_inverse_permutation_matrix(X, perm)
        !!  ```
        !!
        !!  ### Arguments
        !!
        !!  `Q` : Array of vectors derived from the base types defined in the `AbstractVectors`
        !!       module. On entry, it is the original array. On exit, it contains the
        !!       column-permuted version computed in-place. It is an `intent(inout)` argument.
        !!
        !!  `perm` : Rank-1 array of `integer` corresponding to the desired permutation vector.
        !!          It is an `intent(in)` argument.
        module procedure apply_inverse_permutation_matrix_rsp
        module procedure apply_inverse_permutation_matrix_array_rsp
        module procedure apply_inverse_permutation_matrix_rdp
        module procedure apply_inverse_permutation_matrix_array_rdp
        module procedure apply_inverse_permutation_matrix_csp
        module procedure apply_inverse_permutation_matrix_array_csp
        module procedure apply_inverse_permutation_matrix_cdp
        module procedure apply_inverse_permutation_matrix_array_cdp
    end interface

    interface arnoldi
        !!  ### Description
        !!
        !!  Given a square linear operator \( A \), find matrices \( X \) and \( H \) such that
        !!
        !!  \[
        !!      AX_k = X_k H_k + h_{k+1, k} x_{k+1} e_k^T,
        !!  \]
        !!
        !!  where \( X \) is an orthogonal basis and \( H \) is upper Hessenberg.
        !!
        !!  **Algorithmic Features**
        !!
        !!  - The operator \( A \) only needs to be accessed through matrix-vector products.
        !!  - Constructs an orthonormal Krylov basis \( X \) via the Gram-Schmidt process.
        !!  - Constructs an upper Hessenberg matrix \( H \) whose eigenvalues approximates those of \( A \).
        !!  - Checks for convergence and invariant subspaces.
        !   - Block version available (experimental).
        !!
        !!  **References**
        !!
        !!  - Y. Saad. "Iterative methods for sparse linear systems", SIAM 2nd edition, 2003.
        !!    see Chapter 6.3: Arnoldi's method.
        !!
        !!  ### Syntax
        !!
        !!  ```fortran
        !!      call arnoldi(A, X, H, info [, kstart] [, kend] [, tol] [, transpose] [, blksize])
        !!  ```
        !!
        !!  ### Arguments
        !!
        !!  `A` : Linear operator derived from one the base types provided by the `AbstractLinops`
        !!        module. The operator needs to be square, i.e. the dimension of its domain and
        !!        co-domain is the same. It is an `intent(in)` argument.
        !!
        !!  `X` : Array of types derived from one the base types provided by the `AbstractVectors`
        !!        module. It needs to be consistent with the type of `A`. On exit, it contains the
        !!        the computed Krylov vectors. The first entry `X(1)` is the starting vector for
        !!        the Arnoldi factorization. Additionally, the maximum number of Arnoldi steps
        !!        is equal to `size(X) - 1`. It is an `intent(inout)` argument.
        !!
        !!  `H` : `real` or `complex` rank-2 array. On exit, it contains the \( (k+1) \times k\)
        !!         upper Hessenberg matrix computed from the Arnoldi factorization. It is an
        !!         `intent(inout)` argument.
        !!
        !!  `info` : `integer` variable. It is the `LightKrylov` information flag. On exit, if
        !!           `info` > 0, the Arnoldi factorization experienced a lucky breakdown. 
        !!            The array of Krylov vectors `X` spans an \(A\)-invariant subpsace of
        !!            dimension `info`.
        !!
        !!  `kstart` (*optional*): `integer` value determining the index of the first Arnoldi
        !!                          step to be computed. By default, `kstart = 1`.
        !!
        !!  `kend` (*optional*): `integer` value determining the index of the last Arnoldi step
        !!                        to be computed. By default, `kend = size(X) - 1`.
        !!
        !!  `tol` (*optional*): Numerical tolerance below which a subspace is considered
        !!                      to be \( A \)-invariant. By default `tol = atol_sp` or
        !!                      `tol = atol_rp` depending on the kind of `A`.
        !!
        !!  `transpose` (*optional*): `logical` flag determining whether the Arnoldi factorization
        !!                             is applied to \( A \) or \( A^H \). Default `transpose = .false.`
        !!
        !!  `blksize` (*optional*): `integer` value determining the dimension of a block for the
        !!                          block Arnoldi factorization. Default is `blksize=1`.
        module procedure arnoldi_rsp
        module procedure arnoldi_rdp
        module procedure arnoldi_csp
        module procedure arnoldi_cdp
    end interface

    interface initialize_krylov_subspace
        !!  ### Description
        !!
        !!  Utility function to initialize a basis for a Krylov subspace.
        !!
        !!  ### Syntax
        !!
        !!  ```fortran
        !!      call initialize_krylov_subspace(X [, X0])
        !!  ```
        !!
        !!  ### Arguments
        !!
        !!  `X` : Array of vectors that needs to be initialized. It is an `intent(inout)`
        !!        argument. Note that the first action in the subroutine is
        !!        `call zero_basis(X)`, effectively zeroing-out any data stored.
        !!
        !!  `X0` (*optional*) : Collection of vectors which will form the first few
        !!                      Krylov vectors. Note that `X0` need not be an orthonormal
        !!                      basis as this subroutine includes a `call qr(X0)`.
        module procedure initialize_krylov_subspace_rsp
        module procedure initialize_krylov_subspace_rdp
        module procedure initialize_krylov_subspace_csp
        module procedure initialize_krylov_subspace_cdp
    end interface

    interface is_orthonormal
        !!  ### Description
        !!
        !!  Utility function returning a logical `.true.` if the set of vectors stored in \( X \) form
        !!  an orthonormal set of vectors and `.false.` otherwise.
        !!
        !!  ### Syntax
        !!
        !!  ```fortran
        !!      out = is_orthonormal(X)
        !!  ```
        !!
        !!  ### Arguments
        !!
        !!  `X` : Array of derived types extended from the base types provided in the
        !!        `AbstractVectors` module.
        module procedure is_orthonormal_rsp
        module procedure is_orthonormal_rdp
        module procedure is_orthonormal_csp
        module procedure is_orthonormal_cdp
    end interface

    interface orthonormalize_basis
        !!  ### Description
        !!
        !!  Given an array \( X \) of vectors, it computes an orthonormal basis for its
        !!  column-span using the `double_gram_schmidt` process. All computations are done
        !!  in-place.
        !!
        !!  ### Syntax
        !!
        !!  ```fortran
        !!      call orthonormalize_basis(X)
        !!  ```
        !!
        !!  ### Arguments
        !!
        !!  `X` : Array of `abstract_vector` to orthonormalize. Note that this process is done
        !!        in-place. It is an `intent(inout)` argument.
        module procedure orthonormalize_basis_rsp
        module procedure orthonormalize_basis_rdp
        module procedure orthonormalize_basis_csp
        module procedure orthonormalize_basis_cdp
    end interface

    interface orthogonalize_against_basis
        module procedure orthogonalize_vector_against_basis_rsp
        module procedure orthogonalize_basis_against_basis_rsp
        module procedure orthogonalize_vector_against_basis_rdp
        module procedure orthogonalize_basis_against_basis_rdp
        module procedure orthogonalize_vector_against_basis_csp
        module procedure orthogonalize_basis_against_basis_csp
        module procedure orthogonalize_vector_against_basis_cdp
        module procedure orthogonalize_basis_against_basis_cdp
    end interface

    interface double_gram_schmidt_step
        !!  ### Description
        !!
        !!  Given an array \( X \) of `abstract_vector` and an `abstract_vector` 
        !!  (or array of `abstract_vectors`) \( y \), this subroutine returns a modified 
        !!  vector \( y \) orthogonal to all columns of \( X \), i.e.
        !!
        !!  \[
        !!      y \leftarrow \left( I - XX^H \right) y,
        !!  \]
        !!
        !!  using a double Gram-Schmidt process. On exit, \( y \) is orthogonal to \( X \) albeit
        !!  does not have unit-norm. Note moreover that \( X \) is assumed to be an orthonormal 
        !!  set of vectors. The function can also return the projection coefficients 
        !!  \( \beta = X^H y \).
        !!
        !!  ### Syntax
        !!
        !!  ```fortran
        !!      call double_gram_schmidt_step(y, X, info [, if_chk_orthonormal] [, beta])
        !!  ```
        !!
        !!  ### Arguments
        !!
        !!  `y` : `abstract_vector` (or array of `abstract_vector`) that needs to be
        !!        orthogonalize **in-place** against \( X \).
        !!
        !!  `X` : Array of `abstract_vector` against which \( y \) needs to be orthogonalized.
        !!        Note the function assumes that \( X \) is an orthonormal set of vectors, i.e.
        !!        \( X^H X = I \). If it this is not the case, the result are meaningless.
        !!
        !!  `info` : `integer` Information flag.
        !!
        !!  `if_chk_orthonormal` (*optional*) : `logical` flag (default `.true.`) to check
        !!      whether \( X \) is an orthonormal set of vectors or not. If the orthonormality
        !!      returns `.false.`, the function throws an error. Note that this check is however
        !!      computationally expensive and can be disable for the sake of performances.
        !!
        !!  `beta` (*optional*) : `real` or `complex` array containing the coefficients
        !!      \( \beta = X^H y \).
        module procedure DGS_vector_against_basis_rsp
        module procedure DGS_basis_against_basis_rsp
        module procedure DGS_vector_against_basis_rdp
        module procedure DGS_basis_against_basis_rdp
        module procedure DGS_vector_against_basis_csp
        module procedure DGS_basis_against_basis_csp
        module procedure DGS_vector_against_basis_cdp
        module procedure DGS_basis_against_basis_cdp
    end interface

    interface lanczos
        !!  ### Description
        !!
        !!  Given a symmetric or Hermitian linear operator \( A \), find matrices \( X \) and
        !!  \( T \) such that
        !!
        !!  \[
        !!      AX_k = X_k T_k + t_{k+1, k} x_{k+1} e_k^T,
        !!  \]
        !!
        !!  where \( X \) is an orthogonal basis and \( T \) is symmetric tridiagonal.
        !!
        !!  **Algorithmic Features**
        !!
        !!  - The operator \( A \) only needs to be accessed through matrix-vector products.
        !!  - Constructs an orthonormal Krylov basis \( X \) via the Lanczos process with full
        !!    reorthogonalization.
        !!  - Constructs a symmetric tridiagonal matrix \( T \) whose eigenvalues approximates those of \( A \).
        !!  - Checks for convergence and invariant subspaces.
        !!
        !!  **References**
        !!
        !!  - Y. Saad. "Iterative methods for sparse linear systems", SIAM 2nd edition, 2003.
        !!    see Chapter 6.6: The symmetric Lanczos algorithm.
        !!
        !!  ### Syntax
        !!
        !!  ```fortran
        !!      call lanczos(A, X, T, info [, kstart] [, kend] [, tol])
        !!  ```
        !!
        !!  ### Arguments
        !!
        !!  `A` : Symmetric or Hermitian linear operator derived from one the base types 
        !!        provided by the `AbstractLinops` module. It is an `intent(in)` argument.
        !!
        !!  `X` : Array of types derived from one the base types provided by the `AbstractVectors`
        !!        module. It needs to be consistent with the type of `A`. On exit, it contains the
        !!        the computed Krylov vectors. The first entry `X(1)` is the starting vector for
        !!        the Lanczos factorization. Additionally, the maximum number of Lanczos steps
        !!        is equal to `size(X) - 1`. It is an `intent(inout)` argument.
        !!
        !!  `T` : `real` or `complex` rank-2 array. On exit, it contains the \( (k+1) \times k\)
        !!         symmetric tridiagonal matrix computed from the Arnoldi factorization. It is an
        !!         `intent(inout)` argument.
        !!
        !!  `info` : `integer` variable. It is the `LightKrylov` information flag. On exit, if
        !!           `info` > 0, the Lanczos factorization experienced a lucky breakdown. 
        !!            The array of Krylov vectors `X` spans an \(A\)-invariant subpsace of
        !!            dimension `info`.
        !!
        !!  `kstart` (*optional*): `integer` value determining the index of the first Lanczos
        !!                          step to be computed. By default, `kstart = 1`.
        !!
        !!  `kend` (*optional*): `integer` value determining the index of the last Lanczos step
        !!                        to be computed. By default, `kend = size(X) - 1`.
        !!
        !!  `tol` (*optional*): Numerical tolerance below which a subspace is considered
        !!                      to be \( A \)-invariant. By default `tol = atol_sp` or
        !!                      `tol = atol_rp` depending on the kind of `A`.
        module procedure lanczos_tridiagonalization_rsp
        module procedure lanczos_tridiagonalization_rdp
        module procedure lanczos_tridiagonalization_csp
        module procedure lanczos_tridiagonalization_cdp
    end interface

    interface bidiagonalization
        !!  ### Description
        !!
        !!  Given a general linear operator \( A \), find matrices \( U \), \( V \) and
        !!  \( B \) such that
        !!
        !!  \[
        !!      \begin{aligned}
        !!          AV_k & = U_{k+1} B_k, \\
        !!          A^H U_{k+1} & = V_k B_k^T + b_{k+1} v_{k+1} e_{k+1}^T
        !!      \end{aligned}
        !!  \]
        !!
        !!  where \( U \) and \( V \) are orthogonal bases for the column span and row span
        !!  of \( A \), respectively, and \( B \) is a bidiagonal matrix.
        !!
        !!  **Algorithmic Features**
        !!
        !!  - The operator \( A \) only needs to be accessed through matrix-vector products.
        !!  - Constructs an orthonormal Krylov basis \( U \) for the column span of \( A \).
        !!  - Constructs an orthonormal Krylov basis \( V \) for the row span of \( A \).
        !!  - Constructs a bidiagonal matrix \( B \) whose singular values approximates 
        !!    those of \( A \).
        !!  - Checks for convergence and invariant subspaces.
        !!
        !!  **References**
        !!
        !!  - R. M. Larsen. "Lanczos bidiagonalization with partial reorthogonalization." 
        !!    Technical Report, 1998. [(PDF)](http://sun.stanford.edu/~rmunk/PROPACK/paper.pdf)
        !!
        !!  ### Syntax
        !!
        !!  ```fortran
        !!      call bidiagonalization(A, U, V, B, info [, kstart] [, kend] [, tol])
        !!  ```
        !!
        !!  ### Arguments
        !!
        !!  `A` : Linear operator derived from one the base types provided by the 
        !!        `AbstractLinops` module. It is an `intent(in)` argument.
        !!
        !!  `U` : Array of types derived from one the base types provided by the `AbstractVectors`
        !!        module. It needs to be consistent with the type of `A`. On exit, it contains the
        !!        the computed Krylov vectors for the column span of `A`. The first entry `U(1)` 
        !!        is the starting vector for the Lanczos factorization. Additionally, the 
        !!        maximum number of Lanczos steps is equal to `size(X) - 1`. 
        !!        It is an `intent(inout)` argument.
        !!
        !!  `V` : Array of types derived from one the base types provided by the `AbstractVectors`
        !!        module. It needs to be consistent with the type of `A`. On exit, it contains the
        !!        the computed Krylov vectors for the row span of `A`. It is an `intent(inout)` 
        !!        argument.
        !!
        !!  `B` : `real` or `complex` rank-2 array. On exit, it contains the \( (k+1) \times k\)
        !!         bidiagonal matrix computed from the Lanczos factorization. It is an
        !!         `intent(inout)` argument.
        !!
        !!  `info` : `integer` variable. It is the `LightKrylov` information flag. On exit, if
        !!           `info` > 0, the Lanczos factorization experienced a lucky breakdown. 
        !!
        !!  `kstart` (*optional*): `integer` value determining the index of the first Lanczos
        !!                          step to be computed. By default, `kstart = 1`.
        !!
        !!  `kend` (*optional*): `integer` value determining the index of the last Lanczos step
        !!                        to be computed. By default, `kend = size(X) - 1`.
        !!
        !!  `tol` (*optional*): Numerical tolerance below which a subspace is considered
        !!                      to be \( A \)-invariant. By default `tol = atol_sp` or
        !!                      `tol = atol_rp` depending on the kind of `A`.
        module procedure lanczos_bidiagonalization_rsp
        module procedure lanczos_bidiagonalization_rdp
        module procedure lanczos_bidiagonalization_csp
        module procedure lanczos_bidiagonalization_cdp
    end interface

    interface krylov_schur
        !!  ### Description
        !!
        !!  Given a partial Krylov decomposition
        !!
        !!  \[
        !!      AX_k = X_k H_k + h_{k+1, k} x_{k+1} e_k^H,
        !!  \]
        !!
        !!  this subroutine implements the Krylov-Schur restarting strategy proposed by
        !!  Stewart [1].
        !!
        !!  **References**
        !!
        !!  - G. W. Stewart. "A Krylov-Schur algorithm for large eigenproblems".
        !!    SIAM Journal on Matrix Analysis and Applications, vol 23 (3), 2002.
        !!
        !!  ### Syntax
        !!
        !!  ```fortran
        !!      call krylov_schur(n, X, H, select_eigs)
        !!  ```
        !!
        !!  ### Arguments
        !!
        !!  `n` : Number of selected eigenvalues moved to the upper left-block of the 
        !!        Schur matrix. It is an `intent(out)` argument.
        !!
        !!  `X` : On entry, array of `abstract_vector` computed using the Arnoldi process.
        !!        On exit, the first `n` columns form an orthonormal basis for the eigenspace
        !!        associated with eigenvalues moved to the upper left-block of the Schur matrix.
        !!        It is an `intent(inout)` argument.
        !!
        !!  `H` : On entry, `real` of `complex` upper Hessenberg matrix computed using the
        !!        Arnoldi process. On exit, the leading \( n \times n\) block contains the
        !!        \( S_{11} \) block of the re-ordered Schur matrix containing the selected
        !!        eigenvalues. It is an `intent(inout)` argument.
        !!
        !!  `select_eigs` : Procedure to select which eigenvalues to move in the upper-left
        !!                  block. It is an `intent(inout)` argument.
        module procedure krylov_schur_rsp
        module procedure krylov_schur_rdp
        module procedure krylov_schur_csp
        module procedure krylov_schur_cdp
    end interface

    !----------------------------------------------------------
    !-----     ABSTRACT EIGENVALUE SELECTOR INTERFACE     -----
    !----------------------------------------------------------

    abstract interface
        function eigvals_select_sp(lambda) result(selected)
            import sp
            complex(sp), intent(in) :: lambda(:)
            logical                       :: selected(size(lambda))
        end function eigvals_select_sp
        function eigvals_select_dp(lambda) result(selected)
            import dp
            complex(dp), intent(in) :: lambda(:)
            logical                       :: selected(size(lambda))
        end function eigvals_select_dp
    end interface

contains

    !-------------------------------------
    !-----     UTILITY FUNCTIONS     -----
    !-------------------------------------

    logical function is_orthonormal_rsp(X) result(ortho)
        class(abstract_vector_rsp), intent(in) :: X(:)
        real(sp), dimension(size(X), size(X)) :: G
        ortho = .true.
        call innerprod(G, X, X)
        if (norm2(abs(G - eye(size(X)))) > rtol_sp) then
            ! The basis is not orthonormal. Cannot orthonormalize.
            ortho = .false.
        end if
    end function
    logical function is_orthonormal_rdp(X) result(ortho)
        class(abstract_vector_rdp), intent(in) :: X(:)
        real(dp), dimension(size(X), size(X)) :: G
        ortho = .true.
        call innerprod(G, X, X)
        if (norm2(abs(G - eye(size(X)))) > rtol_sp) then
            ! The basis is not orthonormal. Cannot orthonormalize.
            ortho = .false.
        end if
    end function
    logical function is_orthonormal_csp(X) result(ortho)
        class(abstract_vector_csp), intent(in) :: X(:)
        complex(sp), dimension(size(X), size(X)) :: G
        ortho = .true.
        call innerprod(G, X, X)
        if (norm2(abs(G - eye(size(X)))) > rtol_sp) then
            ! The basis is not orthonormal. Cannot orthonormalize.
            ortho = .false.
        end if
    end function
    logical function is_orthonormal_cdp(X) result(ortho)
        class(abstract_vector_cdp), intent(in) :: X(:)
        complex(dp), dimension(size(X), size(X)) :: G
        ortho = .true.
        call innerprod(G, X, X)
        if (norm2(abs(G - eye(size(X)))) > rtol_sp) then
            ! The basis is not orthonormal. Cannot orthonormalize.
            ortho = .false.
        end if
    end function

    subroutine initialize_krylov_subspace_rsp(X, X0)
        class(abstract_vector_rsp), intent(inout) :: X(:)
        class(abstract_vector_rsp), optional, intent(in) :: X0(:)

        ! Internal variables.
        integer :: p

        ! Zero-out X.
        call zero_basis(X)

        ! Deals with optional args.
        if(present(X0)) then
            p = size(X0)
            ! Initialize.
            call copy(X(:p), X0)
            ! Orthonormalize.
            call orthonormalize_basis(X(:p))
        endif

        return
    end subroutine initialize_krylov_subspace_rsp

    subroutine orthonormalize_basis_rsp(X)
      !! Orthonormalizes the `abstract_vector` basis `X`
      class(abstract_vector_rsp), intent(inout) :: X(:)
      !! Input `abstract_vector` basis to orthogonalize against
      
      ! internals
      real(sp) :: R(size(X),size(X))
      integer :: info

      ! internals
      call qr(X, R, info)
      call check_info(info, 'qr', module=this_module, procedure='orthonormalize_basis_rsp')
      
      return
    end subroutine orthonormalize_basis_rsp

    subroutine orthogonalize_vector_against_basis_rsp(y, X, info, if_chk_orthonormal, beta)
      !! Orthogonalizes the `abstract_vector` `y` against a basis `X` of `abstract_vector`.
      class(abstract_vector_rsp), intent(inout) :: y
      !! Input `abstract_vector` to orthogonalize
      class(abstract_vector_rsp), intent(in)    :: X(:)
      !! Input `abstract_vector` basis to orthogonalize against
      integer, intent(out) :: info
      !! Information flag.
      logical,                          optional, intent(in)    :: if_chk_orthonormal
      logical                                                   :: chk_X_orthonormality
      !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
      real(sp),                         optional, intent(out)   :: beta(:)
      !! Projection coefficients if requested

      ! internals
      real(sp) :: proj_coefficients(size(X))

      info = 0

      ! optional input argument
      chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

      ! check for zero vector
      if (y%norm() < atol_sp) info = 1

      if (chk_X_orthonormality) then
         block 
            real(sp), dimension(size(X), size(X)) :: G
            call innerprod(G, X, X)
            if (abs(G(size(X),size(X))) < rtol_sp) then
               ! The last vector in X is zero, it does not impact orthogonalisation
               info = -2
            else if (norm2(abs(G - eye(size(X)))) > rtol_sp) then
               ! The basis is not orthonormal. Cannot orthonormalize.
               info = -1
               return
            end if
         end block
      end if

      ! orthogonalize
      call innerprod(proj_coefficients, X, y)
      block
         class(abstract_vector_rsp), allocatable :: proj
         call linear_combination(proj, X, proj_coefficients)
         call y%sub(proj)
      end block

      if (present(beta)) then
         ! check size
         call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, 'orthogonalize_vector_against_basis_rsp')
         beta = proj_coefficients
      end if
      
      return
    end subroutine orthogonalize_vector_against_basis_rsp

    subroutine orthogonalize_basis_against_basis_rsp(Y, X, info, if_chk_orthonormal, beta)
      !! Orthogonalizes the `abstract_vector` basis `Y` against a basis `X` of `abstract_vector`.
      class(abstract_vector_rsp), intent(inout) :: Y(:)
      !! Input `abstract_vector` basis to orthogonalize
      class(abstract_vector_rsp), intent(in)    :: X(:)
      !! Input `abstract_vector` basis to orthogonalize against
      integer, intent(out) :: info
      !! Information flag.
      logical,                          optional, intent(in)    :: if_chk_orthonormal
      logical                                                   :: chk_X_orthonormality
      !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
      real(sp),                         optional, intent(out)   :: beta(:,:)
      !! Projection coefficients if requested

      ! internals
      real(sp) :: proj_coefficients(size(X), size(Y))
      integer :: i

      info = 0

      ! optional input argument
      chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

      ! check for zero vector
      do i = 1, size(Y)
         if ( Y(i)%norm() < atol_sp) info = i
      end do

      if (chk_X_orthonormality) then
         block 
            real(sp), dimension(size(X), size(X)) :: G
            call innerprod(G, X, X)
            if (abs(G(size(X),size(X))) < rtol_sp) then
               ! The last vector in X is zero, it does not impact orthogonalisation
               info = -2
            else if (norm2(abs(G - eye(size(X)))) > rtol_sp) then
               ! The basis is not orthonormal. Cannot orthonormalize.
               info = -1
               return
            end if
         end block
      end if
      
      ! orthogonalize
      call innerprod(proj_coefficients, X, Y)
      block
         class(abstract_vector_rsp), allocatable :: proj(:)
         call linear_combination(proj, X, proj_coefficients)
         call axpby_basis(Y, one_rsp, proj, -one_rsp)
      end block

      if (present(beta)) then
         ! check size
         call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, 'orthogonalize_basis_against_basis_rsp')
         beta = proj_coefficients
      end if
      
      return
    end subroutine orthogonalize_basis_against_basis_rsp

    subroutine DGS_vector_against_basis_rsp(y, X, info, if_chk_orthonormal, beta)
      !! Computes one step of the double Gram-Schmidt orthogonalization process of the
      !! `abstract_vector` `y` against the `abstract_vector` basis `X`
      class(abstract_vector_rsp), intent(inout) :: y
      !! Input `abstract_vector` to orthogonalize
      class(abstract_vector_rsp), intent(in)    :: X(:)
      !! Input `abstract_vector` basis to orthogonalize against
      integer, intent(out) :: info
      !! Information flag.
      logical,                          optional, intent(in)    :: if_chk_orthonormal
      logical                                                   :: chk_X_orthonormality
      !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
      real(sp),                         optional, intent(out)   :: beta(:)
      !! Projection coefficients if requested

      ! internals
      real(sp), dimension(size(X)) :: proj_coefficients, wrk

      ! optional input argument
      chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

      info = 0

      proj_coefficients = zero_rsp; wrk = zero_rsp

      ! Orthogonalize vector y w.r.t. to Krylov basis X in two passes of GS.
      ! first pass
      call orthogonalize_against_basis(y, X, info, if_chk_orthonormal=.false., beta=proj_coefficients)
      call check_info(info, 'orthogonalize_against_basis_p1', module=this_module, &
                        & procedure='DGS_vector_against_basis_rsp, pass 1')
      ! second pass
      call orthogonalize_against_basis(y, X, info, if_chk_orthonormal=.false., beta=wrk)
      call check_info(info, 'orthogonalize_against_basis_p2', module=this_module, &
                        & procedure='DGS_vector_against_basis_rsp, pass 2')
      ! combine passes
      proj_coefficients = proj_coefficients + wrk

      if (present(beta)) then
         ! check size
         call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, 'DGS_vector_against_basis_rsp')
         beta = proj_coefficients
      end if

    end subroutine DGS_vector_against_basis_rsp

    subroutine DGS_basis_against_basis_rsp(y, X, info, if_chk_orthonormal, beta)
      !! Computes one step of the double Gram-Schmidt orthogonalization process of the
      !! `abstract_vector` `y` against the `abstract_vector` basis `X`
      class(abstract_vector_rsp), intent(inout) :: Y(:)
      !! Input `abstract_vector` basis to orthogonalize
      class(abstract_vector_rsp), intent(in)    :: X(:)
      !! Input `abstract_vector` basis to orthogonalize against
      integer, intent(out) :: info
      !! Information flag.
      logical,                          optional, intent(in)    :: if_chk_orthonormal
      logical                                                   :: chk_X_orthonormality
      !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
      real(sp),                         optional, intent(out)   :: beta(:,:)
      !! Projection coefficients if requested

      ! internals
      real(sp), dimension(size(X),size(Y)) :: proj_coefficients, wrk

      ! optional input argument
      chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

      info = 0

      proj_coefficients = zero_rsp; wrk = zero_rsp

      ! Orthogonalize Krylov basis Y w.r.t. to Krylov basis X in two passes of GS.
      ! first pass
      call orthogonalize_against_basis(Y, X, info, if_chk_orthonormal=.false., beta=proj_coefficients)
      call check_info(info, 'orthogonalize_against_basis_p1', module=this_module, procedure='DGS_basis_against_basis_rsp, first&
          & pass')
      ! second pass
      call orthogonalize_against_basis(Y, X, info, if_chk_orthonormal=.false., beta=wrk)
      call check_info(info, 'orthogonalize_against_basis_p2', module=this_module, procedure='DGS_basis_against_basis_rsp, second&
          & pass')
      ! combine passes
      proj_coefficients = proj_coefficients + wrk

      if (present(beta)) then
         ! check size
         call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, 'DGS_basis_against_basis_rsp')
         beta = proj_coefficients
      end if

    end subroutine DGS_basis_against_basis_rsp  

    subroutine initialize_krylov_subspace_rdp(X, X0)
        class(abstract_vector_rdp), intent(inout) :: X(:)
        class(abstract_vector_rdp), optional, intent(in) :: X0(:)

        ! Internal variables.
        integer :: p

        ! Zero-out X.
        call zero_basis(X)

        ! Deals with optional args.
        if(present(X0)) then
            p = size(X0)
            ! Initialize.
            call copy(X(:p), X0)
            ! Orthonormalize.
            call orthonormalize_basis(X(:p))
        endif

        return
    end subroutine initialize_krylov_subspace_rdp

    subroutine orthonormalize_basis_rdp(X)
      !! Orthonormalizes the `abstract_vector` basis `X`
      class(abstract_vector_rdp), intent(inout) :: X(:)
      !! Input `abstract_vector` basis to orthogonalize against
      
      ! internals
      real(dp) :: R(size(X),size(X))
      integer :: info

      ! internals
      call qr(X, R, info)
      call check_info(info, 'qr', module=this_module, procedure='orthonormalize_basis_rdp')
      
      return
    end subroutine orthonormalize_basis_rdp

    subroutine orthogonalize_vector_against_basis_rdp(y, X, info, if_chk_orthonormal, beta)
      !! Orthogonalizes the `abstract_vector` `y` against a basis `X` of `abstract_vector`.
      class(abstract_vector_rdp), intent(inout) :: y
      !! Input `abstract_vector` to orthogonalize
      class(abstract_vector_rdp), intent(in)    :: X(:)
      !! Input `abstract_vector` basis to orthogonalize against
      integer, intent(out) :: info
      !! Information flag.
      logical,                          optional, intent(in)    :: if_chk_orthonormal
      logical                                                   :: chk_X_orthonormality
      !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
      real(dp),                         optional, intent(out)   :: beta(:)
      !! Projection coefficients if requested

      ! internals
      real(dp) :: proj_coefficients(size(X))

      info = 0

      ! optional input argument
      chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

      ! check for zero vector
      if (y%norm() < atol_dp) info = 1

      if (chk_X_orthonormality) then
         block 
            real(dp), dimension(size(X), size(X)) :: G
            call innerprod(G, X, X)
            if (abs(G(size(X),size(X))) < rtol_dp) then
               ! The last vector in X is zero, it does not impact orthogonalisation
               info = -2
            else if (norm2(abs(G - eye(size(X)))) > rtol_dp) then
               ! The basis is not orthonormal. Cannot orthonormalize.
               info = -1
               return
            end if
         end block
      end if

      ! orthogonalize
      call innerprod(proj_coefficients, X, y)
      block
         class(abstract_vector_rdp), allocatable :: proj
         call linear_combination(proj, X, proj_coefficients)
         call y%sub(proj)
      end block

      if (present(beta)) then
         ! check size
         call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, 'orthogonalize_vector_against_basis_rdp')
         beta = proj_coefficients
      end if
      
      return
    end subroutine orthogonalize_vector_against_basis_rdp

    subroutine orthogonalize_basis_against_basis_rdp(Y, X, info, if_chk_orthonormal, beta)
      !! Orthogonalizes the `abstract_vector` basis `Y` against a basis `X` of `abstract_vector`.
      class(abstract_vector_rdp), intent(inout) :: Y(:)
      !! Input `abstract_vector` basis to orthogonalize
      class(abstract_vector_rdp), intent(in)    :: X(:)
      !! Input `abstract_vector` basis to orthogonalize against
      integer, intent(out) :: info
      !! Information flag.
      logical,                          optional, intent(in)    :: if_chk_orthonormal
      logical                                                   :: chk_X_orthonormality
      !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
      real(dp),                         optional, intent(out)   :: beta(:,:)
      !! Projection coefficients if requested

      ! internals
      real(dp) :: proj_coefficients(size(X), size(Y))
      integer :: i

      info = 0

      ! optional input argument
      chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

      ! check for zero vector
      do i = 1, size(Y)
         if ( Y(i)%norm() < atol_dp) info = i
      end do

      if (chk_X_orthonormality) then
         block 
            real(dp), dimension(size(X), size(X)) :: G
            call innerprod(G, X, X)
            if (abs(G(size(X),size(X))) < rtol_dp) then
               ! The last vector in X is zero, it does not impact orthogonalisation
               info = -2
            else if (norm2(abs(G - eye(size(X)))) > rtol_dp) then
               ! The basis is not orthonormal. Cannot orthonormalize.
               info = -1
               return
            end if
         end block
      end if
      
      ! orthogonalize
      call innerprod(proj_coefficients, X, Y)
      block
         class(abstract_vector_rdp), allocatable :: proj(:)
         call linear_combination(proj, X, proj_coefficients)
         call axpby_basis(Y, one_rdp, proj, -one_rdp)
      end block

      if (present(beta)) then
         ! check size
         call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, 'orthogonalize_basis_against_basis_rdp')
         beta = proj_coefficients
      end if
      
      return
    end subroutine orthogonalize_basis_against_basis_rdp

    subroutine DGS_vector_against_basis_rdp(y, X, info, if_chk_orthonormal, beta)
      !! Computes one step of the double Gram-Schmidt orthogonalization process of the
      !! `abstract_vector` `y` against the `abstract_vector` basis `X`
      class(abstract_vector_rdp), intent(inout) :: y
      !! Input `abstract_vector` to orthogonalize
      class(abstract_vector_rdp), intent(in)    :: X(:)
      !! Input `abstract_vector` basis to orthogonalize against
      integer, intent(out) :: info
      !! Information flag.
      logical,                          optional, intent(in)    :: if_chk_orthonormal
      logical                                                   :: chk_X_orthonormality
      !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
      real(dp),                         optional, intent(out)   :: beta(:)
      !! Projection coefficients if requested

      ! internals
      real(dp), dimension(size(X)) :: proj_coefficients, wrk

      ! optional input argument
      chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

      info = 0

      proj_coefficients = zero_rdp; wrk = zero_rdp

      ! Orthogonalize vector y w.r.t. to Krylov basis X in two passes of GS.
      ! first pass
      call orthogonalize_against_basis(y, X, info, if_chk_orthonormal=.false., beta=proj_coefficients)
      call check_info(info, 'orthogonalize_against_basis_p1', module=this_module, &
                        & procedure='DGS_vector_against_basis_rdp, pass 1')
      ! second pass
      call orthogonalize_against_basis(y, X, info, if_chk_orthonormal=.false., beta=wrk)
      call check_info(info, 'orthogonalize_against_basis_p2', module=this_module, &
                        & procedure='DGS_vector_against_basis_rdp, pass 2')
      ! combine passes
      proj_coefficients = proj_coefficients + wrk

      if (present(beta)) then
         ! check size
         call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, 'DGS_vector_against_basis_rdp')
         beta = proj_coefficients
      end if

    end subroutine DGS_vector_against_basis_rdp

    subroutine DGS_basis_against_basis_rdp(y, X, info, if_chk_orthonormal, beta)
      !! Computes one step of the double Gram-Schmidt orthogonalization process of the
      !! `abstract_vector` `y` against the `abstract_vector` basis `X`
      class(abstract_vector_rdp), intent(inout) :: Y(:)
      !! Input `abstract_vector` basis to orthogonalize
      class(abstract_vector_rdp), intent(in)    :: X(:)
      !! Input `abstract_vector` basis to orthogonalize against
      integer, intent(out) :: info
      !! Information flag.
      logical,                          optional, intent(in)    :: if_chk_orthonormal
      logical                                                   :: chk_X_orthonormality
      !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
      real(dp),                         optional, intent(out)   :: beta(:,:)
      !! Projection coefficients if requested

      ! internals
      real(dp), dimension(size(X),size(Y)) :: proj_coefficients, wrk

      ! optional input argument
      chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

      info = 0

      proj_coefficients = zero_rdp; wrk = zero_rdp

      ! Orthogonalize Krylov basis Y w.r.t. to Krylov basis X in two passes of GS.
      ! first pass
      call orthogonalize_against_basis(Y, X, info, if_chk_orthonormal=.false., beta=proj_coefficients)
      call check_info(info, 'orthogonalize_against_basis_p1', module=this_module, procedure='DGS_basis_against_basis_rdp, first&
          & pass')
      ! second pass
      call orthogonalize_against_basis(Y, X, info, if_chk_orthonormal=.false., beta=wrk)
      call check_info(info, 'orthogonalize_against_basis_p2', module=this_module, procedure='DGS_basis_against_basis_rdp, second&
          & pass')
      ! combine passes
      proj_coefficients = proj_coefficients + wrk

      if (present(beta)) then
         ! check size
         call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, 'DGS_basis_against_basis_rdp')
         beta = proj_coefficients
      end if

    end subroutine DGS_basis_against_basis_rdp  

    subroutine initialize_krylov_subspace_csp(X, X0)
        class(abstract_vector_csp), intent(inout) :: X(:)
        class(abstract_vector_csp), optional, intent(in) :: X0(:)

        ! Internal variables.
        integer :: p

        ! Zero-out X.
        call zero_basis(X)

        ! Deals with optional args.
        if(present(X0)) then
            p = size(X0)
            ! Initialize.
            call copy(X(:p), X0)
            ! Orthonormalize.
            call orthonormalize_basis(X(:p))
        endif

        return
    end subroutine initialize_krylov_subspace_csp

    subroutine orthonormalize_basis_csp(X)
      !! Orthonormalizes the `abstract_vector` basis `X`
      class(abstract_vector_csp), intent(inout) :: X(:)
      !! Input `abstract_vector` basis to orthogonalize against
      
      ! internals
      complex(sp) :: R(size(X),size(X))
      integer :: info

      ! internals
      call qr(X, R, info)
      call check_info(info, 'qr', module=this_module, procedure='orthonormalize_basis_csp')
      
      return
    end subroutine orthonormalize_basis_csp

    subroutine orthogonalize_vector_against_basis_csp(y, X, info, if_chk_orthonormal, beta)
      !! Orthogonalizes the `abstract_vector` `y` against a basis `X` of `abstract_vector`.
      class(abstract_vector_csp), intent(inout) :: y
      !! Input `abstract_vector` to orthogonalize
      class(abstract_vector_csp), intent(in)    :: X(:)
      !! Input `abstract_vector` basis to orthogonalize against
      integer, intent(out) :: info
      !! Information flag.
      logical,                          optional, intent(in)    :: if_chk_orthonormal
      logical                                                   :: chk_X_orthonormality
      !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
      complex(sp),                         optional, intent(out)   :: beta(:)
      !! Projection coefficients if requested

      ! internals
      complex(sp) :: proj_coefficients(size(X))

      info = 0

      ! optional input argument
      chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

      ! check for zero vector
      if (y%norm() < atol_sp) info = 1

      if (chk_X_orthonormality) then
         block 
            complex(sp), dimension(size(X), size(X)) :: G
            call innerprod(G, X, X)
            if (abs(G(size(X),size(X))) < rtol_sp) then
               ! The last vector in X is zero, it does not impact orthogonalisation
               info = -2
            else if (norm2(abs(G - eye(size(X)))) > rtol_sp) then
               ! The basis is not orthonormal. Cannot orthonormalize.
               info = -1
               return
            end if
         end block
      end if

      ! orthogonalize
      call innerprod(proj_coefficients, X, y)
      block
         class(abstract_vector_csp), allocatable :: proj
         call linear_combination(proj, X, proj_coefficients)
         call y%sub(proj)
      end block

      if (present(beta)) then
         ! check size
         call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, 'orthogonalize_vector_against_basis_csp')
         beta = proj_coefficients
      end if
      
      return
    end subroutine orthogonalize_vector_against_basis_csp

    subroutine orthogonalize_basis_against_basis_csp(Y, X, info, if_chk_orthonormal, beta)
      !! Orthogonalizes the `abstract_vector` basis `Y` against a basis `X` of `abstract_vector`.
      class(abstract_vector_csp), intent(inout) :: Y(:)
      !! Input `abstract_vector` basis to orthogonalize
      class(abstract_vector_csp), intent(in)    :: X(:)
      !! Input `abstract_vector` basis to orthogonalize against
      integer, intent(out) :: info
      !! Information flag.
      logical,                          optional, intent(in)    :: if_chk_orthonormal
      logical                                                   :: chk_X_orthonormality
      !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
      complex(sp),                         optional, intent(out)   :: beta(:,:)
      !! Projection coefficients if requested

      ! internals
      complex(sp) :: proj_coefficients(size(X), size(Y))
      integer :: i

      info = 0

      ! optional input argument
      chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

      ! check for zero vector
      do i = 1, size(Y)
         if ( Y(i)%norm() < atol_sp) info = i
      end do

      if (chk_X_orthonormality) then
         block 
            complex(sp), dimension(size(X), size(X)) :: G
            call innerprod(G, X, X)
            if (abs(G(size(X),size(X))) < rtol_sp) then
               ! The last vector in X is zero, it does not impact orthogonalisation
               info = -2
            else if (norm2(abs(G - eye(size(X)))) > rtol_sp) then
               ! The basis is not orthonormal. Cannot orthonormalize.
               info = -1
               return
            end if
         end block
      end if
      
      ! orthogonalize
      call innerprod(proj_coefficients, X, Y)
      block
         class(abstract_vector_csp), allocatable :: proj(:)
         call linear_combination(proj, X, proj_coefficients)
         call axpby_basis(Y, one_csp, proj, -one_csp)
      end block

      if (present(beta)) then
         ! check size
         call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, 'orthogonalize_basis_against_basis_csp')
         beta = proj_coefficients
      end if
      
      return
    end subroutine orthogonalize_basis_against_basis_csp

    subroutine DGS_vector_against_basis_csp(y, X, info, if_chk_orthonormal, beta)
      !! Computes one step of the double Gram-Schmidt orthogonalization process of the
      !! `abstract_vector` `y` against the `abstract_vector` basis `X`
      class(abstract_vector_csp), intent(inout) :: y
      !! Input `abstract_vector` to orthogonalize
      class(abstract_vector_csp), intent(in)    :: X(:)
      !! Input `abstract_vector` basis to orthogonalize against
      integer, intent(out) :: info
      !! Information flag.
      logical,                          optional, intent(in)    :: if_chk_orthonormal
      logical                                                   :: chk_X_orthonormality
      !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
      complex(sp),                         optional, intent(out)   :: beta(:)
      !! Projection coefficients if requested

      ! internals
      complex(sp), dimension(size(X)) :: proj_coefficients, wrk

      ! optional input argument
      chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

      info = 0

      proj_coefficients = zero_csp; wrk = zero_csp

      ! Orthogonalize vector y w.r.t. to Krylov basis X in two passes of GS.
      ! first pass
      call orthogonalize_against_basis(y, X, info, if_chk_orthonormal=.false., beta=proj_coefficients)
      call check_info(info, 'orthogonalize_against_basis_p1', module=this_module, &
                        & procedure='DGS_vector_against_basis_csp, pass 1')
      ! second pass
      call orthogonalize_against_basis(y, X, info, if_chk_orthonormal=.false., beta=wrk)
      call check_info(info, 'orthogonalize_against_basis_p2', module=this_module, &
                        & procedure='DGS_vector_against_basis_csp, pass 2')
      ! combine passes
      proj_coefficients = proj_coefficients + wrk

      if (present(beta)) then
         ! check size
         call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, 'DGS_vector_against_basis_csp')
         beta = proj_coefficients
      end if

    end subroutine DGS_vector_against_basis_csp

    subroutine DGS_basis_against_basis_csp(y, X, info, if_chk_orthonormal, beta)
      !! Computes one step of the double Gram-Schmidt orthogonalization process of the
      !! `abstract_vector` `y` against the `abstract_vector` basis `X`
      class(abstract_vector_csp), intent(inout) :: Y(:)
      !! Input `abstract_vector` basis to orthogonalize
      class(abstract_vector_csp), intent(in)    :: X(:)
      !! Input `abstract_vector` basis to orthogonalize against
      integer, intent(out) :: info
      !! Information flag.
      logical,                          optional, intent(in)    :: if_chk_orthonormal
      logical                                                   :: chk_X_orthonormality
      !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
      complex(sp),                         optional, intent(out)   :: beta(:,:)
      !! Projection coefficients if requested

      ! internals
      complex(sp), dimension(size(X),size(Y)) :: proj_coefficients, wrk

      ! optional input argument
      chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

      info = 0

      proj_coefficients = zero_csp; wrk = zero_csp

      ! Orthogonalize Krylov basis Y w.r.t. to Krylov basis X in two passes of GS.
      ! first pass
      call orthogonalize_against_basis(Y, X, info, if_chk_orthonormal=.false., beta=proj_coefficients)
      call check_info(info, 'orthogonalize_against_basis_p1', module=this_module, procedure='DGS_basis_against_basis_csp, first&
          & pass')
      ! second pass
      call orthogonalize_against_basis(Y, X, info, if_chk_orthonormal=.false., beta=wrk)
      call check_info(info, 'orthogonalize_against_basis_p2', module=this_module, procedure='DGS_basis_against_basis_csp, second&
          & pass')
      ! combine passes
      proj_coefficients = proj_coefficients + wrk

      if (present(beta)) then
         ! check size
         call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, 'DGS_basis_against_basis_csp')
         beta = proj_coefficients
      end if

    end subroutine DGS_basis_against_basis_csp  

    subroutine initialize_krylov_subspace_cdp(X, X0)
        class(abstract_vector_cdp), intent(inout) :: X(:)
        class(abstract_vector_cdp), optional, intent(in) :: X0(:)

        ! Internal variables.
        integer :: p

        ! Zero-out X.
        call zero_basis(X)

        ! Deals with optional args.
        if(present(X0)) then
            p = size(X0)
            ! Initialize.
            call copy(X(:p), X0)
            ! Orthonormalize.
            call orthonormalize_basis(X(:p))
        endif

        return
    end subroutine initialize_krylov_subspace_cdp

    subroutine orthonormalize_basis_cdp(X)
      !! Orthonormalizes the `abstract_vector` basis `X`
      class(abstract_vector_cdp), intent(inout) :: X(:)
      !! Input `abstract_vector` basis to orthogonalize against
      
      ! internals
      complex(dp) :: R(size(X),size(X))
      integer :: info

      ! internals
      call qr(X, R, info)
      call check_info(info, 'qr', module=this_module, procedure='orthonormalize_basis_cdp')
      
      return
    end subroutine orthonormalize_basis_cdp

    subroutine orthogonalize_vector_against_basis_cdp(y, X, info, if_chk_orthonormal, beta)
      !! Orthogonalizes the `abstract_vector` `y` against a basis `X` of `abstract_vector`.
      class(abstract_vector_cdp), intent(inout) :: y
      !! Input `abstract_vector` to orthogonalize
      class(abstract_vector_cdp), intent(in)    :: X(:)
      !! Input `abstract_vector` basis to orthogonalize against
      integer, intent(out) :: info
      !! Information flag.
      logical,                          optional, intent(in)    :: if_chk_orthonormal
      logical                                                   :: chk_X_orthonormality
      !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
      complex(dp),                         optional, intent(out)   :: beta(:)
      !! Projection coefficients if requested

      ! internals
      complex(dp) :: proj_coefficients(size(X))

      info = 0

      ! optional input argument
      chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

      ! check for zero vector
      if (y%norm() < atol_dp) info = 1

      if (chk_X_orthonormality) then
         block 
            complex(dp), dimension(size(X), size(X)) :: G
            call innerprod(G, X, X)
            if (abs(G(size(X),size(X))) < rtol_dp) then
               ! The last vector in X is zero, it does not impact orthogonalisation
               info = -2
            else if (norm2(abs(G - eye(size(X)))) > rtol_dp) then
               ! The basis is not orthonormal. Cannot orthonormalize.
               info = -1
               return
            end if
         end block
      end if

      ! orthogonalize
      call innerprod(proj_coefficients, X, y)
      block
         class(abstract_vector_cdp), allocatable :: proj
         call linear_combination(proj, X, proj_coefficients)
         call y%sub(proj)
      end block

      if (present(beta)) then
         ! check size
         call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, 'orthogonalize_vector_against_basis_cdp')
         beta = proj_coefficients
      end if
      
      return
    end subroutine orthogonalize_vector_against_basis_cdp

    subroutine orthogonalize_basis_against_basis_cdp(Y, X, info, if_chk_orthonormal, beta)
      !! Orthogonalizes the `abstract_vector` basis `Y` against a basis `X` of `abstract_vector`.
      class(abstract_vector_cdp), intent(inout) :: Y(:)
      !! Input `abstract_vector` basis to orthogonalize
      class(abstract_vector_cdp), intent(in)    :: X(:)
      !! Input `abstract_vector` basis to orthogonalize against
      integer, intent(out) :: info
      !! Information flag.
      logical,                          optional, intent(in)    :: if_chk_orthonormal
      logical                                                   :: chk_X_orthonormality
      !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
      complex(dp),                         optional, intent(out)   :: beta(:,:)
      !! Projection coefficients if requested

      ! internals
      complex(dp) :: proj_coefficients(size(X), size(Y))
      integer :: i

      info = 0

      ! optional input argument
      chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

      ! check for zero vector
      do i = 1, size(Y)
         if ( Y(i)%norm() < atol_dp) info = i
      end do

      if (chk_X_orthonormality) then
         block 
            complex(dp), dimension(size(X), size(X)) :: G
            call innerprod(G, X, X)
            if (abs(G(size(X),size(X))) < rtol_dp) then
               ! The last vector in X is zero, it does not impact orthogonalisation
               info = -2
            else if (norm2(abs(G - eye(size(X)))) > rtol_dp) then
               ! The basis is not orthonormal. Cannot orthonormalize.
               info = -1
               return
            end if
         end block
      end if
      
      ! orthogonalize
      call innerprod(proj_coefficients, X, Y)
      block
         class(abstract_vector_cdp), allocatable :: proj(:)
         call linear_combination(proj, X, proj_coefficients)
         call axpby_basis(Y, one_cdp, proj, -one_cdp)
      end block

      if (present(beta)) then
         ! check size
         call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, 'orthogonalize_basis_against_basis_cdp')
         beta = proj_coefficients
      end if
      
      return
    end subroutine orthogonalize_basis_against_basis_cdp

    subroutine DGS_vector_against_basis_cdp(y, X, info, if_chk_orthonormal, beta)
      !! Computes one step of the double Gram-Schmidt orthogonalization process of the
      !! `abstract_vector` `y` against the `abstract_vector` basis `X`
      class(abstract_vector_cdp), intent(inout) :: y
      !! Input `abstract_vector` to orthogonalize
      class(abstract_vector_cdp), intent(in)    :: X(:)
      !! Input `abstract_vector` basis to orthogonalize against
      integer, intent(out) :: info
      !! Information flag.
      logical,                          optional, intent(in)    :: if_chk_orthonormal
      logical                                                   :: chk_X_orthonormality
      !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
      complex(dp),                         optional, intent(out)   :: beta(:)
      !! Projection coefficients if requested

      ! internals
      complex(dp), dimension(size(X)) :: proj_coefficients, wrk

      ! optional input argument
      chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

      info = 0

      proj_coefficients = zero_cdp; wrk = zero_cdp

      ! Orthogonalize vector y w.r.t. to Krylov basis X in two passes of GS.
      ! first pass
      call orthogonalize_against_basis(y, X, info, if_chk_orthonormal=.false., beta=proj_coefficients)
      call check_info(info, 'orthogonalize_against_basis_p1', module=this_module, &
                        & procedure='DGS_vector_against_basis_cdp, pass 1')
      ! second pass
      call orthogonalize_against_basis(y, X, info, if_chk_orthonormal=.false., beta=wrk)
      call check_info(info, 'orthogonalize_against_basis_p2', module=this_module, &
                        & procedure='DGS_vector_against_basis_cdp, pass 2')
      ! combine passes
      proj_coefficients = proj_coefficients + wrk

      if (present(beta)) then
         ! check size
         call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, 'DGS_vector_against_basis_cdp')
         beta = proj_coefficients
      end if

    end subroutine DGS_vector_against_basis_cdp

    subroutine DGS_basis_against_basis_cdp(y, X, info, if_chk_orthonormal, beta)
      !! Computes one step of the double Gram-Schmidt orthogonalization process of the
      !! `abstract_vector` `y` against the `abstract_vector` basis `X`
      class(abstract_vector_cdp), intent(inout) :: Y(:)
      !! Input `abstract_vector` basis to orthogonalize
      class(abstract_vector_cdp), intent(in)    :: X(:)
      !! Input `abstract_vector` basis to orthogonalize against
      integer, intent(out) :: info
      !! Information flag.
      logical,                          optional, intent(in)    :: if_chk_orthonormal
      logical                                                   :: chk_X_orthonormality
      !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
      complex(dp),                         optional, intent(out)   :: beta(:,:)
      !! Projection coefficients if requested

      ! internals
      complex(dp), dimension(size(X),size(Y)) :: proj_coefficients, wrk

      ! optional input argument
      chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

      info = 0

      proj_coefficients = zero_cdp; wrk = zero_cdp

      ! Orthogonalize Krylov basis Y w.r.t. to Krylov basis X in two passes of GS.
      ! first pass
      call orthogonalize_against_basis(Y, X, info, if_chk_orthonormal=.false., beta=proj_coefficients)
      call check_info(info, 'orthogonalize_against_basis_p1', module=this_module, procedure='DGS_basis_against_basis_cdp, first&
          & pass')
      ! second pass
      call orthogonalize_against_basis(Y, X, info, if_chk_orthonormal=.false., beta=wrk)
      call check_info(info, 'orthogonalize_against_basis_p2', module=this_module, procedure='DGS_basis_against_basis_cdp, second&
          & pass')
      ! combine passes
      proj_coefficients = proj_coefficients + wrk

      if (present(beta)) then
         ! check size
         call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, 'DGS_basis_against_basis_cdp')
         beta = proj_coefficients
      end if

    end subroutine DGS_basis_against_basis_cdp  


    !------------------------------------
    !-----     QR FACTORIZATION     -----
    !------------------------------------

    subroutine qr_no_pivoting_rsp(Q, R, info, tol)
        class(abstract_vector_rsp), intent(inout) :: Q(:)
        !! Array of `abstract_vector` to be orthogonalized.
        real(sp), intent(out) :: R(:, :)
        !! Upper triangular matrix \(\mathbf{R}\) resulting from the QR factorization.
        integer, intent(out) :: info
        !! Information flag.
        real(sp), optional, intent(in) :: tol
        !! Tolerance to determine colinearity.
        real(sp)                       :: tolerance

        ! Internal variables.
        real(sp) :: beta
        integer :: j
        logical :: flag
        character(len=128) :: msg

        ! Deals with the optional args.
        tolerance = optval(tol, atol_sp)

        info = 0 ; flag = .false.; R = zero_rsp ; beta = zero_rsp
        do j = 1, size(Q)
            if (j > 1) then
                ! Double Gram-Schmidt orthogonalization
                call double_gram_schmidt_step(Q(j), Q(:j-1), info, if_chk_orthonormal=.false., beta = R(:j-1,j))
                call check_info(info, 'double_gram_schmidt_step', module=this_module, procedure='qr_no_pivoting_rsp')
            end if        

            ! Check for breakdown.
            beta = Q(j)%norm()
            if (abs(beta) < tolerance) then
                if (.not.flag) then
                    flag = .true.
                    info = j
                    write(msg,'(A,I0,A,E15.8)') 'Colinear column detected after ', j, ' steps. beta= ', abs(beta)
                    call logger%log_information(msg, module=this_module, procedure='qr_no_pivoting_rsp')
                end if
                R(j, j) = zero_rsp
                call Q(j)%rand()
                if (j > 1) then
                    call double_gram_schmidt_step(Q(j), Q(:j-1), info, if_chk_orthonormal=.false.)
                    call check_info(info, 'double_gram_schmidt_step', module=this_module, procedure='qr_with_pivoting_rsp')
                end if
                beta = Q(j)%norm()
            else
                R(j, j) = beta
            endif
            ! Normalize column.
            call Q(j)%scal(one_rsp / beta)
        enddo

        return
    end subroutine qr_no_pivoting_rsp

    subroutine qr_with_pivoting_rsp(Q, R, perm, info, tol)
        class(abstract_vector_rsp), intent(inout) :: Q(:)
        !! Array of `abstract_vector` to be orthogonalized.
        real(sp), intent(out) :: R(:, :)
        !! Upper triangular matrix resulting from the QR factorization.
        integer, intent(out) :: perm(size(Q))
        !! Permutation matrix.
        integer, intent(out) :: info
        !! Information flag.
        real(sp), optional, intent(in) :: tol

        !! Tolerance to detect colinearity.
        real(sp) :: tolerance

        ! Internal variables
        real(sp) :: beta
        integer :: idx, i, j, kdim
        integer :: idxv(1)
        real(sp)  :: Rii(size(Q))
        character(len=128) :: msg

        info = 0 ; kdim = size(Q)
        R = zero_rsp ; Rii = zero_rsp
        
        ! Deals with the optional arguments.
        tolerance = optval(tol, atol_sp)

        ! Initialize diagonal entries.
        do i = 1, kdim
            perm(i) = i
            Rii(i) = Q(i)%dot(Q(i))
        enddo

        qr_step: do j = 1, kdim
            idxv = maxloc(abs(Rii)) ; idx = idxv(1)
            if (abs(Rii(idx)) < tolerance) then
                do i = j, kdim
                    call Q(i)%rand()
                    call double_gram_schmidt_step(Q(i), Q(:i-1), info, if_chk_orthonormal=.false.)
                    call check_info(info, 'double_gram_schmidt_step', module=this_module, procedure='qr_with_pivoting_rsp')
                    beta = Q(i)%norm(); call Q(i)%scal(one_rsp / beta)
                enddo
                info = j
                write(msg,'(A,I0,A,E15.8)') 'Breakdown after ', j, ' steps. R_ii= ', abs(Rii(idx))
                call logger%log_information(msg, module=this_module, procedure='qr_with_pivoting_rsp')
                exit qr_step
            endif

            call swap_columns(Q, R, Rii, perm, j, idx)

            ! Check for breakdown.
            beta = Q(j)%norm()
            if (abs(beta) < tolerance) then
                info = j
                R(j, j) = zero_rsp
                call Q(j)%rand()
                call double_gram_schmidt_step(Q(j), Q(:j-1), info, if_chk_orthonormal=.false.)
                call check_info(info, 'double_gram_schmidt_step', module=this_module, procedure='qr_with_pivoting_rsp')
                beta = Q(j)%norm()
            else
                R(j, j) = beta
            endif
            ! Normalize column.
            call Q(j)%scal(one_rsp / beta)

            ! Orthogonalize all columns against new vector.
            do i = j+1, kdim
                beta = Q(j)%dot(Q(i))
                call Q(i)%axpby(one_rsp, Q(j), -beta)
                R(j, i) = beta
            enddo

            ! Update Rii.
            Rii(j) = zero_rsp
            do i = j+1, kdim
                Rii(i) = Rii(i) - R(j, i)**2
            enddo

        enddo qr_step

        return
    end subroutine qr_with_pivoting_rsp

    subroutine swap_columns_rsp(Q, R, Rii, perm, i, j)
        class(abstract_vector_rsp), intent(inout) :: Q(:)
        !! Vector basis whose i-th and j-th columns need swapping.
        real(sp), intent(inout) :: R(:, :)
        !! Upper triangular matrix resulting from QR.
        real(sp), intent(inout) :: Rii(:)
        !! Diagonal entries of R.
        integer, intent(inout) :: perm(:)
        !! Column ordering.
        integer, intent(in) :: i, j
        !! Index of the columns to be swapped.

        ! Internal variables.
        class(abstract_vector_rsp), allocatable :: Qwrk
        real(sp), allocatable :: Rwrk(:)
        integer :: iwrk, m, n

        ! Sanity checks.
        m = size(Q) ; n = min(i, j) - 1

        ! Allocations.
        allocate(Qwrk, source=Q(1)) ; call Qwrk%zero()
        allocate(Rwrk(max(1, n))) ; Rwrk = zero_rsp

        ! Swap columns.
        call copy(Qwrk, Q(j))
        call copy(Q(j), Q(i))
        call copy(Q(i), Qwrk)
        
        Rwrk(1) = Rii(j); Rii(j) = Rii(i); Rii(i) = Rwrk(1)
        iwrk = perm(j); perm(j) = perm(i) ; perm(i) = iwrk

        if (n > 0) then
            Rwrk = R(:n, j) ; R(:n, j) = R(:n, i) ; R(:n, i) = Rwrk
        endif

        return
    end subroutine swap_columns_rsp

    subroutine apply_permutation_matrix_rsp(Q, perm)
        class(abstract_vector_rsp), intent(inout) :: Q(:)
        !! Basis vectors to be permuted.
        integer, intent(in) :: perm(:)
        !! Permutation matrix (vector representation).

        ! Internal variables.
        integer :: i
        class(abstract_vector_rsp), allocatable :: Qwrk(:)

        allocate(Qwrk, source=Q)
        do i = 1, size(perm)
            call copy(Q(i), Qwrk(perm(i)))
        enddo

        return
    end subroutine apply_permutation_matrix_rsp

    subroutine apply_inverse_permutation_matrix_rsp(Q, perm)
        class(abstract_vector_rsp), intent(inout) :: Q(:)
        !! Basis vectors to be (un-) permuted.
        integer, intent(in) :: perm(:)
        !! Permutation matrix (vector representation).

        ! Internal variables.
        integer :: i
        integer :: inv_perm(size(perm))
        class(abstract_vector_rsp), allocatable :: Qwrk(:)

        allocate(Qwrk, source=Q) ; inv_perm = 0

        ! Inverse permutation vector.
        do i = 1, size(perm)
            inv_perm(perm(i)) = i
        enddo

        ! Undo permutation.
        do i = 1, size(perm)
            call copy(Q(i), Qwrk(inv_perm(i)))
        enddo

        return
    end subroutine apply_inverse_permutation_matrix_rsp

    subroutine apply_permutation_matrix_array_rsp(Q, perm)
        real(sp), intent(inout) :: Q(:, :)
        !! Basis vectors to be permuted.
        integer, intent(in) :: perm(:)
        !! Permutation matrix (vector representation).

        ! Internal variables.
        integer :: i
        real(sp), allocatable :: Qwrk(:, :)

        allocate(Qwrk, source=Q)
        do i = 1, size(perm)
            Q(:, i) = Qwrk(:, perm(i))
        enddo

        return
    end subroutine apply_permutation_matrix_array_rsp

    subroutine apply_inverse_permutation_matrix_array_rsp(Q, perm)
        real(sp), intent(inout) :: Q(:, :)
        !! Basis vectors to be (un-) permuted.
        integer, intent(in) :: perm(:)
        !! Permutation matrix (vector representation).

        ! Internal variables.
        integer :: i
        integer :: inv_perm(size(perm))
        real(sp), allocatable :: Qwrk(:, :)

        allocate(Qwrk, source=Q) ; inv_perm = 0

        ! Inverse permutation vector.
        do i = 1, size(perm)
            inv_perm(perm(i)) = i
        enddo

        ! Undo permutation.
        do i = 1, size(perm)
            Q(:, i) = Qwrk(:, inv_perm(i))
        enddo

        return
    end subroutine apply_inverse_permutation_matrix_array_rsp


    subroutine qr_no_pivoting_rdp(Q, R, info, tol)
        class(abstract_vector_rdp), intent(inout) :: Q(:)
        !! Array of `abstract_vector` to be orthogonalized.
        real(dp), intent(out) :: R(:, :)
        !! Upper triangular matrix \(\mathbf{R}\) resulting from the QR factorization.
        integer, intent(out) :: info
        !! Information flag.
        real(dp), optional, intent(in) :: tol
        !! Tolerance to determine colinearity.
        real(dp)                       :: tolerance

        ! Internal variables.
        real(dp) :: beta
        integer :: j
        logical :: flag
        character(len=128) :: msg

        ! Deals with the optional args.
        tolerance = optval(tol, atol_dp)

        info = 0 ; flag = .false.; R = zero_rdp ; beta = zero_rdp
        do j = 1, size(Q)
            if (j > 1) then
                ! Double Gram-Schmidt orthogonalization
                call double_gram_schmidt_step(Q(j), Q(:j-1), info, if_chk_orthonormal=.false., beta = R(:j-1,j))
                call check_info(info, 'double_gram_schmidt_step', module=this_module, procedure='qr_no_pivoting_rdp')
            end if        

            ! Check for breakdown.
            beta = Q(j)%norm()
            if (abs(beta) < tolerance) then
                if (.not.flag) then
                    flag = .true.
                    info = j
                    write(msg,'(A,I0,A,E15.8)') 'Colinear column detected after ', j, ' steps. beta= ', abs(beta)
                    call logger%log_information(msg, module=this_module, procedure='qr_no_pivoting_rdp')
                end if
                R(j, j) = zero_rdp
                call Q(j)%rand()
                if (j > 1) then
                    call double_gram_schmidt_step(Q(j), Q(:j-1), info, if_chk_orthonormal=.false.)
                    call check_info(info, 'double_gram_schmidt_step', module=this_module, procedure='qr_with_pivoting_rdp')
                end if
                beta = Q(j)%norm()
            else
                R(j, j) = beta
            endif
            ! Normalize column.
            call Q(j)%scal(one_rdp / beta)
        enddo

        return
    end subroutine qr_no_pivoting_rdp

    subroutine qr_with_pivoting_rdp(Q, R, perm, info, tol)
        class(abstract_vector_rdp), intent(inout) :: Q(:)
        !! Array of `abstract_vector` to be orthogonalized.
        real(dp), intent(out) :: R(:, :)
        !! Upper triangular matrix resulting from the QR factorization.
        integer, intent(out) :: perm(size(Q))
        !! Permutation matrix.
        integer, intent(out) :: info
        !! Information flag.
        real(dp), optional, intent(in) :: tol

        !! Tolerance to detect colinearity.
        real(dp) :: tolerance

        ! Internal variables
        real(dp) :: beta
        integer :: idx, i, j, kdim
        integer :: idxv(1)
        real(dp)  :: Rii(size(Q))
        character(len=128) :: msg

        info = 0 ; kdim = size(Q)
        R = zero_rdp ; Rii = zero_rdp
        
        ! Deals with the optional arguments.
        tolerance = optval(tol, atol_dp)

        ! Initialize diagonal entries.
        do i = 1, kdim
            perm(i) = i
            Rii(i) = Q(i)%dot(Q(i))
        enddo

        qr_step: do j = 1, kdim
            idxv = maxloc(abs(Rii)) ; idx = idxv(1)
            if (abs(Rii(idx)) < tolerance) then
                do i = j, kdim
                    call Q(i)%rand()
                    call double_gram_schmidt_step(Q(i), Q(:i-1), info, if_chk_orthonormal=.false.)
                    call check_info(info, 'double_gram_schmidt_step', module=this_module, procedure='qr_with_pivoting_rdp')
                    beta = Q(i)%norm(); call Q(i)%scal(one_rdp / beta)
                enddo
                info = j
                write(msg,'(A,I0,A,E15.8)') 'Breakdown after ', j, ' steps. R_ii= ', abs(Rii(idx))
                call logger%log_information(msg, module=this_module, procedure='qr_with_pivoting_rdp')
                exit qr_step
            endif

            call swap_columns(Q, R, Rii, perm, j, idx)

            ! Check for breakdown.
            beta = Q(j)%norm()
            if (abs(beta) < tolerance) then
                info = j
                R(j, j) = zero_rdp
                call Q(j)%rand()
                call double_gram_schmidt_step(Q(j), Q(:j-1), info, if_chk_orthonormal=.false.)
                call check_info(info, 'double_gram_schmidt_step', module=this_module, procedure='qr_with_pivoting_rdp')
                beta = Q(j)%norm()
            else
                R(j, j) = beta
            endif
            ! Normalize column.
            call Q(j)%scal(one_rdp / beta)

            ! Orthogonalize all columns against new vector.
            do i = j+1, kdim
                beta = Q(j)%dot(Q(i))
                call Q(i)%axpby(one_rdp, Q(j), -beta)
                R(j, i) = beta
            enddo

            ! Update Rii.
            Rii(j) = zero_rdp
            do i = j+1, kdim
                Rii(i) = Rii(i) - R(j, i)**2
            enddo

        enddo qr_step

        return
    end subroutine qr_with_pivoting_rdp

    subroutine swap_columns_rdp(Q, R, Rii, perm, i, j)
        class(abstract_vector_rdp), intent(inout) :: Q(:)
        !! Vector basis whose i-th and j-th columns need swapping.
        real(dp), intent(inout) :: R(:, :)
        !! Upper triangular matrix resulting from QR.
        real(dp), intent(inout) :: Rii(:)
        !! Diagonal entries of R.
        integer, intent(inout) :: perm(:)
        !! Column ordering.
        integer, intent(in) :: i, j
        !! Index of the columns to be swapped.

        ! Internal variables.
        class(abstract_vector_rdp), allocatable :: Qwrk
        real(dp), allocatable :: Rwrk(:)
        integer :: iwrk, m, n

        ! Sanity checks.
        m = size(Q) ; n = min(i, j) - 1

        ! Allocations.
        allocate(Qwrk, source=Q(1)) ; call Qwrk%zero()
        allocate(Rwrk(max(1, n))) ; Rwrk = zero_rdp

        ! Swap columns.
        call copy(Qwrk, Q(j))
        call copy(Q(j), Q(i))
        call copy(Q(i), Qwrk)
        
        Rwrk(1) = Rii(j); Rii(j) = Rii(i); Rii(i) = Rwrk(1)
        iwrk = perm(j); perm(j) = perm(i) ; perm(i) = iwrk

        if (n > 0) then
            Rwrk = R(:n, j) ; R(:n, j) = R(:n, i) ; R(:n, i) = Rwrk
        endif

        return
    end subroutine swap_columns_rdp

    subroutine apply_permutation_matrix_rdp(Q, perm)
        class(abstract_vector_rdp), intent(inout) :: Q(:)
        !! Basis vectors to be permuted.
        integer, intent(in) :: perm(:)
        !! Permutation matrix (vector representation).

        ! Internal variables.
        integer :: i
        class(abstract_vector_rdp), allocatable :: Qwrk(:)

        allocate(Qwrk, source=Q)
        do i = 1, size(perm)
            call copy(Q(i), Qwrk(perm(i)))
        enddo

        return
    end subroutine apply_permutation_matrix_rdp

    subroutine apply_inverse_permutation_matrix_rdp(Q, perm)
        class(abstract_vector_rdp), intent(inout) :: Q(:)
        !! Basis vectors to be (un-) permuted.
        integer, intent(in) :: perm(:)
        !! Permutation matrix (vector representation).

        ! Internal variables.
        integer :: i
        integer :: inv_perm(size(perm))
        class(abstract_vector_rdp), allocatable :: Qwrk(:)

        allocate(Qwrk, source=Q) ; inv_perm = 0

        ! Inverse permutation vector.
        do i = 1, size(perm)
            inv_perm(perm(i)) = i
        enddo

        ! Undo permutation.
        do i = 1, size(perm)
            call copy(Q(i), Qwrk(inv_perm(i)))
        enddo

        return
    end subroutine apply_inverse_permutation_matrix_rdp

    subroutine apply_permutation_matrix_array_rdp(Q, perm)
        real(dp), intent(inout) :: Q(:, :)
        !! Basis vectors to be permuted.
        integer, intent(in) :: perm(:)
        !! Permutation matrix (vector representation).

        ! Internal variables.
        integer :: i
        real(dp), allocatable :: Qwrk(:, :)

        allocate(Qwrk, source=Q)
        do i = 1, size(perm)
            Q(:, i) = Qwrk(:, perm(i))
        enddo

        return
    end subroutine apply_permutation_matrix_array_rdp

    subroutine apply_inverse_permutation_matrix_array_rdp(Q, perm)
        real(dp), intent(inout) :: Q(:, :)
        !! Basis vectors to be (un-) permuted.
        integer, intent(in) :: perm(:)
        !! Permutation matrix (vector representation).

        ! Internal variables.
        integer :: i
        integer :: inv_perm(size(perm))
        real(dp), allocatable :: Qwrk(:, :)

        allocate(Qwrk, source=Q) ; inv_perm = 0

        ! Inverse permutation vector.
        do i = 1, size(perm)
            inv_perm(perm(i)) = i
        enddo

        ! Undo permutation.
        do i = 1, size(perm)
            Q(:, i) = Qwrk(:, inv_perm(i))
        enddo

        return
    end subroutine apply_inverse_permutation_matrix_array_rdp


    subroutine qr_no_pivoting_csp(Q, R, info, tol)
        class(abstract_vector_csp), intent(inout) :: Q(:)
        !! Array of `abstract_vector` to be orthogonalized.
        complex(sp), intent(out) :: R(:, :)
        !! Upper triangular matrix \(\mathbf{R}\) resulting from the QR factorization.
        integer, intent(out) :: info
        !! Information flag.
        real(sp), optional, intent(in) :: tol
        !! Tolerance to determine colinearity.
        real(sp)                       :: tolerance

        ! Internal variables.
        complex(sp) :: beta
        integer :: j
        logical :: flag
        character(len=128) :: msg

        ! Deals with the optional args.
        tolerance = optval(tol, atol_sp)

        info = 0 ; flag = .false.; R = zero_rsp ; beta = zero_rsp
        do j = 1, size(Q)
            if (j > 1) then
                ! Double Gram-Schmidt orthogonalization
                call double_gram_schmidt_step(Q(j), Q(:j-1), info, if_chk_orthonormal=.false., beta = R(:j-1,j))
                call check_info(info, 'double_gram_schmidt_step', module=this_module, procedure='qr_no_pivoting_csp')
            end if        

            ! Check for breakdown.
            beta = Q(j)%norm()
            if (abs(beta) < tolerance) then
                if (.not.flag) then
                    flag = .true.
                    info = j
                    write(msg,'(A,I0,A,E15.8)') 'Colinear column detected after ', j, ' steps. beta= ', abs(beta)
                    call logger%log_information(msg, module=this_module, procedure='qr_no_pivoting_csp')
                end if
                R(j, j) = zero_rsp
                call Q(j)%rand()
                if (j > 1) then
                    call double_gram_schmidt_step(Q(j), Q(:j-1), info, if_chk_orthonormal=.false.)
                    call check_info(info, 'double_gram_schmidt_step', module=this_module, procedure='qr_with_pivoting_csp')
                end if
                beta = Q(j)%norm()
            else
                R(j, j) = beta
            endif
            ! Normalize column.
            call Q(j)%scal(one_rsp / beta)
        enddo

        return
    end subroutine qr_no_pivoting_csp

    subroutine qr_with_pivoting_csp(Q, R, perm, info, tol)
        class(abstract_vector_csp), intent(inout) :: Q(:)
        !! Array of `abstract_vector` to be orthogonalized.
        complex(sp), intent(out) :: R(:, :)
        !! Upper triangular matrix resulting from the QR factorization.
        integer, intent(out) :: perm(size(Q))
        !! Permutation matrix.
        integer, intent(out) :: info
        !! Information flag.
        real(sp), optional, intent(in) :: tol

        !! Tolerance to detect colinearity.
        real(sp) :: tolerance

        ! Internal variables
        complex(sp) :: beta
        integer :: idx, i, j, kdim
        integer :: idxv(1)
        complex(sp)  :: Rii(size(Q))
        character(len=128) :: msg

        info = 0 ; kdim = size(Q)
        R = zero_rsp ; Rii = zero_rsp
        
        ! Deals with the optional arguments.
        tolerance = optval(tol, atol_sp)

        ! Initialize diagonal entries.
        do i = 1, kdim
            perm(i) = i
            Rii(i) = Q(i)%dot(Q(i))
        enddo

        qr_step: do j = 1, kdim
            idxv = maxloc(abs(Rii)) ; idx = idxv(1)
            if (abs(Rii(idx)) < tolerance) then
                do i = j, kdim
                    call Q(i)%rand()
                    call double_gram_schmidt_step(Q(i), Q(:i-1), info, if_chk_orthonormal=.false.)
                    call check_info(info, 'double_gram_schmidt_step', module=this_module, procedure='qr_with_pivoting_csp')
                    beta = Q(i)%norm(); call Q(i)%scal(one_csp / beta)
                enddo
                info = j
                write(msg,'(A,I0,A,E15.8)') 'Breakdown after ', j, ' steps. R_ii= ', abs(Rii(idx))
                call logger%log_information(msg, module=this_module, procedure='qr_with_pivoting_csp')
                exit qr_step
            endif

            call swap_columns(Q, R, Rii, perm, j, idx)

            ! Check for breakdown.
            beta = Q(j)%norm()
            if (abs(beta) < tolerance) then
                info = j
                R(j, j) = zero_rsp
                call Q(j)%rand()
                call double_gram_schmidt_step(Q(j), Q(:j-1), info, if_chk_orthonormal=.false.)
                call check_info(info, 'double_gram_schmidt_step', module=this_module, procedure='qr_with_pivoting_csp')
                beta = Q(j)%norm()
            else
                R(j, j) = beta
            endif
            ! Normalize column.
            call Q(j)%scal(one_rsp / beta)

            ! Orthogonalize all columns against new vector.
            do i = j+1, kdim
                beta = Q(j)%dot(Q(i))
                call Q(i)%axpby(one_csp, Q(j), -beta)
                R(j, i) = beta
            enddo

            ! Update Rii.
            Rii(j) = zero_rsp
            do i = j+1, kdim
                Rii(i) = Rii(i) - R(j, i)**2
            enddo

        enddo qr_step

        return
    end subroutine qr_with_pivoting_csp

    subroutine swap_columns_csp(Q, R, Rii, perm, i, j)
        class(abstract_vector_csp), intent(inout) :: Q(:)
        !! Vector basis whose i-th and j-th columns need swapping.
        complex(sp), intent(inout) :: R(:, :)
        !! Upper triangular matrix resulting from QR.
        complex(sp), intent(inout) :: Rii(:)
        !! Diagonal entries of R.
        integer, intent(inout) :: perm(:)
        !! Column ordering.
        integer, intent(in) :: i, j
        !! Index of the columns to be swapped.

        ! Internal variables.
        class(abstract_vector_csp), allocatable :: Qwrk
        complex(sp), allocatable :: Rwrk(:)
        integer :: iwrk, m, n

        ! Sanity checks.
        m = size(Q) ; n = min(i, j) - 1

        ! Allocations.
        allocate(Qwrk, source=Q(1)) ; call Qwrk%zero()
        allocate(Rwrk(max(1, n))) ; Rwrk = zero_rsp

        ! Swap columns.
        call copy(Qwrk, Q(j))
        call copy(Q(j), Q(i))
        call copy(Q(i), Qwrk)
        
        Rwrk(1) = Rii(j); Rii(j) = Rii(i); Rii(i) = Rwrk(1)
        iwrk = perm(j); perm(j) = perm(i) ; perm(i) = iwrk

        if (n > 0) then
            Rwrk = R(:n, j) ; R(:n, j) = R(:n, i) ; R(:n, i) = Rwrk
        endif

        return
    end subroutine swap_columns_csp

    subroutine apply_permutation_matrix_csp(Q, perm)
        class(abstract_vector_csp), intent(inout) :: Q(:)
        !! Basis vectors to be permuted.
        integer, intent(in) :: perm(:)
        !! Permutation matrix (vector representation).

        ! Internal variables.
        integer :: i
        class(abstract_vector_csp), allocatable :: Qwrk(:)

        allocate(Qwrk, source=Q)
        do i = 1, size(perm)
            call copy(Q(i), Qwrk(perm(i)))
        enddo

        return
    end subroutine apply_permutation_matrix_csp

    subroutine apply_inverse_permutation_matrix_csp(Q, perm)
        class(abstract_vector_csp), intent(inout) :: Q(:)
        !! Basis vectors to be (un-) permuted.
        integer, intent(in) :: perm(:)
        !! Permutation matrix (vector representation).

        ! Internal variables.
        integer :: i
        integer :: inv_perm(size(perm))
        class(abstract_vector_csp), allocatable :: Qwrk(:)

        allocate(Qwrk, source=Q) ; inv_perm = 0

        ! Inverse permutation vector.
        do i = 1, size(perm)
            inv_perm(perm(i)) = i
        enddo

        ! Undo permutation.
        do i = 1, size(perm)
            call copy(Q(i), Qwrk(inv_perm(i)))
        enddo

        return
    end subroutine apply_inverse_permutation_matrix_csp

    subroutine apply_permutation_matrix_array_csp(Q, perm)
        complex(sp), intent(inout) :: Q(:, :)
        !! Basis vectors to be permuted.
        integer, intent(in) :: perm(:)
        !! Permutation matrix (vector representation).

        ! Internal variables.
        integer :: i
        complex(sp), allocatable :: Qwrk(:, :)

        allocate(Qwrk, source=Q)
        do i = 1, size(perm)
            Q(:, i) = Qwrk(:, perm(i))
        enddo

        return
    end subroutine apply_permutation_matrix_array_csp

    subroutine apply_inverse_permutation_matrix_array_csp(Q, perm)
        complex(sp), intent(inout) :: Q(:, :)
        !! Basis vectors to be (un-) permuted.
        integer, intent(in) :: perm(:)
        !! Permutation matrix (vector representation).

        ! Internal variables.
        integer :: i
        integer :: inv_perm(size(perm))
        complex(sp), allocatable :: Qwrk(:, :)

        allocate(Qwrk, source=Q) ; inv_perm = 0

        ! Inverse permutation vector.
        do i = 1, size(perm)
            inv_perm(perm(i)) = i
        enddo

        ! Undo permutation.
        do i = 1, size(perm)
            Q(:, i) = Qwrk(:, inv_perm(i))
        enddo

        return
    end subroutine apply_inverse_permutation_matrix_array_csp


    subroutine qr_no_pivoting_cdp(Q, R, info, tol)
        class(abstract_vector_cdp), intent(inout) :: Q(:)
        !! Array of `abstract_vector` to be orthogonalized.
        complex(dp), intent(out) :: R(:, :)
        !! Upper triangular matrix \(\mathbf{R}\) resulting from the QR factorization.
        integer, intent(out) :: info
        !! Information flag.
        real(dp), optional, intent(in) :: tol
        !! Tolerance to determine colinearity.
        real(dp)                       :: tolerance

        ! Internal variables.
        complex(dp) :: beta
        integer :: j
        logical :: flag
        character(len=128) :: msg

        ! Deals with the optional args.
        tolerance = optval(tol, atol_dp)

        info = 0 ; flag = .false.; R = zero_rdp ; beta = zero_rdp
        do j = 1, size(Q)
            if (j > 1) then
                ! Double Gram-Schmidt orthogonalization
                call double_gram_schmidt_step(Q(j), Q(:j-1), info, if_chk_orthonormal=.false., beta = R(:j-1,j))
                call check_info(info, 'double_gram_schmidt_step', module=this_module, procedure='qr_no_pivoting_cdp')
            end if        

            ! Check for breakdown.
            beta = Q(j)%norm()
            if (abs(beta) < tolerance) then
                if (.not.flag) then
                    flag = .true.
                    info = j
                    write(msg,'(A,I0,A,E15.8)') 'Colinear column detected after ', j, ' steps. beta= ', abs(beta)
                    call logger%log_information(msg, module=this_module, procedure='qr_no_pivoting_cdp')
                end if
                R(j, j) = zero_rdp
                call Q(j)%rand()
                if (j > 1) then
                    call double_gram_schmidt_step(Q(j), Q(:j-1), info, if_chk_orthonormal=.false.)
                    call check_info(info, 'double_gram_schmidt_step', module=this_module, procedure='qr_with_pivoting_cdp')
                end if
                beta = Q(j)%norm()
            else
                R(j, j) = beta
            endif
            ! Normalize column.
            call Q(j)%scal(one_rdp / beta)
        enddo

        return
    end subroutine qr_no_pivoting_cdp

    subroutine qr_with_pivoting_cdp(Q, R, perm, info, tol)
        class(abstract_vector_cdp), intent(inout) :: Q(:)
        !! Array of `abstract_vector` to be orthogonalized.
        complex(dp), intent(out) :: R(:, :)
        !! Upper triangular matrix resulting from the QR factorization.
        integer, intent(out) :: perm(size(Q))
        !! Permutation matrix.
        integer, intent(out) :: info
        !! Information flag.
        real(dp), optional, intent(in) :: tol

        !! Tolerance to detect colinearity.
        real(dp) :: tolerance

        ! Internal variables
        complex(dp) :: beta
        integer :: idx, i, j, kdim
        integer :: idxv(1)
        complex(dp)  :: Rii(size(Q))
        character(len=128) :: msg

        info = 0 ; kdim = size(Q)
        R = zero_rdp ; Rii = zero_rdp
        
        ! Deals with the optional arguments.
        tolerance = optval(tol, atol_dp)

        ! Initialize diagonal entries.
        do i = 1, kdim
            perm(i) = i
            Rii(i) = Q(i)%dot(Q(i))
        enddo

        qr_step: do j = 1, kdim
            idxv = maxloc(abs(Rii)) ; idx = idxv(1)
            if (abs(Rii(idx)) < tolerance) then
                do i = j, kdim
                    call Q(i)%rand()
                    call double_gram_schmidt_step(Q(i), Q(:i-1), info, if_chk_orthonormal=.false.)
                    call check_info(info, 'double_gram_schmidt_step', module=this_module, procedure='qr_with_pivoting_cdp')
                    beta = Q(i)%norm(); call Q(i)%scal(one_cdp / beta)
                enddo
                info = j
                write(msg,'(A,I0,A,E15.8)') 'Breakdown after ', j, ' steps. R_ii= ', abs(Rii(idx))
                call logger%log_information(msg, module=this_module, procedure='qr_with_pivoting_cdp')
                exit qr_step
            endif

            call swap_columns(Q, R, Rii, perm, j, idx)

            ! Check for breakdown.
            beta = Q(j)%norm()
            if (abs(beta) < tolerance) then
                info = j
                R(j, j) = zero_rdp
                call Q(j)%rand()
                call double_gram_schmidt_step(Q(j), Q(:j-1), info, if_chk_orthonormal=.false.)
                call check_info(info, 'double_gram_schmidt_step', module=this_module, procedure='qr_with_pivoting_cdp')
                beta = Q(j)%norm()
            else
                R(j, j) = beta
            endif
            ! Normalize column.
            call Q(j)%scal(one_rdp / beta)

            ! Orthogonalize all columns against new vector.
            do i = j+1, kdim
                beta = Q(j)%dot(Q(i))
                call Q(i)%axpby(one_cdp, Q(j), -beta)
                R(j, i) = beta
            enddo

            ! Update Rii.
            Rii(j) = zero_rdp
            do i = j+1, kdim
                Rii(i) = Rii(i) - R(j, i)**2
            enddo

        enddo qr_step

        return
    end subroutine qr_with_pivoting_cdp

    subroutine swap_columns_cdp(Q, R, Rii, perm, i, j)
        class(abstract_vector_cdp), intent(inout) :: Q(:)
        !! Vector basis whose i-th and j-th columns need swapping.
        complex(dp), intent(inout) :: R(:, :)
        !! Upper triangular matrix resulting from QR.
        complex(dp), intent(inout) :: Rii(:)
        !! Diagonal entries of R.
        integer, intent(inout) :: perm(:)
        !! Column ordering.
        integer, intent(in) :: i, j
        !! Index of the columns to be swapped.

        ! Internal variables.
        class(abstract_vector_cdp), allocatable :: Qwrk
        complex(dp), allocatable :: Rwrk(:)
        integer :: iwrk, m, n

        ! Sanity checks.
        m = size(Q) ; n = min(i, j) - 1

        ! Allocations.
        allocate(Qwrk, source=Q(1)) ; call Qwrk%zero()
        allocate(Rwrk(max(1, n))) ; Rwrk = zero_rdp

        ! Swap columns.
        call copy(Qwrk, Q(j))
        call copy(Q(j), Q(i))
        call copy(Q(i), Qwrk)
        
        Rwrk(1) = Rii(j); Rii(j) = Rii(i); Rii(i) = Rwrk(1)
        iwrk = perm(j); perm(j) = perm(i) ; perm(i) = iwrk

        if (n > 0) then
            Rwrk = R(:n, j) ; R(:n, j) = R(:n, i) ; R(:n, i) = Rwrk
        endif

        return
    end subroutine swap_columns_cdp

    subroutine apply_permutation_matrix_cdp(Q, perm)
        class(abstract_vector_cdp), intent(inout) :: Q(:)
        !! Basis vectors to be permuted.
        integer, intent(in) :: perm(:)
        !! Permutation matrix (vector representation).

        ! Internal variables.
        integer :: i
        class(abstract_vector_cdp), allocatable :: Qwrk(:)

        allocate(Qwrk, source=Q)
        do i = 1, size(perm)
            call copy(Q(i), Qwrk(perm(i)))
        enddo

        return
    end subroutine apply_permutation_matrix_cdp

    subroutine apply_inverse_permutation_matrix_cdp(Q, perm)
        class(abstract_vector_cdp), intent(inout) :: Q(:)
        !! Basis vectors to be (un-) permuted.
        integer, intent(in) :: perm(:)
        !! Permutation matrix (vector representation).

        ! Internal variables.
        integer :: i
        integer :: inv_perm(size(perm))
        class(abstract_vector_cdp), allocatable :: Qwrk(:)

        allocate(Qwrk, source=Q) ; inv_perm = 0

        ! Inverse permutation vector.
        do i = 1, size(perm)
            inv_perm(perm(i)) = i
        enddo

        ! Undo permutation.
        do i = 1, size(perm)
            call copy(Q(i), Qwrk(inv_perm(i)))
        enddo

        return
    end subroutine apply_inverse_permutation_matrix_cdp

    subroutine apply_permutation_matrix_array_cdp(Q, perm)
        complex(dp), intent(inout) :: Q(:, :)
        !! Basis vectors to be permuted.
        integer, intent(in) :: perm(:)
        !! Permutation matrix (vector representation).

        ! Internal variables.
        integer :: i
        complex(dp), allocatable :: Qwrk(:, :)

        allocate(Qwrk, source=Q)
        do i = 1, size(perm)
            Q(:, i) = Qwrk(:, perm(i))
        enddo

        return
    end subroutine apply_permutation_matrix_array_cdp

    subroutine apply_inverse_permutation_matrix_array_cdp(Q, perm)
        complex(dp), intent(inout) :: Q(:, :)
        !! Basis vectors to be (un-) permuted.
        integer, intent(in) :: perm(:)
        !! Permutation matrix (vector representation).

        ! Internal variables.
        integer :: i
        integer :: inv_perm(size(perm))
        complex(dp), allocatable :: Qwrk(:, :)

        allocate(Qwrk, source=Q) ; inv_perm = 0

        ! Inverse permutation vector.
        do i = 1, size(perm)
            inv_perm(perm(i)) = i
        enddo

        ! Undo permutation.
        do i = 1, size(perm)
            Q(:, i) = Qwrk(:, inv_perm(i))
        enddo

        return
    end subroutine apply_inverse_permutation_matrix_array_cdp



    !-----------------------------------------
    !-----     ARNOLDI FACTORIZATION     -----
    !-----------------------------------------

    subroutine arnoldi_rsp(A, X, H, info, kstart, kend, tol, transpose, blksize)
        class(abstract_linop_rsp), intent(in) :: A
        !! Linear operator to be factorized.
        class(abstract_vector_rsp), intent(inout) :: X(:)
        !! Orthogonal basis for the generated Krylov subspace.
        real(sp), intent(inout) :: H(:, :)
        !! Upper Hessenberg matrix.
        integer, intent(out) :: info
        !! Information flag.
        integer, optional, intent(in) :: kstart
        !! Starting index for the Arnoldi factorization (default 1).
        integer, optional, intent(in) :: kend
        !! Final index for the Arnoldi factorization (default `size(X)-1`)
        logical, optional, intent(in) :: transpose
        !! Whether \( \mathbf{A} \) is being transposed or not (default `.false.`)
        real(sp), optional, intent(in) :: tol
        !! Tolerance to determine whether an invariant subspace has been computed or not.
        integer, optional, intent(in) :: blksize
        !! Block size for block Arnoldi (default 1).

        ! Internal variables.
        integer :: k_start, k_end, p
        logical :: trans
        real(sp) :: tolerance
        real(sp) :: beta
        real(sp), allocatable :: res(:)
        integer, allocatable :: perm(:)
        integer :: k, i, kdim, kpm, kp, kpp

        ! Deals with optional non-unity blksize and allocations.
        p = optval(blksize, 1) ; allocate(res(p)) ; res = zero_rsp
        allocate(perm(size(H, 2))) ; perm = 0 ; info = 0

        ! Check dimensions.
        kdim = (size(X) - p) / p

        ! Deal with the other optional args.
        k_start = optval(kstart, 1) ; k_end = optval(kend, kdim)
        tolerance = optval (tol, atol_sp)
        trans     = optval(transpose, .false.)

        ! Arnoldi factorization.
        blk_arnoldi: do k = k_start, k_end
            ! Counters
            kpm = (k - 1) * p ; kp = kpm + p ; kpp = kp + p

            ! Matrix-vector product.
            if (trans) then
                do i = 1, p
                    call A%rmatvec(X(kpm+i), X(kp+i))
                enddo
            else
                do i = 1, p
                    call A%matvec(X(kpm+i), X(kp+i))
                enddo
            endif

            ! Update Hessenberg matrix via batch double Gram-Schmidt step.
            call double_gram_schmidt_step(X(kp+1:kpp), X(:kp), info, if_chk_orthonormal=.false., beta=H(:kp, kpm+1:kp))
            call check_info(info, 'double_gram_schmidt_step', module=this_module, procedure='arnoldi_rsp')

            ! Orthogonalize current blk vectors.
            call qr(X(kp+1:kpp), H(kp+1:kpp, kpm+1:kp), info)
            call check_info(info, 'qr', module=this_module, procedure='arnoldi_rsp')

            ! Extract residual norm (smallest diagonal element of H matrix).
            res = zero_rsp
            do i = 1, p
                res(i) = H(kp+i, kpm+i)
            enddo
            beta = minval(abs(res))

            ! Exit Arnoldi loop if needed.
            if (beta < tolerance) then
                ! Dimension of the computed invariant subspace.
                info = kp
                ! Exit the Arnoldi iteration.
                exit blk_arnoldi
            endif

        enddo blk_arnoldi

        return
    end subroutine arnoldi_rsp

    subroutine arnoldi_rdp(A, X, H, info, kstart, kend, tol, transpose, blksize)
        class(abstract_linop_rdp), intent(in) :: A
        !! Linear operator to be factorized.
        class(abstract_vector_rdp), intent(inout) :: X(:)
        !! Orthogonal basis for the generated Krylov subspace.
        real(dp), intent(inout) :: H(:, :)
        !! Upper Hessenberg matrix.
        integer, intent(out) :: info
        !! Information flag.
        integer, optional, intent(in) :: kstart
        !! Starting index for the Arnoldi factorization (default 1).
        integer, optional, intent(in) :: kend
        !! Final index for the Arnoldi factorization (default `size(X)-1`)
        logical, optional, intent(in) :: transpose
        !! Whether \( \mathbf{A} \) is being transposed or not (default `.false.`)
        real(dp), optional, intent(in) :: tol
        !! Tolerance to determine whether an invariant subspace has been computed or not.
        integer, optional, intent(in) :: blksize
        !! Block size for block Arnoldi (default 1).

        ! Internal variables.
        integer :: k_start, k_end, p
        logical :: trans
        real(dp) :: tolerance
        real(dp) :: beta
        real(dp), allocatable :: res(:)
        integer, allocatable :: perm(:)
        integer :: k, i, kdim, kpm, kp, kpp

        ! Deals with optional non-unity blksize and allocations.
        p = optval(blksize, 1) ; allocate(res(p)) ; res = zero_rdp
        allocate(perm(size(H, 2))) ; perm = 0 ; info = 0

        ! Check dimensions.
        kdim = (size(X) - p) / p

        ! Deal with the other optional args.
        k_start = optval(kstart, 1) ; k_end = optval(kend, kdim)
        tolerance = optval (tol, atol_dp)
        trans     = optval(transpose, .false.)

        ! Arnoldi factorization.
        blk_arnoldi: do k = k_start, k_end
            ! Counters
            kpm = (k - 1) * p ; kp = kpm + p ; kpp = kp + p

            ! Matrix-vector product.
            if (trans) then
                do i = 1, p
                    call A%rmatvec(X(kpm+i), X(kp+i))
                enddo
            else
                do i = 1, p
                    call A%matvec(X(kpm+i), X(kp+i))
                enddo
            endif

            ! Update Hessenberg matrix via batch double Gram-Schmidt step.
            call double_gram_schmidt_step(X(kp+1:kpp), X(:kp), info, if_chk_orthonormal=.false., beta=H(:kp, kpm+1:kp))
            call check_info(info, 'double_gram_schmidt_step', module=this_module, procedure='arnoldi_rdp')

            ! Orthogonalize current blk vectors.
            call qr(X(kp+1:kpp), H(kp+1:kpp, kpm+1:kp), info)
            call check_info(info, 'qr', module=this_module, procedure='arnoldi_rdp')

            ! Extract residual norm (smallest diagonal element of H matrix).
            res = zero_rdp
            do i = 1, p
                res(i) = H(kp+i, kpm+i)
            enddo
            beta = minval(abs(res))

            ! Exit Arnoldi loop if needed.
            if (beta < tolerance) then
                ! Dimension of the computed invariant subspace.
                info = kp
                ! Exit the Arnoldi iteration.
                exit blk_arnoldi
            endif

        enddo blk_arnoldi

        return
    end subroutine arnoldi_rdp

    subroutine arnoldi_csp(A, X, H, info, kstart, kend, tol, transpose, blksize)
        class(abstract_linop_csp), intent(in) :: A
        !! Linear operator to be factorized.
        class(abstract_vector_csp), intent(inout) :: X(:)
        !! Orthogonal basis for the generated Krylov subspace.
        complex(sp), intent(inout) :: H(:, :)
        !! Upper Hessenberg matrix.
        integer, intent(out) :: info
        !! Information flag.
        integer, optional, intent(in) :: kstart
        !! Starting index for the Arnoldi factorization (default 1).
        integer, optional, intent(in) :: kend
        !! Final index for the Arnoldi factorization (default `size(X)-1`)
        logical, optional, intent(in) :: transpose
        !! Whether \( \mathbf{A} \) is being transposed or not (default `.false.`)
        real(sp), optional, intent(in) :: tol
        !! Tolerance to determine whether an invariant subspace has been computed or not.
        integer, optional, intent(in) :: blksize
        !! Block size for block Arnoldi (default 1).

        ! Internal variables.
        integer :: k_start, k_end, p
        logical :: trans
        real(sp) :: tolerance
        real(sp) :: beta
        complex(sp), allocatable :: res(:)
        integer, allocatable :: perm(:)
        integer :: k, i, kdim, kpm, kp, kpp

        ! Deals with optional non-unity blksize and allocations.
        p = optval(blksize, 1) ; allocate(res(p)) ; res = zero_rsp
        allocate(perm(size(H, 2))) ; perm = 0 ; info = 0

        ! Check dimensions.
        kdim = (size(X) - p) / p

        ! Deal with the other optional args.
        k_start = optval(kstart, 1) ; k_end = optval(kend, kdim)
        tolerance = optval (tol, atol_sp)
        trans     = optval(transpose, .false.)

        ! Arnoldi factorization.
        blk_arnoldi: do k = k_start, k_end
            ! Counters
            kpm = (k - 1) * p ; kp = kpm + p ; kpp = kp + p

            ! Matrix-vector product.
            if (trans) then
                do i = 1, p
                    call A%rmatvec(X(kpm+i), X(kp+i))
                enddo
            else
                do i = 1, p
                    call A%matvec(X(kpm+i), X(kp+i))
                enddo
            endif

            ! Update Hessenberg matrix via batch double Gram-Schmidt step.
            call double_gram_schmidt_step(X(kp+1:kpp), X(:kp), info, if_chk_orthonormal=.false., beta=H(:kp, kpm+1:kp))
            call check_info(info, 'double_gram_schmidt_step', module=this_module, procedure='arnoldi_csp')

            ! Orthogonalize current blk vectors.
            call qr(X(kp+1:kpp), H(kp+1:kpp, kpm+1:kp), info)
            call check_info(info, 'qr', module=this_module, procedure='arnoldi_csp')

            ! Extract residual norm (smallest diagonal element of H matrix).
            res = zero_rsp
            do i = 1, p
                res(i) = H(kp+i, kpm+i)
            enddo
            beta = minval(abs(res))

            ! Exit Arnoldi loop if needed.
            if (beta < tolerance) then
                ! Dimension of the computed invariant subspace.
                info = kp
                ! Exit the Arnoldi iteration.
                exit blk_arnoldi
            endif

        enddo blk_arnoldi

        return
    end subroutine arnoldi_csp

    subroutine arnoldi_cdp(A, X, H, info, kstart, kend, tol, transpose, blksize)
        class(abstract_linop_cdp), intent(in) :: A
        !! Linear operator to be factorized.
        class(abstract_vector_cdp), intent(inout) :: X(:)
        !! Orthogonal basis for the generated Krylov subspace.
        complex(dp), intent(inout) :: H(:, :)
        !! Upper Hessenberg matrix.
        integer, intent(out) :: info
        !! Information flag.
        integer, optional, intent(in) :: kstart
        !! Starting index for the Arnoldi factorization (default 1).
        integer, optional, intent(in) :: kend
        !! Final index for the Arnoldi factorization (default `size(X)-1`)
        logical, optional, intent(in) :: transpose
        !! Whether \( \mathbf{A} \) is being transposed or not (default `.false.`)
        real(dp), optional, intent(in) :: tol
        !! Tolerance to determine whether an invariant subspace has been computed or not.
        integer, optional, intent(in) :: blksize
        !! Block size for block Arnoldi (default 1).

        ! Internal variables.
        integer :: k_start, k_end, p
        logical :: trans
        real(dp) :: tolerance
        real(dp) :: beta
        complex(dp), allocatable :: res(:)
        integer, allocatable :: perm(:)
        integer :: k, i, kdim, kpm, kp, kpp

        ! Deals with optional non-unity blksize and allocations.
        p = optval(blksize, 1) ; allocate(res(p)) ; res = zero_rdp
        allocate(perm(size(H, 2))) ; perm = 0 ; info = 0

        ! Check dimensions.
        kdim = (size(X) - p) / p

        ! Deal with the other optional args.
        k_start = optval(kstart, 1) ; k_end = optval(kend, kdim)
        tolerance = optval (tol, atol_dp)
        trans     = optval(transpose, .false.)

        ! Arnoldi factorization.
        blk_arnoldi: do k = k_start, k_end
            ! Counters
            kpm = (k - 1) * p ; kp = kpm + p ; kpp = kp + p

            ! Matrix-vector product.
            if (trans) then
                do i = 1, p
                    call A%rmatvec(X(kpm+i), X(kp+i))
                enddo
            else
                do i = 1, p
                    call A%matvec(X(kpm+i), X(kp+i))
                enddo
            endif

            ! Update Hessenberg matrix via batch double Gram-Schmidt step.
            call double_gram_schmidt_step(X(kp+1:kpp), X(:kp), info, if_chk_orthonormal=.false., beta=H(:kp, kpm+1:kp))
            call check_info(info, 'double_gram_schmidt_step', module=this_module, procedure='arnoldi_cdp')

            ! Orthogonalize current blk vectors.
            call qr(X(kp+1:kpp), H(kp+1:kpp, kpm+1:kp), info)
            call check_info(info, 'qr', module=this_module, procedure='arnoldi_cdp')

            ! Extract residual norm (smallest diagonal element of H matrix).
            res = zero_rdp
            do i = 1, p
                res(i) = H(kp+i, kpm+i)
            enddo
            beta = minval(abs(res))

            ! Exit Arnoldi loop if needed.
            if (beta < tolerance) then
                ! Dimension of the computed invariant subspace.
                info = kp
                ! Exit the Arnoldi iteration.
                exit blk_arnoldi
            endif

        enddo blk_arnoldi

        return
    end subroutine arnoldi_cdp



    !---------------------------------------------
    !-----     LANCZOS BIDIAGONALIZATION     -----
    !---------------------------------------------

    subroutine lanczos_bidiagonalization_rsp(A, U, V, B, info, kstart, kend, tol)
        class(abstract_linop_rsp), intent(in) :: A
        !! Linear operator to be factorized.
        class(abstract_vector_rsp), intent(inout) :: U(:)
        !! Orthonormal basis for the column span of \(\mathbf{A}\). On entry, `U(1)` needs to be set to
        !! the starting Krylov vector.
        class(abstract_vector_rsp), intent(inout) :: V(:)
        !! Orthonormal basis for the row span of \(\mathbf{A}\).
        real(sp), intent(inout) :: B(:, :)
        !! Bidiagonal matrix.
        integer, intent(out) :: info
        !! Information flag.
        integer, optional, intent(in) :: kstart
        !! Starting index for the Lanczos factorization (default 1).
        integer, optional, intent(in) :: kend
        !! Final index for the Lanczos factorization (default 1).
        real(sp), optional, intent(in) :: tol
        !! Tolerance to determine whether invariant subspaces have been computed or not.

        ! Internal variables.
        integer :: k_start, k_end
        real(sp) :: tolerance
        real(sp) :: alpha, beta
        integer :: k, kdim

        info = 0

        ! Krylov subspace dimension.
        kdim = size(U) - 1

        ! Deals with the optional args.
        k_start = optval(kstart, 1)
        k_end   = optval(kend, kdim)
        tolerance = optval(tol, atol_sp)

        ! Lanczos bidiagonalization.
        lanczos : do k = k_start, k_end
            ! Transpose matrix-vector product.
            call A%rmatvec(U(k), V(k))

            ! Full re-orthogonalization of the right Krylov subspace.
            if (k > 1) then
                call double_gram_schmidt_step(V(k), V(:k-1), info, if_chk_orthonormal=.false.)
                call check_info(info, 'double_gram_schmidt_step', module=this_module, &
                                    & procedure='lanczos_bidiagonalization_rsp, right basis')
            end if

            ! Normalization step.
            alpha = V(k)%norm() ; B(k, k) = alpha
            if (abs(alpha) > tolerance) then
                call V(k)%scal(one_rsp/alpha)
            else
                info = k
                exit lanczos
            endif

            ! Matrix-vector product.
            call A%matvec(V(k), U(k+1))

            ! Full re-orthogonalization of the left Krylov subspace.
            call double_gram_schmidt_step(U(k+1), U(:k), info, if_chk_orthonormal=.false.)
            call check_info(info, 'double_gram_schmidt_step', module=this_module, &
                                    & procedure='lanczos_bidiagonalization_rsp, left basis')

            ! Normalization step
            beta = U(k+1)%norm() ; B(k+1, k) = beta
            if (abs(beta) > tolerance) then
                call U(k+1)%scal(one_rsp / beta)
            else
                info = k
                exit lanczos
            endif

        enddo lanczos

        return
    end subroutine lanczos_bidiagonalization_rsp

    subroutine lanczos_bidiagonalization_rdp(A, U, V, B, info, kstart, kend, tol)
        class(abstract_linop_rdp), intent(in) :: A
        !! Linear operator to be factorized.
        class(abstract_vector_rdp), intent(inout) :: U(:)
        !! Orthonormal basis for the column span of \(\mathbf{A}\). On entry, `U(1)` needs to be set to
        !! the starting Krylov vector.
        class(abstract_vector_rdp), intent(inout) :: V(:)
        !! Orthonormal basis for the row span of \(\mathbf{A}\).
        real(dp), intent(inout) :: B(:, :)
        !! Bidiagonal matrix.
        integer, intent(out) :: info
        !! Information flag.
        integer, optional, intent(in) :: kstart
        !! Starting index for the Lanczos factorization (default 1).
        integer, optional, intent(in) :: kend
        !! Final index for the Lanczos factorization (default 1).
        real(dp), optional, intent(in) :: tol
        !! Tolerance to determine whether invariant subspaces have been computed or not.

        ! Internal variables.
        integer :: k_start, k_end
        real(dp) :: tolerance
        real(dp) :: alpha, beta
        integer :: k, kdim

        info = 0

        ! Krylov subspace dimension.
        kdim = size(U) - 1

        ! Deals with the optional args.
        k_start = optval(kstart, 1)
        k_end   = optval(kend, kdim)
        tolerance = optval(tol, atol_dp)

        ! Lanczos bidiagonalization.
        lanczos : do k = k_start, k_end
            ! Transpose matrix-vector product.
            call A%rmatvec(U(k), V(k))

            ! Full re-orthogonalization of the right Krylov subspace.
            if (k > 1) then
                call double_gram_schmidt_step(V(k), V(:k-1), info, if_chk_orthonormal=.false.)
                call check_info(info, 'double_gram_schmidt_step', module=this_module, &
                                    & procedure='lanczos_bidiagonalization_rdp, right basis')
            end if

            ! Normalization step.
            alpha = V(k)%norm() ; B(k, k) = alpha
            if (abs(alpha) > tolerance) then
                call V(k)%scal(one_rdp/alpha)
            else
                info = k
                exit lanczos
            endif

            ! Matrix-vector product.
            call A%matvec(V(k), U(k+1))

            ! Full re-orthogonalization of the left Krylov subspace.
            call double_gram_schmidt_step(U(k+1), U(:k), info, if_chk_orthonormal=.false.)
            call check_info(info, 'double_gram_schmidt_step', module=this_module, &
                                    & procedure='lanczos_bidiagonalization_rdp, left basis')

            ! Normalization step
            beta = U(k+1)%norm() ; B(k+1, k) = beta
            if (abs(beta) > tolerance) then
                call U(k+1)%scal(one_rdp / beta)
            else
                info = k
                exit lanczos
            endif

        enddo lanczos

        return
    end subroutine lanczos_bidiagonalization_rdp

    subroutine lanczos_bidiagonalization_csp(A, U, V, B, info, kstart, kend, tol)
        class(abstract_linop_csp), intent(in) :: A
        !! Linear operator to be factorized.
        class(abstract_vector_csp), intent(inout) :: U(:)
        !! Orthonormal basis for the column span of \(\mathbf{A}\). On entry, `U(1)` needs to be set to
        !! the starting Krylov vector.
        class(abstract_vector_csp), intent(inout) :: V(:)
        !! Orthonormal basis for the row span of \(\mathbf{A}\).
        complex(sp), intent(inout) :: B(:, :)
        !! Bidiagonal matrix.
        integer, intent(out) :: info
        !! Information flag.
        integer, optional, intent(in) :: kstart
        !! Starting index for the Lanczos factorization (default 1).
        integer, optional, intent(in) :: kend
        !! Final index for the Lanczos factorization (default 1).
        real(sp), optional, intent(in) :: tol
        !! Tolerance to determine whether invariant subspaces have been computed or not.

        ! Internal variables.
        integer :: k_start, k_end
        real(sp) :: tolerance
        complex(sp) :: alpha, beta
        integer :: k, kdim

        info = 0

        ! Krylov subspace dimension.
        kdim = size(U) - 1

        ! Deals with the optional args.
        k_start = optval(kstart, 1)
        k_end   = optval(kend, kdim)
        tolerance = optval(tol, atol_sp)

        ! Lanczos bidiagonalization.
        lanczos : do k = k_start, k_end
            ! Transpose matrix-vector product.
            call A%rmatvec(U(k), V(k))

            ! Full re-orthogonalization of the right Krylov subspace.
            if (k > 1) then
                call double_gram_schmidt_step(V(k), V(:k-1), info, if_chk_orthonormal=.false.)
                call check_info(info, 'double_gram_schmidt_step', module=this_module, &
                                    & procedure='lanczos_bidiagonalization_csp, right basis')
            end if

            ! Normalization step.
            alpha = V(k)%norm() ; B(k, k) = alpha
            if (abs(alpha) > tolerance) then
                call V(k)%scal(one_csp/alpha)
            else
                info = k
                exit lanczos
            endif

            ! Matrix-vector product.
            call A%matvec(V(k), U(k+1))

            ! Full re-orthogonalization of the left Krylov subspace.
            call double_gram_schmidt_step(U(k+1), U(:k), info, if_chk_orthonormal=.false.)
            call check_info(info, 'double_gram_schmidt_step', module=this_module, &
                                    & procedure='lanczos_bidiagonalization_csp, left basis')

            ! Normalization step
            beta = U(k+1)%norm() ; B(k+1, k) = beta
            if (abs(beta) > tolerance) then
                call U(k+1)%scal(one_csp / beta)
            else
                info = k
                exit lanczos
            endif

        enddo lanczos

        return
    end subroutine lanczos_bidiagonalization_csp

    subroutine lanczos_bidiagonalization_cdp(A, U, V, B, info, kstart, kend, tol)
        class(abstract_linop_cdp), intent(in) :: A
        !! Linear operator to be factorized.
        class(abstract_vector_cdp), intent(inout) :: U(:)
        !! Orthonormal basis for the column span of \(\mathbf{A}\). On entry, `U(1)` needs to be set to
        !! the starting Krylov vector.
        class(abstract_vector_cdp), intent(inout) :: V(:)
        !! Orthonormal basis for the row span of \(\mathbf{A}\).
        complex(dp), intent(inout) :: B(:, :)
        !! Bidiagonal matrix.
        integer, intent(out) :: info
        !! Information flag.
        integer, optional, intent(in) :: kstart
        !! Starting index for the Lanczos factorization (default 1).
        integer, optional, intent(in) :: kend
        !! Final index for the Lanczos factorization (default 1).
        real(dp), optional, intent(in) :: tol
        !! Tolerance to determine whether invariant subspaces have been computed or not.

        ! Internal variables.
        integer :: k_start, k_end
        real(dp) :: tolerance
        complex(dp) :: alpha, beta
        integer :: k, kdim

        info = 0

        ! Krylov subspace dimension.
        kdim = size(U) - 1

        ! Deals with the optional args.
        k_start = optval(kstart, 1)
        k_end   = optval(kend, kdim)
        tolerance = optval(tol, atol_dp)

        ! Lanczos bidiagonalization.
        lanczos : do k = k_start, k_end
            ! Transpose matrix-vector product.
            call A%rmatvec(U(k), V(k))

            ! Full re-orthogonalization of the right Krylov subspace.
            if (k > 1) then
                call double_gram_schmidt_step(V(k), V(:k-1), info, if_chk_orthonormal=.false.)
                call check_info(info, 'double_gram_schmidt_step', module=this_module, &
                                    & procedure='lanczos_bidiagonalization_cdp, right basis')
            end if

            ! Normalization step.
            alpha = V(k)%norm() ; B(k, k) = alpha
            if (abs(alpha) > tolerance) then
                call V(k)%scal(one_cdp/alpha)
            else
                info = k
                exit lanczos
            endif

            ! Matrix-vector product.
            call A%matvec(V(k), U(k+1))

            ! Full re-orthogonalization of the left Krylov subspace.
            call double_gram_schmidt_step(U(k+1), U(:k), info, if_chk_orthonormal=.false.)
            call check_info(info, 'double_gram_schmidt_step', module=this_module, &
                                    & procedure='lanczos_bidiagonalization_cdp, left basis')

            ! Normalization step
            beta = U(k+1)%norm() ; B(k+1, k) = beta
            if (abs(beta) > tolerance) then
                call U(k+1)%scal(one_cdp / beta)
            else
                info = k
                exit lanczos
            endif

        enddo lanczos

        return
    end subroutine lanczos_bidiagonalization_cdp



    !----------------------------------------------
    !-----     LANCZOS TRIDIAGONALIZATION     -----
    !----------------------------------------------
    
    subroutine lanczos_tridiagonalization_rsp(A, X, T, info, kstart, kend, tol)
        class(abstract_sym_linop_rsp), intent(in) :: A
        class(abstract_vector_rsp), intent(inout) :: X(:)
        real(sp), intent(inout) :: T(:, :)
        integer, intent(out) :: info
        integer, optional, intent(in) :: kstart
        integer, optional, intent(in) :: kend
        real(sp), optional, intent(in) :: tol

        ! Internal variables.
        integer :: k_start, k_end
        real(sp) :: tolerance
        real(sp) :: beta
        integer :: k, kdim

        ! Deal with optional args.
        kdim = size(X) - 1
        k_start = optval(kstart, 1)
        k_end = optval(kend, kdim)
        tolerance = optval(tol, atol_sp)
        info = 0

        ! Lanczos tridiagonalization.
        lanczos: do k = k_start, k_end
            ! Matrix-vector product.
            call A%matvec(X(k), X(k+1))
            ! Update tridiagonal matrix.
            call update_tridiag_matrix_rsp(T, X, k)
            beta = X(k+1)%norm() ; T(k+1, k) = beta

            ! Exit Lanczos loop if needed.
            if (beta < tolerance) then
                ! Dimension of the computed invariant subspace.
                info = k
                ! Exit the Lanczos iteration.
                exit lanczos
            else
                ! Normalize the new Krylov vector.
                call X(k+1)%scal(one_rsp / beta)
            endif
        enddo lanczos

        return
    end subroutine lanczos_tridiagonalization_rsp

    subroutine update_tridiag_matrix_rsp(T, X, k)
        integer, intent(in) :: k
        real(sp), intent(inout) :: T(:, :)
        class(abstract_vector_rsp), intent(inout) :: X(:)

        ! Internal variables.
        integer :: i, info

        info = 0

        ! Orthogonalize residual w.r.t. previously computed Krylov vectors to obtain coefficients in tridiag. matrix
        do i = max(1, k-1), k
            T(i, k) = X(i)%dot(X(k+1)) ; call X(k+1)%axpby(one_rsp, X(i), -T(i, k))
        enddo

        ! Full re-orthogonalization against existing basis
        call double_gram_schmidt_step(X(k+1), X(:k), info, if_chk_orthonormal=.false.)
        call check_info(info, 'orthogonalize_against_basis_p1', module=this_module, procedure='update_tridiag_matrix_rsp')

        return
    end subroutine update_tridiag_matrix_rsp

    subroutine lanczos_tridiagonalization_rdp(A, X, T, info, kstart, kend, tol)
        class(abstract_sym_linop_rdp), intent(in) :: A
        class(abstract_vector_rdp), intent(inout) :: X(:)
        real(dp), intent(inout) :: T(:, :)
        integer, intent(out) :: info
        integer, optional, intent(in) :: kstart
        integer, optional, intent(in) :: kend
        real(dp), optional, intent(in) :: tol

        ! Internal variables.
        integer :: k_start, k_end
        real(dp) :: tolerance
        real(dp) :: beta
        integer :: k, kdim

        ! Deal with optional args.
        kdim = size(X) - 1
        k_start = optval(kstart, 1)
        k_end = optval(kend, kdim)
        tolerance = optval(tol, atol_dp)
        info = 0

        ! Lanczos tridiagonalization.
        lanczos: do k = k_start, k_end
            ! Matrix-vector product.
            call A%matvec(X(k), X(k+1))
            ! Update tridiagonal matrix.
            call update_tridiag_matrix_rdp(T, X, k)
            beta = X(k+1)%norm() ; T(k+1, k) = beta

            ! Exit Lanczos loop if needed.
            if (beta < tolerance) then
                ! Dimension of the computed invariant subspace.
                info = k
                ! Exit the Lanczos iteration.
                exit lanczos
            else
                ! Normalize the new Krylov vector.
                call X(k+1)%scal(one_rdp / beta)
            endif
        enddo lanczos

        return
    end subroutine lanczos_tridiagonalization_rdp

    subroutine update_tridiag_matrix_rdp(T, X, k)
        integer, intent(in) :: k
        real(dp), intent(inout) :: T(:, :)
        class(abstract_vector_rdp), intent(inout) :: X(:)

        ! Internal variables.
        integer :: i, info

        info = 0

        ! Orthogonalize residual w.r.t. previously computed Krylov vectors to obtain coefficients in tridiag. matrix
        do i = max(1, k-1), k
            T(i, k) = X(i)%dot(X(k+1)) ; call X(k+1)%axpby(one_rdp, X(i), -T(i, k))
        enddo

        ! Full re-orthogonalization against existing basis
        call double_gram_schmidt_step(X(k+1), X(:k), info, if_chk_orthonormal=.false.)
        call check_info(info, 'orthogonalize_against_basis_p1', module=this_module, procedure='update_tridiag_matrix_rdp')

        return
    end subroutine update_tridiag_matrix_rdp

    subroutine lanczos_tridiagonalization_csp(A, X, T, info, kstart, kend, tol)
        class(abstract_hermitian_linop_csp), intent(in) :: A
        class(abstract_vector_csp), intent(inout) :: X(:)
        complex(sp), intent(inout) :: T(:, :)
        integer, intent(out) :: info
        integer, optional, intent(in) :: kstart
        integer, optional, intent(in) :: kend
        real(sp), optional, intent(in) :: tol

        ! Internal variables.
        integer :: k_start, k_end
        real(sp) :: tolerance
        real(sp) :: beta
        integer :: k, kdim

        ! Deal with optional args.
        kdim = size(X) - 1
        k_start = optval(kstart, 1)
        k_end = optval(kend, kdim)
        tolerance = optval(tol, atol_sp)
        info = 0

        ! Lanczos tridiagonalization.
        lanczos: do k = k_start, k_end
            ! Matrix-vector product.
            call A%matvec(X(k), X(k+1))
            ! Update tridiagonal matrix.
            call update_tridiag_matrix_csp(T, X, k)
            beta = X(k+1)%norm() ; T(k+1, k) = beta

            ! Exit Lanczos loop if needed.
            if (beta < tolerance) then
                ! Dimension of the computed invariant subspace.
                info = k
                ! Exit the Lanczos iteration.
                exit lanczos
            else
                ! Normalize the new Krylov vector.
                call X(k+1)%scal(one_csp / beta)
            endif
        enddo lanczos

        return
    end subroutine lanczos_tridiagonalization_csp

    subroutine update_tridiag_matrix_csp(T, X, k)
        integer, intent(in) :: k
        complex(sp), intent(inout) :: T(:, :)
        class(abstract_vector_csp), intent(inout) :: X(:)

        ! Internal variables.
        integer :: i, info

        info = 0

        ! Orthogonalize residual w.r.t. previously computed Krylov vectors to obtain coefficients in tridiag. matrix
        do i = max(1, k-1), k
            T(i, k) = X(i)%dot(X(k+1)) ; call X(k+1)%axpby(one_csp, X(i), -T(i, k))
        enddo

        ! Full re-orthogonalization against existing basis
        call double_gram_schmidt_step(X(k+1), X(:k), info, if_chk_orthonormal=.false.)
        call check_info(info, 'orthogonalize_against_basis_p1', module=this_module, procedure='update_tridiag_matrix_csp')

        return
    end subroutine update_tridiag_matrix_csp

    subroutine lanczos_tridiagonalization_cdp(A, X, T, info, kstart, kend, tol)
        class(abstract_hermitian_linop_cdp), intent(in) :: A
        class(abstract_vector_cdp), intent(inout) :: X(:)
        complex(dp), intent(inout) :: T(:, :)
        integer, intent(out) :: info
        integer, optional, intent(in) :: kstart
        integer, optional, intent(in) :: kend
        real(dp), optional, intent(in) :: tol

        ! Internal variables.
        integer :: k_start, k_end
        real(dp) :: tolerance
        real(dp) :: beta
        integer :: k, kdim

        ! Deal with optional args.
        kdim = size(X) - 1
        k_start = optval(kstart, 1)
        k_end = optval(kend, kdim)
        tolerance = optval(tol, atol_dp)
        info = 0

        ! Lanczos tridiagonalization.
        lanczos: do k = k_start, k_end
            ! Matrix-vector product.
            call A%matvec(X(k), X(k+1))
            ! Update tridiagonal matrix.
            call update_tridiag_matrix_cdp(T, X, k)
            beta = X(k+1)%norm() ; T(k+1, k) = beta

            ! Exit Lanczos loop if needed.
            if (beta < tolerance) then
                ! Dimension of the computed invariant subspace.
                info = k
                ! Exit the Lanczos iteration.
                exit lanczos
            else
                ! Normalize the new Krylov vector.
                call X(k+1)%scal(one_cdp / beta)
            endif
        enddo lanczos

        return
    end subroutine lanczos_tridiagonalization_cdp

    subroutine update_tridiag_matrix_cdp(T, X, k)
        integer, intent(in) :: k
        complex(dp), intent(inout) :: T(:, :)
        class(abstract_vector_cdp), intent(inout) :: X(:)

        ! Internal variables.
        integer :: i, info

        info = 0

        ! Orthogonalize residual w.r.t. previously computed Krylov vectors to obtain coefficients in tridiag. matrix
        do i = max(1, k-1), k
            T(i, k) = X(i)%dot(X(k+1)) ; call X(k+1)%axpby(one_cdp, X(i), -T(i, k))
        enddo

        ! Full re-orthogonalization against existing basis
        call double_gram_schmidt_step(X(k+1), X(:k), info, if_chk_orthonormal=.false.)
        call check_info(info, 'orthogonalize_against_basis_p1', module=this_module, procedure='update_tridiag_matrix_cdp')

        return
    end subroutine update_tridiag_matrix_cdp


    !----------------------------------------
    !-----     KRYLOV-SCHUR RESTART     -----
    !----------------------------------------

    subroutine krylov_schur_rsp(n, X, H, select_eigs)
        integer, intent(out) :: n
        !! Number eigenvalues that have been moved to the upper
        !! left block of the Schur factorization of `H`.
        class(abstract_vector_rsp), intent(inout) :: X(:)
        !! Krylov basis.
        real(sp), intent(inout) :: H(:, :)
        !! Upper Hessenberg matrix.
        procedure(eigvals_select_sp) :: select_eigs
        !! Procedure to select the eigenvalues to move in the upper left-block.

        !--------------------------------------
        !-----     Internal variables     -----
        !--------------------------------------

        integer :: i, kdim
        
        ! Schur-related.
        real(sp) :: Z(size(H, 2), size(H, 2))
        complex(sp) :: eigvals(size(H, 2))
        logical :: selected(size(H, 2))
       
        ! Krylov subspace dimension.
        kdim = size(X)-1

        ! Allocate and initializes variables.
        eigvals = zero_rsp ; Z = zero_rsp
        selected = .false.

        ! Schur decomposition of the Hessenberg matrix.
        call schur(H(:size(H, 2), :), Z, eigvals)

        ! Eigenvalue selection of the upper left block.
        selected = select_eigs(eigvals) ; n = count(selected)

        ! Re-order the Schur decomposition and Schur basis.
        call ordschur(H(:kdim, :), Z, selected)

        ! Update the Hessenberg matrix and Krylov basis.
        block
        real(sp) :: b(size(H, 2))
        class(abstract_vector_rsp), allocatable :: Xwrk(:)
        
        ! Update the Krylov basis.
        call linear_combination(Xwrk, X(:size(H, 2)), Z(:, :n))
        call copy(X(:n), Xwrk(:n))
        call copy(X(n+1), X(kdim+1))
        call zero_basis(X(n+2:))

        ! Update the Hessenberg matrix.
        b = matmul(H(kdim+1, :), Z)
        H(n+1, :) = b
        H(n+2:, :) = zero_rsp
        H(:, n+1:) = zero_rsp
        end block

        return
    end subroutine krylov_schur_rsp

    subroutine krylov_schur_rdp(n, X, H, select_eigs)
        integer, intent(out) :: n
        !! Number eigenvalues that have been moved to the upper
        !! left block of the Schur factorization of `H`.
        class(abstract_vector_rdp), intent(inout) :: X(:)
        !! Krylov basis.
        real(dp), intent(inout) :: H(:, :)
        !! Upper Hessenberg matrix.
        procedure(eigvals_select_dp) :: select_eigs
        !! Procedure to select the eigenvalues to move in the upper left-block.

        !--------------------------------------
        !-----     Internal variables     -----
        !--------------------------------------

        integer :: i, kdim
        
        ! Schur-related.
        real(dp) :: Z(size(H, 2), size(H, 2))
        complex(dp) :: eigvals(size(H, 2))
        logical :: selected(size(H, 2))
       
        ! Krylov subspace dimension.
        kdim = size(X)-1

        ! Allocate and initializes variables.
        eigvals = zero_rdp ; Z = zero_rdp
        selected = .false.

        ! Schur decomposition of the Hessenberg matrix.
        call schur(H(:size(H, 2), :), Z, eigvals)

        ! Eigenvalue selection of the upper left block.
        selected = select_eigs(eigvals) ; n = count(selected)

        ! Re-order the Schur decomposition and Schur basis.
        call ordschur(H(:kdim, :), Z, selected)

        ! Update the Hessenberg matrix and Krylov basis.
        block
        real(dp) :: b(size(H, 2))
        class(abstract_vector_rdp), allocatable :: Xwrk(:)
        
        ! Update the Krylov basis.
        call linear_combination(Xwrk, X(:size(H, 2)), Z(:, :n))
        call copy(X(:n), Xwrk(:n))
        call copy(X(n+1), X(kdim+1))
        call zero_basis(X(n+2:))

        ! Update the Hessenberg matrix.
        b = matmul(H(kdim+1, :), Z)
        H(n+1, :) = b
        H(n+2:, :) = zero_rdp
        H(:, n+1:) = zero_rdp
        end block

        return
    end subroutine krylov_schur_rdp

    subroutine krylov_schur_csp(n, X, H, select_eigs)
        integer, intent(out) :: n
        !! Number eigenvalues that have been moved to the upper
        !! left block of the Schur factorization of `H`.
        class(abstract_vector_csp), intent(inout) :: X(:)
        !! Krylov basis.
        complex(sp), intent(inout) :: H(:, :)
        !! Upper Hessenberg matrix.
        procedure(eigvals_select_sp) :: select_eigs
        !! Procedure to select the eigenvalues to move in the upper left-block.

        !--------------------------------------
        !-----     Internal variables     -----
        !--------------------------------------

        integer :: i, kdim
        
        ! Schur-related.
        complex(sp) :: Z(size(H, 2), size(H, 2))
        complex(sp) :: eigvals(size(H, 2))
        logical :: selected(size(H, 2))
       
        ! Krylov subspace dimension.
        kdim = size(X)-1

        ! Allocate and initializes variables.
        eigvals = zero_csp ; Z = zero_csp
        selected = .false.

        ! Schur decomposition of the Hessenberg matrix.
        call schur(H(:size(H, 2), :), Z, eigvals)

        ! Eigenvalue selection of the upper left block.
        selected = select_eigs(eigvals) ; n = count(selected)

        ! Re-order the Schur decomposition and Schur basis.
        call ordschur(H(:kdim, :), Z, selected)

        ! Update the Hessenberg matrix and Krylov basis.
        block
        complex(sp) :: b(size(H, 2))
        class(abstract_vector_csp), allocatable :: Xwrk(:)
        
        ! Update the Krylov basis.
        call linear_combination(Xwrk, X(:size(H, 2)), Z(:, :n))
        call copy(X(:n), Xwrk(:n))
        call copy(X(n+1), X(kdim+1))
        call zero_basis(X(n+2:))

        ! Update the Hessenberg matrix.
        b = matmul(H(kdim+1, :), Z)
        H(n+1, :) = b
        H(n+2:, :) = zero_csp
        H(:, n+1:) = zero_csp
        end block

        return
    end subroutine krylov_schur_csp

    subroutine krylov_schur_cdp(n, X, H, select_eigs)
        integer, intent(out) :: n
        !! Number eigenvalues that have been moved to the upper
        !! left block of the Schur factorization of `H`.
        class(abstract_vector_cdp), intent(inout) :: X(:)
        !! Krylov basis.
        complex(dp), intent(inout) :: H(:, :)
        !! Upper Hessenberg matrix.
        procedure(eigvals_select_dp) :: select_eigs
        !! Procedure to select the eigenvalues to move in the upper left-block.

        !--------------------------------------
        !-----     Internal variables     -----
        !--------------------------------------

        integer :: i, kdim
        
        ! Schur-related.
        complex(dp) :: Z(size(H, 2), size(H, 2))
        complex(dp) :: eigvals(size(H, 2))
        logical :: selected(size(H, 2))
       
        ! Krylov subspace dimension.
        kdim = size(X)-1

        ! Allocate and initializes variables.
        eigvals = zero_cdp ; Z = zero_cdp
        selected = .false.

        ! Schur decomposition of the Hessenberg matrix.
        call schur(H(:size(H, 2), :), Z, eigvals)

        ! Eigenvalue selection of the upper left block.
        selected = select_eigs(eigvals) ; n = count(selected)

        ! Re-order the Schur decomposition and Schur basis.
        call ordschur(H(:kdim, :), Z, selected)

        ! Update the Hessenberg matrix and Krylov basis.
        block
        complex(dp) :: b(size(H, 2))
        class(abstract_vector_cdp), allocatable :: Xwrk(:)
        
        ! Update the Krylov basis.
        call linear_combination(Xwrk, X(:size(H, 2)), Z(:, :n))
        call copy(X(:n), Xwrk(:n))
        call copy(X(n+1), X(kdim+1))
        call zero_basis(X(n+2:))

        ! Update the Hessenberg matrix.
        b = matmul(H(kdim+1, :), Z)
        H(n+1, :) = b
        H(n+2:, :) = zero_cdp
        H(:, n+1:) = zero_cdp
        end block

        return
    end subroutine krylov_schur_cdp


end module LightKrylov_BaseKrylov
