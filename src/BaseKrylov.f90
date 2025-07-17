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
    use stdlib_linalg, only: eye, schur, norm, mnorm

    !-------------------------------
    !-----     LightKrylov     -----
    !-------------------------------
    use LightKrylov_Constants
    use LightKrylov_Logger, only: log_warning, log_error, log_message, log_information, &
    &                             stop_error, check_info
    use LightKrylov_Timing, only: timer => global_lightkrylov_timer, time_lightkrylov
    use LightKrylov_Utils
    use LightKrylov_AbstractVectors
    use LightKrylov_AbstractLinops

    implicit none
    private
    
    character(len=*), parameter :: this_module      = 'LK_BKrylov'
    character(len=*), parameter :: this_module_long = 'LightKrylov_BaseKrylov'

    !----- Krylov processes ------
    public :: arnoldi
    public :: bidiagonalization
    public :: lanczos

    !----- Utility functions ------
    public :: qr
    public :: double_gram_schmidt_step
    public :: is_orthonormal
    public :: orthonormalize_basis
    public :: orthogonalize_against_basis
    public :: permcols, invperm

    public :: initialize_krylov_subspace
    public :: krylov_schur

    !------------------------------------
    !-----                          -----
    !-----     KRYLOV PROCESSES     -----
    !-----                          -----
    !------------------------------------

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
        !!  - `A` : Linear operator derived from one the base types provided by the `AbstractLinops`
        !!          module. The operator needs to be square, i.e. the dimension of its domain and
        !!          co-domain is the same. It is an `intent(inout)` argument.
        !!
        !!  - `X` : Array of types derived from one the base types provided by the `AbstractVectors`
        !!          module. It needs to be consistent with the type of `A`. On exit, it contains the
        !!          the computed Krylov vectors. The first entry `X(1)` is the starting vector for
        !!          the Arnoldi factorization. Additionally, the maximum number of Arnoldi steps
        !!          is equal to `size(X) - 1`. It is an `intent(inout)` argument.
        !!
        !!  -`H` : `real` or `complex` rank-2 array. On exit, it contains the \( (k+1) \times k\)
        !!          upper Hessenberg matrix computed from the Arnoldi factorization. It is an
        !!          `intent(inout)` argument.
        !!
        !!  -`info` :   `integer` variable. It is the `LightKrylov` information flag. On exit, if
        !!              `info` > 0, the Arnoldi factorization experienced a lucky breakdown. 
        !!              The array of Krylov vectors `X` spans an \(A\)-invariant subpsace of
        !!              dimension `info`.
        !!
        !!  - `kstart` (*optional*) :   `integer` value determining the index of the first Arnoldi
        !!                              step to be computed. By default, `kstart = 1`.
        !!
        !!  - `kend` (*optional*)   :   `integer` value determining the index of the last Arnoldi step
        !!                              to be computed. By default, `kend = size(X) - 1`.
        !!
        !!  - `tol` (*optional*)    :   Numerical tolerance below which a subspace is considered
        !!                              to be \( A \)-invariant. By default `tol = atol_sp` or
        !!                              `tol = atol_rp` depending on the kind of `A`.
        !!
        !!  - `transpose` (*optional*)  :   `logical` flag determining whether the Arnoldi factorization
        !!                                  is applied to \( A \) or \( A^H \). Default `transpose = .false.`
        !!
        !!  - `blksize` (*optional*)    :   `integer` value determining the dimension of a block for the
        !!                                  block Arnoldi factorization. Default is `blksize=1`.
        module subroutine arnoldi_rsp(A, X, H, info, kstart, kend, tol, transpose, blksize)
            class(abstract_linop_rsp), intent(inout) :: A
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
        end subroutine
        module subroutine arnoldi_rdp(A, X, H, info, kstart, kend, tol, transpose, blksize)
            class(abstract_linop_rdp), intent(inout) :: A
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
        end subroutine
        module subroutine arnoldi_csp(A, X, H, info, kstart, kend, tol, transpose, blksize)
            class(abstract_linop_csp), intent(inout) :: A
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
        end subroutine
        module subroutine arnoldi_cdp(A, X, H, info, kstart, kend, tol, transpose, blksize)
            class(abstract_linop_cdp), intent(inout) :: A
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
        end subroutine
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
        !!  - `A`   :   Symmetric or Hermitian linear operator derived from one the base types 
        !!              provided by the `AbstractLinops` module. It is an `intent(inout)` argument.
        !!
        !!  - `X`   :   Array of types derived from one the base types provided by the `AbstractVectors`
        !!              module. It needs to be consistent with the type of `A`. On exit, it contains the
        !!              the computed Krylov vectors. The first entry `X(1)` is the starting vector for
        !!              the Lanczos factorization. Additionally, the maximum number of Lanczos steps
        !!              is equal to `size(X) - 1`. It is an `intent(inout)` argument.
        !!
        !!  - `T`   :   `real` or `complex` rank-2 array. On exit, it contains the \( (k+1) \times k\)
        !!              symmetric tridiagonal matrix computed from the Arnoldi factorization. It is an
        !!              `intent(inout)` argument.
        !!
        !!  - `info`    :   `integer` variable. It is the `LightKrylov` information flag. On exit, if
        !!                  `info` > 0, the Lanczos factorization experienced a lucky breakdown. 
        !!                  The array of Krylov vectors `X` spans an \(A\)-invariant subpsace of
        !!                  dimension `info`.
        !!
        !!  - `kstart` (*optional*) :   `integer` value determining the index of the first Lanczos
        !!                              step to be computed. By default, `kstart = 1`.
        !!
        !!  - `kend` (*optional*)   :   `integer` value determining the index of the last Lanczos step
        !!                              to be computed. By default, `kend = size(X) - 1`.
        !!
        !!  - `tol` (*optional*)    :   Numerical tolerance below which a subspace is considered
        !!                              to be \( A \)-invariant. By default `tol = atol_sp` or
        !!                              `tol = atol_rp` depending on the kind of `A`.
        module subroutine lanczos_tridiagonalization_rsp(A, X, T, info, kstart, kend, tol)
            class(abstract_sym_linop_rsp), intent(inout) :: A
            class(abstract_vector_rsp), intent(inout) :: X(:)
            real(sp), intent(inout) :: T(:, :)
            integer, intent(out) :: info
            integer, optional, intent(in) :: kstart
            integer, optional, intent(in) :: kend
            real(sp), optional, intent(in) :: tol
        end subroutine  
        module subroutine lanczos_tridiagonalization_rdp(A, X, T, info, kstart, kend, tol)
            class(abstract_sym_linop_rdp), intent(inout) :: A
            class(abstract_vector_rdp), intent(inout) :: X(:)
            real(dp), intent(inout) :: T(:, :)
            integer, intent(out) :: info
            integer, optional, intent(in) :: kstart
            integer, optional, intent(in) :: kend
            real(dp), optional, intent(in) :: tol
        end subroutine  
        module subroutine lanczos_tridiagonalization_csp(A, X, T, info, kstart, kend, tol)
            class(abstract_hermitian_linop_csp), intent(inout) :: A
            class(abstract_vector_csp), intent(inout) :: X(:)
            complex(sp), intent(inout) :: T(:, :)
            integer, intent(out) :: info
            integer, optional, intent(in) :: kstart
            integer, optional, intent(in) :: kend
            real(sp), optional, intent(in) :: tol
        end subroutine  
        module subroutine lanczos_tridiagonalization_cdp(A, X, T, info, kstart, kend, tol)
            class(abstract_hermitian_linop_cdp), intent(inout) :: A
            class(abstract_vector_cdp), intent(inout) :: X(:)
            complex(dp), intent(inout) :: T(:, :)
            integer, intent(out) :: info
            integer, optional, intent(in) :: kstart
            integer, optional, intent(in) :: kend
            real(dp), optional, intent(in) :: tol
        end subroutine  
    end interface

    interface bidiagonalization
        !!  ### Description
        !!
        !!  Given a general linear operator \( A \), find matrices \( U \), \( V \) and
        !!  \( B \) such that
        !!
        !!  \[
        !!      \begin{aligned}
        !!          AV_k & = U_k T_k + t_{k+1, k} v_{k+1} e_k^T, \\
        !!          A^H U_k & = V_k T_k^T + t_{k+1, k} u_{k+1} e_k^T
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
        !!  - `A`   :   Linear operator derived from one the base types provided by the 
        !!              `AbstractLinops` module. It is an `intent(inout)` argument.
        !!
        !!  - `U`   :   Array of types derived from one the base types provided by the `AbstractVectors`
        !!              module. It needs to be consistent with the type of `A`. On exit, it contains the
        !!              the computed Krylov vectors for the column span of `A`. The first entry `U(1)` 
        !!              is the starting vector for the Lanczos factorization. Additionally, the 
        !!              maximum number of Lanczos steps is equal to `size(X) - 1`. 
        !!              It is an `intent(inout)` argument.
        !!
        !!  - `V`   :   Array of types derived from one the base types provided by the `AbstractVectors`
        !!              module. It needs to be consistent with the type of `A`. On exit, it contains the
        !!              the computed Krylov vectors for the row span of `A`. It is an `intent(inout)` 
        !!              argument.
        !!
        !!  - `B`   :   `real` or `complex` rank-2 array. On exit, it contains the \( (k+1) \times k\)
        !!              bidiagonal matrix computed from the Lanczos factorization. It is an
        !!              `intent(inout)` argument.
        !!
        !!  - `info`    :   `integer` variable. It is the `LightKrylov` information flag. On exit, if
        !!                  `info` > 0, the Lanczos factorization experienced a lucky breakdown. 
        !!
        !!  - `kstart` (*optional*) :   `integer` value determining the index of the first Lanczos
        !!                              step to be computed. By default, `kstart = 1`.
        !!
        !!  - `kend` (*optional*)   :   `integer` value determining the index of the last Lanczos step
        !!                              to be computed. By default, `kend = size(X) - 1`.
        !!
        !!  - `tol` (*optional*)    :   Numerical tolerance below which a subspace is considered
        !!                              to be \( A \)-invariant. By default `tol = atol_sp` or
        !!                              `tol = atol_rp` depending on the kind of `A`.
        module subroutine lanczos_bidiagonalization_rsp(A, U, V, B, info, kstart, kend, tol)
            class(abstract_linop_rsp), intent(inout) :: A
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
        end subroutine
        module subroutine lanczos_bidiagonalization_rdp(A, U, V, B, info, kstart, kend, tol)
            class(abstract_linop_rdp), intent(inout) :: A
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
        end subroutine
        module subroutine lanczos_bidiagonalization_csp(A, U, V, B, info, kstart, kend, tol)
            class(abstract_linop_csp), intent(inout) :: A
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
        end subroutine
        module subroutine lanczos_bidiagonalization_cdp(A, U, V, B, info, kstart, kend, tol)
            class(abstract_linop_cdp), intent(inout) :: A
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
        end subroutine
    end interface

    !-------------------------------------
    !-----                           -----
    !-----     UTILITY FUNCTIONS     -----
    !-----                           -----
    !-------------------------------------

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
        !!  - `Q`   :   Array of types derived from one of the base types provided in the
        !!              `AbstractVectors` module. On entry, it contains the original array.
        !!              On exit, it is overwritten by the orthogonal basis for its span.
        !!              It is an `intent(inout)` argument.
        !!
        !!  - `R`   :   `real` or `complex` rank-2 array. On exit, its contains the upper triangular
        !!              matrix resulting from the QR factorization. It is an `intent(out)` argument.
        !!
        !!  - `perm` (*optional*)   :   Rank-1 array of `integer` corresponding to the indices of
        !!                              permuted columns. If `perm` is absent, the naive QR factorization
        !!                              is being computed.
        !!
        !!  - `info`    :   `integer` information flag.
        !!
        !!  - `tol` (*optional*)    :   Numerical tolerance to determine whether two vectors are colinear
        !!                              or not. Default `tol = atol_sp` or `tol = atol_dp`.
        module subroutine qr_no_pivoting_rsp(Q, R, info, tol)
            class(abstract_vector_rsp), intent(inout) :: Q(:)
            !! Array of `abstract_vector` to be orthogonalized.
            real(sp), intent(out) :: R(:, :)
            !! Upper triangular matrix \(\mathbf{R}\) resulting from the QR factorization.
            integer, intent(out) :: info
            !! Information flag.
            real(sp), optional, intent(in) :: tol
        end subroutine

        module subroutine qr_with_pivoting_rsp(Q, R, perm, info, tol)
            class(abstract_vector_rsp), intent(inout) :: Q(:)
            !! Array of `abstract_vector` to be orthogonalized.
            real(sp), intent(out) :: R(:, :)
            !! Upper triangular matrix resulting from the QR factorization.
            integer, intent(out) :: perm(size(Q))
            !! Permutation matrix.
            integer, intent(out) :: info
            !! Information flag.
            real(sp), optional, intent(in) :: tol
        end subroutine 
        module subroutine qr_no_pivoting_rdp(Q, R, info, tol)
            class(abstract_vector_rdp), intent(inout) :: Q(:)
            !! Array of `abstract_vector` to be orthogonalized.
            real(dp), intent(out) :: R(:, :)
            !! Upper triangular matrix \(\mathbf{R}\) resulting from the QR factorization.
            integer, intent(out) :: info
            !! Information flag.
            real(dp), optional, intent(in) :: tol
        end subroutine

        module subroutine qr_with_pivoting_rdp(Q, R, perm, info, tol)
            class(abstract_vector_rdp), intent(inout) :: Q(:)
            !! Array of `abstract_vector` to be orthogonalized.
            real(dp), intent(out) :: R(:, :)
            !! Upper triangular matrix resulting from the QR factorization.
            integer, intent(out) :: perm(size(Q))
            !! Permutation matrix.
            integer, intent(out) :: info
            !! Information flag.
            real(dp), optional, intent(in) :: tol
        end subroutine 
        module subroutine qr_no_pivoting_csp(Q, R, info, tol)
            class(abstract_vector_csp), intent(inout) :: Q(:)
            !! Array of `abstract_vector` to be orthogonalized.
            complex(sp), intent(out) :: R(:, :)
            !! Upper triangular matrix \(\mathbf{R}\) resulting from the QR factorization.
            integer, intent(out) :: info
            !! Information flag.
            real(sp), optional, intent(in) :: tol
        end subroutine

        module subroutine qr_with_pivoting_csp(Q, R, perm, info, tol)
            class(abstract_vector_csp), intent(inout) :: Q(:)
            !! Array of `abstract_vector` to be orthogonalized.
            complex(sp), intent(out) :: R(:, :)
            !! Upper triangular matrix resulting from the QR factorization.
            integer, intent(out) :: perm(size(Q))
            !! Permutation matrix.
            integer, intent(out) :: info
            !! Information flag.
            real(sp), optional, intent(in) :: tol
        end subroutine 
        module subroutine qr_no_pivoting_cdp(Q, R, info, tol)
            class(abstract_vector_cdp), intent(inout) :: Q(:)
            !! Array of `abstract_vector` to be orthogonalized.
            complex(dp), intent(out) :: R(:, :)
            !! Upper triangular matrix \(\mathbf{R}\) resulting from the QR factorization.
            integer, intent(out) :: info
            !! Information flag.
            real(dp), optional, intent(in) :: tol
        end subroutine

        module subroutine qr_with_pivoting_cdp(Q, R, perm, info, tol)
            class(abstract_vector_cdp), intent(inout) :: Q(:)
            !! Array of `abstract_vector` to be orthogonalized.
            complex(dp), intent(out) :: R(:, :)
            !! Upper triangular matrix resulting from the QR factorization.
            integer, intent(out) :: perm(size(Q))
            !! Permutation matrix.
            integer, intent(out) :: info
            !! Information flag.
            real(dp), optional, intent(in) :: tol
        end subroutine 
    end interface

    interface permcols
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
        !!      call permcols(X, perm)
        !!  ```
        !!
        !!  ### Arguments
        !!
        !!  - `X`   :   Array of vectors derived from the base types defined in the `AbstractVectors`
        !!              module. On entry, it is the original array. On exit, it contains the
        !!              column-permuted version computed in-place. It is an `intent(inout)` argument.
        !!
        !!  - `perm`    :   Rank-1 array of `integer` corresponding to the desired permutation vector.
        !!                  It is an `intent(in)` argument.
        module subroutine permcols_basis_rsp(Q, perm)
            class(abstract_vector_rsp), intent(inout) :: Q(:)
            !! Basis vectors to be permuted.
            integer, intent(in) :: perm(:)
        end subroutine
 
        module subroutine permcols_array_rsp(Q, perm)
            real(sp), intent(inout) :: Q(:, :)
            !! Basis vectors to be permuted.
            integer, intent(in) :: perm(:)
        end subroutine
        module subroutine permcols_basis_rdp(Q, perm)
            class(abstract_vector_rdp), intent(inout) :: Q(:)
            !! Basis vectors to be permuted.
            integer, intent(in) :: perm(:)
        end subroutine
 
        module subroutine permcols_array_rdp(Q, perm)
            real(dp), intent(inout) :: Q(:, :)
            !! Basis vectors to be permuted.
            integer, intent(in) :: perm(:)
        end subroutine
        module subroutine permcols_basis_csp(Q, perm)
            class(abstract_vector_csp), intent(inout) :: Q(:)
            !! Basis vectors to be permuted.
            integer, intent(in) :: perm(:)
        end subroutine
 
        module subroutine permcols_array_csp(Q, perm)
            complex(sp), intent(inout) :: Q(:, :)
            !! Basis vectors to be permuted.
            integer, intent(in) :: perm(:)
        end subroutine
        module subroutine permcols_basis_cdp(Q, perm)
            class(abstract_vector_cdp), intent(inout) :: Q(:)
            !! Basis vectors to be permuted.
            integer, intent(in) :: perm(:)
        end subroutine
 
        module subroutine permcols_array_cdp(Q, perm)
            complex(dp), intent(inout) :: Q(:, :)
            !! Basis vectors to be permuted.
            integer, intent(in) :: perm(:)
        end subroutine
    end interface

    interface 
        !!  ### Description
        !!
        !!  Given a permutation vector \( p \), this function computes the vector
        !!  representation of the inverse permutation matrix.
        !!
        !!  ### Syntax
        !!
        !!  ```fortran
        !!      inv_perm = invperm(perm)
        !!  ```
        !!
        !!  ### Arguments
        !!
        !!  - `perm`    :   Rank-1 array of `integer` corresponding to the desired permutation vector.
        !!                  It is an `intent(in)` argument.
        module function invperm(perm) result(inv_perm)
            integer, intent(in) :: perm(:)
            integer, allocatable :: inv_perm(:)
        end function
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
        !!  - `X`   :   Array of vectors that needs to be initialized. It is an `intent(inout)`
        !!              argument. Note that the first action in the subroutine is
        !!              `call zero_basis(X)`, effectively zeroing-out any data stored.
        !!
        !!  - `X0` (*optional*) :   Collection of vectors which will form the first few
        !!                          Krylov vectors. Note that `X0` need not be an orthonormal
        !!                          basis as this subroutine includes a `call qr(X0)`.
        module subroutine initialize_krylov_subspace_rsp(X, X0)
            class(abstract_vector_rsp), intent(inout) :: X(:)
            class(abstract_vector_rsp), optional, intent(in) :: X0(:)
        end subroutine
        module subroutine initialize_krylov_subspace_rdp(X, X0)
            class(abstract_vector_rdp), intent(inout) :: X(:)
            class(abstract_vector_rdp), optional, intent(in) :: X0(:)
        end subroutine
        module subroutine initialize_krylov_subspace_csp(X, X0)
            class(abstract_vector_csp), intent(inout) :: X(:)
            class(abstract_vector_csp), optional, intent(in) :: X0(:)
        end subroutine
        module subroutine initialize_krylov_subspace_cdp(X, X0)
            class(abstract_vector_cdp), intent(inout) :: X(:)
            class(abstract_vector_cdp), optional, intent(in) :: X0(:)
        end subroutine
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
        !!  - `X`   :   Array of derived types extended from the base types provided in the
        !!              `AbstractVectors` module.
        module function is_orthonormal_rsp(X) result(ortho)
            class(abstract_vector_rsp), intent(in) :: X(:)
            logical :: ortho
        end function
        module function is_orthonormal_rdp(X) result(ortho)
            class(abstract_vector_rdp), intent(in) :: X(:)
            logical :: ortho
        end function
        module function is_orthonormal_csp(X) result(ortho)
            class(abstract_vector_csp), intent(in) :: X(:)
            logical :: ortho
        end function
        module function is_orthonormal_cdp(X) result(ortho)
            class(abstract_vector_cdp), intent(in) :: X(:)
            logical :: ortho
        end function
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
        !!  - `X`   :   Array of `abstract_vector` to orthonormalize. Note that this process is done
        !!              in-place. It is an `intent(inout)` argument.
        module subroutine orthonormalize_basis_rsp(X)
            !! Orthonormalizes the `abstract_vector` basis `X`
            class(abstract_vector_rsp), intent(inout) :: X(:)
            !! Input `abstract_vector` basis to orthogonalize against
        end subroutine
        module subroutine orthonormalize_basis_rdp(X)
            !! Orthonormalizes the `abstract_vector` basis `X`
            class(abstract_vector_rdp), intent(inout) :: X(:)
            !! Input `abstract_vector` basis to orthogonalize against
        end subroutine
        module subroutine orthonormalize_basis_csp(X)
            !! Orthonormalizes the `abstract_vector` basis `X`
            class(abstract_vector_csp), intent(inout) :: X(:)
            !! Input `abstract_vector` basis to orthogonalize against
        end subroutine
        module subroutine orthonormalize_basis_cdp(X)
            !! Orthonormalizes the `abstract_vector` basis `X`
            class(abstract_vector_cdp), intent(inout) :: X(:)
            !! Input `abstract_vector` basis to orthogonalize against
        end subroutine
    end interface

    interface orthogonalize_against_basis
        module subroutine orthogonalize_vector_against_basis_rsp(y, X, info, if_chk_orthonormal, beta)
           !! Orthogonalizes the `abstract_vector` `y` against a basis `X` of `abstract_vector`.
            class(abstract_vector_rsp), intent(inout) :: y
            !! Input `abstract_vector` to orthogonalize
            class(abstract_vector_rsp), intent(in)    :: X(:)
            !! Input `abstract_vector` basis to orthogonalize against
            integer, intent(out) :: info
            !! Information flag.
            logical,                          optional, intent(in)    :: if_chk_orthonormal
            !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
            real(sp),                         optional, intent(out)   :: beta(:)
            !! Projection coefficients if requested
        end subroutine

        module subroutine orthogonalize_basis_against_basis_rsp(Y, X, info, if_chk_orthonormal, beta)
            !! Orthogonalizes the `abstract_vector` basis `Y` against a basis `X` of `abstract_vector`.
            class(abstract_vector_rsp), intent(inout) :: Y(:)
            !! Input `abstract_vector` basis to orthogonalize
            class(abstract_vector_rsp), intent(in)    :: X(:)
            !! Input `abstract_vector` basis to orthogonalize against
            integer, intent(out) :: info
            !! Information flag.
            logical,                          optional, intent(in)    :: if_chk_orthonormal
            !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
            real(sp),                         optional, intent(out)   :: beta(:,:)
            !! Projection coefficients if requested
        end subroutine
        module subroutine orthogonalize_vector_against_basis_rdp(y, X, info, if_chk_orthonormal, beta)
           !! Orthogonalizes the `abstract_vector` `y` against a basis `X` of `abstract_vector`.
            class(abstract_vector_rdp), intent(inout) :: y
            !! Input `abstract_vector` to orthogonalize
            class(abstract_vector_rdp), intent(in)    :: X(:)
            !! Input `abstract_vector` basis to orthogonalize against
            integer, intent(out) :: info
            !! Information flag.
            logical,                          optional, intent(in)    :: if_chk_orthonormal
            !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
            real(dp),                         optional, intent(out)   :: beta(:)
            !! Projection coefficients if requested
        end subroutine

        module subroutine orthogonalize_basis_against_basis_rdp(Y, X, info, if_chk_orthonormal, beta)
            !! Orthogonalizes the `abstract_vector` basis `Y` against a basis `X` of `abstract_vector`.
            class(abstract_vector_rdp), intent(inout) :: Y(:)
            !! Input `abstract_vector` basis to orthogonalize
            class(abstract_vector_rdp), intent(in)    :: X(:)
            !! Input `abstract_vector` basis to orthogonalize against
            integer, intent(out) :: info
            !! Information flag.
            logical,                          optional, intent(in)    :: if_chk_orthonormal
            !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
            real(dp),                         optional, intent(out)   :: beta(:,:)
            !! Projection coefficients if requested
        end subroutine
        module subroutine orthogonalize_vector_against_basis_csp(y, X, info, if_chk_orthonormal, beta)
           !! Orthogonalizes the `abstract_vector` `y` against a basis `X` of `abstract_vector`.
            class(abstract_vector_csp), intent(inout) :: y
            !! Input `abstract_vector` to orthogonalize
            class(abstract_vector_csp), intent(in)    :: X(:)
            !! Input `abstract_vector` basis to orthogonalize against
            integer, intent(out) :: info
            !! Information flag.
            logical,                          optional, intent(in)    :: if_chk_orthonormal
            !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
            complex(sp),                         optional, intent(out)   :: beta(:)
            !! Projection coefficients if requested
        end subroutine

        module subroutine orthogonalize_basis_against_basis_csp(Y, X, info, if_chk_orthonormal, beta)
            !! Orthogonalizes the `abstract_vector` basis `Y` against a basis `X` of `abstract_vector`.
            class(abstract_vector_csp), intent(inout) :: Y(:)
            !! Input `abstract_vector` basis to orthogonalize
            class(abstract_vector_csp), intent(in)    :: X(:)
            !! Input `abstract_vector` basis to orthogonalize against
            integer, intent(out) :: info
            !! Information flag.
            logical,                          optional, intent(in)    :: if_chk_orthonormal
            !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
            complex(sp),                         optional, intent(out)   :: beta(:,:)
            !! Projection coefficients if requested
        end subroutine
        module subroutine orthogonalize_vector_against_basis_cdp(y, X, info, if_chk_orthonormal, beta)
           !! Orthogonalizes the `abstract_vector` `y` against a basis `X` of `abstract_vector`.
            class(abstract_vector_cdp), intent(inout) :: y
            !! Input `abstract_vector` to orthogonalize
            class(abstract_vector_cdp), intent(in)    :: X(:)
            !! Input `abstract_vector` basis to orthogonalize against
            integer, intent(out) :: info
            !! Information flag.
            logical,                          optional, intent(in)    :: if_chk_orthonormal
            !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
            complex(dp),                         optional, intent(out)   :: beta(:)
            !! Projection coefficients if requested
        end subroutine

        module subroutine orthogonalize_basis_against_basis_cdp(Y, X, info, if_chk_orthonormal, beta)
            !! Orthogonalizes the `abstract_vector` basis `Y` against a basis `X` of `abstract_vector`.
            class(abstract_vector_cdp), intent(inout) :: Y(:)
            !! Input `abstract_vector` basis to orthogonalize
            class(abstract_vector_cdp), intent(in)    :: X(:)
            !! Input `abstract_vector` basis to orthogonalize against
            integer, intent(out) :: info
            !! Information flag.
            logical,                          optional, intent(in)    :: if_chk_orthonormal
            !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
            complex(dp),                         optional, intent(out)   :: beta(:,:)
            !! Projection coefficients if requested
        end subroutine
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
        !!  - `y`   :   `abstract_vector` (or array of `abstract_vector`) that needs to be
        !!              orthogonalize **in-place** against \( X \).
        !!
        !!  - `X`   :   Array of `abstract_vector` against which \( y \) needs to be orthogonalized.
        !!              Note the function assumes that \( X \) is an orthonormal set of vectors, i.e.
        !!              \( X^H X = I \). If it this is not the case, the result are meaningless.
        !!
        !!  - `info`    :   `integer` Information flag.
        !!
        !!  - `if_chk_orthonormal` (*optional*) :   `logical` flag (default `.true.`) to check
        !!                                          whether \( X \) is an orthonormal set of vectors or not. If the orthonormality
        !!                                          returns `.false.`, the function throws an error. Note that this check is however
        !!                                          computationally expensive and can be disable for the sake of performances.
        !!
        !!  - `beta` (*optional*)   :   `real` or `complex` array containing the coefficients \( \beta = X^H y \).
        module subroutine DGS_vector_against_basis_rsp(y, X, info, if_chk_orthonormal, beta)
          !! Computes one step of the double Gram-Schmidt orthogonalization process of the
          !! `abstract_vector` `y` against the `abstract_vector` basis `X`
          class(abstract_vector_rsp), intent(inout) :: y
          !! Input `abstract_vector` to orthogonalize
          class(abstract_vector_rsp), intent(in)    :: X(:)
          !! Input `abstract_vector` basis to orthogonalize against
          integer, intent(out) :: info
          !! Information flag.
          logical,                          optional, intent(in)    :: if_chk_orthonormal
          !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
          real(sp),                         optional, intent(out)   :: beta(:)
          !! Projection coefficients if requested
        end subroutine

        module subroutine DGS_basis_against_basis_rsp(y, X, info, if_chk_orthonormal, beta)
            !! Computes one step of the double Gram-Schmidt orthogonalization process of the
            !! `abstract_vector` `y` against the `abstract_vector` basis `X`
            class(abstract_vector_rsp), intent(inout) :: Y(:)
            !! Input `abstract_vector` basis to orthogonalize
            class(abstract_vector_rsp), intent(in)    :: X(:)
            !! Input `abstract_vector` basis to orthogonalize against
            integer, intent(out) :: info
            !! Information flag.
            logical,                          optional, intent(in)    :: if_chk_orthonormal
            !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
            real(sp),                         optional, intent(out)   :: beta(:,:)
            !! Projection coefficients if requested
        end subroutine
        module subroutine DGS_vector_against_basis_rdp(y, X, info, if_chk_orthonormal, beta)
          !! Computes one step of the double Gram-Schmidt orthogonalization process of the
          !! `abstract_vector` `y` against the `abstract_vector` basis `X`
          class(abstract_vector_rdp), intent(inout) :: y
          !! Input `abstract_vector` to orthogonalize
          class(abstract_vector_rdp), intent(in)    :: X(:)
          !! Input `abstract_vector` basis to orthogonalize against
          integer, intent(out) :: info
          !! Information flag.
          logical,                          optional, intent(in)    :: if_chk_orthonormal
          !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
          real(dp),                         optional, intent(out)   :: beta(:)
          !! Projection coefficients if requested
        end subroutine

        module subroutine DGS_basis_against_basis_rdp(y, X, info, if_chk_orthonormal, beta)
            !! Computes one step of the double Gram-Schmidt orthogonalization process of the
            !! `abstract_vector` `y` against the `abstract_vector` basis `X`
            class(abstract_vector_rdp), intent(inout) :: Y(:)
            !! Input `abstract_vector` basis to orthogonalize
            class(abstract_vector_rdp), intent(in)    :: X(:)
            !! Input `abstract_vector` basis to orthogonalize against
            integer, intent(out) :: info
            !! Information flag.
            logical,                          optional, intent(in)    :: if_chk_orthonormal
            !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
            real(dp),                         optional, intent(out)   :: beta(:,:)
            !! Projection coefficients if requested
        end subroutine
        module subroutine DGS_vector_against_basis_csp(y, X, info, if_chk_orthonormal, beta)
          !! Computes one step of the double Gram-Schmidt orthogonalization process of the
          !! `abstract_vector` `y` against the `abstract_vector` basis `X`
          class(abstract_vector_csp), intent(inout) :: y
          !! Input `abstract_vector` to orthogonalize
          class(abstract_vector_csp), intent(in)    :: X(:)
          !! Input `abstract_vector` basis to orthogonalize against
          integer, intent(out) :: info
          !! Information flag.
          logical,                          optional, intent(in)    :: if_chk_orthonormal
          !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
          complex(sp),                         optional, intent(out)   :: beta(:)
          !! Projection coefficients if requested
        end subroutine

        module subroutine DGS_basis_against_basis_csp(y, X, info, if_chk_orthonormal, beta)
            !! Computes one step of the double Gram-Schmidt orthogonalization process of the
            !! `abstract_vector` `y` against the `abstract_vector` basis `X`
            class(abstract_vector_csp), intent(inout) :: Y(:)
            !! Input `abstract_vector` basis to orthogonalize
            class(abstract_vector_csp), intent(in)    :: X(:)
            !! Input `abstract_vector` basis to orthogonalize against
            integer, intent(out) :: info
            !! Information flag.
            logical,                          optional, intent(in)    :: if_chk_orthonormal
            !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
            complex(sp),                         optional, intent(out)   :: beta(:,:)
            !! Projection coefficients if requested
        end subroutine
        module subroutine DGS_vector_against_basis_cdp(y, X, info, if_chk_orthonormal, beta)
          !! Computes one step of the double Gram-Schmidt orthogonalization process of the
          !! `abstract_vector` `y` against the `abstract_vector` basis `X`
          class(abstract_vector_cdp), intent(inout) :: y
          !! Input `abstract_vector` to orthogonalize
          class(abstract_vector_cdp), intent(in)    :: X(:)
          !! Input `abstract_vector` basis to orthogonalize against
          integer, intent(out) :: info
          !! Information flag.
          logical,                          optional, intent(in)    :: if_chk_orthonormal
          !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
          complex(dp),                         optional, intent(out)   :: beta(:)
          !! Projection coefficients if requested
        end subroutine

        module subroutine DGS_basis_against_basis_cdp(y, X, info, if_chk_orthonormal, beta)
            !! Computes one step of the double Gram-Schmidt orthogonalization process of the
            !! `abstract_vector` `y` against the `abstract_vector` basis `X`
            class(abstract_vector_cdp), intent(inout) :: Y(:)
            !! Input `abstract_vector` basis to orthogonalize
            class(abstract_vector_cdp), intent(in)    :: X(:)
            !! Input `abstract_vector` basis to orthogonalize against
            integer, intent(out) :: info
            !! Information flag.
            logical,                          optional, intent(in)    :: if_chk_orthonormal
            !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
            complex(dp),                         optional, intent(out)   :: beta(:,:)
            !! Projection coefficients if requested
        end subroutine
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
        !!  - `n`   :   Number of selected eigenvalues moved to the upper left-block of the 
        !!              Schur matrix. It is an `intent(out)` argument.
        !!
        !!  - `X`   :   On entry, array of `abstract_vector` computed using the Arnoldi process.
        !!              On exit, the first `n` columns form an orthonormal basis for the eigenspace
        !!              associated with eigenvalues moved to the upper left-block of the Schur matrix.
        !!              It is an `intent(inout)` argument.
        !!
        !!  - `H`   :   On entry, `real` of `complex` upper Hessenberg matrix computed using the
        !!              Arnoldi process. On exit, the leading \( n \times n\) block contains the
        !!              \( S_{11} \) block of the re-ordered Schur matrix containing the selected
        !!              eigenvalues. It is an `intent(inout)` argument.
        !!
        !!  - `select_eigs` :   Procedure to select which eigenvalues to move in the upper-left
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
            logical, allocatable          :: selected(:)
        end function eigvals_select_sp
        function eigvals_select_dp(lambda) result(selected)
            import dp
            complex(dp), intent(in) :: lambda(:)
            logical, allocatable          :: selected(:)
        end function eigvals_select_dp
    end interface

contains

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

        integer :: kdim
        
        ! Schur-related.
        real(sp) :: Z(size(H, 2), size(H, 2)), T(size(H, 2), size(H, 2))
        complex(sp) :: eigvals(size(H, 2))
        logical :: selected(size(H, 2))
       
        ! Krylov subspace dimension.
        kdim = size(X)-1

        ! Schur decomposition of the Hessenberg matrix.
        call schur(H(:size(H, 2), :), T, Z, eigvals) ; H(:size(H, 2), :) = T

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

        integer :: kdim
        
        ! Schur-related.
        real(dp) :: Z(size(H, 2), size(H, 2)), T(size(H, 2), size(H, 2))
        complex(dp) :: eigvals(size(H, 2))
        logical :: selected(size(H, 2))
       
        ! Krylov subspace dimension.
        kdim = size(X)-1

        ! Schur decomposition of the Hessenberg matrix.
        call schur(H(:size(H, 2), :), T, Z, eigvals) ; H(:size(H, 2), :) = T

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

        integer :: kdim
        
        ! Schur-related.
        complex(sp) :: Z(size(H, 2), size(H, 2)), T(size(H, 2), size(H, 2))
        complex(sp) :: eigvals(size(H, 2))
        logical :: selected(size(H, 2))
       
        ! Krylov subspace dimension.
        kdim = size(X)-1

        ! Schur decomposition of the Hessenberg matrix.
        call schur(H(:size(H, 2), :), T, Z, eigvals) ; H(:size(H, 2), :) = T

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

        integer :: kdim
        
        ! Schur-related.
        complex(dp) :: Z(size(H, 2), size(H, 2)), T(size(H, 2), size(H, 2))
        complex(dp) :: eigvals(size(H, 2))
        logical :: selected(size(H, 2))
       
        ! Krylov subspace dimension.
        kdim = size(X)-1

        ! Schur decomposition of the Hessenberg matrix.
        call schur(H(:size(H, 2), :), T, Z, eigvals) ; H(:size(H, 2), :) = T

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
