! This module is dedicated to iterative solvers for linear equations and eigenvalue problems.
! It imports several modules to handle vectors, linear operations, and Krylov decompositions.
! Additionally, it imports modules from the standard library.
module IterativeSolvers

    use AbstractVector
    use LinearOperator
    use KrylovDecomp

    use stdlib_sorting, only: sort_index, int_size
    use stdlib_optval, only: optval
    use stdlib_io_npy, only: save_npy
    use stdlib_error, only: error_stop
    use stdlib_logger

    implicit none
    include "dtypes.h"

    ! Create a logger instance for structured logging
    type(logger_type) :: logger

    private
    public :: eigs, eighs, gmres, save_eigenspectrum, svds, cg, bicgstab

contains

    !-----------------------------
    !-----     UTILITIES     -----
    !-----------------------------

    !=======================================================================================
    ! Compute Residual for Eigenpairs (compute_residual)
    !=======================================================================================
    !
    ! Purpose and formulation:
    ! --------
    ! Computes the residual associated with an eigenpair, useful for assessing the quality
    ! of the approximated eigenvalues and eigenvectors.
    !
    ! Given the norm of the Krylov residual vector (beta) and the last element of the Ritz
    ! eigenvector (x), the residual is calculated as follows:
    !
    ! residual = |beta * x|
    !
    ! Algorithmic Features:
    ! ----------------------
    ! - Utilizes the absolute value function for residual computation.
    ! - Quick and efficient way to gauge the accuracy of the approximated eigenpair.
    ! - Only computes the residual for a single eigenpair at a time.
    !
    ! Input/Output Parameters:
    ! ------------------------
    ! - beta      : Norm of Krylov residual vector   [Input]
    ! - x         : Last element of Ritz eigenvector [Input]
    ! - residual  : Computed residual                [Output]
    !
    !=======================================================================================
    elemental pure function compute_residual(beta, x) result(residual)
        !> Norm of Krylov residual vector.
        real(kind=wp), intent(in) :: beta
        !> Last element of Ritz eigenvector.
        real(kind=wp), intent(in) :: x
        !> Residual.
        real(kind=wp) :: residual

        ! --> Compute residual.
        residual = abs(beta*x)
        return
    end function compute_residual

    !=======================================================================================
    ! Save Eigenvalues and Residuals to Disk (save_eigenspectrum)
    !=======================================================================================
    !
    ! Purpose and formulation:
    ! -----------------------
    ! Saves the eigenvalues and corresponding residuals to disk for further analysis or
    ! post-processing. The eigenvalues are split into their real and imaginary parts.
    !
    ! Given arrays of real and imaginary parts of eigenvalues (real_part, imag_part)
    ! and residuals, the subroutine saves this data into a structured array:
    !
    ! data = [real_part, imag_part, residuals]
    !
    ! Algorithmic Features:
    ! ----------------------
    ! - Combines real and imaginary parts along with residuals into a single 2D array.
    ! - Utilizes the npy file format for saving the array to disk.
    ! - Facilitates easy storage and retrieval of the eigenvalues and residuals.
    ! - Utilizes a commonly-used file format (npy) for potential compatibility with Python tools.
    ! - Assumes that the real and imaginary parts of the eigenvalues and the residuals are of the same length.
    !
    ! Input/Output Parameters:
    ! ------------------------
    ! - real_part : Real part of the eigenvalues       [Input]
    ! - imag_part : Imaginary part of the eigenvalues  [Input]
    ! - residuals : Residual norms                     [Input]
    ! - filename  : Name of the output file            [Input]
    !
    !=======================================================================================
    subroutine save_eigenspectrum(real_part, imag_part, residuals, filename)
        !> Real and imaginary parts of the eigenvalues.
        real(kind=wp), intent(in) :: real_part(:)
        real(kind=wp), intent(in) :: imag_part(:)
        !> Residual norm computed from the Arnoldi/Lanczos factorization.
        real(kind=wp), intent(in) :: residuals(:)
        !> Name of the output file.
        character(len=*), intent(in) :: filename

        !> Miscellaneous.
        real(kind=wp), dimension(size(real_part), 3) :: data

        ! --> Store the data.
        data(:, 1) = real_part; data(:, 2) = imag_part; data(:, 3) = residuals
        ! --> Save the eigenspectrum to disk using npy file format.
        call save_npy(filename, data)

        return
    end subroutine save_eigenspectrum

    !------------------------------------------
    !-----                                -----
    !-----     EIGENVALUE COMPUTATION     -----
    !-----                                -----
    !------------------------------------------

    !=======================================================================================
    ! Eigenvalue and Eigenvector Solver Subroutine (EIGS)
    !=======================================================================================
    !
    ! Purpose:
    ! --------
    ! Solves the eigenvalue problem for a given square linear operator A using Arnoldi factorization.
    ! The eigenvalues (eigvals) and eigenvectors (eigvecs) are computed within a Krylov subspace (X).
    !
    ! Mathematical Formulation:
    ! -------------------------
    ! Given a square linear operator A of size (n x n), find eigvals and eigvecs such that:
    !
    ! A * eigvecs = eigvals * eigvecs
    !
    ! The Krylov subspace X is formed via Arnoldi factorization, resulting in an upper Hessenberg matrix H.
    ! The eigenvalues of A are approximated by the eigenvalues of H, and the eigenvectors are computed accordingly.
    !
    ! Algorithmic Features:
    ! ----------------------
    ! - Builds a Krylov subspace (X) using Arnoldi factorization.
    ! - Computes eigenpairs of the reduced upper Hessenberg matrix (H).
    ! - Sorts eigvals and associates eigvecs based on magnitude.
    !
    ! Advantages:
    ! -----------
    ! - Suitable for a wide range of square matrices A, including non-symmetric ones.
    ! - Efficient for large-scale problems.
    ! - Allows user-specified initial Krylov vector (X(1)).
    !
    ! Limitations:
    ! ------------
    ! - Accuracy is dependent on the quality of the Krylov subspace (X).
    ! - No preconditioning capabilities in the current implementation.
    !
    ! Input/Output Parameters:
    ! ------------------------
    ! - A          : Linear Operator                [Input]
    ! - X          : Initial Krylov vectors         [Input/Output]
    ! - eigvecs    : Eigenvectors in Krylov basis   [Output]
    ! - eigvals    : Eigenvalues                    [Output]
    ! - residuals  : Residuals for each eigenvalue  [Output]
    ! - info       : Iteration Information flag     [Output]
    ! - verbosity  : Verbosity control flag         [Optional, Input]
    ! - transpose  : Use transpose of A             [Optional, Input]
    !
    ! References:
    ! -----------
    ! - Arnoldi, W. E. (1951). "The Principle of Minimized Iterations in the Solution of the Matrix Eigenvalue Problem,"
    !   Quarterly of Applied Mathematics, 9(1), 17–29.
    !
    !=======================================================================================
    subroutine eigs(A, X, eigvecs, eigvals, residuals, info, verbosity, transpose)
        !> Linear Operator.
        class(abstract_linop), intent(in) :: A
        !> Krylov basis.
        class(abstract_vector), dimension(:), intent(inout) :: X
        !> Coordinates of eigenvectors in Krylov basis, eigenvalues and associated residuals.
        complex(kind=wp), dimension(size(X) - 1, size(X) - 1), intent(out) :: eigvecs
        complex(kind=wp), dimension(size(X) - 1), intent(out) :: eigvals
        real(kind=wp), dimension(size(X) - 1), intent(out) :: residuals
        !> Information flag.
        integer, intent(out) :: info
        !> Verbosity control.
        logical, optional, intent(in) :: verbosity
        logical :: verbose
        !> Transpose operator.
        logical, optional, intent(in) :: transpose
        logical :: trans

        !> Upper Hessenberg matrix.
        real(kind=wp), dimension(size(X), size(X) - 1) :: H
        !> Krylov subspace dimension.
        integer :: kdim
        !> Miscellaneous.
        integer :: i
        integer(int_size), dimension(size(X) - 1) :: indices
        real(kind=wp), dimension(size(X) - 1) :: abs_eigvals

        ! --> Initialize the Krylov subspace dimension.
        kdim = size(X) - 1

        ! --> Handle optional arguments.
        verbose = optval(verbosity, .false.)
        trans = optval(transpose, .false.)

        ! --> Initialize variables and Krylov vectors.
        H = 0.0_wp; residuals = 0.0_wp; eigvals = (0.0_wp, 0.0_wp); eigvecs = (0.0_wp, 0.0_wp)
        call initialize_krylov_basis(X)

        ! --> Compute Arnoldi factorization to build the Krylov subspace.
        call arnoldi_factorization(A, X, H, info, verbosity=verbose, transpose=trans)

        if (info < 0) then
            if (verbose) then
                write (*, *) "INFO : Arnoldi iteration failed. Exiting the eig subroutine."
                write (*, *) "       Arnoldi exit code :", info
            end if
            info = -1
            return
        end if

        ! --> Compute spectral decomposition of the Hessenberg matrix.
        call evd(H(1:kdim, 1:kdim), eigvecs, eigvals, kdim)

        ! --> Sort eigenvalues with decreasing magnitude.
        abs_eigvals = abs(eigvals); call sort_index(abs_eigvals, indices, reverse=.true.)
        eigvals(:) = eigvals(indices); eigvecs = eigvecs(:, indices)

        ! --> Compute the residual associated with each eigenpair.
        residuals = compute_residual(H(kdim + 1, kdim), abs(eigvecs(kdim, :)))

        return
    end subroutine eigs

    !=======================================================================================
    ! Eigenvalue and Eigenvector Solver for Symmetric Positive Definite (SPD) Matrices (EIGHS)
    !=======================================================================================
    !
    ! Purpose:
    ! --------
    ! Solves the eigenvalue problem for SPD matrices (A) using the Lanczos algorithm.
    ! The mathematical formulation is: A * eigvecs = eigvals * eigvecs, where A is of size (n x n).
    !
    ! Algorithmic Features:
    ! ----------------------
    ! - Employs Lanczos tridiagonalization to construct a Krylov subspace represented by matrix T of size (kdim x kdim).
    ! - Computes eigenpairs (eigvals, eigvecs) of T to approximate those of A.
    ! - Sorts eigvals and eigvecs based on magnitude.
    ! - Computes residuals to assess the quality of approximations.
    !
    ! Advantages:
    ! -----------
    ! - Efficient for large, sparse, SPD matrices (A).
    ! - Lower memory requirements compared to full spectrum methods like QR.
    ! - Faster convergence rates for extreme eigvals.
    !
    ! Limitations:
    ! ------------
    ! - Only applicable to SPD matrices (A).
    ! - Eigvals and eigvecs are approximations and may require refinement.
    ! - Not suitable for finding all eigvals in large matrices.
    !
    ! Input/Output Parameters:
    ! ------------------------
    ! - A         : Linear Operator (assumed SPD)  [Input]
    ! - X         : Krylov basis vectors           [Input/Output]
    ! - eigvecs   : Eigenvectors in Krylov basis   [Output]
    ! - eigvals   : Eigenvalues                    [Output]
    ! - residuals : Residuals of eigenpairs        [Output]
    ! - info      : Iteration Information flag     [Output]
    ! - verbosity : Verbosity control flag         [Optional, Input]
    !
    ! References:
    ! -----------
    ! - Lanczos, C. (1950). "An Iteration Method for the Solution of the Eigenvalue Problem of Linear Differential and Integral Operators".
    !   United States Governm. Press Office.
    !
    !=======================================================================================
    subroutine eighs(A, X, eigvecs, eigvals, residuals, info, verbosity)
        !> Linear Operator (assumed SPD)
        class(abstract_spd_linop), intent(in) :: A
        !> Krylov basis.
        class(abstract_vector), dimension(:), intent(inout) :: X
        !> Coordinates of eigenvectors in Krylov basis, eigenvalues, and associated residuals
        real(kind=wp), dimension(size(X) - 1, size(X) - 1), intent(out) :: eigvecs
        real(kind=wp), dimension(size(X) - 1), intent(out) :: eigvals
        real(kind=wp), dimension(size(X) - 1), intent(out) :: residuals
        !> Information flag.
        integer, intent(out) :: info
        !> Verbosity control.
        logical, optional, intent(in) :: verbosity
        logical :: verbose

        !> Tridiagonal matrix.
        real(kind=wp), dimension(size(X), size(X) - 1) :: T
        !> Krylov subspace dimension.
        integer :: kdim
        !> Miscellaneous.
        integer :: i
        integer(int_size), dimension(size(X) - 1) :: indices

        ! --> Dimension of the Krylov subspace.
        kdim = size(X) - 1

        ! --> Deals with the optional argument.
        verbose = optval(verbosity, .false.)

        ! --> Initialize all variables.
        T = 0.0_wp; residuals = 0.0_wp; eigvecs = 0.0_wp; eigvals = 0.0_wp
        call initialize_krylov_basis(X)

        ! --> Compute Lanczos tridiagonalization.
        call lanczos_tridiagonalization(A, X, T, info, verbosity=verbose)

        if (info < 0) then
            if (verbose) then
                write (*, *) "INFO : Lanczos iteration failed. Exiting the eigh subroutine."
                write (*, *) "       Lanczos exit code :", info
            end if
            info = -1
            return
        end if

        ! --> Compute spectral decomposition of the tridiagonal matrix.
        call hevd(T(1:kdim, 1:kdim), eigvecs, eigvals, kdim)

        ! --> Sort eigenvalues in decreasing order.
        call sort_index(eigvals, indices, reverse=.true.)
        eigvecs = eigvecs(:, indices)

        ! --> Compute the residual associated with each eigenpair.
        residuals = compute_residual(T(kdim + 1, kdim), eigvecs(kdim, :))

        return
    end subroutine eighs

    !=======================================================================================
    ! Eigendecomposition Subroutine (evd)
    !=======================================================================================
    !
    ! Purpose:
    ! --------
    ! Performs the eigendecomposition of a given square matrix 'A' and returns the eigenvalues
    ! ('vals') and corresponding eigenvectors ('vecs'). Utilizes LAPACK's dgeev routine for the
    ! numerical operations.
    !
    ! Mathematical Formulation:
    ! -------------------------
    ! Given a square matrix A of dimension n x n, find eigenvalues (lambda) and corresponding
    ! eigenvectors (v) such that:
    !
    ! A * v = lambda * v
    !
    ! Algorithmic Features:
    ! ----------------------
    ! - Calls LAPACK's dgeev function for efficient and robust eigendecomposition.
    ! - Converts real and imaginary components of eigenvalues and eigenvectors to complex(kind=wp) types.
    ! - Extracts the right eigenvectors ('vr') and combines them with corresponding imaginary parts, if any.
    !
    ! Advantages:
    ! -----------
    ! - Relies on LAPACK's established numerical algorithms, ensuring robustness and stability.
    ! - Capable of handling both real and complex eigenvalues and eigenvectors.
    ! - Suitable for dense matrices of moderate dimensions.
    !
    ! Limitations:
    ! ------------
    ! - Operates on dense matrices, making it memory-intensive for large-scale problems.
    ! - Not designed for sparse or structured matrices.
    ! - Requires additional effort for optimal workspace memory ('lwork') setting, defaulted here to 4*n.
    !
    ! Input/Output Parameters:
    ! ------------------------
    ! - A     : Input Matrix (Dimension: n x n)             [Input]
    ! - vecs  : Eigenvectors (Dimension: n x n)             [Output]
    ! - vals  : Eigenvalues (Dimension: n)                  [Output]
    ! - n     : Dimension of the square matrix A            [Input]
    !
    ! References:
    ! -----------
    ! - LAPACK Users' Guide, Third Edition, SIAM. Specifically, the dgeev routine.
    !
    !=======================================================================================
    subroutine evd(A, vecs, vals, n)
        !> Lapack job.
        character*1 :: jobvl = "N", jobvr = "V"
        integer :: n, lwork, info, lda, ldvl, ldvr
        real(kind=wp), dimension(n, n) :: A, A_tilde, vr
        real(kind=wp), dimension(1, n) :: vl
        real(kind=wp), dimension(4*n)  :: work
        real(kind=wp), dimension(n)    :: wr, wi
        complex(kind=wp), dimension(n, n)   :: vecs
        complex(kind=wp), dimension(n)      :: vals
        integer :: i

        interface
            pure subroutine dgeev(fjobvl, fjobvr, fn, fa, flda, fwr, fwi, fvl, fldvl, fvr, fldvr, fwork, flwork, finfo)
                import wp
                character, intent(in) :: fjobvl, fjobvr
                integer, intent(in) :: fn, flda, fldvl, fldvr, flwork, finfo
                real(kind=wp), intent(inout) :: fa(flda, *)
                real(kind=wp), intent(out) :: fwr(fn), fwi(fn), fvl(fldvl, *), fvr(fldvr, *), fwork(flwork)
            end subroutine dgeev
        end interface

        ! --> Compute the eigendecomposition of A.
        lda = n; ldvl = 1; ldvr = n; lwork = 4*n; A_tilde = A
        call dgeev(jobvl, jobvr, n, A_tilde, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)

        ! --> Real to complex arithmetic.
        !     NOTE : Check if a LAPACK function already exists for that purpose.
        vals = wr*(1.0_wp, 0.0_wp) + wi*(0.0_wp, 1.0_wp)
        vecs = vr*(1.0_wp, 0.0_wp)

        do i = 1, n - 1
            if (wi(i) .gt. 0) then
                vecs(:, i) = vr(:, i)*(1.0_wp, 0.0_wp) + vr(:, i + 1)*(0.0_wp, 1.0_wp)
                vecs(:, i + 1) = conjg(vecs(:, i))
            else if (abs(wi(i)) .le. epsilon(wi(i))) then
                vecs(:, i) = vr(:, i)*(1.0_wp, 0.0_wp)
            end if
        end do

        return
    end subroutine evd

    !=======================================================================================
    ! Hermitian Eigendecomposition Subroutine (hevd)
    !=======================================================================================
    !
    ! Purpose:
    ! --------
    ! Performs the eigendecomposition of a given Hermitian (or real symmetric) matrix 'A' and
    ! returns the eigenvalues ('vals') and corresponding eigenvectors ('vecs'). Utilizes LAPACK's
    ! dsyev routine optimized for Hermitian matrices.
    !
    ! Mathematical Formulation:
    ! -------------------------
    ! Given a Hermitian matrix A of dimension n x n, find eigenvalues (lambda) and corresponding
    ! eigenvectors (v) such that:
    !
    ! A * v = lambda * v
    !
    ! Algorithmic Features:
    ! ----------------------
    ! - Calls LAPACK's dsyev function, a specialized subroutine for Hermitian matrices.
    ! - Employs a divide-and-conquer strategy for efficient eigenvalue and eigenvector computation.
    ! - Adapts to both real symmetric and complex Hermitian matrices.
    !
    ! Advantages:
    ! -----------
    ! - Tailored for Hermitian matrices, providing high numerical stability.
    ! - Efficient for dense Hermitian matrices of moderate dimensions.
    ! - Direct linkage with LAPACK ensures reliability and performance.
    !
    ! Limitations:
    ! ------------
    ! - May consume extensive memory for large Hermitian matrices.
    ! - Inefficient for sparse or structured Hermitian matrices (e.g., banded, tridiagonal).
    !
    ! Input/Output Parameters:
    ! ------------------------
    ! - A    : Input Hermitian Matrix (Dimension: n x n)         [Input]
    ! - vecs : Eigenvectors (Dimension: n x n)                   [Output]
    ! - vals : Eigenvalues (Dimension: n)                        [Output]
    ! - n    : Dimension of the square matrix A                   [Input]
    !
    ! References:
    ! -----------
    ! - LAPACK Users' Guide, Third Edition, SIAM.
    ! - Golub, G. H., & Van Loan, C. F. (2012). "Matrix Computations," Fourth Edition.
    !
    !=======================================================================================
    subroutine hevd(A, vecs, vals, n)
        !> Lapack job.
        character :: jobz = "V", uplo = "U"
        integer :: n, lwork, info, lda
        real(kind=wp), dimension(n) :: vals
        real(kind=wp), dimension(n, n) :: A, A_tilde, vecs ! vecs is unused
        real(kind=wp), dimension(3*n - 1) :: work

        interface
            pure subroutine dsyev(fjobz, fuplo, fn, fa, flda, fw, fwork, flwork, finfo)
                import wp
                character, intent(in) :: fjobz, fuplo
                integer, intent(in)  :: fn
                integer, intent(in)  :: flda
                integer, intent(in)  :: flwork
                integer, intent(out) :: finfo
                real(kind=wp), intent(inout) :: fa(flda, *)
                real(kind=wp), intent(out)   :: fw(*)
                real(kind=wp), intent(out)   :: fwork(*)
            end subroutine dsyev
        end interface

        ! --> Compute the eigendecomposition of A.
        lda = n; lwork = 3*n - 1; A_tilde = A
        call dsyev(jobz, uplo, n, A_tilde, lda, vals, work, lwork, info)

        return
    end subroutine hevd

    !----------------------------------------------
    !-----                                    -----
    !-----     SINGULAR VALUE COMPUTATION     -----
    !-----                                    -----
    !----------------------------------------------

    !=======================================================================================
    ! Singular Value Decomposition Subroutine (svds) - lanczos based / matrix-free
    !=======================================================================================
    !
    ! Purpose:
    ! --------
    ! Computes the singular value decomposition of a given linear operator A using Lanczos
    ! bidiagonalization. The subroutine finds matrices U, Sigma, and V such that A = U Sigma V^T.
    !
    ! Mathematical Formulation:
    ! -------------------------
    ! Given a linear operator A, find matrices U, Sigma, and V such that:
    ! A = U Sigma V^T
    !
    ! Algorithmic Features:
    ! ----------------------
    ! - Utilizes Lanczos bidiagonalization to approximate A with a bidiagonal matrix B.
    ! - Computes the singular value decomposition of B to obtain singular values and vectors.
    ! - Calculates residuals for singular triplets as: residual = | A v - sigma u |.
    ! - Optionally verbose, providing runtime diagnostics.
    !
    ! Advantages:
    ! -----------
    ! - Efficient for large, sparse matrices.
    ! - Suitable for problems requiring only a subset of singular values and vectors.
    ! - Can be adapted for preconditioning techniques.
    !
    ! Limitations:
    ! ------------
    ! - Not optimized for dense matrices or those with special structures (e.g., Toeplitz, circulant).
    ! - May require many iterations for ill-conditioned matrices.
    !
    ! Input/Output Parameters:
    ! ------------------------
    ! - A          : Linear Operator                               [Input]
    ! - U, V       : Krylov bases for left (U) and right (V) singular vectors [Input/Output]
    ! - uvecs, vvecs: Coordinates of left (U) and right (V) singular vectors in Krylov bases [Output]
    ! - sigma      : Singular values                               [Output]
    ! - residuals  : Residuals associated with singular triplets  [Output]
    ! - info       : Information flag                              [Output]
    ! - verbosity  : Control for runtime diagnostics               [Optional, Input]
    !
    ! References:
    ! -----------
    ! - Golub, G. H., & Kahan, W. (1965). "Calculating the Singular Values and Pseudo-Inverse of a Matrix."
    ! - Baglama, J., & Reichel, L. (2005). "Augmented implicitly restarted Lanczos bidiagonalization methods."
    !
    !=======================================================================================
    subroutine svds(A, U, V, uvecs, vvecs, sigma, residuals, info, verbosity)
        !> Linear Operator.
        class(abstract_linop), intent(in) :: A
        !> Krylov bases.
        class(abstract_vector), intent(inout) :: U(:) ! Basis for left sing. vectors.
        class(abstract_vector), intent(inout) :: V(:) ! Basis for right sing. vectors.
        !> Coordinates of singular vectors in Krylov bases, singular values, and associated residuals.
        real(kind=wp), intent(out) :: uvecs(size(U) - 1, size(U) - 1), vvecs(size(U) - 1, size(U) - 1)
        real(kind=wp), intent(out) :: sigma(size(U) - 1)
        real(kind=wp), intent(out) :: residuals(size(U) - 1)
        !> Information flag.
        integer, intent(out) :: info
        !> Verbosity control.
        logical, optional, intent(in) :: verbosity
        logical verbose

        !> Bidiagonal matrix.
        real(kind=wp) :: B(size(U), size(U) - 1)
        !> Krylov subspace dimension.
        integer :: kdim
        !> Miscellaneous.
        !integer(int_size) :: indices(size(U) - 1)
        integer :: i

        ! --> Deals with the optional args.
        verbose = optval(verbosity, .false.)
        ! --> Assert size(U) == size(V).
        if (size(U) .ne. size(V)) then
            info = -1
            if (verbose) then
                write (*, *) "INFO : Left and Right Krylov subspaces have different dimensions."
                write (*, *) "       Exiting svds with exit code info =", info
            end if
        else
            kdim = size(U) - 1
        end if

        ! --> Initialize variables.
        B = 0.0_wp; residuals = 0.0_wp; uvecs = 0.0_wp; vvecs = 0.0_wp; sigma = 0.0_wp
        do i = 2, size(U)
            call U(i)%zero(); call V(i)%zero()
        end do

        ! --> Compute the Lanczos bidiagonalization.
        call lanczos_bidiagonalization(A, U, V, B, info, verbosity=verbose)

        ! --> Compute the singular value decomposition of the bidiagonal matrix.
        call svd(B(1:kdim, 1:kdim), uvecs, sigma, vvecs)

        ! --> Compute the residual associated with each singular triplet.
        residuals = compute_residual(B(kdim + 1, kdim), vvecs(kdim, :))

        return
    end subroutine svds

    !=======================================================================================
    ! Singular Value Decomposition (SVD) Subroutine (svd) - dense matrix (full matrix)
    !=======================================================================================
    !
    ! Purpose:
    ! --------
    ! Computes the singular value decomposition of a given matrix A. The subroutine finds
    ! matrices U, S, and V such that A = U S V^T, where U is an m x m orthogonal matrix,
    ! S is an m x n diagonal matrix with non-negative real numbers (singular values),
    ! and V is an n x n orthogonal matrix.
    !
    ! Mathematical Formulation:
    ! -------------------------
    ! Given a matrix A, find matrices U, S, and V such that:
    ! A = U S V^T
    !
    ! Algorithmic Features:
    ! ----------------------
    ! - Utilizes LAPACK's DGESVD routine for the decomposition.
    ! - Employs the "thin" SVD option, storing only the first min(m, n) columns of U and V.
    ! - Operates directly on the input matrix A, which gets overwritten.
    ! - Dynamically allocates workspace based on matrix dimensions.
    !
    ! Advantages:
    ! -----------
    ! - Robust and widely used numerical method with applications in data compression, image processing, etc.
    ! - "Thin" SVD reduces memory requirements.
    !
    ! Limitations:
    ! ------------
    ! - Computationally expensive for large matrices, with time complexity scaling as O(m n^2) or O(n m^2).
    ! - Requires the entire matrix A to be stored in memory.
    !
    ! Input/Output Parameters:
    ! ------------------------
    ! - A     : Matrix to be factorized                         [Input]
    ! - U     : Left singular vectors U                         [Output]
    ! - S     : Singular values S                               [Output]
    ! - V     : Right singular vectors V                        [Output]
    !
    ! References:
    ! -----------
    ! - Golub, G. H., & Reinsch, C. (1970). "Singular value decomposition and least squares solutions,"
    !   Numerische Mathematik, 14(5), 403–420.
    ! - LAPACK Users' Guide, Third Edition.
    !
    !=======================================================================================
    subroutine svd(A, U, S, V)
        !> Matrix to be factorized.
        real(kind=wp), intent(in)  :: A(:, :)
        !> Left singular vectors.
        real(kind=wp), intent(out) :: U(size(A, 1), min(size(A, 1), size(A, 2)))
        !> Singular values.
        real(kind=wp), intent(out) :: S(size(A, 2))
        !> Right singular vectors.
        real(kind=wp), intent(out) :: V(size(A, 2), min(size(A, 1), size(A, 2)))

        !> Lapack job.
        character                  :: jobu = "S", jobvt = "S"
        integer                    :: m, n, lda, ldu, ldvt, lwork, info
        real(kind=wp), allocatable :: work(:)
        real(kind=wp) :: A_tilde(size(A, 1), size(A, 2)), vt(min(size(A, 1), size(A, 2)), size(A, 2))

        interface
            pure subroutine dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
                import wp
                character, intent(in) :: jobu, jobvt
                integer, intent(in) :: m, n, lda, ldu, ldvt, lwork, info
                real(kind=wp), intent(inout) :: a(lda, *)
                real(kind=wp), intent(out)   :: u(ldu, *), s(*), vt(ldvt, *), work(*)
            end subroutine dgesvd
        end interface

        m = size(A, 1); n = size(A, 2)
        lda = size(A, 1); ldu = size(A, 1); ldvt = size(A, 2)
        lwork = max(1, 3*min(m, n), 5*min(m, n)); allocate (work(lwork))

        a_tilde = a
        call dgesvd(jobu, jobvt, m, n, a_tilde, lda, s, u, ldu, vt, ldvt, work, lwork, info)
        v = transpose(vt)

        return
    end subroutine svd

    !--------------------------------------------
    !-----                                  -----
    !-----     ITERATIVE LINEAR SOLVERS     -----
    !-----                                  -----
    !--------------------------------------------

    !=======================================================================================
    ! Generalized Minimal Residual (GMRES) Solver Subroutine
    !=======================================================================================
    !
    ! Purpose:
    ! --------
    ! Implements the classic, unpreconditioned Generalized Minimal Residual (GMRES) algorithm
    ! for solving nonsymmetric, non-Hermitian linear systems of equations.
    !
    ! Mathematical Formulation:
    ! -------------------------
    ! Solves the linear system of equations of the form Ax = b.
    !
    !  A x = b
    !
    ! where,
    !  A ∈ R^{m x n} : Linear Operator          [Input]
    !  b ∈ R^m       : Right-hand side vector   [Input]
    !  x ∈ R^n       : Initial/Updated solution [Input/Output]
    !
    ! Algorithmic Features:
    ! ----------------------
    ! - Constructs a full Krylov subspace without restarts (i.e., not GMRES(m)).
    ! - Utilizes Arnoldi factorization to generate an orthonormal basis for the Krylov subspace.
    ! - Employs a least-squares solve to determine the optimal linear combination of the Krylov vectors.
    ! - Updates the approximate solution based on the least-squares solution.
    !
    ! Advantages:
    ! -----------
    ! - Suitable for nonsymmetric and ill-conditioned matrices.
    ! - Produces monotonically decreasing residuals.
    ! - Fully utilizes the generated Krylov subspace for the solution.
    !
    ! Limitations:
    ! ------------
    ! - Memory-intensive due to the absence of restarts.
    ! - May not be efficient for very large-scale problems.
    ! - No preconditioning capabilities in the current implementation.
    !
    ! Logging:
    ! --------
    ! - Utilizes stdlib_logger for structured logging and debugging.
    ! - Supports different levels of verbosity for finer control over logging.
    !
    ! Input/Output Parameters:
    ! ------------------------
    ! - A        : Linear Operator              [Input]
    ! - b        : Right-hand side vector       [Input]
    ! - x        : Initial/Updated solution     [Input/Output]
    ! - info     : Iteration Information flag   [Output]
    ! - maxiter  : Maximum number of iterations [Optional, Input]
    ! - tol      : Tolerance for convergence    [Optional, Input]
    ! - verbosity: Verbosity control flag       [Optional, Input]
    !
    ! References:
    ! -----------
    ! - Saad, Y., and Schultz, M. H. (1986). "GMRES: A Generalized Minimal Residual Algorithm for Solving Nonsymmetric Linear Systems,"
    !   SIAM Journal on Scientific and Statistical Computing, 7(3), 856–869.
    !
    !=======================================================================================
    subroutine gmres(A, b, x, info, kdim, maxiter, tol, verbosity, transpose)

        !> Input: Linear Operator (implements matvec and rmatvec methods)
        class(abstract_linop), intent(in) :: A

        !> Input: Right-hand side vector (implements norm, zero, etc.)
        class(abstract_vector), intent(in) :: b

        !> Input/Output: Initial guess and updated solution
        class(abstract_vector), intent(inout) :: x

        !> Output: Information flag for the solver's status
        integer, intent(out) :: info

        !> Optional Input: Dimension of Krylov subspace
        integer, optional, intent(in) :: kdim
        integer :: k_dim  ! Actual value used in the algorithm

        !> Optional Input: Maximum number of full GMRES iterations
        integer, optional, intent(in) :: maxiter
        integer :: niter  ! Actual value used in the algorithm

        !> Optional Input: Tolerance for convergence
        real(kind=wp), optional, intent(in) :: tol
        real(kind=wp) :: tolerance  ! Actual value used in the algorithm

        !> Optional Input: Verbosity control flag
        logical, optional, intent(in) :: verbosity
        logical :: verbose  ! Actual value used in the algorithm

        !> Optional Input: Flag to indicate transpose operations
        logical, optional, intent(in) :: transpose
        logical :: trans  ! Actual value used in the algorithm

        !> Declare and allocate internal variables

        ! Krylov subspace vectors
        class(abstract_vector), allocatable :: V(:)

        ! Upper Hessenberg matrix for Arnoldi factorization
        real(kind=wp), allocatable :: H(:, :)

        ! Variables for least-squares solution
        real(kind=wp), allocatable :: y(:)
        real(kind=wp), allocatable :: e(:)
        real(kind=wp) :: beta

        ! Miscellaneous variables
        integer :: i, j, k
        class(abstract_vector), allocatable :: dx
        character(len=128) :: message

        ! Handle optional arguments with stdlib's optval
        k_dim = optval(kdim, 30)
        niter = optval(maxiter, 10)
        tolerance = optval(tol, atol + rtol*b%norm())
        verbose = optval(verbosity, .false.)
        trans = optval(transpose, .false.)

        !Initialize logger with appropriate verbosity level
        call logger%add_log_file('gmres.log')

      !   debug_level
      !   information_level
      !   warning_level
      !   error_level
      !   io_error_level
      !   text_error_level
      !   all_level
      !   none_level

        ! call logger%configure(level=information_level) ! default 
        if (verbose) then
            call logger%configure(level=debug_level)
        end if

        ! --> Initialize Krylov subspace.
        allocate (V(1:k_dim + 1), source=b)
        do i = 1, size(V)
            call V(i)%zero()
        end do
        allocate (H(k_dim + 1, k_dim)); H = 0.0_wp
        allocate (y(1:k_dim)); y = 0.0_wp
        allocate (e(1:k_dim + 1)); e = 0.0_wp

        ! --> Initial Krylov vector.
        if (trans) then
            call A%rmatvec(x, V(1))
        else
            call A%matvec(x, V(1))
        end if
        call V(1)%sub(b); call V(1)%scal(-1.0_wp)
        beta = V(1)%norm(); call V(1)%scal(1.0_wp/beta)

        gmres_iterations: do i = 1, niter
            ! --> Zero-out variables.
            H = 0.0_wp; y = 0.0_wp; e = 0.0_wp; e(1) = beta
            do j = 2, size(V)
                call V(j)%zero()
            end do

            arnoldi: do k = 1, k_dim
                ! --> Step-by-step Arnoldi factorization.
                call arnoldi_factorization(A, V, H, info, kstart=k, kend=k, verbosity=.false., tol=tolerance, transpose=trans)

                ! Error and logging
                if (info < 0) then ! Arnoldi only output info = -1 for error * discuss *
                    call logger%log_error("GMRES - Arnoldi Factorization failed with exit code")
                    call error_stop("GMRES - Arnoldi Factorization failed")
                end if

                ! --> Least-squares problem.
                call lstsq(H(1:k + 1, 1:k), e(1:k + 1), y(1:k))
                ! --> Compute residual.
                beta = norm2(e(1:k + 1) - matmul(H(1:k + 1, 1:k), y(1:k)))
               
                ! --> Check convergence.
                if (beta**2 .lt. tolerance) then
                  exit arnoldi
               endif

            enddo arnoldi

            ! --> Update solution.
            k = min(k, k_dim)
            if (allocated(dx) .eqv. .false.) allocate (dx, source=x)
            call dx%zero(); call get_vec(dx, V(1:k), y(1:k)); call x%add(dx)

            ! --> Recompute residual for sanity check.
            if (trans) then
                call A%rmatvec(x, V(1))
            else
                call A%matvec(x, V(1))
            end if

            call V(1)%sub(b); call V(1)%scal(-1.0_wp)

            ! --> Initialize new starting Krylov vector if needed.
            beta = V(1)%norm(); call V(1)%scal(1.0D+00/beta)

            ! --> Check convergence, informational and debug logs
            if (beta**2 .lt. tolerance) then
                  call logger%log_information("GMRES converged successfully.")
                  info = 0 ! Algorithm converged
                  exit gmres_iterations
            else
                  write (message, '(A, I0, A, E32.5)') "GMRES residual after ",(i-1)*k_dim+k, " iteration: ",beta**2
                  call logger%log_debug(trim(message))
            end if

        end do gmres_iterations

        ! Check if maximum iterations reached without convergence
        if (i >= niter) then
         info = niter ! following lapack info should take the value of the iteration (>0 for error)
         call error_stop("GMRES did not converge within maxiter.")
        endif

        ! --> Deallocate variables.
        deallocate (V, H, y, e)

        return
    end subroutine gmres

    !==============================================================================
    ! Least Squares Solver Subroutine
    !==============================================================================
    !
    ! Mathematical Formulation:
    ! -------------------------
    ! Solves the least squares problem of the form Ax = b, minimizing the residual ||Ax - b||.
    !
    !  minimize ||A x - b||
    !
    ! where,
    !  A ∈ R^{m x n} : Input matrix                 [Input]
    !  b ∈ R^m       : Right-hand side vector       [Input]
    !  x ∈ R^n       : Solution vector              [Output]
    !
    ! Algorithmic Features:
    ! ----------------------
    ! - Utilizes LAPACK's DGELS routine for solving the least squares problem.
    ! - Leverages QR or LQ factorization methods internally, based on the dimensions of A.
    !
    ! Advantages:
    ! -----------
    ! - Efficient and robust, leveraging optimized LAPACK routines.
    ! - Handles both overdetermined and underdetermined systems.
    ! - Can be extended to handle rank-deficient matrices.
    !
    ! Limitations:
    ! ------------
    ! - The matrix A is overwritten during the computation.
    ! - May not be suitable for very large-scale problems due to computational and memory constraints.
    !
    ! Usage:
    ! ------
    !  Given an m x n matrix A and an m-dimensional vector b, finds the n-dimensional
    !  vector x that minimizes the residual ||Ax - b||.
    !
    ! Input Parameters:
    ! -----------------
    !  A : Input matrix                    [Input]
    !  b : Right-hand side vector          [Input]
    !
    ! Output Parameters:
    ! ------------------
    !  x : Solution vector                 [Output]
    !
    ! References:
    ! -----------
    ! - Lawson, C. L., & Hanson, R. J. (1974). "Solving Least Squares Problems," Prentice-Hall, Inc.
    !
    !==============================================================================
    subroutine lstsq(A, b, x)
        !> Input matrix.
        real(kind=wp), dimension(:, :), intent(in)  :: A
        real(kind=wp), dimension(:), intent(in)  :: b
        real(kind=wp), dimension(:), intent(out) :: x

        !> Lapack job.
        character :: trans = "N"
        integer   :: m, n, nrhs, lda, ldb, lwork, info
        real(kind=wp), dimension(size(A, 1), size(A, 2)) :: A_tilde
        real(kind=wp), dimension(size(A, 1))             :: b_tilde
        real(kind=wp), dimension(:), allocatable         :: work

        !> Interface to LAPACK dgels
        interface
            pure subroutine dgels(ftrans, fm, fn, fnrhs, fA, flda, fb, fldb, fwork, flwork, finfo)
                import wp
                character, intent(in)           :: ftrans
                integer, intent(in)           :: fm, fn, fnrhs, flda, fldb, flwork, finfo
                real(kind=wp), intent(inout)    :: fa(flda, *)
                real(kind=wp), intent(inout)    :: fb(flda, *)
                real(kind=wp), intent(out)      :: fwork(*)
            end subroutine dgels
        end interface

        !> Initialize variables.
        m = size(A, 1); n = size(A, 2); nrhs = 1
        lda = m; ldb = m; lwork = max(1, min(m, n) + max(min(m, n), nrhs))
        A_tilde = A; b_tilde = b
        allocate (work(1:lwork)); work = 0.0_wp

        !> Solve the least-squares problem.
        call dgels(trans, m, n, nrhs, A_tilde, lda, b_tilde, ldb, work, lwork, info)

        !> Return solution.
        x = b_tilde(1:n)

        return
    end subroutine lstsq

    !=======================================================================================
    ! Conjugate Gradient (CG) Solver Subroutine
    !=======================================================================================
    !
    ! Purpose:
    ! --------
    ! Solves symmetric positive definite (SPD) linear systems Ax = b using the Conjugate Gradient (CG) method.
    !
    ! Mathematical Description:
    ! -------------------------
    ! Given a SPD matrix A and a right-hand side vector b, the aim is to find x that satisfies Ax = b.
    ! The algorithm starts with an initial guess x and computes the residual r = b - Ax.
    ! It then initializes a direction vector p = r. The residual dot product r_dot_r_old is calculated as r' * r.
    ! The algorithm iteratively updates the solution x and the residual r using the formulas:
    !   alpha = r_dot_r_old / (p' * Ap)
    !   x = x + alpha * p
    !   r = r - alpha * Ap
    !   r_dot_r_new = r' * r
    !   beta = r_dot_r_new / r_dot_r_old
    !   p = r + beta * p
    !
    ! Algorithmic Features:
    ! ----------------------
    ! - Relies on the method of conjugate directions for iteratively updating the solution.
    ! - Utilizes two sets of vectors: residuals (r) and conjugate directions (p) for iterations.
    ! - Each iteration computes new step size alpha and updates direction beta based on residuals.
    !
    ! Advantages:
    ! -----------
    ! - Highly efficient for large, sparse SPD matrices.
    ! - Minimal memory overhead; retains only a few vectors, in contrast to methods like GMRES.
    ! - Guarantees exact solution within 'n' iterations for an 'n'-dimensional SPD matrix under exact arithmetic.
    !
    ! Limitations:
    ! ------------
    ! - Exclusively for SPD matrices; not generalized for other types of matrices.
    ! - No built-in preconditioning.
    ! - Prone to round-off errors, sometimes requiring more than 'n' iterations.
    !
    ! Input/Output Parameters:
    ! ------------------------
    ! - A        : Linear Operator (SPD)                [Input]
    ! - b        : Right-hand side vector               [Input]
    ! - x        : Initial/Updated solution vector      [Input/Output]
    ! - info     : Iteration Information flag           [Output]
    ! - maxiter  : Maximum iterations (Optional)        [Input]
    ! - tol      : Convergence tolerance (Optional)     [Input]
    ! - verbosity: Verbosity control flag (Optional)    [Input]
    !
    ! References:
    ! -----------
    ! - Hestenes, M. R., and Stiefel, E. (1952). "Methods of Conjugate Gradients for Solving Linear Systems,"
    !   Journal of Research of the National Bureau of Standards, 49(6), 409–436.
    !
    !=======================================================================================
    subroutine cg(A, b, x, info, maxiter, tol, verbosity)
        !> Linear problem.
        class(abstract_spd_linop), intent(in) :: A ! Linear Operator.
        class(abstract_vector), intent(in) :: b ! Right-hand side.
        !> Solution vector.
        class(abstract_vector), intent(inout) :: x
        !> Information flag.
        integer, intent(out)   :: info
        !> Optional arguments.
        integer, optional, intent(in) :: maxiter ! Maximum number of CG iterations.
        integer                       :: i, niter
        real(kind=wp), optional, intent(in) :: tol  ! Tolerance for the CG residual.
        real(kind=wp)                       :: tolerance
        logical, optional, intent(in)       :: verbosity ! Verbosity control.
        logical                             :: verbose

        !> Residual and direction vectors.
        class(abstract_vector), allocatable :: r, p, Ap
        !> Scalars used in the CG algorithm.
        real(kind=wp) :: alpha, beta, r_dot_r_old, r_dot_r_new, residual

        ! --> Handle optional arguments.
        niter = optval(maxiter, 100)
        tolerance = optval(tol, atol + rtol*b%norm())
        verbose = optval(verbosity, .false.)

        ! --> Initialize vectors.
        allocate (r, source=b); call r%zero()
        allocate (p, source=b); call p%zero()
        allocate (Ap, source=b); call Ap%zero()

        ! --> Compute initial residual: r = b - Ax.
        call A%matvec(x, r); call r%axpby(-1.0_wp, b, 1.0_wp)

        ! --> Initialize direction vector: p = r.
        p = r

        ! --> Initialize dot product of residual: r_dot_r_old = r' * r.
        r_dot_r_old = r%dot(r)

        ! --> CG Iteration Loop.
        cg_iterations: do i = 1, niter

            ! Compute A * p.
            call A%matvec(p, Ap)

            ! Compute step size alpha = r_dot_r_old / (p' * Ap).
            alpha = r_dot_r_old/p%dot(Ap)

            ! Update solution x = x + alpha * p.
            call x%axpby(1.0_wp, p, alpha)

            ! Update residual r = r - alpha * Ap.
            call r%axpby(1.0_wp, Ap, -alpha)

            ! Compute new dot product of residual r_dot_r_new = r' * r.
            r_dot_r_new = r%dot(r)

            ! Check for convergence.
            residual = sqrt(r_dot_r_new)

            if (verbose) then
                write (*, *) "INFO : CG residual after ", (i), "iterations : ", residual
            end if

            if (residual < tolerance) then
                if (verbose) then
                    write (*, *) "INFO : CG Converged: residual ", residual, "< tolerance: ", tolerance
                end if
                exit cg_iterations
            end if

            ! Compute new direction beta = r_dot_r_new / r_dot_r_old.
            beta = r_dot_r_new/r_dot_r_old

            ! Update direction p = r + beta * p.
            call p%axpby(beta, r, 1.0_wp)

            ! Update r_dot_r_old for next iteration.
            r_dot_r_old = r_dot_r_new

        end do cg_iterations

        ! ! --> Set info flag.
        ! if (residual < tolerance) then
        !     info = 0
        ! else
        !     info = 1
        ! end if
        deallocate (r, p, Ap)

        return
    end subroutine cg

    !=======================================================================================
    ! Biconjugate Gradient Stabilized (BiCGSTAB) Solver Subroutine
    !=======================================================================================
    !
    ! Purpose:
    ! --------
    ! Solves nonsymmetric and potentially ill-conditioned linear systems Ax = b using the Biconjugate Gradient Stabilized (BiCGSTAB) method.
    !
    ! Mathematical Description:
    ! -------------------------
    ! Given a linear operator A and a right-hand side vector b, the goal is to find x that satisfies Ax = b.
    ! The algorithm starts with an initial guess x and computes the residual r = b - Ax. The residual vector r_hat is initialized to r.
    ! The algorithm iteratively performs the following steps:
    !   1. Compute v = A * p
    !   2. Compute alpha = rho / (r_hat' * v)
    !   3. Update s = r - alpha * v
    !   4. Compute t = A * s
    !   5. Compute omega = (t' * s) / (t' * t)
    !   6. Update x = x + omega * s + alpha * p
    !   7. Update r = s - omega * t
    !   8. Compute rho_new = r_hat' * r and beta = (alpha / omega) * (rho_new / rho)
    !   9. Update p = r + beta * (p - omega * v)
    !   10. Update rho = rho_new
    !
    ! Algorithmic Features:
    ! ----------------------
    ! - Extends the BiCG algorithm with a stabilization step, using two search directions and two residuals.
    ! - Iteratively updates the approximate solution x, and the residuals r and s.
    !
    ! Advantages:
    ! -----------
    ! - Effective for nonsymmetric and ill-conditioned matrices.
    ! - Generally faster and more stable than basic BiCG.
    ! - Suitable for large, sparse matrices.
    !
    ! Limitations:
    ! ------------
    ! - May experience stagnation for specific problem types.
    ! - No built-in preconditioning capabilities.
    !
    ! Input/Output Parameters:
    ! ------------------------
    ! - A        : Linear Operator [Input]
    ! - b        : Right-hand side vector [Input]
    ! - x        : Initial/Updated solution vector [Input/Output]
    ! - info     : Iteration Information flag [Output]
    ! - maxiter  : Maximum iterations (Optional) [Input]
    ! - tol      : Convergence tolerance (Optional) [Input]
    ! - verbosity: Verbosity control flag (Optional) [Input]
    ! - transpose: Transpose flag (Optional) [Input]
    !
    ! References:
    ! -----------
    ! - van der Vorst, H. A. (1992). "Bi-CGSTAB: A Fast and Smoothly Converging Variant of Bi-CG for the Solution of Nonsymmetric Linear Systems,"
    !   SIAM Journal on Scientific and Statistical Computing, 13(2), 631–644.
    !
    !=======================================================================================
    subroutine bicgstab(A, b, x, info, maxiter, tol, verbosity, transpose)
        !> Linear problem and initial guess.
        class(abstract_linop), intent(in) :: A
        class(abstract_vector), intent(in) :: b
        class(abstract_vector), intent(inout) :: x
        !> Output and optional input parameters.
        integer, intent(out) :: info
        integer, optional, intent(in) :: maxiter
        real(kind=wp), optional, intent(in) :: tol
        logical, optional, intent(in) :: verbosity
        logical, optional, intent(in) :: transpose

        !> Internal variables.
        integer :: i, niter
        real(kind=wp) :: tolerance, res, alpha, omega, rho, rho_new, beta
        logical :: verbose, trans

        !> BiCGSTAB vectors.
        class(abstract_vector), allocatable :: r, r_hat, p, v, s, t

        ! Initialize optional parameters.
        niter = optval(maxiter, 100)
        tolerance = optval(tol, atol + rtol*b%norm())
        verbose = optval(verbosity, .false.)
        trans = optval(transpose, .false.)
        info = 0 

        ! Initialize vectors.
        allocate (r, source=b)
        allocate (r_hat, source=b)
        allocate (p, source=b)
        allocate (v, source=b)
        allocate (s, source=b)
        allocate (t, source=b)

        ! --> Compute initial residual: r = b - Ax.
        if (trans) then
            call A%rmatvec(x, r)
        else
            call A%matvec(x, r)
        end if
        call r%axpby(-1.0_wp, b, 1.0_wp)

        !call r_hat%copy(r)
        r_hat = r

        rho = r_hat%dot(r)

        !call p%copy(r)
        p = r

        bicgstab_loop: do i = 1, niter

            ! v = A * p
            if (trans) then
                call A%rmatvec(p, v)
            else
                call A%matvec(p, v)
            end if

            alpha = rho/r_hat%dot(v)

            ! s = r - alpha * v
            !call s%copy(r)
            s = r
            call s%axpby(1.0_wp, v, -alpha)

            ! t = A * s
            if (trans) then
                call A%rmatvec(s, t)
            else
                call A%matvec(s, t)
            end if
            omega = t%dot(s)/t%dot(t)

            ! x = x + s * omega + p * alpha
            call x%axpby(1.0_wp, s, omega)
            call x%axpby(1.0_wp, p, alpha)

            ! r = s - t * omega
            !call r%copy(s)
            r = s
            call r%axpby(1.0_wp, t, -omega)

            res = r%norm()
            if (verbose) then
                write (*, *) "INFO : BICGSTAB residual after ", (i), "iterations : ", res
            end if
            if (res < tolerance) then
               info = 0 ! Algorithm converged
               exit bicgstab_loop
            endif
            rho_new = r_hat%dot(r)
            beta = (alpha/omega)*(rho_new/rho)

            ! s = p - v * omega ! reusing s vector
            !call s%copy(p)
            s = p
            call s%axpby(1.0_wp, v, -omega)

            ! p = r + s * beta
            !call p%copy(r)
            p = r
            call p%axpby(1.0_wp, s, beta)

            rho = rho_new

        end do bicgstab_loop

        deallocate (r, r_hat, p, v, s, t)
        return
    end subroutine bicgstab

    ! --> Utility Functions -----

    subroutine initialize_krylov_basis(X)
        ! Initializes Krylov basis vectors to zero, except for
        ! the first vector which is provided by the user
        class(abstract_vector), dimension(:), intent(inout) :: X
        integer :: i
        do i = 2, size(X)
            call X(i)%zero()
        end do
    end subroutine initialize_krylov_basis

end module IterativeSolvers