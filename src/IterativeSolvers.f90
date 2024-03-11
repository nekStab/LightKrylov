! This module is dedicated to iterative solvers for linear equations and eigenvalue problems.
! It imports several modules to handle vectors, linear operations, and Krylov decompositions.
! Additionally, it imports modules from the standard library.
module lightkrylov_IterativeSolvers
   !> LightKrylov modules.
   use lightkrylov_Utils
   use lightkrylov_AbstractVector
   use lightkrylov_LinearOperator
   use lightkrylov_BaseKrylov

   !> Fortran standard library.
   use stdlib_sorting, only: sort_index, int_size
   use stdlib_optval, only: optval
   use stdlib_io_npy, only: save_npy

   implicit none
   include "dtypes.h"

   private
   !> Eigenvalue analysis.
   public :: eigs, eighs, save_eigenspectrum
   !> Singular value decomposition.
   public :: svds
   !> Linear solvers.
   public :: gmres, cg, abstract_linear_solver

   !------------------------------------------------
   !-----                                      -----
   !-----     ABSTRACT PRECONDITIONER TYPE     -----
   !-----                                      -----
   !------------------------------------------------

   ! --> Abstract type definition.
   type, abstract, public :: abstract_preconditioner
   contains
      private
      procedure(abstract_apply_precond), deferred, pass(self), public :: apply
      procedure(abstract_undo_precond), deferred, pass(self), public :: undo
   end type abstract_preconditioner

   ! --> Definition of the type-bound procedures interfaces.
   abstract interface
      !> Apply the preconditionner.
      subroutine abstract_apply_precond(self, vec_inout)
         import abstract_preconditioner, abstract_vector
         class(abstract_preconditioner), intent(in)    :: self
         class(abstract_vector), intent(inout) :: vec_inout
      end subroutine abstract_apply_precond

      !> Undo the action of the preconditioner.
      subroutine abstract_undo_precond(self, vec_inout)
         import abstract_preconditioner, abstract_vector
         class(abstract_preconditioner), intent(in)    :: self
         class(abstract_vector), intent(inout) :: vec_inout
      end subroutine abstract_undo_precond
   end interface

   !--------------------------------------------------------
   !-----                                              -----
   !-----     GENERIC INTERFACE FOR LINEAR SOLVERS     -----
   !-----                                              -----
   !--------------------------------------------------------

   abstract interface
      subroutine abstract_linear_solver(A, b, x, info, preconditioner, options, transpose)
         import abstract_linop, abstract_vector, abstract_opts, abstract_preconditioner
         !> Linear problem.
         class(abstract_linop), intent(in) :: A
         class(abstract_vector), intent(in) :: b
         !> Solution vector.
         class(abstract_vector), intent(inout) :: x
         !> Information flag.
         integer, intent(out) :: info
         !> Solver options.
         class(abstract_opts), optional, intent(in) :: options
         !> Preconditioner.
         class(abstract_preconditioner), optional, intent(in) :: preconditioner
         !> Transposition flag.
         logical, optional, intent(in) :: transpose
      end subroutine abstract_linear_solver
   end interface

contains

   !=======================================================================================
   ! Compute Residual for Eigenpairs (compute_residual)
   !=======================================================================================
   !
   ! Purpose and formulation:
   ! -----------------------
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

   !=======================================================================================
   ! Eigenvalue and Eigenvector Solver Subroutine (EIGS)
   !=======================================================================================
   !
   ! Purpose:
   ! --------
   ! Computes the leading eigenpairs of a linear operator A using its Arnoldi factorization.
   !
   ! Mathematical Formulation:
   ! -------------------------
   ! Given a square linear operator A of size (n x n), find eigvals and eigvecs such that:
   !
   ! A * eigvecs = eigvals * eigvecs
   !
   ! or
   !
   ! transpose(A) * eigvecs = eigvals * eigvecs
   !
   ! The Krylov subspace X is constructed via Arnoldi factorization, resulting in an upper
   ! Hessenberg matrix H. The eigenvalues of A are approximated by the eigenvalues of H,
   ! and the eigenvectors are computed accordingly.
   !
   ! Algorithmic Features:
   ! ---------------------
   ! - Builds a Krylov subspace (X) using Arnoldi factorization.
   ! - Computes eigenpairs of the reduced upper Hessenberg matrix (H).
   ! - Sorts eigvals and associated eigvecs based on magnitude.
   !
   ! Advantages:
   ! -----------
   ! - Suitable for a wide range of square matrices A, including non-symmetric ones.
   ! - Efficient for large-scale problems.
   ! - Allows user-specified initial Krylov vector (X(1)).
   !
   ! Limitations:
   ! ------------
   ! - Accuracy is dependent on the quality of the initial Krylov vector.
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
   ! - Arnoldi, W. E. (1951). "The Principle of Minimized Iterations in the Solution of
   !   the Matrix Eigenvalue Problem." Quarterly of Applied Mathematics, 9(1), 17–29.
   !
   !=======================================================================================
   subroutine eigs(A, X, eigvals, residuals, info, nev, tolerance, verbosity, transpose)
      !> Linear Operator.
      class(abstract_linop), intent(in) :: A
      !> Krylov basis.
      class(abstract_vector), intent(inout) :: X(:)
      class(abstract_vector), allocatable   :: Xwrk(:)
      !> Coordinates of eigenvectors in Krylov basis, eigenvalues and associated residuals.
      real(kind=wp), allocatable :: eigvecs(:, :)
      complex(kind=wp), intent(out) :: eigvals(:)
      real(kind=wp), intent(out) :: residuals(:)
      !> Information flag.
      integer, intent(out) :: info
      !> Number of desired eigenvalues.
      integer, optional, intent(in) :: nev
      integer                       :: nev_, conv
      !> Tolerance control.
      real(kind=wp), optional, intent(in) :: tolerance
      real(kind=wp)                       :: tol
      !> Verbosity control.
      logical, optional, intent(in) :: verbosity
      logical :: verbose
      !> Transpose operator.
      logical, optional, intent(in) :: transpose
      logical :: trans

      !> Upper Hessenberg matrix.
      real(kind=wp) :: H(size(X), size(X) - 1)
      real(kind=wp) :: beta
      !> Krylov subspace dimension.
      integer :: kdim
      !> Miscellaneous.
      integer :: i, j, k
      integer(int_size) :: indices(size(X) - 1)
      real(kind=wp)     :: abs_eigvals(size(X) - 1)
      real(kind=wp)     :: alpha

      ! --> Dimension of the Krylov subspace.
      kdim = size(X) - 1

      ! --> Deals with the optional arguments.
      verbose = optval(verbosity, .false.)
      trans = optval(transpose, .false.)
      nev_ = optval(nev, size(X) - 1)
      tol = optval(tolerance, rtol)

      ! --> Initialize variables.
      allocate(eigvecs(kdim, kdim))
      H = 0.0_wp; residuals = 0.0_wp; eigvals = cmplx(0.0_wp, 0.0_wp, kind=wp); eigvecs = 0.0_wp
      !> Make sure the first Krylov vector has unit-norm.
      alpha = X(1)%norm(); call X(1)%scal(1.0_wp/alpha)
      call initialize_krylov_subspace(X(2:kdim + 1))

      arnoldi: do k = 1, kdim
         !> Arnoldi step.
         call arnoldi_factorization(A, X, H, info, kstart=k, kend=k, verbosity=verbose, transpose=trans)

         if (info < 0) then
            if (verbose) then
               write (*, *) "INFO : Arnoldi iteration failed. Exiting the eig subroutine."
               write (*, *) "       Arnoldi exit code :", info
            end if
            info = -1
            return
         end if

         !> Spectral decomposition of the k x k Hessenberg matrix.
         eigvals = 0.0_wp; eigvecs = cmplx(0.0_wp, 0.0_wp, kind=wp)
         call eig(H(1:k, 1:k), eigvecs(1:k, 1:k), eigvals(1:k))

         !> Sort eigenvalues.
         abs_eigvals = abs(eigvals); call sort_index(abs_eigvals, indices, reverse=.true.)
         eigvals = eigvals(indices); eigvecs = eigvecs(:, indices)

         !> Compute residuals.
         beta = H(k + 1, k) !> Get Krylov residual vector norm.
         residuals(1:k) = compute_residual(beta, abs(eigvecs(k, 1:k)))

         !> Check convergence.
         conv = count(residuals(1:k) < tol)
         if (conv >= nev_) then
            if (verbose) then
               write (*, *) "INFO : The first ", conv, "eigenpairs have converged."
               write (*, *) "       Exiting the computation."
            end if
            exit arnoldi
         end if

      end do arnoldi

      !> Reconstruct the eigenvectors from the Krylov basis.
      allocate(Xwrk, source=X) ; call mat_mult(X(1:kdim), Xwrk(1:kdim), eigvecs)

      return
   end subroutine eigs

   !=======================================================================================
   ! Eigenvalue and Eigenvector Solver for Symmetric Positive Definite (SPD) Matrices
   !=======================================================================================
   !
   ! Purpose:
   ! --------
   ! Computes the leading eigenpairs of a symmetric positive definite operator A using the
   ! iterative Lanczos factorization.
   !
   ! Algorithmic Features:
   ! ---------------------
   ! - Employs Lanczos tridiagonalization to construct a Krylov subspace represented by
   !   matrix T of size (kdim x kdim).
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
   ! - X         : Eigenvectors                   [Input/Output]
   ! - eigvals   : Eigenvalues                    [Output]
   ! - residuals : Residuals of eigenpairs        [Output]
   ! - info      : Iteration Information flag     [Output]
   ! - verbosity : Verbosity control flag         [Optional, Input]
   !
   ! References:
   ! -----------
   ! - Lanczos, C. (1950). "An Iteration Method for the Solution of the Eigenvalue Problem
   !   of Linear Differential and Integral Operators". United States Governm. Press Office.
   !
   !=======================================================================================
   subroutine eighs(A, X, eigvals, residuals, info, nev, tolerance, verbosity)
      !> Linear Operator.
      class(abstract_spd_linop), intent(in) :: A
      !> Krylov basis.
      class(abstract_vector), intent(inout) :: X(:)
      class(abstract_vector), allocatable   :: Xwrk(:)
      !> Coordinates of eigenvectors in Krylov basis, eigenvalues, and associated residuals
      real(kind=wp), allocatable :: eigvecs(:, :)
      real(kind=wp), intent(out) :: eigvals(:)
      real(kind=wp), intent(out) :: residuals(:)
      !> Information flag.
      integer, intent(out) :: info
      !> Number of converged eigenvalues needed.
      integer, optional, intent(in) :: nev
      integer                       :: nev_, conv
      !> Tolerance control.
      real(kind=wp), optional, intent(in) :: tolerance
      real(kind=wp)                       :: tol
      !> Verbosity control.
      logical, optional, intent(in) :: verbosity
      logical :: verbose

      !> Tridiagonal matrix.
      real(kind=wp), dimension(size(X), size(X) - 1) :: T
      real(kind=wp)                                :: beta
      !> Krylov subspace dimension.
      integer :: kdim
      !> Miscellaneous.
      integer :: i, j, k
      integer(int_size), dimension(size(X) - 1) :: indices
      real(kind=wp) :: alpha

      ! --> Dimension of the Krylov subspace.
      kdim = size(X) - 1

      ! --> Deals with the optional argument.
      verbose = optval(verbosity, .false.)
      nev_ = optval(nev, size(X) - 1)
      tol = optval(tolerance, rtol)

      ! --> Initialize all variables.
      allocate(eigvecs(1:kdim, 1:kdim))
      T = 0.0_wp; residuals = 0.0_wp; eigvecs = 0.0_wp; eigvals = 0.0_wp
      !> Make sure the first Krylov vector has unit-norm.
      alpha = X(1)%norm(); call X(1)%scal(1.0_wp/alpha)
      call initialize_krylov_subspace(X(2:kdim + 1))

      lanczos: do k = 1, kdim
         ! --> Compute Lanczos tridiagonalization.
         call lanczos_tridiagonalization(A, X, T, info, kstart=k, kend=k, verbosity=verbose)

         if (info < 0) then
            if (verbose) then
               write (*, *) "INFO : Lanczos iteration failed. Exiting the eigh subroutine."
               write (*, *) "       Lanczos exit code :", info
            end if
            info = -1
            return
         end if

         ! --> Compute spectral decomposition of the tridiagonal matrix.
         eigvals = 0.0_wp; eigvecs = 0.0_wp
         call eigh(T(1:k, 1:k), eigvecs(1:k, 1:k), eigvals(1:k))

         ! --> Sort eigenvalues in decreasing order.
         call sort_index(eigvals, indices, reverse=.true.); eigvecs = eigvecs(:, indices)

         ! --> Compute the residual associated with each eigenpair.
         beta = T(k + 1, k) !> Get Krylov residual vector norm.
         residuals(1:k) = compute_residual(beta, eigvecs(k, 1:k))

         !> Check convergence.
         conv = count(residuals(1:k) < tol)
         if (conv >= nev_) then
            if (verbose) then
               write (*, *) "INFO : The first", conv, "eigenpairs have converged."
               write (*, *) "Exiting the computation."
            end if
            exit lanczos
         end if

      end do lanczos

      !> Compute and returns the eigenvectors constructed from the Krylov basis.
      allocate(Xwrk, source=X) ; call mat_mult(X(1:kdim), Xwrk(1:kdim), eigvecs)

      return
   end subroutine eighs

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
   ! ---------------------
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
   ! - Only a full re-orthogonalization scheme is implemented.
   !
   ! Input/Output Parameters:
   ! ------------------------
   ! - A           : Linear Operator                                                        [Input]
   ! - U, V        : Left (U) and right (V) singular vectors                                [Input/Output]
   ! - sigma       : Singular values                                                        [Output]
   ! - residuals   : Residuals associated with singular triplets                            [Output]
   ! - info        : Information flag                                                       [Output]
   ! - verbosity   : Control for runtime diagnostics                                        [Optional, Input]
   !
   ! References:
   ! -----------
   ! - Golub, G. H., & Kahan, W. (1965). "Calculating the Singular Values and Pseudo-Inverse of a Matrix."
   ! - Baglama, J., & Reichel, L. (2005). "Augmented implicitly restarted Lanczos bidiagonalization methods."
   ! - R. M. Larsen. "Lanczos bidiagonalization with partial reorthogonalization." Technical Report, 1998.
   !   url : http://sun.stanford.edu/~rmunk/PROPACK/paper.pdf
   !
   !=======================================================================================
   subroutine svds(A, U, V, sigma, residuals, info, nev, tolerance, verbosity)
      !> Linear Operator.
      class(abstract_linop), intent(in) :: A
      !> Krylov bases.
      class(abstract_vector), intent(inout) :: U(:) ! Left sing. vectors.
      class(abstract_vector), intent(inout) :: V(:) ! Right sing. vectors.
      class(abstract_vector), allocatable   :: Uwrk(:), Vwrk(:)
      !> Coordinates of singular vectors in Krylov bases, singular values, and associated residuals.
      real(kind=wp), allocatable :: uvecs(:, :), vvecs(:, :)
      real(kind=wp), intent(out) :: sigma(:)
      real(kind=wp), intent(out) :: residuals(:)
      !> Information flag.
      integer, intent(out) :: info
      !> Number of converged singular triplets.
      integer, optional, intent(in) :: nev
      integer                       :: nev_, conv
      !> Tolerance control.
      real(kind=wp), optional, intent(in) :: tolerance
      real(kind=wp)                       :: tol
      !> Verbosity control.
      logical, optional, intent(in) :: verbosity
      logical verbose

      !> Bidiagonal matrix.
      real(kind=wp) :: B(size(U), size(U) - 1)
      real(kind=wp) :: beta
      !> Krylov subspace dimension.
      integer :: kdim
      !> Miscellaneous.
      integer :: i, j, k
      integer(int_size) :: indices(size(U) - 1)

      ! --> Deals with the optional args.
      nev_ = optval(nev, size(U) - 1)
      tol = optval(tolerance, rtol)
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
      B = 0.0_wp; residuals = 0.0_wp; sigma = 0.0_wp
      !> Make sure the first Krylov vector has unit-norm.
      beta = U(1)%norm(); call U(1)%scal(1.0_wp/beta)
      call initialize_krylov_subspace(U(2:kdim + 1))
      call initialize_krylov_subspace(V(2:kdim + 1))

      allocate(uvecs(kdim, kdim)) ; uvecs = 0.0_wp
      allocate(vvecs(kdim, kdim)) ; vvecs = 0.0_wp

      lanczos: do k = 1, kdim
         ! --> Compute the Lanczos bidiagonalization.
         call lanczos_bidiagonalization(A, U, V, B, info, kstart=k, kend=k, verbosity=verbose)

         ! --> Compute the singular value decomposition of the bidiagonal matrix.
         sigma = 0.0_wp; uvecs = 0.0_wp; vvecs = 0.0_wp
         call svd(B(1:k, 1:k), uvecs(1:k, 1:k), sigma(1:k), vvecs(1:k, 1:k))

         ! --> Compute the residual associated with each singular triplet.
         beta = B(k + 1, k) !> Get Krylov residual vector norm.
         residuals(1:k) = compute_residual(beta, vvecs(k, 1:k))

         !> Check convergence.
         conv = count(residuals(1:k) < tol)
         if (conv >= nev_) then
            if (verbose) then
               write (*, *) "INFO : The first ", conv, "singular triplets have converged."
               write (*, *) "       Exiting the computation."
            end if
            exit lanczos
         end if

      end do lanczos

      info = k

      !> Compute and return the low-rank factors from the Krylov bases.
      allocate(Uwrk, source=U) ; allocate(Vwrk, source=V)
      call mat_mult(U(1:kdim), Uwrk(1:kdim), uvecs) ; call mat_mult(V(1:kdim), Vwrk(1:kdim), vvecs)

      return
   end subroutine svds

   !=======================================================================================
   ! Generalized Minimal Residual (GMRES) Solver Subroutine
   !=======================================================================================
   !
   ! Purpose:
   ! --------
   ! Implements the classic, unpreconditioned Generalized Minimal Residual (GMRES) algorithm
   ! for solving nonsymmetric, non-Hermitian linear systems of equations, Ax = b.
   !
   ! Algorithmic Features:
   ! ---------------------
   ! - Constructs a full Krylov subspace without restarts (i.e., not GMRES(m)).
   ! - Utilizes Arnoldi factorization to generate an orthonormal basis for the Krylov subspace.
   ! - Employs a least-squares solve to determine the optimal linear combination of the Krylov vectors.
   ! - Updates the approximate solution based on the least-squares solution.
   !
   ! Advantages:
   ! -----------
   ! - Suitable for nonsymmetric and ill-conditioned matrices.
   ! - Produces monotonically decreasing residuals.
   !
   ! Limitations:
   ! ------------
   ! - Can be memory-intensive due to the absence of restarts.
   ! - Original gmres uses Givens rotations to update the least-squares solution. For the
   !   sake of simplicity, we use directly the lstsq solver from LAPACK.
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
   ! - Saad, Y., and Schultz, M. H. (1986). "GMRES: A Generalized Minimal Residual Algorithm
   !   for Solving Nonsymmetric Linear Systems," SIAM Journal on Scientific and Statistical
   !   Computing, 7(3), 856–869.
   !
   !=======================================================================================
   subroutine gmres(A, b, x, info, preconditioner, options, transpose)
      !> Linear problem.
      class(abstract_linop), intent(in)    :: A ! Linear Operator.
      class(abstract_vector), intent(in)    :: b ! Right-hand side.
      !> Solution vector.
      class(abstract_vector), intent(inout) :: x
      !> Information flag.
      integer, intent(out)   :: info
      !> Preconditioner.
      class(abstract_preconditioner), optional, intent(in) :: preconditioner
      class(abstract_preconditioner), allocatable          :: precond
      logical                                              :: has_precond
      !> Optional arguments.
      class(abstract_opts), optional, intent(in) :: options
      type(gmres_opts)                           :: opts
      logical, optional, intent(in) :: transpose
      logical                                    :: trans

      !> Internal variables.
      integer                                :: k_dim
      integer                                :: maxiter
      real(kind=wp)                          :: tol
      logical                                :: verbose

      !> Krylov subspace.
      class(abstract_vector), allocatable :: V(:)
      !> Upper Hessenberg matrix.
      real(kind=wp), allocatable :: H(:, :)
      !> Least-squares related variables.
      real(kind=wp), allocatable :: y(:)
      real(kind=wp), allocatable :: e(:)
      real(kind=wp)                       :: beta
      !> Miscellaneous.
      integer                             :: i, j, k, l, m
      real(kind=wp)                       :: alpha
      class(abstract_vector), allocatable :: dx, wrk

      ! --> Deals with the optional arguments.
      if (present(options)) then
         select type (options)
         type is (gmres_opts)
            opts = gmres_opts( &
                   kdim=options%kdim, &
                   maxiter=options%maxiter, &
                   atol=options%atol, &
                   rtol=options%rtol, &
                   verbose=options%verbose &
                   )
         end select
      else
         opts = gmres_opts()
      end if
      k_dim = opts%kdim; maxiter = opts%maxiter
      tol = opts%atol + opts%rtol*b%norm(); verbose = opts%verbose
      trans = optval(transpose, .false.)

      if (present(preconditioner)) then
         allocate (precond, source=preconditioner)
         has_precond = .true.
      else
         has_precond = .false.
      end if

      ! --> Initialize Krylov subspace.
      allocate (wrk, source=b); call wrk%zero()
      allocate (V(1:k_dim + 1), source=b)
      call initialize_krylov_subspace(V)
      allocate (H(k_dim + 1, k_dim)); H = 0.0_wp
      allocate (y(1:k_dim)); y = 0.0_wp
      allocate (e(1:k_dim + 1)); e = 0.0_wp

      info = 0

      ! --> Initial Krylov vector.
      if (trans) then
         call A%rmatvec(x, V(1))
      else
         call A%matvec(x, V(1))
      end if
      call V(1)%sub(b); call V(1)%scal(-1.0_wp)
      beta = V(1)%norm(); call V(1)%scal(1.0_wp/beta)
      gmres_iterations: do i = 1, maxiter
         ! --> Zero-out variables.
         H = 0.0_wp; y = 0.0_wp; e = 0.0_wp; e(1) = beta
         call initialize_krylov_subspace(V(2:k_dim + 1))

         ! --> Arnoldi factorization.
         arnoldi: do k = 1, k_dim
            wrk = V(k); if (has_precond) call precond%apply(wrk)

            ! --> Matrix-vector product.
            if (trans) then
               call A%rmatvec(wrk, V(k + 1))
            else
               call A%matvec(wrk, V(k + 1))
            end if

            ! --> Gram-Schmid orthogonalization (twice is enough).
            do j = 1, k
               alpha = V(k + 1)%dot(V(j)); call V(k + 1)%axpby(1.0_wp, V(j), -alpha)
               H(j, k) = alpha
            end do
            do j = 1, k
               alpha = V(k + 1)%dot(V(j)); call V(k + 1)%axpby(1.0_wp, V(j), -alpha)
               H(j, k) = H(j, k) + alpha
            end do

            ! --> Update Hessenberg matrix and normalize new Krylov vector.
            H(k + 1, k) = V(k + 1)%norm()
            if (H(k + 1, k) > tol) then
               call V(k + 1)%scal(1.0_wp/H(k + 1, k))
            end if

            ! --> Least-squares problem.
            call lstsq(H(1:k + 1, 1:k), e(1:k + 1), y(1:k))

            ! --> Compute residual.
            beta = norm2(e(1:k + 1) - matmul(H(1:k + 1, 1:k), y(1:k)))
            if (verbose) then
               write (*, *) "INFO : GMRES residual after ", info + 1, "iteration : ", beta
            end if

            ! --> Current number of iterations performed.
            info = info + 1

            ! --> Check convergence.
            if (beta .lt. tol) then
               exit arnoldi
            end if
         end do arnoldi

         ! --> Update solution.
         k = min(k, k_dim)
         if (allocated(dx) .eqv. .false.) allocate (dx, source=x); call dx%zero()
         call get_vec(dx, V(1:k), y(1:k)); if (has_precond) call precond%apply(dx)
         call x%add(dx)

         ! --> Recompute residual for sanity check.
         if (trans) then
            call A%rmatvec(x, V(1))
         else
            call A%matvec(x, V(1))
         end if

         call V(1)%sub(b); call V(1)%scal(-1.0_wp)

         ! --> Initialize new starting Krylov vector if needed.
         beta = V(1)%norm(); call V(1)%scal(1.0_wp/beta)

         ! --> Exit GMRES if desired accuracy is reached.
         if (beta .lt. tol) exit gmres_iterations

      end do gmres_iterations

      ! --> Convergence information.
      if (verbose) then
         write (*, *) "INFO : GMRES converged with residual ", beta
         write (*, *) "       Computation required ", info, "iterations"
      end if

      return
   end subroutine gmres

   !=======================================================================================
   ! Conjugate Gradient (CG) Solver Subroutine
   !=======================================================================================
   !
   ! Purpose:
   ! --------
   ! Implements the classic Conjugate Gradient (CG) algorithm for solving symmetric positive
   ! definite (SPD) linear systems of equations, Ax = b.
   !
   ! Algorithmic Features:
   ! ---------------------
   ! - Utilizes the method of conjugate directions to iteratively refine the solution.
   ! - Employs two sequences of vectors: residuals (r) and conjugate directions (p).
   ! - Updates the approximate solution based on the computed alpha and beta values.
   !
   ! Advantages:
   ! -----------
   ! - Well-suited for large, sparse, symmetric positive definite (SPD) matrices.
   ! - Memory-efficient, requiring storage for only a few vectors in comparison to GMRES.
   ! - Under exact arithmetic, finds the exact solution within 'n' iterations for an
   !   n-dimensional SPD matrix.
   !
   ! Limitations:
   ! ------------
   ! - Applicability restricted to SPD matrices.
   ! - Might take a large number of iterations for ill-conditioned problems.
   !
   ! Input/Output Parameters:
   ! ------------------------
   ! - A        : Linear Operator (SPD)        [Input]
   ! - b        : Right-hand side vector       [Input]
   ! - x        : Initial/Updated solution     [Input/Output]
   ! - info     : Iteration Information flag   [Output]
   ! - maxiter  : Maximum number of iterations [Optional, Input]
   ! - tol      : Tolerance for convergence    [Optional, Input]
   ! - verbosity: Verbosity control flag       [Optional, Input]
   !
   ! References:
   ! -----------
   ! - Hestenes, M. R., and Stiefel, E. (1952). "Methods of Conjugate Gradients for Solving
   !   Linear Systems," Journal of Research of the National Bureau of Standards, 49(6), 409–436.
   !
   !=======================================================================================
   subroutine cg(A, b, x, info, preconditioner, options)
      !> Linear problem.
      class(abstract_spd_linop), intent(in) :: A ! Linear Operator.
      class(abstract_vector), intent(in) :: b ! Right-hand side.
      !> Solution vector.
      class(abstract_vector), intent(inout) :: x
      !> Information flag.
      integer, intent(out)   :: info
      !> Preconditioner.
      class(abstract_preconditioner), optional, intent(in) :: preconditioner
      class(abstract_preconditioner), allocatable          :: precond
      logical                                              :: has_precond
      !> Optional arguments.
      class(abstract_opts), optional, intent(in) :: options
      type(cg_opts)                              :: opts

      integer :: maxiter ! Maximum number of CG iterations.
      real(kind=wp) :: tol  ! Tolerance for the CG residual.
      logical                             :: verbose

      !> Residual and direction vectors.
      class(abstract_vector), allocatable :: r, p, Ap
      !> Scalars used in the CG algorithm.
      real(kind=wp) :: alpha, beta, r_dot_r_old, r_dot_r_new, residual
      integer :: i, j, k

      ! --> Handle optional arguments.
      if (present(options)) then
         select type (options)
         type is (cg_opts)
            opts = cg_opts( &
                   maxiter=options%maxiter, &
                   atol=options%atol, &
                   rtol=options%rtol, &
                   verbose=options%verbose &
                   )
         end select
      else
         opts = cg_opts()
      end if
      tol = opts%atol + opts%rtol*b%norm(); maxiter = opts%maxiter; verbose = opts%verbose

      if (present(preconditioner)) then
         write (*, *) "INFO: CG does not support preconditioning yet. Precond is thus ignored."
         write (*, *)
      end if

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
      cg_iterations: do i = 1, maxiter

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

         if (residual < tol) then
            if (verbose) then
               write (*, *) "INFO : CG Converged: residual ", residual, "< tolerance: ", tol
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

end module lightkrylov_IterativeSolvers
