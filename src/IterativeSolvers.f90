module lightkrylov_IterativeSolvers
  !! This module provides some of the most important computational routines provided by `LightKrylov`.
  !! These include:
  !!
  !! - `eigs` : Compute the leading eigenpairs of a square linear operator \(\mathbf{A}\).
  !! - `eighs` : Compute the leading eigenpairs of a symmetric positive definite operator \(\mathbf{A}\).
  !! - `svds` : Compute the leading singular triplets of a linear operator \(\mathbf{A}\).
  !! - `gmres` : Solve the linear system \( \mathbf{Ax} = \mathbf{b}\) using the *generalized minimum
  !!    residual method*.
  !! - `cg` : Solve the linear system \( \mathbf{Ax} = \mathbf{b} \) where \(\mathbf{A}\) is symmetric
  !!   positive definite using the *Conjugate Gradient* method.
  !!
  !! It also provides abstract interfaces to pass user-defined linear solvers and preconditioners to
  !! `LightKrylov`. Note that these features are still experimental however.
   use iso_fortran_env, only: output_unit

   ! LightKrylov modules.
   use lightkrylov_Utils
   use lightkrylov_AbstractVector
   use lightkrylov_LinearOperator
   use lightkrylov_BaseKrylov

   ! Fortran standard library.
   use stdlib_sorting, only: sort_index, int_size
   use stdlib_optval, only: optval
   use stdlib_io_npy, only: save_npy

   implicit none
   include "dtypes.h"

   private
   ! Eigenvalue analysis.
   public :: eigs, eighs, save_eigenspectrum
   ! Singular value decomposition.
   public :: svds
   ! Linear solvers.
   public :: gmres, cg, abstract_linear_solver

   !------------------------------------------------
   !-----                                      -----
   !-----     ABSTRACT PRECONDITIONER TYPE     -----
   !-----                                      -----
   !------------------------------------------------

   type, abstract, public :: abstract_preconditioner
      !! Abstract type that needs to be extended to define a preconditioner usable by `LightKrylov`.
      !! Upon extension, the user needs to provide a type-bound procedure `apply`.
   contains
      private
      procedure(abstract_apply_precond), deferred, pass(self), public :: apply
      !! In-place application of the preconditioner.
      procedure(abstract_undo_precond), deferred, pass(self), public :: undo
      !! Deprecated.
   end type abstract_preconditioner

   abstract interface
      subroutine abstract_apply_precond(self, vec_inout)
        !! Abstract interface to apply a preconditioner in `LightKrylov`.
         import abstract_preconditioner, abstract_vector
         class(abstract_preconditioner), intent(in) :: self
         !! Preconditioner.
         class(abstract_vector), intent(inout) :: vec_inout
         !! Input/Output vector.
      end subroutine abstract_apply_precond

      subroutine abstract_undo_precond(self, vec_inout)
        !! Deprecated.
         import abstract_preconditioner, abstract_vector
         class(abstract_preconditioner), intent(in) :: self
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
        !! Abstract interface to use a user-defined linear solver in `LightKrylov`.
         import abstract_linop, abstract_vector, abstract_opts, abstract_preconditioner
         class(abstract_linop), intent(in) :: A
         !! Linear operator to invert.
         class(abstract_vector), intent(in) :: b
         !! Right-hand side vector.
         class(abstract_vector), intent(inout) :: x
         !! Solution vector.
         integer, intent(out) :: info
         !! Information flag.
         class(abstract_opts), optional, intent(in) :: options
         !! Options passed to the linear solver.
         class(abstract_preconditioner), optional, intent(in) :: preconditioner
         !! Preconditioner.
         logical, optional, intent(in) :: transpose
         !! Determine whether \(\mathbf{A}\) (`.false.`) or \(\mathbf{A}^T\) (`.true.`) is being used.
      end subroutine abstract_linear_solver
   end interface

   !----------------------------------------
   !-----                              -----
   !-----     REFINED RITZ VECTORS     -----
   !-----                              -----
   !----------------------------------------

   interface refined_ritz_vector
      module procedure refined_real_ritz_vector
      module procedure refined_cmplx_ritz_vector
   end interface refined_ritz_vector

contains

   elemental pure function compute_residual(beta, x) result(residual)
      !! Computes the residual associated with an eigenpair, useful for assessing the quality
      !! of the approximated eigenvalues and eigenvectors. Given the norm of the Krylov residual
      !! vector (beta) and the last element of the Ritz eigenvector (x), the residual is calculated
      !! as follows \( r = \vert \beta x \vert \).
      real(kind=wp), intent(in) :: beta
      !! Norm of the residual Krylov vector.
      real(kind=wp), intent(in) :: x
      !! Last entry of the low-dimensional Ritz eigenvector.
      real(kind=wp) :: residual
      !! Residual associated to the corresponding Ritz eigenpair.
      residual = abs(beta*x)
      return
   end function compute_residual

   subroutine save_eigenspectrum(real_part, imag_part, residuals, filename)
     !! Saves the eigenvalues and corresponding residuals to disk using the `npy` binary format
     !! to be compatible with `Python`.
      real(kind=wp), intent(in) :: real_part(:)
      !! Real parts of the eigenvalues.
      real(kind=wp), intent(in) :: imag_part(:)
      !! Imaginary parts of the eigenvalues.
      real(kind=wp), intent(in) :: residuals(:)
      !! Residual of the corresponding Ritz eigenpairs.
      character(len=*), intent(in) :: filename
      !! Name of the output file.

      ! Miscellaneous.
      real(kind=wp), dimension(size(real_part), 3) :: data

      ! Store the data.
      data(:, 1) = real_part; data(:, 2) = imag_part; data(:, 3) = residuals
      ! Save the eigenspectrum to disk using npy file format.
      call save_npy(filename, data)

      return
   end subroutine save_eigenspectrum

   function refined_real_ritz_vector(H, theta) result(x)
      real(kind=wp), intent(in) :: H(:, :)
      !! Matrix resulting from the Krylov factorization.
      real(kind=wp), intent(in) :: theta
      !! Eigenvalue for which refinement needs to be done.
      real(kind=wp)             :: x(size(H, 2))
      !! Low-dimensional refined Ritz eigenvector.

      ! Internal variables.
      real(kind=wp) :: A(size(H, 1), size(H, 2))
      real(kind=wp) :: U(size(H, 1), size(H, 1)), V(size(H, 2), size(H, 2))
      real(kind=wp) :: S(size(H, 2))
      integer       :: i

      ! Build H - theta*Id.
      A = H
      do i = 1, size(A, 2)
         A(i, i) = A(i, i) - theta
      end do

      ! USV.T = svd(A).
      call svd(A, U, S, V)

      ! Refined Ritz vector.
      x = V(:, size(H, 2))

      return
   end function refined_real_ritz_vector

   function refined_cmplx_ritz_vector(H, theta) result(x)
      real(kind=wp), intent(in) :: H(:, :)
      !! Matrix resulting from the Krylov factorization.
      complex(kind=wp), intent(in) :: theta
      !! Eigenvalue for which refinement needs to be done.
      complex(kind=wp)             :: x(size(H, 2))
      !! Low-dimensional refined Ritz eigenvector.

      ! Internal variables.
      complex(kind=wp) :: A(size(H, 1), size(H, 2))
      complex(kind=wp) :: U(size(H, 1), size(H, 1)), V(size(H, 2), size(H, 2))
      real(kind=wp)    :: S(size(H, 2))
      integer          :: i

      ! Build H - theta*Id.
      A = H
      do i = 1, size(A, 2)
         A(i, i) = A(i, i) - theta
      end do

      ! USV.H = svd(A).
      call svd(A, U, S, V)

      ! Refined Ritz vector.
      x = V(:, size(H, 2))

      return
   end function refined_cmplx_ritz_vector

   subroutine eigs(A, X, eigvals, residuals, info, nev, tolerance, verbosity, transpose)
     !! Computes the leading eigenpairs of a linear operator A using its Arnoldi factorization.
     !! Given a square linear operator A of size (n x n), find eigvals and eigvecs such that:
     !!
     !! \[
     !!    \mathbf{Ax} = \lambda \mathbf{x}
     !! \]
     !!
     !! or
     !!
     !! \[
     !!    \mathbf{A}^T \mathbf{x} = \lambda \mathbf{x}.
     !! \]
     !!
     !! The Krylov subspace X is constructed via Arnoldi factorization, resulting in an upper
     !! Hessenberg matrix H. The eigenvalues of A are approximated by the eigenvalues of H,
     !! and the eigenvectors are computed accordingly.
     !!
     !! **Algorithmic Features**
     !!
     !! - Builds a Krylov subspace (X) using Arnoldi factorization.
     !! - Computes eigenpairs of the reduced upper Hessenberg matrix (H).
     !! - Sorts eigvals and associated eigvecs based on magnitude.
     !!
     !! **Advantages**
     !!
     !! - Suitable for a wide range of square matrices A, including non-symmetric ones.
     !! - Efficient for large-scale problems.
     !! - Allows user-specified initial Krylov vector (X(1)).
     !!
     !! **Limitations**
     !!
     !! - Accuracy is dependent on the quality of the initial Krylov vector.
     !! - No preconditioning capabilities in the current implementation.
     !!
     !! **References**
     !!
     !! - Arnoldi, W. E. (1951). "The Principle of Minimized Iterations in the Solution of
     !!   the Matrix Eigenvalue Problem." Quarterly of Applied Mathematics, 9(1), 17–29.
      class(abstract_linop), intent(in) :: A
      !! Linear operator whose leading eigenpairs need to be computed.
      class(abstract_vector), intent(inout) :: X(:)
      !! Leading eigenvectors of \(\mathbf{A}\). On entry, `X(1)` needs to be set to the starting Krylov vector.
      !! On exit, storage of the real and imaginary parts of the eigenvectors follows the LAPACK convention.
      complex(kind=wp), intent(out) :: eigvals(:)
      !! Leading eigenvalues of \( \mathbf{A} \).
      real(kind=wp), intent(out) :: residuals(:)
      !! Residual associated with each computed Ritz eigenpair.
      integer, intent(out) :: info
      !! Information flag.
      integer, optional, intent(in) :: nev
      !! Desired number of eigenpairs.
      real(kind=wp), optional, intent(in) :: tolerance
      !! Tolerance to determine whether a Ritz eigenpair has converged or not (default `sqrt(epsilon(1.0_wp))`).
      logical, optional, intent(in) :: verbosity
      !! Verbosity control (default `.false.`)
      logical, optional, intent(in) :: transpose
      !! Determine whether \(\mathbf{A}\) (default `.false.`) or \( \mathbf{A}^T\) (`.true.`) is used.

      ! Internal variables.
      class(abstract_vector), allocatable   :: Xwrk(:)
      real(kind=wp), allocatable :: eigvecs(:, :)
      integer                       :: nev_, conv
      real(kind=wp)                       :: tol
      logical :: verbose
      real(kind=wp) :: H(size(X), size(X) - 1)
      real(kind=wp) :: beta
      integer :: kdim
      integer :: i, j, k
      integer(int_size) :: indices(size(X) - 1)
      real(kind=wp)     :: abs_eigvals(size(X) - 1)
      real(kind=wp)     :: alpha
      logical :: trans

      ! Dimension of the Krylov subspace.
      kdim = size(X) - 1

      ! Deals with the optional arguments.
      verbose = optval(verbosity, .false.)
      trans = optval(transpose, .false.)
      nev_ = optval(nev, size(X) - 1)
      tol = optval(tolerance, rtol)

      ! Initialize variables.
      allocate (eigvecs(kdim, kdim))
      H = 0.0_wp; residuals = 0.0_wp; eigvals = cmplx(0.0_wp, 0.0_wp, kind=wp); eigvecs = 0.0_wp
      ! Make sure the first Krylov vector has unit-norm.
      alpha = X(1)%norm(); call X(1)%scal(1.0_wp/alpha)
      call initialize_krylov_subspace(X(2:kdim + 1))

      arnoldi: do k = 1, kdim
         ! Arnoldi step.
         call arnoldi_factorization(A, X, H, info, kstart=k, kend=k, verbosity=verbose, transpose=trans)

         if (info < 0) then
            if (verbose) then
               write (output_unit, *) "INFO : Arnoldi iteration failed. Exiting the eig subroutine."
               write (output_unit, *) "       Arnoldi exit code :", info
            end if
            info = -1
            return
         end if

         ! Spectral decomposition of the k x k Hessenberg matrix.
         eigvals = 0.0_wp; eigvecs = cmplx(0.0_wp, 0.0_wp, kind=wp)
         call eig(H(1:k, 1:k), eigvecs(1:k, 1:k), eigvals(1:k))

         ! Sort eigenvalues.
         abs_eigvals = abs(eigvals); call sort_index(abs_eigvals, indices, reverse=.true.)
         eigvals = eigvals(indices); eigvecs = eigvecs(:, indices)

         ! Compute residuals.
         beta = H(k + 1, k) ! Get Krylov residual vector norm.
         residuals(1:k) = compute_residual(beta, abs(eigvecs(k, 1:k)))

         ! Check convergence.
         conv = count(residuals(1:k) < tol)
         if (conv >= nev_) then
            if (verbose) then
               write (output_unit, *) "INFO : The first ", conv, "eigenpairs have converged."
               write (output_unit, *) "       Exiting the computation."
            end if
            exit arnoldi
         end if

      end do arnoldi

      ! Reconstruct the eigenvectors from the Krylov basis.
      allocate (Xwrk, source=X); call mat_mult(X(1:kdim), Xwrk(1:kdim), eigvecs)

      return
   end subroutine eigs

   subroutine eighs(A, X, eigvals, residuals, info, nev, tolerance, verbosity)
     !! Computes the leading eigenpairs of a symmetric positive definite operator \( \mathbf{A} \) using the
     !! iterative Lanczos factorization.
     !!
     !! **Algorithmic Features**
     !!
     !! - Employs Lanczos tridiagonalization.
     !! - Sorts eigvals and eigvecs based on magnitude.
     !! - Computes residuals to assess the quality of approximations.
     !!
     !! **Advantages**
     !!
     !! - Efficient for large, sparse, sym. pos. def. matrices (A).
     !! - Lower memory requirements compared to full spectrum methods like QR.
     !! - Faster convergence rates for extreme eigvals.
     !!
     !! **Limitations**
     !!
     !! - Only applicable to sym. pos. def. matrices (A).
     !! - Eigenpairs are approximations and may require refinement.
     !! - Not suitable for finding all eigvals in large dense matrices.
     !!
     !! **References**
     !!
     !! - Lanczos, C. (1950). "An Iteration Method for the Solution of the Eigenvalue Problem
     !!   of Linear Differential and Integral Operators". United States Governm. Press Office.
      class(abstract_spd_linop), intent(in) :: A
      !! Linear operator whose leading eigenpairs need to be computed.
      class(abstract_vector), intent(inout) :: X(:)
      !! Leading eigenvectors of \(\mathbf{A}\). On entry, `X(1)` needs to be set to the starting
      !! Krylov vector.
      real(kind=wp), intent(out) :: eigvals(:)
      !! Leading eigenvalues of \(\mathbf{A}\).
      real(kind=wp), intent(out) :: residuals(:)
      !! Residual associated to each Ritz eigenpair computed.
      integer, intent(out) :: info
      !! Information flag.
      integer, optional, intent(in) :: nev
      !! Desired number of eigenpairs. If not specified, `k = size(X)` Lanczos iterations are done
      !! and all the eigenpairs are reported even if not converged.
      real(kind=wp), optional, intent(in) :: tolerance
      !! Tolerance to check convergence of a given eigenpair (default `sqrt(epsilon(1.0_wp))`).
      logical, optional, intent(in) :: verbosity
      !! Verbosity control (default `.false.`).

      ! Internal variables
      class(abstract_vector), allocatable   :: Xwrk(:)
      real(kind=wp), allocatable :: eigvecs(:, :)
      integer                       :: nev_, conv
      real(kind=wp)                       :: tol
      logical :: verbose
      real(kind=wp), dimension(size(X), size(X) - 1) :: T
      real(kind=wp)                                :: beta
      integer :: kdim
      integer :: i, j, k
      integer(int_size), dimension(size(X) - 1) :: indices
      real(kind=wp) :: alpha

      ! Dimension of the Krylov subspace.
      kdim = size(X) - 1

      ! Deals with the optional argument.
      verbose = optval(verbosity, .false.)
      nev_ = optval(nev, size(X) - 1)
      tol = optval(tolerance, rtol)

      ! Initialize all variables.
      allocate (eigvecs(1:kdim, 1:kdim))
      T = 0.0_wp; residuals = 0.0_wp; eigvecs = 0.0_wp; eigvals = 0.0_wp
      ! Make sure the first Krylov vector has unit-norm.
      alpha = X(1)%norm(); call X(1)%scal(1.0_wp/alpha)
      call initialize_krylov_subspace(X(2:kdim + 1))

      lanczos: do k = 1, kdim
         ! Compute Lanczos tridiagonalization.
         call lanczos_tridiagonalization(A, X, T, info, kstart=k, kend=k, verbosity=verbose)

         if (info < 0) then
            if (verbose) then
               write (output_unit, *) "INFO : Lanczos iteration failed. Exiting the eigh subroutine."
               write (output_unit, *) "       Lanczos exit code :", info
            end if
            info = -1
            return
         end if

         ! Compute spectral decomposition of the tridiagonal matrix.
         eigvals = 0.0_wp; eigvecs = 0.0_wp
         call eigh(T(1:k, 1:k), eigvecs(1:k, 1:k), eigvals(1:k))

         ! Sort eigenvalues in decreasing order.
         call sort_index(eigvals, indices, reverse=.true.); eigvecs = eigvecs(:, indices)

         ! Compute the residual associated with each eigenpair.
         beta = T(k + 1, k) ! Get Krylov residual vector norm.
         residuals(1:k) = compute_residual(beta, eigvecs(k, 1:k))

         ! Check convergence.
         conv = count(residuals(1:k) < tol)
         if (conv >= nev_) then
            if (verbose) then
               write (output_unit, *) "INFO : The first", conv, "eigenpairs have converged."
               write (output_unit, *) "Exiting the computation."
            end if
            exit lanczos
         end if

      end do lanczos

      ! Compute refined Ritz vectors.
      deallocate (eigvecs); allocate (eigvecs(k, conv)); eigvecs = 0.0_wp
      do i = 1, conv
         eigvecs(:, i) = refined_ritz_vector(T(1:k + 1, 1:k), eigvals(i))
      end do

      ! Compute and returns the eigenvectors constructed from the Krylov basis.
      allocate (Xwrk, source=X(1:k)); call initialize_krylov_subspace(X)
      call mat_mult(X(1:conv), Xwrk(1:k), eigvecs)

      return
   end subroutine eighs

   subroutine svds(A, U, V, sigma, residuals, info, nev, tolerance, verbosity)
     !! Computes the singular value decomposition of a given linear operator \(\mathbf{A}\) using Lanczos
     !! bidiagonalization. The subroutine finds matrices \(\mathbf{U}\), \(\boldsymbol{\Sigma}\), and \(\mathbf{V}\)
     !! such that
     !!
     !! \[
     !!    \mathbf{A} = \mathbf{U} \boldsymbol{\Sigma} \mathbf{V}^T.
     !! \]
     !!
     !! **Algorithmic Features**
     !!
     !! - Relies on Lanczos bidiagonalization.
     !! - Computes the singular value decomposition of of a small bidiagonal matrix to obtain approximate
     !!   singular values and vectors of \(\mathbf{A}\).
     !! - Calculates residuals for singular triplets.
     !! - Optionally verbose, providing runtime diagnostics.
     !!
     !! **Advantages**
     !!
     !! - Efficient for large, sparse matrices.
     !! - Suitable for problems requiring only a subset of singular values and vectors.
     !! - Can be adapted for preconditioning techniques.
     !!
     !! **Limitations**
     !!
     !! - Not optimized for dense matrices or those with special structures (e.g., Toeplitz, circulant).
     !! - May require many iterations for ill-conditioned matrices.
     !! - Only a full re-orthogonalization scheme is implemented.
     !!
     !! **References**
     !!
     !! - Golub, G. H., & Kahan, W. (1965). "Calculating the Singular Values and Pseudo-Inverse of a Matrix."
     !! - Baglama, J., & Reichel, L. (2005). "Augmented implicitly restarted Lanczos bidiagonalization methods."
     !! - R. M. Larsen. "Lanczos bidiagonalization with partial reorthogonalization." Technical Report, 1998.
     !!   [(PDF)](http://sun.stanford.edu/~rmunk/PROPACK/paper.pdf)
      class(abstract_linop), intent(in) :: A
      !! Linear operator whose leading singular triplets are to be computed.
      class(abstract_vector), intent(inout) :: U(:)
      !! Leading left singular vectors of \(\mathbf{A}\). On entry, `U(1)` needs to be set to the starting
      !! Krylov vector.
      class(abstract_vector), intent(inout) :: V(:)
      !! Leading right singular vectors of \(\mathbf{A}\).
      real(kind=wp), intent(out) :: sigma(:)
      !! Leading singular values.
      real(kind=wp), intent(out) :: residuals(:)
      !! Residual associated to each Ritz singular triplet.
      integer, intent(out) :: info
      !! Information flag.
      integer, optional, intent(in) :: nev
      !! Desired number of converged singular triplets.
      real(kind=wp), optional, intent(in) :: tolerance
      !! Tolerance to check convergence of each triplet.
      logical, optional, intent(in) :: verbosity
      !! Verbosity control.

      ! Internal variables.
      real(kind=wp) :: B(size(U), size(U) - 1)
      real(kind=wp) :: beta
      integer :: kdim
      integer :: i, j, k
      integer(int_size) :: indices(size(U) - 1)
      class(abstract_vector), allocatable   :: Uwrk(:), Vwrk(:)
      real(kind=wp), allocatable :: uvecs(:, :), vvecs(:, :)
      integer                       :: nev_, conv
      real(kind=wp)                       :: tol
      logical verbose

      ! Deals with the optional args.
      nev_ = optval(nev, size(U) - 1)
      tol = optval(tolerance, rtol)
      verbose = optval(verbosity, .false.)
      ! Assert size(U) == size(V).
      if (size(U) .ne. size(V)) then
         info = -1
         if (verbose) then
            write (output_unit, *) "INFO : Left and Right Krylov subspaces have different dimensions."
            write (output_unit, *) "       Exiting svds with exit code info =", info
         end if
      else
         kdim = size(U) - 1
      end if

      ! Initialize variables.
      B = 0.0_wp; residuals = 0.0_wp; sigma = 0.0_wp
      ! Make sure the first Krylov vector has unit-norm.
      beta = U(1)%norm(); call U(1)%scal(1.0_wp/beta)
      call initialize_krylov_subspace(U(2:kdim + 1))
      call initialize_krylov_subspace(V(2:kdim + 1))

      allocate (uvecs(kdim, kdim)); uvecs = 0.0_wp
      allocate (vvecs(kdim, kdim)); vvecs = 0.0_wp

      lanczos: do k = 1, kdim
         ! Compute the Lanczos bidiagonalization.
         call lanczos_bidiagonalization(A, U, V, B, info, kstart=k, kend=k, verbosity=verbose)

         ! --> Compute the singular value decomposition of the bidiagonal matrix.
         sigma = 0.0_wp; uvecs = 0.0_wp; vvecs = 0.0_wp
         call svd(B(1:k, 1:k), uvecs(1:k, 1:k), sigma(1:k), vvecs(1:k, 1:k))

         ! Compute the residual associated with each singular triplet.
         beta = B(k + 1, k) ! Get Krylov residual vector norm.
         residuals(1:k) = compute_residual(beta, vvecs(k, 1:k))

         ! Check convergence.
         conv = count(residuals(1:k) < tol)
         if (conv >= nev_) then
            if (verbose) then
               write (output_unit, *) "INFO : The first ", conv, "singular triplets have converged."
               write (output_unit, *) "       Exiting the computation."
            end if
            exit lanczos
         end if

      end do lanczos

      info = k

      ! Compute and return the low-rank factors from the Krylov bases.
      allocate (Uwrk, source=U); allocate (Vwrk, source=V)
      call mat_mult(U(1:kdim), Uwrk(1:kdim), uvecs); call mat_mult(V(1:kdim), Vwrk(1:kdim), vvecs)

      return
   end subroutine svds

   subroutine gmres(A, b, x, info, preconditioner, options, transpose)
     !! Implements the classic, unpreconditioned Generalized Minimal Residual (GMRES) algorithm
     !! for solving nonsymmetric linear systems of equations \( \mathbf{Ax} = \mathbf{b}\).
     !!
     !! **Algorithmic Features**
     !!
     !! - Constructs a size-m Krylov subspace.
     !! - Utilizes Arnoldi factorization to generate an orthonormal basis for the Krylov subspace.
     !! - Employs a least-squares solve to determine the optimal linear combination of the Krylov vectors.
     !! - Updates the approximate solution based on the least-squares solution.
     !!
     !! **Advantages*
     !!
     !! - Suitable for nonsymmetric and ill-conditioned matrices.
     !! - Produces monotonically decreasing residuals.
     !!
     !! **Limitations**
     !!
     !! - Original gmres uses Givens rotations to update the least-squares solution. For the
     !!   sake of simplicity, we use directly the lstsq solver from LAPACK.
     !!
     !! **References**
     !!
     !! - Saad, Y., and Schultz, M. H. (1986). "GMRES: A Generalized Minimal Residual Algorithm
     !!   for Solving Nonsymmetric Linear Systems," SIAM Journal on Scientific and Statistical
     !!   Computing, 7(3), 856–869.
      class(abstract_linop), intent(in)     :: A
      !! Linear operator to be inverted.
      class(abstract_vector), intent(in)    :: b
      !! Right-hand side vector.
      class(abstract_vector), intent(inout) :: x
      !! Solution vector.
      integer, intent(out)   :: info
      !! Information flag.
      class(abstract_preconditioner), optional, intent(in) :: preconditioner
      !! Preconditioner.
      class(abstract_opts), optional, intent(in) :: options
      !! Options.
      logical, optional, intent(in) :: transpose
      !! Whether \(\mathbf{A}\) (`.false.`) or \(\mathbf{A}^T\) (`.true.`) is used.

      ! Internal variables.
      integer                                :: k_dim
      integer                                :: maxiter
      real(kind=wp)                          :: tol
      logical                                :: verbose

      ! Krylov subspace.
      class(abstract_vector), allocatable :: V(:)
      ! Upper Hessenberg matrix.
      real(kind=wp), allocatable :: H(:, :)
      ! Least-squares related variables.
      real(kind=wp), allocatable :: y(:)
      real(kind=wp), allocatable :: e(:)
      real(kind=wp)                       :: beta
      ! Miscellaneous.
      integer                             :: i, j, k, l, m
      real(kind=wp)                       :: alpha
      class(abstract_vector), allocatable :: dx, wrk
      logical                                              :: has_precond
      class(abstract_preconditioner), allocatable          :: precond
      type(gmres_opts)                           :: opts
      logical                                    :: trans

      ! Deals with the optional arguments.
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

      ! Initialize Krylov subspace.
      allocate (wrk, source=b); call wrk%zero()
      allocate (V(1:k_dim + 1), source=b)
      call initialize_krylov_subspace(V)
      allocate (H(k_dim + 1, k_dim)); H = 0.0_wp
      allocate (y(1:k_dim)); y = 0.0_wp
      allocate (e(1:k_dim + 1)); e = 0.0_wp

      info = 0

      ! Initial Krylov vector.
      if (trans) then
         call A%rmatvec(x, V(1))
      else
         call A%matvec(x, V(1))
      end if
      call V(1)%sub(b); call V(1)%scal(-1.0_wp)
      beta = V(1)%norm(); call V(1)%scal(1.0_wp/beta)
      gmres_iterations: do i = 1, maxiter
         ! Zero-out variables.
         H = 0.0_wp; y = 0.0_wp; e = 0.0_wp; e(1) = beta
         call initialize_krylov_subspace(V(2:k_dim + 1))

         ! Arnoldi factorization.
         arnoldi: do k = 1, k_dim
            wrk = V(k); if (has_precond) call precond%apply(wrk)

            ! Matrix-vector product.
            if (trans) then
               call A%rmatvec(wrk, V(k + 1))
            else
               call A%matvec(wrk, V(k + 1))
            end if

            ! Gram-Schmid orthogonalization (twice is enough).
            do j = 1, k
               alpha = V(k + 1)%dot(V(j)); call V(k + 1)%axpby(1.0_wp, V(j), -alpha)
               H(j, k) = alpha
            end do
            do j = 1, k
               alpha = V(k + 1)%dot(V(j)); call V(k + 1)%axpby(1.0_wp, V(j), -alpha)
               H(j, k) = H(j, k) + alpha
            end do

            ! Update Hessenberg matrix and normalize new Krylov vector.
            H(k + 1, k) = V(k + 1)%norm()
            if (H(k + 1, k) > tol) then
               call V(k + 1)%scal(1.0_wp/H(k + 1, k))
            end if

            ! Least-squares problem.
            call lstsq(H(1:k + 1, 1:k), e(1:k + 1), y(1:k))

            ! Compute residual.
            beta = norm2(e(1:k + 1) - matmul(H(1:k + 1, 1:k), y(1:k)))
            if (verbose) then
               write (*, *) "INFO : GMRES residual after ", info + 1, "iteration : ", beta
            end if

            ! Current number of iterations performed.
            info = info + 1

            ! Check convergence.
            if (beta .lt. tol) then
               exit arnoldi
            end if
         end do arnoldi

         ! Update solution.
         k = min(k, k_dim)
         if (allocated(dx) .eqv. .false.) allocate (dx, source=x); call dx%zero()
         call get_vec(dx, V(1:k), y(1:k)); if (has_precond) call precond%apply(dx)
         call x%add(dx)

         ! Recompute residual for sanity check.
         if (trans) then
            call A%rmatvec(x, V(1))
         else
            call A%matvec(x, V(1))
         end if

         call V(1)%sub(b); call V(1)%scal(-1.0_wp)

         ! Initialize new starting Krylov vector if needed.
         beta = V(1)%norm(); call V(1)%scal(1.0_wp/beta)

         ! Exit GMRES if desired accuracy is reached.
         if (beta .lt. tol) exit gmres_iterations

      end do gmres_iterations

      ! Convergence information.
      if (verbose) then
         write (output_unit, *) "INFO : GMRES converged with residual ", beta
         write (output_unit, *) "       Computation required ", info, "iterations"
      end if

      return
   end subroutine gmres

   subroutine cg(A, b, x, info, preconditioner, options)
     !! Implements the classic Conjugate Gradient (CG) algorithm for solving symmetric positive
     !! definite linear systems of equations, \(\mathbf{Ax} = \mathbf{b}\).
     !!
     !! **Advantages**
     !!
     !! - Well-suited for large, sparse, symmetric positive definite (SPD) matrices.
     !! - Memory-efficient, requiring storage for only a few vectors in comparison to GMRES.
     !! - Under exact arithmetic, finds the exact solution within 'n' iterations for an
     !!   n-dimensional SPD matrix.
     !!
     !! **Limitations**
     !!
     !! - Applicability restricted to SPD matrices.
     !! - Might take a large number of iterations for ill-conditioned problems.
     !!
     !! **References**
     !!
     !! - Hestenes, M. R., and Stiefel, E. (1952). "Methods of Conjugate Gradients for Solving
     !!   Linear Systems," Journal of Research of the National Bureau of Standards, 49(6), 409–436.
      class(abstract_spd_linop), intent(in) :: A
      !! Linear operator to be inverted.
      class(abstract_vector), intent(in) :: b
      !! Right-hand side vector.
      class(abstract_vector), intent(inout) :: x
      !! Solution vector.
      integer, intent(out)   :: info
      !! Information flag.
      class(abstract_preconditioner), optional, intent(in) :: preconditioner
      !! Preconditioner (not yet supported).
      class(abstract_preconditioner), allocatable          :: precond
      logical                                              :: has_precond
      class(abstract_opts), optional, intent(in) :: options
      !! Options.
      type(cg_opts)                              :: opts

      integer       :: maxiter
      real(kind=wp) :: tol
      logical       :: verbose

      ! Residual and direction vectors.
      class(abstract_vector), allocatable :: r, p, Ap
      ! Scalars used in the CG algorithm.
      real(kind=wp) :: alpha, beta, r_dot_r_old, r_dot_r_new, residual
      integer :: i, j, k

      ! Handle optional arguments.
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
         write (output_unit, *) "INFO: CG does not support preconditioning yet. Precond is thus ignored."
         write (output_unit, *)
      end if

      ! Initialize vectors.
      allocate (r, source=b); call r%zero()
      allocate (p, source=b); call p%zero()
      allocate (Ap, source=b); call Ap%zero()

      ! Compute initial residual: r = b - Ax.
      call A%matvec(x, r); call r%axpby(-1.0_wp, b, 1.0_wp)

      ! Initialize direction vector: p = r.
      p = r

      ! Initialize dot product of residual: r_dot_r_old = r' * r.
      r_dot_r_old = r%dot(r)

      ! CG Iteration Loop.
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
            write (output_unit, *) "INFO : CG residual after ", (i), "iterations : ", residual
         end if

         if (residual < tol) then
            if (verbose) then
               write (output_unit, *) "INFO : CG Converged: residual ", residual, "< tolerance: ", tol
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

      deallocate (r, p, Ap)

      return
   end subroutine cg

end module lightkrylov_IterativeSolvers
