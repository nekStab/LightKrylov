module lightkrylov_BaseKrylov
  !! This module provides a large collection of Krylov-based factorizations forming the computational core
  !! of `LightKrylov`. These include:
  !!
  !! - the Arnoldi factorization for general square matrices,
  !! - the Lanczos factorization for symmetric matrices,
  !! - the Lanczos factorization for non-symmetric matrices,
  !! - the Lanczos bidiagonalization for arbitrary linear operators,
  !!
  !! It also provides a set of utility functions to initialize Krylov subspaces, perform Krylov-Schur restarts
  !! (only for Arnoldi for the moment) or orthogonalization of an array of Krylov vectors using the modified
  !! Gram-Schmid process. Note that the Arnoldi factorization supports both the standard and block
  !! implementations.

   ! Fortran intrinsic module.
   use iso_fortran_env

   ! LightKrylov modules.
   use lightkrylov_Utils
   use lightkrylov_AbstractVector
   use lightkrylov_LinearOperator

   ! Fortran standard library.
   use stdlib_optval, only: optval

   implicit none
   include "dtypes.h"

   private
   ! Factorization for general n x n linear operators.
   public :: arnoldi_factorization, krylov_schur_restart
   ! Factorization for symmetric positive definite operators.
   public :: lanczos_tridiagonalization
   ! Pre-factorization for singular value decomposition.
   public :: lanczos_bidiagonalization
   ! Experimental/Deprecated.
   public :: nonsymmetric_lanczos_tridiagonalization
   public :: two_sided_arnoldi_factorization
   ! Miscellaneous.
   public :: qr_factorization, initialize_krylov_subspace

contains

   subroutine initialize_krylov_subspace(X, B0)
    !! Utility function to initializes an array of `abstract_vector`. By default, it initializes `X` to an array of zero `abstract_vector`s.
    !! If `B0` is passed, it initializes the first `p` columns of `X` with the content of `B0` followed by orthogonalization.
      class(abstract_vector), intent(inout) :: X(:)
    !! Array of `abstract_vector` to initialize.
      class(abstract_vector), optional, intent(in) :: B0(:)
    !! If present, the span of the first `p` columns.

      ! Internal variables.
      class(abstract_vector), allocatable :: B(:)
      real(kind=wp), allocatable :: Rwrk(:, :)
      integer :: i, p, info

      ! Zero out X
      call mat_zero(X)

      ! Deals with the optional starting vector
      if (present(B0)) then
         p = size(B0)
         ! Sanity check
         if (size(X) .lt. p) then
            write (*, *) "ERROR : Mismatch between basis size and size of initial vector."
            stop 1
         end if
         ! Allocate & initialize
         allocate (B(1:p), source=B0(1:p))
         call mat_zero(B); call mat_copy(B, B0)
         ! Orthonormalize
         allocate (Rwrk(1:p, 1:p)); Rwrk = 0.0_wp
         call qr_factorization(B, Rwrk, info)
         ! Set initial vector
         call mat_copy(X(1:p), B)
      end if

      return
   end subroutine initialize_krylov_subspace

   subroutine arnoldi_factorization(A, X, H, info, kstart, kend, verbosity, tol, transpose, block_size)
      !! Given a linear operator \( \mathbf{A} \in \mathbb{R}^{n x n} \), find \(\mathbf{X}\) and \(\mathbf{H}\)
      !! such that:
      !! \[
      !!   \mathbf{AX}_k = \mathbf{X}_k \mathbf{H}_k + h_{k+1, k} \mathbf{x}_{k+1} \mathbf{e}_k
      !! \]
      !!
      !! where \( \mathbf{X} \) is an orthogonal basis and \(\mathbf{H}\) is an upper Hessenberg matrix.
      !!
      !! **Algorithmic Features**
      !!
      !! - The operator \(\mathbf{A}\) only needs to be access through matrix-vector products.
      !! - Constructs an orthonormal Krylov basis X via modified Gram-Schmidt.
      !! - Constructs an upper Hessenberg matrix H.
      !! - Checks for convergence and invariant subspaces.
      !!
      !! **Advantages**
      !!
      !! - Applicable for non-symmetric matrices.
      !! - Basis for many Krylov subspace methods.
      !!
      !! **Limitations**
      !!
      !! - Orthogonality of X may deteriorate for ill-conditioned matrices.
      !! - Not suitable for preconditioning in this implementation.
      !!
      !! **References**
      !!
      !! - Y. Saad. "Iterative methods for sparse linear systems", SIAM 2nd edition, 2003.
      !!   see Chapter 6.3 : Arnoldi's method.

      class(abstract_linop), intent(in) :: A
      !! Linear operator to be factorized.
      class(abstract_vector), intent(inout) :: X(:)
      !! Orthogonal Krylov subspace. Note that on entry `X(1)` needs to be set to the initial Krylov vector.
      real(kind=wp), intent(inout) :: H(:, :)
      !! Upper Hessenberg matrix.
      integer, intent(out) :: info
      !! Information flag. On exit:
      !!
      !! - `info` < 0 : The k-step Arnoldi factorization failed.
      !! - `info` = 0 : The k-step Arnoldi factorization succeeded.
      !! - `info` > 0 : An invariant subspace has been computed after `k=info` steps.
      integer, optional, intent(in) :: kstart
      !! Starting index for the Arnoldi factorization (default 1).
      integer, optional, intent(in) :: kend
      !! Final index for the Arnoldi iteration (default `size(X)-1`).
      logical, optional, intent(in) :: verbosity
      !! Verbosity control (default `.false.`).
      logical, optional, intent(in) :: transpose
      !! Whether \( \mathbf{A} \) is being transposed or not (default `.false.`).
      real(kind=wp), optional, intent(in) :: tol
      !! Tolerance to determine whether an invariant subspace has been computed (default `sqrt(epsilon(1.0_wp))`).
      integer, optional, intent(in) :: block_size
      !! Block size for block Arnoldi (default 1).

      ! Internal variables.
      integer :: k_start, k_end, p
      logical :: verbose, trans
      real(kind=wp) :: tolerance

      ! Internal variables.
      real(kind=wp) :: beta
      real(kind=wp), allocatable :: res(:)
      integer :: k, i, kdim, kpm, kp, kpp

      ! Deals with optional non-unity block size
      p = optval(block_size, 1)
      allocate (res(1:p))

      info = 0

      ! Check dimensions.
      kdim = (size(X) - p)/p; call assert_shape(H, [p*(kdim + 1), p*kdim], "arnoldi_factorization", "H")

      ! Deals with the optional arguments.
      k_start = optval(kstart, 1)
      k_end = optval(kend, kdim)
      verbose = optval(verbosity, .false.)
      tolerance = optval(tol, atol)
      trans = optval(transpose, .false.)

      ! Arnoldi factorization.
      block_arnoldi: do k = k_start, k_end
         ! Counters
         kpm = (k - 1)*p
         kp = kpm + p
         kpp = kp + p
         ! Matrix-vector products.
         if (trans) then
            do i = 1, p
               call A%rmatvec(X(kpm + i), X(kp + i))
            end do
         else
            do i = 1, p
               call A%matvec(X(kpm + i), X(kp + i))
            end do
         end if
         ! Update Hessenberg matrix w.r.t. existing vectors
         call update_hessenberg_matrix(H, X, k, p)
         ! Orthogonalize current vectors
         call qr_factorization(X(kp + 1:kpp), H(kp + 1:kpp, kpm + 1:kp), info)
         ! extract residual norm (smallest diagonal element of R matrix)
         res = 0.0_wp
         do i = 1, p
            res(i) = H(kp + i, kpm + i)
         end do
         beta = minval(res)

         if (verbose) then
            if (p .eq. 1) then
               write (output_unit, *) "--> Arnoldi iteration n°", k, "/", k_end
            else
               write (output_unit, *) "--> Block Arnoldi iteration n°", k, "/", k_end
               write (output_unit, *) "    n° of vectors in basis:", kp
            end if
            write (output_unit, *) "    -----------------"
            write (output_unit, *) "    + Residual norm :", beta
            write (output_unit, *) "    + Elapsed time  :"
            write (output_unit, *) "    + ETA           :"
            write (output_unit, *)
         end if

         ! Exit Arnoldi loop if needed.
         if (beta < tolerance) then
            if (verbose) then
               write (output_unit, *)
               write (output_unit, *) "INFO : An invariant subspace has been computed (beta =", beta, ")."
            end if

            ! Dimension of the computed invariant subspace.
            info = kp

            ! Exit the Arnoldi iteration.
            exit block_arnoldi
         end if

      end do block_arnoldi

      if (verbose) then
         write (output_unit, *) "INFO : Exiting the block Arnoldi factorization with exit code info =", info, "."
         write (output_unit, *)
      end if

      return
   end subroutine arnoldi_factorization

   subroutine update_hessenberg_matrix(H, X, k, block_size)
      integer, intent(in) :: k
      real(kind=wp), intent(inout) :: H(:, :)
      class(abstract_vector) :: X(:)
      ! Optional: size of blocks, default = 1
      integer, optional, intent(in) :: block_size
      ! Internal variables.
      class(abstract_vector), allocatable :: Xwrk(:)
      real(wp), allocatable :: wrk(:, :)
      integer :: p, kpm, kp, kpp

      ! Deals with optional non-unity block size
      p = optval(block_size, 1)

      kpm = (k - 1)*p
      kp = kpm + p
      kpp = kp + p
      allocate (wrk(1:kp, 1:p)); wrk = 0.0_wp
      allocate (Xwrk(1:kp), source=X(1:kp)); call mat_zero(Xwrk)
      ! Orthogonalize residual w.r.t to previously computed Krylov vectors.
      ! - Pass 1
      call mat_mult(H(1:kp, kpm + 1:kp), X(1:kp), X(kp + 1:kpp))
      call mat_mult(Xwrk(1:p), X(1:kp), H(1:kp, kpm + 1:kp))
      ! Project out existing vectors
      call mat_axpby(X(kp + 1:kpp), 1.0_wp, Xwrk(1:p), -1.0_wp)
      ! - Pass 2
      call mat_mult(wrk, X(1:kp), X(kp + 1:kpp))
      call mat_mult(Xwrk(1:p), X(1:kp), wrk)
      ! Project out existing vectors
      call mat_axpby(X(kp + 1:kpp), 1.0_wp, Xwrk(1:p), -1.0_wp)

      ! Update Hessenberg matrix with data from second pass
      call mat_axpby(H(1:kp, kpm + 1:kp), 1.0_wp, wrk, 1.0_wp)
      deallocate (wrk)

      return
   end subroutine update_hessenberg_matrix

   subroutine lanczos_tridiagonalization(A, X, T, info, kstart, kend, verbosity, tol)
     !! Given a symmetric positive definite operator \(\mathbf{A} \in \mathbb{R}^{n \times n}\), find
     !! an orthonormal matrix \(\mathbf{X}\) and symmetric tridiagonal matrix \(\mathbf{T}\) such that
     !! \[
     !!    \mathbf{AX}_k = \mathbf{X}_k \mathbf{T}_k + t_{k+1,k} \mathbf{x}_{k+1} \mathbf{e}_{k}^T.
     !! \]
     !!
     !! **Algorithmic Features**
     !!
     !! - The operator \(\mathbf{A}\) only needs to be accessed through matrix-vector products.
     !! - Constructs an orthonormal Krylov basis \(\mathbf{X}\).
     !! - Constructs a tridiagonal matrix \(\mathbf{T}\).
     !! - Checks for convergence and invariant subspaces.
     !!
     !! **Advantages**
     !!
     !! - Efficient for symmetric positive definite operators.
     !! - Foundation for eigenvalue and linear system solvers.
     !!
     !! **Limitations**
     !!
     !! - Limited to sym. pos. def. matrices.
     !! - Orthogonality may deteriorate for ill-conditioned matrices.
     !!
     !! **References**
     !!
     !! - Y. Saad. "Iterative methods for sparse linear systems", SIAM 2nd edition, 2003.
     !!   see Chapter 6.6 : The symmetric Lanczos algorithm.

      class(abstract_spd_linop), intent(in) :: A
      !! Linear operator to be factorized.
      class(abstract_vector), intent(inout) :: X(:)
      !! Orthogonal Krylov subspace. Note that on entry, `X(1)` needs to be set to the initial Krylov vector.
      real(kind=wp), intent(inout) :: T(:, :)
      !! Tridiagonal matrix.
      integer, intent(out) :: info
      !! Information flag. On exit:
      !!
      !! - `info` < 0 : The k-step Arnoldi factorization failed.
      !! - `info` = 0 : The k-step Arnoldi factorization succeeded.
      !! - `info` > 0 : An invariant subspace has been computed after `k=info` steps.
      integer, optional, intent(in) :: kstart
      !! Starting index for the Lanczos factorization (default 1).
      integer, optional, intent(in) :: kend
      !! Final index for the Lanczos factorization (default 1).
      logical, optional, intent(in) :: verbosity
      !! Verbosity control (default `.false.`).
      real(kind=wp), optional, intent(in) :: tol
      !! Tolerance to determine whether an invariant subspace has been computed (default `sqrt(epsilon(1.0_wp))`).

      ! Internal variables.
      integer :: k_start
      integer :: k_end
      logical :: verbose
      real(kind=wp) :: tolerance
      real(kind=wp) :: beta
      integer       :: k, kdim
      integer       :: i, j

      ! Check dimensions.
      kdim = size(X) - 1; call assert_shape(T, [kdim + 1, kdim], "lanczos_tridiag.", "T")

      ! Deals with the optional arguments.
      k_start = optval(kstart, 1)
      k_end = optval(kend, kdim)
      verbose = optval(verbosity, .false.)
      tolerance = optval(tol, atol)

      ! Lanczos tridiagonalization.
      lanczos: do k = k_start, k_end
         ! Matrix-vector product.
         call A%matvec(X(k), X(k + 1))
         ! Update tridiagonal matrix.
         call update_tridiag_matrix(T, X, k); beta = X(k + 1)%norm(); T(k + 1, k) = beta

         if (verbose) then
            write (output_unit, *) "--> Lanczos iteration n°", k, "/", k_end
            write (output_unit, *) "    -----------------"
            write (output_unit, *) "    + Residual norm :", beta
            write (output_unit, *) "    + Elapsed time  :"
            write (output_unit, *) "    + ETA           :"
         end if

         ! Exit Lanczos loop if needed.
         if (beta < tolerance) then
            if (verbose) then
               write (output_unit, *)
               write (output_unit, *) "INFO : An invariant subspace has been computed (beta =)", beta, ")."
            end if
            ! Dimension of the computed invariant subspace.
            info = k
            ! Exit the Lanczos iteration.
            exit lanczos
         else
            ! Normalize the new Krylov vector.
            call X(k + 1)%scal(1.0_wp/beta)
         end if
      end do lanczos

      if (verbose) then
         write (*, *) "INFO : Exiting the Lanczos factorization with exit code info =", info, "."
      end if

      return
   end subroutine lanczos_tridiagonalization

   subroutine update_tridiag_matrix(T, X, k)
      integer, intent(in) :: k
      real(kind=wp), intent(inout) :: T(:, :)
      class(abstract_vector) :: X(:)
      class(abstract_vector), allocatable :: wrk
      integer :: i
      real(kind=wp) :: alpha

      ! Orthogonalize residual w.r.t to previously computed Krylov vectors.
      do i = max(1, k - 1), k
         alpha = X(k + 1)%dot(X(i)); call X(k + 1)%axpby(1.0_wp, X(i), -alpha)
         ! Update tridiag matrix.
         T(i, k) = alpha
      end do
      ! Full re-orthogonalization.
      do i = 1, k
         alpha = X(k + 1)%dot(X(i)); call X(k + 1)%axpby(1.0_wp, X(i), -alpha)
      end do

      return
   end subroutine update_tridiag_matrix

   subroutine lanczos_bidiagonalization(A, U, V, B, info, kstart, kend, verbosity, tol)
     !! Performs Lanczos bidiagonalization on a given linear operator \(\mathbf{A}\), producing left and
     !! right Krylov bases \(\mathbf{U}\) and \(\mathbf{V}\), and a bidiagonal matrix \(\mathbf{B}\).
     !! Given \(\mathbf{A} \in \mathbb{R}^{m \times n}\), find \(\mathbf{U}\), \(\mathbf{V}\), and \(\mathbf{B}\)
     !! such that:
     !! \[
     !!   \mathbf{AV} = \mathbf{UB}.
     !! \]
     !!
     !! **Algorithmic Features**
     !!
     !! - Constructs orthonormal bases \(\mathbf{U}\) and \(\mathbf{V}\) for the dominant column span
     !!   and row span of \(\mathbf{A}\).
     !! - Constructs a bidiagonal matrix \(\mathbf{B}\).
     !!
     !! **References**
     !!
     !! - R. M. Larsen. "Lanczos bidiagonalization with partial reorthogonalization." Technical Report, 1998.
     !!   [(PDF)](http://sun.stanford.edu/~rmunk/PROPACK/paper.pdf)
      class(abstract_linop), intent(in) :: A
      !! Linear operator to be factorized.
      class(abstract_vector), intent(inout) :: U(:)
      !! Orthonormal basis for the column span of \(\mathbf{A}\). On entry, `U(1)` needs to be set to
      !! the starting Krylov vector.
      class(abstract_vector), intent(inout) :: V(:)
      !! Orthonormal basis for the row span of \(\mathbf{A}\).
      real(kind=wp), intent(inout) :: B(:, :)
      !! Bidiagonal matrix satisfying \( \mathbf{U}^T \mathbf{AV} = \mathbf{B} \).
      integer, intent(out) :: info
      !! Information flag. On exit:
      !!
      !! - `info` < 0 : The k-step Arnoldi factorization failed.
      !! - `info` = 0 : The k-step Arnoldi factorization succeeded.
      !! - `info` > 0 : An invariant subspace has been computed after `k=info` steps.
      integer, optional, intent(in) :: kstart
      !! Starting index for the Lanczos factorization (default 1).
      integer, optional, intent(in) :: kend
      !! Final index for the Lanczos factorization (default `size(U)-1`).
      logical, optional, intent(in) :: verbosity
      !! Verbosity control (default `.false.`).
      real(kind=wp), optional, intent(in) :: tol
      !! Tolerance to determine whether invariant subspaces have been computed (default `sqrt(epsilon(1.0_wp))`).

      ! Internal variables.
      integer                       :: k_start
      integer                       :: k_end
      logical                       :: verbose
      real(kind=wp)                       :: tolerance
      real(kind=wp) :: alpha, beta, gamma
      integer       :: i, j, k, kdim

      ! Check Krylov subspaces dimensions.
      if (size(U) .ne. size(V)) then
         write (output_unit, *) "INFO : Left and right Krylov basis have different dimensions."
         info = -1
         return
      else
         kdim = size(U) - 1
      end if

      ! Check B dimensions.
      call assert_shape(B, [kdim + 1, kdim], "lanczos_bidiag.", "B")

      ! Deals with the optional arguments.
      k_start = optval(kstart, 1)
      k_end = optval(kend, kdim)
      verbose = optval(verbosity, .false.)
      tolerance = optval(tol, atol)

      ! Lanczos bidiagonalization.
      lanczos: do k = k_start, k_end
         ! Transpose matrix-vector product.
         call A%rmatvec(U(k), V(k))
         ! /!\ NOTE : Next lines not needed because already taken care of in the
         !     full re-orthogonalization.
         ! if (k > 1) then
         !    call V(k)%axpby(1.0_wp, V(k-1), -beta)
         ! endif

         ! Full re-orthogonalization of the right Krylov subspace.
         do j = 1, k - 1
            gamma = V(k)%dot(V(j)); call V(k)%axpby(1.0_wp, V(j), -gamma)
         end do

         ! Normalization step.
         alpha = V(k)%norm(); B(k, k) = alpha
         if (alpha > tolerance) then
            call V(k)%scal(1.0_wp/alpha)
         else
            if (verbose) then
               write (output_unit, *) "INFO : alpha = ", alpha
            end if
            info = k
            exit lanczos
         end if

         ! Matrix-vector product.
         call A%matvec(V(k), U(k + 1))
         ! /!\ NOTE : Not needed because taken care of in the full reortho. step.
         ! call U(k+1)%axpby(1.0_wp, U(k), -alpha)

         ! Full re-orthogonalization of the left Krylov subspace.
         do j = 1, k
            gamma = U(k + 1)%dot(U(j)); call U(k + 1)%axpby(1.0_wp, U(j), -gamma)
         end do

         ! Normalization step.
         beta = U(k + 1)%norm(); B(k + 1, k) = beta
         if (beta > tolerance) then
            call U(k + 1)%scal(1.0_wp/beta)
         else
            if (verbose) then
               write (output_unit, *) "INFO : beta = ", beta
            end if
            info = k
            exit lanczos
         end if

      end do lanczos

      return
   end subroutine lanczos_bidiagonalization

   subroutine nonsymmetric_lanczos_tridiagonalization(A, V, W, T, info, kstart, kend, verbosity, tol)
      !! Performs Lanczos tridiagonalization on a given nonsymmetric linear operator A, producing
      !! left and right Krylov bases V and W, and a tridiagonal matrix T satisfying
      !! \[
      !!   \begin{aligned}
      !!   \mathbf{AV} & = \mathbf{WT} \\
      !!    \mathbf{A}^T \mathbf{W} & = \mathbf{VT}^T.
      !! \end{aligned}
      !! \]
      !!
      !! Note additionally that \(\mathbf{V}\) and \(\mathbf{W}\) are bi-orthogonal.
      !!
      !! **Algorithmic Features**
      !!
      !! - Constructs bi-orthogonal bases V and W for the domain and image of A.
      !! - Constructs a tridiagonal matrix T.
      !!
      !! **Advantages**
      !!
      !! - Simultaneously compute a vector basis for the domain and image of A.
      !! - Bi-orthogonality of the two bases is a feature of the algorithm.
      !!
      !! **Limitations**
      !!
      !! - Maintaining the bi-orthogonality for high-dimensional operators is numerically
      !!   unstable and quickly suffers from accumulation of floating point erros.
      !! - It is strongly suggested NOT to use it for production code.
      !!
      !! **References**
      !!
      !! - Y. Saad. "Iterative methods for sparse linear systems", SIAM 2nd edition, 2003.
      !!   see Chapter 7.1 : Lanczos biorthogonalization.
      class(abstract_linop), intent(in)    :: A
      !! Linear operator to be factorized.
      class(abstract_vector), intent(inout) :: V(:)
      !! Bi-orthogonal basis for the column span of \(\mathbf{A}\). On entry `V(1)` needs to be initialized
      !! with the starting Krylov vector for the column span.
      class(abstract_vector), intent(inout) :: W(:)
      !! Bi-orthogonal basis for the row span of \(\mathbf{A}\). On entry `W(1)` needs to be initialized with
      !! the starting Krylov vector for the row span.
      real(kind=wp), intent(inout) :: T(:, :)
      !! Tridiagonal matrix satisfying \(\mathbf{W}^T \mathbf{AV} = \mathbf{T}\).
      integer, intent(out)   :: info
      !! Information flag. On exit:
      !!
      !! - `info` < 0 : The k-step Arnoldi factorization failed.
      !! - `info` = 0 : The k-step Arnoldi factorization succeeded.
      !! - `info` > 0 : An invariant subspace has been computed after `k=info` steps.
      integer, optional, intent(in) :: kstart
      !! Starting index for the Lanczos factorization (default 1).
      integer, optional, intent(in) :: kend
      !! Final index for the Lanczos factorization (default `size(V)-1`).
      logical, optional, intent(in) :: verbosity
      !! Verbosity control (default `.false.`).
      real(kind=wp), optional, intent(in) :: tol
      !! Tolerance to determine whether invariant subspaces have been computed (default `sqrt(epsilon(1.0_wp))`).

      ! Internal variables.
      real(kind=wp) :: alpha, beta, gamma, tmp
      integer       :: i, j, k, kdim
      integer                       :: k_start
      integer                       :: k_end
      logical                       :: verbose
      real(kind=wp)                       :: tolerance

      ! Check Krylov subspaces dimensions.
      if (size(V) .ne. size(W)) then
         write (output_unit, *) "INFO : Left and right Krylov basis have different dimensions."
         info = -1
      else
         kdim = size(V) - 1
      end if

      ! Check T dimensions.
      if (all(shape(T) .ne. [kdim + 1, kdim + 1])) then
         write (output_unit, *) "INFO : Tridiagonal matrix and Krylov subspaces dimensions do not match."
         info = -2
      end if

      ! Deals with the optional arguments.
      k_start = optval(kstart, 1)
      k_end = optval(kend, kdim)
      verbose = optval(verbosity, .false.)
      tolerance = optval(tol, rtol)

      ! Bi-orthogonalize the left and right starting vectors.
      if (k_start .eq. 1) then
         tmp = V(1)%dot(W(1)); beta = sqrt(abs(tmp)); gamma = sign(beta, tmp)
         call V(1)%scal(1.0_wp/beta); call W(1)%scal(1.0_wp/gamma)
      end if

      ! Nonsymmetric Lanczos iterations.
      lanczos: do k = k_start, k_end
         ! --> Matrix-vector product.
         call A%matvec(V(k), V(k + 1))
         call A%rmatvec(W(k), W(k + 1))

         ! Update diagonal entry of the nonsymmetric tridiagonal matrix.
         alpha = V(k + 1)%dot(W(k)); T(k, k) = alpha

         ! Lanczos three term recurrence.
         call V(k + 1)%axpby(1.0_wp, V(k), -alpha)
         if (k > 1) call V(k + 1)%axpby(1.0_wp, V(k - 1), -T(k - 1, k))

         call W(k + 1)%axpby(1.0_wp, W(k), -alpha)
         if (k > 1) call W(k + 1)%axpby(1.0_wp, W(k - 1), -T(k, k - 1))

         ! Full re-biorthogonalization.
         do i = 1, k
            alpha = V(k + 1)%dot(W(i)); call V(k + 1)%axpby(1.0_wp, V(i), -alpha)
            alpha = W(k + 1)%dot(V(i)); call W(k + 1)%axpby(1.0_wp, W(i), -alpha)
         end do

         ! Update the off-diagonal entries of the nonsymmetric tridiagonal matrix.
         tmp = V(k + 1)%dot(W(k + 1)); beta = sqrt(abs(tmp)); gamma = sign(beta, tmp)
         T(k, k + 1) = gamma; T(k + 1, k) = beta

         if (abs(tmp) < tolerance) then
            if (verbose) then
               write (output_unit, *) "INFO : Invariant subspaces have been computed (beta, gamma) = (", beta, ",", gamma, ")."
            end if
         else
            ! Normalization step.
            call V(k + 1)%scal(1.0_wp/beta); call W(k + 1)%scal(1.0_wp/gamma)
         end if

      end do lanczos

      if (verbose) then
         write (output_unit, *) "INFO : Exiting the nonsymmetric Lanczos factorization with exit code info = ", info, "."
      end if

      return
   end subroutine nonsymmetric_lanczos_tridiagonalization

   !-------------------------------------------------------------------------------
   !-----                                                                     -----
   !-----     TWO-SIDED ARNOLDI FACTORIZATION FOR GENERAL SQUARE MATRICES     -----
   !-----                                                                     -----
   !-------------------------------------------------------------------------------

   subroutine two_sided_arnoldi_factorization(A, V, W, H, G, info, kstart, kend, verbosity, tol)
      !> Linear operator to be factorized.
      class(abstract_linop), intent(in)     :: A
      !> Left and right Krylov basis.
      class(abstract_vector), intent(inout)  :: V(:)
      class(abstract_vector), intent(inout)  :: W(:)
      !> Upper Hessenberg matrices.
      real(kind=wp), intent(inout) :: H(:, :), G(:, :)
      !> Information flag.
      integer, intent(out) :: info
      !> Optional arguments.
      integer, optional, intent(in) :: kstart
      integer                             :: k_start
      integer, optional, intent(in) :: kend
      integer                             :: k_end
      logical, optional, intent(in) :: verbosity
      logical :: verbose
      real(kind=wp), optional, intent(in) :: tol
      real(kind=wp)                       :: tolerance
      !> Miscellaneous.
      real(kind=wp)                       :: alpha, beta, gamma, tmp
      real(kind=wp), allocatable          :: M(:, :), invM(:, :) ! Inner-products matrix and its inverse.
      real(kind=wp)                       :: tmp_vec(size(V) - 1)
      class(abstract_vector), allocatable :: tmp_krylov_vec
      integer                             :: i, j, k, kdim

      !> Check Krylov subspace dimensions.
      if (size(V) .ne. size(W)) then
         write (*, *) "INFO : Left and right Krylov bases have different dimensions."
         info = -1
         return
      else
         kdim = size(V) - 1
      end if

      !> Check the dimensions of the Hessenberg matrices.
      call assert_shape(H, [kdim + 1, kdim], "two_sided_arnoldi", "H")
      call assert_shape(G, [kdim + 1, kdim], "two_sided_arnoldi", "G")

      !> Deals with the optional arguments.
      k_start = optval(kstart, 1)
      k_end = optval(kend, kdim)
      verbose = optval(verbosity, .false.)
      tolerance = optval(tol, rtol)

      !-----------------------------------------------
      !-----     Two-sided Arnoldi iteration     -----
      !-----------------------------------------------

      arnoldi: do k = k_start, k_end
         !> Arnoldi step for the image.
         call arnoldi_factorization(A, V, H, info, k, k, verbose, tolerance, .false.)
         !> Arnoldi step for the domain.
         call arnoldi_factorization(A, W, G, info, k, k, verbose, tolerance, .true.)
      end do arnoldi

      !------------------------------------------------------------
      !-----     Computation of the Rayleigh quotient form     -----
      !------------------------------------------------------------

      !> Inner-product matrix.
      allocate (M(kdim + 1, kdim + 1)); allocate (invM(kdim, kdim))

      do i = 1, size(M, 1)
         do j = 1, size(M, 2)
            M(i, j) = W(i)%dot(V(j))
         end do
      end do

      !> Inverse of the inner-product matrix (in-place).
      invM = M(1:kdim, 1:kdim); call inv(invM)

      !> Update the residual vectors.
      call update_residual_vector(V(kdim + 1), V(1:kdim), W(1:kdim), invM)
      call update_residual_vector(W(kdim + 1), W(1:kdim), V(1:kdim), transpose(invM))

      !> Rayleigh Quotient form of the Hessenberg matrices.
      call rayleigh_quotient_form(H, invM, M)
      call rayleigh_quotient_form(G, transpose(invM), transpose(M))

      return

   contains
      subroutine update_residual_vector(x, V, W, M)
         !> Residual vector.
         class(abstract_vector), intent(inout) :: x
         !> Krylov subspaces.
         class(abstract_vector), intent(in)    :: V(:), W(:)
         !> Inverse of the inner-product matrix.
         real(kind=wp), intent(in)    :: M(:, :)
         !> Temporary arrays.
         class(abstract_vector), allocatable :: dx
         real(kind=wp)                       :: z(size(V))
         !> Low-dimensional vector.
         do i = 1, size(W) - 1
            z(i) = W(i)%dot(x)
         end do
         z = matmul(M, z)
         !> Correction vector.
         call get_vec(dx, V, z); call x%axpby(1.0_wp, dx, -1.0_wp)
         return
      end subroutine update_residual_vector

      subroutine rayleigh_quotient_form(H, invM, M)
         !> Hessenberg matrix.
         real(kind=wp), intent(inout)       :: H(:, :)
         !> (Augmented) inner product matrix and its inverse.
         real(kind=wp), intent(in)          :: invM(:, :), M(:, :)
         !> Miscellaneous.
         integer :: k

         !> Column dimension of the Hessenberg matrix.
         k = size(H, 2)

         !> Rayleigh quotient update.
         H(1:k, :) = matmul(invM, matmul(M(1:k, :), H))

         return
      end subroutine rayleigh_quotient_form
   end subroutine two_sided_arnoldi_factorization

   subroutine qr_factorization(Q, R, info, verbosity, tol)
      !! Orthogonalization of an array of `abstract_vector` using the modified Gram-Schmid process.
      !!
      !! **Algorithmic Features**
      !!
      !! - In-place factorization
      !! - Double Gram-Schmidt procedure for stability
      !! - Includes a simple check for premature breakdown
      !!
      !! **Advantages**
      !!
      !! - Suitable for all Krylov subspace dimensions
      !! - Robust with respect to floating point errors
      !!
      !! **Limitations**
      !!
      !! - No pivoting, i.e. the columns are orthonormalized in the natural ordering.
      !!   This may lead to premature breakdown if adjacent
      !!   columns are nearly colinear.
      !! - The current implementation only includes a simple check for breakdown
      !!
      !! **References**
      !!
      !! - G. H. Golub & C. F. Van Loan. "Matrix Computations". 4th edition, The John Hopkins
      !!   University Press, 2013.
      !!   See Chapter 5.2.8 : Modified Gram-Schmidt algorithm.
      class(abstract_vector), intent(inout) :: Q(:)
      !! Array of `abstract_vector` that need to be orthogonalized.
      real(kind=wp), intent(out)   :: R(:, :)
      !! Upper triangular matrix \(\mathbf{R}\) resulting from the QR factoriation.
      integer, intent(out)   :: info
      !! Information flag.
      logical, optional, intent(in)         :: verbosity
      !! Verbosity control (default `.false.`).
      real(kind=wp), optional, intent(in)   :: tol
      !! Tolerance to determine colinearity (default `sqrt(epsilon(1.0_wp))`).

      ! Internal variables.
      real(kind=wp) :: beta
      logical                               :: verbose
      integer       :: i, j, k, kdim
      real(kind=wp)                         :: tolerance

      info = 0

      ! Get basis size.
      kdim = size(Q)

      ! Deal with the optional arguments.
      verbose = optval(verbosity, .false.)
      tolerance = optval(tol, rtol)

      R = 0.0_wp
      ! Double Gram-Schmidt (To avoid stability issues with the classical GS)
      do j = 1, kdim
         ! Orthonormalization against existing columns
         do i = 1, j - 1
            beta = Q(j)%dot(Q(i)); call Q(j)%axpby(1.0_wp, Q(i), -beta)
            ! Update R
            R(i, j) = beta
         end do
         ! Second pass
         do i = 1, j - 1
            beta = Q(j)%dot(Q(i)); call Q(j)%axpby(1.0_wp, Q(i), -beta)
            ! Update R
            R(i, j) = R(i, j) + beta
         end do
         ! Normalize column
         beta = Q(j)%norm(); call Q(j)%scal(1.0_wp/beta)
         R(j, j) = beta
         if (abs(beta) < tolerance) then
            if (verbose) then
               write (output_unit, *) "INFO : Colinear columns detected."
               write (output_unit, *) "       (j, beta) = (", j, ",", beta, ")."
            end if
         end if
      end do
      if (verbose) then
         write (output_unit, *) "INFO : Exiting the QR factorization with exit code info = ", info, "."
      end if

      return
   end subroutine qr_factorization

   subroutine krylov_schur_restart(nblk, X, H, select_eigs, blksize)
      integer, intent(out) :: nblk
      !! Number of blocks (or vectors if `blksize=1`) that have been moved to the upper left block
      !! of the Schur factorization of `H`.
      class(abstract_vector), intent(inout) :: X(:)
      !! Krylov basis.
      class(abstract_vector), allocatable   :: Xwrk(:)
      real(kind=wp), intent(inout) :: H(:, :)
      !! Hessenberg/Schur matrix of the Krylov decomposition.
      real(kind=wp), allocatable   :: b(:, :)
      interface
         function selector(lambda) result(out)
            import wp
            complex(kind=wp), intent(in) :: lambda(:)
            logical                      :: out(size(lambda))
         end function selector
      end interface
      procedure(selector) :: select_eigs
      !! Routine to select the eigenvalues that need to be moved in the upper left block of the
      !! Schur factorization. Its interface needs to be compliant with
      !!```
      !! interface
      !!    function selector(lambda) result(out)
      !!      import wp
      !!      complex(kind=wp), intent(in) :: lambda(:)
      !!      logical                      :: out(size(lambda))
      !!    end function selector
      !! end interface
      !!```
      integer, optional, intent(in) :: blksize
      !! Block-size if block Arnoldi is being used (default `blksize = 1`).
      !!@warning
      !! The routine will return an error if `blksize /= 1`. This is related to an issue with the current
      !! implement if the number of selected eigenvalues is not a multiple of the block size.
      !!@endwarning
      integer                       :: blksize_

      ! Schur-related.
      real(kind=wp), allocatable :: Z(:, :)
      complex(kind=wp), allocatable :: eigvals(:)
      logical, allocatable :: selected(:)

      ! Internal variables.
      integer :: k, i, kdim

      ! Sets up the optional args.
      blksize_ = optval(blksize, 1)
      kdim = (size(X) - blksize_)/blksize_

      if (blksize_ /= 1) then
         call stop_error("Block Krylov-Schur restart is not supported yet.")
      end if

      ! Allocate variables.
      allocate (eigvals(kdim*blksize_)); eigvals = 0.0_wp
      allocate (Z(kdim*blksize_, kdim*blksize_)); Z = 0.0_wp
      allocate (selected(kdim*blksize_)); selected = .false.
      allocate (b(blksize_, size(H, 2))); b = 0.0_wp

      ! Schur decomposition of the Hessenberg matrix.
      call schur(H(1:kdim*blksize_, :), Z, eigvals)
      ! Eigenvalue selection for the upper left block.
      selected = select_eigs(eigvals); nblk = count(selected) ! Number of selected eigs.
      if (mod(nblk, blksize_) /= 0) then
         nblk = nblk/blksize_ + 1
      else
         nblk = nblk/blksize_
      end if

      ! Re-order Schur decomposition and vectors.
      call ordschur(H(1:kdim*blksize_, 1:kdim*blksize_), Z, selected)

      ! Update Krylov basis.
      allocate (Xwrk, source=X); call mat_mult(X(1:kdim*blksize_), Xwrk(1:kdim*blksize_), Z)
      do i = 1, blksize_
         call X(nblk*blksize_ + i)%axpby(0.0_wp, Xwrk(kdim*blksize_ + i), 1.0_wp)
      end do
      call initialize_krylov_subspace(X((nblk + 1)*blksize_ + 1:size(X)))

      ! Update Hessenberg matrix.
      b = matmul(H(kdim*blksize_ + 1:size(H, 1), :), Z); H(nblk*blksize_ + 1:(nblk + 1)*blksize_, :) = b
      H((nblk + 1)*blksize_ + 1:size(H, 1), :) = 0.0_wp; H(:, nblk*blksize_ + 1:size(H, 2)) = 0.0_wp

      return
   end subroutine krylov_schur_restart

end module lightkrylov_BaseKrylov
