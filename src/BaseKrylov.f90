module lightkrylov_BaseKrylov
   !> LightKrylov modules.
   use lightkrylov_Utils
   use lightkrylov_AbstractVector
   use lightkrylov_LinearOperator

   !> Fortran standard library.
   use stdlib_optval, only: optval
   use stdlib_linalg, only: eye

   implicit none
   include "dtypes.h"

   private
   !> Factorization for general n x n linear operators.
   public :: arnoldi_factorization
   !> Factorization for symmetric positive definite operators.
   public :: lanczos_tridiagonalization
   !> Pre-factorization for singular value decomposition.
   public :: lanczos_bidiagonalization
   !> Experimental/Deprecated.
   public :: nonsymmetric_lanczos_tridiagonalization
   public :: two_sided_arnoldi_factorization
   !> Miscellaneous.
   public :: qr_factorization, initialize_krylov_subspace

contains

subroutine initialize_krylov_subspace(X,B0)
   !> Krylov subspace to be initialized
   class(abstract_vector), intent(inout) :: X(:)
   !> Optional: Initial vector/matrix  ::  default is zeros
   class(abstract_vector), optional, intent(in) :: B0(:)
   !> Internal variables.
   class(abstract_vector), allocatable :: B(:)
   real(kind=wp),          allocatable :: Rwrk(:,:)
   real(kind=wp),          allocatable :: Pwrk(:,:)
   integer :: i, p, info

   

   !> zero out X
   call mat_zero(X)

   !> Deals with the optional starting vector
   if (present(B0)) then
      p = size(B0)
      !> Sanity check
      if (size(X) .lt. p) then
         write(*,*) "ERROR : Mismatch between basis size and size of initial vector."
         STOP 1
      endif
      !> allocate & initialize
      allocate(B(1:p), source=B0(1:p))
      call mat_zero(B); call mat_copy(B,B0)
      !> orthonormalize
      allocate(Rwrk(1:p,1:p)); Rwrk = 0.0_wp
      allocate(Pwrk(1:p,1:p)); Pwrk = 0.0_wp
      call qr_factorization(B, Rwrk, Pwrk, info)
      !> Set initial vector
      call mat_copy(X(1:p),B)
   endif      

   return
end subroutine initialize_krylov_subspace

   !=======================================================================================
   ! Arnoldi Factorization Subroutine
   !=======================================================================================
   !
   ! Purpose:
   ! --------
   ! Performs Arnoldi factorization to generate an orthonormal Krylov basis X and an upper
   ! Hessenberg matrix H.
   !
   ! Mathematical Formulation:
   ! -------------------------
   ! Given a linear operator A \in R^(n x n), find X and H such that:
   ! A * X(:, k) = X(:, 1:k) * H(1:k, k) + h_{k+1,k} * X(:, k+1)
   !
   ! Algorithmic Features:
   ! ---------------------
   ! - Constructs an orthonormal Krylov basis X via modified Gram-Schmidt.
   ! - Constructs an upper Hessenberg matrix H.
   ! - Checks for convergence and invariant subspaces.
   !
   ! Advantages:
   ! -----------
   ! - Applicable for non-symmetric matrices.
   ! - Basis for many Krylov subspace methods.
   !
   ! Limitations:
   ! ------------
   ! - Orthogonality of X may deteriorate for ill-conditioned matrices.
   ! - Not suitable for preconditioning in this implementation.
   !
   ! Input/Output Parameters:
   ! ------------------------
   ! A          : Linear Operator, class(abstract_linop), intent(in)
   ! X          : Krylov Basis, class(abstract_vector), dimension(:), intent(inout)
   ! H          : Hessenberg Matrix, double precision, dimension(:, :), intent(inout)
   ! info       : Exit information flag, integer, intent(out)
   ! kstart     : Optional, starting Krylov vector index, integer, intent(in)
   ! kend       : Optional, ending Krylov vector index, integer, intent(in)
   ! verbosity  : Optional, verbosity flag, logical, intent(in)
   ! tol        : Optional, orthogonality tolerance, double precision, intent(in)
   ! transpose  : Optional, transpose flag for A, logical, intent(in)
   !
   ! References:
   ! -----------
   ! - Y. Saad. "Iterative methods for sparse linear systems", SIAM 2nd edition, 2003.
   !   see Chapter 6.3 : Arnoldi's method.
   !
   !=======================================================================================
   subroutine arnoldi_factorization(A, X, H, info, kstart, kend, verbosity, tol, transpose, block_size)

      ! --> Optional arguments (mainly for GMRES)
      integer, optional, intent(in) :: kstart, kend
      logical, optional, intent(in) :: verbosity, transpose
      real(kind=wp), optional, intent(in) :: tol
  
      integer :: k_start, k_end, p
      logical :: verbose, trans
      real(kind=wp) :: tolerance

      ! --> Optional: size of blocks, default = 1
      integer, optional, intent(in) :: block_size
  
      ! --> Linear Operator to be factorized.
      class(abstract_linop), intent(in) :: A
      ! --> Krylov basis.
      class(abstract_vector), intent(inout) :: X(:)
      ! --> Upper Hessenberg matrix.
      real(kind=wp), intent(inout) :: H(:, :)
      ! --> Information.
      integer, intent(out) :: info ! info < 0 : The k-step Arnoldi factorization failed.
      ! info = 0 : The k-step Arnoldi factorization succeeded.
      ! info > 0 : An invariant subspace has been computed after k=info steps.
      ! --> Miscellaneous
      real(kind=wp) :: beta
      real(kind=wp), allocatable :: res(:)
      real(kind=wp), allocatable :: perm(:,:)
      integer :: k, i, kdim, kpm, kp, kpp
      
      ! --> Deals with optional non-unity block size
      p   = optval(block_size, 1)
      allocate(res(1:p))

      allocate(perm(1:size(H,1),1:size(H,2))); perm = 0.0_wp
  
      info = 0

      ! --> Check dimensions.
      kdim = (size(X) - p)/p ; call assert_shape(H, [p*(kdim+1), p*kdim], "arnoldi_factorization", "H")
  
      ! --> Deals with the optional arguments.
      k_start   = optval(kstart, 1)
      k_end     = optval(kend, kdim)
      verbose   = optval(verbosity, .false.)
      tolerance = optval(tol, atol)
      trans     = optval(transpose, .false.)
  
      ! --> Arnoldi factorization.
      block_arnoldi: do k = k_start, k_end
         ! --> Counters
         kpm = (k-1)*p
         kp  = kpm + p
         kpp = kp  + p
         ! --> Matrix-vector products.
         if (trans) then
            do i = 1,p
               call A%rmatvec(X(kpm+i), X(kp+i))
            enddo
         else
            do i = 1,p
               call A%matvec(X(kpm+i), X(kp+i))
            enddo
         endif
         ! --> Update Hessenberg matrix w.r.t. existing vectors
         call update_hessenberg_matrix(H, X, k, p)
         ! --> Orthogonalize current vectors
         call qr_factorization(X(kp+1:kpp), H(kp+1:kpp,kpm+1:kp), perm, info)
         ! --> extract residual norm (smallest diagonal element of R matrix)
         res = 0.0_wp
         do i = 1,p
            res(i) = H(kp+i,kpm+i)
         enddo
         beta = minval(res)
         
         if (verbose) then
            if (p.eq.1) then
               write(*, *) "--> Arnoldi iteration n°", k, "/", k_end
            else
               write(*, *) "--> Block Arnoldi iteration n°", k, "/", k_end
               write(*, *) "    n° of vectors in basis:", kp
            endif
            write(*, *) "    -----------------"
            write(*, *) "    + Residual norm :", beta
            write(*, *) "    + Elapsed time  :"
            write(*, *) "    + ETA           :"
            write(*, *)
         endif
  
         ! --> Exit Arnoldi loop if needed.
         if (beta < tolerance) then
            if (verbose) then
               write(*, *)
               write(*, *) "INFO : An invariant subspace has been computed (beta =", beta, ")."
            endif
  
            ! --> Dimension of the computed invariant subspace.
            info = kp
  
            ! --> Exit the Arnoldi iteration.
            exit block_arnoldi
         endif
  
      enddo block_arnoldi
  
      if(verbose) then
         write(*, *) "INFO : Exiting the block Arnoldi factorization with exit code info =", info, "."
         write(*, *)
      endif
  
      return
    end subroutine arnoldi_factorization

    subroutine update_hessenberg_matrix(H, X, k, block_size)
      integer, intent(in) :: k
      real(kind=wp), intent(inout) :: H(:, :)
      class(abstract_vector) :: X(:)
      ! --> Optional: size of blocks, default = 1
      integer, optional, intent(in) :: block_size
      !> Mist
      class(abstract_vector), allocatable :: Xwrk(:)
      real(wp), allocatable :: wrk(:,:)
      integer :: p, kpm, kp, kpp

      ! --> Deals with optional non-unity block size
      p   = optval(block_size, 1)
    
      kpm = (k-1)*p
      kp  = kpm + p
      kpp = kp  + p
      allocate(wrk(1:kp,1:p)) ; wrk = 0.0_wp
      allocate(Xwrk(1:kp), source=X(1:kp)) ; call mat_zero(Xwrk)
      ! --> Orthogonalize residual w.r.t to previously computed Krylov vectors.
      ! --> pass 1
      call mat_mult(H(1:kp,kpm+1:kp), X(1:kp), X(kp+1:kpp))
      call mat_mult(Xwrk(1:p), X(1:kp), H(1:kp,kpm+1:kp))
      ! --> Project out existing vectors
      call mat_axpby(X(kp+1:kpp), 1.0_wp, Xwrk(1:p), -1.0_wp)
      ! --> pass 2
      call mat_mult(wrk,              X(1:kp), X(kp+1:kpp))
      call mat_mult(Xwrk(1:p), X(1:kp), wrk)
      ! --> Project out existing vectors
      call mat_axpby(X(kp+1:kpp), 1.0_wp, Xwrk(1:p), -1.0_wp)

      ! --> Update Hessenberg matrix with data from second pass
      call mat_axpby(H(1:kp,kpm+1:kp), 1.0_wp, wrk, 1.0_wp)
      deallocate(wrk)
  
      return
    end subroutine update_hessenberg_matrix

   !=======================================================================================
   ! Lanczos Tridiagonalization for Symmetric Positive Definite Matrices
   !=======================================================================================
   !
   ! Purpose:
   ! --------
   ! Performs Lanczos tridiagonalization on a given SPD linear operator A, producing an
   ! orthonormal basis X and a tridiagonal matrix T.
   !
   ! Mathematical Formulation:
   ! -------------------------
   ! Given SPD A \in R^(n x n), find X and T such that:
   ! A * X(:, k) = X(:, 1:k) * T(1:k, k) + t_{k+1,k} * X(:, k+1)
   !
   ! Algorithmic Features:
   ! ---------------------
   ! - Constructs an orthonormal Krylov basis X.
   ! - Constructs a tridiagonal matrix T.
   ! - Checks for convergence and invariant subspaces.
   !
   ! Advantages:
   ! -----------
   ! - Efficient for SPD matrices.
   ! - Foundation for eigenvalue and linear system solvers.
   !
   ! Limitations:
   ! ------------
   ! - Limited to SPD matrices.
   ! - Orthogonality may deteriorate for ill-conditioned matrices.
   !
   ! Input/Output Parameters:
   ! ------------------------
   ! A          : SPD Linear Operator, class(abstract_spd_linop), intent(in)
   ! X          : Krylov Basis, class(abstract_vector), dimension(:), intent(inout)
   ! T          : Tridiagonal Matrix, double precision, dimension(:, :), intent(inout)
   ! info       : Exit information flag, integer, intent(out)
   ! kstart     : Optional, starting Krylov vector index, integer, intent(in)
   ! kend       : Optional, ending Krylov vector index, integer, intent(in)
   ! verbosity  : Optional, verbosity flag, logical, intent(in)
   ! tol        : Optional, orthogonality tolerance, double precision, intent(in)
   !
   ! References:
   ! -----------
   ! - Y. Saad. "Iterative methods for sparse linear systems", SIAM 2nd edition, 2003.
   !   see Chapter 6.6 : The symmetric Lanczos algorithm.
   !
   !=======================================================================================
   subroutine lanczos_tridiagonalization(A, X, T, info, kstart, kend, verbosity, tol)
      !> Linear operator to be factorized;
      class(abstract_spd_linop), intent(in) :: A
      !> Krylov basis.
      class(abstract_vector), intent(inout) :: X(:)
      !> Tri-diagonal matrix.
      real(kind=wp), intent(inout) :: T(:, :)
      !> Information flag.
      integer, intent(out) :: info
      !> Optional arguements.
      integer, optional, intent(in) :: kstart
      integer                       :: k_start
      integer, optional, intent(in) :: kend
      integer                       :: k_end
      logical, optional, intent(in) :: verbosity
      logical                       :: verbose
      real(kind=wp), optional, intent(in) :: tol
      real(kind=wp)                       :: tolerance
      !> Miscellaneous.
      real(kind=wp) :: beta
      integer          :: k, kdim
      integer          :: i, j

      ! --> Check dimensions.
      kdim = size(X) - 1; call assert_shape(T, [kdim + 1, kdim], "lanczos_tridiag.", "T")

      ! --> Deals with the optional arguments.
      k_start = optval(kstart, 1)
      k_end = optval(kend, kdim)
      verbose = optval(verbosity, .false.)
      tolerance = optval(tol, atol)

      ! --> Lanczos tridiagonalization.
      lanczos: do k = k_start, k_end
         ! --> Matrix-vector product.
         call A%matvec(X(k), X(k + 1))
         ! --> Update tridiagonal matrix.
         call update_tridiag_matrix(T, X, k)
         beta = X(k + 1)%norm(); T(k + 1, k) = beta

         if (verbose) then
            write (*, *) "--> Lanczos iteration n°", k, "/", k_end
            write (*, *) "    -----------------"
            write (*, *) "    + Residual norm :", beta
            write (*, *) "    + Elapsed time  :"
            write (*, *) "    + ETA           :"
         end if

         ! --> Exit Lanczos loop if needed.
         if (beta < tolerance) then
            if (verbose) then
               write (*, *)
               write (*, *) "INFO : An invariant subspace has been computed (beta =)", beta, ")."
            end if

            ! --> Dimension of the computed invariant subspace.
            info = k

            ! --> Exit the Lanczos iteration.
            exit lanczos
         else
            ! --> Normalize the new Krylov vector.
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

      ! --> Orthogonalize residual w.r.t to previously computed Krylov vectors.

      do i = max(1, k - 1), k
         alpha = X(k + 1)%dot(X(i)); call X(k + 1)%axpby(1.0_wp, X(i), -alpha)
         ! --> Update tridiag matrix.
         T(i, k) = alpha
      end do
!      ! --> Full re-orthogonalization.
      do i = 1, k
         alpha = X(k + 1)%dot(X(i)); call X(k + 1)%axpby(1.0_wp, X(i), -alpha)
      end do

      return
   end subroutine update_tridiag_matrix

   !=======================================================================================
   ! Lanczos Bidiagonalization for Singular Value Computation
   !=======================================================================================
   !
   ! Purpose:
   ! --------
   ! Performs Lanczos bidiagonalization on a given linear operator A, producing left and
   ! right Krylov bases U and V, and a bidiagonal matrix B.
   !
   ! Mathematical Formulation:
   ! -------------------------
   ! Given A \in R^(m x n), find U, V, and B such that:
   ! A * V = U * B
   !
   ! Algorithmic Features:
   ! ---------------------
   ! - Constructs orthonormal bases U and V.
   ! - Constructs a bidiagonal matrix B.
   !
   ! Input/Output Parameters:
   ! ------------------------
   ! A          : Linear Operator, class(abstract_linop), intent(in)
   ! U          : Left Krylov Basis, class(abstract_vector), dimension(:), intent(inout)
   ! V          : Right Krylov Basis, class(abstract_vector), dimension(:), intent(inout)
   ! B          : Bidiagonal Matrix, real(kind=wp), dimension(:, :), intent(inout)
   ! info       : Exit information flag, integer, intent(out)
   ! kstart     : Optional, starting index, integer, intent(in)
   ! kend       : Optional, ending index, integer, intent(in)
   ! verbosity  : Optional, verbosity flag, logical, intent(in)
   ! tol        : Optional, tolerance, real(kind=wp), intent(in)
   !
   ! References:
   ! -----------
   ! - R. M. Larsen. "Lanczos bidiagonalization with partial reorthogonalization." Technical Report, 1998.
   !   url : http://sun.stanford.edu/~rmunk/PROPACK/paper.pdf
   !
   !=======================================================================================
   subroutine lanczos_bidiagonalization(A, U, V, B, info, kstart, kend, verbosity, tol)
      !> Linear operator to be factorized.
      class(abstract_linop), intent(in) :: A
      !> Left and right Krylov basis.
      class(abstract_vector), intent(inout) :: U(:)
      class(abstract_vector), intent(inout) :: V(:)
      !> Bi-diagonal matrix.
      real(kind=wp), intent(inout) :: B(:, :)
      !> Information flag.
      integer, intent(out) :: info
      !> Optional arguments.
      integer, optional, intent(in) :: kstart
      integer                       :: k_start
      integer, optional, intent(in) :: kend
      integer                       :: k_end
      logical, optional, intent(in) :: verbosity
      logical                       :: verbose
      real(kind=wp), optional, intent(in) :: tol
      real(kind=wp)                       :: tolerance
      !> Miscellaneous.
      real(kind=wp) :: alpha, beta, gamma
      integer       :: i, j, k, kdim

      ! --> Check Krylov subspaces dimensions.
      if (size(U) .ne. size(V)) then
         write (*, *) "INFO : Left and right Krylov basis have different dimensions."
         info = -1
         return
      else
         kdim = size(U) - 1
      end if

      ! --> Check B dimensions.
      call assert_shape(B, [kdim + 1, kdim], "lanczos_bidiag.", "B")

      ! --> Deals with the optional arguments.
      k_start = optval(kstart, 1)
      k_end = optval(kend, kdim)
      verbose = optval(verbosity, .false.)
      tolerance = optval(tol, atol)

      ! --> Lanczos bidiagonalization.
      lanczos: do k = k_start, k_end
         ! --> Transpose matrix-vector product.
         call A%rmatvec(U(k), V(k))
         ! /!\ NOTE : Next lines not needed because already taken care of in the
         !     full re-orthogonalization.
         ! if (k > 1) then
         !    call V(k)%axpby(1.0_wp, V(k-1), -beta)
         ! endif

         ! --> Full re-orthogonalization of the right Krylov subspace.
         do j = 1, k - 1
            gamma = V(k)%dot(V(j)); call V(k)%axpby(1.0_wp, V(j), -gamma)
         end do

         ! --> Normalization step.
         alpha = V(k)%norm(); B(k, k) = alpha
         if (alpha > tolerance) then
            call V(k)%scal(1.0_wp/alpha)
         else
            if (verbose) then
               write (*, *) "INFO : alpha = ", alpha
            end if
            info = k
            exit lanczos
         end if

         ! --> Matrix-vector product.
         call A%matvec(V(k), U(k + 1))
         ! /!\ NOTE : Not needed because taken care of in the full reortho. step.
         ! call U(k+1)%axpby(1.0_wp, U(k), -alpha)

         ! --> Full re-orthogonalization of the left Krylov subspace.
         do j = 1, k
            gamma = U(k + 1)%dot(U(j)); call U(k + 1)%axpby(1.0_wp, U(j), -gamma)
         end do

         ! --> Normalization step.
         beta = U(k + 1)%norm(); B(k + 1, k) = beta
         if (beta > tolerance) then
            call U(k + 1)%scal(1.0_wp/beta)
         else
            if (verbose) then
               write (*, *) "INFO : beta = ", beta
            end if
            info = k
            exit lanczos
         end if

      end do lanczos

      return
   end subroutine lanczos_bidiagonalization

   !=======================================================================================
   ! Nonsymmetric Lanczos Tridiagonalization for General Square Matrices
   !=======================================================================================
   !
   ! Purpose:
   ! --------
   ! Performs Lanczos tridiagonalization on a given nonsymmetric linear operator A, producing
   ! left and right Krylov bases V and W, and a tridiagonal matrix T.
   !
   ! Mathematical Formulation:
   ! -------------------------
   ! Given A \in R^(n x n), find V, W, and T such that:
   ! A * V = W * T
   !
   ! Algorithmic Features:
   ! ---------------------
   ! - Constructs bi-orthogonal bases V and W for the domain and image of A.
   ! - Constructs a tridiagonal matrix T.
   !
   ! Advantages:
   ! -----------
   ! - Simultaneously compute a vector basis for the domain and image of A.
   ! - Bi-orthogonality of the two bases is a feature of the algorithm.
   !
   ! Limitations:
   ! ------------
   ! - Maintaining the bi-orthogonality for high-dimensional operators is numerically
   !   unstable and quickly suffers from accumulation of floating point erros.
   ! - It is strongly suggested NOT to use it for production code.
   !
   ! Input/Output Parameters:
   ! ------------------------
   ! A          : Linear Operator, class(abstract_linop), intent(in)
   ! V          : Left Krylov Basis, class(abstract_vector), dimension(:), intent(inout)
   ! W          : Right Krylov Basis, class(abstract_vector), dimension(:), intent(inout)
   ! T          : Tridiagonal Matrix, real(kind=wp), dimension(:, :), intent(inout)
   ! info       : Exit information flag, integer, intent(out)
   ! kstart     : Optional, starting index, integer, intent(in)
   ! kend       : Optional, ending index, integer, intent(in)
   ! verbosity  : Optional, verbosity flag, logical, intent(in)
   ! tol        : Optional, tolerance, real(kind=wp), intent(in)
   !
   ! References:
   ! -----------
   ! - Y. Saad. "Iterative methods for sparse linear systems", SIAM 2nd edition, 2003.
   !   see Chapter 7.1 : Lanczos biorthogonalization.
   !
   !=======================================================================================
   subroutine nonsymmetric_lanczos_tridiagonalization(A, V, W, T, info, kstart, kend, verbosity, tol)
      !> Linear operator to be factorized.
      class(abstract_linop), intent(in)    :: A
      !> Left and right Krylov basis.
      class(abstract_vector), intent(inout) :: V(:)
      class(abstract_vector), intent(inout) :: W(:)
      !> Tridiagonal matrix.
      real(kind=wp), intent(inout) :: T(:, :)
      !> Information flag.
      integer, intent(out)   :: info
      !> Optional arguments.
      integer, optional, intent(in) :: kstart
      integer                       :: k_start
      integer, optional, intent(in) :: kend
      integer                       :: k_end
      logical, optional, intent(in) :: verbosity
      logical                       :: verbose
      real(kind=wp), optional, intent(in) :: tol
      real(kind=wp)                       :: tolerance
      !> Miscellaneous.
      real(kind=wp) :: alpha, beta, gamma, tmp
      integer       :: i, j, k, kdim

      ! --> Check Krylov subspaces dimensions.
      if (size(V) .ne. size(W)) then
         write (*, *) "INFO : Left and right Krylov basis have different dimensions."
         info = -1
      else
         kdim = size(V) - 1
      end if

      ! --> Check T dimensions.
      if (all(shape(T) .ne. [kdim + 1, kdim + 1])) then
         write (*, *) "INFO : Tridiagonal matrix and Krylov subspaces dimensions do not match."
         info = -2
      end if

      ! --> Deals with the optional arguments.
      k_start = optval(kstart, 1)
      k_end = optval(kend, kdim)
      verbose = optval(verbosity, .false.)
      tolerance = optval(tol, rtol)

      ! --> Bi-orthogonalize the left and right starting vectors.
      if (k_start .eq. 1) then
         tmp = V(1)%dot(W(1)); beta = sqrt(abs(tmp)); gamma = sign(beta, tmp)
         call V(1)%scal(1.0_wp/beta); call W(1)%scal(1.0_wp/gamma)
      end if

      ! --> Nonsymmetric Lanczos iterations.
      lanczos: do k = k_start, k_end
         ! --> Matrix-vector product.
         call A%matvec(V(k), V(k + 1))
         call A%rmatvec(W(k), W(k + 1))

         ! --> Update diagonal entry of the nonsymmetric tridiagonal matrix.
         alpha = V(k + 1)%dot(W(k)); T(k, k) = alpha

         ! --> Lanczos three term recurrence.
         call V(k + 1)%axpby(1.0_wp, V(k), -alpha)
         if (k > 1) call V(k + 1)%axpby(1.0_wp, V(k - 1), -T(k - 1, k))

         call W(k + 1)%axpby(1.0_wp, W(k), -alpha)
         if (k > 1) call W(k + 1)%axpby(1.0_wp, W(k - 1), -T(k, k - 1))

         ! --> Full re-biorthogonalization.
         do i = 1, k
            alpha = V(k + 1)%dot(W(i)); call V(k + 1)%axpby(1.0_wp, V(i), -alpha)
            alpha = W(k + 1)%dot(V(i)); call W(k + 1)%axpby(1.0_wp, W(i), -alpha)
         end do

         ! --> Update the off-diagonal entries of the nonsymmetric tridiagonal matrix.
         tmp = V(k + 1)%dot(W(k + 1)); beta = sqrt(abs(tmp)); gamma = sign(beta, tmp)
         T(k, k + 1) = gamma; T(k + 1, k) = beta

         if (abs(tmp) < tolerance) then
            if (verbose) then
               write (*, *) "INFO : Invariant subspaces have been computed (beta, gamma) = (", beta, ",", gamma, ")."
            end if
         else
            ! --> Normalization step.
            call V(k + 1)%scal(1.0_wp/beta); call W(k + 1)%scal(1.0_wp/gamma)
         end if

      end do lanczos

      if (verbose) then
         write (*, *) "INFO : Exiting the nonsymmetric Lanczos factorization with exit code info = ", info, "."
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

   !=======================================================================================
   ! QR factorization for general n x m matrices
   !=======================================================================================
   !
   ! Purpose:
   ! --------
   ! Simple implementation of the QR factorization of a general (real) Krylov basis.
   !
   ! Algorithmic Features:
   ! ---------------------
   ! - In-place factorization
   ! - Double Gram-Schmidt procedure for stability
   ! - Includes a simple check for premature breakdown
   !
   ! Advantages:
   ! -----------
   ! - Suitable for all Krylov subspace dimensions
   ! - Robust with respect to floating point errors
   !
   ! Limitations:
   ! ------------
   ! - No pivoting, i.e. the columns are orthonormalized in the natural ordering.
   !   This may lead to premature breakdown if adjacent
   !   columns are nearly colinear.
   ! - The current implementation only includes a simple check for breakdown
   !
   ! Input/Output Parameters:
   ! ------------------------
   ! - Q        : Krylov basis to orthonormalize       [Input/Output]
   ! - R        : Upper triangular coefficient matrix  [Output]
   ! - info     : Information flag                     [Output]
   ! - tol      : Tolerance for breakdown detection    [Optional, Input]
   ! - verbosity: Verbosity control flag               [Optional, Input]
   !
   ! References:
   ! -----------
   ! - G. H. Golub & C. F. Van Loan. "Matrix Computations". 4th edition, The John Hopkins
   !   University Press, 2013.
   !   See Chapter 5.2.8 : Modified Gram-Schmidt algorithm.
   !
   !=======================================================================================
   subroutine qr_factorization(Q, R, P, info, ifpivot, verbosity, tol)
      !> Basis to be orthonormalized.
      class(abstract_vector), intent(inout) :: Q(:)
      !> Gram-Schmidt factors
      real(kind=wp), intent(out)            :: R(:, :)
      !> Permutation matrix
      real(kind=wp), intent(out)            :: P(:,:)
      !> Information flag.
      integer, intent(out)                  :: info
      !> Optional arguments.
      logical, optional, intent(in)         :: ifpivot
      logical                               :: pivot
      logical, optional, intent(in)         :: verbosity
      logical                               :: verbose
      real(kind=wp), optional, intent(in)   :: tol
      real(kind=wp)                         :: tolerance
      !> Internal variables.
      real(kind=wp)                         :: beta
      integer                               :: idx, i, j, kdim, iwrk

      integer                               :: idxv(1)
      integer,                allocatable   :: ord(:)
      real(kind=wp),          allocatable   :: Rii(:)
      class(abstract_vector), allocatable   :: Qwrk(1)

      info = 0

      ! --> Get basis size and allocate temporary variables
      kdim = size(Q)
      allocate(ord(1:kdim))
      do i = 1, kdim
         ord(i) = i
      end do
      allocate(Rii(1:kdim)); Rii = 0.0_wp
      allocate(Qwrk(1), source=Q(1))

      ! --> Deal with the optional arguments.
      verbose   = optval(verbosity, .false.)
      pivot     = optval(ifpivot, .false.)
      tolerance = optval(tol, atol)

      R = 0.0_wp
      if (pivot) then

         do i = 1, kdim
            Rii(i) = Q(i)%dot(Q(i))
         end do

         QR_step: do j = 1, kdim
            idxv = maxloc(Rii)
            idx = idxv(1)
            if ( Rii(idx) < tolerance ) then
               if (verbose) then
                  write(*,*) 'INFO : Numerical rank is', j-1
               endif
               exit QR_step
            endif

            call swap_columns(Q, R, Rii, ord, j, idx)
          
            !> Normalize column
            beta = Q(j)%norm();
            !> Check for breakdown
            if (abs(beta) < tolerance) then
               if (verbose) then
                  write (*, *) "INFO : Colinear columns detected."
                  write (*, *) "       (j, beta) = (", j, ",", beta, ")."
               end if
               R(j,j) = 0.0_wp
               call Q(j)%zero()
            else
               call Q(j)%scal(1.0_wp/beta)
               R(j, j) = beta
            end if

            !> Orthonormalize all columns against new vector (MGS)
            do i = j+1, kdim
               beta = Q(j)%dot(Q(i)); call Q(i)%axpby(1.0_wp, Q(j), -beta)
               !> Update R
               R(j, i) = beta
            end do

            ! Update Rii
            Rii(j) = 0.0_wp
            do i = j+1, kdim
               Rii(i) = Rii(i) - R(j,i)**2
            end do

         end do QR_step

         ! compute permutation matrix
         P = 0.0_wp
         do i = 1,kdim
            P(ord(i),i) = 1.0
         end do

      else ! no pivoting

         P = eye(kdim)
         do j = 1, kdim
            !> Orthonormalization against existing columns
            do i = 1, j - 1
               beta = Q(j)%dot(Q(i)); call Q(j)%axpby(1.0_wp, Q(i), -beta)
               !> Update R
               R(i, j) = beta
            end do

            !!> Second pass
            do i = 1, j - 1                                                   ! comment out the
               beta = Q(j)%dot(Q(i)); call Q(j)%axpby(1.0_wp, Q(i), -beta)    ! second pass to
               !> Update R                                                    ! see the simple QR
               R(i, j) = R(i, j) + beta                                       ! factorization fail 
            end do                                                            ! in the test
                                                                
            !> Normalize column
            beta = Q(j)%norm();
            !> Check for breakdown
            if (abs(beta) < tolerance) then
               if (verbose) then
                  write (*, *) "INFO : Colinear columns detected."
                  write (*, *) "       (j, beta) = (", j, ",", beta, ")."
               end if
               R(j,j) = 0.0_wp
               call Q(j)%zero()
            else
               call Q(j)%scal(1.0_wp/beta)
               R(j, j) = beta
            end if
         end do 

      end if
         
      if (verbose) then
         write (*, *) "INFO : Exiting the QR factorization with exit code info = ", info, "."
      end if

      return
   end subroutine qr_factorization

   subroutine swap_columns(Q, R, Rii, ord, i, j)
      implicit none
      !> Krylov Basis
      class(abstract_vector), intent(inout) :: Q(:)
      !> Gram-Schmidt factors
      real(kind=wp),          intent(inout) :: R(:, :)
      !> Diagonal factors
      real(kind=wp),          intent(inout) :: Rii(:)
      !> column ordering 
      integer,                intent(inout) :: ord(:)
      !> column indeces to be swapped
      integer,                intent(in)    :: i, j

      !> internals
      class(abstract_vector), allocatable :: Qwrk
      real(kind=wp),          allocatable :: Rwrk(:)
      integer                             :: iwrk, m, n

      ! sanity check
      m = size(Q);
      n = min(i,j) - 1
      call assert_shape(R, [m, m], "swap_columns", "R")
      if (i .gt. size(Q)) then
         write (*, *) "swap_columns: Index", i, "is out of range."
         STOP 1
      endif
      if (j .gt. size(Q)) then
         write (*, *) "swap_columns: Index", j, "is out of range."
         STOP 1
      endif

      ! allocate memory
      allocate(Qwrk, source=Q(1)); call Qwrk%zero()
      allocate(Rwrk(1:max(1,n))); Rwrk = 0.0_wp

      ! swap columns 
      call Qwrk%axpby(0.0_wp, Q(j), 1.0_wp)
      call Q(j)%axpby(0.0_wp, Q(i), 1.0_wp)
      call Q(i)%axpby(0.0_wp, Qwrk, 1.0_wp)
      
      Rwrk = 0.0_wp
      Rwrk(1) = Rii(j)
      Rii(j)  = Rii(i)
      Rii(i)  = Rwrk(1)

      iwrk   = ord(j)
      ord(j) = ord(i)
      ord(i) = iwrk

      if (n .gt. 0) then
         Rwrk     = R(1:n,j)
         R(1:n,j) = R(1:n,i)
         R(1:n,i) = Rwrk
      endif

      return
   
   end subroutine swap_columns

end module lightkrylov_BaseKrylov
