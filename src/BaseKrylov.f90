module lightkrylov_BaseKrylov
  use lightkrylov_Utils
  use lightkrylov_AbstractVector
  use lightkrylov_LinearOperator
  use stdlib_optval, only : optval
  implicit none
  include "dtypes.h"

  private
  public :: arnoldi_factorization, &
       lanczos_tridiagonalization, &
       lanczos_bidiagonalization,  &
       nonsymmetric_lanczos_tridiagonalization, &
       two_sided_arnoldi_factorization, &
       qr_factorization, &
       initialize_krylov_subspace

contains

  subroutine initialize_krylov_subspace(X)
    !> Krylov subspace to be zeroed-out.
    class(abstract_vector), intent(inout) :: X(:)
    !> Internal variables.
    integer :: i
    do i = 1, size(X)
       call X(i)%zero()
    enddo
    return
  end subroutine initialize_krylov_subspace

  !---------------------------------------------------------------------
  !-----                                                           -----
  !-----     ARNOLDI FACTORIZATION FOR GENERAL SQUARE MATRICES     -----
  !-----                                                           -----
  !---------------------------------------------------------------------

  subroutine arnoldi_factorization(A, X, H, info, kstart, kend, verbosity, tol, transpose)

    ! --> Optional arguments (mainly for GMRES)
    integer, optional, intent(in) :: kstart, kend
    logical, optional, intent(in) :: verbosity, transpose
    real(kind=wp), optional, intent(in) :: tol

    integer :: k_start, k_end
    logical :: verbose, trans
    real(kind=wp) :: tolerance

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
    integer :: k, kdim

    info = 0

    ! --> Check dimensions.
    kdim = size(X) - 1 ; call assert_shape(H, [kdim+1, kdim], "arnoldi_factorization", "H")

    ! --> Deals with the optional arguments.
    k_start   = optval(kstart, 1)
    k_end     = optval(kend, kdim)
    verbose   = optval(verbosity, .false.)
    tolerance = optval(tol, atol)
    trans     = optval(transpose, .false.)

    ! --> Arnoldi factorization.
    arnoldi: do k = k_start, k_end
       ! --> Matrix-vector product.
       if (trans) then
          call A%rmatvec(X(k), X(k+1))
       else
          call A%matvec(X(k), X(k+1))
       endif
       ! --> Update Hessenberg matrix.
       call update_hessenberg_matrix(H, X, k)
       beta = X(k+1)%norm() ; H(k+1, k) = beta

       if (verbose) then
          write(*, *) "--> Arnoldi iteration n°", k, "/", k_end
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
          info = k

          ! --> Exit the Arnoldi iteration.
          exit arnoldi
       else
          ! --> Normalize the new Krylov vector.
          call X(k+1)%scal(1.0D+00 / beta)
       endif

    enddo arnoldi

    if(verbose) then
       write(*, *) "INFO : Exiting the Arnoldi factorization with exit code info =", info, "."
       write(*, *)
    endif

    return
  end subroutine arnoldi_factorization

  subroutine update_hessenberg_matrix(H, X, k)
    integer, intent(in) :: k
    real(kind=wp), intent(inout) :: H(:, :)
    class(abstract_vector) :: X(:)
    class(abstract_vector), allocatable :: wrk
    integer :: i
    real(kind=wp) :: alpha

    ! --> Orthogonalize residual w.r.t to previously computed Krylov vectors.
    do i = 1, k
       alpha = X(k+1)%dot(X(i)) ; call X(k+1)%axpby(1.0_wp, X(i), -alpha)
       ! --> Update Hessenberg matrix.
       H(i, k) = alpha
    enddo

    ! --> Perform full re-orthogonalization (see instability of MGS process)
    do i = 1, k
       alpha = X(k+1)%dot(X(i)) ; call X(k+1)%axpby(1.0_wp, X(i), -alpha)
       ! --> Update Hessenberg matrix.
       H(i, k) = H(i, k) + alpha
    enddo

    return
  end subroutine update_hessenberg_matrix

  !--------------------------------------------------------------------------
  !-----                                                                -----
  !-----     LANCZOS TRIDIAGONALIZATION FOR SYM. POS. DEF. MATRICES     -----
  !-----                                                                -----
  !--------------------------------------------------------------------------

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
    kdim = size(X) - 1 ; call assert_shape(T, [kdim+1, kdim], "lanczos_tridiag.", "T")

    ! --> Deals with the optional arguments.
    k_start   = optval(kstart, 1)
    k_end     = optval(kend, kdim)
    verbose   = optval(verbosity, .false.)
    tolerance = optval(tol, atol)

    ! --> Lanczos tridiagonalization.
    lanczos: do k = k_start, k_end
       ! --> Matrix-vector product.
       call A%matvec(X(k), X(k+1))
       ! --> Update tridiagonal matrix.
       call update_tridiag_matrix(T, X, k)
       beta = X(k+1)%norm() ; T(k+1, k) = beta

       if (verbose) then
          write(*, *) "--> Lanczos iteration n°", k, "/", k_end
          write(*, *) "    -----------------"
          write(*, *) "    + Residual norm :", beta
          write(*, *) "    + Elapsed time  :"
          write(*, *) "    + ETA           :"
       endif

       ! --> Exit Lanczos loop if needed.
       if (beta < tolerance) then
          if (verbose) then
             write(*, *)
             write(*, *) "INFO : An invariant subspace has been computed (beta =)", beta, ")."
          endif

          ! --> Dimension of the computed invariant subspace.
          info = k

          ! --> Exit the Lanczos iteration.
          exit lanczos
       else
          ! --> Normalize the new Krylov vector.
          call X(k+1)%scal(1.0D+00 / beta)
       endif

    enddo lanczos

    if (verbose) then
       write(*, *) "INFO : Exiting the Lanczos factorization with exit code info =", info, "."
    endif

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
    do i = max(1, k-1), k
       alpha = X(k+1)%dot(X(i)) ; call X(k+1)%axpby(1.0_wp, X(i), -alpha)
       ! --> Update tridiag matrix.
       T(i, k) = alpha
    enddo

    ! --> Full re-orthogonalization.
    do i = 1, k
       alpha = X(k+1)%dot(X(i)) ; call X(k+1)%axpby(1.0_wp, X(i), -alpha)
    enddo

    return
  end subroutine update_tridiag_matrix

  !----------------------------------------------------------------------------
  !-----                                                                  -----
  !-----     LANCZOS BIDIAGONALIZATION FOR SINGULAR VALUE COMPUTATION     -----
  !-----                                                                  -----
  !----------------------------------------------------------------------------

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
    class(abstract_vector), allocatable :: wrk

    ! --> Check Krylov subspaces dimensions.
    if (size(U) .ne. size(V)) then
       write(*, *) "INFO : Left and right Krylov basis have different dimensions."
       info = -1
       return
    else
       kdim = size(U) - 1
    endif

    ! --> Check B dimensions.
    call assert_shape(B, [kdim+1, kdim], "lanczos_bidiag.", "B")

    ! --> Deals with the optional arguments.
    k_start = optval(kstart, 1)
    k_end   = optval(kend, kdim)
    verbose = optval(verbosity, .false.)
    tolerance = optval(tol, atol)

    ! --> Lanczos bidiagonalization.
    lanczos : do k = k_start, k_end
       ! --> Transpose matrix-vector product.
       call A%rmatvec(U(k), V(k))
       ! /!\ NOTE : Next lines not needed because already taken care of in the
       !     full re-orthogonalization.
       ! if (k > 1) then
       !    call V(k)%axpby(1.0_wp, V(k-1), -beta)
       ! endif

       ! --> Full re-orthogonalization of the right Krylov subspace.
       do j = 1, k-1
          gamma = V(k)%dot(V(j)) ; call V(k)%axpby(1.0_wp, V(j), -gamma)
       enddo

       ! --> Normalization step.
       alpha = V(k)%norm() ; B(k, k) = alpha
       if (alpha > tolerance) then
          call V(k)%scal(1.0_wp / alpha)
       else
          if (verbose) then
             write(*, *) "INFO : alpha = ", alpha
          endif
          info = k
          exit lanczos
       endif

       ! --> Matrix-vector product.
       call A%matvec(V(k), U(k+1))
       ! /!\ NOTE : Not needed because taken care of in the full reortho. step.
       ! call U(k+1)%axpby(1.0_wp, U(k), -alpha)

       ! --> Full re-orthogonalization of the left Krylov subspace.
       do j = 1, k
          gamma = U(k+1)%dot(U(j)) ; call U(k+1)%axpby(1.0_wp, U(j), -gamma)
       enddo

       ! --> Normalization step.
       beta = U(k+1)%norm() ; B(k+1, k) = beta
       if (beta > tolerance) then
          call U(k+1)%scal(1.0_wp / beta)
       else
          if (verbose) then
             write(*, *) "INFO : beta = ", beta
          endif
          info = k
          exit lanczos
       endif

    enddo lanczos

    return
  end subroutine lanczos_bidiagonalization

  !------------------------------------------------------------------------------------
  !-----                                                                          -----
  !-----     TWO-SIDED LANCZOS TRIDIAGONALIZATION FOR GENERAL SQUARE MATRICES     -----
  !-----                                                                          -----
  !------------------------------------------------------------------------------------

  subroutine nonsymmetric_lanczos_tridiagonalization(A, V, W, T, info, kstart, kend, verbosity, tol)
    !> Linear operator to be factorized.
    class(abstract_linop) , intent(in)    :: A
    !> Left and right Krylov basis.
    class(abstract_vector), intent(inout) :: V(:)
    class(abstract_vector), intent(inout) :: W(:)
    !> Tridiagonal matrix.
    real(kind=wp)         , intent(inout) :: T(:, :)
    !> Information flag.
    integer               , intent(out)   :: info
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
       write(*, *) "INFO : Left and right Krylov basis have different dimensions."
       info = -1
    else
       kdim = size(V) - 1
    endif

    ! --> Check T dimensions.
    if (all(shape(T) .ne. [kdim+1, kdim+1])) then
       write(*, *) "INFO : Tridiagonal matrix and Krylov subspaces dimensions do not match."
       info = -2
    endif

    ! --> Deals with the optional arguments.
    k_start = optval(kstart, 1)
    k_end   = optval(kend, kdim)
    verbose = optval(verbosity, .false.)
    tolerance = optval(tol, rtol)

    ! --> Bi-orthogonalize the left and right starting vectors.
    if (k_start .eq. 1) then
       tmp = V(1)%dot(W(1)) ; beta = sqrt(abs(tmp)) ; gamma = sign(beta, tmp)
       call V(1)%scal(1.0_wp / beta) ; call W(1)%scal(1.0_wp / gamma)
    endif

    ! --> Nonsymmetric Lanczos iterations.
    lanczos : do k = k_start, k_end
       ! --> Matrix-vector product.
       call A%matvec(V(k), V(k+1))
       call A%rmatvec(W(k), W(k+1))

       ! --> Update diagonal entry of the nonsymmetric tridiagonal matrix.
       alpha = V(k+1)%dot(W(k)) ; T(k, k) = alpha

       ! --> Lanczos three term recurrence.
       call V(k+1)%axpby(1.0_wp, V(k), -alpha)
       if (k > 1) call V(k+1)%axpby(1.0_wp, V(k-1), -T(k-1, k))

       call W(k+1)%axpby(1.0_wp, W(k), -alpha)
       if (k > 1) call W(k+1)%axpby(1.0_wp, W(k-1), -T(k, k-1))

       ! --> Full re-biorthogonalization.
       do i = 1, k
          alpha = V(k+1)%dot(W(i)) ; call V(k+1)%axpby(1.0_wp, V(i), -alpha)
          alpha = W(k+1)%dot(V(i)) ; call W(k+1)%axpby(1.0_wp, W(i), -alpha)
       enddo

       ! --> Update the off-diagonal entries of the nonsymmetric tridiagonal matrix.
       tmp = V(k+1)%dot(W(k+1)) ; beta = sqrt(abs(tmp)) ; gamma = sign(beta, tmp)
       T(k, k+1) = gamma ; T(k+1, k) = beta

       if (abs(tmp) < tolerance) then
          if (verbose) then
             write(*, *) "INFO : Invariant subspaces have been computed (beta, gamma) = (", beta, ",", gamma, ")."
          endif
       else
          ! --> Normalization step.
          call V(k+1)%scal(1.0_wp / beta) ; call W(k+1)%scal(1.0_wp / gamma)
       endif

    enddo lanczos

    if (verbose) then
       write(*, *) "INFO : Exiting the nonsymmetric Lanczos factorization with exit code info = ", info, "."
    endif

    return
  end subroutine nonsymmetric_lanczos_tridiagonalization

  !-------------------------------------------------------------------------------
  !-----                                                                     -----
  !-----     TWO-SIDED ARNOLDI FACTORIZATION FOR GENERAL SQUARE MATRICES     -----
  !-----                                                                     -----
  !-------------------------------------------------------------------------------

  subroutine two_sided_arnoldi_factorization(A, V, W, H, G, info, kstart, kend, verbosity, tol)
    !> Linear operator to be factorized.
    class(abstract_linop),  intent(in)     :: A
    !> Left and right Krylov basis.
    class(abstract_vector), intent(inout)  :: V(:)
    class(abstract_vector), intent(inout)  :: W(:)
    !> Upper Hessenberg matrices.
    real(kind=wp),          intent(inout) :: H(:, :), G(:, :)
    !> Information flag.
    integer,                intent(out) :: info
    !> Optional arguments.
    integer, optional,       intent(in) :: kstart
    integer                             :: k_start
    integer, optional,       intent(in) :: kend
    integer                             :: k_end
    logical, optional,       intent(in) :: verbosity
    logical :: verbose
    real(kind=wp), optional, intent(in) :: tol
    real(kind=wp)                       :: tolerance
    !> Miscellaneous.
    real(kind=wp)                       :: alpha, beta, gamma, tmp
    real(kind=wp), allocatable          :: M(:, :), invM(:, :) ! Inner-products matrix and its inverse.
    real(kind=wp)                       :: tmp_vec(size(V)-1)
    class(abstract_vector), allocatable :: tmp_krylov_vec
    integer                             :: i, j, k, kdim

    !> Check Krylov subspace dimensions.
    if (size(V) .ne. size(W)) then
       write(*, *) "INFO : Left and right Krylov bases have different dimensions."
       info = -1
       return
    else
       kdim = size(V) - 1
    endif

    !> Check the dimensions of the Hessenberg matrices.
    call assert_shape(H, [kdim+1, kdim], "two_sided_arnoldi", "H")
    call assert_shape(G, [kdim+1, kdim], "two_sided_arnoldi", "G")


    !> Deals with the optional arguments.
    k_start   = optval(kstart, 1)
    k_end     = optval(kend, kdim)
    verbose   = optval(verbosity, .false.)
    tolerance = optval(tol, rtol)

    !-----------------------------------------------
    !-----     Two-sided Arnoldi iteration     -----
    !-----------------------------------------------

    arnoldi : do k = k_start, k_end
       !> Arnoldi step for the image.
       call arnoldi_factorization(A, V, H, info, k, k, verbose, tolerance, .false.)
       !> Arnoldi step for the domain.
       call arnoldi_factorization(A, W, G, info, k, k, verbose, tolerance, .true.)
    enddo arnoldi

    !------------------------------------------------------------
    !-----     Computation of the Rayleigh quotient form     -----
    !------------------------------------------------------------

    !> Inner-product matrix.
    allocate(M(kdim+1, kdim+1)) ; allocate(invM(kdim, kdim))

    do i = 1, size(M, 1)
       do j = 1, size(M, 2)
          M(i, j) = W(i)%dot(V(j))
       enddo
    enddo

    !> Inverse of the inner-product matrix (in-place).
    invM = M(1:kdim, 1:kdim) ; call inv(invM)

    !> Update the residual vectors.
    call update_residual_vector(V(kdim+1), V(1:kdim), W(1:kdim), invM)
    call update_residual_vector(W(kdim+1), W(1:kdim), V(1:kdim), transpose(invM))

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
      real(kind=wp),          intent(in)    :: M(:, :)
      !> Temporary arrays.
      class(abstract_vector), allocatable :: dx
      real(kind=wp)                       :: z(size(V))
      !> Low-dimensional vector.
      do i = 1, size(W)-1
         z(i) = W(i)%dot(x)
      enddo
      z = matmul(M, z)
      !> Correction vector.
      call get_vec(dx, V, z) ; call x%axpby(1.0_wp, dx, -1.0_wp)
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

  !-------------------------------------------------------------
  !-----                                                   -----
  !-----     QR FACTORIZATION FOR GENERAL NxM MATRICES     -----
  !-----                                                   -----
  !-------------------------------------------------------------
  !
  ! Purpose:
  ! --------
  ! Simple implementation of the QR factorisation of a general (real)
  ! Krylov basis.
  !
  ! Algorithmic Features:
  ! ---------------------
  ! - In-place factorisation
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
  ! - No pivoting, i.e. the columns are orthonormalized in the natural
  !   ordering. This may lead to premature breakdown if adjacent 
  !   columns are nearly colinear.
  ! - The current implementation only includes a simple check for 
  !   breakdown
  !
  ! Input/Output Parameters:
  ! ------------------------
  ! - Q        : Krylov basis to orthonormalize       [Input/Output]
  ! - R        : Upper triangular coefficient matrix  [Output]
  ! - info     : Information flag                     [Output]
  ! - tol      : Tolerance for breakdown detection    [Optional, Input]
  ! - verbosity: Verbosity control flag               [Optional, Input]
  !  
  subroutine qr_factorization(Q, R, info, verbosity, tol)
   !> Basis to be orthonormalized.
   class(abstract_vector), intent(inout) :: Q(:)
   !> Gram-Schmidt factors
   real(kind=wp)         , intent(out)   :: R(:,:)
   !> Information flag.
   integer               , intent(out)   :: info
   !> Optional arguments.
   logical, optional, intent(in)         :: verbosity
   logical                               :: verbose
   real(kind=wp), optional, intent(in)   :: tol
   real(kind=wp)                         :: tolerance 
   !> Internal variables.
   real(kind=wp) :: beta
   integer       :: i, j, k, kdim

   info = 0

   ! --> Get basis size.
   kdim = size(Q)
    
   ! --> Deal with the optional arguments.
   verbose   = optval(verbosity, .false.)
   tolerance = optval(tol, rtol)

   !> Double Gram-Schmidt (To avoid stability issues with the classical GS)
   do j = 1, kdim
      !> Orthonormalization against existing columns
      do i = 1, j-1
         beta = Q(j)%dot(Q(i)); call Q(j)%axpby(1.0_wp, Q(i), -beta)
         !> Update R
         R(i, j) = beta
      enddo
      beta = Q(j)%norm(); call Q(j)%scal(1.0_wp / beta)
      R(j,j) = beta
      if (abs(beta) < tolerance) then
         if (verbose) then
            write(*, *) "INFO : Colinear columns detected."
            write(*, *) "       (j, beta) = (", j,",", beta, ")."
         endif
      endif
   enddo
   !> second pass, update Q & R
   do j = 1, kdim
      !> Orthonormalization against existing columns
      do i = 1, j-1
         beta = Q(j)%dot(Q(i)); call Q(j)%axpby(1.0_wp, Q(i), -beta)
         !> Update R
         R(i, j) = R(i, j) + beta
      enddo
      beta = Q(j)%norm(); call Q(j)%scal(1.0_wp / beta)
      if (abs(beta) < tolerance) then
         if (verbose) then
            write(*, *) "INFO : Colinear columns detected."
            write(*, *) "       (j, beta) = (", j,",", beta, ")."
         endif
      endif
   enddo

   !> Modified Gram-Schmidt
   !do i = 1, kdim
   !   beta = Q(i)%norm(); call Q(i)%scal(1.0_wp / beta)
   !   R(i, i) = beta
   !   do j = i, kdim
   !      beta = Q(i)%dot(Q(j)); call Q(j)%axpby(1.0_wp, Q(i), -beta)
   !      R(i,j) = beta
   !   enddo
   !enddo

   if (verbose) then
      write(*, *) "INFO : Exiting the QR factorization with exit code info = ", info, "."
   endif

   return
   end subroutine qr_factorization

   subroutine block_arnoldi_factorization(A, X, H, p, info, kstart, kend, verbosity, tol, transpose)

      ! --> Optional arguments (mainly for GMRES)
      integer, optional, intent(in) :: kstart, kend
      logical, optional, intent(in) :: verbosity, transpose
      real(kind=wp), optional, intent(in) :: tol
  
      integer :: k_start, k_end
      logical :: verbose, trans
      real(kind=wp) :: tolerance
  
      ! --> Linear Operator to be factorized.
      class(abstract_linop), intent(in) :: A
      ! --> Krylov basis.
      class(abstract_vector), intent(inout) :: X(:)
      ! --> Upper Hessenberg matrix.
      real(kind=wp), intent(inout) :: H(:, :)
      ! --> Size of the blocks
      integer, intent(in) :: p
      ! --> Information.
      integer, intent(out) :: info ! info < 0 : The k-step Arnoldi factorization failed.
      ! info = 0 : The k-step Arnoldi factorization succeeded.
      ! info > 0 : An invariant subspace has been computed after k=info steps.
      ! --> Miscellaneous
      real(kind=wp) :: beta
      integer :: k, kk, i, kdim
      integer :: kpm, kp, kpp
      
      kpm = (k-1)*p
      kp  = kpm + p
      kpp = kp  + p
  
      info = 0
  
      ! --> Check dimensions.
      kdim = size(X) - 1 ; call assert_shape(H, [p*(kdim+1), p*kdim], "arnoldi_factorization", "H")
  
      ! --> Deals with the optional arguments.
      k_start   = optval(kstart, 1)
      k_end     = optval(kend, kdim)
      verbose   = optval(verbosity, .false.)
      tolerance = optval(tol, atol)
      trans     = optval(transpose, .false.)
  
      ! --> Arnoldi factorization.
      block_arnoldi: do k = k_start, k_end
         ! --> Matrix-vector products.
         !> This is a quick-fix that assumes a block arnoldi decomposition with a fixed linear operator.
         !> It would be more general if the user could choose an arbitrary, generic series of rational 
         !> expressions involving A to build the Krylov subspace (like in KPIK alternating between A and Ainv) 
         if (trans) then
            do i = 1,p
               call A%rmatvec(X(kpm+i), X(kp+i))
            enddo
         else
            do i = 1,p
               call A%matvec(X(kpm+i), X(kp+i))
            enddo
         endif
         ! --> Update Hessenberg matrix.
         call block_update_hessenberg_matrix(H, X, k, p)
         call qr_factorization(X(kp+1:kpp),H(kp+1:kpp,kpm+1:kp),info)
  
         !if (verbose) then
         !   write(*, *) "--> Arnoldi iteration n°", k, "/", k_end
         !   write(*, *) "    -----------------"
         !   write(*, *) "    + Residual norm :", beta
         !   write(*, *) "    + Elapsed time  :"
         !   write(*, *) "    + ETA           :"
         !   write(*, *)
         !endif
  
         ! --> Exit Arnoldi loop if needed.
         !if (beta < tolerance) then
         !   if (verbose) then
         !      write(*, *)
         !      write(*, *) "INFO : An invariant subspace has been computed (beta =", beta, ")."
         !   endif
  !
         !   ! --> Dimension of the computed invariant subspace.
         !   info = k
  !
         !   ! --> Exit the Arnoldi iteration.
         !   exit arnoldi
         !else
         !   ! --> Normalize the new Krylov vector.
         !   call X(k+1)%scal(1.0D+00 / beta)
         !endif
  
      enddo block_arnoldi
  
      if(verbose) then
         write(*, *) "INFO : Exiting the block Arnoldi factorization with exit code info =", info, "."
         write(*, *)
      endif
  
      return
    end subroutine block_arnoldi_factorization

    subroutine block_update_hessenberg_matrix(H, X, k, p)
      integer, intent(in) :: k
      ! --> Size of the blocks
      integer, intent(in) :: p
      real(kind=wp), intent(inout) :: H(:, :)
      class(abstract_vector) :: X(:)
      real(wp), allocatable :: wrk(:,:)
      integer :: i, j, kpm, kp, kpp
      real(kind=wp) :: alpha
  
      kpm = (k-1)*p
      kp  = kpm + p
      kpp = kp  + p
      allocate(wrk(1:kp,1:p))
      wrk = 0.0_wp
      ! --> Orthogonalize residual w.r.t to previously computed Krylov vectors.
      call mat_mult(H(1:kp,kpm+1:kp), X(1:kp), X(kp+1:kpp))
        
      ! --> Perform full re-orthogonalization (see instability of MGS process)
      call mat_mult(wrk,              X(1:kp), X(kp+1:kpp))

      ! --> Update Hessenberg matrix
      do i = 1,kp
         do j = 1,p
            H(i,kpm+j) = H(i,kpm+j) + wrk(i,kpm+j)
         enddo
      enddo
  
      return
    end subroutine block_update_hessenberg_matrix

end module lightkrylov_BaseKrylov
