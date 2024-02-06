module BaseKrylov
  use Utils
  use AbstractVector
  use LinearOperator
  use stdlib_optval, only : optval
  implicit none
  include "dtypes.h"

  private
  public :: arnoldi_factorization, &
       lanczos_tridiagonalization, &
       lanczos_bidiagonalization,  &
       nonsymmetric_lanczos_tridiagonalization, &
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
    class(abstract_vector), dimension(:), intent(inout) :: X
    ! --> Upper Hessenberg matrix.
    real(kind=wp), dimension(:, :), intent(inout) :: H
    ! --> Information.
    integer, intent(out) :: info ! info < 0 : The k-step Arnoldi factorization failed.
    ! info = 0 : The k-step Arnoldi factorization succeeded.
    ! info > 0 : An invariant subspace has been computed after k=info steps.
    ! --> Miscellaneous
    real(kind=wp) :: beta
    integer :: k, kdim

    info = 0

    ! --> Check dimensions.
    kdim = size(X) - 1

    if ( all(shape(H) .ne. [kdim+1, kdim]) ) then
       write(*, *) "INFO : Hessenberg matrix and Krylov subspace dimensions do not match."
       info = -1
       return
    endif

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
    real(kind=wp), dimension(:, :), intent(inout) :: H
    class(abstract_vector), dimension(:) :: X
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
    class(abstract_vector), dimension(:), intent(inout) :: X
    !> Tri-diagonal matrix.
    real(kind=wp), dimension(:, :), intent(inout) :: T
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
    kdim = size(X) - 1

    if (all(shape(T) .ne. [kdim+1, kdim])) then
       write(*, *) "INFO : Tridiagonal matrix and Krylov subspace dimensions do not match."
       info = -1
       return
    endif

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
    real(kind=wp), dimension(:, :), intent(inout) :: T
    class(abstract_vector), dimension(:) :: X
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
    if (all(shape(B) .ne. [kdim+1, kdim])) then
       write(*, *) "INFO : Bidiagonal matrix and Krylov subspaces dimensions do not match."
       info = -2
    endif

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

  !-------------------------------------------------------------
  !-----                                                   -----
  !-----     QR FACTORIZATION FOR GENERAL NxM MATRICES     -----
  !-----                                                   -----
  !-------------------------------------------------------------

  subroutine qr_factorisation(A,Q,R, info, verbosity, tol)
   !> Basis to be orthonormalized.
   class(abstract_vector), intent(in)    :: A(:)
   !> Orthonormal subspace basis
   class(abstract_vector), intent(out)   :: Q(:)
   !> Gram-Schmidt factors
   real(kind=wp)         , intent(out)   :: R(:,:)
   !> Information flag.
   integer       , intent(out)   :: info
   !> Optional arguments.
   logical, optional, intent(in) :: verbosity
   logical                       :: verbose
   integer, optional, intent(in) :: tol
   integer                       :: tolerance 
   !> Local variables.
   real(kind=wp) :: beta
   integer       :: i, j, k, kdim

   ! --> Get basis size.
   kdim = size(V)
    
   ! --> Deal with the optional arguments.
   verbose   = optval(verbosity, .false.)
   tolerance = optval(tol, atol)
   
   ! --> Normalize first column and add it to Q.
   call copy(Q(1),A(1)); beta = Q(1)%norm()
   call Q(1)%scal(1.0_wp / beta); R(1, 1) = beta

   !> Double Gram-Schmidt (To avoid stability issues with the classical GS)
   do i = 2, k
      call copy(Q(i),A(i))
      ! --> Orthonormalization against existing columns
      do j = 1, i-1
         beta = Q(i)%dot(Q(j)); call Q(i)%axpby(1.0_wp, Q(j), -beta)
         ! --> Update R.
         R(i, j) = beta
      enddo
      beta = Q(i)%norm(); Q(i)%scal(1.0_wp / beta)
      R(i, j) = beta
      if (abs(beta) < tolerance) then
         if (verbose) then
            write(*, *) "INFO : Invariant subspace has been computed."
            write(*, *) "       (i, beta) = (", i,",", beta, ")."
         endif
      endif
   enddo
   !> second pass
   do i = 2, k
      ! --> Orthonormalization against existing columns
      do j = 1, i-1
         beta = Q(i)%dot(Q(j)); call Q(i)%axpby(1.0_wp, Q(j), -beta)
         ! --> Update R.
         R(i, j) = R(i, j) + beta
      enddo
      beta = Q(i)%norm(); Q(i)%scal(1.0_wp / beta)
      R(i, j) = beta
      if (abs(beta) < tolerance) then
         if (verbose) then
            write(*, *) "INFO : Invariant subspace has been computed."
            write(*, *) "       (i, beta) = (", i,",", beta, ")."
         endif
      endif
   enddo

   !> Modified Gram-Schmidt
   !do i = 1, k
   !   call copy(Q(i),A(i))
   !enddo
!
   !do i = 1, k
   !   beta = Q(i)%norm(); call Q(i)%scal(1.0_wp / beta)
   !   R(i, i) = beta
   !   do j = i, k
   !      beta = Q(i)%dot(Q(j)); call Q(i)%axpby(1.0_wp, Q(j), -beta)
   !      R(i,j) = beta
   !   enddo
   !enddo

   if (verbose) then
      write(*, *) "INFO : Exiting the QR factorization with exit code info = ", info, "."
   endif

   return
   end subroutine qr_factorisation

end module BaseKrylov
