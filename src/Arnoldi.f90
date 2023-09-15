module KrylovDecomp
  use KrylovVector
  use LinearOperator
  implicit none

  private
  public :: arnoldi_factorization

contains

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  !-----                                                           -----
  !-----     ARNOLDI FACTORIZATION FOR GENERAL SQUARE MATRICES     -----
  !-----                                                           -----
  !---------------------------------------------------------------------
  !---------------------------------------------------------------------

  subroutine arnoldi_factorization(A, X, H, kdim, info, kstart, kend, verbosity, tol)

    ! --> Dimension of the desired Krylov factorization.
    integer, intent(in) :: kdim
    ! --> Optional arguments (mainly for GMRES)
    integer, optional, intent(in) :: kstart, kend
    logical, optional, intent(in) :: verbosity
    double precision, optional, intent(in) :: tol

    integer :: k_start, k_end
    logical :: verbose
    double precision :: tolerance

    ! --> Linear Operator to be factorized.
    class(abstract_linop), intent(in) :: A
    ! --> Krylov basis.
    class(abstract_vector), dimension(kdim+1), intent(inout) :: X
    ! --> Upper Hessenberg matrix.
    double precision, dimension(kdim+1, kdim), intent(inout) :: H
    ! --> Information.
    integer, intent(out) :: info ! info < 0 : The k-step Arnoldi factorization failed.
    ! info = 0 : The k-step Arnoldi factorization succeeded.
    ! info > 0 : An invariant subspace has been computed after k=info steps.
    ! --> Miscellaneous
    double precision :: beta
    integer :: k

    ! --> Deals with the optional arguments.
    if (present(kstart)) then
       k_start = kstart
    else
       k_start = 1
    endif

    if (present(kend)) then
       k_end = kend
    else
       k_end = kdim
    endif

    if (present(verbosity)) then
       verbose = verbosity
    else
       verbose = .false.
    endif

    if (present(tol)) then
       tolerance = tol
    else
       tolerance = 1e-12
    endif

    ! --> Allocate variables.
    call X(k_start+1)%zero()

    ! --> Arnoldi factorization.
    arnoldi: do k = k_start, k_end
       ! --> Matrix-vector product.
       call A%matvec(X(k), X(k+1))
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
          write(*, *)
          write(*, *) "INFO : An invariant subspace has been computed (beta =", beta, ")."

          ! --> Dimension of the computed invariant subspace.
          info = k

          ! --> Exit the Arnoldi iteration.
          exit arnoldi
       else
          ! --> Normalize the new Krylov vector.
          call X(k+1)%scalar_mult(1.0D+00 / beta)
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
    double precision, dimension(:, :), intent(inout) :: H
    class(abstract_vector), dimension(:) :: X
    class(abstract_vector), allocatable :: wrk
    integer :: i
    double precision :: alpha

    ! --> Orthogonalize residual w.r.t to previously computed Krylov vectors.
    do i = 1, k
       ! --> Working array.
       wrk = X(i)
       ! --> Update Hessenberg matrix.
       alpha = X(k+1)%dot(wrk) ; H(i, k) = alpha
       ! --> Orthogonalize residual vector.
       call wrk%scalar_mult(alpha) ; call X(k+1)%sub(wrk)
    enddo

    ! --> Perform full re-orthogonalization (see instability of MGS process)
    do i = 1, k
       ! --> Working array.
       wrk = X(i)
       ! --> Update Hessenberg matrix.
       alpha = X(k+1)%dot(wrk) ; H(i, k) = H(i, k) + alpha
       ! --> Orthogonalize residual vectors.
       call wrk%scalar_mult(alpha) ; call X(k+1)%sub(wrk)
    enddo

    return
  end subroutine update_hessenberg_matrix

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  !-----                                                           -----
  !-----     LANCZOS FACTORIZATION FOR SYM. POS. DEF. MATRICES     -----
  !-----                                                           -----
  !---------------------------------------------------------------------
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  !-----                                                                 -----
  !-----      KRYLOV-SCHUR FACTORIZATION FOR GENERAL SQUARE MATRICES     -----
  !-----                                                                 -----
  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------

end module KrylovDecomp
