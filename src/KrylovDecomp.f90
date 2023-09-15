module KrylovDecomp
  use KrylovVector
  use LinearOperator
  implicit none

  private
  public :: arnoldi_factorization

contains

  !-----------------------------------
  !-----                         -----
  !-----     POWER ITERATION     -----
  !-----                         -----
  !-----------------------------------

  subroutine power_iteration(A, x, lambda, niter, info, verbosity, tol)
    !> Linear Operator.
    class(abstract_linop), intent(in) :: A
    !> Starting vector.
    class(abstract_vector), intent(inout) :: x
    !> Estimated eigenvalues.
    double precision, intent(out) :: lambda
    !> Maximum number of power iterations.
    integer, intent(in) :: niter
    !> Information.
    integer, intent(out) :: info

    !> Optional arguments.
    logical, optional, intent(in)          :: verbosity
    double precision, optional, intent(in) :: tol

    !> Miscellaneous.
    class(abstract_vector), allocatable :: wrk
    double precision :: alpha, lambda_old, tolerance
    integer          :: k
    logical          :: verbose

    ! --> Deals with the optional arguments.
    if (present(verbosity)) then
       verbose = verbosity
    else
       verbose = .false.
    end if

    if (present(tol)) then
       tolerance = tol
    else
       tolerance = 1e-6
    endif

    ! --> Initialization.
    lambda_old = 0.0D+00 ; lambda = 0.0D+00 ; info = -1

    ! --> Normalize the starting vector.
    alpha = x%norm() ; call x%scalar_mult(1.0D+00 / alpha)

    ! --> Power iteration
    poweriteration: do k = 1, niter
       allocate(wrk, source=x)
       !> Matrix-vector product.
       call A%matvec(wrk, x)
       !> Estimated eigenvalue.
       lambda = wrk%dot(x) ; deallocate(wrk)
       !> Normalize estimated eigenvectors.
       alpha = x%norm() ; call x%scalar_mult(1.0D+00 / alpha)

       !> Check for convergence.
       if (abs(lambda - lambda_old) < tolerance) then
          if (verbose) then
             write(*, *)
             write(*, *) "INFO : Power iteration has converged to lambda = ", lambda, "."
          endif

          info = 0
          ! --> Exit the loop.
          exit poweriteration
       else
          lambda_old = lambda
       endif
    enddo poweriteration

    if (info >= 0) info = 1
    if (verbose) then
       write(*, *) "INFO : Exiting the Power Iteration with exit code info = ", info, "."
    endif

    return
  end subroutine power_iteration

  !---------------------------------------------------------------------
  !-----                                                           -----
  !-----     ARNOLDI FACTORIZATION FOR GENERAL SQUARE MATRICES     -----
  !-----                                                           -----
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
  !-----                                                           -----
  !-----     LANCZOS FACTORIZATION FOR SYM. POS. DEF. MATRICES     -----
  !-----                                                           -----
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------------
  !-----                                                                 -----
  !-----      KRYLOV-SCHUR FACTORIZATION FOR GENERAL SQUARE MATRICES     -----
  !-----                                                                 -----
  !---------------------------------------------------------------------------

end module KrylovDecomp