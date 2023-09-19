module KrylovDecomp
  use KrylovVector
  use LinearOperator
  use stdlib_optval, only : optval
  implicit none

  private
  public :: power_iteration, &
       arnoldi_factorization, &
       lanczos_factorization

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
    verbose   = optval(verbosity, .false.)
    tolerance = optval(tol, 1.0D-12)

    ! --> Initialization.
    lambda_old = 0.0D+00 ; lambda = 0.0D+00 ; info = -1

    ! --> Normalize the starting vector.
    alpha = x%norm() ; call x%scalar_mult(1.0D+00 / alpha) ; wrk = x

    ! --> Power iteration
    do k = 1, niter
       !> Matrix-vector product.
       call A%matvec(wrk, x)
       !> Estimated eigenvalue.
       lambda = wrk%dot(x)
       !> Normalize estimated eigenvectors.
       alpha = x%norm() ; call x%scalar_mult(1.0D+00 / alpha)

       if (verbose) then
          write(*, *) "--> Power iteration n°", k, "/", niter
          write(*, *) "    ---------------"
          write(*, *) "    + Residual norm :", abs(lambda - lambda_old)
          write(*, *) "    + Elapsed time  :"
          write(*, *) "    + ETA           :"
       endif

       !> Check for convergence.
       if (abs(lambda - lambda_old) < tolerance) then
          if (verbose) then
             write(*, *)
             write(*, *) "INFO : Power iteration has converged to lambda = ", lambda, "."
          endif

          info = 0
          ! --> Exit the loop.
          exit
       else
          !> Update old estimates.
          lambda_old = lambda ; wrk = x
       endif
    enddo

    if (verbose) then
       write(*, *) "INFO : Exiting the Power Iteration with exit code info = ", info, "."
       write(*, *)
    endif

    return
  end subroutine power_iteration

  !---------------------------------------------------------------------
  !-----                                                           -----
  !-----     ARNOLDI FACTORIZATION FOR GENERAL SQUARE MATRICES     -----
  !-----                                                           -----
  !---------------------------------------------------------------------

  subroutine arnoldi_factorization(A, X, H, info, kstart, kend, verbosity, tol)

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
    class(abstract_vector), dimension(:), intent(inout) :: X
    ! --> Upper Hessenberg matrix.
    double precision, dimension(:, :), intent(inout) :: H
    ! --> Information.
    integer, intent(out) :: info ! info < 0 : The k-step Arnoldi factorization failed.
    ! info = 0 : The k-step Arnoldi factorization succeeded.
    ! info > 0 : An invariant subspace has been computed after k=info steps.
    ! --> Miscellaneous
    double precision :: beta
    integer :: k, kdim

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
    tolerance = optval(tol, 1.0D-12)

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

  subroutine lanczos_factorization(A, X, T, info, kstart, kend, verbosity, tol)
    !> Linear operator to be factorized;
    class(abstract_spd_linop), intent(in) :: A
    !> Krylov basis.
    class(abstract_vector), dimension(:), intent(inout) :: X
    !> Tri-diagonal matrix.
    double precision, dimension(:, :), intent(inout) :: T
    !> Information flag.
    integer, intent(out) :: info
    !> Optional arguements.
    integer, optional, intent(in) :: kstart
    integer                       :: k_start
    integer, optional, intent(in) :: kend
    integer                       :: k_end
    logical, optional, intent(in) :: verbosity
    logical                       :: verbose
    double precision, optional, intent(in) :: tol
    double precision                       :: tolerance
    !> Miscellaneous.
    double precision :: beta
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
    tolerance = optval(tol, 1.0D-12)

    ! --> Lanczos factorization.
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
          call X(k+1)%scalar_mult(1.0D+00 / beta)
       endif

    enddo lanczos

    if (verbose) then
       write(*, *) "INFO : Exiting the Lanczos factorization with exit code info =", info, "."
    endif

    return
  end subroutine lanczos_factorization

  subroutine update_tridiag_matrix(T, X, k)
    integer, intent(in) :: k
    double precision, dimension(:, :), intent(inout) :: T
    class(abstract_vector), dimension(:) :: X
    class(abstract_vector), allocatable :: wrk
    integer :: i
    double precision :: alpha

    ! --> Orthogonalize residual w.r.t to previously computed Krylov vectors.
    do i = max(1, k-1), k
       ! --> Working array.
       wrk = X(i)
       ! --> Update tridiag matrix.
       alpha = X(k+1)%dot(wrk) ; T(i, k) = alpha
       call wrk%scalar_mult(alpha) ; call X(k+1)%sub(wrk)
    enddo

    return
  end subroutine update_tridiag_matrix

  !---------------------------------------------------------------------------
  !-----                                                                 -----
  !-----      KRYLOV-SCHUR FACTORIZATION FOR GENERAL SQUARE MATRICES     -----
  !-----                                                                 -----
  !---------------------------------------------------------------------------

end module KrylovDecomp
