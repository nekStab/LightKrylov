module IterativeSolvers
  use KrylovVector
  use LinearOperator
  use KrylovDecomp
  use stdlib_sorting, only: sort_index, int_size
  implicit none

  private
  public :: eigs, eighs

contains

  !------------------------------------------
  !-----                                -----
  !-----     EIGENVALUE COMPUTATION     -----
  !-----                                -----
  !------------------------------------------

  subroutine eigs(A, X, eigvecs, eigvals, residuals, info, verbosity)
    !> Linear Operator.
    class(abstract_linop), intent(in) :: A
    !> Krylov basis.
    class(abstract_vector), dimension(:), intent(inout) :: X
    !> Coordinates of eigenvectors in Krylov basis, eigenvalues and associated residuals.
    double complex  , dimension(size(X)-1, size(X)-1), intent(out) :: eigvecs
    double complex  , dimension(size(X)-1)           , intent(out) :: eigvals
    double precision, dimension(size(X)-1)           , intent(out) :: residuals
    !> Information flag.
    integer, intent(out) :: info
    !> Verbosity control.
    logical, optional, intent(in) :: verbosity
    logical :: verbose

    !> Upper Hessenberg matrix.
    double precision, dimension(size(X), size(X)-1) :: H
    double precision                                :: beta
    !> Krylov subspace dimension.
    integer :: kdim
    !> Miscellaneous.
    integer :: i, j, k
    integer(int_size), dimension(size(X)-1) :: indices
    double precision, dimension(size(X)-1) :: abs_eigvals

    ! --> Dimension of the Krylov subspace.
    kdim = size(X) - 1

    ! --> Deals with the optional arguments.
    if (present(verbosity)) then
       verbose = verbosity
    else
       verbose = .false.
    endif

    ! --> Initialize variables.
    H = 0.0D+00 ; residuals = 0.0D+00 ; eigvals = (0.0D+00, 0.0D+00) ; eigvecs = (0.0D+00, 0.0D+00)
    do i = 2, size(X) ! i=1 is the initial Krylov vector given by the user.
       call X(i)%zero()
    enddo

    ! --> Compute Arnoldi factorization.
    call arnoldi_factorization(A, X, H, info)

    if (info < 0) then
       if (verbose) then
          write(*, *) "INFO : Arnoldi iteration failed. Exiting the eig subroutine."
          write(*, *) "       Arnoldi exit code :", info
       endif
       info = -1
       return
    endif

    ! --> Compute spectral decomposition of the Hessenberg matrix.
    call evd(H(1:kdim, 1:kdim), eigvecs, eigvals, kdim)

    ! --> Sort eigenvalues with decreasing magnitude.
    abs_eigvals = abs(eigvals) ; call sort_index(abs_eigvals, indices, reverse=.true.)
    eigvals(:) = eigvals(indices) ; eigvecs = eigvecs(:, indices)

    ! --> Compute the residual associated with each eigenpair.
    beta = H(kdim+1, kdim) !> Get Krylov residual vector norm.
    do i = 1, kdim
       residuals(i) = abs(beta * eigvecs(kdim, i))
    enddo

    return
  end subroutine eigs

  subroutine eighs(A, X, eigvecs, eigvals, residuals, info, verbosity)
    !> Linear Operator.
    class(abstract_spd_linop), intent(in) :: A
    !> Krylov basis.
    class(abstract_vector), dimension(:), intent(inout) :: X
    !> Coordinates of eigenvectors in Krylov basis, eigenvalues, and associated residuals
    double precision, dimension(size(X)-1, size(X)-1), intent(out) :: eigvecs
    double precision, dimension(size(X)-1)           , intent(out) :: eigvals
    double precision, dimension(size(X)-1)           , intent(out) :: residuals
    !> Information flag.
    integer, intent(out) :: info
    !> Verbosity control.
    logical, optional, intent(in) :: verbosity
    logical :: verbose

    !> Tridiagonal matrix.
    double precision, dimension(size(X), size(X)-1) :: T
    double precision                                :: beta
    !> Krylov subspace dimension.
    integer :: kdim
    !> Miscellaneous.
    integer :: i, j, k
    integer(int_size), dimension(size(X)-1) :: indices

    ! --> Dimension of the Krylov subspace.
    kdim = size(X) - 1

    ! --> Deals with the optional argumentlambdas.
    if (present(verbosity)) then
       verbose = verbosity
    else
       verbose = .false.
    endif

    ! --> Initialize all variables.
    T = 0.0D+00 ; residuals = 0.0D+00 ; eigvecs = 0.0D+00 ; eigvals = 0.0D+00
    do i = 2, size(X) ! i = 1 is the starting Krylov vector provided by the user.
       call X(i)%zero()
    enddo

    ! --> Compute Lanczos factorization.
    call lanczos_factorization(A, X, T, info)

    if (info < 0) then
       if (verbose) then
          write(*, *) "INFO : Lanczos iteration failed. Exiting the eigh subroutine."
          write(*, *) "       Lanczos exit code :", info
       endif
       info = -1
       return
    endif

    ! --> Compute spectral decomposition of the tridiagonal matrix.
    call hevd(T(1:kdim, 1:kdim), eigvecs, eigvals, kdim)

    ! --> Sort eigenvalues in decreasing order.
    call sort_index(eigvals, indices, reverse=.true.) ; eigvecs = eigvecs(:, indices)

    ! --> Compute the residual associated with each eigenpair.
    beta = T(kdim+1, kdim) !> Get Krylov residual vector norm.
    do i = 1, kdim
       residuals(i) = abs(beta * eigvecs(kdim, i))
    enddo

    return
  end subroutine eighs

  subroutine evd(A, vecs, vals, n)
    !> Lapack job.
    character*1 :: jobvl = "N", jobvr = "V"
    integer :: n, lwork, info, lda, ldvl, ldvr
    double precision, dimension(n, n) :: A, A_tilde, vr
    double precision, dimension(1, n) :: vl
    double precision, dimension(4*n)  :: work
    double precision, dimension(n)    :: wr, wi
    double complex, dimension(n, n)   :: vecs
    double complex, dimension(n)      :: vals
    integer :: i
    integer, dimension(n) :: idx

    interface
       pure subroutine dgeev(fjobvl, fjobvr, fn, fa, flda, fwr, fwi, fvl, fldvl, fvr, fldvr, fwork, flwork, finfo)
         character, intent(in) :: fjobvl, fjobvr
         integer  , intent(in) :: fn, flda, fldvl, fldvr, flwork, finfo
         double precision, intent(inout) :: fa(flda, *)
         double precision, intent(out) :: fwr(fn), fwi(fn), fvl(fldvl, *), fvr(fldvr, *), fwork(flwork)
       end subroutine dgeev
    end interface

    ! --> Compute the eigendecomposition of A.
    lda = n ; ldvl = 1 ; ldvr = n ; lwork = 4*n ; A_tilde = A
    call dgeev(jobvl, jobvr, n, A_tilde, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)

    ! --> Real to complex arithmetic.
    !     NOTE : Check if a LAPACK function already exists for that purpose.
    vals = wr*(1.0D+00, 0.0D+00) + wi*(0.0D+00, 1.0D+00)
    vecs = vr*(1.0D+00, 0.0D+00)

    do i = 1, n-1
       if (wi(i) .gt. 0) then
          vecs(:, i) = vr(:, i)*(1.0D+00, 0.0D+00) + vr(:, i+1)*(0.0D+00, 1.0D+00)
          vecs(:, i+1) = conjg(vecs(:, i))
       else if (abs(wi(i)) .le. epsilon(wi(i))) then
          vecs(:, i) = vr(:, i) * (1.0D+00, 0.0D+00)
       endif
    enddo

    return
  end subroutine evd

  subroutine hevd(A, vecs, vals, n)
    !> Lapack job.
    character :: jobz="V", uplo="U"
    integer :: n, lwork, info, lda
    double precision, dimension(n) :: vals
    double precision, dimension(n, n) :: A, A_tilde, vecs
    double precision, dimension(3*n-1) :: work
    integer :: i, j, k

    interface
       pure subroutine dsyev(fjobz, fuplo, fn, fa, flda, fw, fwork, flwork, finfo)
         character, intent(in) :: fjobz, fuplo
         integer, intent(in)  :: fn
         integer, intent(in)  :: flda
         integer, intent(in)  :: flwork
         integer, intent(out) :: finfo
         double precision, intent(inout) :: fa(flda, *)
         double precision, intent(out)   :: fw(*)
         double precision, intent(out)   :: fwork(*)
       end subroutine dsyev
    end interface

    ! --> Compute the eigendecomposition of A.
    lda = n ; lwork = 3*n-1 ; A_tilde = A
    call dsyev(jobz, uplo, n, A_tilde, lda, vals, work, lwork, info)

    return
  end subroutine hevd

  !----------------------------------------------
  !-----                                    -----
  !-----     SINGULAR VALUE COMPUTATION     -----
  !-----                                    -----
  !----------------------------------------------

  ! subroutine svds(A, U, S, V, niter, info, verbosity, tol)
  !   !> Linear Operator.
  !   class(abstract_linop), intent(in) :: A
  !   !> Singular triplets.
  !   class(abstract_vector), dimension(:), allocatable, intent(out) :: U
  !   double precision, dimension(:), allocatable, intent(out)       :: S
  !   class(abstract_vector), dimension(:), allocatable, intent(out) :: V
  !   !> Maximum number of iterations.
  !   integer, intent(in) :: niter
  !   !> Information flag.
  !   integer, intent(out) :: info
  !   !> Verbosity control.
  !   logical, optional, intent(in) :: verbosity
  !   logical :: verbose
  !   !> Tolerance control.
  !   double precision, optional, intent(in) :: tol
  !   double precision :: tolerance

  !   ! --> Deals with the optional arguments.
  !   if (present(verbosity)) then
  !      verbose = verbosity
  !   else
  !      verbose = .false.
  !   endif

  !   if (present(tol)) then
  !      tolerance = tol
  !   else
  !      tolerance = 1e-12
  !   endif

  !   return
  ! end subroutine svds

  !--------------------------------------------
  !-----                                  -----
  !-----     ITERATIVE LINEAR SOLVERS     -----
  !-----                                  -----
  !--------------------------------------------

  ! subroutine gmres(A, b)
  !   !> Linear problem.
  !   class(abstract_linop) , intent(in) :: A ! Linear Operator.
  !   class(abstract_vector), intent(in) :: b ! Right-hand side.
  !   return
  ! end subroutine gmres

  ! subroutine cg(A, b)
  !   !> Linear problem.
  !   class(abstract_spd_linop), intent(in) :: A ! Linear Operator.
  !   class(abstract_vector)   , intent(in) :: b ! Right-hand side.
  !   return
  ! end subroutine cg

end module IterativeSolvers
