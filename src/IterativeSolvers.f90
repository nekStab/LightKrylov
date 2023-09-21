module IterativeSolvers
  use KrylovVector
  use LinearOperator
  use KrylovDecomp
  use stdlib_kinds  , only : dp
  use stdlib_sorting, only : sort_index, int_size
  use stdlib_optval , only : optval
  use stdlib_io_npy , only : save_npy
  implicit none

  private
  public :: eigs, eighs, gmres, save_eigenspectrum

contains

  !-----------------------------
  !-----     UTILITIES     -----
  !-----------------------------

  elemental pure function compute_residual(beta, x) result(residual)
    !> Norm of Krylov residual vector.
    double precision, intent(in) :: beta
    !> Last element of Ritz eigenvector.
    double precision, intent(in) :: x
    !> Residual.
    double precision :: residual

    ! --> Compute residual.
    residual = abs(beta * x)
    return
  end function compute_residual

  subroutine save_eigenspectrum(real_part, imag_part, residuals, filename)
    !> Real and imaginary parts of the eigenvalues.
    real(kind=dp), intent(in) :: real_part(:)
    real(kind=dp), intent(in) :: imag_part(:)
    !> Residual norm computed from the Arnoldi/Lanczos factorization.
    real(kind=dp), intent(in) :: residuals(:)
    !> Name of the output file.
    character(len=*), intent(in) :: filename

    !> Miscellaneous.
    real(kind=dp), dimension(size(real_part), 3) :: data

    ! --> Store the data.
    data(:, 1) = real_part ; data(:, 2) = imag_part ; data(:, 3) = residuals
    ! --> Save the eigenspectrum to disk using npy file format.
    call save_npy(filename, data)

    return
  end subroutine save_eigenspectrum

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
    double precision , dimension(size(X)-1) :: abs_eigvals

    ! --> Dimension of the Krylov subspace.
    kdim = size(X) - 1

    ! --> Deals with the optional arguments.
    verbose = optval(verbosity, .false.)

    ! --> Initialize variables.
    H = 0.0D+00 ; residuals = 0.0D+00 ; eigvals = (0.0D+00, 0.0D+00) ; eigvecs = (0.0D+00, 0.0D+00)
    do i = 2, size(X) ! i=1 is the initial Krylov vector given by the user.
       call X(i)%zero()
    enddo

    ! --> Compute Arnoldi factorization.
    call arnoldi_factorization(A, X, H, info, verbosity=verbose)

    if (info < 0) then
       if (verbose) then
          write(*, *) "INFO : Arnoldi iteration failed. Exiting the eig subroutine."
          write(*, *) "       Arnoldi exit code :", info
       endif
       info = -1
       return
    else
       kdim = info
    endif

    ! --> Compute spectral decomposition of the Hessenberg matrix.
    call evd(H(1:kdim, 1:kdim), eigvecs(1:kdim, 1:kdim), eigvals(1:kdim), kdim)

    ! --> Sort eigenvalues with decreasing magnitude.
    abs_eigvals = abs(eigvals) ; call sort_index(abs_eigvals, indices, reverse=.true.)
    eigvals(:) = eigvals(indices) ; eigvecs = eigvecs(:, indices)

    ! --> Compute the residual associated with each eigenpair.
    beta = H(kdim+1, kdim) !> Get Krylov residual vector norm.
    residuals = compute_residual(beta, abs(eigvecs(kdim, :)))

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

    ! --> Deals with the optional argument.
    verbose = optval(verbosity, .false.)

    ! --> Initialize all variables.
    T = 0.0D+00 ; residuals = 0.0D+00 ; eigvecs = 0.0D+00 ; eigvals = 0.0D+00
    do i = 2, size(X) ! i = 1 is the starting Krylov vector provided by the user.
       call X(i)%zero()
    enddo

    ! --> Compute Lanczos tridiagonalization.
    call lanczos_tridiagonalization(A, X, T, info, verbosity=verbose)

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
    residuals = compute_residual(beta, eigvecs(kdim, :))

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

  subroutine svds(A, U, V, singvecs, singvals, residuals, info, verbosity)
    !> Linear Operator.
    class(abstract_linop), intent(in) :: A
    !> Krylov bases.
    class(abstract_vector), intent(inout) :: U(:) ! Basis for left sing. vectors.
    class(abstract_vector), intent(inout) :: V(:) ! Basis for right sing. vectors.
    !> Coordinates of singular vectors in Krylov bases, singular values, and associated residuals.
    real(kind=dp), intent(out) :: singvecs(size(U)-1, size(U)-1)
    real(kind=dp), intent(out) :: singvals(size(U)-1)
    real(kind=dp), intent(out) :: residuals(size(U)-1)
    !> Information flag.
    integer, intent(out) :: info
    !> Verbosity control.
    logical, optional, intent(in) :: verbosity
    logical verbose

    !> Tridiagonal matrix.
    real(kind=dp) :: T(size(U), size(U)-1)
    real(kind=dp) :: beta
    !> Krylov subspace dimension.
    integer :: kdim
    !> Miscellaneous.
    integer :: i, j, k
    integer(int_size) :: indices(size(U)-1)

    ! --> Deals with the optional args.
    verbose = optval(verbosity, .false.)
    ! --> Assert size(U) == size(V).
    if (size(U) .ne. size(V)) then
       info = -1
       if (verbose) then
          write(*, *) "INFO : Left and Right Krylov subspaces have different dimensions."
          write(*, *) "       Exiting svds with exit code info =", info
       endif
    endif

    ! --> Initialize variables.
    T = 0.0_dp ; residuals = 0.0_dp ; singvecs = 0.0_dp ; singvals = 0.0_dp
    do i = 2, size(U)
       call U(i)%zero() ; call V(i)%zero()
    enddo

    return
  end subroutine svds

  !--------------------------------------------
  !-----                                  -----
  !-----     ITERATIVE LINEAR SOLVERS     -----
  !-----                                  -----
  !--------------------------------------------

  subroutine gmres(A, b, x, info, kdim, maxiter, tol, verbosity)
    !> Linear problem.
    class(abstract_linop) , intent(in) :: A ! Linear Operator.
    class(abstract_vector), intent(in) :: b ! Right-hand side.
    !> Solution vector.
    class(abstract_vector), intent(inout) :: x
    !> Information flag.
    integer                            , intent(out)   :: info
    !> Optional arguments.
    integer, optional, intent(in)          :: kdim      ! Krylov subspace dimension.
    integer                                :: k_dim
    integer, optional, intent(in)          :: maxiter   ! Maximum number full GMRES iterations.
    integer                                :: niter
    double precision, optional, intent(in) :: tol       ! Tolerance for the GMRES residual.
    double precision                       :: tolerance
    logical, optional, intent(in)          :: verbosity ! Verbosity control.
    logical                                :: verbose

    !> Krylov subspace.
    class(abstract_vector), dimension(:)   , allocatable :: V
    !> Upper Hessenberg matrix.
    double precision      , dimension(:, :), allocatable :: H
    !> Least-squares related variables.
    double precision      , dimension(:)   , allocatable :: y
    double precision      , dimension(:)   , allocatable :: e
    double precision                                     :: beta

    !> Miscellaneous.
    integer :: i, j, k, l, m
    double precision :: alpha
    class(abstract_vector), allocatable :: dx

    ! --> Deals with the optional arguments.
    k_dim     = optval(kdim, 30)
    niter     = optval(maxiter, 10)
    tolerance = optval(tol, 1.0D-12)
    verbose   = optval(verbosity, .false.)

    ! --> Initialize Krylov subspace.
    allocate(V(1:k_dim+1), source=b)
    do i = 1, size(V)
       call V(i)%zero()
    enddo
    allocate(H(k_dim+1, k_dim)) ; H = 0.0D+00
    allocate(y(1:k_dim))        ; y = 0.0D+00
    allocate(e(1:k_dim+1))      ; e = 0.0D+00

    ! --> Initial Krylov vector.
    call A%matvec(x, V(1)) ; call V(1)%sub(b) ; call V(1)%scalar_mult(-1.0D+00)
    beta = V(1)%norm() ; call V(1)%scalar_mult(1.0D+00 / beta)

    gmres_iterations : do i = 1, niter
       ! --> Zero-out variables.
       H = 0.0D+00 ; y = 0.0D+00 ; e = 0.0D+00 ; e(1) = beta
       do j = 2, size(V)
          call V(j)%zero()
       enddo

       arnoldi : do k = 1, k_dim
          ! --> Step-by-step Arnoldi factorization.
          call arnoldi_factorization(A, V, H, info, kstart=k, kend=k, verbosity=.false., tol=tolerance)
          if (info < 0) then
             write(*, *) "INFO : Arnoldi Factorization failed with exit code info =", info
             write(*, *) "       Stopping the GMRES computation."
          endif
          ! --> Least-squares problem.
          call lstsq(H(1:k+1, 1:k), e(1:k+1), y(1:k))
          ! --> Compute residual.
          beta = norm2(e(1:k+1) - matmul(H(1:k+1, 1:k), y(1:k)))
          if (verbose) then
             write(*, *) "INFO : GMRES residual after ", (i-1)*k_dim + k, "iteration : ", beta
          endif
          ! --> Check convergence.
          if (beta**2 .lt. tolerance) then
             exit arnoldi
          endif
       enddo arnoldi

       ! --> Update solution.
       k = min(k, k_dim)
       if (allocated(dx) .eqv. .false.) allocate(dx, source=x)
       call dx%zero() ; call get_vec(dx, V(1:k), y(1:k)) ; call x%add(dx)

       ! --> Recompute residual for sanity check.
       call A%matvec(x, V(1)) ; call V(1)%sub(b) ; call V(1)%scalar_mult(-1.0D+00)

       ! --> Initialize new starting Krylov vector if needed.
       beta = V(1)%norm() ; call V(1)%scalar_mult(1.0D+00 / beta)

       if (beta**2 .lt. tolerance) then
          exit gmres_iterations
       endif

    enddo gmres_iterations

    ! --> Deallocate variables.
    deallocate(V, H, y, e)

    return
  end subroutine gmres

  subroutine lstsq(A, b, x)
    !> Input matrix.
    double precision, dimension(:, :), intent(in)  :: A
    double precision, dimension(:)   , intent(in)  :: b
    double precision, dimension(:)   , intent(out) :: x

    !> Lapack job.
    character :: trans = "N"
    integer   :: m, n, nrhs, lda, ldb, lwork, info
    double precision, dimension(size(A, 1), size(A, 2)) :: A_tilde
    double precision, dimension(size(A, 1))             :: b_tilde
    double precision, dimension(:), allocatable         :: work

    !> Interface to LAPACK dgels
    interface
       pure subroutine dgels(ftrans, fm, fn, fnrhs, fA, flda, fb, fldb, fwork, flwork, finfo)
         character, intent(in)           :: ftrans
         integer  , intent(in)           :: fm, fn, fnrhs, flda, fldb, flwork, finfo
         double precision, intent(inout) :: fa(flda, *)
         double precision, intent(inout) :: fb(flda, *)
         double precision, intent(out)   :: fwork(*)
       end subroutine dgels
    end interface

    !> Initialize variables.
    m = size(A, 1) ; n = size(A, 2) ; nrhs = 1
    lda = m ; ldb = m ; lwork = max(1, min(m, n) + max(min(m, n), nrhs))
    A_tilde = A ; b_tilde = b
    allocate(work(1:lwork)) ; work = 0.0D+00

    !> Solve the least-squares problem.
    call dgels(trans, m, n, nrhs, A_tilde, lda, b_tilde, ldb, work, lwork, info)

    !> Return solution.
    x = b_tilde(1:n)

    return
  end subroutine lstsq

  ! subroutine cg(A, b)
  !   !> Linear problem.
  !   class(abstract_spd_linop), intent(in) :: A ! Linear Operator.
  !   class(abstract_vector)   , intent(in) :: b ! Right-hand side.
  !   return
  ! end subroutine cg

end module IterativeSolvers
