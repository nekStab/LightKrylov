module IterativeSolvers
  use AbstractVector
  use LinearOperator
  use KrylovDecomp
  use stdlib_sorting, only : sort_index, int_size
  use stdlib_optval , only : optval
  use stdlib_io_npy , only : save_npy
  implicit none
  include "dtypes.h"

  private
  public :: eigs, eighs, gmres, save_eigenspectrum, svds

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
    real(kind=wp), intent(in) :: real_part(:)
    real(kind=wp), intent(in) :: imag_part(:)
    !> Residual norm computed from the Arnoldi/Lanczos factorization.
    real(kind=wp), intent(in) :: residuals(:)
    !> Name of the output file.
    character(len=*), intent(in) :: filename

    !> Miscellaneous.
    real(kind=wp), dimension(size(real_part), 3) :: data

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
    endif

    ! --> Compute spectral decomposition of the Hessenberg matrix.
    call evd(H(1:kdim, 1:kdim), eigvecs, eigvals, kdim)

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

  subroutine svds(A, U, V, uvecs, vvecs, sigma, residuals, info, verbosity)
    !> Linear Operator.
    class(abstract_linop), intent(in) :: A
    !> Krylov bases.
    class(abstract_vector), intent(inout) :: U(:) ! Basis for left sing. vectors.
    class(abstract_vector), intent(inout) :: V(:) ! Basis for right sing. vectors.
    !> Coordinates of singular vectors in Krylov bases, singular values, and associated residuals.
    double precision, intent(out) :: uvecs(size(U), size(U)-1), vvecs(size(U)-1, size(U)-1)
    double precision, intent(out) :: sigma(size(U)-1)
    double precision, intent(out) :: residuals(size(U)-1)
    !> Information flag.
    integer, intent(out) :: info
    !> Verbosity control.
    logical, optional, intent(in) :: verbosity
    logical verbose

    !> Bidiagonal matrix.
    double precision :: B(size(U), size(U)-1)
    double precision :: beta
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
    else
       kdim = size(U)-1
    endif

    ! --> Initialize variables.
    B = 0.0D+00 ; residuals = 0.0D+00 ; uvecs = 0.0D+00 ; vvecs = 0.0D+00 ; sigma = 0.0D+00
    do i = 2, size(U)
       call U(i)%zero() ; call V(i)%zero()
    enddo

    ! --> Compute the Lanczos bidiagonalization.
    call lanczos_bidiagonalization(A, U, V, B, info, verbosity=verbose)

    ! --> Compute the singular value decomposition of the bidiagonal matrix.
    call svd(B, uvecs, sigma, vvecs)

    ! --> Compute the residual associated with each singular triplet.
    residuals = 0.0D+00

    return
  end subroutine svds

  subroutine svd(A, U, S, V)
    !> Matrix to be factorized.
    double precision, intent(in)  :: A(:, :)
    !> Left singular vectors.
    double precision, intent(out) :: U(size(A, 1), min(size(A, 1), size(A, 2)))
    !> Singular values.
    double precision, intent(out) :: S(size(A, 2))
    !> Right singular vectors.
    double precision, intent(out) :: V(size(A, 2), min(size(A, 1), size(A, 2)))

    !> Lapack job.
    character                  :: jobu = "S", jobvt = "S"
    integer                    :: m, n, lda, ldu, ldvt, lwork, info
    double precision, allocatable :: work(:)
    double precision :: A_tilde(size(A, 1), size(A, 2)), vt(min(size(A, 1), size(A, 2)), size(A, 2))

    interface
       pure subroutine dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
         character, intent(in) :: jobu, jobvt
         integer  , intent(in) :: m, n, lda, ldu, ldvt, lwork, info
         double precision, intent(inout) :: a(lda, *)
         double precision, intent(out)   :: u(ldu, *), s(*), vt(ldvt, *), work(*)
       end subroutine dgesvd
    end interface

    m = size(A, 1) ; n = size(A, 2)
    lda = size(A, 1) ; ldu = size(A, 1) ; ldvt = size(A, 2)
    lwork = max(1, 3*min(m, n), 5*min(m, n)) ; allocate(work(lwork))

    a_tilde = a
    call dgesvd(jobu, jobvt, m, n, a_tilde, lda, s, u, ldu, vt, ldvt, work, lwork, info)
    v = transpose(vt)

    return
  end subroutine svd

  !--------------------------------------------
  !-----                                  -----
  !-----     ITERATIVE LINEAR SOLVERS     -----
  !-----                                  -----
  !--------------------------------------------

  !=======================================================================================
  ! Generalized Minimal Residual (GMRES) Solver Subroutine
  !=======================================================================================
  !
  ! Purpose:
  ! --------
  ! Implements the classic, unpreconditioned Generalized Minimal Residual (GMRES) algorithm
  ! for solving nonsymmetric, non-Hermitian linear systems of equations, Ax = b.
  !
  ! Algorithmic Features:
  ! ----------------------
  ! - Constructs a full Krylov subspace without restarts (i.e., not GMRES(m)).
  ! - Utilizes Arnoldi factorization to generate an orthonormal basis for the Krylov subspace.
  ! - Employs a least-squares solve to determine the optimal linear combination of the Krylov vectors.
  ! - Updates the approximate solution based on the least-squares solution.
  !
  ! Advantages:
  ! -----------
  ! - Suitable for nonsymmetric and ill-conditioned matrices.
  ! - Produces monotonically decreasing residuals.
  ! - Fully utilizes the generated Krylov subspace for the solution.
  !
  ! Limitations:
  ! ------------
  ! - Memory-intensive due to the absence of restarts.
  ! - May not be efficient for very large-scale problems.
  ! - No preconditioning capabilities in the current implementation.
  !
  ! Input/Output Parameters:
  ! ------------------------
  ! - A        : Linear Operator              [Input]
  ! - b        : Right-hand side vector       [Input]
  ! - x        : Initial/Updated solution     [Input/Output]
  ! - info     : Iteration Information flag   [Output]
  ! - maxiter  : Maximum number of iterations [Optional, Input]
  ! - tol      : Tolerance for convergence    [Optional, Input]
  ! - verbosity: Verbosity control flag       [Optional, Input]
  !
  ! References:
  ! -----------
  ! - Saad, Y., and Schultz, M. H. (1986). "GMRES: A Generalized Minimal Residual Algorithm for Solving Nonsymmetric Linear Systems,"
  !   SIAM Journal on Scientific and Statistical Computing, 7(3), 856–869.
  !
  !=======================================================================================
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
    real(kind=wp), optional, intent(in)    :: tol       ! Tolerance for the GMRES residual.
    real(kind=wp)                          :: tolerance
    logical, optional, intent(in)          :: verbosity ! Verbosity control.
    logical                                :: verbose

    !> Krylov subspace.
    class(abstract_vector), dimension(:)   , allocatable :: V
    !> Upper Hessenberg matrix.
    real(kind=wp)         , dimension(:, :), allocatable :: H
    !> Least-squares related variables.
    real(kind=wp)         , dimension(:)   , allocatable :: y
    real(kind=wp)         , dimension(:)   , allocatable :: e
    real(kind=wp)                                        :: beta

    !> Miscellaneous.
    integer                             :: i, j, k, l, m
    real(kind=wp)                       :: alpha
    class(abstract_vector), allocatable :: dx

    ! --> Deals with the optional arguments.
    k_dim     = optval(kdim, 30)
    niter     = optval(maxiter, 10)
    tolerance = optval(tol, atol + rtol*b%norm())
    verbose   = optval(verbosity, .false.)

    ! --> Initialize Krylov subspace.
    allocate(V(1:k_dim+1), source=b)
    do i = 1, size(V)
       call V(i)%zero()
    enddo
    allocate(H(k_dim+1, k_dim)) ; H = 0.0_wp
    allocate(y(1:k_dim))        ; y = 0.0_wp
    allocate(e(1:k_dim+1))      ; e = 0.0_wp

    ! --> Initial Krylov vector.
    call A%matvec(x, V(1)) ; call V(1)%sub(b) ; call V(1)%scal(-1.0_wp)
    beta = V(1)%norm() ; call V(1)%scal(1.0_wp / beta)

    gmres_iterations : do i = 1, niter
       ! --> Zero-out variables.
       H = 0.0_wp ; y = 0.0_wp ; e = 0.0_wp ; e(1) = beta
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
       call A%matvec(x, V(1)) ; call V(1)%sub(b) ; call V(1)%scal(-1.0_wp)

       ! --> Initialize new starting Krylov vector if needed.
       beta = V(1)%norm() ; call V(1)%scal(1.0D+00 / beta)

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
    real(kind=wp), dimension(:, :), intent(in)  :: A
    real(kind=wp), dimension(:)   , intent(in)  :: b
    real(kind=wp), dimension(:)   , intent(out) :: x

    !> Lapack job.
    character :: trans = "N"
    integer   :: m, n, nrhs, lda, ldb, lwork, info
    real(kind=wp), dimension(size(A, 1), size(A, 2)) :: A_tilde
    real(kind=wp), dimension(size(A, 1))             :: b_tilde
    real(kind=wp), dimension(:), allocatable         :: work

    !> Interface to LAPACK dgels
    interface
       pure subroutine dgels(ftrans, fm, fn, fnrhs, fA, flda, fb, fldb, fwork, flwork, finfo)
         import wp
         character, intent(in)           :: ftrans
         integer  , intent(in)           :: fm, fn, fnrhs, flda, fldb, flwork, finfo
         real(kind=wp), intent(inout)    :: fa(flda, *)
         real(kind=wp), intent(inout)    :: fb(flda, *)
         real(kind=wp), intent(out)      :: fwork(*)
       end subroutine dgels
    end interface

    !> Initialize variables.
    m = size(A, 1) ; n = size(A, 2) ; nrhs = 1
    lda = m ; ldb = m ; lwork = max(1, min(m, n) + max(min(m, n), nrhs))
    A_tilde = A ; b_tilde = b
    allocate(work(1:lwork)) ; work = 0.0_wp

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
