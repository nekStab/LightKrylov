module IterativeSolvers
  use Utils
  use AbstractVector
  use LinearOperator
  use BaseKrylov
  use stdlib_sorting, only : sort_index, int_size
  use stdlib_optval , only : optval
  use stdlib_io_npy , only : save_npy
  implicit none
  include "dtypes.h"

  private
  public :: eigs, eighs, gmres, save_eigenspectrum, svds, cg
  public :: abstract_linear_solver

  abstract interface
     subroutine abstract_linear_solver(A, b, x, info, options, transpose)
       import abstract_linop, abstract_vector, abstract_opts
       !> Linear problem.
       class(abstract_linop) , intent(in) :: A
       class(abstract_vector), intent(in) :: b
       !> Solution vector.
       class(abstract_vector), intent(inout) :: x
       !> Information flag.
       integer               , intent(out) :: info
       !> Solver options.
       class(abstract_opts)  , optional, intent(in) :: options
       !> Transposition flag.
       logical               , optional, intent(in) :: transpose
     end subroutine abstract_linear_solver
  end interface

contains

  !-----------------------------
  !-----     UTILITIES     -----
  !-----------------------------

  elemental pure function compute_residual(beta, x) result(residual)
    !> Norm of Krylov residual vector.
    real(kind=wp), intent(in) :: beta
    !> Last element of Ritz eigenvector.
    real(kind=wp), intent(in) :: x
    !> Residual.
    real(kind=wp) :: residual

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

  subroutine eigs(A, X, eigvecs, eigvals, residuals, info, verbosity, transpose)
    !> Linear Operator.
    class(abstract_linop), intent(in) :: A
    !> Krylov basis.
    class(abstract_vector), dimension(:), intent(inout) :: X
    !> Coordinates of eigenvectors in Krylov basis, eigenvalues and associated residuals.
    complex(kind=wp), dimension(size(X)-1, size(X)-1), intent(out) :: eigvecs
    complex(kind=wp), dimension(size(X)-1)           , intent(out) :: eigvals
    real(kind=wp)   , dimension(size(X)-1)           , intent(out) :: residuals
    !> Information flag.
    integer, intent(out) :: info
    !> Verbosity control.
    logical, optional, intent(in) :: verbosity
    logical :: verbose
    !> Transpose operator.
    logical, optional, intent(in) :: transpose
    logical :: trans

    !> Upper Hessenberg matrix.
    real(kind=wp), dimension(size(X), size(X)-1) :: H
    real(kind=wp)                                :: beta
    !> Krylov subspace dimension.
    integer :: kdim
    !> Miscellaneous.
    integer :: i, j, k
    integer(int_size), dimension(size(X)-1) :: indices
    real(kind=wp)    , dimension(size(X)-1) :: abs_eigvals

    ! --> Dimension of the Krylov subspace.
    kdim = size(X) - 1

    ! --> Deals with the optional arguments.
    verbose = optval(verbosity, .false.)
    trans   = optval(transpose, .false.)

    ! --> Initialize variables.
    H = 0.0_wp ; residuals = 0.0_wp ; eigvals = (0.0_wp, 0.0_wp) ; eigvecs = (0.0_wp, 0.0_wp)
    do i = 2, size(X) ! i=1 is the initial Krylov vector given by the user.
       call X(i)%zero()
    enddo

    ! --> Compute Arnoldi factorization.
    call arnoldi_factorization(A, X, H, info, verbosity=verbose, transpose=trans)

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
    real(kind=wp), dimension(size(X)-1, size(X)-1), intent(out) :: eigvecs
    real(kind=wp), dimension(size(X)-1)           , intent(out) :: eigvals
    real(kind=wp), dimension(size(X)-1)           , intent(out) :: residuals
    !> Information flag.
    integer, intent(out) :: info
    !> Verbosity control.
    logical, optional, intent(in) :: verbosity
    logical :: verbose

    !> Tridiagonal matrix.
    real(kind=wp), dimension(size(X), size(X)-1) :: T
    real(kind=wp)                                :: beta
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
    T = 0.0_wp ; residuals = 0.0_wp ; eigvecs = 0.0_wp ; eigvals = 0.0_wp
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
    real(kind=wp), dimension(n, n) :: A, A_tilde, vr
    real(kind=wp), dimension(1, n) :: vl
    real(kind=wp), dimension(4*n)  :: work
    real(kind=wp), dimension(n)    :: wr, wi
    complex(kind=wp), dimension(n, n)   :: vecs
    complex(kind=wp), dimension(n)      :: vals
    integer :: i
    integer, dimension(n) :: idx

    interface
       pure subroutine dgeev(fjobvl, fjobvr, fn, fa, flda, fwr, fwi, fvl, fldvl, fvr, fldvr, fwork, flwork, finfo)
         import wp
         character, intent(in) :: fjobvl, fjobvr
         integer  , intent(in) :: fn, flda, fldvl, fldvr, flwork, finfo
         real(kind=wp), intent(inout) :: fa(flda, *)
         real(kind=wp), intent(out) :: fwr(fn), fwi(fn), fvl(fldvl, *), fvr(fldvr, *), fwork(flwork)
       end subroutine dgeev
    end interface

    ! --> Compute the eigendecomposition of A.
    lda = n ; ldvl = 1 ; ldvr = n ; lwork = 4*n ; A_tilde = A
    call dgeev(jobvl, jobvr, n, A_tilde, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)

    ! --> Real to complex arithmetic.
    !     NOTE : Check if a LAPACK function already exists for that purpose.
    vals = wr*(1.0_wp, 0.0_wp) + wi*(0.0_wp, 1.0_wp)
    vecs = vr*(1.0_wp, 0.0_wp)

    do i = 1, n-1
       if (wi(i) .gt. 0) then
          vecs(:, i) = vr(:, i)*(1.0_wp, 0.0_wp) + vr(:, i+1)*(0.0_wp, 1.0_wp)
          vecs(:, i+1) = conjg(vecs(:, i))
       else if (abs(wi(i)) .le. epsilon(wi(i))) then
          vecs(:, i) = vr(:, i) * (1.0_wp, 0.0_wp)
       endif
    enddo

    return
  end subroutine evd

  subroutine hevd(A, vecs, vals, n)
    !> Lapack job.
    character :: jobz="V", uplo="U"
    integer :: n, lwork, info, lda
    real(kind=wp), dimension(n) :: vals
    real(kind=wp), dimension(n, n) :: A, A_tilde, vecs
    real(kind=wp), dimension(3*n-1) :: work
    integer :: i, j, k

    interface
       pure subroutine dsyev(fjobz, fuplo, fn, fa, flda, fw, fwork, flwork, finfo)
         import wp
         character, intent(in) :: fjobz, fuplo
         integer, intent(in)  :: fn
         integer, intent(in)  :: flda
         integer, intent(in)  :: flwork
         integer, intent(out) :: finfo
         real(kind=wp), intent(inout) :: fa(flda, *)
         real(kind=wp), intent(out)   :: fw(*)
         real(kind=wp), intent(out)   :: fwork(*)
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
    real(kind=wp), intent(out) :: uvecs(size(U)-1, size(U)-1), vvecs(size(U)-1, size(U)-1)
    real(kind=wp), intent(out) :: sigma(size(U)-1)
    real(kind=wp), intent(out) :: residuals(size(U)-1)
    !> Information flag.
    integer, intent(out) :: info
    !> Verbosity control.
    logical, optional, intent(in) :: verbosity
    logical verbose

    !> Bidiagonal matrix.
    real(kind=wp) :: B(size(U), size(U)-1)
    real(kind=wp) :: beta
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
    B = 0.0_wp ; residuals = 0.0_wp ; uvecs = 0.0_wp ; vvecs = 0.0_wp ; sigma = 0.0_wp
    do i = 2, size(U)
       call U(i)%zero() ; call V(i)%zero()
    enddo

    ! --> Compute the Lanczos bidiagonalization.
    call lanczos_bidiagonalization(A, U, V, B, info, verbosity=verbose)

    ! --> Compute the singular value decomposition of the bidiagonal matrix.
    call svd(B(1:kdim, 1:kdim), uvecs, sigma, vvecs)

    ! --> Compute the residual associated with each singular triplet.
    beta = B(kdim+1, kdim) !> Get Krylov residual vector norm.
    residuals = compute_residual(beta, vvecs(kdim, :))

    return
  end subroutine svds

  subroutine svd(A, U, S, V)
    !> Matrix to be factorized.
    real(kind=wp), intent(in)  :: A(:, :)
    !> Left singular vectors.
    real(kind=wp), intent(out) :: U(size(A, 1), min(size(A, 1), size(A, 2)))
    !> Singular values.
    real(kind=wp), intent(out) :: S(size(A, 2))
    !> Right singular vectors.
    real(kind=wp), intent(out) :: V(size(A, 2), min(size(A, 1), size(A, 2)))

    !> Lapack job.
    character                  :: jobu = "S", jobvt = "S"
    integer                    :: m, n, lda, ldu, ldvt, lwork, info
    real(kind=wp), allocatable :: work(:)
    real(kind=wp) :: A_tilde(size(A, 1), size(A, 2)), vt(min(size(A, 1), size(A, 2)), size(A, 2))

    interface
       pure subroutine dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
         import wp
         character, intent(in) :: jobu, jobvt
         integer  , intent(in) :: m, n, lda, ldu, ldvt, lwork, info
         real(kind=wp), intent(inout) :: a(lda, *)
         real(kind=wp), intent(out)   :: u(ldu, *), s(*), vt(ldvt, *), work(*)
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
  ! ---------------------
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
  subroutine gmres(A, b, x, info, options, transpose)
    !> Linear problem.
    class(abstract_linop)  , intent(in)    :: A ! Linear Operator.
    class(abstract_vector) , intent(in)    :: b ! Right-hand side.
    !> Solution vector.
    class(abstract_vector) , intent(inout) :: x
    !> Information flag.
    integer                , intent(out)   :: info
    !> Optional arguments.
    class(abstract_opts), optional, intent(in) :: options
    type(gmres_opts)                       :: opts
    logical, optional      , intent(in)    :: transpose
    logical                                :: trans
    !> Internal variables.
    integer                                :: k_dim
    integer                                :: maxiter
    real(kind=wp)                          :: tol
    logical                                :: verbose
    !> Krylov subspace.
    class(abstract_vector), allocatable :: V(:)
    !> Upper Hessenberg matrix.
    real(kind=wp)         , allocatable :: H(:, :)
    !> Least-squares related variables.
    real(kind=wp)         , allocatable :: y(:)
    real(kind=wp)         , allocatable :: e(:)
    real(kind=wp)                       :: beta
    !> Miscellaneous.
    integer                             :: i, j, k, l, m
    real(kind=wp)                       :: alpha
    class(abstract_vector), allocatable :: dx

    ! --> Deals with the optional arguments.
    if (present(options)) then
       select type(options)
       type is(gmres_opts)
          opts = gmres_opts(              &
               kdim    = options%kdim,    &
               maxiter = options%maxiter, &
               atol    = options%atol,    &
               rtol    = options%rtol,    &
               verbose = options%verbose  &
               )
       end select
    else
       opts = gmres_opts()
    end if
    k_dim = opts%kdim ; maxiter = opts%maxiter
    tol = opts%atol + opts%rtol * b%norm() ; verbose = opts%verbose
    trans = optval(transpose, .false.)

    ! --> Initialize Krylov subspace.
    allocate(V(1:k_dim+1), source=b)
    do i = 1, size(V)
       call V(i)%zero()
    enddo
    allocate(H(k_dim+1, k_dim)) ; H = 0.0_wp
    allocate(y(1:k_dim))        ; y = 0.0_wp
    allocate(e(1:k_dim+1))      ; e = 0.0_wp

    ! --> Initial Krylov vector.
    if (trans) then
       call A%rmatvec(x, V(1))
    else
       call A%matvec(x, V(1))
    endif
    call V(1)%sub(b) ; call V(1)%scal(-1.0_wp)
    beta = V(1)%norm() ; call V(1)%scal(1.0_wp / beta)

    gmres_iterations : do i = 1, maxiter
       ! --> Zero-out variables.
       H = 0.0_wp ; y = 0.0_wp ; e = 0.0_wp ; e(1) = beta
       do j = 2, size(V)
          call V(j)%zero()
       enddo

       arnoldi : do k = 1, k_dim
          ! --> Step-by-step Arnoldi factorization.
          call arnoldi_factorization(A, V, H, info, kstart=k, kend=k, verbosity=.false., tol=tol, transpose=trans)
          if (info < 0) then
             write(*, *) "INFO : Arnoldi Factorization failed with exit code info =", info
             write(*, *) "       Stopping the GMRES computation."
             call exit()
          endif
          ! --> Least-squares problem.
          call lstsq(H(1:k+1, 1:k), e(1:k+1), y(1:k))
          ! --> Compute residual.
          beta = norm2(e(1:k+1) - matmul(H(1:k+1, 1:k), y(1:k)))
          if (verbose) then
             write(*, *) "INFO : GMRES residual after ", (i-1)*k_dim + k, "iteration : ", beta
          endif
          ! --> Check convergence.
          if (beta**2 .lt. tol) then
             exit arnoldi
          endif
       enddo arnoldi

       ! --> Update solution.
       k = min(k, k_dim)
       if (allocated(dx) .eqv. .false.) allocate(dx, source=x)
       call dx%zero() ; call get_vec(dx, V(1:k), y(1:k)) ; call x%add(dx)

       ! --> Recompute residual for sanity check.
       if (trans) then
          call A%rmatvec(x, V(1))
       else
          call A%matvec(x, V(1))
       endif

       call V(1)%sub(b) ; call V(1)%scal(-1.0_wp)

       ! --> Initialize new starting Krylov vector if needed.
       beta = V(1)%norm() ; call V(1)%scal(1.0D+00 / beta)

       if (beta**2 .lt. tol) then
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

  !=======================================================================================
  ! Conjugate Gradient (CG) Solver Subroutine
  !=======================================================================================
  !
  ! Purpose:
  ! --------
  ! Implements the classic Conjugate Gradient (CG) algorithm for solving symmetric positive
  ! definite (SPD) linear systems of equations, Ax = b.
  !
  ! Algorithmic Features:
  ! ----------------------
  ! - Utilizes the method of conjugate directions to iteratively refine the solution.
  ! - Employs two sequences of vectors: residuals (r) and conjugate directions (p).
  ! - Updates the approximate solution based on the computed alpha and beta values.
  !
  ! Advantages:
  ! -----------
  ! - Well-suited for large, sparse, symmetric positive definite (SPD) matrices.
  ! - Memory-efficient, requiring storage for only a few vectors in comparison to GMRES.
  ! - Under exact arithmetic, finds the exact solution within 'n' iterations for an 'n'-dimensional SPD matrix.
  !
  ! Limitations:
  ! ------------
  ! - Applicability restricted to SPD matrices.
  ! - No preconditioning capabilities in the current implementation.
  ! - Subject to numerical rounding errors, which might require more than 'n' iterations in practice.
  !
  ! Input/Output Parameters:
  ! ------------------------
  ! - A        : Linear Operator (SPD)        [Input]
  ! - b        : Right-hand side vector       [Input]
  ! - x        : Initial/Updated solution     [Input/Output]
  ! - info     : Iteration Information flag   [Output]
  ! - maxiter  : Maximum number of iterations [Optional, Input]
  ! - tol      : Tolerance for convergence    [Optional, Input]
  ! - verbosity: Verbosity control flag       [Optional, Input]
  !
  ! References:
  ! -----------
  ! - Hestenes, M. R., and Stiefel, E. (1952). "Methods of Conjugate Gradients for Solving Linear Systems,"
  !   Journal of Research of the National Bureau of Standards, 49(6), 409–436.
  !
  !=======================================================================================
  subroutine cg(A, b, x, info, options)
    !> Linear problem.
    class(abstract_spd_linop), intent(in) :: A ! Linear Operator.
    class(abstract_vector), intent(in) :: b ! Right-hand side.
    !> Solution vector.
    class(abstract_vector), intent(inout) :: x
    !> Information flag.
    integer, intent(out)   :: info
    !> Optional arguments.
    class(abstract_opts), optional, intent(in) :: options
    type(cg_opts)                              :: opts
    integer :: maxiter ! Maximum number of CG iterations.
    real(kind=wp) :: tol  ! Tolerance for the CG residual.
    logical                             :: verbose
    
    !> Residual and direction vectors.
    class(abstract_vector), allocatable :: r, p, Ap
    !> Scalars used in the CG algorithm.
    real(kind=wp) :: alpha, beta, r_dot_r_old, r_dot_r_new, residual
    integer :: i, j, k
    
    ! --> Handle optional arguments.
    if (present(options)) then
       select type(options)
       type is(cg_opts)
          opts = cg_opts(       &
               maxiter = options%maxiter, &
               atol    = options%atol   , &
               rtol    = options%rtol   , &
               verbose = options%verbose  &
               )
       end select
    else
       opts = cg_opts()
    endif
    tol = opts%atol + opts%rtol * b%norm() ; maxiter = opts%maxiter ; verbose = opts%verbose
    
    ! --> Initialize vectors.
    allocate (r, source=b)  ; call r%zero()
    allocate (p, source=b)  ; call p%zero()
    allocate (Ap, source=b) ; call Ap%zero()
    
    ! --> Compute initial residual: r = b - Ax.
    call A%matvec(x, r) ; call r%axpby(-1.0_wp, b, 1.0_wp)
    
    ! --> Initialize direction vector: p = r.
    p = r
    
    ! --> Initialize dot product of residual: r_dot_r_old = r' * r.
    r_dot_r_old = r%dot(r)
    
    ! --> CG Iteration Loop.
    cg_iterations: do i = 1, maxiter
       
       ! Compute A * p.
       call A%matvec(p, Ap)
       
       ! Compute step size alpha = r_dot_r_old / (p' * Ap).
       alpha = r_dot_r_old / p%dot(Ap)
       
       ! Update solution x = x + alpha * p.
       call x%axpby(1.0_wp, p, alpha)
       
       ! Update residual r = r - alpha * Ap.
       call r%axpby(1.0_wp, Ap, -alpha)
       
       ! Compute new dot product of residual r_dot_r_new = r' * r.
       r_dot_r_new = r%dot(r)
       
       ! Check for convergence.
       residual = sqrt(r_dot_r_new)
       
       if (verbose) then
          write (*, *) "INFO : CG residual after ", (i), "iterations : ", residual
       end if
       
       if (residual < tol) then
          if (verbose) then
             write (*, *) "INFO : CG Converged: residual ", residual, "< tolerance: ", tol
          end if
          exit cg_iterations
       end if
       
       ! Compute new direction beta = r_dot_r_new / r_dot_r_old.
       beta = r_dot_r_new/r_dot_r_old
       
       ! Update direction p = r + beta * p.
       call p%axpby(beta, r, 1.0_wp)
       
       ! Update r_dot_r_old for next iteration.
       r_dot_r_old = r_dot_r_new
       
    end do cg_iterations
    
    ! ! --> Set info flag.
    ! if (residual < tolerance) then
    !     info = 0
    ! else
    !     info = 1
    ! end if
    deallocate (r, p, Ap)

    return
  end subroutine cg

  ! !=======================================================================================
  ! ! Biconjugate Gradient Stabilized (BiCGSTAB) Solver Subroutine
  ! !=======================================================================================
  ! !
  ! ! Purpose:
  ! ! --------
  ! ! Implements the Biconjugate Gradient Stabilized (BiCGSTAB) algorithm for solving
  ! ! nonsymmetric and possibly ill-conditioned linear systems Ax = b.
  ! !
  ! ! Algorithmic Features:
  ! ! ----------------------
  ! ! - Extends the BiCG algorithm by stabilizing the iterations.
  ! ! - Utilizes two search directions and two residuals to improve stability.
  ! ! - Iteratively updates both the approximate solution and the residuals.
  ! !
  ! ! Advantages:
  ! ! -----------
  ! ! - Capable of addressing nonsymmetric and ill-conditioned matrices.
  ! ! - Generally more stable and faster compared to the basic BiCG.
  ! ! - Suitable for large and sparse matrices.
  ! !
  ! ! Limitations:
  ! ! ------------
  ! ! - May experience stagnation for certain types of problems.
  ! ! - No preconditioning capabilities in the current implementation.
  ! !
  ! ! Input/Output Parameters:
  ! ! ------------------------
  ! ! - A        : Linear Operator (abstract_linop) [Input]
  ! ! - b        : Right-hand side (abstract_vector) [Input]
  ! ! - x        : Initial/Updated solution (abstract_vector) [Input/Output]
  ! ! - info     : Iteration Information flag (Integer) [Output]
  ! ! - maxiter  : Maximum number of iterations (Integer) [Optional, Input]
  ! ! - tol      : Convergence tolerance (real(kind=dp)) [Optional, Input]
  ! ! - verbosity: Verbosity control flag (Logical) [Optional, Input]
  ! !
  ! ! References:
  ! ! -----------
  ! ! - van der Vorst, H. A. (1992). "Bi-CGSTAB: A Fast and Smoothly Converging Variant of Bi-CG for the Solution of Nonsymmetric Linear Systems,"
  ! !   SIAM Journal on Scientific and Statistical Computing, 13(2), 631–644.
  ! !
  ! !=======================================================================================
  ! subroutine bicgstab(A, b, x, info, options, transpose)
  !   !> Linear problem and initial guess.
  !   class(abstract_linop), intent(in) :: A
  !   class(abstract_vector), intent(in) :: b
  !   class(abstract_vector), intent(inout) :: x
  !   !> Output and optional input parameters.
  !   integer, intent(out) :: info
  !   class(abstract_opts), optional, intent(in) :: options
  !   type(bicgstab_opts)                        :: opts
  !   logical, optional, intent(in) :: transpose
    
  !   !> Internal variables.
  !   integer :: i, maxiter
  !   real(kind=wp) :: tol, res, alpha, omega, rho, rho_new, beta
  !   logical :: verbose, trans
    
  !   !> BiCGSTAB vectors.
  !   class(abstract_vector), allocatable :: r, r_hat, p, p_int, v, s, t

  !   ! --> Deals with the optional arguments.
  !   if (present(options)) then
  !      select type(options)
  !      type is(bicgstab_opts)
  !         opts = bicgstab_opts( &
  !              maxiter = options%maxiter, &
  !              atol    = options%atol,    &
  !              rtol    = options%rtol,    &
  !              verbose = options%verbose  &
  !         )
  !      end select
  !   else
  !      opts = bicgstab_opts()
  !   end if
  !   maxiter = opts%maxiter ; tol = opts%atol + opts%rtol * b%norm()
  !   verbose = opts%verbose ; trans = optval(transpose, .false.)
    
  !   ! Initialize vectors.
  !   allocate (r, source=b)     ; call r%zero()
  !   allocate (r_hat, source=b) ; call r_hat%zero()
  !   allocate (p, source=b)     ; call p%zero()
  !   allocate (v, source=b)     ; call v%zero()
  !   allocate (s, source=b)     ; call s%zero()
  !   allocate (t, source=b)     ; call t%zero()
        
  !   ! --> Compute initial residual: r = b - Ax.
  !   if (trans) then
  !      call A%rmatvec(x, r)
  !   else
  !      call A%matvec(x, r)
  !   endif
  !   call r%axpby(-1.0_wp, b, 1.0_wp)
    
  !   r_hat = r ; rho = r_hat%dot(r) ; p = r 
    
  !   bicgstab_loop: do i = 1, maxiter
  !      if (trans) then
  !         call A%rmatvec(p, v)
  !      else
  !         call A%matvec(p, v)
  !      endif
       
  !      alpha = rho/r_hat%dot(v)
       
  !      ! s = r - alpha * v
  !      s = r ; call s%axpby(1.0_wp, v, -alpha)
       
  !      ! t = A * s
  !      if (trans) then
  !         call A%rmatvec(s, t)
  !      else
  !         call A%matvec(s, t)
  !      endif
  !      omega = t%dot(s)/t%dot(t)
       
  !      ! x = x + s * omega + p * alpha
  !      call x%axpby(1.0_wp, s, omega) ; call x%axpby(1.0_wp, p, alpha)
       
  !      ! r = s - t * omega
  !      r = s ; call r%axpby(1.0_wp, t, -omega)
       
  !      res = r%norm()  
  !      if (verbose) then
  !         write (*, *) "INFO : BICGSTAB residual after ", (i), "iterations : ", res
  !      end if
  !      if (res < tol) exit bicgstab_loop
       
  !      rho_new = r_hat%dot(r)
  !      beta = (alpha/omega) * (rho_new/rho)
       
  !      ! s = p - v * omega ! reusing s vector
  !      s = p ; call s%axpby(1.0_wp, v, -omega)
       
  !      ! p = r + s * beta
  !      p = r ; call p%axpby(1.0_wp, s, beta)
       
  !      rho = rho_new
       
  !   end do bicgstab_loop
    
  !   deallocate (r, r_hat, p, v, s, t)
  !   return
  ! end subroutine bicgstab
  
  ! --> Utility Functions -----
  
  subroutine initialize_krylov_basis(X)
    ! Initializes Krylov basis vectors to zero, except for
    ! the first vector which is provided by the user
    class(abstract_vector), dimension(:), intent(inout) :: X
    integer :: i
    do i = 2, size(X)
       call X(i)%zero()
    end do
  end subroutine initialize_krylov_basis
  
end module IterativeSolvers
