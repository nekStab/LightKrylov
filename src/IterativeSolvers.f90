module lightkrylov_IterativeSolvers
  use lightkrylov_Utils
  use lightkrylov_AbstractVector
  use lightkrylov_LinearOperator
  use lightkrylov_BaseKrylov
  use stdlib_sorting, only : sort_index, int_size
  use stdlib_optval , only : optval
  use stdlib_io_npy , only : save_npy
  implicit none
  include "dtypes.h"

  private
  public :: eigs, eighs, gmres, save_eigenspectrum, svds, cg
  public :: two_sided_eigs
  public :: abstract_linear_solver

  !------------------------------------------------
  !-----                                      -----
  !-----     ABSTRACT PRECONDITIONER TYPE     -----
  !-----                                      -----
  !------------------------------------------------

  ! --> Abstract type definition.
  type, abstract, public :: abstract_preconditioner
   contains
     private
     procedure(abstract_apply_precond), deferred, pass(self), public :: apply
     procedure(abstract_undo_precond) , deferred, pass(self), public :: undo
  end type abstract_preconditioner

  ! --> Definition of the type-bound procedures interfaces.
  abstract interface
     !> Apply the preconditionner.
     subroutine abstract_apply_precond(self, vec_inout)
       import abstract_preconditioner, abstract_vector
       class(abstract_preconditioner), intent(in)    :: self
       class(abstract_vector)        , intent(inout) :: vec_inout
     end subroutine abstract_apply_precond

     !> Undo the action of the preconditioner.
     subroutine abstract_undo_precond(self, vec_inout)
       import abstract_preconditioner, abstract_vector
       class(abstract_preconditioner), intent(in)    :: self
       class(abstract_vector)        , intent(inout) :: vec_inout
     end subroutine abstract_undo_precond
  end interface

  !--------------------------------------------------------
  !-----                                              -----
  !-----     GENERIC INTERFACE FOR LINEAR SOLVERS     -----
  !-----                                              -----
  !--------------------------------------------------------

  abstract interface
     subroutine abstract_linear_solver(A, b, x, info, preconditioner, options, transpose)
       import abstract_linop, abstract_vector, abstract_opts, abstract_preconditioner
       !> Linear problem.
       class(abstract_linop) , intent(in) :: A
       class(abstract_vector), intent(in) :: b
       !> Solution vector.
       class(abstract_vector), intent(inout) :: x
       !> Information flag.
       integer               , intent(out) :: info
       !> Solver options.
       class(abstract_opts)  , optional, intent(in) :: options
       !> Preconditioner.
       class(abstract_preconditioner), optional, intent(in) :: preconditioner
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

  subroutine eigs(A, X, eigvecs, eigvals, residuals, info, nev, tolerance, verbosity, transpose)
    !> Linear Operator.
    class(abstract_linop), intent(in) :: A
    !> Krylov basis.
    class(abstract_vector), intent(inout) :: X(:)
    !> Coordinates of eigenvectors in Krylov basis, eigenvalues and associated residuals.
    complex(kind=wp), intent(out) :: eigvecs(:, :)
    complex(kind=wp), intent(out) :: eigvals(:)
    real(kind=wp)   , intent(out) :: residuals(:)
    !> Information flag.
    integer, intent(out) :: info
    !> Number of desired eigenvalues.
    integer, optional, intent(in) :: nev
    integer                       :: nev_, conv
    !> Tolerance control.
    real(kind=wp), optional, intent(in) :: tolerance
    real(kind=wp)                       :: tol
    !> Verbosity control.
    logical, optional, intent(in) :: verbosity
    logical :: verbose
    !> Transpose operator.
    logical, optional, intent(in) :: transpose
    logical :: trans

    !> Upper Hessenberg matrix.
    real(kind=wp) :: H(size(X), size(X)-1)
    real(kind=wp) :: beta
    !> Krylov subspace dimension.
    integer :: kdim
    !> Miscellaneous.
    integer :: i, j, k
    integer(int_size) :: indices(size(X)-1)
    real(kind=wp)     :: abs_eigvals(size(X)-1)
    real(kind=wp)     :: alpha

    ! --> Dimension of the Krylov subspace.
    kdim = size(X) - 1

    !> Shape assertion.
    call assert_shape(eigvecs, [kdim, kdim], "eigs", "eigvecs")

    ! --> Deals with the optional arguments.
    verbose = optval(verbosity, .false.)
    trans   = optval(transpose, .false.)
    nev_    = optval(nev, size(X)-1)
    tol     = optval(tolerance, rtol)

    ! --> Initialize variables.
    H = 0.0_wp ; residuals = 0.0_wp ; eigvals = (0.0_wp, 0.0_wp) ; eigvecs = (0.0_wp, 0.0_wp)
    !> Make sure the first Krylov vector has unit-norm.
    alpha = X(1)%norm() ; call X(1)%scal(1.0_wp / alpha)
    call initialize_krylov_subspace(X(2:kdim+1))

    arnoldi : do k = 1, kdim
       !> Arnoldi step.
       call arnoldi_factorization(A, X, H, info, kstart=k, kend=k, verbosity=verbose, transpose=trans)

       if (info < 0) then
          if (verbose) then
             write(*, *) "INFO : Arnoldi iteration failed. Exiting the eig subroutine."
             write(*, *) "       Arnoldi exit code :", info
          endif
          info = -1
          return
       endif

       !> Spectral decomposition of the k x k Hessenberg matrix.
       eigvals = (0.0_wp, 0.0_wp) ; eigvecs = (0.0_wp, 0.0_wp)
       call eig(H(1:k, 1:k), eigvecs(1:k, 1:k), eigvals(1:k))

       !> Sort eigenvalues.
       abs_eigvals = abs(eigvals) ; call sort_index(abs_eigvals, indices, reverse=.true.)
       eigvals = eigvals(indices) ; eigvecs = eigvecs(:, indices)

       !> Compute residuals.
       beta = H(k+1, k) !> Get Krylov residual vector norm.
       residuals(1:k) = compute_residual(beta, abs(eigvecs(k, 1:k)))

       !> Check convergence.
       conv = count(residuals(1:k) < tol)
       if (conv >= nev_) then
          if (verbose) then
             write(*, *) "INFO : The first ", conv, "eigenpairs have converged."
             write(*, *) "       Exiting the computation."
          endif
          exit arnoldi
       endif

    enddo arnoldi

    return
  end subroutine eigs

  subroutine eighs(A, X, eigvecs, eigvals, residuals, info, nev, tolerance, verbosity)
    !> Linear Operator.
    class(abstract_spd_linop), intent(in) :: A
    !> Krylov basis.
    class(abstract_vector), intent(inout) :: X(:)
    !> Coordinates of eigenvectors in Krylov basis, eigenvalues, and associated residuals
    real(kind=wp), intent(out) :: eigvecs(:, :)
    real(kind=wp), intent(out) :: eigvals(:)
    real(kind=wp), intent(out) :: residuals(:)
    !> Information flag.
    integer, intent(out) :: info
    !> Number of converged eigenvalues needed.
    integer, optional, intent(in) :: nev
    integer                       :: nev_, conv
    !> Tolerance control.
    real(kind=wp), optional, intent(in) :: tolerance
    real(kind=wp)                       :: tol
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
    real(kind=wp) :: alpha

    ! --> Dimension of the Krylov subspace.
    kdim = size(X) - 1

    ! --> Deals with the optional argument.
    verbose = optval(verbosity, .false.)
    nev_    = optval(nev, size(X)-1)
    tol     = optval(tolerance, rtol)

    ! --> Initialize all variables.
    T = 0.0_wp ; residuals = 0.0_wp ; eigvecs = 0.0_wp ; eigvals = 0.0_wp
    !> Make sure the first Krylov vector has unit-norm.
    alpha = X(1)%norm() ; call X(1)%scal(1.0_wp / alpha)
    call initialize_krylov_subspace(X(2:kdim+1))

    lanczos : do k = 1, kdim
       ! --> Compute Lanczos tridiagonalization.
       call lanczos_tridiagonalization(A, X, T, info, kstart=k, kend=k, verbosity=verbose)
       
       if (info < 0) then
          if (verbose) then
             write(*, *) "INFO : Lanczos iteration failed. Exiting the eigh subroutine."
             write(*, *) "       Lanczos exit code :", info
          endif
          info = -1
          return
       endif
       
       ! --> Compute spectral decomposition of the tridiagonal matrix.
       eigvals = 0.0_wp ; eigvecs = 0.0_wp
       call eigh(T(1:k, 1:k), eigvecs(1:k, 1:k), eigvals(1:k))

       ! --> Sort eigenvalues in decreasing order.
       call sort_index(eigvals, indices, reverse=.true.) ; eigvecs = eigvecs(:, indices)

       ! --> Compute the residual associated with each eigenpair.
       beta = T(k+1, k) !> Get Krylov residual vector norm.
       residuals(1:k) = compute_residual(beta, eigvecs(k, 1:k))

       !> Check convergence.
       conv = count(residuals(1:k) < tol)
       if (conv >= nev_) then
          if (verbose) then
             write(*, *) "INFO : The first", conv, "eigenpairs have converged."
             write(*, *)         "Exiting the computation."
          endif
          exit lanczos
       endif

    enddo lanczos

    return
  end subroutine eighs

  subroutine two_sided_eigs(A, V, W, rvecs, lvecs, eigvals, residuals, info, nev, tolerance, verbosity)
    !> Linear Operator.
    class(abstract_linop), intent(in) :: A
    !> Left and right Krylov basis.
    class(abstract_vector), dimension(:), intent(inout) :: V
    class(abstract_vector), dimension(:), intent(inout) :: W
    !> Coordinates of eigenvectors in Krylov bases, eigenvalues and associated residuals.
    complex(kind=wp), dimension(size(V)-1, size(V)-1), intent(out) :: rvecs
    complex(kind=wp), dimension(size(W)-1, size(W)-1), intent(out) :: lvecs
    complex(kind=wp), dimension(size(V)-1)           , intent(out) :: eigvals
    real(kind=wp)   , dimension(size(V)-1)           , intent(out) :: residuals
    !> Information flag.
    integer, intent(out) :: info
    !> Number of desired eigenvalues.
    integer, optional, intent(in) :: nev
    integer                       :: nev_, conv
    !> Tolerance control.
    real(kind=wp), optional, intent(in) :: tolerance
    real(kind=wp)                       :: tol
    !> Verbosity control.
    logical, optional, intent(in) :: verbosity
    logical                       :: verbose

    !> Tridiagonal matrix.
    real(kind=wp), dimension(size(V), size(V)) :: T
    real(kind=wp)                              :: beta, gamma
    !> Krylov subspace dimension.
    integer :: kdim
    !> Miscellaneous.
    integer :: i, j, k
    integer(int_size), dimension(size(V)-1) :: indices
    real(kind=wp)    , dimension(size(V)-1) :: abs_eigvals
    real(kind=wp)                           :: alpha, tmp

    !> Gram matrix.
    complex(kind=wp), dimension(size(V)-1, size(V)-1) :: G
    complex(kind=wp), dimension(size(V)-1, size(V)-1) :: H
    real(kind=wp) :: direct_residual, adjoint_residual

    !> Dimension of the Krylov subspaces.
    kdim = size(V)-1

    !> Deals with optional args.
    verbose = optval(verbosity, .false.)
    nev_    = optval(nev, kdim)
    tol     = optval(tolerance, rtol)

    !> Initialize variables.
    T = 0.0_wp ; residuals = 0.0_wp ; eigvals = cmplx(0.0_wp, 0.0_wp, kind=wp)
    lvecs = cmplx(0.0_wp, 0.0_wp, kind=wp) ; rvecs = cmplx(0.0_wp, 0.0_wp, kind=wp)

    !> Make sure the first Krylov vectors are bi-orthogonal.
    tmp = V(1)%dot(W(1)) ; beta = sqrt(abs(tmp)) ; gamma = sign(beta, tmp)
    call V(1)%scal(1.0_wp / beta) ; call W(1)%scal(1.0_wp / gamma)

    do i = 2, size(V)
       call V(i)%zero() ; call W(i)%zero()
    enddo

    lanczos : do k = 1, kdim
       !> Lanczos step.
       call nonsymmetric_lanczos_tridiagonalization(A, V, W, T, info, kstart=k, kend=k, verbosity=verbose)
       if (info < 0) then
          if (verbose) then
             write(*, *) "INFO : Lanczos iteration failed. Exiting the eig subroutine."
             write(*, *) "       Lanczos exit code :", info
          endif
          info = -1
          return
       endif

       !> Spectral decomposition of the k x k tridiagonal matrix.
       eigvals = cmplx(0.0_wp, 0.0_wp, kind=wp)
       rvecs = cmplx(0.0_wp, 0.0_wp, kind=wp) ; lvecs = cmplx(0.0_wp, 0.0_wp, kind=wp)
       call eig(T(1:k, 1:k), rvecs(1:k, 1:k), eigvals(1:k))
       lvecs(1:k, 1:k) = rvecs(1:k, 1:k) ; call inv(lvecs(1:k, 1:k)) ; lvecs(1:k, 1:k) = transpose(conjg(lvecs(1:k, 1:k)))

       !> Sort eigenvalues.
       abs_eigvals = abs(eigvals) ; call sort_index(abs_eigvals, indices, reverse=.true.)
       eigvals = eigvals(indices) ; rvecs = rvecs(:, indices) ; lvecs = lvecs(:, indices)

       !> Compute residuals.
       beta = abs(T(k+1, k)) ; gamma = abs(T(k, k+1)) ; alpha = V(k+1)%norm()
       residuals(1:k) = compute_residual(beta*alpha, abs(rvecs(k, 1:k)))

       G = 0.0_wp ; H = 0.0_wp
       do i = 1, k
          G(i, i) = V(i)%dot(V(i)) ; H(i, i) = W(i)%dot(W(i))
          do j = i+1, k
             G(i, j) = V(i)%dot(V(j))
             G(j, i) = G(i, j)

             H(i, j) = W(i)%dot(W(j))
             H(j, i) = H(i, j)
          enddo
       enddo
       residuals = 100.0_wp
       do i = 1, k
          alpha = V(k+1)%norm()
          direct_residual = compute_residual(beta*alpha, abs(rvecs(k, i))) &
               / sqrt( real(dot_product(rvecs(:, k), matmul(G, rvecs(:, k)))) )
          alpha = W(k+1)%norm()
          adjoint_residual = compute_residual(gamma*alpha, abs(lvecs(k, i))) &
               / sqrt( real(dot_product(lvecs(:, k), matmul(H, lvecs(:, k)))) )
          residuals(i) = max(direct_residual, adjoint_residual)
       enddo

       !> Check convergence.
       conv = count(residuals(1:k) < tol)
       write(*, *) "--- Iteration ", k, "---"
       write(*, *) conv, "eigentriplets have converged."
       write(*, *)
       if (conv >= nev_) then
          info = k
          exit lanczos
       endif

    end do lanczos

    call sort_index(residuals, indices, reverse=.false.)
    eigvals = eigvals(indices) ; rvecs = rvecs(:, indices) ; lvecs = lvecs(:, indices)
    return
  end subroutine two_sided_eigs

  !----------------------------------------------
  !-----                                    -----
  !-----     SINGULAR VALUE COMPUTATION     -----
  !-----                                    -----
  !----------------------------------------------

  subroutine svds(A, U, V, uvecs, vvecs, sigma, residuals, info, nev, tolerance, verbosity)
    !> Linear Operator.
    class(abstract_linop), intent(in) :: A
    !> Krylov bases.
    class(abstract_vector), intent(inout) :: U(:) ! Basis for left sing. vectors.
    class(abstract_vector), intent(inout) :: V(:) ! Basis for right sing. vectors.
    !> Coordinates of singular vectors in Krylov bases, singular values, and associated residuals.
    real(kind=wp), intent(out) :: uvecs(:, :), vvecs(:, :)
    real(kind=wp), intent(out) :: sigma(:)
    real(kind=wp), intent(out) :: residuals(:)
    !> Information flag.
    integer, intent(out) :: info
    !> Number of converged singular triplets.
    integer, optional, intent(in) :: nev
    integer                       :: nev_, conv
    !> Tolerance control.
    real(kind=wp), optional, intent(in) :: tolerance
    real(kind=wp)                       :: tol
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
    nev_    = optval(nev, size(U)-1)
    tol     = optval(tolerance, rtol)
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

    call assert_shape(uvecs, [kdim, kdim], "svds", "uvecs")
    call assert_shape(vvecs, [kdim, kdim], "svds", "vvecs")

    ! --> Initialize variables.
    B = 0.0_wp ; residuals = 0.0_wp ; uvecs = 0.0_wp ; vvecs = 0.0_wp ; sigma = 0.0_wp
    !> Make sure the first Krylov vector has unit-norm.
    beta = U(1)%norm() ; call U(1)%scal(1.0_wp / beta)
    call initialize_krylov_subspace(U(2:kdim+1))
    call initialize_krylov_subspace(V(2:kdim+1))

    lanczos : do k = 1, kdim
       ! --> Compute the Lanczos bidiagonalization.
       call lanczos_bidiagonalization(A, U, V, B, info, kstart=k, kend=k, verbosity=verbose)
       
       ! --> Compute the singular value decomposition of the bidiagonal matrix.
       sigma = 0.0_wp ; uvecs = 0.0_wp ; vvecs = 0.0_wp
       call svd(B(1:k, 1:k), uvecs(1:k, 1:k), sigma(1:k), vvecs(1:k, 1:k))
       
       ! --> Compute the residual associated with each singular triplet.
       beta = B(k+1, k) !> Get Krylov residual vector norm.
       residuals(1:k) = compute_residual(beta, vvecs(k, 1:k))

       !> Check convergence.
       conv = count(residuals(1:k) < tol)
       if (conv >= nev_) then
          if (verbose) then
             write(*, *) "INFO : The first ", conv, "singular triplets have converged."
             write(*, *) "       Exiting the computation."
          endif
          exit lanczos
       endif

    enddo lanczos

    info = k

    return
  end subroutine svds

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
  subroutine gmres(A, b, x, info, preconditioner, options, transpose)
    !> Linear problem.
    class(abstract_linop)  , intent(in)    :: A ! Linear Operator.
    class(abstract_vector) , intent(in)    :: b ! Right-hand side.
    !> Solution vector.
    class(abstract_vector) , intent(inout) :: x
    !> Information flag.
    integer                , intent(out)   :: info
    !> Preconditioner.
    class(abstract_preconditioner), optional, intent(in) :: preconditioner
    class(abstract_preconditioner), allocatable          :: precond
    logical                                              :: has_precond
    !> Optional arguments.
    class(abstract_opts), optional, intent(in) :: options
    type(gmres_opts)                           :: opts
    logical             , optional, intent(in) :: transpose
    logical                                    :: trans

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
    class(abstract_vector), allocatable :: dx, wrk

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

    if (present(preconditioner)) then
       allocate(precond, source=preconditioner)
       has_precond = .true.
    else
       has_precond = .false.
    endif

    ! --> Initialize Krylov subspace.
    allocate(wrk, source=b) ; call wrk%zero()
    allocate(V(1:k_dim+1), source=b)
    do i = 1, size(V)
       call V(i)%zero()
    enddo
    allocate(H(k_dim+1, k_dim)) ; H = 0.0_wp
    allocate(y(1:k_dim))        ; y = 0.0_wp
    allocate(e(1:k_dim+1))      ; e = 0.0_wp

    info = 0

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
       
       ! --> Arnoldi factorization.
       arnoldi : do k = 1, k_dim
          wrk = V(k) ; if (has_precond) call precond%apply(wrk)

          ! --> Matrix-vector product.
          if (trans) then
             call A%rmatvec(wrk, V(k+1))
          else
             call A%matvec(wrk, V(k+1))
          endif

          ! --> Gram-Schmid orthogonalization (twice is enough).
          do j = 1, k
             alpha = V(k+1)%dot(V(j)) ; call V(k+1)%axpby(1.0_wp, V(j), -alpha)
             H(j, k) = alpha
          enddo
          do j = 1, k
             alpha = V(k+1)%dot(V(j)) ; call V(k+1)%axpby(1.0_wp, V(j), -alpha)
             H(j, k) = H(j, k) + alpha
          enddo

          ! --> Update Hessenberg matrix and normalize new Krylov vector.
          H(k+1, k) = V(k+1)%norm()
          if (H(k+1, k) > tol) then
             call V(k+1)%scal(1.0_wp / H(k+1, k))
          endif

          ! --> Least-squares problem.
          call lstsq(H(1:k+1, 1:k), e(1:k+1), y(1:k))

          ! --> Compute residual.
          beta = norm2(e(1:k+1) - matmul(H(1:k+1, 1:k), y(1:k)))
          if (verbose) then
             write(*, *) "INFO : GMRES residual after ", info+1, "iteration : ", beta
          endif

          ! --> Current number of iterations performed.
          info = info + 1

          ! --> Check convergence.
          if (beta .lt. tol) then
             exit arnoldi
          endif
       enddo arnoldi

       ! --> Update solution.
       k = min(k, k_dim)
       if (allocated(dx) .eqv. .false.) allocate(dx, source=x) ; call dx%zero()
       call get_vec(dx, V(1:k), y(1:k)) ; if (has_precond) call precond%apply(dx)
       call x%add(dx)

       ! --> Recompute residual for sanity check.
       if (trans) then
          call A%rmatvec(x, V(1))
       else
          call A%matvec(x, V(1))
       endif

       call V(1)%sub(b) ; call V(1)%scal(-1.0_wp)

       ! --> Initialize new starting Krylov vector if needed.
       beta = V(1)%norm() ; call V(1)%scal(1.0_wp / beta)

       ! --> Exit GMRES if desired accuracy is reached.
       if (beta .lt. tol) exit gmres_iterations
    
    enddo gmres_iterations

    ! --> Convergence information.
    write(*, *) "INFO : GMRES converged with residual ", beta
    write(*, *) "       Computation required ", info, "iterations"

    return
  end subroutine gmres

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
  subroutine cg(A, b, x, info, preconditioner, options)
    !> Linear problem.
    class(abstract_spd_linop), intent(in) :: A ! Linear Operator.
    class(abstract_vector), intent(in) :: b ! Right-hand side.
    !> Solution vector.
    class(abstract_vector), intent(inout) :: x
    !> Information flag.
    integer, intent(out)   :: info
    !> Preconditioner.
    class(abstract_preconditioner), optional, intent(in) :: preconditioner
    class(abstract_preconditioner), allocatable          :: precond
    logical                                              :: has_precond
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

    if (present(preconditioner)) then
       write(*, *) "INFO: CG does not support preconditioning yet. Precond is thus ignored."
       write(*, *)
    endif
    
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


  !=======================================================================================
  ! Biconjugate Gradient Stabilized (BiCGSTAB) Solver Subroutine
  !=======================================================================================
  !
  ! Purpose:
  ! --------
  ! Implements the Biconjugate Gradient Stabilized (BiCGSTAB) algorithm for solving
  ! nonsymmetric and possibly ill-conditioned linear systems Ax = b.
  !
  ! Algorithmic Features:
  ! ----------------------
  ! - Extends the BiCG algorithm by stabilizing the iterations.
  ! - Utilizes two search directions and two residuals to improve stability.
  ! - Iteratively updates both the approximate solution and the residuals.
  !
  ! Advantages:
  ! -----------
  ! - Capable of addressing nonsymmetric and ill-conditioned matrices.
  ! - Generally more stable and faster compared to the basic BiCG.
  ! - Suitable for large and sparse matrices.
  !
  ! Limitations:
  ! ------------
  ! - May experience stagnation for certain types of problems.
  ! - No preconditioning capabilities in the current implementation.
  !
  ! Input/Output Parameters:
  ! ------------------------
  ! - A        : Linear Operator (abstract_linop) [Input]
  ! - b        : Right-hand side (abstract_vector) [Input]
  ! - x        : Initial/Updated solution (abstract_vector) [Input/Output]
  ! - info     : Iteration Information flag (Integer) [Output]
  ! - maxiter  : Maximum number of iterations (Integer) [Optional, Input]
  ! - tol      : Convergence tolerance (real(kind=dp)) [Optional, Input]
  ! - verbosity: Verbosity control flag (Logical) [Optional, Input]
  !
  ! References:
  ! -----------
  ! - van der Vorst, H. A. (1992). "Bi-CGSTAB: A Fast and Smoothly Converging Variant of Bi-CG for the Solution of Nonsymmetric Linear Systems,"
  !   SIAM Journal on Scientific and Statistical Computing, 13(2), 631–644.
  !
  !=======================================================================================
  !subroutine bicgstab(A, b, x, info, preconditioner, options, transpose)
  !  !> Linear problem and initial guess.
  !  class(abstract_linop), intent(in) :: A
  !  class(abstract_vector), intent(in) :: b
  !  class(abstract_vector), intent(inout) :: x
  !  !> Output and optional input parameters.
  !  integer, intent(out) :: info
  !  class(abstract_preconditioner), optional, intent(in) :: preconditioner
  !  class(abstract_preconditioner), allocatable          :: precond
  !  logical                                              :: has_precond
  !  class(abstract_opts), optional, intent(in) :: options
  !  type(bicgstab_opts)                        :: opts
  !  logical, optional, intent(in) :: transpose
    
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

  !  if (present(preconditioner)) then
  !    write(*, *) "INFO: BICGSTAB does not support preconditioning yet. Precond is thus ignored."
  !    write(*, *)
  !  endif
  
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
  
end module lightkrylov_IterativeSolvers
