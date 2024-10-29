module LightKrylov_utils
    !!  This module provides a set of utility functions used throughout `LightKrylov`.
    !!  It includes:
    !!
    !!  - `assert_shape`: Assert that the shape of the argument is the expected shape.
    !!  - `eig`: Compute the eigenvalue decomposition of a general matrix.
    !!  - `sqrtm`: Compute the non-negative square root of a symmetric positive definite matrix using its SVD.
    !!  - `sqrtm_eig`: Compute the non-negative square root of a symmetric positive definite matrix using its eigenvalue decomposition.
    !!  - `schur`: Compute the Schur factorization of a general square matrix.
    !!  - `ordschur`: Re-order the Schur factorization to have the selected eigenvalues in the upper left block.
    !!
    !!  Note that as the development of `stdlib` progresses, some of these functions
    !!  will be deprecated in favor of the `stdlib` implementations.

    !--------------------------------------------
    !-----     Standard Fortran Library     -----
    !--------------------------------------------
    use iso_fortran_env, only: output_unit
    use stdlib_linalg, only: is_hermitian, is_symmetric, diag, svd, eigh
    ! Eigenvalue problem.
    use stdlib_linalg_lapack, only: geev
    ! Schur factorization.
    use stdlib_linalg_lapack, only: gees, trsen

    !-------------------------------
    !-----     LightKrylov     -----
    !-------------------------------
    ! Various constants.
    use LightKrylov_Logger
    use LightKrylov_Constants

    implicit none
    private

    character(len=128), parameter :: this_module = 'LightKrylov_Utils'

    public :: assert_shape, norml, log2
    ! Compute AX = XD for general dense matrices.
    public :: eig
    ! Compute matrix sqrt of input symmetric/hermitian positive definite matrix A
    public :: sqrtm
    public :: sqrtm_eig
    ! Compute AX = XS where S is in Schur form.
    public :: schur
    ! Re-orders the Schur factorization of A.
    public :: ordschur

    interface assert_shape
        !! This interface provides methods to assert that the shape of its input vector
        !! or matrix is the expected shape. It throws an error if not.
        module procedure assert_shape_vector_rsp
        module procedure assert_shape_matrix_rsp
        module procedure assert_shape_vector_rdp
        module procedure assert_shape_matrix_rdp
        module procedure assert_shape_vector_csp
        module procedure assert_shape_matrix_csp
        module procedure assert_shape_vector_cdp
        module procedure assert_shape_matrix_cdp
    end interface

    interface norml
        !! This interface provides methods to compute the infinity norm of a matrix.
        !! Note that it'll eventually be superseeded by the `stdlib` implementation.
        module procedure norml_rsp
        module procedure norml_rdp
        module procedure norml_csp
        module procedure norml_cdp
    end interface

    interface log2
        !! Utility function to compute the base-2 logarithm of a real number.
        module procedure log2_rsp
        module procedure log2_rdp
    end interface

    interface eig
        !!  Computes the eigenvalue decomposition of a general square matrix.
        !!
        !!  ### Description
        !!
        !!  This interface provides methods to compute the solution to the eigenproblem
        !!  \( \mathbf{Ax} = \lambda \mathbf{x} \), where $\mathbf{A}$ is a square `real`
        !!  or `complex` matrix.
        !!
        !!  Result array `lambda` returns the eigenvalues of \( \mathbf{A} \), while `vecs`
        !!  returns the corresponding eigenvectors. Note that it follows the LAPACK convention
        !!  when \( \mathbf{A} \) is `real`. The solver is based on LAPACK's `*GEEV` backends.
        !!
        !!  ### Syntax
        !!
        !!  `call eig(A, vecs, lambda)`
        !!
        !!  ### Arguments
        !!
        !!  `A`: `real` or `complex` square array containing the coefficient matrix. It is an `intent(in)` argument.
        !!
        !!  `vecs`: Square array of the same size, type, and kind as `A` containing the eigenvectors
        !!  (following LAPACK's convention for `real` matrices). It is an `intent(out)` argument.
        !!
        !!  `lambda`: `complex` rank-1 array of the same kind as `A` containing the eigenvalues.
        !!  It is an `intent(out)` argument.
        !!
        !!  @note
        !!  Due to the abstrct nature of the vector types defined in `LightKrylov`, it is unlikely
        !!  that this implementation will be superseeded in favor of the `stdlib` one as the latter
        !!  does not follow the LAPACK's convention.
        !!  @endnote
        module procedure eig_rsp
        module procedure eig_rdp
        module procedure eig_csp
        module procedure eig_cdp
    end interface

    interface schur
        !!  Computes the Schur factorization of a general square matrix.
        !!
        !!  ### Description
        !!
        !!  This interface provides methods to compute the Schur factorization of a `real` or
        !!  `complex` square matrix. Note that, if \( \mathbf{A} \) is `real`, it returns the
        !!  real Schur form.
        !!
        !!  Result array `eigvals` returns the eigenvalues of \( \mathbf{A} \) while `Z`
        !!  contains the Schur basis.
        !!
        !!  ### Syntax
        !!
        !!  `call schur(A, Z, eigvals)`
        !!
        !!  ### Arguments
        !!
        !!  `A`: `real` or `complex` square array containing the coefficient matrix. On exit, it
        !!  is overwritten with its (real) Schur factorization. It is an `intent(inout)` argument.
        !!  
        !!  `Z`: Two-dimensional square array of the same size, type and kind as `A`. It contains
        !!  the Schur basis. It is an `intent(out)` argument.
        !!
        !!  `eigvals`: `complex` rank-1 array of the same kind as `A` containing the eigenvalues.
        !!  It is an `intent(out)` arguement.
        module procedure schur_rsp
        module procedure schur_rdp
        module procedure schur_csp
        module procedure schur_cdp
    end interface

    interface ordschur
        !!  Given the Schur factorization and basis of a matrix, reorders it to have the selected
        !!  eigenvalues in the upper left block.
        !!
        !!  ### Description
        !!
        !!  This interface provides methods to re-order the Schur factorization of a `real` or
        !!  `complex` square matrix. Note that, if \( \mathbf{A} \) is `real`, it returns the
        !!  real Schur form.
        !!
        !!  ### Syntax
        !!
        !!  `call ordschur(T, Q, selected)`
        !!
        !!  ### Arguments
        !!
        !!  `T`: `real` or `complex` square array containing the Schur factorization of a matrix. 
        !!  On exit, it is overwritten with its re-ordered counterpart. It is an `intent(inout)` argument.
        !!  
        !!  `Q`: Two-dimensional square array of the same size, type and kind as `A`. It contains
        !!  the original Schur basis on entry and the re-ordered one on exit.
        !!  It is an `intent(inout)` argument.
        !!
        !!  `selected`: `logical` rank-1 array selecting which eigenvalues need to be moved in the
        !!  upper left block of the Schur factorization.
        !!  It is an `intent(in)` arguement.
        module procedure ordschur_rsp
        module procedure ordschur_rdp
        module procedure ordschur_csp
        module procedure ordschur_cdp
    end interface

    interface sqrtm
        !!  Computes the non-negative square root of a symmetric positive definite matrix
        !!  using its singular value decomposition.
        !!
        !!  ### Description
        !!
        !!  This interface provides methods to compute the non-negative square root of a symmetric
        !!  (hermitian) positive definite matrix \( \mathbf{A} \).
        !!
        !!  ### Syntax
        !!
        !!  `call sqrtm(A, sqrtmA, info)`
        !!
        !!  ### Arguments
        !!  
        !!  `A`: Symmetric (hermitian) positive definite matrix whose non-negative square root
        !!  needs to be computed. It is an `intent(in)` argument.
        !!
        !!  `sqrtmA`: Non-negative square root of `A`. It has the same size, kind and type as `A`.
        !!  It is an `intent(out)` argument.
        !!
        !!  `info`: Information flag. It is an `intent(out)` argument.
        module procedure sqrtm_rsp
        module procedure sqrtm_rdp
        module procedure sqrtm_csp
        module procedure sqrtm_cdp
    end interface

    interface sqrtm_eig
        !!  Computes the non-negative square root of a symmetric positive definite matrix
        !!  using its eigenvalue decomposition.
        !!
        !!  ### Description
        !!
        !!  This interface provides methods to compute the non-negative square root of a symmetric
        !!  (hermitian) positive definite matrix \( \mathbf{A} \).
        !!
        !!  ### Syntax
        !!
        !!  `call sqrtm_eig(A, sqrtmA, info)`
        !!
        !!  ### Arguments
        !!  
        !!  `A`: Symmetric (hermitian) positive definite matrix whose non-negative square root
        !!  needs to be computed. It is an `intent(in)` argument.
        !!
        !!  `sqrtmA`: Non-negative square root of `A`. It has the same size, kind and type as `A`.
        !!  It is an `intent(out)` argument.
        !!
        !!  `info`: Information flag. It is an `intent(out)` argument.
        module procedure sqrtm_eig_rsp
        module procedure sqrtm_eig_rdp
        module procedure sqrtm_eig_csp
        module procedure sqrtm_eig_cdp
    end interface

    !------------------------------------------------
    !-----     OPTS TYPE FOR LINEAR SOLVERS     -----
    !------------------------------------------------

    type, abstract, public :: abstract_opts
        !! Abstract type container for options from which all others are being extended.
    end type

    type, extends(abstract_opts), public :: gmres_sp_opts
        !! GMRES options.
        integer :: kdim = 30
        !! Dimension of the Krylov subspace (default: 30).
        integer :: maxiter = 10
        !! Maximum number of `gmres` restarts (default: 10).
    end type

    type, extends(abstract_opts), public :: cg_sp_opts
        !! Conjugate gradient options.
        integer :: maxiter = 100
        !! Maximum number of `cg` iterations (default: 100).
    end type

    type, extends(abstract_opts), public :: newton_sp_opts
        !! Options for Newton-Krylov fixed-point iteration.
        integer :: maxiter = 100
        !! Maximum number of Newton iterations (default = 100)
        logical :: ifbisect = .false.
        !! Bisection toggle to enforce residual reduction (default = .false.)
        integer :: maxstep_bisection = 5
        !! Maximum number of bisections (evaluations of F) for step selection (default = 5)
        !! Ignored if ifbisect = .false.
    end type
    
    type, extends(abstract_opts), public :: gmres_dp_opts
        !! GMRES options.
        integer :: kdim = 30
        !! Dimension of the Krylov subspace (default: 30).
        integer :: maxiter = 10
        !! Maximum number of `gmres` restarts (default: 10).
    end type

    type, extends(abstract_opts), public :: cg_dp_opts
        !! Conjugate gradient options.
        integer :: maxiter = 100
        !! Maximum number of `cg` iterations (default: 100).
    end type

    type, extends(abstract_opts), public :: newton_dp_opts
        !! Options for Newton-Krylov fixed-point iteration.
        integer :: maxiter = 100
        !! Maximum number of Newton iterations (default = 100)
        logical :: ifbisect = .false.
        !! Bisection toggle to enforce residual reduction (default = .false.)
        integer :: maxstep_bisection = 5
        !! Maximum number of bisections (evaluations of F) for step selection (default = 5)
        !! Ignored if ifbisect = .false.
    end type
    

contains

    !-------------------------------------
    !-----     VARIOUS UTILITIES     -----
    !-------------------------------------

    subroutine assert_shape_vector_rsp(v, size, vecname, module, procedure)
        !! Utility function to assert the shape of a vector.
        real(sp), intent(in) :: v(:)
        !! Vector whose dimension need to be asserted.
        integer, intent(in) :: size(:)
        !! Expected dimensions of v.
        character(len=*), intent(in) :: vecname
        !! Name of the asserted vector.
        character(len=*), intent(in) :: module
        !! Name of the module where assertion is done.
        character(len=*), intent(in) :: procedure
        !! Name of the routine where assertion is done.

        if(any(shape(v) /= size)) then
            print *, "Vector "//vecname//" has illegal length ", shape(v), &
                           & ". Expected length is ", size, ". Aborting due to illegal vector length."
            call stop_error('Vector length assertion error', module=module, procedure=procedure)
        endif
        return
    end subroutine assert_shape_vector_rsp

    subroutine assert_shape_matrix_rsp(A, size, matname, module, procedure)
        !! Utility function to assert the shape of a matrix.
        real(sp), intent(in) :: A(:, :)
        !! Matrix whose dimension need to be asserted.
        integer, intent(in) :: size(:)
        !! Expected dimensions of A.
        character(len=*), intent(in) :: matname
        !! Name of the asserted matrix.
        character(len=*), intent(in) :: module
        !! Name of the module where assertion is done.
        character(len=*), intent(in) :: procedure
        !! Name of the routine where assertion is done.

        if(any(shape(A) /= size)) then
            print *, "Matrix "//matname//" has illegal shape ", shape(A), &
                        & ". Expected shape is ", size, ". Aborting due to illegal vector length."
            call stop_error('Matrix shape assertion error', module=module, procedure=procedure)
        endif
        return
    end subroutine assert_shape_matrix_rsp

    subroutine assert_shape_vector_rdp(v, size, vecname, module, procedure)
        !! Utility function to assert the shape of a vector.
        real(dp), intent(in) :: v(:)
        !! Vector whose dimension need to be asserted.
        integer, intent(in) :: size(:)
        !! Expected dimensions of v.
        character(len=*), intent(in) :: vecname
        !! Name of the asserted vector.
        character(len=*), intent(in) :: module
        !! Name of the module where assertion is done.
        character(len=*), intent(in) :: procedure
        !! Name of the routine where assertion is done.

        if(any(shape(v) /= size)) then
            print *, "Vector "//vecname//" has illegal length ", shape(v), &
                           & ". Expected length is ", size, ". Aborting due to illegal vector length."
            call stop_error('Vector length assertion error', module=module, procedure=procedure)
        endif
        return
    end subroutine assert_shape_vector_rdp

    subroutine assert_shape_matrix_rdp(A, size, matname, module, procedure)
        !! Utility function to assert the shape of a matrix.
        real(dp), intent(in) :: A(:, :)
        !! Matrix whose dimension need to be asserted.
        integer, intent(in) :: size(:)
        !! Expected dimensions of A.
        character(len=*), intent(in) :: matname
        !! Name of the asserted matrix.
        character(len=*), intent(in) :: module
        !! Name of the module where assertion is done.
        character(len=*), intent(in) :: procedure
        !! Name of the routine where assertion is done.

        if(any(shape(A) /= size)) then
            print *, "Matrix "//matname//" has illegal shape ", shape(A), &
                        & ". Expected shape is ", size, ". Aborting due to illegal vector length."
            call stop_error('Matrix shape assertion error', module=module, procedure=procedure)
        endif
        return
    end subroutine assert_shape_matrix_rdp

    subroutine assert_shape_vector_csp(v, size, vecname, module, procedure)
        !! Utility function to assert the shape of a vector.
        complex(sp), intent(in) :: v(:)
        !! Vector whose dimension need to be asserted.
        integer, intent(in) :: size(:)
        !! Expected dimensions of v.
        character(len=*), intent(in) :: vecname
        !! Name of the asserted vector.
        character(len=*), intent(in) :: module
        !! Name of the module where assertion is done.
        character(len=*), intent(in) :: procedure
        !! Name of the routine where assertion is done.

        if(any(shape(v) /= size)) then
            print *, "Vector "//vecname//" has illegal length ", shape(v), &
                           & ". Expected length is ", size, ". Aborting due to illegal vector length."
            call stop_error('Vector length assertion error', module=module, procedure=procedure)
        endif
        return
    end subroutine assert_shape_vector_csp

    subroutine assert_shape_matrix_csp(A, size, matname, module, procedure)
        !! Utility function to assert the shape of a matrix.
        complex(sp), intent(in) :: A(:, :)
        !! Matrix whose dimension need to be asserted.
        integer, intent(in) :: size(:)
        !! Expected dimensions of A.
        character(len=*), intent(in) :: matname
        !! Name of the asserted matrix.
        character(len=*), intent(in) :: module
        !! Name of the module where assertion is done.
        character(len=*), intent(in) :: procedure
        !! Name of the routine where assertion is done.

        if(any(shape(A) /= size)) then
            print *, "Matrix "//matname//" has illegal shape ", shape(A), &
                        & ". Expected shape is ", size, ". Aborting due to illegal vector length."
            call stop_error('Matrix shape assertion error', module=module, procedure=procedure)
        endif
        return
    end subroutine assert_shape_matrix_csp

    subroutine assert_shape_vector_cdp(v, size, vecname, module, procedure)
        !! Utility function to assert the shape of a vector.
        complex(dp), intent(in) :: v(:)
        !! Vector whose dimension need to be asserted.
        integer, intent(in) :: size(:)
        !! Expected dimensions of v.
        character(len=*), intent(in) :: vecname
        !! Name of the asserted vector.
        character(len=*), intent(in) :: module
        !! Name of the module where assertion is done.
        character(len=*), intent(in) :: procedure
        !! Name of the routine where assertion is done.

        if(any(shape(v) /= size)) then
            print *, "Vector "//vecname//" has illegal length ", shape(v), &
                           & ". Expected length is ", size, ". Aborting due to illegal vector length."
            call stop_error('Vector length assertion error', module=module, procedure=procedure)
        endif
        return
    end subroutine assert_shape_vector_cdp

    subroutine assert_shape_matrix_cdp(A, size, matname, module, procedure)
        !! Utility function to assert the shape of a matrix.
        complex(dp), intent(in) :: A(:, :)
        !! Matrix whose dimension need to be asserted.
        integer, intent(in) :: size(:)
        !! Expected dimensions of A.
        character(len=*), intent(in) :: matname
        !! Name of the asserted matrix.
        character(len=*), intent(in) :: module
        !! Name of the module where assertion is done.
        character(len=*), intent(in) :: procedure
        !! Name of the routine where assertion is done.

        if(any(shape(A) /= size)) then
            print *, "Matrix "//matname//" has illegal shape ", shape(A), &
                        & ". Expected shape is ", size, ". Aborting due to illegal vector length."
            call stop_error('Matrix shape assertion error', module=module, procedure=procedure)
        endif
        return
    end subroutine assert_shape_matrix_cdp


    !-------------------------------------------
    !-----     LAPACK MATRIX INVERSION     -----
    !-------------------------------------------

    subroutine eig_rsp(A, vecs, vals)
        !! Eigenvalue decomposition of a dense matrix using LAPACK.
        real(sp), intent(in) :: A(:, :)
        !! Matrix to be factorized.
        real(sp), intent(out) :: vecs(:, :)
        !! Eigenvectors.
        complex(sp), intent(out) :: vals(:)

        ! Internal variables
        character :: jobvl = "n", jobvr = "v"
        integer :: n, lwork, info, lda, ldvl, ldvr
        real(sp) :: A_tilde(size(A, 1), size(A, 2)), vl(1, size(A, 2))
        real(sp) :: work(4*size(A, 1)), wr(size(A, 1)), wi(size(A, 1))

        ! Setup variables.
        n = size(A, 1) ; lda = n ; ldvl = 1 ; ldvr = n ; a_tilde = a
        lwork = 4*n

        ! Eigendecomposition.
        call geev(jobvl, jobvr, n, a_tilde, lda, wr, wi, vl, ldvl, vecs, ldvr, work, lwork, info)
        call check_info(info, 'GEEV', module=this_module, procedure='eig_rsp')

        ! Reconstruct eigenvalues
        vals = one_csp*wr + one_im_csp*wi

        return
    end subroutine eig_rsp

    subroutine schur_rsp(A, Z, eigvals)
        !! Compute the Schur form (in-place) and Schur vectors of the matrix `A`.
        real(sp), intent(inout) :: A(:, :)
        !! Matrix to be factorized.
        real(sp), intent(out) :: Z(:, :)
        !! Schur basis.
        complex(sp), intent(out) :: eigvals(:)
        !! Eigenvalues.

        ! Internal variables.
        character :: jobvs = "v", sort = "n"
        integer :: n, lda, sdim, ldvs, lwork, info
        logical, allocatable :: bwork(:)
        real(sp), allocatable :: work(:)
        real(sp), allocatable :: wr(:), wi(:)

        ! Allocate variables.
        n = size(A, 1) ; lda = n ; ldvs = n ; lwork =  3*n 
        allocate(bwork(n)) ; allocate(work(lwork)) ; 

        allocate(wr(size(eigvals)), wi(size(eigvals)))
        call gees(jobvs, sort, dummy_select, n, A, lda, sdim, wr, wi, Z, ldvs, work, lwork, bwork, info)
        call check_info(info, 'GEES', module=this_module, procedure='schur_rsp')

        ! Reconstruct eigenvalues
        eigvals = cmplx(wr, wi, kind=sp)

        return
    contains
        pure function dummy_select(wre, wim) result(out)
            real(sp), intent(in) :: wre
            real(sp), intent(in) :: wim
            logical :: out
            out = .false.
            return
        end function
    end subroutine schur_rsp

    subroutine ordschur_rsp(T, Q, selected)
        !! Re-order the Schur factorization from `schur` such that the selected eigenvalues
        !! are in the upper-left block.
        real(sp), intent(inout) :: T(:, :)
        !! Schur matrix to be re-ordered.
        real(sp), intent(inout) :: Q(:, :)
        !! Schur vectors to be re-ordered.
        logical, intent(in) :: selected(:)
        !! Boolean array defining the selected eigenvalues.

        ! Internal variables
        character :: job="n", compq="v"
        integer info, ldq, ldt, lwork, m, n
        real(sp) :: s, sep
        integer :: iwork(size(T, 1)), liwork
        real(sp) :: wi(size(T, 1)), wr(size(T, 1)), work(size(T, 1))

        ! Setup variables.
        n = size(T, 2) ; ldt = n ; ldq = n ; lwork = max(1, n)

        liwork = 1
        call trsen(job, compq, selected, n, T, ldt, Q, ldq, wr, wi, m, s, sep, work, lwork, iwork, liwork, info)
        call check_info(info, 'TRSEN', module=this_module, procedure='ordschur_rsp')

        return
    end subroutine ordschur_rsp

    subroutine sqrtm_rsp(X, sqrtmX, info)
      !! Matrix-valued sqrt function for dense symmetric/hermitian positive (semi-)definite matrices
      real(sp), intent(inout) :: X(:,:)
      !! Matrix of which to compute the sqrt
      real(sp), intent(out)   :: sqrtmX(size(X,1),size(X,1))
      !! Return matrix
      integer, intent(out) :: info
      !! Information flag

      ! internal
      real(sp) :: S(size(X,1))
      real(sp) :: U(size(X,1), size(X,1)), VT(size(X,1), size(X,1))
      integer :: i
      real(sp) :: symmetry_error
      character(len=256) :: msg

      info = 0

      ! Check if the matrix is symmetric
      symmetry_error = 0.5_sp*maxval(X - transpose(X))
      if (symmetry_error > rtol_sp) then
        write(msg,'(2(A,E9.2))') "Input matrix is not symmetric. 0.5*max(X-X.T) = ", &
            & symmetry_error, ", tol = ", rtol_sp
        call stop_error(msg, module=this_module, procedure='sqrtm_rsp')
      else if (symmetry_error > 10*atol_sp) then
        write(msg,'(A,E9.2)') "Input matrix is not exactly symmetric. 0.5*max(X-X.T) = ", symmetry_error
        call logger%log_warning(msg, module=this_module, procedure='sqrtm_rsp')
      end if

      ! Perform svd
      call svd(X, S, U, VT)

      ! Check if the matrix is positive definite (up to tol)
      do i = 1, size(S)
         if (S(i) .gt. 10*atol_sp ) then
            S(i) = sqrt(S(i))
         else
            S(i) = zero_rsp
            info = 1
         end if
      end do

      ! Reconstruct the square root matrix
      sqrtmX = matmul(U, matmul(diag(S), VT))

      return
    end subroutine

    subroutine sqrtm_eig_rsp(X, sqrtmX, info)
      !! Matrix-valued sqrt function for dense symmetric/hermitian positive (semi-)definite matrices
      real(sp), intent(in)  :: X(:,:)
      !! Matrix of which to compute the sqrt
      real(sp), intent(out) :: sqrtmX(size(X,1),size(X,1))
      !! Return matrix
      integer, intent(out) :: info
      !! Information flag

      ! internals
      real(sp) :: Xtmp(size(X, 1), size(X, 1))
      real(sp) :: lambda(size(X,1))
      real(sp) :: V(size(X,1), size(X,1))
      integer :: i
      character(len=256) :: msg

      info = 0

      ! Check if the matrix is symmetric
      if (.not. is_symmetric(X)) then
        write(msg,'(A)') "Input matrix is not symmetric."
        call stop_error(msg, module=this_module, procedure='sqrtm_rsp')
      end if

      ! Perform eigenvalue decomposition
      Xtmp = X ; call eigh(Xtmp, lambda, vectors=V, overwrite_a=.true.)

      ! Check if the matrix is positive definite (up to tol)
      do i = 1, size(lambda)
         if (abs(lambda(i)) .gt. 10*atol_sp ) then
            if (lambda(i) .gt. zero_rsp) then
               lambda(i) = sqrt(lambda(i))
            else
               lambda(i) = zero_rsp
               info = -1
            end if
         else
            lambda(i) = zero_rsp
            info = 1
         end if
      end do

      ! Reconstruct the square root matrix
      sqrtmX = matmul(V, matmul(diag(lambda), transpose(V)))

      return
    end subroutine
    subroutine eig_rdp(A, vecs, vals)
        !! Eigenvalue decomposition of a dense matrix using LAPACK.
        real(dp), intent(in) :: A(:, :)
        !! Matrix to be factorized.
        real(dp), intent(out) :: vecs(:, :)
        !! Eigenvectors.
        complex(dp), intent(out) :: vals(:)

        ! Internal variables
        character :: jobvl = "n", jobvr = "v"
        integer :: n, lwork, info, lda, ldvl, ldvr
        real(dp) :: A_tilde(size(A, 1), size(A, 2)), vl(1, size(A, 2))
        real(dp) :: work(4*size(A, 1)), wr(size(A, 1)), wi(size(A, 1))

        ! Setup variables.
        n = size(A, 1) ; lda = n ; ldvl = 1 ; ldvr = n ; a_tilde = a
        lwork = 4*n

        ! Eigendecomposition.
        call geev(jobvl, jobvr, n, a_tilde, lda, wr, wi, vl, ldvl, vecs, ldvr, work, lwork, info)
        call check_info(info, 'GEEV', module=this_module, procedure='eig_rdp')

        ! Reconstruct eigenvalues
        vals = one_cdp*wr + one_im_cdp*wi

        return
    end subroutine eig_rdp

    subroutine schur_rdp(A, Z, eigvals)
        !! Compute the Schur form (in-place) and Schur vectors of the matrix `A`.
        real(dp), intent(inout) :: A(:, :)
        !! Matrix to be factorized.
        real(dp), intent(out) :: Z(:, :)
        !! Schur basis.
        complex(dp), intent(out) :: eigvals(:)
        !! Eigenvalues.

        ! Internal variables.
        character :: jobvs = "v", sort = "n"
        integer :: n, lda, sdim, ldvs, lwork, info
        logical, allocatable :: bwork(:)
        real(dp), allocatable :: work(:)
        real(dp), allocatable :: wr(:), wi(:)

        ! Allocate variables.
        n = size(A, 1) ; lda = n ; ldvs = n ; lwork =  3*n 
        allocate(bwork(n)) ; allocate(work(lwork)) ; 

        allocate(wr(size(eigvals)), wi(size(eigvals)))
        call gees(jobvs, sort, dummy_select, n, A, lda, sdim, wr, wi, Z, ldvs, work, lwork, bwork, info)
        call check_info(info, 'GEES', module=this_module, procedure='schur_rdp')

        ! Reconstruct eigenvalues
        eigvals = cmplx(wr, wi, kind=dp)

        return
    contains
        pure function dummy_select(wre, wim) result(out)
            real(dp), intent(in) :: wre
            real(dp), intent(in) :: wim
            logical :: out
            out = .false.
            return
        end function
    end subroutine schur_rdp

    subroutine ordschur_rdp(T, Q, selected)
        !! Re-order the Schur factorization from `schur` such that the selected eigenvalues
        !! are in the upper-left block.
        real(dp), intent(inout) :: T(:, :)
        !! Schur matrix to be re-ordered.
        real(dp), intent(inout) :: Q(:, :)
        !! Schur vectors to be re-ordered.
        logical, intent(in) :: selected(:)
        !! Boolean array defining the selected eigenvalues.

        ! Internal variables
        character :: job="n", compq="v"
        integer info, ldq, ldt, lwork, m, n
        real(dp) :: s, sep
        integer :: iwork(size(T, 1)), liwork
        real(dp) :: wi(size(T, 1)), wr(size(T, 1)), work(size(T, 1))

        ! Setup variables.
        n = size(T, 2) ; ldt = n ; ldq = n ; lwork = max(1, n)

        liwork = 1
        call trsen(job, compq, selected, n, T, ldt, Q, ldq, wr, wi, m, s, sep, work, lwork, iwork, liwork, info)
        call check_info(info, 'TRSEN', module=this_module, procedure='ordschur_rdp')

        return
    end subroutine ordschur_rdp

    subroutine sqrtm_rdp(X, sqrtmX, info)
      !! Matrix-valued sqrt function for dense symmetric/hermitian positive (semi-)definite matrices
      real(dp), intent(inout) :: X(:,:)
      !! Matrix of which to compute the sqrt
      real(dp), intent(out)   :: sqrtmX(size(X,1),size(X,1))
      !! Return matrix
      integer, intent(out) :: info
      !! Information flag

      ! internal
      real(dp) :: S(size(X,1))
      real(dp) :: U(size(X,1), size(X,1)), VT(size(X,1), size(X,1))
      integer :: i
      real(dp) :: symmetry_error
      character(len=256) :: msg

      info = 0

      ! Check if the matrix is symmetric
      symmetry_error = 0.5_dp*maxval(X - transpose(X))
      if (symmetry_error > rtol_dp) then
        write(msg,'(2(A,E9.2))') "Input matrix is not symmetric. 0.5*max(X-X.T) = ", &
            & symmetry_error, ", tol = ", rtol_dp
        call stop_error(msg, module=this_module, procedure='sqrtm_rdp')
      else if (symmetry_error > 10*atol_dp) then
        write(msg,'(A,E9.2)') "Input matrix is not exactly symmetric. 0.5*max(X-X.T) = ", symmetry_error
        call logger%log_warning(msg, module=this_module, procedure='sqrtm_rdp')
      end if

      ! Perform svd
      call svd(X, S, U, VT)

      ! Check if the matrix is positive definite (up to tol)
      do i = 1, size(S)
         if (S(i) .gt. 10*atol_dp ) then
            S(i) = sqrt(S(i))
         else
            S(i) = zero_rdp
            info = 1
         end if
      end do

      ! Reconstruct the square root matrix
      sqrtmX = matmul(U, matmul(diag(S), VT))

      return
    end subroutine

    subroutine sqrtm_eig_rdp(X, sqrtmX, info)
      !! Matrix-valued sqrt function for dense symmetric/hermitian positive (semi-)definite matrices
      real(dp), intent(in)  :: X(:,:)
      !! Matrix of which to compute the sqrt
      real(dp), intent(out) :: sqrtmX(size(X,1),size(X,1))
      !! Return matrix
      integer, intent(out) :: info
      !! Information flag

      ! internals
      real(dp) :: Xtmp(size(X, 1), size(X, 1))
      real(dp) :: lambda(size(X,1))
      real(dp) :: V(size(X,1), size(X,1))
      integer :: i
      character(len=256) :: msg

      info = 0

      ! Check if the matrix is symmetric
      if (.not. is_symmetric(X)) then
        write(msg,'(A)') "Input matrix is not symmetric."
        call stop_error(msg, module=this_module, procedure='sqrtm_rdp')
      end if

      ! Perform eigenvalue decomposition
      Xtmp = X ; call eigh(Xtmp, lambda, vectors=V, overwrite_a=.true.)

      ! Check if the matrix is positive definite (up to tol)
      do i = 1, size(lambda)
         if (abs(lambda(i)) .gt. 10*atol_dp ) then
            if (lambda(i) .gt. zero_rdp) then
               lambda(i) = sqrt(lambda(i))
            else
               lambda(i) = zero_rdp
               info = -1
            end if
         else
            lambda(i) = zero_rdp
            info = 1
         end if
      end do

      ! Reconstruct the square root matrix
      sqrtmX = matmul(V, matmul(diag(lambda), transpose(V)))

      return
    end subroutine
    subroutine eig_csp(A, vecs, vals)
        !! Eigenvalue decomposition of a dense matrix using LAPACK.
        complex(sp), intent(in) :: A(:, :)
        !! Matrix to be factorized.
        complex(sp), intent(out) :: vecs(:, :)
        !! Eigenvectors.
        complex(sp), intent(out) :: vals(:)

        ! Internal variables
        character :: jobvl = "n", jobvr = "v"
        integer :: n, lwork, info, lda, ldvl, ldvr
        complex(sp) :: A_tilde(size(A, 1), size(A, 2)), vl(1, size(A, 1))
        complex(sp) :: work(2*size(A, 1))
        real(sp) :: rwork(2*size(A, 1))

        ! Setup variables.
        n = size(A, 1) ; lda = n ; ldvl = 1 ; ldvr = n ; a_tilde = a
        lwork = 2*n

        ! Eigendecomposition.
        call geev(jobvl, jobvr, n, a_tilde, lda, vals, vl, ldvl, vecs, ldvr, work, lwork, rwork, info)
        call check_info(info, 'GEEV', module=this_module, procedure='eig_csp')


        return
    end subroutine eig_csp

    subroutine schur_csp(A, Z, eigvals)
        !! Compute the Schur form (in-place) and Schur vectors of the matrix `A`.
        complex(sp), intent(inout) :: A(:, :)
        !! Matrix to be factorized.
        complex(sp), intent(out) :: Z(:, :)
        !! Schur basis.
        complex(sp), intent(out) :: eigvals(:)
        !! Eigenvalues.

        ! Internal variables.
        character :: jobvs = "v", sort = "n"
        integer :: n, lda, sdim, ldvs, lwork, info
        logical, allocatable :: bwork(:)
        complex(sp), allocatable :: work(:)
        real(sp), allocatable :: rwork(:)

        ! Allocate variables.
        n = size(A, 1) ; lda = n ; ldvs = n ; lwork =  2*n 
        allocate(bwork(n)) ; allocate(work(lwork)) ;  allocate(rwork(n)) 

        call gees(jobvs, sort, dummy_select, n, A, lda, sdim, eigvals, Z, ldvs, work, lwork, rwork, bwork, info)
        call check_info(info, 'GEES', module=this_module, procedure='schur_csp')


        return
    contains
        pure function dummy_select(w) result(out)
            complex(sp), intent(in) :: w
            logical :: out
            out = .false.
            return
        end function
    end subroutine schur_csp

    subroutine ordschur_csp(T, Q, selected)
        !! Re-order the Schur factorization from `schur` such that the selected eigenvalues
        !! are in the upper-left block.
        complex(sp), intent(inout) :: T(:, :)
        !! Schur matrix to be re-ordered.
        complex(sp), intent(inout) :: Q(:, :)
        !! Schur vectors to be re-ordered.
        logical, intent(in) :: selected(:)
        !! Boolean array defining the selected eigenvalues.

        ! Internal variables
        character :: job="n", compq="v"
        integer info, ldq, ldt, lwork, m, n
        real(sp) :: s, sep
        complex(sp) :: w(size(T, 1)), work(size(T, 1))

        ! Setup variables.
        n = size(T, 2) ; ldt = n ; ldq = n ; lwork = max(1, n)

        call trsen(job, compq, selected, n, T, ldt, Q, ldq, w, m, s, sep, work, lwork, info)
        call check_info(info, 'TRSEN', module=this_module, procedure='ordschur_csp')

        return
    end subroutine ordschur_csp

    subroutine sqrtm_csp(X, sqrtmX, info)
      !! Matrix-valued sqrt function for dense symmetric/hermitian positive (semi-)definite matrices
      complex(sp), intent(inout) :: X(:,:)
      !! Matrix of which to compute the sqrt
      complex(sp), intent(out)   :: sqrtmX(size(X,1),size(X,1))
      !! Return matrix
      integer, intent(out) :: info
      !! Information flag

      ! internal
      real(sp) :: S(size(X,1))
      complex(sp) :: U(size(X,1), size(X,1)), VT(size(X,1), size(X,1))
      integer :: i
      real(sp) :: symmetry_error
      character(len=256) :: msg

      info = 0

      ! Check if the matrix is hermitian
      symmetry_error = 0.5_sp*maxval(abs(X - conjg(transpose(X))))
      if (symmetry_error > rtol_sp) then
        write(msg,'(2(A,E9.2))') "Input matrix is not hermitian. 0.5*max(abs(X-X.H)) = ", &
            & symmetry_error, ", tol = ", rtol_sp
        call stop_error(msg, module=this_module, procedure='sqrtm_csp')
      else if (symmetry_error > 10*atol_sp) then
        write(msg,'(A,E9.2)') "Input matrix is not exactly hermitian. 0.5*max(X-X.T) = ", symmetry_error
        call logger%log_warning(msg, module=this_module, procedure='sqrtm_csp')
      end if

      ! Perform svd
      call svd(X, S, U, VT)

      ! Check if the matrix is positive definite (up to tol)
      do i = 1, size(S)
         if (S(i) .gt. 10*atol_sp ) then
            S(i) = sqrt(S(i))
         else
            S(i) = zero_rsp
            info = 1
         end if
      end do

      ! Reconstruct the square root matrix
      sqrtmX = matmul(U, matmul(diag(S), VT))

      return
    end subroutine

    subroutine sqrtm_eig_csp(X, sqrtmX, info)
      !! Matrix-valued sqrt function for dense symmetric/hermitian positive (semi-)definite matrices
      complex(sp), intent(in)  :: X(:,:)
      !! Matrix of which to compute the sqrt
      complex(sp), intent(out) :: sqrtmX(size(X,1),size(X,1))
      !! Return matrix
      integer, intent(out) :: info
      !! Information flag

      ! internals
      complex(sp) :: Xtmp(size(X, 1), size(X, 1))
      real(sp) :: lambda(size(X,1))
      complex(sp) :: V(size(X,1), size(X,1))
      integer :: i
      character(len=256) :: msg

      info = 0

      ! Check if the matrix is hermitian
      if (.not. is_hermitian(X)) then
        write(msg,'(A)') "Input matrix is not hermitian"
        call stop_error(msg, module=this_module, procedure='sqrtm_csp')
      end if

      ! Perform eigenvalue decomposition
      Xtmp = X ; call eigh(Xtmp, lambda, vectors=V, overwrite_a=.true.)

      ! Check if the matrix is positive definite (up to tol)
      do i = 1, size(lambda)
         if (abs(lambda(i)) .gt. 10*atol_sp ) then
            if (lambda(i) .gt. zero_rsp) then
               lambda(i) = sqrt(lambda(i))
            else
               lambda(i) = zero_rsp
               info = -1
            end if
         else
            lambda(i) = zero_rsp
            info = 1
         end if
      end do

      ! Reconstruct the square root matrix
      sqrtmX = matmul(V, matmul(diag(lambda), conjg(transpose(V))))

      return
    end subroutine
    subroutine eig_cdp(A, vecs, vals)
        !! Eigenvalue decomposition of a dense matrix using LAPACK.
        complex(dp), intent(in) :: A(:, :)
        !! Matrix to be factorized.
        complex(dp), intent(out) :: vecs(:, :)
        !! Eigenvectors.
        complex(dp), intent(out) :: vals(:)

        ! Internal variables
        character :: jobvl = "n", jobvr = "v"
        integer :: n, lwork, info, lda, ldvl, ldvr
        complex(dp) :: A_tilde(size(A, 1), size(A, 2)), vl(1, size(A, 1))
        complex(dp) :: work(2*size(A, 1))
        real(dp) :: rwork(2*size(A, 1))

        ! Setup variables.
        n = size(A, 1) ; lda = n ; ldvl = 1 ; ldvr = n ; a_tilde = a
        lwork = 2*n

        ! Eigendecomposition.
        call geev(jobvl, jobvr, n, a_tilde, lda, vals, vl, ldvl, vecs, ldvr, work, lwork, rwork, info)
        call check_info(info, 'GEEV', module=this_module, procedure='eig_cdp')


        return
    end subroutine eig_cdp

    subroutine schur_cdp(A, Z, eigvals)
        !! Compute the Schur form (in-place) and Schur vectors of the matrix `A`.
        complex(dp), intent(inout) :: A(:, :)
        !! Matrix to be factorized.
        complex(dp), intent(out) :: Z(:, :)
        !! Schur basis.
        complex(dp), intent(out) :: eigvals(:)
        !! Eigenvalues.

        ! Internal variables.
        character :: jobvs = "v", sort = "n"
        integer :: n, lda, sdim, ldvs, lwork, info
        logical, allocatable :: bwork(:)
        complex(dp), allocatable :: work(:)
        real(dp), allocatable :: rwork(:)

        ! Allocate variables.
        n = size(A, 1) ; lda = n ; ldvs = n ; lwork =  2*n 
        allocate(bwork(n)) ; allocate(work(lwork)) ;  allocate(rwork(n)) 

        call gees(jobvs, sort, dummy_select, n, A, lda, sdim, eigvals, Z, ldvs, work, lwork, rwork, bwork, info)
        call check_info(info, 'GEES', module=this_module, procedure='schur_cdp')


        return
    contains
        pure function dummy_select(w) result(out)
            complex(dp), intent(in) :: w
            logical :: out
            out = .false.
            return
        end function
    end subroutine schur_cdp

    subroutine ordschur_cdp(T, Q, selected)
        !! Re-order the Schur factorization from `schur` such that the selected eigenvalues
        !! are in the upper-left block.
        complex(dp), intent(inout) :: T(:, :)
        !! Schur matrix to be re-ordered.
        complex(dp), intent(inout) :: Q(:, :)
        !! Schur vectors to be re-ordered.
        logical, intent(in) :: selected(:)
        !! Boolean array defining the selected eigenvalues.

        ! Internal variables
        character :: job="n", compq="v"
        integer info, ldq, ldt, lwork, m, n
        real(dp) :: s, sep
        complex(dp) :: w(size(T, 1)), work(size(T, 1))

        ! Setup variables.
        n = size(T, 2) ; ldt = n ; ldq = n ; lwork = max(1, n)

        call trsen(job, compq, selected, n, T, ldt, Q, ldq, w, m, s, sep, work, lwork, info)
        call check_info(info, 'TRSEN', module=this_module, procedure='ordschur_cdp')

        return
    end subroutine ordschur_cdp

    subroutine sqrtm_cdp(X, sqrtmX, info)
      !! Matrix-valued sqrt function for dense symmetric/hermitian positive (semi-)definite matrices
      complex(dp), intent(inout) :: X(:,:)
      !! Matrix of which to compute the sqrt
      complex(dp), intent(out)   :: sqrtmX(size(X,1),size(X,1))
      !! Return matrix
      integer, intent(out) :: info
      !! Information flag

      ! internal
      real(dp) :: S(size(X,1))
      complex(dp) :: U(size(X,1), size(X,1)), VT(size(X,1), size(X,1))
      integer :: i
      real(dp) :: symmetry_error
      character(len=256) :: msg

      info = 0

      ! Check if the matrix is hermitian
      symmetry_error = 0.5_dp*maxval(abs(X - conjg(transpose(X))))
      if (symmetry_error > rtol_dp) then
        write(msg,'(2(A,E9.2))') "Input matrix is not hermitian. 0.5*max(abs(X-X.H)) = ", &
            & symmetry_error, ", tol = ", rtol_dp
        call stop_error(msg, module=this_module, procedure='sqrtm_cdp')
      else if (symmetry_error > 10*atol_dp) then
        write(msg,'(A,E9.2)') "Input matrix is not exactly hermitian. 0.5*max(X-X.T) = ", symmetry_error
        call logger%log_warning(msg, module=this_module, procedure='sqrtm_cdp')
      end if

      ! Perform svd
      call svd(X, S, U, VT)

      ! Check if the matrix is positive definite (up to tol)
      do i = 1, size(S)
         if (S(i) .gt. 10*atol_dp ) then
            S(i) = sqrt(S(i))
         else
            S(i) = zero_rdp
            info = 1
         end if
      end do

      ! Reconstruct the square root matrix
      sqrtmX = matmul(U, matmul(diag(S), VT))

      return
    end subroutine

    subroutine sqrtm_eig_cdp(X, sqrtmX, info)
      !! Matrix-valued sqrt function for dense symmetric/hermitian positive (semi-)definite matrices
      complex(dp), intent(in)  :: X(:,:)
      !! Matrix of which to compute the sqrt
      complex(dp), intent(out) :: sqrtmX(size(X,1),size(X,1))
      !! Return matrix
      integer, intent(out) :: info
      !! Information flag

      ! internals
      complex(dp) :: Xtmp(size(X, 1), size(X, 1))
      real(dp) :: lambda(size(X,1))
      complex(dp) :: V(size(X,1), size(X,1))
      integer :: i
      character(len=256) :: msg

      info = 0

      ! Check if the matrix is hermitian
      if (.not. is_hermitian(X)) then
        write(msg,'(A)') "Input matrix is not hermitian"
        call stop_error(msg, module=this_module, procedure='sqrtm_cdp')
      end if

      ! Perform eigenvalue decomposition
      Xtmp = X ; call eigh(Xtmp, lambda, vectors=V, overwrite_a=.true.)

      ! Check if the matrix is positive definite (up to tol)
      do i = 1, size(lambda)
         if (abs(lambda(i)) .gt. 10*atol_dp ) then
            if (lambda(i) .gt. zero_rdp) then
               lambda(i) = sqrt(lambda(i))
            else
               lambda(i) = zero_rdp
               info = -1
            end if
         else
            lambda(i) = zero_rdp
            info = 1
         end if
      end do

      ! Reconstruct the square root matrix
      sqrtmX = matmul(V, matmul(diag(lambda), conjg(transpose(V))))

      return
    end subroutine

    !---------------------------------
    !-----     MISCELLANEOUS     -----
    !---------------------------------

    pure real(sp) function log2_rsp(x) result(y)
        real(sp), intent(in) :: x
        y = log(x) / log(2.0_sp)
    end function

    pure real(sp) function norml_rsp(A) result(norm)
        real(sp), intent(in) :: A(:, :)
        integer :: i, n
        real(sp) :: row_sum

        norm = zero_rsp
        n = size(A, 1)
        do i = 1, n
            row_sum = sum(abs(A(i, :)))
            norm = max(norm, row_sum)
        enddo
    end function

    pure real(dp) function log2_rdp(x) result(y)
        real(dp), intent(in) :: x
        y = log(x) / log(2.0_dp)
    end function

    pure real(dp) function norml_rdp(A) result(norm)
        real(dp), intent(in) :: A(:, :)
        integer :: i, n
        real(dp) :: row_sum

        norm = zero_rdp
        n = size(A, 1)
        do i = 1, n
            row_sum = sum(abs(A(i, :)))
            norm = max(norm, row_sum)
        enddo
    end function


    pure real(sp) function norml_csp(A) result(norm)
        complex(sp), intent(in) :: A(:, :)
        integer :: i, n
        real(sp) :: row_sum

        norm = zero_rsp
        n = size(A, 1)
        do i = 1, n
            row_sum = sum(abs(A(i, :)))
            norm = max(norm, row_sum)
        enddo
    end function


    pure real(dp) function norml_cdp(A) result(norm)
        complex(dp), intent(in) :: A(:, :)
        integer :: i, n
        real(dp) :: row_sum

        norm = zero_rdp
        n = size(A, 1)
        do i = 1, n
            row_sum = sum(abs(A(i, :)))
            norm = max(norm, row_sum)
        enddo
    end function


end module LightKrylov_Utils
