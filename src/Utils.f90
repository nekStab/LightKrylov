module lightkrylov_utils
    !--------------------------------------------
    !-----     Standard Fortran Library     -----
    !--------------------------------------------
    use iso_fortran_env, only: output_unit
    use stdlib_linalg, only: is_hermitian, is_symmetric, diag, svd
    ! Matrix inversion.
    use stdlib_linalg_lapack, only: getrf, getri
    ! Eigenvalue problem (general + symmetric).
    use stdlib_linalg_lapack, only: geev, syev, heev
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
    ! Compute B = inv(A) in-place for dense matrices.
    public :: inv
    ! Compute AX = XD for general dense matrices.
    public :: eig
    ! Compute AX = XD for symmetric/hermitian matrices.
    public :: eigh
    ! Compute matrix sqrt of input symmetric/hermitian positive definite matrix A
    public :: sqrtm
    public :: sqrtm_eig
    ! Compute AX = XS where S is in Schur form.
    public :: schur
    ! Re-orders the Schur factorization of A.
    public :: ordschur

    interface assert_shape
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
        module procedure norml_rsp
        module procedure norml_rdp
        module procedure norml_csp
        module procedure norml_cdp
    end interface

    interface log2
        module procedure log2_rsp
        module procedure log2_rdp
    end interface

    interface inv
        module procedure inv_rsp
        module procedure inv_rdp
        module procedure inv_csp
        module procedure inv_cdp
    end interface

    interface eig
        module procedure eig_rsp
        module procedure eig_rdp
        module procedure eig_csp
        module procedure eig_cdp
    end interface

   interface eigh
        module procedure eigh_rsp
        module procedure eigh_rdp
        module procedure eigh_csp
        module procedure eigh_cdp
    end interface

    interface schur
        module procedure schur_rsp
        module procedure schur_rdp
        module procedure schur_csp
        module procedure schur_cdp
    end interface

    interface ordschur
        module procedure ordschur_rsp
        module procedure ordschur_rdp
        module procedure ordschur_csp
        module procedure ordschur_cdp
    end interface

    interface sqrtm
        module procedure sqrtm_rsp
        module procedure sqrtm_rdp
        module procedure sqrtm_csp
        module procedure sqrtm_cdp
    end interface

    interface sqrtm_eig
        module procedure sqrtm_eig_rsp
        module procedure sqrtm_eig_rdp
        module procedure sqrtm_eig_csp
        module procedure sqrtm_eig_cdp
    end interface

    !------------------------------------------------
    !-----     OPTS TYPE FOR LINEAR SOLVERS     -----
    !------------------------------------------------

    type, abstract, public :: abstract_opts
        !! Abstract type container for options from which all other are being extended.
    end type

    type, extends(abstract_opts), public :: gmres_sp_opts
        !! GMRES options.
        integer :: kdim = 30
        !! Dimension of the Krylov subspace (default: 30).
        integer :: maxiter = 10
        !! Maximum number of `gmres` restarts (default: 10).
        real(sp) :: atol = atol_sp
        !! Absolute tolerance.
        real(sp) :: rtol = rtol_sp
        !! Relative tolerance.
        logical :: verbose = .false.
        !! Verbosity control (default: `.false.`)
    end type

    type, extends(abstract_opts), public :: cg_sp_opts
        !! Conjugate gradient options.
        integer :: maxiter = 100
        !! Maximum number of `cg` iterations (default: 100).
        real(sp) :: atol = atol_sp
        !! Absolute tolerance.
        real(sp) :: rtol = rtol_sp
        !! Relative tolerance.
        logical :: verbose = .false.
        !! Verbosity control (default: `.false.`)
    end type

    type, extends(abstract_opts), public :: gmres_dp_opts
        !! GMRES options.
        integer :: kdim = 30
        !! Dimension of the Krylov subspace (default: 30).
        integer :: maxiter = 10
        !! Maximum number of `gmres` restarts (default: 10).
        real(dp) :: atol = atol_dp
        !! Absolute tolerance.
        real(dp) :: rtol = rtol_dp
        !! Relative tolerance.
        logical :: verbose = .false.
        !! Verbosity control (default: `.false.`)
    end type

    type, extends(abstract_opts), public :: cg_dp_opts
        !! Conjugate gradient options.
        integer :: maxiter = 100
        !! Maximum number of `cg` iterations (default: 100).
        real(dp) :: atol = atol_dp
        !! Absolute tolerance.
        real(dp) :: rtol = rtol_dp
        !! Relative tolerance.
        logical :: verbose = .false.
        !! Verbosity control (default: `.false.`)
    end type


contains

    !-------------------------------------
    !-----     VARIOUS UTILITIES     -----
    !-------------------------------------

    subroutine assert_shape_vector_rsp(v, size, routine, matname)
        !! Utility function to assert the shape of a vector.
        real(sp), intent(in) :: v(:)
        !! Vector whose dimension need to be asserted.
        integer, intent(in) :: size(:)
        !! Expected dimensions of v.
        character(len=*), intent(in) :: routine
        !! Name of the routine where assertion is done.
        character(len=*), intent(in) :: matname
        !! Name of the asserted vector.
        
        ! internals
        character(len=256) :: msg

        if(any(shape(v) /= size)) then
            write(msg, *) "In routine "//routine//" vector "//matname//" has illegal length ", shape(v), &
                           & ". Expected length is ", size, ". Aborting due to illegal vector length."
            call stop_error(msg, module=this_module, procedure='assert_shape_vector_rsp')
        endif
        return
    end subroutine assert_shape_vector_rsp

    subroutine assert_shape_matrix_rsp(A, size, routine, matname)
        !! Utility function to assert the shape of a matrix.
        real(sp), intent(in) :: A(:, :)
        !! Matrix whose dimension need to be asserted.
        integer, intent(in) :: size(:)
        !! Expected dimensions of A.
        character(len=*), intent(in) :: routine
        !! Name of the routine where assertion is done.
        character(len=*), intent(in) :: matname
        !! Name of the asserted matrix.

        ! internals
        character(len=256) :: msg

        if(any(shape(A) /= size)) then
            write(msg, *) "In routine "//routine//" matrix "//matname//" has illegal shape ", shape(A), &
                        & ". Expected shape is ", size, ". Aborting due to illegal vector length."
            call stop_error(msg, module=this_module, procedure='assert_shape_vector_rsp')
        endif
        return
    end subroutine assert_shape_matrix_rsp
    subroutine assert_shape_vector_rdp(v, size, routine, matname)
        !! Utility function to assert the shape of a vector.
        real(dp), intent(in) :: v(:)
        !! Vector whose dimension need to be asserted.
        integer, intent(in) :: size(:)
        !! Expected dimensions of v.
        character(len=*), intent(in) :: routine
        !! Name of the routine where assertion is done.
        character(len=*), intent(in) :: matname
        !! Name of the asserted vector.
        
        ! internals
        character(len=256) :: msg

        if(any(shape(v) /= size)) then
            write(msg, *) "In routine "//routine//" vector "//matname//" has illegal length ", shape(v), &
                           & ". Expected length is ", size, ". Aborting due to illegal vector length."
            call stop_error(msg, module=this_module, procedure='assert_shape_vector_rdp')
        endif
        return
    end subroutine assert_shape_vector_rdp

    subroutine assert_shape_matrix_rdp(A, size, routine, matname)
        !! Utility function to assert the shape of a matrix.
        real(dp), intent(in) :: A(:, :)
        !! Matrix whose dimension need to be asserted.
        integer, intent(in) :: size(:)
        !! Expected dimensions of A.
        character(len=*), intent(in) :: routine
        !! Name of the routine where assertion is done.
        character(len=*), intent(in) :: matname
        !! Name of the asserted matrix.

        ! internals
        character(len=256) :: msg

        if(any(shape(A) /= size)) then
            write(msg, *) "In routine "//routine//" matrix "//matname//" has illegal shape ", shape(A), &
                        & ". Expected shape is ", size, ". Aborting due to illegal vector length."
            call stop_error(msg, module=this_module, procedure='assert_shape_vector_rdp')
        endif
        return
    end subroutine assert_shape_matrix_rdp
    subroutine assert_shape_vector_csp(v, size, routine, matname)
        !! Utility function to assert the shape of a vector.
        complex(sp), intent(in) :: v(:)
        !! Vector whose dimension need to be asserted.
        integer, intent(in) :: size(:)
        !! Expected dimensions of v.
        character(len=*), intent(in) :: routine
        !! Name of the routine where assertion is done.
        character(len=*), intent(in) :: matname
        !! Name of the asserted vector.
        
        ! internals
        character(len=256) :: msg

        if(any(shape(v) /= size)) then
            write(msg, *) "In routine "//routine//" vector "//matname//" has illegal length ", shape(v), &
                           & ". Expected length is ", size, ". Aborting due to illegal vector length."
            call stop_error(msg, module=this_module, procedure='assert_shape_vector_csp')
        endif
        return
    end subroutine assert_shape_vector_csp

    subroutine assert_shape_matrix_csp(A, size, routine, matname)
        !! Utility function to assert the shape of a matrix.
        complex(sp), intent(in) :: A(:, :)
        !! Matrix whose dimension need to be asserted.
        integer, intent(in) :: size(:)
        !! Expected dimensions of A.
        character(len=*), intent(in) :: routine
        !! Name of the routine where assertion is done.
        character(len=*), intent(in) :: matname
        !! Name of the asserted matrix.

        ! internals
        character(len=256) :: msg

        if(any(shape(A) /= size)) then
            write(msg, *) "In routine "//routine//" matrix "//matname//" has illegal shape ", shape(A), &
                        & ". Expected shape is ", size, ". Aborting due to illegal vector length."
            call stop_error(msg, module=this_module, procedure='assert_shape_vector_csp')
        endif
        return
    end subroutine assert_shape_matrix_csp
    subroutine assert_shape_vector_cdp(v, size, routine, matname)
        !! Utility function to assert the shape of a vector.
        complex(dp), intent(in) :: v(:)
        !! Vector whose dimension need to be asserted.
        integer, intent(in) :: size(:)
        !! Expected dimensions of v.
        character(len=*), intent(in) :: routine
        !! Name of the routine where assertion is done.
        character(len=*), intent(in) :: matname
        !! Name of the asserted vector.
        
        ! internals
        character(len=256) :: msg

        if(any(shape(v) /= size)) then
            write(msg, *) "In routine "//routine//" vector "//matname//" has illegal length ", shape(v), &
                           & ". Expected length is ", size, ". Aborting due to illegal vector length."
            call stop_error(msg, module=this_module, procedure='assert_shape_vector_cdp')
        endif
        return
    end subroutine assert_shape_vector_cdp

    subroutine assert_shape_matrix_cdp(A, size, routine, matname)
        !! Utility function to assert the shape of a matrix.
        complex(dp), intent(in) :: A(:, :)
        !! Matrix whose dimension need to be asserted.
        integer, intent(in) :: size(:)
        !! Expected dimensions of A.
        character(len=*), intent(in) :: routine
        !! Name of the routine where assertion is done.
        character(len=*), intent(in) :: matname
        !! Name of the asserted matrix.

        ! internals
        character(len=256) :: msg

        if(any(shape(A) /= size)) then
            write(msg, *) "In routine "//routine//" matrix "//matname//" has illegal shape ", shape(A), &
                        & ". Expected shape is ", size, ". Aborting due to illegal vector length."
            call stop_error(msg, module=this_module, procedure='assert_shape_vector_cdp')
        endif
        return
    end subroutine assert_shape_matrix_cdp

    !-------------------------------------------
    !-----     LAPACK MATRIX INVERSION     -----
    !-------------------------------------------

    subroutine inv_rsp(A)
        !! In-place inversion of A using LAPACK.
        real(sp), intent(inout) :: A(:, :)
        !! Matrix to be inverted (in-place).

        ! Internal variables.
        integer :: n, info
        real(sp) :: work(size(A, 1))
        integer  :: ipiv(size(A, 1))

        ! Compute A = LU (in-place).
        n = size(A, 1) ; call assert_shape(A, [n, n], "inv", "A")
        call getrf(n, n, A, n, ipiv, info)
        call check_info(info, 'GETREF', module=this_module, procedure='inv_rsp')

        ! Compute inv(A) (in-place).
        call getri(n, A, n, ipiv, work, n, info)
        call check_info(info, 'GETRI', module=this_module, procedure='inv_rsp')

        return
    end subroutine inv_rsp

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

    subroutine eigh_rsp(A, vecs, vals)
        !! Eigenvalue decomposition of a dense symmetric/hermitian matrix using LAPACK.
        real(sp), intent(in) :: A(:, :)
        !! Matrix to be factorized.
        real(sp), intent(out) :: vecs(:, :)
        !! Eigenvectors.
        real(sp), intent(out) :: vals(:)
        !! Eigenvalues.

        ! Internal variables.
        character :: jobz = "v", uplo = "u"
        integer :: n, lwork, info, lda
        real(sp) :: A_tilde(size(A, 1), size(A, 2))
        real(sp), allocatable :: work(:)

        ! Setup variables.
        n = size(A, 1) ; lda = n ; a_tilde = a
        lwork = max(1, 3*n-1)
        allocate(work(lwork))

        ! Eigendecomposition.
        call syev(jobz, uplo, n, a_tilde, lda, vals, work, lwork, info)
        call check_info(info, 'SYEV', module=this_module, procedure='eigh_rsp')

        ! Extract eigenvectors
        vecs = a_tilde

        return
    end subroutine eigh_rsp

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
        pure function dummy_select(wr, wi) result(out)
            real(sp), intent(in) :: wr
            real(sp), intent(in) :: wi
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
      character(len=128) :: msg

      info = 0

      ! Check if the matrix is symmetric
      symmetry_error = 0.5*maxval(X - transpose(X))
      print *, "Input matrix is not symmetric. 0.5*max(X-X.T) = ", symmetry_error
      if (symmetry_error > 10*atol_sp) then
        write(msg,*) "Input matrix is not symmetric. 0.5*max(X-X.T) = ", symmetry_error
        call stop_error(msg, module=this_module, procedure='sqrtm_rsp')
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
      real(sp) :: lambda(size(X,1))
      real(sp) :: V(size(X,1), size(X,1))
      integer :: i
      character(len=128) :: msg

      info = 0

      ! Check if the matrix is symmetric
      if (.not. is_symmetric(X)) then
        write(msg,*) "Input matrix is not symmetric."
        call stop_error(msg, module=this_module, procedure='sqrtm_rsp')
      end if

      ! Perform eigenvalue decomposition
      call eigh(X, V, lambda)

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
    subroutine inv_rdp(A)
        !! In-place inversion of A using LAPACK.
        real(dp), intent(inout) :: A(:, :)
        !! Matrix to be inverted (in-place).

        ! Internal variables.
        integer :: n, info
        real(dp) :: work(size(A, 1))
        integer  :: ipiv(size(A, 1))

        ! Compute A = LU (in-place).
        n = size(A, 1) ; call assert_shape(A, [n, n], "inv", "A")
        call getrf(n, n, A, n, ipiv, info)
        call check_info(info, 'GETREF', module=this_module, procedure='inv_rdp')

        ! Compute inv(A) (in-place).
        call getri(n, A, n, ipiv, work, n, info)
        call check_info(info, 'GETRI', module=this_module, procedure='inv_rdp')

        return
    end subroutine inv_rdp

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

    subroutine eigh_rdp(A, vecs, vals)
        !! Eigenvalue decomposition of a dense symmetric/hermitian matrix using LAPACK.
        real(dp), intent(in) :: A(:, :)
        !! Matrix to be factorized.
        real(dp), intent(out) :: vecs(:, :)
        !! Eigenvectors.
        real(dp), intent(out) :: vals(:)
        !! Eigenvalues.

        ! Internal variables.
        character :: jobz = "v", uplo = "u"
        integer :: n, lwork, info, lda
        real(dp) :: A_tilde(size(A, 1), size(A, 2))
        real(dp), allocatable :: work(:)

        ! Setup variables.
        n = size(A, 1) ; lda = n ; a_tilde = a
        lwork = max(1, 3*n-1)
        allocate(work(lwork))

        ! Eigendecomposition.
        call syev(jobz, uplo, n, a_tilde, lda, vals, work, lwork, info)
        call check_info(info, 'SYEV', module=this_module, procedure='eigh_rdp')

        ! Extract eigenvectors
        vecs = a_tilde

        return
    end subroutine eigh_rdp

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
        pure function dummy_select(wr, wi) result(out)
            real(dp), intent(in) :: wr
            real(dp), intent(in) :: wi
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
      character(len=128) :: msg

      info = 0

      ! Check if the matrix is symmetric
      symmetry_error = 0.5*maxval(X - transpose(X))
      print *, "Input matrix is not symmetric. 0.5*max(X-X.T) = ", symmetry_error
      if (symmetry_error > 10*atol_dp) then
        write(msg,*) "Input matrix is not symmetric. 0.5*max(X-X.T) = ", symmetry_error
        call stop_error(msg, module=this_module, procedure='sqrtm_rdp')
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
      real(dp) :: lambda(size(X,1))
      real(dp) :: V(size(X,1), size(X,1))
      integer :: i
      character(len=128) :: msg

      info = 0

      ! Check if the matrix is symmetric
      if (.not. is_symmetric(X)) then
        write(msg,*) "Input matrix is not symmetric."
        call stop_error(msg, module=this_module, procedure='sqrtm_rdp')
      end if

      ! Perform eigenvalue decomposition
      call eigh(X, V, lambda)

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
    subroutine inv_csp(A)
        !! In-place inversion of A using LAPACK.
        complex(sp), intent(inout) :: A(:, :)
        !! Matrix to be inverted (in-place).

        ! Internal variables.
        integer :: n, info
        complex(sp) :: work(size(A, 1))
        integer  :: ipiv(size(A, 1))

        ! Compute A = LU (in-place).
        n = size(A, 1) ; call assert_shape(A, [n, n], "inv", "A")
        call getrf(n, n, A, n, ipiv, info)
        call check_info(info, 'GETREF', module=this_module, procedure='inv_csp')

        ! Compute inv(A) (in-place).
        call getri(n, A, n, ipiv, work, n, info)
        call check_info(info, 'GETRI', module=this_module, procedure='inv_csp')

        return
    end subroutine inv_csp

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

    subroutine eigh_csp(A, vecs, vals)
        !! Eigenvalue decomposition of a dense symmetric/hermitian matrix using LAPACK.
        complex(sp), intent(in) :: A(:, :)
        !! Matrix to be factorized.
        complex(sp), intent(out) :: vecs(:, :)
        !! Eigenvectors.
        real(sp), intent(out) :: vals(:)
        !! Eigenvalues.

        ! Internal variables.
        character :: jobz = "v", uplo = "u"
        integer :: n, lwork, info, lda
        complex(sp) :: A_tilde(size(A, 1), size(A, 2))
        complex(sp), allocatable :: work(:)
        real(sp), allocatable :: rwork(:)

        ! Setup variables.
        n = size(A, 1) ; lda = n ; a_tilde = a
        lwork = max(1, 2*n-1)
        allocate(rwork(max(1, 3*n-2)))
        allocate(work(lwork))

        ! Eigendecomposition.
        call heev(jobz, uplo, n, a_tilde, lda, vals, work, lwork, rwork, info)
        call check_info(info, 'HEEV', module=this_module, procedure='eigh_csp')

        ! Extract eigenvectors
        vecs = a_tilde

        return
    end subroutine eigh_csp

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
      character(len=128) :: msg

      info = 0

      ! Check if the matrix is hermitian
      symmetry_error = 0.5*maxval(abs(X - conjg(transpose(X))))
      if (symmetry_error > 10*atol_sp) then
        write(msg,*) "Input matrix is not hermitian. 0.5*max(abs(X-X.H)) = ", symmetry_error
        call stop_error(msg, module=this_module, procedure='sqrtm_csp')
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
      real(sp) :: lambda(size(X,1))
      complex(sp) :: V(size(X,1), size(X,1))
      integer :: i
      character(len=128) :: msg

      info = 0

      ! Check if the matrix is hermitian
      if (.not. is_hermitian(X)) then
        write(msg,*) "Input matrix is not hermitian"
        call stop_error(msg, module=this_module, procedure='sqrtm_csp')
      end if

      ! Perform eigenvalue decomposition
      call eigh(X, V, lambda)

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
    subroutine inv_cdp(A)
        !! In-place inversion of A using LAPACK.
        complex(dp), intent(inout) :: A(:, :)
        !! Matrix to be inverted (in-place).

        ! Internal variables.
        integer :: n, info
        complex(dp) :: work(size(A, 1))
        integer  :: ipiv(size(A, 1))

        ! Compute A = LU (in-place).
        n = size(A, 1) ; call assert_shape(A, [n, n], "inv", "A")
        call getrf(n, n, A, n, ipiv, info)
        call check_info(info, 'GETREF', module=this_module, procedure='inv_cdp')

        ! Compute inv(A) (in-place).
        call getri(n, A, n, ipiv, work, n, info)
        call check_info(info, 'GETRI', module=this_module, procedure='inv_cdp')

        return
    end subroutine inv_cdp

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

    subroutine eigh_cdp(A, vecs, vals)
        !! Eigenvalue decomposition of a dense symmetric/hermitian matrix using LAPACK.
        complex(dp), intent(in) :: A(:, :)
        !! Matrix to be factorized.
        complex(dp), intent(out) :: vecs(:, :)
        !! Eigenvectors.
        real(dp), intent(out) :: vals(:)
        !! Eigenvalues.

        ! Internal variables.
        character :: jobz = "v", uplo = "u"
        integer :: n, lwork, info, lda
        complex(dp) :: A_tilde(size(A, 1), size(A, 2))
        complex(dp), allocatable :: work(:)
        real(dp), allocatable :: rwork(:)

        ! Setup variables.
        n = size(A, 1) ; lda = n ; a_tilde = a
        lwork = max(1, 2*n-1)
        allocate(rwork(max(1, 3*n-2)))
        allocate(work(lwork))

        ! Eigendecomposition.
        call heev(jobz, uplo, n, a_tilde, lda, vals, work, lwork, rwork, info)
        call check_info(info, 'HEEV', module=this_module, procedure='eigh_cdp')

        ! Extract eigenvectors
        vecs = a_tilde

        return
    end subroutine eigh_cdp

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
      character(len=128) :: msg

      info = 0

      ! Check if the matrix is hermitian
      symmetry_error = 0.5*maxval(abs(X - conjg(transpose(X))))
      if (symmetry_error > 10*atol_dp) then
        write(msg,*) "Input matrix is not hermitian. 0.5*max(abs(X-X.H)) = ", symmetry_error
        call stop_error(msg, module=this_module, procedure='sqrtm_cdp')
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
      real(dp) :: lambda(size(X,1))
      complex(dp) :: V(size(X,1), size(X,1))
      integer :: i
      character(len=128) :: msg

      info = 0

      ! Check if the matrix is hermitian
      if (.not. is_hermitian(X)) then
        write(msg,*) "Input matrix is not hermitian"
        call stop_error(msg, module=this_module, procedure='sqrtm_cdp')
      end if

      ! Perform eigenvalue decomposition
      call eigh(X, V, lambda)

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

    real(sp) function log2_rsp(x) result(y)
        real(sp), intent(in) :: x
        y = log(x) / log(2.0_sp)
    end function

    real(sp) function norml_rsp(A) result(norm)
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

    real(dp) function log2_rdp(x) result(y)
        real(dp), intent(in) :: x
        y = log(x) / log(2.0_dp)
    end function

    real(dp) function norml_rdp(A) result(norm)
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


    real(sp) function norml_csp(A) result(norm)
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


    real(dp) function norml_cdp(A) result(norm)
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


end module lightkrylov_utils
