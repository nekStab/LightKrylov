module lightkrylov_utils
    !--------------------------------------------
    !-----     Standard Fortran Library     -----
    !--------------------------------------------
    use iso_fortran_env, only: output_unit
    ! Check symmetry
    use stdlib_linalg, only: diag, is_symmetric
    ! Matrix inversion.
    use stdlib_linalg_lapack, only: getrf, getri
    ! Singular value decomposition.
    use stdlib_linalg_lapack, only: gesvd
    ! Eigenvalue problem (general + symmetric).
    use stdlib_linalg_lapack, only: geev, syev, heev
    ! Least-squares solver.
    use stdlib_linalg_lapack, only: gels
    ! Schur factorization.
    use stdlib_linalg_lapack, only: gees, trsen

    !-------------------------------
    !-----     LightKrylov     -----
    !-------------------------------
    ! Various constants.
    use lightkrylov_constants

    implicit none

    private

    real(sp), parameter, public :: one_rsp = 1.0_sp
    real(sp), parameter, public :: zero_rsp = 0.0_sp
    real(dp), parameter, public :: one_rdp = 1.0_dp
    real(dp), parameter, public :: zero_rdp = 0.0_dp
    complex(sp), parameter, public :: one_csp = cmplx(1.0_sp, 0.0_sp, kind=sp)
    complex(sp), parameter, public :: zero_csp = cmplx(0.0_sp, 0.0_sp, kind=sp)
    complex(dp), parameter, public :: one_cdp = cmplx(1.0_dp, 0.0_dp, kind=dp)
    complex(dp), parameter, public :: zero_cdp = cmplx(0.0_dp, 0.0_dp, kind=dp)

    public :: stop_error
    public :: assert_shape
    ! Compute B = inv(A) in-place for dense matrices.
    public :: inv
    ! Compute USV^T = svd(A) for dense matrices.
    public :: svd
    ! Compute AX = XD for general dense matrices.
    public :: eig
    ! Compute AX = XD for symmetric/hermitian matrices.
    public :: eigh
    ! Compute matrix sqrt of input SPD matrix A
    public :: sqrtm
    ! Solve min || Ax - b ||_2^2.
    public :: lstsq
    ! Compute AX = XS where S is in Schur form.
    public :: schur
    ! Re-orders the Schur factorization of A.
    public :: ordschur

    public :: abstract_opts
    public :: gmres_sp_opts
    public :: cg_sp_opts
    public :: gmres_dp_opts
    public :: cg_dp_opts

    public :: log2_rsp
    public :: norml_rsp
    public :: log2_rdp
    public :: norml_rdp
    public :: norml_csp
    public :: norml_cdp

    interface assert_shape
        module procedure assert_shape_rsp
        module procedure assert_shape_rdp
        module procedure assert_shape_csp
        module procedure assert_shape_cdp
    end interface

    interface inv
        module procedure inv_rsp
        module procedure inv_rdp
        module procedure inv_csp
        module procedure inv_cdp
    end interface

    interface svd
        module procedure svd_rsp
        module procedure svd_rdp
        module procedure svd_csp
        module procedure svd_cdp
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

    interface lstsq
        module procedure lstsq_rsp
        module procedure lstsq_rdp
        module procedure lstsq_csp
        module procedure lstsq_cdp
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

    subroutine stop_error(msg)
        !! Utility function to print an error message.
        character(len=*), intent(in) :: msg
        !! Error message.
        write(output_unit, *) msg; stop 1
        return
    end subroutine stop_error

    subroutine assert_shape_rsp(A, size, routine, matname)
        !! Utility function to assert the shape of a matrix.
        real(sp), intent(in) :: A(:, :)
        !! Matrix whose dimension need to be asserted.
        integer, intent(in) :: size(:)
        !! Expected dimensions of A.
        character(len=*), intent(in) :: routine
        !! Name of the routine where assertion is done.
        character(len=*), intent(in) :: matname
        !! Name of the asserted matrix.

        if(any(shape(A) /= size)) then
            write(output_unit, *) "In routine "//routine//" matrix "//matname//" has illegal shape ", shape(A)
            write(output_unit, *) "Expected shape is ", size
            call stop_error("Aborting due to illegal matrix size.")
        endif
        return
    end subroutine assert_shape_rsp
    subroutine assert_shape_rdp(A, size, routine, matname)
        !! Utility function to assert the shape of a matrix.
        real(dp), intent(in) :: A(:, :)
        !! Matrix whose dimension need to be asserted.
        integer, intent(in) :: size(:)
        !! Expected dimensions of A.
        character(len=*), intent(in) :: routine
        !! Name of the routine where assertion is done.
        character(len=*), intent(in) :: matname
        !! Name of the asserted matrix.

        if(any(shape(A) /= size)) then
            write(output_unit, *) "In routine "//routine//" matrix "//matname//" has illegal shape ", shape(A)
            write(output_unit, *) "Expected shape is ", size
            call stop_error("Aborting due to illegal matrix size.")
        endif
        return
    end subroutine assert_shape_rdp
    subroutine assert_shape_csp(A, size, routine, matname)
        !! Utility function to assert the shape of a matrix.
        complex(sp), intent(in) :: A(:, :)
        !! Matrix whose dimension need to be asserted.
        integer, intent(in) :: size(:)
        !! Expected dimensions of A.
        character(len=*), intent(in) :: routine
        !! Name of the routine where assertion is done.
        character(len=*), intent(in) :: matname
        !! Name of the asserted matrix.

        if(any(shape(A) /= size)) then
            write(output_unit, *) "In routine "//routine//" matrix "//matname//" has illegal shape ", shape(A)
            write(output_unit, *) "Expected shape is ", size
            call stop_error("Aborting due to illegal matrix size.")
        endif
        return
    end subroutine assert_shape_csp
    subroutine assert_shape_cdp(A, size, routine, matname)
        !! Utility function to assert the shape of a matrix.
        complex(dp), intent(in) :: A(:, :)
        !! Matrix whose dimension need to be asserted.
        integer, intent(in) :: size(:)
        !! Expected dimensions of A.
        character(len=*), intent(in) :: routine
        !! Name of the routine where assertion is done.
        character(len=*), intent(in) :: matname
        !! Name of the asserted matrix.

        if(any(shape(A) /= size)) then
            write(output_unit, *) "In routine "//routine//" matrix "//matname//" has illegal shape ", shape(A)
            write(output_unit, *) "Expected shape is ", size
            call stop_error("Aborting due to illegal matrix size.")
        endif
        return
    end subroutine assert_shape_cdp

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
        if (info /= 0) then
            write(output_unit, *) "GETRF return info =", info
            if (info<0) then
                write(output_unit, *) "The ", -info, "-th argument has illegal value."
            else
                write(output_unit, *) "U(", info, ",", info, ") is exactly zero. The factorization"
                write(output_unit, *) "has been completed but the factor U is exactly singular."
                write(output_unit, *) "Division by zero will occur if used to solve Ax=b."
            endif
            call stop_error("inv: GETREF error")
        endif

        ! Compute inv(A) (in-place).
        call getri(n, A, n, ipiv, work, n, info)
        if (info /= 0) then
            write(output_unit, *) "GETRI return info = ", info
            if (info < 0) then
                write(output_unit, *) "The ", -info, "-th argument has an illegal value."
            else
                write(output_unit, *) "U(", info, ",", info, ") is exactly zero."
                write(output_unit, *) "The matrix is singular and its inverse cannot be computed."
            endif
            call stop_error("in: GETRI error.")
        endif

        return
    end subroutine inv_rsp

    subroutine svd_rsp(A, U, S, V)
        !! Singular value decomposition of a dense matrix.
        real(sp), intent(in) :: A(:, :)
        !! Matrix to be factorized.
        real(sp), intent(out) :: U(:, :)
        !! Left singular vectors.
        real(sp), intent(out) :: S(:)
        !! Singular values.
        real(sp), intent(out) :: V(:, :)
        !! Right singular vectors.

        ! Internal variables.
        character :: jobu = "S", jobvt = "S"
        integer :: m, n, lda, ldu, ldvt, lwork, info
        real(sp), allocatable :: work(:)
        real(sp) :: A_tilde(size(A, 1), size(A, 2)), Vt(min(size(A, 1), size(A, 2)), size(A, 2))

        ! Setup variables.
        m = size(A, 1) ; n = size(A, 2)
        lda = m ; ldu = m ; ldvt = n
        lwork = max(1, 3*min(m, n) + max(m, n), 5*min(m, n)) ; allocate(work(lwork))

        ! Shape assertion.
        call assert_shape(U, [m, m], "svd", "U")
        call assert_shape(V, [n, n], "svd", "V")

        ! SVD computation.
        A_tilde = A
        call gesvd(jobu, jobvt, m, n, A_tilde, lda, S, U, ldu, Vt, ldvt, work, lwork, info)

        return
    end subroutine svd_rsp

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
        real(sp) :: A_tilde(size(A, 1), size(A, 2)), vr(size(A, 1), size(A, 2)), vl(1, size(A, 2))
        real(sp) :: work(4*size(A, 1)), wr(size(A, 1)), wi(size(A, 1))
        integer :: i, idx(size(A, 1))

        ! Setup variables.
        n = size(A, 1) ; lda = n ; ldvl = 1 ; ldvr = n ; a_tilde = a
        lwork = 4*n

        ! Eigendecomposition.
        call geev(jobvl, jobvr, n, a_tilde, lda, wr, wi, vl, ldvl, vecs, ldvr, work, lwork, info)

        if (info /= 0) then
            write(output_unit, *) "GEEV returned info =", info
            if (info < 0) then
                write(output_unit, *) "The ", -info, "-th argument has illegal value."
            else
                write(output_unit, *) "The QR alg. failed to compute all of the eigenvalues."
                write(output_unit, *) "No eigenvector has been computed."
            endif
            call stop_error("eig: geev error.")
        endif

        vals = cmplx(1.0_sp, 0.0_sp, kind=sp)*wr + cmplx(0.0_sp, 1.0_sp, kind=sp)*wi

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

        return
    end subroutine eigh_rsp

    subroutine lstsq_rsp(A, b, x)
        !! Solves a linear least-squares problem \(\min ~ \| \mathbf{Ax} - \mathbf{b} \|_2^2 \) using LAPACK.
        real(sp), intent(in) :: A(:, :)
        !! Matrix to be "pseudo-inversed".
        real(sp), intent(in) :: b(:)
        !! Right-hand side vector.
        real(sp), intent(out) :: x(:)
        !! Solution of the least-squares problem.

        ! Internal variables.
        character :: trans="n"
        integer :: m, n, nrhs, lda, ldb, lwork, info
        real(sp) :: A_tilde(size(A, 1), size(A, 2)), b_tilde(size(A, 1), 1)
        real(sp), allocatable :: work(:)

        ! Setup variables.
        m = size(A, 1) ; n = size(A, 2) ; nrhs = 1
        lda = m ; ldb = m ; lwork = max(1, min(m, n) + max(min(m, n), nrhs))
        a_tilde = a ; b_tilde(:, 1) = b
        allocate(work(lwork)) ; work = 0.0_sp

        ! Solve the least-squares problem.
        call gels(trans, m, n, nrhs, a_tilde, lda, b_tilde, ldb, work, lwork, info)

        if (info /= 0) then
            write(output_unit, *) "The ", -info, "-th argument has illegal value."
            call stop_error("lstsq: dgels error")
        endif

        ! Return solution.
        x = b_tilde(1:n, 1)

        return
    end subroutine lstsq_rsp

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

        return
    end subroutine ordschur_rsp

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
        if (info /= 0) then
            write(output_unit, *) "GETRF return info =", info
            if (info<0) then
                write(output_unit, *) "The ", -info, "-th argument has illegal value."
            else
                write(output_unit, *) "U(", info, ",", info, ") is exactly zero. The factorization"
                write(output_unit, *) "has been completed but the factor U is exactly singular."
                write(output_unit, *) "Division by zero will occur if used to solve Ax=b."
            endif
            call stop_error("inv: GETREF error")
        endif

        ! Compute inv(A) (in-place).
        call getri(n, A, n, ipiv, work, n, info)
        if (info /= 0) then
            write(output_unit, *) "GETRI return info = ", info
            if (info < 0) then
                write(output_unit, *) "The ", -info, "-th argument has an illegal value."
            else
                write(output_unit, *) "U(", info, ",", info, ") is exactly zero."
                write(output_unit, *) "The matrix is singular and its inverse cannot be computed."
            endif
            call stop_error("in: GETRI error.")
        endif

        return
    end subroutine inv_rdp

    subroutine svd_rdp(A, U, S, V)
        !! Singular value decomposition of a dense matrix.
        real(dp), intent(in) :: A(:, :)
        !! Matrix to be factorized.
        real(dp), intent(out) :: U(:, :)
        !! Left singular vectors.
        real(dp), intent(out) :: S(:)
        !! Singular values.
        real(dp), intent(out) :: V(:, :)
        !! Right singular vectors.

        ! Internal variables.
        character :: jobu = "S", jobvt = "S"
        integer :: m, n, lda, ldu, ldvt, lwork, info
        real(dp), allocatable :: work(:)
        real(dp) :: A_tilde(size(A, 1), size(A, 2)), Vt(min(size(A, 1), size(A, 2)), size(A, 2))

        ! Setup variables.
        m = size(A, 1) ; n = size(A, 2)
        lda = m ; ldu = m ; ldvt = n
        lwork = max(1, 3*min(m, n) + max(m, n), 5*min(m, n)) ; allocate(work(lwork))

        ! Shape assertion.
        call assert_shape(U, [m, m], "svd", "U")
        call assert_shape(V, [n, n], "svd", "V")

        ! SVD computation.
        A_tilde = A
        call gesvd(jobu, jobvt, m, n, A_tilde, lda, S, U, ldu, Vt, ldvt, work, lwork, info)

        return
    end subroutine svd_rdp

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
        real(dp) :: A_tilde(size(A, 1), size(A, 2)), vr(size(A, 1), size(A, 2)), vl(1, size(A, 2))
        real(dp) :: work(4*size(A, 1)), wr(size(A, 1)), wi(size(A, 1))
        integer :: i, idx(size(A, 1))

        ! Setup variables.
        n = size(A, 1) ; lda = n ; ldvl = 1 ; ldvr = n ; a_tilde = a
        lwork = 4*n

        ! Eigendecomposition.
        call geev(jobvl, jobvr, n, a_tilde, lda, wr, wi, vl, ldvl, vecs, ldvr, work, lwork, info)

        if (info /= 0) then
            write(output_unit, *) "GEEV returned info =", info
            if (info < 0) then
                write(output_unit, *) "The ", -info, "-th argument has illegal value."
            else
                write(output_unit, *) "The QR alg. failed to compute all of the eigenvalues."
                write(output_unit, *) "No eigenvector has been computed."
            endif
            call stop_error("eig: geev error.")
        endif

        vals = cmplx(1.0_dp, 0.0_dp, kind=dp)*wr + cmplx(0.0_dp, 1.0_dp, kind=dp)*wi

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

        return
    end subroutine eigh_rdp

    subroutine lstsq_rdp(A, b, x)
        !! Solves a linear least-squares problem \(\min ~ \| \mathbf{Ax} - \mathbf{b} \|_2^2 \) using LAPACK.
        real(dp), intent(in) :: A(:, :)
        !! Matrix to be "pseudo-inversed".
        real(dp), intent(in) :: b(:)
        !! Right-hand side vector.
        real(dp), intent(out) :: x(:)
        !! Solution of the least-squares problem.

        ! Internal variables.
        character :: trans="n"
        integer :: m, n, nrhs, lda, ldb, lwork, info
        real(dp) :: A_tilde(size(A, 1), size(A, 2)), b_tilde(size(A, 1), 1)
        real(dp), allocatable :: work(:)

        ! Setup variables.
        m = size(A, 1) ; n = size(A, 2) ; nrhs = 1
        lda = m ; ldb = m ; lwork = max(1, min(m, n) + max(min(m, n), nrhs))
        a_tilde = a ; b_tilde(:, 1) = b
        allocate(work(lwork)) ; work = 0.0_dp

        ! Solve the least-squares problem.
        call gels(trans, m, n, nrhs, a_tilde, lda, b_tilde, ldb, work, lwork, info)

        if (info /= 0) then
            write(output_unit, *) "The ", -info, "-th argument has illegal value."
            call stop_error("lstsq: dgels error")
        endif

        ! Return solution.
        x = b_tilde(1:n, 1)

        return
    end subroutine lstsq_rdp

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

        return
    end subroutine ordschur_rdp

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
        if (info /= 0) then
            write(output_unit, *) "GETRF return info =", info
            if (info<0) then
                write(output_unit, *) "The ", -info, "-th argument has illegal value."
            else
                write(output_unit, *) "U(", info, ",", info, ") is exactly zero. The factorization"
                write(output_unit, *) "has been completed but the factor U is exactly singular."
                write(output_unit, *) "Division by zero will occur if used to solve Ax=b."
            endif
            call stop_error("inv: GETREF error")
        endif

        ! Compute inv(A) (in-place).
        call getri(n, A, n, ipiv, work, n, info)
        if (info /= 0) then
            write(output_unit, *) "GETRI return info = ", info
            if (info < 0) then
                write(output_unit, *) "The ", -info, "-th argument has an illegal value."
            else
                write(output_unit, *) "U(", info, ",", info, ") is exactly zero."
                write(output_unit, *) "The matrix is singular and its inverse cannot be computed."
            endif
            call stop_error("in: GETRI error.")
        endif

        return
    end subroutine inv_csp

    subroutine svd_csp(A, U, S, V)
        !! Singular value decomposition of a dense matrix.
        complex(sp), intent(in) :: A(:, :)
        !! Matrix to be factorized.
        complex(sp), intent(out) :: U(:, :)
        !! Left singular vectors.
        real(sp), intent(out) :: S(:)
        !! Singular values.
        complex(sp), intent(out) :: V(:, :)
        !! Right singular vectors.

        ! Internal variables.
        character :: jobu = "S", jobvt = "S"
        integer :: m, n, lda, ldu, ldvt, lwork, info
        complex(sp), allocatable :: work(:)
        complex(sp) :: A_tilde(size(A, 1), size(A, 2)), Vt(min(size(A, 1), size(A, 2)), size(A, 2))
        real(sp), allocatable :: rwork(:)

        ! Setup variables.
        m = size(A, 1) ; n = size(A, 2)
        lda = m ; ldu = m ; ldvt = n
        lwork = max(1, 3*min(m, n) + max(m, n), 5*min(m, n)) ; allocate(work(lwork))

        ! Shape assertion.
        call assert_shape(U, [m, m], "svd", "U")
        call assert_shape(V, [n, n], "svd", "V")

        ! SVD computation.
        A_tilde = A
        allocate(rwork(5*min(m, n)))
        call gesvd(jobu, jobvt, m, n, A_tilde, lda, S, U, ldu, Vt, ldvt, work, lwork, rwork, info)

        return
    end subroutine svd_csp

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
        complex(sp) :: A_tilde(size(A, 1), size(A, 2)), vr(size(A, 1), size(A, 2)), vl(1, size(A, 1))
        complex(sp) :: work(2*size(A, 1)), w(size(A, 1))
        real(sp) :: rwork(2*size(A, 1))
        integer :: i, idx(size(A, 1))

        ! Setup variables.
        n = size(A, 1) ; lda = n ; ldvl = 1 ; ldvr = n ; a_tilde = a
        lwork = 2*n

        ! Eigendecomposition.
        call geev(jobvl, jobvr, n, a_tilde, lda, vals, vl, ldvl, vecs, ldvr, work, lwork, rwork, info)

        if (info /= 0) then
            write(output_unit, *) "GEEV returned info =", info
            if (info < 0) then
                write(output_unit, *) "The ", -info, "-th argument has illegal value."
            else
                write(output_unit, *) "The QR alg. failed to compute all of the eigenvalues."
                write(output_unit, *) "No eigenvector has been computed."
            endif
            call stop_error("eig: geev error.")
        endif


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

        return
    end subroutine eigh_csp

    subroutine lstsq_csp(A, b, x)
        !! Solves a linear least-squares problem \(\min ~ \| \mathbf{Ax} - \mathbf{b} \|_2^2 \) using LAPACK.
        complex(sp), intent(in) :: A(:, :)
        !! Matrix to be "pseudo-inversed".
        complex(sp), intent(in) :: b(:)
        !! Right-hand side vector.
        complex(sp), intent(out) :: x(:)
        !! Solution of the least-squares problem.

        ! Internal variables.
        character :: trans="n"
        integer :: m, n, nrhs, lda, ldb, lwork, info
        complex(sp) :: A_tilde(size(A, 1), size(A, 2)), b_tilde(size(A, 1), 1)
        complex(sp), allocatable :: work(:)

        ! Setup variables.
        m = size(A, 1) ; n = size(A, 2) ; nrhs = 1
        lda = m ; ldb = m ; lwork = max(1, min(m, n) + max(min(m, n), nrhs))
        a_tilde = a ; b_tilde(:, 1) = b
        allocate(work(lwork)) ; work = 0.0_sp

        ! Solve the least-squares problem.
        call gels(trans, m, n, nrhs, a_tilde, lda, b_tilde, ldb, work, lwork, info)

        if (info /= 0) then
            write(output_unit, *) "The ", -info, "-th argument has illegal value."
            call stop_error("lstsq: dgels error")
        endif

        ! Return solution.
        x = b_tilde(1:n, 1)

        return
    end subroutine lstsq_csp

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

        return
    end subroutine ordschur_csp

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
        if (info /= 0) then
            write(output_unit, *) "GETRF return info =", info
            if (info<0) then
                write(output_unit, *) "The ", -info, "-th argument has illegal value."
            else
                write(output_unit, *) "U(", info, ",", info, ") is exactly zero. The factorization"
                write(output_unit, *) "has been completed but the factor U is exactly singular."
                write(output_unit, *) "Division by zero will occur if used to solve Ax=b."
            endif
            call stop_error("inv: GETREF error")
        endif

        ! Compute inv(A) (in-place).
        call getri(n, A, n, ipiv, work, n, info)
        if (info /= 0) then
            write(output_unit, *) "GETRI return info = ", info
            if (info < 0) then
                write(output_unit, *) "The ", -info, "-th argument has an illegal value."
            else
                write(output_unit, *) "U(", info, ",", info, ") is exactly zero."
                write(output_unit, *) "The matrix is singular and its inverse cannot be computed."
            endif
            call stop_error("in: GETRI error.")
        endif

        return
    end subroutine inv_cdp

    subroutine svd_cdp(A, U, S, V)
        !! Singular value decomposition of a dense matrix.
        complex(dp), intent(in) :: A(:, :)
        !! Matrix to be factorized.
        complex(dp), intent(out) :: U(:, :)
        !! Left singular vectors.
        real(dp), intent(out) :: S(:)
        !! Singular values.
        complex(dp), intent(out) :: V(:, :)
        !! Right singular vectors.

        ! Internal variables.
        character :: jobu = "S", jobvt = "S"
        integer :: m, n, lda, ldu, ldvt, lwork, info
        complex(dp), allocatable :: work(:)
        complex(dp) :: A_tilde(size(A, 1), size(A, 2)), Vt(min(size(A, 1), size(A, 2)), size(A, 2))
        real(dp), allocatable :: rwork(:)

        ! Setup variables.
        m = size(A, 1) ; n = size(A, 2)
        lda = m ; ldu = m ; ldvt = n
        lwork = max(1, 3*min(m, n) + max(m, n), 5*min(m, n)) ; allocate(work(lwork))

        ! Shape assertion.
        call assert_shape(U, [m, m], "svd", "U")
        call assert_shape(V, [n, n], "svd", "V")

        ! SVD computation.
        A_tilde = A
        allocate(rwork(5*min(m, n)))
        call gesvd(jobu, jobvt, m, n, A_tilde, lda, S, U, ldu, Vt, ldvt, work, lwork, rwork, info)

        return
    end subroutine svd_cdp

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
        complex(dp) :: A_tilde(size(A, 1), size(A, 2)), vr(size(A, 1), size(A, 2)), vl(1, size(A, 1))
        complex(dp) :: work(2*size(A, 1)), w(size(A, 1))
        real(dp) :: rwork(2*size(A, 1))
        integer :: i, idx(size(A, 1))

        ! Setup variables.
        n = size(A, 1) ; lda = n ; ldvl = 1 ; ldvr = n ; a_tilde = a
        lwork = 2*n

        ! Eigendecomposition.
        call geev(jobvl, jobvr, n, a_tilde, lda, vals, vl, ldvl, vecs, ldvr, work, lwork, rwork, info)

        if (info /= 0) then
            write(output_unit, *) "GEEV returned info =", info
            if (info < 0) then
                write(output_unit, *) "The ", -info, "-th argument has illegal value."
            else
                write(output_unit, *) "The QR alg. failed to compute all of the eigenvalues."
                write(output_unit, *) "No eigenvector has been computed."
            endif
            call stop_error("eig: geev error.")
        endif


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

        return
    end subroutine eigh_cdp

    subroutine lstsq_cdp(A, b, x)
        !! Solves a linear least-squares problem \(\min ~ \| \mathbf{Ax} - \mathbf{b} \|_2^2 \) using LAPACK.
        complex(dp), intent(in) :: A(:, :)
        !! Matrix to be "pseudo-inversed".
        complex(dp), intent(in) :: b(:)
        !! Right-hand side vector.
        complex(dp), intent(out) :: x(:)
        !! Solution of the least-squares problem.

        ! Internal variables.
        character :: trans="n"
        integer :: m, n, nrhs, lda, ldb, lwork, info
        complex(dp) :: A_tilde(size(A, 1), size(A, 2)), b_tilde(size(A, 1), 1)
        complex(dp), allocatable :: work(:)

        ! Setup variables.
        m = size(A, 1) ; n = size(A, 2) ; nrhs = 1
        lda = m ; ldb = m ; lwork = max(1, min(m, n) + max(min(m, n), nrhs))
        a_tilde = a ; b_tilde(:, 1) = b
        allocate(work(lwork)) ; work = 0.0_dp

        ! Solve the least-squares problem.
        call gels(trans, m, n, nrhs, a_tilde, lda, b_tilde, ldb, work, lwork, info)

        if (info /= 0) then
            write(output_unit, *) "The ", -info, "-th argument has illegal value."
            call stop_error("lstsq: dgels error")
        endif

        ! Return solution.
        x = b_tilde(1:n, 1)

        return
    end subroutine lstsq_cdp

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

        return
    end subroutine ordschur_cdp


    subroutine sqrtm_rsp(X, sqrtmX)
      !! Matrix-valued sqrt function for dense SPD matrices
      real(sp), intent(in)  :: X(:,:)
      !! Matrix of which to compute the sqrt
      real(sp), intent(out) :: sqrtmX(size(X,1),size(X,1))
      !! Return matrix

      ! internals
      real(sp) :: lambda(size(X,1))
      real(sp) :: V(size(X,1), size(X,1))
      logical :: symmetric
      integer :: i

      ! Check if the matrix is symmetric
      if (.not. is_symmetric(X)) then
        write(output_unit,*) "Error: Input matrix is not symmetric"
        STOP
      end if

      ! Perform eigenvalue decomposition
      call eigh(X, V, lambda)

      ! Check if the matrix is positive definite (up to tol)
      do i = 1, size(lambda)
         if (abs(lambda(i)) .gt. 10*atol_sp ) then
            if (lambda(i) .gt. zero_rsp) then
               lambda(i) = sqrt(lambda(i))
            else
               write(output_unit,*) "Error: Input matrix is not positive definite to tolerance"
               STOP
            end if
         else
            lambda(i) = sqrt(abs(lambda(i)))
            write(output_unit,*) "Warning: Input matrix is singular to tolerance"
         end if
      end do

      ! Reconstruct the square root matrix
      sqrtmX = matmul(V, matmul(diag(lambda), transpose(V)))  

      return
    end subroutine
    subroutine sqrtm_rdp(X, sqrtmX)
      !! Matrix-valued sqrt function for dense SPD matrices
      real(dp), intent(in)  :: X(:,:)
      !! Matrix of which to compute the sqrt
      real(dp), intent(out) :: sqrtmX(size(X,1),size(X,1))
      !! Return matrix

      ! internals
      real(dp) :: lambda(size(X,1))
      real(dp) :: V(size(X,1), size(X,1))
      logical :: symmetric
      integer :: i

      ! Check if the matrix is symmetric
      if (.not. is_symmetric(X)) then
        write(output_unit,*) "Error: Input matrix is not symmetric"
        STOP
      end if

      ! Perform eigenvalue decomposition
      call eigh(X, V, lambda)

      ! Check if the matrix is positive definite (up to tol)
      do i = 1, size(lambda)
         if (abs(lambda(i)) .gt. 10*atol_dp ) then
            if (lambda(i) .gt. zero_rdp) then
               lambda(i) = sqrt(lambda(i))
            else
               write(output_unit,*) "Error: Input matrix is not positive definite to tolerance"
               STOP
            end if
         else
            lambda(i) = sqrt(abs(lambda(i)))
            write(output_unit,*) "Warning: Input matrix is singular to tolerance"
         end if
      end do

      ! Reconstruct the square root matrix
      sqrtmX = matmul(V, matmul(diag(lambda), transpose(V)))  

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

        norm = 0.0_sp
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

        norm = 0.0_dp
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

        norm = 0.0_sp
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

        norm = 0.0_dp
        n = size(A, 1)
        do i = 1, n
            row_sum = sum(abs(A(i, :)))
            norm = max(norm, row_sum)
        enddo
    end function


end module lightkrylov_utils
