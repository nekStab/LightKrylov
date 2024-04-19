module lightkrylov_utils
    use iso_fortran_env, only: output_unit
    use stdlib_linalg_lapack, only : getrf, getri, gesvd, geev, syev, gels
    use lightkrylov_constants
    implicit none

    private

    public :: stop_error
    public :: assert_shape
    public :: inv
    public :: svd
    public :: eig
    public :: lstsq

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

    interface lstsq
        module procedure lstsq_rsp
        module procedure lstsq_rdp
        module procedure lstsq_csp
        module procedure lstsq_cdp
    end interface

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

    subroutine lstsq_rsp(A, b, x)
        !! Solves a linear least-squares problem \(\min ~ \| \mathbf{Ax} - \mathbf{b} \|_2^2 \) using LAPACK.
        real(sp), intent(in) :: A(:, :)
        !! Matrix to be "pseudo-inversed".
        real(sp), intent(in) :: b
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
        a_tilde = a ; b_tilde = b
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

    subroutine lstsq_rdp(A, b, x)
        !! Solves a linear least-squares problem \(\min ~ \| \mathbf{Ax} - \mathbf{b} \|_2^2 \) using LAPACK.
        real(dp), intent(in) :: A(:, :)
        !! Matrix to be "pseudo-inversed".
        real(dp), intent(in) :: b
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
        a_tilde = a ; b_tilde = b
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

    subroutine lstsq_csp(A, b, x)
        !! Solves a linear least-squares problem \(\min ~ \| \mathbf{Ax} - \mathbf{b} \|_2^2 \) using LAPACK.
        complex(sp), intent(in) :: A(:, :)
        !! Matrix to be "pseudo-inversed".
        complex(sp), intent(in) :: b
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
        a_tilde = a ; b_tilde = b
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

    subroutine lstsq_cdp(A, b, x)
        !! Solves a linear least-squares problem \(\min ~ \| \mathbf{Ax} - \mathbf{b} \|_2^2 \) using LAPACK.
        complex(dp), intent(in) :: A(:, :)
        !! Matrix to be "pseudo-inversed".
        complex(dp), intent(in) :: b
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
        a_tilde = a ; b_tilde = b
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


end module lightkrylov_utils
