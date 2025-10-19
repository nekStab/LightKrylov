submodule(lightkrylov_utils) utility_functions
    !-------------------------------
    !-----     LightKrylov     -----
    !-------------------------------
    use stdlib_optval, only: optval
    use stdlib_linalg_constants, only: ilp
    use stdlib_linalg_lapack, only: geev, trsen, lartg, lasr
    use stdlib_linalg, only: hermitian, svd, diag, eye, mnorm, inv, norm
    use LightKrylov_Timing, only: timer => global_lightkrylov_timer, time_lightkrylov

    implicit none(type, external)
contains

    module procedure log2_rsp
        y = log(x) / log(2.0_sp)
    end procedure log2_rsp
    module procedure log2_rdp
        y = log(x) / log(2.0_dp)
    end procedure log2_rdp

    !---------------------------------------------
    !-----     Shape Assertion Utilities     -----
    !---------------------------------------------

    module procedure assert_shape_vector_rsp
        if (any(shape(v) /= size)) then
            write(output_unit, *) "Vector "//vecname//" has illegal length", shape(v), &
                                        & ". Expected length is ", size, ". Aborting due to illegal vector length."
            call stop_error("Vector length assertion error", module, procedure)
        endif
    end procedure assert_shape_vector_rsp

    module procedure assert_shape_matrix_rsp
        if (any(shape(A) /= size)) then
            write(output_unit, *) "Matrix "//matname//" has illegal shape", shape(A), &
                                        & ". Expected shape is ", size, ". Aborting due to illegal matrix shape."
            call stop_error("Matrix shape assertion error", module, procedure)
        endif
    end procedure assert_shape_matrix_rsp

    module procedure assert_shape_vector_rdp
        if (any(shape(v) /= size)) then
            write(output_unit, *) "Vector "//vecname//" has illegal length", shape(v), &
                                        & ". Expected length is ", size, ". Aborting due to illegal vector length."
            call stop_error("Vector length assertion error", module, procedure)
        endif
    end procedure assert_shape_vector_rdp

    module procedure assert_shape_matrix_rdp
        if (any(shape(A) /= size)) then
            write(output_unit, *) "Matrix "//matname//" has illegal shape", shape(A), &
                                        & ". Expected shape is ", size, ". Aborting due to illegal matrix shape."
            call stop_error("Matrix shape assertion error", module, procedure)
        endif
    end procedure assert_shape_matrix_rdp

    module procedure assert_shape_vector_csp
        if (any(shape(v) /= size)) then
            write(output_unit, *) "Vector "//vecname//" has illegal length", shape(v), &
                                        & ". Expected length is ", size, ". Aborting due to illegal vector length."
            call stop_error("Vector length assertion error", module, procedure)
        endif
    end procedure assert_shape_vector_csp

    module procedure assert_shape_matrix_csp
        if (any(shape(A) /= size)) then
            write(output_unit, *) "Matrix "//matname//" has illegal shape", shape(A), &
                                        & ". Expected shape is ", size, ". Aborting due to illegal matrix shape."
            call stop_error("Matrix shape assertion error", module, procedure)
        endif
    end procedure assert_shape_matrix_csp

    module procedure assert_shape_vector_cdp
        if (any(shape(v) /= size)) then
            write(output_unit, *) "Vector "//vecname//" has illegal length", shape(v), &
                                        & ". Expected length is ", size, ". Aborting due to illegal vector length."
            call stop_error("Vector length assertion error", module, procedure)
        endif
    end procedure assert_shape_vector_cdp

    module procedure assert_shape_matrix_cdp
        if (any(shape(A) /= size)) then
            write(output_unit, *) "Matrix "//matname//" has illegal shape", shape(A), &
                                        & ". Expected shape is ", size, ". Aborting due to illegal matrix shape."
            call stop_error("Matrix shape assertion error", module, procedure)
        endif
    end procedure assert_shape_matrix_cdp

    !--------------------------------------------
    !-----     Linear Algebra Utilities     -----
    !--------------------------------------------

    !----- Eigenvalue Decomposition -----

    module procedure eig_rsp
        character(len=*), parameter :: this_procedure = 'eig_rsp'
        ! Lapack variables.
        character, parameter :: jobvl = "n", jobvr = "v"
        integer(ilp) :: n, lwork, info, lda, ldvl, ldvr
        real(sp) :: A_tilde(size(A, 1), size(A, 2)), vl(1, size(A, 2))
        real(sp) :: work(4*size(A, 1)), wr(size(A, 1)), wi(size(A, 1))

        ! Setup variables.
        n = size(A, 1) ; lda = n ; ldvl = 1 ; ldvr = n ; a_tilde = a
        lwork =  4*n 

        ! Eigendecomposition.
        call geev(jobvl, jobvr, n, a_tilde, lda, wr, wi, vl, ldvl, vecs, ldvr, work, lwork, info)
        call check_info(info, "GEEV", this_module, "eig_rsp")

        ! Complex eigenvalues.
        vals = one_csp*wr + one_im_csp*wi
    end procedure eig_rsp

    module procedure eig_rdp
        character(len=*), parameter :: this_procedure = 'eig_rdp'
        ! Lapack variables.
        character, parameter :: jobvl = "n", jobvr = "v"
        integer(ilp) :: n, lwork, info, lda, ldvl, ldvr
        real(dp) :: A_tilde(size(A, 1), size(A, 2)), vl(1, size(A, 2))
        real(dp) :: work(4*size(A, 1)), wr(size(A, 1)), wi(size(A, 1))

        ! Setup variables.
        n = size(A, 1) ; lda = n ; ldvl = 1 ; ldvr = n ; a_tilde = a
        lwork =  4*n 

        ! Eigendecomposition.
        call geev(jobvl, jobvr, n, a_tilde, lda, wr, wi, vl, ldvl, vecs, ldvr, work, lwork, info)
        call check_info(info, "GEEV", this_module, "eig_rdp")

        ! Complex eigenvalues.
        vals = one_cdp*wr + one_im_cdp*wi
    end procedure eig_rdp

    module procedure eig_csp
        character(len=*), parameter :: this_procedure = 'eig_csp'
        ! Lapack variables.
        character, parameter :: jobvl = "n", jobvr = "v"
        integer(ilp) :: n, lwork, info, lda, ldvl, ldvr
        complex(sp) :: A_tilde(size(A, 1), size(A, 2)), vl(1, size(A, 2))
        complex(sp) :: work(2*size(A, 1))
        real(sp) :: rwork(2*size(A, 1))

        ! Setup variables.
        n = size(A, 1) ; lda = n ; ldvl = 1 ; ldvr = n ; a_tilde = a
        lwork =  2*n 

        ! Eigendecomposition.
        call geev(jobvl, jobvr, n, a_tilde, lda, vals, vl, ldvl, vecs, ldvr, work, lwork, rwork, &
            & info)
        call check_info(info, "GEEV", this_module, "eig_csp")

    end procedure eig_csp

    module procedure eig_cdp
        character(len=*), parameter :: this_procedure = 'eig_cdp'
        ! Lapack variables.
        character, parameter :: jobvl = "n", jobvr = "v"
        integer(ilp) :: n, lwork, info, lda, ldvl, ldvr
        complex(dp) :: A_tilde(size(A, 1), size(A, 2)), vl(1, size(A, 2))
        complex(dp) :: work(2*size(A, 1))
        real(dp) :: rwork(2*size(A, 1))

        ! Setup variables.
        n = size(A, 1) ; lda = n ; ldvl = 1 ; ldvr = n ; a_tilde = a
        lwork =  2*n 

        ! Eigendecomposition.
        call geev(jobvl, jobvr, n, a_tilde, lda, vals, vl, ldvl, vecs, ldvr, work, lwork, rwork, &
            & info)
        call check_info(info, "GEEV", this_module, "eig_cdp")

    end procedure eig_cdp

    !----- Schur Factorization ------

    !----- OrdSchur Factorization -----
    module procedure ordschur_rsp
        ! Lapack variables.
        character, parameter :: job="n", compq="v"
        integer(ilp) :: info, ldq, ldt, lwork, m, n
        real(sp) :: s, sep
        integer(ilp) :: iwork(size(T, 1)), liwork
        real(sp) :: wi(size(T, 1)), wr(size(T, 1)), work(size(T, 1))

        ! Setup variables.
        n = size(T, 2) ; ldt = n ; ldq = n ; lwork = max(1, n)

        if (time_lightkrylov()) call timer%start('trsen')
        liwork = 1
        call trsen(job, compq, selected, n, T, ldt, Q, ldq, wr, wi, m, s, sep, work, lwork, &
               & iwork, liwork, info)
        call check_info(info, "TRSEN", this_module, "ordschur_rsp")
        if (time_lightkrylov()) call timer%stop('trsen')
    end procedure ordschur_rsp

    module procedure ordschur_rdp
        ! Lapack variables.
        character, parameter :: job="n", compq="v"
        integer(ilp) :: info, ldq, ldt, lwork, m, n
        real(dp) :: s, sep
        integer(ilp) :: iwork(size(T, 1)), liwork
        real(dp) :: wi(size(T, 1)), wr(size(T, 1)), work(size(T, 1))

        ! Setup variables.
        n = size(T, 2) ; ldt = n ; ldq = n ; lwork = max(1, n)

        if (time_lightkrylov()) call timer%start('trsen')
        liwork = 1
        call trsen(job, compq, selected, n, T, ldt, Q, ldq, wr, wi, m, s, sep, work, lwork, &
               & iwork, liwork, info)
        call check_info(info, "TRSEN", this_module, "ordschur_rdp")
        if (time_lightkrylov()) call timer%stop('trsen')
    end procedure ordschur_rdp

    module procedure ordschur_csp
        ! Lapack variables.
        character, parameter :: job="n", compq="v"
        integer(ilp) :: info, ldq, ldt, lwork, m, n
        real(sp) :: s, sep
        complex(sp) :: w(size(T, 1)), work(size(T, 1))

        ! Setup variables.
        n = size(T, 2) ; ldt = n ; ldq = n ; lwork = max(1, n)

        if (time_lightkrylov()) call timer%start('trsen')
        call trsen(job, compq, selected, n, T, ldt, Q, ldq, w, m, s, sep, work, lwork, info)
        call check_info(info, "TRSEN", this_module, "ordschur_csp")
        if (time_lightkrylov()) call timer%stop('trsen')
    end procedure ordschur_csp

    module procedure ordschur_cdp
        ! Lapack variables.
        character, parameter :: job="n", compq="v"
        integer(ilp) :: info, ldq, ldt, lwork, m, n
        real(dp) :: s, sep
        complex(dp) :: w(size(T, 1)), work(size(T, 1))

        ! Setup variables.
        n = size(T, 2) ; ldt = n ; ldq = n ; lwork = max(1, n)

        if (time_lightkrylov()) call timer%start('trsen')
        call trsen(job, compq, selected, n, T, ldt, Q, ldq, w, m, s, sep, work, lwork, info)
        call check_info(info, "TRSEN", this_module, "ordschur_cdp")
        if (time_lightkrylov()) call timer%stop('trsen')
    end procedure ordschur_cdp

    !----- Matrix Square-Root -----

    module procedure sqrtm_rsp
        character(len=*), parameter :: this_procedure = 'sqrtm_rsp'
        ! Singular value decomposition.
        real(sp) :: S(size(A, 1))
        real(sp) :: U(size(A, 1), size(A, 1)), UT(size(A, 1), size(A, 1))
        integer(ilp) :: i
        real(sp) :: symmetry_error
        character(len=256) :: msg
        if (time_lightkrylov()) call timer%start(this_procedure)

        info = 0
        ! Symmetry error.
        symmetry_error = 0.5_sp * maxval( abs(A - hermitian(A)) )
        if (symmetry_error > rtol_sp) then
            write(msg, "(2(A,E9.2))") "Input matrix is not Hermitian. 0.5*max(A - A.H) =", &
                & symmetry_error, ", tol = ", rtol_sp
            call stop_error(msg, this_module, "sqrtm_rsp")
        else if (symmetry_error > 10*atol_sp) then
            write(msg, "(A, E9.2)") "Input matrix is not exactly Hermitian. 0.5*max(A - A.H) =", &
                & symmetry_error
            call log_warning(msg, this_module, "sqrtm_rsp")
        endif

        ! Perform SVD.
        if (time_lightkrylov()) call timer%start('svd')
        call svd(A, S, U, UT)
        if (time_lightkrylov()) call timer%stop('svd')

        ! Check if matrix is pos. def. (up to tol).
        do i = 1, size(S)
            if (S(i) > 10*atol_sp) then
                S(i) = sqrt(S(i))
            else
                S(i) = zero_rsp ; info = 1
            endif
        enddo

        ! Reconstruct the square root matrix.
        sqrtA = matmul(U, matmul(diag(S), hermitian(U)))
        if (time_lightkrylov()) call timer%stop(this_procedure)
    end procedure sqrtm_rsp

    module procedure sqrtm_rdp
        character(len=*), parameter :: this_procedure = 'sqrtm_rdp'
        ! Singular value decomposition.
        real(dp) :: S(size(A, 1))
        real(dp) :: U(size(A, 1), size(A, 1)), UT(size(A, 1), size(A, 1))
        integer(ilp) :: i
        real(dp) :: symmetry_error
        character(len=256) :: msg
        if (time_lightkrylov()) call timer%start(this_procedure)

        info = 0
        ! Symmetry error.
        symmetry_error = 0.5_dp * maxval( abs(A - hermitian(A)) )
        if (symmetry_error > rtol_dp) then
            write(msg, "(2(A,E9.2))") "Input matrix is not Hermitian. 0.5*max(A - A.H) =", &
                & symmetry_error, ", tol = ", rtol_dp
            call stop_error(msg, this_module, "sqrtm_rdp")
        else if (symmetry_error > 10*atol_dp) then
            write(msg, "(A, E9.2)") "Input matrix is not exactly Hermitian. 0.5*max(A - A.H) =", &
                & symmetry_error
            call log_warning(msg, this_module, "sqrtm_rdp")
        endif

        ! Perform SVD.
        if (time_lightkrylov()) call timer%start('svd')
        call svd(A, S, U, UT)
        if (time_lightkrylov()) call timer%stop('svd')

        ! Check if matrix is pos. def. (up to tol).
        do i = 1, size(S)
            if (S(i) > 10*atol_dp) then
                S(i) = sqrt(S(i))
            else
                S(i) = zero_rdp ; info = 1
            endif
        enddo

        ! Reconstruct the square root matrix.
        sqrtA = matmul(U, matmul(diag(S), hermitian(U)))
        if (time_lightkrylov()) call timer%stop(this_procedure)
    end procedure sqrtm_rdp

    module procedure sqrtm_csp
        character(len=*), parameter :: this_procedure = 'sqrtm_csp'
        ! Singular value decomposition.
        real(sp) :: S(size(A, 1))
        complex(sp) :: U(size(A, 1), size(A, 1)), UT(size(A, 1), size(A, 1))
        integer(ilp) :: i
        real(sp) :: symmetry_error
        character(len=256) :: msg
        if (time_lightkrylov()) call timer%start(this_procedure)

        info = 0
        ! Symmetry error.
        symmetry_error = 0.5_sp * maxval( abs(A - hermitian(A)) )
        if (symmetry_error > rtol_sp) then
            write(msg, "(2(A,E9.2))") "Input matrix is not Hermitian. 0.5*max(A - A.H) =", &
                & symmetry_error, ", tol = ", rtol_sp
            call stop_error(msg, this_module, "sqrtm_csp")
        else if (symmetry_error > 10*atol_sp) then
            write(msg, "(A, E9.2)") "Input matrix is not exactly Hermitian. 0.5*max(A - A.H) =", &
                & symmetry_error
            call log_warning(msg, this_module, "sqrtm_csp")
        endif

        ! Perform SVD.
        if (time_lightkrylov()) call timer%start('svd')
        call svd(A, S, U, UT)
        if (time_lightkrylov()) call timer%stop('svd')

        ! Check if matrix is pos. def. (up to tol).
        do i = 1, size(S)
            if (S(i) > 10*atol_sp) then
                S(i) = sqrt(S(i))
            else
                S(i) = zero_rsp ; info = 1
            endif
        enddo

        ! Reconstruct the square root matrix.
        sqrtA = matmul(U, matmul(diag(S), hermitian(U)))
        if (time_lightkrylov()) call timer%stop(this_procedure)
    end procedure sqrtm_csp

    module procedure sqrtm_cdp
        character(len=*), parameter :: this_procedure = 'sqrtm_cdp'
        ! Singular value decomposition.
        real(dp) :: S(size(A, 1))
        complex(dp) :: U(size(A, 1), size(A, 1)), UT(size(A, 1), size(A, 1))
        integer(ilp) :: i
        real(dp) :: symmetry_error
        character(len=256) :: msg
        if (time_lightkrylov()) call timer%start(this_procedure)

        info = 0
        ! Symmetry error.
        symmetry_error = 0.5_dp * maxval( abs(A - hermitian(A)) )
        if (symmetry_error > rtol_dp) then
            write(msg, "(2(A,E9.2))") "Input matrix is not Hermitian. 0.5*max(A - A.H) =", &
                & symmetry_error, ", tol = ", rtol_dp
            call stop_error(msg, this_module, "sqrtm_cdp")
        else if (symmetry_error > 10*atol_dp) then
            write(msg, "(A, E9.2)") "Input matrix is not exactly Hermitian. 0.5*max(A - A.H) =", &
                & symmetry_error
            call log_warning(msg, this_module, "sqrtm_cdp")
        endif

        ! Perform SVD.
        if (time_lightkrylov()) call timer%start('svd')
        call svd(A, S, U, UT)
        if (time_lightkrylov()) call timer%stop('svd')

        ! Check if matrix is pos. def. (up to tol).
        do i = 1, size(S)
            if (S(i) > 10*atol_dp) then
                S(i) = sqrt(S(i))
            else
                S(i) = zero_rdp ; info = 1
            endif
        enddo

        ! Reconstruct the square root matrix.
        sqrtA = matmul(U, matmul(diag(S), hermitian(U)))
        if (time_lightkrylov()) call timer%stop(this_procedure)
    end procedure sqrtm_cdp

    !----- Givens rotations -----

    module procedure givens_rotation_rsp
        g = x / norm(x, 2)
    end procedure givens_rotation_rsp

    module procedure apply_givens_rotation_rsp
        integer(ilp) :: k
        real(sp) :: r
        real(sp), pointer :: hmat(:, :)
        !> Size of the column.
        k = size(h) - 1
        !> Apply previous Givens rotations to this new column.
        hmat(1:k, 1:1) => h(:k)
        call lasr("L", "V", "F", k, 1, c(:k-1), s(:k-1), hmat, k)
        !> Compute the sine and cosine compoennts for the next rotation.
        call lartg(h(k), h(k+1), c(k), s(k), r)
        !> Eliminiate H(k+1, k).
        h(k) = r; h(k+1) = 0.0_sp
    end procedure apply_givens_rotation_rsp
    module procedure givens_rotation_rdp
        g = x / norm(x, 2)
    end procedure givens_rotation_rdp

    module procedure apply_givens_rotation_rdp
        integer(ilp) :: k
        real(dp) :: r
        real(dp), pointer :: hmat(:, :)
        !> Size of the column.
        k = size(h) - 1
        !> Apply previous Givens rotations to this new column.
        hmat(1:k, 1:1) => h(:k)
        call lasr("L", "V", "F", k, 1, c(:k-1), s(:k-1), hmat, k)
        !> Compute the sine and cosine compoennts for the next rotation.
        call lartg(h(k), h(k+1), c(k), s(k), r)
        !> Eliminiate H(k+1, k).
        h(k) = r; h(k+1) = 0.0_dp
    end procedure apply_givens_rotation_rdp
    module procedure givens_rotation_csp
        g = x / norm(x, 2)
    end procedure givens_rotation_csp

    module procedure apply_givens_rotation_csp
        integer(ilp) :: i, k
        complex(sp) :: t, r, g(2)
        !> Size of the column.
        k = size(h) - 1
        !> Apply previous Givens rotations to this new column.
        do i = 1, k-1
            t = c(i)*h(i) + s(i)*h(i+1)
            h(i+1) = -s(i)*h(i) + c(i)*h(i+1)
            h(i) = t
        enddo
        !> Compute the sine and cosine compoennts for the next rotation.
        g = givens_rotation([h(k), h(k+1)]) ; c(k) = g(1) ; s(k) = g(2)
        !> Eliminiate H(k+1, k).
        h(k) = c(k)*h(k) + s(k)*h(k+1) ; h(k+1) = 0.0_sp

    end procedure apply_givens_rotation_csp
    module procedure givens_rotation_cdp
        g = x / norm(x, 2)
    end procedure givens_rotation_cdp

    module procedure apply_givens_rotation_cdp
        integer(ilp) :: i, k
        complex(dp) :: t, r, g(2)
        !> Size of the column.
        k = size(h) - 1
        !> Apply previous Givens rotations to this new column.
        do i = 1, k-1
            t = c(i)*h(i) + s(i)*h(i+1)
            h(i+1) = -s(i)*h(i) + c(i)*h(i+1)
            h(i) = t
        enddo
        !> Compute the sine and cosine compoennts for the next rotation.
        g = givens_rotation([h(k), h(k+1)]) ; c(k) = g(1) ; s(k) = g(2)
        !> Eliminiate H(k+1, k).
        h(k) = c(k)*h(k) + s(k)*h(k+1) ; h(k+1) = 0.0_dp

    end procedure apply_givens_rotation_cdp
end submodule utility_functions

