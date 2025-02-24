submodule(lightkrylov_utils) utility_functions
    !-------------------------------
    !-----     LightKrylov     -----
    !-------------------------------
    use stdlib_optval, only: optval
    use stdlib_linalg_constants, only: ilp, lk
    use stdlib_linalg_lapack, only: geev, trsen
    use stdlib_linalg, only: hermitian, svd, diag, eye, mnorm, inv, norm

    implicit none(type, external)
contains

    module procedure log2_rsp
        y = log(x) / log(2.0_sp)
    end procedure
    module procedure log2_rdp
        y = log(x) / log(2.0_dp)
    end procedure

    !---------------------------------------------
    !-----     Shape Assertion Utilities     -----
    !---------------------------------------------

    module procedure assert_shape_vector_rsp
        if (any(shape(v) /= size)) then
            write(output_unit, *) "Vector "//vecname//" has illegal length", shape(v), &
                                        & ". Expected length is ", size, ". Aborting due to illegal vector length."
            call stop_error("Vector length assertion error", module, procedure)
        endif
    end procedure

    module procedure assert_shape_matrix_rsp
        if (any(shape(A) /= size)) then
            write(output_unit, *) "Matrix "//matname//" has illegal shape", shape(A), &
                                        & ". Expected shape is ", size, ". Aborting due to illegal matrix shape."
            call stop_error("Matrix shape assertion error", module, procedure)
        endif
    end procedure
    module procedure assert_shape_vector_rdp
        if (any(shape(v) /= size)) then
            write(output_unit, *) "Vector "//vecname//" has illegal length", shape(v), &
                                        & ". Expected length is ", size, ". Aborting due to illegal vector length."
            call stop_error("Vector length assertion error", module, procedure)
        endif
    end procedure

    module procedure assert_shape_matrix_rdp
        if (any(shape(A) /= size)) then
            write(output_unit, *) "Matrix "//matname//" has illegal shape", shape(A), &
                                        & ". Expected shape is ", size, ". Aborting due to illegal matrix shape."
            call stop_error("Matrix shape assertion error", module, procedure)
        endif
    end procedure
    module procedure assert_shape_vector_csp
        if (any(shape(v) /= size)) then
            write(output_unit, *) "Vector "//vecname//" has illegal length", shape(v), &
                                        & ". Expected length is ", size, ". Aborting due to illegal vector length."
            call stop_error("Vector length assertion error", module, procedure)
        endif
    end procedure

    module procedure assert_shape_matrix_csp
        if (any(shape(A) /= size)) then
            write(output_unit, *) "Matrix "//matname//" has illegal shape", shape(A), &
                                        & ". Expected shape is ", size, ". Aborting due to illegal matrix shape."
            call stop_error("Matrix shape assertion error", module, procedure)
        endif
    end procedure
    module procedure assert_shape_vector_cdp
        if (any(shape(v) /= size)) then
            write(output_unit, *) "Vector "//vecname//" has illegal length", shape(v), &
                                        & ". Expected length is ", size, ". Aborting due to illegal vector length."
            call stop_error("Vector length assertion error", module, procedure)
        endif
    end procedure

    module procedure assert_shape_matrix_cdp
        if (any(shape(A) /= size)) then
            write(output_unit, *) "Matrix "//matname//" has illegal shape", shape(A), &
                                        & ". Expected shape is ", size, ". Aborting due to illegal matrix shape."
            call stop_error("Matrix shape assertion error", module, procedure)
        endif
    end procedure

    !--------------------------------------------
    !-----     Linear Algebra Utilities     -----
    !--------------------------------------------

    !----- Eigenvalue Decomposition -----

    module procedure eig_rsp
        ! Lapack variables.
        character :: jobvl = "n", jobvr = "v"
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
    end procedure
    module procedure eig_rdp
        ! Lapack variables.
        character :: jobvl = "n", jobvr = "v"
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
    end procedure
    module procedure eig_csp
        ! Lapack variables.
        character :: jobvl = "n", jobvr = "v"
        integer(ilp) :: n, lwork, info, lda, ldvl, ldvr
        complex(sp) :: A_tilde(size(A, 1), size(A, 2)), vl(1, size(A, 2))
        complex(sp) :: work(2*size(A, 1))
        real(sp) :: rwork(2*size(A, 1))

        ! Setup variables.
        n = size(A, 1) ; lda = n ; ldvl = 1 ; ldvr = n ; a_tilde = a
        lwork =  2*n 

        ! Eigendecomposition.
        call geev(jobvl, jobvr, n, a_tilde, lda, vals, vl, ldvl, vecs, ldvr, work, lwork, rwork, info)
        call check_info(info, "GEEV", this_module, "eig_csp")

    end procedure
    module procedure eig_cdp
        ! Lapack variables.
        character :: jobvl = "n", jobvr = "v"
        integer(ilp) :: n, lwork, info, lda, ldvl, ldvr
        complex(dp) :: A_tilde(size(A, 1), size(A, 2)), vl(1, size(A, 2))
        complex(dp) :: work(2*size(A, 1))
        real(dp) :: rwork(2*size(A, 1))

        ! Setup variables.
        n = size(A, 1) ; lda = n ; ldvl = 1 ; ldvr = n ; a_tilde = a
        lwork =  2*n 

        ! Eigendecomposition.
        call geev(jobvl, jobvr, n, a_tilde, lda, vals, vl, ldvl, vecs, ldvr, work, lwork, rwork, info)
        call check_info(info, "GEEV", this_module, "eig_cdp")

    end procedure

    !----- OrdSchur Factorization -----
    module procedure ordschur_rsp
        ! Lapack variables.
        character :: job="n", compq="v"
        integer(ilp) :: info, ldq, ldt, lwork, m, n
        real(sp) :: s, sep
        integer(ilp) :: iwork(size(T, 1)), liwork
        real(sp) :: wi(size(T, 1)), wr(size(T, 1)), work(size(T, 1))

        ! Setup variables.
        n = size(T, 2) ; ldt = n ; ldq = n ; lwork = max(1, n)

        liwork = 1
        call trsen(job, compq, selected, n, T, ldt, Q, ldq, wr, wi, m, s, sep, work, lwork, iwork, liwork, info)
        call check_info(info, "TRSEN", this_module, "ordschur_rsp")
    end procedure
    module procedure ordschur_rdp
        ! Lapack variables.
        character :: job="n", compq="v"
        integer(ilp) :: info, ldq, ldt, lwork, m, n
        real(dp) :: s, sep
        integer(ilp) :: iwork(size(T, 1)), liwork
        real(dp) :: wi(size(T, 1)), wr(size(T, 1)), work(size(T, 1))

        ! Setup variables.
        n = size(T, 2) ; ldt = n ; ldq = n ; lwork = max(1, n)

        liwork = 1
        call trsen(job, compq, selected, n, T, ldt, Q, ldq, wr, wi, m, s, sep, work, lwork, iwork, liwork, info)
        call check_info(info, "TRSEN", this_module, "ordschur_rdp")
    end procedure
    module procedure ordschur_csp
        ! Lapack variables.
        character :: job="n", compq="v"
        integer(ilp) :: info, ldq, ldt, lwork, m, n
        real(sp) :: s, sep
        complex(sp) :: w(size(T, 1)), work(size(T, 1))

        ! Setup variables.
        n = size(T, 2) ; ldt = n ; ldq = n ; lwork = max(1, n)

        call trsen(job, compq, selected, n, T, ldt, Q, ldq, w, m, s, sep, work, lwork, info)
        call check_info(info, "TRSEN", this_module, "ordschur_csp")
    end procedure
    module procedure ordschur_cdp
        ! Lapack variables.
        character :: job="n", compq="v"
        integer(ilp) :: info, ldq, ldt, lwork, m, n
        real(dp) :: s, sep
        complex(dp) :: w(size(T, 1)), work(size(T, 1))

        ! Setup variables.
        n = size(T, 2) ; ldt = n ; ldq = n ; lwork = max(1, n)

        call trsen(job, compq, selected, n, T, ldt, Q, ldq, w, m, s, sep, work, lwork, info)
        call check_info(info, "TRSEN", this_module, "ordschur_cdp")
    end procedure

    !----- Matrix Square-Root -----

    module procedure sqrtm_rsp
        ! Singular value decomposition.
        real(sp) :: S(size(A, 1))
        real(sp) :: U(size(A, 1), size(A, 1)), UT(size(A, 1), size(A, 1))
        integer(ilp) :: i
        real(sp) :: symmetry_error
        character(len=256) :: msg

        info = 0
        ! Symmetry error.
        symmetry_error = 0.5_sp * maxval( abs(A - hermitian(A)) )
        if (symmetry_error > rtol_sp) then
            write(msg, "(2(A,E9.2))") "Input matrix is not Hermitian. 0.5*max(A - A.H) =", &
                & symmetry_error, ", tol = ", rtol_sp
            call stop_error(msg, this_module, "sqrtm_rsp")
        else if (symmetry_error > 10*atol_sp) then
            write(msg, "(A, E9.2)") "Input matrix is not exactly Hermitian. 0.5*max(A - A.H) =", symmetry_error
            call log_warning(msg, this_module, "sqrtm_rsp")
        endif

        ! Perform SVD.
        call svd(A, S, U, UT)

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
    end procedure
    module procedure sqrtm_rdp
        ! Singular value decomposition.
        real(dp) :: S(size(A, 1))
        real(dp) :: U(size(A, 1), size(A, 1)), UT(size(A, 1), size(A, 1))
        integer(ilp) :: i
        real(dp) :: symmetry_error
        character(len=256) :: msg

        info = 0
        ! Symmetry error.
        symmetry_error = 0.5_dp * maxval( abs(A - hermitian(A)) )
        if (symmetry_error > rtol_dp) then
            write(msg, "(2(A,E9.2))") "Input matrix is not Hermitian. 0.5*max(A - A.H) =", &
                & symmetry_error, ", tol = ", rtol_dp
            call stop_error(msg, this_module, "sqrtm_rdp")
        else if (symmetry_error > 10*atol_dp) then
            write(msg, "(A, E9.2)") "Input matrix is not exactly Hermitian. 0.5*max(A - A.H) =", symmetry_error
            call log_warning(msg, this_module, "sqrtm_rdp")
        endif

        ! Perform SVD.
        call svd(A, S, U, UT)

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
    end procedure
    module procedure sqrtm_csp
        ! Singular value decomposition.
        real(sp) :: S(size(A, 1))
        complex(sp) :: U(size(A, 1), size(A, 1)), UT(size(A, 1), size(A, 1))
        integer(ilp) :: i
        real(sp) :: symmetry_error
        character(len=256) :: msg

        info = 0
        ! Symmetry error.
        symmetry_error = 0.5_sp * maxval( abs(A - hermitian(A)) )
        if (symmetry_error > rtol_sp) then
            write(msg, "(2(A,E9.2))") "Input matrix is not Hermitian. 0.5*max(A - A.H) =", &
                & symmetry_error, ", tol = ", rtol_sp
            call stop_error(msg, this_module, "sqrtm_csp")
        else if (symmetry_error > 10*atol_sp) then
            write(msg, "(A, E9.2)") "Input matrix is not exactly Hermitian. 0.5*max(A - A.H) =", symmetry_error
            call log_warning(msg, this_module, "sqrtm_csp")
        endif

        ! Perform SVD.
        call svd(A, S, U, UT)

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
    end procedure
    module procedure sqrtm_cdp
        ! Singular value decomposition.
        real(dp) :: S(size(A, 1))
        complex(dp) :: U(size(A, 1), size(A, 1)), UT(size(A, 1), size(A, 1))
        integer(ilp) :: i
        real(dp) :: symmetry_error
        character(len=256) :: msg

        info = 0
        ! Symmetry error.
        symmetry_error = 0.5_dp * maxval( abs(A - hermitian(A)) )
        if (symmetry_error > rtol_dp) then
            write(msg, "(2(A,E9.2))") "Input matrix is not Hermitian. 0.5*max(A - A.H) =", &
                & symmetry_error, ", tol = ", rtol_dp
            call stop_error(msg, this_module, "sqrtm_cdp")
        else if (symmetry_error > 10*atol_dp) then
            write(msg, "(A, E9.2)") "Input matrix is not exactly Hermitian. 0.5*max(A - A.H) =", symmetry_error
            call log_warning(msg, this_module, "sqrtm_cdp")
        endif

        ! Perform SVD.
        call svd(A, S, U, UT)

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
    end procedure

    !----- Dense Matrix Exponential -----

    module procedure expm_rsp
        real(sp), allocatable :: A2(:, :), Q(:, :), X(:, :)
        real(sp) :: a_norm, c
        integer :: n, ee, k, s
        logical :: p
        integer :: p_order

        ! Deal with optional args.
        p_order = optval(order, 10)

        n = size(A, 1)

        ! Allocate arrays.
        allocate(A2(n, n)) ; allocate(X(n, n)) ; allocate(Q(n, n))

        ! Compute the L-infinity norm.
        a_norm = mnorm(A, "inf")

        ! Determine scaling factor for the matrix.
        ee = int(log2(a_norm)) + 1
        s = max(0, ee+1)

        ! Scale the input matrix & initialize polynomial.
        A2 = A / 2.0_sp**s
        X = A2

        ! Initialize P & Q and add first step.
        c = 0.5_sp
        E = eye(n, mold=1.0_sp) ; E = E + c*A2

        Q = eye(n, mold=1.0_sp) ; Q = Q - c*A2

        ! Iteratively compute the Pade approximation.
        p = .true.
        do k = 2, p_order
            c = c*(p_order - k + 1) / (k * (2*p_order - k + 1))
            X = matmul(A2, X)
            E = E + c*X
            if (p) then
                Q = Q + c*X
            else
                Q = Q - c*X
            endif
            p = .not. p
        enddo

        E = matmul(inv(Q), E)
        do k = 1, s
            E = matmul(E, E)
        enddo

        return
    end procedure
    module procedure expm_rdp
        real(dp), allocatable :: A2(:, :), Q(:, :), X(:, :)
        real(dp) :: a_norm, c
        integer :: n, ee, k, s
        logical :: p
        integer :: p_order

        ! Deal with optional args.
        p_order = optval(order, 10)

        n = size(A, 1)

        ! Allocate arrays.
        allocate(A2(n, n)) ; allocate(X(n, n)) ; allocate(Q(n, n))

        ! Compute the L-infinity norm.
        a_norm = mnorm(A, "inf")

        ! Determine scaling factor for the matrix.
        ee = int(log2(a_norm)) + 1
        s = max(0, ee+1)

        ! Scale the input matrix & initialize polynomial.
        A2 = A / 2.0_dp**s
        X = A2

        ! Initialize P & Q and add first step.
        c = 0.5_dp
        E = eye(n, mold=1.0_dp) ; E = E + c*A2

        Q = eye(n, mold=1.0_dp) ; Q = Q - c*A2

        ! Iteratively compute the Pade approximation.
        p = .true.
        do k = 2, p_order
            c = c*(p_order - k + 1) / (k * (2*p_order - k + 1))
            X = matmul(A2, X)
            E = E + c*X
            if (p) then
                Q = Q + c*X
            else
                Q = Q - c*X
            endif
            p = .not. p
        enddo

        E = matmul(inv(Q), E)
        do k = 1, s
            E = matmul(E, E)
        enddo

        return
    end procedure
    module procedure expm_csp
        complex(sp), allocatable :: A2(:, :), Q(:, :), X(:, :)
        real(sp) :: a_norm, c
        integer :: n, ee, k, s
        logical :: p
        integer :: p_order

        ! Deal with optional args.
        p_order = optval(order, 10)

        n = size(A, 1)

        ! Allocate arrays.
        allocate(A2(n, n)) ; allocate(X(n, n)) ; allocate(Q(n, n))

        ! Compute the L-infinity norm.
        a_norm = mnorm(A, "inf")

        ! Determine scaling factor for the matrix.
        ee = int(log2(a_norm)) + 1
        s = max(0, ee+1)

        ! Scale the input matrix & initialize polynomial.
        A2 = A / 2.0_sp**s
        X = A2

        ! Initialize P & Q and add first step.
        c = 0.5_sp
        E = eye(n, mold=1.0_sp) ; E = E + c*A2

        Q = eye(n, mold=1.0_sp) ; Q = Q - c*A2

        ! Iteratively compute the Pade approximation.
        p = .true.
        do k = 2, p_order
            c = c*(p_order - k + 1) / (k * (2*p_order - k + 1))
            X = matmul(A2, X)
            E = E + c*X
            if (p) then
                Q = Q + c*X
            else
                Q = Q - c*X
            endif
            p = .not. p
        enddo

        E = matmul(inv(Q), E)
        do k = 1, s
            E = matmul(E, E)
        enddo

        return
    end procedure
    module procedure expm_cdp
        complex(dp), allocatable :: A2(:, :), Q(:, :), X(:, :)
        real(dp) :: a_norm, c
        integer :: n, ee, k, s
        logical :: p
        integer :: p_order

        ! Deal with optional args.
        p_order = optval(order, 10)

        n = size(A, 1)

        ! Allocate arrays.
        allocate(A2(n, n)) ; allocate(X(n, n)) ; allocate(Q(n, n))

        ! Compute the L-infinity norm.
        a_norm = mnorm(A, "inf")

        ! Determine scaling factor for the matrix.
        ee = int(log2(a_norm)) + 1
        s = max(0, ee+1)

        ! Scale the input matrix & initialize polynomial.
        A2 = A / 2.0_dp**s
        X = A2

        ! Initialize P & Q and add first step.
        c = 0.5_dp
        E = eye(n, mold=1.0_dp) ; E = E + c*A2

        Q = eye(n, mold=1.0_dp) ; Q = Q - c*A2

        ! Iteratively compute the Pade approximation.
        p = .true.
        do k = 2, p_order
            c = c*(p_order - k + 1) / (k * (2*p_order - k + 1))
            X = matmul(A2, X)
            E = E + c*X
            if (p) then
                Q = Q + c*X
            else
                Q = Q - c*X
            endif
            p = .not. p
        enddo

        E = matmul(inv(Q), E)
        do k = 1, s
            E = matmul(E, E)
        enddo

        return
    end procedure

    !----- Givens rotations -----

    module procedure givens_rotation_rsp
        g = x / norm(x, 2)
    end procedure

    module procedure apply_givens_rotation_rsp
        integer(ilp) :: i, k
        real(sp) :: t, g(2)

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
    end procedure
    module procedure givens_rotation_rdp
        g = x / norm(x, 2)
    end procedure

    module procedure apply_givens_rotation_rdp
        integer(ilp) :: i, k
        real(dp) :: t, g(2)

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
    end procedure
    module procedure givens_rotation_csp
        g = x / norm(x, 2)
    end procedure

    module procedure apply_givens_rotation_csp
        integer(ilp) :: i, k
        complex(sp) :: t, g(2)

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
    end procedure
    module procedure givens_rotation_cdp
        g = x / norm(x, 2)
    end procedure

    module procedure apply_givens_rotation_cdp
        integer(ilp) :: i, k
        complex(dp) :: t, g(2)

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
    end procedure

    !----- Solving triangular systems -----

    module procedure solve_triangular_rsp
        integer(ilp) :: i, j, n
        !> Problem's dimensions.
        n = size(A, 1) ; allocate(x, mold=b) ; x = 0.0_sp
        !> Back-substitution algorithm.
        x(n) = b(n) / A(n, n)
        do i = n-1, 1, -1
            x(i) = (b(i) - sum(A(i, i+1:) * x(i+1:))) / A(i, i)
        enddo
    end procedure
    module procedure solve_triangular_rdp
        integer(ilp) :: i, j, n
        !> Problem's dimensions.
        n = size(A, 1) ; allocate(x, mold=b) ; x = 0.0_dp
        !> Back-substitution algorithm.
        x(n) = b(n) / A(n, n)
        do i = n-1, 1, -1
            x(i) = (b(i) - sum(A(i, i+1:) * x(i+1:))) / A(i, i)
        enddo
    end procedure
    module procedure solve_triangular_csp
        integer(ilp) :: i, j, n
        !> Problem's dimensions.
        n = size(A, 1) ; allocate(x, mold=b) ; x = 0.0_sp
        !> Back-substitution algorithm.
        x(n) = b(n) / A(n, n)
        do i = n-1, 1, -1
            x(i) = (b(i) - sum(A(i, i+1:) * x(i+1:))) / A(i, i)
        enddo
    end procedure
    module procedure solve_triangular_cdp
        integer(ilp) :: i, j, n
        !> Problem's dimensions.
        n = size(A, 1) ; allocate(x, mold=b) ; x = 0.0_dp
        !> Back-substitution algorithm.
        x(n) = b(n) / A(n, n)
        do i = n-1, 1, -1
            x(i) = (b(i) - sum(A(i, i+1:) * x(i+1:))) / A(i, i)
        enddo
    end procedure

    !--------------------------------------------
    !-----     (Modified) Leja Ordering     -----
    !--------------------------------------------

    module procedure leja_ordering_rsp
        logical(lk) :: modified_
        ! Standard or modified Leja ordering.
        modified_ = optval(modified, .false.)
        ! Dispatch to the appropriate driver.
        y = standard_leja_ordering_rsp(x)
    end procedure
    module procedure leja_ordering_rdp
        logical(lk) :: modified_
        ! Standard or modified Leja ordering.
        modified_ = optval(modified, .false.)
        ! Dispatch to the appropriate driver.
        y = standard_leja_ordering_rdp(x)
    end procedure
    module procedure leja_ordering_csp
        logical(lk) :: modified_
        ! Standard or modified Leja ordering.
        modified_ = optval(modified, .false.)
        ! Dispatch to the appropriate driver.
        if (modified_) then
            y = modified_leja_ordering_csp(x)
        else
            y = standard_leja_ordering_csp(x)
        endif
    end procedure
    module procedure leja_ordering_cdp
        logical(lk) :: modified_
        ! Standard or modified Leja ordering.
        modified_ = optval(modified, .false.)
        ! Dispatch to the appropriate driver.
        if (modified_) then
            y = modified_leja_ordering_cdp(x)
        else
            y = standard_leja_ordering_cdp(x)
        endif
    end procedure

    pure function standard_leja_ordering_rsp(x) result(y)
        real(sp), intent(in) :: x(:)
        !! List of points to be sorted.
        real(sp) :: y(size(x))
        !! Sorted list of points.

        !----- Internal variables -----
        logical(lk) :: selected(size(x))
        integer(ilp) :: i, j, k, idx
        real(sp) :: distances(size(x))

        !> Initialize.
        j = 0 ; selected = .false.

        !> Find the first x with maximum absolute value.
        idx = maxloc(abs(x), 1) ; j = j+1
        y(j) = x(idx) ; selected(idx) = .true.

        !> Loop through all the other roots.
        do while (j < size(x))
            ! Distance to the set of points already selected.
            distances = -huge(1.0_sp)
            do concurrent(k=1:size(x))
                if (.not. selected(k)) distances(k) = product(abs(x(k) - y(:j)))
            enddo
            ! Greedily choose the one with the maximum distance.
            idx = maxloc(distances, 1) ; j = j+1
            y(j) = x(idx) ; selected(idx) = .true.
        enddo
    end function
    pure function standard_leja_ordering_rdp(x) result(y)
        real(dp), intent(in) :: x(:)
        !! List of points to be sorted.
        real(dp) :: y(size(x))
        !! Sorted list of points.

        !----- Internal variables -----
        logical(lk) :: selected(size(x))
        integer(ilp) :: i, j, k, idx
        real(dp) :: distances(size(x))

        !> Initialize.
        j = 0 ; selected = .false.

        !> Find the first x with maximum absolute value.
        idx = maxloc(abs(x), 1) ; j = j+1
        y(j) = x(idx) ; selected(idx) = .true.

        !> Loop through all the other roots.
        do while (j < size(x))
            ! Distance to the set of points already selected.
            distances = -huge(1.0_dp)
            do concurrent(k=1:size(x))
                if (.not. selected(k)) distances(k) = product(abs(x(k) - y(:j)))
            enddo
            ! Greedily choose the one with the maximum distance.
            idx = maxloc(distances, 1) ; j = j+1
            y(j) = x(idx) ; selected(idx) = .true.
        enddo
    end function
    pure function standard_leja_ordering_csp(x) result(y)
        complex(sp), intent(in) :: x(:)
        !! List of points to be sorted.
        complex(sp) :: y(size(x))
        !! Sorted list of points.

        !----- Internal variables -----
        logical(lk) :: selected(size(x))
        integer(ilp) :: i, j, k, idx
        real(sp) :: distances(size(x))

        !> Initialize.
        j = 0 ; selected = .false.

        !> Find the first x with maximum absolute value.
        idx = maxloc(abs(x), 1) ; j = j+1
        y(j) = x(idx) ; selected(idx) = .true.

        !> Loop through all the other roots.
        do while (j < size(x))
            ! Distance to the set of points already selected.
            distances = -huge(1.0_sp)
            do concurrent(k=1:size(x))
                if (.not. selected(k)) distances(k) = product(abs(x(k) - y(:j)))
            enddo
            ! Greedily choose the one with the maximum distance.
            idx = maxloc(distances, 1) ; j = j+1
            y(j) = x(idx) ; selected(idx) = .true.
        enddo
    end function
    pure function standard_leja_ordering_cdp(x) result(y)
        complex(dp), intent(in) :: x(:)
        !! List of points to be sorted.
        complex(dp) :: y(size(x))
        !! Sorted list of points.

        !----- Internal variables -----
        logical(lk) :: selected(size(x))
        integer(ilp) :: i, j, k, idx
        real(dp) :: distances(size(x))

        !> Initialize.
        j = 0 ; selected = .false.

        !> Find the first x with maximum absolute value.
        idx = maxloc(abs(x), 1) ; j = j+1
        y(j) = x(idx) ; selected(idx) = .true.

        !> Loop through all the other roots.
        do while (j < size(x))
            ! Distance to the set of points already selected.
            distances = -huge(1.0_dp)
            do concurrent(k=1:size(x))
                if (.not. selected(k)) distances(k) = product(abs(x(k) - y(:j)))
            enddo
            ! Greedily choose the one with the maximum distance.
            idx = maxloc(distances, 1) ; j = j+1
            y(j) = x(idx) ; selected(idx) = .true.
        enddo
    end function

    pure function modified_leja_ordering_csp(x) result(y)
        complex(sp), intent(in) :: x(:)
        !! List of points to be sorted.
        complex(sp) :: y(size(x))
        !! Sorted list of points.

        !----- Internal variables -----
        logical(lk) :: selected(size(x))
        integer(ilp) :: i, j, k, idx
        real(sp) :: distances(size(x)), re, im

        !> Initialize.
        j = 0 ; selected = .false.

        !> Find the first x with maximum absolute value.
        idx = maxloc(abs(x), 1) ; j = j+1 ; re = real(x(idx)) ; im = imag(x(idx))
        y(j) = x(idx) ; selected(idx) = .true.
        !> Add the complex conjugate if needed.
        if (im /= 0.0_sp) then
            j = j+1
            if (im > 0.0_sp) then
                y(j) = x(idx+1) ; selected(idx+1) = .true.
            else
                y(j) = x(idx-1) ; selected(idx-1) = .true.
                y(j-1:j) = conjg(y(j-1:j))
            endif
        endif

        !> Loop through the remaining roots.
        do while (j < size(x))
            !> Compute distances to the set of already selected points.
            distances = -huge(1.0_sp)
            do concurrent(k=1:size(x))
                if (.not. selected(k)) distances(k) = product(abs(x(k) - y(:j)))
            enddo
            !> Greedily select the next point.
            idx = maxloc(distances, 1) ; j = j+1
            y(j) = x(idx) ; selected(idx) = .true.
            re = real(x(idx)) ; im = imag(x(idx))
            !> Add the complex conjugate if needed.
            if (im /= 0.0_sp) then
                j = j+1
                if (im > 0.0_sp) then
                    y(j) = x(idx+1) ; selected(idx+1) = .true.
                else
                    y(j) = x(idx-1) ; selected(idx-1) = .true.
                    y(j-1:j) = conjg(y(j-1:j))
                endif
            endif
        enddo
    end function
    pure function modified_leja_ordering_cdp(x) result(y)
        complex(dp), intent(in) :: x(:)
        !! List of points to be sorted.
        complex(dp) :: y(size(x))
        !! Sorted list of points.

        !----- Internal variables -----
        logical(lk) :: selected(size(x))
        integer(ilp) :: i, j, k, idx
        real(dp) :: distances(size(x)), re, im

        !> Initialize.
        j = 0 ; selected = .false.

        !> Find the first x with maximum absolute value.
        idx = maxloc(abs(x), 1) ; j = j+1 ; re = real(x(idx)) ; im = imag(x(idx))
        y(j) = x(idx) ; selected(idx) = .true.
        !> Add the complex conjugate if needed.
        if (im /= 0.0_dp) then
            j = j+1
            if (im > 0.0_dp) then
                y(j) = x(idx+1) ; selected(idx+1) = .true.
            else
                y(j) = x(idx-1) ; selected(idx-1) = .true.
                y(j-1:j) = conjg(y(j-1:j))
            endif
        endif

        !> Loop through the remaining roots.
        do while (j < size(x))
            !> Compute distances to the set of already selected points.
            distances = -huge(1.0_dp)
            do concurrent(k=1:size(x))
                if (.not. selected(k)) distances(k) = product(abs(x(k) - y(:j)))
            enddo
            !> Greedily select the next point.
            idx = maxloc(distances, 1) ; j = j+1
            y(j) = x(idx) ; selected(idx) = .true.
            re = real(x(idx)) ; im = imag(x(idx))
            !> Add the complex conjugate if needed.
            if (im /= 0.0_dp) then
                j = j+1
                if (im > 0.0_dp) then
                    y(j) = x(idx+1) ; selected(idx+1) = .true.
                else
                    y(j) = x(idx-1) ; selected(idx-1) = .true.
                    y(j-1:j) = conjg(y(j-1:j))
                endif
            endif
        enddo
    end function

    !-------------------------------------------------------------------
    !-----     (Modified) Leja ordering with stability control     -----
    !-------------------------------------------------------------------

    module procedure stabilized_leja_ordering_rsp
        logical(lk) :: modified_
        real(sp), allocatable :: z(:)
        ! Standard or modified Leja ordering.
        modified_ = optval(modified, .false.)
        ! Perform the (modified) Leja ordering.
        z = leja_ordering(x, modified_)
        ! Dispatch to the appropriate stability control driver.
        y = standard_stability_control_rsp(z)
    end procedure
    module procedure stabilized_leja_ordering_rdp
        logical(lk) :: modified_
        real(dp), allocatable :: z(:)
        ! Standard or modified Leja ordering.
        modified_ = optval(modified, .false.)
        ! Perform the (modified) Leja ordering.
        z = leja_ordering(x, modified_)
        ! Dispatch to the appropriate stability control driver.
        y = standard_stability_control_rdp(z)
    end procedure
    module procedure stabilized_leja_ordering_csp
        logical(lk) :: modified_
        complex(sp), allocatable :: z(:)
        ! Standard or modified Leja ordering.
        modified_ = optval(modified, .false.)
        ! Perform the (modified) Leja ordering.
        z = leja_ordering(x, modified_)
        ! Dispatch to the appropriate stability control driver.
        if (modified_) then
            y = modified_stability_control_csp(z)
        else
            y = standard_stability_control_csp(z)
        endif
    end procedure
    module procedure stabilized_leja_ordering_cdp
        logical(lk) :: modified_
        complex(dp), allocatable :: z(:)
        ! Standard or modified Leja ordering.
        modified_ = optval(modified, .false.)
        ! Perform the (modified) Leja ordering.
        z = leja_ordering(x, modified_)
        ! Dispatch to the appropriate stability control driver.
        if (modified_) then
            y = modified_stability_control_cdp(z)
        else
            y = standard_stability_control_cdp(z)
        endif
    end procedure

    pure function standard_stability_control_rsp(x) result(y)
        real(sp), intent(in)  :: x(:)
        !! Roots of the GMRES polynomial.
        real(sp), allocatable :: y(:)
        !! Roots of the stabilized polynomial.

        !----- Internal variables -----
        real(sp) :: pof(size(x))
        integer(ilp) :: i, k
        integer(ilp) :: num_added_roots(size(x))
        logical(lk)  :: mask(size(x))

        !> Product of other factors.
        do k = 1, size(x)
            mask = .true. ; mask(k) = .false.
            pof(k) = product(abs(1.0_sp - x(k)/x), mask)
        enddo

        !> Determine how many additional copies of each root needs to be added.
        do concurrent(k=1:size(x))
            num_added_roots(k) = ceiling( (log10(pof(k)) - 4.0_sp) / 14.0_sp )
        enddo

        !> Append the extra copies at the end.
        y = x
        do k = 1, size(x)
            if (num_added_roots(k) > 0) then
                do i = 1, num_added_roots(k)
                    y = [y, x(k)]
                enddo
            endif
        enddo
    end function
    pure function standard_stability_control_rdp(x) result(y)
        real(dp), intent(in)  :: x(:)
        !! Roots of the GMRES polynomial.
        real(dp), allocatable :: y(:)
        !! Roots of the stabilized polynomial.

        !----- Internal variables -----
        real(dp) :: pof(size(x))
        integer(ilp) :: i, k
        integer(ilp) :: num_added_roots(size(x))
        logical(lk)  :: mask(size(x))

        !> Product of other factors.
        do k = 1, size(x)
            mask = .true. ; mask(k) = .false.
            pof(k) = product(abs(1.0_dp - x(k)/x), mask)
        enddo

        !> Determine how many additional copies of each root needs to be added.
        do concurrent(k=1:size(x))
            num_added_roots(k) = ceiling( (log10(pof(k)) - 4.0_dp) / 14.0_dp )
        enddo

        !> Append the extra copies at the end.
        y = x
        do k = 1, size(x)
            if (num_added_roots(k) > 0) then
                do i = 1, num_added_roots(k)
                    y = [y, x(k)]
                enddo
            endif
        enddo
    end function
    pure function standard_stability_control_csp(x) result(y)
        complex(sp), intent(in)  :: x(:)
        !! Roots of the GMRES polynomial.
        complex(sp), allocatable :: y(:)
        !! Roots of the stabilized polynomial.

        !----- Internal variables -----
        real(sp) :: pof(size(x))
        integer(ilp) :: i, k
        integer(ilp) :: num_added_roots(size(x))
        logical(lk)  :: mask(size(x))

        !> Product of other factors.
        do k = 1, size(x)
            mask = .true. ; mask(k) = .false.
            pof(k) = product(abs(1.0_sp - x(k)/x), mask)
        enddo

        !> Determine how many additional copies of each root needs to be added.
        do concurrent(k=1:size(x))
            num_added_roots(k) = ceiling( (log10(pof(k)) - 4.0_sp) / 14.0_sp )
        enddo

        !> Append the extra copies at the end.
        y = x
        do k = 1, size(x)
            if (num_added_roots(k) > 0) then
                do i = 1, num_added_roots(k)
                    y = [y, x(k)]
                enddo
            endif
        enddo
    end function
    pure function standard_stability_control_cdp(x) result(y)
        complex(dp), intent(in)  :: x(:)
        !! Roots of the GMRES polynomial.
        complex(dp), allocatable :: y(:)
        !! Roots of the stabilized polynomial.

        !----- Internal variables -----
        real(dp) :: pof(size(x))
        integer(ilp) :: i, k
        integer(ilp) :: num_added_roots(size(x))
        logical(lk)  :: mask(size(x))

        !> Product of other factors.
        do k = 1, size(x)
            mask = .true. ; mask(k) = .false.
            pof(k) = product(abs(1.0_dp - x(k)/x), mask)
        enddo

        !> Determine how many additional copies of each root needs to be added.
        do concurrent(k=1:size(x))
            num_added_roots(k) = ceiling( (log10(pof(k)) - 4.0_dp) / 14.0_dp )
        enddo

        !> Append the extra copies at the end.
        y = x
        do k = 1, size(x)
            if (num_added_roots(k) > 0) then
                do i = 1, num_added_roots(k)
                    y = [y, x(k)]
                enddo
            endif
        enddo
    end function

    pure function modified_stability_control_csp(x) result(y)
        complex(sp), intent(in)  :: x(:)
        !! Roots of the GMRES polynomial.
        complex(sp), allocatable :: y(:)
        !! Roots of the stabilized polynomial.

        !----- Internal variables -----
        real(sp) :: pof(size(x)), re, im
        integer(ilp) :: i, k
        integer(ilp) :: num_added_roots(size(x))
        logical(lk)  :: mask(size(x))

        !> Product of other factors.
        do k = 1, size(x)
            mask = .true. ; mask(k) = .false.
            pof(k) = product(abs(1.0_sp - x(k)/x), mask)
        enddo

        !> Determine how many additional copies of each root needs to be added.
        do concurrent(k=1:size(x))
            num_added_roots(k) = ceiling( (log10(pof(k)) - 4.0_sp) / 14.0_sp )
        enddo

        !> Append the extra copies at the end.
        y = x ; k = 1
        do while (k < size(x))
            re = real(x(k)) ; im = imag(x(k))
            if (num_added_roots(k) > 0) then
                do i = 1, num_added_roots(k)
                    if (im == 0.0_sp) then
                        y = [y, x(k)]
                    else
                        y = [y, x(k), x(k+1)]
                    endif
                enddo
            endif
            k = merge(k+1, k+2, im == 0.0_sp)
        enddo
    end function
    pure function modified_stability_control_cdp(x) result(y)
        complex(dp), intent(in)  :: x(:)
        !! Roots of the GMRES polynomial.
        complex(dp), allocatable :: y(:)
        !! Roots of the stabilized polynomial.

        !----- Internal variables -----
        real(dp) :: pof(size(x)), re, im
        integer(ilp) :: i, k
        integer(ilp) :: num_added_roots(size(x))
        logical(lk)  :: mask(size(x))

        !> Product of other factors.
        do k = 1, size(x)
            mask = .true. ; mask(k) = .false.
            pof(k) = product(abs(1.0_dp - x(k)/x), mask)
        enddo

        !> Determine how many additional copies of each root needs to be added.
        do concurrent(k=1:size(x))
            num_added_roots(k) = ceiling( (log10(pof(k)) - 4.0_dp) / 14.0_dp )
        enddo

        !> Append the extra copies at the end.
        y = x ; k = 1
        do while (k < size(x))
            re = real(x(k)) ; im = imag(x(k))
            if (num_added_roots(k) > 0) then
                do i = 1, num_added_roots(k)
                    if (im == 0.0_dp) then
                        y = [y, x(k)]
                    else
                        y = [y, x(k), x(k+1)]
                    endif
                enddo
            endif
            k = merge(k+1, k+2, im == 0.0_dp)
        enddo
    end function
end submodule
