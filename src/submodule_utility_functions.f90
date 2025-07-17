submodule(lightkrylov_utils) utility_functions
    !-------------------------------
    !-----     LightKrylov     -----
    !-------------------------------
    use stdlib_optval, only: optval
    use stdlib_linalg_constants, only: ilp
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

    !----- Schur Factorization ------

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
        integer(ilp) :: i, n
        !> Problem's dimensions.
        n = size(A, 1) ; allocate(x, mold=b) ; x = 0.0_sp
        !> Back-substitution algorithm.
        x(n) = b(n) / A(n, n)
        do i = n-1, 1, -1
            x(i) = (b(i) - sum(A(i, i+1:) * x(i+1:))) / A(i, i)
        enddo
    end procedure
    module procedure solve_triangular_rdp
        integer(ilp) :: i, n
        !> Problem's dimensions.
        n = size(A, 1) ; allocate(x, mold=b) ; x = 0.0_dp
        !> Back-substitution algorithm.
        x(n) = b(n) / A(n, n)
        do i = n-1, 1, -1
            x(i) = (b(i) - sum(A(i, i+1:) * x(i+1:))) / A(i, i)
        enddo
    end procedure
    module procedure solve_triangular_csp
        integer(ilp) :: i, n
        !> Problem's dimensions.
        n = size(A, 1) ; allocate(x, mold=b) ; x = 0.0_sp
        !> Back-substitution algorithm.
        x(n) = b(n) / A(n, n)
        do i = n-1, 1, -1
            x(i) = (b(i) - sum(A(i, i+1:) * x(i+1:))) / A(i, i)
        enddo
    end procedure
    module procedure solve_triangular_cdp
        integer(ilp) :: i, n
        !> Problem's dimensions.
        n = size(A, 1) ; allocate(x, mold=b) ; x = 0.0_dp
        !> Back-substitution algorithm.
        x(n) = b(n) / A(n, n)
        do i = n-1, 1, -1
            x(i) = (b(i) - sum(A(i, i+1:) * x(i+1:))) / A(i, i)
        enddo
    end procedure
end submodule
