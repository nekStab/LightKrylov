submodule (lightkrylov_basekrylov) krylov_utilities
    implicit none(type, external)
contains

    !----------------------------------------
    !-----     Permutation matrices     -----
    !----------------------------------------

    module procedure permcols_basis_rsp
        class(abstract_vector_rsp), allocatable :: Qwrk(:)
        integer :: i
        allocate(Qwrk(size(Q)), mold=Q(1))
        call init_like_basis(Qwrk, Q(1))
        do i = 1, size(Q)
            call copy(Qwrk(i), Q(perm(i)))
        end do
        call copy(Q, Qwrk)
        call free_basis(Qwrk)
    end procedure permcols_basis_rsp

    module procedure permcols_array_rsp
        Q = Q(:, perm)
    end procedure permcols_array_rsp
    module procedure permcols_basis_rdp
        class(abstract_vector_rdp), allocatable :: Qwrk(:)
        integer :: i
        allocate(Qwrk(size(Q)), mold=Q(1))
        call init_like_basis(Qwrk, Q(1))
        do i = 1, size(Q)
            call copy(Qwrk(i), Q(perm(i)))
        end do
        call copy(Q, Qwrk)
        call free_basis(Qwrk)
    end procedure permcols_basis_rdp

    module procedure permcols_array_rdp
        Q = Q(:, perm)
    end procedure permcols_array_rdp
    module procedure permcols_basis_csp
        class(abstract_vector_csp), allocatable :: Qwrk(:)
        integer :: i
        allocate(Qwrk(size(Q)), mold=Q(1))
        call init_like_basis(Qwrk, Q(1))
        do i = 1, size(Q)
            call copy(Qwrk(i), Q(perm(i)))
        end do
        call copy(Q, Qwrk)
        call free_basis(Qwrk)
    end procedure permcols_basis_csp

    module procedure permcols_array_csp
        Q = Q(:, perm)
    end procedure permcols_array_csp
    module procedure permcols_basis_cdp
        class(abstract_vector_cdp), allocatable :: Qwrk(:)
        integer :: i
        allocate(Qwrk(size(Q)), mold=Q(1))
        call init_like_basis(Qwrk, Q(1))
        do i = 1, size(Q)
            call copy(Qwrk(i), Q(perm(i)))
        end do
        call copy(Q, Qwrk)
        call free_basis(Qwrk)
    end procedure permcols_basis_cdp

    module procedure permcols_array_cdp
        Q = Q(:, perm)
    end procedure permcols_array_cdp

    module procedure invperm
        integer :: i, iostat
        character(len=100) :: errmsg
        allocate(inv_perm(size(perm)), source=0, stat=iostat, errmsg=errmsg)
        call check_allocation(iostat, errmsg, this_module, "invperm")
        inv_perm(perm) = [(i, i=1, size(perm))]
    end procedure invperm

    !----------------------------------------------
    !-----     Initialize Krylov subspace     -----
    !----------------------------------------------

    module procedure initialize_krylov_subspace_rsp
        integer :: p

        ! Zero-out X.
        call zero_basis(X)

        ! Deals with optional args.
        if(present(X0)) then
            p = size(X0)
            ! Initialize.
            call copy(X(:p), X0)
            ! Orthonormalize.
            call orthonormalize_basis(X(:p))
        endif
    end procedure initialize_krylov_subspace_rsp
    module procedure initialize_krylov_subspace_rdp
        integer :: p

        ! Zero-out X.
        call zero_basis(X)

        ! Deals with optional args.
        if(present(X0)) then
            p = size(X0)
            ! Initialize.
            call copy(X(:p), X0)
            ! Orthonormalize.
            call orthonormalize_basis(X(:p))
        endif
    end procedure initialize_krylov_subspace_rdp
    module procedure initialize_krylov_subspace_csp
        integer :: p

        ! Zero-out X.
        call zero_basis(X)

        ! Deals with optional args.
        if(present(X0)) then
            p = size(X0)
            ! Initialize.
            call copy(X(:p), X0)
            ! Orthonormalize.
            call orthonormalize_basis(X(:p))
        endif
    end procedure initialize_krylov_subspace_csp
    module procedure initialize_krylov_subspace_cdp
        integer :: p

        ! Zero-out X.
        call zero_basis(X)

        ! Deals with optional args.
        if(present(X0)) then
            p = size(X0)
            ! Initialize.
            call copy(X(:p), X0)
            ! Orthonormalize.
            call orthonormalize_basis(X(:p))
        endif
    end procedure initialize_krylov_subspace_cdp

    !-------------------------------------------------------
    !-----     Initialize random orthonormal basis     -----
    !-------------------------------------------------------

    module procedure initialize_random_orthonormal_basis_rsp
        integer :: p

        ! Fill with random data
        call rand_basis(X)

        ! Orthonormalize
        call orthonormalize_basis(X)
    end procedure initialize_random_orthonormal_basis_rsp
    module procedure initialize_random_orthonormal_basis_rdp
        integer :: p

        ! Fill with random data
        call rand_basis(X)

        ! Orthonormalize
        call orthonormalize_basis(X)
    end procedure initialize_random_orthonormal_basis_rdp
    module procedure initialize_random_orthonormal_basis_csp
        integer :: p

        ! Fill with random data
        call rand_basis(X)

        ! Orthonormalize
        call orthonormalize_basis(X)
    end procedure initialize_random_orthonormal_basis_csp
    module procedure initialize_random_orthonormal_basis_cdp
        integer :: p

        ! Fill with random data
        call rand_basis(X)

        ! Orthonormalize
        call orthonormalize_basis(X)
    end procedure initialize_random_orthonormal_basis_cdp

    !----------------------------------------
    !-----     Orthonormalize basis     -----
    !----------------------------------------

    module procedure orthonormalize_basis_rsp
        character(len=*), parameter :: this_procedure = 'orthonormalize_basis_rsp'
        real(sp) :: R(size(X),size(X))
        integer :: info

        if (time_lightkrylov()) call timer%start(this_procedure)
        ! internals
        call qr(X, R, info)
        call check_info(info, 'qr', this_module, this_procedure)
        if (time_lightkrylov()) call timer%stop(this_procedure)
    end procedure orthonormalize_basis_rsp

    module procedure orthonormalize_basis_rdp
        character(len=*), parameter :: this_procedure = 'orthonormalize_basis_rdp'
        real(dp) :: R(size(X),size(X))
        integer :: info

        if (time_lightkrylov()) call timer%start(this_procedure)
        ! internals
        call qr(X, R, info)
        call check_info(info, 'qr', this_module, this_procedure)
        if (time_lightkrylov()) call timer%stop(this_procedure)
    end procedure orthonormalize_basis_rdp

    module procedure orthonormalize_basis_csp
        character(len=*), parameter :: this_procedure = 'orthonormalize_basis_csp'
        complex(sp) :: R(size(X),size(X))
        integer :: info

        if (time_lightkrylov()) call timer%start(this_procedure)
        ! internals
        call qr(X, R, info)
        call check_info(info, 'qr', this_module, this_procedure)
        if (time_lightkrylov()) call timer%stop(this_procedure)
    end procedure orthonormalize_basis_csp

    module procedure orthonormalize_basis_cdp
        character(len=*), parameter :: this_procedure = 'orthonormalize_basis_cdp'
        complex(dp) :: R(size(X),size(X))
        integer :: info

        if (time_lightkrylov()) call timer%start(this_procedure)
        ! internals
        call qr(X, R, info)
        call check_info(info, 'qr', this_module, this_procedure)
        if (time_lightkrylov()) call timer%stop(this_procedure)
    end procedure orthonormalize_basis_cdp

    !---------------------------------------------------
    !-----     Check orthonormality of a basis     -----
    !---------------------------------------------------

    module procedure is_orthonormal_rsp
        real(sp), dimension(size(X), size(X)) :: G
        ortho = .true.
        G = Gram(X)
        if (mnorm(G - eye(size(X)), "Fro") > rtol_sp) then
            ! The basis is not orthonormal.
            ortho = .false.
        end if
    end procedure is_orthonormal_rsp

    module procedure is_orthonormal_rdp
        real(dp), dimension(size(X), size(X)) :: G
        ortho = .true.
        G = Gram(X)
        if (mnorm(G - eye(size(X)), "Fro") > rtol_sp) then
            ! The basis is not orthonormal.
            ortho = .false.
        end if
    end procedure is_orthonormal_rdp

    module procedure is_orthonormal_csp
        complex(sp), dimension(size(X), size(X)) :: G
        ortho = .true.
        G = Gram(X)
        if (mnorm(G - eye(size(X)), "Fro") > rtol_sp) then
            ! The basis is not orthonormal.
            ortho = .false.
        end if
    end procedure is_orthonormal_csp

    module procedure is_orthonormal_cdp
        complex(dp), dimension(size(X), size(X)) :: G
        ortho = .true.
        G = Gram(X)
        if (mnorm(G - eye(size(X)), "Fro") > rtol_sp) then
            ! The basis is not orthonormal.
            ortho = .false.
        end if
    end procedure is_orthonormal_cdp

    !----------------------------------------------
    !-----     Biorthonormalize two bases     -----
    !----------------------------------------------

    module procedure biorthonormalize_bases_rsp
        character(len=*), parameter :: this_procedure = 'biorthonormalize_bases_rsp'
        !! SVD workspace
        integer  :: n, info_, nretain
        real(sp), allocatable  :: M(:,:), U(:,:), VT(:,:)
        real(sp), allocatable :: S(:)
        real(sp) :: tol_

        if (time_lightkrylov()) call timer%start(this_procedure)

        ! Check sizes.
        if (size(X) /= size(Y)) then
            call stop_error("Krylov bases X and Y have different sizes.", &
                              & this_module, this_procedure)
        endif

        ! handle optional tol
        tol_ = optval(tol, atol_sp)
        n = size(X)
        
        ! compute SVD of inner product matrix
        M = innerprod(Y, X)
        allocate(S(n), U(n, n), VT(n, n))
        call svd(M, S, U, VT)

        ! count + renormalize retained singular values; zero the rest
        nretain = count(S / S(1) > tol_)
        call check_info(merge(-1, 0, nretain == 0), 'biorthonormalize_bases', &
            & this_module, this_procedure)
        S(:nretain) = one_rsp / sqrt(S(:nretain))
        S(nretain+1:) = zero_rsp

        ! apply symmetric transformation in-place
        block
            class(abstract_vector_rsp), allocatable :: Xwrk(:)
            call linear_combination(Xwrk, X, matmul(transpose(VT(:nretain, :)), diag(S(:nretain))))
            call copy(X(:nretain), Xwrk)
            call zero_basis(X(nretain+1:))
            call linear_combination(Xwrk, Y, matmul(U(:, :nretain), diag(S(:nretain))))
            call copy(Y(:nretain), Xwrk)
            call zero_basis(Y(nretain+1:))
        end block

        if (present(info)) info = nretain
        if (time_lightkrylov()) call timer%stop(this_procedure)

    end procedure biorthonormalize_bases_rsp

    module procedure biorthonormalize_bases_rdp
        character(len=*), parameter :: this_procedure = 'biorthonormalize_bases_rdp'
        !! SVD workspace
        integer  :: n, info_, nretain
        real(dp), allocatable  :: M(:,:), U(:,:), VT(:,:)
        real(dp), allocatable :: S(:)
        real(dp) :: tol_

        if (time_lightkrylov()) call timer%start(this_procedure)

        ! Check sizes.
        if (size(X) /= size(Y)) then
            call stop_error("Krylov bases X and Y have different sizes.", &
                              & this_module, this_procedure)
        endif

        ! handle optional tol
        tol_ = optval(tol, atol_dp)
        n = size(X)
        
        ! compute SVD of inner product matrix
        M = innerprod(Y, X)
        allocate(S(n), U(n, n), VT(n, n))
        call svd(M, S, U, VT)

        ! count + renormalize retained singular values; zero the rest
        nretain = count(S / S(1) > tol_)
        call check_info(merge(-1, 0, nretain == 0), 'biorthonormalize_bases', &
            & this_module, this_procedure)
        S(:nretain) = one_rdp / sqrt(S(:nretain))
        S(nretain+1:) = zero_rdp

        ! apply symmetric transformation in-place
        block
            class(abstract_vector_rdp), allocatable :: Xwrk(:)
            call linear_combination(Xwrk, X, matmul(transpose(VT(:nretain, :)), diag(S(:nretain))))
            call copy(X(:nretain), Xwrk)
            call zero_basis(X(nretain+1:))
            call linear_combination(Xwrk, Y, matmul(U(:, :nretain), diag(S(:nretain))))
            call copy(Y(:nretain), Xwrk)
            call zero_basis(Y(nretain+1:))
        end block

        if (present(info)) info = nretain
        if (time_lightkrylov()) call timer%stop(this_procedure)

    end procedure biorthonormalize_bases_rdp

    module procedure biorthonormalize_bases_csp
        character(len=*), parameter :: this_procedure = 'biorthonormalize_bases_csp'
        !! SVD workspace
        integer  :: n, info_, nretain
        complex(sp), allocatable  :: M(:,:), U(:,:), VT(:,:)
        real(sp), allocatable :: S(:)
        real(sp) :: tol_

        if (time_lightkrylov()) call timer%start(this_procedure)

        ! Check sizes.
        if (size(X) /= size(Y)) then
            call stop_error("Krylov bases X and Y have different sizes.", &
                              & this_module, this_procedure)
        endif

        ! handle optional tol
        tol_ = optval(tol, atol_sp)
        n = size(X)
        
        ! compute SVD of inner product matrix
        M = innerprod(Y, X)
        allocate(S(n), U(n, n), VT(n, n))
        call svd(M, S, U, VT)

        ! count + renormalize retained singular values; zero the rest
        nretain = count(S / S(1) > tol_)
        call check_info(merge(-1, 0, nretain == 0), 'biorthonormalize_bases', &
            & this_module, this_procedure)
        S(:nretain) = one_rsp / sqrt(S(:nretain))
        S(nretain+1:) = zero_rsp

        ! apply symmetric transformation in-place
        block
            class(abstract_vector_csp), allocatable :: Xwrk(:)
            call linear_combination(Xwrk, X, matmul(hermitian(VT(:nretain, :)), diag(S(:nretain))))
            call copy(X(:nretain), Xwrk)
            call zero_basis(X(nretain+1:))
            call linear_combination(Xwrk, Y, matmul(U(:, :nretain), diag(S(:nretain))))
            call copy(Y(:nretain), Xwrk)
            call zero_basis(Y(nretain+1:))
        end block

        if (present(info)) info = nretain
        if (time_lightkrylov()) call timer%stop(this_procedure)

    end procedure biorthonormalize_bases_csp

    module procedure biorthonormalize_bases_cdp
        character(len=*), parameter :: this_procedure = 'biorthonormalize_bases_cdp'
        !! SVD workspace
        integer  :: n, info_, nretain
        complex(dp), allocatable  :: M(:,:), U(:,:), VT(:,:)
        real(dp), allocatable :: S(:)
        real(dp) :: tol_

        if (time_lightkrylov()) call timer%start(this_procedure)

        ! Check sizes.
        if (size(X) /= size(Y)) then
            call stop_error("Krylov bases X and Y have different sizes.", &
                              & this_module, this_procedure)
        endif

        ! handle optional tol
        tol_ = optval(tol, atol_dp)
        n = size(X)
        
        ! compute SVD of inner product matrix
        M = innerprod(Y, X)
        allocate(S(n), U(n, n), VT(n, n))
        call svd(M, S, U, VT)

        ! count + renormalize retained singular values; zero the rest
        nretain = count(S / S(1) > tol_)
        call check_info(merge(-1, 0, nretain == 0), 'biorthonormalize_bases', &
            & this_module, this_procedure)
        S(:nretain) = one_rdp / sqrt(S(:nretain))
        S(nretain+1:) = zero_rdp

        ! apply symmetric transformation in-place
        block
            class(abstract_vector_cdp), allocatable :: Xwrk(:)
            call linear_combination(Xwrk, X, matmul(hermitian(VT(:nretain, :)), diag(S(:nretain))))
            call copy(X(:nretain), Xwrk)
            call zero_basis(X(nretain+1:))
            call linear_combination(Xwrk, Y, matmul(U(:, :nretain), diag(S(:nretain))))
            call copy(Y(:nretain), Xwrk)
            call zero_basis(Y(nretain+1:))
        end block

        if (present(info)) info = nretain
        if (time_lightkrylov()) call timer%stop(this_procedure)

    end procedure biorthonormalize_bases_cdp

end submodule krylov_utilities
