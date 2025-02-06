submodule (lightkrylov_basekrylov) krylov_utilities
    implicit none
contains

    !----------------------------------------
    !-----     Permutation matrices     -----
    !----------------------------------------

    module procedure permcols_basis_rsp
        call copy(Q, Q(perm))
    end procedure

    module procedure permcols_array_rsp
        Q = Q(:, perm)
    end procedure
    module procedure permcols_basis_rdp
        call copy(Q, Q(perm))
    end procedure

    module procedure permcols_array_rdp
        Q = Q(:, perm)
    end procedure
    module procedure permcols_basis_csp
        call copy(Q, Q(perm))
    end procedure

    module procedure permcols_array_csp
        Q = Q(:, perm)
    end procedure
    module procedure permcols_basis_cdp
        call copy(Q, Q(perm))
    end procedure

    module procedure permcols_array_cdp
        Q = Q(:, perm)
    end procedure

    module procedure invperm
        integer :: i
        inv_perm(perm) = [(i, i=1, size(perm))]
    end procedure

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

        return
    end procedure
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

        return
    end procedure
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

        return
    end procedure
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

        return
    end procedure
   
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

        return
    end procedure
    module procedure orthonormalize_basis_rdp
        character(len=*), parameter :: this_procedure = 'orthonormalize_basis_rdp'
        real(dp) :: R(size(X),size(X))
        integer :: info

        if (time_lightkrylov()) call timer%start(this_procedure)
        ! internals
        call qr(X, R, info)
        call check_info(info, 'qr', this_module, this_procedure)
        if (time_lightkrylov()) call timer%stop(this_procedure)

        return
    end procedure
    module procedure orthonormalize_basis_csp
        character(len=*), parameter :: this_procedure = 'orthonormalize_basis_csp'
        complex(sp) :: R(size(X),size(X))
        integer :: info

        if (time_lightkrylov()) call timer%start(this_procedure)
        ! internals
        call qr(X, R, info)
        call check_info(info, 'qr', this_module, this_procedure)
        if (time_lightkrylov()) call timer%stop(this_procedure)

        return
    end procedure
    module procedure orthonormalize_basis_cdp
        character(len=*), parameter :: this_procedure = 'orthonormalize_basis_cdp'
        complex(dp) :: R(size(X),size(X))
        integer :: info

        if (time_lightkrylov()) call timer%start(this_procedure)
        ! internals
        call qr(X, R, info)
        call check_info(info, 'qr', this_module, this_procedure)
        if (time_lightkrylov()) call timer%stop(this_procedure)

        return
    end procedure

    !---------------------------------------------------
    !-----     Check orthonormality of a basis     -----
    !---------------------------------------------------

    module procedure is_orthonormal_rsp
        real(sp), dimension(size(X), size(X)) :: G
        ortho = .true.
        G = Gram(X)
        if (mnorm(G - eye(size(X)), "Fro") > rtol_sp) then
            ! The basis is not orthonormal. Cannot orthonormalize.
            ortho = .false.
        end if
    end procedure
    module procedure is_orthonormal_rdp
        real(dp), dimension(size(X), size(X)) :: G
        ortho = .true.
        G = Gram(X)
        if (mnorm(G - eye(size(X)), "Fro") > rtol_sp) then
            ! The basis is not orthonormal. Cannot orthonormalize.
            ortho = .false.
        end if
    end procedure
    module procedure is_orthonormal_csp
        complex(sp), dimension(size(X), size(X)) :: G
        ortho = .true.
        G = Gram(X)
        if (mnorm(G - eye(size(X)), "Fro") > rtol_sp) then
            ! The basis is not orthonormal. Cannot orthonormalize.
            ortho = .false.
        end if
    end procedure
    module procedure is_orthonormal_cdp
        complex(dp), dimension(size(X), size(X)) :: G
        ortho = .true.
        G = Gram(X)
        if (mnorm(G - eye(size(X)), "Fro") > rtol_sp) then
            ! The basis is not orthonormal. Cannot orthonormalize.
            ortho = .false.
        end if
    end procedure

end submodule
