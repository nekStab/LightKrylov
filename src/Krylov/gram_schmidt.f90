submodule (lightkrylov_basekrylov) gram_schmidt_process
    implicit none(type, external)
contains

    !--------------------------------------------
    !-----     DOUBLE GRAM-SCHMIDT STEP     -----
    !--------------------------------------------

    module procedure DGS_vector_against_basis_rsp
        character(len=*), parameter :: this_procedure = 'DGS_vector_against_basis_rsp'
        real(sp), dimension(size(X)) :: proj_coefficients, wrk
        ! internal
        logical :: chk_X_orthonormality
        logical :: is_ortho

        ! optional input argument
        chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

        if (time_lightkrylov()) call timer%start(this_procedure)
        info = 0

        ! optional orthonormality check
        if (chk_X_orthonormality) then
            is_ortho = is_orthonormal_rsp(X)
            if (is_ortho) then
                call log_information("Input basis orthonormal. Remove this check unless necessary for better performance", &
                & this_module, this_procedure)
            else
                call stop_error("Input basis not orthonormal.", this_module, this_procedure)
            end if
        end if

        proj_coefficients = zero_rsp; wrk = zero_rsp

        ! Orthogonalize vector y w.r.t. to Krylov basis X in two passes of GS.
        ! first pass
        call orthogonalize_against_basis(y, X, info, &
                                         if_chk_orthonormal=.false., beta=proj_coefficients)
        call check_info(info, 'orthogonalize_against_basis_p1', &
                        this_module, this_procedure//', pass 1')
        ! second pass
        call orthogonalize_against_basis(y, X, info, if_chk_orthonormal=.false., beta=wrk)
        call check_info(info, 'orthogonalize_against_basis_p2', &
                        this_module, this_procedure//', pass 2')
        ! combine passes
        proj_coefficients = proj_coefficients + wrk

        if (present(beta)) then
            ! check size
            call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, this_procedure)
            beta = proj_coefficients
        end if
        if (time_lightkrylov()) call timer%stop(this_procedure)
    end procedure DGS_vector_against_basis_rsp

    module procedure DGS_basis_against_basis_rsp
        character(len=*), parameter :: this_procedure = 'DGS_basis_against_basis_rsp'
        real(sp), dimension(size(X),size(Y)) :: proj_coefficients, wrk
        ! internal
        logical :: chk_X_orthonormality
        logical :: is_ortho

        ! optional input argument
        chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

        if (time_lightkrylov()) call timer%start(this_procedure)
        info = 0

        ! optional orthonormality check
        if (chk_X_orthonormality) then
            is_ortho = is_orthonormal_rsp(X)
            if (is_ortho) then
                call log_information("Input basis orthonormal. Remove this check unless necessary for better performance", &
                & this_module, this_procedure)
            else
                call stop_error("Input basis not orthonormal.", this_module, this_procedure)
            end if
        end if

        proj_coefficients = zero_rsp; wrk = zero_rsp

        ! Orthogonalize Krylov basis Y w.r.t. to Krylov basis X in two passes of GS.
        ! first pass
        call orthogonalize_against_basis(Y, X, info, &
                                         if_chk_orthonormal=.false., beta=proj_coefficients)
        call check_info(info, 'orthogonalize_against_basis_p1', &
                        this_module, this_procedure//', first pass')
        ! second pass
        call orthogonalize_against_basis(Y, X, info, &
                                         if_chk_orthonormal=.false., beta=wrk)
        call check_info(info, 'orthogonalize_against_basis_p2', &
                        this_module, this_procedure//', second pass')
        ! combine passes
        proj_coefficients = proj_coefficients + wrk

        if (present(beta)) then
            ! check size
            call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, this_procedure)
            beta = proj_coefficients
        end if
        if (time_lightkrylov()) call timer%stop(this_procedure)
    end procedure DGS_basis_against_basis_rsp

    module procedure DGS_vector_against_basis_rdp
        character(len=*), parameter :: this_procedure = 'DGS_vector_against_basis_rdp'
        real(dp), dimension(size(X)) :: proj_coefficients, wrk
        ! internal
        logical :: chk_X_orthonormality
        logical :: is_ortho

        ! optional input argument
        chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

        if (time_lightkrylov()) call timer%start(this_procedure)
        info = 0

        ! optional orthonormality check
        if (chk_X_orthonormality) then
            is_ortho = is_orthonormal_rdp(X)
            if (is_ortho) then
                call log_information("Input basis orthonormal. Remove this check unless necessary for better performance", &
                & this_module, this_procedure)
            else
                call stop_error("Input basis not orthonormal.", this_module, this_procedure)
            end if
        end if

        proj_coefficients = zero_rdp; wrk = zero_rdp

        ! Orthogonalize vector y w.r.t. to Krylov basis X in two passes of GS.
        ! first pass
        call orthogonalize_against_basis(y, X, info, &
                                         if_chk_orthonormal=.false., beta=proj_coefficients)
        call check_info(info, 'orthogonalize_against_basis_p1', &
                        this_module, this_procedure//', pass 1')
        ! second pass
        call orthogonalize_against_basis(y, X, info, if_chk_orthonormal=.false., beta=wrk)
        call check_info(info, 'orthogonalize_against_basis_p2', &
                        this_module, this_procedure//', pass 2')
        ! combine passes
        proj_coefficients = proj_coefficients + wrk

        if (present(beta)) then
            ! check size
            call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, this_procedure)
            beta = proj_coefficients
        end if
        if (time_lightkrylov()) call timer%stop(this_procedure)
    end procedure DGS_vector_against_basis_rdp

    module procedure DGS_basis_against_basis_rdp
        character(len=*), parameter :: this_procedure = 'DGS_basis_against_basis_rdp'
        real(dp), dimension(size(X),size(Y)) :: proj_coefficients, wrk
        ! internal
        logical :: chk_X_orthonormality
        logical :: is_ortho

        ! optional input argument
        chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

        if (time_lightkrylov()) call timer%start(this_procedure)
        info = 0

        ! optional orthonormality check
        if (chk_X_orthonormality) then
            is_ortho = is_orthonormal_rdp(X)
            if (is_ortho) then
                call log_information("Input basis orthonormal. Remove this check unless necessary for better performance", &
                & this_module, this_procedure)
            else
                call stop_error("Input basis not orthonormal.", this_module, this_procedure)
            end if
        end if

        proj_coefficients = zero_rdp; wrk = zero_rdp

        ! Orthogonalize Krylov basis Y w.r.t. to Krylov basis X in two passes of GS.
        ! first pass
        call orthogonalize_against_basis(Y, X, info, &
                                         if_chk_orthonormal=.false., beta=proj_coefficients)
        call check_info(info, 'orthogonalize_against_basis_p1', &
                        this_module, this_procedure//', first pass')
        ! second pass
        call orthogonalize_against_basis(Y, X, info, &
                                         if_chk_orthonormal=.false., beta=wrk)
        call check_info(info, 'orthogonalize_against_basis_p2', &
                        this_module, this_procedure//', second pass')
        ! combine passes
        proj_coefficients = proj_coefficients + wrk

        if (present(beta)) then
            ! check size
            call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, this_procedure)
            beta = proj_coefficients
        end if
        if (time_lightkrylov()) call timer%stop(this_procedure)
    end procedure DGS_basis_against_basis_rdp

    module procedure DGS_vector_against_basis_csp
        character(len=*), parameter :: this_procedure = 'DGS_vector_against_basis_csp'
        complex(sp), dimension(size(X)) :: proj_coefficients, wrk
        ! internal
        logical :: chk_X_orthonormality
        logical :: is_ortho

        ! optional input argument
        chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

        if (time_lightkrylov()) call timer%start(this_procedure)
        info = 0

        ! optional orthonormality check
        if (chk_X_orthonormality) then
            is_ortho = is_orthonormal_csp(X)
            if (is_ortho) then
                call log_information("Input basis orthonormal. Remove this check unless necessary for better performance", &
                & this_module, this_procedure)
            else
                call stop_error("Input basis not orthonormal.", this_module, this_procedure)
            end if
        end if

        proj_coefficients = zero_csp; wrk = zero_csp

        ! Orthogonalize vector y w.r.t. to Krylov basis X in two passes of GS.
        ! first pass
        call orthogonalize_against_basis(y, X, info, &
                                         if_chk_orthonormal=.false., beta=proj_coefficients)
        call check_info(info, 'orthogonalize_against_basis_p1', &
                        this_module, this_procedure//', pass 1')
        ! second pass
        call orthogonalize_against_basis(y, X, info, if_chk_orthonormal=.false., beta=wrk)
        call check_info(info, 'orthogonalize_against_basis_p2', &
                        this_module, this_procedure//', pass 2')
        ! combine passes
        proj_coefficients = proj_coefficients + wrk

        if (present(beta)) then
            ! check size
            call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, this_procedure)
            beta = proj_coefficients
        end if
        if (time_lightkrylov()) call timer%stop(this_procedure)
    end procedure DGS_vector_against_basis_csp

    module procedure DGS_basis_against_basis_csp
        character(len=*), parameter :: this_procedure = 'DGS_basis_against_basis_csp'
        complex(sp), dimension(size(X),size(Y)) :: proj_coefficients, wrk
        ! internal
        logical :: chk_X_orthonormality
        logical :: is_ortho

        ! optional input argument
        chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

        if (time_lightkrylov()) call timer%start(this_procedure)
        info = 0

        ! optional orthonormality check
        if (chk_X_orthonormality) then
            is_ortho = is_orthonormal_csp(X)
            if (is_ortho) then
                call log_information("Input basis orthonormal. Remove this check unless necessary for better performance", &
                & this_module, this_procedure)
            else
                call stop_error("Input basis not orthonormal.", this_module, this_procedure)
            end if
        end if

        proj_coefficients = zero_csp; wrk = zero_csp

        ! Orthogonalize Krylov basis Y w.r.t. to Krylov basis X in two passes of GS.
        ! first pass
        call orthogonalize_against_basis(Y, X, info, &
                                         if_chk_orthonormal=.false., beta=proj_coefficients)
        call check_info(info, 'orthogonalize_against_basis_p1', &
                        this_module, this_procedure//', first pass')
        ! second pass
        call orthogonalize_against_basis(Y, X, info, &
                                         if_chk_orthonormal=.false., beta=wrk)
        call check_info(info, 'orthogonalize_against_basis_p2', &
                        this_module, this_procedure//', second pass')
        ! combine passes
        proj_coefficients = proj_coefficients + wrk

        if (present(beta)) then
            ! check size
            call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, this_procedure)
            beta = proj_coefficients
        end if
        if (time_lightkrylov()) call timer%stop(this_procedure)
    end procedure DGS_basis_against_basis_csp

    module procedure DGS_vector_against_basis_cdp
        character(len=*), parameter :: this_procedure = 'DGS_vector_against_basis_cdp'
        complex(dp), dimension(size(X)) :: proj_coefficients, wrk
        ! internal
        logical :: chk_X_orthonormality
        logical :: is_ortho

        ! optional input argument
        chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

        if (time_lightkrylov()) call timer%start(this_procedure)
        info = 0

        ! optional orthonormality check
        if (chk_X_orthonormality) then
            is_ortho = is_orthonormal_cdp(X)
            if (is_ortho) then
                call log_information("Input basis orthonormal. Remove this check unless necessary for better performance", &
                & this_module, this_procedure)
            else
                call stop_error("Input basis not orthonormal.", this_module, this_procedure)
            end if
        end if

        proj_coefficients = zero_cdp; wrk = zero_cdp

        ! Orthogonalize vector y w.r.t. to Krylov basis X in two passes of GS.
        ! first pass
        call orthogonalize_against_basis(y, X, info, &
                                         if_chk_orthonormal=.false., beta=proj_coefficients)
        call check_info(info, 'orthogonalize_against_basis_p1', &
                        this_module, this_procedure//', pass 1')
        ! second pass
        call orthogonalize_against_basis(y, X, info, if_chk_orthonormal=.false., beta=wrk)
        call check_info(info, 'orthogonalize_against_basis_p2', &
                        this_module, this_procedure//', pass 2')
        ! combine passes
        proj_coefficients = proj_coefficients + wrk

        if (present(beta)) then
            ! check size
            call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, this_procedure)
            beta = proj_coefficients
        end if
        if (time_lightkrylov()) call timer%stop(this_procedure)
    end procedure DGS_vector_against_basis_cdp

    module procedure DGS_basis_against_basis_cdp
        character(len=*), parameter :: this_procedure = 'DGS_basis_against_basis_cdp'
        complex(dp), dimension(size(X),size(Y)) :: proj_coefficients, wrk
        ! internal
        logical :: chk_X_orthonormality
        logical :: is_ortho

        ! optional input argument
        chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

        if (time_lightkrylov()) call timer%start(this_procedure)
        info = 0

        ! optional orthonormality check
        if (chk_X_orthonormality) then
            is_ortho = is_orthonormal_cdp(X)
            if (is_ortho) then
                call log_information("Input basis orthonormal. Remove this check unless necessary for better performance", &
                & this_module, this_procedure)
            else
                call stop_error("Input basis not orthonormal.", this_module, this_procedure)
            end if
        end if

        proj_coefficients = zero_cdp; wrk = zero_cdp

        ! Orthogonalize Krylov basis Y w.r.t. to Krylov basis X in two passes of GS.
        ! first pass
        call orthogonalize_against_basis(Y, X, info, &
                                         if_chk_orthonormal=.false., beta=proj_coefficients)
        call check_info(info, 'orthogonalize_against_basis_p1', &
                        this_module, this_procedure//', first pass')
        ! second pass
        call orthogonalize_against_basis(Y, X, info, &
                                         if_chk_orthonormal=.false., beta=wrk)
        call check_info(info, 'orthogonalize_against_basis_p2', &
                        this_module, this_procedure//', second pass')
        ! combine passes
        proj_coefficients = proj_coefficients + wrk

        if (present(beta)) then
            ! check size
            call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, this_procedure)
            beta = proj_coefficients
        end if
        if (time_lightkrylov()) call timer%stop(this_procedure)
    end procedure DGS_basis_against_basis_cdp

    !-------------------------------------------------------------------------------
    !-----     ORTHOGONALIZE VECTOR/BASIS AGAINST ALREADY ORTHOGONAL BASIS     -----
    !-------------------------------------------------------------------------------

    module procedure orthogonalize_vector_against_basis_rsp
        character(len=*), parameter :: this_procedure = 'orthogonalize_vector_against_basis_rsp'
        real(sp) :: proj_coefficients(size(X))
        ! internal
        logical :: chk_X_orthonormality
        logical :: is_ortho

        if (time_lightkrylov()) call timer%start(this_procedure)
        info = 0

        ! optional input argument
        chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

        ! check for zero vector
        if (y%norm() < atol_sp) info = 1

        ! optional orthonormality check
        if (chk_X_orthonormality) then
            is_ortho = is_orthonormal_rsp(X)
            if (is_ortho) then
                call log_information("Input basis orthonormal. Remove this check unless necessary for better performance", &
                & this_module, this_procedure)
            else
                call stop_error("Input basis not orthonormal.", this_module, this_procedure)
            end if
        end if

        ! orthogonalize
        proj_coefficients = innerprod(X, y)
        ! y <- y - X c, subtracted term by term: no temporary vector is needed.
        block
            integer :: i
            do i = 1, size(X)
                call y%axpby(-proj_coefficients(i), X(i), one_rsp)
            enddo
        end block

        if (present(beta)) then
            ! check size
            call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, this_procedure)
            beta = proj_coefficients
        end if
        if (time_lightkrylov()) call timer%stop(this_procedure)
    end procedure orthogonalize_vector_against_basis_rsp

    module procedure orthogonalize_basis_against_basis_rsp
        character(len=*), parameter :: this_procedure = 'orthogonalize_basis_against_basis_rsp'
        real(sp) :: proj_coefficients(size(X), size(Y))
        ! internal
        integer :: i
        logical :: chk_X_orthonormality
        logical :: is_ortho

        if (time_lightkrylov()) call timer%start(this_procedure)
        info = 0

        ! optional input argument
        chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

        ! check for zero vector
        do i = 1, size(Y)
            if (Y(i)%norm() < atol_sp) info = i
        end do

        ! optional orthonormality check
        if (chk_X_orthonormality) then
            is_ortho = is_orthonormal_rsp(X)
            if (is_ortho) then
                call log_information("Input basis orthonormal. Remove this check unless necessary for better performance", &
                & this_module, this_procedure)
            else
                call stop_error("Input basis not orthonormal.", this_module, this_procedure)
            end if
        end if

        ! orthogonalize
        proj_coefficients = innerprod(X, Y)
        ! Y <- Y - X C, subtracted term by term: no temporary basis is needed.
        block
            integer :: ii, jj
            do jj = 1, size(Y)
                do ii = 1, size(X)
                    call Y(jj)%axpby(-proj_coefficients(ii, jj), X(ii), one_rsp)
                enddo
            enddo
        end block

        if (present(beta)) then
            ! check size
            call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, this_procedure)
            beta = proj_coefficients
        end if
        if (time_lightkrylov()) call timer%stop(this_procedure)
    end procedure orthogonalize_basis_against_basis_rsp

    module procedure orthogonalize_vector_against_basis_rdp
        character(len=*), parameter :: this_procedure = 'orthogonalize_vector_against_basis_rdp'
        real(dp) :: proj_coefficients(size(X))
        ! internal
        logical :: chk_X_orthonormality
        logical :: is_ortho

        if (time_lightkrylov()) call timer%start(this_procedure)
        info = 0

        ! optional input argument
        chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

        ! check for zero vector
        if (y%norm() < atol_dp) info = 1

        ! optional orthonormality check
        if (chk_X_orthonormality) then
            is_ortho = is_orthonormal_rdp(X)
            if (is_ortho) then
                call log_information("Input basis orthonormal. Remove this check unless necessary for better performance", &
                & this_module, this_procedure)
            else
                call stop_error("Input basis not orthonormal.", this_module, this_procedure)
            end if
        end if

        ! orthogonalize
        proj_coefficients = innerprod(X, y)
        ! y <- y - X c, subtracted term by term: no temporary vector is needed.
        block
            integer :: i
            do i = 1, size(X)
                call y%axpby(-proj_coefficients(i), X(i), one_rdp)
            enddo
        end block

        if (present(beta)) then
            ! check size
            call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, this_procedure)
            beta = proj_coefficients
        end if
        if (time_lightkrylov()) call timer%stop(this_procedure)
    end procedure orthogonalize_vector_against_basis_rdp

    module procedure orthogonalize_basis_against_basis_rdp
        character(len=*), parameter :: this_procedure = 'orthogonalize_basis_against_basis_rdp'
        real(dp) :: proj_coefficients(size(X), size(Y))
        ! internal
        integer :: i
        logical :: chk_X_orthonormality
        logical :: is_ortho

        if (time_lightkrylov()) call timer%start(this_procedure)
        info = 0

        ! optional input argument
        chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

        ! check for zero vector
        do i = 1, size(Y)
            if (Y(i)%norm() < atol_dp) info = i
        end do

        ! optional orthonormality check
        if (chk_X_orthonormality) then
            is_ortho = is_orthonormal_rdp(X)
            if (is_ortho) then
                call log_information("Input basis orthonormal. Remove this check unless necessary for better performance", &
                & this_module, this_procedure)
            else
                call stop_error("Input basis not orthonormal.", this_module, this_procedure)
            end if
        end if

        ! orthogonalize
        proj_coefficients = innerprod(X, Y)
        ! Y <- Y - X C, subtracted term by term: no temporary basis is needed.
        block
            integer :: ii, jj
            do jj = 1, size(Y)
                do ii = 1, size(X)
                    call Y(jj)%axpby(-proj_coefficients(ii, jj), X(ii), one_rdp)
                enddo
            enddo
        end block

        if (present(beta)) then
            ! check size
            call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, this_procedure)
            beta = proj_coefficients
        end if
        if (time_lightkrylov()) call timer%stop(this_procedure)
    end procedure orthogonalize_basis_against_basis_rdp

    module procedure orthogonalize_vector_against_basis_csp
        character(len=*), parameter :: this_procedure = 'orthogonalize_vector_against_basis_csp'
        complex(sp) :: proj_coefficients(size(X))
        ! internal
        logical :: chk_X_orthonormality
        logical :: is_ortho

        if (time_lightkrylov()) call timer%start(this_procedure)
        info = 0

        ! optional input argument
        chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

        ! check for zero vector
        if (y%norm() < atol_sp) info = 1

        ! optional orthonormality check
        if (chk_X_orthonormality) then
            is_ortho = is_orthonormal_csp(X)
            if (is_ortho) then
                call log_information("Input basis orthonormal. Remove this check unless necessary for better performance", &
                & this_module, this_procedure)
            else
                call stop_error("Input basis not orthonormal.", this_module, this_procedure)
            end if
        end if

        ! orthogonalize
        proj_coefficients = innerprod(X, y)
        ! y <- y - X c, subtracted term by term: no temporary vector is needed.
        block
            integer :: i
            do i = 1, size(X)
                call y%axpby(-proj_coefficients(i), X(i), one_csp)
            enddo
        end block

        if (present(beta)) then
            ! check size
            call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, this_procedure)
            beta = proj_coefficients
        end if
        if (time_lightkrylov()) call timer%stop(this_procedure)
    end procedure orthogonalize_vector_against_basis_csp

    module procedure orthogonalize_basis_against_basis_csp
        character(len=*), parameter :: this_procedure = 'orthogonalize_basis_against_basis_csp'
        complex(sp) :: proj_coefficients(size(X), size(Y))
        ! internal
        integer :: i
        logical :: chk_X_orthonormality
        logical :: is_ortho

        if (time_lightkrylov()) call timer%start(this_procedure)
        info = 0

        ! optional input argument
        chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

        ! check for zero vector
        do i = 1, size(Y)
            if (Y(i)%norm() < atol_sp) info = i
        end do

        ! optional orthonormality check
        if (chk_X_orthonormality) then
            is_ortho = is_orthonormal_csp(X)
            if (is_ortho) then
                call log_information("Input basis orthonormal. Remove this check unless necessary for better performance", &
                & this_module, this_procedure)
            else
                call stop_error("Input basis not orthonormal.", this_module, this_procedure)
            end if
        end if

        ! orthogonalize
        proj_coefficients = innerprod(X, Y)
        ! Y <- Y - X C, subtracted term by term: no temporary basis is needed.
        block
            integer :: ii, jj
            do jj = 1, size(Y)
                do ii = 1, size(X)
                    call Y(jj)%axpby(-proj_coefficients(ii, jj), X(ii), one_csp)
                enddo
            enddo
        end block

        if (present(beta)) then
            ! check size
            call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, this_procedure)
            beta = proj_coefficients
        end if
        if (time_lightkrylov()) call timer%stop(this_procedure)
    end procedure orthogonalize_basis_against_basis_csp

    module procedure orthogonalize_vector_against_basis_cdp
        character(len=*), parameter :: this_procedure = 'orthogonalize_vector_against_basis_cdp'
        complex(dp) :: proj_coefficients(size(X))
        ! internal
        logical :: chk_X_orthonormality
        logical :: is_ortho

        if (time_lightkrylov()) call timer%start(this_procedure)
        info = 0

        ! optional input argument
        chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

        ! check for zero vector
        if (y%norm() < atol_dp) info = 1

        ! optional orthonormality check
        if (chk_X_orthonormality) then
            is_ortho = is_orthonormal_cdp(X)
            if (is_ortho) then
                call log_information("Input basis orthonormal. Remove this check unless necessary for better performance", &
                & this_module, this_procedure)
            else
                call stop_error("Input basis not orthonormal.", this_module, this_procedure)
            end if
        end if

        ! orthogonalize
        proj_coefficients = innerprod(X, y)
        ! y <- y - X c, subtracted term by term: no temporary vector is needed.
        block
            integer :: i
            do i = 1, size(X)
                call y%axpby(-proj_coefficients(i), X(i), one_cdp)
            enddo
        end block

        if (present(beta)) then
            ! check size
            call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, this_procedure)
            beta = proj_coefficients
        end if
        if (time_lightkrylov()) call timer%stop(this_procedure)
    end procedure orthogonalize_vector_against_basis_cdp

    module procedure orthogonalize_basis_against_basis_cdp
        character(len=*), parameter :: this_procedure = 'orthogonalize_basis_against_basis_cdp'
        complex(dp) :: proj_coefficients(size(X), size(Y))
        ! internal
        integer :: i
        logical :: chk_X_orthonormality
        logical :: is_ortho

        if (time_lightkrylov()) call timer%start(this_procedure)
        info = 0

        ! optional input argument
        chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

        ! check for zero vector
        do i = 1, size(Y)
            if (Y(i)%norm() < atol_dp) info = i
        end do

        ! optional orthonormality check
        if (chk_X_orthonormality) then
            is_ortho = is_orthonormal_cdp(X)
            if (is_ortho) then
                call log_information("Input basis orthonormal. Remove this check unless necessary for better performance", &
                & this_module, this_procedure)
            else
                call stop_error("Input basis not orthonormal.", this_module, this_procedure)
            end if
        end if

        ! orthogonalize
        proj_coefficients = innerprod(X, Y)
        ! Y <- Y - X C, subtracted term by term: no temporary basis is needed.
        block
            integer :: ii, jj
            do jj = 1, size(Y)
                do ii = 1, size(X)
                    call Y(jj)%axpby(-proj_coefficients(ii, jj), X(ii), one_cdp)
                enddo
            enddo
        end block

        if (present(beta)) then
            ! check size
            call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, this_procedure)
            beta = proj_coefficients
        end if
        if (time_lightkrylov()) call timer%stop(this_procedure)
    end procedure orthogonalize_basis_against_basis_cdp

end submodule gram_schmidt_process