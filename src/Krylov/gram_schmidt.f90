submodule (lightkrylov_basekrylov) gram_schmidt_process
    implicit none
contains

    !--------------------------------------------
    !-----     DOUBLE GRAM-SCHMIDT STEP     -----
    !--------------------------------------------

    module procedure DGS_vector_against_basis_rsp
        logical                      :: chk_X_orthonormality
        real(sp), dimension(size(X)) :: proj_coefficients, wrk

        ! optional input argument
        chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

        if (time_lightkrylov()) call timer%start('DGS_vector_against_basis_rsp')
        info = 0

        proj_coefficients = zero_rsp; wrk = zero_rsp

        ! Orthogonalize vector y w.r.t. to Krylov basis X in two passes of GS.
        ! first pass
        call orthogonalize_against_basis(y, X, info, if_chk_orthonormal=.false., beta=proj_coefficients)
        call check_info(info, 'orthogonalize_against_basis_p1', module=this_module, &
                        & procedure='DGS_vector_against_basis_rsp, pass 1')
        ! second pass
        call orthogonalize_against_basis(y, X, info, if_chk_orthonormal=.false., beta=wrk)
        call check_info(info, 'orthogonalize_against_basis_p2', module=this_module, &
                        & procedure='DGS_vector_against_basis_rsp, pass 2')
        ! combine passes
        proj_coefficients = proj_coefficients + wrk

        if (present(beta)) then
            ! check size
            call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, 'DGS_vector_against_basis_rsp')
            beta = proj_coefficients
        end if
        if (time_lightkrylov()) call timer%stop('DGS_vector_against_basis_rsp')

        return
    end procedure

    module procedure DGS_basis_against_basis_rsp
        logical                              :: chk_X_orthonormality
        real(sp), dimension(size(X),size(Y)) :: proj_coefficients, wrk

        ! optional input argument
        chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

        if (time_lightkrylov()) call timer%start('DGS_basis_against_basis_rsp')
        info = 0

        proj_coefficients = zero_rsp; wrk = zero_rsp

        ! Orthogonalize Krylov basis Y w.r.t. to Krylov basis X in two passes of GS.
        ! first pass
        call orthogonalize_against_basis(Y, X, info, if_chk_orthonormal=.false., beta=proj_coefficients)
        call check_info(info, 'orthogonalize_against_basis_p1', module=this_module, procedure='DGS_basis_against_basis_rsp, first&
            & pass')
        ! second pass
        call orthogonalize_against_basis(Y, X, info, if_chk_orthonormal=.false., beta=wrk)
        call check_info(info, 'orthogonalize_against_basis_p2', module=this_module, procedure='DGS_basis_against_basis_rsp, second&
            & pass')
        ! combine passes
        proj_coefficients = proj_coefficients + wrk

        if (present(beta)) then
            ! check size
            call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, 'DGS_basis_against_basis_rsp')
            beta = proj_coefficients
        end if
        if (time_lightkrylov()) call timer%stop('DGS_basis_against_basis_rsp')
        return
    end procedure
    module procedure DGS_vector_against_basis_rdp
        logical                      :: chk_X_orthonormality
        real(dp), dimension(size(X)) :: proj_coefficients, wrk

        ! optional input argument
        chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

        if (time_lightkrylov()) call timer%start('DGS_vector_against_basis_rdp')
        info = 0

        proj_coefficients = zero_rdp; wrk = zero_rdp

        ! Orthogonalize vector y w.r.t. to Krylov basis X in two passes of GS.
        ! first pass
        call orthogonalize_against_basis(y, X, info, if_chk_orthonormal=.false., beta=proj_coefficients)
        call check_info(info, 'orthogonalize_against_basis_p1', module=this_module, &
                        & procedure='DGS_vector_against_basis_rdp, pass 1')
        ! second pass
        call orthogonalize_against_basis(y, X, info, if_chk_orthonormal=.false., beta=wrk)
        call check_info(info, 'orthogonalize_against_basis_p2', module=this_module, &
                        & procedure='DGS_vector_against_basis_rdp, pass 2')
        ! combine passes
        proj_coefficients = proj_coefficients + wrk

        if (present(beta)) then
            ! check size
            call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, 'DGS_vector_against_basis_rdp')
            beta = proj_coefficients
        end if
        if (time_lightkrylov()) call timer%stop('DGS_vector_against_basis_rdp')

        return
    end procedure

    module procedure DGS_basis_against_basis_rdp
        logical                              :: chk_X_orthonormality
        real(dp), dimension(size(X),size(Y)) :: proj_coefficients, wrk

        ! optional input argument
        chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

        if (time_lightkrylov()) call timer%start('DGS_basis_against_basis_rdp')
        info = 0

        proj_coefficients = zero_rdp; wrk = zero_rdp

        ! Orthogonalize Krylov basis Y w.r.t. to Krylov basis X in two passes of GS.
        ! first pass
        call orthogonalize_against_basis(Y, X, info, if_chk_orthonormal=.false., beta=proj_coefficients)
        call check_info(info, 'orthogonalize_against_basis_p1', module=this_module, procedure='DGS_basis_against_basis_rdp, first&
            & pass')
        ! second pass
        call orthogonalize_against_basis(Y, X, info, if_chk_orthonormal=.false., beta=wrk)
        call check_info(info, 'orthogonalize_against_basis_p2', module=this_module, procedure='DGS_basis_against_basis_rdp, second&
            & pass')
        ! combine passes
        proj_coefficients = proj_coefficients + wrk

        if (present(beta)) then
            ! check size
            call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, 'DGS_basis_against_basis_rdp')
            beta = proj_coefficients
        end if
        if (time_lightkrylov()) call timer%stop('DGS_basis_against_basis_rdp')
        return
    end procedure
    module procedure DGS_vector_against_basis_csp
        logical                      :: chk_X_orthonormality
        complex(sp), dimension(size(X)) :: proj_coefficients, wrk

        ! optional input argument
        chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

        if (time_lightkrylov()) call timer%start('DGS_vector_against_basis_csp')
        info = 0

        proj_coefficients = zero_csp; wrk = zero_csp

        ! Orthogonalize vector y w.r.t. to Krylov basis X in two passes of GS.
        ! first pass
        call orthogonalize_against_basis(y, X, info, if_chk_orthonormal=.false., beta=proj_coefficients)
        call check_info(info, 'orthogonalize_against_basis_p1', module=this_module, &
                        & procedure='DGS_vector_against_basis_csp, pass 1')
        ! second pass
        call orthogonalize_against_basis(y, X, info, if_chk_orthonormal=.false., beta=wrk)
        call check_info(info, 'orthogonalize_against_basis_p2', module=this_module, &
                        & procedure='DGS_vector_against_basis_csp, pass 2')
        ! combine passes
        proj_coefficients = proj_coefficients + wrk

        if (present(beta)) then
            ! check size
            call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, 'DGS_vector_against_basis_csp')
            beta = proj_coefficients
        end if
        if (time_lightkrylov()) call timer%stop('DGS_vector_against_basis_csp')

        return
    end procedure

    module procedure DGS_basis_against_basis_csp
        logical                              :: chk_X_orthonormality
        complex(sp), dimension(size(X),size(Y)) :: proj_coefficients, wrk

        ! optional input argument
        chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

        if (time_lightkrylov()) call timer%start('DGS_basis_against_basis_csp')
        info = 0

        proj_coefficients = zero_csp; wrk = zero_csp

        ! Orthogonalize Krylov basis Y w.r.t. to Krylov basis X in two passes of GS.
        ! first pass
        call orthogonalize_against_basis(Y, X, info, if_chk_orthonormal=.false., beta=proj_coefficients)
        call check_info(info, 'orthogonalize_against_basis_p1', module=this_module, procedure='DGS_basis_against_basis_csp, first&
            & pass')
        ! second pass
        call orthogonalize_against_basis(Y, X, info, if_chk_orthonormal=.false., beta=wrk)
        call check_info(info, 'orthogonalize_against_basis_p2', module=this_module, procedure='DGS_basis_against_basis_csp, second&
            & pass')
        ! combine passes
        proj_coefficients = proj_coefficients + wrk

        if (present(beta)) then
            ! check size
            call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, 'DGS_basis_against_basis_csp')
            beta = proj_coefficients
        end if
        if (time_lightkrylov()) call timer%stop('DGS_basis_against_basis_csp')
        return
    end procedure
    module procedure DGS_vector_against_basis_cdp
        logical                      :: chk_X_orthonormality
        complex(dp), dimension(size(X)) :: proj_coefficients, wrk

        ! optional input argument
        chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

        if (time_lightkrylov()) call timer%start('DGS_vector_against_basis_cdp')
        info = 0

        proj_coefficients = zero_cdp; wrk = zero_cdp

        ! Orthogonalize vector y w.r.t. to Krylov basis X in two passes of GS.
        ! first pass
        call orthogonalize_against_basis(y, X, info, if_chk_orthonormal=.false., beta=proj_coefficients)
        call check_info(info, 'orthogonalize_against_basis_p1', module=this_module, &
                        & procedure='DGS_vector_against_basis_cdp, pass 1')
        ! second pass
        call orthogonalize_against_basis(y, X, info, if_chk_orthonormal=.false., beta=wrk)
        call check_info(info, 'orthogonalize_against_basis_p2', module=this_module, &
                        & procedure='DGS_vector_against_basis_cdp, pass 2')
        ! combine passes
        proj_coefficients = proj_coefficients + wrk

        if (present(beta)) then
            ! check size
            call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, 'DGS_vector_against_basis_cdp')
            beta = proj_coefficients
        end if
        if (time_lightkrylov()) call timer%stop('DGS_vector_against_basis_cdp')

        return
    end procedure

    module procedure DGS_basis_against_basis_cdp
        logical                              :: chk_X_orthonormality
        complex(dp), dimension(size(X),size(Y)) :: proj_coefficients, wrk

        ! optional input argument
        chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

        if (time_lightkrylov()) call timer%start('DGS_basis_against_basis_cdp')
        info = 0

        proj_coefficients = zero_cdp; wrk = zero_cdp

        ! Orthogonalize Krylov basis Y w.r.t. to Krylov basis X in two passes of GS.
        ! first pass
        call orthogonalize_against_basis(Y, X, info, if_chk_orthonormal=.false., beta=proj_coefficients)
        call check_info(info, 'orthogonalize_against_basis_p1', module=this_module, procedure='DGS_basis_against_basis_cdp, first&
            & pass')
        ! second pass
        call orthogonalize_against_basis(Y, X, info, if_chk_orthonormal=.false., beta=wrk)
        call check_info(info, 'orthogonalize_against_basis_p2', module=this_module, procedure='DGS_basis_against_basis_cdp, second&
            & pass')
        ! combine passes
        proj_coefficients = proj_coefficients + wrk

        if (present(beta)) then
            ! check size
            call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, 'DGS_basis_against_basis_cdp')
            beta = proj_coefficients
        end if
        if (time_lightkrylov()) call timer%stop('DGS_basis_against_basis_cdp')
        return
    end procedure

    !-------------------------------------------------------------------------------
    !-----     ORTHOGONALIZE VECTOR/BASIS AGAINST ALREADY ORTHOGONAL BASIS     -----
    !-------------------------------------------------------------------------------

    module procedure orthogonalize_vector_against_basis_rsp
        logical  :: chk_X_orthonormality
        real(sp) :: proj_coefficients(size(X))

        if (time_lightkrylov()) call timer%start('orthonormalize_vector_against_basis_rsp')
        info = 0

        ! optional input argument
        chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

        ! check for zero vector
        if (y%norm() < atol_sp) info = 1

        if (chk_X_orthonormality) then
            block 
            real(sp), dimension(size(X), size(X)) :: G
            call innerprod(G, X, X)
            if (abs(G(size(X),size(X))) < rtol_sp) then
                ! The last vector in X is zero, it does not impact orthogonalisation
                info = -2
            else if (mnorm(G - eye(size(X)), "Fro") > rtol_sp) then
                ! The basis is not orthonormal. Cannot orthonormalize.
                info = -1
                return
            end if
            end block
        end if

        ! orthogonalize
        call innerprod(proj_coefficients, X, y)
        block
            class(abstract_vector_rsp), allocatable :: proj
            call linear_combination(proj, X, proj_coefficients)
            call y%sub(proj)
        end block

        if (present(beta)) then
            ! check size
            call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, 'orthogonalize_vector_against_basis_rsp')
            beta = proj_coefficients
        end if
        if (time_lightkrylov()) call timer%stop('orthonormalize_vector_against_basis_rsp')
        return
    end procedure

    module procedure orthogonalize_basis_against_basis_rsp
        logical  :: chk_X_orthonormality
        real(sp) :: proj_coefficients(size(X), size(Y))
        integer  :: i

        if (time_lightkrylov()) call timer%start('orthonormalize_basis_against_basis_rsp')
        info = 0

        ! optional input argument
        chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

        ! check for zero vector
        do i = 1, size(Y)
            if (Y(i)%norm() < atol_sp) info = i
        end do

        if (chk_X_orthonormality) then
            block 
            real(sp), dimension(size(X), size(X)) :: G
            call innerprod(G, X, X)
            if (abs(G(size(X),size(X))) < rtol_sp) then
                ! The last vector in X is zero, it does not impact orthogonalisation
                info = -2
            else if (mnorm(G - eye(size(X)), "Fro") > rtol_sp) then
                ! The basis is not orthonormal. Cannot orthonormalize.
                info = -1
                return
            end if
            end block
        end if

        ! orthogonalize
        call innerprod(proj_coefficients, X, Y)
        block
            class(abstract_vector_rsp), allocatable :: proj(:)
            call linear_combination(proj, X, proj_coefficients)
            call axpby_basis(Y, one_rsp, proj, -one_rsp)
        end block

        if (present(beta)) then
            ! check size
            call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, 'orthogonalize_basis_against_basis_rsp')
            beta = proj_coefficients
        end if
        if (time_lightkrylov()) call timer%stop('orthonormalize_basis_against_basis_rsp')

        return
    end procedure
    module procedure orthogonalize_vector_against_basis_rdp
        logical  :: chk_X_orthonormality
        real(dp) :: proj_coefficients(size(X))

        if (time_lightkrylov()) call timer%start('orthonormalize_vector_against_basis_rdp')
        info = 0

        ! optional input argument
        chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

        ! check for zero vector
        if (y%norm() < atol_dp) info = 1

        if (chk_X_orthonormality) then
            block 
            real(dp), dimension(size(X), size(X)) :: G
            call innerprod(G, X, X)
            if (abs(G(size(X),size(X))) < rtol_dp) then
                ! The last vector in X is zero, it does not impact orthogonalisation
                info = -2
            else if (mnorm(G - eye(size(X)), "Fro") > rtol_dp) then
                ! The basis is not orthonormal. Cannot orthonormalize.
                info = -1
                return
            end if
            end block
        end if

        ! orthogonalize
        call innerprod(proj_coefficients, X, y)
        block
            class(abstract_vector_rdp), allocatable :: proj
            call linear_combination(proj, X, proj_coefficients)
            call y%sub(proj)
        end block

        if (present(beta)) then
            ! check size
            call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, 'orthogonalize_vector_against_basis_rdp')
            beta = proj_coefficients
        end if
        if (time_lightkrylov()) call timer%stop('orthonormalize_vector_against_basis_rdp')
        return
    end procedure

    module procedure orthogonalize_basis_against_basis_rdp
        logical  :: chk_X_orthonormality
        real(dp) :: proj_coefficients(size(X), size(Y))
        integer  :: i

        if (time_lightkrylov()) call timer%start('orthonormalize_basis_against_basis_rdp')
        info = 0

        ! optional input argument
        chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

        ! check for zero vector
        do i = 1, size(Y)
            if (Y(i)%norm() < atol_dp) info = i
        end do

        if (chk_X_orthonormality) then
            block 
            real(dp), dimension(size(X), size(X)) :: G
            call innerprod(G, X, X)
            if (abs(G(size(X),size(X))) < rtol_dp) then
                ! The last vector in X is zero, it does not impact orthogonalisation
                info = -2
            else if (mnorm(G - eye(size(X)), "Fro") > rtol_dp) then
                ! The basis is not orthonormal. Cannot orthonormalize.
                info = -1
                return
            end if
            end block
        end if

        ! orthogonalize
        call innerprod(proj_coefficients, X, Y)
        block
            class(abstract_vector_rdp), allocatable :: proj(:)
            call linear_combination(proj, X, proj_coefficients)
            call axpby_basis(Y, one_rdp, proj, -one_rdp)
        end block

        if (present(beta)) then
            ! check size
            call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, 'orthogonalize_basis_against_basis_rdp')
            beta = proj_coefficients
        end if
        if (time_lightkrylov()) call timer%stop('orthonormalize_basis_against_basis_rdp')

        return
    end procedure
    module procedure orthogonalize_vector_against_basis_csp
        logical  :: chk_X_orthonormality
        complex(sp) :: proj_coefficients(size(X))

        if (time_lightkrylov()) call timer%start('orthonormalize_vector_against_basis_csp')
        info = 0

        ! optional input argument
        chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

        ! check for zero vector
        if (y%norm() < atol_sp) info = 1

        if (chk_X_orthonormality) then
            block 
            complex(sp), dimension(size(X), size(X)) :: G
            call innerprod(G, X, X)
            if (abs(G(size(X),size(X))) < rtol_sp) then
                ! The last vector in X is zero, it does not impact orthogonalisation
                info = -2
            else if (mnorm(G - eye(size(X)), "Fro") > rtol_sp) then
                ! The basis is not orthonormal. Cannot orthonormalize.
                info = -1
                return
            end if
            end block
        end if

        ! orthogonalize
        call innerprod(proj_coefficients, X, y)
        block
            class(abstract_vector_csp), allocatable :: proj
            call linear_combination(proj, X, proj_coefficients)
            call y%sub(proj)
        end block

        if (present(beta)) then
            ! check size
            call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, 'orthogonalize_vector_against_basis_csp')
            beta = proj_coefficients
        end if
        if (time_lightkrylov()) call timer%stop('orthonormalize_vector_against_basis_csp')
        return
    end procedure

    module procedure orthogonalize_basis_against_basis_csp
        logical  :: chk_X_orthonormality
        complex(sp) :: proj_coefficients(size(X), size(Y))
        integer  :: i

        if (time_lightkrylov()) call timer%start('orthonormalize_basis_against_basis_csp')
        info = 0

        ! optional input argument
        chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

        ! check for zero vector
        do i = 1, size(Y)
            if (Y(i)%norm() < atol_sp) info = i
        end do

        if (chk_X_orthonormality) then
            block 
            complex(sp), dimension(size(X), size(X)) :: G
            call innerprod(G, X, X)
            if (abs(G(size(X),size(X))) < rtol_sp) then
                ! The last vector in X is zero, it does not impact orthogonalisation
                info = -2
            else if (mnorm(G - eye(size(X)), "Fro") > rtol_sp) then
                ! The basis is not orthonormal. Cannot orthonormalize.
                info = -1
                return
            end if
            end block
        end if

        ! orthogonalize
        call innerprod(proj_coefficients, X, Y)
        block
            class(abstract_vector_csp), allocatable :: proj(:)
            call linear_combination(proj, X, proj_coefficients)
            call axpby_basis(Y, one_csp, proj, -one_csp)
        end block

        if (present(beta)) then
            ! check size
            call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, 'orthogonalize_basis_against_basis_csp')
            beta = proj_coefficients
        end if
        if (time_lightkrylov()) call timer%stop('orthonormalize_basis_against_basis_csp')

        return
    end procedure
    module procedure orthogonalize_vector_against_basis_cdp
        logical  :: chk_X_orthonormality
        complex(dp) :: proj_coefficients(size(X))

        if (time_lightkrylov()) call timer%start('orthonormalize_vector_against_basis_cdp')
        info = 0

        ! optional input argument
        chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

        ! check for zero vector
        if (y%norm() < atol_dp) info = 1

        if (chk_X_orthonormality) then
            block 
            complex(dp), dimension(size(X), size(X)) :: G
            call innerprod(G, X, X)
            if (abs(G(size(X),size(X))) < rtol_dp) then
                ! The last vector in X is zero, it does not impact orthogonalisation
                info = -2
            else if (mnorm(G - eye(size(X)), "Fro") > rtol_dp) then
                ! The basis is not orthonormal. Cannot orthonormalize.
                info = -1
                return
            end if
            end block
        end if

        ! orthogonalize
        call innerprod(proj_coefficients, X, y)
        block
            class(abstract_vector_cdp), allocatable :: proj
            call linear_combination(proj, X, proj_coefficients)
            call y%sub(proj)
        end block

        if (present(beta)) then
            ! check size
            call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, 'orthogonalize_vector_against_basis_cdp')
            beta = proj_coefficients
        end if
        if (time_lightkrylov()) call timer%stop('orthonormalize_vector_against_basis_cdp')
        return
    end procedure

    module procedure orthogonalize_basis_against_basis_cdp
        logical  :: chk_X_orthonormality
        complex(dp) :: proj_coefficients(size(X), size(Y))
        integer  :: i

        if (time_lightkrylov()) call timer%start('orthonormalize_basis_against_basis_cdp')
        info = 0

        ! optional input argument
        chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

        ! check for zero vector
        do i = 1, size(Y)
            if (Y(i)%norm() < atol_dp) info = i
        end do

        if (chk_X_orthonormality) then
            block 
            complex(dp), dimension(size(X), size(X)) :: G
            call innerprod(G, X, X)
            if (abs(G(size(X),size(X))) < rtol_dp) then
                ! The last vector in X is zero, it does not impact orthogonalisation
                info = -2
            else if (mnorm(G - eye(size(X)), "Fro") > rtol_dp) then
                ! The basis is not orthonormal. Cannot orthonormalize.
                info = -1
                return
            end if
            end block
        end if

        ! orthogonalize
        call innerprod(proj_coefficients, X, Y)
        block
            class(abstract_vector_cdp), allocatable :: proj(:)
            call linear_combination(proj, X, proj_coefficients)
            call axpby_basis(Y, one_cdp, proj, -one_cdp)
        end block

        if (present(beta)) then
            ! check size
            call assert_shape(beta, shape(proj_coefficients), 'beta', this_module, 'orthogonalize_basis_against_basis_cdp')
            beta = proj_coefficients
        end if
        if (time_lightkrylov()) call timer%stop('orthonormalize_basis_against_basis_cdp')

        return
    end procedure

end submodule

