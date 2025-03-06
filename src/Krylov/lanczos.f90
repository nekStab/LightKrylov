submodule (lightkrylov_basekrylov) lanczos_methods
    implicit none
contains
    module procedure lanczos_tridiagonalization_rsp
        character(len=*), parameter :: this_procedure = 'lanczos_tridiagonalization_rsp'
        integer :: k_start, k_end
        real(sp) :: tolerance
        real(sp) :: beta
        integer :: k, kdim

        if (time_lightkrylov()) call timer%start(this_procedure)

        ! Deal with optional args.
        kdim = size(X) - 1
        k_start = optval(kstart, 1)
        k_end = optval(kend, kdim)
        tolerance = optval(tol, atol_sp)
        info = 0

        ! Lanczos tridiagonalization.
        lanczos: do k = k_start, k_end
            ! Matrix-vector product.
            call A%apply_matvec(X(k), X(k+1))
            ! Update tridiagonal matrix.
            call update_tridiag_matrix_rsp(T, X, k)
            beta = X(k+1)%norm() ; T(k+1, k) = beta

            ! Exit Lanczos loop if needed.
            if (beta < tolerance) then
                ! Dimension of the computed invariant subspace.
                info = k
                ! Exit the Lanczos iteration.
                exit lanczos
            else
                ! Normalize the new Krylov vector.
                call X(k+1)%scal(one_rsp / beta)
            endif
        enddo lanczos

        if (time_lightkrylov()) call timer%stop(this_procedure)
        
        return
    end procedure

    subroutine update_tridiag_matrix_rsp(T, X, k)
        integer, intent(in) :: k
        real(sp), intent(inout) :: T(:, :)
        class(abstract_vector_rsp), intent(inout) :: X(:)

        ! Internal variables.
        integer :: i, info

        info = 0

        ! Orthogonalize residual w.r.t. previously computed Krylov vectors to obtain coefficients in tridiag. matrix
        do i = max(1, k-1), k
            T(i, k) = X(i)%dot(X(k+1)) ; call X(k+1)%axpby(-T(i, k), X(i), one_rsp)
        enddo

        ! Full re-orthogonalization against existing basis
        call double_gram_schmidt_step(X(k+1), X(:k), info, if_chk_orthonormal=.false.)
        call check_info(info, 'orthogonalize_against_basis_p1', this_module, 'update_tridiag_matrix_rsp')

        return
    end subroutine update_tridiag_matrix_rsp
    module procedure lanczos_tridiagonalization_rdp
        character(len=*), parameter :: this_procedure = 'lanczos_tridiagonalization_rdp'
        integer :: k_start, k_end
        real(dp) :: tolerance
        real(dp) :: beta
        integer :: k, kdim

        if (time_lightkrylov()) call timer%start(this_procedure)

        ! Deal with optional args.
        kdim = size(X) - 1
        k_start = optval(kstart, 1)
        k_end = optval(kend, kdim)
        tolerance = optval(tol, atol_dp)
        info = 0

        ! Lanczos tridiagonalization.
        lanczos: do k = k_start, k_end
            ! Matrix-vector product.
            call A%apply_matvec(X(k), X(k+1))
            ! Update tridiagonal matrix.
            call update_tridiag_matrix_rdp(T, X, k)
            beta = X(k+1)%norm() ; T(k+1, k) = beta

            ! Exit Lanczos loop if needed.
            if (beta < tolerance) then
                ! Dimension of the computed invariant subspace.
                info = k
                ! Exit the Lanczos iteration.
                exit lanczos
            else
                ! Normalize the new Krylov vector.
                call X(k+1)%scal(one_rdp / beta)
            endif
        enddo lanczos

        if (time_lightkrylov()) call timer%stop(this_procedure)
        
        return
    end procedure

    subroutine update_tridiag_matrix_rdp(T, X, k)
        integer, intent(in) :: k
        real(dp), intent(inout) :: T(:, :)
        class(abstract_vector_rdp), intent(inout) :: X(:)

        ! Internal variables.
        integer :: i, info

        info = 0

        ! Orthogonalize residual w.r.t. previously computed Krylov vectors to obtain coefficients in tridiag. matrix
        do i = max(1, k-1), k
            T(i, k) = X(i)%dot(X(k+1)) ; call X(k+1)%axpby(-T(i, k), X(i), one_rdp)
        enddo

        ! Full re-orthogonalization against existing basis
        call double_gram_schmidt_step(X(k+1), X(:k), info, if_chk_orthonormal=.false.)
        call check_info(info, 'orthogonalize_against_basis_p1', this_module, 'update_tridiag_matrix_rdp')

        return
    end subroutine update_tridiag_matrix_rdp
    module procedure lanczos_tridiagonalization_csp
        character(len=*), parameter :: this_procedure = 'lanczos_tridiagonalization_csp'
        integer :: k_start, k_end
        real(sp) :: tolerance
        real(sp) :: beta
        integer :: k, kdim

        if (time_lightkrylov()) call timer%start(this_procedure)

        ! Deal with optional args.
        kdim = size(X) - 1
        k_start = optval(kstart, 1)
        k_end = optval(kend, kdim)
        tolerance = optval(tol, atol_sp)
        info = 0

        ! Lanczos tridiagonalization.
        lanczos: do k = k_start, k_end
            ! Matrix-vector product.
            call A%apply_matvec(X(k), X(k+1))
            ! Update tridiagonal matrix.
            call update_tridiag_matrix_csp(T, X, k)
            beta = X(k+1)%norm() ; T(k+1, k) = beta

            ! Exit Lanczos loop if needed.
            if (beta < tolerance) then
                ! Dimension of the computed invariant subspace.
                info = k
                ! Exit the Lanczos iteration.
                exit lanczos
            else
                ! Normalize the new Krylov vector.
                call X(k+1)%scal(one_csp / beta)
            endif
        enddo lanczos

        if (time_lightkrylov()) call timer%stop(this_procedure)
        
        return
    end procedure

    subroutine update_tridiag_matrix_csp(T, X, k)
        integer, intent(in) :: k
        complex(sp), intent(inout) :: T(:, :)
        class(abstract_vector_csp), intent(inout) :: X(:)

        ! Internal variables.
        integer :: i, info

        info = 0

        ! Orthogonalize residual w.r.t. previously computed Krylov vectors to obtain coefficients in tridiag. matrix
        do i = max(1, k-1), k
            T(i, k) = X(i)%dot(X(k+1)) ; call X(k+1)%axpby(-T(i, k), X(i), one_csp)
        enddo

        ! Full re-orthogonalization against existing basis
        call double_gram_schmidt_step(X(k+1), X(:k), info, if_chk_orthonormal=.false.)
        call check_info(info, 'orthogonalize_against_basis_p1', this_module, 'update_tridiag_matrix_csp')

        return
    end subroutine update_tridiag_matrix_csp
    module procedure lanczos_tridiagonalization_cdp
        character(len=*), parameter :: this_procedure = 'lanczos_tridiagonalization_cdp'
        integer :: k_start, k_end
        real(dp) :: tolerance
        real(dp) :: beta
        integer :: k, kdim

        if (time_lightkrylov()) call timer%start(this_procedure)

        ! Deal with optional args.
        kdim = size(X) - 1
        k_start = optval(kstart, 1)
        k_end = optval(kend, kdim)
        tolerance = optval(tol, atol_dp)
        info = 0

        ! Lanczos tridiagonalization.
        lanczos: do k = k_start, k_end
            ! Matrix-vector product.
            call A%apply_matvec(X(k), X(k+1))
            ! Update tridiagonal matrix.
            call update_tridiag_matrix_cdp(T, X, k)
            beta = X(k+1)%norm() ; T(k+1, k) = beta

            ! Exit Lanczos loop if needed.
            if (beta < tolerance) then
                ! Dimension of the computed invariant subspace.
                info = k
                ! Exit the Lanczos iteration.
                exit lanczos
            else
                ! Normalize the new Krylov vector.
                call X(k+1)%scal(one_cdp / beta)
            endif
        enddo lanczos

        if (time_lightkrylov()) call timer%stop(this_procedure)
        
        return
    end procedure

    subroutine update_tridiag_matrix_cdp(T, X, k)
        integer, intent(in) :: k
        complex(dp), intent(inout) :: T(:, :)
        class(abstract_vector_cdp), intent(inout) :: X(:)

        ! Internal variables.
        integer :: i, info

        info = 0

        ! Orthogonalize residual w.r.t. previously computed Krylov vectors to obtain coefficients in tridiag. matrix
        do i = max(1, k-1), k
            T(i, k) = X(i)%dot(X(k+1)) ; call X(k+1)%axpby(-T(i, k), X(i), one_cdp)
        enddo

        ! Full re-orthogonalization against existing basis
        call double_gram_schmidt_step(X(k+1), X(:k), info, if_chk_orthonormal=.false.)
        call check_info(info, 'orthogonalize_against_basis_p1', this_module, 'update_tridiag_matrix_cdp')

        return
    end subroutine update_tridiag_matrix_cdp
end submodule

