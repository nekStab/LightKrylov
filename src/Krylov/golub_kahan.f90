submodule (lightkrylov_basekrylov) golub_kahan_methods
    implicit none
contains
    module procedure lanczos_bidiagonalization_rsp
        character(len=*), parameter :: this_procedure = 'lanczos_bidiagonalization_rsp'
        integer :: k_start, k_end
        real(sp) :: tolerance
        real(sp) :: alpha, beta
        integer :: k, kdim

        if (time_lightkrylov()) call timer%start(this_procedure)
        info = 0

        ! Krylov subspace dimension.
        kdim = size(U) - 1

        ! Deals with the optional args.
        k_start = optval(kstart, 1)
        k_end   = optval(kend, kdim)
        tolerance = optval(tol, atol_sp)

        ! Lanczos bidiagonalization.
        lanczos : do k = k_start, k_end
            ! Transpose matrix-vector product.
            call A%apply_rmatvec(U(k), V(k))

            ! Full re-orthogonalization of the right Krylov subspace.
            if (k > 1) then
                call double_gram_schmidt_step(V(k), V(:k-1), info, if_chk_orthonormal=.false.)
                call check_info(info, 'double_gram_schmidt_step', this_module, this_procedure//', right basis')
            end if

            ! Normalization step.
            alpha = V(k)%norm() ; B(k, k) = alpha
            if (abs(alpha) > tolerance) then
                call V(k)%scal(one_rsp/alpha)
            else
                info = k
                exit lanczos
            endif

            ! Matrix-vector product.
            call A%apply_matvec(V(k), U(k+1))

            ! Full re-orthogonalization of the left Krylov subspace.
            call double_gram_schmidt_step(U(k+1), U(:k), info, if_chk_orthonormal=.false.)
            call check_info(info, 'double_gram_schmidt_step', this_module, this_procedure//', left basis')

            ! Normalization step
            beta = U(k+1)%norm() ; B(k+1, k) = beta
            if (abs(beta) > tolerance) then
                call U(k+1)%scal(one_rsp / beta)
            else
                info = k
                exit lanczos
            endif

        enddo lanczos

        if (time_lightkrylov()) call timer%stop(this_procedure)
        
        return
    end procedure
    module procedure lanczos_bidiagonalization_rdp
        character(len=*), parameter :: this_procedure = 'lanczos_bidiagonalization_rdp'
        integer :: k_start, k_end
        real(dp) :: tolerance
        real(dp) :: alpha, beta
        integer :: k, kdim

        if (time_lightkrylov()) call timer%start(this_procedure)
        info = 0

        ! Krylov subspace dimension.
        kdim = size(U) - 1

        ! Deals with the optional args.
        k_start = optval(kstart, 1)
        k_end   = optval(kend, kdim)
        tolerance = optval(tol, atol_dp)

        ! Lanczos bidiagonalization.
        lanczos : do k = k_start, k_end
            ! Transpose matrix-vector product.
            call A%apply_rmatvec(U(k), V(k))

            ! Full re-orthogonalization of the right Krylov subspace.
            if (k > 1) then
                call double_gram_schmidt_step(V(k), V(:k-1), info, if_chk_orthonormal=.false.)
                call check_info(info, 'double_gram_schmidt_step', this_module, this_procedure//', right basis')
            end if

            ! Normalization step.
            alpha = V(k)%norm() ; B(k, k) = alpha
            if (abs(alpha) > tolerance) then
                call V(k)%scal(one_rdp/alpha)
            else
                info = k
                exit lanczos
            endif

            ! Matrix-vector product.
            call A%apply_matvec(V(k), U(k+1))

            ! Full re-orthogonalization of the left Krylov subspace.
            call double_gram_schmidt_step(U(k+1), U(:k), info, if_chk_orthonormal=.false.)
            call check_info(info, 'double_gram_schmidt_step', this_module, this_procedure//', left basis')

            ! Normalization step
            beta = U(k+1)%norm() ; B(k+1, k) = beta
            if (abs(beta) > tolerance) then
                call U(k+1)%scal(one_rdp / beta)
            else
                info = k
                exit lanczos
            endif

        enddo lanczos

        if (time_lightkrylov()) call timer%stop(this_procedure)
        
        return
    end procedure
    module procedure lanczos_bidiagonalization_csp
        character(len=*), parameter :: this_procedure = 'lanczos_bidiagonalization_csp'
        integer :: k_start, k_end
        real(sp) :: tolerance
        complex(sp) :: alpha, beta
        integer :: k, kdim

        if (time_lightkrylov()) call timer%start(this_procedure)
        info = 0

        ! Krylov subspace dimension.
        kdim = size(U) - 1

        ! Deals with the optional args.
        k_start = optval(kstart, 1)
        k_end   = optval(kend, kdim)
        tolerance = optval(tol, atol_sp)

        ! Lanczos bidiagonalization.
        lanczos : do k = k_start, k_end
            ! Transpose matrix-vector product.
            call A%apply_rmatvec(U(k), V(k))

            ! Full re-orthogonalization of the right Krylov subspace.
            if (k > 1) then
                call double_gram_schmidt_step(V(k), V(:k-1), info, if_chk_orthonormal=.false.)
                call check_info(info, 'double_gram_schmidt_step', this_module, this_procedure//', right basis')
            end if

            ! Normalization step.
            alpha = V(k)%norm() ; B(k, k) = alpha
            if (abs(alpha) > tolerance) then
                call V(k)%scal(one_csp/alpha)
            else
                info = k
                exit lanczos
            endif

            ! Matrix-vector product.
            call A%apply_matvec(V(k), U(k+1))

            ! Full re-orthogonalization of the left Krylov subspace.
            call double_gram_schmidt_step(U(k+1), U(:k), info, if_chk_orthonormal=.false.)
            call check_info(info, 'double_gram_schmidt_step', this_module, this_procedure//', left basis')

            ! Normalization step
            beta = U(k+1)%norm() ; B(k+1, k) = beta
            if (abs(beta) > tolerance) then
                call U(k+1)%scal(one_csp / beta)
            else
                info = k
                exit lanczos
            endif

        enddo lanczos

        if (time_lightkrylov()) call timer%stop(this_procedure)
        
        return
    end procedure
    module procedure lanczos_bidiagonalization_cdp
        character(len=*), parameter :: this_procedure = 'lanczos_bidiagonalization_cdp'
        integer :: k_start, k_end
        real(dp) :: tolerance
        complex(dp) :: alpha, beta
        integer :: k, kdim

        if (time_lightkrylov()) call timer%start(this_procedure)
        info = 0

        ! Krylov subspace dimension.
        kdim = size(U) - 1

        ! Deals with the optional args.
        k_start = optval(kstart, 1)
        k_end   = optval(kend, kdim)
        tolerance = optval(tol, atol_dp)

        ! Lanczos bidiagonalization.
        lanczos : do k = k_start, k_end
            ! Transpose matrix-vector product.
            call A%apply_rmatvec(U(k), V(k))

            ! Full re-orthogonalization of the right Krylov subspace.
            if (k > 1) then
                call double_gram_schmidt_step(V(k), V(:k-1), info, if_chk_orthonormal=.false.)
                call check_info(info, 'double_gram_schmidt_step', this_module, this_procedure//', right basis')
            end if

            ! Normalization step.
            alpha = V(k)%norm() ; B(k, k) = alpha
            if (abs(alpha) > tolerance) then
                call V(k)%scal(one_cdp/alpha)
            else
                info = k
                exit lanczos
            endif

            ! Matrix-vector product.
            call A%apply_matvec(V(k), U(k+1))

            ! Full re-orthogonalization of the left Krylov subspace.
            call double_gram_schmidt_step(U(k+1), U(:k), info, if_chk_orthonormal=.false.)
            call check_info(info, 'double_gram_schmidt_step', this_module, this_procedure//', left basis')

            ! Normalization step
            beta = U(k+1)%norm() ; B(k+1, k) = beta
            if (abs(beta) > tolerance) then
                call U(k+1)%scal(one_cdp / beta)
            else
                info = k
                exit lanczos
            endif

        enddo lanczos

        if (time_lightkrylov()) call timer%stop(this_procedure)
        
        return
    end procedure
end submodule
