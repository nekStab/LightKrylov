submodule (lightkrylov_basekrylov) ssy_method
    implicit none (type, external)
contains
    module procedure ssy_rsp
        character(len=*), parameter :: this_procedure = "ssy_rsp"
        character(len=100) :: errmsg
        integer :: k_start, k_end, i, k
        real(sp) :: tolerance
        real(sp) :: alpha, beta, gamma

        if (time_lightkrylov()) call timer%start(this_procedure)
        associate(kdim => size(U)-1)

            ! Deals with optional arguments.
            k_start = optval(kstart, 1) ; k_end = optval(kend, kdim)
            tolerance = optval(tol, atol_sp)

            ! Saunders-Simon-Yip tridiagonalization.
            if (k_start == 1) then
                beta = U(1)%norm()
                call U(1)%scal(one_rsp / beta)

                gamma = V(1)%norm()
                call V(1)%scal(one_rsp / gamma)
            endif

            do k = k_start, k_end
                ! Compute q = A @ v[k] - gamma[k]*u[k-1]
                call A%matvec(V(k), U(k+1))
                if (k > 1) call U(k+1)%axpby(-gamma, U(k-1), one_rsp)

                ! Compute alpha[k] = dot(u[k], q)
                alpha = U(k)%dot(U(k+1))
                T(k, k) = alpha

                ! Compute beta[k+1] * u[k+1] = q - alpha[k]*u[k]
                call double_gram_schmidt_step(U(k+1), U(:k), info, if_chk_orthonormal=.false.)
                beta = U(k+1)%norm() ; T(k+1, k) = beta
                if (abs(beta) > tolerance) then
                    call U(k+1)%scal(one_rsp / beta)
                else
                    info = k
                    exit
                endif

                ! Compute gamma[k+1]*v[k+1] = hermitian(A) @ u[k] - beta[k]*v[k-1] - alpha[k]*v[k]
                call A%rmatvec(U(k), V(k+1))
                call double_gram_schmidt_step(V(k+1), V(:k), info, if_chk_orthonormal=.false.)
                ! if (k > 1) call V(k+1)%axpby(-T(k, k-1), V(k-1), one_rsp)
                ! call V(k+1)%axpby(-T(k, k), V(k), one_rsp)
                gamma = V(k+1)%norm() ; T(k, k+1) = gamma
                if (abs(gamma) > tolerance) then
                    call V(k+1)%scal(one_rsp / gamma)
                else
                    info = k
                    exit
                endif
            enddo

        end associate
        if (time_lightkrylov()) call timer%stop(this_procedure)
    end procedure ssy_rsp
    module procedure ssy_rdp
        character(len=*), parameter :: this_procedure = "ssy_rdp"
        character(len=100) :: errmsg
        integer :: k_start, k_end, i, k
        real(dp) :: tolerance
        real(dp) :: alpha, beta, gamma

        if (time_lightkrylov()) call timer%start(this_procedure)
        associate(kdim => size(U)-1)

            ! Deals with optional arguments.
            k_start = optval(kstart, 1) ; k_end = optval(kend, kdim)
            tolerance = optval(tol, atol_dp)

            ! Saunders-Simon-Yip tridiagonalization.
            if (k_start == 1) then
                beta = U(1)%norm()
                call U(1)%scal(one_rdp / beta)

                gamma = V(1)%norm()
                call V(1)%scal(one_rdp / gamma)
            endif

            do k = k_start, k_end
                ! Compute q = A @ v[k] - gamma[k]*u[k-1]
                call A%matvec(V(k), U(k+1))
                if (k > 1) call U(k+1)%axpby(-gamma, U(k-1), one_rdp)

                ! Compute alpha[k] = dot(u[k], q)
                alpha = U(k)%dot(U(k+1))
                T(k, k) = alpha

                ! Compute beta[k+1] * u[k+1] = q - alpha[k]*u[k]
                call double_gram_schmidt_step(U(k+1), U(:k), info, if_chk_orthonormal=.false.)
                beta = U(k+1)%norm() ; T(k+1, k) = beta
                if (abs(beta) > tolerance) then
                    call U(k+1)%scal(one_rdp / beta)
                else
                    info = k
                    exit
                endif

                ! Compute gamma[k+1]*v[k+1] = hermitian(A) @ u[k] - beta[k]*v[k-1] - alpha[k]*v[k]
                call A%rmatvec(U(k), V(k+1))
                call double_gram_schmidt_step(V(k+1), V(:k), info, if_chk_orthonormal=.false.)
                ! if (k > 1) call V(k+1)%axpby(-T(k, k-1), V(k-1), one_rdp)
                ! call V(k+1)%axpby(-T(k, k), V(k), one_rdp)
                gamma = V(k+1)%norm() ; T(k, k+1) = gamma
                if (abs(gamma) > tolerance) then
                    call V(k+1)%scal(one_rdp / gamma)
                else
                    info = k
                    exit
                endif
            enddo

        end associate
        if (time_lightkrylov()) call timer%stop(this_procedure)
    end procedure ssy_rdp
    module procedure ssy_csp
        character(len=*), parameter :: this_procedure = "ssy_csp"
        character(len=100) :: errmsg
        integer :: k_start, k_end, i, k
        real(sp) :: tolerance
        complex(sp) :: alpha, beta, gamma

        if (time_lightkrylov()) call timer%start(this_procedure)
        associate(kdim => size(U)-1)

            ! Deals with optional arguments.
            k_start = optval(kstart, 1) ; k_end = optval(kend, kdim)
            tolerance = optval(tol, atol_sp)

            ! Saunders-Simon-Yip tridiagonalization.
            if (k_start == 1) then
                beta = U(1)%norm()
                call U(1)%scal(one_csp / beta)

                gamma = V(1)%norm()
                call V(1)%scal(one_csp / gamma)
            endif

            do k = k_start, k_end
                ! Compute q = A @ v[k] - gamma[k]*u[k-1]
                call A%matvec(V(k), U(k+1))
                if (k > 1) call U(k+1)%axpby(-gamma, U(k-1), one_csp)

                ! Compute alpha[k] = dot(u[k], q)
                alpha = U(k)%dot(U(k+1))
                T(k, k) = alpha

                ! Compute beta[k+1] * u[k+1] = q - alpha[k]*u[k]
                call double_gram_schmidt_step(U(k+1), U(:k), info, if_chk_orthonormal=.false.)
                beta = U(k+1)%norm() ; T(k+1, k) = beta
                if (abs(beta) > tolerance) then
                    call U(k+1)%scal(one_csp / beta)
                else
                    info = k
                    exit
                endif

                ! Compute gamma[k+1]*v[k+1] = hermitian(A) @ u[k] - beta[k]*v[k-1] - alpha[k]*v[k]
                call A%rmatvec(U(k), V(k+1))
                call double_gram_schmidt_step(V(k+1), V(:k), info, if_chk_orthonormal=.false.)
                ! if (k > 1) call V(k+1)%axpby(-T(k, k-1), V(k-1), one_csp)
                ! call V(k+1)%axpby(-T(k, k), V(k), one_csp)
                gamma = V(k+1)%norm() ; T(k, k+1) = gamma
                if (abs(gamma) > tolerance) then
                    call V(k+1)%scal(one_csp / gamma)
                else
                    info = k
                    exit
                endif
            enddo

        end associate
        if (time_lightkrylov()) call timer%stop(this_procedure)
    end procedure ssy_csp
    module procedure ssy_cdp
        character(len=*), parameter :: this_procedure = "ssy_cdp"
        character(len=100) :: errmsg
        integer :: k_start, k_end, i, k
        real(dp) :: tolerance
        complex(dp) :: alpha, beta, gamma

        if (time_lightkrylov()) call timer%start(this_procedure)
        associate(kdim => size(U)-1)

            ! Deals with optional arguments.
            k_start = optval(kstart, 1) ; k_end = optval(kend, kdim)
            tolerance = optval(tol, atol_dp)

            ! Saunders-Simon-Yip tridiagonalization.
            if (k_start == 1) then
                beta = U(1)%norm()
                call U(1)%scal(one_cdp / beta)

                gamma = V(1)%norm()
                call V(1)%scal(one_cdp / gamma)
            endif

            do k = k_start, k_end
                ! Compute q = A @ v[k] - gamma[k]*u[k-1]
                call A%matvec(V(k), U(k+1))
                if (k > 1) call U(k+1)%axpby(-gamma, U(k-1), one_cdp)

                ! Compute alpha[k] = dot(u[k], q)
                alpha = U(k)%dot(U(k+1))
                T(k, k) = alpha

                ! Compute beta[k+1] * u[k+1] = q - alpha[k]*u[k]
                call double_gram_schmidt_step(U(k+1), U(:k), info, if_chk_orthonormal=.false.)
                beta = U(k+1)%norm() ; T(k+1, k) = beta
                if (abs(beta) > tolerance) then
                    call U(k+1)%scal(one_cdp / beta)
                else
                    info = k
                    exit
                endif

                ! Compute gamma[k+1]*v[k+1] = hermitian(A) @ u[k] - beta[k]*v[k-1] - alpha[k]*v[k]
                call A%rmatvec(U(k), V(k+1))
                call double_gram_schmidt_step(V(k+1), V(:k), info, if_chk_orthonormal=.false.)
                ! if (k > 1) call V(k+1)%axpby(-T(k, k-1), V(k-1), one_cdp)
                ! call V(k+1)%axpby(-T(k, k), V(k), one_cdp)
                gamma = V(k+1)%norm() ; T(k, k+1) = gamma
                if (abs(gamma) > tolerance) then
                    call V(k+1)%scal(one_cdp / gamma)
                else
                    info = k
                    exit
                endif
            enddo

        end associate
        if (time_lightkrylov()) call timer%stop(this_procedure)
    end procedure ssy_cdp
end submodule ssy_method
