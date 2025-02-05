submodule (lightkrylov_basekrylov) arnoldi_method
    implicit none
contains
    
    module procedure arnoldi_rsp
        character(len=*), parameter :: this_procedure = 'arnoldi_rsp'
        integer :: k_start, k_end, p
        logical :: trans
        real(sp) :: tolerance
        real(sp) :: beta
        real(sp), allocatable :: res(:)
        integer :: k, i, kdim, kpm, kp, kpp

        if (time_lightkrylov()) call timer%start(this_procedure)

        ! Deals with optional non-unity blksize and allocations.
        p = optval(blksize, 1) ; allocate(res(p)) ; res = zero_rsp ; info = 0

        ! Check dimensions.
        kdim = (size(X) - p) / p

        ! Deal with the other optional args.
        k_start = optval(kstart, 1) ; k_end = optval(kend, kdim)
        tolerance = optval (tol, atol_sp)
        trans     = optval(transpose, .false.)

        ! Arnoldi factorization.
        blk_arnoldi: do k = k_start, k_end
            ! Counters
            kpm = (k - 1) * p ; kp = kpm + p ; kpp = kp + p

            ! Matrix-vector product.
            if (trans) then
                do i = 1, p
                    call A%apply_rmatvec(X(kpm+i), X(kp+i))
                enddo
            else
                do i = 1, p
                    call A%apply_matvec(X(kpm+i), X(kp+i))
                enddo
            endif

            ! Update Hessenberg matrix via batch double Gram-Schmidt step.
            call double_gram_schmidt_step(X(kp+1:kpp), X(:kp), info, if_chk_orthonormal=.false., beta=H(:kp, kpm+1:kp))
            call check_info(info, 'double_gram_schmidt_step', this_module, this_procedure)

            ! Orthogonalize current blk vectors.
            call qr(X(kp+1:kpp), H(kp+1:kpp, kpm+1:kp), info)
            call check_info(info, 'qr', this_module, this_procedure)

            ! Extract residual norm (smallest diagonal element of H matrix).
            res = zero_rsp
            do i = 1, p
                res(i) = H(kp+i, kpm+i)
            enddo
            beta = minval(abs(res))

            ! Exit Arnoldi loop if needed.
            if (beta < tolerance) then
                ! Dimension of the computed invariant subspace.
                info = kp
                ! Exit the Arnoldi iteration.
                exit blk_arnoldi
            endif

        enddo blk_arnoldi

        if (time_lightkrylov()) call timer%stop(this_procedure)
        
        return
    end procedure
    module procedure arnoldi_rdp
        character(len=*), parameter :: this_procedure = 'arnoldi_rdp'
        integer :: k_start, k_end, p
        logical :: trans
        real(dp) :: tolerance
        real(dp) :: beta
        real(dp), allocatable :: res(:)
        integer :: k, i, kdim, kpm, kp, kpp

        if (time_lightkrylov()) call timer%start(this_procedure)

        ! Deals with optional non-unity blksize and allocations.
        p = optval(blksize, 1) ; allocate(res(p)) ; res = zero_rdp ; info = 0

        ! Check dimensions.
        kdim = (size(X) - p) / p

        ! Deal with the other optional args.
        k_start = optval(kstart, 1) ; k_end = optval(kend, kdim)
        tolerance = optval (tol, atol_dp)
        trans     = optval(transpose, .false.)

        ! Arnoldi factorization.
        blk_arnoldi: do k = k_start, k_end
            ! Counters
            kpm = (k - 1) * p ; kp = kpm + p ; kpp = kp + p

            ! Matrix-vector product.
            if (trans) then
                do i = 1, p
                    call A%apply_rmatvec(X(kpm+i), X(kp+i))
                enddo
            else
                do i = 1, p
                    call A%apply_matvec(X(kpm+i), X(kp+i))
                enddo
            endif

            ! Update Hessenberg matrix via batch double Gram-Schmidt step.
            call double_gram_schmidt_step(X(kp+1:kpp), X(:kp), info, if_chk_orthonormal=.false., beta=H(:kp, kpm+1:kp))
            call check_info(info, 'double_gram_schmidt_step', this_module, this_procedure)

            ! Orthogonalize current blk vectors.
            call qr(X(kp+1:kpp), H(kp+1:kpp, kpm+1:kp), info)
            call check_info(info, 'qr', this_module, this_procedure)

            ! Extract residual norm (smallest diagonal element of H matrix).
            res = zero_rdp
            do i = 1, p
                res(i) = H(kp+i, kpm+i)
            enddo
            beta = minval(abs(res))

            ! Exit Arnoldi loop if needed.
            if (beta < tolerance) then
                ! Dimension of the computed invariant subspace.
                info = kp
                ! Exit the Arnoldi iteration.
                exit blk_arnoldi
            endif

        enddo blk_arnoldi

        if (time_lightkrylov()) call timer%stop(this_procedure)
        
        return
    end procedure
    module procedure arnoldi_csp
        character(len=*), parameter :: this_procedure = 'arnoldi_csp'
        integer :: k_start, k_end, p
        logical :: trans
        real(sp) :: tolerance
        real(sp) :: beta
        complex(sp), allocatable :: res(:)
        integer :: k, i, kdim, kpm, kp, kpp

        if (time_lightkrylov()) call timer%start(this_procedure)

        ! Deals with optional non-unity blksize and allocations.
        p = optval(blksize, 1) ; allocate(res(p)) ; res = zero_rsp ; info = 0

        ! Check dimensions.
        kdim = (size(X) - p) / p

        ! Deal with the other optional args.
        k_start = optval(kstart, 1) ; k_end = optval(kend, kdim)
        tolerance = optval (tol, atol_sp)
        trans     = optval(transpose, .false.)

        ! Arnoldi factorization.
        blk_arnoldi: do k = k_start, k_end
            ! Counters
            kpm = (k - 1) * p ; kp = kpm + p ; kpp = kp + p

            ! Matrix-vector product.
            if (trans) then
                do i = 1, p
                    call A%apply_rmatvec(X(kpm+i), X(kp+i))
                enddo
            else
                do i = 1, p
                    call A%apply_matvec(X(kpm+i), X(kp+i))
                enddo
            endif

            ! Update Hessenberg matrix via batch double Gram-Schmidt step.
            call double_gram_schmidt_step(X(kp+1:kpp), X(:kp), info, if_chk_orthonormal=.false., beta=H(:kp, kpm+1:kp))
            call check_info(info, 'double_gram_schmidt_step', this_module, this_procedure)

            ! Orthogonalize current blk vectors.
            call qr(X(kp+1:kpp), H(kp+1:kpp, kpm+1:kp), info)
            call check_info(info, 'qr', this_module, this_procedure)

            ! Extract residual norm (smallest diagonal element of H matrix).
            res = zero_rsp
            do i = 1, p
                res(i) = H(kp+i, kpm+i)
            enddo
            beta = minval(abs(res))

            ! Exit Arnoldi loop if needed.
            if (beta < tolerance) then
                ! Dimension of the computed invariant subspace.
                info = kp
                ! Exit the Arnoldi iteration.
                exit blk_arnoldi
            endif

        enddo blk_arnoldi

        if (time_lightkrylov()) call timer%stop(this_procedure)
        
        return
    end procedure
    module procedure arnoldi_cdp
        character(len=*), parameter :: this_procedure = 'arnoldi_cdp'
        integer :: k_start, k_end, p
        logical :: trans
        real(dp) :: tolerance
        real(dp) :: beta
        complex(dp), allocatable :: res(:)
        integer :: k, i, kdim, kpm, kp, kpp

        if (time_lightkrylov()) call timer%start(this_procedure)

        ! Deals with optional non-unity blksize and allocations.
        p = optval(blksize, 1) ; allocate(res(p)) ; res = zero_rdp ; info = 0

        ! Check dimensions.
        kdim = (size(X) - p) / p

        ! Deal with the other optional args.
        k_start = optval(kstart, 1) ; k_end = optval(kend, kdim)
        tolerance = optval (tol, atol_dp)
        trans     = optval(transpose, .false.)

        ! Arnoldi factorization.
        blk_arnoldi: do k = k_start, k_end
            ! Counters
            kpm = (k - 1) * p ; kp = kpm + p ; kpp = kp + p

            ! Matrix-vector product.
            if (trans) then
                do i = 1, p
                    call A%apply_rmatvec(X(kpm+i), X(kp+i))
                enddo
            else
                do i = 1, p
                    call A%apply_matvec(X(kpm+i), X(kp+i))
                enddo
            endif

            ! Update Hessenberg matrix via batch double Gram-Schmidt step.
            call double_gram_schmidt_step(X(kp+1:kpp), X(:kp), info, if_chk_orthonormal=.false., beta=H(:kp, kpm+1:kp))
            call check_info(info, 'double_gram_schmidt_step', this_module, this_procedure)

            ! Orthogonalize current blk vectors.
            call qr(X(kp+1:kpp), H(kp+1:kpp, kpm+1:kp), info)
            call check_info(info, 'qr', this_module, this_procedure)

            ! Extract residual norm (smallest diagonal element of H matrix).
            res = zero_rdp
            do i = 1, p
                res(i) = H(kp+i, kpm+i)
            enddo
            beta = minval(abs(res))

            ! Exit Arnoldi loop if needed.
            if (beta < tolerance) then
                ! Dimension of the computed invariant subspace.
                info = kp
                ! Exit the Arnoldi iteration.
                exit blk_arnoldi
            endif

        enddo blk_arnoldi

        if (time_lightkrylov()) call timer%stop(this_procedure)
        
        return
    end procedure
   
end submodule
