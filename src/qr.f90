submodule (lightkrylov_basekrylov) qr_solvers
    implicit none

    interface swap_columns
        module subroutine swap_columns_rsp(Q, R, Rii, perm, i, j)
            class(abstract_vector_rsp), intent(inout) :: Q(:)
            !! Vector basis whose i-th and j-th columns need swapping.
            real(sp), intent(inout) :: R(:, :)
            !! Upper triangular matrix resulting from QR.
            real(sp), intent(inout) :: Rii(:)
            !! Diagonal entries of R.
            integer, intent(inout) :: perm(:)
            !! Column ordering.
            integer, intent(in) :: i, j
            !! Index of the columns to be swapped.
        end subroutine
        module subroutine swap_columns_rdp(Q, R, Rii, perm, i, j)
            class(abstract_vector_rdp), intent(inout) :: Q(:)
            !! Vector basis whose i-th and j-th columns need swapping.
            real(dp), intent(inout) :: R(:, :)
            !! Upper triangular matrix resulting from QR.
            real(dp), intent(inout) :: Rii(:)
            !! Diagonal entries of R.
            integer, intent(inout) :: perm(:)
            !! Column ordering.
            integer, intent(in) :: i, j
            !! Index of the columns to be swapped.
        end subroutine
        module subroutine swap_columns_csp(Q, R, Rii, perm, i, j)
            class(abstract_vector_csp), intent(inout) :: Q(:)
            !! Vector basis whose i-th and j-th columns need swapping.
            complex(sp), intent(inout) :: R(:, :)
            !! Upper triangular matrix resulting from QR.
            complex(sp), intent(inout) :: Rii(:)
            !! Diagonal entries of R.
            integer, intent(inout) :: perm(:)
            !! Column ordering.
            integer, intent(in) :: i, j
            !! Index of the columns to be swapped.
        end subroutine
        module subroutine swap_columns_cdp(Q, R, Rii, perm, i, j)
            class(abstract_vector_cdp), intent(inout) :: Q(:)
            !! Vector basis whose i-th and j-th columns need swapping.
            complex(dp), intent(inout) :: R(:, :)
            !! Upper triangular matrix resulting from QR.
            complex(dp), intent(inout) :: Rii(:)
            !! Diagonal entries of R.
            integer, intent(inout) :: perm(:)
            !! Column ordering.
            integer, intent(in) :: i, j
            !! Index of the columns to be swapped.
        end subroutine
    end interface

contains

    !------------------------------------
    !-----     QR WITH PIVOTING     -----
    !------------------------------------

    module procedure qr_with_pivoting_rsp
        character(len=*), parameter :: this_procedure = 'qr_with_pivoting_rsp'
        real(sp) :: tolerance
        real(sp) :: beta
        integer :: idx, i, j, kdim
        integer :: idxv(1)
        real(sp)  :: Rii(size(Q))
        character(len=128) :: msg

        if (time_lightkrylov()) call timer%start(this_procedure)
        info = 0 ; kdim = size(Q) ; R = zero_rsp 
        
        ! Deals with the optional arguments.
        tolerance = optval(tol, atol_sp)

        ! Initialize diagonal entries.
        do i = 1, kdim
            perm(i) = i
            Rii(i) = Q(i)%dot(Q(i))
        enddo

        qr_step: do j = 1, kdim
            idxv = maxloc(abs(Rii)) ; idx = idxv(1)
            if (abs(Rii(idx)) < tolerance) then
                do i = j, kdim
                    call Q(i)%rand()
                    call double_gram_schmidt_step(Q(i), Q(:i-1), info, if_chk_orthonormal=.false.)
                    call check_info(info, 'double_gram_schmidt_step', this_module, this_procedure)
                    beta = Q(i)%norm(); call Q(i)%scal(one_rsp / beta)
                enddo
                info = j
                write(msg,'(A,I0,A,E15.8)') 'Breakdown after ', j, ' steps. R_ii= ', abs(Rii(idx))
                call log_information(msg, this_module, this_procedure)
                exit qr_step
            endif

            call swap_columns(Q, R, Rii, perm, j, idx)

            ! Check for breakdown.
            beta = Q(j)%norm()
            if (isnan(beta)) call stop_error('|beta| = NaN detected! Abort', this_module, this_procedure)
            if (abs(beta) < tolerance) then
                info = j
                R(j, j) = zero_rsp
                call Q(j)%rand()
                call double_gram_schmidt_step(Q(j), Q(:j-1), info, if_chk_orthonormal=.false.)
                call check_info(info, 'double_gram_schmidt_step', this_module, this_procedure)
                beta = Q(j)%norm()
            else
                R(j, j) = beta
            endif
            ! Normalize column.
            call Q(j)%scal(one_rsp / beta)

            ! Orthogonalize all columns against new vector.
            do i = j+1, kdim
                beta = Q(j)%dot(Q(i))
                call Q(i)%axpby(-beta, Q(j), one_rsp)   ! Q(i) = Q(i) - beta*Q(j)
                R(j, i) = beta
            enddo

            ! Update Rii.
            Rii(j) = zero_rsp
            do i = j+1, kdim
                Rii(i) = Rii(i) - R(j, i)**2
            enddo

        enddo qr_step
        if (time_lightkrylov()) call timer%stop(this_procedure)
        
        return
    end procedure
    module procedure qr_with_pivoting_rdp
        character(len=*), parameter :: this_procedure = 'qr_with_pivoting_rdp'
        real(dp) :: tolerance
        real(dp) :: beta
        integer :: idx, i, j, kdim
        integer :: idxv(1)
        real(dp)  :: Rii(size(Q))
        character(len=128) :: msg

        if (time_lightkrylov()) call timer%start(this_procedure)
        info = 0 ; kdim = size(Q) ; R = zero_rdp 
        
        ! Deals with the optional arguments.
        tolerance = optval(tol, atol_dp)

        ! Initialize diagonal entries.
        do i = 1, kdim
            perm(i) = i
            Rii(i) = Q(i)%dot(Q(i))
        enddo

        qr_step: do j = 1, kdim
            idxv = maxloc(abs(Rii)) ; idx = idxv(1)
            if (abs(Rii(idx)) < tolerance) then
                do i = j, kdim
                    call Q(i)%rand()
                    call double_gram_schmidt_step(Q(i), Q(:i-1), info, if_chk_orthonormal=.false.)
                    call check_info(info, 'double_gram_schmidt_step', this_module, this_procedure)
                    beta = Q(i)%norm(); call Q(i)%scal(one_rdp / beta)
                enddo
                info = j
                write(msg,'(A,I0,A,E15.8)') 'Breakdown after ', j, ' steps. R_ii= ', abs(Rii(idx))
                call log_information(msg, this_module, this_procedure)
                exit qr_step
            endif

            call swap_columns(Q, R, Rii, perm, j, idx)

            ! Check for breakdown.
            beta = Q(j)%norm()
            if (isnan(beta)) call stop_error('|beta| = NaN detected! Abort', this_module, this_procedure)
            if (abs(beta) < tolerance) then
                info = j
                R(j, j) = zero_rdp
                call Q(j)%rand()
                call double_gram_schmidt_step(Q(j), Q(:j-1), info, if_chk_orthonormal=.false.)
                call check_info(info, 'double_gram_schmidt_step', this_module, this_procedure)
                beta = Q(j)%norm()
            else
                R(j, j) = beta
            endif
            ! Normalize column.
            call Q(j)%scal(one_rdp / beta)

            ! Orthogonalize all columns against new vector.
            do i = j+1, kdim
                beta = Q(j)%dot(Q(i))
                call Q(i)%axpby(-beta, Q(j), one_rdp)   ! Q(i) = Q(i) - beta*Q(j)
                R(j, i) = beta
            enddo

            ! Update Rii.
            Rii(j) = zero_rdp
            do i = j+1, kdim
                Rii(i) = Rii(i) - R(j, i)**2
            enddo

        enddo qr_step
        if (time_lightkrylov()) call timer%stop(this_procedure)
        
        return
    end procedure
    module procedure qr_with_pivoting_csp
        character(len=*), parameter :: this_procedure = 'qr_with_pivoting_csp'
        real(sp) :: tolerance
        complex(sp) :: beta
        integer :: idx, i, j, kdim
        integer :: idxv(1)
        complex(sp)  :: Rii(size(Q))
        character(len=128) :: msg

        if (time_lightkrylov()) call timer%start(this_procedure)
        info = 0 ; kdim = size(Q) ; R = zero_rsp 
        
        ! Deals with the optional arguments.
        tolerance = optval(tol, atol_sp)

        ! Initialize diagonal entries.
        do i = 1, kdim
            perm(i) = i
            Rii(i) = Q(i)%dot(Q(i))
        enddo

        qr_step: do j = 1, kdim
            idxv = maxloc(abs(Rii)) ; idx = idxv(1)
            if (abs(Rii(idx)) < tolerance) then
                do i = j, kdim
                    call Q(i)%rand()
                    call double_gram_schmidt_step(Q(i), Q(:i-1), info, if_chk_orthonormal=.false.)
                    call check_info(info, 'double_gram_schmidt_step', this_module, this_procedure)
                    beta = Q(i)%norm(); call Q(i)%scal(one_csp / beta)
                enddo
                info = j
                write(msg,'(A,I0,A,E15.8)') 'Breakdown after ', j, ' steps. R_ii= ', abs(Rii(idx))
                call log_information(msg, this_module, this_procedure)
                exit qr_step
            endif

            call swap_columns(Q, R, Rii, perm, j, idx)

            ! Check for breakdown.
            beta = Q(j)%norm()
            if (isnan(abs(beta))) call stop_error('|beta| = NaN detected! Abort', this_module, this_procedure)
            if (abs(beta) < tolerance) then
                info = j
                R(j, j) = zero_rsp
                call Q(j)%rand()
                call double_gram_schmidt_step(Q(j), Q(:j-1), info, if_chk_orthonormal=.false.)
                call check_info(info, 'double_gram_schmidt_step', this_module, this_procedure)
                beta = Q(j)%norm()
            else
                R(j, j) = beta
            endif
            ! Normalize column.
            call Q(j)%scal(one_rsp / beta)

            ! Orthogonalize all columns against new vector.
            do i = j+1, kdim
                beta = Q(j)%dot(Q(i))
                call Q(i)%axpby(-beta, Q(j), one_csp)   ! Q(i) = Q(i) - beta*Q(j)
                R(j, i) = beta
            enddo

            ! Update Rii.
            Rii(j) = zero_rsp
            do i = j+1, kdim
                Rii(i) = Rii(i) - R(j, i)**2
            enddo

        enddo qr_step
        if (time_lightkrylov()) call timer%stop(this_procedure)
        
        return
    end procedure
    module procedure qr_with_pivoting_cdp
        character(len=*), parameter :: this_procedure = 'qr_with_pivoting_cdp'
        real(dp) :: tolerance
        complex(dp) :: beta
        integer :: idx, i, j, kdim
        integer :: idxv(1)
        complex(dp)  :: Rii(size(Q))
        character(len=128) :: msg

        if (time_lightkrylov()) call timer%start(this_procedure)
        info = 0 ; kdim = size(Q) ; R = zero_rdp 
        
        ! Deals with the optional arguments.
        tolerance = optval(tol, atol_dp)

        ! Initialize diagonal entries.
        do i = 1, kdim
            perm(i) = i
            Rii(i) = Q(i)%dot(Q(i))
        enddo

        qr_step: do j = 1, kdim
            idxv = maxloc(abs(Rii)) ; idx = idxv(1)
            if (abs(Rii(idx)) < tolerance) then
                do i = j, kdim
                    call Q(i)%rand()
                    call double_gram_schmidt_step(Q(i), Q(:i-1), info, if_chk_orthonormal=.false.)
                    call check_info(info, 'double_gram_schmidt_step', this_module, this_procedure)
                    beta = Q(i)%norm(); call Q(i)%scal(one_cdp / beta)
                enddo
                info = j
                write(msg,'(A,I0,A,E15.8)') 'Breakdown after ', j, ' steps. R_ii= ', abs(Rii(idx))
                call log_information(msg, this_module, this_procedure)
                exit qr_step
            endif

            call swap_columns(Q, R, Rii, perm, j, idx)

            ! Check for breakdown.
            beta = Q(j)%norm()
            if (isnan(abs(beta))) call stop_error('|beta| = NaN detected! Abort', this_module, this_procedure)
            if (abs(beta) < tolerance) then
                info = j
                R(j, j) = zero_rdp
                call Q(j)%rand()
                call double_gram_schmidt_step(Q(j), Q(:j-1), info, if_chk_orthonormal=.false.)
                call check_info(info, 'double_gram_schmidt_step', this_module, this_procedure)
                beta = Q(j)%norm()
            else
                R(j, j) = beta
            endif
            ! Normalize column.
            call Q(j)%scal(one_rdp / beta)

            ! Orthogonalize all columns against new vector.
            do i = j+1, kdim
                beta = Q(j)%dot(Q(i))
                call Q(i)%axpby(-beta, Q(j), one_cdp)   ! Q(i) = Q(i) - beta*Q(j)
                R(j, i) = beta
            enddo

            ! Update Rii.
            Rii(j) = zero_rdp
            do i = j+1, kdim
                Rii(i) = Rii(i) - R(j, i)**2
            enddo

        enddo qr_step
        if (time_lightkrylov()) call timer%stop(this_procedure)
        
        return
    end procedure

    !---------------------------------------------
    !-----     STANDARD QR FACTORIZATION     -----
    !---------------------------------------------

    module procedure qr_no_pivoting_rsp
        character(len=*), parameter :: this_procedure = 'qr_no_pivoting_rsp'
        real(sp) :: tolerance
        real(sp) :: beta
        integer :: j
        logical :: flag
        character(len=128) :: msg

        ! Deals with the optional args.
        tolerance = optval(tol, atol_sp)

        if (time_lightkrylov()) call timer%start(this_procedure)
        info = 0 ; flag = .false.; R = zero_rsp ; beta = zero_rsp
        do j = 1, size(Q)
            if (j > 1) then
                ! Double Gram-Schmidt orthogonalization
                call double_gram_schmidt_step(Q(j), Q(:j-1), info, if_chk_orthonormal=.false., beta = R(:j-1,j))
                call check_info(info, 'double_gram_schmidt_step', this_module, this_procedure)
            end if        

            ! Check for breakdown.
            beta = Q(j)%norm()
            if (isnan(beta)) call stop_error('|beta| = NaN detected! Abort', this_module, this_procedure)
            if (abs(beta) < tolerance) then
                if (.not.flag) then
                    flag = .true.
                    info = j
                    write(msg,'(A,I0,A,E15.8)') 'Colinear column detected after ', j, ' steps. beta= ', abs(beta)
                    call log_information(msg, this_module, this_procedure)
                end if
                R(j, j) = zero_rsp
                call Q(j)%rand()
                if (j > 1) then
                    call double_gram_schmidt_step(Q(j), Q(:j-1), info, if_chk_orthonormal=.false.)
                    call check_info(info, 'double_gram_schmidt_step', this_module, this_procedure)
                end if
                beta = Q(j)%norm()
            else
                R(j, j) = beta
            endif
            ! Normalize column.
            call Q(j)%scal(one_rsp / beta)
        enddo
        if (time_lightkrylov()) call timer%stop(this_procedure)
        
        return
    end procedure
    module procedure qr_no_pivoting_rdp
        character(len=*), parameter :: this_procedure = 'qr_no_pivoting_rdp'
        real(dp) :: tolerance
        real(dp) :: beta
        integer :: j
        logical :: flag
        character(len=128) :: msg

        ! Deals with the optional args.
        tolerance = optval(tol, atol_dp)

        if (time_lightkrylov()) call timer%start(this_procedure)
        info = 0 ; flag = .false.; R = zero_rdp ; beta = zero_rdp
        do j = 1, size(Q)
            if (j > 1) then
                ! Double Gram-Schmidt orthogonalization
                call double_gram_schmidt_step(Q(j), Q(:j-1), info, if_chk_orthonormal=.false., beta = R(:j-1,j))
                call check_info(info, 'double_gram_schmidt_step', this_module, this_procedure)
            end if        

            ! Check for breakdown.
            beta = Q(j)%norm()
            if (isnan(beta)) call stop_error('|beta| = NaN detected! Abort', this_module, this_procedure)
            if (abs(beta) < tolerance) then
                if (.not.flag) then
                    flag = .true.
                    info = j
                    write(msg,'(A,I0,A,E15.8)') 'Colinear column detected after ', j, ' steps. beta= ', abs(beta)
                    call log_information(msg, this_module, this_procedure)
                end if
                R(j, j) = zero_rdp
                call Q(j)%rand()
                if (j > 1) then
                    call double_gram_schmidt_step(Q(j), Q(:j-1), info, if_chk_orthonormal=.false.)
                    call check_info(info, 'double_gram_schmidt_step', this_module, this_procedure)
                end if
                beta = Q(j)%norm()
            else
                R(j, j) = beta
            endif
            ! Normalize column.
            call Q(j)%scal(one_rdp / beta)
        enddo
        if (time_lightkrylov()) call timer%stop(this_procedure)
        
        return
    end procedure
    module procedure qr_no_pivoting_csp
        character(len=*), parameter :: this_procedure = 'qr_no_pivoting_csp'
        real(sp) :: tolerance
        complex(sp) :: beta
        integer :: j
        logical :: flag
        character(len=128) :: msg

        ! Deals with the optional args.
        tolerance = optval(tol, atol_sp)

        if (time_lightkrylov()) call timer%start(this_procedure)
        info = 0 ; flag = .false.; R = zero_rsp ; beta = zero_rsp
        do j = 1, size(Q)
            if (j > 1) then
                ! Double Gram-Schmidt orthogonalization
                call double_gram_schmidt_step(Q(j), Q(:j-1), info, if_chk_orthonormal=.false., beta = R(:j-1,j))
                call check_info(info, 'double_gram_schmidt_step', this_module, this_procedure)
            end if        

            ! Check for breakdown.
            beta = Q(j)%norm()
            if (isnan(abs(beta))) call stop_error('|beta| = NaN detected! Abort', this_module, this_procedure)
            if (abs(beta) < tolerance) then
                if (.not.flag) then
                    flag = .true.
                    info = j
                    write(msg,'(A,I0,A,E15.8)') 'Colinear column detected after ', j, ' steps. beta= ', abs(beta)
                    call log_information(msg, this_module, this_procedure)
                end if
                R(j, j) = zero_rsp
                call Q(j)%rand()
                if (j > 1) then
                    call double_gram_schmidt_step(Q(j), Q(:j-1), info, if_chk_orthonormal=.false.)
                    call check_info(info, 'double_gram_schmidt_step', this_module, this_procedure)
                end if
                beta = Q(j)%norm()
            else
                R(j, j) = beta
            endif
            ! Normalize column.
            call Q(j)%scal(one_rsp / beta)
        enddo
        if (time_lightkrylov()) call timer%stop(this_procedure)
        
        return
    end procedure
    module procedure qr_no_pivoting_cdp
        character(len=*), parameter :: this_procedure = 'qr_no_pivoting_cdp'
        real(dp) :: tolerance
        complex(dp) :: beta
        integer :: j
        logical :: flag
        character(len=128) :: msg

        ! Deals with the optional args.
        tolerance = optval(tol, atol_dp)

        if (time_lightkrylov()) call timer%start(this_procedure)
        info = 0 ; flag = .false.; R = zero_rdp ; beta = zero_rdp
        do j = 1, size(Q)
            if (j > 1) then
                ! Double Gram-Schmidt orthogonalization
                call double_gram_schmidt_step(Q(j), Q(:j-1), info, if_chk_orthonormal=.false., beta = R(:j-1,j))
                call check_info(info, 'double_gram_schmidt_step', this_module, this_procedure)
            end if        

            ! Check for breakdown.
            beta = Q(j)%norm()
            if (isnan(abs(beta))) call stop_error('|beta| = NaN detected! Abort', this_module, this_procedure)
            if (abs(beta) < tolerance) then
                if (.not.flag) then
                    flag = .true.
                    info = j
                    write(msg,'(A,I0,A,E15.8)') 'Colinear column detected after ', j, ' steps. beta= ', abs(beta)
                    call log_information(msg, this_module, this_procedure)
                end if
                R(j, j) = zero_rdp
                call Q(j)%rand()
                if (j > 1) then
                    call double_gram_schmidt_step(Q(j), Q(:j-1), info, if_chk_orthonormal=.false.)
                    call check_info(info, 'double_gram_schmidt_step', this_module, this_procedure)
                end if
                beta = Q(j)%norm()
            else
                R(j, j) = beta
            endif
            ! Normalize column.
            call Q(j)%scal(one_rdp / beta)
        enddo
        if (time_lightkrylov()) call timer%stop(this_procedure)
        
        return
    end procedure

    !-------------------------------------
    !-----     Utility functions     -----
    !-------------------------------------

    module procedure swap_columns_rsp
        class(abstract_vector_rsp), allocatable :: Qwrk
        real(sp), allocatable :: Rwrk(:)
        integer :: iwrk, m, n

        ! Sanity checks.
        m = size(Q) ; n = min(i, j) - 1

        ! Allocations.
        allocate(Qwrk, mold=Q(1)) ; allocate(Rwrk(max(1, n))) ; Rwrk = zero_rsp

        ! Swap columns.
        call copy(Qwrk, Q(j))
        call copy(Q(j), Q(i))
        call copy(Q(i), Qwrk)
        
        Rwrk(1) = Rii(j); Rii(j) = Rii(i); Rii(i) = Rwrk(1)
        iwrk = perm(j); perm(j) = perm(i) ; perm(i) = iwrk

        if (n > 0) then
            Rwrk = R(:n, j) ; R(:n, j) = R(:n, i) ; R(:n, i) = Rwrk
        endif

        return
    end procedure
    module procedure swap_columns_rdp
        class(abstract_vector_rdp), allocatable :: Qwrk
        real(dp), allocatable :: Rwrk(:)
        integer :: iwrk, m, n

        ! Sanity checks.
        m = size(Q) ; n = min(i, j) - 1

        ! Allocations.
        allocate(Qwrk, mold=Q(1)) ; allocate(Rwrk(max(1, n))) ; Rwrk = zero_rdp

        ! Swap columns.
        call copy(Qwrk, Q(j))
        call copy(Q(j), Q(i))
        call copy(Q(i), Qwrk)
        
        Rwrk(1) = Rii(j); Rii(j) = Rii(i); Rii(i) = Rwrk(1)
        iwrk = perm(j); perm(j) = perm(i) ; perm(i) = iwrk

        if (n > 0) then
            Rwrk = R(:n, j) ; R(:n, j) = R(:n, i) ; R(:n, i) = Rwrk
        endif

        return
    end procedure
    module procedure swap_columns_csp
        class(abstract_vector_csp), allocatable :: Qwrk
        complex(sp), allocatable :: Rwrk(:)
        integer :: iwrk, m, n

        ! Sanity checks.
        m = size(Q) ; n = min(i, j) - 1

        ! Allocations.
        allocate(Qwrk, mold=Q(1)) ; allocate(Rwrk(max(1, n))) ; Rwrk = zero_rsp

        ! Swap columns.
        call copy(Qwrk, Q(j))
        call copy(Q(j), Q(i))
        call copy(Q(i), Qwrk)
        
        Rwrk(1) = Rii(j); Rii(j) = Rii(i); Rii(i) = Rwrk(1)
        iwrk = perm(j); perm(j) = perm(i) ; perm(i) = iwrk

        if (n > 0) then
            Rwrk = R(:n, j) ; R(:n, j) = R(:n, i) ; R(:n, i) = Rwrk
        endif

        return
    end procedure
    module procedure swap_columns_cdp
        class(abstract_vector_cdp), allocatable :: Qwrk
        complex(dp), allocatable :: Rwrk(:)
        integer :: iwrk, m, n

        ! Sanity checks.
        m = size(Q) ; n = min(i, j) - 1

        ! Allocations.
        allocate(Qwrk, mold=Q(1)) ; allocate(Rwrk(max(1, n))) ; Rwrk = zero_rdp

        ! Swap columns.
        call copy(Qwrk, Q(j))
        call copy(Q(j), Q(i))
        call copy(Q(i), Qwrk)
        
        Rwrk(1) = Rii(j); Rii(j) = Rii(i); Rii(i) = Rwrk(1)
        iwrk = perm(j); perm(j) = perm(i) ; perm(i) = iwrk

        if (n > 0) then
            Rwrk = R(:n, j) ; R(:n, j) = R(:n, i) ; R(:n, i) = Rwrk
        endif

        return
    end procedure
end submodule
