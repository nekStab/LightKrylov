submodule (lightkrylov_basekrylov) qr_solvers
    implicit none(type, external)

    interface swap_columns
        module subroutine swap_columns_rsp(Q, R, Rii, perm, i, j)
            implicit none(type, external)
            class(abstract_vector_rsp), intent(inout) :: Q(:)
            !! Vector basis whose i-th and j-th columns need swapping.
            real(sp), intent(inout) :: R(:, :)
            !! Upper triangular matrix resulting from QR.
            real(sp), intent(inout) :: Rii(:)
            !! Squared diagonal entries of R.
            integer, intent(inout) :: perm(:)
            !! Permutation vector.
            integer, intent(in) :: i, j
            !! Index of the columns to be swapped.
        end subroutine swap_columns_rsp

        module subroutine swap_columns_rdp(Q, R, Rii, perm, i, j)
            implicit none(type, external)
            class(abstract_vector_rdp), intent(inout) :: Q(:)
            !! Vector basis whose i-th and j-th columns need swapping.
            real(dp), intent(inout) :: R(:, :)
            !! Upper triangular matrix resulting from QR.
            real(dp), intent(inout) :: Rii(:)
            !! Squared diagonal entries of R.
            integer, intent(inout) :: perm(:)
            !! Permutation vector.
            integer, intent(in) :: i, j
            !! Index of the columns to be swapped.
        end subroutine swap_columns_rdp

        module subroutine swap_columns_csp(Q, R, Rii, perm, i, j)
            implicit none(type, external)
            class(abstract_vector_csp), intent(inout) :: Q(:)
            !! Vector basis whose i-th and j-th columns need swapping.
            complex(sp), intent(inout) :: R(:, :)
            !! Upper triangular matrix resulting from QR.
            real(sp), intent(inout) :: Rii(:)
            !! Squared diagonal entries of R.
            integer, intent(inout) :: perm(:)
            !! Permutation vector.
            integer, intent(in) :: i, j
            !! Index of the columns to be swapped.
        end subroutine swap_columns_csp

        module subroutine swap_columns_cdp(Q, R, Rii, perm, i, j)
            implicit none(type, external)
            class(abstract_vector_cdp), intent(inout) :: Q(:)
            !! Vector basis whose i-th and j-th columns need swapping.
            complex(dp), intent(inout) :: R(:, :)
            !! Upper triangular matrix resulting from QR.
            real(dp), intent(inout) :: Rii(:)
            !! Squared diagonal entries of R.
            integer, intent(inout) :: perm(:)
            !! Permutation vector.
            integer, intent(in) :: i, j
            !! Index of the columns to be swapped.
        end subroutine swap_columns_cdp

    end interface

contains

    !------------------------------------
    !-----     QR WITH PIVOTING     -----
    !------------------------------------

    module procedure qr_with_pivoting_rsp
        character(len=*), parameter :: this_procedure = 'qr_with_pivoting_rsp'
        real(sp) :: tolerance
        real(sp) :: alpha, beta
        real(sp) :: gamma
        integer :: idx, i, j, kdim, ierr
        integer :: idxv(1)
        real(sp)  :: Rii(size(Q))
        character(len=128) :: msg

        if (time_lightkrylov()) call timer%start(this_procedure)

        kdim = size(Q) ; R = zero_rsp

        ! Deals with the optional arguments.
        tolerance = optval(tol, atol_sp)

        ! Sanity checks.
        if (size(Q) < 1) then
            info = -1
        else if (size(R, 1) < size(Q) .or. size(R, 2) < size(Q)) then
            info = -2
        else if (size(perm) < size(Q)) then
            info = -3
        else if (tolerance < 0) then
            info = -4
        else
            info = 0
        endif

        if (info /= 0) then
            call check_info(info, 'qr_with_pivoting', this_module, this_procedure)
        else
            ! Initialize diagonal entries.
            do i = 1, kdim
                perm(i) = i
                Rii(i) = real(Q(i)%dot(Q(i)), kind=sp)
            enddo

            qr_step: do j = 1, kdim
                idxv = j - 1 + maxloc(Rii(j:kdim)) ; idx = idxv(1)
                call swap_columns(Q, R, Rii, perm, j, idx)

                ! Check for breakdown.
                beta = Q(j)%norm()
    if (isnan(beta)) call stop_error('|' // "beta" // '| = NaN detected! Abort', this_module, this_procedure)

                if (beta < tolerance) then
                    ! Store old value of beta.
                    alpha = beta
                    ! Remaining columns are numerically in span(Q(1:j-1)): regenerate j..kdim.
                    do i = j, kdim
                        call Q(i)%rand()
                        if (i > 1) then
                            call double_gram_schmidt_step(Q(i), Q(:i-1), ierr, if_chk_orthonormal=.false.)
                            call check_info(ierr, 'double_gram_schmidt_step', this_module, this_procedure)
                        endif
                        beta = Q(i)%norm() ; call Q(i)%scal(one_rsp / beta)
                    enddo
                    info = j
                    write(msg,'(A,I0,A,E15.8)') 'Breakdown after ', j, ' steps. |beta|= ', abs(alpha)
                    call log_information(msg, this_module, this_procedure)
                    exit qr_step
                endif

                R(j, j) = beta
                call Q(j)%scal(one_rsp / beta)

                ! Orthogonalize all columns against new vector and update Rii.
                Rii(j) = zero_rsp
                do i = j+1, kdim
                    gamma = Q(j)%dot(Q(i))
                    call Q(i)%axpby(-gamma, Q(j), one_rsp)   ! Q(i) = Q(i) - gamma*Q(j)
                    R(j, i) = gamma
                    Rii(i) = real(Q(i)%dot(Q(i)), kind=sp)
                enddo

            enddo qr_step
        endif

        if (time_lightkrylov()) call timer%stop(this_procedure)
    end procedure qr_with_pivoting_rsp

    module procedure qr_with_pivoting_rdp
        character(len=*), parameter :: this_procedure = 'qr_with_pivoting_rdp'
        real(dp) :: tolerance
        real(dp) :: alpha, beta
        real(dp) :: gamma
        integer :: idx, i, j, kdim, ierr
        integer :: idxv(1)
        real(dp)  :: Rii(size(Q))
        character(len=128) :: msg

        if (time_lightkrylov()) call timer%start(this_procedure)

        kdim = size(Q) ; R = zero_rdp

        ! Deals with the optional arguments.
        tolerance = optval(tol, atol_dp)

        ! Sanity checks.
        if (size(Q) < 1) then
            info = -1
        else if (size(R, 1) < size(Q) .or. size(R, 2) < size(Q)) then
            info = -2
        else if (size(perm) < size(Q)) then
            info = -3
        else if (tolerance < 0) then
            info = -4
        else
            info = 0
        endif

        if (info /= 0) then
            call check_info(info, 'qr_with_pivoting', this_module, this_procedure)
        else
            ! Initialize diagonal entries.
            do i = 1, kdim
                perm(i) = i
                Rii(i) = real(Q(i)%dot(Q(i)), kind=dp)
            enddo

            qr_step: do j = 1, kdim
                idxv = j - 1 + maxloc(Rii(j:kdim)) ; idx = idxv(1)
                call swap_columns(Q, R, Rii, perm, j, idx)

                ! Check for breakdown.
                beta = Q(j)%norm()
    if (isnan(beta)) call stop_error('|' // "beta" // '| = NaN detected! Abort', this_module, this_procedure)

                if (beta < tolerance) then
                    ! Store old value of beta.
                    alpha = beta
                    ! Remaining columns are numerically in span(Q(1:j-1)): regenerate j..kdim.
                    do i = j, kdim
                        call Q(i)%rand()
                        if (i > 1) then
                            call double_gram_schmidt_step(Q(i), Q(:i-1), ierr, if_chk_orthonormal=.false.)
                            call check_info(ierr, 'double_gram_schmidt_step', this_module, this_procedure)
                        endif
                        beta = Q(i)%norm() ; call Q(i)%scal(one_rdp / beta)
                    enddo
                    info = j
                    write(msg,'(A,I0,A,E15.8)') 'Breakdown after ', j, ' steps. |beta|= ', abs(alpha)
                    call log_information(msg, this_module, this_procedure)
                    exit qr_step
                endif

                R(j, j) = beta
                call Q(j)%scal(one_rdp / beta)

                ! Orthogonalize all columns against new vector and update Rii.
                Rii(j) = zero_rdp
                do i = j+1, kdim
                    gamma = Q(j)%dot(Q(i))
                    call Q(i)%axpby(-gamma, Q(j), one_rdp)   ! Q(i) = Q(i) - gamma*Q(j)
                    R(j, i) = gamma
                    Rii(i) = real(Q(i)%dot(Q(i)), kind=dp)
                enddo

            enddo qr_step
        endif

        if (time_lightkrylov()) call timer%stop(this_procedure)
    end procedure qr_with_pivoting_rdp

    module procedure qr_with_pivoting_csp
        character(len=*), parameter :: this_procedure = 'qr_with_pivoting_csp'
        real(sp) :: tolerance
        real(sp) :: alpha, beta
        complex(sp) :: gamma
        integer :: idx, i, j, kdim, ierr
        integer :: idxv(1)
        real(sp)  :: Rii(size(Q))
        character(len=128) :: msg

        if (time_lightkrylov()) call timer%start(this_procedure)

        kdim = size(Q) ; R = zero_rsp

        ! Deals with the optional arguments.
        tolerance = optval(tol, atol_sp)

        ! Sanity checks.
        if (size(Q) < 1) then
            info = -1
        else if (size(R, 1) < size(Q) .or. size(R, 2) < size(Q)) then
            info = -2
        else if (size(perm) < size(Q)) then
            info = -3
        else if (tolerance < 0) then
            info = -4
        else
            info = 0
        endif

        if (info /= 0) then
            call check_info(info, 'qr_with_pivoting', this_module, this_procedure)
        else
            ! Initialize diagonal entries.
            do i = 1, kdim
                perm(i) = i
                Rii(i) = real(Q(i)%dot(Q(i)), kind=sp)
            enddo

            qr_step: do j = 1, kdim
                idxv = j - 1 + maxloc(Rii(j:kdim)) ; idx = idxv(1)
                call swap_columns(Q, R, Rii, perm, j, idx)

                ! Check for breakdown.
                beta = Q(j)%norm()
    if (isnan(beta)) call stop_error('|' // "beta" // '| = NaN detected! Abort', this_module, this_procedure)

                if (beta < tolerance) then
                    ! Store old value of beta.
                    alpha = beta
                    ! Remaining columns are numerically in span(Q(1:j-1)): regenerate j..kdim.
                    do i = j, kdim
                        call Q(i)%rand()
                        if (i > 1) then
                            call double_gram_schmidt_step(Q(i), Q(:i-1), ierr, if_chk_orthonormal=.false.)
                            call check_info(ierr, 'double_gram_schmidt_step', this_module, this_procedure)
                        endif
                        beta = Q(i)%norm() ; call Q(i)%scal(one_csp / beta)
                    enddo
                    info = j
                    write(msg,'(A,I0,A,E15.8)') 'Breakdown after ', j, ' steps. |beta|= ', abs(alpha)
                    call log_information(msg, this_module, this_procedure)
                    exit qr_step
                endif

                R(j, j) = beta
                call Q(j)%scal(one_csp / beta)

                ! Orthogonalize all columns against new vector and update Rii.
                Rii(j) = zero_rsp
                do i = j+1, kdim
                    gamma = Q(j)%dot(Q(i))
                    call Q(i)%axpby(-gamma, Q(j), one_csp)   ! Q(i) = Q(i) - gamma*Q(j)
                    R(j, i) = gamma
                    Rii(i) = real(Q(i)%dot(Q(i)), kind=sp)
                enddo

            enddo qr_step
        endif

        if (time_lightkrylov()) call timer%stop(this_procedure)
    end procedure qr_with_pivoting_csp

    module procedure qr_with_pivoting_cdp
        character(len=*), parameter :: this_procedure = 'qr_with_pivoting_cdp'
        real(dp) :: tolerance
        real(dp) :: alpha, beta
        complex(dp) :: gamma
        integer :: idx, i, j, kdim, ierr
        integer :: idxv(1)
        real(dp)  :: Rii(size(Q))
        character(len=128) :: msg

        if (time_lightkrylov()) call timer%start(this_procedure)

        kdim = size(Q) ; R = zero_rdp

        ! Deals with the optional arguments.
        tolerance = optval(tol, atol_dp)

        ! Sanity checks.
        if (size(Q) < 1) then
            info = -1
        else if (size(R, 1) < size(Q) .or. size(R, 2) < size(Q)) then
            info = -2
        else if (size(perm) < size(Q)) then
            info = -3
        else if (tolerance < 0) then
            info = -4
        else
            info = 0
        endif

        if (info /= 0) then
            call check_info(info, 'qr_with_pivoting', this_module, this_procedure)
        else
            ! Initialize diagonal entries.
            do i = 1, kdim
                perm(i) = i
                Rii(i) = real(Q(i)%dot(Q(i)), kind=dp)
            enddo

            qr_step: do j = 1, kdim
                idxv = j - 1 + maxloc(Rii(j:kdim)) ; idx = idxv(1)
                call swap_columns(Q, R, Rii, perm, j, idx)

                ! Check for breakdown.
                beta = Q(j)%norm()
    if (isnan(beta)) call stop_error('|' // "beta" // '| = NaN detected! Abort', this_module, this_procedure)

                if (beta < tolerance) then
                    ! Store old value of beta.
                    alpha = beta
                    ! Remaining columns are numerically in span(Q(1:j-1)): regenerate j..kdim.
                    do i = j, kdim
                        call Q(i)%rand()
                        if (i > 1) then
                            call double_gram_schmidt_step(Q(i), Q(:i-1), ierr, if_chk_orthonormal=.false.)
                            call check_info(ierr, 'double_gram_schmidt_step', this_module, this_procedure)
                        endif
                        beta = Q(i)%norm() ; call Q(i)%scal(one_cdp / beta)
                    enddo
                    info = j
                    write(msg,'(A,I0,A,E15.8)') 'Breakdown after ', j, ' steps. |beta|= ', abs(alpha)
                    call log_information(msg, this_module, this_procedure)
                    exit qr_step
                endif

                R(j, j) = beta
                call Q(j)%scal(one_cdp / beta)

                ! Orthogonalize all columns against new vector and update Rii.
                Rii(j) = zero_rdp
                do i = j+1, kdim
                    gamma = Q(j)%dot(Q(i))
                    call Q(i)%axpby(-gamma, Q(j), one_cdp)   ! Q(i) = Q(i) - gamma*Q(j)
                    R(j, i) = gamma
                    Rii(i) = real(Q(i)%dot(Q(i)), kind=dp)
                enddo

            enddo qr_step
        endif

        if (time_lightkrylov()) call timer%stop(this_procedure)
    end procedure qr_with_pivoting_cdp


    !---------------------------------------------
    !-----     STANDARD QR FACTORIZATION     -----
    !---------------------------------------------

    module procedure qr_no_pivoting_rsp
        character(len=*), parameter :: this_procedure = 'qr_no_pivoting_rsp'
        real(sp) :: tolerance
        real(sp) :: beta
        real(sp) :: gamma
        integer :: j
        logical :: flag
        character(len=128) :: msg

        if (time_lightkrylov()) call timer%start(this_procedure)

        ! Deals with the optional args.
        tolerance = optval(tol, atol_sp)

        ! Sanity checks.
        if (size(Q) < 1) then
            info = -1
        else if (size(R, 1) < size(Q) .or. size(R, 2) < size(Q)) then
            info = -2
        else if (tolerance < 0) then
            info = -4
        else
            info = 0
        endif
        if (info /= 0) then
            call check_info(info, 'qr_no_pivoting', this_module, this_procedure)
        else
            flag = .false.
            R = zero_rsp
            beta = zero_rsp
            do j = 1, size(Q)
                if (j > 1) then
                    ! Double Gram-Schmidt orthogonalization
                    call double_gram_schmidt_step(Q(j), Q(:j-1), info, &
                                                  if_chk_orthonormal=.false., beta = R(:j-1,j))
                    call check_info(info, 'double_gram_schmidt_step', this_module, this_procedure)
                end if

                ! Check for breakdown.
                beta = Q(j)%norm()
    if (isnan(beta)) call stop_error('|' // "beta" // '| = NaN detected! Abort', this_module, this_procedure)
                if (beta < tolerance) then
                    if (.not.flag) then
                        flag = .true.
                        info = j
                        write(msg,'(A,I0,A,E15.8)') 'Colinear column detected after ', j, ' steps. beta= ', beta
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
        endif

        if (time_lightkrylov()) call timer%stop(this_procedure)
    end procedure qr_no_pivoting_rsp

    module procedure qr_no_pivoting_rdp
        character(len=*), parameter :: this_procedure = 'qr_no_pivoting_rdp'
        real(dp) :: tolerance
        real(dp) :: beta
        real(dp) :: gamma
        integer :: j
        logical :: flag
        character(len=128) :: msg

        if (time_lightkrylov()) call timer%start(this_procedure)

        ! Deals with the optional args.
        tolerance = optval(tol, atol_dp)

        ! Sanity checks.
        if (size(Q) < 1) then
            info = -1
        else if (size(R, 1) < size(Q) .or. size(R, 2) < size(Q)) then
            info = -2
        else if (tolerance < 0) then
            info = -4
        else
            info = 0
        endif
        if (info /= 0) then
            call check_info(info, 'qr_no_pivoting', this_module, this_procedure)
        else
            flag = .false.
            R = zero_rdp
            beta = zero_rdp
            do j = 1, size(Q)
                if (j > 1) then
                    ! Double Gram-Schmidt orthogonalization
                    call double_gram_schmidt_step(Q(j), Q(:j-1), info, &
                                                  if_chk_orthonormal=.false., beta = R(:j-1,j))
                    call check_info(info, 'double_gram_schmidt_step', this_module, this_procedure)
                end if

                ! Check for breakdown.
                beta = Q(j)%norm()
    if (isnan(beta)) call stop_error('|' // "beta" // '| = NaN detected! Abort', this_module, this_procedure)
                if (beta < tolerance) then
                    if (.not.flag) then
                        flag = .true.
                        info = j
                        write(msg,'(A,I0,A,E15.8)') 'Colinear column detected after ', j, ' steps. beta= ', beta
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
        endif

        if (time_lightkrylov()) call timer%stop(this_procedure)
    end procedure qr_no_pivoting_rdp

    module procedure qr_no_pivoting_csp
        character(len=*), parameter :: this_procedure = 'qr_no_pivoting_csp'
        real(sp) :: tolerance
        real(sp) :: beta
        complex(sp) :: gamma
        integer :: j
        logical :: flag
        character(len=128) :: msg

        if (time_lightkrylov()) call timer%start(this_procedure)

        ! Deals with the optional args.
        tolerance = optval(tol, atol_sp)

        ! Sanity checks.
        if (size(Q) < 1) then
            info = -1
        else if (size(R, 1) < size(Q) .or. size(R, 2) < size(Q)) then
            info = -2
        else if (tolerance < 0) then
            info = -4
        else
            info = 0
        endif
        if (info /= 0) then
            call check_info(info, 'qr_no_pivoting', this_module, this_procedure)
        else
            flag = .false.
            R = zero_rsp
            beta = zero_rsp
            do j = 1, size(Q)
                if (j > 1) then
                    ! Double Gram-Schmidt orthogonalization
                    call double_gram_schmidt_step(Q(j), Q(:j-1), info, &
                                                  if_chk_orthonormal=.false., beta = R(:j-1,j))
                    call check_info(info, 'double_gram_schmidt_step', this_module, this_procedure)
                end if

                ! Check for breakdown.
                beta = Q(j)%norm()
    if (isnan(beta)) call stop_error('|' // "beta" // '| = NaN detected! Abort', this_module, this_procedure)
                if (beta < tolerance) then
                    if (.not.flag) then
                        flag = .true.
                        info = j
                        write(msg,'(A,I0,A,E15.8)') 'Colinear column detected after ', j, ' steps. beta= ', beta
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
                call Q(j)%scal(one_csp / beta)
            enddo
        endif

        if (time_lightkrylov()) call timer%stop(this_procedure)
    end procedure qr_no_pivoting_csp

    module procedure qr_no_pivoting_cdp
        character(len=*), parameter :: this_procedure = 'qr_no_pivoting_cdp'
        real(dp) :: tolerance
        real(dp) :: beta
        complex(dp) :: gamma
        integer :: j
        logical :: flag
        character(len=128) :: msg

        if (time_lightkrylov()) call timer%start(this_procedure)

        ! Deals with the optional args.
        tolerance = optval(tol, atol_dp)

        ! Sanity checks.
        if (size(Q) < 1) then
            info = -1
        else if (size(R, 1) < size(Q) .or. size(R, 2) < size(Q)) then
            info = -2
        else if (tolerance < 0) then
            info = -4
        else
            info = 0
        endif
        if (info /= 0) then
            call check_info(info, 'qr_no_pivoting', this_module, this_procedure)
        else
            flag = .false.
            R = zero_rdp
            beta = zero_rdp
            do j = 1, size(Q)
                if (j > 1) then
                    ! Double Gram-Schmidt orthogonalization
                    call double_gram_schmidt_step(Q(j), Q(:j-1), info, &
                                                  if_chk_orthonormal=.false., beta = R(:j-1,j))
                    call check_info(info, 'double_gram_schmidt_step', this_module, this_procedure)
                end if

                ! Check for breakdown.
                beta = Q(j)%norm()
    if (isnan(beta)) call stop_error('|' // "beta" // '| = NaN detected! Abort', this_module, this_procedure)
                if (beta < tolerance) then
                    if (.not.flag) then
                        flag = .true.
                        info = j
                        write(msg,'(A,I0,A,E15.8)') 'Colinear column detected after ', j, ' steps. beta= ', beta
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
                call Q(j)%scal(one_cdp / beta)
            enddo
        endif

        if (time_lightkrylov()) call timer%stop(this_procedure)
    end procedure qr_no_pivoting_cdp


    !-------------------------------------
    !-----     Utility functions     -----
    !-------------------------------------

    module procedure swap_columns_rsp
        class(abstract_vector_rsp), allocatable :: Qwrk
        real(sp), allocatable :: Rwrk(:)
        integer :: iwrk, m, n, iostat
        character(len=100) :: errmsg

        ! Sanity checks.
        m = size(Q) ; n = min(i, j) - 1

        ! Allocations.
        allocate(Qwrk, mold=Q(1), stat=iostat, errmsg=errmsg)
        call check_allocation(iostat, errmsg, this_module, "swap_columns_rsp")
        allocate(Rwrk(max(1, n)), source=zero_rsp, stat=iostat, errmsg=errmsg)
        call check_allocation(iostat, errmsg, this_module, "swap_columns_rsp")

        ! Swap columns.
        call copy(Qwrk, Q(j))
        call copy(Q(j), Q(i))
        call copy(Q(i), Qwrk)

        Rwrk(1) = Rii(j); Rii(j) = Rii(i); Rii(i) = real(Rwrk(1), kind=sp)
        iwrk = perm(j); perm(j) = perm(i) ; perm(i) = iwrk

        if (n > 0) then
            Rwrk = R(:n, j) ; R(:n, j) = R(:n, i) ; R(:n, i) = Rwrk
        endif
    end procedure swap_columns_rsp

    module procedure swap_columns_rdp
        class(abstract_vector_rdp), allocatable :: Qwrk
        real(dp), allocatable :: Rwrk(:)
        integer :: iwrk, m, n, iostat
        character(len=100) :: errmsg

        ! Sanity checks.
        m = size(Q) ; n = min(i, j) - 1

        ! Allocations.
        allocate(Qwrk, mold=Q(1), stat=iostat, errmsg=errmsg)
        call check_allocation(iostat, errmsg, this_module, "swap_columns_rdp")
        allocate(Rwrk(max(1, n)), source=zero_rdp, stat=iostat, errmsg=errmsg)
        call check_allocation(iostat, errmsg, this_module, "swap_columns_rdp")

        ! Swap columns.
        call copy(Qwrk, Q(j))
        call copy(Q(j), Q(i))
        call copy(Q(i), Qwrk)

        Rwrk(1) = Rii(j); Rii(j) = Rii(i); Rii(i) = real(Rwrk(1), kind=dp)
        iwrk = perm(j); perm(j) = perm(i) ; perm(i) = iwrk

        if (n > 0) then
            Rwrk = R(:n, j) ; R(:n, j) = R(:n, i) ; R(:n, i) = Rwrk
        endif
    end procedure swap_columns_rdp

    module procedure swap_columns_csp
        class(abstract_vector_csp), allocatable :: Qwrk
        complex(sp), allocatable :: Rwrk(:)
        integer :: iwrk, m, n, iostat
        character(len=100) :: errmsg

        ! Sanity checks.
        m = size(Q) ; n = min(i, j) - 1

        ! Allocations.
        allocate(Qwrk, mold=Q(1), stat=iostat, errmsg=errmsg)
        call check_allocation(iostat, errmsg, this_module, "swap_columns_csp")
        allocate(Rwrk(max(1, n)), source=zero_csp, stat=iostat, errmsg=errmsg)
        call check_allocation(iostat, errmsg, this_module, "swap_columns_csp")

        ! Swap columns.
        call copy(Qwrk, Q(j))
        call copy(Q(j), Q(i))
        call copy(Q(i), Qwrk)

        Rwrk(1) = Rii(j); Rii(j) = Rii(i); Rii(i) = real(Rwrk(1), kind=sp)
        iwrk = perm(j); perm(j) = perm(i) ; perm(i) = iwrk

        if (n > 0) then
            Rwrk = R(:n, j) ; R(:n, j) = R(:n, i) ; R(:n, i) = Rwrk
        endif
    end procedure swap_columns_csp

    module procedure swap_columns_cdp
        class(abstract_vector_cdp), allocatable :: Qwrk
        complex(dp), allocatable :: Rwrk(:)
        integer :: iwrk, m, n, iostat
        character(len=100) :: errmsg

        ! Sanity checks.
        m = size(Q) ; n = min(i, j) - 1

        ! Allocations.
        allocate(Qwrk, mold=Q(1), stat=iostat, errmsg=errmsg)
        call check_allocation(iostat, errmsg, this_module, "swap_columns_cdp")
        allocate(Rwrk(max(1, n)), source=zero_cdp, stat=iostat, errmsg=errmsg)
        call check_allocation(iostat, errmsg, this_module, "swap_columns_cdp")

        ! Swap columns.
        call copy(Qwrk, Q(j))
        call copy(Q(j), Q(i))
        call copy(Q(i), Qwrk)

        Rwrk(1) = Rii(j); Rii(j) = Rii(i); Rii(i) = real(Rwrk(1), kind=dp)
        iwrk = perm(j); perm(j) = perm(i) ; perm(i) = iwrk

        if (n > 0) then
            Rwrk = R(:n, j) ; R(:n, j) = R(:n, i) ; R(:n, i) = Rwrk
        endif
    end procedure swap_columns_cdp

end submodule qr_solvers
