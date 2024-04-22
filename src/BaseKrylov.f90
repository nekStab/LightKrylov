module lightkrylov_BaseKrylov
    use iso_fortran_env
    use stdlib_optval, only: optval
    use lightkrylov_constants
    use lightkrylov_AbstractVectors
    use lightkrylov_AbstractLinops
    implicit none
    
    private

    public :: qr
    public :: apply_permutation_matrix
    public :: apply_inverse_permutation_matrix
    public :: arnoldi
    public :: initialize_krylov_subspace
    public :: lanczos_bidiagonalization

    interface qr
        module procedure qr_no_pivoting_rsp
        module procedure qr_with_pivoting_rsp
        module procedure qr_no_pivoting_rdp
        module procedure qr_with_pivoting_rdp
        module procedure qr_no_pivoting_csp
        module procedure qr_with_pivoting_csp
        module procedure qr_no_pivoting_cdp
        module procedure qr_with_pivoting_cdp
    end interface

    interface swap_columns
        module procedure swap_columns_rsp
        module procedure swap_columns_rdp
        module procedure swap_columns_csp
        module procedure swap_columns_cdp
    end interface

    interface apply_permutation_matrix
        module procedure apply_permutation_matrix_rsp
        module procedure apply_permutation_matrix_rdp
        module procedure apply_permutation_matrix_csp
        module procedure apply_permutation_matrix_cdp
    end interface

    interface apply_inverse_permutation_matrix
        module procedure apply_inverse_permutation_matrix_rsp
        module procedure apply_inverse_permutation_matrix_rdp
        module procedure apply_inverse_permutation_matrix_csp
        module procedure apply_inverse_permutation_matrix_cdp
    end interface

    interface arnoldi
        module procedure arnoldi_rsp
        module procedure arnoldi_rdp
        module procedure arnoldi_csp
        module procedure arnoldi_cdp
    end interface

    interface initialize_krylov_subspace
        module procedure initialize_krylov_subspace_rsp
        module procedure initialize_krylov_subspace_rdp
        module procedure initialize_krylov_subspace_csp
        module procedure initialize_krylov_subspace_cdp
    end interface


    interface lanczos_bidiagonalization
        module procedure lanczos_bidiagonalization_rsp
        module procedure lanczos_bidiagonalization_rdp
        module procedure lanczos_bidiagonalization_csp
        module procedure lanczos_bidiagonalization_cdp
    end interface

contains

    !-------------------------------------
    !-----     UTILITY FUNCTIONS     -----
    !-------------------------------------

    subroutine initialize_krylov_subspace_rsp(X, X0)
        class(abstract_vector_rsp), intent(inout) :: X(:)
        class(abstract_vector_rsp), optional, intent(in) :: X0(:)

        ! Internal variables.
        class(abstract_vector_rsp), allocatable :: Q(:)
        real(sp), allocatable :: R(:, :)
        integer, allocatable :: perm(:)
        integer :: i, p, info

        ! Zero-out X.
        do i = 1, size(X)
            call X(i)%zero()
        enddo

        ! Deals with optional args.
        if(present(X0)) then
            p = size(X0)

            ! Initialize.
            allocate(Q(1:p), source=X0(1:p))
            do i = 1, p
                call Q(i)%zero()
                call Q(i)%add(X0(i))
            enddo
    
            ! Orthonormalize.
            allocate(R(1:p, 1:p)) ; R = 0.0_sp
            allocate(perm(1:p)) ; perm = 0
            call qr(Q, R, perm, info)
            do i = 1, p
                call X(i)%add(Q(i))
            enddo
        endif

        return
    end subroutine initialize_krylov_subspace_rsp

    subroutine initialize_krylov_subspace_rdp(X, X0)
        class(abstract_vector_rdp), intent(inout) :: X(:)
        class(abstract_vector_rdp), optional, intent(in) :: X0(:)

        ! Internal variables.
        class(abstract_vector_rdp), allocatable :: Q(:)
        real(dp), allocatable :: R(:, :)
        integer, allocatable :: perm(:)
        integer :: i, p, info

        ! Zero-out X.
        do i = 1, size(X)
            call X(i)%zero()
        enddo

        ! Deals with optional args.
        if(present(X0)) then
            p = size(X0)

            ! Initialize.
            allocate(Q(1:p), source=X0(1:p))
            do i = 1, p
                call Q(i)%zero()
                call Q(i)%add(X0(i))
            enddo
    
            ! Orthonormalize.
            allocate(R(1:p, 1:p)) ; R = 0.0_dp
            allocate(perm(1:p)) ; perm = 0
            call qr(Q, R, perm, info)
            do i = 1, p
                call X(i)%add(Q(i))
            enddo
        endif

        return
    end subroutine initialize_krylov_subspace_rdp

    subroutine initialize_krylov_subspace_csp(X, X0)
        class(abstract_vector_csp), intent(inout) :: X(:)
        class(abstract_vector_csp), optional, intent(in) :: X0(:)

        ! Internal variables.
        class(abstract_vector_csp), allocatable :: Q(:)
        complex(sp), allocatable :: R(:, :)
        integer, allocatable :: perm(:)
        integer :: i, p, info

        ! Zero-out X.
        do i = 1, size(X)
            call X(i)%zero()
        enddo

        ! Deals with optional args.
        if(present(X0)) then
            p = size(X0)

            ! Initialize.
            allocate(Q(1:p), source=X0(1:p))
            do i = 1, p
                call Q(i)%zero()
                call Q(i)%add(X0(i))
            enddo
    
            ! Orthonormalize.
            allocate(R(1:p, 1:p)) ; R = 0.0_sp
            allocate(perm(1:p)) ; perm = 0
            call qr(Q, R, perm, info)
            do i = 1, p
                call X(i)%add(Q(i))
            enddo
        endif

        return
    end subroutine initialize_krylov_subspace_csp

    subroutine initialize_krylov_subspace_cdp(X, X0)
        class(abstract_vector_cdp), intent(inout) :: X(:)
        class(abstract_vector_cdp), optional, intent(in) :: X0(:)

        ! Internal variables.
        class(abstract_vector_cdp), allocatable :: Q(:)
        complex(dp), allocatable :: R(:, :)
        integer, allocatable :: perm(:)
        integer :: i, p, info

        ! Zero-out X.
        do i = 1, size(X)
            call X(i)%zero()
        enddo

        ! Deals with optional args.
        if(present(X0)) then
            p = size(X0)

            ! Initialize.
            allocate(Q(1:p), source=X0(1:p))
            do i = 1, p
                call Q(i)%zero()
                call Q(i)%add(X0(i))
            enddo
    
            ! Orthonormalize.
            allocate(R(1:p, 1:p)) ; R = 0.0_dp
            allocate(perm(1:p)) ; perm = 0
            call qr(Q, R, perm, info)
            do i = 1, p
                call X(i)%add(Q(i))
            enddo
        endif

        return
    end subroutine initialize_krylov_subspace_cdp


    !------------------------------------
    !-----     QR FACTORIZATION     -----
    !------------------------------------

    subroutine qr_no_pivoting_rsp(Q, R, info, verbosity, tol)
        class(abstract_vector_rsp), intent(inout) :: Q(:)
        !! Array of `abstract_vector` to be orthonormalized.
        real(sp), intent(out) :: R(:, :)
        !! Upper triangular matrix \(\mathbf{R}\) resulting from the QR factorization.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        logical                       :: verbose
        real(sp), optional, intent(in) :: tol
        !! Tolerance to determine colinearity.
        real(sp)                       :: tolerance

        ! Internal variables.
        real(sp) :: beta
        integer :: idx, i, j, k, kdim, iwrk

        ! Deals with the optional args.
        verbose   = optval(verbosity, .true.)
        tolerance = optval(tol, rtol_sp)

        info = 0 ; R = 0.0_sp
        do j = 1, size(Q)
            ! First pass
            do i = 1, j-1
                beta = Q(i)%dot(Q(j))
                call Q(j)%axpby(1.0_sp, Q(i), -beta)
                R(i, j) = beta
            enddo

            ! Second pass
            do i = 1, j-1
                beta = Q(i)%dot(Q(j))
                call Q(j)%axpby(1.0_sp, Q(i), -beta)
                R(i, j) = R(i, j) + beta
            enddo

            ! Normalize column.
            beta = Q(j)%norm()

            ! Check for breakdown.
            if (abs(beta) < tolerance) then
                if(verbose) then
                    write(output_unit, *) "INFO: Colinear columns detected."
                    write(output_unit, *) "      (j, beta) = (", j, ", ", beta, ")"
                endif
                R(i, j) = 0.0_sp ; call Q(j)%zero()
            else
                call Q(j)%scal(1.0_sp / beta) ; R(j, j) = beta
            endif
        enddo

        return
    end subroutine qr_no_pivoting_rsp

    subroutine qr_with_pivoting_rsp(Q, R, perm, info, verbosity, tol)
        class(abstract_vector_rsp), intent(inout) :: Q(:)
        !! Array of `abstract_vector` to be orthonormalized.
        real(sp), intent(out) :: R(:, :)
        !! Upper triangular matrix resulting from the QR factorization.
        integer, intent(out) :: perm(size(Q))
        !! Permutation matrix.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        logical :: verbose
        real(sp), optional, intent(in) :: tol

        !! Tolerance to detect colinearity.
        real(sp) :: tolerance

        ! Internal variables
        real(sp) :: beta
        integer :: idx, i, j, k, kdim, iwrk
        integer :: idxv(1)
        real(sp)  :: Rii(size(Q))

        info = 0 ; kdim = size(Q)
        R = 0.0_sp ; Rii = 0.0_sp
        
        ! Deals with the optional arguments.
        verbose = optval(verbosity, .false.)
        tolerance = optval(tol, rtol_sp)

        ! Initialize diagonal entries.
        do i = 1, kdim
            perm(i) = i
            Rii(i) = Q(i)%dot(Q(i))
        enddo

        qr_step: do j = 1, kdim
            idxv = maxloc(abs(Rii)) ; idx = idxv(1)
            if (abs(Rii(idx)) < tolerance) then
                do i = j, kdim
                    call Q(i)%rand(.false.)
                    ! Orthogonalize against existing columns
                    do k = 1, i-1
                        beta = Q(i)%dot(Q(k))
                        call Q(i)%axpby(1.0_sp, Q(k), -beta)
                    enddo
                    beta = Q(i)%norm() ; call Q(i)%scal(1.0_sp / beta)
                enddo
                exit qr_step
            endif

            call swap_columns(Q, R, Rii, perm, j, idx)

            ! Normalize column.
            beta = Q(j)%norm()

            if (abs(beta) < tolerance) then
                R(j, j) = 0.0_sp ; call Q(j)%zero()
            else
                R(j, j) = beta ; call Q(j)%scal(1.0_sp / beta)
            endif

            ! Orthonormalize all columns against new vector (MGS)
            do i = j+1, kdim
                beta = Q(j)%dot(Q(i))
                call Q(i)%axpby(1.0_sp, Q(j), -beta)
                R(j, i) = beta
            enddo

            ! Update Rii
            Rii(j) = 0.0_sp
            do i = j+1, kdim
                Rii(i) = Rii(i) - R(j, i)**2
            enddo

        enddo qr_step

        return
    end subroutine qr_with_pivoting_rsp

    subroutine swap_columns_rsp(Q, R, Rii, perm, i, j)
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

        ! Internal variables.
        class(abstract_vector_rsp), allocatable :: Qwrk
        real(sp), allocatable :: Rwrk(:)
        integer :: iwrk, m, n

        ! Sanity checks.
        m = size(Q) ; n = min(i, j) - 1

        ! Allocations.
        allocate(Qwrk, source=Q(1)) ; call Qwrk%zero()
        allocate(Rwrk(1:max(1, n))) ; Rwrk = 0.0_sp

        ! Swap columns.
        call Qwrk%axpby(0.0_sp, Q(j), 1.0_sp)
        call Q(j)%axpby(0.0_sp, Q(i), 1.0_sp)
        call Q(i)%axpby(0.0_sp, Qwrk, 1.0_sp)
        
        Rwrk(1) = Rii(j); Rii(j) = Rii(i); Rii(i) = Rwrk(1)
        iwrk = perm(j); perm(j) = perm(i) ; perm(i) = iwrk

        if (n > 0) then
            Rwrk = R(1:n, j) ; R(1:n, j) = R(1:n, i) ; R(1:n, i) = Rwrk
        endif

        return
    end subroutine swap_columns_rsp

    subroutine apply_permutation_matrix_rsp(Q, perm)
        class(abstract_vector_rsp), intent(inout) :: Q(:)
        !! Basis vectors to be permuted.
        integer, intent(in) :: perm(:)
        !! Permutation matrix (vector representation).

        ! Internal variables.
        integer :: i
        class(abstract_vector_rsp), allocatable :: Qwrk(:)

        allocate(Qwrk, source=Q)
        do i = 1, size(perm)
            call Q(i)%zero() ; call Q(i)%add(Qwrk(perm(i)))
        enddo

        return
    end subroutine apply_permutation_matrix_rsp

    subroutine apply_inverse_permutation_matrix_rsp(Q, perm)
        class(abstract_vector_rsp), intent(inout) :: Q(:)
        !! Basis vectors to be (un-) permuted.
        integer, intent(in) :: perm(:)
        !! Permutation matrix (vector representation).

        ! Internal variables.
        integer :: i
        integer :: inv_perm(size(perm))
        class(abstract_vector_rsp), allocatable :: Qwrk(:)

        allocate(Qwrk, source=Q) ; inv_perm = 0

        ! Inverse permutation vector.
        do i = 1, size(perm)
            inv_perm(perm(i)) = i
        enddo

        ! Undo permutation.
        do i = 1, size(perm)
            call Q(i)%zero() ; call Q(i)%add(Qwrk(inv_perm(i)))
        enddo

        return
    end subroutine apply_inverse_permutation_matrix_rsp

    subroutine qr_no_pivoting_rdp(Q, R, info, verbosity, tol)
        class(abstract_vector_rdp), intent(inout) :: Q(:)
        !! Array of `abstract_vector` to be orthonormalized.
        real(dp), intent(out) :: R(:, :)
        !! Upper triangular matrix \(\mathbf{R}\) resulting from the QR factorization.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        logical                       :: verbose
        real(dp), optional, intent(in) :: tol
        !! Tolerance to determine colinearity.
        real(dp)                       :: tolerance

        ! Internal variables.
        real(dp) :: beta
        integer :: idx, i, j, k, kdim, iwrk

        ! Deals with the optional args.
        verbose   = optval(verbosity, .true.)
        tolerance = optval(tol, rtol_dp)

        info = 0 ; R = 0.0_dp
        do j = 1, size(Q)
            ! First pass
            do i = 1, j-1
                beta = Q(i)%dot(Q(j))
                call Q(j)%axpby(1.0_dp, Q(i), -beta)
                R(i, j) = beta
            enddo

            ! Second pass
            do i = 1, j-1
                beta = Q(i)%dot(Q(j))
                call Q(j)%axpby(1.0_dp, Q(i), -beta)
                R(i, j) = R(i, j) + beta
            enddo

            ! Normalize column.
            beta = Q(j)%norm()

            ! Check for breakdown.
            if (abs(beta) < tolerance) then
                if(verbose) then
                    write(output_unit, *) "INFO: Colinear columns detected."
                    write(output_unit, *) "      (j, beta) = (", j, ", ", beta, ")"
                endif
                R(i, j) = 0.0_dp ; call Q(j)%zero()
            else
                call Q(j)%scal(1.0_dp / beta) ; R(j, j) = beta
            endif
        enddo

        return
    end subroutine qr_no_pivoting_rdp

    subroutine qr_with_pivoting_rdp(Q, R, perm, info, verbosity, tol)
        class(abstract_vector_rdp), intent(inout) :: Q(:)
        !! Array of `abstract_vector` to be orthonormalized.
        real(dp), intent(out) :: R(:, :)
        !! Upper triangular matrix resulting from the QR factorization.
        integer, intent(out) :: perm(size(Q))
        !! Permutation matrix.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        logical :: verbose
        real(dp), optional, intent(in) :: tol

        !! Tolerance to detect colinearity.
        real(dp) :: tolerance

        ! Internal variables
        real(dp) :: beta
        integer :: idx, i, j, k, kdim, iwrk
        integer :: idxv(1)
        real(dp)  :: Rii(size(Q))

        info = 0 ; kdim = size(Q)
        R = 0.0_dp ; Rii = 0.0_dp
        
        ! Deals with the optional arguments.
        verbose = optval(verbosity, .false.)
        tolerance = optval(tol, rtol_dp)

        ! Initialize diagonal entries.
        do i = 1, kdim
            perm(i) = i
            Rii(i) = Q(i)%dot(Q(i))
        enddo

        qr_step: do j = 1, kdim
            idxv = maxloc(abs(Rii)) ; idx = idxv(1)
            if (abs(Rii(idx)) < tolerance) then
                do i = j, kdim
                    call Q(i)%rand(.false.)
                    ! Orthogonalize against existing columns
                    do k = 1, i-1
                        beta = Q(i)%dot(Q(k))
                        call Q(i)%axpby(1.0_dp, Q(k), -beta)
                    enddo
                    beta = Q(i)%norm() ; call Q(i)%scal(1.0_dp / beta)
                enddo
                exit qr_step
            endif

            call swap_columns(Q, R, Rii, perm, j, idx)

            ! Normalize column.
            beta = Q(j)%norm()

            if (abs(beta) < tolerance) then
                R(j, j) = 0.0_dp ; call Q(j)%zero()
            else
                R(j, j) = beta ; call Q(j)%scal(1.0_dp / beta)
            endif

            ! Orthonormalize all columns against new vector (MGS)
            do i = j+1, kdim
                beta = Q(j)%dot(Q(i))
                call Q(i)%axpby(1.0_dp, Q(j), -beta)
                R(j, i) = beta
            enddo

            ! Update Rii
            Rii(j) = 0.0_dp
            do i = j+1, kdim
                Rii(i) = Rii(i) - R(j, i)**2
            enddo

        enddo qr_step

        return
    end subroutine qr_with_pivoting_rdp

    subroutine swap_columns_rdp(Q, R, Rii, perm, i, j)
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

        ! Internal variables.
        class(abstract_vector_rdp), allocatable :: Qwrk
        real(dp), allocatable :: Rwrk(:)
        integer :: iwrk, m, n

        ! Sanity checks.
        m = size(Q) ; n = min(i, j) - 1

        ! Allocations.
        allocate(Qwrk, source=Q(1)) ; call Qwrk%zero()
        allocate(Rwrk(1:max(1, n))) ; Rwrk = 0.0_dp

        ! Swap columns.
        call Qwrk%axpby(0.0_dp, Q(j), 1.0_dp)
        call Q(j)%axpby(0.0_dp, Q(i), 1.0_dp)
        call Q(i)%axpby(0.0_dp, Qwrk, 1.0_dp)
        
        Rwrk(1) = Rii(j); Rii(j) = Rii(i); Rii(i) = Rwrk(1)
        iwrk = perm(j); perm(j) = perm(i) ; perm(i) = iwrk

        if (n > 0) then
            Rwrk = R(1:n, j) ; R(1:n, j) = R(1:n, i) ; R(1:n, i) = Rwrk
        endif

        return
    end subroutine swap_columns_rdp

    subroutine apply_permutation_matrix_rdp(Q, perm)
        class(abstract_vector_rdp), intent(inout) :: Q(:)
        !! Basis vectors to be permuted.
        integer, intent(in) :: perm(:)
        !! Permutation matrix (vector representation).

        ! Internal variables.
        integer :: i
        class(abstract_vector_rdp), allocatable :: Qwrk(:)

        allocate(Qwrk, source=Q)
        do i = 1, size(perm)
            call Q(i)%zero() ; call Q(i)%add(Qwrk(perm(i)))
        enddo

        return
    end subroutine apply_permutation_matrix_rdp

    subroutine apply_inverse_permutation_matrix_rdp(Q, perm)
        class(abstract_vector_rdp), intent(inout) :: Q(:)
        !! Basis vectors to be (un-) permuted.
        integer, intent(in) :: perm(:)
        !! Permutation matrix (vector representation).

        ! Internal variables.
        integer :: i
        integer :: inv_perm(size(perm))
        class(abstract_vector_rdp), allocatable :: Qwrk(:)

        allocate(Qwrk, source=Q) ; inv_perm = 0

        ! Inverse permutation vector.
        do i = 1, size(perm)
            inv_perm(perm(i)) = i
        enddo

        ! Undo permutation.
        do i = 1, size(perm)
            call Q(i)%zero() ; call Q(i)%add(Qwrk(inv_perm(i)))
        enddo

        return
    end subroutine apply_inverse_permutation_matrix_rdp

    subroutine qr_no_pivoting_csp(Q, R, info, verbosity, tol)
        class(abstract_vector_csp), intent(inout) :: Q(:)
        !! Array of `abstract_vector` to be orthonormalized.
        complex(sp), intent(out) :: R(:, :)
        !! Upper triangular matrix \(\mathbf{R}\) resulting from the QR factorization.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        logical                       :: verbose
        real(sp), optional, intent(in) :: tol
        !! Tolerance to determine colinearity.
        real(sp)                       :: tolerance

        ! Internal variables.
        complex(sp) :: beta
        integer :: idx, i, j, k, kdim, iwrk

        ! Deals with the optional args.
        verbose   = optval(verbosity, .true.)
        tolerance = optval(tol, rtol_sp)

        info = 0 ; R = 0.0_sp
        do j = 1, size(Q)
            ! First pass
            do i = 1, j-1
                beta = Q(i)%dot(Q(j))
                call Q(j)%axpby(cmplx(1.0_sp, 0.0_sp, kind=sp), Q(i), -beta)
                R(i, j) = beta
            enddo

            ! Second pass
            do i = 1, j-1
                beta = Q(i)%dot(Q(j))
                call Q(j)%axpby(cmplx(1.0_sp, 0.0_sp, kind=sp), Q(i), -beta)
                R(i, j) = R(i, j) + beta
            enddo

            ! Normalize column.
            beta = Q(j)%norm()

            ! Check for breakdown.
            if (abs(beta) < tolerance) then
                if(verbose) then
                    write(output_unit, *) "INFO: Colinear columns detected."
                    write(output_unit, *) "      (j, beta) = (", j, ", ", beta, ")"
                endif
                R(i, j) = 0.0_sp ; call Q(j)%zero()
            else
                call Q(j)%scal(1.0_sp / beta) ; R(j, j) = beta
            endif
        enddo

        return
    end subroutine qr_no_pivoting_csp

    subroutine qr_with_pivoting_csp(Q, R, perm, info, verbosity, tol)
        class(abstract_vector_csp), intent(inout) :: Q(:)
        !! Array of `abstract_vector` to be orthonormalized.
        complex(sp), intent(out) :: R(:, :)
        !! Upper triangular matrix resulting from the QR factorization.
        integer, intent(out) :: perm(size(Q))
        !! Permutation matrix.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        logical :: verbose
        real(sp), optional, intent(in) :: tol

        !! Tolerance to detect colinearity.
        real(sp) :: tolerance

        ! Internal variables
        complex(sp) :: beta
        integer :: idx, i, j, k, kdim, iwrk
        integer :: idxv(1)
        complex(sp)  :: Rii(size(Q))

        info = 0 ; kdim = size(Q)
        R = 0.0_sp ; Rii = 0.0_sp
        
        ! Deals with the optional arguments.
        verbose = optval(verbosity, .false.)
        tolerance = optval(tol, rtol_sp)

        ! Initialize diagonal entries.
        do i = 1, kdim
            perm(i) = i
            Rii(i) = Q(i)%dot(Q(i))
        enddo

        qr_step: do j = 1, kdim
            idxv = maxloc(abs(Rii)) ; idx = idxv(1)
            if (abs(Rii(idx)) < tolerance) then
                do i = j, kdim
                    call Q(i)%rand(.false.)
                    ! Orthogonalize against existing columns
                    do k = 1, i-1
                        beta = Q(i)%dot(Q(k))
                        call Q(i)%axpby(cmplx(1.0_sp, 0.0_sp, kind=sp), Q(k), -beta)
                    enddo
                    beta = Q(i)%norm() ; call Q(i)%scal(1.0_sp / beta)
                enddo
                exit qr_step
            endif

            call swap_columns(Q, R, Rii, perm, j, idx)

            ! Normalize column.
            beta = Q(j)%norm()

            if (abs(beta) < tolerance) then
                R(j, j) = 0.0_sp ; call Q(j)%zero()
            else
                R(j, j) = beta ; call Q(j)%scal(1.0_sp / beta)
            endif

            ! Orthonormalize all columns against new vector (MGS)
            do i = j+1, kdim
                beta = Q(j)%dot(Q(i))
                call Q(i)%axpby(cmplx(1.0_sp, 0.0_sp, kind=sp), Q(j), -beta)
                R(j, i) = beta
            enddo

            ! Update Rii
            Rii(j) = 0.0_sp
            do i = j+1, kdim
                Rii(i) = Rii(i) - R(j, i)**2
            enddo

        enddo qr_step

        return
    end subroutine qr_with_pivoting_csp

    subroutine swap_columns_csp(Q, R, Rii, perm, i, j)
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

        ! Internal variables.
        class(abstract_vector_csp), allocatable :: Qwrk
        complex(sp), allocatable :: Rwrk(:)
        integer :: iwrk, m, n

        ! Sanity checks.
        m = size(Q) ; n = min(i, j) - 1

        ! Allocations.
        allocate(Qwrk, source=Q(1)) ; call Qwrk%zero()
        allocate(Rwrk(1:max(1, n))) ; Rwrk = 0.0_sp

        ! Swap columns.
        call Qwrk%axpby(cmplx(0.0_sp, 0.0_sp, kind=sp), Q(j), cmplx(1.0_sp, 0.0_sp, kind=sp))
        call Q(j)%axpby(cmplx(0.0_sp, 0.0_sp, kind=sp), Q(i), cmplx(1.0_sp, 0.0_sp, kind=sp))
        call Q(i)%axpby(cmplx(0.0_sp, 0.0_sp, kind=sp), Qwrk, cmplx(1.0_sp, 0.0_sp, kind=sp))
        
        Rwrk(1) = Rii(j); Rii(j) = Rii(i); Rii(i) = Rwrk(1)
        iwrk = perm(j); perm(j) = perm(i) ; perm(i) = iwrk

        if (n > 0) then
            Rwrk = R(1:n, j) ; R(1:n, j) = R(1:n, i) ; R(1:n, i) = Rwrk
        endif

        return
    end subroutine swap_columns_csp

    subroutine apply_permutation_matrix_csp(Q, perm)
        class(abstract_vector_csp), intent(inout) :: Q(:)
        !! Basis vectors to be permuted.
        integer, intent(in) :: perm(:)
        !! Permutation matrix (vector representation).

        ! Internal variables.
        integer :: i
        class(abstract_vector_csp), allocatable :: Qwrk(:)

        allocate(Qwrk, source=Q)
        do i = 1, size(perm)
            call Q(i)%zero() ; call Q(i)%add(Qwrk(perm(i)))
        enddo

        return
    end subroutine apply_permutation_matrix_csp

    subroutine apply_inverse_permutation_matrix_csp(Q, perm)
        class(abstract_vector_csp), intent(inout) :: Q(:)
        !! Basis vectors to be (un-) permuted.
        integer, intent(in) :: perm(:)
        !! Permutation matrix (vector representation).

        ! Internal variables.
        integer :: i
        integer :: inv_perm(size(perm))
        class(abstract_vector_csp), allocatable :: Qwrk(:)

        allocate(Qwrk, source=Q) ; inv_perm = 0

        ! Inverse permutation vector.
        do i = 1, size(perm)
            inv_perm(perm(i)) = i
        enddo

        ! Undo permutation.
        do i = 1, size(perm)
            call Q(i)%zero() ; call Q(i)%add(Qwrk(inv_perm(i)))
        enddo

        return
    end subroutine apply_inverse_permutation_matrix_csp

    subroutine qr_no_pivoting_cdp(Q, R, info, verbosity, tol)
        class(abstract_vector_cdp), intent(inout) :: Q(:)
        !! Array of `abstract_vector` to be orthonormalized.
        complex(dp), intent(out) :: R(:, :)
        !! Upper triangular matrix \(\mathbf{R}\) resulting from the QR factorization.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        logical                       :: verbose
        real(dp), optional, intent(in) :: tol
        !! Tolerance to determine colinearity.
        real(dp)                       :: tolerance

        ! Internal variables.
        complex(dp) :: beta
        integer :: idx, i, j, k, kdim, iwrk

        ! Deals with the optional args.
        verbose   = optval(verbosity, .true.)
        tolerance = optval(tol, rtol_dp)

        info = 0 ; R = 0.0_dp
        do j = 1, size(Q)
            ! First pass
            do i = 1, j-1
                beta = Q(i)%dot(Q(j))
                call Q(j)%axpby(cmplx(1.0_dp, 0.0_dp, kind=dp), Q(i), -beta)
                R(i, j) = beta
            enddo

            ! Second pass
            do i = 1, j-1
                beta = Q(i)%dot(Q(j))
                call Q(j)%axpby(cmplx(1.0_dp, 0.0_dp, kind=dp), Q(i), -beta)
                R(i, j) = R(i, j) + beta
            enddo

            ! Normalize column.
            beta = Q(j)%norm()

            ! Check for breakdown.
            if (abs(beta) < tolerance) then
                if(verbose) then
                    write(output_unit, *) "INFO: Colinear columns detected."
                    write(output_unit, *) "      (j, beta) = (", j, ", ", beta, ")"
                endif
                R(i, j) = 0.0_dp ; call Q(j)%zero()
            else
                call Q(j)%scal(1.0_dp / beta) ; R(j, j) = beta
            endif
        enddo

        return
    end subroutine qr_no_pivoting_cdp

    subroutine qr_with_pivoting_cdp(Q, R, perm, info, verbosity, tol)
        class(abstract_vector_cdp), intent(inout) :: Q(:)
        !! Array of `abstract_vector` to be orthonormalized.
        complex(dp), intent(out) :: R(:, :)
        !! Upper triangular matrix resulting from the QR factorization.
        integer, intent(out) :: perm(size(Q))
        !! Permutation matrix.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        logical :: verbose
        real(dp), optional, intent(in) :: tol

        !! Tolerance to detect colinearity.
        real(dp) :: tolerance

        ! Internal variables
        complex(dp) :: beta
        integer :: idx, i, j, k, kdim, iwrk
        integer :: idxv(1)
        complex(dp)  :: Rii(size(Q))

        info = 0 ; kdim = size(Q)
        R = 0.0_dp ; Rii = 0.0_dp
        
        ! Deals with the optional arguments.
        verbose = optval(verbosity, .false.)
        tolerance = optval(tol, rtol_dp)

        ! Initialize diagonal entries.
        do i = 1, kdim
            perm(i) = i
            Rii(i) = Q(i)%dot(Q(i))
        enddo

        qr_step: do j = 1, kdim
            idxv = maxloc(abs(Rii)) ; idx = idxv(1)
            if (abs(Rii(idx)) < tolerance) then
                do i = j, kdim
                    call Q(i)%rand(.false.)
                    ! Orthogonalize against existing columns
                    do k = 1, i-1
                        beta = Q(i)%dot(Q(k))
                        call Q(i)%axpby(cmplx(1.0_dp, 0.0_dp, kind=dp), Q(k), -beta)
                    enddo
                    beta = Q(i)%norm() ; call Q(i)%scal(1.0_dp / beta)
                enddo
                exit qr_step
            endif

            call swap_columns(Q, R, Rii, perm, j, idx)

            ! Normalize column.
            beta = Q(j)%norm()

            if (abs(beta) < tolerance) then
                R(j, j) = 0.0_dp ; call Q(j)%zero()
            else
                R(j, j) = beta ; call Q(j)%scal(1.0_dp / beta)
            endif

            ! Orthonormalize all columns against new vector (MGS)
            do i = j+1, kdim
                beta = Q(j)%dot(Q(i))
                call Q(i)%axpby(cmplx(1.0_dp, 0.0_dp, kind=dp), Q(j), -beta)
                R(j, i) = beta
            enddo

            ! Update Rii
            Rii(j) = 0.0_dp
            do i = j+1, kdim
                Rii(i) = Rii(i) - R(j, i)**2
            enddo

        enddo qr_step

        return
    end subroutine qr_with_pivoting_cdp

    subroutine swap_columns_cdp(Q, R, Rii, perm, i, j)
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

        ! Internal variables.
        class(abstract_vector_cdp), allocatable :: Qwrk
        complex(dp), allocatable :: Rwrk(:)
        integer :: iwrk, m, n

        ! Sanity checks.
        m = size(Q) ; n = min(i, j) - 1

        ! Allocations.
        allocate(Qwrk, source=Q(1)) ; call Qwrk%zero()
        allocate(Rwrk(1:max(1, n))) ; Rwrk = 0.0_dp

        ! Swap columns.
        call Qwrk%axpby(cmplx(0.0_dp, 0.0_dp, kind=dp), Q(j), cmplx(1.0_dp, 0.0_dp, kind=dp))
        call Q(j)%axpby(cmplx(0.0_dp, 0.0_dp, kind=dp), Q(i), cmplx(1.0_dp, 0.0_dp, kind=dp))
        call Q(i)%axpby(cmplx(0.0_dp, 0.0_dp, kind=dp), Qwrk, cmplx(1.0_dp, 0.0_dp, kind=dp))
        
        Rwrk(1) = Rii(j); Rii(j) = Rii(i); Rii(i) = Rwrk(1)
        iwrk = perm(j); perm(j) = perm(i) ; perm(i) = iwrk

        if (n > 0) then
            Rwrk = R(1:n, j) ; R(1:n, j) = R(1:n, i) ; R(1:n, i) = Rwrk
        endif

        return
    end subroutine swap_columns_cdp

    subroutine apply_permutation_matrix_cdp(Q, perm)
        class(abstract_vector_cdp), intent(inout) :: Q(:)
        !! Basis vectors to be permuted.
        integer, intent(in) :: perm(:)
        !! Permutation matrix (vector representation).

        ! Internal variables.
        integer :: i
        class(abstract_vector_cdp), allocatable :: Qwrk(:)

        allocate(Qwrk, source=Q)
        do i = 1, size(perm)
            call Q(i)%zero() ; call Q(i)%add(Qwrk(perm(i)))
        enddo

        return
    end subroutine apply_permutation_matrix_cdp

    subroutine apply_inverse_permutation_matrix_cdp(Q, perm)
        class(abstract_vector_cdp), intent(inout) :: Q(:)
        !! Basis vectors to be (un-) permuted.
        integer, intent(in) :: perm(:)
        !! Permutation matrix (vector representation).

        ! Internal variables.
        integer :: i
        integer :: inv_perm(size(perm))
        class(abstract_vector_cdp), allocatable :: Qwrk(:)

        allocate(Qwrk, source=Q) ; inv_perm = 0

        ! Inverse permutation vector.
        do i = 1, size(perm)
            inv_perm(perm(i)) = i
        enddo

        ! Undo permutation.
        do i = 1, size(perm)
            call Q(i)%zero() ; call Q(i)%add(Qwrk(inv_perm(i)))
        enddo

        return
    end subroutine apply_inverse_permutation_matrix_cdp


    !-----------------------------------------
    !-----     ARNOLDI FACTORIZATION     -----
    !-----------------------------------------

    subroutine arnoldi_rsp(A, X, H, info, kstart, kend, verbosity, tol, transpose, blksize)
        class(abstract_linop_rsp), intent(in) :: A
        !! Linear operator to be factorized.
        class(abstract_vector_rsp), intent(inout) :: X(:)
        !! Orthogonal basis for the generated Krylov subspace.
        real(sp), intent(inout) :: H(:, :)
        !! Upper Hessenberg matrix.
        integer, intent(out) :: info
        !! Information flag.
        integer, optional, intent(in) :: kstart
        !! Starting index for the Arnoldi factorization (default 1).
        integer, optional, intent(in) :: kend
        !! Final index for the Arnoldi factorization (default `size(X)-1`)
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        logical, optional, intent(in) :: transpose
        !! Whether \( \mathbf{A} \) is being transposed or not (default `.false.`)
        real(sp), optional, intent(in) :: tol
        !! Tolerance to determine whether an invariant subspace has been computed or not.
        integer, optional, intent(in) :: blksize
        !! Block size for block Arnoldi (default 1).

        ! Internal variables.
        integer :: k_start, k_end, p
        logical :: verbose, trans
        real(sp) :: tolerance
        real(sp) :: beta
        real(sp), allocatable :: res(:)
        integer, allocatable :: perm(:)
        integer :: k, i, kdim, kpm, kp, kpp

        ! Deals with optional non-unity blksize and allocations.
        p = optval(blksize, 1) ; allocate(res(1:p)) ; res = 0.0_sp
        allocate(perm(1:size(H, 2))) ; perm = 0 ; info = 0

        ! Check dimensions.
        kdim = (size(X) - p) / p

        ! Deal with the other optional args.
        k_start = optval(kstart, 1) ; k_end = optval(kend, kdim)
        verbose   = optval(verbosity, .false.)
        tolerance = optval (tol, atol_sp)
        trans     = optval(transpose, .false.)

        ! Arnoldi factorization.
        blk_arnoldi: do k = k_start, k_end
            ! Counters
            kpm = (k - 1) * p ; kp = kpm + p ; kpp = kp + p

            ! Matrix-vector product.
            if (trans) then
                do i = 1, p
                    call A%rmatvec(X(kpm+i), X(kp+i))
                enddo
            else
                do i = 1, p
                    call A%matvec(X(kpm+i), X(kp+i))
                enddo
            endif

            ! Update Hessenberg matrix and orthogonalize w.r.t. previous vectors.
            call update_hessenberg_matrix_rsp(H, X, k, p)

            ! Orthogonalize current blk vectors.
            call qr(X(kp+1:kpp), H(kp+1:kpp, kpm+1:kp), info)

            ! Extract residual norm (smallest diagonal element of H matrix).
            res = 0.0_sp
            do i = 1, p
                res(i) = H(kp+i, kpm+i)
            enddo
            beta = minval(abs(res))
            info = kp
            ! Exit Arnoldi loop if needed.
            if (beta < tolerance) then
                ! Dimension of the computed invariant subspace.
                info = kp
                ! Exit the Arnoldi iteration.
                exit blk_arnoldi
            endif

        enddo blk_arnoldi

        return
    end subroutine arnoldi_rsp

    subroutine update_hessenberg_matrix_rsp(H, X, k, blksize)
        integer, intent(in) :: k
        real(sp), intent(inout) :: H(:, :)
        class(abstract_vector_rsp), intent(inout) :: X(:)
        integer, optional, intent(in) :: blksize

        ! Internal variables.
        real(sp), allocatable :: wrk(:, :)
        integer :: p, kpm, kp, kpp, i, j

        ! Deals with optional non-unity block size.
        p = optval(blksize, 1)

        ! Counters and allocations.
        kpm = (k - 1) * p ; kp = kpm + p ; kpp = kp + p
        allocate(wrk(1:kp, 1:p)) ; wrk = 0.0_sp

        ! Orthogonalize residual vector w.r.t. to previous computed Krylov vectors.
        call innerprod_matrix(H(1:kp, kpm+1:kp), X(1:kp), X(kp+1:kpp))
        do i = 1, p
            do j = 1, kp
                call X(kp+i)%axpby(1.0_sp, X(j), -H(j, kpm+i))
            enddo
        enddo

        call innerprod_matrix(wrk, X(1:kp), X(kp+1:kpp))
        do i = 1, p
            do j = 1, kp
                call X(kp+i)%axpby(1.0_sp, X(j), -wrk(j, i))
            enddo
        enddo

        ! Update Hessenberg matrix with the second pass correction.
        H(1:kp, kpm+1:kp) = H(1:kp, kpm+1:kp) + wrk

        return
    end subroutine update_hessenberg_matrix_rsp

    subroutine arnoldi_rdp(A, X, H, info, kstart, kend, verbosity, tol, transpose, blksize)
        class(abstract_linop_rdp), intent(in) :: A
        !! Linear operator to be factorized.
        class(abstract_vector_rdp), intent(inout) :: X(:)
        !! Orthogonal basis for the generated Krylov subspace.
        real(dp), intent(inout) :: H(:, :)
        !! Upper Hessenberg matrix.
        integer, intent(out) :: info
        !! Information flag.
        integer, optional, intent(in) :: kstart
        !! Starting index for the Arnoldi factorization (default 1).
        integer, optional, intent(in) :: kend
        !! Final index for the Arnoldi factorization (default `size(X)-1`)
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        logical, optional, intent(in) :: transpose
        !! Whether \( \mathbf{A} \) is being transposed or not (default `.false.`)
        real(dp), optional, intent(in) :: tol
        !! Tolerance to determine whether an invariant subspace has been computed or not.
        integer, optional, intent(in) :: blksize
        !! Block size for block Arnoldi (default 1).

        ! Internal variables.
        integer :: k_start, k_end, p
        logical :: verbose, trans
        real(dp) :: tolerance
        real(dp) :: beta
        real(dp), allocatable :: res(:)
        integer, allocatable :: perm(:)
        integer :: k, i, kdim, kpm, kp, kpp

        ! Deals with optional non-unity blksize and allocations.
        p = optval(blksize, 1) ; allocate(res(1:p)) ; res = 0.0_dp
        allocate(perm(1:size(H, 2))) ; perm = 0 ; info = 0

        ! Check dimensions.
        kdim = (size(X) - p) / p

        ! Deal with the other optional args.
        k_start = optval(kstart, 1) ; k_end = optval(kend, kdim)
        verbose   = optval(verbosity, .false.)
        tolerance = optval (tol, atol_dp)
        trans     = optval(transpose, .false.)

        ! Arnoldi factorization.
        blk_arnoldi: do k = k_start, k_end
            ! Counters
            kpm = (k - 1) * p ; kp = kpm + p ; kpp = kp + p

            ! Matrix-vector product.
            if (trans) then
                do i = 1, p
                    call A%rmatvec(X(kpm+i), X(kp+i))
                enddo
            else
                do i = 1, p
                    call A%matvec(X(kpm+i), X(kp+i))
                enddo
            endif

            ! Update Hessenberg matrix and orthogonalize w.r.t. previous vectors.
            call update_hessenberg_matrix_rdp(H, X, k, p)

            ! Orthogonalize current blk vectors.
            call qr(X(kp+1:kpp), H(kp+1:kpp, kpm+1:kp), info)

            ! Extract residual norm (smallest diagonal element of H matrix).
            res = 0.0_dp
            do i = 1, p
                res(i) = H(kp+i, kpm+i)
            enddo
            beta = minval(abs(res))
            info = kp
            ! Exit Arnoldi loop if needed.
            if (beta < tolerance) then
                ! Dimension of the computed invariant subspace.
                info = kp
                ! Exit the Arnoldi iteration.
                exit blk_arnoldi
            endif

        enddo blk_arnoldi

        return
    end subroutine arnoldi_rdp

    subroutine update_hessenberg_matrix_rdp(H, X, k, blksize)
        integer, intent(in) :: k
        real(dp), intent(inout) :: H(:, :)
        class(abstract_vector_rdp), intent(inout) :: X(:)
        integer, optional, intent(in) :: blksize

        ! Internal variables.
        real(dp), allocatable :: wrk(:, :)
        integer :: p, kpm, kp, kpp, i, j

        ! Deals with optional non-unity block size.
        p = optval(blksize, 1)

        ! Counters and allocations.
        kpm = (k - 1) * p ; kp = kpm + p ; kpp = kp + p
        allocate(wrk(1:kp, 1:p)) ; wrk = 0.0_dp

        ! Orthogonalize residual vector w.r.t. to previous computed Krylov vectors.
        call innerprod_matrix(H(1:kp, kpm+1:kp), X(1:kp), X(kp+1:kpp))
        do i = 1, p
            do j = 1, kp
                call X(kp+i)%axpby(1.0_dp, X(j), -H(j, kpm+i))
            enddo
        enddo

        call innerprod_matrix(wrk, X(1:kp), X(kp+1:kpp))
        do i = 1, p
            do j = 1, kp
                call X(kp+i)%axpby(1.0_dp, X(j), -wrk(j, i))
            enddo
        enddo

        ! Update Hessenberg matrix with the second pass correction.
        H(1:kp, kpm+1:kp) = H(1:kp, kpm+1:kp) + wrk

        return
    end subroutine update_hessenberg_matrix_rdp

    subroutine arnoldi_csp(A, X, H, info, kstart, kend, verbosity, tol, transpose, blksize)
        class(abstract_linop_csp), intent(in) :: A
        !! Linear operator to be factorized.
        class(abstract_vector_csp), intent(inout) :: X(:)
        !! Orthogonal basis for the generated Krylov subspace.
        complex(sp), intent(inout) :: H(:, :)
        !! Upper Hessenberg matrix.
        integer, intent(out) :: info
        !! Information flag.
        integer, optional, intent(in) :: kstart
        !! Starting index for the Arnoldi factorization (default 1).
        integer, optional, intent(in) :: kend
        !! Final index for the Arnoldi factorization (default `size(X)-1`)
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        logical, optional, intent(in) :: transpose
        !! Whether \( \mathbf{A} \) is being transposed or not (default `.false.`)
        real(sp), optional, intent(in) :: tol
        !! Tolerance to determine whether an invariant subspace has been computed or not.
        integer, optional, intent(in) :: blksize
        !! Block size for block Arnoldi (default 1).

        ! Internal variables.
        integer :: k_start, k_end, p
        logical :: verbose, trans
        real(sp) :: tolerance
        real(sp) :: beta
        complex(sp), allocatable :: res(:)
        integer, allocatable :: perm(:)
        integer :: k, i, kdim, kpm, kp, kpp

        ! Deals with optional non-unity blksize and allocations.
        p = optval(blksize, 1) ; allocate(res(1:p)) ; res = 0.0_sp
        allocate(perm(1:size(H, 2))) ; perm = 0 ; info = 0

        ! Check dimensions.
        kdim = (size(X) - p) / p

        ! Deal with the other optional args.
        k_start = optval(kstart, 1) ; k_end = optval(kend, kdim)
        verbose   = optval(verbosity, .false.)
        tolerance = optval (tol, atol_sp)
        trans     = optval(transpose, .false.)

        ! Arnoldi factorization.
        blk_arnoldi: do k = k_start, k_end
            ! Counters
            kpm = (k - 1) * p ; kp = kpm + p ; kpp = kp + p

            ! Matrix-vector product.
            if (trans) then
                do i = 1, p
                    call A%rmatvec(X(kpm+i), X(kp+i))
                enddo
            else
                do i = 1, p
                    call A%matvec(X(kpm+i), X(kp+i))
                enddo
            endif

            ! Update Hessenberg matrix and orthogonalize w.r.t. previous vectors.
            call update_hessenberg_matrix_csp(H, X, k, p)

            ! Orthogonalize current blk vectors.
            call qr(X(kp+1:kpp), H(kp+1:kpp, kpm+1:kp), info)

            ! Extract residual norm (smallest diagonal element of H matrix).
            res = 0.0_sp
            do i = 1, p
                res(i) = H(kp+i, kpm+i)
            enddo
            beta = minval(abs(res))
            info = kp
            ! Exit Arnoldi loop if needed.
            if (beta < tolerance) then
                ! Dimension of the computed invariant subspace.
                info = kp
                ! Exit the Arnoldi iteration.
                exit blk_arnoldi
            endif

        enddo blk_arnoldi

        return
    end subroutine arnoldi_csp

    subroutine update_hessenberg_matrix_csp(H, X, k, blksize)
        integer, intent(in) :: k
        complex(sp), intent(inout) :: H(:, :)
        class(abstract_vector_csp), intent(inout) :: X(:)
        integer, optional, intent(in) :: blksize

        ! Internal variables.
        complex(sp), allocatable :: wrk(:, :)
        integer :: p, kpm, kp, kpp, i, j

        ! Deals with optional non-unity block size.
        p = optval(blksize, 1)

        ! Counters and allocations.
        kpm = (k - 1) * p ; kp = kpm + p ; kpp = kp + p
        allocate(wrk(1:kp, 1:p)) ; wrk = 0.0_sp

        ! Orthogonalize residual vector w.r.t. to previous computed Krylov vectors.
        call innerprod_matrix(H(1:kp, kpm+1:kp), X(1:kp), X(kp+1:kpp))
        do i = 1, p
            do j = 1, kp
                call X(kp+i)%axpby(cmplx(1.0_sp, 0.0_sp, kind=sp), X(j), -H(j, kpm+i))
            enddo
        enddo

        call innerprod_matrix(wrk, X(1:kp), X(kp+1:kpp))
        do i = 1, p
            do j = 1, kp
                call X(kp+i)%axpby(cmplx(1.0_sp, 0.0_sp, kind=sp), X(j), -wrk(j, i))
            enddo
        enddo

        ! Update Hessenberg matrix with the second pass correction.
        H(1:kp, kpm+1:kp) = H(1:kp, kpm+1:kp) + wrk

        return
    end subroutine update_hessenberg_matrix_csp

    subroutine arnoldi_cdp(A, X, H, info, kstart, kend, verbosity, tol, transpose, blksize)
        class(abstract_linop_cdp), intent(in) :: A
        !! Linear operator to be factorized.
        class(abstract_vector_cdp), intent(inout) :: X(:)
        !! Orthogonal basis for the generated Krylov subspace.
        complex(dp), intent(inout) :: H(:, :)
        !! Upper Hessenberg matrix.
        integer, intent(out) :: info
        !! Information flag.
        integer, optional, intent(in) :: kstart
        !! Starting index for the Arnoldi factorization (default 1).
        integer, optional, intent(in) :: kend
        !! Final index for the Arnoldi factorization (default `size(X)-1`)
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        logical, optional, intent(in) :: transpose
        !! Whether \( \mathbf{A} \) is being transposed or not (default `.false.`)
        real(dp), optional, intent(in) :: tol
        !! Tolerance to determine whether an invariant subspace has been computed or not.
        integer, optional, intent(in) :: blksize
        !! Block size for block Arnoldi (default 1).

        ! Internal variables.
        integer :: k_start, k_end, p
        logical :: verbose, trans
        real(dp) :: tolerance
        real(dp) :: beta
        complex(dp), allocatable :: res(:)
        integer, allocatable :: perm(:)
        integer :: k, i, kdim, kpm, kp, kpp

        ! Deals with optional non-unity blksize and allocations.
        p = optval(blksize, 1) ; allocate(res(1:p)) ; res = 0.0_dp
        allocate(perm(1:size(H, 2))) ; perm = 0 ; info = 0

        ! Check dimensions.
        kdim = (size(X) - p) / p

        ! Deal with the other optional args.
        k_start = optval(kstart, 1) ; k_end = optval(kend, kdim)
        verbose   = optval(verbosity, .false.)
        tolerance = optval (tol, atol_dp)
        trans     = optval(transpose, .false.)

        ! Arnoldi factorization.
        blk_arnoldi: do k = k_start, k_end
            ! Counters
            kpm = (k - 1) * p ; kp = kpm + p ; kpp = kp + p

            ! Matrix-vector product.
            if (trans) then
                do i = 1, p
                    call A%rmatvec(X(kpm+i), X(kp+i))
                enddo
            else
                do i = 1, p
                    call A%matvec(X(kpm+i), X(kp+i))
                enddo
            endif

            ! Update Hessenberg matrix and orthogonalize w.r.t. previous vectors.
            call update_hessenberg_matrix_cdp(H, X, k, p)

            ! Orthogonalize current blk vectors.
            call qr(X(kp+1:kpp), H(kp+1:kpp, kpm+1:kp), info)

            ! Extract residual norm (smallest diagonal element of H matrix).
            res = 0.0_dp
            do i = 1, p
                res(i) = H(kp+i, kpm+i)
            enddo
            beta = minval(abs(res))
            info = kp
            ! Exit Arnoldi loop if needed.
            if (beta < tolerance) then
                ! Dimension of the computed invariant subspace.
                info = kp
                ! Exit the Arnoldi iteration.
                exit blk_arnoldi
            endif

        enddo blk_arnoldi

        return
    end subroutine arnoldi_cdp

    subroutine update_hessenberg_matrix_cdp(H, X, k, blksize)
        integer, intent(in) :: k
        complex(dp), intent(inout) :: H(:, :)
        class(abstract_vector_cdp), intent(inout) :: X(:)
        integer, optional, intent(in) :: blksize

        ! Internal variables.
        complex(dp), allocatable :: wrk(:, :)
        integer :: p, kpm, kp, kpp, i, j

        ! Deals with optional non-unity block size.
        p = optval(blksize, 1)

        ! Counters and allocations.
        kpm = (k - 1) * p ; kp = kpm + p ; kpp = kp + p
        allocate(wrk(1:kp, 1:p)) ; wrk = 0.0_dp

        ! Orthogonalize residual vector w.r.t. to previous computed Krylov vectors.
        call innerprod_matrix(H(1:kp, kpm+1:kp), X(1:kp), X(kp+1:kpp))
        do i = 1, p
            do j = 1, kp
                call X(kp+i)%axpby(cmplx(1.0_dp, 0.0_dp, kind=dp), X(j), -H(j, kpm+i))
            enddo
        enddo

        call innerprod_matrix(wrk, X(1:kp), X(kp+1:kpp))
        do i = 1, p
            do j = 1, kp
                call X(kp+i)%axpby(cmplx(1.0_dp, 0.0_dp, kind=dp), X(j), -wrk(j, i))
            enddo
        enddo

        ! Update Hessenberg matrix with the second pass correction.
        H(1:kp, kpm+1:kp) = H(1:kp, kpm+1:kp) + wrk

        return
    end subroutine update_hessenberg_matrix_cdp



    !---------------------------------------------
    !-----     LANCZOS BIDIAGONALIZATION     -----
    !---------------------------------------------

    subroutine lanczos_bidiagonalization_rsp(A, U, V, B, info, kstart, kend, verbosity, tol)
        class(abstract_linop_rsp), intent(in) :: A
        !! Linear operator to be factorized.
        class(abstract_vector_rsp), intent(inout) :: U(:)
        !! Orthonormal basis for the column span of \(\mathbf{A}\). On entry, `U(1)` needs to be set to
        !! the starting Krylov vector.
        class(abstract_vector_rsp), intent(inout) :: V(:)
        !! Orthonormal basis for the row span of \(\mathbf{A}\).
        real(sp), intent(inout) :: B(:, :)
        !! Bidiagonal matrix.
        integer, intent(out) :: info
        !! Information flag.
        integer, optional, intent(in) :: kstart
        !! Starting index for the Lanczos factorization (default 1).
        integer, optional, intent(in) :: kend
        !! Final index for the Lanczos factorization (default 1).
        logical, optional, intent(in) :: verbosity
        !! Verbosity control (default `.false.`)
        real(sp), optional, intent(in) :: tol
        !! Tolerance to determine whether invariant subspaces have been computed or not.

        ! Internal variables.
        integer :: k_start, k_end
        logical :: verbose
        real(sp) :: tolerance
        real(sp) :: alpha, beta, gamma
        integer :: i, j, k, kdim


        info = 0

        ! Krylov subspace dimension.
        kdim = size(U) - 1

        ! Deals with the optional args.
        k_start = optval(kstart, 1)
        k_end   = optval(kend, kdim)
        verbose = optval(verbosity, .false.)
        tolerance = optval(tol, atol_sp)

        ! Lanczos bidiagonalization.
        lanczos : do k = k_start, k_end
            ! Transpose matrix-vector product.
            call A%rmatvec(U(k), V(k))

            ! Full reorthogonalization of the right Krylov subspace.
            do j = 1, k-1
                gamma = V(j)%dot(V(k))
                call V(k)%axpby(1.0_sp, V(j), -gamma)
            enddo

            ! Normalization step.
            alpha = V(k)%norm() ; B(k, k) = alpha
            if (abs(alpha) > tolerance) then
                call V(k)%scal(1.0_sp/alpha)
            else
                info = k
                exit lanczos
            endif

            ! Matrix-vector product.
            call A%matvec(V(k), U(k+1))

            ! Full re-orthogonalization of the left Krylov subspace.
            do j = 1, k
                gamma = U(j)%dot(U(k+1))
                call U(k+1)%axpby(1.0_sp, U(j), -gamma)
           enddo

            ! Normalization step
            beta = U(k+1)%norm() ; B(k+1, k) = beta
            if (abs(beta) > tolerance) then
                call U(k+1)%scal(1.0_sp / beta)
            else
                info = k
                exit lanczos
            endif

        enddo lanczos

        return
    end subroutine lanczos_bidiagonalization_rsp

    subroutine lanczos_bidiagonalization_rdp(A, U, V, B, info, kstart, kend, verbosity, tol)
        class(abstract_linop_rdp), intent(in) :: A
        !! Linear operator to be factorized.
        class(abstract_vector_rdp), intent(inout) :: U(:)
        !! Orthonormal basis for the column span of \(\mathbf{A}\). On entry, `U(1)` needs to be set to
        !! the starting Krylov vector.
        class(abstract_vector_rdp), intent(inout) :: V(:)
        !! Orthonormal basis for the row span of \(\mathbf{A}\).
        real(dp), intent(inout) :: B(:, :)
        !! Bidiagonal matrix.
        integer, intent(out) :: info
        !! Information flag.
        integer, optional, intent(in) :: kstart
        !! Starting index for the Lanczos factorization (default 1).
        integer, optional, intent(in) :: kend
        !! Final index for the Lanczos factorization (default 1).
        logical, optional, intent(in) :: verbosity
        !! Verbosity control (default `.false.`)
        real(dp), optional, intent(in) :: tol
        !! Tolerance to determine whether invariant subspaces have been computed or not.

        ! Internal variables.
        integer :: k_start, k_end
        logical :: verbose
        real(dp) :: tolerance
        real(dp) :: alpha, beta, gamma
        integer :: i, j, k, kdim


        info = 0

        ! Krylov subspace dimension.
        kdim = size(U) - 1

        ! Deals with the optional args.
        k_start = optval(kstart, 1)
        k_end   = optval(kend, kdim)
        verbose = optval(verbosity, .false.)
        tolerance = optval(tol, atol_dp)

        ! Lanczos bidiagonalization.
        lanczos : do k = k_start, k_end
            ! Transpose matrix-vector product.
            call A%rmatvec(U(k), V(k))

            ! Full reorthogonalization of the right Krylov subspace.
            do j = 1, k-1
                gamma = V(j)%dot(V(k))
                call V(k)%axpby(1.0_dp, V(j), -gamma)
            enddo

            ! Normalization step.
            alpha = V(k)%norm() ; B(k, k) = alpha
            if (abs(alpha) > tolerance) then
                call V(k)%scal(1.0_dp/alpha)
            else
                info = k
                exit lanczos
            endif

            ! Matrix-vector product.
            call A%matvec(V(k), U(k+1))

            ! Full re-orthogonalization of the left Krylov subspace.
            do j = 1, k
                gamma = U(j)%dot(U(k+1))
                call U(k+1)%axpby(1.0_dp, U(j), -gamma)
           enddo

            ! Normalization step
            beta = U(k+1)%norm() ; B(k+1, k) = beta
            if (abs(beta) > tolerance) then
                call U(k+1)%scal(1.0_dp / beta)
            else
                info = k
                exit lanczos
            endif

        enddo lanczos

        return
    end subroutine lanczos_bidiagonalization_rdp

    subroutine lanczos_bidiagonalization_csp(A, U, V, B, info, kstart, kend, verbosity, tol)
        class(abstract_linop_csp), intent(in) :: A
        !! Linear operator to be factorized.
        class(abstract_vector_csp), intent(inout) :: U(:)
        !! Orthonormal basis for the column span of \(\mathbf{A}\). On entry, `U(1)` needs to be set to
        !! the starting Krylov vector.
        class(abstract_vector_csp), intent(inout) :: V(:)
        !! Orthonormal basis for the row span of \(\mathbf{A}\).
        complex(sp), intent(inout) :: B(:, :)
        !! Bidiagonal matrix.
        integer, intent(out) :: info
        !! Information flag.
        integer, optional, intent(in) :: kstart
        !! Starting index for the Lanczos factorization (default 1).
        integer, optional, intent(in) :: kend
        !! Final index for the Lanczos factorization (default 1).
        logical, optional, intent(in) :: verbosity
        !! Verbosity control (default `.false.`)
        real(sp), optional, intent(in) :: tol
        !! Tolerance to determine whether invariant subspaces have been computed or not.

        ! Internal variables.
        integer :: k_start, k_end
        logical :: verbose
        real(sp) :: tolerance
        complex(sp) :: alpha, beta, gamma
        integer :: i, j, k, kdim


        info = 0

        ! Krylov subspace dimension.
        kdim = size(U) - 1

        ! Deals with the optional args.
        k_start = optval(kstart, 1)
        k_end   = optval(kend, kdim)
        verbose = optval(verbosity, .false.)
        tolerance = optval(tol, atol_sp)

        ! Lanczos bidiagonalization.
        lanczos : do k = k_start, k_end
            ! Transpose matrix-vector product.
            call A%rmatvec(U(k), V(k))

            ! Full reorthogonalization of the right Krylov subspace.
            do j = 1, k-1
                gamma = V(j)%dot(V(k))
                call V(k)%axpby(cmplx(1.0_sp, 0.0_sp, kind=sp), V(j), -gamma)
            enddo

            ! Normalization step.
            alpha = V(k)%norm() ; B(k, k) = alpha
            if (abs(alpha) > tolerance) then
                call V(k)%scal(1.0_sp/alpha)
            else
                info = k
                exit lanczos
            endif

            ! Matrix-vector product.
            call A%matvec(V(k), U(k+1))

            ! Full re-orthogonalization of the left Krylov subspace.
            do j = 1, k
                gamma = U(j)%dot(U(k+1))
                call U(k+1)%axpby(cmplx(1.0_sp, 0.0_sp, kind=sp), U(j), -gamma)
           enddo

            ! Normalization step
            beta = U(k+1)%norm() ; B(k+1, k) = beta
            if (abs(beta) > tolerance) then
                call U(k+1)%scal(1.0_sp / beta)
            else
                info = k
                exit lanczos
            endif

        enddo lanczos

        return
    end subroutine lanczos_bidiagonalization_csp

    subroutine lanczos_bidiagonalization_cdp(A, U, V, B, info, kstart, kend, verbosity, tol)
        class(abstract_linop_cdp), intent(in) :: A
        !! Linear operator to be factorized.
        class(abstract_vector_cdp), intent(inout) :: U(:)
        !! Orthonormal basis for the column span of \(\mathbf{A}\). On entry, `U(1)` needs to be set to
        !! the starting Krylov vector.
        class(abstract_vector_cdp), intent(inout) :: V(:)
        !! Orthonormal basis for the row span of \(\mathbf{A}\).
        complex(dp), intent(inout) :: B(:, :)
        !! Bidiagonal matrix.
        integer, intent(out) :: info
        !! Information flag.
        integer, optional, intent(in) :: kstart
        !! Starting index for the Lanczos factorization (default 1).
        integer, optional, intent(in) :: kend
        !! Final index for the Lanczos factorization (default 1).
        logical, optional, intent(in) :: verbosity
        !! Verbosity control (default `.false.`)
        real(dp), optional, intent(in) :: tol
        !! Tolerance to determine whether invariant subspaces have been computed or not.

        ! Internal variables.
        integer :: k_start, k_end
        logical :: verbose
        real(dp) :: tolerance
        complex(dp) :: alpha, beta, gamma
        integer :: i, j, k, kdim


        info = 0

        ! Krylov subspace dimension.
        kdim = size(U) - 1

        ! Deals with the optional args.
        k_start = optval(kstart, 1)
        k_end   = optval(kend, kdim)
        verbose = optval(verbosity, .false.)
        tolerance = optval(tol, atol_dp)

        ! Lanczos bidiagonalization.
        lanczos : do k = k_start, k_end
            ! Transpose matrix-vector product.
            call A%rmatvec(U(k), V(k))

            ! Full reorthogonalization of the right Krylov subspace.
            do j = 1, k-1
                gamma = V(j)%dot(V(k))
                call V(k)%axpby(cmplx(1.0_dp, 0.0_dp, kind=dp), V(j), -gamma)
            enddo

            ! Normalization step.
            alpha = V(k)%norm() ; B(k, k) = alpha
            if (abs(alpha) > tolerance) then
                call V(k)%scal(1.0_dp/alpha)
            else
                info = k
                exit lanczos
            endif

            ! Matrix-vector product.
            call A%matvec(V(k), U(k+1))

            ! Full re-orthogonalization of the left Krylov subspace.
            do j = 1, k
                gamma = U(j)%dot(U(k+1))
                call U(k+1)%axpby(cmplx(1.0_dp, 0.0_dp, kind=dp), U(j), -gamma)
           enddo

            ! Normalization step
            beta = U(k+1)%norm() ; B(k+1, k) = beta
            if (abs(beta) > tolerance) then
                call U(k+1)%scal(1.0_dp / beta)
            else
                info = k
                exit lanczos
            endif

        enddo lanczos

        return
    end subroutine lanczos_bidiagonalization_cdp


end module lightkrylov_BaseKrylov
