module TestKrylov
    ! Fortran Standard Library.
    use iso_fortran_env
    use stdlib_math, only: is_close, all_close
    use stdlib_linalg, only: eye
    use stdlib_stats, only: median
    ! Testdrive
    use testdrive, only: new_unittest, unittest_type, error_type, check
    ! LightKrylov
    use LightKrylov
    use LightKrylov_Constants
    use LightKrylov_Logger
    use LightKrylov_AbstractVectors
    ! Test Utilities
    use LightKrylov_TestUtils

    implicit none
    
    private

    character(len=*), parameter, private :: this_module      = 'LK_TBKrylov'
    character(len=*), parameter, private :: this_module_long = 'LightKrylov_TestKrylov'

    public :: collect_qr_rsp_testsuite
    public :: collect_arnoldi_rsp_testsuite
    public :: collect_lanczos_bidiag_rsp_testsuite
    public :: collect_lanczos_tridiag_rsp_testsuite

    public :: collect_qr_rdp_testsuite
    public :: collect_arnoldi_rdp_testsuite
    public :: collect_lanczos_bidiag_rdp_testsuite
    public :: collect_lanczos_tridiag_rdp_testsuite

    public :: collect_qr_csp_testsuite
    public :: collect_arnoldi_csp_testsuite
    public :: collect_lanczos_bidiag_csp_testsuite
    public :: collect_lanczos_tridiag_csp_testsuite

    public :: collect_qr_cdp_testsuite
    public :: collect_arnoldi_cdp_testsuite
    public :: collect_lanczos_bidiag_cdp_testsuite
    public :: collect_lanczos_tridiag_cdp_testsuite


contains

    !----------------------------------------------------------------
    !-----     DEFINITIONS OF THE VARIOUS UNIT TESTS FOR QR     -----
    !----------------------------------------------------------------

    subroutine collect_qr_rsp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                        new_unittest("QR factorization", test_qr_factorization_rsp), &
                        new_unittest("Pivoting QR for a rank deficient matrix", test_pivoting_qr_exact_rank_deficiency_rsp) &
                    ]
        return
    end subroutine collect_qr_rsp_testsuite

    subroutine test_qr_factorization_rsp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test Vectors.
        integer, parameter :: kdim = test_size
        type(vector_rsp), allocatable :: A(:)
        ! Upper triangular matrix.
        real(sp) :: R(kdim, kdim) = zero_rsp
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(sp), allocatable :: Adata(:, :), Qdata(:, :)
        real(sp), allocatable :: G(:, :)
        real(sp) :: err
        character(len=256) :: msg

        ! Initialiaze matrix.
        allocate(A(kdim)) ; call init_rand(A)

        ! Get data.
        allocate(Adata(test_size, kdim)) ; call get_data(Adata, A)

        ! In-place QR factorization.
        call qr(A, R, info, tol=atol_sp)
        call check_info(info, 'qr', module=this_module_long, procedure='test_qr_factorization_rsp')

        ! Get data.
        allocate(Qdata(test_size, kdim)) ; call get_data(Qdata, A)

        ! Check correctness.
        err = maxval(abs(Adata - matmul(Qdata, R)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_qr_factorization_rsp', &
                              & info='Factorization', eq='A = Q @ R', context=msg)

        ! Compute Gram matrix associated to the Krylov basis.
        G = Gram(A(:kdim))

        ! Check orthonormality of the computed basis.
        err = norm2(abs(G - eye(kdim, mold=1.0_sp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_qr_factorization_rsp', &
                              & info='Basis orthonormality', eq='Q.H @ Q = I', context=msg)

        return
    end subroutine test_qr_factorization_rsp

    subroutine test_pivoting_qr_exact_rank_deficiency_rsp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test basis.
        type(vector_rsp), allocatable :: A(:)
        ! Krylov subspace dimension.
        integer, parameter :: kdim = 20
        ! Number of zero columns.
        integer, parameter :: nzero = 5
        ! Upper triangular matrix.
        real(sp) :: R(kdim, kdim)
        ! Permutation vector.
        integer :: perm(kdim)
        ! Information flag.
        integer :: info
       
        ! Miscellaneous.
        integer :: k, idx, rk
        real(sp) :: alpha
        logical :: mask(kdim)
        real(sp), allocatable :: Adata(:, :), Qdata(:, :)
        real(sp), allocatable :: G(:, :)
        real(sp) :: err
        character(len=256) :: msg

        ! Effective rank.
        rk = kdim - nzero

        ! Initialize matrix.
        allocate(A(kdim)) ; call init_rand(A)

        ! Add zero vectors at random places.
        mask = .true. ; k = nzero
        do while (k > 0)
            call random_number(alpha)
            idx = 1 + floor(kdim*alpha)
            if (mask(idx)) then
                A(idx)%data = zero_rsp
                mask(idx) = .false.
                k = k-1
            endif
        enddo

        ! Copy data.
        allocate(Adata(test_size, kdim)) ; call get_data(Adata, A)

        ! In-place QR factorization.
        call qr(A, R, perm, info, tol=atol_sp)
        call check_info(info, 'qr_pivot', module=this_module_long, procedure='test_pivoting_qr_exact_rank_deficiency_rsp')

        ! Extract data
        allocate(Qdata(test_size, kdim)) ; call get_data(Qdata, A)
        Adata = Adata(:, perm)

        ! Check correctness.
        err = maxval(abs(Adata - matmul(Qdata, R)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_pivoting_qr_exact_rank_deficiency_rsp', &
                              & info='Factorization', eq='A = Q @ R', context=msg)

        ! Compute Gram matrix associated to the Krylov basis.
        ! allocate(G(kdim, kdim)) ; G = zero_rsp
        G = Gram(A(:kdim))

        ! Check orthonormality of the computed basis.
        err = norm2(abs(G - eye(kdim, mold=1.0_sp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_pivoting_qr_exact_rank_deficiency_rsp', &
                              & info='Basis orthonormality', eq='Q.H @ Q = I', context=msg)

        return
    end subroutine test_pivoting_qr_exact_rank_deficiency_rsp

    subroutine collect_qr_rdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                        new_unittest("QR factorization", test_qr_factorization_rdp), &
                        new_unittest("Pivoting QR for a rank deficient matrix", test_pivoting_qr_exact_rank_deficiency_rdp) &
                    ]
        return
    end subroutine collect_qr_rdp_testsuite

    subroutine test_qr_factorization_rdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test Vectors.
        integer, parameter :: kdim = test_size
        type(vector_rdp), allocatable :: A(:)
        ! Upper triangular matrix.
        real(dp) :: R(kdim, kdim) = zero_rdp
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(dp), allocatable :: Adata(:, :), Qdata(:, :)
        real(dp), allocatable :: G(:, :)
        real(dp) :: err
        character(len=256) :: msg

        ! Initialiaze matrix.
        allocate(A(kdim)) ; call init_rand(A)

        ! Get data.
        allocate(Adata(test_size, kdim)) ; call get_data(Adata, A)

        ! In-place QR factorization.
        call qr(A, R, info, tol=atol_dp)
        call check_info(info, 'qr', module=this_module_long, procedure='test_qr_factorization_rdp')

        ! Get data.
        allocate(Qdata(test_size, kdim)) ; call get_data(Qdata, A)

        ! Check correctness.
        err = maxval(abs(Adata - matmul(Qdata, R)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_qr_factorization_rdp', &
                              & info='Factorization', eq='A = Q @ R', context=msg)

        ! Compute Gram matrix associated to the Krylov basis.
        G = Gram(A(:kdim))

        ! Check orthonormality of the computed basis.
        err = norm2(abs(G - eye(kdim, mold=1.0_dp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_qr_factorization_rdp', &
                              & info='Basis orthonormality', eq='Q.H @ Q = I', context=msg)

        return
    end subroutine test_qr_factorization_rdp

    subroutine test_pivoting_qr_exact_rank_deficiency_rdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test basis.
        type(vector_rdp), allocatable :: A(:)
        ! Krylov subspace dimension.
        integer, parameter :: kdim = 20
        ! Number of zero columns.
        integer, parameter :: nzero = 5
        ! Upper triangular matrix.
        real(dp) :: R(kdim, kdim)
        ! Permutation vector.
        integer :: perm(kdim)
        ! Information flag.
        integer :: info
       
        ! Miscellaneous.
        integer :: k, idx, rk
        real(dp) :: alpha
        logical :: mask(kdim)
        real(dp), allocatable :: Adata(:, :), Qdata(:, :)
        real(dp), allocatable :: G(:, :)
        real(dp) :: err
        character(len=256) :: msg

        ! Effective rank.
        rk = kdim - nzero

        ! Initialize matrix.
        allocate(A(kdim)) ; call init_rand(A)

        ! Add zero vectors at random places.
        mask = .true. ; k = nzero
        do while (k > 0)
            call random_number(alpha)
            idx = 1 + floor(kdim*alpha)
            if (mask(idx)) then
                A(idx)%data = zero_rdp
                mask(idx) = .false.
                k = k-1
            endif
        enddo

        ! Copy data.
        allocate(Adata(test_size, kdim)) ; call get_data(Adata, A)

        ! In-place QR factorization.
        call qr(A, R, perm, info, tol=atol_dp)
        call check_info(info, 'qr_pivot', module=this_module_long, procedure='test_pivoting_qr_exact_rank_deficiency_rdp')

        ! Extract data
        allocate(Qdata(test_size, kdim)) ; call get_data(Qdata, A)
        Adata = Adata(:, perm)

        ! Check correctness.
        err = maxval(abs(Adata - matmul(Qdata, R)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_pivoting_qr_exact_rank_deficiency_rdp', &
                              & info='Factorization', eq='A = Q @ R', context=msg)

        ! Compute Gram matrix associated to the Krylov basis.
        ! allocate(G(kdim, kdim)) ; G = zero_rdp
        G = Gram(A(:kdim))

        ! Check orthonormality of the computed basis.
        err = norm2(abs(G - eye(kdim, mold=1.0_dp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_pivoting_qr_exact_rank_deficiency_rdp', &
                              & info='Basis orthonormality', eq='Q.H @ Q = I', context=msg)

        return
    end subroutine test_pivoting_qr_exact_rank_deficiency_rdp

    subroutine collect_qr_csp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                        new_unittest("QR factorization", test_qr_factorization_csp), &
                        new_unittest("Pivoting QR for a rank deficient matrix", test_pivoting_qr_exact_rank_deficiency_csp) &
                    ]
        return
    end subroutine collect_qr_csp_testsuite

    subroutine test_qr_factorization_csp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test Vectors.
        integer, parameter :: kdim = test_size
        type(vector_csp), allocatable :: A(:)
        ! Upper triangular matrix.
        complex(sp) :: R(kdim, kdim) = zero_csp
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        complex(sp), allocatable :: Adata(:, :), Qdata(:, :)
        complex(sp), allocatable :: G(:, :)
        real(sp) :: err
        character(len=256) :: msg

        ! Initialiaze matrix.
        allocate(A(kdim)) ; call init_rand(A)

        ! Get data.
        allocate(Adata(test_size, kdim)) ; call get_data(Adata, A)

        ! In-place QR factorization.
        call qr(A, R, info, tol=atol_sp)
        call check_info(info, 'qr', module=this_module_long, procedure='test_qr_factorization_csp')

        ! Get data.
        allocate(Qdata(test_size, kdim)) ; call get_data(Qdata, A)

        ! Check correctness.
        err = maxval(abs(Adata - matmul(Qdata, R)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_qr_factorization_csp', &
                              & info='Factorization', eq='A = Q @ R', context=msg)

        ! Compute Gram matrix associated to the Krylov basis.
        G = Gram(A(:kdim))

        ! Check orthonormality of the computed basis.
        err = norm2(abs(G - eye(kdim, mold=1.0_sp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_qr_factorization_csp', &
                              & info='Basis orthonormality', eq='Q.H @ Q = I', context=msg)

        return
    end subroutine test_qr_factorization_csp

    subroutine test_pivoting_qr_exact_rank_deficiency_csp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test basis.
        type(vector_csp), allocatable :: A(:)
        ! Krylov subspace dimension.
        integer, parameter :: kdim = 20
        ! Number of zero columns.
        integer, parameter :: nzero = 5
        ! Upper triangular matrix.
        complex(sp) :: R(kdim, kdim)
        ! Permutation vector.
        integer :: perm(kdim)
        ! Information flag.
        integer :: info
       
        ! Miscellaneous.
        integer :: k, idx, rk
        real(sp) :: alpha
        logical :: mask(kdim)
        complex(sp), allocatable :: Adata(:, :), Qdata(:, :)
        complex(sp), allocatable :: G(:, :)
        real(sp) :: err
        character(len=256) :: msg

        ! Effective rank.
        rk = kdim - nzero

        ! Initialize matrix.
        allocate(A(kdim)) ; call init_rand(A)

        ! Add zero vectors at random places.
        mask = .true. ; k = nzero
        do while (k > 0)
            call random_number(alpha)
            idx = 1 + floor(kdim*alpha)
            if (mask(idx)) then
                A(idx)%data = zero_csp
                mask(idx) = .false.
                k = k-1
            endif
        enddo

        ! Copy data.
        allocate(Adata(test_size, kdim)) ; call get_data(Adata, A)

        ! In-place QR factorization.
        call qr(A, R, perm, info, tol=atol_sp)
        call check_info(info, 'qr_pivot', module=this_module_long, procedure='test_pivoting_qr_exact_rank_deficiency_csp')

        ! Extract data
        allocate(Qdata(test_size, kdim)) ; call get_data(Qdata, A)
        Adata = Adata(:, perm)

        ! Check correctness.
        err = maxval(abs(Adata - matmul(Qdata, R)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_pivoting_qr_exact_rank_deficiency_csp', &
                              & info='Factorization', eq='A = Q @ R', context=msg)

        ! Compute Gram matrix associated to the Krylov basis.
        ! allocate(G(kdim, kdim)) ; G = zero_csp
        G = Gram(A(:kdim))

        ! Check orthonormality of the computed basis.
        err = norm2(abs(G - eye(kdim, mold=1.0_sp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_pivoting_qr_exact_rank_deficiency_csp', &
                              & info='Basis orthonormality', eq='Q.H @ Q = I', context=msg)

        return
    end subroutine test_pivoting_qr_exact_rank_deficiency_csp

    subroutine collect_qr_cdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                        new_unittest("QR factorization", test_qr_factorization_cdp), &
                        new_unittest("Pivoting QR for a rank deficient matrix", test_pivoting_qr_exact_rank_deficiency_cdp) &
                    ]
        return
    end subroutine collect_qr_cdp_testsuite

    subroutine test_qr_factorization_cdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test Vectors.
        integer, parameter :: kdim = test_size
        type(vector_cdp), allocatable :: A(:)
        ! Upper triangular matrix.
        complex(dp) :: R(kdim, kdim) = zero_cdp
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        complex(dp), allocatable :: Adata(:, :), Qdata(:, :)
        complex(dp), allocatable :: G(:, :)
        real(dp) :: err
        character(len=256) :: msg

        ! Initialiaze matrix.
        allocate(A(kdim)) ; call init_rand(A)

        ! Get data.
        allocate(Adata(test_size, kdim)) ; call get_data(Adata, A)

        ! In-place QR factorization.
        call qr(A, R, info, tol=atol_dp)
        call check_info(info, 'qr', module=this_module_long, procedure='test_qr_factorization_cdp')

        ! Get data.
        allocate(Qdata(test_size, kdim)) ; call get_data(Qdata, A)

        ! Check correctness.
        err = maxval(abs(Adata - matmul(Qdata, R)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_qr_factorization_cdp', &
                              & info='Factorization', eq='A = Q @ R', context=msg)

        ! Compute Gram matrix associated to the Krylov basis.
        G = Gram(A(:kdim))

        ! Check orthonormality of the computed basis.
        err = norm2(abs(G - eye(kdim, mold=1.0_dp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_qr_factorization_cdp', &
                              & info='Basis orthonormality', eq='Q.H @ Q = I', context=msg)

        return
    end subroutine test_qr_factorization_cdp

    subroutine test_pivoting_qr_exact_rank_deficiency_cdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test basis.
        type(vector_cdp), allocatable :: A(:)
        ! Krylov subspace dimension.
        integer, parameter :: kdim = 20
        ! Number of zero columns.
        integer, parameter :: nzero = 5
        ! Upper triangular matrix.
        complex(dp) :: R(kdim, kdim)
        ! Permutation vector.
        integer :: perm(kdim)
        ! Information flag.
        integer :: info
       
        ! Miscellaneous.
        integer :: k, idx, rk
        real(dp) :: alpha
        logical :: mask(kdim)
        complex(dp), allocatable :: Adata(:, :), Qdata(:, :)
        complex(dp), allocatable :: G(:, :)
        real(dp) :: err
        character(len=256) :: msg

        ! Effective rank.
        rk = kdim - nzero

        ! Initialize matrix.
        allocate(A(kdim)) ; call init_rand(A)

        ! Add zero vectors at random places.
        mask = .true. ; k = nzero
        do while (k > 0)
            call random_number(alpha)
            idx = 1 + floor(kdim*alpha)
            if (mask(idx)) then
                A(idx)%data = zero_cdp
                mask(idx) = .false.
                k = k-1
            endif
        enddo

        ! Copy data.
        allocate(Adata(test_size, kdim)) ; call get_data(Adata, A)

        ! In-place QR factorization.
        call qr(A, R, perm, info, tol=atol_dp)
        call check_info(info, 'qr_pivot', module=this_module_long, procedure='test_pivoting_qr_exact_rank_deficiency_cdp')

        ! Extract data
        allocate(Qdata(test_size, kdim)) ; call get_data(Qdata, A)
        Adata = Adata(:, perm)

        ! Check correctness.
        err = maxval(abs(Adata - matmul(Qdata, R)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_pivoting_qr_exact_rank_deficiency_cdp', &
                              & info='Factorization', eq='A = Q @ R', context=msg)

        ! Compute Gram matrix associated to the Krylov basis.
        ! allocate(G(kdim, kdim)) ; G = zero_cdp
        G = Gram(A(:kdim))

        ! Check orthonormality of the computed basis.
        err = norm2(abs(G - eye(kdim, mold=1.0_dp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_pivoting_qr_exact_rank_deficiency_cdp', &
                              & info='Basis orthonormality', eq='Q.H @ Q = I', context=msg)

        return
    end subroutine test_pivoting_qr_exact_rank_deficiency_cdp

    
    !--------------------------------------------------------------
    !-----     DEFINITIONS OF THE UNIT-TESTS FOR ARNOLDI      -----
    !--------------------------------------------------------------

    subroutine collect_arnoldi_rsp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("Arnoldi factorization", test_arnoldi_factorization_rsp), &
            new_unittest("Block Arnoldi factorization", test_block_arnoldi_factorization_rsp), &
            new_unittest("Krylov-Schur factorization", test_krylov_schur_rsp) &
                    ]
        return
    end subroutine collect_arnoldi_rsp_testsuite

    subroutine test_arnoldi_factorization_rsp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test linear operator.
        type(linop_rsp), allocatable :: A
        ! Krylov subspace.
        type(vector_rsp), allocatable :: X(:)
        integer, parameter :: kdim = test_size
        ! Hessenberg matrix.
        real(sp), allocatable :: H(:, :)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(sp), allocatable :: Xdata(:, :)
        real(sp), allocatable :: G(:, :)
        real(sp) :: err
        character(len=256) :: msg
           
        ! Initialize linear operator.
        A = linop_rsp() ; call init_rand(A)
        ! Initialize Krylov subspace.
        allocate(X(kdim+1)); call zero_basis(X); call X(1)%rand(ifnorm = .true.)
        allocate(H(kdim+1, kdim)) ; H = zero_rsp
        ! Arnoldi factorization.
        call arnoldi(A, X, H, info, tol=atol_sp)
        call check_info(info, 'arnoldi', module=this_module_long, procedure='test_arnoldi_factorization_rsp')

        ! Check correctness of full factorization.
        allocate(Xdata(test_size, kdim+1)) ; call get_data(Xdata, X)
        err = maxval(abs(matmul(A%data, Xdata(:, :kdim)) - matmul(Xdata, H)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_arnoldi_factorization_rsp', &
                              & info='Factorization', eq='A @ X = X_ @ H_', context=msg)


        ! Compute Gram matrix associated to the Krylov basis.
        ! allocate(G(kdim, kdim)) ; G = zero_rsp
        G = Gram(X(:kdim))

        ! Check orthonormality of the computed basis.
        err = maxval(abs(G - eye(kdim, mold=1.0_sp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_arnoldi_factorization_rsp', &
                              & info='Orthonomality', eq='X.H @ X = I', context=msg)

        return
    end subroutine test_arnoldi_factorization_rsp

    subroutine test_block_arnoldi_factorization_rsp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test Linear operator.
        type(linop_rsp), allocatable :: A
        ! Krylov subspace.
        type(vector_rsp), allocatable :: X(:)
        integer, parameter :: p = 2
        integer, parameter :: kdim = test_size/2
        ! Hessenberg matrix.
        real(sp), allocatable :: H(:, :)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        type(vector_rsp), allocatable :: X0(:)
        real(sp), allocatable :: Xdata(:, :)
        real(sp), allocatable :: G(:, :)
        real(sp) :: err
        character(len=256) :: msg

        ! Initialize linear operator.
        A = linop_rsp() ; call init_rand(A)

        ! Initialize Krylov subspace.
        allocate(X(p*(kdim+1))) ; allocate(X0(p))
        call init_rand(X0) ; call initialize_krylov_subspace(X, X0)
        allocate(H(p*(kdim+1), p*kdim)) ; H = zero_rsp

        ! Arnoldi factorization.
        call arnoldi(A, X, H, info, blksize=p, tol=atol_sp)
        call check_info(info, 'arnoldi', module=this_module_long, procedure='test_block_arnoldi_factorization_rsp')

        ! Check correctness of full factorization.
        allocate(Xdata(test_size, p*(kdim+1))) ; call get_data(Xdata, X)
        err = maxval(abs(matmul(A%data, Xdata(:, :p*kdim)) - matmul(Xdata, H)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_block_arnoldi_factorization_rsp', &
                              & info='Factorization', eq='A @ X = X_ @ H_', context=msg)

        ! Compute Gram matrix associated to the Krylov basis.
        ! allocate(G(p*kdim, p*kdim)) ; G = zero_rsp
        G = Gram(X(:p*kdim))

        ! Check orthonormality of the computed basis.
        err = maxval(abs(G - eye(p*kdim, mold=1.0_sp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_block_arnoldi_factorization_rsp', &
                              & info='Basis orthonormality', eq='X.H @ X = I', context=msg)

        return
    end subroutine test_block_arnoldi_factorization_rsp

    subroutine test_krylov_schur_rsp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test operator.
        type(linop_rsp), allocatable :: A
        ! Krylov subspace.
        type(vector_rsp), allocatable :: X(:)
        ! Krylov subspace dimension.
        integer, parameter :: kdim = 100
        ! Hessenberg matrix.
        real(sp), allocatable :: H(:, :)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        integer :: n
        real(sp), allocatable :: Xdata(:, :)
        real(sp) :: err
        character(len=256) :: msg

        ! Initialize matrix.
        A = linop_rsp() ; call init_rand(A)
        A%data = A%data / norm2(abs(A%data))

        ! Initialize Krylov subspace.
        allocate(X(kdim+1)); call zero_basis(X); call X(1)%rand(ifnorm = .true.)
        allocate(H(kdim+1, kdim)) ; H = zero_rsp

        ! Arnoldi factorization.
        call arnoldi(A, X, H, info, tol=atol_sp)
        call check_info(info, 'arnoldi', module=this_module_long, procedure='test_krylov_schur_rsp')

        ! Krylov-Schur condensation.
        call krylov_schur(n, X, H, select_eigs)

        ! Check correctness.
        allocate(Xdata(test_size, kdim+1)) ; call get_data(Xdata, X)
        err = maxval(abs(matmul(A%data, Xdata(:, :n)) - matmul(Xdata(:, :n+1), H(:n+1, :n))))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_krylov_schur_rsp', &
                              & info='Factorization', eq='A @ X = X_ @ H_', context=msg)

        return
    contains
        function select_eigs(eigvals) result(selected)
            complex(sp), intent(in) :: eigvals(:)
            logical, allocatable :: selected(:)
            selected = abs(eigvals) > median(abs(eigvals))
        end function select_eigs
    end subroutine test_krylov_schur_rsp

    subroutine collect_arnoldi_rdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("Arnoldi factorization", test_arnoldi_factorization_rdp), &
            new_unittest("Block Arnoldi factorization", test_block_arnoldi_factorization_rdp), &
            new_unittest("Krylov-Schur factorization", test_krylov_schur_rdp) &
                    ]
        return
    end subroutine collect_arnoldi_rdp_testsuite

    subroutine test_arnoldi_factorization_rdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test linear operator.
        type(linop_rdp), allocatable :: A
        ! Krylov subspace.
        type(vector_rdp), allocatable :: X(:)
        integer, parameter :: kdim = test_size
        ! Hessenberg matrix.
        real(dp), allocatable :: H(:, :)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(dp), allocatable :: Xdata(:, :)
        real(dp), allocatable :: G(:, :)
        real(dp) :: err
        character(len=256) :: msg
           
        ! Initialize linear operator.
        A = linop_rdp() ; call init_rand(A)
        ! Initialize Krylov subspace.
        allocate(X(kdim+1)); call zero_basis(X); call X(1)%rand(ifnorm = .true.)
        allocate(H(kdim+1, kdim)) ; H = zero_rdp
        ! Arnoldi factorization.
        call arnoldi(A, X, H, info, tol=atol_dp)
        call check_info(info, 'arnoldi', module=this_module_long, procedure='test_arnoldi_factorization_rdp')

        ! Check correctness of full factorization.
        allocate(Xdata(test_size, kdim+1)) ; call get_data(Xdata, X)
        err = maxval(abs(matmul(A%data, Xdata(:, :kdim)) - matmul(Xdata, H)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_arnoldi_factorization_rdp', &
                              & info='Factorization', eq='A @ X = X_ @ H_', context=msg)


        ! Compute Gram matrix associated to the Krylov basis.
        ! allocate(G(kdim, kdim)) ; G = zero_rdp
        G = Gram(X(:kdim))

        ! Check orthonormality of the computed basis.
        err = maxval(abs(G - eye(kdim, mold=1.0_dp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_arnoldi_factorization_rdp', &
                              & info='Orthonomality', eq='X.H @ X = I', context=msg)

        return
    end subroutine test_arnoldi_factorization_rdp

    subroutine test_block_arnoldi_factorization_rdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test Linear operator.
        type(linop_rdp), allocatable :: A
        ! Krylov subspace.
        type(vector_rdp), allocatable :: X(:)
        integer, parameter :: p = 2
        integer, parameter :: kdim = test_size/2
        ! Hessenberg matrix.
        real(dp), allocatable :: H(:, :)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        type(vector_rdp), allocatable :: X0(:)
        real(dp), allocatable :: Xdata(:, :)
        real(dp), allocatable :: G(:, :)
        real(dp) :: err
        character(len=256) :: msg

        ! Initialize linear operator.
        A = linop_rdp() ; call init_rand(A)

        ! Initialize Krylov subspace.
        allocate(X(p*(kdim+1))) ; allocate(X0(p))
        call init_rand(X0) ; call initialize_krylov_subspace(X, X0)
        allocate(H(p*(kdim+1), p*kdim)) ; H = zero_rdp

        ! Arnoldi factorization.
        call arnoldi(A, X, H, info, blksize=p, tol=atol_dp)
        call check_info(info, 'arnoldi', module=this_module_long, procedure='test_block_arnoldi_factorization_rdp')

        ! Check correctness of full factorization.
        allocate(Xdata(test_size, p*(kdim+1))) ; call get_data(Xdata, X)
        err = maxval(abs(matmul(A%data, Xdata(:, :p*kdim)) - matmul(Xdata, H)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_block_arnoldi_factorization_rdp', &
                              & info='Factorization', eq='A @ X = X_ @ H_', context=msg)

        ! Compute Gram matrix associated to the Krylov basis.
        ! allocate(G(p*kdim, p*kdim)) ; G = zero_rdp
        G = Gram(X(:p*kdim))

        ! Check orthonormality of the computed basis.
        err = maxval(abs(G - eye(p*kdim, mold=1.0_dp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_block_arnoldi_factorization_rdp', &
                              & info='Basis orthonormality', eq='X.H @ X = I', context=msg)

        return
    end subroutine test_block_arnoldi_factorization_rdp

    subroutine test_krylov_schur_rdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test operator.
        type(linop_rdp), allocatable :: A
        ! Krylov subspace.
        type(vector_rdp), allocatable :: X(:)
        ! Krylov subspace dimension.
        integer, parameter :: kdim = 100
        ! Hessenberg matrix.
        real(dp), allocatable :: H(:, :)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        integer :: n
        real(dp), allocatable :: Xdata(:, :)
        real(dp) :: err
        character(len=256) :: msg

        ! Initialize matrix.
        A = linop_rdp() ; call init_rand(A)
        A%data = A%data / norm2(abs(A%data))

        ! Initialize Krylov subspace.
        allocate(X(kdim+1)); call zero_basis(X); call X(1)%rand(ifnorm = .true.)
        allocate(H(kdim+1, kdim)) ; H = zero_rdp

        ! Arnoldi factorization.
        call arnoldi(A, X, H, info, tol=atol_dp)
        call check_info(info, 'arnoldi', module=this_module_long, procedure='test_krylov_schur_rdp')

        ! Krylov-Schur condensation.
        call krylov_schur(n, X, H, select_eigs)

        ! Check correctness.
        allocate(Xdata(test_size, kdim+1)) ; call get_data(Xdata, X)
        err = maxval(abs(matmul(A%data, Xdata(:, :n)) - matmul(Xdata(:, :n+1), H(:n+1, :n))))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_krylov_schur_rdp', &
                              & info='Factorization', eq='A @ X = X_ @ H_', context=msg)

        return
    contains
        function select_eigs(eigvals) result(selected)
            complex(dp), intent(in) :: eigvals(:)
            logical, allocatable :: selected(:)
            selected = abs(eigvals) > median(abs(eigvals))
        end function select_eigs
    end subroutine test_krylov_schur_rdp

    subroutine collect_arnoldi_csp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("Arnoldi factorization", test_arnoldi_factorization_csp), &
            new_unittest("Block Arnoldi factorization", test_block_arnoldi_factorization_csp), &
            new_unittest("Krylov-Schur factorization", test_krylov_schur_csp) &
                    ]
        return
    end subroutine collect_arnoldi_csp_testsuite

    subroutine test_arnoldi_factorization_csp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test linear operator.
        type(linop_csp), allocatable :: A
        ! Krylov subspace.
        type(vector_csp), allocatable :: X(:)
        integer, parameter :: kdim = test_size
        ! Hessenberg matrix.
        complex(sp), allocatable :: H(:, :)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        complex(sp), allocatable :: Xdata(:, :)
        complex(sp), allocatable :: G(:, :)
        real(sp) :: err
        character(len=256) :: msg
           
        ! Initialize linear operator.
        A = linop_csp() ; call init_rand(A)
        ! Initialize Krylov subspace.
        allocate(X(kdim+1)); call zero_basis(X); call X(1)%rand(ifnorm = .true.)
        allocate(H(kdim+1, kdim)) ; H = zero_csp
        ! Arnoldi factorization.
        call arnoldi(A, X, H, info, tol=atol_sp)
        call check_info(info, 'arnoldi', module=this_module_long, procedure='test_arnoldi_factorization_csp')

        ! Check correctness of full factorization.
        allocate(Xdata(test_size, kdim+1)) ; call get_data(Xdata, X)
        err = maxval(abs(matmul(A%data, Xdata(:, :kdim)) - matmul(Xdata, H)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_arnoldi_factorization_csp', &
                              & info='Factorization', eq='A @ X = X_ @ H_', context=msg)


        ! Compute Gram matrix associated to the Krylov basis.
        ! allocate(G(kdim, kdim)) ; G = zero_csp
        G = Gram(X(:kdim))

        ! Check orthonormality of the computed basis.
        err = maxval(abs(G - eye(kdim, mold=1.0_sp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_arnoldi_factorization_csp', &
                              & info='Orthonomality', eq='X.H @ X = I', context=msg)

        return
    end subroutine test_arnoldi_factorization_csp

    subroutine test_block_arnoldi_factorization_csp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test Linear operator.
        type(linop_csp), allocatable :: A
        ! Krylov subspace.
        type(vector_csp), allocatable :: X(:)
        integer, parameter :: p = 2
        integer, parameter :: kdim = test_size/2
        ! Hessenberg matrix.
        complex(sp), allocatable :: H(:, :)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        type(vector_csp), allocatable :: X0(:)
        complex(sp), allocatable :: Xdata(:, :)
        complex(sp), allocatable :: G(:, :)
        real(sp) :: err
        character(len=256) :: msg

        ! Initialize linear operator.
        A = linop_csp() ; call init_rand(A)

        ! Initialize Krylov subspace.
        allocate(X(p*(kdim+1))) ; allocate(X0(p))
        call init_rand(X0) ; call initialize_krylov_subspace(X, X0)
        allocate(H(p*(kdim+1), p*kdim)) ; H = zero_csp

        ! Arnoldi factorization.
        call arnoldi(A, X, H, info, blksize=p, tol=atol_sp)
        call check_info(info, 'arnoldi', module=this_module_long, procedure='test_block_arnoldi_factorization_csp')

        ! Check correctness of full factorization.
        allocate(Xdata(test_size, p*(kdim+1))) ; call get_data(Xdata, X)
        err = maxval(abs(matmul(A%data, Xdata(:, :p*kdim)) - matmul(Xdata, H)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_block_arnoldi_factorization_csp', &
                              & info='Factorization', eq='A @ X = X_ @ H_', context=msg)

        ! Compute Gram matrix associated to the Krylov basis.
        ! allocate(G(p*kdim, p*kdim)) ; G = zero_csp
        G = Gram(X(:p*kdim))

        ! Check orthonormality of the computed basis.
        err = maxval(abs(G - eye(p*kdim, mold=1.0_sp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_block_arnoldi_factorization_csp', &
                              & info='Basis orthonormality', eq='X.H @ X = I', context=msg)

        return
    end subroutine test_block_arnoldi_factorization_csp

    subroutine test_krylov_schur_csp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test operator.
        type(linop_csp), allocatable :: A
        ! Krylov subspace.
        type(vector_csp), allocatable :: X(:)
        ! Krylov subspace dimension.
        integer, parameter :: kdim = 100
        ! Hessenberg matrix.
        complex(sp), allocatable :: H(:, :)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        integer :: n
        complex(sp), allocatable :: Xdata(:, :)
        real(sp) :: err
        character(len=256) :: msg

        ! Initialize matrix.
        A = linop_csp() ; call init_rand(A)
        A%data = A%data / norm2(abs(A%data))

        ! Initialize Krylov subspace.
        allocate(X(kdim+1)); call zero_basis(X); call X(1)%rand(ifnorm = .true.)
        allocate(H(kdim+1, kdim)) ; H = zero_csp

        ! Arnoldi factorization.
        call arnoldi(A, X, H, info, tol=atol_sp)
        call check_info(info, 'arnoldi', module=this_module_long, procedure='test_krylov_schur_csp')

        ! Krylov-Schur condensation.
        call krylov_schur(n, X, H, select_eigs)

        ! Check correctness.
        allocate(Xdata(test_size, kdim+1)) ; call get_data(Xdata, X)
        err = maxval(abs(matmul(A%data, Xdata(:, :n)) - matmul(Xdata(:, :n+1), H(:n+1, :n))))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_krylov_schur_csp', &
                              & info='Factorization', eq='A @ X = X_ @ H_', context=msg)

        return
    contains
        function select_eigs(eigvals) result(selected)
            complex(sp), intent(in) :: eigvals(:)
            logical, allocatable :: selected(:)
            selected = abs(eigvals) > median(abs(eigvals))
        end function select_eigs
    end subroutine test_krylov_schur_csp

    subroutine collect_arnoldi_cdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("Arnoldi factorization", test_arnoldi_factorization_cdp), &
            new_unittest("Block Arnoldi factorization", test_block_arnoldi_factorization_cdp), &
            new_unittest("Krylov-Schur factorization", test_krylov_schur_cdp) &
                    ]
        return
    end subroutine collect_arnoldi_cdp_testsuite

    subroutine test_arnoldi_factorization_cdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test linear operator.
        type(linop_cdp), allocatable :: A
        ! Krylov subspace.
        type(vector_cdp), allocatable :: X(:)
        integer, parameter :: kdim = test_size
        ! Hessenberg matrix.
        complex(dp), allocatable :: H(:, :)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        complex(dp), allocatable :: Xdata(:, :)
        complex(dp), allocatable :: G(:, :)
        real(dp) :: err
        character(len=256) :: msg
           
        ! Initialize linear operator.
        A = linop_cdp() ; call init_rand(A)
        ! Initialize Krylov subspace.
        allocate(X(kdim+1)); call zero_basis(X); call X(1)%rand(ifnorm = .true.)
        allocate(H(kdim+1, kdim)) ; H = zero_cdp
        ! Arnoldi factorization.
        call arnoldi(A, X, H, info, tol=atol_dp)
        call check_info(info, 'arnoldi', module=this_module_long, procedure='test_arnoldi_factorization_cdp')

        ! Check correctness of full factorization.
        allocate(Xdata(test_size, kdim+1)) ; call get_data(Xdata, X)
        err = maxval(abs(matmul(A%data, Xdata(:, :kdim)) - matmul(Xdata, H)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_arnoldi_factorization_cdp', &
                              & info='Factorization', eq='A @ X = X_ @ H_', context=msg)


        ! Compute Gram matrix associated to the Krylov basis.
        ! allocate(G(kdim, kdim)) ; G = zero_cdp
        G = Gram(X(:kdim))

        ! Check orthonormality of the computed basis.
        err = maxval(abs(G - eye(kdim, mold=1.0_dp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_arnoldi_factorization_cdp', &
                              & info='Orthonomality', eq='X.H @ X = I', context=msg)

        return
    end subroutine test_arnoldi_factorization_cdp

    subroutine test_block_arnoldi_factorization_cdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test Linear operator.
        type(linop_cdp), allocatable :: A
        ! Krylov subspace.
        type(vector_cdp), allocatable :: X(:)
        integer, parameter :: p = 2
        integer, parameter :: kdim = test_size/2
        ! Hessenberg matrix.
        complex(dp), allocatable :: H(:, :)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        type(vector_cdp), allocatable :: X0(:)
        complex(dp), allocatable :: Xdata(:, :)
        complex(dp), allocatable :: G(:, :)
        real(dp) :: err
        character(len=256) :: msg

        ! Initialize linear operator.
        A = linop_cdp() ; call init_rand(A)

        ! Initialize Krylov subspace.
        allocate(X(p*(kdim+1))) ; allocate(X0(p))
        call init_rand(X0) ; call initialize_krylov_subspace(X, X0)
        allocate(H(p*(kdim+1), p*kdim)) ; H = zero_cdp

        ! Arnoldi factorization.
        call arnoldi(A, X, H, info, blksize=p, tol=atol_dp)
        call check_info(info, 'arnoldi', module=this_module_long, procedure='test_block_arnoldi_factorization_cdp')

        ! Check correctness of full factorization.
        allocate(Xdata(test_size, p*(kdim+1))) ; call get_data(Xdata, X)
        err = maxval(abs(matmul(A%data, Xdata(:, :p*kdim)) - matmul(Xdata, H)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_block_arnoldi_factorization_cdp', &
                              & info='Factorization', eq='A @ X = X_ @ H_', context=msg)

        ! Compute Gram matrix associated to the Krylov basis.
        ! allocate(G(p*kdim, p*kdim)) ; G = zero_cdp
        G = Gram(X(:p*kdim))

        ! Check orthonormality of the computed basis.
        err = maxval(abs(G - eye(p*kdim, mold=1.0_dp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_block_arnoldi_factorization_cdp', &
                              & info='Basis orthonormality', eq='X.H @ X = I', context=msg)

        return
    end subroutine test_block_arnoldi_factorization_cdp

    subroutine test_krylov_schur_cdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test operator.
        type(linop_cdp), allocatable :: A
        ! Krylov subspace.
        type(vector_cdp), allocatable :: X(:)
        ! Krylov subspace dimension.
        integer, parameter :: kdim = 100
        ! Hessenberg matrix.
        complex(dp), allocatable :: H(:, :)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        integer :: n
        complex(dp), allocatable :: Xdata(:, :)
        real(dp) :: err
        character(len=256) :: msg

        ! Initialize matrix.
        A = linop_cdp() ; call init_rand(A)
        A%data = A%data / norm2(abs(A%data))

        ! Initialize Krylov subspace.
        allocate(X(kdim+1)); call zero_basis(X); call X(1)%rand(ifnorm = .true.)
        allocate(H(kdim+1, kdim)) ; H = zero_cdp

        ! Arnoldi factorization.
        call arnoldi(A, X, H, info, tol=atol_dp)
        call check_info(info, 'arnoldi', module=this_module_long, procedure='test_krylov_schur_cdp')

        ! Krylov-Schur condensation.
        call krylov_schur(n, X, H, select_eigs)

        ! Check correctness.
        allocate(Xdata(test_size, kdim+1)) ; call get_data(Xdata, X)
        err = maxval(abs(matmul(A%data, Xdata(:, :n)) - matmul(Xdata(:, :n+1), H(:n+1, :n))))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_krylov_schur_cdp', &
                              & info='Factorization', eq='A @ X = X_ @ H_', context=msg)

        return
    contains
        function select_eigs(eigvals) result(selected)
            complex(dp), intent(in) :: eigvals(:)
            logical, allocatable :: selected(:)
            selected = abs(eigvals) > median(abs(eigvals))
        end function select_eigs
    end subroutine test_krylov_schur_cdp


    !------------------------------------------------------------------------------
    !-----     DEFINITION OF THE UNIT-TESTS FOR LANCZOS BIDIAGONALIZATION     -----
    !------------------------------------------------------------------------------

    subroutine collect_lanczos_bidiag_rsp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("Lanczos Bidiagonalization", test_lanczos_bidiag_factorization_rsp) &
                ]
        return
    end subroutine collect_lanczos_bidiag_rsp_testsuite

    subroutine test_lanczos_bidiag_factorization_rsp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test Linear Operator
        type(linop_rsp), allocatable :: A
        ! Left and right Krylov bases.
        type(vector_rsp), allocatable :: U(:), V(:)
        ! Krylov subspace dimension.
        integer, parameter :: kdim = test_size
        ! Bidiagonal matrix.
        real(sp), allocatable :: B(:, :)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(sp), allocatable :: Udata(:, :), Vdata(:, :)
        real(sp), allocatable :: G(:, :)
        real(sp) :: err
        character(len=256) :: msg

        ! Initialize linear operator.
        A = linop_rsp() ; call init_rand(A)

        ! Initialize Krylov subspaces.
        allocate(U(kdim+1), V(kdim+1), B(kdim+1,kdim))
        call zero_basis(U); call U(1)%rand(ifnorm = .true.)
        call zero_basis(V)
        B = zero_rsp

        ! Lanczos bidiagonalization.
        call bidiagonalization(A, U, V, B, info, tol=atol_sp)
        call check_info(info, 'bidiagonalization', module=this_module_long, &
                        & procedure='test_lanczos_bidiag_factorization_rsp')

        ! Check correctness.
        allocate(Udata(test_size, kdim+1)) ; call get_data(Udata, U)
        allocate(Vdata(test_size, kdim+1)) ; call get_data(Vdata, V)

        err = maxval(abs(matmul(A%data, Vdata(:, :kdim)) - matmul(Udata, B)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_lanczos_bidiag_factorization_rsp', &
                              & info='Factorization', eq='A @ V = U_ @ B_', context=msg)

        ! Compute Gram matrix associated to the left Krylov basis.
        G = Gram(U(:kdim))

        ! Check orthonormality of the left basis.
        err = maxval(abs(G - eye(kdim, mold=1.0_sp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_lanczos_bidiag_factorization_rsp', &
                              & info='Basis orthonormality (left)', eq='U.H @ U = I', context=msg)

        ! Compute Gram matrix associated to the right Krylov basis.
        G = Gram(V(:kdim))

        ! Check orthonormality of the right basis.
        err = maxval(abs(G - eye(kdim, mold=1.0_sp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_lanczos_bidiag_factorization_rsp', &
                              & info='Basis orthonormality (right)', eq='V.H @ V = I', context=msg)

        return
    end subroutine test_lanczos_bidiag_factorization_rsp

    subroutine collect_lanczos_bidiag_rdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("Lanczos Bidiagonalization", test_lanczos_bidiag_factorization_rdp) &
                ]
        return
    end subroutine collect_lanczos_bidiag_rdp_testsuite

    subroutine test_lanczos_bidiag_factorization_rdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test Linear Operator
        type(linop_rdp), allocatable :: A
        ! Left and right Krylov bases.
        type(vector_rdp), allocatable :: U(:), V(:)
        ! Krylov subspace dimension.
        integer, parameter :: kdim = test_size
        ! Bidiagonal matrix.
        real(dp), allocatable :: B(:, :)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(dp), allocatable :: Udata(:, :), Vdata(:, :)
        real(dp), allocatable :: G(:, :)
        real(dp) :: err
        character(len=256) :: msg

        ! Initialize linear operator.
        A = linop_rdp() ; call init_rand(A)

        ! Initialize Krylov subspaces.
        allocate(U(kdim+1), V(kdim+1), B(kdim+1,kdim))
        call zero_basis(U); call U(1)%rand(ifnorm = .true.)
        call zero_basis(V)
        B = zero_rdp

        ! Lanczos bidiagonalization.
        call bidiagonalization(A, U, V, B, info, tol=atol_dp)
        call check_info(info, 'bidiagonalization', module=this_module_long, &
                        & procedure='test_lanczos_bidiag_factorization_rdp')

        ! Check correctness.
        allocate(Udata(test_size, kdim+1)) ; call get_data(Udata, U)
        allocate(Vdata(test_size, kdim+1)) ; call get_data(Vdata, V)

        err = maxval(abs(matmul(A%data, Vdata(:, :kdim)) - matmul(Udata, B)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_lanczos_bidiag_factorization_rdp', &
                              & info='Factorization', eq='A @ V = U_ @ B_', context=msg)

        ! Compute Gram matrix associated to the left Krylov basis.
        G = Gram(U(:kdim))

        ! Check orthonormality of the left basis.
        err = maxval(abs(G - eye(kdim, mold=1.0_dp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_lanczos_bidiag_factorization_rdp', &
                              & info='Basis orthonormality (left)', eq='U.H @ U = I', context=msg)

        ! Compute Gram matrix associated to the right Krylov basis.
        G = Gram(V(:kdim))

        ! Check orthonormality of the right basis.
        err = maxval(abs(G - eye(kdim, mold=1.0_dp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_lanczos_bidiag_factorization_rdp', &
                              & info='Basis orthonormality (right)', eq='V.H @ V = I', context=msg)

        return
    end subroutine test_lanczos_bidiag_factorization_rdp

    subroutine collect_lanczos_bidiag_csp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("Lanczos Bidiagonalization", test_lanczos_bidiag_factorization_csp) &
                ]
        return
    end subroutine collect_lanczos_bidiag_csp_testsuite

    subroutine test_lanczos_bidiag_factorization_csp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test Linear Operator
        type(linop_csp), allocatable :: A
        ! Left and right Krylov bases.
        type(vector_csp), allocatable :: U(:), V(:)
        ! Krylov subspace dimension.
        integer, parameter :: kdim = test_size
        ! Bidiagonal matrix.
        complex(sp), allocatable :: B(:, :)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        complex(sp), allocatable :: Udata(:, :), Vdata(:, :)
        complex(sp), allocatable :: G(:, :)
        real(sp) :: err
        character(len=256) :: msg

        ! Initialize linear operator.
        A = linop_csp() ; call init_rand(A)

        ! Initialize Krylov subspaces.
        allocate(U(kdim+1), V(kdim+1), B(kdim+1,kdim))
        call zero_basis(U); call U(1)%rand(ifnorm = .true.)
        call zero_basis(V)
        B = zero_csp

        ! Lanczos bidiagonalization.
        call bidiagonalization(A, U, V, B, info, tol=atol_sp)
        call check_info(info, 'bidiagonalization', module=this_module_long, &
                        & procedure='test_lanczos_bidiag_factorization_csp')

        ! Check correctness.
        allocate(Udata(test_size, kdim+1)) ; call get_data(Udata, U)
        allocate(Vdata(test_size, kdim+1)) ; call get_data(Vdata, V)

        err = maxval(abs(matmul(A%data, Vdata(:, :kdim)) - matmul(Udata, B)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_lanczos_bidiag_factorization_csp', &
                              & info='Factorization', eq='A @ V = U_ @ B_', context=msg)

        ! Compute Gram matrix associated to the left Krylov basis.
        G = Gram(U(:kdim))

        ! Check orthonormality of the left basis.
        err = maxval(abs(G - eye(kdim, mold=1.0_sp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_lanczos_bidiag_factorization_csp', &
                              & info='Basis orthonormality (left)', eq='U.H @ U = I', context=msg)

        ! Compute Gram matrix associated to the right Krylov basis.
        G = Gram(V(:kdim))

        ! Check orthonormality of the right basis.
        err = maxval(abs(G - eye(kdim, mold=1.0_sp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_lanczos_bidiag_factorization_csp', &
                              & info='Basis orthonormality (right)', eq='V.H @ V = I', context=msg)

        return
    end subroutine test_lanczos_bidiag_factorization_csp

    subroutine collect_lanczos_bidiag_cdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("Lanczos Bidiagonalization", test_lanczos_bidiag_factorization_cdp) &
                ]
        return
    end subroutine collect_lanczos_bidiag_cdp_testsuite

    subroutine test_lanczos_bidiag_factorization_cdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test Linear Operator
        type(linop_cdp), allocatable :: A
        ! Left and right Krylov bases.
        type(vector_cdp), allocatable :: U(:), V(:)
        ! Krylov subspace dimension.
        integer, parameter :: kdim = test_size
        ! Bidiagonal matrix.
        complex(dp), allocatable :: B(:, :)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        complex(dp), allocatable :: Udata(:, :), Vdata(:, :)
        complex(dp), allocatable :: G(:, :)
        real(dp) :: err
        character(len=256) :: msg

        ! Initialize linear operator.
        A = linop_cdp() ; call init_rand(A)

        ! Initialize Krylov subspaces.
        allocate(U(kdim+1), V(kdim+1), B(kdim+1,kdim))
        call zero_basis(U); call U(1)%rand(ifnorm = .true.)
        call zero_basis(V)
        B = zero_cdp

        ! Lanczos bidiagonalization.
        call bidiagonalization(A, U, V, B, info, tol=atol_dp)
        call check_info(info, 'bidiagonalization', module=this_module_long, &
                        & procedure='test_lanczos_bidiag_factorization_cdp')

        ! Check correctness.
        allocate(Udata(test_size, kdim+1)) ; call get_data(Udata, U)
        allocate(Vdata(test_size, kdim+1)) ; call get_data(Vdata, V)

        err = maxval(abs(matmul(A%data, Vdata(:, :kdim)) - matmul(Udata, B)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_lanczos_bidiag_factorization_cdp', &
                              & info='Factorization', eq='A @ V = U_ @ B_', context=msg)

        ! Compute Gram matrix associated to the left Krylov basis.
        G = Gram(U(:kdim))

        ! Check orthonormality of the left basis.
        err = maxval(abs(G - eye(kdim, mold=1.0_dp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_lanczos_bidiag_factorization_cdp', &
                              & info='Basis orthonormality (left)', eq='U.H @ U = I', context=msg)

        ! Compute Gram matrix associated to the right Krylov basis.
        G = Gram(V(:kdim))

        ! Check orthonormality of the right basis.
        err = maxval(abs(G - eye(kdim, mold=1.0_dp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_lanczos_bidiag_factorization_cdp', &
                              & info='Basis orthonormality (right)', eq='V.H @ V = I', context=msg)

        return
    end subroutine test_lanczos_bidiag_factorization_cdp


    !-------------------------------------------------------------------------------
    !-----     DEFINITION OF THE UNIT TESTS FOR LANCZOS TRIDIAGONALIZATION     -----
    !-------------------------------------------------------------------------------

    subroutine collect_lanczos_tridiag_rsp_testsuite(testsuite)
        ! Collection of unit tests.
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
             new_unittest("Lanczos Tridiagonalization", test_lanczos_tridiag_factorization_rsp) &
            ]

        return
    end subroutine collect_lanczos_tridiag_rsp_testsuite

    subroutine test_lanczos_tridiag_factorization_rsp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test matrix.
        type(spd_linop_rsp), allocatable :: A
        ! Krylov subspace.
        type(vector_rsp), allocatable :: X(:)
        ! Krylov subspace dimension.
        integer, parameter :: kdim = test_size
        ! Tridiagonal matrix.
        real(sp), allocatable :: T(:, :)
        ! Information flag.
        integer :: info

        ! Internal variables.
        real(sp), allocatable :: Xdata(:, :)
        real(sp), allocatable :: G(:, :)
        real(sp) :: err
        character(len=256) :: msg

        ! Initialize tridiagonal matrix.
        allocate(T(kdim+1, kdim)) ; T = zero_rsp

        ! Initialize operator.
        A = spd_linop_rsp()
        call init_rand(A)

        ! Initialize Krylov subspace.
        allocate(X(kdim+1)); call zero_basis(X); call X(1)%rand(ifnorm = .true.)

        ! Lanczos factorization.
        call lanczos(A, X, T, info, tol=atol_sp)
        call check_info(info, 'lanczos', module=this_module_long, & 
                        & procedure='test_lanczos_tridiag_factorization_rsp')

        ! Check correctness.
        allocate(Xdata(test_size, kdim+1)) ; call get_data(Xdata, X)

        ! Infinity-norm check.
        err = maxval(abs(matmul(A%data, Xdata(:, :kdim)) - matmul(Xdata, T)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_lanczos_tridiag_factorization_rsp', &
                                 & info='Factorization', eq='A @ X = X_ @ T_', context=msg)

        ! Compute Gram matrix associated to the right Krylov basis.
        ! allocate(G(kdim, kdim)) ; G = zero_rsp
        G = Gram(X(:kdim))

        ! Check orthonormality of the Krylov basis.
        err = maxval(abs(G - eye(kdim, mold=1.0_sp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_lanczos_tridiag_factorization_rsp', &
                                 & info='Orthonomality', eq='X.H @ X = I', context=msg)

        return
    end subroutine test_lanczos_tridiag_factorization_rsp

    subroutine collect_lanczos_tridiag_rdp_testsuite(testsuite)
        ! Collection of unit tests.
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
             new_unittest("Lanczos Tridiagonalization", test_lanczos_tridiag_factorization_rdp) &
            ]

        return
    end subroutine collect_lanczos_tridiag_rdp_testsuite

    subroutine test_lanczos_tridiag_factorization_rdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test matrix.
        type(spd_linop_rdp), allocatable :: A
        ! Krylov subspace.
        type(vector_rdp), allocatable :: X(:)
        ! Krylov subspace dimension.
        integer, parameter :: kdim = test_size
        ! Tridiagonal matrix.
        real(dp), allocatable :: T(:, :)
        ! Information flag.
        integer :: info

        ! Internal variables.
        real(dp), allocatable :: Xdata(:, :)
        real(dp), allocatable :: G(:, :)
        real(dp) :: err
        character(len=256) :: msg

        ! Initialize tridiagonal matrix.
        allocate(T(kdim+1, kdim)) ; T = zero_rdp

        ! Initialize operator.
        A = spd_linop_rdp()
        call init_rand(A)

        ! Initialize Krylov subspace.
        allocate(X(kdim+1)); call zero_basis(X); call X(1)%rand(ifnorm = .true.)

        ! Lanczos factorization.
        call lanczos(A, X, T, info, tol=atol_dp)
        call check_info(info, 'lanczos', module=this_module_long, & 
                        & procedure='test_lanczos_tridiag_factorization_rdp')

        ! Check correctness.
        allocate(Xdata(test_size, kdim+1)) ; call get_data(Xdata, X)

        ! Infinity-norm check.
        err = maxval(abs(matmul(A%data, Xdata(:, :kdim)) - matmul(Xdata, T)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_lanczos_tridiag_factorization_rdp', &
                                 & info='Factorization', eq='A @ X = X_ @ T_', context=msg)

        ! Compute Gram matrix associated to the right Krylov basis.
        ! allocate(G(kdim, kdim)) ; G = zero_rdp
        G = Gram(X(:kdim))

        ! Check orthonormality of the Krylov basis.
        err = maxval(abs(G - eye(kdim, mold=1.0_dp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_lanczos_tridiag_factorization_rdp', &
                                 & info='Orthonomality', eq='X.H @ X = I', context=msg)

        return
    end subroutine test_lanczos_tridiag_factorization_rdp

    subroutine collect_lanczos_tridiag_csp_testsuite(testsuite)
        ! Collection of unit tests.
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
             new_unittest("Lanczos Tridiagonalization", test_lanczos_tridiag_factorization_csp) &
            ]

        return
    end subroutine collect_lanczos_tridiag_csp_testsuite

    subroutine test_lanczos_tridiag_factorization_csp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test matrix.
        type(hermitian_linop_csp), allocatable :: A
        ! Krylov subspace.
        type(vector_csp), allocatable :: X(:)
        ! Krylov subspace dimension.
        integer, parameter :: kdim = test_size
        ! Tridiagonal matrix.
        complex(sp), allocatable :: T(:, :)
        ! Information flag.
        integer :: info

        ! Internal variables.
        complex(sp), allocatable :: Xdata(:, :)
        complex(sp), allocatable :: G(:, :)
        real(sp) :: err
        character(len=256) :: msg

        ! Initialize tridiagonal matrix.
        allocate(T(kdim+1, kdim)) ; T = zero_csp

        ! Initialize operator.
        A = hermitian_linop_csp()
        call init_rand(A)

        ! Initialize Krylov subspace.
        allocate(X(kdim+1)); call zero_basis(X); call X(1)%rand(ifnorm = .true.)

        ! Lanczos factorization.
        call lanczos(A, X, T, info, tol=atol_sp)
        call check_info(info, 'lanczos', module=this_module_long, & 
                        & procedure='test_lanczos_tridiag_factorization_csp')

        ! Check correctness.
        allocate(Xdata(test_size, kdim+1)) ; call get_data(Xdata, X)

        ! Infinity-norm check.
        err = maxval(abs(matmul(A%data, Xdata(:, :kdim)) - matmul(Xdata, T)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_lanczos_tridiag_factorization_csp', &
                                 & info='Factorization', eq='A @ X = X_ @ T_', context=msg)

        ! Compute Gram matrix associated to the right Krylov basis.
        ! allocate(G(kdim, kdim)) ; G = zero_csp
        G = Gram(X(:kdim))

        ! Check orthonormality of the Krylov basis.
        err = maxval(abs(G - eye(kdim, mold=1.0_sp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_lanczos_tridiag_factorization_csp', &
                                 & info='Orthonomality', eq='X.H @ X = I', context=msg)

        return
    end subroutine test_lanczos_tridiag_factorization_csp

    subroutine collect_lanczos_tridiag_cdp_testsuite(testsuite)
        ! Collection of unit tests.
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
             new_unittest("Lanczos Tridiagonalization", test_lanczos_tridiag_factorization_cdp) &
            ]

        return
    end subroutine collect_lanczos_tridiag_cdp_testsuite

    subroutine test_lanczos_tridiag_factorization_cdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test matrix.
        type(hermitian_linop_cdp), allocatable :: A
        ! Krylov subspace.
        type(vector_cdp), allocatable :: X(:)
        ! Krylov subspace dimension.
        integer, parameter :: kdim = test_size
        ! Tridiagonal matrix.
        complex(dp), allocatable :: T(:, :)
        ! Information flag.
        integer :: info

        ! Internal variables.
        complex(dp), allocatable :: Xdata(:, :)
        complex(dp), allocatable :: G(:, :)
        real(dp) :: err
        character(len=256) :: msg

        ! Initialize tridiagonal matrix.
        allocate(T(kdim+1, kdim)) ; T = zero_cdp

        ! Initialize operator.
        A = hermitian_linop_cdp()
        call init_rand(A)

        ! Initialize Krylov subspace.
        allocate(X(kdim+1)); call zero_basis(X); call X(1)%rand(ifnorm = .true.)

        ! Lanczos factorization.
        call lanczos(A, X, T, info, tol=atol_dp)
        call check_info(info, 'lanczos', module=this_module_long, & 
                        & procedure='test_lanczos_tridiag_factorization_cdp')

        ! Check correctness.
        allocate(Xdata(test_size, kdim+1)) ; call get_data(Xdata, X)

        ! Infinity-norm check.
        err = maxval(abs(matmul(A%data, Xdata(:, :kdim)) - matmul(Xdata, T)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_lanczos_tridiag_factorization_cdp', &
                                 & info='Factorization', eq='A @ X = X_ @ T_', context=msg)

        ! Compute Gram matrix associated to the right Krylov basis.
        ! allocate(G(kdim, kdim)) ; G = zero_cdp
        G = Gram(X(:kdim))

        ! Check orthonormality of the Krylov basis.
        err = maxval(abs(G - eye(kdim, mold=1.0_dp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_lanczos_tridiag_factorization_cdp', &
                                 & info='Orthonomality', eq='X.H @ X = I', context=msg)

        return
    end subroutine test_lanczos_tridiag_factorization_cdp

end module TestKrylov
