module TestKrylov
    ! Fortran Standard Library.
    use iso_fortran_env, only: output_unit
    use stdlib_math, only: is_close, all_close
    use stdlib_linalg, only: eye, hermitian, is_hessenberg
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
    use TestUtils
    use stdlib_io_npy, only: save_npy

    implicit none (type, external)

    private

    character(len=*), parameter, private :: this_module      = 'LK_TBKrylov'
    character(len=*), parameter, private :: this_module_long = 'LightKrylov_TestKrylov'

    public :: collect_qr_rsp_testsuite
    public :: collect_arnoldi_rsp_testsuite
    public :: collect_lanczos_bidiag_rsp_testsuite
    public :: collect_lanczos_tridiag_rsp_testsuite
    public :: collect_ssy_tridiag_rsp_testsuite
    public :: collect_krylov_utilities_rsp_testsuite

    public :: collect_qr_rdp_testsuite
    public :: collect_arnoldi_rdp_testsuite
    public :: collect_lanczos_bidiag_rdp_testsuite
    public :: collect_lanczos_tridiag_rdp_testsuite
    public :: collect_ssy_tridiag_rdp_testsuite
    public :: collect_krylov_utilities_rdp_testsuite

    public :: collect_qr_csp_testsuite
    public :: collect_arnoldi_csp_testsuite
    public :: collect_lanczos_bidiag_csp_testsuite
    public :: collect_lanczos_tridiag_csp_testsuite
    public :: collect_ssy_tridiag_csp_testsuite
    public :: collect_krylov_utilities_csp_testsuite

    public :: collect_qr_cdp_testsuite
    public :: collect_arnoldi_cdp_testsuite
    public :: collect_lanczos_bidiag_cdp_testsuite
    public :: collect_lanczos_tridiag_cdp_testsuite
    public :: collect_ssy_tridiag_cdp_testsuite
    public :: collect_krylov_utilities_cdp_testsuite


contains

    !----------------------------------------------------------------
    !-----     DEFINITIONS OF THE VARIOUS UNIT TESTS FOR QR     -----
    !----------------------------------------------------------------

    subroutine collect_qr_rsp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                        new_unittest("QR invalid inputs", test_qr_invalid_inputs_rsp), &
                        new_unittest("QR factorization", test_qr_factorization_rsp), &
                        new_unittest("QR rank deficient", test_qr_rank_deficient_rsp), &
                        new_unittest("QR single vector", test_qr_single_vector_rsp) &
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
        real(sp) :: R(kdim, kdim)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(sp), allocatable :: Adata(:, :), Qdata(:, :)
        real(sp), allocatable :: G(:, :)
        real(sp) :: err
        integer :: perm(kdim), i
        character(len=256) :: msg

        ! Initialiaze matrix.
        allocate(A(kdim)) ; call init_rand(A); R = zero_rsp

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


        call init_rand(A); R = zero_rsp ; call get_data(Adata, A)

        ! In-place Pivoted QR factorization.
        call qr(A, R, perm, info, tol=atol_sp)
        call check_info(info, 'qr_pivot', module=this_module_long, procedure='test_qr_factorization_rsp')

        ! Get data.
        call get_data(Qdata, A)

        ! Check correctness.
        err = maxval(abs(Adata(:, perm) - matmul(Qdata, R)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_qr_factorization_rsp', &
                              & info='Pivoted factorization', eq='AP = Q @ R', context=msg)

        ! Compute Gram matrix associated to the Krylov basis.
        G = Gram(A(:kdim))

        ! Check orthonormality of the computed basis.
        err = norm2(abs(G - eye(kdim, mold=1.0_sp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_qr_factorization_rsp', &
                              & info='Pivoted basis orthonormality', eq='Q.H @ Q = I', context=msg)

        ! Check that diagonal entries are non-increasing in magnitude.
        ! |R(1,1)| >= |R(2,2)| >= ... >= |R(kdim,kdim)|
        do i = 1, kdim-1
            err = abs(R(i+1, i+1)) - abs(R(i, i))
            call check(error, err <= rtol_sp)
            if (allocated(error)) exit
        end do
        call get_err_str(msg, "max deviation: ", err)
        call check_test(error, 'test_qr_factorization_rsp', &
                              & info='Diagonal ordering', &
                              & eq='|R(1,1)| >= |R(2,2)| >= ...', context=msg)

        return
    end subroutine test_qr_factorization_rsp

    subroutine test_qr_rank_deficient_rsp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test Vectors.
        integer, parameter :: kdim = 10
        type(vector_rsp), allocatable :: A(:)
        ! Upper triangular matrix.
        real(sp) :: R(kdim, kdim)
        ! Information flag.
        integer :: info
        ! Column to make collinear.
        integer, parameter :: j_col = 3
        ! Miscellaneous.
        real(sp), allocatable :: Adata(:, :), Qdata(:, :)
        real(sp), allocatable :: G(:, :)
        real(sp) :: err
        character(len=256) :: msg
        real(sp) :: large_tol, beta, alpha
        integer :: i, k, nzero, rk, idx, perm(kdim)
        logical :: mask(kdim)

        ! Use a large tolerance to trigger collinearity detection.
        large_tol = sqrt(epsilon(1.0_sp))

        ! Initialize matrix.
        allocate(A(kdim)) ; call init_rand(A)

        ! Build orthonormal basis: orthogonalize each vector against previous ones.
        call orthonormalize_basis(A)

        ! Make column j_col exactly collinear with column 1: A(j_col) = A(1)
        call copy(A(j_col), A(1))

        ! Save data before QR factorization.
        allocate(Adata(test_size, kdim)) ; call get_data(Adata, A)
        R = zero_rsp

        do
            ! In-place QR factorization.
            call qr(A, R, info, tol=large_tol)

            ! Check correct column has been flagged.
            call check(error, info == j_col)
            if (allocated(error)) exit

            ! Check corresponding entry in R is zero.
            call check(error, abs(R(j_col, j_col)) == 0)
            if (allocated(error)) exit

            ! Get Q data after factorization.
            allocate(Qdata(test_size, kdim)) ; call get_data(Qdata, A)

            ! Check correctness: A = Q @ R must hold even with collinear column.
            err = maxval(abs(Adata - matmul(Qdata, R)))
            call get_err_str(msg, "max err: ", err)
            call check(error, err < rtol_sp)
            exit
        end do
        call check_test(error, 'test_qr_rank_deficient_rsp', &
                              & info='Factorization', eq='A = Q @ R', context=msg)

        ! Compute Gram matrix associated to the Krylov basis.
        G = Gram(A(:kdim))

        ! Check orthonormality of the computed basis.
        err = norm2(abs(G - eye(kdim, mold=1.0_sp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_qr_rank_deficient_rsp', &
                              & info='Basis orthonormality', eq='Q.H @ Q = I', context=msg)

        ! Effective rank.
        rk = kdim - nzero

        ! Initialize matrix.
        call init_rand(A)

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
        call get_data(Adata, A)

        ! In-place QR factorization.
        call qr(A, R, perm, info, tol=atol_sp)
        call check_info(info, 'qr_pivot', module=this_module_long, procedure='test_qr_rank_deficient_rsp')

        ! Extract data
        call get_data(Qdata, A)
        Adata = Adata(:, perm)

        ! Check correctness.
        err = maxval(abs(Adata - matmul(Qdata, R)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_qr_rank_deficient_rsp', &
                              & info='Pivoted factorization', eq='AP = Q @ R', context=msg)

        ! Compute Gram matrix associated to the Krylov basis.
        G = Gram(A(:kdim))

        ! Check orthonormality of the computed basis.
        err = norm2(abs(G - eye(kdim, mold=1.0_sp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_qr_rank_deficient_rsp', &
                              & info='Pivoted basis orthonormality', eq='Q.H @ Q = I', context=msg)

        return
    end subroutine test_qr_rank_deficient_rsp

    subroutine test_qr_invalid_inputs_rsp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        type(vector_rsp), allocatable :: A(:)
        ! Upper triangular matrix.
        real(sp) :: R(5, 5), Rsmall(4, 4)
        ! Permutation vector (too small).
        integer :: perm(5), perm_small(4)
        ! Information flag.
        integer :: info
        character(len=256) :: msg

        ! Test with empty set of vectors.
        test_loop: do
            allocate(A(0))
            call qr(A, R, info) ! Standard QR
            call check(error, info == -1)
            if (allocated(error)) exit test_loop
            call qr(A, R, perm, info)   ! Pivoting QR
            call check(error, info == -1)
            if (allocated(error)) exit test_loop

            ! Test matrix R with inconsistent dimensions.
            deallocate(A) ; allocate(A(5)) ; call init_rand(A)
            call qr(A, Rsmall, info)    ! Standard QR
            call check(error, info==-2)
            if (allocated(error)) exit test_loop
            call qr(A, Rsmall, perm, info)  ! Pivoting QR
            call check(error, info==-2)
            if (allocated(error)) exit test_loop

            ! Test perm too small on qr_with_pivoting.
            call qr(A, R, perm_small, info, tol=atol_sp)
            call check(error, info == -3)
            if (allocated(error)) exit test_loop

            ! Test negative tolerance.
            call qr(A, R, info, tol=-1.0_sp)  ! Standard QR
            call check(error, info == -4)
            if (allocated(error)) exit test_loop
            call qr(A, R, perm, info, tol=-1.0_sp)    ! Pivoting QR
            call check(error, info == -4)
            exit test_loop
        end do test_loop
        call check_test(error, 'test_qr_invalid_inputs_rsp', &
                              & info='Invalid input parameters', eq='', context=msg)

        return
    end subroutine test_qr_invalid_inputs_rsp

    subroutine test_qr_single_vector_rsp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test Vectors.
        integer, parameter :: kdim = 1
        type(vector_rsp), allocatable :: A(:)
        ! Upper triangular matrix.
        real(sp) :: R(kdim, kdim)
        ! Information flag.
        integer :: info, perm(kdim)
        ! Miscellaneous.
        real(sp), allocatable :: Adata(:, :), Qdata(:, :)
        real(sp) :: err
        character(len=256) :: msg

        ! Initialize single vector.
        allocate(A(kdim)) ; call init_rand(A)
        allocate(Adata(test_size, kdim)) ; call get_data(Adata, A)
        R = zero_rsp

        ! In-place QR factorization.
        call qr(A, R, info, tol=atol_sp)
        call check_info(info, 'qr', module=this_module_long, &
            & procedure='test_qr_single_vector_rsp')

        ! Get Q data.
        allocate(Qdata(test_size, kdim)) ; call get_data(Qdata, A)

        ! Check correctness.
        err = maxval(abs(Adata - matmul(Qdata, R)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)

        if (.not. allocated(error)) then    ! Pivoting QR.
            ! Initialize single vector.
            call init_rand(A) ; call get_data(Adata, A)
            R = zero_rsp

            ! In-place QR factorization.
            call qr(A, R, perm, info, tol=atol_sp)
            call check_info(info, 'qr', module=this_module_long, &
                & procedure='test_qr_single_vector_rsp')

            ! Check perm(1) = 1.
            call check(error, perm(1) == 1)

            ! Get Q data.
            call get_data(Qdata, A)

            ! Check correctness.
            err = maxval(abs(Adata - matmul(Qdata, R)))
            call get_err_str(msg, "max err: ", err)
            call check(error, err < rtol_sp)
        endif
        call check_test(error, 'test_qr_single_vector_rsp', &
                              & info='Factorization', eq='A = Q @ R', context=msg)

        return
    end subroutine test_qr_single_vector_rsp

    subroutine collect_qr_rdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                        new_unittest("QR invalid inputs", test_qr_invalid_inputs_rdp), &
                        new_unittest("QR factorization", test_qr_factorization_rdp), &
                        new_unittest("QR rank deficient", test_qr_rank_deficient_rdp), &
                        new_unittest("QR single vector", test_qr_single_vector_rdp) &
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
        real(dp) :: R(kdim, kdim)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(dp), allocatable :: Adata(:, :), Qdata(:, :)
        real(dp), allocatable :: G(:, :)
        real(dp) :: err
        integer :: perm(kdim), i
        character(len=256) :: msg

        ! Initialiaze matrix.
        allocate(A(kdim)) ; call init_rand(A); R = zero_rdp

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


        call init_rand(A); R = zero_rdp ; call get_data(Adata, A)

        ! In-place Pivoted QR factorization.
        call qr(A, R, perm, info, tol=atol_dp)
        call check_info(info, 'qr_pivot', module=this_module_long, procedure='test_qr_factorization_rdp')

        ! Get data.
        call get_data(Qdata, A)

        ! Check correctness.
        err = maxval(abs(Adata(:, perm) - matmul(Qdata, R)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_qr_factorization_rdp', &
                              & info='Pivoted factorization', eq='AP = Q @ R', context=msg)

        ! Compute Gram matrix associated to the Krylov basis.
        G = Gram(A(:kdim))

        ! Check orthonormality of the computed basis.
        err = norm2(abs(G - eye(kdim, mold=1.0_dp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_qr_factorization_rdp', &
                              & info='Pivoted basis orthonormality', eq='Q.H @ Q = I', context=msg)

        ! Check that diagonal entries are non-increasing in magnitude.
        ! |R(1,1)| >= |R(2,2)| >= ... >= |R(kdim,kdim)|
        do i = 1, kdim-1
            err = abs(R(i+1, i+1)) - abs(R(i, i))
            call check(error, err <= rtol_dp)
            if (allocated(error)) exit
        end do
        call get_err_str(msg, "max deviation: ", err)
        call check_test(error, 'test_qr_factorization_rdp', &
                              & info='Diagonal ordering', &
                              & eq='|R(1,1)| >= |R(2,2)| >= ...', context=msg)

        return
    end subroutine test_qr_factorization_rdp

    subroutine test_qr_rank_deficient_rdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test Vectors.
        integer, parameter :: kdim = 10
        type(vector_rdp), allocatable :: A(:)
        ! Upper triangular matrix.
        real(dp) :: R(kdim, kdim)
        ! Information flag.
        integer :: info
        ! Column to make collinear.
        integer, parameter :: j_col = 3
        ! Miscellaneous.
        real(dp), allocatable :: Adata(:, :), Qdata(:, :)
        real(dp), allocatable :: G(:, :)
        real(dp) :: err
        character(len=256) :: msg
        real(dp) :: large_tol, beta, alpha
        integer :: i, k, nzero, rk, idx, perm(kdim)
        logical :: mask(kdim)

        ! Use a large tolerance to trigger collinearity detection.
        large_tol = sqrt(epsilon(1.0_dp))

        ! Initialize matrix.
        allocate(A(kdim)) ; call init_rand(A)

        ! Build orthonormal basis: orthogonalize each vector against previous ones.
        call orthonormalize_basis(A)

        ! Make column j_col exactly collinear with column 1: A(j_col) = A(1)
        call copy(A(j_col), A(1))

        ! Save data before QR factorization.
        allocate(Adata(test_size, kdim)) ; call get_data(Adata, A)
        R = zero_rdp

        do
            ! In-place QR factorization.
            call qr(A, R, info, tol=large_tol)

            ! Check correct column has been flagged.
            call check(error, info == j_col)
            if (allocated(error)) exit

            ! Check corresponding entry in R is zero.
            call check(error, abs(R(j_col, j_col)) == 0)
            if (allocated(error)) exit

            ! Get Q data after factorization.
            allocate(Qdata(test_size, kdim)) ; call get_data(Qdata, A)

            ! Check correctness: A = Q @ R must hold even with collinear column.
            err = maxval(abs(Adata - matmul(Qdata, R)))
            call get_err_str(msg, "max err: ", err)
            call check(error, err < rtol_dp)
            exit
        end do
        call check_test(error, 'test_qr_rank_deficient_rdp', &
                              & info='Factorization', eq='A = Q @ R', context=msg)

        ! Compute Gram matrix associated to the Krylov basis.
        G = Gram(A(:kdim))

        ! Check orthonormality of the computed basis.
        err = norm2(abs(G - eye(kdim, mold=1.0_dp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_qr_rank_deficient_rdp', &
                              & info='Basis orthonormality', eq='Q.H @ Q = I', context=msg)

        ! Effective rank.
        rk = kdim - nzero

        ! Initialize matrix.
        call init_rand(A)

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
        call get_data(Adata, A)

        ! In-place QR factorization.
        call qr(A, R, perm, info, tol=atol_dp)
        call check_info(info, 'qr_pivot', module=this_module_long, procedure='test_qr_rank_deficient_rdp')

        ! Extract data
        call get_data(Qdata, A)
        Adata = Adata(:, perm)

        ! Check correctness.
        err = maxval(abs(Adata - matmul(Qdata, R)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_qr_rank_deficient_rdp', &
                              & info='Pivoted factorization', eq='AP = Q @ R', context=msg)

        ! Compute Gram matrix associated to the Krylov basis.
        G = Gram(A(:kdim))

        ! Check orthonormality of the computed basis.
        err = norm2(abs(G - eye(kdim, mold=1.0_dp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_qr_rank_deficient_rdp', &
                              & info='Pivoted basis orthonormality', eq='Q.H @ Q = I', context=msg)

        return
    end subroutine test_qr_rank_deficient_rdp

    subroutine test_qr_invalid_inputs_rdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        type(vector_rdp), allocatable :: A(:)
        ! Upper triangular matrix.
        real(dp) :: R(5, 5), Rsmall(4, 4)
        ! Permutation vector (too small).
        integer :: perm(5), perm_small(4)
        ! Information flag.
        integer :: info
        character(len=256) :: msg

        ! Test with empty set of vectors.
        test_loop: do
            allocate(A(0))
            call qr(A, R, info) ! Standard QR
            call check(error, info == -1)
            if (allocated(error)) exit test_loop
            call qr(A, R, perm, info)   ! Pivoting QR
            call check(error, info == -1)
            if (allocated(error)) exit test_loop

            ! Test matrix R with inconsistent dimensions.
            deallocate(A) ; allocate(A(5)) ; call init_rand(A)
            call qr(A, Rsmall, info)    ! Standard QR
            call check(error, info==-2)
            if (allocated(error)) exit test_loop
            call qr(A, Rsmall, perm, info)  ! Pivoting QR
            call check(error, info==-2)
            if (allocated(error)) exit test_loop

            ! Test perm too small on qr_with_pivoting.
            call qr(A, R, perm_small, info, tol=atol_dp)
            call check(error, info == -3)
            if (allocated(error)) exit test_loop

            ! Test negative tolerance.
            call qr(A, R, info, tol=-1.0_dp)  ! Standard QR
            call check(error, info == -4)
            if (allocated(error)) exit test_loop
            call qr(A, R, perm, info, tol=-1.0_dp)    ! Pivoting QR
            call check(error, info == -4)
            exit test_loop
        end do test_loop
        call check_test(error, 'test_qr_invalid_inputs_rdp', &
                              & info='Invalid input parameters', eq='', context=msg)

        return
    end subroutine test_qr_invalid_inputs_rdp

    subroutine test_qr_single_vector_rdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test Vectors.
        integer, parameter :: kdim = 1
        type(vector_rdp), allocatable :: A(:)
        ! Upper triangular matrix.
        real(dp) :: R(kdim, kdim)
        ! Information flag.
        integer :: info, perm(kdim)
        ! Miscellaneous.
        real(dp), allocatable :: Adata(:, :), Qdata(:, :)
        real(dp) :: err
        character(len=256) :: msg

        ! Initialize single vector.
        allocate(A(kdim)) ; call init_rand(A)
        allocate(Adata(test_size, kdim)) ; call get_data(Adata, A)
        R = zero_rdp

        ! In-place QR factorization.
        call qr(A, R, info, tol=atol_dp)
        call check_info(info, 'qr', module=this_module_long, &
            & procedure='test_qr_single_vector_rdp')

        ! Get Q data.
        allocate(Qdata(test_size, kdim)) ; call get_data(Qdata, A)

        ! Check correctness.
        err = maxval(abs(Adata - matmul(Qdata, R)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)

        if (.not. allocated(error)) then    ! Pivoting QR.
            ! Initialize single vector.
            call init_rand(A) ; call get_data(Adata, A)
            R = zero_rdp

            ! In-place QR factorization.
            call qr(A, R, perm, info, tol=atol_dp)
            call check_info(info, 'qr', module=this_module_long, &
                & procedure='test_qr_single_vector_rdp')

            ! Check perm(1) = 1.
            call check(error, perm(1) == 1)

            ! Get Q data.
            call get_data(Qdata, A)

            ! Check correctness.
            err = maxval(abs(Adata - matmul(Qdata, R)))
            call get_err_str(msg, "max err: ", err)
            call check(error, err < rtol_dp)
        endif
        call check_test(error, 'test_qr_single_vector_rdp', &
                              & info='Factorization', eq='A = Q @ R', context=msg)

        return
    end subroutine test_qr_single_vector_rdp

    subroutine collect_qr_csp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                        new_unittest("QR invalid inputs", test_qr_invalid_inputs_csp), &
                        new_unittest("QR factorization", test_qr_factorization_csp), &
                        new_unittest("QR rank deficient", test_qr_rank_deficient_csp), &
                        new_unittest("QR single vector", test_qr_single_vector_csp) &
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
        complex(sp) :: R(kdim, kdim)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        complex(sp), allocatable :: Adata(:, :), Qdata(:, :)
        complex(sp), allocatable :: G(:, :)
        real(sp) :: err
        integer :: perm(kdim), i
        character(len=256) :: msg

        ! Initialiaze matrix.
        allocate(A(kdim)) ; call init_rand(A); R = zero_csp

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


        call init_rand(A); R = zero_csp ; call get_data(Adata, A)

        ! In-place Pivoted QR factorization.
        call qr(A, R, perm, info, tol=atol_sp)
        call check_info(info, 'qr_pivot', module=this_module_long, procedure='test_qr_factorization_csp')

        ! Get data.
        call get_data(Qdata, A)

        ! Check correctness.
        err = maxval(abs(Adata(:, perm) - matmul(Qdata, R)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_qr_factorization_csp', &
                              & info='Pivoted factorization', eq='AP = Q @ R', context=msg)

        ! Compute Gram matrix associated to the Krylov basis.
        G = Gram(A(:kdim))

        ! Check orthonormality of the computed basis.
        err = norm2(abs(G - eye(kdim, mold=1.0_sp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_qr_factorization_csp', &
                              & info='Pivoted basis orthonormality', eq='Q.H @ Q = I', context=msg)

        ! Check that diagonal entries are non-increasing in magnitude.
        ! |R(1,1)| >= |R(2,2)| >= ... >= |R(kdim,kdim)|
        do i = 1, kdim-1
            err = abs(R(i+1, i+1)) - abs(R(i, i))
            call check(error, err <= rtol_sp)
            if (allocated(error)) exit
        end do
        call get_err_str(msg, "max deviation: ", err)
        call check_test(error, 'test_qr_factorization_csp', &
                              & info='Diagonal ordering', &
                              & eq='|R(1,1)| >= |R(2,2)| >= ...', context=msg)

        return
    end subroutine test_qr_factorization_csp

    subroutine test_qr_rank_deficient_csp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test Vectors.
        integer, parameter :: kdim = 10
        type(vector_csp), allocatable :: A(:)
        ! Upper triangular matrix.
        complex(sp) :: R(kdim, kdim)
        ! Information flag.
        integer :: info
        ! Column to make collinear.
        integer, parameter :: j_col = 3
        ! Miscellaneous.
        complex(sp), allocatable :: Adata(:, :), Qdata(:, :)
        complex(sp), allocatable :: G(:, :)
        real(sp) :: err
        character(len=256) :: msg
        real(sp) :: large_tol, beta, alpha
        integer :: i, k, nzero, rk, idx, perm(kdim)
        logical :: mask(kdim)

        ! Use a large tolerance to trigger collinearity detection.
        large_tol = sqrt(epsilon(1.0_sp))

        ! Initialize matrix.
        allocate(A(kdim)) ; call init_rand(A)

        ! Build orthonormal basis: orthogonalize each vector against previous ones.
        call orthonormalize_basis(A)

        ! Make column j_col exactly collinear with column 1: A(j_col) = A(1)
        call copy(A(j_col), A(1))

        ! Save data before QR factorization.
        allocate(Adata(test_size, kdim)) ; call get_data(Adata, A)
        R = zero_csp

        do
            ! In-place QR factorization.
            call qr(A, R, info, tol=large_tol)

            ! Check correct column has been flagged.
            call check(error, info == j_col)
            if (allocated(error)) exit

            ! Check corresponding entry in R is zero.
            call check(error, abs(R(j_col, j_col)) == 0)
            if (allocated(error)) exit

            ! Get Q data after factorization.
            allocate(Qdata(test_size, kdim)) ; call get_data(Qdata, A)

            ! Check correctness: A = Q @ R must hold even with collinear column.
            err = maxval(abs(Adata - matmul(Qdata, R)))
            call get_err_str(msg, "max err: ", err)
            call check(error, err < rtol_sp)
            exit
        end do
        call check_test(error, 'test_qr_rank_deficient_csp', &
                              & info='Factorization', eq='A = Q @ R', context=msg)

        ! Compute Gram matrix associated to the Krylov basis.
        G = Gram(A(:kdim))

        ! Check orthonormality of the computed basis.
        err = norm2(abs(G - eye(kdim, mold=1.0_sp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_qr_rank_deficient_csp', &
                              & info='Basis orthonormality', eq='Q.H @ Q = I', context=msg)

        ! Effective rank.
        rk = kdim - nzero

        ! Initialize matrix.
        call init_rand(A)

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
        call get_data(Adata, A)

        ! In-place QR factorization.
        call qr(A, R, perm, info, tol=atol_sp)
        call check_info(info, 'qr_pivot', module=this_module_long, procedure='test_qr_rank_deficient_csp')

        ! Extract data
        call get_data(Qdata, A)
        Adata = Adata(:, perm)

        ! Check correctness.
        err = maxval(abs(Adata - matmul(Qdata, R)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_qr_rank_deficient_csp', &
                              & info='Pivoted factorization', eq='AP = Q @ R', context=msg)

        ! Compute Gram matrix associated to the Krylov basis.
        G = Gram(A(:kdim))

        ! Check orthonormality of the computed basis.
        err = norm2(abs(G - eye(kdim, mold=1.0_sp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_qr_rank_deficient_csp', &
                              & info='Pivoted basis orthonormality', eq='Q.H @ Q = I', context=msg)

        return
    end subroutine test_qr_rank_deficient_csp

    subroutine test_qr_invalid_inputs_csp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        type(vector_csp), allocatable :: A(:)
        ! Upper triangular matrix.
        complex(sp) :: R(5, 5), Rsmall(4, 4)
        ! Permutation vector (too small).
        integer :: perm(5), perm_small(4)
        ! Information flag.
        integer :: info
        character(len=256) :: msg

        ! Test with empty set of vectors.
        test_loop: do
            allocate(A(0))
            call qr(A, R, info) ! Standard QR
            call check(error, info == -1)
            if (allocated(error)) exit test_loop
            call qr(A, R, perm, info)   ! Pivoting QR
            call check(error, info == -1)
            if (allocated(error)) exit test_loop

            ! Test matrix R with inconsistent dimensions.
            deallocate(A) ; allocate(A(5)) ; call init_rand(A)
            call qr(A, Rsmall, info)    ! Standard QR
            call check(error, info==-2)
            if (allocated(error)) exit test_loop
            call qr(A, Rsmall, perm, info)  ! Pivoting QR
            call check(error, info==-2)
            if (allocated(error)) exit test_loop

            ! Test perm too small on qr_with_pivoting.
            call qr(A, R, perm_small, info, tol=atol_sp)
            call check(error, info == -3)
            if (allocated(error)) exit test_loop

            ! Test negative tolerance.
            call qr(A, R, info, tol=-1.0_sp)  ! Standard QR
            call check(error, info == -4)
            if (allocated(error)) exit test_loop
            call qr(A, R, perm, info, tol=-1.0_sp)    ! Pivoting QR
            call check(error, info == -4)
            exit test_loop
        end do test_loop
        call check_test(error, 'test_qr_invalid_inputs_csp', &
                              & info='Invalid input parameters', eq='', context=msg)

        return
    end subroutine test_qr_invalid_inputs_csp

    subroutine test_qr_single_vector_csp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test Vectors.
        integer, parameter :: kdim = 1
        type(vector_csp), allocatable :: A(:)
        ! Upper triangular matrix.
        complex(sp) :: R(kdim, kdim)
        ! Information flag.
        integer :: info, perm(kdim)
        ! Miscellaneous.
        complex(sp), allocatable :: Adata(:, :), Qdata(:, :)
        real(sp) :: err
        character(len=256) :: msg

        ! Initialize single vector.
        allocate(A(kdim)) ; call init_rand(A)
        allocate(Adata(test_size, kdim)) ; call get_data(Adata, A)
        R = zero_csp

        ! In-place QR factorization.
        call qr(A, R, info, tol=atol_sp)
        call check_info(info, 'qr', module=this_module_long, &
            & procedure='test_qr_single_vector_csp')

        ! Get Q data.
        allocate(Qdata(test_size, kdim)) ; call get_data(Qdata, A)

        ! Check correctness.
        err = maxval(abs(Adata - matmul(Qdata, R)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)

        if (.not. allocated(error)) then    ! Pivoting QR.
            ! Initialize single vector.
            call init_rand(A) ; call get_data(Adata, A)
            R = zero_csp

            ! In-place QR factorization.
            call qr(A, R, perm, info, tol=atol_sp)
            call check_info(info, 'qr', module=this_module_long, &
                & procedure='test_qr_single_vector_csp')

            ! Check perm(1) = 1.
            call check(error, perm(1) == 1)

            ! Get Q data.
            call get_data(Qdata, A)

            ! Check correctness.
            err = maxval(abs(Adata - matmul(Qdata, R)))
            call get_err_str(msg, "max err: ", err)
            call check(error, err < rtol_sp)
        endif
        call check_test(error, 'test_qr_single_vector_csp', &
                              & info='Factorization', eq='A = Q @ R', context=msg)

        return
    end subroutine test_qr_single_vector_csp

    subroutine collect_qr_cdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                        new_unittest("QR invalid inputs", test_qr_invalid_inputs_cdp), &
                        new_unittest("QR factorization", test_qr_factorization_cdp), &
                        new_unittest("QR rank deficient", test_qr_rank_deficient_cdp), &
                        new_unittest("QR single vector", test_qr_single_vector_cdp) &
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
        complex(dp) :: R(kdim, kdim)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        complex(dp), allocatable :: Adata(:, :), Qdata(:, :)
        complex(dp), allocatable :: G(:, :)
        real(dp) :: err
        integer :: perm(kdim), i
        character(len=256) :: msg

        ! Initialiaze matrix.
        allocate(A(kdim)) ; call init_rand(A); R = zero_cdp

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


        call init_rand(A); R = zero_cdp ; call get_data(Adata, A)

        ! In-place Pivoted QR factorization.
        call qr(A, R, perm, info, tol=atol_dp)
        call check_info(info, 'qr_pivot', module=this_module_long, procedure='test_qr_factorization_cdp')

        ! Get data.
        call get_data(Qdata, A)

        ! Check correctness.
        err = maxval(abs(Adata(:, perm) - matmul(Qdata, R)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_qr_factorization_cdp', &
                              & info='Pivoted factorization', eq='AP = Q @ R', context=msg)

        ! Compute Gram matrix associated to the Krylov basis.
        G = Gram(A(:kdim))

        ! Check orthonormality of the computed basis.
        err = norm2(abs(G - eye(kdim, mold=1.0_dp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_qr_factorization_cdp', &
                              & info='Pivoted basis orthonormality', eq='Q.H @ Q = I', context=msg)

        ! Check that diagonal entries are non-increasing in magnitude.
        ! |R(1,1)| >= |R(2,2)| >= ... >= |R(kdim,kdim)|
        do i = 1, kdim-1
            err = abs(R(i+1, i+1)) - abs(R(i, i))
            call check(error, err <= rtol_dp)
            if (allocated(error)) exit
        end do
        call get_err_str(msg, "max deviation: ", err)
        call check_test(error, 'test_qr_factorization_cdp', &
                              & info='Diagonal ordering', &
                              & eq='|R(1,1)| >= |R(2,2)| >= ...', context=msg)

        return
    end subroutine test_qr_factorization_cdp

    subroutine test_qr_rank_deficient_cdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test Vectors.
        integer, parameter :: kdim = 10
        type(vector_cdp), allocatable :: A(:)
        ! Upper triangular matrix.
        complex(dp) :: R(kdim, kdim)
        ! Information flag.
        integer :: info
        ! Column to make collinear.
        integer, parameter :: j_col = 3
        ! Miscellaneous.
        complex(dp), allocatable :: Adata(:, :), Qdata(:, :)
        complex(dp), allocatable :: G(:, :)
        real(dp) :: err
        character(len=256) :: msg
        real(dp) :: large_tol, beta, alpha
        integer :: i, k, nzero, rk, idx, perm(kdim)
        logical :: mask(kdim)

        ! Use a large tolerance to trigger collinearity detection.
        large_tol = sqrt(epsilon(1.0_dp))

        ! Initialize matrix.
        allocate(A(kdim)) ; call init_rand(A)

        ! Build orthonormal basis: orthogonalize each vector against previous ones.
        call orthonormalize_basis(A)

        ! Make column j_col exactly collinear with column 1: A(j_col) = A(1)
        call copy(A(j_col), A(1))

        ! Save data before QR factorization.
        allocate(Adata(test_size, kdim)) ; call get_data(Adata, A)
        R = zero_cdp

        do
            ! In-place QR factorization.
            call qr(A, R, info, tol=large_tol)

            ! Check correct column has been flagged.
            call check(error, info == j_col)
            if (allocated(error)) exit

            ! Check corresponding entry in R is zero.
            call check(error, abs(R(j_col, j_col)) == 0)
            if (allocated(error)) exit

            ! Get Q data after factorization.
            allocate(Qdata(test_size, kdim)) ; call get_data(Qdata, A)

            ! Check correctness: A = Q @ R must hold even with collinear column.
            err = maxval(abs(Adata - matmul(Qdata, R)))
            call get_err_str(msg, "max err: ", err)
            call check(error, err < rtol_dp)
            exit
        end do
        call check_test(error, 'test_qr_rank_deficient_cdp', &
                              & info='Factorization', eq='A = Q @ R', context=msg)

        ! Compute Gram matrix associated to the Krylov basis.
        G = Gram(A(:kdim))

        ! Check orthonormality of the computed basis.
        err = norm2(abs(G - eye(kdim, mold=1.0_dp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_qr_rank_deficient_cdp', &
                              & info='Basis orthonormality', eq='Q.H @ Q = I', context=msg)

        ! Effective rank.
        rk = kdim - nzero

        ! Initialize matrix.
        call init_rand(A)

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
        call get_data(Adata, A)

        ! In-place QR factorization.
        call qr(A, R, perm, info, tol=atol_dp)
        call check_info(info, 'qr_pivot', module=this_module_long, procedure='test_qr_rank_deficient_cdp')

        ! Extract data
        call get_data(Qdata, A)
        Adata = Adata(:, perm)

        ! Check correctness.
        err = maxval(abs(Adata - matmul(Qdata, R)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_qr_rank_deficient_cdp', &
                              & info='Pivoted factorization', eq='AP = Q @ R', context=msg)

        ! Compute Gram matrix associated to the Krylov basis.
        G = Gram(A(:kdim))

        ! Check orthonormality of the computed basis.
        err = norm2(abs(G - eye(kdim, mold=1.0_dp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_qr_rank_deficient_cdp', &
                              & info='Pivoted basis orthonormality', eq='Q.H @ Q = I', context=msg)

        return
    end subroutine test_qr_rank_deficient_cdp

    subroutine test_qr_invalid_inputs_cdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        type(vector_cdp), allocatable :: A(:)
        ! Upper triangular matrix.
        complex(dp) :: R(5, 5), Rsmall(4, 4)
        ! Permutation vector (too small).
        integer :: perm(5), perm_small(4)
        ! Information flag.
        integer :: info
        character(len=256) :: msg

        ! Test with empty set of vectors.
        test_loop: do
            allocate(A(0))
            call qr(A, R, info) ! Standard QR
            call check(error, info == -1)
            if (allocated(error)) exit test_loop
            call qr(A, R, perm, info)   ! Pivoting QR
            call check(error, info == -1)
            if (allocated(error)) exit test_loop

            ! Test matrix R with inconsistent dimensions.
            deallocate(A) ; allocate(A(5)) ; call init_rand(A)
            call qr(A, Rsmall, info)    ! Standard QR
            call check(error, info==-2)
            if (allocated(error)) exit test_loop
            call qr(A, Rsmall, perm, info)  ! Pivoting QR
            call check(error, info==-2)
            if (allocated(error)) exit test_loop

            ! Test perm too small on qr_with_pivoting.
            call qr(A, R, perm_small, info, tol=atol_dp)
            call check(error, info == -3)
            if (allocated(error)) exit test_loop

            ! Test negative tolerance.
            call qr(A, R, info, tol=-1.0_dp)  ! Standard QR
            call check(error, info == -4)
            if (allocated(error)) exit test_loop
            call qr(A, R, perm, info, tol=-1.0_dp)    ! Pivoting QR
            call check(error, info == -4)
            exit test_loop
        end do test_loop
        call check_test(error, 'test_qr_invalid_inputs_cdp', &
                              & info='Invalid input parameters', eq='', context=msg)

        return
    end subroutine test_qr_invalid_inputs_cdp

    subroutine test_qr_single_vector_cdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test Vectors.
        integer, parameter :: kdim = 1
        type(vector_cdp), allocatable :: A(:)
        ! Upper triangular matrix.
        complex(dp) :: R(kdim, kdim)
        ! Information flag.
        integer :: info, perm(kdim)
        ! Miscellaneous.
        complex(dp), allocatable :: Adata(:, :), Qdata(:, :)
        real(dp) :: err
        character(len=256) :: msg

        ! Initialize single vector.
        allocate(A(kdim)) ; call init_rand(A)
        allocate(Adata(test_size, kdim)) ; call get_data(Adata, A)
        R = zero_cdp

        ! In-place QR factorization.
        call qr(A, R, info, tol=atol_dp)
        call check_info(info, 'qr', module=this_module_long, &
            & procedure='test_qr_single_vector_cdp')

        ! Get Q data.
        allocate(Qdata(test_size, kdim)) ; call get_data(Qdata, A)

        ! Check correctness.
        err = maxval(abs(Adata - matmul(Qdata, R)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)

        if (.not. allocated(error)) then    ! Pivoting QR.
            ! Initialize single vector.
            call init_rand(A) ; call get_data(Adata, A)
            R = zero_cdp

            ! In-place QR factorization.
            call qr(A, R, perm, info, tol=atol_dp)
            call check_info(info, 'qr', module=this_module_long, &
                & procedure='test_qr_single_vector_cdp')

            ! Check perm(1) = 1.
            call check(error, perm(1) == 1)

            ! Get Q data.
            call get_data(Qdata, A)

            ! Check correctness.
            err = maxval(abs(Adata - matmul(Qdata, R)))
            call get_err_str(msg, "max err: ", err)
            call check(error, err < rtol_dp)
        endif
        call check_test(error, 'test_qr_single_vector_cdp', &
                              & info='Factorization', eq='A = Q @ R', context=msg)

        return
    end subroutine test_qr_single_vector_cdp


    !--------------------------------------------------------------
    !-----     DEFINITIONS OF THE UNIT-TESTS FOR ARNOLDI      -----
    !--------------------------------------------------------------

    subroutine collect_arnoldi_rsp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("Arnoldi invalid parameters", test_arnoldi_invalid_params_rsp), &
            new_unittest("Arnoldi factorization", test_arnoldi_factorization_rsp), &
            new_unittest("Arnoldi shifted matrix", test_arnoldi_shifted_matrix_rsp), &
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

        ! Check H is indeed Hessenberg.
        call check(error, is_hessenberg(H, uplo='u'))
        call check_test(error, 'test_arnoldi_factorization_rsp', &
                               info='Upper Hessenberg', eq='Is H Hessenberg?', context=msg)

        ! Check correctness of full factorization.
        allocate(Xdata(test_size, kdim+1)) ; call get_data(Xdata, X)
        err = maxval(abs(matmul(A%data, Xdata(:, :kdim)) - matmul(Xdata, H)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_arnoldi_factorization_rsp', &
                              & info='Factorization', eq='A @ X = X_ @ H_', context=msg)

        ! Compute Gram matrix associated to the Krylov basis.
        G = Gram(X(:kdim))

        ! Check orthonormality of the computed basis.
        err = maxval(abs(G - eye(kdim, mold=1.0_sp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_arnoldi_factorization_rsp', &
                              & info='Orthonomality', eq='X.H @ X = I', context=msg)

        block
        ! Krylov subspaces.
        type(vector_rsp), allocatable :: Xfull(:), Xrestart(:)
        integer, parameter :: kstart = kdim/2
        ! Hessenberg matrix.
        real(sp), allocatable :: Hfull(:, :), Hrestart(:, :)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(sp), allocatable :: Xfull_data(:, :), Xrestart_data(:, :)

        ! Initialize Krylov subspace.
        allocate(Xfull(kdim+1)); call zero_basis(Xfull); call Xfull(1)%rand(ifnorm=.true.)
        allocate(Hfull(kdim+1, kdim), Hrestart(kdim+1, kdim), source=zero_rsp)

        ! Full Arnoldi factorization.
        call arnoldi(A, Xfull, Hfull, info, tol=atol_sp)
        call check_info(info, 'arnoldi', module=this_module_long, procedure='test_arnoldi_factorization_rsp')

        ! Copy data for restart.
        allocate(Xrestart(kdim+1)) ; call zero_basis(Xrestart)
        call copy(Xrestart(:kstart), Xfull(:kstart))
        Hrestart(:kstart, :kstart-1) = Hfull(:kstart, :kstart-1)

        ! Restart Arnoldi factorization.
        call arnoldi(A, Xrestart, Hrestart, info, kstart=kstart, tol=atol_sp)

        ! Compute inner product between the two bases.
        G = innerprod(Xfull(:kdim), Xrestart(:kdim))
        err = maxval(abs(G - eye(kdim, mold=1.0_sp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_arnoldi_factorization_rsp', &
                                info='Restart', eq='Xfull = Xrestart', context=msg)

        ! Check Hessenberg matrices.
        err = maxval(abs(Hfull - Hrestart))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_arnoldi_factorization_rsp', &
                                info='Restart', eq='Hfull = Hrestart', context=msg)

        end block

        block
        integer :: k_inv, k_max

        k_inv = 10
        k_max = test_size
        ! Create block lower triangular matrix.
        call init_rand(A)
        A%data(k_inv+1:, :k_inv) = zero_rsp

        ! Initial vector v1 must be in the invariant subspace
        deallocate(X) ; allocate(X(k_max+1))
        call zero_basis(X)
        call X(1)%rand(ifnorm=.false.)
        X(1)%data(k_inv+1:) = zero_rsp
        err = X(1)%norm()
        call X(1)%scal(one_rsp/err)

        deallocate(H) ; allocate(H(k_max+1, k_max), source=zero_rsp)
        call arnoldi(A, X, H, info, tol=atol_sp)

        ! 1. Check if Arnoldi detected the invariant subspace dimension
        call check(error, info == k_inv)
        call check_test(error, 'test_arnoldi_factorization_rsp', &
                              & info='Subspace Dim', eq='info == k_inv', context='Invariant detection')

        ! 2. Check AX_k = X_k H_k
        deallocate(Xdata) ; allocate(Xdata(test_size, info))
        call get_data(Xdata, X(:info))

        err = maxval(abs(matmul(A%data, Xdata) - matmul(Xdata, H(1:info, 1:info))))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_arnoldi_factorization_rsp', &
                              & info='Invariant Property', eq='A @ X = X @ H', context=msg)

        end block

        return
    end subroutine test_arnoldi_factorization_rsp

    subroutine test_arnoldi_invalid_params_rsp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test linear operator.
        type(linop_rsp), allocatable :: A
        ! Krylov subspace.
        type(vector_rsp), allocatable :: X(:)
        ! Hessenberg matrix.
        real(sp), allocatable :: H(:, :)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        integer, parameter :: kdim = test_size
        integer :: p

        ! Common initialization for each case.
        do
            A = linop_rsp() ; call init_rand(A)
            allocate(H(kdim+1, kdim)) ; H = zero_rsp
            p = 1

            ! --- Case 1: k_start < 1 or k_start > k_end (info = -5) ---
            allocate(X(kdim+1)) ; call zero_basis(X)
            call arnoldi(A, X, H, info, kstart=0, tol=atol_sp)
            call check(error, info == -5)
            if (allocated(error)) exit
            call arnoldi(A, X, H, info, kstart=kdim, kend=1, tol=atol_sp)
            call check(error, info == -5)
            if (allocated(error)) exit
            ! --- Case 2: k_end > kdim (info = -6) ---
            call arnoldi(A, X, H, info, kend=kdim+1, tol=atol_sp)
            call check(error, info == -6)
            if (allocated(error)) exit
            ! --- Case 3: tolerance < 0 (info = -7) ---
            call arnoldi(A, X, H, info, tol=-1.0_sp)
            call check(error, info == -7)
            if (allocated(error)) exit
            ! --- Case 4: mod(size(X), p) /= 0 (info = -9) ---
            ! If kdim=20, size(X)=21. p=2 does not divide 21.
            p = 2
            call arnoldi(A, X, H, info, blksize=p, tol=atol_sp)
            call check(error, info == -9)
            if (allocated(error)) exit
            ! --- Case 5: blksize <= 0 (info = -9) ---
            p = 0
            call arnoldi(A, X, H, info, blksize=p, tol=atol_sp)
            call check(error, info == -9)
            if (allocated(error)) exit
            ! --- Case 6: H has wrong first dimension (info = -3) ---
            p = 1
            deallocate(H)
            allocate(H(kdim, kdim)) ; H = zero_rsp
            call arnoldi(A, X, H, info, blksize=p, tol=atol_sp)
            call check(error, info == -3)
            if (allocated(error)) exit
            ! --- Case 7: H has wrong second dimension (info = -3) ---
            deallocate(H)
            allocate(H(kdim+1, kdim-1)) ; H = zero_rsp
            call arnoldi(A, X, H, info, blksize=p, tol=atol_sp)
            call check(error, info == -3)
            if (allocated(error)) exit
            ! --- Case 8: block Arnoldi with too-small H (info = -3) ---
            block
                integer, parameter :: p_block = 2
                integer, parameter :: kdim_block = test_size/2
                type(vector_rsp), allocatable :: X0(:)

                deallocate(X, H)
                allocate(X(p_block*(kdim_block+1))) ; allocate(X0(p_block))
                call init_rand(X0) ; call initialize_krylov_subspace(X, X0)
                allocate(H(p_block*(kdim_block+1) - 1, p_block*kdim_block)) ; H = zero_rsp
                call arnoldi(A, X, H, info, blksize=p_block, tol=atol_sp)
                call check(error, info == -3)
           end block
            exit
        enddo
        call check_test(error, 'test_arnoldi_invalid_params_rsp', &
                          & info='Invalid parameters', eq='', context='block p=2')

        return
    end subroutine test_arnoldi_invalid_params_rsp

    subroutine test_arnoldi_shifted_matrix_rsp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Linear operators.
        type(linop_rsp), allocatable :: A
        type(axpby_linop_rsp), allocatable :: A_shift
        ! Krylov subspaces.
        type(vector_rsp), allocatable :: X(:), X_shift(:)
        ! Hessenberg matrices.
        real(sp), allocatable :: H(:, :), H_shift(:, :)
        ! Information flags.
        integer :: info, info_shift, i
        ! Miscellaneous.
        integer, parameter :: kdim = test_size
        real(sp) :: sigma
        real(sp), allocatable :: G_basis(:, :)
        real(sp), allocatable :: H_diff(:, :)
        real(sp) :: err
        character(len=256) :: msg

        ! 1. Initialize linear operator A.
        A = linop_rsp() ; call init_rand(A)

        ! 2. Arnoldi factorization for A.
        allocate(X(kdim+1)); call zero_basis(X); call X(1)%rand(ifnorm = .true.)
        allocate(H(kdim+1, kdim)) ; H = zero_rsp
        call arnoldi(A, X, H, info, tol=atol_sp)
        call check_info(info, 'arnoldi', module=this_module_long, procedure='test_arnoldi_shifted_matrix_rsp')

        ! 3. Construct shifted operator A_shift = A + sigma*I.
        call random_number(sigma)
        allocate(A_shift)
        A_shift%A = A
        A_shift%B = Id_rsp()
        A_shift%alpha = 1.0_sp
        A_shift%beta = sigma
        A_shift%transA = .false.
        A_shift%transB = .false.

        ! 4. Arnoldi factorization for A_shift using the same starting vector.
        allocate(X_shift(kdim+1)); call zero_basis(X_shift)
        call copy(X_shift(1), X(1))
        allocate(H_shift(kdim+1, kdim)) ; H_shift = zero_rsp
        call arnoldi(A_shift, X_shift, H_shift, info_shift, tol=atol_sp)
        call check_info(info_shift, 'arnoldi', module=this_module_long, procedure='test_arnoldi_shifted_matrix_rsp')

        ! 5. Verify that bases are the same.
        G_basis = innerprod(X(:kdim), X_shift(:kdim))
        err = maxval(abs(G_basis - eye(kdim, mold=1.0_sp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_arnoldi_shifted_matrix_rsp', &
                              & info='Basis Invariance', eq='X = X_shift', context=msg)

        ! 6. Verify that H_shift = H + sigma*I.
        allocate(H_diff(kdim+1, kdim))
        H_diff = H_shift - H
        do i = 1, kdim
            H_diff(i, i) = H_diff(i, i) - sigma
        end do
        err = maxval(abs(H_diff))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_arnoldi_shifted_matrix_rsp', &
                              & info='Hessenberg Shift', eq='H_shift = H + sigma*I', context=msg)

        return
    end subroutine test_arnoldi_shifted_matrix_rsp

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

        block
        type(vector_rsp), allocatable :: Xfull(:), Xrestart(:)
        integer :: kdim_, kstart
        ! Hessenberg matrix.
        real(sp), allocatable :: Hfull(:, :), Hrestart(:, :)
        ! Miscellaneous.
        real(sp), allocatable :: Xfull_data(:, :), Xrestart_data(:, :)
        real(sp), allocatable :: G(:, :)
        real(sp) :: err

        ! Block Arnoldi parameters.
        kdim_ = test_size/p - 1
        kstart = kdim_/2

        ! Initialize full block Krylov subspace.
        deallocate(X0) ; allocate(Xfull(p*(kdim_+1)), X0(p))
        call init_rand(X0) ; call initialize_krylov_subspace(Xfull, X0)
        allocate(Hfull(p*(kdim_+1), p*kdim_), source=zero_rsp)
        allocate(Hrestart(p*(kdim_+1), p*kdim_), source=zero_rsp)

        ! Full block Arnoldi factorization.
        call arnoldi(A, Xfull, Hfull, info, blksize=p, tol=atol_sp)
        call check_info(info, 'arnoldi', module=this_module_long, procedure='test_block_arnoldi_factorization_rsp')

        ! Copy data for restart.
        allocate(Xrestart(p*(kdim_+1))) ; call zero_basis(Xrestart)
        call copy(Xrestart(:kstart*p), Xfull(:kstart*p))
        Hrestart(:kstart*p, :kstart*p-1) = Hfull(:kstart*p, :kstart*p-1)

        ! Restart block Arnoldi factorization.
        call arnoldi(A, Xrestart, Hrestart, info, kstart=kstart, blksize=p, tol=atol_sp)

        ! Compute inner product between the two bases.
        G = innerprod(Xfull(:p*kdim_), Xrestart(:p*kdim_))
        err = maxval(abs(G - eye(p*kdim_, mold=1.0_sp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_block_arnoldi_factorization_rsp', &
                                info='Restart', eq='Xfull = Xrestart', context=msg)

        ! Check Hessenberg matrices (compare the first kstart*p rows/cols).
        err = maxval(abs(Hfull(:kstart*p+1, :kstart*p-1) - Hrestart(:kstart*p+1, :kstart*p-1)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_block_arnoldi_factorization_rsp', &
                                info='Restart', eq='Hfull = Hrestart', context=msg)

        end block

        block
         integer :: kdim_, k_inv, i

        k_inv = 2  ! Multiple of block size p=2
        kdim_ = test_size/p - 1

        ! Create block lower triangular matrix.
        A%data(k_inv+1:, :k_inv) = zero_rsp

        ! Initialize starting vectors confined to the invariant subspace.
        deallocate(X, X0) ; allocate(X(p*(kdim_+1)), X0(p))
        call init_rand(X0)
        ! Zero out the components outside the invariant subspace.
        do i = 1, p
            call X0(i)%rand(ifnorm=.false.)
            X0(i)%data(k_inv+1:) = zero_rsp
        enddo
        call initialize_krylov_subspace(X, X0)

        deallocate(H) ; allocate(H(p*(kdim_+1), p*kdim_), source=zero_rsp)
        call arnoldi(A, X, H, info, blksize=p, tol=atol_sp)

        ! 1. Check if block Arnoldi detected the invariant subspace dimension.
        ! For block Arnoldi with p=2 and k_inv=2, we expect info = k_inv = 2.
        call check(error, info == k_inv)
        call check_test(error, 'test_block_arnoldi_factorization_rsp', &
                              & info='Subspace Dim', eq='info == k_inv', context='Invariant detection')

        ! 2. Check A @ X = X @ H for the computed invariant subspace.
        deallocate(Xdata) ; allocate(Xdata(test_size, info))
        call get_data(Xdata, X(:info))

        err = maxval(abs(matmul(A%data, Xdata) - matmul(Xdata, H(1:info, 1:info))))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_block_arnoldi_factorization_rsp', &
                              & info='Invariant Property', eq='A @ X = X @ H', context=msg)
        end block

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
            new_unittest("Arnoldi invalid parameters", test_arnoldi_invalid_params_rdp), &
            new_unittest("Arnoldi factorization", test_arnoldi_factorization_rdp), &
            new_unittest("Arnoldi shifted matrix", test_arnoldi_shifted_matrix_rdp), &
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

        ! Check H is indeed Hessenberg.
        call check(error, is_hessenberg(H, uplo='u'))
        call check_test(error, 'test_arnoldi_factorization_rdp', &
                               info='Upper Hessenberg', eq='Is H Hessenberg?', context=msg)

        ! Check correctness of full factorization.
        allocate(Xdata(test_size, kdim+1)) ; call get_data(Xdata, X)
        err = maxval(abs(matmul(A%data, Xdata(:, :kdim)) - matmul(Xdata, H)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_arnoldi_factorization_rdp', &
                              & info='Factorization', eq='A @ X = X_ @ H_', context=msg)

        ! Compute Gram matrix associated to the Krylov basis.
        G = Gram(X(:kdim))

        ! Check orthonormality of the computed basis.
        err = maxval(abs(G - eye(kdim, mold=1.0_dp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_arnoldi_factorization_rdp', &
                              & info='Orthonomality', eq='X.H @ X = I', context=msg)

        block
        ! Krylov subspaces.
        type(vector_rdp), allocatable :: Xfull(:), Xrestart(:)
        integer, parameter :: kstart = kdim/2
        ! Hessenberg matrix.
        real(dp), allocatable :: Hfull(:, :), Hrestart(:, :)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(dp), allocatable :: Xfull_data(:, :), Xrestart_data(:, :)

        ! Initialize Krylov subspace.
        allocate(Xfull(kdim+1)); call zero_basis(Xfull); call Xfull(1)%rand(ifnorm=.true.)
        allocate(Hfull(kdim+1, kdim), Hrestart(kdim+1, kdim), source=zero_rdp)

        ! Full Arnoldi factorization.
        call arnoldi(A, Xfull, Hfull, info, tol=atol_dp)
        call check_info(info, 'arnoldi', module=this_module_long, procedure='test_arnoldi_factorization_rdp')

        ! Copy data for restart.
        allocate(Xrestart(kdim+1)) ; call zero_basis(Xrestart)
        call copy(Xrestart(:kstart), Xfull(:kstart))
        Hrestart(:kstart, :kstart-1) = Hfull(:kstart, :kstart-1)

        ! Restart Arnoldi factorization.
        call arnoldi(A, Xrestart, Hrestart, info, kstart=kstart, tol=atol_dp)

        ! Compute inner product between the two bases.
        G = innerprod(Xfull(:kdim), Xrestart(:kdim))
        err = maxval(abs(G - eye(kdim, mold=1.0_dp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_arnoldi_factorization_rdp', &
                                info='Restart', eq='Xfull = Xrestart', context=msg)

        ! Check Hessenberg matrices.
        err = maxval(abs(Hfull - Hrestart))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_arnoldi_factorization_rdp', &
                                info='Restart', eq='Hfull = Hrestart', context=msg)

        end block

        block
        integer :: k_inv, k_max

        k_inv = 10
        k_max = test_size
        ! Create block lower triangular matrix.
        call init_rand(A)
        A%data(k_inv+1:, :k_inv) = zero_rdp

        ! Initial vector v1 must be in the invariant subspace
        deallocate(X) ; allocate(X(k_max+1))
        call zero_basis(X)
        call X(1)%rand(ifnorm=.false.)
        X(1)%data(k_inv+1:) = zero_rdp
        err = X(1)%norm()
        call X(1)%scal(one_rdp/err)

        deallocate(H) ; allocate(H(k_max+1, k_max), source=zero_rdp)
        call arnoldi(A, X, H, info, tol=atol_dp)

        ! 1. Check if Arnoldi detected the invariant subspace dimension
        call check(error, info == k_inv)
        call check_test(error, 'test_arnoldi_factorization_rdp', &
                              & info='Subspace Dim', eq='info == k_inv', context='Invariant detection')

        ! 2. Check AX_k = X_k H_k
        deallocate(Xdata) ; allocate(Xdata(test_size, info))
        call get_data(Xdata, X(:info))

        err = maxval(abs(matmul(A%data, Xdata) - matmul(Xdata, H(1:info, 1:info))))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_arnoldi_factorization_rdp', &
                              & info='Invariant Property', eq='A @ X = X @ H', context=msg)

        end block

        return
    end subroutine test_arnoldi_factorization_rdp

    subroutine test_arnoldi_invalid_params_rdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test linear operator.
        type(linop_rdp), allocatable :: A
        ! Krylov subspace.
        type(vector_rdp), allocatable :: X(:)
        ! Hessenberg matrix.
        real(dp), allocatable :: H(:, :)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        integer, parameter :: kdim = test_size
        integer :: p

        ! Common initialization for each case.
        do
            A = linop_rdp() ; call init_rand(A)
            allocate(H(kdim+1, kdim)) ; H = zero_rdp
            p = 1

            ! --- Case 1: k_start < 1 or k_start > k_end (info = -5) ---
            allocate(X(kdim+1)) ; call zero_basis(X)
            call arnoldi(A, X, H, info, kstart=0, tol=atol_dp)
            call check(error, info == -5)
            if (allocated(error)) exit
            call arnoldi(A, X, H, info, kstart=kdim, kend=1, tol=atol_dp)
            call check(error, info == -5)
            if (allocated(error)) exit
            ! --- Case 2: k_end > kdim (info = -6) ---
            call arnoldi(A, X, H, info, kend=kdim+1, tol=atol_dp)
            call check(error, info == -6)
            if (allocated(error)) exit
            ! --- Case 3: tolerance < 0 (info = -7) ---
            call arnoldi(A, X, H, info, tol=-1.0_dp)
            call check(error, info == -7)
            if (allocated(error)) exit
            ! --- Case 4: mod(size(X), p) /= 0 (info = -9) ---
            ! If kdim=20, size(X)=21. p=2 does not divide 21.
            p = 2
            call arnoldi(A, X, H, info, blksize=p, tol=atol_dp)
            call check(error, info == -9)
            if (allocated(error)) exit
            ! --- Case 5: blksize <= 0 (info = -9) ---
            p = 0
            call arnoldi(A, X, H, info, blksize=p, tol=atol_dp)
            call check(error, info == -9)
            if (allocated(error)) exit
            ! --- Case 6: H has wrong first dimension (info = -3) ---
            p = 1
            deallocate(H)
            allocate(H(kdim, kdim)) ; H = zero_rdp
            call arnoldi(A, X, H, info, blksize=p, tol=atol_dp)
            call check(error, info == -3)
            if (allocated(error)) exit
            ! --- Case 7: H has wrong second dimension (info = -3) ---
            deallocate(H)
            allocate(H(kdim+1, kdim-1)) ; H = zero_rdp
            call arnoldi(A, X, H, info, blksize=p, tol=atol_dp)
            call check(error, info == -3)
            if (allocated(error)) exit
            ! --- Case 8: block Arnoldi with too-small H (info = -3) ---
            block
                integer, parameter :: p_block = 2
                integer, parameter :: kdim_block = test_size/2
                type(vector_rdp), allocatable :: X0(:)

                deallocate(X, H)
                allocate(X(p_block*(kdim_block+1))) ; allocate(X0(p_block))
                call init_rand(X0) ; call initialize_krylov_subspace(X, X0)
                allocate(H(p_block*(kdim_block+1) - 1, p_block*kdim_block)) ; H = zero_rdp
                call arnoldi(A, X, H, info, blksize=p_block, tol=atol_dp)
                call check(error, info == -3)
           end block
            exit
        enddo
        call check_test(error, 'test_arnoldi_invalid_params_rdp', &
                          & info='Invalid parameters', eq='', context='block p=2')

        return
    end subroutine test_arnoldi_invalid_params_rdp

    subroutine test_arnoldi_shifted_matrix_rdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Linear operators.
        type(linop_rdp), allocatable :: A
        type(axpby_linop_rdp), allocatable :: A_shift
        ! Krylov subspaces.
        type(vector_rdp), allocatable :: X(:), X_shift(:)
        ! Hessenberg matrices.
        real(dp), allocatable :: H(:, :), H_shift(:, :)
        ! Information flags.
        integer :: info, info_shift, i
        ! Miscellaneous.
        integer, parameter :: kdim = test_size
        real(dp) :: sigma
        real(dp), allocatable :: G_basis(:, :)
        real(dp), allocatable :: H_diff(:, :)
        real(dp) :: err
        character(len=256) :: msg

        ! 1. Initialize linear operator A.
        A = linop_rdp() ; call init_rand(A)

        ! 2. Arnoldi factorization for A.
        allocate(X(kdim+1)); call zero_basis(X); call X(1)%rand(ifnorm = .true.)
        allocate(H(kdim+1, kdim)) ; H = zero_rdp
        call arnoldi(A, X, H, info, tol=atol_dp)
        call check_info(info, 'arnoldi', module=this_module_long, procedure='test_arnoldi_shifted_matrix_rdp')

        ! 3. Construct shifted operator A_shift = A + sigma*I.
        call random_number(sigma)
        allocate(A_shift)
        A_shift%A = A
        A_shift%B = Id_rdp()
        A_shift%alpha = 1.0_dp
        A_shift%beta = sigma
        A_shift%transA = .false.
        A_shift%transB = .false.

        ! 4. Arnoldi factorization for A_shift using the same starting vector.
        allocate(X_shift(kdim+1)); call zero_basis(X_shift)
        call copy(X_shift(1), X(1))
        allocate(H_shift(kdim+1, kdim)) ; H_shift = zero_rdp
        call arnoldi(A_shift, X_shift, H_shift, info_shift, tol=atol_dp)
        call check_info(info_shift, 'arnoldi', module=this_module_long, procedure='test_arnoldi_shifted_matrix_rdp')

        ! 5. Verify that bases are the same.
        G_basis = innerprod(X(:kdim), X_shift(:kdim))
        err = maxval(abs(G_basis - eye(kdim, mold=1.0_dp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_arnoldi_shifted_matrix_rdp', &
                              & info='Basis Invariance', eq='X = X_shift', context=msg)

        ! 6. Verify that H_shift = H + sigma*I.
        allocate(H_diff(kdim+1, kdim))
        H_diff = H_shift - H
        do i = 1, kdim
            H_diff(i, i) = H_diff(i, i) - sigma
        end do
        err = maxval(abs(H_diff))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_arnoldi_shifted_matrix_rdp', &
                              & info='Hessenberg Shift', eq='H_shift = H + sigma*I', context=msg)

        return
    end subroutine test_arnoldi_shifted_matrix_rdp

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

        block
        type(vector_rdp), allocatable :: Xfull(:), Xrestart(:)
        integer :: kdim_, kstart
        ! Hessenberg matrix.
        real(dp), allocatable :: Hfull(:, :), Hrestart(:, :)
        ! Miscellaneous.
        real(dp), allocatable :: Xfull_data(:, :), Xrestart_data(:, :)
        real(dp), allocatable :: G(:, :)
        real(dp) :: err

        ! Block Arnoldi parameters.
        kdim_ = test_size/p - 1
        kstart = kdim_/2

        ! Initialize full block Krylov subspace.
        deallocate(X0) ; allocate(Xfull(p*(kdim_+1)), X0(p))
        call init_rand(X0) ; call initialize_krylov_subspace(Xfull, X0)
        allocate(Hfull(p*(kdim_+1), p*kdim_), source=zero_rdp)
        allocate(Hrestart(p*(kdim_+1), p*kdim_), source=zero_rdp)

        ! Full block Arnoldi factorization.
        call arnoldi(A, Xfull, Hfull, info, blksize=p, tol=atol_dp)
        call check_info(info, 'arnoldi', module=this_module_long, procedure='test_block_arnoldi_factorization_rdp')

        ! Copy data for restart.
        allocate(Xrestart(p*(kdim_+1))) ; call zero_basis(Xrestart)
        call copy(Xrestart(:kstart*p), Xfull(:kstart*p))
        Hrestart(:kstart*p, :kstart*p-1) = Hfull(:kstart*p, :kstart*p-1)

        ! Restart block Arnoldi factorization.
        call arnoldi(A, Xrestart, Hrestart, info, kstart=kstart, blksize=p, tol=atol_dp)

        ! Compute inner product between the two bases.
        G = innerprod(Xfull(:p*kdim_), Xrestart(:p*kdim_))
        err = maxval(abs(G - eye(p*kdim_, mold=1.0_dp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_block_arnoldi_factorization_rdp', &
                                info='Restart', eq='Xfull = Xrestart', context=msg)

        ! Check Hessenberg matrices (compare the first kstart*p rows/cols).
        err = maxval(abs(Hfull(:kstart*p+1, :kstart*p-1) - Hrestart(:kstart*p+1, :kstart*p-1)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_block_arnoldi_factorization_rdp', &
                                info='Restart', eq='Hfull = Hrestart', context=msg)

        end block

        block
         integer :: kdim_, k_inv, i

        k_inv = 2  ! Multiple of block size p=2
        kdim_ = test_size/p - 1

        ! Create block lower triangular matrix.
        A%data(k_inv+1:, :k_inv) = zero_rdp

        ! Initialize starting vectors confined to the invariant subspace.
        deallocate(X, X0) ; allocate(X(p*(kdim_+1)), X0(p))
        call init_rand(X0)
        ! Zero out the components outside the invariant subspace.
        do i = 1, p
            call X0(i)%rand(ifnorm=.false.)
            X0(i)%data(k_inv+1:) = zero_rdp
        enddo
        call initialize_krylov_subspace(X, X0)

        deallocate(H) ; allocate(H(p*(kdim_+1), p*kdim_), source=zero_rdp)
        call arnoldi(A, X, H, info, blksize=p, tol=atol_dp)

        ! 1. Check if block Arnoldi detected the invariant subspace dimension.
        ! For block Arnoldi with p=2 and k_inv=2, we expect info = k_inv = 2.
        call check(error, info == k_inv)
        call check_test(error, 'test_block_arnoldi_factorization_rdp', &
                              & info='Subspace Dim', eq='info == k_inv', context='Invariant detection')

        ! 2. Check A @ X = X @ H for the computed invariant subspace.
        deallocate(Xdata) ; allocate(Xdata(test_size, info))
        call get_data(Xdata, X(:info))

        err = maxval(abs(matmul(A%data, Xdata) - matmul(Xdata, H(1:info, 1:info))))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_block_arnoldi_factorization_rdp', &
                              & info='Invariant Property', eq='A @ X = X @ H', context=msg)
        end block

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
            new_unittest("Arnoldi invalid parameters", test_arnoldi_invalid_params_csp), &
            new_unittest("Arnoldi factorization", test_arnoldi_factorization_csp), &
            new_unittest("Arnoldi shifted matrix", test_arnoldi_shifted_matrix_csp), &
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

        ! Check H is indeed Hessenberg.
        call check(error, is_hessenberg(H, uplo='u'))
        call check_test(error, 'test_arnoldi_factorization_csp', &
                               info='Upper Hessenberg', eq='Is H Hessenberg?', context=msg)

        ! Check correctness of full factorization.
        allocate(Xdata(test_size, kdim+1)) ; call get_data(Xdata, X)
        err = maxval(abs(matmul(A%data, Xdata(:, :kdim)) - matmul(Xdata, H)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_arnoldi_factorization_csp', &
                              & info='Factorization', eq='A @ X = X_ @ H_', context=msg)

        ! Compute Gram matrix associated to the Krylov basis.
        G = Gram(X(:kdim))

        ! Check orthonormality of the computed basis.
        err = maxval(abs(G - eye(kdim, mold=1.0_sp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_arnoldi_factorization_csp', &
                              & info='Orthonomality', eq='X.H @ X = I', context=msg)

        block
        ! Krylov subspaces.
        type(vector_csp), allocatable :: Xfull(:), Xrestart(:)
        integer, parameter :: kstart = kdim/2
        ! Hessenberg matrix.
        complex(sp), allocatable :: Hfull(:, :), Hrestart(:, :)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        complex(sp), allocatable :: Xfull_data(:, :), Xrestart_data(:, :)

        ! Initialize Krylov subspace.
        allocate(Xfull(kdim+1)); call zero_basis(Xfull); call Xfull(1)%rand(ifnorm=.true.)
        allocate(Hfull(kdim+1, kdim), Hrestart(kdim+1, kdim), source=zero_csp)

        ! Full Arnoldi factorization.
        call arnoldi(A, Xfull, Hfull, info, tol=atol_sp)
        call check_info(info, 'arnoldi', module=this_module_long, procedure='test_arnoldi_factorization_csp')

        ! Copy data for restart.
        allocate(Xrestart(kdim+1)) ; call zero_basis(Xrestart)
        call copy(Xrestart(:kstart), Xfull(:kstart))
        Hrestart(:kstart, :kstart-1) = Hfull(:kstart, :kstart-1)

        ! Restart Arnoldi factorization.
        call arnoldi(A, Xrestart, Hrestart, info, kstart=kstart, tol=atol_sp)

        ! Compute inner product between the two bases.
        G = innerprod(Xfull(:kdim), Xrestart(:kdim))
        err = maxval(abs(G - eye(kdim, mold=1.0_sp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_arnoldi_factorization_csp', &
                                info='Restart', eq='Xfull = Xrestart', context=msg)

        ! Check Hessenberg matrices.
        err = maxval(abs(Hfull - Hrestart))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_arnoldi_factorization_csp', &
                                info='Restart', eq='Hfull = Hrestart', context=msg)

        end block

        block
        integer :: k_inv, k_max

        k_inv = 10
        k_max = test_size
        ! Create block lower triangular matrix.
        call init_rand(A)
        A%data(k_inv+1:, :k_inv) = zero_csp

        ! Initial vector v1 must be in the invariant subspace
        deallocate(X) ; allocate(X(k_max+1))
        call zero_basis(X)
        call X(1)%rand(ifnorm=.false.)
        X(1)%data(k_inv+1:) = zero_csp
        err = X(1)%norm()
        call X(1)%scal(one_csp/err)

        deallocate(H) ; allocate(H(k_max+1, k_max), source=zero_csp)
        call arnoldi(A, X, H, info, tol=atol_sp)

        ! 1. Check if Arnoldi detected the invariant subspace dimension
        call check(error, info == k_inv)
        call check_test(error, 'test_arnoldi_factorization_csp', &
                              & info='Subspace Dim', eq='info == k_inv', context='Invariant detection')

        ! 2. Check AX_k = X_k H_k
        deallocate(Xdata) ; allocate(Xdata(test_size, info))
        call get_data(Xdata, X(:info))

        err = maxval(abs(matmul(A%data, Xdata) - matmul(Xdata, H(1:info, 1:info))))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_arnoldi_factorization_csp', &
                              & info='Invariant Property', eq='A @ X = X @ H', context=msg)

        end block

        return
    end subroutine test_arnoldi_factorization_csp

    subroutine test_arnoldi_invalid_params_csp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test linear operator.
        type(linop_csp), allocatable :: A
        ! Krylov subspace.
        type(vector_csp), allocatable :: X(:)
        ! Hessenberg matrix.
        complex(sp), allocatable :: H(:, :)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        integer, parameter :: kdim = test_size
        integer :: p

        ! Common initialization for each case.
        do
            A = linop_csp() ; call init_rand(A)
            allocate(H(kdim+1, kdim)) ; H = zero_csp
            p = 1

            ! --- Case 1: k_start < 1 or k_start > k_end (info = -5) ---
            allocate(X(kdim+1)) ; call zero_basis(X)
            call arnoldi(A, X, H, info, kstart=0, tol=atol_sp)
            call check(error, info == -5)
            if (allocated(error)) exit
            call arnoldi(A, X, H, info, kstart=kdim, kend=1, tol=atol_sp)
            call check(error, info == -5)
            if (allocated(error)) exit
            ! --- Case 2: k_end > kdim (info = -6) ---
            call arnoldi(A, X, H, info, kend=kdim+1, tol=atol_sp)
            call check(error, info == -6)
            if (allocated(error)) exit
            ! --- Case 3: tolerance < 0 (info = -7) ---
            call arnoldi(A, X, H, info, tol=-1.0_sp)
            call check(error, info == -7)
            if (allocated(error)) exit
            ! --- Case 4: mod(size(X), p) /= 0 (info = -9) ---
            ! If kdim=20, size(X)=21. p=2 does not divide 21.
            p = 2
            call arnoldi(A, X, H, info, blksize=p, tol=atol_sp)
            call check(error, info == -9)
            if (allocated(error)) exit
            ! --- Case 5: blksize <= 0 (info = -9) ---
            p = 0
            call arnoldi(A, X, H, info, blksize=p, tol=atol_sp)
            call check(error, info == -9)
            if (allocated(error)) exit
            ! --- Case 6: H has wrong first dimension (info = -3) ---
            p = 1
            deallocate(H)
            allocate(H(kdim, kdim)) ; H = zero_csp
            call arnoldi(A, X, H, info, blksize=p, tol=atol_sp)
            call check(error, info == -3)
            if (allocated(error)) exit
            ! --- Case 7: H has wrong second dimension (info = -3) ---
            deallocate(H)
            allocate(H(kdim+1, kdim-1)) ; H = zero_csp
            call arnoldi(A, X, H, info, blksize=p, tol=atol_sp)
            call check(error, info == -3)
            if (allocated(error)) exit
            ! --- Case 8: block Arnoldi with too-small H (info = -3) ---
            block
                integer, parameter :: p_block = 2
                integer, parameter :: kdim_block = test_size/2
                type(vector_csp), allocatable :: X0(:)

                deallocate(X, H)
                allocate(X(p_block*(kdim_block+1))) ; allocate(X0(p_block))
                call init_rand(X0) ; call initialize_krylov_subspace(X, X0)
                allocate(H(p_block*(kdim_block+1) - 1, p_block*kdim_block)) ; H = zero_csp
                call arnoldi(A, X, H, info, blksize=p_block, tol=atol_sp)
                call check(error, info == -3)
           end block
            exit
        enddo
        call check_test(error, 'test_arnoldi_invalid_params_csp', &
                          & info='Invalid parameters', eq='', context='block p=2')

        return
    end subroutine test_arnoldi_invalid_params_csp

    subroutine test_arnoldi_shifted_matrix_csp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Linear operators.
        type(linop_csp), allocatable :: A
        type(axpby_linop_csp), allocatable :: A_shift
        ! Krylov subspaces.
        type(vector_csp), allocatable :: X(:), X_shift(:)
        ! Hessenberg matrices.
        complex(sp), allocatable :: H(:, :), H_shift(:, :)
        ! Information flags.
        integer :: info, info_shift, i
        ! Miscellaneous.
        integer, parameter :: kdim = test_size
        complex(sp) :: sigma
        real(sp) :: sigma_r, sigma_i
        complex(sp), allocatable :: G_basis(:, :)
        complex(sp), allocatable :: H_diff(:, :)
        real(sp) :: err
        character(len=256) :: msg

        ! 1. Initialize linear operator A.
        A = linop_csp() ; call init_rand(A)

        ! 2. Arnoldi factorization for A.
        allocate(X(kdim+1)); call zero_basis(X); call X(1)%rand(ifnorm = .true.)
        allocate(H(kdim+1, kdim)) ; H = zero_csp
        call arnoldi(A, X, H, info, tol=atol_sp)
        call check_info(info, 'arnoldi', module=this_module_long, procedure='test_arnoldi_shifted_matrix_csp')

        ! 3. Construct shifted operator A_shift = A + sigma*I.
        call random_number(sigma_r)
        call random_number(sigma_i)
        sigma = cmplx(sigma_r, sigma_i, kind=sp)
        allocate(A_shift)
        A_shift%A = A
        A_shift%B = Id_csp()
        A_shift%alpha = 1.0_sp
        A_shift%beta = sigma
        A_shift%transA = .false.
        A_shift%transB = .false.

        ! 4. Arnoldi factorization for A_shift using the same starting vector.
        allocate(X_shift(kdim+1)); call zero_basis(X_shift)
        call copy(X_shift(1), X(1))
        allocate(H_shift(kdim+1, kdim)) ; H_shift = zero_csp
        call arnoldi(A_shift, X_shift, H_shift, info_shift, tol=atol_sp)
        call check_info(info_shift, 'arnoldi', module=this_module_long, procedure='test_arnoldi_shifted_matrix_csp')

        ! 5. Verify that bases are the same.
        G_basis = innerprod(X(:kdim), X_shift(:kdim))
        err = maxval(abs(G_basis - eye(kdim, mold=1.0_sp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_arnoldi_shifted_matrix_csp', &
                              & info='Basis Invariance', eq='X = X_shift', context=msg)

        ! 6. Verify that H_shift = H + sigma*I.
        allocate(H_diff(kdim+1, kdim))
        H_diff = H_shift - H
        do i = 1, kdim
            H_diff(i, i) = H_diff(i, i) - sigma
        end do
        err = maxval(abs(H_diff))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_arnoldi_shifted_matrix_csp', &
                              & info='Hessenberg Shift', eq='H_shift = H + sigma*I', context=msg)

        return
    end subroutine test_arnoldi_shifted_matrix_csp

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

        block
        type(vector_csp), allocatable :: Xfull(:), Xrestart(:)
        integer :: kdim_, kstart
        ! Hessenberg matrix.
        complex(sp), allocatable :: Hfull(:, :), Hrestart(:, :)
        ! Miscellaneous.
        complex(sp), allocatable :: Xfull_data(:, :), Xrestart_data(:, :)
        complex(sp), allocatable :: G(:, :)
        real(sp) :: err

        ! Block Arnoldi parameters.
        kdim_ = test_size/p - 1
        kstart = kdim_/2

        ! Initialize full block Krylov subspace.
        deallocate(X0) ; allocate(Xfull(p*(kdim_+1)), X0(p))
        call init_rand(X0) ; call initialize_krylov_subspace(Xfull, X0)
        allocate(Hfull(p*(kdim_+1), p*kdim_), source=zero_csp)
        allocate(Hrestart(p*(kdim_+1), p*kdim_), source=zero_csp)

        ! Full block Arnoldi factorization.
        call arnoldi(A, Xfull, Hfull, info, blksize=p, tol=atol_sp)
        call check_info(info, 'arnoldi', module=this_module_long, procedure='test_block_arnoldi_factorization_csp')

        ! Copy data for restart.
        allocate(Xrestart(p*(kdim_+1))) ; call zero_basis(Xrestart)
        call copy(Xrestart(:kstart*p), Xfull(:kstart*p))
        Hrestart(:kstart*p, :kstart*p-1) = Hfull(:kstart*p, :kstart*p-1)

        ! Restart block Arnoldi factorization.
        call arnoldi(A, Xrestart, Hrestart, info, kstart=kstart, blksize=p, tol=atol_sp)

        ! Compute inner product between the two bases.
        G = innerprod(Xfull(:p*kdim_), Xrestart(:p*kdim_))
        err = maxval(abs(G - eye(p*kdim_, mold=1.0_sp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_block_arnoldi_factorization_csp', &
                                info='Restart', eq='Xfull = Xrestart', context=msg)

        ! Check Hessenberg matrices (compare the first kstart*p rows/cols).
        err = maxval(abs(Hfull(:kstart*p+1, :kstart*p-1) - Hrestart(:kstart*p+1, :kstart*p-1)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_block_arnoldi_factorization_csp', &
                                info='Restart', eq='Hfull = Hrestart', context=msg)

        end block

        block
         integer :: kdim_, k_inv, i

        k_inv = 2  ! Multiple of block size p=2
        kdim_ = test_size/p - 1

        ! Create block lower triangular matrix.
        A%data(k_inv+1:, :k_inv) = zero_csp

        ! Initialize starting vectors confined to the invariant subspace.
        deallocate(X, X0) ; allocate(X(p*(kdim_+1)), X0(p))
        call init_rand(X0)
        ! Zero out the components outside the invariant subspace.
        do i = 1, p
            call X0(i)%rand(ifnorm=.false.)
            X0(i)%data(k_inv+1:) = zero_csp
        enddo
        call initialize_krylov_subspace(X, X0)

        deallocate(H) ; allocate(H(p*(kdim_+1), p*kdim_), source=zero_csp)
        call arnoldi(A, X, H, info, blksize=p, tol=atol_sp)

        ! 1. Check if block Arnoldi detected the invariant subspace dimension.
        ! For block Arnoldi with p=2 and k_inv=2, we expect info = k_inv = 2.
        call check(error, info == k_inv)
        call check_test(error, 'test_block_arnoldi_factorization_csp', &
                              & info='Subspace Dim', eq='info == k_inv', context='Invariant detection')

        ! 2. Check A @ X = X @ H for the computed invariant subspace.
        deallocate(Xdata) ; allocate(Xdata(test_size, info))
        call get_data(Xdata, X(:info))

        err = maxval(abs(matmul(A%data, Xdata) - matmul(Xdata, H(1:info, 1:info))))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_block_arnoldi_factorization_csp', &
                              & info='Invariant Property', eq='A @ X = X @ H', context=msg)
        end block

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
            new_unittest("Arnoldi invalid parameters", test_arnoldi_invalid_params_cdp), &
            new_unittest("Arnoldi factorization", test_arnoldi_factorization_cdp), &
            new_unittest("Arnoldi shifted matrix", test_arnoldi_shifted_matrix_cdp), &
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

        ! Check H is indeed Hessenberg.
        call check(error, is_hessenberg(H, uplo='u'))
        call check_test(error, 'test_arnoldi_factorization_cdp', &
                               info='Upper Hessenberg', eq='Is H Hessenberg?', context=msg)

        ! Check correctness of full factorization.
        allocate(Xdata(test_size, kdim+1)) ; call get_data(Xdata, X)
        err = maxval(abs(matmul(A%data, Xdata(:, :kdim)) - matmul(Xdata, H)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_arnoldi_factorization_cdp', &
                              & info='Factorization', eq='A @ X = X_ @ H_', context=msg)

        ! Compute Gram matrix associated to the Krylov basis.
        G = Gram(X(:kdim))

        ! Check orthonormality of the computed basis.
        err = maxval(abs(G - eye(kdim, mold=1.0_dp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_arnoldi_factorization_cdp', &
                              & info='Orthonomality', eq='X.H @ X = I', context=msg)

        block
        ! Krylov subspaces.
        type(vector_cdp), allocatable :: Xfull(:), Xrestart(:)
        integer, parameter :: kstart = kdim/2
        ! Hessenberg matrix.
        complex(dp), allocatable :: Hfull(:, :), Hrestart(:, :)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        complex(dp), allocatable :: Xfull_data(:, :), Xrestart_data(:, :)

        ! Initialize Krylov subspace.
        allocate(Xfull(kdim+1)); call zero_basis(Xfull); call Xfull(1)%rand(ifnorm=.true.)
        allocate(Hfull(kdim+1, kdim), Hrestart(kdim+1, kdim), source=zero_cdp)

        ! Full Arnoldi factorization.
        call arnoldi(A, Xfull, Hfull, info, tol=atol_dp)
        call check_info(info, 'arnoldi', module=this_module_long, procedure='test_arnoldi_factorization_cdp')

        ! Copy data for restart.
        allocate(Xrestart(kdim+1)) ; call zero_basis(Xrestart)
        call copy(Xrestart(:kstart), Xfull(:kstart))
        Hrestart(:kstart, :kstart-1) = Hfull(:kstart, :kstart-1)

        ! Restart Arnoldi factorization.
        call arnoldi(A, Xrestart, Hrestart, info, kstart=kstart, tol=atol_dp)

        ! Compute inner product between the two bases.
        G = innerprod(Xfull(:kdim), Xrestart(:kdim))
        err = maxval(abs(G - eye(kdim, mold=1.0_dp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_arnoldi_factorization_cdp', &
                                info='Restart', eq='Xfull = Xrestart', context=msg)

        ! Check Hessenberg matrices.
        err = maxval(abs(Hfull - Hrestart))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_arnoldi_factorization_cdp', &
                                info='Restart', eq='Hfull = Hrestart', context=msg)

        end block

        block
        integer :: k_inv, k_max

        k_inv = 10
        k_max = test_size
        ! Create block lower triangular matrix.
        call init_rand(A)
        A%data(k_inv+1:, :k_inv) = zero_cdp

        ! Initial vector v1 must be in the invariant subspace
        deallocate(X) ; allocate(X(k_max+1))
        call zero_basis(X)
        call X(1)%rand(ifnorm=.false.)
        X(1)%data(k_inv+1:) = zero_cdp
        err = X(1)%norm()
        call X(1)%scal(one_cdp/err)

        deallocate(H) ; allocate(H(k_max+1, k_max), source=zero_cdp)
        call arnoldi(A, X, H, info, tol=atol_dp)

        ! 1. Check if Arnoldi detected the invariant subspace dimension
        call check(error, info == k_inv)
        call check_test(error, 'test_arnoldi_factorization_cdp', &
                              & info='Subspace Dim', eq='info == k_inv', context='Invariant detection')

        ! 2. Check AX_k = X_k H_k
        deallocate(Xdata) ; allocate(Xdata(test_size, info))
        call get_data(Xdata, X(:info))

        err = maxval(abs(matmul(A%data, Xdata) - matmul(Xdata, H(1:info, 1:info))))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_arnoldi_factorization_cdp', &
                              & info='Invariant Property', eq='A @ X = X @ H', context=msg)

        end block

        return
    end subroutine test_arnoldi_factorization_cdp

    subroutine test_arnoldi_invalid_params_cdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test linear operator.
        type(linop_cdp), allocatable :: A
        ! Krylov subspace.
        type(vector_cdp), allocatable :: X(:)
        ! Hessenberg matrix.
        complex(dp), allocatable :: H(:, :)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        integer, parameter :: kdim = test_size
        integer :: p

        ! Common initialization for each case.
        do
            A = linop_cdp() ; call init_rand(A)
            allocate(H(kdim+1, kdim)) ; H = zero_cdp
            p = 1

            ! --- Case 1: k_start < 1 or k_start > k_end (info = -5) ---
            allocate(X(kdim+1)) ; call zero_basis(X)
            call arnoldi(A, X, H, info, kstart=0, tol=atol_dp)
            call check(error, info == -5)
            if (allocated(error)) exit
            call arnoldi(A, X, H, info, kstart=kdim, kend=1, tol=atol_dp)
            call check(error, info == -5)
            if (allocated(error)) exit
            ! --- Case 2: k_end > kdim (info = -6) ---
            call arnoldi(A, X, H, info, kend=kdim+1, tol=atol_dp)
            call check(error, info == -6)
            if (allocated(error)) exit
            ! --- Case 3: tolerance < 0 (info = -7) ---
            call arnoldi(A, X, H, info, tol=-1.0_dp)
            call check(error, info == -7)
            if (allocated(error)) exit
            ! --- Case 4: mod(size(X), p) /= 0 (info = -9) ---
            ! If kdim=20, size(X)=21. p=2 does not divide 21.
            p = 2
            call arnoldi(A, X, H, info, blksize=p, tol=atol_dp)
            call check(error, info == -9)
            if (allocated(error)) exit
            ! --- Case 5: blksize <= 0 (info = -9) ---
            p = 0
            call arnoldi(A, X, H, info, blksize=p, tol=atol_dp)
            call check(error, info == -9)
            if (allocated(error)) exit
            ! --- Case 6: H has wrong first dimension (info = -3) ---
            p = 1
            deallocate(H)
            allocate(H(kdim, kdim)) ; H = zero_cdp
            call arnoldi(A, X, H, info, blksize=p, tol=atol_dp)
            call check(error, info == -3)
            if (allocated(error)) exit
            ! --- Case 7: H has wrong second dimension (info = -3) ---
            deallocate(H)
            allocate(H(kdim+1, kdim-1)) ; H = zero_cdp
            call arnoldi(A, X, H, info, blksize=p, tol=atol_dp)
            call check(error, info == -3)
            if (allocated(error)) exit
            ! --- Case 8: block Arnoldi with too-small H (info = -3) ---
            block
                integer, parameter :: p_block = 2
                integer, parameter :: kdim_block = test_size/2
                type(vector_cdp), allocatable :: X0(:)

                deallocate(X, H)
                allocate(X(p_block*(kdim_block+1))) ; allocate(X0(p_block))
                call init_rand(X0) ; call initialize_krylov_subspace(X, X0)
                allocate(H(p_block*(kdim_block+1) - 1, p_block*kdim_block)) ; H = zero_cdp
                call arnoldi(A, X, H, info, blksize=p_block, tol=atol_dp)
                call check(error, info == -3)
           end block
            exit
        enddo
        call check_test(error, 'test_arnoldi_invalid_params_cdp', &
                          & info='Invalid parameters', eq='', context='block p=2')

        return
    end subroutine test_arnoldi_invalid_params_cdp

    subroutine test_arnoldi_shifted_matrix_cdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Linear operators.
        type(linop_cdp), allocatable :: A
        type(axpby_linop_cdp), allocatable :: A_shift
        ! Krylov subspaces.
        type(vector_cdp), allocatable :: X(:), X_shift(:)
        ! Hessenberg matrices.
        complex(dp), allocatable :: H(:, :), H_shift(:, :)
        ! Information flags.
        integer :: info, info_shift, i
        ! Miscellaneous.
        integer, parameter :: kdim = test_size
        complex(dp) :: sigma
        real(dp) :: sigma_r, sigma_i
        complex(dp), allocatable :: G_basis(:, :)
        complex(dp), allocatable :: H_diff(:, :)
        real(dp) :: err
        character(len=256) :: msg

        ! 1. Initialize linear operator A.
        A = linop_cdp() ; call init_rand(A)

        ! 2. Arnoldi factorization for A.
        allocate(X(kdim+1)); call zero_basis(X); call X(1)%rand(ifnorm = .true.)
        allocate(H(kdim+1, kdim)) ; H = zero_cdp
        call arnoldi(A, X, H, info, tol=atol_dp)
        call check_info(info, 'arnoldi', module=this_module_long, procedure='test_arnoldi_shifted_matrix_cdp')

        ! 3. Construct shifted operator A_shift = A + sigma*I.
        call random_number(sigma_r)
        call random_number(sigma_i)
        sigma = cmplx(sigma_r, sigma_i, kind=dp)
        allocate(A_shift)
        A_shift%A = A
        A_shift%B = Id_cdp()
        A_shift%alpha = 1.0_dp
        A_shift%beta = sigma
        A_shift%transA = .false.
        A_shift%transB = .false.

        ! 4. Arnoldi factorization for A_shift using the same starting vector.
        allocate(X_shift(kdim+1)); call zero_basis(X_shift)
        call copy(X_shift(1), X(1))
        allocate(H_shift(kdim+1, kdim)) ; H_shift = zero_cdp
        call arnoldi(A_shift, X_shift, H_shift, info_shift, tol=atol_dp)
        call check_info(info_shift, 'arnoldi', module=this_module_long, procedure='test_arnoldi_shifted_matrix_cdp')

        ! 5. Verify that bases are the same.
        G_basis = innerprod(X(:kdim), X_shift(:kdim))
        err = maxval(abs(G_basis - eye(kdim, mold=1.0_dp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_arnoldi_shifted_matrix_cdp', &
                              & info='Basis Invariance', eq='X = X_shift', context=msg)

        ! 6. Verify that H_shift = H + sigma*I.
        allocate(H_diff(kdim+1, kdim))
        H_diff = H_shift - H
        do i = 1, kdim
            H_diff(i, i) = H_diff(i, i) - sigma
        end do
        err = maxval(abs(H_diff))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_arnoldi_shifted_matrix_cdp', &
                              & info='Hessenberg Shift', eq='H_shift = H + sigma*I', context=msg)

        return
    end subroutine test_arnoldi_shifted_matrix_cdp

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

        block
        type(vector_cdp), allocatable :: Xfull(:), Xrestart(:)
        integer :: kdim_, kstart
        ! Hessenberg matrix.
        complex(dp), allocatable :: Hfull(:, :), Hrestart(:, :)
        ! Miscellaneous.
        complex(dp), allocatable :: Xfull_data(:, :), Xrestart_data(:, :)
        complex(dp), allocatable :: G(:, :)
        real(dp) :: err

        ! Block Arnoldi parameters.
        kdim_ = test_size/p - 1
        kstart = kdim_/2

        ! Initialize full block Krylov subspace.
        deallocate(X0) ; allocate(Xfull(p*(kdim_+1)), X0(p))
        call init_rand(X0) ; call initialize_krylov_subspace(Xfull, X0)
        allocate(Hfull(p*(kdim_+1), p*kdim_), source=zero_cdp)
        allocate(Hrestart(p*(kdim_+1), p*kdim_), source=zero_cdp)

        ! Full block Arnoldi factorization.
        call arnoldi(A, Xfull, Hfull, info, blksize=p, tol=atol_dp)
        call check_info(info, 'arnoldi', module=this_module_long, procedure='test_block_arnoldi_factorization_cdp')

        ! Copy data for restart.
        allocate(Xrestart(p*(kdim_+1))) ; call zero_basis(Xrestart)
        call copy(Xrestart(:kstart*p), Xfull(:kstart*p))
        Hrestart(:kstart*p, :kstart*p-1) = Hfull(:kstart*p, :kstart*p-1)

        ! Restart block Arnoldi factorization.
        call arnoldi(A, Xrestart, Hrestart, info, kstart=kstart, blksize=p, tol=atol_dp)

        ! Compute inner product between the two bases.
        G = innerprod(Xfull(:p*kdim_), Xrestart(:p*kdim_))
        err = maxval(abs(G - eye(p*kdim_, mold=1.0_dp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_block_arnoldi_factorization_cdp', &
                                info='Restart', eq='Xfull = Xrestart', context=msg)

        ! Check Hessenberg matrices (compare the first kstart*p rows/cols).
        err = maxval(abs(Hfull(:kstart*p+1, :kstart*p-1) - Hrestart(:kstart*p+1, :kstart*p-1)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_block_arnoldi_factorization_cdp', &
                                info='Restart', eq='Hfull = Hrestart', context=msg)

        end block

        block
         integer :: kdim_, k_inv, i

        k_inv = 2  ! Multiple of block size p=2
        kdim_ = test_size/p - 1

        ! Create block lower triangular matrix.
        A%data(k_inv+1:, :k_inv) = zero_cdp

        ! Initialize starting vectors confined to the invariant subspace.
        deallocate(X, X0) ; allocate(X(p*(kdim_+1)), X0(p))
        call init_rand(X0)
        ! Zero out the components outside the invariant subspace.
        do i = 1, p
            call X0(i)%rand(ifnorm=.false.)
            X0(i)%data(k_inv+1:) = zero_cdp
        enddo
        call initialize_krylov_subspace(X, X0)

        deallocate(H) ; allocate(H(p*(kdim_+1), p*kdim_), source=zero_cdp)
        call arnoldi(A, X, H, info, blksize=p, tol=atol_dp)

        ! 1. Check if block Arnoldi detected the invariant subspace dimension.
        ! For block Arnoldi with p=2 and k_inv=2, we expect info = k_inv = 2.
        call check(error, info == k_inv)
        call check_test(error, 'test_block_arnoldi_factorization_cdp', &
                              & info='Subspace Dim', eq='info == k_inv', context='Invariant detection')

        ! 2. Check A @ X = X @ H for the computed invariant subspace.
        deallocate(Xdata) ; allocate(Xdata(test_size, info))
        call get_data(Xdata, X(:info))

        err = maxval(abs(matmul(A%data, Xdata) - matmul(Xdata, H(1:info, 1:info))))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_block_arnoldi_factorization_cdp', &
                              & info='Invariant Property', eq='A @ X = X @ H', context=msg)
        end block

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


    !-----------------------------------------------------------------------------------------
    !-----     DEFINITION OF THE UNIT TESTS FOR SAUNDER-SIMON-YIP TRIDIAGONALIZATION     -----
    !-----------------------------------------------------------------------------------------

    subroutine collect_ssy_tridiag_rsp_testsuite(testsuite)
        ! Collection of unit tests.
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
             new_unittest("Saunder-Simon-Yip Tridiagonalization", test_ssy_tridiag_factorization_rsp) &
            ]

        return
    end subroutine collect_ssy_tridiag_rsp_testsuite

     subroutine test_ssy_tridiag_factorization_rsp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test matrix.
        type(linop_rsp), allocatable :: A
        ! Krylov subspace.
        type(vector_rsp), allocatable :: U(:), V(:)
        ! Krylov subspace dimension.
        integer, parameter :: kdim = test_size
        ! Tridiagonal matrix.
        real(sp), allocatable :: T(:, :)
        ! Information flag.
        integer :: info

        ! Internal variables.
        real(sp), allocatable :: Udata(:, :), Vdata(:, :)
        real(sp), allocatable :: G(:, :)
        real(sp) :: err
        character(len=256) :: msg

        ! Initialize tridiagonal matrix.
        allocate(T(kdim+1, kdim+1)) ; T = zero_rsp

        ! Initialize operator.
        A = linop_rsp()
        call init_rand(A)

        ! Initialize Krylov subspace.
        allocate(U(kdim+1)); call zero_basis(U); call U(1)%rand(ifnorm = .false.)
        allocate(V(kdim+1)); call zero_basis(V); call V(1)%rand(ifnorm = .false.)

        ! Saunders-Simon-Yip tridiagonalization.
        call ssy(A, U, V, T, info, tol=atol_sp)
        call check_info(info, "ssy", module=this_module_long, &
                        procedure="test_ssy_tridiag_factorization_rsp")

        ! Orthogonality of the column-span basis.
        G = Gram(U(:kdim)) ; call save_npy("UG_matrix.npy", G)
        err = maxval(abs(G - eye(kdim, mold=1.0_sp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_ssy_tridiag_factorization_rsp', &
                                 & info='Orthonomality', eq='U.H @ U = I', context=msg)

        ! Orthogonality of the row-span basis.
        G = Gram(V(:kdim))
        err = maxval(abs(G - eye(kdim, mold=1.0_sp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_ssy_tridiag_factorization_rsp', &
                                 & info='Orthonomality', eq='V.H @ V = I', context=msg)

        ! Check correctness.
        allocate(Udata(test_size, kdim+1)) ; call get_data(Udata, U)
        allocate(Vdata(test_size, kdim+1)) ; call get_data(Vdata, V)

        ! Infinity-norm check.
        err = maxval(abs(matmul(A%data, Vdata(:, :kdim)) - matmul(Udata, T(:kdim+1, :kdim))))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_ssy_tridiag_factorization_rsp', &
                                 & info='Factorization', eq='A @ V = U_ @ T_', context=msg)

        err = maxval(abs(matmul(hermitian(A%data), Udata(:, :kdim)) - matmul(Vdata, hermitian(T(:kdim, :kdim+1)))))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_ssy_tridiag_factorization_rsp', &
                                 & info='Factorization', eq='A.H @ U = V_ @ T_.H', context=msg)
        return
    end subroutine test_ssy_tridiag_factorization_rsp
    subroutine collect_ssy_tridiag_rdp_testsuite(testsuite)
        ! Collection of unit tests.
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
             new_unittest("Saunder-Simon-Yip Tridiagonalization", test_ssy_tridiag_factorization_rdp) &
            ]

        return
    end subroutine collect_ssy_tridiag_rdp_testsuite

     subroutine test_ssy_tridiag_factorization_rdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test matrix.
        type(linop_rdp), allocatable :: A
        ! Krylov subspace.
        type(vector_rdp), allocatable :: U(:), V(:)
        ! Krylov subspace dimension.
        integer, parameter :: kdim = test_size
        ! Tridiagonal matrix.
        real(dp), allocatable :: T(:, :)
        ! Information flag.
        integer :: info

        ! Internal variables.
        real(dp), allocatable :: Udata(:, :), Vdata(:, :)
        real(dp), allocatable :: G(:, :)
        real(dp) :: err
        character(len=256) :: msg

        ! Initialize tridiagonal matrix.
        allocate(T(kdim+1, kdim+1)) ; T = zero_rdp

        ! Initialize operator.
        A = linop_rdp()
        call init_rand(A)

        ! Initialize Krylov subspace.
        allocate(U(kdim+1)); call zero_basis(U); call U(1)%rand(ifnorm = .false.)
        allocate(V(kdim+1)); call zero_basis(V); call V(1)%rand(ifnorm = .false.)

        ! Saunders-Simon-Yip tridiagonalization.
        call ssy(A, U, V, T, info, tol=atol_dp)
        call check_info(info, "ssy", module=this_module_long, &
                        procedure="test_ssy_tridiag_factorization_rdp")

        ! Orthogonality of the column-span basis.
        G = Gram(U(:kdim)) ; call save_npy("UG_matrix.npy", G)
        err = maxval(abs(G - eye(kdim, mold=1.0_dp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_ssy_tridiag_factorization_rdp', &
                                 & info='Orthonomality', eq='U.H @ U = I', context=msg)

        ! Orthogonality of the row-span basis.
        G = Gram(V(:kdim))
        err = maxval(abs(G - eye(kdim, mold=1.0_dp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_ssy_tridiag_factorization_rdp', &
                                 & info='Orthonomality', eq='V.H @ V = I', context=msg)

        ! Check correctness.
        allocate(Udata(test_size, kdim+1)) ; call get_data(Udata, U)
        allocate(Vdata(test_size, kdim+1)) ; call get_data(Vdata, V)

        ! Infinity-norm check.
        err = maxval(abs(matmul(A%data, Vdata(:, :kdim)) - matmul(Udata, T(:kdim+1, :kdim))))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_ssy_tridiag_factorization_rdp', &
                                 & info='Factorization', eq='A @ V = U_ @ T_', context=msg)

        err = maxval(abs(matmul(hermitian(A%data), Udata(:, :kdim)) - matmul(Vdata, hermitian(T(:kdim, :kdim+1)))))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_ssy_tridiag_factorization_rdp', &
                                 & info='Factorization', eq='A.H @ U = V_ @ T_.H', context=msg)
        return
    end subroutine test_ssy_tridiag_factorization_rdp
    subroutine collect_ssy_tridiag_csp_testsuite(testsuite)
        ! Collection of unit tests.
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
             new_unittest("Saunder-Simon-Yip Tridiagonalization", test_ssy_tridiag_factorization_csp) &
            ]

        return
    end subroutine collect_ssy_tridiag_csp_testsuite

     subroutine test_ssy_tridiag_factorization_csp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test matrix.
        type(linop_csp), allocatable :: A
        ! Krylov subspace.
        type(vector_csp), allocatable :: U(:), V(:)
        ! Krylov subspace dimension.
        integer, parameter :: kdim = test_size
        ! Tridiagonal matrix.
        complex(sp), allocatable :: T(:, :)
        ! Information flag.
        integer :: info

        ! Internal variables.
        complex(sp), allocatable :: Udata(:, :), Vdata(:, :)
        complex(sp), allocatable :: G(:, :)
        real(sp) :: err
        character(len=256) :: msg

        ! Initialize tridiagonal matrix.
        allocate(T(kdim+1, kdim+1)) ; T = zero_csp

        ! Initialize operator.
        A = linop_csp()
        call init_rand(A)

        ! Initialize Krylov subspace.
        allocate(U(kdim+1)); call zero_basis(U); call U(1)%rand(ifnorm = .false.)
        allocate(V(kdim+1)); call zero_basis(V); call V(1)%rand(ifnorm = .false.)

        ! Saunders-Simon-Yip tridiagonalization.
        call ssy(A, U, V, T, info, tol=atol_sp)
        call check_info(info, "ssy", module=this_module_long, &
                        procedure="test_ssy_tridiag_factorization_csp")

        ! Orthogonality of the column-span basis.
        G = Gram(U(:kdim)) ; call save_npy("UG_matrix.npy", G)
        err = maxval(abs(G - eye(kdim, mold=1.0_sp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_ssy_tridiag_factorization_csp', &
                                 & info='Orthonomality', eq='U.H @ U = I', context=msg)

        ! Orthogonality of the row-span basis.
        G = Gram(V(:kdim))
        err = maxval(abs(G - eye(kdim, mold=1.0_sp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_ssy_tridiag_factorization_csp', &
                                 & info='Orthonomality', eq='V.H @ V = I', context=msg)

        ! Check correctness.
        allocate(Udata(test_size, kdim+1)) ; call get_data(Udata, U)
        allocate(Vdata(test_size, kdim+1)) ; call get_data(Vdata, V)

        ! Infinity-norm check.
        err = maxval(abs(matmul(A%data, Vdata(:, :kdim)) - matmul(Udata, T(:kdim+1, :kdim))))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_ssy_tridiag_factorization_csp', &
                                 & info='Factorization', eq='A @ V = U_ @ T_', context=msg)

        err = maxval(abs(matmul(hermitian(A%data), Udata(:, :kdim)) - matmul(Vdata, hermitian(T(:kdim, :kdim+1)))))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_ssy_tridiag_factorization_csp', &
                                 & info='Factorization', eq='A.H @ U = V_ @ T_.H', context=msg)
        return
    end subroutine test_ssy_tridiag_factorization_csp
    subroutine collect_ssy_tridiag_cdp_testsuite(testsuite)
        ! Collection of unit tests.
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
             new_unittest("Saunder-Simon-Yip Tridiagonalization", test_ssy_tridiag_factorization_cdp) &
            ]

        return
    end subroutine collect_ssy_tridiag_cdp_testsuite

     subroutine test_ssy_tridiag_factorization_cdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test matrix.
        type(linop_cdp), allocatable :: A
        ! Krylov subspace.
        type(vector_cdp), allocatable :: U(:), V(:)
        ! Krylov subspace dimension.
        integer, parameter :: kdim = test_size
        ! Tridiagonal matrix.
        complex(dp), allocatable :: T(:, :)
        ! Information flag.
        integer :: info

        ! Internal variables.
        complex(dp), allocatable :: Udata(:, :), Vdata(:, :)
        complex(dp), allocatable :: G(:, :)
        real(dp) :: err
        character(len=256) :: msg

        ! Initialize tridiagonal matrix.
        allocate(T(kdim+1, kdim+1)) ; T = zero_cdp

        ! Initialize operator.
        A = linop_cdp()
        call init_rand(A)

        ! Initialize Krylov subspace.
        allocate(U(kdim+1)); call zero_basis(U); call U(1)%rand(ifnorm = .false.)
        allocate(V(kdim+1)); call zero_basis(V); call V(1)%rand(ifnorm = .false.)

        ! Saunders-Simon-Yip tridiagonalization.
        call ssy(A, U, V, T, info, tol=atol_dp)
        call check_info(info, "ssy", module=this_module_long, &
                        procedure="test_ssy_tridiag_factorization_cdp")

        ! Orthogonality of the column-span basis.
        G = Gram(U(:kdim)) ; call save_npy("UG_matrix.npy", G)
        err = maxval(abs(G - eye(kdim, mold=1.0_dp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_ssy_tridiag_factorization_cdp', &
                                 & info='Orthonomality', eq='U.H @ U = I', context=msg)

        ! Orthogonality of the row-span basis.
        G = Gram(V(:kdim))
        err = maxval(abs(G - eye(kdim, mold=1.0_dp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_ssy_tridiag_factorization_cdp', &
                                 & info='Orthonomality', eq='V.H @ V = I', context=msg)

        ! Check correctness.
        allocate(Udata(test_size, kdim+1)) ; call get_data(Udata, U)
        allocate(Vdata(test_size, kdim+1)) ; call get_data(Vdata, V)

        ! Infinity-norm check.
        err = maxval(abs(matmul(A%data, Vdata(:, :kdim)) - matmul(Udata, T(:kdim+1, :kdim))))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_ssy_tridiag_factorization_cdp', &
                                 & info='Factorization', eq='A @ V = U_ @ T_', context=msg)

        err = maxval(abs(matmul(hermitian(A%data), Udata(:, :kdim)) - matmul(Vdata, hermitian(T(:kdim, :kdim+1)))))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_ssy_tridiag_factorization_cdp', &
                                 & info='Factorization', eq='A.H @ U = V_ @ T_.H', context=msg)
        return
    end subroutine test_ssy_tridiag_factorization_cdp

    !-----------------------------------------------------------------------
    !-----     DEFINITIONS OF THE VARIOUS UNIT TESTS FOR UTILITIES     -----
    !-----------------------------------------------------------------------

    subroutine collect_krylov_utilities_rsp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [ &
            new_unittest("Orthonormalize basis", test_orthonormalize_basis_rsp), &
            new_unittest("Biorthonormalize bases", test_biorthonormalize_bases_rsp), &
            new_unittest("Biorthonormalize bases rank deficient", test_biorthonormalize_bases_rank_deficient_rsp) &
        ]
        return
    end subroutine collect_krylov_utilities_rsp_testsuite

    !-------------------------------------------------
    !-----     TEST ORTHONORMALIZE_BASIS          -----
    !-------------------------------------------------

    subroutine test_orthonormalize_basis_rsp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        integer, parameter :: kdim = test_size
        type(vector_rsp), allocatable :: X(:)
        ! Gram matrix.
        real(sp), allocatable :: G(:,:)
        ! Miscellaneous.
        real(sp) :: err
        character(len=256) :: msg

        ! Initialize random basis.
        allocate(X(kdim)); call init_rand(X)

        ! Orthonormalize in-place.
        call orthonormalize_basis(X)

        ! Check orthonormality via Gram matrix.
        G = Gram(X)
        err = norm2(abs(G - eye(kdim, mold=1.0_sp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_orthonormalize_basis_rsp', &
            & info='Basis orthonormality', eq='Q.H @ Q = I', context=msg)
        return
    end subroutine test_orthonormalize_basis_rsp

    !-------------------------------------------------
    !-----     TEST BIORTHONORMALIZE_BASES        -----
    !-------------------------------------------------

    subroutine test_biorthonormalize_bases_rsp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        integer, parameter :: kdim = test_size
        type(vector_rsp), allocatable :: X(:), Y(:)
        ! Cross-Gram matrix.
        real(sp), allocatable :: G(:,:)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(sp) :: err
        character(len=256) :: msg

        ! Initialize random bases.
        allocate(X(kdim), Y(kdim))
        call init_rand(X); call init_rand(Y)
        call orthonormalize_basis(X)
        call orthonormalize_basis(Y)

        ! Biorthonormalize in-place.
        call biorthonormalize_bases(X, Y, info=info)
        call check_info(info, 'biorthonormalize_bases', &
            & module=this_module_long, &
            & procedure='test_biorthonormalize_bases_rsp')

        ! Check biorthonormality: Y.H @ X = I
        G = innerprod(Y, X)
        err = maxval(abs(G - eye(info, mold=1.0_sp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_biorthonormalize_bases_rsp', &
            & info='Biorthonormality', eq='Y.H @ X = I', context=msg)
        return
    end subroutine test_biorthonormalize_bases_rsp

    !-----------------------------------------------------------
    !-----     TEST BIORTHONORMALIZE_BASES RANK DEFICIENT  -----
    !-----------------------------------------------------------

    subroutine test_biorthonormalize_bases_rank_deficient_rsp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        integer, parameter :: kdim = test_size
        integer, parameter :: rank = test_size / 2
        type(vector_rsp), allocatable :: X(:), Y(:)
        ! Cross-Gram matrix.
        real(sp), allocatable :: G(:,:)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(sp) :: err
        character(len=256) :: msg
        integer :: i

        ! Build rank-deficient bases: last kdim-rank vectors are linear combinations
        ! of the first rank vectors.
        allocate(X(kdim), Y(kdim))
        call zero_basis(X); call zero_basis(Y)
        call init_rand(X(:rank)); call init_rand(Y(:rank))

        ! Biorthonormalize: should detect rank deficiency and return nretain < kdim.
        call biorthonormalize_bases(X, Y, tol=atol_sp, info=info)

        ! info should equal rank (number of non-negligible singular values).
        call check(error, info == rank)
        call get_err_str(msg, "retained rank: ", real(info, sp))
        call check_test(error, 'test_biorthonormalize_bases_rank_deficient_rsp', &
            & info='Rank detection', eq='nretain == rank', context=msg)

        ! Check biorthonormality of the retained subspace.
        G = innerprod(Y(:info), X(:info))
        err = maxval(abs(G - eye(info, mold=1.0_sp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_biorthonormalize_bases_rank_deficient_rsp', &
            & info='Retained subspace biorth.', eq='Y.H @ X = I', context=msg)
        return
    end subroutine test_biorthonormalize_bases_rank_deficient_rsp

    subroutine collect_krylov_utilities_rdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [ &
            new_unittest("Orthonormalize basis", test_orthonormalize_basis_rdp), &
            new_unittest("Biorthonormalize bases", test_biorthonormalize_bases_rdp), &
            new_unittest("Biorthonormalize bases rank deficient", test_biorthonormalize_bases_rank_deficient_rdp) &
        ]
        return
    end subroutine collect_krylov_utilities_rdp_testsuite

    !-------------------------------------------------
    !-----     TEST ORTHONORMALIZE_BASIS          -----
    !-------------------------------------------------

    subroutine test_orthonormalize_basis_rdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        integer, parameter :: kdim = test_size
        type(vector_rdp), allocatable :: X(:)
        ! Gram matrix.
        real(dp), allocatable :: G(:,:)
        ! Miscellaneous.
        real(dp) :: err
        character(len=256) :: msg

        ! Initialize random basis.
        allocate(X(kdim)); call init_rand(X)

        ! Orthonormalize in-place.
        call orthonormalize_basis(X)

        ! Check orthonormality via Gram matrix.
        G = Gram(X)
        err = norm2(abs(G - eye(kdim, mold=1.0_dp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_orthonormalize_basis_rdp', &
            & info='Basis orthonormality', eq='Q.H @ Q = I', context=msg)
        return
    end subroutine test_orthonormalize_basis_rdp

    !-------------------------------------------------
    !-----     TEST BIORTHONORMALIZE_BASES        -----
    !-------------------------------------------------

    subroutine test_biorthonormalize_bases_rdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        integer, parameter :: kdim = test_size
        type(vector_rdp), allocatable :: X(:), Y(:)
        ! Cross-Gram matrix.
        real(dp), allocatable :: G(:,:)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(dp) :: err
        character(len=256) :: msg

        ! Initialize random bases.
        allocate(X(kdim), Y(kdim))
        call init_rand(X); call init_rand(Y)
        call orthonormalize_basis(X)
        call orthonormalize_basis(Y)

        ! Biorthonormalize in-place.
        call biorthonormalize_bases(X, Y, info=info)
        call check_info(info, 'biorthonormalize_bases', &
            & module=this_module_long, &
            & procedure='test_biorthonormalize_bases_rdp')

        ! Check biorthonormality: Y.H @ X = I
        G = innerprod(Y, X)
        err = maxval(abs(G - eye(info, mold=1.0_dp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_biorthonormalize_bases_rdp', &
            & info='Biorthonormality', eq='Y.H @ X = I', context=msg)
        return
    end subroutine test_biorthonormalize_bases_rdp

    !-----------------------------------------------------------
    !-----     TEST BIORTHONORMALIZE_BASES RANK DEFICIENT  -----
    !-----------------------------------------------------------

    subroutine test_biorthonormalize_bases_rank_deficient_rdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        integer, parameter :: kdim = test_size
        integer, parameter :: rank = test_size / 2
        type(vector_rdp), allocatable :: X(:), Y(:)
        ! Cross-Gram matrix.
        real(dp), allocatable :: G(:,:)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(dp) :: err
        character(len=256) :: msg
        integer :: i

        ! Build rank-deficient bases: last kdim-rank vectors are linear combinations
        ! of the first rank vectors.
        allocate(X(kdim), Y(kdim))
        call zero_basis(X); call zero_basis(Y)
        call init_rand(X(:rank)); call init_rand(Y(:rank))

        ! Biorthonormalize: should detect rank deficiency and return nretain < kdim.
        call biorthonormalize_bases(X, Y, tol=atol_dp, info=info)

        ! info should equal rank (number of non-negligible singular values).
        call check(error, info == rank)
        call get_err_str(msg, "retained rank: ", real(info, dp))
        call check_test(error, 'test_biorthonormalize_bases_rank_deficient_rdp', &
            & info='Rank detection', eq='nretain == rank', context=msg)

        ! Check biorthonormality of the retained subspace.
        G = innerprod(Y(:info), X(:info))
        err = maxval(abs(G - eye(info, mold=1.0_dp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_biorthonormalize_bases_rank_deficient_rdp', &
            & info='Retained subspace biorth.', eq='Y.H @ X = I', context=msg)
        return
    end subroutine test_biorthonormalize_bases_rank_deficient_rdp

    subroutine collect_krylov_utilities_csp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [ &
            new_unittest("Orthonormalize basis", test_orthonormalize_basis_csp), &
            new_unittest("Biorthonormalize bases", test_biorthonormalize_bases_csp), &
            new_unittest("Biorthonormalize bases rank deficient", test_biorthonormalize_bases_rank_deficient_csp) &
        ]
        return
    end subroutine collect_krylov_utilities_csp_testsuite

    !-------------------------------------------------
    !-----     TEST ORTHONORMALIZE_BASIS          -----
    !-------------------------------------------------

    subroutine test_orthonormalize_basis_csp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        integer, parameter :: kdim = test_size
        type(vector_csp), allocatable :: X(:)
        ! Gram matrix.
        complex(sp), allocatable :: G(:,:)
        ! Miscellaneous.
        real(sp) :: err
        character(len=256) :: msg

        ! Initialize random basis.
        allocate(X(kdim)); call init_rand(X)

        ! Orthonormalize in-place.
        call orthonormalize_basis(X)

        ! Check orthonormality via Gram matrix.
        G = Gram(X)
        err = norm2(abs(G - eye(kdim, mold=1.0_sp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_orthonormalize_basis_csp', &
            & info='Basis orthonormality', eq='Q.H @ Q = I', context=msg)
        return
    end subroutine test_orthonormalize_basis_csp

    !-------------------------------------------------
    !-----     TEST BIORTHONORMALIZE_BASES        -----
    !-------------------------------------------------

    subroutine test_biorthonormalize_bases_csp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        integer, parameter :: kdim = test_size
        type(vector_csp), allocatable :: X(:), Y(:)
        ! Cross-Gram matrix.
        complex(sp), allocatable :: G(:,:)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(sp) :: err
        character(len=256) :: msg

        ! Initialize random bases.
        allocate(X(kdim), Y(kdim))
        call init_rand(X); call init_rand(Y)
        call orthonormalize_basis(X)
        call orthonormalize_basis(Y)

        ! Biorthonormalize in-place.
        call biorthonormalize_bases(X, Y, info=info)
        call check_info(info, 'biorthonormalize_bases', &
            & module=this_module_long, &
            & procedure='test_biorthonormalize_bases_csp')

        ! Check biorthonormality: Y.H @ X = I
        G = innerprod(Y, X)
        err = maxval(abs(G - eye(info, mold=1.0_sp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_biorthonormalize_bases_csp', &
            & info='Biorthonormality', eq='Y.H @ X = I', context=msg)
        return
    end subroutine test_biorthonormalize_bases_csp

    !-----------------------------------------------------------
    !-----     TEST BIORTHONORMALIZE_BASES RANK DEFICIENT  -----
    !-----------------------------------------------------------

    subroutine test_biorthonormalize_bases_rank_deficient_csp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        integer, parameter :: kdim = test_size
        integer, parameter :: rank = test_size / 2
        type(vector_csp), allocatable :: X(:), Y(:)
        ! Cross-Gram matrix.
        complex(sp), allocatable :: G(:,:)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(sp) :: err
        character(len=256) :: msg
        integer :: i

        ! Build rank-deficient bases: last kdim-rank vectors are linear combinations
        ! of the first rank vectors.
        allocate(X(kdim), Y(kdim))
        call zero_basis(X); call zero_basis(Y)
        call init_rand(X(:rank)); call init_rand(Y(:rank))

        ! Biorthonormalize: should detect rank deficiency and return nretain < kdim.
        call biorthonormalize_bases(X, Y, tol=atol_sp, info=info)

        ! info should equal rank (number of non-negligible singular values).
        call check(error, info == rank)
        call get_err_str(msg, "retained rank: ", real(info, sp))
        call check_test(error, 'test_biorthonormalize_bases_rank_deficient_csp', &
            & info='Rank detection', eq='nretain == rank', context=msg)

        ! Check biorthonormality of the retained subspace.
        G = innerprod(Y(:info), X(:info))
        err = maxval(abs(G - eye(info, mold=1.0_sp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_biorthonormalize_bases_rank_deficient_csp', &
            & info='Retained subspace biorth.', eq='Y.H @ X = I', context=msg)
        return
    end subroutine test_biorthonormalize_bases_rank_deficient_csp

    subroutine collect_krylov_utilities_cdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [ &
            new_unittest("Orthonormalize basis", test_orthonormalize_basis_cdp), &
            new_unittest("Biorthonormalize bases", test_biorthonormalize_bases_cdp), &
            new_unittest("Biorthonormalize bases rank deficient", test_biorthonormalize_bases_rank_deficient_cdp) &
        ]
        return
    end subroutine collect_krylov_utilities_cdp_testsuite

    !-------------------------------------------------
    !-----     TEST ORTHONORMALIZE_BASIS          -----
    !-------------------------------------------------

    subroutine test_orthonormalize_basis_cdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        integer, parameter :: kdim = test_size
        type(vector_cdp), allocatable :: X(:)
        ! Gram matrix.
        complex(dp), allocatable :: G(:,:)
        ! Miscellaneous.
        real(dp) :: err
        character(len=256) :: msg

        ! Initialize random basis.
        allocate(X(kdim)); call init_rand(X)

        ! Orthonormalize in-place.
        call orthonormalize_basis(X)

        ! Check orthonormality via Gram matrix.
        G = Gram(X)
        err = norm2(abs(G - eye(kdim, mold=1.0_dp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_orthonormalize_basis_cdp', &
            & info='Basis orthonormality', eq='Q.H @ Q = I', context=msg)
        return
    end subroutine test_orthonormalize_basis_cdp

    !-------------------------------------------------
    !-----     TEST BIORTHONORMALIZE_BASES        -----
    !-------------------------------------------------

    subroutine test_biorthonormalize_bases_cdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        integer, parameter :: kdim = test_size
        type(vector_cdp), allocatable :: X(:), Y(:)
        ! Cross-Gram matrix.
        complex(dp), allocatable :: G(:,:)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(dp) :: err
        character(len=256) :: msg

        ! Initialize random bases.
        allocate(X(kdim), Y(kdim))
        call init_rand(X); call init_rand(Y)
        call orthonormalize_basis(X)
        call orthonormalize_basis(Y)

        ! Biorthonormalize in-place.
        call biorthonormalize_bases(X, Y, info=info)
        call check_info(info, 'biorthonormalize_bases', &
            & module=this_module_long, &
            & procedure='test_biorthonormalize_bases_cdp')

        ! Check biorthonormality: Y.H @ X = I
        G = innerprod(Y, X)
        err = maxval(abs(G - eye(info, mold=1.0_dp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_biorthonormalize_bases_cdp', &
            & info='Biorthonormality', eq='Y.H @ X = I', context=msg)
        return
    end subroutine test_biorthonormalize_bases_cdp

    !-----------------------------------------------------------
    !-----     TEST BIORTHONORMALIZE_BASES RANK DEFICIENT  -----
    !-----------------------------------------------------------

    subroutine test_biorthonormalize_bases_rank_deficient_cdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        integer, parameter :: kdim = test_size
        integer, parameter :: rank = test_size / 2
        type(vector_cdp), allocatable :: X(:), Y(:)
        ! Cross-Gram matrix.
        complex(dp), allocatable :: G(:,:)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(dp) :: err
        character(len=256) :: msg
        integer :: i

        ! Build rank-deficient bases: last kdim-rank vectors are linear combinations
        ! of the first rank vectors.
        allocate(X(kdim), Y(kdim))
        call zero_basis(X); call zero_basis(Y)
        call init_rand(X(:rank)); call init_rand(Y(:rank))

        ! Biorthonormalize: should detect rank deficiency and return nretain < kdim.
        call biorthonormalize_bases(X, Y, tol=atol_dp, info=info)

        ! info should equal rank (number of non-negligible singular values).
        call check(error, info == rank)
        call get_err_str(msg, "retained rank: ", real(info, dp))
        call check_test(error, 'test_biorthonormalize_bases_rank_deficient_cdp', &
            & info='Rank detection', eq='nretain == rank', context=msg)

        ! Check biorthonormality of the retained subspace.
        G = innerprod(Y(:info), X(:info))
        err = maxval(abs(G - eye(info, mold=1.0_dp)))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_biorthonormalize_bases_rank_deficient_cdp', &
            & info='Retained subspace biorth.', eq='Y.H @ X = I', context=msg)
        return
    end subroutine test_biorthonormalize_bases_rank_deficient_cdp

end module TestKrylov
