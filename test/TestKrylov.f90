module TestKrylov
    ! Fortran Standard Library.
    use iso_fortran_env
    use stdlib_math, only: is_close, all_close
    use stdlib_linalg, only: eye

    ! LightKrylov
    use LightKrylov

    ! Testdrive
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use TestVectors
    use TestLinops
    use TestUtils

    implicit none
    private

    public :: collect_qr_rsp_testsuite
    public :: collect_qr_rdp_testsuite
    public :: collect_qr_csp_testsuite
    public :: collect_qr_cdp_testsuite

contains

    !--------------------------------------------------------
    !-----     DEFINITIONS OF THE VARIOUS UNIT TEST     -----
    !--------------------------------------------------------

    subroutine collect_qr_rsp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                        new_unittest("QR factorization", test_qr_factorization_rsp), &
                        new_unittest("Pivoting QR (exact rank def.)", test_pivoting_qr_exact_rank_deficiency_rsp) &
                    ]
        return
    end subroutine collect_qr_rsp_testsuite

    subroutine test_qr_factorization_rsp(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        !> Test Vectors.
        integer, parameter :: kdim = test_size
        type(vector_rsp), allocatable :: A(:)
        !> Upper triangular matrix.
        real(sp) :: R(kdim, kdim) = 0.0_sp
        !> Information flag.
        integer :: info
        !> Miscellaneous.
        real(sp) :: Adata(test_size, kdim), Qdata(test_size, kdim)

        ! Initiliaze matrix.
        allocate(A(1:kdim)) ; call init_rand(A)

        ! Get data.
        call get_data(Adata, A)

        ! In-place QR factorization.
        call qr(A, R, info)

        ! Get data.
        call get_data(Qdata, A)

        ! Check correctness.
        call check(error, maxval(abs(Adata - matmul(Qdata, R))) < rtol_sp)

        return
    end subroutine test_qr_factorization_rsp

    subroutine test_pivoting_qr_exact_rank_deficiency_rsp(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        !> Test basis.
        type(vector_rsp), allocatable :: A(:)
        !> Krylov subspace dimension.
        integer, parameter :: kdim = 20
        !> Number of zero columns.
        integer, parameter :: nzero = 5
        !> Upper triangular matrix.
        real(sp) :: R(kdim, kdim)
        !> Permutation vector.
        integer :: perm(kdim)
        real(sp) :: Id(kdim, kdim)
        !> Information flag.
        integer :: info
       
        !> Miscellaneous.
        integer :: k, idx, rk
        real(sp) :: Adata(test_size, kdim), Qdata(test_size, kdim)
        real(sp) :: alpha
        logical :: mask(kdim)

        ! Effective rank.
        rk = kdim - nzero

        ! Initialize matrix.
        allocate(A(1:kdim)) ; call init_rand(A)

        ! Add zero vectors at random places.
        mask = .true. ; k = nzero
        do while (k > 0)
            call random_number(alpha)
            idx = 1 + floor(kdim*alpha)
            if (mask(idx)) then
                A(idx)%data = 0.0_sp
                mask(idx) = .false.
                k = k-1
            endif
        enddo

        ! Copy data.
        call get_data(Adata, A)

        ! In-place QR factorization.
        call qr(A, R, perm, info)

        ! Extract data
        call get_data(Qdata, A)
        Adata = Adata(:, perm)

        ! Check correctness.
        call check(error, maxval(abs(Adata - matmul(Qdata, R))) < rtol_sp)

        return
    end subroutine test_pivoting_qr_exact_rank_deficiency_rsp

    subroutine collect_qr_rdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                        new_unittest("QR factorization", test_qr_factorization_rdp), &
                        new_unittest("Pivoting QR (exact rank def.)", test_pivoting_qr_exact_rank_deficiency_rdp) &
                    ]
        return
    end subroutine collect_qr_rdp_testsuite

    subroutine test_qr_factorization_rdp(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        !> Test Vectors.
        integer, parameter :: kdim = test_size
        type(vector_rdp), allocatable :: A(:)
        !> Upper triangular matrix.
        real(dp) :: R(kdim, kdim) = 0.0_dp
        !> Information flag.
        integer :: info
        !> Miscellaneous.
        real(dp) :: Adata(test_size, kdim), Qdata(test_size, kdim)

        ! Initiliaze matrix.
        allocate(A(1:kdim)) ; call init_rand(A)

        ! Get data.
        call get_data(Adata, A)

        ! In-place QR factorization.
        call qr(A, R, info)

        ! Get data.
        call get_data(Qdata, A)

        ! Check correctness.
        call check(error, maxval(abs(Adata - matmul(Qdata, R))) < rtol_dp)

        return
    end subroutine test_qr_factorization_rdp

    subroutine test_pivoting_qr_exact_rank_deficiency_rdp(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        !> Test basis.
        type(vector_rdp), allocatable :: A(:)
        !> Krylov subspace dimension.
        integer, parameter :: kdim = 20
        !> Number of zero columns.
        integer, parameter :: nzero = 5
        !> Upper triangular matrix.
        real(dp) :: R(kdim, kdim)
        !> Permutation vector.
        integer :: perm(kdim)
        real(dp) :: Id(kdim, kdim)
        !> Information flag.
        integer :: info
       
        !> Miscellaneous.
        integer :: k, idx, rk
        real(dp) :: Adata(test_size, kdim), Qdata(test_size, kdim)
        real(dp) :: alpha
        logical :: mask(kdim)

        ! Effective rank.
        rk = kdim - nzero

        ! Initialize matrix.
        allocate(A(1:kdim)) ; call init_rand(A)

        ! Add zero vectors at random places.
        mask = .true. ; k = nzero
        do while (k > 0)
            call random_number(alpha)
            idx = 1 + floor(kdim*alpha)
            if (mask(idx)) then
                A(idx)%data = 0.0_dp
                mask(idx) = .false.
                k = k-1
            endif
        enddo

        ! Copy data.
        call get_data(Adata, A)

        ! In-place QR factorization.
        call qr(A, R, perm, info)

        ! Extract data
        call get_data(Qdata, A)
        Adata = Adata(:, perm)

        ! Check correctness.
        call check(error, maxval(abs(Adata - matmul(Qdata, R))) < rtol_dp)

        return
    end subroutine test_pivoting_qr_exact_rank_deficiency_rdp

    subroutine collect_qr_csp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                        new_unittest("QR factorization", test_qr_factorization_csp), &
                        new_unittest("Pivoting QR (exact rank def.)", test_pivoting_qr_exact_rank_deficiency_csp) &
                    ]
        return
    end subroutine collect_qr_csp_testsuite

    subroutine test_qr_factorization_csp(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        !> Test Vectors.
        integer, parameter :: kdim = test_size
        type(vector_csp), allocatable :: A(:)
        !> Upper triangular matrix.
        complex(sp) :: R(kdim, kdim) = 0.0_sp
        !> Information flag.
        integer :: info
        !> Miscellaneous.
        complex(sp) :: Adata(test_size, kdim), Qdata(test_size, kdim)

        ! Initiliaze matrix.
        allocate(A(1:kdim)) ; call init_rand(A)

        ! Get data.
        call get_data(Adata, A)

        ! In-place QR factorization.
        call qr(A, R, info)

        ! Get data.
        call get_data(Qdata, A)

        ! Check correctness.
        call check(error, maxval(abs(Adata - matmul(Qdata, R))) < rtol_sp)

        return
    end subroutine test_qr_factorization_csp

    subroutine test_pivoting_qr_exact_rank_deficiency_csp(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        !> Test basis.
        type(vector_csp), allocatable :: A(:)
        !> Krylov subspace dimension.
        integer, parameter :: kdim = 20
        !> Number of zero columns.
        integer, parameter :: nzero = 5
        !> Upper triangular matrix.
        complex(sp) :: R(kdim, kdim)
        !> Permutation vector.
        integer :: perm(kdim)
        complex(sp) :: Id(kdim, kdim)
        !> Information flag.
        integer :: info
       
        !> Miscellaneous.
        integer :: k, idx, rk
        complex(sp) :: Adata(test_size, kdim), Qdata(test_size, kdim)
        real(sp) :: alpha
        logical :: mask(kdim)

        ! Effective rank.
        rk = kdim - nzero

        ! Initialize matrix.
        allocate(A(1:kdim)) ; call init_rand(A)

        ! Add zero vectors at random places.
        mask = .true. ; k = nzero
        do while (k > 0)
            call random_number(alpha)
            idx = 1 + floor(kdim*alpha)
            if (mask(idx)) then
                A(idx)%data = 0.0_sp
                mask(idx) = .false.
                k = k-1
            endif
        enddo

        ! Copy data.
        call get_data(Adata, A)

        ! In-place QR factorization.
        call qr(A, R, perm, info)

        ! Extract data
        call get_data(Qdata, A)
        Adata = Adata(:, perm)

        ! Check correctness.
        call check(error, maxval(abs(Adata - matmul(Qdata, R))) < rtol_sp)

        return
    end subroutine test_pivoting_qr_exact_rank_deficiency_csp

    subroutine collect_qr_cdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                        new_unittest("QR factorization", test_qr_factorization_cdp), &
                        new_unittest("Pivoting QR (exact rank def.)", test_pivoting_qr_exact_rank_deficiency_cdp) &
                    ]
        return
    end subroutine collect_qr_cdp_testsuite

    subroutine test_qr_factorization_cdp(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        !> Test Vectors.
        integer, parameter :: kdim = test_size
        type(vector_cdp), allocatable :: A(:)
        !> Upper triangular matrix.
        complex(dp) :: R(kdim, kdim) = 0.0_dp
        !> Information flag.
        integer :: info
        !> Miscellaneous.
        complex(dp) :: Adata(test_size, kdim), Qdata(test_size, kdim)

        ! Initiliaze matrix.
        allocate(A(1:kdim)) ; call init_rand(A)

        ! Get data.
        call get_data(Adata, A)

        ! In-place QR factorization.
        call qr(A, R, info)

        ! Get data.
        call get_data(Qdata, A)

        ! Check correctness.
        call check(error, maxval(abs(Adata - matmul(Qdata, R))) < rtol_dp)

        return
    end subroutine test_qr_factorization_cdp

    subroutine test_pivoting_qr_exact_rank_deficiency_cdp(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        !> Test basis.
        type(vector_cdp), allocatable :: A(:)
        !> Krylov subspace dimension.
        integer, parameter :: kdim = 20
        !> Number of zero columns.
        integer, parameter :: nzero = 5
        !> Upper triangular matrix.
        complex(dp) :: R(kdim, kdim)
        !> Permutation vector.
        integer :: perm(kdim)
        complex(dp) :: Id(kdim, kdim)
        !> Information flag.
        integer :: info
       
        !> Miscellaneous.
        integer :: k, idx, rk
        complex(dp) :: Adata(test_size, kdim), Qdata(test_size, kdim)
        real(dp) :: alpha
        logical :: mask(kdim)

        ! Effective rank.
        rk = kdim - nzero

        ! Initialize matrix.
        allocate(A(1:kdim)) ; call init_rand(A)

        ! Add zero vectors at random places.
        mask = .true. ; k = nzero
        do while (k > 0)
            call random_number(alpha)
            idx = 1 + floor(kdim*alpha)
            if (mask(idx)) then
                A(idx)%data = 0.0_dp
                mask(idx) = .false.
                k = k-1
            endif
        enddo

        ! Copy data.
        call get_data(Adata, A)

        ! In-place QR factorization.
        call qr(A, R, perm, info)

        ! Extract data
        call get_data(Qdata, A)
        Adata = Adata(:, perm)

        ! Check correctness.
        call check(error, maxval(abs(Adata - matmul(Qdata, R))) < rtol_dp)

        return
    end subroutine test_pivoting_qr_exact_rank_deficiency_cdp


end module TestKrylov
