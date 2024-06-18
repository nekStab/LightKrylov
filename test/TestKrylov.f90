module TestKrylov
    ! Fortran Standard Library.
    use iso_fortran_env
    use stdlib_math, only: is_close, all_close
    use stdlib_linalg, only: eye
    use stdlib_stats, only: median

    ! LightKrylov
    use LightKrylov
    use LightKrylov_Constants
    use LightKrylov_Logger
    use LightKrylov_AbstractVectors

    ! Testdrive
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use TestVectors
    use TestLinops
    use TestUtils

    implicit none
    
    private

    character*128, parameter, private :: this_module = 'LightKrylov_TestKrylov'

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
                        new_unittest("QR factorization (correctness & orthonormality)", test_qr_factorization_rsp), &
                        new_unittest("Pivoting QR for a rank deficient matrix (correctness & orthonormality))",&
                            & test_pivoting_qr_exact_rank_deficiency_rsp) &
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
        real(sp), dimension(test_size, kdim) :: Adata, Qdata
        real(sp), dimension(kdim, kdim) :: G, Id

        ! Initiliaze matrix.
        allocate(A(1:kdim)) ; call init_rand(A)

        ! Get data.
        call get_data(Adata, A)

        ! In-place QR factorization.
        call qr(A, R, info)
        call check_info(info, 'qr', module=this_module, procedure='test_qr_factorization_rsp')

        ! Get data.
        call get_data(Qdata, A)

        ! Check correctness.
        call check(error, maxval(abs(Adata - matmul(Qdata, R))) < rtol_sp)
        call check_test(error, 'test_qr_factorization_rsp', &
                        & 'QR factorization not correct: A /= Q @ R')

        ! Compute Gram matrix associated to the Krylov basis.
        G = zero_rsp
        call innerprod(G, A(1:kdim), A(1:kdim))

        ! Check orthonormality of the computed basis.
        Id = eye(kdim)
        call check(error, norm2(abs(G - Id)) < rtol_sp)
        call check_test(error, 'test_qr_factorization_rsp', &
                        & 'Basis Q not orthonormal: Q.H @ Q /= I')

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
        real(sp), dimension(test_size, kdim) :: Adata, Qdata
        real(sp), dimension(kdim, kdim) :: G, Id

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
                A(idx)%data = zero_rsp
                mask(idx) = .false.
                k = k-1
            endif
        enddo

        ! Copy data.
        call get_data(Adata, A)

        ! In-place QR factorization.
        call qr(A, R, perm, info)
        call check_info(info, 'qr_pivot', module=this_module, procedure='test_pivoting_qr_exact_rank_deficiency_rsp')

        ! Extract data
        call get_data(Qdata, A)
        Adata = Adata(:, perm)

        ! Check correctness.
        call check(error, maxval(abs(Adata - matmul(Qdata, R))) < rtol_sp)
        call check_test(error, 'test_pivoting_qr_exact_rank_deficiency_rsp', &
                        & 'QR factorization not correct: A /= Q @ R')

        ! Compute Gram matrix associated to the Krylov basis.
        G = zero_rsp
        call innerprod(G, A(1:kdim), A(1:kdim))

        ! Check orthonormality of the computed basis.
        Id = eye(kdim)
        call check(error, norm2(abs(G - Id)) < rtol_sp)
        call check_test(error, 'test_pivoting_qr_exact_rank_deficiency_rsp', &
                        & 'Basis Q not orthonormal: Q.H @ Q /= I')

        return
    end subroutine test_pivoting_qr_exact_rank_deficiency_rsp

    subroutine collect_qr_rdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                        new_unittest("QR factorization (correctness & orthonormality)", test_qr_factorization_rdp), &
                        new_unittest("Pivoting QR for a rank deficient matrix (correctness & orthonormality))",&
                            & test_pivoting_qr_exact_rank_deficiency_rdp) &
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
        real(dp), dimension(test_size, kdim) :: Adata, Qdata
        real(dp), dimension(kdim, kdim) :: G, Id

        ! Initiliaze matrix.
        allocate(A(1:kdim)) ; call init_rand(A)

        ! Get data.
        call get_data(Adata, A)

        ! In-place QR factorization.
        call qr(A, R, info)
        call check_info(info, 'qr', module=this_module, procedure='test_qr_factorization_rdp')

        ! Get data.
        call get_data(Qdata, A)

        ! Check correctness.
        call check(error, maxval(abs(Adata - matmul(Qdata, R))) < rtol_dp)
        call check_test(error, 'test_qr_factorization_rdp', &
                        & 'QR factorization not correct: A /= Q @ R')

        ! Compute Gram matrix associated to the Krylov basis.
        G = zero_rdp
        call innerprod(G, A(1:kdim), A(1:kdim))

        ! Check orthonormality of the computed basis.
        Id = eye(kdim)
        call check(error, norm2(abs(G - Id)) < rtol_dp)
        call check_test(error, 'test_qr_factorization_rdp', &
                        & 'Basis Q not orthonormal: Q.H @ Q /= I')

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
        real(dp), dimension(test_size, kdim) :: Adata, Qdata
        real(dp), dimension(kdim, kdim) :: G, Id

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
                A(idx)%data = zero_rdp
                mask(idx) = .false.
                k = k-1
            endif
        enddo

        ! Copy data.
        call get_data(Adata, A)

        ! In-place QR factorization.
        call qr(A, R, perm, info)
        call check_info(info, 'qr_pivot', module=this_module, procedure='test_pivoting_qr_exact_rank_deficiency_rdp')

        ! Extract data
        call get_data(Qdata, A)
        Adata = Adata(:, perm)

        ! Check correctness.
        call check(error, maxval(abs(Adata - matmul(Qdata, R))) < rtol_dp)
        call check_test(error, 'test_pivoting_qr_exact_rank_deficiency_rdp', &
                        & 'QR factorization not correct: A /= Q @ R')

        ! Compute Gram matrix associated to the Krylov basis.
        G = zero_rdp
        call innerprod(G, A(1:kdim), A(1:kdim))

        ! Check orthonormality of the computed basis.
        Id = eye(kdim)
        call check(error, norm2(abs(G - Id)) < rtol_dp)
        call check_test(error, 'test_pivoting_qr_exact_rank_deficiency_rdp', &
                        & 'Basis Q not orthonormal: Q.H @ Q /= I')

        return
    end subroutine test_pivoting_qr_exact_rank_deficiency_rdp

    subroutine collect_qr_csp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                        new_unittest("QR factorization (correctness & orthonormality)", test_qr_factorization_csp), &
                        new_unittest("Pivoting QR for a rank deficient matrix (correctness & orthonormality))",&
                            & test_pivoting_qr_exact_rank_deficiency_csp) &
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
        complex(sp), dimension(test_size, kdim) :: Adata, Qdata
        complex(sp), dimension(kdim, kdim) :: G, Id

        ! Initiliaze matrix.
        allocate(A(1:kdim)) ; call init_rand(A)

        ! Get data.
        call get_data(Adata, A)

        ! In-place QR factorization.
        call qr(A, R, info)
        call check_info(info, 'qr', module=this_module, procedure='test_qr_factorization_csp')

        ! Get data.
        call get_data(Qdata, A)

        ! Check correctness.
        call check(error, maxval(abs(Adata - matmul(Qdata, R))) < rtol_sp)
        call check_test(error, 'test_qr_factorization_csp', &
                        & 'QR factorization not correct: A /= Q @ R')

        ! Compute Gram matrix associated to the Krylov basis.
        G = zero_csp
        call innerprod(G, A(1:kdim), A(1:kdim))

        ! Check orthonormality of the computed basis.
        Id = eye(kdim)
        call check(error, norm2(abs(G - Id)) < rtol_sp)
        call check_test(error, 'test_qr_factorization_csp', &
                        & 'Basis Q not orthonormal: Q.H @ Q /= I')

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
        complex(sp), dimension(test_size, kdim) :: Adata, Qdata
        complex(sp), dimension(kdim, kdim) :: G, Id

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
                A(idx)%data = zero_csp
                mask(idx) = .false.
                k = k-1
            endif
        enddo

        ! Copy data.
        call get_data(Adata, A)

        ! In-place QR factorization.
        call qr(A, R, perm, info)
        call check_info(info, 'qr_pivot', module=this_module, procedure='test_pivoting_qr_exact_rank_deficiency_csp')

        ! Extract data
        call get_data(Qdata, A)
        Adata = Adata(:, perm)

        ! Check correctness.
        call check(error, maxval(abs(Adata - matmul(Qdata, R))) < rtol_sp)
        call check_test(error, 'test_pivoting_qr_exact_rank_deficiency_csp', &
                        & 'QR factorization not correct: A /= Q @ R')

        ! Compute Gram matrix associated to the Krylov basis.
        G = zero_csp
        call innerprod(G, A(1:kdim), A(1:kdim))

        ! Check orthonormality of the computed basis.
        Id = eye(kdim)
        call check(error, norm2(abs(G - Id)) < rtol_sp)
        call check_test(error, 'test_pivoting_qr_exact_rank_deficiency_csp', &
                        & 'Basis Q not orthonormal: Q.H @ Q /= I')

        return
    end subroutine test_pivoting_qr_exact_rank_deficiency_csp

    subroutine collect_qr_cdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                        new_unittest("QR factorization (correctness & orthonormality)", test_qr_factorization_cdp), &
                        new_unittest("Pivoting QR for a rank deficient matrix (correctness & orthonormality))",&
                            & test_pivoting_qr_exact_rank_deficiency_cdp) &
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
        complex(dp), dimension(test_size, kdim) :: Adata, Qdata
        complex(dp), dimension(kdim, kdim) :: G, Id

        ! Initiliaze matrix.
        allocate(A(1:kdim)) ; call init_rand(A)

        ! Get data.
        call get_data(Adata, A)

        ! In-place QR factorization.
        call qr(A, R, info)
        call check_info(info, 'qr', module=this_module, procedure='test_qr_factorization_cdp')

        ! Get data.
        call get_data(Qdata, A)

        ! Check correctness.
        call check(error, maxval(abs(Adata - matmul(Qdata, R))) < rtol_dp)
        call check_test(error, 'test_qr_factorization_cdp', &
                        & 'QR factorization not correct: A /= Q @ R')

        ! Compute Gram matrix associated to the Krylov basis.
        G = zero_cdp
        call innerprod(G, A(1:kdim), A(1:kdim))

        ! Check orthonormality of the computed basis.
        Id = eye(kdim)
        call check(error, norm2(abs(G - Id)) < rtol_dp)
        call check_test(error, 'test_qr_factorization_cdp', &
                        & 'Basis Q not orthonormal: Q.H @ Q /= I')

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
        complex(dp), dimension(test_size, kdim) :: Adata, Qdata
        complex(dp), dimension(kdim, kdim) :: G, Id

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
                A(idx)%data = zero_cdp
                mask(idx) = .false.
                k = k-1
            endif
        enddo

        ! Copy data.
        call get_data(Adata, A)

        ! In-place QR factorization.
        call qr(A, R, perm, info)
        call check_info(info, 'qr_pivot', module=this_module, procedure='test_pivoting_qr_exact_rank_deficiency_cdp')

        ! Extract data
        call get_data(Qdata, A)
        Adata = Adata(:, perm)

        ! Check correctness.
        call check(error, maxval(abs(Adata - matmul(Qdata, R))) < rtol_dp)
        call check_test(error, 'test_pivoting_qr_exact_rank_deficiency_cdp', &
                        & 'QR factorization not correct: A /= Q @ R')

        ! Compute Gram matrix associated to the Krylov basis.
        G = zero_cdp
        call innerprod(G, A(1:kdim), A(1:kdim))

        ! Check orthonormality of the computed basis.
        Id = eye(kdim)
        call check(error, norm2(abs(G - Id)) < rtol_dp)
        call check_test(error, 'test_pivoting_qr_exact_rank_deficiency_cdp', &
                        & 'Basis Q not orthonormal: Q.H @ Q /= I')

        return
    end subroutine test_pivoting_qr_exact_rank_deficiency_cdp

    
    !--------------------------------------------------------------
    !-----     DEFINITIONS OF THE UNIT-TESTS FOR ARNOLDI      -----
    !--------------------------------------------------------------

    subroutine collect_arnoldi_rsp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("Arnoldi factorization (correctness & orthonormality)", test_arnoldi_factorization_rsp), &
            new_unittest("Block Arnoldi factorization (correctness & orthonormality)", test_block_arnoldi_factorization_rsp), &
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
        type(vector_rsp), dimension(:), allocatable :: X0
        integer, parameter :: kdim = test_size
        ! Hessenberg matrix.
        real(sp) :: H(kdim+1, kdim)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(sp) :: Xdata(test_size, kdim+1)
        real(sp), dimension(kdim, kdim) :: G, Id
           
        ! Initialize linear operator.
        A = linop_rsp() ; call init_rand(A)
        ! Initialize Krylov subspace.
        allocate(X(1:kdim+1)) ; allocate (X0(1))
        call init_rand(X0) ; call initialize_krylov_subspace(X, X0)
        H = zero_rsp
        ! Arnoldi factorization.
        call arnoldi(A, X, H, info)
        call check_info(info, 'arnoldi', module=this_module, procedure='test_arnoldi_factorization_rsp')

        ! Check correctness of full factorization.
        call get_data(Xdata, X)
        call check(error, maxval(abs(matmul(A%data, Xdata(:, 1:kdim)) - matmul(Xdata, H))) < rtol_sp)
        call check_test(error, 'test_arnoldi_factorization_rsp', &
                        & 'Arnoldi factorization not correct: A @ X[:,1:k] /= X[:,1:k+1] @ H[1:k+1,1:k]')


        ! Compute Gram matrix associated to the Krylov basis.
        G = zero_rsp
        call innerprod(G, X(1:kdim), X(1:kdim))

        ! Check orthonormality of the computed basis.
        Id = eye(kdim)
        call check(error, norm2(abs(G - Id)) < rtol_sp)
        call check_test(error, 'test_arnoldi_factorization_rsp', &
                        & 'Basis X not orthonormal: X[:,1:k].H @ X[:,1:k] /= I')

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
        real(sp) :: H(p*(kdim+1), p*kdim)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        type(vector_rsp), allocatable :: X0(:)
        real(sp), dimension(test_size, p*(kdim+1)) :: Xdata
        real(sp), dimension(p*kdim, p*kdim) :: G
        real(sp), dimension(p*kdim, p*kdim) :: Id

        ! Initialize linear operator.
        A = linop_rsp() ; call init_rand(A)

        ! Initialize Krylov subspace.
        allocate(X(1:p*(kdim+1))) ; allocate(X0(1:p))
        call init_rand(X0) ; call initialize_krylov_subspace(X, X0)
        H = zero_rsp

        ! Arnoldi factorization.
        call arnoldi(A, X, H, info, blksize=p)
        call check_info(info, 'arnoldi', module=this_module, procedure='test_block_arnoldi_factorization_rsp')

        ! Check correctness of full factorization.
        call get_data(Xdata, X)
        call check(error, maxval(abs(matmul(A%data, Xdata(:, 1:p*kdim)) - matmul(Xdata, H))) < rtol_sp)
        call check_test(error, 'test_block_arnoldi_factorization_rsp', &
                        & 'Arnoldi factorization not correct: A @ X[:,1:p*k] /= X[:,1:p*(k+1)] @ H[1:p*(k+1),1:k]')

        ! Compute Gram matrix associated to the Krylov basis.
        G = zero_rsp
        call innerprod(G, X(1:p*kdim), X(1:p*kdim))

        ! Check orthonormality of the computed basis.
        Id = eye(p*kdim)
        call check(error, norm2(abs(G - Id)) < rtol_sp)
        call check_test(error, 'test_block_arnoldi_factorization_rsp', &
                        & 'Basis X not orthonormal: X[:,1:k].H @ X[:,1:k] /= I')

        return
    end subroutine test_block_arnoldi_factorization_rsp

    subroutine test_krylov_schur_rsp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test operator.
        type(linop_rsp), allocatable :: A
        ! Krylov subspace.
        type(vector_rsp), allocatable :: X(:), X0(:)
        ! Krylov subspace dimension.
        integer, parameter :: kdim = 100
        ! Hessenberg matrix.
        real(sp) :: H(kdim+1, kdim)
        ! Information flag.
        integer :: info

        ! Miscellaneous.
        integer :: n
        real(sp) :: Xdata(test_size, kdim+1)
        real(sp) :: alpha

        ! Initialize matrix.
        A = linop_rsp() ; call init_rand(A)
        A%data = A%data / norm2(abs(A%data))

        ! Initialize Krylov subspace.
        allocate(X(kdim+1)) ; allocate(X0(1))
        call init_rand(X0) ; call initialize_krylov_subspace(X, X0)
        H = zero_rsp

        ! Arnoldi factorization.
        call arnoldi(A, X, H, info)
        call check_info(info, 'arnoldi', module=this_module, procedure='test_krylov_schur_rsp')

        ! Krylov-Schur condensation.
        call krylov_schur(n, X, H, select_eigs)

        ! Check correctness.
        call get_data(Xdata, X)
        alpha = maxval(abs(matmul(A%data, Xdata(:, :n)) - matmul(Xdata(:, :n+1), H(:n+1, :n))))
        call check(error, alpha < rtol_sp)
        call check_test(error, 'test_krylov_schur_rsp', &
                        & 'Arnoldi factorization not correct: A @ X[:,1:k] /= X[:,1:k+1] @ H[1:k+1,1:k]')

        return
    contains
        function select_eigs(eigvals) result(selected)
            complex(sp), intent(in) :: eigvals(:)
            logical                       :: selected(size(eigvals))
           
            selected = abs(eigvals) > median(abs(eigvals))
            return
        end function select_eigs
    end subroutine test_krylov_schur_rsp

    subroutine collect_arnoldi_rdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("Arnoldi factorization (correctness & orthonormality)", test_arnoldi_factorization_rdp), &
            new_unittest("Block Arnoldi factorization (correctness & orthonormality)", test_block_arnoldi_factorization_rdp), &
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
        type(vector_rdp), dimension(:), allocatable :: X0
        integer, parameter :: kdim = test_size
        ! Hessenberg matrix.
        real(dp) :: H(kdim+1, kdim)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(dp) :: Xdata(test_size, kdim+1)
        real(dp), dimension(kdim, kdim) :: G, Id
           
        ! Initialize linear operator.
        A = linop_rdp() ; call init_rand(A)
        ! Initialize Krylov subspace.
        allocate(X(1:kdim+1)) ; allocate (X0(1))
        call init_rand(X0) ; call initialize_krylov_subspace(X, X0)
        H = zero_rdp
        ! Arnoldi factorization.
        call arnoldi(A, X, H, info)
        call check_info(info, 'arnoldi', module=this_module, procedure='test_arnoldi_factorization_rdp')

        ! Check correctness of full factorization.
        call get_data(Xdata, X)
        call check(error, maxval(abs(matmul(A%data, Xdata(:, 1:kdim)) - matmul(Xdata, H))) < rtol_dp)
        call check_test(error, 'test_arnoldi_factorization_rdp', &
                        & 'Arnoldi factorization not correct: A @ X[:,1:k] /= X[:,1:k+1] @ H[1:k+1,1:k]')


        ! Compute Gram matrix associated to the Krylov basis.
        G = zero_rdp
        call innerprod(G, X(1:kdim), X(1:kdim))

        ! Check orthonormality of the computed basis.
        Id = eye(kdim)
        call check(error, norm2(abs(G - Id)) < rtol_dp)
        call check_test(error, 'test_arnoldi_factorization_rdp', &
                        & 'Basis X not orthonormal: X[:,1:k].H @ X[:,1:k] /= I')

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
        real(dp) :: H(p*(kdim+1), p*kdim)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        type(vector_rdp), allocatable :: X0(:)
        real(dp), dimension(test_size, p*(kdim+1)) :: Xdata
        real(dp), dimension(p*kdim, p*kdim) :: G
        real(dp), dimension(p*kdim, p*kdim) :: Id

        ! Initialize linear operator.
        A = linop_rdp() ; call init_rand(A)

        ! Initialize Krylov subspace.
        allocate(X(1:p*(kdim+1))) ; allocate(X0(1:p))
        call init_rand(X0) ; call initialize_krylov_subspace(X, X0)
        H = zero_rdp

        ! Arnoldi factorization.
        call arnoldi(A, X, H, info, blksize=p)
        call check_info(info, 'arnoldi', module=this_module, procedure='test_block_arnoldi_factorization_rdp')

        ! Check correctness of full factorization.
        call get_data(Xdata, X)
        call check(error, maxval(abs(matmul(A%data, Xdata(:, 1:p*kdim)) - matmul(Xdata, H))) < rtol_dp)
        call check_test(error, 'test_block_arnoldi_factorization_rdp', &
                        & 'Arnoldi factorization not correct: A @ X[:,1:p*k] /= X[:,1:p*(k+1)] @ H[1:p*(k+1),1:k]')

        ! Compute Gram matrix associated to the Krylov basis.
        G = zero_rdp
        call innerprod(G, X(1:p*kdim), X(1:p*kdim))

        ! Check orthonormality of the computed basis.
        Id = eye(p*kdim)
        call check(error, norm2(abs(G - Id)) < rtol_dp)
        call check_test(error, 'test_block_arnoldi_factorization_rdp', &
                        & 'Basis X not orthonormal: X[:,1:k].H @ X[:,1:k] /= I')

        return
    end subroutine test_block_arnoldi_factorization_rdp

    subroutine test_krylov_schur_rdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test operator.
        type(linop_rdp), allocatable :: A
        ! Krylov subspace.
        type(vector_rdp), allocatable :: X(:), X0(:)
        ! Krylov subspace dimension.
        integer, parameter :: kdim = 100
        ! Hessenberg matrix.
        real(dp) :: H(kdim+1, kdim)
        ! Information flag.
        integer :: info

        ! Miscellaneous.
        integer :: n
        real(dp) :: Xdata(test_size, kdim+1)
        real(dp) :: alpha

        ! Initialize matrix.
        A = linop_rdp() ; call init_rand(A)
        A%data = A%data / norm2(abs(A%data))

        ! Initialize Krylov subspace.
        allocate(X(kdim+1)) ; allocate(X0(1))
        call init_rand(X0) ; call initialize_krylov_subspace(X, X0)
        H = zero_rdp

        ! Arnoldi factorization.
        call arnoldi(A, X, H, info)
        call check_info(info, 'arnoldi', module=this_module, procedure='test_krylov_schur_rdp')

        ! Krylov-Schur condensation.
        call krylov_schur(n, X, H, select_eigs)

        ! Check correctness.
        call get_data(Xdata, X)
        alpha = maxval(abs(matmul(A%data, Xdata(:, :n)) - matmul(Xdata(:, :n+1), H(:n+1, :n))))
        call check(error, alpha < rtol_dp)
        call check_test(error, 'test_krylov_schur_rdp', &
                        & 'Arnoldi factorization not correct: A @ X[:,1:k] /= X[:,1:k+1] @ H[1:k+1,1:k]')

        return
    contains
        function select_eigs(eigvals) result(selected)
            complex(dp), intent(in) :: eigvals(:)
            logical                       :: selected(size(eigvals))
           
            selected = abs(eigvals) > median(abs(eigvals))
            return
        end function select_eigs
    end subroutine test_krylov_schur_rdp

    subroutine collect_arnoldi_csp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("Arnoldi factorization (correctness & orthonormality)", test_arnoldi_factorization_csp), &
            new_unittest("Block Arnoldi factorization (correctness & orthonormality)", test_block_arnoldi_factorization_csp), &
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
        type(vector_csp), dimension(:), allocatable :: X0
        integer, parameter :: kdim = test_size
        ! Hessenberg matrix.
        complex(sp) :: H(kdim+1, kdim)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        complex(sp) :: Xdata(test_size, kdim+1)
        complex(sp), dimension(kdim, kdim) :: G, Id
           
        ! Initialize linear operator.
        A = linop_csp() ; call init_rand(A)
        ! Initialize Krylov subspace.
        allocate(X(1:kdim+1)) ; allocate (X0(1))
        call init_rand(X0) ; call initialize_krylov_subspace(X, X0)
        H = zero_csp
        ! Arnoldi factorization.
        call arnoldi(A, X, H, info)
        call check_info(info, 'arnoldi', module=this_module, procedure='test_arnoldi_factorization_csp')

        ! Check correctness of full factorization.
        call get_data(Xdata, X)
        call check(error, maxval(abs(matmul(A%data, Xdata(:, 1:kdim)) - matmul(Xdata, H))) < rtol_sp)
        call check_test(error, 'test_arnoldi_factorization_csp', &
                        & 'Arnoldi factorization not correct: A @ X[:,1:k] /= X[:,1:k+1] @ H[1:k+1,1:k]')


        ! Compute Gram matrix associated to the Krylov basis.
        G = zero_csp
        call innerprod(G, X(1:kdim), X(1:kdim))

        ! Check orthonormality of the computed basis.
        Id = eye(kdim)
        call check(error, norm2(abs(G - Id)) < rtol_sp)
        call check_test(error, 'test_arnoldi_factorization_csp', &
                        & 'Basis X not orthonormal: X[:,1:k].H @ X[:,1:k] /= I')

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
        complex(sp) :: H(p*(kdim+1), p*kdim)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        type(vector_csp), allocatable :: X0(:)
        complex(sp), dimension(test_size, p*(kdim+1)) :: Xdata
        complex(sp), dimension(p*kdim, p*kdim) :: G
        real(sp), dimension(p*kdim, p*kdim) :: Id

        ! Initialize linear operator.
        A = linop_csp() ; call init_rand(A)

        ! Initialize Krylov subspace.
        allocate(X(1:p*(kdim+1))) ; allocate(X0(1:p))
        call init_rand(X0) ; call initialize_krylov_subspace(X, X0)
        H = zero_csp

        ! Arnoldi factorization.
        call arnoldi(A, X, H, info, blksize=p)
        call check_info(info, 'arnoldi', module=this_module, procedure='test_block_arnoldi_factorization_csp')

        ! Check correctness of full factorization.
        call get_data(Xdata, X)
        call check(error, maxval(abs(matmul(A%data, Xdata(:, 1:p*kdim)) - matmul(Xdata, H))) < rtol_sp)
        call check_test(error, 'test_block_arnoldi_factorization_csp', &
                        & 'Arnoldi factorization not correct: A @ X[:,1:p*k] /= X[:,1:p*(k+1)] @ H[1:p*(k+1),1:k]')

        ! Compute Gram matrix associated to the Krylov basis.
        G = zero_csp
        call innerprod(G, X(1:p*kdim), X(1:p*kdim))

        ! Check orthonormality of the computed basis.
        Id = eye(p*kdim)
        call check(error, norm2(abs(G - Id)) < rtol_sp)
        call check_test(error, 'test_block_arnoldi_factorization_csp', &
                        & 'Basis X not orthonormal: X[:,1:k].H @ X[:,1:k] /= I')

        return
    end subroutine test_block_arnoldi_factorization_csp

    subroutine test_krylov_schur_csp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test operator.
        type(linop_csp), allocatable :: A
        ! Krylov subspace.
        type(vector_csp), allocatable :: X(:), X0(:)
        ! Krylov subspace dimension.
        integer, parameter :: kdim = 100
        ! Hessenberg matrix.
        complex(sp) :: H(kdim+1, kdim)
        ! Information flag.
        integer :: info

        ! Miscellaneous.
        integer :: n
        complex(sp) :: Xdata(test_size, kdim+1)
        real(sp) :: alpha

        ! Initialize matrix.
        A = linop_csp() ; call init_rand(A)
        A%data = A%data / norm2(abs(A%data))

        ! Initialize Krylov subspace.
        allocate(X(kdim+1)) ; allocate(X0(1))
        call init_rand(X0) ; call initialize_krylov_subspace(X, X0)
        H = zero_csp

        ! Arnoldi factorization.
        call arnoldi(A, X, H, info)
        call check_info(info, 'arnoldi', module=this_module, procedure='test_krylov_schur_csp')

        ! Krylov-Schur condensation.
        call krylov_schur(n, X, H, select_eigs)

        ! Check correctness.
        call get_data(Xdata, X)
        alpha = maxval(abs(matmul(A%data, Xdata(:, :n)) - matmul(Xdata(:, :n+1), H(:n+1, :n))))
        call check(error, alpha < rtol_sp)
        call check_test(error, 'test_krylov_schur_csp', &
                        & 'Arnoldi factorization not correct: A @ X[:,1:k] /= X[:,1:k+1] @ H[1:k+1,1:k]')

        return
    contains
        function select_eigs(eigvals) result(selected)
            complex(sp), intent(in) :: eigvals(:)
            logical                       :: selected(size(eigvals))
           
            selected = abs(eigvals) > median(abs(eigvals))
            return
        end function select_eigs
    end subroutine test_krylov_schur_csp

    subroutine collect_arnoldi_cdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("Arnoldi factorization (correctness & orthonormality)", test_arnoldi_factorization_cdp), &
            new_unittest("Block Arnoldi factorization (correctness & orthonormality)", test_block_arnoldi_factorization_cdp), &
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
        type(vector_cdp), dimension(:), allocatable :: X0
        integer, parameter :: kdim = test_size
        ! Hessenberg matrix.
        complex(dp) :: H(kdim+1, kdim)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        complex(dp) :: Xdata(test_size, kdim+1)
        complex(dp), dimension(kdim, kdim) :: G, Id
           
        ! Initialize linear operator.
        A = linop_cdp() ; call init_rand(A)
        ! Initialize Krylov subspace.
        allocate(X(1:kdim+1)) ; allocate (X0(1))
        call init_rand(X0) ; call initialize_krylov_subspace(X, X0)
        H = zero_cdp
        ! Arnoldi factorization.
        call arnoldi(A, X, H, info)
        call check_info(info, 'arnoldi', module=this_module, procedure='test_arnoldi_factorization_cdp')

        ! Check correctness of full factorization.
        call get_data(Xdata, X)
        call check(error, maxval(abs(matmul(A%data, Xdata(:, 1:kdim)) - matmul(Xdata, H))) < rtol_dp)
        call check_test(error, 'test_arnoldi_factorization_cdp', &
                        & 'Arnoldi factorization not correct: A @ X[:,1:k] /= X[:,1:k+1] @ H[1:k+1,1:k]')


        ! Compute Gram matrix associated to the Krylov basis.
        G = zero_cdp
        call innerprod(G, X(1:kdim), X(1:kdim))

        ! Check orthonormality of the computed basis.
        Id = eye(kdim)
        call check(error, norm2(abs(G - Id)) < rtol_dp)
        call check_test(error, 'test_arnoldi_factorization_cdp', &
                        & 'Basis X not orthonormal: X[:,1:k].H @ X[:,1:k] /= I')

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
        complex(dp) :: H(p*(kdim+1), p*kdim)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        type(vector_cdp), allocatable :: X0(:)
        complex(dp), dimension(test_size, p*(kdim+1)) :: Xdata
        complex(dp), dimension(p*kdim, p*kdim) :: G
        real(dp), dimension(p*kdim, p*kdim) :: Id

        ! Initialize linear operator.
        A = linop_cdp() ; call init_rand(A)

        ! Initialize Krylov subspace.
        allocate(X(1:p*(kdim+1))) ; allocate(X0(1:p))
        call init_rand(X0) ; call initialize_krylov_subspace(X, X0)
        H = zero_cdp

        ! Arnoldi factorization.
        call arnoldi(A, X, H, info, blksize=p)
        call check_info(info, 'arnoldi', module=this_module, procedure='test_block_arnoldi_factorization_cdp')

        ! Check correctness of full factorization.
        call get_data(Xdata, X)
        call check(error, maxval(abs(matmul(A%data, Xdata(:, 1:p*kdim)) - matmul(Xdata, H))) < rtol_dp)
        call check_test(error, 'test_block_arnoldi_factorization_cdp', &
                        & 'Arnoldi factorization not correct: A @ X[:,1:p*k] /= X[:,1:p*(k+1)] @ H[1:p*(k+1),1:k]')

        ! Compute Gram matrix associated to the Krylov basis.
        G = zero_cdp
        call innerprod(G, X(1:p*kdim), X(1:p*kdim))

        ! Check orthonormality of the computed basis.
        Id = eye(p*kdim)
        call check(error, norm2(abs(G - Id)) < rtol_dp)
        call check_test(error, 'test_block_arnoldi_factorization_cdp', &
                        & 'Basis X not orthonormal: X[:,1:k].H @ X[:,1:k] /= I')

        return
    end subroutine test_block_arnoldi_factorization_cdp

    subroutine test_krylov_schur_cdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test operator.
        type(linop_cdp), allocatable :: A
        ! Krylov subspace.
        type(vector_cdp), allocatable :: X(:), X0(:)
        ! Krylov subspace dimension.
        integer, parameter :: kdim = 100
        ! Hessenberg matrix.
        complex(dp) :: H(kdim+1, kdim)
        ! Information flag.
        integer :: info

        ! Miscellaneous.
        integer :: n
        complex(dp) :: Xdata(test_size, kdim+1)
        real(dp) :: alpha

        ! Initialize matrix.
        A = linop_cdp() ; call init_rand(A)
        A%data = A%data / norm2(abs(A%data))

        ! Initialize Krylov subspace.
        allocate(X(kdim+1)) ; allocate(X0(1))
        call init_rand(X0) ; call initialize_krylov_subspace(X, X0)
        H = zero_cdp

        ! Arnoldi factorization.
        call arnoldi(A, X, H, info)
        call check_info(info, 'arnoldi', module=this_module, procedure='test_krylov_schur_cdp')

        ! Krylov-Schur condensation.
        call krylov_schur(n, X, H, select_eigs)

        ! Check correctness.
        call get_data(Xdata, X)
        alpha = maxval(abs(matmul(A%data, Xdata(:, :n)) - matmul(Xdata(:, :n+1), H(:n+1, :n))))
        call check(error, alpha < rtol_dp)
        call check_test(error, 'test_krylov_schur_cdp', &
                        & 'Arnoldi factorization not correct: A @ X[:,1:k] /= X[:,1:k+1] @ H[1:k+1,1:k]')

        return
    contains
        function select_eigs(eigvals) result(selected)
            complex(dp), intent(in) :: eigvals(:)
            logical                       :: selected(size(eigvals))
           
            selected = abs(eigvals) > median(abs(eigvals))
            return
        end function select_eigs
    end subroutine test_krylov_schur_cdp


    !------------------------------------------------------------------------------
    !-----     DEFINITION OF THE UNIT-TESTS FOR LANCZOS BIDIAGONALIZATION     -----
    !------------------------------------------------------------------------------

    subroutine collect_lanczos_bidiag_rsp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("Lanczos Bidiagonalization (correctness and left/right orthonormality)",&
                & test_lanczos_bidiag_factorization_rsp) &
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
        real(sp) :: B(kdim+1, kdim)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(sp) :: alpha
        real(sp), dimension(test_size, kdim+1) :: Udata, Vdata
        real(sp), dimension(kdim, kdim) :: G, Id
        type(vector_rsp), allocatable :: X0(:)

        ! Initialize linear operator.
        A = linop_rsp() ; call init_rand(A)

        ! Initialize Krylov subspaces.
        allocate(U(1:kdim+1)) ; allocate(V(1:kdim+1)) ; allocate(X0(1))
        call init_rand(X0) ; call initialize_krylov_subspace(U, X0)
        call initialize_krylov_subspace(V) ; B = zero_rsp

        ! Lanczos bidiagonalization.
        call lanczos_bidiagonalization(A, U, V, B, info)
        call check_info(info, 'lanczos_bidiagonalization', module=this_module, &
                        & procedure='test_lanczos_bidiag_factorization_rsp')

        ! Check correctness.
        call get_data(Udata, U)
        call get_data(Vdata, V)

        alpha = maxval(abs(matmul(A%data, Vdata(:, 1:kdim)) - matmul(Udata, B)))
        call check(error, alpha < rtol_sp)
        call check_test(error, 'test_lanczos_bidiag_factorization_rsp', &
                              & 'Bidiagonal Lanczos factorization not correct: A @ V[:,1:k] /= U[:,1:k+1] @ B[1:k+1,1:k]')

        ! Compute Gram matrix associated to the left Krylov basis.
        G = zero_rsp
        call innerprod(G, U(1:kdim), U(1:kdim))

        ! Check orthonormality of the left basis.
        Id = eye(kdim)
        call check(error, norm2(abs(G - Id)) < rtol_sp)
        call check_test(error, 'test_lanczos_bidiag_factorization_rsp', &
                              & 'Left basis U not orthonormal: U[:,1:k].H @ U[:,1:k] /= I')

        ! Compute Gram matrix associated to the right Krylov basis.
        G = zero_rsp
        call innerprod(G, V(1:kdim), V(1:kdim))

        ! Check orthonormality of the right basis.
        call check(error, norm2(abs(G - Id)) < rtol_sp)
        call check_test(error, 'test_lanczos_bidiag_factorization_rsp', &
                              & 'Right basis V not orthonormal: V[:,1:k].H @ V[:,1:k] /= I')

        return
    end subroutine test_lanczos_bidiag_factorization_rsp

    subroutine collect_lanczos_bidiag_rdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("Lanczos Bidiagonalization (correctness and left/right orthonormality)",&
                & test_lanczos_bidiag_factorization_rdp) &
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
        real(dp) :: B(kdim+1, kdim)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(dp) :: alpha
        real(dp), dimension(test_size, kdim+1) :: Udata, Vdata
        real(dp), dimension(kdim, kdim) :: G, Id
        type(vector_rdp), allocatable :: X0(:)

        ! Initialize linear operator.
        A = linop_rdp() ; call init_rand(A)

        ! Initialize Krylov subspaces.
        allocate(U(1:kdim+1)) ; allocate(V(1:kdim+1)) ; allocate(X0(1))
        call init_rand(X0) ; call initialize_krylov_subspace(U, X0)
        call initialize_krylov_subspace(V) ; B = zero_rdp

        ! Lanczos bidiagonalization.
        call lanczos_bidiagonalization(A, U, V, B, info)
        call check_info(info, 'lanczos_bidiagonalization', module=this_module, &
                        & procedure='test_lanczos_bidiag_factorization_rdp')

        ! Check correctness.
        call get_data(Udata, U)
        call get_data(Vdata, V)

        alpha = maxval(abs(matmul(A%data, Vdata(:, 1:kdim)) - matmul(Udata, B)))
        call check(error, alpha < rtol_dp)
        call check_test(error, 'test_lanczos_bidiag_factorization_rdp', &
                              & 'Bidiagonal Lanczos factorization not correct: A @ V[:,1:k] /= U[:,1:k+1] @ B[1:k+1,1:k]')

        ! Compute Gram matrix associated to the left Krylov basis.
        G = zero_rdp
        call innerprod(G, U(1:kdim), U(1:kdim))

        ! Check orthonormality of the left basis.
        Id = eye(kdim)
        call check(error, norm2(abs(G - Id)) < rtol_dp)
        call check_test(error, 'test_lanczos_bidiag_factorization_rdp', &
                              & 'Left basis U not orthonormal: U[:,1:k].H @ U[:,1:k] /= I')

        ! Compute Gram matrix associated to the right Krylov basis.
        G = zero_rdp
        call innerprod(G, V(1:kdim), V(1:kdim))

        ! Check orthonormality of the right basis.
        call check(error, norm2(abs(G - Id)) < rtol_dp)
        call check_test(error, 'test_lanczos_bidiag_factorization_rdp', &
                              & 'Right basis V not orthonormal: V[:,1:k].H @ V[:,1:k] /= I')

        return
    end subroutine test_lanczos_bidiag_factorization_rdp

    subroutine collect_lanczos_bidiag_csp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("Lanczos Bidiagonalization (correctness and left/right orthonormality)",&
                & test_lanczos_bidiag_factorization_csp) &
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
        complex(sp) :: B(kdim+1, kdim)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(sp) :: alpha
        complex(sp), dimension(test_size, kdim+1) :: Udata, Vdata
        complex(sp), dimension(kdim, kdim) :: G, Id
        type(vector_csp), allocatable :: X0(:)

        ! Initialize linear operator.
        A = linop_csp() ; call init_rand(A)

        ! Initialize Krylov subspaces.
        allocate(U(1:kdim+1)) ; allocate(V(1:kdim+1)) ; allocate(X0(1))
        call init_rand(X0) ; call initialize_krylov_subspace(U, X0)
        call initialize_krylov_subspace(V) ; B = zero_csp

        ! Lanczos bidiagonalization.
        call lanczos_bidiagonalization(A, U, V, B, info)
        call check_info(info, 'lanczos_bidiagonalization', module=this_module, &
                        & procedure='test_lanczos_bidiag_factorization_csp')

        ! Check correctness.
        call get_data(Udata, U)
        call get_data(Vdata, V)

        alpha = maxval(abs(matmul(A%data, Vdata(:, 1:kdim)) - matmul(Udata, B)))
        call check(error, alpha < rtol_sp)
        call check_test(error, 'test_lanczos_bidiag_factorization_csp', &
                              & 'Bidiagonal Lanczos factorization not correct: A @ V[:,1:k] /= U[:,1:k+1] @ B[1:k+1,1:k]')

        ! Compute Gram matrix associated to the left Krylov basis.
        G = zero_csp
        call innerprod(G, U(1:kdim), U(1:kdim))

        ! Check orthonormality of the left basis.
        Id = eye(kdim)
        call check(error, norm2(abs(G - Id)) < rtol_sp)
        call check_test(error, 'test_lanczos_bidiag_factorization_csp', &
                              & 'Left basis U not orthonormal: U[:,1:k].H @ U[:,1:k] /= I')

        ! Compute Gram matrix associated to the right Krylov basis.
        G = zero_csp
        call innerprod(G, V(1:kdim), V(1:kdim))

        ! Check orthonormality of the right basis.
        call check(error, norm2(abs(G - Id)) < rtol_sp)
        call check_test(error, 'test_lanczos_bidiag_factorization_csp', &
                              & 'Right basis V not orthonormal: V[:,1:k].H @ V[:,1:k] /= I')

        return
    end subroutine test_lanczos_bidiag_factorization_csp

    subroutine collect_lanczos_bidiag_cdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("Lanczos Bidiagonalization (correctness and left/right orthonormality)",&
                & test_lanczos_bidiag_factorization_cdp) &
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
        complex(dp) :: B(kdim+1, kdim)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(dp) :: alpha
        complex(dp), dimension(test_size, kdim+1) :: Udata, Vdata
        complex(dp), dimension(kdim, kdim) :: G, Id
        type(vector_cdp), allocatable :: X0(:)

        ! Initialize linear operator.
        A = linop_cdp() ; call init_rand(A)

        ! Initialize Krylov subspaces.
        allocate(U(1:kdim+1)) ; allocate(V(1:kdim+1)) ; allocate(X0(1))
        call init_rand(X0) ; call initialize_krylov_subspace(U, X0)
        call initialize_krylov_subspace(V) ; B = zero_cdp

        ! Lanczos bidiagonalization.
        call lanczos_bidiagonalization(A, U, V, B, info)
        call check_info(info, 'lanczos_bidiagonalization', module=this_module, &
                        & procedure='test_lanczos_bidiag_factorization_cdp')

        ! Check correctness.
        call get_data(Udata, U)
        call get_data(Vdata, V)

        alpha = maxval(abs(matmul(A%data, Vdata(:, 1:kdim)) - matmul(Udata, B)))
        call check(error, alpha < rtol_dp)
        call check_test(error, 'test_lanczos_bidiag_factorization_cdp', &
                              & 'Bidiagonal Lanczos factorization not correct: A @ V[:,1:k] /= U[:,1:k+1] @ B[1:k+1,1:k]')

        ! Compute Gram matrix associated to the left Krylov basis.
        G = zero_cdp
        call innerprod(G, U(1:kdim), U(1:kdim))

        ! Check orthonormality of the left basis.
        Id = eye(kdim)
        call check(error, norm2(abs(G - Id)) < rtol_dp)
        call check_test(error, 'test_lanczos_bidiag_factorization_cdp', &
                              & 'Left basis U not orthonormal: U[:,1:k].H @ U[:,1:k] /= I')

        ! Compute Gram matrix associated to the right Krylov basis.
        G = zero_cdp
        call innerprod(G, V(1:kdim), V(1:kdim))

        ! Check orthonormality of the right basis.
        call check(error, norm2(abs(G - Id)) < rtol_dp)
        call check_test(error, 'test_lanczos_bidiag_factorization_cdp', &
                              & 'Right basis V not orthonormal: V[:,1:k].H @ V[:,1:k] /= I')

        return
    end subroutine test_lanczos_bidiag_factorization_cdp


    !-------------------------------------------------------------------------------
    !-----     DEFINITION OF THE UNIT TESTS FOR LANCZOS TRIDIAGONALIZATION     -----
    !-------------------------------------------------------------------------------

    subroutine collect_lanczos_tridiag_rsp_testsuite(testsuite)
        ! Collection of unit tests.
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
             new_unittest("Lanczos Tridiagonalization (correctness & orthonormality)", test_lanczos_tridiag_factorization_rsp) &
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
        real(sp) :: T(kdim+1, kdim)
        ! Information flag.
        integer :: info

        ! Internal variables.
        real(sp) :: Xdata(test_size, kdim+1)
        real(sp) :: alpha
        class(vector_rsp), allocatable :: X0(:)
        real(sp), dimension(kdim, kdim) :: G, Id

        ! Initialize tridiagonal matrix.
        T = zero_rsp

        ! Initialize operator.
        A = spd_linop_rsp()
        call init_rand(A)

        ! Initialize Krylov subspace.
        allocate(X(kdim+1)) ; allocate(X0(1))
        call init_rand(X0) ; call initialize_krylov_subspace(X, X0)

        ! Lanczos factorization.
        call lanczos_tridiagonalization(A, X, T, info)
        call check_info(info, 'lanczos_tridiagonalization', module=this_module, & 
                        & procedure='test_lanczos_tridiag_factorization_rsp')

        ! Check correctness.
        call get_data(Xdata, X)

        ! Infinity-norm check.
        alpha = maxval(abs(matmul(A%data, Xdata(:, 1:kdim)) - matmul(Xdata, T)))
        call check(error, alpha < rtol_sp)
        call check_test(error, 'test_lanczos_tridiag_factorization_rsp', &
                                 & 'Tridiagonal Lanczos factorization not correct: A @ X[:,1:k] /= X[:,1:k+1] @ T[1:k+1,1:k]')

        ! Compute Gram matrix associated to the right Krylov basis.
        G = zero_rsp
        call innerprod(G, X(1:kdim), X(1:kdim))

        ! Check orthonormality of the Krylov basis.
        Id = eye(kdim)
        call check(error, norm2(abs(G - Id)) < rtol_sp)
        call check_test(error, 'test_lanczos_tridiag_factorization_rsp', &
                                 &'Basis X not orthonormal: X[:,1:k].H @ X[:,1:k] /= I')

        return
    end subroutine test_lanczos_tridiag_factorization_rsp

    subroutine collect_lanczos_tridiag_rdp_testsuite(testsuite)
        ! Collection of unit tests.
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
             new_unittest("Lanczos Tridiagonalization (correctness & orthonormality)", test_lanczos_tridiag_factorization_rdp) &
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
        real(dp) :: T(kdim+1, kdim)
        ! Information flag.
        integer :: info

        ! Internal variables.
        real(dp) :: Xdata(test_size, kdim+1)
        real(dp) :: alpha
        class(vector_rdp), allocatable :: X0(:)
        real(dp), dimension(kdim, kdim) :: G, Id

        ! Initialize tridiagonal matrix.
        T = zero_rdp

        ! Initialize operator.
        A = spd_linop_rdp()
        call init_rand(A)

        ! Initialize Krylov subspace.
        allocate(X(kdim+1)) ; allocate(X0(1))
        call init_rand(X0) ; call initialize_krylov_subspace(X, X0)

        ! Lanczos factorization.
        call lanczos_tridiagonalization(A, X, T, info)
        call check_info(info, 'lanczos_tridiagonalization', module=this_module, & 
                        & procedure='test_lanczos_tridiag_factorization_rdp')

        ! Check correctness.
        call get_data(Xdata, X)

        ! Infinity-norm check.
        alpha = maxval(abs(matmul(A%data, Xdata(:, 1:kdim)) - matmul(Xdata, T)))
        call check(error, alpha < rtol_dp)
        call check_test(error, 'test_lanczos_tridiag_factorization_rdp', &
                                 & 'Tridiagonal Lanczos factorization not correct: A @ X[:,1:k] /= X[:,1:k+1] @ T[1:k+1,1:k]')

        ! Compute Gram matrix associated to the right Krylov basis.
        G = zero_rdp
        call innerprod(G, X(1:kdim), X(1:kdim))

        ! Check orthonormality of the Krylov basis.
        Id = eye(kdim)
        call check(error, norm2(abs(G - Id)) < rtol_dp)
        call check_test(error, 'test_lanczos_tridiag_factorization_rdp', &
                                 &'Basis X not orthonormal: X[:,1:k].H @ X[:,1:k] /= I')

        return
    end subroutine test_lanczos_tridiag_factorization_rdp

    subroutine collect_lanczos_tridiag_csp_testsuite(testsuite)
        ! Collection of unit tests.
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
             new_unittest("Lanczos Tridiagonalization (correctness & orthonormality)", test_lanczos_tridiag_factorization_csp) &
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
        complex(sp) :: T(kdim+1, kdim)
        ! Information flag.
        integer :: info

        ! Internal variables.
        complex(sp) :: Xdata(test_size, kdim+1)
        real(sp) :: alpha
        class(vector_csp), allocatable :: X0(:)
        complex(sp), dimension(kdim, kdim) :: G, Id

        ! Initialize tridiagonal matrix.
        T = zero_csp

        ! Initialize operator.
        A = hermitian_linop_csp()
        call init_rand(A)

        ! Initialize Krylov subspace.
        allocate(X(kdim+1)) ; allocate(X0(1))
        call init_rand(X0) ; call initialize_krylov_subspace(X, X0)

        ! Lanczos factorization.
        call lanczos_tridiagonalization(A, X, T, info)
        call check_info(info, 'lanczos_tridiagonalization', module=this_module, & 
                        & procedure='test_lanczos_tridiag_factorization_csp')

        ! Check correctness.
        call get_data(Xdata, X)

        ! Infinity-norm check.
        alpha = maxval(abs(matmul(A%data, Xdata(:, 1:kdim)) - matmul(Xdata, T)))
        call check(error, alpha < rtol_sp)
        call check_test(error, 'test_lanczos_tridiag_factorization_csp', &
                                 & 'Tridiagonal Lanczos factorization not correct: A @ X[:,1:k] /= X[:,1:k+1] @ T[1:k+1,1:k]')

        ! Compute Gram matrix associated to the right Krylov basis.
        G = zero_csp
        call innerprod(G, X(1:kdim), X(1:kdim))

        ! Check orthonormality of the Krylov basis.
        Id = eye(kdim)
        call check(error, norm2(abs(G - Id)) < rtol_sp)
        call check_test(error, 'test_lanczos_tridiag_factorization_csp', &
                                 &'Basis X not orthonormal: X[:,1:k].H @ X[:,1:k] /= I')

        return
    end subroutine test_lanczos_tridiag_factorization_csp

    subroutine collect_lanczos_tridiag_cdp_testsuite(testsuite)
        ! Collection of unit tests.
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
             new_unittest("Lanczos Tridiagonalization (correctness & orthonormality)", test_lanczos_tridiag_factorization_cdp) &
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
        complex(dp) :: T(kdim+1, kdim)
        ! Information flag.
        integer :: info

        ! Internal variables.
        complex(dp) :: Xdata(test_size, kdim+1)
        real(dp) :: alpha
        class(vector_cdp), allocatable :: X0(:)
        complex(dp), dimension(kdim, kdim) :: G, Id

        ! Initialize tridiagonal matrix.
        T = zero_cdp

        ! Initialize operator.
        A = hermitian_linop_cdp()
        call init_rand(A)

        ! Initialize Krylov subspace.
        allocate(X(kdim+1)) ; allocate(X0(1))
        call init_rand(X0) ; call initialize_krylov_subspace(X, X0)

        ! Lanczos factorization.
        call lanczos_tridiagonalization(A, X, T, info)
        call check_info(info, 'lanczos_tridiagonalization', module=this_module, & 
                        & procedure='test_lanczos_tridiag_factorization_cdp')

        ! Check correctness.
        call get_data(Xdata, X)

        ! Infinity-norm check.
        alpha = maxval(abs(matmul(A%data, Xdata(:, 1:kdim)) - matmul(Xdata, T)))
        call check(error, alpha < rtol_dp)
        call check_test(error, 'test_lanczos_tridiag_factorization_cdp', &
                                 & 'Tridiagonal Lanczos factorization not correct: A @ X[:,1:k] /= X[:,1:k+1] @ T[1:k+1,1:k]')

        ! Compute Gram matrix associated to the right Krylov basis.
        G = zero_cdp
        call innerprod(G, X(1:kdim), X(1:kdim))

        ! Check orthonormality of the Krylov basis.
        Id = eye(kdim)
        call check(error, norm2(abs(G - Id)) < rtol_dp)
        call check_test(error, 'test_lanczos_tridiag_factorization_cdp', &
                                 &'Basis X not orthonormal: X[:,1:k].H @ X[:,1:k] /= I')

        return
    end subroutine test_lanczos_tridiag_factorization_cdp

end module TestKrylov
