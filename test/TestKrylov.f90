module TestKrylov
    ! Fortran Standard Library.
    use iso_fortran_env
    use stdlib_math, only: is_close, all_close
    use stdlib_linalg, only: eye
    use stdlib_stats, only: median

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
                        new_unittest("Pivoting QR (exact rank def.)", test_pivoting_qr_exact_rank_deficiency_rsp) &
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
        real(sp) :: R(kdim, kdim) = 0.0_sp
        ! Information flag.
        integer :: info
        ! Miscellaneous.
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
        real(sp) :: Id(kdim, kdim)
        ! Information flag.
        integer :: info
       
        ! Miscellaneous.
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
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test Vectors.
        integer, parameter :: kdim = test_size
        type(vector_rdp), allocatable :: A(:)
        ! Upper triangular matrix.
        real(dp) :: R(kdim, kdim) = 0.0_dp
        ! Information flag.
        integer :: info
        ! Miscellaneous.
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
        real(dp) :: Id(kdim, kdim)
        ! Information flag.
        integer :: info
       
        ! Miscellaneous.
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
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test Vectors.
        integer, parameter :: kdim = test_size
        type(vector_csp), allocatable :: A(:)
        ! Upper triangular matrix.
        complex(sp) :: R(kdim, kdim) = 0.0_sp
        ! Information flag.
        integer :: info
        ! Miscellaneous.
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
        complex(sp) :: Id(kdim, kdim)
        ! Information flag.
        integer :: info
       
        ! Miscellaneous.
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
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test Vectors.
        integer, parameter :: kdim = test_size
        type(vector_cdp), allocatable :: A(:)
        ! Upper triangular matrix.
        complex(dp) :: R(kdim, kdim) = 0.0_dp
        ! Information flag.
        integer :: info
        ! Miscellaneous.
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
        complex(dp) :: Id(kdim, kdim)
        ! Information flag.
        integer :: info
       
        ! Miscellaneous.
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

    
    !--------------------------------------------------------------
    !-----     DEFINITIONS OF THE UNIT-TESTS FOR ARNOLDI      -----
    !--------------------------------------------------------------

    subroutine collect_arnoldi_rsp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("Arnoldi factorization", test_arnoldi_factorization_rsp), &
            new_unittest("Arnoldi orthogonality", test_arnoldi_basis_orthogonality_rsp), &
            new_unittest("Block Arnoldi factorization", test_block_arnoldi_factorization_rsp), &
            new_unittest("Block Arnoldi orthogonality", test_block_arnoldi_basis_orthogonality_rsp), &
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
        real(sp) :: H(kdim+1, kdim)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(sp) :: Xdata(test_size, kdim+1), alpha
    
        ! Initialize linear operator.
        A = linop_rsp() ; call init_rand(A)
        ! Initialize Krylov subspace.
        allocate(X(1:kdim+1)) ;
        call initialize_krylov_subspace(X)
        call X(1)%rand(ifnorm=.true.) ; alpha = X(1)%norm()
        call X(1)%scal(1.0_sp / alpha)
        H = 0.0_sp
        ! Arnoldi factorization.
        call arnoldi(A, X, H, info)

        ! Check correctness of full factorization.
        call get_data(Xdata, X)

        call check(error, maxval(abs(matmul(A%data, Xdata(:, 1:kdim)) - matmul(Xdata, H))) < rtol_sp)

        return
    end subroutine test_arnoldi_factorization_rsp

    subroutine test_arnoldi_basis_orthogonality_rsp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test matrix.
        type(linop_rsp), allocatable :: A
        ! Krylov subspace.
        type(vector_rsp), dimension(:), allocatable :: X
        type(vector_rsp), dimension(:), allocatable :: X0
        ! Krylov subspace dimension.
        integer, parameter :: kdim = test_size
        ! Hessenberg matrix.
        real(sp), dimension(kdim + 1, kdim) :: H
        ! Information flag.
        integer :: info
        ! Misc.
        real(sp), dimension(kdim, kdim) :: G, Id
        integer :: i, j

        ! Initialize random matrix.
        A = linop_rsp(); call init_rand(A)
        ! Initialize Krylov subspace.
        allocate (X(1:kdim + 1)); allocate (X0(1))
        call init_rand(X0) ; call initialize_krylov_subspace(X, X0)
        H = 0.0_sp

        ! Arnoldi factorization.
        call arnoldi(A, X, H, info)

        ! Compute Gram matrix associated to the Krylov basis.
        G = 0.0_sp
        do j = 1, kdim
            do i = 1, kdim
                G(i, j) = X(i)%dot(X(j))
            enddo
        enddo
        ! call innerprod_matrix(G, X(1:kdim), X(1:kdim))

        ! Check result.
        Id = eye(kdim)
        call check(error, norm2(abs(G - Id)) < rtol_sp)

        return
    end subroutine test_arnoldi_basis_orthogonality_rsp

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
        real(sp) :: Xdata(test_size, p*(kdim+1)), G(p*(kdim+1), p*(kdim+1))
        integer :: i, j

        ! Initialize linear operator.
        A = linop_rsp() ; call init_rand(A)

        ! Initialize Krylov subspace.
        allocate(X(1:p*(kdim+1))) ; allocate(X0(1:p))
        call init_rand(X0) ; call initialize_krylov_subspace(X, X0)
        H = 0.0_sp

        ! Arnoldi factorization.
        call arnoldi(A, X, H, info, blksize=p)

        G = 0.0_sp
        do j = 1, size(G, 2)
            do i = 1, size(G, 1)
                G(i, j) = X(i)%dot(X(j))
            enddo
        enddo

        ! Check correctness.
        call get_data(Xdata, X)
        
        call check(error, maxval(abs(matmul(A%data, Xdata(:, 1:p*kdim)) - matmul(Xdata, H))) < rtol_sp)

        return
    end subroutine test_block_arnoldi_factorization_rsp

    subroutine test_block_arnoldi_basis_orthogonality_rsp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test matrix.
        class(linop_rsp), allocatable :: A
        ! Krylov subspace.
        type(vector_rsp), dimension(:), allocatable :: X
        type(vector_rsp), dimension(:), allocatable :: X0
        ! Krylov subspace dimension.
        integer, parameter :: p = 2
        integer, parameter :: kdim = test_size/2
        ! Hessenberg matrix.
        real(sp), dimension(p*(kdim + 1), p*kdim) :: H
        ! Information flag.
        integer :: info
        ! Misc.
        real(sp), dimension(p*kdim, p*kdim) :: G, Id
        integer :: i, j

        ! Initialize random matrix.
        A = linop_rsp(); call init_rand(A)

        ! Initialize Krylov subspace.
        allocate (X(1:p*(kdim + 1))); allocate (X0(1:p)); 
        call init_rand(X0)
        call initialize_krylov_subspace(X, X0)
        H = 0.0_sp

        ! Arnoldi factorization.
        call arnoldi(A, X, H, info, blksize=p)

        ! Compute Gram matrix associated to the Krylov basis.
        G = 0.0_sp
        do j = 1, size(G, 2)
            do i = 1, size(G, 1)
                G(i, j) = X(i)%dot(X(j))
            enddo
        enddo

        ! Check result.
        Id = eye(p*kdim)
        call check(error, norm2(abs(G - Id)) < rtol_sp)

        return
    end subroutine test_block_arnoldi_basis_orthogonality_rsp

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
        H = 0.0_sp

        ! Arnoldi factorization.
        call arnoldi(A, X, H, info)

        ! Krylov-Schur condensation.
        call krylov_schur(n, X, H, select_eigs)

        ! Check correctness.
        call get_data(Xdata, X)
        alpha = maxval(abs(matmul(A%data, Xdata(:, :n)) - matmul(Xdata(:, :n+1), H(:n+1, :n))))
        call check(error, alpha < rtol_sp)

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
            new_unittest("Arnoldi factorization", test_arnoldi_factorization_rdp), &
            new_unittest("Arnoldi orthogonality", test_arnoldi_basis_orthogonality_rdp), &
            new_unittest("Block Arnoldi factorization", test_block_arnoldi_factorization_rdp), &
            new_unittest("Block Arnoldi orthogonality", test_block_arnoldi_basis_orthogonality_rdp), &
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
        real(dp) :: H(kdim+1, kdim)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(dp) :: Xdata(test_size, kdim+1), alpha
    
        ! Initialize linear operator.
        A = linop_rdp() ; call init_rand(A)
        ! Initialize Krylov subspace.
        allocate(X(1:kdim+1)) ;
        call initialize_krylov_subspace(X)
        call X(1)%rand(ifnorm=.true.) ; alpha = X(1)%norm()
        call X(1)%scal(1.0_dp / alpha)
        H = 0.0_dp
        ! Arnoldi factorization.
        call arnoldi(A, X, H, info)

        ! Check correctness of full factorization.
        call get_data(Xdata, X)

        call check(error, maxval(abs(matmul(A%data, Xdata(:, 1:kdim)) - matmul(Xdata, H))) < rtol_dp)

        return
    end subroutine test_arnoldi_factorization_rdp

    subroutine test_arnoldi_basis_orthogonality_rdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test matrix.
        type(linop_rdp), allocatable :: A
        ! Krylov subspace.
        type(vector_rdp), dimension(:), allocatable :: X
        type(vector_rdp), dimension(:), allocatable :: X0
        ! Krylov subspace dimension.
        integer, parameter :: kdim = test_size
        ! Hessenberg matrix.
        real(dp), dimension(kdim + 1, kdim) :: H
        ! Information flag.
        integer :: info
        ! Misc.
        real(dp), dimension(kdim, kdim) :: G, Id
        integer :: i, j

        ! Initialize random matrix.
        A = linop_rdp(); call init_rand(A)
        ! Initialize Krylov subspace.
        allocate (X(1:kdim + 1)); allocate (X0(1))
        call init_rand(X0) ; call initialize_krylov_subspace(X, X0)
        H = 0.0_dp

        ! Arnoldi factorization.
        call arnoldi(A, X, H, info)

        ! Compute Gram matrix associated to the Krylov basis.
        G = 0.0_dp
        do j = 1, kdim
            do i = 1, kdim
                G(i, j) = X(i)%dot(X(j))
            enddo
        enddo
        ! call innerprod_matrix(G, X(1:kdim), X(1:kdim))

        ! Check result.
        Id = eye(kdim)
        call check(error, norm2(abs(G - Id)) < rtol_dp)

        return
    end subroutine test_arnoldi_basis_orthogonality_rdp

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
        real(dp) :: Xdata(test_size, p*(kdim+1)), G(p*(kdim+1), p*(kdim+1))
        integer :: i, j

        ! Initialize linear operator.
        A = linop_rdp() ; call init_rand(A)

        ! Initialize Krylov subspace.
        allocate(X(1:p*(kdim+1))) ; allocate(X0(1:p))
        call init_rand(X0) ; call initialize_krylov_subspace(X, X0)
        H = 0.0_dp

        ! Arnoldi factorization.
        call arnoldi(A, X, H, info, blksize=p)

        G = 0.0_dp
        do j = 1, size(G, 2)
            do i = 1, size(G, 1)
                G(i, j) = X(i)%dot(X(j))
            enddo
        enddo

        ! Check correctness.
        call get_data(Xdata, X)
        
        call check(error, maxval(abs(matmul(A%data, Xdata(:, 1:p*kdim)) - matmul(Xdata, H))) < rtol_dp)

        return
    end subroutine test_block_arnoldi_factorization_rdp

    subroutine test_block_arnoldi_basis_orthogonality_rdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test matrix.
        class(linop_rdp), allocatable :: A
        ! Krylov subspace.
        type(vector_rdp), dimension(:), allocatable :: X
        type(vector_rdp), dimension(:), allocatable :: X0
        ! Krylov subspace dimension.
        integer, parameter :: p = 2
        integer, parameter :: kdim = test_size/2
        ! Hessenberg matrix.
        real(dp), dimension(p*(kdim + 1), p*kdim) :: H
        ! Information flag.
        integer :: info
        ! Misc.
        real(dp), dimension(p*kdim, p*kdim) :: G, Id
        integer :: i, j

        ! Initialize random matrix.
        A = linop_rdp(); call init_rand(A)

        ! Initialize Krylov subspace.
        allocate (X(1:p*(kdim + 1))); allocate (X0(1:p)); 
        call init_rand(X0)
        call initialize_krylov_subspace(X, X0)
        H = 0.0_dp

        ! Arnoldi factorization.
        call arnoldi(A, X, H, info, blksize=p)

        ! Compute Gram matrix associated to the Krylov basis.
        G = 0.0_dp
        do j = 1, size(G, 2)
            do i = 1, size(G, 1)
                G(i, j) = X(i)%dot(X(j))
            enddo
        enddo

        ! Check result.
        Id = eye(p*kdim)
        call check(error, norm2(abs(G - Id)) < rtol_dp)

        return
    end subroutine test_block_arnoldi_basis_orthogonality_rdp

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
        H = 0.0_dp

        ! Arnoldi factorization.
        call arnoldi(A, X, H, info)

        ! Krylov-Schur condensation.
        call krylov_schur(n, X, H, select_eigs)

        ! Check correctness.
        call get_data(Xdata, X)
        alpha = maxval(abs(matmul(A%data, Xdata(:, :n)) - matmul(Xdata(:, :n+1), H(:n+1, :n))))
        call check(error, alpha < rtol_dp)

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
            new_unittest("Arnoldi factorization", test_arnoldi_factorization_csp), &
            new_unittest("Arnoldi orthogonality", test_arnoldi_basis_orthogonality_csp), &
            new_unittest("Block Arnoldi factorization", test_block_arnoldi_factorization_csp), &
            new_unittest("Block Arnoldi orthogonality", test_block_arnoldi_basis_orthogonality_csp), &
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
        complex(sp) :: H(kdim+1, kdim)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        complex(sp) :: Xdata(test_size, kdim+1), alpha
    
        ! Initialize linear operator.
        A = linop_csp() ; call init_rand(A)
        ! Initialize Krylov subspace.
        allocate(X(1:kdim+1)) ;
        call initialize_krylov_subspace(X)
        call X(1)%rand(ifnorm=.true.) ; alpha = X(1)%norm()
        call X(1)%scal(1.0_sp / alpha)
        H = 0.0_sp
        ! Arnoldi factorization.
        call arnoldi(A, X, H, info)

        ! Check correctness of full factorization.
        call get_data(Xdata, X)

        call check(error, maxval(abs(matmul(A%data, Xdata(:, 1:kdim)) - matmul(Xdata, H))) < rtol_sp)

        return
    end subroutine test_arnoldi_factorization_csp

    subroutine test_arnoldi_basis_orthogonality_csp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test matrix.
        type(linop_csp), allocatable :: A
        ! Krylov subspace.
        type(vector_csp), dimension(:), allocatable :: X
        type(vector_csp), dimension(:), allocatable :: X0
        ! Krylov subspace dimension.
        integer, parameter :: kdim = test_size
        ! Hessenberg matrix.
        complex(sp), dimension(kdim + 1, kdim) :: H
        ! Information flag.
        integer :: info
        ! Misc.
        complex(sp), dimension(kdim, kdim) :: G, Id
        integer :: i, j

        ! Initialize random matrix.
        A = linop_csp(); call init_rand(A)
        ! Initialize Krylov subspace.
        allocate (X(1:kdim + 1)); allocate (X0(1))
        call init_rand(X0) ; call initialize_krylov_subspace(X, X0)
        H = 0.0_sp

        ! Arnoldi factorization.
        call arnoldi(A, X, H, info)

        ! Compute Gram matrix associated to the Krylov basis.
        G = 0.0_sp
        do j = 1, kdim
            do i = 1, kdim
                G(i, j) = X(i)%dot(X(j))
            enddo
        enddo
        ! call innerprod_matrix(G, X(1:kdim), X(1:kdim))

        ! Check result.
        Id = eye(kdim)
        call check(error, norm2(abs(G - Id)) < rtol_sp)

        return
    end subroutine test_arnoldi_basis_orthogonality_csp

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
        complex(sp) :: Xdata(test_size, p*(kdim+1)), G(p*(kdim+1), p*(kdim+1))
        integer :: i, j

        ! Initialize linear operator.
        A = linop_csp() ; call init_rand(A)

        ! Initialize Krylov subspace.
        allocate(X(1:p*(kdim+1))) ; allocate(X0(1:p))
        call init_rand(X0) ; call initialize_krylov_subspace(X, X0)
        H = 0.0_sp

        ! Arnoldi factorization.
        call arnoldi(A, X, H, info, blksize=p)

        G = 0.0_sp
        do j = 1, size(G, 2)
            do i = 1, size(G, 1)
                G(i, j) = X(i)%dot(X(j))
            enddo
        enddo

        ! Check correctness.
        call get_data(Xdata, X)
        
        call check(error, maxval(abs(matmul(A%data, Xdata(:, 1:p*kdim)) - matmul(Xdata, H))) < rtol_sp)

        return
    end subroutine test_block_arnoldi_factorization_csp

    subroutine test_block_arnoldi_basis_orthogonality_csp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test matrix.
        class(linop_csp), allocatable :: A
        ! Krylov subspace.
        type(vector_csp), dimension(:), allocatable :: X
        type(vector_csp), dimension(:), allocatable :: X0
        ! Krylov subspace dimension.
        integer, parameter :: p = 2
        integer, parameter :: kdim = test_size/2
        ! Hessenberg matrix.
        complex(sp), dimension(p*(kdim + 1), p*kdim) :: H
        ! Information flag.
        integer :: info
        ! Misc.
        complex(sp), dimension(p*kdim, p*kdim) :: G, Id
        integer :: i, j

        ! Initialize random matrix.
        A = linop_csp(); call init_rand(A)

        ! Initialize Krylov subspace.
        allocate (X(1:p*(kdim + 1))); allocate (X0(1:p)); 
        call init_rand(X0)
        call initialize_krylov_subspace(X, X0)
        H = 0.0_sp

        ! Arnoldi factorization.
        call arnoldi(A, X, H, info, blksize=p)

        ! Compute Gram matrix associated to the Krylov basis.
        G = 0.0_sp
        do j = 1, size(G, 2)
            do i = 1, size(G, 1)
                G(i, j) = X(i)%dot(X(j))
            enddo
        enddo

        ! Check result.
        Id = eye(p*kdim)
        call check(error, norm2(abs(G - Id)) < rtol_sp)

        return
    end subroutine test_block_arnoldi_basis_orthogonality_csp

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
        H = 0.0_sp

        ! Arnoldi factorization.
        call arnoldi(A, X, H, info)

        ! Krylov-Schur condensation.
        call krylov_schur(n, X, H, select_eigs)

        ! Check correctness.
        call get_data(Xdata, X)
        alpha = maxval(abs(matmul(A%data, Xdata(:, :n)) - matmul(Xdata(:, :n+1), H(:n+1, :n))))
        call check(error, alpha < rtol_sp)

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
            new_unittest("Arnoldi factorization", test_arnoldi_factorization_cdp), &
            new_unittest("Arnoldi orthogonality", test_arnoldi_basis_orthogonality_cdp), &
            new_unittest("Block Arnoldi factorization", test_block_arnoldi_factorization_cdp), &
            new_unittest("Block Arnoldi orthogonality", test_block_arnoldi_basis_orthogonality_cdp), &
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
        complex(dp) :: H(kdim+1, kdim)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        complex(dp) :: Xdata(test_size, kdim+1), alpha
    
        ! Initialize linear operator.
        A = linop_cdp() ; call init_rand(A)
        ! Initialize Krylov subspace.
        allocate(X(1:kdim+1)) ;
        call initialize_krylov_subspace(X)
        call X(1)%rand(ifnorm=.true.) ; alpha = X(1)%norm()
        call X(1)%scal(1.0_dp / alpha)
        H = 0.0_dp
        ! Arnoldi factorization.
        call arnoldi(A, X, H, info)

        ! Check correctness of full factorization.
        call get_data(Xdata, X)

        call check(error, maxval(abs(matmul(A%data, Xdata(:, 1:kdim)) - matmul(Xdata, H))) < rtol_dp)

        return
    end subroutine test_arnoldi_factorization_cdp

    subroutine test_arnoldi_basis_orthogonality_cdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test matrix.
        type(linop_cdp), allocatable :: A
        ! Krylov subspace.
        type(vector_cdp), dimension(:), allocatable :: X
        type(vector_cdp), dimension(:), allocatable :: X0
        ! Krylov subspace dimension.
        integer, parameter :: kdim = test_size
        ! Hessenberg matrix.
        complex(dp), dimension(kdim + 1, kdim) :: H
        ! Information flag.
        integer :: info
        ! Misc.
        complex(dp), dimension(kdim, kdim) :: G, Id
        integer :: i, j

        ! Initialize random matrix.
        A = linop_cdp(); call init_rand(A)
        ! Initialize Krylov subspace.
        allocate (X(1:kdim + 1)); allocate (X0(1))
        call init_rand(X0) ; call initialize_krylov_subspace(X, X0)
        H = 0.0_dp

        ! Arnoldi factorization.
        call arnoldi(A, X, H, info)

        ! Compute Gram matrix associated to the Krylov basis.
        G = 0.0_dp
        do j = 1, kdim
            do i = 1, kdim
                G(i, j) = X(i)%dot(X(j))
            enddo
        enddo
        ! call innerprod_matrix(G, X(1:kdim), X(1:kdim))

        ! Check result.
        Id = eye(kdim)
        call check(error, norm2(abs(G - Id)) < rtol_dp)

        return
    end subroutine test_arnoldi_basis_orthogonality_cdp

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
        complex(dp) :: Xdata(test_size, p*(kdim+1)), G(p*(kdim+1), p*(kdim+1))
        integer :: i, j

        ! Initialize linear operator.
        A = linop_cdp() ; call init_rand(A)

        ! Initialize Krylov subspace.
        allocate(X(1:p*(kdim+1))) ; allocate(X0(1:p))
        call init_rand(X0) ; call initialize_krylov_subspace(X, X0)
        H = 0.0_dp

        ! Arnoldi factorization.
        call arnoldi(A, X, H, info, blksize=p)

        G = 0.0_dp
        do j = 1, size(G, 2)
            do i = 1, size(G, 1)
                G(i, j) = X(i)%dot(X(j))
            enddo
        enddo

        ! Check correctness.
        call get_data(Xdata, X)
        
        call check(error, maxval(abs(matmul(A%data, Xdata(:, 1:p*kdim)) - matmul(Xdata, H))) < rtol_dp)

        return
    end subroutine test_block_arnoldi_factorization_cdp

    subroutine test_block_arnoldi_basis_orthogonality_cdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test matrix.
        class(linop_cdp), allocatable :: A
        ! Krylov subspace.
        type(vector_cdp), dimension(:), allocatable :: X
        type(vector_cdp), dimension(:), allocatable :: X0
        ! Krylov subspace dimension.
        integer, parameter :: p = 2
        integer, parameter :: kdim = test_size/2
        ! Hessenberg matrix.
        complex(dp), dimension(p*(kdim + 1), p*kdim) :: H
        ! Information flag.
        integer :: info
        ! Misc.
        complex(dp), dimension(p*kdim, p*kdim) :: G, Id
        integer :: i, j

        ! Initialize random matrix.
        A = linop_cdp(); call init_rand(A)

        ! Initialize Krylov subspace.
        allocate (X(1:p*(kdim + 1))); allocate (X0(1:p)); 
        call init_rand(X0)
        call initialize_krylov_subspace(X, X0)
        H = 0.0_dp

        ! Arnoldi factorization.
        call arnoldi(A, X, H, info, blksize=p)

        ! Compute Gram matrix associated to the Krylov basis.
        G = 0.0_dp
        do j = 1, size(G, 2)
            do i = 1, size(G, 1)
                G(i, j) = X(i)%dot(X(j))
            enddo
        enddo

        ! Check result.
        Id = eye(p*kdim)
        call check(error, norm2(abs(G - Id)) < rtol_dp)

        return
    end subroutine test_block_arnoldi_basis_orthogonality_cdp

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
        H = 0.0_dp

        ! Arnoldi factorization.
        call arnoldi(A, X, H, info)

        ! Krylov-Schur condensation.
        call krylov_schur(n, X, H, select_eigs)

        ! Check correctness.
        call get_data(Xdata, X)
        alpha = maxval(abs(matmul(A%data, Xdata(:, :n)) - matmul(Xdata(:, :n+1), H(:n+1, :n))))
        call check(error, alpha < rtol_dp)

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
            new_unittest("Lanczos Bidiagonalization", test_lanczos_bidiag_factorization_rsp), &
            new_unittest("Lanczos left orthogonality", test_lanczos_bidiag_left_orthogonality_rsp), &
            new_unittest("Lanczos right orthogonality", test_lanczos_bidiag_right_orthogonality_rsp) &
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
        real(sp) :: Udata(test_size, kdim+1), Vdata(test_size, kdim+1)
        type(vector_rsp), allocatable :: X0(:)

        ! Initialize linear operator.
        A = linop_rsp() ; call init_rand(A)

        ! Initialize Krylov subspaces.
        allocate(U(1:kdim+1)) ; allocate(V(1:kdim+1)) ; allocate(X0(1))
        call init_rand(X0) ; call initialize_krylov_subspace(U, X0)
        call initialize_krylov_subspace(V) ; B = 0.0_sp

        ! Lanczos bidiagonalization.
        call lanczos_bidiagonalization(A, U, V, B, info)

        ! Check correctness.
        call get_data(Udata, U)
        call get_data(Vdata, V)

        alpha = maxval(abs(matmul(A%data, Vdata(:, 1:kdim)) - matmul(Udata, B)))
        call check(error, alpha < rtol_sp)

        return
    end subroutine test_lanczos_bidiag_factorization_rsp

    subroutine test_lanczos_bidiag_left_orthogonality_rsp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test Linear Operator
        type(linop_rsp), allocatable :: A
        ! Left and right Krylov bases.
        type(vector_rsp), allocatable :: U(:), V(:)
        ! Krylov subspace dimension.
        integer, parameter :: kdim = 3!test_size
        ! Bidiagonal matrix.
        real(sp) :: B(kdim+1, kdim)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(sp) :: Id(kdim, kdim), G(kdim, kdim)
        type(vector_rsp), allocatable :: X0(:)
        integer :: i, j

        ! Initialize linear operator.
        A = linop_rsp() ; call init_rand(A)

        ! Initialize Krylov subspaces.
        allocate(U(1:kdim+1)) ; allocate(V(1:kdim+1)) ; allocate(X0(1))
        call init_rand(X0) ; call initialize_krylov_subspace(U, X0)
        call initialize_krylov_subspace(V) ; B = 0.0_sp

        ! Lanczos bidiagonalization.
        call lanczos_bidiagonalization(A, U, V, B, info)

        ! Check correctness.
        Id = eye(kdim)
        G = 0.0_sp
        do j = 1, kdim
            do i = 1, kdim
                G(i, j) = U(i)%dot(U(j))
            enddo
        enddo

        call check(error, maxval(abs(Id - G)) < rtol_sp)

        return
    end subroutine test_lanczos_bidiag_left_orthogonality_rsp

    subroutine test_lanczos_bidiag_right_orthogonality_rsp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test Linear Operator
        type(linop_rsp), allocatable :: A
        ! Left and right Krylov bases.
        type(vector_rsp), allocatable :: U(:), V(:)
        ! Krylov subspace dimension.
        integer, parameter :: kdim = 3!test_size
        ! Bidiagonal matrix.
        real(sp) :: B(kdim+1, kdim)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(sp) :: Id(kdim, kdim), G(kdim, kdim)
        type(vector_rsp), allocatable :: X0(:)
        integer :: i, j

        ! Initialize linear operator.
        A = linop_rsp() ; call init_rand(A)

        ! Initialize Krylov subspaces.
        allocate(U(1:kdim+1)) ; allocate(V(1:kdim+1)) ; allocate(X0(1))
        call init_rand(X0) ; call initialize_krylov_subspace(U, X0)
        call initialize_krylov_subspace(V) ; B = 0.0_sp

        ! Lanczos bidiagonalization.
        call lanczos_bidiagonalization(A, U, V, B, info)

        ! Check correctness.
        Id = eye(kdim)
        G = 0.0_sp
        do j = 1, kdim
            do i = 1, kdim
                G(i, j) = V(i)%dot(V(j))
            enddo
        enddo

        call check(error, maxval(abs(Id - G)) < rtol_sp)

        return
    end subroutine test_lanczos_bidiag_right_orthogonality_rsp

    subroutine collect_lanczos_bidiag_rdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("Lanczos Bidiagonalization", test_lanczos_bidiag_factorization_rdp), &
            new_unittest("Lanczos left orthogonality", test_lanczos_bidiag_left_orthogonality_rdp), &
            new_unittest("Lanczos right orthogonality", test_lanczos_bidiag_right_orthogonality_rdp) &
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
        real(dp) :: Udata(test_size, kdim+1), Vdata(test_size, kdim+1)
        type(vector_rdp), allocatable :: X0(:)

        ! Initialize linear operator.
        A = linop_rdp() ; call init_rand(A)

        ! Initialize Krylov subspaces.
        allocate(U(1:kdim+1)) ; allocate(V(1:kdim+1)) ; allocate(X0(1))
        call init_rand(X0) ; call initialize_krylov_subspace(U, X0)
        call initialize_krylov_subspace(V) ; B = 0.0_dp

        ! Lanczos bidiagonalization.
        call lanczos_bidiagonalization(A, U, V, B, info)

        ! Check correctness.
        call get_data(Udata, U)
        call get_data(Vdata, V)

        alpha = maxval(abs(matmul(A%data, Vdata(:, 1:kdim)) - matmul(Udata, B)))
        call check(error, alpha < rtol_dp)

        return
    end subroutine test_lanczos_bidiag_factorization_rdp

    subroutine test_lanczos_bidiag_left_orthogonality_rdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test Linear Operator
        type(linop_rdp), allocatable :: A
        ! Left and right Krylov bases.
        type(vector_rdp), allocatable :: U(:), V(:)
        ! Krylov subspace dimension.
        integer, parameter :: kdim = 3!test_size
        ! Bidiagonal matrix.
        real(dp) :: B(kdim+1, kdim)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(dp) :: Id(kdim, kdim), G(kdim, kdim)
        type(vector_rdp), allocatable :: X0(:)
        integer :: i, j

        ! Initialize linear operator.
        A = linop_rdp() ; call init_rand(A)

        ! Initialize Krylov subspaces.
        allocate(U(1:kdim+1)) ; allocate(V(1:kdim+1)) ; allocate(X0(1))
        call init_rand(X0) ; call initialize_krylov_subspace(U, X0)
        call initialize_krylov_subspace(V) ; B = 0.0_dp

        ! Lanczos bidiagonalization.
        call lanczos_bidiagonalization(A, U, V, B, info)

        ! Check correctness.
        Id = eye(kdim)
        G = 0.0_dp
        do j = 1, kdim
            do i = 1, kdim
                G(i, j) = U(i)%dot(U(j))
            enddo
        enddo

        call check(error, maxval(abs(Id - G)) < rtol_dp)

        return
    end subroutine test_lanczos_bidiag_left_orthogonality_rdp

    subroutine test_lanczos_bidiag_right_orthogonality_rdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test Linear Operator
        type(linop_rdp), allocatable :: A
        ! Left and right Krylov bases.
        type(vector_rdp), allocatable :: U(:), V(:)
        ! Krylov subspace dimension.
        integer, parameter :: kdim = 3!test_size
        ! Bidiagonal matrix.
        real(dp) :: B(kdim+1, kdim)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(dp) :: Id(kdim, kdim), G(kdim, kdim)
        type(vector_rdp), allocatable :: X0(:)
        integer :: i, j

        ! Initialize linear operator.
        A = linop_rdp() ; call init_rand(A)

        ! Initialize Krylov subspaces.
        allocate(U(1:kdim+1)) ; allocate(V(1:kdim+1)) ; allocate(X0(1))
        call init_rand(X0) ; call initialize_krylov_subspace(U, X0)
        call initialize_krylov_subspace(V) ; B = 0.0_dp

        ! Lanczos bidiagonalization.
        call lanczos_bidiagonalization(A, U, V, B, info)

        ! Check correctness.
        Id = eye(kdim)
        G = 0.0_dp
        do j = 1, kdim
            do i = 1, kdim
                G(i, j) = V(i)%dot(V(j))
            enddo
        enddo

        call check(error, maxval(abs(Id - G)) < rtol_dp)

        return
    end subroutine test_lanczos_bidiag_right_orthogonality_rdp

    subroutine collect_lanczos_bidiag_csp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("Lanczos Bidiagonalization", test_lanczos_bidiag_factorization_csp), &
            new_unittest("Lanczos left orthogonality", test_lanczos_bidiag_left_orthogonality_csp), &
            new_unittest("Lanczos right orthogonality", test_lanczos_bidiag_right_orthogonality_csp) &
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
        complex(sp) :: Udata(test_size, kdim+1), Vdata(test_size, kdim+1)
        type(vector_csp), allocatable :: X0(:)

        ! Initialize linear operator.
        A = linop_csp() ; call init_rand(A)

        ! Initialize Krylov subspaces.
        allocate(U(1:kdim+1)) ; allocate(V(1:kdim+1)) ; allocate(X0(1))
        call init_rand(X0) ; call initialize_krylov_subspace(U, X0)
        call initialize_krylov_subspace(V) ; B = 0.0_sp

        ! Lanczos bidiagonalization.
        call lanczos_bidiagonalization(A, U, V, B, info)

        ! Check correctness.
        call get_data(Udata, U)
        call get_data(Vdata, V)

        alpha = maxval(abs(matmul(A%data, Vdata(:, 1:kdim)) - matmul(Udata, B)))
        call check(error, alpha < rtol_sp)

        return
    end subroutine test_lanczos_bidiag_factorization_csp

    subroutine test_lanczos_bidiag_left_orthogonality_csp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test Linear Operator
        type(linop_csp), allocatable :: A
        ! Left and right Krylov bases.
        type(vector_csp), allocatable :: U(:), V(:)
        ! Krylov subspace dimension.
        integer, parameter :: kdim = 3!test_size
        ! Bidiagonal matrix.
        complex(sp) :: B(kdim+1, kdim)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(sp) :: Id(kdim, kdim), G(kdim, kdim)
        type(vector_csp), allocatable :: X0(:)
        integer :: i, j

        ! Initialize linear operator.
        A = linop_csp() ; call init_rand(A)

        ! Initialize Krylov subspaces.
        allocate(U(1:kdim+1)) ; allocate(V(1:kdim+1)) ; allocate(X0(1))
        call init_rand(X0) ; call initialize_krylov_subspace(U, X0)
        call initialize_krylov_subspace(V) ; B = 0.0_sp

        ! Lanczos bidiagonalization.
        call lanczos_bidiagonalization(A, U, V, B, info)

        ! Check correctness.
        Id = eye(kdim)
        G = 0.0_sp
        do j = 1, kdim
            do i = 1, kdim
                G(i, j) = U(i)%dot(U(j))
            enddo
        enddo

        call check(error, maxval(abs(Id - G)) < rtol_sp)

        return
    end subroutine test_lanczos_bidiag_left_orthogonality_csp

    subroutine test_lanczos_bidiag_right_orthogonality_csp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test Linear Operator
        type(linop_csp), allocatable :: A
        ! Left and right Krylov bases.
        type(vector_csp), allocatable :: U(:), V(:)
        ! Krylov subspace dimension.
        integer, parameter :: kdim = 3!test_size
        ! Bidiagonal matrix.
        complex(sp) :: B(kdim+1, kdim)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(sp) :: Id(kdim, kdim), G(kdim, kdim)
        type(vector_csp), allocatable :: X0(:)
        integer :: i, j

        ! Initialize linear operator.
        A = linop_csp() ; call init_rand(A)

        ! Initialize Krylov subspaces.
        allocate(U(1:kdim+1)) ; allocate(V(1:kdim+1)) ; allocate(X0(1))
        call init_rand(X0) ; call initialize_krylov_subspace(U, X0)
        call initialize_krylov_subspace(V) ; B = 0.0_sp

        ! Lanczos bidiagonalization.
        call lanczos_bidiagonalization(A, U, V, B, info)

        ! Check correctness.
        Id = eye(kdim)
        G = 0.0_sp
        do j = 1, kdim
            do i = 1, kdim
                G(i, j) = V(i)%dot(V(j))
            enddo
        enddo

        call check(error, maxval(abs(Id - G)) < rtol_sp)

        return
    end subroutine test_lanczos_bidiag_right_orthogonality_csp

    subroutine collect_lanczos_bidiag_cdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("Lanczos Bidiagonalization", test_lanczos_bidiag_factorization_cdp), &
            new_unittest("Lanczos left orthogonality", test_lanczos_bidiag_left_orthogonality_cdp), &
            new_unittest("Lanczos right orthogonality", test_lanczos_bidiag_right_orthogonality_cdp) &
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
        complex(dp) :: Udata(test_size, kdim+1), Vdata(test_size, kdim+1)
        type(vector_cdp), allocatable :: X0(:)

        ! Initialize linear operator.
        A = linop_cdp() ; call init_rand(A)

        ! Initialize Krylov subspaces.
        allocate(U(1:kdim+1)) ; allocate(V(1:kdim+1)) ; allocate(X0(1))
        call init_rand(X0) ; call initialize_krylov_subspace(U, X0)
        call initialize_krylov_subspace(V) ; B = 0.0_dp

        ! Lanczos bidiagonalization.
        call lanczos_bidiagonalization(A, U, V, B, info)

        ! Check correctness.
        call get_data(Udata, U)
        call get_data(Vdata, V)

        alpha = maxval(abs(matmul(A%data, Vdata(:, 1:kdim)) - matmul(Udata, B)))
        call check(error, alpha < rtol_dp)

        return
    end subroutine test_lanczos_bidiag_factorization_cdp

    subroutine test_lanczos_bidiag_left_orthogonality_cdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test Linear Operator
        type(linop_cdp), allocatable :: A
        ! Left and right Krylov bases.
        type(vector_cdp), allocatable :: U(:), V(:)
        ! Krylov subspace dimension.
        integer, parameter :: kdim = 3!test_size
        ! Bidiagonal matrix.
        complex(dp) :: B(kdim+1, kdim)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(dp) :: Id(kdim, kdim), G(kdim, kdim)
        type(vector_cdp), allocatable :: X0(:)
        integer :: i, j

        ! Initialize linear operator.
        A = linop_cdp() ; call init_rand(A)

        ! Initialize Krylov subspaces.
        allocate(U(1:kdim+1)) ; allocate(V(1:kdim+1)) ; allocate(X0(1))
        call init_rand(X0) ; call initialize_krylov_subspace(U, X0)
        call initialize_krylov_subspace(V) ; B = 0.0_dp

        ! Lanczos bidiagonalization.
        call lanczos_bidiagonalization(A, U, V, B, info)

        ! Check correctness.
        Id = eye(kdim)
        G = 0.0_dp
        do j = 1, kdim
            do i = 1, kdim
                G(i, j) = U(i)%dot(U(j))
            enddo
        enddo

        call check(error, maxval(abs(Id - G)) < rtol_dp)

        return
    end subroutine test_lanczos_bidiag_left_orthogonality_cdp

    subroutine test_lanczos_bidiag_right_orthogonality_cdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test Linear Operator
        type(linop_cdp), allocatable :: A
        ! Left and right Krylov bases.
        type(vector_cdp), allocatable :: U(:), V(:)
        ! Krylov subspace dimension.
        integer, parameter :: kdim = 3!test_size
        ! Bidiagonal matrix.
        complex(dp) :: B(kdim+1, kdim)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(dp) :: Id(kdim, kdim), G(kdim, kdim)
        type(vector_cdp), allocatable :: X0(:)
        integer :: i, j

        ! Initialize linear operator.
        A = linop_cdp() ; call init_rand(A)

        ! Initialize Krylov subspaces.
        allocate(U(1:kdim+1)) ; allocate(V(1:kdim+1)) ; allocate(X0(1))
        call init_rand(X0) ; call initialize_krylov_subspace(U, X0)
        call initialize_krylov_subspace(V) ; B = 0.0_dp

        ! Lanczos bidiagonalization.
        call lanczos_bidiagonalization(A, U, V, B, info)

        ! Check correctness.
        Id = eye(kdim)
        G = 0.0_dp
        do j = 1, kdim
            do i = 1, kdim
                G(i, j) = V(i)%dot(V(j))
            enddo
        enddo

        call check(error, maxval(abs(Id - G)) < rtol_dp)

        return
    end subroutine test_lanczos_bidiag_right_orthogonality_cdp


    !-------------------------------------------------------------------------------
    !-----     DEFINITION OF THE UNIT TESTS FOR LANCZOS TRIDIAGONALIZATION     -----
    !-------------------------------------------------------------------------------

    subroutine collect_lanczos_tridiag_rsp_testsuite(testsuite)
        ! Collection of unit tests.
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
             new_unittest("Lanczos Tridiagonalization", test_lanczos_tridiag_full_factorization_rsp), &
             new_unittest("Lanczos Tridiagonalization orthogonality", test_lanczos_tridiag_orthogonality_rsp) &
                    ]

        return
    end subroutine collect_lanczos_tridiag_rsp_testsuite

    subroutine test_lanczos_tridiag_full_factorization_rsp(error)
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

        ! Initialize tridiagonal matrix.
        T = 0.0_sp

        ! Initialize operator.
        A = spd_linop_rsp()
        call init_rand(A)

        ! Initialize Krylov subspace.
        allocate(X(kdim+1)) ; allocate(X0(1))
        call init_rand(X0) ; call initialize_krylov_subspace(X, X0)

        ! Lanczos factorization.
        call lanczos_tridiagonalization(A, X, T, info)

        ! Check correctness.
        call get_data(Xdata, X)

        ! Infinity-norm check.
        alpha = maxval(abs(matmul(A%data, Xdata(:, 1:kdim)) - matmul(Xdata, T)))
        call check(error, alpha < rtol_sp)

        return
    end subroutine test_lanczos_tridiag_full_factorization_rsp

    subroutine test_lanczos_tridiag_orthogonality_rsp(error)
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
        real(sp) :: alpha
        real(sp) :: Xdata(test_size, kdim+1)
        class(vector_rsp), allocatable :: X0(:)
        real(sp) :: Id(kdim+1, kdim+1)
        real(sp) :: G(kdim+1, kdim+1)

        ! Initialize tridiagonal matrix.
        T = 0.0_sp

        ! Initialize operator.
        A = spd_linop_rsp()
        call init_rand(A)

        ! Initialize Krylov subspace.
        allocate(X(kdim+1)) ; allocate(X0(1))
        call init_rand(X0) ; call initialize_krylov_subspace(X, X0)

        ! Lanczos factorization.
        call lanczos_tridiagonalization(A, X, T, info)

        ! Check correctness.
        call get_data(Xdata, X)
        G = matmul(transpose(Xdata), Xdata)
        Id = 0.0_sp ; Id(1:kdim, 1:kdim) = eye(kdim)
        call check(error, norm2(abs(G-Id)) < rtol_sp)

        return
    end subroutine test_lanczos_tridiag_orthogonality_rsp

    subroutine collect_lanczos_tridiag_rdp_testsuite(testsuite)
        ! Collection of unit tests.
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
             new_unittest("Lanczos Tridiagonalization", test_lanczos_tridiag_full_factorization_rdp), &
             new_unittest("Lanczos Tridiagonalization orthogonality", test_lanczos_tridiag_orthogonality_rdp) &
                    ]

        return
    end subroutine collect_lanczos_tridiag_rdp_testsuite

    subroutine test_lanczos_tridiag_full_factorization_rdp(error)
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

        ! Initialize tridiagonal matrix.
        T = 0.0_dp

        ! Initialize operator.
        A = spd_linop_rdp()
        call init_rand(A)

        ! Initialize Krylov subspace.
        allocate(X(kdim+1)) ; allocate(X0(1))
        call init_rand(X0) ; call initialize_krylov_subspace(X, X0)

        ! Lanczos factorization.
        call lanczos_tridiagonalization(A, X, T, info)

        ! Check correctness.
        call get_data(Xdata, X)

        ! Infinity-norm check.
        alpha = maxval(abs(matmul(A%data, Xdata(:, 1:kdim)) - matmul(Xdata, T)))
        call check(error, alpha < rtol_dp)

        return
    end subroutine test_lanczos_tridiag_full_factorization_rdp

    subroutine test_lanczos_tridiag_orthogonality_rdp(error)
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
        real(dp) :: alpha
        real(dp) :: Xdata(test_size, kdim+1)
        class(vector_rdp), allocatable :: X0(:)
        real(dp) :: Id(kdim+1, kdim+1)
        real(dp) :: G(kdim+1, kdim+1)

        ! Initialize tridiagonal matrix.
        T = 0.0_dp

        ! Initialize operator.
        A = spd_linop_rdp()
        call init_rand(A)

        ! Initialize Krylov subspace.
        allocate(X(kdim+1)) ; allocate(X0(1))
        call init_rand(X0) ; call initialize_krylov_subspace(X, X0)

        ! Lanczos factorization.
        call lanczos_tridiagonalization(A, X, T, info)

        ! Check correctness.
        call get_data(Xdata, X)
        G = matmul(transpose(Xdata), Xdata)
        Id = 0.0_dp ; Id(1:kdim, 1:kdim) = eye(kdim)
        call check(error, norm2(abs(G-Id)) < rtol_dp)

        return
    end subroutine test_lanczos_tridiag_orthogonality_rdp

    subroutine collect_lanczos_tridiag_csp_testsuite(testsuite)
        ! Collection of unit tests.
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
             new_unittest("Lanczos Tridiagonalization", test_lanczos_tridiag_full_factorization_csp), &
             new_unittest("Lanczos Tridiagonalization orthogonality", test_lanczos_tridiag_orthogonality_csp) &
                    ]

        return
    end subroutine collect_lanczos_tridiag_csp_testsuite

    subroutine test_lanczos_tridiag_full_factorization_csp(error)
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

        ! Initialize tridiagonal matrix.
        T = 0.0_sp

        ! Initialize operator.
        A = hermitian_linop_csp()
        call init_rand(A)

        ! Initialize Krylov subspace.
        allocate(X(kdim+1)) ; allocate(X0(1))
        call init_rand(X0) ; call initialize_krylov_subspace(X, X0)

        ! Lanczos factorization.
        call lanczos_tridiagonalization(A, X, T, info)

        ! Check correctness.
        call get_data(Xdata, X)

        ! Infinity-norm check.
        alpha = maxval(abs(matmul(A%data, Xdata(:, 1:kdim)) - matmul(Xdata, T)))
        call check(error, alpha < rtol_sp)

        return
    end subroutine test_lanczos_tridiag_full_factorization_csp

    subroutine test_lanczos_tridiag_orthogonality_csp(error)
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
        real(sp) :: alpha
        complex(sp) :: Xdata(test_size, kdim+1)
        class(vector_csp), allocatable :: X0(:)
        real(sp) :: Id(kdim+1, kdim+1)
        complex(sp) :: G(kdim+1, kdim+1)

        ! Initialize tridiagonal matrix.
        T = 0.0_sp

        ! Initialize operator.
        A = hermitian_linop_csp()
        call init_rand(A)

        ! Initialize Krylov subspace.
        allocate(X(kdim+1)) ; allocate(X0(1))
        call init_rand(X0) ; call initialize_krylov_subspace(X, X0)

        ! Lanczos factorization.
        call lanczos_tridiagonalization(A, X, T, info)

        ! Check correctness.
        call get_data(Xdata, X)
        G = matmul(transpose(conjg(Xdata)), Xdata)
        Id = 0.0_sp ; Id(1:kdim, 1:kdim) = eye(kdim)
        call check(error, norm2(abs(G-Id)) < rtol_sp)

        return
    end subroutine test_lanczos_tridiag_orthogonality_csp

    subroutine collect_lanczos_tridiag_cdp_testsuite(testsuite)
        ! Collection of unit tests.
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
             new_unittest("Lanczos Tridiagonalization", test_lanczos_tridiag_full_factorization_cdp), &
             new_unittest("Lanczos Tridiagonalization orthogonality", test_lanczos_tridiag_orthogonality_cdp) &
                    ]

        return
    end subroutine collect_lanczos_tridiag_cdp_testsuite

    subroutine test_lanczos_tridiag_full_factorization_cdp(error)
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

        ! Initialize tridiagonal matrix.
        T = 0.0_dp

        ! Initialize operator.
        A = hermitian_linop_cdp()
        call init_rand(A)

        ! Initialize Krylov subspace.
        allocate(X(kdim+1)) ; allocate(X0(1))
        call init_rand(X0) ; call initialize_krylov_subspace(X, X0)

        ! Lanczos factorization.
        call lanczos_tridiagonalization(A, X, T, info)

        ! Check correctness.
        call get_data(Xdata, X)

        ! Infinity-norm check.
        alpha = maxval(abs(matmul(A%data, Xdata(:, 1:kdim)) - matmul(Xdata, T)))
        call check(error, alpha < rtol_dp)

        return
    end subroutine test_lanczos_tridiag_full_factorization_cdp

    subroutine test_lanczos_tridiag_orthogonality_cdp(error)
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
        real(dp) :: alpha
        complex(dp) :: Xdata(test_size, kdim+1)
        class(vector_cdp), allocatable :: X0(:)
        real(dp) :: Id(kdim+1, kdim+1)
        complex(dp) :: G(kdim+1, kdim+1)

        ! Initialize tridiagonal matrix.
        T = 0.0_dp

        ! Initialize operator.
        A = hermitian_linop_cdp()
        call init_rand(A)

        ! Initialize Krylov subspace.
        allocate(X(kdim+1)) ; allocate(X0(1))
        call init_rand(X0) ; call initialize_krylov_subspace(X, X0)

        ! Lanczos factorization.
        call lanczos_tridiagonalization(A, X, T, info)

        ! Check correctness.
        call get_data(Xdata, X)
        G = matmul(transpose(conjg(Xdata)), Xdata)
        Id = 0.0_dp ; Id(1:kdim, 1:kdim) = eye(kdim)
        call check(error, norm2(abs(G-Id)) < rtol_dp)

        return
    end subroutine test_lanczos_tridiag_orthogonality_cdp


end module TestKrylov
