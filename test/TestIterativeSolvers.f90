module TestIterativeSolvers
    ! Fortran Standard library.
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
    use TestKrylov

    implicit none
    private

    public :: collect_eig_rsp_testsuite
    public :: collect_eig_rdp_testsuite
    public :: collect_eig_csp_testsuite
    public :: collect_eig_cdp_testsuite

contains

    !---------------------------------------------------------
    !-----     DEFINITION OF THE UNIT TESTS FOR EIGS     -----
    !---------------------------------------------------------

    subroutine collect_eig_rsp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("Eigs computation", test_evp_rsp) &
                    ]
        return
    end subroutine collect_eig_rsp_testsuite

    subroutine test_evp_rsp(error)
        !> Error type.
        type(error_type), allocatable, intent(out) :: error
        !> Test linear operator.
        type(linop_rsp), allocatable :: A
        !> Eigenvectors.
        type(vector_rsp), allocatable :: X(:)
        !> Eigenvalues.
        complex(sp), allocatable :: eigvals(:)
        !> Residuals.
        real(sp), allocatable :: residuals(:)
        !> Information flag.
        integer :: info
        !> Miscellaneous.
        integer :: i, k, n
        complex(sp) :: true_eigvals(test_size)
        real(sp) :: pi = 4.0_sp * atan(1.0_sp)
        real(sp) :: a_, b_

        ! Allocate eigenvectors.
        allocate(X(test_size)) ; call initialize_krylov_subspace(X)
        
        ! Initialize linear operator with random tridiagonal Toeplitz matrix.
        A = linop_rsp() ; A%data = 0.0_sp ; n = size(A%data, 1)

        call random_number(a_) ; call random_number(b_) ; b_ = abs(b_)
        do i = 1, n
            ! Diagonal entry.
            A%data(i, i) = a_
            ! Upper diagonal entry.
            if (i < n) then
                A%data(i, i+1) = b_
                A%data(i+1, i) = -b_
            endif
        enddo

        ! Compute spectral decomposition.
        call eigs(A, X, eigvals, residuals, info)

        ! Analytical eigenvalues.
        true_eigvals = cmplx(0.0_sp, 0.0_sp, kind=sp) ; k = 1
        do i = 1, test_size, 2
            true_eigvals(i) = a_*cmplx(1.0_sp, 0.0_sp, kind=sp) + (2*b_*cos(k*pi/(test_size+1)))*cmplx(0.0_sp, 1.0_sp, kind=sp)
            true_eigvals(i+1) = a_*cmplx(1.0_sp, 0.0_sp, kind=sp) - (2*b_*cos(k*pi/(test_size+1)))*cmplx(0.0_sp, 1.0_sp, kind=sp)
            k = k+1
        enddo

        call check(error, norm2(abs(eigvals - true_eigvals)) < rtol_sp)

        return
    end subroutine test_evp_rsp

    subroutine collect_eig_rdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("Eigs computation", test_evp_rdp) &
                    ]
        return
    end subroutine collect_eig_rdp_testsuite

    subroutine test_evp_rdp(error)
        !> Error type.
        type(error_type), allocatable, intent(out) :: error
        !> Test linear operator.
        type(linop_rdp), allocatable :: A
        !> Eigenvectors.
        type(vector_rdp), allocatable :: X(:)
        !> Eigenvalues.
        complex(dp), allocatable :: eigvals(:)
        !> Residuals.
        real(dp), allocatable :: residuals(:)
        !> Information flag.
        integer :: info
        !> Miscellaneous.
        integer :: i, k, n
        complex(dp) :: true_eigvals(test_size)
        real(dp) :: pi = 4.0_dp * atan(1.0_dp)
        real(dp) :: a_, b_

        ! Allocate eigenvectors.
        allocate(X(test_size)) ; call initialize_krylov_subspace(X)
        
        ! Initialize linear operator with random tridiagonal Toeplitz matrix.
        A = linop_rdp() ; A%data = 0.0_dp ; n = size(A%data, 1)

        call random_number(a_) ; call random_number(b_) ; b_ = abs(b_)
        do i = 1, n
            ! Diagonal entry.
            A%data(i, i) = a_
            ! Upper diagonal entry.
            if (i < n) then
                A%data(i, i+1) = b_
                A%data(i+1, i) = -b_
            endif
        enddo

        ! Compute spectral decomposition.
        call eigs(A, X, eigvals, residuals, info)

        ! Analytical eigenvalues.
        true_eigvals = cmplx(0.0_dp, 0.0_dp, kind=dp) ; k = 1
        do i = 1, test_size, 2
            true_eigvals(i) = a_*cmplx(1.0_dp, 0.0_dp, kind=dp) + (2*b_*cos(k*pi/(test_size+1)))*cmplx(0.0_dp, 1.0_dp, kind=dp)
            true_eigvals(i+1) = a_*cmplx(1.0_dp, 0.0_dp, kind=dp) - (2*b_*cos(k*pi/(test_size+1)))*cmplx(0.0_dp, 1.0_dp, kind=dp)
            k = k+1
        enddo

        call check(error, norm2(abs(eigvals - true_eigvals)) < rtol_dp)

        return
    end subroutine test_evp_rdp

    subroutine collect_eig_csp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("Eigs computation", test_evp_csp) &
                    ]
        return
    end subroutine collect_eig_csp_testsuite

    subroutine test_evp_csp(error)
        !> Error type.
        type(error_type), allocatable, intent(out) :: error
        !> Test linear operator.
        type(linop_csp), allocatable :: A
        !> Eigenvectors.
        type(vector_csp), allocatable :: X(:)
        !> Eigenvalues.
        complex(sp), allocatable :: eigvals(:)
        !> Residuals.
        real(sp), allocatable :: residuals(:)
        !> Information flag.
        integer :: info
        !> Miscellaneous.
        integer :: i, k, n
        complex(sp) :: true_eigvals(test_size)
        real(sp) :: pi = 4.0_sp * atan(1.0_sp)

        return
    end subroutine test_evp_csp

    subroutine collect_eig_cdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("Eigs computation", test_evp_cdp) &
                    ]
        return
    end subroutine collect_eig_cdp_testsuite

    subroutine test_evp_cdp(error)
        !> Error type.
        type(error_type), allocatable, intent(out) :: error
        !> Test linear operator.
        type(linop_cdp), allocatable :: A
        !> Eigenvectors.
        type(vector_cdp), allocatable :: X(:)
        !> Eigenvalues.
        complex(dp), allocatable :: eigvals(:)
        !> Residuals.
        real(dp), allocatable :: residuals(:)
        !> Information flag.
        integer :: info
        !> Miscellaneous.
        integer :: i, k, n
        complex(dp) :: true_eigvals(test_size)
        real(dp) :: pi = 4.0_dp * atan(1.0_dp)

        return
    end subroutine test_evp_cdp


end module TestIterativeSolvers

