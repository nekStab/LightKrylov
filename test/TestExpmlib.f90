module TestExpmlib
    ! Fortran Standard Library.
    use iso_fortran_env
    use stdlib_math, only: is_close, all_close
    use stdlib_linalg, only: eye
    use stdlib_io_npy, only: save_npy

    ! LightKrylov
    use LightKrylov

    ! Testdrive
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use TestVectors
    use TestLinops
    use TestUtils
    use TestKrylov

    public :: collect_expm_rsp_testsuite
    public :: collect_expm_rdp_testsuite
    public :: collect_expm_csp_testsuite
    public :: collect_expm_cdp_testsuite

contains

    !--------------------------------------------------------------
    !-----     UNIT TESTS FOR DENSE MATRIX EXPONENTIATION     -----
    !--------------------------------------------------------------

    subroutine collect_expm_rsp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                        new_unittest("Dense expm.", test_dense_expm_rsp) &
                    ]

        return
    end subroutine collect_expm_rsp_testsuite

    subroutine test_dense_expm_rsp(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        !> Problem dimension.
        integer, parameter :: n = 5, m = 6
        !> Test matrix.
        real(sp) :: A(n, n), E(n, n), Eref(n, n)
        integer :: i, j

        ! Initialize matrix.
        A = 0.0_sp
        do i = 1, n-1
            A(i, i+1) = m*1.0_sp
        enddo

        ! Reference with analytical exponential.
        Eref = eye(n)
        do i = 1, n-1
            do j = 1, n-i
                Eref(i, i+j) = Eref(i, i+j-1)*m/j
            enddo
        enddo

        ! Compute matrix exponential.
        call expm(E, A)

        ! Check result.
        call check(error, maxval(abs(E-Eref)) < rtol_sp)

        return
    end subroutine test_dense_expm_rsp

    subroutine collect_expm_rdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                        new_unittest("Dense expm.", test_dense_expm_rdp) &
                    ]

        return
    end subroutine collect_expm_rdp_testsuite

    subroutine test_dense_expm_rdp(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        !> Problem dimension.
        integer, parameter :: n = 5, m = 6
        !> Test matrix.
        real(dp) :: A(n, n), E(n, n), Eref(n, n)
        integer :: i, j

        ! Initialize matrix.
        A = 0.0_dp
        do i = 1, n-1
            A(i, i+1) = m*1.0_dp
        enddo

        ! Reference with analytical exponential.
        Eref = eye(n)
        do i = 1, n-1
            do j = 1, n-i
                Eref(i, i+j) = Eref(i, i+j-1)*m/j
            enddo
        enddo

        ! Compute matrix exponential.
        call expm(E, A)

        ! Check result.
        call check(error, maxval(abs(E-Eref)) < rtol_dp)

        return
    end subroutine test_dense_expm_rdp

    subroutine collect_expm_csp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                        new_unittest("Dense expm.", test_dense_expm_csp) &
                    ]

        return
    end subroutine collect_expm_csp_testsuite

    subroutine test_dense_expm_csp(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        !> Problem dimension.
        integer, parameter :: n = 5, m = 6
        !> Test matrix.
        complex(sp) :: A(n, n), E(n, n), Eref(n, n)
        integer :: i, j

        ! Initialize matrix.
        A = 0.0_sp
        do i = 1, n-1
            A(i, i+1) = m*1.0_sp
        enddo

        ! Reference with analytical exponential.
        Eref = eye(n)
        do i = 1, n-1
            do j = 1, n-i
                Eref(i, i+j) = Eref(i, i+j-1)*m/j
            enddo
        enddo

        ! Compute matrix exponential.
        call expm(E, A)

        ! Check result.
        call check(error, maxval(abs(E-Eref)) < rtol_sp)

        return
    end subroutine test_dense_expm_csp

    subroutine collect_expm_cdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                        new_unittest("Dense expm.", test_dense_expm_cdp) &
                    ]

        return
    end subroutine collect_expm_cdp_testsuite

    subroutine test_dense_expm_cdp(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        !> Problem dimension.
        integer, parameter :: n = 5, m = 6
        !> Test matrix.
        complex(dp) :: A(n, n), E(n, n), Eref(n, n)
        integer :: i, j

        ! Initialize matrix.
        A = 0.0_dp
        do i = 1, n-1
            A(i, i+1) = m*1.0_dp
        enddo

        ! Reference with analytical exponential.
        Eref = eye(n)
        do i = 1, n-1
            do j = 1, n-i
                Eref(i, i+j) = Eref(i, i+j-1)*m/j
            enddo
        enddo

        ! Compute matrix exponential.
        call expm(E, A)

        ! Check result.
        call check(error, maxval(abs(E-Eref)) < rtol_dp)

        return
    end subroutine test_dense_expm_cdp


end module

