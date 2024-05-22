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
                        new_unittest("Dense expm.", test_dense_expm_rsp), &
                        new_unittest("Krylov expm.", test_kexptA_rsp) &
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

    subroutine test_kexptA_rsp(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        !> Test matrix.
        class(linop_rsp), allocatable :: A
        !> Basis vectors.
        class(vector_rsp), allocatable :: Q, Xref, Xkryl
        !> Krylov subspace dimension.
        integer, parameter :: kdim = test_size
    
        ! ----- Internal variables -----
        real(sp) :: E(kdim, kdim)
        real(sp), parameter :: tau = 0.1_sp
        logical, parameter :: verb = .true.
        integer, parameter :: nkmax = 64
        real(sp) :: err
        integer :: info
        
        ! Initialize data.
        A = linop_rsp() ; call init_rand(A)
        allocate(Q) ; call init_rand(Q)
        allocate(Xref) ; call Xref%zero()
        allocate(XKryl) ; call Xkryl%zero()

        ! Dense computation.
        call expm(E, tau*A%data)
        Xref%data = matmul(E, Q%data)

        ! Krylov exponential.
        call kexpm(Xkryl, A, Q, tau, rtol_sp, info, verbosity=verb, kdim=nkmax)

        call save_npy("test_krylov_expm_operator.npy", A%data)
        call save_npy("test_krylov_expm_rhs.npy", Q%data)
        call save_npy("test_krylov_expm_ref.npy", Xref%data)
        call save_npy("test_krylov_expm_kexpm.npy", Xkryl%data)
 
        ! Check result.
        call Xkryl%sub(Xref) ; err = Xkryl%norm()
        if (verb) write(output_unit, *) "     True error: ||error||_2 =", err

       call check(error, err < rtol_sp)

        return
    end subroutine test_kexptA_rsp

    subroutine collect_expm_rdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                        new_unittest("Dense expm.", test_dense_expm_rdp), &
                        new_unittest("Krylov expm.", test_kexptA_rdp) &
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

    subroutine test_kexptA_rdp(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        !> Test matrix.
        class(linop_rdp), allocatable :: A
        !> Basis vectors.
        class(vector_rdp), allocatable :: Q, Xref, Xkryl
        !> Krylov subspace dimension.
        integer, parameter :: kdim = test_size
    
        ! ----- Internal variables -----
        real(dp) :: E(kdim, kdim)
        real(dp), parameter :: tau = 0.1_dp
        logical, parameter :: verb = .true.
        integer, parameter :: nkmax = 64
        real(dp) :: err
        integer :: info
        
        ! Initialize data.
        A = linop_rdp() ; call init_rand(A)
        allocate(Q) ; call init_rand(Q)
        allocate(Xref) ; call Xref%zero()
        allocate(XKryl) ; call Xkryl%zero()

        ! Dense computation.
        call expm(E, tau*A%data)
        Xref%data = matmul(E, Q%data)

        ! Krylov exponential.
        call kexpm(Xkryl, A, Q, tau, rtol_dp, info, verbosity=verb, kdim=nkmax)

        call save_npy("test_krylov_expm_operator.npy", A%data)
        call save_npy("test_krylov_expm_rhs.npy", Q%data)
        call save_npy("test_krylov_expm_ref.npy", Xref%data)
        call save_npy("test_krylov_expm_kexpm.npy", Xkryl%data)
 
        ! Check result.
        call Xkryl%sub(Xref) ; err = Xkryl%norm()
        if (verb) write(output_unit, *) "     True error: ||error||_2 =", err

       call check(error, err < rtol_dp)

        return
    end subroutine test_kexptA_rdp

    subroutine collect_expm_csp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                        new_unittest("Dense expm.", test_dense_expm_csp), &
                        new_unittest("Krylov expm.", test_kexptA_csp) &
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

    subroutine test_kexptA_csp(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        !> Test matrix.
        class(linop_csp), allocatable :: A
        !> Basis vectors.
        class(vector_csp), allocatable :: Q, Xref, Xkryl
        !> Krylov subspace dimension.
        integer, parameter :: kdim = test_size
    
        ! ----- Internal variables -----
        complex(sp) :: E(kdim, kdim)
        real(sp), parameter :: tau = 0.1_sp
        logical, parameter :: verb = .true.
        integer, parameter :: nkmax = 64
        real(sp) :: err
        integer :: info
        
        ! Initialize data.
        A = linop_csp() ; call init_rand(A)
        allocate(Q) ; call init_rand(Q)
        allocate(Xref) ; call Xref%zero()
        allocate(XKryl) ; call Xkryl%zero()

        ! Dense computation.
        call expm(E, tau*A%data)
        Xref%data = matmul(E, Q%data)

        ! Krylov exponential.
        call kexpm(Xkryl, A, Q, tau, rtol_sp, info, verbosity=verb, kdim=nkmax)

        call save_npy("test_krylov_expm_operator.npy", A%data)
        call save_npy("test_krylov_expm_rhs.npy", Q%data)
        call save_npy("test_krylov_expm_ref.npy", Xref%data)
        call save_npy("test_krylov_expm_kexpm.npy", Xkryl%data)
 
        ! Check result.
        call Xkryl%sub(Xref) ; err = Xkryl%norm()
        if (verb) write(output_unit, *) "     True error: ||error||_2 =", err

       call check(error, err < rtol_sp)

        return
    end subroutine test_kexptA_csp

    subroutine collect_expm_cdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                        new_unittest("Dense expm.", test_dense_expm_cdp), &
                        new_unittest("Krylov expm.", test_kexptA_cdp) &
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

    subroutine test_kexptA_cdp(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        !> Test matrix.
        class(linop_cdp), allocatable :: A
        !> Basis vectors.
        class(vector_cdp), allocatable :: Q, Xref, Xkryl
        !> Krylov subspace dimension.
        integer, parameter :: kdim = test_size
    
        ! ----- Internal variables -----
        complex(dp) :: E(kdim, kdim)
        real(dp), parameter :: tau = 0.1_dp
        logical, parameter :: verb = .true.
        integer, parameter :: nkmax = 64
        real(dp) :: err
        integer :: info
        
        ! Initialize data.
        A = linop_cdp() ; call init_rand(A)
        allocate(Q) ; call init_rand(Q)
        allocate(Xref) ; call Xref%zero()
        allocate(XKryl) ; call Xkryl%zero()

        ! Dense computation.
        call expm(E, tau*A%data)
        Xref%data = matmul(E, Q%data)

        ! Krylov exponential.
        call kexpm(Xkryl, A, Q, tau, rtol_dp, info, verbosity=verb, kdim=nkmax)

        call save_npy("test_krylov_expm_operator.npy", A%data)
        call save_npy("test_krylov_expm_rhs.npy", Q%data)
        call save_npy("test_krylov_expm_ref.npy", Xref%data)
        call save_npy("test_krylov_expm_kexpm.npy", Xkryl%data)
 
        ! Check result.
        call Xkryl%sub(Xref) ; err = Xkryl%norm()
        if (verb) write(output_unit, *) "     True error: ||error||_2 =", err

       call check(error, err < rtol_dp)

        return
    end subroutine test_kexptA_cdp


end module

