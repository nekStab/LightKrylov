module TestExpmlib
    ! Fortran Standard Library.
    use iso_fortran_env
    use stdlib_math, only: is_close, all_close
    use stdlib_linalg, only: eye, diag
    use stdlib_io_npy, only: save_npy

    ! LightKrylov
    use LightKrylov
    use LightKrylov_Constants
    use LightKrylov_Logger
    use LightKrylov_Utils, only : eig, sqrtm

    ! Testdrive
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use TestVectors
    use TestLinops
    use TestUtils
    use TestKrylov

    implicit none

    private

    character(len=128), parameter, private :: this_module = 'LightKrylov_TestExpmLib'

    public :: collect_expm_rsp_testsuite
    public :: collect_expm_rdp_testsuite
    public :: collect_expm_csp_testsuite
    public :: collect_expm_cdp_testsuite

    public :: collect_sqrtm_rsp_testsuite
    public :: collect_sqrtm_rdp_testsuite
    public :: collect_sqrtm_csp_testsuite
    public :: collect_sqrtm_cdp_testsuite

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
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Problem dimension.
        integer, parameter :: n = 5, m = 6
        ! Test matrix.
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
        call check_test(error, 'test_dense_expm_rsp', 'maxval(abs(E-Eref)) < rtol')

        return
    end subroutine test_dense_expm_rsp

    subroutine test_kexptA_rsp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test matrix.
        class(linop_rsp), allocatable :: A
        ! Basis vectors.
        class(vector_rsp), allocatable :: Q, Xref, Xkryl
        ! Krylov subspace dimension.
        integer, parameter :: kdim = test_size
    
        ! ----- Internal variables -----
        real(sp), allocatable :: E(:, :)
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
        allocate(E(kdim, kdim)) ; call expm(E, tau*A%data)
        Xref%data = matmul(E, Q%data)

        ! Krylov exponential.
        call kexpm(Xkryl, A, Q, tau, rtol_sp, info, verbosity=verb, kdim=nkmax)
        call check_info(info, 'kexpm', module=this_module, procedure='test_kexptA_rsp')

        ! Check result.
        call Xkryl%sub(Xref) ; err = Xkryl%norm()
        if (verb) write(output_unit, *) "     True error: ||error||_2 =", err

        call check(error, err < rtol_sp)
        call check_test(error, 'test_kexptA_rsp', 'norm2(x - x_true) < rol')

        return
    end subroutine test_kexptA_rsp

    subroutine test_block_kexptA_rsp(error)
        !! This function tests the Krylov based approximation of the action of the exponential
        !! propagator against the dense computation for a random operator, a random RHS and a 
        !! typical value of tau.

        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        class(linop_rsp), allocatable :: A
        ! Basis vectors.
        class(vector_rsp), allocatable :: B(:)
        class(vector_rsp), allocatable :: Cref(:)
        class(vector_rsp), allocatable :: C(:), Cblk(:)
        ! Krylov subspace dimension.
        integer, parameter :: kdim = test_size
        ! Test matrix.
        real(sp), allocatable :: Adata(:, :), Edata(:, :)
        ! Information flag.
        integer :: info
        ! Test parameters
        integer       , parameter :: nkmax = 15
        integer       , parameter :: p = 3
        real(sp), parameter :: tau = 0.1_sp
        real(sp), parameter :: tol = rtol_sp
        logical       , parameter :: verb = .true.
        ! Misc.
        integer  :: i, j
        real(sp), allocatable :: Xdata(:, :), Qdata(:, :)
        real(sp) :: alpha
        real(sp) :: err(p, p)

        allocate(Adata(kdim, kdim)) ; Adata = 0.0_sp
        allocate(Edata(kdim, kdim)) ; Edata = 0.0_sp
        allocate(Xdata(test_size, p)) ; Xdata = 0.0_sp

        allocate(Cref(p)) ; call initialize_krylov_subspace(Cref)
        allocate(C(p))    ; call initialize_krylov_subspace(C)
        allocate(Cblk(p)) ; call initialize_krylov_subspace(Cblk)

        ! --> Initialize operator.
        A = linop_rsp() ; call init_rand(A) ; call get_data(Adata, A)
       
        ! --> Initialize rhs.
        allocate(B(1:p)) ; call init_rand(B)
        allocate(Qdata(test_size, p)) ; call get_data(Qdata, B)

        ! Comparison is dense computation (10th order Pade approximation)
        call expm(Edata, tau*Adata) ; Xdata = matmul(Edata,Qdata) ; call put_data(Cref, Xdata)

        ! Compute Krylov matrix exponential using sequential arnoldi method for each input column
        if (verb) write(*,*) 'SEQUENTIAL ARNOLDI'
        do i = 1,p
            if (verb) write(*,*) '    column',i
            call kexpm(C(i), A, B(i), tau, tol, info, verbosity=verb, kdim=nkmax)
            call check_info(info, 'kexpm', module=this_module, procedure='test_block_kexptA_rsp, 1')
        end do
        
        ! Compute Krylov matrix exponential using block-arnoldi method
        if (verb) write(*,*) 'BLOCK-ARNOLDI'
        call kexpm(Cblk, A, B, tau, tol, info, verbosity=verb, kdim=nkmax)
        call check_info(info, 'kexpm', module=this_module, procedure='test_block_kexptA_rsp, 2')
    
        do i = 1, p
            write(output_unit, *) C(i)%norm(), Cblk(i)%norm()
            call C(i)%sub(Cref(i)) ; call Cblk(i)%sub(Cref(i))
        end do

        ! Compute 2-norm of the error
        if (verb) then
            do i = 1, size(C)
                do j = 1, size(C)
                    err(i, j) = C(i)%dot(C(j))
                enddo
            enddo
            alpha = norm2(abs(err))
            write(*,*) '--------------------------------------------------------------------'
            write(*, *) '    true error (seq.):   ||error||_2 = ', alpha
        endif
 
        do i = 1, size(Cblk)
            do j = 1, size(Cblk)
                err(i, j) = Cblk(i)%dot(Cblk(j))
            enddo
        enddo
        
        alpha = norm2(abs(err))
        if (verb) write(*, *) '    true error (block):  ||error||_2 = ', alpha

        call check(error, alpha < rtol_sp)
        call check_test(error, 'test_block_kexptA_rsp', 'norm2(x - x_true) < rol')

        return
    end subroutine test_block_kexptA_rsp

    subroutine collect_expm_rdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                        new_unittest("Dense expm.", test_dense_expm_rdp), &
                        new_unittest("Krylov expm.", test_kexptA_rdp) &
                    ]

        return
    end subroutine collect_expm_rdp_testsuite

    subroutine test_dense_expm_rdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Problem dimension.
        integer, parameter :: n = 5, m = 6
        ! Test matrix.
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
        call check_test(error, 'test_dense_expm_rdp', 'maxval(abs(E-Eref)) < rtol')

        return
    end subroutine test_dense_expm_rdp

    subroutine test_kexptA_rdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test matrix.
        class(linop_rdp), allocatable :: A
        ! Basis vectors.
        class(vector_rdp), allocatable :: Q, Xref, Xkryl
        ! Krylov subspace dimension.
        integer, parameter :: kdim = test_size
    
        ! ----- Internal variables -----
        real(dp), allocatable :: E(:, :)
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
        allocate(E(kdim, kdim)) ; call expm(E, tau*A%data)
        Xref%data = matmul(E, Q%data)

        ! Krylov exponential.
        call kexpm(Xkryl, A, Q, tau, rtol_dp, info, verbosity=verb, kdim=nkmax)
        call check_info(info, 'kexpm', module=this_module, procedure='test_kexptA_rdp')

        ! Check result.
        call Xkryl%sub(Xref) ; err = Xkryl%norm()
        if (verb) write(output_unit, *) "     True error: ||error||_2 =", err

        call check(error, err < rtol_dp)
        call check_test(error, 'test_kexptA_rdp', 'norm2(x - x_true) < rol')

        return
    end subroutine test_kexptA_rdp

    subroutine test_block_kexptA_rdp(error)
        !! This function tests the Krylov based approximation of the action of the exponential
        !! propagator against the dense computation for a random operator, a random RHS and a 
        !! typical value of tau.

        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        class(linop_rdp), allocatable :: A
        ! Basis vectors.
        class(vector_rdp), allocatable :: B(:)
        class(vector_rdp), allocatable :: Cref(:)
        class(vector_rdp), allocatable :: C(:), Cblk(:)
        ! Krylov subspace dimension.
        integer, parameter :: kdim = test_size
        ! Test matrix.
        real(dp), allocatable :: Adata(:, :), Edata(:, :)
        ! Information flag.
        integer :: info
        ! Test parameters
        integer       , parameter :: nkmax = 15
        integer       , parameter :: p = 3
        real(dp), parameter :: tau = 0.1_dp
        real(dp), parameter :: tol = rtol_dp
        logical       , parameter :: verb = .true.
        ! Misc.
        integer  :: i, j
        real(dp), allocatable :: Xdata(:, :), Qdata(:, :)
        real(dp) :: alpha
        real(dp) :: err(p, p)

        allocate(Adata(kdim, kdim)) ; Adata = 0.0_dp
        allocate(Edata(kdim, kdim)) ; Edata = 0.0_dp
        allocate(Xdata(test_size, p)) ; Xdata = 0.0_dp

        allocate(Cref(p)) ; call initialize_krylov_subspace(Cref)
        allocate(C(p))    ; call initialize_krylov_subspace(C)
        allocate(Cblk(p)) ; call initialize_krylov_subspace(Cblk)

        ! --> Initialize operator.
        A = linop_rdp() ; call init_rand(A) ; call get_data(Adata, A)
       
        ! --> Initialize rhs.
        allocate(B(1:p)) ; call init_rand(B)
        allocate(Qdata(test_size, p)) ; call get_data(Qdata, B)

        ! Comparison is dense computation (10th order Pade approximation)
        call expm(Edata, tau*Adata) ; Xdata = matmul(Edata,Qdata) ; call put_data(Cref, Xdata)

        ! Compute Krylov matrix exponential using sequential arnoldi method for each input column
        if (verb) write(*,*) 'SEQUENTIAL ARNOLDI'
        do i = 1,p
            if (verb) write(*,*) '    column',i
            call kexpm(C(i), A, B(i), tau, tol, info, verbosity=verb, kdim=nkmax)
            call check_info(info, 'kexpm', module=this_module, procedure='test_block_kexptA_rdp, 1')
        end do
        
        ! Compute Krylov matrix exponential using block-arnoldi method
        if (verb) write(*,*) 'BLOCK-ARNOLDI'
        call kexpm(Cblk, A, B, tau, tol, info, verbosity=verb, kdim=nkmax)
        call check_info(info, 'kexpm', module=this_module, procedure='test_block_kexptA_rdp, 2')
    
        do i = 1, p
            write(output_unit, *) C(i)%norm(), Cblk(i)%norm()
            call C(i)%sub(Cref(i)) ; call Cblk(i)%sub(Cref(i))
        end do

        ! Compute 2-norm of the error
        if (verb) then
            do i = 1, size(C)
                do j = 1, size(C)
                    err(i, j) = C(i)%dot(C(j))
                enddo
            enddo
            alpha = norm2(abs(err))
            write(*,*) '--------------------------------------------------------------------'
            write(*, *) '    true error (seq.):   ||error||_2 = ', alpha
        endif
 
        do i = 1, size(Cblk)
            do j = 1, size(Cblk)
                err(i, j) = Cblk(i)%dot(Cblk(j))
            enddo
        enddo
        
        alpha = norm2(abs(err))
        if (verb) write(*, *) '    true error (block):  ||error||_2 = ', alpha

        call check(error, alpha < rtol_dp)
        call check_test(error, 'test_block_kexptA_rdp', 'norm2(x - x_true) < rol')

        return
    end subroutine test_block_kexptA_rdp

    subroutine collect_expm_csp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                        new_unittest("Dense expm.", test_dense_expm_csp), &
                        new_unittest("Krylov expm.", test_kexptA_csp) &
                    ]

        return
    end subroutine collect_expm_csp_testsuite

    subroutine test_dense_expm_csp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Problem dimension.
        integer, parameter :: n = 5, m = 6
        ! Test matrix.
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
        call check_test(error, 'test_dense_expm_csp', 'maxval(abs(E-Eref)) < rtol')

        return
    end subroutine test_dense_expm_csp

    subroutine test_kexptA_csp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test matrix.
        class(linop_csp), allocatable :: A
        ! Basis vectors.
        class(vector_csp), allocatable :: Q, Xref, Xkryl
        ! Krylov subspace dimension.
        integer, parameter :: kdim = test_size
    
        ! ----- Internal variables -----
        complex(sp), allocatable :: E(:, :)
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
        allocate(E(kdim, kdim)) ; call expm(E, tau*A%data)
        Xref%data = matmul(E, Q%data)

        ! Krylov exponential.
        call kexpm(Xkryl, A, Q, tau, rtol_sp, info, verbosity=verb, kdim=nkmax)
        call check_info(info, 'kexpm', module=this_module, procedure='test_kexptA_csp')

        ! Check result.
        call Xkryl%sub(Xref) ; err = Xkryl%norm()
        if (verb) write(output_unit, *) "     True error: ||error||_2 =", err

        call check(error, err < rtol_sp)
        call check_test(error, 'test_kexptA_csp', 'norm2(x - x_true) < rol')

        return
    end subroutine test_kexptA_csp

    subroutine test_block_kexptA_csp(error)
        !! This function tests the Krylov based approximation of the action of the exponential
        !! propagator against the dense computation for a random operator, a random RHS and a 
        !! typical value of tau.

        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        class(linop_csp), allocatable :: A
        ! Basis vectors.
        class(vector_csp), allocatable :: B(:)
        class(vector_csp), allocatable :: Cref(:)
        class(vector_csp), allocatable :: C(:), Cblk(:)
        ! Krylov subspace dimension.
        integer, parameter :: kdim = test_size
        ! Test matrix.
        complex(sp), allocatable :: Adata(:, :), Edata(:, :)
        ! Information flag.
        integer :: info
        ! Test parameters
        integer       , parameter :: nkmax = 15
        integer       , parameter :: p = 3
        real(sp), parameter :: tau = 0.1_sp
        real(sp), parameter :: tol = rtol_sp
        logical       , parameter :: verb = .true.
        ! Misc.
        integer  :: i, j
        complex(sp), allocatable :: Xdata(:, :), Qdata(:, :)
        real(sp) :: alpha
        complex(sp) :: err(p, p)

        allocate(Adata(kdim, kdim)) ; Adata = 0.0_sp
        allocate(Edata(kdim, kdim)) ; Edata = 0.0_sp
        allocate(Xdata(test_size, p)) ; Xdata = 0.0_sp

        allocate(Cref(p)) ; call initialize_krylov_subspace(Cref)
        allocate(C(p))    ; call initialize_krylov_subspace(C)
        allocate(Cblk(p)) ; call initialize_krylov_subspace(Cblk)

        ! --> Initialize operator.
        A = linop_csp() ; call init_rand(A) ; call get_data(Adata, A)
       
        ! --> Initialize rhs.
        allocate(B(1:p)) ; call init_rand(B)
        allocate(Qdata(test_size, p)) ; call get_data(Qdata, B)

        ! Comparison is dense computation (10th order Pade approximation)
        call expm(Edata, tau*Adata) ; Xdata = matmul(Edata,Qdata) ; call put_data(Cref, Xdata)

        ! Compute Krylov matrix exponential using sequential arnoldi method for each input column
        if (verb) write(*,*) 'SEQUENTIAL ARNOLDI'
        do i = 1,p
            if (verb) write(*,*) '    column',i
            call kexpm(C(i), A, B(i), tau, tol, info, verbosity=verb, kdim=nkmax)
            call check_info(info, 'kexpm', module=this_module, procedure='test_block_kexptA_csp, 1')
        end do
        
        ! Compute Krylov matrix exponential using block-arnoldi method
        if (verb) write(*,*) 'BLOCK-ARNOLDI'
        call kexpm(Cblk, A, B, tau, tol, info, verbosity=verb, kdim=nkmax)
        call check_info(info, 'kexpm', module=this_module, procedure='test_block_kexptA_csp, 2')
    
        do i = 1, p
            write(output_unit, *) C(i)%norm(), Cblk(i)%norm()
            call C(i)%sub(Cref(i)) ; call Cblk(i)%sub(Cref(i))
        end do

        ! Compute 2-norm of the error
        if (verb) then
            do i = 1, size(C)
                do j = 1, size(C)
                    err(i, j) = C(i)%dot(C(j))
                enddo
            enddo
            alpha = norm2(abs(err))
            write(*,*) '--------------------------------------------------------------------'
            write(*, *) '    true error (seq.):   ||error||_2 = ', alpha
        endif
 
        do i = 1, size(Cblk)
            do j = 1, size(Cblk)
                err(i, j) = Cblk(i)%dot(Cblk(j))
            enddo
        enddo
        
        alpha = norm2(abs(err))
        if (verb) write(*, *) '    true error (block):  ||error||_2 = ', alpha

        call check(error, alpha < rtol_sp)
        call check_test(error, 'test_block_kexptA_csp', 'norm2(x - x_true) < rol')

        return
    end subroutine test_block_kexptA_csp

    subroutine collect_expm_cdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                        new_unittest("Dense expm.", test_dense_expm_cdp), &
                        new_unittest("Krylov expm.", test_kexptA_cdp) &
                    ]

        return
    end subroutine collect_expm_cdp_testsuite

    subroutine test_dense_expm_cdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Problem dimension.
        integer, parameter :: n = 5, m = 6
        ! Test matrix.
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
        call check_test(error, 'test_dense_expm_cdp', 'maxval(abs(E-Eref)) < rtol')

        return
    end subroutine test_dense_expm_cdp

    subroutine test_kexptA_cdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test matrix.
        class(linop_cdp), allocatable :: A
        ! Basis vectors.
        class(vector_cdp), allocatable :: Q, Xref, Xkryl
        ! Krylov subspace dimension.
        integer, parameter :: kdim = test_size
    
        ! ----- Internal variables -----
        complex(dp), allocatable :: E(:, :)
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
        allocate(E(kdim, kdim)) ; call expm(E, tau*A%data)
        Xref%data = matmul(E, Q%data)

        ! Krylov exponential.
        call kexpm(Xkryl, A, Q, tau, rtol_dp, info, verbosity=verb, kdim=nkmax)
        call check_info(info, 'kexpm', module=this_module, procedure='test_kexptA_cdp')

        ! Check result.
        call Xkryl%sub(Xref) ; err = Xkryl%norm()
        if (verb) write(output_unit, *) "     True error: ||error||_2 =", err

        call check(error, err < rtol_dp)
        call check_test(error, 'test_kexptA_cdp', 'norm2(x - x_true) < rol')

        return
    end subroutine test_kexptA_cdp

    subroutine test_block_kexptA_cdp(error)
        !! This function tests the Krylov based approximation of the action of the exponential
        !! propagator against the dense computation for a random operator, a random RHS and a 
        !! typical value of tau.

        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        class(linop_cdp), allocatable :: A
        ! Basis vectors.
        class(vector_cdp), allocatable :: B(:)
        class(vector_cdp), allocatable :: Cref(:)
        class(vector_cdp), allocatable :: C(:), Cblk(:)
        ! Krylov subspace dimension.
        integer, parameter :: kdim = test_size
        ! Test matrix.
        complex(dp), allocatable :: Adata(:, :), Edata(:, :)
        ! Information flag.
        integer :: info
        ! Test parameters
        integer       , parameter :: nkmax = 15
        integer       , parameter :: p = 3
        real(dp), parameter :: tau = 0.1_dp
        real(dp), parameter :: tol = rtol_dp
        logical       , parameter :: verb = .true.
        ! Misc.
        integer  :: i, j
        complex(dp), allocatable :: Xdata(:, :), Qdata(:, :)
        real(dp) :: alpha
        complex(dp) :: err(p, p)

        allocate(Adata(kdim, kdim)) ; Adata = 0.0_dp
        allocate(Edata(kdim, kdim)) ; Edata = 0.0_dp
        allocate(Xdata(test_size, p)) ; Xdata = 0.0_dp

        allocate(Cref(p)) ; call initialize_krylov_subspace(Cref)
        allocate(C(p))    ; call initialize_krylov_subspace(C)
        allocate(Cblk(p)) ; call initialize_krylov_subspace(Cblk)

        ! --> Initialize operator.
        A = linop_cdp() ; call init_rand(A) ; call get_data(Adata, A)
       
        ! --> Initialize rhs.
        allocate(B(1:p)) ; call init_rand(B)
        allocate(Qdata(test_size, p)) ; call get_data(Qdata, B)

        ! Comparison is dense computation (10th order Pade approximation)
        call expm(Edata, tau*Adata) ; Xdata = matmul(Edata,Qdata) ; call put_data(Cref, Xdata)

        ! Compute Krylov matrix exponential using sequential arnoldi method for each input column
        if (verb) write(*,*) 'SEQUENTIAL ARNOLDI'
        do i = 1,p
            if (verb) write(*,*) '    column',i
            call kexpm(C(i), A, B(i), tau, tol, info, verbosity=verb, kdim=nkmax)
            call check_info(info, 'kexpm', module=this_module, procedure='test_block_kexptA_cdp, 1')
        end do
        
        ! Compute Krylov matrix exponential using block-arnoldi method
        if (verb) write(*,*) 'BLOCK-ARNOLDI'
        call kexpm(Cblk, A, B, tau, tol, info, verbosity=verb, kdim=nkmax)
        call check_info(info, 'kexpm', module=this_module, procedure='test_block_kexptA_cdp, 2')
    
        do i = 1, p
            write(output_unit, *) C(i)%norm(), Cblk(i)%norm()
            call C(i)%sub(Cref(i)) ; call Cblk(i)%sub(Cref(i))
        end do

        ! Compute 2-norm of the error
        if (verb) then
            do i = 1, size(C)
                do j = 1, size(C)
                    err(i, j) = C(i)%dot(C(j))
                enddo
            enddo
            alpha = norm2(abs(err))
            write(*,*) '--------------------------------------------------------------------'
            write(*, *) '    true error (seq.):   ||error||_2 = ', alpha
        endif
 
        do i = 1, size(Cblk)
            do j = 1, size(Cblk)
                err(i, j) = Cblk(i)%dot(Cblk(j))
            enddo
        enddo
        
        alpha = norm2(abs(err))
        if (verb) write(*, *) '    true error (block):  ||error||_2 = ', alpha

        call check(error, alpha < rtol_dp)
        call check_test(error, 'test_block_kexptA_cdp', 'norm2(x - x_true) < rol')

        return
    end subroutine test_block_kexptA_cdp

    !-----------------------------------------------------------
    !-----     UNIT TESTS FOR DENSE MATRIX SQUARE ROOT     -----
    !-----------------------------------------------------------

    subroutine collect_sqrtm_rsp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                        new_unittest("Dense sqrtm for positive definite matrices.", test_dense_sqrtm_pos_def_rsp), &
                        new_unittest("Dense sqrtm for positive semi-definite matrices.", test_dense_sqrtm_pos_semi_def_rsp) &
                    ]

        return
    end subroutine collect_sqrtm_rsp_testsuite

    subroutine test_dense_sqrtm_pos_def_rsp(error)
       !> This function tests the matrix version of the sqrt function for the case of
       ! a symmetric/hermitian positive definite matrix
    
       !> Error type to be returned.
       type(error_type), allocatable, intent(out) :: error
       !> Problem dimension.
       integer, parameter :: n = 5
       !> Test matrix.
       real(sp) :: A(n, n)
       real(sp) :: sqrtmA(n, n)
       complex(sp) :: lambda(n)
       integer :: i, info
    
       ! --> Initialize matrix.
       call random_number(A)
       ! make symmetric/hermitian positive definite
       A = 0.5_sp*(A + transpose(A))
       call eig(A, sqrtmA, lambda)
       do i = 1,n
          lambda(i) = abs(lambda(i)) + 0.1_sp
       end do
       ! reconstruct matrix
       A = matmul(sqrtmA, matmul(diag(abs(lambda)), transpose(sqrtmA)))
       ! ensure it is exactly symmetric/hermitian
       A = 0.5_sp*(A + transpose(A))
     
       ! compute matrix square root
       call sqrtm(A, sqrtmA, info)
       call check_info(info, 'sqrtm', module=this_module, procedure='test_dense_sqrtm_pos_def_rsp')
    
       write(*,*) 'max err: ', maxval(matmul(sqrtmA, sqrtmA) - A)
       call check(error, maxval(matmul(sqrtmA, sqrtmA) - A) < rtol_sp)
       call check_test(error, 'test_dense_sqrtm_pos_def_rsp', 'sqrt(A)**2 /= A')
    
       return
    end subroutine test_dense_sqrtm_pos_def_rsp
    
    subroutine test_dense_sqrtm_pos_semi_def_rsp(error)
       !> This function tests the matrix version of the sqrt function for the case 
       ! of a symmetric semi-definite matrix
    
       !> Error type to be returned.
       type(error_type), allocatable, intent(out) :: error
       !> Problem dimension.
       integer, parameter :: n = 5
       !> Test matrix.
       real(sp) :: A(n, n)
       real(sp) :: sqrtmA(n, n)
       complex(sp) :: lambda(n)
       integer :: i, info
    
       ! --> Initialize matrix.
       call random_number(A)
       ! make symmetric/hermitian positive semi-definite
       A = 0.5_sp*(A + transpose(A))
       call eig(A, sqrtmA, lambda)
       do i = 1,n-1
          lambda(i) = abs(lambda(i)) + 0.1_sp
       end do
       lambda(n) = zero_rsp
       ! reconstruct matrix
       A = matmul(sqrtmA, matmul(diag(abs(lambda)), transpose(sqrtmA)))
       ! ensure it is exactly symmetric/hermitian
       A = 0.5_sp*(A + transpose(A))
    
       ! compute matrix square root
       call sqrtm(A, sqrtmA, info)
       call check_info(info, 'sqrtm', module=this_module, procedure='test_dense_sqrtm_pos_semi_def_rsp')
    
       write(*,*) 'max err: ', maxval(matmul(sqrtmA, sqrtmA) - A)
       call check(error, maxval(matmul(sqrtmA, sqrtmA) - A) < rtol_sp)
       call check_test(error, 'test_dense_sqrtm_pos_semi_def_rsp', 'sqrt(A)**2 /= A')
    
       return
    end subroutine test_dense_sqrtm_pos_semi_def_rsp

    subroutine collect_sqrtm_rdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                        new_unittest("Dense sqrtm for positive definite matrices.", test_dense_sqrtm_pos_def_rdp), &
                        new_unittest("Dense sqrtm for positive semi-definite matrices.", test_dense_sqrtm_pos_semi_def_rdp) &
                    ]

        return
    end subroutine collect_sqrtm_rdp_testsuite

    subroutine test_dense_sqrtm_pos_def_rdp(error)
       !> This function tests the matrix version of the sqrt function for the case of
       ! a symmetric/hermitian positive definite matrix
    
       !> Error type to be returned.
       type(error_type), allocatable, intent(out) :: error
       !> Problem dimension.
       integer, parameter :: n = 5
       !> Test matrix.
       real(dp) :: A(n, n)
       real(dp) :: sqrtmA(n, n)
       complex(dp) :: lambda(n)
       integer :: i, info
    
       ! --> Initialize matrix.
       call random_number(A)
       ! make symmetric/hermitian positive definite
       A = 0.5_dp*(A + transpose(A))
       call eig(A, sqrtmA, lambda)
       do i = 1,n
          lambda(i) = abs(lambda(i)) + 0.1_dp
       end do
       ! reconstruct matrix
       A = matmul(sqrtmA, matmul(diag(abs(lambda)), transpose(sqrtmA)))
       ! ensure it is exactly symmetric/hermitian
       A = 0.5_dp*(A + transpose(A))
     
       ! compute matrix square root
       call sqrtm(A, sqrtmA, info)
       call check_info(info, 'sqrtm', module=this_module, procedure='test_dense_sqrtm_pos_def_rdp')
    
       write(*,*) 'max err: ', maxval(matmul(sqrtmA, sqrtmA) - A)
       call check(error, maxval(matmul(sqrtmA, sqrtmA) - A) < rtol_dp)
       call check_test(error, 'test_dense_sqrtm_pos_def_rdp', 'sqrt(A)**2 /= A')
    
       return
    end subroutine test_dense_sqrtm_pos_def_rdp
    
    subroutine test_dense_sqrtm_pos_semi_def_rdp(error)
       !> This function tests the matrix version of the sqrt function for the case 
       ! of a symmetric semi-definite matrix
    
       !> Error type to be returned.
       type(error_type), allocatable, intent(out) :: error
       !> Problem dimension.
       integer, parameter :: n = 5
       !> Test matrix.
       real(dp) :: A(n, n)
       real(dp) :: sqrtmA(n, n)
       complex(dp) :: lambda(n)
       integer :: i, info
    
       ! --> Initialize matrix.
       call random_number(A)
       ! make symmetric/hermitian positive semi-definite
       A = 0.5_dp*(A + transpose(A))
       call eig(A, sqrtmA, lambda)
       do i = 1,n-1
          lambda(i) = abs(lambda(i)) + 0.1_dp
       end do
       lambda(n) = zero_rdp
       ! reconstruct matrix
       A = matmul(sqrtmA, matmul(diag(abs(lambda)), transpose(sqrtmA)))
       ! ensure it is exactly symmetric/hermitian
       A = 0.5_dp*(A + transpose(A))
    
       ! compute matrix square root
       call sqrtm(A, sqrtmA, info)
       call check_info(info, 'sqrtm', module=this_module, procedure='test_dense_sqrtm_pos_semi_def_rdp')
    
       write(*,*) 'max err: ', maxval(matmul(sqrtmA, sqrtmA) - A)
       call check(error, maxval(matmul(sqrtmA, sqrtmA) - A) < rtol_dp)
       call check_test(error, 'test_dense_sqrtm_pos_semi_def_rdp', 'sqrt(A)**2 /= A')
    
       return
    end subroutine test_dense_sqrtm_pos_semi_def_rdp

    subroutine collect_sqrtm_csp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                        new_unittest("Dense sqrtm for positive definite matrices.", test_dense_sqrtm_pos_def_csp), &
                        new_unittest("Dense sqrtm for positive semi-definite matrices.", test_dense_sqrtm_pos_semi_def_csp) &
                    ]

        return
    end subroutine collect_sqrtm_csp_testsuite

    subroutine test_dense_sqrtm_pos_def_csp(error)
       !> This function tests the matrix version of the sqrt function for the case of
       ! a symmetric/hermitian positive definite matrix
    
       !> Error type to be returned.
       type(error_type), allocatable, intent(out) :: error
       !> Problem dimension.
       integer, parameter :: n = 5
       !> Test matrix.
       complex(sp) :: A(n, n)
       complex(sp) :: sqrtmA(n, n)
       complex(sp) :: lambda(n)
       integer :: i, info
    
       ! --> Initialize matrix.
       call random_number(A%re)
       call random_number(A%im)
       ! make symmetric/hermitian positive definite
       A = 0.5_sp*(A + conjg(transpose(A)))
       call eig(A, sqrtmA, lambda)
       do i = 1,n
          lambda(i) = abs(lambda(i)) + 0.1_sp
       end do
       ! reconstruct matrix
       A = matmul(sqrtmA, matmul(diag(abs(lambda)), conjg(transpose(sqrtmA))))
       ! ensure it is exactly symmetric/hermitian
       A = 0.5_sp*(A + conjg(transpose(A)))
     
       ! compute matrix square root
       call sqrtm(A, sqrtmA, info)
       call check_info(info, 'sqrtm', module=this_module, procedure='test_dense_sqrtm_pos_def_csp')
    
       write(*,*) 'max err: ', maxval(abs(matmul(sqrtmA, sqrtmA) - A))
       call check(error, maxval(abs(matmul(sqrtmA, sqrtmA) - A)) < rtol_sp)
       call check_test(error, 'test_dense_sqrtm_pos_def_csp', 'sqrt(A)**2 /= A')
    
       return
    end subroutine test_dense_sqrtm_pos_def_csp
    
    subroutine test_dense_sqrtm_pos_semi_def_csp(error)
       !> This function tests the matrix version of the sqrt function for the case 
       ! of a symmetric semi-definite matrix
    
       !> Error type to be returned.
       type(error_type), allocatable, intent(out) :: error
       !> Problem dimension.
       integer, parameter :: n = 5
       !> Test matrix.
       complex(sp) :: A(n, n)
       complex(sp) :: sqrtmA(n, n)
       complex(sp) :: lambda(n)
       integer :: i, info
    
       ! --> Initialize matrix.
       call random_number(A%re)
       call random_number(A%im)
       ! make symmetric/hermitian positive semi-definite
       A = 0.5_sp*(A + conjg(transpose(A)))
       call eig(A, sqrtmA, lambda)
       do i = 1,n-1
          lambda(i) = abs(lambda(i)) + 0.1_sp
       end do
       lambda(n) = zero_rsp
       ! reconstruct matrix
       A = matmul(sqrtmA, matmul(diag(abs(lambda)), conjg(transpose(sqrtmA))))
       ! ensure it is exactly symmetric/hermitian
       A = 0.5_sp*(A + conjg(transpose(A)))
    
       ! compute matrix square root
       call sqrtm(A, sqrtmA, info)
       call check_info(info, 'sqrtm', module=this_module, procedure='test_dense_sqrtm_pos_semi_def_csp')
    
       write(*,*) 'max err: ', maxval(abs(matmul(sqrtmA, sqrtmA) - A))
       call check(error, maxval(abs(matmul(sqrtmA, sqrtmA) - A)) < rtol_sp)
       call check_test(error, 'test_dense_sqrtm_pos_semi_def_csp', 'sqrt(A)**2 /= A')
    
       return
    end subroutine test_dense_sqrtm_pos_semi_def_csp

    subroutine collect_sqrtm_cdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                        new_unittest("Dense sqrtm for positive definite matrices.", test_dense_sqrtm_pos_def_cdp), &
                        new_unittest("Dense sqrtm for positive semi-definite matrices.", test_dense_sqrtm_pos_semi_def_cdp) &
                    ]

        return
    end subroutine collect_sqrtm_cdp_testsuite

    subroutine test_dense_sqrtm_pos_def_cdp(error)
       !> This function tests the matrix version of the sqrt function for the case of
       ! a symmetric/hermitian positive definite matrix
    
       !> Error type to be returned.
       type(error_type), allocatable, intent(out) :: error
       !> Problem dimension.
       integer, parameter :: n = 5
       !> Test matrix.
       complex(dp) :: A(n, n)
       complex(dp) :: sqrtmA(n, n)
       complex(dp) :: lambda(n)
       integer :: i, info
    
       ! --> Initialize matrix.
       call random_number(A%re)
       call random_number(A%im)
       ! make symmetric/hermitian positive definite
       A = 0.5_dp*(A + conjg(transpose(A)))
       call eig(A, sqrtmA, lambda)
       do i = 1,n
          lambda(i) = abs(lambda(i)) + 0.1_dp
       end do
       ! reconstruct matrix
       A = matmul(sqrtmA, matmul(diag(abs(lambda)), conjg(transpose(sqrtmA))))
       ! ensure it is exactly symmetric/hermitian
       A = 0.5_dp*(A + conjg(transpose(A)))
     
       ! compute matrix square root
       call sqrtm(A, sqrtmA, info)
       call check_info(info, 'sqrtm', module=this_module, procedure='test_dense_sqrtm_pos_def_cdp')
    
       write(*,*) 'max err: ', maxval(abs(matmul(sqrtmA, sqrtmA) - A))
       call check(error, maxval(abs(matmul(sqrtmA, sqrtmA) - A)) < rtol_dp)
       call check_test(error, 'test_dense_sqrtm_pos_def_cdp', 'sqrt(A)**2 /= A')
    
       return
    end subroutine test_dense_sqrtm_pos_def_cdp
    
    subroutine test_dense_sqrtm_pos_semi_def_cdp(error)
       !> This function tests the matrix version of the sqrt function for the case 
       ! of a symmetric semi-definite matrix
    
       !> Error type to be returned.
       type(error_type), allocatable, intent(out) :: error
       !> Problem dimension.
       integer, parameter :: n = 5
       !> Test matrix.
       complex(dp) :: A(n, n)
       complex(dp) :: sqrtmA(n, n)
       complex(dp) :: lambda(n)
       integer :: i, info
    
       ! --> Initialize matrix.
       call random_number(A%re)
       call random_number(A%im)
       ! make symmetric/hermitian positive semi-definite
       A = 0.5_dp*(A + conjg(transpose(A)))
       call eig(A, sqrtmA, lambda)
       do i = 1,n-1
          lambda(i) = abs(lambda(i)) + 0.1_dp
       end do
       lambda(n) = zero_rdp
       ! reconstruct matrix
       A = matmul(sqrtmA, matmul(diag(abs(lambda)), conjg(transpose(sqrtmA))))
       ! ensure it is exactly symmetric/hermitian
       A = 0.5_dp*(A + conjg(transpose(A)))
    
       ! compute matrix square root
       call sqrtm(A, sqrtmA, info)
       call check_info(info, 'sqrtm', module=this_module, procedure='test_dense_sqrtm_pos_semi_def_cdp')
    
       write(*,*) 'max err: ', maxval(abs(matmul(sqrtmA, sqrtmA) - A))
       call check(error, maxval(abs(matmul(sqrtmA, sqrtmA) - A)) < rtol_dp)
       call check_test(error, 'test_dense_sqrtm_pos_semi_def_cdp', 'sqrt(A)**2 /= A')
    
       return
    end subroutine test_dense_sqrtm_pos_semi_def_cdp


end module

