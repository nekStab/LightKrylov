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

        ! Check result.
        call Xkryl%sub(Xref) ; err = Xkryl%norm()
        if (verb) write(output_unit, *) "     True error: ||error||_2 =", err

       call check(error, err < rtol_sp)

        return
    end subroutine test_kexptA_rsp

    subroutine test_block_kexptA_rsp(error)
        !> This function tests the Krylov based approximation of the action of the exponential
        ! propagator against the dense computation for a random operator, a random RHS and a 
        ! typical value of tau.

        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        class(linop_rsp), allocatable :: A
        !> Basis vectors.
        class(vector_rsp), allocatable :: B(:)
        class(vector_rsp), allocatable :: Cref(:)
        class(vector_rsp), allocatable :: C(:), Cblk(:)
        !> Krylov subspace dimension.
        integer, parameter :: kdim = test_size
        !> Test matrix.
        real(sp) :: Adata(kdim, kdim)
        real(sp) :: Edata(kdim, kdim)
        !> GS factors.
        real(sp) :: R(kdim, kdim)
        real(sp) :: Id(kdim, kdim)
        !> Information flag.
        integer :: info
        !> Test parameters
        integer       , parameter :: nkmax = 15
        integer       , parameter :: p = 3
        real(sp), parameter :: tau = 0.1_sp
        real(sp), parameter :: tol = rtol_sp
        logical       , parameter :: verb = .true.
        !> Misc.
        integer  :: i, j, k
        real(sp) :: Xdata(test_size,p), Qdata(test_size,p)
        real(sp) :: alpha
        real(sp) :: err(p, p)

        Adata = 0.0_sp ; Edata = 0.0_sp ; Xdata = 0.0_sp

        allocate(Cref(p)) ; call initialize_krylov_subspace(Cref)
        allocate(C(p))    ; call initialize_krylov_subspace(C)
        allocate(Cblk(p)) ; call initialize_krylov_subspace(Cblk)

        ! --> Initialize operator.
        A = linop_rsp() ; call init_rand(A) ; call get_data(Adata, A)
       
        ! --> Initialize rhs.
        allocate(B(1:p)) ; call init_rand(B) ; call get_data(Qdata, B)

        !> Comparison is dense computation (10th order Pade approximation)
        call expm(Edata, tau*Adata) ; Xdata = matmul(Edata,Qdata) ; call put_data(Cref, Xdata)

        !> Compute Krylov matrix exponential using sequential arnoldi method for each input column
        if (verb) write(*,*) 'SEQUENTIAL ARNOLDI'
        do i = 1,p
            if (verb) write(*,*) '    column',i
            call kexpm(C(i), A, B(i), tau, tol, info, verbosity=verb, kdim=nkmax)
        end do
        
        !> Compute Krylov matrix exponential using block-arnoldi method
        if (verb) write(*,*) 'BLOCK-ARNOLDI'
        call kexpm(Cblk, A, B, tau, tol, info, verbosity=verb, kdim=nkmax)
        do i = 1, p
            write(output_unit, *) C(i)%norm(), Cblk(i)%norm()
            call C(i)%sub(Cref(i)) ; call Cblk(i)%sub(Cref(i))
        end do

        !> Compute 2-norm of the error
        if (verb) then
            do i = 1, size(C)
                do j = 1, size(C)
                    err(i, j) = C(i)%dot(C(j))
                enddo
            enddo
            alpha = sqrt(norm2(abs(err)))
            write(*,*) '--------------------------------------------------------------------'
            write(*, *) '    true error (seq.):   ||error||_2 = ', alpha
        endif
 
        do i = 1, size(Cblk)
            do j = 1, size(Cblk)
                err(i, j) = Cblk(i)%dot(Cblk(j))
            enddo
        enddo
        
        alpha = sqrt(norm2(abs(err)))
        if (verb) write(*, *) '    true error (block):  ||error||_2 = ', alpha

        call check(error, alpha < rtol_sp)

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

        ! Check result.
        call Xkryl%sub(Xref) ; err = Xkryl%norm()
        if (verb) write(output_unit, *) "     True error: ||error||_2 =", err

       call check(error, err < rtol_dp)

        return
    end subroutine test_kexptA_rdp

    subroutine test_block_kexptA_rdp(error)
        !> This function tests the Krylov based approximation of the action of the exponential
        ! propagator against the dense computation for a random operator, a random RHS and a 
        ! typical value of tau.

        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        class(linop_rdp), allocatable :: A
        !> Basis vectors.
        class(vector_rdp), allocatable :: B(:)
        class(vector_rdp), allocatable :: Cref(:)
        class(vector_rdp), allocatable :: C(:), Cblk(:)
        !> Krylov subspace dimension.
        integer, parameter :: kdim = test_size
        !> Test matrix.
        real(dp) :: Adata(kdim, kdim)
        real(dp) :: Edata(kdim, kdim)
        !> GS factors.
        real(dp) :: R(kdim, kdim)
        real(dp) :: Id(kdim, kdim)
        !> Information flag.
        integer :: info
        !> Test parameters
        integer       , parameter :: nkmax = 15
        integer       , parameter :: p = 3
        real(dp), parameter :: tau = 0.1_dp
        real(dp), parameter :: tol = rtol_dp
        logical       , parameter :: verb = .true.
        !> Misc.
        integer  :: i, j, k
        real(dp) :: Xdata(test_size,p), Qdata(test_size,p)
        real(dp) :: alpha
        real(dp) :: err(p, p)

        Adata = 0.0_dp ; Edata = 0.0_dp ; Xdata = 0.0_dp

        allocate(Cref(p)) ; call initialize_krylov_subspace(Cref)
        allocate(C(p))    ; call initialize_krylov_subspace(C)
        allocate(Cblk(p)) ; call initialize_krylov_subspace(Cblk)

        ! --> Initialize operator.
        A = linop_rdp() ; call init_rand(A) ; call get_data(Adata, A)
       
        ! --> Initialize rhs.
        allocate(B(1:p)) ; call init_rand(B) ; call get_data(Qdata, B)

        !> Comparison is dense computation (10th order Pade approximation)
        call expm(Edata, tau*Adata) ; Xdata = matmul(Edata,Qdata) ; call put_data(Cref, Xdata)

        !> Compute Krylov matrix exponential using sequential arnoldi method for each input column
        if (verb) write(*,*) 'SEQUENTIAL ARNOLDI'
        do i = 1,p
            if (verb) write(*,*) '    column',i
            call kexpm(C(i), A, B(i), tau, tol, info, verbosity=verb, kdim=nkmax)
        end do
        
        !> Compute Krylov matrix exponential using block-arnoldi method
        if (verb) write(*,*) 'BLOCK-ARNOLDI'
        call kexpm(Cblk, A, B, tau, tol, info, verbosity=verb, kdim=nkmax)
        do i = 1, p
            write(output_unit, *) C(i)%norm(), Cblk(i)%norm()
            call C(i)%sub(Cref(i)) ; call Cblk(i)%sub(Cref(i))
        end do

        !> Compute 2-norm of the error
        if (verb) then
            do i = 1, size(C)
                do j = 1, size(C)
                    err(i, j) = C(i)%dot(C(j))
                enddo
            enddo
            alpha = sqrt(norm2(abs(err)))
            write(*,*) '--------------------------------------------------------------------'
            write(*, *) '    true error (seq.):   ||error||_2 = ', alpha
        endif
 
        do i = 1, size(Cblk)
            do j = 1, size(Cblk)
                err(i, j) = Cblk(i)%dot(Cblk(j))
            enddo
        enddo
        
        alpha = sqrt(norm2(abs(err)))
        if (verb) write(*, *) '    true error (block):  ||error||_2 = ', alpha

        call check(error, alpha < rtol_dp)

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

        ! Check result.
        call Xkryl%sub(Xref) ; err = Xkryl%norm()
        if (verb) write(output_unit, *) "     True error: ||error||_2 =", err

       call check(error, err < rtol_sp)

        return
    end subroutine test_kexptA_csp

    subroutine test_block_kexptA_csp(error)
        !> This function tests the Krylov based approximation of the action of the exponential
        ! propagator against the dense computation for a random operator, a random RHS and a 
        ! typical value of tau.

        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        class(linop_csp), allocatable :: A
        !> Basis vectors.
        class(vector_csp), allocatable :: B(:)
        class(vector_csp), allocatable :: Cref(:)
        class(vector_csp), allocatable :: C(:), Cblk(:)
        !> Krylov subspace dimension.
        integer, parameter :: kdim = test_size
        !> Test matrix.
        complex(sp) :: Adata(kdim, kdim)
        complex(sp) :: Edata(kdim, kdim)
        !> GS factors.
        complex(sp) :: R(kdim, kdim)
        complex(sp) :: Id(kdim, kdim)
        !> Information flag.
        integer :: info
        !> Test parameters
        integer       , parameter :: nkmax = 15
        integer       , parameter :: p = 3
        real(sp), parameter :: tau = 0.1_sp
        real(sp), parameter :: tol = rtol_sp
        logical       , parameter :: verb = .true.
        !> Misc.
        integer  :: i, j, k
        complex(sp) :: Xdata(test_size,p), Qdata(test_size,p)
        real(sp) :: alpha
        complex(sp) :: err(p, p)

        Adata = 0.0_sp ; Edata = 0.0_sp ; Xdata = 0.0_sp

        allocate(Cref(p)) ; call initialize_krylov_subspace(Cref)
        allocate(C(p))    ; call initialize_krylov_subspace(C)
        allocate(Cblk(p)) ; call initialize_krylov_subspace(Cblk)

        ! --> Initialize operator.
        A = linop_csp() ; call init_rand(A) ; call get_data(Adata, A)
       
        ! --> Initialize rhs.
        allocate(B(1:p)) ; call init_rand(B) ; call get_data(Qdata, B)

        !> Comparison is dense computation (10th order Pade approximation)
        call expm(Edata, tau*Adata) ; Xdata = matmul(Edata,Qdata) ; call put_data(Cref, Xdata)

        !> Compute Krylov matrix exponential using sequential arnoldi method for each input column
        if (verb) write(*,*) 'SEQUENTIAL ARNOLDI'
        do i = 1,p
            if (verb) write(*,*) '    column',i
            call kexpm(C(i), A, B(i), tau, tol, info, verbosity=verb, kdim=nkmax)
        end do
        
        !> Compute Krylov matrix exponential using block-arnoldi method
        if (verb) write(*,*) 'BLOCK-ARNOLDI'
        call kexpm(Cblk, A, B, tau, tol, info, verbosity=verb, kdim=nkmax)
        do i = 1, p
            write(output_unit, *) C(i)%norm(), Cblk(i)%norm()
            call C(i)%sub(Cref(i)) ; call Cblk(i)%sub(Cref(i))
        end do

        !> Compute 2-norm of the error
        if (verb) then
            do i = 1, size(C)
                do j = 1, size(C)
                    err(i, j) = C(i)%dot(C(j))
                enddo
            enddo
            alpha = sqrt(norm2(abs(err)))
            write(*,*) '--------------------------------------------------------------------'
            write(*, *) '    true error (seq.):   ||error||_2 = ', alpha
        endif
 
        do i = 1, size(Cblk)
            do j = 1, size(Cblk)
                err(i, j) = Cblk(i)%dot(Cblk(j))
            enddo
        enddo
        
        alpha = sqrt(norm2(abs(err)))
        if (verb) write(*, *) '    true error (block):  ||error||_2 = ', alpha

        call check(error, alpha < rtol_sp)

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

        ! Check result.
        call Xkryl%sub(Xref) ; err = Xkryl%norm()
        if (verb) write(output_unit, *) "     True error: ||error||_2 =", err

       call check(error, err < rtol_dp)

        return
    end subroutine test_kexptA_cdp

    subroutine test_block_kexptA_cdp(error)
        !> This function tests the Krylov based approximation of the action of the exponential
        ! propagator against the dense computation for a random operator, a random RHS and a 
        ! typical value of tau.

        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        class(linop_cdp), allocatable :: A
        !> Basis vectors.
        class(vector_cdp), allocatable :: B(:)
        class(vector_cdp), allocatable :: Cref(:)
        class(vector_cdp), allocatable :: C(:), Cblk(:)
        !> Krylov subspace dimension.
        integer, parameter :: kdim = test_size
        !> Test matrix.
        complex(dp) :: Adata(kdim, kdim)
        complex(dp) :: Edata(kdim, kdim)
        !> GS factors.
        complex(dp) :: R(kdim, kdim)
        complex(dp) :: Id(kdim, kdim)
        !> Information flag.
        integer :: info
        !> Test parameters
        integer       , parameter :: nkmax = 15
        integer       , parameter :: p = 3
        real(dp), parameter :: tau = 0.1_dp
        real(dp), parameter :: tol = rtol_dp
        logical       , parameter :: verb = .true.
        !> Misc.
        integer  :: i, j, k
        complex(dp) :: Xdata(test_size,p), Qdata(test_size,p)
        real(dp) :: alpha
        complex(dp) :: err(p, p)

        Adata = 0.0_dp ; Edata = 0.0_dp ; Xdata = 0.0_dp

        allocate(Cref(p)) ; call initialize_krylov_subspace(Cref)
        allocate(C(p))    ; call initialize_krylov_subspace(C)
        allocate(Cblk(p)) ; call initialize_krylov_subspace(Cblk)

        ! --> Initialize operator.
        A = linop_cdp() ; call init_rand(A) ; call get_data(Adata, A)
       
        ! --> Initialize rhs.
        allocate(B(1:p)) ; call init_rand(B) ; call get_data(Qdata, B)

        !> Comparison is dense computation (10th order Pade approximation)
        call expm(Edata, tau*Adata) ; Xdata = matmul(Edata,Qdata) ; call put_data(Cref, Xdata)

        !> Compute Krylov matrix exponential using sequential arnoldi method for each input column
        if (verb) write(*,*) 'SEQUENTIAL ARNOLDI'
        do i = 1,p
            if (verb) write(*,*) '    column',i
            call kexpm(C(i), A, B(i), tau, tol, info, verbosity=verb, kdim=nkmax)
        end do
        
        !> Compute Krylov matrix exponential using block-arnoldi method
        if (verb) write(*,*) 'BLOCK-ARNOLDI'
        call kexpm(Cblk, A, B, tau, tol, info, verbosity=verb, kdim=nkmax)
        do i = 1, p
            write(output_unit, *) C(i)%norm(), Cblk(i)%norm()
            call C(i)%sub(Cref(i)) ; call Cblk(i)%sub(Cref(i))
        end do

        !> Compute 2-norm of the error
        if (verb) then
            do i = 1, size(C)
                do j = 1, size(C)
                    err(i, j) = C(i)%dot(C(j))
                enddo
            enddo
            alpha = sqrt(norm2(abs(err)))
            write(*,*) '--------------------------------------------------------------------'
            write(*, *) '    true error (seq.):   ||error||_2 = ', alpha
        endif
 
        do i = 1, size(Cblk)
            do j = 1, size(Cblk)
                err(i, j) = Cblk(i)%dot(Cblk(j))
            enddo
        enddo
        
        alpha = sqrt(norm2(abs(err)))
        if (verb) write(*, *) '    true error (block):  ||error||_2 = ', alpha

        call check(error, alpha < rtol_dp)

        return
    end subroutine test_block_kexptA_cdp


end module
