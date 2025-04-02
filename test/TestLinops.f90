module TestLinops
    ! Fortran Standard Library
    use stdlib_math, only: is_close, all_close
    use stdlib_linalg, only: norm, hermitian
    ! Testdrive
    use testdrive, only: new_unittest, unittest_type, error_type, check
    ! LightKrylov
    use LightKrylov
    use LightKrylov_Constants
    use LightKrylov_Logger
    ! Test Utilities
    use LightKrylov_TestUtils

    implicit none
    
    private

    character(len=*), parameter, private :: this_module      = 'LK_TLinops'
    character(len=*), parameter, private :: this_module_long = 'LightKrylov_TestLinops'
    integer, parameter :: n = 128

    public :: collect_linop_rsp_testsuite
    public :: collect_linop_rdp_testsuite
    public :: collect_linop_csp_testsuite
    public :: collect_linop_cdp_testsuite

contains

    !---------------------------------------------------------
    !-----     DEFINITIONS OF THE VARIOUS UNIT TESTS     -----
    !---------------------------------------------------------

    subroutine collect_linop_rsp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("Matrix-vector product", test_matvec_rsp), &
                    new_unittest("Tranpose Matrix-vector product", test_rmatvec_rsp),  &
                    new_unittest("Adjoint Matrix-vector product", test_adjoint_matvec_rsp),  &
                    new_unittest("Tranpose Adjoint Matrix-vector product", test_adjoint_rmatvec_rsp) &
                    ]
        return
    end subroutine collect_linop_rsp_testsuite

    subroutine test_matvec_rsp(error)
        ! Error to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        type(dense_vector_rsp) :: x, y
        real(sp) :: x_(n), y_(n)
        ! Test LinOp.
        type(dense_linop_rsp) :: A

        ! Initialize vectors.
        x = dense_vector(x_) ; call x%rand()
        y = dense_vector(y_) ; call y%zero()

        ! Initialize matrix.
        A = dense_linop_rsp() ; allocate(A%data(n, n))
        call random_number(A%data)

        ! Compute matrix-vector product.
        call A%apply_matvec(x, y)

        ! Check result.
        call check(error, norm(y%data - matmul(A%data, x%data), 2) < rtol_sp)
        call check_test(error, 'test_matvec_rsp', eq='norm(y - A @ x)')
        
        return
    end subroutine test_matvec_rsp

    subroutine test_rmatvec_rsp(error)
        ! Error to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        type(dense_vector_rsp) :: x, y
        real(sp) :: x_(n), y_(n)
        ! Test LinOp.
        type(dense_linop_rsp) :: A

        ! Initialize vectors.
        x = dense_vector(x_) ; call x%rand()
        y = dense_vector(y_) ; call y%zero()

        ! Initialize matrix.
        allocate(A%data(n, n))
        call random_number(A%data)

        ! Compute matrix-vector product.
        call A%apply_rmatvec(x, y)

        ! Check result.
        call check(error, norm(y%data - matmul(hermitian(A%data), x%data), 2) < rtol_sp)
        call check_test(error, 'test_rmatvec_rsp', eq='norm(y - A.H @ x)')
       
        return
    end subroutine test_rmatvec_rsp

    subroutine test_adjoint_matvec_rsp(error)
        ! Error to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        type(dense_vector_rsp) :: x, y
        real(sp) :: x_(n), y_(n)
        ! Test LinOp.
        type(dense_linop_rsp) :: A
        type(adjoint_linop_rsp), allocatable :: B

        ! Initialize vectors.
        x = dense_vector(x_) ; call x%rand()
        y = dense_vector(y_) ; call y%zero()

        ! Initialize matrix.
        allocate(A%data(n, n))
        call random_number(A%data)

        ! Compute matrix-vector product.
        B = adjoint(A) ; call B%apply_matvec(x, y)

        ! Check result.
        call check(error, norm(y%data - matmul(hermitian(A%data), x%data), 2) < rtol_sp)
        call check_test(error, 'test_adjoint_matvec_rsp', eq='norm(y - adjoint(A) @ x)')

       return
    end subroutine test_adjoint_matvec_rsp

    subroutine test_adjoint_rmatvec_rsp(error)
        ! Error to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        type(dense_vector_rsp) :: x, y
        real(sp) :: x_(n), y_(n)
        ! Test LinOp.
        type(dense_linop_rsp) :: A
        type(adjoint_linop_rsp), allocatable :: B

        ! Initialize vectors.
        x = dense_vector(x_) ; call x%rand()
        y = dense_vector(y_) ; call y%zero()

        ! Initialize matrix.
        allocate(A%data(n, n))
        call random_number(A%data)

        ! Compute matrix-vector product.
        B = adjoint(A) ; call B%apply_rmatvec(x, y)

        ! Check result.
        call check(error, norm(y%data - matmul(A%data, x%data), 2) < rtol_sp)
        call check_test(error, 'test_adjoint_rmatvec_rsp', eq='norm(y - A @ x)')

       return
    end subroutine test_adjoint_rmatvec_rsp

    subroutine collect_linop_rdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("Matrix-vector product", test_matvec_rdp), &
                    new_unittest("Tranpose Matrix-vector product", test_rmatvec_rdp),  &
                    new_unittest("Adjoint Matrix-vector product", test_adjoint_matvec_rdp),  &
                    new_unittest("Tranpose Adjoint Matrix-vector product", test_adjoint_rmatvec_rdp) &
                    ]
        return
    end subroutine collect_linop_rdp_testsuite

    subroutine test_matvec_rdp(error)
        ! Error to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        type(dense_vector_rdp) :: x, y
        real(dp) :: x_(n), y_(n)
        ! Test LinOp.
        type(dense_linop_rdp) :: A

        ! Initialize vectors.
        x = dense_vector(x_) ; call x%rand()
        y = dense_vector(y_) ; call y%zero()

        ! Initialize matrix.
        A = dense_linop_rdp() ; allocate(A%data(n, n))
        call random_number(A%data)

        ! Compute matrix-vector product.
        call A%apply_matvec(x, y)

        ! Check result.
        call check(error, norm(y%data - matmul(A%data, x%data), 2) < rtol_dp)
        call check_test(error, 'test_matvec_rdp', eq='norm(y - A @ x)')
        
        return
    end subroutine test_matvec_rdp

    subroutine test_rmatvec_rdp(error)
        ! Error to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        type(dense_vector_rdp) :: x, y
        real(dp) :: x_(n), y_(n)
        ! Test LinOp.
        type(dense_linop_rdp) :: A

        ! Initialize vectors.
        x = dense_vector(x_) ; call x%rand()
        y = dense_vector(y_) ; call y%zero()

        ! Initialize matrix.
        allocate(A%data(n, n))
        call random_number(A%data)

        ! Compute matrix-vector product.
        call A%apply_rmatvec(x, y)

        ! Check result.
        call check(error, norm(y%data - matmul(hermitian(A%data), x%data), 2) < rtol_dp)
        call check_test(error, 'test_rmatvec_rdp', eq='norm(y - A.H @ x)')
       
        return
    end subroutine test_rmatvec_rdp

    subroutine test_adjoint_matvec_rdp(error)
        ! Error to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        type(dense_vector_rdp) :: x, y
        real(dp) :: x_(n), y_(n)
        ! Test LinOp.
        type(dense_linop_rdp) :: A
        type(adjoint_linop_rdp), allocatable :: B

        ! Initialize vectors.
        x = dense_vector(x_) ; call x%rand()
        y = dense_vector(y_) ; call y%zero()

        ! Initialize matrix.
        allocate(A%data(n, n))
        call random_number(A%data)

        ! Compute matrix-vector product.
        B = adjoint(A) ; call B%apply_matvec(x, y)

        ! Check result.
        call check(error, norm(y%data - matmul(hermitian(A%data), x%data), 2) < rtol_dp)
        call check_test(error, 'test_adjoint_matvec_rdp', eq='norm(y - adjoint(A) @ x)')

       return
    end subroutine test_adjoint_matvec_rdp

    subroutine test_adjoint_rmatvec_rdp(error)
        ! Error to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        type(dense_vector_rdp) :: x, y
        real(dp) :: x_(n), y_(n)
        ! Test LinOp.
        type(dense_linop_rdp) :: A
        type(adjoint_linop_rdp), allocatable :: B

        ! Initialize vectors.
        x = dense_vector(x_) ; call x%rand()
        y = dense_vector(y_) ; call y%zero()

        ! Initialize matrix.
        allocate(A%data(n, n))
        call random_number(A%data)

        ! Compute matrix-vector product.
        B = adjoint(A) ; call B%apply_rmatvec(x, y)

        ! Check result.
        call check(error, norm(y%data - matmul(A%data, x%data), 2) < rtol_dp)
        call check_test(error, 'test_adjoint_rmatvec_rdp', eq='norm(y - A @ x)')

       return
    end subroutine test_adjoint_rmatvec_rdp

    subroutine collect_linop_csp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("Matrix-vector product", test_matvec_csp), &
                    new_unittest("Tranpose Matrix-vector product", test_rmatvec_csp),  &
                    new_unittest("Adjoint Matrix-vector product", test_adjoint_matvec_csp),  &
                    new_unittest("Tranpose Adjoint Matrix-vector product", test_adjoint_rmatvec_csp) &
                    ]
        return
    end subroutine collect_linop_csp_testsuite

    subroutine test_matvec_csp(error)
        ! Error to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        type(dense_vector_csp) :: x, y
        complex(sp) :: x_(n), y_(n)
        ! Test LinOp.
        type(dense_linop_csp) :: A
        real(sp) :: Adata(n, n, 2)

        ! Initialize vectors.
        x = dense_vector(x_) ; call x%rand()
        y = dense_vector(y_) ; call y%zero()

        ! Initialize matrix.
        A = dense_linop_csp() ; allocate(A%data(n, n))
        call random_number(Adata) ; A%data%re = Adata(:, :, 1) ; A%data%im = Adata(:, :, 2)

        ! Compute matrix-vector product.
        call A%apply_matvec(x, y)

        ! Check result.
        call check(error, norm(y%data - matmul(A%data, x%data), 2) < rtol_sp)
        call check_test(error, 'test_matvec_csp', eq='norm(y - A @ x)')
        
        return
    end subroutine test_matvec_csp

    subroutine test_rmatvec_csp(error)
        ! Error to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        type(dense_vector_csp) :: x, y
        complex(sp) :: x_(n), y_(n)
        ! Test LinOp.
        type(dense_linop_csp) :: A
        real(sp) :: Adata(n, n, 2)

        ! Initialize vectors.
        x = dense_vector(x_) ; call x%rand()
        y = dense_vector(y_) ; call y%zero()

        ! Initialize matrix.
        allocate(A%data(n, n))
        call random_number(Adata) ; A%data%re = Adata(:, :, 1) ; A%data%im = Adata(:, :, 2)

        ! Compute matrix-vector product.
        call A%apply_rmatvec(x, y)

        ! Check result.
        call check(error, norm(y%data - matmul(hermitian(A%data), x%data), 2) < rtol_sp)
        call check_test(error, 'test_rmatvec_csp', eq='norm(y - A.H @ x)')
       
        return
    end subroutine test_rmatvec_csp

    subroutine test_adjoint_matvec_csp(error)
        ! Error to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        type(dense_vector_csp) :: x, y
        complex(sp) :: x_(n), y_(n)
        ! Test LinOp.
        type(dense_linop_csp) :: A
        type(adjoint_linop_csp), allocatable :: B
        real(sp) :: Adata(n, n, 2)

        ! Initialize vectors.
        x = dense_vector(x_) ; call x%rand()
        y = dense_vector(y_) ; call y%zero()

        ! Initialize matrix.
        allocate(A%data(n, n))
        call random_number(Adata) ; A%data%re = Adata(:, :, 1) ; A%data%im = Adata(:, :, 2)

        ! Compute matrix-vector product.
        B = adjoint(A) ; call B%apply_matvec(x, y)

        ! Check result.
        call check(error, norm(y%data - matmul(hermitian(A%data), x%data), 2) < rtol_sp)
        call check_test(error, 'test_adjoint_matvec_csp', eq='norm(y - adjoint(A) @ x)')

       return
    end subroutine test_adjoint_matvec_csp

    subroutine test_adjoint_rmatvec_csp(error)
        ! Error to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        type(dense_vector_csp) :: x, y
        complex(sp) :: x_(n), y_(n)
        ! Test LinOp.
        type(dense_linop_csp) :: A
        type(adjoint_linop_csp), allocatable :: B
        real(sp) :: Adata(n, n, 2)

        ! Initialize vectors.
        x = dense_vector(x_) ; call x%rand()
        y = dense_vector(y_) ; call y%zero()

        ! Initialize matrix.
        allocate(A%data(n, n))
        call random_number(Adata) ; A%data%re = Adata(:, :, 1) ; A%data%im = Adata(:, :, 2)

        ! Compute matrix-vector product.
        B = adjoint(A) ; call B%apply_rmatvec(x, y)

        ! Check result.
        call check(error, norm(y%data - matmul(A%data, x%data), 2) < rtol_sp)
        call check_test(error, 'test_adjoint_rmatvec_csp', eq='norm(y - A @ x)')

       return
    end subroutine test_adjoint_rmatvec_csp

    subroutine collect_linop_cdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("Matrix-vector product", test_matvec_cdp), &
                    new_unittest("Tranpose Matrix-vector product", test_rmatvec_cdp),  &
                    new_unittest("Adjoint Matrix-vector product", test_adjoint_matvec_cdp),  &
                    new_unittest("Tranpose Adjoint Matrix-vector product", test_adjoint_rmatvec_cdp) &
                    ]
        return
    end subroutine collect_linop_cdp_testsuite

    subroutine test_matvec_cdp(error)
        ! Error to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        type(dense_vector_cdp) :: x, y
        complex(dp) :: x_(n), y_(n)
        ! Test LinOp.
        type(dense_linop_cdp) :: A
        real(dp) :: Adata(n, n, 2)

        ! Initialize vectors.
        x = dense_vector(x_) ; call x%rand()
        y = dense_vector(y_) ; call y%zero()

        ! Initialize matrix.
        A = dense_linop_cdp() ; allocate(A%data(n, n))
        call random_number(Adata) ; A%data%re = Adata(:, :, 1) ; A%data%im = Adata(:, :, 2)

        ! Compute matrix-vector product.
        call A%apply_matvec(x, y)

        ! Check result.
        call check(error, norm(y%data - matmul(A%data, x%data), 2) < rtol_dp)
        call check_test(error, 'test_matvec_cdp', eq='norm(y - A @ x)')
        
        return
    end subroutine test_matvec_cdp

    subroutine test_rmatvec_cdp(error)
        ! Error to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        type(dense_vector_cdp) :: x, y
        complex(dp) :: x_(n), y_(n)
        ! Test LinOp.
        type(dense_linop_cdp) :: A
        real(dp) :: Adata(n, n, 2)

        ! Initialize vectors.
        x = dense_vector(x_) ; call x%rand()
        y = dense_vector(y_) ; call y%zero()

        ! Initialize matrix.
        allocate(A%data(n, n))
        call random_number(Adata) ; A%data%re = Adata(:, :, 1) ; A%data%im = Adata(:, :, 2)

        ! Compute matrix-vector product.
        call A%apply_rmatvec(x, y)

        ! Check result.
        call check(error, norm(y%data - matmul(hermitian(A%data), x%data), 2) < rtol_dp)
        call check_test(error, 'test_rmatvec_cdp', eq='norm(y - A.H @ x)')
       
        return
    end subroutine test_rmatvec_cdp

    subroutine test_adjoint_matvec_cdp(error)
        ! Error to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        type(dense_vector_cdp) :: x, y
        complex(dp) :: x_(n), y_(n)
        ! Test LinOp.
        type(dense_linop_cdp) :: A
        type(adjoint_linop_cdp), allocatable :: B
        real(dp) :: Adata(n, n, 2)

        ! Initialize vectors.
        x = dense_vector(x_) ; call x%rand()
        y = dense_vector(y_) ; call y%zero()

        ! Initialize matrix.
        allocate(A%data(n, n))
        call random_number(Adata) ; A%data%re = Adata(:, :, 1) ; A%data%im = Adata(:, :, 2)

        ! Compute matrix-vector product.
        B = adjoint(A) ; call B%apply_matvec(x, y)

        ! Check result.
        call check(error, norm(y%data - matmul(hermitian(A%data), x%data), 2) < rtol_dp)
        call check_test(error, 'test_adjoint_matvec_cdp', eq='norm(y - adjoint(A) @ x)')

       return
    end subroutine test_adjoint_matvec_cdp

    subroutine test_adjoint_rmatvec_cdp(error)
        ! Error to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        type(dense_vector_cdp) :: x, y
        complex(dp) :: x_(n), y_(n)
        ! Test LinOp.
        type(dense_linop_cdp) :: A
        type(adjoint_linop_cdp), allocatable :: B
        real(dp) :: Adata(n, n, 2)

        ! Initialize vectors.
        x = dense_vector(x_) ; call x%rand()
        y = dense_vector(y_) ; call y%zero()

        ! Initialize matrix.
        allocate(A%data(n, n))
        call random_number(Adata) ; A%data%re = Adata(:, :, 1) ; A%data%im = Adata(:, :, 2)

        ! Compute matrix-vector product.
        B = adjoint(A) ; call B%apply_rmatvec(x, y)

        ! Check result.
        call check(error, norm(y%data - matmul(A%data, x%data), 2) < rtol_dp)
        call check_test(error, 'test_adjoint_rmatvec_cdp', eq='norm(y - A @ x)')

       return
    end subroutine test_adjoint_rmatvec_cdp

end module TestLinops

