module TestLinops
    ! Fortran Standard Library
    use stdlib_math, only: is_close, all_close
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

    character*128, parameter, private :: this_module = 'LightKrylov_TestLinops'

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
        type(vector_rsp), allocatable :: x, y
        ! Test LinOp.
        type(linop_rsp), allocatable :: A

        ! Initialize vectors.
        x = vector_rsp() ; call x%rand()
        y = vector_rsp() ; call y%rand()

        ! Initialize matrix.
        A = linop_rsp()
        call random_number(A%data)

        ! Compute matrix-vector product.
        call A%matvec(x, y)

        ! Check result.
        call check(error, norm2(y%data - matmul(A%data, x%data)) < rtol_sp)
        call check_test(error, 'test_matvec_rsp', eq='norm2(y - A @ x)')
        
        return
    end subroutine test_matvec_rsp

    subroutine test_rmatvec_rsp(error)
        ! Error to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        type(vector_rsp), allocatable :: x, y
        ! Test LinOp.
        type(linop_rsp), allocatable :: A

        ! Initialize vectors.
        x = vector_rsp() ; call x%rand()
        y = vector_rsp() ; call y%rand()

        ! Initialize matrix.
        A = linop_rsp()
        call random_number(A%data)

        ! Compute matrix-vector product.
        call A%rmatvec(x, y)

        ! Check result.
        call check(error, norm2(y%data - matmul(transpose(A%data), x%data)) < rtol_sp)
        call check_test(error, 'test_rmatvec_rsp', eq='norm2(y - A.T @ x)')
        
        return
    end subroutine test_rmatvec_rsp

    subroutine test_adjoint_matvec_rsp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        type(vector_rsp), allocatable :: x, y
        ! Test LinOp.
        type(linop_rsp), allocatable :: A
        type(adjoint_linop_rsp), allocatable :: B

        ! Internal variable.
        real(sp) :: alpha

        ! Initialize matrix.
        A = linop_rsp()
        call random_number(A%data)

        ! Adjoint operator.
        allocate(B) ; B%A = A

        ! Initialize vectors.
        x = vector_rsp() ; call x%rand()
        y = vector_rsp() ; call y%zero()

        ! Compute adjoint matrix-vector product
        call B%matvec(x, y)

        ! Check result.
        alpha = norm2(y%data - matmul(transpose(A%data), x%data))
        call check(error, alpha < rtol_sp)
        call check_test(error, 'test_adjoint_matvec_rsp', eq='norm2(y - A.T @ x)')

        

       return
    end subroutine test_adjoint_matvec_rsp

    subroutine test_adjoint_rmatvec_rsp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        type(vector_rsp), allocatable :: x, y
        ! Test LinOp.
        type(linop_rsp), allocatable :: A
        type(adjoint_linop_rsp), allocatable :: B

        ! Internal variable.
        real(sp) :: alpha

        ! Initialize matrix.
        A = linop_rsp()
        call random_number(A%data)

        ! Adjoint operator.
        allocate(B) ; B%A = A

        ! Initialize vectors.
        x = vector_rsp() ; call x%rand()
        y = vector_rsp() ; call y%zero()

        ! Compute adjoint matrix-vector product
        call B%rmatvec(x, y)

        ! Check result.
        alpha = norm2(y%data - matmul(A%data, x%data))
        call check(error, alpha < rtol_sp)
        call check_test(error, 'test_adjoint_rmatvec_rsp', eq='norm2(y - A @ x)')

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
        type(vector_rdp), allocatable :: x, y
        ! Test LinOp.
        type(linop_rdp), allocatable :: A

        ! Initialize vectors.
        x = vector_rdp() ; call x%rand()
        y = vector_rdp() ; call y%rand()

        ! Initialize matrix.
        A = linop_rdp()
        call random_number(A%data)

        ! Compute matrix-vector product.
        call A%matvec(x, y)

        ! Check result.
        call check(error, norm2(y%data - matmul(A%data, x%data)) < rtol_dp)
        call check_test(error, 'test_matvec_rdp', eq='norm2(y - A @ x)')
        
        return
    end subroutine test_matvec_rdp

    subroutine test_rmatvec_rdp(error)
        ! Error to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        type(vector_rdp), allocatable :: x, y
        ! Test LinOp.
        type(linop_rdp), allocatable :: A

        ! Initialize vectors.
        x = vector_rdp() ; call x%rand()
        y = vector_rdp() ; call y%rand()

        ! Initialize matrix.
        A = linop_rdp()
        call random_number(A%data)

        ! Compute matrix-vector product.
        call A%rmatvec(x, y)

        ! Check result.
        call check(error, norm2(y%data - matmul(transpose(A%data), x%data)) < rtol_dp)
        call check_test(error, 'test_rmatvec_rdp', eq='norm2(y - A.T @ x)')
        
        return
    end subroutine test_rmatvec_rdp

    subroutine test_adjoint_matvec_rdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        type(vector_rdp), allocatable :: x, y
        ! Test LinOp.
        type(linop_rdp), allocatable :: A
        type(adjoint_linop_rdp), allocatable :: B

        ! Internal variable.
        real(dp) :: alpha

        ! Initialize matrix.
        A = linop_rdp()
        call random_number(A%data)

        ! Adjoint operator.
        allocate(B) ; B%A = A

        ! Initialize vectors.
        x = vector_rdp() ; call x%rand()
        y = vector_rdp() ; call y%zero()

        ! Compute adjoint matrix-vector product
        call B%matvec(x, y)

        ! Check result.
        alpha = norm2(y%data - matmul(transpose(A%data), x%data))
        call check(error, alpha < rtol_dp)
        call check_test(error, 'test_adjoint_matvec_rdp', eq='norm2(y - A.T @ x)')

        

       return
    end subroutine test_adjoint_matvec_rdp

    subroutine test_adjoint_rmatvec_rdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        type(vector_rdp), allocatable :: x, y
        ! Test LinOp.
        type(linop_rdp), allocatable :: A
        type(adjoint_linop_rdp), allocatable :: B

        ! Internal variable.
        real(dp) :: alpha

        ! Initialize matrix.
        A = linop_rdp()
        call random_number(A%data)

        ! Adjoint operator.
        allocate(B) ; B%A = A

        ! Initialize vectors.
        x = vector_rdp() ; call x%rand()
        y = vector_rdp() ; call y%zero()

        ! Compute adjoint matrix-vector product
        call B%rmatvec(x, y)

        ! Check result.
        alpha = norm2(y%data - matmul(A%data, x%data))
        call check(error, alpha < rtol_dp)
        call check_test(error, 'test_adjoint_rmatvec_rdp', eq='norm2(y - A @ x)')

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
        type(vector_csp), allocatable :: x, y
        ! Test LinOp.
        type(linop_csp), allocatable :: A
        real(sp), allocatable :: Adata(:, :, :)

        ! Initialize vectors.
        x = vector_csp() ; call x%rand()
        y = vector_csp() ; call y%rand()

        ! Initialize matrix.
        A = linop_csp()
        allocate(Adata(test_size, test_size, 2))
        call random_number(Adata) ; A%data%re = Adata(:, :, 1) ; A%data%im = Adata(:, :, 2)

        ! Compute matrix-vector product.
        call A%matvec(x, y)

        ! Check result.
        call check(error, norm2(abs(y%data - matmul(A%data, x%data))) < rtol_sp)
        call check_test(error, 'test_matvec_csp', eq='norm2(|y - A @ x|)')
        
        return
    end subroutine test_matvec_csp

    subroutine test_rmatvec_csp(error)
        ! Error to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        type(vector_csp), allocatable :: x, y
        ! Test LinOp.
        type(linop_csp), allocatable :: A
        real(sp), allocatable :: Adata(:, :, :)

        ! Initialize vectors.
        x = vector_csp() ; call x%rand()
        y = vector_csp() ; call y%rand()

        ! Initialize matrix.
        A = linop_csp()
        allocate(Adata(test_size, test_size, 2))
        call random_number(Adata) ; A%data%re = Adata(:, :, 1) ; A%data%im = Adata(:, :, 2)

        ! Compute matrix-vector product.
        call A%rmatvec(x, y)

        ! Check result.
        call check(error, norm2(abs(y%data - matmul(transpose(conjg(A%data)), x%data))) < rtol_sp)
        call check_test(error, 'test_rmatvec_csp', eq='norm2(|y - A.H @ x|)')
        
        return
    end subroutine test_rmatvec_csp

    subroutine test_adjoint_matvec_csp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        type(vector_csp), allocatable :: x, y
        ! Test LinOp.
        type(linop_csp), allocatable :: A
        type(adjoint_linop_csp), allocatable :: B
        real(sp), allocatable :: Adata(:, :, :)

        ! Internal variable.
        real(sp) :: alpha

        ! Initialize matrix.
        A = linop_csp()
        allocate(Adata(test_size, test_size, 2))
        call random_number(Adata) ; A%data%re = Adata(:, :, 1) ; A%data%im = Adata(:, :, 2)

        ! Adjoint operator.
        allocate(B) ; B%A = A

        ! Initialize vectors.
        x = vector_csp() ; call x%rand()
        y = vector_csp() ; call y%zero()

        ! Compute adjoint matrix-vector product
        call B%matvec(x, y)

        ! Check result.
        alpha = norm2(abs(y%data - matmul(transpose(conjg(A%data)), x%data)))
        call check(error, alpha < rtol_sp)
        call check_test(error, 'test_adjoint_matvec_csp', eq='norm2(|y - A.H @ x|)')

        

       return
    end subroutine test_adjoint_matvec_csp

    subroutine test_adjoint_rmatvec_csp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        type(vector_csp), allocatable :: x, y
        ! Test LinOp.
        type(linop_csp), allocatable :: A
        type(adjoint_linop_csp), allocatable :: B
        real(sp), allocatable :: Adata(:, :, :)

        ! Internal variable.
        real(sp) :: alpha

        ! Initialize matrix.
        A = linop_csp()
        allocate(Adata(test_size, test_size, 2))
        call random_number(Adata) ; A%data%re = Adata(:, :, 1) ; A%data%im = Adata(:, :, 2)

        ! Adjoint operator.
        allocate(B) ; B%A = A

        ! Initialize vectors.
        x = vector_csp() ; call x%rand()
        y = vector_csp() ; call y%zero()

        ! Compute adjoint matrix-vector product
        call B%rmatvec(x, y)

        ! Check result.
        alpha = norm2(abs(y%data - matmul(A%data, x%data)))
        call check(error, alpha < rtol_sp)
        call check_test(error, 'test_adjoint_rmatvec_csp', eq='norm2(|y - A @ x|)')

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
        type(vector_cdp), allocatable :: x, y
        ! Test LinOp.
        type(linop_cdp), allocatable :: A
        real(dp), allocatable :: Adata(:, :, :)

        ! Initialize vectors.
        x = vector_cdp() ; call x%rand()
        y = vector_cdp() ; call y%rand()

        ! Initialize matrix.
        A = linop_cdp()
        allocate(Adata(test_size, test_size, 2))
        call random_number(Adata) ; A%data%re = Adata(:, :, 1) ; A%data%im = Adata(:, :, 2)

        ! Compute matrix-vector product.
        call A%matvec(x, y)

        ! Check result.
        call check(error, norm2(abs(y%data - matmul(A%data, x%data))) < rtol_dp)
        call check_test(error, 'test_matvec_cdp', eq='norm2(|y - A @ x|)')
        
        return
    end subroutine test_matvec_cdp

    subroutine test_rmatvec_cdp(error)
        ! Error to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        type(vector_cdp), allocatable :: x, y
        ! Test LinOp.
        type(linop_cdp), allocatable :: A
        real(dp), allocatable :: Adata(:, :, :)

        ! Initialize vectors.
        x = vector_cdp() ; call x%rand()
        y = vector_cdp() ; call y%rand()

        ! Initialize matrix.
        A = linop_cdp()
        allocate(Adata(test_size, test_size, 2))
        call random_number(Adata) ; A%data%re = Adata(:, :, 1) ; A%data%im = Adata(:, :, 2)

        ! Compute matrix-vector product.
        call A%rmatvec(x, y)

        ! Check result.
        call check(error, norm2(abs(y%data - matmul(transpose(conjg(A%data)), x%data))) < rtol_dp)
        call check_test(error, 'test_rmatvec_cdp', eq='norm2(|y - A.H @ x|)')
        
        return
    end subroutine test_rmatvec_cdp

    subroutine test_adjoint_matvec_cdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        type(vector_cdp), allocatable :: x, y
        ! Test LinOp.
        type(linop_cdp), allocatable :: A
        type(adjoint_linop_cdp), allocatable :: B
        real(dp), allocatable :: Adata(:, :, :)

        ! Internal variable.
        real(dp) :: alpha

        ! Initialize matrix.
        A = linop_cdp()
        allocate(Adata(test_size, test_size, 2))
        call random_number(Adata) ; A%data%re = Adata(:, :, 1) ; A%data%im = Adata(:, :, 2)

        ! Adjoint operator.
        allocate(B) ; B%A = A

        ! Initialize vectors.
        x = vector_cdp() ; call x%rand()
        y = vector_cdp() ; call y%zero()

        ! Compute adjoint matrix-vector product
        call B%matvec(x, y)

        ! Check result.
        alpha = norm2(abs(y%data - matmul(transpose(conjg(A%data)), x%data)))
        call check(error, alpha < rtol_dp)
        call check_test(error, 'test_adjoint_matvec_cdp', eq='norm2(|y - A.H @ x|)')

        

       return
    end subroutine test_adjoint_matvec_cdp

    subroutine test_adjoint_rmatvec_cdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vectors.
        type(vector_cdp), allocatable :: x, y
        ! Test LinOp.
        type(linop_cdp), allocatable :: A
        type(adjoint_linop_cdp), allocatable :: B
        real(dp), allocatable :: Adata(:, :, :)

        ! Internal variable.
        real(dp) :: alpha

        ! Initialize matrix.
        A = linop_cdp()
        allocate(Adata(test_size, test_size, 2))
        call random_number(Adata) ; A%data%re = Adata(:, :, 1) ; A%data%im = Adata(:, :, 2)

        ! Adjoint operator.
        allocate(B) ; B%A = A

        ! Initialize vectors.
        x = vector_cdp() ; call x%rand()
        y = vector_cdp() ; call y%zero()

        ! Compute adjoint matrix-vector product
        call B%rmatvec(x, y)

        ! Check result.
        alpha = norm2(abs(y%data - matmul(A%data, x%data)))
        call check(error, alpha < rtol_dp)
        call check_test(error, 'test_adjoint_rmatvec_cdp', eq='norm2(|y - A @ x|)')

       return
    end subroutine test_adjoint_rmatvec_cdp

end module TestLinops

