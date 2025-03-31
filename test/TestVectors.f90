module TestVectors
    ! Fortran Standard Library
    use stdlib_math, only: is_close, all_close
    use stdlib_linalg, only: norm
    use stdlib_stats_distribution_normal, only: normal => rvs_normal
    use stdlib_optval, only: optval
    ! Testdrive
    use testdrive, only: new_unittest, unittest_type, error_type, check
    ! LightKrylov
    use LightKrylov
    use LightKrylov_Logger
    ! Test Utilities
    use LightKrylov_TestUtils
    
    implicit none
    
    private

    character(len=*), parameter, private :: this_module      = 'LK_TVectors'
    character(len=*), parameter, private :: this_module_long = 'LightKrylov_TestVectors'
    integer, parameter :: n = 128

    public :: collect_vector_rsp_testsuite
    public :: collect_vector_rdp_testsuite
    public :: collect_vector_csp_testsuite
    public :: collect_vector_cdp_testsuite

contains
        
    !---------------------------------------------------------
    !-----     DEFINITIONS OF THE VARIOUS UNIT TESTS     -----
    !---------------------------------------------------------

    subroutine collect_vector_rsp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("Vector norm", test_vector_rsp_norm)      , &
                    new_unittest("Vector scale", test_vector_rsp_scal)     , &
                    new_unittest("Vector addition", test_vector_rsp_add)   , &
                    new_unittest("Vector subtraction", test_vector_rsp_sub), &
                    new_unittest("Vector dot product", test_vector_rsp_dot)  &
                    ]
        return
    end subroutine collect_vector_rsp_testsuite

    subroutine test_vector_rsp_norm(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vector.
        type(dense_vector_rsp) :: x
        real(sp) :: alpha

        ! Initialize vector.
        x = dense_vector_rsp(n) ; call x%rand()
        
        ! Compute its norm.
        alpha = x%norm()

        ! Check result.
        call check(error, is_close(alpha, norm(x%data, 2)))
        call check_test(error, 'test_vector_rsp_norm', eq='is_close(x%norm, norm(x, 2))')
        
        return
    end subroutine test_vector_rsp_norm

    subroutine test_vector_rsp_add(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error

        ! Test vectors.
        type(dense_vector_rsp), allocatable :: x, y, z

        ! Initialize vectors.
        x = dense_vector_rsp(n) ; call x%rand()
        y = dense_vector_rsp(n) ; call y%rand()
        z = x

        ! Vector addition.
        call z%add(y)

        ! Check correctness.
        call check(error, norm(z%data - x%data - y%data, 2) < rtol_sp)
        call check_test(error, 'test_vector_rsp_norm', eq='is_close(x%norm, norm(z - (x+y), 2))')

        return
    end subroutine test_vector_rsp_add
 
    subroutine test_vector_rsp_sub(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error

        ! Test vectors.
        type(dense_vector_rsp), allocatable :: x, y, z

        ! Initialize vectors.
        x = dense_vector_rsp(n) ; call x%rand()
        y = dense_vector_rsp(n) ; call y%rand()
        z = x

        ! Vector addition.
        call z%sub(y)

        ! Check correctness.
        call check(error, norm(z%data - (x%data - y%data), 2) < rtol_sp)
        call check_test(error, 'test_vector_rsp_norm', eq='is_close(x%norm, norm(z - (x-y), 2))')

        return
    end subroutine test_vector_rsp_sub

    subroutine test_vector_rsp_dot(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error

        ! Test vectors.
        type(dense_vector_rsp), allocatable :: x, y
        real(sp) :: alpha

        ! Initialize vectors.
        x = dense_vector_rsp(n) ; call x%rand()
        y = dense_vector_rsp(n) ; call y%rand()

        ! Compute inner-product.
        alpha = x%dot(y)

        ! Check correctness.
        call check(error, abs(alpha - dot_product(x%data, y%data)) < rtol_sp)
        call check_test(error, 'test_vector_rsp_dot', eq='abs(x%dot(y) - dot_product(x, y))')

        return
    end subroutine test_vector_rsp_dot

    subroutine test_vector_rsp_scal(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error

        ! Test vector.
        type(dense_vector_rsp), allocatable :: x, y
        real(sp) :: alpha

        ! Initialize vector.
        x = dense_vector_rsp(n) ; call x%rand(ifnorm=.true.)
        y = dense_vector_rsp(x%data)
        call random_number(alpha)
        print *, "alpha : ", alpha
        print *
        print *, "x     : ", x%data
        print *
        print *, "y     : ", y%data
        print *
        print *, "x-y     : ", x%data - y%data
        
        ! Scale the vector.
        call x%scal(alpha)
        print *
        print *, "scal(x) :", x%data
        print *
        print *, "scal(y) :", alpha*y%data
        print *
        print *, "diff    :", abs(x%data - alpha*y%data)
        print *
        print *, "norm(diff) :", norm(x%data - alpha*y%data, 2) 

        ! Check correctness.
        call check(error, norm(x%data - alpha*y%data, 2) < rtol_sp)
        call check_test(error, 'test_vector_rsp_scal', eq='norm(x - alpha*y, 2)')

        return
    end subroutine test_vector_rsp_scal

    subroutine collect_vector_rdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("Vector norm", test_vector_rdp_norm)      , &
                    new_unittest("Vector scale", test_vector_rdp_scal)     , &
                    new_unittest("Vector addition", test_vector_rdp_add)   , &
                    new_unittest("Vector subtraction", test_vector_rdp_sub), &
                    new_unittest("Vector dot product", test_vector_rdp_dot)  &
                    ]
        return
    end subroutine collect_vector_rdp_testsuite

    subroutine test_vector_rdp_norm(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vector.
        type(dense_vector_rdp) :: x
        real(dp) :: alpha

        ! Initialize vector.
        x = dense_vector_rdp(n) ; call x%rand()
        
        ! Compute its norm.
        alpha = x%norm()

        ! Check result.
        call check(error, is_close(alpha, norm(x%data, 2)))
        call check_test(error, 'test_vector_rdp_norm', eq='is_close(x%norm, norm(x, 2))')
        
        return
    end subroutine test_vector_rdp_norm

    subroutine test_vector_rdp_add(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error

        ! Test vectors.
        type(dense_vector_rdp), allocatable :: x, y, z

        ! Initialize vectors.
        x = dense_vector_rdp(n) ; call x%rand()
        y = dense_vector_rdp(n) ; call y%rand()
        z = x

        ! Vector addition.
        call z%add(y)

        ! Check correctness.
        call check(error, norm(z%data - x%data - y%data, 2) < rtol_dp)
        call check_test(error, 'test_vector_rdp_norm', eq='is_close(x%norm, norm(z - (x+y), 2))')

        return
    end subroutine test_vector_rdp_add
 
    subroutine test_vector_rdp_sub(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error

        ! Test vectors.
        type(dense_vector_rdp), allocatable :: x, y, z

        ! Initialize vectors.
        x = dense_vector_rdp(n) ; call x%rand()
        y = dense_vector_rdp(n) ; call y%rand()
        z = x

        ! Vector addition.
        call z%sub(y)

        ! Check correctness.
        call check(error, norm(z%data - (x%data - y%data), 2) < rtol_dp)
        call check_test(error, 'test_vector_rdp_norm', eq='is_close(x%norm, norm(z - (x-y), 2))')

        return
    end subroutine test_vector_rdp_sub

    subroutine test_vector_rdp_dot(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error

        ! Test vectors.
        type(dense_vector_rdp), allocatable :: x, y
        real(dp) :: alpha

        ! Initialize vectors.
        x = dense_vector_rdp(n) ; call x%rand()
        y = dense_vector_rdp(n) ; call y%rand()

        ! Compute inner-product.
        alpha = x%dot(y)

        ! Check correctness.
        call check(error, abs(alpha - dot_product(x%data, y%data)) < rtol_dp)
        call check_test(error, 'test_vector_rdp_dot', eq='abs(x%dot(y) - dot_product(x, y))')

        return
    end subroutine test_vector_rdp_dot

    subroutine test_vector_rdp_scal(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error

        ! Test vector.
        type(dense_vector_rdp), allocatable :: x, y
        real(dp) :: alpha

        ! Initialize vector.
        x = dense_vector_rdp(n) ; call x%rand(ifnorm=.true.)
        y = dense_vector_rdp(x%data)
        call random_number(alpha)
        print *, "alpha : ", alpha
        print *
        print *, "x     : ", x%data
        print *
        print *, "y     : ", y%data
        print *
        print *, "x-y     : ", x%data - y%data
        
        ! Scale the vector.
        call x%scal(alpha)
        print *
        print *, "scal(x) :", x%data
        print *
        print *, "scal(y) :", alpha*y%data
        print *
        print *, "diff    :", abs(x%data - alpha*y%data)
        print *
        print *, "norm(diff) :", norm(x%data - alpha*y%data, 2) 

        ! Check correctness.
        call check(error, norm(x%data - alpha*y%data, 2) < rtol_dp)
        call check_test(error, 'test_vector_rdp_scal', eq='norm(x - alpha*y, 2)')

        return
    end subroutine test_vector_rdp_scal

    subroutine collect_vector_csp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("Vector norm", test_vector_csp_norm)      , &
                    new_unittest("Vector scale", test_vector_csp_scal)     , &
                    new_unittest("Vector addition", test_vector_csp_add)   , &
                    new_unittest("Vector subtraction", test_vector_csp_sub), &
                    new_unittest("Vector dot product", test_vector_csp_dot)  &
                    ]
        return
    end subroutine collect_vector_csp_testsuite

    subroutine test_vector_csp_norm(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vector.
        type(dense_vector_csp) :: x
        real(sp) :: alpha

        ! Initialize vector.
        x = dense_vector_csp(n) ; call x%rand()
        
        ! Compute its norm.
        alpha = x%norm()

        ! Check result.
        call check(error, is_close(alpha, norm(x%data, 2)))
        call check_test(error, 'test_vector_csp_norm', eq='is_close(x%norm, norm(x, 2))')
        
        return
    end subroutine test_vector_csp_norm

    subroutine test_vector_csp_add(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error

        ! Test vectors.
        type(dense_vector_csp), allocatable :: x, y, z

        ! Initialize vectors.
        x = dense_vector_csp(n) ; call x%rand()
        y = dense_vector_csp(n) ; call y%rand()
        z = x

        ! Vector addition.
        call z%add(y)

        ! Check correctness.
        call check(error, norm(z%data - x%data - y%data, 2) < rtol_sp)
        call check_test(error, 'test_vector_csp_norm', eq='is_close(x%norm, norm(z - (x+y), 2))')

        return
    end subroutine test_vector_csp_add
 
    subroutine test_vector_csp_sub(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error

        ! Test vectors.
        type(dense_vector_csp), allocatable :: x, y, z

        ! Initialize vectors.
        x = dense_vector_csp(n) ; call x%rand()
        y = dense_vector_csp(n) ; call y%rand()
        z = x

        ! Vector addition.
        call z%sub(y)

        ! Check correctness.
        call check(error, norm(z%data - (x%data - y%data), 2) < rtol_sp)
        call check_test(error, 'test_vector_csp_norm', eq='is_close(x%norm, norm(z - (x-y), 2))')

        return
    end subroutine test_vector_csp_sub

    subroutine test_vector_csp_dot(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error

        ! Test vectors.
        type(dense_vector_csp), allocatable :: x, y
        complex(sp) :: alpha

        ! Initialize vectors.
        x = dense_vector_csp(n) ; call x%rand()
        y = dense_vector_csp(n) ; call y%rand()

        ! Compute inner-product.
        alpha = x%dot(y)

        ! Check correctness.
        call check(error, abs(alpha - dot_product(x%data, y%data)) < rtol_sp)
        call check_test(error, 'test_vector_csp_dot', eq='abs(x%dot(y) - dot_product(x, y))')

        return
    end subroutine test_vector_csp_dot

    subroutine test_vector_csp_scal(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error

        ! Test vector.
        type(dense_vector_csp), allocatable :: x, y
        complex(sp) :: alpha

        ! Initialize vector.
        x = dense_vector_csp(n) ; call x%rand(ifnorm=.true.)
        y = dense_vector_csp(x%data)
        alpha = 0.0_sp ; call random_number(alpha%re) ; call random_number(alpha%im)
        alpha = alpha / abs(alpha)
        print *, "alpha : ", alpha
        print *
        print *, "x     : ", x%data
        print *
        print *, "y     : ", y%data
        print *
        print *, "x-y     : ", x%data - y%data
        
        ! Scale the vector.
        call x%scal(alpha)
        print *
        print *, "scal(x) :", x%data
        print *
        print *, "scal(y) :", alpha*y%data
        print *
        print *, "diff    :", abs(x%data - alpha*y%data)
        print *
        print *, "norm(diff) :", norm(x%data - alpha*y%data, 2) 

        ! Check correctness.
        call check(error, norm(x%data - alpha*y%data, 2) < rtol_sp)
        call check_test(error, 'test_vector_csp_scal', eq='norm(x - alpha*y, 2)')

        return
    end subroutine test_vector_csp_scal

    subroutine collect_vector_cdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("Vector norm", test_vector_cdp_norm)      , &
                    new_unittest("Vector scale", test_vector_cdp_scal)     , &
                    new_unittest("Vector addition", test_vector_cdp_add)   , &
                    new_unittest("Vector subtraction", test_vector_cdp_sub), &
                    new_unittest("Vector dot product", test_vector_cdp_dot)  &
                    ]
        return
    end subroutine collect_vector_cdp_testsuite

    subroutine test_vector_cdp_norm(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vector.
        type(dense_vector_cdp) :: x
        real(dp) :: alpha

        ! Initialize vector.
        x = dense_vector_cdp(n) ; call x%rand()
        
        ! Compute its norm.
        alpha = x%norm()

        ! Check result.
        call check(error, is_close(alpha, norm(x%data, 2)))
        call check_test(error, 'test_vector_cdp_norm', eq='is_close(x%norm, norm(x, 2))')
        
        return
    end subroutine test_vector_cdp_norm

    subroutine test_vector_cdp_add(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error

        ! Test vectors.
        type(dense_vector_cdp), allocatable :: x, y, z

        ! Initialize vectors.
        x = dense_vector_cdp(n) ; call x%rand()
        y = dense_vector_cdp(n) ; call y%rand()
        z = x

        ! Vector addition.
        call z%add(y)

        ! Check correctness.
        call check(error, norm(z%data - x%data - y%data, 2) < rtol_dp)
        call check_test(error, 'test_vector_cdp_norm', eq='is_close(x%norm, norm(z - (x+y), 2))')

        return
    end subroutine test_vector_cdp_add
 
    subroutine test_vector_cdp_sub(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error

        ! Test vectors.
        type(dense_vector_cdp), allocatable :: x, y, z

        ! Initialize vectors.
        x = dense_vector_cdp(n) ; call x%rand()
        y = dense_vector_cdp(n) ; call y%rand()
        z = x

        ! Vector addition.
        call z%sub(y)

        ! Check correctness.
        call check(error, norm(z%data - (x%data - y%data), 2) < rtol_dp)
        call check_test(error, 'test_vector_cdp_norm', eq='is_close(x%norm, norm(z - (x-y), 2))')

        return
    end subroutine test_vector_cdp_sub

    subroutine test_vector_cdp_dot(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error

        ! Test vectors.
        type(dense_vector_cdp), allocatable :: x, y
        complex(dp) :: alpha

        ! Initialize vectors.
        x = dense_vector_cdp(n) ; call x%rand()
        y = dense_vector_cdp(n) ; call y%rand()

        ! Compute inner-product.
        alpha = x%dot(y)

        ! Check correctness.
        call check(error, abs(alpha - dot_product(x%data, y%data)) < rtol_dp)
        call check_test(error, 'test_vector_cdp_dot', eq='abs(x%dot(y) - dot_product(x, y))')

        return
    end subroutine test_vector_cdp_dot

    subroutine test_vector_cdp_scal(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error

        ! Test vector.
        type(dense_vector_cdp), allocatable :: x, y
        complex(dp) :: alpha

        ! Initialize vector.
        x = dense_vector_cdp(n) ; call x%rand(ifnorm=.true.)
        y = dense_vector_cdp(x%data)
        alpha = 0.0_dp ; call random_number(alpha%re) ; call random_number(alpha%im)
        alpha = alpha / abs(alpha)
        print *, "alpha : ", alpha
        print *
        print *, "x     : ", x%data
        print *
        print *, "y     : ", y%data
        print *
        print *, "x-y     : ", x%data - y%data
        
        ! Scale the vector.
        call x%scal(alpha)
        print *
        print *, "scal(x) :", x%data
        print *
        print *, "scal(y) :", alpha*y%data
        print *
        print *, "diff    :", abs(x%data - alpha*y%data)
        print *
        print *, "norm(diff) :", norm(x%data - alpha*y%data, 2) 

        ! Check correctness.
        call check(error, norm(x%data - alpha*y%data, 2) < rtol_dp)
        call check_test(error, 'test_vector_cdp_scal', eq='norm(x - alpha*y, 2)')

        return
    end subroutine test_vector_cdp_scal


end module TestVectors
