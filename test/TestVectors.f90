module TestVectors
    ! Fortran Standard Library
    use stdlib_math, only: is_close, all_close
    use stdlib_linalg, only: norm, eye
    use stdlib_stats_distribution_normal, only: normal => rvs_normal
    use stdlib_optval, only: optval
    ! Testdrive
    use testdrive, only: new_unittest, unittest_type, error_type, check
    ! LightKrylov
    use LightKrylov
    use LightKrylov_Logger
    ! Test Utilities
    use LightKrylov_Constants
    use LightKrylov_TestUtils
    use TestUtils
    implicit none (type, external)

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
                    new_unittest("Vector dot product", test_vector_rsp_dot), &
                    new_unittest("Vector space axioms", test_vector_axioms_rsp), &
                    new_unittest("Vector get_size", test_vector_rsp_get_size), &
                    new_unittest("Vector chsgn", test_vector_rsp_chsgn), &
                    new_unittest("Linear combination vector", test_linear_combination_vector_rsp), &
                    ! new_unittest("Linear combination matrix", test_linear_combination_matrix_rsp), &
                    new_unittest("Gram matrix", test_gram_rsp), &
                    new_unittest("Innerprod vector", test_innerprod_vector_rsp), &
                    new_unittest("Innerprod matrix", test_innerprod_matrix_rsp), &
                    new_unittest("axpby_basis", test_axpby_basis_rsp), &
                    new_unittest("zero_basis", test_zero_basis_rsp), &
                    new_unittest("copy_basis", test_copy_basis_rsp), &
                    new_unittest("rand_basis", test_rand_basis_rsp) &
                    ]
        return
    end subroutine collect_vector_rsp_testsuite

    subroutine test_vector_axioms_rsp(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_rsp) :: x
        real(sp) :: x_(n)
        logical :: success
        ! Initialize vector.
        x_ = 0.0_sp ; x = dense_vector(x_)
        success = verify_vector_axioms(x)
        call check(error, success .eqv. .true.)
        call check_test(error, 'test_vector_axioms_rsp', eq='Vector space axioms')
    end subroutine test_vector_axioms_rsp

    subroutine test_vector_rsp_norm(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vector.
        type(dense_vector_rsp) :: x
        real(sp) :: x_(n)
        real(sp) :: alpha

        ! Initialize vector.
        x_ = 0.0_sp ; x = dense_vector(x_) ; call x%rand()
        
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
        real(sp) :: x_(n), y_(n), z_(n)

        ! Initialize vectors.
        x = dense_vector(x_) ; call x%rand()
        y = dense_vector(y_) ; call y%rand()
        z = x

        ! Vector addition.
        call z%add(y)

        ! Check correctness.
        call check(error, norm(z%data - x%data - y%data, 2) < rtol_sp)
        call check_test(error, 'test_vector_rsp_add', eq='is_close(x%norm, norm(z - (x+y), 2))')

        return
    end subroutine test_vector_rsp_add
 
    subroutine test_vector_rsp_sub(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error

        ! Test vectors.
        type(dense_vector_rsp), allocatable :: x, y, z
        real(sp) :: x_(n), y_(n), z_(n)

        ! Initialize vectors.
        x = dense_vector(x_) ; call x%rand()
        y = dense_vector(y_) ; call y%rand()
        z = x

        ! Vector addition.
        call z%sub(y)

        ! Check correctness.
        call check(error, norm(z%data - (x%data - y%data), 2) < rtol_sp)
        call check_test(error, 'test_vector_rsp_sub', eq='is_close(x%norm, norm(z - (x-y), 2))')

        return
    end subroutine test_vector_rsp_sub

    subroutine test_vector_rsp_dot(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error

        ! Test vectors.
        type(dense_vector_rsp), allocatable :: x, y
        real(sp) :: x_(n), y_(n)
        real(sp) :: alpha

        ! Initialize vectors.
        x = dense_vector(x_) ; call x%rand()
        y = dense_vector(y_) ; call y%rand()

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
        real(sp) :: x_(n), y_(n)
        real(sp) :: alpha

        ! Initialize vector.
        x = dense_vector(x_) ; call x%rand(ifnorm=.true.)
        y = x
        call random_number(alpha)
        
        ! Scale the vector.
        call x%scal(alpha)

        ! Check correctness.
        call check(error, norm(x%data - alpha*y%data, 2) < rtol_sp)
        call check_test(error, 'test_vector_rsp_scal', eq='norm(x - alpha*y, 2)')

        return
    end subroutine test_vector_rsp_scal

    subroutine test_vector_rsp_get_size(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_rsp) :: x
        real(sp) :: x_(n)
        x_ = 0.0_sp ; x = dense_vector(x_)
        call check(error, x%get_size() == n)
        call check_test(error, 'test_vector_rsp_get_size', eq='x%get_size() == n')
    end subroutine test_vector_rsp_get_size

    subroutine test_vector_rsp_chsgn(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_rsp) :: x, y
        real(sp) :: x_(n)
        x = dense_vector(x_) ; call x%rand()
        y = x
        call x%chsgn()
        call y%scal(-one_rsp)
        call check(error, norm(x%data - y%data, 2) < rtol_sp)
        call check_test(error, 'test_vector_rsp_chsgn', eq='norm(x - (-y))')
    end subroutine test_vector_rsp_chsgn

    subroutine test_linear_combination_vector_rsp(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_rsp), allocatable :: X(:), y
        real(sp) :: x_(n), v(5)
        integer :: k, i
        k = 5
        allocate(X(k))
        x_ = zero_rsp
        do i = 1, k
            X(i) = dense_vector(x_)
            call X(i)%rand()
        end do
        call random_number(v)
        call linear_combination(y, X, v)
        block
            real(sp) :: expected(n)
            expected = zero_rsp
            do i = 1, k
                expected = expected + v(i) * X(i)%data
            end do
            call check(error, norm(y%data - expected, 2) < rtol_sp)
        end block
        call check_test(error, 'test_linear_combination_vector_rsp', eq='y = Xv')
    end subroutine test_linear_combination_vector_rsp

    ! subroutine test_linear_combination_matrix_rsp(error)
    !     type(error_type), allocatable, intent(out) :: error
    !     type(dense_vector_rsp), allocatable :: X(:), Y(:)
    !     real(sp) :: x_(n)
    !     real(sp), allocatable :: B(:,:)
    !     integer :: k, m, i, j
    !     k = 5 ; m = 3
    !     allocate(X(k), Y(m), B(k, m))
    !     x_ = zero_rsp
    !     do i = 1, k
    !         X(i) = dense_vector(x_)
    !         call X(i)%rand()
    !     end do
    !     #:if type[0] == "c"
    !     call random_number(B%re)
    !     call random_number(B%im)
    !     #:else
    !     call random_number(B)
    !     #:endif
    !     call linear_combination(Y, X, B)
    !     block
    !         real(sp) :: expected(n, m)
    !         expected = zero_rsp
    !         do j = 1, m
    !             do i = 1, k
    !                 expected(:, j) = expected(:, j) + B(i, j) * X(i)%data
    !             end do
    !         end do
    !         call check(error, norm(Y(1)%data - expected(:,1), 2) < rtol_sp .and. &
    !                        norm(Y(m)%data - expected(:,m), 2) < rtol_sp)
    !     end block
    !     call check_test(error, 'test_linear_combination_matrix_rsp', eq='Y = XB')
    ! end subroutine test_linear_combination_matrix_rsp

    subroutine test_gram_rsp(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_rsp), allocatable :: X(:)
        real(sp) :: x_(n)
        real(sp), allocatable :: G(:,:)
        integer :: k, i
        k = 5
        allocate(X(k))
        x_ = zero_rsp
        do i = 1, k
            X(i) = dense_vector(x_)
            call X(i)%rand()
        end do
        ! Orthonormalize X
        call orthonormalize_basis(X)
        G = Gram(X)
        call check(error, norm(G - eye(k, mold=1.0_sp), 2) < rtol_sp)
        call check_test(error, 'test_gram_rsp', eq='Gram(X) = I')
    end subroutine test_gram_rsp

    subroutine test_innerprod_vector_rsp(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_rsp), allocatable :: X(:), y
        real(sp) :: x_(n)
        real(sp), allocatable :: v(:)
        integer :: k, i
        k = 5
        allocate(X(k), y, v(k))
        x_ = zero_rsp
        do i = 1, k
            X(i) = dense_vector(x_)
            call X(i)%rand()
        end do
        y = dense_vector(x_)
        call y%rand()
        v = innerprod(X, y)
        block
            real(sp) :: expected(k)
            do i = 1, k
                expected(i) = X(i)%dot(y)
            end do
            call check(error, norm(v - expected, 2) < rtol_sp)
        end block
        call check_test(error, 'test_innerprod_vector_rsp', eq='v = X.dot(y)')
    end subroutine test_innerprod_vector_rsp

    subroutine test_innerprod_matrix_rsp(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_rsp), allocatable :: X(:), Y(:)
        real(sp) :: x_(n)
        real(sp), allocatable :: M(:,:)
        integer :: k, l, i, j
        k = 5 ; l = 3
        allocate(X(k), Y(l), M(k, l))
        x_ = zero_rsp
        do i = 1, k
            X(i) = dense_vector(x_)
            call X(i)%rand()
        end do
        do i = 1, l
            Y(i) = dense_vector(x_)
            call Y(i)%rand()
        end do
        M = innerprod(X, Y)
        block
            real(sp) :: expected(k, l)
            do j = 1, l
                do i = 1, k
                    expected(i, j) = X(i)%dot(Y(j))
                end do
            end do
            call check(error, norm(M - expected, 2) < rtol_sp)
        end block
        call check_test(error, 'test_innerprod_matrix_rsp', eq='M = X.dot(Y)')
    end subroutine test_innerprod_matrix_rsp

    subroutine test_axpby_basis_rsp(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_rsp), allocatable :: X(:), Y(:), Y_orig(:)
        real(sp) :: x_(n), alpha, beta
        integer :: k, i
        k = 5
        allocate(X(k), Y(k), Y_orig(k))
        x_ = zero_rsp
        do i = 1, k
            X(i) = dense_vector(x_)
            Y(i) = dense_vector(x_)
            Y_orig(i) = dense_vector(x_)
            call X(i)%rand()
            call Y(i)%rand()
            call copy(Y(i), Y_orig(i))
        end do
        call random_number(alpha)
        call random_number(beta)
        call axpby_basis(alpha, X, beta, Y)
        do i = 1, k
            call Y_orig(i)%axpby(alpha, X(i), beta)
        end do
        call check(error, norm(Y(1)%data - Y_orig(1)%data, 2) < rtol_sp)
        call check_test(error, 'test_axpby_basis_rsp', eq='Y = alpha*X + beta*Y')
    end subroutine test_axpby_basis_rsp

    subroutine test_zero_basis_rsp(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_rsp), allocatable :: X(:)
        real(sp) :: x_(n)
        integer :: k, i
        k = 5
        allocate(X(k))
        x_ = zero_rsp
        do i = 1, k
            X(i) = dense_vector(x_)
            call X(i)%rand()
        end do
        call zero_basis(X)
        call check(error, norm(X(1)%data, 2) <= atol_sp)
        call check_test(error, 'test_zero_basis_rsp', eq='X == 0')
    end subroutine test_zero_basis_rsp

    subroutine test_copy_basis_rsp(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_rsp), allocatable :: X(:), Y(:)
        real(sp) :: x_(n)
        integer :: k, i
        k = 5
        allocate(X(k), Y(k))
        x_ = zero_rsp
        do i = 1, k
            X(i) = dense_vector(x_)
            Y(i) = dense_vector(x_)
            call X(i)%rand()
        end do
        call copy(Y, X)
        call check(error, norm(X(1)%data - Y(1)%data, 2) < rtol_sp)
        call check_test(error, 'test_copy_basis_rsp', eq='Y == X')
    end subroutine test_copy_basis_rsp

    subroutine test_rand_basis_rsp(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_rsp), allocatable :: X(:)
        real(sp) :: x_(n)
        integer :: k, i
        real(sp) :: var1, var2
        k = 5
        allocate(X(k))
        x_ = zero_rsp
        do i = 1, k
            X(i) = dense_vector(x_)
        end do
        ! Test without normalization
        call rand_basis(X, ifnorm=.false.)
        var1 = sum(abs(X(1)%data))
        call X(1)%rand()
        var2 = sum(abs(X(2)%data))
        call check(error, abs(var1 - var2) > 0.0_sp .or. abs(var1) > 0.0_sp)
        call check_test(error, 'test_rand_basis_rsp', eq='rand vectors')
        ! Test with normalization
        call rand_basis(X, ifnorm=.true.)
        call check(error, abs(X(1)%norm() - 1.0_sp) < rtol_sp)
        call check_test(error, 'test_rand_basis_rsp', eq='norm == 1')
    end subroutine test_rand_basis_rsp

    subroutine collect_vector_rdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("Vector norm", test_vector_rdp_norm)      , &
                    new_unittest("Vector scale", test_vector_rdp_scal)     , &
                    new_unittest("Vector addition", test_vector_rdp_add)   , &
                    new_unittest("Vector subtraction", test_vector_rdp_sub), &
                    new_unittest("Vector dot product", test_vector_rdp_dot), &
                    new_unittest("Vector space axioms", test_vector_axioms_rdp), &
                    new_unittest("Vector get_size", test_vector_rdp_get_size), &
                    new_unittest("Vector chsgn", test_vector_rdp_chsgn), &
                    new_unittest("Linear combination vector", test_linear_combination_vector_rdp), &
                    ! new_unittest("Linear combination matrix", test_linear_combination_matrix_rdp), &
                    new_unittest("Gram matrix", test_gram_rdp), &
                    new_unittest("Innerprod vector", test_innerprod_vector_rdp), &
                    new_unittest("Innerprod matrix", test_innerprod_matrix_rdp), &
                    new_unittest("axpby_basis", test_axpby_basis_rdp), &
                    new_unittest("zero_basis", test_zero_basis_rdp), &
                    new_unittest("copy_basis", test_copy_basis_rdp), &
                    new_unittest("rand_basis", test_rand_basis_rdp) &
                    ]
        return
    end subroutine collect_vector_rdp_testsuite

    subroutine test_vector_axioms_rdp(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_rdp) :: x
        real(dp) :: x_(n)
        logical :: success
        ! Initialize vector.
        x_ = 0.0_dp ; x = dense_vector(x_)
        success = verify_vector_axioms(x)
        call check(error, success .eqv. .true.)
        call check_test(error, 'test_vector_axioms_rdp', eq='Vector space axioms')
    end subroutine test_vector_axioms_rdp

    subroutine test_vector_rdp_norm(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vector.
        type(dense_vector_rdp) :: x
        real(dp) :: x_(n)
        real(dp) :: alpha

        ! Initialize vector.
        x_ = 0.0_dp ; x = dense_vector(x_) ; call x%rand()
        
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
        real(dp) :: x_(n), y_(n), z_(n)

        ! Initialize vectors.
        x = dense_vector(x_) ; call x%rand()
        y = dense_vector(y_) ; call y%rand()
        z = x

        ! Vector addition.
        call z%add(y)

        ! Check correctness.
        call check(error, norm(z%data - x%data - y%data, 2) < rtol_dp)
        call check_test(error, 'test_vector_rdp_add', eq='is_close(x%norm, norm(z - (x+y), 2))')

        return
    end subroutine test_vector_rdp_add
 
    subroutine test_vector_rdp_sub(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error

        ! Test vectors.
        type(dense_vector_rdp), allocatable :: x, y, z
        real(dp) :: x_(n), y_(n), z_(n)

        ! Initialize vectors.
        x = dense_vector(x_) ; call x%rand()
        y = dense_vector(y_) ; call y%rand()
        z = x

        ! Vector addition.
        call z%sub(y)

        ! Check correctness.
        call check(error, norm(z%data - (x%data - y%data), 2) < rtol_dp)
        call check_test(error, 'test_vector_rdp_sub', eq='is_close(x%norm, norm(z - (x-y), 2))')

        return
    end subroutine test_vector_rdp_sub

    subroutine test_vector_rdp_dot(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error

        ! Test vectors.
        type(dense_vector_rdp), allocatable :: x, y
        real(dp) :: x_(n), y_(n)
        real(dp) :: alpha

        ! Initialize vectors.
        x = dense_vector(x_) ; call x%rand()
        y = dense_vector(y_) ; call y%rand()

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
        real(dp) :: x_(n), y_(n)
        real(dp) :: alpha

        ! Initialize vector.
        x = dense_vector(x_) ; call x%rand(ifnorm=.true.)
        y = x
        call random_number(alpha)
        
        ! Scale the vector.
        call x%scal(alpha)

        ! Check correctness.
        call check(error, norm(x%data - alpha*y%data, 2) < rtol_dp)
        call check_test(error, 'test_vector_rdp_scal', eq='norm(x - alpha*y, 2)')

        return
    end subroutine test_vector_rdp_scal

    subroutine test_vector_rdp_get_size(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_rdp) :: x
        real(dp) :: x_(n)
        x_ = 0.0_dp ; x = dense_vector(x_)
        call check(error, x%get_size() == n)
        call check_test(error, 'test_vector_rdp_get_size', eq='x%get_size() == n')
    end subroutine test_vector_rdp_get_size

    subroutine test_vector_rdp_chsgn(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_rdp) :: x, y
        real(dp) :: x_(n)
        x = dense_vector(x_) ; call x%rand()
        y = x
        call x%chsgn()
        call y%scal(-one_rdp)
        call check(error, norm(x%data - y%data, 2) < rtol_dp)
        call check_test(error, 'test_vector_rdp_chsgn', eq='norm(x - (-y))')
    end subroutine test_vector_rdp_chsgn

    subroutine test_linear_combination_vector_rdp(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_rdp), allocatable :: X(:), y
        real(dp) :: x_(n), v(5)
        integer :: k, i
        k = 5
        allocate(X(k))
        x_ = zero_rdp
        do i = 1, k
            X(i) = dense_vector(x_)
            call X(i)%rand()
        end do
        call random_number(v)
        call linear_combination(y, X, v)
        block
            real(dp) :: expected(n)
            expected = zero_rdp
            do i = 1, k
                expected = expected + v(i) * X(i)%data
            end do
            call check(error, norm(y%data - expected, 2) < rtol_dp)
        end block
        call check_test(error, 'test_linear_combination_vector_rdp', eq='y = Xv')
    end subroutine test_linear_combination_vector_rdp

    ! subroutine test_linear_combination_matrix_rdp(error)
    !     type(error_type), allocatable, intent(out) :: error
    !     type(dense_vector_rdp), allocatable :: X(:), Y(:)
    !     real(dp) :: x_(n)
    !     real(dp), allocatable :: B(:,:)
    !     integer :: k, m, i, j
    !     k = 5 ; m = 3
    !     allocate(X(k), Y(m), B(k, m))
    !     x_ = zero_rdp
    !     do i = 1, k
    !         X(i) = dense_vector(x_)
    !         call X(i)%rand()
    !     end do
    !     #:if type[0] == "c"
    !     call random_number(B%re)
    !     call random_number(B%im)
    !     #:else
    !     call random_number(B)
    !     #:endif
    !     call linear_combination(Y, X, B)
    !     block
    !         real(dp) :: expected(n, m)
    !         expected = zero_rdp
    !         do j = 1, m
    !             do i = 1, k
    !                 expected(:, j) = expected(:, j) + B(i, j) * X(i)%data
    !             end do
    !         end do
    !         call check(error, norm(Y(1)%data - expected(:,1), 2) < rtol_dp .and. &
    !                        norm(Y(m)%data - expected(:,m), 2) < rtol_dp)
    !     end block
    !     call check_test(error, 'test_linear_combination_matrix_rdp', eq='Y = XB')
    ! end subroutine test_linear_combination_matrix_rdp

    subroutine test_gram_rdp(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_rdp), allocatable :: X(:)
        real(dp) :: x_(n)
        real(dp), allocatable :: G(:,:)
        integer :: k, i
        k = 5
        allocate(X(k))
        x_ = zero_rdp
        do i = 1, k
            X(i) = dense_vector(x_)
            call X(i)%rand()
        end do
        ! Orthonormalize X
        call orthonormalize_basis(X)
        G = Gram(X)
        call check(error, norm(G - eye(k, mold=1.0_dp), 2) < rtol_dp)
        call check_test(error, 'test_gram_rdp', eq='Gram(X) = I')
    end subroutine test_gram_rdp

    subroutine test_innerprod_vector_rdp(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_rdp), allocatable :: X(:), y
        real(dp) :: x_(n)
        real(dp), allocatable :: v(:)
        integer :: k, i
        k = 5
        allocate(X(k), y, v(k))
        x_ = zero_rdp
        do i = 1, k
            X(i) = dense_vector(x_)
            call X(i)%rand()
        end do
        y = dense_vector(x_)
        call y%rand()
        v = innerprod(X, y)
        block
            real(dp) :: expected(k)
            do i = 1, k
                expected(i) = X(i)%dot(y)
            end do
            call check(error, norm(v - expected, 2) < rtol_dp)
        end block
        call check_test(error, 'test_innerprod_vector_rdp', eq='v = X.dot(y)')
    end subroutine test_innerprod_vector_rdp

    subroutine test_innerprod_matrix_rdp(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_rdp), allocatable :: X(:), Y(:)
        real(dp) :: x_(n)
        real(dp), allocatable :: M(:,:)
        integer :: k, l, i, j
        k = 5 ; l = 3
        allocate(X(k), Y(l), M(k, l))
        x_ = zero_rdp
        do i = 1, k
            X(i) = dense_vector(x_)
            call X(i)%rand()
        end do
        do i = 1, l
            Y(i) = dense_vector(x_)
            call Y(i)%rand()
        end do
        M = innerprod(X, Y)
        block
            real(dp) :: expected(k, l)
            do j = 1, l
                do i = 1, k
                    expected(i, j) = X(i)%dot(Y(j))
                end do
            end do
            call check(error, norm(M - expected, 2) < rtol_dp)
        end block
        call check_test(error, 'test_innerprod_matrix_rdp', eq='M = X.dot(Y)')
    end subroutine test_innerprod_matrix_rdp

    subroutine test_axpby_basis_rdp(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_rdp), allocatable :: X(:), Y(:), Y_orig(:)
        real(dp) :: x_(n), alpha, beta
        integer :: k, i
        k = 5
        allocate(X(k), Y(k), Y_orig(k))
        x_ = zero_rdp
        do i = 1, k
            X(i) = dense_vector(x_)
            Y(i) = dense_vector(x_)
            Y_orig(i) = dense_vector(x_)
            call X(i)%rand()
            call Y(i)%rand()
            call copy(Y(i), Y_orig(i))
        end do
        call random_number(alpha)
        call random_number(beta)
        call axpby_basis(alpha, X, beta, Y)
        do i = 1, k
            call Y_orig(i)%axpby(alpha, X(i), beta)
        end do
        call check(error, norm(Y(1)%data - Y_orig(1)%data, 2) < rtol_dp)
        call check_test(error, 'test_axpby_basis_rdp', eq='Y = alpha*X + beta*Y')
    end subroutine test_axpby_basis_rdp

    subroutine test_zero_basis_rdp(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_rdp), allocatable :: X(:)
        real(dp) :: x_(n)
        integer :: k, i
        k = 5
        allocate(X(k))
        x_ = zero_rdp
        do i = 1, k
            X(i) = dense_vector(x_)
            call X(i)%rand()
        end do
        call zero_basis(X)
        call check(error, norm(X(1)%data, 2) <= atol_dp)
        call check_test(error, 'test_zero_basis_rdp', eq='X == 0')
    end subroutine test_zero_basis_rdp

    subroutine test_copy_basis_rdp(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_rdp), allocatable :: X(:), Y(:)
        real(dp) :: x_(n)
        integer :: k, i
        k = 5
        allocate(X(k), Y(k))
        x_ = zero_rdp
        do i = 1, k
            X(i) = dense_vector(x_)
            Y(i) = dense_vector(x_)
            call X(i)%rand()
        end do
        call copy(Y, X)
        call check(error, norm(X(1)%data - Y(1)%data, 2) < rtol_dp)
        call check_test(error, 'test_copy_basis_rdp', eq='Y == X')
    end subroutine test_copy_basis_rdp

    subroutine test_rand_basis_rdp(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_rdp), allocatable :: X(:)
        real(dp) :: x_(n)
        integer :: k, i
        real(dp) :: var1, var2
        k = 5
        allocate(X(k))
        x_ = zero_rdp
        do i = 1, k
            X(i) = dense_vector(x_)
        end do
        ! Test without normalization
        call rand_basis(X, ifnorm=.false.)
        var1 = sum(abs(X(1)%data))
        call X(1)%rand()
        var2 = sum(abs(X(2)%data))
        call check(error, abs(var1 - var2) > 0.0_dp .or. abs(var1) > 0.0_dp)
        call check_test(error, 'test_rand_basis_rdp', eq='rand vectors')
        ! Test with normalization
        call rand_basis(X, ifnorm=.true.)
        call check(error, abs(X(1)%norm() - 1.0_dp) < rtol_dp)
        call check_test(error, 'test_rand_basis_rdp', eq='norm == 1')
    end subroutine test_rand_basis_rdp

    subroutine collect_vector_csp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("Vector norm", test_vector_csp_norm)      , &
                    new_unittest("Vector scale", test_vector_csp_scal)     , &
                    new_unittest("Vector addition", test_vector_csp_add)   , &
                    new_unittest("Vector subtraction", test_vector_csp_sub), &
                    new_unittest("Vector dot product", test_vector_csp_dot), &
                    new_unittest("Vector space axioms", test_vector_axioms_csp), &
                    new_unittest("Vector get_size", test_vector_csp_get_size), &
                    new_unittest("Vector chsgn", test_vector_csp_chsgn), &
                    new_unittest("Linear combination vector", test_linear_combination_vector_csp), &
                    ! new_unittest("Linear combination matrix", test_linear_combination_matrix_csp), &
                    new_unittest("Gram matrix", test_gram_csp), &
                    new_unittest("Innerprod vector", test_innerprod_vector_csp), &
                    new_unittest("Innerprod matrix", test_innerprod_matrix_csp), &
                    new_unittest("axpby_basis", test_axpby_basis_csp), &
                    new_unittest("zero_basis", test_zero_basis_csp), &
                    new_unittest("copy_basis", test_copy_basis_csp), &
                    new_unittest("rand_basis", test_rand_basis_csp) &
                    ]
        return
    end subroutine collect_vector_csp_testsuite

    subroutine test_vector_axioms_csp(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_csp) :: x
        complex(sp) :: x_(n)
        logical :: success
        ! Initialize vector.
        x_ = 0.0_sp ; x = dense_vector(x_)
        success = verify_vector_axioms(x)
        call check(error, success .eqv. .true.)
        call check_test(error, 'test_vector_axioms_csp', eq='Vector space axioms')
    end subroutine test_vector_axioms_csp

    subroutine test_vector_csp_norm(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vector.
        type(dense_vector_csp) :: x
        complex(sp) :: x_(n)
        real(sp) :: alpha

        ! Initialize vector.
        x_ = 0.0_sp ; x = dense_vector(x_) ; call x%rand()
        
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
        complex(sp) :: x_(n), y_(n), z_(n)

        ! Initialize vectors.
        x = dense_vector(x_) ; call x%rand()
        y = dense_vector(y_) ; call y%rand()
        z = x

        ! Vector addition.
        call z%add(y)

        ! Check correctness.
        call check(error, norm(z%data - x%data - y%data, 2) < rtol_sp)
        call check_test(error, 'test_vector_csp_add', eq='is_close(x%norm, norm(z - (x+y), 2))')

        return
    end subroutine test_vector_csp_add
 
    subroutine test_vector_csp_sub(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error

        ! Test vectors.
        type(dense_vector_csp), allocatable :: x, y, z
        complex(sp) :: x_(n), y_(n), z_(n)

        ! Initialize vectors.
        x = dense_vector(x_) ; call x%rand()
        y = dense_vector(y_) ; call y%rand()
        z = x

        ! Vector addition.
        call z%sub(y)

        ! Check correctness.
        call check(error, norm(z%data - (x%data - y%data), 2) < rtol_sp)
        call check_test(error, 'test_vector_csp_sub', eq='is_close(x%norm, norm(z - (x-y), 2))')

        return
    end subroutine test_vector_csp_sub

    subroutine test_vector_csp_dot(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error

        ! Test vectors.
        type(dense_vector_csp), allocatable :: x, y
        complex(sp) :: x_(n), y_(n)
        complex(sp) :: alpha

        ! Initialize vectors.
        x = dense_vector(x_) ; call x%rand()
        y = dense_vector(y_) ; call y%rand()

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
        complex(sp) :: x_(n), y_(n)
        complex(sp) :: alpha

        ! Initialize vector.
        x = dense_vector(x_) ; call x%rand(ifnorm=.true.)
        y = x
        alpha = 0.0_sp ; call random_number(alpha%re) ; call random_number(alpha%im)
        alpha = alpha / abs(alpha)
        
        ! Scale the vector.
        call x%scal(alpha)

        ! Check correctness.
        call check(error, norm(x%data - alpha*y%data, 2) < rtol_sp)
        call check_test(error, 'test_vector_csp_scal', eq='norm(x - alpha*y, 2)')

        return
    end subroutine test_vector_csp_scal

    subroutine test_vector_csp_get_size(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_csp) :: x
        complex(sp) :: x_(n)
        x_ = 0.0_sp ; x = dense_vector(x_)
        call check(error, x%get_size() == n)
        call check_test(error, 'test_vector_csp_get_size', eq='x%get_size() == n')
    end subroutine test_vector_csp_get_size

    subroutine test_vector_csp_chsgn(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_csp) :: x, y
        complex(sp) :: x_(n)
        x = dense_vector(x_) ; call x%rand()
        y = x
        call x%chsgn()
        call y%scal(-one_csp)
        call check(error, norm(x%data - y%data, 2) < rtol_sp)
        call check_test(error, 'test_vector_csp_chsgn', eq='norm(x - (-y))')
    end subroutine test_vector_csp_chsgn

    subroutine test_linear_combination_vector_csp(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_csp), allocatable :: X(:), y
        complex(sp) :: x_(n), v(5)
        integer :: k, i
        k = 5
        allocate(X(k))
        x_ = zero_csp
        do i = 1, k
            X(i) = dense_vector(x_)
            call X(i)%rand()
        end do
        call random_number(v%re)
        call random_number(v%im)
        call linear_combination(y, X, v)
        block
            complex(sp) :: expected(n)
            expected = zero_csp
            do i = 1, k
                expected = expected + v(i) * X(i)%data
            end do
            call check(error, norm(y%data - expected, 2) < rtol_sp)
        end block
        call check_test(error, 'test_linear_combination_vector_csp', eq='y = Xv')
    end subroutine test_linear_combination_vector_csp

    ! subroutine test_linear_combination_matrix_csp(error)
    !     type(error_type), allocatable, intent(out) :: error
    !     type(dense_vector_csp), allocatable :: X(:), Y(:)
    !     complex(sp) :: x_(n)
    !     complex(sp), allocatable :: B(:,:)
    !     integer :: k, m, i, j
    !     k = 5 ; m = 3
    !     allocate(X(k), Y(m), B(k, m))
    !     x_ = zero_csp
    !     do i = 1, k
    !         X(i) = dense_vector(x_)
    !         call X(i)%rand()
    !     end do
    !     #:if type[0] == "c"
    !     call random_number(B%re)
    !     call random_number(B%im)
    !     #:else
    !     call random_number(B)
    !     #:endif
    !     call linear_combination(Y, X, B)
    !     block
    !         complex(sp) :: expected(n, m)
    !         expected = zero_csp
    !         do j = 1, m
    !             do i = 1, k
    !                 expected(:, j) = expected(:, j) + B(i, j) * X(i)%data
    !             end do
    !         end do
    !         call check(error, norm(Y(1)%data - expected(:,1), 2) < rtol_sp .and. &
    !                        norm(Y(m)%data - expected(:,m), 2) < rtol_sp)
    !     end block
    !     call check_test(error, 'test_linear_combination_matrix_csp', eq='Y = XB')
    ! end subroutine test_linear_combination_matrix_csp

    subroutine test_gram_csp(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_csp), allocatable :: X(:)
        complex(sp) :: x_(n)
        complex(sp), allocatable :: G(:,:)
        integer :: k, i
        k = 5
        allocate(X(k))
        x_ = zero_csp
        do i = 1, k
            X(i) = dense_vector(x_)
            call X(i)%rand()
        end do
        ! Orthonormalize X
        call orthonormalize_basis(X)
        G = Gram(X)
        call check(error, norm(G - eye(k, mold=1.0_sp), 2) < rtol_sp)
        call check_test(error, 'test_gram_csp', eq='Gram(X) = I')
    end subroutine test_gram_csp

    subroutine test_innerprod_vector_csp(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_csp), allocatable :: X(:), y
        complex(sp) :: x_(n)
        complex(sp), allocatable :: v(:)
        integer :: k, i
        k = 5
        allocate(X(k), y, v(k))
        x_ = zero_csp
        do i = 1, k
            X(i) = dense_vector(x_)
            call X(i)%rand()
        end do
        y = dense_vector(x_)
        call y%rand()
        v = innerprod(X, y)
        block
            complex(sp) :: expected(k)
            do i = 1, k
                expected(i) = X(i)%dot(y)
            end do
            call check(error, norm(v - expected, 2) < rtol_sp)
        end block
        call check_test(error, 'test_innerprod_vector_csp', eq='v = X.dot(y)')
    end subroutine test_innerprod_vector_csp

    subroutine test_innerprod_matrix_csp(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_csp), allocatable :: X(:), Y(:)
        complex(sp) :: x_(n)
        complex(sp), allocatable :: M(:,:)
        integer :: k, l, i, j
        k = 5 ; l = 3
        allocate(X(k), Y(l), M(k, l))
        x_ = zero_csp
        do i = 1, k
            X(i) = dense_vector(x_)
            call X(i)%rand()
        end do
        do i = 1, l
            Y(i) = dense_vector(x_)
            call Y(i)%rand()
        end do
        M = innerprod(X, Y)
        block
            complex(sp) :: expected(k, l)
            do j = 1, l
                do i = 1, k
                    expected(i, j) = X(i)%dot(Y(j))
                end do
            end do
            call check(error, norm(M - expected, 2) < rtol_sp)
        end block
        call check_test(error, 'test_innerprod_matrix_csp', eq='M = X.dot(Y)')
    end subroutine test_innerprod_matrix_csp

    subroutine test_axpby_basis_csp(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_csp), allocatable :: X(:), Y(:), Y_orig(:)
        complex(sp) :: x_(n), alpha, beta
        integer :: k, i
        k = 5
        allocate(X(k), Y(k), Y_orig(k))
        x_ = zero_csp
        do i = 1, k
            X(i) = dense_vector(x_)
            Y(i) = dense_vector(x_)
            Y_orig(i) = dense_vector(x_)
            call X(i)%rand()
            call Y(i)%rand()
            call copy(Y(i), Y_orig(i))
        end do
        call random_number(alpha%re)
        call random_number(alpha%im)
        call random_number(beta%re)
        call random_number(beta%im)
        call axpby_basis(alpha, X, beta, Y)
        do i = 1, k
            call Y_orig(i)%axpby(alpha, X(i), beta)
        end do
        call check(error, norm(Y(1)%data - Y_orig(1)%data, 2) < rtol_sp)
        call check_test(error, 'test_axpby_basis_csp', eq='Y = alpha*X + beta*Y')
    end subroutine test_axpby_basis_csp

    subroutine test_zero_basis_csp(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_csp), allocatable :: X(:)
        complex(sp) :: x_(n)
        integer :: k, i
        k = 5
        allocate(X(k))
        x_ = zero_csp
        do i = 1, k
            X(i) = dense_vector(x_)
            call X(i)%rand()
        end do
        call zero_basis(X)
        call check(error, norm(X(1)%data, 2) <= atol_sp)
        call check_test(error, 'test_zero_basis_csp', eq='X == 0')
    end subroutine test_zero_basis_csp

    subroutine test_copy_basis_csp(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_csp), allocatable :: X(:), Y(:)
        complex(sp) :: x_(n)
        integer :: k, i
        k = 5
        allocate(X(k), Y(k))
        x_ = zero_csp
        do i = 1, k
            X(i) = dense_vector(x_)
            Y(i) = dense_vector(x_)
            call X(i)%rand()
        end do
        call copy(Y, X)
        call check(error, norm(X(1)%data - Y(1)%data, 2) < rtol_sp)
        call check_test(error, 'test_copy_basis_csp', eq='Y == X')
    end subroutine test_copy_basis_csp

    subroutine test_rand_basis_csp(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_csp), allocatable :: X(:)
        complex(sp) :: x_(n)
        integer :: k, i
        real(sp) :: var1, var2
        k = 5
        allocate(X(k))
        x_ = zero_csp
        do i = 1, k
            X(i) = dense_vector(x_)
        end do
        ! Test without normalization
        call rand_basis(X, ifnorm=.false.)
        var1 = sum(abs(X(1)%data))
        call X(1)%rand()
        var2 = sum(abs(X(2)%data))
        call check(error, abs(var1 - var2) > 0.0_sp .or. abs(var1) > 0.0_sp)
        call check_test(error, 'test_rand_basis_csp', eq='rand vectors')
        ! Test with normalization
        call rand_basis(X, ifnorm=.true.)
        call check(error, abs(X(1)%norm() - 1.0_sp) < rtol_sp)
        call check_test(error, 'test_rand_basis_csp', eq='norm == 1')
    end subroutine test_rand_basis_csp

    subroutine collect_vector_cdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("Vector norm", test_vector_cdp_norm)      , &
                    new_unittest("Vector scale", test_vector_cdp_scal)     , &
                    new_unittest("Vector addition", test_vector_cdp_add)   , &
                    new_unittest("Vector subtraction", test_vector_cdp_sub), &
                    new_unittest("Vector dot product", test_vector_cdp_dot), &
                    new_unittest("Vector space axioms", test_vector_axioms_cdp), &
                    new_unittest("Vector get_size", test_vector_cdp_get_size), &
                    new_unittest("Vector chsgn", test_vector_cdp_chsgn), &
                    new_unittest("Linear combination vector", test_linear_combination_vector_cdp), &
                    ! new_unittest("Linear combination matrix", test_linear_combination_matrix_cdp), &
                    new_unittest("Gram matrix", test_gram_cdp), &
                    new_unittest("Innerprod vector", test_innerprod_vector_cdp), &
                    new_unittest("Innerprod matrix", test_innerprod_matrix_cdp), &
                    new_unittest("axpby_basis", test_axpby_basis_cdp), &
                    new_unittest("zero_basis", test_zero_basis_cdp), &
                    new_unittest("copy_basis", test_copy_basis_cdp), &
                    new_unittest("rand_basis", test_rand_basis_cdp) &
                    ]
        return
    end subroutine collect_vector_cdp_testsuite

    subroutine test_vector_axioms_cdp(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_cdp) :: x
        complex(dp) :: x_(n)
        logical :: success
        ! Initialize vector.
        x_ = 0.0_dp ; x = dense_vector(x_)
        success = verify_vector_axioms(x)
        call check(error, success .eqv. .true.)
        call check_test(error, 'test_vector_axioms_cdp', eq='Vector space axioms')
    end subroutine test_vector_axioms_cdp

    subroutine test_vector_cdp_norm(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test vector.
        type(dense_vector_cdp) :: x
        complex(dp) :: x_(n)
        real(dp) :: alpha

        ! Initialize vector.
        x_ = 0.0_dp ; x = dense_vector(x_) ; call x%rand()
        
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
        complex(dp) :: x_(n), y_(n), z_(n)

        ! Initialize vectors.
        x = dense_vector(x_) ; call x%rand()
        y = dense_vector(y_) ; call y%rand()
        z = x

        ! Vector addition.
        call z%add(y)

        ! Check correctness.
        call check(error, norm(z%data - x%data - y%data, 2) < rtol_dp)
        call check_test(error, 'test_vector_cdp_add', eq='is_close(x%norm, norm(z - (x+y), 2))')

        return
    end subroutine test_vector_cdp_add
 
    subroutine test_vector_cdp_sub(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error

        ! Test vectors.
        type(dense_vector_cdp), allocatable :: x, y, z
        complex(dp) :: x_(n), y_(n), z_(n)

        ! Initialize vectors.
        x = dense_vector(x_) ; call x%rand()
        y = dense_vector(y_) ; call y%rand()
        z = x

        ! Vector addition.
        call z%sub(y)

        ! Check correctness.
        call check(error, norm(z%data - (x%data - y%data), 2) < rtol_dp)
        call check_test(error, 'test_vector_cdp_sub', eq='is_close(x%norm, norm(z - (x-y), 2))')

        return
    end subroutine test_vector_cdp_sub

    subroutine test_vector_cdp_dot(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error

        ! Test vectors.
        type(dense_vector_cdp), allocatable :: x, y
        complex(dp) :: x_(n), y_(n)
        complex(dp) :: alpha

        ! Initialize vectors.
        x = dense_vector(x_) ; call x%rand()
        y = dense_vector(y_) ; call y%rand()

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
        complex(dp) :: x_(n), y_(n)
        complex(dp) :: alpha

        ! Initialize vector.
        x = dense_vector(x_) ; call x%rand(ifnorm=.true.)
        y = x
        alpha = 0.0_dp ; call random_number(alpha%re) ; call random_number(alpha%im)
        alpha = alpha / abs(alpha)
        
        ! Scale the vector.
        call x%scal(alpha)

        ! Check correctness.
        call check(error, norm(x%data - alpha*y%data, 2) < rtol_dp)
        call check_test(error, 'test_vector_cdp_scal', eq='norm(x - alpha*y, 2)')

        return
    end subroutine test_vector_cdp_scal

    subroutine test_vector_cdp_get_size(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_cdp) :: x
        complex(dp) :: x_(n)
        x_ = 0.0_dp ; x = dense_vector(x_)
        call check(error, x%get_size() == n)
        call check_test(error, 'test_vector_cdp_get_size', eq='x%get_size() == n')
    end subroutine test_vector_cdp_get_size

    subroutine test_vector_cdp_chsgn(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_cdp) :: x, y
        complex(dp) :: x_(n)
        x = dense_vector(x_) ; call x%rand()
        y = x
        call x%chsgn()
        call y%scal(-one_cdp)
        call check(error, norm(x%data - y%data, 2) < rtol_dp)
        call check_test(error, 'test_vector_cdp_chsgn', eq='norm(x - (-y))')
    end subroutine test_vector_cdp_chsgn

    subroutine test_linear_combination_vector_cdp(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_cdp), allocatable :: X(:), y
        complex(dp) :: x_(n), v(5)
        integer :: k, i
        k = 5
        allocate(X(k))
        x_ = zero_cdp
        do i = 1, k
            X(i) = dense_vector(x_)
            call X(i)%rand()
        end do
        call random_number(v%re)
        call random_number(v%im)
        call linear_combination(y, X, v)
        block
            complex(dp) :: expected(n)
            expected = zero_cdp
            do i = 1, k
                expected = expected + v(i) * X(i)%data
            end do
            call check(error, norm(y%data - expected, 2) < rtol_dp)
        end block
        call check_test(error, 'test_linear_combination_vector_cdp', eq='y = Xv')
    end subroutine test_linear_combination_vector_cdp

    ! subroutine test_linear_combination_matrix_cdp(error)
    !     type(error_type), allocatable, intent(out) :: error
    !     type(dense_vector_cdp), allocatable :: X(:), Y(:)
    !     complex(dp) :: x_(n)
    !     complex(dp), allocatable :: B(:,:)
    !     integer :: k, m, i, j
    !     k = 5 ; m = 3
    !     allocate(X(k), Y(m), B(k, m))
    !     x_ = zero_cdp
    !     do i = 1, k
    !         X(i) = dense_vector(x_)
    !         call X(i)%rand()
    !     end do
    !     #:if type[0] == "c"
    !     call random_number(B%re)
    !     call random_number(B%im)
    !     #:else
    !     call random_number(B)
    !     #:endif
    !     call linear_combination(Y, X, B)
    !     block
    !         complex(dp) :: expected(n, m)
    !         expected = zero_cdp
    !         do j = 1, m
    !             do i = 1, k
    !                 expected(:, j) = expected(:, j) + B(i, j) * X(i)%data
    !             end do
    !         end do
    !         call check(error, norm(Y(1)%data - expected(:,1), 2) < rtol_dp .and. &
    !                        norm(Y(m)%data - expected(:,m), 2) < rtol_dp)
    !     end block
    !     call check_test(error, 'test_linear_combination_matrix_cdp', eq='Y = XB')
    ! end subroutine test_linear_combination_matrix_cdp

    subroutine test_gram_cdp(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_cdp), allocatable :: X(:)
        complex(dp) :: x_(n)
        complex(dp), allocatable :: G(:,:)
        integer :: k, i
        k = 5
        allocate(X(k))
        x_ = zero_cdp
        do i = 1, k
            X(i) = dense_vector(x_)
            call X(i)%rand()
        end do
        ! Orthonormalize X
        call orthonormalize_basis(X)
        G = Gram(X)
        call check(error, norm(G - eye(k, mold=1.0_dp), 2) < rtol_dp)
        call check_test(error, 'test_gram_cdp', eq='Gram(X) = I')
    end subroutine test_gram_cdp

    subroutine test_innerprod_vector_cdp(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_cdp), allocatable :: X(:), y
        complex(dp) :: x_(n)
        complex(dp), allocatable :: v(:)
        integer :: k, i
        k = 5
        allocate(X(k), y, v(k))
        x_ = zero_cdp
        do i = 1, k
            X(i) = dense_vector(x_)
            call X(i)%rand()
        end do
        y = dense_vector(x_)
        call y%rand()
        v = innerprod(X, y)
        block
            complex(dp) :: expected(k)
            do i = 1, k
                expected(i) = X(i)%dot(y)
            end do
            call check(error, norm(v - expected, 2) < rtol_dp)
        end block
        call check_test(error, 'test_innerprod_vector_cdp', eq='v = X.dot(y)')
    end subroutine test_innerprod_vector_cdp

    subroutine test_innerprod_matrix_cdp(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_cdp), allocatable :: X(:), Y(:)
        complex(dp) :: x_(n)
        complex(dp), allocatable :: M(:,:)
        integer :: k, l, i, j
        k = 5 ; l = 3
        allocate(X(k), Y(l), M(k, l))
        x_ = zero_cdp
        do i = 1, k
            X(i) = dense_vector(x_)
            call X(i)%rand()
        end do
        do i = 1, l
            Y(i) = dense_vector(x_)
            call Y(i)%rand()
        end do
        M = innerprod(X, Y)
        block
            complex(dp) :: expected(k, l)
            do j = 1, l
                do i = 1, k
                    expected(i, j) = X(i)%dot(Y(j))
                end do
            end do
            call check(error, norm(M - expected, 2) < rtol_dp)
        end block
        call check_test(error, 'test_innerprod_matrix_cdp', eq='M = X.dot(Y)')
    end subroutine test_innerprod_matrix_cdp

    subroutine test_axpby_basis_cdp(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_cdp), allocatable :: X(:), Y(:), Y_orig(:)
        complex(dp) :: x_(n), alpha, beta
        integer :: k, i
        k = 5
        allocate(X(k), Y(k), Y_orig(k))
        x_ = zero_cdp
        do i = 1, k
            X(i) = dense_vector(x_)
            Y(i) = dense_vector(x_)
            Y_orig(i) = dense_vector(x_)
            call X(i)%rand()
            call Y(i)%rand()
            call copy(Y(i), Y_orig(i))
        end do
        call random_number(alpha%re)
        call random_number(alpha%im)
        call random_number(beta%re)
        call random_number(beta%im)
        call axpby_basis(alpha, X, beta, Y)
        do i = 1, k
            call Y_orig(i)%axpby(alpha, X(i), beta)
        end do
        call check(error, norm(Y(1)%data - Y_orig(1)%data, 2) < rtol_dp)
        call check_test(error, 'test_axpby_basis_cdp', eq='Y = alpha*X + beta*Y')
    end subroutine test_axpby_basis_cdp

    subroutine test_zero_basis_cdp(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_cdp), allocatable :: X(:)
        complex(dp) :: x_(n)
        integer :: k, i
        k = 5
        allocate(X(k))
        x_ = zero_cdp
        do i = 1, k
            X(i) = dense_vector(x_)
            call X(i)%rand()
        end do
        call zero_basis(X)
        call check(error, norm(X(1)%data, 2) <= atol_dp)
        call check_test(error, 'test_zero_basis_cdp', eq='X == 0')
    end subroutine test_zero_basis_cdp

    subroutine test_copy_basis_cdp(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_cdp), allocatable :: X(:), Y(:)
        complex(dp) :: x_(n)
        integer :: k, i
        k = 5
        allocate(X(k), Y(k))
        x_ = zero_cdp
        do i = 1, k
            X(i) = dense_vector(x_)
            Y(i) = dense_vector(x_)
            call X(i)%rand()
        end do
        call copy(Y, X)
        call check(error, norm(X(1)%data - Y(1)%data, 2) < rtol_dp)
        call check_test(error, 'test_copy_basis_cdp', eq='Y == X')
    end subroutine test_copy_basis_cdp

    subroutine test_rand_basis_cdp(error)
        type(error_type), allocatable, intent(out) :: error
        type(dense_vector_cdp), allocatable :: X(:)
        complex(dp) :: x_(n)
        integer :: k, i
        real(dp) :: var1, var2
        k = 5
        allocate(X(k))
        x_ = zero_cdp
        do i = 1, k
            X(i) = dense_vector(x_)
        end do
        ! Test without normalization
        call rand_basis(X, ifnorm=.false.)
        var1 = sum(abs(X(1)%data))
        call X(1)%rand()
        var2 = sum(abs(X(2)%data))
        call check(error, abs(var1 - var2) > 0.0_dp .or. abs(var1) > 0.0_dp)
        call check_test(error, 'test_rand_basis_cdp', eq='rand vectors')
        ! Test with normalization
        call rand_basis(X, ifnorm=.true.)
        call check(error, abs(X(1)%norm() - 1.0_dp) < rtol_dp)
        call check_test(error, 'test_rand_basis_cdp', eq='norm == 1')
    end subroutine test_rand_basis_cdp


end module TestVectors
