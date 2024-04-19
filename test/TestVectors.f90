module TestVectors
    use LightKrylov
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use stdlib_math, only: is_close, all_close
    use stdlib_optval, only: optval
    implicit none
    private

    integer, parameter, public :: test_size = 128

    public :: collect_vector_rsp_testsuite
    public :: collect_vector_rdp_testsuite
    public :: collect_vector_csp_testsuite
    public :: collect_vector_cdp_testsuite

    type, extends(abstract_vector_rsp), public :: vector_rsp
        real(sp), dimension(test_size) :: data = 0.0_sp
    contains
        private
        procedure, pass(self), public :: zero => zero_rsp
        procedure, pass(self), public :: dot => dot_rsp
        procedure, pass(self), public :: scal => scal_rsp
        procedure, pass(self), public :: axpby => axpby_rsp
        procedure, pass(self), public :: rand => rand_rsp
    end type vector_rsp

    type, extends(abstract_vector_rdp), public :: vector_rdp
        real(dp), dimension(test_size) :: data = 0.0_dp
    contains
        private
        procedure, pass(self), public :: zero => zero_rdp
        procedure, pass(self), public :: dot => dot_rdp
        procedure, pass(self), public :: scal => scal_rdp
        procedure, pass(self), public :: axpby => axpby_rdp
        procedure, pass(self), public :: rand => rand_rdp
    end type vector_rdp

    type, extends(abstract_vector_csp), public :: vector_csp
        complex(sp), dimension(test_size) :: data = 0.0_sp
    contains
        private
        procedure, pass(self), public :: zero => zero_csp
        procedure, pass(self), public :: dot => dot_csp
        procedure, pass(self), public :: scal => scal_csp
        procedure, pass(self), public :: axpby => axpby_csp
        procedure, pass(self), public :: rand => rand_csp
    end type vector_csp

    type, extends(abstract_vector_cdp), public :: vector_cdp
        complex(dp), dimension(test_size) :: data = 0.0_dp
    contains
        private
        procedure, pass(self), public :: zero => zero_cdp
        procedure, pass(self), public :: dot => dot_cdp
        procedure, pass(self), public :: scal => scal_cdp
        procedure, pass(self), public :: axpby => axpby_cdp
        procedure, pass(self), public :: rand => rand_cdp
    end type vector_cdp

contains
    
    !--------------------------------------------------------------------
    !-----     DEFINITIONS OF THE VARIOUS TYPE-BOUND PROCEDURES     -----
    !--------------------------------------------------------------------

    subroutine zero_rsp(self)
        class(vector_rsp), intent(inout) :: self
        self%data = 0.0_sp
        return
    end subroutine zero_rsp

    function dot_rsp(self, vec) result(alpha)
        class(vector_rsp), intent(in) :: self
        class(abstract_vector_rsp), intent(in) :: vec
        real(sp) :: alpha

        select type(vec)
        type is(vector_rsp)
            alpha = dot_product(self%data, vec%data)
        end select
    end function dot_rsp

    subroutine scal_rsp(self, alpha)
        class(vector_rsp), intent(inout) :: self
        real(sp), intent(in) :: alpha
        self%data = alpha * self%data
        return
    end subroutine scal_rsp

    subroutine axpby_rsp(self, alpha, vec, beta)
        class(vector_rsp), intent(inout) :: self
        class(abstract_vector_rsp), intent(in) :: vec
        real(sp), intent(in) :: alpha, beta

        select type(vec)
        type is(vector_rsp)
            self%data = alpha*self%data + beta*vec%data
        end select
        return
    end subroutine

    subroutine rand_rsp(self, ifnorm)
        class(vector_rsp), intent(inout) :: self
        logical, optional, intent(in) :: ifnorm
        logical :: normalized
        real(sp) :: alpha
        call random_number(self%data)

        normalized = optval(ifnorm, .false.)
        if (normalized) then
            alpha = self%norm()
            call self%scal(1.0_sp/alpha)
        endif
    end subroutine rand_rsp

    subroutine zero_rdp(self)
        class(vector_rdp), intent(inout) :: self
        self%data = 0.0_dp
        return
    end subroutine zero_rdp

    function dot_rdp(self, vec) result(alpha)
        class(vector_rdp), intent(in) :: self
        class(abstract_vector_rdp), intent(in) :: vec
        real(dp) :: alpha

        select type(vec)
        type is(vector_rdp)
            alpha = dot_product(self%data, vec%data)
        end select
    end function dot_rdp

    subroutine scal_rdp(self, alpha)
        class(vector_rdp), intent(inout) :: self
        real(dp), intent(in) :: alpha
        self%data = alpha * self%data
        return
    end subroutine scal_rdp

    subroutine axpby_rdp(self, alpha, vec, beta)
        class(vector_rdp), intent(inout) :: self
        class(abstract_vector_rdp), intent(in) :: vec
        real(dp), intent(in) :: alpha, beta

        select type(vec)
        type is(vector_rdp)
            self%data = alpha*self%data + beta*vec%data
        end select
        return
    end subroutine

    subroutine rand_rdp(self, ifnorm)
        class(vector_rdp), intent(inout) :: self
        logical, optional, intent(in) :: ifnorm
        logical :: normalized
        real(dp) :: alpha
        call random_number(self%data)

        normalized = optval(ifnorm, .false.)
        if (normalized) then
            alpha = self%norm()
            call self%scal(1.0_dp/alpha)
        endif
    end subroutine rand_rdp

    subroutine zero_csp(self)
        class(vector_csp), intent(inout) :: self
        self%data = 0.0_sp
        return
    end subroutine zero_csp

    function dot_csp(self, vec) result(alpha)
        class(vector_csp), intent(in) :: self
        class(abstract_vector_csp), intent(in) :: vec
        complex(sp) :: alpha

        select type(vec)
        type is(vector_csp)
            alpha = dot_product(self%data, vec%data)
        end select
    end function dot_csp

    subroutine scal_csp(self, alpha)
        class(vector_csp), intent(inout) :: self
        complex(sp), intent(in) :: alpha
        self%data = alpha * self%data
        return
    end subroutine scal_csp

    subroutine axpby_csp(self, alpha, vec, beta)
        class(vector_csp), intent(inout) :: self
        class(abstract_vector_csp), intent(in) :: vec
        complex(sp), intent(in) :: alpha, beta

        select type(vec)
        type is(vector_csp)
            self%data = alpha*self%data + beta*vec%data
        end select
        return
    end subroutine

    subroutine rand_csp(self, ifnorm)
        class(vector_csp), intent(inout) :: self
        logical, optional, intent(in) :: ifnorm
        logical :: normalized
        complex(sp) :: alpha
        real(sp), dimension(test_size, 2) :: data

        call random_number(data)
        self%data%re = data(:, 1)
        self%data%im = data(:, 2)

        normalized = optval(ifnorm, .false.)
        if (normalized) then
            alpha = self%norm()
            call self%scal(1.0_sp/alpha)
        endif
    end subroutine rand_csp

    subroutine zero_cdp(self)
        class(vector_cdp), intent(inout) :: self
        self%data = 0.0_dp
        return
    end subroutine zero_cdp

    function dot_cdp(self, vec) result(alpha)
        class(vector_cdp), intent(in) :: self
        class(abstract_vector_cdp), intent(in) :: vec
        complex(dp) :: alpha

        select type(vec)
        type is(vector_cdp)
            alpha = dot_product(self%data, vec%data)
        end select
    end function dot_cdp

    subroutine scal_cdp(self, alpha)
        class(vector_cdp), intent(inout) :: self
        complex(dp), intent(in) :: alpha
        self%data = alpha * self%data
        return
    end subroutine scal_cdp

    subroutine axpby_cdp(self, alpha, vec, beta)
        class(vector_cdp), intent(inout) :: self
        class(abstract_vector_cdp), intent(in) :: vec
        complex(dp), intent(in) :: alpha, beta

        select type(vec)
        type is(vector_cdp)
            self%data = alpha*self%data + beta*vec%data
        end select
        return
    end subroutine

    subroutine rand_cdp(self, ifnorm)
        class(vector_cdp), intent(inout) :: self
        logical, optional, intent(in) :: ifnorm
        logical :: normalized
        complex(dp) :: alpha
        real(dp), dimension(test_size, 2) :: data

        call random_number(data)
        self%data%re = data(:, 1)
        self%data%im = data(:, 2)

        normalized = optval(ifnorm, .false.)
        if (normalized) then
            alpha = self%norm()
            call self%scal(1.0_dp/alpha)
        endif
    end subroutine rand_cdp

    
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
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        !> Test vector.
        type(vector_rsp), allocatable :: x
        real(sp) :: alpha

        ! Initialize vector.
        x = vector_rsp() ; call x%rand()
        
        ! Compute its norm.
        alpha = x%norm()

        ! Check result.
        call check(error, is_close(alpha, norm2(x%data)))
        
        return
    end subroutine test_vector_rsp_norm

    subroutine test_vector_rsp_add(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error

        !> Test vectors.
        type(vector_rsp), allocatable :: x, y, z

        ! Initialize vectors.
        x = vector_rsp() ; call x%rand()
        y = vector_rsp() ; call y%rand()
        z = x

        ! Vector addition.
        call z%add(y)

        ! Check correctness.
        call check(error, norm2(z%data - x%data - y%data) < rtol_sp)

        return
    end subroutine test_vector_rsp_add
 
    subroutine test_vector_rsp_sub(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error

        !> Test vectors.
        type(vector_rsp), allocatable :: x, y, z

        ! Initialize vectors.
        x = vector_rsp() ; call x%rand()
        y = vector_rsp() ; call y%rand()
        z = x

        ! Vector addition.
        call z%sub(y)

        ! Check correctness.
        call check(error, norm2(z%data - x%data + y%data) < rtol_sp)

        return
    end subroutine test_vector_rsp_sub

    subroutine test_vector_rsp_dot(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error

        !> Test vectors.
        type(vector_rsp), allocatable :: x, y
        real(sp) :: alpha

        ! Initialize vectors.
        x = vector_rsp() ; call x%rand()
        y = vector_rsp() ; call y%rand()

        ! Compute inner-product.
        alpha = x%dot(y)

        ! Check correctness.
        call check(error, abs(alpha - dot_product(x%data, y%data)) < rtol_sp)

        return
    end subroutine test_vector_rsp_dot

    subroutine test_vector_rsp_scal(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error

        !> Test vector.
        type(vector_rsp), allocatable :: x, y
        real(sp) :: alpha

        ! Initialize vector.
        x = vector_rsp() ; call x%rand()
        y = x
        call random_number(alpha)
        
        ! Scale the vector.
        call x%scal(alpha)

        ! Check correctness.
        call check(error, norm2(x%data - alpha*y%data) < rtol_sp)

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
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        !> Test vector.
        type(vector_rdp), allocatable :: x
        real(dp) :: alpha

        ! Initialize vector.
        x = vector_rdp() ; call x%rand()
        
        ! Compute its norm.
        alpha = x%norm()

        ! Check result.
        call check(error, is_close(alpha, norm2(x%data)))
        
        return
    end subroutine test_vector_rdp_norm

    subroutine test_vector_rdp_add(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error

        !> Test vectors.
        type(vector_rdp), allocatable :: x, y, z

        ! Initialize vectors.
        x = vector_rdp() ; call x%rand()
        y = vector_rdp() ; call y%rand()
        z = x

        ! Vector addition.
        call z%add(y)

        ! Check correctness.
        call check(error, norm2(z%data - x%data - y%data) < rtol_dp)

        return
    end subroutine test_vector_rdp_add
 
    subroutine test_vector_rdp_sub(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error

        !> Test vectors.
        type(vector_rdp), allocatable :: x, y, z

        ! Initialize vectors.
        x = vector_rdp() ; call x%rand()
        y = vector_rdp() ; call y%rand()
        z = x

        ! Vector addition.
        call z%sub(y)

        ! Check correctness.
        call check(error, norm2(z%data - x%data + y%data) < rtol_dp)

        return
    end subroutine test_vector_rdp_sub

    subroutine test_vector_rdp_dot(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error

        !> Test vectors.
        type(vector_rdp), allocatable :: x, y
        real(dp) :: alpha

        ! Initialize vectors.
        x = vector_rdp() ; call x%rand()
        y = vector_rdp() ; call y%rand()

        ! Compute inner-product.
        alpha = x%dot(y)

        ! Check correctness.
        call check(error, abs(alpha - dot_product(x%data, y%data)) < rtol_dp)

        return
    end subroutine test_vector_rdp_dot

    subroutine test_vector_rdp_scal(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error

        !> Test vector.
        type(vector_rdp), allocatable :: x, y
        real(dp) :: alpha

        ! Initialize vector.
        x = vector_rdp() ; call x%rand()
        y = x
        call random_number(alpha)
        
        ! Scale the vector.
        call x%scal(alpha)

        ! Check correctness.
        call check(error, norm2(x%data - alpha*y%data) < rtol_dp)

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
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        !> Test vector.
        type(vector_csp), allocatable :: x
        real(sp) :: alpha

        ! Initialize vector.
        x = vector_csp() ; call x%rand()
        
        ! Compute its norm.
        alpha = x%norm()

        ! Check result.
        call check(error, is_close(alpha, sqrt(sum(x%data%re**2 + x%data%im**2))))
        
        return
    end subroutine test_vector_csp_norm

    subroutine test_vector_csp_add(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error

        !> Test vectors.
        type(vector_csp), allocatable :: x, y, z

        ! Initialize vectors.
        x = vector_csp() ; call x%rand()
        y = vector_csp() ; call y%rand()
        z = x

        ! Vector addition.
        call z%add(y)

        ! Check correctness.
        call check(error, sqrt(sum((z%data%re - x%data%re - y%data%re)**2 + (z%data%im - x%data%im - y%data%im)**2)) < rtol_sp)

        return
    end subroutine test_vector_csp_add
 
    subroutine test_vector_csp_sub(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error

        !> Test vectors.
        type(vector_csp), allocatable :: x, y, z

        ! Initialize vectors.
        x = vector_csp() ; call x%rand()
        y = vector_csp() ; call y%rand()
        z = x

        ! Vector addition.
        call z%sub(y)

        ! Check correctness.
        call check(error, sqrt(sum((z%data%re - x%data%re + y%data%re)**2 + (z%data%im - x%data%im + y%data%im)**2)) < rtol_sp)

        return
    end subroutine test_vector_csp_sub

    subroutine test_vector_csp_dot(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error

        !> Test vectors.
        type(vector_csp), allocatable :: x, y
        complex(sp) :: alpha

        ! Initialize vectors.
        x = vector_csp() ; call x%rand()
        y = vector_csp() ; call y%rand()

        ! Compute inner-product.
        alpha = x%dot(y)

        ! Check correctness.
        call check(error, abs(alpha - dot_product(x%data, y%data)) < rtol_sp)

        return
    end subroutine test_vector_csp_dot

    subroutine test_vector_csp_scal(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error

        !> Test vector.
        type(vector_csp), allocatable :: x, y
        complex(sp) :: alpha
        complex(sp) :: tmp(test_size)

        ! Initialize vector.
        x = vector_csp() ; call x%rand()
        y = x
        call random_number(alpha%re) ; call random_number(alpha%im)
        
        ! Scale the vector.
        call x%scal(alpha)

        ! Check correctness.
        tmp = x%data - alpha*y%data
        call check(error, sqrt(sum(tmp%re**2 + tmp%im**2)) < rtol_sp)

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
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        !> Test vector.
        type(vector_cdp), allocatable :: x
        real(dp) :: alpha

        ! Initialize vector.
        x = vector_cdp() ; call x%rand()
        
        ! Compute its norm.
        alpha = x%norm()

        ! Check result.
        call check(error, is_close(alpha, sqrt(sum(x%data%re**2 + x%data%im**2))))
        
        return
    end subroutine test_vector_cdp_norm

    subroutine test_vector_cdp_add(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error

        !> Test vectors.
        type(vector_cdp), allocatable :: x, y, z

        ! Initialize vectors.
        x = vector_cdp() ; call x%rand()
        y = vector_cdp() ; call y%rand()
        z = x

        ! Vector addition.
        call z%add(y)

        ! Check correctness.
        call check(error, sqrt(sum((z%data%re - x%data%re - y%data%re)**2 + (z%data%im - x%data%im - y%data%im)**2)) < rtol_dp)

        return
    end subroutine test_vector_cdp_add
 
    subroutine test_vector_cdp_sub(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error

        !> Test vectors.
        type(vector_cdp), allocatable :: x, y, z

        ! Initialize vectors.
        x = vector_cdp() ; call x%rand()
        y = vector_cdp() ; call y%rand()
        z = x

        ! Vector addition.
        call z%sub(y)

        ! Check correctness.
        call check(error, sqrt(sum((z%data%re - x%data%re + y%data%re)**2 + (z%data%im - x%data%im + y%data%im)**2)) < rtol_dp)

        return
    end subroutine test_vector_cdp_sub

    subroutine test_vector_cdp_dot(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error

        !> Test vectors.
        type(vector_cdp), allocatable :: x, y
        complex(dp) :: alpha

        ! Initialize vectors.
        x = vector_cdp() ; call x%rand()
        y = vector_cdp() ; call y%rand()

        ! Compute inner-product.
        alpha = x%dot(y)

        ! Check correctness.
        call check(error, abs(alpha - dot_product(x%data, y%data)) < rtol_dp)

        return
    end subroutine test_vector_cdp_dot

    subroutine test_vector_cdp_scal(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error

        !> Test vector.
        type(vector_cdp), allocatable :: x, y
        complex(dp) :: alpha
        complex(dp) :: tmp(test_size)

        ! Initialize vector.
        x = vector_cdp() ; call x%rand()
        y = x
        call random_number(alpha%re) ; call random_number(alpha%im)
        
        ! Scale the vector.
        call x%scal(alpha)

        ! Check correctness.
        tmp = x%data - alpha*y%data
        call check(error, sqrt(sum(tmp%re**2 + tmp%im**2)) < rtol_dp)

        return
    end subroutine test_vector_cdp_scal


end module TestVectors
