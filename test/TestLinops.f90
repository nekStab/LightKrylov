module TestLinops
    ! Fortran Standard Library
    use stdlib_math, only: is_close, all_close

    ! LightKrylov
    use LightKrylov
    
    ! Testdrive
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use TestVectors

    implicit none
    private

    public :: collect_linop_rsp_testsuite
    public :: collect_linop_rdp_testsuite
    public :: collect_linop_csp_testsuite
    public :: collect_linop_cdp_testsuite

    type, extends(abstract_linop_rsp), public :: linop_rsp
        real(sp), dimension(test_size, test_size) :: data = 0.0_sp
    contains
        private
        procedure, pass(self), public :: matvec  => matvec_rsp
        procedure, pass(self), public :: rmatvec => rmatvec_rsp
    end type

    type, extends(abstract_spd_linop_rsp), public :: spd_linop_rsp
        real(sp), dimension(test_size, test_size) :: data = 0.0_sp
    contains
        private
        procedure, pass(self), public :: matvec => sdp_matvec_rsp
        procedure, pass(self), public :: rmatvec => sdp_matvec_rsp
    end type

    type, extends(abstract_linop_rdp), public :: linop_rdp
        real(dp), dimension(test_size, test_size) :: data = 0.0_dp
    contains
        private
        procedure, pass(self), public :: matvec  => matvec_rdp
        procedure, pass(self), public :: rmatvec => rmatvec_rdp
    end type

    type, extends(abstract_spd_linop_rdp), public :: spd_linop_rdp
        real(dp), dimension(test_size, test_size) :: data = 0.0_dp
    contains
        private
        procedure, pass(self), public :: matvec => sdp_matvec_rdp
        procedure, pass(self), public :: rmatvec => sdp_matvec_rdp
    end type

    type, extends(abstract_linop_csp), public :: linop_csp
        complex(sp), dimension(test_size, test_size) :: data = 0.0_sp
    contains
        private
        procedure, pass(self), public :: matvec  => matvec_csp
        procedure, pass(self), public :: rmatvec => rmatvec_csp
    end type

    type, extends(abstract_hermitian_linop_csp), public :: hermitian_linop_csp
        complex(sp), dimension(test_size, test_size) :: data = 0.0_sp
    contains
        private
        procedure, pass(self), public :: matvec => hermitian_matvec_csp
        procedure, pass(self), public :: rmatvec => hermitian_matvec_csp
    end type

    type, extends(abstract_linop_cdp), public :: linop_cdp
        complex(dp), dimension(test_size, test_size) :: data = 0.0_dp
    contains
        private
        procedure, pass(self), public :: matvec  => matvec_cdp
        procedure, pass(self), public :: rmatvec => rmatvec_cdp
    end type

    type, extends(abstract_hermitian_linop_cdp), public :: hermitian_linop_cdp
        complex(dp), dimension(test_size, test_size) :: data = 0.0_dp
    contains
        private
        procedure, pass(self), public :: matvec => hermitian_matvec_cdp
        procedure, pass(self), public :: rmatvec => hermitian_matvec_cdp
    end type

contains

    !--------------------------------------------------------------------
    !-----     DEFINITIONS OF THE VARIOUS TYPE-BOUND PROCEDURES     -----
    !--------------------------------------------------------------------

    subroutine matvec_rsp(self, vec_in, vec_out)
        class(linop_rsp), intent(in)  :: self
        class(abstract_vector_rsp)       , intent(in)  :: vec_in
        class(abstract_vector_rsp)       , intent(out) :: vec_out

        select type(vec_in)
        type is(vector_rsp)
            select type(vec_out)
            type is(vector_rsp)

            vec_out%data = matmul(self%data, vec_in%data)

            end select
        end select

        return
    end subroutine matvec_rsp

    subroutine rmatvec_rsp(self, vec_in, vec_out)
        class(linop_rsp), intent(in)  :: self
        class(abstract_vector_rsp)       , intent(in)  :: vec_in
        class(abstract_vector_rsp)       , intent(out) :: vec_out

        select type(vec_in)
        type is(vector_rsp)
            select type(vec_out)
            type is(vector_rsp)

            vec_out%data = matmul(transpose(self%data), vec_in%data)

            end select
        end select

        return
    end subroutine rmatvec_rsp

    subroutine sdp_matvec_rsp(self, vec_in, vec_out)
        class(spd_linop_rsp), intent(in)  :: self
        class(abstract_vector_rsp)       , intent(in)  :: vec_in
        class(abstract_vector_rsp)       , intent(out) :: vec_out

        select type(vec_in)
        type is(vector_rsp)
            select type(vec_out)
            type is(vector_rsp)

            vec_out%data = matmul(self%data, vec_in%data)

            end select
        end select

        return
    end subroutine sdp_matvec_rsp

    subroutine matvec_rdp(self, vec_in, vec_out)
        class(linop_rdp), intent(in)  :: self
        class(abstract_vector_rdp)       , intent(in)  :: vec_in
        class(abstract_vector_rdp)       , intent(out) :: vec_out

        select type(vec_in)
        type is(vector_rdp)
            select type(vec_out)
            type is(vector_rdp)

            vec_out%data = matmul(self%data, vec_in%data)

            end select
        end select

        return
    end subroutine matvec_rdp

    subroutine rmatvec_rdp(self, vec_in, vec_out)
        class(linop_rdp), intent(in)  :: self
        class(abstract_vector_rdp)       , intent(in)  :: vec_in
        class(abstract_vector_rdp)       , intent(out) :: vec_out

        select type(vec_in)
        type is(vector_rdp)
            select type(vec_out)
            type is(vector_rdp)

            vec_out%data = matmul(transpose(self%data), vec_in%data)

            end select
        end select

        return
    end subroutine rmatvec_rdp

    subroutine sdp_matvec_rdp(self, vec_in, vec_out)
        class(spd_linop_rdp), intent(in)  :: self
        class(abstract_vector_rdp)       , intent(in)  :: vec_in
        class(abstract_vector_rdp)       , intent(out) :: vec_out

        select type(vec_in)
        type is(vector_rdp)
            select type(vec_out)
            type is(vector_rdp)

            vec_out%data = matmul(self%data, vec_in%data)

            end select
        end select

        return
    end subroutine sdp_matvec_rdp

    subroutine matvec_csp(self, vec_in, vec_out)
        class(linop_csp), intent(in)  :: self
        class(abstract_vector_csp)       , intent(in)  :: vec_in
        class(abstract_vector_csp)       , intent(out) :: vec_out

        select type(vec_in)
        type is(vector_csp)
            select type(vec_out)
            type is(vector_csp)

            vec_out%data = matmul(self%data, vec_in%data)

            end select
        end select

        return
    end subroutine matvec_csp

    subroutine rmatvec_csp(self, vec_in, vec_out)
        class(linop_csp), intent(in)  :: self
        class(abstract_vector_csp)       , intent(in)  :: vec_in
        class(abstract_vector_csp)       , intent(out) :: vec_out

        select type(vec_in)
        type is(vector_csp)
            select type(vec_out)
            type is(vector_csp)

            vec_out%data = matmul(transpose(conjg(self%data)), vec_in%data)

            end select
        end select

        return
    end subroutine rmatvec_csp

    subroutine hermitian_matvec_csp(self, vec_in, vec_out)
        class(hermitian_linop_csp), intent(in)  :: self
        class(abstract_vector_csp)       , intent(in)  :: vec_in
        class(abstract_vector_csp)       , intent(out) :: vec_out

        select type(vec_in)
        type is(vector_csp)
            select type(vec_out)
            type is(vector_csp)

            vec_out%data = matmul(self%data, vec_in%data)

            end select
        end select

        return
    end subroutine hermitian_matvec_csp

    subroutine matvec_cdp(self, vec_in, vec_out)
        class(linop_cdp), intent(in)  :: self
        class(abstract_vector_cdp)       , intent(in)  :: vec_in
        class(abstract_vector_cdp)       , intent(out) :: vec_out

        select type(vec_in)
        type is(vector_cdp)
            select type(vec_out)
            type is(vector_cdp)

            vec_out%data = matmul(self%data, vec_in%data)

            end select
        end select

        return
    end subroutine matvec_cdp

    subroutine rmatvec_cdp(self, vec_in, vec_out)
        class(linop_cdp), intent(in)  :: self
        class(abstract_vector_cdp)       , intent(in)  :: vec_in
        class(abstract_vector_cdp)       , intent(out) :: vec_out

        select type(vec_in)
        type is(vector_cdp)
            select type(vec_out)
            type is(vector_cdp)

            vec_out%data = matmul(transpose(conjg(self%data)), vec_in%data)

            end select
        end select

        return
    end subroutine rmatvec_cdp

    subroutine hermitian_matvec_cdp(self, vec_in, vec_out)
        class(hermitian_linop_cdp), intent(in)  :: self
        class(abstract_vector_cdp)       , intent(in)  :: vec_in
        class(abstract_vector_cdp)       , intent(out) :: vec_out

        select type(vec_in)
        type is(vector_cdp)
            select type(vec_out)
            type is(vector_cdp)

            vec_out%data = matmul(self%data, vec_in%data)

            end select
        end select

        return
    end subroutine hermitian_matvec_cdp


    !---------------------------------------------------------
    !-----     DEFINITIONS OF THE VARIOUS UNIT TESTS     -----
    !---------------------------------------------------------

    subroutine collect_linop_rsp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("Matrix-vector product", test_matvec_rsp), &
                    new_unittest("Tranpose Matrix-vector product", test_rmatvec_rsp)  &
                    ]
        return
    end subroutine collect_linop_rsp_testsuite

    subroutine test_matvec_rsp(error)
        !> Error to be returned.
        type(error_type), allocatable, intent(out) :: error
        !> Test vectors.
        type(vector_rsp), allocatable :: x, y
        !> Test LinOp.
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

        return
    end subroutine test_matvec_rsp

    subroutine test_rmatvec_rsp(error)
        !> Error to be returned.
        type(error_type), allocatable, intent(out) :: error
        !> Test vectors.
        type(vector_rsp), allocatable :: x, y
        !> Test LinOp.
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
        
        return
    end subroutine test_rmatvec_rsp

    subroutine collect_linop_rdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("Matrix-vector product", test_matvec_rdp), &
                    new_unittest("Tranpose Matrix-vector product", test_rmatvec_rdp)  &
                    ]
        return
    end subroutine collect_linop_rdp_testsuite

    subroutine test_matvec_rdp(error)
        !> Error to be returned.
        type(error_type), allocatable, intent(out) :: error
        !> Test vectors.
        type(vector_rdp), allocatable :: x, y
        !> Test LinOp.
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

        return
    end subroutine test_matvec_rdp

    subroutine test_rmatvec_rdp(error)
        !> Error to be returned.
        type(error_type), allocatable, intent(out) :: error
        !> Test vectors.
        type(vector_rdp), allocatable :: x, y
        !> Test LinOp.
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
        
        return
    end subroutine test_rmatvec_rdp

    subroutine collect_linop_csp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("Matrix-vector product", test_matvec_csp), &
                    new_unittest("Tranpose Matrix-vector product", test_rmatvec_csp)  &
                    ]
        return
    end subroutine collect_linop_csp_testsuite

    subroutine test_matvec_csp(error)
        !> Error to be returned.
        type(error_type), allocatable, intent(out) :: error
        !> Test vectors.
        type(vector_csp), allocatable :: x, y
        !> Test LinOp.
        type(linop_csp), allocatable :: A

        ! Initialize vectors.
        x = vector_csp() ; call x%rand()
        y = vector_csp() ; call y%rand()

        ! Initialize matrix.
        A = linop_csp()
        call random_number(A%data%re) ; call random_number(A%data%im)

        ! Compute matrix-vector product.
        call A%matvec(x, y)

        ! Check result.
        call check(error, norm2(abs(y%data - matmul(A%data, x%data))) < rtol_sp)

        return
    end subroutine test_matvec_csp

    subroutine test_rmatvec_csp(error)
        !> Error to be returned.
        type(error_type), allocatable, intent(out) :: error
        !> Test vectors.
        type(vector_csp), allocatable :: x, y
        !> Test LinOp.
        type(linop_csp), allocatable :: A

        ! Initialize vectors.
        x = vector_csp() ; call x%rand()
        y = vector_csp() ; call y%rand()

        ! Initialize matrix.
        A = linop_csp()
        call random_number(A%data%re) ; call random_number(A%data%im)

        ! Compute matrix-vector product.
        call A%rmatvec(x, y)

        ! Check result.
        call check(error, norm2(abs(y%data - matmul(transpose(conjg(A%data)), x%data))) < rtol_sp)
        
        return
    end subroutine test_rmatvec_csp

    subroutine collect_linop_cdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("Matrix-vector product", test_matvec_cdp), &
                    new_unittest("Tranpose Matrix-vector product", test_rmatvec_cdp)  &
                    ]
        return
    end subroutine collect_linop_cdp_testsuite

    subroutine test_matvec_cdp(error)
        !> Error to be returned.
        type(error_type), allocatable, intent(out) :: error
        !> Test vectors.
        type(vector_cdp), allocatable :: x, y
        !> Test LinOp.
        type(linop_cdp), allocatable :: A

        ! Initialize vectors.
        x = vector_cdp() ; call x%rand()
        y = vector_cdp() ; call y%rand()

        ! Initialize matrix.
        A = linop_cdp()
        call random_number(A%data%re) ; call random_number(A%data%im)

        ! Compute matrix-vector product.
        call A%matvec(x, y)

        ! Check result.
        call check(error, norm2(abs(y%data - matmul(A%data, x%data))) < rtol_dp)

        return
    end subroutine test_matvec_cdp

    subroutine test_rmatvec_cdp(error)
        !> Error to be returned.
        type(error_type), allocatable, intent(out) :: error
        !> Test vectors.
        type(vector_cdp), allocatable :: x, y
        !> Test LinOp.
        type(linop_cdp), allocatable :: A

        ! Initialize vectors.
        x = vector_cdp() ; call x%rand()
        y = vector_cdp() ; call y%rand()

        ! Initialize matrix.
        A = linop_cdp()
        call random_number(A%data%re) ; call random_number(A%data%im)

        ! Compute matrix-vector product.
        call A%rmatvec(x, y)

        ! Check result.
        call check(error, norm2(abs(y%data - matmul(transpose(conjg(A%data)), x%data))) < rtol_dp)
        
        return
    end subroutine test_rmatvec_cdp

end module TestLinops

