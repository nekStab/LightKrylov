module LightKrylov_TestUtils
    ! Fortran Standard Library
    use stdlib_stats_distribution_normal, only: normal => rvs_normal
    use stdlib_optval, only: optval
    use stdlib_linalg, only: eye
    ! LightKrylov
    use LightKrylov
    use LightKrylov_Constants
    
    implicit none
    
    private

    character(len=128), parameter, private :: this_module = 'LightKrylov_TestUtils'

    integer, parameter, public :: test_size = 128

    public :: get_data
    public :: put_data
    public :: init_rand
    public :: get_err_str

    !-----------------------------------------------
    !-----     TEST VECTOR TYPE DEFINITION     -----
    !-----------------------------------------------

    type, extends(abstract_vector_rsp), public :: vector_rsp
        real(sp), dimension(test_size) :: data = 0.0_sp
    contains
        private
        procedure, pass(self), public :: zero => init_zero_rsp
        procedure, pass(self), public :: dot => dot_rsp
        procedure, pass(self), public :: scal => scal_rsp
        procedure, pass(self), public :: axpby => axpby_rsp
        procedure, pass(self), public :: rand => rand_rsp
        procedure, pass(self), public :: get_size => get_size_rsp
    end type vector_rsp

    type, extends(abstract_vector_rdp), public :: vector_rdp
        real(dp), dimension(test_size) :: data = 0.0_dp
    contains
        private
        procedure, pass(self), public :: zero => init_zero_rdp
        procedure, pass(self), public :: dot => dot_rdp
        procedure, pass(self), public :: scal => scal_rdp
        procedure, pass(self), public :: axpby => axpby_rdp
        procedure, pass(self), public :: rand => rand_rdp
        procedure, pass(self), public :: get_size => get_size_rdp
    end type vector_rdp

    type, extends(abstract_vector_csp), public :: vector_csp
        complex(sp), dimension(test_size) :: data = 0.0_sp
    contains
        private
        procedure, pass(self), public :: zero => init_zero_csp
        procedure, pass(self), public :: dot => dot_csp
        procedure, pass(self), public :: scal => scal_csp
        procedure, pass(self), public :: axpby => axpby_csp
        procedure, pass(self), public :: rand => rand_csp
        procedure, pass(self), public :: get_size => get_size_csp
    end type vector_csp

    type, extends(abstract_vector_cdp), public :: vector_cdp
        complex(dp), dimension(test_size) :: data = 0.0_dp
    contains
        private
        procedure, pass(self), public :: zero => init_zero_cdp
        procedure, pass(self), public :: dot => dot_cdp
        procedure, pass(self), public :: scal => scal_cdp
        procedure, pass(self), public :: axpby => axpby_cdp
        procedure, pass(self), public :: rand => rand_cdp
        procedure, pass(self), public :: get_size => get_size_cdp
    end type vector_cdp


    !----------------------------------------------
    !-----     TEST LINOP TYPE DEFINITION     -----
    !----------------------------------------------

    type, extends(abstract_linop_rsp), public :: linop_rsp
        real(sp), dimension(test_size, test_size) :: data = 0.0_sp
    contains
        private
        procedure, pass(self), public :: matvec  => matvec_rsp
        procedure, pass(self), public :: rmatvec => rmatvec_rsp
    end type

    type, extends(abstract_sym_linop_rsp), public :: spd_linop_rsp
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

    type, extends(abstract_sym_linop_rdp), public :: spd_linop_rdp
        real(dp), dimension(test_size, test_size) :: data = 0.0_dp
    contains
        private
        procedure, pass(self), public :: matvec => sdp_matvec_rdp
        procedure, pass(self), public :: rmatvec => sdp_matvec_rdp
    end type

    type, extends(abstract_linop_csp), public :: linop_csp
        complex(sp), dimension(test_size, test_size) :: data = cmplx(0.0_sp, 0.0_sp, kind=sp)
    contains
        private
        procedure, pass(self), public :: matvec  => matvec_csp
        procedure, pass(self), public :: rmatvec => rmatvec_csp
    end type

    type, extends(abstract_hermitian_linop_csp), public :: hermitian_linop_csp
        complex(sp), dimension(test_size, test_size) :: data = cmplx(0.0_sp, 0.0_sp, kind=sp)
    contains
        private
        procedure, pass(self), public :: matvec => hermitian_matvec_csp
        procedure, pass(self), public :: rmatvec => hermitian_matvec_csp
    end type

    type, extends(abstract_linop_cdp), public :: linop_cdp
        complex(dp), dimension(test_size, test_size) :: data = cmplx(0.0_dp, 0.0_dp, kind=dp)
    contains
        private
        procedure, pass(self), public :: matvec  => matvec_cdp
        procedure, pass(self), public :: rmatvec => rmatvec_cdp
    end type

    type, extends(abstract_hermitian_linop_cdp), public :: hermitian_linop_cdp
        complex(dp), dimension(test_size, test_size) :: data = cmplx(0.0_dp, 0.0_dp, kind=dp)
    contains
        private
        procedure, pass(self), public :: matvec => hermitian_matvec_cdp
        procedure, pass(self), public :: rmatvec => hermitian_matvec_cdp
    end type


    interface get_data
        module procedure get_data_vec_rsp
        module procedure get_data_vec_basis_rsp
        module procedure get_data_linop_rsp
        module procedure get_data_vec_rdp
        module procedure get_data_vec_basis_rdp
        module procedure get_data_linop_rdp
        module procedure get_data_vec_csp
        module procedure get_data_vec_basis_csp
        module procedure get_data_linop_csp
        module procedure get_data_vec_cdp
        module procedure get_data_vec_basis_cdp
        module procedure get_data_linop_cdp
    end interface

    interface put_data
        module procedure put_data_vec_rsp
        module procedure put_data_vec_basis_rsp
        module procedure put_data_linop_rsp
        module procedure put_data_vec_rdp
        module procedure put_data_vec_basis_rdp
        module procedure put_data_linop_rdp
        module procedure put_data_vec_csp
        module procedure put_data_vec_basis_csp
        module procedure put_data_linop_csp
        module procedure put_data_vec_cdp
        module procedure put_data_vec_basis_cdp
        module procedure put_data_linop_cdp
    end interface

    interface init_rand
        module procedure init_rand_vec_rsp
        module procedure init_rand_basis_rsp
        module procedure init_rand_linop_rsp
        module procedure init_rand_spd_linop_rsp
        module procedure init_rand_vec_rdp
        module procedure init_rand_basis_rdp
        module procedure init_rand_linop_rdp
        module procedure init_rand_spd_linop_rdp
        module procedure init_rand_vec_csp
        module procedure init_rand_basis_csp
        module procedure init_rand_linop_csp
        module procedure init_rand_hermitian_linop_csp
        module procedure init_rand_vec_cdp
        module procedure init_rand_basis_cdp
        module procedure init_rand_linop_cdp
        module procedure init_rand_hermitian_linop_cdp
    end interface

    interface get_err_str
        module procedure get_err_str_sp
        module procedure get_err_str_dp
    end interface

contains
    
    !----------------------------------------------------------
    !-----     TYPE-BOUND PROCEDURES FOR TEST VECTORS     -----
    !----------------------------------------------------------

    subroutine init_zero_rsp(self)
        class(vector_rsp), intent(inout) :: self
        self%data = 0.0_sp
        return
    end subroutine init_zero_rsp

    function dot_rsp(self, vec) result(alpha)
        class(vector_rsp), intent(in) :: self
        class(abstract_vector_rsp), intent(in) :: vec
        real(sp) :: alpha

        select type(vec)
        type is(vector_rsp)
            alpha = dot_product(self%data, vec%data)
        end select
    end function dot_rsp

    integer function get_size_rsp(self) result(N)
        class(vector_rsp), intent(in) :: self
        N = test_size
        return
    end function get_size_rsp

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
        real(sp) :: mu(test_size), var(test_size)
        real(sp) :: alpha
 
        mu = 0.0_sp
        var = 1.0_sp
        self%data = normal(mu, var)
 
        normalized = optval(ifnorm, .false.)
        if (normalized) then
            alpha = self%norm()
            call self%scal(1.0_sp/alpha)
        endif
    end subroutine rand_rsp

    subroutine init_zero_rdp(self)
        class(vector_rdp), intent(inout) :: self
        self%data = 0.0_dp
        return
    end subroutine init_zero_rdp

    function dot_rdp(self, vec) result(alpha)
        class(vector_rdp), intent(in) :: self
        class(abstract_vector_rdp), intent(in) :: vec
        real(dp) :: alpha

        select type(vec)
        type is(vector_rdp)
            alpha = dot_product(self%data, vec%data)
        end select
    end function dot_rdp

    integer function get_size_rdp(self) result(N)
        class(vector_rdp), intent(in) :: self
        N = test_size
        return
    end function get_size_rdp

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
        real(dp) :: mu(test_size), var(test_size)
        real(dp) :: alpha
 
        mu = 0.0_dp
        var = 1.0_dp
        self%data = normal(mu, var)
 
        normalized = optval(ifnorm, .false.)
        if (normalized) then
            alpha = self%norm()
            call self%scal(1.0_dp/alpha)
        endif
    end subroutine rand_rdp

    subroutine init_zero_csp(self)
        class(vector_csp), intent(inout) :: self
        self%data = 0.0_sp
        return
    end subroutine init_zero_csp

    function dot_csp(self, vec) result(alpha)
        class(vector_csp), intent(in) :: self
        class(abstract_vector_csp), intent(in) :: vec
        complex(sp) :: alpha

        select type(vec)
        type is(vector_csp)
            alpha = dot_product(self%data, vec%data)
        end select
    end function dot_csp

    integer function get_size_csp(self) result(N)
        class(vector_csp), intent(in) :: self
        N = test_size
        return
    end function get_size_csp

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
        complex(sp) :: mu(test_size), var(test_size)
        complex(sp) :: alpha
 
        mu = 0.0_sp
        var = cmplx(1.0_sp, 1.0_sp, kind=sp)
        self%data = normal(mu, var)
 
        normalized = optval(ifnorm, .false.)
        if (normalized) then
            alpha = self%norm()
            call self%scal(1.0_sp/alpha)
        endif
    end subroutine rand_csp

    subroutine init_zero_cdp(self)
        class(vector_cdp), intent(inout) :: self
        self%data = 0.0_dp
        return
    end subroutine init_zero_cdp

    function dot_cdp(self, vec) result(alpha)
        class(vector_cdp), intent(in) :: self
        class(abstract_vector_cdp), intent(in) :: vec
        complex(dp) :: alpha

        select type(vec)
        type is(vector_cdp)
            alpha = dot_product(self%data, vec%data)
        end select
    end function dot_cdp

    integer function get_size_cdp(self) result(N)
        class(vector_cdp), intent(in) :: self
        N = test_size
        return
    end function get_size_cdp

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
        complex(dp) :: mu(test_size), var(test_size)
        complex(dp) :: alpha
 
        mu = 0.0_dp
        var = cmplx(1.0_dp, 1.0_dp, kind=dp)
        self%data = normal(mu, var)
 
        normalized = optval(ifnorm, .false.)
        if (normalized) then
            alpha = self%norm()
            call self%scal(1.0_dp/alpha)
        endif
    end subroutine rand_cdp


    !---------------------------------------------------------
    !-----     TYPE-BOUND PROCEDURES FOR TEST LINOPS     -----
    !---------------------------------------------------------

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


    !----------------------------------------------------
    !-----     EXTRACT DATA FROM ABSTRACT TYPES     -----
    !----------------------------------------------------

    subroutine get_data_vec_rsp(vec_out, vec_in)
        real(sp), intent(out) :: vec_out(:)
        type(vector_rsp), intent(in) :: vec_in
        vec_out = vec_in%data
        return
    end subroutine get_data_vec_rsp

    subroutine get_data_vec_basis_rsp(basis_out, basis_in)
        real(sp), intent(out) :: basis_out(:, :)
        type(vector_rsp), intent(in) :: basis_in(:)
        ! Internal variables.
        integer :: k
        do k = 1, size(basis_in)
            basis_out(:, k) = basis_in(k)%data
        enddo
        return
    end subroutine get_data_vec_basis_rsp

    subroutine get_data_linop_rsp(mat_out, linop_in)
        real(sp), intent(out) :: mat_out(:, :)
        type(linop_rsp), intent(in) :: linop_in
        mat_out = linop_in%data
        return
    end subroutine get_data_linop_rsp

    subroutine get_data_vec_rdp(vec_out, vec_in)
        real(dp), intent(out) :: vec_out(:)
        type(vector_rdp), intent(in) :: vec_in
        vec_out = vec_in%data
        return
    end subroutine get_data_vec_rdp

    subroutine get_data_vec_basis_rdp(basis_out, basis_in)
        real(dp), intent(out) :: basis_out(:, :)
        type(vector_rdp), intent(in) :: basis_in(:)
        ! Internal variables.
        integer :: k
        do k = 1, size(basis_in)
            basis_out(:, k) = basis_in(k)%data
        enddo
        return
    end subroutine get_data_vec_basis_rdp

    subroutine get_data_linop_rdp(mat_out, linop_in)
        real(dp), intent(out) :: mat_out(:, :)
        type(linop_rdp), intent(in) :: linop_in
        mat_out = linop_in%data
        return
    end subroutine get_data_linop_rdp

    subroutine get_data_vec_csp(vec_out, vec_in)
        complex(sp), intent(out) :: vec_out(:)
        type(vector_csp), intent(in) :: vec_in
        vec_out = vec_in%data
        return
    end subroutine get_data_vec_csp

    subroutine get_data_vec_basis_csp(basis_out, basis_in)
        complex(sp), intent(out) :: basis_out(:, :)
        type(vector_csp), intent(in) :: basis_in(:)
        ! Internal variables.
        integer :: k
        do k = 1, size(basis_in)
            basis_out(:, k) = basis_in(k)%data
        enddo
        return
    end subroutine get_data_vec_basis_csp

    subroutine get_data_linop_csp(mat_out, linop_in)
        complex(sp), intent(out) :: mat_out(:, :)
        type(linop_csp), intent(in) :: linop_in
        mat_out = linop_in%data
        return
    end subroutine get_data_linop_csp

    subroutine get_data_vec_cdp(vec_out, vec_in)
        complex(dp), intent(out) :: vec_out(:)
        type(vector_cdp), intent(in) :: vec_in
        vec_out = vec_in%data
        return
    end subroutine get_data_vec_cdp

    subroutine get_data_vec_basis_cdp(basis_out, basis_in)
        complex(dp), intent(out) :: basis_out(:, :)
        type(vector_cdp), intent(in) :: basis_in(:)
        ! Internal variables.
        integer :: k
        do k = 1, size(basis_in)
            basis_out(:, k) = basis_in(k)%data
        enddo
        return
    end subroutine get_data_vec_basis_cdp

    subroutine get_data_linop_cdp(mat_out, linop_in)
        complex(dp), intent(out) :: mat_out(:, :)
        type(linop_cdp), intent(in) :: linop_in
        mat_out = linop_in%data
        return
    end subroutine get_data_linop_cdp


    !----------------------------------------------
    !-----     PUT DATA TO ABSTRACT TYPES     -----
    !----------------------------------------------

    subroutine put_data_vec_rsp(vec_out, vec_in)
        type(vector_rsp), intent(out) :: vec_out
        real(sp), intent(in) :: vec_in
        vec_out%data = vec_in
        return
    end subroutine put_data_vec_rsp

    subroutine put_data_vec_basis_rsp(basis_out, basis_in)
        type(vector_rsp), intent(out) :: basis_out(:)
        real(sp), intent(in) :: basis_in(:, :)
        ! Internal variables.
        integer :: k
        do k = 1, size(basis_out)
            basis_out(k)%data = basis_in(:, k)
        enddo
        return
    end subroutine put_data_vec_basis_rsp

    subroutine put_data_linop_rsp(linop_out, mat_in)
        type(linop_rsp), intent(out) :: linop_out
        real(sp), intent(in) :: mat_in(:, :)
        ! Internal variables.
        linop_out%data = mat_in
        return
    end subroutine put_data_linop_rsp

    subroutine put_data_vec_rdp(vec_out, vec_in)
        type(vector_rdp), intent(out) :: vec_out
        real(dp), intent(in) :: vec_in
        vec_out%data = vec_in
        return
    end subroutine put_data_vec_rdp

    subroutine put_data_vec_basis_rdp(basis_out, basis_in)
        type(vector_rdp), intent(out) :: basis_out(:)
        real(dp), intent(in) :: basis_in(:, :)
        ! Internal variables.
        integer :: k
        do k = 1, size(basis_out)
            basis_out(k)%data = basis_in(:, k)
        enddo
        return
    end subroutine put_data_vec_basis_rdp

    subroutine put_data_linop_rdp(linop_out, mat_in)
        type(linop_rdp), intent(out) :: linop_out
        real(dp), intent(in) :: mat_in(:, :)
        ! Internal variables.
        linop_out%data = mat_in
        return
    end subroutine put_data_linop_rdp

    subroutine put_data_vec_csp(vec_out, vec_in)
        type(vector_csp), intent(out) :: vec_out
        complex(sp), intent(in) :: vec_in
        vec_out%data = vec_in
        return
    end subroutine put_data_vec_csp

    subroutine put_data_vec_basis_csp(basis_out, basis_in)
        type(vector_csp), intent(out) :: basis_out(:)
        complex(sp), intent(in) :: basis_in(:, :)
        ! Internal variables.
        integer :: k
        do k = 1, size(basis_out)
            basis_out(k)%data = basis_in(:, k)
        enddo
        return
    end subroutine put_data_vec_basis_csp

    subroutine put_data_linop_csp(linop_out, mat_in)
        type(linop_csp), intent(out) :: linop_out
        complex(sp), intent(in) :: mat_in(:, :)
        ! Internal variables.
        linop_out%data = mat_in
        return
    end subroutine put_data_linop_csp

    subroutine put_data_vec_cdp(vec_out, vec_in)
        type(vector_cdp), intent(out) :: vec_out
        complex(dp), intent(in) :: vec_in
        vec_out%data = vec_in
        return
    end subroutine put_data_vec_cdp

    subroutine put_data_vec_basis_cdp(basis_out, basis_in)
        type(vector_cdp), intent(out) :: basis_out(:)
        complex(dp), intent(in) :: basis_in(:, :)
        ! Internal variables.
        integer :: k
        do k = 1, size(basis_out)
            basis_out(k)%data = basis_in(:, k)
        enddo
        return
    end subroutine put_data_vec_basis_cdp

    subroutine put_data_linop_cdp(linop_out, mat_in)
        type(linop_cdp), intent(out) :: linop_out
        complex(dp), intent(in) :: mat_in(:, :)
        ! Internal variables.
        linop_out%data = mat_in
        return
    end subroutine put_data_linop_cdp


    !--------------------------------------------------------------
    !-----     INITIALIZE ABSTRACT TYPES WITH RANDOM DATA     -----
    !--------------------------------------------------------------

    subroutine init_rand_vec_rsp(x)
        type(vector_rsp), intent(inout) :: x
        call x%rand()
        return
    end subroutine init_rand_vec_rsp

    subroutine init_rand_basis_rsp(X)
        type(vector_rsp), intent(inout) :: X(:)
        integer :: i
        do i = 1, size(X)
            call X(i)%rand()
        enddo
        return
    end subroutine init_rand_basis_rsp

    subroutine init_rand_linop_rsp(linop)
        type(linop_rsp), intent(inout) :: linop
        real(sp), allocatable :: mu(:, :), var(:, :)
        allocate(mu(test_size, test_size)) ; mu = 0.0_sp
        allocate(var(test_size, test_size))
        var = 1.0_sp
        linop%data = normal(mu, var)
        return
    end subroutine init_rand_linop_rsp

    subroutine init_rand_spd_linop_rsp(linop)
        type(spd_linop_rsp), intent(inout) :: linop
        real(sp), allocatable :: mu(:, :), var(:, :)
        real(sp), allocatable :: data(:, :)
        allocate(mu(test_size, test_size)) ; mu = zero_rsp
        allocate(var(test_size, test_size)) ; var = one_rsp

        data = normal(mu, var)
        linop%data = matmul(data, transpose(data))/test_size + 0.01*eye(test_size)

        return
    end subroutine init_rand_spd_linop_rsp

    subroutine init_rand_vec_rdp(x)
        type(vector_rdp), intent(inout) :: x
        call x%rand()
        return
    end subroutine init_rand_vec_rdp

    subroutine init_rand_basis_rdp(X)
        type(vector_rdp), intent(inout) :: X(:)
        integer :: i
        do i = 1, size(X)
            call X(i)%rand()
        enddo
        return
    end subroutine init_rand_basis_rdp

    subroutine init_rand_linop_rdp(linop)
        type(linop_rdp), intent(inout) :: linop
        real(dp), allocatable :: mu(:, :), var(:, :)
        allocate(mu(test_size, test_size)) ; mu = 0.0_dp
        allocate(var(test_size, test_size))
        var = 1.0_dp
        linop%data = normal(mu, var)
        return
    end subroutine init_rand_linop_rdp

    subroutine init_rand_spd_linop_rdp(linop)
        type(spd_linop_rdp), intent(inout) :: linop
        real(dp), allocatable :: mu(:, :), var(:, :)
        real(dp), allocatable :: data(:, :)
        allocate(mu(test_size, test_size)) ; mu = zero_rdp
        allocate(var(test_size, test_size)) ; var = one_rdp

        data = normal(mu, var)
        linop%data = matmul(data, transpose(data))/test_size + 0.01*eye(test_size)

        return
    end subroutine init_rand_spd_linop_rdp

    subroutine init_rand_vec_csp(x)
        type(vector_csp), intent(inout) :: x
        call x%rand()
        return
    end subroutine init_rand_vec_csp

    subroutine init_rand_basis_csp(X)
        type(vector_csp), intent(inout) :: X(:)
        integer :: i
        do i = 1, size(X)
            call X(i)%rand()
        enddo
        return
    end subroutine init_rand_basis_csp

    subroutine init_rand_linop_csp(linop)
        type(linop_csp), intent(inout) :: linop
        complex(sp), allocatable :: mu(:, :), var(:, :)
        allocate(mu(test_size, test_size)) ; mu = 0.0_sp
        allocate(var(test_size, test_size))
        var = cmplx(1.0_sp, 1.0_sp, kind=sp)
        linop%data = normal(mu, var)
        return
    end subroutine init_rand_linop_csp

    subroutine init_rand_hermitian_linop_csp(linop)
        type(hermitian_linop_csp), intent(inout) :: linop
        complex(sp), allocatable :: data(:, :)
        complex(sp), allocatable :: mu(:, :), var(:, :)

        allocate(mu(test_size, test_size)) ; mu = 0.0_sp
        allocate(var(test_size, test_size)) ; var = cmplx(1.0_sp, 1.0_sp, kind=sp)

        data = normal(mu, var)
        data = matmul(data, transpose(conjg(data)))/test_size + 0.01*eye(test_size)
        linop%data = data

        return
    end subroutine init_rand_hermitian_linop_csp

    subroutine init_rand_vec_cdp(x)
        type(vector_cdp), intent(inout) :: x
        call x%rand()
        return
    end subroutine init_rand_vec_cdp

    subroutine init_rand_basis_cdp(X)
        type(vector_cdp), intent(inout) :: X(:)
        integer :: i
        do i = 1, size(X)
            call X(i)%rand()
        enddo
        return
    end subroutine init_rand_basis_cdp

    subroutine init_rand_linop_cdp(linop)
        type(linop_cdp), intent(inout) :: linop
        complex(dp), allocatable :: mu(:, :), var(:, :)
        allocate(mu(test_size, test_size)) ; mu = 0.0_dp
        allocate(var(test_size, test_size))
        var = cmplx(1.0_dp, 1.0_dp, kind=dp)
        linop%data = normal(mu, var)
        return
    end subroutine init_rand_linop_cdp

    subroutine init_rand_hermitian_linop_cdp(linop)
        type(hermitian_linop_cdp), intent(inout) :: linop
        complex(dp), allocatable :: data(:, :)
        complex(dp), allocatable :: mu(:, :), var(:, :)

        allocate(mu(test_size, test_size)) ; mu = 0.0_dp
        allocate(var(test_size, test_size)) ; var = cmplx(1.0_dp, 1.0_dp, kind=dp)

        data = normal(mu, var)
        data = matmul(data, transpose(conjg(data)))/test_size + 0.01*eye(test_size)
        linop%data = data

        return
    end subroutine init_rand_hermitian_linop_cdp


    subroutine get_err_str_sp(msg, info, err)
      character(len=*), intent(inout) :: msg
      character(len=*), intent(in)    :: info
      real(sp) :: err

      ! internals
      character*8 :: value_str
      character(len=*), parameter :: indent = repeat(" ", 4)

      write(value_str, '(E8.2)') err
      msg = indent // info // value_str // achar(10)
       
    end subroutine get_err_str_sp
    subroutine get_err_str_dp(msg, info, err)
      character(len=*), intent(inout) :: msg
      character(len=*), intent(in)    :: info
      real(dp) :: err

      ! internals
      character*8 :: value_str
      character(len=*), parameter :: indent = repeat(" ", 4)

      write(value_str, '(E8.2)') err
      msg = indent // info // value_str // achar(10)
       
    end subroutine get_err_str_dp
    
end module LightKrylov_TestUtils
