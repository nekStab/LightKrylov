module LightKrylov_TestUtils
    ! Fortran Standard Library
    use stdlib_stats_distribution_normal, only: normal => rvs_normal
    use stdlib_optval, only: optval
    use stdlib_linalg, only: eye, hermitian
    ! LightKrylov
    use LightKrylov
    use LightKrylov_Logger
    use LightKrylov_Constants
    
    implicit none
    
    private

    character(len=*), parameter :: this_module      = 'LK_TUtils'
    character(len=*), parameter :: this_module_long = 'LightKrylov_TestUtils'

    integer, parameter, public :: test_size = 128

    public :: get_data
    public :: put_data
    public :: init_rand
    public :: get_err_str

    ! Roessler
    public :: get_state_rsp
    public :: roessler_analytical_fp_rsp
    public :: get_state_rdp
    public :: roessler_analytical_fp_rdp

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
    interface linop_rsp
        pure module function construct_linop_rsp(data) result(A)
            real(sp), dimension(test_size, test_size), intent(in) :: data
            type(linop_rsp) :: A
        end function
    end interface

    type, extends(abstract_sym_linop_rsp), public :: spd_linop_rsp
        real(sp), dimension(test_size, test_size) :: data = 0.0_sp
    contains
        private
        procedure, pass(self), public :: matvec => sdp_matvec_rsp
        procedure, pass(self), public :: rmatvec => sdp_matvec_rsp
    end type
    interface spd_linop_rsp
        pure module function construct_spd_linop_rsp(data) result(A)
            real(sp), dimension(test_size, test_size), intent(in) :: data
            type(spd_linop_rsp) :: A
        end function
    end interface

    interface hermitian_linop_rsp
        pure module function construct_hermitian_linop_rsp(data) result(A)
            real(sp), dimension(test_size, test_size), intent(in) :: data
            type(linop_rsp) :: A
        end function
    end interface

    type, extends(abstract_linop_rdp), public :: linop_rdp
        real(dp), dimension(test_size, test_size) :: data = 0.0_dp
    contains
        private
        procedure, pass(self), public :: matvec  => matvec_rdp
        procedure, pass(self), public :: rmatvec => rmatvec_rdp
    end type
    interface linop_rdp
        pure module function construct_linop_rdp(data) result(A)
            real(dp), dimension(test_size, test_size), intent(in) :: data
            type(linop_rdp) :: A
        end function
    end interface

    type, extends(abstract_sym_linop_rdp), public :: spd_linop_rdp
        real(dp), dimension(test_size, test_size) :: data = 0.0_dp
    contains
        private
        procedure, pass(self), public :: matvec => sdp_matvec_rdp
        procedure, pass(self), public :: rmatvec => sdp_matvec_rdp
    end type
    interface spd_linop_rdp
        pure module function construct_spd_linop_rdp(data) result(A)
            real(dp), dimension(test_size, test_size), intent(in) :: data
            type(spd_linop_rdp) :: A
        end function
    end interface

    interface hermitian_linop_rdp
        pure module function construct_hermitian_linop_rdp(data) result(A)
            real(dp), dimension(test_size, test_size), intent(in) :: data
            type(linop_rdp) :: A
        end function
    end interface

    type, extends(abstract_linop_csp), public :: linop_csp
        complex(sp), dimension(test_size, test_size) :: data = zero_csp
    contains
        private
        procedure, pass(self), public :: matvec  => matvec_csp
        procedure, pass(self), public :: rmatvec => rmatvec_csp
    end type
    interface linop_csp
        pure module function construct_linop_csp(data) result(A)
            complex(sp), dimension(test_size, test_size), intent(in) :: data
            type(linop_csp) :: A
        end function
    end interface

    type, extends(abstract_hermitian_linop_csp), public :: hermitian_linop_csp
        complex(sp), dimension(test_size, test_size) :: data = zero_csp
    contains
        private
        procedure, pass(self), public :: matvec => hermitian_matvec_csp
        procedure, pass(self), public :: rmatvec => hermitian_matvec_csp
    end type
    interface hermitian_linop_csp
        pure module function construct_hermitian_linop_csp(data) result(A)
            complex(sp), dimension(test_size, test_size), intent(in) :: data
            type(linop_csp) :: A
        end function
    end interface

    type, extends(abstract_linop_cdp), public :: linop_cdp
        complex(dp), dimension(test_size, test_size) :: data = zero_cdp
    contains
        private
        procedure, pass(self), public :: matvec  => matvec_cdp
        procedure, pass(self), public :: rmatvec => rmatvec_cdp
    end type
    interface linop_cdp
        pure module function construct_linop_cdp(data) result(A)
            complex(dp), dimension(test_size, test_size), intent(in) :: data
            type(linop_cdp) :: A
        end function
    end interface

    type, extends(abstract_hermitian_linop_cdp), public :: hermitian_linop_cdp
        complex(dp), dimension(test_size, test_size) :: data = zero_cdp
    contains
        private
        procedure, pass(self), public :: matvec => hermitian_matvec_cdp
        procedure, pass(self), public :: rmatvec => hermitian_matvec_cdp
    end type
    interface hermitian_linop_cdp
        pure module function construct_hermitian_linop_cdp(data) result(A)
            complex(dp), dimension(test_size, test_size), intent(in) :: data
            type(linop_cdp) :: A
        end function
    end interface


    ! ROESSLER SYSTEM

    real(sp), parameter :: a_sp = 0.2_sp
    real(sp), parameter :: b_sp = 0.2_sp
    real(sp), parameter :: c_sp = 5.7_sp

    type, extends(abstract_vector_rsp), public :: state_vector_rsp
       real(sp) :: x = 0.0_sp
       real(sp) :: y = 0.0_sp
       real(sp) :: z = 0.0_sp
    contains
       private
       procedure, pass(self), public :: zero => zero_state_rsp
       procedure, pass(self), public :: dot => dot_state_rsp
       procedure, pass(self), public :: scal => scal_state_rsp
       procedure, pass(self), public :: axpby => axpby_state_rsp
       procedure, pass(self), public :: rand => rand_state_rsp
       procedure, pass(self), public :: get_size => get_size_state_rsp
    end type state_vector_rsp

    type, extends(abstract_system_rsp), public :: roessler_rsp
    contains
       private
       procedure, pass(self), public :: response => eval_roessler_rsp
    end type roessler_rsp

    type, extends(abstract_jacobian_linop_rsp), public :: jacobian_rsp
    contains
       private
       procedure, pass(self), public :: matvec => lin_roessler_rsp
       procedure, pass(self), public :: rmatvec => adj_lin_roessler_rsp
    end type jacobian_rsp
    real(dp), parameter :: a_dp = 0.2_dp
    real(dp), parameter :: b_dp = 0.2_dp
    real(dp), parameter :: c_dp = 5.7_dp

    type, extends(abstract_vector_rdp), public :: state_vector_rdp
       real(dp) :: x = 0.0_dp
       real(dp) :: y = 0.0_dp
       real(dp) :: z = 0.0_dp
    contains
       private
       procedure, pass(self), public :: zero => zero_state_rdp
       procedure, pass(self), public :: dot => dot_state_rdp
       procedure, pass(self), public :: scal => scal_state_rdp
       procedure, pass(self), public :: axpby => axpby_state_rdp
       procedure, pass(self), public :: rand => rand_state_rdp
       procedure, pass(self), public :: get_size => get_size_state_rdp
    end type state_vector_rdp

    type, extends(abstract_system_rdp), public :: roessler_rdp
    contains
       private
       procedure, pass(self), public :: response => eval_roessler_rdp
    end type roessler_rdp

    type, extends(abstract_jacobian_linop_rdp), public :: jacobian_rdp
    contains
       private
       procedure, pass(self), public :: matvec => lin_roessler_rdp
       procedure, pass(self), public :: rmatvec => adj_lin_roessler_rdp
    end type jacobian_rdp

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

    !--------------------------------
    !-----     CONSTRUCTORS     -----
    !--------------------------------

    module procedure construct_linop_rsp
    A%data = data
    end procedure

    module procedure construct_spd_linop_rsp
    A%data = data
    end procedure
    module procedure construct_linop_rdp
    A%data = data
    end procedure

    module procedure construct_spd_linop_rdp
    A%data = data
    end procedure
    module procedure construct_linop_csp
    A%data = data
    end procedure

    module procedure construct_hermitian_linop_csp
    A%data = data
    end procedure
    module procedure construct_linop_cdp
    A%data = data
    end procedure

    module procedure construct_hermitian_linop_cdp
    A%data = data
    end procedure

    !----------------------------------------------------------
    !-----     TYPE-BOUND PROCEDURES FOR TEST VECTORS     -----
    !----------------------------------------------------------

    subroutine init_zero_rsp(self)
        class(vector_rsp), intent(inout) :: self
        self%data = 0.0_sp
    end subroutine init_zero_rsp

    function dot_rsp(self, vec) result(alpha)
        class(vector_rsp), intent(in) :: self
        class(abstract_vector_rsp), intent(in) :: vec
        real(sp) :: alpha
        select type(vec)
        type is(vector_rsp)
            alpha = dot_product(self%data, vec%data)
        class default
            call type_error('vec','vector_rsp','IN',this_module,'dot_rsp')
        end select
    end function dot_rsp

    integer function get_size_rsp(self) result(N)
        class(vector_rsp), intent(in) :: self
        N = test_size
    end function get_size_rsp

    subroutine scal_rsp(self, alpha)
        class(vector_rsp), intent(inout) :: self
        real(sp), intent(in) :: alpha
        self%data = alpha * self%data
    end subroutine scal_rsp

    subroutine axpby_rsp(alpha, vec, beta, self)
        class(vector_rsp), intent(inout) :: self
        class(abstract_vector_rsp), intent(in) :: vec
        real(sp), intent(in) :: alpha, beta
        select type(vec)
        type is(vector_rsp)
            self%data = alpha*vec%data + beta*self%data
        class default
            call type_error('vec','vector_rsp','IN',this_module,'axpby_rsp')
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
    end subroutine init_zero_rdp

    function dot_rdp(self, vec) result(alpha)
        class(vector_rdp), intent(in) :: self
        class(abstract_vector_rdp), intent(in) :: vec
        real(dp) :: alpha
        select type(vec)
        type is(vector_rdp)
            alpha = dot_product(self%data, vec%data)
        class default
            call type_error('vec','vector_rdp','IN',this_module,'dot_rdp')
        end select
    end function dot_rdp

    integer function get_size_rdp(self) result(N)
        class(vector_rdp), intent(in) :: self
        N = test_size
    end function get_size_rdp

    subroutine scal_rdp(self, alpha)
        class(vector_rdp), intent(inout) :: self
        real(dp), intent(in) :: alpha
        self%data = alpha * self%data
    end subroutine scal_rdp

    subroutine axpby_rdp(alpha, vec, beta, self)
        class(vector_rdp), intent(inout) :: self
        class(abstract_vector_rdp), intent(in) :: vec
        real(dp), intent(in) :: alpha, beta
        select type(vec)
        type is(vector_rdp)
            self%data = alpha*vec%data + beta*self%data
        class default
            call type_error('vec','vector_rdp','IN',this_module,'axpby_rdp')
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
    end subroutine init_zero_csp

    function dot_csp(self, vec) result(alpha)
        class(vector_csp), intent(in) :: self
        class(abstract_vector_csp), intent(in) :: vec
        complex(sp) :: alpha
        select type(vec)
        type is(vector_csp)
            alpha = dot_product(self%data, vec%data)
        class default
            call type_error('vec','vector_csp','IN',this_module,'dot_csp')
        end select
    end function dot_csp

    integer function get_size_csp(self) result(N)
        class(vector_csp), intent(in) :: self
        N = test_size
    end function get_size_csp

    subroutine scal_csp(self, alpha)
        class(vector_csp), intent(inout) :: self
        complex(sp), intent(in) :: alpha
        self%data = alpha * self%data
    end subroutine scal_csp

    subroutine axpby_csp(alpha, vec, beta, self)
        class(vector_csp), intent(inout) :: self
        class(abstract_vector_csp), intent(in) :: vec
        complex(sp), intent(in) :: alpha, beta
        select type(vec)
        type is(vector_csp)
            self%data = alpha*vec%data + beta*self%data
        class default
            call type_error('vec','vector_csp','IN',this_module,'axpby_csp')
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
    end subroutine init_zero_cdp

    function dot_cdp(self, vec) result(alpha)
        class(vector_cdp), intent(in) :: self
        class(abstract_vector_cdp), intent(in) :: vec
        complex(dp) :: alpha
        select type(vec)
        type is(vector_cdp)
            alpha = dot_product(self%data, vec%data)
        class default
            call type_error('vec','vector_cdp','IN',this_module,'dot_cdp')
        end select
    end function dot_cdp

    integer function get_size_cdp(self) result(N)
        class(vector_cdp), intent(in) :: self
        N = test_size
    end function get_size_cdp

    subroutine scal_cdp(self, alpha)
        class(vector_cdp), intent(inout) :: self
        complex(dp), intent(in) :: alpha
        self%data = alpha * self%data
    end subroutine scal_cdp

    subroutine axpby_cdp(alpha, vec, beta, self)
        class(vector_cdp), intent(inout) :: self
        class(abstract_vector_cdp), intent(in) :: vec
        complex(dp), intent(in) :: alpha, beta
        select type(vec)
        type is(vector_cdp)
            self%data = alpha*vec%data + beta*self%data
        class default
            call type_error('vec','vector_cdp','IN',this_module,'axpby_cdp')
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
        class(linop_rsp), intent(inout)  :: self
        class(abstract_vector_rsp)       , intent(in)  :: vec_in
        class(abstract_vector_rsp)       , intent(out) :: vec_out

        select type(vec_in)
        type is(vector_rsp)
            select type(vec_out)
            type is(vector_rsp)

                vec_out%data = matmul(self%data, vec_in%data)

            class default
                call type_error('vec_out','vector_rsp','OUT',this_module,'matvec_rsp')
            end select
        class default
            call type_error('vec_in','vector_rsp','IN',this_module,'matvec_rsp')
        end select
    end subroutine matvec_rsp

    subroutine rmatvec_rsp(self, vec_in, vec_out)
        class(linop_rsp), intent(inout)  :: self
        class(abstract_vector_rsp)       , intent(in)  :: vec_in
        class(abstract_vector_rsp)       , intent(out) :: vec_out

        select type(vec_in)
        type is(vector_rsp)
            select type(vec_out)
            type is(vector_rsp)
                vec_out%data = matmul(hermitian(self%data), vec_in%data)
            class default
                call type_error('vec_out','vector_rsp','OUT',this_module,'rmatvec_rsp')
            end select
        class default
            call type_error('vec_in','vector_rsp','IN',this_module,'rmatvec_rsp')
        end select
    end subroutine rmatvec_rsp

    subroutine sdp_matvec_rsp(self, vec_in, vec_out)
        class(spd_linop_rsp), intent(inout)  :: self
        class(abstract_vector_rsp)       , intent(in)  :: vec_in
        class(abstract_vector_rsp)       , intent(out) :: vec_out

        select type(vec_in)
        type is(vector_rsp)
            select type(vec_out)
            type is(vector_rsp)

                vec_out%data = matmul(self%data, vec_in%data)

            class default
                call type_error('vec_out','vector_rsp','OUT',this_module,'spd_matvec_rsp')
            end select
        class default
            call type_error('vec_in','vector_rsp','IN',this_module,'spd_matvec_rsp')
        end select
    end subroutine sdp_matvec_rsp

    subroutine matvec_rdp(self, vec_in, vec_out)
        class(linop_rdp), intent(inout)  :: self
        class(abstract_vector_rdp)       , intent(in)  :: vec_in
        class(abstract_vector_rdp)       , intent(out) :: vec_out

        select type(vec_in)
        type is(vector_rdp)
            select type(vec_out)
            type is(vector_rdp)

                vec_out%data = matmul(self%data, vec_in%data)

            class default
                call type_error('vec_out','vector_rdp','OUT',this_module,'matvec_rdp')
            end select
        class default
            call type_error('vec_in','vector_rdp','IN',this_module,'matvec_rdp')
        end select
    end subroutine matvec_rdp

    subroutine rmatvec_rdp(self, vec_in, vec_out)
        class(linop_rdp), intent(inout)  :: self
        class(abstract_vector_rdp)       , intent(in)  :: vec_in
        class(abstract_vector_rdp)       , intent(out) :: vec_out

        select type(vec_in)
        type is(vector_rdp)
            select type(vec_out)
            type is(vector_rdp)
                vec_out%data = matmul(hermitian(self%data), vec_in%data)
            class default
                call type_error('vec_out','vector_rdp','OUT',this_module,'rmatvec_rdp')
            end select
        class default
            call type_error('vec_in','vector_rdp','IN',this_module,'rmatvec_rdp')
        end select
    end subroutine rmatvec_rdp

    subroutine sdp_matvec_rdp(self, vec_in, vec_out)
        class(spd_linop_rdp), intent(inout)  :: self
        class(abstract_vector_rdp)       , intent(in)  :: vec_in
        class(abstract_vector_rdp)       , intent(out) :: vec_out

        select type(vec_in)
        type is(vector_rdp)
            select type(vec_out)
            type is(vector_rdp)

                vec_out%data = matmul(self%data, vec_in%data)

            class default
                call type_error('vec_out','vector_rdp','OUT',this_module,'spd_matvec_rdp')
            end select
        class default
            call type_error('vec_in','vector_rdp','IN',this_module,'spd_matvec_rdp')
        end select
    end subroutine sdp_matvec_rdp

    subroutine matvec_csp(self, vec_in, vec_out)
        class(linop_csp), intent(inout)  :: self
        class(abstract_vector_csp)       , intent(in)  :: vec_in
        class(abstract_vector_csp)       , intent(out) :: vec_out

        select type(vec_in)
        type is(vector_csp)
            select type(vec_out)
            type is(vector_csp)

                vec_out%data = matmul(self%data, vec_in%data)

            class default
                call type_error('vec_out','vector_csp','OUT',this_module,'matvec_csp')
            end select
        class default
            call type_error('vec_in','vector_csp','IN',this_module,'matvec_csp')
        end select
    end subroutine matvec_csp

    subroutine rmatvec_csp(self, vec_in, vec_out)
        class(linop_csp), intent(inout)  :: self
        class(abstract_vector_csp)       , intent(in)  :: vec_in
        class(abstract_vector_csp)       , intent(out) :: vec_out

        select type(vec_in)
        type is(vector_csp)
            select type(vec_out)
            type is(vector_csp)
                vec_out%data = matmul(hermitian(self%data), vec_in%data)
            class default
                call type_error('vec_out','vector_csp','OUT',this_module,'rmatvec_csp')
            end select
        class default
            call type_error('vec_in','vector_csp','IN',this_module,'rmatvec_csp')
        end select
    end subroutine rmatvec_csp

    subroutine hermitian_matvec_csp(self, vec_in, vec_out)
        class(hermitian_linop_csp), intent(inout)  :: self
        class(abstract_vector_csp)       , intent(in)  :: vec_in
        class(abstract_vector_csp)       , intent(out) :: vec_out

        select type(vec_in)
        type is(vector_csp)
            select type(vec_out)
            type is(vector_csp)

                vec_out%data = matmul(self%data, vec_in%data)

            class default
                call type_error('vec_out','vector_csp','OUT',this_module,'hermitian_matvec_csp')
            end select
        class default
            call type_error('vec_in','vector_csp','IN',this_module,'hermitian_matvec_csp')
        end select
    end subroutine hermitian_matvec_csp

    subroutine matvec_cdp(self, vec_in, vec_out)
        class(linop_cdp), intent(inout)  :: self
        class(abstract_vector_cdp)       , intent(in)  :: vec_in
        class(abstract_vector_cdp)       , intent(out) :: vec_out

        select type(vec_in)
        type is(vector_cdp)
            select type(vec_out)
            type is(vector_cdp)

                vec_out%data = matmul(self%data, vec_in%data)

            class default
                call type_error('vec_out','vector_cdp','OUT',this_module,'matvec_cdp')
            end select
        class default
            call type_error('vec_in','vector_cdp','IN',this_module,'matvec_cdp')
        end select
    end subroutine matvec_cdp

    subroutine rmatvec_cdp(self, vec_in, vec_out)
        class(linop_cdp), intent(inout)  :: self
        class(abstract_vector_cdp)       , intent(in)  :: vec_in
        class(abstract_vector_cdp)       , intent(out) :: vec_out

        select type(vec_in)
        type is(vector_cdp)
            select type(vec_out)
            type is(vector_cdp)
                vec_out%data = matmul(hermitian(self%data), vec_in%data)
            class default
                call type_error('vec_out','vector_cdp','OUT',this_module,'rmatvec_cdp')
            end select
        class default
            call type_error('vec_in','vector_cdp','IN',this_module,'rmatvec_cdp')
        end select
    end subroutine rmatvec_cdp

    subroutine hermitian_matvec_cdp(self, vec_in, vec_out)
        class(hermitian_linop_cdp), intent(inout)  :: self
        class(abstract_vector_cdp)       , intent(in)  :: vec_in
        class(abstract_vector_cdp)       , intent(out) :: vec_out

        select type(vec_in)
        type is(vector_cdp)
            select type(vec_out)
            type is(vector_cdp)

                vec_out%data = matmul(self%data, vec_in%data)

            class default
                call type_error('vec_out','vector_cdp','OUT',this_module,'hermitian_matvec_cdp')
            end select
        class default
            call type_error('vec_in','vector_cdp','IN',this_module,'hermitian_matvec_cdp')
        end select
    end subroutine hermitian_matvec_cdp


    !----------------------------------------------------
    !-----     EXTRACT DATA FROM ABSTRACT TYPES     -----
    !----------------------------------------------------

    subroutine get_data_vec_rsp(vec_out, vec_in)
        real(sp), intent(out) :: vec_out(:)
        type(vector_rsp), intent(in) :: vec_in
        vec_out = vec_in%data
    end subroutine get_data_vec_rsp

    subroutine get_data_vec_basis_rsp(basis_out, basis_in)
        real(sp), intent(out) :: basis_out(:, :)
        type(vector_rsp), intent(in) :: basis_in(:)
        ! Internal variables.
        integer :: k
        do k = 1, size(basis_in)
            basis_out(:, k) = basis_in(k)%data
        enddo
    end subroutine get_data_vec_basis_rsp

    subroutine get_data_linop_rsp(mat_out, linop_in)
        real(sp), intent(out) :: mat_out(:, :)
        type(linop_rsp), intent(in) :: linop_in
        mat_out = linop_in%data
    end subroutine get_data_linop_rsp

    subroutine get_data_vec_rdp(vec_out, vec_in)
        real(dp), intent(out) :: vec_out(:)
        type(vector_rdp), intent(in) :: vec_in
        vec_out = vec_in%data
    end subroutine get_data_vec_rdp

    subroutine get_data_vec_basis_rdp(basis_out, basis_in)
        real(dp), intent(out) :: basis_out(:, :)
        type(vector_rdp), intent(in) :: basis_in(:)
        ! Internal variables.
        integer :: k
        do k = 1, size(basis_in)
            basis_out(:, k) = basis_in(k)%data
        enddo
    end subroutine get_data_vec_basis_rdp

    subroutine get_data_linop_rdp(mat_out, linop_in)
        real(dp), intent(out) :: mat_out(:, :)
        type(linop_rdp), intent(in) :: linop_in
        mat_out = linop_in%data
    end subroutine get_data_linop_rdp

    subroutine get_data_vec_csp(vec_out, vec_in)
        complex(sp), intent(out) :: vec_out(:)
        type(vector_csp), intent(in) :: vec_in
        vec_out = vec_in%data
    end subroutine get_data_vec_csp

    subroutine get_data_vec_basis_csp(basis_out, basis_in)
        complex(sp), intent(out) :: basis_out(:, :)
        type(vector_csp), intent(in) :: basis_in(:)
        ! Internal variables.
        integer :: k
        do k = 1, size(basis_in)
            basis_out(:, k) = basis_in(k)%data
        enddo
    end subroutine get_data_vec_basis_csp

    subroutine get_data_linop_csp(mat_out, linop_in)
        complex(sp), intent(out) :: mat_out(:, :)
        type(linop_csp), intent(in) :: linop_in
        mat_out = linop_in%data
    end subroutine get_data_linop_csp

    subroutine get_data_vec_cdp(vec_out, vec_in)
        complex(dp), intent(out) :: vec_out(:)
        type(vector_cdp), intent(in) :: vec_in
        vec_out = vec_in%data
    end subroutine get_data_vec_cdp

    subroutine get_data_vec_basis_cdp(basis_out, basis_in)
        complex(dp), intent(out) :: basis_out(:, :)
        type(vector_cdp), intent(in) :: basis_in(:)
        ! Internal variables.
        integer :: k
        do k = 1, size(basis_in)
            basis_out(:, k) = basis_in(k)%data
        enddo
    end subroutine get_data_vec_basis_cdp

    subroutine get_data_linop_cdp(mat_out, linop_in)
        complex(dp), intent(out) :: mat_out(:, :)
        type(linop_cdp), intent(in) :: linop_in
        mat_out = linop_in%data
    end subroutine get_data_linop_cdp


    !----------------------------------------------
    !-----     PUT DATA TO ABSTRACT TYPES     -----
    !----------------------------------------------

    subroutine put_data_vec_rsp(vec_out, vec_in)
        type(vector_rsp), intent(out) :: vec_out
        real(sp), intent(in) :: vec_in
        vec_out%data = vec_in
    end subroutine put_data_vec_rsp

    subroutine put_data_vec_basis_rsp(basis_out, basis_in)
        type(vector_rsp), intent(out) :: basis_out(:)
        real(sp), intent(in) :: basis_in(:, :)
        ! Internal variables.
        integer :: k
        do k = 1, size(basis_out)
            basis_out(k)%data = basis_in(:, k)
        enddo
    end subroutine put_data_vec_basis_rsp

    subroutine put_data_linop_rsp(linop_out, mat_in)
        type(linop_rsp), intent(out) :: linop_out
        real(sp), intent(in) :: mat_in(:, :)
        ! Internal variables.
        linop_out%data = mat_in
    end subroutine put_data_linop_rsp

    subroutine put_data_vec_rdp(vec_out, vec_in)
        type(vector_rdp), intent(out) :: vec_out
        real(dp), intent(in) :: vec_in
        vec_out%data = vec_in
    end subroutine put_data_vec_rdp

    subroutine put_data_vec_basis_rdp(basis_out, basis_in)
        type(vector_rdp), intent(out) :: basis_out(:)
        real(dp), intent(in) :: basis_in(:, :)
        ! Internal variables.
        integer :: k
        do k = 1, size(basis_out)
            basis_out(k)%data = basis_in(:, k)
        enddo
    end subroutine put_data_vec_basis_rdp

    subroutine put_data_linop_rdp(linop_out, mat_in)
        type(linop_rdp), intent(out) :: linop_out
        real(dp), intent(in) :: mat_in(:, :)
        ! Internal variables.
        linop_out%data = mat_in
    end subroutine put_data_linop_rdp

    subroutine put_data_vec_csp(vec_out, vec_in)
        type(vector_csp), intent(out) :: vec_out
        complex(sp), intent(in) :: vec_in
        vec_out%data = vec_in
    end subroutine put_data_vec_csp

    subroutine put_data_vec_basis_csp(basis_out, basis_in)
        type(vector_csp), intent(out) :: basis_out(:)
        complex(sp), intent(in) :: basis_in(:, :)
        ! Internal variables.
        integer :: k
        do k = 1, size(basis_out)
            basis_out(k)%data = basis_in(:, k)
        enddo
    end subroutine put_data_vec_basis_csp

    subroutine put_data_linop_csp(linop_out, mat_in)
        type(linop_csp), intent(out) :: linop_out
        complex(sp), intent(in) :: mat_in(:, :)
        ! Internal variables.
        linop_out%data = mat_in
    end subroutine put_data_linop_csp

    subroutine put_data_vec_cdp(vec_out, vec_in)
        type(vector_cdp), intent(out) :: vec_out
        complex(dp), intent(in) :: vec_in
        vec_out%data = vec_in
    end subroutine put_data_vec_cdp

    subroutine put_data_vec_basis_cdp(basis_out, basis_in)
        type(vector_cdp), intent(out) :: basis_out(:)
        complex(dp), intent(in) :: basis_in(:, :)
        ! Internal variables.
        integer :: k
        do k = 1, size(basis_out)
            basis_out(k)%data = basis_in(:, k)
        enddo
    end subroutine put_data_vec_basis_cdp

    subroutine put_data_linop_cdp(linop_out, mat_in)
        type(linop_cdp), intent(out) :: linop_out
        complex(dp), intent(in) :: mat_in(:, :)
        ! Internal variables.
        linop_out%data = mat_in
    end subroutine put_data_linop_cdp


    !--------------------------------------------------------------
    !-----     INITIALIZE ABSTRACT TYPES WITH RANDOM DATA     -----
    !--------------------------------------------------------------

    subroutine init_rand_vec_rsp(x)
        type(vector_rsp), intent(inout) :: x
        call x%rand()
    end subroutine init_rand_vec_rsp

    subroutine init_rand_basis_rsp(X)
        type(vector_rsp), intent(inout) :: X(:)
        integer :: i
        do i = 1, size(X)
            call X(i)%rand()
        enddo
    end subroutine init_rand_basis_rsp

    subroutine init_rand_linop_rsp(linop)
        type(linop_rsp), intent(inout) :: linop
        real(sp), allocatable :: mu(:, :), var(:, :)
        allocate(mu(test_size, test_size)) ; mu = 0.0_sp
        allocate(var(test_size, test_size))
        var = 1.0_sp
        linop%data = normal(mu, var)
    end subroutine init_rand_linop_rsp

    subroutine init_rand_spd_linop_rsp(linop)
        type(spd_linop_rsp), intent(inout) :: linop
        real(sp), allocatable :: mu(:, :), var(:, :)
        real(sp), allocatable :: data(:, :)
        allocate(mu(test_size, test_size)) ; mu = zero_rsp
        allocate(var(test_size, test_size)) ; var = one_rsp

        data = normal(mu, var)
        linop%data = matmul(data, transpose(data))/test_size + 0.01_sp*eye(test_size, mold=1.0_sp)

    end subroutine init_rand_spd_linop_rsp

    subroutine init_rand_vec_rdp(x)
        type(vector_rdp), intent(inout) :: x
        call x%rand()
    end subroutine init_rand_vec_rdp

    subroutine init_rand_basis_rdp(X)
        type(vector_rdp), intent(inout) :: X(:)
        integer :: i
        do i = 1, size(X)
            call X(i)%rand()
        enddo
    end subroutine init_rand_basis_rdp

    subroutine init_rand_linop_rdp(linop)
        type(linop_rdp), intent(inout) :: linop
        real(dp), allocatable :: mu(:, :), var(:, :)
        allocate(mu(test_size, test_size)) ; mu = 0.0_dp
        allocate(var(test_size, test_size))
        var = 1.0_dp
        linop%data = normal(mu, var)
    end subroutine init_rand_linop_rdp

    subroutine init_rand_spd_linop_rdp(linop)
        type(spd_linop_rdp), intent(inout) :: linop
        real(dp), allocatable :: mu(:, :), var(:, :)
        real(dp), allocatable :: data(:, :)
        allocate(mu(test_size, test_size)) ; mu = zero_rdp
        allocate(var(test_size, test_size)) ; var = one_rdp

        data = normal(mu, var)
        linop%data = matmul(data, transpose(data))/test_size + 0.01_dp*eye(test_size, mold=1.0_dp)

    end subroutine init_rand_spd_linop_rdp

    subroutine init_rand_vec_csp(x)
        type(vector_csp), intent(inout) :: x
        call x%rand()
    end subroutine init_rand_vec_csp

    subroutine init_rand_basis_csp(X)
        type(vector_csp), intent(inout) :: X(:)
        integer :: i
        do i = 1, size(X)
            call X(i)%rand()
        enddo
    end subroutine init_rand_basis_csp

    subroutine init_rand_linop_csp(linop)
        type(linop_csp), intent(inout) :: linop
        complex(sp), allocatable :: mu(:, :), var(:, :)
        allocate(mu(test_size, test_size)) ; mu = 0.0_sp
        allocate(var(test_size, test_size))
        var = cmplx(1.0_sp, 1.0_sp, kind=sp)
        linop%data = normal(mu, var)
    end subroutine init_rand_linop_csp

    subroutine init_rand_hermitian_linop_csp(linop)
        type(hermitian_linop_csp), intent(inout) :: linop
        complex(sp), allocatable :: data(:, :)
        complex(sp), allocatable :: mu(:, :), var(:, :)

        allocate(mu(test_size, test_size)) ; mu = 0.0_sp
        allocate(var(test_size, test_size)) ; var = cmplx(1.0_sp, 1.0_sp, kind=sp)

        data = normal(mu, var)
        data = matmul(data, hermitian(data))/test_size + 0.01_sp*eye(test_size, mold=1.0_sp)
        linop%data = data

    end subroutine init_rand_hermitian_linop_csp

    subroutine init_rand_vec_cdp(x)
        type(vector_cdp), intent(inout) :: x
        call x%rand()
    end subroutine init_rand_vec_cdp

    subroutine init_rand_basis_cdp(X)
        type(vector_cdp), intent(inout) :: X(:)
        integer :: i
        do i = 1, size(X)
            call X(i)%rand()
        enddo
    end subroutine init_rand_basis_cdp

    subroutine init_rand_linop_cdp(linop)
        type(linop_cdp), intent(inout) :: linop
        complex(dp), allocatable :: mu(:, :), var(:, :)
        allocate(mu(test_size, test_size)) ; mu = 0.0_dp
        allocate(var(test_size, test_size))
        var = cmplx(1.0_dp, 1.0_dp, kind=dp)
        linop%data = normal(mu, var)
    end subroutine init_rand_linop_cdp

    subroutine init_rand_hermitian_linop_cdp(linop)
        type(hermitian_linop_cdp), intent(inout) :: linop
        complex(dp), allocatable :: data(:, :)
        complex(dp), allocatable :: mu(:, :), var(:, :)

        allocate(mu(test_size, test_size)) ; mu = 0.0_dp
        allocate(var(test_size, test_size)) ; var = cmplx(1.0_dp, 1.0_dp, kind=dp)

        data = normal(mu, var)
        data = matmul(data, hermitian(data))/test_size + 0.01_dp*eye(test_size, mold=1.0_dp)
        linop%data = data

    end subroutine init_rand_hermitian_linop_cdp


    subroutine get_err_str_sp(msg, info, err)
        character(len=*), intent(inout) :: msg
        character(len=*), intent(in)    :: info
        real(sp) :: err

        ! internals
        character(len=9) :: value_str
        character(len=*), parameter :: indent = repeat(" ", 4)

        write(value_str, '(E9.2)') err
        msg = indent // info // value_str // achar(10)
       
    end subroutine get_err_str_sp
    subroutine get_err_str_dp(msg, info, err)
        character(len=*), intent(inout) :: msg
        character(len=*), intent(in)    :: info
        real(dp) :: err

        ! internals
        character(len=9) :: value_str
        character(len=*), parameter :: indent = repeat(" ", 4)

        write(value_str, '(E9.2)') err
        msg = indent // info // value_str // achar(10)
       
    end subroutine get_err_str_dp

    !-----------------------------------------------
    !-----     ROESSLER SYSTEM TYPE DEFINITION     -----
    !-----------------------------------------------
 
    subroutine zero_state_rsp(self)
        class(state_vector_rsp), intent(inout) :: self
        self%x = 0.0_sp
        self%y = 0.0_sp
        self%z = 0.0_sp
    end subroutine zero_state_rsp
  
    real(sp) function dot_state_rsp(self, vec) result(alpha)
        class(state_vector_rsp)   , intent(in) :: self
        class(abstract_vector_rsp), intent(in) :: vec
        select type(vec)
        type is(state_vector_rsp)
            alpha = self%x*vec%x + self%y*vec%y + self%z*vec%z
        class default
            call type_error('vec','state_vector_rsp','IN',this_module,'dot_state_rsp')
        end select
    end function dot_state_rsp
  
    subroutine scal_state_rsp(self, alpha)
        class(state_vector_rsp), intent(inout) :: self
        real(sp)               , intent(in)    :: alpha
        self%x = self%x * alpha
        self%y = self%y * alpha
        self%z = self%z * alpha
    end subroutine scal_state_rsp
  
    subroutine axpby_state_rsp(alpha, vec, beta, self)
        class(state_vector_rsp)   , intent(inout) :: self
        class(abstract_vector_rsp), intent(in)    :: vec
        real(sp)                  , intent(in)    :: alpha, beta
        select type(vec)
        type is(state_vector_rsp)
            self%x = beta*self%x + alpha*vec%x
            self%y = beta*self%y + alpha*vec%y
            self%z = beta*self%z + alpha*vec%z
        class default
            call type_error('vec','state_vector_rsp','IN',this_module,'axpby_state_rsp')
        end select
    end subroutine axpby_state_rsp
  
    integer function get_size_state_rsp(self) result(N)
        class(state_vector_rsp), intent(in) :: self
        N = 3
    end function get_size_state_rsp
  
    subroutine rand_state_rsp(self, ifnorm)
        class(state_vector_rsp), intent(inout) :: self
        logical, optional,   intent(in)        :: ifnorm
        logical :: normalized
        real(sp) :: mu, var
        real(sp) :: alpha
    
        mu = zero_rsp
        var = one_rsp
        self%x = normal(mu, var)
        self%y = normal(mu, var)
        self%z = normal(mu, var)
    
        normalized = optval(ifnorm, .false.)
        if (normalized) then
            alpha = self%norm()
            call self%scal(one_rsp/alpha)
        endif
    end subroutine rand_state_rsp

    subroutine zero_state_rdp(self)
        class(state_vector_rdp), intent(inout) :: self
        self%x = 0.0_dp
        self%y = 0.0_dp
        self%z = 0.0_dp
    end subroutine zero_state_rdp
  
    real(dp) function dot_state_rdp(self, vec) result(alpha)
        class(state_vector_rdp)   , intent(in) :: self
        class(abstract_vector_rdp), intent(in) :: vec
        select type(vec)
        type is(state_vector_rdp)
            alpha = self%x*vec%x + self%y*vec%y + self%z*vec%z
        class default
            call type_error('vec','state_vector_rdp','IN',this_module,'dot_state_rdp')
        end select
    end function dot_state_rdp
  
    subroutine scal_state_rdp(self, alpha)
        class(state_vector_rdp), intent(inout) :: self
        real(dp)               , intent(in)    :: alpha
        self%x = self%x * alpha
        self%y = self%y * alpha
        self%z = self%z * alpha
    end subroutine scal_state_rdp
  
    subroutine axpby_state_rdp(alpha, vec, beta, self)
        class(state_vector_rdp)   , intent(inout) :: self
        class(abstract_vector_rdp), intent(in)    :: vec
        real(dp)                  , intent(in)    :: alpha, beta
        select type(vec)
        type is(state_vector_rdp)
            self%x = beta*self%x + alpha*vec%x
            self%y = beta*self%y + alpha*vec%y
            self%z = beta*self%z + alpha*vec%z
        class default
            call type_error('vec','state_vector_rdp','IN',this_module,'axpby_state_rdp')
        end select
    end subroutine axpby_state_rdp
  
    integer function get_size_state_rdp(self) result(N)
        class(state_vector_rdp), intent(in) :: self
        N = 3
    end function get_size_state_rdp
  
    subroutine rand_state_rdp(self, ifnorm)
        class(state_vector_rdp), intent(inout) :: self
        logical, optional,   intent(in)        :: ifnorm
        logical :: normalized
        real(dp) :: mu, var
        real(dp) :: alpha
    
        mu = zero_rdp
        var = one_rdp
        self%x = normal(mu, var)
        self%y = normal(mu, var)
        self%z = normal(mu, var)
    
        normalized = optval(ifnorm, .false.)
        if (normalized) then
            alpha = self%norm()
            call self%scal(one_rdp/alpha)
        endif
    end subroutine rand_state_rdp


    subroutine eval_roessler_rsp(self, vec_in, vec_out, atol)
        class(roessler_rsp),            intent(inout)  :: self
        class(abstract_vector_rsp), intent(in)  :: vec_in
        class(abstract_vector_rsp), intent(out) :: vec_out
        real(sp),                    intent(in)  :: atol

        select type(vec_in)
        type is(state_vector_rsp)
            select type(vec_out)
            type is(state_vector_rsp)

                vec_out%x = -vec_in%y - vec_in%z
                vec_out%y = vec_in%x + a_sp * vec_in%y
                vec_out%z = b_sp + vec_in%z * (vec_in%x - c_sp)

            class default
                call type_error('vec_out','state_vector_rsp','OUT',this_module,'eval_roessler_rsp')
            end select
        class default
            call type_error('vec_in','state_vector_rsp','IN',this_module,'eval_roessler_rsp')
        end select
    end subroutine eval_roessler_rsp

    subroutine get_state_rsp(state, X, Y, Z)
        class(abstract_vector_rsp),   intent(in)  :: state
        real(sp),                 intent(out) :: X, Y, z

        select type (state)
        type is (state_vector_rsp)
            X = state%x
            Y = state%y
            Z = state%z
        class default
            call type_error('state','state_vector_rsp','IN',this_module,'get_state_rsp')
        end select
    end subroutine get_state_rsp

    subroutine lin_roessler_rsp(self, vec_in, vec_out)
        class(jacobian_rsp),            intent(inout)  :: self
        class(abstract_vector_rsp), intent(in)  :: vec_in
        class(abstract_vector_rsp), intent(out) :: vec_out

        real(sp) :: X, Y, Z

        call get_state_rsp(self%X, X, Y, Z)

        select type(vec_in)
        type is(state_vector_rsp)
            select type(vec_out)
            type is(state_vector_rsp)

                vec_out%x = -vec_in%y - vec_in%z
                vec_out%y =  vec_in%x + a_sp*vec_in%y
                vec_out%z =  vec_in%x*Z + vec_in%z*(X - c_sp)

            class default
                call type_error('vec_out','state_vector_rsp','OUT',this_module,'lin_roessler_rsp')
            end select
        class default
            call type_error('vec_in','state_vector_rsp','IN',this_module,'lin_roessler_rsp')
        end select
    end subroutine lin_roessler_rsp

    subroutine adj_lin_roessler_rsp(self, vec_in, vec_out)
        class(jacobian_rsp),            intent(inout)  :: self
        class(abstract_vector_rsp), intent(in)  :: vec_in
        class(abstract_vector_rsp), intent(out) :: vec_out

        real(sp) :: X, Y, Z
        
        call get_state_rsp(self%X, X, Y, Z)

        select type(vec_in)
        type is(state_vector_rsp)
            select type(vec_out)
            type is(state_vector_rsp) 

                vec_out%x =  vec_in%y + vec_in%z*Z
                vec_out%y = -vec_in%x + a_sp*vec_in%y
                vec_out%z = -vec_in%x + vec_in%z*(X - c_sp)

            class default
                call type_error('vec_out','state_vector_rsp','OUT',this_module,'adj_lin_roessler_rsp')
            end select
        class default
            call type_error('vec_in','state_vector_rsp','IN',this_module,'adj_lin_roessler_rsp')
        end select
    end subroutine adj_lin_roessler_rsp

    subroutine roessler_analytical_fp_rsp(fp1, fp2)
        class(state_vector_rsp), intent(out) :: fp1, fp2

        real(sp) :: d

        d = sqrt(c_sp**2 - 4*a_sp*b_sp)

        fp1%x = ( c_sp - d)/ 2
        fp1%y = (-c_sp + d)/(2*a_sp)
        fp1%z = ( c_sp - d)/(2*a_sp)

        fp2%x = ( c_sp + d)/ 2
        fp2%y = (-c_sp - d)/(2*a_sp)
        fp2%z = ( c_sp + d)/(2*a_sp)
    end subroutine roessler_analytical_fp_rsp

    subroutine eval_roessler_rdp(self, vec_in, vec_out, atol)
        class(roessler_rdp),            intent(inout)  :: self
        class(abstract_vector_rdp), intent(in)  :: vec_in
        class(abstract_vector_rdp), intent(out) :: vec_out
        real(dp),                    intent(in)  :: atol

        select type(vec_in)
        type is(state_vector_rdp)
            select type(vec_out)
            type is(state_vector_rdp)

                vec_out%x = -vec_in%y - vec_in%z
                vec_out%y = vec_in%x + a_dp * vec_in%y
                vec_out%z = b_dp + vec_in%z * (vec_in%x - c_dp)

            class default
                call type_error('vec_out','state_vector_rdp','OUT',this_module,'eval_roessler_rdp')
            end select
        class default
            call type_error('vec_in','state_vector_rdp','IN',this_module,'eval_roessler_rdp')
        end select
    end subroutine eval_roessler_rdp

    subroutine get_state_rdp(state, X, Y, Z)
        class(abstract_vector_rdp),   intent(in)  :: state
        real(dp),                 intent(out) :: X, Y, z

        select type (state)
        type is (state_vector_rdp)
            X = state%x
            Y = state%y
            Z = state%z
        class default
            call type_error('state','state_vector_rdp','IN',this_module,'get_state_rdp')
        end select
    end subroutine get_state_rdp

    subroutine lin_roessler_rdp(self, vec_in, vec_out)
        class(jacobian_rdp),            intent(inout)  :: self
        class(abstract_vector_rdp), intent(in)  :: vec_in
        class(abstract_vector_rdp), intent(out) :: vec_out

        real(dp) :: X, Y, Z

        call get_state_rdp(self%X, X, Y, Z)

        select type(vec_in)
        type is(state_vector_rdp)
            select type(vec_out)
            type is(state_vector_rdp)

                vec_out%x = -vec_in%y - vec_in%z
                vec_out%y =  vec_in%x + a_dp*vec_in%y
                vec_out%z =  vec_in%x*Z + vec_in%z*(X - c_dp)

            class default
                call type_error('vec_out','state_vector_rdp','OUT',this_module,'lin_roessler_rdp')
            end select
        class default
            call type_error('vec_in','state_vector_rdp','IN',this_module,'lin_roessler_rdp')
        end select
    end subroutine lin_roessler_rdp

    subroutine adj_lin_roessler_rdp(self, vec_in, vec_out)
        class(jacobian_rdp),            intent(inout)  :: self
        class(abstract_vector_rdp), intent(in)  :: vec_in
        class(abstract_vector_rdp), intent(out) :: vec_out

        real(dp) :: X, Y, Z
        
        call get_state_rdp(self%X, X, Y, Z)

        select type(vec_in)
        type is(state_vector_rdp)
            select type(vec_out)
            type is(state_vector_rdp) 

                vec_out%x =  vec_in%y + vec_in%z*Z
                vec_out%y = -vec_in%x + a_dp*vec_in%y
                vec_out%z = -vec_in%x + vec_in%z*(X - c_dp)

            class default
                call type_error('vec_out','state_vector_rdp','OUT',this_module,'adj_lin_roessler_rdp')
            end select
        class default
            call type_error('vec_in','state_vector_rdp','IN',this_module,'adj_lin_roessler_rdp')
        end select
    end subroutine adj_lin_roessler_rdp

    subroutine roessler_analytical_fp_rdp(fp1, fp2)
        class(state_vector_rdp), intent(out) :: fp1, fp2

        real(dp) :: d

        d = sqrt(c_dp**2 - 4*a_dp*b_dp)

        fp1%x = ( c_dp - d)/ 2
        fp1%y = (-c_dp + d)/(2*a_dp)
        fp1%z = ( c_dp - d)/(2*a_dp)

        fp2%x = ( c_dp + d)/ 2
        fp2%y = (-c_dp - d)/(2*a_dp)
        fp2%z = ( c_dp + d)/(2*a_dp)
    end subroutine roessler_analytical_fp_rdp

end module LightKrylov_TestUtils
