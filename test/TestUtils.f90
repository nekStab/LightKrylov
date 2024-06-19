module TestUtils
    use stdlib_io_npy, only: save_npy
    use stdlib_linalg, only: eye
    use stdlib_stats_distribution_normal, only: normal => rvs_normal
    use LightKrylov
    use LightKrylov_Constants
    use TestVectors
    use TestLinops
    
    implicit none
    
    private

    character(len=128), parameter, private :: this_module = 'LightKrylov_TestUtils'

    public :: get_data
    public :: put_data
    public :: init_rand

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

contains

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
        real(sp), allocatable :: mu(:, :), loc(:, :)
        allocate(mu(test_size, test_size)) ; mu = zero_rsp
        allocate(loc(test_size, test_size)) ; loc = one_rsp
        linop%data = normal(mu, loc)
        return
    end subroutine init_rand_linop_rsp

    subroutine init_rand_spd_linop_rsp(linop)
        type(spd_linop_rsp), intent(inout) :: linop
        real(sp), allocatable :: mu(:, :), loc(:, :)
        real(sp), allocatable :: data(:, :)
        allocate(mu(test_size, test_size)) ; mu = zero_rsp
        allocate(loc(test_size, test_size)) ; loc = one_rsp

        data = normal(mu, loc)
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
        real(dp), allocatable :: mu(:, :), loc(:, :)
        allocate(mu(test_size, test_size)) ; mu = zero_rdp
        allocate(loc(test_size, test_size)) ; loc = one_rdp
        linop%data = normal(mu, loc)
        return
    end subroutine init_rand_linop_rdp

    subroutine init_rand_spd_linop_rdp(linop)
        type(spd_linop_rdp), intent(inout) :: linop
        real(dp), allocatable :: mu(:, :), loc(:, :)
        real(dp), allocatable :: data(:, :)
        allocate(mu(test_size, test_size)) ; mu = zero_rdp
        allocate(loc(test_size, test_size)) ; loc = one_rdp

        data = normal(mu, loc)
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
        complex(sp), allocatable :: mu(:, :), loc(:, :)
        allocate(mu(test_size, test_size)) ; mu = zero_csp
        allocate(loc(test_size, test_size)) ; loc = one_csp
        linop%data = normal(mu, loc)
        return
    end subroutine init_rand_linop_csp

    subroutine init_rand_hermitian_linop_csp(linop)
        type(hermitian_linop_csp), intent(inout) :: linop
        complex(sp), allocatable :: data(:, :)
        complex(sp), allocatable :: mu(:, :), loc(:, :)

        allocate(mu(test_size, test_size)) ; mu = zero_csp
        allocate(loc(test_size, test_size)) ; loc = one_csp

        data = normal(loc, mu)
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
        complex(dp), allocatable :: mu(:, :), loc(:, :)
        allocate(mu(test_size, test_size)) ; mu = zero_cdp
        allocate(loc(test_size, test_size)) ; loc = one_cdp
        linop%data = normal(mu, loc)
        return
    end subroutine init_rand_linop_cdp

    subroutine init_rand_hermitian_linop_cdp(linop)
        type(hermitian_linop_cdp), intent(inout) :: linop
        complex(dp), allocatable :: data(:, :)
        complex(dp), allocatable :: mu(:, :), loc(:, :)

        allocate(mu(test_size, test_size)) ; mu = zero_cdp
        allocate(loc(test_size, test_size)) ; loc = one_cdp

        data = normal(loc, mu)
        data = matmul(data, transpose(conjg(data)))/test_size + 0.01*eye(test_size)
        linop%data = data

        return
    end subroutine init_rand_hermitian_linop_cdp


end module
