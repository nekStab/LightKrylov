module TestSpecialMatrices
    ! Fortran Standard Library.
    use iso_fortran_env
    use stdlib_math, only: is_close, all_close
    use stdlib_io_npy, only: save_npy
    ! Testdrive
    use testdrive, only: new_unittest, unittest_type, error_type, check
    ! LightKrylov
    use LightKrylov
    use LightKrylov_Logger
    use LightKrylov_TestUtils
    ! SpecialMatrices
    use SpecialMatrices

    implicit none
    private

    character(len=*), parameter, private :: this_module      = "LK_TSpecialMatrices"
    character(len=*), parameter, private :: this_module_long = "LightKrylov_TestSpecialMatrices"
    integer, parameter :: nx = 16, ny = 8
    real(dp), parameter :: dx = 1.0_dp/(nx+1), dy = 1.0_dp/(ny+1)

    public :: collect_specialmatrices_rdp_testsuite

    !-----------------------------------------------
    !-----     ABSTRACT POISSON2D OPERATOR     -----
    !-----------------------------------------------

    type, extends(abstract_sym_linop_rdp), public :: AbstractPoisson2D
        type(Poisson2D), allocatable :: data
    contains
        private
        procedure, pass(self), public :: matvec => matvec_rdp
        procedure, pass(self), public :: rmatvec => rmatvec_rdp
    end type

    type, extends(abstract_precond_rdp) :: BlockJacobiPreconditioner
        type(SymTridiagonal), allocatable :: D
    contains
        private
        procedure, pass(self), public :: apply => apply_blockjacobi
    end type

contains

    !---------------------------------------------------------------
    !-----     TYPE-BOUND PROCEDURES FOR AbstractPoisson2D     -----
    !---------------------------------------------------------------

    subroutine matvec_rdp(self, vec_in, vec_out)
        class(AbstractPoisson2D), intent(inout) :: self
        class(abstract_vector_rdp), intent(in) :: vec_in
        class(abstract_vector_rdp), intent(out) :: vec_out

        select type(vec_in)
        type is(vector_rdp)
            select type(vec_out)
            type is(vector_rdp)
                vec_out%data = -matmul(self%data, vec_in%data)
            class default
                call type_error('vec_out','state_vector','OUT',this_module,'Poisson2D_matvec')
            end select
        class default
            call type_error('vec_in','state_vector','IN',this_module,'Poisson2D_matvec')
        end select
    end subroutine

    subroutine rmatvec_rdp(self, vec_in, vec_out)
        class(AbstractPoisson2D), intent(inout) :: self
        class(abstract_vector_rdp), intent(in) :: vec_in
        class(abstract_vector_rdp), intent(out) :: vec_out
        call matvec_rdp(self, vec_in, vec_out)
    end subroutine

    function construct_preconditioner(A) result(M)
        type(Poisson2D), intent(in) :: A
        type(BlockJacobiPreconditioner) :: M
        type(SymTridiagonal) :: D
        real(dp), allocatable :: dv(:), ev(:)
        integer :: i

        dv = [(2*(1.0_dp/dx**2 + 1.0_dp/dy**2), i=1, nx)]
        ev = [(-1.0_dp / dx**2, i=1, nx-1)]
        D = SymTridiagonal(dv, ev)
        M = BlockJacobiPreconditioner() ; M%D = D
    end function

    subroutine apply_blockjacobi(self, vec, iter, current_residual, target_residual)
        class(BlockJacobiPreconditioner), intent(inout) :: self
        class(abstract_vector_rdp), intent(inout) :: vec
        integer, optional, intent(in) :: iter
        real(dp), optional, intent(in) :: current_residual
        real(dp), optional, intent(in) :: target_residual
        
        real(dp), allocatable, target :: x(:), y(:)
        real(dp), pointer :: xmat(:, :), ymat(:, :)
        real(dp), allocatable :: z(:), Dmat(:, :)
        integer :: i, j

        select type(vec)
        type is(vector_rdp)
            x = vec%data ; allocate(y(nx*ny)) ; y = 0.0_dp
            xmat(1:nx, 1:ny) => x ; ymat(1:nx, 1:ny) => y
            do i = 1, ny
                ymat(:, i) = solve(self%D, xmat(:, i))
            enddo
            vec%data = y
        end select
    end subroutine

    !--------------------------------------------------------------------
    !-----     DEFINITION OF THE UNIT TESTS FOR SPECIALMATRICES     -----
    !--------------------------------------------------------------------

    subroutine collect_specialmatrices_rdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [ &
                    new_unittest("Precond. Conj. Gradient (Poisson2D)", test_pcg_poisson_rdp) &
                    ]
    end subroutine

    subroutine test_pcg_poisson_rdp(error)
        ! Error type.
        type(error_type), allocatable, intent(out) :: error
        ! Linear problem.
        type(Poisson2D) :: Amat
        type(AbstractPoisson2D) :: A
        type(BlockJacobiPreconditioner) :: P
        type(vector_rdp), allocatable :: x, b, r
        ! CG options and metadata.
        type(cg_dp_opts) :: opts
        type(cg_dp_metadata) :: meta
        ! Information flag.
        integer :: info
        ! Misc.
        real(dp) :: err
        character(len=256) :: msg

        ! Initialize linear problem.
        A = AbstractPoisson2D() ; A%data = Poisson2D(nx, ny)
        P = construct_preconditioner(A%data)
        b = vector_rdp() ; call init_rand(b)
        x = vector_rdp() ; call x%zero()
        r = vector_rdp() ; call r%zero()

        ! CG solver.
        opts = cg_dp_opts(maxiter=2*nx*ny, if_print_metadata=.true.)
        call cg(A, b, x, info, preconditioner=P, rtol=rtol_dp, atol=atol_dp, options=opts, meta=meta)
        call check_info(info, 'cg', module=this_module_long, procedure='test_pcg_poisson_rdp')

        ! Check convergence.
        call A%matvec(x, r) ; call r%sub(b) ; call P%apply(r)
        err = r%norm()
        call get_err_str(msg, "max_err:", err)
        call check(error, err < b%norm()*rtol_dp)
        call check_test(error, 'test_pcg_poisson_rdp', eq='A @ x = b', context=msg)
    end subroutine

end module
