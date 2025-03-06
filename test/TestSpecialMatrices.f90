module TestSpecialMatrices
    ! Fortran Standard Library.
    use iso_fortran_env
    use stdlib_math, only: is_close, all_close
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
                vec_out%data = matmul(self%data, vec_in%data)
            class default
                call stop_error("The intent [OUT] 'vec_out' must be of type 'vector_rdp'.", &
                                & this_module, "Poisson2D_matvec")
            end select
        class default
            call stop_error("The intent [IN] 'vec_in' must be of type 'vector_rdp'.", &
                            & this_module, "Poisson2D_matvec")
        end select
    end subroutine

    subroutine rmatvec_rdp(self, vec_in, vec_out)
        class(AbstractPoisson2D), intent(inout) :: self
        class(abstract_vector_rdp), intent(in) :: vec_in
        class(abstract_vector_rdp), intent(out) :: vec_out
        call matvec_rdp(self, vec_in, vec_out)
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
        type(vector_rdp), allocatable :: x, b
        ! CG options and metadata.
        type(cg_dp_opts) :: opts
        type(cg_dp_metadata) :: meta
        ! Information flag.
        integer :: info
        ! Misc.
        real(dp) :: err
        integer, parameter :: nx=16, ny=8
        character(len=256) :: msg

        ! Initialize linear problem.
        A = AbstractPoisson2D() ; A%data = Poisson2D(nx, ny)
        b = vector_rdp() ; call init_rand(b)
        x = vector_rdp() ; call x%zero()

        ! CG solver.
        opts = cg_dp_opts(maxiter=2*nx*ny, if_print_metadata=.true.)
        call cg(A, b, x, info, rtol=rtol_dp, atol=atol_dp, options=opts, meta=meta)
        call check_info(info, 'cg', module=this_module_long, procedure='test_pcg_poisson_rdp')

        ! Check convergence.
        err = norm2(matmul(A%data, x%data) - b%data)
        call get_err_str(msg, "max_err:", err)
        call check(error, err < b%norm()*rtol_dp)
        call check_test(error, 'test_pcg_poisson_rdp', eq='A @ x = b', context=msg)
    end subroutine

end module
