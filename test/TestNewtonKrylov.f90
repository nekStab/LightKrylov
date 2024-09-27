module TestNewtonKrylov
    ! Fortran Standard library.
    use iso_fortran_env
    use stdlib_io_npy, only: save_npy
    use stdlib_math, only: is_close, all_close
    use stdlib_linalg, only: eye, diag
    use stdlib_stats, only : median
    ! Testdrive
    use testdrive, only: new_unittest, unittest_type, error_type, check
    ! LightKrylov
    use LightKrylov
    use LightKrylov_Constants
    use LightKrylov_Logger
    use LightKrylov_AbstractVectors
    use LightKrylov_NewtonKrylov
    ! Test Utilities
    use LightKrylov_TestUtils

    implicit none
    
    private

    character(len=128), parameter, private :: this_module = 'LightKrylov_TestNewtonKrylov'

    public :: collect_newton_rdp_testsuite

contains

    !----------------------------------------------------------
    !-----     DEFINITION OF THE UNIT TESTS FOR NEWTON    -----
    !----------------------------------------------------------

    subroutine collect_newton_rdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

         testsuite = [ &
                    new_unittest("Fixed point calculation", test_fixedp_rdp) &
                    ]
        return
    end subroutine collect_newton_rdp_testsuite

   subroutine test_fixedp_rdp(error)
       ! Error type.
       type(error_type), allocatable, intent(out) :: error
       ! Roessler system
       type(roessler), allocatable :: sys
       ! Jacobian
       type(jacobian), allocatable :: J
       ! Initial guess
       type(state_vector), allocatable :: X, Xtest
       ! Solver tolerance
       real(dp) :: tol = 1e-10_dp
       ! Verbosity
       logical :: verb = .true.
       ! Information flag.
       integer :: info
       ! Misc
       real(dp) :: err
       character(len=256) :: msg

       ! Allocate and set initial guess
       allocate(X);
       X%x = 0.0_dp
       X%y = 0.0_dp
       X%z = 0.0_dp

       ! Allocate and set roessler system
       sys = roessler()

       J = jacobian(X)
       sys%jacobian = J

       call newton(sys, X, tol, verb, info)

       ! check eigenvalues
       err = 0.0_dp !maxval(abs(eigvals - true_eigvals(:nev)))
       call get_err_str(msg, "max err: ", err)
       call check(error, err < rtol_dp)
       call check_test(error, 'test_fixedp_rdp', info='Fixed point value.', context=msg)

       return
   end subroutine test_fixedp_rdp

end module TestNewtonKrylov

