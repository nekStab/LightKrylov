module TestNewtonKrylov
    ! Fortran Standard library.
    use iso_fortran_env
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
       type(state_vector), allocatable :: X, fp1, fp2
       ! Solver tolerance
       real(dp) :: tol = 1e-10_dp
       ! Verbosity
       logical :: verb = .false.
       ! Information flag.
       integer :: info, i
       ! Misc
       real(dp) :: err
       character(len=256) :: msg, infomsg

       ! Allocate and set initial guess
       allocate(X, fp1, fp2);

       ! Allocate and set Roessler system and Jacobian
       sys = roessler()
       J   = jacobian()
       sys%jacobian = J
       
       call roessler_analytical_fp(fp1, fp2)

       X%x = 0.0_dp
       X%y = 0.0_dp
       X%z = 0.0_dp
       call newton(sys, X, tol, verb, info)
       call X%sub(fp1)

       ! check fixed point 1
       write(infomsg, '(A1,E8.2,A1,E9.2,A1,E8.2,A1)') '(',fp1%x,',',fp1%y,',',fp1%z,')'
       !write(infomsg, *) '|| X_newton - fp1 ||_2'
       err = X%norm()
       call get_err_str(msg, "max err: ", err)
       call check(error, err < rtol_dp)
       call check_test(error, 'test_fixedp_rdp', info=infomsg, context=msg)
       
       X%x = 10.0_dp
       X%y = -5.0_dp
       X%z = 20.0_dp
       call newton(sys, X, tol, verb, info)
       call X%sub(fp2)

       ! check fixed point 2
       write(infomsg, '(A1,E8.2,A1,E9.2,A1,E8.2,A1)') '(',fp2%x,',',fp2%y,',',fp2%z,')'
       !write(infomsg, *) '|| X_newton - fp2 ||_2'
       err = X%norm()
       call get_err_str(msg, "max err: ", err)
       call check(error, err < rtol_dp)
       call check_test(error, 'test_fixedp_rdp', info=infomsg, context=msg)
       
       return
   end subroutine test_fixedp_rdp

end module TestNewtonKrylov

