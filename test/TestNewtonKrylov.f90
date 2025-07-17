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
    use LightKrylov_Utils
    ! Test Utilities
    use LightKrylov_TestUtils

    implicit none
    
    private

    character(len=*), parameter, private :: this_module      = 'LK_TNwtKryl'
    character(len=*), parameter, private :: this_module_long = 'LightKrylov_TestNewtonKrylov'

    public :: collect_newton_rsp_testsuite
    public :: collect_newton_rdp_testsuite
    public :: collect_newton_csp_testsuite
    public :: collect_newton_cdp_testsuite

contains

    !----------------------------------------------------------
    !-----     DEFINITION OF THE UNIT TESTS FOR NEWTON    -----
    !----------------------------------------------------------

    subroutine collect_newton_rsp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

         testsuite = [ &
                    new_unittest("Fixed point calculation", test_fixedp_rsp) &
                    ]
        return
    end subroutine collect_newton_rsp_testsuite

    subroutine test_fixedp_rsp(error)
       ! Error type.
       type(error_type), allocatable, intent(out) :: error
       ! Roessler system
       type(roessler_rsp), allocatable :: sys
       ! Initial guess
       type(state_vector_rsp), allocatable :: X, fp1, fp2
       ! Newton options
       type(newton_sp_opts) :: opts
       ! Information flag.
       integer :: info
       ! Misc
       real(sp) :: err
       character(len=256) :: msg, infomsg

       ! Allocate solution variables and reference values
       allocate(X, fp1, fp2);

       ! Newton opts
       opts = newton_sp_opts(maxiter=60, ifbisect=.false., if_print_metadata=.true.) 

       ! Allocate and set Roessler system and Jacobian
       sys = roessler_rsp()
       sys%jacobian = jacobian_rsp()
       
       call roessler_analytical_fp_rsp(fp1, fp2)

       X%x = zero_rsp
       X%y = zero_rsp
       X%z = zero_rsp
       call newton(sys, X, gmres_rsp, info, rtol=10*atol_sp, atol=10*atol_sp, options=opts)
       call X%sub(fp1)

       ! check fixed point 1
       !write(infomsg, *) '|| X_newton - fp1 ||_2'
       err = X%norm()
       call get_err_str(msg, "max err: ", err)
       call check(error, err < rtol_sp)
       write(infomsg, '(A1,E8.2,A1,E9.2,A1,E8.2,A1)') '(',fp1%x,',',fp1%y,',',fp1%z,')'
       call check_test(error, 'test_fixedp_rsp', info=infomsg, context=msg)

       X%x = zero_rsp
       X%y = zero_rsp
       X%z = zero_rsp
       opts%ifbisect = .true.
       call newton(sys, X, gmres_rsp, info, rtol=10*atol_sp, atol=10*atol_sp, options=opts)
       call X%sub(fp1)

       ! check fixed point 1 with bisection (if necessary)
       err = X%norm()
       call get_err_str(msg, "max err: ", err)
       call check(error, err < rtol_sp)
       call check_test(error, 'test_fixedp_rsp', info='Newton with step bisection', context=msg)

       return
   end subroutine test_fixedp_rsp
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
       type(roessler_rdp), allocatable :: sys
       ! Initial guess
       type(state_vector_rdp), allocatable :: X, fp1, fp2
       ! Newton options
       type(newton_dp_opts) :: opts
       ! Information flag.
       integer :: info
       ! Misc
       real(dp) :: err
       character(len=256) :: msg, infomsg

       ! Allocate solution variables and reference values
       allocate(X, fp1, fp2);

       ! Newton opts
       opts = newton_dp_opts(maxiter=60, ifbisect=.false., if_print_metadata=.true.) 

       ! Allocate and set Roessler system and Jacobian
       sys = roessler_rdp()
       sys%jacobian = jacobian_rdp()
       
       call roessler_analytical_fp_rdp(fp1, fp2)

       X%x = zero_rdp
       X%y = zero_rdp
       X%z = zero_rdp
       call newton(sys, X, gmres_rdp, info, rtol=10*atol_dp, atol=10*atol_dp, options=opts)
       call X%sub(fp1)

       ! check fixed point 1
       !write(infomsg, *) '|| X_newton - fp1 ||_2'
       err = X%norm()
       call get_err_str(msg, "max err: ", err)
       call check(error, err < rtol_dp)
       write(infomsg, '(A1,E8.2,A1,E9.2,A1,E8.2,A1)') '(',fp1%x,',',fp1%y,',',fp1%z,')'
       call check_test(error, 'test_fixedp_rdp', info=infomsg, context=msg)

       X%x = zero_rdp
       X%y = zero_rdp
       X%z = zero_rdp
       opts%ifbisect = .true.
       call newton(sys, X, gmres_rdp, info, rtol=10*atol_dp, atol=10*atol_dp, options=opts)
       call X%sub(fp1)

       ! check fixed point 1 with bisection (if necessary)
       err = X%norm()
       call get_err_str(msg, "max err: ", err)
       call check(error, err < rtol_dp)
       call check_test(error, 'test_fixedp_rdp', info='Newton with step bisection', context=msg)

       return
   end subroutine test_fixedp_rdp
    subroutine collect_newton_csp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

         testsuite = [ &
                    new_unittest("Fixed point calculation", test_fixedp_csp) &
                    ]
        return
    end subroutine collect_newton_csp_testsuite

    subroutine test_fixedp_csp(error)
       ! Error type.
       type(error_type), allocatable, intent(out) :: error
       ! Roessler system
       type(roessler_csp), allocatable :: sys
       ! Initial guess
       type(state_vector_csp), allocatable :: X, fp1, fp2
       ! Newton options
       type(newton_sp_opts) :: opts
       ! Information flag.
       integer :: info
       ! Misc
       real(sp) :: err
       character(len=256) :: msg, infomsg

       ! Allocate solution variables and reference values
       allocate(X, fp1, fp2);

       ! Newton opts
       opts = newton_sp_opts(maxiter=60, ifbisect=.false., if_print_metadata=.true.) 

       ! Allocate and set Roessler system and Jacobian
       sys = roessler_csp()
       sys%jacobian = jacobian_csp()
       
       call roessler_analytical_fp_csp(fp1, fp2)

       X%x = zero_csp
       X%y = zero_csp
       X%z = zero_csp
       call newton(sys, X, gmres_csp, info, rtol=10*atol_sp, atol=10*atol_sp, options=opts)
       call X%sub(fp1)

       ! check fixed point 1
       !write(infomsg, *) '|| X_newton - fp1 ||_2'
       err = X%norm()
       call get_err_str(msg, "max err: ", err)
       call check(error, err < rtol_sp)
       write(infomsg, '(A1,E8.2,A1,E9.2,A1,E8.2,A1)') '(',fp1%x%re,',',fp1%y%re,',',fp1%z%re,')'
       call check_test(error, 'test_fixedp_csp', info=infomsg, context=msg)
       write(infomsg, '(A1,E8.2,A1,E9.2,A1,E8.2,A1)') '(',fp1%x%im,',',fp1%y%im,',',fp1%z%im,')'
       call check_test(error, 'test_fixedp_csp', info=infomsg, context=msg)

       X%x = zero_csp
       X%y = zero_csp
       X%z = zero_csp
       opts%ifbisect = .true.
       call newton(sys, X, gmres_csp, info, rtol=10*atol_sp, atol=10*atol_sp, options=opts)
       call X%sub(fp1)

       ! check fixed point 1 with bisection (if necessary)
       err = X%norm()
       call get_err_str(msg, "max err: ", err)
       call check(error, err < rtol_sp)
       call check_test(error, 'test_fixedp_csp', info='Newton with step bisection', context=msg)

       return
   end subroutine test_fixedp_csp
    subroutine collect_newton_cdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

         testsuite = [ &
                    new_unittest("Fixed point calculation", test_fixedp_cdp) &
                    ]
        return
    end subroutine collect_newton_cdp_testsuite

    subroutine test_fixedp_cdp(error)
       ! Error type.
       type(error_type), allocatable, intent(out) :: error
       ! Roessler system
       type(roessler_cdp), allocatable :: sys
       ! Initial guess
       type(state_vector_cdp), allocatable :: X, fp1, fp2
       ! Newton options
       type(newton_dp_opts) :: opts
       ! Information flag.
       integer :: info
       ! Misc
       real(dp) :: err
       character(len=256) :: msg, infomsg

       ! Allocate solution variables and reference values
       allocate(X, fp1, fp2);

       ! Newton opts
       opts = newton_dp_opts(maxiter=60, ifbisect=.false., if_print_metadata=.true.) 

       ! Allocate and set Roessler system and Jacobian
       sys = roessler_cdp()
       sys%jacobian = jacobian_cdp()
       
       call roessler_analytical_fp_cdp(fp1, fp2)

       X%x = zero_cdp
       X%y = zero_cdp
       X%z = zero_cdp
       call newton(sys, X, gmres_cdp, info, rtol=10*atol_dp, atol=10*atol_dp, options=opts)
       call X%sub(fp1)

       ! check fixed point 1
       !write(infomsg, *) '|| X_newton - fp1 ||_2'
       err = X%norm()
       call get_err_str(msg, "max err: ", err)
       call check(error, err < rtol_dp)
       write(infomsg, '(A1,E8.2,A1,E9.2,A1,E8.2,A1)') '(',fp1%x%re,',',fp1%y%re,',',fp1%z%re,')'
       call check_test(error, 'test_fixedp_cdp', info=infomsg, context=msg)
       write(infomsg, '(A1,E8.2,A1,E9.2,A1,E8.2,A1)') '(',fp1%x%im,',',fp1%y%im,',',fp1%z%im,')'
       call check_test(error, 'test_fixedp_cdp', info=infomsg, context=msg)

       X%x = zero_cdp
       X%y = zero_cdp
       X%z = zero_cdp
       opts%ifbisect = .true.
       call newton(sys, X, gmres_cdp, info, rtol=10*atol_dp, atol=10*atol_dp, options=opts)
       call X%sub(fp1)

       ! check fixed point 1 with bisection (if necessary)
       err = X%norm()
       call get_err_str(msg, "max err: ", err)
       call check(error, err < rtol_dp)
       call check_test(error, 'test_fixedp_cdp', info='Newton with step bisection', context=msg)

       return
   end subroutine test_fixedp_cdp

end module TestNewtonKrylov

