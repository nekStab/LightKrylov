program demo
   use stdlib_io_npy, only : save_npy
   use stdlib_logger, only : information_level, warning_level, debug_level, error_level, none_level
   use LightKrylov
   use LightKrylov, only: wp => dp
   use LightKrylov_Logger
   use LightKrylov_Utils
   use Roessler
   implicit none

   character(len=128), parameter :: this_module = 'Example Roessler'

   ! Roessler system.
   type(roessler_upo), allocatable :: sys
   ! State vectors
   type(state_vector), allocatable :: bf, dx, residual
   ! Integration time.
   real(wp), parameter :: Tend = 10.0_wp

   ! Misc
   integer :: i, info
   type(newton_dp_opts) :: opts
   type(gmres_dp_opts) :: gmres_opts
   real(wp) :: rnorm
   character(len=20) :: fmt
   
   write(fmt,*) '(A20,4(1X,F18.6))'

   ! Set up logging
   call logger_setup()
   call logger%configure(level=warning_level, time_stamp=.false.)

   ! Initialize baseflow and perturbation state vectors
   allocate(bf, dx, residual)
   call bf%zero()
   call dx%zero()
   call residual%zero()
   print '(A20,4(18X,A))', '         ', 'X', 'Y', 'Z', 'T'

   call set_position((/ 1.0_wp, 1.0_wp, 0.0_wp/), bf)  ! some initial guess
   bf%T = Tend ! period guess
   print fmt, 'Initial guess PO:', bf%x, bf%y, bf%z, bf%T
   print *,''

   ! Initialize system and Jacobian
   sys = roessler_upo()
   ! Set Jacobian and baseflow
   sys%jacobian = jacobian()
   sys%jacobian%X = bf

   opts       = newton_dp_opts(maxiter=30, ifbisect=.false., verbose=.true.)
   gmres_opts = gmres_dp_opts(atol=1e-12, verbose=.true.)
   call newton(sys, bf, info, opts, linear_solver_options=gmres_opts, scheduler=constant_atol_dp)

   call sys%eval(bf, residual, 0)
   print *,''
   print fmt, ' PO(0):  ', bf%x, bf%y, bf%z, bf%T
   print *, 'Compute residual of newton solution:'
   print fmt, ' res:    ', residual%x, residual%y, residual%z, residual%T
   print *,''

   call set_position((/ 1.0_wp, 1.0_wp, 0.0_wp/), bf)  ! some initial guess
   bf%T = Tend ! period guess
   sys%jacobian%X = bf

   call newton(sys, bf, info, opts, linear_solver_options=gmres_opts, scheduler=dynamic_tol_dp)

   call sys%eval(bf, residual, 0)
   print *,''
   print fmt, ' PO(0):  ', bf%x, bf%y, bf%z, bf%T
   print *, 'Compute residual of newton solution:'
   print fmt, ' res:    ', residual%x, residual%y, residual%z, residual%T
   print *,''

   return
end program demo
