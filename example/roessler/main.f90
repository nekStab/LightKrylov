program demo
   use stdlib_io_npy, only : save_npy
   use LightKrylov
   use LightKrylov, only: wp => dp
   use LightKrylov_Logger
   use Roessler
   implicit none

   character(len=128), parameter :: this_module = 'Example Roessler'

   ! Roessler system.
   type(roessler_upo), allocatable :: sys
   ! Jacobian.
   type(jacobian), allocatable :: J
   ! State vectors
   type(state_vector), allocatable :: bf, dx
   ! Integration time.
   real(wp), parameter :: Tend = 10.0_wp
   ! Misc
   integer :: iunit, iostat, info
   type(newton_opts) :: opts
   type(gmres_dp_opts) :: gmres_opts

   ! NEWTON
   type(state_vector), allocatable :: residual, xtest
   real(wp) :: rnorm
   integer :: i
   character(len=128) :: fmt

   write(fmt,*) '(A,4(1X,F18.6))'

   ! Set up logging
   call logger_setup()

   ! Initialize system
   sys = roessler_upo() !tau=Tend)

   ! Initialize baseflow and perturbation state vectors
   allocate(bf, dx, residual, xtest)
   call bf%zero(); call dx%zero(); call residual%zero(); call xtest%zero()

   call set_position((/ 1.0_wp, 1.0_wp, 0.0_wp/), bf)  ! some initial guess
   bf%T = Tend ! period guess

   J = jacobian() !tau=Tend)
   sys%jacobian = J
   sys%jacobian%X = bf

   ! compute residual
   call sys%eval(bf, residual)
   call residual%chsgn()
   print fmt, '-res:', residual%x, residual%y, residual%z, residual%T
   
   call gmres(sys%jacobian, residual, dx, info)

   print fmt, 'dX  :', dx%x, dx%y, dx%z, dx%T

   call sys%jacobian%matvec(dx, xtest)
   call xtest%add(residual)
   print fmt, 'J@dx + res:', xtest%x, xtest%y, xtest%z, xtest%T

   STOP 4

   opts       = newton_opts(maxiter=10, ifbisect=.false., verbose=.true.)
   gmres_opts = gmres_dp_opts(atol=1e-6, verbose=.true.)
   call newton(sys, bf, info, opts, linear_solver_options=gmres_opts)

   return
end program demo
