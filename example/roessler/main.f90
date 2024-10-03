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
   type(state_vector), allocatable :: bfin, bfout, xpin, xpout
   ! Integration time.
   real(wp), parameter :: Tend = 1.0_wp
   ! Misc
   integer :: iunit, iostat, info
   type(newton_opts) :: opts
   type(gmres_dp_opts) :: gmres_opts

   ! NEWTON
   type(state_vector), allocatable :: residual
   real(wp) :: rnorm
   integer :: i
   character(len=128) :: fmt

   write(fmt,*) '(A,4(1X,F18.6))'

   ! Set up logging
   call logger_setup()

   ! Initialize system
   sys = roessler_upo(tau=Tend)

   ! Initialize baseflow and perturbation state vectors
   allocate(xpin, xpout); call xpin%zero(); call xpout%zero()
   allocate(bfin, bfout, residual); call bfin%zero(); call bfout%zero(); call residual%zero()

   call set_position((/ 1.0_wp, 1.0_wp, 0.0_wp/), bfin)
   !call set_position((/ 0.1_wp, 6.0917_wp, 1.2997_wp/), bfin)
   call set_position((/ -1.0_wp, 0.0_wp, 0.0_wp/), xpin)
   bfin%T = Tend

   !open(newunit=iunit, file='new_orbit.txt', status='new', action='write')
   !close(iunit)

   call sys%eval(bfin, bfout)

   !call rename('new_orbit.txt', 'BF_orbit.txt', iostat)

   J = jacobian(tau=Tend)
   sys%jacobian = J
   sys%jacobian%X = bfin
   
   !call sys%jacobian%matvec(xpin, xpout)

   !! NEWTON

   !print fmt, 'BFIN:', bfin%x, bfin%y, bfin%z, bfin%T

   i = 1
   !print fmt, 'RES: ', residual%x, residual%y, residual%z, residual%T
   !rnorm = residual%norm()
   !write(*,*) "Iteration", i, ": Residual norm = ", rnorm




   !call rename('new_orbit.txt', 'BF_orbit.txt', iostat)




   !open(newunit=iunit, file='new_orbit.txt', status='new', action='write')
   !close(iunit)

   ! Define the Jacobian
   sys%jacobian%X = bfin

   !print fmt, 'XPIN:',xpin%x, xpin%y, xpin%z, xpin%T
   call sys%jacobian%matvec(xpin, xpout)

   !call set_position((/ 1.0_wp, 0.0_wp, 0.0_wp/), xpin)
   !call sys%jacobian%matvec(xpin, xpout)
   !call set_position((/ 0.0_wp, 1.0_wp, 0.0_wp/), xpin)
   !call sys%jacobian%matvec(xpin, xpout)
   !call set_position((/ 0.0_wp, 0.0_wp, 1.0_wp/), xpin)
   !call sys%jacobian%matvec(xpin, xpout)

   !call rename('new_orbit.txt', 'PT_orbit.txt', iostat)
   !STOP 6
   print *, 'compute residual:'
   call sys%eval(bfin, residual)
   call residual%chsgn()
   print fmt, '-RES:',residual%x, residual%y, residual%z, residual%T
   !call sys%jacobian%matvec(xpin, xpout)
   !print fmt, 'J@-r:', xpout%x, xpout%y, xpout%z, xpout%T

   !call sys%jacobian%matvec(residual, xpout)
   !print fmt, 'J@-r:',xpout%x, xpout%y, xpout%z, xpout%T
   !call sys%jacobian%matvec(xpout, residual)
   !print fmt, 'J^2@-r:', residual%x, residual%y, residual%z, residual%T
   !call sys%jacobian%matvec(residual, xpout)
   !print fmt, 'J^3@-r:',xpout%x, xpout%y, xpout%z, xpout%T

   call gmres(sys%jacobian, residual, xpout, info)

   print fmt, 'dX  :',xpout%x, xpout%y, xpout%z, xpout%T

   call sys%jacobian%matvec(xpout, xpin)
   call xpin%add(residual)
   print fmt, 'J@dx + r:',xpin%x, xpin%y, xpin%z, xpin%T

   STOP 5

   opts = newton_opts(maxiter=10, ifbisect=.false., verbose=.true.)
   call newton(sys, bfin, info, opts)

   return
end program demo
