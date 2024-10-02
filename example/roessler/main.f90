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
   !type(jacobian), allocatable :: J
   ! State vectors
   type(state_vector), allocatable :: Xin, Xout
   !> Integration time.
   real(wp), parameter :: Tend = 1.0_wp

   !> Set up logging
   call logger_setup()

   !> Initialize system
   sys = roessler_upo(tau=Tend)

   !> Initialize state vector
   allocate(Xin, Xout)
   call Xin%zero()
   call Xout%zero()

   call sys%eval(Xin, Xout)

   print *, Xout%x, Xout%y, Xout%z

end program demo
