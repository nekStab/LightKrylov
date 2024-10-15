program demo
   use stdlib_linalg, only : eye, eigvals
   use stdlib_io_npy, only : save_npy
   use stdlib_logger, only : information_level, warning_level, debug_level, error_level, none_level
   use LightKrylov
   use LightKrylov, only: wp => dp
   use LightKrylov_Logger
   use lightkrylov_IterativeSolvers, only: gmres_rdp
   use LightKrylov_Utils
   ! Roessler
   use Roessler
   use Roessler_OTD
   implicit none

   character(len=128), parameter :: this_module = 'Example Roessler'

   ! Roessler system.
   type(roessler_upo), allocatable :: sys
   ! State vectors
   type(state_vector), allocatable :: bf, dx, residual
   ! OTD basis
   type(pos_vector), allocatable :: OTD_in(:), OTD_out(:)
   ! Integration time.
   real(wp) :: Tend, TGS
   ! Report OTD?
   logical, parameter :: if_report_OTD    = .true.
   logical, parameter :: if_report_stdout = .true.

   ! Misc
   integer :: i, info, stat, iunit
   type(newton_dp_opts) :: opts
   type(gmres_dp_opts) :: gmres_opts
   real(wp) :: rnorm, tol
   character(len=20) :: fmt, file
   real(wp) :: M(npts,npts), Id(npts, npts)
   complex(wp) :: floquet_exponents(npts)
   logical :: exist
   
   write(fmt,*) '(A22,4(1X,F18.6))'
   file = 'roessler_OTD.txt'

   ! Set up logging
   call logger_setup()
   call logger%configure(level=error_level, time_stamp=.false.)

   ! Initialize baseflow and perturbation state vectors
   allocate(bf, dx, residual)
   call bf%zero()
   call dx%zero()
   call residual%zero()

   ! Set tolerace
   tol = 1e-12_wp

   print *, '########################################################################################'
   print '(A,E9.2,A)',' #             Newton iteration with constant tolerance (tol=', tol, ')                 #'
   print *, '########################################################################################'
   print *, ''

   call set_position((/ 0.0_wp, 6.1_wp, 1.3_wp /), bf)  ! some initial guess
   bf%T = 6.0_wp ! period guess
   print '(A22,4(16X,A,2X))', '         ', 'X', 'Y', 'Z', 'T'
   print fmt, 'Initial guess PO:   ', bf%x, bf%y, bf%z, bf%T
   print *,''

   ! Initialize system and Jacobian
   sys = roessler_upo()
   ! Set Jacobian and baseflow
   sys%jacobian = jacobian()
   sys%jacobian%X = bf

   opts = newton_dp_opts(maxiter=30, ifbisect=.false.)
   call newton(sys, bf, info, tolerance=tol, options=opts, linear_solver=gmres_rdp, scheduler=constant_atol_dp)

   call sys%eval(bf, residual, tol)
   print *, ''
   print '(A22,4(16X,A,2X))', '       ', 'X', 'Y', 'Z', 'T'
   print fmt, 'Solution:           ', bf%x, bf%y, bf%z, bf%T
   print fmt, 'Solution residual:  ', residual%x, residual%y, residual%z, residual%T
   print *,''

   print *, '########################################################################################'
   print '(A,E9.2,A)',' #             Newton iteration with dynamic tolerances (target=', tol, ')              #'
   print *, '########################################################################################'
   print *, ''

   call set_position((/ 0.0_wp, 6.1_wp, 1.3_wp /), bf)  ! some initial guess
   bf%T = 6.0_wp ! period guess
   print '(A22,4(16X,A,2X))', '         ', 'X', 'Y', 'Z', 'T'
   print fmt, 'Initial guess PO:  ', bf%x, bf%y, bf%z, bf%T
   print *,''
   sys%jacobian%X = bf
   
   call newton(sys, bf, info, tolerance=tol, options=opts, linear_solver=gmres_rdp, scheduler=dynamic_tol_dp)

   call sys%eval(bf, residual, tol)
   print *, ''
   print '(A22,4(16X,A,2X))', '       ', 'X', 'Y', 'Z', 'T'
   print fmt, 'Solution:           ', bf%x, bf%y, bf%z, bf%T
   print fmt, 'Solution residual:  ', residual%x, residual%y, residual%z, residual%T
   print *,''

   ! Compute the stability of the orbit
   sys%jacobian = floquet_operator()
   sys%jacobian%X = bf  ! <- periodic orbit

   M = 0.0_wp
   Id = eye(npts)
   do i = 1, npts
      call set_position(Id(:,i), dx)
      call sys%jacobian%matvec(dx, residual)
      call get_position(residual, M(:,i))
   end do
   floquet_exponents = eigvals(M)
   print *, 'Real part of the Floquet multipliers exp(T*mu) along the PO:'
   print *, ''
   do i = 1, npts
      print '(4X,I1,": ",E14.6)', i, real(floquet_exponents(i))
   end do
   print *, ''

   ! Compute OTD modes
   allocate(OTD_in(r), OTD_out(r))
   call zero_basis(OTD_in); call zero_basis(OTD_out)

   call rand_basis(OTD_in, ifnorm=.false.)
   call orthonormalize_basis(OTD_in)

   Tend = 10.0_wp
   if (if_report_OTD) then
      inquire(file=file, exist=exist)
      if (exist) open(unit=1234, file=file, status='old'); close(1234, status='delete')
      call write_header(file)
      call OTD_map(bf, OTD_in, Tend, OTD_out, if_report_stdout)
   else
      call OTD_map(bf, OTD_in, Tend, OTD_out)
   end if

end program demo
