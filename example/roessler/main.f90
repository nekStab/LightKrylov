program demo
   use stdlib_linalg, only : eye, eigvals
   use stdlib_io_npy, only : save_npy
   use stdlib_sorting, only : sort
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
   type(state_vector), allocatable :: bf, dx, residual, fp
   ! Position vectors
   type(pos_vector), allocatable   :: bfp
   ! OTD basis
   type(pos_vector), allocatable   :: OTD_in(:), OTD_out(:)

   ! Misc
   type(newton_dp_opts)            :: opts
   type(gmres_dp_opts)             :: gmres_opts
   integer                         :: i, j, info
   real(wp)                        :: rnorm, tol, Tend, t_FTLE, d
   real(wp), dimension(npts, npts) :: M, Id
   real(wp), dimension(npts)       :: eval, vec
   real(wp), dimension(npts, r)    :: u, Lu
   real(wp), dimension(r, r)       :: Lr
   ! IO
   character(len=20)    :: fmt
   
   write(fmt,*) '(A22,4(1X,F18.6))'

   ! Set up logging
   call logger_setup()
   call logger%configure(level=error_level, time_stamp=.false.)

   ! Initialize baseflow and perturbation state vectors
   allocate(bf, dx, residual)
   call bf%zero(); call dx%zero(); call residual%zero()

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

   print *, '########################################################################################'
   print *, '#                        Monodromy matrix and floquet exponents                        #'
   print *, '########################################################################################'
   print *, ''

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
   eval = real(eigvals(M))
   call sort(eval, reverse=.true.)
   print *, 'Real part of the Floquet multipliers exp(T*mu) along the PO:'
   print *, ''
   do i = 1, npts
      print '(4X,I1,": ",F15.12)', i, eval(i)
   end do
   print *, ''

   print *, '########################################################################################'
   print *, '#                  Optimally Time-Dependent (OTD) modes on fixed point                 #'
   print *, '########################################################################################'
   print *, ''
   allocate(bfp); call bfp%zero()
   ! Set the baseflow to a fixed point
   d = sqrt(c**2 - 4*a*b)
   bfp%x = ( c - d)/ 2
   bfp%y = (-c + d)/(2*a)
   bfp%z = ( c - d)/(2*a)

   ! Compute OTD modes on the fixed point
   allocate(OTD_in(r), OTD_out(r))
   call zero_basis(OTD_in); call zero_basis(OTD_out)

   ! Initialize basis
   call rand_basis(OTD_in, ifnorm=.false.)
   call orthonormalize_basis(OTD_in)

   ! We need long enough to converge to the invariant tangent space
   Tend = 5.0_wp
   t_FTLE = 5.0_wp
   call write_header()
   call OTD_map(bfp, OTD_in, Tend, OTD_out, t_FTLE)
   call rename(file, 'example/roessler/FP_OTD.txt')
   ! get baseflow
   call get_pos(bfp, vec)
   ! get OTD basis vectors
   u  = 0.0_wp
   Lu = 0.0_wp
   do i = 1, r
      call get_pos(OTD_out(i), u(:,i))
      call linear_roessler(u(:,i), vec, Lu(:,i))
   end do
   ! compute Lr
   Lr = 0.0_wp
   do i = 1, r
      do j = 1, r
         Lr(i,j) = dot_product(Lu(:,i), u(:,j))
      end do
   end do
   eval = 0.0_wp
   eval(1:r) = eigvals(Lr)
   print '(*(A16,1X))', ' ', 'lambda_1', 'lambda_2'
   print '(A16,1X,*(F16.9,1X))', 'Reference   ', EV_ref
   print *, ''
   print '(A10,F6.3,1X,*(F16.9,1X))', 'OTD:  t=', Tend, eval(1:r)
   print *, ''
   print *, '########################################################################################'
   print *, '#                  Optimally Time-Dependent (OTD) modes on periodic orbit              #'
   print *, '########################################################################################'
   print *, ''
   ! Now move to the periodic orbit
   call get_position(bf, vec)
   call set_pos(vec, bfp)

   ! Reinitialize basis
   call rand_basis(OTD_in, ifnorm=.false.)
   call orthonormalize_basis(OTD_in)

   ! We need long enough to converge to the invariant periodic tangent space
   Tend = 30.0_wp*bf%T
   t_FTLE = bf%T
   call write_header(); call write_header_LE()
   print '(*(A16,1X))', ' ', 'LE_1', 'LE_2'
   print '(A16,1X,2(F16.9,1X),A16,1X,F16.9)', 'Reference   ', LE_ref, 'Period T=', bf%T
   call OTD_map(bfp, OTD_in, Tend, OTD_out, t_FTLE, if_rst=.true.)
   call rename(file,    'example/roessler/PO_OTD.txt')
   call rename(file_LE, 'example/roessler/PO_LE.txt')
   print *, ''

   print *, ''
   print *, '########################################################################################'
   print *, '#                  Optimally Time-Dependent (OTD) modes on Route to Chaos              #'
   print *, '########################################################################################'
   print *, ''

   ! We use the old converged basis
   call orthonormalize_basis(OTD_in)

   ! We need long enough for the orbit to return to the chaotic attractor
   Tend = 60.0_wp*bf%T
   t_FTLE = bf%T
   call write_header(); call write_header_LE()
   print '(*(A16,1X))', ' ', 'FTLE_1', 'FTLE_2'
   print '(A16,1X,*(F16.9,1X))', 'Reference   ', LE_ref
   print *, ''
   call OTD_map(bfp, OTD_in, Tend, OTD_out, t_FTLE, if_rst=.false.) ! we do not reset the bf!
   call rename(file,    'example/roessler/PO-chaos_OTD.txt')
   call rename(file_LE, 'example/roessler/PO-chaos_LE.txt')
   print *, ''
   
end program demo