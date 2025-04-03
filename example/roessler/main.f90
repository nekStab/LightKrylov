program demo
   use stdlib_linalg, only: eye, eigvals
   use stdlib_io_npy, only: save_npy
   use stdlib_sorting, only: sort
   use stdlib_strings, only: padl
   use stdlib_logger, only: information_level, warning_level, debug_level, error_level, none_level
   ! RKLIB module for time integration.
   use rklib_module
   ! LightKrylov for linear algebra
   use LightKrylov
   use LightKrylov_Constants
   use LightKrylov_Logger
   use LightKrylov_Timing, only: timer => global_lightkrylov_timer
   use LightKrylov_Utils
   use lightkrylov_IterativeSolvers, only: gmres_rdp
   ! Roessler system
   use Roessler
   use Roessler_OTD
   implicit none

   ! character(len=*), parameter :: this_module = 'Example Roessler'

   ! Roessler system.
   type(roessler_upo) :: sys_jac, sys_floquet
   ! State vectors
   type(state_vector) :: bf, dx, residual
   ! Position vectors
   type(pos_vector)   :: bfp
   ! OTD basis
   type(pos_vector)   :: OTD_in(r), OTD_out(r)
   ! Time-integrator
   type(rks54_class)  :: nonlinear_integrator

   ! Misc
   type(newton_dp_opts)            :: opts
   integer                         :: i, j, info
   real(dp)                        :: tol, Tend, t_FTLE, d
   real(dp), dimension(npts, npts) :: M, Id
   real(dp), dimension(npts)       :: eval, vec
   real(dp), dimension(npts, r)    :: u, Lu
   real(dp), dimension(r, r)       :: Lr

   ! IO
   character(len=20)    :: data_fmt, header_fmt

   write (header_fmt, *) '(22X,*(A,2X))'
   write (data_fmt, *) '(A22,*(1X,F15.6))'

   ! Set up logging
   call logger_setup()
   call logger%configure(level=error_level, time_stamp=.false.)

   ! Set up timing
   call timer%initialize()
   call timer%add_timer('Roessler example (total)', start=.true.)
   call timer%add_timer('Chaotic attractor', start=.true.)

   ! Initialize baseflow and perturbation state vectors
   call bf%zero(); call dx%zero(); call residual%zero()

   print *, '########################################################################################'
   print *, '#                         Roessler system chaotic attractor                            #'
   print *, '########################################################################################'
   print *, ''

   vec = (/0.0_dp, -5.0_dp, 0.05_dp/) ! some intial point
   Tend = 300.0_dp ! Integration time
   ! Integrate equations
   call write_report_header
   call nonlinear_integrator%initialize(n=npts, f=NL_rhs, report=roessler_report_file, report_rate=20)
   call nonlinear_integrator%integrate(0.0_dp, vec, 1.0_dp, Tend, eval)
   call rename(report_file, 'example/roessler/roessler_attractor.txt')

   print header_fmt, padl('X', 14), padl('Y', 14), padl('Z', 14), padl('time', 14)
   print data_fmt, 'Initial position :', vec(1), vec(2), vec(3), 0.0_dp
   print data_fmt, 'Final position :', eval(1), eval(2), eval(3), Tend
   print *, ''

   call timer%stop('Chaotic attractor')
   call timer%add_timer('Newton iteration (const. tol)', start=.true.)

   print *, '########################################################################################'
   print '(A,E9.2,A)', ' #             Newton iteration with constant tolerance (tol=', tol, ')                 #'
   print *, '########################################################################################'
   print *, ''

   call set_position((/0.0_dp, 6.1_dp, 1.3_dp/), bf)  ! initial guess
   bf%T = 6.0_dp ! period guess
   print header_fmt, padl('X', 14), padl('Y', 14), padl('Z', 14), padl('T', 14)
   print data_fmt, 'Initial guess PO:', bf%x, bf%y, bf%z, bf%T
   print *, ''

   ! Initialize system and Jacobian
   sys_jac = roessler_upo()
   ! Set Jacobian and baseflow
   sys_jac%jacobian = jacobian()
   sys_jac%jacobian%X = bf

   ! Set tolerance
   tol = 1e-12_dp

   opts = newton_dp_opts(maxiter=30, ifbisect=.false.)
   call newton(sys_jac, bf, gmres_rdp, info, rtol=tol, atol=tol, options=opts, scheduler=constant_tol_dp)

   call sys_jac%eval(bf, residual, tol)
   print *, ''
   print header_fmt, padl('X', 14), padl('Y', 14), padl('Z', 14), padl('T', 14)
   print data_fmt, 'Solution:         ', bf%x, bf%y, bf%z, bf%T
   print data_fmt, 'Solution residual:', residual%x, residual%y, residual%z, residual%T
   print *, ''

   ! Reset timers
   call timer%stop('Newton iteration (const. tol)')
   call sys_jac%reset_timer()
   call sys_jac%jacobian%reset_timer()
   call sys_jac%jacobian%finalize_timer()
   call timer%reset_all()
   call timer%add_timer('Newton iteration (dyn. tol)', start=.true.)

   print *, '########################################################################################'
   print '(A,E9.2,A)', ' #             Newton iteration with dynamic tolerances (target=', tol, ')              #'
   print *, '########################################################################################'
   print *, ''

   call set_position((/0.0_dp, 6.1_dp, 1.3_dp/), bf)  ! some initial guess
   bf%T = 6.0_dp ! period guess
   print header_fmt, padl('X', 14), padl('Y', 14), padl('Z', 14), padl('T', 14)
   print data_fmt, 'Initial guess PO:  ', bf%x, bf%y, bf%z, bf%T
   print *, ''
   sys_jac%jacobian%X = bf

   call newton(sys_jac, bf, gmres_rdp, info, rtol=tol, atol=tol, options=opts, scheduler=dynamic_tol_dp)

   call sys_jac%eval(bf, residual, tol)
   print *, ''
   print header_fmt, padl('X', 14), padl('Y', 14), padl('Z', 14), padl('T', 14)
   print data_fmt, 'Solution:         ', bf%x, bf%y, bf%z, bf%T
   print data_fmt, 'Solution residual:', residual%x, residual%y, residual%z, residual%T
   print *, ''

   ! Reset timers
   call timer%stop('Newton iteration (dyn. tol)')
   call sys_jac%reset_timer()
   call sys_jac%jacobian%reset_timer()
   call timer%reset_all()
   call timer%add_timer('Monodromy matrix & Floquet exp.', start=.true.)

   print *, '########################################################################################'
   print *, '#                        Monodromy matrix and floquet exponents                        #'
   print *, '########################################################################################'
   print *, ''

   ! Initialize system and Jacobian
   sys_floquet = roessler_upo()
   ! Compute the stability of the orbit
   sys_floquet%jacobian = floquet_operator()
   sys_floquet%jacobian%X = bf  ! <- periodic orbit

   M = 0.0_dp
   Id = eye(npts)
   do i = 1, npts
      call set_position(Id(:, i), dx)
      call sys_floquet%jacobian%apply_matvec(dx, residual)
      call get_position(residual, M(:, i))
   end do
   eval = real(eigvals(M))
   call sort(eval, reverse=.true.)
   print *, 'Real part of the Floquet multipliers exp(T*mu) along the PO:'
   print *, ''
   do i = 1, npts
      print '(4X,I1,": ",F15.12)', i, eval(i)
   end do
   print *, ''

   ! Reset timers
   call timer%stop('Monodromy matrix & Floquet exp.')
   call sys_floquet%reset_timer()
   call sys_floquet%jacobian%reset_timer()
   call timer%reset_all()
   call timer%add_timer('OTD modes - fixed-point', start=.true.)

   print *, '########################################################################################'
   print *, '#                  Optimally Time-Dependent (OTD) modes on fixed point                 #'
   print *, '########################################################################################'
   print *, ''
   ! Initialize fixed point
   call bfp%zero()
   ! Set the baseflow to a fixed point
   d = sqrt(c**2 - 4*a*b)
   bfp%x = (c - d)/2
   bfp%y = (-c + d)/(2*a)
   bfp%z = (c - d)/(2*a)

   ! Compute OTD modes on the fixed point
   call zero_basis(OTD_in); call zero_basis(OTD_out)

   ! Initialize basis
   call rand_basis(OTD_in, ifnorm=.false.)
   call orthonormalize_basis(OTD_in)

   ! We need long enough to converge to the invariant tangent space
   Tend = 5.0_dp
   t_FTLE = 5.0_dp
   call write_header()
   call OTD_map(bfp, OTD_in, Tend, OTD_out, t_FTLE)
   call rename(report_file_OTD, 'example/roessler/FP_OTD.txt')
   ! get baseflow
   call get_pos(bfp, vec)
   ! get OTD basis vectors
   u = 0.0_dp
   Lu = 0.0_dp
   do i = 1, r
      call get_pos(OTD_out(i), u(:, i))
      call linear_roessler(u(:, i), vec, Lu(:, i))
   end do
   ! compute Lr
   Lr = 0.0_dp
   do i = 1, r
      do j = 1, r
         Lr(i, j) = dot_product(Lu(:, i), u(:, j))
      end do
   end do
   eval = 0.0_dp
   eval(1:r) = real(eigvals(Lr), kind=dp)
   print '(*(A16,1X))', ' ', 'lambda_1', 'lambda_2'
   print '(A16,1X,*(F16.9,1X))', 'Reference   ', EV_ref
   print *, ''
   print '(A10,F6.3,1X,*(F16.9,1X))', 'OTD:  t=', Tend, eval(1:r)

   ! Reset timers
   call timer%stop('OTD modes - fixed-point')
   call sys_floquet%reset_timer()
   call sys_floquet%jacobian%reset_timer()
   call timer%reset_all()
   call timer%add_timer('OTD modes - periodic orbit', start=.true.)

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
   Tend = 30.0_dp*bf%T
   t_FTLE = bf%T
   call write_header(); call write_header_LE()
   print '(*(A16,1X))', ' ', 'LE_1', 'LE_2'
   print '(A16,1X,2(F16.9,1X),A16,1X,F16.9)', 'Reference   ', LE_ref, 'Period T=', bf%T
   call OTD_map(bfp, OTD_in, Tend, OTD_out, t_FTLE, if_rst=.true.)
   call rename(report_file_OTD, 'example/roessler/PO_OTD.txt')
   call rename(report_file_OTD_LE, 'example/roessler/PO_LE.txt')
   print *, ''

   ! Reset timers
   call timer%stop('OTD modes - periodic orbit')
   call sys_floquet%reset_timer()
   call sys_floquet%jacobian%reset_timer()
   call timer%reset_all()
   call timer%add_timer('OTD modes - route to chaos', start=.true.)

   print *, ''
   print *, '########################################################################################'
   print *, '#                  Optimally Time-Dependent (OTD) modes on Route to Chaos              #'
   print *, '########################################################################################'
   print *, ''

   ! We use the old converged basis
   call orthonormalize_basis(OTD_in)

   ! We need long enough for the orbit to return to the chaotic attractor
   Tend = 60.0_dp*bf%T
   t_FTLE = bf%T
   call write_header(); call write_header_LE()
   print '(*(A16,1X))', ' ', 'FTLE_1', 'FTLE_2'
   print '(A16,1X,*(F16.9,1X))', 'Reference   ', LE_ref
   print *, ''
   call OTD_map(bfp, OTD_in, Tend, OTD_out, t_FTLE, if_rst=.false.) ! we do not reset the bf!
   call rename(report_file_OTD, 'example/roessler/PO-chaos_OTD.txt')
   call rename(report_file_OTD_LE, 'example/roessler/PO-chaos_LE.txt')
   print *, ''

   ! Print timing info for system evaulations
   print *, ''
   print '(A)', 'Timing information for sys_jac:'
   print *, ''
   call sys_jac%finalize_timer()
   call sys_jac%jacobian%finalize_timer()
   ! Print timing info for system evaulations
   print *, ''
   print '(A)', 'Timing information for sys_floquet:'
   print *, ''
   call sys_floquet%finalize_timer()
   call sys_floquet%jacobian%finalize_timer()
   ! Finalize timing
   call timer%finalize()

end program demo
