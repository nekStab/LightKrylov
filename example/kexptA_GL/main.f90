program demo
   use stdlib_math, only: logspace
   use stdlib_logger, only: error_level
   use stdlib_stats, only: mean, var
   use LightKrylov
   use LightKrylov_Constants
   use LightKrylov_Logger
   use LightKrylov_Timing, only: timer => global_lightkrylov_timer
   ! GInzburg-Landau
   use Ginzburg_Landau_Base
   use Ginzburg_Landau_Operators
   implicit none

   character(len=128), parameter :: this_module = 'Example Ginzburg_Landau'

   ! Linear Ginzburg-Landau operator
   type(GL_operator) :: GL
   ! Exponential propagator using RKlib
   type(exponential_prop_rk) :: exptA_rk
   ! Exponential propagator using kexpm
   type(exponential_prop_krylov) :: exptA_krylov

   ! State vectors
   type(state_vector) :: vec, vec_rk, vec_krylov

   ! Ingration horizon
   real(dp), dimension(:), allocatable :: tauv, tolv
   real(dp) :: tau, tol

   ! Misc
   integer :: i, j, k, nrun, ntau, ntol, nk, kdim, info
   integer :: nacc, nrej, kp
   real(dp) :: kp_avg, err_avg, time
   real(dp) :: t_avg_rk, t_avg_krylov
   real(dp) :: stddev_rk, stddev_krylov
   real(dp) :: norm_rk, norm_krylov, ip, err_norm, nacc_avg, nrej_avg

   ! Parameters
   integer, parameter :: nruns = 100     ! number of tests
   logical, parameter :: trans = .false.

   real(dp) :: ttime(nruns, 2)

   ! timer
   integer :: clock_rate, clock_start, clock_stop

   call system_clock(count_rate=clock_rate)

   ! Set up logging
   call logger_setup(log_level=error_level)

   ! Initialize GL params
   call initialize_parameters()

   ! Define set of integration horizons and tolerances
   ntau = 9
   ntol = 5
   tauv = logspace(-5.0_wp, -1.0_wp, ntau, 10)
   tolv = logspace(-14.0_wp, -6.0_wp, ntol, 10)
   tolv = tolv(ntol:1:-1) ! reverse vector

   kdim = 250

   ! Set the Ginzburg-Landau operator for the Krylov exponential operator
   exptA_krylov = exponential_prop_krylov(GL)

   ! compute exp(tA) v examples
   tau = 0.01_dp
   time = 0
   call vec%rand(ifnorm=.true.)
   print *, 'exptA_example', 0, 'time', time, vec%norm()
   do nrun = 1, 5000
      time = time + tau
      call exptA(vec_krylov, exptA_krylov, vec, tau, info, trans)
      vec = vec_krylov
      print *, 'exptA_example', nrun, 'time', time, vec%norm()
   end do

   print '(3(A6,1x),6x,5(A8,1x),A16,2(A25),A15)', '# runs', 'nx', 'kdim', 'dt', 'tol', & 
         & 'n_acc', 'n_rej', 'kp', '||err||_2', &
         & 'timings (RK)', 'timings (krylov)', 'Elapsed time'
   do i = 1, ntau
      tau = tauv(i)
      do j = 1, ntol
         tol = tolv(j)

         nacc = 0
         nrej = 0
         kp   = 0
         err_avg = 0.0_dp
         ttime = 0.0_dp
         do nrun = 1, nruns
            ! create a random input vector
            call vec%rand(ifnorm=.true.)
            
            ! integrate with RKlib
            exptA_rk%tol = tol
            call system_clock(count=clock_start)     ! Start Timer
            call exptA(vec_rk, exptA_rk, vec, tau, info, trans)
            call system_clock(count=clock_stop)      ! Stop Timer
            ttime(nrun,1) = real(clock_stop-clock_start)/real(clock_rate)
            nacc = nacc + exptA_rk%n_accepted
            nrej = nrej + exptA_rk%n_rejected
            
            norm_rk = vec_rk%norm()
               
            ! integration with kexpm
            exptA_krylov%tol  = tol
            exptA_krylov%kdim = kdim
            call system_clock(count=clock_start)     ! Start Timer
            call exptA(vec_krylov, exptA_krylov, vec, tau, info, trans)
            call system_clock(count=clock_stop)      ! Stop Timer
            ttime(nrun,2) = real(clock_stop-clock_start)/real(clock_rate)

            ! compare outputs
            norm_krylov = vec_krylov%norm()
            ip          = vec_rk%dot(vec_krylov)/(norm_rk*norm_krylov)
            call vec_krylov%sub(vec_rk)

            ! cumulate data
            err_avg = err_avg + vec_krylov%norm()
            kp      = kp + exptA_krylov%info
         end do
         nacc_avg = 1.0_dp*nacc/nruns
         nrej_avg = 1.0_dp*nrej/nruns
         kp_avg   = 1.0_dp*kp/nruns
         err_avg  = err_avg/nruns
         t_avg_rk      = mean(ttime(:,1))
         t_avg_krylov  = mean(ttime(:,2))
         stddev_rk     = var(ttime(:,1))
         stddev_krylov = var(ttime(:,2))
         
         print '(3(I6,1x),6x,F8.6,1x,E8.1,1x,3(F8.2,1x),E16.8,2(4x,F8.6,1x,F12.10),3x,F12.4)', &
            & nruns, nx, kdim, tau, tol, &
            & nacc_avg, nrej_avg, kp_avg, err_avg, &
            & t_avg_rk, stddev_rk, t_avg_krylov, stddev_krylov, sum(ttime)
      end do
      !print *, ''
   end do

end program demo
