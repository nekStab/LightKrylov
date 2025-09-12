program demo
   use stdlib_math, only: logspace
   use stdlib_logger, only: error_level
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
   real(dp) :: kp_avg, err_avg
   real(dp) :: norm_rk, norm_krylov, ip, err_norm, nacc_avg, nrej_avg

   ! Parameters
   integer, parameter :: nruns = 100     ! number of tests
   logical, parameter :: trans = .false.

   ! Set up logging
   call logger_setup(log_level=error_level)

   ! Initialize GL params
   call initialize_parameters()

   ! Define set of integration horizons and tolerances
   ntau = 3
   ntol = 5
   tauv = logspace(-2.0_wp, 0.0_wp, ntau, 10)
   tolv = logspace(-14.0_wp, -6.0_wp, ntol, 10)
   tolv = tolv(ntol:1:-1) ! reverse vector

   kdim = 150

   ! Set the Ginzburg-Landau operator for the Krylov exponential operator
   exptA_krylov = exponential_prop_krylov(GL)

   print '(3x,3(A6,1x),6x,5(A8,1x),A16)', '# runs', 'nx', 'kdim', 'dt', 'tol', 'n_acc', 'n_rej', 'kp', '||err||_2'
   do i = 1, ntau
      tau = tauv(i)
      do j = 1, ntol
         tol = tolv(j)

         nacc = 0
         nrej = 0
         kp   = 0
         err_avg = 0.0_dp
         do nrun = 1, nruns
            ! create a random input vector
            call vec%rand(ifnorm=.true.)
            
            ! integrate with RKlib
            exptA_rk%tol = tol
            call exptA(vec_rk, exptA_rk, vec, tau, info, trans)
            nacc = nacc + exptA_rk%n_accepted
            nrej = nrej + exptA_rk%n_rejected
            
            norm_rk = vec_rk%norm()
               
            ! integration with kexpm
            exptA_krylov%tol  = tol
            exptA_krylov%kdim = kdim
            call exptA(vec_krylov, exptA_krylov, vec, tau, info, trans)

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
         
         print '(3x,3(I6,1x),6x,F8.6,1x,E8.1,1x,3(F8.2,1x),E16.8))', nruns, nx, kdim, tau, tol, nacc_avg, nrej_avg, kp_avg, err_avg
      end do
      print *, ''
   end do

end program demo
