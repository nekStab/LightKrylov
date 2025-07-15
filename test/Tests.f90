program Tester
   ! Fortran best practice.
   use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
   use stdlib_logger, only: information_level, warning_level, debug_level, error_level, none_level
   ! Unit-test utility.
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   ! Abstract implementation of Krylov-based techniques.
   use LightKrylov
   use LightKrylov_Logger
   ! Implementation of the vector types and associated matrices.
   use TestVectors
   use TestLinops
   use TestKrylov
   use TestIterativeSolvers
   use TestExpmlib
   use TestNewtonKrylov
   use TestSpecialMatrices

   implicit none

   ! Unit-test related.
   integer :: status, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   ! Display information about the version of LightKrylov being tested.
   call greetings()

   ! Turn off logging during tests (unless you REALLY want it)
   call logger_setup(log_level=error_level, log_timestamp=.false.); write(*,*) 'Logging set to error_level.'
   write(*,*) ""; write(*,*) ""

   status = 0

   !----------------------------------------------------
   !-----     REAL SINGLE-PRECISION TEST SUITE     -----
   !----------------------------------------------------

   testsuites = [ &
                new_testsuite("Real Vector (sp) Test Suite", collect_vector_rsp_testsuite) , &
                new_testsuite("Real Linops (sp) Test Suite", collect_linop_rsp_testsuite), &
                new_testsuite("Real QR (sp) Test Suite", collect_qr_rsp_testsuite), &
                new_testsuite("Real Arnoldi (sp) Test Suite", collect_arnoldi_rsp_testsuite), &
                new_testsuite("Real Lanczos bidiagonalization (sp) Test Suite", collect_lanczos_bidiag_rsp_testsuite), &
                new_testsuite("Real Lanczos tridiagonalization (sp) Test Suite", collect_lanczos_tridiag_rsp_testsuite), &
                new_testsuite("Real EVP (sp) Test Suite", collect_eig_rsp_testsuite), &
                new_testsuite("Real SVD (sp) Test Suite", collect_svd_rsp_testsuite), &
                new_testsuite("Real GMRES (sp) Test Suite", collect_gmres_rsp_testsuite), &
                new_testsuite("Real CG (sp) Test Suite", collect_cg_rsp_testsuite), &
                ! new_testsuite("Real Expm (sp) Test Suite", collect_expm_rsp_testsuite), & 
                new_testsuite("Real Sqrtm (sp) Test Suite", collect_sqrtm_rsp_testsuite), &
                new_testsuite("Real Newton-Krylov fixed-point iteration (sp) Test Suite", collect_newton_rsp_testsuite) &
                ]

   write(output_unit, *) "----------------------------------------------------------------"
   write(output_unit, *) "-----                                                      -----"
   write(output_unit, *) "-----     RUNNING THE REAL SINGLE-PRECISION TEST SUITE     -----"
   write(output_unit, *) "-----                                                      -----"
   write(output_unit, *) "----------------------------------------------------------------"
   write(output_unit, *)

   do is = 1, size(testsuites)
      write (*, *) "-------------------------------"
      write (error_unit, fmt) "Testing :", testsuites(is)%name
      write (*, *) "-------------------------------"
      write (*, *)
      call run_testsuite(testsuites(is)%collect, error_unit, status)
      write (*, *)
   end do

   if (status > 0) then
      write (error_unit, '(i0, 1x, a)') status, "test(s) failed!"
      error stop
   else if (status == 0) then
      write (*, *) "All tests successfully passed!"
      write (*, *)
   end if

   !----------------------------------------------------
   !-----     REAL DOUBLE-PRECISION TEST SUITE     -----
   !----------------------------------------------------

   testsuites = [ &
                new_testsuite("Real Vector (dp) Test Suite", collect_vector_rdp_testsuite) , &
                new_testsuite("Real Linops (dp) Test Suite", collect_linop_rdp_testsuite), &
                new_testsuite("Real QR (dp) Test Suite", collect_qr_rdp_testsuite), &
                new_testsuite("Real Arnoldi (dp) Test Suite", collect_arnoldi_rdp_testsuite), &
                new_testsuite("Real Lanczos bidiagonalization (dp) Test Suite", collect_lanczos_bidiag_rdp_testsuite), &
                new_testsuite("Real Lanczos tridiagonalization (dp) Test Suite", collect_lanczos_tridiag_rdp_testsuite), &
                new_testsuite("Real EVP (dp) Test Suite", collect_eig_rdp_testsuite), &
                new_testsuite("Real SVD (dp) Test Suite", collect_svd_rdp_testsuite), &
                new_testsuite("Real GMRES (dp) Test Suite", collect_gmres_rdp_testsuite), &
                new_testsuite("Real CG (dp) Test Suite", collect_cg_rdp_testsuite), &
                new_testsuite("Real Expm (dp) Test Suite", collect_expm_rdp_testsuite), &
                new_testsuite("Real Sqrtm (dp) Test Suite", collect_sqrtm_rdp_testsuite), &
                new_testsuite("Real Newton-Krylov fixed-point iteration (dp) Test Suite", collect_newton_rdp_testsuite), &
                new_testsuite("Special Matrices (dp) Test Suite", collect_specialmatrices_rdp_testsuite) &
                ]

   write(output_unit, *) "----------------------------------------------------------------"
   write(output_unit, *) "-----                                                      -----"
   write(output_unit, *) "-----     RUNNING THE REAL DOUBLE-PRECISION TEST SUITE     -----"
   write(output_unit, *) "-----                                                      -----"
   write(output_unit, *) "----------------------------------------------------------------"
   write(output_unit, *)

   do is = 1, size(testsuites)
      write (*, *) "-------------------------------"
      write (error_unit, fmt) "Testing :", testsuites(is)%name
      write (*, *) "-------------------------------"
      write (*, *)
      call run_testsuite(testsuites(is)%collect, error_unit, status)
      write (*, *)
   end do

   if (status > 0) then
      write (error_unit, '(i0, 1x, a)') status, "test(s) failed!"
      error stop
   else if (status == 0) then
      write (*, *) "All tests successfully passed!"
      write (*, *)
   end if

   !-----------------------------------------------------
   !-----     COMPLEX SING-PRECISION TEST SUITE     -----
   !-----------------------------------------------------

   testsuites = [ &
                new_testsuite("Complex Vector (sp) Test Suite", collect_vector_csp_testsuite),  &
                new_testsuite("Complex Linops (sp) Test Suite", collect_linop_csp_testsuite),  &
                new_testsuite("Complex QR (sp) Test Suite", collect_qr_csp_testsuite), &
                new_testsuite("Complex Arnoldi (sp) Test Suite", collect_arnoldi_csp_testsuite), &
                ! new_testsuite("Complex Lanczos bidiagonalization (sp) Test Suite", collect_lanczos_bidiag_csp_testsuite), &
                new_testsuite("Complex Lanczos tridiagonalization (sp) Test Suite", collect_lanczos_tridiag_csp_testsuite), &
                new_testsuite("Complex GMRES (sp) Test Suite", collect_gmres_csp_testsuite), &
                new_testsuite("Complex CG (sp) Test Suite", collect_cg_csp_testsuite), &
                new_testsuite("Complex Expm. (sp) Test Suite", collect_expm_csp_testsuite), &
                new_testsuite("Complex Sqrtm (sp) Test Suite", collect_sqrtm_csp_testsuite) &
                ]

   write(output_unit, *) "-------------------------------------------------------------------"
   write(output_unit, *) "-----                                                         -----"
   write(output_unit, *) "-----     RUNNING THE COMPLEX SINGLE-PRECISION TEST SUITE     -----"
   write(output_unit, *) "-----                                                         -----"
   write(output_unit, *) "-------------------------------------------------------------------"
   write(output_unit, *)

   do is = 1, size(testsuites)
      write (*, *) "-------------------------------"
      write (error_unit, fmt) "Testing :", testsuites(is)%name
      write (*, *) "-------------------------------"
      write (*, *)
      call run_testsuite(testsuites(is)%collect, error_unit, status)
      write (*, *)
   end do

   if (status > 0) then
      write (error_unit, '(i0, 1x, a)') status, "test(s) failed!"
      error stop
   else if (status == 0) then
      write (*, *) "All tests successfully passed!"
      write (*, *)
   end if

   !-------------------------------------------------------
   !-----     COMPLEX DOUBLE-PRECISION TEST SUITE     -----
   !-------------------------------------------------------

   testsuites = [ &
                new_testsuite("Complex Vector (dp) Test Suite", collect_vector_cdp_testsuite),  &
                new_testsuite("Complex Linops (dp) Test Suite", collect_linop_cdp_testsuite),  &
                new_testsuite("Complex QR (dp) Test Suite", collect_qr_cdp_testsuite), &
                new_testsuite("Complex Arnoldi (dp) Test Suite", collect_arnoldi_cdp_testsuite), &
                ! new_testsuite("Complex Lanczos bidiagonalization (dp) Test Suite", collect_lanczos_bidiag_cdp_testsuite), &
                new_testsuite("Complex Lanczos tridiagonalization (dp) Test Suite", collect_lanczos_tridiag_cdp_testsuite), &
                new_testsuite("Complex GMRES (dp) Test Suite", collect_gmres_cdp_testsuite), &
                new_testsuite("Complex CG (dp) Test Suite", collect_cg_cdp_testsuite), &
                new_testsuite("Complex Expm. (dp) Test Suite", collect_expm_cdp_testsuite), &
                new_testsuite("Complex Sqrtm (dp) Test Suite", collect_sqrtm_cdp_testsuite) &
                ]

   write(output_unit, *) "-------------------------------------------------------------------"
   write(output_unit, *) "-----                                                         -----"
   write(output_unit, *) "-----     RUNNING THE COMPLEX DOUBLE-PRECISION TEST SUITE     -----"
   write(output_unit, *) "-----                                                         -----"
   write(output_unit, *) "-------------------------------------------------------------------"
   write(output_unit, *)

   do is = 1, size(testsuites)
      write (*, *) "-------------------------------"
      write (error_unit, fmt) "Testing :", testsuites(is)%name
      write (*, *) "-------------------------------"
      write (*, *)
      call run_testsuite(testsuites(is)%collect, error_unit, status)
      write (*, *)
   end do

   if (status > 0) then
      write (error_unit, '(i0, 1x, a)') status, "test(s) failed!"
      error stop
   else if (status == 0) then
      write (*, *) "All tests successfully passed!"
      write (*, *)
   end if


end program Tester
