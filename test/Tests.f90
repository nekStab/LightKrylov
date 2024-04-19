program Tester
   !> Fortran best practice.
   use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
   !> Unit-test utility.
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   !> Abstract implementation of Krylov-based techniques.
   use LightKrylov
   !> Implementation of the vector types and associated matrices.
   use TestVectors
   use TestLinops
   use TestKrylov

   implicit none

   !> Unit-test related.
   integer :: status, is, num_tests
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("+", *(1x, a))'

   ! --> Display information about the version of LightKrylov being tested.
   call greetings()

   status = 0

   !----------------------------------------------------
   !-----     REAL SINGLE-PRECISION TEST SUITE     -----
   !----------------------------------------------------

   testsuites = [ &
                new_testsuite("Real Vector (sp) Test Suite", collect_vector_rsp_testsuite), &
                new_testsuite("Real Linops (sp) Test Suite", collect_linop_rsp_testsuite), &
                new_testsuite("Real QR (sp) Test Suite", collect_qr_rsp_testsuite), &
                new_testsuite("Real Arnoldi (sp) Test Suite", collect_arnoldi_rsp_testsuite) &
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
                new_testsuite("Real Vector (dp) Test Suite", collect_vector_rdp_testsuite), &
                new_testsuite("Real Linops (dp) Test Suite", collect_linop_rdp_testsuite), &
                new_testsuite("Real QR (dp) Test Suite", collect_qr_rdp_testsuite), &
                new_testsuite("Real Arnoldi (dp) Test Suite", collect_arnoldi_rdp_testsuite) &
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
                new_testsuite("Complex Vector (sp) Test Suite", collect_vector_csp_testsuite), &
                new_testsuite("Complex Linops (sp) Test Suite", collect_linop_csp_testsuite), &
                new_testsuite("Complex QR (sp) Test Suite", collect_qr_csp_testsuite), &
                new_testsuite("Complex Arnoldi (sp) Test Suite", collect_arnoldi_csp_testsuite) &
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

   !----------------------------------------------------
   !-----     COMPLEX DOUBLE-PRECISION TEST SUITE     -----
   !----------------------------------------------------

   testsuites = [ &
                new_testsuite("Complex Vector (dp) Test Suite", collect_vector_cdp_testsuite),  &
                new_testsuite("Complex Linops (dp) Test Suite", collect_linop_cdp_testsuite),  &
                new_testsuite("Complex QR (dp) Test Suite", collect_qr_cdp_testsuite), &
                new_testsuite("Complex QR (dp) Test Suite", collect_arnoldi_cdp_testsuite) &
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
