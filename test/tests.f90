program Tester
   !> Fortran best practice.
   use, intrinsic :: iso_fortran_env, only: error_unit
   !> Unit-test utility.
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   !> Abstract implementation of Krylov-based techniques.
   use LightKrylov
   !> Implementation of simple 3-dimensional vector types and associated matrices.
   use TestVector
   use TestMatrices
   use TestKrylov
   use TestExpm
   use TestIterativeSolvers

   implicit none

   !> Unit-test related.
   integer :: status, is, num_tests
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("+", *(1x, a))'

   ! --> Display information about the version of LightKrylov being tested.
   call greetings()

   status = 0

   testsuites = [ &
                new_testsuite("Real Vector Test Suite", collect_real_vector_testsuite), &
                new_testsuite("Real Matrix Test Suite", collect_real_matrix_testsuite), &
                new_testsuite("Operations on Abstract Lin. Op.", collect_abstract_linop_operations_testsuite), &
                new_testsuite("Arnoldi Test Suite", collect_arnoldi_testsuite), &
                new_testsuite("Two-sided Arnoldi Test Suite", collect_two_sided_arnoldi_testsuite), &
                !new_testsuite("Rational Arnoldi Test Suite", collect_rational_arnoldi_testsuite),         &
                new_testsuite("Lanczos tridiagonalization Test Suite", collect_lanczos_tridiag_testsuite), &
                new_testsuite("Lanczos bidiagonalization Test Suite", collect_lanczos_bidiag_testsuite), &
                new_testsuite("Eigenvalues Test Suite", collect_evp_testsuite), &
                new_testsuite("GMRES Test Suite", collect_gmres_testsuite), &
                new_testsuite("SVD Test Suite", collect_svd_testsuite), &
                new_testsuite("CG Test Suite", collect_cg_testsuite), &
                !new_testsuite("BICGSTAB Test Suite", collect_bicgstab_testsuite)                          &
                !new_testsuite("Non-symetric Lanczos Test Suite", collect_nonsymmetric_lanczos_testsuite),  &
                new_testsuite("QR factorization Test Suite", collect_qr_testsuite), &
                new_testsuite("Dense Matrix Functions Test Suite", collect_dense_matrix_functions_testsuite), &
                new_testsuite("Krylov Matrix Exponential Test Suite", collect_kexpm_testsuite) &
                ]

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
