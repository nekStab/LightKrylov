program Tester
  !> Fortran best practice.
  use, intrinsic :: iso_fortran_env, only : error_unit
  !> Unit-test utility.
  use testdrive                    , only : run_testsuite, new_testsuite, testsuite_type
  !> Abstract implementation of Krylov-based techniques.
  use LightKrylov
  !> Implementation of simple 3-dimensional vector types and associated matrices.
  use TestVector                   , only : collect_real_vector_testsuite
  use TestMatrices                 , only : collect_real_matrix_testsuite
  use TestKrylov                   , only : collect_power_iteration_testsuite, collect_arnoldi_testsuite, &
       collect_lanczos_testsuite
  use TestIterativeSolvers         , only : collect_evp_testsuite

  implicit none

  !> Unit-test related.
  integer :: status, is
  type(testsuite_type), allocatable :: testsuites(:)
  character(len=*), parameter :: fmt = '("+", *(1x, a))'

  ! --> Display information about the version of LightKrylov being tested.
  call greetings()


  status = 0

  testsuites = [&
       new_testsuite("Real Vector Test Suite", collect_real_vector_testsuite),         &
       new_testsuite("Real Matrix Test Suite", collect_real_matrix_testsuite),         &
       new_testsuite("Power iteration Test Suite", collect_power_iteration_testsuite), &
       new_testsuite("Arnoldi Test Suite", collect_arnoldi_testsuite)                , &
       new_testsuite("Lanczos Test Suite", collect_lanczos_testsuite)                , &
       new_testsuite("Eigenvalues Test Suite", collect_evp_testsuite)                  &
       ]

  do is = 1, size(testsuites)
     write(*, *)            "-------------------------------"
     write(error_unit, fmt) "Testing :", testsuites(is)%name
     write(*, *)            "-------------------------------"
     write(*, *)
     call run_testsuite(testsuites(is)%collect, error_unit, status)
     write(*, *)
  enddo

  if (status > 0) then
     write(error_unit, '(i0, 1x, a)') status, "test(s) failed!"
     error stop
  else if (status == 0) then
     write(*, *) "All tests successfully passed!"
     write(*, *)
  endif

end program Tester
