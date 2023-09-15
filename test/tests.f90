program Tester
  use, intrinsic :: iso_fortran_env, only : error_unit
  use testdrive, only : run_testsuite, new_testsuite, testsuite_type
  use LightKrylov
  use TestVector, only : collect_vector_testsuite
  implicit none

  integer :: status, is
  type(testsuite_type), allocatable :: testsuites(:)
  character(len=*), parameter :: fmt = '("+", *(1x, a))'

  call greetings()

  status = 0

  testsuites = [new_testsuite("Real Vector Test Suite", collect_vector_testsuite)]

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
