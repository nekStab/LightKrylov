module TestKrylov
  use LightKrylov
  use TestVector
  use TestMatrices
  use testdrive, only : new_unittest, unittest_type, error_type, check
  implicit none

  private

  public :: collect_power_iteration_testsuite, &
            collect_arnoldi_testsuite,         &
            collect_lanczos_testsuite,         &
            collect_krylov_schur_testsuite

contains

  !------------------------------------------------------
  !-----                                            -----
  !-----     TEST SUITE FOR THE POWER ITERATION     -----
  !-----                                            -----
  !------------------------------------------------------

  subroutine collect_power_iteration_testsuite(testsuite)
    !> Collection of tests.
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    return
  end subroutine collect_power_iteration_testsuite

  !------------------------------------------------------------
  !-----                                                  -----
  !-----     TEST SUITE FOR THE ARNOLDI FACTORIZATION     -----
  !-----                                                  -----
  !------------------------------------------------------------

  subroutine collect_arnoldi_testsuite(testsuite)
    !> Collection of tests.
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    return
  end subroutine collect_arnoldi_testsuite

  !------------------------------------------------------------
  !-----                                                  -----
  !-----     TEST SUITE FOR THE LANCZOS FACTORIZATION     -----
  !-----                                                  -----
  !------------------------------------------------------------

  subroutine collect_lanczos_testsuite(testsuite)
    !> Collection of tests.
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    return
  end subroutine collect_lanczos_testsuite

  !-----------------------------------------------------------------
  !-----                                                       -----
  !-----     TEST SUITE FOR THE KRYLOV-SCHUR FACTORIZATION     -----
  !-----                                                       -----
  !-----------------------------------------------------------------

  subroutine collect_krylov_schur_testsuite(testsuite)
    !> Collection of tests.
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    return
  end subroutine collect_krylov_schur_testsuite

end module TestKrylov
