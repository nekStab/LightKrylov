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

    testsuite = [&
         new_unittest("Power Iteration w/ real 3x3 matrix", test_real_matrix_power_iteration), &
         new_unittest("Power Iteration w/ 3x3 Strang matrix", test_spd_matrix_power_iteration) &
         ]
    return
  end subroutine collect_power_iteration_testsuite

  subroutine test_real_matrix_power_iteration(error)
    !> Error type to be returned.
    type(error_type), allocatable, intent(out) :: error
    !> Test matrix.
    class(rmatrix), allocatable :: A
    !> Starting vector.
    class(rvector), allocatable :: x
    !> Estimated eigenvalue.
    double precision :: lambda
    !> Maximum number of iterations.
    integer :: niter = 50
    !> Information flag.
    integer :: info

    ! --> Initialize matrix.
    A = rmatrix(reshape([1, 2, 0, &
                        -2, 1, 2, &
                         1, 3, 1], shape=[3, 3], order=[2, 1])) ! order=[2, 1] -> fill the matrix row-wise.
    ! --> Initialize vector.
    x = rvector([1, 1, 1])
    ! --> Power iteration method.
    call power_iteration(A, x, lambda, niter, info)
    ! --> Check result.
    call check(error, info == 0)

    return
  end subroutine test_real_matrix_power_iteration

  subroutine test_spd_matrix_power_iteration(error)
    !> Error type to be returned.
    type(error_type), allocatable, intent(out) :: error
    !> Test matrix.
    class(spd_matrix), allocatable :: A
    !> Starting vector.
    class(rvector), allocatable :: x
    !> Estimated eigenvalue.
    double precision :: lambda
    !> Maximum number of iterations.
    integer :: niter = 50
    !> Information flag.
    integer :: info

    ! --> Initialize matrix.
    A = spd_matrix(reshape([2, -1,  0, &
                           -1,  2, -1, &
                            0, -1,  2 ], shape=[3, 3]))
    ! --> Initialize vector.
    x = rvector([1, 1, 1])
    ! --> Power iteration method.
    call power_iteration(A, x, lambda, niter, info)
    ! --> Check results.
    call check(error, info == 0)

    return
  end subroutine test_spd_matrix_power_iteration

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
