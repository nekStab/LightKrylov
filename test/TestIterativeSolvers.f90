module TestIterativeSolvers
  use LightKrylov
  use TestVector
  use TestMatrices
  use testdrive  , only : new_unittest, unittest_type, error_type, check
  use stdlib_math, only : all_close
  implicit none

  private

  public :: collect_gmres_testsuite,    &
       collect_cg_testsuite,       &
       collect_bicgstab_testsuite

contains

  !-------------------------------------------------------------
  !-----                                                   -----
  !-----     TEST SUITE FOR GENERAL EIGENVALUE PROBLEM     -----
  !-----                                                   -----
  !-------------------------------------------------------------

  !----------------------------------------------------------
  !-----                                                -----
  !-----     TEST FOR SYMMETRIC EIGENVALUE PROBLEMS     -----
  !-----                                                -----
  !----------------------------------------------------------


  !------------------------------------
  !-----                          -----
  !-----     GMRES TEST SUITE     -----
  !-----                          -----
  !------------------------------------

  subroutine collect_gmres_testsuite(testsuite)
    !> Collection of tests.
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [&
         new_unittest("GMRES full computation", test_gmres_full_computation),                                    &
         new_unittest("GMRES full computation w. sym. pos. def. matrix", test_gmres_full_computation_spd_matrix) &
         ]
    return
  end subroutine collect_gmres_testsuite

  subroutine test_gmres_full_computation(error)
    !> Error type to be returned.
    type(error_type), allocatable, intent(out) :: error
    !> Linear Problem.
    class(rmatrix), allocatable :: A ! Linear Operator.
    class(rvector), allocatable :: b ! Right-hand side vector.
    class(rvector), allocatable :: x ! Solution vector.
    !> GMRES options.
    type(gmres_opts) :: opts
    !> Information flag.
    integer :: info

    ! --> Initialize linear problem.
    A = rmatrix() ; call random_number(A%data)
    b = rvector() ; call random_number(b%data)
    x = rvector() ; call x%zero()
    ! --> GMRES solver.
    opts = gmres_opts(kdim=test_size, verbose=.false.)
    call gmres(A, b, x, info, opts)
    ! --> Check convergence.
    call check(error, norm2(matmul(A%data, x%data) - b%data)**2 < rtol)

    return
  end subroutine test_gmres_full_computation

  subroutine test_gmres_full_computation_spd_matrix(error)
    !> Error type to be returned.
    type(error_type), allocatable, intent(out) :: error
    !> Linear Problem.
    class(spd_matrix), allocatable :: A ! Linear Operator.
    class(rvector), allocatable :: b ! Right-hand side vector.
    class(rvector), allocatable :: x ! Solution vector.
    !> GMRES options.
    type(gmres_opts) :: opts
    !> Information flag.
    integer :: info

    ! --> Initialize linear problem.
    A = spd_matrix() ; call random_number(A%data) ; A%data = matmul(A%data, transpose(A%data))
    b = rvector()    ; call random_number(b%data)
    x = rvector()    ; call x%zero()
    ! --> GMRES solver.
    opts = gmres_opts(kdim=test_size, verbose=.false.)
    call gmres(A, b, x, info, opts)
    ! --> Check convergence.
    call check(error, norm2(matmul(A%data, x%data) - b%data)**2 < rtol)

    return
  end subroutine test_gmres_full_computation_spd_matrix

  !----------------------------------
  !-----                        -----
  !-----     SVD TEST SUITE     -----
  !-----                        -----
  !----------------------------------

  
  !------------------------------------
  !-----                          -----
  !-----     CG TEST SUITE        -----
  !-----                          -----
  !------------------------------------
  
  subroutine collect_cg_testsuite(testsuite)
    !> Collection of tests.
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    
    testsuite = [ &
         new_unittest("CG full computation w. s.p.d. matrix", test_cg_full_computation_spd_matrix) &
         ]
    return
  end subroutine collect_cg_testsuite
  
  subroutine test_cg_full_computation_spd_matrix(error)
    !> Error type to be returned.
    type(error_type), allocatable, intent(out) :: error
    !> Linear Problem.
    real(kind=wp) :: Q(test_size, 10_wp*test_size)
    class(spd_matrix), allocatable :: A ! Linear Operator.
    class(rvector), allocatable :: b ! Right-hand side vector.
    class(rvector), allocatable :: x ! Solution vector.
    type(cg_opts)               :: opts
    !> Information flag.
    integer :: info

    ! --> Initialize linear problem.
    A = spd_matrix()
    call random_number(Q) ; A%data = matmul(Q, transpose(Q))
    A%data = 0.5 * (A%data + transpose(A%data))
    b = rvector(); call random_number(b%data)
    x = rvector(); call x%zero()
    ! --> CG solver.
    opts = cg_opts(verbose=.false., atol=1e-12_wp, rtol=0.0_wp)
    call cg(A, b, x, info, options=opts)
    ! --> Check convergence.
    call check(error, norm2(matmul(A%data, x%data) - b%data)**2 < rtol)
    
    return
  end subroutine test_cg_full_computation_spd_matrix
  
  !--------------------------------------
  !-----                            -----
  !-----     BICGSTAB TEST SUITE    -----
  !-----                            -----
  !--------------------------------------
  
  subroutine collect_bicgstab_testsuite(testsuite)
    !> Collection of tests.
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    
    testsuite = [ &
         new_unittest("BiCGSTAB full computation", test_bicgstab_full_computation) &
         ]
    return
  end subroutine collect_bicgstab_testsuite
  
  subroutine test_bicgstab_full_computation(error)
    !> Error type to be returned.
    type(error_type), allocatable, intent(out) :: error
    !> Linear Problem.
    class(rmatrix), allocatable :: A ! Linear Operator.
    class(rvector), allocatable :: b ! Right-hand side vector.
    class(rvector), allocatable :: x ! Solution vector.
    !> Information flag.
    integer :: info
    
    ! --> Initialize linear problem.
    A = rmatrix(); call random_number(A%data)
    b = rvector(); call random_number(b%data)
    x = rvector(); call x%zero()
    ! --> BiCGSTAB solver.
    call bicgstab(A, b, x, info)
    ! --> Check convergence.
    call check(error, norm2(matmul(A%data, x%data) - b%data)**2 < rtol)
    
    return
  end subroutine test_bicgstab_full_computation

end module TestIterativeSolvers
