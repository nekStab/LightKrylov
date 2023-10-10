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
       collect_bicgstab_testsuite, &
       collect_evp_testsuite

contains

  !-------------------------------------------------------------
  !-----                                                   -----
  !-----     TEST SUITE FOR GENERAL EIGENVALUE PROBLEM     -----
  !-----                                                   -----
  !-------------------------------------------------------------

  subroutine collect_evp_testsuite(testsuite)
    !> Collection of tests.
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
         new_unittest("Sym. EVP computation", test_spd_evp_problem), &
         new_unittest("EVP computation", test_evp_problem)           &
         ]
    return
  end subroutine collect_evp_testsuite

  subroutine test_spd_evp_problem(error)
    !> Error type to be returned.
    type(error_type), allocatable, intent(out) :: error
    !> Test matrix.
    class(spd_matrix), allocatable :: A
    !> Krylov subspace.
    class(rvector), allocatable :: X(:)
    !> Eigenvalues and eigenvectors.
    real(kind=wp) :: evals(test_size), evecs(test_size, test_size)
    !> Residuals.
    real(kind=wp) :: residuals(test_size)
    !> Information flag.
    integer :: info
    !> Toeplitz matrix.
    real(kind=wp) :: T(test_size, test_size), a_, b_
    !> Miscellaneous.
    integer :: i, j, k
    real(kind=wp) :: alpha
    real(kind=wp) :: true_evals(test_size), pi

    !> Create the sym. pos. def. Toeplitz matrix.
    call random_number(a_) ; call random_number(b_) ; b_ = -abs(b_)
    do i = 1, test_size
       !> Diagonal entry.
       T(i, i) = a_
       !> Upper diagonal entry.
       if (i < test_size) T(i, i+1) = b_
       !> Lower diagonal entry.
       if (i < test_size) T(i+1, i) = b_
    enddo

    !> Test matrix.
    A = spd_matrix(T)

    !> Initialize Krylov subspace.
    allocate(X(1:test_size+1))
    call random_number(X(1)%data) ; alpha = X(1)%norm() ; call X(1)%scal(1.0_wp / alpha)
    do i = 2, size(X)
       call X(i)%zero()
    enddo

    !> Initialize internal variables.
    evals = 0.0_wp ; evecs = 0.0_wp ; residuals = 0.0_wp

    !> Compute spectral decomposition.
    call eighs(A, X, evecs, evals, residuals, info)

    !> Analytical eigenvalues.
    true_evals = 0.0_wp ; pi = 4.0_wp * atan(1.0_wp)
    do i = 1, test_size
       true_evals(i) = a_ + 2*abs(b_) * cos(i*pi / (test_size+1))
    enddo

    !> Check correctness.
    call check(error, all_close(evals, true_evals, rtol, atol))

    return
  end subroutine test_spd_evp_problem

  subroutine test_evp_problem(error)
    !> Error type to be returned.
    type(error_type), allocatable, intent(out) :: error
    !> Test matrix.
    class(rmatrix), allocatable :: A
    !> Krylov subspace.
    class(rvector), allocatable :: X(:)
    !> Eigenvalues and eigenvectors.
    complex(kind=wp) :: evals(test_size), evecs(test_size, test_size)
    !> Residuals.
    real(kind=wp) :: residuals(test_size)
    !> Information flag.
    integer :: info
    !> Toeplitz matrix.
    real(kind=wp) :: T(test_size, test_size), a_, b_
    !> Miscellaneous.
    integer :: i, j, k
    real(kind=wp) :: alpha
    complex(kind=wp) :: true_evals(test_size), pi

    !> Create the sym. pos. def. Toeplitz matrix.
    call random_number(a_) ; call random_number(b_) ; b_ = abs(b_)
    do i = 1, test_size
       !> Diagonal entry.
       T(i, i) = a_
       !> Upper diagonal entry.
       if (i < test_size) T(i, i+1) = b_
       !> Lower diagonal entry.
       if (i < test_size) T(i+1, i) = -b_
    enddo

    !> Test matrix.
    A = rmatrix(T)

    !> Initialize Krylov subspace.
    allocate(X(1:test_size+1))
    call random_number(X(1)%data) ; alpha = X(1)%norm() ; call X(1)%scal(1.0_wp / alpha)
    do i = 2, size(X)
       call X(i)%zero()
    enddo

    !> Initialize internal variables.
    evals = cmplx(0.0_wp, 0.0_wp, kind=wp) ; evecs = cmplx(0.0_wp, 0.0_wp, kind=wp) ; residuals = 0.0_wp

    !> Compute spectral decomposition.
    call eigs(A, X, evecs, evals, residuals, info)

    !> Analytical eigenvalues.
    true_evals = cmplx(0.0_wp, 0.0_wp, kind=wp) ; pi = 4.0_wp * atan(1.0_wp)
    k = 1
    do i = 1, test_size, 2
       true_evals(i)   = a_*cmplx(1.0_wp, 0.0_wp, kind=wp) + (2_wp*b_*cos(k*pi / (test_size+1)))*cmplx(0.0_wp, 1.0_wp, kind=wp)
       true_evals(i+1)   = a_*cmplx(1.0_wp, 0.0_wp, kind=wp) - (2_wp*b_*cos(k*pi / (test_size+1)))*cmplx(0.0_wp, 1.0_wp, kind=wp)
       k = k +1
    enddo
    
    !> Check correctness.
    call check(error, norm2(abs(evals - true_evals)) < rtol)

    return
  end subroutine test_evp_problem

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
    type(bicgstab_opts) :: opts
    !> Information flag.
    integer :: info
    
    ! --> Initialize linear problem.
    A = rmatrix(); call random_number(A%data)
    b = rvector(); call random_number(b%data)
    x = rvector(); call x%zero()
    ! --> BiCGSTAB solver.
    opts = bicgstab_opts(verbose=.true., atol=1e-12_wp, rtol=0.0_wp)
    call bicgstab(A, b, x, info, options=opts)
    ! --> Check convergence.
    call check(error, norm2(matmul(A%data, x%data) - b%data)**2 < rtol)
    
    return
  end subroutine test_bicgstab_full_computation

end module TestIterativeSolvers
