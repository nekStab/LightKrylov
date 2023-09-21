module TestIterativeSolvers
  use LightKrylov
  use TestVector
  use TestMatrices
  use testdrive  , only : new_unittest, unittest_type, error_type, check
  use stdlib_math, only : all_close
  implicit none

  private

  public :: collect_evp_testsuite, collect_gmres_testsuite

contains

  !-------------------------------------------------------------
  !-----                                                   -----
  !-----     TEST SUITE FOR GENERAL EIGENVALUE PROBLEM     -----
  !-----                                                   -----
  !-------------------------------------------------------------

  subroutine collect_evp_testsuite(testsuite)
    !> Collection of tests.
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [&
         new_unittest("Triangular matrix eigenvalues", test_eigvals_triangular_matrix), &
         new_unittest("Sym. Pos. Def. matrix eigenvalues", test_eigvals_spd_matrix)     &
         ]
    return
  end subroutine collect_evp_testsuite

  subroutine test_eigvals_triangular_matrix(error)
    !> Error type to be returned.
    type(error_type), allocatable, intent(out) :: error
    !> Test matrix.
    class(rmatrix), allocatable :: A
    !> Krylov subspace.
    class(rvector), dimension(:), allocatable :: X
    !> Krylov subspace dimension.
    integer, parameter :: kdim = 3
    !> Eigenvectors, eigenvalues, and residuals.
    double complex  , dimension(kdim, kdim) :: eigvecs
    double complex  , dimension(kdim)       :: eigvals
    double precision, dimension(kdim)       :: residuals
    !> Information flag.
    integer :: info

    !> Miscellaneous.
    double precision :: alpha
    double precision, dimension(3), parameter :: true_eigs = [6, -5, 3]

    ! --> Initialize matrix.
    A = rmatrix() ; A%data = reshape([-2, -2, 4, -4, 1, 2, 2, 2, 5], shape=[3, 3])
    ! --> Initialize Krylov subspace.
    allocate(X(1:kdim+1)) ; call random_number(X(1)%data)
    alpha = X(1)%norm() ; call X(1)%scalar_mult(1.0D+00 / alpha)
    ! --> Compute eigenpairs.
    call eigs(A, X, eigvecs, eigvals, residuals, info)
    ! --> Check results.
    call check(error, all_close(real(eigvals), true_eigs), .true.)
    return
  end subroutine test_eigvals_triangular_matrix

  !----------------------------------------------------------
  !-----                                                -----
  !-----     TEST FOR SYMMETRIC EIGENVALUE PROBLEMS     -----
  !-----                                                -----
  !----------------------------------------------------------

  subroutine test_eigvals_spd_matrix(error)
    !> Error type to be returned
    type(error_type), allocatable, intent(out) :: error
    !> Test matrix.
    class(spd_matrix), allocatable :: A
    !> Krylov subspace.
    class(rvector), dimension(:), allocatable :: X
    !> Krylov subspace dimension.
    integer, parameter :: kdim = 3
    !> Eigenvectors, eigenvalues, and residuals.
    double precision, dimension(kdim, kdim) :: eigvecs
    double precision, dimension(kdim)       :: eigvals
    double precision, dimension(kdim)       :: residuals
    !> Information flag.
    integer :: info

    !> Miscellaneous.
    integer :: i, j, k
    double precision :: alpha
    double precision, dimension(3), parameter :: true_eigs = [1.5D+00, 1.0D+00, 0.5D+00]

    ! --> Initialize matrix.
    A = spd_matrix() ;
    A%data(1, 1) = 1.0D+00    ; A%data(1, 2) = -0.5D+00 ; A%data(1, 3) = 0.0D+00
    A%data(2, 1) = -0.5D+00   ; A%data(2, 2) = 1.0D+00  ; A%data(2, 3) = 0.0D+00
    A%data(3, 1) = 0.0D+00    ; A%data(3, 2) = 0.0D+00  ; A%data(3, 3) = 1.0D+00
    ! --> Initialize Krylov subspace.
    allocate(X(1:kdim+1)) ; call random_number(X(1)%data)
    X(1)%data = [1.0D+00, 2.0D+00, 1.0D+00]
    alpha = X(1)%norm() ; call X(1)%scalar_mult(1.0D+00 / alpha)
    ! --> Compute eigenpairs.
    call eighs(A, X, eigvecs, eigvals, residuals, info)
    ! --> Check results.
    call check(error, all_close(eigvals, true_eigs), .true.)
    return
  end subroutine test_eigvals_spd_matrix

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
    !> Information flag.
    integer :: info

    ! --> Initialize linear problem.
    A = rmatrix() ; call random_number(A%data)
    b = rvector() ; call random_number(b%data)
    x = rvector() ; call x%zero()
    ! --> GMRES solver.
    call gmres(A, b, x, info, kdim=3)
    ! --> Check convergence.
    call check(error, norm2(matmul(A%data, x%data) - b%data) < 1e-12)

    return
  end subroutine test_gmres_full_computation

  subroutine test_gmres_full_computation_spd_matrix(error)
    !> Error type to be returned.
    type(error_type), allocatable, intent(out) :: error
    !> Linear Problem.
    class(spd_matrix), allocatable :: A ! Linear Operator.
    class(rvector), allocatable :: b ! Right-hand side vector.
    class(rvector), allocatable :: x ! Solution vector.
    !> Information flag.
    integer :: info

    ! --> Initialize linear problem.
    A = spd_matrix() ; call random_number(A%data) ; A%data = matmul(A%data, transpose(A%data))
    b = rvector()    ; call random_number(b%data)
    x = rvector()    ; call x%zero()
    ! --> GMRES solver.
    call gmres(A, b, x, info, kdim=3)
    ! --> Check convergence.
    call check(error, norm2(matmul(A%data, x%data) - b%data) < 1e-12)

    return
  end subroutine test_gmres_full_computation_spd_matrix

  ! subroutine test_restarted_gmres_computation(error)
  !   !> Eror type to be returned.
  !   type(error_type), allocatable, intent(out) :: error
  !   !> Linear Problem.
  !   class(rmatrix), allocatable :: A ! Linear Operator.
  !   class(rvector), allocatable :: b ! Right-hand side vector.
  !   class(rvector), allocatable :: x ! Solution vector.
  !   !> Information flag.
  !   integer :: info

  !   ! --> Initialize linear problem.
  !   A = rmatrix() ; call random_number(A%data)
  !   b = rvector() ; call random_number(b%data)
  !   x = rvector() ; call x%zero()
  !   ! --> GMRES solver.
  !   call gmres(A, b, x, info, kdim=2, verbosity=.true., maxiter=100)
  !   ! --> Check convergence.
  !   call check(error, norm2(matmul(A%data, x%data) - b%data)**2 < 1e-12)
  !   return
  ! end subroutine test_restarted_gmres_computation

end module TestIterativeSolvers