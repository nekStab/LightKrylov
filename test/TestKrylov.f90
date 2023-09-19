module TestKrylov
  use LightKrylov
  use TestVector
  use TestMatrices
  use testdrive  , only : new_unittest, unittest_type, error_type, check
  use stdlib_math, only : all_close
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

    testsuite = [&
         new_unittest("Arnoldi full factorization", test_arnoldi_full_factorization),            &
         new_unittest("Arnoldi basis orthogonality", test_arnoldi_basis_orthogonality),          &
         new_unittest("Arnoldi invariant subspace computation", test_arnoldi_invariant_subspace) &
         ]

    return
  end subroutine collect_arnoldi_testsuite

  subroutine test_arnoldi_full_factorization(error)
    ! This function checks the correctness of the Arnoldi implementation by
    ! verifying that the full factorization is correct, i.e. A @ X[1:k] = X[1:k+1] @ H.
    ! A random 3x3 matrix is used for testing.

    !> Error type to be returned.
    type(error_type), allocatable, intent(out) :: error
    !> Test matrix.
    class(rmatrix), allocatable :: A
    double precision, dimension(3, 3) :: Adata
    !> Krylov subspace.
    class(rvector), dimension(:), allocatable :: X
    !> Krylov subspace dimension.
    integer, parameter :: kdim=3
    !> Hessenberg matrix.
    double precision, dimension(kdim+1, kdim) :: H
    !> Information flag.
    integer :: info
    !> Misc.
    integer :: k
    double precision, dimension(3, kdim+1) :: Xdata
    double precision :: alpha

    ! --> Initialize matrix.
    call random_number(Adata) ; A = rmatrix(Adata)
    ! --> Initialize Krylov subspace.
    allocate(X(1:kdim+1)) ; call random_number(X(1)%data)
    alpha = X(1)%norm() ; call X(1)%scalar_mult(1.0D+00 / alpha)
    do k = 2, size(X)
       call X(k)%zero()
    enddo
    H = 0.0D+00
    ! --> Arnoldi factorization.
    call arnoldi_factorization(A, X, H, info)
    ! --> Check correctness of full factorization.
    do k = 1, kdim+1
       Xdata(:, k) = X(k)%data
    enddo
    call check(error, norm2(matmul(Adata, Xdata(:, 1:kdim)) - matmul(Xdata, H)) < 1e-12 )

    return
  end subroutine test_arnoldi_full_factorization

  subroutine test_arnoldi_invariant_subspace(error)
    ! This function checks that the Arnoldi iteration stops whenever an invariant
    ! has been computed. An arbitrary 3x3 matrix with a 2-dimensional invariant subspace
    ! is being used for testing. The starting vector for Arnoldi has component only
    ! within this 2-dimensional subspace.

    !> Error type to be returned.
    type(error_type), allocatable, intent(out) :: error
    !> Test matrix.
    class(rmatrix), allocatable :: A
    !> Krylov subspace.
    class(rvector), dimension(:), allocatable :: X
    !> Krylov subspace dimension.
    integer, parameter :: kdim = 3
    !> Hessenberg matrix.
    double precision, dimension(kdim+1, kdim) :: H
    !> Information flag.
    integer :: info
    !> Misc.
    double precision, dimension(3, 3) :: Adata
    double precision :: alpha
    integer :: i, j, k

    ! --> Initialize matrix with a 2-dimensional invariant subspace.
    A = rmatrix() ; call random_number(A%data)
    A%data(3, 1:2) = 0.0D+00 ; A%data(1:2, 3) = 0.0D+00
    ! --> Initialize Krylov subspace.
    allocate(X(1:kdim+1)) ; call random_number(X(1)%data)
    X(1)%data(3) = 0.0D+00 ! Makes sure X(1) has component only within the invariant subspace.
    alpha = X(1)%norm() ; call X(1)%scalar_mult(1.0D+00 / alpha)
    do k = 2, size(X)
       call X(k)%zero()
    enddo
    H = 0.0D+00
    ! --> Arnoldi factorization.
    call arnoldi_factorization(A, X, H, info)
    ! --> Check result.
    call check(error, info == 2)

    return
  end subroutine test_arnoldi_invariant_subspace

  subroutine test_arnoldi_basis_orthogonality(error)
    !> Error type to be returned.
    type(error_type), allocatable, intent(out) :: error
    !> Test matrix.
    class(rmatrix), allocatable :: A
    !> Krylov subspace.
    class(rvector), dimension(:), allocatable :: X
    !> Krylov subspace dimension.
    integer, parameter :: kdim = 2
    !> Hessenberg matrix.
    double precision, dimension(kdim+1, kdim) :: H
    !> Information flag.
    integer :: info
    !> Misc.
    double precision, dimension(3, 3) :: G, Id
    double precision :: alpha
    integer :: i, j, k

    ! --> Initialize random matrix.
    A = rmatrix() ; call random_number(A%data)
    ! --> Initialize Krylov subspace.
    allocate(X(1:kdim+1)) ; call random_number(X(1)%data)
    alpha = X(1)%norm() ; call X(1)%scalar_mult(1.0D+00 / alpha)
    do k = 2, size(X)
       call X(k)%zero()
    enddo
    H = 0.0D+00
    ! --> Arnoldi factorization.
    call arnoldi_factorization(A, X, H, info)
    ! --> Compute Gram matrix associated to the Krylov basis.
    do i = 1, kdim+1
       do j = 1, kdim+1
          G(i, j) = X(i)%dot(X(j))
       enddo
    enddo
    ! --> Check result.
    Id = reshape([1, 0, 0, 0, 1, 0, 0, 0, 1], shape=[3, 3]) !> Identity matrix.
    call check(error, norm2(G - Id) < 1e-12)

    return
  end subroutine test_arnoldi_basis_orthogonality

  !------------------------------------------------------------
  !-----                                                  -----
  !-----     TEST SUITE FOR THE LANCZOS FACTORIZATION     -----
  !-----                                                  -----
  !------------------------------------------------------------

  subroutine collect_lanczos_testsuite(testsuite)
    !> Collection of tests.
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [&
         new_unittest("Lanczos full factorization", test_lanczos_full_factorization),            &
         new_unittest("Lanczos basis orthogonality", test_lanczos_basis_orthogonality),          &
         new_unittest("Lanczos invariant subspace computation", test_lanczos_invariant_subspace) &
         ]
    return
  end subroutine collect_lanczos_testsuite

  subroutine test_lanczos_full_factorization(error)

    !> Error type to be returned.
    type(error_type), allocatable, intent(out) :: error
    !> Test matrix.
    class(spd_matrix), allocatable :: A
    double precision, dimension(3, 3) :: Adata
    !> Krylov subspace.
    class(rvector), dimension(:), allocatable :: X
    !> Krylov subspace dimension.
    integer, parameter :: kdim=3
    !> Tridiagonal matrix.
    double precision, dimension(kdim+1, kdim) :: T
    !> Information flag.
    integer :: info
    !> Misc.
    integer :: i, j, k
    double precision, dimension(3, kdim+1) :: Xdata
    double precision :: alpha
    double precision, dimension(3, 3) :: B, C

    ! --> Initialize matrix.
    A = spd_matrix()
    A%data(1, 1) = 1.0D+00    ; A%data(1, 2) = -0.5D+00 ; A%data(1, 3) = 0.0D+00
    A%data(2, 1) = -0.5D+00   ; A%data(2, 2) = 1.0D+00  ; A%data(2, 3) = 0.0D+00
    A%data(3, 1) = 0.0D+00    ; A%data(3, 2) = 0.0D+00  ; A%data(3, 3) = 1.0D+00
    Adata = A%data

    ! --> Initialize Krylov subspace.
    allocate(X(1:kdim+1)) ; call random_number(X(1)%data)
    X(1)%data = [1.0D+00, 2.0D+00, 1.0D+00]
    alpha = X(1)%norm() ; call X(1)%scalar_mult(1.0D+00 / alpha)
    do k = 2, size(X)
       call X(k)%zero()
    enddo
    T = 0.0D+00
    ! --> Lanczos factorization.
    call lanczos_factorization(A, X, T, info)
    ! --> Check correctness of full factorization.
    do k = 1, kdim+1
       Xdata(:, k) = X(k)%data
    enddo
    alpha = norm2( matmul(Adata, Xdata(:, 1:kdim)) - matmul(Xdata, T) )
    call check(error, alpha < 1e-12)

    return
  end subroutine test_lanczos_full_factorization

  subroutine test_lanczos_invariant_subspace(error)
    !> Error-type to be returned.
    type(error_type), allocatable, intent(out) :: error
    !> Test matrix.
    class(spd_matrix), allocatable :: A
    !> Krylov subspace.
    class(rvector), dimension(:), allocatable :: X
    !> Krylov subspace dimension.
    integer, parameter :: kdim=3
    !> Tridiagonal matrix.
    double precision, dimension(kdim+1, kdim) :: T
    !> Information flag.
    integer :: info
    !> Misc.
    double precision, dimension(3, 3) :: Adata
    double precision :: alpha
    integer :: i, j, k

    ! --> Initialize matrix with a 2-dimensional invariant subspace.
    A = spd_matrix()
    A%data(1, 1) = 1.0D+00    ; A%data(1, 2) = -0.5D+00 ; A%data(1, 3) = 0.0D+00
    A%data(2, 1) = -0.5D+00   ; A%data(2, 2) = 1.0D+00  ; A%data(2, 3) = 0.0D+00
    A%data(3, 1) = 0.0D+00    ; A%data(3, 2) = 0.0D+00  ; A%data(3, 3) = 1.0D+00

    ! --> Initialize Krylov subspace.
    allocate(X(1:kdim+1)) ; call random_number(X(1)%data)
    X(1)%data(3) = 0.0D+00 ! Makes sure X(1) has component only within the invariant subspace.
    alpha = X(1)%norm() ; call X(1)%scalar_mult(1.0D+00 / alpha)
    do k = 2, size(X)
       call X(k)%zero()
    enddo
    T = 0.0D+00
    ! --> Lanczos factorization.
    call lanczos_factorization(A, X, T, info)
    ! --> Check result.
    call check(error, info == 2)
    return
  end subroutine test_lanczos_invariant_subspace

  subroutine test_lanczos_basis_orthogonality(error)
    !> Error type to be returned.
    type(error_type), allocatable, intent(out) :: error
    !> Test matrix.
    class(spd_matrix), allocatable :: A
    !> Krylov subspace.
    class(rvector), dimension(:), allocatable :: X
    !> Krylov subspace dimension.
    integer, parameter :: kdim = 3
    !> Tridiagonal matrix.
    double precision, dimension(kdim+1, kdim) :: T
    !> Information flag.
    integer :: info
    !> Misc.
    double precision, dimension(kdim, kdim) :: G, Id
    double precision :: alpha
    integer :: i, j, k

    ! --> Initialize random spd matrix.
    A = spd_matrix() ; call random_number(A%data) ; A%data = matmul(A%data, transpose(A%data))
    ! --> Initialize Krylov subspace.
    allocate(X(1:kdim+1)) ; call random_number(X(1)%data)
    alpha = X(1)%norm() ; call X(1)%scalar_mult(1.0D+00 / alpha)
    do k = 2, size(X)
       call X(k)%zero()
    enddo
    T = 0.0D+00
    ! --> Lanczos factorization.
    call lanczos_factorization(A, X, T, info)
    ! --> Compute Gram matrix associated to the Krylov basis.
    do i = 1, kdim
       do j = 1, kdim
          G(i, j) = X(i)%dot(X(j))
       enddo
    enddo
    ! --> Check result.
    Id = reshape([1, 0, 0, 0, 1, 0, 0, 0, 1], shape=[3, 3])
    call check(error, norm2(G - Id) < 1e-12)
    return
  end subroutine test_lanczos_basis_orthogonality

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
