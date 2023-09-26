module TestKrylov
  use LightKrylov
  use TestVector
  use TestMatrices
  use testdrive  , only : new_unittest, unittest_type, error_type, check
  use stdlib_math, only : all_close
  implicit none

  private

  public :: collect_arnoldi_testsuite,              &
       collect_lanczos_tridiag_testsuite,      &
       collect_lanczos_bidiag_testsuite,       &
       collect_krylov_schur_testsuite

contains

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
    alpha = X(1)%norm() ; call X(1)%scal(1.0D+00 / alpha)
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
    alpha = X(1)%norm() ; call X(1)%scal(1.0D+00 / alpha)
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
    alpha = X(1)%norm() ; call X(1)%scal(1.0D+00 / alpha)
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

  !-----------------------------------------------------------------
  !-----                                                       -----
  !-----     TEST SUITE FOR THE LANCZOS TRIDIAGONALIZATION     -----
  !-----                                                       -----
  !-----------------------------------------------------------------

  subroutine collect_lanczos_tridiag_testsuite(testsuite)
    !> Collection of tests.
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [&
         new_unittest("Lanczos tridiag. full factorization", test_lanczos_tridiag_full_factorization),            &
         new_unittest("Lanczos tridiag. basis orthogonality", test_lanczos_tridiag_basis_orthogonality),          &
         new_unittest("Lanczos tridiag. invariant subspace computation", test_lanczos_tridiag_invariant_subspace) &
         ]
    return
  end subroutine collect_lanczos_tridiag_testsuite

  subroutine test_lanczos_tridiag_full_factorization(error)

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
    alpha = X(1)%norm() ; call X(1)%scal(1.0D+00 / alpha)
    do k = 2, size(X)
       call X(k)%zero()
    enddo
    T = 0.0D+00
    ! --> Lanczos factorization.
    call lanczos_tridiagonalization(A, X, T, info)
    ! --> Check correctness of full factorization.
    do k = 1, kdim+1
       Xdata(:, k) = X(k)%data
    enddo
    alpha = norm2( matmul(Adata, Xdata(:, 1:kdim)) - matmul(Xdata, T) )
    call check(error, alpha < 1e-12)

    return
  end subroutine test_lanczos_tridiag_full_factorization

  subroutine test_lanczos_tridiag_invariant_subspace(error)
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
    alpha = X(1)%norm() ; call X(1)%scal(1.0D+00 / alpha)
    do k = 2, size(X)
       call X(k)%zero()
    enddo
    T = 0.0D+00
    ! --> Lanczos factorization.
    call lanczos_tridiagonalization(A, X, T, info)
    ! --> Check result.
    call check(error, info == 2)
    return
  end subroutine test_lanczos_tridiag_invariant_subspace

  subroutine test_lanczos_tridiag_basis_orthogonality(error)
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
    alpha = X(1)%norm() ; call X(1)%scal(1.0D+00 / alpha)
    do k = 2, size(X)
       call X(k)%zero()
    enddo
    T = 0.0D+00
    ! --> Lanczos factorization.
    call lanczos_tridiagonalization(A, X, T, info)
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
  end subroutine test_lanczos_tridiag_basis_orthogonality

  !----------------------------------------------------------------
  !-----                                                      -----
  !-----     TEST SUITE FOR THE LANCZOS BIDIAGONALIZATION     -----
  !-----                                                      -----
  !----------------------------------------------------------------

  subroutine collect_lanczos_bidiag_testsuite(testsuite)
    !> Collection of tests.
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [&
         new_unittest("Lanczos bidiag. full factorization", test_lanczos_bidiag_full_factorization) &
         ]

    return
  end subroutine collect_lanczos_bidiag_testsuite

  subroutine test_lanczos_bidiag_full_factorization(error)
    !> Error type to be returned.
    type(error_type), allocatable, intent(out) :: error
    !> Test matrix.
    class(rmatrix), allocatable :: A
    !> Left and right Krylov subspaces.
    class(rvector), allocatable :: U(:)
    class(rvector), allocatable :: V(:)
    !> Krylov subspace dimension.
    integer, parameter :: kdim = 3
    !> Bidiagonal matrix.
    real(kind=wp) :: B(kdim+1, kdim)
    !> Information flag.
    integer :: info
    !> Miscellaneous.
    integer :: i, j, k
    real(kind=wp) :: alpha
    real(kind=wp) :: Udata(3, kdim+1), Vdata(3, kdim+1)

    ! --> Initialize matrix.
    A = rmatrix() ; call random_number(A%data)

    ! --> Initialize Krylov subspaces.
    allocate(U(1:kdim+1)) ; allocate(V(1:kdim+1))
    do k = 1, size(U)
       call U(k)%zero() ; call V(k)%zero()
    enddo
    call random_number(U(1)%data)
    alpha = U(1)%norm() ; call U(1)%scal(1.0_wp / alpha)
    B = 0.0_wp
    ! --> Lanczos bidiagonalization.
    call lanczos_bidiagonalization(A, U, V, B, info)
    ! --> Check correctness of full factorization.
    do k = 1, size(U)
       Udata(:, k) = U(k)%data ; Vdata(:, k) = V(k)%data
    enddo
    alpha = norm2(matmul(A%data, Vdata(:, 1:kdim)) - matmul(Udata, B))
    call check(error, alpha < 1e-12)
    return
  end subroutine test_lanczos_bidiag_full_factorization

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
