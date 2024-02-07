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
            collect_nonsymmetric_lanczos_testsuite, &
            !collect_rational_arnoldi_testsuite,     &
            collect_two_sided_arnoldi_testsuite,    &
            collect_qr_testsuite

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
         new_unittest("Arnoldi full factorization", test_arnoldi_full_factorization),  &
         new_unittest("Arnoldi basis orthogonality", test_arnoldi_basis_orthogonality) &
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
    !> Krylov subspace.
    class(rvector), dimension(:), allocatable :: X
    !> Krylov subspace dimension.
    integer, parameter :: kdim=test_size
    !> Hessenberg matrix.
    real(kind=wp) :: H(kdim+1, kdim)
    !> Information flag.
    integer :: info
    !> Misc.
    integer :: k
    real(kind=wp) :: Xdata(test_size, kdim+1)
    real(kind=wp) :: alpha

    ! --> Initialize matrix.
    A = rmatrix() ; call random_number(A%data)
    ! --> Initialize Krylov subspace.
    allocate(X(1:kdim+1)) ; call random_number(X(1)%data)
    alpha = X(1)%norm() ; call X(1)%scal(1.0D+00 / alpha)
    do k = 2, size(X)
       call X(k)%zero()
    enddo
    H = 0.0_wp
    ! --> Arnoldi factorization.
    call arnoldi_factorization(A, X, H, info)
    ! --> Check correctness of full factorization.
    do k = 1, kdim+1
       Xdata(:, k) = X(k)%data
    enddo
    call check(error, all_close(matmul(A%data, Xdata(:, 1:kdim)),  matmul(Xdata, H), rtol, atol) )

    return
  end subroutine test_arnoldi_full_factorization

  subroutine test_arnoldi_basis_orthogonality(error)
    !> Error type to be returned.
    type(error_type), allocatable, intent(out) :: error
    !> Test matrix.
    class(rmatrix), allocatable :: A
    !> Krylov subspace.
    class(rvector), dimension(:), allocatable :: X
    !> Krylov subspace dimension.
    integer, parameter :: kdim = test_size
    !> Hessenberg matrix.
    double precision, dimension(kdim+1, kdim) :: H
    !> Information flag.
    integer :: info
    !> Misc.
    double precision, dimension(kdim, kdim) :: G, Id
    double precision :: alpha
    integer :: i, j, k

    ! --> Initialize random matrix.
    A = rmatrix() ; call random_number(A%data)
    ! --> Initialize Krylov subspace.
    allocate(X(1:kdim+1)) ; call random_number(X(1)%data)
    alpha = X(1)%norm() ; call X(1)%scal(1.0_wp / alpha)
    do k = 2, size(X)
       call X(k)%zero()
    enddo
    H = 0.0_wp
    ! --> Arnoldi factorization.
    call arnoldi_factorization(A, X, H, info)
    ! --> Compute Gram matrix associated to the Krylov basis.
    do i = 1, kdim
       do j = 1, kdim
          G(i, j) = X(i)%dot(X(j))
       enddo
    enddo

    ! --> Check result.
    Id = 0.0_wp ; forall(i=1:kdim) Id(i, i) = 1.0_wp
    call check(error, norm2(G - Id) < rtol)

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
         new_unittest("Lanczos tridiag. full factorization", test_lanczos_tridiag_full_factorization),  &
         new_unittest("Lanczos tridiag. basis orthogonality", test_lanczos_tridiag_basis_orthogonality) &
         ]
    return
  end subroutine collect_lanczos_tridiag_testsuite

  subroutine test_lanczos_tridiag_full_factorization(error)

    !> Error type to be returned.
    type(error_type), allocatable, intent(out) :: error
    !> Test matrix.
    class(spd_matrix), allocatable :: A
    !> Krylov subspace.
    class(rvector), dimension(:), allocatable :: X
    !> Krylov subspace dimension.
    integer, parameter :: kdim=test_size
    !> Tridiagonal matrix.
    real(kind=wp), dimension(kdim+1, kdim) :: T
    !> Information flag.
    integer :: info
    !> Misc.
    integer :: i, j, k
    real(kind=wp), dimension(test_size, kdim+1) :: Xdata
    real(kind=wp) :: alpha

    ! --> Initialize matrix.
    A = spd_matrix() ; call random_number(A%data) ; A%data = matmul(A%data, transpose(A%data))

    ! --> Initialize Krylov subspace.
    allocate(X(1:kdim+1)) ; call random_number(X(1)%data)
    alpha = X(1)%norm() ; call X(1)%scal(1.0D+00 / alpha)
    do k = 2, size(X)
       call X(k)%zero()
    enddo
    T = 0.0_wp
    ! --> Lanczos factorization.
    call lanczos_tridiagonalization(A, X, T, info)
    ! --> Check correctness of full factorization.
    do k = 1, kdim+1
       Xdata(:, k) = X(k)%data
    enddo
    ! --> Infinity-norm check.
    alpha = maxval(abs(matmul(A%data, Xdata(:, 1:kdim)) - matmul(Xdata, T)))
    write(*, *) "Infinity-norm      :", alpha
    write(*, *) "Relative tolerance :", rtol
    call check(error, alpha < rtol)

    return
  end subroutine test_lanczos_tridiag_full_factorization

  subroutine test_lanczos_tridiag_basis_orthogonality(error)
    !> Error type to be returned.
    type(error_type), allocatable, intent(out) :: error
    !> Test matrix.
    class(spd_matrix), allocatable :: A
    !> Krylov subspace.
    class(rvector), dimension(:), allocatable :: X
    !> Krylov subspace dimension.
    integer, parameter :: kdim = test_size
    !> Tridiagonal matrix.
    double precision, dimension(kdim+1, kdim) :: T
    !> Information flag.
    integer :: info
    !> Misc.
    double precision, dimension(kdim+1, kdim+1) :: G, Id
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
    T = 0.0_wp
    ! --> Lanczos factorization.
    call lanczos_tridiagonalization(A, X, T, info)
    ! --> Compute Gram matrix associated to the Krylov basis.
    do i = 1, kdim
       do j = 1, kdim
          G(i, j) = X(i)%dot(X(j))
       enddo
    enddo
    ! --> Check result.
    Id = 0.0_wp ; forall(i=1:kdim) Id(i, i) = 1.0_wp
    call check(error, norm2(G - Id) < rtol)
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
    integer, parameter :: kdim = test_size
    !> Bidiagonal matrix.
    real(kind=wp) :: B(kdim+1, kdim)
    !> Information flag.
    integer :: info
    !> Miscellaneous.
    integer :: i, j, k
    real(kind=wp) :: alpha
    real(kind=wp) :: Udata(test_size, kdim+1), Vdata(test_size, kdim+1)

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
    ! --> Infinity-norm check.
    alpha = maxval(abs(matmul(A%data, Vdata(:, 1:kdim)) - matmul(Udata, B)))
    write(*, *) "Infinity norm      :", alpha
    write(*, *) "Relative tolerance :", rtol
    call check(error, alpha < rtol)
    return
  end subroutine test_lanczos_bidiag_full_factorization

  !--------------------------------------------------------------------------
  !-----                                                                -----
  !-----     TEST SUITE FOR THE NON-SYMMETRIC LANCZOS FACTORIZATION     -----
  !-----                                                                -----
  !--------------------------------------------------------------------------

  subroutine collect_nonsymmetric_lanczos_testsuite(testsuite)
    !> Collection of tests.
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
         new_unittest("Non-symmetric Lanczos tridiag. full factorization", test_nonsym_lanczos_full_factorization) &
         ]

    return
  end subroutine collect_nonsymmetric_lanczos_testsuite

  subroutine test_nonsym_lanczos_full_factorization(error)
    !> Error type to be returned.
    type(error_type), allocatable, intent(out) :: error
    !> Test matrix.
    class(rmatrix), allocatable :: A
    !> Left and right Krylov subspaces.
    class(rvector), allocatable :: V(:), W(:)
    !> Krylov subspace dimenion.
    integer, parameter :: kdim = test_size
    !> Tridiagonal matrix.
    real(kind=wp) :: T(kdim+1, kdim+1)
    !> Information flag.
    integer :: info
    !> Miscellaneous.
    integer :: i, j, k
    real(kind=wp) :: alpha, beta
    real(kind=wp) :: Vdata(test_size, kdim+1), Wdata(test_size, kdim+1)

    ! --> Initialize matrix.
    A = rmatrix() ; call random_number(A%data)

    ! --> Initialize Krylov subspaces.
    allocate(V(1:kdim+1)) ; allocate(W(1:kdim+1))
    call random_number(V(1)%data) ; call random_number(W(1)%data)
    alpha = V(1)%norm() ; call V(1)%scal(1.0_wp / alpha)
    alpha = W(1)%norm() ; call W(1)%scal(1.0_wp / alpha)
    do k = 2, size(V)
       call V(k)%zero() ; call W(k)%zero()
    enddo
    T = 0.0_wp

    ! --> Nonsymmetric Lanczos factorization.
    call nonsymmetric_lanczos_tridiagonalization(A, V, W, T, info, verbosity=.true.)

    ! --> Check correctness of the factorization.
    do k = 1, size(V)
       Vdata(:, k) = V(k)%data ; Wdata(:, k) = W(k)%data
    enddo

    ! --> Infinity-norm check.
    alpha = maxval( abs(matmul(A%data, Vdata(:, 1:kdim)) - matmul(Vdata, T(1:kdim+1, 1:kdim))) )
    write(*, *) "Infinity norm      :", alpha
    write(*, *) "Relative tolerance :", rtol
    call check(error, alpha < rtol)

    return
  end subroutine test_nonsym_lanczos_full_factorization

  !----------------------------------------------------------------------
  !-----                                                            -----
  !-----     TEST SUITE FOR THE TWO-SIDED ARNOLDI FACTORIZATION     -----
  !-----                                                            -----
  !----------------------------------------------------------------------

  subroutine collect_two_sided_arnoldi_testsuite(testsuite)
    !> Collection of tests.
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [&
         new_unittest("Two-sided Arnoldi full factorization", test_two_sided_arnoldi_full_factorization),              &
         new_unittest("Two-sided Arnoldi full factorization adjoint", test_two_sided_arnoldi_full_factorization_bis),  &
         new_unittest("Two-sided Arnoldi basis orthogonality", test_two_sided_arnoldi_basis_orthogonality),            &
         new_unittest("Two-sided Arnoldi basis orthogonality adjoint", test_two_sided_arnoldi_basis_orthogonality_bis) &
         ]

    return
  end subroutine collect_two_sided_arnoldi_testsuite

  subroutine test_two_sided_arnoldi_full_factorization(error)
    !> Error type to be returned.
    type(error_type), allocatable, intent(out) :: error
    !> Test matrix.
    class(rmatrix), allocatable :: A
    !> Krylov subspaces.
    class(rvector), allocatable :: V(:), W(:)
    !> Krylov subspace dimension.
    integer, parameter :: kdim=test_size
    !> Hessenberg matrices.
    real(kind=wp) :: H(kdim+1, kdim), G(kdim+1, kdim)
    !> Information flag.
    integer :: info
    !> Miscellaneous.
    integer :: k
    real(kind=wp) :: Vdata(test_size, kdim+1), Wdata(test_size, kdim+1)
    real(kind=wp) :: alpha

    !> Initialize matrix.
    A = rmatrix() ; call random_number(A%data)
    !> Initialize Krylov subspaces.
    allocate(V(1:kdim+1)) ; allocate(W(1:kdim+1))
    call random_number(V(1)%data) ; call random_number(W(1)%data)
    alpha = V(1)%norm() ; call V(1)%scal(1.0_wp / alpha)
    alpha = W(1)%norm() ; call W(1)%scal(1.0_wp / alpha)
    do k = 2, size(V)
       call V(k)%zero() ; call W(k)%zero()
    enddo

    !> Initialize Hessenberg matrices.
    H = 0.0_wp ; G = 0.0_wp

    !> Two-sided Arnoldi factoriztion.
    call two_sided_arnoldi_factorization(A, V, W, H, G, info)

    !> Check correctness of the full factorization.
    do k = 1, kdim+1
       Vdata(:, k) = V(k)%data
    enddo
    call check(error, all_close(matmul(A%data, Vdata(:, 1:kdim)),  matmul(Vdata, H), rtol, atol) )

    return
  end subroutine test_two_sided_arnoldi_full_factorization

  subroutine test_two_sided_arnoldi_full_factorization_bis(error)
    !> Error type to be returned.
    type(error_type), allocatable, intent(out) :: error
    !> Test matrix.
    class(rmatrix), allocatable :: A
    !> Krylov subspaces.
    class(rvector), allocatable :: V(:), W(:)
    !> Krylov subspace dimension.
    integer, parameter :: kdim=test_size
    !> Hessenberg matrices.
    real(kind=wp) :: H(kdim+1, kdim), G(kdim+1, kdim)
    !> Information flag.
    integer :: info
    !> Miscellaneous.
    integer :: k
    real(kind=wp) :: Vdata(test_size, kdim+1), Wdata(test_size, kdim+1)
    real(kind=wp) :: alpha

    !> Initialize matrix.
    A = rmatrix() ; call random_number(A%data)
    !> Initialize Krylov subspaces.
    allocate(V(1:kdim+1)) ; allocate(W(1:kdim+1))
    call random_number(V(1)%data) ; call random_number(W(1)%data)
    alpha = V(1)%norm() ; call V(1)%scal(1.0_wp / alpha)
    alpha = W(1)%norm() ; call W(1)%scal(1.0_wp / alpha)
    do k = 2, size(V)
       call V(k)%zero() ; call W(k)%zero()
    enddo

    !> Initialize Hessenberg matrices.
    H = 0.0_wp ; G = 0.0_wp

    !> Two-sided Arnoldi factoriztion.
    call two_sided_arnoldi_factorization(A, V, W, H, G, info)

    !> Check correctness of the full factorization.
    do k = 1, kdim+1
       Wdata(:, k) = W(k)%data
    enddo
    call check(error, all_close(matmul(transpose(A%data), Wdata(:, 1:kdim)),  matmul(Wdata, G), rtol, atol) )

    return
  end subroutine test_two_sided_arnoldi_full_factorization_bis

   subroutine test_two_sided_arnoldi_basis_orthogonality(error)
    !> Error type to be returned.
    type(error_type), allocatable, intent(out) :: error
    !> Test matrix.
    class(rmatrix), allocatable :: A
    !> Krylov subspaces.
    class(rvector), allocatable :: V(:), W(:)
    !> Krylov subspace dimension.
    integer, parameter :: kdim=test_size
    !> Hessenberg matrices.
    real(kind=wp) :: H(kdim+1, kdim), G(kdim+1, kdim)
    !> Information flag.
    integer :: info
    !> Miscellaneous.
    real(kind=wp) :: M(kdim, kdim), Id(kdim, kdim)
    integer :: i, j, k
    real(kind=wp) :: alpha

    !> Initialize matrix.
    A = rmatrix() ; call random_number(A%data)
    !> Initialize Krylov subspaces.
    allocate(V(1:kdim+1)) ; allocate(W(1:kdim+1))
    call random_number(V(1)%data) ; call random_number(W(1)%data)
    alpha = V(1)%norm() ; call V(1)%scal(1.0_wp / alpha)
    alpha = W(1)%norm() ; call W(1)%scal(1.0_wp / alpha)
    do k = 2, size(V)
       call V(k)%zero() ; call W(k)%zero()
    enddo

    !> Initialize Hessenberg matrices.
    H = 0.0_wp ; G = 0.0_wp

    !> Two-sided Arnoldi factoriztion.
    call two_sided_arnoldi_factorization(A, V, W, H, G, info)

    !> Inner-product matrix.
    M = 0.0_wp
    do i = 1, kdim
       do j = 1, kdim
          M(i, j) = V(i)%dot(V(j))
       enddo
    enddo

    !> Check results.
    Id = 0.0_wp ; forall (i=1:kdim) Id(i, i) = 1.0_wp
    call check(error, norm2(M-Id) < rtol)

    return
  end subroutine test_two_sided_arnoldi_basis_orthogonality

 subroutine test_two_sided_arnoldi_basis_orthogonality_bis(error)
    !> Error type to be returned.
    type(error_type), allocatable, intent(out) :: error
    !> Test matrix.
    class(rmatrix), allocatable :: A
    !> Krylov subspaces.
    class(rvector), allocatable :: V(:), W(:)
    !> Krylov subspace dimension.
    integer, parameter :: kdim=test_size
    !> Hessenberg matrices.
    real(kind=wp) :: H(kdim+1, kdim), G(kdim+1, kdim)
    !> Information flag.
    integer :: info
    !> Miscellaneous.
    real(kind=wp) :: M(kdim, kdim), Id(kdim, kdim)
    integer :: i, j, k
    real(kind=wp) :: alpha

    !> Initialize matrix.
    A = rmatrix() ; call random_number(A%data)
    !> Initialize Krylov subspaces.
    allocate(V(1:kdim+1)) ; allocate(W(1:kdim+1))
    call random_number(V(1)%data) ; call random_number(W(1)%data)
    alpha = V(1)%norm() ; call V(1)%scal(1.0_wp / alpha)
    alpha = W(1)%norm() ; call W(1)%scal(1.0_wp / alpha)
    do k = 2, size(V)
       call V(k)%zero() ; call W(k)%zero()
    enddo

    !> Initialize Hessenberg matrices.
    H = 0.0_wp ; G = 0.0_wp

    !> Two-sided Arnoldi factoriztion.
    call two_sided_arnoldi_factorization(A, V, W, H, G, info)

    !> Inner-product matrix.
    M = 0.0_wp
    do i = 1, kdim
       do j = 1, kdim
          M(i, j) = W(i)%dot(W(j))
       enddo
    enddo

    !> Check results.
    Id = 0.0_wp ; forall (i=1:kdim) Id(i, i) = 1.0_wp
    call check(error, norm2(M-Id) < rtol)

    return
  end subroutine test_two_sided_arnoldi_basis_orthogonality_bis

  !-------------------------------------------------------
  !-----                                             -----
  !-----     TEST SUITE FOR THE QR FACTORIZATION     -----
  !-----                                             -----
  !-------------------------------------------------------

   subroutine collect_qr_testsuite(testsuite)
      !> Collection of tests.
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [&
        new_unittest("QR factorization", test_qr_factorization),  &
        new_unittest("Q orthonormality", test_qr_basis_orthonormality)  &
        ]

      return
   end subroutine collect_qr_testsuite

   subroutine test_qr_factorization(error)
     ! This function checks the correctness of the QR implementation by
     ! verifying that the factorization is correct, i.e. A = Q @ R.
     ! A random matrix is used for testing.

     !> Error type to be returned.
     type(error_type), allocatable, intent(out) :: error
     !> Test matrix.
     class(rvector), dimension(:), allocatable :: A
     !> Basis vectors.
     class(rvector), dimension(:), allocatable :: Q
     !> Krylov subspace dimension.
     integer, parameter :: kdim = 3
     !> GS factors.
     real(kind=wp) :: R(kdim, kdim)
     !> Information flag.
     integer :: info
     !> Misc.
     integer :: i,j,k
     real(kind=wp) :: Amat(test_size, kdim), Qmat(test_size, kdim)
     real(kind=wp) :: alpha

     ! --> Initialize matrix.
     allocate(A(1:kdim));
     do k = 1, size(A)
        call random_number(A(k)%data)
     enddo
     ! --> Initialize Krylov subspace.
     allocate(Q(1:kdim));
     do k = 1, size(Q)
        call Q(k)%zero()
     enddo
     R = 0.0_wp
     ! --> QR factorization.
     call qr_factorization(A, Q, R, info)
     ! --> Extract data 
     do k = 1, kdim
        Amat(:, k) = A(k)%data
        Qmat(:, k) = Q(k)%data
     enddo
     ! --> Check correctness of QR factorization.
     call check(error, all_close(Amat, matmul(Qmat, R), rtol, atol) )

     return
   end subroutine test_qr_factorization

   subroutine test_qr_basis_orthonormality(error)
     ! This function checks the correctness of the QR implementation by
     ! verifying that the obtained basis is orthonormal, i.e. Q.T @ Q = I.
     ! A random matrix is used for testing.
      
     !> Error type to be returned.
     type(error_type), allocatable, intent(out) :: error
     !> Test matrix.
     class(rvector), dimension(:), allocatable :: A
     !> Basis vectors.
     class(rvector), dimension(:), allocatable :: Q
     !> Krylov subspace dimension.
     integer, parameter :: kdim = 3
     !> GS factors.
     real(kind=wp) :: R(kdim, kdim)
     real(kind=wp) :: Id(kdim, kdim)
     !> Information flag.
     integer :: info
     !> Misc.
     integer :: i,k
     real(kind=wp) :: Qtmat(kdim, test_size), Qmat(test_size, kdim)
     real(kind=wp) :: alpha

     ! --> Initialize matrix.
     allocate(A(1:kdim));
     do k = 1, size(A)
        call random_number(A(k)%data)
     enddo
     ! --> Initialize Krylov subspace.
     allocate(Q(1:kdim));
     do k = 1, size(Q)
        call Q(k)%zero()
     enddo
     R = 0.0_wp
     ! --> QR factorization.
     call qr_factorization(A, Q, R, info)
     ! --> Extract data 
     do k = 1, kdim
        Qmat(:, k) = Q(k)%data
        do i = 1, test_size
           Qtmat(k,i) = Qmat(i,k)
        enddo
     enddo
     ! --> Identity
     Id = 0.0_wp
     do k = 1, kdim
        Id(k,k) = 1.0_wp
     enddo
     ! --> Check correctness of QR factorization.
     call check(error, all_close(Id, matmul(Qtmat, Qmat), rtol, atol) )

   end subroutine test_qr_basis_orthonormality

end module TestKrylov
