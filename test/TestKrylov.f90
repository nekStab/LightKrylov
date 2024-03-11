module TestKrylov
   use LightKrylov
   use TestVector
   use TestMatrices
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use stdlib_math, only: all_close, is_close
   implicit none

   private

   public :: collect_arnoldi_testsuite, &
             collect_lanczos_tridiag_testsuite, &
             collect_lanczos_bidiag_testsuite, &
             collect_nonsymmetric_lanczos_testsuite, &
             !collect_rational_arnoldi_testsuite,     &
             collect_two_sided_arnoldi_testsuite, &
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

      testsuite = [ &
           new_unittest("Arnoldi full factorization", test_arnoldi_full_factorization), &
           new_unittest("Arnoldi basis orthogonality", test_arnoldi_basis_orthogonality), &
           new_unittest("Block Arnoldi full factorization", test_block_arnoldi_full_factorization), &
           new_unittest("Block Arnoldi basis orthogonality", test_block_arnoldi_basis_orthogonality), &
           new_unittest("Krylov-Schur restart", test_krylov_schur) &
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
      integer, parameter :: kdim = test_size
      !> Hessenberg matrix.
      real(kind=wp) :: H(kdim + 1, kdim)
      !> Information flag.
      integer :: info
      !> Misc.
      integer :: k
      real(kind=wp) :: Xdata(test_size, kdim + 1)
      real(kind=wp) :: alpha
      class(rvector), allocatable :: X0(1)

      ! --> Initialize matrix.
      A = rmatrix(); call random_number(A%data)
      ! --> Initialize Krylov subspace.
      allocate (X(1:kdim + 1)); allocate (X0(1))
      call random_number(X0(1)%data);
      call initialize_krylov_subspace(X, X0)
      H = 0.0_wp

      ! --> Arnoldi factorization.
      call arnoldi_factorization(A, X, H, info)

      ! --> Check correctness of full factorization.
      do k = 1, kdim + 1
         Xdata(:, k) = X(k)%data
      end do
      call check(error, all_close(matmul(A%data, Xdata(:, 1:kdim)), matmul(Xdata, H), rtol, atol))

      return
   end subroutine test_arnoldi_full_factorization

   subroutine test_arnoldi_basis_orthogonality(error)
      !> Error type to be returned.
      type(error_type), allocatable, intent(out) :: error
      !> Test matrix.
      class(rmatrix), allocatable :: A
      !> Krylov subspace.
      class(rvector), dimension(:), allocatable :: X
      class(rvector), dimension(:), allocatable :: X0
      !> Krylov subspace dimension.
      integer, parameter :: kdim = test_size
      !> Hessenberg matrix.
      double precision, dimension(kdim + 1, kdim) :: H
      !> Information flag.
      integer :: info
      !> Misc.
      double precision, dimension(kdim, kdim) :: G, Id
      double precision :: alpha
      integer :: i, j, k

      ! --> Initialize random matrix.
      A = rmatrix(); call random_number(A%data)
      ! --> Initialize Krylov subspace.
      allocate (X(1:kdim + 1)); allocate (X0(1))
      call random_number(X0(1)%data);
      call initialize_krylov_subspace(X, X0)
      H = 0.0_wp

      ! --> Arnoldi factorization.
      call arnoldi_factorization(A, X, H, info)
      
      ! --> Compute Gram matrix associated to the Krylov basis.
      G = 0.0_wp
      call mat_mult(G,X(1:kdim),X(1:kdim))
      
      ! --> Check result.
      Id = 0.0_wp; forall (i=1:kdim) Id(i, i) = 1.0_wp
      call check(error, norm2(G - Id) < rtol)

      return
   end subroutine test_arnoldi_basis_orthogonality

   subroutine test_block_arnoldi_full_factorization(error)
      ! This function checks the correctness of the block Arnoldi implementation by
      ! verifying that the full factorization is correct, i.e. A @ X[1:p*k] = X[1:p*(k+1)] @ H.

      !> Error type to be returned.
      type(error_type), allocatable, intent(out) :: error
      !> Test matrix.
      class(rmatrix), allocatable :: A
      !> Krylov subspace.
      class(rvector), dimension(:), allocatable :: X
      !> Krylov subspace dimension.
      integer, parameter :: p = 2
      integer, parameter :: kdim = test_size/p
      !> Hessenberg matrix.
      real(kind=wp) :: H(p*(kdim + 1), p*kdim)
      !> Information flag.
      integer :: info
      !> Misc.
      integer :: k
      real(kind=wp) :: Xdata(test_size, p*(kdim + 1))
      real(kind=wp) :: Rwrk(p,p)
      double precision, dimension(p*kdim, p*kdim) :: G, Id
      class(rvector), dimension(:), allocatable :: X0
      
      ! --> Initialize matrix.
      A = rmatrix(); call random_number(A%data)
      ! --> Initialize Krylov subspace.
      allocate (X(1:p*(kdim + 1))); allocate (X0(1:p)); 
      do k = 1,p
         call random_number(X0(k)%data)
      enddo
      call initialize_krylov_subspace(X, X0)
      H = 0.0_wp

      ! --> Arnoldi factorization.
      call arnoldi_factorization(A, X, H, info, block_size = p) 
      
      ! --> Check correctness of full factorization.
      do k = 1, p*(kdim + 1)
         Xdata(:, k) = X(k)%data
      end do
      call check(error, all_close(matmul(A%data, Xdata(:, 1:p*kdim)), matmul(Xdata, H), rtol, atol))

      return
   end subroutine test_block_arnoldi_full_factorization

   subroutine test_block_arnoldi_basis_orthogonality(error)
      !> Error type to be returned.
      type(error_type), allocatable, intent(out) :: error
      !> Test matrix.
      class(rmatrix), allocatable :: A
      !> Krylov subspace.
      class(rvector), dimension(:), allocatable :: X
      class(rvector), dimension(:), allocatable :: X0
      !> Krylov subspace dimension.
      integer, parameter :: p = 2
      integer, parameter :: kdim = test_size/p
      !> Hessenberg matrix.
      double precision, dimension(p*(kdim + 1), p*kdim) :: H
      !> Information flag.
      integer :: info
      !> Misc.
      double precision, dimension(p*kdim, p*kdim) :: G, Id
      double precision :: Rwrk(p,p), alpha
      integer :: i, j, k

      !! --> Initialize random matrix.
      A = rmatrix(); call random_number(A%data)
      !! --> Initialize Krylov subspace.
      allocate (X(1:p*(kdim + 1))); allocate (X0(1:p)); 
      do k = 1,p
         call random_number(X0(k)%data)
      enddo
      call initialize_krylov_subspace(X, X0)
      H = 0.0_wp

      ! --> Arnoldi factorization.
      call arnoldi_factorization(A, X, H, info, block_size = p)

      ! --> Compute Gram matrix associated to the Krylov basis.
      G = 0.0_wp
      call mat_mult(G,X(1:p*kdim),X(1:p*kdim))

      ! --> Check result.
      Id = 0.0_wp; forall (i=1:p*kdim) Id(i, i) = 1.0_wp
      call check(error, norm2(G - Id) < rtol)

      return
   end subroutine test_block_arnoldi_basis_orthogonality

   subroutine test_krylov_schur(error)
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
     integer, parameter :: kdim = 10
     !> Hessenberg matrix.
     real(kind=wp) :: H(kdim + 1, kdim)
     !> Information flag.
     integer :: info
     !> Misc.
     integer :: k, nblk
     real(kind=wp) :: Xdata(test_size, kdim + 1)
     real(kind=wp) :: alpha
     
     ! --> Initialize matrix.
     A = rmatrix(); call random_number(A%data)
     ! --> Initialize Krylov subspace.
     allocate (X(1:kdim + 1)) ; call initialize_krylov_subspace(X)
     call random_number(X(1)%data); alpha = X(1)%norm() ; call X(1)%scal(1.0_wp / alpha)
     H = 0.0_wp
     
     ! --> Arnoldi factorization.
     call arnoldi_factorization(A, X, H, info)
     
     ! --> Krylov-Schur condensation.
     call krylov_schur_restart(nblk, X, H, select_eigs)
     
     ! --> Check correctness of full factorization.
     do k = 1, kdim + 1
        Xdata(:, k) = X(k)%data
     end do
     
     !> Infinity-norm of the error.
     alpha = maxval(abs(matmul(A%data, Xdata(:, 1:nblk)) - matmul(Xdata(:, 1:nblk+1), H(1:nblk+1, 1:nblk))))
     
     !> Check correctness.
     call check(error, alpha < rtol)
     
     return
   contains
     function select_eigs(eigvals) result(selected)
       complex(kind=wp), intent(in) :: eigvals(:)
       logical                      :: selected(size(eigvals))
       selected = (eigvals%re > 0.5_wp)
       return
     end function select_eigs
   end subroutine test_krylov_schur
   
   subroutine test_block_krylov_schur(error)
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
     integer, parameter :: p = 3
     integer, parameter :: kdim = 7
     !> Hessenberg matrix.
     real(kind=wp) :: H(p*(kdim + 1), p*kdim)
     !> Information flag.
     integer :: info
     !> Misc.
     integer :: k, nblk
     real(kind=wp) :: Xdata(test_size, p*(kdim + 1))
     real(kind=wp) :: alpha, m(test_size)
     class(rvector), allocatable :: X0(:)
     
     ! --> Initialize matrix.
     A = rmatrix(); call random_number(A%data) ; m = sum(A%data, 2)/test_size
     do k = 1, test_size
        A%data(:, k) = A%data(:, k) - m
     enddo
     !A%data = 0.5 * (A%data + transpose(A%data))
     ! --> Initialize Krylov subspace.
     allocate (X(p*(kdim + 1))) ; allocate(X0(p))
     do k = 1, p
        call random_number(X0(k)%data)
     enddo
     call initialize_krylov_subspace(X, X0)
     H = 0.0_wp
     
     ! --> Arnoldi factorization.
     call arnoldi_factorization(A, X, H, info, block_size=p)
     
     ! --> Krylov-Schur condensation.
     call krylov_schur_restart(nblk, X, H, select_eigs, p)
     
     ! --> Check correctness of full factorization.
     do k = 1, p*(kdim + 1)
        Xdata(:, k) = X(k)%data
     end do
     
     !> Infinity-norm of the error.
     alpha = maxval(abs(matmul(A%data, Xdata(:, 1:p*nblk)) - matmul(Xdata(:, 1:p*(nblk+1)), H(1:p*(nblk+1), 1:p*nblk))))
     
     !> Check correctness.
     call check(error, alpha < rtol)
     
     return
   contains
     function select_eigs(eigvals) result(selected)
       complex(kind=wp), intent(in) :: eigvals(:)
       logical                      :: selected(size(eigvals))
       selected = (abs(eigvals) > 2.0_wp)
       return
     end function select_eigs
   end subroutine test_block_krylov_schur
   
   subroutine test_krylov_schur_basis_orthogonality(error)
     !> Error type to be returned.
     type(error_type), allocatable, intent(out) :: error
     !> Test matrix.
     class(rmatrix), allocatable :: A
     !> Krylov subspace.
     class(rvector), dimension(:), allocatable :: X
     class(rvector), dimension(:), allocatable :: X0
     !> Krylov subspace dimension.
     integer, parameter :: kdim = 10
     !> Hessenberg matrix.
     double precision, dimension(kdim + 1, kdim) :: H
     !> Information flag.
     integer :: info
     !> Misc.
     double precision, dimension(kdim, kdim) :: G, Id
     double precision :: alpha
     integer :: i, j, k, nblk
     
     ! --> Initialize random matrix.
     A = rmatrix(); call random_number(A%data)
     ! --> Initialize Krylov subspace.
     allocate(X(kdim + 1)) ; call initialize_krylov_subspace(X)
     call random_number(X(1)%data) ; alpha = X(1)%norm() ; call X(1)%scal(1.0_wp / alpha)
     H = 0.0_wp
     
     ! --> Arnoldi factorization.
     call arnoldi_factorization(A, X, H, info)
     
     ! --> Krylov-Schur condensation.
     call krylov_schur_restart(nblk, X, H, select_eigs)
     
     ! --> Compute Gram matrix associated to the Krylov basis.
     G = 0.0_wp
     call mat_mult(G(1:nblk, 1:nblk),X(1:nblk),X(1:nblk))
     
     ! --> Check result.
     Id = 0.0_wp; forall (i=1:nblk) Id(i, i) = 1.0_wp
     call check(error, norm2(G - Id) < rtol)
     
     return
   contains
     function select_eigs(eigvals) result(selected)
       complex(kind=wp), intent(in) :: eigvals(:)
       logical                      :: selected(size(eigvals))
       selected = (abs(eigvals) > 0.5_wp)
       return
     end function select_eigs
   end subroutine test_krylov_schur_basis_orthogonality
   
   !-----------------------------------------------------------------
   !-----                                                       -----
   !-----     TEST SUITE FOR THE LANCZOS TRIDIAGONALIZATION     -----
   !-----                                                       -----
   !-----------------------------------------------------------------

   subroutine collect_lanczos_tridiag_testsuite(testsuite)
      !> Collection of tests.
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("Lanczos tridiag. full factorization", test_lanczos_tridiag_full_factorization), &
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
      integer, parameter :: kdim = test_size
      !> Tridiagonal matrix.
      real(kind=wp), dimension(kdim + 1, kdim) :: T
      !> Information flag.
      integer :: info
      !> Misc.
      integer :: i, j, k
      real(kind=wp), dimension(test_size, kdim + 1) :: Xdata
      real(kind=wp) :: alpha
      class(rvector), allocatable :: X0(1)

      ! --> Initialize matrix.
      A = spd_matrix(); call random_number(A%data); A%data = matmul(A%data, transpose(A%data))
      ! --> Initialize Krylov subspace.
      allocate (X(1:kdim + 1)); allocate (X0(1))
      call random_number(X0(1)%data);
      call initialize_krylov_subspace(X, X0)
      T = 0.0_wp

      ! --> Lanczos factorization.
      call lanczos_tridiagonalization(A, X, T, info)

      ! --> Check correctness of full factorization.
      do k = 1, kdim + 1
         Xdata(:, k) = X(k)%data
      end do
      ! --> Infinity-norm check.
      alpha = maxval(abs(matmul(A%data, Xdata(:, 1:kdim)) - matmul(Xdata, T)))
      write (*, *) "Infinity-norm      :", alpha
      write (*, *) "Relative tolerance :", rtol
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
      double precision, dimension(kdim + 1, kdim) :: T
      !> Information flag.
      integer :: info
      !> Misc.
      double precision, dimension(kdim + 1, kdim + 1) :: G, Id
      double precision :: alpha
      class(rvector), allocatable :: X0(1)
      integer :: i, j, k

      ! --> Initialize random spd matrix.
      A = spd_matrix(); call random_number(A%data); A%data = matmul(A%data, transpose(A%data))
      ! --> Initialize Krylov subspace.
      allocate (X(1:kdim + 1)); allocate (X0(1))
      call random_number(X0(1)%data);
      call initialize_krylov_subspace(X, X0)
      T = 0.0_wp

      ! --> Lanczos factorization.
      call lanczos_tridiagonalization(A, X, T, info)

      ! --> Compute Gram matrix associated to the Krylov basis.
      G = 0.0_wp
      call mat_mult(G,X,X)

      ! --> Check result.
      Id = 0.0_wp; forall (i=1:kdim) Id(i, i) = 1.0_wp
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

      testsuite = [ &
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
      class(rvector), allocatable :: U(:), V(:)
      !> Krylov subspace dimension.
      integer, parameter :: kdim = test_size
      !> Bidiagonal matrix.
      real(kind=wp) :: B(kdim + 1, kdim)
      !> Information flag.
      integer :: info
      !> Miscellaneous.
      integer :: k
      real(kind=wp) :: alpha
      real(kind=wp) :: Udata(test_size, kdim + 1), Vdata(test_size, kdim + 1)
      class(rvector), allocatable :: X0(1)

      ! --> Initialize matrix.
      A = rmatrix(); call random_number(A%data)
      ! --> Initialize Krylov subspaces.
      allocate (U(1:kdim + 1)); allocate (V(1:kdim + 1)); allocate (X0(1))
      call random_number(X0(1)%data);
      call initialize_krylov_subspace(U, X0)
      call random_number(X0(1)%data); ! new random number for V
      call initialize_krylov_subspace(V, X0)
      B = 0.0_wp

      ! --> Lanczos bidiagonalization.
      call lanczos_bidiagonalization(A, U, V, B, info)

      ! --> Check correctness of full factorization.
      do k = 1, size(U)
         Udata(:, k) = U(k)%data; Vdata(:, k) = V(k)%data
      end do
      ! --> Infinity-norm check.
      alpha = maxval(abs(matmul(A%data, Vdata(:, 1:kdim)) - matmul(Udata, B)))
      write (*, *) "Infinity norm      :", alpha
      write (*, *) "Relative tolerance :", rtol
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
      real(kind=wp) :: T(kdim + 1, kdim + 1)
      !> Information flag.
      integer :: info
      !> Miscellaneous.
      integer :: i, j, k
      real(kind=wp) :: alpha, beta
      real(kind=wp) :: Vdata(test_size, kdim + 1), Wdata(test_size, kdim + 1)
      class(rvector), allocatable :: X0(1)

      ! --> Initialize matrix.
      A = rmatrix(); call random_number(A%data)

      ! --> Initialize Krylov subspaces.
      allocate (V(1:kdim + 1)); allocate (W(1:kdim + 1)); allocate (X0(1))
      call random_number(X0(1)%data);
      call initialize_krylov_subspace(V, X0)
      call random_number(X0(1)%data); ! new random number for W
      call initialize_krylov_subspace(W, X0)
      T = 0.0_wp

      ! --> Nonsymmetric Lanczos factorization.
      call nonsymmetric_lanczos_tridiagonalization(A, V, W, T, info, verbosity=.true.)

      ! --> Check correctness of the factorization.
      do k = 1, size(V)
         Vdata(:, k) = V(k)%data; Wdata(:, k) = W(k)%data
      end do
      ! --> Infinity-norm check.
      alpha = maxval(abs(matmul(A%data, Vdata(:, 1:kdim)) - matmul(Vdata, T(1:kdim + 1, 1:kdim))))
      write (*, *) "Infinity norm      :", alpha
      write (*, *) "Relative tolerance :", rtol
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

      testsuite = [ &
                  new_unittest("Two-sided Arnoldi full factorization", test_two_sided_arnoldi_full_factorization), &
                  new_unittest("Two-sided Arnoldi full factorization adjoint", test_two_sided_arnoldi_full_factorization_bis), &
                  new_unittest("Two-sided Arnoldi basis orthogonality", test_two_sided_arnoldi_basis_orthogonality), &
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
      integer, parameter :: kdim = test_size
      !> Hessenberg matrices.
      real(kind=wp) :: H(kdim + 1, kdim), G(kdim + 1, kdim)
      !> Information flag.
      integer :: info
      !> Miscellaneous.
      integer :: k
      real(kind=wp)  :: Vdata(test_size, kdim + 1), Wdata(test_size, kdim + 1)
      real(kind=wp)  :: alpha
      class(rvector), allocatable :: X0(1)

      !> Initialize matrix.
      A = rmatrix(); call random_number(A%data)
      !> Initialize Krylov subspaces.
      allocate (V(1:kdim + 1)); allocate (W(1:kdim + 1));  allocate (X0(1))
      call random_number(X0(1)%data);
      call initialize_krylov_subspace(V, X0)
      call random_number(X0(1)%data); ! new random number for W
      call initialize_krylov_subspace(W, X0)
      H = 0.0_wp; G = 0.0_wp

      !> Two-sided Arnoldi factoriztion.
      call two_sided_arnoldi_factorization(A, V, W, H, G, info)

      !> Check correctness of the full factorization.
      do k = 1, kdim + 1
         Vdata(:, k) = V(k)%data
      end do
      call check(error, all_close(matmul(A%data, Vdata(:, 1:kdim)), matmul(Vdata, H), rtol, atol))

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
      integer, parameter :: kdim = test_size
      !> Hessenberg matrices.
      real(kind=wp) :: H(kdim + 1, kdim), G(kdim + 1, kdim)
      !> Information flag.
      integer :: info
      !> Miscellaneous.
      integer :: k
      real(kind=wp) :: Vdata(test_size, kdim + 1), Wdata(test_size, kdim + 1)
      real(kind=wp) :: alpha
      class(rvector), allocatable :: X0(1)

      !> Initialize matrix.
      A = rmatrix(); call random_number(A%data)
      !> Initialize Krylov subspaces.
      allocate (V(1:kdim + 1)); allocate (W(1:kdim + 1)); allocate (X0(1))
      call random_number(X0(1)%data);
      call initialize_krylov_subspace(V, X0)
      call random_number(X0(1)%data); ! new random number for W
      call initialize_krylov_subspace(W, X0)
      H = 0.0_wp; G = 0.0_wp

      !> Two-sided Arnoldi factoriztion.
      call two_sided_arnoldi_factorization(A, V, W, H, G, info)

      !> Check correctness of the full factorization.
      do k = 1, kdim + 1
         Wdata(:, k) = W(k)%data
      end do
      call check(error, all_close(matmul(transpose(A%data), Wdata(:, 1:kdim)), matmul(Wdata, G), rtol, atol))

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
      integer, parameter :: kdim = test_size
      !> Hessenberg matrices.
      real(kind=wp) :: H(kdim + 1, kdim), G(kdim + 1, kdim)
      !> Information flag.
      integer :: info
      !> Miscellaneous.
      real(kind=wp) :: M(kdim, kdim), Id(kdim, kdim)
      integer :: i, j, k
      real(kind=wp) :: alpha
      class(rvector), allocatable :: X0(1)

      !> Initialize matrix.
      A = rmatrix(); call random_number(A%data)
      !> Initialize Krylov subspaces.
      allocate (V(1:kdim + 1)); allocate (W(1:kdim + 1)); allocate (X0(1)); 
      call random_number(X0(1)%data);
      call initialize_krylov_subspace(V, X0)
      call random_number(X0(1)%data); ! new random number for W
      call initialize_krylov_subspace(W, X0)
      H = 0.0_wp; G = 0.0_wp

      !> Two-sided Arnoldi factoriztion.
      call two_sided_arnoldi_factorization(A, V, W, H, G, info)

      !> Inner-product matrix.
      M = 0.0_wp
      call mat_mult(M,V(1:kdim),V(1:kdim))

      !> Check results.
      Id = 0.0_wp; forall (i=1:kdim) Id(i, i) = 1.0_wp
      call check(error, norm2(M - Id) < rtol)

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
      integer, parameter :: kdim = test_size
      !> Hessenberg matrices.
      real(kind=wp) :: H(kdim + 1, kdim), G(kdim + 1, kdim)
      !> Information flag.
      integer :: info
      !> Miscellaneous.
      real(kind=wp) :: M(kdim, kdim), Id(kdim, kdim)
      integer :: i, j, k
      real(kind=wp) :: alpha
      class(rvector), allocatable :: X0(1)

      !> Initialize matrix.
      A = rmatrix(); call random_number(A%data)
      !> Initialize Krylov subspaces.
      allocate (V(1:kdim + 1)); allocate (W(1:kdim + 1)); allocate (X0(1)); 
      call random_number(X0(1)%data);
      call initialize_krylov_subspace(V, X0)
      call random_number(X0(1)%data); ! new random number for W
      call initialize_krylov_subspace(W, X0)
      H = 0.0_wp; G = 0.0_wp

      !> Two-sided Arnoldi factoriztion.
      call two_sided_arnoldi_factorization(A, V, W, H, G, info)

      !> Inner-product matrix.
      M = 0.0_wp
      call mat_mult(M,V(1:kdim),V(1:kdim))
      
      !> Check results.
      Id = 0.0_wp; forall (i=1:kdim) Id(i, i) = 1.0_wp
      call check(error, norm2(M - Id) < rtol)

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

      testsuite = [ &
                  new_unittest("QR factorization", test_qr_factorization), &
                  new_unittest("Q orthonormality", test_qr_basis_orthonormality), &
                  new_unittest("QR worst case breakdown", test_qr_breakdown) &
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
      !> Krylov subspace dimension.
      integer, parameter :: kdim = 3
      !> GS factors.
      real(kind=wp) :: R(kdim, kdim)
      !> Information flag.
      integer :: info
      !> Misc.
      integer :: i, j, k
      real(kind=wp) :: Amat(test_size, kdim), Qmat(test_size, kdim)
      real(kind=wp) :: alpha

      ! --> Initialize matrix.
      allocate (A(1:kdim)); 
      do k = 1, size(A)
         call random_number(A(k)%data)
      end do
      ! --> Copy input matrix data for comparison
      do k = 1, kdim
         Amat(:, k) = A(k)%data
      end do
      R = 0.0_wp
      ! --> In-place QR factorization.
      call qr_factorization(A, R, info)
      ! --> Extract data
      do k = 1, kdim
         Qmat(:, k) = A(k)%data
      end do
      ! --> Check correctness of QR factorization.
      call check(error, all_close(Amat, matmul(Qmat, R), rtol, atol))

      return
   end subroutine test_qr_factorization

   subroutine test_qr_basis_orthonormality(error)
      ! This function checks the correctness of the QR implementation by
      ! verifying that the obtained basis is orthonormal, i.e. Q.T @ Q = I.
      ! A random matrix is used for testing.

      !> Error type to be returned.
      type(error_type), allocatable, intent(out) :: error
      !> Test matrix.
      class(rvector), dimension(:), allocatable  :: A
      !> Krylov subspace dimension.
      integer, parameter :: kdim = 3
      !> GS factors.
      real(kind=wp) :: R(kdim, kdim)
      real(kind=wp) :: Id(kdim, kdim)
      !> Information flag.
      integer :: info
      !> Misc.
      integer :: i, k
      real(kind=wp) :: Qmat(test_size, kdim)
      real(kind=wp) :: alpha

      ! --> Initialize matrix.
      allocate (A(1:kdim)); 
      do k = 1, size(A)
         call random_number(A(k)%data)
      enddo
      R = 0.0_wp
      ! --> In-place QR factorization.
      call qr_factorization(A, R, info)
      ! --> Extract data
      do k = 1, kdim
         Qmat(:, k) = A(k)%data
      enddo
      ! --> Identity
      Id = 0.0_wp; forall (i=1:kdim) Id(i, i) = 1.0_wp
      ! --> Check correctness of QR factorization.
      call check(error, all_close(Id, matmul(transpose(Qmat), Qmat), rtol, atol))

   end subroutine test_qr_basis_orthonormality

   subroutine test_qr_breakdown(error)
      ! This function checks the correctness of the QR implementation in a worst
      ! case scenario where the basis vectors are nearly linearly dependent.
      ! The snaller the value of eps, the closer the columns are to linear dependence.

      !> Error type to be returned.
      type(error_type), allocatable, intent(out) :: error
      !> Test matrix.
      class(rvector), dimension(:), allocatable  :: A
      !> Krylov subspace dimension.
      integer,       parameter :: kdim = 6
      real(kind=wp), parameter :: eps = 1e-10
      !> GS factors.
      real(kind=wp) :: R(kdim, kdim)
      real(kind=wp) :: Id(kdim, kdim)
      !> Information flag.
      integer :: info
      !> Misc.
      integer :: i, k
      real(kind=wp) :: Qmat(test_size, kdim)
      class(rvector) , dimension(:), allocatable :: wrk

      ! --> Initialize matrix with worst case scenario
      allocate (A(1:kdim)); call mat_zero(A)
      allocate (wrk(1))
      call random_number(A(1)%data)
      do k = 2, size(A)
         ! each column is different from the previous one by eps
         call random_number(wrk(1)%data)
         call A(k)%axpby(0.0_wp, A(k-1), 1.0_wp)
         call A(k)%axpby(1.0_wp, wrk(1), eps)
      end do
      R = 0.0_wp
      ! --> In-place QR factorization.
      call qr_factorization(A, R, info)
      ! --> Extract data
      do k = 1, kdim
         Qmat(:, k) = A(k)%data
      end do
      ! --> Identity
      Id = 0.0_wp; forall (i=1:kdim) Id(i, i) = 1.0_wp
      ! --> Check correctness of QR factorization.
      call check(error, all_close(Id, matmul(transpose(Qmat), Qmat), rtol, atol))

   end subroutine test_qr_breakdown

end module TestKrylov
