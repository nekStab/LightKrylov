module TestKrylov
   use LightKrylov
   use TestVector
   use TestMatrices
   use TestUtils
   use lightkrylov_Utils
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use stdlib_math, only: all_close, is_close
   use stdlib_linalg, only: eye
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
      real(kind=wp) :: Xdata(test_size, kdim + 1)
      class(rvector), allocatable :: X0(1)

      ! --> Initialize matrix.
      A = rmatrix(); call init_rand(A)
      ! --> Initialize Krylov subspace.
      allocate (X(1:kdim + 1)); allocate (X0(1))
      call init_rand(X0)
      call initialize_krylov_subspace(X, X0)
      H = 0.0_wp

      ! --> Arnoldi factorization.
      call arnoldi_factorization(A, X, H, info)

      ! --> Check correctness of full factorization.
      call get_data(Xdata, X)

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

      ! --> Initialize random matrix.
      A = rmatrix(); call init_rand(A)
      ! --> Initialize Krylov subspace.
      allocate (X(1:kdim + 1)); allocate (X0(1))
      call init_rand(X0)
      call initialize_krylov_subspace(X, X0)
      H = 0.0_wp

      ! --> Arnoldi factorization.
      call arnoldi_factorization(A, X, H, info)
      
      ! --> Compute Gram matrix associated to the Krylov basis.
      G = 0.0_wp
      call mat_mult(G,X(1:kdim),X(1:kdim))
      
      ! --> Check result.
      Id = eye(kdim)
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
      real(kind=wp) :: Xdata(test_size, p*(kdim + 1))
      double precision, dimension(p*kdim, p*kdim) :: G, Id
      class(rvector), dimension(:), allocatable :: X0
      
      ! --> Initialize matrix.
      A = rmatrix(); call init_rand(A)
      ! --> Initialize Krylov subspace.
      allocate (X(1:p*(kdim + 1))); allocate (X0(1:p)); 
      call init_rand(X0)
      call initialize_krylov_subspace(X, X0)
      H = 0.0_wp

      ! --> Arnoldi factorization.
      call arnoldi_factorization(A, X, H, info, block_size = p) 
      
      ! --> Check correctness of full factorization.
      call get_data(Xdata, X)

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

      !! --> Initialize random matrix.
      A = rmatrix(); call init_rand(A)
      !! --> Initialize Krylov subspace.
      allocate (X(1:p*(kdim + 1))); allocate (X0(1:p)); 
      call init_rand(X0)
      call initialize_krylov_subspace(X, X0)
      H = 0.0_wp

      ! --> Arnoldi factorization.
      call arnoldi_factorization(A, X, H, info, block_size = p)

      ! --> Compute Gram matrix associated to the Krylov basis.
      G = 0.0_wp
      call mat_mult(G,X(1:p*kdim),X(1:p*kdim))

      ! --> Check result.
      Id = eye(p*kdim)
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
     class(rvector), dimension(:), allocatable :: X0
     !> Krylov subspace dimension.
     integer, parameter :: kdim = 10
     !> Hessenberg matrix.
     real(kind=wp) :: H(kdim + 1, kdim)
     !> Information flag.
     integer :: info
     !> Misc.
     integer :: nblk
     real(kind=wp) :: Xdata(test_size, kdim + 1)
     real(kind=wp) :: alpha
     
     ! --> Initialize matrix.
     A = rmatrix(); call init_rand(A)
     ! --> Initialize Krylov subspace.
     allocate (X(1:kdim + 1)); allocate (X0(1))
     call init_rand(X0)
     call initialize_krylov_subspace(X, X0)
     H = 0.0_wp
     
     ! --> Arnoldi factorization.
     call arnoldi_factorization(A, X, H, info)
     
     ! --> Krylov-Schur condensation.
     call krylov_schur_restart(nblk, X, H, select_eigs)
     
     ! --> Check correctness of full factorization.
     call get_data(Xdata, X)
     
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
     A = rmatrix(); call init_rand(A) ; m = sum(A%data, 2)/test_size
     do k = 1, test_size
        A%data(:, k) = A%data(:, k) - m
     enddo
     !A%data = 0.5 * (A%data + transpose(A%data))
     ! --> Initialize Krylov subspace.
     allocate (X(p*(kdim + 1))) ; allocate(X0(p))
     call init_rand(X0)
     call initialize_krylov_subspace(X, X0)
     H = 0.0_wp
     
     ! --> Arnoldi factorization.
     call arnoldi_factorization(A, X, H, info, block_size=p)
     
     ! --> Krylov-Schur condensation.
     call krylov_schur_restart(nblk, X, H, select_eigs, p)
     
     ! --> Check correctness of full factorization.
     call get_data(Xdata, X)
     
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
     integer :: nblk
     
     ! --> Initialize random matrix.
     A = rmatrix(); call init_rand(A)
     ! --> Initialize Krylov subspace.
     allocate (X(1:kdim + 1)); allocate (X0(1))
     call init_rand(X0)
     call initialize_krylov_subspace(X, X0)
     H = 0.0_wp
     
     ! --> Arnoldi factorization.
     call arnoldi_factorization(A, X, H, info)
     
     ! --> Krylov-Schur condensation.
     call krylov_schur_restart(nblk, X, H, select_eigs)
     
     ! --> Compute Gram matrix associated to the Krylov basis.
     G = 0.0_wp
     call mat_mult(G(1:nblk, 1:nblk),X(1:nblk),X(1:nblk))
     
     ! --> Check result.
     Id = eye(nblk)
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
      real(kind=wp), dimension(test_size, kdim + 1) :: Xdata
      real(kind=wp) :: alpha
      class(rvector), allocatable :: X0(1)

      ! --> Initialize matrix.
      A = spd_matrix(); call init_rand(A)
      ! --> Initialize Krylov subspace.
      allocate (X(1:kdim + 1)); allocate (X0(1))
      call init_rand(X0)
      call initialize_krylov_subspace(X, X0)
      T = 0.0_wp

      ! --> Lanczos factorization.
      call lanczos_tridiagonalization(A, X, T, info)

      ! --> Check correctness of full factorization.
      call get_data(Xdata, X)

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
      class(rvector), allocatable :: X0(1)

      ! --> Initialize random spd matrix.
      A = spd_matrix(); call init_rand(A)
      ! --> Initialize Krylov subspace.
      allocate (X(1:kdim + 1)); allocate (X0(1))
      call init_rand(X0)
      call initialize_krylov_subspace(X, X0)
      T = 0.0_wp

      ! --> Lanczos factorization.
      call lanczos_tridiagonalization(A, X, T, info)

      ! --> Compute Gram matrix associated to the Krylov basis.
      G = 0.0_wp
      call mat_mult(G,X,X)

      ! --> Check result.
      Id = 0.0_wp; Id(1:kdim,1:kdim) = eye(kdim)

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
      real(kind=wp) :: alpha
      real(kind=wp) :: Udata(test_size, kdim + 1), Vdata(test_size, kdim + 1)
      class(rvector), allocatable :: X0(1)

      ! --> Initialize matrix.
      A = rmatrix(); call init_rand(A)
      ! --> Initialize Krylov subspaces.
      allocate (U(1:kdim + 1)); allocate (V(1:kdim + 1)); allocate (X0(1))
      call init_rand(X0)
      call initialize_krylov_subspace(U, X0)
      call init_rand(X0)
      call initialize_krylov_subspace(V, X0)
      B = 0.0_wp

      ! --> Lanczos bidiagonalization.
      call lanczos_bidiagonalization(A, U, V, B, info)

      ! --> Check correctness of full factorization.
      call get_data(Udata, U)
      call get_data(Vdata, V)

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
      real(kind=wp) :: alpha
      real(kind=wp) :: Vdata(test_size, kdim + 1), Wdata(test_size, kdim + 1)
      class(rvector), allocatable :: X0(1)

      ! --> Initialize matrix.
      A = rmatrix(); call init_rand(A)
      ! --> Initialize Krylov subspaces.
      allocate (V(1:kdim + 1)); allocate (W(1:kdim + 1)); allocate (X0(1))
      call init_rand(X0)
      call initialize_krylov_subspace(V, X0)
      call init_rand(X0)
      call initialize_krylov_subspace(W, X0)
      T = 0.0_wp

      ! --> Nonsymmetric Lanczos factorization.
      call nonsymmetric_lanczos_tridiagonalization(A, V, W, T, info, verbosity=.true.)

      ! --> Check correctness of the factorization.
      call get_data(Vdata, V)
      call get_data(Wdata, W)
      
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
      real(kind=wp)  :: Vdata(test_size, kdim + 1), Wdata(test_size, kdim + 1)
      class(rvector), allocatable :: X0(1)

      !> Initialize matrix.
      A = rmatrix(); call init_rand(A)
      !> Initialize Krylov subspaces.
      allocate (V(1:kdim + 1)); allocate (W(1:kdim + 1));  allocate (X0(1))
      call init_rand(X0)
      call initialize_krylov_subspace(V, X0)
      call init_rand(X0)
      call initialize_krylov_subspace(W, X0)
      H = 0.0_wp; G = 0.0_wp

      !> Two-sided Arnoldi factoriztion.
      call two_sided_arnoldi_factorization(A, V, W, H, G, info)

      !> Check correctness of the full factorization.
      call get_data(Vdata, V)

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
      real(kind=wp) :: Vdata(test_size, kdim + 1), Wdata(test_size, kdim + 1)
      class(rvector), allocatable :: X0(1)

      !> Initialize matrix.
      A = rmatrix(); call init_rand(A)
      !> Initialize Krylov subspaces.
      allocate (V(1:kdim + 1)); allocate (W(1:kdim + 1)); allocate (X0(1))
      call init_rand(X0)
      call initialize_krylov_subspace(V, X0)
      call init_rand(X0) ! new random number for W
      call initialize_krylov_subspace(W, X0)
      H = 0.0_wp; G = 0.0_wp

      !> Two-sided Arnoldi factoriztion.
      call two_sided_arnoldi_factorization(A, V, W, H, G, info)

      !> Check correctness of the full factorization.
      call get_data(Wdata, W)

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
      class(rvector), allocatable :: X0(1)

      !> Initialize matrix.
      A = rmatrix(); call init_rand(A)
      !> Initialize Krylov subspaces.
      allocate (V(1:kdim + 1)); allocate (W(1:kdim + 1)); allocate (X0(1)); 
      call init_rand(X0)
      call initialize_krylov_subspace(V, X0)
      call init_rand(X0) ! new random number for W
      call initialize_krylov_subspace(W, X0)
      H = 0.0_wp; G = 0.0_wp

      !> Two-sided Arnoldi factoriztion.
      call two_sided_arnoldi_factorization(A, V, W, H, G, info)

      !> Inner-product matrix.
      M = 0.0_wp
      call mat_mult(M,V(1:kdim),V(1:kdim))

      !> Check results.
      Id = eye(kdim)
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
      class(rvector), allocatable :: X0(1)

      !> Initialize matrix.
      A = rmatrix(); call init_rand(A)
      !> Initialize Krylov subspaces.
      allocate (V(1:kdim + 1)); allocate (W(1:kdim + 1)); allocate (X0(1)); 
      call init_rand(X0)
      call initialize_krylov_subspace(V, X0)
      call init_rand(X0)
      call initialize_krylov_subspace(W, X0)
      H = 0.0_wp; G = 0.0_wp

      !> Two-sided Arnoldi factoriztion.
      call two_sided_arnoldi_factorization(A, V, W, H, G, info)

      !> Inner-product matrix.
      M = 0.0_wp
      call mat_mult(M,V(1:kdim),V(1:kdim))
      
      !> Check results.
      Id = eye(kdim)
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
                  new_unittest("QR worst case breakdown", test_qr_breakdown), &
                  new_unittest("Pivoted QR for exactly rank deficient matrices", test_piv_qr_absolute_rank_deficiency), &
                  new_unittest("Pivoted QR for numerically rank deficient matrices", test_piv_qr_num_rank_deficiency) &
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
      integer, parameter :: kdim = test_size
      !> GS factors.
      real(kind=wp) :: R(kdim, kdim)
      !> Permutation matrix.
      real(kind=wp) :: P(kdim, kdim)
      !> Information flag.
      integer :: info
      !> Misc.
      real(kind=wp) :: Adata(test_size, kdim), Qdata(test_size, kdim)

      ! --> Initialize matrix.
      allocate (A(1:kdim)); call init_rand(A)

      ! --> Copy input matrix data for comparison
      call get_data(Adata, A)
      
      ! --> In-place QR factorization.
      call qr_factorization(A, R, P, info)

      ! --> Extract data
      call get_data(Qdata, A)

      ! --> Check correctness of QR factorization.
      call check(error, all_close(Adata, matmul(Qdata, R), rtol, atol))

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
      integer, parameter :: kdim = test_size
      !> GS factors.
      real(kind=wp) :: R(kdim, kdim)
      !> Permutation matrix.
      real(kind=wp) :: P(kdim, kdim)
      real(kind=wp) :: Id(kdim, kdim)
      !> Information flag.
      integer :: info
      !> Misc.
      real(kind=wp) :: Qdata(test_size, kdim)

      ! --> Initialize matrix.
      allocate (A(1:kdim)); call init_rand(A)
      
      ! --> In-place QR factorization.
      call qr_factorization(A, R, P, info)

      ! --> Extract data
      call get_data(Qdata, A)

      ! --> Check correctness of QR factorization.
      Id = eye(kdim)
      call check(error, all_close(Id, matmul(transpose(Qdata), Qdata), rtol, atol))

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
      !> Permutation matrix.
      real(kind=wp) :: P(kdim, kdim)
      real(kind=wp) :: Id(kdim, kdim)
      !> Information flag.
      integer :: info
      !> Misc.
      integer :: k
      real(kind=wp) :: Qdata(test_size, kdim)
      class(rvector), allocatable :: wrk

      ! --> Initialize matrix with worst case scenario
      allocate (A(1:kdim)); call init_rand(A)
      allocate (wrk)
      do k = 2, size(A)
         ! each column is different from the previous one by eps
         call init_rand(wrk)
         call A(k)%axpby(0.0_wp, A(k-1), 1.0_wp)
         call A(k)%axpby(1.0_wp, wrk, eps)
      end do
      
      ! --> In-place QR factorization.
      call qr_factorization(A, R, P, info)

      ! --> Extract data
      call get_data(Qdata, A)

      ! --> Identity
      Id = eye(kdim)
      ! --> Check correctness of QR factorization.
      call check(error, all_close(Id, matmul(transpose(Qdata), Qdata), rtol, atol))

   end subroutine test_qr_breakdown

   subroutine test_piv_qr_absolute_rank_deficiency(error)
      ! This function checks the correctness of the pivoted  QR implementation 
      ! by testing it on a rank deficient matrix

      !> Error type to be returned.
      type(error_type), allocatable, intent(out) :: error
      !> Test matrix.
      class(rvector), dimension(:), allocatable  :: A
      !> Krylov subspace dimension.
      integer, parameter :: kdim = 20
      !> Number of zero columns
      integer, parameter :: nzero = 5
      !> GS factors.
      real(kind=wp) :: R(kdim, kdim)
      !> Permutation matrix.
      real(kind=wp) :: P(kdim, kdim)
      real(kind=wp) :: Id(kdim, kdim)
      !> Information flag.
      integer :: info
      !> Misc.
      integer :: k, idx, rk
      real(kind=wp) :: Adata(test_size, kdim), Qdata(test_size, kdim)
      real(kind=wp) :: alpha
      logical       :: mask(kdim)

      ! Effective rank 
      rk = kdim - nzero

      ! --> Initialize matrix.
      allocate (A(1:kdim)); call init_rand(A)

      ! add zero vectors at random places
      mask = .true.; k = nzero
      do while ( k .gt. 0 )
         call random_number(alpha)
         idx = 1 + floor(kdim*alpha)
         if (mask(idx)) then
            A(idx)%data = 0.0_wp
            mask(idx) = .false.
            k = k - 1
         end if
      end do

      ! copy data
      call get_data(Adata, A)

      ! --> In-place QR factorization.
      call qr_factorization(A, R, P, info,  ifpivot = .true.)

      ! --> Extract data
      call get_data(Qdata, A)

      ! --> Check correctness of QR factorization.
      call check(error, all_close(matmul(Adata,P), matmul(Qdata, R), rtol, atol))

   end subroutine test_piv_qr_absolute_rank_deficiency

   subroutine test_piv_qr_num_rank_deficiency(error)
      ! This function checks the correctness of the pivoted  QR implementation 
      ! by testing it on a rank deficient matrix

      !> Error type to be returned.
      type(error_type), allocatable, intent(out) :: error
      !> Test matrix.
      class(rvector), dimension(:), allocatable  :: A
      !> Krylov subspace dimension.
      integer, parameter :: kdim = 20
      !> Number of zero columns
      integer, parameter :: nzero = 5
      !> GS factors.
      real(kind=wp) :: R(kdim, kdim)
      !> Permutation matrix.
      real(kind=wp) :: P(kdim, kdim)
      real(kind=wp) :: Id(kdim, kdim)
      !> Information flag.
      integer :: info
      !> Misc.
      integer :: k, idx, rk
      real(kind=wp) :: Adata(test_size, kdim), Qdata(test_size, kdim)
      real(kind=wp) :: rnd(test_size,2)
      real(kind=wp) :: alpha
      logical       :: mask(kdim)

      ! Effective rank 
      rk = kdim - nzero

      ! --> Initialize matrix.
      allocate (A(1:kdim)); call init_rand(A)

      ! add zero vectors at random places
      mask = .true.; k = 1
      do while ( k .le. nzero )
         call random_number(alpha)
         idx = 1 + floor(kdim*alpha)
         if (mask(idx)) then
            call random_number(rnd)
            A(idx)%data = rnd(:,1)*A(k)%data + atol*rnd(:,2)
            mask(idx) = .false.
            k = k + 1
         end if
      end do

      ! copy data
      call get_data(Adata, A)

      ! --> In-place QR factorization.
      call qr_factorization(A, R, P, info, ifpivot = .true.)
     
      ! --> Extract data
      call get_data(Qdata, A)

      ! --> Check correctness of QR factorization.
      call check(error, all_close(matmul(Adata,P), matmul(Qdata, R), rtol, atol))

   end subroutine test_piv_qr_num_rank_deficiency

end module TestKrylov
