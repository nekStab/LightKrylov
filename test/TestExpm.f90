module TestExpm
   use LightKrylov
   use TestVector
   use TestMatrices
   use TestUtils
   use LightKrylov_expmlib
   use LightKrylov_utils
   use testdrive  , only : new_unittest, unittest_type, error_type, check
   use stdlib_linalg, only: diag
   use stdlib_math, only : all_close
   use stdlib_io_npy, only : save_npy
   implicit none
 
   private
 
   public :: collect_dense_matrix_functions_testsuite, collect_kexpm_testsuite

contains

   !---------------------------------------------------------
   !-----                                               -----
   !-----     TEST SUITE FOR DENSE MATRIX FUNCTIONS     -----
   !-----                                               -----
   !---------------------------------------------------------
 
  subroutine collect_dense_matrix_functions_testsuite(testsuite)
     !> Collection of tests.
     type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [&
           new_unittest("Dense Matrix Exponential", test_dense_matrix_exponential), &
           new_unittest("Dense Matrix Square Root - SPD", test_dense_sqrtm_SPD), &
           new_unittest("Dense Matrix Square Root - symmetric positive semi-definite", test_dense_sqrtm_pos_semidefinite) &
           ]

      return
   end subroutine collect_dense_matrix_functions_testsuite

   subroutine test_dense_matrix_exponential(error)
      !> This function tests the scaling and squaring followed by rational Pade approximation
      ! of the matrix exponential for a matrix for which the exponential propagator is known
      ! analytically

      !> Error type to be returned.
      type(error_type), allocatable, intent(out) :: error
      !> Problem dimension.
      integer, parameter :: n = 5
      integer, parameter :: m = 6
      !> Test matrix.
      real(kind=wp) :: A(n, n)
      real(kind=wp) :: E(n, n)
      real(kind=wp) :: Eref(n, n)
      integer :: i, j

      ! --> Initialize matrix.
      A = 0.0_wp
      do i = 1, n-1
         A(i,i+1) = m*1.0_wp
      end do
      ! --> Reference with analytical exponential
      Eref = 0.0_wp
      forall (i=1:n) Eref(i, i) = 1.0_wp
      do i = 1, n-1
         do j = 1, n-i
            Eref(i,i+j) = Eref(i,i+j-1)*m/j
         end do
      end do
      ! --> Compute exponential numerically
      E = 0.0_wp
      call expm(E, A)

      call check(error, maxval(E-Eref) < rtol)

      return
   end subroutine test_dense_matrix_exponential

   subroutine test_dense_sqrtm_SPD(error)
      !> This function tests the matrix version of the sqrt function for the case of
      ! a SPD matrix

      !> Error type to be returned.
      type(error_type), allocatable, intent(out) :: error
      !> Problem dimension.
      integer, parameter :: n = 5
      !> Test matrix.
      real(kind=wp) :: A(n, n)
      real(kind=wp) :: sqrtmA(n, n)
      complex(kind=wp) :: lambda(n)
      integer :: i, j

      ! --> Initialize matrix.
      call random_number(A)
      ! make SPD
      A = 0.5_wp*(A + transpose(A))
      call eig(A, sqrtmA, lambda)
      do i = 1,n
         lambda(i) = abs(lambda(i)) + 0.1_wp
      end do
      ! reconstruct matrix
      A = matmul(sqrtmA, matmul(diag(lambda), transpose(sqrtmA)))
      
      ! compute matrix square root
      call sqrtm(sqrtmA, A)

      call check(error, maxval(matmul(sqrtmA, sqrtmA) - A) < 10*atol)

      return
   end subroutine test_dense_sqrtm_SPD

   subroutine test_dense_sqrtm_pos_semidefinite(error)
      !> This function tests the matrix version of the sqrt function for the case 
      ! of a symmetric semi-definite matrix

      !> Error type to be returned.
      type(error_type), allocatable, intent(out) :: error
      !> Problem dimension.
      integer, parameter :: n = 5
      !> Test matrix.
      real(kind=wp) :: A(n, n)
      real(kind=wp) :: sqrtmA(n, n)
      complex(kind=wp) :: lambda(n)
      integer :: i, j

      ! --> Initialize matrix.
      call random_number(A)
      ! make positive semi-definite
      A = 0.5_wp*(A + transpose(A))
      call eig(A, sqrtmA, lambda)
      do i = 1,n-1
         lambda(i) = abs(lambda(i)) + 0.1_wp
      end do
      lambda(n) = 0.0_wp
      ! reconstruct matrix
      A = matmul(sqrtmA, matmul(diag(lambda), transpose(sqrtmA)))

      ! compute matrix square root
      call sqrtm(sqrtmA, A)

      call check(error, abs(maxval(matmul(sqrtmA, sqrtmA) - A)) < 10*atol)

      return
   end subroutine test_dense_sqrtm_pos_semidefinite
   
   !---------------------------------------------------------
   !-----                                               -----
   !-----     TEST SUITE FOR THE MATRIX EXPONENTIAL     -----
   !-----                                               -----
   !---------------------------------------------------------

   subroutine collect_kexpm_testsuite(testsuite)
      !> Collection of tests.
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [&
            new_unittest("Krylov Matrix Exponential", test_krylov_matrix_exponential), &
            new_unittest("Block Krylov Matrix Exponential", test_block_krylov_matrix_exponential) &
            ]

       return
   end subroutine collect_kexpm_testsuite

   subroutine test_krylov_matrix_exponential(error)
      !> This function tests the Krylov based approximation of the action of the exponential
      ! propagator against the dense computation for a random operator, a random RHS and a 
      ! typical value of tau.
      !> Error type to be returned.
      type(error_type), allocatable, intent(out) :: error
      class(rmatrix), allocatable :: A
      !> Basis vectors.
      class(rvector), allocatable :: Q
      class(rvector), allocatable :: Xref
      class(rvector), allocatable :: Xkryl
      !> Krylov subspace dimension.
      integer, parameter :: kdim = test_size
      !> Test matrix.
      real(kind=wp) :: Adata(kdim, kdim)
      real(kind=wp) :: Edata(kdim, kdim)
      !> GS factors.
      real(kind=wp) :: R(kdim, kdim)
      real(kind=wp) :: Id(kdim, kdim)
      !> Information flag.
      integer :: info
      !> Test parameters
      integer, parameter         :: nkmax = 15
      real(kind=wp), parameter   :: tau   = 0.1_wp
      real(kind=wp), parameter   :: tol   = 1e-10_wp
      logical, parameter         :: verb  = .true.
      !> Misc.
      integer :: i,j,k
      real(kind=wp) :: Xdata(test_size), Qdata(test_size)
      real(kind=wp) :: err

!#define DEBUG

      Adata = 0.0_wp; Edata = 0.0_wp; Xdata = 0.0_wp
      allocate(Q); allocate(Xref); allocate(Xkryl)
      call Xref%zero()
      call Xkryl%zero()

      ! --> Initialize operator.
      A = rmatrix()
      call init_rand(A)
      call get_data(Adata, A)   
      ! --> Initialize rhs.
      call init_rand(Q)
      call get_data(Qdata, Q)

      !> Comparison is dense computation (10th order Pade approximation)
      call expm(Edata, tau*Adata)
      Xdata = matmul(Edata,Qdata)

      !> Copy reference data into Krylov vector
      call put_data(Xref, Xdata)

      !> Compute Krylov matrix exponential using the arnoldi method
      call kexpm(Xkryl, A, Q, tau, tol, info, verbosity = verb, kdim = nkmax)

#ifdef DEBUG
      !> Save test data to disk.
      call save_npy("debug/test_krylov_expm_operator.npy", Adata)
      call save_npy("debug/test_krylov_expm_rhs.npy", Qdata)
      call save_npy("debug/test_krylov_expm_ref.npy", Xdata)
      call get_data(Xdata, Xkryl)
      call save_npy("debug/test_krylov_expm_kexpm.npy", Xdata)
#endif
#undef DEBUG

      call Xkryl%axpby(1.0_wp, Xref, -1.0_wp)

      !> Compute 2-norm of the error
      err = Xkryl%norm()
      if (verb) write(*, *) '    true error:          ||error||_2 = ', err

      call check(error, err < rtol)

      return
   end subroutine test_krylov_matrix_exponential

   subroutine test_block_krylov_matrix_exponential(error)
      !> This function tests the Krylov based approximation of the action of the exponential
      ! propagator against the dense computation for a random operator, a random RHS and a 
      ! typical value of tau.

      !> Error type to be returned.
      type(error_type), allocatable, intent(out) :: error
      class(rmatrix), allocatable :: A
      !> Basis vectors.
      class(rvector), allocatable :: Q(:)
      class(rvector), allocatable :: Xref(:)
      class(rvector), allocatable :: Xkryl(:)
      class(rvector), allocatable :: Xkryl_block(:)
      !> Krylov subspace dimension.
      integer, parameter :: kdim = test_size
      !> Test matrix.
      real(kind=wp) :: Adata(kdim, kdim)
      real(kind=wp) :: Edata(kdim, kdim)
      !> GS factors.
      real(kind=wp) :: R(kdim, kdim)
      real(kind=wp) :: Id(kdim, kdim)
      !> Information flag.
      integer :: info
      !> Test parameters
      integer, parameter         :: nkmax = 15
      integer, parameter         :: p     = 3
      real(kind=wp), parameter   :: tau   = 0.1_wp
      real(kind=wp), parameter   :: tol   = 1e-10_wp
      logical, parameter         :: verb  = .true.
      !> Misc.
      integer :: i,j,k
      real(kind=wp) :: Xdata(test_size,p), Qdata(test_size,p)
      real(kind=wp) :: alpha
      real(kind=wp) :: err(p,p)

!#define DEBUG

      Adata = 0.0_wp; Edata = 0.0_wp; Xdata = 0.0_wp
      allocate(Xref(1:p)) ; call mat_zero(Xref)
      allocate(Xkryl(1:p)) ; call mat_zero(Xkryl)
      allocate(Xkryl_block(1:p)) ; call mat_zero(Xkryl_block)

      ! --> Initialize operator.
      A = rmatrix()
      call init_rand(A)
      call get_data(Adata, A)
      ! --> Initialize rhs.
      allocate(Q(1:p))
      call init_rand(Q) 
      call get_data(Qdata, Q)

      !> Comparison is dense computation (10th order Pade approximation)
      call expm(Edata, tau*Adata)
      Xdata = matmul(Edata,Qdata)
      !> Copy reference data into Krylov vector
      call put_data(Xref, Xdata)

      !> Compute Krylov matrix exponential using sequential arnoldi method for each input column
      if (verb) write(*,*) 'SEQUENTIAL ARNOLDI'
      do i = 1,p
         if (verb) write(*,*) '    column',i
         call kexpm(Xkryl(i:i), A, Q(i:i), tau, tol, info, verbosity = verb, kdim = nkmax)
      end do
      if (verb) write(*,*) 'BLOCK-ARNOLDI'
      !> Compute Krylov matrix exponential using block-arnoldi method
      call kexpm(Xkryl_block(1:p), A, Q(1:p), tau, tol, info, verbosity = verb, kdim = nkmax)

      do i = 1,p
         call Xkryl(i)%axpby(1.0_wp, Xref(i), -1.0_wp)
         call Xkryl_block(i)%axpby(1.0_wp, Xref(i), -1.0_wp)
      end do

#ifdef DEBUG
      !> Save test data to disk.
      call save_npy("debug/test_block_krylov_expm_operator.npy", Adata)
      call save_npy("debug/test_block_krylov_expm_rhs.npy", Qdata)
      call save_npy("debug/test_block_krylov_expm_ref.npy", Xdata)
      call get_data(Xdata, Xkryl)
      call save_npy("debug/test_block_krylov_expm_kexpm_seq.npy", Xdata)
      call get_data(Xdata, Xkryl_block)
      call save_npy("debug/test_block_krylov_expm_kexpm_blk.npy", Xdata)
#endif
#undef DEBUG

      !> Compute 2-norm of the error
      if (verb) then
         call mat_mult(err,Xkryl(1:p),Xkryl(1:p))
         alpha = sqrt(norm2(err))
         write(*,*) '--------------------------------------------------------------------'
         write(*, *) '    true error (seq.):   ||error||_2 = ', alpha
      endif
      call mat_mult(err,Xkryl_block(1:p),Xkryl_block(1:p))
      alpha = sqrt(norm2(err))
      if (verb) write(*, *) '    true error (block):  ||error||_2 = ', alpha

      call check(error, alpha < rtol)

      return
   end subroutine test_block_krylov_matrix_exponential

end module TestExpm
