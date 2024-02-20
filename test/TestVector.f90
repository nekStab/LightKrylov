module TestVector
   use LightKrylov
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use stdlib_math, only: is_close, all_close

   implicit none

   private

   public :: collect_real_vector_testsuite, test_size

   integer, parameter :: test_size = 100

   type, extends(abstract_vector), public :: rvector
      real(kind=wp), dimension(test_size) :: data = 0.0_wp
   contains
      private
      procedure, pass(self), public :: zero
      procedure, pass(self), public :: dot
      procedure, pass(self), public :: scal
      procedure, pass(self), public :: axpby
   end type rvector

contains

   !-----------------------------------------------------------
   !-----                                                 -----
   !-----     DEFINITION OF THE TYPE-BOUND PROCEDURES     -----
   !-----                                                 -----
   !-----------------------------------------------------------

   !--> Zero-out a vector.
   subroutine zero(self)
      class(rvector), intent(inout) :: self
      self%data = 0.0_wp
      return
   end subroutine zero

   double precision function dot(self, vec) result(alpha)
      class(rvector), intent(in)         :: self
      class(abstract_vector), intent(in) :: vec

      select type (vec)
      type is (rvector)
         alpha = dot_product(self%data, vec%data)
      end select
      return
   end function dot

   ! --> In-place scalar multiplication.
   subroutine scal(self, alpha)
      class(rvector), intent(inout) :: self
      real(kind=wp), intent(in) :: alpha
      self%data = self%data*alpha
      return
   end subroutine scal

   ! --> axpby interface
   subroutine axpby(self, alpha, vec, beta)
      class(rvector), intent(inout) :: self
      class(abstract_vector), intent(in) :: vec
      real(kind=wp), intent(in) :: alpha, beta

      select type (vec)
      type is (rvector)
         self%data = alpha*self%data + beta*vec%data
      end select
      return
   end subroutine axpby

   !-------------------------------------
   !-----                           -----
   !-----     VECTOR TEST SUITE     -----
   !-----                           -----
   !-------------------------------------

   subroutine collect_real_vector_testsuite(testsuite)
      !> Collection of tests.
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("Vector norm", test_vector_norm), &
                  new_unittest("Vector dot product", test_vector_dot), &
                  new_unittest("Vector-scalar multiplication", test_vector_mult), &
                  new_unittest("Vector addition", test_vector_add), &
                  new_unittest("Vector subtraction", test_vector_sub), &
                  new_unittest("Vector zeroing", test_vector_zero), &
                  new_unittest("Matrix product direct", test_direct_krylov_matrix_product), &
                  new_unittest("Matrix product transpose", test_transpose_krylov_matrix_product), &
                  new_unittest("Matrix axpby (real matrices)", test_real_matrix_axpby), &
                  new_unittest("Matrix axpby (krylov matrices)", test_krylov_matrix_axpby) & 
                  ]
      return
   end subroutine collect_real_vector_testsuite

   subroutine test_vector_norm(error)
      ! --> Error type to be returned.
      type(error_type), allocatable, intent(out) :: error

      ! --> Test vector.
      class(rvector), allocatable :: x
      real(kind=wp) :: alpha

      ! --> Initialize vector.
      x = rvector(); call random_number(x%data)
      ! --> Compute the vector norm.
      alpha = x%norm()
      ! --> Check the correctness.
      call check(error, is_close(alpha, norm2(x%data)))

      return
   end subroutine test_vector_norm

   subroutine test_vector_add(error)
      !> Error type to be returned.
      type(error_type), allocatable, intent(out) :: error
      !> Test vectors.
      class(rvector), allocatable :: x, y, z

      ! --> Initialize vectors.
      x = rvector(); call random_number(x%data)
      y = rvector(); call random_number(y%data)
      z = rvector(); z = x
      ! --> Vector addition.
      call z%add(y)
      ! --> Check result.
      call check(error, norm2(z%data - (x%data + y%data)) < rtol)

      return
   end subroutine test_vector_add

   subroutine test_vector_sub(error)
      !> Error type to be returned.
      type(error_type), allocatable, intent(out) :: error
      !> Test vectors.
      class(rvector), allocatable :: x, y, z

      ! --> Initialize vectors.
      x = rvector(); call random_number(x%data)
      y = rvector(); call random_number(y%data)
      z = rvector(); z = x
      ! --> Vector addition.
      call z%sub(y)
      ! --> Check result.
      call check(error, norm2(z%data - (x%data - y%data)) < rtol)

      return
   end subroutine test_vector_sub

   subroutine test_vector_zero(error)
      !> Error type to be returned.
      type(error_type), allocatable, intent(out) :: error
      !> Test vector.
      class(rvector), allocatable :: x

      ! --> Initialize vector.
      x = rvector(); call random_number(x%data)
      ! --> Zero-out vector.
      call x%zero()
      ! --> Check result.
      call check(error, x%norm() < rtol)

      return
   end subroutine test_vector_zero

   subroutine test_vector_dot(error)
      !> Error type to be returned.
      type(error_type), allocatable, intent(out) :: error

      !> Test vectors.
      class(rvector), allocatable :: x, y
      real(kind=wp) :: alpha

      ! --> Initialize vectors.
      x = rvector(); call random_number(x%data)
      y = rvector(); call random_number(y%data)
      ! --> Compute dot product.
      alpha = x%dot(y)
      ! --> Check result.
      call check(error, is_close(alpha, dot_product(x%data, y%data)))

      return
   end subroutine test_vector_dot

   subroutine test_vector_mult(error)
      !> Error type to be returned.
      type(error_type), allocatable, intent(out) :: error
      !> Test vector.
      class(rvector), allocatable :: x, y

      ! --> Random data.
      x = rvector(); call random_number(x%data)
      y = rvector(); y = x
      ! --> Scalar multiplication.
      call x%scal(2.0_wp)
      ! --> Check result.
      call check(error, norm2(x%data - 2.0_wp*y%data) < rtol)

      return
   end subroutine test_vector_mult

   subroutine test_direct_krylov_matrix_product(error)
      !> Error type to be returned.
      type(error_type), allocatable, intent(out) :: error
      !> Test bases.
      class(rvector), dimension(:), allocatable :: A(:)
      class(rvector), dimension(:), allocatable :: C(:)
      !> Krylov subspace dimension.
      integer, parameter :: kdim1 = 3
      !> Number of columns in coefficien matrix
      integer, parameter :: kdim2 = 4
      !> Test matrices.
      real(kind=wp)               :: B(kdim1, kdim2)
      real(kind=wp)               :: Amat(test_size, kdim1)
      real(kind=wp)               :: Cmat(test_size, kdim2)
      !> Misc.
      integer :: i,j

      !> Initialize basis and copy data to matrix
      allocate(A(1:kdim1))
      do i = 1, size(A)
         call random_number(A(i)%data)
         Amat(:, i) = A(i)%data
      enddo
      allocate(C(1:kdim2))
      B = 0.0_wp
      do i = 1, size(A)
         do j = 1, size(C)
            call random_number(B(i,j))
         enddo
      enddo
      call mat_zero(C)
      !> Compute product
      call mat_mult(C,A,B)
      !> Copy data
      do i = 1, kdim2
         Cmat(:, i) = C(i)%data
      enddo
      call check(error, all_close(matmul(Amat, B), Cmat, rtol, atol) )
      return
   end subroutine test_direct_krylov_matrix_product

   subroutine test_transpose_krylov_matrix_product(error)
      !> Error type to be returned.
      type(error_type), allocatable, intent(out) :: error
      !> Test bases.
      class(rvector), dimension(:), allocatable :: A(:)
      class(rvector), dimension(:), allocatable :: B(:)
      !> Krylov subspace dimension.
      integer, parameter :: kdim1 = 3
      integer, parameter :: kdim2 = 4
      !> Test matrices.
      real(kind=wp)               :: C(kdim1, kdim2)
      real(kind=wp)               :: Amat(test_size, kdim1)
      real(kind=wp)               :: Bmat(test_size, kdim2)
      !> Misc.
      integer :: k

      !> Initialize bases and copy data to matrices
      allocate(A(1:kdim1))
      Amat = 0.0_wp
      do k = 1, size(A)
         call random_number(A(k)%data)
         Amat(:, k) = A(k)%data
      enddo
      allocate(B(1:kdim2))
      Bmat = 0.0_wp
      do k = 1, size(B)
         call random_number(B(k)%data)
         Bmat(:, k) = B(k)%data
      enddo
      C = 0.0_wp
      !> Compute product
      call mat_mult(C,A,B)
      call check(error, all_close(matmul(transpose(Amat), Bmat), C, rtol, atol) )
      return
   end subroutine test_transpose_krylov_matrix_product

   subroutine test_real_matrix_axpby(error)
      !> Error type to be returned.
      type(error_type), allocatable, intent(out) :: error
      !> Test matrices.
      real(kind=wp) , allocatable :: A(:,:)
      real(kind=wp) , allocatable :: B(:,:)
      ! factors
      real(kind=wp) :: alpha
      real(kind=wp) :: beta   
      !> Size
      integer, parameter :: kdim = 3
      !> Comparison.
      real(kind=wp) :: Z(test_size, kdim)
      !> Misc.
      integer :: i,j

      !> Initialize matrices
      allocate(A(1:test_size, 1:kdim))
      allocate(B(1:test_size, 1:kdim))
      do i = 1, test_size
         do j = 1, kdim
            call random_number(A(i,j))
            B(i,j) = -2.0*A(i,j)
         enddo
      enddo
      Z = 0.0_wp
      !> Compute sum
      call mat_axpby(A,2.0_wp,B,1.0_wp)
      call check(error, all_close(A, Z, rtol, atol) )
      return
   end subroutine test_real_matrix_axpby

   subroutine test_krylov_matrix_axpby(error)
      !> Error type to be returned.
      type(error_type), allocatable, intent(out) :: error
      !> Test matrices.
      class(rvector) , allocatable :: A(:)
      class(rvector) , allocatable :: B(:)
      ! factors
      real(kind=wp) :: alpha
      real(kind=wp) :: beta   
      !> Size
      integer, parameter :: kdim = 3
      !> Comparison.
      real(kind=wp) :: Amat(test_size, kdim)
      real(kind=wp) :: Zmat(test_size, kdim)
      !> Misc.
      integer :: i

      !> Initialize bases and copy data to matrices
      allocate(A(1:kdim))
      allocate(B(1:kdim))
      do i = 1, kdim
         call random_number(A(i)%data)
         call B(i)%axpby(0.0_wp,A(i),-2.0_wp)
      enddo
      Zmat = 0.0_wp
      !> Compute sum
      call mat_axpby(A,4.0_wp,B,2.0_wp)
      Amat = 0.0_wp
      !> Copy data to matrix
      do i = 1, kdim
         Amat(:, i) = A(i)%data
      enddo
      call check(error, all_close(Amat, Zmat, rtol, atol) )
      return
   end subroutine test_krylov_matrix_axpby

end module TestVector
