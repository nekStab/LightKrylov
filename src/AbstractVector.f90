module lightkrylov_AbstractVector
   implicit none
   include "dtypes.h"

   private
   public :: get_vec
   !> Krylov Basis utilities
   public :: mat_mult, mat_axpby, mat_zero, mat_copy

   interface mat_mult
      !> Krylov Basis multiplication
      module procedure mat_mult_direct
      module procedure mat_mult_transpose
   end interface mat_mult
      
   interface mat_axpby
      !> Krylov Basis operations
      module procedure mat_axpby_realmat
      module procedure mat_axpby_avecmat
   end interface mat_axpby

   !---------------------------------------------------
   !-----     ABSTRACT VECTOR TYPE DEFINITION     -----
   !---------------------------------------------------

   type, abstract, public :: abstract_vector
   contains
      private

      !> Basic operations.
      procedure, pass(from), public :: copy
      generic, public             :: assignment(=) => copy

      procedure(abstract_zero), deferred, public :: zero

      !> Scalar-vector product.
      procedure(abstract_scal), deferred, public :: scal

      !> Vector-vector operations.
      procedure(axpby_interface), deferred, pass(self), public :: axpby
      procedure, pass(self), public :: add
      procedure, pass(self), public :: sub

      !> Reduction operations.
      procedure(abstract_dot), pass(self), deferred, public :: dot
      procedure, pass(self), public                         :: norm
   end type abstract_vector

   !> Abstract interfaces for the type-bound procedures.
   abstract interface
      !> Basic operations.
      subroutine abstract_zero(self)
         import abstract_vector
         class(abstract_vector), intent(inout) :: self
      end subroutine abstract_zero

      !> Scalar-vector product.
      subroutine abstract_scal(self, alpha)
         import abstract_vector, wp
         class(abstract_vector), intent(inout) :: self
         real(kind=wp), intent(in)             :: alpha
      end subroutine abstract_scal

      !> Vector-vector operations.
      subroutine axpby_interface(self, alpha, vec, beta)
         import abstract_vector, wp
         class(abstract_vector), intent(inout) :: self
         class(abstract_vector), intent(in)    :: vec
         real(kind=wp), intent(in)    :: alpha, beta
      end subroutine axpby_interface

      !> Reduction operations.
      real(kind=wp) function abstract_dot(self, vec) result(alpha)
         import abstract_vector, wp
         class(abstract_vector), intent(in) :: self, vec
      end function abstract_dot
   end interface

contains

   subroutine copy(out, from)
      class(abstract_vector), intent(in)               :: from
      class(abstract_vector), allocatable, intent(out) :: out
      if (allocated(out)) deallocate (out)
      allocate (out, source=from)
      return
   end subroutine copy

   subroutine add(self, y)
      class(abstract_vector), intent(inout) :: self
      class(abstract_vector), intent(in)    :: y
      !> Vector addition.
      call self%axpby(1.0_wp, y, 1.0_wp)
      return
   end subroutine add

   subroutine sub(self, y)
      class(abstract_vector), intent(inout) :: self
      class(abstract_vector), intent(in)    :: y
      !> Vector subtraction.
      call self%axpby(1.0_wp, y, -1.0_wp)
      return
   end subroutine sub

   real(kind=wp) function norm(self) result(alpha)
      class(abstract_vector), intent(in) :: self
      alpha = sqrt(self%dot(self))
   end function norm

   !------------------------------------------
   !-----                                -----
   !-----     MATRIX VECTOR PRODUCTS     -----
   !-----                                -----
   !------------------------------------------

   subroutine get_vec(y, X, v)
      !> Output vector.
      class(abstract_vector), allocatable, intent(out) :: y
      !> Krylov basis.
      class(abstract_vector), intent(in) :: X(:)
      !> Coordinates of the output vector y in the Krylov basis X.
      real(kind=wp), intent(in) :: v(:)
      !> Temporary Krylov vector.
      class(abstract_vector), allocatable :: wrk
      !> Miscellaneous.
      integer :: i
      !> Check sizes.
      if (size(X) .ne. size(v)) then
         write (*, *) "INFO : Krylov basis X and low-dimension vector v have different sizes."
         return
      end if
      !> Initialize output vector.
      if (allocated(y) .eqv. .false.) allocate (y, source=X(1)); call y%zero()
      !> Compute output vector.
      if (.not. allocated(wrk)) allocate (wrk, source=X(1))
      do i = 1, size(X)
         wrk = X(i)
         call wrk%scal(v(i))
         call y%add(wrk)
      end do
      if (allocated(wrk)) deallocate (wrk)
      return
   end subroutine get_vec

   !------------------------------------------
   !-----                                -----
   !-----     MATRIX MATRIX PRODUCTS     -----
   !-----                                -----
   !------------------------------------------

   subroutine mat_mult_direct(C,A,B)
      ! Compute the matrix product C = A @ B with
      !     C: abstract vector type Krylov basis :: size nxq
      !     A: abstract vector type Krylov basis :: size nxr
      !     B: real matrix                       :: size rxq
      class(abstract_vector) , intent(out) :: C(:)   ! result
      class(abstract_vector) , intent(in)  :: A(:)   ! krylov basis
      real(kind=wp)          , intent(in)  :: B(:,:) ! coefficient matrix
      !> Local variables
      class(abstract_vector) , allocatable :: wrk
      integer :: i
    
      !> Check sizes.
      if (size(A) .ne. size(B,1)) then
         write(*,*) "INFO : Abstract vector basis A and coefficient matrix B have incompatible sizes for the product A @ B."
         return
      endif
      allocate(wrk, source=A(1))
      !> Compute product column-wise
      do i = 1, size(B,2)
         call get_vec(wrk, A, B(:, i))
         call C(i)%axpby(0.0_wp, wrk, 1.0_wp)
      enddo
      deallocate(wrk)
      return
    end subroutine mat_mult_direct

    subroutine mat_mult_transpose(C,A,B)
       ! Compute the matrix product C = A.T @ B with
       !     C: real matrix                       :: size rxq
       !     A: abstract vector type Krylov basis :: size nxr
       !     B: abstract vector type Krylov basis :: size nxq
       class(abstract_vector) , intent(in)  :: A(:)   ! krylov basis
       class(abstract_vector) , intent(in)  :: B(:)   ! krylov basis
       real(kind=wp)          , intent(out) :: C(size(A),size(B)) ! result
       !> Local variables
       integer :: i, j

       !> Compute product column-wise
       C = 0.0_wp
       do i = 1, size(A)
          do j = 1, size(B)
             C(i,j) = A(i)%dot(B(j))
          enddo
       enddo
       return
    end subroutine mat_mult_transpose

    !--------------------------------
    !-----                      -----
    !-----     MATRIX UTILS     -----
    !-----                      -----
    !--------------------------------
 
    subroutine mat_axpby_realmat(A,alpha,B,beta)
       ! Compute the scaled sum of two matrices in-place
       !     alpha*A + beta*B
       ! with
       !     A,B: real matrices :: size nxr
       real(kind=wp) , intent(inout)  :: A(:,:)
       real(kind=wp) , intent(in)     :: B(:,:)
       real(kind=wp) , intent(in)     :: alpha
       real(kind=wp) , intent(in)     :: beta
       !> local variables
       integer :: i,j
       !> Check size 
       if (any(shape(A) .ne. shape(B))) then
          write(*, *) "INFO : Array sizes incompatible for summation. "
          stop 1
       endif
       do i = 1, size(A,1)
          do j = 1, size(A,2)
             A(i,j) = alpha*A(i,j) + beta*B(i,j)
          enddo
       enddo
    end subroutine mat_axpby_realmat
    
    subroutine mat_axpby_avecmat(A,alpha,B,beta)
       ! Compute the scaled sum of two matrices in-place
       !     alpha*A + beta*B
       ! with
       !     A,B: abstract_vector type Krylov bases :: size nxr
       class(abstract_vector) , intent(inout)  :: A(:)
       class(abstract_vector) , intent(in)     :: B(:)
       real(kind=wp)          , intent(in)     :: alpha
       real(kind=wp)          , intent(in)     :: beta
       !> local variables
       integer :: i
       !> Check size 
       if (size(A) .ne. size(B)) then
          write(*, *) "INFO : Array sizes incompatible for summation. "
          stop 1
       endif
       do i = 1, size(A)
          call A(i)%axpby(alpha, B(i), beta)
       enddo
    end subroutine mat_axpby_avecmat

    subroutine mat_zero(A)
       ! Initialise Krylov basis with zeros
       class(abstract_vector) , intent(inout)  :: A(:)
       !> local variables
       integer :: i
       do i = 1, size(A)
          call A(i)%zero()
       enddo
    end subroutine mat_zero

    subroutine mat_copy(A,B)
       ! Copy data from B to A
       class(abstract_vector) , intent(out)  :: A(:)
       class(abstract_vector) , intent(in)   :: B(:)
       call mat_axpby(A,0.0_wp,B,1.0_wp)
    end subroutine mat_copy

end module lightkrylov_AbstractVector
