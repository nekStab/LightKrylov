module lightkrylov_AbstractVector
  !! This module provides the base class `abstract_vector` from which all Krylov vectors need to be derived.
  !! In order to use `LightKrylov`, you need to extend this class for your own definition of a vector and
  !! provide implementations for the following type-bound procedures:
  !!
  !! - `zero(self)` : A function zeroing-out the vector.
  !! - `scal(self, alpha)` : A function computing **in-place** the product \( \mathbf{x} = \alpha \mathbf{x} \) with \(\alpha \in \mathbb{R} \).
  !! - `axpby(self, alpha, vec, beta)` : A function computing **in-place** the result of \( \mathbf{x} = \alpha \mathbf{x} + \beta \mathbf{y} \) with \( \mathbf{x} \) and \( \mathbf{y} \) two similar instances derived from `abstract_vector` and \( \alpha, \beta \in \mathbb{R} \).
  !! - `dot(self, vec)` : A function computing the inner product between two similar instances derived from `abstract_vector`, i.e. \( \alpha = \mathbf{x}^T \mathbf{y} \) (with \( \alpha \in \mathbb{R} \)).
  !!
  !! Once these type-bound procedures have been provided by the user, they will automatically be used to define the vector addition `add(self, v) = axpby(self, 1.0_wp, vec, 1.0_wp)`, subtraction `sub(self, vec) = axpby(self, 1.0_wp, vec, -1.0_wp)` and `norm(self) = sqrt(dot(self, self))`.
  !!
  !! This module also provides the following utility subroutines:
  !!
  !! - `get_vec(y, X, v)` : return \( \mathbf{y} = \mathbf{Xv} \) where
  !!     + `y` is an `abstract_vector`.
  !!     + `X` is an array of `abstract_vector`
  !!     + `v` is a standard vector in \( \mathbb{R}^n \) represented as a double precision fortran array.
  !! - `mat_mult_direct(C, A, B)` : return \( \mathbf{C} = \mathbf{AB} \) where
  !!     + `C` is an array of `abstract_vector` of size \(q\).
  !!     + `A` is an array of `abstract_vector` of size \(r\).
  !!     + `B` is a two-dimensional Fortran array of size \(r \times q\).
  !! - `mat_mult_transpose(C, A, B)` : return \( \mathbf{C} = \mathbf{A}^T \mathbf{B} \) where
  !!     + `C` is a two-dimensional Fortran array of size \(r \times q\).
  !!     + `A` is an array of `abstract_vector` of size \(r\).
  !!     + `B` is an array of `abstract_vector` of size \(q\).
  !! - `mat_axpby_realmat(A, alpha, B, beta)` : return **in-place** the sum \( \mathbf{A} = \alpha \mathbf{A} + \beta \mathbf{B} \) where
  !!     + `A` and `B` are two-dimensional Fortran array of compatible dimensions.
  !!     + `alpha` and `beta` are double precision floating point numbers.
  !! - `mat_axpby_avecmat(A, alpha, B, beta)` : return **in-place** the sum \( \mathbf{A} = \alpha \mathbf{A} + \beta \mathbf{B} \) where
  !!     + `A` and `B` are two arrays for `abstract_vector` with the same size (and type).
  !!     + `alpha` and `beta` are two double precision floating point numbers.
  !! - `mat_zero(A)` : return \( \mathbf{A} = \mathbf{0} \) where
  !!     + `A` is an array of `abstract_vector`.
  !! - `mat_copy(A, B)` : return \( \mathbf{A} = \mathbf{B} \) where
  !!     + `A` and `B` are two arrays of `abstract_vector` with the same size.
  !!
  !! @warning
  !! For the moment `LightKrylov` does not handle natively complex-valued vectors.
  !! If you need to work with complex vectors in \( \mathbb{C}^n \), you will have to rewrite your operations in terms of vectors in \( \mathbb{R}^{2n} \).
  !! @endwarning
   implicit none
   include "dtypes.h"

   private
   public :: get_vec
   ! Krylov Basis utilities
   public :: mat_mult, mat_axpby, mat_zero, mat_copy

   interface mat_mult
      ! Krylov Basis multiplication
      module procedure mat_mult_direct
      module procedure mat_mult_transpose
   end interface mat_mult

   interface mat_axpby
      ! Krylov Basis operations
      module procedure mat_axpby_realmat
      module procedure mat_axpby_avecmat
   end interface mat_axpby

   !---------------------------------------------------
   !-----     ABSTRACT VECTOR TYPE DEFINITION     -----
   !---------------------------------------------------

   type, abstract, public :: abstract_vector
      !! Abstract type from which users need to extend in order to define the vectors used by `LightKrylov`.
      !! Upon extension, the user needs to provide a handful of type-bound procedures:
      !!
      !! - `zero(self)`
      !! - `scal(self, alpha)`
      !! - `axpby(self, alpha, vec, beta)`
      !! - `dot(self, vec)`
   contains
      private
      procedure, pass(from), public :: copy
      generic, public               :: assignment(=) => copy
      !! Overload the assignment operator.
      procedure(abstract_zero), deferred, public :: zero
      !! Sets an `abstract_vector` to zero.
      procedure(abstract_rand), deferred, public :: rand
      !! Create a random `æbstract_vector.
      procedure(abstract_scal), deferred, public :: scal
      !! Compute the scalar-vector product.
      procedure(axpby_interface), deferred, pass(self), public :: axpby
      !! In-place computation of \( \mathbf{x} = \alpha \mathbf{x} + \beta \mathbf{y} \).
      procedure, pass(self), public :: add
      !! In-place vector addition.
      procedure, pass(self), public :: sub
      !! In-place vector subtraction.
      procedure(abstract_dot), pass(self), deferred, public :: dot
      !! Computes the dot product between two `abstract_vector`.
      procedure, pass(self), public                         :: norm
      !! Compute the norm of an `abstract_vector`.
   end type abstract_vector

   ! Abstract interfaces for the type-bound procedures.
   abstract interface
      subroutine abstract_zero(self)
        !! Abstract interface to zero-out a vector in-place.
         import abstract_vector
         class(abstract_vector), intent(inout) :: self
         !! Vector to be zeroed-out.
      end subroutine abstract_zero

      subroutine abstract_rand(self, ifnorm)
         !! Abstract interface to set vector to (normalized) random vector
         import abstract_vector
         class(abstract_vector), intent(inout) :: self
         logical, optional,      intent(in)    :: ifnorm
      end subroutine abstract_rand

      subroutine abstract_scal(self, alpha)
        !! Abstract interface for in-place scalar-vector multiplication.
         import abstract_vector, wp
         class(abstract_vector), intent(inout) :: self
         !! Input/Output vector.
         real(kind=wp), intent(in)             :: alpha
         !! Scalar multiplier.
      end subroutine abstract_scal

      subroutine axpby_interface(self, alpha, vec, beta)
        !! Abstract interface to compute \( \mathbf{x} = \alpha \mathbf{x} + \beta \mathbf{y} \)
         import abstract_vector, wp
         class(abstract_vector), intent(inout) :: self
         !! Input/Output vector.
         class(abstract_vector), intent(in)    :: vec
         !! Vector to be add/subtracted.
         real(kind=wp), intent(in)    :: alpha, beta
         !! Scalar multipliers.
      end subroutine axpby_interface

      real(kind=wp) function abstract_dot(self, vec) result(alpha)
        !! Abstract interface to compute the dot product between two `abstract_vector`.
         import abstract_vector, wp
         class(abstract_vector), intent(in) :: self, vec
         !! Vectors whose dot product will be computed.
      end function abstract_dot

   end interface

contains

   subroutine copy(out, from)
     !! This function copies the content in `from` to `out`.
      class(abstract_vector), intent(in)               :: from
     !! Vector we wish to copy from.
      class(abstract_vector), allocatable, intent(out) :: out
     !! Vector we wish to copy to.
      if (allocated(out)) deallocate (out)
      allocate (out, source=from)
      return
   end subroutine copy

   subroutine add(self, y)
     !! Add two `abstract_vector` in-place.
      class(abstract_vector), intent(inout) :: self
     !! Input/Output vector.
      class(abstract_vector), intent(in)    :: y
     !! Vector to add to `self`.
      call self%axpby(1.0_wp, y, 1.0_wp)
      return
   end subroutine add

   subroutine sub(self, y)
     !! Subtract two `abstract_vector` in-place.
      class(abstract_vector), intent(inout) :: self
     !! Input/Output vector.
      class(abstract_vector), intent(in)    :: y
     !! Vector to subtract from `self`.
      call self%axpby(1.0_wp, y, -1.0_wp)
      return
   end subroutine sub

   real(kind=wp) function norm(self) result(alpha)
     !! Compute the norm of an `abstract_vector`. Note that it is based
     !! on the subroutine `dot` which needs to be provided by the user.
      class(abstract_vector), intent(in) :: self
     !! Vector whose norm is computed.
      alpha = sqrt(self%dot(self))
   end function norm

   !------------------------------------------
   !-----                                -----
   !-----     MATRIX VECTOR PRODUCTS     -----
   !-----                                -----
   !------------------------------------------

   subroutine get_vec(y, X, v)
     !! Given `X` and `v`, this function returns \( \mathbf{y} = \mathbf{Xv} \) where
     !! `y` is an `abstract_vector`, `X` an array of `abstract_vector` and `v` a Fortran
     !! array containing the coefficients of the linear combination.
      class(abstract_vector), allocatable, intent(out) :: y
     !! Output vector.
      class(abstract_vector), dimension(:), intent(in) :: X
     !! Krylov basis.
      real(kind=wp), dimension(:), intent(in) :: v
     !! Coordinates of the output vector y in the Krylov basis X.
      integer :: i
      ! Miscellaneous.

      ! Check sizes.
      if (size(X) .ne. size(v)) then
         write (*, *) "INFO : Krylov basis X and low-dimension vector v have different sizes."
         return
      end if
      ! Initialize output vector.
      if (allocated(y) .eqv. .false.) allocate (y, source=X(1)); call y%zero()
      ! Compute output vector.
      do i = 1, size(X)
         call y%axpby(1.0_wp, X(i), v(i))
      end do
      return
   end subroutine get_vec

   !------------------------------------------
   !-----                                -----
   !-----     MATRIX MATRIX PRODUCTS     -----
   !-----                                -----
   !------------------------------------------

   subroutine mat_mult_direct(C, A, B)
     !! Utility function to compute the product \( \mathbf{C} = \mathbf{AB} \) where
     !! `A` and `C` are arrays of `abstract_vector` and `B` a real-valued matrix.
      class(abstract_vector), intent(out) :: C(:)
     !! Array of `abstract_vector` containing the result of the operation.
      class(abstract_vector), intent(in)  :: A(:)
     !! Array of `abstract_vector` to be right-multiplied by `B`.
      real(kind=wp), intent(in)  :: B(:, :)
     !! Real-valued matrix to be left-multiplied by `A`.

      ! Local variables
      integer :: i, j

      ! Check sizes.
      if (size(A) .ne. size(B, 1)) then
         write (*, *) "INFO : mat_mult dimension error"
         write (*, *) "Abstract vector basis A and coefficient matrix B have incompatible sizes for the product C = A @ B."
         stop 1
      end if
      if (size(C) .ne. size(B, 2)) then
         write (*, *) "INFO : mat_mult dimension error"
         write (*, *) "Coefficient matrix B does not have the same amout of colums as output basis C for the product C = A @ B."
         stop 1
      end if

      ! Compute product column-wise
      do j = 1, size(B, 2)
         call C(j)%zero()
         do i = 1, size(B, 1)
            call C(j)%axpby(1.0_wp, A(i), B(i, j))
         enddo
      end do

      return
   end subroutine mat_mult_direct

   subroutine mat_mult_transpose(C, A, B)
      !! Utility function to compute the product \(\mathbf{C} = \mathbf{A}^T \mathbf{B} \) where
      !! `C` is a real-valued matrix while `A` and `B` are arrays of `abstract_vector`.
      class(abstract_vector), intent(in)  :: A(:)
       !! Array of `abstract_vector`.
      class(abstract_vector), intent(in)  :: B(:)
       !! Array of `abstract_vector`.
      real(kind=wp), intent(out) :: C(size(A), size(B))
       !! Inner-product matrix.

      ! Local variables
      integer :: i, j

      ! Compute product column-wise
      C = 0.0_wp
      do j = 1, size(B)
         do i = 1, size(A)
            C(i, j) = A(i)%dot(B(j))
         end do
      end do
      return
   end subroutine mat_mult_transpose

   !--------------------------------
   !-----                      -----
   !-----     MATRIX UTILS     -----
   !-----                      -----
   !--------------------------------

   subroutine mat_axpby_realmat(A, alpha, B, beta)
      !! Utility function to compute in-place the result \( \mathbf{A} = \alpha \mathbf{A} + \beta \mathbf{B} \)
      !! where `A` and `B` are real-valued matrices.
      !!@note
      !! This function might probably be replaced by a call to an appropriate `blas` function.
      !!@endnote
      real(kind=wp), intent(inout)  :: A(:, :)
      !! Input/Ouput matrix.
      real(kind=wp), intent(in)     :: B(:, :)
      !! Matrix to be added/subtracted to `A`.
      real(kind=wp), intent(in)     :: alpha, beta
      !! Scalar multipliers.

      ! local variables
      integer :: i, j

      ! Check size
      if (any(shape(A) .ne. shape(B))) then
         write (*, *) "INFO : Array sizes incompatible for summation. "
         stop 1
      end if
      do j = 1, size(A, 2)
         do i = 1, size(A, 1)
            A(i, j) = alpha*A(i, j) + beta*B(i, j)
         end do
      end do
   end subroutine mat_axpby_realmat

   subroutine mat_axpby_avecmat(A, alpha, B, beta)
      !! Utility function to compute in-place the result \( \mathbf{A} = \alpha \mathbf{A} + \beta \mathbf{B} \)
      !! where `A` and `B` are arrays of `abstract_vector`.
      class(abstract_vector), intent(inout)  :: A(:)
       !! Input/Output array of `abstract_vector`.
      class(abstract_vector), intent(in)     :: B(:)
       !! Array of `abstract_vector` to be added/subtracted to `A`.
      real(kind=wp), intent(in)     :: alpha, beta
       !! Scalar multipliers.

      ! local variables
      integer :: i

      ! Check size
      if (size(A) .ne. size(B)) then
         write (*, *) "INFO : Array sizes incompatible for summation. "
         stop 1
      end if
      do i = 1, size(A)
         call A(i)%axpby(alpha, B(i), beta)
      end do
   end subroutine mat_axpby_avecmat

   subroutine mat_zero(A)
      !! Initialize an array of `abstract_vector` to zero by repeatedly calling
      !! `A(i)%zero()` for `i = 1, size(A)`.
      class(abstract_vector), intent(inout)  :: A(:)
      !! Array of `abstract_vector` to be zeroed-out.

      ! local variables
      integer :: i

      do i = 1, size(A)
         call A(i)%zero()
      end do
   end subroutine mat_zero

   subroutine mat_copy(A, B)
      !! Utility function to copy an array of `abstract_vector` into another one.
      class(abstract_vector), intent(out)  :: A(:)
      !! Output array.
      class(abstract_vector), intent(in)   :: B(:)
      !! Input array from which we wish to copy the content.
      call mat_axpby(A, 0.0_wp, B, 1.0_wp)
   end subroutine mat_copy

end module lightkrylov_AbstractVector
