module KrylovVector
  implicit none

  private
  public :: get_vec

  type, abstract, public :: abstract_vector
   contains
     private
     ! ℓ₂ norm of the vector.
     procedure(abstract_norm), pass(self), deferred, public :: norm

     ! Copy vector.
     procedure, pass(from), public :: copy
     generic, public :: assignment(=) => copy

     ! Vector addition.
     procedure(abstract_add), pass(self), deferred, public :: add

     ! Vector subtraction.
     procedure(abstract_sub), pass(self), deferred, public :: sub

     ! Zero vector.
     procedure(abstract_zero), pass(self), deferred, public :: zero

     ! Dot product.
     procedure(abstract_dot_product), pass(self), deferred, public :: dot

     ! --> Scalar multiplication.
     procedure(abstract_scalar_multiplication), pass(self), deferred, public :: scalar_mult

  end type abstract_vector

  abstract interface
     ! Interface for the ℓ₂ norm.
     double precision function abstract_norm(self) result(out)
       import abstract_vector
       class(abstract_vector), intent(in) :: self
     end function abstract_norm

     ! Interface for vector addition.
     subroutine abstract_add(self, vec)
       import abstract_vector
       class(abstract_vector), intent(inout) :: self
       class(abstract_vector), intent(in)    :: vec
     end subroutine abstract_add

     ! Interface for vector subtraction.
     subroutine abstract_sub(self, vec)
       import abstract_vector
       class(abstract_vector), intent(inout) :: self
       class(abstract_vector), intent(in)    :: vec
     end subroutine abstract_sub

     ! Interface to zero-out the vector.
     subroutine abstract_zero(self)
       import abstract_vector
       class(abstract_vector), intent(inout) :: self
     end subroutine abstract_zero

     ! Interface for the dot product definition.
     double precision function abstract_dot_product(self, vec) result(alpha)
       import abstract_vector
       class(abstract_vector), intent(in) :: self
       class(abstract_vector), intent(in) :: vec
     end function abstract_dot_product

     ! Interface for scalar multiplication.
     subroutine abstract_scalar_multiplication(self, alpha)
       import abstract_vector
       class(abstract_vector), intent(inout) :: self
       double precision, intent(in)          :: alpha
     end subroutine abstract_scalar_multiplication

  end interface

contains

  subroutine copy(out, from)
    class(abstract_vector), intent(in) :: from
    class(abstract_vector), intent(out), allocatable :: out
    allocate(out, source=from)
  end subroutine copy

  subroutine get_vec(y, X, v)
    !> Output Krylov vector.
    class(abstract_vector), allocatable, intent(inout) :: y
    !> Krylov basis.
    class(abstract_vector), dimension(:), intent(in) :: X
    !> Coordinates of the output vector y in the Krylov basis X.
    double precision, dimension(:), intent(in) :: v
    !> Temporary Krylov vector.
    class(abstract_vector), allocatable :: wrk
    !> Miscellaneous.
    integer :: i

    ! --> Check sizes.
    if (size(X) .ne. size(v)) then
       write(*, *) "INFO : Krylov basis X and low-dimension vector v have different sizes."
       call exit()
    endif
    ! --> Initialize output vector.
    allocate(y, source=X(1)) ; call y%zero()
    ! --> Compute output vector.
    do i = 1, size(X)
       wrk = X(i) ; call wrk%scalar_mult(v(i))
       call y%add(wrk)
    enddo
    return
  end subroutine get_vec


end module KrylovVector
