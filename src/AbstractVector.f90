module AbstractVector
  implicit none
  include "dtypes.h"

  private
  public :: get_vec

  !---------------------------------------------------
  !-----     ABSTRACT VECTOR TYPE DEFINITION     -----
  !---------------------------------------------------

  type, abstract, public :: abstract_vector
   contains
     private

  !> Basic operations.
     procedure, pass(from), public :: copy
     generic  , public             :: assignment(=) => copy

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
       real(kind=wp)         , intent(in)    :: alpha, beta
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
    if (allocated(out)) deallocate(out)
    allocate(out, source=from)
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
       write(*, *) "INFO : Krylov basis X and low-dimension vector v have different sizes."
       return
    endif
  !> Initialize output vector.
    if (allocated(y) .eqv. .false.) allocate(y, source=X(1)) ; call y%zero()
  !> Compute output vector.
    if (.not. allocated(wrk)) allocate(wrk, source=X(1))
    do i = 1, size(X)
       wrk = X(i)
       call wrk%scal(v(i))
       call y%add(wrk)
    enddo
    if (allocated(wrk)) deallocate(wrk)
    return
  end subroutine get_vec

end module AbstractVector