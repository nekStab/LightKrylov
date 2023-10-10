module TestVector
  use LightKrylov
  use testdrive, only: new_unittest, unittest_type, error_type, check
  use stdlib_math, only : is_close

  implicit none

  private

  public :: collect_real_vector_testsuite, test_size

  integer, parameter :: test_size = 100

  type, extends(abstract_vector), public :: rvector
     double precision, dimension(test_size) :: data = 0.0D+00
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
    self%data = 0.0D+00
    return
  end subroutine zero

  double precision function dot(self, vec) result(alpha)
    class(rvector), intent(in)         :: self
    class(abstract_vector), intent(in) :: vec

    select type(vec)
    type is(rvector)
       alpha = dot_product(self%data, vec%data)
    end select
    return
  end function dot

  ! --> In-place scalar multiplication.
  subroutine scal(self, alpha)
    class(rvector), intent(inout) :: self
    double precision, intent(in) :: alpha
    self%data = self%data * alpha
    return
  end subroutine scal

  ! --> axpby interface
  subroutine axpby(self, alpha, vec, beta)
    class(rvector), intent(inout) :: self
    class(abstract_vector), intent(in) :: vec
    real(kind=wp)         , intent(in) :: alpha, beta

    select type(vec)
    type is(rvector)
       self%data = alpha * self%data + beta*vec%data
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

    testsuite = [&
         new_unittest("Vector norm", test_vector_norm), &
         new_unittest("Vector dot product", test_vector_dot), &
         new_unittest("Vector-scalar multiplication", test_vector_mult), &
         new_unittest("Vector addition", test_vector_add), &
         new_unittest("Vector subtraction", test_vector_sub), &
         new_unittest("Vector zeroing", test_vector_zero) &
         ]
    return
  end subroutine collect_real_vector_testsuite

  subroutine test_vector_norm(error)
    ! --> Error type to be returned.
    type(error_type), allocatable, intent(out) :: error

    ! --> Test vector.
    class(rvector), allocatable :: x
    double precision :: alpha

    ! --> Initialize vector.
    x = rvector() ; call random_number(x%data)
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
    x = rvector() ; call random_number(x%data)
    y = rvector() ; call random_number(y%data)
    z = rvector() ; z = x
    ! --> Vector addition.
    call z%add(y)
    ! --> Check result.
    call check(error, norm2(z%data - (x%data+y%data)) < rtol)

    return
  end subroutine test_vector_add

  subroutine test_vector_sub(error)
    !> Error type to be returned.
    type(error_type), allocatable, intent(out) :: error
    !> Test vectors.
    class(rvector), allocatable :: x, y, z

    ! --> Initialize vectors.
    x = rvector() ; call random_number(x%data)
    y = rvector() ; call random_number(y%data)
    z = rvector() ; z = x
    ! --> Vector addition.
    call z%sub(y)
    ! --> Check result.
    call check(error, norm2(z%data - (x%data-y%data)) < rtol)

    return
  end subroutine test_vector_sub

  subroutine test_vector_zero(error)
    !> Error type to be returned.
    type(error_type), allocatable, intent(out) :: error
    !> Test vector.
    class(rvector), allocatable :: x

    ! --> Initialize vector.
    x = rvector() ; call random_number(x%data)
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
    double precision :: alpha

    ! --> Initialize vectors.
    x = rvector() ; call random_number(x%data)
    y = rvector() ; call random_number(y%data)
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
    x = rvector() ; call random_number(x%data)
    y = rvector() ; y = x
    ! --> Scalar multiplication.
    call x%scal(2.0_wp)
    ! --> Check result.
    call check(error, norm2(x%data-2.0_wp*y%data) < rtol)

    return
  end subroutine test_vector_mult

end module TestVector
