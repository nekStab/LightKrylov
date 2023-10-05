module TestMatrices
  use LightKrylov
  use TestVector
  use testdrive  , only : new_unittest, unittest_type, error_type, check
  use stdlib_math, only : all_close
  implicit none

  private

  public :: collect_real_matrix_testsuite, collect_abstract_linop_operations_testsuite

  !---------------------------------------
  !-----     GENERAL REAL MATRIX     -----
  !---------------------------------------
  type, extends(abstract_linop), public :: rmatrix
     double precision, dimension(3, 3) :: data = 0.0D+00
   contains
     private
     procedure, pass(self), public :: matvec => general_matvec
     procedure, pass(self), public :: rmatvec => general_rmatvec
  end type rmatrix

  !----------------------------------------
  !-----     SYM. POS. DEF MATRIX     -----
  !----------------------------------------
  type, extends(abstract_spd_linop), public :: spd_matrix
     double precision, dimension(3, 3) :: data = 0.0D+00
   contains
     private
     procedure, pass(self), public :: matvec => spd_matvec
     procedure, pass(self), public :: rmatvec => spd_matvec
  end type spd_matrix

contains

  !-------------------------------------------------------------------------------------
  !-----                                                                           -----
  !-----     DEFINITION OF THE TYPE-BOUND PROCEDURES FOR GENERAL REAL MATRICES     -----
  !-----                                                                           -----
  !-------------------------------------------------------------------------------------

  subroutine general_matvec(self, vec_in, vec_out)
    class(rmatrix)        , intent(in)  :: self
    class(abstract_vector), intent(in)  :: vec_in
    class(abstract_vector), intent(out) :: vec_out

    select type(vec_in)
    type is(rvector)
       select type(vec_out)
       type is(rvector)
          vec_out%data = matmul(self%data, vec_in%data)
       end select
    end select
  end subroutine general_matvec

  subroutine general_rmatvec(self, vec_in, vec_out)
    class(rmatrix), intent(in)          :: self
    class(abstract_vector), intent(in)  :: vec_in
    class(abstract_vector), intent(out) :: vec_out

    select type(vec_in)
    type is(rvector)
       select type(vec_out)
       type is(rvector)
          vec_out%data = matmul(transpose(self%data), vec_in%data)
       end select
    end select
  end subroutine general_rmatvec

  !-----------------------------------------------------------------------------------------------
  !-----                                                                                     -----
  !-----     DEFINITION OF THE TYPE-BOUND PROCEDURES FOR SYM. POS. DEF. LINEAR OPERATORS     -----
  !-----                                                                                     -----
  !-----------------------------------------------------------------------------------------------

  subroutine spd_matvec(self, vec_in, vec_out)
    class(spd_matrix)     , intent(in)  :: self
    class(abstract_vector), intent(in)  :: vec_in
    class(abstract_vector), intent(out) :: vec_out

    select type(vec_in)
    type is(rvector)
       select type(vec_out)
       type is(rvector)
          vec_out%data = matmul(self%data, vec_in%data)
       end select
    end select
    return
  end subroutine spd_matvec

  !----------------------------------------------------------------------------------
  !-----                                                                        -----
  !-----     DEFINITION OF THE TYPE-BOUND PROCEDURES FOR HERMITIAN MATRICES     -----
  !-----                                                                        -----
  !----------------------------------------------------------------------------------

  !--------------------------------------------------------
  !-----                                              -----
  !-----     TEST SUITE FOR GENERAL REAL MATRICES     -----
  !-----                                              -----
  !--------------------------------------------------------

  subroutine collect_real_matrix_testsuite(testsuite)
    !> Collection of tests.
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [&
         new_unittest("Matrix-vector product",  test_real_matvec), &
         new_unittest("Tranpose matrix-vector product", test_real_rmatvec) &
         ]

    return
  end subroutine collect_real_matrix_testsuite

  subroutine test_real_matvec(error)
    !> Error type to be returned.
    type(error_type), allocatable, intent(out) :: error
    !> Test vector.
    class(rvector), allocatable :: x, y
    !> Test matrix.
    class(rmatrix), allocatable :: A
    !> Misc.
    double precision, dimension(3)    :: xdata, ydata
    double precision, dimension(3, 3) :: Adata

    ! --> Initialize vector.
    call random_number(xdata) ; x = rvector(xdata)
    y = rvector()
    ! --> Initialize matrix.
    call random_number(Adata) ; A = rmatrix(Adata)
    ! --> Compute matrix-vector product.
    call A%matvec(x, y) ; ydata = matmul(Adata, xdata)
    ! --> Check result.
    call check(error, all_close(y%data, ydata), .true.)

    return
  end subroutine test_real_matvec

  subroutine test_real_rmatvec(error)
    !> Error type to be returned.
    type(error_type), allocatable, intent(out) :: error
    !> Test vector.
    class(rvector), allocatable :: x, y
    !> Test matrix.
    class(rmatrix), allocatable :: A
    !> Misc.
    double precision, dimension(3)    :: xdata, ydata
    double precision, dimension(3, 3) :: Adata

    ! --> Initialize vector.
    call random_number(xdata) ; x = rvector(xdata)
    y = rvector()
    ! --> Initialize matrix.
    call random_number(Adata) ; A = rmatrix(Adata)
    ! --> Compute transpose matrix-vector product.
    call A%rmatvec(x, y) ; ydata = matmul(transpose(Adata), xdata)
    ! --> Check result.
    call check(error, all_close(y%data, ydata), .true.)
    return
  end subroutine test_real_rmatvec

  !----------------------------------------------------------
  !-----                                                -----
  !-----     TEST SUITE FOR SYM. POS. DEF. MATRICES     -----
  !-----                                                -----
  !----------------------------------------------------------

  subroutine collect_spd_matrix_testsuite(testsuite)
    !> Collection of tests.
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    return
  end subroutine collect_spd_matrix_testsuite

  !----------------------------------------------------------------
  !-----                                                      -----
  !-----     TEST SUITE FOR OPERATIONS ON ABSTRACT LINOPS     -----
  !-----                                                      -----
  !----------------------------------------------------------------

  subroutine collect_abstract_linop_operations_testsuite(testsuite)
    !> Collection of tests.
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
         new_unittest("Scaled linear operator", test_scaled_linop) &
    ]
    return
  end subroutine collect_abstract_linop_operations_testsuite

  subroutine test_scaled_linop(error)
    !> Error type to be returned.
    type(error_type), allocatable, intent(out) :: error
    !> Test matrix.
    class(rmatrix), allocatable :: A
    !> Test vector.
    class(rvector), allocatable :: x
    class(rvector), allocatable :: y
    !> Scaling factor.
    real(kind=wp) :: sigma
    !> Scaled linear operator.
    class(scaled_linop), allocatable :: B

    ! --> Initialize test matrix and vector.
    A = rmatrix() ; call random_number(A%data)
    x = rvector() ; call random_number(x%data)
    y = rvector() ; call y%zero()
    ! --> Scaled operator.
    sigma = 2.0D+00 ; B = scaled_linop(A, sigma)
    ! --> Compute scaled matrix-vector product.
    call B%matvec(x, y)
    write(*, *) y%data, sigma*matmul(A%data, x%data)
    ! --> Check error.
    call check(error, all_close(y%data, sigma*matmul(A%data, x%data)), .true.)

    return
  end subroutine test_scaled_linop

end module TestMatrices
