module TestMatrices
  use LightKrylov
  use TestVector
  use testdrive, only: new_unittest, unittest_type, error_type, check
  implicit none

  private

  type, extends(abstract_linop), public :: rmatrix
     double precision, dimension(3, 3) :: data = 0.0D+00
   contains
     private
     procedure, pass(self), public :: matvec => general_matvec
     procedure, pass(self), public :: rmatvec => general_rmatvec
  end type rmatrix

contains

  !--------------------------------------------------------------------------------
  !-----                                                                      -----
  !-----     DEFINITION OF THE TYPE-BOUND PROCEDURES FOR GENERAL MATRICES     -----
  !-----                                                                      -----
  !--------------------------------------------------------------------------------

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

end module TestMatrices
