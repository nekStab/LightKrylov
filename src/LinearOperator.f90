module LinearOperator
  implicit none

  private

  ! -------------------------------------------------------
  ! -----                                            ------
  ! -----     ABSTRACT TYPE FOR GENERAL MATRICES     ------
  ! -----                                            ------
  ! -------------------------------------------------------

  type, abstract, public :: abstract_linop
   contains
     private
     ! Matrix-vector product.
     procedure(abstract_matvec), pass(self), deferred, public  :: matvec
     procedure(abstract_rmatvec), pass(self), deferred, public :: rmatvec
  end type abstract_linop

  abstract interface
     ! Interface for the matrix-vector product.
     subroutine abstract_matvec(self, vec_in, vec_out)
       use KrylovVector
       import abstract_linop
       class(abstract_linop) , intent(in)  :: self
       class(abstract_vector), intent(in)  :: vec_in
       class(abstract_vector), intent(out) :: vec_out
     end subroutine abstract_matvec

     ! Interface for the vector-matrix product.
     subroutine abstract_rmatvec(self, vec_in, vec_out)
       use KrylovVector
       import abstract_linop
       class(abstract_linop) , intent(in)  :: self
       class(abstract_vector), intent(in)  :: vec_in
       class(abstract_vector), intent(out) :: vec_out
     end subroutine abstract_rmatvec
  end interface

  !--------------------------------------------------------------
  !-----                                                    -----
  !-----     ABSTRACT TYPE FOR SYM. POS. DEF. OPERATORS     -----
  !-----                                                    -----
  !--------------------------------------------------------------

  type, extends(abstract_linop), abstract, public :: abstract_spd_linop
   contains
     private
  end type abstract_spd_linop

end module LinearOperator
