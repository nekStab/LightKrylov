module LinearOperator
  use AbstractVector
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
       use AbstractVector
       import abstract_linop
       class(abstract_linop) , intent(in)  :: self
       class(abstract_vector), intent(in)  :: vec_in
       class(abstract_vector), intent(out) :: vec_out
     end subroutine abstract_matvec

     ! Interface for the vector-matrix product.
     subroutine abstract_rmatvec(self, vec_in, vec_out)
       use AbstractVector
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

  !-------------------------------------------------------
  !-----                                             -----
  !-----     ABSTRACT TYPE FOR LOW-RANK MATRICES     -----
  !-----                                             -----
  !-------------------------------------------------------

  type, extends(abstract_linop), abstract, public :: abstract_lowrank_linop
   contains
     private
  end type abstract_lowrank_linop

  !-------------------------------------
  !-----                           -----
  !-----     IDENTITY OPERATOR     -----
  !-----                           -----
  !-------------------------------------

  type, extends(abstract_linop), public :: Identity_operator
   contains
     private
     procedure, pass(self), public :: matvec  => identity_matvec
     procedure, pass(self), public :: rmatvec => identity_matvec
  end type Identity_operator

contains

  !-------------------------------------------------------------------
  !-----                                                         -----
  !-----     TYPE-BOUND PROCEDURES FOR THE IDENTITY OPERATOR     -----
  !-----                                                         -----
  !-------------------------------------------------------------------

  subroutine identity_matvec(self, vec_in, vec_out)
    class(identity_operator), intent(in)  :: self
    class(abstract_vector)  , intent(in)  :: vec_in
    class(abstract_vector)  , intent(out) :: vec_out
    ! /!\ NOTE : This needs to be improved. It is a simple hack but I ain't happy with it.
    call vec_out%zero() ; call vec_out%add(vec_in)
    return
  end subroutine identity_matvec
  
end module LinearOperator
