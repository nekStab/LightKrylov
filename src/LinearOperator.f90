module lightkrylov_LinearOperator
  use lightkrylov_AbstractVector
  implicit none
  include "dtypes.h"

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
       use lightkrylov_AbstractVector
       import abstract_linop
       class(abstract_linop) , intent(in)  :: self
       class(abstract_vector), intent(in)  :: vec_in
       class(abstract_vector), intent(out) :: vec_out
     end subroutine abstract_matvec

     ! Interface for the vector-matrix product.
     subroutine abstract_rmatvec(self, vec_in, vec_out)
       use lightkrylov_AbstractVector
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

  type, extends(abstract_linop), public :: identity_linop
   contains
     private
     procedure, pass(self), public :: matvec  => identity_matvec
     procedure, pass(self), public :: rmatvec => identity_matvec
  end type identity_linop

  !-----------------------------------
  !-----                         -----
  !-----     SCALED OPERATOR     -----
  !-----                         -----
  !-----------------------------------

  type, extends(abstract_linop), public :: scaled_linop
     !> Original operator.
     class(abstract_linop), allocatable :: A
     !> Scaling factor.
     real(kind=wp)         :: sigma
     contains
       private
       procedure, pass(self), public :: matvec  => scaled_matvec
       procedure, pass(self), public :: rmatvec => scaled_rmatvec
  end type scaled_linop

  !----------------------------------
  !-----                        -----
  !-----     AXBPY OPERATOR     -----
  !-----                        -----
  !----------------------------------

  type, extends(abstract_linop), public :: axpby_linop
     !> First operator.
     class(abstract_linop), allocatable :: A
     !> Second operator.
     class(abstract_linop), allocatable :: B
     !> Constants.
     real(kind=wp) :: alpha, beta
     !> Logicals.
     logical :: transA = .false., transB = .false.
   contains
     private
     procedure, pass(self), public :: matvec  => axpby_matvec
     procedure, pass(self), public :: rmatvec => axpby_rmatvec
  end type axpby_linop

contains

  !-------------------------------------------------------------------
  !-----                                                         -----
  !-----     TYPE-BOUND PROCEDURES FOR THE IDENTITY OPERATOR     -----
  !-----                                                         -----
  !-------------------------------------------------------------------

  subroutine identity_matvec(self, vec_in, vec_out)
    class(identity_linop) , intent(in)  :: self
    class(abstract_vector), intent(in)  :: vec_in
    class(abstract_vector), intent(out) :: vec_out
    ! /!\ NOTE : This needs to be improved. It is a simple hack but I ain't happy with it.
    call vec_out%axpby(0.0_wp, vec_in, 1.0_wp)
    return
  end subroutine identity_matvec

  !------------------------------------------------------------------
  !-----                                                        -----
  !-----      TYPE-BOUND PROCEDURES FOR THE SCALED OPERATOR     -----
  !-----                                                        -----
  !------------------------------------------------------------------

  subroutine scaled_matvec(self, vec_in, vec_out)
    !> Arguments.
    class(scaled_linop), intent(in)     :: self
    class(abstract_vector), intent(in)  :: vec_in
    class(abstract_vector), intent(out) :: vec_out
    !> Original matrix-vector product.
    call self%A%matvec(vec_in, vec_out)
    !> Scale the result.
    call vec_out%scal(self%sigma)
    return
  end subroutine scaled_matvec
  
  subroutine scaled_rmatvec(self, vec_in, vec_out)
    !> Arguments.
    class(scaled_linop), intent(in)     :: self
    class(abstract_vector), intent(in)  :: vec_in
    class(abstract_vector), intent(out) :: vec_out
    !> Original matrix-vector product.
    call self%A%rmatvec(vec_in, vec_out)
    !> Scale the result.
    call vec_out%scal(self%sigma)
    return
  end subroutine scaled_rmatvec

  !-------------------------------------------------------------------
  !-----                                                         -----
  !-----     TYPE-BOUND PROCEDURES FOR AXPBY LINEAR OPERATOR     -----
  !-----                                                         -----
  !-------------------------------------------------------------------

  subroutine axpby_matvec(self, vec_in, vec_out)
    !> Arguments.
    class(axpby_linop)    , intent(in)  :: self
    class(abstract_vector), intent(in)  :: vec_in
    class(abstract_vector), intent(out) :: vec_out
    !> Working array.
    class(abstract_vector), allocatable :: wrk

    ! --> Allocate working array.
    allocate(wrk, source=vec_in)

    ! --> w = A @ x
    if (self%transA) then
       call self%A%rmatvec(vec_in, wrk)
    else
       call self%A%matvec(vec_in, wrk)
    endif

    ! --> y = B @ x
    if (self%transB) then
       call self%B%rmatvec(vec_in, vec_out)
    else
       call self%B%matvec(vec_in, vec_out)
    endif

    ! --> y = alpha*w + beta*y
    call vec_out%axpby(self%beta, wrk, self%alpha)

    return
  end subroutine axpby_matvec

  subroutine axpby_rmatvec(self, vec_in, vec_out)
    !> Arguments.
    class(axpby_linop)    , intent(in)  :: self
    class(abstract_vector), intent(in)  :: vec_in
    class(abstract_vector), intent(out) :: vec_out
    !> Working array.
    class(abstract_vector), allocatable :: wrk

    ! --> w = A @ x
    if (self%transA) then
       call self%A%matvec(vec_in, wrk)
    else
       call self%A%rmatvec(vec_in, wrk)
    endif

    ! --> y = B @ x
    if (self%transB) then
       call self%B%matvec(vec_in, vec_out)
    else
       call self%B%rmatvec(vec_in, vec_out)
    endif

    ! --> y = alpha*w + beta*y
    call vec_out%axpby(self%beta, wrk, self%alpha)

    return
  end subroutine axpby_rmatvec

end module lightkrylov_LinearOperator
