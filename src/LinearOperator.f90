module lightkrylov_LinearOperator
  !! This module provides the base class `abstract_linop` from which all linear operators used in `LightKrylov` needs to be extended.
  !! To extend `abstract_linop`, you simply need to provide two type-bound procedures:
  !!
  !! - `matvec(self, vec_in, vec_out)` : Computes the matrix-vector product \(\mathbf{y}  = \mathbf{Ax} \).
  !! - `rmatvec(self, vec_in, vec_out)` : Compute the transpose matrix-vector product \( \mathbf{y} = \mathbf{A}^T \mathbf{x} \).
  !!
  !! It also provides extended types to define the identity operator (`identity_linop`), low-rank linear operator (`abstract_low_rank`), symmetric positive definite operators (`abstract_spd_linop`) and a few others.
  !!
  !!@note
  !! `LightKrylov` primarily aims at applications where the operator \(\mathbf{A}\) is not readily available but
  !! can only be accessed via a procedure computing the matrix-vector product \(\mathbf{Ax}\). As such, if you do
  !! have access to the matrix \(\mathbf{A}\) and can leverage its special structure to solve various problems,
  !! `LightKrylov` is very likely **not** the tool you want to use for your problem.
  !!@endnote
  !!
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
      !! Abstract type from which users need to extend in order to define an arbitrary linear operator to
      !! be used by `LightKrylov`. Upon extension, the user needs to provide two type-bound procedures:
      !!
      !! - `matvec(self, x)`
      !! - `rmatvec(self, x)`
   contains
      private
      procedure(abstract_matvec), pass(self), deferred, public  :: matvec
      !! Definition of the matrix-vector product \(\mathbf{y} = \mathbf{Ax}\).
      procedure(abstract_rmatvec), pass(self), deferred, public :: rmatvec
      !! Definition of the matrix-vector product \(\mathbf{y} = \mathbf{A}^T \mathbf{x}\).
   end type abstract_linop

   abstract interface
      subroutine abstract_matvec(self, vec_in, vec_out)
        !! Abstract interface defining the matrix-vector \(\mathbf{y} = \mathbf{Ax}\).
         use lightkrylov_AbstractVector
         import abstract_linop
         class(abstract_linop), intent(in)  :: self
         !! Linear operator `A`.
         class(abstract_vector), intent(in)  :: vec_in
         !! Vector to be multiplied by \(\mathbf{A}\).
         class(abstract_vector), intent(out) :: vec_out
         !! Result of the matrix-vector product.
      end subroutine abstract_matvec

      subroutine abstract_rmatvec(self, vec_in, vec_out)
        !! Abstract interface defining the transpose matrix-vector product \(\mathbf{y} = \mathbf{A}^T \mathbf{x}\).
         use lightkrylov_AbstractVector
         import abstract_linop
         class(abstract_linop), intent(in)  :: self
         !! Linear operator `A`.
         class(abstract_vector), intent(in)  :: vec_in
         !! Vector to be multiplied by \(\mathbf{A}^T\).
         class(abstract_vector), intent(out) :: vec_out
         !! Result of the transpose matrix-vector product.
      end subroutine abstract_rmatvec
   end interface

   !--------------------------------------------------------------
   !-----                                                    -----
   !-----     ABSTRACT TYPE FOR SYM. POS. DEF. OPERATORS     -----
   !-----                                                    -----
   !--------------------------------------------------------------

   type, extends(abstract_linop), abstract, public :: abstract_spd_linop
      !! Extended abstract type to define symmetric positive definite linear operators.
   contains
      private
   end type abstract_spd_linop

   !-------------------------------------------------------
   !-----                                             -----
   !-----     ABSTRACT TYPE FOR LOW-RANK MATRICES     -----
   !-----                                             -----
   !-------------------------------------------------------

   type, extends(abstract_linop), abstract, public :: abstract_lowrank_linop
      !! Extended abstract type to define low-rank linear operators.
   contains
      private
   end type abstract_lowrank_linop

   !-------------------------------------
   !-----                           -----
   !-----     IDENTITY OPERATOR     -----
   !-----                           -----
   !-------------------------------------

   type, extends(abstract_linop), public :: identity_linop
      !! Definition of the identity map.
   contains
      private
      procedure, pass(self), public :: matvec => identity_matvec
      !! Returns \(\mathbf{x} = \mathbf{I} \mathbf{x}\).
      procedure, pass(self), public :: rmatvec => identity_matvec
      !! Returns \(\mathbf{x} = \mathbf{I} \mathbf{x}\).
   end type identity_linop

   !-----------------------------------
   !-----                         -----
   !-----     SCALED OPERATOR     -----
   !-----                         -----
   !-----------------------------------

   type, extends(abstract_linop), public :: scaled_linop
      !! Defines a scaled linear operator \( \mathbf{B} = \sigma \mathbf{A} \) with \( \sigma \in \mathbb{R} \).
      class(abstract_linop), allocatable :: A
      !! Base linear operator to be scaled.
      real(kind=wp)         :: sigma
      !! Scaling factor.
   contains
      private
      procedure, pass(self), public :: matvec => scaled_matvec
      !! Compute the scaled matrix-vector product \( \mathbf{y} = \sigma \mathbf{Ax} \).
      procedure, pass(self), public :: rmatvec => scaled_rmatvec
      !! Compute the scaled transpose matrix-vector product \( \mathbf{y} = \sigma \mathbf{A}^T \mathbf{x} \).
   end type scaled_linop

   !----------------------------------
   !-----                        -----
   !-----     AXBPY OPERATOR     -----
   !-----                        -----
   !----------------------------------

   type, extends(abstract_linop), public :: axpby_linop
      !! Extended type to define addition of linear operators \( \mathbf{C} = \alpha \mathrm{op}(\mathbf{A}) + \beta \mathrm{op}(\mathbf{B})\).
      class(abstract_linop), allocatable :: A
      !! Linear operator.
      class(abstract_linop), allocatable :: B
      !! Linear operator.
      real(kind=wp) :: alpha, beta
      !! Scalar multipliers.
      logical :: transA = .false., transB = .false.
      !! Logical controlling whether \(\mathrm{op}(\mathbf{A}) = \mathbf{A}\) or \(\mathrm{op}(\mathbf{A}) = \mathbf{A}^T\).
      !! Likewise for \(\mathbf{B}\).
   contains
      private
      procedure, pass(self), public :: matvec => axpby_matvec
      !! Definition of the matrix-vector product \( \mathbf{y} = \left( \alpha \mathbf{A} + \beta \mathbf{B} \right) \mathbf{x} \).
      procedure, pass(self), public :: rmatvec => axpby_rmatvec
      !! Definition of the matrix-vector product \( \mathbf{y} = \left( \alpha \mathbf{A}^T + \beta \mathbf{B}^T \right) \mathbf{x} \)  end type axpby_linop
   end type axpby_linop
contains

   !-------------------------------------------------------------------
   !-----                                                         -----
   !-----     TYPE-BOUND PROCEDURES FOR THE IDENTITY OPERATOR     -----
   !-----                                                         -----
   !-------------------------------------------------------------------

   subroutine identity_matvec(self, vec_in, vec_out)
      class(identity_linop), intent(in)  :: self
      class(abstract_vector), intent(in)  :: vec_in
      class(abstract_vector), intent(out) :: vec_out
      call vec_out%axpby(0.0_wp, vec_in, 1.0_wp)
      return
   end subroutine identity_matvec

   !------------------------------------------------------------------
   !-----                                                        -----
   !-----      TYPE-BOUND PROCEDURES FOR THE SCALED OPERATOR     -----
   !-----                                                        -----
   !------------------------------------------------------------------

   subroutine scaled_matvec(self, vec_in, vec_out)
      class(scaled_linop), intent(in)     :: self
      !! Linear operator \( \mathbf{B} = \sigma \mathbf{A} \).
      class(abstract_vector), intent(in)  :: vec_in
      !! Vector to be multiplied by \(\mathbf{B}\).
      class(abstract_vector), intent(out) :: vec_out
      !! Result of the matrix-vector product \(\mathbf{y} = \sigma \mathbf{Ax}\).

      call self%A%matvec(vec_in, vec_out); call vec_out%scal(self%sigma)
      return
   end subroutine scaled_matvec

   subroutine scaled_rmatvec(self, vec_in, vec_out)
      class(scaled_linop), intent(in)     :: self
      !! Linear operator \( \mathbf{B} = \sigma \mathbf{A} \).
      class(abstract_vector), intent(in)  :: vec_in
      !! Vector to be multiplied by \( \mathbf{B} \).
      class(abstract_vector), intent(out) :: vec_out
      !! Result of the transpose matrix-vector product \(\mathbf{y} = \sigma \mathbf{A}^T \mathbf{x}\).

      call self%A%rmatvec(vec_in, vec_out); call vec_out%scal(self%sigma)
      return
   end subroutine scaled_rmatvec

   !-------------------------------------------------------------------
   !-----                                                         -----
   !-----     TYPE-BOUND PROCEDURES FOR AXPBY LINEAR OPERATOR     -----
   !-----                                                         -----
   !-------------------------------------------------------------------

   subroutine axpby_matvec(self, vec_in, vec_out)
      class(axpby_linop), intent(in)  :: self
      !! Linear operator \( \mathbf{C} = \alpha \mathbf{A} + \beta \mathbf{B} \).
      class(abstract_vector), intent(in)  :: vec_in
      !! Vector to be multiplied by \( \mathbf{C} \).
      class(abstract_vector), intent(out) :: vec_out
      !! Result of the matrix-vector product \(\mathbf{y} = \mathbf{Cx}\).

      ! Working array.
      class(abstract_vector), allocatable :: wrk

      ! --> Allocate working array.
      allocate (wrk, source=vec_in)

      ! --> w = A @ x
      if (self%transA) then
         call self%A%rmatvec(vec_in, wrk)
      else
         call self%A%matvec(vec_in, wrk)
      end if

      ! --> y = B @ x
      if (self%transB) then
         call self%B%rmatvec(vec_in, vec_out)
      else
         call self%B%matvec(vec_in, vec_out)
      end if

      ! --> y = alpha*w + beta*y
      call vec_out%axpby(self%beta, wrk, self%alpha)

      return
   end subroutine axpby_matvec

   subroutine axpby_rmatvec(self, vec_in, vec_out)
      class(axpby_linop), intent(in)  :: self
      !! Linear operator \( \mathbf{C} = \alpha \mathbf{A} + \beta \mathbf{B} \).
      class(abstract_vector), intent(in)  :: vec_in
      !! Vector to be multiplied by \(\mathbf{C}\).
      class(abstract_vector), intent(out) :: vec_out
      !! Result of the transpose matrix-vector product \(\mathbf{y} = \mathbf{C}^T \mathbf{x}\).

      ! Working array.
      class(abstract_vector), allocatable :: wrk

      ! --> w = A @ x
      if (self%transA) then
         call self%A%matvec(vec_in, wrk)
      else
         call self%A%rmatvec(vec_in, wrk)
      end if

      ! --> y = B @ x
      if (self%transB) then
         call self%B%matvec(vec_in, vec_out)
      else
         call self%B%rmatvec(vec_in, vec_out)
      end if

      ! --> y = alpha*w + beta*y
      call vec_out%axpby(self%beta, wrk, self%alpha)

      return
   end subroutine axpby_rmatvec

end module lightkrylov_LinearOperator
