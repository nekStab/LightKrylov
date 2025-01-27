module lightkrylov_utils
    !!  This module provides a set of utility functions used throughout `LightKrylov`.
    !!  It includes:
    !!
    !!  - `assert_shape`: Assert that the shape of the argument is the expected shape.
    !!  - `eig`: Compute the eigenvalue decomposition of a general matrix.
    !!  - `sqrtm`: Compute the non-negative square root of a symmetric positive definite matrix using its SVD.
    !!  - `sqrtm_eig`: Compute the non-negative square root of a symmetric positive definite matrix using its eigenvalue decomposition.
    !!  - `schur`: Compute the Schur factorization of a general square matrix.
    !!  - `ordschur`: Re-order the Schur factorization to have the selected eigenvalues in the upper left block.
    !!
    !!  Note that as the development of `stdlib` progresses, some of these functions
    !!  will be deprecated in favor of the `stdlib` implementations.

    !--------------------------------------------
    !-----     Standard Fortran Library     -----
    !--------------------------------------------
    use iso_fortran_env, only: output_unit
    use stdlib_optval, only: optval
    use stdlib_strings, only: padr
    use stdlib_linalg, only: is_hermitian, is_symmetric, diag, svd, eigh, hermitian
    ! Eigenvalue problem.
    use stdlib_linalg_lapack, only: geev
    ! Schur factorization.
    use stdlib_linalg_lapack, only: gees, trsen

    !-------------------------------
    !-----     LightKrylov     -----
    !-------------------------------
    ! Various constants.
    use lightkrylov_utils_bis
    use LightKrylov_Logger
    use LightKrylov_Constants

    implicit none
    private

    character(len=*), parameter :: this_module      = 'LK_Utils'
    character(len=*), parameter :: this_module_long = 'LightKrylov_Utils'

    public :: assert_shape, log2, abstract_metadata, abstract_opts
    ! Compute AX = XD for general dense matrices.
    public :: eig
    ! Compute matrix sqrt of input symmetric/hermitian positive definite matrix A
    public :: sqrtm
    public :: sqrtm_eig
    ! Re-orders the Schur factorization of A.
    public :: ordschur

    interface sqrtm_eig
        !!  Computes the non-negative square root of a symmetric positive definite matrix
        !!  using its eigenvalue decomposition.
        !!
        !!  ### Description
        !!
        !!  This interface provides methods to compute the non-negative square root of a symmetric
        !!  (hermitian) positive definite matrix \( \mathbf{A} \).
        !!
        !!  ### Syntax
        !!
        !!  `call sqrtm_eig(A, sqrtmA, info)`
        !!
        !!  ### Arguments
        !!  
        !!  `A`: Symmetric (hermitian) positive definite matrix whose non-negative square root
        !!  needs to be computed. It is an `intent(in)` argument.
        !!
        !!  `sqrtmA`: Non-negative square root of `A`. It has the same size, kind and type as `A`.
        !!  It is an `intent(out)` argument.
        !!
        !!  `info`: Information flag. It is an `intent(out)` argument.
        module procedure sqrtm_eig_rsp
        module procedure sqrtm_eig_rdp
        module procedure sqrtm_eig_csp
        module procedure sqrtm_eig_cdp
    end interface

contains

    ! NOTE: This function will be deprecated soon.
    subroutine sqrtm_eig_rsp(X, sqrtmX, info)
      !! Matrix-valued sqrt function for dense symmetric/hermitian positive (semi-)definite matrices
      real(sp), intent(in)  :: X(:,:)
      !! Matrix of which to compute the sqrt
      real(sp), intent(out) :: sqrtmX(size(X,1),size(X,1))
      !! Return matrix
      integer, intent(out) :: info
      !! Information flag

      ! internals
      real(sp) :: Xtmp(size(X, 1), size(X, 1))
      real(sp) :: lambda(size(X,1))
      real(sp) :: V(size(X,1), size(X,1))
      integer :: i
      character(len=256) :: msg

      info = 0

      ! Check if the matrix is symmetric
      if (.not. is_symmetric(X)) then
        write(msg,'(A)') "Input matrix is not symmetric."
        call stop_error(msg, module=this_module, procedure='sqrtm_rsp')
      end if

      ! Perform eigenvalue decomposition
      Xtmp = X ; call eigh(Xtmp, lambda, vectors=V, overwrite_a=.true.)

      ! Check if the matrix is positive definite (up to tol)
      do i = 1, size(lambda)
         if (abs(lambda(i)) .gt. 10*atol_sp ) then
            if (lambda(i) .gt. zero_rsp) then
               lambda(i) = sqrt(lambda(i))
            else
               lambda(i) = zero_rsp
               info = -1
            end if
         else
            lambda(i) = zero_rsp
            info = 1
         end if
      end do

      ! Reconstruct the square root matrix
      sqrtmX = matmul(V, matmul(diag(lambda), hermitian(V)))

      return
    end subroutine
    subroutine sqrtm_eig_rdp(X, sqrtmX, info)
      !! Matrix-valued sqrt function for dense symmetric/hermitian positive (semi-)definite matrices
      real(dp), intent(in)  :: X(:,:)
      !! Matrix of which to compute the sqrt
      real(dp), intent(out) :: sqrtmX(size(X,1),size(X,1))
      !! Return matrix
      integer, intent(out) :: info
      !! Information flag

      ! internals
      real(dp) :: Xtmp(size(X, 1), size(X, 1))
      real(dp) :: lambda(size(X,1))
      real(dp) :: V(size(X,1), size(X,1))
      integer :: i
      character(len=256) :: msg

      info = 0

      ! Check if the matrix is symmetric
      if (.not. is_symmetric(X)) then
        write(msg,'(A)') "Input matrix is not symmetric."
        call stop_error(msg, module=this_module, procedure='sqrtm_rdp')
      end if

      ! Perform eigenvalue decomposition
      Xtmp = X ; call eigh(Xtmp, lambda, vectors=V, overwrite_a=.true.)

      ! Check if the matrix is positive definite (up to tol)
      do i = 1, size(lambda)
         if (abs(lambda(i)) .gt. 10*atol_dp ) then
            if (lambda(i) .gt. zero_rdp) then
               lambda(i) = sqrt(lambda(i))
            else
               lambda(i) = zero_rdp
               info = -1
            end if
         else
            lambda(i) = zero_rdp
            info = 1
         end if
      end do

      ! Reconstruct the square root matrix
      sqrtmX = matmul(V, matmul(diag(lambda), hermitian(V)))

      return
    end subroutine
    subroutine sqrtm_eig_csp(X, sqrtmX, info)
      !! Matrix-valued sqrt function for dense symmetric/hermitian positive (semi-)definite matrices
      complex(sp), intent(in)  :: X(:,:)
      !! Matrix of which to compute the sqrt
      complex(sp), intent(out) :: sqrtmX(size(X,1),size(X,1))
      !! Return matrix
      integer, intent(out) :: info
      !! Information flag

      ! internals
      complex(sp) :: Xtmp(size(X, 1), size(X, 1))
      real(sp) :: lambda(size(X,1))
      complex(sp) :: V(size(X,1), size(X,1))
      integer :: i
      character(len=256) :: msg

      info = 0

      ! Check if the matrix is hermitian
      if (.not. is_hermitian(X)) then
        write(msg,'(A)') "Input matrix is not hermitian"
        call stop_error(msg, module=this_module, procedure='sqrtm_csp')
      end if

      ! Perform eigenvalue decomposition
      Xtmp = X ; call eigh(Xtmp, lambda, vectors=V, overwrite_a=.true.)

      ! Check if the matrix is positive definite (up to tol)
      do i = 1, size(lambda)
         if (abs(lambda(i)) .gt. 10*atol_sp ) then
            if (lambda(i) .gt. zero_rsp) then
               lambda(i) = sqrt(lambda(i))
            else
               lambda(i) = zero_rsp
               info = -1
            end if
         else
            lambda(i) = zero_rsp
            info = 1
         end if
      end do

      ! Reconstruct the square root matrix
      sqrtmX = matmul(V, matmul(diag(lambda), hermitian(V)))

      return
    end subroutine
    subroutine sqrtm_eig_cdp(X, sqrtmX, info)
      !! Matrix-valued sqrt function for dense symmetric/hermitian positive (semi-)definite matrices
      complex(dp), intent(in)  :: X(:,:)
      !! Matrix of which to compute the sqrt
      complex(dp), intent(out) :: sqrtmX(size(X,1),size(X,1))
      !! Return matrix
      integer, intent(out) :: info
      !! Information flag

      ! internals
      complex(dp) :: Xtmp(size(X, 1), size(X, 1))
      real(dp) :: lambda(size(X,1))
      complex(dp) :: V(size(X,1), size(X,1))
      integer :: i
      character(len=256) :: msg

      info = 0

      ! Check if the matrix is hermitian
      if (.not. is_hermitian(X)) then
        write(msg,'(A)') "Input matrix is not hermitian"
        call stop_error(msg, module=this_module, procedure='sqrtm_cdp')
      end if

      ! Perform eigenvalue decomposition
      Xtmp = X ; call eigh(Xtmp, lambda, vectors=V, overwrite_a=.true.)

      ! Check if the matrix is positive definite (up to tol)
      do i = 1, size(lambda)
         if (abs(lambda(i)) .gt. 10*atol_dp ) then
            if (lambda(i) .gt. zero_rdp) then
               lambda(i) = sqrt(lambda(i))
            else
               lambda(i) = zero_rdp
               info = -1
            end if
         else
            lambda(i) = zero_rdp
            info = 1
         end if
      end do

      ! Reconstruct the square root matrix
      sqrtmX = matmul(V, matmul(diag(lambda), hermitian(V)))

      return
    end subroutine

end module lightkrylov_utils
