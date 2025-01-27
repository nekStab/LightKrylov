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
    use stdlib_linalg, only: is_hermitian, is_symmetric, diag, svd, eigh
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

    !------------------------------------------------
    !-----     OPTS TYPE FOR LINEAR SOLVERS     -----
    !------------------------------------------------

    type, extends(abstract_opts), public :: newton_sp_opts
        !! Options for Newton-Krylov fixed-point iteration.
        integer :: maxiter = 100
        !! Maximum number of Newton iterations (default = 100)
        logical :: ifbisect = .false.
        !! Bisection toggle to enforce residual reduction (default = .false.)
        integer :: maxstep_bisection = 5
        !! Maximum number of bisections (evaluations of F) for step selection (default = 5)
        !! Ignored if ifbisect = .false.
        logical :: if_print_metadata = .false.
        !! Print interation metadata on exit (default = .false.)
    end type
    
    type, extends(abstract_opts), public :: newton_dp_opts
        !! Options for Newton-Krylov fixed-point iteration.
        integer :: maxiter = 100
        !! Maximum number of Newton iterations (default = 100)
        logical :: ifbisect = .false.
        !! Bisection toggle to enforce residual reduction (default = .false.)
        integer :: maxstep_bisection = 5
        !! Maximum number of bisections (evaluations of F) for step selection (default = 5)
        !! Ignored if ifbisect = .false.
        logical :: if_print_metadata = .false.
        !! Print interation metadata on exit (default = .false.)
    end type
    

    type, extends(abstract_metadata), public :: newton_sp_metadata
        !! Metadata for Newton-Krylov fixed-point iteration.
        integer :: n_iter = 0
        !! Iteration counter
        integer :: eval_counter_record = 0
        !! System response evaluation counter:
        !! N.B.: For each of these evals the current residual and tolerance are recorded.
        real(sp), dimension(:), allocatable :: res
        !! Residual history
        real(sp), dimension(:), allocatable :: tol
        !! Tolerance history
        logical :: converged = .false.
        !! Convergence flag
        integer :: info = 0
        !! Copy of the information flag for completeness
    contains
        procedure, pass(self), public :: print => print_newton_sp
        procedure, pass(self), public :: reset => reset_newton_sp
        procedure, pass(self), public :: record => record_data_sp
    end type
   
    type, extends(abstract_metadata), public :: newton_dp_metadata
        !! Metadata for Newton-Krylov fixed-point iteration.
        integer :: n_iter = 0
        !! Iteration counter
        integer :: eval_counter_record = 0
        !! System response evaluation counter:
        !! N.B.: For each of these evals the current residual and tolerance are recorded.
        real(dp), dimension(:), allocatable :: res
        !! Residual history
        real(dp), dimension(:), allocatable :: tol
        !! Tolerance history
        logical :: converged = .false.
        !! Convergence flag
        integer :: info = 0
        !! Copy of the information flag for completeness
    contains
        procedure, pass(self), public :: print => print_newton_dp
        procedure, pass(self), public :: reset => reset_newton_dp
        procedure, pass(self), public :: record => record_data_dp
    end type
   

contains

    !------------------------------------------------------
    !-----     TYPE BOUND PROCEDURES FOR METADATA     -----
    !------------------------------------------------------

    subroutine print_newton_sp(self, reset_counters, verbose)
        class(newton_sp_metadata), intent(inout) :: self
        logical, optional, intent(in) :: reset_counters
        !! Reset all counters to zero after printing?
        logical, optional, intent(in) :: verbose
        !! Print the residual full residual history?
        ! internals
        integer :: i
        logical :: ifreset, ifverbose
        character(len=128) :: msg

        ifreset   = optval(reset_counters, .false.)
        ifverbose = optval(verbose, .false.)

        write(msg,'(A30,I20)') padr('Iterations: ', 30), self%n_iter
        call logger%log_message(msg, module=this_module, procedure='newton_metadata')
        if (ifverbose) then
            write(msg,'(14X,A15,2X,A15)') 'Residual', 'Tolerance'
            call logger%log_message(msg, module=this_module, procedure='newton_metadata')
            write(msg,'(A14,E15.8,2X,E15.8)') '   INIT:', self%res(1), self%tol(1)
            call logger%log_message(msg, module=this_module, procedure='newton_metadata')
            do i = 2, size(self%res) - 1
               write(msg,'(A,I4,A,E15.8,2X,E15.8)') '   Step ', i-1, ': ', self%res(i), self%tol(i)
               call logger%log_message(msg, module=this_module, procedure='newton_metadata')
            end do
            i = size(self%res)
            write(msg,'(A14,E15.8,2X,E15.8)') '   FINAL:', self%res(i), self%tol(i)
            call logger%log_message(msg, module=this_module, procedure='newton_metadata')
        else
            write(msg,'(A30,E20.8)') padr('Residual: ', 30), self%res(size(self%res))
            call logger%log_message(msg, module=this_module, procedure='newton_metadata')
        end if
        if (self%converged) then
            call logger%log_message('Status: CONVERGED', module=this_module, procedure='newton_metadata')
        else
            call logger%log_message('Status: NOT CONVERGED', module=this_module, procedure='newton_metadata')
        end if
        if (ifreset) call self%reset()
        return
    end subroutine print_newton_sp

    subroutine reset_newton_sp(self)
        class(newton_sp_metadata), intent(inout) :: self
        self%n_iter = 0
        self%eval_counter_record = 0
        self%converged = .false.
        self%info = 0
        if (allocated(self%res)) deallocate(self%res)
        if (allocated(self%tol)) deallocate(self%tol)
        return
    end subroutine reset_newton_sp

    subroutine record_data_sp(self, res, tol)
        class(newton_sp_metadata), intent(inout) :: self
        real(sp) :: res
        !! Residual of the current evaluation
        real(sp) :: tol
        !! Tolerance of the current evaluation
        if (.not.allocated(self%res)) then
            allocate(self%res(1))
            self%res(1) = res
        else
            self%res = [ self%res, res ]
        end if
        if (.not.allocated(self%tol)) then
            allocate(self%tol(1))
            self%tol(1) = tol
        else
            self%tol = [ self%tol, tol ]
        end if
        self%eval_counter_record = self%eval_counter_record + 1
        return
    end subroutine record_data_sp

    subroutine print_newton_dp(self, reset_counters, verbose)
        class(newton_dp_metadata), intent(inout) :: self
        logical, optional, intent(in) :: reset_counters
        !! Reset all counters to zero after printing?
        logical, optional, intent(in) :: verbose
        !! Print the residual full residual history?
        ! internals
        integer :: i
        logical :: ifreset, ifverbose
        character(len=128) :: msg

        ifreset   = optval(reset_counters, .false.)
        ifverbose = optval(verbose, .false.)

        write(msg,'(A30,I20)') padr('Iterations: ', 30), self%n_iter
        call logger%log_message(msg, module=this_module, procedure='newton_metadata')
        if (ifverbose) then
            write(msg,'(14X,A15,2X,A15)') 'Residual', 'Tolerance'
            call logger%log_message(msg, module=this_module, procedure='newton_metadata')
            write(msg,'(A14,E15.8,2X,E15.8)') '   INIT:', self%res(1), self%tol(1)
            call logger%log_message(msg, module=this_module, procedure='newton_metadata')
            do i = 2, size(self%res) - 1
               write(msg,'(A,I4,A,E15.8,2X,E15.8)') '   Step ', i-1, ': ', self%res(i), self%tol(i)
               call logger%log_message(msg, module=this_module, procedure='newton_metadata')
            end do
            i = size(self%res)
            write(msg,'(A14,E15.8,2X,E15.8)') '   FINAL:', self%res(i), self%tol(i)
            call logger%log_message(msg, module=this_module, procedure='newton_metadata')
        else
            write(msg,'(A30,E20.8)') padr('Residual: ', 30), self%res(size(self%res))
            call logger%log_message(msg, module=this_module, procedure='newton_metadata')
        end if
        if (self%converged) then
            call logger%log_message('Status: CONVERGED', module=this_module, procedure='newton_metadata')
        else
            call logger%log_message('Status: NOT CONVERGED', module=this_module, procedure='newton_metadata')
        end if
        if (ifreset) call self%reset()
        return
    end subroutine print_newton_dp

    subroutine reset_newton_dp(self)
        class(newton_dp_metadata), intent(inout) :: self
        self%n_iter = 0
        self%eval_counter_record = 0
        self%converged = .false.
        self%info = 0
        if (allocated(self%res)) deallocate(self%res)
        if (allocated(self%tol)) deallocate(self%tol)
        return
    end subroutine reset_newton_dp

    subroutine record_data_dp(self, res, tol)
        class(newton_dp_metadata), intent(inout) :: self
        real(dp) :: res
        !! Residual of the current evaluation
        real(dp) :: tol
        !! Tolerance of the current evaluation
        if (.not.allocated(self%res)) then
            allocate(self%res(1))
            self%res(1) = res
        else
            self%res = [ self%res, res ]
        end if
        if (.not.allocated(self%tol)) then
            allocate(self%tol(1))
            self%tol(1) = tol
        else
            self%tol = [ self%tol, tol ]
        end if
        self%eval_counter_record = self%eval_counter_record + 1
        return
    end subroutine record_data_dp


    !-------------------------------------------
    !-----     LAPACK MATRIX INVERSION     -----
    !-------------------------------------------

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
      sqrtmX = matmul(V, matmul(diag(lambda), transpose(V)))

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
      sqrtmX = matmul(V, matmul(diag(lambda), transpose(V)))

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
      sqrtmX = matmul(V, matmul(diag(lambda), conjg(transpose(V))))

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
      sqrtmX = matmul(V, matmul(diag(lambda), conjg(transpose(V))))

      return
    end subroutine

end module lightkrylov_utils
