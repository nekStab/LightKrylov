module LightKrylov_Logger
   use stdlib_optval, only : optval
   use stdlib_logger
   use stdlib_logger, only: logger => global_logger
   use stdlib_ascii, only : to_lower
   use stdlib_strings, only : chomp, replace_all
   ! Testdrive
   use testdrive, only: error_type, test_failed
   implicit none
   private

   logical, parameter :: exit_on_error = .true.
   logical, parameter :: exit_on_test_error = .true.

   public :: stop_error
   public :: check_info
   public :: check_test
   public :: logger

contains

   subroutine stop_error(msg, module, procedure)
      character(len=*),            intent(in)  :: msg
      !! The name of the procedure in which the call happens
      character(len=*), optional,  intent(in)  :: module
      !! The name of the module in which the call happens
      character(len=*), optional,  intent(in)  :: procedure
      !! The name of the procedure in which the call happens
      call check_info(-1, origin='', module=module, procedure=procedure, info_msg=msg)
      return
   end subroutine stop_error

   subroutine check_info(info, origin, module, procedure, info_msg)
      integer,                     intent(in)  :: info
      !! Informaion flag
      character(len=*),            intent(in)  :: origin
      !! The name of the subroutine from which the flag originates
      character(len=*), optional,  intent(in)  :: module
      !! The name of the module in which the call happens
      character(len=*), optional,  intent(in)  :: procedure
      !! The name of the procedure in which the call happens
      character(len=*), optional,  intent(in)  :: info_msg
      character*128                            :: str
      !! Optional extra message

      ! internals
      character*256 :: msg
      integer :: ierr

      str = optval(info_msg, '')

      ierr = 0
      if (info == 0) then
         ! Successful exit --> only log on debug
         write(msg, *) 'The subroutine "'//trim(origin)//'" returned successfully. ', trim(str)
         call logger%log_debug(trim(msg), module=module, procedure=procedure)
      else
         !
         !   LAPACK
         !
         if (trim(to_lower(origin)) == 'getref') then
            ! GETREF
            if (info < 0 ) then
               write(msg, *) "The ", -info, "-th argument has illegal value. ", trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            else
               write(msg, *) "U(", info, ",", info, ") is exactly zero. The factorization ", &
                           & "has been completed but the factor U is exactly singular. ", &
                           & "Division by zero will occur if used to solve Ax=b. ", trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'getri') then
            ! GETRI
            if (info < 0 ) then
               write(msg, *) "The ", -info, "-th argument has illegal value. ", trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            else
               write(msg, *) "U(", info, ",", info, ") is exactly zero. ", &
                           & "The matrix is singular and its inverse cannot be computed. ", trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'geev') then
            ! GEEV
            if (info < 0) then
               write(msg, *) "The ", -info, "-th argument has illegal value. ", trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            else 
               write(msg, *) "The QR alg. failed to compute all of the eigenvalues.", &
                           & "No eigenvector has been computed. ", trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'syev') then
            ! SYEV
            if (info < 0) then
               write(msg, *) "The ", -info, "-th argument has illegal value. ", trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            else 
               write(msg, *) "The QR alg. failed to compute all of the eigenvalues.", &
                           & "No eigenvector has been computed. ", trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'heev') then
            ! HEEV
            if (info < 0) then
               write(msg, *) "The ", -info, "-th argument has illegal value. ", trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            else 
               write(msg, *) "The QR alg. failed to compute all of the eigenvalues.", &
                           & "No eigenvector has been computed. ", trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'gels') then
            ! GELS
            if (info < 0) then
               write(msg, *) "The ", -info, "-th argument has illegal value. ", trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            else
               write(msg, *) "Undocumented error. ", trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'gees') then
            ! GEES
            if (info < 0) then
               write(msg, *) "The ", -info, "-th argument has illegal value. ", trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            else
               write(msg, *) "The QR alg. failed to compute all of the eigenvalues.", &
                           & "No eigenvector has been computed. ", trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'trsen') then
            ! GEES
            if (info < 0) then
               write(msg, *) "The ", -info, "-th argument has illegal value. ", trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            else if (info == 1) then
               write(msg, *) "The reordering of T failed because some eigenvalues are too", &
                           & "close to separate (the problem is very ill-conditioned); ", &
                           & "T may have been partially reordered, and WR and WI ", &
                           & "contain the eigenvalues in the same order as in T; S and", &
                           & "SEP (if requested) are set to zero. ", trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            else
               write(msg, *) "Undocumented error. ", trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
        !
         !   LightKrylov_Utils
         !
         else if (trim(to_lower(origin)) == 'sqrtm') then
            if (info == 1) then
               write(msg, *) 'The input matrix is singular to tolerance. The singular eigenvalues are set to zero. '
               call logger%log_warning(trim(msg), module=module, procedure=procedure)
            else if (info == -1) then
               write(msg, *) "The input matrix is not positive (semi-)definite. "
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            else
               write(msg, *) "Undocumented error. ", trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         !
         !   LightKrylov_BaseKrylov
         !
         else if (trim(to_lower(origin)) == 'orthogonalize_against_basis') then
            ! orthogonalization
            if (info == 1) then
               write(msg, *) 'Orthogonalization: The ', info, 'th input vector is numerically zero.'
               call logger%log_information(trim(msg), module=module, procedure=procedure)
            else if (info == -1) then
               write(msg, *) 'The input Krylov basis is not orthonormal.'
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            else if (info == -2) then
               write(msg, *) 'Orthogonalization: The last column of the input basis is zero.'
               call logger%log_warning(trim(msg), module=module, procedure=procedure)
            else
               write(msg, *) "Undocumented error. ", trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'double_gram_schmidt_step') then
            ! orthogonalization
            if (info == 1) then
               write(msg, *) 'Orthogonalization: The ', info, 'th input vector is numerically zero.'
               call logger%log_information(trim(msg), module=module, procedure=procedure)
            else if (info == -1) then
               write(msg, *) 'The input Krylov basis is not orthonormal.'
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            else if (info == -2) then
               write(msg, *) 'Orthogonalization: The last column of the input basis is zero.'
               call logger%log_warning(trim(msg), module=module, procedure=procedure)
            else
               write(msg, *) "Undocumented error. ", trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'qr') then
            ! qr
            if (info > 0) then
               write(msg, *) 'QR factorization: Colinear column detected in column', info
               call logger%log_information(trim(msg), module=module, procedure=procedure)
            else
               write(msg, *) "Undocumented error. ", trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'qr_pivot') then
            ! qr_pivot
            if (info > 0) then
               write(msg, *) 'QR factorization: Invariant subspace found after', info, 'steps.'
               call logger%log_information(trim(msg), module=module, procedure=procedure)
            else
               write(msg, *) "Undocumented error. ", trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'arnoldi') then
            ! arnoldi
            if (info > 0) then
               write(msg, *) 'Arnoldi factorization: Invariant subspace computed after', info, 'iterations.'
               call logger%log_information(trim(msg), module=module, procedure=procedure)
            else
               write(msg, *) "Undocumented error. ", trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'lanczos_bidiagonalization') then
            ! lanczos_bidiagonalization
            if (info > 0) then
               write(msg, *) 'Lanczos Bidiagonalisation: Invariant subspace found after', info, 'steps.'
               call logger%log_information(trim(msg), module=module, procedure=procedure)
            else
               write(msg, *) "Undocumented error. ", trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'lanczos_tridiagonalization') then
            ! lanczos_tridiagonalization
            if (info > 0) then
               write(msg, *) 'Lanczos Tridiagonalisation: Invariant subspace found after', info, 'steps.'
               call logger%log_information(trim(msg), module=module, procedure=procedure)
            else
               write(msg, *) "Undocumented error. ", trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         !
         !   LightKrylov_IterativeSolvers
         !
         else if (trim(to_lower(origin)) == 'eigs') then
            ! GMRES
            if (info > 0) then
               write(msg, *) 'eigs iteration converged after', info, 'iterations'
               call logger%log_information(trim(msg), module=module, procedure=procedure)
            else
               write(msg, *) "Undocumented error. ", trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'eighs') then
            ! GMRES
            if (info > 0) then
               write(msg, *) 'eigs iteration converged after', info, 'iterations'
               call logger%log_information(trim(msg), module=module, procedure=procedure)
            else
               write(msg, *) "Undocumented error. ", trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'svds') then
            ! GMRES
            if (info > 0) then
               write(msg, *) 'svds iteration converged after', info, 'iterations'
               call logger%log_information(trim(msg), module=module, procedure=procedure)
            else
               write(msg, *) "Undocumented error. ", trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'gmres') then
            ! GMRES
            if (info > 0) then
               write(msg, *) 'GMRES iteration converged after', info, 'iterations'
               call logger%log_information(trim(msg), module=module, procedure=procedure)
            else
               write(msg, *) "Undocumented error. ", trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'cg') then
            ! CG
            if (info > 0) then
               write(msg, *) 'CG iteration converged after', info, 'iterations'
               call logger%log_information(trim(msg), module=module, procedure=procedure)
            else
               write(msg, *) "Undocumented error. ", trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         !
         !   LightKrylov_ExpmLib
         !
         else if (trim(to_lower(origin)) == 'kexpm') then
            ! Krylov Matrix Exponential
            if (info > 0) then
               write(msg, *) 'kexpm converged. Estimated error below tolerance using', info, 'Krylov vectors.'
               call logger%log_information(trim(msg), module=module, procedure=procedure)
            else if (info == -2) then
               write(msg, *) 'kexpm converged. Arnoldi iteration breakdown. Approximation is exact to arnoldi tolerance.'
               call logger%log_information(trim(msg), module=module, procedure=procedure)
            else if (info == -1) then
               write(msg, *) 'kexpm did not converge. Maximum number of Krylov vectors reached.'
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            write(msg, *) "Undocumented error. ", trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         !
         !   Default
         !
         else
            write(msg,*) 'subroutine "'//trim(origin)//'" returned with a non-zero error flag.'
            call logger%log_error(trim(msg), module=module, procedure=procedure, stat=info, errmsg=trim(str))
            ierr = -1
         end if
      end if ! info /= 0

      call error_handler(ierr)      
      
   end subroutine check_info

   subroutine error_handler(ierr)
      integer, intent(in) :: ierr

      if (ierr == 0) then
         return
      else
         if (exit_on_error) then
            write(*,*)
            write(*,*) 'A fatal error was encountered. Aborting calculation as per user directive.'
            write(*,*)
            STOP 1
         end if
      end if
      
   end subroutine error_handler

   subroutine check_test(error, test_name, info, eq, context)
      use face 
      type(error_type), allocatable, intent(inout) :: error
      character(len=*),              intent(in)    :: test_name
      character(len=*),    optional, intent(in)    :: info
      character(len=*),    optional, intent(in)    :: eq
      character(len=*),    optional, intent(in)    :: context
      character*128                                :: name
      
      ! internals
      character*128                  :: msg, info_, eq_
      character(len=*), parameter :: indent = repeat(" ", 7)
      character(len=4), dimension(4) :: substrings
      integer :: i

      info_ = optval(info, '')
      eq_   = optval(eq, '')

      name = trim(to_lower(test_name))
      substrings = ["_rsp", "_rdp", "_csp", "_cdp"]
      do i = 1, size(substrings)
         name = replace_all(name, substrings(i), "")
      end do
      name = replace_all(name, "test_", "")

      write(*, '(A33)', ADVANCE='NO') name
      write(*, '(A3)',  ADVANCE='NO') ' % '

      if (len(trim(info_)) == 0) then
         msg = eq_
      else
         if (len(info_) > 30) then
            msg = info_(:30) // eq_
         else
            msg = info_ // repeat(' ', 30 - len(trim(info_))) // eq_
         end if 
      end if
      write(*, '(A62)', ADVANCE='NO') msg

      if (allocated(error)) then
         print *, colorize('FAILED', color_fg='red')
         if (present(context)) then
            write(*, '(A)', ADVANCE='NO') trim(context)
         end if
         write(*,*)
         write(*,*) 'The most recent test failed. Aborting calculation as per user directive.'
         write(*,*)
         STOP 1
      else
         print *, colorize('PASSED', color_fg='green')
      end if

   end subroutine check_test

end module LightKrylov_Logger
