module LightKrylov_Logger
   use stdlib_optval, only : optval
   use stdlib_logger
   use stdlib_logger, only: logger => global_logger
   use stdlib_ascii, only : to_lower
   implicit none
   private

   logical, parameter :: exit_on_error = .true.

   public :: check_info
   public :: logger

contains

   subroutine check_info(info, origin, module, procedure)
      integer,                     intent(in)  :: info
      !! Informaion flag
      character(len=*),            intent(in)  :: origin
      !! The name of the subroutine from which the flag originates
      character(len=*), optional,  intent(in)  :: module
      !! The name of the module in which the call happens
      character(len=*), optional,  intent(in)  :: procedure
      !! The name of the procedure in which the call happens

      ! internals
      character*256 :: msg
      integer :: ierr

      ierr = 0
      if (info == 0) then
         ! Successful exit --> only log on debug
         write(msg, *) 'The subroutine "'//trim(origin)//'" returned successfully.'
         call logger%log_debug(trim(msg), module=module, procedure=procedure)
      else
         !
         !   LAPACK
         !
         if (trim(to_lower(origin)) == 'getref') then
            ! GETREF
            if (info < 0 ) then
               write(msg, *) "The ", -info, "-th argument has illegal value."
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            else
               write(msg, *) "U(", info, ",", info, ") is exactly zero. The factorization ", &
                           & "has been completed but the factor U is exactly singular. ", &
                           & "Division by zero will occur if used to solve Ax=b."
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'getri') then
            ! GETRI
            if (info < 0 ) then
               write(msg, *) "The ", -info, "-th argument has illegal value."
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            else
               write(msg, *) "U(", info, ",", info, ") is exactly zero. ", &
                           & "The matrix is singular and its inverse cannot be computed."
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'geev') then
            ! GEEV
            if (info < 0) then
               write(msg, *) "The ", -info, "-th argument has illegal value."
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            else 
               write(msg, *) "The QR alg. failed to compute all of the eigenvalues.", &
                           & "No eigenvector has been computed."
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'syev') then
            ! SYEV
            if (info < 0) then
               write(msg, *) "The ", -info, "-th argument has illegal value."
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            else 
               write(msg, *) "The QR alg. failed to compute all of the eigenvalues.", &
                           & "No eigenvector has been computed."
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'heev') then
            ! HEEV
            if (info < 0) then
               write(msg, *) "The ", -info, "-th argument has illegal value."
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            else 
               write(msg, *) "The QR alg. failed to compute all of the eigenvalues.", &
                           & "No eigenvector has been computed."
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'gels') then
            ! GELS
            if (info < 0) then
               write(msg, *) "The ", -info, "-th argument has illegal value."
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            else
               write(msg, *) "Undocumented error."
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'gees') then
            ! GEES
            if (info < 0) then
               write(msg, *) "The ", -info, "-th argument has illegal value."
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            else
               write(msg, *) "The QR alg. failed to compute all of the eigenvalues.", &
                           & "No eigenvector has been computed."
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'trsen') then
            ! GEES
            if (info < 0) then
               write(msg, *) "The ", -info, "-th argument has illegal value."
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            else if (info == 1) then
               write(msg, *) "The reordering of T failed because some eigenvalues are too", &
                           & "close to separate (the problem is very ill-conditioned); ", &
                           & "T may have been partially reordered, and WR and WI ", &
                           & "contain the eigenvalues in the same order as in T; S and", &
                           & "SEP (if requested) are set to zero."
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            else
               write(msg, *) "Undocumented error."
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
               write(msg, *) "Undocumented error."
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
               write(msg, *) "Undocumented error."
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
               write(msg, *) "Undocumented error."
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'qr') then
            ! qr
            if (info > 0) then
               write(msg, *) 'QR factorization: Colinear column detected in column', info
               call logger%log_information(trim(msg), module=module, procedure=procedure)
            else
               write(msg, *) "Undocumented error."
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'qr_pivot') then
            ! qr_pivot
            if (info > 0) then
               write(msg, *) 'QR factorization: Invariant subspace found after', info, 'steps.'
               call logger%log_information(trim(msg), module=module, procedure=procedure)
            else
               write(msg, *) "Undocumented error."
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'arnoldi') then
            ! arnoldi
            if (info > 0) then
               write(msg, *) 'Arnoldi factorization: Invariant subspace computed after', info, 'iterations.'
               call logger%log_information(trim(msg), module=module, procedure=procedure)
            else
               write(msg, *) "Undocumented error."
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'lanczos_bidiagonalization') then
            ! lanczos_bidiagonalization
            if (info > 0) then
               write(msg, *) 'Lanczos Bidiagonalisation: Invariant subspace found after', info, 'steps.'
               call logger%log_information(trim(msg), module=module, procedure=procedure)
            else
               write(msg, *) "Undocumented error."
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'lanczos_tridiagonalization') then
            ! lanczos_tridiagonalization
            if (info > 0) then
               write(msg, *) 'Lanczos Tridiagonalisation: Invariant subspace found after', info, 'steps.'
               call logger%log_information(trim(msg), module=module, procedure=procedure)
            else
               write(msg, *) "Undocumented error."
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
               write(msg, *) "Undocumented error."
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'eighs') then
            ! GMRES
            if (info > 0) then
               write(msg, *) 'eigs iteration converged after', info, 'iterations'
               call logger%log_information(trim(msg), module=module, procedure=procedure)
            else
               write(msg, *) "Undocumented error."
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'svds') then
            ! GMRES
            if (info > 0) then
               write(msg, *) 'svds iteration converged after', info, 'iterations'
               call logger%log_information(trim(msg), module=module, procedure=procedure)
            else
               write(msg, *) "Undocumented error."
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'gmres') then
            ! GMRES
            if (info > 0) then
               write(msg, *) 'GMRES iteration converged after', info, 'iterations'
               call logger%log_information(trim(msg), module=module, procedure=procedure)
            else
               write(msg, *) "Undocumented error."
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'cg') then
            ! CG
            if (info > 0) then
               write(msg, *) 'CG iteration converged after', info, 'iterations'
               call logger%log_information(trim(msg), module=module, procedure=procedure)
            else
               write(msg, *) "Undocumented error."
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
            write(msg, *) "Undocumented error."
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         !
         !   Default
         !
         else 
            write(msg,*) 'subroutine "'//trim(origin)//'" returned with flag: ', info
            call logger%log_error(trim(msg), module=module, procedure=procedure)
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

end module LightKrylov_Logger
