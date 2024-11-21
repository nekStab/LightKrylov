module LightKrylov_Logger
#ifdef MPI
   use mpi_f08
#endif
   ! Fortran Standard Library
   use stdlib_optval, only : optval
   use stdlib_logger, only: logger => global_logger
   use stdlib_ascii, only : to_lower
   use stdlib_strings, only : chomp, replace_all
   ! Testdrive
   use testdrive, only: error_type, test_failed
   ! LightKrylov
   use LightKrylov_Constants

   implicit none
   private

   character(len=128), parameter, private :: this_module = 'LightKrylov_Logger'

   logical, parameter, private :: exit_on_error = .true.
   logical, parameter, private :: exit_on_test_error = .true.

   public :: stop_error
   public :: check_info
   public :: check_test
   public :: logger

   public :: logger_setup

   ! MPI subroutines
   public :: comm_setup
   public :: comm_close

contains

   subroutine logger_setup(logfile, nio, log_level, log_stdout, log_timestamp, iunit)
      !! Wrapper to set up MPI if needed and initialize log files
      character(len=*), optional, intent(in) :: logfile
      !! name of the dedicated LightKrylov logfile
      integer, optional, intent(in)          :: nio
      !! I/O rank for logging
      integer, optional, intent(in)          :: log_level
      !! set logging level
      !! 0   : all_level
      !! 10  : debug_level
      !! 20  : information_level
      !! 30  : warning_level	
      !! 40  : error_level	
      !! 100 : none_level
      logical, optional, intent(in)          :: log_stdout
      !! duplicate log messages to stdout?
      logical, optional, intent(in)          :: log_timestamp
      !! add timestamp to log messages
      integer, optional, intent(out)         :: iunit
      !! log unit identifier

      ! internals
      character(len=:), allocatable :: logfile_
      integer                       :: nio_
      integer                       :: log_level_
      logical                       :: log_stdout_
      logical                       :: log_timestamp_
      integer                       :: iunit_
      ! misc
      integer :: stat

      logfile_       = optval(logfile, 'lightkrylov.log')
      nio_           = optval(nio, 0)
      log_level_     = optval(log_level, 20)
      log_level_     = max(0, min(log_level_, 100))
      log_stdout_    = optval(log_stdout, .true.)
      log_timestamp_ = optval(log_timestamp, .true.)

      ! set log level
      call logger%configure(level=log_level_, time_stamp=log_timestamp_) 

      ! set up LightKrylov log file
      call logger%add_log_file(logfile_, unit=iunit, stat=stat)
      if (stat /= 0) call stop_error('Unable to open logfile '//trim(logfile_)//'.', module=this_module, procedure='logger_setup')

      ! Set up comms
      call comm_setup()

      ! Set I/O rank
      if (nio_/=0) call set_io_rank(nio_)

      ! log to stdout
      if (log_stdout_) then
         call logger%add_log_unit(6, stat=stat)
         if (stat /= 0) call stop_error('Unable to add stdout to logger.', module=this_module, procedure='logger_setup')
      end if

      ! return unit if requested
      if (present(iunit)) iunit = iunit_
      return
   end subroutine logger_setup

   subroutine comm_setup()
      ! internal
      integer :: ierr, nid, comm_size
      logical :: mpi_is_initialized
      character(len=128) :: msg
#ifdef MPI
      ! check if MPI has already been initialized and if not, initialize
      call MPI_Initialized(mpi_is_initialized, ierr)
      if (.not. mpi_is_initialized) then
         call logger%log_debug('Set up parallel run with MPI.', module='LightKrylov', procedure='comm_setup')
         call MPI_Init(ierr)
         if (ierr /= MPI_SUCCESS) call stop_error("Error initializing MPI", module='LightKrylov',procedure='mpi_init')
      else
         call logger%log_debug('MPI already initialized.', module='LightKrylov', procedure='comm_setup')
      end if
      call MPI_Comm_rank(MPI_COMM_WORLD, nid, ierr); call set_rank(nid)
      call MPI_Comm_size(MPI_COMM_WORLD, comm_size, ierr); call set_comm_size(comm_size)
      write(msg,'(A,I0,A,I0)') 'rank', nid, ', comm_size = ', comm_size
      call logger%log_debug(trim(msg), module='LightKrylov', procedure='comm_setup')
#else
      write(msg,'(A)') 'Setup serial run'
      call set_rank(0)
      call set_comm_size(1)
      call logger%log_debug(trim(msg), module='LightKrylov', procedure='comm_setup')
#endif
      return
   end subroutine comm_setup

   subroutine comm_close()
      integer :: ierr
#ifdef MPI
      character(len=128) :: msg
      ! Finalize MPI
      call MPI_Finalize(ierr)
      if (ierr /= MPI_SUCCESS) call stop_error("Error finalizing MPI", module='LightKrylov',procedure='comm_close')
#else
      ierr = 0
#endif
      return
   end subroutine comm_close

   subroutine stop_error(msg, module, procedure)
      character(len=*),            intent(in)  :: msg
      !! The name of the procedure in which the call happens
      character(len=*), optional,  intent(in)  :: module
      !! The name of the module in which the call happens
      character(len=*), optional,  intent(in)  :: procedure
      !! The name of the procedure in which the call happens
      call check_info(-1, origin="STOP_ERROR", module=module, procedure=procedure, info_msg=msg)
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
         write(msg,'(A)') 'The subroutine "'//trim(origin)//'" returned successfully. '//trim(str)
         call logger%log_debug(trim(msg), module=module, procedure=procedure)
      else
         !
         !   LAPACK
         !
         if (trim(to_lower(origin)) == 'getref') then
            ! GETREF
            if (info < 0 ) then
               write(msg,'(A,I0,A)') "The ", -info, "-th argument has illegal value. "//trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            else
               write(msg,'(A,I0,A,I0,A)') "U(", info, ",", info, ") is exactly zero. The factorization ", &
                           & "has been completed but the factor U is exactly singular. ", &
                           & "Division by zero will occur if used to solve Ax=b. "//trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'getri') then
            ! GETRI
            if (info < 0 ) then
               write(msg,'(A,I0,A)') "The ", -info, "-th argument has illegal value. "//trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            else
               write(msg,'(A,I0,A)') "U(", info, ",", info, ") is exactly zero. ", &
                           & "The matrix is singular and its inverse cannot be computed. "//trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'geev') then
            ! GEEV
            if (info < 0) then
               write(msg,'(A,I0,A)') "The ", -info, "-th argument has illegal value. "//trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            else 
               write(msg,'(A,I0,A)') "The QR alg. failed to compute all of the eigenvalues.", &
                           & "No eigenvector has been computed. "//trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'syev') then
            ! SYEV
            if (info < 0) then
               write(msg,'(A,I0,A)') "The ", -info, "-th argument has illegal value. "//trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            else 
               write(msg,'(A)') "The QR alg. failed to compute all of the eigenvalues.", &
                           & "No eigenvector has been computed. "//trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'heev') then
            ! HEEV
            if (info < 0) then
               write(msg,'(A,I0,A)') "The ", -info, "-th argument has illegal value. "//trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            else 
               write(msg,'(A)') "The QR alg. failed to compute all of the eigenvalues.", &
                           & "No eigenvector has been computed. "//trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'gels') then
            ! GELS
            if (info < 0) then
               write(msg,'(A,I0,A)') "The ", -info, "-th argument has illegal value. "//trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            else
               write(msg,'(A)') "Undocumented error. "//trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'gees') then
            ! GEES
            if (info < 0) then
               write(msg,'(A,I0,A)') "The ", -info, "-th argument has illegal value. "//trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            else
               write(msg,'(A)') "The QR alg. failed to compute all of the eigenvalues.", &
                           & "No eigenvector has been computed. "//trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'trsen') then
            ! GEES
            if (info < 0) then
               write(msg,'(A,I0,A)') "The ", -info, "-th argument has illegal value. "//trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            else if (info == 1) then
               write(msg,'(A)') "The reordering of T failed because some eigenvalues are too", &
                           & "close to separate (the problem is very ill-conditioned); ", &
                           & "T may have been partially reordered, and WR and WI ", &
                           & "contain the eigenvalues in the same order as in T; S and", &
                           & "SEP (if requested) are set to zero. "//trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            else
               write(msg,'(A)') "Undocumented error. "//trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         !
         !   LightKrylov_Utils
         !
         else if (trim(to_lower(origin)) == 'sqrtm') then
            if (info == 1) then
               write(msg,'(A)') 'The input matrix is singular to tolerance. The singular eigenvalues are set to zero. '
               call logger%log_warning(trim(msg), module=module, procedure=procedure)
            else if (info == -1) then
               write(msg,'(A)') "The input matrix is not positive (semi-)definite. "
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            else
               write(msg,'(A)') "Undocumented error. "//trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         !
         !   LightKrylov_BaseKrylov
         !
         else if (trim(to_lower(origin)) == 'orthogonalize_against_basis') then ! the regular case
            ! orthogonalization
            if (info > 0) then
               write(msg,'(A,I0,A)') 'Orthogonalization: The ', info, 'th input vector is numerically zero.'
               call logger%log_debug(trim(msg), module=module, procedure=procedure)
            else if (info == -1) then
               write(msg,'(A)') 'The input Krylov basis is not orthonormal.'
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            else if (info == -2) then
               write(msg,'(A)') 'Orthogonalization: The last column of the input basis is zero.'
               call logger%log_debug(trim(msg), module=module, procedure=procedure)
            else
               write(msg,'(A)') "Undocumented error. "//trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'orthogonalize_against_basis_p1') then
            ! orthogonalization
            if (info > 0) then
               write(msg,'(A,I0,A)') 'Orthogonalization: The ', info, 'th input vector is numerically zero.'
               call logger%log_debug(trim(msg), module=module, procedure=procedure)
            else if (info == -1) then
               write(msg,'(A)') 'The input Krylov basis is not orthonormal.'
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            else if (info == -2) then
               write(msg,'(A)') 'Orthogonalization: The last column of the input basis is zero.'
               call logger%log_warning(trim(msg), module=module, procedure=procedure)
            else
               write(msg,'(A)') "Undocumented error. "//trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'orthogonalize_against_basis_p2') then
            ! orthogonalization
            if (info > 0) then ! show this information only for debugging
               write(msg,'(A,I0,A)') 'Orthogonalization: The ', info, 'th input vector is numerically zero.'
               call logger%log_debug(trim(msg), module=module, procedure=procedure)
            else if (info == -1) then
               write(msg,'(A)') 'The input Krylov basis is not orthonormal.'
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            else if (info == -2) then
               write(msg,'(A)') 'Orthogonalization: The last column of the input basis is zero.'
               call logger%log_warning(trim(msg), module=module, procedure=procedure)
            else
               write(msg,'(A)') "Undocumented error. "//trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'double_gram_schmidt_step') then
            ! orthogonalization
            if (info > 0) then
               write(msg,'(A,I0,A)') 'Orthogonalization: The ', info, 'th input vector is numerically zero.'
               call logger%log_debug(trim(msg), module=module, procedure=procedure)
            else if (info == -1) then
               write(msg,'(A)') 'The input Krylov basis is not orthonormal.'
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            else if (info == -2) then
               write(msg,'(A)') 'Orthogonalization: The last column of the input basis is zero.'
               call logger%log_warning(trim(msg), module=module, procedure=procedure)
            else
               write(msg,'(A)') "Undocumented error. "//trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'qr') then
            ! qr
            if (info > 0) then
               write(msg,'(A,I0,A)') 'QR factorization: Colinear column detected in column ', info,  &
                           & '. NOTE: Other subsequent columns may also be colinear.'
               call logger%log_debug(trim(msg), module=module, procedure=procedure)
            else
               write(msg,'(A)') "Undocumented error. "//trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'qr_pivot') then
            ! qr_pivot
            if (info > 0) then
               write(msg,'(A,I0,A)') 'QR factorization: Invariant subspace found after ', info, ' steps.'
               call logger%log_debug(trim(msg), module=module, procedure=procedure)
            else
               write(msg,'(A)') "Undocumented error. "//trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'arnoldi') then
            ! arnoldi
            if (info > 0) then
               write(msg,'(A,I0,A)') 'Arnoldi factorization: Invariant subspace computed after ', info, ' iterations.'
               call logger%log_debug(trim(msg), module=module, procedure=procedure)
            else
               write(msg,'(A)') "Undocumented error. "//trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'bidiagonalization') then
            ! lanczos_bidiagonalization
            if (info > 0) then
               write(msg,'(A,I0,A)') 'Lanczos Bidiagonalisation: Invariant subspace found after ', info, ' steps.'
               call logger%log_debug(trim(msg), module=module, procedure=procedure)
            else
               write(msg,'(A)') "Undocumented error. "//trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'lanczos') then
            ! lanczos_tridiagonalization
            if (info > 0) then
               write(msg,'(A,I0,A)') 'Lanczos Tridiagonalisation: Invariant subspace found after ', info, ' steps.'
               call logger%log_debug(trim(msg), module=module, procedure=procedure)
            else
               write(msg,'(A)') "Undocumented error. "//trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         !
         !   LightKrylov_IterativeSolvers
         !
         else if (trim(to_lower(origin)) == 'eigs') then
            ! GMRES
            if (info > 0) then
               write(msg,'(A,I0,A)') 'eigs iteration converged after ', info, ' iterations'
               call logger%log_information(trim(msg), module=module, procedure=procedure)
            else
               write(msg,'(A)') "Undocumented error. "//trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'eighs') then
            ! GMRES
            if (info > 0) then
               write(msg,'(A,I0,A)') 'eigs iteration converged after ', info, ' iterations'
               call logger%log_information(trim(msg), module=module, procedure=procedure)
            else
               write(msg,'(A)') "Undocumented error. "//trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'svds') then
            ! GMRES
            if (info > 0) then
               write(msg,'(A,I0,A)') 'svds iteration converged after ', info, ' iterations'
               call logger%log_information(trim(msg), module=module, procedure=procedure)
            else
               write(msg,'(A)') "Undocumented error. "//trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'gmres') then
            ! GMRES
            if (info > 0) then
               write(msg,'(A,I0,A)') 'GMRES iteration converged after ', info, ' iterations'
               call logger%log_message(trim(msg), module=module, procedure=procedure)
            else
               write(msg,'(A)') "Undocumented error. "//trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'cg') then
            ! CG
            if (info > 0) then
               write(msg,'(A,I0,A)') 'CG iteration converged after ', info, ' iterations'
               call logger%log_message(trim(msg), module=module, procedure=procedure)
            else
               write(msg,'(A)') "Undocumented error. "//trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         else if (trim(to_lower(origin)) == 'linear_solver') then
            ! Abstract linear solver
            if (info > 0) then
               write(msg,'(A,I0,A)') 'The linear solver converged after ', info, ' iterations'
               call logger%log_message(trim(msg), module=module, procedure=procedure)
            else
               write(msg,'(A)') "Undocumented error. "//trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         !
         !   LightKrylov_ExpmLib
         !
         else if (trim(to_lower(origin)) == 'kexpm') then
            ! Krylov Matrix Exponential
            if (info > 0) then
               write(msg,'(A,I0,A)') 'kexpm converged. Estimated error below tolerance using ', info, ' Krylov vectors.'
               call logger%log_debug(trim(msg), module=module, procedure=procedure)
            else if (info == -2) then
               write(msg,'(A)') 'kexpm converged. Arnoldi iteration breakdown. Approximation is exact to arnoldi tolerance.'
               call logger%log_debug(trim(msg), module=module, procedure=procedure)
            else if (info == -1) then
               write(msg,'(A)') 'kexpm did not converge. Maximum number of Krylov vectors reached.'
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            write(msg,'(A)') "Undocumented error. "//trim(str)
               call logger%log_error(origin, module=module, procedure=procedure, stat=info, errmsg=trim(msg))
               ierr = -1
            end if
         !
         !  stop error
         !
         else if (trim(origin) == 'STOP_ERROR') then
            call logger%log_error(trim(origin), module=module, procedure=procedure, stat=info, errmsg=trim(str))
            ierr = -1
         !
         !   Default
         !
         else
            write(msg,'(A)') 'subroutine "'//trim(origin)//'" returned with a non-zero error flag.'
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
