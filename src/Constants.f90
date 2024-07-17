module LightKrylov_Constants
    use stdlib_logger, only: logger => global_logger
    use LightKrylov_Logger
#ifdef MPI
    use mpi_f08
#endif
    implicit none
    private

    integer, private :: nid = 0
    integer, private :: comm_size = 1
    integer, private :: nio = 0

    integer , parameter, public :: sp = selected_real_kind(6, 37)
    !! Definition of the single precision data type.
    real(sp), parameter, public :: atol_sp = 10.0_sp ** -precision(1.0_sp)
    !! Definition of the absolute tolerance for single precision computations.
    real(sp), parameter, public :: rtol_sp = sqrt(atol_sp)
    !! Definition of the relative tolerance for single precision computations.

    integer , parameter, public :: dp = selected_real_kind(15, 307)
    !! Definition of the double precision data type.
    real(dp), parameter, public :: atol_dp = 10.0_dp ** -precision(1.0_dp)
    !! Definition of the absolute tolerance for double precision computations.
    real(dp), parameter, public :: rtol_dp = sqrt(atol_dp)
    !! Definition of the relative tolerance for double precision computations.

    real(sp), parameter, public :: one_rsp  = 1.0_sp
    real(sp), parameter, public :: zero_rsp = 0.0_sp
    real(dp), parameter, public :: one_rdp  = 1.0_dp
    real(dp), parameter, public :: zero_rdp = 0.0_dp
    complex(sp), parameter, public :: one_csp    = cmplx(1.0_sp, 0.0_sp, kind=sp)
    complex(sp), parameter, public :: one_im_csp = cmplx(0.0_sp, 1.0_sp, kind=sp)
    complex(sp), parameter, public :: zero_csp   = cmplx(0.0_sp, 0.0_sp, kind=sp)
    complex(dp), parameter, public :: one_cdp    = cmplx(1.0_dp, 0.0_dp, kind=dp)
    complex(sp), parameter, public :: one_im_cdp = cmplx(0.0_dp, 1.0_dp, kind=dp)
    complex(dp), parameter, public :: zero_cdp   = cmplx(0.0_dp, 0.0_dp, kind=dp)

    ! MPI subroutines
    public :: comm_setup, comm_close
    
    ! Getter/setter
    public :: set_io_rank
    public :: io_rank
    public :: get_rank
    
contains

   subroutine comm_setup()
      integer :: ierr
      character(len=128) :: msg
#ifdef MPI
      ! Initialize MPI
      call MPI_Init(ierr)
      if (ierr /= MPI_SUCCESS) call stop_error("Error initializing MPI", module='LightKrylov',procedure='mpi_init')
      call MPI_Comm_rank(MPI_COMM_WORLD, nid, ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, comm_size, ierr)
      write(msg, '(A,I4,A,I4)') 'Setup parallel run: rank', nid, ', comm_size = ', comm_size
      call logger%log_message(trim(msg), module='LightKrylov', procedure='mpi_setup')
#else
      write(msg, *) 'Setup serial run'
      call logger%log_message(trim(msg), module='LightKrylov', procedure='mpi_setup')
#endif
      call set_io_rank(0)
      return
   end subroutine comm_setup

   subroutine comm_close
      integer :: ierr
#ifdef MPI
      character(len=128) :: msg
      ! Finalize MPI
      call MPI_Finalize(ierr)
      if (ierr /= MPI_SUCCESS) call stop_error("Error finalizing MPI", module='LightKrylov',procedure='mpi_finalize')
#else
      ierr = 0
#endif
      return
   end subroutine comm_close

   subroutine set_io_rank(rk)
      integer, intent(in) :: rk
      character(len=128) :: msg
      if (rk > comm_size) then
         write(msg, *) 'Invalid I/O rank specified!'
         if (io_rank()) call logger%log_message(trim(msg), module='LightKrylov', procedure='set_io_rank')
      else
         nio = rk
         write(msg, '(A,I4)') 'I/O rank --> rank ', nio
         if (io_rank()) call logger%log_message(trim(msg), module='LightKrylov', procedure='set_io_rank')
      end if
   end 

   logical function io_rank() result(is_io)
      is_io = .false.      
      if (nid == nio) is_io = .true.
   end function io_rank

   integer function get_rank() result(rank)
      rank = nid
   end function get_rank

end module LightKrylov_Constants
