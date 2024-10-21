module LightKrylov_Constants
   !! This module defines a list of simple constants used throughout `LightKrylov`.
   !! It also provides some utility functions related to let `LightKrylov` be aware
   !! of any MPI-related information (e.g. rank of the current MPI process, dimension
   !! of the MPI communicator, etc).
   implicit none
   private

   integer, private :: nid = 0
   !! Rank of the current process (local MPI variable).
   integer, private :: comm_size = 1
   !! Dimension of the MPI communicator.
   integer, private :: nio = 0
   !! Rank of the processor for logging information.

   integer , parameter, public :: sp = selected_real_kind(6, 37)
   !! Single precision data type.
   real(sp), parameter, public :: atol_sp = 10.0_sp ** -precision(1.0_sp)
   !! Absolute tolerance for single precision computations.
   real(sp), parameter, public :: rtol_sp = sqrt(atol_sp)
   !! Relative tolerance for single precision computations.
   real(sp), parameter, public :: one_rsp  = 1.0_sp
   !! Real-valued single precision one. 
   real(sp), parameter, public :: zero_rsp = 0.0_sp
   !! Real-valued single precision zero.
   complex(sp), parameter, public :: one_csp    = cmplx(1.0_sp, 0.0_sp, kind=sp)
   !! Complex-valued single precision one.
   complex(sp), parameter, public :: one_im_csp = cmplx(0.0_sp, 1.0_sp, kind=sp)
   !! Complex-valued single precision imaginary unit.
   complex(sp), parameter, public :: zero_csp   = cmplx(0.0_sp, 0.0_sp, kind=sp)
   !! Complex-value single precision zero.

   integer , parameter, public :: dp = selected_real_kind(15, 307)
   !! Double precision data type.
   real(dp), parameter, public :: atol_dp = 10.0_dp ** -precision(1.0_dp)
   !! Absolute tolerance for double precision computations.
   real(dp), parameter, public :: rtol_dp = sqrt(atol_dp)
   !! Relative tolerance for double precision computations.
   real(dp), parameter, public :: one_rdp  = 1.0_dp
   !! Real-valued double precision one.
   real(dp), parameter, public :: zero_rdp = 0.0_dp
   !! Real-valued double precision zero.
   complex(dp), parameter, public :: one_cdp    = cmplx(1.0_dp, 0.0_dp, kind=dp)
   !! Complex-valued double precision one.
   complex(sp), parameter, public :: one_im_cdp = cmplx(0.0_dp, 1.0_dp, kind=dp)
   !! Complex-valued double precision imaginary unit.
   complex(dp), parameter, public :: zero_cdp   = cmplx(0.0_dp, 0.0_dp, kind=dp)
   !! Complex-valued double precision zero.

   ! Getter/setter routines
   public :: set_comm_size
   public :: set_rank
   public :: set_io_rank
   public :: get_rank
   public :: get_comm_size
   public :: io_rank
    
contains

   subroutine set_rank(rank)
      !! Utility function to set the rank of an MPI process.
      integer :: rank
      !! Desired rank identification.
      nid = rank
   end subroutine set_rank

   subroutine set_comm_size(c_size)
      !! Utility function to inform `LightKrylov` of the MPI-communicator's dimension.
      integer :: c_size
      !! Dimension of the MPI communicator.
      comm_size = c_size
   end subroutine set_comm_size

   subroutine set_io_rank(rk)
      !! Utility function to set the rank of the process doing I/O.
      integer, intent(in) :: rk
      !! Desired rank for the IO process.
      if (rk > comm_size .or. rk < 0) then
         if (io_rank()) print *, 'Invalid I/O rank specified!', rk
      else
         nio = rk
         if (io_rank()) print *, 'I/O rank --> rank ', nio
      end if
   end

   pure integer function get_rank() result(rank)
      !! Utility function to get the rank of the current MPI process.
      rank = nid
   end function get_rank
   
   pure integer function get_comm_size() result(c_size)
      !! Utility function to get the dimension of the communicator known to `LightKrylov`.
      c_size = comm_size
   end function get_comm_size

   pure logical function io_rank() result(is_io)
      !! Utility function to determine whether the current MPI process can do I/O.
      is_io = .false.
      if (nid == nio) is_io = .true.
   end function io_rank

end module LightKrylov_Constants
