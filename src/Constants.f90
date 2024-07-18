module LightKrylov_Constants
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

   ! Getter/setter routines
   public :: get_rank
   public :: set_io_rank
   public :: io_rank
    
contains

   subroutine set_io_rank(rk)
      integer, intent(in) :: rk
      if (rk > comm_size .or. rk < 0) then
         if (io_rank()) print *, 'Invalid I/O rank specified!', rk
      else
         nio = rk
         if (io_rank()) print *, 'I/O rank --> rank ', nio
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
