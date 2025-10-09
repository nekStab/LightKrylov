module TestUtils
   use stdlib_optval, only: optval
   use stdlib_ascii, only: to_lower
   use stdlib_strings, only: replace_all
   use testdrive, only: error_type
   implicit none(type, external)
   private
   interface
      module subroutine check_test(error, test_name, info, eq, context)
         use face
         implicit none(type, external)
         type(error_type), allocatable, intent(inout) :: error
         character(len=*), intent(in)    :: test_name
         character(len=*), optional, intent(in)    :: info
         character(len=*), optional, intent(in)    :: eq
         character(len=*), optional, intent(in)    :: context
      end subroutine check_test
   end interface
   public :: check_test
contains
   module subroutine check_test(error, test_name, info, eq, context)
      use face
      implicit none(type, external)
      type(error_type), allocatable, intent(inout) :: error
      character(len=*), intent(in)    :: test_name
      character(len=*), optional, intent(in)    :: info
      character(len=*), optional, intent(in)    :: eq
      character(len=*), optional, intent(in)    :: context
      character(len=128)                           :: name

      ! internals
      character(len=128)                           :: msg, info_, eq_
      ! character(len=*), parameter :: indent = repeat(" ", 7)
      character(len=4), dimension(4) :: substrings
      integer :: i

      info_ = optval(info, '')
      eq_ = optval(eq, '')

      name = trim(to_lower(test_name))
      substrings = ["_rsp", "_rdp", "_csp", "_cdp"]
      do i = 1, size(substrings)
         name = replace_all(name, substrings(i), "")
      end do
      name = replace_all(name, "test_", "")

      write (*, '(A33)', ADVANCE='NO') name
      write (*, '(A3)', ADVANCE='NO') ' % '

      if (len(trim(info_)) == 0) then
         msg = eq_
      else
         if (len(info_) > 30) then
            msg = info_(:30)//eq_
         else
            msg = info_//repeat(' ', 30 - len(trim(info_)))//eq_
         end if
      end if
      write (*, '(A62)', ADVANCE='NO') msg

      if (allocated(error)) then
         print *, colorize('FAILED', color_fg='red')
         if (present(context)) then
            write (*, '(A)', ADVANCE='NO') trim(context)
         end if
         write (*, *)
         write (*, *) 'The most recent test failed. Aborting calculation as per user directive.'
         write (*, *)
         STOP 1
      else
         print *, colorize('PASSED', color_fg='green')
      end if

   end subroutine check_test
end module TestUtils
