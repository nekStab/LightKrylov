module LightKrylov_Timer_Utils
   use stdlib_optval, only: optval
   use stdlib_ascii, only: to_lower
   use LightKrylov_Constants, only: dp
   use LightKrylov_Logger
   implicit none
   private
   character(len=*), parameter :: this_module      = 'LK_TmrUtils'
   character(len=*), parameter :: this_module_long = 'LightKrylov_Timer_Utils'
   
   ! Timer type
   type, public :: lightkrylov_timer
      !! Individual timer.
      !! Atomic timer that is associated to a particular 'event' by name which may be a procedure or a 
      !! user-defined string at instantiation.
      !!
      !! The timing information in gathered for each timer independently. The individual timers are gathered
      !! into groups (relevant only for timing output) and managed by a central watch that is derived from 
      !! the `abstract_watch` type. The individual timers are rarely user independently but all timing 
      !! actions are typically performed via procedures in the central timer.
      !!
      !! A notable exception are the `matvec`/`rmatvec` as well as `reponse` timers associated with the types
      !! `abstract_linop` and `abstract_system`, respectively, which are managed via their parent types and the
      !! corresponding type-bound procedures only.
      private
      character(len=128), public :: name = 'default_timer'
      !! Timer name.
      real(dp) :: etime       = 0.0_dp
      !! Elapsed time since reset.
      real(dp) :: etime_pause = 0.0_dp
      !! Elapsed time up until most recent pause.
      real(dp) :: start_time  = 0.0_dp
      !! Start time for comparison.
      real(dp) :: etime_max   = 0.0_dp
      !! Maximum elapsed time since reset.
      real(dp) :: etime_min   = huge(1.0_dp)
      !! Minimum elapsed time since reset.
      integer :: local_count  = 0
      !! Call counter since reset.
      integer :: reset_count  = 0
      !! Reset counter.
      integer :: count        = 0
      !! Global counter (only reset when data is flushed).
      logical :: running      = .false.
      !! Protection against repeated starts.
      logical :: is_finalized = .false.
      !! Switch for printing.
      real(dp), dimension(:), allocatable :: etime_data
      real(dp), dimension(:), allocatable :: etavg_data
      real(dp), dimension(:), allocatable :: etmin_data
      real(dp), dimension(:), allocatable :: etmax_data
      integer,  dimension(:), allocatable :: count_data
   contains
      private
      procedure, pass(self), public :: start => start_timer
      procedure, pass(self), public :: stop => stop_timer
      procedure, pass(self), public :: pause => pause_timer
      procedure, pass(self), public :: reset => reset_timer
      !! Reset timing data (soft/hard, clear/save data).
      procedure, pass(self), public :: finalize => finalize_timer
      !! Aggregate data and prepare summary.
      procedure, pass(self), public :: get_time => get_timer_etime
      !! Getter routine to access self%etime only.
      procedure, pass(self), public :: get_data => get_timer_data
      !! Getter routine to access to all local data: etime, counts, etime_max, etime_min, etime_pause.
      procedure, pass(self), public :: print_info => print_timer_info
      !! Print timing data.
      procedure, pass(self), public :: save_timer_data
      !! Transfer timing data to arrays.
   end type lightkrylov_timer

   ! Timer group type
   type, public :: lightkrylov_timer_group
      !! Simple type to allow for some structure in the timer output.
      private
      character(len=128), public :: name = 'default_group'
      !! group name.
      integer :: istart = 0
      !! global index of first timer in group.
      integer :: iend = 0
      !! global index of last timer in group.
   end type lightkrylov_timer_group

   ! Abstract watch type
   type, abstract, public :: abstract_watch
      !! Base type to define a global timer.
      !! All watches within LightKrylov and associated tools are derived from this type.
      !!
      !! Within LightKrylov, the derived type `global_lightkrylov_timer` is used to manage all
      !! atomic timers associated with both internal routines (private) as well as user-defined
      !! (public) timers that can be added and removed as necessary. In order to protect the 
      !! private timers, they are defined and set only during the initialisation of the derived 
      !! type via the deferred procedure `set_private_timers_and_name` if timing is requested 
      !! and cannot be deleted.
      !! Once the global watch is initalized, the user can define and subsequently remove 
      !! user-defined timers at any point that will be managed by the global watch in the same
      !! way as private timers.
      !! Upon finalization, the user-defined timers with the associated timing information is 
      !! presented together with that of the private timers.
      !!
      !! The type-bound procedures of the abstract_watch type allow individual access to each
      !! managed timer individually using the timer name as reference.
      private
      character(len=128) :: name = 'default_watch'
      type(lightkrylov_timer), dimension(:), allocatable :: timers
      !! Array of timers contained in the watch.
      integer :: timer_count    = 0
      !! Number of timers managed by watch.
      type(lightkrylov_timer_group), dimension(:), allocatable :: groups
      !! Array of timer groups contained in the watch.
      integer :: group_count    = 0
      !! Number of timer groups managed by watch.
      integer :: private_count  = 0
      !! Number of private timers (immutable by user).
      integer :: user_count     = 0
      !! Number of user defined timers.
      logical :: user_mode      = .false.
      !! Number of user defined timers.
      logical :: is_initialized = .false.
   contains
      private
      procedure, pass(self), public :: add_timer
      !! Add new timer to the watch.
      procedure, pass(self), public :: remove_timer
      !! Remove existing timer from the watch.
      procedure, pass(self), public :: add_group
      !! Add new timer group to the watch.
      ! Getter/Setter and helper routines for watch components.
      procedure, pass(self), public :: get_timer_id
      procedure, pass(self), public :: get_group_id
      procedure, pass(self), public :: set_watch_name
      ! Wrappers for the basic timing routines/access to timer components.
      procedure, pass(self), public :: start => start_timer_by_name
      procedure, pass(self), public :: stop => stop_timer_by_name
      procedure, pass(self), public :: pause => pause_timer_by_name
      procedure, pass(self), public :: reset => reset_timer_by_name
      procedure, pass(self), public :: get_time => get_timer_etime_by_name
      procedure, pass(self), public :: get_data => get_timer_data_by_name
      procedure, pass(self), public :: print_info => print_timer_info_by_name
      ! Global manipulation routines.
      procedure, pass(self), public :: reset_all
      !! Reset all timers in watch.
      procedure, pass(self), public :: enumerate
      !! Print summary of registered timers and their current status.
      procedure, pass(self), public :: initialize
      !! Set up private timers, flags and counters. Switch on timing.
      procedure, pass(self), public :: finalize
      !! Gather timing information and print it to screen/logfile.
      procedure(abstract_set_timers), pass(self), deferred, public :: set_private_timers_and_name
      !! Define private timers that cannot be removed by the user.
   end type abstract_watch

   abstract interface
      subroutine abstract_set_timers(self)
         !! Interface for defining timers and groups.
         import abstract_watch
         class(abstract_watch), intent(inout) :: self
      end subroutine abstract_set_timers
   end interface

   ! format strings for uniform printing
   character(len=128), parameter :: fmt_h = '(2X,A30," : ",  10X,A9,4(A15))'      ! headers
   character(len=128), parameter :: fmt_t = '(2X,A30," : ",A6,4X,I9,4(1X,F14.6))' ! data total
   character(len=128), parameter :: fmt_1 = '(2X,A30," : ",A6,4X,I9,1X,F14.6,3(A15))' ! data 1 call
   character(len=128), parameter :: fmt_r = '(2X,30X,3X,   A6,I4,I9,4(1X,F14.6))' ! data reset
   character(len=128), parameter :: fmt_n = '(2X,30X,3X,   A6,I4,I9,A60)'         ! not called

contains

   !--------------------------------------------------------------
   !  Type-bound procedures for lightkrylov_timer type
   !--------------------------------------------------------------

   subroutine start_timer(self)
      !! Type-bound to lightkrylov_timer: Start timer.
      class(lightkrylov_timer), intent(inout) :: self
      if (.not. self%running) then
         call cpu_time(self%start_time)
         self%running = .true.
         self%count = self%count + 1
         self%local_count = self%local_count + 1
      end if
   end subroutine start_timer

   subroutine stop_timer(self)
      !! Type-bound to lightkrylov_timer: Stop timer.
      class(lightkrylov_timer), intent(inout) :: self
      ! internal
      real(dp) :: t_now, etime
      call cpu_time(t_now)
      if (self%running) then
         etime            = t_now - self%start_time
         self%etime       = self%etime + etime + self%etime_pause
         self%etime_pause = 0.0_dp
         self%etime_min   = min(self%etime_min, etime)
         self%etime_max   = max(self%etime_max, etime)
         self%running     = .false.
      end if
   end subroutine stop_timer

   subroutine pause_timer(self)
      !! Type-bound to lightkrylov_timer: Pause timer.
      class(lightkrylov_timer), intent(inout) :: self
      ! internal
      real(dp) :: t_now
      call cpu_time(t_now)
      if (self%running) then
         self%etime_pause = self%etime_pause + (t_now - self%start_time)
         self%running = .false.
      end if
   end subroutine pause_timer

   subroutine save_timer_data(self)
      !! Type-bound to lightkrylov_timer: Save current timing information. 
      !! Note: This is done irrespective of the call/run status of the timer.
      class(lightkrylov_timer), intent(inout) :: self
      if (self%reset_count == 0) then
         allocate(self%etime_data(1))
         allocate(self%etmin_data(1))
         allocate(self%etmax_data(1))
         allocate(self%etavg_data(1))
         allocate(self%count_data(1))
         if (self%local_count > 0) then
            self%etime_data(1) = self%etime
            self%etmin_data(1) = self%etime_min
            self%etmax_data(1) = self%etime_max
            self%etavg_data(1) = self%etime/self%local_count
            self%count_data(1) = self%local_count
         else
            self%etime_data(1) = 0.0_dp
            self%etavg_data(1) = 0.0_dp
            self%etmin_data(1) = 0.0_dp
            self%etmax_data(1) = 0.0_dp
            self%count_data(1) = 0
         end if
         self%reset_count = 1
      else
         if (self%local_count > 0) then
            self%etime_data = [ self%etime_data, self%etime ]
            self%etmin_data = [ self%etmin_data, self%etime_min ]
            self%etmax_data = [ self%etmax_data, self%etime_max ]
            self%etavg_data = [ self%etavg_data, self%etime/self%local_count ]
            self%count_data = [ self%count_data, self%local_count ]
         else
            self%etime_data = [ self%etime_data, 0.0_dp ]
            self%etmin_data = [ self%etmin_data, 0.0_dp ]
            self%etmax_data = [ self%etmax_data, 0.0_dp ]
            self%etavg_data = [ self%etavg_data, 0.0_dp ]
            self%count_data = [ self%count_data, 0 ]
         end if
         self%reset_count = self%reset_count + 1
      end if
   end subroutine save_timer_data

   subroutine reset_timer(self, soft, clean, verbose)
      !! Type-bound to lightkrylov_timer: Reset timer.
      class(lightkrylov_timer), intent(inout) :: self
      logical, optional, intent(in) :: soft
      !! Save timing data and reset only if data was collected (i.e. timer was called), default = .true.
      logical, optional, intent(in) :: clean
      !! Flush timing data as well as previously saved timing data, default = .false.
      logical, optional, intent(in) :: verbose
      !! Always print information about the reset process.
      ! internal
      logical :: save_data, flush_timer, print_info
      character(len=128) :: msg
      save_data = optval(soft, .true.)
      flush_timer  = optval(clean, .false.)
      print_info   = optval(verbose, .false.)
      if (self%running) then
         write(msg,'(A)') 'Timer "'//trim(self%name)//'" is curently running. Timer not reset.'
         call logger%log_message(msg, this_module, 'reset_timer')
      else
         write(msg,'(A,L1,3X,A,L1)') 'soft reset: ', save_data, 'flush timers: ', flush_timer
         if (print_info) then
            call logger%log_message(msg, this_module, self%name)
         else
            call logger%log_debug(msg, this_module, self%name)
         end if
         if (save_data .and. .not. flush_timer) then
            if (self%local_count > 0) then
               call self%save_timer_data()
               self%etime       = 0.0_dp
               self%etime_pause = 0.0_dp
               self%start_time  = 0.0_dp
               self%etime_min   = huge(1.0_dp)
               self%etime_max   = 0.0_dp
               self%running     = .false.
               self%local_count = 0
            end if
         else
            ! hard reset
            self%etime       = 0.0_dp
            self%etime_pause = 0.0_dp
            self%etime_min   = huge(1.0_dp)
            self%etime_max   = 0.0_dp
            self%start_time  = 0.0_dp
            self%running     = .false.
            self%local_count = 0
            self%reset_count = 0
            if(allocated(self%etime_data)) deallocate(self%etime_data)
            if(allocated(self%etmin_data)) deallocate(self%etmin_data)
            if(allocated(self%etmax_data)) deallocate(self%etmax_data)
            if(allocated(self%etavg_data)) deallocate(self%etavg_data)
            if(allocated(self%count_data)) deallocate(self%count_data)
         end if
         if (flush_timer) then
            self%count = 0
            self%is_finalized = .false.
         end if
      end if
   end subroutine reset_timer

   real(dp) function get_timer_etime(self, restart) result(etime)
      !! Type-bound to lightkrylov_timer: Getter routine to return the current timer etime.
      !! Note: If it is running, the timer is stopped.
      class(lightkrylov_timer) :: self
      logical, optional :: restart
      ! internal
      logical :: restart_timer
      restart_timer = optval(restart, .false.)
      if (self%running) call self%stop()
      etime = self%etime
      if (restart_timer) call self%start()
   end function

   subroutine get_timer_data(self, restart, etime, etmin, etmax, etimp, lcount, rcount, gcount)
      !! Type-bound to lightkrylov_timer: Getter routine to return the timer data.
      !! Note: If it is running, the timer is stopped.
      class(lightkrylov_timer), intent(inout) :: self
      logical, optional, intent(in) :: restart
      real(dp), optional, intent(out) :: etime
      real(dp), optional, intent(out) :: etmin
      real(dp), optional, intent(out) :: etmax
      real(dp), optional, intent(out) :: etimp
      integer, optional, intent(out) :: lcount
      integer, optional, intent(out) :: rcount
      integer, optional, intent(out) :: gcount
      ! internal
      logical :: restart_timer
      restart_timer = optval(restart, .false.)
      if (self%running) call self%stop()
      if (present(etime)) etime = self%etime
      if (present(etmin)) etmin = self%etime_min
      if (present(etmax)) etmax = self%etime_max
      if (present(etimp)) etimp = self%etime_pause
      if (present(lcount)) lcount = self%local_count
      if (present(rcount)) rcount = self%reset_count
      if (present(gcount)) gcount = self%count
      if (restart_timer) call self%start()
   end subroutine get_timer_data

   subroutine print_timer_info(self)
      !! Type-bound to lightkrylov_timer: Compute spimple statistics and print timing information to screen.
      class(lightkrylov_timer), intent(inout) :: self
      ! internal
      character(len=128) :: msg
      call logger%log_message('#########        Timer info        ##########', this_module)
      if (self%count == 0) then
         write(msg, '(A)') 'No timing data available for "'//trim(self%name)//'": Timer not called.'
         call logger%log_message(msg, this_module)
      else
         call print_summary_header('Summary', self%name)
         if (.not.self%is_finalized) then
            call self%save_timer_data()
         else ! is_finalized
            if (self%reset_count == 0) then
               call stop_error(trim(self%name)//': reset_count inconsistent!', this_module, 'finalize_timer')
            end if
         end if
         call print_summary(self)
      end if
   end subroutine print_timer_info

   subroutine finalize_timer(self, if_silent)
      !! Type-bound to lightkrylov_timer: Prepare timer summary.
      class(lightkrylov_timer), intent(inout) :: self
      logical, optional, intent(in) :: if_silent
      !! No output
      ! internal
      logical :: silent
      character(len=128) :: msg
      silent = optval(if_silent, .false.)
      call self%stop()
      call self%save_timer_data()
      self%is_finalized = .true.
      if (.not. silent) then
         write(msg,'(*(A))') trim(self%name), ' finalization complete.'
         call logger%log_message(msg, this_module)
         call self%print_info()
      end if
   end subroutine finalize_timer

   !--------------------------------------------------------------
   !  Type-bound procedures for abstract_watch type
   !--------------------------------------------------------------

   subroutine add_timer(self, name, start, count)
      !! Type-bound to abstract_watch: Add timer to watch and optionally start it immediately.
      !! Note: The new timer name must be unique.
      class(abstract_watch), intent(inout) :: self
      character(len=*), intent(in) :: name
      logical, optional, intent(in) :: start
      integer, optional, intent(out) :: count
      ! internal
      logical :: start_
      character(len=128) :: msg, tname
      start_ = optval(start, .false.)
      tname = to_lower(name)
      if (self%timer_count == 0) then
         allocate(self%timers(1))
         self%timers(1) = lightkrylov_timer(tname)
         self%timer_count = 1
      else
         if (self%get_timer_id(name) > 0) call element_exists(tname, 'Timer', 'add_timer')
         self%timers = [ self%timers, lightkrylov_timer(tname) ]
         self%timer_count = self%timer_count + 1
         if (self%user_mode) self%user_count = self%user_count + 1
      end if
      write(msg,'(A,I0)') 'Timer "'//trim(tname)//'" added: timer_count: ', self%timer_count
      call logger%log_debug(msg, this_module)
      if (present(count)) count = self%timer_count
      if (start_) call self%start(tname)
   end subroutine add_timer

   subroutine remove_timer(self, name, count)
      !! Type-bound to abstract_watch: Remove timer from watch.
      !! Note: Timers considered private (defined during initialisation) cannot be removed.
      class(abstract_watch), intent(inout) :: self
      character(len=*), intent(in) :: name
      integer, optional, intent(out) :: count
      ! internal
      type(lightkrylov_timer), dimension(:), allocatable :: new_timers
      character(len=128) :: msg, tname
      integer :: id
      tname = to_lower(name)
      id = self%get_timer_id(tname)
      if (id == 0) call timer_not_found(tname, 'remove_timer')
      if (id <= self%private_count) then
         write(msg,'(A)') 'Cannot remove private timer "'//trim(tname)//'". Do nothing.'
         call logger%log_message(msg, this_module, 'remove_timer')
      else
         self%timer_count = self%timer_count - 1
         allocate(new_timers(self%timer_count))
         new_timers(1:id-1) = self%timers(1:id-1)
         new_timers(id:)    = self%timers(id+1:)
         deallocate(self%timers)
         self%timers = new_timers
         write(msg,'(A,I0)') 'Timer "'//trim(tname)//'" removed: timer_count: ', self%timer_count
         call logger%log_debug(msg, this_module, 'remove_timer')
      end if
      if (present(count)) count = self%timer_count
   end subroutine remove_timer

   subroutine add_group(self, name, istart, iend, count)
      !! Type-bound to abstract_watch: Add timer group to watch.
      !! Note: The new group name must be unique. This is a quick hack and should be done better.
      class(abstract_watch), intent(inout) :: self
      character(len=*), intent(in) :: name
      integer, intent(in) :: istart
      integer, intent(in) :: iend
      integer, optional, intent(out) :: count
      ! internal
      character(len=*), parameter :: this_procedure = 'add_group'
      character(len=128) :: msg, gname
      ! Sanity checks
      if (istart < 1 .or. iend < 1) then
         call stop_error('Inconsistent input for istart, iend.', this_module, this_procedure)
      else if (istart > iend) then
         call stop_error('istart > iend.', this_module, this_procedure)
      else if (iend > self%timer_count) then
         call stop_error('iend > timer_count.', this_module, this_procedure)
      end if
      gname = to_lower(name)
      if (self%group_count == 0) then
         allocate(self%groups(1))
         self%groups(1) = lightkrylov_timer_group(name=gname, istart=istart, iend=iend)
         self%group_count = 1
      else
         if (self%get_group_id(name) > 0) call element_exists(gname, 'Group', this_module)
         self%groups = [ self%groups, lightkrylov_timer_group(name=gname, istart=istart, iend=iend) ]
         self%group_count = self%group_count + 1
      end if
      write(msg,'(A,I0)') 'Timer group "'//trim(gname)//'" added: group_count: ', self%group_count
      call logger%log_debug(msg, this_module)
      if (present(count)) count = self%group_count
   end subroutine add_group

   integer function get_timer_id(self, name) result(id)
      !! Type-bound to abstract_watch: Getter routine to return the timer id based on name.
      class(abstract_watch) :: self
      character(len=*)  :: name
      !! Timer name
      ! internal
      integer :: i
      id = 0
      do i = 1, self%timer_count
         if (self%timers(i)%name == to_lower(name)) id = i
      end do
   end function get_timer_id

   integer function get_group_id(self, name) result(id)
      !! Type-bound to abstract_watch: Getter routine to return the group id based on name.
      class(abstract_watch) :: self
      character(len=*) :: name
      !! Timer name
      ! internal
      integer :: i
      id = 0
      do i = 1, self%group_count
         if (self%groups(i)%name == to_lower(name)) id = i
      end do
   end function get_group_id

   subroutine set_watch_name(self, name)
      !! Type-bound to abstract_watch: Set name of watch
      class(abstract_watch), intent(inout) :: self
      character(len=*), intent(in) :: name
      !! Watch name
      self%name = name
   end subroutine set_watch_name
   
   subroutine start_timer_by_name(self, name)
      !! Type-bound to abstract_watch: Start timer referenced by name.
      !! Note: Wrapper of the corresponding routine from lightkrylov_timer.
      class(abstract_watch), intent(inout) :: self
      character(len=*), intent(in) :: name
      ! internal
      integer :: id
      character(len=128) :: tname
      tname = to_lower(name)
      id = self%get_timer_id(tname)
      if (id == 0) call timer_not_found(tname, 'start_timer_by_name')
      call self%timers(id)%start()
      call logger%log_debug('Timer "'//trim(tname)//'" started.', this_module, self%name)
   end subroutine start_timer_by_name

   subroutine stop_timer_by_name(self, name)
      !! Type-bound to abstract_watch: Stop timer referenced by name.
      !! Note: Wrapper of the corresponding routine from lightkrylov_timer.
      class(abstract_watch), intent(inout) :: self
      character(len=*), intent(in) :: name
      ! internal
      integer :: id
      character(len=128) :: tname
      tname = to_lower(name)
      id = self%get_timer_id(tname)
      if (id == 0) call timer_not_found(tname, 'stop_timer_by_name')
      call self%timers(id)%stop()
      call logger%log_debug('Timer "'//trim(tname)//'" stopped.', this_module, self%name)
   end subroutine stop_timer_by_name

   subroutine pause_timer_by_name(self, name)
      !! Type-bound to abstract_watch: Pause timer referenced by name.
      !! Note: Wrapper of the corresponding routine from lightkrylov_timer.
      class(abstract_watch), intent(inout) :: self
      character(len=*), intent(in) :: name
      ! internal
      integer :: id
      character(len=128) :: tname
      tname = to_lower(name)
      id = self%get_timer_id(tname)
      if (id == 0) call timer_not_found(tname, 'pause_timer_by_name')
      call self%timers(id)%pause()
      call logger%log_debug('Timer "'//trim(tname)//'" paused.', this_module, self%name)
   end subroutine

   subroutine reset_timer_by_name(self, name, soft, clean)
      !! Type-bound to abstract_watch: Reset timer referenced by name.
      !! Note: Wrapper of the corresponding routine from lightkrylov_timer.
      class(abstract_watch), intent(inout) :: self
      character(len=*), intent(in) :: name
      logical, optional, intent(in) :: soft
      logical, optional, intent(in) :: clean
      ! internal
      integer :: id
      character(len=128) :: tname
      tname = to_lower(name)
      id = self%get_timer_id(tname)
      if (id == 0) call timer_not_found(tname, 'reset_timer_by_name')
      call self%timers(id)%reset(soft, clean)
   end subroutine

   real(dp) function get_timer_etime_by_name(self, name, restart) result(etime)
      !! Type-bound to abstract_watch: Getter routine to return the current etime for timer referenced by name.
      !! Notes: Wrapper of the corresponding routine from lightkrylov_timer.
      !!        If it is running, the timer is stopped.
      class(abstract_watch) :: self
      character(len=*) :: name
      logical, optional :: restart
      ! internal
      integer :: id
      character(len=128) :: tname
      tname = to_lower(name)
      id = self%get_timer_id(tname)
      if (id == 0) call timer_not_found(tname, 'get_timer_etime_by_name')
      etime = self%timers(id)%get_time(restart)
   end function get_timer_etime_by_name

   subroutine get_timer_data_by_name(self, name, restart, etime, etmin, etmax, etimp, lcount, rcount, gcount)
      !! Type-bound to abstract_watch: Getter routine to return the data for timer referenced by name.
      !! Notes: Wrapper of the corresponding routine from lightkrylov_timer.
      !!        If it is running, the timer is stopped.
      class(abstract_watch), intent(inout) :: self
      character(len=*), intent(in) :: name
      logical, optional, intent(in) :: restart
      real(dp), optional, intent(out) :: etime
      real(dp), optional, intent(out) :: etmin
      real(dp), optional, intent(out) :: etmax
      real(dp), optional, intent(out) :: etimp
      integer, optional, intent(out) :: lcount
      integer, optional, intent(out) :: rcount
      integer, optional, intent(out) :: gcount
      ! internal
      integer :: id
      character(len=128) :: tname
      tname = to_lower(name)
      id = self%get_timer_id(tname)
      if (id == 0) then 
         call timer_not_found(tname, 'get_timer_data_by_name')
      else
         call self%timers(id)%get_data(restart, etime, etmin, etmax, etimp, lcount, rcount, gcount)
      end if
   end subroutine get_timer_data_by_name

   subroutine print_timer_info_by_name(self, name)
      !! Type-bound to abstract_watch: Print timing information for timer referenced by name.
      !! Note: Wrapper of the corresponding routine from lightkrylov_timer.
      class(abstract_watch), intent(inout) :: self
      character(len=*), intent(in) :: name
      ! internal
      integer :: id
      character(len=128) :: tname
      tname = to_lower(name)
      id = self%get_timer_id(tname)
      if (id == 0) then 
         call timer_not_found(tname, 'print_timer_info_by_name')
      else
         call self%timers(id)%print_info()
      end if
   end subroutine

   subroutine reset_all(self, soft, clean)
      !! Type-bound to abstract_watch: Utility function to reset all timers at once.
      !! Note: Wrapper of the corresponding routine from lightkrylov_timer.
      class(abstract_watch), intent(inout) :: self
      logical, optional, intent(in) :: soft
      logical, optional, intent(in) :: clean
      ! internal
      integer :: i
      logical :: soft_
      logical :: clean_
      character(len=128) :: msg
      soft_  = optval(soft, .true.)
      clean_ = optval(clean, .false.)
      do i = 1, self%timer_count
         call self%timers(i)%reset(soft, clean, verbose=.false.)
      end do
      write(msg,'(A,2(A,I0))') 'All timers reset: ', 'private: ', self%private_count, ', user: ', self%user_count
      call logger%log_message(msg, this_module, self%name)
      write(msg,'(A,L,3X,A,L)') 'soft reset: ', soft_, 'flush timers: ', clean_
      call logger%log_message(msg, this_module, self%name)
   end subroutine reset_all

   subroutine enumerate(self, only_user)
      !! Type-bound to abstract_watch: Summarize registered timers and their status.
      class(abstract_watch), intent(in) :: self
      logical, optional, intent(in) :: only_user
      !! Summarize only user defined timers? default: .false.
      ! internal
      integer :: i, j
      logical :: only_user_
      character(len=128) :: msg, fmt_e
      fmt_e = '(2X,I3,A50," :",3(1X,I0))'
      only_user_ = optval(only_user, .false.)
      if (.not. only_user_) then
         call logger%log_message('Registered timers: all', this_module, self%name)
         do i = 1, self%group_count
            call logger%log_message(trim(self%groups(i)%name)//":", this_module)
            do j = self%groups(i)%istart, self%groups(i)%iend
               associate (t => self%timers(j))
                  write(msg,fmt_e) j, trim(t%name), t%count, t%local_count, t%reset_count
                  call logger%log_message(msg, this_module)
               end associate
            end do
         end do
      end if
      if (self%user_count > 0) then
         call logger%log_message('Registered timers: user', this_module, self%name)
         do i = self%private_count+1, self%timer_count
            associate (t => self%timers(i))
               write(msg,fmt_e) i, trim(t%name), t%count, t%local_count, t%reset_count
               call logger%log_message(msg, this_module)
            end associate
         end do
      end if
   end subroutine enumerate

   subroutine initialize(self)
      !! Initialize global watch within LightKrylov and define private system timers.
      class(abstract_watch), intent(inout) :: self
      ! internal
      integer :: i, count
      character(len=128) :: msg
      if (.not. self%is_initialized) then
         call self%set_private_timers_and_name()
         self%private_count = self%timer_count
         write(msg,'(2(I0,A))') self%private_count, ' private timers registered in ', self%group_count, ' groups:'
         call logger%log_information(msg, this_module, self%name)
         do i = 1, self%group_count
            count = self%groups(i)%iend - self%groups(i)%istart + 1
            write(msg,'(3X,A20," : ",I3," timers.")') self%groups(i)%name, count
            call logger%log_information(msg, this_module, self%name)
         end do
         self%is_initialized = .true.
      else
         ! If the system timers have already been defined, we want to flush the data.
         call self%reset_all(soft = .false.)
         write(msg,'(3X,I4,A)') self%private_count, ' private timers registered and fully reset.'
         call logger%log_information(msg, this_module, self%name)
         if (self%user_count > 0) then
            write(msg,'(3X,I4,A)') self%user_count, ' user defined timers registered and fully reset.'
            call logger%log_information(msg, this_module, self%name)
         end if
      end if
      ! All subsequent timers that are added are user defined.
      self%user_mode = .true.
      call logger%log_message('Private timer initialization complete.', this_module, self%name)
   end subroutine initialize

   subroutine finalize(self, write_to_file)
      !! Finalize global watch within LightKrylov and print used timers.
      class(abstract_watch), intent(inout) :: self
      logical, optional, intent(in) :: write_to_file
      !! Print timer summary to a dedicated timer logfile as well as the general log file? default: .true.
      ! internal
      integer :: i, j, icalled, ic_user, tmr_unit
      integer, dimension(:), allocatable :: ic
      logical :: to_file
      character(len=128) :: msg, logfile_tmr
      to_file = optval(write_to_file, .true.)
      if (to_file) then
         write(logfile_tmr, '(A,A)') trim(to_lower(self%name)), '.log'
         call logger%add_log_file(logfile_tmr, tmr_unit, status='replace', action='write')
      end if
      icalled = 0
      allocate(ic(self%group_count))
      do i = 1, self%timer_count
         call self%timers(i)%finalize(if_silent=.true.)
         if (self%timers(i)%count > 0) icalled = icalled + 1
         do j = 1, self%group_count
            if (i == self%groups(j)%iend) then
               ic(j) = icalled
               icalled = 0
            end if
         end do
      end do
      ic_user = icalled
      call logger%log_message('Timer finalization complete.', this_module, self%name)
      call logger%log_message('#########   Global timer summary   ##########', this_module)
      call logger%log_message('Overview:', this_module, self%name)
      write(msg, '(2X,A60,I5)') 'Total active timers:', self%timer_count
      call logger%log_message(msg, this_module)
      write(msg, '(2X,A60,I5)') 'User defined:', self%user_count
      call logger%log_message(msg, this_module)
      write(msg, '(2X,A60,I5)') 'Called timers:', sum(ic) + ic_user
      call logger%log_message(msg, this_module)
      do i = 1, self%group_count
         if (ic(i) > 0) then
            associate(g => self%groups(i))
               call print_summary_header(g%name, self%name)
               do j = g%istart, g%iend
                  call print_summary(self%timers(j))
               end do
            end associate
         end if
      end do
      if (self%user_count > 0 .and. ic_user > 0) then
         call print_summary_header('User-defined', self%name)
         do i = self%private_count + 1, self%timer_count
            call print_summary(self%timers(i))
         end do
      end if
      call logger%log_message('#########   Global timer summary   ##########', this_module)
      if (to_file) call logger%remove_log_unit(tmr_unit, close_unit=.true.)
   end subroutine finalize

   !--------------------------------------------------------------
   !  Helper subroutines for pretty output
   !--------------------------------------------------------------

   subroutine print_summary_header(section_name, watch_name)
      !! Print section headers for the private and user defined timers.
      character(len=*), intent(in) :: section_name
      character(len=*), intent(in) :: watch_name
      ! internal
      character(len=128) :: msg
      call logger%log_message(trim(section_name)//':', this_module, watch_name)
      write(msg, fmt_h) 'name', 'calls', 'total (s)', 'avg (s)', 'min (s)', 'max (s)'
      call logger%log_message(msg, this_module)
      call logger%log_message('_____________________________________________', this_module)
   end subroutine print_summary_header

   subroutine print_summary(t)
      !! Print the full timer summary
      class(lightkrylov_timer), intent(in) :: t
      ! internal
      integer  :: i, count, count2
      real(dp) :: etime, etavg, etmin, etmax
      character(len=128) :: msg
      count  = sum(t%count_data)
      count2 = 0
      etmin  = huge(0.0_dp)
      etmax  = 0.0_dp
      do i = 1, t%reset_count
         if (t%count_data(i) > 0) then
            count2 = count2 + 1
            etmin = min(etmin, t%etmin_data(i))
            etmax = max(etmax, t%etmax_data(i))
         end if
      end do
      if (count == 1) then
         write(msg,fmt_1) trim(t%name), 'total', count, t%etime_data(1), '-', '-', '-'
         call logger%log_message(msg, this_module)
      else if (count > 1) then
         etime = sum(t%etime_data)
         etavg = sum(t%etavg_data)/count2
         write(msg,fmt_t) trim(t%name), 'total', count, etime, etavg, etmin, etmax
         call logger%log_message(msg, this_module)
         if (t%reset_count > 1) then
            do i = 1, t%reset_count
               etime = t%etime_data(i)
               etmin = t%etmin_data(i)
               etmax = t%etmax_data(i)
               etavg = t%etavg_data(i)
               count = t%count_data(i)
               if (count > 0) then
                  write(msg,fmt_r) 'reset', i, count, etime, etavg, etmin, etmax
               else
                  write(msg,fmt_n) 'reset', i, count, 'not called'
               end if
               call logger%log_message(msg, this_module)
            end do
         end if
      end if
   end subroutine print_summary

   !--------------------------------------------------------------
   !  Helper subroutines for error handling
   !--------------------------------------------------------------

   subroutine timer_not_found(name, procedure)
      character(len=*), intent(in) :: name
      character(len=*), optional, intent(in) :: procedure
      call stop_error('Timer "'//trim(name)//'" not found!', this_module, procedure)
   end subroutine timer_not_found

   subroutine element_exists(name, element, procedure)
      character(len=*), intent(in) :: name
      character(len=*), intent(in) :: element
      character(len=*), optional, intent(in) :: procedure
      call stop_error(trim(element)//' "'//trim(name)//'" already defined!', this_module, procedure)
   end subroutine element_exists

end module LightKrylov_Timer_Utils
