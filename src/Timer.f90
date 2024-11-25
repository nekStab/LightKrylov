module LightKrylov_Timing
   use stdlib_optval, only: optval
   use stdlib_ascii, only: to_lower
   use LightKrylov_Constants, only: dp
   use LightKrylov_Logger
   implicit none
   private
   character(len=*), parameter :: this_module      = 'LK_Timer'
   character(len=*), parameter :: this_module_long = 'LightKrylov_Timer'
   logical :: if_time = .false.

   public :: time_lightkrylov
   public :: global_lightkrylov_timer

   ! Timer type
   type, public :: lightkrylov_timer
      private
      character(len=128), public :: name = 'default_timer'
      real(dp) :: elapsed_time = 0.0_dp
      real(dp) :: start_time = 0.0_dp
      real(dp), dimension(:), allocatable :: etime_history
      real(dp), dimension(:), allocatable :: etavg_history
      integer,  dimension(:), allocatable :: count_history
      logical :: running = .false.
      integer :: count = 0
      integer :: reset_counter = 0
   contains
      private
      procedure, pass(self), public :: start => start_timer
      procedure, pass(self), public :: stop => stop_timer
      procedure, pass(self), public :: pause => pause_timer
      procedure, pass(self), public :: reset => reset_timer
      procedure, pass(self), public :: finalize => finalize_timer
      procedure, pass(self), public :: get_time => get_timer_time
      procedure, pass(self), public :: print_info => print_timer_info
      procedure, pass(self), public :: save_timer_history
   end type lightkrylov_timer

   ! Abstract watch type
   type, abstract, public :: abstract_watch
      !! Base type to define a global timer.
      private
      type(lightkrylov_timer), dimension(:), allocatable :: timers
      integer :: timer_count = 0
      integer :: private_count = 0
      logical :: user_mode = .false.
      integer :: user_count = 0
   contains
      private
      procedure, pass(self), public :: add_timer
      procedure, pass(self), public :: remove_timer
      procedure, pass(self), public :: get_timer_id
      procedure, pass(self), public :: enumerate
      procedure, pass(self), public :: reset_all
      procedure, pass(self), public :: start => start_timer_by_name
      procedure, pass(self), public :: stop => stop_timer_by_name
      procedure, pass(self), public :: pause => pause_timer_by_name
      procedure, pass(self), public :: reset => reset_timer_by_name
      procedure(abstract_watch_init), pass(self), deferred, public :: initialize
      procedure(abstract_watch_exit), pass(self), deferred, public :: finalize
   end type abstract_watch

   abstract interface
      subroutine abstract_watch_init(self)
         !! Interface for the initialization of the structure.
         import abstract_watch
         class(abstract_watch), intent(inout) :: self
      end subroutine abstract_watch_init
      subroutine abstract_watch_exit(self)
         !! Interface for the finalization of the structure including the printing of the results
         import abstract_watch
         class(abstract_watch), intent(inout) :: self
      end subroutine abstract_watch_exit
   end interface

   ! LightKrylov_watch type
   type, extends(abstract_watch), public :: lightkrylov_watch
      !! Global timing structure to contain all timers within Lightkrylov
      character(len=128) :: name = 'lightkrylov_timer'
      integer :: basekrylov_count = 0
      integer :: iterativesolvers_count = 0
      integer :: newtonkrylov_count = 0
   contains
      private
      procedure, pass(self), public :: initialize => initialize_lightkrylov_watch
      procedure, pass(self), public :: finalize => finalize_lightkrylov_watch
   end type lightkrylov_watch

   type(lightkrylov_watch) :: global_lightkrylov_timer
contains

   logical function time_lightkrylov() result(if_time_lightkrylov)
      if_time_lightkrylov = if_time
   end function time_lightkrylov

   !--------------------------------------------------------------
   !  Type-bound procedures for lightkrylov_timer type
   !--------------------------------------------------------------

   subroutine start_timer(self)
      class(lightkrylov_timer), intent(inout) :: self
      if (.not. self%running) then
         call cpu_time(self%start_time)
         self%running = .true.
         self%count = self%count + 1
      end if
   end subroutine start_timer

   subroutine stop_timer(self)
      class(lightkrylov_timer), intent(inout) :: self
      ! internal
      real(dp) :: t_now
      call cpu_time(t_now)
      if (self%running) then
         self%elapsed_time = self%elapsed_time + (t_now - self%start_time)
         self%running = .false.
      end if
   end subroutine stop_timer

   subroutine pause_timer(self)
      class(lightkrylov_timer), intent(inout) :: self
      ! internal
      real(dp) :: t_now
      call cpu_time(t_now)
      if (self%running) then
         self%elapsed_time = self%elapsed_time + (t_now - self%start_time)
         self%running = .false.
      end if
   end subroutine pause_timer

   subroutine save_timer_history(self)
      class(lightkrylov_timer), intent(inout) :: self
      if (self%reset_counter == 0) then
         allocate(self%etime_history(1))
         allocate(self%etavg_history(1))
         allocate(self%count_history(1))
         if (self%count > 0) then
            self%etime_history(1) = self%elapsed_time
            self%etavg_history(1) = self%elapsed_time/self%count
            self%count_history(1) = self%count
         else
            self%etime_history(1) = 0.0_dp
            self%etavg_history(1) = 0.0_dp
            self%count_history(1) = self%count
         end if
         self%reset_counter = 1
      else
         if (self%count > 0) then
            self%etime_history = [ self%etime_history, self%elapsed_time ]
            self%etavg_history = [ self%etavg_history, self%elapsed_time/self%count ]
            self%count_history = [ self%count_history, self%count ]
         else
            self%etime_history(1) = 0.0_dp
            self%etavg_history(1) = 0.0_dp
            self%count_history(1) = self%count
         end if
         self%reset_counter = self%reset_counter + 1
      end if
   end subroutine save_timer_history

   subroutine reset_timer(self, save_history)
      class(lightkrylov_timer), intent(inout) :: self
      logical, optional, intent(in) :: save_history
      ! internal
      logical :: ifsave
      ifsave = optval(save_history, .true.)
      if (self%count > 0) then
         if (ifsave) call self%save_timer_history()
         self%elapsed_time = 0.0_dp
         self%start_time = 0.0_dp
         self%running = .false.
         self%count = 0
      end if
   end subroutine reset_timer

   real(dp) function get_timer_time(self) result(etime)
      class(lightkrylov_timer), intent(inout) :: self
      if (self%running) then
         call self%stop()
      end if
      etime = self%elapsed_time
   end function

   subroutine print_timer_info(self)
      class(lightkrylov_timer), intent(inout) :: self
      ! internal
      integer :: i
      real(dp) :: etime, etavg
      character(len=128) :: msg, timer_fmt
      timer_fmt       = '(2X,A30," : ",I7,2(1X,F12.6))'
      call logger%log_message('###        Timer info                  #######################################', & 
                              & module=this_module)
      write(msg, '(A32," : ",A7,2(1X,A12))') 'name', 'calls', 'total (s)', 'avg (s)'
      call logger%log_message(msg, module=this_module)
      etime = 0.0_dp
      etavg = 0.0_dp
      etime = self%get_time()
      if (self%count > 0) etavg = etime/self%count
      write(msg,timer_fmt) trim(self%name), self%count, etime, etavg
      call logger%log_message(msg, module=this_module)
      call logger%log_message('###        Timer info                  #######################################', & 
                              & module=this_module)
   end subroutine print_timer_info

   subroutine finalize_timer(self)
      class(lightkrylov_timer), intent(inout) :: self
      ! internal
      integer :: i
      integer :: ic_bk, ic_is, ic_nk, ic_user, count
      real(dp) :: etime, etavg
      character(len=128) :: msg, timer_fmt, timer_fmt_reset
      timer_fmt       = '(2X,A30," : ",A6,1X,I7,2(1X,F12.6))'
      timer_fmt_reset = '(2X,33X,A6,I3,1X,I7,2(1X,F12.6))'
      call logger%log_message('###        Timer summary               #######################################', & 
                              & module=this_module)
      write(msg, '(A32," : ",7X,A7,2(1X,A12))') 'name', 'calls', 'total (s)', 'avg (s)'
      call logger%log_message(msg, module=this_module)
      call logger%log_message('______________________________________________________________________________', & 
                              & module=this_module)
      etavg = 0.0_dp
      if (self%count > 0 .or. self%reset_counter > 0) then
         if (self%reset_counter == 0) call self%save_timer_history()
         etime = sum(self%etime_history)
         count = sum(self%count_history)
         etavg = sum(self%etavg_history)/self%reset_counter
         write(msg,timer_fmt) trim(self%name), 'total', count, etime, etavg
         call logger%log_message(msg, module=this_module)
         if (self%reset_counter > 1) then
            do i = 1, self%reset_counter
               etime = self%etime_history(i)
               etavg = self%etavg_history(i)
               count = self%count_history(i)
               write(msg,timer_fmt_reset) 'reset', i, count, etime, etavg
               call logger%log_message(msg, module=this_module)
            end do
         end if
      end if
      call logger%log_message('###        Timer summary               #######################################', & 
                              & module=this_module)
   end subroutine finalize_timer

   !--------------------------------------------------------------
   !  Type-bound procedures for abstract_watch type
   !--------------------------------------------------------------

   integer function get_timer_id(self, name) result(id)
      class(abstract_watch) :: self
      character(len=*)  :: name
      ! internal
      integer :: i
      id = 0
      do i = 1, self%timer_count
         if (self%timers(i)%name == to_lower(name)) then
            id = i
         end if
      end do
   end function get_timer_id

   subroutine add_timer(self, name)
      class(abstract_watch), intent(inout) :: self
      character(len=*), intent(in) :: name
      if (self%timer_count == 0) then
         allocate(self%timers(1))
         self%timers(1) = lightkrylov_timer(to_lower(name))
         self%timer_count = 1
      else
         if (self%get_timer_id(name) > 0) then
            call stop_error('Timer "'//to_lower(trim(name))//'" already defined!', & 
                              & module=this_module, procedure='add_timer')
         end if
         self%timers = [ self%timers, lightkrylov_timer(name) ]
         self%timer_count = self%timer_count + 1
         if (self%user_mode) self%user_count = self%user_count + 1
      end if
      call logger%log_debug('Timer "'//to_lower(trim(name))//'" added.', module=this_module)
   end subroutine add_timer

   subroutine remove_timer(self, name)
      class(abstract_watch), intent(inout) :: self
      character(len=*), intent(in) :: name
      ! internal
      type(lightkrylov_timer), dimension(:), allocatable :: new_timers
      integer :: id
      id = self%get_timer_id(name)
      if (id == 0) then
         call stop_error('Timer "'//to_lower(trim(name))//'" not defined!', & 
                              & module=this_module, procedure='remove_timer')
      else
         if (id <= self%private_count) then
            call logger%log_message('Cannot remove private timer "'//to_lower(trim(name))//'".', & 
                              & module=this_module, procedure='remove_timer')
         else
            self%timer_count = self%timer_count - 1
            allocate(new_timers(self%timer_count))
            new_timers(1:id-1) = self%timers(1:id-1)
            new_timers(id:)    = self%timers(id+1:)
            deallocate(self%timers)
            self%timers = new_timers
         end if
      end if
      call logger%log_debug('Timer "'//to_lower(trim(name))//'" removed.', module=this_module)
   end subroutine remove_timer
   
   subroutine start_timer_by_name(self, name)
      class(abstract_watch), intent(inout) :: self
      character(len=*), intent(in) :: name
      ! internal
      integer :: id
      id = self%get_timer_id(name)
      if (id == 0) then 
         call stop_error('Timer "'//to_lower(trim(name))//'" not found!', & 
                              & module=this_module, procedure='start_timer_by_name')
      else
         call self%timers(id)%start()
      end if
      call logger%log_debug('Timer "'//to_lower(trim(name))//'" started.', module=this_module)
   end subroutine start_timer_by_name

   subroutine stop_timer_by_name(self, name)
      class(abstract_watch), intent(inout) :: self
      character(len=*), intent(in) :: name
      ! internal
      integer :: id
      id = self%get_timer_id(name)
      if (id == 0) then 
         call stop_error('Timer "'//to_lower(trim(name))//'" not found!', & 
                              & module=this_module, procedure='start_timer_by_name')
      else
         call self%timers(id)%stop()
      end if
      call logger%log_debug('Timer "'//to_lower(trim(name))//'" stopped.', module=this_module)
   end subroutine stop_timer_by_name

   subroutine pause_timer_by_name(self, name)
      class(abstract_watch), intent(inout) :: self
      character(len=*), intent(in) :: name
      ! internal
      integer :: id
      id = self%get_timer_id(name)
      if (id == 0) then 
         call stop_error('Timer "'//to_lower(trim(name))//'" not found!', & 
                              & module=this_module, procedure='start_timer_by_name')
      else
         call self%timers(id)%pause()
      end if
      call logger%log_debug('Timer "'//to_lower(trim(name))//'" paused.', module=this_module)
   end subroutine

   subroutine reset_timer_by_name(self, name, save_history)
      class(abstract_watch), intent(inout) :: self
      character(len=*), intent(in) :: name
      logical, optional, intent(in) :: save_history
      ! internal
      integer :: id
      id = self%get_timer_id(name)
      if (id == 0) then 
         call stop_error('Timer "'//to_lower(trim(name))//'" not found!', & 
                              & module=this_module, procedure='start_timer_by_name')
      else
         call self%timers(id)%reset(save_history)
      end if
   end subroutine

   subroutine enumerate(self, only_user)
      class(abstract_watch), intent(in) :: self
      logical, optional, intent(in) :: only_user
      ! internal
      integer :: i
      logical :: only_user_
      character(len=128) :: msg
      only_user_ = optval(only_user, .true.)
      if (.not. only_user_) then
         call logger%log_message('Registered timers: all', module=this_module)
         do i = 1, self%private_count
            write(msg,'(4X,I4,A,A)') i, ' : ', trim(self%timers(i)%name)
            call logger%log_message(msg, module=this_module)
         end do
      end if
      if (self%user_count > 0) then
         call logger%log_message('Registered timers: user', module=this_module)
         do i = self%private_count+1, self%timer_count
            write(msg,'(4X,I4,A,A)') i, ' : ', trim(self%timers(i)%name)
            call logger%log_message(msg, module=this_module)
         end do
      end if
   end subroutine enumerate

   subroutine reset_all(self, save_history)
      class(abstract_watch), intent(inout) :: self
      logical, optional, intent(in) :: save_history
      ! internal
      integer :: i
      character(len=128) :: msg
      do i = 1, self%timer_count
         call self%timers(i)%reset(save_history)
      end do
   end subroutine reset_all

   !--------------------------------------------------------------
   !  Concrete implementations for the lightkrylov_watch type
   !--------------------------------------------------------------

   subroutine initialize_lightkrylov_watch(self)
      class(lightkrylov_watch), intent(inout) :: self
      ! timers for LightKrylov_BaseKrylov
      ! rsp
      call self%add_timer('qr_with_pivoting_rsp')
      call self%add_timer('qr_no_pivoting_rsp')
      call self%add_timer('orthonormalize_basis_rsp')
      call self%add_timer('orthonormalize_vector_against_basis_rsp')
      call self%add_timer('orthonormalize_basis_against_basis_rsp')
      call self%add_timer('dgs_vector_against_basis_rsp')
      call self%add_timer('dgs_basis_against_basis_rsp')
      call self%add_timer('arnoldi_rsp')
      call self%add_timer('lanczos_bidiagonalization_rsp')
      call self%add_timer('lanczos_tridiagonalization_rsp')
      self%basekrylov_count = self%timer_count
      ! rdp
      call self%add_timer('qr_with_pivoting_rdp')
      call self%add_timer('qr_no_pivoting_rdp')
      call self%add_timer('orthonormalize_basis_rdp')
      call self%add_timer('orthonormalize_vector_against_basis_rdp')
      call self%add_timer('orthonormalize_basis_against_basis_rdp')
      call self%add_timer('dgs_vector_against_basis_rdp')
      call self%add_timer('dgs_basis_against_basis_rdp')
      call self%add_timer('arnoldi_rdp')
      call self%add_timer('lanczos_bidiagonalization_rdp')
      call self%add_timer('lanczos_tridiagonalization_rdp')
      self%basekrylov_count = self%timer_count
      ! csp
      call self%add_timer('qr_with_pivoting_csp')
      call self%add_timer('qr_no_pivoting_csp')
      call self%add_timer('orthonormalize_basis_csp')
      call self%add_timer('orthonormalize_vector_against_basis_csp')
      call self%add_timer('orthonormalize_basis_against_basis_csp')
      call self%add_timer('dgs_vector_against_basis_csp')
      call self%add_timer('dgs_basis_against_basis_csp')
      call self%add_timer('arnoldi_csp')
      call self%add_timer('lanczos_bidiagonalization_csp')
      call self%add_timer('lanczos_tridiagonalization_csp')
      self%basekrylov_count = self%timer_count
      ! cdp
      call self%add_timer('qr_with_pivoting_cdp')
      call self%add_timer('qr_no_pivoting_cdp')
      call self%add_timer('orthonormalize_basis_cdp')
      call self%add_timer('orthonormalize_vector_against_basis_cdp')
      call self%add_timer('orthonormalize_basis_against_basis_cdp')
      call self%add_timer('dgs_vector_against_basis_cdp')
      call self%add_timer('dgs_basis_against_basis_cdp')
      call self%add_timer('arnoldi_cdp')
      call self%add_timer('lanczos_bidiagonalization_cdp')
      call self%add_timer('lanczos_tridiagonalization_cdp')
      self%basekrylov_count = self%timer_count
      ! timers for LightKrylov_IterativeSolvers
      ! rsp
      call self%add_timer('eigs_rsp')
      call self%add_timer('eighs_rsp')
      call self%add_timer('svds_rsp')
      call self%add_timer('gmres_rsp')
      call self%add_timer('fgmres_rsp')
      call self%add_timer('cg_rsp')
      self%iterativesolvers_count = self%timer_count
      ! rdp
      call self%add_timer('eigs_rdp')
      call self%add_timer('eighs_rdp')
      call self%add_timer('svds_rdp')
      call self%add_timer('gmres_rdp')
      call self%add_timer('fgmres_rdp')
      call self%add_timer('cg_rdp')
      self%iterativesolvers_count = self%timer_count
      ! csp
      call self%add_timer('eigs_csp')
      call self%add_timer('eighs_csp')
      call self%add_timer('svds_csp')
      call self%add_timer('gmres_csp')
      call self%add_timer('fgmres_csp')
      call self%add_timer('cg_csp')
      self%iterativesolvers_count = self%timer_count
      ! cdp
      call self%add_timer('eigs_cdp')
      call self%add_timer('eighs_cdp')
      call self%add_timer('svds_cdp')
      call self%add_timer('gmres_cdp')
      call self%add_timer('fgmres_cdp')
      call self%add_timer('cg_cdp')
      self%iterativesolvers_count = self%timer_count
      ! timers for LightKrylov_NewtonKrylov
      ! rsp
      call self%add_timer('newton_rsp')
      ! rdp
      call self%add_timer('newton_rdp')
      ! csp
      call self%add_timer('newton_csp')
      ! cdp
      call self%add_timer('newton_cdp')
      self%newtonkrylov_count = self%timer_count
      self%private_count = self%timer_count
      self%user_mode = .true.
      if_time = .true.
      call logger%log_message('LightKrylov system timer initialization complete.', module=this_module)
   end subroutine initialize_lightkrylov_watch

   subroutine finalize_lightkrylov_watch(self)
      class(lightkrylov_watch), intent(inout) :: self
      ! internal
      integer :: i, j, icalled
      integer :: ic_bk, ic_is, ic_nk, ic_user, count, rcount
      real(dp) :: etime, etavg
      character(len=128) :: msg, timer_fmt, timer_fmt_reset
      icalled = 0
      do i = 1, self%timer_count
         call self%timers(i)%stop()
         if (self%timers(i)%count > 0) icalled = icalled + 1
         call self%timers(i)%save_timer_history()
         if (i == self%basekrylov_count) then
            ic_bk = icalled
         else if (i == self%iterativesolvers_count) then
            ic_is = icalled - ic_bk
         else if (i == self%newtonkrylov_count) then
            ic_nk = icalled - ic_is - ic_bk
         end if
      end do
      ic_user = icalled - ic_nk - ic_is - ic_bk
      if_time = .false.
      call logger%log_message('LightKrylov timer finalization complete.', module=this_module)
      call logger%log_message('###        Global timer summary        #######################################', & 
                              & module=this_module)
      call logger%log_message('____________________', module=this_module)
      call logger%log_message('Overview:', module=this_module)
      write(msg, '(2X,A40,I5)') 'Total active timers:', self%timer_count
      call logger%log_message(msg, module=this_module)
      write(msg, '(2X,A40,I5)') 'User defined:', self%user_count
      call logger%log_message(msg, module=this_module)
      write(msg, '(2X,A40,I5)') 'Called timers:', icalled
      call logger%log_message(msg, module=this_module)
      timer_fmt = '(2X,A30," : ",A6,1X,I7,2(1X,F12.6))'
      timer_fmt_reset = '(2X,33X,A6,I3,1X,I7,2(1X,F12.6))'
      if (ic_bk > 0) then
         call logger%log_message('____________________', module=this_module)
         call logger%log_message('BaseKrylov:', module=this_module)
         write(msg, '(A32," : ",7X,A7,2(1X,A12))') 'name', 'calls', 'total (s)', 'avg (s)'
         call logger%log_message(msg, module=this_module)
         call logger%log_message('______________________________________________________________________________', & 
                              & module=this_module)
         do i = 1, self%basekrylov_count
            associate(t => self%timers(i))
               rcount = t%reset_counter
               if (t%count_history(rcount) > 0) then
                  etime = sum(t%etime_history)
                  etavg = sum(t%etavg_history)/rcount
                  count = sum(t%count_history)
                  write(msg,timer_fmt) trim(t%name), 'total', count, etime, etavg
                  call logger%log_message(msg, module=this_module)
                  if (rcount > 1) then
                     do j = 1, rcount
                        etime = t%etime_history(j)
                        etavg = t%etavg_history(j)
                        count = t%count_history(j)
                        write(msg,timer_fmt_reset) 'reset', j, count, etime, etavg
                        call logger%log_message(msg, module=this_module)
                     end do
                  end if
               end if
            end associate
         end do
      end if
      j = self%basekrylov_count
      if (ic_is > 0) then
         call logger%log_message('____________________', module=this_module)
         call logger%log_message('IterativeSolvers:', module=this_module)
         write(msg, '(A32," : ",7X,A7,2(1X,A12))') 'name', 'calls', 'total (s)', 'avg (s)'
         call logger%log_message(msg, module=this_module)
         call logger%log_message('______________________________________________________________________________', & 
                              & module=this_module)
         do i = j, self%iterativesolvers_count
            associate(t => self%timers(i))
               rcount = t%reset_counter
               if (t%count_history(rcount) > 0) then
                  etime = sum(t%etime_history)
                  etavg = sum(t%etavg_history)/rcount
                  count = sum(t%count_history)
                  write(msg,timer_fmt) trim(t%name), 'total', count, etime, etavg
                  call logger%log_message(msg, module=this_module)
                  if (rcount > 1) then
                     do j = 1, rcount
                        etime = t%etime_history(j)
                        etavg = t%etavg_history(j)
                        count = t%count_history(j)
                        write(msg,timer_fmt_reset) 'reset', j, count, etime, etavg
                        call logger%log_message(msg, module=this_module)
                     end do
                  end if
               end if
            end associate
         end do
      end if
      j = self%iterativesolvers_count
      if (ic_nk > 0) then
         call logger%log_message('____________________', module=this_module)
         call logger%log_message('NewtonKrylov:', module=this_module)
         write(msg, '(A32," : ",7X,A7,2(1X,A12))') 'name', 'calls', 'total (s)', 'avg (s)'
         call logger%log_message(msg, module=this_module)
         call logger%log_message('______________________________________________________________________________', & 
                              & module=this_module)
         do i = j, self%newtonkrylov_count
            associate(t => self%timers(i))
               rcount = t%reset_counter
               if (t%count_history(rcount) > 0) then
                  etime = sum(t%etime_history)
                  etavg = sum(t%etavg_history)/rcount
                  count = sum(t%count_history)
                  write(msg,timer_fmt) trim(t%name), 'total', count, etime, etavg
                  call logger%log_message(msg, module=this_module)
                  if (rcount > 1) then
                     do j = 1, rcount
                        etime = t%etime_history(j)
                        etavg = t%etavg_history(j)
                        count = t%count_history(j)
                        write(msg,timer_fmt_reset) 'reset', j, count, etime, etavg
                        call logger%log_message(msg, module=this_module)
                     end do
                  end if
               end if
            end associate
         end do
      end if
      if (self%user_count > 0 .and. ic_user > 0) then
         j = self%private_count
         call logger%log_message('____________________', module=this_module)
         call logger%log_message('User-defined:', module=this_module)
         write(msg, '(A32," : ",7X,A7,2(1X,A12))') 'name', 'calls', 'total (s)', 'avg (s)'
         call logger%log_message(msg, module=this_module)
         call logger%log_message('______________________________________________________________________________', & 
                              & module=this_module)
         do i = j, self%timer_count
            associate(t => self%timers(i))
               rcount = t%reset_counter
               if (t%count_history(rcount) > 0) then
                  etime = sum(t%etime_history)
                  etavg = sum(t%etavg_history)/rcount
                  count = sum(t%count_history)
                  write(msg,timer_fmt) trim(t%name), 'total', count, etime, etavg
                  call logger%log_message(msg, module=this_module)
                  if (rcount > 1) then
                     do j = 1, rcount
                        etime = t%etime_history(j)
                        etavg = t%etavg_history(j)
                        count = t%count_history(j)
                        write(msg,timer_fmt_reset) 'reset', j, count, etime, etavg
                        call logger%log_message(msg, module=this_module)
                     end do
                  end if
               end if
            end associate
         end do
      end if
      call logger%log_message('###        Global timer summary        #######################################',  & 
                              & module=this_module)
   end subroutine finalize_lightkrylov_watch

end module LightKrylov_Timing