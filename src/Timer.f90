module LightKrylov_Timer
   use stdlib_optval, only: optval
   use LightKrylov_Constants, only: dp
   use LightKrylov_Logger
   implicit none
   private
   character(len=128), parameter :: this_module = 'LightKrylov_Timer'
   logical :: if_time = .false.

   public :: time_lightkrylov
   public :: global_timer

   ! Timer type
   type, public :: timer_type
      private
      character(len=128) :: name
      real(dp) :: elapsed_time = 0.0_dp
      real(dp) :: start_time = 0.0_dp
      logical :: running = .false.
      integer :: count = 0
   contains
      procedure, pass(self), public :: start => start_timer
      procedure, pass(self), public :: stop => stop_timer
      procedure, pass(self), public :: pause => pause_timer
      procedure, pass(self), public :: reset => reset_timer
      procedure, pass(self), public :: get_time => get_timer_time
   end type timer_type

   ! Watch type
   type :: watch_type
      private
      type(timer_type), dimension(:), allocatable :: timers
      integer :: timer_count = 0
      integer :: basekrylov_count = 0
      integer :: iterativesolvers_count = 0
      integer :: newtonkrylov_count = 0
      integer :: private_count = 0
      logical :: user_mode = .false.
      integer :: user_count = 0
   contains
      procedure, pass(self), public :: add_timer
      procedure, pass(self), public :: remove_timer
      procedure, pass(self), public :: get_timer_id
      procedure, pass(self), public :: initialize
      procedure, pass(self), public :: enumerate
      procedure, pass(self), public :: finalize
      procedure, pass(self), public :: start => start_timer_by_name
      procedure, pass(self), public :: stop => stop_timer_by_name
      procedure, pass(self), public :: pause => pause_timer_by_name
      procedure, pass(self), public :: reset => reset_timer_by_name
   end type watch_type

   type(watch_type) :: global_timer

contains

   function time_lightkrylov() result(if_time_lightkrylov)
      logical :: if_time_lightkrylov
      if_time_lightkrylov = if_time
   end function time_lightkrylov

   subroutine start_timer(self)
      class(timer_type), intent(inout) :: self
      if (.not. self%running) then
         call cpu_time(self%start_time)
         self%running = .true.
         self%count = self%count + 1
      end if
   end subroutine start_timer

   subroutine stop_timer(self)
      class(timer_type), intent(inout) :: self
      ! internal
      real(dp) :: t_now
      call cpu_time(t_now)
      if (self%running) then
         self%elapsed_time = self%elapsed_time + (t_now - self%start_time)
         self%running = .false.
      end if
   end subroutine stop_timer

   subroutine pause_timer(self)
      class(timer_type), intent(inout) :: self
      ! internal
      real(dp) :: t_now
      call cpu_time(t_now)
      if (self%running) then
         self%elapsed_time = self%elapsed_time + (t_now - self%start_time)
         self%running = .false.
      end if
   end subroutine pause_timer

   subroutine reset_timer(self)
      class(timer_type), intent(inout) :: self
      self%elapsed_time = 0.0_dp
      self%start_time = 0.0_dp
      self%running = .false.
      self%count = 0
   end subroutine reset_timer

   real(dp) function get_timer_time(self) result(etime)
      class(timer_type), intent(inout) :: self
      if (self%running) then
         call self%stop()
      end if
      etime = self%elapsed_time
   end function

   integer function get_timer_id(self, name) result(id)
      class(watch_type) :: self
      character(len=*)  :: name
      ! internal
      integer :: i
      id = 0
      do i = 1, self%timer_count
         if (self%timers(i)%name == name) then
            id = i
         end if
      end do
   end function get_timer_id

   subroutine add_timer(self, name)
      class(watch_type), intent(inout) :: self
      character(len=*), intent(in) :: name
      if (self%timer_count == 0) then
         allocate(self%timers(1))
         self%timers(1) = timer_type(name)
         self%timer_count = 1
      else
         if (self%get_timer_id(name) > 0) then
            call stop_error('Timer "'//trim(name)//'" already defined!', module=this_module, procedure='add_timer')
         end if
         self%timers = [ self%timers, timer_type(name) ]
         self%timer_count = self%timer_count + 1
         if (self%user_mode) self%user_count = self%user_count + 1
      end if
   end subroutine add_timer

   subroutine remove_timer(self, name)
      class(watch_type), intent(inout) :: self
      character(len=*), intent(in) :: name
      ! internal
      type(timer_type), dimension(:), allocatable :: new_timers
      integer :: id
      id = self%get_timer_id(name)
      if (id == 0) then
         call stop_error('Timer "'//trim(name)//'" not defined!', module=this_module, procedure='remove_timer')
      else
         if (id <= self%private_count) then
            call logger%log_message('Cannot remove private timer "'//trim(name)//'".', module=this_module, procedure='remove_timer')
         else
            self%timer_count = self%timer_count - 1
            allocate(new_timers(self%timer_count))
            new_timers(1:id-1) = self%timers(1:id-1)
            new_timers(id:)    = self%timers(id+1:)
            deallocate(self%timers)
            self%timers = new_timers
         end if
      end if
   end subroutine remove_timer
   
   subroutine start_timer_by_name(self, name)
      class(watch_type), intent(inout) :: self
      character(len=*), intent(in) :: name
      ! internal
      integer :: id
      id = self%get_timer_id(name)
      if (id == 0) then 
         call stop_error('Timer "'//trim(name)//'" not found!', module=this_module, procedure='start_timer_by_name')
      else
         call self%timers(id)%start()
      end if
   end subroutine start_timer_by_name

   subroutine stop_timer_by_name(self, name)
      class(watch_type), intent(inout) :: self
      character(len=*), intent(in) :: name
      ! internal
      integer :: id
      id = self%get_timer_id(name)
      if (id == 0) then 
         call stop_error('Timer "'//trim(name)//'" not found!', module=this_module, procedure='start_timer_by_name')
      else
         call self%timers(id)%stop()
      end if
   end subroutine stop_timer_by_name

   subroutine pause_timer_by_name(self, name)
      class(watch_type), intent(inout) :: self
      character(len=*), intent(in) :: name
      ! internal
      integer :: id
      id = self%get_timer_id(name)
      if (id == 0) then 
         call stop_error('Timer "'//trim(name)//'" not found!', module=this_module, procedure='start_timer_by_name')
      else
         call self%timers(id)%start()
      end if
   end subroutine

   subroutine reset_timer_by_name(self, name)
      class(watch_type), intent(inout) :: self
      character(len=*), intent(in) :: name
      ! internal
      integer :: id
      id = self%get_timer_id(name)
      if (id == 0) then 
         call stop_error('Timer "'//trim(name)//'" not found!', module=this_module, procedure='start_timer_by_name')
      else
         call self%timers(id)%start()
      end if
   end subroutine

   subroutine enumerate(self, only_user)
      class(watch_type), intent(in) :: self
      logical, optional, intent(in) :: only_user
      ! internal
      integer :: i
      logical :: only_user_
      character(len=128) :: msg
      only_user_ = optval(only_user, .true.)
      if (.not. only_user_) then
         call logger%log_message('Registered timers: default', module=this_module)
         do i = 1, self%private_count
            write(msg,'(4X,I4,A,A)') i, ' : ', trim(self%timers(i)%name)
            call logger%log_message(msg, module=this_module)
         end do
      end if
      call logger%log_message('Registered timers: user', module=this_module)
      do i = self%private_count+1, self%timer_count
         write(msg,'(4X,I4,A,A)') i, ' : ', trim(self%timers(i)%name)
         call logger%log_message(msg, module=this_module)
      end do
   end subroutine enumerate

   subroutine initialize(self)
      class(watch_type), intent(inout) :: self
      ! timers for LightKrylov_BaseKrylov
      ! rsp
      call self%add_timer('QR_with_pivoting_rsp')
      call self%add_timer('QR_no_pivoting_rsp')
      call self%add_timer('Orthonormalize_basis_rsp')
      call self%add_timer('Orthonormalize_vector_against_basis_rsp')
      call self%add_timer('Orthonormalize_basis_against_basis_rsp')
      call self%add_timer('DGS_vector_against_basis_rsp')
      call self%add_timer('DGS_basis_against_basis_rsp')
      self%basekrylov_count = self%timer_count
      ! rdp
      call self%add_timer('QR_with_pivoting_rdp')
      call self%add_timer('QR_no_pivoting_rdp')
      call self%add_timer('Orthonormalize_basis_rdp')
      call self%add_timer('Orthonormalize_vector_against_basis_rdp')
      call self%add_timer('Orthonormalize_basis_against_basis_rdp')
      call self%add_timer('DGS_vector_against_basis_rdp')
      call self%add_timer('DGS_basis_against_basis_rdp')
      self%basekrylov_count = self%timer_count
      ! csp
      call self%add_timer('QR_with_pivoting_csp')
      call self%add_timer('QR_no_pivoting_csp')
      call self%add_timer('Orthonormalize_basis_csp')
      call self%add_timer('Orthonormalize_vector_against_basis_csp')
      call self%add_timer('Orthonormalize_basis_against_basis_csp')
      call self%add_timer('DGS_vector_against_basis_csp')
      call self%add_timer('DGS_basis_against_basis_csp')
      self%basekrylov_count = self%timer_count
      ! cdp
      call self%add_timer('QR_with_pivoting_cdp')
      call self%add_timer('QR_no_pivoting_cdp')
      call self%add_timer('Orthonormalize_basis_cdp')
      call self%add_timer('Orthonormalize_vector_against_basis_cdp')
      call self%add_timer('Orthonormalize_basis_against_basis_cdp')
      call self%add_timer('DGS_vector_against_basis_cdp')
      call self%add_timer('DGS_basis_against_basis_cdp')
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
   end subroutine initialize

   subroutine finalize(self)
      class(watch_type), intent(inout) :: self
      ! internal
      integer :: i, j, icalled
      integer :: ic_bk, ic_is, ic_nk, ic_user
      character(len=128) :: msg
      icalled = 0
      ic_bk = 0
      do i = 1, self%timer_count
         call self%timers(i)%stop()
         if (self%timers(i)%count > 0) icalled = icalled + 1
         if (i == self%basekrylov_count) then
            ic_bk = icalled
         else if (i == self%iterativesolvers_count) then
            ic_is = icalled - ic_bk
         else if (i == self%newtonkrylov_count) then
            ic_nk = icalled - ic_is
         end if
      end do
      ic_user = icalled - ic_nk
      if_time = .false.
      call logger%log_message('##############        Timer summary        ##############', module=this_module)
      write(msg, '(2X,A20,I3)') 'Total active timers:', self%timer_count
      call logger%log_message(msg, module=this_module)
      write(msg, '(2X,A20,I3)') 'User defined:', self%user_count
      call logger%log_message(msg, module=this_module)
      write(msg, '(2X,A20,I3)') 'Called timers:', icalled
      call logger%log_message(msg, module=this_module)
      write(msg, '(A32,A5,2(1X,A10))') 'name', 'cnt', 'total', 'avg'
      call logger%log_message(msg, module=this_module)
      if (ic_bk > 0) then
         call logger%log_message('BaseKrylov:', module=this_module)
         do i = 1, self%basekrylov_count
            associate(t => self%timers(i))
               if (t%count > 0) then
                  write(msg,'(2X,A30,I5,2(1X,F10.6))') trim(t%name), t%count, t%get_time(), t%get_time()/t%count
                  call logger%log_message(msg, module=this_module)
               end if
            end associate
         end do
      end if
      j = self%basekrylov_count
      if (ic_is > 0) then
         call logger%log_message('IterativeSolvers:', module=this_module)
         do i = j, self%iterativesolvers_count
            associate(t => self%timers(i))
               if (t%count > 0) then
                  write(msg,'(2X,A30,I5,2(1X,F10.6))') trim(t%name), t%count, t%get_time(), t%get_time()/t%count
                  call logger%log_message(msg, module=this_module)
               end if
            end associate
         end do
      end if
      j = self%iterativesolvers_count
      if (ic_nk > 0) then
         call logger%log_message('NewtonKrylov:', module=this_module)
         do i = j, self%newtonkrylov_count
            associate(t => self%timers(i))
               if (t%count > 0) then
                  write(msg,'(2X,A30,I5,2(1X,F10.6))') trim(t%name), t%count, t%get_time(), t%get_time()/t%count
                  call logger%log_message(msg, module=this_module)
               end if
            end associate
         end do
      end if
      if (self%user_count > 0 .and. ic_user > 0) then
         j = self%private_count
         call logger%log_message('User-defined:', module=this_module)
         do i = j, self%timer_count
            associate(t => self%timers(i))
               if (t%count > 0) then
                  write(msg,'(2X,A30,I5,2(1X,F10.6))') trim(t%name), t%count, t%get_time(), t%get_time()/t%count
                  call logger%log_message(msg, module=this_module)
               end if
            end associate
         end do
      end if
      call logger%log_message('##############        Timer summary        ##############', module=this_module)
   end subroutine finalize

end module LightKrylov_Timer