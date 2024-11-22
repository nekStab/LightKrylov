module LightKrylov_Timer
   use LightKrylov_Constants, only: dp
   use LightKrylov_Logger
   implicit none
   private
   character(len=128), parameter :: this_module = 'LightKrylov_Timer'
   logical :: if_time = .false.

   public :: initialize_timers, enumerate_timers, finalize_timers
   public :: timeit

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
      integer :: private_count = 0
   contains
      procedure, pass(self), public :: add_timer
      procedure, pass(self), public :: remove_timer
      procedure, pass(self), public :: enumerate_timers
      procedure, pass(self), public :: start => start_timer_by_name
      procedure, pass(self), public :: stop => stop_timer_by_name
      procedure, pass(self), public :: pause => pause_timer_by_name
      procedure, pass(self), public :: reset => reset_timer_by_name
   end type watch_type

   type(watch_type), public :: global_timer

contains

   function timeit() result(if_time_lightkrylov)
      logical :: if_time_lightkrylov
      if_time_lightkrylov = if_time
   end function timeit

   subroutine initialize_timers()
      ! timers for LightKrylov_BaseKrylov
      ! rsp
      !call global_timer%add_timer('QR_with_pivoting_rsp')
      !call global_timer%add_timer('QR_no_pivoting_rsp')
      !call global_timer%add_timer('Orthonormalize_basis_rsp')
      !call global_timer%add_timer('Orthonormalize_vector_against_basis_rsp')
      !call global_timer%add_timer('Orthonormalize_basis_against_basis_rsp')
      !call global_timer%add_timer('DGS_vector_against_basis_rsp')
      !call global_timer%add_timer('DGS_basis_against_basis_rsp')
      !#:endfor
      ! timers for LightKrylov_IterativeSolvers
      !#:for kind, type in RC_KINDS_TYPES
      ! rsp
      !call global_timer%add_timer('eigs_rsp')
      !call global_timer%add_timer('eighs_rsp')
      call global_timer%add_timer('svds_rsp')
      call global_timer%add_timer('gmres_rsp')
      !call global_timer%add_timer('fgmres_rsp')
      !call global_timer%add_timer('cg_rsp')
      ! rdp
      !call global_timer%add_timer('QR_with_pivoting_rdp')
      !call global_timer%add_timer('QR_no_pivoting_rdp')
      !call global_timer%add_timer('Orthonormalize_basis_rdp')
      !call global_timer%add_timer('Orthonormalize_vector_against_basis_rdp')
      !call global_timer%add_timer('Orthonormalize_basis_against_basis_rdp')
      !call global_timer%add_timer('DGS_vector_against_basis_rdp')
      !call global_timer%add_timer('DGS_basis_against_basis_rdp')
      !#:endfor
      ! timers for LightKrylov_IterativeSolvers
      !#:for kind, type in RC_KINDS_TYPES
      ! rdp
      !call global_timer%add_timer('eigs_rdp')
      !call global_timer%add_timer('eighs_rdp')
      call global_timer%add_timer('svds_rdp')
      call global_timer%add_timer('gmres_rdp')
      !call global_timer%add_timer('fgmres_rdp')
      !call global_timer%add_timer('cg_rdp')
      ! csp
      !call global_timer%add_timer('QR_with_pivoting_csp')
      !call global_timer%add_timer('QR_no_pivoting_csp')
      !call global_timer%add_timer('Orthonormalize_basis_csp')
      !call global_timer%add_timer('Orthonormalize_vector_against_basis_csp')
      !call global_timer%add_timer('Orthonormalize_basis_against_basis_csp')
      !call global_timer%add_timer('DGS_vector_against_basis_csp')
      !call global_timer%add_timer('DGS_basis_against_basis_csp')
      !#:endfor
      ! timers for LightKrylov_IterativeSolvers
      !#:for kind, type in RC_KINDS_TYPES
      ! csp
      !call global_timer%add_timer('eigs_csp')
      !call global_timer%add_timer('eighs_csp')
      call global_timer%add_timer('svds_csp')
      call global_timer%add_timer('gmres_csp')
      !call global_timer%add_timer('fgmres_csp')
      !call global_timer%add_timer('cg_csp')
      ! cdp
      !call global_timer%add_timer('QR_with_pivoting_cdp')
      !call global_timer%add_timer('QR_no_pivoting_cdp')
      !call global_timer%add_timer('Orthonormalize_basis_cdp')
      !call global_timer%add_timer('Orthonormalize_vector_against_basis_cdp')
      !call global_timer%add_timer('Orthonormalize_basis_against_basis_cdp')
      !call global_timer%add_timer('DGS_vector_against_basis_cdp')
      !call global_timer%add_timer('DGS_basis_against_basis_cdp')
      !#:endfor
      ! timers for LightKrylov_IterativeSolvers
      !#:for kind, type in RC_KINDS_TYPES
      ! cdp
      !call global_timer%add_timer('eigs_cdp')
      !call global_timer%add_timer('eighs_cdp')
      call global_timer%add_timer('svds_cdp')
      call global_timer%add_timer('gmres_cdp')
      !call global_timer%add_timer('fgmres_cdp')
      !call global_timer%add_timer('cg_cdp')
      ! timers for LightKrylov_NewtonKrylov
      ! rsp
      call global_timer%add_timer('newton_rsp')
      ! rdp
      call global_timer%add_timer('newton_rdp')
      ! csp
      call global_timer%add_timer('newton_csp')
      ! cdp
      call global_timer%add_timer('newton_cdp')
      global_timer%private_count = global_timer%timer_count
      if_time = .true.
   end subroutine initialize_timers

   subroutine finalize_timers()
      ! internal
      integer :: i
      do i = 1, global_timer%timer_count
         call global_timer%timers(i)%stop()
      end do
      if_time = .false.
      do i = 1, global_timer%timer_count
         print *, trim(global_timer%timers(i)%name), ':', global_timer%timers(i)%count, global_timer%timers(i)%get_time()
      end do
   end subroutine finalize_timers

   subroutine start_timer(self)
      class(timer_type), intent(inout) :: self
      if (.not. self%running) then
         call cpu_time(self%start_time)
         self%running = .true.
      end if
   end subroutine start_timer

   subroutine stop_timer(self)
      class(timer_type), intent(inout) :: self
      ! internal
      real(dp) :: t_now
      call cpu_time(t_now)
      if (self%running) then
         self%elapsed_time = self%elapsed_time + (t_now - self%start_time)
         self%count = self%count + 1
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

   subroutine add_timer(self, name)
      class(watch_type), intent(inout) :: self
      character(len=*), intent(in) :: name
      ! internal
      integer :: i
      if (self%timer_count == 0) then
         allocate(self%timers(1))
         self%timers(1) = timer_type(name)
         self%timer_count = 1
      else
         do i = 1, self%timer_count
            if (self%timers(i)%name == name) then
               call stop_error('Timer "'//trim(name)//'" already defined!', module=this_module, procedure='add_timer')
            end if
         end do
         self%timers = [ self%timers, timer_type(name) ]
         self%timer_count = self%timer_count + 1
      end if
   end subroutine add_timer

   subroutine remove_timer(self, name)
      class(watch_type), intent(inout) :: self
      character(len=*), intent(in) :: name
      ! internal
      type(timer_type), dimension(:), allocatable :: new_timers
      integer :: i, i_rm
      logical :: found
      found = .false.
      do i = 1, self%timer_count
         if (self%timers(i)%name == name) then
            i_rm = i
            found = .true.
         end if
      end do
      if (.not. found) then
         call stop_error('Timer "'//trim(name)//'" not defined!', module=this_module, procedure='remove_timer')
      else
         if (i_rm <= self%private_count) then
            call logger%log_message('Cannot remove private timer "'//trim(name)//'".', module=this_module, procedure='remove_timer')
         else
            self%timer_count = self%timer_count - 1
            allocate(new_timers(self%timer_count))
            new_timers(1:i_rm-1) = self%timers(1:i_rm-1)
            new_timers(i_rm:)    = self%timers(i_rm+1:)
            deallocate(self%timers)
            self%timers = new_timers
         end if
      end if
   end subroutine remove_timer

   subroutine enumerate_timers(self)
      class(watch_type), intent(in) :: self
      ! internal
      integer :: i
      do i = 1, self%timer_count
         print *, 'Timer', i, ':', trim(self%timers(i)%name)
      end do
   end subroutine
   
   subroutine start_timer_by_name(self, name)
      class(watch_type), intent(inout) :: self
      character(len=*), intent(in) :: name
      ! internal
      integer :: i
      logical :: found
      found = .false.
      do i = 1, self%timer_count
         if (self%timers(i)%name == name) then
            call self%timers(i)%start()
            found = .true.
         end if
      end do
      if (.not. found) call stop_error('Timer "'//trim(name)//'" not found!', module=this_module, procedure='start_timer_by_name')
   end subroutine start_timer_by_name

   subroutine stop_timer_by_name(self, name)
      class(watch_type), intent(inout) :: self
      character(len=*), intent(in) :: name
      ! internal
      integer :: i
      logical :: found
      found = .false.
      do i = 1, self%timer_count
         if (self%timers(i)%name == name) then
            call self%timers(i)%stop()
            found = .true.
         end if
      end do
      if (.not. found) call stop_error('Timer "'//trim(name)//'" not found!', module=this_module, procedure='stop_timer_by_name')
   end subroutine stop_timer_by_name

   subroutine pause_timer_by_name(self, name)
      class(watch_type), intent(inout) :: self
      character(len=*), intent(in) :: name
      ! internal
      integer :: i
      logical :: found
      found = .false.
      do i = 1, self%timer_count
         if (self%timers(i)%name == name) then
            call self%timers(i)%pause()
            found = .true.
         end if
      end do
      if (.not. found) call stop_error('Timer "'//trim(name)//'" not found!', module=this_module, procedure='pause_timer_by_name')
   end subroutine

   subroutine reset_timer_by_name(self, name)
      class(watch_type), intent(inout) :: self
      character(len=*), intent(in) :: name
      ! internal
      integer :: i
      logical :: found
      found = .false.
      do i = 1, self%timer_count
         if (self%timers(i)%name == name) then
            call self%timers(i)%reset()
            found = .true.
         end if
      end do
      if (.not. found) call stop_error('Timer "'//trim(name)//'" not found!', module=this_module, procedure='reset_timer_by_name')
   end subroutine

end module LightKrylov_Timer