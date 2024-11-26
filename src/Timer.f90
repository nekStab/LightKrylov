module LightKrylov_Timing
   use stdlib_optval, only: optval
   use stdlib_ascii, only: to_lower
   use LightKrylov_Constants, only: dp
   use LightKrylov_Logger
   use LightKrylov_Timer_Utils
   implicit none
   private
   character(len=*), parameter :: this_module      = 'LK_Timer'
   character(len=*), parameter :: this_module_long = 'LightKrylov_Timer'
   logical :: if_time = .false.

   public :: global_lightkrylov_timer

   ! LightKrylov_watch type
   type, extends(abstract_watch), public :: lightkrylov_watch
      !! Global timing structure to contain all timers within Lightkrylov
      character(len=128) :: name = 'lightkrylov_timer'
   contains
      private
      procedure, pass(self), public :: set_private_timers => set_lightkrylov_timers
   end type lightkrylov_watch

   type(lightkrylov_watch) :: global_lightkrylov_timer

contains

   !--------------------------------------------------------------
   !  Concrete implementations for the lightkrylov_watch type
   !--------------------------------------------------------------

   subroutine set_lightkrylov_timers(self)
      !! Initialize global watch within LightKrylov and define private system timers.
      class(lightkrylov_watch), intent(inout) :: self
      ! internal
      integer :: count_old
      ! timers for LightKrylov_BaseKrylov
      count_old = self%get_timer_count()
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
      call self%add_group('BaseKrylov', istart=count_old+1, iend=self%get_timer_count())
      ! timers for LightKrylov_IterativeSolvers
      count_old = self%get_timer_count()
      ! rsp
      call self%add_timer('eigs_rsp')
      call self%add_timer('eighs_rsp')
      call self%add_timer('svds_rsp')
      call self%add_timer('gmres_rsp')
      call self%add_timer('fgmres_rsp')
      call self%add_timer('cg_rsp')
      ! rdp
      call self%add_timer('eigs_rdp')
      call self%add_timer('eighs_rdp')
      call self%add_timer('svds_rdp')
      call self%add_timer('gmres_rdp')
      call self%add_timer('fgmres_rdp')
      call self%add_timer('cg_rdp')
      ! csp
      call self%add_timer('eigs_csp')
      call self%add_timer('eighs_csp')
      call self%add_timer('svds_csp')
      call self%add_timer('gmres_csp')
      call self%add_timer('fgmres_csp')
      call self%add_timer('cg_csp')
      ! cdp
      call self%add_timer('eigs_cdp')
      call self%add_timer('eighs_cdp')
      call self%add_timer('svds_cdp')
      call self%add_timer('gmres_cdp')
      call self%add_timer('fgmres_cdp')
      call self%add_timer('cg_cdp')
      call self%add_group('IterativeSolvers', istart=count_old+1, iend=self%get_timer_count())
      ! timers for LightKrylov_NewtonKrylov
      count_old = self%get_timer_count()
      ! rsp
      call self%add_timer('newton_rsp')
      ! rdp
      call self%add_timer('newton_rdp')
      ! csp
      call self%add_timer('newton_csp')
      ! cdp
      call self%add_timer('newton_cdp')
      call self%add_group('NewtonKrylov', istart=count_old+1, iend=self%get_timer_count())
   end subroutine set_lightkrylov_timers

end module LightKrylov_Timing