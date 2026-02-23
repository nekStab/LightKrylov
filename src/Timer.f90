module LightKrylov_Timing
   !!  This module provides a set of utility functions to define default timers within `LightKrylov`.

   !--------------------------------------------
   !-----     Standard Fortran Library     -----
   !--------------------------------------------
   use stdlib_optval, only: optval
   use stdlib_ascii, only: to_lower

   !-------------------------------
   !-----     LightKrylov     -----
   !-------------------------------
   use LightKrylov_Logger
   use LightKrylov_Timer_Utils

   implicit none(type, external)
   private

   character(len=*), parameter :: this_module      = 'LK_Timer'
   character(len=*), parameter :: this_module_long = 'LightKrylov_Timer'
   logical :: if_time = .false. 

   !----------------------------------
   !-----     Public exports     -----
   !----------------------------------

   public :: time_lightkrylov
   public :: global_lightkrylov_timer

   ! LightKrylov_watch type
   type, extends(abstract_watch), public :: lightkrylov_watch
      !! Global timing structure to contain all timers within Lightkrylov
   contains
      private
      procedure, pass(self), public :: set_private_timers_and_name => set_lightkrylov_timers
   end type lightkrylov_watch

   type(lightkrylov_watch) :: global_lightkrylov_timer

contains

   pure logical function time_lightkrylov() result(if_time_lightkrylov)
      if_time_lightkrylov = if_time
   end function time_lightkrylov

   subroutine set_lightkrylov_timer_switch(value)
      logical, intent(in) :: value     
      if (if_time .neqv. value) then
         if_time = value
         if (if_time) then
            call log_message('LightKrylov timing enabled.', module=this_module)
         else
            call log_message('LightKrylov timing disabled.', module=this_module)
         end if
      else
         call log_debug('LightKrylov timing switched unchanged.', module=this_module)
      end if      
   end subroutine set_lightkrylov_timer_switch

   !--------------------------------------------------------------
   !  Concrete implementations for the lightkrylov_watch type
   !--------------------------------------------------------------

   subroutine set_lightkrylov_timers(self)
      !! Initialize global watch within LightKrylov and define private system timers.
      class(lightkrylov_watch), intent(inout) :: self
      ! Internal
      integer :: istart, iend
      call self%set_watch_name('LightKrylov_timer')

      ! Timers for LAPACK
      call self%add_timer('eig', count=istart)
      call self%add_timer('eigh')
      call self%add_timer('svd')
      call self%add_timer('trsen', count=iend)
      ! Define LAPACK group
      call self%add_group('LAPACK', istart=istart, iend=iend)

      ! Timers for Utils
      ! rsp
      call self%add_timer('sqrtm_rsp', count=istart)
      call self%add_timer('expm_rsp')
      ! rdp
      call self%add_timer('sqrtm_rdp')
      call self%add_timer('expm_rdp')
      ! csp
      call self%add_timer('sqrtm_csp')
      call self%add_timer('expm_csp')
      ! cdp
      call self%add_timer('sqrtm_cdp')
      call self%add_timer('expm_cdp', count=iend)
      ! Define Utils group
      call self%add_group('Utils', istart=istart, iend=iend)

      ! Timers for BaseKrylov
      ! rsp
      call self%add_timer('qr_with_pivoting_rsp', count=istart)
      call self%add_timer('qr_no_pivoting_rsp')
      call self%add_timer('orthonormalize_basis_rsp')
      call self%add_timer('orthogonalize_vector_against_basis_rsp')
      call self%add_timer('orthogonalize_basis_against_basis_rsp')
      call self%add_timer('dgs_vector_against_basis_rsp')
      call self%add_timer('dgs_basis_against_basis_rsp')
      call self%add_timer('arnoldi_rsp')
      call self%add_timer('lanczos_bidiagonalization_rsp')
      call self%add_timer('lanczos_tridiagonalization_rsp')
      call self%add_timer('krylov_schur_rsp')
      ! rdp
      call self%add_timer('qr_with_pivoting_rdp')
      call self%add_timer('qr_no_pivoting_rdp')
      call self%add_timer('orthonormalize_basis_rdp')
      call self%add_timer('orthogonalize_vector_against_basis_rdp')
      call self%add_timer('orthogonalize_basis_against_basis_rdp')
      call self%add_timer('dgs_vector_against_basis_rdp')
      call self%add_timer('dgs_basis_against_basis_rdp')
      call self%add_timer('arnoldi_rdp')
      call self%add_timer('lanczos_bidiagonalization_rdp')
      call self%add_timer('lanczos_tridiagonalization_rdp')
      call self%add_timer('krylov_schur_rdp')
      ! csp
      call self%add_timer('qr_with_pivoting_csp')
      call self%add_timer('qr_no_pivoting_csp')
      call self%add_timer('orthonormalize_basis_csp')
      call self%add_timer('orthogonalize_vector_against_basis_csp')
      call self%add_timer('orthogonalize_basis_against_basis_csp')
      call self%add_timer('dgs_vector_against_basis_csp')
      call self%add_timer('dgs_basis_against_basis_csp')
      call self%add_timer('arnoldi_csp')
      call self%add_timer('lanczos_bidiagonalization_csp')
      call self%add_timer('lanczos_tridiagonalization_csp')
      call self%add_timer('krylov_schur_csp')
      ! cdp
      call self%add_timer('qr_with_pivoting_cdp')
      call self%add_timer('qr_no_pivoting_cdp')
      call self%add_timer('orthonormalize_basis_cdp')
      call self%add_timer('orthogonalize_vector_against_basis_cdp')
      call self%add_timer('orthogonalize_basis_against_basis_cdp')
      call self%add_timer('dgs_vector_against_basis_cdp')
      call self%add_timer('dgs_basis_against_basis_cdp')
      call self%add_timer('arnoldi_cdp')
      call self%add_timer('lanczos_bidiagonalization_cdp')
      call self%add_timer('lanczos_tridiagonalization_cdp')
      call self%add_timer('krylov_schur_cdp', count=iend)
      ! Define BaseKrylov group
      call self%add_group('BaseKrylov', istart=istart, iend=iend)

      ! Timers for IterativeSolvers
      ! rsp
      call self%add_timer('eigs_rsp', count=istart)
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
      call self%add_timer('cg_cdp', count=iend)
      ! Define IterativeSolvers group
      call self%add_group('IterativeSolvers', istart=istart, iend=iend)

      ! Timers for NewtonKrylov
      ! rsp
      call self%add_timer('newton_rsp', count=istart)
      ! rdp
      call self%add_timer('newton_rdp')
      ! csp
      call self%add_timer('newton_csp')
      ! cdp
      call self%add_timer('newton_cdp', count=iend)
      ! Define NewtonKrylov group
      call self%add_group('NewtonKrylov', istart=istart, iend=iend)

      ! Enable timing
      call set_lightkrylov_timer_switch(.true.)
   end subroutine set_lightkrylov_timers

end module LightKrylov_Timing
