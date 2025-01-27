module lightkrylov_gmres
    !--------------------------------------------
    !-----     Fortran Standard Library     -----
    !--------------------------------------------
    use iso_fortran_env, only: output_unit
    use stdlib_strings, only: padr
    use stdlib_optval, only: optval

    !-------------------------------
    !-----     LightKrylov     -----
    !-------------------------------
    use LightKrylov_Constants
    use LightKrylov_Logger
    use lightkrylov_Utils

    implicit none
    private
 
    character(len=*), parameter :: this_module      = 'LK_GMRES'
    character(len=*), parameter :: this_module_long = 'LightKrylov_GMRES'

    !----------------------------------------
    !-----     OPTIONS AND METADATA     -----
    !----------------------------------------

    type, extends(abstract_opts), public :: gmres_sp_opts
        integer :: kdim = 30
        integer :: maxiter = 10
        logical :: if_print_metadata = .false.
        logical :: sanity_check = .true.
    end type

    type, extends(abstract_metadata), public :: gmres_sp_metadata
        integer :: n_iter = 0
        integer :: n_inner = 0
        integer :: n_outer = 0
        real(sp), dimension(:), allocatable :: res
        logical :: converged = .false.
        integer :: info = 0
    contains
        procedure, pass(self), public :: print => print_gmres_sp
        procedure, pass(self), public :: reset => reset_gmres_sp
    end type
    type, extends(abstract_opts), public :: gmres_dp_opts
        integer :: kdim = 30
        integer :: maxiter = 10
        logical :: if_print_metadata = .false.
        logical :: sanity_check = .true.
    end type

    type, extends(abstract_metadata), public :: gmres_dp_metadata
        integer :: n_iter = 0
        integer :: n_inner = 0
        integer :: n_outer = 0
        real(dp), dimension(:), allocatable :: res
        logical :: converged = .false.
        integer :: info = 0
    contains
        procedure, pass(self), public :: print => print_gmres_dp
        procedure, pass(self), public :: reset => reset_gmres_dp
    end type

contains
    subroutine print_gmres_sp(self, reset_counters, verbose)
        class(gmres_sp_metadata), intent(inout) :: self
        logical, optional, intent(in) :: reset_counters
        !! Reset all counters to zero after printing?
        logical, optional, intent(in) :: verbose
        !! Print the residual full residual history?
        ! internals
        integer :: i
        logical :: ifreset, ifverbose
        character(len=128) :: msg

        ifreset   = optval(reset_counters, .false.)
        ifverbose = optval(verbose, .false.)
  
        write(msg,'(A30,I6,"  (",I6,"/",I3,")")') padr('Iterations   (inner/outer): ', 30), &
                  & self%n_iter, self%n_inner, self%n_outer
        call logger%log_message(msg, module=this_module, procedure='gmres_metadata')
        if (ifverbose) then
            write(msg,'(14X,A15)') 'Residual'
            call logger%log_message(msg, module=this_module, procedure='gmres_metadata')
            write(msg,'(A14,E15.8)') '   INIT:', self%res(1)
            call logger%log_message(msg, module=this_module, procedure='gmres_metadata')
            do i = 1, self%n_iter
               write(msg,'(A,I3,A,E20.8)') '   Step ', i, ': ', self%res(i)
               call logger%log_message(msg, module=this_module, procedure='gmres_metadata')
            end do
        else
            write(msg,'(A30,I20)') padr('Number of records: ', 30), size(self%res)
            call logger%log_message(msg, module=this_module, procedure='gmres_metadata')
            write(msg,'(A30,E20.8)') padr('Residual: ', 30), self%res(size(self%res))
            call logger%log_message(msg, module=this_module, procedure='gmres_metadata')
        end if
        if (self%converged) then
            call logger%log_message('Status: CONVERGED', module=this_module, procedure='gmres_metadata')
        else
            call logger%log_message('Status: NOT CONVERGED', module=this_module, procedure='gmres_metadata')
        end if
        if (ifreset) call self%reset()
        return
    end subroutine print_gmres_sp

    subroutine reset_gmres_sp(self)
        class(gmres_sp_metadata), intent(inout) :: self
        self%n_iter = 0
        self%n_inner = 0
        self%n_outer = 0
        self%converged = .false.
        self%info = 0
        if (allocated(self%res)) deallocate(self%res)
        return
    end subroutine reset_gmres_sp
    subroutine print_gmres_dp(self, reset_counters, verbose)
        class(gmres_dp_metadata), intent(inout) :: self
        logical, optional, intent(in) :: reset_counters
        !! Reset all counters to zero after printing?
        logical, optional, intent(in) :: verbose
        !! Print the residual full residual history?
        ! internals
        integer :: i
        logical :: ifreset, ifverbose
        character(len=128) :: msg

        ifreset   = optval(reset_counters, .false.)
        ifverbose = optval(verbose, .false.)
  
        write(msg,'(A30,I6,"  (",I6,"/",I3,")")') padr('Iterations   (inner/outer): ', 30), &
                  & self%n_iter, self%n_inner, self%n_outer
        call logger%log_message(msg, module=this_module, procedure='gmres_metadata')
        if (ifverbose) then
            write(msg,'(14X,A15)') 'Residual'
            call logger%log_message(msg, module=this_module, procedure='gmres_metadata')
            write(msg,'(A14,E15.8)') '   INIT:', self%res(1)
            call logger%log_message(msg, module=this_module, procedure='gmres_metadata')
            do i = 1, self%n_iter
               write(msg,'(A,I3,A,E20.8)') '   Step ', i, ': ', self%res(i)
               call logger%log_message(msg, module=this_module, procedure='gmres_metadata')
            end do
        else
            write(msg,'(A30,I20)') padr('Number of records: ', 30), size(self%res)
            call logger%log_message(msg, module=this_module, procedure='gmres_metadata')
            write(msg,'(A30,E20.8)') padr('Residual: ', 30), self%res(size(self%res))
            call logger%log_message(msg, module=this_module, procedure='gmres_metadata')
        end if
        if (self%converged) then
            call logger%log_message('Status: CONVERGED', module=this_module, procedure='gmres_metadata')
        else
            call logger%log_message('Status: NOT CONVERGED', module=this_module, procedure='gmres_metadata')
        end if
        if (ifreset) call self%reset()
        return
    end subroutine print_gmres_dp

    subroutine reset_gmres_dp(self)
        class(gmres_dp_metadata), intent(inout) :: self
        self%n_iter = 0
        self%n_inner = 0
        self%n_outer = 0
        self%converged = .false.
        self%info = 0
        if (allocated(self%res)) deallocate(self%res)
        return
    end subroutine reset_gmres_dp
end module
