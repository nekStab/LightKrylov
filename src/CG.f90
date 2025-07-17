submodule (lightkrylov_iterativesolvers) cg_solver
    use stdlib_strings, only: padr
    implicit none
contains

    !----------------------------------------
    !-----     Options and Metadata     -----
    !----------------------------------------

    module procedure print_cg_sp
        character(len=*), parameter :: this_procedure = 'print_cg_sp'
        integer :: i
        logical :: ifreset, ifverbose
        character(len=128) :: msg

        ifreset   = optval(reset_counters, .false.)
        ifverbose = optval(verbose, .false.)

        write(msg,'(A30,I20)') padr('Iterations: ', 30), self%n_iter
        call log_message(msg, this_module, this_procedure)
        if (ifverbose) then
            write(msg,'(14X,A15)') 'Residual'
            call log_message(msg, this_module, this_procedure)
            call log_message('Residual history:', this_module, this_procedure)
            write(msg,'(A14,E15.8)') '   INIT:', self%res(1)
            call log_message(msg, this_module, this_procedure)
            do i = 2, self%n_iter+1
               write(msg,'(A,I4,A,E15.8)') '   Step ', i-1, ': ', self%res(i)
               call log_message(msg, this_module, this_procedure)
            end do
        else
            write(msg,'(A30,I20)') padr('Number of records: ', 30), size(self%res)
            call log_message(msg, this_module, this_procedure)
            write(msg,'(A30,E20.8)') padr('Residual: ', 30), self%res(size(self%res))
            call log_message(msg, this_module, this_procedure)
        end if
        if (self%converged) then
            call log_message('Status: CONVERGED', this_module, this_procedure)
        else
            call log_message('Status: NOT CONVERGED', this_module, this_procedure)
        end if
        if (ifreset) call self%reset()
        return
    end procedure 

    module procedure reset_cg_sp
        self%n_iter = 0
        self%converged = .false.
        self%info = 0
        if (allocated(self%res)) deallocate(self%res)
    end procedure
    module procedure print_cg_dp
        character(len=*), parameter :: this_procedure = 'print_cg_dp'
        integer :: i
        logical :: ifreset, ifverbose
        character(len=128) :: msg

        ifreset   = optval(reset_counters, .false.)
        ifverbose = optval(verbose, .false.)

        write(msg,'(A30,I20)') padr('Iterations: ', 30), self%n_iter
        call log_message(msg, this_module, this_procedure)
        if (ifverbose) then
            write(msg,'(14X,A15)') 'Residual'
            call log_message(msg, this_module, this_procedure)
            call log_message('Residual history:', this_module, this_procedure)
            write(msg,'(A14,E15.8)') '   INIT:', self%res(1)
            call log_message(msg, this_module, this_procedure)
            do i = 2, self%n_iter+1
               write(msg,'(A,I4,A,E15.8)') '   Step ', i-1, ': ', self%res(i)
               call log_message(msg, this_module, this_procedure)
            end do
        else
            write(msg,'(A30,I20)') padr('Number of records: ', 30), size(self%res)
            call log_message(msg, this_module, this_procedure)
            write(msg,'(A30,E20.8)') padr('Residual: ', 30), self%res(size(self%res))
            call log_message(msg, this_module, this_procedure)
        end if
        if (self%converged) then
            call log_message('Status: CONVERGED', this_module, this_procedure)
        else
            call log_message('Status: NOT CONVERGED', this_module, this_procedure)
        end if
        if (ifreset) call self%reset()
        return
    end procedure 

    module procedure reset_cg_dp
        self%n_iter = 0
        self%converged = .false.
        self%info = 0
        if (allocated(self%res)) deallocate(self%res)
    end procedure

    !-------------------------------------------------
    !-----     CG SOLVERS FOR ABSTRACT TYPES     -----
    !-------------------------------------------------

    module procedure cg_rsp
        ! Options.
        integer :: maxiter
        real(sp) :: tol, rtol_, atol_
        type(cg_sp_opts)     :: opts
        type(cg_sp_metadata) :: cg_meta

        ! Working variables.
        class(abstract_vector_rsp), allocatable :: r, p, Ap, z
        real(sp) :: alpha, beta, r_dot_r_old, r_dot_r_new
        real(sp) :: residual

        ! Miscellaneous.
        character(len=*), parameter :: this_procedure = 'cg_rsp'
        integer :: i
        character(len=256) :: msg

        call log_debug('start', this_module, this_procedure)
        if (time_lightkrylov()) call timer%start(this_procedure)
        ! Deals with the optional args.
        rtol_ = optval(rtol, rtol_sp)
        atol_ = optval(atol, atol_sp)
        if (present(options)) then
            opts = options
        else
            opts = cg_sp_opts()
        endif
        tol = atol_ + rtol_ * b%norm() ; maxiter = opts%maxiter

        ! Initialize vectors.
        allocate(r, mold=b)  ; call r%zero()
        allocate(p, mold=b)  ; call p%zero()
        allocate(Ap, mold=b) ; call Ap%zero()

         ! Initialize meta & reset matvec counter
        cg_meta = cg_sp_metadata()
        call A%reset_counter(.false., 'cg%init')

        info = 0

        associate(ifprecond => present(preconditioner))
        ! Compute initial residual r = b - Ax.
        if (x%norm() > 0) call A%apply_matvec(x, r)
        call r%sub(b) ; call r%chsgn()

        ! Deal with the preconditioner (if available).
        if (ifprecond) then
            z = r ; call preconditioner%apply(z) ; p = z
            r_dot_r_old = r%dot(z)
        else
            p = r ; r_dot_r_old = r%dot(r)
        endif

        allocate(cg_meta%res(1)); cg_meta%res(1) = sqrt(abs(r_dot_r_old))

        ! Conjugate gradient iteration.
        cg_loop: do i = 1, maxiter
            ! Compute A @ p
            call A%apply_matvec(p, Ap)
            ! Compute step size.
            alpha = r_dot_r_old / p%dot(Ap)
            ! Update solution x = x + alpha*p
            call x%axpby(alpha, p, one_rsp)
            ! Update residual r = r - alpha*Ap
            call r%axpby(-alpha, Ap, one_rsp)

            if(ifprecond) then
                z = r ; call preconditioner%apply(z) ; r_dot_r_new = r%dot(z)
            else
                ! Compute new dot product of residual r_dot_r_new = r' * r.
                r_dot_r_new = r%dot(r)
            endif

            ! Check for convergence.
            residual = sqrt(r_dot_r_new)

            ! Save metadata.
            cg_meta%n_iter = cg_meta%n_iter + 1
            cg_meta%res = [ cg_meta%res, residual ]

            if (residual < tol) then
               cg_meta%converged = .true.
               exit cg_loop
            end if

            ! Compute new direction beta = r_dot_r_new / r_dot_r_old.
            beta = r_dot_r_new / r_dot_r_old
            ! Update direction p = beta*p + r

            if (ifprecond) then
                call p%axpby(one_rsp, z, beta)
            else
                call p%axpby(one_rsp, r, beta)
            endif

            ! Update r_dot_r_old for next iteration.
            r_dot_r_old = r_dot_r_new

            write(msg,'(A,I3,2(A,E9.2))') 'CG step ', i, ': res= ', residual, ', tol= ', tol
            call log_information(msg, this_module, this_procedure)
        enddo cg_loop
        end associate

        ! Set and copy info flag for completeness
        info = cg_meta%n_iter
        cg_meta%info = info

        if (opts%if_print_metadata) call cg_meta%print()

        ! Set metadata output
        if (present(meta)) then
            select type(meta)
            type is (cg_sp_metadata)
                meta = cg_meta
            class default
                call type_error('meta','cg_sp_metadata','OUT',this_module,this_procedure)
            end select
        end if

        call A%reset_counter(.false., 'cg%post')
        if (time_lightkrylov()) call timer%stop(this_procedure)
        call log_debug('end', this_module, this_procedure)

        return
    end procedure
    module procedure cg_rdp
        ! Options.
        integer :: maxiter
        real(dp) :: tol, rtol_, atol_
        type(cg_dp_opts)     :: opts
        type(cg_dp_metadata) :: cg_meta

        ! Working variables.
        class(abstract_vector_rdp), allocatable :: r, p, Ap, z
        real(dp) :: alpha, beta, r_dot_r_old, r_dot_r_new
        real(dp) :: residual

        ! Miscellaneous.
        character(len=*), parameter :: this_procedure = 'cg_rdp'
        integer :: i
        character(len=256) :: msg

        call log_debug('start', this_module, this_procedure)
        if (time_lightkrylov()) call timer%start(this_procedure)
        ! Deals with the optional args.
        rtol_ = optval(rtol, rtol_dp)
        atol_ = optval(atol, atol_dp)
        if (present(options)) then
            opts = options
        else
            opts = cg_dp_opts()
        endif
        tol = atol_ + rtol_ * b%norm() ; maxiter = opts%maxiter

        ! Initialize vectors.
        allocate(r, mold=b)  ; call r%zero()
        allocate(p, mold=b)  ; call p%zero()
        allocate(Ap, mold=b) ; call Ap%zero()

         ! Initialize meta & reset matvec counter
        cg_meta = cg_dp_metadata()
        call A%reset_counter(.false., 'cg%init')

        info = 0

        associate(ifprecond => present(preconditioner))
        ! Compute initial residual r = b - Ax.
        if (x%norm() > 0) call A%apply_matvec(x, r)
        call r%sub(b) ; call r%chsgn()

        ! Deal with the preconditioner (if available).
        if (ifprecond) then
            z = r ; call preconditioner%apply(z) ; p = z
            r_dot_r_old = r%dot(z)
        else
            p = r ; r_dot_r_old = r%dot(r)
        endif

        allocate(cg_meta%res(1)); cg_meta%res(1) = sqrt(abs(r_dot_r_old))

        ! Conjugate gradient iteration.
        cg_loop: do i = 1, maxiter
            ! Compute A @ p
            call A%apply_matvec(p, Ap)
            ! Compute step size.
            alpha = r_dot_r_old / p%dot(Ap)
            ! Update solution x = x + alpha*p
            call x%axpby(alpha, p, one_rdp)
            ! Update residual r = r - alpha*Ap
            call r%axpby(-alpha, Ap, one_rdp)

            if(ifprecond) then
                z = r ; call preconditioner%apply(z) ; r_dot_r_new = r%dot(z)
            else
                ! Compute new dot product of residual r_dot_r_new = r' * r.
                r_dot_r_new = r%dot(r)
            endif

            ! Check for convergence.
            residual = sqrt(r_dot_r_new)

            ! Save metadata.
            cg_meta%n_iter = cg_meta%n_iter + 1
            cg_meta%res = [ cg_meta%res, residual ]

            if (residual < tol) then
               cg_meta%converged = .true.
               exit cg_loop
            end if

            ! Compute new direction beta = r_dot_r_new / r_dot_r_old.
            beta = r_dot_r_new / r_dot_r_old
            ! Update direction p = beta*p + r

            if (ifprecond) then
                call p%axpby(one_rdp, z, beta)
            else
                call p%axpby(one_rdp, r, beta)
            endif

            ! Update r_dot_r_old for next iteration.
            r_dot_r_old = r_dot_r_new

            write(msg,'(A,I3,2(A,E9.2))') 'CG step ', i, ': res= ', residual, ', tol= ', tol
            call log_information(msg, this_module, this_procedure)
        enddo cg_loop
        end associate

        ! Set and copy info flag for completeness
        info = cg_meta%n_iter
        cg_meta%info = info

        if (opts%if_print_metadata) call cg_meta%print()

        ! Set metadata output
        if (present(meta)) then
            select type(meta)
            type is (cg_dp_metadata)
                meta = cg_meta
            class default
                call type_error('meta','cg_dp_metadata','OUT',this_module,this_procedure)
            end select
        end if

        call A%reset_counter(.false., 'cg%post')
        if (time_lightkrylov()) call timer%stop(this_procedure)
        call log_debug('end', this_module, this_procedure)

        return
    end procedure
    module procedure cg_csp
        ! Options.
        integer :: maxiter
        real(sp) :: tol, rtol_, atol_
        type(cg_sp_opts)     :: opts
        type(cg_sp_metadata) :: cg_meta

        ! Working variables.
        class(abstract_vector_csp), allocatable :: r, p, Ap, z
        complex(sp) :: alpha, beta, r_dot_r_old, r_dot_r_new
        real(sp) :: residual

        ! Miscellaneous.
        character(len=*), parameter :: this_procedure = 'cg_csp'
        integer :: i
        character(len=256) :: msg

        call log_debug('start', this_module, this_procedure)
        if (time_lightkrylov()) call timer%start(this_procedure)
        ! Deals with the optional args.
        rtol_ = optval(rtol, rtol_sp)
        atol_ = optval(atol, atol_sp)
        if (present(options)) then
            opts = options
        else
            opts = cg_sp_opts()
        endif
        tol = atol_ + rtol_ * b%norm() ; maxiter = opts%maxiter

        ! Initialize vectors.
        allocate(r, mold=b)  ; call r%zero()
        allocate(p, mold=b)  ; call p%zero()
        allocate(Ap, mold=b) ; call Ap%zero()

         ! Initialize meta & reset matvec counter
        cg_meta = cg_sp_metadata()
        call A%reset_counter(.false., 'cg%init')

        info = 0

        associate(ifprecond => present(preconditioner))
        ! Compute initial residual r = b - Ax.
        if (x%norm() > 0) call A%apply_matvec(x, r)
        call r%sub(b) ; call r%chsgn()

        ! Deal with the preconditioner (if available).
        if (ifprecond) then
            z = r ; call preconditioner%apply(z) ; p = z
            r_dot_r_old = r%dot(z)
        else
            p = r ; r_dot_r_old = r%dot(r)
        endif

        allocate(cg_meta%res(1)); cg_meta%res(1) = sqrt(abs(r_dot_r_old))

        ! Conjugate gradient iteration.
        cg_loop: do i = 1, maxiter
            ! Compute A @ p
            call A%apply_matvec(p, Ap)
            ! Compute step size.
            alpha = r_dot_r_old / p%dot(Ap)
            ! Update solution x = x + alpha*p
            call x%axpby(alpha, p, one_csp)
            ! Update residual r = r - alpha*Ap
            call r%axpby(-alpha, Ap, one_csp)

            if(ifprecond) then
                z = r ; call preconditioner%apply(z) ; r_dot_r_new = r%dot(z)
            else
                ! Compute new dot product of residual r_dot_r_new = r' * r.
                r_dot_r_new = r%dot(r)
            endif

            ! Check for convergence.
            residual = sqrt(abs(r_dot_r_new))

            ! Save metadata.
            cg_meta%n_iter = cg_meta%n_iter + 1
            cg_meta%res = [ cg_meta%res, residual ]

            if (residual < tol) then
               cg_meta%converged = .true.
               exit cg_loop
            end if

            ! Compute new direction beta = r_dot_r_new / r_dot_r_old.
            beta = r_dot_r_new / r_dot_r_old
            ! Update direction p = beta*p + r

            if (ifprecond) then
                call p%axpby(one_csp, z, beta)
            else
                call p%axpby(one_csp, r, beta)
            endif

            ! Update r_dot_r_old for next iteration.
            r_dot_r_old = r_dot_r_new

            write(msg,'(A,I3,2(A,E9.2))') 'CG step ', i, ': res= ', residual, ', tol= ', tol
            call log_information(msg, this_module, this_procedure)
        enddo cg_loop
        end associate

        ! Set and copy info flag for completeness
        info = cg_meta%n_iter
        cg_meta%info = info

        if (opts%if_print_metadata) call cg_meta%print()

        ! Set metadata output
        if (present(meta)) then
            select type(meta)
            type is (cg_sp_metadata)
                meta = cg_meta
            class default
                call type_error('meta','cg_sp_metadata','OUT',this_module,this_procedure)
            end select
        end if

        call A%reset_counter(.false., 'cg%post')
        if (time_lightkrylov()) call timer%stop(this_procedure)
        call log_debug('end', this_module, this_procedure)

        return
    end procedure
    module procedure cg_cdp
        ! Options.
        integer :: maxiter
        real(dp) :: tol, rtol_, atol_
        type(cg_dp_opts)     :: opts
        type(cg_dp_metadata) :: cg_meta

        ! Working variables.
        class(abstract_vector_cdp), allocatable :: r, p, Ap, z
        complex(dp) :: alpha, beta, r_dot_r_old, r_dot_r_new
        real(dp) :: residual

        ! Miscellaneous.
        character(len=*), parameter :: this_procedure = 'cg_cdp'
        integer :: i
        character(len=256) :: msg

        call log_debug('start', this_module, this_procedure)
        if (time_lightkrylov()) call timer%start(this_procedure)
        ! Deals with the optional args.
        rtol_ = optval(rtol, rtol_dp)
        atol_ = optval(atol, atol_dp)
        if (present(options)) then
            opts = options
        else
            opts = cg_dp_opts()
        endif
        tol = atol_ + rtol_ * b%norm() ; maxiter = opts%maxiter

        ! Initialize vectors.
        allocate(r, mold=b)  ; call r%zero()
        allocate(p, mold=b)  ; call p%zero()
        allocate(Ap, mold=b) ; call Ap%zero()

         ! Initialize meta & reset matvec counter
        cg_meta = cg_dp_metadata()
        call A%reset_counter(.false., 'cg%init')

        info = 0

        associate(ifprecond => present(preconditioner))
        ! Compute initial residual r = b - Ax.
        if (x%norm() > 0) call A%apply_matvec(x, r)
        call r%sub(b) ; call r%chsgn()

        ! Deal with the preconditioner (if available).
        if (ifprecond) then
            z = r ; call preconditioner%apply(z) ; p = z
            r_dot_r_old = r%dot(z)
        else
            p = r ; r_dot_r_old = r%dot(r)
        endif

        allocate(cg_meta%res(1)); cg_meta%res(1) = sqrt(abs(r_dot_r_old))

        ! Conjugate gradient iteration.
        cg_loop: do i = 1, maxiter
            ! Compute A @ p
            call A%apply_matvec(p, Ap)
            ! Compute step size.
            alpha = r_dot_r_old / p%dot(Ap)
            ! Update solution x = x + alpha*p
            call x%axpby(alpha, p, one_cdp)
            ! Update residual r = r - alpha*Ap
            call r%axpby(-alpha, Ap, one_cdp)

            if(ifprecond) then
                z = r ; call preconditioner%apply(z) ; r_dot_r_new = r%dot(z)
            else
                ! Compute new dot product of residual r_dot_r_new = r' * r.
                r_dot_r_new = r%dot(r)
            endif

            ! Check for convergence.
            residual = sqrt(abs(r_dot_r_new))

            ! Save metadata.
            cg_meta%n_iter = cg_meta%n_iter + 1
            cg_meta%res = [ cg_meta%res, residual ]

            if (residual < tol) then
               cg_meta%converged = .true.
               exit cg_loop
            end if

            ! Compute new direction beta = r_dot_r_new / r_dot_r_old.
            beta = r_dot_r_new / r_dot_r_old
            ! Update direction p = beta*p + r

            if (ifprecond) then
                call p%axpby(one_cdp, z, beta)
            else
                call p%axpby(one_cdp, r, beta)
            endif

            ! Update r_dot_r_old for next iteration.
            r_dot_r_old = r_dot_r_new

            write(msg,'(A,I3,2(A,E9.2))') 'CG step ', i, ': res= ', residual, ', tol= ', tol
            call log_information(msg, this_module, this_procedure)
        enddo cg_loop
        end associate

        ! Set and copy info flag for completeness
        info = cg_meta%n_iter
        cg_meta%info = info

        if (opts%if_print_metadata) call cg_meta%print()

        ! Set metadata output
        if (present(meta)) then
            select type(meta)
            type is (cg_dp_metadata)
                meta = cg_meta
            class default
                call type_error('meta','cg_dp_metadata','OUT',this_module,this_procedure)
            end select
        end if

        call A%reset_counter(.false., 'cg%post')
        if (time_lightkrylov()) call timer%stop(this_procedure)
        call log_debug('end', this_module, this_procedure)

        return
    end procedure

end submodule

