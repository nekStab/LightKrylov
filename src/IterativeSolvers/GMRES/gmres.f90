submodule (lightkrylov_iterativesolvers) gmres_solver
    use stdlib_strings, only: padr
    use stdlib_linalg, only: lstsq, norm
    implicit none
contains

    !----------------------------------------
    !-----     Options and Metadata     -----
    !----------------------------------------

    module procedure print_gmres_sp
        ! internals
        integer :: i
        logical :: ifreset, ifverbose
        character(len=128) :: msg

        ifreset   = optval(reset_counters, .false.)
        ifverbose = optval(verbose, .false.)
  
        write(msg,'(A30,I6,"  (",I6,"/",I3,")")') padr('Iterations   (inner/outer): ', 30), &
                  & self%n_iter, self%n_inner, self%n_outer
        call log_message(msg, module=this_module, procedure='gmres_metadata')
        if (ifverbose) then
            write(msg,'(14X,A15)') 'Residual'
            call log_message(msg, module=this_module, procedure='gmres_metadata')
            write(msg,'(A14,E15.8)') '   INIT:', self%res(1)
            call log_message(msg, module=this_module, procedure='gmres_metadata')
            do i = 1, self%n_iter
               write(msg,'(A,I3,A,E20.8)') '   Step ', i, ': ', self%res(i)
               call log_message(msg, module=this_module, procedure='gmres_metadata')
            end do
        else
            write(msg,'(A30,I20)') padr('Number of records: ', 30), size(self%res)
            call log_message(msg, module=this_module, procedure='gmres_metadata')
            write(msg,'(A30,E20.8)') padr('Residual: ', 30), self%res(size(self%res))
            call log_message(msg, module=this_module, procedure='gmres_metadata')
        end if
        if (self%converged) then
            call log_message('Status: CONVERGED', module=this_module, procedure='gmres_metadata')
        else
            call log_message('Status: NOT CONVERGED', module=this_module, procedure='gmres_metadata')
        end if
        if (ifreset) call self%reset()
        return
    end procedure

    module procedure reset_gmres_sp
        self%n_iter = 0
        self%n_inner = 0
        self%n_outer = 0
        self%converged = .false.
        self%info = 0
        if (allocated(self%res)) deallocate(self%res)
        return
    end procedure
    module procedure print_gmres_dp
        ! internals
        integer :: i
        logical :: ifreset, ifverbose
        character(len=128) :: msg

        ifreset   = optval(reset_counters, .false.)
        ifverbose = optval(verbose, .false.)
  
        write(msg,'(A30,I6,"  (",I6,"/",I3,")")') padr('Iterations   (inner/outer): ', 30), &
                  & self%n_iter, self%n_inner, self%n_outer
        call log_message(msg, module=this_module, procedure='gmres_metadata')
        if (ifverbose) then
            write(msg,'(14X,A15)') 'Residual'
            call log_message(msg, module=this_module, procedure='gmres_metadata')
            write(msg,'(A14,E15.8)') '   INIT:', self%res(1)
            call log_message(msg, module=this_module, procedure='gmres_metadata')
            do i = 1, self%n_iter
               write(msg,'(A,I3,A,E20.8)') '   Step ', i, ': ', self%res(i)
               call log_message(msg, module=this_module, procedure='gmres_metadata')
            end do
        else
            write(msg,'(A30,I20)') padr('Number of records: ', 30), size(self%res)
            call log_message(msg, module=this_module, procedure='gmres_metadata')
            write(msg,'(A30,E20.8)') padr('Residual: ', 30), self%res(size(self%res))
            call log_message(msg, module=this_module, procedure='gmres_metadata')
        end if
        if (self%converged) then
            call log_message('Status: CONVERGED', module=this_module, procedure='gmres_metadata')
        else
            call log_message('Status: NOT CONVERGED', module=this_module, procedure='gmres_metadata')
        end if
        if (ifreset) call self%reset()
        return
    end procedure

    module procedure reset_gmres_dp
        self%n_iter = 0
        self%n_inner = 0
        self%n_outer = 0
        self%converged = .false.
        self%info = 0
        if (allocated(self%res)) deallocate(self%res)
        return
    end procedure

    !----------------------------------------------------
    !-----     GMRES SOLVERS FOR ABSTRACT TYPES     -----
    !----------------------------------------------------

    module procedure gmres_rsp
       ! Options.
        integer :: kdim, maxiter
        real(sp) :: tol, rtol_, atol_
        logical :: trans
        type(gmres_sp_opts)     :: opts
        type(gmres_sp_metadata) :: gmres_meta

        ! Krylov subspace
        class(abstract_vector_rsp), allocatable :: V(:)
        ! Hessenberg matrix.
        real(sp), allocatable :: H(:, :)
        ! Least-squares variables.
        real(sp), allocatable :: y(:), e(:)
        real(sp) :: beta

        ! Preconditioner
        logical :: has_precond
        class(abstract_precond_rsp), allocatable :: precond

        ! Miscellaneous.
        integer :: i, k
        class(abstract_vector_rsp), allocatable :: dx, wrk
        character(len=256) :: msg

        if (time_lightkrylov()) call timer%start('gmres_rsp')
        ! Deals with the optional args.
        rtol_ = optval(rtol, rtol_sp)
        atol_ = optval(atol, atol_sp)
        if (present(options)) then
            select type (options)
            type is (gmres_sp_opts)
                opts = options
            class default
                call stop_error("The optional intent [IN] argument 'options' must be of type 'gmres_sp_opts'", module=this_module,&
                    & procedure='gmres_rsp')
            end select
        else
            opts = gmres_sp_opts()
        endif

        kdim = opts%kdim ; maxiter = opts%maxiter
        tol = atol_ + rtol_ * b%norm()
        trans = optval(transpose, .false.)

        ! Deals with the preconditioner.
        has_precond = optval(present(preconditioner), .false.)
        if (has_precond) allocate(precond, source=preconditioner)

        ! Initialize working variables.
        allocate(wrk, mold=b) ; call wrk%zero()
        allocate(V(kdim+1), mold=b) ; call zero_basis(V)
        allocate(H(kdim+1, kdim)) ; H = 0.0_sp
        allocate(y(kdim)) ; y = 0.0_sp
        allocate(e(kdim+1)) ; e = 0.0_sp

        ! Initialize metadata and & reset matvec counter
        gmres_meta = gmres_sp_metadata()
        call A%reset_counter(trans, 'gmres%init')

        info = 0

        ! Initial Krylov vector.
        if (x%norm() > 0) then
            if (trans) then
                call A%apply_rmatvec(x, V(1))
            else
                call A%apply_matvec(x, V(1))
            endif
        endif

        call V(1)%sub(b) ; call V(1)%chsgn()
        beta = V(1)%norm() ; call V(1)%scal(one_rsp/beta)
        allocate(gmres_meta%res(1)); gmres_meta%res(1) = abs(beta)

        write(msg,'(2(A,E11.4))') 'GMRES(k)   init step     : |res|= ', &
                    & abs(beta), ', tol= ', tol
        call log_information(msg, module=this_module, procedure='gmres_rsp')

        ! Iterative solver.
        gmres_iter : do i = 1, maxiter
            ! Zero-out variables.
            H = 0.0_sp ; y = 0.0_sp ; e = 0.0_sp ; e(1) = beta

            ! Arnoldi factorization.
            arnoldi_fact: do k = 1, kdim
                ! Preconditioner.
                wrk = V(k) ; if (has_precond) call precond%apply(wrk, k, beta, tol)

                ! Matrix-vector product.
                if (trans) then
                    call A%apply_rmatvec(wrk, V(k+1))
                else
                    call A%apply_matvec(wrk, V(k+1))
                endif

                ! Double Gram-Schmid orthogonalization
                call double_gram_schmidt_step(V(k+1), V(:k), info, if_chk_orthonormal=.false., beta=H(:k, k))
                call check_info(info, 'double_gram_schmidt_step', module=this_module, procedure='gmres_rsp')

                ! Update Hessenberg matrix and normalize residual Krylov vector.
                H(k+1, k) = V(k+1)%norm()
                if (abs(H(k+1, k)) > tol) call V(k+1)%scal(one_rsp / H(k+1, k))

                ! Least-squares problem.
                y(:k) = lstsq(H(:k+1, :k), e(:k+1))

                ! Compute residual.
                beta = norm(e(:k+1) - matmul(H(:k+1, :k), y(:k)), 2)

                ! Save metadata.
                gmres_meta%n_iter  = gmres_meta%n_iter + 1
                gmres_meta%n_inner = gmres_meta%n_inner + 1
                gmres_meta%res = [ gmres_meta%res, abs(beta) ]

                ! Check convergence.
                write(msg,'(A,I3,2(A,E11.4))') 'GMRES(k)   inner step ', k, ': |res|= ', &
                            & abs(beta), ', tol= ', tol
                call log_information(msg, module=this_module, procedure='gmres_rsp')
                if (abs(beta) <= tol) then
                    gmres_meta%converged = .true.
                    exit arnoldi_fact
                endif
            enddo arnoldi_fact

            ! Update solution.
            k = min(k, kdim) ; call linear_combination(dx, V(:k), y(:k))
            if (has_precond) call precond%apply(dx) ; call x%add(dx)

            ! Recompute residual for sanity check.
            if (trans) then
                call A%apply_rmatvec(x, v(1))
            else
                call A%apply_matvec(x, v(1))
            endif
            call v(1)%sub(b) ; call v(1)%chsgn()

            ! Initialize new starting Krylov vector if needed.
            beta = v(1)%norm() ; call v(1)%scal(one_rsp / beta)

            ! Save metadata.
            gmres_meta%n_iter  = gmres_meta%n_iter + 1
            gmres_meta%n_outer = gmres_meta%n_outer + 1
            gmres_meta%res = [ gmres_meta%res, abs(beta) ]

            write(msg,'(A,I3,2(A,E11.4))') 'GMRES(k) outer step   ', i, ': |res|= ', &
                            & abs(beta), ', tol= ', tol
            call log_information(msg, module=this_module, procedure='gmres_rsp')

            ! Exit gmres if desired accuracy is reached.
            if (abs(beta) <= tol) then
               gmres_meta%converged = .true.
               exit gmres_iter
            end if
        enddo gmres_iter

        ! Returns the number of iterations.
        info = gmres_meta%n_iter
        gmres_meta%info = info

        if (opts%if_print_metadata) call gmres_meta%print()

        ! Set metadata output
        if (present(meta)) then
            select type(meta)
            type is (gmres_sp_metadata)
                meta = gmres_meta
            class default
                call stop_error("The optional intent [OUT] argument 'meta' must be of type 'gmres_sp_metadata'",&
                    & module=this_module, procedure='gmres_rsp')
            end select
        end if

        call A%reset_counter(trans, 'gmres%post')
        if (time_lightkrylov()) call timer%stop('gmres_rsp')
        
        return
    end procedure 
    module procedure gmres_rdp
       ! Options.
        integer :: kdim, maxiter
        real(dp) :: tol, rtol_, atol_
        logical :: trans
        type(gmres_dp_opts)     :: opts
        type(gmres_dp_metadata) :: gmres_meta

        ! Krylov subspace
        class(abstract_vector_rdp), allocatable :: V(:)
        ! Hessenberg matrix.
        real(dp), allocatable :: H(:, :)
        ! Least-squares variables.
        real(dp), allocatable :: y(:), e(:)
        real(dp) :: beta

        ! Preconditioner
        logical :: has_precond
        class(abstract_precond_rdp), allocatable :: precond

        ! Miscellaneous.
        integer :: i, k
        class(abstract_vector_rdp), allocatable :: dx, wrk
        character(len=256) :: msg

        if (time_lightkrylov()) call timer%start('gmres_rdp')
        ! Deals with the optional args.
        rtol_ = optval(rtol, rtol_dp)
        atol_ = optval(atol, atol_dp)
        if (present(options)) then
            select type (options)
            type is (gmres_dp_opts)
                opts = options
            class default
                call stop_error("The optional intent [IN] argument 'options' must be of type 'gmres_dp_opts'", module=this_module,&
                    & procedure='gmres_rdp')
            end select
        else
            opts = gmres_dp_opts()
        endif

        kdim = opts%kdim ; maxiter = opts%maxiter
        tol = atol_ + rtol_ * b%norm()
        trans = optval(transpose, .false.)

        ! Deals with the preconditioner.
        has_precond = optval(present(preconditioner), .false.)
        if (has_precond) allocate(precond, source=preconditioner)

        ! Initialize working variables.
        allocate(wrk, mold=b) ; call wrk%zero()
        allocate(V(kdim+1), mold=b) ; call zero_basis(V)
        allocate(H(kdim+1, kdim)) ; H = 0.0_dp
        allocate(y(kdim)) ; y = 0.0_dp
        allocate(e(kdim+1)) ; e = 0.0_dp

        ! Initialize metadata and & reset matvec counter
        gmres_meta = gmres_dp_metadata()
        call A%reset_counter(trans, 'gmres%init')

        info = 0

        ! Initial Krylov vector.
        if (x%norm() > 0) then
            if (trans) then
                call A%apply_rmatvec(x, V(1))
            else
                call A%apply_matvec(x, V(1))
            endif
        endif

        call V(1)%sub(b) ; call V(1)%chsgn()
        beta = V(1)%norm() ; call V(1)%scal(one_rdp/beta)
        allocate(gmres_meta%res(1)); gmres_meta%res(1) = abs(beta)

        write(msg,'(2(A,E11.4))') 'GMRES(k)   init step     : |res|= ', &
                    & abs(beta), ', tol= ', tol
        call log_information(msg, module=this_module, procedure='gmres_rdp')

        ! Iterative solver.
        gmres_iter : do i = 1, maxiter
            ! Zero-out variables.
            H = 0.0_dp ; y = 0.0_dp ; e = 0.0_dp ; e(1) = beta

            ! Arnoldi factorization.
            arnoldi_fact: do k = 1, kdim
                ! Preconditioner.
                wrk = V(k) ; if (has_precond) call precond%apply(wrk, k, beta, tol)

                ! Matrix-vector product.
                if (trans) then
                    call A%apply_rmatvec(wrk, V(k+1))
                else
                    call A%apply_matvec(wrk, V(k+1))
                endif

                ! Double Gram-Schmid orthogonalization
                call double_gram_schmidt_step(V(k+1), V(:k), info, if_chk_orthonormal=.false., beta=H(:k, k))
                call check_info(info, 'double_gram_schmidt_step', module=this_module, procedure='gmres_rdp')

                ! Update Hessenberg matrix and normalize residual Krylov vector.
                H(k+1, k) = V(k+1)%norm()
                if (abs(H(k+1, k)) > tol) call V(k+1)%scal(one_rdp / H(k+1, k))

                ! Least-squares problem.
                y(:k) = lstsq(H(:k+1, :k), e(:k+1))

                ! Compute residual.
                beta = norm(e(:k+1) - matmul(H(:k+1, :k), y(:k)), 2)

                ! Save metadata.
                gmres_meta%n_iter  = gmres_meta%n_iter + 1
                gmres_meta%n_inner = gmres_meta%n_inner + 1
                gmres_meta%res = [ gmres_meta%res, abs(beta) ]

                ! Check convergence.
                write(msg,'(A,I3,2(A,E11.4))') 'GMRES(k)   inner step ', k, ': |res|= ', &
                            & abs(beta), ', tol= ', tol
                call log_information(msg, module=this_module, procedure='gmres_rdp')
                if (abs(beta) <= tol) then
                    gmres_meta%converged = .true.
                    exit arnoldi_fact
                endif
            enddo arnoldi_fact

            ! Update solution.
            k = min(k, kdim) ; call linear_combination(dx, V(:k), y(:k))
            if (has_precond) call precond%apply(dx) ; call x%add(dx)

            ! Recompute residual for sanity check.
            if (trans) then
                call A%apply_rmatvec(x, v(1))
            else
                call A%apply_matvec(x, v(1))
            endif
            call v(1)%sub(b) ; call v(1)%chsgn()

            ! Initialize new starting Krylov vector if needed.
            beta = v(1)%norm() ; call v(1)%scal(one_rdp / beta)

            ! Save metadata.
            gmres_meta%n_iter  = gmres_meta%n_iter + 1
            gmres_meta%n_outer = gmres_meta%n_outer + 1
            gmres_meta%res = [ gmres_meta%res, abs(beta) ]

            write(msg,'(A,I3,2(A,E11.4))') 'GMRES(k) outer step   ', i, ': |res|= ', &
                            & abs(beta), ', tol= ', tol
            call log_information(msg, module=this_module, procedure='gmres_rdp')

            ! Exit gmres if desired accuracy is reached.
            if (abs(beta) <= tol) then
               gmres_meta%converged = .true.
               exit gmres_iter
            end if
        enddo gmres_iter

        ! Returns the number of iterations.
        info = gmres_meta%n_iter
        gmres_meta%info = info

        if (opts%if_print_metadata) call gmres_meta%print()

        ! Set metadata output
        if (present(meta)) then
            select type(meta)
            type is (gmres_dp_metadata)
                meta = gmres_meta
            class default
                call stop_error("The optional intent [OUT] argument 'meta' must be of type 'gmres_dp_metadata'",&
                    & module=this_module, procedure='gmres_rdp')
            end select
        end if

        call A%reset_counter(trans, 'gmres%post')
        if (time_lightkrylov()) call timer%stop('gmres_rdp')
        
        return
    end procedure 
    module procedure gmres_csp
       ! Options.
        integer :: kdim, maxiter
        real(sp) :: tol, rtol_, atol_
        logical :: trans
        type(gmres_sp_opts)     :: opts
        type(gmres_sp_metadata) :: gmres_meta

        ! Krylov subspace
        class(abstract_vector_csp), allocatable :: V(:)
        ! Hessenberg matrix.
        complex(sp), allocatable :: H(:, :)
        ! Least-squares variables.
        complex(sp), allocatable :: y(:), e(:)
        real(sp) :: beta

        ! Preconditioner
        logical :: has_precond
        class(abstract_precond_csp), allocatable :: precond

        ! Miscellaneous.
        integer :: i, k
        class(abstract_vector_csp), allocatable :: dx, wrk
        character(len=256) :: msg

        if (time_lightkrylov()) call timer%start('gmres_csp')
        ! Deals with the optional args.
        rtol_ = optval(rtol, rtol_sp)
        atol_ = optval(atol, atol_sp)
        if (present(options)) then
            select type (options)
            type is (gmres_sp_opts)
                opts = options
            class default
                call stop_error("The optional intent [IN] argument 'options' must be of type 'gmres_sp_opts'", module=this_module,&
                    & procedure='gmres_csp')
            end select
        else
            opts = gmres_sp_opts()
        endif

        kdim = opts%kdim ; maxiter = opts%maxiter
        tol = atol_ + rtol_ * b%norm()
        trans = optval(transpose, .false.)

        ! Deals with the preconditioner.
        has_precond = optval(present(preconditioner), .false.)
        if (has_precond) allocate(precond, source=preconditioner)

        ! Initialize working variables.
        allocate(wrk, mold=b) ; call wrk%zero()
        allocate(V(kdim+1), mold=b) ; call zero_basis(V)
        allocate(H(kdim+1, kdim)) ; H = 0.0_sp
        allocate(y(kdim)) ; y = 0.0_sp
        allocate(e(kdim+1)) ; e = 0.0_sp

        ! Initialize metadata and & reset matvec counter
        gmres_meta = gmres_sp_metadata()
        call A%reset_counter(trans, 'gmres%init')

        info = 0

        ! Initial Krylov vector.
        if (x%norm() > 0) then
            if (trans) then
                call A%apply_rmatvec(x, V(1))
            else
                call A%apply_matvec(x, V(1))
            endif
        endif

        call V(1)%sub(b) ; call V(1)%chsgn()
        beta = V(1)%norm() ; call V(1)%scal(one_csp/beta)
        allocate(gmres_meta%res(1)); gmres_meta%res(1) = abs(beta)

        write(msg,'(2(A,E11.4))') 'GMRES(k)   init step     : |res|= ', &
                    & abs(beta), ', tol= ', tol
        call log_information(msg, module=this_module, procedure='gmres_csp')

        ! Iterative solver.
        gmres_iter : do i = 1, maxiter
            ! Zero-out variables.
            H = 0.0_sp ; y = 0.0_sp ; e = 0.0_sp ; e(1) = beta

            ! Arnoldi factorization.
            arnoldi_fact: do k = 1, kdim
                ! Preconditioner.
                wrk = V(k) ; if (has_precond) call precond%apply(wrk, k, beta, tol)

                ! Matrix-vector product.
                if (trans) then
                    call A%apply_rmatvec(wrk, V(k+1))
                else
                    call A%apply_matvec(wrk, V(k+1))
                endif

                ! Double Gram-Schmid orthogonalization
                call double_gram_schmidt_step(V(k+1), V(:k), info, if_chk_orthonormal=.false., beta=H(:k, k))
                call check_info(info, 'double_gram_schmidt_step', module=this_module, procedure='gmres_csp')

                ! Update Hessenberg matrix and normalize residual Krylov vector.
                H(k+1, k) = V(k+1)%norm()
                if (abs(H(k+1, k)) > tol) call V(k+1)%scal(one_csp / H(k+1, k))

                ! Least-squares problem.
                y(:k) = lstsq(H(:k+1, :k), e(:k+1))

                ! Compute residual.
                beta = norm(e(:k+1) - matmul(H(:k+1, :k), y(:k)), 2)

                ! Save metadata.
                gmres_meta%n_iter  = gmres_meta%n_iter + 1
                gmres_meta%n_inner = gmres_meta%n_inner + 1
                gmres_meta%res = [ gmres_meta%res, abs(beta) ]

                ! Check convergence.
                write(msg,'(A,I3,2(A,E11.4))') 'GMRES(k)   inner step ', k, ': |res|= ', &
                            & abs(beta), ', tol= ', tol
                call log_information(msg, module=this_module, procedure='gmres_csp')
                if (abs(beta) <= tol) then
                    gmres_meta%converged = .true.
                    exit arnoldi_fact
                endif
            enddo arnoldi_fact

            ! Update solution.
            k = min(k, kdim) ; call linear_combination(dx, V(:k), y(:k))
            if (has_precond) call precond%apply(dx) ; call x%add(dx)

            ! Recompute residual for sanity check.
            if (trans) then
                call A%apply_rmatvec(x, v(1))
            else
                call A%apply_matvec(x, v(1))
            endif
            call v(1)%sub(b) ; call v(1)%chsgn()

            ! Initialize new starting Krylov vector if needed.
            beta = v(1)%norm() ; call v(1)%scal(one_csp / beta)

            ! Save metadata.
            gmres_meta%n_iter  = gmres_meta%n_iter + 1
            gmres_meta%n_outer = gmres_meta%n_outer + 1
            gmres_meta%res = [ gmres_meta%res, abs(beta) ]

            write(msg,'(A,I3,2(A,E11.4))') 'GMRES(k) outer step   ', i, ': |res|= ', &
                            & abs(beta), ', tol= ', tol
            call log_information(msg, module=this_module, procedure='gmres_csp')

            ! Exit gmres if desired accuracy is reached.
            if (abs(beta) <= tol) then
               gmres_meta%converged = .true.
               exit gmres_iter
            end if
        enddo gmres_iter

        ! Returns the number of iterations.
        info = gmres_meta%n_iter
        gmres_meta%info = info

        if (opts%if_print_metadata) call gmres_meta%print()

        ! Set metadata output
        if (present(meta)) then
            select type(meta)
            type is (gmres_sp_metadata)
                meta = gmres_meta
            class default
                call stop_error("The optional intent [OUT] argument 'meta' must be of type 'gmres_sp_metadata'",&
                    & module=this_module, procedure='gmres_csp')
            end select
        end if

        call A%reset_counter(trans, 'gmres%post')
        if (time_lightkrylov()) call timer%stop('gmres_csp')
        
        return
    end procedure 
    module procedure gmres_cdp
       ! Options.
        integer :: kdim, maxiter
        real(dp) :: tol, rtol_, atol_
        logical :: trans
        type(gmres_dp_opts)     :: opts
        type(gmres_dp_metadata) :: gmres_meta

        ! Krylov subspace
        class(abstract_vector_cdp), allocatable :: V(:)
        ! Hessenberg matrix.
        complex(dp), allocatable :: H(:, :)
        ! Least-squares variables.
        complex(dp), allocatable :: y(:), e(:)
        real(dp) :: beta

        ! Preconditioner
        logical :: has_precond
        class(abstract_precond_cdp), allocatable :: precond

        ! Miscellaneous.
        integer :: i, k
        class(abstract_vector_cdp), allocatable :: dx, wrk
        character(len=256) :: msg

        if (time_lightkrylov()) call timer%start('gmres_cdp')
        ! Deals with the optional args.
        rtol_ = optval(rtol, rtol_dp)
        atol_ = optval(atol, atol_dp)
        if (present(options)) then
            select type (options)
            type is (gmres_dp_opts)
                opts = options
            class default
                call stop_error("The optional intent [IN] argument 'options' must be of type 'gmres_dp_opts'", module=this_module,&
                    & procedure='gmres_cdp')
            end select
        else
            opts = gmres_dp_opts()
        endif

        kdim = opts%kdim ; maxiter = opts%maxiter
        tol = atol_ + rtol_ * b%norm()
        trans = optval(transpose, .false.)

        ! Deals with the preconditioner.
        has_precond = optval(present(preconditioner), .false.)
        if (has_precond) allocate(precond, source=preconditioner)

        ! Initialize working variables.
        allocate(wrk, mold=b) ; call wrk%zero()
        allocate(V(kdim+1), mold=b) ; call zero_basis(V)
        allocate(H(kdim+1, kdim)) ; H = 0.0_dp
        allocate(y(kdim)) ; y = 0.0_dp
        allocate(e(kdim+1)) ; e = 0.0_dp

        ! Initialize metadata and & reset matvec counter
        gmres_meta = gmres_dp_metadata()
        call A%reset_counter(trans, 'gmres%init')

        info = 0

        ! Initial Krylov vector.
        if (x%norm() > 0) then
            if (trans) then
                call A%apply_rmatvec(x, V(1))
            else
                call A%apply_matvec(x, V(1))
            endif
        endif

        call V(1)%sub(b) ; call V(1)%chsgn()
        beta = V(1)%norm() ; call V(1)%scal(one_cdp/beta)
        allocate(gmres_meta%res(1)); gmres_meta%res(1) = abs(beta)

        write(msg,'(2(A,E11.4))') 'GMRES(k)   init step     : |res|= ', &
                    & abs(beta), ', tol= ', tol
        call log_information(msg, module=this_module, procedure='gmres_cdp')

        ! Iterative solver.
        gmres_iter : do i = 1, maxiter
            ! Zero-out variables.
            H = 0.0_dp ; y = 0.0_dp ; e = 0.0_dp ; e(1) = beta

            ! Arnoldi factorization.
            arnoldi_fact: do k = 1, kdim
                ! Preconditioner.
                wrk = V(k) ; if (has_precond) call precond%apply(wrk, k, beta, tol)

                ! Matrix-vector product.
                if (trans) then
                    call A%apply_rmatvec(wrk, V(k+1))
                else
                    call A%apply_matvec(wrk, V(k+1))
                endif

                ! Double Gram-Schmid orthogonalization
                call double_gram_schmidt_step(V(k+1), V(:k), info, if_chk_orthonormal=.false., beta=H(:k, k))
                call check_info(info, 'double_gram_schmidt_step', module=this_module, procedure='gmres_cdp')

                ! Update Hessenberg matrix and normalize residual Krylov vector.
                H(k+1, k) = V(k+1)%norm()
                if (abs(H(k+1, k)) > tol) call V(k+1)%scal(one_cdp / H(k+1, k))

                ! Least-squares problem.
                y(:k) = lstsq(H(:k+1, :k), e(:k+1))

                ! Compute residual.
                beta = norm(e(:k+1) - matmul(H(:k+1, :k), y(:k)), 2)

                ! Save metadata.
                gmres_meta%n_iter  = gmres_meta%n_iter + 1
                gmres_meta%n_inner = gmres_meta%n_inner + 1
                gmres_meta%res = [ gmres_meta%res, abs(beta) ]

                ! Check convergence.
                write(msg,'(A,I3,2(A,E11.4))') 'GMRES(k)   inner step ', k, ': |res|= ', &
                            & abs(beta), ', tol= ', tol
                call log_information(msg, module=this_module, procedure='gmres_cdp')
                if (abs(beta) <= tol) then
                    gmres_meta%converged = .true.
                    exit arnoldi_fact
                endif
            enddo arnoldi_fact

            ! Update solution.
            k = min(k, kdim) ; call linear_combination(dx, V(:k), y(:k))
            if (has_precond) call precond%apply(dx) ; call x%add(dx)

            ! Recompute residual for sanity check.
            if (trans) then
                call A%apply_rmatvec(x, v(1))
            else
                call A%apply_matvec(x, v(1))
            endif
            call v(1)%sub(b) ; call v(1)%chsgn()

            ! Initialize new starting Krylov vector if needed.
            beta = v(1)%norm() ; call v(1)%scal(one_cdp / beta)

            ! Save metadata.
            gmres_meta%n_iter  = gmres_meta%n_iter + 1
            gmres_meta%n_outer = gmres_meta%n_outer + 1
            gmres_meta%res = [ gmres_meta%res, abs(beta) ]

            write(msg,'(A,I3,2(A,E11.4))') 'GMRES(k) outer step   ', i, ': |res|= ', &
                            & abs(beta), ', tol= ', tol
            call log_information(msg, module=this_module, procedure='gmres_cdp')

            ! Exit gmres if desired accuracy is reached.
            if (abs(beta) <= tol) then
               gmres_meta%converged = .true.
               exit gmres_iter
            end if
        enddo gmres_iter

        ! Returns the number of iterations.
        info = gmres_meta%n_iter
        gmres_meta%info = info

        if (opts%if_print_metadata) call gmres_meta%print()

        ! Set metadata output
        if (present(meta)) then
            select type(meta)
            type is (gmres_dp_metadata)
                meta = gmres_meta
            class default
                call stop_error("The optional intent [OUT] argument 'meta' must be of type 'gmres_dp_metadata'",&
                    & module=this_module, procedure='gmres_cdp')
            end select
        end if

        call A%reset_counter(trans, 'gmres%post')
        if (time_lightkrylov()) call timer%stop('gmres_cdp')
        
        return
    end procedure 

end submodule
