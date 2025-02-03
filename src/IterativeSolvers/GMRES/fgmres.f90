submodule (lightkrylov_iterativesolvers) fgmres_solver
    use stdlib_strings, only: padr
    use stdlib_linalg, only: lstsq, norm
    implicit none
contains

    !----------------------------------------
    !-----     Options and Metadata     -----
    !----------------------------------------

    module procedure print_fgmres_sp
        ! internals
        integer :: i
        logical :: ifreset, ifverbose
        character(len=128) :: msg

        ifreset   = optval(reset_counters, .false.)
        ifverbose = optval(verbose, .false.)
  
        write(msg,'(A30,I6,"  (",I6,"/",I3,")")') padr('Iterations   (inner/outer): ', 30), &
                  & self%n_iter, self%n_inner, self%n_outer
        call log_message(msg, module=this_module, procedure='fgmres_metadata')
        if (ifverbose) then
            write(msg,'(14X,A15)') 'Residual'
            call log_message(msg, module=this_module, procedure='fgmres_metadata')
            write(msg,'(A14,E15.8)') '   INIT:', self%res(1)
            call log_message(msg, module=this_module, procedure='fgmres_metadata')
            do i = 1, self%n_iter
               write(msg,'(A,I3,A,E20.8)') '   Step ', i, ': ', self%res(i)
               call log_message(msg, module=this_module, procedure='fgmres_metadata')
            end do
        else
            write(msg,'(A30,I20)') padr('Number of records: ', 30), size(self%res)
            call log_message(msg, module=this_module, procedure='fgmres_metadata')
            write(msg,'(A30,E20.8)') padr('Residual: ', 30), self%res(size(self%res))
            call log_message(msg, module=this_module, procedure='fgmres_metadata')
        end if
        if (self%converged) then
            call log_message('Status: CONVERGED', module=this_module, procedure='fgmres_metadata')
        else
            call log_message('Status: NOT CONVERGED', module=this_module, procedure='fgmres_metadata')
        end if
        if (ifreset) call self%reset()
        return
    end procedure

    module procedure reset_fgmres_sp
        self%n_iter = 0
        self%n_inner = 0
        self%n_outer = 0
        self%converged = .false.
        self%info = 0
        if (allocated(self%res)) deallocate(self%res)
        return
    end procedure
    module procedure print_fgmres_dp
        ! internals
        integer :: i
        logical :: ifreset, ifverbose
        character(len=128) :: msg

        ifreset   = optval(reset_counters, .false.)
        ifverbose = optval(verbose, .false.)
  
        write(msg,'(A30,I6,"  (",I6,"/",I3,")")') padr('Iterations   (inner/outer): ', 30), &
                  & self%n_iter, self%n_inner, self%n_outer
        call log_message(msg, module=this_module, procedure='fgmres_metadata')
        if (ifverbose) then
            write(msg,'(14X,A15)') 'Residual'
            call log_message(msg, module=this_module, procedure='fgmres_metadata')
            write(msg,'(A14,E15.8)') '   INIT:', self%res(1)
            call log_message(msg, module=this_module, procedure='fgmres_metadata')
            do i = 1, self%n_iter
               write(msg,'(A,I3,A,E20.8)') '   Step ', i, ': ', self%res(i)
               call log_message(msg, module=this_module, procedure='fgmres_metadata')
            end do
        else
            write(msg,'(A30,I20)') padr('Number of records: ', 30), size(self%res)
            call log_message(msg, module=this_module, procedure='fgmres_metadata')
            write(msg,'(A30,E20.8)') padr('Residual: ', 30), self%res(size(self%res))
            call log_message(msg, module=this_module, procedure='fgmres_metadata')
        end if
        if (self%converged) then
            call log_message('Status: CONVERGED', module=this_module, procedure='fgmres_metadata')
        else
            call log_message('Status: NOT CONVERGED', module=this_module, procedure='fgmres_metadata')
        end if
        if (ifreset) call self%reset()
        return
    end procedure

    module procedure reset_fgmres_dp
        self%n_iter = 0
        self%n_inner = 0
        self%n_outer = 0
        self%converged = .false.
        self%info = 0
        if (allocated(self%res)) deallocate(self%res)
        return
    end procedure

    !-------------------------------------------------------------
    !-----     FLEXIBLE GMRES SOLVERS FOR ABSTRACT TYPES     -----
    !-------------------------------------------------------------

    module procedure fgmres_rsp
        ! Options.
        integer :: kdim, maxiter
        real(sp) :: tol, rtol_, atol_
        logical :: trans
        type(fgmres_sp_opts)     :: opts
        type(fgmres_sp_metadata) :: fgmres_meta

        ! Krylov subspace
        class(abstract_vector_rsp), allocatable :: V(:)
        class(abstract_vector_rsp), allocatable :: Z(:)
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
        class(abstract_vector_rsp), allocatable :: dx
        character(len=256) :: msg

        call log_debug('start', module=this_module, procedure='fgmres_rsp')
        if (time_lightkrylov()) call timer%start('fgmres_rsp')
        ! Deals with the optional args.
        rtol_ = optval(rtol, rtol_sp)
        atol_ = optval(atol, atol_sp)
        if (present(options)) then
            select type (options)
            type is (fgmres_sp_opts)
                opts = options
            end select
        else
            opts = fgmres_sp_opts()
        endif

        kdim = opts%kdim ; maxiter = opts%maxiter
        tol = atol_ + rtol_ * b%norm()
        trans = optval(transpose, .false.)

        ! Deals with the preconditioner.
        has_precond = optval(present(preconditioner), .false.)
        if (has_precond) allocate(precond, source=preconditioner)

        ! Initialize working variables.
        allocate(V(kdim+1), mold=b) ; call zero_basis(V)
        allocate(Z(kdim+1), mold=b) ; call zero_basis(Z)
        allocate(H(kdim+1, kdim)) ; H = 0.0_sp
        allocate(y(kdim)) ; y = 0.0_sp
        allocate(e(kdim+1)) ; e = 0.0_sp

        ! Initialize metadata and & reset matvec counter
        fgmres_meta = fgmres_sp_metadata()
        call A%reset_counter(trans, 'fgmres%init')

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
        allocate(fgmres_meta%res(1)); fgmres_meta%res(1) = abs(beta)

        write(msg,'(A,I3,2(A,E9.2))') 'FGMRES(k)   inner step ', 0, ': |res|= ', &
                    & abs(beta), ', tol= ', tol
        call log_information(msg, module=this_module, procedure='fgmres_rsp')

        ! Iterative solver.
        fgmres_iter : do i = 1, maxiter
            ! Zero-out variables.
            H = 0.0_sp ; y = 0.0_sp ; e = 0.0_sp ; e(1) = beta

            ! Arnoldi factorization.
            arnoldi_fact: do k = 1, kdim
                ! Preconditioner.
                call copy(Z(k), V(k)); if (has_precond) call precond%apply(Z(k), k, beta, tol)

                ! Matrix-vector product.
                if (trans) then
                    call A%apply_rmatvec(Z(k), V(k+1))
                else
                    call A%apply_matvec(Z(k), V(k+1))
                endif

                ! Double Gram-Schmid orthogonalization
                call double_gram_schmidt_step(V(k+1), V(:k), info, if_chk_orthonormal=.false., beta=H(:k, k))
                call check_info(info, 'double_gram_schmidt_step', module=this_module, procedure='fgmres_rsp')

                ! Update Hessenberg matrix and normalize residual Krylov vector.
                H(k+1, k) = V(k+1)%norm()
                if (abs(H(k+1, k)) > tol) call V(k+1)%scal(one_rsp / H(k+1, k))

                ! Least-squares problem.
                y(:k) = lstsq(H(:k+1, :k), e(:k+1))

                ! Compute residual.
                beta = norm(e(:k+1) - matmul(H(:k+1, :k), y(:k)), 2)

                ! Save metadata.
                fgmres_meta%n_iter  = fgmres_meta%n_iter + 1
                fgmres_meta%n_inner = fgmres_meta%n_inner + 1
                fgmres_meta%res = [ fgmres_meta%res, abs(beta) ]

                ! Check convergence.
                write(msg,'(A,I3,2(A,E9.2))') 'FGMRES(k)   inner step ', k, ': |res|= ', &
                            & abs(beta), ', tol= ', tol
                call log_information(msg, module=this_module, procedure='fgmres_rsp')
                if (abs(beta) <= tol) then
                    fgmres_meta%converged = .true.
                    exit arnoldi_fact
                endif
            enddo arnoldi_fact

            ! Update solution.
            k = min(k, kdim) ; call linear_combination(dx, Z(:k), y(:k)); call x%add(dx)

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
            fgmres_meta%n_iter  = fgmres_meta%n_iter + 1
            fgmres_meta%n_outer = fgmres_meta%n_outer + 1
            fgmres_meta%res = [ fgmres_meta%res, abs(beta) ]

            write(msg,'(A,I3,2(A,E9.2))') 'FGMRES(k) outer step   ', i, ': |res|= ', &
                            & abs(beta), ', tol= ', tol
            call log_information(msg, module=this_module, procedure='fgmres_rsp')

            ! Exit gmres if desired accuracy is reached.
            if (abs(beta) <= tol) then
               fgmres_meta%converged = .true.
               exit fgmres_iter
            end if
        enddo fgmres_iter

        ! Returns the number of iterations.
        info = fgmres_meta%n_iter
        fgmres_meta%info = info

        if (opts%if_print_metadata) call fgmres_meta%print()

        ! Set metadata output
        if (present(meta)) then
           select type(meta)
                 type is (fgmres_sp_metadata)
                    meta = fgmres_meta
           end select
        end if

        call A%reset_counter(trans, 'fgmres%post')
        if (time_lightkrylov()) call timer%stop('fgmres_rsp')
        call log_debug('end', module=this_module, procedure='fgmres_rsp')
        
        return
    end procedure
    module procedure fgmres_rdp
        ! Options.
        integer :: kdim, maxiter
        real(dp) :: tol, rtol_, atol_
        logical :: trans
        type(fgmres_dp_opts)     :: opts
        type(fgmres_dp_metadata) :: fgmres_meta

        ! Krylov subspace
        class(abstract_vector_rdp), allocatable :: V(:)
        class(abstract_vector_rdp), allocatable :: Z(:)
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
        class(abstract_vector_rdp), allocatable :: dx
        character(len=256) :: msg

        call log_debug('start', module=this_module, procedure='fgmres_rdp')
        if (time_lightkrylov()) call timer%start('fgmres_rdp')
        ! Deals with the optional args.
        rtol_ = optval(rtol, rtol_dp)
        atol_ = optval(atol, atol_dp)
        if (present(options)) then
            select type (options)
            type is (fgmres_dp_opts)
                opts = options
            end select
        else
            opts = fgmres_dp_opts()
        endif

        kdim = opts%kdim ; maxiter = opts%maxiter
        tol = atol_ + rtol_ * b%norm()
        trans = optval(transpose, .false.)

        ! Deals with the preconditioner.
        has_precond = optval(present(preconditioner), .false.)
        if (has_precond) allocate(precond, source=preconditioner)

        ! Initialize working variables.
        allocate(V(kdim+1), mold=b) ; call zero_basis(V)
        allocate(Z(kdim+1), mold=b) ; call zero_basis(Z)
        allocate(H(kdim+1, kdim)) ; H = 0.0_dp
        allocate(y(kdim)) ; y = 0.0_dp
        allocate(e(kdim+1)) ; e = 0.0_dp

        ! Initialize metadata and & reset matvec counter
        fgmres_meta = fgmres_dp_metadata()
        call A%reset_counter(trans, 'fgmres%init')

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
        allocate(fgmres_meta%res(1)); fgmres_meta%res(1) = abs(beta)

        write(msg,'(A,I3,2(A,E9.2))') 'FGMRES(k)   inner step ', 0, ': |res|= ', &
                    & abs(beta), ', tol= ', tol
        call log_information(msg, module=this_module, procedure='fgmres_rdp')

        ! Iterative solver.
        fgmres_iter : do i = 1, maxiter
            ! Zero-out variables.
            H = 0.0_dp ; y = 0.0_dp ; e = 0.0_dp ; e(1) = beta

            ! Arnoldi factorization.
            arnoldi_fact: do k = 1, kdim
                ! Preconditioner.
                call copy(Z(k), V(k)); if (has_precond) call precond%apply(Z(k), k, beta, tol)

                ! Matrix-vector product.
                if (trans) then
                    call A%apply_rmatvec(Z(k), V(k+1))
                else
                    call A%apply_matvec(Z(k), V(k+1))
                endif

                ! Double Gram-Schmid orthogonalization
                call double_gram_schmidt_step(V(k+1), V(:k), info, if_chk_orthonormal=.false., beta=H(:k, k))
                call check_info(info, 'double_gram_schmidt_step', module=this_module, procedure='fgmres_rdp')

                ! Update Hessenberg matrix and normalize residual Krylov vector.
                H(k+1, k) = V(k+1)%norm()
                if (abs(H(k+1, k)) > tol) call V(k+1)%scal(one_rdp / H(k+1, k))

                ! Least-squares problem.
                y(:k) = lstsq(H(:k+1, :k), e(:k+1))

                ! Compute residual.
                beta = norm(e(:k+1) - matmul(H(:k+1, :k), y(:k)), 2)

                ! Save metadata.
                fgmres_meta%n_iter  = fgmres_meta%n_iter + 1
                fgmres_meta%n_inner = fgmres_meta%n_inner + 1
                fgmres_meta%res = [ fgmres_meta%res, abs(beta) ]

                ! Check convergence.
                write(msg,'(A,I3,2(A,E9.2))') 'FGMRES(k)   inner step ', k, ': |res|= ', &
                            & abs(beta), ', tol= ', tol
                call log_information(msg, module=this_module, procedure='fgmres_rdp')
                if (abs(beta) <= tol) then
                    fgmres_meta%converged = .true.
                    exit arnoldi_fact
                endif
            enddo arnoldi_fact

            ! Update solution.
            k = min(k, kdim) ; call linear_combination(dx, Z(:k), y(:k)); call x%add(dx)

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
            fgmres_meta%n_iter  = fgmres_meta%n_iter + 1
            fgmres_meta%n_outer = fgmres_meta%n_outer + 1
            fgmres_meta%res = [ fgmres_meta%res, abs(beta) ]

            write(msg,'(A,I3,2(A,E9.2))') 'FGMRES(k) outer step   ', i, ': |res|= ', &
                            & abs(beta), ', tol= ', tol
            call log_information(msg, module=this_module, procedure='fgmres_rdp')

            ! Exit gmres if desired accuracy is reached.
            if (abs(beta) <= tol) then
               fgmres_meta%converged = .true.
               exit fgmres_iter
            end if
        enddo fgmres_iter

        ! Returns the number of iterations.
        info = fgmres_meta%n_iter
        fgmres_meta%info = info

        if (opts%if_print_metadata) call fgmres_meta%print()

        ! Set metadata output
        if (present(meta)) then
           select type(meta)
                 type is (fgmres_dp_metadata)
                    meta = fgmres_meta
           end select
        end if

        call A%reset_counter(trans, 'fgmres%post')
        if (time_lightkrylov()) call timer%stop('fgmres_rdp')
        call log_debug('end', module=this_module, procedure='fgmres_rdp')
        
        return
    end procedure
    module procedure fgmres_csp
        ! Options.
        integer :: kdim, maxiter
        real(sp) :: tol, rtol_, atol_
        logical :: trans
        type(fgmres_sp_opts)     :: opts
        type(fgmres_sp_metadata) :: fgmres_meta

        ! Krylov subspace
        class(abstract_vector_csp), allocatable :: V(:)
        class(abstract_vector_csp), allocatable :: Z(:)
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
        class(abstract_vector_csp), allocatable :: dx
        character(len=256) :: msg

        call log_debug('start', module=this_module, procedure='fgmres_csp')
        if (time_lightkrylov()) call timer%start('fgmres_csp')
        ! Deals with the optional args.
        rtol_ = optval(rtol, rtol_sp)
        atol_ = optval(atol, atol_sp)
        if (present(options)) then
            select type (options)
            type is (fgmres_sp_opts)
                opts = options
            end select
        else
            opts = fgmres_sp_opts()
        endif

        kdim = opts%kdim ; maxiter = opts%maxiter
        tol = atol_ + rtol_ * b%norm()
        trans = optval(transpose, .false.)

        ! Deals with the preconditioner.
        has_precond = optval(present(preconditioner), .false.)
        if (has_precond) allocate(precond, source=preconditioner)

        ! Initialize working variables.
        allocate(V(kdim+1), mold=b) ; call zero_basis(V)
        allocate(Z(kdim+1), mold=b) ; call zero_basis(Z)
        allocate(H(kdim+1, kdim)) ; H = 0.0_sp
        allocate(y(kdim)) ; y = 0.0_sp
        allocate(e(kdim+1)) ; e = 0.0_sp

        ! Initialize metadata and & reset matvec counter
        fgmres_meta = fgmres_sp_metadata()
        call A%reset_counter(trans, 'fgmres%init')

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
        allocate(fgmres_meta%res(1)); fgmres_meta%res(1) = abs(beta)

        write(msg,'(A,I3,2(A,E9.2))') 'FGMRES(k)   inner step ', 0, ': |res|= ', &
                    & abs(beta), ', tol= ', tol
        call log_information(msg, module=this_module, procedure='fgmres_csp')

        ! Iterative solver.
        fgmres_iter : do i = 1, maxiter
            ! Zero-out variables.
            H = 0.0_sp ; y = 0.0_sp ; e = 0.0_sp ; e(1) = beta

            ! Arnoldi factorization.
            arnoldi_fact: do k = 1, kdim
                ! Preconditioner.
                call copy(Z(k), V(k)); if (has_precond) call precond%apply(Z(k), k, beta, tol)

                ! Matrix-vector product.
                if (trans) then
                    call A%apply_rmatvec(Z(k), V(k+1))
                else
                    call A%apply_matvec(Z(k), V(k+1))
                endif

                ! Double Gram-Schmid orthogonalization
                call double_gram_schmidt_step(V(k+1), V(:k), info, if_chk_orthonormal=.false., beta=H(:k, k))
                call check_info(info, 'double_gram_schmidt_step', module=this_module, procedure='fgmres_csp')

                ! Update Hessenberg matrix and normalize residual Krylov vector.
                H(k+1, k) = V(k+1)%norm()
                if (abs(H(k+1, k)) > tol) call V(k+1)%scal(one_csp / H(k+1, k))

                ! Least-squares problem.
                y(:k) = lstsq(H(:k+1, :k), e(:k+1))

                ! Compute residual.
                beta = norm(e(:k+1) - matmul(H(:k+1, :k), y(:k)), 2)

                ! Save metadata.
                fgmres_meta%n_iter  = fgmres_meta%n_iter + 1
                fgmres_meta%n_inner = fgmres_meta%n_inner + 1
                fgmres_meta%res = [ fgmres_meta%res, abs(beta) ]

                ! Check convergence.
                write(msg,'(A,I3,2(A,E9.2))') 'FGMRES(k)   inner step ', k, ': |res|= ', &
                            & abs(beta), ', tol= ', tol
                call log_information(msg, module=this_module, procedure='fgmres_csp')
                if (abs(beta) <= tol) then
                    fgmres_meta%converged = .true.
                    exit arnoldi_fact
                endif
            enddo arnoldi_fact

            ! Update solution.
            k = min(k, kdim) ; call linear_combination(dx, Z(:k), y(:k)); call x%add(dx)

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
            fgmres_meta%n_iter  = fgmres_meta%n_iter + 1
            fgmres_meta%n_outer = fgmres_meta%n_outer + 1
            fgmres_meta%res = [ fgmres_meta%res, abs(beta) ]

            write(msg,'(A,I3,2(A,E9.2))') 'FGMRES(k) outer step   ', i, ': |res|= ', &
                            & abs(beta), ', tol= ', tol
            call log_information(msg, module=this_module, procedure='fgmres_csp')

            ! Exit gmres if desired accuracy is reached.
            if (abs(beta) <= tol) then
               fgmres_meta%converged = .true.
               exit fgmres_iter
            end if
        enddo fgmres_iter

        ! Returns the number of iterations.
        info = fgmres_meta%n_iter
        fgmres_meta%info = info

        if (opts%if_print_metadata) call fgmres_meta%print()

        ! Set metadata output
        if (present(meta)) then
           select type(meta)
                 type is (fgmres_sp_metadata)
                    meta = fgmres_meta
           end select
        end if

        call A%reset_counter(trans, 'fgmres%post')
        if (time_lightkrylov()) call timer%stop('fgmres_csp')
        call log_debug('end', module=this_module, procedure='fgmres_csp')
        
        return
    end procedure
    module procedure fgmres_cdp
        ! Options.
        integer :: kdim, maxiter
        real(dp) :: tol, rtol_, atol_
        logical :: trans
        type(fgmres_dp_opts)     :: opts
        type(fgmres_dp_metadata) :: fgmres_meta

        ! Krylov subspace
        class(abstract_vector_cdp), allocatable :: V(:)
        class(abstract_vector_cdp), allocatable :: Z(:)
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
        class(abstract_vector_cdp), allocatable :: dx
        character(len=256) :: msg

        call log_debug('start', module=this_module, procedure='fgmres_cdp')
        if (time_lightkrylov()) call timer%start('fgmres_cdp')
        ! Deals with the optional args.
        rtol_ = optval(rtol, rtol_dp)
        atol_ = optval(atol, atol_dp)
        if (present(options)) then
            select type (options)
            type is (fgmres_dp_opts)
                opts = options
            end select
        else
            opts = fgmres_dp_opts()
        endif

        kdim = opts%kdim ; maxiter = opts%maxiter
        tol = atol_ + rtol_ * b%norm()
        trans = optval(transpose, .false.)

        ! Deals with the preconditioner.
        has_precond = optval(present(preconditioner), .false.)
        if (has_precond) allocate(precond, source=preconditioner)

        ! Initialize working variables.
        allocate(V(kdim+1), mold=b) ; call zero_basis(V)
        allocate(Z(kdim+1), mold=b) ; call zero_basis(Z)
        allocate(H(kdim+1, kdim)) ; H = 0.0_dp
        allocate(y(kdim)) ; y = 0.0_dp
        allocate(e(kdim+1)) ; e = 0.0_dp

        ! Initialize metadata and & reset matvec counter
        fgmres_meta = fgmres_dp_metadata()
        call A%reset_counter(trans, 'fgmres%init')

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
        allocate(fgmres_meta%res(1)); fgmres_meta%res(1) = abs(beta)

        write(msg,'(A,I3,2(A,E9.2))') 'FGMRES(k)   inner step ', 0, ': |res|= ', &
                    & abs(beta), ', tol= ', tol
        call log_information(msg, module=this_module, procedure='fgmres_cdp')

        ! Iterative solver.
        fgmres_iter : do i = 1, maxiter
            ! Zero-out variables.
            H = 0.0_dp ; y = 0.0_dp ; e = 0.0_dp ; e(1) = beta

            ! Arnoldi factorization.
            arnoldi_fact: do k = 1, kdim
                ! Preconditioner.
                call copy(Z(k), V(k)); if (has_precond) call precond%apply(Z(k), k, beta, tol)

                ! Matrix-vector product.
                if (trans) then
                    call A%apply_rmatvec(Z(k), V(k+1))
                else
                    call A%apply_matvec(Z(k), V(k+1))
                endif

                ! Double Gram-Schmid orthogonalization
                call double_gram_schmidt_step(V(k+1), V(:k), info, if_chk_orthonormal=.false., beta=H(:k, k))
                call check_info(info, 'double_gram_schmidt_step', module=this_module, procedure='fgmres_cdp')

                ! Update Hessenberg matrix and normalize residual Krylov vector.
                H(k+1, k) = V(k+1)%norm()
                if (abs(H(k+1, k)) > tol) call V(k+1)%scal(one_cdp / H(k+1, k))

                ! Least-squares problem.
                y(:k) = lstsq(H(:k+1, :k), e(:k+1))

                ! Compute residual.
                beta = norm(e(:k+1) - matmul(H(:k+1, :k), y(:k)), 2)

                ! Save metadata.
                fgmres_meta%n_iter  = fgmres_meta%n_iter + 1
                fgmres_meta%n_inner = fgmres_meta%n_inner + 1
                fgmres_meta%res = [ fgmres_meta%res, abs(beta) ]

                ! Check convergence.
                write(msg,'(A,I3,2(A,E9.2))') 'FGMRES(k)   inner step ', k, ': |res|= ', &
                            & abs(beta), ', tol= ', tol
                call log_information(msg, module=this_module, procedure='fgmres_cdp')
                if (abs(beta) <= tol) then
                    fgmres_meta%converged = .true.
                    exit arnoldi_fact
                endif
            enddo arnoldi_fact

            ! Update solution.
            k = min(k, kdim) ; call linear_combination(dx, Z(:k), y(:k)); call x%add(dx)

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
            fgmres_meta%n_iter  = fgmres_meta%n_iter + 1
            fgmres_meta%n_outer = fgmres_meta%n_outer + 1
            fgmres_meta%res = [ fgmres_meta%res, abs(beta) ]

            write(msg,'(A,I3,2(A,E9.2))') 'FGMRES(k) outer step   ', i, ': |res|= ', &
                            & abs(beta), ', tol= ', tol
            call log_information(msg, module=this_module, procedure='fgmres_cdp')

            ! Exit gmres if desired accuracy is reached.
            if (abs(beta) <= tol) then
               fgmres_meta%converged = .true.
               exit fgmres_iter
            end if
        enddo fgmres_iter

        ! Returns the number of iterations.
        info = fgmres_meta%n_iter
        fgmres_meta%info = info

        if (opts%if_print_metadata) call fgmres_meta%print()

        ! Set metadata output
        if (present(meta)) then
           select type(meta)
                 type is (fgmres_dp_metadata)
                    meta = fgmres_meta
           end select
        end if

        call A%reset_counter(trans, 'fgmres%post')
        if (time_lightkrylov()) call timer%stop('fgmres_cdp')
        call log_debug('end', module=this_module, procedure='fgmres_cdp')
        
        return
    end procedure

end submodule
