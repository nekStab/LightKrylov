submodule (lightkrylov_iterativesolvers) fgmres_solver
    use stdlib_strings, only: padr
    use stdlib_linalg, only: lstsq, norm
    use stdlib_linalg_lapack, only: trtrs
    implicit none(type, external)
contains

    !----------------------------------------
    !-----     Options and Metadata     -----
    !----------------------------------------

    module procedure print_fgmres_sp
        ! internals
        character(len=*), parameter :: this_procedure = 'print_fgmres_sp'
        integer :: i
        logical :: ifreset, ifverbose
        character(len=128) :: msg

        ifreset   = optval(reset_counters, .false.)
        ifverbose = optval(verbose, .false.)
  
        write(msg,'(A30,I6,"  (",I6,"/",I3,")")') padr('Iterations   (inner/outer): ', 30), &
                  & self%n_iter, self%n_inner, self%n_outer
        call log_message(msg, this_module, this_procedure)
        if (ifverbose) then
            write(msg,'(14X,A15)') 'Residual'
            call log_message(msg, this_module, this_procedure)
            write(msg,'(A14,E15.8)') '   INIT:', self%res(1)
            call log_message(msg, this_module, this_procedure)
            do i = 1, self%n_iter
               write(msg,'(A,I3,A,E20.8)') '   Step ', i, ': ', self%res(i)
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
    end procedure

    module procedure reset_fgmres_sp
        self%n_iter = 0
        self%n_inner = 0
        self%n_outer = 0
        self%converged = .false.
        self%info = 0
        if (allocated(self%res)) deallocate(self%res)
    end procedure

    module procedure print_fgmres_dp
        ! internals
        character(len=*), parameter :: this_procedure = 'print_fgmres_dp'
        integer :: i
        logical :: ifreset, ifverbose
        character(len=128) :: msg

        ifreset   = optval(reset_counters, .false.)
        ifverbose = optval(verbose, .false.)
  
        write(msg,'(A30,I6,"  (",I6,"/",I3,")")') padr('Iterations   (inner/outer): ', 30), &
                  & self%n_iter, self%n_inner, self%n_outer
        call log_message(msg, this_module, this_procedure)
        if (ifverbose) then
            write(msg,'(14X,A15)') 'Residual'
            call log_message(msg, this_module, this_procedure)
            write(msg,'(A14,E15.8)') '   INIT:', self%res(1)
            call log_message(msg, this_module, this_procedure)
            do i = 1, self%n_iter
               write(msg,'(A,I3,A,E20.8)') '   Step ', i, ': ', self%res(i)
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
    end procedure

    module procedure reset_fgmres_dp
        self%n_iter = 0
        self%n_inner = 0
        self%n_outer = 0
        self%converged = .false.
        self%info = 0
        if (allocated(self%res)) deallocate(self%res)
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
        class(abstract_vector_rsp), allocatable :: V(:), Z(:)
        ! Hessenberg matrix.
        real(sp), allocatable :: H(:, :)
        ! Least-squares variables.
        real(sp), target, allocatable :: e(:)
        real(sp), pointer :: y(:, :)
        real(sp) :: beta
        ! Givens rotations.
        real(sp), allocatable :: c(:), s(:)

        ! Miscellaneous.
        character(len=*), parameter :: this_procedure = 'fgmres_rsp'
        integer :: k, iostat
        class(abstract_vector_rsp), allocatable :: dx, wrk
        character(len=256) :: msg

        if (time_lightkrylov()) call timer%start(this_procedure)
        ! Deals with the optional args.
        rtol_ = optval(rtol, rtol_sp)
        atol_ = optval(atol, atol_sp)
        if (present(options)) then
            select type (options)
            type is (fgmres_sp_opts)
                opts = options
            class default
                call type_error('options','fgmres_sp_opts','IN',this_module,this_procedure)
            end select
        else
            opts = fgmres_sp_opts()
        endif

        kdim  = opts%kdim ; maxiter = opts%maxiter
        tol   = atol_ + rtol_ * b%norm()
        trans = optval(transpose, .false.)

        ! Initialize working variables.
        allocate(wrk, source=b, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)
        call wrk%zero()
        allocate(V(kdim+1), source=b, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)
        call zero_basis(V)
        allocate(Z(kdim+1), source=b, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)
        call zero_basis(Z)
        allocate(H(kdim+1, kdim), source=zero_rsp, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)
        allocate(e(kdim+1), source=zero_rsp, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)
        allocate(c(kdim), source=zero_rsp, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)
        allocate(s(kdim), source=zero_rsp, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)

        ! Initialize metadata and & reset matvec counter
        fgmres_meta = fgmres_sp_metadata() ; fgmres_meta%converged = .false.
        call A%reset_counter(trans, 'fgmres%init')

        info = 0
        associate(ifprecond => present(preconditioner))
        do while ((.not. fgmres_meta%converged) .and. (fgmres_meta%n_outer <= maxiter))
            !> Initialize data
            H = 0.0_sp ; call zero_basis(V)
            if (x%norm() /= 0.0_sp) then
                if (trans) then
                    call A%apply_rmatvec(x, V(1))
                else
                    call A%apply_matvec(x, V(1))
                endif
            endif
            call V(1)%sub(b) ; call V(1)%chsgn()
            e = 0.0_sp ; beta = V(1)%norm() ; e(1) = beta
            call V(1)%scal(one_rsp/beta)
            c = 0.0_sp ; s = 0.0_sp
            if (fgmres_meta%n_outer == 0) then
               allocate(fgmres_meta%res(1), source=abs(beta), stat=iostat, errmsg=msg)
               call check_allocation(iostat, msg, this_module, this_procedure)
               write(msg,'(2(A,E11.4))') 'FGMRES(k)   init step     : |res|= ', &
                        & abs(beta), ', tol= ', tol
               call log_information(msg, this_module, this_procedure)
            end if

            gmres_iter: do k = 1, kdim
                !> Preconditioner.
                call copy(Z(k), V(k)) ; if (ifprecond) call preconditioner%apply(Z(k), k, beta, tol)

                !-----------------------------------------
                !-----     Arnoldi factorization     -----
                !-----------------------------------------
                !> Matrix vector product.
                if (trans) then
                    call A%apply_rmatvec(Z(k), V(k+1))
                else
                    call A%apply_matvec(Z(k), V(k+1))
                endif
                !> Orthogonalization + Hessenberg update.
                call double_gram_schmidt_step(V(k+1), V(:k), info, if_chk_orthonormal=.false., beta=H(:k, k))
                call check_info(info, 'double_gram_schmidt_step', this_module, this_procedure)
                !> Update Hessenberg matrix and normalize residual Krylov vector.
                H(k+1, k) = V(k+1)%norm() 
                if (abs(H(k+1, k)) > tol) call V(k+1)%scal(one_rsp / H(k+1, k))

                !-----------------------------------------
                !-----     Least-Squares Problem     -----
                !-----------------------------------------
                !> Apply Givens rotations to the Hessenberg matrix.
                call apply_givens_rotation(H(:k+1, k), c(:k), s(:k))
                !> Update the right-hand side vector accordingly.
                e(k+1) = -s(k)*e(k) ; e(k) = c(k)*e(k)
                !> Least-squares residual.
                beta = abs(e(k+1))
 
                ! Save metadata.
                fgmres_meta%n_iter  = fgmres_meta%n_iter + 1
                fgmres_meta%n_inner = fgmres_meta%n_inner + 1
                fgmres_meta%res     = [ fgmres_meta%res, abs(beta) ]

                ! Check convergence.
                write(msg,'(A,I3,2(A,E11.4))') 'FGMRES(k)   inner step ', k, ': |res|= ', &
                            & abs(beta), ', tol= ', tol
                call log_information(msg, this_module, this_procedure)

                if (abs(beta) < tol) fgmres_meta%converged = .true.
                if (fgmres_meta%converged) exit gmres_iter
            enddo gmres_iter

            ! Update solution.
            k = min(k, kdim)
            y(1:k, 1:1) => e(:k) ; call trtrs("u", "n", "n", k, 1, H(:k, :k), k, y, k, info)
            call linear_combination(dx, Z(:k), e(:k)) ; call x%add(dx)

            ! Recompute residual for sanity check.
            if (trans) then
                call A%apply_rmatvec(x, v(1))
            else
                call A%apply_matvec(x, v(1))
            endif
            call v(1)%sub(b) ; call v(1)%chsgn()

            ! Initialize new starting Krylov vector if needed.
            beta = v(1)%norm() ; if (abs(beta) > 0.0_sp) call v(1)%scal(one_rsp / beta)

            ! Save metadata.
            fgmres_meta%n_iter  = fgmres_meta%n_iter + 1
            fgmres_meta%n_outer = fgmres_meta%n_outer + 1
            fgmres_meta%res = [ fgmres_meta%res, abs(beta) ]

            write(msg,'(A,I3,2(A,E11.4))') 'FGMRES(k) outer step   ', fgmres_meta%n_outer, ': |res|= ', &
                            & abs(beta), ', tol= ', tol
            call log_information(msg, this_module, this_procedure)

            ! Exit gmres if desired accuracy is reached.
            if (abs(beta) < tol) then
               fgmres_meta%converged = .true.
               exit 
            end if
        enddo
        end associate

        ! Returns the number of iterations if converged
        if (fgmres_meta%converged) then
            info = fgmres_meta%n_iter
        else
            info = -fgmres_meta%n_iter
        end if
        fgmres_meta%info = info

        if (opts%if_print_metadata) call fgmres_meta%print()

        ! Set metadata output
        if (present(meta)) then
            select type(meta)
            type is (fgmres_sp_metadata)
                meta = fgmres_meta
            class default
                call type_error('meta','fgmres_sp_metadata','OUT',this_module,this_procedure)
            end select
        end if

        call A%reset_counter(trans, 'fgmres%post')
        if (time_lightkrylov()) call timer%stop(this_procedure)
    end procedure

    module procedure fgmres_rdp
       ! Options.
        integer :: kdim, maxiter
        real(dp) :: tol, rtol_, atol_
        logical :: trans
        type(fgmres_dp_opts)     :: opts
        type(fgmres_dp_metadata) :: fgmres_meta

        ! Krylov subspace
        class(abstract_vector_rdp), allocatable :: V(:), Z(:)
        ! Hessenberg matrix.
        real(dp), allocatable :: H(:, :)
        ! Least-squares variables.
        real(dp), target, allocatable :: e(:)
        real(dp), pointer :: y(:, :)
        real(dp) :: beta
        ! Givens rotations.
        real(dp), allocatable :: c(:), s(:)

        ! Miscellaneous.
        character(len=*), parameter :: this_procedure = 'fgmres_rdp'
        integer :: k, iostat
        class(abstract_vector_rdp), allocatable :: dx, wrk
        character(len=256) :: msg

        if (time_lightkrylov()) call timer%start(this_procedure)
        ! Deals with the optional args.
        rtol_ = optval(rtol, rtol_dp)
        atol_ = optval(atol, atol_dp)
        if (present(options)) then
            select type (options)
            type is (fgmres_dp_opts)
                opts = options
            class default
                call type_error('options','fgmres_dp_opts','IN',this_module,this_procedure)
            end select
        else
            opts = fgmres_dp_opts()
        endif

        kdim  = opts%kdim ; maxiter = opts%maxiter
        tol   = atol_ + rtol_ * b%norm()
        trans = optval(transpose, .false.)

        ! Initialize working variables.
        allocate(wrk, source=b, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)
        call wrk%zero()
        allocate(V(kdim+1), source=b, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)
        call zero_basis(V)
        allocate(Z(kdim+1), source=b, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)
        call zero_basis(Z)
        allocate(H(kdim+1, kdim), source=zero_rdp, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)
        allocate(e(kdim+1), source=zero_rdp, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)
        allocate(c(kdim), source=zero_rdp, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)
        allocate(s(kdim), source=zero_rdp, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)

        ! Initialize metadata and & reset matvec counter
        fgmres_meta = fgmres_dp_metadata() ; fgmres_meta%converged = .false.
        call A%reset_counter(trans, 'fgmres%init')

        info = 0
        associate(ifprecond => present(preconditioner))
        do while ((.not. fgmres_meta%converged) .and. (fgmres_meta%n_outer <= maxiter))
            !> Initialize data
            H = 0.0_dp ; call zero_basis(V)
            if (x%norm() /= 0.0_dp) then
                if (trans) then
                    call A%apply_rmatvec(x, V(1))
                else
                    call A%apply_matvec(x, V(1))
                endif
            endif
            call V(1)%sub(b) ; call V(1)%chsgn()
            e = 0.0_dp ; beta = V(1)%norm() ; e(1) = beta
            call V(1)%scal(one_rdp/beta)
            c = 0.0_dp ; s = 0.0_dp
            if (fgmres_meta%n_outer == 0) then
               allocate(fgmres_meta%res(1), source=abs(beta), stat=iostat, errmsg=msg)
               call check_allocation(iostat, msg, this_module, this_procedure)
               write(msg,'(2(A,E11.4))') 'FGMRES(k)   init step     : |res|= ', &
                        & abs(beta), ', tol= ', tol
               call log_information(msg, this_module, this_procedure)
            end if

            gmres_iter: do k = 1, kdim
                !> Preconditioner.
                call copy(Z(k), V(k)) ; if (ifprecond) call preconditioner%apply(Z(k), k, beta, tol)

                !-----------------------------------------
                !-----     Arnoldi factorization     -----
                !-----------------------------------------
                !> Matrix vector product.
                if (trans) then
                    call A%apply_rmatvec(Z(k), V(k+1))
                else
                    call A%apply_matvec(Z(k), V(k+1))
                endif
                !> Orthogonalization + Hessenberg update.
                call double_gram_schmidt_step(V(k+1), V(:k), info, if_chk_orthonormal=.false., beta=H(:k, k))
                call check_info(info, 'double_gram_schmidt_step', this_module, this_procedure)
                !> Update Hessenberg matrix and normalize residual Krylov vector.
                H(k+1, k) = V(k+1)%norm() 
                if (abs(H(k+1, k)) > tol) call V(k+1)%scal(one_rdp / H(k+1, k))

                !-----------------------------------------
                !-----     Least-Squares Problem     -----
                !-----------------------------------------
                !> Apply Givens rotations to the Hessenberg matrix.
                call apply_givens_rotation(H(:k+1, k), c(:k), s(:k))
                !> Update the right-hand side vector accordingly.
                e(k+1) = -s(k)*e(k) ; e(k) = c(k)*e(k)
                !> Least-squares residual.
                beta = abs(e(k+1))
 
                ! Save metadata.
                fgmres_meta%n_iter  = fgmres_meta%n_iter + 1
                fgmres_meta%n_inner = fgmres_meta%n_inner + 1
                fgmres_meta%res     = [ fgmres_meta%res, abs(beta) ]

                ! Check convergence.
                write(msg,'(A,I3,2(A,E11.4))') 'FGMRES(k)   inner step ', k, ': |res|= ', &
                            & abs(beta), ', tol= ', tol
                call log_information(msg, this_module, this_procedure)

                if (abs(beta) < tol) fgmres_meta%converged = .true.
                if (fgmres_meta%converged) exit gmres_iter
            enddo gmres_iter

            ! Update solution.
            k = min(k, kdim)
            y(1:k, 1:1) => e(:k) ; call trtrs("u", "n", "n", k, 1, H(:k, :k), k, y, k, info)
            call linear_combination(dx, Z(:k), e(:k)) ; call x%add(dx)

            ! Recompute residual for sanity check.
            if (trans) then
                call A%apply_rmatvec(x, v(1))
            else
                call A%apply_matvec(x, v(1))
            endif
            call v(1)%sub(b) ; call v(1)%chsgn()

            ! Initialize new starting Krylov vector if needed.
            beta = v(1)%norm() ; if (abs(beta) > 0.0_dp) call v(1)%scal(one_rdp / beta)

            ! Save metadata.
            fgmres_meta%n_iter  = fgmres_meta%n_iter + 1
            fgmres_meta%n_outer = fgmres_meta%n_outer + 1
            fgmres_meta%res = [ fgmres_meta%res, abs(beta) ]

            write(msg,'(A,I3,2(A,E11.4))') 'FGMRES(k) outer step   ', fgmres_meta%n_outer, ': |res|= ', &
                            & abs(beta), ', tol= ', tol
            call log_information(msg, this_module, this_procedure)

            ! Exit gmres if desired accuracy is reached.
            if (abs(beta) < tol) then
               fgmres_meta%converged = .true.
               exit 
            end if
        enddo
        end associate

        ! Returns the number of iterations if converged
        if (fgmres_meta%converged) then
            info = fgmres_meta%n_iter
        else
            info = -fgmres_meta%n_iter
        end if
        fgmres_meta%info = info

        if (opts%if_print_metadata) call fgmres_meta%print()

        ! Set metadata output
        if (present(meta)) then
            select type(meta)
            type is (fgmres_dp_metadata)
                meta = fgmres_meta
            class default
                call type_error('meta','fgmres_dp_metadata','OUT',this_module,this_procedure)
            end select
        end if

        call A%reset_counter(trans, 'fgmres%post')
        if (time_lightkrylov()) call timer%stop(this_procedure)
    end procedure

    module procedure fgmres_csp
       ! Options.
        integer :: kdim, maxiter
        real(sp) :: tol, rtol_, atol_
        logical :: trans
        type(fgmres_sp_opts)     :: opts
        type(fgmres_sp_metadata) :: fgmres_meta

        ! Krylov subspace
        class(abstract_vector_csp), allocatable :: V(:), Z(:)
        ! Hessenberg matrix.
        complex(sp), allocatable :: H(:, :)
        ! Least-squares variables.
        complex(sp), target, allocatable :: e(:)
        complex(sp), pointer :: y(:, :)
        real(sp) :: beta
        ! Givens rotations.
        complex(sp), allocatable :: c(:), s(:)

        ! Miscellaneous.
        character(len=*), parameter :: this_procedure = 'fgmres_csp'
        integer :: k, iostat
        class(abstract_vector_csp), allocatable :: dx, wrk
        character(len=256) :: msg

        if (time_lightkrylov()) call timer%start(this_procedure)
        ! Deals with the optional args.
        rtol_ = optval(rtol, rtol_sp)
        atol_ = optval(atol, atol_sp)
        if (present(options)) then
            select type (options)
            type is (fgmres_sp_opts)
                opts = options
            class default
                call type_error('options','fgmres_sp_opts','IN',this_module,this_procedure)
            end select
        else
            opts = fgmres_sp_opts()
        endif

        kdim  = opts%kdim ; maxiter = opts%maxiter
        tol   = atol_ + rtol_ * b%norm()
        trans = optval(transpose, .false.)

        ! Initialize working variables.
        allocate(wrk, source=b, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)
        call wrk%zero()
        allocate(V(kdim+1), source=b, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)
        call zero_basis(V)
        allocate(Z(kdim+1), source=b, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)
        call zero_basis(Z)
        allocate(H(kdim+1, kdim), source=zero_csp, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)
        allocate(e(kdim+1), source=zero_csp, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)
        allocate(c(kdim), source=zero_csp, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)
        allocate(s(kdim), source=zero_csp, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)

        ! Initialize metadata and & reset matvec counter
        fgmres_meta = fgmres_sp_metadata() ; fgmres_meta%converged = .false.
        call A%reset_counter(trans, 'fgmres%init')

        info = 0
        associate(ifprecond => present(preconditioner))
        do while ((.not. fgmres_meta%converged) .and. (fgmres_meta%n_outer <= maxiter))
            !> Initialize data
            H = 0.0_sp ; call zero_basis(V)
            if (x%norm() /= 0.0_sp) then
                if (trans) then
                    call A%apply_rmatvec(x, V(1))
                else
                    call A%apply_matvec(x, V(1))
                endif
            endif
            call V(1)%sub(b) ; call V(1)%chsgn()
            e = 0.0_sp ; beta = V(1)%norm() ; e(1) = beta
            call V(1)%scal(one_csp/beta)
            c = 0.0_sp ; s = 0.0_sp
            if (fgmres_meta%n_outer == 0) then
               allocate(fgmres_meta%res(1), source=abs(beta), stat=iostat, errmsg=msg)
               call check_allocation(iostat, msg, this_module, this_procedure)
               write(msg,'(2(A,E11.4))') 'FGMRES(k)   init step     : |res|= ', &
                        & abs(beta), ', tol= ', tol
               call log_information(msg, this_module, this_procedure)
            end if

            gmres_iter: do k = 1, kdim
                !> Preconditioner.
                call copy(Z(k), V(k)) ; if (ifprecond) call preconditioner%apply(Z(k), k, beta, tol)

                !-----------------------------------------
                !-----     Arnoldi factorization     -----
                !-----------------------------------------
                !> Matrix vector product.
                if (trans) then
                    call A%apply_rmatvec(Z(k), V(k+1))
                else
                    call A%apply_matvec(Z(k), V(k+1))
                endif
                !> Orthogonalization + Hessenberg update.
                call double_gram_schmidt_step(V(k+1), V(:k), info, if_chk_orthonormal=.false., beta=H(:k, k))
                call check_info(info, 'double_gram_schmidt_step', this_module, this_procedure)
                !> Update Hessenberg matrix and normalize residual Krylov vector.
                H(k+1, k) = V(k+1)%norm() 
                if (abs(H(k+1, k)) > tol) call V(k+1)%scal(one_csp / H(k+1, k))

                !-----------------------------------------
                !-----     Least-Squares Problem     -----
                !-----------------------------------------
                !> Apply Givens rotations to the Hessenberg matrix.
                call apply_givens_rotation(H(:k+1, k), c(:k), s(:k))
                !> Update the right-hand side vector accordingly.
                e(k+1) = -s(k)*e(k) ; e(k) = c(k)*e(k)
                !> Least-squares residual.
                beta = abs(e(k+1))
 
                ! Save metadata.
                fgmres_meta%n_iter  = fgmres_meta%n_iter + 1
                fgmres_meta%n_inner = fgmres_meta%n_inner + 1
                fgmres_meta%res     = [ fgmres_meta%res, abs(beta) ]

                ! Check convergence.
                write(msg,'(A,I3,2(A,E11.4))') 'FGMRES(k)   inner step ', k, ': |res|= ', &
                            & abs(beta), ', tol= ', tol
                call log_information(msg, this_module, this_procedure)

                if (abs(beta) < tol) fgmres_meta%converged = .true.
                if (fgmres_meta%converged) exit gmres_iter
            enddo gmres_iter

            ! Update solution.
            k = min(k, kdim)
            y(1:k, 1:1) => e(:k) ; call trtrs("u", "n", "n", k, 1, H(:k, :k), k, y, k, info)
            call linear_combination(dx, Z(:k), e(:k)) ; call x%add(dx)

            ! Recompute residual for sanity check.
            if (trans) then
                call A%apply_rmatvec(x, v(1))
            else
                call A%apply_matvec(x, v(1))
            endif
            call v(1)%sub(b) ; call v(1)%chsgn()

            ! Initialize new starting Krylov vector if needed.
            beta = v(1)%norm() ; if (abs(beta) > 0.0_sp) call v(1)%scal(one_csp / beta)

            ! Save metadata.
            fgmres_meta%n_iter  = fgmres_meta%n_iter + 1
            fgmres_meta%n_outer = fgmres_meta%n_outer + 1
            fgmres_meta%res = [ fgmres_meta%res, abs(beta) ]

            write(msg,'(A,I3,2(A,E11.4))') 'FGMRES(k) outer step   ', fgmres_meta%n_outer, ': |res|= ', &
                            & abs(beta), ', tol= ', tol
            call log_information(msg, this_module, this_procedure)

            ! Exit gmres if desired accuracy is reached.
            if (abs(beta) < tol) then
               fgmres_meta%converged = .true.
               exit 
            end if
        enddo
        end associate

        ! Returns the number of iterations if converged
        if (fgmres_meta%converged) then
            info = fgmres_meta%n_iter
        else
            info = -fgmres_meta%n_iter
        end if
        fgmres_meta%info = info

        if (opts%if_print_metadata) call fgmres_meta%print()

        ! Set metadata output
        if (present(meta)) then
            select type(meta)
            type is (fgmres_sp_metadata)
                meta = fgmres_meta
            class default
                call type_error('meta','fgmres_sp_metadata','OUT',this_module,this_procedure)
            end select
        end if

        call A%reset_counter(trans, 'fgmres%post')
        if (time_lightkrylov()) call timer%stop(this_procedure)
    end procedure

    module procedure fgmres_cdp
       ! Options.
        integer :: kdim, maxiter
        real(dp) :: tol, rtol_, atol_
        logical :: trans
        type(fgmres_dp_opts)     :: opts
        type(fgmres_dp_metadata) :: fgmres_meta

        ! Krylov subspace
        class(abstract_vector_cdp), allocatable :: V(:), Z(:)
        ! Hessenberg matrix.
        complex(dp), allocatable :: H(:, :)
        ! Least-squares variables.
        complex(dp), target, allocatable :: e(:)
        complex(dp), pointer :: y(:, :)
        real(dp) :: beta
        ! Givens rotations.
        complex(dp), allocatable :: c(:), s(:)

        ! Miscellaneous.
        character(len=*), parameter :: this_procedure = 'fgmres_cdp'
        integer :: k, iostat
        class(abstract_vector_cdp), allocatable :: dx, wrk
        character(len=256) :: msg

        if (time_lightkrylov()) call timer%start(this_procedure)
        ! Deals with the optional args.
        rtol_ = optval(rtol, rtol_dp)
        atol_ = optval(atol, atol_dp)
        if (present(options)) then
            select type (options)
            type is (fgmres_dp_opts)
                opts = options
            class default
                call type_error('options','fgmres_dp_opts','IN',this_module,this_procedure)
            end select
        else
            opts = fgmres_dp_opts()
        endif

        kdim  = opts%kdim ; maxiter = opts%maxiter
        tol   = atol_ + rtol_ * b%norm()
        trans = optval(transpose, .false.)

        ! Initialize working variables.
        allocate(wrk, source=b, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)
        call wrk%zero()
        allocate(V(kdim+1), source=b, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)
        call zero_basis(V)
        allocate(Z(kdim+1), source=b, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)
        call zero_basis(Z)
        allocate(H(kdim+1, kdim), source=zero_cdp, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)
        allocate(e(kdim+1), source=zero_cdp, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)
        allocate(c(kdim), source=zero_cdp, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)
        allocate(s(kdim), source=zero_cdp, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)

        ! Initialize metadata and & reset matvec counter
        fgmres_meta = fgmres_dp_metadata() ; fgmres_meta%converged = .false.
        call A%reset_counter(trans, 'fgmres%init')

        info = 0
        associate(ifprecond => present(preconditioner))
        do while ((.not. fgmres_meta%converged) .and. (fgmres_meta%n_outer <= maxiter))
            !> Initialize data
            H = 0.0_dp ; call zero_basis(V)
            if (x%norm() /= 0.0_dp) then
                if (trans) then
                    call A%apply_rmatvec(x, V(1))
                else
                    call A%apply_matvec(x, V(1))
                endif
            endif
            call V(1)%sub(b) ; call V(1)%chsgn()
            e = 0.0_dp ; beta = V(1)%norm() ; e(1) = beta
            call V(1)%scal(one_cdp/beta)
            c = 0.0_dp ; s = 0.0_dp
            if (fgmres_meta%n_outer == 0) then
               allocate(fgmres_meta%res(1), source=abs(beta), stat=iostat, errmsg=msg)
               call check_allocation(iostat, msg, this_module, this_procedure)
               write(msg,'(2(A,E11.4))') 'FGMRES(k)   init step     : |res|= ', &
                        & abs(beta), ', tol= ', tol
               call log_information(msg, this_module, this_procedure)
            end if

            gmres_iter: do k = 1, kdim
                !> Preconditioner.
                call copy(Z(k), V(k)) ; if (ifprecond) call preconditioner%apply(Z(k), k, beta, tol)

                !-----------------------------------------
                !-----     Arnoldi factorization     -----
                !-----------------------------------------
                !> Matrix vector product.
                if (trans) then
                    call A%apply_rmatvec(Z(k), V(k+1))
                else
                    call A%apply_matvec(Z(k), V(k+1))
                endif
                !> Orthogonalization + Hessenberg update.
                call double_gram_schmidt_step(V(k+1), V(:k), info, if_chk_orthonormal=.false., beta=H(:k, k))
                call check_info(info, 'double_gram_schmidt_step', this_module, this_procedure)
                !> Update Hessenberg matrix and normalize residual Krylov vector.
                H(k+1, k) = V(k+1)%norm() 
                if (abs(H(k+1, k)) > tol) call V(k+1)%scal(one_cdp / H(k+1, k))

                !-----------------------------------------
                !-----     Least-Squares Problem     -----
                !-----------------------------------------
                !> Apply Givens rotations to the Hessenberg matrix.
                call apply_givens_rotation(H(:k+1, k), c(:k), s(:k))
                !> Update the right-hand side vector accordingly.
                e(k+1) = -s(k)*e(k) ; e(k) = c(k)*e(k)
                !> Least-squares residual.
                beta = abs(e(k+1))
 
                ! Save metadata.
                fgmres_meta%n_iter  = fgmres_meta%n_iter + 1
                fgmres_meta%n_inner = fgmres_meta%n_inner + 1
                fgmres_meta%res     = [ fgmres_meta%res, abs(beta) ]

                ! Check convergence.
                write(msg,'(A,I3,2(A,E11.4))') 'FGMRES(k)   inner step ', k, ': |res|= ', &
                            & abs(beta), ', tol= ', tol
                call log_information(msg, this_module, this_procedure)

                if (abs(beta) < tol) fgmres_meta%converged = .true.
                if (fgmres_meta%converged) exit gmres_iter
            enddo gmres_iter

            ! Update solution.
            k = min(k, kdim)
            y(1:k, 1:1) => e(:k) ; call trtrs("u", "n", "n", k, 1, H(:k, :k), k, y, k, info)
            call linear_combination(dx, Z(:k), e(:k)) ; call x%add(dx)

            ! Recompute residual for sanity check.
            if (trans) then
                call A%apply_rmatvec(x, v(1))
            else
                call A%apply_matvec(x, v(1))
            endif
            call v(1)%sub(b) ; call v(1)%chsgn()

            ! Initialize new starting Krylov vector if needed.
            beta = v(1)%norm() ; if (abs(beta) > 0.0_dp) call v(1)%scal(one_cdp / beta)

            ! Save metadata.
            fgmres_meta%n_iter  = fgmres_meta%n_iter + 1
            fgmres_meta%n_outer = fgmres_meta%n_outer + 1
            fgmres_meta%res = [ fgmres_meta%res, abs(beta) ]

            write(msg,'(A,I3,2(A,E11.4))') 'FGMRES(k) outer step   ', fgmres_meta%n_outer, ': |res|= ', &
                            & abs(beta), ', tol= ', tol
            call log_information(msg, this_module, this_procedure)

            ! Exit gmres if desired accuracy is reached.
            if (abs(beta) < tol) then
               fgmres_meta%converged = .true.
               exit 
            end if
        enddo
        end associate

        ! Returns the number of iterations if converged
        if (fgmres_meta%converged) then
            info = fgmres_meta%n_iter
        else
            info = -fgmres_meta%n_iter
        end if
        fgmres_meta%info = info

        if (opts%if_print_metadata) call fgmres_meta%print()

        ! Set metadata output
        if (present(meta)) then
            select type(meta)
            type is (fgmres_dp_metadata)
                meta = fgmres_meta
            class default
                call type_error('meta','fgmres_dp_metadata','OUT',this_module,this_procedure)
            end select
        end if

        call A%reset_counter(trans, 'fgmres%post')
        if (time_lightkrylov()) call timer%stop(this_procedure)
    end procedure


    module procedure dense_fgmres_rsp
    type(dense_vector_rsp) :: b_, x_
    type(dense_linop_rsp)  :: A_
    ! Wrap data into convenience types.
    A_ = dense_linop(A)
    b_ = dense_vector(b)
    x_ = dense_vector(x)
    ! Call abstract gmres.
    call fgmres(A_, b_, x_, info, rtol, atol, preconditioner, options, transpose, meta)
    ! Extract solution.
    x = x_%data
    end procedure

    module procedure dense_fgmres_rdp
    type(dense_vector_rdp) :: b_, x_
    type(dense_linop_rdp)  :: A_
    ! Wrap data into convenience types.
    A_ = dense_linop(A)
    b_ = dense_vector(b)
    x_ = dense_vector(x)
    ! Call abstract gmres.
    call fgmres(A_, b_, x_, info, rtol, atol, preconditioner, options, transpose, meta)
    ! Extract solution.
    x = x_%data
    end procedure

    module procedure dense_fgmres_csp
    type(dense_vector_csp) :: b_, x_
    type(dense_linop_csp)  :: A_
    ! Wrap data into convenience types.
    A_ = dense_linop(A)
    b_ = dense_vector(b)
    x_ = dense_vector(x)
    ! Call abstract gmres.
    call fgmres(A_, b_, x_, info, rtol, atol, preconditioner, options, transpose, meta)
    ! Extract solution.
    x = x_%data
    end procedure

    module procedure dense_fgmres_cdp
    type(dense_vector_cdp) :: b_, x_
    type(dense_linop_cdp)  :: A_
    ! Wrap data into convenience types.
    A_ = dense_linop(A)
    b_ = dense_vector(b)
    x_ = dense_vector(x)
    ! Call abstract gmres.
    call fgmres(A_, b_, x_, info, rtol, atol, preconditioner, options, transpose, meta)
    ! Extract solution.
    x = x_%data
    end procedure

end submodule
