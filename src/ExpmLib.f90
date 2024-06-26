module lightkrylov_expmlib

    ! Iso Fortran.
    use iso_fortran_env, only: output_unit

    ! Fortran standard library.
    use stdlib_optval, only: optval
    use stdlib_linalg, only: eye

    ! LightKrylov.
    use Lightkrylov_Constants
    use LightKrylov_Logger
    use LightKrylov_Utils
    use LightKrylov_AbstractVectors
    use LightKrylov_AbstractLinops
    use LightKrylov_BaseKrylov

    implicit none
    
    character(len=128), parameter, private :: this_module = 'LightKrylov_ExpmLib'

    public :: abstract_exptA_rsp
    public :: abstract_exptA_rdp
    public :: abstract_exptA_csp
    public :: abstract_exptA_cdp
    public :: expm
    public :: kexpm
    public :: kexpm_var_dt
    public :: k_exptA

    abstract interface
        subroutine abstract_exptA_rsp(vec_out, A, vec_in, tau, dt, info, trans)
            import sp
            import abstract_vector_rsp
            import abstract_linop_rsp
            !! Abstract interface to define the matrix exponential-vector product.
            class(abstract_vector_rsp), intent(out) :: vec_out
            !! Solution vector.
            class(abstract_linop_rsp), intent(inout) :: A
            !! Linear operator to be exponentiated.
            class(abstract_vector_rsp), intent(in) :: vec_in
            !! Input vector.
            real(sp), intent(in) :: tau
            !! Time horizon for integration.
            real(sp), intent(inout) :: dt
            !! Current time step.
            integer, intent(out) :: info
            !! Information flag.
            logical, optional, intent(in) :: trans
            !! Use transpose ?
        end subroutine abstract_exptA_rsp

        subroutine abstract_exptA_rdp(vec_out, A, vec_in, tau, dt, info, trans)
            import dp
            import abstract_vector_rdp
            import abstract_linop_rdp
            !! Abstract interface to define the matrix exponential-vector product.
            class(abstract_vector_rdp), intent(out) :: vec_out
            !! Solution vector.
            class(abstract_linop_rdp), intent(inout) :: A
            !! Linear operator to be exponentiated.
            class(abstract_vector_rdp), intent(in) :: vec_in
            !! Input vector.
            real(dp), intent(in) :: tau
            !! Time horizon for integration.
            real(dp), intent(inout) :: dt
            !! Current time step.
            integer, intent(out) :: info
            !! Information flag.
            logical, optional, intent(in) :: trans
            !! Use transpose ?
        end subroutine abstract_exptA_rdp

        subroutine abstract_exptA_csp(vec_out, A, vec_in, tau, dt, info, trans)
            import sp
            import abstract_vector_csp
            import abstract_linop_csp
            !! Abstract interface to define the matrix exponential-vector product.
            class(abstract_vector_csp), intent(out) :: vec_out
            !! Solution vector.
            class(abstract_linop_csp), intent(inout) :: A
            !! Linear operator to be exponentiated.
            class(abstract_vector_csp), intent(in) :: vec_in
            !! Input vector.
            real(sp), intent(in) :: tau
            !! Time horizon for integration.
            real(sp), intent(inout) :: dt
            !! Current time step.
            integer, intent(out) :: info
            !! Information flag.
            logical, optional, intent(in) :: trans
            !! Use transpose ?
        end subroutine abstract_exptA_csp

        subroutine abstract_exptA_cdp(vec_out, A, vec_in, tau, dt, info, trans)
            import dp
            import abstract_vector_cdp
            import abstract_linop_cdp
            !! Abstract interface to define the matrix exponential-vector product.
            class(abstract_vector_cdp), intent(out) :: vec_out
            !! Solution vector.
            class(abstract_linop_cdp), intent(inout) :: A
            !! Linear operator to be exponentiated.
            class(abstract_vector_cdp), intent(in) :: vec_in
            !! Input vector.
            real(dp), intent(in) :: tau
            !! Time horizon for integration.
            real(dp), intent(inout) :: dt
            !! Current time step.
            integer, intent(out) :: info
            !! Information flag.
            logical, optional, intent(in) :: trans
            !! Use transpose ?
        end subroutine abstract_exptA_cdp

    end interface

    interface expm
        module procedure expm_rsp
        module procedure expm_rdp
        module procedure expm_csp
        module procedure expm_cdp
    end interface

    interface kexpm
        module procedure kexpm_vec_rsp
        module procedure kexpm_mat_rsp
        module procedure kexpm_vec_rdp
        module procedure kexpm_mat_rdp
        module procedure kexpm_vec_csp
        module procedure kexpm_mat_csp
        module procedure kexpm_vec_cdp
        module procedure kexpm_mat_cdp
    end interface

    interface kexpm_var_dt
        module procedure kexpm_vec_variable_dt_rsp
        module procedure kexpm_mat_variable_dt_rsp
        module procedure kexpm_vec_variable_dt_rdp
        module procedure kexpm_mat_variable_dt_rdp
        module procedure kexpm_vec_variable_dt_csp
        module procedure kexpm_mat_variable_dt_csp
        module procedure kexpm_vec_variable_dt_cdp
        module procedure kexpm_mat_variable_dt_cdp
    end interface

    interface k_exptA
        module procedure k_exptA_rsp
        module procedure k_exptA_rdp
        module procedure k_exptA_csp
        module procedure k_exptA_cdp
    end interface

contains

    !--------------------------------------------
    !-----     DENSE MATRIX EXPONENTIAL     -----
    !--------------------------------------------

    subroutine expm_rsp(E, A, order)
        real(sp), intent(in) :: A(:, :)
        !! Matrix to be exponentiated.
        real(sp), intent(out) :: E(:, :)
        !! Output matrix E = exp(tA).
        integer, optional :: order
        !! Order of the Pade approximation.

        !----- Internal variables -----
        real(sp), allocatable :: A2(:, :), Q(:, :), X(:, :), invQ(:, :), wrk(:)
        real(sp) :: a_norm, c
        integer :: n, ee, k, s
        logical :: p
        integer :: p_order

        ! Deal with optional args.
        p_order = optval(order, 10)

        n = size(A, 1)

        ! Allocate arrays.
        allocate(A2(n, n)) ; allocate(X(n, n))
        allocate(Q(n, n)) ; allocate(invQ(n, n))
        allocate(wrk(n))

        ! Compute the L-infinity norm.
        a_norm = norml_rsp(A)

        ! Determine scaling factor for the matrix.
        ee = int(log2_rsp(a_norm)) + 1
        s = max(0, ee+1)

        ! Scale the input matrix & initialize polynomial.
        A2 = A / 2.0_sp**s
        X = A2

        ! Initialize P & Q and add first step.
        c = 0.5_sp
        E = eye(n) ; E = E + c*A2

        Q = eye(n) ; Q = Q - c*A2

        ! Iteratively compute the Pade approximation.
        p = .true.
        do k = 2, p_order
            c = c*(p_order - k + 1) / (k * (2*p_order - k + 1))
            X = matmul(A2, X)
            E = E + c*X
            if (p) then
                Q = Q + c*X
            else
                Q = Q - c*X
            endif
            p = .not. p
        enddo

        invQ = Q ; call inv(invQ)
        E = matmul(invQ, E)

        do k = 1, s
            E = matmul(E, E)
        enddo

        return
    end subroutine expm_rsp

    subroutine kexpm_vec_rsp(c, A, b, tau, tol, info, trans, verbosity, kdim)
        class(abstract_vector_rsp), intent(out) :: c
        !! Best approximation of \( \exp(\tau \mathbf{A}) \mathbf{b} \) in the computed Krylov subspace
        class(abstract_linop_rsp), intent(in) :: A
        !! Linear operator to be exponentiated.
        class(abstract_vector_rsp), intent(in) :: b
        !! Input vector on which to apply \( \exp(\tau \mathbf{A}) \).
        real(sp), intent(in) :: tau
        !! Time horizon for the exponentiation.
        real(sp), intent(in) :: tol
        !! Solution tolerance based on error estimates.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: trans
        !! Use transpose?
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        integer, optional, intent(in) :: kdim
        !! Maximum size of the Krylov subspace.

        ! ----- Internal variables -----
        integer, parameter :: kmax = 100
        integer :: k, km, kp, nk
        ! Arnoldi factorization.
        class(abstract_vector_rsp), allocatable :: X(:)
        real(sp), allocatable :: H(:, :)
        ! Normaliztion & temporary arrays.
        real(sp), allocatable :: E(:, :)
        class(abstract_vector_rsp), allocatable :: Xwrk
        real(sp) :: err_est, beta
        character(len=512) :: msg

        ! Optional arguments.
        logical :: transpose
        logical :: verbose
        integer :: nsteps
    
        ! Deals with optional args.
        transpose = optval(trans, .false.)
        verbose   = optval(verbosity, .false.)
        nsteps    = optval(kdim, kmax)
        nk        = nsteps

        info = 0

        ! Allocate arrays.
        allocate(X(nk+1), source=b) ; allocate(H(nk+1, nk+1))
        allocate(E(nk+1, nk+1)) ; allocate(Xwrk, source=b)

        ! Normalize input vector and initialize Krylov subspace.
        beta = b%norm()
        if (beta == 0.0_sp) then
            ! Input is zero => Output is zero.
            call c%zero()
            err_est = 0.0_sp
            kp = 1
        else
            call initialize_krylov_subspace(X)
            call X(1)%add(b) ; call X(1)%scal(one_rsp / beta)
            H = 0.0_sp

            expm_arnoldi: do k = 1, nk
                km = k-1 ; kp = k+1

                ! Reset work arrays.
                E = 0.0_sp

                ! Compute k-th step Arnoldi factorization.
                call arnoldi(A, X, H, info, kstart=k, kend=k, transpose=transpose)
                call check_info(info, 'arnoldi', module=this_module, procedure='kexpm_vec_rsp')

                ! Compute approximation.
                if (info == k) then
                    ! Arnoldi breakdown, do not consider extended matrix.
                    kp = k
                    info = -2
                endif

                ! Compute the (dense) matrix exponential of the extended Hessenberg matrix.
                call expm(E(:kp, :kp), tau*H(:kp, :kp))

                ! Project back into original space.
                call linear_combination(Xwrk, X(:kp), E(:kp, 1))
                call c%axpby(zero_rsp, Xwrk, one_rsp*beta)

                ! Cheap error esimate (this actually is the magnitude of the included correction
                ! and thus is very conservative).
                if (info == k) then
                    ! Approximation is exact.
                    err_est = 0.0_sp
                else
                    err_est = abs(E(kp, 1) * beta)
                endif

                ! Check convergence.
                if (err_est <= tol) exit expm_arnoldi

            enddo expm_arnoldi
        endif

        if (err_est <= tol) then
            info = kp
            write(msg, *) new_line('A'),' Arnoldi approxmation of the action of the exp. propagator converged!'
            write(msg, '(A,4X,A,I3)') trim(msg)//new_line('A'), 'n° of vectors:', kp
            write(msg, '(A,4X,A,E10.4)') trim(msg)//new_line('A'), 'estimated error:   ||err_est||_2 = ', err_est
            write(msg, '(A,4X,A,E10.4)') trim(msg)//new_line('A'), 'desired tolerance:           tol = ', tol
        else
            info = -1
            write(msg, *) new_line('A'),' Arnoldi approxmation of the action of the exp. propagator did not converge!'
            write(msg, '(A,4X,A,I3)') trim(msg)//new_line('A'), 'Maximum n° of vectors reached: ', nk + 1
            write(msg, '(A,4X,A,E10.4)') trim(msg)//new_line('A'), 'Estimated error              :   ||err_est||_2 = ', err_est
            write(msg, '(A,4X,A,E10.4)') trim(msg)//new_line('A'), 'Desired tolerance            :           tol = ', tol
        endif
        if (verbose) call logger%log_message(trim(msg), module=this_module, procedure='kexpm_vec_rsp')

        return
    end subroutine kexpm_vec_rsp

    subroutine kexpm_mat_rsp(C, A, B, tau, tol, info, trans, verbosity, kdim)
        class(abstract_vector_rsp), intent(out) :: C(:)
        !! Best Krylov approximation of \( \mathbf{C} = \exp(\tau \mathbf{A}) \mathbf{B} \).
        class(abstract_linop_rsp), intent(in) :: A
        !! Linear operator to be exponentiated.
        class(abstract_vector_rsp), intent(in) :: B(:)
        !! Input matrix on which to apply \( \exp(\tau \mathbf{A}) \).
        real(sp), intent(in) :: tau
        !! Time horizon for the exponentiation.
        real(sp), intent(in) :: tol
        !! Solution toleance based on error estimates.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: trans
        !! Use transpose ?
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        integer, optional, intent(in) :: kdim
        !! Maximum size of the Krylov subspace.

        ! ----- Internal variables -----
        integer, parameter :: kmax = 100
        integer :: i, j, k, p, kpm, kp, kpp, nk
        ! Block-Arnoldi factorization.
        class(abstract_vector_rsp), allocatable :: X(:)
        real(sp), allocatable :: H(:, :)
        ! Normalization & temporary arrays.
        real(sp), allocatable :: R(:, :), E(:, :), em(:, :)
        integer, allocatable :: perm(:), ptrans(:)
        class(abstract_vector_rsp), allocatable :: Xwrk(:), Cwrk(:)
        real(sp) :: err_est
        character(len=512) :: msg

        ! Optional arguments.
        logical :: transpose
        logical :: verbose
        integer :: nsteps

        ! Determine block size.
        p = size(B)

        ! Deals with the optional args.
        transpose = optval(trans, .false.)
        verbose   = optval(verbosity, .false.)
        nsteps    = optval(kdim, kmax)
        nk        = nsteps*p
        
        info = 0

        ! Allocate arrays.
        allocate(R(p, p)) ; allocate(perm(p)) ; allocate(ptrans(p))
        allocate(X(p*(nk+1)), source=B(1)) ; allocate(H(p*(nk+1), p*(nk+1)))
        allocate(E(p*(nk+1), p*(nk+1))) ; allocate(em(p, p))

        ! Scratch arrays.
        allocate(Xwrk(p), source=B) ; allocate(Cwrk(p), source=B(1))

        ! Normalize input matrix and initialize Krylov subspace.
        R = 0.0_sp
        call qr(Xwrk, R, perm, info) ; call apply_inverse_permutation_matrix(R, perm)

        if (norm2(abs(R)) == 0.0_sp) then
            ! Input matrix is zero.
            do i = 1, size(C)
                call C(i)%zero()
            enddo
            err_est = 0.0_sp ; k = 0 ; kpp = p
        else
            call initialize_krylov_subspace(X, Xwrk) ; H = 0.0_sp

            expm_arnoldi: do k = 1, nk
                ! Set counters.
                kpm = (k-1)*p ; kp = kpm + p ; kpp = kp + p

                ! Reset working arrays.
                E = 0.0_sp
                do i = 1, size(Xwrk)
                    call Xwrk(i)%zero()
                enddo

                ! Compute the k-th step of the Arnoldi factorization.
                call arnoldi(A, X, H, info, kstart=k, kend=k, transpose=transpose, blksize=p)
                call check_info(info, 'arnoldi', module=this_module, procedure='kexpm_mat_rsp')


                if (info == kp) then
                    ! Arnoldi breakdown. Do not consider extended matrix.
                    kpp = kp
                    info = -2
                endif

                ! Compute the (dense) matrix exponential of the extended Hessenberg matrix.
                call expm(E(:kpp, :kpp), tau*H(:kpp, :kpp))

                ! Project back to original space.
                do i = 1, size(Xwrk)
                    call Xwrk(i)%zero()
                    do j = 1, kpp
                        call Xwrk(i)%axpby(one_rsp, X(j), E(j, i))
                    enddo
                enddo

                do i = 1, p
                    call C(i)%zero()
                    do j = 1, p
                        call C(i)%axpby(one_rsp, Xwrk(j), R(j, i))
                    enddo
                enddo

                ! Cheap error estimate.
                if (info == kp) then
                    ! Approximation is exact.
                    err_est = 0.0_sp
                else
                    em = matmul(E(kp+1:kpp, :p), R(:p, :p))
                    err_est = norm2(abs(em))
                endif

                if (err_est <= tol) exit expm_arnoldi

            enddo expm_arnoldi
        endif

        if (err_est .le. tol) then
            info = kpp
            if (p.eq.1) then
               write(msg, *) new_line('A'),' Arnoldi approxmation of the action of the exp. propagator converged!'
            else
               write(msg, *) new_line('A'),' Block Arnoldi approxmation of the action of the exp. propagator converged!'
            endif 
            write(msg, '(A,2(A,I3))') trim(msg)//new_line('A'), 'n° of vectors:', k+1, ' per input vector, total:', kpp
            write(msg, '(A,4X,A,E10.4)') trim(msg)//new_line('A'), 'estimated error:   ||err_est||_2 = ', err_est
            write(msg, '(A,4X,A,E10.4)') trim(msg)//new_line('A'), 'desired tolerance:           tol = ', tol
        else
            info = -1
            if (p.eq.1) then
                write(msg, *) 'Arnoldi approxmation of the action of the exp. propagator did not converge!'
            else
                write(msg, *) 'Block Arnoldi approxmation of the action of the exp. propagator did not converge!'
            endif
            write(msg, '(A,4X,2(A,I3))') trim(msg)//new_line('A'), 'maximum n° of vectors reached: ', nsteps+1, &
                                   & ' per input vector, total:', kpp
            write(msg, '(A,4X,A,E10.4)') trim(msg)//new_line('A'), 'estimated error:   ||err_est||_2 = ', err_est
            write(msg, '(A,4X,A,E10.4)') trim(msg)//new_line('A'), 'desired tolerance:           tol = ', tol
        endif

        if (verbose) call logger%log_message(trim(msg), module=this_module, procedure='kexpm_mat_rsp')

        return
    end subroutine kexpm_mat_rsp

    subroutine kexpm_vec_variable_dt_rsp(c, A, b, tau, dt, tol, info, trans, verbosity, kdim)
        class(abstract_vector_rsp), intent(out) :: c
        !! Best approximation of \( \exp(\tau \mathbf{A}) \mathbf{b} \) in the computed Krylov subspace
        class(abstract_linop_rsp), intent(in) :: A
        !! Linear operator to be exponentiated.
        class(abstract_vector_rsp), intent(in) :: b
        !! Input vector on which to apply \( \exp(\tau \mathbf{A}) \).
        real(sp), intent(in) :: tau
        !! Time horizon for the exponentiation.
        real(sp), intent(inout) :: dt
        !! Current time step.
        real(sp), intent(in) :: tol
        !! Solution tolerance based on error estimates.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: trans
        !! Use transpose?
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        integer, optional, intent(in) :: kdim
        !! Maximum size of the Krylov subspace.

        ! ----- Internal variables -----
        integer, parameter :: kmax = 100
        integer :: k, km, kp, nk, tstep, rstep
        ! Arnoldi factorization.
        class(abstract_vector_rsp), allocatable :: X(:)
        real(sp), allocatable :: H(:, :)
        ! Normaliztion & temporary arrays.
        real(sp), allocatable :: E(:, :)
        class(abstract_vector_rsp), allocatable :: b_
        real(sp) :: err_est, beta, Tend, dt_, tol_, scale, sfact
        character(len=128) :: msg

        ! Optional arguments.
        logical :: transpose
        logical :: verbose
        integer :: nsteps
    
        ! Deals with optional args.
        transpose = optval(trans, .false.)
        verbose   = optval(verbosity, .false.)
        nsteps    = optval(kdim, kmax)
        nk        = nsteps

        info = 0
        scale = 1.5_sp
        sfact = 0.5_sp

        ! Allocate arrays.
        allocate(X(nk+1), b_, source=b)
        allocate(H(nk+1, nk+1), E(nk+1, nk+1))

        ! Normalize input vector and initialize Krylov subspace.
        beta = b%norm()
        if (beta == 0.0_sp) then
            ! Input is zero => Output is zero.
            call c%zero()
        else
            tstep = 1
            Tend = 0.0_sp
            if (verbose) then
               write(msg, '(3(A,E12.6))') 'tau = ', tau, ', dt_0 = ', dt, ', tol = ', tol
               call logger%log_message(trim(msg), module=this_module, procedure='kexpm_vec_variable_dt_rsp')
            end if
            do while ( Tend < tau )

               ! (Re-)initialization of X
               call initialize_krylov_subspace(X) ; call X(1)%add(b_) ;
               beta = b_%norm() ; call X(1)%scal(one_rsp / beta)
               H = 0.0_sp

               ! choose time step
               dt  = min(tau, scale*dt)    ! Reuse step size of last call (and try to increase)
               if (dt > tau-Tend) then
                  dt_ = tau-Tend       ! If a small final step is needed to reach Tend, ignore it on restart
                  if (verbose) then
                     write(msg, '(A,E12.6)') 'Filler step, set dt_ = ', dt_
                     call logger%log_message(trim(msg), module=this_module, procedure='kexpm_vec_variable_dt_rsp')
                  end if
               else
                  dt_ = dt
               end if
               tol_ = dt_/tau*tol*sfact
               if (verbose) then
                  write(msg, '(A,I3,4(A,E8.2))') 'step ', tstep, ': T= ', Tend, ', dt= ', dt, ', dt_= ', dt_, ', tol_= ', tol_
                  call logger%log_message(trim(msg), module=this_module, procedure='kexpm_vec_variable_dt_rsp')
               end if

               expm_arnoldi: do k = 1, nk
                   km = k-1 ; kp = k+1

                   ! Reset work arrays.
                   E = 0.0_sp

                   ! Compute k-th step Arnoldi factorization.
                   call arnoldi(A, X, H, info, kstart=k, kend=k, transpose=transpose)
                   call check_info(info, 'arnoldi', module=this_module, procedure='kexpm_vec_variable_dt_rsp')

                   ! Compute approximation.
                   if (info == k) then
                       ! Arnoldi breakdown, do not consider extended matrix.
                       kp = k
                       info = -2
                   endif

                   ! Compute the (dense) matrix exponential of the extended Hessenberg matrix.
                   call expm(E(:kp, :kp), dt_*H(:kp, :kp))

                   ! Cheap error esimate (this actually is the magnitude of the included correction
                   ! and thus is very conservative).
                   if (info == k) then
                       ! Approximation is exact.
                       err_est = 0.0_sp
                   else
                       err_est = abs(E(kp, 1) * beta)
                   endif

                   ! Check convergence.
                   if (err_est <= tol_) exit expm_arnoldi

               enddo expm_arnoldi

               if (verbose) then
                  write(msg, '(34X,A,E12.6)') 'err_est = ', err_est
                  call logger%log_message(trim(msg), module=this_module, procedure='kexpm_vec_variable_dt_rsp')
               end if

               if (err_est > tol_) then ! decrease dt_ until error estimate is below tolerance
                                       ! NOTE: We do not need to recompute X,H!
                  rstep = 1
                  do while (err_est > tol_)

                     !hfact = (tol/err_est)**(one_rsp/(kdim+1))
                     !dt_ = 0.5_sp*sfact_sp*hfact*dt_ 
                     dt_ = dt_/scale ! this is a quick hack, should choose step better

                     tol_ = dt_/tau*tol*sfact ! update accepted tol based on new time step

                     if (dt_ < atol_sp) then
                        write(msg, *) 'dt < atol_sp. Cannot compute action of the matrix exponential!'
                        call stop_error(trim(msg), module=this_module, procedure='kexpm_vec_variable_dt_rsp')
                     end if

                     ! Compute the (dense) matrix exponential of the extended Hessenberg matrix.
                     call expm(E(:kp, :kp), dt_*H(:kp, :kp))

                     ! Compute new error estimate
                     err_est = abs(E(kp, 1) * beta)

                     if (verbose) then
                        write(msg, '(6X,3(A,E10.4))') 'v dt_ = ', dt_, ', err_est = ', err_est,',  tol_ = ', tol_
                        call logger%log_message(trim(msg), module=this_module, procedure='kexpm_vec_variable_dt_rsp')
                     end if

                     rstep = rstep + 1

                  end do

                  ! save dt for reuse
                  dt = dt_

               end if

               ! Project back into original space.
               block
                  class(abstract_vector_rsp), allocatable :: xwrk
                  call linear_combination(xwrk, X(:kp), E(:kp, 1))
                  call b_%axpby(zero_rsp, xwrk, one_rsp*beta)
               end block

               ! move forward in time
               Tend = Tend + dt_

               tstep = tstep + 1

               if (verbose) then
                  write(msg, '(A,E10.4,2(A,I3),A,E10.4)') '    post: T = ', Tend, ', k = ', k, '/', kdim+1, ', dt => ', dt
                  call logger%log_message(trim(msg), module=this_module, procedure='kexpm_vec_variable_dt_rsp')
               end if
           
            end do ! integration

            ! copy result into output vector
            call c%axpby(zero_rsp, b_, one_rsp)

        endif      

        return
    end subroutine kexpm_vec_variable_dt_rsp

    subroutine kexpm_mat_variable_dt_rsp(C, A, B, tau, dt, tol, info, trans, verbosity, kdim)
        class(abstract_vector_rsp), intent(out) :: C(:)
        !! Best Krylov approximation of \( \mathbf{C} = \exp(\tau \mathbf{A}) \mathbf{B} \).
        class(abstract_linop_rsp), intent(in) :: A
        !! Linear operator to be exponentiated.
        class(abstract_vector_rsp), intent(in) :: B(:)
        !! Input matrix on which to apply \( \exp(\tau \mathbf{A}) \).
        real(sp), intent(in) :: tau
        !! Time horizon for the exponentiation.
        real(sp), intent(inout) :: dt
        !! Current time step.
        real(sp), intent(in) :: tol
        !! Solution toleance based on error estimates.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: trans
        !! Use transpose ?
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        integer, optional, intent(in) :: kdim
        !! Maximum size of the Krylov subspace.

        ! ----- Internal variables -----
        integer, parameter :: kmax = 100
        integer :: i, j, k, p, kpm, kp, kpp, nk, tstep
        ! Block-Arnoldi factorization.
        class(abstract_vector_rsp), allocatable :: X(:)
        real(sp), allocatable :: H(:, :)
        ! Normalization & temporary arrays.
        real(sp), allocatable :: R(:, :), E(:, :), em(:, :)
        integer, allocatable :: perm(:), ptrans(:)
        class(abstract_vector_rsp), allocatable :: Xwrk(:), Cwrk(:), B_(:)
        real(sp) :: err_est, Tend, dt_, tol_, scale
        character(len=128) :: msg

        ! Optional arguments.
        logical :: transpose
        logical :: verbose
        integer :: nsteps

        ! Determine block size.
        p = size(B)

        ! Deals with the optional args.
        transpose = optval(trans, .false.)
        verbose   = optval(verbosity, .false.)
        nsteps    = optval(kdim, kmax)
        nk        = nsteps*p
        
        info = 0
        scale = 1.5

        ! Allocate arrays.
        allocate(perm(p), ptrans(p))
        allocate(H(p*(nk+1), p*(nk+1)), E(p*(nk+1), p*(nk+1)), R(p, p), em(p, p))
        allocate(X(p*(nk+1)), source=B(1))

        ! Scratch arrays.
        allocate(B_(p), source=B) ; allocate(Cwrk(p), source=B(1))

        if (norm2(abs(R)) == 0.0_sp) then
            ! Input matrix is zero.
            call zero_basis(C)
        else
            tstep = 1
            Tend = 0.0_sp
            do while ( Tend < tau )

               ! Normalize B_ and (re-)initialize Krylov subspace.
               R = 0.0_sp
               call qr(B_, R, perm, info) ; call apply_inverse_permutation_matrix(R, perm)
               call initialize_krylov_subspace(X, B_)
               H = 0.0_sp

               ! choose time step
               dt  = min(tau, scale*dt)    ! Reuse step size of last call (and try to increase)
               if (dt > tau-Tend) then
                  dt_ = tau-Tend       ! If a small final step is needed to reach Tend, ignore it on restart
                  if (verbose) then
                     write(msg, '(A,E12.6)') 'Filler step, set dt_ = ', dt_
                     call logger%log_message(trim(msg), module=this_module, procedure='kexpm_mat_variable_dt_rsp')
                  end if
               else
                  dt_ = dt
               end if
               tol_ = dt_/tau*tol
               if (verbose) then
                  write(msg, '(A,I3,4(A,E8.2))') 'step ', tstep, ': T= ', Tend, ', dt= ', dt, ', dt_= ', dt_, ', tol_= ', tol_
                  call logger%log_message(trim(msg), module=this_module, procedure='kexpm_mat_variable_dt_rsp')
               end if

               expm_arnoldi: do k = 1, nk
                   ! Set counters.
                   kpm = (k-1)*p ; kp = kpm + p ; kpp = kp + p

                   ! Reset working arrays.
                   E = 0.0_sp

                   ! Compute the k-th step of the Arnoldi factorization.
                   call arnoldi(A, X, H, info, kstart=k, kend=k, transpose=transpose, blksize=p)
                   call check_info(info, 'arnoldi', module=this_module, procedure='kexpm_mat_variable_dt_rsp')

                   if (info == kp) then
                       ! Arnoldi breakdown. Do not consider extended matrix.
                       kpp = kp
                       info = -2
                   endif

                   ! Compute the (dense) matrix exponential of the extended Hessenberg matrix.
                   call expm(E(:kpp, :kpp), dt_*H(:kpp, :kpp))

                   ! Cheap error estimate.
                   if (info == kp) then
                       ! Approximation is exact.
                       err_est = 0.0_sp
                   else
                       em = matmul(E(kp+1:kpp, :p), R(:p, :p))
                       err_est = norm2(abs(em))
                   endif

                   if (err_est <= tol_) exit expm_arnoldi

               enddo expm_arnoldi

               if (verbose) then
                  write(msg, '(34X,A,E10.4)') 'err_est = ', err_est
                  call logger%log_message(trim(msg), module=this_module, procedure='kexpm_mat_variable_dt_rsp')
               end if

               if (err_est > tol_) then ! decrease dt_ until error estimate is below tolerance
                                       ! NOTE: We do not need to recompute X,H!
                  do while (err_est > tol_)

                     !hfact = (tol/err_est)**(one_rsp/kdim)
                     !dt_ = sfact_sp*hfact*dt_
                     dt_ = dt_/scale ! this is a quick hack, should choose step better

                     tol_ = dt_/tau*tol ! update accepted tol based on new time step

                     if (dt_ < atol_sp) then
                        write(msg, *) 'dt < atol_sp. Cannot compute action of the matrix exponential!'
                        call stop_error(trim(msg), module=this_module, procedure='kexpm_mat_variable_dt_rsp')
                     end if

                     ! Compute the (dense) matrix exponential of the extended Hessenberg matrix.
                     call expm(E(:kpp, :kpp), dt_*H(:kpp, :kpp))

                     ! Compute new error estimate
                     em = matmul(E(kp+1:kpp, :p), R(:p, :p))
                     err_est = norm2(abs(em))

                     if (verbose) then
                        write(msg, '(6X,3(A,E12.6))') 'v dt_ = ', dt_, ', err_est = ', err_est, ', tol_ = ', tol_
                        call logger%log_message(trim(msg), module=this_module, procedure='kexpm_mat_variable_dt_rsp')
                     end if

                  end do

                  ! save dt for reuse
                  dt = dt_
               
               end if

               ! Project back into original space.
               block
                  class(abstract_vector_rsp), allocatable :: Xwrk(:)
                  call linear_combination(Xwrk, X(:kp), matmul(E(:kp, :p), R))
               end block

               ! move forward in time
               Tend = Tend + dt_

               tstep = tstep + 1

               if (verbose) then
                  write(msg, '(A,E10.4,2(A,I3),A,E10.4)') '    post: T = ', Tend, ', k = ', k, '/', kdim+1, ', dt => ', dt
                  call logger%log_message(trim(msg), module=this_module, procedure='kexpm_mat_variable_dt_rsp')
               end if
           
            end do ! integration

            ! copy result into output vector
            call copy_basis(C, B_)
        endif

        return
    end subroutine kexpm_mat_variable_dt_rsp

    subroutine k_exptA_rsp(vec_out, A, vec_in, tau, info, trans)
        class(abstract_vector_rsp), intent(out) :: vec_out
        !! Solution vector.
        class(abstract_linop_rsp), intent(inout) :: A
        !! Linear operator to be exponentiated.
        class(abstract_vector_rsp), intent(in) :: vec_in
        !! Input vector to be multiplied by \( \exp(\tau \mathbf{A}) \).
        real(sp), intent(in) :: tau
        !! Time horizon for the exponentiation.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: trans
        !! Use adjoint ?

        ! ----- Internal variables -----
        real(sp) :: tol
        integer :: kdim
        logical :: verbose

        tol = rtol_sp
        kdim = 30
        verbose = .false.

        call kexpm(vec_out, A, vec_in, tau, tol, info, trans=trans, verbosity=verbose, kdim=kdim)
        call check_info(info, 'kexpm', module=this_module, procedure='k_exptA_rsp')

        return
    end subroutine k_exptA_rsp

    subroutine expm_rdp(E, A, order)
        real(dp), intent(in) :: A(:, :)
        !! Matrix to be exponentiated.
        real(dp), intent(out) :: E(:, :)
        !! Output matrix E = exp(tA).
        integer, optional :: order
        !! Order of the Pade approximation.

        !----- Internal variables -----
        real(dp), allocatable :: A2(:, :), Q(:, :), X(:, :), invQ(:, :), wrk(:)
        real(dp) :: a_norm, c
        integer :: n, ee, k, s
        logical :: p
        integer :: p_order

        ! Deal with optional args.
        p_order = optval(order, 10)

        n = size(A, 1)

        ! Allocate arrays.
        allocate(A2(n, n)) ; allocate(X(n, n))
        allocate(Q(n, n)) ; allocate(invQ(n, n))
        allocate(wrk(n))

        ! Compute the L-infinity norm.
        a_norm = norml_rdp(A)

        ! Determine scaling factor for the matrix.
        ee = int(log2_rdp(a_norm)) + 1
        s = max(0, ee+1)

        ! Scale the input matrix & initialize polynomial.
        A2 = A / 2.0_dp**s
        X = A2

        ! Initialize P & Q and add first step.
        c = 0.5_dp
        E = eye(n) ; E = E + c*A2

        Q = eye(n) ; Q = Q - c*A2

        ! Iteratively compute the Pade approximation.
        p = .true.
        do k = 2, p_order
            c = c*(p_order - k + 1) / (k * (2*p_order - k + 1))
            X = matmul(A2, X)
            E = E + c*X
            if (p) then
                Q = Q + c*X
            else
                Q = Q - c*X
            endif
            p = .not. p
        enddo

        invQ = Q ; call inv(invQ)
        E = matmul(invQ, E)

        do k = 1, s
            E = matmul(E, E)
        enddo

        return
    end subroutine expm_rdp

    subroutine kexpm_vec_rdp(c, A, b, tau, tol, info, trans, verbosity, kdim)
        class(abstract_vector_rdp), intent(out) :: c
        !! Best approximation of \( \exp(\tau \mathbf{A}) \mathbf{b} \) in the computed Krylov subspace
        class(abstract_linop_rdp), intent(in) :: A
        !! Linear operator to be exponentiated.
        class(abstract_vector_rdp), intent(in) :: b
        !! Input vector on which to apply \( \exp(\tau \mathbf{A}) \).
        real(dp), intent(in) :: tau
        !! Time horizon for the exponentiation.
        real(dp), intent(in) :: tol
        !! Solution tolerance based on error estimates.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: trans
        !! Use transpose?
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        integer, optional, intent(in) :: kdim
        !! Maximum size of the Krylov subspace.

        ! ----- Internal variables -----
        integer, parameter :: kmax = 100
        integer :: k, km, kp, nk
        ! Arnoldi factorization.
        class(abstract_vector_rdp), allocatable :: X(:)
        real(dp), allocatable :: H(:, :)
        ! Normaliztion & temporary arrays.
        real(dp), allocatable :: E(:, :)
        class(abstract_vector_rdp), allocatable :: Xwrk
        real(dp) :: err_est, beta
        character(len=512) :: msg

        ! Optional arguments.
        logical :: transpose
        logical :: verbose
        integer :: nsteps
    
        ! Deals with optional args.
        transpose = optval(trans, .false.)
        verbose   = optval(verbosity, .false.)
        nsteps    = optval(kdim, kmax)
        nk        = nsteps

        info = 0

        ! Allocate arrays.
        allocate(X(nk+1), source=b) ; allocate(H(nk+1, nk+1))
        allocate(E(nk+1, nk+1)) ; allocate(Xwrk, source=b)

        ! Normalize input vector and initialize Krylov subspace.
        beta = b%norm()
        if (beta == 0.0_dp) then
            ! Input is zero => Output is zero.
            call c%zero()
            err_est = 0.0_dp
            kp = 1
        else
            call initialize_krylov_subspace(X)
            call X(1)%add(b) ; call X(1)%scal(one_rdp / beta)
            H = 0.0_dp

            expm_arnoldi: do k = 1, nk
                km = k-1 ; kp = k+1

                ! Reset work arrays.
                E = 0.0_dp

                ! Compute k-th step Arnoldi factorization.
                call arnoldi(A, X, H, info, kstart=k, kend=k, transpose=transpose)
                call check_info(info, 'arnoldi', module=this_module, procedure='kexpm_vec_rdp')

                ! Compute approximation.
                if (info == k) then
                    ! Arnoldi breakdown, do not consider extended matrix.
                    kp = k
                    info = -2
                endif

                ! Compute the (dense) matrix exponential of the extended Hessenberg matrix.
                call expm(E(:kp, :kp), tau*H(:kp, :kp))

                ! Project back into original space.
                call linear_combination(Xwrk, X(:kp), E(:kp, 1))
                call c%axpby(zero_rdp, Xwrk, one_rdp*beta)

                ! Cheap error esimate (this actually is the magnitude of the included correction
                ! and thus is very conservative).
                if (info == k) then
                    ! Approximation is exact.
                    err_est = 0.0_dp
                else
                    err_est = abs(E(kp, 1) * beta)
                endif

                ! Check convergence.
                if (err_est <= tol) exit expm_arnoldi

            enddo expm_arnoldi
        endif

        if (err_est <= tol) then
            info = kp
            write(msg, *) new_line('A'),' Arnoldi approxmation of the action of the exp. propagator converged!'
            write(msg, '(A,4X,A,I3)') trim(msg)//new_line('A'), 'n° of vectors:', kp
            write(msg, '(A,4X,A,E10.4)') trim(msg)//new_line('A'), 'estimated error:   ||err_est||_2 = ', err_est
            write(msg, '(A,4X,A,E10.4)') trim(msg)//new_line('A'), 'desired tolerance:           tol = ', tol
        else
            info = -1
            write(msg, *) new_line('A'),' Arnoldi approxmation of the action of the exp. propagator did not converge!'
            write(msg, '(A,4X,A,I3)') trim(msg)//new_line('A'), 'Maximum n° of vectors reached: ', nk + 1
            write(msg, '(A,4X,A,E10.4)') trim(msg)//new_line('A'), 'Estimated error              :   ||err_est||_2 = ', err_est
            write(msg, '(A,4X,A,E10.4)') trim(msg)//new_line('A'), 'Desired tolerance            :           tol = ', tol
        endif
        if (verbose) call logger%log_message(trim(msg), module=this_module, procedure='kexpm_vec_rdp')

        return
    end subroutine kexpm_vec_rdp

    subroutine kexpm_mat_rdp(C, A, B, tau, tol, info, trans, verbosity, kdim)
        class(abstract_vector_rdp), intent(out) :: C(:)
        !! Best Krylov approximation of \( \mathbf{C} = \exp(\tau \mathbf{A}) \mathbf{B} \).
        class(abstract_linop_rdp), intent(in) :: A
        !! Linear operator to be exponentiated.
        class(abstract_vector_rdp), intent(in) :: B(:)
        !! Input matrix on which to apply \( \exp(\tau \mathbf{A}) \).
        real(dp), intent(in) :: tau
        !! Time horizon for the exponentiation.
        real(dp), intent(in) :: tol
        !! Solution toleance based on error estimates.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: trans
        !! Use transpose ?
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        integer, optional, intent(in) :: kdim
        !! Maximum size of the Krylov subspace.

        ! ----- Internal variables -----
        integer, parameter :: kmax = 100
        integer :: i, j, k, p, kpm, kp, kpp, nk
        ! Block-Arnoldi factorization.
        class(abstract_vector_rdp), allocatable :: X(:)
        real(dp), allocatable :: H(:, :)
        ! Normalization & temporary arrays.
        real(dp), allocatable :: R(:, :), E(:, :), em(:, :)
        integer, allocatable :: perm(:), ptrans(:)
        class(abstract_vector_rdp), allocatable :: Xwrk(:), Cwrk(:)
        real(dp) :: err_est
        character(len=512) :: msg

        ! Optional arguments.
        logical :: transpose
        logical :: verbose
        integer :: nsteps

        ! Determine block size.
        p = size(B)

        ! Deals with the optional args.
        transpose = optval(trans, .false.)
        verbose   = optval(verbosity, .false.)
        nsteps    = optval(kdim, kmax)
        nk        = nsteps*p
        
        info = 0

        ! Allocate arrays.
        allocate(R(p, p)) ; allocate(perm(p)) ; allocate(ptrans(p))
        allocate(X(p*(nk+1)), source=B(1)) ; allocate(H(p*(nk+1), p*(nk+1)))
        allocate(E(p*(nk+1), p*(nk+1))) ; allocate(em(p, p))

        ! Scratch arrays.
        allocate(Xwrk(p), source=B) ; allocate(Cwrk(p), source=B(1))

        ! Normalize input matrix and initialize Krylov subspace.
        R = 0.0_dp
        call qr(Xwrk, R, perm, info) ; call apply_inverse_permutation_matrix(R, perm)

        if (norm2(abs(R)) == 0.0_dp) then
            ! Input matrix is zero.
            do i = 1, size(C)
                call C(i)%zero()
            enddo
            err_est = 0.0_dp ; k = 0 ; kpp = p
        else
            call initialize_krylov_subspace(X, Xwrk) ; H = 0.0_dp

            expm_arnoldi: do k = 1, nk
                ! Set counters.
                kpm = (k-1)*p ; kp = kpm + p ; kpp = kp + p

                ! Reset working arrays.
                E = 0.0_dp
                do i = 1, size(Xwrk)
                    call Xwrk(i)%zero()
                enddo

                ! Compute the k-th step of the Arnoldi factorization.
                call arnoldi(A, X, H, info, kstart=k, kend=k, transpose=transpose, blksize=p)
                call check_info(info, 'arnoldi', module=this_module, procedure='kexpm_mat_rdp')


                if (info == kp) then
                    ! Arnoldi breakdown. Do not consider extended matrix.
                    kpp = kp
                    info = -2
                endif

                ! Compute the (dense) matrix exponential of the extended Hessenberg matrix.
                call expm(E(:kpp, :kpp), tau*H(:kpp, :kpp))

                ! Project back to original space.
                do i = 1, size(Xwrk)
                    call Xwrk(i)%zero()
                    do j = 1, kpp
                        call Xwrk(i)%axpby(one_rdp, X(j), E(j, i))
                    enddo
                enddo

                do i = 1, p
                    call C(i)%zero()
                    do j = 1, p
                        call C(i)%axpby(one_rdp, Xwrk(j), R(j, i))
                    enddo
                enddo

                ! Cheap error estimate.
                if (info == kp) then
                    ! Approximation is exact.
                    err_est = 0.0_dp
                else
                    em = matmul(E(kp+1:kpp, :p), R(:p, :p))
                    err_est = norm2(abs(em))
                endif

                if (err_est <= tol) exit expm_arnoldi

            enddo expm_arnoldi
        endif

        if (err_est .le. tol) then
            info = kpp
            if (p.eq.1) then
               write(msg, *) new_line('A'),' Arnoldi approxmation of the action of the exp. propagator converged!'
            else
               write(msg, *) new_line('A'),' Block Arnoldi approxmation of the action of the exp. propagator converged!'
            endif 
            write(msg, '(A,2(A,I3))') trim(msg)//new_line('A'), 'n° of vectors:', k+1, ' per input vector, total:', kpp
            write(msg, '(A,4X,A,E10.4)') trim(msg)//new_line('A'), 'estimated error:   ||err_est||_2 = ', err_est
            write(msg, '(A,4X,A,E10.4)') trim(msg)//new_line('A'), 'desired tolerance:           tol = ', tol
        else
            info = -1
            if (p.eq.1) then
                write(msg, *) 'Arnoldi approxmation of the action of the exp. propagator did not converge!'
            else
                write(msg, *) 'Block Arnoldi approxmation of the action of the exp. propagator did not converge!'
            endif
            write(msg, '(A,4X,2(A,I3))') trim(msg)//new_line('A'), 'maximum n° of vectors reached: ', nsteps+1, &
                                   & ' per input vector, total:', kpp
            write(msg, '(A,4X,A,E10.4)') trim(msg)//new_line('A'), 'estimated error:   ||err_est||_2 = ', err_est
            write(msg, '(A,4X,A,E10.4)') trim(msg)//new_line('A'), 'desired tolerance:           tol = ', tol
        endif

        if (verbose) call logger%log_message(trim(msg), module=this_module, procedure='kexpm_mat_rdp')

        return
    end subroutine kexpm_mat_rdp

    subroutine kexpm_vec_variable_dt_rdp(c, A, b, tau, dt, tol, info, trans, verbosity, kdim)
        class(abstract_vector_rdp), intent(out) :: c
        !! Best approximation of \( \exp(\tau \mathbf{A}) \mathbf{b} \) in the computed Krylov subspace
        class(abstract_linop_rdp), intent(in) :: A
        !! Linear operator to be exponentiated.
        class(abstract_vector_rdp), intent(in) :: b
        !! Input vector on which to apply \( \exp(\tau \mathbf{A}) \).
        real(dp), intent(in) :: tau
        !! Time horizon for the exponentiation.
        real(dp), intent(inout) :: dt
        !! Current time step.
        real(dp), intent(in) :: tol
        !! Solution tolerance based on error estimates.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: trans
        !! Use transpose?
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        integer, optional, intent(in) :: kdim
        !! Maximum size of the Krylov subspace.

        ! ----- Internal variables -----
        integer, parameter :: kmax = 100
        integer :: k, km, kp, nk, tstep, rstep
        ! Arnoldi factorization.
        class(abstract_vector_rdp), allocatable :: X(:)
        real(dp), allocatable :: H(:, :)
        ! Normaliztion & temporary arrays.
        real(dp), allocatable :: E(:, :)
        class(abstract_vector_rdp), allocatable :: b_
        real(dp) :: err_est, beta, Tend, dt_, tol_, scale, sfact
        character(len=128) :: msg

        ! Optional arguments.
        logical :: transpose
        logical :: verbose
        integer :: nsteps
    
        ! Deals with optional args.
        transpose = optval(trans, .false.)
        verbose   = optval(verbosity, .false.)
        nsteps    = optval(kdim, kmax)
        nk        = nsteps

        info = 0
        scale = 1.5_dp
        sfact = 0.5_dp

        ! Allocate arrays.
        allocate(X(nk+1), b_, source=b)
        allocate(H(nk+1, nk+1), E(nk+1, nk+1))

        ! Normalize input vector and initialize Krylov subspace.
        beta = b%norm()
        if (beta == 0.0_dp) then
            ! Input is zero => Output is zero.
            call c%zero()
        else
            tstep = 1
            Tend = 0.0_dp
            if (verbose) then
               write(msg, '(3(A,E12.6))') 'tau = ', tau, ', dt_0 = ', dt, ', tol = ', tol
               call logger%log_message(trim(msg), module=this_module, procedure='kexpm_vec_variable_dt_rdp')
            end if
            do while ( Tend < tau )

               ! (Re-)initialization of X
               call initialize_krylov_subspace(X) ; call X(1)%add(b_) ;
               beta = b_%norm() ; call X(1)%scal(one_rdp / beta)
               H = 0.0_dp

               ! choose time step
               dt  = min(tau, scale*dt)    ! Reuse step size of last call (and try to increase)
               if (dt > tau-Tend) then
                  dt_ = tau-Tend       ! If a small final step is needed to reach Tend, ignore it on restart
                  if (verbose) then
                     write(msg, '(A,E12.6)') 'Filler step, set dt_ = ', dt_
                     call logger%log_message(trim(msg), module=this_module, procedure='kexpm_vec_variable_dt_rdp')
                  end if
               else
                  dt_ = dt
               end if
               tol_ = dt_/tau*tol*sfact
               if (verbose) then
                  write(msg, '(A,I3,4(A,E8.2))') 'step ', tstep, ': T= ', Tend, ', dt= ', dt, ', dt_= ', dt_, ', tol_= ', tol_
                  call logger%log_message(trim(msg), module=this_module, procedure='kexpm_vec_variable_dt_rdp')
               end if

               expm_arnoldi: do k = 1, nk
                   km = k-1 ; kp = k+1

                   ! Reset work arrays.
                   E = 0.0_dp

                   ! Compute k-th step Arnoldi factorization.
                   call arnoldi(A, X, H, info, kstart=k, kend=k, transpose=transpose)
                   call check_info(info, 'arnoldi', module=this_module, procedure='kexpm_vec_variable_dt_rdp')

                   ! Compute approximation.
                   if (info == k) then
                       ! Arnoldi breakdown, do not consider extended matrix.
                       kp = k
                       info = -2
                   endif

                   ! Compute the (dense) matrix exponential of the extended Hessenberg matrix.
                   call expm(E(:kp, :kp), dt_*H(:kp, :kp))

                   ! Cheap error esimate (this actually is the magnitude of the included correction
                   ! and thus is very conservative).
                   if (info == k) then
                       ! Approximation is exact.
                       err_est = 0.0_dp
                   else
                       err_est = abs(E(kp, 1) * beta)
                   endif

                   ! Check convergence.
                   if (err_est <= tol_) exit expm_arnoldi

               enddo expm_arnoldi

               if (verbose) then
                  write(msg, '(34X,A,E12.6)') 'err_est = ', err_est
                  call logger%log_message(trim(msg), module=this_module, procedure='kexpm_vec_variable_dt_rdp')
               end if

               if (err_est > tol_) then ! decrease dt_ until error estimate is below tolerance
                                       ! NOTE: We do not need to recompute X,H!
                  rstep = 1
                  do while (err_est > tol_)

                     !hfact = (tol/err_est)**(one_rdp/(kdim+1))
                     !dt_ = 0.5_dp*sfact_dp*hfact*dt_ 
                     dt_ = dt_/scale ! this is a quick hack, should choose step better

                     tol_ = dt_/tau*tol*sfact ! update accepted tol based on new time step

                     if (dt_ < atol_dp) then
                        write(msg, *) 'dt < atol_dp. Cannot compute action of the matrix exponential!'
                        call stop_error(trim(msg), module=this_module, procedure='kexpm_vec_variable_dt_rdp')
                     end if

                     ! Compute the (dense) matrix exponential of the extended Hessenberg matrix.
                     call expm(E(:kp, :kp), dt_*H(:kp, :kp))

                     ! Compute new error estimate
                     err_est = abs(E(kp, 1) * beta)

                     if (verbose) then
                        write(msg, '(6X,3(A,E10.4))') 'v dt_ = ', dt_, ', err_est = ', err_est,',  tol_ = ', tol_
                        call logger%log_message(trim(msg), module=this_module, procedure='kexpm_vec_variable_dt_rdp')
                     end if

                     rstep = rstep + 1

                  end do

                  ! save dt for reuse
                  dt = dt_

               end if

               ! Project back into original space.
               block
                  class(abstract_vector_rdp), allocatable :: xwrk
                  call linear_combination(xwrk, X(:kp), E(:kp, 1))
                  call b_%axpby(zero_rdp, xwrk, one_rdp*beta)
               end block

               ! move forward in time
               Tend = Tend + dt_

               tstep = tstep + 1

               if (verbose) then
                  write(msg, '(A,E10.4,2(A,I3),A,E10.4)') '    post: T = ', Tend, ', k = ', k, '/', kdim+1, ', dt => ', dt
                  call logger%log_message(trim(msg), module=this_module, procedure='kexpm_vec_variable_dt_rdp')
               end if
           
            end do ! integration

            ! copy result into output vector
            call c%axpby(zero_rdp, b_, one_rdp)

        endif      

        return
    end subroutine kexpm_vec_variable_dt_rdp

    subroutine kexpm_mat_variable_dt_rdp(C, A, B, tau, dt, tol, info, trans, verbosity, kdim)
        class(abstract_vector_rdp), intent(out) :: C(:)
        !! Best Krylov approximation of \( \mathbf{C} = \exp(\tau \mathbf{A}) \mathbf{B} \).
        class(abstract_linop_rdp), intent(in) :: A
        !! Linear operator to be exponentiated.
        class(abstract_vector_rdp), intent(in) :: B(:)
        !! Input matrix on which to apply \( \exp(\tau \mathbf{A}) \).
        real(dp), intent(in) :: tau
        !! Time horizon for the exponentiation.
        real(dp), intent(inout) :: dt
        !! Current time step.
        real(dp), intent(in) :: tol
        !! Solution toleance based on error estimates.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: trans
        !! Use transpose ?
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        integer, optional, intent(in) :: kdim
        !! Maximum size of the Krylov subspace.

        ! ----- Internal variables -----
        integer, parameter :: kmax = 100
        integer :: i, j, k, p, kpm, kp, kpp, nk, tstep
        ! Block-Arnoldi factorization.
        class(abstract_vector_rdp), allocatable :: X(:)
        real(dp), allocatable :: H(:, :)
        ! Normalization & temporary arrays.
        real(dp), allocatable :: R(:, :), E(:, :), em(:, :)
        integer, allocatable :: perm(:), ptrans(:)
        class(abstract_vector_rdp), allocatable :: Xwrk(:), Cwrk(:), B_(:)
        real(dp) :: err_est, Tend, dt_, tol_, scale
        character(len=128) :: msg

        ! Optional arguments.
        logical :: transpose
        logical :: verbose
        integer :: nsteps

        ! Determine block size.
        p = size(B)

        ! Deals with the optional args.
        transpose = optval(trans, .false.)
        verbose   = optval(verbosity, .false.)
        nsteps    = optval(kdim, kmax)
        nk        = nsteps*p
        
        info = 0
        scale = 1.5

        ! Allocate arrays.
        allocate(perm(p), ptrans(p))
        allocate(H(p*(nk+1), p*(nk+1)), E(p*(nk+1), p*(nk+1)), R(p, p), em(p, p))
        allocate(X(p*(nk+1)), source=B(1))

        ! Scratch arrays.
        allocate(B_(p), source=B) ; allocate(Cwrk(p), source=B(1))

        if (norm2(abs(R)) == 0.0_dp) then
            ! Input matrix is zero.
            call zero_basis(C)
        else
            tstep = 1
            Tend = 0.0_dp
            do while ( Tend < tau )

               ! Normalize B_ and (re-)initialize Krylov subspace.
               R = 0.0_dp
               call qr(B_, R, perm, info) ; call apply_inverse_permutation_matrix(R, perm)
               call initialize_krylov_subspace(X, B_)
               H = 0.0_dp

               ! choose time step
               dt  = min(tau, scale*dt)    ! Reuse step size of last call (and try to increase)
               if (dt > tau-Tend) then
                  dt_ = tau-Tend       ! If a small final step is needed to reach Tend, ignore it on restart
                  if (verbose) then
                     write(msg, '(A,E12.6)') 'Filler step, set dt_ = ', dt_
                     call logger%log_message(trim(msg), module=this_module, procedure='kexpm_mat_variable_dt_rdp')
                  end if
               else
                  dt_ = dt
               end if
               tol_ = dt_/tau*tol
               if (verbose) then
                  write(msg, '(A,I3,4(A,E8.2))') 'step ', tstep, ': T= ', Tend, ', dt= ', dt, ', dt_= ', dt_, ', tol_= ', tol_
                  call logger%log_message(trim(msg), module=this_module, procedure='kexpm_mat_variable_dt_rdp')
               end if

               expm_arnoldi: do k = 1, nk
                   ! Set counters.
                   kpm = (k-1)*p ; kp = kpm + p ; kpp = kp + p

                   ! Reset working arrays.
                   E = 0.0_dp

                   ! Compute the k-th step of the Arnoldi factorization.
                   call arnoldi(A, X, H, info, kstart=k, kend=k, transpose=transpose, blksize=p)
                   call check_info(info, 'arnoldi', module=this_module, procedure='kexpm_mat_variable_dt_rdp')

                   if (info == kp) then
                       ! Arnoldi breakdown. Do not consider extended matrix.
                       kpp = kp
                       info = -2
                   endif

                   ! Compute the (dense) matrix exponential of the extended Hessenberg matrix.
                   call expm(E(:kpp, :kpp), dt_*H(:kpp, :kpp))

                   ! Cheap error estimate.
                   if (info == kp) then
                       ! Approximation is exact.
                       err_est = 0.0_dp
                   else
                       em = matmul(E(kp+1:kpp, :p), R(:p, :p))
                       err_est = norm2(abs(em))
                   endif

                   if (err_est <= tol_) exit expm_arnoldi

               enddo expm_arnoldi

               if (verbose) then
                  write(msg, '(34X,A,E10.4)') 'err_est = ', err_est
                  call logger%log_message(trim(msg), module=this_module, procedure='kexpm_mat_variable_dt_rdp')
               end if

               if (err_est > tol_) then ! decrease dt_ until error estimate is below tolerance
                                       ! NOTE: We do not need to recompute X,H!
                  do while (err_est > tol_)

                     !hfact = (tol/err_est)**(one_rdp/kdim)
                     !dt_ = sfact_dp*hfact*dt_
                     dt_ = dt_/scale ! this is a quick hack, should choose step better

                     tol_ = dt_/tau*tol ! update accepted tol based on new time step

                     if (dt_ < atol_dp) then
                        write(msg, *) 'dt < atol_dp. Cannot compute action of the matrix exponential!'
                        call stop_error(trim(msg), module=this_module, procedure='kexpm_mat_variable_dt_rdp')
                     end if

                     ! Compute the (dense) matrix exponential of the extended Hessenberg matrix.
                     call expm(E(:kpp, :kpp), dt_*H(:kpp, :kpp))

                     ! Compute new error estimate
                     em = matmul(E(kp+1:kpp, :p), R(:p, :p))
                     err_est = norm2(abs(em))

                     if (verbose) then
                        write(msg, '(6X,3(A,E12.6))') 'v dt_ = ', dt_, ', err_est = ', err_est, ', tol_ = ', tol_
                        call logger%log_message(trim(msg), module=this_module, procedure='kexpm_mat_variable_dt_rdp')
                     end if

                  end do

                  ! save dt for reuse
                  dt = dt_
               
               end if

               ! Project back into original space.
               block
                  class(abstract_vector_rdp), allocatable :: Xwrk(:)
                  call linear_combination(Xwrk, X(:kp), matmul(E(:kp, :p), R))
               end block

               ! move forward in time
               Tend = Tend + dt_

               tstep = tstep + 1

               if (verbose) then
                  write(msg, '(A,E10.4,2(A,I3),A,E10.4)') '    post: T = ', Tend, ', k = ', k, '/', kdim+1, ', dt => ', dt
                  call logger%log_message(trim(msg), module=this_module, procedure='kexpm_mat_variable_dt_rdp')
               end if
           
            end do ! integration

            ! copy result into output vector
            call copy_basis(C, B_)
        endif

        return
    end subroutine kexpm_mat_variable_dt_rdp

    subroutine k_exptA_rdp(vec_out, A, vec_in, tau, info, trans)
        class(abstract_vector_rdp), intent(out) :: vec_out
        !! Solution vector.
        class(abstract_linop_rdp), intent(inout) :: A
        !! Linear operator to be exponentiated.
        class(abstract_vector_rdp), intent(in) :: vec_in
        !! Input vector to be multiplied by \( \exp(\tau \mathbf{A}) \).
        real(dp), intent(in) :: tau
        !! Time horizon for the exponentiation.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: trans
        !! Use adjoint ?

        ! ----- Internal variables -----
        real(dp) :: tol
        integer :: kdim
        logical :: verbose

        tol = rtol_dp
        kdim = 30
        verbose = .false.

        call kexpm(vec_out, A, vec_in, tau, tol, info, trans=trans, verbosity=verbose, kdim=kdim)
        call check_info(info, 'kexpm', module=this_module, procedure='k_exptA_rdp')

        return
    end subroutine k_exptA_rdp

    subroutine expm_csp(E, A, order)
        complex(sp), intent(in) :: A(:, :)
        !! Matrix to be exponentiated.
        complex(sp), intent(out) :: E(:, :)
        !! Output matrix E = exp(tA).
        integer, optional :: order
        !! Order of the Pade approximation.

        !----- Internal variables -----
        complex(sp), allocatable :: A2(:, :), Q(:, :), X(:, :), invQ(:, :), wrk(:)
        real(sp) :: a_norm, c
        integer :: n, ee, k, s
        logical :: p
        integer :: p_order

        ! Deal with optional args.
        p_order = optval(order, 10)

        n = size(A, 1)

        ! Allocate arrays.
        allocate(A2(n, n)) ; allocate(X(n, n))
        allocate(Q(n, n)) ; allocate(invQ(n, n))
        allocate(wrk(n))

        ! Compute the L-infinity norm.
        a_norm = norml_csp(A)

        ! Determine scaling factor for the matrix.
        ee = int(log2_rsp(a_norm)) + 1
        s = max(0, ee+1)

        ! Scale the input matrix & initialize polynomial.
        A2 = A / 2.0_sp**s
        X = A2

        ! Initialize P & Q and add first step.
        c = 0.5_sp
        E = eye(n) ; E = E + c*A2

        Q = eye(n) ; Q = Q - c*A2

        ! Iteratively compute the Pade approximation.
        p = .true.
        do k = 2, p_order
            c = c*(p_order - k + 1) / (k * (2*p_order - k + 1))
            X = matmul(A2, X)
            E = E + c*X
            if (p) then
                Q = Q + c*X
            else
                Q = Q - c*X
            endif
            p = .not. p
        enddo

        invQ = Q ; call inv(invQ)
        E = matmul(invQ, E)

        do k = 1, s
            E = matmul(E, E)
        enddo

        return
    end subroutine expm_csp

    subroutine kexpm_vec_csp(c, A, b, tau, tol, info, trans, verbosity, kdim)
        class(abstract_vector_csp), intent(out) :: c
        !! Best approximation of \( \exp(\tau \mathbf{A}) \mathbf{b} \) in the computed Krylov subspace
        class(abstract_linop_csp), intent(in) :: A
        !! Linear operator to be exponentiated.
        class(abstract_vector_csp), intent(in) :: b
        !! Input vector on which to apply \( \exp(\tau \mathbf{A}) \).
        real(sp), intent(in) :: tau
        !! Time horizon for the exponentiation.
        real(sp), intent(in) :: tol
        !! Solution tolerance based on error estimates.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: trans
        !! Use transpose?
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        integer, optional, intent(in) :: kdim
        !! Maximum size of the Krylov subspace.

        ! ----- Internal variables -----
        integer, parameter :: kmax = 100
        integer :: k, km, kp, nk
        ! Arnoldi factorization.
        class(abstract_vector_csp), allocatable :: X(:)
        complex(sp), allocatable :: H(:, :)
        ! Normaliztion & temporary arrays.
        complex(sp), allocatable :: E(:, :)
        class(abstract_vector_csp), allocatable :: Xwrk
        real(sp) :: err_est, beta
        character(len=512) :: msg

        ! Optional arguments.
        logical :: transpose
        logical :: verbose
        integer :: nsteps
    
        ! Deals with optional args.
        transpose = optval(trans, .false.)
        verbose   = optval(verbosity, .false.)
        nsteps    = optval(kdim, kmax)
        nk        = nsteps

        info = 0

        ! Allocate arrays.
        allocate(X(nk+1), source=b) ; allocate(H(nk+1, nk+1))
        allocate(E(nk+1, nk+1)) ; allocate(Xwrk, source=b)

        ! Normalize input vector and initialize Krylov subspace.
        beta = b%norm()
        if (beta == 0.0_sp) then
            ! Input is zero => Output is zero.
            call c%zero()
            err_est = 0.0_sp
            kp = 1
        else
            call initialize_krylov_subspace(X)
            call X(1)%add(b) ; call X(1)%scal(one_csp / beta)
            H = 0.0_sp

            expm_arnoldi: do k = 1, nk
                km = k-1 ; kp = k+1

                ! Reset work arrays.
                E = 0.0_sp

                ! Compute k-th step Arnoldi factorization.
                call arnoldi(A, X, H, info, kstart=k, kend=k, transpose=transpose)
                call check_info(info, 'arnoldi', module=this_module, procedure='kexpm_vec_csp')

                ! Compute approximation.
                if (info == k) then
                    ! Arnoldi breakdown, do not consider extended matrix.
                    kp = k
                    info = -2
                endif

                ! Compute the (dense) matrix exponential of the extended Hessenberg matrix.
                call expm(E(:kp, :kp), tau*H(:kp, :kp))

                ! Project back into original space.
                call linear_combination(Xwrk, X(:kp), E(:kp, 1))
                call c%axpby(zero_csp, Xwrk, one_csp*beta)

                ! Cheap error esimate (this actually is the magnitude of the included correction
                ! and thus is very conservative).
                if (info == k) then
                    ! Approximation is exact.
                    err_est = 0.0_sp
                else
                    err_est = abs(E(kp, 1) * beta)
                endif

                ! Check convergence.
                if (err_est <= tol) exit expm_arnoldi

            enddo expm_arnoldi
        endif

        if (err_est <= tol) then
            info = kp
            write(msg, *) new_line('A'),' Arnoldi approxmation of the action of the exp. propagator converged!'
            write(msg, '(A,4X,A,I3)') trim(msg)//new_line('A'), 'n° of vectors:', kp
            write(msg, '(A,4X,A,E10.4)') trim(msg)//new_line('A'), 'estimated error:   ||err_est||_2 = ', err_est
            write(msg, '(A,4X,A,E10.4)') trim(msg)//new_line('A'), 'desired tolerance:           tol = ', tol
        else
            info = -1
            write(msg, *) new_line('A'),' Arnoldi approxmation of the action of the exp. propagator did not converge!'
            write(msg, '(A,4X,A,I3)') trim(msg)//new_line('A'), 'Maximum n° of vectors reached: ', nk + 1
            write(msg, '(A,4X,A,E10.4)') trim(msg)//new_line('A'), 'Estimated error              :   ||err_est||_2 = ', err_est
            write(msg, '(A,4X,A,E10.4)') trim(msg)//new_line('A'), 'Desired tolerance            :           tol = ', tol
        endif
        if (verbose) call logger%log_message(trim(msg), module=this_module, procedure='kexpm_vec_csp')

        return
    end subroutine kexpm_vec_csp

    subroutine kexpm_mat_csp(C, A, B, tau, tol, info, trans, verbosity, kdim)
        class(abstract_vector_csp), intent(out) :: C(:)
        !! Best Krylov approximation of \( \mathbf{C} = \exp(\tau \mathbf{A}) \mathbf{B} \).
        class(abstract_linop_csp), intent(in) :: A
        !! Linear operator to be exponentiated.
        class(abstract_vector_csp), intent(in) :: B(:)
        !! Input matrix on which to apply \( \exp(\tau \mathbf{A}) \).
        real(sp), intent(in) :: tau
        !! Time horizon for the exponentiation.
        real(sp), intent(in) :: tol
        !! Solution toleance based on error estimates.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: trans
        !! Use transpose ?
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        integer, optional, intent(in) :: kdim
        !! Maximum size of the Krylov subspace.

        ! ----- Internal variables -----
        integer, parameter :: kmax = 100
        integer :: i, j, k, p, kpm, kp, kpp, nk
        ! Block-Arnoldi factorization.
        class(abstract_vector_csp), allocatable :: X(:)
        complex(sp), allocatable :: H(:, :)
        ! Normalization & temporary arrays.
        complex(sp), allocatable :: R(:, :), E(:, :), em(:, :)
        integer, allocatable :: perm(:), ptrans(:)
        class(abstract_vector_csp), allocatable :: Xwrk(:), Cwrk(:)
        real(sp) :: err_est
        character(len=512) :: msg

        ! Optional arguments.
        logical :: transpose
        logical :: verbose
        integer :: nsteps

        ! Determine block size.
        p = size(B)

        ! Deals with the optional args.
        transpose = optval(trans, .false.)
        verbose   = optval(verbosity, .false.)
        nsteps    = optval(kdim, kmax)
        nk        = nsteps*p
        
        info = 0

        ! Allocate arrays.
        allocate(R(p, p)) ; allocate(perm(p)) ; allocate(ptrans(p))
        allocate(X(p*(nk+1)), source=B(1)) ; allocate(H(p*(nk+1), p*(nk+1)))
        allocate(E(p*(nk+1), p*(nk+1))) ; allocate(em(p, p))

        ! Scratch arrays.
        allocate(Xwrk(p), source=B) ; allocate(Cwrk(p), source=B(1))

        ! Normalize input matrix and initialize Krylov subspace.
        R = 0.0_sp
        call qr(Xwrk, R, perm, info) ; call apply_inverse_permutation_matrix(R, perm)

        if (norm2(abs(R)) == 0.0_sp) then
            ! Input matrix is zero.
            do i = 1, size(C)
                call C(i)%zero()
            enddo
            err_est = 0.0_sp ; k = 0 ; kpp = p
        else
            call initialize_krylov_subspace(X, Xwrk) ; H = 0.0_sp

            expm_arnoldi: do k = 1, nk
                ! Set counters.
                kpm = (k-1)*p ; kp = kpm + p ; kpp = kp + p

                ! Reset working arrays.
                E = 0.0_sp
                do i = 1, size(Xwrk)
                    call Xwrk(i)%zero()
                enddo

                ! Compute the k-th step of the Arnoldi factorization.
                call arnoldi(A, X, H, info, kstart=k, kend=k, transpose=transpose, blksize=p)
                call check_info(info, 'arnoldi', module=this_module, procedure='kexpm_mat_csp')


                if (info == kp) then
                    ! Arnoldi breakdown. Do not consider extended matrix.
                    kpp = kp
                    info = -2
                endif

                ! Compute the (dense) matrix exponential of the extended Hessenberg matrix.
                call expm(E(:kpp, :kpp), tau*H(:kpp, :kpp))

                ! Project back to original space.
                do i = 1, size(Xwrk)
                    call Xwrk(i)%zero()
                    do j = 1, kpp
                        call Xwrk(i)%axpby(one_csp, X(j), E(j, i))
                    enddo
                enddo

                do i = 1, p
                    call C(i)%zero()
                    do j = 1, p
                        call C(i)%axpby(one_csp, Xwrk(j), R(j, i))
                    enddo
                enddo

                ! Cheap error estimate.
                if (info == kp) then
                    ! Approximation is exact.
                    err_est = 0.0_sp
                else
                    em = matmul(E(kp+1:kpp, :p), R(:p, :p))
                    err_est = norm2(abs(em))
                endif

                if (err_est <= tol) exit expm_arnoldi

            enddo expm_arnoldi
        endif

        if (err_est .le. tol) then
            info = kpp
            if (p.eq.1) then
               write(msg, *) new_line('A'),' Arnoldi approxmation of the action of the exp. propagator converged!'
            else
               write(msg, *) new_line('A'),' Block Arnoldi approxmation of the action of the exp. propagator converged!'
            endif 
            write(msg, '(A,2(A,I3))') trim(msg)//new_line('A'), 'n° of vectors:', k+1, ' per input vector, total:', kpp
            write(msg, '(A,4X,A,E10.4)') trim(msg)//new_line('A'), 'estimated error:   ||err_est||_2 = ', err_est
            write(msg, '(A,4X,A,E10.4)') trim(msg)//new_line('A'), 'desired tolerance:           tol = ', tol
        else
            info = -1
            if (p.eq.1) then
                write(msg, *) 'Arnoldi approxmation of the action of the exp. propagator did not converge!'
            else
                write(msg, *) 'Block Arnoldi approxmation of the action of the exp. propagator did not converge!'
            endif
            write(msg, '(A,4X,2(A,I3))') trim(msg)//new_line('A'), 'maximum n° of vectors reached: ', nsteps+1, &
                                   & ' per input vector, total:', kpp
            write(msg, '(A,4X,A,E10.4)') trim(msg)//new_line('A'), 'estimated error:   ||err_est||_2 = ', err_est
            write(msg, '(A,4X,A,E10.4)') trim(msg)//new_line('A'), 'desired tolerance:           tol = ', tol
        endif

        if (verbose) call logger%log_message(trim(msg), module=this_module, procedure='kexpm_mat_csp')

        return
    end subroutine kexpm_mat_csp

    subroutine kexpm_vec_variable_dt_csp(c, A, b, tau, dt, tol, info, trans, verbosity, kdim)
        class(abstract_vector_csp), intent(out) :: c
        !! Best approximation of \( \exp(\tau \mathbf{A}) \mathbf{b} \) in the computed Krylov subspace
        class(abstract_linop_csp), intent(in) :: A
        !! Linear operator to be exponentiated.
        class(abstract_vector_csp), intent(in) :: b
        !! Input vector on which to apply \( \exp(\tau \mathbf{A}) \).
        real(sp), intent(in) :: tau
        !! Time horizon for the exponentiation.
        real(sp), intent(inout) :: dt
        !! Current time step.
        real(sp), intent(in) :: tol
        !! Solution tolerance based on error estimates.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: trans
        !! Use transpose?
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        integer, optional, intent(in) :: kdim
        !! Maximum size of the Krylov subspace.

        ! ----- Internal variables -----
        integer, parameter :: kmax = 100
        integer :: k, km, kp, nk, tstep, rstep
        ! Arnoldi factorization.
        class(abstract_vector_csp), allocatable :: X(:)
        complex(sp), allocatable :: H(:, :)
        ! Normaliztion & temporary arrays.
        complex(sp), allocatable :: E(:, :)
        class(abstract_vector_csp), allocatable :: b_
        real(sp) :: err_est, beta, Tend, dt_, tol_, scale, sfact
        character(len=128) :: msg

        ! Optional arguments.
        logical :: transpose
        logical :: verbose
        integer :: nsteps
    
        ! Deals with optional args.
        transpose = optval(trans, .false.)
        verbose   = optval(verbosity, .false.)
        nsteps    = optval(kdim, kmax)
        nk        = nsteps

        info = 0
        scale = 1.5_sp
        sfact = 0.5_sp

        ! Allocate arrays.
        allocate(X(nk+1), b_, source=b)
        allocate(H(nk+1, nk+1), E(nk+1, nk+1))

        ! Normalize input vector and initialize Krylov subspace.
        beta = b%norm()
        if (beta == 0.0_sp) then
            ! Input is zero => Output is zero.
            call c%zero()
        else
            tstep = 1
            Tend = 0.0_sp
            if (verbose) then
               write(msg, '(3(A,E12.6))') 'tau = ', tau, ', dt_0 = ', dt, ', tol = ', tol
               call logger%log_message(trim(msg), module=this_module, procedure='kexpm_vec_variable_dt_csp')
            end if
            do while ( Tend < tau )

               ! (Re-)initialization of X
               call initialize_krylov_subspace(X) ; call X(1)%add(b_) ;
               beta = b_%norm() ; call X(1)%scal(one_csp / beta)
               H = 0.0_sp

               ! choose time step
               dt  = min(tau, scale*dt)    ! Reuse step size of last call (and try to increase)
               if (dt > tau-Tend) then
                  dt_ = tau-Tend       ! If a small final step is needed to reach Tend, ignore it on restart
                  if (verbose) then
                     write(msg, '(A,E12.6)') 'Filler step, set dt_ = ', dt_
                     call logger%log_message(trim(msg), module=this_module, procedure='kexpm_vec_variable_dt_csp')
                  end if
               else
                  dt_ = dt
               end if
               tol_ = dt_/tau*tol*sfact
               if (verbose) then
                  write(msg, '(A,I3,4(A,E8.2))') 'step ', tstep, ': T= ', Tend, ', dt= ', dt, ', dt_= ', dt_, ', tol_= ', tol_
                  call logger%log_message(trim(msg), module=this_module, procedure='kexpm_vec_variable_dt_csp')
               end if

               expm_arnoldi: do k = 1, nk
                   km = k-1 ; kp = k+1

                   ! Reset work arrays.
                   E = 0.0_sp

                   ! Compute k-th step Arnoldi factorization.
                   call arnoldi(A, X, H, info, kstart=k, kend=k, transpose=transpose)
                   call check_info(info, 'arnoldi', module=this_module, procedure='kexpm_vec_variable_dt_csp')

                   ! Compute approximation.
                   if (info == k) then
                       ! Arnoldi breakdown, do not consider extended matrix.
                       kp = k
                       info = -2
                   endif

                   ! Compute the (dense) matrix exponential of the extended Hessenberg matrix.
                   call expm(E(:kp, :kp), dt_*H(:kp, :kp))

                   ! Cheap error esimate (this actually is the magnitude of the included correction
                   ! and thus is very conservative).
                   if (info == k) then
                       ! Approximation is exact.
                       err_est = 0.0_sp
                   else
                       err_est = abs(E(kp, 1) * beta)
                   endif

                   ! Check convergence.
                   if (err_est <= tol_) exit expm_arnoldi

               enddo expm_arnoldi

               if (verbose) then
                  write(msg, '(34X,A,E12.6)') 'err_est = ', err_est
                  call logger%log_message(trim(msg), module=this_module, procedure='kexpm_vec_variable_dt_csp')
               end if

               if (err_est > tol_) then ! decrease dt_ until error estimate is below tolerance
                                       ! NOTE: We do not need to recompute X,H!
                  rstep = 1
                  do while (err_est > tol_)

                     !hfact = (tol/err_est)**(one_csp/(kdim+1))
                     !dt_ = 0.5_sp*sfact_sp*hfact*dt_ 
                     dt_ = dt_/scale ! this is a quick hack, should choose step better

                     tol_ = dt_/tau*tol*sfact ! update accepted tol based on new time step

                     if (dt_ < atol_sp) then
                        write(msg, *) 'dt < atol_sp. Cannot compute action of the matrix exponential!'
                        call stop_error(trim(msg), module=this_module, procedure='kexpm_vec_variable_dt_csp')
                     end if

                     ! Compute the (dense) matrix exponential of the extended Hessenberg matrix.
                     call expm(E(:kp, :kp), dt_*H(:kp, :kp))

                     ! Compute new error estimate
                     err_est = abs(E(kp, 1) * beta)

                     if (verbose) then
                        write(msg, '(6X,3(A,E10.4))') 'v dt_ = ', dt_, ', err_est = ', err_est,',  tol_ = ', tol_
                        call logger%log_message(trim(msg), module=this_module, procedure='kexpm_vec_variable_dt_csp')
                     end if

                     rstep = rstep + 1

                  end do

                  ! save dt for reuse
                  dt = dt_

               end if

               ! Project back into original space.
               block
                  class(abstract_vector_csp), allocatable :: xwrk
                  call linear_combination(xwrk, X(:kp), E(:kp, 1))
                  call b_%axpby(zero_csp, xwrk, one_csp*beta)
               end block

               ! move forward in time
               Tend = Tend + dt_

               tstep = tstep + 1

               if (verbose) then
                  write(msg, '(A,E10.4,2(A,I3),A,E10.4)') '    post: T = ', Tend, ', k = ', k, '/', kdim+1, ', dt => ', dt
                  call logger%log_message(trim(msg), module=this_module, procedure='kexpm_vec_variable_dt_csp')
               end if
           
            end do ! integration

            ! copy result into output vector
            call c%axpby(zero_csp, b_, one_csp)

        endif      

        return
    end subroutine kexpm_vec_variable_dt_csp

    subroutine kexpm_mat_variable_dt_csp(C, A, B, tau, dt, tol, info, trans, verbosity, kdim)
        class(abstract_vector_csp), intent(out) :: C(:)
        !! Best Krylov approximation of \( \mathbf{C} = \exp(\tau \mathbf{A}) \mathbf{B} \).
        class(abstract_linop_csp), intent(in) :: A
        !! Linear operator to be exponentiated.
        class(abstract_vector_csp), intent(in) :: B(:)
        !! Input matrix on which to apply \( \exp(\tau \mathbf{A}) \).
        real(sp), intent(in) :: tau
        !! Time horizon for the exponentiation.
        real(sp), intent(inout) :: dt
        !! Current time step.
        real(sp), intent(in) :: tol
        !! Solution toleance based on error estimates.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: trans
        !! Use transpose ?
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        integer, optional, intent(in) :: kdim
        !! Maximum size of the Krylov subspace.

        ! ----- Internal variables -----
        integer, parameter :: kmax = 100
        integer :: i, j, k, p, kpm, kp, kpp, nk, tstep
        ! Block-Arnoldi factorization.
        class(abstract_vector_csp), allocatable :: X(:)
        complex(sp), allocatable :: H(:, :)
        ! Normalization & temporary arrays.
        complex(sp), allocatable :: R(:, :), E(:, :), em(:, :)
        integer, allocatable :: perm(:), ptrans(:)
        class(abstract_vector_csp), allocatable :: Xwrk(:), Cwrk(:), B_(:)
        real(sp) :: err_est, Tend, dt_, tol_, scale
        character(len=128) :: msg

        ! Optional arguments.
        logical :: transpose
        logical :: verbose
        integer :: nsteps

        ! Determine block size.
        p = size(B)

        ! Deals with the optional args.
        transpose = optval(trans, .false.)
        verbose   = optval(verbosity, .false.)
        nsteps    = optval(kdim, kmax)
        nk        = nsteps*p
        
        info = 0
        scale = 1.5

        ! Allocate arrays.
        allocate(perm(p), ptrans(p))
        allocate(H(p*(nk+1), p*(nk+1)), E(p*(nk+1), p*(nk+1)), R(p, p), em(p, p))
        allocate(X(p*(nk+1)), source=B(1))

        ! Scratch arrays.
        allocate(B_(p), source=B) ; allocate(Cwrk(p), source=B(1))

        if (norm2(abs(R)) == 0.0_sp) then
            ! Input matrix is zero.
            call zero_basis(C)
        else
            tstep = 1
            Tend = 0.0_sp
            do while ( Tend < tau )

               ! Normalize B_ and (re-)initialize Krylov subspace.
               R = 0.0_sp
               call qr(B_, R, perm, info) ; call apply_inverse_permutation_matrix(R, perm)
               call initialize_krylov_subspace(X, B_)
               H = 0.0_sp

               ! choose time step
               dt  = min(tau, scale*dt)    ! Reuse step size of last call (and try to increase)
               if (dt > tau-Tend) then
                  dt_ = tau-Tend       ! If a small final step is needed to reach Tend, ignore it on restart
                  if (verbose) then
                     write(msg, '(A,E12.6)') 'Filler step, set dt_ = ', dt_
                     call logger%log_message(trim(msg), module=this_module, procedure='kexpm_mat_variable_dt_csp')
                  end if
               else
                  dt_ = dt
               end if
               tol_ = dt_/tau*tol
               if (verbose) then
                  write(msg, '(A,I3,4(A,E8.2))') 'step ', tstep, ': T= ', Tend, ', dt= ', dt, ', dt_= ', dt_, ', tol_= ', tol_
                  call logger%log_message(trim(msg), module=this_module, procedure='kexpm_mat_variable_dt_csp')
               end if

               expm_arnoldi: do k = 1, nk
                   ! Set counters.
                   kpm = (k-1)*p ; kp = kpm + p ; kpp = kp + p

                   ! Reset working arrays.
                   E = 0.0_sp

                   ! Compute the k-th step of the Arnoldi factorization.
                   call arnoldi(A, X, H, info, kstart=k, kend=k, transpose=transpose, blksize=p)
                   call check_info(info, 'arnoldi', module=this_module, procedure='kexpm_mat_variable_dt_csp')

                   if (info == kp) then
                       ! Arnoldi breakdown. Do not consider extended matrix.
                       kpp = kp
                       info = -2
                   endif

                   ! Compute the (dense) matrix exponential of the extended Hessenberg matrix.
                   call expm(E(:kpp, :kpp), dt_*H(:kpp, :kpp))

                   ! Cheap error estimate.
                   if (info == kp) then
                       ! Approximation is exact.
                       err_est = 0.0_sp
                   else
                       em = matmul(E(kp+1:kpp, :p), R(:p, :p))
                       err_est = norm2(abs(em))
                   endif

                   if (err_est <= tol_) exit expm_arnoldi

               enddo expm_arnoldi

               if (verbose) then
                  write(msg, '(34X,A,E10.4)') 'err_est = ', err_est
                  call logger%log_message(trim(msg), module=this_module, procedure='kexpm_mat_variable_dt_csp')
               end if

               if (err_est > tol_) then ! decrease dt_ until error estimate is below tolerance
                                       ! NOTE: We do not need to recompute X,H!
                  do while (err_est > tol_)

                     !hfact = (tol/err_est)**(one_csp/kdim)
                     !dt_ = sfact_sp*hfact*dt_
                     dt_ = dt_/scale ! this is a quick hack, should choose step better

                     tol_ = dt_/tau*tol ! update accepted tol based on new time step

                     if (dt_ < atol_sp) then
                        write(msg, *) 'dt < atol_sp. Cannot compute action of the matrix exponential!'
                        call stop_error(trim(msg), module=this_module, procedure='kexpm_mat_variable_dt_csp')
                     end if

                     ! Compute the (dense) matrix exponential of the extended Hessenberg matrix.
                     call expm(E(:kpp, :kpp), dt_*H(:kpp, :kpp))

                     ! Compute new error estimate
                     em = matmul(E(kp+1:kpp, :p), R(:p, :p))
                     err_est = norm2(abs(em))

                     if (verbose) then
                        write(msg, '(6X,3(A,E12.6))') 'v dt_ = ', dt_, ', err_est = ', err_est, ', tol_ = ', tol_
                        call logger%log_message(trim(msg), module=this_module, procedure='kexpm_mat_variable_dt_csp')
                     end if

                  end do

                  ! save dt for reuse
                  dt = dt_
               
               end if

               ! Project back into original space.
               block
                  class(abstract_vector_csp), allocatable :: Xwrk(:)
                  call linear_combination(Xwrk, X(:kp), matmul(E(:kp, :p), R))
               end block

               ! move forward in time
               Tend = Tend + dt_

               tstep = tstep + 1

               if (verbose) then
                  write(msg, '(A,E10.4,2(A,I3),A,E10.4)') '    post: T = ', Tend, ', k = ', k, '/', kdim+1, ', dt => ', dt
                  call logger%log_message(trim(msg), module=this_module, procedure='kexpm_mat_variable_dt_csp')
               end if
           
            end do ! integration

            ! copy result into output vector
            call copy_basis(C, B_)
        endif

        return
    end subroutine kexpm_mat_variable_dt_csp

    subroutine k_exptA_csp(vec_out, A, vec_in, tau, info, trans)
        class(abstract_vector_csp), intent(out) :: vec_out
        !! Solution vector.
        class(abstract_linop_csp), intent(inout) :: A
        !! Linear operator to be exponentiated.
        class(abstract_vector_csp), intent(in) :: vec_in
        !! Input vector to be multiplied by \( \exp(\tau \mathbf{A}) \).
        real(sp), intent(in) :: tau
        !! Time horizon for the exponentiation.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: trans
        !! Use adjoint ?

        ! ----- Internal variables -----
        real(sp) :: tol
        integer :: kdim
        logical :: verbose

        tol = rtol_sp
        kdim = 30
        verbose = .false.

        call kexpm(vec_out, A, vec_in, tau, tol, info, trans=trans, verbosity=verbose, kdim=kdim)
        call check_info(info, 'kexpm', module=this_module, procedure='k_exptA_csp')

        return
    end subroutine k_exptA_csp

    subroutine expm_cdp(E, A, order)
        complex(dp), intent(in) :: A(:, :)
        !! Matrix to be exponentiated.
        complex(dp), intent(out) :: E(:, :)
        !! Output matrix E = exp(tA).
        integer, optional :: order
        !! Order of the Pade approximation.

        !----- Internal variables -----
        complex(dp), allocatable :: A2(:, :), Q(:, :), X(:, :), invQ(:, :), wrk(:)
        real(dp) :: a_norm, c
        integer :: n, ee, k, s
        logical :: p
        integer :: p_order

        ! Deal with optional args.
        p_order = optval(order, 10)

        n = size(A, 1)

        ! Allocate arrays.
        allocate(A2(n, n)) ; allocate(X(n, n))
        allocate(Q(n, n)) ; allocate(invQ(n, n))
        allocate(wrk(n))

        ! Compute the L-infinity norm.
        a_norm = norml_cdp(A)

        ! Determine scaling factor for the matrix.
        ee = int(log2_rdp(a_norm)) + 1
        s = max(0, ee+1)

        ! Scale the input matrix & initialize polynomial.
        A2 = A / 2.0_dp**s
        X = A2

        ! Initialize P & Q and add first step.
        c = 0.5_dp
        E = eye(n) ; E = E + c*A2

        Q = eye(n) ; Q = Q - c*A2

        ! Iteratively compute the Pade approximation.
        p = .true.
        do k = 2, p_order
            c = c*(p_order - k + 1) / (k * (2*p_order - k + 1))
            X = matmul(A2, X)
            E = E + c*X
            if (p) then
                Q = Q + c*X
            else
                Q = Q - c*X
            endif
            p = .not. p
        enddo

        invQ = Q ; call inv(invQ)
        E = matmul(invQ, E)

        do k = 1, s
            E = matmul(E, E)
        enddo

        return
    end subroutine expm_cdp

    subroutine kexpm_vec_cdp(c, A, b, tau, tol, info, trans, verbosity, kdim)
        class(abstract_vector_cdp), intent(out) :: c
        !! Best approximation of \( \exp(\tau \mathbf{A}) \mathbf{b} \) in the computed Krylov subspace
        class(abstract_linop_cdp), intent(in) :: A
        !! Linear operator to be exponentiated.
        class(abstract_vector_cdp), intent(in) :: b
        !! Input vector on which to apply \( \exp(\tau \mathbf{A}) \).
        real(dp), intent(in) :: tau
        !! Time horizon for the exponentiation.
        real(dp), intent(in) :: tol
        !! Solution tolerance based on error estimates.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: trans
        !! Use transpose?
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        integer, optional, intent(in) :: kdim
        !! Maximum size of the Krylov subspace.

        ! ----- Internal variables -----
        integer, parameter :: kmax = 100
        integer :: k, km, kp, nk
        ! Arnoldi factorization.
        class(abstract_vector_cdp), allocatable :: X(:)
        complex(dp), allocatable :: H(:, :)
        ! Normaliztion & temporary arrays.
        complex(dp), allocatable :: E(:, :)
        class(abstract_vector_cdp), allocatable :: Xwrk
        real(dp) :: err_est, beta
        character(len=512) :: msg

        ! Optional arguments.
        logical :: transpose
        logical :: verbose
        integer :: nsteps
    
        ! Deals with optional args.
        transpose = optval(trans, .false.)
        verbose   = optval(verbosity, .false.)
        nsteps    = optval(kdim, kmax)
        nk        = nsteps

        info = 0

        ! Allocate arrays.
        allocate(X(nk+1), source=b) ; allocate(H(nk+1, nk+1))
        allocate(E(nk+1, nk+1)) ; allocate(Xwrk, source=b)

        ! Normalize input vector and initialize Krylov subspace.
        beta = b%norm()
        if (beta == 0.0_dp) then
            ! Input is zero => Output is zero.
            call c%zero()
            err_est = 0.0_dp
            kp = 1
        else
            call initialize_krylov_subspace(X)
            call X(1)%add(b) ; call X(1)%scal(one_cdp / beta)
            H = 0.0_dp

            expm_arnoldi: do k = 1, nk
                km = k-1 ; kp = k+1

                ! Reset work arrays.
                E = 0.0_dp

                ! Compute k-th step Arnoldi factorization.
                call arnoldi(A, X, H, info, kstart=k, kend=k, transpose=transpose)
                call check_info(info, 'arnoldi', module=this_module, procedure='kexpm_vec_cdp')

                ! Compute approximation.
                if (info == k) then
                    ! Arnoldi breakdown, do not consider extended matrix.
                    kp = k
                    info = -2
                endif

                ! Compute the (dense) matrix exponential of the extended Hessenberg matrix.
                call expm(E(:kp, :kp), tau*H(:kp, :kp))

                ! Project back into original space.
                call linear_combination(Xwrk, X(:kp), E(:kp, 1))
                call c%axpby(zero_cdp, Xwrk, one_cdp*beta)

                ! Cheap error esimate (this actually is the magnitude of the included correction
                ! and thus is very conservative).
                if (info == k) then
                    ! Approximation is exact.
                    err_est = 0.0_dp
                else
                    err_est = abs(E(kp, 1) * beta)
                endif

                ! Check convergence.
                if (err_est <= tol) exit expm_arnoldi

            enddo expm_arnoldi
        endif

        if (err_est <= tol) then
            info = kp
            write(msg, *) new_line('A'),' Arnoldi approxmation of the action of the exp. propagator converged!'
            write(msg, '(A,4X,A,I3)') trim(msg)//new_line('A'), 'n° of vectors:', kp
            write(msg, '(A,4X,A,E10.4)') trim(msg)//new_line('A'), 'estimated error:   ||err_est||_2 = ', err_est
            write(msg, '(A,4X,A,E10.4)') trim(msg)//new_line('A'), 'desired tolerance:           tol = ', tol
        else
            info = -1
            write(msg, *) new_line('A'),' Arnoldi approxmation of the action of the exp. propagator did not converge!'
            write(msg, '(A,4X,A,I3)') trim(msg)//new_line('A'), 'Maximum n° of vectors reached: ', nk + 1
            write(msg, '(A,4X,A,E10.4)') trim(msg)//new_line('A'), 'Estimated error              :   ||err_est||_2 = ', err_est
            write(msg, '(A,4X,A,E10.4)') trim(msg)//new_line('A'), 'Desired tolerance            :           tol = ', tol
        endif
        if (verbose) call logger%log_message(trim(msg), module=this_module, procedure='kexpm_vec_cdp')

        return
    end subroutine kexpm_vec_cdp

    subroutine kexpm_mat_cdp(C, A, B, tau, tol, info, trans, verbosity, kdim)
        class(abstract_vector_cdp), intent(out) :: C(:)
        !! Best Krylov approximation of \( \mathbf{C} = \exp(\tau \mathbf{A}) \mathbf{B} \).
        class(abstract_linop_cdp), intent(in) :: A
        !! Linear operator to be exponentiated.
        class(abstract_vector_cdp), intent(in) :: B(:)
        !! Input matrix on which to apply \( \exp(\tau \mathbf{A}) \).
        real(dp), intent(in) :: tau
        !! Time horizon for the exponentiation.
        real(dp), intent(in) :: tol
        !! Solution toleance based on error estimates.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: trans
        !! Use transpose ?
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        integer, optional, intent(in) :: kdim
        !! Maximum size of the Krylov subspace.

        ! ----- Internal variables -----
        integer, parameter :: kmax = 100
        integer :: i, j, k, p, kpm, kp, kpp, nk
        ! Block-Arnoldi factorization.
        class(abstract_vector_cdp), allocatable :: X(:)
        complex(dp), allocatable :: H(:, :)
        ! Normalization & temporary arrays.
        complex(dp), allocatable :: R(:, :), E(:, :), em(:, :)
        integer, allocatable :: perm(:), ptrans(:)
        class(abstract_vector_cdp), allocatable :: Xwrk(:), Cwrk(:)
        real(dp) :: err_est
        character(len=512) :: msg

        ! Optional arguments.
        logical :: transpose
        logical :: verbose
        integer :: nsteps

        ! Determine block size.
        p = size(B)

        ! Deals with the optional args.
        transpose = optval(trans, .false.)
        verbose   = optval(verbosity, .false.)
        nsteps    = optval(kdim, kmax)
        nk        = nsteps*p
        
        info = 0

        ! Allocate arrays.
        allocate(R(p, p)) ; allocate(perm(p)) ; allocate(ptrans(p))
        allocate(X(p*(nk+1)), source=B(1)) ; allocate(H(p*(nk+1), p*(nk+1)))
        allocate(E(p*(nk+1), p*(nk+1))) ; allocate(em(p, p))

        ! Scratch arrays.
        allocate(Xwrk(p), source=B) ; allocate(Cwrk(p), source=B(1))

        ! Normalize input matrix and initialize Krylov subspace.
        R = 0.0_dp
        call qr(Xwrk, R, perm, info) ; call apply_inverse_permutation_matrix(R, perm)

        if (norm2(abs(R)) == 0.0_dp) then
            ! Input matrix is zero.
            do i = 1, size(C)
                call C(i)%zero()
            enddo
            err_est = 0.0_dp ; k = 0 ; kpp = p
        else
            call initialize_krylov_subspace(X, Xwrk) ; H = 0.0_dp

            expm_arnoldi: do k = 1, nk
                ! Set counters.
                kpm = (k-1)*p ; kp = kpm + p ; kpp = kp + p

                ! Reset working arrays.
                E = 0.0_dp
                do i = 1, size(Xwrk)
                    call Xwrk(i)%zero()
                enddo

                ! Compute the k-th step of the Arnoldi factorization.
                call arnoldi(A, X, H, info, kstart=k, kend=k, transpose=transpose, blksize=p)
                call check_info(info, 'arnoldi', module=this_module, procedure='kexpm_mat_cdp')


                if (info == kp) then
                    ! Arnoldi breakdown. Do not consider extended matrix.
                    kpp = kp
                    info = -2
                endif

                ! Compute the (dense) matrix exponential of the extended Hessenberg matrix.
                call expm(E(:kpp, :kpp), tau*H(:kpp, :kpp))

                ! Project back to original space.
                do i = 1, size(Xwrk)
                    call Xwrk(i)%zero()
                    do j = 1, kpp
                        call Xwrk(i)%axpby(one_cdp, X(j), E(j, i))
                    enddo
                enddo

                do i = 1, p
                    call C(i)%zero()
                    do j = 1, p
                        call C(i)%axpby(one_cdp, Xwrk(j), R(j, i))
                    enddo
                enddo

                ! Cheap error estimate.
                if (info == kp) then
                    ! Approximation is exact.
                    err_est = 0.0_dp
                else
                    em = matmul(E(kp+1:kpp, :p), R(:p, :p))
                    err_est = norm2(abs(em))
                endif

                if (err_est <= tol) exit expm_arnoldi

            enddo expm_arnoldi
        endif

        if (err_est .le. tol) then
            info = kpp
            if (p.eq.1) then
               write(msg, *) new_line('A'),' Arnoldi approxmation of the action of the exp. propagator converged!'
            else
               write(msg, *) new_line('A'),' Block Arnoldi approxmation of the action of the exp. propagator converged!'
            endif 
            write(msg, '(A,2(A,I3))') trim(msg)//new_line('A'), 'n° of vectors:', k+1, ' per input vector, total:', kpp
            write(msg, '(A,4X,A,E10.4)') trim(msg)//new_line('A'), 'estimated error:   ||err_est||_2 = ', err_est
            write(msg, '(A,4X,A,E10.4)') trim(msg)//new_line('A'), 'desired tolerance:           tol = ', tol
        else
            info = -1
            if (p.eq.1) then
                write(msg, *) 'Arnoldi approxmation of the action of the exp. propagator did not converge!'
            else
                write(msg, *) 'Block Arnoldi approxmation of the action of the exp. propagator did not converge!'
            endif
            write(msg, '(A,4X,2(A,I3))') trim(msg)//new_line('A'), 'maximum n° of vectors reached: ', nsteps+1, &
                                   & ' per input vector, total:', kpp
            write(msg, '(A,4X,A,E10.4)') trim(msg)//new_line('A'), 'estimated error:   ||err_est||_2 = ', err_est
            write(msg, '(A,4X,A,E10.4)') trim(msg)//new_line('A'), 'desired tolerance:           tol = ', tol
        endif

        if (verbose) call logger%log_message(trim(msg), module=this_module, procedure='kexpm_mat_cdp')

        return
    end subroutine kexpm_mat_cdp

    subroutine kexpm_vec_variable_dt_cdp(c, A, b, tau, dt, tol, info, trans, verbosity, kdim)
        class(abstract_vector_cdp), intent(out) :: c
        !! Best approximation of \( \exp(\tau \mathbf{A}) \mathbf{b} \) in the computed Krylov subspace
        class(abstract_linop_cdp), intent(in) :: A
        !! Linear operator to be exponentiated.
        class(abstract_vector_cdp), intent(in) :: b
        !! Input vector on which to apply \( \exp(\tau \mathbf{A}) \).
        real(dp), intent(in) :: tau
        !! Time horizon for the exponentiation.
        real(dp), intent(inout) :: dt
        !! Current time step.
        real(dp), intent(in) :: tol
        !! Solution tolerance based on error estimates.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: trans
        !! Use transpose?
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        integer, optional, intent(in) :: kdim
        !! Maximum size of the Krylov subspace.

        ! ----- Internal variables -----
        integer, parameter :: kmax = 100
        integer :: k, km, kp, nk, tstep, rstep
        ! Arnoldi factorization.
        class(abstract_vector_cdp), allocatable :: X(:)
        complex(dp), allocatable :: H(:, :)
        ! Normaliztion & temporary arrays.
        complex(dp), allocatable :: E(:, :)
        class(abstract_vector_cdp), allocatable :: b_
        real(dp) :: err_est, beta, Tend, dt_, tol_, scale, sfact
        character(len=128) :: msg

        ! Optional arguments.
        logical :: transpose
        logical :: verbose
        integer :: nsteps
    
        ! Deals with optional args.
        transpose = optval(trans, .false.)
        verbose   = optval(verbosity, .false.)
        nsteps    = optval(kdim, kmax)
        nk        = nsteps

        info = 0
        scale = 1.5_dp
        sfact = 0.5_dp

        ! Allocate arrays.
        allocate(X(nk+1), b_, source=b)
        allocate(H(nk+1, nk+1), E(nk+1, nk+1))

        ! Normalize input vector and initialize Krylov subspace.
        beta = b%norm()
        if (beta == 0.0_dp) then
            ! Input is zero => Output is zero.
            call c%zero()
        else
            tstep = 1
            Tend = 0.0_dp
            if (verbose) then
               write(msg, '(3(A,E12.6))') 'tau = ', tau, ', dt_0 = ', dt, ', tol = ', tol
               call logger%log_message(trim(msg), module=this_module, procedure='kexpm_vec_variable_dt_cdp')
            end if
            do while ( Tend < tau )

               ! (Re-)initialization of X
               call initialize_krylov_subspace(X) ; call X(1)%add(b_) ;
               beta = b_%norm() ; call X(1)%scal(one_cdp / beta)
               H = 0.0_dp

               ! choose time step
               dt  = min(tau, scale*dt)    ! Reuse step size of last call (and try to increase)
               if (dt > tau-Tend) then
                  dt_ = tau-Tend       ! If a small final step is needed to reach Tend, ignore it on restart
                  if (verbose) then
                     write(msg, '(A,E12.6)') 'Filler step, set dt_ = ', dt_
                     call logger%log_message(trim(msg), module=this_module, procedure='kexpm_vec_variable_dt_cdp')
                  end if
               else
                  dt_ = dt
               end if
               tol_ = dt_/tau*tol*sfact
               if (verbose) then
                  write(msg, '(A,I3,4(A,E8.2))') 'step ', tstep, ': T= ', Tend, ', dt= ', dt, ', dt_= ', dt_, ', tol_= ', tol_
                  call logger%log_message(trim(msg), module=this_module, procedure='kexpm_vec_variable_dt_cdp')
               end if

               expm_arnoldi: do k = 1, nk
                   km = k-1 ; kp = k+1

                   ! Reset work arrays.
                   E = 0.0_dp

                   ! Compute k-th step Arnoldi factorization.
                   call arnoldi(A, X, H, info, kstart=k, kend=k, transpose=transpose)
                   call check_info(info, 'arnoldi', module=this_module, procedure='kexpm_vec_variable_dt_cdp')

                   ! Compute approximation.
                   if (info == k) then
                       ! Arnoldi breakdown, do not consider extended matrix.
                       kp = k
                       info = -2
                   endif

                   ! Compute the (dense) matrix exponential of the extended Hessenberg matrix.
                   call expm(E(:kp, :kp), dt_*H(:kp, :kp))

                   ! Cheap error esimate (this actually is the magnitude of the included correction
                   ! and thus is very conservative).
                   if (info == k) then
                       ! Approximation is exact.
                       err_est = 0.0_dp
                   else
                       err_est = abs(E(kp, 1) * beta)
                   endif

                   ! Check convergence.
                   if (err_est <= tol_) exit expm_arnoldi

               enddo expm_arnoldi

               if (verbose) then
                  write(msg, '(34X,A,E12.6)') 'err_est = ', err_est
                  call logger%log_message(trim(msg), module=this_module, procedure='kexpm_vec_variable_dt_cdp')
               end if

               if (err_est > tol_) then ! decrease dt_ until error estimate is below tolerance
                                       ! NOTE: We do not need to recompute X,H!
                  rstep = 1
                  do while (err_est > tol_)

                     !hfact = (tol/err_est)**(one_cdp/(kdim+1))
                     !dt_ = 0.5_dp*sfact_dp*hfact*dt_ 
                     dt_ = dt_/scale ! this is a quick hack, should choose step better

                     tol_ = dt_/tau*tol*sfact ! update accepted tol based on new time step

                     if (dt_ < atol_dp) then
                        write(msg, *) 'dt < atol_dp. Cannot compute action of the matrix exponential!'
                        call stop_error(trim(msg), module=this_module, procedure='kexpm_vec_variable_dt_cdp')
                     end if

                     ! Compute the (dense) matrix exponential of the extended Hessenberg matrix.
                     call expm(E(:kp, :kp), dt_*H(:kp, :kp))

                     ! Compute new error estimate
                     err_est = abs(E(kp, 1) * beta)

                     if (verbose) then
                        write(msg, '(6X,3(A,E10.4))') 'v dt_ = ', dt_, ', err_est = ', err_est,',  tol_ = ', tol_
                        call logger%log_message(trim(msg), module=this_module, procedure='kexpm_vec_variable_dt_cdp')
                     end if

                     rstep = rstep + 1

                  end do

                  ! save dt for reuse
                  dt = dt_

               end if

               ! Project back into original space.
               block
                  class(abstract_vector_cdp), allocatable :: xwrk
                  call linear_combination(xwrk, X(:kp), E(:kp, 1))
                  call b_%axpby(zero_cdp, xwrk, one_cdp*beta)
               end block

               ! move forward in time
               Tend = Tend + dt_

               tstep = tstep + 1

               if (verbose) then
                  write(msg, '(A,E10.4,2(A,I3),A,E10.4)') '    post: T = ', Tend, ', k = ', k, '/', kdim+1, ', dt => ', dt
                  call logger%log_message(trim(msg), module=this_module, procedure='kexpm_vec_variable_dt_cdp')
               end if
           
            end do ! integration

            ! copy result into output vector
            call c%axpby(zero_cdp, b_, one_cdp)

        endif      

        return
    end subroutine kexpm_vec_variable_dt_cdp

    subroutine kexpm_mat_variable_dt_cdp(C, A, B, tau, dt, tol, info, trans, verbosity, kdim)
        class(abstract_vector_cdp), intent(out) :: C(:)
        !! Best Krylov approximation of \( \mathbf{C} = \exp(\tau \mathbf{A}) \mathbf{B} \).
        class(abstract_linop_cdp), intent(in) :: A
        !! Linear operator to be exponentiated.
        class(abstract_vector_cdp), intent(in) :: B(:)
        !! Input matrix on which to apply \( \exp(\tau \mathbf{A}) \).
        real(dp), intent(in) :: tau
        !! Time horizon for the exponentiation.
        real(dp), intent(inout) :: dt
        !! Current time step.
        real(dp), intent(in) :: tol
        !! Solution toleance based on error estimates.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: trans
        !! Use transpose ?
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        integer, optional, intent(in) :: kdim
        !! Maximum size of the Krylov subspace.

        ! ----- Internal variables -----
        integer, parameter :: kmax = 100
        integer :: i, j, k, p, kpm, kp, kpp, nk, tstep
        ! Block-Arnoldi factorization.
        class(abstract_vector_cdp), allocatable :: X(:)
        complex(dp), allocatable :: H(:, :)
        ! Normalization & temporary arrays.
        complex(dp), allocatable :: R(:, :), E(:, :), em(:, :)
        integer, allocatable :: perm(:), ptrans(:)
        class(abstract_vector_cdp), allocatable :: Xwrk(:), Cwrk(:), B_(:)
        real(dp) :: err_est, Tend, dt_, tol_, scale
        character(len=128) :: msg

        ! Optional arguments.
        logical :: transpose
        logical :: verbose
        integer :: nsteps

        ! Determine block size.
        p = size(B)

        ! Deals with the optional args.
        transpose = optval(trans, .false.)
        verbose   = optval(verbosity, .false.)
        nsteps    = optval(kdim, kmax)
        nk        = nsteps*p
        
        info = 0
        scale = 1.5

        ! Allocate arrays.
        allocate(perm(p), ptrans(p))
        allocate(H(p*(nk+1), p*(nk+1)), E(p*(nk+1), p*(nk+1)), R(p, p), em(p, p))
        allocate(X(p*(nk+1)), source=B(1))

        ! Scratch arrays.
        allocate(B_(p), source=B) ; allocate(Cwrk(p), source=B(1))

        if (norm2(abs(R)) == 0.0_dp) then
            ! Input matrix is zero.
            call zero_basis(C)
        else
            tstep = 1
            Tend = 0.0_dp
            do while ( Tend < tau )

               ! Normalize B_ and (re-)initialize Krylov subspace.
               R = 0.0_dp
               call qr(B_, R, perm, info) ; call apply_inverse_permutation_matrix(R, perm)
               call initialize_krylov_subspace(X, B_)
               H = 0.0_dp

               ! choose time step
               dt  = min(tau, scale*dt)    ! Reuse step size of last call (and try to increase)
               if (dt > tau-Tend) then
                  dt_ = tau-Tend       ! If a small final step is needed to reach Tend, ignore it on restart
                  if (verbose) then
                     write(msg, '(A,E12.6)') 'Filler step, set dt_ = ', dt_
                     call logger%log_message(trim(msg), module=this_module, procedure='kexpm_mat_variable_dt_cdp')
                  end if
               else
                  dt_ = dt
               end if
               tol_ = dt_/tau*tol
               if (verbose) then
                  write(msg, '(A,I3,4(A,E8.2))') 'step ', tstep, ': T= ', Tend, ', dt= ', dt, ', dt_= ', dt_, ', tol_= ', tol_
                  call logger%log_message(trim(msg), module=this_module, procedure='kexpm_mat_variable_dt_cdp')
               end if

               expm_arnoldi: do k = 1, nk
                   ! Set counters.
                   kpm = (k-1)*p ; kp = kpm + p ; kpp = kp + p

                   ! Reset working arrays.
                   E = 0.0_dp

                   ! Compute the k-th step of the Arnoldi factorization.
                   call arnoldi(A, X, H, info, kstart=k, kend=k, transpose=transpose, blksize=p)
                   call check_info(info, 'arnoldi', module=this_module, procedure='kexpm_mat_variable_dt_cdp')

                   if (info == kp) then
                       ! Arnoldi breakdown. Do not consider extended matrix.
                       kpp = kp
                       info = -2
                   endif

                   ! Compute the (dense) matrix exponential of the extended Hessenberg matrix.
                   call expm(E(:kpp, :kpp), dt_*H(:kpp, :kpp))

                   ! Cheap error estimate.
                   if (info == kp) then
                       ! Approximation is exact.
                       err_est = 0.0_dp
                   else
                       em = matmul(E(kp+1:kpp, :p), R(:p, :p))
                       err_est = norm2(abs(em))
                   endif

                   if (err_est <= tol_) exit expm_arnoldi

               enddo expm_arnoldi

               if (verbose) then
                  write(msg, '(34X,A,E10.4)') 'err_est = ', err_est
                  call logger%log_message(trim(msg), module=this_module, procedure='kexpm_mat_variable_dt_cdp')
               end if

               if (err_est > tol_) then ! decrease dt_ until error estimate is below tolerance
                                       ! NOTE: We do not need to recompute X,H!
                  do while (err_est > tol_)

                     !hfact = (tol/err_est)**(one_cdp/kdim)
                     !dt_ = sfact_dp*hfact*dt_
                     dt_ = dt_/scale ! this is a quick hack, should choose step better

                     tol_ = dt_/tau*tol ! update accepted tol based on new time step

                     if (dt_ < atol_dp) then
                        write(msg, *) 'dt < atol_dp. Cannot compute action of the matrix exponential!'
                        call stop_error(trim(msg), module=this_module, procedure='kexpm_mat_variable_dt_cdp')
                     end if

                     ! Compute the (dense) matrix exponential of the extended Hessenberg matrix.
                     call expm(E(:kpp, :kpp), dt_*H(:kpp, :kpp))

                     ! Compute new error estimate
                     em = matmul(E(kp+1:kpp, :p), R(:p, :p))
                     err_est = norm2(abs(em))

                     if (verbose) then
                        write(msg, '(6X,3(A,E12.6))') 'v dt_ = ', dt_, ', err_est = ', err_est, ', tol_ = ', tol_
                        call logger%log_message(trim(msg), module=this_module, procedure='kexpm_mat_variable_dt_cdp')
                     end if

                  end do

                  ! save dt for reuse
                  dt = dt_
               
               end if

               ! Project back into original space.
               block
                  class(abstract_vector_cdp), allocatable :: Xwrk(:)
                  call linear_combination(Xwrk, X(:kp), matmul(E(:kp, :p), R))
               end block

               ! move forward in time
               Tend = Tend + dt_

               tstep = tstep + 1

               if (verbose) then
                  write(msg, '(A,E10.4,2(A,I3),A,E10.4)') '    post: T = ', Tend, ', k = ', k, '/', kdim+1, ', dt => ', dt
                  call logger%log_message(trim(msg), module=this_module, procedure='kexpm_mat_variable_dt_cdp')
               end if
           
            end do ! integration

            ! copy result into output vector
            call copy_basis(C, B_)
        endif

        return
    end subroutine kexpm_mat_variable_dt_cdp

    subroutine k_exptA_cdp(vec_out, A, vec_in, tau, info, trans)
        class(abstract_vector_cdp), intent(out) :: vec_out
        !! Solution vector.
        class(abstract_linop_cdp), intent(inout) :: A
        !! Linear operator to be exponentiated.
        class(abstract_vector_cdp), intent(in) :: vec_in
        !! Input vector to be multiplied by \( \exp(\tau \mathbf{A}) \).
        real(dp), intent(in) :: tau
        !! Time horizon for the exponentiation.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: trans
        !! Use adjoint ?

        ! ----- Internal variables -----
        real(dp) :: tol
        integer :: kdim
        logical :: verbose

        tol = rtol_dp
        kdim = 30
        verbose = .false.

        call kexpm(vec_out, A, vec_in, tau, tol, info, trans=trans, verbosity=verbose, kdim=kdim)
        call check_info(info, 'kexpm', module=this_module, procedure='k_exptA_cdp')

        return
    end subroutine k_exptA_cdp


end module lightkrylov_expmlib
