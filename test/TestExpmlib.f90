module TestExpmlib
    ! Fortran Standard Library.
    use stdlib_linalg, only: eye, diag, expm, hermitian
    ! Testdrive
    use testdrive, only: new_unittest, unittest_type, error_type, check
    ! LightKrylov
    use LightKrylov
    use LightKrylov_Constants
    use LightKrylov_Logger
    use LightKrylov_Utils, only : eig, sqrtm
    ! Test Utilities
    use LightKrylov_TestUtils
    use TestUtils

    implicit none (type, external)

    private

    character(len=*), parameter, private :: this_module      = 'LK_TExpmLib'
    character(len=*), parameter, private :: this_module_long = 'LightKrylov_TestExpmLib'

    public :: collect_expm_rsp_testsuite
    public :: collect_expm_rdp_testsuite
    public :: collect_expm_csp_testsuite
    public :: collect_expm_cdp_testsuite

    public :: collect_sqrtm_rsp_testsuite
    public :: collect_sqrtm_rdp_testsuite
    public :: collect_sqrtm_csp_testsuite
    public :: collect_sqrtm_cdp_testsuite

contains

    !--------------------------------------------------------------
    !-----     UNIT TESTS FOR DENSE MATRIX EXPONENTIATION     -----
    !--------------------------------------------------------------

    subroutine collect_expm_rsp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                        new_unittest("Krylov expm invalid inputs", test_kexptA_invalid_inputs_rsp), &
                        new_unittest("Dense expm", test_dense_expm_rsp), &
                        new_unittest("Krylov expm", test_kexptA_rsp), &
                        new_unittest("Krylov expm adjoint", test_kexptA_adjoint_rsp), &
                        new_unittest("Krylov expm zero input", test_kexptA_zero_input_rsp), &
                        new_unittest("Krylov expm invariant subspace", test_kexptA_breakdown_rsp), &
                        new_unittest("Krylov expm stiff (sub-stepping)", test_kexptA_stiff_rsp), &
                        new_unittest("Krylov expm warm start", test_kexptA_warm_start_rsp), &
                        new_unittest("Block Krylov expm", test_block_kexptA_rsp) &
                    ]

        return
    end subroutine collect_expm_rsp_testsuite

    subroutine test_kexptA_invalid_inputs_rsp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test linear operator.
        type(linop_rsp), allocatable :: A
        ! Test vectors.
        type(vector_rsp), allocatable :: b, c, Bblk(:), Cblk(:), X(:)
        ! Sub-stepping options and metadata.
        type(kexpm_sp_opts) :: opts
        type(kexpm_sp_metadata) :: meta
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        integer, parameter :: p = 2
        real(sp), parameter :: tau = 0.1_sp
        integer :: icase
        character(len=256) :: msg

        ! Common initialization.
        A = linop_rsp() ; call init_rand(A)
        allocate(b, c) ; call init_rand(b) ; call c%zero()
        allocate(Bblk(p)) ; call init_rand(Bblk)

        test_loop: do
            !-----     Single vector kexpm     -----

            ! --- Case 1: workspace with fewer than 2 vectors (info = -3) ---
            icase = 1
            allocate(X(1)) ; call zero_basis(X)
            call kexpm(c, A, b, tau, rtol_sp, info, X=X)
            call check(error, info == -3)
            if (allocated(error)) exit test_loop

            ! --- Case 2: kdim < 1 (info = -4) ---
            icase = 2
            call kexpm(c, A, b, tau, rtol_sp, info, kdim=0)
            call check(error, info == -4)
            if (allocated(error)) exit test_loop

            ! --- Case 3: negative tolerance (info = -5) ---
            icase = 3
            call kexpm(c, A, b, tau, -1.0_sp, info)
            call check(error, info == -5)
            if (allocated(error)) exit test_loop

            !-----     Block kexpm     -----

            ! --- Case 4: size(C) /= size(B) (info = -6) ---
            icase = 4
            allocate(Cblk(p+1)) ; call zero_basis(Cblk)
            call kexpm(Cblk, A, Bblk, tau, rtol_sp, info)
            call check(error, info == -6)
            if (allocated(error)) exit test_loop

            ! --- Case 5: empty block (info = -6) ---
            icase = 5
            block
                type(vector_rsp), allocatable :: B0(:), C0(:)
                allocate(B0(0), C0(0))
                call kexpm(C0, A, B0, tau, rtol_sp, info)
            end block
            call check(error, info == -6)
            if (allocated(error)) exit test_loop

            ! --- Case 6: workspace size not a multiple of the block size (info = -3) ---
            icase = 6
            deallocate(Cblk, X) ; allocate(Cblk(p), X(2*p+1)) ; call zero_basis(Cblk) ; call zero_basis(X)
            call kexpm(Cblk, A, Bblk, tau, rtol_sp, info, X=X)
            call check(error, info == -3)
            if (allocated(error)) exit test_loop

            ! --- Case 7: workspace with fewer than 2*p vectors (info = -3) ---
            icase = 7
            deallocate(X) ; allocate(X(p)) ; call zero_basis(X)
            call kexpm(Cblk, A, Bblk, tau, rtol_sp, info, X=X)
            call check(error, info == -3)
            if (allocated(error)) exit test_loop

            ! --- Case 8: kdim < 1 (info = -4) ---
            icase = 8
            call kexpm(Cblk, A, Bblk, tau, rtol_sp, info, kdim=0)
            call check(error, info == -4)
            if (allocated(error)) exit test_loop

            ! --- Case 9: negative tolerance (info = -5) ---
            icase = 9
            call kexpm(Cblk, A, Bblk, tau, -1.0_sp, info)
            call check(error, info == -5)
            if (allocated(error)) exit test_loop

            !-----     Sub-stepped kexpm_substep     -----

            ! --- Case 10: opts%kdim < 1 (info = -4) ---
            icase = 10
            opts = kexpm_sp_opts(kdim=0)
            call kexpm_substep(c, A, b, tau, rtol_sp, info, options=opts)
            call check(error, info == -4)
            if (allocated(error)) exit test_loop

            ! --- Case 11: negative tolerance (info = -5) ---
            icase = 11
            call kexpm_substep(c, A, b, tau, -1.0_sp, info)
            call check(error, info == -5)
            if (allocated(error)) exit test_loop

            ! --- Case 12: opts%max_rejections < 1 (info = -7) ---
            icase = 12
            opts = kexpm_sp_opts(max_rejections=0)
            call kexpm_substep(c, A, b, tau, rtol_sp, info, options=opts, meta=meta)
            call check(error, info == -7)
            if (allocated(error)) exit test_loop

            ! --- Case 13: opts%max_substeps < 1 (info = -7), metadata left untouched ---
            icase = 13
            opts = kexpm_sp_opts(max_substeps=0)
            call kexpm_substep(c, A, b, tau, rtol_sp, info, options=opts, meta=meta)
            call check(error, info == -7 .and. meta%n_calls == 0 .and. .not. meta%learned)
            exit test_loop
        end do test_loop
        write(msg, '(A,I0)') 'last case: ', icase
        call check_test(error, 'test_kexptA_invalid_inputs_rsp', &
                              & info='Invalid input parameters', eq='', context=msg)

        return
    end subroutine test_kexptA_invalid_inputs_rsp

    subroutine test_dense_expm_rsp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Problem dimension.
        integer, parameter :: n = 5, m = 6
        ! Test matrix.
        real(sp) :: A(n, n), E(n, n), Eref(n, n)
        ! Miscellaneous.
        integer :: i, j
        real(sp) :: err
        character(len=256) :: msg

        ! Initialize nilpotent matrix (the exponential is a finite series).
        A = 0.0_sp
        do i = 1, n-1
            A(i, i+1) = m*1.0_sp
        enddo

        ! Reference with analytical exponential: Eref(i, i+j) = m**j / j!.
        Eref = eye(n, mold=1.0_sp)
        do i = 1, n-1
            do j = 1, n-i
                Eref(i, i+j) = Eref(i, i+j-1)*m/j
            enddo
        enddo

        ! Compute matrix exponential.
        E = expm(A)

        ! Check result.
        err = maxval(abs(E-Eref))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_dense_expm_rsp', &
                              & info='Nilpotent matrix', eq='expm(A) = Eref', context=msg)

        return
    end subroutine test_dense_expm_rsp

    subroutine test_kexptA_rsp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test linear operator.
        type(linop_rsp), allocatable :: A
        ! Test vectors.
        type(vector_rsp), allocatable :: b, c, cref
        ! Krylov subspace (user-supplied workspace).
        integer, parameter :: kdim = 64
        type(vector_rsp), allocatable :: X(:)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(sp), parameter :: tau = 0.1_sp
        real(sp) :: err
        character(len=256) :: msg

        ! Initialize linear operator and vectors.
        A = linop_rsp() ; call init_rand(A)
        allocate(b, c, cref) ; call init_rand(b) ; call c%zero()
        allocate(X(kdim+1)) ; call zero_basis(X)

        ! Reference: dense computation.
        cref%data = matmul(expm(tau*A%data), b%data)

        ! kexpm with internal workspace.
        call kexpm(c, A, b, tau, rtol_sp, info, kdim=kdim)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_rsp')
        call c%sub(cref) ; err = c%norm() / cref%norm()
        call get_err_str(msg, "max err: ", err)
        call check(error, info > 0 .and. err < rtol_sp)
        call check_test(error, 'test_kexptA_rsp', &
                              & info='kexpm', eq='c = expm(tau A) @ b', context=msg)

        ! krylov_exptA (adaptive sub-stepping).
        call krylov_exptA(c, A, b, tau, info)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_rsp')
        call c%sub(cref) ; err = c%norm() / cref%norm()
        call get_err_str(msg, "max err: ", err)
        call check(error, info > 0 .and. err < rtol_sp)
        call check_test(error, 'test_kexptA_rsp', &
                              & info='krylov_exptA', eq='c = expm(tau A) @ b', context=msg)

        ! kexpm with user-supplied workspace.
        call kexpm(c, A, b, tau, rtol_sp, info, X=X)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_rsp')
        call c%sub(cref) ; err = c%norm() / cref%norm()
        call get_err_str(msg, "max err: ", err)
        call check(error, info > 0 .and. err < rtol_sp)
        call check_test(error, 'test_kexptA_rsp', &
                              & info='kexpm, workspace', eq='c = expm(tau A) @ b', context=msg)

        ! Re-used (non-zeroed) workspace with kdim smaller than the workspace.
        call kexpm(c, A, b, tau, rtol_sp, info, kdim=kdim/2, X=X)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_rsp')
        call c%sub(cref) ; err = c%norm() / cref%norm()
        call get_err_str(msg, "max err: ", err)
        call check(error, info > 0 .and. info <= kdim/2 + 1 .and. err < rtol_sp)
        call check_test(error, 'test_kexptA_rsp', &
                              & info='kexpm, workspace + kdim', eq='c = expm(tau A) @ b', context=msg)

        return
    end subroutine test_kexptA_rsp

    subroutine test_kexptA_adjoint_rsp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test linear operator.
        type(linop_rsp), allocatable :: A
        ! Test vectors.
        type(vector_rsp), allocatable :: b, c, cref
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(sp), parameter :: tau = 0.1_sp
        real(sp) :: err
        character(len=256) :: msg

        ! Initialize linear operator and vectors.
        A = linop_rsp() ; call init_rand(A)
        allocate(b, c, cref) ; call init_rand(b) ; call c%zero()

        ! Reference: dense computation with the adjoint operator.
        cref%data = matmul(expm(tau*hermitian(A%data)), b%data)

        ! kexpm.
        call kexpm(c, A, b, tau, rtol_sp, info, trans=.true.)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_adjoint_rsp')
        call c%sub(cref) ; err = c%norm() / cref%norm()
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_kexptA_adjoint_rsp', &
                              & info='kexpm', eq='c = expm(tau A.H) @ b', context=msg)

        ! krylov_exptA.
        call krylov_exptA(c, A, b, tau, info, trans=.true.)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_adjoint_rsp')
        call c%sub(cref) ; err = c%norm() / cref%norm()
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_kexptA_adjoint_rsp', &
                              & info='krylov_exptA', eq='c = expm(tau A.H) @ b', context=msg)

        return
    end subroutine test_kexptA_adjoint_rsp

    subroutine test_kexptA_zero_input_rsp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test linear operator.
        type(linop_rsp), allocatable :: A
        ! Test vectors.
        integer, parameter :: p = 2
        type(vector_rsp), allocatable :: b, c, Bblk(:), Cblk(:)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(sp), parameter :: tau = 0.1_sp
        real(sp) :: err
        integer :: i
        character(len=256) :: msg

        ! Initialize linear operator, zero inputs and non-zero outputs.
        A = linop_rsp() ; call init_rand(A)
        allocate(b, c) ; call b%zero()
        allocate(Bblk(p), Cblk(p)) ; call zero_basis(Bblk)

        do
            ! kexpm.
            call c%rand()
            call kexpm(c, A, b, tau, rtol_sp, info)
            err = c%norm()
            call check(error, info > 0 .and. err == 0.0_sp)
            if (allocated(error)) exit

            ! krylov_exptA.
            call c%rand()
            call krylov_exptA(c, A, b, tau, info)
            err = c%norm()
            call check(error, info > 0 .and. err == 0.0_sp)
            if (allocated(error)) exit

            ! Block kexpm.
            call rand_basis(Cblk)
            call kexpm(Cblk, A, Bblk, tau, rtol_sp, info)
            err = 0.0_sp
            do i = 1, p
                err = err + Cblk(i)%norm()
            enddo
            call check(error, info > 0 .and. err == 0.0_sp)
            exit
        enddo
        call get_err_str(msg, "norm: ", err)
        call check_test(error, 'test_kexptA_zero_input_rsp', &
                              & info='Zero input', eq='expm(tau A) @ 0 = 0', context=msg)

        return
    end subroutine test_kexptA_zero_input_rsp

    subroutine test_kexptA_breakdown_rsp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test linear operator.
        type(linop_rsp), allocatable :: A
        ! Test vectors.
        type(vector_rsp), allocatable :: b, c, cref
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        integer, parameter :: k_inv = 10
        real(sp), parameter :: tau = 0.1_sp
        real(sp) :: err
        character(len=256) :: msg

        ! Create block lower triangular matrix: span(e_1, ..., e_k_inv) is invariant.
        A = linop_rsp() ; call init_rand(A)
        A%data(k_inv+1:, :k_inv) = zero_rsp

        ! Input vector confined to the invariant subspace.
        allocate(b, c, cref) ; call init_rand(b) ; call c%zero()
        b%data(k_inv+1:) = zero_rsp

        ! Reference: dense computation.
        cref%data = matmul(expm(tau*A%data), b%data)

        ! kexpm with tol = 0: convergence is only possible through the Arnoldi breakdown,
        ! which must happen after exactly k_inv steps (the approximation is then exact).
        call kexpm(c, A, b, tau, 0.0_sp, info)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_breakdown_rsp')
        call check(error, info == k_inv)
        call check_test(error, 'test_kexptA_breakdown_rsp', &
                              & info='Subspace Dim', eq='info == k_inv', context='Invariant detection')

        call c%sub(cref) ; err = c%norm() / cref%norm()
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_kexptA_breakdown_rsp', &
                              & info='Invariant Property', eq='c = expm(tau A) @ b', context=msg)

        ! krylov_exptA: a single sub-step, never more than k_inv matrix-vector products.
        call krylov_exptA(c, A, b, tau, info)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_breakdown_rsp')
        call c%sub(cref) ; err = c%norm() / cref%norm()
        call get_err_str(msg, "max err: ", err)
        call check(error, info <= k_inv .and. err < rtol_sp)
        call check_test(error, 'test_kexptA_breakdown_rsp', &
                              & info='krylov_exptA', eq='<= k_inv matvecs, c = expm(tau A) @ b', context=msg)

        return
    end subroutine test_kexptA_breakdown_rsp

    subroutine test_kexptA_stiff_rsp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test linear operator.
        type(linop_rsp), allocatable :: A
        ! Test vectors.
        type(vector_rsp), allocatable :: b, c, cref
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(sp), parameter :: tau = 1.0_sp, s = 30.0_sp
        real(sp) :: err
        character(len=256) :: msg

        ! Stable random operator with || tau A || ~ 60: far beyond what a single
        ! 30-vector Krylov basis can resolve.
        A = linop_rsp() ; call init_rand(A)
        A%data = s*A%data/sqrt(real(test_size, sp)) - s*eye(test_size, mold=1.0_sp)
        allocate(b, c, cref) ; call init_rand(b) ; call c%zero()

        ! Reference: dense computation.
        cref%data = matmul(expm(tau*A%data), b%data)

        ! krylov_exptA must converge thanks to the adaptive sub-stepping (> 30 matvecs).
        call krylov_exptA(c, A, b, tau, info)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_stiff_rsp')
        call c%sub(cref) ; err = c%norm() / cref%norm()
        call get_err_str(msg, "max err: ", err)
        call check(error, info > 30 .and. err < rtol_sp)
        call check_test(error, 'test_kexptA_stiff_rsp', &
                              & info='Sub-stepping', eq='c = expm(tau A) @ b', context=msg)

        return
    end subroutine test_kexptA_stiff_rsp

    subroutine test_kexptA_warm_start_rsp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test linear operator.
        type(linop_rsp), allocatable :: A
        ! Test vectors.
        type(vector_rsp), allocatable :: b, c, cref
        ! Sub-stepping options and metadata.
        type(kexpm_sp_opts) :: opts
        type(kexpm_sp_metadata) :: meta
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        integer, parameter :: ncalls = 10
        real(sp), parameter :: tau = 1.0_sp, s = 30.0_sp
        real(sp), allocatable :: E(:, :)
        real(sp) :: err, tol
        integer :: icall, cold, last
        character(len=256) :: msg

        ! Stiff stable random operator.
        A = linop_rsp() ; call init_rand(A)
        A%data = s*A%data/sqrt(real(test_size, sp)) - s*eye(test_size, mold=1.0_sp)
        E = expm(tau*A%data)
        allocate(b, c, cref) ; call c%zero()

        ! Repeated calls with a persistent metadata object and new right-hand sides.
        err = 0.0_sp
        do icall = 1, ncalls
            call init_rand(b) ; cref%data = matmul(E, b%data) ; tol = atol_sp*b%norm()
            call kexpm_substep(c, A, b, tau, tol, info, meta=meta)
            call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_warm_start_rsp')
            if (icall == 1) cold = info
            last = info
            call c%sub(cref) ; err = max(err, c%norm() / cref%norm())
        enddo
        write(msg, '(A,E9.2,A,I0,A,I0,A,I0)') 'max err: ', err, ', cold: ', cold, ', last: ', last, ', nsub: ', meta%nsub
        call check(error, err < rtol_sp)
        call check_test(error, 'test_kexptA_warm_start_rsp', &
                              & info='Accuracy', eq='c = expm(tau A) @ b', context=msg)

        ! Settings learned, frozen, and not more expensive than the cold call (up to the
        ! call-to-call variation of the adaptive cost with the right-hand side).
        call check(error, meta%learned .and. meta%frozen .and. meta%n_calls == ncalls &
                        & .and. last <= cold + max(2, cold/10))
        call check_test(error, 'test_kexptA_warm_start_rsp', &
                              & info='Learned settings', eq='frozen, cost <= cold cost', context=msg)

        ! reset() keeps the learned settings, forget() clears them.
        do
            call meta%reset()
            call check(error, meta%learned .and. meta%n_calls == 0 .and. meta%n_matvec == 0)
            if (allocated(error)) exit
            call meta%forget()
            call check(error, .not. meta%learned .and. meta%nsub == 0)
            exit
        enddo
        call check_test(error, 'test_kexptA_warm_start_rsp', &
                              & info='reset / forget', eq='reset keeps, forget clears', context='')

        ! tune = .false. : the sub-step count of the cold call is kept.
        opts = kexpm_sp_opts(tune=.false.)
        call init_rand(b) ; cref%data = matmul(E, b%data) ; tol = atol_sp*b%norm()
        call kexpm_substep(c, A, b, tau, tol, info, options=opts, meta=meta)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_warm_start_rsp')
        call kexpm_substep(c, A, b, tau, tol, info, options=opts, meta=meta)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_warm_start_rsp')
        call c%sub(cref) ; err = c%norm() / cref%norm()
        write(msg, '(A,E9.2,A,I0,A,I0)') 'max err: ', err, ', nsub: ', meta%nsub, ', nsub0: ', meta%nsub0
        call check(error, err < rtol_sp .and. meta%frozen .and. meta%nsub == meta%nsub0)
        call check_test(error, 'test_kexptA_warm_start_rsp', &
                              & info='No tuning', eq='nsub == nsub0', context=msg)

        ! A different horizon invalidates the learned settings.
        cref%data = matmul(expm(0.5_sp*tau*A%data), b%data)
        call kexpm_substep(c, A, b, 0.5_sp*tau, tol, info, meta=meta)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_warm_start_rsp')
        call c%sub(cref) ; err = c%norm() / cref%norm()
        write(msg, '(A,E9.2,A,E9.2)') 'max err: ', err, ', tau_learned: ', meta%tau_learned
        call check(error, err < rtol_sp .and. abs(meta%tau_learned - 0.5_sp*tau) <= epsilon(tau))
        call check_test(error, 'test_kexptA_warm_start_rsp', &
                              & info='New horizon', eq='settings re-learned', context=msg)

        return
    end subroutine test_kexptA_warm_start_rsp

    subroutine test_block_kexptA_rsp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test linear operator.
        type(linop_rsp), allocatable :: A
        ! Test vectors.
        integer, parameter :: p = 3
        type(vector_rsp), allocatable :: B(:), C(:)
        ! Krylov subspace (user-supplied workspace).
        integer, parameter :: kdim = 15
        type(vector_rsp), allocatable :: X(:)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(sp), parameter :: tau = 0.1_sp
        real(sp), allocatable :: Bdata(:, :), Cdata(:, :), Cref(:, :)
        real(sp) :: err
        integer :: i
        character(len=256) :: msg

        ! Initialize linear operator and input block.
        A = linop_rsp() ; call init_rand(A)
        allocate(B(p), C(p), X(p*(kdim+1))) ; call init_rand(B) ; call zero_basis(C) ; call zero_basis(X)
        allocate(Bdata(test_size, p), Cdata(test_size, p)) ; call get_data(Bdata, B)

        ! Reference: dense computation.
        Cref = matmul(expm(tau*A%data), Bdata)

        ! Column-by-column kexpm.
        do i = 1, p
            call kexpm(C(i), A, B(i), tau, rtol_sp, info, kdim=kdim)
            call check_info(info, 'kexpm', module=this_module_long, procedure='test_block_kexptA_rsp')
        enddo
        call get_data(Cdata, C)
        err = norm2(abs(Cdata - Cref)) / norm2(abs(Cref))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_block_kexptA_rsp', &
                              & info='Column-wise', eq='C = expm(tau A) @ B', context=msg)

        ! Block kexpm.
        call kexpm(C, A, B, tau, rtol_sp, info, kdim=kdim)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_block_kexptA_rsp')
        call get_data(Cdata, C)
        err = norm2(abs(Cdata - Cref)) / norm2(abs(Cref))
        call get_err_str(msg, "max err: ", err)
        call check(error, info > 0 .and. info <= p*(kdim+1) .and. err < rtol_sp)
        call check_test(error, 'test_block_kexptA_rsp', &
                              & info='Block', eq='C = expm(tau A) @ B', context=msg)

        ! Block kexpm with user-supplied workspace.
        call kexpm(C, A, B, tau, rtol_sp, info, X=X)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_block_kexptA_rsp')
        call get_data(Cdata, C)
        err = norm2(abs(Cdata - Cref)) / norm2(abs(Cref))
        call get_err_str(msg, "max err: ", err)
        call check(error, info > 0 .and. err < rtol_sp)
        call check_test(error, 'test_block_kexptA_rsp', &
                              & info='Block, workspace', eq='C = expm(tau A) @ B', context=msg)

        block
        integer, parameter :: pinv = 2, k_inv = 2*pinv
        type(vector_rsp), allocatable :: B2(:), C2(:)
        real(sp), allocatable :: B2data(:, :), C2data(:, :), C2ref(:, :)

        ! Create block lower triangular matrix: span(e_1, ..., e_k_inv) is invariant.
        call init_rand(A)
        A%data(k_inv+1:, :k_inv) = zero_rsp

        ! Input block confined to the invariant subspace.
        allocate(B2(pinv), C2(pinv)) ; call init_rand(B2) ; call zero_basis(C2)
        do i = 1, pinv
            B2(i)%data(k_inv+1:) = zero_rsp
        enddo
        allocate(B2data(test_size, pinv), C2data(test_size, pinv)) ; call get_data(B2data, B2)
        C2ref = matmul(expm(tau*A%data), B2data)

        ! tol = 0: block Arnoldi must break down after k_inv/pinv block steps (exact approximation).
        call kexpm(C2, A, B2, tau, 0.0_sp, info)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_block_kexptA_rsp')
        call check(error, info == k_inv)
        call check_test(error, 'test_block_kexptA_rsp', &
                              & info='Subspace Dim', eq='info == k_inv', context='Invariant detection')

        call get_data(C2data, C2)
        err = norm2(abs(C2data - C2ref)) / norm2(abs(C2ref))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_block_kexptA_rsp', &
                              & info='Invariant Property', eq='C = expm(tau A) @ B', context=msg)
        end block

        return
    end subroutine test_block_kexptA_rsp

    subroutine collect_expm_rdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                        new_unittest("Krylov expm invalid inputs", test_kexptA_invalid_inputs_rdp), &
                        new_unittest("Dense expm", test_dense_expm_rdp), &
                        new_unittest("Krylov expm", test_kexptA_rdp), &
                        new_unittest("Krylov expm adjoint", test_kexptA_adjoint_rdp), &
                        new_unittest("Krylov expm zero input", test_kexptA_zero_input_rdp), &
                        new_unittest("Krylov expm invariant subspace", test_kexptA_breakdown_rdp), &
                        new_unittest("Krylov expm stiff (sub-stepping)", test_kexptA_stiff_rdp), &
                        new_unittest("Krylov expm warm start", test_kexptA_warm_start_rdp), &
                        new_unittest("Block Krylov expm", test_block_kexptA_rdp) &
                    ]

        return
    end subroutine collect_expm_rdp_testsuite

    subroutine test_kexptA_invalid_inputs_rdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test linear operator.
        type(linop_rdp), allocatable :: A
        ! Test vectors.
        type(vector_rdp), allocatable :: b, c, Bblk(:), Cblk(:), X(:)
        ! Sub-stepping options and metadata.
        type(kexpm_dp_opts) :: opts
        type(kexpm_dp_metadata) :: meta
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        integer, parameter :: p = 2
        real(dp), parameter :: tau = 0.1_dp
        integer :: icase
        character(len=256) :: msg

        ! Common initialization.
        A = linop_rdp() ; call init_rand(A)
        allocate(b, c) ; call init_rand(b) ; call c%zero()
        allocate(Bblk(p)) ; call init_rand(Bblk)

        test_loop: do
            !-----     Single vector kexpm     -----

            ! --- Case 1: workspace with fewer than 2 vectors (info = -3) ---
            icase = 1
            allocate(X(1)) ; call zero_basis(X)
            call kexpm(c, A, b, tau, rtol_dp, info, X=X)
            call check(error, info == -3)
            if (allocated(error)) exit test_loop

            ! --- Case 2: kdim < 1 (info = -4) ---
            icase = 2
            call kexpm(c, A, b, tau, rtol_dp, info, kdim=0)
            call check(error, info == -4)
            if (allocated(error)) exit test_loop

            ! --- Case 3: negative tolerance (info = -5) ---
            icase = 3
            call kexpm(c, A, b, tau, -1.0_dp, info)
            call check(error, info == -5)
            if (allocated(error)) exit test_loop

            !-----     Block kexpm     -----

            ! --- Case 4: size(C) /= size(B) (info = -6) ---
            icase = 4
            allocate(Cblk(p+1)) ; call zero_basis(Cblk)
            call kexpm(Cblk, A, Bblk, tau, rtol_dp, info)
            call check(error, info == -6)
            if (allocated(error)) exit test_loop

            ! --- Case 5: empty block (info = -6) ---
            icase = 5
            block
                type(vector_rdp), allocatable :: B0(:), C0(:)
                allocate(B0(0), C0(0))
                call kexpm(C0, A, B0, tau, rtol_dp, info)
            end block
            call check(error, info == -6)
            if (allocated(error)) exit test_loop

            ! --- Case 6: workspace size not a multiple of the block size (info = -3) ---
            icase = 6
            deallocate(Cblk, X) ; allocate(Cblk(p), X(2*p+1)) ; call zero_basis(Cblk) ; call zero_basis(X)
            call kexpm(Cblk, A, Bblk, tau, rtol_dp, info, X=X)
            call check(error, info == -3)
            if (allocated(error)) exit test_loop

            ! --- Case 7: workspace with fewer than 2*p vectors (info = -3) ---
            icase = 7
            deallocate(X) ; allocate(X(p)) ; call zero_basis(X)
            call kexpm(Cblk, A, Bblk, tau, rtol_dp, info, X=X)
            call check(error, info == -3)
            if (allocated(error)) exit test_loop

            ! --- Case 8: kdim < 1 (info = -4) ---
            icase = 8
            call kexpm(Cblk, A, Bblk, tau, rtol_dp, info, kdim=0)
            call check(error, info == -4)
            if (allocated(error)) exit test_loop

            ! --- Case 9: negative tolerance (info = -5) ---
            icase = 9
            call kexpm(Cblk, A, Bblk, tau, -1.0_dp, info)
            call check(error, info == -5)
            if (allocated(error)) exit test_loop

            !-----     Sub-stepped kexpm_substep     -----

            ! --- Case 10: opts%kdim < 1 (info = -4) ---
            icase = 10
            opts = kexpm_dp_opts(kdim=0)
            call kexpm_substep(c, A, b, tau, rtol_dp, info, options=opts)
            call check(error, info == -4)
            if (allocated(error)) exit test_loop

            ! --- Case 11: negative tolerance (info = -5) ---
            icase = 11
            call kexpm_substep(c, A, b, tau, -1.0_dp, info)
            call check(error, info == -5)
            if (allocated(error)) exit test_loop

            ! --- Case 12: opts%max_rejections < 1 (info = -7) ---
            icase = 12
            opts = kexpm_dp_opts(max_rejections=0)
            call kexpm_substep(c, A, b, tau, rtol_dp, info, options=opts, meta=meta)
            call check(error, info == -7)
            if (allocated(error)) exit test_loop

            ! --- Case 13: opts%max_substeps < 1 (info = -7), metadata left untouched ---
            icase = 13
            opts = kexpm_dp_opts(max_substeps=0)
            call kexpm_substep(c, A, b, tau, rtol_dp, info, options=opts, meta=meta)
            call check(error, info == -7 .and. meta%n_calls == 0 .and. .not. meta%learned)
            exit test_loop
        end do test_loop
        write(msg, '(A,I0)') 'last case: ', icase
        call check_test(error, 'test_kexptA_invalid_inputs_rdp', &
                              & info='Invalid input parameters', eq='', context=msg)

        return
    end subroutine test_kexptA_invalid_inputs_rdp

    subroutine test_dense_expm_rdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Problem dimension.
        integer, parameter :: n = 5, m = 6
        ! Test matrix.
        real(dp) :: A(n, n), E(n, n), Eref(n, n)
        ! Miscellaneous.
        integer :: i, j
        real(dp) :: err
        character(len=256) :: msg

        ! Initialize nilpotent matrix (the exponential is a finite series).
        A = 0.0_dp
        do i = 1, n-1
            A(i, i+1) = m*1.0_dp
        enddo

        ! Reference with analytical exponential: Eref(i, i+j) = m**j / j!.
        Eref = eye(n, mold=1.0_dp)
        do i = 1, n-1
            do j = 1, n-i
                Eref(i, i+j) = Eref(i, i+j-1)*m/j
            enddo
        enddo

        ! Compute matrix exponential.
        E = expm(A)

        ! Check result.
        err = maxval(abs(E-Eref))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_dense_expm_rdp', &
                              & info='Nilpotent matrix', eq='expm(A) = Eref', context=msg)

        return
    end subroutine test_dense_expm_rdp

    subroutine test_kexptA_rdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test linear operator.
        type(linop_rdp), allocatable :: A
        ! Test vectors.
        type(vector_rdp), allocatable :: b, c, cref
        ! Krylov subspace (user-supplied workspace).
        integer, parameter :: kdim = 64
        type(vector_rdp), allocatable :: X(:)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(dp), parameter :: tau = 0.1_dp
        real(dp) :: err
        character(len=256) :: msg

        ! Initialize linear operator and vectors.
        A = linop_rdp() ; call init_rand(A)
        allocate(b, c, cref) ; call init_rand(b) ; call c%zero()
        allocate(X(kdim+1)) ; call zero_basis(X)

        ! Reference: dense computation.
        cref%data = matmul(expm(tau*A%data), b%data)

        ! kexpm with internal workspace.
        call kexpm(c, A, b, tau, rtol_dp, info, kdim=kdim)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_rdp')
        call c%sub(cref) ; err = c%norm() / cref%norm()
        call get_err_str(msg, "max err: ", err)
        call check(error, info > 0 .and. err < rtol_dp)
        call check_test(error, 'test_kexptA_rdp', &
                              & info='kexpm', eq='c = expm(tau A) @ b', context=msg)

        ! krylov_exptA (adaptive sub-stepping).
        call krylov_exptA(c, A, b, tau, info)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_rdp')
        call c%sub(cref) ; err = c%norm() / cref%norm()
        call get_err_str(msg, "max err: ", err)
        call check(error, info > 0 .and. err < rtol_dp)
        call check_test(error, 'test_kexptA_rdp', &
                              & info='krylov_exptA', eq='c = expm(tau A) @ b', context=msg)

        ! kexpm with user-supplied workspace.
        call kexpm(c, A, b, tau, rtol_dp, info, X=X)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_rdp')
        call c%sub(cref) ; err = c%norm() / cref%norm()
        call get_err_str(msg, "max err: ", err)
        call check(error, info > 0 .and. err < rtol_dp)
        call check_test(error, 'test_kexptA_rdp', &
                              & info='kexpm, workspace', eq='c = expm(tau A) @ b', context=msg)

        ! Re-used (non-zeroed) workspace with kdim smaller than the workspace.
        call kexpm(c, A, b, tau, rtol_dp, info, kdim=kdim/2, X=X)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_rdp')
        call c%sub(cref) ; err = c%norm() / cref%norm()
        call get_err_str(msg, "max err: ", err)
        call check(error, info > 0 .and. info <= kdim/2 + 1 .and. err < rtol_dp)
        call check_test(error, 'test_kexptA_rdp', &
                              & info='kexpm, workspace + kdim', eq='c = expm(tau A) @ b', context=msg)

        return
    end subroutine test_kexptA_rdp

    subroutine test_kexptA_adjoint_rdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test linear operator.
        type(linop_rdp), allocatable :: A
        ! Test vectors.
        type(vector_rdp), allocatable :: b, c, cref
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(dp), parameter :: tau = 0.1_dp
        real(dp) :: err
        character(len=256) :: msg

        ! Initialize linear operator and vectors.
        A = linop_rdp() ; call init_rand(A)
        allocate(b, c, cref) ; call init_rand(b) ; call c%zero()

        ! Reference: dense computation with the adjoint operator.
        cref%data = matmul(expm(tau*hermitian(A%data)), b%data)

        ! kexpm.
        call kexpm(c, A, b, tau, rtol_dp, info, trans=.true.)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_adjoint_rdp')
        call c%sub(cref) ; err = c%norm() / cref%norm()
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_kexptA_adjoint_rdp', &
                              & info='kexpm', eq='c = expm(tau A.H) @ b', context=msg)

        ! krylov_exptA.
        call krylov_exptA(c, A, b, tau, info, trans=.true.)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_adjoint_rdp')
        call c%sub(cref) ; err = c%norm() / cref%norm()
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_kexptA_adjoint_rdp', &
                              & info='krylov_exptA', eq='c = expm(tau A.H) @ b', context=msg)

        return
    end subroutine test_kexptA_adjoint_rdp

    subroutine test_kexptA_zero_input_rdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test linear operator.
        type(linop_rdp), allocatable :: A
        ! Test vectors.
        integer, parameter :: p = 2
        type(vector_rdp), allocatable :: b, c, Bblk(:), Cblk(:)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(dp), parameter :: tau = 0.1_dp
        real(dp) :: err
        integer :: i
        character(len=256) :: msg

        ! Initialize linear operator, zero inputs and non-zero outputs.
        A = linop_rdp() ; call init_rand(A)
        allocate(b, c) ; call b%zero()
        allocate(Bblk(p), Cblk(p)) ; call zero_basis(Bblk)

        do
            ! kexpm.
            call c%rand()
            call kexpm(c, A, b, tau, rtol_dp, info)
            err = c%norm()
            call check(error, info > 0 .and. err == 0.0_dp)
            if (allocated(error)) exit

            ! krylov_exptA.
            call c%rand()
            call krylov_exptA(c, A, b, tau, info)
            err = c%norm()
            call check(error, info > 0 .and. err == 0.0_dp)
            if (allocated(error)) exit

            ! Block kexpm.
            call rand_basis(Cblk)
            call kexpm(Cblk, A, Bblk, tau, rtol_dp, info)
            err = 0.0_dp
            do i = 1, p
                err = err + Cblk(i)%norm()
            enddo
            call check(error, info > 0 .and. err == 0.0_dp)
            exit
        enddo
        call get_err_str(msg, "norm: ", err)
        call check_test(error, 'test_kexptA_zero_input_rdp', &
                              & info='Zero input', eq='expm(tau A) @ 0 = 0', context=msg)

        return
    end subroutine test_kexptA_zero_input_rdp

    subroutine test_kexptA_breakdown_rdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test linear operator.
        type(linop_rdp), allocatable :: A
        ! Test vectors.
        type(vector_rdp), allocatable :: b, c, cref
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        integer, parameter :: k_inv = 10
        real(dp), parameter :: tau = 0.1_dp
        real(dp) :: err
        character(len=256) :: msg

        ! Create block lower triangular matrix: span(e_1, ..., e_k_inv) is invariant.
        A = linop_rdp() ; call init_rand(A)
        A%data(k_inv+1:, :k_inv) = zero_rdp

        ! Input vector confined to the invariant subspace.
        allocate(b, c, cref) ; call init_rand(b) ; call c%zero()
        b%data(k_inv+1:) = zero_rdp

        ! Reference: dense computation.
        cref%data = matmul(expm(tau*A%data), b%data)

        ! kexpm with tol = 0: convergence is only possible through the Arnoldi breakdown,
        ! which must happen after exactly k_inv steps (the approximation is then exact).
        call kexpm(c, A, b, tau, 0.0_dp, info)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_breakdown_rdp')
        call check(error, info == k_inv)
        call check_test(error, 'test_kexptA_breakdown_rdp', &
                              & info='Subspace Dim', eq='info == k_inv', context='Invariant detection')

        call c%sub(cref) ; err = c%norm() / cref%norm()
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_kexptA_breakdown_rdp', &
                              & info='Invariant Property', eq='c = expm(tau A) @ b', context=msg)

        ! krylov_exptA: a single sub-step, never more than k_inv matrix-vector products.
        call krylov_exptA(c, A, b, tau, info)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_breakdown_rdp')
        call c%sub(cref) ; err = c%norm() / cref%norm()
        call get_err_str(msg, "max err: ", err)
        call check(error, info <= k_inv .and. err < rtol_dp)
        call check_test(error, 'test_kexptA_breakdown_rdp', &
                              & info='krylov_exptA', eq='<= k_inv matvecs, c = expm(tau A) @ b', context=msg)

        return
    end subroutine test_kexptA_breakdown_rdp

    subroutine test_kexptA_stiff_rdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test linear operator.
        type(linop_rdp), allocatable :: A
        ! Test vectors.
        type(vector_rdp), allocatable :: b, c, cref
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(dp), parameter :: tau = 1.0_dp, s = 30.0_dp
        real(dp) :: err
        character(len=256) :: msg

        ! Stable random operator with || tau A || ~ 60: far beyond what a single
        ! 30-vector Krylov basis can resolve.
        A = linop_rdp() ; call init_rand(A)
        A%data = s*A%data/sqrt(real(test_size, dp)) - s*eye(test_size, mold=1.0_dp)
        allocate(b, c, cref) ; call init_rand(b) ; call c%zero()

        ! Reference: dense computation.
        cref%data = matmul(expm(tau*A%data), b%data)

        ! krylov_exptA must converge thanks to the adaptive sub-stepping (> 30 matvecs).
        call krylov_exptA(c, A, b, tau, info)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_stiff_rdp')
        call c%sub(cref) ; err = c%norm() / cref%norm()
        call get_err_str(msg, "max err: ", err)
        call check(error, info > 30 .and. err < rtol_dp)
        call check_test(error, 'test_kexptA_stiff_rdp', &
                              & info='Sub-stepping', eq='c = expm(tau A) @ b', context=msg)

        return
    end subroutine test_kexptA_stiff_rdp

    subroutine test_kexptA_warm_start_rdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test linear operator.
        type(linop_rdp), allocatable :: A
        ! Test vectors.
        type(vector_rdp), allocatable :: b, c, cref
        ! Sub-stepping options and metadata.
        type(kexpm_dp_opts) :: opts
        type(kexpm_dp_metadata) :: meta
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        integer, parameter :: ncalls = 10
        real(dp), parameter :: tau = 1.0_dp, s = 30.0_dp
        real(dp), allocatable :: E(:, :)
        real(dp) :: err, tol
        integer :: icall, cold, last
        character(len=256) :: msg

        ! Stiff stable random operator.
        A = linop_rdp() ; call init_rand(A)
        A%data = s*A%data/sqrt(real(test_size, dp)) - s*eye(test_size, mold=1.0_dp)
        E = expm(tau*A%data)
        allocate(b, c, cref) ; call c%zero()

        ! Repeated calls with a persistent metadata object and new right-hand sides.
        err = 0.0_dp
        do icall = 1, ncalls
            call init_rand(b) ; cref%data = matmul(E, b%data) ; tol = atol_dp*b%norm()
            call kexpm_substep(c, A, b, tau, tol, info, meta=meta)
            call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_warm_start_rdp')
            if (icall == 1) cold = info
            last = info
            call c%sub(cref) ; err = max(err, c%norm() / cref%norm())
        enddo
        write(msg, '(A,E9.2,A,I0,A,I0,A,I0)') 'max err: ', err, ', cold: ', cold, ', last: ', last, ', nsub: ', meta%nsub
        call check(error, err < rtol_dp)
        call check_test(error, 'test_kexptA_warm_start_rdp', &
                              & info='Accuracy', eq='c = expm(tau A) @ b', context=msg)

        ! Settings learned, frozen, and not more expensive than the cold call (up to the
        ! call-to-call variation of the adaptive cost with the right-hand side).
        call check(error, meta%learned .and. meta%frozen .and. meta%n_calls == ncalls &
                        & .and. last <= cold + max(2, cold/10))
        call check_test(error, 'test_kexptA_warm_start_rdp', &
                              & info='Learned settings', eq='frozen, cost <= cold cost', context=msg)

        ! reset() keeps the learned settings, forget() clears them.
        do
            call meta%reset()
            call check(error, meta%learned .and. meta%n_calls == 0 .and. meta%n_matvec == 0)
            if (allocated(error)) exit
            call meta%forget()
            call check(error, .not. meta%learned .and. meta%nsub == 0)
            exit
        enddo
        call check_test(error, 'test_kexptA_warm_start_rdp', &
                              & info='reset / forget', eq='reset keeps, forget clears', context='')

        ! tune = .false. : the sub-step count of the cold call is kept.
        opts = kexpm_dp_opts(tune=.false.)
        call init_rand(b) ; cref%data = matmul(E, b%data) ; tol = atol_dp*b%norm()
        call kexpm_substep(c, A, b, tau, tol, info, options=opts, meta=meta)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_warm_start_rdp')
        call kexpm_substep(c, A, b, tau, tol, info, options=opts, meta=meta)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_warm_start_rdp')
        call c%sub(cref) ; err = c%norm() / cref%norm()
        write(msg, '(A,E9.2,A,I0,A,I0)') 'max err: ', err, ', nsub: ', meta%nsub, ', nsub0: ', meta%nsub0
        call check(error, err < rtol_dp .and. meta%frozen .and. meta%nsub == meta%nsub0)
        call check_test(error, 'test_kexptA_warm_start_rdp', &
                              & info='No tuning', eq='nsub == nsub0', context=msg)

        ! A different horizon invalidates the learned settings.
        cref%data = matmul(expm(0.5_dp*tau*A%data), b%data)
        call kexpm_substep(c, A, b, 0.5_dp*tau, tol, info, meta=meta)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_warm_start_rdp')
        call c%sub(cref) ; err = c%norm() / cref%norm()
        write(msg, '(A,E9.2,A,E9.2)') 'max err: ', err, ', tau_learned: ', meta%tau_learned
        call check(error, err < rtol_dp .and. abs(meta%tau_learned - 0.5_dp*tau) <= epsilon(tau))
        call check_test(error, 'test_kexptA_warm_start_rdp', &
                              & info='New horizon', eq='settings re-learned', context=msg)

        return
    end subroutine test_kexptA_warm_start_rdp

    subroutine test_block_kexptA_rdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test linear operator.
        type(linop_rdp), allocatable :: A
        ! Test vectors.
        integer, parameter :: p = 3
        type(vector_rdp), allocatable :: B(:), C(:)
        ! Krylov subspace (user-supplied workspace).
        integer, parameter :: kdim = 15
        type(vector_rdp), allocatable :: X(:)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(dp), parameter :: tau = 0.1_dp
        real(dp), allocatable :: Bdata(:, :), Cdata(:, :), Cref(:, :)
        real(dp) :: err
        integer :: i
        character(len=256) :: msg

        ! Initialize linear operator and input block.
        A = linop_rdp() ; call init_rand(A)
        allocate(B(p), C(p), X(p*(kdim+1))) ; call init_rand(B) ; call zero_basis(C) ; call zero_basis(X)
        allocate(Bdata(test_size, p), Cdata(test_size, p)) ; call get_data(Bdata, B)

        ! Reference: dense computation.
        Cref = matmul(expm(tau*A%data), Bdata)

        ! Column-by-column kexpm.
        do i = 1, p
            call kexpm(C(i), A, B(i), tau, rtol_dp, info, kdim=kdim)
            call check_info(info, 'kexpm', module=this_module_long, procedure='test_block_kexptA_rdp')
        enddo
        call get_data(Cdata, C)
        err = norm2(abs(Cdata - Cref)) / norm2(abs(Cref))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_block_kexptA_rdp', &
                              & info='Column-wise', eq='C = expm(tau A) @ B', context=msg)

        ! Block kexpm.
        call kexpm(C, A, B, tau, rtol_dp, info, kdim=kdim)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_block_kexptA_rdp')
        call get_data(Cdata, C)
        err = norm2(abs(Cdata - Cref)) / norm2(abs(Cref))
        call get_err_str(msg, "max err: ", err)
        call check(error, info > 0 .and. info <= p*(kdim+1) .and. err < rtol_dp)
        call check_test(error, 'test_block_kexptA_rdp', &
                              & info='Block', eq='C = expm(tau A) @ B', context=msg)

        ! Block kexpm with user-supplied workspace.
        call kexpm(C, A, B, tau, rtol_dp, info, X=X)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_block_kexptA_rdp')
        call get_data(Cdata, C)
        err = norm2(abs(Cdata - Cref)) / norm2(abs(Cref))
        call get_err_str(msg, "max err: ", err)
        call check(error, info > 0 .and. err < rtol_dp)
        call check_test(error, 'test_block_kexptA_rdp', &
                              & info='Block, workspace', eq='C = expm(tau A) @ B', context=msg)

        block
        integer, parameter :: pinv = 2, k_inv = 2*pinv
        type(vector_rdp), allocatable :: B2(:), C2(:)
        real(dp), allocatable :: B2data(:, :), C2data(:, :), C2ref(:, :)

        ! Create block lower triangular matrix: span(e_1, ..., e_k_inv) is invariant.
        call init_rand(A)
        A%data(k_inv+1:, :k_inv) = zero_rdp

        ! Input block confined to the invariant subspace.
        allocate(B2(pinv), C2(pinv)) ; call init_rand(B2) ; call zero_basis(C2)
        do i = 1, pinv
            B2(i)%data(k_inv+1:) = zero_rdp
        enddo
        allocate(B2data(test_size, pinv), C2data(test_size, pinv)) ; call get_data(B2data, B2)
        C2ref = matmul(expm(tau*A%data), B2data)

        ! tol = 0: block Arnoldi must break down after k_inv/pinv block steps (exact approximation).
        call kexpm(C2, A, B2, tau, 0.0_dp, info)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_block_kexptA_rdp')
        call check(error, info == k_inv)
        call check_test(error, 'test_block_kexptA_rdp', &
                              & info='Subspace Dim', eq='info == k_inv', context='Invariant detection')

        call get_data(C2data, C2)
        err = norm2(abs(C2data - C2ref)) / norm2(abs(C2ref))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_block_kexptA_rdp', &
                              & info='Invariant Property', eq='C = expm(tau A) @ B', context=msg)
        end block

        return
    end subroutine test_block_kexptA_rdp

    subroutine collect_expm_csp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                        new_unittest("Krylov expm invalid inputs", test_kexptA_invalid_inputs_csp), &
                        new_unittest("Dense expm", test_dense_expm_csp), &
                        new_unittest("Krylov expm", test_kexptA_csp), &
                        new_unittest("Krylov expm adjoint", test_kexptA_adjoint_csp), &
                        new_unittest("Krylov expm zero input", test_kexptA_zero_input_csp), &
                        new_unittest("Krylov expm invariant subspace", test_kexptA_breakdown_csp), &
                        new_unittest("Krylov expm stiff (sub-stepping)", test_kexptA_stiff_csp), &
                        new_unittest("Krylov expm warm start", test_kexptA_warm_start_csp), &
                        new_unittest("Block Krylov expm", test_block_kexptA_csp) &
                    ]

        return
    end subroutine collect_expm_csp_testsuite

    subroutine test_kexptA_invalid_inputs_csp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test linear operator.
        type(linop_csp), allocatable :: A
        ! Test vectors.
        type(vector_csp), allocatable :: b, c, Bblk(:), Cblk(:), X(:)
        ! Sub-stepping options and metadata.
        type(kexpm_sp_opts) :: opts
        type(kexpm_sp_metadata) :: meta
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        integer, parameter :: p = 2
        real(sp), parameter :: tau = 0.1_sp
        integer :: icase
        character(len=256) :: msg

        ! Common initialization.
        A = linop_csp() ; call init_rand(A)
        allocate(b, c) ; call init_rand(b) ; call c%zero()
        allocate(Bblk(p)) ; call init_rand(Bblk)

        test_loop: do
            !-----     Single vector kexpm     -----

            ! --- Case 1: workspace with fewer than 2 vectors (info = -3) ---
            icase = 1
            allocate(X(1)) ; call zero_basis(X)
            call kexpm(c, A, b, tau, rtol_sp, info, X=X)
            call check(error, info == -3)
            if (allocated(error)) exit test_loop

            ! --- Case 2: kdim < 1 (info = -4) ---
            icase = 2
            call kexpm(c, A, b, tau, rtol_sp, info, kdim=0)
            call check(error, info == -4)
            if (allocated(error)) exit test_loop

            ! --- Case 3: negative tolerance (info = -5) ---
            icase = 3
            call kexpm(c, A, b, tau, -1.0_sp, info)
            call check(error, info == -5)
            if (allocated(error)) exit test_loop

            !-----     Block kexpm     -----

            ! --- Case 4: size(C) /= size(B) (info = -6) ---
            icase = 4
            allocate(Cblk(p+1)) ; call zero_basis(Cblk)
            call kexpm(Cblk, A, Bblk, tau, rtol_sp, info)
            call check(error, info == -6)
            if (allocated(error)) exit test_loop

            ! --- Case 5: empty block (info = -6) ---
            icase = 5
            block
                type(vector_csp), allocatable :: B0(:), C0(:)
                allocate(B0(0), C0(0))
                call kexpm(C0, A, B0, tau, rtol_sp, info)
            end block
            call check(error, info == -6)
            if (allocated(error)) exit test_loop

            ! --- Case 6: workspace size not a multiple of the block size (info = -3) ---
            icase = 6
            deallocate(Cblk, X) ; allocate(Cblk(p), X(2*p+1)) ; call zero_basis(Cblk) ; call zero_basis(X)
            call kexpm(Cblk, A, Bblk, tau, rtol_sp, info, X=X)
            call check(error, info == -3)
            if (allocated(error)) exit test_loop

            ! --- Case 7: workspace with fewer than 2*p vectors (info = -3) ---
            icase = 7
            deallocate(X) ; allocate(X(p)) ; call zero_basis(X)
            call kexpm(Cblk, A, Bblk, tau, rtol_sp, info, X=X)
            call check(error, info == -3)
            if (allocated(error)) exit test_loop

            ! --- Case 8: kdim < 1 (info = -4) ---
            icase = 8
            call kexpm(Cblk, A, Bblk, tau, rtol_sp, info, kdim=0)
            call check(error, info == -4)
            if (allocated(error)) exit test_loop

            ! --- Case 9: negative tolerance (info = -5) ---
            icase = 9
            call kexpm(Cblk, A, Bblk, tau, -1.0_sp, info)
            call check(error, info == -5)
            if (allocated(error)) exit test_loop

            !-----     Sub-stepped kexpm_substep     -----

            ! --- Case 10: opts%kdim < 1 (info = -4) ---
            icase = 10
            opts = kexpm_sp_opts(kdim=0)
            call kexpm_substep(c, A, b, tau, rtol_sp, info, options=opts)
            call check(error, info == -4)
            if (allocated(error)) exit test_loop

            ! --- Case 11: negative tolerance (info = -5) ---
            icase = 11
            call kexpm_substep(c, A, b, tau, -1.0_sp, info)
            call check(error, info == -5)
            if (allocated(error)) exit test_loop

            ! --- Case 12: opts%max_rejections < 1 (info = -7) ---
            icase = 12
            opts = kexpm_sp_opts(max_rejections=0)
            call kexpm_substep(c, A, b, tau, rtol_sp, info, options=opts, meta=meta)
            call check(error, info == -7)
            if (allocated(error)) exit test_loop

            ! --- Case 13: opts%max_substeps < 1 (info = -7), metadata left untouched ---
            icase = 13
            opts = kexpm_sp_opts(max_substeps=0)
            call kexpm_substep(c, A, b, tau, rtol_sp, info, options=opts, meta=meta)
            call check(error, info == -7 .and. meta%n_calls == 0 .and. .not. meta%learned)
            exit test_loop
        end do test_loop
        write(msg, '(A,I0)') 'last case: ', icase
        call check_test(error, 'test_kexptA_invalid_inputs_csp', &
                              & info='Invalid input parameters', eq='', context=msg)

        return
    end subroutine test_kexptA_invalid_inputs_csp

    subroutine test_dense_expm_csp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Problem dimension.
        integer, parameter :: n = 5, m = 6
        ! Test matrix.
        complex(sp) :: A(n, n), E(n, n), Eref(n, n)
        ! Miscellaneous.
        integer :: i, j
        real(sp) :: err
        character(len=256) :: msg

        ! Initialize nilpotent matrix (the exponential is a finite series).
        A = 0.0_sp
        do i = 1, n-1
            A(i, i+1) = m*1.0_sp
        enddo

        ! Reference with analytical exponential: Eref(i, i+j) = m**j / j!.
        Eref = eye(n, mold=1.0_sp)
        do i = 1, n-1
            do j = 1, n-i
                Eref(i, i+j) = Eref(i, i+j-1)*m/j
            enddo
        enddo

        ! Compute matrix exponential.
        E = expm(A)

        ! Check result.
        err = maxval(abs(E-Eref))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_dense_expm_csp', &
                              & info='Nilpotent matrix', eq='expm(A) = Eref', context=msg)

        return
    end subroutine test_dense_expm_csp

    subroutine test_kexptA_csp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test linear operator.
        type(linop_csp), allocatable :: A
        ! Test vectors.
        type(vector_csp), allocatable :: b, c, cref
        ! Krylov subspace (user-supplied workspace).
        integer, parameter :: kdim = 64
        type(vector_csp), allocatable :: X(:)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(sp), parameter :: tau = 0.1_sp
        real(sp) :: err
        character(len=256) :: msg

        ! Initialize linear operator and vectors.
        A = linop_csp() ; call init_rand(A)
        allocate(b, c, cref) ; call init_rand(b) ; call c%zero()
        allocate(X(kdim+1)) ; call zero_basis(X)

        ! Reference: dense computation.
        cref%data = matmul(expm(tau*A%data), b%data)

        ! kexpm with internal workspace.
        call kexpm(c, A, b, tau, rtol_sp, info, kdim=kdim)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_csp')
        call c%sub(cref) ; err = c%norm() / cref%norm()
        call get_err_str(msg, "max err: ", err)
        call check(error, info > 0 .and. err < rtol_sp)
        call check_test(error, 'test_kexptA_csp', &
                              & info='kexpm', eq='c = expm(tau A) @ b', context=msg)

        ! krylov_exptA (adaptive sub-stepping).
        call krylov_exptA(c, A, b, tau, info)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_csp')
        call c%sub(cref) ; err = c%norm() / cref%norm()
        call get_err_str(msg, "max err: ", err)
        call check(error, info > 0 .and. err < rtol_sp)
        call check_test(error, 'test_kexptA_csp', &
                              & info='krylov_exptA', eq='c = expm(tau A) @ b', context=msg)

        ! kexpm with user-supplied workspace.
        call kexpm(c, A, b, tau, rtol_sp, info, X=X)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_csp')
        call c%sub(cref) ; err = c%norm() / cref%norm()
        call get_err_str(msg, "max err: ", err)
        call check(error, info > 0 .and. err < rtol_sp)
        call check_test(error, 'test_kexptA_csp', &
                              & info='kexpm, workspace', eq='c = expm(tau A) @ b', context=msg)

        ! Re-used (non-zeroed) workspace with kdim smaller than the workspace.
        call kexpm(c, A, b, tau, rtol_sp, info, kdim=kdim/2, X=X)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_csp')
        call c%sub(cref) ; err = c%norm() / cref%norm()
        call get_err_str(msg, "max err: ", err)
        call check(error, info > 0 .and. info <= kdim/2 + 1 .and. err < rtol_sp)
        call check_test(error, 'test_kexptA_csp', &
                              & info='kexpm, workspace + kdim', eq='c = expm(tau A) @ b', context=msg)

        return
    end subroutine test_kexptA_csp

    subroutine test_kexptA_adjoint_csp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test linear operator.
        type(linop_csp), allocatable :: A
        ! Test vectors.
        type(vector_csp), allocatable :: b, c, cref
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(sp), parameter :: tau = 0.1_sp
        real(sp) :: err
        character(len=256) :: msg

        ! Initialize linear operator and vectors.
        A = linop_csp() ; call init_rand(A)
        allocate(b, c, cref) ; call init_rand(b) ; call c%zero()

        ! Reference: dense computation with the adjoint operator.
        cref%data = matmul(expm(tau*hermitian(A%data)), b%data)

        ! kexpm.
        call kexpm(c, A, b, tau, rtol_sp, info, trans=.true.)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_adjoint_csp')
        call c%sub(cref) ; err = c%norm() / cref%norm()
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_kexptA_adjoint_csp', &
                              & info='kexpm', eq='c = expm(tau A.H) @ b', context=msg)

        ! krylov_exptA.
        call krylov_exptA(c, A, b, tau, info, trans=.true.)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_adjoint_csp')
        call c%sub(cref) ; err = c%norm() / cref%norm()
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_kexptA_adjoint_csp', &
                              & info='krylov_exptA', eq='c = expm(tau A.H) @ b', context=msg)

        return
    end subroutine test_kexptA_adjoint_csp

    subroutine test_kexptA_zero_input_csp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test linear operator.
        type(linop_csp), allocatable :: A
        ! Test vectors.
        integer, parameter :: p = 2
        type(vector_csp), allocatable :: b, c, Bblk(:), Cblk(:)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(sp), parameter :: tau = 0.1_sp
        real(sp) :: err
        integer :: i
        character(len=256) :: msg

        ! Initialize linear operator, zero inputs and non-zero outputs.
        A = linop_csp() ; call init_rand(A)
        allocate(b, c) ; call b%zero()
        allocate(Bblk(p), Cblk(p)) ; call zero_basis(Bblk)

        do
            ! kexpm.
            call c%rand()
            call kexpm(c, A, b, tau, rtol_sp, info)
            err = c%norm()
            call check(error, info > 0 .and. err == 0.0_sp)
            if (allocated(error)) exit

            ! krylov_exptA.
            call c%rand()
            call krylov_exptA(c, A, b, tau, info)
            err = c%norm()
            call check(error, info > 0 .and. err == 0.0_sp)
            if (allocated(error)) exit

            ! Block kexpm.
            call rand_basis(Cblk)
            call kexpm(Cblk, A, Bblk, tau, rtol_sp, info)
            err = 0.0_sp
            do i = 1, p
                err = err + Cblk(i)%norm()
            enddo
            call check(error, info > 0 .and. err == 0.0_sp)
            exit
        enddo
        call get_err_str(msg, "norm: ", err)
        call check_test(error, 'test_kexptA_zero_input_csp', &
                              & info='Zero input', eq='expm(tau A) @ 0 = 0', context=msg)

        return
    end subroutine test_kexptA_zero_input_csp

    subroutine test_kexptA_breakdown_csp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test linear operator.
        type(linop_csp), allocatable :: A
        ! Test vectors.
        type(vector_csp), allocatable :: b, c, cref
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        integer, parameter :: k_inv = 10
        real(sp), parameter :: tau = 0.1_sp
        real(sp) :: err
        character(len=256) :: msg

        ! Create block lower triangular matrix: span(e_1, ..., e_k_inv) is invariant.
        A = linop_csp() ; call init_rand(A)
        A%data(k_inv+1:, :k_inv) = zero_csp

        ! Input vector confined to the invariant subspace.
        allocate(b, c, cref) ; call init_rand(b) ; call c%zero()
        b%data(k_inv+1:) = zero_csp

        ! Reference: dense computation.
        cref%data = matmul(expm(tau*A%data), b%data)

        ! kexpm with tol = 0: convergence is only possible through the Arnoldi breakdown,
        ! which must happen after exactly k_inv steps (the approximation is then exact).
        call kexpm(c, A, b, tau, 0.0_sp, info)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_breakdown_csp')
        call check(error, info == k_inv)
        call check_test(error, 'test_kexptA_breakdown_csp', &
                              & info='Subspace Dim', eq='info == k_inv', context='Invariant detection')

        call c%sub(cref) ; err = c%norm() / cref%norm()
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_kexptA_breakdown_csp', &
                              & info='Invariant Property', eq='c = expm(tau A) @ b', context=msg)

        ! krylov_exptA: a single sub-step, never more than k_inv matrix-vector products.
        call krylov_exptA(c, A, b, tau, info)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_breakdown_csp')
        call c%sub(cref) ; err = c%norm() / cref%norm()
        call get_err_str(msg, "max err: ", err)
        call check(error, info <= k_inv .and. err < rtol_sp)
        call check_test(error, 'test_kexptA_breakdown_csp', &
                              & info='krylov_exptA', eq='<= k_inv matvecs, c = expm(tau A) @ b', context=msg)

        return
    end subroutine test_kexptA_breakdown_csp

    subroutine test_kexptA_stiff_csp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test linear operator.
        type(linop_csp), allocatable :: A
        ! Test vectors.
        type(vector_csp), allocatable :: b, c, cref
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(sp), parameter :: tau = 1.0_sp, s = 30.0_sp
        real(sp) :: err
        character(len=256) :: msg

        ! Stable random operator with || tau A || ~ 60: far beyond what a single
        ! 30-vector Krylov basis can resolve.
        A = linop_csp() ; call init_rand(A)
        A%data = s*A%data/sqrt(real(test_size, sp)) - s*eye(test_size, mold=1.0_sp)
        allocate(b, c, cref) ; call init_rand(b) ; call c%zero()

        ! Reference: dense computation.
        cref%data = matmul(expm(tau*A%data), b%data)

        ! krylov_exptA must converge thanks to the adaptive sub-stepping (> 30 matvecs).
        call krylov_exptA(c, A, b, tau, info)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_stiff_csp')
        call c%sub(cref) ; err = c%norm() / cref%norm()
        call get_err_str(msg, "max err: ", err)
        call check(error, info > 30 .and. err < rtol_sp)
        call check_test(error, 'test_kexptA_stiff_csp', &
                              & info='Sub-stepping', eq='c = expm(tau A) @ b', context=msg)

        return
    end subroutine test_kexptA_stiff_csp

    subroutine test_kexptA_warm_start_csp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test linear operator.
        type(linop_csp), allocatable :: A
        ! Test vectors.
        type(vector_csp), allocatable :: b, c, cref
        ! Sub-stepping options and metadata.
        type(kexpm_sp_opts) :: opts
        type(kexpm_sp_metadata) :: meta
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        integer, parameter :: ncalls = 10
        real(sp), parameter :: tau = 1.0_sp, s = 30.0_sp
        complex(sp), allocatable :: E(:, :)
        real(sp) :: err, tol
        integer :: icall, cold, last
        character(len=256) :: msg

        ! Stiff stable random operator.
        A = linop_csp() ; call init_rand(A)
        A%data = s*A%data/sqrt(real(test_size, sp)) - s*eye(test_size, mold=1.0_sp)
        E = expm(tau*A%data)
        allocate(b, c, cref) ; call c%zero()

        ! Repeated calls with a persistent metadata object and new right-hand sides.
        err = 0.0_sp
        do icall = 1, ncalls
            call init_rand(b) ; cref%data = matmul(E, b%data) ; tol = atol_sp*b%norm()
            call kexpm_substep(c, A, b, tau, tol, info, meta=meta)
            call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_warm_start_csp')
            if (icall == 1) cold = info
            last = info
            call c%sub(cref) ; err = max(err, c%norm() / cref%norm())
        enddo
        write(msg, '(A,E9.2,A,I0,A,I0,A,I0)') 'max err: ', err, ', cold: ', cold, ', last: ', last, ', nsub: ', meta%nsub
        call check(error, err < rtol_sp)
        call check_test(error, 'test_kexptA_warm_start_csp', &
                              & info='Accuracy', eq='c = expm(tau A) @ b', context=msg)

        ! Settings learned, frozen, and not more expensive than the cold call (up to the
        ! call-to-call variation of the adaptive cost with the right-hand side).
        call check(error, meta%learned .and. meta%frozen .and. meta%n_calls == ncalls &
                        & .and. last <= cold + max(2, cold/10))
        call check_test(error, 'test_kexptA_warm_start_csp', &
                              & info='Learned settings', eq='frozen, cost <= cold cost', context=msg)

        ! reset() keeps the learned settings, forget() clears them.
        do
            call meta%reset()
            call check(error, meta%learned .and. meta%n_calls == 0 .and. meta%n_matvec == 0)
            if (allocated(error)) exit
            call meta%forget()
            call check(error, .not. meta%learned .and. meta%nsub == 0)
            exit
        enddo
        call check_test(error, 'test_kexptA_warm_start_csp', &
                              & info='reset / forget', eq='reset keeps, forget clears', context='')

        ! tune = .false. : the sub-step count of the cold call is kept.
        opts = kexpm_sp_opts(tune=.false.)
        call init_rand(b) ; cref%data = matmul(E, b%data) ; tol = atol_sp*b%norm()
        call kexpm_substep(c, A, b, tau, tol, info, options=opts, meta=meta)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_warm_start_csp')
        call kexpm_substep(c, A, b, tau, tol, info, options=opts, meta=meta)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_warm_start_csp')
        call c%sub(cref) ; err = c%norm() / cref%norm()
        write(msg, '(A,E9.2,A,I0,A,I0)') 'max err: ', err, ', nsub: ', meta%nsub, ', nsub0: ', meta%nsub0
        call check(error, err < rtol_sp .and. meta%frozen .and. meta%nsub == meta%nsub0)
        call check_test(error, 'test_kexptA_warm_start_csp', &
                              & info='No tuning', eq='nsub == nsub0', context=msg)

        ! A different horizon invalidates the learned settings.
        cref%data = matmul(expm(0.5_sp*tau*A%data), b%data)
        call kexpm_substep(c, A, b, 0.5_sp*tau, tol, info, meta=meta)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_warm_start_csp')
        call c%sub(cref) ; err = c%norm() / cref%norm()
        write(msg, '(A,E9.2,A,E9.2)') 'max err: ', err, ', tau_learned: ', meta%tau_learned
        call check(error, err < rtol_sp .and. abs(meta%tau_learned - 0.5_sp*tau) <= epsilon(tau))
        call check_test(error, 'test_kexptA_warm_start_csp', &
                              & info='New horizon', eq='settings re-learned', context=msg)

        return
    end subroutine test_kexptA_warm_start_csp

    subroutine test_block_kexptA_csp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test linear operator.
        type(linop_csp), allocatable :: A
        ! Test vectors.
        integer, parameter :: p = 3
        type(vector_csp), allocatable :: B(:), C(:)
        ! Krylov subspace (user-supplied workspace).
        integer, parameter :: kdim = 15
        type(vector_csp), allocatable :: X(:)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(sp), parameter :: tau = 0.1_sp
        complex(sp), allocatable :: Bdata(:, :), Cdata(:, :), Cref(:, :)
        real(sp) :: err
        integer :: i
        character(len=256) :: msg

        ! Initialize linear operator and input block.
        A = linop_csp() ; call init_rand(A)
        allocate(B(p), C(p), X(p*(kdim+1))) ; call init_rand(B) ; call zero_basis(C) ; call zero_basis(X)
        allocate(Bdata(test_size, p), Cdata(test_size, p)) ; call get_data(Bdata, B)

        ! Reference: dense computation.
        Cref = matmul(expm(tau*A%data), Bdata)

        ! Column-by-column kexpm.
        do i = 1, p
            call kexpm(C(i), A, B(i), tau, rtol_sp, info, kdim=kdim)
            call check_info(info, 'kexpm', module=this_module_long, procedure='test_block_kexptA_csp')
        enddo
        call get_data(Cdata, C)
        err = norm2(abs(Cdata - Cref)) / norm2(abs(Cref))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_block_kexptA_csp', &
                              & info='Column-wise', eq='C = expm(tau A) @ B', context=msg)

        ! Block kexpm.
        call kexpm(C, A, B, tau, rtol_sp, info, kdim=kdim)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_block_kexptA_csp')
        call get_data(Cdata, C)
        err = norm2(abs(Cdata - Cref)) / norm2(abs(Cref))
        call get_err_str(msg, "max err: ", err)
        call check(error, info > 0 .and. info <= p*(kdim+1) .and. err < rtol_sp)
        call check_test(error, 'test_block_kexptA_csp', &
                              & info='Block', eq='C = expm(tau A) @ B', context=msg)

        ! Block kexpm with user-supplied workspace.
        call kexpm(C, A, B, tau, rtol_sp, info, X=X)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_block_kexptA_csp')
        call get_data(Cdata, C)
        err = norm2(abs(Cdata - Cref)) / norm2(abs(Cref))
        call get_err_str(msg, "max err: ", err)
        call check(error, info > 0 .and. err < rtol_sp)
        call check_test(error, 'test_block_kexptA_csp', &
                              & info='Block, workspace', eq='C = expm(tau A) @ B', context=msg)

        block
        integer, parameter :: pinv = 2, k_inv = 2*pinv
        type(vector_csp), allocatable :: B2(:), C2(:)
        complex(sp), allocatable :: B2data(:, :), C2data(:, :), C2ref(:, :)

        ! Create block lower triangular matrix: span(e_1, ..., e_k_inv) is invariant.
        call init_rand(A)
        A%data(k_inv+1:, :k_inv) = zero_csp

        ! Input block confined to the invariant subspace.
        allocate(B2(pinv), C2(pinv)) ; call init_rand(B2) ; call zero_basis(C2)
        do i = 1, pinv
            B2(i)%data(k_inv+1:) = zero_csp
        enddo
        allocate(B2data(test_size, pinv), C2data(test_size, pinv)) ; call get_data(B2data, B2)
        C2ref = matmul(expm(tau*A%data), B2data)

        ! tol = 0: block Arnoldi must break down after k_inv/pinv block steps (exact approximation).
        call kexpm(C2, A, B2, tau, 0.0_sp, info)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_block_kexptA_csp')
        call check(error, info == k_inv)
        call check_test(error, 'test_block_kexptA_csp', &
                              & info='Subspace Dim', eq='info == k_inv', context='Invariant detection')

        call get_data(C2data, C2)
        err = norm2(abs(C2data - C2ref)) / norm2(abs(C2ref))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_sp)
        call check_test(error, 'test_block_kexptA_csp', &
                              & info='Invariant Property', eq='C = expm(tau A) @ B', context=msg)
        end block

        return
    end subroutine test_block_kexptA_csp

    subroutine collect_expm_cdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                        new_unittest("Krylov expm invalid inputs", test_kexptA_invalid_inputs_cdp), &
                        new_unittest("Dense expm", test_dense_expm_cdp), &
                        new_unittest("Krylov expm", test_kexptA_cdp), &
                        new_unittest("Krylov expm adjoint", test_kexptA_adjoint_cdp), &
                        new_unittest("Krylov expm zero input", test_kexptA_zero_input_cdp), &
                        new_unittest("Krylov expm invariant subspace", test_kexptA_breakdown_cdp), &
                        new_unittest("Krylov expm stiff (sub-stepping)", test_kexptA_stiff_cdp), &
                        new_unittest("Krylov expm warm start", test_kexptA_warm_start_cdp), &
                        new_unittest("Block Krylov expm", test_block_kexptA_cdp) &
                    ]

        return
    end subroutine collect_expm_cdp_testsuite

    subroutine test_kexptA_invalid_inputs_cdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test linear operator.
        type(linop_cdp), allocatable :: A
        ! Test vectors.
        type(vector_cdp), allocatable :: b, c, Bblk(:), Cblk(:), X(:)
        ! Sub-stepping options and metadata.
        type(kexpm_dp_opts) :: opts
        type(kexpm_dp_metadata) :: meta
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        integer, parameter :: p = 2
        real(dp), parameter :: tau = 0.1_dp
        integer :: icase
        character(len=256) :: msg

        ! Common initialization.
        A = linop_cdp() ; call init_rand(A)
        allocate(b, c) ; call init_rand(b) ; call c%zero()
        allocate(Bblk(p)) ; call init_rand(Bblk)

        test_loop: do
            !-----     Single vector kexpm     -----

            ! --- Case 1: workspace with fewer than 2 vectors (info = -3) ---
            icase = 1
            allocate(X(1)) ; call zero_basis(X)
            call kexpm(c, A, b, tau, rtol_dp, info, X=X)
            call check(error, info == -3)
            if (allocated(error)) exit test_loop

            ! --- Case 2: kdim < 1 (info = -4) ---
            icase = 2
            call kexpm(c, A, b, tau, rtol_dp, info, kdim=0)
            call check(error, info == -4)
            if (allocated(error)) exit test_loop

            ! --- Case 3: negative tolerance (info = -5) ---
            icase = 3
            call kexpm(c, A, b, tau, -1.0_dp, info)
            call check(error, info == -5)
            if (allocated(error)) exit test_loop

            !-----     Block kexpm     -----

            ! --- Case 4: size(C) /= size(B) (info = -6) ---
            icase = 4
            allocate(Cblk(p+1)) ; call zero_basis(Cblk)
            call kexpm(Cblk, A, Bblk, tau, rtol_dp, info)
            call check(error, info == -6)
            if (allocated(error)) exit test_loop

            ! --- Case 5: empty block (info = -6) ---
            icase = 5
            block
                type(vector_cdp), allocatable :: B0(:), C0(:)
                allocate(B0(0), C0(0))
                call kexpm(C0, A, B0, tau, rtol_dp, info)
            end block
            call check(error, info == -6)
            if (allocated(error)) exit test_loop

            ! --- Case 6: workspace size not a multiple of the block size (info = -3) ---
            icase = 6
            deallocate(Cblk, X) ; allocate(Cblk(p), X(2*p+1)) ; call zero_basis(Cblk) ; call zero_basis(X)
            call kexpm(Cblk, A, Bblk, tau, rtol_dp, info, X=X)
            call check(error, info == -3)
            if (allocated(error)) exit test_loop

            ! --- Case 7: workspace with fewer than 2*p vectors (info = -3) ---
            icase = 7
            deallocate(X) ; allocate(X(p)) ; call zero_basis(X)
            call kexpm(Cblk, A, Bblk, tau, rtol_dp, info, X=X)
            call check(error, info == -3)
            if (allocated(error)) exit test_loop

            ! --- Case 8: kdim < 1 (info = -4) ---
            icase = 8
            call kexpm(Cblk, A, Bblk, tau, rtol_dp, info, kdim=0)
            call check(error, info == -4)
            if (allocated(error)) exit test_loop

            ! --- Case 9: negative tolerance (info = -5) ---
            icase = 9
            call kexpm(Cblk, A, Bblk, tau, -1.0_dp, info)
            call check(error, info == -5)
            if (allocated(error)) exit test_loop

            !-----     Sub-stepped kexpm_substep     -----

            ! --- Case 10: opts%kdim < 1 (info = -4) ---
            icase = 10
            opts = kexpm_dp_opts(kdim=0)
            call kexpm_substep(c, A, b, tau, rtol_dp, info, options=opts)
            call check(error, info == -4)
            if (allocated(error)) exit test_loop

            ! --- Case 11: negative tolerance (info = -5) ---
            icase = 11
            call kexpm_substep(c, A, b, tau, -1.0_dp, info)
            call check(error, info == -5)
            if (allocated(error)) exit test_loop

            ! --- Case 12: opts%max_rejections < 1 (info = -7) ---
            icase = 12
            opts = kexpm_dp_opts(max_rejections=0)
            call kexpm_substep(c, A, b, tau, rtol_dp, info, options=opts, meta=meta)
            call check(error, info == -7)
            if (allocated(error)) exit test_loop

            ! --- Case 13: opts%max_substeps < 1 (info = -7), metadata left untouched ---
            icase = 13
            opts = kexpm_dp_opts(max_substeps=0)
            call kexpm_substep(c, A, b, tau, rtol_dp, info, options=opts, meta=meta)
            call check(error, info == -7 .and. meta%n_calls == 0 .and. .not. meta%learned)
            exit test_loop
        end do test_loop
        write(msg, '(A,I0)') 'last case: ', icase
        call check_test(error, 'test_kexptA_invalid_inputs_cdp', &
                              & info='Invalid input parameters', eq='', context=msg)

        return
    end subroutine test_kexptA_invalid_inputs_cdp

    subroutine test_dense_expm_cdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Problem dimension.
        integer, parameter :: n = 5, m = 6
        ! Test matrix.
        complex(dp) :: A(n, n), E(n, n), Eref(n, n)
        ! Miscellaneous.
        integer :: i, j
        real(dp) :: err
        character(len=256) :: msg

        ! Initialize nilpotent matrix (the exponential is a finite series).
        A = 0.0_dp
        do i = 1, n-1
            A(i, i+1) = m*1.0_dp
        enddo

        ! Reference with analytical exponential: Eref(i, i+j) = m**j / j!.
        Eref = eye(n, mold=1.0_dp)
        do i = 1, n-1
            do j = 1, n-i
                Eref(i, i+j) = Eref(i, i+j-1)*m/j
            enddo
        enddo

        ! Compute matrix exponential.
        E = expm(A)

        ! Check result.
        err = maxval(abs(E-Eref))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_dense_expm_cdp', &
                              & info='Nilpotent matrix', eq='expm(A) = Eref', context=msg)

        return
    end subroutine test_dense_expm_cdp

    subroutine test_kexptA_cdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test linear operator.
        type(linop_cdp), allocatable :: A
        ! Test vectors.
        type(vector_cdp), allocatable :: b, c, cref
        ! Krylov subspace (user-supplied workspace).
        integer, parameter :: kdim = 64
        type(vector_cdp), allocatable :: X(:)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(dp), parameter :: tau = 0.1_dp
        real(dp) :: err
        character(len=256) :: msg

        ! Initialize linear operator and vectors.
        A = linop_cdp() ; call init_rand(A)
        allocate(b, c, cref) ; call init_rand(b) ; call c%zero()
        allocate(X(kdim+1)) ; call zero_basis(X)

        ! Reference: dense computation.
        cref%data = matmul(expm(tau*A%data), b%data)

        ! kexpm with internal workspace.
        call kexpm(c, A, b, tau, rtol_dp, info, kdim=kdim)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_cdp')
        call c%sub(cref) ; err = c%norm() / cref%norm()
        call get_err_str(msg, "max err: ", err)
        call check(error, info > 0 .and. err < rtol_dp)
        call check_test(error, 'test_kexptA_cdp', &
                              & info='kexpm', eq='c = expm(tau A) @ b', context=msg)

        ! krylov_exptA (adaptive sub-stepping).
        call krylov_exptA(c, A, b, tau, info)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_cdp')
        call c%sub(cref) ; err = c%norm() / cref%norm()
        call get_err_str(msg, "max err: ", err)
        call check(error, info > 0 .and. err < rtol_dp)
        call check_test(error, 'test_kexptA_cdp', &
                              & info='krylov_exptA', eq='c = expm(tau A) @ b', context=msg)

        ! kexpm with user-supplied workspace.
        call kexpm(c, A, b, tau, rtol_dp, info, X=X)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_cdp')
        call c%sub(cref) ; err = c%norm() / cref%norm()
        call get_err_str(msg, "max err: ", err)
        call check(error, info > 0 .and. err < rtol_dp)
        call check_test(error, 'test_kexptA_cdp', &
                              & info='kexpm, workspace', eq='c = expm(tau A) @ b', context=msg)

        ! Re-used (non-zeroed) workspace with kdim smaller than the workspace.
        call kexpm(c, A, b, tau, rtol_dp, info, kdim=kdim/2, X=X)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_cdp')
        call c%sub(cref) ; err = c%norm() / cref%norm()
        call get_err_str(msg, "max err: ", err)
        call check(error, info > 0 .and. info <= kdim/2 + 1 .and. err < rtol_dp)
        call check_test(error, 'test_kexptA_cdp', &
                              & info='kexpm, workspace + kdim', eq='c = expm(tau A) @ b', context=msg)

        return
    end subroutine test_kexptA_cdp

    subroutine test_kexptA_adjoint_cdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test linear operator.
        type(linop_cdp), allocatable :: A
        ! Test vectors.
        type(vector_cdp), allocatable :: b, c, cref
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(dp), parameter :: tau = 0.1_dp
        real(dp) :: err
        character(len=256) :: msg

        ! Initialize linear operator and vectors.
        A = linop_cdp() ; call init_rand(A)
        allocate(b, c, cref) ; call init_rand(b) ; call c%zero()

        ! Reference: dense computation with the adjoint operator.
        cref%data = matmul(expm(tau*hermitian(A%data)), b%data)

        ! kexpm.
        call kexpm(c, A, b, tau, rtol_dp, info, trans=.true.)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_adjoint_cdp')
        call c%sub(cref) ; err = c%norm() / cref%norm()
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_kexptA_adjoint_cdp', &
                              & info='kexpm', eq='c = expm(tau A.H) @ b', context=msg)

        ! krylov_exptA.
        call krylov_exptA(c, A, b, tau, info, trans=.true.)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_adjoint_cdp')
        call c%sub(cref) ; err = c%norm() / cref%norm()
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_kexptA_adjoint_cdp', &
                              & info='krylov_exptA', eq='c = expm(tau A.H) @ b', context=msg)

        return
    end subroutine test_kexptA_adjoint_cdp

    subroutine test_kexptA_zero_input_cdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test linear operator.
        type(linop_cdp), allocatable :: A
        ! Test vectors.
        integer, parameter :: p = 2
        type(vector_cdp), allocatable :: b, c, Bblk(:), Cblk(:)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(dp), parameter :: tau = 0.1_dp
        real(dp) :: err
        integer :: i
        character(len=256) :: msg

        ! Initialize linear operator, zero inputs and non-zero outputs.
        A = linop_cdp() ; call init_rand(A)
        allocate(b, c) ; call b%zero()
        allocate(Bblk(p), Cblk(p)) ; call zero_basis(Bblk)

        do
            ! kexpm.
            call c%rand()
            call kexpm(c, A, b, tau, rtol_dp, info)
            err = c%norm()
            call check(error, info > 0 .and. err == 0.0_dp)
            if (allocated(error)) exit

            ! krylov_exptA.
            call c%rand()
            call krylov_exptA(c, A, b, tau, info)
            err = c%norm()
            call check(error, info > 0 .and. err == 0.0_dp)
            if (allocated(error)) exit

            ! Block kexpm.
            call rand_basis(Cblk)
            call kexpm(Cblk, A, Bblk, tau, rtol_dp, info)
            err = 0.0_dp
            do i = 1, p
                err = err + Cblk(i)%norm()
            enddo
            call check(error, info > 0 .and. err == 0.0_dp)
            exit
        enddo
        call get_err_str(msg, "norm: ", err)
        call check_test(error, 'test_kexptA_zero_input_cdp', &
                              & info='Zero input', eq='expm(tau A) @ 0 = 0', context=msg)

        return
    end subroutine test_kexptA_zero_input_cdp

    subroutine test_kexptA_breakdown_cdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test linear operator.
        type(linop_cdp), allocatable :: A
        ! Test vectors.
        type(vector_cdp), allocatable :: b, c, cref
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        integer, parameter :: k_inv = 10
        real(dp), parameter :: tau = 0.1_dp
        real(dp) :: err
        character(len=256) :: msg

        ! Create block lower triangular matrix: span(e_1, ..., e_k_inv) is invariant.
        A = linop_cdp() ; call init_rand(A)
        A%data(k_inv+1:, :k_inv) = zero_cdp

        ! Input vector confined to the invariant subspace.
        allocate(b, c, cref) ; call init_rand(b) ; call c%zero()
        b%data(k_inv+1:) = zero_cdp

        ! Reference: dense computation.
        cref%data = matmul(expm(tau*A%data), b%data)

        ! kexpm with tol = 0: convergence is only possible through the Arnoldi breakdown,
        ! which must happen after exactly k_inv steps (the approximation is then exact).
        call kexpm(c, A, b, tau, 0.0_dp, info)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_breakdown_cdp')
        call check(error, info == k_inv)
        call check_test(error, 'test_kexptA_breakdown_cdp', &
                              & info='Subspace Dim', eq='info == k_inv', context='Invariant detection')

        call c%sub(cref) ; err = c%norm() / cref%norm()
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_kexptA_breakdown_cdp', &
                              & info='Invariant Property', eq='c = expm(tau A) @ b', context=msg)

        ! krylov_exptA: a single sub-step, never more than k_inv matrix-vector products.
        call krylov_exptA(c, A, b, tau, info)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_breakdown_cdp')
        call c%sub(cref) ; err = c%norm() / cref%norm()
        call get_err_str(msg, "max err: ", err)
        call check(error, info <= k_inv .and. err < rtol_dp)
        call check_test(error, 'test_kexptA_breakdown_cdp', &
                              & info='krylov_exptA', eq='<= k_inv matvecs, c = expm(tau A) @ b', context=msg)

        return
    end subroutine test_kexptA_breakdown_cdp

    subroutine test_kexptA_stiff_cdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test linear operator.
        type(linop_cdp), allocatable :: A
        ! Test vectors.
        type(vector_cdp), allocatable :: b, c, cref
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(dp), parameter :: tau = 1.0_dp, s = 30.0_dp
        real(dp) :: err
        character(len=256) :: msg

        ! Stable random operator with || tau A || ~ 60: far beyond what a single
        ! 30-vector Krylov basis can resolve.
        A = linop_cdp() ; call init_rand(A)
        A%data = s*A%data/sqrt(real(test_size, dp)) - s*eye(test_size, mold=1.0_dp)
        allocate(b, c, cref) ; call init_rand(b) ; call c%zero()

        ! Reference: dense computation.
        cref%data = matmul(expm(tau*A%data), b%data)

        ! krylov_exptA must converge thanks to the adaptive sub-stepping (> 30 matvecs).
        call krylov_exptA(c, A, b, tau, info)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_stiff_cdp')
        call c%sub(cref) ; err = c%norm() / cref%norm()
        call get_err_str(msg, "max err: ", err)
        call check(error, info > 30 .and. err < rtol_dp)
        call check_test(error, 'test_kexptA_stiff_cdp', &
                              & info='Sub-stepping', eq='c = expm(tau A) @ b', context=msg)

        return
    end subroutine test_kexptA_stiff_cdp

    subroutine test_kexptA_warm_start_cdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test linear operator.
        type(linop_cdp), allocatable :: A
        ! Test vectors.
        type(vector_cdp), allocatable :: b, c, cref
        ! Sub-stepping options and metadata.
        type(kexpm_dp_opts) :: opts
        type(kexpm_dp_metadata) :: meta
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        integer, parameter :: ncalls = 10
        real(dp), parameter :: tau = 1.0_dp, s = 30.0_dp
        complex(dp), allocatable :: E(:, :)
        real(dp) :: err, tol
        integer :: icall, cold, last
        character(len=256) :: msg

        ! Stiff stable random operator.
        A = linop_cdp() ; call init_rand(A)
        A%data = s*A%data/sqrt(real(test_size, dp)) - s*eye(test_size, mold=1.0_dp)
        E = expm(tau*A%data)
        allocate(b, c, cref) ; call c%zero()

        ! Repeated calls with a persistent metadata object and new right-hand sides.
        err = 0.0_dp
        do icall = 1, ncalls
            call init_rand(b) ; cref%data = matmul(E, b%data) ; tol = atol_dp*b%norm()
            call kexpm_substep(c, A, b, tau, tol, info, meta=meta)
            call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_warm_start_cdp')
            if (icall == 1) cold = info
            last = info
            call c%sub(cref) ; err = max(err, c%norm() / cref%norm())
        enddo
        write(msg, '(A,E9.2,A,I0,A,I0,A,I0)') 'max err: ', err, ', cold: ', cold, ', last: ', last, ', nsub: ', meta%nsub
        call check(error, err < rtol_dp)
        call check_test(error, 'test_kexptA_warm_start_cdp', &
                              & info='Accuracy', eq='c = expm(tau A) @ b', context=msg)

        ! Settings learned, frozen, and not more expensive than the cold call (up to the
        ! call-to-call variation of the adaptive cost with the right-hand side).
        call check(error, meta%learned .and. meta%frozen .and. meta%n_calls == ncalls &
                        & .and. last <= cold + max(2, cold/10))
        call check_test(error, 'test_kexptA_warm_start_cdp', &
                              & info='Learned settings', eq='frozen, cost <= cold cost', context=msg)

        ! reset() keeps the learned settings, forget() clears them.
        do
            call meta%reset()
            call check(error, meta%learned .and. meta%n_calls == 0 .and. meta%n_matvec == 0)
            if (allocated(error)) exit
            call meta%forget()
            call check(error, .not. meta%learned .and. meta%nsub == 0)
            exit
        enddo
        call check_test(error, 'test_kexptA_warm_start_cdp', &
                              & info='reset / forget', eq='reset keeps, forget clears', context='')

        ! tune = .false. : the sub-step count of the cold call is kept.
        opts = kexpm_dp_opts(tune=.false.)
        call init_rand(b) ; cref%data = matmul(E, b%data) ; tol = atol_dp*b%norm()
        call kexpm_substep(c, A, b, tau, tol, info, options=opts, meta=meta)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_warm_start_cdp')
        call kexpm_substep(c, A, b, tau, tol, info, options=opts, meta=meta)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_warm_start_cdp')
        call c%sub(cref) ; err = c%norm() / cref%norm()
        write(msg, '(A,E9.2,A,I0,A,I0)') 'max err: ', err, ', nsub: ', meta%nsub, ', nsub0: ', meta%nsub0
        call check(error, err < rtol_dp .and. meta%frozen .and. meta%nsub == meta%nsub0)
        call check_test(error, 'test_kexptA_warm_start_cdp', &
                              & info='No tuning', eq='nsub == nsub0', context=msg)

        ! A different horizon invalidates the learned settings.
        cref%data = matmul(expm(0.5_dp*tau*A%data), b%data)
        call kexpm_substep(c, A, b, 0.5_dp*tau, tol, info, meta=meta)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexptA_warm_start_cdp')
        call c%sub(cref) ; err = c%norm() / cref%norm()
        write(msg, '(A,E9.2,A,E9.2)') 'max err: ', err, ', tau_learned: ', meta%tau_learned
        call check(error, err < rtol_dp .and. abs(meta%tau_learned - 0.5_dp*tau) <= epsilon(tau))
        call check_test(error, 'test_kexptA_warm_start_cdp', &
                              & info='New horizon', eq='settings re-learned', context=msg)

        return
    end subroutine test_kexptA_warm_start_cdp

    subroutine test_block_kexptA_cdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test linear operator.
        type(linop_cdp), allocatable :: A
        ! Test vectors.
        integer, parameter :: p = 3
        type(vector_cdp), allocatable :: B(:), C(:)
        ! Krylov subspace (user-supplied workspace).
        integer, parameter :: kdim = 15
        type(vector_cdp), allocatable :: X(:)
        ! Information flag.
        integer :: info
        ! Miscellaneous.
        real(dp), parameter :: tau = 0.1_dp
        complex(dp), allocatable :: Bdata(:, :), Cdata(:, :), Cref(:, :)
        real(dp) :: err
        integer :: i
        character(len=256) :: msg

        ! Initialize linear operator and input block.
        A = linop_cdp() ; call init_rand(A)
        allocate(B(p), C(p), X(p*(kdim+1))) ; call init_rand(B) ; call zero_basis(C) ; call zero_basis(X)
        allocate(Bdata(test_size, p), Cdata(test_size, p)) ; call get_data(Bdata, B)

        ! Reference: dense computation.
        Cref = matmul(expm(tau*A%data), Bdata)

        ! Column-by-column kexpm.
        do i = 1, p
            call kexpm(C(i), A, B(i), tau, rtol_dp, info, kdim=kdim)
            call check_info(info, 'kexpm', module=this_module_long, procedure='test_block_kexptA_cdp')
        enddo
        call get_data(Cdata, C)
        err = norm2(abs(Cdata - Cref)) / norm2(abs(Cref))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_block_kexptA_cdp', &
                              & info='Column-wise', eq='C = expm(tau A) @ B', context=msg)

        ! Block kexpm.
        call kexpm(C, A, B, tau, rtol_dp, info, kdim=kdim)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_block_kexptA_cdp')
        call get_data(Cdata, C)
        err = norm2(abs(Cdata - Cref)) / norm2(abs(Cref))
        call get_err_str(msg, "max err: ", err)
        call check(error, info > 0 .and. info <= p*(kdim+1) .and. err < rtol_dp)
        call check_test(error, 'test_block_kexptA_cdp', &
                              & info='Block', eq='C = expm(tau A) @ B', context=msg)

        ! Block kexpm with user-supplied workspace.
        call kexpm(C, A, B, tau, rtol_dp, info, X=X)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_block_kexptA_cdp')
        call get_data(Cdata, C)
        err = norm2(abs(Cdata - Cref)) / norm2(abs(Cref))
        call get_err_str(msg, "max err: ", err)
        call check(error, info > 0 .and. err < rtol_dp)
        call check_test(error, 'test_block_kexptA_cdp', &
                              & info='Block, workspace', eq='C = expm(tau A) @ B', context=msg)

        block
        integer, parameter :: pinv = 2, k_inv = 2*pinv
        type(vector_cdp), allocatable :: B2(:), C2(:)
        complex(dp), allocatable :: B2data(:, :), C2data(:, :), C2ref(:, :)

        ! Create block lower triangular matrix: span(e_1, ..., e_k_inv) is invariant.
        call init_rand(A)
        A%data(k_inv+1:, :k_inv) = zero_cdp

        ! Input block confined to the invariant subspace.
        allocate(B2(pinv), C2(pinv)) ; call init_rand(B2) ; call zero_basis(C2)
        do i = 1, pinv
            B2(i)%data(k_inv+1:) = zero_cdp
        enddo
        allocate(B2data(test_size, pinv), C2data(test_size, pinv)) ; call get_data(B2data, B2)
        C2ref = matmul(expm(tau*A%data), B2data)

        ! tol = 0: block Arnoldi must break down after k_inv/pinv block steps (exact approximation).
        call kexpm(C2, A, B2, tau, 0.0_dp, info)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_block_kexptA_cdp')
        call check(error, info == k_inv)
        call check_test(error, 'test_block_kexptA_cdp', &
                              & info='Subspace Dim', eq='info == k_inv', context='Invariant detection')

        call get_data(C2data, C2)
        err = norm2(abs(C2data - C2ref)) / norm2(abs(C2ref))
        call get_err_str(msg, "max err: ", err)
        call check(error, err < rtol_dp)
        call check_test(error, 'test_block_kexptA_cdp', &
                              & info='Invariant Property', eq='C = expm(tau A) @ B', context=msg)
        end block

        return
    end subroutine test_block_kexptA_cdp

    !-----------------------------------------------------------
    !-----     UNIT TESTS FOR DENSE MATRIX SQUARE ROOT     -----
    !-----------------------------------------------------------

    subroutine collect_sqrtm_rsp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                        new_unittest("Dense sqrtm for positive definite matrices.", test_dense_sqrtm_pos_def_rsp), &
                        new_unittest("Dense sqrtm for positive semi-definite matrices.", test_dense_sqrtm_pos_semi_def_rsp) &
                    ]

        return
    end subroutine collect_sqrtm_rsp_testsuite

    subroutine test_dense_sqrtm_pos_def_rsp(error)
       !> This function tests the matrix version of the sqrt function for the case of
       ! a symmetric/hermitian positive definite matrix
    
       !> Error type to be returned.
       type(error_type), allocatable, intent(out) :: error
       !> Problem dimension.
       integer, parameter :: n = 5
       !> Test matrix.
       real(sp) :: A(n, n)
       real(sp) :: sqrtmA(n, n)
       complex(sp) :: lambda(n)
       real(sp) :: err
       integer :: i, info
       character(len=256) :: msg
    
       ! --> Initialize matrix.
       call random_number(A)
       ! make symmetric/hermitian positive definite
       A = 0.5_sp*(A + transpose(A))
       call eig(A, sqrtmA, lambda)
       do i = 1,n
          lambda(i) = abs(lambda(i)) + 0.1_sp
       end do
       ! reconstruct matrix
       A = matmul(sqrtmA, matmul(diag(abs(lambda)), transpose(sqrtmA)))
     
       ! compute matrix square root
       call sqrtm(A, sqrtmA, info)
       call check_info(info, 'sqrtm', module=this_module_long, procedure='test_dense_sqrtm_pos_def_rsp')
    
       err = maxval(abs(matmul(sqrtmA, sqrtmA) - A))
       call get_err_str(msg, "max err: ", err)
       call check(error, err < rtol_sp)
       call check_test(error, 'test_dense_sqrtm_pos_def_rsp', eq='sqrt(A)**2 = A', context=msg)
    
       return
    end subroutine test_dense_sqrtm_pos_def_rsp
    
    subroutine test_dense_sqrtm_pos_semi_def_rsp(error)
       !> This function tests the matrix version of the sqrt function for the case 
       ! of a symmetric semi-definite matrix
    
       !> Error type to be returned.
       type(error_type), allocatable, intent(out) :: error
       !> Problem dimension.
       integer, parameter :: n = 25
       !> Test matrix.
       real(sp) :: A(n, n)
       real(sp) :: sqrtmA(n, n)
       complex(sp) :: lambda(n)
       real(sp) :: err
       integer :: i, info
       character(len=256) :: msg
    
       ! --> Initialize matrix.
       call random_number(A)
       ! make symmetric/hermitian positive semi-definite
       A = 0.5_sp*(A + transpose(A))
       call eig(A, sqrtmA, lambda)
       do i = 1,n-1
          lambda(i) = abs(lambda(i)) + 0.1_sp
       end do
       lambda(n) = zero_rsp
       ! reconstruct matrix
       A = matmul(sqrtmA, matmul(diag(abs(lambda)), transpose(sqrtmA)))
    
       ! compute matrix square root
       call sqrtm(A, sqrtmA, info)
       call check_info(info, 'sqrtm', module=this_module_long, procedure='test_dense_sqrtm_pos_semi_def_rsp')
    
       err = maxval(abs(matmul(sqrtmA, sqrtmA) - A))
       call get_err_str(msg, "max err: ", err)
       call check(error, err < rtol_sp)
       call check_test(error, 'test_dense_sqrtm_pos_semi_def_rsp', eq='sqrt(A)**2 = A', context=trim(msg))
    
       return
    end subroutine test_dense_sqrtm_pos_semi_def_rsp

    subroutine collect_sqrtm_rdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                        new_unittest("Dense sqrtm for positive definite matrices.", test_dense_sqrtm_pos_def_rdp), &
                        new_unittest("Dense sqrtm for positive semi-definite matrices.", test_dense_sqrtm_pos_semi_def_rdp) &
                    ]

        return
    end subroutine collect_sqrtm_rdp_testsuite

    subroutine test_dense_sqrtm_pos_def_rdp(error)
       !> This function tests the matrix version of the sqrt function for the case of
       ! a symmetric/hermitian positive definite matrix
    
       !> Error type to be returned.
       type(error_type), allocatable, intent(out) :: error
       !> Problem dimension.
       integer, parameter :: n = 5
       !> Test matrix.
       real(dp) :: A(n, n)
       real(dp) :: sqrtmA(n, n)
       complex(dp) :: lambda(n)
       real(dp) :: err
       integer :: i, info
       character(len=256) :: msg
    
       ! --> Initialize matrix.
       call random_number(A)
       ! make symmetric/hermitian positive definite
       A = 0.5_dp*(A + transpose(A))
       call eig(A, sqrtmA, lambda)
       do i = 1,n
          lambda(i) = abs(lambda(i)) + 0.1_dp
       end do
       ! reconstruct matrix
       A = matmul(sqrtmA, matmul(diag(abs(lambda)), transpose(sqrtmA)))
     
       ! compute matrix square root
       call sqrtm(A, sqrtmA, info)
       call check_info(info, 'sqrtm', module=this_module_long, procedure='test_dense_sqrtm_pos_def_rdp')
    
       err = maxval(abs(matmul(sqrtmA, sqrtmA) - A))
       call get_err_str(msg, "max err: ", err)
       call check(error, err < rtol_dp)
       call check_test(error, 'test_dense_sqrtm_pos_def_rdp', eq='sqrt(A)**2 = A', context=msg)
    
       return
    end subroutine test_dense_sqrtm_pos_def_rdp
    
    subroutine test_dense_sqrtm_pos_semi_def_rdp(error)
       !> This function tests the matrix version of the sqrt function for the case 
       ! of a symmetric semi-definite matrix
    
       !> Error type to be returned.
       type(error_type), allocatable, intent(out) :: error
       !> Problem dimension.
       integer, parameter :: n = 25
       !> Test matrix.
       real(dp) :: A(n, n)
       real(dp) :: sqrtmA(n, n)
       complex(dp) :: lambda(n)
       real(dp) :: err
       integer :: i, info
       character(len=256) :: msg
    
       ! --> Initialize matrix.
       call random_number(A)
       ! make symmetric/hermitian positive semi-definite
       A = 0.5_dp*(A + transpose(A))
       call eig(A, sqrtmA, lambda)
       do i = 1,n-1
          lambda(i) = abs(lambda(i)) + 0.1_dp
       end do
       lambda(n) = zero_rdp
       ! reconstruct matrix
       A = matmul(sqrtmA, matmul(diag(abs(lambda)), transpose(sqrtmA)))
    
       ! compute matrix square root
       call sqrtm(A, sqrtmA, info)
       call check_info(info, 'sqrtm', module=this_module_long, procedure='test_dense_sqrtm_pos_semi_def_rdp')
    
       err = maxval(abs(matmul(sqrtmA, sqrtmA) - A))
       call get_err_str(msg, "max err: ", err)
       call check(error, err < rtol_dp)
       call check_test(error, 'test_dense_sqrtm_pos_semi_def_rdp', eq='sqrt(A)**2 = A', context=trim(msg))
    
       return
    end subroutine test_dense_sqrtm_pos_semi_def_rdp

    subroutine collect_sqrtm_csp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                        new_unittest("Dense sqrtm for positive definite matrices.", test_dense_sqrtm_pos_def_csp), &
                        new_unittest("Dense sqrtm for positive semi-definite matrices.", test_dense_sqrtm_pos_semi_def_csp) &
                    ]

        return
    end subroutine collect_sqrtm_csp_testsuite

    subroutine test_dense_sqrtm_pos_def_csp(error)
       !> This function tests the matrix version of the sqrt function for the case of
       ! a symmetric/hermitian positive definite matrix
    
       !> Error type to be returned.
       type(error_type), allocatable, intent(out) :: error
       !> Problem dimension.
       integer, parameter :: n = 5
       !> Test matrix.
       complex(sp) :: A(n, n)
       complex(sp) :: sqrtmA(n, n)
       complex(sp) :: lambda(n)
       real(sp) :: err
       integer :: i, info
       character(len=256) :: msg
    
       ! --> Initialize matrix.
       call random_number(A%re)
       call random_number(A%im)
       ! make symmetric/hermitian positive definite
       A = 0.5_sp*(A + conjg(transpose(A)))
       call eig(A, sqrtmA, lambda)
       do i = 1,n
          lambda(i) = abs(lambda(i)) + 0.1_sp
       end do
       ! reconstruct matrix
       A = matmul(sqrtmA, matmul(diag(abs(lambda)), conjg(transpose(sqrtmA))))
     
       ! compute matrix square root
       call sqrtm(A, sqrtmA, info)
       call check_info(info, 'sqrtm', module=this_module_long, procedure='test_dense_sqrtm_pos_def_csp')
    
       err =  maxval(abs(matmul(sqrtmA, sqrtmA) - A))
       call get_err_str(msg, "max err: ", err)
       call check(error, err < rtol_sp)
       call check_test(error, 'test_dense_sqrtm_pos_def_csp', eq='sqrt(A)**2 = A', context=msg)
    
       return
    end subroutine test_dense_sqrtm_pos_def_csp
    
    subroutine test_dense_sqrtm_pos_semi_def_csp(error)
       !> This function tests the matrix version of the sqrt function for the case 
       ! of a symmetric semi-definite matrix
    
       !> Error type to be returned.
       type(error_type), allocatable, intent(out) :: error
       !> Problem dimension.
       integer, parameter :: n = 25
       !> Test matrix.
       complex(sp) :: A(n, n)
       complex(sp) :: sqrtmA(n, n)
       complex(sp) :: lambda(n)
       real(sp) :: err
       integer :: i, info
       character(len=256) :: msg
    
       ! --> Initialize matrix.
       call random_number(A%re)
       call random_number(A%im)
       ! make symmetric/hermitian positive semi-definite
       A = 0.5_sp*(A + conjg(transpose(A)))
       call eig(A, sqrtmA, lambda)
       do i = 1,n-1
          lambda(i) = abs(lambda(i)) + 0.1_sp
       end do
       lambda(n) = zero_rsp
       ! reconstruct matrix
       A = matmul(sqrtmA, matmul(diag(abs(lambda)), conjg(transpose(sqrtmA))))
    
       ! compute matrix square root
       call sqrtm(A, sqrtmA, info)
       call check_info(info, 'sqrtm', module=this_module_long, procedure='test_dense_sqrtm_pos_semi_def_csp')
    
       err = maxval(abs(matmul(sqrtmA, sqrtmA) - A))
       call get_err_str(msg, "max err: ", err)
       call check(error, err < rtol_sp)
       call check_test(error, 'test_dense_sqrtm_pos_semi_def_csp', eq='sqrt(A)**2 = A', context=trim(msg))
    
       return
    end subroutine test_dense_sqrtm_pos_semi_def_csp

    subroutine collect_sqrtm_cdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                        new_unittest("Dense sqrtm for positive definite matrices.", test_dense_sqrtm_pos_def_cdp), &
                        new_unittest("Dense sqrtm for positive semi-definite matrices.", test_dense_sqrtm_pos_semi_def_cdp) &
                    ]

        return
    end subroutine collect_sqrtm_cdp_testsuite

    subroutine test_dense_sqrtm_pos_def_cdp(error)
       !> This function tests the matrix version of the sqrt function for the case of
       ! a symmetric/hermitian positive definite matrix
    
       !> Error type to be returned.
       type(error_type), allocatable, intent(out) :: error
       !> Problem dimension.
       integer, parameter :: n = 5
       !> Test matrix.
       complex(dp) :: A(n, n)
       complex(dp) :: sqrtmA(n, n)
       complex(dp) :: lambda(n)
       real(dp) :: err
       integer :: i, info
       character(len=256) :: msg
    
       ! --> Initialize matrix.
       call random_number(A%re)
       call random_number(A%im)
       ! make symmetric/hermitian positive definite
       A = 0.5_dp*(A + conjg(transpose(A)))
       call eig(A, sqrtmA, lambda)
       do i = 1,n
          lambda(i) = abs(lambda(i)) + 0.1_dp
       end do
       ! reconstruct matrix
       A = matmul(sqrtmA, matmul(diag(abs(lambda)), conjg(transpose(sqrtmA))))
     
       ! compute matrix square root
       call sqrtm(A, sqrtmA, info)
       call check_info(info, 'sqrtm', module=this_module_long, procedure='test_dense_sqrtm_pos_def_cdp')
    
       err =  maxval(abs(matmul(sqrtmA, sqrtmA) - A))
       call get_err_str(msg, "max err: ", err)
       call check(error, err < rtol_dp)
       call check_test(error, 'test_dense_sqrtm_pos_def_cdp', eq='sqrt(A)**2 = A', context=msg)
    
       return
    end subroutine test_dense_sqrtm_pos_def_cdp
    
    subroutine test_dense_sqrtm_pos_semi_def_cdp(error)
       !> This function tests the matrix version of the sqrt function for the case 
       ! of a symmetric semi-definite matrix
    
       !> Error type to be returned.
       type(error_type), allocatable, intent(out) :: error
       !> Problem dimension.
       integer, parameter :: n = 25
       !> Test matrix.
       complex(dp) :: A(n, n)
       complex(dp) :: sqrtmA(n, n)
       complex(dp) :: lambda(n)
       real(dp) :: err
       integer :: i, info
       character(len=256) :: msg
    
       ! --> Initialize matrix.
       call random_number(A%re)
       call random_number(A%im)
       ! make symmetric/hermitian positive semi-definite
       A = 0.5_dp*(A + conjg(transpose(A)))
       call eig(A, sqrtmA, lambda)
       do i = 1,n-1
          lambda(i) = abs(lambda(i)) + 0.1_dp
       end do
       lambda(n) = zero_rdp
       ! reconstruct matrix
       A = matmul(sqrtmA, matmul(diag(abs(lambda)), conjg(transpose(sqrtmA))))
    
       ! compute matrix square root
       call sqrtm(A, sqrtmA, info)
       call check_info(info, 'sqrtm', module=this_module_long, procedure='test_dense_sqrtm_pos_semi_def_cdp')
    
       err = maxval(abs(matmul(sqrtmA, sqrtmA) - A))
       call get_err_str(msg, "max err: ", err)
       call check(error, err < rtol_dp)
       call check_test(error, 'test_dense_sqrtm_pos_semi_def_cdp', eq='sqrt(A)**2 = A', context=trim(msg))
    
       return
    end subroutine test_dense_sqrtm_pos_semi_def_cdp


end module TestExpmLib