module TestWorkspace
    !! Tests for the caller-provided workspace (`work` arguments) of the solvers.
    !!
    !! Each test first runs a routine letting it allocate its own workspace (reference), then
    !! runs it twice more on identical inputs with the same caller-provided `work` array:
    !!
    !! - trial 1/2, fresh workspace: `work` is freshly allocated and unshaped, so the routine
    !!   must shape it itself (`prepare_workspace` -> `init_like`) and slice it correctly.
    !! - trial 2/2, reused workspace: `work` is already shaped and still holds stale data from
    !!   trial 1. This checks that `init_like` takes its no-op path, that the routine does not
    !!   depend on the previous content of the workspace, and that the first call left it usable
    !!   (in particular, did not `free` it).
    !!
    !! Both trials must reproduce the reference result. The log line of each check names the
    !! trial (`info` column) and the quantity compared (`eq` column).
    !!
    !! These tests do not detect hidden internal allocations (the test vector types manage
    !! their own memory); that is checked separately with a resource-counting vector type.
    ! Testdrive
    use testdrive, only: new_unittest, unittest_type, error_type, check
    ! LightKrylov
    use LightKrylov
    use LightKrylov_Constants
    use LightKrylov_Logger
    use LightKrylov_AbstractVectors
    ! Test Utilities
    use LightKrylov_TestUtils
    use TestUtils

    implicit none (type, external)

    private

    character(len=*), parameter, private :: this_module      = 'LK_TWork'
    character(len=*), parameter, private :: this_module_long = 'LightKrylov_TestWorkspace'

    public :: collect_workspace_rsp_testsuite
    public :: collect_workspace_rdp_testsuite
    public :: collect_workspace_csp_testsuite
    public :: collect_workspace_cdp_testsuite

contains

    pure function trial_info(trial) result(info)
        !! Label of a workspace trial in the test log (at most 30 characters).
        integer, intent(in) :: trial
        character(len=30) :: info
        if (trial == 1) then
            info = 'fresh (unshaped) work'
        else
            info = 'reused (stale) work  '
        endif
    end function trial_info

    pure function newton_trial_info(trial) result(info)
        !! Label of a Newton workspace trial in the test log (at most 30 characters).
        integer, intent(in) :: trial
        character(len=30) :: info
        if (trial == 1) then
            info = 'fresh work, no bisection'
        else
            info = 'reused work + bisection'
        endif
    end function newton_trial_info

    subroutine collect_workspace_rsp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [ &
            new_unittest("GMRES with workspace", test_gmres_work_rsp), &
            new_unittest("FGMRES with workspace", test_fgmres_work_rsp), &
            new_unittest("CG with workspace", test_cg_work_rsp), &
            new_unittest("KS eigs with workspace", test_eigs_work_rsp), &
            new_unittest("Sym. eighs with workspace", test_eighs_work_rsp), &
            new_unittest("SVDS with workspace", test_svds_work_rsp), &
            new_unittest("Krylov expm with workspace", test_kexpm_work_rsp), &
            new_unittest("Newton-Krylov with workspace", test_newton_work_rsp) &
        ]
    end subroutine collect_workspace_rsp_testsuite

    !-----------------------------------------------------
    !-----     LINEAR SOLVERS WITH A WORKSPACE         -----
    !-----------------------------------------------------

    subroutine test_gmres_work_rsp(error)
        type(error_type), allocatable, intent(out) :: error
        type(linop_rsp), allocatable :: A
        type(vector_rsp), allocatable :: b, xref, x
        type(vector_rsp), allocatable :: work(:)
        type(gmres_sp_opts) :: opts
        integer :: info, trial
        real(sp) :: err
        character(len=256) :: msg

        A = linop_rsp() ; call init_rand(A)
        b = vector_rsp() ; call init_rand(b)
        xref = vector_rsp() ; call xref%zero()
        opts = gmres_sp_opts(kdim=test_size)

        ! Reference: workspace allocated by the solver.
        call gmres(A, b, xref, info, rtol=rtol_sp, atol=atol_sp, options=opts)
        call check_info(info, 'gmres', module=this_module_long, procedure='test_gmres_work_rsp')

        ! Caller-provided workspace, used for two consecutive solves.
        allocate(work(gmres_work_size(test_size)))
        do trial = 1, 2
            x = vector_rsp() ; call x%zero()
            call gmres(A, b, x, info, rtol=rtol_sp, atol=atol_sp, options=opts, work=work)
            call check_info(info, 'gmres', module=this_module_long, procedure='test_gmres_work_rsp')
            call x%sub(xref) ; err = x%norm() / xref%norm()
            call get_err_str(msg, "max err: ", err)
            call check(error, err < rtol_sp)
            call check_test(error, 'test_gmres_work_rsp', &
                          & info=trial_info(trial), eq='x(work) == x(internal)', context=msg)
        enddo
    end subroutine test_gmres_work_rsp

    subroutine test_fgmres_work_rsp(error)
        type(error_type), allocatable, intent(out) :: error
        type(linop_rsp), allocatable :: A
        type(vector_rsp), allocatable :: b, xref, x
        type(vector_rsp), allocatable :: work(:)
        type(fgmres_sp_opts) :: opts
        integer :: info, trial
        real(sp) :: err
        character(len=256) :: msg

        A = linop_rsp() ; call init_rand(A)
        b = vector_rsp() ; call init_rand(b)
        xref = vector_rsp() ; call xref%zero()
        opts = fgmres_sp_opts(kdim=test_size)

        call fgmres(A, b, xref, info, rtol=rtol_sp, atol=atol_sp, options=opts)
        call check_info(info, 'fgmres', module=this_module_long, procedure='test_fgmres_work_rsp')

        allocate(work(fgmres_work_size(test_size)))
        do trial = 1, 2
            x = vector_rsp() ; call x%zero()
            call fgmres(A, b, x, info, rtol=rtol_sp, atol=atol_sp, options=opts, work=work)
            call check_info(info, 'fgmres', module=this_module_long, procedure='test_fgmres_work_rsp')
            call x%sub(xref) ; err = x%norm() / xref%norm()
            call get_err_str(msg, "max err: ", err)
            call check(error, err < rtol_sp)
            call check_test(error, 'test_fgmres_work_rsp', &
                          & info=trial_info(trial), eq='x(work) == x(internal)', context=msg)
        enddo
    end subroutine test_fgmres_work_rsp

    subroutine test_cg_work_rsp(error)
        type(error_type), allocatable, intent(out) :: error
        type(spd_linop_rsp), allocatable :: A
        type(vector_rsp), allocatable :: b, xref, x
        type(vector_rsp), allocatable :: work(:)
        type(cg_sp_opts) :: opts
        integer :: info, trial
        real(sp) :: err
        character(len=256) :: msg

        A = spd_linop_rsp() ; call init_rand(A)
        b = vector_rsp() ; call init_rand(b)
        xref = vector_rsp() ; call xref%zero()
        opts = cg_sp_opts(maxiter=2*test_size)

        call cg(A, b, xref, info, rtol=rtol_sp, atol=atol_sp, options=opts)
        call check_info(info, 'cg', module=this_module_long, procedure='test_cg_work_rsp')

        allocate(work(cg_work_size()))
        do trial = 1, 2
            x = vector_rsp() ; call x%zero()
            call cg(A, b, x, info, rtol=rtol_sp, atol=atol_sp, options=opts, work=work)
            call check_info(info, 'cg', module=this_module_long, procedure='test_cg_work_rsp')
            call x%sub(xref) ; err = x%norm() / xref%norm()
            call get_err_str(msg, "max err: ", err)
            call check(error, err < rtol_sp)
            call check_test(error, 'test_cg_work_rsp', &
                          & info=trial_info(trial), eq='x(work) == x(internal)', context=msg)
        enddo
    end subroutine test_cg_work_rsp

    !--------------------------------------------------------
    !-----     EIGENVALUE / SVD SOLVERS WITH A WORKSPACE     -----
    !--------------------------------------------------------

    subroutine test_eigs_work_rsp(error)
        !! Same Toeplitz problem as the Krylov-Schur eigs test (nev < n, so restarts happen
        !! and the Krylov-Schur part of the workspace is used).
        type(error_type), allocatable, intent(out) :: error
        type(linop_rsp), allocatable :: A
        type(vector_rsp), allocatable :: X(:), x0
        type(vector_rsp), allocatable :: work(:)
        complex(sp), allocatable :: eigvals_ref(:), eigvals(:)
        real(sp), allocatable :: residuals(:)
        integer, parameter :: nev = 8
        integer :: i, n, info, trial
        real(sp) :: a_, b_
        real(sp) :: err
        character(len=256) :: msg

        A = linop_rsp() ; A%data = 0.0_sp ; n = size(A%data, 1)
        call random_number(a_) ; call random_number(b_) ; b_ = abs(b_)
        do i = 1, n
            A%data(i, i) = a_
            if (i < n) then
                A%data(i, i+1) = b_
                A%data(i+1, i) = -b_
            endif
        enddo
        x0 = vector_rsp() ; call init_rand(x0)

        allocate(X(nev)) ; call zero_basis(X)
        call eigs(A, X, eigvals_ref, residuals, info, x0=x0, tolerance=atol_sp)
        call check_info(info, 'eigs', module=this_module_long, procedure='test_eigs_work_rsp')

        allocate(work(eigs_work_size(nev)))
        do trial = 1, 2
            call zero_basis(X)
            call eigs(A, X, eigvals, residuals, info, x0=x0, tolerance=atol_sp, work=work)
            call check_info(info, 'eigs', module=this_module_long, procedure='test_eigs_work_rsp')
            err = maxval(abs(eigvals - eigvals_ref))
            call get_err_str(msg, "max err: ", err)
            call check(error, err < rtol_sp)
            call check_test(error, 'test_eigs_work_rsp', &
                          & info=trial_info(trial), eq='evals(work) == evals(internal)', context=msg)
        enddo
    end subroutine test_eigs_work_rsp

    subroutine test_eighs_work_rsp(error)
        type(error_type), allocatable, intent(out) :: error
        type(spd_linop_rsp), allocatable :: A
        type(vector_rsp), allocatable :: X(:), x0
        type(vector_rsp), allocatable :: work(:)
        real(sp), allocatable :: eigvals_ref(:), eigvals(:), residuals(:)
        integer :: i, info, trial
        real(sp) :: err
        character(len=256) :: msg

        ! Strang matrix (symmetric positive definite, real-valued).
        A = spd_linop_rsp()
        A%data = zero_rsp
        do i = 1, test_size
            A%data(i, i) = 2*one_rsp
            if (i < test_size) then
                A%data(i, i+1) = -one_rsp
                A%data(i+1, i) = -one_rsp
            endif
        enddo
        x0 = vector_rsp() ; call init_rand(x0)

        allocate(X(test_size)) ; call zero_basis(X)
        call eighs(A, X, eigvals_ref, residuals, info, x0=x0, tolerance=atol_sp)
        call check_info(info, 'eighs', module=this_module_long, procedure='test_eighs_work_rsp')

        allocate(work(eighs_work_size(test_size)))
        do trial = 1, 2
            call zero_basis(X)
            call eighs(A, X, eigvals, residuals, info, x0=x0, tolerance=atol_sp, work=work)
            call check_info(info, 'eighs', module=this_module_long, procedure='test_eighs_work_rsp')
            err = maxval(abs(eigvals - eigvals_ref))
            call get_err_str(msg, "max err: ", err)
            call check(error, err < rtol_sp)
            call check_test(error, 'test_eighs_work_rsp', &
                          & info=trial_info(trial), eq='evals(work) == evals(internal)', context=msg)
        enddo
    end subroutine test_eighs_work_rsp

    subroutine test_svds_work_rsp(error)
        type(error_type), allocatable, intent(out) :: error
        type(linop_rsp), allocatable :: A
        type(vector_rsp), allocatable :: U(:), V(:), u0
        type(vector_rsp), allocatable :: uwork(:), vwork(:)
        real(sp), allocatable :: S_ref(:), S(:), residuals(:)
        integer :: i, info, trial
        real(sp) :: err
        character(len=256) :: msg

        ! Strang matrix.
        A = linop_rsp() ; A%data = zero_rsp
        do i = 1, test_size
            A%data(i, i) = 2*one_rsp
            if (i < test_size) then
                A%data(i, i+1) = -one_rsp
                A%data(i+1, i) = -one_rsp
            endif
        enddo
        u0 = vector_rsp() ; call init_rand(u0)

        allocate(U(test_size), V(test_size)) ; call zero_basis(U) ; call zero_basis(V)
        call svds(A, U, S_ref, V, residuals, info, u0=u0, tolerance=atol_sp)
        call check_info(info, 'svds', module=this_module_long, procedure='test_svds_work_rsp')

        allocate(uwork(svds_work_size(test_size)), vwork(svds_work_size(test_size)))
        do trial = 1, 2
            call zero_basis(U) ; call zero_basis(V)
            call svds(A, U, S, V, residuals, info, u0=u0, tolerance=atol_sp, uwork=uwork, vwork=vwork)
            call check_info(info, 'svds', module=this_module_long, procedure='test_svds_work_rsp')
            err = maxval(abs(S - S_ref))
            call get_err_str(msg, "max err: ", err)
            call check(error, err < rtol_sp)
            call check_test(error, 'test_svds_work_rsp', &
                          & info=trial_info(trial), eq='S(work) == S(internal)', context=msg)
        enddo
    end subroutine test_svds_work_rsp

    !-------------------------------------------------------------
    !-----     MATRIX EXPONENTIAL AND NEWTON WITH A WORKSPACE     -----
    !-------------------------------------------------------------

    subroutine test_kexpm_work_rsp(error)
        type(error_type), allocatable, intent(out) :: error
        type(linop_rsp), allocatable :: A
        type(vector_rsp), allocatable :: Q, Xref, X
        type(vector_rsp), allocatable :: work(:)
        real(sp), parameter :: tau = 0.1_sp
        integer, parameter :: nkmax = 64
        integer :: info, trial
        real(sp) :: err
        character(len=256) :: msg

        A = linop_rsp() ; call init_rand(A)
        allocate(Q) ; call init_rand(Q)
        allocate(Xref) ; call Xref%zero()

        call kexpm(Xref, A, Q, tau, rtol_sp, info, kdim=nkmax)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexpm_work_rsp')

        allocate(work(kexpm_work_size(nkmax)))
        do trial = 1, 2
            X = vector_rsp() ; call X%zero()
            call kexpm(X, A, Q, tau, rtol_sp, info, kdim=nkmax, work=work)
            call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexpm_work_rsp')
            call X%sub(Xref) ; err = X%norm() / Xref%norm()
            call get_err_str(msg, "max err: ", err)
            call check(error, err < rtol_sp)
            call check_test(error, 'test_kexpm_work_rsp', &
                          & info=trial_info(trial), eq='expm(work) == expm(internal)', context=msg)
        enddo
    end subroutine test_kexpm_work_rsp

    subroutine test_newton_work_rsp(error)
        !! Newton-Krylov on the Roessler system, the workspace holding both the Newton work
        !! vectors and the GMRES workspace forwarded to the linear solver.
        type(error_type), allocatable, intent(out) :: error
        type(roessler_rsp), allocatable :: sys
        type(state_vector_rsp), allocatable :: X, fp1, fp2
        type(state_vector_rsp), allocatable :: work(:)
        type(newton_sp_opts) :: opts
        integer :: info, trial
        real(sp) :: err
        character(len=256) :: msg

        allocate(X, fp1, fp2)
        sys = roessler_rsp()
        sys%jacobian = jacobian_rsp()
        call roessler_analytical_fp_rsp(fp1, fp2)

        allocate(work(newton_work_size(gmres_work_size())))
        do trial = 1, 2
            ! Second trial also exercises the step bisection work vectors.
            opts = newton_sp_opts(maxiter=60, ifbisect=(trial == 2))
            X%x = zero_rsp ; X%y = zero_rsp ; X%z = zero_rsp
            call newton(sys, X, gmres_rsp, info, rtol=10*atol_sp, atol=10*atol_sp, &
                      & options=opts, work=work)
            call X%sub(fp1) ; err = X%norm()
            call get_err_str(msg, "max err: ", err)
            call check(error, err < rtol_sp)
            call check_test(error, 'test_newton_work_rsp', &
                          & info=newton_trial_info(trial), eq='X == analytical fixed point', context=msg)
        enddo
    end subroutine test_newton_work_rsp

    subroutine collect_workspace_rdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [ &
            new_unittest("GMRES with workspace", test_gmres_work_rdp), &
            new_unittest("FGMRES with workspace", test_fgmres_work_rdp), &
            new_unittest("CG with workspace", test_cg_work_rdp), &
            new_unittest("KS eigs with workspace", test_eigs_work_rdp), &
            new_unittest("Sym. eighs with workspace", test_eighs_work_rdp), &
            new_unittest("SVDS with workspace", test_svds_work_rdp), &
            new_unittest("Krylov expm with workspace", test_kexpm_work_rdp), &
            new_unittest("Newton-Krylov with workspace", test_newton_work_rdp) &
        ]
    end subroutine collect_workspace_rdp_testsuite

    !-----------------------------------------------------
    !-----     LINEAR SOLVERS WITH A WORKSPACE         -----
    !-----------------------------------------------------

    subroutine test_gmres_work_rdp(error)
        type(error_type), allocatable, intent(out) :: error
        type(linop_rdp), allocatable :: A
        type(vector_rdp), allocatable :: b, xref, x
        type(vector_rdp), allocatable :: work(:)
        type(gmres_dp_opts) :: opts
        integer :: info, trial
        real(dp) :: err
        character(len=256) :: msg

        A = linop_rdp() ; call init_rand(A)
        b = vector_rdp() ; call init_rand(b)
        xref = vector_rdp() ; call xref%zero()
        opts = gmres_dp_opts(kdim=test_size)

        ! Reference: workspace allocated by the solver.
        call gmres(A, b, xref, info, rtol=rtol_dp, atol=atol_dp, options=opts)
        call check_info(info, 'gmres', module=this_module_long, procedure='test_gmres_work_rdp')

        ! Caller-provided workspace, used for two consecutive solves.
        allocate(work(gmres_work_size(test_size)))
        do trial = 1, 2
            x = vector_rdp() ; call x%zero()
            call gmres(A, b, x, info, rtol=rtol_dp, atol=atol_dp, options=opts, work=work)
            call check_info(info, 'gmres', module=this_module_long, procedure='test_gmres_work_rdp')
            call x%sub(xref) ; err = x%norm() / xref%norm()
            call get_err_str(msg, "max err: ", err)
            call check(error, err < rtol_dp)
            call check_test(error, 'test_gmres_work_rdp', &
                          & info=trial_info(trial), eq='x(work) == x(internal)', context=msg)
        enddo
    end subroutine test_gmres_work_rdp

    subroutine test_fgmres_work_rdp(error)
        type(error_type), allocatable, intent(out) :: error
        type(linop_rdp), allocatable :: A
        type(vector_rdp), allocatable :: b, xref, x
        type(vector_rdp), allocatable :: work(:)
        type(fgmres_dp_opts) :: opts
        integer :: info, trial
        real(dp) :: err
        character(len=256) :: msg

        A = linop_rdp() ; call init_rand(A)
        b = vector_rdp() ; call init_rand(b)
        xref = vector_rdp() ; call xref%zero()
        opts = fgmres_dp_opts(kdim=test_size)

        call fgmres(A, b, xref, info, rtol=rtol_dp, atol=atol_dp, options=opts)
        call check_info(info, 'fgmres', module=this_module_long, procedure='test_fgmres_work_rdp')

        allocate(work(fgmres_work_size(test_size)))
        do trial = 1, 2
            x = vector_rdp() ; call x%zero()
            call fgmres(A, b, x, info, rtol=rtol_dp, atol=atol_dp, options=opts, work=work)
            call check_info(info, 'fgmres', module=this_module_long, procedure='test_fgmres_work_rdp')
            call x%sub(xref) ; err = x%norm() / xref%norm()
            call get_err_str(msg, "max err: ", err)
            call check(error, err < rtol_dp)
            call check_test(error, 'test_fgmres_work_rdp', &
                          & info=trial_info(trial), eq='x(work) == x(internal)', context=msg)
        enddo
    end subroutine test_fgmres_work_rdp

    subroutine test_cg_work_rdp(error)
        type(error_type), allocatable, intent(out) :: error
        type(spd_linop_rdp), allocatable :: A
        type(vector_rdp), allocatable :: b, xref, x
        type(vector_rdp), allocatable :: work(:)
        type(cg_dp_opts) :: opts
        integer :: info, trial
        real(dp) :: err
        character(len=256) :: msg

        A = spd_linop_rdp() ; call init_rand(A)
        b = vector_rdp() ; call init_rand(b)
        xref = vector_rdp() ; call xref%zero()
        opts = cg_dp_opts(maxiter=2*test_size)

        call cg(A, b, xref, info, rtol=rtol_dp, atol=atol_dp, options=opts)
        call check_info(info, 'cg', module=this_module_long, procedure='test_cg_work_rdp')

        allocate(work(cg_work_size()))
        do trial = 1, 2
            x = vector_rdp() ; call x%zero()
            call cg(A, b, x, info, rtol=rtol_dp, atol=atol_dp, options=opts, work=work)
            call check_info(info, 'cg', module=this_module_long, procedure='test_cg_work_rdp')
            call x%sub(xref) ; err = x%norm() / xref%norm()
            call get_err_str(msg, "max err: ", err)
            call check(error, err < rtol_dp)
            call check_test(error, 'test_cg_work_rdp', &
                          & info=trial_info(trial), eq='x(work) == x(internal)', context=msg)
        enddo
    end subroutine test_cg_work_rdp

    !--------------------------------------------------------
    !-----     EIGENVALUE / SVD SOLVERS WITH A WORKSPACE     -----
    !--------------------------------------------------------

    subroutine test_eigs_work_rdp(error)
        !! Same Toeplitz problem as the Krylov-Schur eigs test (nev < n, so restarts happen
        !! and the Krylov-Schur part of the workspace is used).
        type(error_type), allocatable, intent(out) :: error
        type(linop_rdp), allocatable :: A
        type(vector_rdp), allocatable :: X(:), x0
        type(vector_rdp), allocatable :: work(:)
        complex(dp), allocatable :: eigvals_ref(:), eigvals(:)
        real(dp), allocatable :: residuals(:)
        integer, parameter :: nev = 8
        integer :: i, n, info, trial
        real(dp) :: a_, b_
        real(dp) :: err
        character(len=256) :: msg

        A = linop_rdp() ; A%data = 0.0_dp ; n = size(A%data, 1)
        call random_number(a_) ; call random_number(b_) ; b_ = abs(b_)
        do i = 1, n
            A%data(i, i) = a_
            if (i < n) then
                A%data(i, i+1) = b_
                A%data(i+1, i) = -b_
            endif
        enddo
        x0 = vector_rdp() ; call init_rand(x0)

        allocate(X(nev)) ; call zero_basis(X)
        call eigs(A, X, eigvals_ref, residuals, info, x0=x0, tolerance=atol_dp)
        call check_info(info, 'eigs', module=this_module_long, procedure='test_eigs_work_rdp')

        allocate(work(eigs_work_size(nev)))
        do trial = 1, 2
            call zero_basis(X)
            call eigs(A, X, eigvals, residuals, info, x0=x0, tolerance=atol_dp, work=work)
            call check_info(info, 'eigs', module=this_module_long, procedure='test_eigs_work_rdp')
            err = maxval(abs(eigvals - eigvals_ref))
            call get_err_str(msg, "max err: ", err)
            call check(error, err < rtol_dp)
            call check_test(error, 'test_eigs_work_rdp', &
                          & info=trial_info(trial), eq='evals(work) == evals(internal)', context=msg)
        enddo
    end subroutine test_eigs_work_rdp

    subroutine test_eighs_work_rdp(error)
        type(error_type), allocatable, intent(out) :: error
        type(spd_linop_rdp), allocatable :: A
        type(vector_rdp), allocatable :: X(:), x0
        type(vector_rdp), allocatable :: work(:)
        real(dp), allocatable :: eigvals_ref(:), eigvals(:), residuals(:)
        integer :: i, info, trial
        real(dp) :: err
        character(len=256) :: msg

        ! Strang matrix (symmetric positive definite, real-valued).
        A = spd_linop_rdp()
        A%data = zero_rdp
        do i = 1, test_size
            A%data(i, i) = 2*one_rdp
            if (i < test_size) then
                A%data(i, i+1) = -one_rdp
                A%data(i+1, i) = -one_rdp
            endif
        enddo
        x0 = vector_rdp() ; call init_rand(x0)

        allocate(X(test_size)) ; call zero_basis(X)
        call eighs(A, X, eigvals_ref, residuals, info, x0=x0, tolerance=atol_dp)
        call check_info(info, 'eighs', module=this_module_long, procedure='test_eighs_work_rdp')

        allocate(work(eighs_work_size(test_size)))
        do trial = 1, 2
            call zero_basis(X)
            call eighs(A, X, eigvals, residuals, info, x0=x0, tolerance=atol_dp, work=work)
            call check_info(info, 'eighs', module=this_module_long, procedure='test_eighs_work_rdp')
            err = maxval(abs(eigvals - eigvals_ref))
            call get_err_str(msg, "max err: ", err)
            call check(error, err < rtol_dp)
            call check_test(error, 'test_eighs_work_rdp', &
                          & info=trial_info(trial), eq='evals(work) == evals(internal)', context=msg)
        enddo
    end subroutine test_eighs_work_rdp

    subroutine test_svds_work_rdp(error)
        type(error_type), allocatable, intent(out) :: error
        type(linop_rdp), allocatable :: A
        type(vector_rdp), allocatable :: U(:), V(:), u0
        type(vector_rdp), allocatable :: uwork(:), vwork(:)
        real(dp), allocatable :: S_ref(:), S(:), residuals(:)
        integer :: i, info, trial
        real(dp) :: err
        character(len=256) :: msg

        ! Strang matrix.
        A = linop_rdp() ; A%data = zero_rdp
        do i = 1, test_size
            A%data(i, i) = 2*one_rdp
            if (i < test_size) then
                A%data(i, i+1) = -one_rdp
                A%data(i+1, i) = -one_rdp
            endif
        enddo
        u0 = vector_rdp() ; call init_rand(u0)

        allocate(U(test_size), V(test_size)) ; call zero_basis(U) ; call zero_basis(V)
        call svds(A, U, S_ref, V, residuals, info, u0=u0, tolerance=atol_dp)
        call check_info(info, 'svds', module=this_module_long, procedure='test_svds_work_rdp')

        allocate(uwork(svds_work_size(test_size)), vwork(svds_work_size(test_size)))
        do trial = 1, 2
            call zero_basis(U) ; call zero_basis(V)
            call svds(A, U, S, V, residuals, info, u0=u0, tolerance=atol_dp, uwork=uwork, vwork=vwork)
            call check_info(info, 'svds', module=this_module_long, procedure='test_svds_work_rdp')
            err = maxval(abs(S - S_ref))
            call get_err_str(msg, "max err: ", err)
            call check(error, err < rtol_dp)
            call check_test(error, 'test_svds_work_rdp', &
                          & info=trial_info(trial), eq='S(work) == S(internal)', context=msg)
        enddo
    end subroutine test_svds_work_rdp

    !-------------------------------------------------------------
    !-----     MATRIX EXPONENTIAL AND NEWTON WITH A WORKSPACE     -----
    !-------------------------------------------------------------

    subroutine test_kexpm_work_rdp(error)
        type(error_type), allocatable, intent(out) :: error
        type(linop_rdp), allocatable :: A
        type(vector_rdp), allocatable :: Q, Xref, X
        type(vector_rdp), allocatable :: work(:)
        real(dp), parameter :: tau = 0.1_dp
        integer, parameter :: nkmax = 64
        integer :: info, trial
        real(dp) :: err
        character(len=256) :: msg

        A = linop_rdp() ; call init_rand(A)
        allocate(Q) ; call init_rand(Q)
        allocate(Xref) ; call Xref%zero()

        call kexpm(Xref, A, Q, tau, rtol_dp, info, kdim=nkmax)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexpm_work_rdp')

        allocate(work(kexpm_work_size(nkmax)))
        do trial = 1, 2
            X = vector_rdp() ; call X%zero()
            call kexpm(X, A, Q, tau, rtol_dp, info, kdim=nkmax, work=work)
            call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexpm_work_rdp')
            call X%sub(Xref) ; err = X%norm() / Xref%norm()
            call get_err_str(msg, "max err: ", err)
            call check(error, err < rtol_dp)
            call check_test(error, 'test_kexpm_work_rdp', &
                          & info=trial_info(trial), eq='expm(work) == expm(internal)', context=msg)
        enddo
    end subroutine test_kexpm_work_rdp

    subroutine test_newton_work_rdp(error)
        !! Newton-Krylov on the Roessler system, the workspace holding both the Newton work
        !! vectors and the GMRES workspace forwarded to the linear solver.
        type(error_type), allocatable, intent(out) :: error
        type(roessler_rdp), allocatable :: sys
        type(state_vector_rdp), allocatable :: X, fp1, fp2
        type(state_vector_rdp), allocatable :: work(:)
        type(newton_dp_opts) :: opts
        integer :: info, trial
        real(dp) :: err
        character(len=256) :: msg

        allocate(X, fp1, fp2)
        sys = roessler_rdp()
        sys%jacobian = jacobian_rdp()
        call roessler_analytical_fp_rdp(fp1, fp2)

        allocate(work(newton_work_size(gmres_work_size())))
        do trial = 1, 2
            ! Second trial also exercises the step bisection work vectors.
            opts = newton_dp_opts(maxiter=60, ifbisect=(trial == 2))
            X%x = zero_rdp ; X%y = zero_rdp ; X%z = zero_rdp
            call newton(sys, X, gmres_rdp, info, rtol=10*atol_dp, atol=10*atol_dp, &
                      & options=opts, work=work)
            call X%sub(fp1) ; err = X%norm()
            call get_err_str(msg, "max err: ", err)
            call check(error, err < rtol_dp)
            call check_test(error, 'test_newton_work_rdp', &
                          & info=newton_trial_info(trial), eq='X == analytical fixed point', context=msg)
        enddo
    end subroutine test_newton_work_rdp

    subroutine collect_workspace_csp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [ &
            new_unittest("GMRES with workspace", test_gmres_work_csp), &
            new_unittest("FGMRES with workspace", test_fgmres_work_csp), &
            new_unittest("CG with workspace", test_cg_work_csp), &
            new_unittest("Sym. eighs with workspace", test_eighs_work_csp), &
            new_unittest("SVDS with workspace", test_svds_work_csp), &
            new_unittest("Krylov expm with workspace", test_kexpm_work_csp), &
            new_unittest("Newton-Krylov with workspace", test_newton_work_csp) &
        ]
    end subroutine collect_workspace_csp_testsuite

    !-----------------------------------------------------
    !-----     LINEAR SOLVERS WITH A WORKSPACE         -----
    !-----------------------------------------------------

    subroutine test_gmres_work_csp(error)
        type(error_type), allocatable, intent(out) :: error
        type(linop_csp), allocatable :: A
        type(vector_csp), allocatable :: b, xref, x
        type(vector_csp), allocatable :: work(:)
        type(gmres_sp_opts) :: opts
        integer :: info, trial
        real(sp) :: err
        character(len=256) :: msg

        A = linop_csp() ; call init_rand(A)
        b = vector_csp() ; call init_rand(b)
        xref = vector_csp() ; call xref%zero()
        opts = gmres_sp_opts(kdim=test_size)

        ! Reference: workspace allocated by the solver.
        call gmres(A, b, xref, info, rtol=rtol_sp, atol=atol_sp, options=opts)
        call check_info(info, 'gmres', module=this_module_long, procedure='test_gmres_work_csp')

        ! Caller-provided workspace, used for two consecutive solves.
        allocate(work(gmres_work_size(test_size)))
        do trial = 1, 2
            x = vector_csp() ; call x%zero()
            call gmres(A, b, x, info, rtol=rtol_sp, atol=atol_sp, options=opts, work=work)
            call check_info(info, 'gmres', module=this_module_long, procedure='test_gmres_work_csp')
            call x%sub(xref) ; err = x%norm() / xref%norm()
            call get_err_str(msg, "max err: ", err)
            call check(error, err < rtol_sp)
            call check_test(error, 'test_gmres_work_csp', &
                          & info=trial_info(trial), eq='x(work) == x(internal)', context=msg)
        enddo
    end subroutine test_gmres_work_csp

    subroutine test_fgmres_work_csp(error)
        type(error_type), allocatable, intent(out) :: error
        type(linop_csp), allocatable :: A
        type(vector_csp), allocatable :: b, xref, x
        type(vector_csp), allocatable :: work(:)
        type(fgmres_sp_opts) :: opts
        integer :: info, trial
        real(sp) :: err
        character(len=256) :: msg

        A = linop_csp() ; call init_rand(A)
        b = vector_csp() ; call init_rand(b)
        xref = vector_csp() ; call xref%zero()
        opts = fgmres_sp_opts(kdim=test_size)

        call fgmres(A, b, xref, info, rtol=rtol_sp, atol=atol_sp, options=opts)
        call check_info(info, 'fgmres', module=this_module_long, procedure='test_fgmres_work_csp')

        allocate(work(fgmres_work_size(test_size)))
        do trial = 1, 2
            x = vector_csp() ; call x%zero()
            call fgmres(A, b, x, info, rtol=rtol_sp, atol=atol_sp, options=opts, work=work)
            call check_info(info, 'fgmres', module=this_module_long, procedure='test_fgmres_work_csp')
            call x%sub(xref) ; err = x%norm() / xref%norm()
            call get_err_str(msg, "max err: ", err)
            call check(error, err < rtol_sp)
            call check_test(error, 'test_fgmres_work_csp', &
                          & info=trial_info(trial), eq='x(work) == x(internal)', context=msg)
        enddo
    end subroutine test_fgmres_work_csp

    subroutine test_cg_work_csp(error)
        type(error_type), allocatable, intent(out) :: error
        type(hermitian_linop_csp), allocatable :: A
        type(vector_csp), allocatable :: b, xref, x
        type(vector_csp), allocatable :: work(:)
        type(cg_sp_opts) :: opts
        integer :: info, trial
        real(sp) :: err
        character(len=256) :: msg

        A = hermitian_linop_csp() ; call init_rand(A)
        b = vector_csp() ; call init_rand(b)
        xref = vector_csp() ; call xref%zero()
        opts = cg_sp_opts(maxiter=2*test_size)

        call cg(A, b, xref, info, rtol=rtol_sp, atol=atol_sp, options=opts)
        call check_info(info, 'cg', module=this_module_long, procedure='test_cg_work_csp')

        allocate(work(cg_work_size()))
        do trial = 1, 2
            x = vector_csp() ; call x%zero()
            call cg(A, b, x, info, rtol=rtol_sp, atol=atol_sp, options=opts, work=work)
            call check_info(info, 'cg', module=this_module_long, procedure='test_cg_work_csp')
            call x%sub(xref) ; err = x%norm() / xref%norm()
            call get_err_str(msg, "max err: ", err)
            call check(error, err < rtol_sp)
            call check_test(error, 'test_cg_work_csp', &
                          & info=trial_info(trial), eq='x(work) == x(internal)', context=msg)
        enddo
    end subroutine test_cg_work_csp

    !--------------------------------------------------------
    !-----     EIGENVALUE / SVD SOLVERS WITH A WORKSPACE     -----
    !--------------------------------------------------------


    subroutine test_eighs_work_csp(error)
        type(error_type), allocatable, intent(out) :: error
        type(hermitian_linop_csp), allocatable :: A
        type(vector_csp), allocatable :: X(:), x0
        type(vector_csp), allocatable :: work(:)
        real(sp), allocatable :: eigvals_ref(:), eigvals(:), residuals(:)
        integer :: i, info, trial
        real(sp) :: err
        character(len=256) :: msg

        ! Strang matrix (symmetric positive definite, real-valued).
        A = hermitian_linop_csp()
        A%data = zero_csp
        do i = 1, test_size
            A%data(i, i) = 2*one_csp
            if (i < test_size) then
                A%data(i, i+1) = -one_csp
                A%data(i+1, i) = -one_csp
            endif
        enddo
        x0 = vector_csp() ; call init_rand(x0)

        allocate(X(test_size)) ; call zero_basis(X)
        call eighs(A, X, eigvals_ref, residuals, info, x0=x0, tolerance=atol_sp)
        call check_info(info, 'eighs', module=this_module_long, procedure='test_eighs_work_csp')

        allocate(work(eighs_work_size(test_size)))
        do trial = 1, 2
            call zero_basis(X)
            call eighs(A, X, eigvals, residuals, info, x0=x0, tolerance=atol_sp, work=work)
            call check_info(info, 'eighs', module=this_module_long, procedure='test_eighs_work_csp')
            err = maxval(abs(eigvals - eigvals_ref))
            call get_err_str(msg, "max err: ", err)
            call check(error, err < rtol_sp)
            call check_test(error, 'test_eighs_work_csp', &
                          & info=trial_info(trial), eq='evals(work) == evals(internal)', context=msg)
        enddo
    end subroutine test_eighs_work_csp

    subroutine test_svds_work_csp(error)
        type(error_type), allocatable, intent(out) :: error
        type(linop_csp), allocatable :: A
        type(vector_csp), allocatable :: U(:), V(:), u0
        type(vector_csp), allocatable :: uwork(:), vwork(:)
        real(sp), allocatable :: S_ref(:), S(:), residuals(:)
        integer :: i, info, trial
        real(sp) :: err
        character(len=256) :: msg

        ! Strang matrix.
        A = linop_csp() ; A%data = zero_csp
        do i = 1, test_size
            A%data(i, i) = 2*one_csp
            if (i < test_size) then
                A%data(i, i+1) = -one_csp
                A%data(i+1, i) = -one_csp
            endif
        enddo
        u0 = vector_csp() ; call init_rand(u0)

        allocate(U(test_size), V(test_size)) ; call zero_basis(U) ; call zero_basis(V)
        call svds(A, U, S_ref, V, residuals, info, u0=u0, tolerance=atol_sp)
        call check_info(info, 'svds', module=this_module_long, procedure='test_svds_work_csp')

        allocate(uwork(svds_work_size(test_size)), vwork(svds_work_size(test_size)))
        do trial = 1, 2
            call zero_basis(U) ; call zero_basis(V)
            call svds(A, U, S, V, residuals, info, u0=u0, tolerance=atol_sp, uwork=uwork, vwork=vwork)
            call check_info(info, 'svds', module=this_module_long, procedure='test_svds_work_csp')
            err = maxval(abs(S - S_ref))
            call get_err_str(msg, "max err: ", err)
            call check(error, err < rtol_sp)
            call check_test(error, 'test_svds_work_csp', &
                          & info=trial_info(trial), eq='S(work) == S(internal)', context=msg)
        enddo
    end subroutine test_svds_work_csp

    !-------------------------------------------------------------
    !-----     MATRIX EXPONENTIAL AND NEWTON WITH A WORKSPACE     -----
    !-------------------------------------------------------------

    subroutine test_kexpm_work_csp(error)
        type(error_type), allocatable, intent(out) :: error
        type(linop_csp), allocatable :: A
        type(vector_csp), allocatable :: Q, Xref, X
        type(vector_csp), allocatable :: work(:)
        real(sp), parameter :: tau = 0.1_sp
        integer, parameter :: nkmax = 64
        integer :: info, trial
        real(sp) :: err
        character(len=256) :: msg

        A = linop_csp() ; call init_rand(A)
        allocate(Q) ; call init_rand(Q)
        allocate(Xref) ; call Xref%zero()

        call kexpm(Xref, A, Q, tau, rtol_sp, info, kdim=nkmax)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexpm_work_csp')

        allocate(work(kexpm_work_size(nkmax)))
        do trial = 1, 2
            X = vector_csp() ; call X%zero()
            call kexpm(X, A, Q, tau, rtol_sp, info, kdim=nkmax, work=work)
            call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexpm_work_csp')
            call X%sub(Xref) ; err = X%norm() / Xref%norm()
            call get_err_str(msg, "max err: ", err)
            call check(error, err < rtol_sp)
            call check_test(error, 'test_kexpm_work_csp', &
                          & info=trial_info(trial), eq='expm(work) == expm(internal)', context=msg)
        enddo
    end subroutine test_kexpm_work_csp

    subroutine test_newton_work_csp(error)
        !! Newton-Krylov on the Roessler system, the workspace holding both the Newton work
        !! vectors and the GMRES workspace forwarded to the linear solver.
        type(error_type), allocatable, intent(out) :: error
        type(roessler_csp), allocatable :: sys
        type(state_vector_csp), allocatable :: X, fp1, fp2
        type(state_vector_csp), allocatable :: work(:)
        type(newton_sp_opts) :: opts
        integer :: info, trial
        real(sp) :: err
        character(len=256) :: msg

        allocate(X, fp1, fp2)
        sys = roessler_csp()
        sys%jacobian = jacobian_csp()
        call roessler_analytical_fp_csp(fp1, fp2)

        allocate(work(newton_work_size(gmres_work_size())))
        do trial = 1, 2
            ! Second trial also exercises the step bisection work vectors.
            opts = newton_sp_opts(maxiter=60, ifbisect=(trial == 2))
            X%x = zero_csp ; X%y = zero_csp ; X%z = zero_csp
            call newton(sys, X, gmres_csp, info, rtol=10*atol_sp, atol=10*atol_sp, &
                      & options=opts, work=work)
            call X%sub(fp1) ; err = X%norm()
            call get_err_str(msg, "max err: ", err)
            call check(error, err < rtol_sp)
            call check_test(error, 'test_newton_work_csp', &
                          & info=newton_trial_info(trial), eq='X == analytical fixed point', context=msg)
        enddo
    end subroutine test_newton_work_csp

    subroutine collect_workspace_cdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [ &
            new_unittest("GMRES with workspace", test_gmres_work_cdp), &
            new_unittest("FGMRES with workspace", test_fgmres_work_cdp), &
            new_unittest("CG with workspace", test_cg_work_cdp), &
            new_unittest("Sym. eighs with workspace", test_eighs_work_cdp), &
            new_unittest("SVDS with workspace", test_svds_work_cdp), &
            new_unittest("Krylov expm with workspace", test_kexpm_work_cdp), &
            new_unittest("Newton-Krylov with workspace", test_newton_work_cdp) &
        ]
    end subroutine collect_workspace_cdp_testsuite

    !-----------------------------------------------------
    !-----     LINEAR SOLVERS WITH A WORKSPACE         -----
    !-----------------------------------------------------

    subroutine test_gmres_work_cdp(error)
        type(error_type), allocatable, intent(out) :: error
        type(linop_cdp), allocatable :: A
        type(vector_cdp), allocatable :: b, xref, x
        type(vector_cdp), allocatable :: work(:)
        type(gmres_dp_opts) :: opts
        integer :: info, trial
        real(dp) :: err
        character(len=256) :: msg

        A = linop_cdp() ; call init_rand(A)
        b = vector_cdp() ; call init_rand(b)
        xref = vector_cdp() ; call xref%zero()
        opts = gmres_dp_opts(kdim=test_size)

        ! Reference: workspace allocated by the solver.
        call gmres(A, b, xref, info, rtol=rtol_dp, atol=atol_dp, options=opts)
        call check_info(info, 'gmres', module=this_module_long, procedure='test_gmres_work_cdp')

        ! Caller-provided workspace, used for two consecutive solves.
        allocate(work(gmres_work_size(test_size)))
        do trial = 1, 2
            x = vector_cdp() ; call x%zero()
            call gmres(A, b, x, info, rtol=rtol_dp, atol=atol_dp, options=opts, work=work)
            call check_info(info, 'gmres', module=this_module_long, procedure='test_gmres_work_cdp')
            call x%sub(xref) ; err = x%norm() / xref%norm()
            call get_err_str(msg, "max err: ", err)
            call check(error, err < rtol_dp)
            call check_test(error, 'test_gmres_work_cdp', &
                          & info=trial_info(trial), eq='x(work) == x(internal)', context=msg)
        enddo
    end subroutine test_gmres_work_cdp

    subroutine test_fgmres_work_cdp(error)
        type(error_type), allocatable, intent(out) :: error
        type(linop_cdp), allocatable :: A
        type(vector_cdp), allocatable :: b, xref, x
        type(vector_cdp), allocatable :: work(:)
        type(fgmres_dp_opts) :: opts
        integer :: info, trial
        real(dp) :: err
        character(len=256) :: msg

        A = linop_cdp() ; call init_rand(A)
        b = vector_cdp() ; call init_rand(b)
        xref = vector_cdp() ; call xref%zero()
        opts = fgmres_dp_opts(kdim=test_size)

        call fgmres(A, b, xref, info, rtol=rtol_dp, atol=atol_dp, options=opts)
        call check_info(info, 'fgmres', module=this_module_long, procedure='test_fgmres_work_cdp')

        allocate(work(fgmres_work_size(test_size)))
        do trial = 1, 2
            x = vector_cdp() ; call x%zero()
            call fgmres(A, b, x, info, rtol=rtol_dp, atol=atol_dp, options=opts, work=work)
            call check_info(info, 'fgmres', module=this_module_long, procedure='test_fgmres_work_cdp')
            call x%sub(xref) ; err = x%norm() / xref%norm()
            call get_err_str(msg, "max err: ", err)
            call check(error, err < rtol_dp)
            call check_test(error, 'test_fgmres_work_cdp', &
                          & info=trial_info(trial), eq='x(work) == x(internal)', context=msg)
        enddo
    end subroutine test_fgmres_work_cdp

    subroutine test_cg_work_cdp(error)
        type(error_type), allocatable, intent(out) :: error
        type(hermitian_linop_cdp), allocatable :: A
        type(vector_cdp), allocatable :: b, xref, x
        type(vector_cdp), allocatable :: work(:)
        type(cg_dp_opts) :: opts
        integer :: info, trial
        real(dp) :: err
        character(len=256) :: msg

        A = hermitian_linop_cdp() ; call init_rand(A)
        b = vector_cdp() ; call init_rand(b)
        xref = vector_cdp() ; call xref%zero()
        opts = cg_dp_opts(maxiter=2*test_size)

        call cg(A, b, xref, info, rtol=rtol_dp, atol=atol_dp, options=opts)
        call check_info(info, 'cg', module=this_module_long, procedure='test_cg_work_cdp')

        allocate(work(cg_work_size()))
        do trial = 1, 2
            x = vector_cdp() ; call x%zero()
            call cg(A, b, x, info, rtol=rtol_dp, atol=atol_dp, options=opts, work=work)
            call check_info(info, 'cg', module=this_module_long, procedure='test_cg_work_cdp')
            call x%sub(xref) ; err = x%norm() / xref%norm()
            call get_err_str(msg, "max err: ", err)
            call check(error, err < rtol_dp)
            call check_test(error, 'test_cg_work_cdp', &
                          & info=trial_info(trial), eq='x(work) == x(internal)', context=msg)
        enddo
    end subroutine test_cg_work_cdp

    !--------------------------------------------------------
    !-----     EIGENVALUE / SVD SOLVERS WITH A WORKSPACE     -----
    !--------------------------------------------------------


    subroutine test_eighs_work_cdp(error)
        type(error_type), allocatable, intent(out) :: error
        type(hermitian_linop_cdp), allocatable :: A
        type(vector_cdp), allocatable :: X(:), x0
        type(vector_cdp), allocatable :: work(:)
        real(dp), allocatable :: eigvals_ref(:), eigvals(:), residuals(:)
        integer :: i, info, trial
        real(dp) :: err
        character(len=256) :: msg

        ! Strang matrix (symmetric positive definite, real-valued).
        A = hermitian_linop_cdp()
        A%data = zero_cdp
        do i = 1, test_size
            A%data(i, i) = 2*one_cdp
            if (i < test_size) then
                A%data(i, i+1) = -one_cdp
                A%data(i+1, i) = -one_cdp
            endif
        enddo
        x0 = vector_cdp() ; call init_rand(x0)

        allocate(X(test_size)) ; call zero_basis(X)
        call eighs(A, X, eigvals_ref, residuals, info, x0=x0, tolerance=atol_dp)
        call check_info(info, 'eighs', module=this_module_long, procedure='test_eighs_work_cdp')

        allocate(work(eighs_work_size(test_size)))
        do trial = 1, 2
            call zero_basis(X)
            call eighs(A, X, eigvals, residuals, info, x0=x0, tolerance=atol_dp, work=work)
            call check_info(info, 'eighs', module=this_module_long, procedure='test_eighs_work_cdp')
            err = maxval(abs(eigvals - eigvals_ref))
            call get_err_str(msg, "max err: ", err)
            call check(error, err < rtol_dp)
            call check_test(error, 'test_eighs_work_cdp', &
                          & info=trial_info(trial), eq='evals(work) == evals(internal)', context=msg)
        enddo
    end subroutine test_eighs_work_cdp

    subroutine test_svds_work_cdp(error)
        type(error_type), allocatable, intent(out) :: error
        type(linop_cdp), allocatable :: A
        type(vector_cdp), allocatable :: U(:), V(:), u0
        type(vector_cdp), allocatable :: uwork(:), vwork(:)
        real(dp), allocatable :: S_ref(:), S(:), residuals(:)
        integer :: i, info, trial
        real(dp) :: err
        character(len=256) :: msg

        ! Strang matrix.
        A = linop_cdp() ; A%data = zero_cdp
        do i = 1, test_size
            A%data(i, i) = 2*one_cdp
            if (i < test_size) then
                A%data(i, i+1) = -one_cdp
                A%data(i+1, i) = -one_cdp
            endif
        enddo
        u0 = vector_cdp() ; call init_rand(u0)

        allocate(U(test_size), V(test_size)) ; call zero_basis(U) ; call zero_basis(V)
        call svds(A, U, S_ref, V, residuals, info, u0=u0, tolerance=atol_dp)
        call check_info(info, 'svds', module=this_module_long, procedure='test_svds_work_cdp')

        allocate(uwork(svds_work_size(test_size)), vwork(svds_work_size(test_size)))
        do trial = 1, 2
            call zero_basis(U) ; call zero_basis(V)
            call svds(A, U, S, V, residuals, info, u0=u0, tolerance=atol_dp, uwork=uwork, vwork=vwork)
            call check_info(info, 'svds', module=this_module_long, procedure='test_svds_work_cdp')
            err = maxval(abs(S - S_ref))
            call get_err_str(msg, "max err: ", err)
            call check(error, err < rtol_dp)
            call check_test(error, 'test_svds_work_cdp', &
                          & info=trial_info(trial), eq='S(work) == S(internal)', context=msg)
        enddo
    end subroutine test_svds_work_cdp

    !-------------------------------------------------------------
    !-----     MATRIX EXPONENTIAL AND NEWTON WITH A WORKSPACE     -----
    !-------------------------------------------------------------

    subroutine test_kexpm_work_cdp(error)
        type(error_type), allocatable, intent(out) :: error
        type(linop_cdp), allocatable :: A
        type(vector_cdp), allocatable :: Q, Xref, X
        type(vector_cdp), allocatable :: work(:)
        real(dp), parameter :: tau = 0.1_dp
        integer, parameter :: nkmax = 64
        integer :: info, trial
        real(dp) :: err
        character(len=256) :: msg

        A = linop_cdp() ; call init_rand(A)
        allocate(Q) ; call init_rand(Q)
        allocate(Xref) ; call Xref%zero()

        call kexpm(Xref, A, Q, tau, rtol_dp, info, kdim=nkmax)
        call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexpm_work_cdp')

        allocate(work(kexpm_work_size(nkmax)))
        do trial = 1, 2
            X = vector_cdp() ; call X%zero()
            call kexpm(X, A, Q, tau, rtol_dp, info, kdim=nkmax, work=work)
            call check_info(info, 'kexpm', module=this_module_long, procedure='test_kexpm_work_cdp')
            call X%sub(Xref) ; err = X%norm() / Xref%norm()
            call get_err_str(msg, "max err: ", err)
            call check(error, err < rtol_dp)
            call check_test(error, 'test_kexpm_work_cdp', &
                          & info=trial_info(trial), eq='expm(work) == expm(internal)', context=msg)
        enddo
    end subroutine test_kexpm_work_cdp

    subroutine test_newton_work_cdp(error)
        !! Newton-Krylov on the Roessler system, the workspace holding both the Newton work
        !! vectors and the GMRES workspace forwarded to the linear solver.
        type(error_type), allocatable, intent(out) :: error
        type(roessler_cdp), allocatable :: sys
        type(state_vector_cdp), allocatable :: X, fp1, fp2
        type(state_vector_cdp), allocatable :: work(:)
        type(newton_dp_opts) :: opts
        integer :: info, trial
        real(dp) :: err
        character(len=256) :: msg

        allocate(X, fp1, fp2)
        sys = roessler_cdp()
        sys%jacobian = jacobian_cdp()
        call roessler_analytical_fp_cdp(fp1, fp2)

        allocate(work(newton_work_size(gmres_work_size())))
        do trial = 1, 2
            ! Second trial also exercises the step bisection work vectors.
            opts = newton_dp_opts(maxiter=60, ifbisect=(trial == 2))
            X%x = zero_cdp ; X%y = zero_cdp ; X%z = zero_cdp
            call newton(sys, X, gmres_cdp, info, rtol=10*atol_dp, atol=10*atol_dp, &
                      & options=opts, work=work)
            call X%sub(fp1) ; err = X%norm()
            call get_err_str(msg, "max err: ", err)
            call check(error, err < rtol_dp)
            call check_test(error, 'test_newton_work_cdp', &
                          & info=newton_trial_info(trial), eq='X == analytical fixed point', context=msg)
        enddo
    end subroutine test_newton_work_cdp

end module TestWorkspace