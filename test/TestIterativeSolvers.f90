module TestIterativeSolvers
    ! Fortran Standard library.
    use iso_fortran_env
    use stdlib_math, only: is_close, all_close
    use stdlib_linalg, only: eye

    ! LightKrylov
    use LightKrylov

    ! Testdrive
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use TestVectors
    use TestLinops
    use TestUtils
    use TestKrylov

    implicit none
    private

    public :: collect_eig_rsp_testsuite
    public :: collect_svd_rsp_testsuite
    public :: collect_gmres_rsp_testsuite
    public :: collect_cg_rsp_testsuite
    public :: collect_eig_rdp_testsuite
    public :: collect_svd_rdp_testsuite
    public :: collect_gmres_rdp_testsuite
    public :: collect_cg_rdp_testsuite
    public :: collect_eig_csp_testsuite
    public :: collect_svd_csp_testsuite
    public :: collect_gmres_csp_testsuite
    public :: collect_cg_csp_testsuite
    public :: collect_eig_cdp_testsuite
    public :: collect_svd_cdp_testsuite
    public :: collect_gmres_cdp_testsuite
    public :: collect_cg_cdp_testsuite

contains

    !---------------------------------------------------------
    !-----     DEFINITION OF THE UNIT TESTS FOR EIGS     -----
    !---------------------------------------------------------

    subroutine collect_eig_rsp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

         testsuite = [ &
                    new_unittest("Eigs computation", test_evp_rsp), &
                    new_unittest("Sym. eigs computation", test_sym_evp_rsp) &
                    ]
        return
    end subroutine collect_eig_rsp_testsuite

    subroutine test_evp_rsp(error)
        !> Error type.
        type(error_type), allocatable, intent(out) :: error
        !> Test linear operator.
        type(linop_rsp), allocatable :: A
        !> Eigenvectors.
        type(vector_rsp), allocatable :: X(:)
        !> Eigenvalues.
        complex(sp), allocatable :: eigvals(:)
        !> Residuals.
        real(sp), allocatable :: residuals(:)
        !> Information flag.
        integer :: info
        !> Miscellaneous.
        integer :: i, k, n
        complex(sp) :: true_eigvals(test_size)
        real(sp) :: pi = 4.0_sp * atan(1.0_sp)
        real(sp) :: a_, b_

        ! Allocate eigenvectors.
        allocate(X(test_size)) ; call initialize_krylov_subspace(X)
        
        ! Initialize linear operator with random tridiagonal Toeplitz matrix.
        A = linop_rsp() ; A%data = 0.0_sp ; n = size(A%data, 1)

        call random_number(a_) ; call random_number(b_) ; b_ = abs(b_)
        do i = 1, n
            ! Diagonal entry.
            A%data(i, i) = a_
            ! Upper diagonal entry.
            if (i < n) then
                A%data(i, i+1) = b_
                A%data(i+1, i) = -b_
            endif
        enddo

        ! Compute spectral decomposition.
        call eigs(A, X, eigvals, residuals, info)

        ! Analytical eigenvalues.
        true_eigvals = cmplx(0.0_sp, 0.0_sp, kind=sp) ; k = 1
        do i = 1, test_size, 2
            true_eigvals(i) = a_*cmplx(1.0_sp, 0.0_sp, kind=sp) + (2*b_*cos(k*pi/(test_size+1)))*cmplx(0.0_sp, 1.0_sp, kind=sp)
            true_eigvals(i+1) = a_*cmplx(1.0_sp, 0.0_sp, kind=sp) - (2*b_*cos(k*pi/(test_size+1)))*cmplx(0.0_sp, 1.0_sp, kind=sp)
            k = k+1
        enddo

        call check(error, norm2(abs(eigvals - true_eigvals)) < rtol_sp)

        return
    end subroutine test_evp_rsp

    subroutine test_sym_evp_rsp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test matrix.
        type(spd_linop_rsp), allocatable :: A
        ! Eigenvectors.
        type(vector_rsp), allocatable :: X(:)
        ! Eigenvalues.
        real(sp), allocatable :: evals(:)
        ! Residuals.
        real(sp), allocatable :: residuals(:)
        ! Information flag.
        integer :: info
        ! Toeplitz matrix.
        real(sp) :: T(test_size, test_size), a_, b_
        ! Miscellaneous.
        integer :: i
        real(sp) :: alpha, true_evals(test_size)
        real(sp), parameter :: pi = 4.0_sp * atan(1.0_sp)

        ! Create the sym. pos. def. Toeplitz matrix.
        call random_number(a_) ; call random_number(b_) ; b_ = -abs(b_)
        T = 0.0_sp
        do i = 1, test_size
            ! Diagonal entry.
            T(i, i) = a_
            if (i < test_size) then
                ! Upper diagonal entry.
                T(i, i+1) = b_
                ! Lower diagonal entry.
                T(i+1, i) = b_
            endif
        enddo

        ! Allocations.
        A = spd_linop_rsp(T)
        allocate(X(test_size)) ; call initialize_krylov_subspace(X)

        ! Spectral decomposition.
        call eighs(A, X, evals, residuals, info, kdim=test_size)

        ! Analytical eigenvalues.
        true_evals = 0.0_sp
        do i = 1, test_size
            true_evals(i) = a_ + 2*abs(b_) * cos(i*pi/(test_size+1))
        enddo

        ! Check error.
        call check(error, all_close(evals, true_evals, rtol_sp, atol_sp))

        return
    end subroutine test_sym_evp_rsp

    subroutine collect_eig_rdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

         testsuite = [ &
                    new_unittest("Eigs computation", test_evp_rdp), &
                    new_unittest("Sym. eigs computation", test_sym_evp_rdp) &
                    ]
        return
    end subroutine collect_eig_rdp_testsuite

    subroutine test_evp_rdp(error)
        !> Error type.
        type(error_type), allocatable, intent(out) :: error
        !> Test linear operator.
        type(linop_rdp), allocatable :: A
        !> Eigenvectors.
        type(vector_rdp), allocatable :: X(:)
        !> Eigenvalues.
        complex(dp), allocatable :: eigvals(:)
        !> Residuals.
        real(dp), allocatable :: residuals(:)
        !> Information flag.
        integer :: info
        !> Miscellaneous.
        integer :: i, k, n
        complex(dp) :: true_eigvals(test_size)
        real(dp) :: pi = 4.0_dp * atan(1.0_dp)
        real(dp) :: a_, b_

        ! Allocate eigenvectors.
        allocate(X(test_size)) ; call initialize_krylov_subspace(X)
        
        ! Initialize linear operator with random tridiagonal Toeplitz matrix.
        A = linop_rdp() ; A%data = 0.0_dp ; n = size(A%data, 1)

        call random_number(a_) ; call random_number(b_) ; b_ = abs(b_)
        do i = 1, n
            ! Diagonal entry.
            A%data(i, i) = a_
            ! Upper diagonal entry.
            if (i < n) then
                A%data(i, i+1) = b_
                A%data(i+1, i) = -b_
            endif
        enddo

        ! Compute spectral decomposition.
        call eigs(A, X, eigvals, residuals, info)

        ! Analytical eigenvalues.
        true_eigvals = cmplx(0.0_dp, 0.0_dp, kind=dp) ; k = 1
        do i = 1, test_size, 2
            true_eigvals(i) = a_*cmplx(1.0_dp, 0.0_dp, kind=dp) + (2*b_*cos(k*pi/(test_size+1)))*cmplx(0.0_dp, 1.0_dp, kind=dp)
            true_eigvals(i+1) = a_*cmplx(1.0_dp, 0.0_dp, kind=dp) - (2*b_*cos(k*pi/(test_size+1)))*cmplx(0.0_dp, 1.0_dp, kind=dp)
            k = k+1
        enddo

        call check(error, norm2(abs(eigvals - true_eigvals)) < rtol_dp)

        return
    end subroutine test_evp_rdp

    subroutine test_sym_evp_rdp(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test matrix.
        type(spd_linop_rdp), allocatable :: A
        ! Eigenvectors.
        type(vector_rdp), allocatable :: X(:)
        ! Eigenvalues.
        real(dp), allocatable :: evals(:)
        ! Residuals.
        real(dp), allocatable :: residuals(:)
        ! Information flag.
        integer :: info
        ! Toeplitz matrix.
        real(dp) :: T(test_size, test_size), a_, b_
        ! Miscellaneous.
        integer :: i
        real(dp) :: alpha, true_evals(test_size)
        real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)

        ! Create the sym. pos. def. Toeplitz matrix.
        call random_number(a_) ; call random_number(b_) ; b_ = -abs(b_)
        T = 0.0_dp
        do i = 1, test_size
            ! Diagonal entry.
            T(i, i) = a_
            if (i < test_size) then
                ! Upper diagonal entry.
                T(i, i+1) = b_
                ! Lower diagonal entry.
                T(i+1, i) = b_
            endif
        enddo

        ! Allocations.
        A = spd_linop_rdp(T)
        allocate(X(test_size)) ; call initialize_krylov_subspace(X)

        ! Spectral decomposition.
        call eighs(A, X, evals, residuals, info, kdim=test_size)

        ! Analytical eigenvalues.
        true_evals = 0.0_dp
        do i = 1, test_size
            true_evals(i) = a_ + 2*abs(b_) * cos(i*pi/(test_size+1))
        enddo

        ! Check error.
        call check(error, all_close(evals, true_evals, rtol_dp, atol_dp))

        return
    end subroutine test_sym_evp_rdp

    subroutine collect_eig_csp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("Eigs computation", test_evp_csp) &
                    ]
        return
    end subroutine collect_eig_csp_testsuite

    subroutine test_evp_csp(error)
        !> Error type.
        type(error_type), allocatable, intent(out) :: error
        !> Test linear operator.
        type(linop_csp), allocatable :: A
        !> Eigenvectors.
        type(vector_csp), allocatable :: X(:)
        !> Eigenvalues.
        complex(sp), allocatable :: eigvals(:)
        !> Residuals.
        real(sp), allocatable :: residuals(:)
        !> Information flag.
        integer :: info
        !> Miscellaneous.
        integer :: i, k, n
        complex(sp) :: true_eigvals(test_size)
        real(sp) :: pi = 4.0_sp * atan(1.0_sp)

        return
    end subroutine test_evp_csp


    subroutine collect_eig_cdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("Eigs computation", test_evp_cdp) &
                    ]
        return
    end subroutine collect_eig_cdp_testsuite

    subroutine test_evp_cdp(error)
        !> Error type.
        type(error_type), allocatable, intent(out) :: error
        !> Test linear operator.
        type(linop_cdp), allocatable :: A
        !> Eigenvectors.
        type(vector_cdp), allocatable :: X(:)
        !> Eigenvalues.
        complex(dp), allocatable :: eigvals(:)
        !> Residuals.
        real(dp), allocatable :: residuals(:)
        !> Information flag.
        integer :: info
        !> Miscellaneous.
        integer :: i, k, n
        complex(dp) :: true_eigvals(test_size)
        real(dp) :: pi = 4.0_dp * atan(1.0_dp)

        return
    end subroutine test_evp_cdp




    !---------------------------------------------------------
    !-----     DEFINITION OF THE UNIT TESTS FOR SVDS     -----
    !---------------------------------------------------------

    subroutine collect_svd_rsp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("SVDS computation", test_svd_rsp) &
                    ]
        return
    end subroutine collect_svd_rsp_testsuite

    subroutine test_svd_rsp(error)
        !> Error type.
        type(error_type), allocatable, intent(out) :: error
        !> Test linear operator.
        type(linop_rsp), allocatable :: A
        !> Singular vectors.
        type(vector_rsp), allocatable :: U(:), V(:)
        !> Singular values.
        real(sp), allocatable :: S(:)
        !> Residuals.
        real(sp), allocatable :: residuals(:)
        !> Information flag.
        integer :: info
        !> Miscellaneous.
        integer :: i, k, n
        real(sp) :: true_svdvals(test_size)
        real(sp) :: pi = 4.0_sp * atan(1.0_sp)

        ! Allocate eigenvectors.
        allocate(U(test_size)) ; call initialize_krylov_subspace(U)
        allocate(V(test_size)) ; call initialize_krylov_subspace(V)
        
        ! Initialize linear operator with the Strang matrix.
        A = linop_rsp() ; A%data = 0.0_sp ; n = size(A%data, 1)

        do i = 1, n
            ! Diagonal entry.
            A%data(i, i) = 2.0_sp
            ! Upper diagonal entry.
            if (i < n) then
                A%data(i, i+1) = -1.0_sp
                A%data(i+1, i) = -1.0_sp
            endif
        enddo

        ! Compute spectral decomposition.
        call svds(A, U, S, V, residuals, info)

        ! Analytical singular values.
        do i = 1, test_size
            true_svdvals(i) = 2.0_sp * (1.0_sp + cos(i*pi/(test_size+1)))
        enddo

        call check(error, norm2(s - true_svdvals)**2 < rtol_sp)

        return
    end subroutine test_svd_rsp

    subroutine collect_svd_rdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("SVDS computation", test_svd_rdp) &
                    ]
        return
    end subroutine collect_svd_rdp_testsuite

    subroutine test_svd_rdp(error)
        !> Error type.
        type(error_type), allocatable, intent(out) :: error
        !> Test linear operator.
        type(linop_rdp), allocatable :: A
        !> Singular vectors.
        type(vector_rdp), allocatable :: U(:), V(:)
        !> Singular values.
        real(dp), allocatable :: S(:)
        !> Residuals.
        real(dp), allocatable :: residuals(:)
        !> Information flag.
        integer :: info
        !> Miscellaneous.
        integer :: i, k, n
        real(dp) :: true_svdvals(test_size)
        real(dp) :: pi = 4.0_dp * atan(1.0_dp)

        ! Allocate eigenvectors.
        allocate(U(test_size)) ; call initialize_krylov_subspace(U)
        allocate(V(test_size)) ; call initialize_krylov_subspace(V)
        
        ! Initialize linear operator with the Strang matrix.
        A = linop_rdp() ; A%data = 0.0_dp ; n = size(A%data, 1)

        do i = 1, n
            ! Diagonal entry.
            A%data(i, i) = 2.0_dp
            ! Upper diagonal entry.
            if (i < n) then
                A%data(i, i+1) = -1.0_dp
                A%data(i+1, i) = -1.0_dp
            endif
        enddo

        ! Compute spectral decomposition.
        call svds(A, U, S, V, residuals, info)

        ! Analytical singular values.
        do i = 1, test_size
            true_svdvals(i) = 2.0_dp * (1.0_dp + cos(i*pi/(test_size+1)))
        enddo

        call check(error, norm2(s - true_svdvals)**2 < rtol_dp)

        return
    end subroutine test_svd_rdp

    subroutine collect_svd_csp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("SVDS computation", test_svd_csp) &
                    ]
        return
    end subroutine collect_svd_csp_testsuite

    subroutine test_svd_csp(error)
        !> Error type.
        type(error_type), allocatable, intent(out) :: error
        !> Test linear operator.
        type(linop_csp), allocatable :: A
        !> Singular vectors.
        type(vector_csp), allocatable :: U(:), V(:)
        !> Singular values.
        real(sp), allocatable :: S(:)
        !> Residuals.
        real(sp), allocatable :: residuals(:)
        !> Information flag.
        integer :: info
        !> Miscellaneous.
        integer :: i, k, n
        real(sp) :: true_svdvals(test_size)
        real(sp) :: pi = 4.0_sp * atan(1.0_sp)

        return
    end subroutine test_svd_csp

    subroutine collect_svd_cdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("SVDS computation", test_svd_cdp) &
                    ]
        return
    end subroutine collect_svd_cdp_testsuite

    subroutine test_svd_cdp(error)
        !> Error type.
        type(error_type), allocatable, intent(out) :: error
        !> Test linear operator.
        type(linop_cdp), allocatable :: A
        !> Singular vectors.
        type(vector_cdp), allocatable :: U(:), V(:)
        !> Singular values.
        real(dp), allocatable :: S(:)
        !> Residuals.
        real(dp), allocatable :: residuals(:)
        !> Information flag.
        integer :: info
        !> Miscellaneous.
        integer :: i, k, n
        real(dp) :: true_svdvals(test_size)
        real(dp) :: pi = 4.0_dp * atan(1.0_dp)

        return
    end subroutine test_svd_cdp


    !----------------------------------------------------------
    !-----     DEFINITION OF THE UNIT TESTS FOR GMRES     -----
    !----------------------------------------------------------

    subroutine collect_gmres_rsp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [ &
                    new_unittest("Full GMRES", test_gmres_rsp), &
                    new_unittest("Full (SPD) GMRES", test_gmres_rsp) &
                    ]
        return
    end subroutine collect_gmres_rsp_testsuite

    subroutine test_gmres_rsp(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        !> Linear problem.
        type(linop_rsp) , allocatable :: A ! Linear Operator.
        type(vector_rsp), allocatable :: b ! Right-hand side vector.
        type(vector_rsp), allocatable :: x ! Solution vector.
        !> GMRES options.
        type(gmres_sp_opts) :: opts
        !> Information flag.
        integer :: info

        ! Initialize linear problem.
        A = linop_rsp()  ; call init_rand(A)
        b = vector_rsp() ; call init_rand(b)
        x = vector_rsp() ; call x%zero()

        ! GMRES solver.
        opts = gmres_sp_opts(kdim=test_size, verbose=.false.)
        call gmres(A, b, x, info, options=opts)

        ! Check convergence.
        call check(error, norm2(abs(matmul(A%data, x%data) - b%data)) < b%norm() * rtol_sp)

        return
    end subroutine test_gmres_rsp
    
    subroutine test_gmres_spd_rsp(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        !> Linear problem.
        type(spd_linop_rsp) , allocatable :: A ! Linear Operator.
        type(vector_rsp), allocatable :: b ! Right-hand side vector.
        type(vector_rsp), allocatable :: x ! Solution vector.
        !> GMRES options.
        type(gmres_sp_opts) :: opts
        !> Information flag.
        integer :: info

        ! Initialize linear problem.
        A = spd_linop_rsp()  ; call init_rand(A)
        b = vector_rsp() ; call init_rand(b)
        x = vector_rsp() ; call x%zero()

        ! GMRES solver.
        opts = gmres_sp_opts(kdim=test_size, verbose=.false.)
        call gmres(A, b, x, info, options=opts)

        ! Check convergence.
        call check(error, norm2(abs(matmul(A%data, x%data) - b%data)) < b%norm() * rtol_sp)

        return
    end subroutine test_gmres_spd_rsp

    subroutine collect_gmres_rdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [ &
                    new_unittest("Full GMRES", test_gmres_rdp), &
                    new_unittest("Full (SPD) GMRES", test_gmres_rdp) &
                    ]
        return
    end subroutine collect_gmres_rdp_testsuite

    subroutine test_gmres_rdp(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        !> Linear problem.
        type(linop_rdp) , allocatable :: A ! Linear Operator.
        type(vector_rdp), allocatable :: b ! Right-hand side vector.
        type(vector_rdp), allocatable :: x ! Solution vector.
        !> GMRES options.
        type(gmres_dp_opts) :: opts
        !> Information flag.
        integer :: info

        ! Initialize linear problem.
        A = linop_rdp()  ; call init_rand(A)
        b = vector_rdp() ; call init_rand(b)
        x = vector_rdp() ; call x%zero()

        ! GMRES solver.
        opts = gmres_dp_opts(kdim=test_size, verbose=.false.)
        call gmres(A, b, x, info, options=opts)

        ! Check convergence.
        call check(error, norm2(abs(matmul(A%data, x%data) - b%data)) < b%norm() * rtol_dp)

        return
    end subroutine test_gmres_rdp
    
    subroutine test_gmres_spd_rdp(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        !> Linear problem.
        type(spd_linop_rdp) , allocatable :: A ! Linear Operator.
        type(vector_rdp), allocatable :: b ! Right-hand side vector.
        type(vector_rdp), allocatable :: x ! Solution vector.
        !> GMRES options.
        type(gmres_dp_opts) :: opts
        !> Information flag.
        integer :: info

        ! Initialize linear problem.
        A = spd_linop_rdp()  ; call init_rand(A)
        b = vector_rdp() ; call init_rand(b)
        x = vector_rdp() ; call x%zero()

        ! GMRES solver.
        opts = gmres_dp_opts(kdim=test_size, verbose=.false.)
        call gmres(A, b, x, info, options=opts)

        ! Check convergence.
        call check(error, norm2(abs(matmul(A%data, x%data) - b%data)) < b%norm() * rtol_dp)

        return
    end subroutine test_gmres_spd_rdp

    subroutine collect_gmres_csp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [ &
                    new_unittest("Full GMRES", test_gmres_csp), &
                    new_unittest("Full (SPD) GMRES", test_gmres_csp) &
                    ]
        return
    end subroutine collect_gmres_csp_testsuite

    subroutine test_gmres_csp(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        !> Linear problem.
        type(linop_csp) , allocatable :: A ! Linear Operator.
        type(vector_csp), allocatable :: b ! Right-hand side vector.
        type(vector_csp), allocatable :: x ! Solution vector.
        !> GMRES options.
        type(gmres_sp_opts) :: opts
        !> Information flag.
        integer :: info

        ! Initialize linear problem.
        A = linop_csp()  ; call init_rand(A)
        b = vector_csp() ; call init_rand(b)
        x = vector_csp() ; call x%zero()

        ! GMRES solver.
        opts = gmres_sp_opts(kdim=test_size, verbose=.false.)
        call gmres(A, b, x, info, options=opts)

        ! Check convergence.
        call check(error, norm2(abs(matmul(A%data, x%data) - b%data)) < b%norm() * rtol_sp)

        return
    end subroutine test_gmres_csp
    
    subroutine test_gmres_spd_csp(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        !> Linear problem.
        type(hermitian_linop_csp), allocatable :: A
        type(vector_csp), allocatable :: b ! Right-hand side vector.
        type(vector_csp), allocatable :: x ! Solution vector.
        !> GMRES options.
        type(gmres_sp_opts) :: opts
        !> Information flag.
        integer :: info

        ! Initialize linear problem.
        A = hermitian_linop_csp() ; call init_rand(A)
        b = vector_csp() ; call init_rand(b)
        x = vector_csp() ; call x%zero()

        ! GMRES solver.
        opts = gmres_sp_opts(kdim=test_size, verbose=.false.)
        call gmres(A, b, x, info, options=opts)

        ! Check convergence.
        call check(error, norm2(abs(matmul(A%data, x%data) - b%data)) < b%norm() * rtol_sp)

        return
    end subroutine test_gmres_spd_csp

    subroutine collect_gmres_cdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [ &
                    new_unittest("Full GMRES", test_gmres_cdp), &
                    new_unittest("Full (SPD) GMRES", test_gmres_cdp) &
                    ]
        return
    end subroutine collect_gmres_cdp_testsuite

    subroutine test_gmres_cdp(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        !> Linear problem.
        type(linop_cdp) , allocatable :: A ! Linear Operator.
        type(vector_cdp), allocatable :: b ! Right-hand side vector.
        type(vector_cdp), allocatable :: x ! Solution vector.
        !> GMRES options.
        type(gmres_dp_opts) :: opts
        !> Information flag.
        integer :: info

        ! Initialize linear problem.
        A = linop_cdp()  ; call init_rand(A)
        b = vector_cdp() ; call init_rand(b)
        x = vector_cdp() ; call x%zero()

        ! GMRES solver.
        opts = gmres_dp_opts(kdim=test_size, verbose=.false.)
        call gmres(A, b, x, info, options=opts)

        ! Check convergence.
        call check(error, norm2(abs(matmul(A%data, x%data) - b%data)) < b%norm() * rtol_dp)

        return
    end subroutine test_gmres_cdp
    
    subroutine test_gmres_spd_cdp(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        !> Linear problem.
        type(hermitian_linop_cdp), allocatable :: A
        type(vector_cdp), allocatable :: b ! Right-hand side vector.
        type(vector_cdp), allocatable :: x ! Solution vector.
        !> GMRES options.
        type(gmres_dp_opts) :: opts
        !> Information flag.
        integer :: info

        ! Initialize linear problem.
        A = hermitian_linop_cdp() ; call init_rand(A)
        b = vector_cdp() ; call init_rand(b)
        x = vector_cdp() ; call x%zero()

        ! GMRES solver.
        opts = gmres_dp_opts(kdim=test_size, verbose=.false.)
        call gmres(A, b, x, info, options=opts)

        ! Check convergence.
        call check(error, norm2(abs(matmul(A%data, x%data) - b%data)) < b%norm() * rtol_dp)

        return
    end subroutine test_gmres_spd_cdp


    !-----------------------------------------------------------------------
    !-----     DEFINITION OF THE UNIT TESTS FOR CONJUGATE GRADIENT     -----
    !-----------------------------------------------------------------------

    subroutine collect_cg_rsp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [ &
                    new_unittest("Cong. Gradient", test_cg_rsp) &
                    ]
        return
    end subroutine collect_cg_rsp_testsuite

    subroutine test_cg_rsp(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        !> Linear problem.
        type(spd_linop_rsp), allocatable :: A
        type(vector_rsp), allocatable :: b
        type(vector_rsp), allocatable :: x
        type(cg_sp_opts) :: opts
        !> Information flag
        integer :: info, i

        ! Initialize linear problem.
        A = spd_linop_rsp() ; call init_rand(A)
        b = vector_rsp() ; call init_rand(b)
        x = vector_rsp() ; call x%zero()

        ! CG solver.
        opts = cg_sp_opts()
        call cg(A, b, x, info, options=opts)

        write(*, *) norm2(abs(matmul(A%data, x%data) - b%data)), b%norm() * rtol_sp 

        ! Check convergence.
        call check(error, norm2(abs(matmul(A%data, x%data) - b%data)) < b%norm() * rtol_sp)

        return
    end subroutine test_cg_rsp

    subroutine collect_cg_rdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [ &
                    new_unittest("Cong. Gradient", test_cg_rdp) &
                    ]
        return
    end subroutine collect_cg_rdp_testsuite

    subroutine test_cg_rdp(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        !> Linear problem.
        type(spd_linop_rdp), allocatable :: A
        type(vector_rdp), allocatable :: b
        type(vector_rdp), allocatable :: x
        type(cg_dp_opts) :: opts
        !> Information flag
        integer :: info, i

        ! Initialize linear problem.
        A = spd_linop_rdp() ; call init_rand(A)
        b = vector_rdp() ; call init_rand(b)
        x = vector_rdp() ; call x%zero()

        ! CG solver.
        opts = cg_dp_opts()
        call cg(A, b, x, info, options=opts)

        write(*, *) norm2(abs(matmul(A%data, x%data) - b%data)), b%norm() * rtol_dp 

        ! Check convergence.
        call check(error, norm2(abs(matmul(A%data, x%data) - b%data)) < b%norm() * rtol_dp)

        return
    end subroutine test_cg_rdp

    subroutine collect_cg_csp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [ &
                    new_unittest("Cong. Gradient", test_cg_csp) &
                    ]
        return
    end subroutine collect_cg_csp_testsuite

    subroutine test_cg_csp(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        !> Linear problem.
        type(hermitian_linop_csp), allocatable :: A
        type(vector_csp), allocatable :: b
        type(vector_csp), allocatable :: x
        type(cg_sp_opts) :: opts
        !> Information flag
        integer :: info, i

        ! Initialize linear problem.
        A = hermitian_linop_csp() ; call init_rand(A)
        b = vector_csp() ; call init_rand(b)
        x = vector_csp() ; call x%zero()

        ! CG solver.
        opts = cg_sp_opts()
        call cg(A, b, x, info, options=opts)

        write(*, *) norm2(abs(matmul(A%data, x%data) - b%data)), b%norm() * rtol_sp 

        ! Check convergence.
        call check(error, norm2(abs(matmul(A%data, x%data) - b%data)) < b%norm() * rtol_sp)

        return
    end subroutine test_cg_csp

    subroutine collect_cg_cdp_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [ &
                    new_unittest("Cong. Gradient", test_cg_cdp) &
                    ]
        return
    end subroutine collect_cg_cdp_testsuite

    subroutine test_cg_cdp(error)
        !> Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        !> Linear problem.
        type(hermitian_linop_cdp), allocatable :: A
        type(vector_cdp), allocatable :: b
        type(vector_cdp), allocatable :: x
        type(cg_dp_opts) :: opts
        !> Information flag
        integer :: info, i

        ! Initialize linear problem.
        A = hermitian_linop_cdp() ; call init_rand(A)
        b = vector_cdp() ; call init_rand(b)
        x = vector_cdp() ; call x%zero()

        ! CG solver.
        opts = cg_dp_opts()
        call cg(A, b, x, info, options=opts)

        write(*, *) norm2(abs(matmul(A%data, x%data) - b%data)), b%norm() * rtol_dp 

        ! Check convergence.
        call check(error, norm2(abs(matmul(A%data, x%data) - b%data)) < b%norm() * rtol_dp)

        return
    end subroutine test_cg_cdp


end module TestIterativeSolvers

