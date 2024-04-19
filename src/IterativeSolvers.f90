module lightkrylov_IterativeSolvers

    use iso_fortran_env, only: output_unit
    
    use stdlib_sorting, only: sort_index, int_size
    use stdlib_optval, only: optval
    use stdlib_io_npy, only: save_npy

    use lightkrylov_constants
    use lightkrylov_Utils
    use lightkrylov_AbstractVectors
    use lightkrylov_AbstractLinops
    use lightkrylov_BaseKrylov

    implicit none

    private

    public :: save_eigenspectrum
    public :: eigs

    interface save_eigenspectrum
        module procedure save_eigenspectrum_sp
        module procedure save_eigenspectrum_dp
    end interface

    interface eigs
        module procedure eigs_rsp
        module procedure eigs_rdp
        module procedure eigs_csp
        module procedure eigs_cdp
    end interface

contains

    !-------------------------------------
    !-----     UTILITY FUNCTIONS     -----
    !-------------------------------------

    elemental pure function compute_residual_rsp(beta, x) result(residual)
        !! Computes the residual associated with a Ritz eigenpair.
        real(sp), intent(in) :: beta
        !! Norm of the residual Krylov vector.
        real(sp), intent(in) :: x
        !! Last entry of the low-dimensional Ritz eigenvector.
        real(sp) :: residual
        !! Residual associated to the corresponding Ritz eigenpair.
        residual = abs(beta*x)
        return
    end function compute_residual_rsp

    elemental pure function compute_residual_rdp(beta, x) result(residual)
        !! Computes the residual associated with a Ritz eigenpair.
        real(dp), intent(in) :: beta
        !! Norm of the residual Krylov vector.
        real(dp), intent(in) :: x
        !! Last entry of the low-dimensional Ritz eigenvector.
        real(dp) :: residual
        !! Residual associated to the corresponding Ritz eigenpair.
        residual = abs(beta*x)
        return
    end function compute_residual_rdp

    elemental pure function compute_residual_csp(beta, x) result(residual)
        !! Computes the residual associated with a Ritz eigenpair.
        complex(sp), intent(in) :: beta
        !! Norm of the residual Krylov vector.
        complex(sp), intent(in) :: x
        !! Last entry of the low-dimensional Ritz eigenvector.
        real(sp) :: residual
        !! Residual associated to the corresponding Ritz eigenpair.
        residual = abs(beta*x)
        return
    end function compute_residual_csp

    elemental pure function compute_residual_cdp(beta, x) result(residual)
        !! Computes the residual associated with a Ritz eigenpair.
        complex(dp), intent(in) :: beta
        !! Norm of the residual Krylov vector.
        complex(dp), intent(in) :: x
        !! Last entry of the low-dimensional Ritz eigenvector.
        real(dp) :: residual
        !! Residual associated to the corresponding Ritz eigenpair.
        residual = abs(beta*x)
        return
    end function compute_residual_cdp


    subroutine save_eigenspectrum_sp(eigvals, residuals, fname)
        !! Saves the eigenspectrum and corresponding residuals to disk use the `npy` binary format.
        complex(sp), intent(in) :: eigvals(:)
        !! Eigenalues.
        real(sp), intent(in) :: residuals(:)
        !! Residual of the corresponding Ritz eigenpairs.
        character(len=*), intent(in) :: fname
        !! Name of the output file.

        ! Internal variables.
        real(sp) :: data(size(eigvals), 3)

        ! Store data.
        data(:, 1) = eigvals%re ; data(:, 2) = eigvals%im ; data(:, 3) = residuals
        ! Save the eigenspectrum to disk.
        call save_npy(fname, data)

        return
    end subroutine save_eigenspectrum_sp

    subroutine save_eigenspectrum_dp(eigvals, residuals, fname)
        !! Saves the eigenspectrum and corresponding residuals to disk use the `npy` binary format.
        complex(dp), intent(in) :: eigvals(:)
        !! Eigenalues.
        real(dp), intent(in) :: residuals(:)
        !! Residual of the corresponding Ritz eigenpairs.
        character(len=*), intent(in) :: fname
        !! Name of the output file.

        ! Internal variables.
        real(dp) :: data(size(eigvals), 3)

        ! Store data.
        data(:, 1) = eigvals%re ; data(:, 2) = eigvals%im ; data(:, 3) = residuals
        ! Save the eigenspectrum to disk.
        call save_npy(fname, data)

        return
    end subroutine save_eigenspectrum_dp


    !---------------------------------------------------
    !-----     GENERAL EIGENVALUE COMPUTATIONS     -----
    !---------------------------------------------------

    subroutine eigs_rsp(A, X, eigvals, residuals, info, kdim, tolerance, verbosity, transpose)
        class(abstract_linop_rsp), intent(in) :: A
        !! Linear operator whose leading eigenpairs need to be computed.
        class(abstract_vector_rsp), intent(out) :: X(:)
        !! Leading eigenvectors of \(\mathbf{A}\).
        complex(sp), allocatable, intent(out) :: eigvals(:)
        !! Leading eigenvalues of \(\mathbf{A}\).
        real(sp), allocatable, intent(out) :: residuals(:)
        !! Residuals associated to each Ritz eigenpair.
        integer, intent(out) :: info
        !! Information flag.
        integer, optional, intent(in) :: kdim
        !! Desired number of eigenpairs.
        real(sp), optional, intent(in) :: tolerance
        !! Tolerance.
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        logical, optional, intent(in) :: transpose
        !! Determine whether \(\mathbf{A}\) or \(\mathbf{A}^H\) is being used.

        !--------------------------------------
        !-----     Internal variables     -----
        !--------------------------------------

        !> Krylov subspace and Krylov subspace dimension.
        class(abstract_vector_rsp), allocatable :: Xwrk(:)
        integer :: kdim_
        !> Hessenberg matrix.
        real(sp), allocatable :: H(:, :)
        !> Working arrays for the eigenvectors and eigenvalues.
        real(sp), allocatable :: eigvecs_wrk(:, :)
        complex(sp), allocatable :: eigvals_wrk(:)
        real(sp), allocatable :: residuals_wrk(:)
        !> Miscellaneous.
        integer :: nev, conv
        integer :: i, j, k
        logical :: verbose, trans
        real(sp) :: tol
        real(sp) :: beta
        real(sp) :: alpha

        ! Deals with optional parameters.
        nev = size(X)
        kdim_   = optval(kdim, 4*nev)
        verbose = optval(verbosity, .false.)
        tol     = optval(tolerance, rtol_sp)

        ! Allocate eigenvalues.
        allocate(eigvals(nev)) ; eigvals = 0.0_sp

        ! Allocate working variables.
        allocate(Xwrk(kdim_+1), source=X(1)) ; call initialize_krylov_subspace(Xwrk) ; call Xwrk(1)%rand(.true.)
        allocate(H(kdim_+1, kdim_)) ; H = 0.0_sp
        allocate(eigvecs_wrk(kdim_, kdim_)) ; eigvecs_wrk = 0.0_sp
        allocate(eigvals_wrk(kdim_)) ; eigvals_wrk = 0.0_sp
        allocate(residuals_wrk(kdim_)) ; residuals_wrk = 0.0_sp

        ! Ritz eigenpairs computation.
        arnoldi_factorization: do k = 1, kdim_
            ! Arnoldi step.
            call arnoldi(A, Xwrk, H, info, kstart=k, kend=k, verbosity=verbose, transpose=transpose)
            if (info < 0) then
                call stop_error("eigs: Error in Arnoldi iteration.")
            endif

            ! Spectral decomposition of the k x k Hessenberg matrix.
            eigvals_wrk = 0.0_sp ; eigvecs_wrk = 0.0_sp
            call eig(H(1:k, 1:k), eigvecs_wrk(1:k, 1:k), eigvals_wrk(1:k))

            ! Compute residuals.
            beta = H(k+1, k)
            do i = 1, k
                if (eigvals_wrk(i)%im > 0) then
                    alpha = abs(cmplx(eigvecs_wrk(k, i), eigvecs_wrk(k, i+1), kind=sp))
                else if (eigvals_wrk(i)%im < 0) then
                    alpha = abs(cmplx(eigvecs_wrk(k, i-1), eigvecs_wrk(k, i), kind=sp))
                else
                    alpha = abs(eigvecs_wrk(k, i))
                endif
                residuals_wrk(i) = compute_residual_rsp(beta, alpha)
            enddo

            ! Check convergence.
            conv = count(residuals_wrk(1:k) < tol)
            if (conv >= nev) then
                exit arnoldi_factorization
            endif

        enddo arnoldi_factorization

        !--------------------------------
        !-----     POST-PROCESS     -----
        !--------------------------------

        block
        integer(int_size) :: indices(kdim_)
        real(sp) :: abs_eigvals(kdim_)
       
        ! Sort eigenvalues.
        abs_eigvals = abs(eigvals_wrk) ; call sort_index(abs_eigvals, indices, reverse=.true.)
        eigvals_wrk = eigvals_wrk(indices) ; eigvecs_wrk = eigvecs_wrk(:, indices)

        ! Store converged eigenvalues.
        eigvals = eigvals_wrk(1:nev)
        end block

        ! Construct eigenvectors.
        k = min(k, kdim_)
        do i = 1, nev
            call X(i)%zero()
            do j = 1, k
                call X(i)%axpby(1.0_sp, Xwrk(j), eigvecs_wrk(j, i))
            enddo
        enddo

        return
    end subroutine eigs_rsp
    subroutine eigs_rdp(A, X, eigvals, residuals, info, kdim, tolerance, verbosity, transpose)
        class(abstract_linop_rdp), intent(in) :: A
        !! Linear operator whose leading eigenpairs need to be computed.
        class(abstract_vector_rdp), intent(out) :: X(:)
        !! Leading eigenvectors of \(\mathbf{A}\).
        complex(dp), allocatable, intent(out) :: eigvals(:)
        !! Leading eigenvalues of \(\mathbf{A}\).
        real(dp), allocatable, intent(out) :: residuals(:)
        !! Residuals associated to each Ritz eigenpair.
        integer, intent(out) :: info
        !! Information flag.
        integer, optional, intent(in) :: kdim
        !! Desired number of eigenpairs.
        real(dp), optional, intent(in) :: tolerance
        !! Tolerance.
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        logical, optional, intent(in) :: transpose
        !! Determine whether \(\mathbf{A}\) or \(\mathbf{A}^H\) is being used.

        !--------------------------------------
        !-----     Internal variables     -----
        !--------------------------------------

        !> Krylov subspace and Krylov subspace dimension.
        class(abstract_vector_rdp), allocatable :: Xwrk(:)
        integer :: kdim_
        !> Hessenberg matrix.
        real(dp), allocatable :: H(:, :)
        !> Working arrays for the eigenvectors and eigenvalues.
        real(dp), allocatable :: eigvecs_wrk(:, :)
        complex(dp), allocatable :: eigvals_wrk(:)
        real(dp), allocatable :: residuals_wrk(:)
        !> Miscellaneous.
        integer :: nev, conv
        integer :: i, j, k
        logical :: verbose, trans
        real(dp) :: tol
        real(dp) :: beta
        real(dp) :: alpha

        ! Deals with optional parameters.
        nev = size(X)
        kdim_   = optval(kdim, 4*nev)
        verbose = optval(verbosity, .false.)
        tol     = optval(tolerance, rtol_dp)

        ! Allocate eigenvalues.
        allocate(eigvals(nev)) ; eigvals = 0.0_dp

        ! Allocate working variables.
        allocate(Xwrk(kdim_+1), source=X(1)) ; call initialize_krylov_subspace(Xwrk) ; call Xwrk(1)%rand(.true.)
        allocate(H(kdim_+1, kdim_)) ; H = 0.0_dp
        allocate(eigvecs_wrk(kdim_, kdim_)) ; eigvecs_wrk = 0.0_dp
        allocate(eigvals_wrk(kdim_)) ; eigvals_wrk = 0.0_dp
        allocate(residuals_wrk(kdim_)) ; residuals_wrk = 0.0_dp

        ! Ritz eigenpairs computation.
        arnoldi_factorization: do k = 1, kdim_
            ! Arnoldi step.
            call arnoldi(A, Xwrk, H, info, kstart=k, kend=k, verbosity=verbose, transpose=transpose)
            if (info < 0) then
                call stop_error("eigs: Error in Arnoldi iteration.")
            endif

            ! Spectral decomposition of the k x k Hessenberg matrix.
            eigvals_wrk = 0.0_dp ; eigvecs_wrk = 0.0_dp
            call eig(H(1:k, 1:k), eigvecs_wrk(1:k, 1:k), eigvals_wrk(1:k))

            ! Compute residuals.
            beta = H(k+1, k)
            do i = 1, k
                if (eigvals_wrk(i)%im > 0) then
                    alpha = abs(cmplx(eigvecs_wrk(k, i), eigvecs_wrk(k, i+1), kind=dp))
                else if (eigvals_wrk(i)%im < 0) then
                    alpha = abs(cmplx(eigvecs_wrk(k, i-1), eigvecs_wrk(k, i), kind=dp))
                else
                    alpha = abs(eigvecs_wrk(k, i))
                endif
                residuals_wrk(i) = compute_residual_rdp(beta, alpha)
            enddo

            ! Check convergence.
            conv = count(residuals_wrk(1:k) < tol)
            if (conv >= nev) then
                exit arnoldi_factorization
            endif

        enddo arnoldi_factorization

        !--------------------------------
        !-----     POST-PROCESS     -----
        !--------------------------------

        block
        integer(int_size) :: indices(kdim_)
        real(dp) :: abs_eigvals(kdim_)
       
        ! Sort eigenvalues.
        abs_eigvals = abs(eigvals_wrk) ; call sort_index(abs_eigvals, indices, reverse=.true.)
        eigvals_wrk = eigvals_wrk(indices) ; eigvecs_wrk = eigvecs_wrk(:, indices)

        ! Store converged eigenvalues.
        eigvals = eigvals_wrk(1:nev)
        end block

        ! Construct eigenvectors.
        k = min(k, kdim_)
        do i = 1, nev
            call X(i)%zero()
            do j = 1, k
                call X(i)%axpby(1.0_dp, Xwrk(j), eigvecs_wrk(j, i))
            enddo
        enddo

        return
    end subroutine eigs_rdp
    subroutine eigs_csp(A, X, eigvals, residuals, info, kdim, tolerance, verbosity, transpose)
        class(abstract_linop_csp), intent(in) :: A
        !! Linear operator whose leading eigenpairs need to be computed.
        class(abstract_vector_csp), intent(out) :: X(:)
        !! Leading eigenvectors of \(\mathbf{A}\).
        complex(sp), allocatable, intent(out) :: eigvals(:)
        !! Leading eigenvalues of \(\mathbf{A}\).
        real(sp), allocatable, intent(out) :: residuals(:)
        !! Residuals associated to each Ritz eigenpair.
        integer, intent(out) :: info
        !! Information flag.
        integer, optional, intent(in) :: kdim
        !! Desired number of eigenpairs.
        real(sp), optional, intent(in) :: tolerance
        !! Tolerance.
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        logical, optional, intent(in) :: transpose
        !! Determine whether \(\mathbf{A}\) or \(\mathbf{A}^H\) is being used.

        !--------------------------------------
        !-----     Internal variables     -----
        !--------------------------------------

        !> Krylov subspace and Krylov subspace dimension.
        class(abstract_vector_csp), allocatable :: Xwrk(:)
        integer :: kdim_
        !> Hessenberg matrix.
        complex(sp), allocatable :: H(:, :)
        !> Working arrays for the eigenvectors and eigenvalues.
        complex(sp), allocatable :: eigvecs_wrk(:, :)
        complex(sp), allocatable :: eigvals_wrk(:)
        real(sp), allocatable :: residuals_wrk(:)
        !> Miscellaneous.
        integer :: nev, conv
        integer :: i, j, k
        logical :: verbose, trans
        real(sp) :: tol
        complex(sp) :: beta

        ! Deals with optional parameters.
        nev = size(X)
        kdim_   = optval(kdim, 4*nev)
        verbose = optval(verbosity, .false.)
        tol     = optval(tolerance, rtol_sp)

        ! Allocate eigenvalues.
        allocate(eigvals(nev)) ; eigvals = 0.0_sp

        ! Allocate working variables.
        allocate(Xwrk(kdim_+1), source=X(1)) ; call initialize_krylov_subspace(Xwrk) ; call Xwrk(1)%rand(.true.)
        allocate(H(kdim_+1, kdim_)) ; H = 0.0_sp
        allocate(eigvecs_wrk(kdim_, kdim_)) ; eigvecs_wrk = 0.0_sp
        allocate(eigvals_wrk(kdim_)) ; eigvals_wrk = 0.0_sp
        allocate(residuals_wrk(kdim_)) ; residuals_wrk = 0.0_sp

        ! Ritz eigenpairs computation.
        arnoldi_factorization: do k = 1, kdim_
            ! Arnoldi step.
            call arnoldi(A, Xwrk, H, info, kstart=k, kend=k, verbosity=verbose, transpose=transpose)
            if (info < 0) then
                call stop_error("eigs: Error in Arnoldi iteration.")
            endif

            ! Spectral decomposition of the k x k Hessenberg matrix.
            eigvals_wrk = 0.0_sp ; eigvecs_wrk = 0.0_sp
            call eig(H(1:k, 1:k), eigvecs_wrk(1:k, 1:k), eigvals_wrk(1:k))

            ! Compute residuals.
            beta = H(k+1, k)
            residuals_wrk(1:k) = compute_residual_csp(beta, eigvecs_wrk(k,1:k))

            ! Check convergence.
            conv = count(residuals_wrk(1:k) < tol)
            if (conv >= nev) then
                exit arnoldi_factorization
            endif

        enddo arnoldi_factorization

        !--------------------------------
        !-----     POST-PROCESS     -----
        !--------------------------------

        block
        integer(int_size) :: indices(kdim_)
        real(sp) :: abs_eigvals(kdim_)
       
        ! Sort eigenvalues.
        abs_eigvals = abs(eigvals_wrk) ; call sort_index(abs_eigvals, indices, reverse=.true.)
        eigvals_wrk = eigvals_wrk(indices) ; eigvecs_wrk = eigvecs_wrk(:, indices)

        ! Store converged eigenvalues.
        eigvals = eigvals_wrk(1:nev)
        end block

        ! Construct eigenvectors.
        k = min(k, kdim_)
        do i = 1, nev
            call X(i)%zero()
            do j = 1, k
                call X(i)%axpby(cmplx(1.0_sp, 0.0_sp, kind=sp), Xwrk(j), eigvecs_wrk(j, i))
            enddo
        enddo

        return
    end subroutine eigs_csp
    subroutine eigs_cdp(A, X, eigvals, residuals, info, kdim, tolerance, verbosity, transpose)
        class(abstract_linop_cdp), intent(in) :: A
        !! Linear operator whose leading eigenpairs need to be computed.
        class(abstract_vector_cdp), intent(out) :: X(:)
        !! Leading eigenvectors of \(\mathbf{A}\).
        complex(dp), allocatable, intent(out) :: eigvals(:)
        !! Leading eigenvalues of \(\mathbf{A}\).
        real(dp), allocatable, intent(out) :: residuals(:)
        !! Residuals associated to each Ritz eigenpair.
        integer, intent(out) :: info
        !! Information flag.
        integer, optional, intent(in) :: kdim
        !! Desired number of eigenpairs.
        real(dp), optional, intent(in) :: tolerance
        !! Tolerance.
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        logical, optional, intent(in) :: transpose
        !! Determine whether \(\mathbf{A}\) or \(\mathbf{A}^H\) is being used.

        !--------------------------------------
        !-----     Internal variables     -----
        !--------------------------------------

        !> Krylov subspace and Krylov subspace dimension.
        class(abstract_vector_cdp), allocatable :: Xwrk(:)
        integer :: kdim_
        !> Hessenberg matrix.
        complex(dp), allocatable :: H(:, :)
        !> Working arrays for the eigenvectors and eigenvalues.
        complex(dp), allocatable :: eigvecs_wrk(:, :)
        complex(dp), allocatable :: eigvals_wrk(:)
        real(dp), allocatable :: residuals_wrk(:)
        !> Miscellaneous.
        integer :: nev, conv
        integer :: i, j, k
        logical :: verbose, trans
        real(dp) :: tol
        complex(dp) :: beta

        ! Deals with optional parameters.
        nev = size(X)
        kdim_   = optval(kdim, 4*nev)
        verbose = optval(verbosity, .false.)
        tol     = optval(tolerance, rtol_dp)

        ! Allocate eigenvalues.
        allocate(eigvals(nev)) ; eigvals = 0.0_dp

        ! Allocate working variables.
        allocate(Xwrk(kdim_+1), source=X(1)) ; call initialize_krylov_subspace(Xwrk) ; call Xwrk(1)%rand(.true.)
        allocate(H(kdim_+1, kdim_)) ; H = 0.0_dp
        allocate(eigvecs_wrk(kdim_, kdim_)) ; eigvecs_wrk = 0.0_dp
        allocate(eigvals_wrk(kdim_)) ; eigvals_wrk = 0.0_dp
        allocate(residuals_wrk(kdim_)) ; residuals_wrk = 0.0_dp

        ! Ritz eigenpairs computation.
        arnoldi_factorization: do k = 1, kdim_
            ! Arnoldi step.
            call arnoldi(A, Xwrk, H, info, kstart=k, kend=k, verbosity=verbose, transpose=transpose)
            if (info < 0) then
                call stop_error("eigs: Error in Arnoldi iteration.")
            endif

            ! Spectral decomposition of the k x k Hessenberg matrix.
            eigvals_wrk = 0.0_dp ; eigvecs_wrk = 0.0_dp
            call eig(H(1:k, 1:k), eigvecs_wrk(1:k, 1:k), eigvals_wrk(1:k))

            ! Compute residuals.
            beta = H(k+1, k)
            residuals_wrk(1:k) = compute_residual_cdp(beta, eigvecs_wrk(k,1:k))

            ! Check convergence.
            conv = count(residuals_wrk(1:k) < tol)
            if (conv >= nev) then
                exit arnoldi_factorization
            endif

        enddo arnoldi_factorization

        !--------------------------------
        !-----     POST-PROCESS     -----
        !--------------------------------

        block
        integer(int_size) :: indices(kdim_)
        real(dp) :: abs_eigvals(kdim_)
       
        ! Sort eigenvalues.
        abs_eigvals = abs(eigvals_wrk) ; call sort_index(abs_eigvals, indices, reverse=.true.)
        eigvals_wrk = eigvals_wrk(indices) ; eigvecs_wrk = eigvecs_wrk(:, indices)

        ! Store converged eigenvalues.
        eigvals = eigvals_wrk(1:nev)
        end block

        ! Construct eigenvectors.
        k = min(k, kdim_)
        do i = 1, nev
            call X(i)%zero()
            do j = 1, k
                call X(i)%axpby(cmplx(1.0_dp, 0.0_dp, kind=dp), Xwrk(j), eigvecs_wrk(j, i))
            enddo
        enddo

        return
    end subroutine eigs_cdp

end module lightkrylov_IterativeSolvers
