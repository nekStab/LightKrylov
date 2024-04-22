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
    public :: svds
    public :: gmres

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

    interface svds
        module procedure svds_rsp
        module procedure svds_rdp
        module procedure svds_csp
        module procedure svds_cdp
    end interface

    interface gmres
        module procedure gmres_rsp
        module procedure gmres_rdp
        module procedure gmres_csp
        module procedure gmres_cdp
    end interface

    !------------------------------------------------------
    !-----     ABSTRACT PRECONDITIONER DEFINITION     -----
    !------------------------------------------------------

    type, abstract, public :: abstract_precond_rsp
    contains
        private
        procedure(abstract_apply_rsp), pass(self), deferred :: apply
    end type

    abstract interface
        subroutine abstract_apply_rsp(self, vec)
            !! Abstract interface to apply a preconditioner in `LightKrylov`.
            import abstract_precond_rsp, abstract_vector_rsp
            class(abstract_precond_rsp), intent(in) :: self
            !! Preconditioner.
            class(abstract_vector_rsp), intent(inout) :: vec
            !! Input/Output vector.
        end subroutine abstract_apply_rsp
    end interface
    
    type, abstract, public :: abstract_precond_rdp
    contains
        private
        procedure(abstract_apply_rdp), pass(self), deferred :: apply
    end type

    abstract interface
        subroutine abstract_apply_rdp(self, vec)
            !! Abstract interface to apply a preconditioner in `LightKrylov`.
            import abstract_precond_rdp, abstract_vector_rdp
            class(abstract_precond_rdp), intent(in) :: self
            !! Preconditioner.
            class(abstract_vector_rdp), intent(inout) :: vec
            !! Input/Output vector.
        end subroutine abstract_apply_rdp
    end interface
    
    type, abstract, public :: abstract_precond_csp
    contains
        private
        procedure(abstract_apply_csp), pass(self), deferred :: apply
    end type

    abstract interface
        subroutine abstract_apply_csp(self, vec)
            !! Abstract interface to apply a preconditioner in `LightKrylov`.
            import abstract_precond_csp, abstract_vector_csp
            class(abstract_precond_csp), intent(in) :: self
            !! Preconditioner.
            class(abstract_vector_csp), intent(inout) :: vec
            !! Input/Output vector.
        end subroutine abstract_apply_csp
    end interface
    
    type, abstract, public :: abstract_precond_cdp
    contains
        private
        procedure(abstract_apply_cdp), pass(self), deferred :: apply
    end type

    abstract interface
        subroutine abstract_apply_cdp(self, vec)
            !! Abstract interface to apply a preconditioner in `LightKrylov`.
            import abstract_precond_cdp, abstract_vector_cdp
            class(abstract_precond_cdp), intent(in) :: self
            !! Preconditioner.
            class(abstract_vector_cdp), intent(inout) :: vec
            !! Input/Output vector.
        end subroutine abstract_apply_cdp
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
        H = 0.0_sp
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
        H = 0.0_dp
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
        H = 0.0_sp
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
        H = 0.0_dp
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


    !------------------------------------------------
    !-----     SINGULAR VALUE DECOMPOSITION     -----
    !------------------------------------------------

    subroutine svds_rsp(A, U, S, V, residuals, info, kdim, tolerance, verbosity)
        class(abstract_linop_rsp), intent(in) :: A
        !! Linear operator whose leading singular triplets need to be computed.
        class(abstract_vector_rsp), intent(out) :: U(:)
        !! Leading left singular vectors.
        real(sp), allocatable, intent(out) :: S(:)
        !! Leading singular values.
        class(abstract_vector_rsp), intent(out) :: V(:)
        !! Leading right singular vectors.
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

        !--------------------------------------
        !-----     Internal variables     -----
        !--------------------------------------
        ! Left and right Krylov subspaces.
        integer :: kdim_
        class(abstract_vector_rsp), allocatable :: Uwrk(:), Vwrk(:)
        ! Bidiagonal matrix.
        real(sp), allocatable :: B(:, :)
        ! Working arrays for the singular vectors and singular values.
        real(sp), allocatable :: svdvals_wrk(:)
        real(sp), allocatable :: umat(:, :), vmat(:, :)
        real(sp), allocatable :: residuals_wrk(:)
        ! Miscellaneous.
        integer :: nsv, conv
        integer :: i, j, k
        logical :: verbose
        real(sp) :: tol
        real(sp) :: beta, alpha

        ! Deals with the optional arguments.
        nsv = size(U)
        kdim_ = optval(kdim, 4*nsv)
        verbose = optval(verbosity, .false.)
        tol     = optval(tolerance, rtol_sp)

        ! Allocate working variables.
        allocate(Uwrk(kdim_+1), source=U(1)) ; call initialize_krylov_subspace(Uwrk) ; call Uwrk(1)%rand(.true.)
        allocate(Vwrk(kdim_+1), source=V(1)) ; call initialize_krylov_subspace(Vwrk)
        allocate(svdvals_wrk(kdim_)) ; svdvals_wrk = 0.0_sp
        allocate(umat(kdim_, kdim_)) ; umat = 0.0_sp
        allocate(vmat(kdim_, kdim_)) ; vmat = 0.0_sp
        allocate(residuals_wrk(kdim_)) ; residuals_wrk = 0.0_sp
        allocate(B(kdim_+1, kdim_)) ; B = 0.0_sp

        ! Ritz singular triplets computation.
        lanczos : do k = 1, kdim_
            ! Lanczos bidiag. step.
            call lanczos_bidiagonalization(A, Uwrk, Vwrk, B, info, kstart=k, kend=k, verbosity=verbosity, tol=tol)
            if (info < 0) then
                call stop_error("svds: Error in Lanczos bidiagonalization.")
            endif

            ! SVD of the k x k bidiagonal matrix.
            svdvals_wrk = 0.0_sp ; umat = 0.0_sp ; vmat = 0.0_sp
            call svd(B(1:k, 1:k), umat(1:k, 1:k), svdvals_wrk(1:k), vmat(1:k, 1:k))

            ! Compute residuals.
            beta = B(k+1, k)
            residuals_wrk(1:k) = compute_residual_rsp(beta, vmat(k, 1:k))

            ! Check for convergence.
            conv = count(residuals_wrk(1:k) < tol)
            if (conv >= nsv) then
                exit lanczos
            endif
        enddo lanczos

        !--------------------------------
        !-----     POST-PROCESS     -----
        !--------------------------------

        ! Singular values.
        S = svdvals_wrk(1:nsv)

        ! Singular vectors.
        k = min(k, kdim_)
        do i = 1, nsv
            call U(i)%zero() ; call V(i)%zero()
            do j = 1, k
                call U(i)%axpby(1.0_sp, Uwrk(j), umat(j, i))
            enddo
        enddo

        return
    end subroutine svds_rsp

    subroutine svds_rdp(A, U, S, V, residuals, info, kdim, tolerance, verbosity)
        class(abstract_linop_rdp), intent(in) :: A
        !! Linear operator whose leading singular triplets need to be computed.
        class(abstract_vector_rdp), intent(out) :: U(:)
        !! Leading left singular vectors.
        real(dp), allocatable, intent(out) :: S(:)
        !! Leading singular values.
        class(abstract_vector_rdp), intent(out) :: V(:)
        !! Leading right singular vectors.
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

        !--------------------------------------
        !-----     Internal variables     -----
        !--------------------------------------
        ! Left and right Krylov subspaces.
        integer :: kdim_
        class(abstract_vector_rdp), allocatable :: Uwrk(:), Vwrk(:)
        ! Bidiagonal matrix.
        real(dp), allocatable :: B(:, :)
        ! Working arrays for the singular vectors and singular values.
        real(dp), allocatable :: svdvals_wrk(:)
        real(dp), allocatable :: umat(:, :), vmat(:, :)
        real(dp), allocatable :: residuals_wrk(:)
        ! Miscellaneous.
        integer :: nsv, conv
        integer :: i, j, k
        logical :: verbose
        real(dp) :: tol
        real(dp) :: beta, alpha

        ! Deals with the optional arguments.
        nsv = size(U)
        kdim_ = optval(kdim, 4*nsv)
        verbose = optval(verbosity, .false.)
        tol     = optval(tolerance, rtol_dp)

        ! Allocate working variables.
        allocate(Uwrk(kdim_+1), source=U(1)) ; call initialize_krylov_subspace(Uwrk) ; call Uwrk(1)%rand(.true.)
        allocate(Vwrk(kdim_+1), source=V(1)) ; call initialize_krylov_subspace(Vwrk)
        allocate(svdvals_wrk(kdim_)) ; svdvals_wrk = 0.0_dp
        allocate(umat(kdim_, kdim_)) ; umat = 0.0_dp
        allocate(vmat(kdim_, kdim_)) ; vmat = 0.0_dp
        allocate(residuals_wrk(kdim_)) ; residuals_wrk = 0.0_dp
        allocate(B(kdim_+1, kdim_)) ; B = 0.0_dp

        ! Ritz singular triplets computation.
        lanczos : do k = 1, kdim_
            ! Lanczos bidiag. step.
            call lanczos_bidiagonalization(A, Uwrk, Vwrk, B, info, kstart=k, kend=k, verbosity=verbosity, tol=tol)
            if (info < 0) then
                call stop_error("svds: Error in Lanczos bidiagonalization.")
            endif

            ! SVD of the k x k bidiagonal matrix.
            svdvals_wrk = 0.0_dp ; umat = 0.0_dp ; vmat = 0.0_dp
            call svd(B(1:k, 1:k), umat(1:k, 1:k), svdvals_wrk(1:k), vmat(1:k, 1:k))

            ! Compute residuals.
            beta = B(k+1, k)
            residuals_wrk(1:k) = compute_residual_rdp(beta, vmat(k, 1:k))

            ! Check for convergence.
            conv = count(residuals_wrk(1:k) < tol)
            if (conv >= nsv) then
                exit lanczos
            endif
        enddo lanczos

        !--------------------------------
        !-----     POST-PROCESS     -----
        !--------------------------------

        ! Singular values.
        S = svdvals_wrk(1:nsv)

        ! Singular vectors.
        k = min(k, kdim_)
        do i = 1, nsv
            call U(i)%zero() ; call V(i)%zero()
            do j = 1, k
                call U(i)%axpby(1.0_dp, Uwrk(j), umat(j, i))
            enddo
        enddo

        return
    end subroutine svds_rdp

    subroutine svds_csp(A, U, S, V, residuals, info, kdim, tolerance, verbosity)
        class(abstract_linop_csp), intent(in) :: A
        !! Linear operator whose leading singular triplets need to be computed.
        class(abstract_vector_csp), intent(out) :: U(:)
        !! Leading left singular vectors.
        real(sp), allocatable, intent(out) :: S(:)
        !! Leading singular values.
        class(abstract_vector_csp), intent(out) :: V(:)
        !! Leading right singular vectors.
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

        !--------------------------------------
        !-----     Internal variables     -----
        !--------------------------------------
        ! Left and right Krylov subspaces.
        integer :: kdim_
        class(abstract_vector_csp), allocatable :: Uwrk(:), Vwrk(:)
        ! Bidiagonal matrix.
        complex(sp), allocatable :: B(:, :)
        ! Working arrays for the singular vectors and singular values.
        real(sp), allocatable :: svdvals_wrk(:)
        complex(sp), allocatable :: umat(:, :), vmat(:, :)
        real(sp), allocatable :: residuals_wrk(:)
        ! Miscellaneous.
        integer :: nsv, conv
        integer :: i, j, k
        logical :: verbose
        real(sp) :: tol
        complex(sp) :: beta, alpha

        ! Deals with the optional arguments.
        nsv = size(U)
        kdim_ = optval(kdim, 4*nsv)
        verbose = optval(verbosity, .false.)
        tol     = optval(tolerance, rtol_sp)

        ! Allocate working variables.
        allocate(Uwrk(kdim_+1), source=U(1)) ; call initialize_krylov_subspace(Uwrk) ; call Uwrk(1)%rand(.true.)
        allocate(Vwrk(kdim_+1), source=V(1)) ; call initialize_krylov_subspace(Vwrk)
        allocate(svdvals_wrk(kdim_)) ; svdvals_wrk = 0.0_sp
        allocate(umat(kdim_, kdim_)) ; umat = 0.0_sp
        allocate(vmat(kdim_, kdim_)) ; vmat = 0.0_sp
        allocate(residuals_wrk(kdim_)) ; residuals_wrk = 0.0_sp
        allocate(B(kdim_+1, kdim_)) ; B = 0.0_sp

        ! Ritz singular triplets computation.
        lanczos : do k = 1, kdim_
            ! Lanczos bidiag. step.
            call lanczos_bidiagonalization(A, Uwrk, Vwrk, B, info, kstart=k, kend=k, verbosity=verbosity, tol=tol)
            if (info < 0) then
                call stop_error("svds: Error in Lanczos bidiagonalization.")
            endif

            ! SVD of the k x k bidiagonal matrix.
            svdvals_wrk = 0.0_sp ; umat = 0.0_sp ; vmat = 0.0_sp
            call svd(B(1:k, 1:k), umat(1:k, 1:k), svdvals_wrk(1:k), vmat(1:k, 1:k))

            ! Compute residuals.
            beta = B(k+1, k)
            residuals_wrk(1:k) = compute_residual_csp(beta, vmat(k, 1:k))

            ! Check for convergence.
            conv = count(residuals_wrk(1:k) < tol)
            if (conv >= nsv) then
                exit lanczos
            endif
        enddo lanczos

        !--------------------------------
        !-----     POST-PROCESS     -----
        !--------------------------------

        ! Singular values.
        S = svdvals_wrk(1:nsv)

        ! Singular vectors.
        k = min(k, kdim_)
        do i = 1, nsv
            call U(i)%zero() ; call V(i)%zero()
            do j = 1, k
                call U(i)%axpby(cmplx(1.0_sp, 0.0_sp, kind=sp), Uwrk(j), umat(j, i))
                call V(i)%axpby(cmplx(1.0_sp, 0.0_sp, kind=sp), Vwrk(j), Vmat(j, i))
            enddo
        enddo

        return
    end subroutine svds_csp

    subroutine svds_cdp(A, U, S, V, residuals, info, kdim, tolerance, verbosity)
        class(abstract_linop_cdp), intent(in) :: A
        !! Linear operator whose leading singular triplets need to be computed.
        class(abstract_vector_cdp), intent(out) :: U(:)
        !! Leading left singular vectors.
        real(dp), allocatable, intent(out) :: S(:)
        !! Leading singular values.
        class(abstract_vector_cdp), intent(out) :: V(:)
        !! Leading right singular vectors.
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

        !--------------------------------------
        !-----     Internal variables     -----
        !--------------------------------------
        ! Left and right Krylov subspaces.
        integer :: kdim_
        class(abstract_vector_cdp), allocatable :: Uwrk(:), Vwrk(:)
        ! Bidiagonal matrix.
        complex(dp), allocatable :: B(:, :)
        ! Working arrays for the singular vectors and singular values.
        real(dp), allocatable :: svdvals_wrk(:)
        complex(dp), allocatable :: umat(:, :), vmat(:, :)
        real(dp), allocatable :: residuals_wrk(:)
        ! Miscellaneous.
        integer :: nsv, conv
        integer :: i, j, k
        logical :: verbose
        real(dp) :: tol
        complex(dp) :: beta, alpha

        ! Deals with the optional arguments.
        nsv = size(U)
        kdim_ = optval(kdim, 4*nsv)
        verbose = optval(verbosity, .false.)
        tol     = optval(tolerance, rtol_dp)

        ! Allocate working variables.
        allocate(Uwrk(kdim_+1), source=U(1)) ; call initialize_krylov_subspace(Uwrk) ; call Uwrk(1)%rand(.true.)
        allocate(Vwrk(kdim_+1), source=V(1)) ; call initialize_krylov_subspace(Vwrk)
        allocate(svdvals_wrk(kdim_)) ; svdvals_wrk = 0.0_dp
        allocate(umat(kdim_, kdim_)) ; umat = 0.0_dp
        allocate(vmat(kdim_, kdim_)) ; vmat = 0.0_dp
        allocate(residuals_wrk(kdim_)) ; residuals_wrk = 0.0_dp
        allocate(B(kdim_+1, kdim_)) ; B = 0.0_dp

        ! Ritz singular triplets computation.
        lanczos : do k = 1, kdim_
            ! Lanczos bidiag. step.
            call lanczos_bidiagonalization(A, Uwrk, Vwrk, B, info, kstart=k, kend=k, verbosity=verbosity, tol=tol)
            if (info < 0) then
                call stop_error("svds: Error in Lanczos bidiagonalization.")
            endif

            ! SVD of the k x k bidiagonal matrix.
            svdvals_wrk = 0.0_dp ; umat = 0.0_dp ; vmat = 0.0_dp
            call svd(B(1:k, 1:k), umat(1:k, 1:k), svdvals_wrk(1:k), vmat(1:k, 1:k))

            ! Compute residuals.
            beta = B(k+1, k)
            residuals_wrk(1:k) = compute_residual_cdp(beta, vmat(k, 1:k))

            ! Check for convergence.
            conv = count(residuals_wrk(1:k) < tol)
            if (conv >= nsv) then
                exit lanczos
            endif
        enddo lanczos

        !--------------------------------
        !-----     POST-PROCESS     -----
        !--------------------------------

        ! Singular values.
        S = svdvals_wrk(1:nsv)

        ! Singular vectors.
        k = min(k, kdim_)
        do i = 1, nsv
            call U(i)%zero() ; call V(i)%zero()
            do j = 1, k
                call U(i)%axpby(cmplx(1.0_dp, 0.0_dp, kind=dp), Uwrk(j), umat(j, i))
                call V(i)%axpby(cmplx(1.0_dp, 0.0_dp, kind=dp), Vwrk(j), Vmat(j, i))
            enddo
        enddo

        return
    end subroutine svds_cdp


    !-------------------------------------------------------
    !-----     GENERALIZED MINIMUM RESIDUAL METHOD     -----
    !-------------------------------------------------------

    subroutine gmres_rsp(A, b, x, info, preconditioner, options, transpose)
        class(abstract_linop_rsp), intent(in) :: A
        !! Linear operator to be inverted.
        class(abstract_vector_rsp), intent(in) :: b
        !! Right-hand side vector.
        class(abstract_vector_rsp), intent(out) :: x
        !! Solution vector.
        integer, intent(out) :: info
        !! Information flag.
        class(abstract_precond_rsp), optional, intent(in) :: preconditioner
        !! Preconditioner (optional).
        type(gmres_sp_opts), optional, intent(in) :: options
        !! GMRES options.   
        logical, optional, intent(in) :: transpose
        !! Whether \(\mathbf{A}\) or \(\mathbf{A}^H\) is being used.

        !--------------------------------------
        !-----     Internal variables     -----
        !--------------------------------------

        ! Options.
        integer :: kdim, maxiter
        real(sp) :: tol
        logical :: verbose, trans
        type(gmres_sp_opts) :: opts

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
        integer :: i, j, k
        real(sp) :: alpha
        class(abstract_vector_rsp), allocatable :: dx, wrk

        ! Deals with the optional args.
        if (present(options)) then
            opts = gmres_sp_opts( &
                        kdim    = options%kdim, &
                        maxiter = options%maxiter, &
                        atol    = options%atol, &
                        rtol    = options%rtol, &
                        verbose = options%verbose &
                    )
        else
            opts = gmres_sp_opts()
        endif

        kdim = opts%kdim ; maxiter = opts%maxiter
        tol = opts%atol + opts%rtol * b%norm() ; verbose = opts%verbose
        trans = optval(transpose, .false.)

        ! Deals with the preconditioner.
        if (present(preconditioner)) then
            has_precond = .true.
            allocate(precond, source=preconditioner)
        else
            has_precond = .false.
        endif

        ! Initialize working variables.
        allocate(wrk, source=b) ; call wrk%zero()
        allocate(V(1:kdim+1), source=b) ; call initialize_krylov_subspace(V)
        allocate(H(kdim+1, kdim)) ; H = 0.0_sp
        allocate(y(kdim)) ; y = 0.0_sp
        allocate(e(kdim+1)) ; e = 0.0_sp

        info = 0

        ! Initial Krylov vector.
        if (x%norm() > 0) then
            if (trans) then
                call A%rmatvec(x, V(1))
            else
                call A%matvec(x, V(1))
            endif
        endif

        call V(1)%sub(b) ; call V(1)%chsgn()
        beta = V(1)%norm() ; call V(1)%scal(1.0_sp/beta)

        ! Iterative solver.
        gmres_iter : do i = 1, maxiter
            ! Zero-out variables.
            H = 0.0_sp ; y = 0.0_sp ; e = 0.0_sp ; e(1) = beta
            call initialize_krylov_subspace(V(2:kdim+1))

            ! Arnoldi factorization.
            arnoldi_fact: do k = 1, kdim
                ! Preconditioner.
                wrk = V(k) ; if (has_precond) call precond%apply(wrk)

                ! Matrix-vector product.
                if (trans) then
                    call A%rmatvec(wrk, V(k+1))
                else
                    call A%matvec(wrk, V(k+1))
                endif

                ! Gram-Schmid orthogonalization (twice is enough).
                do j = 1, k
                    alpha = V(j)%dot(V(k+1))
                    call V(k+1)%axpby(one_rsp, V(j), -alpha)
                    H(j, k) = alpha
                enddo
                do j = 1, k
                    alpha = V(j)%dot(V(k+1))
                    call V(k+1)%axpby(one_rsp, V(j), -alpha)
                    H(j, k) = H(j, k) + alpha
                enddo

                ! Update Hessenberg matrix and normalize residual Krylov vector.
                H(k+1, k) = V(k+1)%norm()
                if (abs(H(k+1, k)) > tol) then
                    call V(k+1)%scal(1.0_sp / H(k+1, k))
                endif

                ! Least-squares problem.
                call lstsq(H(1:k+1, 1:k), e(1:k+1), y(1:k))

                ! Compute residual.
                beta = norm2(abs(e(1:k+1) - matmul(H(1:k+1, 1:k), y(1:k))))

                ! Current number of iterations performed.
                info = info + 1

                ! Check convergence.
                if (abs(beta) <= tol) then
                    exit arnoldi_fact
                endif
            enddo arnoldi_fact

            ! Update solution.
            k = min(k, kdim)
            if (allocated(dx) .eqv. .false.) allocate(dx, source=x); call dx%zero()
            do j = 1, k
                call dx%axpby(one_rsp, V(j), y(j))
            enddo
            if (has_precond) call precond%apply(dx)
            call x%add(dx)

            ! Recompute residual for sanity check.
            if (trans) then
                call A%rmatvec(x, v(1))
            else
                call A%matvec(x, v(1))
            endif
            call v(1)%sub(b) ; call v(1)%chsgn()

            ! Initialize new starting Krylov vector if needed.
            beta = v(1)%norm() ; call v(1)%scal(1.0_sp / beta)

            ! Exit gmres if desired accuracy is reached.
            if (abs(beta) <= tol) exit gmres_iter

        enddo gmres_iter

        return
    end subroutine gmres_rsp

    subroutine gmres_rdp(A, b, x, info, preconditioner, options, transpose)
        class(abstract_linop_rdp), intent(in) :: A
        !! Linear operator to be inverted.
        class(abstract_vector_rdp), intent(in) :: b
        !! Right-hand side vector.
        class(abstract_vector_rdp), intent(out) :: x
        !! Solution vector.
        integer, intent(out) :: info
        !! Information flag.
        class(abstract_precond_rdp), optional, intent(in) :: preconditioner
        !! Preconditioner (optional).
        type(gmres_dp_opts), optional, intent(in) :: options
        !! GMRES options.   
        logical, optional, intent(in) :: transpose
        !! Whether \(\mathbf{A}\) or \(\mathbf{A}^H\) is being used.

        !--------------------------------------
        !-----     Internal variables     -----
        !--------------------------------------

        ! Options.
        integer :: kdim, maxiter
        real(dp) :: tol
        logical :: verbose, trans
        type(gmres_dp_opts) :: opts

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
        integer :: i, j, k
        real(dp) :: alpha
        class(abstract_vector_rdp), allocatable :: dx, wrk

        ! Deals with the optional args.
        if (present(options)) then
            opts = gmres_dp_opts( &
                        kdim    = options%kdim, &
                        maxiter = options%maxiter, &
                        atol    = options%atol, &
                        rtol    = options%rtol, &
                        verbose = options%verbose &
                    )
        else
            opts = gmres_dp_opts()
        endif

        kdim = opts%kdim ; maxiter = opts%maxiter
        tol = opts%atol + opts%rtol * b%norm() ; verbose = opts%verbose
        trans = optval(transpose, .false.)

        ! Deals with the preconditioner.
        if (present(preconditioner)) then
            has_precond = .true.
            allocate(precond, source=preconditioner)
        else
            has_precond = .false.
        endif

        ! Initialize working variables.
        allocate(wrk, source=b) ; call wrk%zero()
        allocate(V(1:kdim+1), source=b) ; call initialize_krylov_subspace(V)
        allocate(H(kdim+1, kdim)) ; H = 0.0_dp
        allocate(y(kdim)) ; y = 0.0_dp
        allocate(e(kdim+1)) ; e = 0.0_dp

        info = 0

        ! Initial Krylov vector.
        if (x%norm() > 0) then
            if (trans) then
                call A%rmatvec(x, V(1))
            else
                call A%matvec(x, V(1))
            endif
        endif

        call V(1)%sub(b) ; call V(1)%chsgn()
        beta = V(1)%norm() ; call V(1)%scal(1.0_dp/beta)

        ! Iterative solver.
        gmres_iter : do i = 1, maxiter
            ! Zero-out variables.
            H = 0.0_dp ; y = 0.0_dp ; e = 0.0_dp ; e(1) = beta
            call initialize_krylov_subspace(V(2:kdim+1))

            ! Arnoldi factorization.
            arnoldi_fact: do k = 1, kdim
                ! Preconditioner.
                wrk = V(k) ; if (has_precond) call precond%apply(wrk)

                ! Matrix-vector product.
                if (trans) then
                    call A%rmatvec(wrk, V(k+1))
                else
                    call A%matvec(wrk, V(k+1))
                endif

                ! Gram-Schmid orthogonalization (twice is enough).
                do j = 1, k
                    alpha = V(j)%dot(V(k+1))
                    call V(k+1)%axpby(one_rdp, V(j), -alpha)
                    H(j, k) = alpha
                enddo
                do j = 1, k
                    alpha = V(j)%dot(V(k+1))
                    call V(k+1)%axpby(one_rdp, V(j), -alpha)
                    H(j, k) = H(j, k) + alpha
                enddo

                ! Update Hessenberg matrix and normalize residual Krylov vector.
                H(k+1, k) = V(k+1)%norm()
                if (abs(H(k+1, k)) > tol) then
                    call V(k+1)%scal(1.0_dp / H(k+1, k))
                endif

                ! Least-squares problem.
                call lstsq(H(1:k+1, 1:k), e(1:k+1), y(1:k))

                ! Compute residual.
                beta = norm2(abs(e(1:k+1) - matmul(H(1:k+1, 1:k), y(1:k))))

                ! Current number of iterations performed.
                info = info + 1

                ! Check convergence.
                if (abs(beta) <= tol) then
                    exit arnoldi_fact
                endif
            enddo arnoldi_fact

            ! Update solution.
            k = min(k, kdim)
            if (allocated(dx) .eqv. .false.) allocate(dx, source=x); call dx%zero()
            do j = 1, k
                call dx%axpby(one_rdp, V(j), y(j))
            enddo
            if (has_precond) call precond%apply(dx)
            call x%add(dx)

            ! Recompute residual for sanity check.
            if (trans) then
                call A%rmatvec(x, v(1))
            else
                call A%matvec(x, v(1))
            endif
            call v(1)%sub(b) ; call v(1)%chsgn()

            ! Initialize new starting Krylov vector if needed.
            beta = v(1)%norm() ; call v(1)%scal(1.0_dp / beta)

            ! Exit gmres if desired accuracy is reached.
            if (abs(beta) <= tol) exit gmres_iter

        enddo gmres_iter

        return
    end subroutine gmres_rdp

    subroutine gmres_csp(A, b, x, info, preconditioner, options, transpose)
        class(abstract_linop_csp), intent(in) :: A
        !! Linear operator to be inverted.
        class(abstract_vector_csp), intent(in) :: b
        !! Right-hand side vector.
        class(abstract_vector_csp), intent(out) :: x
        !! Solution vector.
        integer, intent(out) :: info
        !! Information flag.
        class(abstract_precond_csp), optional, intent(in) :: preconditioner
        !! Preconditioner (optional).
        type(gmres_sp_opts), optional, intent(in) :: options
        !! GMRES options.   
        logical, optional, intent(in) :: transpose
        !! Whether \(\mathbf{A}\) or \(\mathbf{A}^H\) is being used.

        !--------------------------------------
        !-----     Internal variables     -----
        !--------------------------------------

        ! Options.
        integer :: kdim, maxiter
        real(sp) :: tol
        logical :: verbose, trans
        type(gmres_sp_opts) :: opts

        ! Krylov subspace
        class(abstract_vector_csp), allocatable :: V(:)
        ! Hessenberg matrix.
        complex(sp), allocatable :: H(:, :)
        ! Least-squares variables.
        complex(sp), allocatable :: y(:), e(:)
        complex(sp) :: beta

        ! Preconditioner
        logical :: has_precond
        class(abstract_precond_csp), allocatable :: precond

        ! Miscellaneous.
        integer :: i, j, k
        complex(sp) :: alpha
        class(abstract_vector_csp), allocatable :: dx, wrk

        ! Deals with the optional args.
        if (present(options)) then
            opts = gmres_sp_opts( &
                        kdim    = options%kdim, &
                        maxiter = options%maxiter, &
                        atol    = options%atol, &
                        rtol    = options%rtol, &
                        verbose = options%verbose &
                    )
        else
            opts = gmres_sp_opts()
        endif

        kdim = opts%kdim ; maxiter = opts%maxiter
        tol = opts%atol + opts%rtol * b%norm() ; verbose = opts%verbose
        trans = optval(transpose, .false.)

        ! Deals with the preconditioner.
        if (present(preconditioner)) then
            has_precond = .true.
            allocate(precond, source=preconditioner)
        else
            has_precond = .false.
        endif

        ! Initialize working variables.
        allocate(wrk, source=b) ; call wrk%zero()
        allocate(V(1:kdim+1), source=b) ; call initialize_krylov_subspace(V)
        allocate(H(kdim+1, kdim)) ; H = 0.0_sp
        allocate(y(kdim)) ; y = 0.0_sp
        allocate(e(kdim+1)) ; e = 0.0_sp

        info = 0

        ! Initial Krylov vector.
        if (x%norm() > 0) then
            if (trans) then
                call A%rmatvec(x, V(1))
            else
                call A%matvec(x, V(1))
            endif
        endif

        call V(1)%sub(b) ; call V(1)%chsgn()
        beta = V(1)%norm() ; call V(1)%scal(1.0_sp/beta)

        ! Iterative solver.
        gmres_iter : do i = 1, maxiter
            ! Zero-out variables.
            H = 0.0_sp ; y = 0.0_sp ; e = 0.0_sp ; e(1) = beta
            call initialize_krylov_subspace(V(2:kdim+1))

            ! Arnoldi factorization.
            arnoldi_fact: do k = 1, kdim
                ! Preconditioner.
                wrk = V(k) ; if (has_precond) call precond%apply(wrk)

                ! Matrix-vector product.
                if (trans) then
                    call A%rmatvec(wrk, V(k+1))
                else
                    call A%matvec(wrk, V(k+1))
                endif

                ! Gram-Schmid orthogonalization (twice is enough).
                do j = 1, k
                    alpha = V(j)%dot(V(k+1))
                    call V(k+1)%axpby(one_csp, V(j), -alpha)
                    H(j, k) = alpha
                enddo
                do j = 1, k
                    alpha = V(j)%dot(V(k+1))
                    call V(k+1)%axpby(one_csp, V(j), -alpha)
                    H(j, k) = H(j, k) + alpha
                enddo

                ! Update Hessenberg matrix and normalize residual Krylov vector.
                H(k+1, k) = V(k+1)%norm()
                if (abs(H(k+1, k)) > tol) then
                    call V(k+1)%scal(1.0_sp / H(k+1, k))
                endif

                ! Least-squares problem.
                call lstsq(H(1:k+1, 1:k), e(1:k+1), y(1:k))

                ! Compute residual.
                beta = norm2(abs(e(1:k+1) - matmul(H(1:k+1, 1:k), y(1:k))))

                ! Current number of iterations performed.
                info = info + 1

                ! Check convergence.
                if (abs(beta) <= tol) then
                    exit arnoldi_fact
                endif
            enddo arnoldi_fact

            ! Update solution.
            k = min(k, kdim)
            if (allocated(dx) .eqv. .false.) allocate(dx, source=x); call dx%zero()
            do j = 1, k
                call dx%axpby(one_csp, V(j), y(j))
            enddo
            if (has_precond) call precond%apply(dx)
            call x%add(dx)

            ! Recompute residual for sanity check.
            if (trans) then
                call A%rmatvec(x, v(1))
            else
                call A%matvec(x, v(1))
            endif
            call v(1)%sub(b) ; call v(1)%chsgn()

            ! Initialize new starting Krylov vector if needed.
            beta = v(1)%norm() ; call v(1)%scal(1.0_sp / beta)

            ! Exit gmres if desired accuracy is reached.
            if (abs(beta) <= tol) exit gmres_iter

        enddo gmres_iter

        return
    end subroutine gmres_csp

    subroutine gmres_cdp(A, b, x, info, preconditioner, options, transpose)
        class(abstract_linop_cdp), intent(in) :: A
        !! Linear operator to be inverted.
        class(abstract_vector_cdp), intent(in) :: b
        !! Right-hand side vector.
        class(abstract_vector_cdp), intent(out) :: x
        !! Solution vector.
        integer, intent(out) :: info
        !! Information flag.
        class(abstract_precond_cdp), optional, intent(in) :: preconditioner
        !! Preconditioner (optional).
        type(gmres_dp_opts), optional, intent(in) :: options
        !! GMRES options.   
        logical, optional, intent(in) :: transpose
        !! Whether \(\mathbf{A}\) or \(\mathbf{A}^H\) is being used.

        !--------------------------------------
        !-----     Internal variables     -----
        !--------------------------------------

        ! Options.
        integer :: kdim, maxiter
        real(dp) :: tol
        logical :: verbose, trans
        type(gmres_dp_opts) :: opts

        ! Krylov subspace
        class(abstract_vector_cdp), allocatable :: V(:)
        ! Hessenberg matrix.
        complex(dp), allocatable :: H(:, :)
        ! Least-squares variables.
        complex(dp), allocatable :: y(:), e(:)
        complex(dp) :: beta

        ! Preconditioner
        logical :: has_precond
        class(abstract_precond_cdp), allocatable :: precond

        ! Miscellaneous.
        integer :: i, j, k
        complex(dp) :: alpha
        class(abstract_vector_cdp), allocatable :: dx, wrk

        ! Deals with the optional args.
        if (present(options)) then
            opts = gmres_dp_opts( &
                        kdim    = options%kdim, &
                        maxiter = options%maxiter, &
                        atol    = options%atol, &
                        rtol    = options%rtol, &
                        verbose = options%verbose &
                    )
        else
            opts = gmres_dp_opts()
        endif

        kdim = opts%kdim ; maxiter = opts%maxiter
        tol = opts%atol + opts%rtol * b%norm() ; verbose = opts%verbose
        trans = optval(transpose, .false.)

        ! Deals with the preconditioner.
        if (present(preconditioner)) then
            has_precond = .true.
            allocate(precond, source=preconditioner)
        else
            has_precond = .false.
        endif

        ! Initialize working variables.
        allocate(wrk, source=b) ; call wrk%zero()
        allocate(V(1:kdim+1), source=b) ; call initialize_krylov_subspace(V)
        allocate(H(kdim+1, kdim)) ; H = 0.0_dp
        allocate(y(kdim)) ; y = 0.0_dp
        allocate(e(kdim+1)) ; e = 0.0_dp

        info = 0

        ! Initial Krylov vector.
        if (x%norm() > 0) then
            if (trans) then
                call A%rmatvec(x, V(1))
            else
                call A%matvec(x, V(1))
            endif
        endif

        call V(1)%sub(b) ; call V(1)%chsgn()
        beta = V(1)%norm() ; call V(1)%scal(1.0_dp/beta)

        ! Iterative solver.
        gmres_iter : do i = 1, maxiter
            ! Zero-out variables.
            H = 0.0_dp ; y = 0.0_dp ; e = 0.0_dp ; e(1) = beta
            call initialize_krylov_subspace(V(2:kdim+1))

            ! Arnoldi factorization.
            arnoldi_fact: do k = 1, kdim
                ! Preconditioner.
                wrk = V(k) ; if (has_precond) call precond%apply(wrk)

                ! Matrix-vector product.
                if (trans) then
                    call A%rmatvec(wrk, V(k+1))
                else
                    call A%matvec(wrk, V(k+1))
                endif

                ! Gram-Schmid orthogonalization (twice is enough).
                do j = 1, k
                    alpha = V(j)%dot(V(k+1))
                    call V(k+1)%axpby(one_cdp, V(j), -alpha)
                    H(j, k) = alpha
                enddo
                do j = 1, k
                    alpha = V(j)%dot(V(k+1))
                    call V(k+1)%axpby(one_cdp, V(j), -alpha)
                    H(j, k) = H(j, k) + alpha
                enddo

                ! Update Hessenberg matrix and normalize residual Krylov vector.
                H(k+1, k) = V(k+1)%norm()
                if (abs(H(k+1, k)) > tol) then
                    call V(k+1)%scal(1.0_dp / H(k+1, k))
                endif

                ! Least-squares problem.
                call lstsq(H(1:k+1, 1:k), e(1:k+1), y(1:k))

                ! Compute residual.
                beta = norm2(abs(e(1:k+1) - matmul(H(1:k+1, 1:k), y(1:k))))

                ! Current number of iterations performed.
                info = info + 1

                ! Check convergence.
                if (abs(beta) <= tol) then
                    exit arnoldi_fact
                endif
            enddo arnoldi_fact

            ! Update solution.
            k = min(k, kdim)
            if (allocated(dx) .eqv. .false.) allocate(dx, source=x); call dx%zero()
            do j = 1, k
                call dx%axpby(one_cdp, V(j), y(j))
            enddo
            if (has_precond) call precond%apply(dx)
            call x%add(dx)

            ! Recompute residual for sanity check.
            if (trans) then
                call A%rmatvec(x, v(1))
            else
                call A%matvec(x, v(1))
            endif
            call v(1)%sub(b) ; call v(1)%chsgn()

            ! Initialize new starting Krylov vector if needed.
            beta = v(1)%norm() ; call v(1)%scal(1.0_dp / beta)

            ! Exit gmres if desired accuracy is reached.
            if (abs(beta) <= tol) exit gmres_iter

        enddo gmres_iter

        return
    end subroutine gmres_cdp


end module lightkrylov_IterativeSolvers
