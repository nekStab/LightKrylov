submodule (lightkrylov_iterativesolvers) hermitian_eigensolvers
    use stdlib_strings, only: padr
    use stdlib_linalg, only: eigh
    implicit none
    character(len=*), parameter :: eighs_output = 'eighs_output.txt'
contains

    !----- Utility functions -----
    elemental pure function eigenvalue_residual_rsp(beta, x) result(residual)
        !! Computes the residual associated with a Ritz eigenpair.
        real(sp), intent(in) :: beta
        !! Norm of the residual Krylov vector.
        real(sp), intent(in) :: x
        !! Last entry of the low-dimensional Ritz eigenvector.
        real(sp) :: residual
        !! Residual associated to the corresponding Ritz eigenpair.
        residual = abs(beta*x)
        return
    end function eigenvalue_residual_rsp
    elemental pure function eigenvalue_residual_rdp(beta, x) result(residual)
        !! Computes the residual associated with a Ritz eigenpair.
        real(dp), intent(in) :: beta
        !! Norm of the residual Krylov vector.
        real(dp), intent(in) :: x
        !! Last entry of the low-dimensional Ritz eigenvector.
        real(dp) :: residual
        !! Residual associated to the corresponding Ritz eigenpair.
        residual = abs(beta*x)
        return
    end function eigenvalue_residual_rdp
    elemental pure function eigenvalue_residual_csp(beta, x) result(residual)
        !! Computes the residual associated with a Ritz eigenpair.
        complex(sp), intent(in) :: beta
        !! Norm of the residual Krylov vector.
        complex(sp), intent(in) :: x
        !! Last entry of the low-dimensional Ritz eigenvector.
        real(sp) :: residual
        !! Residual associated to the corresponding Ritz eigenpair.
        residual = abs(beta*x)
        return
    end function eigenvalue_residual_csp
    elemental pure function eigenvalue_residual_cdp(beta, x) result(residual)
        !! Computes the residual associated with a Ritz eigenpair.
        complex(dp), intent(in) :: beta
        !! Norm of the residual Krylov vector.
        complex(dp), intent(in) :: x
        !! Last entry of the low-dimensional Ritz eigenvector.
        real(dp) :: residual
        !! Residual associated to the corresponding Ritz eigenpair.
        residual = abs(beta*x)
        return
    end function eigenvalue_residual_cdp
   
    !--------------------------------------------------
    !-----     ABSTRACT HERMITIAN EIGENSOLVER     -----
    !--------------------------------------------------

    module procedure eighs_rsp
        class(abstract_vector_rsp), allocatable :: Xwrk(:)
        ! Krylov subspace.
        integer :: kdim_
        ! Krylov subspace dimension.
        real(sp), allocatable :: T(:, :)
        ! Tridiagonal matrix.
        real(sp), allocatable :: eigvecs_wrk(:, :)
        ! Working array for the Ritz eigenvectors.
        real(sp), allocatable :: eigvals_wrk(:)
        ! Working array for the Ritz eigenvalues.
        real(sp), allocatable :: residuals_wrk(:)
        ! Working array for the Ritz residuals.
        real(sp) :: x0_norm

        ! Miscellaneous.
        character(len=*), parameter :: this_procedure = 'eighs_rsp'
        integer :: i, j, k, nev, conv
        real(sp) :: tol
        real(sp) :: beta
        logical :: outpost
        character(len=256) :: msg

        if (time_lightkrylov()) call timer%start(this_procedure)
        ! Deaks with the optional args.
        nev = size(X)
        kdim_   = optval(kdim, 4*nev)
        tol     = optval(tolerance, rtol_sp)
        outpost = optval(write_intermediate, .false.)

        ! Allocate working variables.
        allocate(Xwrk(kdim_+1), mold=X(1)) ; call zero_basis(Xwrk)
        if (present(x0)) then
            call copy(Xwrk(1), x0)
            x0_norm = x0%norm(); call Xwrk(1)%scal(one_rsp/x0_norm)
        else
            call Xwrk(1)%rand(.true.)
        endif
        allocate(T(kdim_+1, kdim_)) ; T = zero_rsp
        allocate(eigvecs_wrk(kdim_, kdim_)) ; eigvecs_wrk = zero_rsp
        allocate(eigvals_wrk(kdim_)) ; eigvals_wrk = 0.0_sp
        allocate(residuals_wrk(kdim_)) ; residuals_wrk = 0.0_sp

        ! Ritz eigenpairs computation.
        lanczos_iter : do k = 1, kdim_
            ! Symmetric Lanczos step.
            call lanczos(A, Xwrk, T, info, kstart=k, kend=k)
            call check_info(info, 'lanczos', this_module, this_procedure)

            ! Spectral decomposition of the k x k tridiagonal matrix.
            eigvals_wrk = 0.0_sp ; eigvecs_wrk = zero_rsp
            call eigh(T(:k, :k), eigvals_wrk(:k), vectors=eigvecs_wrk(:k, :k))

            ! Compute residuals.
            beta = T(k+1, k)
            residuals_wrk(:k) = eigenvalue_residual_rsp(beta, eigvecs_wrk(k, :k))

            ! Check convergence.
            conv = count(residuals_wrk(:k) < tol)
            write(msg,'(I0,A,I0,A,I0,A)') conv, '/', nev, ' eigenvalues converged after ', k, &
                            & ' iterations of the Lanczos process.'
            call log_information(msg, this_module, this_procedure)
            if (outpost) call write_results_rsp(eighs_output, eigvals_wrk(:k), residuals_wrk(:k), tol)
            if (conv >= nev) exit lanczos_iter
        enddo lanczos_iter

        !--------------------------------
        !-----     POST-PROCESS     -----
        !--------------------------------

        block
            integer :: indices(kdim_)
            call sort_index(eigvals_wrk, indices, reverse=.true.)
            eigvecs_wrk = eigvecs_wrk(:, indices)
            residuals_wrk = residuals_wrk(indices)
            ! Store converged eigenvalues.
            eigvals = eigvals_wrk(:nev) ; residuals = residuals_wrk(:nev)
        end block

        ! Construct eigenvectors.
        k = min(k, kdim_)
        do i = 1, nev
            call X(i)%zero()
            do j = 1, k
                call X(i)%axpby(eigvecs_wrk(j, i), Xwrk(j), one_rsp)
            enddo
        enddo
        
        info = k
        if (time_lightkrylov()) call timer%stop(this_procedure)
        
        return
    end procedure 

    module procedure eighs_rdp
        class(abstract_vector_rdp), allocatable :: Xwrk(:)
        ! Krylov subspace.
        integer :: kdim_
        ! Krylov subspace dimension.
        real(dp), allocatable :: T(:, :)
        ! Tridiagonal matrix.
        real(dp), allocatable :: eigvecs_wrk(:, :)
        ! Working array for the Ritz eigenvectors.
        real(dp), allocatable :: eigvals_wrk(:)
        ! Working array for the Ritz eigenvalues.
        real(dp), allocatable :: residuals_wrk(:)
        ! Working array for the Ritz residuals.
        real(dp) :: x0_norm

        ! Miscellaneous.
        character(len=*), parameter :: this_procedure = 'eighs_rdp'
        integer :: i, j, k, nev, conv
        real(dp) :: tol
        real(dp) :: beta
        logical :: outpost
        character(len=256) :: msg

        if (time_lightkrylov()) call timer%start(this_procedure)
        ! Deaks with the optional args.
        nev = size(X)
        kdim_   = optval(kdim, 4*nev)
        tol     = optval(tolerance, rtol_dp)
        outpost = optval(write_intermediate, .false.)

        ! Allocate working variables.
        allocate(Xwrk(kdim_+1), mold=X(1)) ; call zero_basis(Xwrk)
        if (present(x0)) then
            call copy(Xwrk(1), x0)
            x0_norm = x0%norm(); call Xwrk(1)%scal(one_rdp/x0_norm)
        else
            call Xwrk(1)%rand(.true.)
        endif
        allocate(T(kdim_+1, kdim_)) ; T = zero_rdp
        allocate(eigvecs_wrk(kdim_, kdim_)) ; eigvecs_wrk = zero_rdp
        allocate(eigvals_wrk(kdim_)) ; eigvals_wrk = 0.0_dp
        allocate(residuals_wrk(kdim_)) ; residuals_wrk = 0.0_dp

        ! Ritz eigenpairs computation.
        lanczos_iter : do k = 1, kdim_
            ! Symmetric Lanczos step.
            call lanczos(A, Xwrk, T, info, kstart=k, kend=k)
            call check_info(info, 'lanczos', this_module, this_procedure)

            ! Spectral decomposition of the k x k tridiagonal matrix.
            eigvals_wrk = 0.0_dp ; eigvecs_wrk = zero_rdp
            call eigh(T(:k, :k), eigvals_wrk(:k), vectors=eigvecs_wrk(:k, :k))

            ! Compute residuals.
            beta = T(k+1, k)
            residuals_wrk(:k) = eigenvalue_residual_rdp(beta, eigvecs_wrk(k, :k))

            ! Check convergence.
            conv = count(residuals_wrk(:k) < tol)
            write(msg,'(I0,A,I0,A,I0,A)') conv, '/', nev, ' eigenvalues converged after ', k, &
                            & ' iterations of the Lanczos process.'
            call log_information(msg, this_module, this_procedure)
            if (outpost) call write_results_rdp(eighs_output, eigvals_wrk(:k), residuals_wrk(:k), tol)
            if (conv >= nev) exit lanczos_iter
        enddo lanczos_iter

        !--------------------------------
        !-----     POST-PROCESS     -----
        !--------------------------------

        block
            integer :: indices(kdim_)
            call sort_index(eigvals_wrk, indices, reverse=.true.)
            eigvecs_wrk = eigvecs_wrk(:, indices)
            residuals_wrk = residuals_wrk(indices)
            ! Store converged eigenvalues.
            eigvals = eigvals_wrk(:nev) ; residuals = residuals_wrk(:nev)
        end block

        ! Construct eigenvectors.
        k = min(k, kdim_)
        do i = 1, nev
            call X(i)%zero()
            do j = 1, k
                call X(i)%axpby(eigvecs_wrk(j, i), Xwrk(j), one_rdp)
            enddo
        enddo
        
        info = k
        if (time_lightkrylov()) call timer%stop(this_procedure)
        
        return
    end procedure 

    module procedure eighs_csp
        class(abstract_vector_csp), allocatable :: Xwrk(:)
        ! Krylov subspace.
        integer :: kdim_
        ! Krylov subspace dimension.
        complex(sp), allocatable :: T(:, :)
        ! Tridiagonal matrix.
        complex(sp), allocatable :: eigvecs_wrk(:, :)
        ! Working array for the Ritz eigenvectors.
        real(sp), allocatable :: eigvals_wrk(:)
        ! Working array for the Ritz eigenvalues.
        real(sp), allocatable :: residuals_wrk(:)
        ! Working array for the Ritz residuals.
        real(sp) :: x0_norm

        ! Miscellaneous.
        character(len=*), parameter :: this_procedure = 'eighs_csp'
        integer :: i, j, k, nev, conv
        real(sp) :: tol
        complex(sp) :: beta
        logical :: outpost
        character(len=256) :: msg

        if (time_lightkrylov()) call timer%start(this_procedure)
        ! Deaks with the optional args.
        nev = size(X)
        kdim_   = optval(kdim, 4*nev)
        tol     = optval(tolerance, rtol_sp)
        outpost = optval(write_intermediate, .false.)

        ! Allocate working variables.
        allocate(Xwrk(kdim_+1), mold=X(1)) ; call zero_basis(Xwrk)
        if (present(x0)) then
            call copy(Xwrk(1), x0)
            x0_norm = x0%norm(); call Xwrk(1)%scal(one_csp/x0_norm)
        else
            call Xwrk(1)%rand(.true.)
        endif
        allocate(T(kdim_+1, kdim_)) ; T = zero_csp
        allocate(eigvecs_wrk(kdim_, kdim_)) ; eigvecs_wrk = zero_csp
        allocate(eigvals_wrk(kdim_)) ; eigvals_wrk = 0.0_sp
        allocate(residuals_wrk(kdim_)) ; residuals_wrk = 0.0_sp

        ! Ritz eigenpairs computation.
        lanczos_iter : do k = 1, kdim_
            ! Symmetric Lanczos step.
            call lanczos(A, Xwrk, T, info, kstart=k, kend=k)
            call check_info(info, 'lanczos', this_module, this_procedure)

            ! Spectral decomposition of the k x k tridiagonal matrix.
            eigvals_wrk = 0.0_sp ; eigvecs_wrk = zero_csp
            call eigh(T(:k, :k), eigvals_wrk(:k), vectors=eigvecs_wrk(:k, :k))

            ! Compute residuals.
            beta = T(k+1, k)
            residuals_wrk(:k) = eigenvalue_residual_csp(beta, eigvecs_wrk(k, :k))

            ! Check convergence.
            conv = count(residuals_wrk(:k) < tol)
            write(msg,'(I0,A,I0,A,I0,A)') conv, '/', nev, ' eigenvalues converged after ', k, &
                            & ' iterations of the Lanczos process.'
            call log_information(msg, this_module, this_procedure)
            if (outpost) call write_results_rsp(eighs_output, eigvals_wrk(:k), residuals_wrk(:k), tol)
            if (conv >= nev) exit lanczos_iter
        enddo lanczos_iter

        !--------------------------------
        !-----     POST-PROCESS     -----
        !--------------------------------

        block
            integer :: indices(kdim_)
            call sort_index(eigvals_wrk, indices, reverse=.true.)
            eigvecs_wrk = eigvecs_wrk(:, indices)
            residuals_wrk = residuals_wrk(indices)
            ! Store converged eigenvalues.
            eigvals = eigvals_wrk(:nev) ; residuals = residuals_wrk(:nev)
        end block

        ! Construct eigenvectors.
        k = min(k, kdim_)
        do i = 1, nev
            call X(i)%zero()
            do j = 1, k
                call X(i)%axpby(eigvecs_wrk(j, i), Xwrk(j), one_csp)
            enddo
        enddo
        
        info = k
        if (time_lightkrylov()) call timer%stop(this_procedure)
        
        return
    end procedure 

    module procedure eighs_cdp
        class(abstract_vector_cdp), allocatable :: Xwrk(:)
        ! Krylov subspace.
        integer :: kdim_
        ! Krylov subspace dimension.
        complex(dp), allocatable :: T(:, :)
        ! Tridiagonal matrix.
        complex(dp), allocatable :: eigvecs_wrk(:, :)
        ! Working array for the Ritz eigenvectors.
        real(dp), allocatable :: eigvals_wrk(:)
        ! Working array for the Ritz eigenvalues.
        real(dp), allocatable :: residuals_wrk(:)
        ! Working array for the Ritz residuals.
        real(dp) :: x0_norm

        ! Miscellaneous.
        character(len=*), parameter :: this_procedure = 'eighs_cdp'
        integer :: i, j, k, nev, conv
        real(dp) :: tol
        complex(dp) :: beta
        logical :: outpost
        character(len=256) :: msg

        if (time_lightkrylov()) call timer%start(this_procedure)
        ! Deaks with the optional args.
        nev = size(X)
        kdim_   = optval(kdim, 4*nev)
        tol     = optval(tolerance, rtol_dp)
        outpost = optval(write_intermediate, .false.)

        ! Allocate working variables.
        allocate(Xwrk(kdim_+1), mold=X(1)) ; call zero_basis(Xwrk)
        if (present(x0)) then
            call copy(Xwrk(1), x0)
            x0_norm = x0%norm(); call Xwrk(1)%scal(one_cdp/x0_norm)
        else
            call Xwrk(1)%rand(.true.)
        endif
        allocate(T(kdim_+1, kdim_)) ; T = zero_cdp
        allocate(eigvecs_wrk(kdim_, kdim_)) ; eigvecs_wrk = zero_cdp
        allocate(eigvals_wrk(kdim_)) ; eigvals_wrk = 0.0_dp
        allocate(residuals_wrk(kdim_)) ; residuals_wrk = 0.0_dp

        ! Ritz eigenpairs computation.
        lanczos_iter : do k = 1, kdim_
            ! Symmetric Lanczos step.
            call lanczos(A, Xwrk, T, info, kstart=k, kend=k)
            call check_info(info, 'lanczos', this_module, this_procedure)

            ! Spectral decomposition of the k x k tridiagonal matrix.
            eigvals_wrk = 0.0_dp ; eigvecs_wrk = zero_cdp
            call eigh(T(:k, :k), eigvals_wrk(:k), vectors=eigvecs_wrk(:k, :k))

            ! Compute residuals.
            beta = T(k+1, k)
            residuals_wrk(:k) = eigenvalue_residual_cdp(beta, eigvecs_wrk(k, :k))

            ! Check convergence.
            conv = count(residuals_wrk(:k) < tol)
            write(msg,'(I0,A,I0,A,I0,A)') conv, '/', nev, ' eigenvalues converged after ', k, &
                            & ' iterations of the Lanczos process.'
            call log_information(msg, this_module, this_procedure)
            if (outpost) call write_results_rdp(eighs_output, eigvals_wrk(:k), residuals_wrk(:k), tol)
            if (conv >= nev) exit lanczos_iter
        enddo lanczos_iter

        !--------------------------------
        !-----     POST-PROCESS     -----
        !--------------------------------

        block
            integer :: indices(kdim_)
            call sort_index(eigvals_wrk, indices, reverse=.true.)
            eigvecs_wrk = eigvecs_wrk(:, indices)
            residuals_wrk = residuals_wrk(indices)
            ! Store converged eigenvalues.
            eigvals = eigvals_wrk(:nev) ; residuals = residuals_wrk(:nev)
        end block

        ! Construct eigenvectors.
        k = min(k, kdim_)
        do i = 1, nev
            call X(i)%zero()
            do j = 1, k
                call X(i)%axpby(eigvecs_wrk(j, i), Xwrk(j), one_cdp)
            enddo
        enddo
        
        info = k
        if (time_lightkrylov()) call timer%stop(this_procedure)
        
        return
    end procedure 

   

end submodule

