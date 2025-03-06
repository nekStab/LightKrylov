submodule (lightkrylov_iterativesolvers) svds_solver
    use stdlib_linalg, only: hermitian, svd
    implicit none
    character(len=*), parameter :: svds_output = 'svds_output.txt'
contains

    !----- Utility functions -----
    elemental pure function svd_residual_rsp(beta, x) result(residual)
        !! Computes the residual associated with a Ritz eigenpair.
        real(sp), intent(in) :: beta
        !! Norm of the residual Krylov vector.
        real(sp), intent(in) :: x
        !! Last entry of the low-dimensional Ritz eigenvector.
        real(sp) :: residual
        !! Residual associated to the corresponding Ritz eigenpair.
        residual = abs(beta*x)
        return
    end function svd_residual_rsp
    elemental pure function svd_residual_rdp(beta, x) result(residual)
        !! Computes the residual associated with a Ritz eigenpair.
        real(dp), intent(in) :: beta
        !! Norm of the residual Krylov vector.
        real(dp), intent(in) :: x
        !! Last entry of the low-dimensional Ritz eigenvector.
        real(dp) :: residual
        !! Residual associated to the corresponding Ritz eigenpair.
        residual = abs(beta*x)
        return
    end function svd_residual_rdp
    elemental pure function svd_residual_csp(beta, x) result(residual)
        !! Computes the residual associated with a Ritz eigenpair.
        complex(sp), intent(in) :: beta
        !! Norm of the residual Krylov vector.
        complex(sp), intent(in) :: x
        !! Last entry of the low-dimensional Ritz eigenvector.
        real(sp) :: residual
        !! Residual associated to the corresponding Ritz eigenpair.
        residual = abs(beta*x)
        return
    end function svd_residual_csp
    elemental pure function svd_residual_cdp(beta, x) result(residual)
        !! Computes the residual associated with a Ritz eigenpair.
        complex(dp), intent(in) :: beta
        !! Norm of the residual Krylov vector.
        complex(dp), intent(in) :: x
        !! Last entry of the low-dimensional Ritz eigenvector.
        real(dp) :: residual
        !! Residual associated to the corresponding Ritz eigenpair.
        residual = abs(beta*x)
        return
    end function svd_residual_cdp

    !----------------------------------------
    !-----     ITERATIVE SVD SOLVER     -----
    !----------------------------------------

    module procedure svds_rsp
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
        character(len=*), parameter :: this_procedure = 'svds_rsp'
        integer :: nsv, conv
        integer :: i, j, k
        real(sp) :: tol, u0_norm
        logical :: outpost
        character(len=256) :: msg

        if (time_lightkrylov()) call timer%start(this_procedure)
        ! Deals with the optional arguments.
        nsv = size(U)
        kdim_   = optval(kdim, 4*nsv)
        tol     = optval(tolerance, rtol_sp)
        outpost = optval(write_intermediate, .false.)

        ! Allocate working variables.
        allocate(Uwrk(kdim_+1), mold=U(1)) ; call zero_basis(Uwrk)
        if (present(u0)) then
            call copy(Uwrk(1), u0)
            u0_norm = u0%norm(); call Uwrk(1)%scal(one_rsp/u0_norm)
        else
            call Uwrk(1)%rand(.true.)
        endif
        allocate(Vwrk(kdim_+1), mold=V(1)) ; call zero_basis(Vwrk)
        allocate(svdvals_wrk(kdim_)) ; svdvals_wrk = 0.0_sp
        allocate(umat(kdim_, kdim_)) ; umat = 0.0_sp
        allocate(vmat(kdim_, kdim_)) ; vmat = 0.0_sp
        allocate(residuals_wrk(kdim_)) ; residuals_wrk = 0.0_sp
        allocate(B(kdim_+1, kdim_)) ; B = 0.0_sp

        info = 0

        ! Ritz singular triplets computation.
        lanczos_iter : do k = 1, kdim_
            ! Lanczos bidiag. step.
            call bidiagonalization(A, Uwrk, Vwrk, B, info, kstart=k, kend=k, tol=tol)
            call check_info(info, 'bidiagonalization', this_module, this_procedure)

            ! SVD of the k x k bidiagonal matrix and residual computation.
            svdvals_wrk = 0.0_sp ; umat = 0.0_sp ; vmat = 0.0_sp

            call svd(B(:k, :k), svdvals_wrk(:k), umat(:k, :k), vmat(:k, :k))
            vmat(:k, :k) = hermitian(vmat(:k, :k))

            residuals_wrk(:k) = svd_residual_rsp(B(k+1, k), vmat(k, :k))

            ! Check for convergence.
            conv = count(residuals_wrk(:k) < tol)
            write(msg,'(I0,A,I0,A,I0,A)') conv, '/', nsv, ' singular values converged after ', k, &
                            & ' iterations of the Lanczos process.'
            call log_information(msg, this_module, this_procedure)
            if (outpost) call write_results_rsp(svds_output, svdvals_wrk(:k), residuals_wrk(:k), tol)
            if (conv >= nsv) exit lanczos_iter
        enddo lanczos_iter

        !--------------------------------
        !-----     POST-PROCESS     -----
        !--------------------------------

        ! Singular values.
        S = svdvals_wrk(:nsv) ; residuals = residuals_wrk(:nsv)

        ! Singular vectors.
        k = min(k, kdim_) ; info = k
        call zero_basis(U) ; call zero_basis(V)
        do i = 1, nsv
            do j = 1, k
                call U(i)%axpby(umat(j, i), Uwrk(j), one_rsp)
                call V(i)%axpby(vmat(j, i), Vwrk(j), one_rsp)
            enddo
        enddo
        if (time_lightkrylov()) call timer%stop(this_procedure)
        
        return
    end procedure

    module procedure svds_rdp
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
        character(len=*), parameter :: this_procedure = 'svds_rdp'
        integer :: nsv, conv
        integer :: i, j, k
        real(dp) :: tol, u0_norm
        logical :: outpost
        character(len=256) :: msg

        if (time_lightkrylov()) call timer%start(this_procedure)
        ! Deals with the optional arguments.
        nsv = size(U)
        kdim_   = optval(kdim, 4*nsv)
        tol     = optval(tolerance, rtol_dp)
        outpost = optval(write_intermediate, .false.)

        ! Allocate working variables.
        allocate(Uwrk(kdim_+1), mold=U(1)) ; call zero_basis(Uwrk)
        if (present(u0)) then
            call copy(Uwrk(1), u0)
            u0_norm = u0%norm(); call Uwrk(1)%scal(one_rdp/u0_norm)
        else
            call Uwrk(1)%rand(.true.)
        endif
        allocate(Vwrk(kdim_+1), mold=V(1)) ; call zero_basis(Vwrk)
        allocate(svdvals_wrk(kdim_)) ; svdvals_wrk = 0.0_dp
        allocate(umat(kdim_, kdim_)) ; umat = 0.0_dp
        allocate(vmat(kdim_, kdim_)) ; vmat = 0.0_dp
        allocate(residuals_wrk(kdim_)) ; residuals_wrk = 0.0_dp
        allocate(B(kdim_+1, kdim_)) ; B = 0.0_dp

        info = 0

        ! Ritz singular triplets computation.
        lanczos_iter : do k = 1, kdim_
            ! Lanczos bidiag. step.
            call bidiagonalization(A, Uwrk, Vwrk, B, info, kstart=k, kend=k, tol=tol)
            call check_info(info, 'bidiagonalization', this_module, this_procedure)

            ! SVD of the k x k bidiagonal matrix and residual computation.
            svdvals_wrk = 0.0_dp ; umat = 0.0_dp ; vmat = 0.0_dp

            call svd(B(:k, :k), svdvals_wrk(:k), umat(:k, :k), vmat(:k, :k))
            vmat(:k, :k) = hermitian(vmat(:k, :k))

            residuals_wrk(:k) = svd_residual_rdp(B(k+1, k), vmat(k, :k))

            ! Check for convergence.
            conv = count(residuals_wrk(:k) < tol)
            write(msg,'(I0,A,I0,A,I0,A)') conv, '/', nsv, ' singular values converged after ', k, &
                            & ' iterations of the Lanczos process.'
            call log_information(msg, this_module, this_procedure)
            if (outpost) call write_results_rdp(svds_output, svdvals_wrk(:k), residuals_wrk(:k), tol)
            if (conv >= nsv) exit lanczos_iter
        enddo lanczos_iter

        !--------------------------------
        !-----     POST-PROCESS     -----
        !--------------------------------

        ! Singular values.
        S = svdvals_wrk(:nsv) ; residuals = residuals_wrk(:nsv)

        ! Singular vectors.
        k = min(k, kdim_) ; info = k
        call zero_basis(U) ; call zero_basis(V)
        do i = 1, nsv
            do j = 1, k
                call U(i)%axpby(umat(j, i), Uwrk(j), one_rdp)
                call V(i)%axpby(vmat(j, i), Vwrk(j), one_rdp)
            enddo
        enddo
        if (time_lightkrylov()) call timer%stop(this_procedure)
        
        return
    end procedure

    module procedure svds_csp
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
        character(len=*), parameter :: this_procedure = 'svds_csp'
        integer :: nsv, conv
        integer :: i, j, k
        real(sp) :: tol, u0_norm
        logical :: outpost
        character(len=256) :: msg

        if (time_lightkrylov()) call timer%start(this_procedure)
        ! Deals with the optional arguments.
        nsv = size(U)
        kdim_   = optval(kdim, 4*nsv)
        tol     = optval(tolerance, rtol_sp)
        outpost = optval(write_intermediate, .false.)

        ! Allocate working variables.
        allocate(Uwrk(kdim_+1), mold=U(1)) ; call zero_basis(Uwrk)
        if (present(u0)) then
            call copy(Uwrk(1), u0)
            u0_norm = u0%norm(); call Uwrk(1)%scal(one_csp/u0_norm)
        else
            call Uwrk(1)%rand(.true.)
        endif
        allocate(Vwrk(kdim_+1), mold=V(1)) ; call zero_basis(Vwrk)
        allocate(svdvals_wrk(kdim_)) ; svdvals_wrk = 0.0_sp
        allocate(umat(kdim_, kdim_)) ; umat = 0.0_sp
        allocate(vmat(kdim_, kdim_)) ; vmat = 0.0_sp
        allocate(residuals_wrk(kdim_)) ; residuals_wrk = 0.0_sp
        allocate(B(kdim_+1, kdim_)) ; B = 0.0_sp

        info = 0

        ! Ritz singular triplets computation.
        lanczos_iter : do k = 1, kdim_
            ! Lanczos bidiag. step.
            call bidiagonalization(A, Uwrk, Vwrk, B, info, kstart=k, kend=k, tol=tol)
            call check_info(info, 'bidiagonalization', this_module, this_procedure)

            ! SVD of the k x k bidiagonal matrix and residual computation.
            svdvals_wrk = 0.0_sp ; umat = 0.0_sp ; vmat = 0.0_sp

            call svd(B(:k, :k), svdvals_wrk(:k), umat(:k, :k), vmat(:k, :k))
            vmat(:k, :k) = hermitian(vmat(:k, :k))

            residuals_wrk(:k) = svd_residual_csp(B(k+1, k), vmat(k, :k))

            ! Check for convergence.
            conv = count(residuals_wrk(:k) < tol)
            write(msg,'(I0,A,I0,A,I0,A)') conv, '/', nsv, ' singular values converged after ', k, &
                            & ' iterations of the Lanczos process.'
            call log_information(msg, this_module, this_procedure)
            if (outpost) call write_results_rsp(svds_output, svdvals_wrk(:k), residuals_wrk(:k), tol)
            if (conv >= nsv) exit lanczos_iter
        enddo lanczos_iter

        !--------------------------------
        !-----     POST-PROCESS     -----
        !--------------------------------

        ! Singular values.
        S = svdvals_wrk(:nsv) ; residuals = residuals_wrk(:nsv)

        ! Singular vectors.
        k = min(k, kdim_) ; info = k
        call zero_basis(U) ; call zero_basis(V)
        do i = 1, nsv
            do j = 1, k
                call U(i)%axpby(umat(j, i), Uwrk(j), one_csp)
                call V(i)%axpby(vmat(j, i), Vwrk(j), one_csp)
            enddo
        enddo
        if (time_lightkrylov()) call timer%stop(this_procedure)
        
        return
    end procedure

    module procedure svds_cdp
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
        character(len=*), parameter :: this_procedure = 'svds_cdp'
        integer :: nsv, conv
        integer :: i, j, k
        real(dp) :: tol, u0_norm
        logical :: outpost
        character(len=256) :: msg

        if (time_lightkrylov()) call timer%start(this_procedure)
        ! Deals with the optional arguments.
        nsv = size(U)
        kdim_   = optval(kdim, 4*nsv)
        tol     = optval(tolerance, rtol_dp)
        outpost = optval(write_intermediate, .false.)

        ! Allocate working variables.
        allocate(Uwrk(kdim_+1), mold=U(1)) ; call zero_basis(Uwrk)
        if (present(u0)) then
            call copy(Uwrk(1), u0)
            u0_norm = u0%norm(); call Uwrk(1)%scal(one_cdp/u0_norm)
        else
            call Uwrk(1)%rand(.true.)
        endif
        allocate(Vwrk(kdim_+1), mold=V(1)) ; call zero_basis(Vwrk)
        allocate(svdvals_wrk(kdim_)) ; svdvals_wrk = 0.0_dp
        allocate(umat(kdim_, kdim_)) ; umat = 0.0_dp
        allocate(vmat(kdim_, kdim_)) ; vmat = 0.0_dp
        allocate(residuals_wrk(kdim_)) ; residuals_wrk = 0.0_dp
        allocate(B(kdim_+1, kdim_)) ; B = 0.0_dp

        info = 0

        ! Ritz singular triplets computation.
        lanczos_iter : do k = 1, kdim_
            ! Lanczos bidiag. step.
            call bidiagonalization(A, Uwrk, Vwrk, B, info, kstart=k, kend=k, tol=tol)
            call check_info(info, 'bidiagonalization', this_module, this_procedure)

            ! SVD of the k x k bidiagonal matrix and residual computation.
            svdvals_wrk = 0.0_dp ; umat = 0.0_dp ; vmat = 0.0_dp

            call svd(B(:k, :k), svdvals_wrk(:k), umat(:k, :k), vmat(:k, :k))
            vmat(:k, :k) = hermitian(vmat(:k, :k))

            residuals_wrk(:k) = svd_residual_cdp(B(k+1, k), vmat(k, :k))

            ! Check for convergence.
            conv = count(residuals_wrk(:k) < tol)
            write(msg,'(I0,A,I0,A,I0,A)') conv, '/', nsv, ' singular values converged after ', k, &
                            & ' iterations of the Lanczos process.'
            call log_information(msg, this_module, this_procedure)
            if (outpost) call write_results_rdp(svds_output, svdvals_wrk(:k), residuals_wrk(:k), tol)
            if (conv >= nsv) exit lanczos_iter
        enddo lanczos_iter

        !--------------------------------
        !-----     POST-PROCESS     -----
        !--------------------------------

        ! Singular values.
        S = svdvals_wrk(:nsv) ; residuals = residuals_wrk(:nsv)

        ! Singular vectors.
        k = min(k, kdim_) ; info = k
        call zero_basis(U) ; call zero_basis(V)
        do i = 1, nsv
            do j = 1, k
                call U(i)%axpby(umat(j, i), Uwrk(j), one_cdp)
                call V(i)%axpby(vmat(j, i), Vwrk(j), one_cdp)
            enddo
        enddo
        if (time_lightkrylov()) call timer%stop(this_procedure)
        
        return
    end procedure

end submodule

