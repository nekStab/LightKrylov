submodule (lightkrylov_iterativesolvers) general_eigensolvers
    use stdlib_stats, only: median
    implicit none(type, external)
    character(len=*), parameter :: eigs_output = 'eigs_output.txt'
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
    end function eigenvalue_residual_cdp

   !------------------------------------------------------------
   !-----     ABSTRACT GENERAL EIGENVALUE COMPUTATIONS     -----
   !------------------------------------------------------------

   module procedure eigs_rsp
      ! Krylov subspace and Krylov subspace dimension.
      class(abstract_vector_rsp), allocatable :: Xwrk(:)
      integer :: kdim_, kstart
      ! Hessenberg matrix.
      real(sp), allocatable :: H(:, :)
      ! Working arrays for the eigenvectors and eigenvalues.
      real(sp), allocatable :: eigvecs_wrk(:, :)
      complex(sp), allocatable :: eigvals_wrk(:)
      real(sp), allocatable :: residuals_wrk(:)
      ! Miscellaneous.
      character(len=*), parameter :: this_procedure = 'eigs_rsp'
      integer :: nev, conv
      integer :: i, j, k, niter, krst, iostat
      real(sp) :: tol, x0_norm
      real(sp) :: beta
      real(sp) :: alpha
      logical :: outpost
      character(len=256) :: msg

      if (time_lightkrylov()) call timer%start(this_procedure)
      ! Deals with optional parameters.
      nev = size(X)
      kdim_   = optval(kdim, 4*nev)
      tol     = optval(tolerance, rtol_sp)
      outpost = optval(write_intermediate, .true.)

      ! Allocate eigenvalues.
      allocate(eigvals(nev), source=zero_csp, stat=iostat, errmsg=msg)
      call check_allocation(iostat, msg, this_module, this_procedure)

      ! Allocate working variables.
      allocate(Xwrk(kdim_+1), source=X(1), stat=iostat, errmsg=msg)
      call check_allocation(iostat, msg, this_module, this_procedure)
      call zero_basis(Xwrk)

      if (present(x0)) then
         call copy(Xwrk(1), x0)
         x0_norm = x0%norm(); call Xwrk(1)%scal(one_rsp/x0_norm)
      else
         call Xwrk(1)%rand(.true.)
      endif
      allocate(H(kdim_+1, kdim_), source=zero_rsp, stat=iostat, errmsg=msg)
      call check_allocation(iostat, msg, this_module, this_procedure)
      allocate(eigvecs_wrk(kdim_, kdim_), source=zero_rsp, stat=iostat, errmsg=msg)
      call check_allocation(iostat, msg, this_module, this_procedure)
      allocate(eigvals_wrk(kdim_), source=zero_csp, stat=iostat, errmsg=msg)
      call check_allocation(iostat, msg, this_module, this_procedure)
      allocate(residuals_wrk(kdim_), source=zero_rsp, stat=iostat, errmsg=msg)
      call check_allocation(iostat, msg, this_module, this_procedure)

      ! Ritz eigenpairs computation.
      H = 0.0_sp

      kstart = 1 ; conv = 0 ; niter = 0 ; krst = 1
      krylovschur: do while (conv < nev)

         arnoldi_factorization: do k = kstart, kdim_
               ! Arnoldi step.
               call arnoldi(A, Xwrk, H, info, kstart=k, kend=k, transpose=transpose)
               call check_info(info, 'arnoldi', this_module, this_procedure)

               ! Spectral decomposition of the k x k Hessenberg matrix.
               eigvals_wrk = 0.0_sp ; eigvecs_wrk = 0.0_sp
               if (time_lightkrylov()) call timer%start('eig')
               call eig(H(:k, :k), eigvecs_wrk(:k, :k), eigvals_wrk(:k))
               if (time_lightkrylov()) call timer%stop('eig')

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
                  residuals_wrk(i) = eigenvalue_residual_rsp(beta, alpha)
               enddo

               ! Check convergence.
               niter = niter + 1
               conv = count(residuals_wrk(:k) < tol)
               write(msg,'(I0,A,I0,A,I0,A)') conv, '/', nev, ' eigenvalues converged after ', niter, &
                           & ' steps of the Arnoldi process.'
               call log_information(msg, this_module, this_procedure)
               if (outpost) call write_results_csp(eigs_output, eigvals_wrk(:k), residuals_wrk(:k), tol)
               if (conv >= nev) exit arnoldi_factorization
         enddo arnoldi_factorization

         write(msg,'(I0,A,I0,A,I0,A)') conv, '/', nev, ' eigenvalues converged after ', krst, &
                           & ' Krylov-Schur restarts of the Arnoldi process.'
         call log_information(msg, this_module, this_procedure)
         ! Krylov-Schur restarting procedure.
         krst  = krst + 1
         call krylov_schur(kstart, Xwrk, H, median_selector) ; kstart = kstart + 1

      end do krylovschur

      !--------------------------------
      !-----     POST-PROCESS     -----
      !--------------------------------

      block
      integer :: indices(kdim_)
      real(sp) :: abs_eigvals(kdim_)

      ! Re-compute eigenvalues and eigenvectors.
      k = min(k, kdim_)
      if (time_lightkrylov()) call timer%start('eig')
      call eig(H(:k, :k), eigvecs_wrk(:k, :k), eigvals_wrk(:k))
      if (time_lightkrylov()) call timer%stop('eig')
      ! Sort eigenvalues.
      abs_eigvals = abs(eigvals_wrk) ; call sort_index(abs_eigvals, indices, reverse=.true.)
      eigvals_wrk = eigvals_wrk(indices) ; eigvecs_wrk = eigvecs_wrk(:, indices)
      residuals_wrk = residuals_wrk(indices)

      ! Store converged eigenvalues.
      eigvals = eigvals_wrk(:nev) ; residuals = residuals_wrk(:nev)
      end block

      ! Construct eigenvectors.
      do i = 1, nev
         call X(i)%zero()
         do j = 1, k
               call X(i)%axpby(eigvecs_wrk(j, i), Xwrk(j), one_rsp)
         enddo
      enddo

      info = niter
      if (time_lightkrylov()) call timer%stop(this_procedure)
   contains
      function median_selector(lambda) result(selected)
         implicit none(type, external)
         complex(sp), intent(in) :: lambda(:)
         logical, allocatable :: selected(:)
         selected = abs(lambda) > median(abs(lambda))
      end function median_selector
   end procedure eigs_rsp

   module procedure eigs_rdp
      ! Krylov subspace and Krylov subspace dimension.
      class(abstract_vector_rdp), allocatable :: Xwrk(:)
      integer :: kdim_, kstart
      ! Hessenberg matrix.
      real(dp), allocatable :: H(:, :)
      ! Working arrays for the eigenvectors and eigenvalues.
      real(dp), allocatable :: eigvecs_wrk(:, :)
      complex(dp), allocatable :: eigvals_wrk(:)
      real(dp), allocatable :: residuals_wrk(:)
      ! Miscellaneous.
      character(len=*), parameter :: this_procedure = 'eigs_rdp'
      integer :: nev, conv
      integer :: i, j, k, niter, krst, iostat
      real(dp) :: tol, x0_norm
      real(dp) :: beta
      real(dp) :: alpha
      logical :: outpost
      character(len=256) :: msg

      if (time_lightkrylov()) call timer%start(this_procedure)
      ! Deals with optional parameters.
      nev = size(X)
      kdim_   = optval(kdim, 4*nev)
      tol     = optval(tolerance, rtol_dp)
      outpost = optval(write_intermediate, .true.)

      ! Allocate eigenvalues.
      allocate(eigvals(nev), source=zero_cdp, stat=iostat, errmsg=msg)
      call check_allocation(iostat, msg, this_module, this_procedure)

      ! Allocate working variables.
      allocate(Xwrk(kdim_+1), source=X(1), stat=iostat, errmsg=msg)
      call check_allocation(iostat, msg, this_module, this_procedure)
      call zero_basis(Xwrk)

      if (present(x0)) then
         call copy(Xwrk(1), x0)
         x0_norm = x0%norm(); call Xwrk(1)%scal(one_rdp/x0_norm)
      else
         call Xwrk(1)%rand(.true.)
      endif
      allocate(H(kdim_+1, kdim_), source=zero_rdp, stat=iostat, errmsg=msg)
      call check_allocation(iostat, msg, this_module, this_procedure)
      allocate(eigvecs_wrk(kdim_, kdim_), source=zero_rdp, stat=iostat, errmsg=msg)
      call check_allocation(iostat, msg, this_module, this_procedure)
      allocate(eigvals_wrk(kdim_), source=zero_cdp, stat=iostat, errmsg=msg)
      call check_allocation(iostat, msg, this_module, this_procedure)
      allocate(residuals_wrk(kdim_), source=zero_rdp, stat=iostat, errmsg=msg)
      call check_allocation(iostat, msg, this_module, this_procedure)

      ! Ritz eigenpairs computation.
      H = 0.0_dp

      kstart = 1 ; conv = 0 ; niter = 0 ; krst = 1
      krylovschur: do while (conv < nev)

         arnoldi_factorization: do k = kstart, kdim_
               ! Arnoldi step.
               call arnoldi(A, Xwrk, H, info, kstart=k, kend=k, transpose=transpose)
               call check_info(info, 'arnoldi', this_module, this_procedure)

               ! Spectral decomposition of the k x k Hessenberg matrix.
               eigvals_wrk = 0.0_dp ; eigvecs_wrk = 0.0_dp
               if (time_lightkrylov()) call timer%start('eig')
               call eig(H(:k, :k), eigvecs_wrk(:k, :k), eigvals_wrk(:k))
               if (time_lightkrylov()) call timer%stop('eig')

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
                  residuals_wrk(i) = eigenvalue_residual_rdp(beta, alpha)
               enddo

               ! Check convergence.
               niter = niter + 1
               conv = count(residuals_wrk(:k) < tol)
               write(msg,'(I0,A,I0,A,I0,A)') conv, '/', nev, ' eigenvalues converged after ', niter, &
                           & ' steps of the Arnoldi process.'
               call log_information(msg, this_module, this_procedure)
               if (outpost) call write_results_cdp(eigs_output, eigvals_wrk(:k), residuals_wrk(:k), tol)
               if (conv >= nev) exit arnoldi_factorization
         enddo arnoldi_factorization

         write(msg,'(I0,A,I0,A,I0,A)') conv, '/', nev, ' eigenvalues converged after ', krst, &
                           & ' Krylov-Schur restarts of the Arnoldi process.'
         call log_information(msg, this_module, this_procedure)
         ! Krylov-Schur restarting procedure.
         krst  = krst + 1
         call krylov_schur(kstart, Xwrk, H, median_selector) ; kstart = kstart + 1

      end do krylovschur

      !--------------------------------
      !-----     POST-PROCESS     -----
      !--------------------------------

      block
      integer :: indices(kdim_)
      real(dp) :: abs_eigvals(kdim_)

      ! Re-compute eigenvalues and eigenvectors.
      k = min(k, kdim_)
      if (time_lightkrylov()) call timer%start('eig')
      call eig(H(:k, :k), eigvecs_wrk(:k, :k), eigvals_wrk(:k))
      if (time_lightkrylov()) call timer%stop('eig')
      ! Sort eigenvalues.
      abs_eigvals = abs(eigvals_wrk) ; call sort_index(abs_eigvals, indices, reverse=.true.)
      eigvals_wrk = eigvals_wrk(indices) ; eigvecs_wrk = eigvecs_wrk(:, indices)
      residuals_wrk = residuals_wrk(indices)

      ! Store converged eigenvalues.
      eigvals = eigvals_wrk(:nev) ; residuals = residuals_wrk(:nev)
      end block

      ! Construct eigenvectors.
      do i = 1, nev
         call X(i)%zero()
         do j = 1, k
               call X(i)%axpby(eigvecs_wrk(j, i), Xwrk(j), one_rdp)
         enddo
      enddo

      info = niter
      if (time_lightkrylov()) call timer%stop(this_procedure)
   contains
      function median_selector(lambda) result(selected)
         implicit none(type, external)
         complex(dp), intent(in) :: lambda(:)
         logical, allocatable :: selected(:)
         selected = abs(lambda) > median(abs(lambda))
      end function median_selector
   end procedure eigs_rdp

   module procedure eigs_csp
      ! Krylov subspace and Krylov subspace dimension.
      class(abstract_vector_csp), allocatable :: Xwrk(:)
      integer :: kdim_, kstart
      ! Hessenberg matrix.
      complex(sp), allocatable :: H(:, :)
      ! Working arrays for the eigenvectors and eigenvalues.
      complex(sp), allocatable :: eigvecs_wrk(:, :)
      complex(sp), allocatable :: eigvals_wrk(:)
      real(sp), allocatable :: residuals_wrk(:)
      ! Miscellaneous.
      character(len=*), parameter :: this_procedure = 'eigs_csp'
      integer :: nev, conv
      integer :: i, j, k, niter, krst, iostat
      real(sp) :: tol, x0_norm
      complex(sp) :: beta
      logical :: outpost
      character(len=256) :: msg

      if (time_lightkrylov()) call timer%start(this_procedure)
      ! Deals with optional parameters.
      nev = size(X)
      kdim_   = optval(kdim, 4*nev)
      tol     = optval(tolerance, rtol_sp)
      outpost = optval(write_intermediate, .true.)

      ! Allocate eigenvalues.
      allocate(eigvals(nev), source=zero_csp, stat=iostat, errmsg=msg)
      call check_allocation(iostat, msg, this_module, this_procedure)

      ! Allocate working variables.
      allocate(Xwrk(kdim_+1), source=X(1), stat=iostat, errmsg=msg)
      call check_allocation(iostat, msg, this_module, this_procedure)
      call zero_basis(Xwrk)

      if (present(x0)) then
         call copy(Xwrk(1), x0)
         x0_norm = x0%norm(); call Xwrk(1)%scal(one_csp/x0_norm)
      else
         call Xwrk(1)%rand(.true.)
      endif
      allocate(H(kdim_+1, kdim_), source=zero_csp, stat=iostat, errmsg=msg)
      call check_allocation(iostat, msg, this_module, this_procedure)
      allocate(eigvecs_wrk(kdim_, kdim_), source=zero_csp, stat=iostat, errmsg=msg)
      call check_allocation(iostat, msg, this_module, this_procedure)
      allocate(eigvals_wrk(kdim_), source=zero_csp, stat=iostat, errmsg=msg)
      call check_allocation(iostat, msg, this_module, this_procedure)
      allocate(residuals_wrk(kdim_), source=zero_rsp, stat=iostat, errmsg=msg)
      call check_allocation(iostat, msg, this_module, this_procedure)

      ! Ritz eigenpairs computation.
      H = 0.0_sp

      kstart = 1 ; conv = 0 ; niter = 0 ; krst = 1
      krylovschur: do while (conv < nev)

         arnoldi_factorization: do k = kstart, kdim_
               ! Arnoldi step.
               call arnoldi(A, Xwrk, H, info, kstart=k, kend=k, transpose=transpose)
               call check_info(info, 'arnoldi', this_module, this_procedure)

               ! Spectral decomposition of the k x k Hessenberg matrix.
               eigvals_wrk = 0.0_sp ; eigvecs_wrk = 0.0_sp
               if (time_lightkrylov()) call timer%start('eig')
               call eig(H(:k, :k), eigvecs_wrk(:k, :k), eigvals_wrk(:k))
               if (time_lightkrylov()) call timer%stop('eig')

               ! Compute residuals.
               beta = H(k+1, k)
               residuals_wrk(:k) = eigenvalue_residual_csp(beta, eigvecs_wrk(k,:k))

               ! Check convergence.
               niter = niter + 1
               conv = count(residuals_wrk(:k) < tol)
               write(msg,'(I0,A,I0,A,I0,A)') conv, '/', nev, ' eigenvalues converged after ', niter, &
                           & ' steps of the Arnoldi process.'
               call log_information(msg, this_module, this_procedure)
               if (outpost) call write_results_csp(eigs_output, eigvals_wrk(:k), residuals_wrk(:k), tol)
               if (conv >= nev) exit arnoldi_factorization
         enddo arnoldi_factorization

         write(msg,'(I0,A,I0,A,I0,A)') conv, '/', nev, ' eigenvalues converged after ', krst, &
                           & ' Krylov-Schur restarts of the Arnoldi process.'
         call log_information(msg, this_module, this_procedure)
         ! Krylov-Schur restarting procedure.
         krst  = krst + 1
         call krylov_schur(kstart, Xwrk, H, median_selector) ; kstart = kstart + 1

      end do krylovschur

      !--------------------------------
      !-----     POST-PROCESS     -----
      !--------------------------------

      block
      integer :: indices(kdim_)
      real(sp) :: abs_eigvals(kdim_)

      ! Re-compute eigenvalues and eigenvectors.
      k = min(k, kdim_)
      if (time_lightkrylov()) call timer%start('eig')
      call eig(H(:k, :k), eigvecs_wrk(:k, :k), eigvals_wrk(:k))
      if (time_lightkrylov()) call timer%stop('eig')
      ! Sort eigenvalues.
      abs_eigvals = abs(eigvals_wrk) ; call sort_index(abs_eigvals, indices, reverse=.true.)
      eigvals_wrk = eigvals_wrk(indices) ; eigvecs_wrk = eigvecs_wrk(:, indices)
      residuals_wrk = residuals_wrk(indices)

      ! Store converged eigenvalues.
      eigvals = eigvals_wrk(:nev) ; residuals = residuals_wrk(:nev)
      end block

      ! Construct eigenvectors.
      do i = 1, nev
         call X(i)%zero()
         do j = 1, k
               call X(i)%axpby(eigvecs_wrk(j, i), Xwrk(j), one_csp)
         enddo
      enddo

      info = niter
      if (time_lightkrylov()) call timer%stop(this_procedure)
   contains
      function median_selector(lambda) result(selected)
         implicit none(type, external)
         complex(sp), intent(in) :: lambda(:)
         logical, allocatable :: selected(:)
         selected = abs(lambda) > median(abs(lambda))
      end function median_selector
   end procedure eigs_csp

   module procedure eigs_cdp
      ! Krylov subspace and Krylov subspace dimension.
      class(abstract_vector_cdp), allocatable :: Xwrk(:)
      integer :: kdim_, kstart
      ! Hessenberg matrix.
      complex(dp), allocatable :: H(:, :)
      ! Working arrays for the eigenvectors and eigenvalues.
      complex(dp), allocatable :: eigvecs_wrk(:, :)
      complex(dp), allocatable :: eigvals_wrk(:)
      real(dp), allocatable :: residuals_wrk(:)
      ! Miscellaneous.
      character(len=*), parameter :: this_procedure = 'eigs_cdp'
      integer :: nev, conv
      integer :: i, j, k, niter, krst, iostat
      real(dp) :: tol, x0_norm
      complex(dp) :: beta
      logical :: outpost
      character(len=256) :: msg

      if (time_lightkrylov()) call timer%start(this_procedure)
      ! Deals with optional parameters.
      nev = size(X)
      kdim_   = optval(kdim, 4*nev)
      tol     = optval(tolerance, rtol_dp)
      outpost = optval(write_intermediate, .true.)

      ! Allocate eigenvalues.
      allocate(eigvals(nev), source=zero_cdp, stat=iostat, errmsg=msg)
      call check_allocation(iostat, msg, this_module, this_procedure)

      ! Allocate working variables.
      allocate(Xwrk(kdim_+1), source=X(1), stat=iostat, errmsg=msg)
      call check_allocation(iostat, msg, this_module, this_procedure)
      call zero_basis(Xwrk)

      if (present(x0)) then
         call copy(Xwrk(1), x0)
         x0_norm = x0%norm(); call Xwrk(1)%scal(one_cdp/x0_norm)
      else
         call Xwrk(1)%rand(.true.)
      endif
      allocate(H(kdim_+1, kdim_), source=zero_cdp, stat=iostat, errmsg=msg)
      call check_allocation(iostat, msg, this_module, this_procedure)
      allocate(eigvecs_wrk(kdim_, kdim_), source=zero_cdp, stat=iostat, errmsg=msg)
      call check_allocation(iostat, msg, this_module, this_procedure)
      allocate(eigvals_wrk(kdim_), source=zero_cdp, stat=iostat, errmsg=msg)
      call check_allocation(iostat, msg, this_module, this_procedure)
      allocate(residuals_wrk(kdim_), source=zero_rdp, stat=iostat, errmsg=msg)
      call check_allocation(iostat, msg, this_module, this_procedure)

      ! Ritz eigenpairs computation.
      H = 0.0_dp

      kstart = 1 ; conv = 0 ; niter = 0 ; krst = 1
      krylovschur: do while (conv < nev)

         arnoldi_factorization: do k = kstart, kdim_
               ! Arnoldi step.
               call arnoldi(A, Xwrk, H, info, kstart=k, kend=k, transpose=transpose)
               call check_info(info, 'arnoldi', this_module, this_procedure)

               ! Spectral decomposition of the k x k Hessenberg matrix.
               eigvals_wrk = 0.0_dp ; eigvecs_wrk = 0.0_dp
               if (time_lightkrylov()) call timer%start('eig')
               call eig(H(:k, :k), eigvecs_wrk(:k, :k), eigvals_wrk(:k))
               if (time_lightkrylov()) call timer%stop('eig')

               ! Compute residuals.
               beta = H(k+1, k)
               residuals_wrk(:k) = eigenvalue_residual_cdp(beta, eigvecs_wrk(k,:k))

               ! Check convergence.
               niter = niter + 1
               conv = count(residuals_wrk(:k) < tol)
               write(msg,'(I0,A,I0,A,I0,A)') conv, '/', nev, ' eigenvalues converged after ', niter, &
                           & ' steps of the Arnoldi process.'
               call log_information(msg, this_module, this_procedure)
               if (outpost) call write_results_cdp(eigs_output, eigvals_wrk(:k), residuals_wrk(:k), tol)
               if (conv >= nev) exit arnoldi_factorization
         enddo arnoldi_factorization

         write(msg,'(I0,A,I0,A,I0,A)') conv, '/', nev, ' eigenvalues converged after ', krst, &
                           & ' Krylov-Schur restarts of the Arnoldi process.'
         call log_information(msg, this_module, this_procedure)
         ! Krylov-Schur restarting procedure.
         krst  = krst + 1
         call krylov_schur(kstart, Xwrk, H, median_selector) ; kstart = kstart + 1

      end do krylovschur

      !--------------------------------
      !-----     POST-PROCESS     -----
      !--------------------------------

      block
      integer :: indices(kdim_)
      real(dp) :: abs_eigvals(kdim_)

      ! Re-compute eigenvalues and eigenvectors.
      k = min(k, kdim_)
      if (time_lightkrylov()) call timer%start('eig')
      call eig(H(:k, :k), eigvecs_wrk(:k, :k), eigvals_wrk(:k))
      if (time_lightkrylov()) call timer%stop('eig')
      ! Sort eigenvalues.
      abs_eigvals = abs(eigvals_wrk) ; call sort_index(abs_eigvals, indices, reverse=.true.)
      eigvals_wrk = eigvals_wrk(indices) ; eigvecs_wrk = eigvecs_wrk(:, indices)
      residuals_wrk = residuals_wrk(indices)

      ! Store converged eigenvalues.
      eigvals = eigvals_wrk(:nev) ; residuals = residuals_wrk(:nev)
      end block

      ! Construct eigenvectors.
      do i = 1, nev
         call X(i)%zero()
         do j = 1, k
               call X(i)%axpby(eigvecs_wrk(j, i), Xwrk(j), one_cdp)
         enddo
      enddo

      info = niter
      if (time_lightkrylov()) call timer%stop(this_procedure)
   contains
      function median_selector(lambda) result(selected)
         implicit none(type, external)
         complex(dp), intent(in) :: lambda(:)
         logical, allocatable :: selected(:)
         selected = abs(lambda) > median(abs(lambda))
      end function median_selector
   end procedure eigs_cdp

end submodule general_eigensolvers
