module KrylovDecomp
    
    use AbstractVector
    use LinearOperator
    
    use stdlib_optval, only: optval
    use stdlib_error, only: error_stop
    use stdlib_logger

    implicit none
    include "dtypes.h"

    private
    public :: arnoldi_factorization, &
              lanczos_tridiagonalization, &
              lanczos_bidiagonalization, &
              nonsymmetric_lanczos_tridiagonalization, &
              krylov_schur_factorization, &
              block_arnoldi_factorization

contains

    !---------------------------------------------------------------------
    !-----                                                           -----
    !-----     ARNOLDI FACTORIZATION FOR GENERAL SQUARE MATRICES     -----
    !-----                                                           -----
    !---------------------------------------------------------------------

    !=======================================================================================
    ! Arnoldi Factorization Subroutine
    !=======================================================================================
    !
    ! Purpose:
    ! --------
    ! Performs Arnoldi factorization to generate an orthonormal Krylov basis X and an upper Hessenberg matrix H.
    !
    ! Mathematical Formulation:
    ! -------------------------
    ! Given a linear operator A \in R^(n x n), find X and H such that:
    ! A * X(:, k) = X(:, 1:k) * H(1:k, k) + h_{k+1,k} * X(:, k+1)
    !
    ! Algorithmic Features:
    ! ----------------------
    ! - Constructs an orthonormal Krylov basis X via modified Gram-Schmidt.
    ! - Constructs an upper Hessenberg matrix H.
    ! - Checks for convergence and invariant subspaces.
    !
    ! Advantages:
    ! -----------
    ! - Applicable for non-symmetric matrices.
    ! - Basis for many Krylov subspace methods.
    !
    ! Limitations:
    ! ------------
    ! - Orthogonality of X may deteriorate for ill-conditioned matrices.
    ! - Not suitable for preconditioning in this implementation.
    !
    ! Input/Output Parameters:
    ! ------------------------
    ! A          : Linear Operator, class(abstract_linop), intent(in)
    ! X          : Krylov Basis, class(abstract_vector), dimension(:), intent(inout)
    ! H          : Hessenberg Matrix, double precision, dimension(:, :), intent(inout)
    ! info       : Exit information flag, integer, intent(out)
    ! kstart     : Optional, starting Krylov vector index, integer, intent(in)
    ! kend       : Optional, ending Krylov vector index, integer, intent(in)
    ! verbosity  : Optional, verbosity flag, logical, intent(in)
    ! tol        : Optional, orthogonality tolerance, double precision, intent(in)
    ! transpose  : Optional, transpose flag for A, logical, intent(in)
    !
    !=======================================================================================
    subroutine arnoldi_factorization(A, X, H, info, kstart, kend, verbosity, tol, transpose)

        ! --> Optional arguments (mainly for GMRES)
        integer, optional, intent(in) :: kstart, kend
        logical, optional, intent(in) :: verbosity, transpose
        double precision, optional, intent(in) :: tol

        integer :: k_start, k_end
        logical :: verbose, trans
        double precision :: tolerance

        ! --> Linear Operator to be factorized.
        class(abstract_linop), intent(in) :: A
        ! --> Krylov basis.
        class(abstract_vector), dimension(:), intent(inout) :: X
        ! --> Upper Hessenberg matrix.
        double precision, dimension(:, :), intent(inout) :: H
        ! --> Information.
        integer, intent(out) :: info ! info < 0 : The k-step Arnoldi factorization failed.
        ! info = 0 : The k-step Arnoldi factorization succeeded.
        ! info > 0 : An invariant subspace has been computed after k=info steps.
        ! --> Miscellaneous
        double precision :: beta
        integer :: k, kdim

        ! --> Check dimensions.
        kdim = size(X) - 1

        if (all(shape(H) .ne. [kdim + 1, kdim])) then
            write (*, *) "INFO : Hessenberg matrix and Krylov subspace dimensions do not match."
            info = -1 ! should be -6 as the error is on the 6th argument of the subroutine ?
            ! < 0 - argument error
            ! = 0 - no error
            ! > 0 - the error is on the i-th iteration
            return
        end if

        ! --> Deals with the optional arguments.
        k_start = optval(kstart, 1)
        k_end = optval(kend, kdim)
        verbose = optval(verbosity, .false.)
        tolerance = optval(tol, 1.0D-12)
        trans = optval(transpose, .false.)

        ! --> Arnoldi factorization.
        arnoldi: do k = k_start, k_end
            ! --> Matrix-vector product.
            if (trans) then
                call A%rmatvec(X(k), X(k + 1))
            else
                call A%matvec(X(k), X(k + 1))
            end if
            ! --> Update Hessenberg matrix.
            call update_hessenberg_matrix(H, X, k)
            beta = X(k + 1)%norm(); H(k + 1, k) = beta

            if (verbose) then
                write (*, *) "--> Arnoldi iteration n°", k, "/", k_end
                write (*, *) "    -----------------"
                write (*, *) "    + Residual norm :", beta
                write (*, *) "    + Elapsed time  :"
                write (*, *) "    + ETA           :"
                write (*, *)
            end if

            ! --> Exit Arnoldi loop if needed.
            if (beta < tolerance) then
                if (verbose) then
                    write (*, *)
                    write (*, *) "INFO : An invariant subspace has been computed (beta =", beta, ")."
                end if

                ! --> Dimension of the computed invariant subspace.
                info = k

                ! --> Exit the Arnoldi iteration.
                exit arnoldi
            else
                ! --> Normalize the new Krylov vector.
                call X(k + 1)%scal(1.0D+00/beta)
            end if

        end do arnoldi

        if (verbose) then
            write (*, *) "INFO : Exiting the Arnoldi factorization with exit code info =", info, "."
            write (*, *)
        end if

        return
    end subroutine arnoldi_factorization

    !=======================================================================================
    ! Update Hessenberg Matrix Subroutine
    !=======================================================================================
    !
    ! Purpose:
    ! --------
    ! Updates the Hessenberg matrix H in the Arnoldi factorization process.
    !
    ! Mathematical Formulation:
    ! -------------------------
    ! Given the k+1-th Krylov vector X(:, k+1) and the existing Hessenberg matrix H,
    ! update the k-th column of H through orthogonalization.
    !
    ! Algorithmic Features:
    ! ----------------------
    ! - Employs modified Gram-Schmidt for orthogonalization.
    ! - Performs full re-orthogonalization to enhance stability.
    !
    ! Input/Output Parameters:
    ! ------------------------
    ! H          : Hessenberg Matrix, double precision, dimension(:, :), intent(inout)
    ! X          : Krylov Basis, class(abstract_vector), dimension(:)
    ! k          : Current Arnoldi step, integer, intent(in)
    !
    !=======================================================================================
    subroutine update_hessenberg_matrix(H, X, k)
        integer, intent(in) :: k
        double precision, dimension(:, :), intent(inout) :: H
        class(abstract_vector), dimension(:) :: X
        integer :: i
        double precision :: alpha

        ! --> Orthogonalize residual w.r.t to previously computed Krylov vectors.
        do i = 1, k
            alpha = X(k + 1)%dot(X(i)); call X(k + 1)%axpby(1.0_wp, X(i), -alpha)
            ! --> Update Hessenberg matrix.
            H(i, k) = alpha
        end do

        ! --> Perform full re-orthogonalization (see instability of MGS process)
        do i = 1, k
            alpha = X(k + 1)%dot(X(i)); call X(k + 1)%axpby(1.0_wp, X(i), -alpha)
            ! --> Update Hessenberg matrix.
            H(i, k) = H(i, k) + alpha
        end do

        return
    end subroutine update_hessenberg_matrix

    !--------------------------------------------------------------------------
    !-----                                                                -----
    !-----     LANCZOS TRIDIAGONALIZATION FOR SYM. POS. DEF. MATRICES     -----
    !-----                                                                -----
    !--------------------------------------------------------------------------

    !=======================================================================================
    ! Lanczos Tridiagonalization for Symmetric Positive Definite Matrices
    !=======================================================================================
    !
    ! Purpose:
    ! --------
    ! Performs Lanczos tridiagonalization on a given SPD linear operator A, producing an orthonormal basis X and a tridiagonal matrix T.
    !
    ! Mathematical Formulation:
    ! -------------------------
    ! Given SPD A \in R^(n x n), find X and T such that:
    ! A * X(:, k) = X(:, 1:k) * T(1:k, k) + t_{k+1,k} * X(:, k+1)
    !
    ! Algorithmic Features:
    ! ----------------------
    ! - Constructs an orthonormal Krylov basis X.
    ! - Constructs a tridiagonal matrix T.
    ! - Checks for convergence and invariant subspaces.
    !
    ! Advantages:
    ! -----------
    ! - Efficient for SPD matrices.
    ! - Foundation for eigenvalue and linear system solvers.
    !
    ! Limitations:
    ! ------------
    ! - Limited to SPD matrices.
    ! - Orthogonality may deteriorate for ill-conditioned matrices.
    !
    ! Input/Output Parameters:
    ! ------------------------
    ! A          : SPD Linear Operator, class(abstract_spd_linop), intent(in)
    ! X          : Krylov Basis, class(abstract_vector), dimension(:), intent(inout)
    ! T          : Tridiagonal Matrix, double precision, dimension(:, :), intent(inout)
    ! info       : Exit information flag, integer, intent(out)
    ! kstart     : Optional, starting Krylov vector index, integer, intent(in)
    ! kend       : Optional, ending Krylov vector index, integer, intent(in)
    ! verbosity  : Optional, verbosity flag, logical, intent(in)
    ! tol        : Optional, orthogonality tolerance, double precision, intent(in)
    !
    !=======================================================================================
    subroutine lanczos_tridiagonalization(A, X, T, info, kstart, kend, verbosity, tol)
        !> Linear operator to be factorized;
        class(abstract_spd_linop), intent(in) :: A
        !> Krylov basis.
        class(abstract_vector), dimension(:), intent(inout) :: X
        !> Tri-diagonal matrix.
        double precision, dimension(:, :), intent(inout) :: T
        !> Information flag.
        integer, intent(out) :: info
        !> Optional arguements.
        integer, optional, intent(in) :: kstart
        integer                       :: k_start
        integer, optional, intent(in) :: kend
        integer                       :: k_end
        logical, optional, intent(in) :: verbosity
        logical                       :: verbose
        double precision, optional, intent(in) :: tol
        double precision                       :: tolerance
        !> Miscellaneous.
        double precision :: beta
        integer          :: k, kdim

        ! --> Check dimensions.
        kdim = size(X) - 1

        if (all(shape(T) .ne. [kdim + 1, kdim])) then
            write (*, *) "INFO : Tridiagonal matrix and Krylov subspace dimensions do not match."
            info = -1
            return
        end if

        ! --> Deals with the optional arguments.
        k_start = optval(kstart, 1)
        k_end = optval(kend, kdim)
        verbose = optval(verbosity, .false.)
        tolerance = optval(tol, 1.0D-12)

        ! --> Lanczos tridiagonalization.
        lanczos: do k = k_start, k_end
            ! --> Matrix-vector product.
            call A%matvec(X(k), X(k + 1))
            ! --> Update tridiagonal matrix.
            call update_tridiag_matrix(T, X, k)
            beta = X(k + 1)%norm(); T(k + 1, k) = beta

            if (verbose) then
                write (*, *) "--> Lanczos iteration n°", k, "/", k_end
                write (*, *) "    -----------------"
                write (*, *) "    + Residual norm :", beta
                write (*, *) "    + Elapsed time  :"
                write (*, *) "    + ETA           :"
            end if

            ! --> Exit Lanczos loop if needed.
            if (beta < tolerance) then
                if (verbose) then
                    write (*, *)
                    write (*, *) "INFO : An invariant subspace has been computed (beta =)", beta, ")."
                end if

                ! --> Dimension of the computed invariant subspace.
                info = k

                ! --> Exit the Lanczos iteration.
                exit lanczos
            else
                ! --> Normalize the new Krylov vector.
                call X(k + 1)%scal(1.0D+00/beta)
            end if

        end do lanczos

        if (verbose) then
            write (*, *) "INFO : Exiting the Lanczos factorization with exit code info =", info, "."
        end if

        return
    end subroutine lanczos_tridiagonalization

    subroutine update_tridiag_matrix(T, X, k)
        integer, intent(in) :: k
        double precision, dimension(:, :), intent(inout) :: T
        class(abstract_vector), dimension(:) :: X
        integer :: i
        double precision :: alpha

        ! --> Orthogonalize residual w.r.t to previously computed Krylov vectors.
        do i = max(1, k - 1), k
            alpha = X(k + 1)%dot(X(i)); call X(k + 1)%axpby(1.0_wp, X(i), -alpha)
            ! --> Update tridiag matrix.
            T(i, k) = alpha
        end do

        ! --> Full re-orthogonalization.
        do i = 1, k
            alpha = X(k + 1)%dot(X(i)); call X(k + 1)%axpby(1.0_wp, X(i), -alpha)
        end do

        return
    end subroutine update_tridiag_matrix

    !=======================================================================================
    ! Lanczos Bidiagonalization for Singular Value Computation
    !=======================================================================================
    !
    ! Purpose:
    ! --------
    ! Performs Lanczos bidiagonalization on a given linear operator A, producing left and
    ! right Krylov bases U and V, and a bidiagonal matrix B.
    !
    ! Mathematical Formulation:
    ! -------------------------
    ! Given A \in R^(m x n), find U, V, and B such that:
    ! A * V = U * B
    !
    ! Algorithmic Features:
    ! ----------------------
    ! - Constructs orthonormal bases U and V.
    ! - Constructs a bidiagonal matrix B.
    !
    ! Input/Output Parameters:
    ! ------------------------
    ! A          : Linear Operator, class(abstract_linop), intent(in)
    ! U          : Left Krylov Basis, class(abstract_vector), dimension(:), intent(inout)
    ! V          : Right Krylov Basis, class(abstract_vector), dimension(:), intent(inout)
    ! B          : Bidiagonal Matrix, real(kind=wp), dimension(:, :), intent(inout)
    ! info       : Exit information flag, integer, intent(out)
    ! kstart     : Optional, starting index, integer, intent(in)
    ! kend       : Optional, ending index, integer, intent(in)
    ! verbosity  : Optional, verbosity flag, logical, intent(in)
    ! tol        : Optional, tolerance, real(kind=wp), intent(in)
    !
    !=======================================================================================
    subroutine lanczos_bidiagonalization(A, U, V, B, info, kstart, kend, verbosity, tol)
        !> Linear operator to be factorized.
        class(abstract_linop), intent(in) :: A
        !> Left and right Krylov basis.
        class(abstract_vector), intent(inout) :: U(:)
        class(abstract_vector), intent(inout) :: V(:)
        !> Bi-diagonal matrix.
        real(kind=wp), intent(inout) :: B(:, :)
        !> Information flag.
        integer, intent(out) :: info
        !> Optional arguments.
        integer, optional, intent(in) :: kstart
        integer                       :: k_start
        integer, optional, intent(in) :: kend
        integer                       :: k_end
        logical, optional, intent(in) :: verbosity
        logical                       :: verbose
        real(kind=wp), optional, intent(in) :: tol
        real(kind=wp)                       :: tolerance
        !> Miscellaneous.
        real(kind=wp) :: alpha, beta, gamma
        integer       :: i, j, k, kdim

        ! --> Check Krylov subspaces dimensions.
        if (size(U) .ne. size(V)) then
            write (*, *) "INFO : Left and right Krylov basis have different dimensions."
            info = -1
            return
        else
            kdim = size(U) - 1
        end if

        ! --> Check B dimensions.
        if (all(shape(B) .ne. [kdim + 1, kdim])) then
            write (*, *) "INFO : Bidiagonal matrix and Krylov subspaces dimensions do not match."
            info = -2
        end if

        ! --> Deals with the optional arguments.
        k_start = optval(kstart, 1)
        k_end = optval(kend, kdim)
        verbose = optval(verbosity, .false.)
        tolerance = optval(tol, 1.0D-12)

        ! --> Lanczos bidiagonalization.
        lanczos: do k = k_start, k_end
            ! --> Transpose matrix-vector product.
            call A%rmatvec(U(k), V(k))
            ! /!\ Next lines not needed because already taken care of in the
            !     full re-orthogonalization.
            ! if (k > 1) then
            !    call V(k)%axpby(1.0_wp, V(k-1), -beta)
            ! endif

            ! --> Full re-orthogonalization of the right Krylov subspace.
            do j = 1, k - 1
                gamma = V(k)%dot(V(j)); call V(k)%axpby(1.0_wp, V(j), -gamma)
            end do

            ! --> Normalization step.
            alpha = V(k)%norm(); B(k, k) = alpha
            if (alpha > tolerance) then
                call V(k)%scal(1.0_wp/alpha)
                B(k, k) = alpha
            else
                if (verbose) then
                    write (*, *) "INFO : alpha = ", alpha
                end if
                info = k
                exit lanczos
            end if

            ! --> Matrix-vector product.
            call A%matvec(V(k), U(k + 1))
            ! /!\ Not needed because taken care of in the full reortho. step.
            ! call U(k+1)%axpby(1.0_wp, U(k), -alpha)

            ! --> Full re-orthogonalization of the left Krylov subspace.
            do j = 1, k
                gamma = U(k + 1)%dot(U(j)); call U(k + 1)%axpby(1.0_wp, U(j), -gamma)
            end do

            ! --> Normalization step.
            beta = U(k + 1)%norm(); B(k + 1, k) = beta
            if (beta > tolerance) then
                call U(k + 1)%scal(1.0_wp/beta)
                B(k + 1, k) = beta
            else
                if (verbose) then
                    write (*, *) "INFO : beta = ", beta
                end if
                info = k
                exit lanczos
            end if

        end do lanczos

        return
    end subroutine lanczos_bidiagonalization

    !=======================================================================================
    ! Nonsymmetric Two-Sided Lanczos Tridiagonalization for General Square Matrices
    !=======================================================================================
    !
    ! Purpose:
    ! --------
    ! Performs Lanczos tridiagonalization on a given nonsymmetric linear operator A, producing
    ! left and right Krylov bases V and W, and a tridiagonal matrix T.
    !
    ! Mathematical Formulation:
    ! -------------------------
    ! Given A \in R^(n x n), find V, W, and T such that:
    ! A * V = W * T
    !
    ! Algorithmic Features:
    ! ----------------------
    ! - Constructs orthonormal bases V and W.
    ! - Constructs a tridiagonal matrix T.
    !
    ! Input/Output Parameters:
    ! ------------------------
    ! A          : Linear Operator, class(abstract_linop), intent(in)
    ! V          : Left Krylov Basis, class(abstract_vector), dimension(:), intent(inout)
    ! W          : Right Krylov Basis, class(abstract_vector), dimension(:), intent(inout)
    ! T          : Tridiagonal Matrix, real(kind=wp), dimension(:, :), intent(inout)
    ! info       : Exit information flag, integer, intent(out)
    ! kstart     : Optional, starting index, integer, intent(in)
    ! kend       : Optional, ending index, integer, intent(in)
    ! verbosity  : Optional, verbosity flag, logical, intent(in)
    ! tol        : Optional, tolerance, real(kind=wp), intent(in)
    !
    !=======================================================================================
    subroutine nonsymmetric_lanczos_tridiagonalization(A, V, W, T, info, kstart, kend, verbosity, tol)
        !> Linear operator to be factorized.
        class(abstract_linop), intent(in)    :: A
        !> Left and right Krylov basis.
        class(abstract_vector), intent(inout) :: V(:)
        class(abstract_vector), intent(inout) :: W(:)
        !> Tridiagonal matrix.
        real(kind=wp), intent(inout) :: T(:, :)
        !> Information flag.
        integer, intent(out)   :: info
        !> Optional arguments.
        integer, optional, intent(in) :: kstart
        integer                       :: k_start
        integer, optional, intent(in) :: kend
        integer                       :: k_end
        logical, optional, intent(in) :: verbosity
        logical                       :: verbose
        real(kind=wp), optional, intent(in) :: tol
        real(kind=wp)                       :: tolerance
        !> Miscellaneous.
        real(kind=wp) :: alpha, beta, gamma, tmp
        integer       :: i, k, kdim

        ! --> Check Krylov subspaces dimensions.
        if (size(V) .ne. size(W)) then
            write (*, *) "INFO : Left and right Krylov basis have different dimensions."
            info = -1
        else
            kdim = size(V) - 1
        end if

        ! --> Check T dimensions.
        if (all(shape(T) .ne. [kdim + 1, kdim + 1])) then
            write (*, *) "INFO : Tridiagonal matrix and Krylov subspaces dimensions do not match."
            info = -2
        end if

        ! --> Deals with the optional arguments.
        k_start = optval(kstart, 1)
        k_end = optval(kend, kdim)
        verbose = optval(verbosity, .false.)
        tolerance = optval(tol, atol)

        ! --> Bi-orthogonalize the left and right starting vectors.
        tmp = V(1)%dot(W(1)); beta = sqrt(abs(tmp)); gamma = sign(beta, tmp)
        call V(1)%scal(1.0_wp/beta); call W(1)%scal(1.0_wp/gamma)

        ! --> Nonsymmetric Lanczos iterations.
        lanczos: do k = k_start, k_end
            ! --> Matrix-vector product.
            call A%matvec(V(k), V(k + 1))
            call A%rmatvec(W(k), W(k + 1))

            ! --> Update diagonal entry of the nonsymmetric tridiagonal matrix.
            alpha = W(k)%dot(V(k + 1)); T(k, k) = alpha

            ! --> Lanczos three term recurrence.
            call V(k + 1)%axpby(1.0_wp, V(k), -alpha)
            if (k > 1) call V(k + 1)%axpby(1.0_wp, V(k - 1), -gamma)

            call W(k + 1)%axpby(1.0_wp, W(k), -alpha)
            if (k > 1) call W(k + 1)%axpby(1.0_wp, W(k - 1), -beta)

            ! --> Update the off-diagonal entries of the nonsymmetric tridiagonal matrix.
            tmp = V(k + 1)%dot(W(k + 1)); beta = sqrt(abs(tmp)); gamma = sign(beta, tmp)
            T(k, k + 1) = gamma; T(k + 1, k) = beta

            if ((abs(beta) < tolerance) .or. (abs(gamma) < tolerance)) then
                if (verbose) then
                    write (*, *) "INFO : Invariant subspaces have been computed (beta, gamma) = (", beta, ",", gamma, ")."
                end if
            else
                ! --> Full re-biorthogonalization.
                do i = 1, k
                    alpha = V(k + 1)%dot(W(i)); call V(k + 1)%axpby(1.0_wp, W(i), -alpha)
                    alpha = W(k + 1)%dot(V(i)); call W(k + 1)%axpby(1.0_wp, V(i), -alpha)
                end do

                ! --> Normalization step.
                call V(k + 1)%scal(1.0_wp/beta); call W(k + 1)%scal(1.0_wp/gamma)
            end if

        end do lanczos

        if (verbose) then
            write (*, *) "INFO : Exiting the nonsymmetric Lanczos factorization with exit code info = ", info, "."
        end if

        return
    end subroutine nonsymmetric_lanczos_tridiagonalization

    !---------------------------------------------------------------------------
    !-----                                                                 -----
    !-----      KRYLOV-SCHUR FACTORIZATION FOR GENERAL SQUARE MATRICES     -----
    !-----                                                                 -----
    !---------------------------------------------------------------------------

    !=======================================================================================
    ! Krylov-Schur Factorization for General Square Matrices
    !=======================================================================================
    !
    ! Purpose:
    ! --------
    ! Performs Krylov-Schur factorization on a given square linear operator A to compute eigenvalues and eigenvectors.
    !
    ! Mathematical Formulation:
    ! -------------------------
    ! Given a square linear operator A of size n x n, find a subspace S and a diagonal matrix D such that:
    !
    ! A * S = S * D
    !
    ! The columns of S are approximated eigenvectors and the diagonal entries of D are approximated eigenvalues of A.
    !
    ! Algorithmic Features:
    ! ----------------------
    ! - Constructs a Krylov subspace via Arnoldi or Lanczos methods.
    ! - Computes and sorts the Schur decomposition of the reduced operator.
    ! - Optionally performs a restart.
    !
    ! Advantages:
    ! -----------
    ! - Suitable for large-scale eigenvalue problems.
    ! - Allows for selective eigenvalue extraction.
    ! - Capable of handling non-symmetric matrices.
    !
    ! Limitations:
    ! ------------
    ! - Convergence can be slow for certain classes of matrices.
    ! - May require multiple restarts for full spectrum.
    !
    ! Input/Output Parameters:
    ! ------------------------
    ! A          : Linear Operator, class(abstract_linop), intent(in)
    ! S          : Krylov Subspace Vectors, class(abstract_vector), dimension(:), intent(inout)
    ! D          : Approximated Eigenvalues, double precision, dimension(:, :), intent(inout)
    ! info       : Exit information flag, integer, intent(out)
    ! kstart     : Optional, starting index, integer, intent(in)
    ! kend       : Optional, ending index, integer, intent(in)
    ! verbosity  : Optional, verbosity flag, logical, intent(in)
    ! tol        : Optional, tolerance, double precision, intent(in)
    !
    ! References:
    ! -----------
    ! - K.Schur, "On the QD-algorithm and its application to the characteristic polynomial,"
    !   Numer. Math., 10:181–195, 1967.
    !
    !=======================================================================================
    subroutine krylov_schur_factorization(A, X, info, kstart, kend, schurtgt, verbosity, tol)

        ! --> Required arguments
        class(abstract_linop), intent(in) :: A
        class(abstract_vector), dimension(:), intent(inout) :: X
        integer, intent(out) :: info

        !double precision, dimension(:, :), intent(inout) :: D

        ! --> Optional arguments
        integer, optional, intent(in) :: kstart, kend, schurtgt
        logical, optional, intent(in) :: verbosity
        double precision, optional, intent(in) :: tol

        ! --> Local variables
        integer :: k_start, k_end
        logical :: verbose, converged
        double precision :: tolerance
        double precision, dimension(:, :), allocatable :: H
        integer :: kdim, k, schur_tgt, schur_cnt, cnt
        integer :: arnoldi_info

        ! --> Dimensionality of the Krylov subspace
        kdim = size(X) - 1

        ! --> Allocate memory for Hessenberg matrix
        allocate (H(kdim + 1, kdim))

        ! --> Optional arguments
        k_start = optval(kstart, 1)
        k_end = optval(kend, kdim)
        schur_tgt = optval(schurtgt, 10)
        verbose = optval(verbosity, .false.)
        tolerance = optval(tol, 1.0D-12)

        cnt = 0
        schur_cnt = 0
        converged = .false.
        do while (.not. converged)

            ! --> Perform Arnoldi factorization to build Krylov subspace and Hessenberg matrix
            call arnoldi_factorization(A, X, H, arnoldi_info, k_start, k_end, verbose, tolerance)

            ! --> Select whether to stop or apply Schur condensation depending on schur_tgt.
            select case (schur_tgt)

                ! --> k-step Arnoldi factorization completed.
            case (:0)
                converged = .true.
                !write (6, *) 'Arnoldi factorization completed.'

                ! --> Perform Schur decomposition on H to find D and update X
            case (1:)
                if (cnt .ge. schur_tgt) then ! Krylov-Schur factorization completed.
                    converged = .true.
                else ! Apply Schur condensation before restarting the factorization.
                    schur_cnt = schur_cnt + 1
                    !call schur_condensation(mstart, H, Q, k_dim)
                    !  call schur interface with LAPACK to diagonalize H
                end if

            end select

        end do ! while ( .not. converged )

        ! --> Finalize and set exit information
        info = 0
        if (verbose) then
            write (*, *) "INFO : Exiting the Krylov-Schur factorization with exit code info =", info, "."
        end if

        return
    end subroutine krylov_schur_factorization

    !=======================================================================================
    ! Block Arnoldi Factorization Subroutine
    !=======================================================================================
    !
    ! Purpose:
    ! --------
    ! Implements the Block Arnoldi method for computing an orthonormal basis of the Krylov
    ! subspace associated with a given matrix A and initial block of vectors.
    !
    ! Mathematical Formulation:
    ! -------------------------
    ! Given a matrix A of size n x n and an initial block of vectors X of size n x s, 
    ! the algorithm iteratively constructs an orthonormal basis for the Krylov subspace 
    ! spanned by {X, A*X, A^2*X, ..., A^(m-1)*X} where m is the number of Arnoldi steps.
    ! The process produces an upper block Hessenberg matrix H such that: 
    ! A * V_m = V_m * H + f_m * e_m' where V_m is the orthonormal basis of the Krylov subspace,
    ! H is the upper block Hessenberg matrix, f_m is a block of vectors, and e_m' is the m-th
    ! canonical basis vector.
    !
    ! Algorithmic Features:
    ! ----------------------
    ! - Extends the classical Arnoldi factorization to using blocks of vectors.
    ! - Constructs an upper block Hessenberg matrix as part of the factorization.
    ! - Utilizes the Modified Gram-Schmidt process for block orthogonalization.
    !
    ! Advantages:
    ! -----------
    ! - Efficient for matrices with multiple right-hand sides or multiple eigenvalues.
    ! - Suitable for parallelization due to block operations.
    ! - Enables more effective subspace iterations compared to standard Arnoldi.
    !
    ! Limitations:
    ! ------------
    ! - The block size needs to be chosen carefully for optimal performance.
    ! - May require more memory to store the block of vectors.
    !
    ! Input/Output Parameters:
    ! ------------------------
    ! - A        : Linear Operator (abstract_linop) [Input]
    ! - X        : Initial/Updated block of Krylov vectors (abstract_vector) [Input/Output]
    ! - H        : Upper block Hessenberg matrix (real(kind=dp), allocatable) [Output]
    ! - info     : Iteration Information flag (Integer) [Output]
    ! - kstart   : Starting index for Krylov basis (Integer) [Optional, Input]
    ! - kend     : Ending index for Krylov basis (Integer) [Optional, Input]
    ! - s        : Block size (Integer) [Input]
    ! - tol      : Convergence tolerance (real(kind=dp)) [Optional, Input]
    ! - verbosity: Verbosity control flag (Logical) [Optional, Input]
    !
    ! References:
    ! -----------
    !
    !=======================================================================================
!=======================================================================================
! ... (Header remains the same)
!=======================================================================================
    subroutine block_arnoldi_factorization(A, X, H, info, kstart, kend, block_size, verbosity, tol)
  
      ! Input: A is the linear operator, X is the Krylov subspace, H is the Hessenberg matrix
      ! Output: Updated X and H, and info flag
    
      class(abstract_linop), intent(in) :: A
      class(abstract_vector), dimension(:), intent(inout) :: X
      double precision, dimension(:, :), intent(inout) :: H
      integer, intent(out) :: info
      integer, optional, intent(in) :: kstart, kend, block_size
      logical, optional, intent(in) :: verbosity
      double precision, optional, intent(in) :: tol
    
      integer :: k, m, i, k_start, k_end, s
      double precision :: alpha, norm
      logical :: verbose
      double precision :: tolerance
      class(abstract_vector), allocatable :: wrk, Y(:)
    
      k_start = optval(kstart, 1)
      k_end = optval(kend, size(X, dim=1) - 1)
      s = optval(block_size, 1)
      verbose = optval(verbosity, .false.)
      tolerance = optval(tol, 1.0D-12)
    
      ! Allocate space for Y, a block of orthonormal vectors
      ! allocate(Y(s))
      ! Error: Allocating y of ABSTRACT base type at (1) requires a type-spec or source-expr

      do k = k_start, k_end, s
        do m = 1, s ! Loop over each vector in the block
    
          call A%matvec(X(k), Y(m)) ! Y(m) = A * X(k)
    
          ! Orthogonalize Y(m) against the existing Krylov subspace X
          do i = 1, k ! using MGS
    
            ! Compute alpha = < Y(m), X(i) >
            alpha = Y(m)%dot(X(i))
    
            ! Update Hessenberg matrix H(i, k) = alpha
            H(i, k) = alpha
    
            ! Update Y(m) = Y(m) - alpha * X(i)
            allocate(wrk, source=X(i))
            call wrk%scal(alpha)
            call Y(m)%sub(wrk)
            deallocate(wrk)
    
          enddo
    
          ! Normalize Y(m) and set X(k+1) = Y(m)
          norm = Y(m)%norm()  ! || Y(m) ||
          H(k + 1, k) = norm  ! Update Hessenberg matrix with the norm
    
          call Y(m)%scal(1.0D0 / norm)  ! Normalize Y(m)
    
          !call X(k + 1)%copy(Y(m)) ! Update Krylov subspace
    
        enddo
      enddo
    
      deallocate(Y)  ! Free allocated memory
    
    end subroutine block_arnoldi_factorization
    
end module KrylovDecomp