module KrylovDecomp
  use AbstractVector
  use LinearOperator
  use stdlib_optval, only : optval
  implicit none
  include "dtypes.h"

  private
  public :: arnoldi_factorization, &
       lanczos_tridiagonalization, &
       lanczos_bidiagonalization,  &
       nonsymmetric_lanczos_tridiagonalization, &
       rational_arnoldi_factorization

contains

  !---------------------------------------------------------------------
  !-----                                                           -----
  !-----     ARNOLDI FACTORIZATION FOR GENERAL SQUARE MATRICES     -----
  !-----                                                           -----
  !---------------------------------------------------------------------

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

    if ( all(shape(H) .ne. [kdim+1, kdim]) ) then
       write(*, *) "INFO : Hessenberg matrix and Krylov subspace dimensions do not match."
       info = -1
       return
    endif

    ! --> Deals with the optional arguments.
    k_start   = optval(kstart, 1)
    k_end     = optval(kend, kdim)
    verbose   = optval(verbosity, .false.)
    tolerance = optval(tol, 1.0D-12)
    trans     = optval(transpose, .false.)

    ! --> Arnoldi factorization.
    arnoldi: do k = k_start, k_end
       ! --> Matrix-vector product.
       if (trans) then
          call A%rmatvec(X(k), X(k+1))
       else
          call A%matvec(X(k), X(k+1))
       endif
       ! --> Update Hessenberg matrix.
       call update_hessenberg_matrix(H, X, k)
       beta = X(k+1)%norm() ; H(k+1, k) = beta

       if (verbose) then
          write(*, *) "--> Arnoldi iteration n°", k, "/", k_end
          write(*, *) "    -----------------"
          write(*, *) "    + Residual norm :", beta
          write(*, *) "    + Elapsed time  :"
          write(*, *) "    + ETA           :"
          write(*, *)
       endif

       ! --> Exit Arnoldi loop if needed.
       if (beta < tolerance) then
          if (verbose) then
             write(*, *)
             write(*, *) "INFO : An invariant subspace has been computed (beta =", beta, ")."
          endif

          ! --> Dimension of the computed invariant subspace.
          info = k

          ! --> Exit the Arnoldi iteration.
          exit arnoldi
       else
          ! --> Normalize the new Krylov vector.
          call X(k+1)%scal(1.0D+00 / beta)
       endif

    enddo arnoldi

    if(verbose) then
       write(*, *) "INFO : Exiting the Arnoldi factorization with exit code info =", info, "."
       write(*, *)
    endif

    return
  end subroutine arnoldi_factorization

  subroutine update_hessenberg_matrix(H, X, k)
    integer, intent(in) :: k
    double precision, dimension(:, :), intent(inout) :: H
    class(abstract_vector), dimension(:) :: X
    class(abstract_vector), allocatable :: wrk
    integer :: i
    double precision :: alpha

    ! --> Orthogonalize residual w.r.t to previously computed Krylov vectors.
    do i = 1, k
       alpha = X(k+1)%dot(X(i)) ; call X(k+1)%axpby(1.0_wp, X(i), -alpha)
       ! --> Update Hessenberg matrix.
       H(i, k) = alpha
    enddo

    ! --> Perform full re-orthogonalization (see instability of MGS process)
    do i = 1, k
       alpha = X(k+1)%dot(X(i)) ; call X(k+1)%axpby(1.0_wp, X(i), -alpha)
       ! --> Update Hessenberg matrix.
       H(i, k) = H(i, k) + alpha
    enddo

    return
  end subroutine update_hessenberg_matrix

  !--------------------------------------------------------------------------
  !-----                                                                -----
  !-----     LANCZOS TRIDIAGONALIZATION FOR SYM. POS. DEF. MATRICES     -----
  !-----                                                                -----
  !--------------------------------------------------------------------------

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
    integer          :: i, j

    ! --> Check dimensions.
    kdim = size(X) - 1

    if (all(shape(T) .ne. [kdim+1, kdim])) then
       write(*, *) "INFO : Tridiagonal matrix and Krylov subspace dimensions do not match."
       info = -1
       return
    endif

    ! --> Deals with the optional arguments.
    k_start   = optval(kstart, 1)
    k_end     = optval(kend, kdim)
    verbose   = optval(verbosity, .false.)
    tolerance = optval(tol, 1.0D-12)

    ! --> Lanczos tridiagonalization.
    lanczos: do k = k_start, k_end
       ! --> Matrix-vector product.
       call A%matvec(X(k), X(k+1))
       ! --> Update tridiagonal matrix.
       call update_tridiag_matrix(T, X, k)
       beta = X(k+1)%norm() ; T(k+1, k) = beta

       if (verbose) then
          write(*, *) "--> Lanczos iteration n°", k, "/", k_end
          write(*, *) "    -----------------"
          write(*, *) "    + Residual norm :", beta
          write(*, *) "    + Elapsed time  :"
          write(*, *) "    + ETA           :"
       endif

       ! --> Exit Lanczos loop if needed.
       if (beta < tolerance) then
          if (verbose) then
             write(*, *)
             write(*, *) "INFO : An invariant subspace has been computed (beta =)", beta, ")."
          endif

          ! --> Dimension of the computed invariant subspace.
          info = k

          ! --> Exit the Lanczos iteration.
          exit lanczos
       else
          ! --> Normalize the new Krylov vector.
          call X(k+1)%scal(1.0D+00 / beta)
       endif

    enddo lanczos

    if (verbose) then
       write(*, *) "INFO : Exiting the Lanczos factorization with exit code info =", info, "."
    endif

    return
  end subroutine lanczos_tridiagonalization

  subroutine update_tridiag_matrix(T, X, k)
    integer, intent(in) :: k
    double precision, dimension(:, :), intent(inout) :: T
    class(abstract_vector), dimension(:) :: X
    class(abstract_vector), allocatable :: wrk
    integer :: i
    double precision :: alpha

    ! --> Orthogonalize residual w.r.t to previously computed Krylov vectors.
    do i = max(1, k-1), k
       alpha = X(k+1)%dot(X(i)) ; call X(k+1)%axpby(1.0_wp, X(i), -alpha)
       ! --> Update tridiag matrix.
       T(i, k) = alpha
    enddo

    ! --> Full re-orthogonalization.
    do i = 1, k
       alpha = X(k+1)%dot(X(i)) ; call X(k+1)%axpby(1.0_wp, X(i), -alpha)
    enddo

    return
  end subroutine update_tridiag_matrix

  !----------------------------------------------------------------------------
  !-----                                                                  -----
  !-----     LANCZOS BIDIAGONALIZATION FOR SINGULAR VALUE COMPUTATION     -----
  !-----                                                                  -----
  !----------------------------------------------------------------------------

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
    class(abstract_vector), allocatable :: wrk

    ! --> Check Krylov subspaces dimensions.
    if (size(U) .ne. size(V)) then
       write(*, *) "INFO : Left and right Krylov basis have different dimensions."
       info = -1
       return
    else
       kdim = size(U) - 1
    endif

    ! --> Check B dimensions.
    if (all(shape(B) .ne. [kdim+1, kdim])) then
       write(*, *) "INFO : Bidiagonal matrix and Krylov subspaces dimensions do not match."
       info = -2
    endif

    ! --> Deals with the optional arguments.
    k_start = optval(kstart, 1)
    k_end   = optval(kend, kdim)
    verbose = optval(verbosity, .false.)
    tolerance = optval(tol, 1.0D-12)

    ! --> Lanczos bidiagonalization.
    lanczos : do k = k_start, k_end
       ! --> Transpose matrix-vector product.
       call A%rmatvec(U(k), V(k))
       ! /!\ NOTE : Next lines not needed because already taken care of in the
       !     full re-orthogonalization.
       ! if (k > 1) then
       !    call V(k)%axpby(1.0_wp, V(k-1), -beta)
       ! endif

       ! --> Full re-orthogonalization of the right Krylov subspace.
       do j = 1, k-1
          gamma = V(k)%dot(V(j)) ; call V(k)%axpby(1.0_wp, V(j), -gamma)
       enddo

       ! --> Normalization step.
       alpha = V(k)%norm() ; B(k, k) = alpha
       if (alpha > tolerance) then
          call V(k)%scal(1.0_wp / alpha)
          B(k, k) = alpha
       else
          if (verbose) then
             write(*, *) "INFO : alpha = ", alpha
          endif
          info = k
          exit lanczos
       endif

       ! --> Matrix-vector product.
       call A%matvec(V(k), U(k+1))
       ! /!\ NOTE : Not needed because taken care of in the full reortho. step.
       ! call U(k+1)%axpby(1.0_wp, U(k), -alpha)

       ! --> Full re-orthogonalization of the left Krylov subspace.
       do j = 1, k
          gamma = U(k+1)%dot(U(j)) ; call U(k+1)%axpby(1.0_wp, U(j), -gamma)
       enddo

       ! --> Normalization step.
       beta = U(k+1)%norm() ; B(k+1, k) = beta
       if (beta > tolerance) then
          call U(k+1)%scal(1.0_wp / beta)
          B(k+1, k) = beta
       else
          if (verbose) then
             write(*, *) "INFO : beta = ", beta
          endif
          info = k
          exit lanczos
       endif

    enddo lanczos

    return
  end subroutine lanczos_bidiagonalization

  !------------------------------------------------------------------------------------
  !-----                                                                          -----
  !-----     TWO-SIDED LANCZOS TRIDIAGONALIZATION FOR GENERAL SQUARE MATRICES     -----
  !-----                                                                          -----
  !------------------------------------------------------------------------------------

  subroutine nonsymmetric_lanczos_tridiagonalization(A, V, W, T, info, kstart, kend, verbosity, tol)
    !> Linear operator to be factorized.
    class(abstract_linop) , intent(in)    :: A
    !> Left and right Krylov basis.
    class(abstract_vector), intent(inout) :: V(:)
    class(abstract_vector), intent(inout) :: W(:)
    !> Tridiagonal matrix.
    real(kind=wp)         , intent(inout) :: T(:, :)
    !> Information flag.
    integer               , intent(out)   :: info
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
    integer       :: i, j, k, kdim

    ! --> Check Krylov subspaces dimensions.
    if (size(V) .ne. size(W)) then
       write(*, *) "INFO : Left and right Krylov basis have different dimensions."
       info = -1
    else
       kdim = size(V) - 1
    endif

    ! --> Check T dimensions.
    if (all(shape(T) .ne. [kdim+1, kdim+1])) then
       write(*, *) "INFO : Tridiagonal matrix and Krylov subspaces dimensions do not match."
       info = -2
    endif

    ! --> Deals with the optional arguments.
    k_start = optval(kstart, 1)
    k_end   = optval(kend, kdim)
    verbose = optval(verbosity, .false.)
    tolerance = optval(tol, atol)

    ! --> Bi-orthogonalize the left and right starting vectors.
    tmp = V(1)%dot(W(1)) ; beta = sqrt(abs(tmp)) ; gamma = sign(beta, tmp)
    call V(1)%scal(1.0_wp / beta) ; call W(1)%scal(1.0_wp / gamma)

    ! --> Nonsymmetric Lanczos iterations.
    lanczos : do k = k_start, k_end
       ! --> Matrix-vector product.
       call A%matvec(V(k), V(k+1))
       call A%rmatvec(W(k), W(k+1))

       ! --> Update diagonal entry of the nonsymmetric tridiagonal matrix.
       alpha = W(k)%dot(V(k+1)) ; T(k, k) = alpha

       ! --> Lanczos three term recurrence.
       call V(k+1)%axpby(1.0_wp, V(k), -alpha)
       if (k > 1) call V(k+1)%axpby(1.0_wp, V(k-1), -gamma)

       call W(k+1)%axpby(1.0_wp, W(k), -alpha)
       if (k > 1) call W(k+1)%axpby(1.0_wp, W(k-1), -beta)

       ! --> Update the off-diagonal entries of the nonsymmetric tridiagonal matrix.
       tmp = V(k+1)%dot(W(k+1)) ; beta = sqrt(abs(tmp)) ; gamma = sign(beta, tmp)
       T(k, k+1) = gamma ; T(k+1, k) = beta

       if ((abs(beta) < tolerance) .or. (abs(gamma) < tolerance)) then
          if (verbose) then
             write(*, *) "INFO : Invariant subspaces have been computed (beta, gamma) = (", beta, ",", gamma, ")."
          endif
       else
          ! --> Full re-biorthogonalization.
          do i = 1, k
             alpha = V(k+1)%dot(W(i)) ; call V(k+1)%axpby(1.0_wp, W(i), -alpha)
             alpha = W(k+1)%dot(V(i)) ; call W(k+1)%axpby(1.0_wp, V(i), -alpha)
          enddo

          ! --> Normalization step.
          call V(k+1)%scal(1.0_wp / beta) ; call W(k+1)%scal(1.0_wp / gamma)
       endif
       
    enddo lanczos

    if (verbose) then
       write(*, *) "INFO : Exiting the nonsymmetric Lanczos factorization with exit code info = ", info, "."
    endif

    return
  end subroutine nonsymmetric_lanczos_tridiagonalization

  !======================================================================================  
  ! Rational Arnoldi Factorization with single pole Subroutine
  !======================================================================================
  !
  ! Purpose:
  ! -------
  ! Implements the Rational Arnoldi Factorization for arbitrary square linear operators.
  ! Only the single pole case is considered with the exception of the last iteration (k=kend)
  ! where the pole is set to +infinity to define a proper Rational Arnoldi factorization.
  !
  ! Algorithmic Features:
  ! ---------------------
  ! - Constructs a rational Krylov decomposition A*X*(I + HD) + v = XH + w where H is an upper
  !   Hessenberg matrix, D a diagonal matrix constructed from the inverse shifts, and v and w
  !   are residual vectors.
  ! - The orthogonal basis X for the Krylov subspace is iteratively constructed through the
  !   Arnoldi procedure.
  ! - At the last iteration, the shift is chosen to be +infinity to ensure that we obtain a
  !   well defined Arnoldi decomposition.
  !
  ! Advantages:
  ! -----------
  ! - Rational Arnoldi produces better (rational) approximation for matrix-valued function f(A)
  !   than its polynomial counterpart for identical dimension of the generated Krylov subspace.
  !
  ! Limitations:
  ! ------------
  ! - By default, this implementation only handles the case of a single real repeated pole (with the
  !   exception that the last pole is set to + infinity).
  ! - Computational cost is mostly dominated by solving (I - A/sigma) x = b. An efficient linear solver
  !   passed via the solve argument is thus needed. At the current time, we do not support preconditioning yet.
  ! - Choosing an appropriate shift sigma has a huge influence on the capabilities of the rational
  !   Arnoldi method to provide good approximations of matrix-valued functions f(A)*x. No good default choice
  !   can be recommended as this is highly problem-dependent.
  !
  ! Input/Output Parameters:
  ! ------------------------
  ! - A         : Linear operator [Input]
  ! - X         : Krylov subspace, assumed to be initialized                [Input/Output]
  ! - H         : Upper Hessenberg matrix                                   [Input/Output]
  ! - G         : Upper Hessenberg matrix                                   [Input/Output]
  ! - sigma     : real shift                                                [Input]
  ! - info      : Information flag                                          [Output]
  ! - solve     : Linear solver for (I - A/sigma)*x=b                       [Input]
  ! - kstart    : Starting index for the Krylov iteration                   [Optional, Input]
  ! - kend      : Last index for the Krylov iteration                       [Optional, Input]
  ! - verbosity : Verbosity control flag                                    [Optional, Input]
  ! - tol       : Tolerance for convergence                                 [Optional, Input]
  ! - transpose : Logical flag to control whether A or transpose(A) is used [Optional, Input]
  ! References:
  ! -----------
  ! - S. Guttel. "Rational Krylov Methods for Operator Functions", PhD Thesis, 2010.
  !
  !======================================================================================
  subroutine rational_arnoldi_factorization(A, X, H, G, sigma, info, solve, kstart, kend, verbosity, tol, transpose)
    !> Linear operator to be factorized
    class(abstract_linop), intent(in) :: A
    !> Krylov basis.
    class(abstract_vector), intent(inout) :: X(:)
    !> Upper Hessenberg matrices.
    real(kind=wp), intent(inout) :: H(:, :), G(:, :)
    !> Real shift for the rational approximation.
    real(kind=wp), intent(in) :: sigma
    !> Information flag.
    integer, intent(out) :: info
    !> Linear solver.
    interface
       subroutine solve(A, b, x, info, kdim, maxiter, tol, verbosity, transpose)
         import abstract_linop, abstract_vector, wp
         class(abstract_linop)  , intent(in)    :: A
         class(abstract_vector) , intent(in)    :: b
         class(abstract_vector) , intent(inout) :: x
         integer                , intent(out)   :: info
         integer, optional      , intent(in)    :: kdim
         integer, optional      , intent(in)    :: maxiter
         real(kind=wp), optional, intent(in)    :: tol
         logical, optional      , intent(in)    :: verbosity
         logical, optional      , intent(in)    :: transpose
       end subroutine solve
    end interface
    !> Optional arguments.
    integer, optional, intent(in) :: kstart
    integer                       :: k_start
    integer, optional, intent(in) :: kend
    integer                       :: k_end
    logical, optional, intent(in) :: verbosity
    logical                       :: verbose
    real(kind=wp), optional, intent(in) :: tol
    real(kind=wp)                       :: tolerance
    logical, optional, intent(in) :: transpose
    logical                       :: trans

    !> Internal variables.
    integer :: kdim
    integer :: i, j, k
    real(kind=wp) :: alpha, beta
    class(identity_linop), allocatable  :: Id
    class(axpby_linop)   , allocatable  :: S
    class(abstract_vector), allocatable :: wrk

    !> Krylov subspace dimension.
    kdim = size(X) - 1

    !------------------------------------
    !-----     Check dimensions     -----
    !------------------------------------

    if ((size(X) .ne. size(H, 1)) .or. (size(X)-1 .ne. size(H, 2))) then
       info = -3
    endif

    if ((size(X) .ne. size(G, 1)) .or. (size(X)-1 .ne. size(G, 2))) then
       info = -4
    endif

    !---------------------------------
    !-----     OPTIONAL ARGS     -----
    !---------------------------------

    k_start   = optval(kstart, 1)
    k_end     = optval(kend, kdim)
    verbose   = optval(verbosity, .false.)
    tolerance = optval(tol, atol + rtol)
    trans     = optval(transpose, .false.)

    !-------------------------------------------
    !-----     RATIONAL ARNOLDI METHOD     -----
    !-------------------------------------------

    !> Definition of the shifted operator to be inverted.
    Id = identity_linop()
    S  = axpby_linop(Id, A, 1.0_wp, -1.0_wp/sigma) ! S = I - A/sigma

    !> Initialize working array.
    allocate(wrk, source=X(1))

    !> (Rational) Arnoldi factorization.
    arnoldi: do k = k_start, k_end
       ! --> Zero-out working array (sanity).
       call wrk%zero()

       ! --> Matrix vector product.
       if (trans) then
          call A%rmatvec(X(k), wrk)
       else
          call A%matvec(X(k), wrk)
       endif

       if (k < k_end) then
          call solve(S, wrk, X(k+1), info, kdim, transpose=trans)
       else !> Last pole is set to +infinity.
         call Id%matvec(wrk, X(k+1))
       endif

       ! --> Orthogonalize residual vector w.r.t. to previously computed Krylov vectors.
       do i = 1, k
          alpha = X(k+1)%dot(X(i)) ; call X(k+1)%axpby(1.0_wp, X(i), -alpha)
          ! --> Update Hessenberg matrix H.
          H(i, k) = alpha
       enddo
       ! --> Perform full re-orthogonaliation (see instability of MGS process)
       do i = 1, k
          alpha = X(k+1)%dot(X(i)) ; call X(k+1)%axpby(1.0_wp, X(i), -alpha)
          ! --> Update Hessenberg matrix H.
          H(i, k) = H(i, k) + alpha
       enddo

       ! --> Normalize residual and update the last row of the Hessenberg matrix.
       beta = X(k+1)%norm() ; H(k+1, k) = beta

       ! --> Update the second Hessenberg matrix.
       if (k < k_end) then
          G(1:k, k) = H(1:k, k) / sigma ; G(k, k) = G(k, k) + 1.0_wp ; G(k+1, k) = beta / sigma
       else !> Last pole is set to +infinity.
          G(k, k) = 1.0_wp ; G(k+1, k) = 0.0_wp
       endif

       ! --> Exit Arnoldi loop if needed.
       if (beta < tolerance) then
          !> Logging information.
          !> Dimension of the computed invariant subspace.
          info = k
          !> Exit the loop.
          exit arnoldi
       else
          call X(k+1)%scal(1.0_wp / beta)
       endif

    enddo arnoldi

    return
  end subroutine rational_arnoldi_factorization

end module KrylovDecomp
