module lightkrylov_RationalKrylov
  !!@warning
  !! The content of this module still is experimental and is not considered ready for production.
  !!@endwarning

   ! LightKrylov modules.
   use lightkrylov_Utils
   use lightkrylov_AbstractVector
   use lightkrylov_LinearOperator
   use lightkrylov_IterativeSolvers

   ! Fortran standard library.
   use stdlib_optval, only: optval

   implicit none
   include "dtypes.h"

   private

   public :: rational_arnoldi_factorization

contains

   subroutine rational_arnoldi_factorization( &
      A, X, H, G, sigma, info, &
      solve, linsolve_opts, preconditioner, &
      kstart, kend, verbosity, tol, transpose)
      !! Implements the Rational Arnoldi Factorization for arbitrary square linear operators.
      !! Only the single pole case is considered with the exception of the last iteration (k=kend)
      !! where the pole is set to +infinity to define a proper Rational Arnoldi factorization.
      !!
      !! **Algorithmic Features**
      !!
      !! - Constructs a rational Krylov decomposition A*X*(I + HD) + v = XH + w where H is an upper
      !!   Hessenberg matrix, D a diagonal matrix constructed from the inverse shifts, and v and w
      !!   are residual vectors.
      !! - The orthogonal basis X for the Krylov subspace is iteratively constructed through the
      !!   Arnoldi procedure.
      !! - At the last iteration, the shift is chosen to be +infinity to ensure that we obtain a
      !!   well defined Arnoldi decomposition.
      !!
      !! **Advantages**
      !!
      !! - Rational Arnoldi produces better (rational) approximation for matrix-valued function f(A)
      !!   than its polynomial counterpart for identical dimension of the generated Krylov subspace.
      !!
      !! **Limitations**
      !!
      !! - By default, this implementation only handles the case of a single real repeated pole (with the
      !!   exception that the last pole is set to + infinity).
      !! - Computational cost is mostly dominated by solving (I - A/sigma) x = b. An efficient linear solver
      !!   passed via the solve argument is thus needed. At the current time, we do not support preconditioning yet.
      !! - Choosing an appropriate shift sigma has a huge influence on the capabilities of the rational
      !!   Arnoldi method to provide good approximations of matrix-valued functions f(A)*x. No good default choice
      !!   can be recommended as this is highly problem-dependent.
      !!
      !! **References**
      !!
      !! - S. Guttel. "Rational Krylov Methods for Operator Functions", PhD Thesis, 2010.
      class(abstract_linop), intent(in) :: A
      !! Linear operator to be factorized.
      class(abstract_vector), intent(inout) :: X(:)
      !! Resulting rational Krylov basis. On entry, `X(1)` needs to contain the starting Krylov vector.
      real(kind=wp), intent(inout) :: H(:, :), G(:, :)
      !! Upper Hessenberg matrices defining the rational factorization.
      real(kind=wp), intent(in) :: sigma
      !! Real-value shift for the shift inverse system.
      integer, intent(out) :: info
      !! Information flag.
      procedure(abstract_linear_solver) :: solve
      !! Linear solver to be used. It needs to conform to the `abstract_linear_solver` interface.
      class(abstract_opts), optional, intent(in) :: linsolve_opts
      !! Options to be passed to the linear solver.
      class(abstract_preconditioner), optional, intent(in) :: preconditioner
      !! Right-preconditioner used to solve the linear system.
      integer, optional, intent(in) :: kstart
      !! Starting index for the Krylov iteration.
      integer                       :: k_start
      integer, optional, intent(in) :: kend
      !! Final index for the Krylov iteration.
      integer                       :: k_end
      logical, optional, intent(in) :: verbosity
      !! Verbosity control.
      logical                       :: verbose
      real(kind=wp), optional, intent(in) :: tol
      !! Tolerance for invariant subspace detection.
      real(kind=wp)                       :: tolerance
      logical, optional, intent(in) :: transpose
      !! Determine whether \(\mathbf{A}\) (`.false.`) or \(\mathbf{A}^T\) (`.true.`) is being used.
      logical                       :: trans

      ! Internal variables.
      integer :: kdim
      integer :: i, j, k
      real(kind=wp) :: alpha, beta
      class(identity_linop), allocatable  :: Id
      class(axpby_linop), allocatable  :: S
      class(abstract_vector), allocatable :: wrk

      ! Krylov subspace dimension.
      kdim = size(X) - 1

      !------------------------------------
      !-----     Check dimensions     -----
      !------------------------------------
      call assert_shape(H, [kdim + 1, kdim], "rational_arnoldi", "H")
      call assert_shape(G, [kdim + 1, kdim], "rational_arnoldi", "G")

      !---------------------------------
      !-----     OPTIONAL ARGS     -----
      !---------------------------------

      k_start = optval(kstart, 1)
      k_end = optval(kend, kdim)
      verbose = optval(verbosity, .false.)
      tolerance = optval(tol, atol + rtol)
      trans = optval(transpose, .false.)

      !-------------------------------------------
      !-----     RATIONAL ARNOLDI METHOD     -----
      !-------------------------------------------

      ! Definition of the shifted operator to be inverted.
      Id = identity_linop()
      S = axpby_linop(Id, A, 1.0_wp, -1.0_wp/sigma) ! S = I - A/sigma

      ! Initialize working array.
      allocate (wrk, source=X(1))

      ! (Rational) Arnoldi factorization.
      arnoldi: do k = k_start, k_end
         ! Zero-out working array (sanity).
         call wrk%zero()

         ! Matrix vector product.
         if (trans) then
            call A%rmatvec(X(k), wrk)
         else
            call A%matvec(X(k), wrk)
         end if

         if (k < k_end) then
            call solve(S, wrk, X(k + 1), info, options=linsolve_opts, transpose=trans, preconditioner=preconditioner)
         else ! Last pole is set to +infinity.
            call Id%matvec(wrk, X(k + 1))
         end if

         ! Orthogonalize residual vector w.r.t. to previously computed Krylov vectors.
         do i = 1, k
            alpha = X(k + 1)%dot(X(i)); call X(k + 1)%axpby(1.0_wp, X(i), -alpha)
            ! Update Hessenberg matrix H.
            H(i, k) = alpha
         end do
         ! Perform full re-orthogonaliation (see instability of MGS process)
         do i = 1, k
            alpha = X(k + 1)%dot(X(i)); call X(k + 1)%axpby(1.0_wp, X(i), -alpha)
            ! Update Hessenberg matrix H.
            H(i, k) = H(i, k) + alpha
         end do

         ! Normalize residual and update the last row of the Hessenberg matrix.
         beta = X(k + 1)%norm(); H(k + 1, k) = beta

         ! Update the second Hessenberg matrix.
         if (k < k_end) then
            G(1:k, k) = H(1:k, k)/sigma; G(k, k) = G(k, k) + 1.0_wp; G(k + 1, k) = beta/sigma
         else
            ! Last pole is set to +infinity.
            G(k, k) = 1.0_wp; G(k + 1, k) = 0.0_wp
         end if

         ! Exit Arnoldi loop if needed.
         if (beta < tolerance) then
            ! Dimension of the computed invariant subspace.
            info = k
            ! Exit the loop.
            exit arnoldi
         else
            call X(k + 1)%scal(1.0_wp/beta)
         end if

      end do arnoldi

      return
   end subroutine rational_arnoldi_factorization

end module lightkrylov_RationalKrylov
