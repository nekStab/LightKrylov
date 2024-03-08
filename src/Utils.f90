module lightkrylov_utils
   implicit none
   include "dtypes.h"

   private
   !> General-purpose utilities.
   public :: assert_shape
   !> Linear Algebra Utilities.
   public :: inv, svd, eig, eigh, lstsq

   !-------------------------------------------------------
   !-----                                             -----
   !-----     OPTS TYPE OBJECT FOR LINEAR SOLVERS     -----
   !-----                                             -----
   !-------------------------------------------------------

   ! --> Base type.
   type, abstract, public :: abstract_opts
   end type abstract_opts

   ! --> GMRES options.
   type, extends(abstract_opts), public :: gmres_opts
      !> Default dimension of the Krylov subspace.
      integer :: kdim = 30
      !> Default maximum number of gmres restarts.
      integer :: maxiter = 10
      !> Default tolerance.
      real(kind=wp) :: atol = atol
      real(kind=wp) :: rtol = rtol
      !> Default verbosity control.
      logical :: verbose = .false.
   end type gmres_opts

   ! --> BICGSTAB options.
   type, extends(abstract_opts), public :: bicgstab_opts
      !> Default maximum number of iterations.
      integer :: maxiter = 100
      !> Default tolerance.
      real(kind=wp) :: atol = atol
      real(kind=wp) :: rtol = rtol
      !> Default verbosity control.
      logical :: verbose = .false.
   end type bicgstab_opts

   ! --> Conjugate Gradient options.
   type, extends(abstract_opts), public :: cg_opts
      !> Default maximum number of iterations.
      integer :: maxiter = 100
      !> Default tolerance.
      real(kind=wp) :: rtol = rtol
      real(kind=wp) :: atol = atol
      !> Default verbosity control.
      logical :: verbose = .false.
   end type cg_opts

   !------------------------------
   !-----     INTERFACES     -----
   !------------------------------

   !> Check dimensions of an array.
   interface assert_shape
      module procedure dassert_shape
      module procedure zassert_shape
   end interface assert_shape

   !> Compute inv(A).
   interface inv
      module procedure dinv
      module procedure zinv
   end interface inv

   !> Compute U, S, V = svd(A)
   interface svd
      module procedure dsvd
      module procedure zsvd
   end interface svd

   !> Compute EVD(A).
   interface eig
      module procedure deig
   end interface eig

   !> Compute EVD(A) with A sym. pos. def / Hermitian.
   interface eigh
      module procedure deigh
   end interface eigh

contains

   !-------------------------------------
   !-----                           -----
   !-----     VARIOUS UTILITIES     -----
   !-----                           -----
   !-------------------------------------

   subroutine stop_error(msg)
      !> Error message.
      character(len=*), intent(in) :: msg

      write (*, *) msg
      stop 1

      return
   end subroutine stop_error

   subroutine dassert_shape(A, size, routine, matname)
      !> Input matrix and expected dimensions.
      real(kind=wp), intent(in) :: A(:, :)
      integer, intent(in) :: size(:)
      character(len=*), intent(in) :: routine, matname

      if (any(shape(A) /= size)) then
         write (*, *) "In routine "//routine//" matrix "//matname//" has illegal shape ", shape(A)
         write (*, *) "Expected shape is ", size
         call stop_error("Aborting due to illegal matrix operation.")
      end if

      return
   end subroutine dassert_shape

   subroutine zassert_shape(A, size, routine, matname)
      !> Input matrix and expected dimensions.
      complex(kind=wp), intent(in) :: A(:, :)
      integer, intent(in) :: size(:)
      character(len=*), intent(in) :: routine, matname

      if (any(shape(A) /= size)) then
         write (*, *) "In routine "//routine//" matrix "//matname//" has illegal shape ", shape(A)
         write (*, *) "Expected shape is ", size
         call stop_error("Aborting due to illegal matrix operation.")
      end if

      return
   end subroutine zassert_shape

   !-------------------------------------------
   !-----                                 -----
   !-----     LAPACK MATRIX INVERSION     -----
   !-----                                 -----
   !-------------------------------------------

   subroutine dinv(A)
      !> Matrix to invert (in-place)
      real(kind=wp), intent(inout) :: A(:, :)
      !> Lapack-related.
      integer :: n, info
      real(kind=wp) :: work(size(A, 1))
      integer       :: ipiv(size(A, 1))

      !> Compute A = LU (in-place)
      n = size(A, 1); call assert_shape(A, [n, n], "inv", "A")
      call dgetrf(n, n, A, n, ipiv, info)
      if (info /= 0) then
         write (*, *) "DGETRF returned info = ", info
         if (info < 0) then
            write (*, *) "The ", -info, "-th argument has an illegal value."
         else
            write (*, *) "U(", info, ",", info, ") is exactly zero. The factorization"
            write (*, *) "has been completed but the factor U is exactly singular."
            write (*, *) "Division by zero will occur if used to solve Ax = b."
         end if
         call stop_error("inv: DGETRF error.")
      end if

      !> Compute inv(A).
      call dgetri(n, A, n, ipiv, work, n, info)
      if (info /= 0) then
         write (*, *) "DGETRI returned info =", info
         if (info < 0) then
            write (*, *) "The ", -info, "-th argument has an illegal value."
         else
            write (*, *) "U(", info, ",", info, ") is exactly zero."
            write (*, *) "The matrix is singular and its inverse cannot be computed."
         end if
         call stop_error("inv: DGETRI error.")
      end if

      return
   end subroutine dinv

   subroutine zinv(A)
      !> Matrix to invert (in-place)
      complex(kind=wp), intent(inout) :: A(:, :)
      !> Lapack-related.
      integer :: n, info
      complex(kind=wp) :: work(size(A, 1))
      integer          :: ipiv(size(A, 1))

      !> Compute A = LU (in-place).
      n = size(A, 1); call assert_shape(A, [n, n], "inv", "A")
      call zgetrf(n, n, A, n, ipiv, info)
      if (info /= 0) then
         write (*, *) "ZGETRF returned info = ", info
         if (info < 0) then
            write (*, *) "The ", -info, "-th argument has an illegal value."
         else
            write (*, *) "U(", info, ",", info, ") is exactly zero. The factorization"
            write (*, *) "has been completed but the factor U is exactly singular."
            write (*, *) "Division by zero will occur if used to solve Ax = b."
         end if
         call stop_error("inv: ZGETRF error.")
      end if

      !> Compute inv(A).
      call zgetri(n, A, n, ipiv, work, n, info)
      if (info /= 0) then
         write (*, *) "ZGETRI returned info =", info
         if (info < 0) then
            write (*, *) "The ", -info, "-th argument has an illegal value."
         else
            write (*, *) "U(", info, ",", info, ") is exactly zero."
            write (*, *) "The matrix is singular and its inverse cannot be computed."
         end if
         call stop_error("inv: ZGETRI error.")
      end if

      return
   end subroutine zinv

   !------------------------------------------
   !-----                                -----
   !-----     LAPACK SVD COMPUTATION     -----
   !-----                                -----
   !------------------------------------------

   subroutine dsvd(A, U, S, V)
      !> Matrix to be factorized
      real(kind=wp), intent(in)  :: A(:, :)
      !> Left singular vectors.
      real(kind=wp), intent(out) :: U(:, :)
      !> Singular values.
      real(kind=wp), intent(out) :: S(:)
      !> Right singular vectors.
      real(kind=wp), intent(out) :: V(:, :)

      !> Lapack-related.
      character :: jobu = "S", jobvt = "S"
      integer   :: m, n, lda, ldu, ldvt, lwork, info
      real(kind=wp), allocatable :: work(:)
      real(kind=wp) :: A_tilde(size(A, 1), size(A, 2)), vt(min(size(A, 1), size(A, 2)), size(A, 2))

      !> Setup variables.
      m = size(A, 1); n = size(A, 2)
      lda = m; ldu = m; ldvt = n
      lwork = max(1, 3*min(m, n), 5*min(m, n)); allocate (work(lwork))

      !> Shape assertions.
      call assert_shape(U, [m, m], "svd", "U")
      call assert_shape(V, [n, n], "svd", "V")

      !> SVD computation.
      a_tilde = a
      call dgesvd(jobu, jobvt, m, n, a_tilde, lda, s, u, ldu, vt, ldvt, work, lwork, info)
      if (info /= 0) then
         write (*, *) "DGESVD returned info = ", info
         if (info < 0) then
            write (*, *) "The ", -info, "-th argument has an illegal value."
         else
            write (*, *) "DBSQR did not converge. There are ", info, "superdiagonals"
            write (*, *) "of an intermediate bidiagonal matrix form B which did not"
            write (*, *) "converge to zero. See Lapack documentation for more details."
         end if
         call stop_error("svd: dgesvd error")
      end if

      !> Return the transpose of V.
      v = transpose(vt)

      return
   end subroutine dsvd

   subroutine zsvd(A, U, S, V)
     !> Matrix to be factorized
     complex(kind=wp), intent(in)  :: A(:, :)
     !> Left singular vectors.
     complex(kind=wp), intent(out) :: U(:, :)
     !> Singular values.
     real(kind=wp), intent(out) :: S(:)
     !> Right singular vectors.
     complex(kind=wp), intent(out) :: V(:, :)

     !> Lapack-related.
     character :: jobu = "S", jobvt = "S"
     integer   :: m, n, lda, ldu, ldvt, lwork, info
     complex(kind=wp), allocatable :: work(:)
     real(kind=wp), allocatable :: rwork(:)
     complex(kind=wp) :: A_tilde(size(A, 1), size(A, 2)), vt(min(size(A, 1), size(A, 2)), size(A, 2))

     !> Setup variables.
     m = size(A, 1); n = size(A, 2)
     lda = m; ldu = m; ldvt = n
     lwork = max(1, 3*min(m, n), 5*min(m, n)); allocate (work(lwork)) ; allocate(rwork(5*min(m, n)))

     !> Shape assertion.
     call assert_shape(U, [m, m], "svd", "U")
     call assert_shape(V, [n, n], "svd", "V")

     !> SVD computation.
     a_tilde = a
     call zgesvd(jobu, jobvt, m, n, a_tilde, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info)
     if (info /= 0) then
        write (*, *) "ZGESVD returned info = ", info
        if (info < 0) then
           write (*, *) "The ", -info, "-th argument has an illegal value."
        else
           write (*, *) "ZBSQR did not converge. There are ", info, "superdiagonals"
           write (*, *) "of an intermediate bidiagonal matrix form B which did not"
           write (*, *) "converge to zero. See Lapack documentation for more details."
        end if
        call stop_error("svd: zgesvd error")
     end if

    !> Return the transpose of V.
     v = transpose(vt) ; v = conjg(v)

     return
   end subroutine zsvd

   !-------------------------------------------
   !-----                                 -----
   !-----     LAPACK EVD COMPUTATIONS     -----
   !-----                                 -----
   !-------------------------------------------

   subroutine deig(A, vecs, vals)
      !> Matrix to be factorized.
      real(kind=wp), intent(in) :: A(:, :)
      !> Eigenvectors.
      real(kind=wp), intent(out) :: vecs(:, :)
      !> Eigenvalues.
      complex(kind=wp), intent(out) :: vals(:)

      !> Lapack-related.
      character :: jobvl = "n", jobvr = "v"
      integer   :: n, lwork, info, lda, ldvl, ldvr
      real(kind=wp) :: A_tilde(size(A, 1), size(A, 2)), vr(size(A, 1), size(A, 2))
      real(kind=wp) :: vl(1, size(A, 1))
      real(kind=wp) :: work(4*size(A, 1))
      real(kind=wp) :: wr(size(A, 1)), wi(size(A, 1))
      integer :: i, idx(size(A, 1))

      !> Setup variables.
      n = size(A, 1); lda = n; ldvl = 1; ldvr = n; lwork = 4*n; a_tilde = a

      !> Shape assertion.
      call assert_shape(A, [n, n], "eig", "A")
      call assert_shape(vecs, [n, n], "eig", "vecs")

      !> Eigendecomposition.
      call dgeev(jobvl, jobvr, n, a_tilde, lda, wr, wi, vl, ldvl, vecs, ldvr, work, lwork, info)

      if (info /= 0) then
         write (*, *) "DGEEV returned info = ", info
         if (info < 0) then
            write (*, *) "The ", -info, "-th argument has an illegal value."
         else
            write (*, *) "The QR alg. failed to compute all of the eigenvalues."
            write (*, *) "No eigenvector has been computed."
         end if
         call stop_error("eig: dgeev error")
      end if

      ! !> Real to complex arithmetic.
      ! !> NOTE : Check if a LAPACK function already exists for that purpose.
      vals = cmplx(1.0_wp, 0.0_wp, kind=wp)*wr + cmplx(0.0_wp, 1.0_wp, kind=wp)*wi
      ! vecs = cmplx(0.0_wp, 0.0_wp, kind=wp)*vr

      ! do i = 1, n
      !    if (wi(i) > 0.0_wp) then
      !       vecs(:, i) = cmplx(1.0_wp, 0.0_wp, kind=wp)*vr(:, i) + cmplx(0.0_wp, 1.0_wp, kind=wp)*vr(:, i + 1)
      !    else if (wi(i) < 0.0_wp) then
      !       vecs(:, i) = cmplx(1.0_wp, 0.0_wp, kind=wp)*vr(:, i - 1) - cmplx(0.0_wp, 1.0_wp, kind=wp)*vr(:, i)
      !    else
      !       vecs(:, i) = vr(:, i)
      !    end if
      ! end do

      return
   end subroutine deig

   subroutine deigh(A, vecs, vals)
      !> Matrix to be factorized.
      real(kind=wp), intent(in)  :: A(:, :)
      !> Eigenvectors.
      real(kind=wp), intent(out) :: vecs(:, :)
      !> Eigenvalues.
      real(kind=wp), intent(out) :: vals(:)

      !> Lapack-related.
      character :: jobz = "v", uplo = "u"
      integer   :: n, lwork, info, lda
      real(kind=wp) :: a_tilde(size(A, 1), size(A, 2))
      real(kind=wp) :: work(3*size(A, 1) - 1)

      !> Setup variables.
      n = size(A, 1); lda = n; lwork = 3*n - 1; a_tilde = a

      !> Shape assertion.
      call assert_shape(A, [n, n], "eigh", "A")
      call assert_shape(vecs, [n, n], "eigh", "vecs")

      !> Eigendecomposition.
      call dsyev(jobz, uplo, n, a_tilde, lda, vals, work, lwork, info)

      if (info /= 0) then
         write (*, *) "DSYEV returned info = ", info
         if (info < 0) then
            write (*, *) "The ", -info, "-th argument has an illegal value."
         else
            write (*, *) "The computation failed. See lapack documentation for more details."
         end if
         call stop_error("eigh: dsyev error")
      end if

      return
   end subroutine deigh

   !-----------------------------------------------------
   !-----                                           -----
   !-----     LAPACK LEAST-SQUARES COMPUTATIONS     -----
   !-----                                           -----
   !-----------------------------------------------------

   subroutine lstsq(A, b, x)
      !> Input matrix.
      real(kind=wp), dimension(:, :), intent(in)  :: A
      real(kind=wp), dimension(:), intent(in)  :: b
      real(kind=wp), dimension(:), intent(out) :: x

      !> Lapack job.
      character :: trans = "N"
      integer   :: m, n, nrhs, lda, ldb, lwork, info
      real(kind=wp), dimension(size(A, 1), size(A, 2)) :: A_tilde
      real(kind=wp), dimension(size(A, 1))             :: b_tilde
      real(kind=wp), dimension(:), allocatable         :: work

      !> Initialize variables.
      m = size(A, 1); n = size(A, 2); nrhs = 1
      lda = m; ldb = m; lwork = max(1, min(m, n) + max(min(m, n), nrhs))
      A_tilde = A; b_tilde = b
      allocate (work(1:lwork)); work = 0.0_wp

      !> Solve the least-squares problem.
      call dgels(trans, m, n, nrhs, A_tilde, lda, b_tilde, ldb, work, lwork, info)

      if (info /= 0) then
         write (*, *) "The ", -info, "-th argument has an illegal value."
         call stop_error("lstsq: dgels error")
      end if

      !> Return solution.
      x = b_tilde(1:n)

      return
   end subroutine lstsq

   !----------------------------------------------
   !-----                                    -----
   !-----     LAPACK SCHUR DECOMPOSITION     -----
   !-----                                    -----
   !----------------------------------------------

   subroutine schur(A, Z, eigvals)
     !> Matrix to be factorize.
     real(kind=wp)   , intent(inout) :: A(:, :)
     !> Schur basis.
     real(kind=wp)   , intent(out)   :: Z(:, :)
     !> Eigenvalues.
     complex(kind=wp), intent(out)   :: eigvals(:)

     !> LAPACK-related.
     character :: jobvs="v", sort="n"
     integer   :: n, lda, sdim, ldvs, lwork, info
     real(kind=wp) :: wr(size(A, 1)), wi(size(A, 1))
     real(kind=wp) :: work(3*size(A, 1))
     logical       :: bwork(size(A, 1))

     !> Setup lapack variables.
     n = size(A, 1); lda = max(1, n); ldvs = max(1, n); lwork = max(1, 3*n)

     !> Perform Schur decomposition.
     call dgees(jobvs, sort, dummy_select, n, A, lda, sdim, wr, wi, Z, ldvs, work, lwork, bwork, info)

     return
   end subroutine schur

   subroutine ordschur(T, Q, selected)
     !> Schur matrix to be reordered.
     real(kind=wp), intent(inout) :: T(:, :)
     !> Schur basis to be reordered.
     real(kind=wp), intent(inout) :: Q(:, :)
     !> Array of selected eigenvalues.
     logical      , intent(in)    :: selected(:)

     !> LAPACK-related.
     character :: job="n", compq="v"
     integer   :: info, ldq, ldt, liwork, lwork, m, n, iwork(size(T, 1))
     real(kind=wp) :: s, sep
     real(kind=wp) :: work(size(T, 1)), wr(size(T, 1)), wi(size(T, 1))

     !> Setup variables.
     n = size(T, 2) ; ldt = n ; ldq = n ; lwork = max(1, n) ; liwork = 1

     !> Re-order Schur.
     call dtrsen(job, compq, selected, n, T, ldt, Q, ldq, wr, wi, m, s, sep, work, lwork, iwork, liwork, info)

     return
   end subroutine ordschur

   logical function dummy_select(wr, wi) result(out)
     real(kind=wp), intent(in) :: wr, wi
     return
   end function dummy_select

end module lightkrylov_utils
