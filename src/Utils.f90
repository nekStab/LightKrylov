module lightkrylov_utils
  !! This module provides a set of utilities used throughout `LightKrylov`.
  !! It also provides a selection wrapper around LAPACK to perform standard linear algebra computations.
   use iso_fortran_env, only: output_unit
   implicit none
   include "dtypes.h"

   private
   ! General-purpose utilities.
   public :: assert_shape, stop_error, iargsort, norml, log2
   ! Print matrices
   public :: print_mat
   ! Linear Algebra Utilities.
   public :: inv, svd, eig, eigh, lstsq, schur, ordschur

   !-------------------------------------------------------
   !-----                                             -----
   !-----     OPTS TYPE OBJECT FOR LINEAR SOLVERS     -----
   !-----                                             -----
   !-------------------------------------------------------

   type, abstract, public :: abstract_opts
      !! Abstract type container for options to be passed to the various iterative solvers.
   end type abstract_opts

   type, extends(abstract_opts), public :: gmres_opts
      !! Extended `abstract_opts` type to pass options to the `gmres` solver.
      integer :: kdim = 30
      !! Dimension of the Krylov subspace (default: 30).
      integer :: maxiter = 10
      !! Maximum number of `gmres` restarts (default: 10)
      real(kind=wp) :: atol = atol
      !! Absolute tolerance (default: `epsilon(1.0_wp)`).
      real(kind=wp) :: rtol = rtol
      !! Relative tolerance (default: `sqrt(atol)`).
      logical :: verbose = .false.
      !! Verbosity control (default: `.false.`).
   end type gmres_opts

   type, extends(abstract_opts), public :: cg_opts
      !! Extended `abstract_opts` type to pass options to the `cg` solver.
      integer :: maxiter = 100
      !! Maximum number of `cg` iterations (default: 100).
      real(kind=wp) :: rtol = rtol
      !! Relative tolerance (default: `sqrt(atol)`).
      real(kind=wp) :: atol = atol
      !! Absolute tolerance (default: `epsilon(1.0_wp)`).
      logical :: verbose = .false.
      !! Verbosity control (default: `.false.`).
   end type cg_opts

   !------------------------------
   !-----     INTERFACES     -----
   !------------------------------

   interface assert_shape
      !! Interface to assert the shape of a matrix.
      module procedure dassert_shape
      module procedure zassert_shape
   end interface assert_shape

   interface inv
      !! Interface to compute the inverse of a matrix (in-place).
      module procedure dinv
      module procedure zinv
   end interface inv

   interface svd
      !! Interface to compute the SVD of a matrix.
      module procedure dsvd
      module procedure zsvd
   end interface svd

   interface eig
      !! Interface to compute the EVD of a matrix.
      module procedure deig
   end interface eig

   interface eigh
      !! Interface to compute the EVD of a sym. pos. def. matrix.
      module procedure deigh
   end interface eigh

contains

   !-------------------------------------
   !-----                           -----
   !-----     VARIOUS UTILITIES     -----
   !-----                           -----
   !-------------------------------------

   subroutine stop_error(msg)
     !! Utility function to print an error message.
      character(len=*), intent(in) :: msg
     !! Error message.
      write (output_unit, *) msg; stop 1
      return
   end subroutine stop_error

   subroutine dassert_shape(A, size, routine, matname)
     !! Utility function to assert the shape of a real-valued matrix.
      real(kind=wp), intent(in) :: A(:, :)
     !! Matrix whose dimensions need to be asserted.
      integer, intent(in) :: size(:)
     !! Expected dimensions of A.
      character(len=*), intent(in) :: routine
     !! Name of the routine where assertion is done.
      character(len=*), intent(in) :: matname
     !! Name of the asserted matrix.

      if (any(shape(A) /= size)) then
         write (output_unit, *) "In routine "//routine//" matrix "//matname//" has illegal shape ", shape(A)
         write (output_unit, *) "Expected shape is ", size
         call stop_error("Aborting due to illegal matrix operation.")
      end if

      return
   end subroutine dassert_shape

   subroutine zassert_shape(A, size, routine, matname)
     !! Utility function to assert the shape of a complex-valued matrix.
      complex(kind=wp), intent(in) :: A(:, :)
     !! Matrix whose dimensions need to be asserted.
      integer, intent(in) :: size(:)
     !! Expected dimensions of A.
      character(len=*), intent(in) :: routine
     !! Name of the routine where assertion is done.
      character(len=*), intent(in) :: matname
     !! Name of the asserted matrix.

      if (any(shape(A) /= size)) then
         write (output_unit, *) "In routine "//routine//" matrix "//matname//" has illegal shape ", shape(A)
         write (output_unit, *) "Expected shape is ", size
         call stop_error("Aborting due to illegal matrix operation.")
      end if

      return
   end subroutine zassert_shape

   subroutine print_mat(m, n, A, name)
      !! Utility function to print a m-by-n matrix with a title
      integer , intent(in)                         :: n
      integer , intent(in)                         :: m
      real (kind=wp)  , intent(in)                 :: a(m,n)
      character(*) , optional,  intent(in)         :: name

      ! internal variables
      real (kind=wp) :: amin, amax
      character ( len = 10 ) iform
      integer :: i, ihi, ilo, j, jhi, jlo, lmax, npline
      logical :: integ
      
      write(*,*)
      if (present(name)) then
         write(*,*) 'Output matrix: ', trim(name)
      endif
      
      ! Check if all entries are integral. 
      integ = .true.

      do i = 1, m
         do j = 1, n
            if ( integ ) then
               if ( a(i,j) /= real ( int ( a(i,j) ),kind=wp) ) then
                  integ = .false.
               end if
            end if
         end do
      end do

      ! Find the maximum and minimum entries.
      amax = maxval ( a(1:m,1:n) )
      amin = minval ( a(1:m,1:n) )

      ! Use the information about the maximum size of an entry to
      ! compute an intelligent format for use with integer entries.      
      if ( amax .lt. 1e-12) amax = 1e-12 
      ! to avoid problems with zero matrix
      lmax = int ( log10 ( amax ) )
      if ( lmax .lt. -2 ) lmax = 0
   
      if ( integ ) then
         npline = 79 / ( lmax + 3 )
         write ( iform, '(''('',i2,''I'',i2,'')'')' ) npline, lmax+3
      else
         npline = 8
         iform = ' '
      end if

      ! Print a scalar quantity.
      if ( m == 1 .and. n == 1 ) then
         if ( integ ) then
            write ( *, iform ) int ( a(1,1) )
         else
            write ( *, '(2x,g10.2)' ) a(1,1)
         end if

      ! Column vector of length M,
      else if ( n == 1 ) then
         do ilo = 1, m, npline
            ihi = min ( ilo+npline-1, m )
            if ( integ ) then
               write ( *, iform ) ( int ( a(i,1) ), i = ilo, ihi )
            else
               write ( *, '(2x,8g10.2)' ) a(ilo:ihi,1)
            end if
         end do

      ! Row vector of length N,
      else if ( m == 1 ) then
         do jlo = 1, n, npline
            jhi = min ( jlo+npline-1, n )
            if ( integ ) then
               write ( *, iform ) int ( a(1,jlo:jhi) )
            else
               write ( *, '(2x,8g10.2)' ) a(1,jlo:jhi)
            end if
         end do
      
      ! M by N Array
      else
         do jlo = 1, n, npline
            jhi = min ( jlo+npline-1, n )
            if ( npline < n ) then
               write ( *, '(a)' ) ' '
               write ( *, '(a,i8,a,i8)' ) 'Matrix columns ', jlo, ' to ', jhi
               write ( *, '(a)' ) ' '
            end if
            do i = 1, m
               if ( integ ) then
                  write ( *, iform ) int ( a(i,jlo:jhi) )
               else
                  write ( *, '(2x,8g14.6)' ) a(i,jlo:jhi)
               end if
            end do
         end do
      end if
      write(*,*)
   
      return
   end subroutine print_mat

   !-------------------------------------
   !-----                           -----
   !-----     VARIOUS FUNCTIONS     -----
   !-----                           -----
   !-------------------------------------

   function log2(x) result(y)
      !! compute the base-2 logarithm of the input
      implicit none
      real(kind=wp), intent(in) :: x
      real(kind=wp) :: y
      y = log(x) / log(2.0_wp)
   end function log2
 
   function norml(A) result(norm)
      !! compute the infinity norm of the real-valued input matrix A
      implicit none   
      real(kind=wp), intent(in) :: A(:,:)
      ! Internal variables
      integer       :: i, n
      real(kind=wp) :: row_sum, norm
      
      norm = 0.0_wp
      n = size(A,1)
      do i = 1, n
      row_sum = sum ( abs ( A(i,1:n) ) )
      norm = max ( norm, row_sum )
      end do
   end function norml

   function iargsort(a) result(idx)
      integer, intent(in):: a(:)    ! array of numbers
      integer :: idx(size(a))       ! indices into the array 'a' that sort it
      integer :: N                  ! number of numbers/vectors
      integer :: i,imin             ! indices: i, i of smallest
      integer :: temp               ! temporary
      integer :: a2(size(a))
      a2 = a
      N=size(a)
      do i = 1, N
          idx(i) = i
      end do
      do i = 1, N-1
          ! find ith smallest in 'a'
          imin = minloc(a2(i:),1) + i - 1
          ! swap to position i in 'a' and 'b', if not already there
          if (imin /= i) then
              temp = a2(i);a2(i) = a2(imin); a2(imin) = temp
              temp = idx(i); idx(i) = idx(imin); idx(imin) = temp
          end if
      end do
   end function

   !-------------------------------------------
   !-----                                 -----
   !-----     LAPACK MATRIX INVERSION     -----
   !-----                                 -----
   !-------------------------------------------

   subroutine dinv(A)
     !! In-place inversion of a real-valued matrix A using LAPACK.
      real(kind=wp), intent(inout) :: A(:, :)
      !! Matrix to be inverted (in-place).

      ! Internal variables.
      integer :: n, info
      real(kind=wp) :: work(size(A, 1))
      integer       :: ipiv(size(A, 1))

      ! Compute A = LU (in-place)
      n = size(A, 1); call assert_shape(A, [n, n], "inv", "A")
      call dgetrf(n, n, A, n, ipiv, info)
      if (info /= 0) then
         write (output_unit, *) "DGETRF returned info = ", info
         if (info < 0) then
            write (output_unit, *) "The ", -info, "-th argument has an illegal value."
         else
            write (output_unit, *) "U(", info, ",", info, ") is exactly zero. The factorization"
            write (output_unit, *) "has been completed but the factor U is exactly singular."
            write (output_unit, *) "Division by zero will occur if used to solve Ax = b."
         end if
         call stop_error("inv: DGETRF error.")
      end if

      ! Compute inv(A).
      call dgetri(n, A, n, ipiv, work, n, info)
      if (info /= 0) then
         write (output_unit, *) "DGETRI returned info =", info
         if (info < 0) then
            write (output_unit, *) "The ", -info, "-th argument has an illegal value."
         else
            write (output_unit, *) "U(", info, ",", info, ") is exactly zero."
            write (output_unit, *) "The matrix is singular and its inverse cannot be computed."
         end if
         call stop_error("inv: DGETRI error.")
      end if

      return
   end subroutine dinv

   subroutine zinv(A)
     !! In-place inversion of a complex-valued matrix using LAPACK.
      complex(kind=wp), intent(inout) :: A(:, :)
      !! Matrix be inverted (in-place).

      ! Internal variables.
      integer :: n, info
      complex(kind=wp) :: work(size(A, 1))
      integer          :: ipiv(size(A, 1))

      ! Compute A = LU (in-place).
      n = size(A, 1); call assert_shape(A, [n, n], "inv", "A")
      call zgetrf(n, n, A, n, ipiv, info)
      if (info /= 0) then
         write (output_unit, *) "ZGETRF returned info = ", info
         if (info < 0) then
            write (output_unit, *) "The ", -info, "-th argument has an illegal value."
         else
            write (output_unit, *) "U(", info, ",", info, ") is exactly zero. The factorization"
            write (output_unit, *) "has been completed but the factor U is exactly singular."
            write (output_unit, *) "Division by zero will occur if used to solve Ax = b."
         end if
         call stop_error("inv: ZGETRF error.")
      end if

      ! Compute inv(A).
      call zgetri(n, A, n, ipiv, work, n, info)
      if (info /= 0) then
         write (output_unit, *) "ZGETRI returned info =", info
         if (info < 0) then
            write (output_unit, *) "The ", -info, "-th argument has an illegal value."
         else
            write (output_unit, *) "U(", info, ",", info, ") is exactly zero."
            write (output_unit, *) "The matrix is singular and its inverse cannot be computed."
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
     !! Singular Value Decomposition of a real-valued matrix using LAPACK.
      real(kind=wp), intent(in)  :: A(:, :)
      !! Matrix to be factorized.
      real(kind=wp), intent(out) :: U(:, :)
      !! Left singular vectors.
      real(kind=wp), intent(out) :: S(:)
      !! Singular values.
      real(kind=wp), intent(out) :: V(:, :)
      !! Right singular vectors.

      ! Lapack-related.
      character :: jobu = "S", jobvt = "S"
      integer   :: m, n, lda, ldu, ldvt, lwork, info
      real(kind=wp), allocatable :: work(:)
      real(kind=wp) :: A_tilde(size(A, 1), size(A, 2)), vt(min(size(A, 1), size(A, 2)), size(A, 2))

      ! Setup variables.
      m = size(A, 1); n = size(A, 2)
      lda = m; ldu = m; ldvt = n
      lwork = max(1, 3*min(m, n) + max(m, n), 5*min(m, n)); allocate (work(lwork))

      ! Shape assertions.
      call assert_shape(U, [m, m], "svd", "U")
      call assert_shape(V, [n, n], "svd", "V")

      ! SVD computation.
      a_tilde = a
      call dgesvd(jobu, jobvt, m, n, a_tilde, lda, s, u, ldu, vt, ldvt, work, lwork, info)
      if (info /= 0) then
         write (output_unit, *) "DGESVD returned info = ", info
         if (info < 0) then
            write (output_unit, *) "The ", -info, "-th argument has an illegal value."
         else
            write (output_unit, *) "DBSQR did not converge. There are ", info, "superdiagonals"
            write (output_unit, *) "of an intermediate bidiagonal matrix form B which did not"
            write (output_unit, *) "converge to zero. See Lapack documentation for more details."
         end if
         call stop_error("svd: dgesvd error")
      end if

      ! Return the transpose of V.
      v = transpose(vt)

      return
   end subroutine dsvd

   subroutine zsvd(A, U, S, V)
     !! Singular Value Decomposition of a complex-valued matrix using LAPACK.
      complex(kind=wp), intent(in)  :: A(:, :)
     !! Matrix to be factorized.
      complex(kind=wp), intent(out) :: U(:, :)
     !! Left singular vectors.
      real(kind=wp), intent(out) :: S(:)
     !! Singular values.
      complex(kind=wp), intent(out) :: V(:, :)
     !! Right singular vectors.

      ! Lapack-related.
      character :: jobu = "S", jobvt = "S"
      integer   :: m, n, lda, ldu, ldvt, lwork, info
      complex(kind=wp), allocatable :: work(:)
      real(kind=wp), allocatable :: rwork(:)
      complex(kind=wp) :: A_tilde(size(A, 1), size(A, 2)), vt(min(size(A, 1), size(A, 2)), size(A, 2))

      ! Setup variables.
      m = size(A, 1); n = size(A, 2)
      lda = m; ldu = m; ldvt = n
      lwork = max(1, 3*min(m, n), 5*min(m, n)); allocate (work(lwork)); allocate (rwork(5*min(m, n)))

      ! Shape assertion.
      call assert_shape(U, [m, m], "svd", "U")
      call assert_shape(V, [n, n], "svd", "V")

      ! SVD computation.
      a_tilde = a
      call zgesvd(jobu, jobvt, m, n, a_tilde, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info)
      if (info /= 0) then
         write (output_unit, *) "ZGESVD returned info = ", info
         if (info < 0) then
            write (output_unit, *) "The ", -info, "-th argument has an illegal value."
         else
            write (output_unit, *) "ZBSQR did not converge. There are ", info, "superdiagonals"
            write (output_unit, *) "of an intermediate bidiagonal matrix form B which did not"
            write (output_unit, *) "converge to zero. See Lapack documentation for more details."
         end if
         call stop_error("svd: zgesvd error")
      end if

      ! Return the transpose of V.
      v = transpose(vt); v = conjg(v)

      return
   end subroutine zsvd

   !-------------------------------------------
   !-----                                 -----
   !-----     LAPACK EVD COMPUTATIONS     -----
   !-----                                 -----
   !-------------------------------------------

   subroutine deig(A, vecs, vals)
     !! Eigenvalue decomposition of a real-valued matrix using LAPACK.
      real(kind=wp), intent(in) :: A(:, :)
      !! Matrix to be factorized.
      real(kind=wp), intent(out) :: vecs(:, :)
      !! Eigenvectors.
      complex(kind=wp), intent(out) :: vals(:)
      !! Eigenvalues.

      ! Lapack-related.
      character :: jobvl = "n", jobvr = "v"
      integer   :: n, lwork, info, lda, ldvl, ldvr
      real(kind=wp) :: A_tilde(size(A, 1), size(A, 2)), vr(size(A, 1), size(A, 2))
      real(kind=wp) :: vl(1, size(A, 1))
      real(kind=wp) :: work(4*size(A, 1))
      real(kind=wp) :: wr(size(A, 1)), wi(size(A, 1))
      integer :: i, idx(size(A, 1))

      ! Setup variables.
      n = size(A, 1); lda = n; ldvl = 1; ldvr = n; lwork = 4*n; a_tilde = a

      ! Shape assertion.
      call assert_shape(A, [n, n], "eig", "A")
      call assert_shape(vecs, [n, n], "eig", "vecs")

      ! Eigendecomposition.
      call dgeev(jobvl, jobvr, n, a_tilde, lda, wr, wi, vl, ldvl, vecs, ldvr, work, lwork, info)

      if (info /= 0) then
         write (output_unit, *) "DGEEV returned info = ", info
         if (info < 0) then
            write (output_unit, *) "The ", -info, "-th argument has an illegal value."
         else
            write (output_unit, *) "The QR alg. failed to compute all of the eigenvalues."
            write (output_unit, *) "No eigenvector has been computed."
         end if
         call stop_error("eig: dgeev error")
      end if

      ! Real to complex arithmetic.
      ! NOTE : Check if a LAPACK function already exists for that purpose.
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
     !! Eigenvalue decomposition of a real-valued sym. pos. def. matrix using LAPACK.
      real(kind=wp), intent(in)  :: A(:, :)
      !! Matrix to be factorized.
      real(kind=wp), intent(out) :: vecs(:, :)
      !! Eigenvectors.
      real(kind=wp), intent(out) :: vals(:)
      !! Eigenvalues.

      ! Lapack-related.
      character :: jobz = "v", uplo = "u"
      integer   :: n, lwork, info, lda
      real(kind=wp) :: a_tilde(size(A, 1), size(A, 2))
      real(kind=wp) :: work(3*size(A, 1) - 1)

      ! Setup variables.
      n = size(A, 1); lda = n; lwork = 3*n - 1; a_tilde = a

      ! Shape assertion.
      call assert_shape(A, [n, n], "eigh", "A")
      call assert_shape(vecs, [n, n], "eigh", "vecs")

      ! Eigendecomposition.
      call dsyev(jobz, uplo, n, a_tilde, lda, vals, work, lwork, info)

      if (info /= 0) then
         write (output_unit, *) "DSYEV returned info = ", info
         if (info < 0) then
            write (output_unit, *) "The ", -info, "-th argument has an illegal value."
         else
            write (output_unit, *) "The computation failed. See lapack documentation for more details."
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
     !! Solves a linear least-squares problem \( \min~\| \mathbf{Ax} - \mathbf{b} \|_2^2 \) using LAPACK.
      real(kind=wp), dimension(:, :), intent(in)  :: A
      !! Matrix to be "pseudo-inverted".
      real(kind=wp), dimension(:), intent(in)  :: b
      !! Right-hand side vector.
      real(kind=wp), dimension(:), intent(out) :: x
      !! Solution of the least-squares problem.

      ! Lapack job.
      character :: trans = "N"
      integer   :: m, n, nrhs, lda, ldb, lwork, info
      real(kind=wp), dimension(size(A, 1), size(A, 2)) :: A_tilde
      real(kind=wp), dimension(size(A, 1))             :: b_tilde
      real(kind=wp), dimension(:), allocatable         :: work

      ! Initialize variables.
      m = size(A, 1); n = size(A, 2); nrhs = 1
      lda = m; ldb = m; lwork = max(1, min(m, n) + max(min(m, n), nrhs))
      A_tilde = A; b_tilde = b
      allocate (work(1:lwork)); work = 0.0_wp

      ! Solve the least-squares problem.
      call dgels(trans, m, n, nrhs, A_tilde, lda, b_tilde, ldb, work, lwork, info)

      if (info /= 0) then
         write (output_unit, *) "The ", -info, "-th argument has an illegal value."
         call stop_error("lstsq: dgels error")
      end if

      ! Return solution.
      x = b_tilde(1:n)

      return
   end subroutine lstsq

   !----------------------------------------------
   !-----                                    -----
   !-----     LAPACK SCHUR DECOMPOSITION     -----
   !-----                                    -----
   !----------------------------------------------

   subroutine schur(A, Z, eigvals)
     !! Compute the Schur form (in-place) and Schur vectors of a real-valued matrix.
      real(kind=wp), intent(inout) :: A(:, :)
     !! Matrix to be factorized.
      real(kind=wp), intent(out)   :: Z(:, :)
     !! Schur basis.
      complex(kind=wp), intent(out)   :: eigvals(:)
     !! Eigenvalues.

      ! LAPACK-related.
      character :: jobvs = "v", sort = "n"
      integer   :: n, lda, sdim, ldvs, lwork, info
      real(kind=wp) :: wr(size(A, 1)), wi(size(A, 1))
      real(kind=wp) :: work(3*size(A, 1))
      logical       :: bwork(size(A, 1))

      ! Setup lapack variables.
      n = size(A, 1); lda = max(1, n); ldvs = max(1, n); lwork = max(1, 3*n)

      ! Perform Schur decomposition.
      call dgees(jobvs, sort, dummy_select, n, A, lda, sdim, wr, wi, Z, ldvs, work, lwork, bwork, info)

      ! Eigenvalues.
      eigvals = cmplx(wr, wi, kind=wp)

      return
   end subroutine schur

   subroutine ordschur(T, Q, selected)
     !! Re-order the Schur factorization of a real-valued matrix by moving the selected eigenvalues
     !! in the upper-left block.
      real(kind=wp), intent(inout) :: T(:, :)
     !! Schur matrix to be reordered.
      real(kind=wp), intent(inout) :: Q(:, :)
     !! Schur vectors to be reordered.
      logical, intent(in)    :: selected(:)
     !! Boolean array defining the selected eigenvalues.

      ! LAPACK-related.
      character :: job = "n", compq = "v"
      integer   :: info, ldq, ldt, liwork, lwork, m, n, iwork(size(T, 1))
      real(kind=wp) :: s, sep
      real(kind=wp) :: work(size(T, 1)), wr(size(T, 1)), wi(size(T, 1))

      ! Setup variables.
      n = size(T, 2); ldt = n; ldq = n; lwork = max(1, n); liwork = 1

      ! Re-order Schur.
      call dtrsen(job, compq, selected, n, T, ldt, Q, ldq, wr, wi, m, s, sep, work, lwork, iwork, liwork, info)

      return
   end subroutine ordschur

   pure logical function dummy_select(wr, wi) result(out)
      real(kind=wp), intent(in) :: wr, wi
      return
   end function dummy_select

end module lightkrylov_utils
