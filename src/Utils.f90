module Utils
  implicit none
  include "dtypes.h"

  private
  public :: rinv, cinv, svd, evd, hevd, lstsq

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
     integer :: kdim    = 30
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

contains

  !-------------------------------------------
  !-----                                 -----
  !-----     LAPACK MATRIX INVERSION     -----
  !-----                                 -----
  !-------------------------------------------

  subroutine rinv(A)
    !> Matrix to invert (in-place)
    real(kind=wp), intent(inout) :: A(:, :)
    !> Lapack-related.
    integer :: n, info
    real(kind=wp) :: work(size(A, 1))
    integer       :: ipiv(size(A, 1))

    !> Compute A = LU (in-place)
    n = size(A, 1) ; call dgetrf(n, n, A, n, ipiv, info)
    !> Compute inv(A).
    call dgetri(n, A, n, ipiv, work, n, info)

    return
  end subroutine rinv

  subroutine cinv(A)
    !> Matrix to invert (in-place)
    complex(kind=wp), intent(inout) :: A(:, :)
    !> Lapack-related.
    integer :: n, info
    complex(kind=wp) :: work(size(A, 1))
    integer          :: ipiv(size(A, 1))

    !> Compute A = LU (in-place).
    n = size(A, 1) ; call zgetrf(n, n, A, n, ipiv, info)
    !> Compute inv(A).
    call zgetri(n, A, n, ipiv, work, n, info)

    return
  end subroutine cinv

  !------------------------------------------
  !-----                                -----
  !-----     LAPACK SVD COMPUTATION     -----
  !-----                                -----
  !------------------------------------------

  subroutine svd(A, U, S, V)
    !> Matrix to be factorized
    real(kind=wp), intent(in)  :: A(:, :)
    !> Left singular vectors.
    real(kind=wp), intent(out) :: U(size(A, 1), min(size(A, 1), size(A, 2)))
    !> Singular values.
    real(kind=wp), intent(out) :: S(size(A, 2))
    !> Right singular vectors.
    real(kind=wp), intent(out) :: V(size(A, 2), min(size(A, 1), size(A, 2)))

    !> Lapack-related.
    character :: jobu="S", jobvt="S"
    integer   :: m, n, lda, ldu, ldvt, lwork, info
    real(kind=wp), allocatable :: work(:)
    real(kind=wp) :: A_tilde(size(A, 1), size(A, 2)), vt(min(size(A, 1), size(A, 2)), size(A, 2))

    !> Setup variables.
    m = size(A, 1) ; n = size(A, 2)
    lda = m ; ldu = m ; ldvt = n
    lwork = max(1, 3*min(m, n), 5*min(m, n)) ; allocate(work(lwork))

    !> SVD computation.
    a_tilde = a
    call dgesvd(jobu, jobvt, m, n, a_tilde, lda, s, u, ldu, vt, ldvt, work, lwork, info)
    v = transpose(vt)

    return
  end subroutine svd

  !-------------------------------------------
  !-----                                 -----
  !-----     LAPACK EVD COMPUTATIONS     -----
  !-----                                 -----
  !-------------------------------------------

  subroutine evd(A, vecs, vals)
    !> Matrix to be factorized.
    real(kind=wp), intent(in) :: A(:, :)
    !> Eigenvectors.
    complex(kind=wp), intent(out) :: vecs(size(A, 1), size(A, 2))
    !> Eigenvalues.
    complex(kind=wp), intent(out) :: vals(size(A, 1))

    !> Lapack-related.
    character :: jobvl="n", jobvr="v"
    integer   :: n, lwork, info, lda, ldvl, ldvr
    real(kind=wp) :: A_tilde(size(A, 1), size(A, 2)), vr(size(A, 1), size(A, 2))
    real(kind=wp) :: vl(1, size(A, 1))
    real(kind=wp) :: work(4*size(A, 1))
    real(kind=wp) :: wr(size(A, 1)), wi(size(A, 1))
    integer :: i, idx(size(A, 1))

    !> Setup variables.
    n = size(A, 1) ; lda = n ; ldvl = 1 ; ldvr = n ; lwork = 4*n ; a_tilde = a
    !> Eigendecomposition.
    call dgeev(jobvl, jobvr, n, a_tilde, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
    !> Real to complex arithmetic.
    !> NOTE : Check if a LAPACK function already exists for that purpose.
    vals = cmplx(1.0_wp, 0.0_wp, kind=wp)*wr + cmplx(0.0_wp, 1.0_wp, kind=wp)*wi
    vecs = cmplx(1.0_wp, 0.0_wp, kind=wp)*vr

    do i = 1, n-1
       if (wi(i) .gt. 0) then
          vecs(:, i)   = cmplx(1.0_wp, 0.0_wp, kind=wp)*vr(:, i) + cmplx(0.0_wp, 1.0_wp, kind=wp)*vr(:, i+1)
          vecs(:, i+1) = conjg(vecs(:, i))
       else if (abs(wi(i)) .le. epsilon(wi(i))) then
          vecs(:, i) = cmplx(1.0_wp, 0.0_wp, kind=wp)*vr(:, i)
       end if
    enddo

    return
  end subroutine evd

  subroutine hevd(A, vecs, vals)
    !> Matrix to be factorized.
    real(kind=wp), intent(in)  :: A(:, :)
    !> Eigenvectors.
    real(kind=wp), intent(out) :: vecs(size(A, 1), size(A, 2))
    !> Eigenvalues.
    real(kind=wp), intent(out) :: vals(size(A, 1))

    !> Lapack-related.
    character :: jobz="v", uplo="u"
    integer   :: n, lwork, info, lda
    real(kind=wp) :: a_tilde(size(A, 1), size(A, 2))
    real(kind=wp) :: work(3*size(A, 1) - 1)

    !> Setup variables.
    n = size(A, 1) ; lda = n ; lwork = 3*n-1 ; a_tilde = a
    !> Eigendecomposition.
    call dsyev(jobz, uplo, n, a_tilde, lda, vals, work, lwork, info)

    return
  end subroutine hevd

  !-----------------------------------------------------
  !-----                                           -----
  !-----     LAPACK LEAST-SQUARES COMPUTATIONS     -----
  !-----                                           -----
  !-----------------------------------------------------

  subroutine lstsq(A, b, x)
    !> Input matrix.
    real(kind=wp), dimension(:, :), intent(in)  :: A
    real(kind=wp), dimension(:)   , intent(in)  :: b
    real(kind=wp), dimension(:)   , intent(out) :: x

    !> Lapack job.
    character :: trans = "N"
    integer   :: m, n, nrhs, lda, ldb, lwork, info
    real(kind=wp), dimension(size(A, 1), size(A, 2)) :: A_tilde
    real(kind=wp), dimension(size(A, 1))             :: b_tilde
    real(kind=wp), dimension(:), allocatable         :: work

    !> Initialize variables.
    m = size(A, 1) ; n = size(A, 2) ; nrhs = 1
    lda = m ; ldb = m ; lwork = max(1, min(m, n) + max(min(m, n), nrhs))
    A_tilde = A ; b_tilde = b
    allocate(work(1:lwork)) ; work = 0.0_wp

    !> Solve the least-squares problem.
    call dgels(trans, m, n, nrhs, A_tilde, lda, b_tilde, ldb, work, lwork, info)

    if (info > 0) then
       write(*, *) "INFO: Error in LAPACK least-squares solver. Stopping the computation."
       call exit()
    endif

    !> Return solution.
    x = b_tilde(1:n)

    return
  end subroutine lstsq



end module Utils
