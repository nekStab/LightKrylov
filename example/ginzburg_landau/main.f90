program demo
  use LightKrylov
  use LightKrylov, only: wp => dp
  use Ginzburg_Landau
  use stdlib_io_npy, only : save_npy
  implicit none

  !------------------------------------------------
  !-----     LINEAR OPERATOR INVESTIGATED     -----
  !------------------------------------------------

  !> Exponential propagator.
  type(exponential_prop), allocatable :: A
  !> Sampling time.
  real(kind=wp), parameter :: tau = 0.1_wp

  !---------------------------------------------------
  !-----     KRYLOV-BASED EIGENDECOMPOSITION     -----
  !---------------------------------------------------

  !> Number of eigenvalues we wish to converge.
  integer, parameter :: nev = 32
  !> Krylov subspace.
  type(state_vector), allocatable :: X(:)
  !> Eigenvalues.
  complex(kind=wp), allocatable :: lambda(:)
  !> Residual.
  real(kind=wp), allocatable    :: residuals(:)
  !> Information flag.
  integer          :: info

  !> Miscellaneous.
  integer       :: i, j, k
  real(kind=wp) :: alpha
  complex(wp) :: eigenvectors(nx, nev)

  !=============================================================================

  !----------------------------------
  !-----     INITIALIZATION     -----
  !----------------------------------

  !> Initialize physical parameters.
  call initialize_parameters()

  !> Initialize exponential propagator.
  A = exponential_prop(tau)

  !> Initialize Krylov subspace.
  allocate(X(nev)) ; call initialize_krylov_subspace(X)

  !------------------------------------------
  !-----     EIGENVALUE COMPUTATION     -----
  !------------------------------------------

  !> Call to LightKrylov.
  call eigs(A, X, lambda, residuals, info)

  !> Transform eigenspectrum from unit-disk to standard complex plane.
  lambda = log(lambda) / tau

  !--------------------------------
  !-----     SAVE TO DISK     -----
  !--------------------------------

  !> Save the eigenspectrum.
  call save_eigenspectrum(lambda, residuals, "example/ginzburg_landau/eigenspectrum.npy")

  !> Reconstruct the leading eigenvectors from the Krylov basis.
  do i = 1, nev
    eigenvectors(:, i) = X(i)%state
  enddo

  !> Save eigenvectors to disk.
  call save_npy("example/ginzburg_landau/eigenvectors.npy", eigenvectors)

end program demo
