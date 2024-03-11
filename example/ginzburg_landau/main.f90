program demo
  use LightKrylov
  use Ginzburg_Landau
  use stdlib_io_npy, only : save_npy
  implicit none

  !------------------------------------------------
  !-----     LINEAR OPERATOR INVESTIGATED     -----
  !------------------------------------------------

  !> Exponential propagator.
  type(exponential_prop), allocatable :: A
  !> Sampling time.
  real(kind=wp), parameter :: tau = 1.0_wp

  !---------------------------------------------------
  !-----     KRYLOV-BASED EIGENDECOMPOSITION     -----
  !---------------------------------------------------

  !> Number of eigenvalues we wish to converge.
  integer, parameter :: nev = 128
  !> Krylov subspace dimension.
  integer, parameter :: kdim = 2*nx
  !> Krylov subspace.
  type(state_vector), allocatable :: X(:)
  !> Eigenvalues.
  complex(kind=wp) :: lambda(kdim)
  !> Residual.
  real(kind=wp)    :: residuals(kdim)
  !> Information flag.
  integer          :: info

  !> Miscellaneous.
  integer       :: i, j, k
  real(kind=wp) :: alpha
  class(abstract_vector), allocatable :: wrk
  complex(kind=wp)                    :: eigenvectors(nx, nev)

  !=============================================================================

  !----------------------------------
  !-----     INITIALIZATION     -----
  !----------------------------------

  !> Initialize physical parameters.
  call initialize_parameters()

  !> Initialize exponential propagator.
  A = exponential_prop(tau)

  !> Initialize Krylov subspace.
  allocate(X(1:kdim+1)) ; call initialize_krylov_subspace(X)

  !> Random initial Krylov vector.
  call random_number(X(1)%state) ; alpha = X(1)%norm() ; call X(1)%scal(1.0_wp / alpha)

  !------------------------------------------
  !-----     EIGENVALUE COMPUTATION     -----
  !------------------------------------------

  !> Call to LightKrylov.
  call eigs(A, X, lambda, residuals, info, nev=nev)

  !> Transform eigenspectrum from unit-disk to standard complex plane.
  lambda = log(lambda) / tau

  !--------------------------------
  !-----     SAVE TO DISK     -----
  !--------------------------------

  !> Save the eigenspectrum.
  call save_eigenspectrum(lambda%re, lambda%im, residuals, "example/ginzburg_landau/eigenspectrum.npy")

  !> Reconstruct the leading eigenvectors from the Krylov basis.
  do i = 1, nev
     !> Real part.
     eigenvectors(:, i)%re = X(i)%state(1:nx)
     !> Imaginary part.
     eigenvectors(:, i)%im = X(i+1)%state(1:nx)
  enddo

  !> Save eigenvectors to disk.
  call save_npy("example/ginzburg_landau/eigenvectors.npy", eigenvectors)

end program demo
