program demo
   use stdlib_io_npy, only: save_npy
   use LightKrylov
   use LightKrylov_Constants
   use LightKrylov_Logger
   use LightKrylov_Timing, only: timer => global_lightkrylov_timer
   use Ginzburg_Landau
   implicit none

   character(len=128), parameter :: this_module = 'Example Ginzburg_Landau'

   !------------------------------------------------
   !-----     LINEAR OPERATOR INVESTIGATED     -----
   !------------------------------------------------

   !> Exponential propagator.
   type(exponential_prop), allocatable :: A
   !> Sampling time.
   real(kind=dp), parameter :: tau = 0.1_dp

   !---------------------------------------------------
   !-----     KRYLOV-BASED EIGENDECOMPOSITION     -----
   !---------------------------------------------------

   !> Number of eigenvalues we wish to converge.
   integer, parameter :: nev = 32
   !> Krylov subspace.
   type(state_vector), allocatable :: X(:)
   !> Eigenvalues.
   complex(kind=dp), allocatable :: lambda(:)
   !> Residual.
   real(kind=dp), allocatable    :: residuals(:)
   !> Information flag.
   integer          :: info

   !> Miscellaneous.
   integer       :: i
   complex(dp) :: eigenvectors(nx, nev)

   !=============================================================================

   !----------------------------------
   !-----     INITIALIZATION     -----
   !----------------------------------

   !> Set up logging
   call logger_setup()

   !> Set up timing
   call timer%initialize()
   call timer%add_timer('Ginzburg-Landau example', start=.true.)

   !> Initialize physical parameters.
   call initialize_parameters()

   !> Initialize exponential propagator.
   A = exponential_prop(tau)

   !> Initialize Krylov subspace.
   allocate (X(nev)); call zero_basis(X)

   !------------------------------------------
   !-----     EIGENVALUE COMPUTATION     -----
   !------------------------------------------

   !> Call to LightKrylov.
   call eigs(A, X, lambda, residuals, info)
   call check_info(info, 'eigs', module=this_module, procedure='main')

   !> Transform eigenspectrum from unit-disk to standard complex plane.
   lambda = log(lambda)/tau

   !--------------------------------
   !-----     SAVE TO DISK     -----
   !--------------------------------

   !> Save the eigenspectrum.
   call save_eigenspectrum(lambda, residuals, "example/ginzburg_landau/eigenspectrum.npy")

   !> Reconstruct the leading eigenvectors from the Krylov basis.
   do i = 1, nev
      eigenvectors(:, i) = X(i)%state
   end do

   !> Save eigenvectors to disk.
   call save_npy("example/ginzburg_landau/eigenvectors.npy", eigenvectors)

   ! Print timing info for exponential propagator
   call A%finalize_timer()
   ! Finalize timing
   call timer%finalize()

end program demo
