module Ginzburg_Landau_Base
   ! Standard Library.
   use stdlib_stats_distribution_normal, only: normal => rvs_normal
   use stdlib_optval, only : optval
   use stdlib_math, only : linspace
   use stdlib_io_npy, only: save_npy
   use stdlib_strings, only: replace_all
   use stdlib_linalg, only: svdvals, eye
   ! LightKrylov for linear algebra.
   use LightKrylov
   use LightKrylov, only: wp => dp
   use LightKrylov_Logger
   use LightKrylov_Utils, only : assert_shape
   use LightKrylov_AbstractVectors
   implicit none

   private :: this_module
   character(len=*), parameter :: this_module = 'Ginzburg_Landau_Base'
   
   public  :: initialize_parameters
   
   public  :: L, nx, N, dx
   public  :: nu, gamma, mu_0, c_mu, mu_2, mu
   public  :: weight
   !public  :: rk_b, x_b, s_b, rk_c, x_c, s_c
   !public  :: B, CT, weight, weight_mat
   !public  :: BBTW, CTCW
   !public  :: Qc, Rinv, CTQcCW, BRinvBTW
   
   !-------------------------------
   !-----     PARAMETERS 1    -----
   !-------------------------------

   ! Mesh related parameters.
   real(wp), parameter :: L  = 50.0_wp ! Domain length
   integer,  parameter :: nx = 256      ! Number of grid points (excluding boundaries).
   real(wp)            :: dx           ! Grid size.

   !-------------------------------------------
   !-----     LIGHTKRYLOV VECTOR TYPE     -----
   !-------------------------------------------

   type, extends(abstract_vector_rdp), public :: state_vector
      real(wp) :: state(2*nx) = 0.0_wp
   contains
      private
      procedure, pass(self), public :: zero
      procedure, pass(self), public :: dot
      procedure, pass(self), public :: scal
      procedure, pass(self), public :: axpby
      procedure, pass(self), public :: rand
      procedure, pass(self), public :: get_size
   end type state_vector

   !-------------------------------
   !-----     PARAMETERS 2    -----
   !-------------------------------

   ! Physical parameters.
   complex(wp), parameter :: nu    = cmplx(2.0_wp, 0.2_wp, wp)
   complex(wp), parameter :: gamma = cmplx(1.0_wp, -1.0_wp, wp)
   real(wp),    parameter :: mu_0  = 0.38_wp
   real(wp),    parameter :: c_mu  = 0.2_wp
   real(wp),    parameter :: mu_2  = -0.01_wp
   real(wp)               :: mu(nx)

   ! Input-Output system parameters
   real(wp)               :: weight(2*nx)       ! integration weights
   !integer,  parameter    :: rk_b = 2           ! number of inputs to the system
   !real(wp), parameter    :: x_b = -11.0_wp     ! location of input Gaussian
   !real(wp), parameter    :: s_b = 1.0_wp       ! variance of input Gaussian
   !type(state_vector)     :: B(rk_b)
   !real(wp), parameter    :: x_c = sqrt(-2.0_wp*(mu_0 - c_mu**2)/mu_2) ! location of input Gaussian
   !real(wp), parameter    :: s_c = 1.0_wp       ! variance of input Gaussian
   !integer,  parameter    :: rk_c = 2           ! number of outputs to the system
   !type(state_vector)     :: CT(rk_c)
   !real(wp)               :: Qc(rk_c,rk_c)
   !real(wp)               :: Rinv(rk_b,rk_b)

   ! Data matrices for RK lyap
   integer,  parameter    :: N = 2*nx           ! Number of grid points (excluding boundaries).
   !real(wp)               :: weight_mat(N,N)    ! integration weights matrix
   !real(wp)               :: weight_flat(N**2)    ! integration weights flat
   !real(wp)               :: BBTW(N,N)
   !real(wp)               :: CTCW(N,N)
   ! Data matrices for Riccati
   !real(wp)               :: CTQcCW(N,N)
   !real(wp)               :: BRinvBTW(N,N)

contains

subroutine initialize_parameters()
   implicit none
   ! Mesh array.
   real(wp), allocatable :: x(:)
   integer               :: i

   ! Construct mesh.
   x = linspace(-L/2, L/2, nx+2)
   dx = x(2)-x(1)

   ! Construct mu(x)
   mu(:) = (mu_0 - c_mu**2) + (mu_2 / 2.0_wp) * x(2:nx+1)**2

   ! Define integration weights
   weight          = dx

   print '(A)', ' ----------------------------------------'
   print '(A)', '    LINEAR GINZBURG LANDAU PARAMETERS'
   print '(A)', ' ----------------------------------------'
   print '(4X,A,F10.6," + ",F10.6," i")', 'nu    = ', nu
   print '(4X,A,F10.6," + ",F10.6," i")', 'gamma = ', gamma
   print '(4X,A,F10.6)', 'mu_0  = ', mu_0
   print '(4X,A,F10.6)', 'c_mu  = ', c_mu
   print '(4X,A,F10.6)', 'mu_2  = ', mu_2
   print '(A)', ' ----------------------------------------'

   return
end subroutine initialize_parameters

   !=========================================================
   !=========================================================
   !=====                                               =====
   !=====     LIGHTKRYLOV MANDATORY IMPLEMENTATIONS     =====
   !=====                                               =====
   !=========================================================
   !=========================================================

   !----------------------------------------------------
   !-----     TYPE-BOUND PROCEDURE FOR VECTORS     -----
   !----------------------------------------------------

   subroutine zero(self)
      class(state_vector), intent(inout) :: self
      self%state = 0.0_wp
   end subroutine zero

   real(wp) function dot(self, vec) result(alpha)
      ! weighted inner product
      class(state_vector),        intent(in) :: self
      class(abstract_vector_rdp), intent(in) :: vec
      select type(vec)
      type is(state_vector)
         alpha = dot_product(self%state, weight*vec%state)
      class default
         call stop_error('vec must be a state_vector', this_module, 'dot')
      end select
   end function dot
   
   integer function get_size(self) result(N)
     class(state_vector), intent(in) :: self
     N = 2*nx
   end function get_size

   subroutine scal(self, alpha)
      class(state_vector), intent(inout) :: self
      real(wp),            intent(in)    :: alpha
      self%state = self%state * alpha
   end subroutine scal

   subroutine axpby(alpha, vec, beta, self)
      class(state_vector),        intent(inout) :: self
      class(abstract_vector_rdp), intent(in)    :: vec
      real(wp),                   intent(in)    :: alpha, beta
      select type(vec)
      type is(state_vector)
         self%state = beta*self%state + alpha*vec%state
      class default
         call stop_error('vec must be a state_vector', this_module, 'axpby')
      end select
   end subroutine axpby

   subroutine rand(self, ifnorm)
      class(state_vector), intent(inout) :: self
      logical, optional,   intent(in)    :: ifnorm
      ! internals
      logical :: normalize
      real(wp) :: alpha
      real(wp), dimension(2*nx) :: mean, std
      normalize = optval(ifnorm,.true.)
      mean = 0.0_wp
      std  = 1.0_wp
      self%state = normal(mean,std)
      if (normalize) then
         alpha = self%norm()
         call self%scal(1.0/alpha)
      endif
   end subroutine rand

end module Ginzburg_Landau_Base
