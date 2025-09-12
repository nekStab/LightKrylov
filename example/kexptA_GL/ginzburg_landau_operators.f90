module Ginzburg_Landau_Operators
   ! Standard Library.
   use stdlib_optval, only : optval
   ! RKLIB module for time integration.
   use rklib_module
   ! LightKrylov for linear algebra.
   use LightKrylov
   use LightKrylov, only : wp => dp
   use LightKrylov_utils, only : assert_shape
   ! Ginzburg Landau
   use Ginzburg_Landau_Base
   implicit none

   private :: this_module
   character(len=*), parameter :: this_module = 'Ginzburg_Landau_Operators'
   public  :: direct_GL, adjoint_GL
   public  :: exptA

   !--------------------------------------
   !-----     LINEAR GL OPERATOR     -----
   !--------------------------------------

   type, extends(abstract_linop_rdp), public :: GL_operator
   contains
      private
      procedure, pass(self), public :: matvec  => direct_matvec_GL
      procedure, pass(self), public :: rmatvec => adjoint_matvec_GL
   end type GL_operator

   !-------------------------------------------
   !-----     EXPONENTIAL PROPAGATORS     -----
   !-------------------------------------------

   type, extends(abstract_expta_linop_rdp), public :: exponential_prop_rk
      real(wp) :: tol = 1e-8_wp
      integer :: n_accepted
      integer :: n_rejected
      real(wp) :: last_dt
   contains
      private
      procedure, pass(self), public :: matvec => direct_rk_solver
      procedure, pass(self), public :: rmatvec => adjoint_rk_solver
   end type exponential_prop_rk

   type, extends(abstract_expta_linop_rdp), public :: exponential_prop_krylov
      class(abstract_linop_rdp), allocatable :: A
      real(wp) :: tol = 1e-8_wp
      integer :: kdim = 30
      integer :: info
   contains
      private
      procedure, pass(self), public :: matvec => direct_krylov_solver
      procedure, pass(self), public :: rmatvec => adjoint_krylov_solver
   end type exponential_prop_krylov

   interface exponential_prop_krylov
      module function construct_krylov_prop(A) result(prop)
         class(abstract_linop_rdp), intent(in) :: A
         type(exponential_prop_krylov) :: prop
      end function
   end interface

contains

   module procedure construct_krylov_prop
      prop%A = A
   end procedure

   !========================================================================
   !========================================================================
   !=====                                                              =====
   !=====     PHYSICAL MODEL : LINEARIZED GINZBURG-LANDAU EQUATION     =====
   !=====                                                              =====
   !========================================================================
   !========================================================================

   !---------------------------------------------------------
   !-----      LINEARIZED GINZBURG-LANDAU EQUATIONS     -----
   !---------------------------------------------------------

   subroutine direct_GL(vec_in, vec_out)

      !> State vector.
      real(wp), dimension(:), intent(in)  :: vec_in
      !> Time-derivative.
      real(wp), dimension(:), intent(out) :: vec_out

      !> Internal variables.
      integer                 :: i
      real(wp), dimension(nx) :: u, v, du, dv
      real(wp)                :: d2u, d2v, cu, cv

      u = vec_in(1:nx)     
      v = vec_in(nx+1:2*nx)

      !---------------------------------------------------
      !-----     Linear Ginzburg Landau Equation     -----
      !---------------------------------------------------

      cu = u(2) / (2*dx) ; cv = v(2) / (2*dx)
      du(1) = -(real(nu)*cu - aimag(nu)*cv) ! Convective term.
      dv(1) = -(aimag(nu)*cu + real(nu)*cv) ! Convective term.

      d2u = (u(2) - 2*u(1)) / dx**2 ; d2v = (v(2) - 2*v(1)) / dx**2
      du(1) = du(1) + real(gamma)*d2u - aimag(gamma)*d2v ! Diffusion term.
      dv(1) = dv(1) + aimag(gamma)*d2u + real(gamma)*d2v ! Diffusion term.

      du(1) = du(1) + mu(1)*u(1) ! Non-parallel term.
      dv(1) = dv(1) + mu(1)*v(1) ! Non-parallel term.

      ! Interior nodes.
      do i = 2, nx-1
         ! Convective term.
         cu = (u(i+1) - u(i-1)) / (2*dx)
         cv = (v(i+1) - v(i-1)) / (2*dx)
         du(i) = -(real(nu)*cu - aimag(nu)*cv)
         dv(i) = -(aimag(nu)*cu + real(nu)*cv)

         ! Diffusion term.
         d2u = (u(i+1) - 2*u(i) + u(i-1)) / dx**2
         d2v = (v(i+1) - 2*v(i) + v(i-1)) / dx**2
         du(i) = du(i) + real(gamma)*d2u - aimag(gamma)*d2v
         dv(i) = dv(i) + aimag(gamma)*d2u + real(gamma)*d2v

         ! Non-parallel term.
         du(i) = du(i) + mu(i)*u(i)
         dv(i) = dv(i) + mu(i)*v(i)
      enddo

      ! Right most boundary points.
      cu = -u(nx-1) / (2*dx) ; cv = -v(nx-1) / (2*dx)
      du(nx) = -(real(nu)*cu - aimag(nu)*cv) ! Convective term.
      dv(nx) = -(aimag(nu)*cu + real(nu)*cv) ! Convective term.

      d2u = (-2*u(nx) + u(nx-1)) / dx**2 ; d2v = (-2*v(nx) + v(nx-1)) / dx**2
      du(nx) = du(nx) + real(gamma)*d2u - aimag(gamma)*d2v ! Diffusion term.
      dv(nx) = dv(nx) + aimag(gamma)*d2u + real(gamma)*d2v ! Diffusion term.

      du(nx) = du(nx) + mu(nx)*u(nx) ! Non-parallel term.
      dv(nx) = dv(nx) + mu(nx)*v(nx) ! Non-parallel term.

      vec_out(1:nx)      = du
      vec_out(nx+1:2*nx) = dv

   end subroutine direct_GL

   !-----------------------------------------------------------
   !-----     Adjoint linear Ginzburg-Landau equation     -----
   !-----------------------------------------------------------

   subroutine adjoint_GL(vec_in, vec_out)
      !> State vector.
      real(wp), dimension(:), intent(in)  :: vec_in
      !> Time-derivative.
      real(wp), dimension(:), intent(out) :: vec_out

      ! Internal variables.
      integer :: i
      real(wp), dimension(nx) :: u, du
      real(wp), dimension(nx) :: v, dv
      real(wp)                :: d2u, d2v, cu, cv

      ! Sets the internal variables.
      u = vec_in(1:nx)     
      v = vec_in(nx+1:2*nx)

      !---------------------------------------------------
      !-----     Linear Ginzburg Landau Equation     -----
      !---------------------------------------------------

      ! Left most boundary points.
      cu = u(2) / (2*dx) ; cv = v(2) / (2*dx)
      du(1) = (real(nu)*cu + aimag(nu)*cv) ! Convective term.
      dv(1) = (-aimag(nu)*cu + real(nu)*cv) ! Convective term.

      d2u = (u(2) - 2*u(1)) / dx**2 ; d2v = (v(2) - 2*v(1)) / dx**2
      du(1) = du(1) + real(gamma)*d2u + aimag(gamma)*d2v ! Diffusion term.
      dv(1) = dv(1) - aimag(gamma)*d2u + real(gamma)*d2v ! Diffusion term.

      du(1) = du(1) + mu(1)*u(1) ! Non-parallel term.
      dv(1) = dv(1) + mu(1)*v(1) ! Non-parallel term.

      ! Interior nodes.
      do i = 2, nx-1
         ! Convective term.
         cu = (u(i+1) - u(i-1)) / (2*dx)
         cv = (v(i+1) - v(i-1)) / (2*dx)
         du(i) = (real(nu)*cu + aimag(nu)*cv)
         dv(i) = (-aimag(nu)*cu + real(nu)*cv)

         ! Diffusion term.
         d2u = (u(i+1) - 2*u(i) + u(i-1)) / dx**2
         d2v = (v(i+1) - 2*v(i) + v(i-1)) / dx**2
         du(i) = du(i) + real(gamma)*d2u + aimag(gamma)*d2v
         dv(i) = dv(i) - aimag(gamma)*d2u + real(gamma)*d2v

         ! Non-parallel term.
         du(i) = du(i) + mu(i)*u(i)
         dv(i) = dv(i) + mu(i)*v(i)
      enddo

      ! Right most boundary points.
      cu = -u(nx-1) / (2*dx) ; cv = -v(nx-1) / (2*dx)
      du(nx) = (real(nu)*cu + aimag(nu)*cv) ! Convective term.
      dv(nx) = (-aimag(nu)*cu + real(nu)*cv) ! Convective term.

      d2u = (-2*u(nx) + u(nx-1)) / dx**2 ; d2v = (-2*v(nx) + v(nx-1)) / dx**2
      du(nx) = du(nx) + real(gamma)*d2u + aimag(gamma)*d2v ! Diffusion term.
      dv(nx) = dv(nx) - aimag(gamma)*d2u + real(gamma)*d2v ! Diffusion term.

      du(nx) = du(nx) + mu(nx)*u(nx) ! Non-parallel term.
      dv(nx) = dv(nx) + mu(nx)*v(nx) ! Non-parallel term.

      ! Copy results to the output array.
      vec_out(1:nx)      = du
      vec_out(nx+1:2*nx) = dv

   end subroutine adjoint_GL

   !--------------------------------------
   !-----     WRAPPERS FOR RKLIB     -----
   !--------------------------------------

   subroutine rhs(me, t, x, f)
      ! Time-integrator.
      class(rk_class),               intent(inout) :: me
      ! Current time.
      real(wp),                      intent(in)    :: t
      ! State vector.
      real(wp),        dimension(:), intent(in)    :: x
      ! Time-derivative.
      real(wp),        dimension(:), intent(out)   :: f

      f = 0.0_wp
      call direct_GL(x, f)

   end subroutine rhs

   subroutine adjoint_rhs(me, t, x, f)
      ! Time-integrator.
      class(rk_class),               intent(inout) :: me
      ! Current time.
      real(wp),                      intent(in)    :: t
      ! State vector.
      real(wp),        dimension(:), intent(in)    :: x
      ! Time-derivative.
      real(wp),        dimension(:), intent(out)   :: f

      f = 0.0_wp
      call adjoint_GL(x, f)

   end subroutine adjoint_rhs

   !-------------------------------------------------------------
   !-----     TYPE-BOUND PROCEDURES FOR THE GL OPERATOR     -----
   !-------------------------------------------------------------

   subroutine direct_matvec_GL(self, vec_in, vec_out)
      !> Linear Operator.
      class(GL_operator),          intent(inout)  :: self
      !> Input vector.
      class(abstract_vector_rdp), intent(in)  :: vec_in
      !> Output vector.
      class(abstract_vector_rdp), intent(out) :: vec_out
      select type(vec_in)
      type is (state_vector)
         select type(vec_out)
         type is (state_vector)
            call direct_GL(vec_in%state, vec_out%state)
         class default
            call type_error('vec_out', 'state_vector', this_module, 'direct_matvec_GL')
         end select
      class default
         call type_error('vec_in', 'state_vector', this_module, 'direct_matvec_GL')
      end select
   end subroutine direct_matvec_GL

   subroutine adjoint_matvec_GL(self, vec_in, vec_out)
      !> Linear Operator.
      class(GL_operator),         intent(inout)  :: self
      !> Input vector.
      class(abstract_vector_rdp), intent(in)  :: vec_in
      !> Output vector.
      class(abstract_vector_rdp), intent(out) :: vec_out
      select type(vec_in)
      type is (state_vector)
         select type(vec_out)
         type is (state_vector)
            call adjoint_GL(vec_in%state, vec_out%state)
         class default
            call type_error('vec_out', 'state_vector', this_module, 'adjoint_matvec_GL')
         end select
      class default
         call type_error('vec_in', 'state_vector', this_module, 'adjoint_matvec_GL')
      end select
   end subroutine adjoint_matvec_GL

   !------------------------------------------------------------------------
   !-----     TYPE-BOUND PROCEDURES FOR THE EXPONENTIAL PROPAGATOR     -----
   !------------------------------------------------------------------------

   subroutine direct_rk_solver(self, vec_in, vec_out)
      ! Linear Operator.
      class(exponential_prop_rk),  intent(inout)  :: self
      ! Input vector.
      class(abstract_vector_rdp),  intent(in)  :: vec_in
      ! Output vector.
      class(abstract_vector_rdp),  intent(out) :: vec_out

      ! Time-integrator.
      type(rks54_class) :: prop
      real(wp)          :: dt = 1.0_wp

      select type(vec_in)
      type is(state_vector)
         select type(vec_out)
         type is(state_vector)

            ! Initialize propagator.
            call prop%initialize(n=2*nx, f=rhs, rtol=[self%tol], atol=[self%tol])
            ! Integrate forward in time.
            call prop%integrate(0.0_wp, vec_in%state, dt, self%tau, vec_out%state)
            ! get info
            call prop%info(self%n_accepted, self%n_rejected, self%last_dt)

         class default
            call type_error('vec_out', 'state_vector', this_module, 'direct_rk_solver')
         end select
      class default
         call type_error('vec_in', 'state_vector', this_module, 'direct_rk_solver')
      end select
   end subroutine direct_rk_solver

   subroutine adjoint_rk_solver(self, vec_in, vec_out)
      ! Linear Operator.
      class(exponential_prop_rk),  intent(inout)  :: self
      ! Input vector.
      class(abstract_vector_rdp),  intent(in)  :: vec_in
      ! Output vector.
      class(abstract_vector_rdp),  intent(out) :: vec_out

      ! Time-integrator.
      type(rks54_class) :: prop
      real(wp)          :: dt = 1.0_wp

      select type(vec_in)
      type is(state_vector)
         select type(vec_out)
         type is(state_vector)

            ! Initialize propagator.
            call prop%initialize(n=2*nx, f=adjoint_rhs, rtol=[self%tol], atol=[self%tol])
            ! Integrate forward in time.
            call prop%integrate(0.0_wp, vec_in%state, dt, self%tau, vec_out%state)
            ! get info
            call prop%info(self%n_accepted, self%n_rejected, self%last_dt)

         class default
            call type_error('vec_out', 'state_vector', this_module, 'adjoint_rk_solver')
         end select
      class default
         call type_error('vec_in', 'state_vector', this_module, 'adjoint_rk_solver')
      end select
   end subroutine adjoint_rk_solver

   subroutine direct_krylov_solver(self, vec_in, vec_out)
      ! Linear Operator.
      class(exponential_prop_krylov), intent(inout)  :: self
      ! Input vector.
      class(abstract_vector_rdp),     intent(in)  :: vec_in
      ! Output vector.
      class(abstract_vector_rdp),     intent(out) :: vec_out

      ! Time-integrator.
      select type(vec_in)
      type is(state_vector)
         select type(vec_out)
         type is(state_vector)

            call kexpm(vec_out, self%A, vec_in, self%tau, self%tol, self%info, trans=.false., kdim=self%kdim)

         class default
            call type_error('vec_out', 'state_vector', this_module, 'direct_krylov_solver')
         end select
      class default
         call type_error('vec_in', 'state_vector', this_module, 'direct_krylov_solver')
      end select
   end subroutine direct_krylov_solver

   subroutine adjoint_krylov_solver(self, vec_in, vec_out)
      ! Linear Operator.
      class(exponential_prop_krylov), intent(inout)  :: self
      ! Input vector.
      class(abstract_vector_rdp),     intent(in)  :: vec_in
      ! Output vector.
      class(abstract_vector_rdp),     intent(out) :: vec_out

      ! Time-integrator.
      select type(vec_in)
      type is(state_vector)
         select type(vec_out)
         type is(state_vector)

            call kexpm(vec_out, self%A, vec_in, self%tau, self%tol, self%info, trans=.true., kdim=self%kdim)

         class default
            call type_error('vec_out', 'state_vector', this_module, 'adjoint_krylov_solver')
         end select
      class default
         call type_error('vec_in', 'state_vector', this_module, 'adjoint_krylov_solver')
      end select
   end subroutine adjoint_krylov_solver

   !--------------------------------------
   !-----     EXP(tA) SUBROUTINE     -----
   !--------------------------------------

   subroutine exptA(vec_out, prop, vec_in, tau, info, trans)
      !! Subroutine for the exponential propagator that conforms with the abstract interface
      !! defined in expmlib.f90
      class(abstract_vector_rdp),      intent(out)   :: vec_out
      !! Output vector
      class(abstract_expta_linop_rdp), intent(inout) :: prop
      !! Linear operator
      class(abstract_vector_rdp),      intent(in)    :: vec_in
      !! Input vector.
      real(wp),                        intent(in)    :: tau
      !! Integration horizon
      integer,                         intent(out)   :: info
      !! Information flag
      logical, optional,               intent(in)    :: trans
      logical                                        :: transpose
      !! Direct or Adjoint?

      ! optional argument
      transpose = optval(trans, .false.)

      ! time integrator
      select type (vec_in)
      type is (state_vector)
         select type (vec_out)
         type is (state_vector)
            prop%tau = tau
            if (transpose) then
               call prop%rmatvec(vec_in, vec_out)
            else
               call prop%matvec(vec_in, vec_out)
            end if 
         class default
            call type_error('vec_out', 'state_vector', this_module, 'exptA')
         end select
      class default
         call type_error('vec_in', 'state_vector', this_module, 'exptA')
      end select
   end subroutine exptA

end module Ginzburg_Landau_Operators