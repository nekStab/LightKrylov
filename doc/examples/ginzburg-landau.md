---
title: Ginzburg-Landau
---

# **Leading eigenpairs of the linearized Ginzburg-Landau equation**

The linearized complex Ginzburg-Landau equation is a one-dimensional partial differential equation of convection-diffusion type.
It is given by
$$
\dfrac{\partial u}{\partial t} = \left( -\nu \dfrac{\partial}{\partial x} + \gamma \dfrac{\partial^2}{\partial x^2} + \mu(x) \right) u.
$$
The spatial domain considered is formaly infinite, hence there are no boundary conditions.
Yet, the complex-valued solution \(u(x, t)\) is required to remain finite as \(x \to \pm \infty\).
The complex convective velocity is defined as \(\nu = U + 2i c_u\), while the complex diffusion parameter is \(\gamma = 1 + i c_d\).
The spatially varying parameter \(\mu(x)\) is defined as
$$
\mu(x) = \left( \mu_0 - c_u^2 \right) + \mu_2 \dfrac{x^2}{2}.
$$
In the region given by
$$
-\sqrt{-2 \left( \mu_0 - c_u^2 \right)/2 } \leq x \leq \sqrt{-2 \left( \mu_0 - c_u^2 \right)/2 }
$$
perturbations are spatially amplified.
Outside of this region, they are exponentially attenuated.
This simple partial differential equation has been extensively used in the hydrodynamic stability community as a model for spatially developing flows, as well as a prototypical example for reduced-order modeling and controller design.
For more details, interested readers are refered to [[1](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.60.25), [2](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.78.4387), [3](https://www.annualreviews.org/content/journals/10.1146/annurev.fluid.37.061903.175810), [4](https://epubs.siam.org/doi/10.1137/100787350)] and references therein.
Here, we use this seemingly simple test case to illustrate how one can compute the leading eigenpairs of this linear operator with `LightKrylov` using a *time-stepper* approach.

**Time-stepping and eigenvalue analysis -** Suppose the equations have been discretized in space such that we now consider the finite-dimensional continuous-time linear time invariant dynamical system
$$
\dot{\mathbf{x}} = \mathbf{Ax},
$$
where \(\mathbf{x} \in \mathbb{C}^n\) is the state vector of the system, and \(\mathbf{A} \in \mathbb{C}^{n \times n}\) its dynamics matrix.
When studying the linear stability of such systems, one is often interested in the eigenvalues of \(\mathbf{A}\) having small real parts (i.e. in the vicinity of the imaginary axis).
While direct solvers such as `QZ` can be used when the dimension `n` is relatively small, one often needs to turn to iterative solvers for large-scale systems.
Unfortunately, techniques such as the Arnoldi iteration tend to converge first the eigenvalues having the largest magnitude.
In order to overcome this limitation, a possibly remedy is to use a *shift-invert* strategy.
Such a strategy however requires the inversion of a possibly very ill-conditioned linear system at each iteration.
Another alternative, pioneered by [Laurette Tuckerman](https://blog.espci.fr/laurette/), is to consider the following eigenvalue problem
$$
\exp(\tau \mathbf{A}) \mathbf{v}_i = \mu_i \mathbf{v}_i,
$$
where \(\mathbf{v}_i\) is an eigenvector of both \(\mathbf{A}\) and \(\exp(\tau \mathbf{A})\) and the eigenvalues \(\mu_i\) are related to those of \(\mathbf{A}\) via a simple logarithmic transformation.
Such a reformulation of the problem is advantageous for two main reasons:

- The eigenvalues of \( \mathbf{A} \) close to the imaginary axis get mapped along the unit-circle while those with very negative real parts get mapped close to zero. As a consequence, in the transformed complex plane, the eigenvalues of interest become the outer-most ones of the exponential matrix and thus the first to converge when using an iterative procedure.
- When using an interative procedure, only the matrix-vector product is needed. For a matrix exponential, this product can be easily and efficiently computed using a standard time-integration scheme. Consequently, it requires only a minimal set of modifications to adapt an existing simulation code for it to perform linear stability analyses.

Over the past two decades, this time-stepper approach to stability analysis has been instrumental in the field of hydrodynamic stability analysis.
For more details, please refer to [5, 6, 7, 8, 9].

## **Computing the leading eigenpairs of the Ginzburg-Landau equation via time-stepping**

@note
Owing to the small size of the problem, computing the leading eigenpairs of the linearized Ginzburg-Landau operator using a time-stepper approach is not the recommended way to proceed. The code is only used for illustratory purposes.
@endnote

The performance-critical part in a time-stepper approach for eigenvalue computation is the evaluation of the exponential matrix \(\exp(\tau \mathbf{A})\) applied to a given vector \(\mathbf{x}\).
Note that the exponential matrix does not have to be computed using classical dens elinear algebra.
Instead, the \(\mathbf{y} = \exp(\tau \mathbf{A}) \mathbf{x} \) can be computed by simply time-marching the linearized equations.

**Numerical setup -** The computational domain consider is \(-100 \leq x \leq 100\) along with homogeneous Dirichlet boundary conditions. This domain is large enough such that finite-length effects are negligible.
The spatial derivatives are discretized using second-order accurate finite differences on a uniform grid. The routine below then uses the Runge-Kutta class provided by [`rklib`](https://github.com/jacobwilliams/rklib) to implement the right-hand side of the linearized discretized system.

```fortran
subroutine rhs(me, t, x, f)
    ! Time-integrator.
    class(rk_class), intent(inout)           :: me
    ! Current time.
    real(kind=dp), intent(in)                :: t
    ! State vector.
    real(kind=dp), dimension(:), intent(in)  :: x
    ! Time-derivative.
    real(kind=dp), dimension(:), intent(out) :: f

    ! Internal variables.
    integer :: i
    real(kind=dp), dimension(nx) :: u, du, v, dv
    real(kind=dp)                :: d2u, d2v, cu, cv

    ! Sets the internal variables.
    f = 0.0_dp ; u = x(1:nx); du = 0.0_dp
    v = x(nx + 1:2*nx) ; dv = 0.0_dp

    !---------------------------------------------------
    !-----     Linear Ginzburg Landau Equation     -----
    !---------------------------------------------------

    ! Left most boundary points.
    cu = u(2)/(2*dx); cv = v(2)/(2*dx)
    du(1) = -(real(nu)*cu - aimag(nu)*cv) ! Convective term.
    dv(1) = -(aimag(nu)*cu + real(nu)*cv) ! Convective term.

    d2u = (u(2) - 2*u(1))/dx**2; d2v = (v(2) - 2*v(1))/dx**2
    du(1) = du(1) + real(gamma)*d2u - aimag(gamma)*d2v ! Diffusion term.
    dv(1) = dv(1) + aimag(gamma)*d2u + real(gamma)*d2v ! Diffusion term.

    du(1) = du(1) + mu(1)*u(1) ! Non-parallel term.
    dv(1) = dv(1) + mu(1)*v(1) ! Non-parallel term.

    ! Interior nodes.
    do i = 2, nx - 1
        ! Convective term.
        cu = (u(i + 1) - u(i - 1))/(2*dx)
        cv = (v(i + 1) - v(i - 1))/(2*dx)
        du(i) = -(real(nu)*cu - aimag(nu)*cv)
        dv(i) = -(aimag(nu)*cu + real(nu)*cv)

        ! Diffusion term.
        d2u = (u(i + 1) - 2*u(i) + u(i - 1))/dx**2
        d2v = (v(i + 1) - 2*v(i) + v(i - 1))/dx**2
        du(i) = du(i) + real(gamma)*d2u - aimag(gamma)*d2v
        dv(i) = dv(i) + aimag(gamma)*d2u + real(gamma)*d2v

        ! Non-parallel term.
        du(i) = du(i) + mu(i)*u(i)
        dv(i) = dv(i) + mu(i)*v(i)
    end do

    ! Right most boundary points.
    cu = -u(nx - 1)/(2*dx); cv = -v(nx - 1)/(2*dx)
    du(nx) = -(real(nu)*cu - aimag(nu)*cv) ! Convective term.
    dv(nx) = -(aimag(nu)*cu + real(nu)*cv) ! Convective term.

    d2u = (-2*u(nx) + u(nx - 1))/dx**2; d2v = (-2*v(nx) + v(nx - 1))/dx**2
    du(nx) = du(nx) + real(gamma)*d2u - aimag(gamma)*d2v ! Diffusion term.
    dv(nx) = dv(nx) + aimag(gamma)*d2u + real(gamma)*d2v ! Diffusion term.

    du(nx) = du(nx) + mu(nx)*u(nx) ! Non-parallel term.
    dv(nx) = dv(nx) + mu(nx)*v(nx) ! Non-parallel term.

    ! Copy results to the output array.
    f(1:nx) = du ; f(nx + 1:2*nx) = dv
end subroutine rhs
```

It can then be embeded into another subroutine having the interface required for the definition of an abstract linear operator (see [[LightKrylov_AbstractLinops(module):matvec(interface)]] for more details).
The definition the `state_vector` derived-type extended from `abstract_vector_cdp` is not presented here for the sake of conciseness.
It amounts to a simple wrapper around the standard Fortran rank-1 arrays.

```fortran
subroutine direct_solver(self, vec_in, vec_out)
    ! Linear Operator.
    class(exponential_prop), intent(inout)  :: self
    ! Input vector.
    class(abstract_vector_cdp), intent(in)  :: vec_in
    ! Output vector.
    class(abstract_vector_cdp), intent(out) :: vec_out

    ! Time-integrator.
    type(rks54_class) :: prop
    real(kind=dp)     :: dt = 1.0_dp
    real(kind=dp)     :: state_ic(2*nx), state_fc(2*nx)

    select type (vec_in)
    type is (state_vector)
        select type (vec_out)
        type is (state_vector)
            ! Get state vector.
            state_ic(:nx) = vec_in%state%re
            state_ic(nx + 1:) = vec_in%state%im
            ! Initialize propagator.
            call prop%initialize(n=2*nx, f=rhs)
            ! Integrate forward in time.
            call prop%integrate(0.0_dp, state_ic, dt, self%tau, state_fc)
            ! Pass-back the state vector.
            vec_out%state%re = state_fc(:nx)
            vec_out%state%im = state_fc(nx + 1:)
        end select
    end select
end subroutine direct_solver
```

Provided a similar set of routines have been implemented for the adjoint system, one can then extend the `abstract_linop_cdp` type to define the exponential propagator as illustrated below.

```fortran
type, extends(abstract_linop_cdp), public :: exponential_prop
    real(kind=dp), public :: tau ! Integration time.
contains
    private
    procedure, pass(self), public :: matvec => direct_solver
    procedure, pass(self), public :: rmatvec => adjoint_solver
end type exponential_prop
```

This is pretty much everything you need to implement in order to benefit from the iterative eigenvalue and singular value solvers provided by `LightKrylov`.
The Fortran program below uses the derived-types we've just defined and pass them to [[LightKrylov_IterativeSolvers(module):eigs(interface)]] in order to compute the leading eigenpairs of the linearized Ginzburg-Landau operator.

```fortran
program demo
    use stdlib_io_npy, only: save_npy
    use LightKrylov, only: wp => dp
    use LightKrylov
    use Ginzburg_Landau
    implicit none

    !> Exponential propagator.
    type(exponential_prop), allocatable :: A
    !> Sampling time.
    real(kind=wp), parameter :: tau = 0.1_wp

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
    integer       :: i
    complex(wp) :: eigenvectors(nx, nev)

    !> Initialize physical parameters.
    call initialize_parameters()

    !> Initialize exponential propagator.
    A = exponential_prop(tau)

    !> Initialize Krylov subspace.
    allocate (X(nev)); call zero_basis(X)

    !> Call to LightKrylov.
    call eigs(A, X, lambda, residuals, info)

    !> Transform eigenspectrum from unit-disk to standard complex plane.
    lambda = log(lambda)/tau

    !> Extract the leading eigenvectors from the Krylov basis into a standard array.
    do i = 1, nev
        eigenvectors(:, i) = X(i)%state
    end do

    !> Save eigenvalues and eigenvectors to disk.
    call save_eigenspectrum(lambda, residuals, "example/ginzburg_landau/eigenspectrum.npy")
    call save_npy("example/ginzburg_landau/eigenvectors.npy", eigenvectors)
end program demo
```

From the root directory of `LightKrylov`, you can run this example by typing the following command:

```bash
fpm run --example Ginzburg-Landau
```

Then `python eigenplots.py` in the `example/ginzburg-landau` directory will produce the following figure depicting the eigenspectrum and leading eigenvector of the discretized linearized Ginzburg-Landau operator.

<div style="text-align:center;">
    <img src="ginzburg_landau_modal_analysis.png" style="width:80%;">
    </br>
    <b>Figure:</b> Left panel depicts the eigenspectrum computed with `LightKrylov` and compares against the analytical expressions. The branch splitting observed for very stable eigenvalues result from the discretization of the operator. Right panel depicts the real and imaginary parts of the leading eigenvector.
</div>

## **References**

[1] J.-M. Chomaz, P. Huerre and L. G. Redekopp. [*Bifurcations to local and global modes in spatially developping flows*](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.60.25). Phys. Rev. Lett., 60 (25), 1988.

[2] C. Cossu and J.-M. Chomaz. [*Global measures of local convective instabilities*](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.78.4387). Phys. Rev. Lett., 78 (4387), 1997.

[3] J.-M. Chomaz. [*Global instabilities in spatially developing flow: non-normality and nonlinearity*](https://www.annualreviews.org/content/journals/10.1146/annurev.fluid.37.061903.175810). Ann. Rev. Fluid Mech., 37 (357-392), 2005.

[4] M. Ilak, S. Bagheri, L. Brandt, C. W. Rowley, and D. S. Henningson. [*Model reduction of the nonlinear complex Ginzburg-Landau equation*](https://epubs.siam.org/doi/10.1137/100787350). SIAM Journal on Applied Dynamical Sysmte, 9 (4), 2010.

[5] W. S. Edward, L. S. Tuckerman, R. A. Friesner and D. C. Sorensen. [*Krylov methods for the incompressible Navier-Stokes equations](https://www.sciencedirect.com/science/article/pii/S0021999184710072). Journal of Computational Physics, 110 (1), 1994.

[6] Tuckerman, Laurette S., and Dwight Barkley. *Bifurcation analysis for timesteppers.* Numerical methods for bifurcation problems and large-scale dynamical systems, 2000.

[7] S. Bagheri, E. Akervik, L. Brandt, and D. S. Henningson. [*Matrix-free methods for stability and control of boundary layers*](https://arc.aiaa.org/doi/abs/10.2514/1.41365). AIAA Journal, 47 (5), 2009.

[8] J.-Ch. Loiseau, J.-Ch. Robinet, S. Cherubini and E. Leriche. [*Investigation of the roughness-induced transition: global stability analyses and direct numerical similations*](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/investigation-of-the-roughnessinduced-transition-global-stability-analyses-and-direct-numerical-simulations/0D62A95BD1FD2D2C664A8EDABCFAFC46). Journal of Fluid Mechanics, 2014.

[9] R. A. S. Frantz, J.-Ch. Loiseau and J.-Ch. Robinet. [*Krylov methods for large-scale dynamical systems: applications in fluid dynamics*](https://asmedigitalcollection.asme.org/appliedmechanicsreviews/article-abstract/75/3/030802/1156502/Krylov-Methods-for-Large-Scale-Dynamical-Systems). Applied Mechanics Review, 75 (3), 2023.
