---
title: '`LightKrylov`: Lightweight implementation of Krylov subspace techniques in modern Fortran'
tags:
  - Fortran
  - Numerical linear algebra
  - Krylov methods
  - Sparse linear systems
  - Eigenvalues and singular values
authors:
  - name: J. Simon Kern
    orcid: 0000-0002-2460-578X
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 1
  - name: Ricardo S. Frantz
    orcid: 0000-0001-8219-3368
    equal-contrib: true
    affiliation: 1
  - name: Jean-Christophe Loiseau
    orcid: 0000-0002-7244-8416
    corresponding: true
    equal-contrib: true
    affiliation: 1 # (Multiple affiliations must be quoted)
affiliations:
 - name: Arts et Métiers Institute of Technology
   index: 1
date: 13 August 2017
bibliography: paper.bib
---

# Summary

[`LightKrylov`](https://github.com/nekStab/LightKrylov) is a Fortran package with minimal dependencies providing a collection of Krylov-based methods for large-scale linear algebra.
Direct methods tend to have a computational cost scaling as $\mathcal{O}(n^3)$ making them impractical for large-scale problems.
Even when $A$ is sparse, computing matrix factorizations can incur an $\mathcal{O}(n^2)$ storage cost in the worst case by introducing fill-in.
In contrast, Krylov methods only need a function computing the matrix-vector product $u \leftarrow Av$ (or $u \leftarrow A^H v$) to iteratively construct a *Krylov subspace* [@krylov-1931] from which the solutions can be extracted, see @ipsen-1998, @saad-2003 and to @frantz-2023.

# Statement of need

## A collection of Krylov-based algorithms in pure modern Fortran

`LightKrylov` provides Fortran users with SciPy-inspired interfaces to widely used Krylov techniques, including:

- **Krylov processes:** Arnoldi factorization, Golub-Kahan bidiagonalization, as well as Lanczos tridiagonalization for Hermitian operators.
- **Krylov methods:**
    - **Linear systems -** Conjugate Gradient (CG), Generalized Minimal Residual method (GMRES) and Flexible GMRES.
    - **Spectral decomposition -** Arnoldi iteration with Krylov-Schur restart for non-Hermitian operators, Lanczos tridiagonalization with thick restart for Hermitian ones.
    - **SVD -** restarted Golub-Kahan bidiagonalization.

Libraries exposing more methods exist in other languages, e.g., `Krylov.jl` [@montoison-2023] in Julia or `PETSc` [@petsc-web-page].
Integrating multi-language libraries into existing Fortran codes often presents significant challenges (e.g. dependency management or potential performance overhead for instance).
In contrast, `LightKrylov` provides a pure Fortran alternative, compliant with the Fortran 2018 standard, and only requires the Fortran standard library `stdlib` as an external runtime dependency.

Additionally, using `fypp` [@fypp-webpage] as a build-time dependency significantly reduces code duplication and maintenance effort, ensuring consistency across different data types (`single`, `double`, `extended`, and `quadruple precision`) both for `real` and `complex` numbers while keeping the core algorithmic logic centralized.
Finally, its build process relies on the Fortran package manager `fpm`, greatly facilitating its installation and integration with the modern Fortran ecosystem.

## A focus on abstract linear operators and abstract vectors

Krylov methods can be implemented without explicit reference to the data structure used to represent a vector or a linear operator, nor to how the actual matrix-vector product is implemented.
To do so, `LightKrylov` uses modern Fortran `abstract` type constructs.
A stripped-down version of the abstract vector type is shown below.

```fortran
type, abstract :: abstract_vector_rdp
  contains
    ! Abstract procedure to compute the scalar-vector product.
    procedure(abstract_scal_rdp), pass(self), deferred :: scal
    ! Abstract procedure to compute y = alpha*x + beta*y.
    procedure(abstract_axpby_rdp), pass(self), deferred :: axpby
    ! Abstract procedure to compute the vector dot product.
    procedure(abstract_dot_rdp), pass(self), deferred :: dot
end type
```

These type-bound procedures cover the basic operations on vectors: scalar-vector product, linear combination of two vectors, and the dot product. These operations are the essential building blocks required by Krylov algorithms.
Their signatures follow, to the extent possible, the BLAS standard for a more familiar use.
For instance, the `abstract_axpby_rdp` interface reads

```fortran
abstract interface
  subroutine abstract_axpby_rdp(alpha, vec, beta, self)
    double precision, intent(in)    :: alpha, beta
    class(abstract_vector_rdp), intent(in)    :: vec
    class(abstract_vector_rdp), intent(inout) :: self
  end subroutine
end interface
```
mimicking the signature of the (extended) BLAS-1 subroutine `axpby`.

Similarly, a stripped-down version of an abstract linear operator type is shown below.

```fortran
type, abstract :: abstract_linop_rdp
  contains
    ! Abstract procedure to compute the matrix-vector product y = A * x
    procedure(abstract_matvec_rdp), pass(self), deferred :: matvec
    ! Abstract procedure to compute the matrix-vector product y = A^H * x
    procedure(abstract_matvec_rdp), pass(self), deferred :: rmatvec
end type
```

The two type-bound procedures need to implement the matrix-vector product and the transposed matrix-vector product.
The corresponding abstract interface reads

```fortran
abstract interface
  subroutine abstract_matvec_rdp(self, vec_in, vec_out)
    class(abstract_linop_rdp), intent(inout) :: self
    class(abstract_vector_rdp), intent(in)    :: vec_in
    class(abstract_vector_rdp), intent(out)   :: vec_out
  end subroutine
end interface
```

Using such `abstract` types enables us to focus on the high-level implementation of the different algorithms while leaving the performance-critical details to the users.
After extending these `abstract` types for their particular applications, users can solve linear systems or compute eigenvalues as easily as

```fortran
! Solve linear system.
call gmres(A, b, x, info)

! Compute eigenvalues and eigenvectors.
call eigs(A, V, lambda, residuals, info)
```

In addition, `LightKrylov` exposes abstract interfaces for preconditioners, as well as a Newton-GMRES solver.
Examples illustrating how to extend these `abstract` types for computing the leading eigenpairs of the linearized Ginzburg-Landau equation or computing unstable periodic orbits of the Rössler system and studying their Lyapunov exponents can be found [here](https://github.com/nekStab/LightKrylov/tree/main/example).

# Performance in a production-ready open-source codebase

`LightKrylov` was successfully integrated into [`neklab`](https://github.com/nekStab/neklab), a toolbox for stability and bifurcation analysis using the high-performance spectral element solver `Nek5000`. The abstract vector interface allows direct use of `Nek5000`'s distributed data structures, and the pure-Fortran nature facilitated integration with its existing build system, demonstrating the library's suitability for large-scale HPC applications.

## Hydrodynamic stability of an unstable fixed point of the nonlinear Navier-Stokes equations

Computing (stable or unstable) fixed points of the nonlinear governing equations is a necessary step when studying the stability properties of a nonlinear system. In hydrodynamic instability applications, this amounts to finding a stationary solution of the nonlinear Navier-Stokes equations. Investigating the spectral properties of the linearized Navier-Stokes operator about such fixed points is critical to understand the transition to unsteadiness and turbulence in many flow applications. 

Using the supercritical two-dimensional flow past a circular cylinder at Reynolds number of 100 as an example, we showcase the efficient integration of `LightKrylov` and `Nek5000` and validate the results with algorithms provided by the [`KTH Framework`](https://github.com/KTH-Nek5000/KTH_Framework) toolbox [@kth-framework] based on the same Navier-Stokes solver. In all cases, the governing equations are discretized with 1996 spectral elements using a 5th-order polynomial approximation in each direction, resulting in 71,856 grid points and roughly 175,000 degrees of freedom (i.e., the two velocity components and pressure). All computations were run in parallel on 12 Intel Core Ultra 7 processors and the numerical settings, in particular those internal to `Nek5000`, are identical for all runs.

The unstable fixed point of the nonlinear Navier-Stokes equations is obtained using both `LightKrylov`'s *time-stepper*-based Newton-GMRES solver [@frantz-2023] and the `KTH Framework` implementation of the selective frequency damping algorithm [@pof-sfd] (often considered as the gold standard in this community). The top row in \autoref{fig:timings} depicts the streamwise (*u*) velocity distribution at the fixed point (left) and the pointwise squared difference between the solutions obtained with the two solution techniques and a tolerance of $\varepsilon = 10^{-8}$ (right). Note that, while the definition of the error norm slightly differs between the methods, the one used by the Newton-GMRES solver is stricter. The code to run this example can be found [here](?).

For the linear stability analysis, the first complex-conjugate eigenpair is computed with the Krylov-Schur algorithm using `LightKrylov`'s abstract interfaces. In the lower-left panel of \autoref{fig:timings}, it is compared to the spectrum computed using `KTH Framework` relying under the hood on `ARPACK`. Both solvers effectively use a *time-stepper* approach, approximating the leading eigenvalues of the exponential propagator $\exp(\tau \mathbf{A})$ rather than $\mathbf{A}$ directly. The code to run this example can be found [here](?).

The table in the lower-right panel of \autoref{fig:timings} summarizes the wall-clock times of the `neklab` computations. Isolating the intrinsic cost of the algorithms in `LightKrylov` from the cost of the calls to LAPACK and the linear and nonlinear Navier-Stokes solvers (`matvec` and `response`, respectively) conclusively shows that the use of extended `abstract` types and object-oriented programming in Fortran incurs a negligible computational overhead for such large-scale applications.

![Validation of `LightKrylov`: Newton-GMRES and `eigs`.](LK_newton_eigs_timings.png){#fig:timings}

# Perspectives

Despite being in its early development stage, `LightKrylov` has already been used in hydrodynamic stability production runs on dozens of processors. It has also been interfaced with [`neko`](https://github.com/ExtremeFLOW/neko), a modernized implementation of `Nek5000` running on GPUs, and is currently being interfaced with [`dNami`](https://github.com/dNamiLab/dNami), a high-performance source code generator for hyperbolic partial differential equations.

Ongoing development efforts include:

- **Build system:** `LightKrylov` relies on the Fortran package manager `fpm` for its build. While it ensures nice integration into the modern Fortran ecosystem, it also limits its usage to codes built using `fpm`. An additional build system based on `CMake` for easier integration into non-`fpm` codes will be proposed shortly.
- **Krylov processes:** As of summer 2025, `LightKrylov` provides implementations for three of the most important Krylov processes. Future releases will include the non-Hermitian Lanczos and two-sided Arnoldi factorizations useful for reduced-order modeling of linear time-invariant (LTI) systems as well as the Saunders-Simon-Yip process [@ssy-krylov].
- **Iterative solvers:** Beyond `CG` and `GMRES`, the integration of new Krylov processes will enable us to provide an extended list of iterative solvers, including `MINRES` (indefinite Hermitian systems), `CGNE` and `LSQR` (least-squares problems), and `LSMR` (least-norm problems).

### Acknowledgements

We acknowledge the financial support of the French National Agency for Research (ANR) through the ANR-33-CE46-0008-CONMAN grant agreement. We would also like to thank the [fortran-lang](https://fortran-lang.org/) community for the development of [`stdlib`](https://stdlib.fortran-lang.org/), and in particular Frederico Perini, Jeremie Vandenplas, and Jose Alvez for their work on the `stdlib_linalg` module which `LightKrylov` relies on.

# References
