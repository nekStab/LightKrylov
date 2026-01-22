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

Direct solvers for linear algebraic systems scale cubically in the problem's dimension, rapidly becoming intractable for large-scale problems, while sparse factorization may still require quadratic storage due to fill-in.
Krylov techniques [@krylov-1931] avoid these costs by needing only a routine that computes a matrix-vector product, iteratively building a subspace from which the solution is obtained, see @ipsen-1998, @saad-2003 and @frantz-2023.
[`LightKrylov`](https://github.com/nekStab/LightKrylov) is a Fortran package providing a suite of such Krylov methods along with an easy-to-use high level API based on `abstract` types.
It is primarily intended for applications where the linear operator of interest is only available implicitly via a matrix-vector subroutine and enables users to maximally re-use existing components of their code base (including parallelisation), thus requiring a minimal set of changes without sacrificing computational performance.

# Statement of need

## A collection of Krylov-based algorithms in pure modern Fortran

`LightKrylov` provides Fortran users with SciPy-inspired interfaces to widely used Krylov techniques, including:

- **Linear systems -** Conjugate Gradient (CG), Generalized Minimal Residual method (GMRES) and Flexible GMRES [@saad:fgmres:siam].
- **Spectral decomposition -** Arnoldi method (with Krylov-Schur restart) for non-Hermitian operators, Lanczos tridiagonalisation for Hermitian ones.
- **SVD -** Golub-Kahan bidiagonalisation.

It is a pure Fortran package, compliant with the 2018 standard, and requiring only the community-led Fortran standard library [`stdlib`](https://stdlib.fortran-lang.org/) [@stdlib:ieee] as dependency.
Moreover, its build process relies on the Fortran package manager `fpm`, facilitating its integration with the modern Fortran ecosystem.

## A focus on abstract linear operators and abstract vectors

Krylov methods can be implemented without explicit reference to the data structure used to represent a vector or linear operator, nor to how the matrix-vector product is implemented.
To do so, `LightKrylov` uses modern Fortran `abstract` types.
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

The type-bound procedures cover the basic operations on vectors: scalar-vector product, linear combination of two vectors, and the dot product. These operations are the essential building blocks required by Krylov algorithms.
Their signatures follow, to the extent possible, the BLAS standard.
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
Abstract linear operators are defined similarly, with two type-bound procedures required to implement the matrix-vector and transpose matrix-vector product.
Using `abstract` types enables us to focus on the high-level implementation of the different algorithms while leaving the performance-critical details to the users.
In addition, `LightKrylov` exposes abstract types for preconditioners, as well as a Newton-GMRES solver for nonlinear systems.
After extending these `abstract` types for their application, one can solve linear systems or compute eigenvalues as easily as `call gmres(A, b, x, info)` or `call eigs(A, V, lambda, residuals, info)`.

## High-level comparison with other libraries

`PETSc` [@petsc-web-page] is a widely used library for large-scale linear algebra problems, especially those resulting from PDE discretisation.
While it offers more than Krylov methods, its many data structures and build process can make integration difficult into an already existing large code base when only linear solvers are needed.

`LightKrylov` is thus closer to `Krylov.jl` [@montoison-2023] in Julia: a minimal package with a high level of abstraction specialised for Krylov methods only.
While the latter offers a broader collection of methods, we are actively working to bridge the gap.
Additionally, calling Julia code from Fortran remains a delicate process, burdened by the *two-language problem*.
In that regard, the fact that `LightKrylov` is written in pure, standard‑compliant Fortran makes it an ideal candidate for integration into existing Fortran codebases.
Moreover, its high level of abstraction enables users to re-use existing components of their codebase (including parallelisation), requiring only minimal changes while preserving computational performance.

# Performance in a production-ready open-source codebase

`LightKrylov` has been integrated into [`neklab`](https://github.com/nekStab/neklab), a toolbox for stability and bifurcation analysis for the spectral element solver `Nek5000`. The abstract vector interface allows direct use of `Nek5000`'s distributed data structures, and the pure-Fortran nature facilitated integration with its existing build system, demonstrating the library's suitability for large-scale HPC applications.

## Hydrodynamic stability of an unstable fixed point of the nonlinear Navier-Stokes equations

Using the two-dimensional flow past a circular cylinder at Reynolds number of 100, we showcase the efficient integration of `LightKrylov` and `Nek5000` and validate the results with algorithms provided by the [`KTH Framework`](https://github.com/KTH-Nek5000/KTH_Framework) toolbox [@kth-framework] based on the same solver.
Discretisation of the governing equations leads to systems with approximately 175,000 degrees of freedom.
All computations were run in parallel on 12 Intel Core Ultra 7 processors and the numerical settings are identical for both libraries.

The unstable fixed point of the nonlinear Navier-Stokes equations is computed using both `LightKrylov`'s *time-stepper*-based Newton-GMRES solver and the selected frequency damping implementation from `KTH Framework` .
Likewise, the leading eigenpair of the corresponding linearised Navier-Stokes operator is computed using `LightKrylov`'s implementation of the Krylov-Schur algorithm and the `KTH Framework`'s wrapper for `ARPACK` [@lehoucq:arpack:siam].

A visual comparison is provided in \autoref{fig:timings} showing excellent agreement.
The table in the lower-right panel of \autoref{fig:timings} summarizes the wall-clock times of the `neklab` computations. Isolating the intrinsic cost of the algorithms in `LightKrylov` from the cost of the calls to LAPACK and the linear and nonlinear Navier-Stokes solvers (`matvec` and `response`, respectively) shows that extended `abstract` types and object-oriented programming in Fortran incurs a negligible computational overhead for such large-scale applications.

![Validation of `LightKrylov`: Newton-GMRES and `eigs`. The top row depicts the streamwise velocity of the unstable solution computed using `LightKrylov` and the pointwise difference with the reference one computed with `KTH Framework`.](LK_newton_eigs_timings.pdf){#fig:timings}

# Perspectives

Despite being in its early development stage, `LightKrylov` has already been used in production runs on dozens of processors. It has also been interfaced with [`neko`](https://github.com/ExtremeFLOW/neko), a modernized implementation of `Nek5000` running on GPUs, and is currently being interfaced with [`dNami`](https://github.com/dNamiLab/dNami), a high-performance source code generator for hyperbolic partial differential equations.

Ongoing development efforts include:

- **Build system:** `LightKrylov` relies on the Fortran package manager `fpm` for its build. While it ensures nice integration into the modern Fortran ecosystem, it also limits its usage to codes built using `fpm`. A build system based on `CMake` for easier integration into non-`fpm` codes is in development.
- **Krylov processes:** `LightKrylov` currently provides implementations for three of the most important Krylov processes. Future releases will include the non-Hermitian Lanczos and two-sided Arnoldi factorizations, as well as the Saunders-Simon-Yip process [@ssy-krylov].
- **Iterative solvers:** The integration of new Krylov processes will enable us to provide an extended list of solvers, including `MINRES`, `CGNE`, `LSQR`, and `LSMR`.

### Acknowledgements

We acknowledge the financial support of the French National Agency for Research (ANR) through the ANR-33-CE46-0008-CONMAN grant agreement. We also thank the [fortran-lang](https://fortran-lang.org/) community for the development of [`stdlib`](https://stdlib.fortran-lang.org/), and in particular [Frederico Perini](https://github.com/perazz), [Jeremie Vandenplas](https://github.com/jvdp1), and [José Alvez](https://github.com/jalvesz) for their work on the `stdlib_linalg` module.

# References
