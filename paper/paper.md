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

[`LightKrylov`](https://github.com/nekStab/LightKrylov) is a Fortran package with minimal dependencies providing a collection of Krylov-based methods for solving linear systems, computing the leading eigenpairs of a linear operator as well as its leading singular triplets.

| Square systems   | Eigenvalue problems | Singular value decomposition           |
| :--------------: | :-----------------: | :------------------------------------: |
| $Ax = b$         | $Ax = \lambda x$    | $Av = \sigma u \quad A^H u = \sigma v$ |


When solving such problems, direct methods tend to have a computational cost scaling as $\mathcal{O}(n^3)$ (where $n$ is the leading dimension of $A$) making them impractical for very large-scale problems.
Even when $A$ is sparse, care needs to be taken when computing matrix factorizations such as LU or Cholesky since they can incur an $\mathcal{O}(n^2)$ storage cost in the worst case scenarios by introducing fill-in, where the factorization becomes significantly less sparse than the original matrix.
In contrast, Krylov methods only need a function computing the matrix-vector product $u \leftarrow Av$ (or $u \leftarrow A^H v$) to iteratively construct a *Krylov subspace* [@krylov-1931] given by $\mathcal{K}_k(A, v) = \text{span}\{v, Av, A^2v, \ldots, A^{k-1}v\}$ for an initial vector $v$.
We refer interested readers to @ipsen-1998 for an introduction to Krylov methods, to @saad-2003 for technical details and to @frantz-2023 for examples of their usage in the field of computational fluid dynamics.

# Statement of need

## A collection of Krylov-based algorithms in pure modern Fortran

`LightKrylov` aims to provide Fortran users with familiar `scipy`-inspired interfaces to a collection of widely used Krylov techniques.
These include:

- **Krylov processes:** Arnoldi factorization, Golub-Kahan bidiagonalization, as well as Lanczos tridiagonalization for Hermitian operators.
- **Krylov methods:**
    - **Linear systems -** Conjugate Gradient (CG), Generalized Minimal Residual method (GMRES) and Flexible GMRES.
    - **Spectral decomposition -** Arnoldi iteration with Krylov-Schur restart for non-Hermitian operators, Lanczos tridiagonalization with thick restart for Hermitian ones.
    - **SVD -** restarted Golub-Kahan bidiagonalization.

A block version of the Arnoldi iterative process is also provided.
Libraries exposing more methods exist in other languages, e.g. `Krylov.jl` [@montoison-2023] in `Julia` or `PETSc` [@petsc-web-page] for instance.
Integrating large multi-language libraries like PETSc into existing Fortran codes often presents significant challenges related to build systems, dependency management, and potential performance overhead from language interfacing ('two-language problem'). In contrast, `LightKrylov` provides a pure Fortran alternative, fully compliant with the Fortran 2018 standard, and only requires `stdlib` as external run-time dependency.
It leverages `stdlib` not only for basic utilities but specifically builds upon `stdlib_linalg` for foundational linear algebra operations, ensuring adherence to community standards and benefiting from its rapid ongoing developments.

Additionally, using `fypp` [@fypp-webpage] as build-time dependency, significantly reduces code duplication and maintenance effort, ensuring consistency across different data types (`single`, `double`, `extended` and `quadruple precision`) both for `real` and `complex` numbers while keeping the core algorithmic logic centralized. `fypp` is a Python-based Fortran preprocessor that automatically generates the necessary code for different data types.
Finally, its build process relies on the Fortran package manager `fpm`, greatly facilitating its installation and its incorporation into the modern Fortran ecosystem.

**Licensing :** `LightKrylov` is distributed under the permissive BSD-3 license, encouraging broad adoption and contribution.

**Testing/CI :** The library includes a suite of unit and integration tests, automatically executed via GitHub Actions for each new commit. The set of compilers currently tested in CI include `gfortran` (versions 12 to 14) as well as the Intel `ifort` and `ifx` compilers. The operating systems include the latest Ubuntu and Windows OS as well as MacOS 12.

## A focus on abstract linear operators and abstract vectors

From a mathematical point of view, Krylov methods can be implemented without making explicit reference to the particular data structure used to represent a vector or a linear operator nor to how the actual matrix-vector product is being implemented.
To do so, `LightKrylov` uses modern Fortran `abstract` type constructs.
A stripped-down version of the abstract vector type is shown below.

```fortran
type, abstract :: abstract_vector_rdp
  ! Abstract type defining a (double precision) vector which can be extended
  ! by users to accomodate their particular data structure and associated
  ! computational routines.
  contains
    ! Abstract procedure to compute the scalar-vector product.
    procedure(abstract_scal_rdp) , pass(self), deferred :: scal
    ! Abstract procedure to compute y = alpha*x + beta*y
    procedure(abstract_axpby_rdp), pass(self), deferred:: axpby
    ! Abstract procedure to compute the vector dot product.
    procedure(abstract_dot_rdp)  , pass(self), deferred :: dot
end type
```

These three type-bound procedures correspond to the basic set of operations on vectors, namely scalar-vector product, linear combination of two vectors, and the dot product. These operations correspond to the essential building blocks required by Krylov algorithms.
The signatures of these type-bound procedures follow, to the extent possible, the standard signatures of the corresponding `blas` functions for a more familiar use.
For instance, the `abstract_axpby_rdp` interface is defined as

```fortran
abstract interface
  subroutine abstract_axpby_rdp(alpha, vec, beta, self)
    double precision          , intent(in)    :: alpha, beta
    class(abstract_vector_rdp), intent(in)    :: vec
    class(abstract_vector_rdp), intent(inout) :: self
  end subroutine
end interface
```
mimicking the signature of the (extended) BLAS-2 subroutine `axpby`.
Note that, in practice, `LightKrylov` requires the definition of two additional procedures, one to set a vector to zero and the second to fill it with random data.

Similarly, a stripped-down version of an abstract linear operator type is shown below.

```fortran
type, abstract: abstract_linop_rdp
  ! Abstract type defining a (real, double precision) linear operator which
  ! can be extended by users to accomodate their particular data structure
  ! and associated kernels for the matrix-vector product.
  contains
    ! Abstract procedure to compute the matrix-vector product y = A * x
    procedure(abstract_matvec_rdp), pass(self), deferred :: matvec
    ! Abstract procedure to compute the matrix-vector product y = A^H * x
    procedure(abstract_matvec_rdp), pass(self), deferred :: rmatvec
end type
```

The two type-bound procedures need to implement the matrix-vector product and the transposed matrix-vector product, respectively.
The corresponding abstract interface reads

```fortran
abstract interface
  subroutine abstract_matvec_rdp(self, vec_in, vec_out)
    class(abstract_linop_rdp) , intent(inout) :: self
    class(abstract_vector_rdp), intent(in)    :: vec_in
    class(abstract_vector_rdp), intent(out)   :: vec_out
  end subroutine
end interface
```
Using such `abstract` types enables us to focus on the high-level implementation of the different Krylov-based algorithms while leaving the performance-critical details of the different vector and matrix-vector operations to the users.
After having extended these `abstract` types for their particular applications, users can solve linear systems or compute eigenvalues as easily as

```fortran

! Solve linear system.
call gmres(A, b, x, info)

! Compute eigenvalues and eigenvectors.
call eigs(A, V, lambda, residuals, info)
```

In addition, `LightKrylov` exposes abstract interfaces to define preconditioners, which can improve significantly the convergence of iterative solvers, as well as a Newton-GMRES solver to find the roots of a multivariate function $f : \mathbb{K}^n \to \mathbb{K}^n$ (where $\mathbb{K}$ is a placeholder for $\mathbb{R}$ or $\mathbb{C}$).
Examples illustrating how to extend these `abstract` types and interfaces for computing the leading eigenpairs of the linearized Ginzburg-Landau equation or to compute unstable periodic orbits of the Rössler system and study its Lyapunov exponents can be found [here](https://github.com/nekStab/LightKrylov/tree/main/example).

# Performances in a production-ready open-source code

`LightKrylov` was successfully integrated into [`neklab`](https://github.com/nekStab/neklab), a toolbox for stability and bifurcation analysis using the high-performance spectral element solver `Nek5000`. The abstract vector interface allowed direct use of `Nek5000`'s distributed data structures, and the pure-Fortran nature facilitated integration with its existing build system, demonstrating the library's suitability for large-scale HPC applications.

## Computing a fixed point of the nonlinear Navier-Stokes equations

Computing (stable or unstable) fixed points of the nonlinear governing equations is a often a necessary step when studying the stability properties of a nonlinear system.
In hydrodynamic instability applications, this amount to finding a stationary solution of the nonlinear Navier-Stokes equations.
Using the two-dimensional flow past a circular cylinder as an example, the figure below depicts the convergence history for two fixed point computation techniques, namely `LightKrylov`'s Newton-GMRES solver and the selective frequency damping algorithm [@pof-sfd] (often considered as the gold standard in this community).

**ADD FIGURE**

Using the spectral element solver `Nek5000`, the governing equations are discretized with 1996 elements using a 7-th order polynomial approximation in each direction resulting in 127 744 grid points and roughly 350 000 degrees of freedom (i.e. the two velocity components as well as the pressure).
Both solvers make use of the same initial condition and computations were run on some numbers of Intel something-something processors.
Note that the Newton-GMRES solver make use of the *time-stepper* approach recently reviewed in @frantz-2023 while the selective frequency damping solver uses the implementation provided by the [`KTH Framework`](https://github.com/KTH-Nek5000/KTH_Framework) toolbox [@kth-framework].
Leveraging the abstract interfaces provided by `LightKrylov`, the implementation of the Newton-GMRES solver in `neklab` required a very minimal set of modification of the `Nek5000` solver. The code to run this example can found online [here](???).

## Computing the leading eigenpair of the linearized Navier-Stokes operator

Investigating the spectral properties of the linearized Navier-Stokes operator is critical to better understand the transition to unsteadiness and turbulence in many flow applications. Using once again the two-dimensional cylinder flow as an example, the table below compares the wall-clock time and number of matrix-vector products needed to converge the first pair of complex-conjugate eigenvalues. The two solvers considered include once again `neklab` and `KTH Framework`. While `neklab` uses `LightKrylov`'s abstract interfaces and implementation of the Krylov-Schur algorithm, `KTH Framework` relies under the hood on `ARPACK` and its reversed communication methodology.

**ADD TABLE**

Both solvers effectively use a *time-stepper* approach, approximating the leading eigenvalues of the exponential propagator $\exp(\tau \mathbf{A})$ rather than $\mathbf{A}$ directly.
For both cases, the sampling time is set to $\tau = ??$ and the maximum dimension of the Krylov subspace is set to $\mathrm{dim}(\mathcal{K}) = 64$. All other parameters, including those internal to `Nek5000`, are set to the same values for both solvers. The code to run this example can be found [here](?).
From table ?, it can be seen that `LightKrylov`'s eigenvalue solver is competitive with `ARPACK` and that the use of extended `abstract` types and object-oriented programming Fortran lead to no computational performances degradation for such large-scale applications.

<!-- Figure sizes can be customized by adding an optional second parameter: -->
<!-- ![Caption for example figure.](figure.png){ width=20% } -->

# Perspectives

Despite being in its early development stage, `LightKrylov` has already been used in hydrodynamic stability production runs on ??? processors. It has also been successfully interfaced with [`neko`](https://github.com/ExtremeFLOW/neko), a modernized implementation of `Nek5000` running on GPU and is currently in the process of being interfaced with [`dNami`](https://github.com/dNamiLab/dNami), a high-performance source code generator for hyperbolic partial differential equations.

When it comes to the development of the package itself, on-going efforts include:

- **Build system :** Currently, `LightKrylov` relies exclusively on the Fortran package manager `fpm` for its build. While it ensures nice integration into the modern Fortran ecosystem, it also limits its usage to codes build using `fpm` as well. Following what is being done with the Fortran standard library `stdlib`, plans are being made to propose an alternative build system based on `Cmake` for easier integration into non-`fpm` codes.
- **Krylov processes :** As of summer 2025, `LightKrylov` provides implementations for three of the most important Krylov processes. Future releases will extend this list by including the non-Hermitian Lanczos and two-sided Arnoldi factorizations useful for reduced-order modeling of LTI systems based on moment-matching as well as the Saunders-Simon-Yip process [@ssy-krylov].
- **Iterative solvers :** Beyond `CG` and `GMRES`, the integration of new Krylov processes will enable us to provide implementations for an extended list of iterative linear solvers including `MINRES` (indefinite Hermitian system), `CGNE`, `LSQR` (least-squares problems) and `LSMR` (least-norm problems) to begin with.

It needs to be emphasized finally that the development of `LightKrylov` sprang a community-driven effort in the fortran-lang ecosystem to integrate some of its features (namely the abstraction layer and the iterative solvers) into the standard library `stdlib`.

### Acknowledgements

We acknowledge the financial support of the French National Agency for Research (ANR) through the ANR-33-CE46-0008-CONMAN grant agreement. We also would like to thank the [fortran-lang](https://fortran-lang.org/) community for the development of [`stdlib`](https://stdlib.fortran-lang.org/), and in particular Frederico Perini, Jeremie Vandenplas, and Jose Alvez for their work on the `stdlib_linalg` module on top which `LightKrylov` heavily relies on.

# References
