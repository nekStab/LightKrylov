---
project: LightKrylov
summary: Lightweight implementation of Krylov methods using modern Fortran.
author: Jean-Christophe Loiseau, Simon Kern, and Ricardo Frantz
src_dir: src
include: src
exclude_dir: test
output_dir: API-doc
page_dir: doc
media_dir: doc/media
display: public
         protected
source: true
proc_internals: true
md_extensions: markdown.extensions.toc
sort: permission-alpha
project_github: https://github.com/nekStab/LightKrylov
github: https://github.com/nekStab/LightKrylov
project_website: https://nekstab.github.io/LightKrylov/
favicon: doc/media/favicon.ico
license: bsd
css: ./bootstrap.min.css
---

@warning
   This API documentation for the `LightKrylov` pacakge is a work in progress.
   It is build from the source code in the `main` branch and does not track the current development in `dev` or any other branches.
   If you use another branch, please refer to the in-code documentation.
   Use the navigation bar at the top of the screen to browse `modules`, `procedures`, `source files`, etc.
   The listings near the bootom of the page are incomplete.
@endwarning

This is the main API documentation landing page generated by [FORD](https://github.com/Fortran-FOSS-Programmers/ford#readme). This documentation is released under the [`CC-BY-SA`](https://creativecommons.org/licenses/by-sa/4.0/) license while the `LightKrylov` source code is distribution under the [`BSD-3 Clause`](https://opensource.org/license/bsd-3-clause) one.

## Scope

The goal of `LightKrylov` is to provide a lightweight implementation of many standard Krylov techniques using modern `Fortran` features, including: 

- Linear solvers for \( \mathbf{Ax} = \mathbf{b} \)
      + [`cg`](https://en.wikipedia.org/wiki/Conjugate_gradient_method) when \( \mathbf{A} \) is a symmetric (hermitian) positive definite matrix.
      + [`gmres`](https://en.wikipedia.org/wiki/Generalized_minimal_residual_method) when \( \mathbf{A} \) is a non-symmetric square linear operator.
- Krylov-based matrix factorizations
      + `arnoldi` - Construct an upper Hessenberg matrix \( \mathbf{H} \) and an orthonormal (unitary) basis \( \mathbf{X} \) capturing the dominant eigenspace of a square non-symmetric matrix \( \mathbf{A} \) using the [`Arnoldi`](https://en.wikipedia.org/wiki/Arnoldi_iteration) iterative process.
      + `lanczos` - Construct a symmetric (hermitian) tridiagonal matrix \( \mathbf{T} \) and an orthonormal (unitary) basis \( \mathbf{X} \) capturing the dominant eigenspace of a symmetric (hermitian) linear operator \( \mathbf{A} \) using the [`Lanczos`](https://en.wikipedia.org/wiki/Lanczos_algorithm) iterative process.
      + `bidiagonalization` - Construct a bidiagonal matrix \( \mathbf{B} \) and orthonormal bases \( \mathbf{U} \) and \( \mathbf{V} \) for the dominant column (resp. row) span of a general linear operator \( \mathbf{A} \) using the  [`Lanczos bidiagonalization`](https://en.wikipedia.org/wiki/Bidiagonalization) iterative process.
- Eigenvalue solvers for \( \mathbf{Ax} = \lambda \mathbf{x} \)
      + `eigs` to compute the largest eigenvalues and associated eigenvectors of a general square linear operator \( \mathbf{A} \) using the [`Arnoldi`](https://en.wikipedia.org/wiki/Arnoldi_iteration) iterative process.
      + `eighs` to compute the largest eigenvalues and associated eigenvectors of a symmetric (hermitian) linear operator \( \mathbf{A} \) using the [`Lanczos`](https://en.wikipedia.org/wiki/Lanczos_algorithm) iterative process.
- Singular value decomposition \( \mathbf{A} = \mathbf{U} \boldsymbol{\Sigma} \mathbf{V}^H \)
      + `svds` to compute the leading singular triplets of a general linear operator \( \mathbf{A} \) using the [`Lanczos bidiagonalization`](https://en.wikipedia.org/wiki/Bidiagonalization) iterative process.
- Solving a system of nonlinear equations \( F(\mathbf{x}) = \mathbf{0} \)
      + `newton` to find the solution to \( F(\mathbf{x}) = \mathbf{0} \) using a Newton-Krylov solver with optimal step size found by a simple bisection method.

While similar and more feature-complete packages exist (e.g. [Arpack](https://www.arpack.org), [SLEPC](https://slepc.upv.es/) or [Trilinos](https://trilinos.github.io/)), the use of `abstract_type` in `LightKrylov` and its nearly non-existant list of dependencies makes it far easier to incorporate into an existing code base. Preliminary benchmark results moreover show that it is on par with `Arpack` in terms of accuracy and computational performances.

## Acknowledgment

The development of `LightKrylov` is part of an on-going research project funded by [Agence Nationale pour la Recherche](https://anr.fr/en/) (ANR) under the grant agreement ANR-22-CE46-0008. The project started in January 2023 and will run until December 2026.
We are also very grateful to the [fortran-lang](https://fortran-lang.org/) community and the maintainers of [`stdlib`](https://github.com/fortran-lang/stdlib).

## Related projects

`LightKrylov` is the base package of our ecosystem. If you like it, you may also be interested in:

- [`LightROM`](https://github.com/nekStab/LightROM) : a lightweight Fortran package providing a set of functions for reduced-order modeling, control and estimation of large-scale linear time invariant dynamical systems.
- [`neklab`](https://github.com/nekStab/neklab) : a bifurcation and stability analysis toolbox based on `LightKrylov` for the massively parallel spectral element solver [`Nek5000`](https://github.com/Nek5000/Nek5000).
