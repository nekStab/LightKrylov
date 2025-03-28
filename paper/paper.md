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

[`LightKrylov`](https://github.com/nekStab/LightKrylov) is a Fortran package providing a collection of Krylov-based methods for solving linear systems, computing the leading eigenpairs of a linear operator as well as its leading singular triplets.

| Square systems   | Eigenvalue problems | Singular value decomposition           |
| :--------------: | :-----------------: | :------------------------------------: |
| $Ax = b$         | $Ax = \lambda x$    | $Av = \sigma u \quad A^H u = \sigma v$ |

When solving problems such as the ones shown above, direct methods tend to have a computational cost scaling as $\mathcal{O}(n^3)$ (where $n$ is the leading dimension of $A$) making them impractical for very large-scale problems.
Even when $A$ is sparse, care needs to be taken when computing matrix factorizations such as LU or Cholesky since they can incur an $\mathcal{O}(n^2)$ storage cost in the worst case scenarios.
In contrast, Krylov methods only need a function computing the matrix-vector product $u \leftarrow Av$ (and possibly $u \leftarrow A^H v$) to iteratively construct a *Krylov subspace* [@krylov-1931] from which an approximate solution can be extracted.
We refer interested readers to @ipsen-1998 for an introduction to Krylov methods, to @saad-2003 for technical details and to @frantz-2023 for examples of their usage in the field of computational fluid dynamics.


# Statement of need

## A collection of Krylov-based algorithms

`LightKrylov` aims to provide Fortran users with familiar `scipy`-inspired interfaces to a collection of widely used Krylov techniques.
These include:

- **Krylov processes:** Arnoldi factorization, Golub-Kahan bidiagonalization, as well as Lanczos tridiagonalization for Hermitian operators.
- **Krylov methods:**
    - **Linear systems -** Conjugate Gradient (CG), Generalized Minimal Residual method (GMRES) and Flexible GMRES.
    - **Spectral decomposition -** Arnoldi iteration with Krylov-Schur restart for non-Hermitian operators, Lanczos tridiagonalization with thick restart for Hermitian ones.
    - **SVD -** restarted Golub-Kahan bidiagonalization.

A block version of the Arnoldi iterative process is also provided.
Libraries exposing more methods exist in other languages, e.g. `Krylov.jl` [@montoison-2023] in `Julia` or `PETSc` [@petsc-web-page] for instance.
Integrating these into an existing Fortran code base might however prove challenging, either because of the two languages problem or the need of special privileges for installation on a cluster and managing the dependencies.
In contrast, `LightKrylov` is written in pure Fortran, fully compliant with the Fortran 2018 standard, and only requires `stdlib` as external dependency.
It is being developed with `fypp` [@fypp-webpage], a Python-based Fortran preprocesor enabling the automatic generation of `single`, `double`, `extended` and `quadruple precision` versions of the methods, both for `real` and `complex` numbers.
Finally, its build process relies on the Fortran package manager `fpm`, greatly facilitating its installation and its incorporation into the modern Fortran ecosystem.

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

The three abstract type-bound procedures correspond to the basic set of operations on vectors, namely scalar-vector product, linear combination of two vectors, and the dot product.
The signatures of these type-bound procedures follow, to the extent possible, the standard signatures of the corresponding `blas` functions.
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

# Examples

<!-- Although there exist libraries exposing more methods in other languages, e.g. `Krylov.jl` [@montoison-2023] in `Julia`, it needs to be emphasized that `LightKrylov` is written in pure Fortran and relies on a minimalistic set of dependencies which can all be taken care of using the Fortran package manager `fpm` [@?]. -->
<!-- Hence, it makes it a suitable choice for integration into existing code bases used in numerous areas of high-performance scientific computing. -->
<!-- As an example, `LightKrylov` has recently been integrated into `neklab`, an open-source toolbox for stability and bifurcation analysis built around the massively parallel spectral element solver `Nek5000`. -->
<!---->
<!-- ## Support for `abstract` vectors and linear operators -->
<!---->
<!-- # Example -->

<!-- # Citations -->
<!---->
<!-- Citations to entries in paper.bib should be in -->
<!-- [rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html) -->
<!-- format. -->
<!---->
<!-- If you want to cite a software repository URL (e.g. something on GitHub without a preferred -->
<!-- citation) then you can do it with the example BibTeX entry below for @fidgit. -->
<!---->
<!-- For a quick reference, the following citation commands can be used: -->
<!-- - `@author:2001`  ->  "Author et al. (2001)" -->
<!-- - `[@author:2001]` -> "(Author et al., 2001)" -->
<!-- - `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)" -->
<!---->
<!-- # Figures -->
<!---->
<!-- Figures can be included like this: -->
<!-- ![Caption for example figure.\label{fig:example}](figure.png) -->
<!-- and referenced from text using \autoref{fig:example}. -->
<!---->
<!-- Figure sizes can be customized by adding an optional second parameter: -->
<!-- ![Caption for example figure.](figure.png){ width=20% } -->

# Acknowledgements

We acknowledge the financial support of the French National Agency for Research (ANR) through the ANR-33-CE46-0008-CONMAN grant agreement. We also would like to thank the [fortran-lang](https://fortran-lang.org/) community for the development of [`stdlib`](https://stdlib.fortran-lang.org/), and in particular Frederico Perini, Jeremie Vandenplas, and Jose ??? for their work on the `stdlib_linalg` module on top which `LightKrylov` heavily relies.

# References
