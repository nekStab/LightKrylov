---
title: '`LightKrylov`: Lightweight implementation of Krylov subspace techniques in modern Fortran'
tags:
  - Fortran
  - Numerical linear algebra
  - Krylov methods
  - Sparse linear systems
  - Eigenvalues and singular values
authors:
  - name: Jean-Christophe Loiseau
    orcid: 0000-0002-7244-8416
    corresponding: true
    equal-contrib: true
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: J. Simon Kern
    orcid: 0000-0002-2460-578X
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 1
  - name: Ricardo S. Frantz
    orcid: 0000-0001-8219-3368
    equal-contrib: true
    affiliation: 1
affiliations:
 - name: Arts et Métiers Institute of Technology
   index: 1
date: 13 August 2017
bibliography: paper.bib
---

# Summary

[`LightKrylov`](https://github.com/nekStab/LightKrylov) is a modern Fortran package implementing a collection of Krylov subspace methods for solving a variety of linear problems:

| Square systems   | Eigenvalue problems | Singular value decomposition           |
| :--------------: | :-----------------: | :------------------------------------: |
| $Ax = b$         | $Ax = \lambda x$    | $Av = \sigma u \quad A^H u = \sigma v$ |

Here, $A^H$ denotes the conjugate transpose of $A$, coinciding with the standard transpose if $A$ is real-valued.
Direct methods often require explicit storage of $A$.
Even if $A$ is sparse, computing its QR or Cholesky factorization to invert a linear system might lead to prohibitve storage requirements or computational costs.
In contrast, Krylov methods only need a function computing the matrix-vector product $u \leftarrow Av$ (and possibly $u \leftarrow A^H v$) to iteratively construct the *Krylov subspace* [@krylov-1931].
Over the past decades, Krylov methods have become a critical element of high-performance computing.
We refer interested readers to [@ipsen-1998] for an introduction to Krylov methods, to [@saad-2003] for technical details and to [@frantz-2023] for examples of their usage in the field of computational fluid dynamics.

# Statement of need

## A collection of Krylov-based algorithms

## Support for `real` and `complex` types for standard floating-point systems supported by Fortran

## Support for `abstract` vectors and linear operators

# Example

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

We acknowledge the financial support of the French National Agency for Research (ANR) through the ANR-????-CONMAN grant agreement. We also would like to thank the [fortran-lang](https://fortran-lang.org/) community for the development of [`stdlib`](https://stdlib.fortran-lang.org/), and in particular Frederico Perini, Jeremie Vandenplas, and Jose ??? for their work on the `stdlib_linalg` module on top which `LightKrylov` heavily relies.

# References
