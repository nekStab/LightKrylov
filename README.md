<img src="imgs/logo-white.png" style="align:center; width:512px" />



|                         **License**                          |                       **Build Status**                       | **Documentation** |
| :----------------------------------------------------------: | :----------------------------------------------------------: | :---------------: |
| [![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause) | [![Github actions](https://github.com/nekStab/LightKrylov/actions/workflows/gcc.yml/badge.svg?event=push)](https://github.com/nekStab/LightKrylov/actions) |                   |

# LightKrylov

Targeting large-scale linear algebra applications where the matrix $\mathbf{A}$ is only defined implicitly (e.g. through a call to a `matvec` subroutine), this package provides lightweight Fortran implementations of certain of the most useful Krylov methods to solve a variety of problems, among which:

1. Eigenvalue Decomposition
   $$\mathbf{A} \mathbf{x} = \lambda \mathbf{x}$$

2. Singular Value Decomposition
   $$\mathbf{A} = \mathbf{U} \boldsymbol{\Sigma} \mathbf{V}^T$$


3. Linear system of equations
   $$\mathbf{Ax} = \mathbf{b}$$

Krylov methods are particularly appropriate in situations where such problems must be solved but factorizing the matrix $\mathbf{A}$ is not possible because:

- $\mathbf{A}$ is not available explicitly but only implicitly through a `matvec` subroutine computing the matrix-vector product $\mathbf{Ax}$.
- $\mathbf{A}$ or its factors (e.g. `LU` or `Cholesky`) are dense and would consume an excessive amount of memory.

Krylov methods are *iterative methods*, i.e. they iteratively refine the solution of the problem until a desired accuracy is reached. While they are not recommended when a machine-precision solution is needed, they can nonetheless provide highly accurate approximations of the solution after a relatively small number of iterations. Krylov methods form the workhorses of large-scale numerical linear algebra.

## Capabilities

## Installation

## Help and support

## Contributing

### Current developers

### Bug reports and contributions

## References

### How to cite

### Related projects
