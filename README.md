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

`LightKrylov` leverages Fortran's `abstract type` feature to provide generic implementations of the various Krylov methods.
The only requirement from the user to benefit from the capabilities of `LightKrylov` is to extend the `abstract_vector` and `abstract_linop` types to define their notion of vectors and linear operators. `LightKrylov` then provides the following functionalities:

- Krylov factorizations : `arnoldi_factorization`, `lanczos_tridiagonalization`, `lanczos_bidiagonalization`.
- Spectral analysis : `eigs`, `eighs`, `svds`.
- Linear systems : `gmres`, `cg`.

At the present time, none of the algorithms support Krylov-Schur restarting procedure although this is part of our plan for `LightKrylov v2.0`.

### Known limitations

For the sake of simplicity, `LightKrylov` only works with `real` or `double precision` data. While this might seem restrictive at first, consider that a complex-valued $n \times n$ linear system $\mathbf{Ax} = \mathbf{b}$ can always be rewritten using only real arithmetic as a $2n \times 2n$ real-valued system.
The primary reason to develop `LightKrylov` is to couple it with high-performance solvers in computational mechanics which often use exclusively real-valued data types. As such, we do not have any plan to natively support complex-valued linear operators or vectors.

### Examples

Several examples can be found in the `example` folder. These include:
- [Ginzburg-Landau]() : Serial computation of the leading eigenpairs of a complex-valued linear operator via time-stepping.
- [Laplace operator]() : Parallel computation of the leading eigenpairs of the Laplace operator defined on the unit-square.

Alternatively, you can also look at [`neklab`](), a bifurcation and stability analysis toolbox based on `LightKrylov` and designed to augment the functionalities of the massively parallel spectral element solver [`Nek5000`]().

| [**Ginzburg-Landau**]() | [**Laplace operator**]() |
| :---------------------: | :----------------------: |
|                         |                          |

## Installation

Provided you have `git` installed, getting the code is as simple as:

```
git clone https://github.com/nekStab/LightKrylov
```

Alternatively, using `gh-cli`, you can type

```
gh repo clone nekStab/LightKrylov
```

### Dependencies

`LightKrylov` has a very minimal set of dependencies. These only include:

- a Fortran compiler,
- [LAPACK]() (or similar),
- [`fpm`](https://github.com/fortran-lang/fpm) or `make` for building the code.

And that's all of it. To date, the tested compilers include:

- `gfortran 12.0.3`

### Building with `fpm`

Provided you have cloned the repo, installing `LightKrylov` with `fpm` is as simple as

```
fpm build --profile release
```

### Building with `make`

N/A

### Running the tests

To see if the library has been compiled correctly, a set of unit tests are provided in [test](). If you use `fpm`, running these tests is as simple as

```
fpm test
```

If everything went fine, you should see

```
All tests successfully passed!
```

If not, please feel free to open an Issue.

### Running the examples

To run the examples:

```
fpm run --example
```

This command will run all of the examples sequentially. You can alternatively run a specific example using e.g.

```
fpm run --example Ginzburg-Landau
```

For more details, please refer to each of the examples.

## Contributing

### Current developers

`LightKrylov` is currently developed and maintained by a team of three:
- [Jean-Christophe Loiseau](https://loiseaujc.github.io/) : Assistant Professor of Applied maths and Fluid dynamics at [DynFluid](https://dynfluid.ensam.eu/), Arts et Métiers Institute of Technology, Paris, France.
- [Ricardo Frantz](https://github.com/ricardofrantz) : PhD in Fluid dynamics (Arts et Métiers, France, 2022) and currently postdoctoral researcher at DynFluid.
- [Simon Kern](https://github.com/Simkern/) : PhD in Fluid dynamics (KTH, Sweden, 2023) and currently postdoctoral researcher at DynFluid.

Anyone else interested in contributing is obviously most welcomed!

## Acknowledgment

The development of `LightKrylov` is part of an on-going research project funded by [Agence Nationale pour la Recherche](https://anr.fr/en/) (ANR) under the grant agreement ANR-22-CE46-0008.

### Related projects

`LightKrylov` is the base package of our ecosystem. If you like it, you may also be interested in :
- [`LightROM`](https://github.com/nekStab/LightROM) : a lightweight Fortran package providing a set of functions for reduced-order modeling, control and estimation of large-scale linear time invariant dynamical systems.
- [`neklab`]() : a bifurcation and stability analysis toolbox based on `LightKrylov` for the massively parallel spectral element solver [`Nek5000`]().
