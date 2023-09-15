# LightKrylov

[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

**LightKrylov** is a lightweight fortran implementation of Krylov subspace techniques. It serves as the base library for [**nekStab**](https://github.com/nekStab/nekStab), a toolbox for performing bifurcation analysis using the spectral element CFD solver [**Nek5000**](https://github.com/Nek5000/Nek5000).
Its primary use case if for bifurcation and stability analysis of large-scale nonlinear dynamical systems of the form

$$
\dot{\mathbf{X}} = \mathcal{F}(\mathbf{X}, \boldsymbol{\mu}),
$$

where $\mathbf{X} \in \mathbb{R}^{n}$ is the state vector of the system, $\boldsymbol{\mu} \in \mathbb{R}^p$ the parameter vector and $\mathcal{F} : \mathbb{R}^n \times \mathbb{R}^p \to \mathbb{R}^n$ the dynamics.
Given a fixed point $\mathbf{X}_*$ of the system, the dynamics of infinitesimal perturbations evolving in its vicinity are governed by the *linearized* equations

$$
\dot{\mathbf{x}} = \mathbf{A} \mathbf{x},
$$

where $\mathbf{A}$ is the *Jacobian* matrix of the system.
