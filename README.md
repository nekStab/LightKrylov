# LightKrylov

[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![Github actions](https://github.com/nekStab/LightKrylov/actions/workflows/gcc.yml/badge.svg?event=push)](https://github.com/nekStab/LightKrylov/actions)

**LightKrylov** is a lightweight fortran implementation of Krylov subspace techniques. It serves as the base library for [**nekStab**](https://github.com/nekStab/nekStab), a toolbox for performing bifurcation analysis using the spectral element CFD solver [**Nek5000**](https://github.com/Nek5000/Nek5000).
Its primary use case is for bifurcation and stability analysis of large-scale nonlinear dynamical systems of the form

$$
\dot{\mathbf{X}} = \mathcal{F}(\mathbf{X}, \boldsymbol{\mu}),
$$

where $\mathbf{X} \in \mathbb{R}^{n}$ is the state vector of the system, $\boldsymbol{\mu} \in \mathbb{R}^p$ the parameter vector and $\mathcal{F} : \mathbb{R}^n \times \mathbb{R}^p \to \mathbb{R}^n$ the dynamics.
Given a fixed point $\mathbf{X}_*$ of the system, the dynamics of infinitesimal perturbations evolving in its vicinity are governed by the *linearized* equations

$$
\dot{\mathbf{x}} = \mathbf{A} \mathbf{x},
$$

where $\mathbf{A}$ is the *Jacobian* matrix of the system.
Its solution is given by

$$
\mathbf{x}(\tau) = \exp\left( \tau \mathbf{A} \right) \mathbf{x}_0.
$$

**LightKrylov** aims at providing a simple set of Krylov-based techniques to study the spectral properties of the *exponential propagator* $\exp \left( \tau \mathbf{A} \right)$.
This operator is represented through an *abstract type* and its implementation is deferred to the user, typically through the use of a *matrix-free approach* (e.g. *time-stepping*).
More details about the *time-stepper* approach to stability analysis of large-scale dynamical systems can be found in [1].

### References

[1] R. S. Frantz, J.-Ch. Loiseau, and J.-Ch. Robinet. Krylov methods for large-scale dynamical systems: applications in fluid dynamics. *Appl. Mech. Rev.*, 2023. [[arXiv]](https://arxiv.org/abs/2301.12940)
