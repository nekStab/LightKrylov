module LightKrylov
   ! --> Utilities.
   use lightkrylov_Utils
   ! --> Definition of the abstract vector type.
   use lightkrylov_AbstractVector
   ! --> Definition of the abstract linear operator type.
   use lightkrylov_LinearOperator
   ! --> Implementation of the various Krylov decompositions.
   use lightkrylov_BaseKrylov
   use lightkrylov_RationalKrylov
   ! --> Iterative Solvers.
   use lightkrylov_IterativeSolvers
   implicit none
   include "dtypes.h"

   private

   !> Global variables.
   public :: greetings, wp, atol, rtol
   !> Abstract vectors.
   public :: abstract_vector, get_vec, mat_mult, mat_axpby, mat_zero, mat_copy
   !> Abstract linear operators.
   public :: abstract_linop, abstract_spd_linop, identity_linop, scaled_linop, axpby_linop
   !> Krylov factorization for general matrix.
   public :: initialize_krylov_subspace
   public :: arnoldi_factorization, nonsymmetric_lanczos_tridiagonalization, two_sided_arnoldi_factorization
   public :: lanczos_bidiagonalization, rational_arnoldi_factorization
   !> Krylov factorization for sym. pos. def. matrices.
   public :: lanczos_tridiagonalization
   !> Linear solvers.
   public :: gmres, cg
   public :: abstract_opts, gmres_opts, cg_opts
   public :: abstract_linear_solver, abstract_preconditioner
   !> Matrix factorization.
   public :: eigs, eighs, svds, two_sided_eigs
   public :: qr_factorization
   public :: save_eigenspectrum

contains

   subroutine greetings()
      write (*, *)
      write (*, *)
      write (*, *) "-----------------------------------------------------------------"
      write (*, *) "-----------------------------------------------------------------"
      write (*, *)

      write (*, *) "      _     _       _     _   _  __           _            "
      write (*, *) "     | |   (_) __ _| |__ | |_| |/ /_ __ _   _| | _____   __"
      write (*, *) "     | |   | |/ _` | '_ \| __| ' /| '__| | | | |/ _ \ \ / /"
      write (*, *) "     | |___| | (_| | | | | |_| . \| |  | |_| | | (_) \ V / "
      write (*, *) "     |_____|_|\__, |_| |_|\__|_|\_\_|   \__, |_|\___/ \_/  "
      write (*, *) "              |___/                     |___/              "

      write (*, *)
      write (*, *) "Developped by: Jean-Christophe Loiseau."
      write (*, *) "               Arts & Métiers Institute of Technology, 2023."
      write (*, *) "               jean-christophe.loiseau@ensam.eu"
      write (*, *)

      write (*, *) "Version -- 0.1.0"
      write (*, *) "License -- BSD 3-Clause"
      write (*, *)

      write (*, *) "-----------------------------------------------------------------"
      write (*, *) "-----------------------------------------------------------------"
      write (*, *)
      write (*, *)
   end subroutine greetings

end module LightKrylov
