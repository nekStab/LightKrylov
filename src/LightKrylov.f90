module LightKrylov
  ! --> Utilities.
  use Utils
  ! --> Definition of the abstract vector type.
  use AbstractVector
  ! --> Definition of the abstract linear operator type.
  use LinearOperator
  ! --> Implementation of the various Krylov decompositions.
  use BaseKrylov
  use RationalKrylov
  ! --> Iterative Solvers.
  use IterativeSolvers
  implicit none
  include "dtypes.h"

  private

  public :: greetings, wp, atol, rtol,                                     &
       abstract_vector, get_vec,                                           &
       abstract_linop, abstract_spd_linop, identity_linop, scaled_linop, axpby_linop, &
       arnoldi_factorization, lanczos_tridiagonalization, lanczos_bidiagonalization, &
       nonsymmetric_lanczos_tridiagonalization, rational_arnoldi_factorization, &
       eigs, eighs, gmres, save_eigenspectrum, svds, cg, bicgstab, &
       abstract_opts, gmres_opts, bicgstab_opts, cg_opts

contains

  subroutine greetings()
    write(*, *)
    write(*, *)
    write(*, *) "-----------------------------------------------------------------"
    write(*, *) "-----------------------------------------------------------------"
    write(*, *)

    write(*, *) "      _     _       _     _   _  __           _            "
    write(*, *) "     | |   (_) __ _| |__ | |_| |/ /_ __ _   _| | _____   __"
    write(*, *) "     | |   | |/ _` | '_ \| __| ' /| '__| | | | |/ _ \ \ / /"
    write(*, *) "     | |___| | (_| | | | | |_| . \| |  | |_| | | (_) \ V / "
    write(*, *) "     |_____|_|\__, |_| |_|\__|_|\_\_|   \__, |_|\___/ \_/  "
    write(*, *) "              |___/                     |___/              "


    write(*, *)
    write(*, *) "Developped by: Jean-Christophe Loiseau."
    write(*, *) "               Arts & Métiers Institute of Technology, 2023."
    write(*, *) "               jean-christophe.loiseau@ensam.eu"
    write(*, *)

    write(*, *) "Version -- 0.1.0"
    write(*, *) "License -- BSD 3-Clause"
    write(*, *)

    write(*, *) "-----------------------------------------------------------------"
    write(*, *) "-----------------------------------------------------------------"
    write(*, *)
    write(*, *)
  end subroutine greetings

end module LightKrylov
