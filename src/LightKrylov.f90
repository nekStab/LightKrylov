module LightKrylov
  ! --> Definition of the abstract vector type.
  use AbstractVector
  ! --> Definition of the abstract linear operator type.
  use LinearOperator
  ! --> Implementation of the various Krylov decompositions.
  use KrylovDecomp
  ! --> Iterative Solvers.
  use IterativeSolvers
  implicit none
  include "dtypes.h"

  private

  public :: greetings, wp,                                                 &
       abstract_vector, get_vec,                                           &
       abstract_linop, abstract_spd_linop,                                 &
       arnoldi_factorization, lanczos_tridiagonalization, lanczos_bidiagonalization, &
       eigs, eighs, gmres, save_eigenspectrum, svds

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
