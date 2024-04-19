module LightKrylov
    ! --> Definitions of various constants.
    use lightkrylov_constants
    ! --> Set of utility functions.
    use lightkrylov_utils
    ! --> Definitions of the abstract vector types.
    use lightkrylov_AbstractVectors
    ! --> Definitions of the abstract linear operators.
    use lightkrylov_AbstractLinops
    ! --> Standard Krylov techniques.
    use lightkrylov_BaseKrylov
    ! --> Iterative solvers.
    use lightkrylov_IterativeSolvers
    implicit none
    private

    ! LightKrylov exports.
    public :: greetings

    ! Constants exports.
    public :: sp, atol_sp, rtol_sp
    public :: dp, atol_dp, rtol_dp

    ! AbstractVectors exports.
    public :: abstract_vector
    public :: abstract_vector_rsp
    public :: abstract_vector_rdp
    public :: abstract_vector_csp
    public :: abstract_vector_cdp
    
    ! AbstractLinops exports.
    public :: abstract_linop
    public :: abstract_linop_rsp
    public :: abstract_linop_rdp
    public :: abstract_linop_csp
    public :: abstract_linop_cdp
    
    ! BaseKrylov exports.
    public :: qr
    public :: apply_permutation_matrix, apply_inverse_permutation_matrix
    public :: arnoldi
    public :: initialize_krylov_subspace
    public :: lanczos_bidiagonalization

    ! IterativeSolvers exports.
    public :: eigs
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
