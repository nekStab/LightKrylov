module LightKrylov
    ! --> Definitions of various constants.
    use LightKrylov_Constants
    ! --> Set of utility functions.
    use LightKrylov_Utils
    ! --> Definitions of the abstract vector types.
    use LightKrylov_AbstractVectors
    ! --> Definitions of the abstract linear operators.
    use LightKrylov_AbstractLinops
    ! --> Definitions of the abstract dynamical systems.
    use LightKrylov_AbstractSystems
    ! --> Standard Krylov techniques.
    use LightKrylov_BaseKrylov
    ! --> Iterative solvers.
    use LightKrylov_NewtonKrylov
    ! --> Newton-Krylov fixed-point iteration.
    use LightKrylov_IterativeSolvers
    ! --> Expmlib
    use LightKrylov_Expmlib
    ! --> TestTypes
    implicit none
    private

    ! LightKrylov exports.
    public :: greetings

    ! Constants exports.
    public :: sp, atol_sp, rtol_sp
    public :: dp, atol_dp, rtol_dp

    ! Utils exports.
    public :: gmres_sp_opts
    public :: cg_sp_opts
    public :: newton_sp_opts
    public :: gmres_dp_opts
    public :: cg_dp_opts
    public :: newton_dp_opts

    ! AbstractVectors exports.
    public :: abstract_vector
    public :: abstract_vector_rsp
    public :: abstract_vector_rdp
    public :: abstract_vector_csp
    public :: abstract_vector_cdp
    public :: innerprod
    public :: linear_combination
    public :: axpby_basis
    public :: zero_basis
    public :: copy
    public :: rand_basis
    
    ! AbstractLinops exports.
    public :: abstract_linop
    public :: abstract_linop_rsp
    public :: adjoint_linop_rsp
    public :: Id_rsp
    public :: scaled_linop_rsp
    public :: axpby_linop_rsp
    public :: abstract_sym_linop_rsp
    public :: abstract_linop_rdp
    public :: adjoint_linop_rdp
    public :: Id_rdp
    public :: scaled_linop_rdp
    public :: axpby_linop_rdp
    public :: abstract_sym_linop_rdp
    public :: abstract_linop_csp
    public :: adjoint_linop_csp
    public :: Id_csp
    public :: scaled_linop_csp
    public :: axpby_linop_csp
    public :: abstract_hermitian_linop_csp
    public :: abstract_linop_cdp
    public :: adjoint_linop_cdp
    public :: Id_cdp
    public :: scaled_linop_cdp
    public :: axpby_linop_cdp
    public :: abstract_hermitian_linop_cdp

    ! AbstractSystems exports.
    public :: abstract_system_rsp
    public :: abstract_jacobian_linop_rsp
    public :: abstract_system_rdp
    public :: abstract_jacobian_linop_rdp
    public :: abstract_system_csp
    public :: abstract_jacobian_linop_csp
    public :: abstract_system_cdp
    public :: abstract_jacobian_linop_cdp
    
    ! BaseKrylov exports.
    public :: qr
    public :: apply_permutation_matrix, apply_inverse_permutation_matrix
    public :: arnoldi
    public :: initialize_krylov_subspace
    public :: orthogonalize_against_basis
    public :: orthonormalize_basis
    public :: bidiagonalization
    public :: lanczos
    public :: krylov_schur

    ! IterativeSolvers exports.
    public :: abstract_precond_rsp
    public :: abstract_precond_rdp
    public :: abstract_precond_csp
    public :: abstract_precond_cdp
    public :: abstract_linear_solver_rsp
    public :: abstract_linear_solver_rdp
    public :: abstract_linear_solver_csp
    public :: abstract_linear_solver_cdp
    public :: eigs, eighs, save_eigenspectrum
    public :: svds
    public :: gmres
    public :: gmres_rsp
    public :: gmres_rdp
    public :: gmres_csp
    public :: gmres_cdp
    public :: cg

    ! Newton-Krylov exports.
    public :: newton
    public :: constant_atol_sp, dynamic_tol_sp
    public :: constant_atol_dp, dynamic_tol_dp

    ! ExpmLib exports.
    public :: abstract_exptA_rsp
    public :: abstract_exptA_rdp
    public :: abstract_exptA_csp
    public :: abstract_exptA_cdp
    public :: expm
    public :: kexpm
    public :: k_exptA

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
      write (*, *) "Developed by: Jean-Christophe Loiseau"
      write (*, *) "              J. Simon Kern"
      write (*, *) "              Arts & Métiers Institute of Technology, 2023."
      write (*, *) "              jean-christophe.loiseau@ensam.eu"
      write (*, *)

      write (*, *) "Version -- beta 0.1.0"
      write (*, *) "License -- BSD 3-Clause"
      write (*, *)

      write (*, *) "-----------------------------------------------------------------"
      write (*, *) "-----------------------------------------------------------------"
      write (*, *)
      write (*, *)
   end subroutine greetings

end module LightKrylov
