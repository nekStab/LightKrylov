module lightkrylov_IterativeSolvers
    !!  This module provides some of the most important computational routines provided by
    !!  `LightKrylov`. These include:
    !!
    !!  - `eigs` : Compute the leading eigenpairs of a square linear operator \( A \).
    !!  - `eighs` : Compute the leading eigenpairs of a symmetric positive definite 
    !!              operator \( A \).
    !!  - `svds` : Compute the leading singular triplets of a linear operator \( A \).
    !!  - `gmres` : Solve the linear system \( Ax = b \) using the *generalized minimum
    !!              residual method*.
    !!  - `cg` : Solve the linear system \( Ax = b \) where \( A \) is symmetric
    !!           positive definite using the *Conjugate Gradient* method.
    !!
    !!  It also provides abstract interfaces to pass user-defined solvers and preconditioners
    !!  to `LightKrylov`. Note that these features are still experimental however.

    !--------------------------------------------
    !-----     Fortran Standard Library     -----
    !--------------------------------------------
    use iso_fortran_env, only: output_unit
    
    use stdlib_sorting, only: sort_index
    use stdlib_optval, only: optval
    use stdlib_io_npy, only: save_npy
    use stdlib_linalg, only: lstsq, svd, eigh
    use stdlib_stats, only: median

    !-------------------------------
    !-----     LightKrylov     -----
    !-------------------------------
    use LightKrylov_Constants
    Use LightKrylov_Logger
    use LightKrylov_Utils
    use LightKrylov_Timing, only: timer => global_lightkrylov_timer, time_lightkrylov
    use LightKrylov_AbstractVectors
    use LightKrylov_AbstractLinops
    use LightKrylov_BaseKrylov

    implicit none
    private

    character(len=*), parameter :: this_module      = 'LK_Solvers'
    character(len=*), parameter :: this_module_long = 'LightKrylov_IterativeSolvers'

    public :: abstract_linear_solver_rsp
    public :: abstract_linear_solver_rdp
    public :: abstract_linear_solver_csp
    public :: abstract_linear_solver_cdp
    public :: save_eigenspectrum
    public :: eigs
    public :: eighs
    public :: svds
    public :: gmres
    public :: gmres_rsp
    public :: gmres_rdp
    public :: gmres_csp
    public :: gmres_cdp
    public :: fgmres
    public :: fgmres_rsp
    public :: fgmres_rdp
    public :: fgmres_csp
    public :: fgmres_cdp
    public :: cg
    public :: cg_rsp
    public :: cg_rdp
    public :: cg_csp
    public :: cg_cdp

    !------------------------------------------------------
    !-----     ABSTRACT PRECONDITIONER DEFINITION     -----
    !------------------------------------------------------

    type, abstract, public :: abstract_precond_rsp
    contains
        private
        procedure(abstract_apply_rsp), pass(self), public, deferred :: apply
    end type

    abstract interface
        subroutine abstract_apply_rsp(self, vec, iter, current_residual, target_residual)
            !! Abstract interface to apply a preconditioner in `LightKrylov`.
            import abstract_precond_rsp, abstract_vector_rsp
            import sp
            class(abstract_precond_rsp), intent(inout) :: self
            !! Preconditioner.
            class(abstract_vector_rsp), intent(inout) :: vec
            !! Input/Output vector.
            integer, optional, intent(in) :: iter
            !! Current iteration number.
            real(sp), optional, intent(in) :: current_residual
            real(sp), optional, intent(in) :: target_residual
        end subroutine abstract_apply_rsp
    end interface

    type, abstract, public :: abstract_precond_rdp
    contains
        private
        procedure(abstract_apply_rdp), pass(self), public, deferred :: apply
    end type

    abstract interface
        subroutine abstract_apply_rdp(self, vec, iter, current_residual, target_residual)
            !! Abstract interface to apply a preconditioner in `LightKrylov`.
            import abstract_precond_rdp, abstract_vector_rdp
            import dp
            class(abstract_precond_rdp), intent(inout) :: self
            !! Preconditioner.
            class(abstract_vector_rdp), intent(inout) :: vec
            !! Input/Output vector.
            integer, optional, intent(in) :: iter
            !! Current iteration number.
            real(dp), optional, intent(in) :: current_residual
            real(dp), optional, intent(in) :: target_residual
        end subroutine abstract_apply_rdp
    end interface

    type, abstract, public :: abstract_precond_csp
    contains
        private
        procedure(abstract_apply_csp), pass(self), public, deferred :: apply
    end type

    abstract interface
        subroutine abstract_apply_csp(self, vec, iter, current_residual, target_residual)
            !! Abstract interface to apply a preconditioner in `LightKrylov`.
            import abstract_precond_csp, abstract_vector_csp
            import sp
            class(abstract_precond_csp), intent(inout) :: self
            !! Preconditioner.
            class(abstract_vector_csp), intent(inout) :: vec
            !! Input/Output vector.
            integer, optional, intent(in) :: iter
            !! Current iteration number.
            real(sp), optional, intent(in) :: current_residual
            real(sp), optional, intent(in) :: target_residual
        end subroutine abstract_apply_csp
    end interface

    type, abstract, public :: abstract_precond_cdp
    contains
        private
        procedure(abstract_apply_cdp), pass(self), public, deferred :: apply
    end type

    abstract interface
        subroutine abstract_apply_cdp(self, vec, iter, current_residual, target_residual)
            !! Abstract interface to apply a preconditioner in `LightKrylov`.
            import abstract_precond_cdp, abstract_vector_cdp
            import dp
            class(abstract_precond_cdp), intent(inout) :: self
            !! Preconditioner.
            class(abstract_vector_cdp), intent(inout) :: vec
            !! Input/Output vector.
            integer, optional, intent(in) :: iter
            !! Current iteration number.
            real(dp), optional, intent(in) :: current_residual
            real(dp), optional, intent(in) :: target_residual
        end subroutine abstract_apply_cdp
    end interface


    !-------------------------------------------------------
    !-----                                             -----
    !-----     GENERALIZED MINIMUM RESIDUAL METHOD     -----
    !-----                                             -----
    !--------------------------------------------------------

    !----- Options and Metadata -----
    type, extends(abstract_opts), public :: gmres_sp_opts
        !! GMRES options.
        integer :: kdim = 30
        !! Dimension of the Krylov subspace (default: 30).
        integer :: maxiter = 10
        !! Maximum number of `gmres` restarts (default: 10).
        logical :: if_print_metadata = .false.
        !! Print iteration metadata on exit (default: .false.).
        logical :: sanity_check = .true.
        !! Performs extra matrix-vector product for sanity check.
    end type

    type, extends(abstract_metadata), public :: gmres_sp_metadata
        !! GMRES metadata.
        integer :: n_iter = 0
        !! Total iteration counter.
        integer :: n_inner = 0
        !! Number of inner iterations.
        integer :: n_outer = 0
        !! Number of outer iterations (i.e. restarts)
        real(sp), dimension(:), allocatable :: res
        !! Residual history.
        logical :: converged = .false.
        !! Convergence flag.
        integer :: info = 0
        !! Copy of the information flag for completeness.
    contains
        procedure, pass(self), public :: print => print_gmres_sp
        procedure, pass(self), public :: reset => reset_gmres_sp
    end type

    interface
        module subroutine print_gmres_sp(self, reset_counters, verbose)
            class(gmres_sp_metadata), intent(inout) :: self
            logical, optional, intent(in) :: reset_counters
            !! Reset all counters to zero after printing?
            logical, optional, intent(in) :: verbose
            !! Print the full residual history?
        end subroutine

        module subroutine reset_gmres_sp(self)
            class(gmres_sp_metadata), intent(inout) :: self
        end subroutine
    end interface

    type, extends(abstract_opts), public :: gmres_dp_opts
        !! GMRES options.
        integer :: kdim = 30
        !! Dimension of the Krylov subspace (default: 30).
        integer :: maxiter = 10
        !! Maximum number of `gmres` restarts (default: 10).
        logical :: if_print_metadata = .false.
        !! Print iteration metadata on exit (default: .false.).
        logical :: sanity_check = .true.
        !! Performs extra matrix-vector product for sanity check.
    end type

    type, extends(abstract_metadata), public :: gmres_dp_metadata
        !! GMRES metadata.
        integer :: n_iter = 0
        !! Total iteration counter.
        integer :: n_inner = 0
        !! Number of inner iterations.
        integer :: n_outer = 0
        !! Number of outer iterations (i.e. restarts)
        real(dp), dimension(:), allocatable :: res
        !! Residual history.
        logical :: converged = .false.
        !! Convergence flag.
        integer :: info = 0
        !! Copy of the information flag for completeness.
    contains
        procedure, pass(self), public :: print => print_gmres_dp
        procedure, pass(self), public :: reset => reset_gmres_dp
    end type

    interface
        module subroutine print_gmres_dp(self, reset_counters, verbose)
            class(gmres_dp_metadata), intent(inout) :: self
            logical, optional, intent(in) :: reset_counters
            !! Reset all counters to zero after printing?
            logical, optional, intent(in) :: verbose
            !! Print the full residual history?
        end subroutine

        module subroutine reset_gmres_dp(self)
            class(gmres_dp_metadata), intent(inout) :: self
        end subroutine
    end interface


    !----- Interfaces for the GMRES solvers -----
    interface gmres
        !!  ### Description
        !!
        !!  Solve a square linear system of equations
        !!
        !!  \[
        !!      Ax = b
        !!  \]
        !!
        !!  using the *Generalized Minimum RESidual* (GMRES) method.
        !!
        !!  **References**
        !!
        !!  - Saad Y. and Schultz M. H. "GMRES: A generalized minimal residual algorithm for
        !!  solving nonsymmetric linear systems." SIAM Journal on Scientific and Statistical
        !!  Computing, 7(3), 1986.
        !!
        !!  ### Syntax
        !!
        !!  ```fortran
        !!      call gmres(A, b, x, info [, rtol] [, atol] [, preconditioner] [, options] [, transpose])
        !!  ```
        !!
        !!  ### Arguments
        !!
        !!  `A` : Linear operator derived from one of the `abstract_linop` provided by the
        !!  `AbstractLinops` module. It is an `intent(inout)` argument.
        !!
        !!  `b` : Right-hand side vector derived from one the `abstract_vector` types provided
        !!  by the `AbstractVectors` module. It needs to have the same type and kind as `A`.
        !!  It is an `intent(in)` argument.
        !!
        !!  `x` : On entry, initial guess for the solution. On exit, the solution computed by
        !!  gmres. It is a vector derived from one the `abstract_vector` types provided by the
        !!  `AbstractVectors` module. It needs to have the same type and kind as `A`. It is
        !!  an `intent(inout)` argument.
        !!
        !!  `info` : `integer` information flag.
        !!
        !!  `rtol` (optional) : `real` relative tolerance for the solver.
        !!
        !!  `atol` (optional) : `real` absolute tolerance for the solver.
        !!
        !!  `preconditioner` (optional) : Right preconditioner used to solve the system. It needs
        !!  to be consistent with the `abstract_preconditioner` interface. It is an `intent(in)`
        !!  argument.
        !!
        !!  `options` (optional) : Container for the gmres options given by the `gmres_opts` type.
        !!  It is an `intent(in)` argument.
        !!
        !!  `transpose` (optional) : `logical` flag controlling whether \( Ax = b\) or
        !!  \( A^H x = b \) is being solver.

        ! --- Interface for GMRES with Abstract Linops and Abstract Vectors ---
        module subroutine gmres_rsp(A, b, x, info, rtol, atol, preconditioner, options, transpose, meta)
            class(abstract_linop_rsp), intent(inout) :: A
            !! Linear operator to be inverted.
            class(abstract_vector_rsp), intent(in) :: b
            !! Right-hand side vector.
            class(abstract_vector_rsp), intent(inout) :: x
            !! Solution vector.
            integer, intent(out) :: info
            !! Information flag.
            real(sp), optional, intent(in) :: rtol
            !! Relative solver tolerance
            real(sp), optional, intent(in) :: atol
            !! Absolute solver tolerance
            class(abstract_precond_rsp), optional, intent(in) :: preconditioner
            !! Preconditioner (optional).
            class(abstract_opts), optional, intent(in) :: options
            !! GMRES options.   
            logical, optional, intent(in) :: transpose
            !! Whether \(\mathbf{A}\) or \(\mathbf{A}^H\) is being used.
            class(abstract_metadata), optional, intent(out) :: meta
            !! Metadata.
        end subroutine
        module subroutine gmres_rdp(A, b, x, info, rtol, atol, preconditioner, options, transpose, meta)
            class(abstract_linop_rdp), intent(inout) :: A
            !! Linear operator to be inverted.
            class(abstract_vector_rdp), intent(in) :: b
            !! Right-hand side vector.
            class(abstract_vector_rdp), intent(inout) :: x
            !! Solution vector.
            integer, intent(out) :: info
            !! Information flag.
            real(dp), optional, intent(in) :: rtol
            !! Relative solver tolerance
            real(dp), optional, intent(in) :: atol
            !! Absolute solver tolerance
            class(abstract_precond_rdp), optional, intent(in) :: preconditioner
            !! Preconditioner (optional).
            class(abstract_opts), optional, intent(in) :: options
            !! GMRES options.   
            logical, optional, intent(in) :: transpose
            !! Whether \(\mathbf{A}\) or \(\mathbf{A}^H\) is being used.
            class(abstract_metadata), optional, intent(out) :: meta
            !! Metadata.
        end subroutine
        module subroutine gmres_csp(A, b, x, info, rtol, atol, preconditioner, options, transpose, meta)
            class(abstract_linop_csp), intent(inout) :: A
            !! Linear operator to be inverted.
            class(abstract_vector_csp), intent(in) :: b
            !! Right-hand side vector.
            class(abstract_vector_csp), intent(inout) :: x
            !! Solution vector.
            integer, intent(out) :: info
            !! Information flag.
            real(sp), optional, intent(in) :: rtol
            !! Relative solver tolerance
            real(sp), optional, intent(in) :: atol
            !! Absolute solver tolerance
            class(abstract_precond_csp), optional, intent(in) :: preconditioner
            !! Preconditioner (optional).
            class(abstract_opts), optional, intent(in) :: options
            !! GMRES options.   
            logical, optional, intent(in) :: transpose
            !! Whether \(\mathbf{A}\) or \(\mathbf{A}^H\) is being used.
            class(abstract_metadata), optional, intent(out) :: meta
            !! Metadata.
        end subroutine
        module subroutine gmres_cdp(A, b, x, info, rtol, atol, preconditioner, options, transpose, meta)
            class(abstract_linop_cdp), intent(inout) :: A
            !! Linear operator to be inverted.
            class(abstract_vector_cdp), intent(in) :: b
            !! Right-hand side vector.
            class(abstract_vector_cdp), intent(inout) :: x
            !! Solution vector.
            integer, intent(out) :: info
            !! Information flag.
            real(dp), optional, intent(in) :: rtol
            !! Relative solver tolerance
            real(dp), optional, intent(in) :: atol
            !! Absolute solver tolerance
            class(abstract_precond_cdp), optional, intent(in) :: preconditioner
            !! Preconditioner (optional).
            class(abstract_opts), optional, intent(in) :: options
            !! GMRES options.   
            logical, optional, intent(in) :: transpose
            !! Whether \(\mathbf{A}\) or \(\mathbf{A}^H\) is being used.
            class(abstract_metadata), optional, intent(out) :: meta
            !! Metadata.
        end subroutine
    end interface

    !----------------------------------------------------------------
    !-----                                                      -----
    !-----     FLEXIBLE GENERALIZED MINIMUM RESIDUAL METHOD     -----
    !-----                                                      -----
    !----------------------------------------------------------------

    !----- Options and Metadata -----
    type, extends(abstract_opts), public :: fgmres_sp_opts
        !! FGMRES options.
        integer :: kdim = 30
        !! Dimension of the Krylov subspace (default: 30).
        integer :: maxiter = 10
        !! Maximum number of `fgmres` restarts (default: 10).
        logical :: if_print_metadata = .false.
        !! Print iteration metadata on exit (default: .false.).
        logical :: sanity_check = .true.
        !! Performs extra matrix-vector product for sanity check.
    end type

    type, extends(abstract_metadata), public :: fgmres_sp_metadata
        !! FGMRES metadata.
        integer :: n_iter = 0
        !! Total iteration counter.
        integer :: n_inner = 0
        !! Number of inner iterations.
        integer :: n_outer = 0
        !! Number of outer iterations (i.e. restarts)
        real(sp), dimension(:), allocatable :: res
        !! Residual history.
        logical :: converged = .false.
        !! Convergence flag.
        integer :: info = 0
        !! Copy of the information flag for completeness.
    contains
        procedure, pass(self), public :: print => print_fgmres_sp
        procedure, pass(self), public :: reset => reset_fgmres_sp
    end type

    interface
        module subroutine print_fgmres_sp(self, reset_counters, verbose)
            class(fgmres_sp_metadata), intent(inout) :: self
            logical, optional, intent(in) :: reset_counters
            !! Reset all counters to zero after printing?
            logical, optional, intent(in) :: verbose
            !! Print the full residual history?
        end subroutine

        module subroutine reset_fgmres_sp(self)
            class(fgmres_sp_metadata), intent(inout) :: self
        end subroutine
    end interface

    type, extends(abstract_opts), public :: fgmres_dp_opts
        !! FGMRES options.
        integer :: kdim = 30
        !! Dimension of the Krylov subspace (default: 30).
        integer :: maxiter = 10
        !! Maximum number of `fgmres` restarts (default: 10).
        logical :: if_print_metadata = .false.
        !! Print iteration metadata on exit (default: .false.).
        logical :: sanity_check = .true.
        !! Performs extra matrix-vector product for sanity check.
    end type

    type, extends(abstract_metadata), public :: fgmres_dp_metadata
        !! FGMRES metadata.
        integer :: n_iter = 0
        !! Total iteration counter.
        integer :: n_inner = 0
        !! Number of inner iterations.
        integer :: n_outer = 0
        !! Number of outer iterations (i.e. restarts)
        real(dp), dimension(:), allocatable :: res
        !! Residual history.
        logical :: converged = .false.
        !! Convergence flag.
        integer :: info = 0
        !! Copy of the information flag for completeness.
    contains
        procedure, pass(self), public :: print => print_fgmres_dp
        procedure, pass(self), public :: reset => reset_fgmres_dp
    end type

    interface
        module subroutine print_fgmres_dp(self, reset_counters, verbose)
            class(fgmres_dp_metadata), intent(inout) :: self
            logical, optional, intent(in) :: reset_counters
            !! Reset all counters to zero after printing?
            logical, optional, intent(in) :: verbose
            !! Print the full residual history?
        end subroutine

        module subroutine reset_fgmres_dp(self)
            class(fgmres_dp_metadata), intent(inout) :: self
        end subroutine
    end interface


    !----- Interfaces for the FGMRES solvers -----
    interface fgmres
        !!  ### Description
        !!
        !!  Solve a square linear system of equations
        !!
        !!  \[
        !!      Ax = b
        !!  \]
        !!
        !!  using the *Flexible Generalized Minimum RESidual* (FGMRES) method.
        !!
        !!  **References**
        !!
        !!  - Saad Y. and Schultz M. H. "GMRES: A generalized minimal residual algorithm for
        !!  solving nonsymmetric linear systems." SIAM Journal on Scientific and Statistical
        !!  Computing, 7(3), 1986.
        !!
        !!  ### Syntax
        !!
        !!  ```fortran
        !!      call fgmres(A, b, x, info [, rtol] [, atol] [, preconditioner] [, options] [, transpose])
        !!  ```
        !!
        !!  ### Arguments
        !!
        !!  `A` : Linear operator derived from one of the `abstract_linop` provided by the
        !!  `AbstractLinops` module. It is an `intent(inout)` argument.
        !!
        !!  `b` : Right-hand side vector derived from one the `abstract_vector` types provided
        !!  by the `AbstractVectors` module. It needs to have the same type and kind as `A`.
        !!  It is an `intent(in)` argument.
        !!
        !!  `x` : On entry, initial guess for the solution. On exit, the solution computed by
        !!  gmres. It is a vector derived from one the `abstract_vector` types provided by the
        !!  `AbstractVectors` module. It needs to have the same type and kind as `A`. It is
        !!  an `intent(inout)` argument.
        !!
        !!  `info` : `integer` information flag.
        !!
        !!  `rtol` (optional) : `real` relative tolerance for the solver.
        !!
        !!  `atol` (optional) : `real` absolute tolerance for the solver.
        !!
        !!  `preconditioner` (optional) : Right preconditioner used to solve the system. It needs
        !!  to be consistent with the `abstract_preconditioner` interface. It is an `intent(in)`
        !!  argument.
        !!
        !!  `options` (optional) : Container for the gmres options given by the `gmres_opts` type.
        !!  It is an `intent(in)` argument.
        !!
        !!  `transpose` (optional) : `logical` flag controlling whether \( Ax = b\) or
        !!  \( A^H x = b \) is being solved.
        module subroutine fgmres_rsp(A, b, x, info, rtol, atol, preconditioner, options, transpose, meta)
            class(abstract_linop_rsp), intent(inout) :: A
            !! Linear operator to be inverted.
            class(abstract_vector_rsp), intent(in) :: b
            !! Right-hand side vector.
            class(abstract_vector_rsp), intent(inout) :: x
            !! Solution vector.
            integer, intent(out) :: info
            !! Information flag.
            real(sp), optional, intent(in) :: rtol
            !! Relative solver tolerance
            real(sp), optional, intent(in) :: atol
            !! Absolute solver tolerance
            class(abstract_precond_rsp), optional, intent(in) :: preconditioner
            !! Preconditioner (optional).
            class(abstract_opts), optional, intent(in) :: options
            !! GMRES options.   
            logical, optional, intent(in) :: transpose
            !! Whether \(\mathbf{A}\) or \(\mathbf{A}^H\) is being used.
            class(abstract_metadata), optional, intent(out) :: meta
            !! Metadata.
        end subroutine
        module subroutine fgmres_rdp(A, b, x, info, rtol, atol, preconditioner, options, transpose, meta)
            class(abstract_linop_rdp), intent(inout) :: A
            !! Linear operator to be inverted.
            class(abstract_vector_rdp), intent(in) :: b
            !! Right-hand side vector.
            class(abstract_vector_rdp), intent(inout) :: x
            !! Solution vector.
            integer, intent(out) :: info
            !! Information flag.
            real(dp), optional, intent(in) :: rtol
            !! Relative solver tolerance
            real(dp), optional, intent(in) :: atol
            !! Absolute solver tolerance
            class(abstract_precond_rdp), optional, intent(in) :: preconditioner
            !! Preconditioner (optional).
            class(abstract_opts), optional, intent(in) :: options
            !! GMRES options.   
            logical, optional, intent(in) :: transpose
            !! Whether \(\mathbf{A}\) or \(\mathbf{A}^H\) is being used.
            class(abstract_metadata), optional, intent(out) :: meta
            !! Metadata.
        end subroutine
        module subroutine fgmres_csp(A, b, x, info, rtol, atol, preconditioner, options, transpose, meta)
            class(abstract_linop_csp), intent(inout) :: A
            !! Linear operator to be inverted.
            class(abstract_vector_csp), intent(in) :: b
            !! Right-hand side vector.
            class(abstract_vector_csp), intent(inout) :: x
            !! Solution vector.
            integer, intent(out) :: info
            !! Information flag.
            real(sp), optional, intent(in) :: rtol
            !! Relative solver tolerance
            real(sp), optional, intent(in) :: atol
            !! Absolute solver tolerance
            class(abstract_precond_csp), optional, intent(in) :: preconditioner
            !! Preconditioner (optional).
            class(abstract_opts), optional, intent(in) :: options
            !! GMRES options.   
            logical, optional, intent(in) :: transpose
            !! Whether \(\mathbf{A}\) or \(\mathbf{A}^H\) is being used.
            class(abstract_metadata), optional, intent(out) :: meta
            !! Metadata.
        end subroutine
        module subroutine fgmres_cdp(A, b, x, info, rtol, atol, preconditioner, options, transpose, meta)
            class(abstract_linop_cdp), intent(inout) :: A
            !! Linear operator to be inverted.
            class(abstract_vector_cdp), intent(in) :: b
            !! Right-hand side vector.
            class(abstract_vector_cdp), intent(inout) :: x
            !! Solution vector.
            integer, intent(out) :: info
            !! Information flag.
            real(dp), optional, intent(in) :: rtol
            !! Relative solver tolerance
            real(dp), optional, intent(in) :: atol
            !! Absolute solver tolerance
            class(abstract_precond_cdp), optional, intent(in) :: preconditioner
            !! Preconditioner (optional).
            class(abstract_opts), optional, intent(in) :: options
            !! GMRES options.   
            logical, optional, intent(in) :: transpose
            !! Whether \(\mathbf{A}\) or \(\mathbf{A}^H\) is being used.
            class(abstract_metadata), optional, intent(out) :: meta
            !! Metadata.
        end subroutine
    end interface

    !---------------------------------------------
    !-----                                   -----
    !-----     CONJUGATE GRADIENT METHOD     -----
    !-----                                   -----
    !---------------------------------------------

    !----- Options and Metadata -----
    type, extends(abstract_opts), public :: cg_sp_opts
        !! Conjugate gradient options.
        integer :: maxiter = 100
        !! Maximum number of `cg` iterations (default: 100).
        logical :: if_print_metadata = .false.
        !! Print interation metadata on exit (default = .false.)
    end type

    type, extends(abstract_metadata), public :: cg_sp_metadata
        !! Conjugate gradient metadata.
        integer :: n_iter = 0
        !! Iteration counter
        real(sp), dimension(:), allocatable :: res
        !! Residual history
        logical :: converged = .false.
        !! Convergence flag
        integer :: info = 0
        !! Copy of the information flag for completeness
    contains
        procedure, pass(self), public :: print => print_cg_sp
        procedure, pass(self), public :: reset => reset_cg_sp
    end type

    interface
        module subroutine print_cg_sp(self, reset_counters, verbose)
            class(cg_sp_metadata), intent(inout) :: self
            logical, optional, intent(in) :: reset_counters
            !! Reset all counters to zero after printing?
            logical, optional, intent(in) :: verbose
            !! Print the residual full residual history?
        end subroutine

        module subroutine reset_cg_sp(self)
            class(cg_sp_metadata), intent(inout) :: self
        end subroutine
    end interface
    type, extends(abstract_opts), public :: cg_dp_opts
        !! Conjugate gradient options.
        integer :: maxiter = 100
        !! Maximum number of `cg` iterations (default: 100).
        logical :: if_print_metadata = .false.
        !! Print interation metadata on exit (default = .false.)
    end type

    type, extends(abstract_metadata), public :: cg_dp_metadata
        !! Conjugate gradient metadata.
        integer :: n_iter = 0
        !! Iteration counter
        real(dp), dimension(:), allocatable :: res
        !! Residual history
        logical :: converged = .false.
        !! Convergence flag
        integer :: info = 0
        !! Copy of the information flag for completeness
    contains
        procedure, pass(self), public :: print => print_cg_dp
        procedure, pass(self), public :: reset => reset_cg_dp
    end type

    interface
        module subroutine print_cg_dp(self, reset_counters, verbose)
            class(cg_dp_metadata), intent(inout) :: self
            logical, optional, intent(in) :: reset_counters
            !! Reset all counters to zero after printing?
            logical, optional, intent(in) :: verbose
            !! Print the residual full residual history?
        end subroutine

        module subroutine reset_cg_dp(self)
            class(cg_dp_metadata), intent(inout) :: self
        end subroutine
    end interface

    !----- Interfaces for the Conjugate Gradient solvers -----
    interface cg
        !!  ### Description
        !!
        !!  Given a symmetric (positive definite) matrix \( A \), solves the linear system
        !!
        !!  \[
        !!      Ax = b
        !!  \]
        !!
        !!  using the *Conjugate Gradient* method.
        !!
        !!  **References**
        !!
        !!  - Hestenes, M. R., and Stiefel, E. (1952). "Methods of Conjugate Gradients for Solving
        !!  Linear Systems," Journal of Research of the National Bureau of Standards,
        !!  49(6), 409–436.
        !!
        !!  ### Syntax
        !!
        !!  ```fortran
        !!      call cg(A, b, x, info [, rtol] [, atol] [, preconditioner] [, options])
        !!  ```
        !!
        !!  ### Arguments
        !!
        !!  `A` : Linear operator derived from one of the `abstract_sym_linop` or
        !!  `abstract_hermitian_linop` types provided by the `AbstractLinops` module. It is an
        !!  `intent(inout)` argument.
        !!
        !!  `b` : Right-hand side vector derived from one the `abstract_vector` types provided
        !!  by the `AbstractVectors` module. It needs to have the same type and kind as `A`.
        !!  It is an `intent(in)` argument.
        !!
        !!  `x` : On entry, initial guess for the solution. On exit, the solution computed by
        !!  cg. It is a vector derived from one the `abstract_vector` types provided by the
        !!  `AbstractVectors` module. It needs to have the same type and kind as `A`. It is
        !!  an `intent(inout)` argument.
        !!
        !!  `info` : `integer` information flag.
        !!
        !!  `rtol` (optional) : `real` relative tolerance for the solver.
        !!
        !!  `atol` (optional) : `real` absolute tolerance for the solver.
        !!
        !!  `preconditioner` (optional) : Right preconditioner used to solve the system. It needs
        !!  to be consistent with the `abstract_preconditioner` interface. It is an `intent(in)`
        !!  argument.
        !!
        !!  `options` (optional) : Container for the gmres options given by the `cg_opts` type.
        !!  It is an `intent(in)` argument.
        !!
        !!  @note
        !!  Although the interface to pass a preconditioner is exposed, it is not currently
        !!  implemented.
        !!  @endnote
        module subroutine cg_rsp(A, b, x, info, rtol, atol, preconditioner, options, meta)
            class(abstract_sym_linop_rsp), intent(inout) :: A
            !! Linear operator to be inverted.
            class(abstract_vector_rsp), intent(in) :: b
            !! Right-hand side vector.
            class(abstract_vector_rsp), intent(inout) :: x
            !! Solution vector.
            integer, intent(out) :: info
            !! Information flag.
            real(sp), optional, intent(in) :: rtol
            !! Relative solver tolerance
            real(sp), optional, intent(in) :: atol
            !! Absolute solver tolerance
            class(abstract_precond_rsp), optional, intent(in) :: preconditioner
            !! Preconditioner (not yet supported).
            type(cg_sp_opts), optional, intent(in) :: options
            !! Options for the conjugate gradient solver.
            class(abstract_metadata), optional, intent(out) :: meta
            !! Metadata.
        end subroutine
        module subroutine cg_rdp(A, b, x, info, rtol, atol, preconditioner, options, meta)
            class(abstract_sym_linop_rdp), intent(inout) :: A
            !! Linear operator to be inverted.
            class(abstract_vector_rdp), intent(in) :: b
            !! Right-hand side vector.
            class(abstract_vector_rdp), intent(inout) :: x
            !! Solution vector.
            integer, intent(out) :: info
            !! Information flag.
            real(dp), optional, intent(in) :: rtol
            !! Relative solver tolerance
            real(dp), optional, intent(in) :: atol
            !! Absolute solver tolerance
            class(abstract_precond_rdp), optional, intent(in) :: preconditioner
            !! Preconditioner (not yet supported).
            type(cg_dp_opts), optional, intent(in) :: options
            !! Options for the conjugate gradient solver.
            class(abstract_metadata), optional, intent(out) :: meta
            !! Metadata.
        end subroutine
        module subroutine cg_csp(A, b, x, info, rtol, atol, preconditioner, options, meta)
            class(abstract_hermitian_linop_csp), intent(inout) :: A
            !! Linear operator to be inverted.
            class(abstract_vector_csp), intent(in) :: b
            !! Right-hand side vector.
            class(abstract_vector_csp), intent(inout) :: x
            !! Solution vector.
            integer, intent(out) :: info
            !! Information flag.
            real(sp), optional, intent(in) :: rtol
            !! Relative solver tolerance
            real(sp), optional, intent(in) :: atol
            !! Absolute solver tolerance
            class(abstract_precond_csp), optional, intent(in) :: preconditioner
            !! Preconditioner (not yet supported).
            type(cg_sp_opts), optional, intent(in) :: options
            !! Options for the conjugate gradient solver.
            class(abstract_metadata), optional, intent(out) :: meta
            !! Metadata.
        end subroutine
        module subroutine cg_cdp(A, b, x, info, rtol, atol, preconditioner, options, meta)
            class(abstract_hermitian_linop_cdp), intent(inout) :: A
            !! Linear operator to be inverted.
            class(abstract_vector_cdp), intent(in) :: b
            !! Right-hand side vector.
            class(abstract_vector_cdp), intent(inout) :: x
            !! Solution vector.
            integer, intent(out) :: info
            !! Information flag.
            real(dp), optional, intent(in) :: rtol
            !! Relative solver tolerance
            real(dp), optional, intent(in) :: atol
            !! Absolute solver tolerance
            class(abstract_precond_cdp), optional, intent(in) :: preconditioner
            !! Preconditioner (not yet supported).
            type(cg_dp_opts), optional, intent(in) :: options
            !! Options for the conjugate gradient solver.
            class(abstract_metadata), optional, intent(out) :: meta
            !! Metadata.
        end subroutine
    end interface

    !------------------------------------------
    !-----                                -----
    !-----     SINGULAR VALUE SOLVERS     -----
    !-----                                -----
    !------------------------------------------

    interface svds
        !!  ### Description
        !!
        !!  Computes the leading singular triplets of an arbitrary linear operator \(A\)
        !!  using the Lanczos iterative process. Given a linear operator \(A\), it finds
        !!  the leading singular values and singular vectors such that:
        !!
        !!  \[
        !!      \begin{aligned}
        !!      Av & = \sigma u \\
        !!      A^H u & = \sigma v.
        !!      \end{aligned}
        !!  \]
        !!
        !!  The subspaces \(U\) and \(V\) are constructed via Lanczos factorization, resulting in
        !!  a bidiagonal matrix \(B\). The singular values of \(A\) are approximated by those of
        !!  \(B\) and the singular vectors are computed accordingly.
        !!
        !!  **References**
        !!
        !!  - Golub, G. H., & Kahan, W. (1965). "Calculating the Singular Values and
        !!   Pseudo-Inverse of a Matrix."
        !!  - Baglama, J., & Reichel, L. (2005). "Augmented implicitly restarted Lanczos
        !!   bidiagonalization methods."
        !!  - R. M. Larsen. "Lanczos bidiagonalization with partial reorthogonalization."
        !!   Technical Report, 1998.
        !!
        !!  ### Syntax
        !!
        !!  ```fortran
        !!      call svds(A, U, S, V, residuals, info [, kdim] [,tolerance])
        !!  ```
        !!
        !!  ### Arguments
        !!
        !!  `A` : Linear operator derived from `abstract_sym_linop_rsp`, `abstract_sym_linop_rdp`,
        !!        `abstract_hermitian_linop_csp` or `abstract_hermitian_linop_cdp` whose leading
        !!        eigenpairs need to be computed. It is an `intent(inout)` argument.
        !!
        !!  `U` : Array of `abstract_vectors` with the same type and kind as `A`. On exit, it
        !!        contains the left singular vectors of `A`. Note that the dimension of `U` fixes
        !!        the number of eigenpairs computed.
        !!
        !!  `S` : Rank-1 array of `real` numbers. On exit, it contains the leading
        !!        singular values of `A`. It is an `intent(out)` argument.
        !!
        !!  `V` : Array of `abstract_vectors` with the same type and kind as `A`. On exit, it
        !!        contains the left singular vectors of `A`. Note that the dimension of `U` fixes
        !!        the number of eigenpairs computed.
        !!
        !!  `residuals` : Rank-1 array of `real` numbers. On exit, it contains the residuals
        !!                associated with each singular triplet. It is an `intent(out)` argument.
        !!
        !!  `info` : `integer` Information flag.
        !!
        !!  `kdim` (*optional*) : `integer`, maximum dimension of the Krylov subspace used to
        !!                        approximate the leading singular triplets. It is an `intent(in)`
        !!                        argument. By default, `kdim = 4*size(X)`.
        !!
        !!  `tolerance` (*optional*) : `real` tolerance below which a triplet is considered as
        !!                             being converged. It is an `intent(in)` agument. By default,
        !!                             `tolerance = rtol_sp` or `tolerance = rtol_dp`.
        module subroutine svds_rsp(A, U, S, V, residuals, info, u0, kdim, tolerance)
        class(abstract_linop_rsp), intent(inout) :: A
        !! Linear operator whose leading singular triplets need to be computed.
        class(abstract_vector_rsp), intent(out) :: U(:)
        !! Leading left singular vectors.
        real(sp), allocatable, intent(out) :: S(:)
        !! Leading singular values.
        class(abstract_vector_rsp), intent(out) :: V(:)
        !! Leading right singular vectors.
        real(sp), allocatable, intent(out) :: residuals(:)
        !! Residuals associated to each Ritz eigenpair.
        integer, intent(out) :: info
        !! Information flag.
        class(abstract_vector_rsp), optional, intent(in) :: u0
        integer, optional, intent(in) :: kdim
        !! Desired number of eigenpairs.
        real(sp), optional, intent(in) :: tolerance
        !! Tolerance.
        end subroutine
        module subroutine svds_rdp(A, U, S, V, residuals, info, u0, kdim, tolerance)
        class(abstract_linop_rdp), intent(inout) :: A
        !! Linear operator whose leading singular triplets need to be computed.
        class(abstract_vector_rdp), intent(out) :: U(:)
        !! Leading left singular vectors.
        real(dp), allocatable, intent(out) :: S(:)
        !! Leading singular values.
        class(abstract_vector_rdp), intent(out) :: V(:)
        !! Leading right singular vectors.
        real(dp), allocatable, intent(out) :: residuals(:)
        !! Residuals associated to each Ritz eigenpair.
        integer, intent(out) :: info
        !! Information flag.
        class(abstract_vector_rdp), optional, intent(in) :: u0
        integer, optional, intent(in) :: kdim
        !! Desired number of eigenpairs.
        real(dp), optional, intent(in) :: tolerance
        !! Tolerance.
        end subroutine
        module subroutine svds_csp(A, U, S, V, residuals, info, u0, kdim, tolerance)
        class(abstract_linop_csp), intent(inout) :: A
        !! Linear operator whose leading singular triplets need to be computed.
        class(abstract_vector_csp), intent(out) :: U(:)
        !! Leading left singular vectors.
        real(sp), allocatable, intent(out) :: S(:)
        !! Leading singular values.
        class(abstract_vector_csp), intent(out) :: V(:)
        !! Leading right singular vectors.
        real(sp), allocatable, intent(out) :: residuals(:)
        !! Residuals associated to each Ritz eigenpair.
        integer, intent(out) :: info
        !! Information flag.
        class(abstract_vector_csp), optional, intent(in) :: u0
        integer, optional, intent(in) :: kdim
        !! Desired number of eigenpairs.
        real(sp), optional, intent(in) :: tolerance
        !! Tolerance.
        end subroutine
        module subroutine svds_cdp(A, U, S, V, residuals, info, u0, kdim, tolerance)
        class(abstract_linop_cdp), intent(inout) :: A
        !! Linear operator whose leading singular triplets need to be computed.
        class(abstract_vector_cdp), intent(out) :: U(:)
        !! Leading left singular vectors.
        real(dp), allocatable, intent(out) :: S(:)
        !! Leading singular values.
        class(abstract_vector_cdp), intent(out) :: V(:)
        !! Leading right singular vectors.
        real(dp), allocatable, intent(out) :: residuals(:)
        !! Residuals associated to each Ritz eigenpair.
        integer, intent(out) :: info
        !! Information flag.
        class(abstract_vector_cdp), optional, intent(in) :: u0
        integer, optional, intent(in) :: kdim
        !! Desired number of eigenpairs.
        real(dp), optional, intent(in) :: tolerance
        !! Tolerance.
        end subroutine
    end interface
































    interface save_eigenspectrum
        !!  ### Description
        !!
        !!  Utility function to save the eigenspectrum computed from the Arnoldi factorization.
        !!  It outpost a .npy file.
        !!
        !!  ### Syntax
        !!
        !!  ```fortran
        !!      call save_eigenspectrum(eigvals, residuals, fname)
        !!  ```
        !!
        !!  ### Arguments
        !!
        !!  `eigvals` : `complex` rank-1 array containing the eigenvalues.
        !!
        !!  `residuals` : `real` rank-1 array containing the residuals associated to each
        !!                eigenvalues.
        !!
        !!  `fname` : Name of the file to save the eigenspectrum.
        module procedure save_eigenspectrum_sp
        module procedure save_eigenspectrum_dp
    end interface

    interface eighs
        !!  ### Description
        !!
        !!  Computes the leading eigenpairs of a symmetric operator \(A\)
        !!  using the Lanczos iterative process. Given a square linear operator \(A\), it finds
        !!  the leading eigvalues and eigvectors such that:
        !!
        !!  \[
        !!      Ax = \lambda x
        !!  \]
        !!
        !!  The subspace \(X\) is constructed via Lanczos factorization, resulting in a symmetric
        !!  tridiagonal matrix \(T\). The eigenvalues of \(A\) are approximated by those of \(T\)
        !!  and the eigenvectors are computed accordingly.
        !!
        !!  **References**
        !!
        !!  - Lanczos, C. (1950). "An Iteration Method for the Solution of the Eigenvalue Problem
        !!  of Linear Differential and Integral Operators". United States Governm. Press Office.
        !!
        !!  ### Syntax
        !!
        !!  ```fortran
        !!      call eighs(A, X, eigvals, residuals, info [, kdim] [,tolerance])
        !!  ```
        !!
        !!  ### Arguments
        !!
        !!  `A` : Linear operator derived from `abstract_sym_linop_rsp`, `abstract_sym_linop_rdp`,
        !!        `abstract_hermitian_linop_csp` or `abstract_hermitian_linop_cdp` whose leading
        !!        eigenpairs need to be computed. It is an `intent(inout)` argument.
        !!
        !!  `X` : Array of `abstract_vectors` with the same type and kind as `A`. On exit, it
        !!        contains the leading eigenvectors of `A`. Note that the dimension of `X` fixes
        !!        the number of eigenpairs computed.
        !!
        !!  `eigvals` : Rank-1 array of `real` numbers. On exit, it contains the leading
        !!              eigenvalues of `A`. It is an `intent(out)` argument.
        !!
        !!  `residuals` : Rank-1 array of `real` numbers. On exit, it contains the residuals
        !!                associated with each eigenpairs. It is an `intent(out)` argument.
        !!
        !!  `info` : `integer` Information flag.
        !!
        !!  `kdim` (*optional*) : `integer`, maximum dimension of the Krylov subspace used to
        !!                        approximate the leading eigenpairs. It is an `intent(in)`
        !!                        argument. By default, `kdim = 4*size(X)`.
        !!
        !!  `tolerance` (*optional*) : `real` tolerance below which an eigenpair is considered as
        !!                             being converged. It is an `intent(in)` agument. By default,
        !!                             `tolerance = rtol_sp` or `tolerance = rtol_dp`.
        module procedure eighs_rsp
        module procedure eighs_rdp
        module procedure eighs_csp
        module procedure eighs_cdp
    end interface

    interface eigs
        !!  ### Description
        !!
        !!  Computes the leading eigenpairs of a square linear operator \(A\)
        !!  using the Arnoldi iterative process. Given a square linear operator \(A\), it finds
        !!  the leading eigvalues and eigvectorss such that:
        !!
        !!  \[
        !!      Ax = \lambda x
        !!  \]
        !!
        !!  or
        !!
        !!  \[
        !!      A^H x = \lambda x.
        !!  \]
        !!
        !!  The subspace \(X\) is constructed via Arnoldi factorization, resulting in an upper
        !!  Hessenberg matrix \(H\). The eigenvalues of \(A\) are approximated by those of \(H\)
        !!  and the eigenvectors are computed accordingly.
        !!
        !!  **References**
        !!
        !!  - Arnoldi, W. E. (1951). "The Principle of Minimized Iterations in the Solution of
        !!    the Matrix Eigenvalue Problem." Quarterly of Applied Mathematics, 9(1), 17–29.
        !!
        !!  ### Syntax
        !!
        !!  ```fortran
        !!      call eigs(A, X, eigvals, residuals, info [, kdim] [, select] [,tolerance] [, transpose])
        !!  ```
        !!
        !!  ### Arguments
        !!
        !!  `A` : Linear operator derived from `abstract_sym_linop_rsp`, `abstract_sym_linop_rdp`,
        !!        `abstract_hermitian_linop_csp` or `abstract_hermitian_linop_cdp` whose leading
        !!        eigenpairs need to be computed. It is an `intent(inout)` argument.
        !!
        !!  `X` : Array of `abstract_vectors` with the same type and kind as `A`. On exit, it
        !!        contains the leading eigenvectors of `A`. Note that the dimension of `X` fixes
        !!        the number of eigenpairs computed.
        !!
        !!  `eigvals` : Rank-1 array of `real` numbers. On exit, it contains the leading
        !!              eigenvalues of `A`. It is an `intent(out)` argument.
        !!
        !!  `residuals` : Rank-1 array of `real` numbers. On exit, it contains the residuals
        !!                associated with each eigenpairs. It is an `intent(out)` argument.
        !!
        !!  `info` : `integer` Information flag.
        !!
        !!  `kdim` (*optional*) : `integer`, maximum dimension of the Krylov subspace used to
        !!                        approximate the leading eigenpairs. It is an `intent(in)`
        !!                        argument. By default, `kdim = 4*size(X)`.
        !!
        !!  `select` (*optional*) : Function to select which eigenvalues to compute.
        !!
        !!  `tolerance` (*optional*) : `real` tolerance below which an eigenpair is considered as
        !!                             being converged. It is an `intent(in)` agument. By default,
        !!                             `tolerance = rtol_sp` or `tolerance = rtol_dp`.
        !!
        !!  `transpose` (*optional*) : `logical` flag determining whether the eigenvalues of \(A\)
        !!                             or \(A^H\) need to be computed.
        module procedure eigs_rsp
        module procedure eigs_rdp
        module procedure eigs_csp
        module procedure eigs_cdp
    end interface





    !--------------------------------------------------------
    !-----                                              -----
    !-----     GENERIC INTERFACE FOR LINEAR SOLVERS     -----
    !-----                                              -----
    !--------------------------------------------------------

    abstract interface
        subroutine abstract_linear_solver_rsp(A, b, x, info, rtol, atol, preconditioner, options, transpose, meta)
            !! Abstract interface to use a user-defined linear solver in `LightKrylov`.
            import abstract_linop_rsp, abstract_vector_rsp, abstract_opts, abstract_metadata, abstract_precond_rsp, sp
            class(abstract_linop_rsp), intent(inout) :: A
            !! Linear operator to invert.
            class(abstract_vector_rsp), intent(in) :: b
            !! Right-hand side vector.
            class(abstract_vector_rsp), intent(inout) :: x
            !! Solution vector.
            integer, intent(out) :: info
            !! Information flag. In case of successful exit, the flag should return the number of iterations required for convergence.
            real(sp), optional, intent(in) :: rtol
            !! Relative solver tolerance
            real(sp), optional, intent(in) :: atol
            !! Absolute solver tolerance
            class(abstract_precond_rsp), optional, intent(in) :: preconditioner
            !! Preconditioner.
            class(abstract_opts), optional, intent(in) :: options
            !! Options passed to the linear solver.
            logical, optional, intent(in) :: transpose
            !! Determine whether \(\mathbf{A}\) (`.false.`) or \(\mathbf{A}^T\) (`.true.`) is being used.
            class(abstract_metadata), optional, intent(out) :: meta
            !! Metadata.
        end subroutine abstract_linear_solver_rsp

        subroutine abstract_linear_solver_rdp(A, b, x, info, rtol, atol, preconditioner, options, transpose, meta)
            !! Abstract interface to use a user-defined linear solver in `LightKrylov`.
            import abstract_linop_rdp, abstract_vector_rdp, abstract_opts, abstract_metadata, abstract_precond_rdp, dp
            class(abstract_linop_rdp), intent(inout) :: A
            !! Linear operator to invert.
            class(abstract_vector_rdp), intent(in) :: b
            !! Right-hand side vector.
            class(abstract_vector_rdp), intent(inout) :: x
            !! Solution vector.
            integer, intent(out) :: info
            !! Information flag. In case of successful exit, the flag should return the number of iterations required for convergence.
            real(dp), optional, intent(in) :: rtol
            !! Relative solver tolerance
            real(dp), optional, intent(in) :: atol
            !! Absolute solver tolerance
            class(abstract_precond_rdp), optional, intent(in) :: preconditioner
            !! Preconditioner.
            class(abstract_opts), optional, intent(in) :: options
            !! Options passed to the linear solver.
            logical, optional, intent(in) :: transpose
            !! Determine whether \(\mathbf{A}\) (`.false.`) or \(\mathbf{A}^T\) (`.true.`) is being used.
            class(abstract_metadata), optional, intent(out) :: meta
            !! Metadata.
        end subroutine abstract_linear_solver_rdp

        subroutine abstract_linear_solver_csp(A, b, x, info, rtol, atol, preconditioner, options, transpose, meta)
            !! Abstract interface to use a user-defined linear solver in `LightKrylov`.
            import abstract_linop_csp, abstract_vector_csp, abstract_opts, abstract_metadata, abstract_precond_csp, sp
            class(abstract_linop_csp), intent(inout) :: A
            !! Linear operator to invert.
            class(abstract_vector_csp), intent(in) :: b
            !! Right-hand side vector.
            class(abstract_vector_csp), intent(inout) :: x
            !! Solution vector.
            integer, intent(out) :: info
            !! Information flag. In case of successful exit, the flag should return the number of iterations required for convergence.
            real(sp), optional, intent(in) :: rtol
            !! Relative solver tolerance
            real(sp), optional, intent(in) :: atol
            !! Absolute solver tolerance
            class(abstract_precond_csp), optional, intent(in) :: preconditioner
            !! Preconditioner.
            class(abstract_opts), optional, intent(in) :: options
            !! Options passed to the linear solver.
            logical, optional, intent(in) :: transpose
            !! Determine whether \(\mathbf{A}\) (`.false.`) or \(\mathbf{A}^T\) (`.true.`) is being used.
            class(abstract_metadata), optional, intent(out) :: meta
            !! Metadata.
        end subroutine abstract_linear_solver_csp

        subroutine abstract_linear_solver_cdp(A, b, x, info, rtol, atol, preconditioner, options, transpose, meta)
            !! Abstract interface to use a user-defined linear solver in `LightKrylov`.
            import abstract_linop_cdp, abstract_vector_cdp, abstract_opts, abstract_metadata, abstract_precond_cdp, dp
            class(abstract_linop_cdp), intent(inout) :: A
            !! Linear operator to invert.
            class(abstract_vector_cdp), intent(in) :: b
            !! Right-hand side vector.
            class(abstract_vector_cdp), intent(inout) :: x
            !! Solution vector.
            integer, intent(out) :: info
            !! Information flag. In case of successful exit, the flag should return the number of iterations required for convergence.
            real(dp), optional, intent(in) :: rtol
            !! Relative solver tolerance
            real(dp), optional, intent(in) :: atol
            !! Absolute solver tolerance
            class(abstract_precond_cdp), optional, intent(in) :: preconditioner
            !! Preconditioner.
            class(abstract_opts), optional, intent(in) :: options
            !! Options passed to the linear solver.
            logical, optional, intent(in) :: transpose
            !! Determine whether \(\mathbf{A}\) (`.false.`) or \(\mathbf{A}^T\) (`.true.`) is being used.
            class(abstract_metadata), optional, intent(out) :: meta
            !! Metadata.
        end subroutine abstract_linear_solver_cdp

    end interface

contains

    !-------------------------------------
    !-----     UTILITY FUNCTIONS     -----
    !-------------------------------------

    elemental pure function compute_residual_rsp(beta, x) result(residual)
        !! Computes the residual associated with a Ritz eigenpair.
        real(sp), intent(in) :: beta
        !! Norm of the residual Krylov vector.
        real(sp), intent(in) :: x
        !! Last entry of the low-dimensional Ritz eigenvector.
        real(sp) :: residual
        !! Residual associated to the corresponding Ritz eigenpair.
        residual = abs(beta*x)
        return
    end function compute_residual_rsp

    elemental pure function compute_residual_rdp(beta, x) result(residual)
        !! Computes the residual associated with a Ritz eigenpair.
        real(dp), intent(in) :: beta
        !! Norm of the residual Krylov vector.
        real(dp), intent(in) :: x
        !! Last entry of the low-dimensional Ritz eigenvector.
        real(dp) :: residual
        !! Residual associated to the corresponding Ritz eigenpair.
        residual = abs(beta*x)
        return
    end function compute_residual_rdp

    elemental pure function compute_residual_csp(beta, x) result(residual)
        !! Computes the residual associated with a Ritz eigenpair.
        complex(sp), intent(in) :: beta
        !! Norm of the residual Krylov vector.
        complex(sp), intent(in) :: x
        !! Last entry of the low-dimensional Ritz eigenvector.
        real(sp) :: residual
        !! Residual associated to the corresponding Ritz eigenpair.
        residual = abs(beta*x)
        return
    end function compute_residual_csp

    elemental pure function compute_residual_cdp(beta, x) result(residual)
        !! Computes the residual associated with a Ritz eigenpair.
        complex(dp), intent(in) :: beta
        !! Norm of the residual Krylov vector.
        complex(dp), intent(in) :: x
        !! Last entry of the low-dimensional Ritz eigenvector.
        real(dp) :: residual
        !! Residual associated to the corresponding Ritz eigenpair.
        residual = abs(beta*x)
        return
    end function compute_residual_cdp


    subroutine save_eigenspectrum_sp(eigvals, residuals, fname)
        !! Saves the eigenspectrum and corresponding residuals to disk use the `npy` binary format.
        complex(sp), intent(in) :: eigvals(:)
        !! Eigenalues.
        real(sp), intent(in) :: residuals(:)
        !! Residual of the corresponding Ritz eigenpairs.
        character(len=*), intent(in) :: fname
        !! Name of the output file.

        ! Internal variables.
        real(sp) :: data(size(eigvals), 3)

        ! Store data.
        data(:, 1) = eigvals%re ; data(:, 2) = eigvals%im ; data(:, 3) = residuals
        ! Save the eigenspectrum to disk.
        call save_npy(fname, data)

        return
    end subroutine save_eigenspectrum_sp

    function median_eigvals_selector_sp(lambda) result(selected)
        complex(sp), intent(in) :: lambda(:)
        logical :: selected(size(lambda))
        selected = abs(lambda) > median(abs(lambda))
        return
    end function median_eigvals_selector_sp

    subroutine save_eigenspectrum_dp(eigvals, residuals, fname)
        !! Saves the eigenspectrum and corresponding residuals to disk use the `npy` binary format.
        complex(dp), intent(in) :: eigvals(:)
        !! Eigenalues.
        real(dp), intent(in) :: residuals(:)
        !! Residual of the corresponding Ritz eigenpairs.
        character(len=*), intent(in) :: fname
        !! Name of the output file.

        ! Internal variables.
        real(dp) :: data(size(eigvals), 3)

        ! Store data.
        data(:, 1) = eigvals%re ; data(:, 2) = eigvals%im ; data(:, 3) = residuals
        ! Save the eigenspectrum to disk.
        call save_npy(fname, data)

        return
    end subroutine save_eigenspectrum_dp

    function median_eigvals_selector_dp(lambda) result(selected)
        complex(dp), intent(in) :: lambda(:)
        logical :: selected(size(lambda))
        selected = abs(lambda) > median(abs(lambda))
        return
    end function median_eigvals_selector_dp


    !---------------------------------------------------
    !-----     GENERAL EIGENVALUE COMPUTATIONS     -----
    !---------------------------------------------------

    subroutine eigs_rsp(A, X, eigvals, residuals, info, x0, kdim, select, tolerance, transpose)
        class(abstract_linop_rsp), intent(inout) :: A
        !! Linear operator whose leading eigenpairs need to be computed.
        class(abstract_vector_rsp), intent(out) :: X(:)
        !! Leading eigenvectors of \(\mathbf{A}\).
        complex(sp), allocatable, intent(out) :: eigvals(:)
        !! Leading eigenvalues of \(\mathbf{A}\).
        real(sp), allocatable, intent(out) :: residuals(:)
        !! Residuals associated to each Ritz eigenpair.
        integer, intent(out) :: info
        !! Information flag.
        class(abstract_vector_rsp), optional, intent(in) :: x0
        !! Optional starting vector for generating the Krylov subspace.
        integer, optional, intent(in) :: kdim
        procedure(eigvals_select_sp), optional :: select
        !! Desired number of eigenpairs.
        real(sp), optional, intent(in) :: tolerance
        !! Tolerance.
        logical, optional, intent(in) :: transpose
        !! Determine whether \(\mathbf{A}\) or \(\mathbf{A}^H\) is being used.

        !--------------------------------------
        !-----     Internal variables     -----
        !--------------------------------------

        ! Krylov subspace and Krylov subspace dimension.
        class(abstract_vector_rsp), allocatable :: Xwrk(:)
        integer :: kdim_, kstart
        ! Hessenberg matrix.
        real(sp), allocatable :: H(:, :)
        ! Working arrays for the eigenvectors and eigenvalues.
        real(sp), allocatable :: eigvecs_wrk(:, :)
        complex(sp), allocatable :: eigvals_wrk(:)
        real(sp), allocatable :: residuals_wrk(:)
        ! Miscellaneous.
        integer :: nev, conv
        integer :: i, j, k, niter, krst
        real(sp) :: tol, x0_norm
        real(sp) :: beta
        real(sp) :: alpha
        character(len=256) :: msg
        ! Eigenvalue selection.
        procedure(eigvals_select_sp), pointer :: select_

        if (time_lightkrylov()) call timer%start('eigs_rsp')
        ! Deals with optional parameters.
        nev = size(X)
        kdim_   = optval(kdim, 4*nev)
        tol     = optval(tolerance, rtol_sp)

        if (present(select)) then
            select_ => select
        else
            select_ => median_eigvals_selector_sp
        endif

        ! Allocate eigenvalues.
        allocate(eigvals(nev)) ; eigvals = 0.0_sp

        ! Allocate working variables.
        allocate(Xwrk(kdim_+1), source=X(1)) ; call zero_basis(Xwrk)
        if (present(x0)) then
            call copy(Xwrk(1), x0)
            x0_norm = x0%norm(); call Xwrk(1)%scal(one_rsp/x0_norm)
        else
            call Xwrk(1)%rand(.true.)
        endif
        allocate(H(kdim_+1, kdim_)) ; H = 0.0_sp
        allocate(eigvecs_wrk(kdim_, kdim_)) ; eigvecs_wrk = 0.0_sp
        allocate(eigvals_wrk(kdim_)) ; eigvals_wrk = 0.0_sp
        allocate(residuals_wrk(kdim_)) ; residuals_wrk = 0.0_sp

        ! Ritz eigenpairs computation.
        H = 0.0_sp

        kstart = 1 ; conv = 0 ; niter = 0 ; krst = 1
        krylovschur: do while (conv < nev)

           arnoldi_factorization: do k = kstart, kdim_
                ! Arnoldi step.
                call arnoldi(A, Xwrk, H, info, kstart=k, kend=k, transpose=transpose)
                call check_info(info, 'arnoldi', module=this_module, procedure='eigs_rsp')

                ! Spectral decomposition of the k x k Hessenberg matrix.
                eigvals_wrk = 0.0_sp ; eigvecs_wrk = 0.0_sp
                call eig(H(:k, :k), eigvecs_wrk(:k, :k), eigvals_wrk(:k))

                ! Compute residuals.
                beta = H(k+1, k)
                do i = 1, k
                    if (eigvals_wrk(i)%im > 0) then
                        alpha = abs(cmplx(eigvecs_wrk(k, i), eigvecs_wrk(k, i+1), kind=sp))
                    else if (eigvals_wrk(i)%im < 0) then
                        alpha = abs(cmplx(eigvecs_wrk(k, i-1), eigvecs_wrk(k, i), kind=sp))
                    else
                        alpha = abs(eigvecs_wrk(k, i))
                    endif
                    residuals_wrk(i) = compute_residual_rsp(beta, alpha)
                enddo

                ! Check convergence.
                niter = niter + 1
                conv = count(residuals_wrk(:k) < tol)
                write(msg,'(I0,A,I0,A,I0,A)') conv, '/', nev, ' eigenvalues converged after ', niter, &
                            & ' steps of the Arnoldi process.'
                call logger%log_information(msg, module=this_module, procedure='eigs_rsp')
                if (conv >= nev) exit arnoldi_factorization
            enddo arnoldi_factorization

            write(msg,'(I0,A,I0,A,I0,A)') conv, '/', nev, ' eigenvalues converged after ', krst, &
                            & ' Krylov-Schur restarts of the Arnoldi process.'
            call logger%log_information(msg, module=this_module, procedure='eigs_rsp')
            ! Krylov-Schur restarting procedure.
            krst  = krst + 1
            call krylov_schur(kstart, Xwrk, H, select_) ; kstart = kstart + 1
            
        end do krylovschur

        !--------------------------------
        !-----     POST-PROCESS     -----
        !--------------------------------

        block
        integer :: indices(kdim_)
        real(sp) :: abs_eigvals(kdim_)
       
        ! Re-compute eigenvalues and eigenvectors.
        k = min(k, kdim_) ; call eig(H(:k, :k), eigvecs_wrk(:k, :k), eigvals_wrk(:k))
        ! Sort eigenvalues.
        abs_eigvals = abs(eigvals_wrk) ; call sort_index(abs_eigvals, indices, reverse=.true.)
        eigvals_wrk = eigvals_wrk(indices) ; eigvecs_wrk = eigvecs_wrk(:, indices)
        residuals_wrk = residuals_wrk(indices)

        ! Store converged eigenvalues.
        eigvals = eigvals_wrk(:nev) ; residuals = residuals_wrk(:nev)
        end block

        ! Construct eigenvectors.
        do i = 1, nev
            call X(i)%zero()
            do j = 1, k
                call X(i)%axpby(one_rsp, Xwrk(j), eigvecs_wrk(j, i))
            enddo
        enddo

        info = niter
        if (time_lightkrylov()) call timer%stop('eigs_rsp')
        
        return
    end subroutine eigs_rsp

    subroutine eigs_rdp(A, X, eigvals, residuals, info, x0, kdim, select, tolerance, transpose)
        class(abstract_linop_rdp), intent(inout) :: A
        !! Linear operator whose leading eigenpairs need to be computed.
        class(abstract_vector_rdp), intent(out) :: X(:)
        !! Leading eigenvectors of \(\mathbf{A}\).
        complex(dp), allocatable, intent(out) :: eigvals(:)
        !! Leading eigenvalues of \(\mathbf{A}\).
        real(dp), allocatable, intent(out) :: residuals(:)
        !! Residuals associated to each Ritz eigenpair.
        integer, intent(out) :: info
        !! Information flag.
        class(abstract_vector_rdp), optional, intent(in) :: x0
        !! Optional starting vector for generating the Krylov subspace.
        integer, optional, intent(in) :: kdim
        procedure(eigvals_select_dp), optional :: select
        !! Desired number of eigenpairs.
        real(dp), optional, intent(in) :: tolerance
        !! Tolerance.
        logical, optional, intent(in) :: transpose
        !! Determine whether \(\mathbf{A}\) or \(\mathbf{A}^H\) is being used.

        !--------------------------------------
        !-----     Internal variables     -----
        !--------------------------------------

        ! Krylov subspace and Krylov subspace dimension.
        class(abstract_vector_rdp), allocatable :: Xwrk(:)
        integer :: kdim_, kstart
        ! Hessenberg matrix.
        real(dp), allocatable :: H(:, :)
        ! Working arrays for the eigenvectors and eigenvalues.
        real(dp), allocatable :: eigvecs_wrk(:, :)
        complex(dp), allocatable :: eigvals_wrk(:)
        real(dp), allocatable :: residuals_wrk(:)
        ! Miscellaneous.
        integer :: nev, conv
        integer :: i, j, k, niter, krst
        real(dp) :: tol, x0_norm
        real(dp) :: beta
        real(dp) :: alpha
        character(len=256) :: msg
        ! Eigenvalue selection.
        procedure(eigvals_select_dp), pointer :: select_

        if (time_lightkrylov()) call timer%start('eigs_rdp')
        ! Deals with optional parameters.
        nev = size(X)
        kdim_   = optval(kdim, 4*nev)
        tol     = optval(tolerance, rtol_dp)

        if (present(select)) then
            select_ => select
        else
            select_ => median_eigvals_selector_dp
        endif

        ! Allocate eigenvalues.
        allocate(eigvals(nev)) ; eigvals = 0.0_dp

        ! Allocate working variables.
        allocate(Xwrk(kdim_+1), source=X(1)) ; call zero_basis(Xwrk)
        if (present(x0)) then
            call copy(Xwrk(1), x0)
            x0_norm = x0%norm(); call Xwrk(1)%scal(one_rdp/x0_norm)
        else
            call Xwrk(1)%rand(.true.)
        endif
        allocate(H(kdim_+1, kdim_)) ; H = 0.0_dp
        allocate(eigvecs_wrk(kdim_, kdim_)) ; eigvecs_wrk = 0.0_dp
        allocate(eigvals_wrk(kdim_)) ; eigvals_wrk = 0.0_dp
        allocate(residuals_wrk(kdim_)) ; residuals_wrk = 0.0_dp

        ! Ritz eigenpairs computation.
        H = 0.0_dp

        kstart = 1 ; conv = 0 ; niter = 0 ; krst = 1
        krylovschur: do while (conv < nev)

           arnoldi_factorization: do k = kstart, kdim_
                ! Arnoldi step.
                call arnoldi(A, Xwrk, H, info, kstart=k, kend=k, transpose=transpose)
                call check_info(info, 'arnoldi', module=this_module, procedure='eigs_rdp')

                ! Spectral decomposition of the k x k Hessenberg matrix.
                eigvals_wrk = 0.0_dp ; eigvecs_wrk = 0.0_dp
                call eig(H(:k, :k), eigvecs_wrk(:k, :k), eigvals_wrk(:k))

                ! Compute residuals.
                beta = H(k+1, k)
                do i = 1, k
                    if (eigvals_wrk(i)%im > 0) then
                        alpha = abs(cmplx(eigvecs_wrk(k, i), eigvecs_wrk(k, i+1), kind=dp))
                    else if (eigvals_wrk(i)%im < 0) then
                        alpha = abs(cmplx(eigvecs_wrk(k, i-1), eigvecs_wrk(k, i), kind=dp))
                    else
                        alpha = abs(eigvecs_wrk(k, i))
                    endif
                    residuals_wrk(i) = compute_residual_rdp(beta, alpha)
                enddo

                ! Check convergence.
                niter = niter + 1
                conv = count(residuals_wrk(:k) < tol)
                write(msg,'(I0,A,I0,A,I0,A)') conv, '/', nev, ' eigenvalues converged after ', niter, &
                            & ' steps of the Arnoldi process.'
                call logger%log_information(msg, module=this_module, procedure='eigs_rdp')
                if (conv >= nev) exit arnoldi_factorization
            enddo arnoldi_factorization

            write(msg,'(I0,A,I0,A,I0,A)') conv, '/', nev, ' eigenvalues converged after ', krst, &
                            & ' Krylov-Schur restarts of the Arnoldi process.'
            call logger%log_information(msg, module=this_module, procedure='eigs_rdp')
            ! Krylov-Schur restarting procedure.
            krst  = krst + 1
            call krylov_schur(kstart, Xwrk, H, select_) ; kstart = kstart + 1
            
        end do krylovschur

        !--------------------------------
        !-----     POST-PROCESS     -----
        !--------------------------------

        block
        integer :: indices(kdim_)
        real(dp) :: abs_eigvals(kdim_)
       
        ! Re-compute eigenvalues and eigenvectors.
        k = min(k, kdim_) ; call eig(H(:k, :k), eigvecs_wrk(:k, :k), eigvals_wrk(:k))
        ! Sort eigenvalues.
        abs_eigvals = abs(eigvals_wrk) ; call sort_index(abs_eigvals, indices, reverse=.true.)
        eigvals_wrk = eigvals_wrk(indices) ; eigvecs_wrk = eigvecs_wrk(:, indices)
        residuals_wrk = residuals_wrk(indices)

        ! Store converged eigenvalues.
        eigvals = eigvals_wrk(:nev) ; residuals = residuals_wrk(:nev)
        end block

        ! Construct eigenvectors.
        do i = 1, nev
            call X(i)%zero()
            do j = 1, k
                call X(i)%axpby(one_rdp, Xwrk(j), eigvecs_wrk(j, i))
            enddo
        enddo

        info = niter
        if (time_lightkrylov()) call timer%stop('eigs_rdp')
        
        return
    end subroutine eigs_rdp

    subroutine eigs_csp(A, X, eigvals, residuals, info, x0, kdim, select, tolerance, transpose)
        class(abstract_linop_csp), intent(inout) :: A
        !! Linear operator whose leading eigenpairs need to be computed.
        class(abstract_vector_csp), intent(out) :: X(:)
        !! Leading eigenvectors of \(\mathbf{A}\).
        complex(sp), allocatable, intent(out) :: eigvals(:)
        !! Leading eigenvalues of \(\mathbf{A}\).
        real(sp), allocatable, intent(out) :: residuals(:)
        !! Residuals associated to each Ritz eigenpair.
        integer, intent(out) :: info
        !! Information flag.
        class(abstract_vector_csp), optional, intent(in) :: x0
        !! Optional starting vector for generating the Krylov subspace.
        integer, optional, intent(in) :: kdim
        procedure(eigvals_select_sp), optional :: select
        !! Desired number of eigenpairs.
        real(sp), optional, intent(in) :: tolerance
        !! Tolerance.
        logical, optional, intent(in) :: transpose
        !! Determine whether \(\mathbf{A}\) or \(\mathbf{A}^H\) is being used.

        !--------------------------------------
        !-----     Internal variables     -----
        !--------------------------------------

        ! Krylov subspace and Krylov subspace dimension.
        class(abstract_vector_csp), allocatable :: Xwrk(:)
        integer :: kdim_, kstart
        ! Hessenberg matrix.
        complex(sp), allocatable :: H(:, :)
        ! Working arrays for the eigenvectors and eigenvalues.
        complex(sp), allocatable :: eigvecs_wrk(:, :)
        complex(sp), allocatable :: eigvals_wrk(:)
        real(sp), allocatable :: residuals_wrk(:)
        ! Miscellaneous.
        integer :: nev, conv
        integer :: i, j, k, niter, krst
        real(sp) :: tol, x0_norm
        complex(sp) :: beta
        character(len=256) :: msg
        ! Eigenvalue selection.
        procedure(eigvals_select_sp), pointer :: select_

        if (time_lightkrylov()) call timer%start('eigs_csp')
        ! Deals with optional parameters.
        nev = size(X)
        kdim_   = optval(kdim, 4*nev)
        tol     = optval(tolerance, rtol_sp)

        if (present(select)) then
            select_ => select
        else
            select_ => median_eigvals_selector_sp
        endif

        ! Allocate eigenvalues.
        allocate(eigvals(nev)) ; eigvals = 0.0_sp

        ! Allocate working variables.
        allocate(Xwrk(kdim_+1), source=X(1)) ; call zero_basis(Xwrk)
        if (present(x0)) then
            call copy(Xwrk(1), x0)
            x0_norm = x0%norm(); call Xwrk(1)%scal(one_csp/x0_norm)
        else
            call Xwrk(1)%rand(.true.)
        endif
        allocate(H(kdim_+1, kdim_)) ; H = 0.0_sp
        allocate(eigvecs_wrk(kdim_, kdim_)) ; eigvecs_wrk = 0.0_sp
        allocate(eigvals_wrk(kdim_)) ; eigvals_wrk = 0.0_sp
        allocate(residuals_wrk(kdim_)) ; residuals_wrk = 0.0_sp

        ! Ritz eigenpairs computation.
        H = 0.0_sp

        kstart = 1 ; conv = 0 ; niter = 0 ; krst = 1
        krylovschur: do while (conv < nev)

           arnoldi_factorization: do k = kstart, kdim_
                ! Arnoldi step.
                call arnoldi(A, Xwrk, H, info, kstart=k, kend=k, transpose=transpose)
                call check_info(info, 'arnoldi', module=this_module, procedure='eigs_csp')

                ! Spectral decomposition of the k x k Hessenberg matrix.
                eigvals_wrk = 0.0_sp ; eigvecs_wrk = 0.0_sp
                call eig(H(:k, :k), eigvecs_wrk(:k, :k), eigvals_wrk(:k))

                ! Compute residuals.
                beta = H(k+1, k)
                residuals_wrk(:k) = compute_residual_csp(beta, eigvecs_wrk(k,:k))

                ! Check convergence.
                niter = niter + 1
                conv = count(residuals_wrk(:k) < tol)
                write(msg,'(I0,A,I0,A,I0,A)') conv, '/', nev, ' eigenvalues converged after ', niter, &
                            & ' steps of the Arnoldi process.'
                call logger%log_information(msg, module=this_module, procedure='eigs_csp')
                if (conv >= nev) exit arnoldi_factorization
            enddo arnoldi_factorization

            write(msg,'(I0,A,I0,A,I0,A)') conv, '/', nev, ' eigenvalues converged after ', krst, &
                            & ' Krylov-Schur restarts of the Arnoldi process.'
            call logger%log_information(msg, module=this_module, procedure='eigs_csp')
            ! Krylov-Schur restarting procedure.
            krst  = krst + 1
            call krylov_schur(kstart, Xwrk, H, select_) ; kstart = kstart + 1
            
        end do krylovschur

        !--------------------------------
        !-----     POST-PROCESS     -----
        !--------------------------------

        block
        integer :: indices(kdim_)
        real(sp) :: abs_eigvals(kdim_)
       
        ! Re-compute eigenvalues and eigenvectors.
        k = min(k, kdim_) ; call eig(H(:k, :k), eigvecs_wrk(:k, :k), eigvals_wrk(:k))
        ! Sort eigenvalues.
        abs_eigvals = abs(eigvals_wrk) ; call sort_index(abs_eigvals, indices, reverse=.true.)
        eigvals_wrk = eigvals_wrk(indices) ; eigvecs_wrk = eigvecs_wrk(:, indices)
        residuals_wrk = residuals_wrk(indices)

        ! Store converged eigenvalues.
        eigvals = eigvals_wrk(:nev) ; residuals = residuals_wrk(:nev)
        end block

        ! Construct eigenvectors.
        do i = 1, nev
            call X(i)%zero()
            do j = 1, k
                call X(i)%axpby(one_csp, Xwrk(j), eigvecs_wrk(j, i))
            enddo
        enddo

        info = niter
        if (time_lightkrylov()) call timer%stop('eigs_csp')
        
        return
    end subroutine eigs_csp

    subroutine eigs_cdp(A, X, eigvals, residuals, info, x0, kdim, select, tolerance, transpose)
        class(abstract_linop_cdp), intent(inout) :: A
        !! Linear operator whose leading eigenpairs need to be computed.
        class(abstract_vector_cdp), intent(out) :: X(:)
        !! Leading eigenvectors of \(\mathbf{A}\).
        complex(dp), allocatable, intent(out) :: eigvals(:)
        !! Leading eigenvalues of \(\mathbf{A}\).
        real(dp), allocatable, intent(out) :: residuals(:)
        !! Residuals associated to each Ritz eigenpair.
        integer, intent(out) :: info
        !! Information flag.
        class(abstract_vector_cdp), optional, intent(in) :: x0
        !! Optional starting vector for generating the Krylov subspace.
        integer, optional, intent(in) :: kdim
        procedure(eigvals_select_dp), optional :: select
        !! Desired number of eigenpairs.
        real(dp), optional, intent(in) :: tolerance
        !! Tolerance.
        logical, optional, intent(in) :: transpose
        !! Determine whether \(\mathbf{A}\) or \(\mathbf{A}^H\) is being used.

        !--------------------------------------
        !-----     Internal variables     -----
        !--------------------------------------

        ! Krylov subspace and Krylov subspace dimension.
        class(abstract_vector_cdp), allocatable :: Xwrk(:)
        integer :: kdim_, kstart
        ! Hessenberg matrix.
        complex(dp), allocatable :: H(:, :)
        ! Working arrays for the eigenvectors and eigenvalues.
        complex(dp), allocatable :: eigvecs_wrk(:, :)
        complex(dp), allocatable :: eigvals_wrk(:)
        real(dp), allocatable :: residuals_wrk(:)
        ! Miscellaneous.
        integer :: nev, conv
        integer :: i, j, k, niter, krst
        real(dp) :: tol, x0_norm
        complex(dp) :: beta
        character(len=256) :: msg
        ! Eigenvalue selection.
        procedure(eigvals_select_dp), pointer :: select_

        if (time_lightkrylov()) call timer%start('eigs_cdp')
        ! Deals with optional parameters.
        nev = size(X)
        kdim_   = optval(kdim, 4*nev)
        tol     = optval(tolerance, rtol_dp)

        if (present(select)) then
            select_ => select
        else
            select_ => median_eigvals_selector_dp
        endif

        ! Allocate eigenvalues.
        allocate(eigvals(nev)) ; eigvals = 0.0_dp

        ! Allocate working variables.
        allocate(Xwrk(kdim_+1), source=X(1)) ; call zero_basis(Xwrk)
        if (present(x0)) then
            call copy(Xwrk(1), x0)
            x0_norm = x0%norm(); call Xwrk(1)%scal(one_cdp/x0_norm)
        else
            call Xwrk(1)%rand(.true.)
        endif
        allocate(H(kdim_+1, kdim_)) ; H = 0.0_dp
        allocate(eigvecs_wrk(kdim_, kdim_)) ; eigvecs_wrk = 0.0_dp
        allocate(eigvals_wrk(kdim_)) ; eigvals_wrk = 0.0_dp
        allocate(residuals_wrk(kdim_)) ; residuals_wrk = 0.0_dp

        ! Ritz eigenpairs computation.
        H = 0.0_dp

        kstart = 1 ; conv = 0 ; niter = 0 ; krst = 1
        krylovschur: do while (conv < nev)

           arnoldi_factorization: do k = kstart, kdim_
                ! Arnoldi step.
                call arnoldi(A, Xwrk, H, info, kstart=k, kend=k, transpose=transpose)
                call check_info(info, 'arnoldi', module=this_module, procedure='eigs_cdp')

                ! Spectral decomposition of the k x k Hessenberg matrix.
                eigvals_wrk = 0.0_dp ; eigvecs_wrk = 0.0_dp
                call eig(H(:k, :k), eigvecs_wrk(:k, :k), eigvals_wrk(:k))

                ! Compute residuals.
                beta = H(k+1, k)
                residuals_wrk(:k) = compute_residual_cdp(beta, eigvecs_wrk(k,:k))

                ! Check convergence.
                niter = niter + 1
                conv = count(residuals_wrk(:k) < tol)
                write(msg,'(I0,A,I0,A,I0,A)') conv, '/', nev, ' eigenvalues converged after ', niter, &
                            & ' steps of the Arnoldi process.'
                call logger%log_information(msg, module=this_module, procedure='eigs_cdp')
                if (conv >= nev) exit arnoldi_factorization
            enddo arnoldi_factorization

            write(msg,'(I0,A,I0,A,I0,A)') conv, '/', nev, ' eigenvalues converged after ', krst, &
                            & ' Krylov-Schur restarts of the Arnoldi process.'
            call logger%log_information(msg, module=this_module, procedure='eigs_cdp')
            ! Krylov-Schur restarting procedure.
            krst  = krst + 1
            call krylov_schur(kstart, Xwrk, H, select_) ; kstart = kstart + 1
            
        end do krylovschur

        !--------------------------------
        !-----     POST-PROCESS     -----
        !--------------------------------

        block
        integer :: indices(kdim_)
        real(dp) :: abs_eigvals(kdim_)
       
        ! Re-compute eigenvalues and eigenvectors.
        k = min(k, kdim_) ; call eig(H(:k, :k), eigvecs_wrk(:k, :k), eigvals_wrk(:k))
        ! Sort eigenvalues.
        abs_eigvals = abs(eigvals_wrk) ; call sort_index(abs_eigvals, indices, reverse=.true.)
        eigvals_wrk = eigvals_wrk(indices) ; eigvecs_wrk = eigvecs_wrk(:, indices)
        residuals_wrk = residuals_wrk(indices)

        ! Store converged eigenvalues.
        eigvals = eigvals_wrk(:nev) ; residuals = residuals_wrk(:nev)
        end block

        ! Construct eigenvectors.
        do i = 1, nev
            call X(i)%zero()
            do j = 1, k
                call X(i)%axpby(one_cdp, Xwrk(j), eigvecs_wrk(j, i))
            enddo
        enddo

        info = niter
        if (time_lightkrylov()) call timer%stop('eigs_cdp')
        
        return
    end subroutine eigs_cdp


    !-----------------------------------------------------------------------------
    !-----      EIGENVALUE COMPUTATIONS FOR SYMMETRIC/HERMITIAN MATRICES     -----
    !-----------------------------------------------------------------------------

    subroutine eighs_rsp(A, X, eigvals, residuals, info, x0, kdim, tolerance)
        class(abstract_sym_linop_rsp), intent(inout) :: A
        !! Linear operator whose leading eigenpairs need to be computed.
        class(abstract_vector_rsp), intent(out) :: X(:)
        !! Leading eigevectors of \( \mathbf{A} \).
        real(sp), allocatable, intent(out) :: eigvals(:)
        !! Leading eigenvalues of \( \mathbf{A} \).
        real(sp), allocatable, intent(out) :: residuals(:)
        !! Residuals associated to each Ritz eigenpairs.
        integer, intent(out) :: info
        !! Information flag.
        class(abstract_vector_rsp), optional, intent(in) :: x0
        !! Optional starting vector to generate the Krylov subspace.
        integer, optional, intent(in) :: kdim
        !! Desired number of eigenpairs.
        real(sp), optional, intent(in) :: tolerance
        !! Tolerance

        !--------------------------------------
        !-----     Internal variables     -----
        !--------------------------------------

        class(abstract_vector_rsp), allocatable :: Xwrk(:)
        ! Krylov subspace.
        integer :: kdim_
        ! Krylov subspace dimension.
        real(sp), allocatable :: T(:, :)
        ! Tridiagonal matrix.
        real(sp), allocatable :: eigvecs_wrk(:, :)
        ! Working array for the Ritz eigenvectors.
        real(sp), allocatable :: eigvals_wrk(:)
        ! Working array for the Ritz eigenvalues.
        real(sp), allocatable :: residuals_wrk(:)
        ! Working array for the Ritz residuals.
        real(sp) :: x0_norm

        ! Miscellaneous.
        integer :: i, j, k, nev, conv
        real(sp) :: tol
        real(sp) :: beta
        character(len=256) :: msg

        if (time_lightkrylov()) call timer%start('eighs_rsp')
        ! Deaks with the optional args.
        nev = size(X)
        kdim_ = optval(kdim, 4*nev)
        tol = optval(tolerance, rtol_sp)

        ! Allocate working variables.
        allocate(Xwrk(kdim_+1), source=X(1)) ; call zero_basis(Xwrk)
        if (present(x0)) then
            call copy(Xwrk(1), x0)
            x0_norm = x0%norm(); call Xwrk(1)%scal(one_rsp/x0_norm)
        else
            call Xwrk(1)%rand(.true.)
        endif
        allocate(T(kdim_+1, kdim_)) ; T = zero_rsp
        allocate(eigvecs_wrk(kdim_, kdim_)) ; eigvecs_wrk = zero_rsp
        allocate(eigvals_wrk(kdim_)) ; eigvals_wrk = 0.0_sp
        allocate(residuals_wrk(kdim_)) ; residuals_wrk = 0.0_sp

        ! Ritz eigenpairs computation.
        lanczos_iter : do k = 1, kdim_
            ! Symmetric Lanczos step.
            call lanczos(A, Xwrk, T, info, kstart=k, kend=k)
            call check_info(info, 'lanczos', module=this_module, procedure='eighs_rsp')

            ! Spectral decomposition of the k x k tridiagonal matrix.
            eigvals_wrk = 0.0_sp ; eigvecs_wrk = zero_rsp
            call eigh(T(:k, :k), eigvals_wrk(:k), vectors=eigvecs_wrk(:k, :k))

            ! Compute residuals.
            beta = T(k+1, k)
            residuals_wrk(:k) = compute_residual_rsp(beta, eigvecs_wrk(k, :k))

            ! Check convergence.
            conv = count(residuals_wrk(:k) < tol)
            write(msg,'(I0,A,I0,A,I0,A)') conv, '/', nev, ' eigenvalues converged after ', k, &
                            & ' iterations of the Lanczos process.'
            call logger%log_information(msg, module=this_module, procedure='eighs_rsp')
            if (conv >= nev) exit lanczos_iter
        enddo lanczos_iter

        !--------------------------------
        !-----     POST-PROCESS     -----
        !--------------------------------

        block
            integer :: indices(kdim_)
            call sort_index(eigvals_wrk, indices, reverse=.true.)
            !eigvals_wrk = eigvals_wrk(indices) ; 
            eigvecs_wrk = eigvecs_wrk(:, indices)
            ! Store converged eigenvalues.
            eigvals = eigvals_wrk(:nev) ; residuals = residuals_wrk(:nev)
        end block

        ! Construct eigenvectors.
        k = min(k, kdim_)
        do i = 1, nev
            call X(i)%zero()
            do j = 1, k
                call X(i)%axpby(one_rsp, Xwrk(j), eigvecs_wrk(j, i))
            enddo
        enddo
        
        info = k
        if (time_lightkrylov()) call timer%stop('eighs_rsp')
        
        return
    end subroutine eighs_rsp

    subroutine eighs_rdp(A, X, eigvals, residuals, info, x0, kdim, tolerance)
        class(abstract_sym_linop_rdp), intent(inout) :: A
        !! Linear operator whose leading eigenpairs need to be computed.
        class(abstract_vector_rdp), intent(out) :: X(:)
        !! Leading eigevectors of \( \mathbf{A} \).
        real(dp), allocatable, intent(out) :: eigvals(:)
        !! Leading eigenvalues of \( \mathbf{A} \).
        real(dp), allocatable, intent(out) :: residuals(:)
        !! Residuals associated to each Ritz eigenpairs.
        integer, intent(out) :: info
        !! Information flag.
        class(abstract_vector_rdp), optional, intent(in) :: x0
        !! Optional starting vector to generate the Krylov subspace.
        integer, optional, intent(in) :: kdim
        !! Desired number of eigenpairs.
        real(dp), optional, intent(in) :: tolerance
        !! Tolerance

        !--------------------------------------
        !-----     Internal variables     -----
        !--------------------------------------

        class(abstract_vector_rdp), allocatable :: Xwrk(:)
        ! Krylov subspace.
        integer :: kdim_
        ! Krylov subspace dimension.
        real(dp), allocatable :: T(:, :)
        ! Tridiagonal matrix.
        real(dp), allocatable :: eigvecs_wrk(:, :)
        ! Working array for the Ritz eigenvectors.
        real(dp), allocatable :: eigvals_wrk(:)
        ! Working array for the Ritz eigenvalues.
        real(dp), allocatable :: residuals_wrk(:)
        ! Working array for the Ritz residuals.
        real(dp) :: x0_norm

        ! Miscellaneous.
        integer :: i, j, k, nev, conv
        real(dp) :: tol
        real(dp) :: beta
        character(len=256) :: msg

        if (time_lightkrylov()) call timer%start('eighs_rdp')
        ! Deaks with the optional args.
        nev = size(X)
        kdim_ = optval(kdim, 4*nev)
        tol = optval(tolerance, rtol_dp)

        ! Allocate working variables.
        allocate(Xwrk(kdim_+1), source=X(1)) ; call zero_basis(Xwrk)
        if (present(x0)) then
            call copy(Xwrk(1), x0)
            x0_norm = x0%norm(); call Xwrk(1)%scal(one_rdp/x0_norm)
        else
            call Xwrk(1)%rand(.true.)
        endif
        allocate(T(kdim_+1, kdim_)) ; T = zero_rdp
        allocate(eigvecs_wrk(kdim_, kdim_)) ; eigvecs_wrk = zero_rdp
        allocate(eigvals_wrk(kdim_)) ; eigvals_wrk = 0.0_dp
        allocate(residuals_wrk(kdim_)) ; residuals_wrk = 0.0_dp

        ! Ritz eigenpairs computation.
        lanczos_iter : do k = 1, kdim_
            ! Symmetric Lanczos step.
            call lanczos(A, Xwrk, T, info, kstart=k, kend=k)
            call check_info(info, 'lanczos', module=this_module, procedure='eighs_rdp')

            ! Spectral decomposition of the k x k tridiagonal matrix.
            eigvals_wrk = 0.0_dp ; eigvecs_wrk = zero_rdp
            call eigh(T(:k, :k), eigvals_wrk(:k), vectors=eigvecs_wrk(:k, :k))

            ! Compute residuals.
            beta = T(k+1, k)
            residuals_wrk(:k) = compute_residual_rdp(beta, eigvecs_wrk(k, :k))

            ! Check convergence.
            conv = count(residuals_wrk(:k) < tol)
            write(msg,'(I0,A,I0,A,I0,A)') conv, '/', nev, ' eigenvalues converged after ', k, &
                            & ' iterations of the Lanczos process.'
            call logger%log_information(msg, module=this_module, procedure='eighs_rdp')
            if (conv >= nev) exit lanczos_iter
        enddo lanczos_iter

        !--------------------------------
        !-----     POST-PROCESS     -----
        !--------------------------------

        block
            integer :: indices(kdim_)
            call sort_index(eigvals_wrk, indices, reverse=.true.)
            !eigvals_wrk = eigvals_wrk(indices) ; 
            eigvecs_wrk = eigvecs_wrk(:, indices)
            ! Store converged eigenvalues.
            eigvals = eigvals_wrk(:nev) ; residuals = residuals_wrk(:nev)
        end block

        ! Construct eigenvectors.
        k = min(k, kdim_)
        do i = 1, nev
            call X(i)%zero()
            do j = 1, k
                call X(i)%axpby(one_rdp, Xwrk(j), eigvecs_wrk(j, i))
            enddo
        enddo
        
        info = k
        if (time_lightkrylov()) call timer%stop('eighs_rdp')
        
        return
    end subroutine eighs_rdp

    subroutine eighs_csp(A, X, eigvals, residuals, info, x0, kdim, tolerance)
        class(abstract_hermitian_linop_csp), intent(inout) :: A
        !! Linear operator whose leading eigenpairs need to be computed.
        class(abstract_vector_csp), intent(out) :: X(:)
        !! Leading eigevectors of \( \mathbf{A} \).
        real(sp), allocatable, intent(out) :: eigvals(:)
        !! Leading eigenvalues of \( \mathbf{A} \).
        real(sp), allocatable, intent(out) :: residuals(:)
        !! Residuals associated to each Ritz eigenpairs.
        integer, intent(out) :: info
        !! Information flag.
        class(abstract_vector_csp), optional, intent(in) :: x0
        !! Optional starting vector to generate the Krylov subspace.
        integer, optional, intent(in) :: kdim
        !! Desired number of eigenpairs.
        real(sp), optional, intent(in) :: tolerance
        !! Tolerance

        !--------------------------------------
        !-----     Internal variables     -----
        !--------------------------------------

        class(abstract_vector_csp), allocatable :: Xwrk(:)
        ! Krylov subspace.
        integer :: kdim_
        ! Krylov subspace dimension.
        complex(sp), allocatable :: T(:, :)
        ! Tridiagonal matrix.
        complex(sp), allocatable :: eigvecs_wrk(:, :)
        ! Working array for the Ritz eigenvectors.
        real(sp), allocatable :: eigvals_wrk(:)
        ! Working array for the Ritz eigenvalues.
        real(sp), allocatable :: residuals_wrk(:)
        ! Working array for the Ritz residuals.
        real(sp) :: x0_norm

        ! Miscellaneous.
        integer :: i, j, k, nev, conv
        real(sp) :: tol
        complex(sp) :: beta
        character(len=256) :: msg

        if (time_lightkrylov()) call timer%start('eighs_csp')
        ! Deaks with the optional args.
        nev = size(X)
        kdim_ = optval(kdim, 4*nev)
        tol = optval(tolerance, rtol_sp)

        ! Allocate working variables.
        allocate(Xwrk(kdim_+1), source=X(1)) ; call zero_basis(Xwrk)
        if (present(x0)) then
            call copy(Xwrk(1), x0)
            x0_norm = x0%norm(); call Xwrk(1)%scal(one_csp/x0_norm)
        else
            call Xwrk(1)%rand(.true.)
        endif
        allocate(T(kdim_+1, kdim_)) ; T = zero_csp
        allocate(eigvecs_wrk(kdim_, kdim_)) ; eigvecs_wrk = zero_csp
        allocate(eigvals_wrk(kdim_)) ; eigvals_wrk = 0.0_sp
        allocate(residuals_wrk(kdim_)) ; residuals_wrk = 0.0_sp

        ! Ritz eigenpairs computation.
        lanczos_iter : do k = 1, kdim_
            ! Symmetric Lanczos step.
            call lanczos(A, Xwrk, T, info, kstart=k, kend=k)
            call check_info(info, 'lanczos', module=this_module, procedure='eighs_csp')

            ! Spectral decomposition of the k x k tridiagonal matrix.
            eigvals_wrk = 0.0_sp ; eigvecs_wrk = zero_csp
            call eigh(T(:k, :k), eigvals_wrk(:k), vectors=eigvecs_wrk(:k, :k))

            ! Compute residuals.
            beta = T(k+1, k)
            residuals_wrk(:k) = compute_residual_csp(beta, eigvecs_wrk(k, :k))

            ! Check convergence.
            conv = count(residuals_wrk(:k) < tol)
            write(msg,'(I0,A,I0,A,I0,A)') conv, '/', nev, ' eigenvalues converged after ', k, &
                            & ' iterations of the Lanczos process.'
            call logger%log_information(msg, module=this_module, procedure='eighs_csp')
            if (conv >= nev) exit lanczos_iter
        enddo lanczos_iter

        !--------------------------------
        !-----     POST-PROCESS     -----
        !--------------------------------

        block
            integer :: indices(kdim_)
            call sort_index(eigvals_wrk, indices, reverse=.true.)
            !eigvals_wrk = eigvals_wrk(indices) ; 
            eigvecs_wrk = eigvecs_wrk(:, indices)
            ! Store converged eigenvalues.
            eigvals = eigvals_wrk(:nev) ; residuals = residuals_wrk(:nev)
        end block

        ! Construct eigenvectors.
        k = min(k, kdim_)
        do i = 1, nev
            call X(i)%zero()
            do j = 1, k
                call X(i)%axpby(one_csp, Xwrk(j), eigvecs_wrk(j, i))
            enddo
        enddo
        
        info = k
        if (time_lightkrylov()) call timer%stop('eighs_csp')
        
        return
    end subroutine eighs_csp

    subroutine eighs_cdp(A, X, eigvals, residuals, info, x0, kdim, tolerance)
        class(abstract_hermitian_linop_cdp), intent(inout) :: A
        !! Linear operator whose leading eigenpairs need to be computed.
        class(abstract_vector_cdp), intent(out) :: X(:)
        !! Leading eigevectors of \( \mathbf{A} \).
        real(dp), allocatable, intent(out) :: eigvals(:)
        !! Leading eigenvalues of \( \mathbf{A} \).
        real(dp), allocatable, intent(out) :: residuals(:)
        !! Residuals associated to each Ritz eigenpairs.
        integer, intent(out) :: info
        !! Information flag.
        class(abstract_vector_cdp), optional, intent(in) :: x0
        !! Optional starting vector to generate the Krylov subspace.
        integer, optional, intent(in) :: kdim
        !! Desired number of eigenpairs.
        real(dp), optional, intent(in) :: tolerance
        !! Tolerance

        !--------------------------------------
        !-----     Internal variables     -----
        !--------------------------------------

        class(abstract_vector_cdp), allocatable :: Xwrk(:)
        ! Krylov subspace.
        integer :: kdim_
        ! Krylov subspace dimension.
        complex(dp), allocatable :: T(:, :)
        ! Tridiagonal matrix.
        complex(dp), allocatable :: eigvecs_wrk(:, :)
        ! Working array for the Ritz eigenvectors.
        real(dp), allocatable :: eigvals_wrk(:)
        ! Working array for the Ritz eigenvalues.
        real(dp), allocatable :: residuals_wrk(:)
        ! Working array for the Ritz residuals.
        real(dp) :: x0_norm

        ! Miscellaneous.
        integer :: i, j, k, nev, conv
        real(dp) :: tol
        complex(dp) :: beta
        character(len=256) :: msg

        if (time_lightkrylov()) call timer%start('eighs_cdp')
        ! Deaks with the optional args.
        nev = size(X)
        kdim_ = optval(kdim, 4*nev)
        tol = optval(tolerance, rtol_dp)

        ! Allocate working variables.
        allocate(Xwrk(kdim_+1), source=X(1)) ; call zero_basis(Xwrk)
        if (present(x0)) then
            call copy(Xwrk(1), x0)
            x0_norm = x0%norm(); call Xwrk(1)%scal(one_cdp/x0_norm)
        else
            call Xwrk(1)%rand(.true.)
        endif
        allocate(T(kdim_+1, kdim_)) ; T = zero_cdp
        allocate(eigvecs_wrk(kdim_, kdim_)) ; eigvecs_wrk = zero_cdp
        allocate(eigvals_wrk(kdim_)) ; eigvals_wrk = 0.0_dp
        allocate(residuals_wrk(kdim_)) ; residuals_wrk = 0.0_dp

        ! Ritz eigenpairs computation.
        lanczos_iter : do k = 1, kdim_
            ! Symmetric Lanczos step.
            call lanczos(A, Xwrk, T, info, kstart=k, kend=k)
            call check_info(info, 'lanczos', module=this_module, procedure='eighs_cdp')

            ! Spectral decomposition of the k x k tridiagonal matrix.
            eigvals_wrk = 0.0_dp ; eigvecs_wrk = zero_cdp
            call eigh(T(:k, :k), eigvals_wrk(:k), vectors=eigvecs_wrk(:k, :k))

            ! Compute residuals.
            beta = T(k+1, k)
            residuals_wrk(:k) = compute_residual_cdp(beta, eigvecs_wrk(k, :k))

            ! Check convergence.
            conv = count(residuals_wrk(:k) < tol)
            write(msg,'(I0,A,I0,A,I0,A)') conv, '/', nev, ' eigenvalues converged after ', k, &
                            & ' iterations of the Lanczos process.'
            call logger%log_information(msg, module=this_module, procedure='eighs_cdp')
            if (conv >= nev) exit lanczos_iter
        enddo lanczos_iter

        !--------------------------------
        !-----     POST-PROCESS     -----
        !--------------------------------

        block
            integer :: indices(kdim_)
            call sort_index(eigvals_wrk, indices, reverse=.true.)
            !eigvals_wrk = eigvals_wrk(indices) ; 
            eigvecs_wrk = eigvecs_wrk(:, indices)
            ! Store converged eigenvalues.
            eigvals = eigvals_wrk(:nev) ; residuals = residuals_wrk(:nev)
        end block

        ! Construct eigenvectors.
        k = min(k, kdim_)
        do i = 1, nev
            call X(i)%zero()
            do j = 1, k
                call X(i)%axpby(one_cdp, Xwrk(j), eigvecs_wrk(j, i))
            enddo
        enddo
        
        info = k
        if (time_lightkrylov()) call timer%stop('eighs_cdp')
        
        return
    end subroutine eighs_cdp


end module lightkrylov_IterativeSolvers
