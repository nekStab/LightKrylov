module LightKrylov_ExpmLib
    !!  This module implements the evaluation of the "matrix-exponential times vector" procedure
    !!  using Krylov methods.
    !!
    !!  - `kexpm` computes \( \exp(\tau \mathbf{A}) \mathbf{b} \) (or its block version) by projection
    !!    onto a single (block) Krylov subspace whose dimension is increased until the error estimate
    !!    drops below the user-defined tolerance.
    !!  - `kexpm_substep` splits the time horizon into adaptively chosen sub-steps (Expokit-like) so
    !!    that convergence no longer depends on \( \| \tau \mathbf{A} \| \) being small enough for a
    !!    fixed Krylov dimension. With a persistent metadata object it learns, over successive calls,
    !!    the cheapest number of uniform sub-steps (warm start) while keeping the error check.
    !!  - `krylov_exptA` is a wrapper of `kexpm_substep` conforming to the `abstract_exptA` interface.

    ! Fortran standard library.
    use stdlib_optval, only: optval
    use stdlib_linalg, only: mnorm, matrix_exp

    ! LightKrylov.
    use LightKrylov_Constants
    use LightKrylov_Logger
    use LightKrylov_Utils
    use LightKrylov_AbstractVectors
    use LightKrylov_AbstractLinops
    use LightKrylov_BaseKrylov

    implicit none(type, external)
    private

    character(len=*), parameter :: this_module      = 'LK_ExpmLib'
    character(len=*), parameter :: this_module_long = 'LightKrylov_ExpmLib'

    integer, parameter :: kmax = 100
    !! Default maximum number of Arnoldi steps in `kexpm`.

    public :: abstract_exptA_rsp
    public :: abstract_exptA_rdp
    public :: abstract_exptA_csp
    public :: abstract_exptA_cdp
    public :: kexpm
    public :: kexpm_substep
    public :: krylov_exptA
    public :: krylov_exptA_rsp
    public :: krylov_exptA_rdp
    public :: krylov_exptA_csp
    public :: krylov_exptA_cdp

    !-----------------------------------------------------------
    !-----     Options and metadata for `kexpm_substep`     -----
    !-----------------------------------------------------------

    type, extends(abstract_opts), public :: kexpm_sp_opts
        !! Options for the sub-stepped Krylov exponential `kexpm_substep`.
        integer :: kdim = 30
        !! Maximum number of Krylov vectors per sub-step (default: 30).
        logical :: warm_start = .true.
        !! Re-use the sub-step count learned during previous calls (stored in the metadata).
        !! Only active when a `kexpm_sp_metadata` object is passed (default: .true.).
        logical :: tune = .true.
        !! Over successive calls, decrease the number of uniform sub-steps while this lowers the
        !! number of matrix-vector products, and freeze it at the last setting that did not need
        !! the safety net (default: .true.). With `.false.`, the count of the cold call is kept.
        integer :: max_rejections = 25
        !! Maximum number of consecutive step-size reductions within one sub-step (default: 25).
        integer :: max_substeps = 10000
        !! Safeguard on the number of sub-steps per call (default: 10000).
        logical :: if_print_metadata = .false.
        !! Print the metadata on exit (default: .false.).
    end type kexpm_sp_opts

    type, extends(abstract_metadata), public :: kexpm_sp_metadata
        !! Counters and learned settings of `kexpm_substep`. Pass the same object to successive calls
        !! to benefit from the warm start. `reset()` clears the counters only, `forget()` clears the
        !! learned settings (the next call is then a cold, fully adaptive one).
        ! ----- Counters -----
        integer :: n_calls = 0
        !! Number of calls.
        integer :: n_matvec = 0
        !! Total number of matrix-vector products.
        integer :: n_substeps = 0
        !! Total number of sub-steps.
        integer :: n_rejections = 0
        !! Total number of rejected step sizes (each costs one small dense exponential).
        integer :: n_fallbacks = 0
        !! Number of calls in which the learned uniform step was too large (safety net triggered).
        integer :: info = 0
        !! Information flag of the last call.
        logical :: converged = .false.
        !! Convergence flag of the last call.
        ! ----- Learned settings (warm start) -----
        logical :: learned = .false.
        !! Whether settings have been learned.
        integer :: nsub = 0
        !! Number of uniform sub-steps used by the next call (0: fully adaptive).
        integer :: nsub_best = 0
        !! Best number of uniform sub-steps found so far (0: adaptive was best).
        integer :: mv_best = 0
        !! Matrix-vector products of the best setting.
        integer :: nsub0 = 0
        !! Number of sub-steps taken by the cold adaptive call (starting point of the search).
        logical :: frozen = .false.
        !! Whether the search is over.
        real(sp) :: tau_learned = 0.0_sp
        real(sp) :: tol_learned = 0.0_sp
        integer :: kdim_learned = 0
        !! Problem settings for which the learned values are valid: |tau|, tol/||b|| (within a
        !! factor 2) and kdim. Any change triggers `forget`.
    contains
        procedure, pass(self), public :: print => print_kexpm_sp
        procedure, pass(self), public :: reset => reset_kexpm_sp
        procedure, pass(self), public :: forget => forget_kexpm_sp
    end type kexpm_sp_metadata

    type, extends(abstract_opts), public :: kexpm_dp_opts
        !! Options for the sub-stepped Krylov exponential `kexpm_substep`.
        integer :: kdim = 30
        !! Maximum number of Krylov vectors per sub-step (default: 30).
        logical :: warm_start = .true.
        !! Re-use the sub-step count learned during previous calls (stored in the metadata).
        !! Only active when a `kexpm_dp_metadata` object is passed (default: .true.).
        logical :: tune = .true.
        !! Over successive calls, decrease the number of uniform sub-steps while this lowers the
        !! number of matrix-vector products, and freeze it at the last setting that did not need
        !! the safety net (default: .true.). With `.false.`, the count of the cold call is kept.
        integer :: max_rejections = 25
        !! Maximum number of consecutive step-size reductions within one sub-step (default: 25).
        integer :: max_substeps = 10000
        !! Safeguard on the number of sub-steps per call (default: 10000).
        logical :: if_print_metadata = .false.
        !! Print the metadata on exit (default: .false.).
    end type kexpm_dp_opts

    type, extends(abstract_metadata), public :: kexpm_dp_metadata
        !! Counters and learned settings of `kexpm_substep`. Pass the same object to successive calls
        !! to benefit from the warm start. `reset()` clears the counters only, `forget()` clears the
        !! learned settings (the next call is then a cold, fully adaptive one).
        ! ----- Counters -----
        integer :: n_calls = 0
        !! Number of calls.
        integer :: n_matvec = 0
        !! Total number of matrix-vector products.
        integer :: n_substeps = 0
        !! Total number of sub-steps.
        integer :: n_rejections = 0
        !! Total number of rejected step sizes (each costs one small dense exponential).
        integer :: n_fallbacks = 0
        !! Number of calls in which the learned uniform step was too large (safety net triggered).
        integer :: info = 0
        !! Information flag of the last call.
        logical :: converged = .false.
        !! Convergence flag of the last call.
        ! ----- Learned settings (warm start) -----
        logical :: learned = .false.
        !! Whether settings have been learned.
        integer :: nsub = 0
        !! Number of uniform sub-steps used by the next call (0: fully adaptive).
        integer :: nsub_best = 0
        !! Best number of uniform sub-steps found so far (0: adaptive was best).
        integer :: mv_best = 0
        !! Matrix-vector products of the best setting.
        integer :: nsub0 = 0
        !! Number of sub-steps taken by the cold adaptive call (starting point of the search).
        logical :: frozen = .false.
        !! Whether the search is over.
        real(dp) :: tau_learned = 0.0_dp
        real(dp) :: tol_learned = 0.0_dp
        integer :: kdim_learned = 0
        !! Problem settings for which the learned values are valid: |tau|, tol/||b|| (within a
        !! factor 2) and kdim. Any change triggers `forget`.
    contains
        procedure, pass(self), public :: print => print_kexpm_dp
        procedure, pass(self), public :: reset => reset_kexpm_dp
        procedure, pass(self), public :: forget => forget_kexpm_dp
    end type kexpm_dp_metadata


    abstract interface
        subroutine abstract_exptA_rsp(vec_out, A, vec_in, tau, info, trans)
            !! Abstract interface to define the matrix exponential-vector product.
            import sp
            import abstract_vector_rsp
            import abstract_linop_rsp
            implicit none(type, external)
            class(abstract_vector_rsp), intent(out) :: vec_out
            !! Solution vector.
            class(abstract_linop_rsp), intent(inout) :: A
            !! Linear operator to be exponentiated.
            class(abstract_vector_rsp), intent(in) :: vec_in
            !! Input vector.
            real(sp), intent(in) :: tau
            !! Time horizon for integration.
            integer, intent(out) :: info
            !! Information flag.
            logical, optional, intent(in) :: trans
            !! Use transpose ?
        end subroutine abstract_exptA_rsp

        subroutine abstract_exptA_rdp(vec_out, A, vec_in, tau, info, trans)
            !! Abstract interface to define the matrix exponential-vector product.
            import dp
            import abstract_vector_rdp
            import abstract_linop_rdp
            implicit none(type, external)
            class(abstract_vector_rdp), intent(out) :: vec_out
            !! Solution vector.
            class(abstract_linop_rdp), intent(inout) :: A
            !! Linear operator to be exponentiated.
            class(abstract_vector_rdp), intent(in) :: vec_in
            !! Input vector.
            real(dp), intent(in) :: tau
            !! Time horizon for integration.
            integer, intent(out) :: info
            !! Information flag.
            logical, optional, intent(in) :: trans
            !! Use transpose ?
        end subroutine abstract_exptA_rdp

        subroutine abstract_exptA_csp(vec_out, A, vec_in, tau, info, trans)
            !! Abstract interface to define the matrix exponential-vector product.
            import sp
            import abstract_vector_csp
            import abstract_linop_csp
            implicit none(type, external)
            class(abstract_vector_csp), intent(out) :: vec_out
            !! Solution vector.
            class(abstract_linop_csp), intent(inout) :: A
            !! Linear operator to be exponentiated.
            class(abstract_vector_csp), intent(in) :: vec_in
            !! Input vector.
            real(sp), intent(in) :: tau
            !! Time horizon for integration.
            integer, intent(out) :: info
            !! Information flag.
            logical, optional, intent(in) :: trans
            !! Use transpose ?
        end subroutine abstract_exptA_csp

        subroutine abstract_exptA_cdp(vec_out, A, vec_in, tau, info, trans)
            !! Abstract interface to define the matrix exponential-vector product.
            import dp
            import abstract_vector_cdp
            import abstract_linop_cdp
            implicit none(type, external)
            class(abstract_vector_cdp), intent(out) :: vec_out
            !! Solution vector.
            class(abstract_linop_cdp), intent(inout) :: A
            !! Linear operator to be exponentiated.
            class(abstract_vector_cdp), intent(in) :: vec_in
            !! Input vector.
            real(dp), intent(in) :: tau
            !! Time horizon for integration.
            integer, intent(out) :: info
            !! Information flag.
            logical, optional, intent(in) :: trans
            !! Use transpose ?
        end subroutine abstract_exptA_cdp

    end interface

    interface kexpm
        !!  ### Description
        !!
        !!  This interface provides methods to evaluate the matrix-vector product
        !!  \( c = \exp(\tau A) b \) based on the Arnoldi method.
        !!
        !!  ### Syntax
        !!
        !!  ```fortran
        !!      call kexpm(c, A, b, tau, tol, info [, trans] [, kdim] [, X])
        !!  ```
        !!
        !!  ### Arguments
        !!
        !!  - `c`   :   Output vector (or vectors). It is an `intent(out)` argument.
        !!
        !!  - `A`   :   Linear operator to be exponentiated. It is an `intent(inout)` argument.
        !!
        !!  - `b`   :   Vector (or vectors) to be multiplied by \( \exp(\tau A) \). It is an `intent(in)` argument.
        !!
        !!  - `tau` :   `real` (single or double) time over which the matrix exponential needs to
        !!              be computed. It is an `intent(in)` argument.
        !!
        !!  - `tol` :   Absolute tolerance on the (conservative) error estimate.
        !!
        !!  - `info`    :   `integer` Information flag. On exit:
        !!                  - `info > 0`  : converged, number of Krylov vectors used.
        !!                  - `info = -1` : the error estimate did not drop below `tol` within the
        !!                                  allowed Krylov dimension.
        !!                  - `info = -3` : invalid workspace `X` (fewer than `2p` vectors, or not a
        !!                                  multiple of the block size `p`).
        !!                  - `info = -4` : `kdim < 1`.
        !!                  - `info = -5` : `tol < 0`.
        !!                  - `info = -6` : empty `B`, or `size(C) /= size(B)` (block version).
        !!                  For invalid inputs (`info <= -3`) nothing is computed.
        !!
        !!  - `trans` (optional)    :   Whether \( A \) or \( A^H \) is being used.
        !!                              (default `trans=.false.`)
        !!
        !!  - `kdim` (optional)     :   Maximum number of (block) Arnoldi steps (default 100). The Krylov
        !!                              basis then holds `p*(kdim+1)` vectors, `p` being the block size.
        !!
        !!  - `X` (optional)        :   Pre-allocated workspace for the Krylov basis. Its size must be a
        !!                              multiple of the block size `p` and at least `2p`. When present,
        !!                              no Krylov vector is allocated and the number of steps is
        !!                              `min(kdim, size(X)/p - 1)`. On exit it contains the Krylov basis.
        module procedure kexpm_vec_rsp
        module procedure kexpm_mat_rsp
        module procedure kexpm_vec_rdp
        module procedure kexpm_mat_rdp
        module procedure kexpm_vec_csp
        module procedure kexpm_mat_csp
        module procedure kexpm_vec_cdp
        module procedure kexpm_mat_cdp
    end interface

    interface kexpm_substep
        !!  ### Description
        !!
        !!  Evaluates \( c = \exp(\tau A) b \) as a product of sub-steps \( \exp(\Delta t_j A) \).
        !!  Each sub-step grows a Krylov basis of at most `opts%kdim` vectors until the error estimate
        !!  is below `tol * dt / |tau|`. When the basis is exhausted, `dt` is reduced using the
        !!  Hessenberg matrix already computed (no extra matrix-vector product).
        !!
        !!  **Warm start.** When the same `kexpm_*_metadata` object is passed to successive calls with
        !!  the same `|tau|`, `kdim` and relative tolerance `tol/||b||` (within a factor 2), the first call is fully adaptive and records the number
        !!  of sub-steps. Later calls use *uniform* sub-steps and, with `opts%tune`, reduce their
        !!  number by one per call while this lowers the cost, then freeze it. The error check is
        !!  always kept: if a learned step is too large, the rest of that call falls back to the
        !!  adaptive controller (the result stays within tolerance), and a frozen count is raised
        !!  and re-tuned.
        !!  Call `meta%forget()` if the operator, the horizon or the tolerance change substantially.
        !!
        !!  ### Syntax
        !!
        !!  ```fortran
        !!      call kexpm_substep(c, A, b, tau, tol, info [, trans] [, options] [, meta])
        !!  ```
        !!
        !!  ### Arguments
        !!
        !!  - `c`       :   Output vector. `intent(out)`.
        !!  - `A`       :   Linear operator to be exponentiated. `intent(inout)`.
        !!  - `b`       :   Input vector. `intent(in)`.
        !!  - `tau`     :   Time horizon (may be negative). `intent(in)`.
        !!  - `tol`     :   Absolute tolerance on the accumulated error estimate. `intent(in)`.
        !!  - `info`    :   Total number of matrix-vector products on success, -1 if not converged,
        !!                  -4 if `opts%kdim < 1`, -5 if `tol < 0`, -7 if `opts%max_rejections < 1`
        !!                  or `opts%max_substeps < 1` (nothing is computed for invalid inputs).
        !!  - `trans` (optional)   :   Use \( A^H \) instead of \( A \) (default `.false.`).
        !!  - `options` (optional) :   `kexpm_sp_opts` / `kexpm_dp_opts`.
        !!  - `meta` (optional)    :   `kexpm_sp_metadata` / `kexpm_dp_metadata`, `intent(inout)`:
        !!                             counters and learned settings, kept between calls.
        module procedure kexpm_substep_rsp
        module procedure kexpm_substep_rdp
        module procedure kexpm_substep_csp
        module procedure kexpm_substep_cdp
    end interface

    interface krylov_exptA
        !!  ### Description
        !!
        !!  Utility function to evaluate the matrix-exponential times vector. The time horizon is split
        !!  into adaptively chosen sub-steps (step-size control based on the Krylov error estimate,
        !!  in the spirit of Sidje's Expokit `expv`), each one built from a Krylov basis of at most
        !!  30 vectors. Rejected step sizes only cost a small dense exponential, never extra
        !!  matrix-vector products.
        !!
        !!  ### Syntax
        !!
        !!  ```fortran
        !!      call krylov_exptA(vec_out, A, vec_in, tau, info, trans)
        !!  ```
        !!
        !!  ### Arguments
        !!
        !!  - `vec_out` :   Output vector.
        !!
        !!  - `A`       :   Matrix to be exponentiated.
        !!
        !!  - `vec_in`  :   Input vector.
        !!
        !!  - `tau`     :   Integration time.
        !!
        !!  - `info`    :   Information flag (total number of matrix-vector products on success).
        !!
        !!  - `trans`   :   Whether \( A \) or \( A^H \) is being used.
        module procedure krylov_exptA_rsp
        module procedure krylov_exptA_rdp
        module procedure krylov_exptA_csp
        module procedure krylov_exptA_cdp
    end interface

contains
    !-------------------------------------------------------------------
    !-----     Public driver : single vector, fixed time horizon     -----
    !-------------------------------------------------------------------

    subroutine kexpm_vec_rsp(c, A, b, tau, tol, info, trans, kdim, X)
        !! Best approximation of \( \exp(\tau \mathbf{A}) \mathbf{b} \) in the computed Krylov subspace
        implicit none(type, external)
        class(abstract_vector_rsp), intent(out) :: c
        !! Output vector.
        class(abstract_linop_rsp), intent(inout) :: A
        !! Linear operator to be exponentiated.
        class(abstract_vector_rsp), intent(in) :: b
        !! Input vector on which to apply \( \exp(\tau \mathbf{A}) \).
        real(sp), intent(in) :: tau
        !! Time horizon for the exponentiation.
        real(sp), intent(in) :: tol
        !! Solution tolerance based on error estimates.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: trans
        !! Use transpose?
        integer, optional, intent(in) :: kdim
        !! Maximum number of Arnoldi steps.
        class(abstract_vector_rsp), optional, intent(inout) :: X(:)
        !! Optional workspace for the Krylov basis (avoids allocating it at every call).

        ! ----- Internal variables -----
        character(len=*), parameter :: this_procedure = 'kexpm_vec_rsp'
        class(abstract_vector_rsp), allocatable :: Xloc(:)
        integer :: nsteps, iostat
        character(len=256) :: msg

        ! Sanity checks (no work is done for invalid inputs).
        info = 0
        if (present(X)) then
            if (size(X) < 2) info = -3
        endif
        if (info == 0 .and. present(kdim)) then
            if (kdim < 1) info = -4
        endif
        if (info == 0 .and. tol < 0.0_sp) info = -5
        if (info /= 0) then
            call log_kexpm_invalid_input(info, this_procedure)
            return
        endif

        if (present(X)) then
            nsteps = min(optval(kdim, size(X)-1), size(X)-1)
            call kexpm_vec_kernel_rsp(c, A, b, tau, tol, info, optval(trans, .false.), X(:nsteps+1))
        else
            nsteps = optval(kdim, kmax)
            ! `mold` : no copy of b into every basis vector (the kernel zeroes the basis anyway).
            allocate(Xloc(nsteps+1), mold=b, stat=iostat, errmsg=msg)
            call check_allocation(iostat, msg, this_module, this_procedure)
            call kexpm_vec_kernel_rsp(c, A, b, tau, tol, info, optval(trans, .false.), Xloc)
            ! Xloc is allocatable: automatically deallocated on return (no leak).
        endif
    end subroutine kexpm_vec_rsp

    subroutine kexpm_vec_kernel_rsp(c, A, b, tau, tol, info, transpose, X)
        !! Worker for `kexpm_vec`: the Krylov basis `X` is a plain (non-optional) dummy argument.
        implicit none(type, external)
        class(abstract_vector_rsp), intent(out) :: c
        class(abstract_linop_rsp), intent(inout) :: A
        class(abstract_vector_rsp), intent(in) :: b
        real(sp), intent(in) :: tau, tol
        integer, intent(out) :: info
        logical, intent(in) :: transpose
        class(abstract_vector_rsp), intent(inout) :: X(:)

        ! ----- Internal variables -----
        character(len=*), parameter :: this_procedure = 'kexpm_vec_rsp'
        real(sp), allocatable :: H(:, :), E(:, :)
        real(sp) :: beta, err_est
        integer :: nk, kp, ksteps, iostat
        logical :: breakdown
        character(len=256) :: msg

        nk = size(X) - 1

        beta = b%norm()
        if (beta == 0.0_sp) then
            ! Input is zero => Output is zero.
            call c%zero()
            info = 1
            call log_information('Zero input vector, zero output.', this_module, this_procedure)
            return
        endif

        allocate(H(nk+1, nk+1), E(nk+1, nk+1), source=zero_rsp, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)

        ! Initialize Krylov subspace with the normalized input vector.
        call zero_basis(X)
        call X(1)%add(b) ; call X(1)%scal(one_rsp / beta)

        ! Grow the Krylov subspace until the error estimate is below tol.
        call arnoldi_expm_rsp(A, X, H, E, beta, tau, tol, transpose, &
                                              & ksteps, kp, err_est, breakdown)

        ! Single projection back onto the full space: c = beta * X(:kp) @ E(:kp, 1).
        call project_rsp(c, X(:kp), beta*E(:kp, 1))

        if (err_est <= tol) then
            info = kp
            write(msg,'(A,I0,2(A,E9.2),A,L1)') 'Converged. kp= ', kp, ', err_est= ', err_est, &
                                             & ', tol= ', tol, ', breakdown= ', breakdown
            call log_information(msg, this_module, this_procedure)
        else
            info = -1
            write(msg,'(A,I0,2(A,E9.2))') 'Not converged. kp= ', kp, ', err_est= ', err_est, ', tol= ', tol
            call log_information(msg, this_module, this_procedure)
        endif
    end subroutine kexpm_vec_kernel_rsp

    !----------------------------------------------------------------
    !-----     Public driver : block version, fixed horizon     -----
    !----------------------------------------------------------------

    subroutine kexpm_mat_rsp(C, A, B, tau, tol, info, trans, kdim, X)
        !! Best Krylov approximation of \( \mathbf{C} = \exp(\tau \mathbf{A}) \mathbf{B} \) in the computed Krylov subspace.
        implicit none(type, external)
        class(abstract_vector_rsp), intent(out) :: C(:)
        !! Output vectors.
        class(abstract_linop_rsp), intent(inout) :: A
        !! Linear operator to be exponentiated.
        class(abstract_vector_rsp), intent(in) :: B(:)
        !! Input matrix on which to apply \( \exp(\tau \mathbf{A}) \).
        real(sp), intent(in) :: tau
        !! Time horizon for the exponentiation.
        real(sp), intent(in) :: tol
        !! Solution tolerance based on error estimates.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: trans
        !! Use transpose ?
        integer, optional, intent(in) :: kdim
        !! Maximum number of block-Arnoldi steps.
        class(abstract_vector_rsp), optional, intent(inout) :: X(:)
        !! Optional workspace for the block Krylov basis (size multiple of size(B)).

        ! ----- Internal variables -----
        character(len=*), parameter :: this_procedure = 'kexpm_mat_rsp'
        class(abstract_vector_rsp), allocatable :: Xloc(:)
        integer :: p, nsteps, iostat
        character(len=256) :: msg

        p = size(B)

        ! Sanity checks (no work is done for invalid inputs).
        info = 0
        if (p < 1 .or. size(C) /= p) then
            info = -6
        else if (present(X)) then
            if (mod(size(X), p) /= 0 .or. size(X) < 2*p) info = -3
        endif
        if (info == 0 .and. present(kdim)) then
            if (kdim < 1) info = -4
        endif
        if (info == 0 .and. tol < 0.0_sp) info = -5
        if (info /= 0) then
            call log_kexpm_invalid_input(info, this_procedure)
            return
        endif

        if (present(X)) then
            nsteps = min(optval(kdim, size(X)/p - 1), size(X)/p - 1)
            call kexpm_mat_kernel_rsp(C, A, B, tau, tol, info, optval(trans, .false.), X(:p*(nsteps+1)))
        else
            ! NOTE: kdim counts block steps. The basis holds p*(kdim+1) vectors
            !       (previously p*(p*kdim+1) because the block size was counted twice).
            nsteps = optval(kdim, kmax)
            allocate(Xloc(p*(nsteps+1)), mold=B(1), stat=iostat, errmsg=msg)
            call check_allocation(iostat, msg, this_module, this_procedure)
            call kexpm_mat_kernel_rsp(C, A, B, tau, tol, info, optval(trans, .false.), Xloc)
        endif
    end subroutine kexpm_mat_rsp

    subroutine kexpm_mat_kernel_rsp(C, A, B, tau, tol, info, transpose, X)
        !! Worker for `kexpm_mat`.
        implicit none(type, external)
        class(abstract_vector_rsp), intent(out) :: C(:)
        class(abstract_linop_rsp), intent(inout) :: A
        class(abstract_vector_rsp), intent(in) :: B(:)
        real(sp), intent(in) :: tau, tol
        integer, intent(out) :: info
        logical, intent(in) :: transpose
        class(abstract_vector_rsp), intent(inout) :: X(:)

        ! ----- Internal variables -----
        character(len=*), parameter :: this_procedure = 'kexpm_mat_rsp'
        real(sp), allocatable :: H(:, :), E(:, :), R(:, :), coef(:, :)
        integer, allocatable :: perm(:)
        real(sp) :: err_est
        integer :: p, nk, k, kpm, kp, kpp, iostat
        logical :: breakdown, partial_breakdown
        character(len=256) :: msg

        p  = size(B)
        nk = size(X)/p - 1

        allocate(R(p, p), H(p*(nk+1), p*(nk+1)), E(p*(nk+1), p*(nk+1)), &
                 source=zero_rsp, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)
        allocate(perm(p), source=0, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)

        ! Orthonormalize the input block directly inside the Krylov basis: B = X(:p) @ R.
        call zero_basis(X)
        call copy(X(:p), B)
        call qr(X(:p), R, perm, info)
        call check_info(info, 'qr_pivot', this_module, this_procedure)
        call permcols(R, invperm(perm))

        if (mnorm(R, "fro") == 0.0_sp) then
            ! Input matrix is zero.
            call zero_basis(C)
            info = p
            call log_information('Zero input block, zero output.', this_module, this_procedure)
            return
        endif

        err_est = huge(1.0_sp) ; kp = p ; kpp = p
        breakdown = .false. ; partial_breakdown = .false.

        expm_arnoldi: do k = 1, nk
            ! Set counters.
            kpm = (k-1)*p ; kp = kpm + p ; kpp = kp + p

            ! Compute the k-th step of the block Arnoldi factorization.
            call arnoldi(A, X, H, info, kstart=k, kend=k, transpose=transpose, blksize=p)
            call check_info(info, 'arnoldi', this_module, this_procedure)

            if (info == kp) then
                ! Block Arnoldi reports breakdown as soon as the new block is rank deficient.
                ! The subspace is invariant only if the whole coupling block vanishes.
                breakdown = mnorm(H(kp+1:kpp, kpm+1:kp), "fro") < atol_sp
                ! Otherwise, the QR completed the block with random vectors that are not
                ! orthogonal to the previous basis: we cannot safely continue (no deflation).
                partial_breakdown = .not. breakdown
            endif

            if (breakdown) then
                ! Exact invariant subspace: do not consider the extended matrix.
                kpp = kp
                call small_expm_rsp(E, H, kpp, kpp, tau)
                err_est = 0.0_sp
                exit expm_arnoldi
            endif

            ! Dense exponential of the extended block-Hessenberg matrix.
            call small_expm_rsp(E, H, kpp, kp, tau)

            ! Cheap (conservative) error estimate.
            err_est = mnorm(matmul(E(kp+1:kpp, :p), R), "fro")

            if (err_est <= tol .or. partial_breakdown) exit expm_arnoldi
        enddo expm_arnoldi

        ! Single projection back to the full space: C = X(:kpp) @ (E(:kpp, :p) @ R).
        coef = matmul(E(:kpp, :p), R)
        call project_block_rsp(C, X(:kpp), coef)

        if (err_est <= tol) then
            info = kpp
            write(msg,'(A,I0,2(A,E9.2),A,L1)') 'Converged. kp= ', kpp, ', err_est= ', err_est, &
                                             & ', tol= ', tol, ', breakdown= ', breakdown
            call log_information(msg, this_module, this_procedure)
        else
            info = -1
            if (partial_breakdown) then
                write(msg,'(A,I0,2(A,E9.2))') 'Not converged (partial block breakdown, deflation not supported). kp= ', &
                                            & kpp, ', err_est= ', err_est, ', tol= ', tol
            else
                write(msg,'(A,I0,2(A,E9.2))') 'Not converged. kp= ', kpp, ', err_est= ', err_est, ', tol= ', tol
            endif
            call log_information(msg, this_module, this_procedure)
        endif
    end subroutine kexpm_mat_kernel_rsp

    !--------------------------------------------------------------------
    !-----     Public driver : adaptive sub-stepping (Expokit-like)     -----
    !--------------------------------------------------------------------

    subroutine krylov_exptA_rsp(vec_out, A, vec_in, tau, info, trans)
        !! Wrapper for the Krylov-based evaluation of the action of the matrix exponential operator
        !! on a vector that conforms to the `abstract_exptA_rsp` interface.
        !! It uses `kexpm_substep` with default options and no warm start; call `kexpm_substep`
        !! directly with a persistent metadata object to re-use learned settings.
        implicit none(type, external)
        class(abstract_vector_rsp), intent(out) :: vec_out
        !! Solution vector.
        class(abstract_linop_rsp), intent(inout) :: A
        !! Linear operator to be exponentiated.
        class(abstract_vector_rsp), intent(in) :: vec_in
        !! Input vector to be multiplied by \( \exp(\tau \mathbf{A}) \).
        real(sp), intent(in) :: tau
        !! Time horizon for the exponentiation.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: trans
        !! Use adjoint ?

        ! ----- Internal variables -----
        character(len=*), parameter :: this_procedure = 'krylov_exptA_rsp'
        real(sp) :: tol

        ! Tolerance relative to the norm of the input (identical to the former atol for unit vectors).
        tol = atol_sp * vec_in%norm()

        call kexpm_substep_rsp(vec_out, A, vec_in, tau, tol, info, trans=trans)
        call check_info(info, 'kexpm', this_module, this_procedure)
    end subroutine krylov_exptA_rsp

    subroutine kexpm_substep_rsp(c, A, b, tau, tol, info, trans, options, meta)
        !! Sub-stepped Krylov approximation of \( \exp(\tau \mathbf{A}) \mathbf{b} \) with optional
        !! warm start (see the `kexpm_substep` interface).
        implicit none(type, external)
        class(abstract_vector_rsp), intent(out) :: c
        !! Output vector.
        class(abstract_linop_rsp), intent(inout) :: A
        !! Linear operator to be exponentiated.
        class(abstract_vector_rsp), intent(in) :: b
        !! Input vector.
        real(sp), intent(in) :: tau
        !! Time horizon (may be negative).
        real(sp), intent(in) :: tol
        !! Absolute tolerance on the accumulated error estimate.
        integer, intent(out) :: info
        !! Total number of matvecs on success, -1 otherwise.
        logical, optional, intent(in) :: trans
        !! Use adjoint?
        class(abstract_opts), optional, intent(in) :: options
        !! `kexpm_sp_opts`.
        class(abstract_metadata), optional, intent(inout) :: meta
        !! `kexpm_sp_metadata` (counters and learned settings, kept between calls).

        ! ----- Internal variables -----
        character(len=*), parameter :: this_procedure = 'kexpm_substep_rsp'
        real(sp), parameter :: safety = 0.9_sp, fac_min = 0.2_sp, fac_max = 5.0_sp
        type(kexpm_sp_opts) :: opts
        type(kexpm_sp_metadata) :: lmeta
        class(abstract_vector_rsp), allocatable :: X(:)
        real(sp), allocatable :: H(:, :), E(:, :)
        real(sp) :: t_end, sgn, t, dt, beta, err, tol_loc, fac, floor_tol, tol_rel
        integer :: m, ksteps, kp, nsub, nmatvec, nrej, nrej_tot, n_uniform, iostat
        logical :: transpose, breakdown, last, uniform, used_uniform, fallback
        character(len=256) :: msg

        ! ----- Options and metadata -----
        opts = kexpm_sp_opts()
        if (present(options)) then
            select type (options)
            type is (kexpm_sp_opts)
                opts = options
            class default
                call type_error('options', 'kexpm_sp_opts', 'IN', this_module, this_procedure)
            end select
        endif
        if (present(meta)) then
            select type (meta)
            type is (kexpm_sp_metadata)
                lmeta = meta
            class default
                call type_error('meta', 'kexpm_sp_metadata', 'INOUT', this_module, this_procedure)
            end select
        endif

        ! Sanity checks (no work is done and the metadata is left untouched for invalid inputs).
        info = 0
        if (opts%kdim < 1) then
            info = -4
        else if (tol < 0.0_sp) then
            info = -5
        else if (opts%max_rejections < 1 .or. opts%max_substeps < 1) then
            info = -7
        endif
        if (info /= 0) then
            call log_kexpm_invalid_input(info, this_procedure)
            return
        endif

        transpose = optval(trans, .false.)
        m = opts%kdim
        t_end = abs(tau) ; sgn = sign(1.0_sp, tau)
        lmeta%n_calls = lmeta%n_calls + 1

        ! Learned settings are only valid for the same horizon, tolerance and Krylov dimension.
        ! The tolerance is compared relative to ||b|| (within a factor 2), so that calls with
        ! tol = rtol*||b|| and different right-hand sides share the learned settings.
        beta = b%norm()
        tol_rel = merge(tol/beta, 0.0_sp, beta > 0.0_sp)
        if (lmeta%learned) then
            if (abs(t_end - lmeta%tau_learned) > 10*epsilon(1.0_sp)*t_end &
                & .or. tol_rel > 2*lmeta%tol_learned .or. 2*tol_rel < lmeta%tol_learned &
                & .or. m /= lmeta%kdim_learned) call lmeta%forget()
        endif
        used_uniform = present(meta) .and. opts%warm_start .and. lmeta%learned .and. lmeta%nsub > 0
        n_uniform = merge(lmeta%nsub, 0, used_uniform)
        uniform = used_uniform ; fallback = .false.

        ! c holds the running solution w(t).
        call c%zero() ; call c%add(b)
        beta = c%norm()
        t = 0.0_sp ; nsub = 0 ; nmatvec = 0 ; nrej_tot = 0 ; info = 0

        if (beta /= 0.0_sp .and. t_end /= 0.0_sp) then
            ! Workspace (allocated once, reused by every sub-step).
            allocate(X(m+1), mold=b, stat=iostat, errmsg=msg)
            call check_allocation(iostat, msg, this_module, this_procedure)
            allocate(H(m+1, m+1), E(m+1, m+1), source=zero_rsp, stat=iostat, errmsg=msg)
            call check_allocation(iostat, msg, this_module, this_procedure)

            dt = merge(t_end/max(n_uniform, 1), t_end, uniform)

            substeps: do while (t < t_end)
                nsub = nsub + 1
                if (nsub > opts%max_substeps) then
                    info = -1 ; exit substeps
                endif

                ! Last step (with a small slack so that round-off does not create a tiny extra step).
                last = (dt >= (t_end - t) - 10*epsilon(1.0_sp)*t_end) ; if (last) dt = t_end - t

                beta = c%norm()
                if (beta == 0.0_sp) exit substeps   ! Solution is exactly zero from now on.

                ! Roundoff floor on the local tolerance to avoid endless rejections.
                floor_tol = 2*epsilon(1.0_sp)*beta
                tol_loc = max(tol*dt/t_end, floor_tol)

                ! Fresh Krylov basis started from the current solution.
                call zero_basis(X) ; H = zero_rsp
                call X(1)%add(c) ; call X(1)%scal(one_rsp / beta)

                call arnoldi_expm_rsp(A, X, H, E, beta, sgn*dt, tol_loc, transpose, &
                                                      & ksteps, kp, err, breakdown)
                nmatvec = nmatvec + ksteps

                if (breakdown) then
                    ! Invariant subspace: exact for any dt, finish the horizon in one go.
                    dt = t_end - t ; last = .true.
                    call small_expm_rsp(E, H, kp, kp, sgn*dt)
                    err = 0.0_sp
                else if (err > tol_loc) then
                    ! Krylov basis exhausted. Safety net: a learned uniform step that is too large
                    ! switches the rest of this call to the adaptive controller.
                    if (uniform) then
                        uniform = .false. ; fallback = .true.
                    endif
                    ! Shrink dt reusing H (only small dense exponentials).
                    last = .false.
                    do nrej = 1, opts%max_rejections
                        fac = min(safety, max(0.1_sp, safety*(tol_loc/err)**(1.0_sp/m)))
                        dt  = fac*dt
                        tol_loc = max(tol*dt/t_end, floor_tol)
                        call small_expm_rsp(E, H, kp, ksteps, sgn*dt)
                        err = beta*abs(E(kp, 1))
                        nrej_tot = nrej_tot + 1
                        if (err <= tol_loc) exit
                    enddo
                    if (err > tol_loc) then
                        info = -1 ; exit substeps
                    endif
                endif

                ! Advance: w(t+dt) = beta * X(:kp) @ E(:kp, 1).
                call project_rsp(c, X(:kp), beta*E(:kp, 1))
                if (last) then
                    t = t_end
                else
                    t = t + dt
                endif

                ! Next step size.
                if (uniform) then
                    dt = t_end/n_uniform
                else
                    ! Error controller (local error ~ dt^ksteps).
                    if (err > 0.0_sp) then
                        fac = safety*(tol_loc/err)**(1.0_sp/max(ksteps, 1))
                    else
                        fac = fac_max
                    endif
                    dt = dt*min(fac_max, max(fac_min, fac))
                endif
            enddo substeps
        endif

        if (info == 0) info = max(nmatvec, 1)

        ! ----- Update metadata (counters, learned settings) -----
        lmeta%n_matvec     = lmeta%n_matvec + nmatvec
        lmeta%n_substeps   = lmeta%n_substeps + nsub
        lmeta%n_rejections = lmeta%n_rejections + nrej_tot
        if (fallback) lmeta%n_fallbacks = lmeta%n_fallbacks + 1
        lmeta%info = info ; lmeta%converged = (info > 0)
        if (present(meta) .and. opts%warm_start .and. info > 0 .and. nsub > 0) then
            call learn_kexpm_sp(lmeta, used_uniform, n_uniform, nsub, nmatvec, fallback, &
                                      & opts%tune, t_end, tol_rel, m)
        endif

        if (info > 0) then
            write(msg,'(A,I0,A,I0,A,I0,A,L1,A,L1)') 'Converged. substeps= ', nsub, ', matvecs= ', nmatvec, &
                & ', rejections= ', nrej_tot, ', warm start= ', used_uniform, ', fallback= ', fallback
        else
            write(msg,'(A,I0,A,E9.2,A,E9.2)') 'Not converged. substeps= ', nsub, ', t= ', sgn*t, ', tau= ', tau
        endif
        call log_information(msg, this_module, this_procedure)

        if (opts%if_print_metadata) call lmeta%print()
        if (present(meta)) then
            select type (meta)
            type is (kexpm_sp_metadata)
                meta = lmeta
            end select
        endif
    end subroutine kexpm_substep_rsp

    !----------------------------------------
    !-----     Private helper kernels     -----
    !----------------------------------------

    subroutine arnoldi_expm_rsp(A, X, H, E, beta, tau, tol, transpose, ksteps, kp, err_est, breakdown)
        !! Grows the Arnoldi basis `X` (with `X(1)` already normalized) one vector at a time and
        !! evaluates the exponential of the extended Hessenberg matrix after each step, until
        !! `beta*|E(kp,1)| <= tol`, breakdown, or `size(X)-1` steps.
        !! On exit `E(:kp, :kp)` holds the dense exponential used for the projection.
        implicit none(type, external)
        class(abstract_linop_rsp), intent(inout) :: A
        class(abstract_vector_rsp), intent(inout) :: X(:)
        real(sp), intent(inout) :: H(:, :), E(:, :)
        real(sp), intent(in) :: beta, tau, tol
        logical, intent(in) :: transpose
        integer, intent(out) :: ksteps, kp
        real(sp), intent(out) :: err_est
        logical, intent(out) :: breakdown

        character(len=*), parameter :: this_procedure = 'arnoldi_expm_rsp'
        integer :: k, info

        err_est = huge(1.0_sp) ; breakdown = .false. ; ksteps = 0 ; kp = 1

        do k = 1, size(X) - 1
            call arnoldi(A, X, H, info, kstart=k, kend=k, transpose=transpose)
            call check_info(info, 'arnoldi', this_module, this_procedure)
            ksteps = k
            breakdown = (info == k)

            if (breakdown) then
                ! Invariant subspace: the approximation is exact, drop the extended row.
                kp = k
                call small_expm_rsp(E, H, kp, kp, tau)
                err_est = 0.0_sp
                return
            endif

            ! Exponential of the (k+1)x(k+1) extended Hessenberg matrix [H_k 0; h e_k^T 0].
            kp = k + 1
            call small_expm_rsp(E, H, kp, k, tau)

            ! Cheap error estimate (magnitude of the last included correction, conservative).
            err_est = beta*abs(E(kp, 1))
            if (err_est <= tol) return
        enddo
    end subroutine arnoldi_expm_rsp

    subroutine small_expm_rsp(E, H, n, ncols, tau)
        !! E(:n, :n) = exp(tau * Hbar) where Hbar holds the first `ncols` columns of H(:n, :)
        !! and zeros elsewhere. Only the leading n x n block of E is touched (no full reset),
        !! and stale columns of H (e.g. from a previous sub-step) are never used.
        implicit none(type, external)
        real(sp), intent(inout) :: E(:, :)
        real(sp), intent(in) :: H(:, :)
        integer, intent(in) :: n, ncols
        real(sp), intent(in) :: tau

        E(:n, :ncols) = tau*H(:n, :ncols)
        if (ncols < n) E(:n, ncols+1:n) = zero_rsp
        call matrix_exp(E(:n, :n))
    end subroutine small_expm_rsp

    subroutine project_rsp(c, X, v)
        !! c = X @ v, written directly into the (non-allocatable) output vector.
        implicit none(type, external)
        class(abstract_vector_rsp), intent(inout) :: c
        class(abstract_vector_rsp), intent(in) :: X(:)
        real(sp), intent(in) :: v(:)
        integer :: i
        call c%zero()
        do i = 1, size(X)
            call c%axpby(v(i), X(i), one_rsp)
        enddo
    end subroutine project_rsp

    subroutine project_block_rsp(C, X, V)
        !! C = X @ V for a block of output vectors.
        implicit none(type, external)
        class(abstract_vector_rsp), intent(inout) :: C(:)
        class(abstract_vector_rsp), intent(in) :: X(:)
        real(sp), intent(in) :: V(:, :)
        integer :: j
        do j = 1, size(C)
            call project_rsp(C(j), X, V(:, j))
        enddo
    end subroutine project_block_rsp

    !-------------------------------------------------------------------
    !-----     Public driver : single vector, fixed time horizon     -----
    !-------------------------------------------------------------------

    subroutine kexpm_vec_rdp(c, A, b, tau, tol, info, trans, kdim, X)
        !! Best approximation of \( \exp(\tau \mathbf{A}) \mathbf{b} \) in the computed Krylov subspace
        implicit none(type, external)
        class(abstract_vector_rdp), intent(out) :: c
        !! Output vector.
        class(abstract_linop_rdp), intent(inout) :: A
        !! Linear operator to be exponentiated.
        class(abstract_vector_rdp), intent(in) :: b
        !! Input vector on which to apply \( \exp(\tau \mathbf{A}) \).
        real(dp), intent(in) :: tau
        !! Time horizon for the exponentiation.
        real(dp), intent(in) :: tol
        !! Solution tolerance based on error estimates.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: trans
        !! Use transpose?
        integer, optional, intent(in) :: kdim
        !! Maximum number of Arnoldi steps.
        class(abstract_vector_rdp), optional, intent(inout) :: X(:)
        !! Optional workspace for the Krylov basis (avoids allocating it at every call).

        ! ----- Internal variables -----
        character(len=*), parameter :: this_procedure = 'kexpm_vec_rdp'
        class(abstract_vector_rdp), allocatable :: Xloc(:)
        integer :: nsteps, iostat
        character(len=256) :: msg

        ! Sanity checks (no work is done for invalid inputs).
        info = 0
        if (present(X)) then
            if (size(X) < 2) info = -3
        endif
        if (info == 0 .and. present(kdim)) then
            if (kdim < 1) info = -4
        endif
        if (info == 0 .and. tol < 0.0_dp) info = -5
        if (info /= 0) then
            call log_kexpm_invalid_input(info, this_procedure)
            return
        endif

        if (present(X)) then
            nsteps = min(optval(kdim, size(X)-1), size(X)-1)
            call kexpm_vec_kernel_rdp(c, A, b, tau, tol, info, optval(trans, .false.), X(:nsteps+1))
        else
            nsteps = optval(kdim, kmax)
            ! `mold` : no copy of b into every basis vector (the kernel zeroes the basis anyway).
            allocate(Xloc(nsteps+1), mold=b, stat=iostat, errmsg=msg)
            call check_allocation(iostat, msg, this_module, this_procedure)
            call kexpm_vec_kernel_rdp(c, A, b, tau, tol, info, optval(trans, .false.), Xloc)
            ! Xloc is allocatable: automatically deallocated on return (no leak).
        endif
    end subroutine kexpm_vec_rdp

    subroutine kexpm_vec_kernel_rdp(c, A, b, tau, tol, info, transpose, X)
        !! Worker for `kexpm_vec`: the Krylov basis `X` is a plain (non-optional) dummy argument.
        implicit none(type, external)
        class(abstract_vector_rdp), intent(out) :: c
        class(abstract_linop_rdp), intent(inout) :: A
        class(abstract_vector_rdp), intent(in) :: b
        real(dp), intent(in) :: tau, tol
        integer, intent(out) :: info
        logical, intent(in) :: transpose
        class(abstract_vector_rdp), intent(inout) :: X(:)

        ! ----- Internal variables -----
        character(len=*), parameter :: this_procedure = 'kexpm_vec_rdp'
        real(dp), allocatable :: H(:, :), E(:, :)
        real(dp) :: beta, err_est
        integer :: nk, kp, ksteps, iostat
        logical :: breakdown
        character(len=256) :: msg

        nk = size(X) - 1

        beta = b%norm()
        if (beta == 0.0_dp) then
            ! Input is zero => Output is zero.
            call c%zero()
            info = 1
            call log_information('Zero input vector, zero output.', this_module, this_procedure)
            return
        endif

        allocate(H(nk+1, nk+1), E(nk+1, nk+1), source=zero_rdp, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)

        ! Initialize Krylov subspace with the normalized input vector.
        call zero_basis(X)
        call X(1)%add(b) ; call X(1)%scal(one_rdp / beta)

        ! Grow the Krylov subspace until the error estimate is below tol.
        call arnoldi_expm_rdp(A, X, H, E, beta, tau, tol, transpose, &
                                              & ksteps, kp, err_est, breakdown)

        ! Single projection back onto the full space: c = beta * X(:kp) @ E(:kp, 1).
        call project_rdp(c, X(:kp), beta*E(:kp, 1))

        if (err_est <= tol) then
            info = kp
            write(msg,'(A,I0,2(A,E9.2),A,L1)') 'Converged. kp= ', kp, ', err_est= ', err_est, &
                                             & ', tol= ', tol, ', breakdown= ', breakdown
            call log_information(msg, this_module, this_procedure)
        else
            info = -1
            write(msg,'(A,I0,2(A,E9.2))') 'Not converged. kp= ', kp, ', err_est= ', err_est, ', tol= ', tol
            call log_information(msg, this_module, this_procedure)
        endif
    end subroutine kexpm_vec_kernel_rdp

    !----------------------------------------------------------------
    !-----     Public driver : block version, fixed horizon     -----
    !----------------------------------------------------------------

    subroutine kexpm_mat_rdp(C, A, B, tau, tol, info, trans, kdim, X)
        !! Best Krylov approximation of \( \mathbf{C} = \exp(\tau \mathbf{A}) \mathbf{B} \) in the computed Krylov subspace.
        implicit none(type, external)
        class(abstract_vector_rdp), intent(out) :: C(:)
        !! Output vectors.
        class(abstract_linop_rdp), intent(inout) :: A
        !! Linear operator to be exponentiated.
        class(abstract_vector_rdp), intent(in) :: B(:)
        !! Input matrix on which to apply \( \exp(\tau \mathbf{A}) \).
        real(dp), intent(in) :: tau
        !! Time horizon for the exponentiation.
        real(dp), intent(in) :: tol
        !! Solution tolerance based on error estimates.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: trans
        !! Use transpose ?
        integer, optional, intent(in) :: kdim
        !! Maximum number of block-Arnoldi steps.
        class(abstract_vector_rdp), optional, intent(inout) :: X(:)
        !! Optional workspace for the block Krylov basis (size multiple of size(B)).

        ! ----- Internal variables -----
        character(len=*), parameter :: this_procedure = 'kexpm_mat_rdp'
        class(abstract_vector_rdp), allocatable :: Xloc(:)
        integer :: p, nsteps, iostat
        character(len=256) :: msg

        p = size(B)

        ! Sanity checks (no work is done for invalid inputs).
        info = 0
        if (p < 1 .or. size(C) /= p) then
            info = -6
        else if (present(X)) then
            if (mod(size(X), p) /= 0 .or. size(X) < 2*p) info = -3
        endif
        if (info == 0 .and. present(kdim)) then
            if (kdim < 1) info = -4
        endif
        if (info == 0 .and. tol < 0.0_dp) info = -5
        if (info /= 0) then
            call log_kexpm_invalid_input(info, this_procedure)
            return
        endif

        if (present(X)) then
            nsteps = min(optval(kdim, size(X)/p - 1), size(X)/p - 1)
            call kexpm_mat_kernel_rdp(C, A, B, tau, tol, info, optval(trans, .false.), X(:p*(nsteps+1)))
        else
            ! NOTE: kdim counts block steps. The basis holds p*(kdim+1) vectors
            !       (previously p*(p*kdim+1) because the block size was counted twice).
            nsteps = optval(kdim, kmax)
            allocate(Xloc(p*(nsteps+1)), mold=B(1), stat=iostat, errmsg=msg)
            call check_allocation(iostat, msg, this_module, this_procedure)
            call kexpm_mat_kernel_rdp(C, A, B, tau, tol, info, optval(trans, .false.), Xloc)
        endif
    end subroutine kexpm_mat_rdp

    subroutine kexpm_mat_kernel_rdp(C, A, B, tau, tol, info, transpose, X)
        !! Worker for `kexpm_mat`.
        implicit none(type, external)
        class(abstract_vector_rdp), intent(out) :: C(:)
        class(abstract_linop_rdp), intent(inout) :: A
        class(abstract_vector_rdp), intent(in) :: B(:)
        real(dp), intent(in) :: tau, tol
        integer, intent(out) :: info
        logical, intent(in) :: transpose
        class(abstract_vector_rdp), intent(inout) :: X(:)

        ! ----- Internal variables -----
        character(len=*), parameter :: this_procedure = 'kexpm_mat_rdp'
        real(dp), allocatable :: H(:, :), E(:, :), R(:, :), coef(:, :)
        integer, allocatable :: perm(:)
        real(dp) :: err_est
        integer :: p, nk, k, kpm, kp, kpp, iostat
        logical :: breakdown, partial_breakdown
        character(len=256) :: msg

        p  = size(B)
        nk = size(X)/p - 1

        allocate(R(p, p), H(p*(nk+1), p*(nk+1)), E(p*(nk+1), p*(nk+1)), &
                 source=zero_rdp, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)
        allocate(perm(p), source=0, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)

        ! Orthonormalize the input block directly inside the Krylov basis: B = X(:p) @ R.
        call zero_basis(X)
        call copy(X(:p), B)
        call qr(X(:p), R, perm, info)
        call check_info(info, 'qr_pivot', this_module, this_procedure)
        call permcols(R, invperm(perm))

        if (mnorm(R, "fro") == 0.0_dp) then
            ! Input matrix is zero.
            call zero_basis(C)
            info = p
            call log_information('Zero input block, zero output.', this_module, this_procedure)
            return
        endif

        err_est = huge(1.0_dp) ; kp = p ; kpp = p
        breakdown = .false. ; partial_breakdown = .false.

        expm_arnoldi: do k = 1, nk
            ! Set counters.
            kpm = (k-1)*p ; kp = kpm + p ; kpp = kp + p

            ! Compute the k-th step of the block Arnoldi factorization.
            call arnoldi(A, X, H, info, kstart=k, kend=k, transpose=transpose, blksize=p)
            call check_info(info, 'arnoldi', this_module, this_procedure)

            if (info == kp) then
                ! Block Arnoldi reports breakdown as soon as the new block is rank deficient.
                ! The subspace is invariant only if the whole coupling block vanishes.
                breakdown = mnorm(H(kp+1:kpp, kpm+1:kp), "fro") < atol_dp
                ! Otherwise, the QR completed the block with random vectors that are not
                ! orthogonal to the previous basis: we cannot safely continue (no deflation).
                partial_breakdown = .not. breakdown
            endif

            if (breakdown) then
                ! Exact invariant subspace: do not consider the extended matrix.
                kpp = kp
                call small_expm_rdp(E, H, kpp, kpp, tau)
                err_est = 0.0_dp
                exit expm_arnoldi
            endif

            ! Dense exponential of the extended block-Hessenberg matrix.
            call small_expm_rdp(E, H, kpp, kp, tau)

            ! Cheap (conservative) error estimate.
            err_est = mnorm(matmul(E(kp+1:kpp, :p), R), "fro")

            if (err_est <= tol .or. partial_breakdown) exit expm_arnoldi
        enddo expm_arnoldi

        ! Single projection back to the full space: C = X(:kpp) @ (E(:kpp, :p) @ R).
        coef = matmul(E(:kpp, :p), R)
        call project_block_rdp(C, X(:kpp), coef)

        if (err_est <= tol) then
            info = kpp
            write(msg,'(A,I0,2(A,E9.2),A,L1)') 'Converged. kp= ', kpp, ', err_est= ', err_est, &
                                             & ', tol= ', tol, ', breakdown= ', breakdown
            call log_information(msg, this_module, this_procedure)
        else
            info = -1
            if (partial_breakdown) then
                write(msg,'(A,I0,2(A,E9.2))') 'Not converged (partial block breakdown, deflation not supported). kp= ', &
                                            & kpp, ', err_est= ', err_est, ', tol= ', tol
            else
                write(msg,'(A,I0,2(A,E9.2))') 'Not converged. kp= ', kpp, ', err_est= ', err_est, ', tol= ', tol
            endif
            call log_information(msg, this_module, this_procedure)
        endif
    end subroutine kexpm_mat_kernel_rdp

    !--------------------------------------------------------------------
    !-----     Public driver : adaptive sub-stepping (Expokit-like)     -----
    !--------------------------------------------------------------------

    subroutine krylov_exptA_rdp(vec_out, A, vec_in, tau, info, trans)
        !! Wrapper for the Krylov-based evaluation of the action of the matrix exponential operator
        !! on a vector that conforms to the `abstract_exptA_rdp` interface.
        !! It uses `kexpm_substep` with default options and no warm start; call `kexpm_substep`
        !! directly with a persistent metadata object to re-use learned settings.
        implicit none(type, external)
        class(abstract_vector_rdp), intent(out) :: vec_out
        !! Solution vector.
        class(abstract_linop_rdp), intent(inout) :: A
        !! Linear operator to be exponentiated.
        class(abstract_vector_rdp), intent(in) :: vec_in
        !! Input vector to be multiplied by \( \exp(\tau \mathbf{A}) \).
        real(dp), intent(in) :: tau
        !! Time horizon for the exponentiation.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: trans
        !! Use adjoint ?

        ! ----- Internal variables -----
        character(len=*), parameter :: this_procedure = 'krylov_exptA_rdp'
        real(dp) :: tol

        ! Tolerance relative to the norm of the input (identical to the former atol for unit vectors).
        tol = atol_dp * vec_in%norm()

        call kexpm_substep_rdp(vec_out, A, vec_in, tau, tol, info, trans=trans)
        call check_info(info, 'kexpm', this_module, this_procedure)
    end subroutine krylov_exptA_rdp

    subroutine kexpm_substep_rdp(c, A, b, tau, tol, info, trans, options, meta)
        !! Sub-stepped Krylov approximation of \( \exp(\tau \mathbf{A}) \mathbf{b} \) with optional
        !! warm start (see the `kexpm_substep` interface).
        implicit none(type, external)
        class(abstract_vector_rdp), intent(out) :: c
        !! Output vector.
        class(abstract_linop_rdp), intent(inout) :: A
        !! Linear operator to be exponentiated.
        class(abstract_vector_rdp), intent(in) :: b
        !! Input vector.
        real(dp), intent(in) :: tau
        !! Time horizon (may be negative).
        real(dp), intent(in) :: tol
        !! Absolute tolerance on the accumulated error estimate.
        integer, intent(out) :: info
        !! Total number of matvecs on success, -1 otherwise.
        logical, optional, intent(in) :: trans
        !! Use adjoint?
        class(abstract_opts), optional, intent(in) :: options
        !! `kexpm_dp_opts`.
        class(abstract_metadata), optional, intent(inout) :: meta
        !! `kexpm_dp_metadata` (counters and learned settings, kept between calls).

        ! ----- Internal variables -----
        character(len=*), parameter :: this_procedure = 'kexpm_substep_rdp'
        real(dp), parameter :: safety = 0.9_dp, fac_min = 0.2_dp, fac_max = 5.0_dp
        type(kexpm_dp_opts) :: opts
        type(kexpm_dp_metadata) :: lmeta
        class(abstract_vector_rdp), allocatable :: X(:)
        real(dp), allocatable :: H(:, :), E(:, :)
        real(dp) :: t_end, sgn, t, dt, beta, err, tol_loc, fac, floor_tol, tol_rel
        integer :: m, ksteps, kp, nsub, nmatvec, nrej, nrej_tot, n_uniform, iostat
        logical :: transpose, breakdown, last, uniform, used_uniform, fallback
        character(len=256) :: msg

        ! ----- Options and metadata -----
        opts = kexpm_dp_opts()
        if (present(options)) then
            select type (options)
            type is (kexpm_dp_opts)
                opts = options
            class default
                call type_error('options', 'kexpm_dp_opts', 'IN', this_module, this_procedure)
            end select
        endif
        if (present(meta)) then
            select type (meta)
            type is (kexpm_dp_metadata)
                lmeta = meta
            class default
                call type_error('meta', 'kexpm_dp_metadata', 'INOUT', this_module, this_procedure)
            end select
        endif

        ! Sanity checks (no work is done and the metadata is left untouched for invalid inputs).
        info = 0
        if (opts%kdim < 1) then
            info = -4
        else if (tol < 0.0_dp) then
            info = -5
        else if (opts%max_rejections < 1 .or. opts%max_substeps < 1) then
            info = -7
        endif
        if (info /= 0) then
            call log_kexpm_invalid_input(info, this_procedure)
            return
        endif

        transpose = optval(trans, .false.)
        m = opts%kdim
        t_end = abs(tau) ; sgn = sign(1.0_dp, tau)
        lmeta%n_calls = lmeta%n_calls + 1

        ! Learned settings are only valid for the same horizon, tolerance and Krylov dimension.
        ! The tolerance is compared relative to ||b|| (within a factor 2), so that calls with
        ! tol = rtol*||b|| and different right-hand sides share the learned settings.
        beta = b%norm()
        tol_rel = merge(tol/beta, 0.0_dp, beta > 0.0_dp)
        if (lmeta%learned) then
            if (abs(t_end - lmeta%tau_learned) > 10*epsilon(1.0_dp)*t_end &
                & .or. tol_rel > 2*lmeta%tol_learned .or. 2*tol_rel < lmeta%tol_learned &
                & .or. m /= lmeta%kdim_learned) call lmeta%forget()
        endif
        used_uniform = present(meta) .and. opts%warm_start .and. lmeta%learned .and. lmeta%nsub > 0
        n_uniform = merge(lmeta%nsub, 0, used_uniform)
        uniform = used_uniform ; fallback = .false.

        ! c holds the running solution w(t).
        call c%zero() ; call c%add(b)
        beta = c%norm()
        t = 0.0_dp ; nsub = 0 ; nmatvec = 0 ; nrej_tot = 0 ; info = 0

        if (beta /= 0.0_dp .and. t_end /= 0.0_dp) then
            ! Workspace (allocated once, reused by every sub-step).
            allocate(X(m+1), mold=b, stat=iostat, errmsg=msg)
            call check_allocation(iostat, msg, this_module, this_procedure)
            allocate(H(m+1, m+1), E(m+1, m+1), source=zero_rdp, stat=iostat, errmsg=msg)
            call check_allocation(iostat, msg, this_module, this_procedure)

            dt = merge(t_end/max(n_uniform, 1), t_end, uniform)

            substeps: do while (t < t_end)
                nsub = nsub + 1
                if (nsub > opts%max_substeps) then
                    info = -1 ; exit substeps
                endif

                ! Last step (with a small slack so that round-off does not create a tiny extra step).
                last = (dt >= (t_end - t) - 10*epsilon(1.0_dp)*t_end) ; if (last) dt = t_end - t

                beta = c%norm()
                if (beta == 0.0_dp) exit substeps   ! Solution is exactly zero from now on.

                ! Roundoff floor on the local tolerance to avoid endless rejections.
                floor_tol = 2*epsilon(1.0_dp)*beta
                tol_loc = max(tol*dt/t_end, floor_tol)

                ! Fresh Krylov basis started from the current solution.
                call zero_basis(X) ; H = zero_rdp
                call X(1)%add(c) ; call X(1)%scal(one_rdp / beta)

                call arnoldi_expm_rdp(A, X, H, E, beta, sgn*dt, tol_loc, transpose, &
                                                      & ksteps, kp, err, breakdown)
                nmatvec = nmatvec + ksteps

                if (breakdown) then
                    ! Invariant subspace: exact for any dt, finish the horizon in one go.
                    dt = t_end - t ; last = .true.
                    call small_expm_rdp(E, H, kp, kp, sgn*dt)
                    err = 0.0_dp
                else if (err > tol_loc) then
                    ! Krylov basis exhausted. Safety net: a learned uniform step that is too large
                    ! switches the rest of this call to the adaptive controller.
                    if (uniform) then
                        uniform = .false. ; fallback = .true.
                    endif
                    ! Shrink dt reusing H (only small dense exponentials).
                    last = .false.
                    do nrej = 1, opts%max_rejections
                        fac = min(safety, max(0.1_dp, safety*(tol_loc/err)**(1.0_dp/m)))
                        dt  = fac*dt
                        tol_loc = max(tol*dt/t_end, floor_tol)
                        call small_expm_rdp(E, H, kp, ksteps, sgn*dt)
                        err = beta*abs(E(kp, 1))
                        nrej_tot = nrej_tot + 1
                        if (err <= tol_loc) exit
                    enddo
                    if (err > tol_loc) then
                        info = -1 ; exit substeps
                    endif
                endif

                ! Advance: w(t+dt) = beta * X(:kp) @ E(:kp, 1).
                call project_rdp(c, X(:kp), beta*E(:kp, 1))
                if (last) then
                    t = t_end
                else
                    t = t + dt
                endif

                ! Next step size.
                if (uniform) then
                    dt = t_end/n_uniform
                else
                    ! Error controller (local error ~ dt^ksteps).
                    if (err > 0.0_dp) then
                        fac = safety*(tol_loc/err)**(1.0_dp/max(ksteps, 1))
                    else
                        fac = fac_max
                    endif
                    dt = dt*min(fac_max, max(fac_min, fac))
                endif
            enddo substeps
        endif

        if (info == 0) info = max(nmatvec, 1)

        ! ----- Update metadata (counters, learned settings) -----
        lmeta%n_matvec     = lmeta%n_matvec + nmatvec
        lmeta%n_substeps   = lmeta%n_substeps + nsub
        lmeta%n_rejections = lmeta%n_rejections + nrej_tot
        if (fallback) lmeta%n_fallbacks = lmeta%n_fallbacks + 1
        lmeta%info = info ; lmeta%converged = (info > 0)
        if (present(meta) .and. opts%warm_start .and. info > 0 .and. nsub > 0) then
            call learn_kexpm_dp(lmeta, used_uniform, n_uniform, nsub, nmatvec, fallback, &
                                      & opts%tune, t_end, tol_rel, m)
        endif

        if (info > 0) then
            write(msg,'(A,I0,A,I0,A,I0,A,L1,A,L1)') 'Converged. substeps= ', nsub, ', matvecs= ', nmatvec, &
                & ', rejections= ', nrej_tot, ', warm start= ', used_uniform, ', fallback= ', fallback
        else
            write(msg,'(A,I0,A,E9.2,A,E9.2)') 'Not converged. substeps= ', nsub, ', t= ', sgn*t, ', tau= ', tau
        endif
        call log_information(msg, this_module, this_procedure)

        if (opts%if_print_metadata) call lmeta%print()
        if (present(meta)) then
            select type (meta)
            type is (kexpm_dp_metadata)
                meta = lmeta
            end select
        endif
    end subroutine kexpm_substep_rdp

    !----------------------------------------
    !-----     Private helper kernels     -----
    !----------------------------------------

    subroutine arnoldi_expm_rdp(A, X, H, E, beta, tau, tol, transpose, ksteps, kp, err_est, breakdown)
        !! Grows the Arnoldi basis `X` (with `X(1)` already normalized) one vector at a time and
        !! evaluates the exponential of the extended Hessenberg matrix after each step, until
        !! `beta*|E(kp,1)| <= tol`, breakdown, or `size(X)-1` steps.
        !! On exit `E(:kp, :kp)` holds the dense exponential used for the projection.
        implicit none(type, external)
        class(abstract_linop_rdp), intent(inout) :: A
        class(abstract_vector_rdp), intent(inout) :: X(:)
        real(dp), intent(inout) :: H(:, :), E(:, :)
        real(dp), intent(in) :: beta, tau, tol
        logical, intent(in) :: transpose
        integer, intent(out) :: ksteps, kp
        real(dp), intent(out) :: err_est
        logical, intent(out) :: breakdown

        character(len=*), parameter :: this_procedure = 'arnoldi_expm_rdp'
        integer :: k, info

        err_est = huge(1.0_dp) ; breakdown = .false. ; ksteps = 0 ; kp = 1

        do k = 1, size(X) - 1
            call arnoldi(A, X, H, info, kstart=k, kend=k, transpose=transpose)
            call check_info(info, 'arnoldi', this_module, this_procedure)
            ksteps = k
            breakdown = (info == k)

            if (breakdown) then
                ! Invariant subspace: the approximation is exact, drop the extended row.
                kp = k
                call small_expm_rdp(E, H, kp, kp, tau)
                err_est = 0.0_dp
                return
            endif

            ! Exponential of the (k+1)x(k+1) extended Hessenberg matrix [H_k 0; h e_k^T 0].
            kp = k + 1
            call small_expm_rdp(E, H, kp, k, tau)

            ! Cheap error estimate (magnitude of the last included correction, conservative).
            err_est = beta*abs(E(kp, 1))
            if (err_est <= tol) return
        enddo
    end subroutine arnoldi_expm_rdp

    subroutine small_expm_rdp(E, H, n, ncols, tau)
        !! E(:n, :n) = exp(tau * Hbar) where Hbar holds the first `ncols` columns of H(:n, :)
        !! and zeros elsewhere. Only the leading n x n block of E is touched (no full reset),
        !! and stale columns of H (e.g. from a previous sub-step) are never used.
        implicit none(type, external)
        real(dp), intent(inout) :: E(:, :)
        real(dp), intent(in) :: H(:, :)
        integer, intent(in) :: n, ncols
        real(dp), intent(in) :: tau

        E(:n, :ncols) = tau*H(:n, :ncols)
        if (ncols < n) E(:n, ncols+1:n) = zero_rdp
        call matrix_exp(E(:n, :n))
    end subroutine small_expm_rdp

    subroutine project_rdp(c, X, v)
        !! c = X @ v, written directly into the (non-allocatable) output vector.
        implicit none(type, external)
        class(abstract_vector_rdp), intent(inout) :: c
        class(abstract_vector_rdp), intent(in) :: X(:)
        real(dp), intent(in) :: v(:)
        integer :: i
        call c%zero()
        do i = 1, size(X)
            call c%axpby(v(i), X(i), one_rdp)
        enddo
    end subroutine project_rdp

    subroutine project_block_rdp(C, X, V)
        !! C = X @ V for a block of output vectors.
        implicit none(type, external)
        class(abstract_vector_rdp), intent(inout) :: C(:)
        class(abstract_vector_rdp), intent(in) :: X(:)
        real(dp), intent(in) :: V(:, :)
        integer :: j
        do j = 1, size(C)
            call project_rdp(C(j), X, V(:, j))
        enddo
    end subroutine project_block_rdp

    !-------------------------------------------------------------------
    !-----     Public driver : single vector, fixed time horizon     -----
    !-------------------------------------------------------------------

    subroutine kexpm_vec_csp(c, A, b, tau, tol, info, trans, kdim, X)
        !! Best approximation of \( \exp(\tau \mathbf{A}) \mathbf{b} \) in the computed Krylov subspace
        implicit none(type, external)
        class(abstract_vector_csp), intent(out) :: c
        !! Output vector.
        class(abstract_linop_csp), intent(inout) :: A
        !! Linear operator to be exponentiated.
        class(abstract_vector_csp), intent(in) :: b
        !! Input vector on which to apply \( \exp(\tau \mathbf{A}) \).
        real(sp), intent(in) :: tau
        !! Time horizon for the exponentiation.
        real(sp), intent(in) :: tol
        !! Solution tolerance based on error estimates.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: trans
        !! Use transpose?
        integer, optional, intent(in) :: kdim
        !! Maximum number of Arnoldi steps.
        class(abstract_vector_csp), optional, intent(inout) :: X(:)
        !! Optional workspace for the Krylov basis (avoids allocating it at every call).

        ! ----- Internal variables -----
        character(len=*), parameter :: this_procedure = 'kexpm_vec_csp'
        class(abstract_vector_csp), allocatable :: Xloc(:)
        integer :: nsteps, iostat
        character(len=256) :: msg

        ! Sanity checks (no work is done for invalid inputs).
        info = 0
        if (present(X)) then
            if (size(X) < 2) info = -3
        endif
        if (info == 0 .and. present(kdim)) then
            if (kdim < 1) info = -4
        endif
        if (info == 0 .and. tol < 0.0_sp) info = -5
        if (info /= 0) then
            call log_kexpm_invalid_input(info, this_procedure)
            return
        endif

        if (present(X)) then
            nsteps = min(optval(kdim, size(X)-1), size(X)-1)
            call kexpm_vec_kernel_csp(c, A, b, tau, tol, info, optval(trans, .false.), X(:nsteps+1))
        else
            nsteps = optval(kdim, kmax)
            ! `mold` : no copy of b into every basis vector (the kernel zeroes the basis anyway).
            allocate(Xloc(nsteps+1), mold=b, stat=iostat, errmsg=msg)
            call check_allocation(iostat, msg, this_module, this_procedure)
            call kexpm_vec_kernel_csp(c, A, b, tau, tol, info, optval(trans, .false.), Xloc)
            ! Xloc is allocatable: automatically deallocated on return (no leak).
        endif
    end subroutine kexpm_vec_csp

    subroutine kexpm_vec_kernel_csp(c, A, b, tau, tol, info, transpose, X)
        !! Worker for `kexpm_vec`: the Krylov basis `X` is a plain (non-optional) dummy argument.
        implicit none(type, external)
        class(abstract_vector_csp), intent(out) :: c
        class(abstract_linop_csp), intent(inout) :: A
        class(abstract_vector_csp), intent(in) :: b
        real(sp), intent(in) :: tau, tol
        integer, intent(out) :: info
        logical, intent(in) :: transpose
        class(abstract_vector_csp), intent(inout) :: X(:)

        ! ----- Internal variables -----
        character(len=*), parameter :: this_procedure = 'kexpm_vec_csp'
        complex(sp), allocatable :: H(:, :), E(:, :)
        real(sp) :: beta, err_est
        integer :: nk, kp, ksteps, iostat
        logical :: breakdown
        character(len=256) :: msg

        nk = size(X) - 1

        beta = b%norm()
        if (beta == 0.0_sp) then
            ! Input is zero => Output is zero.
            call c%zero()
            info = 1
            call log_information('Zero input vector, zero output.', this_module, this_procedure)
            return
        endif

        allocate(H(nk+1, nk+1), E(nk+1, nk+1), source=zero_csp, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)

        ! Initialize Krylov subspace with the normalized input vector.
        call zero_basis(X)
        call X(1)%add(b) ; call X(1)%scal(one_csp / beta)

        ! Grow the Krylov subspace until the error estimate is below tol.
        call arnoldi_expm_csp(A, X, H, E, beta, tau, tol, transpose, &
                                              & ksteps, kp, err_est, breakdown)

        ! Single projection back onto the full space: c = beta * X(:kp) @ E(:kp, 1).
        call project_csp(c, X(:kp), beta*E(:kp, 1))

        if (err_est <= tol) then
            info = kp
            write(msg,'(A,I0,2(A,E9.2),A,L1)') 'Converged. kp= ', kp, ', err_est= ', err_est, &
                                             & ', tol= ', tol, ', breakdown= ', breakdown
            call log_information(msg, this_module, this_procedure)
        else
            info = -1
            write(msg,'(A,I0,2(A,E9.2))') 'Not converged. kp= ', kp, ', err_est= ', err_est, ', tol= ', tol
            call log_information(msg, this_module, this_procedure)
        endif
    end subroutine kexpm_vec_kernel_csp

    !----------------------------------------------------------------
    !-----     Public driver : block version, fixed horizon     -----
    !----------------------------------------------------------------

    subroutine kexpm_mat_csp(C, A, B, tau, tol, info, trans, kdim, X)
        !! Best Krylov approximation of \( \mathbf{C} = \exp(\tau \mathbf{A}) \mathbf{B} \) in the computed Krylov subspace.
        implicit none(type, external)
        class(abstract_vector_csp), intent(out) :: C(:)
        !! Output vectors.
        class(abstract_linop_csp), intent(inout) :: A
        !! Linear operator to be exponentiated.
        class(abstract_vector_csp), intent(in) :: B(:)
        !! Input matrix on which to apply \( \exp(\tau \mathbf{A}) \).
        real(sp), intent(in) :: tau
        !! Time horizon for the exponentiation.
        real(sp), intent(in) :: tol
        !! Solution tolerance based on error estimates.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: trans
        !! Use transpose ?
        integer, optional, intent(in) :: kdim
        !! Maximum number of block-Arnoldi steps.
        class(abstract_vector_csp), optional, intent(inout) :: X(:)
        !! Optional workspace for the block Krylov basis (size multiple of size(B)).

        ! ----- Internal variables -----
        character(len=*), parameter :: this_procedure = 'kexpm_mat_csp'
        class(abstract_vector_csp), allocatable :: Xloc(:)
        integer :: p, nsteps, iostat
        character(len=256) :: msg

        p = size(B)

        ! Sanity checks (no work is done for invalid inputs).
        info = 0
        if (p < 1 .or. size(C) /= p) then
            info = -6
        else if (present(X)) then
            if (mod(size(X), p) /= 0 .or. size(X) < 2*p) info = -3
        endif
        if (info == 0 .and. present(kdim)) then
            if (kdim < 1) info = -4
        endif
        if (info == 0 .and. tol < 0.0_sp) info = -5
        if (info /= 0) then
            call log_kexpm_invalid_input(info, this_procedure)
            return
        endif

        if (present(X)) then
            nsteps = min(optval(kdim, size(X)/p - 1), size(X)/p - 1)
            call kexpm_mat_kernel_csp(C, A, B, tau, tol, info, optval(trans, .false.), X(:p*(nsteps+1)))
        else
            ! NOTE: kdim counts block steps. The basis holds p*(kdim+1) vectors
            !       (previously p*(p*kdim+1) because the block size was counted twice).
            nsteps = optval(kdim, kmax)
            allocate(Xloc(p*(nsteps+1)), mold=B(1), stat=iostat, errmsg=msg)
            call check_allocation(iostat, msg, this_module, this_procedure)
            call kexpm_mat_kernel_csp(C, A, B, tau, tol, info, optval(trans, .false.), Xloc)
        endif
    end subroutine kexpm_mat_csp

    subroutine kexpm_mat_kernel_csp(C, A, B, tau, tol, info, transpose, X)
        !! Worker for `kexpm_mat`.
        implicit none(type, external)
        class(abstract_vector_csp), intent(out) :: C(:)
        class(abstract_linop_csp), intent(inout) :: A
        class(abstract_vector_csp), intent(in) :: B(:)
        real(sp), intent(in) :: tau, tol
        integer, intent(out) :: info
        logical, intent(in) :: transpose
        class(abstract_vector_csp), intent(inout) :: X(:)

        ! ----- Internal variables -----
        character(len=*), parameter :: this_procedure = 'kexpm_mat_csp'
        complex(sp), allocatable :: H(:, :), E(:, :), R(:, :), coef(:, :)
        integer, allocatable :: perm(:)
        real(sp) :: err_est
        integer :: p, nk, k, kpm, kp, kpp, iostat
        logical :: breakdown, partial_breakdown
        character(len=256) :: msg

        p  = size(B)
        nk = size(X)/p - 1

        allocate(R(p, p), H(p*(nk+1), p*(nk+1)), E(p*(nk+1), p*(nk+1)), &
                 source=zero_csp, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)
        allocate(perm(p), source=0, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)

        ! Orthonormalize the input block directly inside the Krylov basis: B = X(:p) @ R.
        call zero_basis(X)
        call copy(X(:p), B)
        call qr(X(:p), R, perm, info)
        call check_info(info, 'qr_pivot', this_module, this_procedure)
        call permcols(R, invperm(perm))

        if (mnorm(R, "fro") == 0.0_sp) then
            ! Input matrix is zero.
            call zero_basis(C)
            info = p
            call log_information('Zero input block, zero output.', this_module, this_procedure)
            return
        endif

        err_est = huge(1.0_sp) ; kp = p ; kpp = p
        breakdown = .false. ; partial_breakdown = .false.

        expm_arnoldi: do k = 1, nk
            ! Set counters.
            kpm = (k-1)*p ; kp = kpm + p ; kpp = kp + p

            ! Compute the k-th step of the block Arnoldi factorization.
            call arnoldi(A, X, H, info, kstart=k, kend=k, transpose=transpose, blksize=p)
            call check_info(info, 'arnoldi', this_module, this_procedure)

            if (info == kp) then
                ! Block Arnoldi reports breakdown as soon as the new block is rank deficient.
                ! The subspace is invariant only if the whole coupling block vanishes.
                breakdown = mnorm(H(kp+1:kpp, kpm+1:kp), "fro") < atol_sp
                ! Otherwise, the QR completed the block with random vectors that are not
                ! orthogonal to the previous basis: we cannot safely continue (no deflation).
                partial_breakdown = .not. breakdown
            endif

            if (breakdown) then
                ! Exact invariant subspace: do not consider the extended matrix.
                kpp = kp
                call small_expm_csp(E, H, kpp, kpp, tau)
                err_est = 0.0_sp
                exit expm_arnoldi
            endif

            ! Dense exponential of the extended block-Hessenberg matrix.
            call small_expm_csp(E, H, kpp, kp, tau)

            ! Cheap (conservative) error estimate.
            err_est = mnorm(matmul(E(kp+1:kpp, :p), R), "fro")

            if (err_est <= tol .or. partial_breakdown) exit expm_arnoldi
        enddo expm_arnoldi

        ! Single projection back to the full space: C = X(:kpp) @ (E(:kpp, :p) @ R).
        coef = matmul(E(:kpp, :p), R)
        call project_block_csp(C, X(:kpp), coef)

        if (err_est <= tol) then
            info = kpp
            write(msg,'(A,I0,2(A,E9.2),A,L1)') 'Converged. kp= ', kpp, ', err_est= ', err_est, &
                                             & ', tol= ', tol, ', breakdown= ', breakdown
            call log_information(msg, this_module, this_procedure)
        else
            info = -1
            if (partial_breakdown) then
                write(msg,'(A,I0,2(A,E9.2))') 'Not converged (partial block breakdown, deflation not supported). kp= ', &
                                            & kpp, ', err_est= ', err_est, ', tol= ', tol
            else
                write(msg,'(A,I0,2(A,E9.2))') 'Not converged. kp= ', kpp, ', err_est= ', err_est, ', tol= ', tol
            endif
            call log_information(msg, this_module, this_procedure)
        endif
    end subroutine kexpm_mat_kernel_csp

    !--------------------------------------------------------------------
    !-----     Public driver : adaptive sub-stepping (Expokit-like)     -----
    !--------------------------------------------------------------------

    subroutine krylov_exptA_csp(vec_out, A, vec_in, tau, info, trans)
        !! Wrapper for the Krylov-based evaluation of the action of the matrix exponential operator
        !! on a vector that conforms to the `abstract_exptA_csp` interface.
        !! It uses `kexpm_substep` with default options and no warm start; call `kexpm_substep`
        !! directly with a persistent metadata object to re-use learned settings.
        implicit none(type, external)
        class(abstract_vector_csp), intent(out) :: vec_out
        !! Solution vector.
        class(abstract_linop_csp), intent(inout) :: A
        !! Linear operator to be exponentiated.
        class(abstract_vector_csp), intent(in) :: vec_in
        !! Input vector to be multiplied by \( \exp(\tau \mathbf{A}) \).
        real(sp), intent(in) :: tau
        !! Time horizon for the exponentiation.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: trans
        !! Use adjoint ?

        ! ----- Internal variables -----
        character(len=*), parameter :: this_procedure = 'krylov_exptA_csp'
        real(sp) :: tol

        ! Tolerance relative to the norm of the input (identical to the former atol for unit vectors).
        tol = atol_sp * vec_in%norm()

        call kexpm_substep_csp(vec_out, A, vec_in, tau, tol, info, trans=trans)
        call check_info(info, 'kexpm', this_module, this_procedure)
    end subroutine krylov_exptA_csp

    subroutine kexpm_substep_csp(c, A, b, tau, tol, info, trans, options, meta)
        !! Sub-stepped Krylov approximation of \( \exp(\tau \mathbf{A}) \mathbf{b} \) with optional
        !! warm start (see the `kexpm_substep` interface).
        implicit none(type, external)
        class(abstract_vector_csp), intent(out) :: c
        !! Output vector.
        class(abstract_linop_csp), intent(inout) :: A
        !! Linear operator to be exponentiated.
        class(abstract_vector_csp), intent(in) :: b
        !! Input vector.
        real(sp), intent(in) :: tau
        !! Time horizon (may be negative).
        real(sp), intent(in) :: tol
        !! Absolute tolerance on the accumulated error estimate.
        integer, intent(out) :: info
        !! Total number of matvecs on success, -1 otherwise.
        logical, optional, intent(in) :: trans
        !! Use adjoint?
        class(abstract_opts), optional, intent(in) :: options
        !! `kexpm_sp_opts`.
        class(abstract_metadata), optional, intent(inout) :: meta
        !! `kexpm_sp_metadata` (counters and learned settings, kept between calls).

        ! ----- Internal variables -----
        character(len=*), parameter :: this_procedure = 'kexpm_substep_csp'
        real(sp), parameter :: safety = 0.9_sp, fac_min = 0.2_sp, fac_max = 5.0_sp
        type(kexpm_sp_opts) :: opts
        type(kexpm_sp_metadata) :: lmeta
        class(abstract_vector_csp), allocatable :: X(:)
        complex(sp), allocatable :: H(:, :), E(:, :)
        real(sp) :: t_end, sgn, t, dt, beta, err, tol_loc, fac, floor_tol, tol_rel
        integer :: m, ksteps, kp, nsub, nmatvec, nrej, nrej_tot, n_uniform, iostat
        logical :: transpose, breakdown, last, uniform, used_uniform, fallback
        character(len=256) :: msg

        ! ----- Options and metadata -----
        opts = kexpm_sp_opts()
        if (present(options)) then
            select type (options)
            type is (kexpm_sp_opts)
                opts = options
            class default
                call type_error('options', 'kexpm_sp_opts', 'IN', this_module, this_procedure)
            end select
        endif
        if (present(meta)) then
            select type (meta)
            type is (kexpm_sp_metadata)
                lmeta = meta
            class default
                call type_error('meta', 'kexpm_sp_metadata', 'INOUT', this_module, this_procedure)
            end select
        endif

        ! Sanity checks (no work is done and the metadata is left untouched for invalid inputs).
        info = 0
        if (opts%kdim < 1) then
            info = -4
        else if (tol < 0.0_sp) then
            info = -5
        else if (opts%max_rejections < 1 .or. opts%max_substeps < 1) then
            info = -7
        endif
        if (info /= 0) then
            call log_kexpm_invalid_input(info, this_procedure)
            return
        endif

        transpose = optval(trans, .false.)
        m = opts%kdim
        t_end = abs(tau) ; sgn = sign(1.0_sp, tau)
        lmeta%n_calls = lmeta%n_calls + 1

        ! Learned settings are only valid for the same horizon, tolerance and Krylov dimension.
        ! The tolerance is compared relative to ||b|| (within a factor 2), so that calls with
        ! tol = rtol*||b|| and different right-hand sides share the learned settings.
        beta = b%norm()
        tol_rel = merge(tol/beta, 0.0_sp, beta > 0.0_sp)
        if (lmeta%learned) then
            if (abs(t_end - lmeta%tau_learned) > 10*epsilon(1.0_sp)*t_end &
                & .or. tol_rel > 2*lmeta%tol_learned .or. 2*tol_rel < lmeta%tol_learned &
                & .or. m /= lmeta%kdim_learned) call lmeta%forget()
        endif
        used_uniform = present(meta) .and. opts%warm_start .and. lmeta%learned .and. lmeta%nsub > 0
        n_uniform = merge(lmeta%nsub, 0, used_uniform)
        uniform = used_uniform ; fallback = .false.

        ! c holds the running solution w(t).
        call c%zero() ; call c%add(b)
        beta = c%norm()
        t = 0.0_sp ; nsub = 0 ; nmatvec = 0 ; nrej_tot = 0 ; info = 0

        if (beta /= 0.0_sp .and. t_end /= 0.0_sp) then
            ! Workspace (allocated once, reused by every sub-step).
            allocate(X(m+1), mold=b, stat=iostat, errmsg=msg)
            call check_allocation(iostat, msg, this_module, this_procedure)
            allocate(H(m+1, m+1), E(m+1, m+1), source=zero_csp, stat=iostat, errmsg=msg)
            call check_allocation(iostat, msg, this_module, this_procedure)

            dt = merge(t_end/max(n_uniform, 1), t_end, uniform)

            substeps: do while (t < t_end)
                nsub = nsub + 1
                if (nsub > opts%max_substeps) then
                    info = -1 ; exit substeps
                endif

                ! Last step (with a small slack so that round-off does not create a tiny extra step).
                last = (dt >= (t_end - t) - 10*epsilon(1.0_sp)*t_end) ; if (last) dt = t_end - t

                beta = c%norm()
                if (beta == 0.0_sp) exit substeps   ! Solution is exactly zero from now on.

                ! Roundoff floor on the local tolerance to avoid endless rejections.
                floor_tol = 2*epsilon(1.0_sp)*beta
                tol_loc = max(tol*dt/t_end, floor_tol)

                ! Fresh Krylov basis started from the current solution.
                call zero_basis(X) ; H = zero_csp
                call X(1)%add(c) ; call X(1)%scal(one_csp / beta)

                call arnoldi_expm_csp(A, X, H, E, beta, sgn*dt, tol_loc, transpose, &
                                                      & ksteps, kp, err, breakdown)
                nmatvec = nmatvec + ksteps

                if (breakdown) then
                    ! Invariant subspace: exact for any dt, finish the horizon in one go.
                    dt = t_end - t ; last = .true.
                    call small_expm_csp(E, H, kp, kp, sgn*dt)
                    err = 0.0_sp
                else if (err > tol_loc) then
                    ! Krylov basis exhausted. Safety net: a learned uniform step that is too large
                    ! switches the rest of this call to the adaptive controller.
                    if (uniform) then
                        uniform = .false. ; fallback = .true.
                    endif
                    ! Shrink dt reusing H (only small dense exponentials).
                    last = .false.
                    do nrej = 1, opts%max_rejections
                        fac = min(safety, max(0.1_sp, safety*(tol_loc/err)**(1.0_sp/m)))
                        dt  = fac*dt
                        tol_loc = max(tol*dt/t_end, floor_tol)
                        call small_expm_csp(E, H, kp, ksteps, sgn*dt)
                        err = beta*abs(E(kp, 1))
                        nrej_tot = nrej_tot + 1
                        if (err <= tol_loc) exit
                    enddo
                    if (err > tol_loc) then
                        info = -1 ; exit substeps
                    endif
                endif

                ! Advance: w(t+dt) = beta * X(:kp) @ E(:kp, 1).
                call project_csp(c, X(:kp), beta*E(:kp, 1))
                if (last) then
                    t = t_end
                else
                    t = t + dt
                endif

                ! Next step size.
                if (uniform) then
                    dt = t_end/n_uniform
                else
                    ! Error controller (local error ~ dt^ksteps).
                    if (err > 0.0_sp) then
                        fac = safety*(tol_loc/err)**(1.0_sp/max(ksteps, 1))
                    else
                        fac = fac_max
                    endif
                    dt = dt*min(fac_max, max(fac_min, fac))
                endif
            enddo substeps
        endif

        if (info == 0) info = max(nmatvec, 1)

        ! ----- Update metadata (counters, learned settings) -----
        lmeta%n_matvec     = lmeta%n_matvec + nmatvec
        lmeta%n_substeps   = lmeta%n_substeps + nsub
        lmeta%n_rejections = lmeta%n_rejections + nrej_tot
        if (fallback) lmeta%n_fallbacks = lmeta%n_fallbacks + 1
        lmeta%info = info ; lmeta%converged = (info > 0)
        if (present(meta) .and. opts%warm_start .and. info > 0 .and. nsub > 0) then
            call learn_kexpm_sp(lmeta, used_uniform, n_uniform, nsub, nmatvec, fallback, &
                                      & opts%tune, t_end, tol_rel, m)
        endif

        if (info > 0) then
            write(msg,'(A,I0,A,I0,A,I0,A,L1,A,L1)') 'Converged. substeps= ', nsub, ', matvecs= ', nmatvec, &
                & ', rejections= ', nrej_tot, ', warm start= ', used_uniform, ', fallback= ', fallback
        else
            write(msg,'(A,I0,A,E9.2,A,E9.2)') 'Not converged. substeps= ', nsub, ', t= ', sgn*t, ', tau= ', tau
        endif
        call log_information(msg, this_module, this_procedure)

        if (opts%if_print_metadata) call lmeta%print()
        if (present(meta)) then
            select type (meta)
            type is (kexpm_sp_metadata)
                meta = lmeta
            end select
        endif
    end subroutine kexpm_substep_csp

    !----------------------------------------
    !-----     Private helper kernels     -----
    !----------------------------------------

    subroutine arnoldi_expm_csp(A, X, H, E, beta, tau, tol, transpose, ksteps, kp, err_est, breakdown)
        !! Grows the Arnoldi basis `X` (with `X(1)` already normalized) one vector at a time and
        !! evaluates the exponential of the extended Hessenberg matrix after each step, until
        !! `beta*|E(kp,1)| <= tol`, breakdown, or `size(X)-1` steps.
        !! On exit `E(:kp, :kp)` holds the dense exponential used for the projection.
        implicit none(type, external)
        class(abstract_linop_csp), intent(inout) :: A
        class(abstract_vector_csp), intent(inout) :: X(:)
        complex(sp), intent(inout) :: H(:, :), E(:, :)
        real(sp), intent(in) :: beta, tau, tol
        logical, intent(in) :: transpose
        integer, intent(out) :: ksteps, kp
        real(sp), intent(out) :: err_est
        logical, intent(out) :: breakdown

        character(len=*), parameter :: this_procedure = 'arnoldi_expm_csp'
        integer :: k, info

        err_est = huge(1.0_sp) ; breakdown = .false. ; ksteps = 0 ; kp = 1

        do k = 1, size(X) - 1
            call arnoldi(A, X, H, info, kstart=k, kend=k, transpose=transpose)
            call check_info(info, 'arnoldi', this_module, this_procedure)
            ksteps = k
            breakdown = (info == k)

            if (breakdown) then
                ! Invariant subspace: the approximation is exact, drop the extended row.
                kp = k
                call small_expm_csp(E, H, kp, kp, tau)
                err_est = 0.0_sp
                return
            endif

            ! Exponential of the (k+1)x(k+1) extended Hessenberg matrix [H_k 0; h e_k^T 0].
            kp = k + 1
            call small_expm_csp(E, H, kp, k, tau)

            ! Cheap error estimate (magnitude of the last included correction, conservative).
            err_est = beta*abs(E(kp, 1))
            if (err_est <= tol) return
        enddo
    end subroutine arnoldi_expm_csp

    subroutine small_expm_csp(E, H, n, ncols, tau)
        !! E(:n, :n) = exp(tau * Hbar) where Hbar holds the first `ncols` columns of H(:n, :)
        !! and zeros elsewhere. Only the leading n x n block of E is touched (no full reset),
        !! and stale columns of H (e.g. from a previous sub-step) are never used.
        implicit none(type, external)
        complex(sp), intent(inout) :: E(:, :)
        complex(sp), intent(in) :: H(:, :)
        integer, intent(in) :: n, ncols
        real(sp), intent(in) :: tau

        E(:n, :ncols) = tau*H(:n, :ncols)
        if (ncols < n) E(:n, ncols+1:n) = zero_csp
        call matrix_exp(E(:n, :n))
    end subroutine small_expm_csp

    subroutine project_csp(c, X, v)
        !! c = X @ v, written directly into the (non-allocatable) output vector.
        implicit none(type, external)
        class(abstract_vector_csp), intent(inout) :: c
        class(abstract_vector_csp), intent(in) :: X(:)
        complex(sp), intent(in) :: v(:)
        integer :: i
        call c%zero()
        do i = 1, size(X)
            call c%axpby(v(i), X(i), one_csp)
        enddo
    end subroutine project_csp

    subroutine project_block_csp(C, X, V)
        !! C = X @ V for a block of output vectors.
        implicit none(type, external)
        class(abstract_vector_csp), intent(inout) :: C(:)
        class(abstract_vector_csp), intent(in) :: X(:)
        complex(sp), intent(in) :: V(:, :)
        integer :: j
        do j = 1, size(C)
            call project_csp(C(j), X, V(:, j))
        enddo
    end subroutine project_block_csp

    !-------------------------------------------------------------------
    !-----     Public driver : single vector, fixed time horizon     -----
    !-------------------------------------------------------------------

    subroutine kexpm_vec_cdp(c, A, b, tau, tol, info, trans, kdim, X)
        !! Best approximation of \( \exp(\tau \mathbf{A}) \mathbf{b} \) in the computed Krylov subspace
        implicit none(type, external)
        class(abstract_vector_cdp), intent(out) :: c
        !! Output vector.
        class(abstract_linop_cdp), intent(inout) :: A
        !! Linear operator to be exponentiated.
        class(abstract_vector_cdp), intent(in) :: b
        !! Input vector on which to apply \( \exp(\tau \mathbf{A}) \).
        real(dp), intent(in) :: tau
        !! Time horizon for the exponentiation.
        real(dp), intent(in) :: tol
        !! Solution tolerance based on error estimates.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: trans
        !! Use transpose?
        integer, optional, intent(in) :: kdim
        !! Maximum number of Arnoldi steps.
        class(abstract_vector_cdp), optional, intent(inout) :: X(:)
        !! Optional workspace for the Krylov basis (avoids allocating it at every call).

        ! ----- Internal variables -----
        character(len=*), parameter :: this_procedure = 'kexpm_vec_cdp'
        class(abstract_vector_cdp), allocatable :: Xloc(:)
        integer :: nsteps, iostat
        character(len=256) :: msg

        ! Sanity checks (no work is done for invalid inputs).
        info = 0
        if (present(X)) then
            if (size(X) < 2) info = -3
        endif
        if (info == 0 .and. present(kdim)) then
            if (kdim < 1) info = -4
        endif
        if (info == 0 .and. tol < 0.0_dp) info = -5
        if (info /= 0) then
            call log_kexpm_invalid_input(info, this_procedure)
            return
        endif

        if (present(X)) then
            nsteps = min(optval(kdim, size(X)-1), size(X)-1)
            call kexpm_vec_kernel_cdp(c, A, b, tau, tol, info, optval(trans, .false.), X(:nsteps+1))
        else
            nsteps = optval(kdim, kmax)
            ! `mold` : no copy of b into every basis vector (the kernel zeroes the basis anyway).
            allocate(Xloc(nsteps+1), mold=b, stat=iostat, errmsg=msg)
            call check_allocation(iostat, msg, this_module, this_procedure)
            call kexpm_vec_kernel_cdp(c, A, b, tau, tol, info, optval(trans, .false.), Xloc)
            ! Xloc is allocatable: automatically deallocated on return (no leak).
        endif
    end subroutine kexpm_vec_cdp

    subroutine kexpm_vec_kernel_cdp(c, A, b, tau, tol, info, transpose, X)
        !! Worker for `kexpm_vec`: the Krylov basis `X` is a plain (non-optional) dummy argument.
        implicit none(type, external)
        class(abstract_vector_cdp), intent(out) :: c
        class(abstract_linop_cdp), intent(inout) :: A
        class(abstract_vector_cdp), intent(in) :: b
        real(dp), intent(in) :: tau, tol
        integer, intent(out) :: info
        logical, intent(in) :: transpose
        class(abstract_vector_cdp), intent(inout) :: X(:)

        ! ----- Internal variables -----
        character(len=*), parameter :: this_procedure = 'kexpm_vec_cdp'
        complex(dp), allocatable :: H(:, :), E(:, :)
        real(dp) :: beta, err_est
        integer :: nk, kp, ksteps, iostat
        logical :: breakdown
        character(len=256) :: msg

        nk = size(X) - 1

        beta = b%norm()
        if (beta == 0.0_dp) then
            ! Input is zero => Output is zero.
            call c%zero()
            info = 1
            call log_information('Zero input vector, zero output.', this_module, this_procedure)
            return
        endif

        allocate(H(nk+1, nk+1), E(nk+1, nk+1), source=zero_cdp, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)

        ! Initialize Krylov subspace with the normalized input vector.
        call zero_basis(X)
        call X(1)%add(b) ; call X(1)%scal(one_cdp / beta)

        ! Grow the Krylov subspace until the error estimate is below tol.
        call arnoldi_expm_cdp(A, X, H, E, beta, tau, tol, transpose, &
                                              & ksteps, kp, err_est, breakdown)

        ! Single projection back onto the full space: c = beta * X(:kp) @ E(:kp, 1).
        call project_cdp(c, X(:kp), beta*E(:kp, 1))

        if (err_est <= tol) then
            info = kp
            write(msg,'(A,I0,2(A,E9.2),A,L1)') 'Converged. kp= ', kp, ', err_est= ', err_est, &
                                             & ', tol= ', tol, ', breakdown= ', breakdown
            call log_information(msg, this_module, this_procedure)
        else
            info = -1
            write(msg,'(A,I0,2(A,E9.2))') 'Not converged. kp= ', kp, ', err_est= ', err_est, ', tol= ', tol
            call log_information(msg, this_module, this_procedure)
        endif
    end subroutine kexpm_vec_kernel_cdp

    !----------------------------------------------------------------
    !-----     Public driver : block version, fixed horizon     -----
    !----------------------------------------------------------------

    subroutine kexpm_mat_cdp(C, A, B, tau, tol, info, trans, kdim, X)
        !! Best Krylov approximation of \( \mathbf{C} = \exp(\tau \mathbf{A}) \mathbf{B} \) in the computed Krylov subspace.
        implicit none(type, external)
        class(abstract_vector_cdp), intent(out) :: C(:)
        !! Output vectors.
        class(abstract_linop_cdp), intent(inout) :: A
        !! Linear operator to be exponentiated.
        class(abstract_vector_cdp), intent(in) :: B(:)
        !! Input matrix on which to apply \( \exp(\tau \mathbf{A}) \).
        real(dp), intent(in) :: tau
        !! Time horizon for the exponentiation.
        real(dp), intent(in) :: tol
        !! Solution tolerance based on error estimates.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: trans
        !! Use transpose ?
        integer, optional, intent(in) :: kdim
        !! Maximum number of block-Arnoldi steps.
        class(abstract_vector_cdp), optional, intent(inout) :: X(:)
        !! Optional workspace for the block Krylov basis (size multiple of size(B)).

        ! ----- Internal variables -----
        character(len=*), parameter :: this_procedure = 'kexpm_mat_cdp'
        class(abstract_vector_cdp), allocatable :: Xloc(:)
        integer :: p, nsteps, iostat
        character(len=256) :: msg

        p = size(B)

        ! Sanity checks (no work is done for invalid inputs).
        info = 0
        if (p < 1 .or. size(C) /= p) then
            info = -6
        else if (present(X)) then
            if (mod(size(X), p) /= 0 .or. size(X) < 2*p) info = -3
        endif
        if (info == 0 .and. present(kdim)) then
            if (kdim < 1) info = -4
        endif
        if (info == 0 .and. tol < 0.0_dp) info = -5
        if (info /= 0) then
            call log_kexpm_invalid_input(info, this_procedure)
            return
        endif

        if (present(X)) then
            nsteps = min(optval(kdim, size(X)/p - 1), size(X)/p - 1)
            call kexpm_mat_kernel_cdp(C, A, B, tau, tol, info, optval(trans, .false.), X(:p*(nsteps+1)))
        else
            ! NOTE: kdim counts block steps. The basis holds p*(kdim+1) vectors
            !       (previously p*(p*kdim+1) because the block size was counted twice).
            nsteps = optval(kdim, kmax)
            allocate(Xloc(p*(nsteps+1)), mold=B(1), stat=iostat, errmsg=msg)
            call check_allocation(iostat, msg, this_module, this_procedure)
            call kexpm_mat_kernel_cdp(C, A, B, tau, tol, info, optval(trans, .false.), Xloc)
        endif
    end subroutine kexpm_mat_cdp

    subroutine kexpm_mat_kernel_cdp(C, A, B, tau, tol, info, transpose, X)
        !! Worker for `kexpm_mat`.
        implicit none(type, external)
        class(abstract_vector_cdp), intent(out) :: C(:)
        class(abstract_linop_cdp), intent(inout) :: A
        class(abstract_vector_cdp), intent(in) :: B(:)
        real(dp), intent(in) :: tau, tol
        integer, intent(out) :: info
        logical, intent(in) :: transpose
        class(abstract_vector_cdp), intent(inout) :: X(:)

        ! ----- Internal variables -----
        character(len=*), parameter :: this_procedure = 'kexpm_mat_cdp'
        complex(dp), allocatable :: H(:, :), E(:, :), R(:, :), coef(:, :)
        integer, allocatable :: perm(:)
        real(dp) :: err_est
        integer :: p, nk, k, kpm, kp, kpp, iostat
        logical :: breakdown, partial_breakdown
        character(len=256) :: msg

        p  = size(B)
        nk = size(X)/p - 1

        allocate(R(p, p), H(p*(nk+1), p*(nk+1)), E(p*(nk+1), p*(nk+1)), &
                 source=zero_cdp, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)
        allocate(perm(p), source=0, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)

        ! Orthonormalize the input block directly inside the Krylov basis: B = X(:p) @ R.
        call zero_basis(X)
        call copy(X(:p), B)
        call qr(X(:p), R, perm, info)
        call check_info(info, 'qr_pivot', this_module, this_procedure)
        call permcols(R, invperm(perm))

        if (mnorm(R, "fro") == 0.0_dp) then
            ! Input matrix is zero.
            call zero_basis(C)
            info = p
            call log_information('Zero input block, zero output.', this_module, this_procedure)
            return
        endif

        err_est = huge(1.0_dp) ; kp = p ; kpp = p
        breakdown = .false. ; partial_breakdown = .false.

        expm_arnoldi: do k = 1, nk
            ! Set counters.
            kpm = (k-1)*p ; kp = kpm + p ; kpp = kp + p

            ! Compute the k-th step of the block Arnoldi factorization.
            call arnoldi(A, X, H, info, kstart=k, kend=k, transpose=transpose, blksize=p)
            call check_info(info, 'arnoldi', this_module, this_procedure)

            if (info == kp) then
                ! Block Arnoldi reports breakdown as soon as the new block is rank deficient.
                ! The subspace is invariant only if the whole coupling block vanishes.
                breakdown = mnorm(H(kp+1:kpp, kpm+1:kp), "fro") < atol_dp
                ! Otherwise, the QR completed the block with random vectors that are not
                ! orthogonal to the previous basis: we cannot safely continue (no deflation).
                partial_breakdown = .not. breakdown
            endif

            if (breakdown) then
                ! Exact invariant subspace: do not consider the extended matrix.
                kpp = kp
                call small_expm_cdp(E, H, kpp, kpp, tau)
                err_est = 0.0_dp
                exit expm_arnoldi
            endif

            ! Dense exponential of the extended block-Hessenberg matrix.
            call small_expm_cdp(E, H, kpp, kp, tau)

            ! Cheap (conservative) error estimate.
            err_est = mnorm(matmul(E(kp+1:kpp, :p), R), "fro")

            if (err_est <= tol .or. partial_breakdown) exit expm_arnoldi
        enddo expm_arnoldi

        ! Single projection back to the full space: C = X(:kpp) @ (E(:kpp, :p) @ R).
        coef = matmul(E(:kpp, :p), R)
        call project_block_cdp(C, X(:kpp), coef)

        if (err_est <= tol) then
            info = kpp
            write(msg,'(A,I0,2(A,E9.2),A,L1)') 'Converged. kp= ', kpp, ', err_est= ', err_est, &
                                             & ', tol= ', tol, ', breakdown= ', breakdown
            call log_information(msg, this_module, this_procedure)
        else
            info = -1
            if (partial_breakdown) then
                write(msg,'(A,I0,2(A,E9.2))') 'Not converged (partial block breakdown, deflation not supported). kp= ', &
                                            & kpp, ', err_est= ', err_est, ', tol= ', tol
            else
                write(msg,'(A,I0,2(A,E9.2))') 'Not converged. kp= ', kpp, ', err_est= ', err_est, ', tol= ', tol
            endif
            call log_information(msg, this_module, this_procedure)
        endif
    end subroutine kexpm_mat_kernel_cdp

    !--------------------------------------------------------------------
    !-----     Public driver : adaptive sub-stepping (Expokit-like)     -----
    !--------------------------------------------------------------------

    subroutine krylov_exptA_cdp(vec_out, A, vec_in, tau, info, trans)
        !! Wrapper for the Krylov-based evaluation of the action of the matrix exponential operator
        !! on a vector that conforms to the `abstract_exptA_cdp` interface.
        !! It uses `kexpm_substep` with default options and no warm start; call `kexpm_substep`
        !! directly with a persistent metadata object to re-use learned settings.
        implicit none(type, external)
        class(abstract_vector_cdp), intent(out) :: vec_out
        !! Solution vector.
        class(abstract_linop_cdp), intent(inout) :: A
        !! Linear operator to be exponentiated.
        class(abstract_vector_cdp), intent(in) :: vec_in
        !! Input vector to be multiplied by \( \exp(\tau \mathbf{A}) \).
        real(dp), intent(in) :: tau
        !! Time horizon for the exponentiation.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: trans
        !! Use adjoint ?

        ! ----- Internal variables -----
        character(len=*), parameter :: this_procedure = 'krylov_exptA_cdp'
        real(dp) :: tol

        ! Tolerance relative to the norm of the input (identical to the former atol for unit vectors).
        tol = atol_dp * vec_in%norm()

        call kexpm_substep_cdp(vec_out, A, vec_in, tau, tol, info, trans=trans)
        call check_info(info, 'kexpm', this_module, this_procedure)
    end subroutine krylov_exptA_cdp

    subroutine kexpm_substep_cdp(c, A, b, tau, tol, info, trans, options, meta)
        !! Sub-stepped Krylov approximation of \( \exp(\tau \mathbf{A}) \mathbf{b} \) with optional
        !! warm start (see the `kexpm_substep` interface).
        implicit none(type, external)
        class(abstract_vector_cdp), intent(out) :: c
        !! Output vector.
        class(abstract_linop_cdp), intent(inout) :: A
        !! Linear operator to be exponentiated.
        class(abstract_vector_cdp), intent(in) :: b
        !! Input vector.
        real(dp), intent(in) :: tau
        !! Time horizon (may be negative).
        real(dp), intent(in) :: tol
        !! Absolute tolerance on the accumulated error estimate.
        integer, intent(out) :: info
        !! Total number of matvecs on success, -1 otherwise.
        logical, optional, intent(in) :: trans
        !! Use adjoint?
        class(abstract_opts), optional, intent(in) :: options
        !! `kexpm_dp_opts`.
        class(abstract_metadata), optional, intent(inout) :: meta
        !! `kexpm_dp_metadata` (counters and learned settings, kept between calls).

        ! ----- Internal variables -----
        character(len=*), parameter :: this_procedure = 'kexpm_substep_cdp'
        real(dp), parameter :: safety = 0.9_dp, fac_min = 0.2_dp, fac_max = 5.0_dp
        type(kexpm_dp_opts) :: opts
        type(kexpm_dp_metadata) :: lmeta
        class(abstract_vector_cdp), allocatable :: X(:)
        complex(dp), allocatable :: H(:, :), E(:, :)
        real(dp) :: t_end, sgn, t, dt, beta, err, tol_loc, fac, floor_tol, tol_rel
        integer :: m, ksteps, kp, nsub, nmatvec, nrej, nrej_tot, n_uniform, iostat
        logical :: transpose, breakdown, last, uniform, used_uniform, fallback
        character(len=256) :: msg

        ! ----- Options and metadata -----
        opts = kexpm_dp_opts()
        if (present(options)) then
            select type (options)
            type is (kexpm_dp_opts)
                opts = options
            class default
                call type_error('options', 'kexpm_dp_opts', 'IN', this_module, this_procedure)
            end select
        endif
        if (present(meta)) then
            select type (meta)
            type is (kexpm_dp_metadata)
                lmeta = meta
            class default
                call type_error('meta', 'kexpm_dp_metadata', 'INOUT', this_module, this_procedure)
            end select
        endif

        ! Sanity checks (no work is done and the metadata is left untouched for invalid inputs).
        info = 0
        if (opts%kdim < 1) then
            info = -4
        else if (tol < 0.0_dp) then
            info = -5
        else if (opts%max_rejections < 1 .or. opts%max_substeps < 1) then
            info = -7
        endif
        if (info /= 0) then
            call log_kexpm_invalid_input(info, this_procedure)
            return
        endif

        transpose = optval(trans, .false.)
        m = opts%kdim
        t_end = abs(tau) ; sgn = sign(1.0_dp, tau)
        lmeta%n_calls = lmeta%n_calls + 1

        ! Learned settings are only valid for the same horizon, tolerance and Krylov dimension.
        ! The tolerance is compared relative to ||b|| (within a factor 2), so that calls with
        ! tol = rtol*||b|| and different right-hand sides share the learned settings.
        beta = b%norm()
        tol_rel = merge(tol/beta, 0.0_dp, beta > 0.0_dp)
        if (lmeta%learned) then
            if (abs(t_end - lmeta%tau_learned) > 10*epsilon(1.0_dp)*t_end &
                & .or. tol_rel > 2*lmeta%tol_learned .or. 2*tol_rel < lmeta%tol_learned &
                & .or. m /= lmeta%kdim_learned) call lmeta%forget()
        endif
        used_uniform = present(meta) .and. opts%warm_start .and. lmeta%learned .and. lmeta%nsub > 0
        n_uniform = merge(lmeta%nsub, 0, used_uniform)
        uniform = used_uniform ; fallback = .false.

        ! c holds the running solution w(t).
        call c%zero() ; call c%add(b)
        beta = c%norm()
        t = 0.0_dp ; nsub = 0 ; nmatvec = 0 ; nrej_tot = 0 ; info = 0

        if (beta /= 0.0_dp .and. t_end /= 0.0_dp) then
            ! Workspace (allocated once, reused by every sub-step).
            allocate(X(m+1), mold=b, stat=iostat, errmsg=msg)
            call check_allocation(iostat, msg, this_module, this_procedure)
            allocate(H(m+1, m+1), E(m+1, m+1), source=zero_cdp, stat=iostat, errmsg=msg)
            call check_allocation(iostat, msg, this_module, this_procedure)

            dt = merge(t_end/max(n_uniform, 1), t_end, uniform)

            substeps: do while (t < t_end)
                nsub = nsub + 1
                if (nsub > opts%max_substeps) then
                    info = -1 ; exit substeps
                endif

                ! Last step (with a small slack so that round-off does not create a tiny extra step).
                last = (dt >= (t_end - t) - 10*epsilon(1.0_dp)*t_end) ; if (last) dt = t_end - t

                beta = c%norm()
                if (beta == 0.0_dp) exit substeps   ! Solution is exactly zero from now on.

                ! Roundoff floor on the local tolerance to avoid endless rejections.
                floor_tol = 2*epsilon(1.0_dp)*beta
                tol_loc = max(tol*dt/t_end, floor_tol)

                ! Fresh Krylov basis started from the current solution.
                call zero_basis(X) ; H = zero_cdp
                call X(1)%add(c) ; call X(1)%scal(one_cdp / beta)

                call arnoldi_expm_cdp(A, X, H, E, beta, sgn*dt, tol_loc, transpose, &
                                                      & ksteps, kp, err, breakdown)
                nmatvec = nmatvec + ksteps

                if (breakdown) then
                    ! Invariant subspace: exact for any dt, finish the horizon in one go.
                    dt = t_end - t ; last = .true.
                    call small_expm_cdp(E, H, kp, kp, sgn*dt)
                    err = 0.0_dp
                else if (err > tol_loc) then
                    ! Krylov basis exhausted. Safety net: a learned uniform step that is too large
                    ! switches the rest of this call to the adaptive controller.
                    if (uniform) then
                        uniform = .false. ; fallback = .true.
                    endif
                    ! Shrink dt reusing H (only small dense exponentials).
                    last = .false.
                    do nrej = 1, opts%max_rejections
                        fac = min(safety, max(0.1_dp, safety*(tol_loc/err)**(1.0_dp/m)))
                        dt  = fac*dt
                        tol_loc = max(tol*dt/t_end, floor_tol)
                        call small_expm_cdp(E, H, kp, ksteps, sgn*dt)
                        err = beta*abs(E(kp, 1))
                        nrej_tot = nrej_tot + 1
                        if (err <= tol_loc) exit
                    enddo
                    if (err > tol_loc) then
                        info = -1 ; exit substeps
                    endif
                endif

                ! Advance: w(t+dt) = beta * X(:kp) @ E(:kp, 1).
                call project_cdp(c, X(:kp), beta*E(:kp, 1))
                if (last) then
                    t = t_end
                else
                    t = t + dt
                endif

                ! Next step size.
                if (uniform) then
                    dt = t_end/n_uniform
                else
                    ! Error controller (local error ~ dt^ksteps).
                    if (err > 0.0_dp) then
                        fac = safety*(tol_loc/err)**(1.0_dp/max(ksteps, 1))
                    else
                        fac = fac_max
                    endif
                    dt = dt*min(fac_max, max(fac_min, fac))
                endif
            enddo substeps
        endif

        if (info == 0) info = max(nmatvec, 1)

        ! ----- Update metadata (counters, learned settings) -----
        lmeta%n_matvec     = lmeta%n_matvec + nmatvec
        lmeta%n_substeps   = lmeta%n_substeps + nsub
        lmeta%n_rejections = lmeta%n_rejections + nrej_tot
        if (fallback) lmeta%n_fallbacks = lmeta%n_fallbacks + 1
        lmeta%info = info ; lmeta%converged = (info > 0)
        if (present(meta) .and. opts%warm_start .and. info > 0 .and. nsub > 0) then
            call learn_kexpm_dp(lmeta, used_uniform, n_uniform, nsub, nmatvec, fallback, &
                                      & opts%tune, t_end, tol_rel, m)
        endif

        if (info > 0) then
            write(msg,'(A,I0,A,I0,A,I0,A,L1,A,L1)') 'Converged. substeps= ', nsub, ', matvecs= ', nmatvec, &
                & ', rejections= ', nrej_tot, ', warm start= ', used_uniform, ', fallback= ', fallback
        else
            write(msg,'(A,I0,A,E9.2,A,E9.2)') 'Not converged. substeps= ', nsub, ', t= ', sgn*t, ', tau= ', tau
        endif
        call log_information(msg, this_module, this_procedure)

        if (opts%if_print_metadata) call lmeta%print()
        if (present(meta)) then
            select type (meta)
            type is (kexpm_dp_metadata)
                meta = lmeta
            end select
        endif
    end subroutine kexpm_substep_cdp

    !----------------------------------------
    !-----     Private helper kernels     -----
    !----------------------------------------

    subroutine arnoldi_expm_cdp(A, X, H, E, beta, tau, tol, transpose, ksteps, kp, err_est, breakdown)
        !! Grows the Arnoldi basis `X` (with `X(1)` already normalized) one vector at a time and
        !! evaluates the exponential of the extended Hessenberg matrix after each step, until
        !! `beta*|E(kp,1)| <= tol`, breakdown, or `size(X)-1` steps.
        !! On exit `E(:kp, :kp)` holds the dense exponential used for the projection.
        implicit none(type, external)
        class(abstract_linop_cdp), intent(inout) :: A
        class(abstract_vector_cdp), intent(inout) :: X(:)
        complex(dp), intent(inout) :: H(:, :), E(:, :)
        real(dp), intent(in) :: beta, tau, tol
        logical, intent(in) :: transpose
        integer, intent(out) :: ksteps, kp
        real(dp), intent(out) :: err_est
        logical, intent(out) :: breakdown

        character(len=*), parameter :: this_procedure = 'arnoldi_expm_cdp'
        integer :: k, info

        err_est = huge(1.0_dp) ; breakdown = .false. ; ksteps = 0 ; kp = 1

        do k = 1, size(X) - 1
            call arnoldi(A, X, H, info, kstart=k, kend=k, transpose=transpose)
            call check_info(info, 'arnoldi', this_module, this_procedure)
            ksteps = k
            breakdown = (info == k)

            if (breakdown) then
                ! Invariant subspace: the approximation is exact, drop the extended row.
                kp = k
                call small_expm_cdp(E, H, kp, kp, tau)
                err_est = 0.0_dp
                return
            endif

            ! Exponential of the (k+1)x(k+1) extended Hessenberg matrix [H_k 0; h e_k^T 0].
            kp = k + 1
            call small_expm_cdp(E, H, kp, k, tau)

            ! Cheap error estimate (magnitude of the last included correction, conservative).
            err_est = beta*abs(E(kp, 1))
            if (err_est <= tol) return
        enddo
    end subroutine arnoldi_expm_cdp

    subroutine small_expm_cdp(E, H, n, ncols, tau)
        !! E(:n, :n) = exp(tau * Hbar) where Hbar holds the first `ncols` columns of H(:n, :)
        !! and zeros elsewhere. Only the leading n x n block of E is touched (no full reset),
        !! and stale columns of H (e.g. from a previous sub-step) are never used.
        implicit none(type, external)
        complex(dp), intent(inout) :: E(:, :)
        complex(dp), intent(in) :: H(:, :)
        integer, intent(in) :: n, ncols
        real(dp), intent(in) :: tau

        E(:n, :ncols) = tau*H(:n, :ncols)
        if (ncols < n) E(:n, ncols+1:n) = zero_cdp
        call matrix_exp(E(:n, :n))
    end subroutine small_expm_cdp

    subroutine project_cdp(c, X, v)
        !! c = X @ v, written directly into the (non-allocatable) output vector.
        implicit none(type, external)
        class(abstract_vector_cdp), intent(inout) :: c
        class(abstract_vector_cdp), intent(in) :: X(:)
        complex(dp), intent(in) :: v(:)
        integer :: i
        call c%zero()
        do i = 1, size(X)
            call c%axpby(v(i), X(i), one_cdp)
        enddo
    end subroutine project_cdp

    subroutine project_block_cdp(C, X, V)
        !! C = X @ V for a block of output vectors.
        implicit none(type, external)
        class(abstract_vector_cdp), intent(inout) :: C(:)
        class(abstract_vector_cdp), intent(in) :: X(:)
        complex(dp), intent(in) :: V(:, :)
        integer :: j
        do j = 1, size(C)
            call project_cdp(C(j), X, V(:, j))
        enddo
    end subroutine project_block_cdp

    !------------------------------------------------------------------
    !-----     Metadata of kexpm_substep : learning and printing     -----
    !------------------------------------------------------------------

    subroutine learn_kexpm_sp(self, used_uniform, n_uniform, nsub, mv, fallback, tune, t_end, tol_rel, m)
        !! Updates the learned number of uniform sub-steps after a successful call.
        !!
        !! The cost (matvecs) of N uniform sub-steps increases with N, while too small an N makes
        !! the error check fail (safety net). The cheapest setting is thus the smallest N that does
        !! not need the safety net. The search is:
        !!
        !! - Cold (adaptive) call: record its cost and number of sub-steps `nsub0`; next call
        !!   tries `nsub0 - 1` uniform sub-steps.
        !! - Uniform call, no fallback, cheaper than the best so far: keep it, try one sub-step less.
        !! - Uniform call with fallback, or not cheaper: stop and freeze on the best setting
        !!   (`nsub = 0`, i.e. fully adaptive, if no uniform setting was cheaper).
        !! - Fallback while frozen (e.g. the operator became stiffer): raise the count to
        !!   `max(nsub+1, sub-steps taken by this call)` and, with `tune`, restart the descent.
        implicit none(type, external)
        class(kexpm_sp_metadata), intent(inout) :: self
        logical, intent(in) :: used_uniform, fallback, tune
        integer, intent(in) :: n_uniform, nsub, mv, m
        real(sp), intent(in) :: t_end, tol_rel

        if (.not. self%learned) then
            ! Cold adaptive call.
            self%learned = .true.
            self%tau_learned = t_end ; self%tol_learned = tol_rel ; self%kdim_learned = m
            self%nsub0 = nsub ; self%nsub_best = 0 ; self%mv_best = mv
            if (tune .and. nsub > 1) then
                self%nsub = nsub - 1 ; self%frozen = .false.
            else
                self%nsub = nsub ; self%frozen = .true.
            endif
            return
        endif

        if (.not. used_uniform) return

        if (self%frozen) then
            if (fallback) then
                ! The frozen count became too small: restart the descent from a safe count.
                self%nsub = max(n_uniform + 1, nsub)
                self%nsub_best = 0 ; self%mv_best = huge(1)
                self%frozen = .not. tune
            endif
            return
        endif

        if ((.not. fallback) .and. mv < self%mv_best) then
            self%nsub_best = n_uniform ; self%mv_best = mv
            if (n_uniform > 1) then
                self%nsub = n_uniform - 1
                return
            endif
        endif

        ! Search over: freeze on the best setting (0 = fully adaptive).
        self%frozen = .true. ; self%nsub = self%nsub_best
    end subroutine learn_kexpm_sp

    subroutine print_kexpm_sp(self, reset_counters, verbose)
        implicit none(type, external)
        class(kexpm_sp_metadata), intent(inout) :: self
        logical, optional, intent(in) :: reset_counters
        !! Reset the counters (not the learned settings) after printing?
        logical, optional, intent(in) :: verbose
        !! Also print the internal state of the step-count search?
        character(len=*), parameter :: this_procedure = 'print_kexpm_sp'
        character(len=128) :: msg

        write(msg,'(A,I0,A,I0,A,I0)') 'Calls: ', self%n_calls, ', matvecs: ', self%n_matvec, &
                                    & ', sub-steps: ', self%n_substeps
        call log_message(msg, this_module, this_procedure)
        write(msg,'(A,I0,A,I0)') 'Rejected step sizes: ', self%n_rejections, ', fallbacks: ', self%n_fallbacks
        call log_message(msg, this_module, this_procedure)
        if (self%learned) then
            write(msg,'(A,I0,A,L1,A,I0)') 'Learned uniform sub-steps (0 = adaptive): ', self%nsub, &
                                        & ', frozen: ', self%frozen, ', best cost: ', self%mv_best
        else
            msg = 'No settings learned yet.'
        endif
        call log_message(msg, this_module, this_procedure)
        if (optval(verbose, .false.)) then
            write(msg,'(A,I0,A,I0)') 'Search: nsub0= ', self%nsub0, ', nsub_best= ', self%nsub_best
            call log_message(msg, this_module, this_procedure)
        endif
        if (self%converged) then
            call log_message('Status: CONVERGED', this_module, this_procedure)
        else
            call log_message('Status: NOT CONVERGED', this_module, this_procedure)
        endif
        if (optval(reset_counters, .false.)) call self%reset()
    end subroutine print_kexpm_sp

    subroutine reset_kexpm_sp(self)
        !! Resets the counters. The learned settings are kept (use `forget` to clear them).
        implicit none(type, external)
        class(kexpm_sp_metadata), intent(inout) :: self
        self%n_calls = 0 ; self%n_matvec = 0 ; self%n_substeps = 0
        self%n_rejections = 0 ; self%n_fallbacks = 0
        self%info = 0 ; self%converged = .false.
    end subroutine reset_kexpm_sp

    subroutine forget_kexpm_sp(self)
        !! Clears the learned settings: the next call is a cold, fully adaptive one.
        implicit none(type, external)
        class(kexpm_sp_metadata), intent(inout) :: self
        self%learned = .false. ; self%frozen = .false.
        self%nsub = 0 ; self%nsub_best = 0 ; self%nsub0 = 0 ; self%mv_best = 0
        self%tau_learned = 0.0_sp ; self%tol_learned = 0.0_sp ; self%kdim_learned = 0
    end subroutine forget_kexpm_sp

    !------------------------------------------------------------------
    !-----     Metadata of kexpm_substep : learning and printing     -----
    !------------------------------------------------------------------

    subroutine learn_kexpm_dp(self, used_uniform, n_uniform, nsub, mv, fallback, tune, t_end, tol_rel, m)
        !! Updates the learned number of uniform sub-steps after a successful call.
        !!
        !! The cost (matvecs) of N uniform sub-steps increases with N, while too small an N makes
        !! the error check fail (safety net). The cheapest setting is thus the smallest N that does
        !! not need the safety net. The search is:
        !!
        !! - Cold (adaptive) call: record its cost and number of sub-steps `nsub0`; next call
        !!   tries `nsub0 - 1` uniform sub-steps.
        !! - Uniform call, no fallback, cheaper than the best so far: keep it, try one sub-step less.
        !! - Uniform call with fallback, or not cheaper: stop and freeze on the best setting
        !!   (`nsub = 0`, i.e. fully adaptive, if no uniform setting was cheaper).
        !! - Fallback while frozen (e.g. the operator became stiffer): raise the count to
        !!   `max(nsub+1, sub-steps taken by this call)` and, with `tune`, restart the descent.
        implicit none(type, external)
        class(kexpm_dp_metadata), intent(inout) :: self
        logical, intent(in) :: used_uniform, fallback, tune
        integer, intent(in) :: n_uniform, nsub, mv, m
        real(dp), intent(in) :: t_end, tol_rel

        if (.not. self%learned) then
            ! Cold adaptive call.
            self%learned = .true.
            self%tau_learned = t_end ; self%tol_learned = tol_rel ; self%kdim_learned = m
            self%nsub0 = nsub ; self%nsub_best = 0 ; self%mv_best = mv
            if (tune .and. nsub > 1) then
                self%nsub = nsub - 1 ; self%frozen = .false.
            else
                self%nsub = nsub ; self%frozen = .true.
            endif
            return
        endif

        if (.not. used_uniform) return

        if (self%frozen) then
            if (fallback) then
                ! The frozen count became too small: restart the descent from a safe count.
                self%nsub = max(n_uniform + 1, nsub)
                self%nsub_best = 0 ; self%mv_best = huge(1)
                self%frozen = .not. tune
            endif
            return
        endif

        if ((.not. fallback) .and. mv < self%mv_best) then
            self%nsub_best = n_uniform ; self%mv_best = mv
            if (n_uniform > 1) then
                self%nsub = n_uniform - 1
                return
            endif
        endif

        ! Search over: freeze on the best setting (0 = fully adaptive).
        self%frozen = .true. ; self%nsub = self%nsub_best
    end subroutine learn_kexpm_dp

    subroutine print_kexpm_dp(self, reset_counters, verbose)
        implicit none(type, external)
        class(kexpm_dp_metadata), intent(inout) :: self
        logical, optional, intent(in) :: reset_counters
        !! Reset the counters (not the learned settings) after printing?
        logical, optional, intent(in) :: verbose
        !! Also print the internal state of the step-count search?
        character(len=*), parameter :: this_procedure = 'print_kexpm_dp'
        character(len=128) :: msg

        write(msg,'(A,I0,A,I0,A,I0)') 'Calls: ', self%n_calls, ', matvecs: ', self%n_matvec, &
                                    & ', sub-steps: ', self%n_substeps
        call log_message(msg, this_module, this_procedure)
        write(msg,'(A,I0,A,I0)') 'Rejected step sizes: ', self%n_rejections, ', fallbacks: ', self%n_fallbacks
        call log_message(msg, this_module, this_procedure)
        if (self%learned) then
            write(msg,'(A,I0,A,L1,A,I0)') 'Learned uniform sub-steps (0 = adaptive): ', self%nsub, &
                                        & ', frozen: ', self%frozen, ', best cost: ', self%mv_best
        else
            msg = 'No settings learned yet.'
        endif
        call log_message(msg, this_module, this_procedure)
        if (optval(verbose, .false.)) then
            write(msg,'(A,I0,A,I0)') 'Search: nsub0= ', self%nsub0, ', nsub_best= ', self%nsub_best
            call log_message(msg, this_module, this_procedure)
        endif
        if (self%converged) then
            call log_message('Status: CONVERGED', this_module, this_procedure)
        else
            call log_message('Status: NOT CONVERGED', this_module, this_procedure)
        endif
        if (optval(reset_counters, .false.)) call self%reset()
    end subroutine print_kexpm_dp

    subroutine reset_kexpm_dp(self)
        !! Resets the counters. The learned settings are kept (use `forget` to clear them).
        implicit none(type, external)
        class(kexpm_dp_metadata), intent(inout) :: self
        self%n_calls = 0 ; self%n_matvec = 0 ; self%n_substeps = 0
        self%n_rejections = 0 ; self%n_fallbacks = 0
        self%info = 0 ; self%converged = .false.
    end subroutine reset_kexpm_dp

    subroutine forget_kexpm_dp(self)
        !! Clears the learned settings: the next call is a cold, fully adaptive one.
        implicit none(type, external)
        class(kexpm_dp_metadata), intent(inout) :: self
        self%learned = .false. ; self%frozen = .false.
        self%nsub = 0 ; self%nsub_best = 0 ; self%nsub0 = 0 ; self%mv_best = 0
        self%tau_learned = 0.0_dp ; self%tol_learned = 0.0_dp ; self%kdim_learned = 0
    end subroutine forget_kexpm_dp

    !--------------------------------------------
    !-----     Invalid-input diagnostics     -----
    !--------------------------------------------

    subroutine log_kexpm_invalid_input(info, procedure)
        !! Logs the reason for an invalid-input `info` code (the caller returns immediately).
        implicit none(type, external)
        integer, intent(in) :: info
        character(len=*), intent(in) :: procedure
        character(len=128) :: msg
        select case (info)
        case (-3)
            msg = 'Invalid workspace X (fewer than 2*blksize vectors or not a multiple of blksize).'
        case (-4)
            msg = 'Invalid Krylov dimension (kdim < 1).'
        case (-5)
            msg = 'Invalid tolerance (tol < 0).'
        case (-6)
            msg = 'Invalid block size (B is empty or size(C) /= size(B)).'
        case (-7)
            msg = 'Invalid options (max_rejections < 1 or max_substeps < 1).'
        case default
            write(msg, '(A,I0)') 'Invalid input, info = ', info
        end select
        call log_warning(trim(msg), this_module, procedure)
    end subroutine log_kexpm_invalid_input

end module LightKrylov_ExpmLib