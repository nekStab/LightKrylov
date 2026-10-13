module LightKrylov_ExpmLib
    !!  This module implements the evaluation of the "matrix-exponential times vector" procedure
    !!  using Krylov methods.

    ! Iso Fortran.
    use, intrinsic :: iso_fortran_env, only: output_unit

    ! Fortran standard library.
    use stdlib_optval, only: optval
    use stdlib_linalg, only: mnorm, norm, expm

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

    public :: abstract_exptA_rsp
    public :: abstract_exptA_rdp
    public :: abstract_exptA_csp
    public :: abstract_exptA_cdp
    public :: kexpm
    public :: kexpm_work_size
    public :: krylov_exptA
    public :: krylov_exptA_rsp
    public :: krylov_exptA_rdp
    public :: krylov_exptA_csp
    public :: krylov_exptA_cdp

    abstract interface
        subroutine abstract_exptA_rsp(vec_out, A, vec_in, tau, info, trans)
            !! Abstract interface to define the matrix exponential-vector product.
            import sp
            import abstract_vector_rsp
            import abstract_linop_rsp
            implicit none(type, external)
            class(abstract_vector_rsp), intent(inout) :: vec_out
            !! Solution vector. Allocated by the caller; its content is overwritten.
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
            class(abstract_vector_rdp), intent(inout) :: vec_out
            !! Solution vector. Allocated by the caller; its content is overwritten.
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
            class(abstract_vector_csp), intent(inout) :: vec_out
            !! Solution vector. Allocated by the caller; its content is overwritten.
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
            class(abstract_vector_cdp), intent(inout) :: vec_out
            !! Solution vector. Allocated by the caller; its content is overwritten.
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
        !!      call kexpm(c, A, b, tau, tol, info [, trans] [, kdim] [, work])
        !!  ```
        !!
        !!  ### Arguments
        !!
        !!  - `c`   :   Output vector (or vectors), allocated by the caller. Its content is
        !!              overwritten. It is an `intent(inout)` argument.
        !!
        !!  - `A`   :   Linear operator to be exponentiated. It is an `intent(inout)` argument.
        !!
        !!  - `b`   :   Vector to be multiplied by \( \exp(\tau A) \). It is an `intent(in)` argument.
        !!
        !!  - `tau` :   `real` (singe or double) time over which the matrix exponential needs to
        !!              be computed. It is an `intent(in)` argument.
        !!
        !!  - `info`    :   `integer` Information flag.
        !!
        !!  - `trans` (optional)    :   Whether \( A \) or \( A^H \) is being used.
        !!                              (default `trans=.false.`)
        !!
        !!  - `kdim` (optional)     :   Dimension of the Krylov subspace used in the Arnoldi method.
        !!
        !!  - `work` (optional)     :   Caller-provided workspace of at least `kexpm_work_size(kdim [, size(B)])`
        !!                              vectors. It is shaped with `init_like`, its content is undefined on
        !!                              exit, and it is never freed. If absent, it is allocated internally.
        module procedure kexpm_vec_rsp
        module procedure kexpm_mat_rsp
        module procedure kexpm_vec_rdp
        module procedure kexpm_mat_rdp
        module procedure kexpm_vec_csp
        module procedure kexpm_mat_csp
        module procedure kexpm_vec_cdp
        module procedure kexpm_mat_cdp
    end interface

    interface krylov_exptA
        !!  ### Description
        !!
        !!  Utility function to evaluate the matrix-exponential times vector.
        !!
        !!  ### Syntax
        !!
        !!  ```fortran
        !!      call k_exptA(vec_out, A, vec_in, tau, info, trans)
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
        !!  - `info`    :   Information flag.
        !!
        !!  - `trans`   :   Whether \( A \) or \( A^H \) is being used.
        module procedure krylov_exptA_rsp
        module procedure krylov_exptA_rdp
        module procedure krylov_exptA_csp
        module procedure krylov_exptA_cdp
    end interface

contains

    pure integer function kexpm_work_size(kdim, blksize) result(nwork)
        !! Number of vectors required in the `work` argument of `kexpm`: `kdim + 1` for a single
        !! vector, and `p*(kdim*p + 1)` for a block of `p = size(B)` vectors (the block variant
        !! performs up to `kdim*p` block-Arnoldi steps).
        integer, optional, intent(in) :: kdim
        !! Maximum dimension of the Krylov subspace (default 100, as in `kexpm`).
        integer, optional, intent(in) :: blksize
        !! Number of input vectors, i.e. `size(B)` for the block variant (default 1).
        integer :: k, p
        k = optval(kdim, 100) ; p = optval(blksize, 1)
        nwork = p*(k*p + 1)
    end function kexpm_work_size

    subroutine kexpm_vec_rsp(c, A, b, tau, tol, info, trans, kdim, work)
        !! Best approximation of \( \exp(\tau \mathbf{A}) \mathbf{b} \) in the computed Krylov subspace
        implicit none(type, external)
        class(abstract_vector_rsp), intent(inout) :: c
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
        !! Maximum size of the Krylov subspace.
        class(abstract_vector_rsp), optional, intent(inout), target :: work(:)
        !! Workspace (optional) of at least `kexpm_work_size(kdim [, blksize])` vectors.

        ! ----- Internal variables -----
        character(len=*), parameter :: this_procedure = 'kexpm_vec_rsp'
        integer, parameter :: kmax = 100
        integer :: k, km, kp, nk
        ! Arnoldi factorization (Krylov basis points into `work` or into `wloc`).
        class(abstract_vector_rsp), allocatable, target :: wloc(:)
        class(abstract_vector_rsp), pointer :: X(:)
        real(sp), allocatable :: H(:, :)
        integer :: nwork
        ! Normaliztion & temporary arrays.
        real(sp), allocatable :: E(:, :)
        real(sp) :: err_est, beta
        ! Optional arguments.
        logical :: transpose
        integer :: nsteps, iostat
        character(len=256) :: msg

        ! Deals with optional args.
        transpose = optval(trans, .false.)
        nsteps    = optval(kdim, kmax)
        nk        = nsteps

        info = 0

        ! Workspace: user-provided, or allocated here and released at exit.
        nwork = kexpm_work_size(nk)
        if (present(work)) then
            call prepare_workspace(work, nwork, b, this_module, this_procedure)
            X => work(1:nwork)
        else
            allocate(wloc(nwork), mold=b, stat=iostat, errmsg=msg)
            call check_allocation(iostat, msg, this_module, this_procedure)
            call init_like_basis(wloc, b)
            X => wloc
        endif
        allocate(H(nk+1, nk+1), E(nk+1, nk+1), source=zero_rsp, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)
        ! inout arguments
        call c%init_like(b)

        ! Normalize input vector and initialize Krylov subspace.
        beta = b%norm()
        if (beta == 0.0_sp) then
            ! Input is zero => Output is zero.
            call c%zero()
            err_est = 0.0_sp
            kp = 1
        else
            call zero_basis(X)
            call X(1)%add(b) ; call X(1)%scal(one_rsp / beta)

            expm_arnoldi: do k = 1, nk
                km = k-1 ; kp = k+1

                ! Reset work arrays.
                E = 0.0_sp

                ! Compute k-th step Arnoldi factorization.
                call arnoldi(A, X, H, info, kstart=k, kend=k, transpose=transpose)
                call check_info(info, 'arnoldi', this_module, this_procedure)

                ! Compute approximation.
                if (info == k) then
                    ! Arnoldi breakdown, do not consider extended matrix.
                    kp = k
                    info = -2
                endif

                ! Compute the (dense) matrix exponential of the extended Hessenberg matrix.
                E(:kp, :kp) = expm(tau*H(:kp, :kp))

                ! Cheap error esimate (this actually is the magnitude of the included correction
                ! and thus is very conservative).
                err_est = merge(0.0_sp, abs(E(kp, 1)*beta), info==k)

                ! Check convergence.
                if (err_est <= tol) exit expm_arnoldi

            enddo expm_arnoldi

            ! Project back into original space.
            ! c = beta X(:kp) E(:kp, 1), computed directly into c (no temporary vector).
            call linear_combination(c, X(:kp), beta*E(:kp, 1))
        endif

        if (err_est <= tol) then
            info = kp
            write(msg,'(A,I0,2(A,E9.2))') 'Converged. kp= ', kp, ', err_est= ', err_est, ', tol= ', tol
            call log_information(msg, this_module, this_procedure)
        else
            write(msg,'(A,I0,2(A,E9.2))') 'Not converged. kp= ', nk+1, ', err_est= ', err_est, ', tol= ', tol
            call log_information(msg, this_module, this_procedure)
            info = -1
        endif

        !> Cleanup (a caller-provided workspace is never freed here).
        if (.not. present(work)) call free_basis(wloc)
        nullify(X)
    end subroutine kexpm_vec_rsp

    subroutine kexpm_mat_rsp(C, A, B, tau, tol, info, trans, kdim, work)
        !! Best Krylov approximation of \( \mathbf{C} = \exp(\tau \mathbf{A}) \mathbf{B} \) in the computed Krylov subspace.
        implicit none(type, external)
        class(abstract_vector_rsp), intent(inout) :: C(:)
        class(abstract_linop_rsp), intent(inout) :: A
        !! Linear operator to be exponentiated.
        class(abstract_vector_rsp), intent(in) :: B(:)
        !! Input matrix on which to apply \( \exp(\tau \mathbf{A}) \).
        real(sp), intent(in) :: tau
        !! Time horizon for the exponentiation.
        real(sp), intent(in) :: tol
        !! Solution toleance based on error estimates.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: trans
        !! Use transpose ?
        integer, optional, intent(in) :: kdim
        !! Maximum size of the Krylov subspace.
        class(abstract_vector_rsp), optional, intent(inout), target :: work(:)
        !! Workspace (optional) of at least `kexpm_work_size(kdim [, blksize])` vectors.

        ! ----- Internal variables -----
        character(len=*), parameter :: this_procedure = 'kexpm_mat_rsp'
        integer, parameter :: kmax = 100
        integer :: k, p, kpm, kp, kpp, nk
        ! Block-Arnoldi factorization (Krylov basis points into `work` or into `wloc`).
        class(abstract_vector_rsp), allocatable, target :: wloc(:)
        class(abstract_vector_rsp), pointer :: X(:)
        integer :: nwork
        real(sp), allocatable :: H(:, :)
        ! Normalization & temporary arrays.
        real(sp), allocatable :: R(:, :), E(:, :)
        integer, allocatable :: perm(:), ptrans(:)
        real(sp) :: err_est
        ! Optional arguments.
        logical :: transpose
        integer :: nsteps, iostat
        character(len=256) :: msg

        ! Determine block size.
        p = size(B)

        ! Deals with the optional args.
        transpose = optval(trans, .false.)
        nsteps    = optval(kdim, kmax)
        nk        = nsteps*p

        info = 0

        ! Allocate arrays.
        allocate(R(p, p), H(p*(nk+1), p*(nk+1)), E(p*(nk+1), p*(nk+1)), &
                 source=zero_rsp, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)

        allocate(perm(p), ptrans(p), source=0, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)

        ! Workspace: user-provided, or allocated here and released at exit.
        nwork = kexpm_work_size(nsteps, p)
        if (present(work)) then
            call prepare_workspace(work, nwork, B(1), this_module, this_procedure)
            X => work(1:nwork)
        else
            allocate(wloc(nwork), mold=B(1), stat=iostat, errmsg=msg)
            call check_allocation(iostat, msg, this_module, this_procedure)
            call init_like_basis(wloc, B(1))
            X => wloc
        endif

        ! Initialize Krylov subspace.
        call copy(X(:p), B)

        ! Normalize input matrix and initialize Krylov subspace.
        call qr(X(:p), R, perm, info) ; call permcols(R, invperm(perm))

        if (mnorm(R, "fro") == 0.0_sp) then
            call zero_basis(C)
            ! Input matrix is zero.
            err_est = 0.0_sp ; k = 0 ; kpp = p
        else
            H = 0.0_sp

            expm_arnoldi: do k = 1, nk
                ! Set counters.
                kpm = (k-1)*p ; kp = kpm + p ; kpp = kp + p

                ! Reset working arrays.
                E = 0.0_sp

                ! Compute the k-th step of the Arnoldi factorization.
                call arnoldi(A, X, H, info, kstart=k, kend=k, transpose=transpose, blksize=p)
                call check_info(info, 'arnoldi', this_module, this_procedure)

                if (info == kp) then
                    ! Arnoldi breakdown. Do not consider extended matrix.
                    kpp = kp
                    info = -2
                endif

                ! Compute the (dense) matrix exponential of the extended Hessenberg matrix.
                E(:kpp, :kpp) = expm(tau*H(:kpp, :kpp))

                ! Cheap error estimate.
                if (info == kp) then
                    ! Approximation is exact.
                    err_est = 0.0_sp
                else
                    err_est = norm(matmul(E(kp+1:kpp, :p), R(:p, :p)), 2)
                endif

                if (err_est <= tol) exit expm_arnoldi

            enddo expm_arnoldi

            ! Project back to original space.
            ! C = (X E) R = X (E R), computed directly into C (no temporary basis).
            call linear_combination(C, X(:kpp), matmul(E(:kpp, :p), R(:p, :p)))
        endif

        if (err_est <= tol) then
            info = kpp
            write(msg,'(A,I0,2(A,E9.2))') 'Converged. kp= ', kpp, ', err_est= ', err_est, ', tol= ', tol
            call log_information(msg, this_module, this_procedure)
        else
            write(msg,'(A,I0,2(A,E9.2))') 'Not converged. kp= ', kpp, ', err_est= ', err_est, ', tol= ', tol
            call log_information(msg, this_module, this_procedure)
            info = -1
        endif

        ! Cleanup (a caller-provided workspace is never freed here).
        if (.not. present(work)) call free_basis(wloc)
        nullify(X)
    end subroutine kexpm_mat_rsp

    subroutine krylov_exptA_rsp(vec_out, A, vec_in, tau, info, trans)
        !! Wrapper for the Krylov-based evaluation of the action of the matrix exponential operator
        !! on a vector that conforms to the `abstract_exptA_rsp` interface.
        implicit none(type, external)
        class(abstract_vector_rsp), intent(inout) :: vec_out
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
        integer :: kdim

        tol = atol_sp
        kdim = 30

        call kexpm(vec_out, A, vec_in, tau, tol, info, trans=trans, kdim=kdim)
        call check_info(info, 'kexpm', this_module, this_procedure)
    end subroutine krylov_exptA_rsp

    subroutine kexpm_vec_rdp(c, A, b, tau, tol, info, trans, kdim, work)
        !! Best approximation of \( \exp(\tau \mathbf{A}) \mathbf{b} \) in the computed Krylov subspace
        implicit none(type, external)
        class(abstract_vector_rdp), intent(inout) :: c
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
        !! Maximum size of the Krylov subspace.
        class(abstract_vector_rdp), optional, intent(inout), target :: work(:)
        !! Workspace (optional) of at least `kexpm_work_size(kdim [, blksize])` vectors.

        ! ----- Internal variables -----
        character(len=*), parameter :: this_procedure = 'kexpm_vec_rdp'
        integer, parameter :: kmax = 100
        integer :: k, km, kp, nk
        ! Arnoldi factorization (Krylov basis points into `work` or into `wloc`).
        class(abstract_vector_rdp), allocatable, target :: wloc(:)
        class(abstract_vector_rdp), pointer :: X(:)
        real(dp), allocatable :: H(:, :)
        integer :: nwork
        ! Normaliztion & temporary arrays.
        real(dp), allocatable :: E(:, :)
        real(dp) :: err_est, beta
        ! Optional arguments.
        logical :: transpose
        integer :: nsteps, iostat
        character(len=256) :: msg

        ! Deals with optional args.
        transpose = optval(trans, .false.)
        nsteps    = optval(kdim, kmax)
        nk        = nsteps

        info = 0

        ! Workspace: user-provided, or allocated here and released at exit.
        nwork = kexpm_work_size(nk)
        if (present(work)) then
            call prepare_workspace(work, nwork, b, this_module, this_procedure)
            X => work(1:nwork)
        else
            allocate(wloc(nwork), mold=b, stat=iostat, errmsg=msg)
            call check_allocation(iostat, msg, this_module, this_procedure)
            call init_like_basis(wloc, b)
            X => wloc
        endif
        allocate(H(nk+1, nk+1), E(nk+1, nk+1), source=zero_rdp, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)
        ! inout arguments
        call c%init_like(b)

        ! Normalize input vector and initialize Krylov subspace.
        beta = b%norm()
        if (beta == 0.0_dp) then
            ! Input is zero => Output is zero.
            call c%zero()
            err_est = 0.0_dp
            kp = 1
        else
            call zero_basis(X)
            call X(1)%add(b) ; call X(1)%scal(one_rdp / beta)

            expm_arnoldi: do k = 1, nk
                km = k-1 ; kp = k+1

                ! Reset work arrays.
                E = 0.0_dp

                ! Compute k-th step Arnoldi factorization.
                call arnoldi(A, X, H, info, kstart=k, kend=k, transpose=transpose)
                call check_info(info, 'arnoldi', this_module, this_procedure)

                ! Compute approximation.
                if (info == k) then
                    ! Arnoldi breakdown, do not consider extended matrix.
                    kp = k
                    info = -2
                endif

                ! Compute the (dense) matrix exponential of the extended Hessenberg matrix.
                E(:kp, :kp) = expm(tau*H(:kp, :kp))

                ! Cheap error esimate (this actually is the magnitude of the included correction
                ! and thus is very conservative).
                err_est = merge(0.0_dp, abs(E(kp, 1)*beta), info==k)

                ! Check convergence.
                if (err_est <= tol) exit expm_arnoldi

            enddo expm_arnoldi

            ! Project back into original space.
            ! c = beta X(:kp) E(:kp, 1), computed directly into c (no temporary vector).
            call linear_combination(c, X(:kp), beta*E(:kp, 1))
        endif

        if (err_est <= tol) then
            info = kp
            write(msg,'(A,I0,2(A,E9.2))') 'Converged. kp= ', kp, ', err_est= ', err_est, ', tol= ', tol
            call log_information(msg, this_module, this_procedure)
        else
            write(msg,'(A,I0,2(A,E9.2))') 'Not converged. kp= ', nk+1, ', err_est= ', err_est, ', tol= ', tol
            call log_information(msg, this_module, this_procedure)
            info = -1
        endif

        !> Cleanup (a caller-provided workspace is never freed here).
        if (.not. present(work)) call free_basis(wloc)
        nullify(X)
    end subroutine kexpm_vec_rdp

    subroutine kexpm_mat_rdp(C, A, B, tau, tol, info, trans, kdim, work)
        !! Best Krylov approximation of \( \mathbf{C} = \exp(\tau \mathbf{A}) \mathbf{B} \) in the computed Krylov subspace.
        implicit none(type, external)
        class(abstract_vector_rdp), intent(inout) :: C(:)
        class(abstract_linop_rdp), intent(inout) :: A
        !! Linear operator to be exponentiated.
        class(abstract_vector_rdp), intent(in) :: B(:)
        !! Input matrix on which to apply \( \exp(\tau \mathbf{A}) \).
        real(dp), intent(in) :: tau
        !! Time horizon for the exponentiation.
        real(dp), intent(in) :: tol
        !! Solution toleance based on error estimates.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: trans
        !! Use transpose ?
        integer, optional, intent(in) :: kdim
        !! Maximum size of the Krylov subspace.
        class(abstract_vector_rdp), optional, intent(inout), target :: work(:)
        !! Workspace (optional) of at least `kexpm_work_size(kdim [, blksize])` vectors.

        ! ----- Internal variables -----
        character(len=*), parameter :: this_procedure = 'kexpm_mat_rdp'
        integer, parameter :: kmax = 100
        integer :: k, p, kpm, kp, kpp, nk
        ! Block-Arnoldi factorization (Krylov basis points into `work` or into `wloc`).
        class(abstract_vector_rdp), allocatable, target :: wloc(:)
        class(abstract_vector_rdp), pointer :: X(:)
        integer :: nwork
        real(dp), allocatable :: H(:, :)
        ! Normalization & temporary arrays.
        real(dp), allocatable :: R(:, :), E(:, :)
        integer, allocatable :: perm(:), ptrans(:)
        real(dp) :: err_est
        ! Optional arguments.
        logical :: transpose
        integer :: nsteps, iostat
        character(len=256) :: msg

        ! Determine block size.
        p = size(B)

        ! Deals with the optional args.
        transpose = optval(trans, .false.)
        nsteps    = optval(kdim, kmax)
        nk        = nsteps*p

        info = 0

        ! Allocate arrays.
        allocate(R(p, p), H(p*(nk+1), p*(nk+1)), E(p*(nk+1), p*(nk+1)), &
                 source=zero_rdp, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)

        allocate(perm(p), ptrans(p), source=0, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)

        ! Workspace: user-provided, or allocated here and released at exit.
        nwork = kexpm_work_size(nsteps, p)
        if (present(work)) then
            call prepare_workspace(work, nwork, B(1), this_module, this_procedure)
            X => work(1:nwork)
        else
            allocate(wloc(nwork), mold=B(1), stat=iostat, errmsg=msg)
            call check_allocation(iostat, msg, this_module, this_procedure)
            call init_like_basis(wloc, B(1))
            X => wloc
        endif

        ! Initialize Krylov subspace.
        call copy(X(:p), B)

        ! Normalize input matrix and initialize Krylov subspace.
        call qr(X(:p), R, perm, info) ; call permcols(R, invperm(perm))

        if (mnorm(R, "fro") == 0.0_dp) then
            call zero_basis(C)
            ! Input matrix is zero.
            err_est = 0.0_dp ; k = 0 ; kpp = p
        else
            H = 0.0_dp

            expm_arnoldi: do k = 1, nk
                ! Set counters.
                kpm = (k-1)*p ; kp = kpm + p ; kpp = kp + p

                ! Reset working arrays.
                E = 0.0_dp

                ! Compute the k-th step of the Arnoldi factorization.
                call arnoldi(A, X, H, info, kstart=k, kend=k, transpose=transpose, blksize=p)
                call check_info(info, 'arnoldi', this_module, this_procedure)

                if (info == kp) then
                    ! Arnoldi breakdown. Do not consider extended matrix.
                    kpp = kp
                    info = -2
                endif

                ! Compute the (dense) matrix exponential of the extended Hessenberg matrix.
                E(:kpp, :kpp) = expm(tau*H(:kpp, :kpp))

                ! Cheap error estimate.
                if (info == kp) then
                    ! Approximation is exact.
                    err_est = 0.0_dp
                else
                    err_est = norm(matmul(E(kp+1:kpp, :p), R(:p, :p)), 2)
                endif

                if (err_est <= tol) exit expm_arnoldi

            enddo expm_arnoldi

            ! Project back to original space.
            ! C = (X E) R = X (E R), computed directly into C (no temporary basis).
            call linear_combination(C, X(:kpp), matmul(E(:kpp, :p), R(:p, :p)))
        endif

        if (err_est <= tol) then
            info = kpp
            write(msg,'(A,I0,2(A,E9.2))') 'Converged. kp= ', kpp, ', err_est= ', err_est, ', tol= ', tol
            call log_information(msg, this_module, this_procedure)
        else
            write(msg,'(A,I0,2(A,E9.2))') 'Not converged. kp= ', kpp, ', err_est= ', err_est, ', tol= ', tol
            call log_information(msg, this_module, this_procedure)
            info = -1
        endif

        ! Cleanup (a caller-provided workspace is never freed here).
        if (.not. present(work)) call free_basis(wloc)
        nullify(X)
    end subroutine kexpm_mat_rdp

    subroutine krylov_exptA_rdp(vec_out, A, vec_in, tau, info, trans)
        !! Wrapper for the Krylov-based evaluation of the action of the matrix exponential operator
        !! on a vector that conforms to the `abstract_exptA_rdp` interface.
        implicit none(type, external)
        class(abstract_vector_rdp), intent(inout) :: vec_out
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
        integer :: kdim

        tol = atol_dp
        kdim = 30

        call kexpm(vec_out, A, vec_in, tau, tol, info, trans=trans, kdim=kdim)
        call check_info(info, 'kexpm', this_module, this_procedure)
    end subroutine krylov_exptA_rdp

    subroutine kexpm_vec_csp(c, A, b, tau, tol, info, trans, kdim, work)
        !! Best approximation of \( \exp(\tau \mathbf{A}) \mathbf{b} \) in the computed Krylov subspace
        implicit none(type, external)
        class(abstract_vector_csp), intent(inout) :: c
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
        !! Maximum size of the Krylov subspace.
        class(abstract_vector_csp), optional, intent(inout), target :: work(:)
        !! Workspace (optional) of at least `kexpm_work_size(kdim [, blksize])` vectors.

        ! ----- Internal variables -----
        character(len=*), parameter :: this_procedure = 'kexpm_vec_csp'
        integer, parameter :: kmax = 100
        integer :: k, km, kp, nk
        ! Arnoldi factorization (Krylov basis points into `work` or into `wloc`).
        class(abstract_vector_csp), allocatable, target :: wloc(:)
        class(abstract_vector_csp), pointer :: X(:)
        complex(sp), allocatable :: H(:, :)
        integer :: nwork
        ! Normaliztion & temporary arrays.
        complex(sp), allocatable :: E(:, :)
        real(sp) :: err_est, beta
        ! Optional arguments.
        logical :: transpose
        integer :: nsteps, iostat
        character(len=256) :: msg

        ! Deals with optional args.
        transpose = optval(trans, .false.)
        nsteps    = optval(kdim, kmax)
        nk        = nsteps

        info = 0

        ! Workspace: user-provided, or allocated here and released at exit.
        nwork = kexpm_work_size(nk)
        if (present(work)) then
            call prepare_workspace(work, nwork, b, this_module, this_procedure)
            X => work(1:nwork)
        else
            allocate(wloc(nwork), mold=b, stat=iostat, errmsg=msg)
            call check_allocation(iostat, msg, this_module, this_procedure)
            call init_like_basis(wloc, b)
            X => wloc
        endif
        allocate(H(nk+1, nk+1), E(nk+1, nk+1), source=zero_csp, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)
        ! inout arguments
        call c%init_like(b)

        ! Normalize input vector and initialize Krylov subspace.
        beta = b%norm()
        if (beta == 0.0_sp) then
            ! Input is zero => Output is zero.
            call c%zero()
            err_est = 0.0_sp
            kp = 1
        else
            call zero_basis(X)
            call X(1)%add(b) ; call X(1)%scal(one_csp / beta)

            expm_arnoldi: do k = 1, nk
                km = k-1 ; kp = k+1

                ! Reset work arrays.
                E = 0.0_sp

                ! Compute k-th step Arnoldi factorization.
                call arnoldi(A, X, H, info, kstart=k, kend=k, transpose=transpose)
                call check_info(info, 'arnoldi', this_module, this_procedure)

                ! Compute approximation.
                if (info == k) then
                    ! Arnoldi breakdown, do not consider extended matrix.
                    kp = k
                    info = -2
                endif

                ! Compute the (dense) matrix exponential of the extended Hessenberg matrix.
                E(:kp, :kp) = expm(tau*H(:kp, :kp))

                ! Cheap error esimate (this actually is the magnitude of the included correction
                ! and thus is very conservative).
                err_est = merge(0.0_sp, abs(E(kp, 1)*beta), info==k)

                ! Check convergence.
                if (err_est <= tol) exit expm_arnoldi

            enddo expm_arnoldi

            ! Project back into original space.
            ! c = beta X(:kp) E(:kp, 1), computed directly into c (no temporary vector).
            call linear_combination(c, X(:kp), beta*E(:kp, 1))
        endif

        if (err_est <= tol) then
            info = kp
            write(msg,'(A,I0,2(A,E9.2))') 'Converged. kp= ', kp, ', err_est= ', err_est, ', tol= ', tol
            call log_information(msg, this_module, this_procedure)
        else
            write(msg,'(A,I0,2(A,E9.2))') 'Not converged. kp= ', nk+1, ', err_est= ', err_est, ', tol= ', tol
            call log_information(msg, this_module, this_procedure)
            info = -1
        endif

        !> Cleanup (a caller-provided workspace is never freed here).
        if (.not. present(work)) call free_basis(wloc)
        nullify(X)
    end subroutine kexpm_vec_csp

    subroutine kexpm_mat_csp(C, A, B, tau, tol, info, trans, kdim, work)
        !! Best Krylov approximation of \( \mathbf{C} = \exp(\tau \mathbf{A}) \mathbf{B} \) in the computed Krylov subspace.
        implicit none(type, external)
        class(abstract_vector_csp), intent(inout) :: C(:)
        class(abstract_linop_csp), intent(inout) :: A
        !! Linear operator to be exponentiated.
        class(abstract_vector_csp), intent(in) :: B(:)
        !! Input matrix on which to apply \( \exp(\tau \mathbf{A}) \).
        real(sp), intent(in) :: tau
        !! Time horizon for the exponentiation.
        real(sp), intent(in) :: tol
        !! Solution toleance based on error estimates.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: trans
        !! Use transpose ?
        integer, optional, intent(in) :: kdim
        !! Maximum size of the Krylov subspace.
        class(abstract_vector_csp), optional, intent(inout), target :: work(:)
        !! Workspace (optional) of at least `kexpm_work_size(kdim [, blksize])` vectors.

        ! ----- Internal variables -----
        character(len=*), parameter :: this_procedure = 'kexpm_mat_csp'
        integer, parameter :: kmax = 100
        integer :: k, p, kpm, kp, kpp, nk
        ! Block-Arnoldi factorization (Krylov basis points into `work` or into `wloc`).
        class(abstract_vector_csp), allocatable, target :: wloc(:)
        class(abstract_vector_csp), pointer :: X(:)
        integer :: nwork
        complex(sp), allocatable :: H(:, :)
        ! Normalization & temporary arrays.
        complex(sp), allocatable :: R(:, :), E(:, :)
        integer, allocatable :: perm(:), ptrans(:)
        real(sp) :: err_est
        ! Optional arguments.
        logical :: transpose
        integer :: nsteps, iostat
        character(len=256) :: msg

        ! Determine block size.
        p = size(B)

        ! Deals with the optional args.
        transpose = optval(trans, .false.)
        nsteps    = optval(kdim, kmax)
        nk        = nsteps*p

        info = 0

        ! Allocate arrays.
        allocate(R(p, p), H(p*(nk+1), p*(nk+1)), E(p*(nk+1), p*(nk+1)), &
                 source=zero_csp, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)

        allocate(perm(p), ptrans(p), source=0, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)

        ! Workspace: user-provided, or allocated here and released at exit.
        nwork = kexpm_work_size(nsteps, p)
        if (present(work)) then
            call prepare_workspace(work, nwork, B(1), this_module, this_procedure)
            X => work(1:nwork)
        else
            allocate(wloc(nwork), mold=B(1), stat=iostat, errmsg=msg)
            call check_allocation(iostat, msg, this_module, this_procedure)
            call init_like_basis(wloc, B(1))
            X => wloc
        endif

        ! Initialize Krylov subspace.
        call copy(X(:p), B)

        ! Normalize input matrix and initialize Krylov subspace.
        call qr(X(:p), R, perm, info) ; call permcols(R, invperm(perm))

        if (mnorm(R, "fro") == 0.0_sp) then
            call zero_basis(C)
            ! Input matrix is zero.
            err_est = 0.0_sp ; k = 0 ; kpp = p
        else
            H = 0.0_sp

            expm_arnoldi: do k = 1, nk
                ! Set counters.
                kpm = (k-1)*p ; kp = kpm + p ; kpp = kp + p

                ! Reset working arrays.
                E = 0.0_sp

                ! Compute the k-th step of the Arnoldi factorization.
                call arnoldi(A, X, H, info, kstart=k, kend=k, transpose=transpose, blksize=p)
                call check_info(info, 'arnoldi', this_module, this_procedure)

                if (info == kp) then
                    ! Arnoldi breakdown. Do not consider extended matrix.
                    kpp = kp
                    info = -2
                endif

                ! Compute the (dense) matrix exponential of the extended Hessenberg matrix.
                E(:kpp, :kpp) = expm(tau*H(:kpp, :kpp))

                ! Cheap error estimate.
                if (info == kp) then
                    ! Approximation is exact.
                    err_est = 0.0_sp
                else
                    err_est = norm(matmul(E(kp+1:kpp, :p), R(:p, :p)), 2)
                endif

                if (err_est <= tol) exit expm_arnoldi

            enddo expm_arnoldi

            ! Project back to original space.
            ! C = (X E) R = X (E R), computed directly into C (no temporary basis).
            call linear_combination(C, X(:kpp), matmul(E(:kpp, :p), R(:p, :p)))
        endif

        if (err_est <= tol) then
            info = kpp
            write(msg,'(A,I0,2(A,E9.2))') 'Converged. kp= ', kpp, ', err_est= ', err_est, ', tol= ', tol
            call log_information(msg, this_module, this_procedure)
        else
            write(msg,'(A,I0,2(A,E9.2))') 'Not converged. kp= ', kpp, ', err_est= ', err_est, ', tol= ', tol
            call log_information(msg, this_module, this_procedure)
            info = -1
        endif

        ! Cleanup (a caller-provided workspace is never freed here).
        if (.not. present(work)) call free_basis(wloc)
        nullify(X)
    end subroutine kexpm_mat_csp

    subroutine krylov_exptA_csp(vec_out, A, vec_in, tau, info, trans)
        !! Wrapper for the Krylov-based evaluation of the action of the matrix exponential operator
        !! on a vector that conforms to the `abstract_exptA_csp` interface.
        implicit none(type, external)
        class(abstract_vector_csp), intent(inout) :: vec_out
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
        integer :: kdim

        tol = atol_sp
        kdim = 30

        call kexpm(vec_out, A, vec_in, tau, tol, info, trans=trans, kdim=kdim)
        call check_info(info, 'kexpm', this_module, this_procedure)
    end subroutine krylov_exptA_csp

    subroutine kexpm_vec_cdp(c, A, b, tau, tol, info, trans, kdim, work)
        !! Best approximation of \( \exp(\tau \mathbf{A}) \mathbf{b} \) in the computed Krylov subspace
        implicit none(type, external)
        class(abstract_vector_cdp), intent(inout) :: c
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
        !! Maximum size of the Krylov subspace.
        class(abstract_vector_cdp), optional, intent(inout), target :: work(:)
        !! Workspace (optional) of at least `kexpm_work_size(kdim [, blksize])` vectors.

        ! ----- Internal variables -----
        character(len=*), parameter :: this_procedure = 'kexpm_vec_cdp'
        integer, parameter :: kmax = 100
        integer :: k, km, kp, nk
        ! Arnoldi factorization (Krylov basis points into `work` or into `wloc`).
        class(abstract_vector_cdp), allocatable, target :: wloc(:)
        class(abstract_vector_cdp), pointer :: X(:)
        complex(dp), allocatable :: H(:, :)
        integer :: nwork
        ! Normaliztion & temporary arrays.
        complex(dp), allocatable :: E(:, :)
        real(dp) :: err_est, beta
        ! Optional arguments.
        logical :: transpose
        integer :: nsteps, iostat
        character(len=256) :: msg

        ! Deals with optional args.
        transpose = optval(trans, .false.)
        nsteps    = optval(kdim, kmax)
        nk        = nsteps

        info = 0

        ! Workspace: user-provided, or allocated here and released at exit.
        nwork = kexpm_work_size(nk)
        if (present(work)) then
            call prepare_workspace(work, nwork, b, this_module, this_procedure)
            X => work(1:nwork)
        else
            allocate(wloc(nwork), mold=b, stat=iostat, errmsg=msg)
            call check_allocation(iostat, msg, this_module, this_procedure)
            call init_like_basis(wloc, b)
            X => wloc
        endif
        allocate(H(nk+1, nk+1), E(nk+1, nk+1), source=zero_cdp, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)
        ! inout arguments
        call c%init_like(b)

        ! Normalize input vector and initialize Krylov subspace.
        beta = b%norm()
        if (beta == 0.0_dp) then
            ! Input is zero => Output is zero.
            call c%zero()
            err_est = 0.0_dp
            kp = 1
        else
            call zero_basis(X)
            call X(1)%add(b) ; call X(1)%scal(one_cdp / beta)

            expm_arnoldi: do k = 1, nk
                km = k-1 ; kp = k+1

                ! Reset work arrays.
                E = 0.0_dp

                ! Compute k-th step Arnoldi factorization.
                call arnoldi(A, X, H, info, kstart=k, kend=k, transpose=transpose)
                call check_info(info, 'arnoldi', this_module, this_procedure)

                ! Compute approximation.
                if (info == k) then
                    ! Arnoldi breakdown, do not consider extended matrix.
                    kp = k
                    info = -2
                endif

                ! Compute the (dense) matrix exponential of the extended Hessenberg matrix.
                E(:kp, :kp) = expm(tau*H(:kp, :kp))

                ! Cheap error esimate (this actually is the magnitude of the included correction
                ! and thus is very conservative).
                err_est = merge(0.0_dp, abs(E(kp, 1)*beta), info==k)

                ! Check convergence.
                if (err_est <= tol) exit expm_arnoldi

            enddo expm_arnoldi

            ! Project back into original space.
            ! c = beta X(:kp) E(:kp, 1), computed directly into c (no temporary vector).
            call linear_combination(c, X(:kp), beta*E(:kp, 1))
        endif

        if (err_est <= tol) then
            info = kp
            write(msg,'(A,I0,2(A,E9.2))') 'Converged. kp= ', kp, ', err_est= ', err_est, ', tol= ', tol
            call log_information(msg, this_module, this_procedure)
        else
            write(msg,'(A,I0,2(A,E9.2))') 'Not converged. kp= ', nk+1, ', err_est= ', err_est, ', tol= ', tol
            call log_information(msg, this_module, this_procedure)
            info = -1
        endif

        !> Cleanup (a caller-provided workspace is never freed here).
        if (.not. present(work)) call free_basis(wloc)
        nullify(X)
    end subroutine kexpm_vec_cdp

    subroutine kexpm_mat_cdp(C, A, B, tau, tol, info, trans, kdim, work)
        !! Best Krylov approximation of \( \mathbf{C} = \exp(\tau \mathbf{A}) \mathbf{B} \) in the computed Krylov subspace.
        implicit none(type, external)
        class(abstract_vector_cdp), intent(inout) :: C(:)
        class(abstract_linop_cdp), intent(inout) :: A
        !! Linear operator to be exponentiated.
        class(abstract_vector_cdp), intent(in) :: B(:)
        !! Input matrix on which to apply \( \exp(\tau \mathbf{A}) \).
        real(dp), intent(in) :: tau
        !! Time horizon for the exponentiation.
        real(dp), intent(in) :: tol
        !! Solution toleance based on error estimates.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: trans
        !! Use transpose ?
        integer, optional, intent(in) :: kdim
        !! Maximum size of the Krylov subspace.
        class(abstract_vector_cdp), optional, intent(inout), target :: work(:)
        !! Workspace (optional) of at least `kexpm_work_size(kdim [, blksize])` vectors.

        ! ----- Internal variables -----
        character(len=*), parameter :: this_procedure = 'kexpm_mat_cdp'
        integer, parameter :: kmax = 100
        integer :: k, p, kpm, kp, kpp, nk
        ! Block-Arnoldi factorization (Krylov basis points into `work` or into `wloc`).
        class(abstract_vector_cdp), allocatable, target :: wloc(:)
        class(abstract_vector_cdp), pointer :: X(:)
        integer :: nwork
        complex(dp), allocatable :: H(:, :)
        ! Normalization & temporary arrays.
        complex(dp), allocatable :: R(:, :), E(:, :)
        integer, allocatable :: perm(:), ptrans(:)
        real(dp) :: err_est
        ! Optional arguments.
        logical :: transpose
        integer :: nsteps, iostat
        character(len=256) :: msg

        ! Determine block size.
        p = size(B)

        ! Deals with the optional args.
        transpose = optval(trans, .false.)
        nsteps    = optval(kdim, kmax)
        nk        = nsteps*p

        info = 0

        ! Allocate arrays.
        allocate(R(p, p), H(p*(nk+1), p*(nk+1)), E(p*(nk+1), p*(nk+1)), &
                 source=zero_cdp, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)

        allocate(perm(p), ptrans(p), source=0, stat=iostat, errmsg=msg)
        call check_allocation(iostat, msg, this_module, this_procedure)

        ! Workspace: user-provided, or allocated here and released at exit.
        nwork = kexpm_work_size(nsteps, p)
        if (present(work)) then
            call prepare_workspace(work, nwork, B(1), this_module, this_procedure)
            X => work(1:nwork)
        else
            allocate(wloc(nwork), mold=B(1), stat=iostat, errmsg=msg)
            call check_allocation(iostat, msg, this_module, this_procedure)
            call init_like_basis(wloc, B(1))
            X => wloc
        endif

        ! Initialize Krylov subspace.
        call copy(X(:p), B)

        ! Normalize input matrix and initialize Krylov subspace.
        call qr(X(:p), R, perm, info) ; call permcols(R, invperm(perm))

        if (mnorm(R, "fro") == 0.0_dp) then
            call zero_basis(C)
            ! Input matrix is zero.
            err_est = 0.0_dp ; k = 0 ; kpp = p
        else
            H = 0.0_dp

            expm_arnoldi: do k = 1, nk
                ! Set counters.
                kpm = (k-1)*p ; kp = kpm + p ; kpp = kp + p

                ! Reset working arrays.
                E = 0.0_dp

                ! Compute the k-th step of the Arnoldi factorization.
                call arnoldi(A, X, H, info, kstart=k, kend=k, transpose=transpose, blksize=p)
                call check_info(info, 'arnoldi', this_module, this_procedure)

                if (info == kp) then
                    ! Arnoldi breakdown. Do not consider extended matrix.
                    kpp = kp
                    info = -2
                endif

                ! Compute the (dense) matrix exponential of the extended Hessenberg matrix.
                E(:kpp, :kpp) = expm(tau*H(:kpp, :kpp))

                ! Cheap error estimate.
                if (info == kp) then
                    ! Approximation is exact.
                    err_est = 0.0_dp
                else
                    err_est = norm(matmul(E(kp+1:kpp, :p), R(:p, :p)), 2)
                endif

                if (err_est <= tol) exit expm_arnoldi

            enddo expm_arnoldi

            ! Project back to original space.
            ! C = (X E) R = X (E R), computed directly into C (no temporary basis).
            call linear_combination(C, X(:kpp), matmul(E(:kpp, :p), R(:p, :p)))
        endif

        if (err_est <= tol) then
            info = kpp
            write(msg,'(A,I0,2(A,E9.2))') 'Converged. kp= ', kpp, ', err_est= ', err_est, ', tol= ', tol
            call log_information(msg, this_module, this_procedure)
        else
            write(msg,'(A,I0,2(A,E9.2))') 'Not converged. kp= ', kpp, ', err_est= ', err_est, ', tol= ', tol
            call log_information(msg, this_module, this_procedure)
            info = -1
        endif

        ! Cleanup (a caller-provided workspace is never freed here).
        if (.not. present(work)) call free_basis(wloc)
        nullify(X)
    end subroutine kexpm_mat_cdp

    subroutine krylov_exptA_cdp(vec_out, A, vec_in, tau, info, trans)
        !! Wrapper for the Krylov-based evaluation of the action of the matrix exponential operator
        !! on a vector that conforms to the `abstract_exptA_cdp` interface.
        implicit none(type, external)
        class(abstract_vector_cdp), intent(inout) :: vec_out
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
        integer :: kdim

        tol = atol_dp
        kdim = 30

        call kexpm(vec_out, A, vec_in, tau, tol, info, trans=trans, kdim=kdim)
        call check_info(info, 'kexpm', this_module, this_procedure)
    end subroutine krylov_exptA_cdp

end module LightKrylov_ExpmLib