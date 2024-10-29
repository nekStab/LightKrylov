module LightKrylov_ExpmLib
    !!  This module implements the evaluation of the "matrix-exponential times vector" procedure
    !!  using Krylov methods.

    ! Iso Fortran.
    use iso_fortran_env, only: output_unit

    ! Fortran standard library.
    use stdlib_optval, only: optval
    use stdlib_linalg, only: eye, inv

    ! LightKrylov.
    use LightKrylov_Constants
    use LightKrylov_Logger
    use LightKrylov_Utils
    use LightKrylov_AbstractVectors
    use LightKrylov_AbstractLinops
    use LightKrylov_BaseKrylov

    implicit none
    private
    
    character(len=128), parameter, private :: this_module= 'LightKrylov_ExpmLib'
    public :: abstract_exptA_rsp
    public :: abstract_exptA_rdp
    public :: abstract_exptA_csp
    public :: abstract_exptA_cdp
    public :: expm
    public :: kexpm
    public :: k_exptA

    abstract interface
        subroutine abstract_exptA_rsp(vec_out, A, vec_in, tau, info, trans)
            import sp
            import abstract_vector_rsp
            import abstract_linop_rsp
            !! Abstract interface to define the matrix exponential-vector product.
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
            import dp
            import abstract_vector_rdp
            import abstract_linop_rdp
            !! Abstract interface to define the matrix exponential-vector product.
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
            import sp
            import abstract_vector_csp
            import abstract_linop_csp
            !! Abstract interface to define the matrix exponential-vector product.
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
            import dp
            import abstract_vector_cdp
            import abstract_linop_cdp
            !! Abstract interface to define the matrix exponential-vector product.
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

    interface expm
        !!  ### Description
        !!
        !!  Evaluate the exponential of a dense matrix using Pade approximations.
        !!
        !!  ### Syntax
        !!
        !!  ```fortran
        !!      E = expm(A, order)
        !!  ```
        !!
        !!  ### Arguments
        !!
        !!  `E` : `real` or `complex` rank-2 array with \( E = \exp(A) \).
        !!
        !!  `A` : `real` or `complex` matrix that needs to be exponentiated.
        !!
        !!  `order` (optional) : Order of the Pade approximation. By default `order = 10`.
        module procedure expm_rsp
        module procedure expm_rdp
        module procedure expm_csp
        module procedure expm_cdp
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
        !!      call kexpm(c, A, b, tau, tol, info [, trans] [, kdim])
        !!  ```
        !!
        !!  ### Arguments
        !!
        !!  `c` : Output vector (or vectors). It is an `intent(out)` argument.
        !!
        !!  `A` : Linear operator to be exponentiated. It is an `intent(in)` argument.
        !!
        !!  `b` : Vector to be multiplied by \( \exp(\tau A) \). It is an `intent(in)` argument.
        !!
        !!  `tau` : `real` (singe or double) time over which the matrix exponential needs to
        !!          be computed. It is an `intent(in)` argument.
        !!
        !!  `info` : `integer` Information flag.
        !!
        !!  `trans` (optional) : Whether \( A \) or \( A^H \) is being used.
        !!                      (default `trans=.false.`)
        !!
        !!  `kdim` (optional) : Dimension of the Krylov subspace used in the Arnoldi method.
        module procedure kexpm_vec_rsp
        module procedure kexpm_mat_rsp
        module procedure kexpm_vec_rdp
        module procedure kexpm_mat_rdp
        module procedure kexpm_vec_csp
        module procedure kexpm_mat_csp
        module procedure kexpm_vec_cdp
        module procedure kexpm_mat_cdp
    end interface

    interface k_exptA
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
        !!  `vec_out` : Output vector.
        !!
        !!  `A` : Matrix to be exponentiated.
        !!
        !!  `vec_in` : Input vector.
        !!
        !!  `tau` : Integration time.
        !!
        !!  `info` : Information flag.
        !!
        !!  `trans` : Whether \( A \) or \( A^H \) is being used.
        module procedure k_exptA_rsp
        module procedure k_exptA_rdp
        module procedure k_exptA_csp
        module procedure k_exptA_cdp
    end interface

contains

    !--------------------------------------------
    !-----     DENSE MATRIX EXPONENTIAL     -----
    !--------------------------------------------

    function expm_rsp(A, order) result(E)
        real(sp), intent(in) :: A(:, :)
        !! Matrix to be exponentiated.
        real(sp) :: E(size(A, 1), size(A, 1))
        !! Output matrix E = exp(tA).
        integer, intent(in), optional :: order
        !! Order of the Pade approximation.

        !----- Internal variables -----
        real(sp), allocatable :: A2(:, :), Q(:, :), X(:, :)
        real(sp) :: a_norm, c
        integer :: n, ee, k, s
        logical :: p
        integer :: p_order

        ! Deal with optional args.
        p_order = optval(order, 10)

        n = size(A, 1)

        ! Allocate arrays.
        allocate(A2(n, n)) ; allocate(X(n, n)) ; allocate(Q(n, n))

        ! Compute the L-infinity norm.
        a_norm = norml(A)

        ! Determine scaling factor for the matrix.
        ee = int(log2(a_norm)) + 1
        s = max(0, ee+1)

        ! Scale the input matrix & initialize polynomial.
        A2 = A / 2.0_sp**s
        X = A2

        ! Initialize P & Q and add first step.
        c = 0.5_sp
        E = eye(n) ; E = E + c*A2

        Q = eye(n) ; Q = Q - c*A2

        ! Iteratively compute the Pade approximation.
        p = .true.
        do k = 2, p_order
            c = c*(p_order - k + 1) / (k * (2*p_order - k + 1))
            X = matmul(A2, X)
            E = E + c*X
            if (p) then
                Q = Q + c*X
            else
                Q = Q - c*X
            endif
            p = .not. p
        enddo

        E = matmul(inv(Q), E)
        do k = 1, s
            E = matmul(E, E)
        enddo

        return
    end function expm_rsp

    subroutine kexpm_vec_rsp(c, A, b, tau, tol, info, trans, kdim)
        class(abstract_vector_rsp), intent(out) :: c
        !! Best approximation of \( \exp(\tau \mathbf{A}) \mathbf{b} \) in the computed Krylov subspace
        class(abstract_linop_rsp), intent(in) :: A
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

        ! ----- Internal variables -----
        integer, parameter :: kmax = 100
        integer :: k, km, kp, nk
        ! Arnoldi factorization.
        class(abstract_vector_rsp), allocatable :: X(:)
        real(sp), allocatable :: H(:, :)
        ! Normaliztion & temporary arrays.
        real(sp), allocatable :: E(:, :)
        class(abstract_vector_rsp), allocatable :: Xwrk
        real(sp) :: err_est, beta
        ! Optional arguments.
        logical :: transpose
        integer :: nsteps
        character(len=256) :: msg
    
        ! Deals with optional args.
        transpose = optval(trans, .false.)
        nsteps    = optval(kdim, kmax)
        nk        = nsteps

        info = 0

        ! Allocate arrays.
        allocate(X(nk+1), source=b) ; allocate(H(nk+1, nk+1))
        allocate(E(nk+1, nk+1)) ; allocate(Xwrk, source=b)

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
            H = 0.0_sp

            expm_arnoldi: do k = 1, nk
                km = k-1 ; kp = k+1

                ! Reset work arrays.
                E = 0.0_sp

                ! Compute k-th step Arnoldi factorization.
                call arnoldi(A, X, H, info, kstart=k, kend=k, transpose=transpose)
                call check_info(info, 'arnoldi', module=this_module, procedure='kexpm_vec_rsp')

                ! Compute approximation.
                if (info == k) then
                    ! Arnoldi breakdown, do not consider extended matrix.
                    kp = k
                    info = -2
                endif

                ! Compute the (dense) matrix exponential of the extended Hessenberg matrix.
                E(:kp, :kp) = expm(tau*H(:kp, :kp))

                ! Project back into original space.
                call linear_combination(Xwrk, X(:kp), E(:kp, 1))
                call c%axpby(zero_rsp, Xwrk, one_rsp*beta)

                ! Cheap error esimate (this actually is the magnitude of the included correction
                ! and thus is very conservative).
                if (info == k) then
                    ! Approximation is exact.
                    err_est = 0.0_sp
                else
                    err_est = abs(E(kp, 1) * beta)
                endif

                ! Check convergence.
                if (err_est <= tol) exit expm_arnoldi

            enddo expm_arnoldi
        endif

        if (err_est <= tol) then
            info = kp
            write(msg,'(A,I0,2(A,E9.2))') 'Converged. kp= ', kp, ', err_est= ', err_est, ', tol= ', tol
            call logger%log_information(msg, module=this_module, procedure='kexpm_vec_rsp')
        else
            write(msg,'(A,I0,2(A,E9.2))') 'Not converged. kp= ', nk+1, ', err_est= ', err_est, ', tol= ', tol
            call logger%log_information(msg, module=this_module, procedure='kexpm_vec_rsp')
            info = -1
        endif

        return
    end subroutine kexpm_vec_rsp

    subroutine kexpm_mat_rsp(C, A, B, tau, tol, info, trans, kdim)
        class(abstract_vector_rsp), intent(out) :: C(:)
        !! Best Krylov approximation of \( \mathbf{C} = \exp(\tau \mathbf{A}) \mathbf{B} \).
        class(abstract_linop_rsp), intent(in) :: A
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

        ! ----- Internal variables -----
        integer, parameter :: kmax = 100
        integer :: i, j, k, p, kpm, kp, kpp, nk
        ! Block-Arnoldi factorization.
        class(abstract_vector_rsp), allocatable :: X(:)
        real(sp), allocatable :: H(:, :)
        ! Normalization & temporary arrays.
        real(sp), allocatable :: R(:, :), E(:, :), em(:, :)
        integer, allocatable :: perm(:), ptrans(:)
        class(abstract_vector_rsp), allocatable :: Xwrk(:), Cwrk(:)
        real(sp) :: err_est
        ! Optional arguments.
        logical :: transpose
        integer :: nsteps
        character(len=256) :: msg

        ! Determine block size.
        p = size(B)

        ! Deals with the optional args.
        transpose = optval(trans, .false.)
        nsteps    = optval(kdim, kmax)
        nk        = nsteps*p
        
        info = 0

        ! Allocate arrays.
        allocate(R(p, p)) ; allocate(perm(p)) ; allocate(ptrans(p))
        allocate(X(p*(nk+1)), source=B(1)) ; allocate(H(p*(nk+1), p*(nk+1)))
        allocate(E(p*(nk+1), p*(nk+1))) ; allocate(em(p, p))

        ! Scratch arrays.
        allocate(Xwrk(p), source=B) ; allocate(Cwrk(p), source=B(1))

        ! Normalize input matrix and initialize Krylov subspace.
        R = 0.0_sp
        call qr(Xwrk, R, perm, info) ; call apply_inverse_permutation_matrix(R, perm)

        if (norm2(abs(R)) == 0.0_sp) then
            ! Input matrix is zero.
            call zero_basis(C)
            err_est = 0.0_sp ; k = 0 ; kpp = p
        else
            call initialize_krylov_subspace(X, Xwrk) ; H = 0.0_sp

            expm_arnoldi: do k = 1, nk
                ! Set counters.
                kpm = (k-1)*p ; kp = kpm + p ; kpp = kp + p

                ! Reset working arrays.
                E = 0.0_sp ; call zero_basis(Xwrk)

                ! Compute the k-th step of the Arnoldi factorization.
                call arnoldi(A, X, H, info, kstart=k, kend=k, transpose=transpose, blksize=p)
                call check_info(info, 'arnoldi', module=this_module, procedure='kexpm_mat_rsp')

                if (info == kp) then
                    ! Arnoldi breakdown. Do not consider extended matrix.
                    kpp = kp
                    info = -2
                endif

                ! Compute the (dense) matrix exponential of the extended Hessenberg matrix.
                E(:kpp, :kpp) = expm(tau*H(:kpp, :kpp))

                ! Project back to original space.
                do i = 1, size(Xwrk)
                    call Xwrk(i)%zero()
                    do j = 1, kpp
                        call Xwrk(i)%axpby(one_rsp, X(j), E(j, i))
                    enddo
                enddo

                do i = 1, p
                    call C(i)%zero()
                    do j = 1, p
                        call C(i)%axpby(one_rsp, Xwrk(j), R(j, i))
                    enddo
                enddo

                ! Cheap error estimate.
                if (info == kp) then
                    ! Approximation is exact.
                    err_est = 0.0_sp
                else
                    em = matmul(E(kp+1:kpp, :p), R(:p, :p))
                    err_est = norm2(abs(em))
                endif

                if (err_est <= tol) exit expm_arnoldi

            enddo expm_arnoldi
        endif

        if (err_est .le. tol) then
            info = kpp
            write(msg,'(A,I0,2(A,E9.2))') 'Converged. kp= ', kpp, ', err_est= ', err_est, ', tol= ', tol
            call logger%log_information(msg, module=this_module, procedure='kexpm_mat_rsp')
        else
            write(msg,'(A,I0,2(A,E9.2))') 'Not converged. kp= ', kpp, ', err_est= ', err_est, ', tol= ', tol
            call logger%log_information(msg, module=this_module, procedure='kexpm_mat_rsp')
            info = -1
        endif

        return
    end subroutine kexpm_mat_rsp

    subroutine k_exptA_rsp(vec_out, A, vec_in, tau, info, trans)
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
        real(sp) :: tol
        integer :: kdim

        tol = atol_sp
        kdim = 30

        call kexpm(vec_out, A, vec_in, tau, tol, info, trans=trans, kdim=kdim)
        call check_info(info, 'kexpm', module=this_module, procedure='k_exptA_rsp')

        return
    end subroutine k_exptA_rsp

    function expm_rdp(A, order) result(E)
        real(dp), intent(in) :: A(:, :)
        !! Matrix to be exponentiated.
        real(dp) :: E(size(A, 1), size(A, 1))
        !! Output matrix E = exp(tA).
        integer, intent(in), optional :: order
        !! Order of the Pade approximation.

        !----- Internal variables -----
        real(dp), allocatable :: A2(:, :), Q(:, :), X(:, :)
        real(dp) :: a_norm, c
        integer :: n, ee, k, s
        logical :: p
        integer :: p_order

        ! Deal with optional args.
        p_order = optval(order, 10)

        n = size(A, 1)

        ! Allocate arrays.
        allocate(A2(n, n)) ; allocate(X(n, n)) ; allocate(Q(n, n))

        ! Compute the L-infinity norm.
        a_norm = norml(A)

        ! Determine scaling factor for the matrix.
        ee = int(log2(a_norm)) + 1
        s = max(0, ee+1)

        ! Scale the input matrix & initialize polynomial.
        A2 = A / 2.0_dp**s
        X = A2

        ! Initialize P & Q and add first step.
        c = 0.5_dp
        E = eye(n) ; E = E + c*A2

        Q = eye(n) ; Q = Q - c*A2

        ! Iteratively compute the Pade approximation.
        p = .true.
        do k = 2, p_order
            c = c*(p_order - k + 1) / (k * (2*p_order - k + 1))
            X = matmul(A2, X)
            E = E + c*X
            if (p) then
                Q = Q + c*X
            else
                Q = Q - c*X
            endif
            p = .not. p
        enddo

        E = matmul(inv(Q), E)
        do k = 1, s
            E = matmul(E, E)
        enddo

        return
    end function expm_rdp

    subroutine kexpm_vec_rdp(c, A, b, tau, tol, info, trans, kdim)
        class(abstract_vector_rdp), intent(out) :: c
        !! Best approximation of \( \exp(\tau \mathbf{A}) \mathbf{b} \) in the computed Krylov subspace
        class(abstract_linop_rdp), intent(in) :: A
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

        ! ----- Internal variables -----
        integer, parameter :: kmax = 100
        integer :: k, km, kp, nk
        ! Arnoldi factorization.
        class(abstract_vector_rdp), allocatable :: X(:)
        real(dp), allocatable :: H(:, :)
        ! Normaliztion & temporary arrays.
        real(dp), allocatable :: E(:, :)
        class(abstract_vector_rdp), allocatable :: Xwrk
        real(dp) :: err_est, beta
        ! Optional arguments.
        logical :: transpose
        integer :: nsteps
        character(len=256) :: msg
    
        ! Deals with optional args.
        transpose = optval(trans, .false.)
        nsteps    = optval(kdim, kmax)
        nk        = nsteps

        info = 0

        ! Allocate arrays.
        allocate(X(nk+1), source=b) ; allocate(H(nk+1, nk+1))
        allocate(E(nk+1, nk+1)) ; allocate(Xwrk, source=b)

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
            H = 0.0_dp

            expm_arnoldi: do k = 1, nk
                km = k-1 ; kp = k+1

                ! Reset work arrays.
                E = 0.0_dp

                ! Compute k-th step Arnoldi factorization.
                call arnoldi(A, X, H, info, kstart=k, kend=k, transpose=transpose)
                call check_info(info, 'arnoldi', module=this_module, procedure='kexpm_vec_rdp')

                ! Compute approximation.
                if (info == k) then
                    ! Arnoldi breakdown, do not consider extended matrix.
                    kp = k
                    info = -2
                endif

                ! Compute the (dense) matrix exponential of the extended Hessenberg matrix.
                E(:kp, :kp) = expm(tau*H(:kp, :kp))

                ! Project back into original space.
                call linear_combination(Xwrk, X(:kp), E(:kp, 1))
                call c%axpby(zero_rdp, Xwrk, one_rdp*beta)

                ! Cheap error esimate (this actually is the magnitude of the included correction
                ! and thus is very conservative).
                if (info == k) then
                    ! Approximation is exact.
                    err_est = 0.0_dp
                else
                    err_est = abs(E(kp, 1) * beta)
                endif

                ! Check convergence.
                if (err_est <= tol) exit expm_arnoldi

            enddo expm_arnoldi
        endif

        if (err_est <= tol) then
            info = kp
            write(msg,'(A,I0,2(A,E9.2))') 'Converged. kp= ', kp, ', err_est= ', err_est, ', tol= ', tol
            call logger%log_information(msg, module=this_module, procedure='kexpm_vec_rdp')
        else
            write(msg,'(A,I0,2(A,E9.2))') 'Not converged. kp= ', nk+1, ', err_est= ', err_est, ', tol= ', tol
            call logger%log_information(msg, module=this_module, procedure='kexpm_vec_rdp')
            info = -1
        endif

        return
    end subroutine kexpm_vec_rdp

    subroutine kexpm_mat_rdp(C, A, B, tau, tol, info, trans, kdim)
        class(abstract_vector_rdp), intent(out) :: C(:)
        !! Best Krylov approximation of \( \mathbf{C} = \exp(\tau \mathbf{A}) \mathbf{B} \).
        class(abstract_linop_rdp), intent(in) :: A
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

        ! ----- Internal variables -----
        integer, parameter :: kmax = 100
        integer :: i, j, k, p, kpm, kp, kpp, nk
        ! Block-Arnoldi factorization.
        class(abstract_vector_rdp), allocatable :: X(:)
        real(dp), allocatable :: H(:, :)
        ! Normalization & temporary arrays.
        real(dp), allocatable :: R(:, :), E(:, :), em(:, :)
        integer, allocatable :: perm(:), ptrans(:)
        class(abstract_vector_rdp), allocatable :: Xwrk(:), Cwrk(:)
        real(dp) :: err_est
        ! Optional arguments.
        logical :: transpose
        integer :: nsteps
        character(len=256) :: msg

        ! Determine block size.
        p = size(B)

        ! Deals with the optional args.
        transpose = optval(trans, .false.)
        nsteps    = optval(kdim, kmax)
        nk        = nsteps*p
        
        info = 0

        ! Allocate arrays.
        allocate(R(p, p)) ; allocate(perm(p)) ; allocate(ptrans(p))
        allocate(X(p*(nk+1)), source=B(1)) ; allocate(H(p*(nk+1), p*(nk+1)))
        allocate(E(p*(nk+1), p*(nk+1))) ; allocate(em(p, p))

        ! Scratch arrays.
        allocate(Xwrk(p), source=B) ; allocate(Cwrk(p), source=B(1))

        ! Normalize input matrix and initialize Krylov subspace.
        R = 0.0_dp
        call qr(Xwrk, R, perm, info) ; call apply_inverse_permutation_matrix(R, perm)

        if (norm2(abs(R)) == 0.0_dp) then
            ! Input matrix is zero.
            call zero_basis(C)
            err_est = 0.0_dp ; k = 0 ; kpp = p
        else
            call initialize_krylov_subspace(X, Xwrk) ; H = 0.0_dp

            expm_arnoldi: do k = 1, nk
                ! Set counters.
                kpm = (k-1)*p ; kp = kpm + p ; kpp = kp + p

                ! Reset working arrays.
                E = 0.0_dp ; call zero_basis(Xwrk)

                ! Compute the k-th step of the Arnoldi factorization.
                call arnoldi(A, X, H, info, kstart=k, kend=k, transpose=transpose, blksize=p)
                call check_info(info, 'arnoldi', module=this_module, procedure='kexpm_mat_rdp')

                if (info == kp) then
                    ! Arnoldi breakdown. Do not consider extended matrix.
                    kpp = kp
                    info = -2
                endif

                ! Compute the (dense) matrix exponential of the extended Hessenberg matrix.
                E(:kpp, :kpp) = expm(tau*H(:kpp, :kpp))

                ! Project back to original space.
                do i = 1, size(Xwrk)
                    call Xwrk(i)%zero()
                    do j = 1, kpp
                        call Xwrk(i)%axpby(one_rdp, X(j), E(j, i))
                    enddo
                enddo

                do i = 1, p
                    call C(i)%zero()
                    do j = 1, p
                        call C(i)%axpby(one_rdp, Xwrk(j), R(j, i))
                    enddo
                enddo

                ! Cheap error estimate.
                if (info == kp) then
                    ! Approximation is exact.
                    err_est = 0.0_dp
                else
                    em = matmul(E(kp+1:kpp, :p), R(:p, :p))
                    err_est = norm2(abs(em))
                endif

                if (err_est <= tol) exit expm_arnoldi

            enddo expm_arnoldi
        endif

        if (err_est .le. tol) then
            info = kpp
            write(msg,'(A,I0,2(A,E9.2))') 'Converged. kp= ', kpp, ', err_est= ', err_est, ', tol= ', tol
            call logger%log_information(msg, module=this_module, procedure='kexpm_mat_rdp')
        else
            write(msg,'(A,I0,2(A,E9.2))') 'Not converged. kp= ', kpp, ', err_est= ', err_est, ', tol= ', tol
            call logger%log_information(msg, module=this_module, procedure='kexpm_mat_rdp')
            info = -1
        endif

        return
    end subroutine kexpm_mat_rdp

    subroutine k_exptA_rdp(vec_out, A, vec_in, tau, info, trans)
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
        real(dp) :: tol
        integer :: kdim

        tol = atol_dp
        kdim = 30

        call kexpm(vec_out, A, vec_in, tau, tol, info, trans=trans, kdim=kdim)
        call check_info(info, 'kexpm', module=this_module, procedure='k_exptA_rdp')

        return
    end subroutine k_exptA_rdp

    function expm_csp(A, order) result(E)
        complex(sp), intent(in) :: A(:, :)
        !! Matrix to be exponentiated.
        complex(sp) :: E(size(A, 1), size(A, 1))
        !! Output matrix E = exp(tA).
        integer, intent(in), optional :: order
        !! Order of the Pade approximation.

        !----- Internal variables -----
        complex(sp), allocatable :: A2(:, :), Q(:, :), X(:, :)
        real(sp) :: a_norm, c
        integer :: n, ee, k, s
        logical :: p
        integer :: p_order

        ! Deal with optional args.
        p_order = optval(order, 10)

        n = size(A, 1)

        ! Allocate arrays.
        allocate(A2(n, n)) ; allocate(X(n, n)) ; allocate(Q(n, n))

        ! Compute the L-infinity norm.
        a_norm = norml(A)

        ! Determine scaling factor for the matrix.
        ee = int(log2(a_norm)) + 1
        s = max(0, ee+1)

        ! Scale the input matrix & initialize polynomial.
        A2 = A / 2.0_sp**s
        X = A2

        ! Initialize P & Q and add first step.
        c = 0.5_sp
        E = eye(n) ; E = E + c*A2

        Q = eye(n) ; Q = Q - c*A2

        ! Iteratively compute the Pade approximation.
        p = .true.
        do k = 2, p_order
            c = c*(p_order - k + 1) / (k * (2*p_order - k + 1))
            X = matmul(A2, X)
            E = E + c*X
            if (p) then
                Q = Q + c*X
            else
                Q = Q - c*X
            endif
            p = .not. p
        enddo

        E = matmul(inv(Q), E)
        do k = 1, s
            E = matmul(E, E)
        enddo

        return
    end function expm_csp

    subroutine kexpm_vec_csp(c, A, b, tau, tol, info, trans, kdim)
        class(abstract_vector_csp), intent(out) :: c
        !! Best approximation of \( \exp(\tau \mathbf{A}) \mathbf{b} \) in the computed Krylov subspace
        class(abstract_linop_csp), intent(in) :: A
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

        ! ----- Internal variables -----
        integer, parameter :: kmax = 100
        integer :: k, km, kp, nk
        ! Arnoldi factorization.
        class(abstract_vector_csp), allocatable :: X(:)
        complex(sp), allocatable :: H(:, :)
        ! Normaliztion & temporary arrays.
        complex(sp), allocatable :: E(:, :)
        class(abstract_vector_csp), allocatable :: Xwrk
        real(sp) :: err_est, beta
        ! Optional arguments.
        logical :: transpose
        integer :: nsteps
        character(len=256) :: msg
    
        ! Deals with optional args.
        transpose = optval(trans, .false.)
        nsteps    = optval(kdim, kmax)
        nk        = nsteps

        info = 0

        ! Allocate arrays.
        allocate(X(nk+1), source=b) ; allocate(H(nk+1, nk+1))
        allocate(E(nk+1, nk+1)) ; allocate(Xwrk, source=b)

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
            H = 0.0_sp

            expm_arnoldi: do k = 1, nk
                km = k-1 ; kp = k+1

                ! Reset work arrays.
                E = 0.0_sp

                ! Compute k-th step Arnoldi factorization.
                call arnoldi(A, X, H, info, kstart=k, kend=k, transpose=transpose)
                call check_info(info, 'arnoldi', module=this_module, procedure='kexpm_vec_csp')

                ! Compute approximation.
                if (info == k) then
                    ! Arnoldi breakdown, do not consider extended matrix.
                    kp = k
                    info = -2
                endif

                ! Compute the (dense) matrix exponential of the extended Hessenberg matrix.
                E(:kp, :kp) = expm(tau*H(:kp, :kp))

                ! Project back into original space.
                call linear_combination(Xwrk, X(:kp), E(:kp, 1))
                call c%axpby(zero_csp, Xwrk, one_csp*beta)

                ! Cheap error esimate (this actually is the magnitude of the included correction
                ! and thus is very conservative).
                if (info == k) then
                    ! Approximation is exact.
                    err_est = 0.0_sp
                else
                    err_est = abs(E(kp, 1) * beta)
                endif

                ! Check convergence.
                if (err_est <= tol) exit expm_arnoldi

            enddo expm_arnoldi
        endif

        if (err_est <= tol) then
            info = kp
            write(msg,'(A,I0,2(A,E9.2))') 'Converged. kp= ', kp, ', err_est= ', err_est, ', tol= ', tol
            call logger%log_information(msg, module=this_module, procedure='kexpm_vec_csp')
        else
            write(msg,'(A,I0,2(A,E9.2))') 'Not converged. kp= ', nk+1, ', err_est= ', err_est, ', tol= ', tol
            call logger%log_information(msg, module=this_module, procedure='kexpm_vec_csp')
            info = -1
        endif

        return
    end subroutine kexpm_vec_csp

    subroutine kexpm_mat_csp(C, A, B, tau, tol, info, trans, kdim)
        class(abstract_vector_csp), intent(out) :: C(:)
        !! Best Krylov approximation of \( \mathbf{C} = \exp(\tau \mathbf{A}) \mathbf{B} \).
        class(abstract_linop_csp), intent(in) :: A
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

        ! ----- Internal variables -----
        integer, parameter :: kmax = 100
        integer :: i, j, k, p, kpm, kp, kpp, nk
        ! Block-Arnoldi factorization.
        class(abstract_vector_csp), allocatable :: X(:)
        complex(sp), allocatable :: H(:, :)
        ! Normalization & temporary arrays.
        complex(sp), allocatable :: R(:, :), E(:, :), em(:, :)
        integer, allocatable :: perm(:), ptrans(:)
        class(abstract_vector_csp), allocatable :: Xwrk(:), Cwrk(:)
        real(sp) :: err_est
        ! Optional arguments.
        logical :: transpose
        integer :: nsteps
        character(len=256) :: msg

        ! Determine block size.
        p = size(B)

        ! Deals with the optional args.
        transpose = optval(trans, .false.)
        nsteps    = optval(kdim, kmax)
        nk        = nsteps*p
        
        info = 0

        ! Allocate arrays.
        allocate(R(p, p)) ; allocate(perm(p)) ; allocate(ptrans(p))
        allocate(X(p*(nk+1)), source=B(1)) ; allocate(H(p*(nk+1), p*(nk+1)))
        allocate(E(p*(nk+1), p*(nk+1))) ; allocate(em(p, p))

        ! Scratch arrays.
        allocate(Xwrk(p), source=B) ; allocate(Cwrk(p), source=B(1))

        ! Normalize input matrix and initialize Krylov subspace.
        R = 0.0_sp
        call qr(Xwrk, R, perm, info) ; call apply_inverse_permutation_matrix(R, perm)

        if (norm2(abs(R)) == 0.0_sp) then
            ! Input matrix is zero.
            call zero_basis(C)
            err_est = 0.0_sp ; k = 0 ; kpp = p
        else
            call initialize_krylov_subspace(X, Xwrk) ; H = 0.0_sp

            expm_arnoldi: do k = 1, nk
                ! Set counters.
                kpm = (k-1)*p ; kp = kpm + p ; kpp = kp + p

                ! Reset working arrays.
                E = 0.0_sp ; call zero_basis(Xwrk)

                ! Compute the k-th step of the Arnoldi factorization.
                call arnoldi(A, X, H, info, kstart=k, kend=k, transpose=transpose, blksize=p)
                call check_info(info, 'arnoldi', module=this_module, procedure='kexpm_mat_csp')

                if (info == kp) then
                    ! Arnoldi breakdown. Do not consider extended matrix.
                    kpp = kp
                    info = -2
                endif

                ! Compute the (dense) matrix exponential of the extended Hessenberg matrix.
                E(:kpp, :kpp) = expm(tau*H(:kpp, :kpp))

                ! Project back to original space.
                do i = 1, size(Xwrk)
                    call Xwrk(i)%zero()
                    do j = 1, kpp
                        call Xwrk(i)%axpby(one_csp, X(j), E(j, i))
                    enddo
                enddo

                do i = 1, p
                    call C(i)%zero()
                    do j = 1, p
                        call C(i)%axpby(one_csp, Xwrk(j), R(j, i))
                    enddo
                enddo

                ! Cheap error estimate.
                if (info == kp) then
                    ! Approximation is exact.
                    err_est = 0.0_sp
                else
                    em = matmul(E(kp+1:kpp, :p), R(:p, :p))
                    err_est = norm2(abs(em))
                endif

                if (err_est <= tol) exit expm_arnoldi

            enddo expm_arnoldi
        endif

        if (err_est .le. tol) then
            info = kpp
            write(msg,'(A,I0,2(A,E9.2))') 'Converged. kp= ', kpp, ', err_est= ', err_est, ', tol= ', tol
            call logger%log_information(msg, module=this_module, procedure='kexpm_mat_csp')
        else
            write(msg,'(A,I0,2(A,E9.2))') 'Not converged. kp= ', kpp, ', err_est= ', err_est, ', tol= ', tol
            call logger%log_information(msg, module=this_module, procedure='kexpm_mat_csp')
            info = -1
        endif

        return
    end subroutine kexpm_mat_csp

    subroutine k_exptA_csp(vec_out, A, vec_in, tau, info, trans)
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
        real(sp) :: tol
        integer :: kdim

        tol = atol_sp
        kdim = 30

        call kexpm(vec_out, A, vec_in, tau, tol, info, trans=trans, kdim=kdim)
        call check_info(info, 'kexpm', module=this_module, procedure='k_exptA_csp')

        return
    end subroutine k_exptA_csp

    function expm_cdp(A, order) result(E)
        complex(dp), intent(in) :: A(:, :)
        !! Matrix to be exponentiated.
        complex(dp) :: E(size(A, 1), size(A, 1))
        !! Output matrix E = exp(tA).
        integer, intent(in), optional :: order
        !! Order of the Pade approximation.

        !----- Internal variables -----
        complex(dp), allocatable :: A2(:, :), Q(:, :), X(:, :)
        real(dp) :: a_norm, c
        integer :: n, ee, k, s
        logical :: p
        integer :: p_order

        ! Deal with optional args.
        p_order = optval(order, 10)

        n = size(A, 1)

        ! Allocate arrays.
        allocate(A2(n, n)) ; allocate(X(n, n)) ; allocate(Q(n, n))

        ! Compute the L-infinity norm.
        a_norm = norml(A)

        ! Determine scaling factor for the matrix.
        ee = int(log2(a_norm)) + 1
        s = max(0, ee+1)

        ! Scale the input matrix & initialize polynomial.
        A2 = A / 2.0_dp**s
        X = A2

        ! Initialize P & Q and add first step.
        c = 0.5_dp
        E = eye(n) ; E = E + c*A2

        Q = eye(n) ; Q = Q - c*A2

        ! Iteratively compute the Pade approximation.
        p = .true.
        do k = 2, p_order
            c = c*(p_order - k + 1) / (k * (2*p_order - k + 1))
            X = matmul(A2, X)
            E = E + c*X
            if (p) then
                Q = Q + c*X
            else
                Q = Q - c*X
            endif
            p = .not. p
        enddo

        E = matmul(inv(Q), E)
        do k = 1, s
            E = matmul(E, E)
        enddo

        return
    end function expm_cdp

    subroutine kexpm_vec_cdp(c, A, b, tau, tol, info, trans, kdim)
        class(abstract_vector_cdp), intent(out) :: c
        !! Best approximation of \( \exp(\tau \mathbf{A}) \mathbf{b} \) in the computed Krylov subspace
        class(abstract_linop_cdp), intent(in) :: A
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

        ! ----- Internal variables -----
        integer, parameter :: kmax = 100
        integer :: k, km, kp, nk
        ! Arnoldi factorization.
        class(abstract_vector_cdp), allocatable :: X(:)
        complex(dp), allocatable :: H(:, :)
        ! Normaliztion & temporary arrays.
        complex(dp), allocatable :: E(:, :)
        class(abstract_vector_cdp), allocatable :: Xwrk
        real(dp) :: err_est, beta
        ! Optional arguments.
        logical :: transpose
        integer :: nsteps
        character(len=256) :: msg
    
        ! Deals with optional args.
        transpose = optval(trans, .false.)
        nsteps    = optval(kdim, kmax)
        nk        = nsteps

        info = 0

        ! Allocate arrays.
        allocate(X(nk+1), source=b) ; allocate(H(nk+1, nk+1))
        allocate(E(nk+1, nk+1)) ; allocate(Xwrk, source=b)

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
            H = 0.0_dp

            expm_arnoldi: do k = 1, nk
                km = k-1 ; kp = k+1

                ! Reset work arrays.
                E = 0.0_dp

                ! Compute k-th step Arnoldi factorization.
                call arnoldi(A, X, H, info, kstart=k, kend=k, transpose=transpose)
                call check_info(info, 'arnoldi', module=this_module, procedure='kexpm_vec_cdp')

                ! Compute approximation.
                if (info == k) then
                    ! Arnoldi breakdown, do not consider extended matrix.
                    kp = k
                    info = -2
                endif

                ! Compute the (dense) matrix exponential of the extended Hessenberg matrix.
                E(:kp, :kp) = expm(tau*H(:kp, :kp))

                ! Project back into original space.
                call linear_combination(Xwrk, X(:kp), E(:kp, 1))
                call c%axpby(zero_cdp, Xwrk, one_cdp*beta)

                ! Cheap error esimate (this actually is the magnitude of the included correction
                ! and thus is very conservative).
                if (info == k) then
                    ! Approximation is exact.
                    err_est = 0.0_dp
                else
                    err_est = abs(E(kp, 1) * beta)
                endif

                ! Check convergence.
                if (err_est <= tol) exit expm_arnoldi

            enddo expm_arnoldi
        endif

        if (err_est <= tol) then
            info = kp
            write(msg,'(A,I0,2(A,E9.2))') 'Converged. kp= ', kp, ', err_est= ', err_est, ', tol= ', tol
            call logger%log_information(msg, module=this_module, procedure='kexpm_vec_cdp')
        else
            write(msg,'(A,I0,2(A,E9.2))') 'Not converged. kp= ', nk+1, ', err_est= ', err_est, ', tol= ', tol
            call logger%log_information(msg, module=this_module, procedure='kexpm_vec_cdp')
            info = -1
        endif

        return
    end subroutine kexpm_vec_cdp

    subroutine kexpm_mat_cdp(C, A, B, tau, tol, info, trans, kdim)
        class(abstract_vector_cdp), intent(out) :: C(:)
        !! Best Krylov approximation of \( \mathbf{C} = \exp(\tau \mathbf{A}) \mathbf{B} \).
        class(abstract_linop_cdp), intent(in) :: A
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

        ! ----- Internal variables -----
        integer, parameter :: kmax = 100
        integer :: i, j, k, p, kpm, kp, kpp, nk
        ! Block-Arnoldi factorization.
        class(abstract_vector_cdp), allocatable :: X(:)
        complex(dp), allocatable :: H(:, :)
        ! Normalization & temporary arrays.
        complex(dp), allocatable :: R(:, :), E(:, :), em(:, :)
        integer, allocatable :: perm(:), ptrans(:)
        class(abstract_vector_cdp), allocatable :: Xwrk(:), Cwrk(:)
        real(dp) :: err_est
        ! Optional arguments.
        logical :: transpose
        integer :: nsteps
        character(len=256) :: msg

        ! Determine block size.
        p = size(B)

        ! Deals with the optional args.
        transpose = optval(trans, .false.)
        nsteps    = optval(kdim, kmax)
        nk        = nsteps*p
        
        info = 0

        ! Allocate arrays.
        allocate(R(p, p)) ; allocate(perm(p)) ; allocate(ptrans(p))
        allocate(X(p*(nk+1)), source=B(1)) ; allocate(H(p*(nk+1), p*(nk+1)))
        allocate(E(p*(nk+1), p*(nk+1))) ; allocate(em(p, p))

        ! Scratch arrays.
        allocate(Xwrk(p), source=B) ; allocate(Cwrk(p), source=B(1))

        ! Normalize input matrix and initialize Krylov subspace.
        R = 0.0_dp
        call qr(Xwrk, R, perm, info) ; call apply_inverse_permutation_matrix(R, perm)

        if (norm2(abs(R)) == 0.0_dp) then
            ! Input matrix is zero.
            call zero_basis(C)
            err_est = 0.0_dp ; k = 0 ; kpp = p
        else
            call initialize_krylov_subspace(X, Xwrk) ; H = 0.0_dp

            expm_arnoldi: do k = 1, nk
                ! Set counters.
                kpm = (k-1)*p ; kp = kpm + p ; kpp = kp + p

                ! Reset working arrays.
                E = 0.0_dp ; call zero_basis(Xwrk)

                ! Compute the k-th step of the Arnoldi factorization.
                call arnoldi(A, X, H, info, kstart=k, kend=k, transpose=transpose, blksize=p)
                call check_info(info, 'arnoldi', module=this_module, procedure='kexpm_mat_cdp')

                if (info == kp) then
                    ! Arnoldi breakdown. Do not consider extended matrix.
                    kpp = kp
                    info = -2
                endif

                ! Compute the (dense) matrix exponential of the extended Hessenberg matrix.
                E(:kpp, :kpp) = expm(tau*H(:kpp, :kpp))

                ! Project back to original space.
                do i = 1, size(Xwrk)
                    call Xwrk(i)%zero()
                    do j = 1, kpp
                        call Xwrk(i)%axpby(one_cdp, X(j), E(j, i))
                    enddo
                enddo

                do i = 1, p
                    call C(i)%zero()
                    do j = 1, p
                        call C(i)%axpby(one_cdp, Xwrk(j), R(j, i))
                    enddo
                enddo

                ! Cheap error estimate.
                if (info == kp) then
                    ! Approximation is exact.
                    err_est = 0.0_dp
                else
                    em = matmul(E(kp+1:kpp, :p), R(:p, :p))
                    err_est = norm2(abs(em))
                endif

                if (err_est <= tol) exit expm_arnoldi

            enddo expm_arnoldi
        endif

        if (err_est .le. tol) then
            info = kpp
            write(msg,'(A,I0,2(A,E9.2))') 'Converged. kp= ', kpp, ', err_est= ', err_est, ', tol= ', tol
            call logger%log_information(msg, module=this_module, procedure='kexpm_mat_cdp')
        else
            write(msg,'(A,I0,2(A,E9.2))') 'Not converged. kp= ', kpp, ', err_est= ', err_est, ', tol= ', tol
            call logger%log_information(msg, module=this_module, procedure='kexpm_mat_cdp')
            info = -1
        endif

        return
    end subroutine kexpm_mat_cdp

    subroutine k_exptA_cdp(vec_out, A, vec_in, tau, info, trans)
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
        real(dp) :: tol
        integer :: kdim

        tol = atol_dp
        kdim = 30

        call kexpm(vec_out, A, vec_in, tau, tol, info, trans=trans, kdim=kdim)
        call check_info(info, 'kexpm', module=this_module, procedure='k_exptA_cdp')

        return
    end subroutine k_exptA_cdp


end module LightKrylov_expmlib
