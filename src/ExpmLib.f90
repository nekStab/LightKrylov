module lightkrylov_expmlib

    ! Iso Fortran.
    use iso_fortran_env, only: output_unit

    ! Fortran standard library.
    use stdlib_optval, only: optval
    use stdlib_linalg, only: eye

    ! LightKrylov.
    use lightkrylov_constants
    use lightkrylov_utils
    use lightkrylov_AbstractVectors
    use lightkrylov_AbstractLinops
    use lightkrylov_BaseKrylov

    implicit none
    private

    public :: abstract_exptA_rsp
    public :: abstract_exptA_rdp
    public :: abstract_exptA_csp
    public :: abstract_exptA_cdp
    public :: expm

    abstract interface
        subroutine abstract_exptA_rsp(vec_out, A, vec_in, tau, info, trans)
            import sp
            import abstract_vector_rsp
            import abstract_linop_rsp
            !! Abstract interface to define the matrix exponential-vector product.
            class(abstract_vector_rsp), intent(out) :: vec_out
            !! Solution vector.
            class(abstract_linop_rsp), intent(in) :: A
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
            class(abstract_linop_rdp), intent(in) :: A
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
            class(abstract_linop_csp), intent(in) :: A
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
            class(abstract_linop_cdp), intent(in) :: A
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
        module procedure expm_rsp
        module procedure expm_rdp
        module procedure expm_csp
        module procedure expm_cdp
    end interface

contains

    !--------------------------------------------
    !-----     DENSE MATRIX EXPONENTIAL     -----
    !--------------------------------------------

    subroutine expm_rsp(E, A, order)
        real(sp), intent(in) :: A(:, :)
        !! Matrix to be exponentiated.
        real(sp), intent(out) :: E(:, :)
        !! Output matrix E = exp(tA).
        integer, optional :: order
        !! Order of the Pade approximation.

        !----- Internal variables -----
        real(sp), allocatable :: A2(:, :), Q(:, :), X(:, :), invQ(:, :), wrk(:)
        real(sp) :: a_norm, c
        integer :: n, ee, k, s
        logical :: p
        integer :: p_order

        ! Deal with optional args.
        p_order = optval(order, 10)

        n = size(A, 1)

        ! Allocate arrays.
        allocate(A2(n, n)) ; allocate(X(n, n))
        allocate(Q(n, n)) ; allocate(invQ(n, n))
        allocate(wrk(n))

        ! Compute the L-infinity norm.
        a_norm = norml_rsp(A)

        ! Determine scaling factor for the matrix.
        ee = int(log2_rsp(a_norm)) + 1
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

        invQ = Q ; call inv(invQ)
        E = matmul(invQ, E)

        do k = 1, s
            E = matmul(E, E)
        enddo

        return
    end subroutine expm_rsp
    subroutine expm_rdp(E, A, order)
        real(dp), intent(in) :: A(:, :)
        !! Matrix to be exponentiated.
        real(dp), intent(out) :: E(:, :)
        !! Output matrix E = exp(tA).
        integer, optional :: order
        !! Order of the Pade approximation.

        !----- Internal variables -----
        real(dp), allocatable :: A2(:, :), Q(:, :), X(:, :), invQ(:, :), wrk(:)
        real(dp) :: a_norm, c
        integer :: n, ee, k, s
        logical :: p
        integer :: p_order

        ! Deal with optional args.
        p_order = optval(order, 10)

        n = size(A, 1)

        ! Allocate arrays.
        allocate(A2(n, n)) ; allocate(X(n, n))
        allocate(Q(n, n)) ; allocate(invQ(n, n))
        allocate(wrk(n))

        ! Compute the L-infinity norm.
        a_norm = norml_rdp(A)

        ! Determine scaling factor for the matrix.
        ee = int(log2_rdp(a_norm)) + 1
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

        invQ = Q ; call inv(invQ)
        E = matmul(invQ, E)

        do k = 1, s
            E = matmul(E, E)
        enddo

        return
    end subroutine expm_rdp
    subroutine expm_csp(E, A, order)
        complex(sp), intent(in) :: A(:, :)
        !! Matrix to be exponentiated.
        complex(sp), intent(out) :: E(:, :)
        !! Output matrix E = exp(tA).
        integer, optional :: order
        !! Order of the Pade approximation.

        !----- Internal variables -----
        complex(sp), allocatable :: A2(:, :), Q(:, :), X(:, :), invQ(:, :), wrk(:)
        real(sp) :: a_norm, c
        integer :: n, ee, k, s
        logical :: p
        integer :: p_order

        ! Deal with optional args.
        p_order = optval(order, 10)

        n = size(A, 1)

        ! Allocate arrays.
        allocate(A2(n, n)) ; allocate(X(n, n))
        allocate(Q(n, n)) ; allocate(invQ(n, n))
        allocate(wrk(n))

        ! Compute the L-infinity norm.
        a_norm = norml_csp(A)

        ! Determine scaling factor for the matrix.
        ee = int(log2_rsp(a_norm)) + 1
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

        invQ = Q ; call inv(invQ)
        E = matmul(invQ, E)

        do k = 1, s
            E = matmul(E, E)
        enddo

        return
    end subroutine expm_csp
    subroutine expm_cdp(E, A, order)
        complex(dp), intent(in) :: A(:, :)
        !! Matrix to be exponentiated.
        complex(dp), intent(out) :: E(:, :)
        !! Output matrix E = exp(tA).
        integer, optional :: order
        !! Order of the Pade approximation.

        !----- Internal variables -----
        complex(dp), allocatable :: A2(:, :), Q(:, :), X(:, :), invQ(:, :), wrk(:)
        real(dp) :: a_norm, c
        integer :: n, ee, k, s
        logical :: p
        integer :: p_order

        ! Deal with optional args.
        p_order = optval(order, 10)

        n = size(A, 1)

        ! Allocate arrays.
        allocate(A2(n, n)) ; allocate(X(n, n))
        allocate(Q(n, n)) ; allocate(invQ(n, n))
        allocate(wrk(n))

        ! Compute the L-infinity norm.
        a_norm = norml_cdp(A)

        ! Determine scaling factor for the matrix.
        ee = int(log2_rdp(a_norm)) + 1
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

        invQ = Q ; call inv(invQ)
        E = matmul(invQ, E)

        do k = 1, s
            E = matmul(E, E)
        enddo

        return
    end subroutine expm_cdp

end module lightkrylov_expmlib


