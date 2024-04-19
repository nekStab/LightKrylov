module lightkrylov_AbstractVectors
    use lightkrylov_constants
    implicit none
    private

    public :: innerprod_matrix

    interface innerprod_matrix
        module procedure innerprod_matrix_rsp
        module procedure innerprod_matrix_rdp
        module procedure innerprod_matrix_csp
        module procedure innerprod_matrix_cdp
    end interface


    type, abstract, public :: abstract_vector
    contains
    private
        procedure, pass(from), public :: copy
        generic, public :: assignment(=) => copy
    end type abstract_vector




    !----------------------------------------------------------------------------
    !-----     Definition of an abstract real(sp) vector with kind=sp     -----
    !----------------------------------------------------------------------------

    type, abstract, extends(abstract_vector), public :: abstract_vector_rsp
    contains
        private
        procedure(abstract_zero_rsp), pass(self), deferred, public :: zero
        !! Sets and `abstract_vector_rsp` to zero.
        procedure(abstract_rand_rsp), pass(self), deferred, public :: rand
        !! Creates a random `abstract_vector_rsp.
        procedure(abstract_scal_rsp), pass(self), deferred, public :: scal
        !! Compute the scalar-vector product.
        procedure(abstract_axpby_rsp), pass(self), deferred, public :: axpby
        !! In-place computation of \( \mathbf{x} = \alpha \mathbf{x} + \beta \mathbf{y} \).
        procedure(abstract_dot_rsp), pass(self), deferred, public :: dot
        !! Computes the dot product between two `abstract_vector_rsp`.
        procedure, pass(self), public :: norm => norm_rsp
        !! Computes the norm of the `abstract_vector`.
        procedure, pass(self), public :: add => add_rsp
        !! Adds two `abstract_vector`.
        procedure, pass(self), public :: sub => sub_rsp
        !! Subtracts two `abstract_vector`.
        procedure, pass(self), public :: chsgn => chsgn_rsp
    end type

    abstract interface
        subroutine abstract_zero_rsp(self)
            !! Abstract interface to zero-out a vector in-place.
            import abstract_vector_rsp, sp

            class(abstract_vector_rsp), intent(inout) :: self
            !! Vector to be zeroed-out.
        end subroutine abstract_zero_rsp

        subroutine abstract_rand_rsp(self, ifnorm)
            !! Abstract interface to generate a random (normalized) vector.
            import abstract_vector_rsp, sp

            class(abstract_vector_rsp), intent(inout) :: self
            logical, optional, intent(in) :: ifnorm
        end subroutine abstract_rand_rsp

        subroutine abstract_scal_rsp(self, alpha)
            !! Abstract interface to scale a vector.
            import abstract_vector_rsp, sp

            class(abstract_vector_rsp), intent(inout) :: self
            !! Input/Output vector.
            real(sp), intent(in) :: alpha
        end subroutine abstract_scal_rsp

        subroutine abstract_axpby_rsp(self, alpha, vec, beta)
            !! Abstract interface to add/scale two vectors in-place.
            import abstract_vector_rsp, sp

            class(abstract_vector_rsp), intent(inout) :: self
            !! Input/Output vector.
            class(abstract_vector_rsp), intent(in) :: vec
            !! Vector to be added/subtracted.
            real(sp), intent(in) :: alpha, beta
        end subroutine abstract_axpby_rsp

        function abstract_dot_rsp(self, vec) result(alpha)
            !! Abstract interface to compute the dot product.
            import abstract_vector_rsp, sp

            class(abstract_vector_rsp), intent(in) :: self, vec
            !! Vectors whose dot product will be computed.
            real(sp) :: alpha
        end function abstract_dot_rsp
    end interface





    !----------------------------------------------------------------------------
    !-----     Definition of an abstract real(dp) vector with kind=dp     -----
    !----------------------------------------------------------------------------

    type, abstract, extends(abstract_vector), public :: abstract_vector_rdp
    contains
        private
        procedure(abstract_zero_rdp), pass(self), deferred, public :: zero
        !! Sets and `abstract_vector_rdp` to zero.
        procedure(abstract_rand_rdp), pass(self), deferred, public :: rand
        !! Creates a random `abstract_vector_rdp.
        procedure(abstract_scal_rdp), pass(self), deferred, public :: scal
        !! Compute the scalar-vector product.
        procedure(abstract_axpby_rdp), pass(self), deferred, public :: axpby
        !! In-place computation of \( \mathbf{x} = \alpha \mathbf{x} + \beta \mathbf{y} \).
        procedure(abstract_dot_rdp), pass(self), deferred, public :: dot
        !! Computes the dot product between two `abstract_vector_rdp`.
        procedure, pass(self), public :: norm => norm_rdp
        !! Computes the norm of the `abstract_vector`.
        procedure, pass(self), public :: add => add_rdp
        !! Adds two `abstract_vector`.
        procedure, pass(self), public :: sub => sub_rdp
        !! Subtracts two `abstract_vector`.
        procedure, pass(self), public :: chsgn => chsgn_rdp
    end type

    abstract interface
        subroutine abstract_zero_rdp(self)
            !! Abstract interface to zero-out a vector in-place.
            import abstract_vector_rdp, dp

            class(abstract_vector_rdp), intent(inout) :: self
            !! Vector to be zeroed-out.
        end subroutine abstract_zero_rdp

        subroutine abstract_rand_rdp(self, ifnorm)
            !! Abstract interface to generate a random (normalized) vector.
            import abstract_vector_rdp, dp

            class(abstract_vector_rdp), intent(inout) :: self
            logical, optional, intent(in) :: ifnorm
        end subroutine abstract_rand_rdp

        subroutine abstract_scal_rdp(self, alpha)
            !! Abstract interface to scale a vector.
            import abstract_vector_rdp, dp

            class(abstract_vector_rdp), intent(inout) :: self
            !! Input/Output vector.
            real(dp), intent(in) :: alpha
        end subroutine abstract_scal_rdp

        subroutine abstract_axpby_rdp(self, alpha, vec, beta)
            !! Abstract interface to add/scale two vectors in-place.
            import abstract_vector_rdp, dp

            class(abstract_vector_rdp), intent(inout) :: self
            !! Input/Output vector.
            class(abstract_vector_rdp), intent(in) :: vec
            !! Vector to be added/subtracted.
            real(dp), intent(in) :: alpha, beta
        end subroutine abstract_axpby_rdp

        function abstract_dot_rdp(self, vec) result(alpha)
            !! Abstract interface to compute the dot product.
            import abstract_vector_rdp, dp

            class(abstract_vector_rdp), intent(in) :: self, vec
            !! Vectors whose dot product will be computed.
            real(dp) :: alpha
        end function abstract_dot_rdp
    end interface





    !----------------------------------------------------------------------------
    !-----     Definition of an abstract complex(sp) vector with kind=sp     -----
    !----------------------------------------------------------------------------

    type, abstract, extends(abstract_vector), public :: abstract_vector_csp
    contains
        private
        procedure(abstract_zero_csp), pass(self), deferred, public :: zero
        !! Sets and `abstract_vector_csp` to zero.
        procedure(abstract_rand_csp), pass(self), deferred, public :: rand
        !! Creates a random `abstract_vector_csp.
        procedure(abstract_scal_csp), pass(self), deferred, public :: scal
        !! Compute the scalar-vector product.
        procedure(abstract_axpby_csp), pass(self), deferred, public :: axpby
        !! In-place computation of \( \mathbf{x} = \alpha \mathbf{x} + \beta \mathbf{y} \).
        procedure(abstract_dot_csp), pass(self), deferred, public :: dot
        !! Computes the dot product between two `abstract_vector_csp`.
        procedure, pass(self), public :: norm => norm_csp
        !! Computes the norm of the `abstract_vector`.
        procedure, pass(self), public :: add => add_csp
        !! Adds two `abstract_vector`.
        procedure, pass(self), public :: sub => sub_csp
        !! Subtracts two `abstract_vector`.
        procedure, pass(self), public :: chsgn => chsgn_csp
    end type

    abstract interface
        subroutine abstract_zero_csp(self)
            !! Abstract interface to zero-out a vector in-place.
            import abstract_vector_csp, sp

            class(abstract_vector_csp), intent(inout) :: self
            !! Vector to be zeroed-out.
        end subroutine abstract_zero_csp

        subroutine abstract_rand_csp(self, ifnorm)
            !! Abstract interface to generate a random (normalized) vector.
            import abstract_vector_csp, sp

            class(abstract_vector_csp), intent(inout) :: self
            logical, optional, intent(in) :: ifnorm
        end subroutine abstract_rand_csp

        subroutine abstract_scal_csp(self, alpha)
            !! Abstract interface to scale a vector.
            import abstract_vector_csp, sp

            class(abstract_vector_csp), intent(inout) :: self
            !! Input/Output vector.
            complex(sp), intent(in) :: alpha
        end subroutine abstract_scal_csp

        subroutine abstract_axpby_csp(self, alpha, vec, beta)
            !! Abstract interface to add/scale two vectors in-place.
            import abstract_vector_csp, sp

            class(abstract_vector_csp), intent(inout) :: self
            !! Input/Output vector.
            class(abstract_vector_csp), intent(in) :: vec
            !! Vector to be added/subtracted.
            complex(sp), intent(in) :: alpha, beta
        end subroutine abstract_axpby_csp

        function abstract_dot_csp(self, vec) result(alpha)
            !! Abstract interface to compute the dot product.
            import abstract_vector_csp, sp

            class(abstract_vector_csp), intent(in) :: self, vec
            !! Vectors whose dot product will be computed.
            complex(sp) :: alpha
        end function abstract_dot_csp
    end interface





    !----------------------------------------------------------------------------
    !-----     Definition of an abstract complex(dp) vector with kind=dp     -----
    !----------------------------------------------------------------------------

    type, abstract, extends(abstract_vector), public :: abstract_vector_cdp
    contains
        private
        procedure(abstract_zero_cdp), pass(self), deferred, public :: zero
        !! Sets and `abstract_vector_cdp` to zero.
        procedure(abstract_rand_cdp), pass(self), deferred, public :: rand
        !! Creates a random `abstract_vector_cdp.
        procedure(abstract_scal_cdp), pass(self), deferred, public :: scal
        !! Compute the scalar-vector product.
        procedure(abstract_axpby_cdp), pass(self), deferred, public :: axpby
        !! In-place computation of \( \mathbf{x} = \alpha \mathbf{x} + \beta \mathbf{y} \).
        procedure(abstract_dot_cdp), pass(self), deferred, public :: dot
        !! Computes the dot product between two `abstract_vector_cdp`.
        procedure, pass(self), public :: norm => norm_cdp
        !! Computes the norm of the `abstract_vector`.
        procedure, pass(self), public :: add => add_cdp
        !! Adds two `abstract_vector`.
        procedure, pass(self), public :: sub => sub_cdp
        !! Subtracts two `abstract_vector`.
        procedure, pass(self), public :: chsgn => chsgn_cdp
    end type

    abstract interface
        subroutine abstract_zero_cdp(self)
            !! Abstract interface to zero-out a vector in-place.
            import abstract_vector_cdp, dp

            class(abstract_vector_cdp), intent(inout) :: self
            !! Vector to be zeroed-out.
        end subroutine abstract_zero_cdp

        subroutine abstract_rand_cdp(self, ifnorm)
            !! Abstract interface to generate a random (normalized) vector.
            import abstract_vector_cdp, dp

            class(abstract_vector_cdp), intent(inout) :: self
            logical, optional, intent(in) :: ifnorm
        end subroutine abstract_rand_cdp

        subroutine abstract_scal_cdp(self, alpha)
            !! Abstract interface to scale a vector.
            import abstract_vector_cdp, dp

            class(abstract_vector_cdp), intent(inout) :: self
            !! Input/Output vector.
            complex(dp), intent(in) :: alpha
        end subroutine abstract_scal_cdp

        subroutine abstract_axpby_cdp(self, alpha, vec, beta)
            !! Abstract interface to add/scale two vectors in-place.
            import abstract_vector_cdp, dp

            class(abstract_vector_cdp), intent(inout) :: self
            !! Input/Output vector.
            class(abstract_vector_cdp), intent(in) :: vec
            !! Vector to be added/subtracted.
            complex(dp), intent(in) :: alpha, beta
        end subroutine abstract_axpby_cdp

        function abstract_dot_cdp(self, vec) result(alpha)
            !! Abstract interface to compute the dot product.
            import abstract_vector_cdp, dp

            class(abstract_vector_cdp), intent(in) :: self, vec
            !! Vectors whose dot product will be computed.
            complex(dp) :: alpha
        end function abstract_dot_cdp
    end interface





contains

    subroutine copy(out, from)
        class(abstract_vector), intent(in)  :: from
        class(abstract_vector), allocatable, intent(out) :: out
        if (allocated(out)) deallocate(out)
        allocate(out, source=from)
        return
    end subroutine copy

    function norm_rsp(self) result(alpha)
        class(abstract_vector_rsp), intent(in) :: self
        real(sp) :: alpha
        alpha = abs(self%dot(self)) ; alpha = sqrt(alpha)
    end function norm_rsp

    subroutine sub_rsp(self, vec)
        !! Subtract two `abstract_vector` in-place.
        class(abstract_vector_rsp), intent(inout) :: self
        !! Input/Output vector.
        class(abstract_vector_rsp), intent(in) :: vec
        !! Vector to be added.
        call self%axpby(1.0_sp, vec, -1.0_sp)
    end subroutine sub_rsp

    subroutine add_rsp(self, vec)
        !! Add two `abstract_vector` in-place.
        class(abstract_vector_rsp), intent(inout) :: self
        !! Input/Output vector.
        class(abstract_vector_rsp), intent(in) :: vec
        !! Vector to be added.
        call self%axpby(1.0_sp, vec, 1.0_sp)
    end subroutine add_rsp

    subroutine chsgn_rsp(self)
        !! Changes the sign of the `abstract_vector`.
        class(abstract_vector_rsp), intent(inout) :: self
        call self%scal(-1.0_sp)
    end subroutine chsgn_rsp

    function norm_rdp(self) result(alpha)
        class(abstract_vector_rdp), intent(in) :: self
        real(dp) :: alpha
        alpha = abs(self%dot(self)) ; alpha = sqrt(alpha)
    end function norm_rdp

    subroutine sub_rdp(self, vec)
        !! Subtract two `abstract_vector` in-place.
        class(abstract_vector_rdp), intent(inout) :: self
        !! Input/Output vector.
        class(abstract_vector_rdp), intent(in) :: vec
        !! Vector to be added.
        call self%axpby(1.0_dp, vec, -1.0_dp)
    end subroutine sub_rdp

    subroutine add_rdp(self, vec)
        !! Add two `abstract_vector` in-place.
        class(abstract_vector_rdp), intent(inout) :: self
        !! Input/Output vector.
        class(abstract_vector_rdp), intent(in) :: vec
        !! Vector to be added.
        call self%axpby(1.0_dp, vec, 1.0_dp)
    end subroutine add_rdp

    subroutine chsgn_rdp(self)
        !! Changes the sign of the `abstract_vector`.
        class(abstract_vector_rdp), intent(inout) :: self
        call self%scal(-1.0_dp)
    end subroutine chsgn_rdp

    function norm_csp(self) result(alpha)
        class(abstract_vector_csp), intent(in) :: self
        real(sp) :: alpha
        alpha = abs(self%dot(self)) ; alpha = sqrt(alpha)
    end function norm_csp

    subroutine sub_csp(self, vec)
        !! Subtract two `abstract_vector` in-place.
        class(abstract_vector_csp), intent(inout) :: self
        !! Input/Output vector.
        class(abstract_vector_csp), intent(in) :: vec
        !! Vector to be added.
        call self%axpby(cmplx(1.0_sp, 0.0_sp, kind=sp), vec, cmplx(-1.0_sp, 0.0_sp, kind=sp))
    end subroutine sub_csp

    subroutine add_csp(self, vec)
        !! Add two `abstract_vector` in-place.
        class(abstract_vector_csp), intent(inout) :: self
        !! Input/Output vector.
        class(abstract_vector_csp), intent(in) :: vec
        !! Vector to be added.
        call self%axpby(cmplx(1.0_sp, 0.0_sp, kind=sp), vec, cmplx(1.0_sp, 0.0_sp, kind=sp))
    end subroutine add_csp

    subroutine chsgn_csp(self)
        !! Changes the sign of the `abstract_vector`.
        class(abstract_vector_csp), intent(inout) :: self
        call self%scal(cmplx(-1.0_sp, 0.0_sp, kind=sp))
    end subroutine chsgn_csp

    function norm_cdp(self) result(alpha)
        class(abstract_vector_cdp), intent(in) :: self
        real(dp) :: alpha
        alpha = abs(self%dot(self)) ; alpha = sqrt(alpha)
    end function norm_cdp

    subroutine sub_cdp(self, vec)
        !! Subtract two `abstract_vector` in-place.
        class(abstract_vector_cdp), intent(inout) :: self
        !! Input/Output vector.
        class(abstract_vector_cdp), intent(in) :: vec
        !! Vector to be added.
        call self%axpby(cmplx(1.0_dp, 0.0_dp, kind=dp), vec, cmplx(-1.0_dp, 0.0_dp, kind=dp))
    end subroutine sub_cdp

    subroutine add_cdp(self, vec)
        !! Add two `abstract_vector` in-place.
        class(abstract_vector_cdp), intent(inout) :: self
        !! Input/Output vector.
        class(abstract_vector_cdp), intent(in) :: vec
        !! Vector to be added.
        call self%axpby(cmplx(1.0_dp, 0.0_dp, kind=dp), vec, cmplx(1.0_dp, 0.0_dp, kind=dp))
    end subroutine add_cdp

    subroutine chsgn_cdp(self)
        !! Changes the sign of the `abstract_vector`.
        class(abstract_vector_cdp), intent(inout) :: self
        call self%scal(cmplx(-1.0_dp, 0.0_dp, kind=dp))
    end subroutine chsgn_cdp

    
    !--------------------------------------
    !-----      UTILITY FUNCTIONS     -----
    !--------------------------------------

    subroutine innerprod_matrix_rsp(M, X, Y)
        class(abstract_vector_rsp), intent(in) :: X(:)
        class(abstract_vector_rsp), intent(in) :: Y(:)
        real(sp), intent(out) :: M(size(X), size(Y))

        ! Local variables.
        integer :: i, j

        M = 0.0_sp
        do j = 1, size(Y)
            do i = 1, size(X)
                M(i, j) = X(i)%dot(Y(j))
            enddo
        enddo

        return
    end subroutine innerprod_matrix_rsp
    
    subroutine innerprod_matrix_rdp(M, X, Y)
        class(abstract_vector_rdp), intent(in) :: X(:)
        class(abstract_vector_rdp), intent(in) :: Y(:)
        real(dp), intent(out) :: M(size(X), size(Y))

        ! Local variables.
        integer :: i, j

        M = 0.0_dp
        do j = 1, size(Y)
            do i = 1, size(X)
                M(i, j) = X(i)%dot(Y(j))
            enddo
        enddo

        return
    end subroutine innerprod_matrix_rdp
    
    subroutine innerprod_matrix_csp(M, X, Y)
        class(abstract_vector_csp), intent(in) :: X(:)
        class(abstract_vector_csp), intent(in) :: Y(:)
        complex(sp), intent(out) :: M(size(X), size(Y))

        ! Local variables.
        integer :: i, j

        M = 0.0_sp
        do j = 1, size(Y)
            do i = 1, size(X)
                M(i, j) = X(i)%dot(Y(j))
            enddo
        enddo

        return
    end subroutine innerprod_matrix_csp
    
    subroutine innerprod_matrix_cdp(M, X, Y)
        class(abstract_vector_cdp), intent(in) :: X(:)
        class(abstract_vector_cdp), intent(in) :: Y(:)
        complex(dp), intent(out) :: M(size(X), size(Y))

        ! Local variables.
        integer :: i, j

        M = 0.0_dp
        do j = 1, size(Y)
            do i = 1, size(X)
                M(i, j) = X(i)%dot(Y(j))
            enddo
        enddo

        return
    end subroutine innerprod_matrix_cdp
    

end module lightkrylov_AbstractVectors
