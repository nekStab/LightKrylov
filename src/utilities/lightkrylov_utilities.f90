module lightkrylov_utils_bis
    !--------------------------------------------
    !-----     Standard Fortran Library     -----
    !--------------------------------------------
    use iso_fortran_env, only: output_unit

    !-------------------------------
    !-----     LightKrylov     -----
    !-------------------------------
    use LightKrylov_Constants
    use LightKrylov_Logger

    implicit none
    private

    character(len=*), parameter :: this_module      = 'LK_Utils'
    character(len=*), parameter :: this_module_long = 'LightKrylov_Utils'

    !----------------------------------
    !-----     Public exports     -----
    !----------------------------------

    public :: assert_shape
    public :: log2
    public :: eig
    public :: ordschur
    public :: sqrtm

    !-------------------------------------------------
    !-----     Options for iterative solvers     -----
    !-------------------------------------------------

    type, abstract, public :: abstract_opts
    end type

    type, abstract, public :: abstract_metadata
        private
        contains
        procedure(abstract_print_metadata), pass(self), deferred, public :: print
        procedure(abstract_reset_metadata), pass(self), deferred, public :: reset
    end type

    abstract interface
        subroutine abstract_print_metadata(self, reset_counters, verbose)
            import abstract_metadata
            class(abstract_metadata), intent(inout) :: self
            logical, optional, intent(in) :: reset_counters
            logical, optional, intent(in) :: verbose
        end subroutine

        subroutine abstract_reset_metadata(self)
            import abstract_metadata
            class(abstract_metadata), intent(inout) :: self
        end subroutine
    end interface

    !-------------------------------------
    !-----     Utility functions     -----
    !-------------------------------------

    ! NOTE : Most of these functions will gradually disappear as more stable
    !        versions make their ways into the Fortran stdlib library.

    interface assert_shape
        module subroutine assert_shape_vector_rsp(v, size, vecname, module, procedure)
            real(sp), intent(in) :: v(:)
            integer, intent(in) :: size(:)
            character(len=*), intent(in) :: vecname
            character(len=*), intent(in) :: module
            character(len=*), intent(in) :: procedure
        end subroutine

        module subroutine assert_shape_matrix_rsp(A, size, matname, module, procedure)
            real(sp), intent(in) :: A(:, :)
            integer, intent(in) :: size(:)
            character(len=*), intent(in) :: matname
            character(len=*), intent(in) :: module
            character(len=*), intent(in) :: procedure
        end subroutine
        module subroutine assert_shape_vector_rdp(v, size, vecname, module, procedure)
            real(dp), intent(in) :: v(:)
            integer, intent(in) :: size(:)
            character(len=*), intent(in) :: vecname
            character(len=*), intent(in) :: module
            character(len=*), intent(in) :: procedure
        end subroutine

        module subroutine assert_shape_matrix_rdp(A, size, matname, module, procedure)
            real(dp), intent(in) :: A(:, :)
            integer, intent(in) :: size(:)
            character(len=*), intent(in) :: matname
            character(len=*), intent(in) :: module
            character(len=*), intent(in) :: procedure
        end subroutine
        module subroutine assert_shape_vector_csp(v, size, vecname, module, procedure)
            complex(sp), intent(in) :: v(:)
            integer, intent(in) :: size(:)
            character(len=*), intent(in) :: vecname
            character(len=*), intent(in) :: module
            character(len=*), intent(in) :: procedure
        end subroutine

        module subroutine assert_shape_matrix_csp(A, size, matname, module, procedure)
            complex(sp), intent(in) :: A(:, :)
            integer, intent(in) :: size(:)
            character(len=*), intent(in) :: matname
            character(len=*), intent(in) :: module
            character(len=*), intent(in) :: procedure
        end subroutine
        module subroutine assert_shape_vector_cdp(v, size, vecname, module, procedure)
            complex(dp), intent(in) :: v(:)
            integer, intent(in) :: size(:)
            character(len=*), intent(in) :: vecname
            character(len=*), intent(in) :: module
            character(len=*), intent(in) :: procedure
        end subroutine

        module subroutine assert_shape_matrix_cdp(A, size, matname, module, procedure)
            complex(dp), intent(in) :: A(:, :)
            integer, intent(in) :: size(:)
            character(len=*), intent(in) :: matname
            character(len=*), intent(in) :: module
            character(len=*), intent(in) :: procedure
        end subroutine
    end interface

    interface log2
        elemental real(sp) module function log2_rsp(x) result(y)
            real(sp), intent(in) :: x
        end function
        elemental real(dp) module function log2_rdp(x) result(y)
            real(dp), intent(in) :: x
        end function
    end interface

    interface eig
        module subroutine eig_rsp(A, vecs, vals)
            real(sp), intent(in) :: A(:, :)
            real(sp), intent(out) :: vecs(:, :)
            complex(sp), intent(out) :: vals(:)
        end subroutine
        module subroutine eig_rdp(A, vecs, vals)
            real(dp), intent(in) :: A(:, :)
            real(dp), intent(out) :: vecs(:, :)
            complex(dp), intent(out) :: vals(:)
        end subroutine
        module subroutine eig_csp(A, vecs, vals)
            complex(sp), intent(in) :: A(:, :)
            complex(sp), intent(out) :: vecs(:, :)
            complex(sp), intent(out) :: vals(:)
        end subroutine
        module subroutine eig_cdp(A, vecs, vals)
            complex(dp), intent(in) :: A(:, :)
            complex(dp), intent(out) :: vecs(:, :)
            complex(dp), intent(out) :: vals(:)
        end subroutine
    end interface

    interface ordschur
        module subroutine ordschur_rsp(T, Q, selected)
            real(sp), intent(inout) :: T(:, :)
            real(sp), intent(inout) :: Q(:, :)
            logical, intent(in) :: selected(:)
        end subroutine
        module subroutine ordschur_rdp(T, Q, selected)
            real(dp), intent(inout) :: T(:, :)
            real(dp), intent(inout) :: Q(:, :)
            logical, intent(in) :: selected(:)
        end subroutine
        module subroutine ordschur_csp(T, Q, selected)
            complex(sp), intent(inout) :: T(:, :)
            complex(sp), intent(inout) :: Q(:, :)
            logical, intent(in) :: selected(:)
        end subroutine
        module subroutine ordschur_cdp(T, Q, selected)
            complex(dp), intent(inout) :: T(:, :)
            complex(dp), intent(inout) :: Q(:, :)
            logical, intent(in) :: selected(:)
        end subroutine
    end interface

    interface sqrtm
        module subroutine sqrtm_rsp(A, sqrtA, info)
            real(sp), intent(inout) :: A(:, :)
            real(sp), intent(out) :: sqrtA(:, :)
            integer, intent(out) :: info
        end subroutine
        module subroutine sqrtm_rdp(A, sqrtA, info)
            real(dp), intent(inout) :: A(:, :)
            real(dp), intent(out) :: sqrtA(:, :)
            integer, intent(out) :: info
        end subroutine
        module subroutine sqrtm_csp(A, sqrtA, info)
            complex(sp), intent(inout) :: A(:, :)
            complex(sp), intent(out) :: sqrtA(:, :)
            integer, intent(out) :: info
        end subroutine
        module subroutine sqrtm_cdp(A, sqrtA, info)
            complex(dp), intent(inout) :: A(:, :)
            complex(dp), intent(out) :: sqrtA(:, :)
            integer, intent(out) :: info
        end subroutine
    end interface
contains
end module
