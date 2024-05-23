module lightkrylov_constants
    implicit none
    private

    integer , parameter, public :: sp = selected_real_kind(6, 37)
    !! Definition of the single precision data type.
    real(sp), parameter, public :: atol_sp = 10*epsilon(1.0_sp)
    !! Definition of the absolute tolerance for single precision computations.
    real(sp), parameter, public :: rtol_sp = sqrt(atol_sp)
    !! Definition of the relative tolerance for single precision computations.

    integer , parameter, public :: dp = selected_real_kind(15, 307)
    !! Definition of the double precision data type.
    real(dp), parameter, public :: atol_dp = epsilon(1.0_dp)
    !! Definition of the absolute tolerance for double precision computations.
    real(dp), parameter, public :: rtol_dp = sqrt(atol_dp)
    !! Definition of the relative tolerance for double precision computations.

end module lightkrylov_constants
