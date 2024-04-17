module lightkrylov_constants
    implicit none
    private

    integer , parameter, public :: sp = selected_real_kind(6, 37)
    real(sp), parameter, public :: atol_sp = epsilon(1.0_sp)
    real(sp), parameter, public :: rtol_sp = sqrt(atol_sp)

    integer , parameter, public :: dp = selected_real_kind(15, 307)
    real(dp), parameter, public :: atol_dp = epsilon(1.0_dp)
    real(dp), parameter, public :: rtol_dp = sqrt(atol_dp)

end module lightkrylov_constants
