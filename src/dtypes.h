  integer, parameter :: wp = selected_real_kind(15, 307)
  real(kind=wp), parameter :: atol = epsilon(1.0_wp)
  real(kind=wp), parameter :: rtol = sqrt(atol)
