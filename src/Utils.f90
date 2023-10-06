module Utils
  implicit none
  include "dtypes.h"

  private

  !-------------------------------------------------------
  !-----                                             -----
  !-----     OPTS TYPE OBJECT FOR LINEAR SOLVERS     -----
  !-----                                             -----
  !-------------------------------------------------------

  ! --> Base type.
  type, abstract, public :: abstract_opts
  end type abstract_opts

  ! --> GMRES options.
  type, extends(abstract_opts), public :: gmres_opts
     !> Default dimension of the Krylov subspace.
     integer :: kdim    = 30
     !> Default maximum number of gmres restarts.
     integer :: maxiter = 10
     !> Default tolerance.
     real(kind=wp) :: tol = rtol
     !> Default verbosity control.
     logical :: verbose = .false.
  end type gmres_opts

  ! --> BICGSTAB options.
  type, extends(abstract_opts), public :: bicgstab_opts
     !> Default maximum number of iterations.
     integer :: maxiter = 100
     !> Default tolerance.
     real(kind=wp) :: tol = rtol
     !> Default verbosity control.
     logical :: verbose = .false.
  end type bicgstab_opts

  ! --> Conjugate Gradient options.
  type, extends(abstract_opts), public :: cg_opts
     !> Default maximum number of iterations.
     integer :: maxiter = 100
     !> Default tolerance.
     real(kind=wp) :: tol = rtol
     !> Default verbosity control.
     logical :: verbose = .false.
  end type cg_opts

contains

end module Utils
