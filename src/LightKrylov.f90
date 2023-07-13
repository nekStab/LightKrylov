module LightKrylov
  implicit none
  private

  public :: say_hello
contains
  subroutine say_hello
    print *, "Hello, LightKrylov!"
  end subroutine say_hello
end module LightKrylov
