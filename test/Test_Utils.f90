module TestUtils
   !! This module provides a set of utilities used for testing
   use LightKrylov
   use TestVector
   use TestMatrices
   implicit none

   private
   public :: get_data, put_data, init_rand

   !---------------------------------------------
   !-----                                   -----
   !-----     INTERFACES FOR DATA TYPES     -----
   !-----                                   -----
   !---------------------------------------------
   
   interface get_data
      module procedure get_data_vec
      module procedure get_data_vec_basis
      module procedure get_data_mat
   end interface

   interface put_data
      module procedure put_data_vec
      module procedure put_data_vec_basis
      module procedure put_data_mat
   end interface

   interface init_rand
      module procedure init_rand_vec
      module procedure init_rand_vec_basis
      module procedure init_rand_mat
   end interface

contains

   !----------------------------------------------------
   !-----                                          -----
   !-----     EXTRACT DATA FROM ABSTRACT TYPES     -----
   !-----                                          -----
   !----------------------------------------------------

   subroutine get_data_vec(vec_out, vec_in)
      !! Utility function to transfer data from a rvector to a real array
      real(kind=wp),          intent(out) :: vec_out(:)
      class(abstract_vector), intent(in)  :: vec_in
      ! internal variables
      integer :: k, kdim
      vec_out = 0.0_wp
      select type (vec_in)
      type is (rvector)
         vec_out(:) = vec_in%data
      end select
      return
   end subroutine get_data_vec

   subroutine get_data_vec_basis(mat_out, basis_in)
      !! Utility function to transfer data from a rvector basis to a real array
      real(kind=wp),          intent(out) :: mat_out(:,:)
      class(abstract_vector), intent(in)  :: basis_in(:)
      ! internal variables
      integer :: k, kdim
      mat_out = 0.0_wp
      select type (basis_in)
      type is (rvector)
         kdim = size(basis_in)
         do k = 1, kdim
            mat_out(:,k) = basis_in(k)%data
         end do
      end select
      return
   end subroutine get_data_vec_basis

   subroutine get_data_mat(mat_out, mat_in)
      !! Utility function to transfer data from a rmatrix to a real array
      real(kind=wp),         intent(out) :: mat_out(:,:)
      class(abstract_linop), intent(in)  :: mat_in
      mat_out = 0.0_wp
      select type (mat_in)
      type is (rmatrix)
         mat_out = mat_in%data
      type is (spd_matrix)
         mat_out = mat_in%data
      end select
      return
   end subroutine get_data_mat

   !------------------------------------------------
   !-----                                      -----
   !-----     INPUT DATA TO ABSTRACT TYPES     -----
   !-----                                      -----
   !------------------------------------------------

   subroutine put_data_vec(vec_out, vec_in)
      !! Utility function to transfer data from a real array to a rvector
      class(abstract_vector), intent(out) :: vec_out
      real(kind=wp),          intent(in)  :: vec_in(:)
      ! internal variables
      select type (vec_out)
      type is (rvector)
         call vec_out%zero()
         vec_out%data = vec_in(:)
      end select
      return
   end subroutine put_data_vec

   subroutine put_data_vec_basis(basis_out, mat_in)
      !! Utility function to transfer data from a real array to a rvector basis
      class(abstract_vector), intent(out) :: basis_out(:)
      real(kind=wp),          intent(in)  :: mat_in(:,:)
      ! internal variables
      integer :: k, kdim
      select type (basis_out)
      type is (rvector)
         call mat_zero(basis_out)
         kdim = size(basis_out)
         do k = 1, kdim
            basis_out(k)%data = mat_in(:,k)
         end do
      end select
      return
   end subroutine put_data_vec_basis

   subroutine put_data_mat(mat_out, mat_in)
      !! Utility function to transfer data from a real array to a rmatrix
      class(abstract_linop), intent(out) :: mat_out
      real(kind=wp),         intent(in)  :: mat_in(:,:)
      select type (mat_out)
      type is (rmatrix)
         mat_out%data = mat_in
      end select
      return
   end subroutine put_data_mat

   !-------------------------------------------------------------
   !-----                                                   -----
   !-----     INITIALIZE ABSTRACT TYPES WITH RANDOM DATA    -----
   !-----                                                   -----
   !-------------------------------------------------------------

   subroutine init_rand_vec(vec)
      !! Utility function to initialize a rvector
      class(abstract_vector), intent(inout) :: vec
      select type (vec)
      type is (rvector)
         call vec%rand()
      end select
      return
   end subroutine init_rand_vec

   subroutine init_rand_vec_basis (basis)
      !! Utility function to initialize a rvector basis
      class(abstract_vector), intent(inout) :: basis (:)
      ! internal variables
      integer :: k, kdim
      select type (basis)
      type is (rvector)
         call mat_zero(basis)
         kdim = size(basis)
         do k = 1, kdim
            call basis(k)%rand()
         end do
      end select
      return
   end subroutine init_rand_vec_basis

   subroutine init_rand_mat(mat)
      !! Utility function to initialize a rmatrix
      class(abstract_linop), intent(inout) :: mat
      select type (mat)
      type is (rmatrix)
         call random_number(mat%data)
      type is (spd_matrix)
         call random_number(mat%data)
         mat%data = matmul(mat%data, transpose(mat%data))
      end select
      return
   end subroutine init_rand_mat

end module TestUtils