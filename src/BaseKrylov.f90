module lightkrylov_BaseKrylov
    use iso_fortran_env
    use stdlib_optval, only: optval
    use stdlib_linalg, only: eye
    use LightKrylov_constants
    use LightKrylov_Logger
    use LightKrylov_Utils
    use LightKrylov_AbstractVectors
    use LightKrylov_AbstractLinops
    implicit none
    private
    
    character*128, parameter :: this_module = 'LightKrylov_BaseKrylov'

    public :: qr
    public :: apply_permutation_matrix
    public :: apply_inverse_permutation_matrix
    public :: arnoldi
    public :: initialize_krylov_subspace
    public :: orthogonalize_against_basis
    public :: lanczos_bidiagonalization
    public :: lanczos_tridiagonalization
    public :: krylov_schur
    public :: double_gram_schmidt_step
    public :: eigvals_select_sp
    public :: eigvals_select_dp

    interface qr
        module procedure qr_no_pivoting_rsp
        module procedure qr_with_pivoting_rsp
        module procedure qr_no_pivoting_rdp
        module procedure qr_with_pivoting_rdp
        module procedure qr_no_pivoting_csp
        module procedure qr_with_pivoting_csp
        module procedure qr_no_pivoting_cdp
        module procedure qr_with_pivoting_cdp
    end interface

    interface swap_columns
        module procedure swap_columns_rsp
        module procedure swap_columns_rdp
        module procedure swap_columns_csp
        module procedure swap_columns_cdp
    end interface

    interface apply_permutation_matrix
        module procedure apply_permutation_matrix_rsp
        module procedure apply_permutation_matrix_array_rsp
        module procedure apply_permutation_matrix_rdp
        module procedure apply_permutation_matrix_array_rdp
        module procedure apply_permutation_matrix_csp
        module procedure apply_permutation_matrix_array_csp
        module procedure apply_permutation_matrix_cdp
        module procedure apply_permutation_matrix_array_cdp
    end interface

    interface apply_inverse_permutation_matrix
        module procedure apply_inverse_permutation_matrix_rsp
        module procedure apply_inverse_permutation_matrix_array_rsp
        module procedure apply_inverse_permutation_matrix_rdp
        module procedure apply_inverse_permutation_matrix_array_rdp
        module procedure apply_inverse_permutation_matrix_csp
        module procedure apply_inverse_permutation_matrix_array_csp
        module procedure apply_inverse_permutation_matrix_cdp
        module procedure apply_inverse_permutation_matrix_array_cdp
    end interface

    interface arnoldi
        module procedure arnoldi_rsp
        module procedure arnoldi_rdp
        module procedure arnoldi_csp
        module procedure arnoldi_cdp
    end interface

    interface initialize_krylov_subspace
        module procedure initialize_krylov_subspace_rsp
        module procedure initialize_krylov_subspace_rdp
        module procedure initialize_krylov_subspace_csp
        module procedure initialize_krylov_subspace_cdp
    end interface

    interface orthogonalize_against_basis
        module procedure orthogonalize_vector_against_basis_rsp
        module procedure orthogonalize_basis_against_basis_rsp
        module procedure orthogonalize_vector_against_basis_rdp
        module procedure orthogonalize_basis_against_basis_rdp
        module procedure orthogonalize_vector_against_basis_csp
        module procedure orthogonalize_basis_against_basis_csp
        module procedure orthogonalize_vector_against_basis_cdp
        module procedure orthogonalize_basis_against_basis_cdp
    end interface

    interface double_gram_schmidt_step
        module procedure DGS_vector_against_basis_rsp
        module procedure DGS_basis_against_basis_rsp
        module procedure DGS_vector_against_basis_rdp
        module procedure DGS_basis_against_basis_rdp
        module procedure DGS_vector_against_basis_csp
        module procedure DGS_basis_against_basis_csp
        module procedure DGS_vector_against_basis_cdp
        module procedure DGS_basis_against_basis_cdp
    end interface

    interface lanczos_tridiagonalization
        module procedure lanczos_tridiagonalization_rsp
        module procedure lanczos_tridiagonalization_rdp
        module procedure lanczos_tridiagonalization_csp
        module procedure lanczos_tridiagonalization_cdp
    end interface

    interface lanczos_bidiagonalization
        module procedure lanczos_bidiagonalization_rsp
        module procedure lanczos_bidiagonalization_rdp
        module procedure lanczos_bidiagonalization_csp
        module procedure lanczos_bidiagonalization_cdp
    end interface

    interface krylov_schur
        module procedure krylov_schur_rsp
        module procedure krylov_schur_rdp
        module procedure krylov_schur_csp
        module procedure krylov_schur_cdp
    end interface

    !----------------------------------------------------------
    !-----     ABSTRACT EIGENVALUE SELECTOR INTERFACE     -----
    !----------------------------------------------------------

    abstract interface
        function eigvals_select_sp(lambda) result(selected)
            import sp
            complex(sp), intent(in) :: lambda(:)
            logical                       :: selected(size(lambda))
        end function eigvals_select_sp
        function eigvals_select_dp(lambda) result(selected)
            import dp
            complex(dp), intent(in) :: lambda(:)
            logical                       :: selected(size(lambda))
        end function eigvals_select_dp
    end interface

contains

    !-------------------------------------
    !-----     UTILITY FUNCTIONS     -----
    !-------------------------------------

    subroutine initialize_krylov_subspace_rsp(X, X0)
        class(abstract_vector_rsp), intent(inout) :: X(:)
        class(abstract_vector_rsp), optional, intent(in) :: X0(:)

        ! Internal variables.
        class(abstract_vector_rsp), allocatable :: Q(:)
        real(sp), allocatable :: R(:, :)
        integer, allocatable :: perm(:)
        integer :: i, p, info

        ! Zero-out X.
        call zero_basis(X)

        ! Deals with optional args.
        if(present(X0)) then
            p = size(X0)

            ! Initialize.
            allocate(Q(1:p), source=X0(1:p))
            do i = 1, p
                call Q(i)%zero()
                call Q(i)%add(X0(i))
            enddo
    
            ! Orthonormalize.
            allocate(R(p, p)) ; R = zero_rsp
            allocate(perm(p)) ; perm = 0
            call qr(Q, R, perm, info)
            do i = 1, p
                call X(i)%add(Q(i))
            enddo
        endif

        return
    end subroutine initialize_krylov_subspace_rsp

    subroutine orthogonalize_vector_against_basis_rsp(y, X, info, if_chk_orthonormal, beta)
      !! Orthonormalizes the `abstract_vector` `y` against a basis `X` of `abstract_vector`.
      class(abstract_vector_rsp), intent(inout) :: y
      !! Input `abstract_vector` to orthogonalize
      class(abstract_vector_rsp), intent(in)    :: X(:)
      !! Input `abstract_vector` basis to orthogonalize against
      integer, intent(out) :: info
      !! Information flag.
      logical,                          optional, intent(in)    :: if_chk_orthonormal
      logical                                                   :: chk_X_orthonormality
      !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
      real(sp),                         optional, intent(out)   :: beta(:)
      !! Projection coefficients if requested

      ! internals
      real(sp) :: proj_coefficients(size(X))

      info = 0

      ! optional input argument
      chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

      ! check for zero vector
      if (y%norm() < atol_sp) info = 1

      if (chk_X_orthonormality) then
         block 
            real(sp), dimension(size(X), size(X)) :: G
            call innerprod(G, X, X)
            if (abs(G(size(X),size(X))) < rtol_sp) then
               ! The last vector in X is zero, it does not impact orthogonalisation
               info = -2
            else if (norm2(abs(G - eye(size(X)))) > rtol_sp) then
               ! The basis is not orthonormal. Cannot orthonormalize.
               info = -1
               return
            end if
         end block
      end if

      ! orthogonalize
      call innerprod(proj_coefficients, X, y)
      block
         class(abstract_vector_rsp), allocatable :: proj
         call linear_combination(proj, X, proj_coefficients)
         call y%sub(proj)
      end block

      if (present(beta)) then
         ! check size
         call assert_shape(beta, shape(proj_coefficients), 'orthogonalize_vector_against_basis_rsp', 'beta')
         beta = proj_coefficients
      end if
      
      return
    end subroutine orthogonalize_vector_against_basis_rsp

    subroutine orthogonalize_basis_against_basis_rsp(Y, X, info, if_chk_orthonormal, beta)
      !! Orthonormalizes the `abstract_vector` basis `Y` against a basis `X` of `abstract_vector`.
      class(abstract_vector_rsp), intent(inout) :: Y(:)
      !! Input `abstract_vector` basis to orthogonalize
      class(abstract_vector_rsp), intent(in)    :: X(:)
      !! Input `abstract_vector` basis to orthogonalize against
      integer, intent(out) :: info
      !! Information flag.
      logical,                          optional, intent(in)    :: if_chk_orthonormal
      logical                                                   :: chk_X_orthonormality
      !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
      real(sp),                         optional, intent(out)   :: beta(:,:)
      !! Projection coefficients if requested

      ! internals
      real(sp) :: proj_coefficients(size(X), size(Y))
      integer :: i

      info = 0

      ! optional input argument
      chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

      ! check for zero vector
      do i = 1, size(Y)
         if ( Y(i)%norm() < atol_sp) info = i
      end do

      if (chk_X_orthonormality) then
         block 
            real(sp), dimension(size(X), size(X)) :: G
            call innerprod(G, X, X)
            if (abs(G(size(X),size(X))) < rtol_sp) then
               ! The last vector in X is zero, it does not impact orthogonalisation
               info = -2
            else if (norm2(abs(G - eye(size(X)))) > rtol_sp) then
               ! The basis is not orthonormal. Cannot orthonormalize.
               info = -1
               return
            end if
         end block
      end if
      
      ! orthogonalize
      call innerprod(proj_coefficients, X, Y)
      block
         class(abstract_vector_rsp), allocatable :: proj(:)
         call linear_combination(proj, X, proj_coefficients)
         call axpby_basis(Y, one_rsp, proj, -one_rsp)
      end block

      if (present(beta)) then
         ! check size
         call assert_shape(beta, shape(proj_coefficients), 'orthogonalize_basis_against_basis_rsp', 'beta')
         beta = proj_coefficients
      end if
      
      return
    end subroutine orthogonalize_basis_against_basis_rsp

    subroutine DGS_vector_against_basis_rsp(y, X, info, if_chk_orthonormal, beta)
      !! Computes one step of the double Gram-Schmidt orthogonalization process of the
      !! `abstract_vector` `y` against the `abstract_vector` basis `X`
      class(abstract_vector_rsp), intent(inout) :: y
      !! Input `abstract_vector` to orthogonalize
      class(abstract_vector_rsp), intent(in)    :: X(:)
      !! Input `abstract_vector` basis to orthogonalize against
      integer, intent(out) :: info
      !! Information flag.
      logical,                          optional, intent(in)    :: if_chk_orthonormal
      logical                                                   :: chk_X_orthonormality
      !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
      real(sp),                         optional, intent(out)   :: beta(:)
      !! Projection coefficients if requested

      ! internals
      real(sp), dimension(size(X)) :: proj_coefficients, wrk

      ! optional input argument
      chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

      info = 0

      proj_coefficients = zero_rsp; wrk = zero_rsp

      ! Orthogonalize vector y w.r.t. to Krylov basis X in two passes of GS.
      ! first pass
      call orthogonalize_against_basis(y, X, info, if_chk_orthonormal=.false., beta=proj_coefficients)
      call check_info(info, 'orthogonalize_against_basis', module=this_module, procedure='DGS_vector_against_basis, first pass')
      ! second pass
      call orthogonalize_against_basis(y, X, info, if_chk_orthonormal=.false., beta=wrk)
      call check_info(info, 'orthogonalize_against_basis', module=this_module, procedure='DGS_vector_against_basis_rsp, second&
          & pass')
      ! combine passes
      proj_coefficients = proj_coefficients + wrk

      if (present(beta)) then
         ! check size
         call assert_shape(beta, shape(proj_coefficients), 'DGS_vector_against_basis_rsp', 'beta')
         beta = proj_coefficients
      end if

    end subroutine DGS_vector_against_basis_rsp

    subroutine DGS_basis_against_basis_rsp(y, X, info, if_chk_orthonormal, beta)
      !! Computes one step of the double Gram-Schmidt orthogonalization process of the
      !! `abstract_vector` `y` against the `abstract_vector` basis `X`
      class(abstract_vector_rsp), intent(inout) :: Y(:)
      !! Input `abstract_vector` basis to orthogonalize
      class(abstract_vector_rsp), intent(in)    :: X(:)
      !! Input `abstract_vector` basis to orthogonalize against
      integer, intent(out) :: info
      !! Information flag.
      logical,                          optional, intent(in)    :: if_chk_orthonormal
      logical                                                   :: chk_X_orthonormality
      !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
      real(sp),                         optional, intent(out)   :: beta(:,:)
      !! Projection coefficients if requested

      ! internals
      real(sp), dimension(size(X),size(Y)) :: proj_coefficients, wrk

      ! optional input argument
      chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

      info = 0

      proj_coefficients = zero_rsp; wrk = zero_rsp

      ! Orthogonalize Krylov basis Y w.r.t. to Krylov basis X in two passes of GS.
      ! first pass
      call orthogonalize_against_basis(Y, X, info, if_chk_orthonormal=.false., beta=proj_coefficients)
      call check_info(info, 'orthogonalize_against_basis', module=this_module, procedure='DGS_basis_against_basis_rsp, first pass')
      ! second pass
      call orthogonalize_against_basis(Y, X, info, if_chk_orthonormal=.false., beta=wrk)
      call check_info(info, 'orthogonalize_against_basis', module=this_module, procedure='DGS_basis_against_basis_rsp, second pass')
      ! combine passes
      proj_coefficients = proj_coefficients + wrk

      if (present(beta)) then
         ! check size
         call assert_shape(beta, shape(proj_coefficients), 'DGS_basis_against_basis_rsp', 'beta')
         beta = proj_coefficients
      end if

    end subroutine DGS_basis_against_basis_rsp  

    subroutine initialize_krylov_subspace_rdp(X, X0)
        class(abstract_vector_rdp), intent(inout) :: X(:)
        class(abstract_vector_rdp), optional, intent(in) :: X0(:)

        ! Internal variables.
        class(abstract_vector_rdp), allocatable :: Q(:)
        real(dp), allocatable :: R(:, :)
        integer, allocatable :: perm(:)
        integer :: i, p, info

        ! Zero-out X.
        call zero_basis(X)

        ! Deals with optional args.
        if(present(X0)) then
            p = size(X0)

            ! Initialize.
            allocate(Q(1:p), source=X0(1:p))
            do i = 1, p
                call Q(i)%zero()
                call Q(i)%add(X0(i))
            enddo
    
            ! Orthonormalize.
            allocate(R(p, p)) ; R = zero_rdp
            allocate(perm(p)) ; perm = 0
            call qr(Q, R, perm, info)
            do i = 1, p
                call X(i)%add(Q(i))
            enddo
        endif

        return
    end subroutine initialize_krylov_subspace_rdp

    subroutine orthogonalize_vector_against_basis_rdp(y, X, info, if_chk_orthonormal, beta)
      !! Orthonormalizes the `abstract_vector` `y` against a basis `X` of `abstract_vector`.
      class(abstract_vector_rdp), intent(inout) :: y
      !! Input `abstract_vector` to orthogonalize
      class(abstract_vector_rdp), intent(in)    :: X(:)
      !! Input `abstract_vector` basis to orthogonalize against
      integer, intent(out) :: info
      !! Information flag.
      logical,                          optional, intent(in)    :: if_chk_orthonormal
      logical                                                   :: chk_X_orthonormality
      !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
      real(dp),                         optional, intent(out)   :: beta(:)
      !! Projection coefficients if requested

      ! internals
      real(dp) :: proj_coefficients(size(X))

      info = 0

      ! optional input argument
      chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

      ! check for zero vector
      if (y%norm() < atol_dp) info = 1

      if (chk_X_orthonormality) then
         block 
            real(dp), dimension(size(X), size(X)) :: G
            call innerprod(G, X, X)
            if (abs(G(size(X),size(X))) < rtol_dp) then
               ! The last vector in X is zero, it does not impact orthogonalisation
               info = -2
            else if (norm2(abs(G - eye(size(X)))) > rtol_dp) then
               ! The basis is not orthonormal. Cannot orthonormalize.
               info = -1
               return
            end if
         end block
      end if

      ! orthogonalize
      call innerprod(proj_coefficients, X, y)
      block
         class(abstract_vector_rdp), allocatable :: proj
         call linear_combination(proj, X, proj_coefficients)
         call y%sub(proj)
      end block

      if (present(beta)) then
         ! check size
         call assert_shape(beta, shape(proj_coefficients), 'orthogonalize_vector_against_basis_rdp', 'beta')
         beta = proj_coefficients
      end if
      
      return
    end subroutine orthogonalize_vector_against_basis_rdp

    subroutine orthogonalize_basis_against_basis_rdp(Y, X, info, if_chk_orthonormal, beta)
      !! Orthonormalizes the `abstract_vector` basis `Y` against a basis `X` of `abstract_vector`.
      class(abstract_vector_rdp), intent(inout) :: Y(:)
      !! Input `abstract_vector` basis to orthogonalize
      class(abstract_vector_rdp), intent(in)    :: X(:)
      !! Input `abstract_vector` basis to orthogonalize against
      integer, intent(out) :: info
      !! Information flag.
      logical,                          optional, intent(in)    :: if_chk_orthonormal
      logical                                                   :: chk_X_orthonormality
      !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
      real(dp),                         optional, intent(out)   :: beta(:,:)
      !! Projection coefficients if requested

      ! internals
      real(dp) :: proj_coefficients(size(X), size(Y))
      integer :: i

      info = 0

      ! optional input argument
      chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

      ! check for zero vector
      do i = 1, size(Y)
         if ( Y(i)%norm() < atol_dp) info = i
      end do

      if (chk_X_orthonormality) then
         block 
            real(dp), dimension(size(X), size(X)) :: G
            call innerprod(G, X, X)
            if (abs(G(size(X),size(X))) < rtol_dp) then
               ! The last vector in X is zero, it does not impact orthogonalisation
               info = -2
            else if (norm2(abs(G - eye(size(X)))) > rtol_dp) then
               ! The basis is not orthonormal. Cannot orthonormalize.
               info = -1
               return
            end if
         end block
      end if
      
      ! orthogonalize
      call innerprod(proj_coefficients, X, Y)
      block
         class(abstract_vector_rdp), allocatable :: proj(:)
         call linear_combination(proj, X, proj_coefficients)
         call axpby_basis(Y, one_rdp, proj, -one_rdp)
      end block

      if (present(beta)) then
         ! check size
         call assert_shape(beta, shape(proj_coefficients), 'orthogonalize_basis_against_basis_rdp', 'beta')
         beta = proj_coefficients
      end if
      
      return
    end subroutine orthogonalize_basis_against_basis_rdp

    subroutine DGS_vector_against_basis_rdp(y, X, info, if_chk_orthonormal, beta)
      !! Computes one step of the double Gram-Schmidt orthogonalization process of the
      !! `abstract_vector` `y` against the `abstract_vector` basis `X`
      class(abstract_vector_rdp), intent(inout) :: y
      !! Input `abstract_vector` to orthogonalize
      class(abstract_vector_rdp), intent(in)    :: X(:)
      !! Input `abstract_vector` basis to orthogonalize against
      integer, intent(out) :: info
      !! Information flag.
      logical,                          optional, intent(in)    :: if_chk_orthonormal
      logical                                                   :: chk_X_orthonormality
      !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
      real(dp),                         optional, intent(out)   :: beta(:)
      !! Projection coefficients if requested

      ! internals
      real(dp), dimension(size(X)) :: proj_coefficients, wrk

      ! optional input argument
      chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

      info = 0

      proj_coefficients = zero_rdp; wrk = zero_rdp

      ! Orthogonalize vector y w.r.t. to Krylov basis X in two passes of GS.
      ! first pass
      call orthogonalize_against_basis(y, X, info, if_chk_orthonormal=.false., beta=proj_coefficients)
      call check_info(info, 'orthogonalize_against_basis', module=this_module, procedure='DGS_vector_against_basis, first pass')
      ! second pass
      call orthogonalize_against_basis(y, X, info, if_chk_orthonormal=.false., beta=wrk)
      call check_info(info, 'orthogonalize_against_basis', module=this_module, procedure='DGS_vector_against_basis_rdp, second&
          & pass')
      ! combine passes
      proj_coefficients = proj_coefficients + wrk

      if (present(beta)) then
         ! check size
         call assert_shape(beta, shape(proj_coefficients), 'DGS_vector_against_basis_rdp', 'beta')
         beta = proj_coefficients
      end if

    end subroutine DGS_vector_against_basis_rdp

    subroutine DGS_basis_against_basis_rdp(y, X, info, if_chk_orthonormal, beta)
      !! Computes one step of the double Gram-Schmidt orthogonalization process of the
      !! `abstract_vector` `y` against the `abstract_vector` basis `X`
      class(abstract_vector_rdp), intent(inout) :: Y(:)
      !! Input `abstract_vector` basis to orthogonalize
      class(abstract_vector_rdp), intent(in)    :: X(:)
      !! Input `abstract_vector` basis to orthogonalize against
      integer, intent(out) :: info
      !! Information flag.
      logical,                          optional, intent(in)    :: if_chk_orthonormal
      logical                                                   :: chk_X_orthonormality
      !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
      real(dp),                         optional, intent(out)   :: beta(:,:)
      !! Projection coefficients if requested

      ! internals
      real(dp), dimension(size(X),size(Y)) :: proj_coefficients, wrk

      ! optional input argument
      chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

      info = 0

      proj_coefficients = zero_rdp; wrk = zero_rdp

      ! Orthogonalize Krylov basis Y w.r.t. to Krylov basis X in two passes of GS.
      ! first pass
      call orthogonalize_against_basis(Y, X, info, if_chk_orthonormal=.false., beta=proj_coefficients)
      call check_info(info, 'orthogonalize_against_basis', module=this_module, procedure='DGS_basis_against_basis_rdp, first pass')
      ! second pass
      call orthogonalize_against_basis(Y, X, info, if_chk_orthonormal=.false., beta=wrk)
      call check_info(info, 'orthogonalize_against_basis', module=this_module, procedure='DGS_basis_against_basis_rdp, second pass')
      ! combine passes
      proj_coefficients = proj_coefficients + wrk

      if (present(beta)) then
         ! check size
         call assert_shape(beta, shape(proj_coefficients), 'DGS_basis_against_basis_rdp', 'beta')
         beta = proj_coefficients
      end if

    end subroutine DGS_basis_against_basis_rdp  

    subroutine initialize_krylov_subspace_csp(X, X0)
        class(abstract_vector_csp), intent(inout) :: X(:)
        class(abstract_vector_csp), optional, intent(in) :: X0(:)

        ! Internal variables.
        class(abstract_vector_csp), allocatable :: Q(:)
        complex(sp), allocatable :: R(:, :)
        integer, allocatable :: perm(:)
        integer :: i, p, info

        ! Zero-out X.
        call zero_basis(X)

        ! Deals with optional args.
        if(present(X0)) then
            p = size(X0)

            ! Initialize.
            allocate(Q(1:p), source=X0(1:p))
            do i = 1, p
                call Q(i)%zero()
                call Q(i)%add(X0(i))
            enddo
    
            ! Orthonormalize.
            allocate(R(p, p)) ; R = zero_rsp
            allocate(perm(p)) ; perm = 0
            call qr(Q, R, perm, info)
            do i = 1, p
                call X(i)%add(Q(i))
            enddo
        endif

        return
    end subroutine initialize_krylov_subspace_csp

    subroutine orthogonalize_vector_against_basis_csp(y, X, info, if_chk_orthonormal, beta)
      !! Orthonormalizes the `abstract_vector` `y` against a basis `X` of `abstract_vector`.
      class(abstract_vector_csp), intent(inout) :: y
      !! Input `abstract_vector` to orthogonalize
      class(abstract_vector_csp), intent(in)    :: X(:)
      !! Input `abstract_vector` basis to orthogonalize against
      integer, intent(out) :: info
      !! Information flag.
      logical,                          optional, intent(in)    :: if_chk_orthonormal
      logical                                                   :: chk_X_orthonormality
      !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
      complex(sp),                         optional, intent(out)   :: beta(:)
      !! Projection coefficients if requested

      ! internals
      complex(sp) :: proj_coefficients(size(X))

      info = 0

      ! optional input argument
      chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

      ! check for zero vector
      if (y%norm() < atol_sp) info = 1

      if (chk_X_orthonormality) then
         block 
            complex(sp), dimension(size(X), size(X)) :: G
            call innerprod(G, X, X)
            if (abs(G(size(X),size(X))) < rtol_sp) then
               ! The last vector in X is zero, it does not impact orthogonalisation
               info = -2
            else if (norm2(abs(G - eye(size(X)))) > rtol_sp) then
               ! The basis is not orthonormal. Cannot orthonormalize.
               info = -1
               return
            end if
         end block
      end if

      ! orthogonalize
      call innerprod(proj_coefficients, X, y)
      block
         class(abstract_vector_csp), allocatable :: proj
         call linear_combination(proj, X, proj_coefficients)
         call y%sub(proj)
      end block

      if (present(beta)) then
         ! check size
         call assert_shape(beta, shape(proj_coefficients), 'orthogonalize_vector_against_basis_csp', 'beta')
         beta = proj_coefficients
      end if
      
      return
    end subroutine orthogonalize_vector_against_basis_csp

    subroutine orthogonalize_basis_against_basis_csp(Y, X, info, if_chk_orthonormal, beta)
      !! Orthonormalizes the `abstract_vector` basis `Y` against a basis `X` of `abstract_vector`.
      class(abstract_vector_csp), intent(inout) :: Y(:)
      !! Input `abstract_vector` basis to orthogonalize
      class(abstract_vector_csp), intent(in)    :: X(:)
      !! Input `abstract_vector` basis to orthogonalize against
      integer, intent(out) :: info
      !! Information flag.
      logical,                          optional, intent(in)    :: if_chk_orthonormal
      logical                                                   :: chk_X_orthonormality
      !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
      complex(sp),                         optional, intent(out)   :: beta(:,:)
      !! Projection coefficients if requested

      ! internals
      complex(sp) :: proj_coefficients(size(X), size(Y))
      integer :: i

      info = 0

      ! optional input argument
      chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

      ! check for zero vector
      do i = 1, size(Y)
         if ( Y(i)%norm() < atol_sp) info = i
      end do

      if (chk_X_orthonormality) then
         block 
            complex(sp), dimension(size(X), size(X)) :: G
            call innerprod(G, X, X)
            if (abs(G(size(X),size(X))) < rtol_sp) then
               ! The last vector in X is zero, it does not impact orthogonalisation
               info = -2
            else if (norm2(abs(G - eye(size(X)))) > rtol_sp) then
               ! The basis is not orthonormal. Cannot orthonormalize.
               info = -1
               return
            end if
         end block
      end if
      
      ! orthogonalize
      call innerprod(proj_coefficients, X, Y)
      block
         class(abstract_vector_csp), allocatable :: proj(:)
         call linear_combination(proj, X, proj_coefficients)
         call axpby_basis(Y, one_csp, proj, -one_csp)
      end block

      if (present(beta)) then
         ! check size
         call assert_shape(beta, shape(proj_coefficients), 'orthogonalize_basis_against_basis_csp', 'beta')
         beta = proj_coefficients
      end if
      
      return
    end subroutine orthogonalize_basis_against_basis_csp

    subroutine DGS_vector_against_basis_csp(y, X, info, if_chk_orthonormal, beta)
      !! Computes one step of the double Gram-Schmidt orthogonalization process of the
      !! `abstract_vector` `y` against the `abstract_vector` basis `X`
      class(abstract_vector_csp), intent(inout) :: y
      !! Input `abstract_vector` to orthogonalize
      class(abstract_vector_csp), intent(in)    :: X(:)
      !! Input `abstract_vector` basis to orthogonalize against
      integer, intent(out) :: info
      !! Information flag.
      logical,                          optional, intent(in)    :: if_chk_orthonormal
      logical                                                   :: chk_X_orthonormality
      !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
      complex(sp),                         optional, intent(out)   :: beta(:)
      !! Projection coefficients if requested

      ! internals
      complex(sp), dimension(size(X)) :: proj_coefficients, wrk

      ! optional input argument
      chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

      info = 0

      proj_coefficients = zero_csp; wrk = zero_csp

      ! Orthogonalize vector y w.r.t. to Krylov basis X in two passes of GS.
      ! first pass
      call orthogonalize_against_basis(y, X, info, if_chk_orthonormal=.false., beta=proj_coefficients)
      call check_info(info, 'orthogonalize_against_basis', module=this_module, procedure='DGS_vector_against_basis, first pass')
      ! second pass
      call orthogonalize_against_basis(y, X, info, if_chk_orthonormal=.false., beta=wrk)
      call check_info(info, 'orthogonalize_against_basis', module=this_module, procedure='DGS_vector_against_basis_csp, second&
          & pass')
      ! combine passes
      proj_coefficients = proj_coefficients + wrk

      if (present(beta)) then
         ! check size
         call assert_shape(beta, shape(proj_coefficients), 'DGS_vector_against_basis_csp', 'beta')
         beta = proj_coefficients
      end if

    end subroutine DGS_vector_against_basis_csp

    subroutine DGS_basis_against_basis_csp(y, X, info, if_chk_orthonormal, beta)
      !! Computes one step of the double Gram-Schmidt orthogonalization process of the
      !! `abstract_vector` `y` against the `abstract_vector` basis `X`
      class(abstract_vector_csp), intent(inout) :: Y(:)
      !! Input `abstract_vector` basis to orthogonalize
      class(abstract_vector_csp), intent(in)    :: X(:)
      !! Input `abstract_vector` basis to orthogonalize against
      integer, intent(out) :: info
      !! Information flag.
      logical,                          optional, intent(in)    :: if_chk_orthonormal
      logical                                                   :: chk_X_orthonormality
      !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
      complex(sp),                         optional, intent(out)   :: beta(:,:)
      !! Projection coefficients if requested

      ! internals
      complex(sp), dimension(size(X),size(Y)) :: proj_coefficients, wrk

      ! optional input argument
      chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

      info = 0

      proj_coefficients = zero_csp; wrk = zero_csp

      ! Orthogonalize Krylov basis Y w.r.t. to Krylov basis X in two passes of GS.
      ! first pass
      call orthogonalize_against_basis(Y, X, info, if_chk_orthonormal=.false., beta=proj_coefficients)
      call check_info(info, 'orthogonalize_against_basis', module=this_module, procedure='DGS_basis_against_basis_csp, first pass')
      ! second pass
      call orthogonalize_against_basis(Y, X, info, if_chk_orthonormal=.false., beta=wrk)
      call check_info(info, 'orthogonalize_against_basis', module=this_module, procedure='DGS_basis_against_basis_csp, second pass')
      ! combine passes
      proj_coefficients = proj_coefficients + wrk

      if (present(beta)) then
         ! check size
         call assert_shape(beta, shape(proj_coefficients), 'DGS_basis_against_basis_csp', 'beta')
         beta = proj_coefficients
      end if

    end subroutine DGS_basis_against_basis_csp  

    subroutine initialize_krylov_subspace_cdp(X, X0)
        class(abstract_vector_cdp), intent(inout) :: X(:)
        class(abstract_vector_cdp), optional, intent(in) :: X0(:)

        ! Internal variables.
        class(abstract_vector_cdp), allocatable :: Q(:)
        complex(dp), allocatable :: R(:, :)
        integer, allocatable :: perm(:)
        integer :: i, p, info

        ! Zero-out X.
        call zero_basis(X)

        ! Deals with optional args.
        if(present(X0)) then
            p = size(X0)

            ! Initialize.
            allocate(Q(1:p), source=X0(1:p))
            do i = 1, p
                call Q(i)%zero()
                call Q(i)%add(X0(i))
            enddo
    
            ! Orthonormalize.
            allocate(R(p, p)) ; R = zero_rdp
            allocate(perm(p)) ; perm = 0
            call qr(Q, R, perm, info)
            do i = 1, p
                call X(i)%add(Q(i))
            enddo
        endif

        return
    end subroutine initialize_krylov_subspace_cdp

    subroutine orthogonalize_vector_against_basis_cdp(y, X, info, if_chk_orthonormal, beta)
      !! Orthonormalizes the `abstract_vector` `y` against a basis `X` of `abstract_vector`.
      class(abstract_vector_cdp), intent(inout) :: y
      !! Input `abstract_vector` to orthogonalize
      class(abstract_vector_cdp), intent(in)    :: X(:)
      !! Input `abstract_vector` basis to orthogonalize against
      integer, intent(out) :: info
      !! Information flag.
      logical,                          optional, intent(in)    :: if_chk_orthonormal
      logical                                                   :: chk_X_orthonormality
      !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
      complex(dp),                         optional, intent(out)   :: beta(:)
      !! Projection coefficients if requested

      ! internals
      complex(dp) :: proj_coefficients(size(X))

      info = 0

      ! optional input argument
      chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

      ! check for zero vector
      if (y%norm() < atol_dp) info = 1

      if (chk_X_orthonormality) then
         block 
            complex(dp), dimension(size(X), size(X)) :: G
            call innerprod(G, X, X)
            if (abs(G(size(X),size(X))) < rtol_dp) then
               ! The last vector in X is zero, it does not impact orthogonalisation
               info = -2
            else if (norm2(abs(G - eye(size(X)))) > rtol_dp) then
               ! The basis is not orthonormal. Cannot orthonormalize.
               info = -1
               return
            end if
         end block
      end if

      ! orthogonalize
      call innerprod(proj_coefficients, X, y)
      block
         class(abstract_vector_cdp), allocatable :: proj
         call linear_combination(proj, X, proj_coefficients)
         call y%sub(proj)
      end block

      if (present(beta)) then
         ! check size
         call assert_shape(beta, shape(proj_coefficients), 'orthogonalize_vector_against_basis_cdp', 'beta')
         beta = proj_coefficients
      end if
      
      return
    end subroutine orthogonalize_vector_against_basis_cdp

    subroutine orthogonalize_basis_against_basis_cdp(Y, X, info, if_chk_orthonormal, beta)
      !! Orthonormalizes the `abstract_vector` basis `Y` against a basis `X` of `abstract_vector`.
      class(abstract_vector_cdp), intent(inout) :: Y(:)
      !! Input `abstract_vector` basis to orthogonalize
      class(abstract_vector_cdp), intent(in)    :: X(:)
      !! Input `abstract_vector` basis to orthogonalize against
      integer, intent(out) :: info
      !! Information flag.
      logical,                          optional, intent(in)    :: if_chk_orthonormal
      logical                                                   :: chk_X_orthonormality
      !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
      complex(dp),                         optional, intent(out)   :: beta(:,:)
      !! Projection coefficients if requested

      ! internals
      complex(dp) :: proj_coefficients(size(X), size(Y))
      integer :: i

      info = 0

      ! optional input argument
      chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

      ! check for zero vector
      do i = 1, size(Y)
         if ( Y(i)%norm() < atol_dp) info = i
      end do

      if (chk_X_orthonormality) then
         block 
            complex(dp), dimension(size(X), size(X)) :: G
            call innerprod(G, X, X)
            if (abs(G(size(X),size(X))) < rtol_dp) then
               ! The last vector in X is zero, it does not impact orthogonalisation
               info = -2
            else if (norm2(abs(G - eye(size(X)))) > rtol_dp) then
               ! The basis is not orthonormal. Cannot orthonormalize.
               info = -1
               return
            end if
         end block
      end if
      
      ! orthogonalize
      call innerprod(proj_coefficients, X, Y)
      block
         class(abstract_vector_cdp), allocatable :: proj(:)
         call linear_combination(proj, X, proj_coefficients)
         call axpby_basis(Y, one_cdp, proj, -one_cdp)
      end block

      if (present(beta)) then
         ! check size
         call assert_shape(beta, shape(proj_coefficients), 'orthogonalize_basis_against_basis_cdp', 'beta')
         beta = proj_coefficients
      end if
      
      return
    end subroutine orthogonalize_basis_against_basis_cdp

    subroutine DGS_vector_against_basis_cdp(y, X, info, if_chk_orthonormal, beta)
      !! Computes one step of the double Gram-Schmidt orthogonalization process of the
      !! `abstract_vector` `y` against the `abstract_vector` basis `X`
      class(abstract_vector_cdp), intent(inout) :: y
      !! Input `abstract_vector` to orthogonalize
      class(abstract_vector_cdp), intent(in)    :: X(:)
      !! Input `abstract_vector` basis to orthogonalize against
      integer, intent(out) :: info
      !! Information flag.
      logical,                          optional, intent(in)    :: if_chk_orthonormal
      logical                                                   :: chk_X_orthonormality
      !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
      complex(dp),                         optional, intent(out)   :: beta(:)
      !! Projection coefficients if requested

      ! internals
      complex(dp), dimension(size(X)) :: proj_coefficients, wrk

      ! optional input argument
      chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

      info = 0

      proj_coefficients = zero_cdp; wrk = zero_cdp

      ! Orthogonalize vector y w.r.t. to Krylov basis X in two passes of GS.
      ! first pass
      call orthogonalize_against_basis(y, X, info, if_chk_orthonormal=.false., beta=proj_coefficients)
      call check_info(info, 'orthogonalize_against_basis', module=this_module, procedure='DGS_vector_against_basis, first pass')
      ! second pass
      call orthogonalize_against_basis(y, X, info, if_chk_orthonormal=.false., beta=wrk)
      call check_info(info, 'orthogonalize_against_basis', module=this_module, procedure='DGS_vector_against_basis_cdp, second&
          & pass')
      ! combine passes
      proj_coefficients = proj_coefficients + wrk

      if (present(beta)) then
         ! check size
         call assert_shape(beta, shape(proj_coefficients), 'DGS_vector_against_basis_cdp', 'beta')
         beta = proj_coefficients
      end if

    end subroutine DGS_vector_against_basis_cdp

    subroutine DGS_basis_against_basis_cdp(y, X, info, if_chk_orthonormal, beta)
      !! Computes one step of the double Gram-Schmidt orthogonalization process of the
      !! `abstract_vector` `y` against the `abstract_vector` basis `X`
      class(abstract_vector_cdp), intent(inout) :: Y(:)
      !! Input `abstract_vector` basis to orthogonalize
      class(abstract_vector_cdp), intent(in)    :: X(:)
      !! Input `abstract_vector` basis to orthogonalize against
      integer, intent(out) :: info
      !! Information flag.
      logical,                          optional, intent(in)    :: if_chk_orthonormal
      logical                                                   :: chk_X_orthonormality
      !! Check that input Krylov vectors `X` form an orthonormal basis (expensive!)
      complex(dp),                         optional, intent(out)   :: beta(:,:)
      !! Projection coefficients if requested

      ! internals
      complex(dp), dimension(size(X),size(Y)) :: proj_coefficients, wrk

      ! optional input argument
      chk_X_orthonormality = optval(if_chk_orthonormal, .true.) ! default to true!

      info = 0

      proj_coefficients = zero_cdp; wrk = zero_cdp

      ! Orthogonalize Krylov basis Y w.r.t. to Krylov basis X in two passes of GS.
      ! first pass
      call orthogonalize_against_basis(Y, X, info, if_chk_orthonormal=.false., beta=proj_coefficients)
      call check_info(info, 'orthogonalize_against_basis', module=this_module, procedure='DGS_basis_against_basis_cdp, first pass')
      ! second pass
      call orthogonalize_against_basis(Y, X, info, if_chk_orthonormal=.false., beta=wrk)
      call check_info(info, 'orthogonalize_against_basis', module=this_module, procedure='DGS_basis_against_basis_cdp, second pass')
      ! combine passes
      proj_coefficients = proj_coefficients + wrk

      if (present(beta)) then
         ! check size
         call assert_shape(beta, shape(proj_coefficients), 'DGS_basis_against_basis_cdp', 'beta')
         beta = proj_coefficients
      end if

    end subroutine DGS_basis_against_basis_cdp  


    !------------------------------------
    !-----     QR FACTORIZATION     -----
    !------------------------------------

    subroutine qr_no_pivoting_rsp(Q, R, info, verbosity, tol)
        class(abstract_vector_rsp), intent(inout) :: Q(:)
        !! Array of `abstract_vector` to be orthogonalized.
        real(sp), intent(out) :: R(:, :)
        !! Upper triangular matrix \(\mathbf{R}\) resulting from the QR factorization.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        logical                       :: verbose
        real(sp), optional, intent(in) :: tol
        !! Tolerance to determine colinearity.
        real(sp)                       :: tolerance

        ! Internal variables.
        real(sp) :: alpha
        real(sp) :: beta(size(Q))
        integer :: idx, j, kdim, iwrk

        ! Deals with the optional args.
        verbose   = optval(verbosity, .false.)
        tolerance = optval(tol, rtol_sp)

        info = 0 ; R = zero_rsp ; beta = zero_rsp
        do j = 1, size(Q)
            if (j > 1) then
                ! Double Gram-Schmidt orthogonalization
                call double_gram_schmidt_step(Q(j), Q(1:j-1), info, if_chk_orthonormal=.false., beta=R(1:j-1,j))
                call check_info(info, 'double_gram_schmidt_step', module=this_module, procedure='qr_no_pivoting_rsp')
            end if

            ! Normalize column.
            alpha = Q(j)%norm()

            ! Check for breakdown.
            if (abs(alpha) < tolerance) then
                if(verbose) then
                    write(output_unit, *) "INFO: Colinear columns detected."
                    write(output_unit, *) "      (j, alpha) = (", j, ", ", alpha, ")"
                endif
                info = j
                R(j, j) = zero_rsp ; call Q(j)%zero()
            else
                call Q(j)%scal(one_rsp / alpha) ; R(j, j) = alpha
            endif
        enddo

        return
    end subroutine qr_no_pivoting_rsp

    subroutine qr_with_pivoting_rsp(Q, R, perm, info, verbosity, tol)
        class(abstract_vector_rsp), intent(inout) :: Q(:)
        !! Array of `abstract_vector` to be orthogonalized.
        real(sp), intent(out) :: R(:, :)
        !! Upper triangular matrix resulting from the QR factorization.
        integer, intent(out) :: perm(size(Q))
        !! Permutation matrix.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        logical :: verbose
        real(sp), optional, intent(in) :: tol

        !! Tolerance to detect colinearity.
        real(sp) :: tolerance

        ! Internal variables
        real(sp) :: beta
        integer :: idx, i, j, kdim, iwrk
        integer :: idxv(1)
        real(sp)  :: Rii(size(Q))

        info = 0 ; kdim = size(Q)
        R = zero_rsp ; Rii = zero_rsp
        
        ! Deals with the optional arguments.
        verbose = optval(verbosity, .false.)
        tolerance = optval(tol, rtol_sp)

        ! Initialize diagonal entries.
        do i = 1, kdim
            perm(i) = i
            Rii(i) = Q(i)%dot(Q(i))
        enddo

        qr_step: do j = 1, kdim
            idxv = maxloc(abs(Rii)) ; idx = idxv(1)
            if (abs(Rii(idx)) < tolerance) then
                do i = j, kdim
                    call Q(i)%rand(.false.)
                    call orthogonalize_against_basis(Q(i), Q(1:i-1), info, if_chk_orthonormal=.false.)
                    call check_info(info, 'orthogonalize_against_basis', module=this_module, procedure='qr_with_pivoting_rsp')
                    beta = Q(i)%norm(); call Q(i)%scal(one_rsp / beta)
                enddo
                info = j
                exit qr_step
            endif

            call swap_columns(Q, R, Rii, perm, j, idx)

            ! Normalize column.
            beta = Q(j)%norm()

            if (abs(beta) < tolerance) then
               info = j
                R(j, j) = zero_rsp ; call Q(j)%zero()
            else
                R(j, j) = beta ; call Q(j)%scal(one_rsp / beta)
            endif

            ! Orthonormalize all columns against new vector (MGS)
            do i = j+1, kdim
                beta = Q(j)%dot(Q(i))
                call Q(i)%axpby(one_rsp, Q(j), -beta)
                R(j, i) = beta
            enddo

            ! Update Rii
            Rii(j) = zero_rsp
            do i = j+1, kdim
                Rii(i) = Rii(i) - R(j, i)**2
            enddo

        enddo qr_step

        return
    end subroutine qr_with_pivoting_rsp

    subroutine swap_columns_rsp(Q, R, Rii, perm, i, j)
        class(abstract_vector_rsp), intent(inout) :: Q(:)
        !! Vector basis whose i-th and j-th columns need swapping.
        real(sp), intent(inout) :: R(:, :)
        !! Upper triangular matrix resulting from QR.
        real(sp), intent(inout) :: Rii(:)
        !! Diagonal entries of R.
        integer, intent(inout) :: perm(:)
        !! Column ordering.
        integer, intent(in) :: i, j
        !! Index of the columns to be swapped.

        ! Internal variables.
        class(abstract_vector_rsp), allocatable :: Qwrk
        real(sp), allocatable :: Rwrk(:)
        integer :: iwrk, m, n

        ! Sanity checks.
        m = size(Q) ; n = min(i, j) - 1

        ! Allocations.
        allocate(Qwrk, source=Q(1)) ; call Qwrk%zero()
        allocate(Rwrk(1:max(1, n))) ; Rwrk = zero_rsp

        ! Swap columns.
        call Qwrk%axpby(zero_rsp, Q(j), one_rsp)
        call Q(j)%axpby(zero_rsp, Q(i), one_rsp)
        call Q(i)%axpby(zero_rsp, Qwrk, one_rsp)
        
        Rwrk(1) = Rii(j); Rii(j) = Rii(i); Rii(i) = Rwrk(1)
        iwrk = perm(j); perm(j) = perm(i) ; perm(i) = iwrk

        if (n > 0) then
            Rwrk = R(1:n, j) ; R(1:n, j) = R(1:n, i) ; R(1:n, i) = Rwrk
        endif

        return
    end subroutine swap_columns_rsp

    subroutine apply_permutation_matrix_rsp(Q, perm)
        class(abstract_vector_rsp), intent(inout) :: Q(:)
        !! Basis vectors to be permuted.
        integer, intent(in) :: perm(:)
        !! Permutation matrix (vector representation).

        ! Internal variables.
        integer :: i
        class(abstract_vector_rsp), allocatable :: Qwrk(:)

        allocate(Qwrk, source=Q)
        do i = 1, size(perm)
            call Q(i)%axpby(zero_rsp, Qwrk(perm(i)), one_rsp)
        enddo

        return
    end subroutine apply_permutation_matrix_rsp

    subroutine apply_inverse_permutation_matrix_rsp(Q, perm)
        class(abstract_vector_rsp), intent(inout) :: Q(:)
        !! Basis vectors to be (un-) permuted.
        integer, intent(in) :: perm(:)
        !! Permutation matrix (vector representation).

        ! Internal variables.
        integer :: i
        integer :: inv_perm(size(perm))
        class(abstract_vector_rsp), allocatable :: Qwrk(:)

        allocate(Qwrk, source=Q) ; inv_perm = 0

        ! Inverse permutation vector.
        do i = 1, size(perm)
            inv_perm(perm(i)) = i
        enddo

        ! Undo permutation.
        do i = 1, size(perm)
            call Q(i)%axpby(zero_rsp, Qwrk(inv_perm(i)), one_rsp)
        enddo

        return
    end subroutine apply_inverse_permutation_matrix_rsp

    subroutine apply_permutation_matrix_array_rsp(Q, perm)
        real(sp), intent(inout) :: Q(:, :)
        !! Basis vectors to be permuted.
        integer, intent(in) :: perm(:)
        !! Permutation matrix (vector representation).

        ! Internal variables.
        integer :: i
        real(sp), allocatable :: Qwrk(:, :)

        Qwrk = Q
        do i = 1, size(perm)
            Q(:, i) = Qwrk(:, perm(i))
        enddo

        return
    end subroutine apply_permutation_matrix_array_rsp

    subroutine apply_inverse_permutation_matrix_array_rsp(Q, perm)
        real(sp), intent(inout) :: Q(:, :)
        !! Basis vectors to be (un-) permuted.
        integer, intent(in) :: perm(:)
        !! Permutation matrix (vector representation).

        ! Internal variables.
        integer :: i
        integer :: inv_perm(size(perm))
        real(sp), allocatable :: Qwrk(:, :)

        Qwrk = Q ; inv_perm = 0

        ! Inverse permutation vector.
        do i = 1, size(perm)
            inv_perm(perm(i)) = i
        enddo

        ! Undo permutation.
        do i = 1, size(perm)
            Q(:, i) = Qwrk(:, inv_perm(i))
        enddo

        return
    end subroutine apply_inverse_permutation_matrix_array_rsp


    subroutine qr_no_pivoting_rdp(Q, R, info, verbosity, tol)
        class(abstract_vector_rdp), intent(inout) :: Q(:)
        !! Array of `abstract_vector` to be orthogonalized.
        real(dp), intent(out) :: R(:, :)
        !! Upper triangular matrix \(\mathbf{R}\) resulting from the QR factorization.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        logical                       :: verbose
        real(dp), optional, intent(in) :: tol
        !! Tolerance to determine colinearity.
        real(dp)                       :: tolerance

        ! Internal variables.
        real(dp) :: alpha
        real(dp) :: beta(size(Q))
        integer :: idx, j, kdim, iwrk

        ! Deals with the optional args.
        verbose   = optval(verbosity, .false.)
        tolerance = optval(tol, rtol_dp)

        info = 0 ; R = zero_rdp ; beta = zero_rdp
        do j = 1, size(Q)
            if (j > 1) then
                ! Double Gram-Schmidt orthogonalization
                call double_gram_schmidt_step(Q(j), Q(1:j-1), info, if_chk_orthonormal=.false., beta=R(1:j-1,j))
                call check_info(info, 'double_gram_schmidt_step', module=this_module, procedure='qr_no_pivoting_rdp')
            end if

            ! Normalize column.
            alpha = Q(j)%norm()

            ! Check for breakdown.
            if (abs(alpha) < tolerance) then
                if(verbose) then
                    write(output_unit, *) "INFO: Colinear columns detected."
                    write(output_unit, *) "      (j, alpha) = (", j, ", ", alpha, ")"
                endif
                info = j
                R(j, j) = zero_rdp ; call Q(j)%zero()
            else
                call Q(j)%scal(one_rdp / alpha) ; R(j, j) = alpha
            endif
        enddo

        return
    end subroutine qr_no_pivoting_rdp

    subroutine qr_with_pivoting_rdp(Q, R, perm, info, verbosity, tol)
        class(abstract_vector_rdp), intent(inout) :: Q(:)
        !! Array of `abstract_vector` to be orthogonalized.
        real(dp), intent(out) :: R(:, :)
        !! Upper triangular matrix resulting from the QR factorization.
        integer, intent(out) :: perm(size(Q))
        !! Permutation matrix.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        logical :: verbose
        real(dp), optional, intent(in) :: tol

        !! Tolerance to detect colinearity.
        real(dp) :: tolerance

        ! Internal variables
        real(dp) :: beta
        integer :: idx, i, j, kdim, iwrk
        integer :: idxv(1)
        real(dp)  :: Rii(size(Q))

        info = 0 ; kdim = size(Q)
        R = zero_rdp ; Rii = zero_rdp
        
        ! Deals with the optional arguments.
        verbose = optval(verbosity, .false.)
        tolerance = optval(tol, rtol_dp)

        ! Initialize diagonal entries.
        do i = 1, kdim
            perm(i) = i
            Rii(i) = Q(i)%dot(Q(i))
        enddo

        qr_step: do j = 1, kdim
            idxv = maxloc(abs(Rii)) ; idx = idxv(1)
            if (abs(Rii(idx)) < tolerance) then
                do i = j, kdim
                    call Q(i)%rand(.false.)
                    call orthogonalize_against_basis(Q(i), Q(1:i-1), info, if_chk_orthonormal=.false.)
                    call check_info(info, 'orthogonalize_against_basis', module=this_module, procedure='qr_with_pivoting_rdp')
                    beta = Q(i)%norm(); call Q(i)%scal(one_rdp / beta)
                enddo
                info = j
                exit qr_step
            endif

            call swap_columns(Q, R, Rii, perm, j, idx)

            ! Normalize column.
            beta = Q(j)%norm()

            if (abs(beta) < tolerance) then
               info = j
                R(j, j) = zero_rdp ; call Q(j)%zero()
            else
                R(j, j) = beta ; call Q(j)%scal(one_rdp / beta)
            endif

            ! Orthonormalize all columns against new vector (MGS)
            do i = j+1, kdim
                beta = Q(j)%dot(Q(i))
                call Q(i)%axpby(one_rdp, Q(j), -beta)
                R(j, i) = beta
            enddo

            ! Update Rii
            Rii(j) = zero_rdp
            do i = j+1, kdim
                Rii(i) = Rii(i) - R(j, i)**2
            enddo

        enddo qr_step

        return
    end subroutine qr_with_pivoting_rdp

    subroutine swap_columns_rdp(Q, R, Rii, perm, i, j)
        class(abstract_vector_rdp), intent(inout) :: Q(:)
        !! Vector basis whose i-th and j-th columns need swapping.
        real(dp), intent(inout) :: R(:, :)
        !! Upper triangular matrix resulting from QR.
        real(dp), intent(inout) :: Rii(:)
        !! Diagonal entries of R.
        integer, intent(inout) :: perm(:)
        !! Column ordering.
        integer, intent(in) :: i, j
        !! Index of the columns to be swapped.

        ! Internal variables.
        class(abstract_vector_rdp), allocatable :: Qwrk
        real(dp), allocatable :: Rwrk(:)
        integer :: iwrk, m, n

        ! Sanity checks.
        m = size(Q) ; n = min(i, j) - 1

        ! Allocations.
        allocate(Qwrk, source=Q(1)) ; call Qwrk%zero()
        allocate(Rwrk(1:max(1, n))) ; Rwrk = zero_rdp

        ! Swap columns.
        call Qwrk%axpby(zero_rdp, Q(j), one_rdp)
        call Q(j)%axpby(zero_rdp, Q(i), one_rdp)
        call Q(i)%axpby(zero_rdp, Qwrk, one_rdp)
        
        Rwrk(1) = Rii(j); Rii(j) = Rii(i); Rii(i) = Rwrk(1)
        iwrk = perm(j); perm(j) = perm(i) ; perm(i) = iwrk

        if (n > 0) then
            Rwrk = R(1:n, j) ; R(1:n, j) = R(1:n, i) ; R(1:n, i) = Rwrk
        endif

        return
    end subroutine swap_columns_rdp

    subroutine apply_permutation_matrix_rdp(Q, perm)
        class(abstract_vector_rdp), intent(inout) :: Q(:)
        !! Basis vectors to be permuted.
        integer, intent(in) :: perm(:)
        !! Permutation matrix (vector representation).

        ! Internal variables.
        integer :: i
        class(abstract_vector_rdp), allocatable :: Qwrk(:)

        allocate(Qwrk, source=Q)
        do i = 1, size(perm)
            call Q(i)%axpby(zero_rdp, Qwrk(perm(i)), one_rdp)
        enddo

        return
    end subroutine apply_permutation_matrix_rdp

    subroutine apply_inverse_permutation_matrix_rdp(Q, perm)
        class(abstract_vector_rdp), intent(inout) :: Q(:)
        !! Basis vectors to be (un-) permuted.
        integer, intent(in) :: perm(:)
        !! Permutation matrix (vector representation).

        ! Internal variables.
        integer :: i
        integer :: inv_perm(size(perm))
        class(abstract_vector_rdp), allocatable :: Qwrk(:)

        allocate(Qwrk, source=Q) ; inv_perm = 0

        ! Inverse permutation vector.
        do i = 1, size(perm)
            inv_perm(perm(i)) = i
        enddo

        ! Undo permutation.
        do i = 1, size(perm)
            call Q(i)%axpby(zero_rdp, Qwrk(inv_perm(i)), one_rdp)
        enddo

        return
    end subroutine apply_inverse_permutation_matrix_rdp

    subroutine apply_permutation_matrix_array_rdp(Q, perm)
        real(dp), intent(inout) :: Q(:, :)
        !! Basis vectors to be permuted.
        integer, intent(in) :: perm(:)
        !! Permutation matrix (vector representation).

        ! Internal variables.
        integer :: i
        real(dp), allocatable :: Qwrk(:, :)

        Qwrk = Q
        do i = 1, size(perm)
            Q(:, i) = Qwrk(:, perm(i))
        enddo

        return
    end subroutine apply_permutation_matrix_array_rdp

    subroutine apply_inverse_permutation_matrix_array_rdp(Q, perm)
        real(dp), intent(inout) :: Q(:, :)
        !! Basis vectors to be (un-) permuted.
        integer, intent(in) :: perm(:)
        !! Permutation matrix (vector representation).

        ! Internal variables.
        integer :: i
        integer :: inv_perm(size(perm))
        real(dp), allocatable :: Qwrk(:, :)

        Qwrk = Q ; inv_perm = 0

        ! Inverse permutation vector.
        do i = 1, size(perm)
            inv_perm(perm(i)) = i
        enddo

        ! Undo permutation.
        do i = 1, size(perm)
            Q(:, i) = Qwrk(:, inv_perm(i))
        enddo

        return
    end subroutine apply_inverse_permutation_matrix_array_rdp


    subroutine qr_no_pivoting_csp(Q, R, info, verbosity, tol)
        class(abstract_vector_csp), intent(inout) :: Q(:)
        !! Array of `abstract_vector` to be orthogonalized.
        complex(sp), intent(out) :: R(:, :)
        !! Upper triangular matrix \(\mathbf{R}\) resulting from the QR factorization.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        logical                       :: verbose
        real(sp), optional, intent(in) :: tol
        !! Tolerance to determine colinearity.
        real(sp)                       :: tolerance

        ! Internal variables.
        complex(sp) :: alpha
        complex(sp) :: beta(size(Q))
        integer :: idx, j, kdim, iwrk

        ! Deals with the optional args.
        verbose   = optval(verbosity, .false.)
        tolerance = optval(tol, rtol_sp)

        info = 0 ; R = zero_rsp ; beta = zero_rsp
        do j = 1, size(Q)
            if (j > 1) then
                ! Double Gram-Schmidt orthogonalization
                call double_gram_schmidt_step(Q(j), Q(1:j-1), info, if_chk_orthonormal=.false., beta=R(1:j-1,j))
                call check_info(info, 'double_gram_schmidt_step', module=this_module, procedure='qr_no_pivoting_csp')
            end if

            ! Normalize column.
            alpha = Q(j)%norm()

            ! Check for breakdown.
            if (abs(alpha) < tolerance) then
                if(verbose) then
                    write(output_unit, *) "INFO: Colinear columns detected."
                    write(output_unit, *) "      (j, alpha) = (", j, ", ", alpha, ")"
                endif
                info = j
                R(j, j) = zero_rsp ; call Q(j)%zero()
            else
                call Q(j)%scal(one_rsp / alpha) ; R(j, j) = alpha
            endif
        enddo

        return
    end subroutine qr_no_pivoting_csp

    subroutine qr_with_pivoting_csp(Q, R, perm, info, verbosity, tol)
        class(abstract_vector_csp), intent(inout) :: Q(:)
        !! Array of `abstract_vector` to be orthogonalized.
        complex(sp), intent(out) :: R(:, :)
        !! Upper triangular matrix resulting from the QR factorization.
        integer, intent(out) :: perm(size(Q))
        !! Permutation matrix.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        logical :: verbose
        real(sp), optional, intent(in) :: tol

        !! Tolerance to detect colinearity.
        real(sp) :: tolerance

        ! Internal variables
        complex(sp) :: beta
        integer :: idx, i, j, kdim, iwrk
        integer :: idxv(1)
        complex(sp)  :: Rii(size(Q))

        info = 0 ; kdim = size(Q)
        R = zero_rsp ; Rii = zero_rsp
        
        ! Deals with the optional arguments.
        verbose = optval(verbosity, .false.)
        tolerance = optval(tol, rtol_sp)

        ! Initialize diagonal entries.
        do i = 1, kdim
            perm(i) = i
            Rii(i) = Q(i)%dot(Q(i))
        enddo

        qr_step: do j = 1, kdim
            idxv = maxloc(abs(Rii)) ; idx = idxv(1)
            if (abs(Rii(idx)) < tolerance) then
                do i = j, kdim
                    call Q(i)%rand(.false.)
                    call orthogonalize_against_basis(Q(i), Q(1:i-1), info, if_chk_orthonormal=.false.)
                    call check_info(info, 'orthogonalize_against_basis', module=this_module, procedure='qr_with_pivoting_csp')
                    beta = Q(i)%norm(); call Q(i)%scal(one_csp / beta)
                enddo
                info = j
                exit qr_step
            endif

            call swap_columns(Q, R, Rii, perm, j, idx)

            ! Normalize column.
            beta = Q(j)%norm()

            if (abs(beta) < tolerance) then
               info = j
                R(j, j) = zero_rsp ; call Q(j)%zero()
            else
                R(j, j) = beta ; call Q(j)%scal(one_rsp / beta)
            endif

            ! Orthonormalize all columns against new vector (MGS)
            do i = j+1, kdim
                beta = Q(j)%dot(Q(i))
                call Q(i)%axpby(one_csp, Q(j), -beta)
                R(j, i) = beta
            enddo

            ! Update Rii
            Rii(j) = zero_rsp
            do i = j+1, kdim
                Rii(i) = Rii(i) - R(j, i)**2
            enddo

        enddo qr_step

        return
    end subroutine qr_with_pivoting_csp

    subroutine swap_columns_csp(Q, R, Rii, perm, i, j)
        class(abstract_vector_csp), intent(inout) :: Q(:)
        !! Vector basis whose i-th and j-th columns need swapping.
        complex(sp), intent(inout) :: R(:, :)
        !! Upper triangular matrix resulting from QR.
        complex(sp), intent(inout) :: Rii(:)
        !! Diagonal entries of R.
        integer, intent(inout) :: perm(:)
        !! Column ordering.
        integer, intent(in) :: i, j
        !! Index of the columns to be swapped.

        ! Internal variables.
        class(abstract_vector_csp), allocatable :: Qwrk
        complex(sp), allocatable :: Rwrk(:)
        integer :: iwrk, m, n

        ! Sanity checks.
        m = size(Q) ; n = min(i, j) - 1

        ! Allocations.
        allocate(Qwrk, source=Q(1)) ; call Qwrk%zero()
        allocate(Rwrk(1:max(1, n))) ; Rwrk = zero_rsp

        ! Swap columns.
        call Qwrk%axpby(zero_csp, Q(j), one_csp)
        call Q(j)%axpby(zero_csp, Q(i), one_csp)
        call Q(i)%axpby(zero_csp, Qwrk, one_csp)
        
        Rwrk(1) = Rii(j); Rii(j) = Rii(i); Rii(i) = Rwrk(1)
        iwrk = perm(j); perm(j) = perm(i) ; perm(i) = iwrk

        if (n > 0) then
            Rwrk = R(1:n, j) ; R(1:n, j) = R(1:n, i) ; R(1:n, i) = Rwrk
        endif

        return
    end subroutine swap_columns_csp

    subroutine apply_permutation_matrix_csp(Q, perm)
        class(abstract_vector_csp), intent(inout) :: Q(:)
        !! Basis vectors to be permuted.
        integer, intent(in) :: perm(:)
        !! Permutation matrix (vector representation).

        ! Internal variables.
        integer :: i
        class(abstract_vector_csp), allocatable :: Qwrk(:)

        allocate(Qwrk, source=Q)
        do i = 1, size(perm)
            call Q(i)%axpby(zero_csp, Qwrk(perm(i)), one_csp)
        enddo

        return
    end subroutine apply_permutation_matrix_csp

    subroutine apply_inverse_permutation_matrix_csp(Q, perm)
        class(abstract_vector_csp), intent(inout) :: Q(:)
        !! Basis vectors to be (un-) permuted.
        integer, intent(in) :: perm(:)
        !! Permutation matrix (vector representation).

        ! Internal variables.
        integer :: i
        integer :: inv_perm(size(perm))
        class(abstract_vector_csp), allocatable :: Qwrk(:)

        allocate(Qwrk, source=Q) ; inv_perm = 0

        ! Inverse permutation vector.
        do i = 1, size(perm)
            inv_perm(perm(i)) = i
        enddo

        ! Undo permutation.
        do i = 1, size(perm)
            call Q(i)%axpby(zero_csp, Qwrk(inv_perm(i)), one_csp)
        enddo

        return
    end subroutine apply_inverse_permutation_matrix_csp

    subroutine apply_permutation_matrix_array_csp(Q, perm)
        complex(sp), intent(inout) :: Q(:, :)
        !! Basis vectors to be permuted.
        integer, intent(in) :: perm(:)
        !! Permutation matrix (vector representation).

        ! Internal variables.
        integer :: i
        complex(sp), allocatable :: Qwrk(:, :)

        Qwrk = Q
        do i = 1, size(perm)
            Q(:, i) = Qwrk(:, perm(i))
        enddo

        return
    end subroutine apply_permutation_matrix_array_csp

    subroutine apply_inverse_permutation_matrix_array_csp(Q, perm)
        complex(sp), intent(inout) :: Q(:, :)
        !! Basis vectors to be (un-) permuted.
        integer, intent(in) :: perm(:)
        !! Permutation matrix (vector representation).

        ! Internal variables.
        integer :: i
        integer :: inv_perm(size(perm))
        complex(sp), allocatable :: Qwrk(:, :)

        Qwrk = Q ; inv_perm = 0

        ! Inverse permutation vector.
        do i = 1, size(perm)
            inv_perm(perm(i)) = i
        enddo

        ! Undo permutation.
        do i = 1, size(perm)
            Q(:, i) = Qwrk(:, inv_perm(i))
        enddo

        return
    end subroutine apply_inverse_permutation_matrix_array_csp


    subroutine qr_no_pivoting_cdp(Q, R, info, verbosity, tol)
        class(abstract_vector_cdp), intent(inout) :: Q(:)
        !! Array of `abstract_vector` to be orthogonalized.
        complex(dp), intent(out) :: R(:, :)
        !! Upper triangular matrix \(\mathbf{R}\) resulting from the QR factorization.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        logical                       :: verbose
        real(dp), optional, intent(in) :: tol
        !! Tolerance to determine colinearity.
        real(dp)                       :: tolerance

        ! Internal variables.
        complex(dp) :: alpha
        complex(dp) :: beta(size(Q))
        integer :: idx, j, kdim, iwrk

        ! Deals with the optional args.
        verbose   = optval(verbosity, .false.)
        tolerance = optval(tol, rtol_dp)

        info = 0 ; R = zero_rdp ; beta = zero_rdp
        do j = 1, size(Q)
            if (j > 1) then
                ! Double Gram-Schmidt orthogonalization
                call double_gram_schmidt_step(Q(j), Q(1:j-1), info, if_chk_orthonormal=.false., beta=R(1:j-1,j))
                call check_info(info, 'double_gram_schmidt_step', module=this_module, procedure='qr_no_pivoting_cdp')
            end if

            ! Normalize column.
            alpha = Q(j)%norm()

            ! Check for breakdown.
            if (abs(alpha) < tolerance) then
                if(verbose) then
                    write(output_unit, *) "INFO: Colinear columns detected."
                    write(output_unit, *) "      (j, alpha) = (", j, ", ", alpha, ")"
                endif
                info = j
                R(j, j) = zero_rdp ; call Q(j)%zero()
            else
                call Q(j)%scal(one_rdp / alpha) ; R(j, j) = alpha
            endif
        enddo

        return
    end subroutine qr_no_pivoting_cdp

    subroutine qr_with_pivoting_cdp(Q, R, perm, info, verbosity, tol)
        class(abstract_vector_cdp), intent(inout) :: Q(:)
        !! Array of `abstract_vector` to be orthogonalized.
        complex(dp), intent(out) :: R(:, :)
        !! Upper triangular matrix resulting from the QR factorization.
        integer, intent(out) :: perm(size(Q))
        !! Permutation matrix.
        integer, intent(out) :: info
        !! Information flag.
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        logical :: verbose
        real(dp), optional, intent(in) :: tol

        !! Tolerance to detect colinearity.
        real(dp) :: tolerance

        ! Internal variables
        complex(dp) :: beta
        integer :: idx, i, j, kdim, iwrk
        integer :: idxv(1)
        complex(dp)  :: Rii(size(Q))

        info = 0 ; kdim = size(Q)
        R = zero_rdp ; Rii = zero_rdp
        
        ! Deals with the optional arguments.
        verbose = optval(verbosity, .false.)
        tolerance = optval(tol, rtol_dp)

        ! Initialize diagonal entries.
        do i = 1, kdim
            perm(i) = i
            Rii(i) = Q(i)%dot(Q(i))
        enddo

        qr_step: do j = 1, kdim
            idxv = maxloc(abs(Rii)) ; idx = idxv(1)
            if (abs(Rii(idx)) < tolerance) then
                do i = j, kdim
                    call Q(i)%rand(.false.)
                    call orthogonalize_against_basis(Q(i), Q(1:i-1), info, if_chk_orthonormal=.false.)
                    call check_info(info, 'orthogonalize_against_basis', module=this_module, procedure='qr_with_pivoting_cdp')
                    beta = Q(i)%norm(); call Q(i)%scal(one_cdp / beta)
                enddo
                info = j
                exit qr_step
            endif

            call swap_columns(Q, R, Rii, perm, j, idx)

            ! Normalize column.
            beta = Q(j)%norm()

            if (abs(beta) < tolerance) then
               info = j
                R(j, j) = zero_rdp ; call Q(j)%zero()
            else
                R(j, j) = beta ; call Q(j)%scal(one_rdp / beta)
            endif

            ! Orthonormalize all columns against new vector (MGS)
            do i = j+1, kdim
                beta = Q(j)%dot(Q(i))
                call Q(i)%axpby(one_cdp, Q(j), -beta)
                R(j, i) = beta
            enddo

            ! Update Rii
            Rii(j) = zero_rdp
            do i = j+1, kdim
                Rii(i) = Rii(i) - R(j, i)**2
            enddo

        enddo qr_step

        return
    end subroutine qr_with_pivoting_cdp

    subroutine swap_columns_cdp(Q, R, Rii, perm, i, j)
        class(abstract_vector_cdp), intent(inout) :: Q(:)
        !! Vector basis whose i-th and j-th columns need swapping.
        complex(dp), intent(inout) :: R(:, :)
        !! Upper triangular matrix resulting from QR.
        complex(dp), intent(inout) :: Rii(:)
        !! Diagonal entries of R.
        integer, intent(inout) :: perm(:)
        !! Column ordering.
        integer, intent(in) :: i, j
        !! Index of the columns to be swapped.

        ! Internal variables.
        class(abstract_vector_cdp), allocatable :: Qwrk
        complex(dp), allocatable :: Rwrk(:)
        integer :: iwrk, m, n

        ! Sanity checks.
        m = size(Q) ; n = min(i, j) - 1

        ! Allocations.
        allocate(Qwrk, source=Q(1)) ; call Qwrk%zero()
        allocate(Rwrk(1:max(1, n))) ; Rwrk = zero_rdp

        ! Swap columns.
        call Qwrk%axpby(zero_cdp, Q(j), one_cdp)
        call Q(j)%axpby(zero_cdp, Q(i), one_cdp)
        call Q(i)%axpby(zero_cdp, Qwrk, one_cdp)
        
        Rwrk(1) = Rii(j); Rii(j) = Rii(i); Rii(i) = Rwrk(1)
        iwrk = perm(j); perm(j) = perm(i) ; perm(i) = iwrk

        if (n > 0) then
            Rwrk = R(1:n, j) ; R(1:n, j) = R(1:n, i) ; R(1:n, i) = Rwrk
        endif

        return
    end subroutine swap_columns_cdp

    subroutine apply_permutation_matrix_cdp(Q, perm)
        class(abstract_vector_cdp), intent(inout) :: Q(:)
        !! Basis vectors to be permuted.
        integer, intent(in) :: perm(:)
        !! Permutation matrix (vector representation).

        ! Internal variables.
        integer :: i
        class(abstract_vector_cdp), allocatable :: Qwrk(:)

        allocate(Qwrk, source=Q)
        do i = 1, size(perm)
            call Q(i)%axpby(zero_cdp, Qwrk(perm(i)), one_cdp)
        enddo

        return
    end subroutine apply_permutation_matrix_cdp

    subroutine apply_inverse_permutation_matrix_cdp(Q, perm)
        class(abstract_vector_cdp), intent(inout) :: Q(:)
        !! Basis vectors to be (un-) permuted.
        integer, intent(in) :: perm(:)
        !! Permutation matrix (vector representation).

        ! Internal variables.
        integer :: i
        integer :: inv_perm(size(perm))
        class(abstract_vector_cdp), allocatable :: Qwrk(:)

        allocate(Qwrk, source=Q) ; inv_perm = 0

        ! Inverse permutation vector.
        do i = 1, size(perm)
            inv_perm(perm(i)) = i
        enddo

        ! Undo permutation.
        do i = 1, size(perm)
            call Q(i)%axpby(zero_cdp, Qwrk(inv_perm(i)), one_cdp)
        enddo

        return
    end subroutine apply_inverse_permutation_matrix_cdp

    subroutine apply_permutation_matrix_array_cdp(Q, perm)
        complex(dp), intent(inout) :: Q(:, :)
        !! Basis vectors to be permuted.
        integer, intent(in) :: perm(:)
        !! Permutation matrix (vector representation).

        ! Internal variables.
        integer :: i
        complex(dp), allocatable :: Qwrk(:, :)

        Qwrk = Q
        do i = 1, size(perm)
            Q(:, i) = Qwrk(:, perm(i))
        enddo

        return
    end subroutine apply_permutation_matrix_array_cdp

    subroutine apply_inverse_permutation_matrix_array_cdp(Q, perm)
        complex(dp), intent(inout) :: Q(:, :)
        !! Basis vectors to be (un-) permuted.
        integer, intent(in) :: perm(:)
        !! Permutation matrix (vector representation).

        ! Internal variables.
        integer :: i
        integer :: inv_perm(size(perm))
        complex(dp), allocatable :: Qwrk(:, :)

        Qwrk = Q ; inv_perm = 0

        ! Inverse permutation vector.
        do i = 1, size(perm)
            inv_perm(perm(i)) = i
        enddo

        ! Undo permutation.
        do i = 1, size(perm)
            Q(:, i) = Qwrk(:, inv_perm(i))
        enddo

        return
    end subroutine apply_inverse_permutation_matrix_array_cdp



    !-----------------------------------------
    !-----     ARNOLDI FACTORIZATION     -----
    !-----------------------------------------

    subroutine arnoldi_rsp(A, X, H, info, kstart, kend, verbosity, tol, transpose, blksize)
        class(abstract_linop_rsp), intent(in) :: A
        !! Linear operator to be factorized.
        class(abstract_vector_rsp), intent(inout) :: X(:)
        !! Orthogonal basis for the generated Krylov subspace.
        real(sp), intent(inout) :: H(:, :)
        !! Upper Hessenberg matrix.
        integer, intent(out) :: info
        !! Information flag.
        integer, optional, intent(in) :: kstart
        !! Starting index for the Arnoldi factorization (default 1).
        integer, optional, intent(in) :: kend
        !! Final index for the Arnoldi factorization (default `size(X)-1`)
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        logical, optional, intent(in) :: transpose
        !! Whether \( \mathbf{A} \) is being transposed or not (default `.false.`)
        real(sp), optional, intent(in) :: tol
        !! Tolerance to determine whether an invariant subspace has been computed or not.
        integer, optional, intent(in) :: blksize
        !! Block size for block Arnoldi (default 1).

        ! Internal variables.
        integer :: k_start, k_end, p
        logical :: verbose, trans
        real(sp) :: tolerance
        real(sp) :: beta
        real(sp), allocatable :: res(:)
        integer, allocatable :: perm(:)
        integer :: k, i, kdim, kpm, kp, kpp

        ! Deals with optional non-unity blksize and allocations.
        p = optval(blksize, 1) ; allocate(res(1:p)) ; res = zero_rsp
        allocate(perm(1:size(H, 2))) ; perm = 0 ; info = 0

        ! Check dimensions.
        kdim = (size(X) - p) / p

        ! Deal with the other optional args.
        k_start = optval(kstart, 1) ; k_end = optval(kend, kdim)
        verbose   = optval(verbosity, .false.)
        tolerance = optval (tol, atol_sp)
        trans     = optval(transpose, .false.)

        ! Arnoldi factorization.
        blk_arnoldi: do k = k_start, k_end
            ! Counters
            kpm = (k - 1) * p ; kp = kpm + p ; kpp = kp + p

            ! Matrix-vector product.
            if (trans) then
                do i = 1, p
                    call A%rmatvec(X(kpm+i), X(kp+i))
                enddo
            else
                do i = 1, p
                    call A%matvec(X(kpm+i), X(kp+i))
                enddo
            endif

            ! Update Hessenberg matrix via batch double Gram-Schmidt step.
            call double_gram_schmidt_step(X(kp+1:kpp), X(1:kp), info, if_chk_orthonormal=.false., beta=H(1:kp, kpm+1:kp))
            call check_info(info, 'double_gram_schmidt_step', module=this_module, procedure='arnoldi_rsp')

            ! Orthogonalize current blk vectors.
            call qr(X(kp+1:kpp), H(kp+1:kpp, kpm+1:kp), info)
            call check_info(info, 'qr', module=this_module, procedure='arnoldi_rsp')

            ! Extract residual norm (smallest diagonal element of H matrix).
            res = zero_rsp
            do i = 1, p
                res(i) = H(kp+i, kpm+i)
            enddo
            beta = minval(abs(res))

            ! Exit Arnoldi loop if needed.
            if (beta < tolerance) then
                ! Dimension of the computed invariant subspace.
                info = kp
                ! Exit the Arnoldi iteration.
                exit blk_arnoldi
            endif

        enddo blk_arnoldi

        return
    end subroutine arnoldi_rsp

    subroutine arnoldi_rdp(A, X, H, info, kstart, kend, verbosity, tol, transpose, blksize)
        class(abstract_linop_rdp), intent(in) :: A
        !! Linear operator to be factorized.
        class(abstract_vector_rdp), intent(inout) :: X(:)
        !! Orthogonal basis for the generated Krylov subspace.
        real(dp), intent(inout) :: H(:, :)
        !! Upper Hessenberg matrix.
        integer, intent(out) :: info
        !! Information flag.
        integer, optional, intent(in) :: kstart
        !! Starting index for the Arnoldi factorization (default 1).
        integer, optional, intent(in) :: kend
        !! Final index for the Arnoldi factorization (default `size(X)-1`)
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        logical, optional, intent(in) :: transpose
        !! Whether \( \mathbf{A} \) is being transposed or not (default `.false.`)
        real(dp), optional, intent(in) :: tol
        !! Tolerance to determine whether an invariant subspace has been computed or not.
        integer, optional, intent(in) :: blksize
        !! Block size for block Arnoldi (default 1).

        ! Internal variables.
        integer :: k_start, k_end, p
        logical :: verbose, trans
        real(dp) :: tolerance
        real(dp) :: beta
        real(dp), allocatable :: res(:)
        integer, allocatable :: perm(:)
        integer :: k, i, kdim, kpm, kp, kpp

        ! Deals with optional non-unity blksize and allocations.
        p = optval(blksize, 1) ; allocate(res(1:p)) ; res = zero_rdp
        allocate(perm(1:size(H, 2))) ; perm = 0 ; info = 0

        ! Check dimensions.
        kdim = (size(X) - p) / p

        ! Deal with the other optional args.
        k_start = optval(kstart, 1) ; k_end = optval(kend, kdim)
        verbose   = optval(verbosity, .false.)
        tolerance = optval (tol, atol_dp)
        trans     = optval(transpose, .false.)

        ! Arnoldi factorization.
        blk_arnoldi: do k = k_start, k_end
            ! Counters
            kpm = (k - 1) * p ; kp = kpm + p ; kpp = kp + p

            ! Matrix-vector product.
            if (trans) then
                do i = 1, p
                    call A%rmatvec(X(kpm+i), X(kp+i))
                enddo
            else
                do i = 1, p
                    call A%matvec(X(kpm+i), X(kp+i))
                enddo
            endif

            ! Update Hessenberg matrix via batch double Gram-Schmidt step.
            call double_gram_schmidt_step(X(kp+1:kpp), X(1:kp), info, if_chk_orthonormal=.false., beta=H(1:kp, kpm+1:kp))
            call check_info(info, 'double_gram_schmidt_step', module=this_module, procedure='arnoldi_rdp')

            ! Orthogonalize current blk vectors.
            call qr(X(kp+1:kpp), H(kp+1:kpp, kpm+1:kp), info)
            call check_info(info, 'qr', module=this_module, procedure='arnoldi_rdp')

            ! Extract residual norm (smallest diagonal element of H matrix).
            res = zero_rdp
            do i = 1, p
                res(i) = H(kp+i, kpm+i)
            enddo
            beta = minval(abs(res))

            ! Exit Arnoldi loop if needed.
            if (beta < tolerance) then
                ! Dimension of the computed invariant subspace.
                info = kp
                ! Exit the Arnoldi iteration.
                exit blk_arnoldi
            endif

        enddo blk_arnoldi

        return
    end subroutine arnoldi_rdp

    subroutine arnoldi_csp(A, X, H, info, kstart, kend, verbosity, tol, transpose, blksize)
        class(abstract_linop_csp), intent(in) :: A
        !! Linear operator to be factorized.
        class(abstract_vector_csp), intent(inout) :: X(:)
        !! Orthogonal basis for the generated Krylov subspace.
        complex(sp), intent(inout) :: H(:, :)
        !! Upper Hessenberg matrix.
        integer, intent(out) :: info
        !! Information flag.
        integer, optional, intent(in) :: kstart
        !! Starting index for the Arnoldi factorization (default 1).
        integer, optional, intent(in) :: kend
        !! Final index for the Arnoldi factorization (default `size(X)-1`)
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        logical, optional, intent(in) :: transpose
        !! Whether \( \mathbf{A} \) is being transposed or not (default `.false.`)
        real(sp), optional, intent(in) :: tol
        !! Tolerance to determine whether an invariant subspace has been computed or not.
        integer, optional, intent(in) :: blksize
        !! Block size for block Arnoldi (default 1).

        ! Internal variables.
        integer :: k_start, k_end, p
        logical :: verbose, trans
        real(sp) :: tolerance
        real(sp) :: beta
        complex(sp), allocatable :: res(:)
        integer, allocatable :: perm(:)
        integer :: k, i, kdim, kpm, kp, kpp

        ! Deals with optional non-unity blksize and allocations.
        p = optval(blksize, 1) ; allocate(res(1:p)) ; res = zero_rsp
        allocate(perm(1:size(H, 2))) ; perm = 0 ; info = 0

        ! Check dimensions.
        kdim = (size(X) - p) / p

        ! Deal with the other optional args.
        k_start = optval(kstart, 1) ; k_end = optval(kend, kdim)
        verbose   = optval(verbosity, .false.)
        tolerance = optval (tol, atol_sp)
        trans     = optval(transpose, .false.)

        ! Arnoldi factorization.
        blk_arnoldi: do k = k_start, k_end
            ! Counters
            kpm = (k - 1) * p ; kp = kpm + p ; kpp = kp + p

            ! Matrix-vector product.
            if (trans) then
                do i = 1, p
                    call A%rmatvec(X(kpm+i), X(kp+i))
                enddo
            else
                do i = 1, p
                    call A%matvec(X(kpm+i), X(kp+i))
                enddo
            endif

            ! Update Hessenberg matrix via batch double Gram-Schmidt step.
            call double_gram_schmidt_step(X(kp+1:kpp), X(1:kp), info, if_chk_orthonormal=.false., beta=H(1:kp, kpm+1:kp))
            call check_info(info, 'double_gram_schmidt_step', module=this_module, procedure='arnoldi_csp')

            ! Orthogonalize current blk vectors.
            call qr(X(kp+1:kpp), H(kp+1:kpp, kpm+1:kp), info)
            call check_info(info, 'qr', module=this_module, procedure='arnoldi_csp')

            ! Extract residual norm (smallest diagonal element of H matrix).
            res = zero_rsp
            do i = 1, p
                res(i) = H(kp+i, kpm+i)
            enddo
            beta = minval(abs(res))

            ! Exit Arnoldi loop if needed.
            if (beta < tolerance) then
                ! Dimension of the computed invariant subspace.
                info = kp
                ! Exit the Arnoldi iteration.
                exit blk_arnoldi
            endif

        enddo blk_arnoldi

        return
    end subroutine arnoldi_csp

    subroutine arnoldi_cdp(A, X, H, info, kstart, kend, verbosity, tol, transpose, blksize)
        class(abstract_linop_cdp), intent(in) :: A
        !! Linear operator to be factorized.
        class(abstract_vector_cdp), intent(inout) :: X(:)
        !! Orthogonal basis for the generated Krylov subspace.
        complex(dp), intent(inout) :: H(:, :)
        !! Upper Hessenberg matrix.
        integer, intent(out) :: info
        !! Information flag.
        integer, optional, intent(in) :: kstart
        !! Starting index for the Arnoldi factorization (default 1).
        integer, optional, intent(in) :: kend
        !! Final index for the Arnoldi factorization (default `size(X)-1`)
        logical, optional, intent(in) :: verbosity
        !! Verbosity control.
        logical, optional, intent(in) :: transpose
        !! Whether \( \mathbf{A} \) is being transposed or not (default `.false.`)
        real(dp), optional, intent(in) :: tol
        !! Tolerance to determine whether an invariant subspace has been computed or not.
        integer, optional, intent(in) :: blksize
        !! Block size for block Arnoldi (default 1).

        ! Internal variables.
        integer :: k_start, k_end, p
        logical :: verbose, trans
        real(dp) :: tolerance
        real(dp) :: beta
        complex(dp), allocatable :: res(:)
        integer, allocatable :: perm(:)
        integer :: k, i, kdim, kpm, kp, kpp

        ! Deals with optional non-unity blksize and allocations.
        p = optval(blksize, 1) ; allocate(res(1:p)) ; res = zero_rdp
        allocate(perm(1:size(H, 2))) ; perm = 0 ; info = 0

        ! Check dimensions.
        kdim = (size(X) - p) / p

        ! Deal with the other optional args.
        k_start = optval(kstart, 1) ; k_end = optval(kend, kdim)
        verbose   = optval(verbosity, .false.)
        tolerance = optval (tol, atol_dp)
        trans     = optval(transpose, .false.)

        ! Arnoldi factorization.
        blk_arnoldi: do k = k_start, k_end
            ! Counters
            kpm = (k - 1) * p ; kp = kpm + p ; kpp = kp + p

            ! Matrix-vector product.
            if (trans) then
                do i = 1, p
                    call A%rmatvec(X(kpm+i), X(kp+i))
                enddo
            else
                do i = 1, p
                    call A%matvec(X(kpm+i), X(kp+i))
                enddo
            endif

            ! Update Hessenberg matrix via batch double Gram-Schmidt step.
            call double_gram_schmidt_step(X(kp+1:kpp), X(1:kp), info, if_chk_orthonormal=.false., beta=H(1:kp, kpm+1:kp))
            call check_info(info, 'double_gram_schmidt_step', module=this_module, procedure='arnoldi_cdp')

            ! Orthogonalize current blk vectors.
            call qr(X(kp+1:kpp), H(kp+1:kpp, kpm+1:kp), info)
            call check_info(info, 'qr', module=this_module, procedure='arnoldi_cdp')

            ! Extract residual norm (smallest diagonal element of H matrix).
            res = zero_rdp
            do i = 1, p
                res(i) = H(kp+i, kpm+i)
            enddo
            beta = minval(abs(res))

            ! Exit Arnoldi loop if needed.
            if (beta < tolerance) then
                ! Dimension of the computed invariant subspace.
                info = kp
                ! Exit the Arnoldi iteration.
                exit blk_arnoldi
            endif

        enddo blk_arnoldi

        return
    end subroutine arnoldi_cdp



    !---------------------------------------------
    !-----     LANCZOS BIDIAGONALIZATION     -----
    !---------------------------------------------

    subroutine lanczos_bidiagonalization_rsp(A, U, V, B, info, kstart, kend, verbosity, tol)
        class(abstract_linop_rsp), intent(in) :: A
        !! Linear operator to be factorized.
        class(abstract_vector_rsp), intent(inout) :: U(:)
        !! Orthonormal basis for the column span of \(\mathbf{A}\). On entry, `U(1)` needs to be set to
        !! the starting Krylov vector.
        class(abstract_vector_rsp), intent(inout) :: V(:)
        !! Orthonormal basis for the row span of \(\mathbf{A}\).
        real(sp), intent(inout) :: B(:, :)
        !! Bidiagonal matrix.
        integer, intent(out) :: info
        !! Information flag.
        integer, optional, intent(in) :: kstart
        !! Starting index for the Lanczos factorization (default 1).
        integer, optional, intent(in) :: kend
        !! Final index for the Lanczos factorization (default 1).
        logical, optional, intent(in) :: verbosity
        !! Verbosity control (default `.false.`)
        real(sp), optional, intent(in) :: tol
        !! Tolerance to determine whether invariant subspaces have been computed or not.

        ! Internal variables.
        integer :: k_start, k_end
        logical :: verbose
        real(sp) :: tolerance
        real(sp) :: alpha, beta
        integer :: i, j, k, kdim


        info = 0

        ! Krylov subspace dimension.
        kdim = size(U) - 1

        ! Deals with the optional args.
        k_start = optval(kstart, 1)
        k_end   = optval(kend, kdim)
        verbose = optval(verbosity, .false.)
        tolerance = optval(tol, atol_sp)

        ! Lanczos bidiagonalization.
        lanczos : do k = k_start, k_end
            ! Transpose matrix-vector product.
            call A%rmatvec(U(k), V(k))

            ! Full re-orthogonalization of the right Krylov subspace.
            if (k > 1 ) then
            
                call double_gram_schmidt_step(V(k), V(:k-1), info, if_chk_orthonormal=.true.)
                call check_info(info, 'double_gram_schmidt_step', module=this_module, &
                                    & procedure='lanczos_bidiagonalization_rsp, first pass')
            end if

            ! Normalization step.
            alpha = V(k)%norm() ; B(k, k) = alpha
            if (abs(alpha) > tolerance) then
                call V(k)%scal(one_rsp/alpha)
            else
                info = k
                exit lanczos
            endif

            ! Matrix-vector product.
            call A%matvec(V(k), U(k+1))

            ! Full re-orthogonalization of the left Krylov subspace.
            call double_gram_schmidt_step(U(k+1), U(:k), info, if_chk_orthonormal=.true.)
            call check_info(info, 'orthogonalize_against_basis', module=this_module, &
                                    & procedure='lanczos_bidiagonalization_rsp, second pass')

            ! Normalization step
            beta = U(k+1)%norm() ; B(k+1, k) = beta
            if (abs(beta) > tolerance) then
                call U(k+1)%scal(one_rsp / beta)
            else
                info = k
                exit lanczos
            endif

        enddo lanczos

        return
    end subroutine lanczos_bidiagonalization_rsp

    subroutine lanczos_bidiagonalization_rdp(A, U, V, B, info, kstart, kend, verbosity, tol)
        class(abstract_linop_rdp), intent(in) :: A
        !! Linear operator to be factorized.
        class(abstract_vector_rdp), intent(inout) :: U(:)
        !! Orthonormal basis for the column span of \(\mathbf{A}\). On entry, `U(1)` needs to be set to
        !! the starting Krylov vector.
        class(abstract_vector_rdp), intent(inout) :: V(:)
        !! Orthonormal basis for the row span of \(\mathbf{A}\).
        real(dp), intent(inout) :: B(:, :)
        !! Bidiagonal matrix.
        integer, intent(out) :: info
        !! Information flag.
        integer, optional, intent(in) :: kstart
        !! Starting index for the Lanczos factorization (default 1).
        integer, optional, intent(in) :: kend
        !! Final index for the Lanczos factorization (default 1).
        logical, optional, intent(in) :: verbosity
        !! Verbosity control (default `.false.`)
        real(dp), optional, intent(in) :: tol
        !! Tolerance to determine whether invariant subspaces have been computed or not.

        ! Internal variables.
        integer :: k_start, k_end
        logical :: verbose
        real(dp) :: tolerance
        real(dp) :: alpha, beta
        integer :: i, j, k, kdim


        info = 0

        ! Krylov subspace dimension.
        kdim = size(U) - 1

        ! Deals with the optional args.
        k_start = optval(kstart, 1)
        k_end   = optval(kend, kdim)
        verbose = optval(verbosity, .false.)
        tolerance = optval(tol, atol_dp)

        ! Lanczos bidiagonalization.
        lanczos : do k = k_start, k_end
            ! Transpose matrix-vector product.
            call A%rmatvec(U(k), V(k))

            ! Full re-orthogonalization of the right Krylov subspace.
            if (k > 1 ) then
            
                call double_gram_schmidt_step(V(k), V(:k-1), info, if_chk_orthonormal=.true.)
                call check_info(info, 'double_gram_schmidt_step', module=this_module, &
                                    & procedure='lanczos_bidiagonalization_rdp, first pass')
            end if

            ! Normalization step.
            alpha = V(k)%norm() ; B(k, k) = alpha
            if (abs(alpha) > tolerance) then
                call V(k)%scal(one_rdp/alpha)
            else
                info = k
                exit lanczos
            endif

            ! Matrix-vector product.
            call A%matvec(V(k), U(k+1))

            ! Full re-orthogonalization of the left Krylov subspace.
            call double_gram_schmidt_step(U(k+1), U(:k), info, if_chk_orthonormal=.true.)
            call check_info(info, 'orthogonalize_against_basis', module=this_module, &
                                    & procedure='lanczos_bidiagonalization_rdp, second pass')

            ! Normalization step
            beta = U(k+1)%norm() ; B(k+1, k) = beta
            if (abs(beta) > tolerance) then
                call U(k+1)%scal(one_rdp / beta)
            else
                info = k
                exit lanczos
            endif

        enddo lanczos

        return
    end subroutine lanczos_bidiagonalization_rdp

    subroutine lanczos_bidiagonalization_csp(A, U, V, B, info, kstart, kend, verbosity, tol)
        class(abstract_linop_csp), intent(in) :: A
        !! Linear operator to be factorized.
        class(abstract_vector_csp), intent(inout) :: U(:)
        !! Orthonormal basis for the column span of \(\mathbf{A}\). On entry, `U(1)` needs to be set to
        !! the starting Krylov vector.
        class(abstract_vector_csp), intent(inout) :: V(:)
        !! Orthonormal basis for the row span of \(\mathbf{A}\).
        complex(sp), intent(inout) :: B(:, :)
        !! Bidiagonal matrix.
        integer, intent(out) :: info
        !! Information flag.
        integer, optional, intent(in) :: kstart
        !! Starting index for the Lanczos factorization (default 1).
        integer, optional, intent(in) :: kend
        !! Final index for the Lanczos factorization (default 1).
        logical, optional, intent(in) :: verbosity
        !! Verbosity control (default `.false.`)
        real(sp), optional, intent(in) :: tol
        !! Tolerance to determine whether invariant subspaces have been computed or not.

        ! Internal variables.
        integer :: k_start, k_end
        logical :: verbose
        real(sp) :: tolerance
        complex(sp) :: alpha, beta
        integer :: i, j, k, kdim


        info = 0

        ! Krylov subspace dimension.
        kdim = size(U) - 1

        ! Deals with the optional args.
        k_start = optval(kstart, 1)
        k_end   = optval(kend, kdim)
        verbose = optval(verbosity, .false.)
        tolerance = optval(tol, atol_sp)

        ! Lanczos bidiagonalization.
        lanczos : do k = k_start, k_end
            ! Transpose matrix-vector product.
            call A%rmatvec(U(k), V(k))

            ! Full re-orthogonalization of the right Krylov subspace.
            if (k > 1 ) then
            
                call double_gram_schmidt_step(V(k), V(:k-1), info, if_chk_orthonormal=.true.)
                call check_info(info, 'double_gram_schmidt_step', module=this_module, &
                                    & procedure='lanczos_bidiagonalization_csp, first pass')
            end if

            ! Normalization step.
            alpha = V(k)%norm() ; B(k, k) = alpha
            if (abs(alpha) > tolerance) then
                call V(k)%scal(one_csp/alpha)
            else
                info = k
                exit lanczos
            endif

            ! Matrix-vector product.
            call A%matvec(V(k), U(k+1))

            ! Full re-orthogonalization of the left Krylov subspace.
            call double_gram_schmidt_step(U(k+1), U(:k), info, if_chk_orthonormal=.true.)
            call check_info(info, 'orthogonalize_against_basis', module=this_module, &
                                    & procedure='lanczos_bidiagonalization_csp, second pass')

            ! Normalization step
            beta = U(k+1)%norm() ; B(k+1, k) = beta
            if (abs(beta) > tolerance) then
                call U(k+1)%scal(one_csp / beta)
            else
                info = k
                exit lanczos
            endif

        enddo lanczos

        return
    end subroutine lanczos_bidiagonalization_csp

    subroutine lanczos_bidiagonalization_cdp(A, U, V, B, info, kstart, kend, verbosity, tol)
        class(abstract_linop_cdp), intent(in) :: A
        !! Linear operator to be factorized.
        class(abstract_vector_cdp), intent(inout) :: U(:)
        !! Orthonormal basis for the column span of \(\mathbf{A}\). On entry, `U(1)` needs to be set to
        !! the starting Krylov vector.
        class(abstract_vector_cdp), intent(inout) :: V(:)
        !! Orthonormal basis for the row span of \(\mathbf{A}\).
        complex(dp), intent(inout) :: B(:, :)
        !! Bidiagonal matrix.
        integer, intent(out) :: info
        !! Information flag.
        integer, optional, intent(in) :: kstart
        !! Starting index for the Lanczos factorization (default 1).
        integer, optional, intent(in) :: kend
        !! Final index for the Lanczos factorization (default 1).
        logical, optional, intent(in) :: verbosity
        !! Verbosity control (default `.false.`)
        real(dp), optional, intent(in) :: tol
        !! Tolerance to determine whether invariant subspaces have been computed or not.

        ! Internal variables.
        integer :: k_start, k_end
        logical :: verbose
        real(dp) :: tolerance
        complex(dp) :: alpha, beta
        integer :: i, j, k, kdim


        info = 0

        ! Krylov subspace dimension.
        kdim = size(U) - 1

        ! Deals with the optional args.
        k_start = optval(kstart, 1)
        k_end   = optval(kend, kdim)
        verbose = optval(verbosity, .false.)
        tolerance = optval(tol, atol_dp)

        ! Lanczos bidiagonalization.
        lanczos : do k = k_start, k_end
            ! Transpose matrix-vector product.
            call A%rmatvec(U(k), V(k))

            ! Full re-orthogonalization of the right Krylov subspace.
            if (k > 1 ) then
            
                call double_gram_schmidt_step(V(k), V(:k-1), info, if_chk_orthonormal=.true.)
                call check_info(info, 'double_gram_schmidt_step', module=this_module, &
                                    & procedure='lanczos_bidiagonalization_cdp, first pass')
            end if

            ! Normalization step.
            alpha = V(k)%norm() ; B(k, k) = alpha
            if (abs(alpha) > tolerance) then
                call V(k)%scal(one_cdp/alpha)
            else
                info = k
                exit lanczos
            endif

            ! Matrix-vector product.
            call A%matvec(V(k), U(k+1))

            ! Full re-orthogonalization of the left Krylov subspace.
            call double_gram_schmidt_step(U(k+1), U(:k), info, if_chk_orthonormal=.true.)
            call check_info(info, 'orthogonalize_against_basis', module=this_module, &
                                    & procedure='lanczos_bidiagonalization_cdp, second pass')

            ! Normalization step
            beta = U(k+1)%norm() ; B(k+1, k) = beta
            if (abs(beta) > tolerance) then
                call U(k+1)%scal(one_cdp / beta)
            else
                info = k
                exit lanczos
            endif

        enddo lanczos

        return
    end subroutine lanczos_bidiagonalization_cdp



    !----------------------------------------------
    !-----     LANCZOS TRIDIAGONALIZATION     -----
    !----------------------------------------------
    
    subroutine lanczos_tridiagonalization_rsp(A, X, T, info, kstart, kend, verbosity, tol)
        class(abstract_sym_linop_rsp), intent(in) :: A
        class(abstract_vector_rsp), intent(inout) :: X(:)
        real(sp), intent(inout) :: T(:, :)
        integer, intent(out) :: info
        integer, optional, intent(in) :: kstart
        integer, optional, intent(in) :: kend
        logical, optional, intent(in) :: verbosity
        real(sp), optional, intent(in) :: tol

        ! Internal variables.
        integer :: k_start, k_end
        logical :: verbose
        real(sp) :: tolerance
        real(sp) :: beta
        integer :: i, j, k, kdim

        ! Deal with optional args.
        kdim = size(X) - 1
        k_start = optval(kstart, 1)
        k_end = optval(kend, kdim)
        verbose = optval(verbosity, .false.)
        tolerance = optval(tol, atol_sp)

        ! Lanczos tridiagonalization.
        lanczos: do k = k_start, k_end
            ! Matrix-vector product.
            call A%matvec(X(k), X(k+1))
            ! Update tridiagonal matrix.
            call update_tridiag_matrix_rsp(T, X, k)
            beta = X(k+1)%norm() ; T(k+1, k) = beta

            ! Exit Lanczos loop if needed.
            if (beta < tolerance) then
                ! Dimension of the computed invariant subspace.
                info = k
                ! Exit the Lanczos iteration.
                exit lanczos
            else
                ! Normalize the new Krylov vector.
                call X(k+1)%scal(one_rsp / beta)
            endif
        enddo lanczos

        return
    end subroutine lanczos_tridiagonalization_rsp

    subroutine update_tridiag_matrix_rsp(T, X, k)
        integer, intent(in) :: k
        real(sp), intent(inout) :: T(:, :)
        class(abstract_vector_rsp), intent(inout) :: X(:)

        ! Internal variables.
        class(abstract_vector_rsp), allocatable :: wrk
        integer :: i, info

        info = 0

        ! Orthogonalize residual w.r.t. previously computed Krylov vectors to obtain coefficients in tridiag. matrix
        do i = max(1, k-1), k
            T(i, k) = X(i)%dot(X(k+1)) ; call X(k+1)%axpby(one_rsp, X(i), -T(i, k))
        enddo

        ! Full re-orthogonalization against existing basis
        call orthogonalize_against_basis(X(k+1), X(1:k), info, if_chk_orthonormal=.false.)
        call check_info(info, 'orthogonalize_against_basis', module=this_module, procedure='update_tridiag_matrix_rsp')

        return
    end subroutine update_tridiag_matrix_rsp

    subroutine lanczos_tridiagonalization_rdp(A, X, T, info, kstart, kend, verbosity, tol)
        class(abstract_sym_linop_rdp), intent(in) :: A
        class(abstract_vector_rdp), intent(inout) :: X(:)
        real(dp), intent(inout) :: T(:, :)
        integer, intent(out) :: info
        integer, optional, intent(in) :: kstart
        integer, optional, intent(in) :: kend
        logical, optional, intent(in) :: verbosity
        real(dp), optional, intent(in) :: tol

        ! Internal variables.
        integer :: k_start, k_end
        logical :: verbose
        real(dp) :: tolerance
        real(dp) :: beta
        integer :: i, j, k, kdim

        ! Deal with optional args.
        kdim = size(X) - 1
        k_start = optval(kstart, 1)
        k_end = optval(kend, kdim)
        verbose = optval(verbosity, .false.)
        tolerance = optval(tol, atol_dp)

        ! Lanczos tridiagonalization.
        lanczos: do k = k_start, k_end
            ! Matrix-vector product.
            call A%matvec(X(k), X(k+1))
            ! Update tridiagonal matrix.
            call update_tridiag_matrix_rdp(T, X, k)
            beta = X(k+1)%norm() ; T(k+1, k) = beta

            ! Exit Lanczos loop if needed.
            if (beta < tolerance) then
                ! Dimension of the computed invariant subspace.
                info = k
                ! Exit the Lanczos iteration.
                exit lanczos
            else
                ! Normalize the new Krylov vector.
                call X(k+1)%scal(one_rdp / beta)
            endif
        enddo lanczos

        return
    end subroutine lanczos_tridiagonalization_rdp

    subroutine update_tridiag_matrix_rdp(T, X, k)
        integer, intent(in) :: k
        real(dp), intent(inout) :: T(:, :)
        class(abstract_vector_rdp), intent(inout) :: X(:)

        ! Internal variables.
        class(abstract_vector_rdp), allocatable :: wrk
        integer :: i, info

        info = 0

        ! Orthogonalize residual w.r.t. previously computed Krylov vectors to obtain coefficients in tridiag. matrix
        do i = max(1, k-1), k
            T(i, k) = X(i)%dot(X(k+1)) ; call X(k+1)%axpby(one_rdp, X(i), -T(i, k))
        enddo

        ! Full re-orthogonalization against existing basis
        call orthogonalize_against_basis(X(k+1), X(1:k), info, if_chk_orthonormal=.false.)
        call check_info(info, 'orthogonalize_against_basis', module=this_module, procedure='update_tridiag_matrix_rdp')

        return
    end subroutine update_tridiag_matrix_rdp

    subroutine lanczos_tridiagonalization_csp(A, X, T, info, kstart, kend, verbosity, tol)
        class(abstract_hermitian_linop_csp), intent(in) :: A
        class(abstract_vector_csp), intent(inout) :: X(:)
        complex(sp), intent(inout) :: T(:, :)
        integer, intent(out) :: info
        integer, optional, intent(in) :: kstart
        integer, optional, intent(in) :: kend
        logical, optional, intent(in) :: verbosity
        real(sp), optional, intent(in) :: tol

        ! Internal variables.
        integer :: k_start, k_end
        logical :: verbose
        real(sp) :: tolerance
        real(sp) :: beta
        integer :: i, j, k, kdim

        ! Deal with optional args.
        kdim = size(X) - 1
        k_start = optval(kstart, 1)
        k_end = optval(kend, kdim)
        verbose = optval(verbosity, .false.)
        tolerance = optval(tol, atol_sp)

        ! Lanczos tridiagonalization.
        lanczos: do k = k_start, k_end
            ! Matrix-vector product.
            call A%matvec(X(k), X(k+1))
            ! Update tridiagonal matrix.
            call update_tridiag_matrix_csp(T, X, k)
            beta = X(k+1)%norm() ; T(k+1, k) = beta

            ! Exit Lanczos loop if needed.
            if (beta < tolerance) then
                ! Dimension of the computed invariant subspace.
                info = k
                ! Exit the Lanczos iteration.
                exit lanczos
            else
                ! Normalize the new Krylov vector.
                call X(k+1)%scal(one_csp / beta)
            endif
        enddo lanczos

        return
    end subroutine lanczos_tridiagonalization_csp

    subroutine update_tridiag_matrix_csp(T, X, k)
        integer, intent(in) :: k
        complex(sp), intent(inout) :: T(:, :)
        class(abstract_vector_csp), intent(inout) :: X(:)

        ! Internal variables.
        class(abstract_vector_csp), allocatable :: wrk
        integer :: i, info

        info = 0

        ! Orthogonalize residual w.r.t. previously computed Krylov vectors to obtain coefficients in tridiag. matrix
        do i = max(1, k-1), k
            T(i, k) = X(i)%dot(X(k+1)) ; call X(k+1)%axpby(one_csp, X(i), -T(i, k))
        enddo

        ! Full re-orthogonalization against existing basis
        call orthogonalize_against_basis(X(k+1), X(1:k), info, if_chk_orthonormal=.false.)
        call check_info(info, 'orthogonalize_against_basis', module=this_module, procedure='update_tridiag_matrix_csp')

        return
    end subroutine update_tridiag_matrix_csp

    subroutine lanczos_tridiagonalization_cdp(A, X, T, info, kstart, kend, verbosity, tol)
        class(abstract_hermitian_linop_cdp), intent(in) :: A
        class(abstract_vector_cdp), intent(inout) :: X(:)
        complex(dp), intent(inout) :: T(:, :)
        integer, intent(out) :: info
        integer, optional, intent(in) :: kstart
        integer, optional, intent(in) :: kend
        logical, optional, intent(in) :: verbosity
        real(dp), optional, intent(in) :: tol

        ! Internal variables.
        integer :: k_start, k_end
        logical :: verbose
        real(dp) :: tolerance
        real(dp) :: beta
        integer :: i, j, k, kdim

        ! Deal with optional args.
        kdim = size(X) - 1
        k_start = optval(kstart, 1)
        k_end = optval(kend, kdim)
        verbose = optval(verbosity, .false.)
        tolerance = optval(tol, atol_dp)

        ! Lanczos tridiagonalization.
        lanczos: do k = k_start, k_end
            ! Matrix-vector product.
            call A%matvec(X(k), X(k+1))
            ! Update tridiagonal matrix.
            call update_tridiag_matrix_cdp(T, X, k)
            beta = X(k+1)%norm() ; T(k+1, k) = beta

            ! Exit Lanczos loop if needed.
            if (beta < tolerance) then
                ! Dimension of the computed invariant subspace.
                info = k
                ! Exit the Lanczos iteration.
                exit lanczos
            else
                ! Normalize the new Krylov vector.
                call X(k+1)%scal(one_cdp / beta)
            endif
        enddo lanczos

        return
    end subroutine lanczos_tridiagonalization_cdp

    subroutine update_tridiag_matrix_cdp(T, X, k)
        integer, intent(in) :: k
        complex(dp), intent(inout) :: T(:, :)
        class(abstract_vector_cdp), intent(inout) :: X(:)

        ! Internal variables.
        class(abstract_vector_cdp), allocatable :: wrk
        integer :: i, info

        info = 0

        ! Orthogonalize residual w.r.t. previously computed Krylov vectors to obtain coefficients in tridiag. matrix
        do i = max(1, k-1), k
            T(i, k) = X(i)%dot(X(k+1)) ; call X(k+1)%axpby(one_cdp, X(i), -T(i, k))
        enddo

        ! Full re-orthogonalization against existing basis
        call orthogonalize_against_basis(X(k+1), X(1:k), info, if_chk_orthonormal=.false.)
        call check_info(info, 'orthogonalize_against_basis', module=this_module, procedure='update_tridiag_matrix_cdp')

        return
    end subroutine update_tridiag_matrix_cdp


    !----------------------------------------
    !-----     KRYLOV-SCHUR RESTART     -----
    !----------------------------------------

    subroutine krylov_schur_rsp(n, X, H, select_eigs)
        integer, intent(out) :: n
        !! Number eigenvalues that have been moved to the upper
        !! left block of the Schur factorization of `H`.
        class(abstract_vector_rsp), intent(inout) :: X(:)
        !! Krylov basis.
        real(sp), intent(inout) :: H(:, :)
        !! Upper Hessenberg matrix.
        procedure(eigvals_select_sp) :: select_eigs
        !! Procedure to select the eigenvalues to move in the upper left-block.

        !--------------------------------------
        !-----     Internal variables     -----
        !--------------------------------------

        integer :: i, j, kdim
        
        ! Schur-related.
        real(sp) :: Z(size(H, 2), size(H, 2))
        complex(sp) :: eigvals(size(H, 2))
        logical :: selected(size(H, 2))
       
        ! Krylov subspace dimension.
        kdim = size(X)-1

        ! Allocate and initializes variables.
        eigvals = zero_rsp ; Z = zero_rsp
        selected = .false.

        ! Schur decomposition of the Hessenberg matrix.
        call schur(H(:size(H, 2), :), Z, eigvals)

        ! Eigenvalue selection of the upper left block.
        selected = select_eigs(eigvals) ; n = count(selected)

        ! Re-order the Schur decomposition and Schur basis.
        call ordschur(H(:kdim, :), Z, selected)

        ! Update the Hessenberg matrix and Krylov basis.
        block
        real(sp) :: b(size(H, 2))
        class(abstract_vector_rsp), allocatable :: Xwrk(:)
        
        ! Update the Krylov basis.
        call linear_combination(Xwrk, X(:size(H, 2)), Z(:, :n))

        do i = 1, n
            call X(i)%axpby(zero_rsp, Xwrk(i), one_rsp)
        enddo

        call X(n+1)%axpby(zero_rsp, X(kdim+1), one_rsp)

        call zero_basis(X(n+2:))

        ! Update the Hessenberg matrix.
        b = matmul(H(kdim+1, :), Z)
        H(n+1, :) = b
        H(n+2:, :) = zero_rsp
        H(:, n+1:) = zero_rsp
        end block

        return
    end subroutine krylov_schur_rsp

    subroutine krylov_schur_rdp(n, X, H, select_eigs)
        integer, intent(out) :: n
        !! Number eigenvalues that have been moved to the upper
        !! left block of the Schur factorization of `H`.
        class(abstract_vector_rdp), intent(inout) :: X(:)
        !! Krylov basis.
        real(dp), intent(inout) :: H(:, :)
        !! Upper Hessenberg matrix.
        procedure(eigvals_select_dp) :: select_eigs
        !! Procedure to select the eigenvalues to move in the upper left-block.

        !--------------------------------------
        !-----     Internal variables     -----
        !--------------------------------------

        integer :: i, j, kdim
        
        ! Schur-related.
        real(dp) :: Z(size(H, 2), size(H, 2))
        complex(dp) :: eigvals(size(H, 2))
        logical :: selected(size(H, 2))
       
        ! Krylov subspace dimension.
        kdim = size(X)-1

        ! Allocate and initializes variables.
        eigvals = zero_rdp ; Z = zero_rdp
        selected = .false.

        ! Schur decomposition of the Hessenberg matrix.
        call schur(H(:size(H, 2), :), Z, eigvals)

        ! Eigenvalue selection of the upper left block.
        selected = select_eigs(eigvals) ; n = count(selected)

        ! Re-order the Schur decomposition and Schur basis.
        call ordschur(H(:kdim, :), Z, selected)

        ! Update the Hessenberg matrix and Krylov basis.
        block
        real(dp) :: b(size(H, 2))
        class(abstract_vector_rdp), allocatable :: Xwrk(:)
        
        ! Update the Krylov basis.
        call linear_combination(Xwrk, X(:size(H, 2)), Z(:, :n))

        do i = 1, n
            call X(i)%axpby(zero_rdp, Xwrk(i), one_rdp)
        enddo

        call X(n+1)%axpby(zero_rdp, X(kdim+1), one_rdp)

        call zero_basis(X(n+2:))

        ! Update the Hessenberg matrix.
        b = matmul(H(kdim+1, :), Z)
        H(n+1, :) = b
        H(n+2:, :) = zero_rdp
        H(:, n+1:) = zero_rdp
        end block

        return
    end subroutine krylov_schur_rdp

    subroutine krylov_schur_csp(n, X, H, select_eigs)
        integer, intent(out) :: n
        !! Number eigenvalues that have been moved to the upper
        !! left block of the Schur factorization of `H`.
        class(abstract_vector_csp), intent(inout) :: X(:)
        !! Krylov basis.
        complex(sp), intent(inout) :: H(:, :)
        !! Upper Hessenberg matrix.
        procedure(eigvals_select_sp) :: select_eigs
        !! Procedure to select the eigenvalues to move in the upper left-block.

        !--------------------------------------
        !-----     Internal variables     -----
        !--------------------------------------

        integer :: i, j, kdim
        
        ! Schur-related.
        complex(sp) :: Z(size(H, 2), size(H, 2))
        complex(sp) :: eigvals(size(H, 2))
        logical :: selected(size(H, 2))
       
        ! Krylov subspace dimension.
        kdim = size(X)-1

        ! Allocate and initializes variables.
        eigvals = zero_csp ; Z = zero_csp
        selected = .false.

        ! Schur decomposition of the Hessenberg matrix.
        call schur(H(:size(H, 2), :), Z, eigvals)

        ! Eigenvalue selection of the upper left block.
        selected = select_eigs(eigvals) ; n = count(selected)

        ! Re-order the Schur decomposition and Schur basis.
        call ordschur(H(:kdim, :), Z, selected)

        ! Update the Hessenberg matrix and Krylov basis.
        block
        complex(sp) :: b(size(H, 2))
        class(abstract_vector_csp), allocatable :: Xwrk(:)
        
        ! Update the Krylov basis.
        call linear_combination(Xwrk, X(:size(H, 2)), Z(:, :n))

        do i = 1, n
            call X(i)%axpby(zero_csp, Xwrk(i), one_csp)
        enddo

        call X(n+1)%axpby(zero_csp, X(kdim+1), one_csp)

        call zero_basis(X(n+2:))

        ! Update the Hessenberg matrix.
        b = matmul(H(kdim+1, :), Z)
        H(n+1, :) = b
        H(n+2:, :) = zero_csp
        H(:, n+1:) = zero_csp
        end block

        return
    end subroutine krylov_schur_csp

    subroutine krylov_schur_cdp(n, X, H, select_eigs)
        integer, intent(out) :: n
        !! Number eigenvalues that have been moved to the upper
        !! left block of the Schur factorization of `H`.
        class(abstract_vector_cdp), intent(inout) :: X(:)
        !! Krylov basis.
        complex(dp), intent(inout) :: H(:, :)
        !! Upper Hessenberg matrix.
        procedure(eigvals_select_dp) :: select_eigs
        !! Procedure to select the eigenvalues to move in the upper left-block.

        !--------------------------------------
        !-----     Internal variables     -----
        !--------------------------------------

        integer :: i, j, kdim
        
        ! Schur-related.
        complex(dp) :: Z(size(H, 2), size(H, 2))
        complex(dp) :: eigvals(size(H, 2))
        logical :: selected(size(H, 2))
       
        ! Krylov subspace dimension.
        kdim = size(X)-1

        ! Allocate and initializes variables.
        eigvals = zero_cdp ; Z = zero_cdp
        selected = .false.

        ! Schur decomposition of the Hessenberg matrix.
        call schur(H(:size(H, 2), :), Z, eigvals)

        ! Eigenvalue selection of the upper left block.
        selected = select_eigs(eigvals) ; n = count(selected)

        ! Re-order the Schur decomposition and Schur basis.
        call ordschur(H(:kdim, :), Z, selected)

        ! Update the Hessenberg matrix and Krylov basis.
        block
        complex(dp) :: b(size(H, 2))
        class(abstract_vector_cdp), allocatable :: Xwrk(:)
        
        ! Update the Krylov basis.
        call linear_combination(Xwrk, X(:size(H, 2)), Z(:, :n))

        do i = 1, n
            call X(i)%axpby(zero_cdp, Xwrk(i), one_cdp)
        enddo

        call X(n+1)%axpby(zero_cdp, X(kdim+1), one_cdp)

        call zero_basis(X(n+2:))

        ! Update the Hessenberg matrix.
        b = matmul(H(kdim+1, :), Z)
        H(n+1, :) = b
        H(n+2:, :) = zero_cdp
        H(:, n+1:) = zero_cdp
        end block

        return
    end subroutine krylov_schur_cdp


end module lightkrylov_BaseKrylov
