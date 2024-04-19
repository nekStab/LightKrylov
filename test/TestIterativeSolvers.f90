module TestIterativeSolvers
   use LightKrylov
   use TestVector
   use TestMatrices
   use TestUtils
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use stdlib_math, only: all_close
   implicit none

   private

   public :: collect_gmres_testsuite, &
             collect_cg_testsuite, &
             collect_evp_testsuite, &
             collect_svd_testsuite

   !-----------------------------------------------------------
   !-----                                                 -----
   !-----     JACOBI-BASED PRECONDITIONER FOR TESTING     -----
   !-----                                                 -----
   !-----------------------------------------------------------

   type, extends(abstract_preconditioner), public :: jacobi_preconditioner
      real(kind=wp), dimension(test_size) :: data
   contains
      private
      procedure, pass(self), public :: apply
      procedure, pass(self), public :: undo
   end type jacobi_preconditioner

contains

   !--------------------------------------------------------------------
   !-----                                                          -----
   !-----     TYPE-BOUND PROCEDURES FOR JACOBI-PRECONDITIONING     -----
   !-----                                                          -----
   !--------------------------------------------------------------------

   subroutine apply(self, vec_inout)
      ! Preconditioner.
      class(jacobi_preconditioner), intent(in)    :: self
      ! Input/output vector.
      class(abstract_vector), intent(inout) :: vec_inout

      ! Diagonal scaling.
      select type (vec_inout)
      type is (rvector)
         vec_inout%data = (1.0_wp/self%data)*vec_inout%data
      end select
      return
   end subroutine apply

   subroutine undo(self, vec_inout)
      ! Preconditioner.
      class(jacobi_preconditioner), intent(in)    :: self
      ! Input/Output vector.
      class(abstract_vector), intent(inout) :: vec_inout

      ! Undo diagonal scaling.
      select type (vec_inout)
      type is (rvector)
         vec_inout%data = self%data*vec_inout%data
      end select
      return
   end subroutine undo

   !-------------------------------------------------------------
   !-----                                                   -----
   !-----     TEST SUITE FOR GENERAL EIGENVALUE PROBLEM     -----
   !-----                                                   -----
   !-------------------------------------------------------------

   subroutine collect_evp_testsuite(testsuite)
      ! Collection of tests.
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("Sym. EVP computation", test_spd_evp_problem), &
                  new_unittest("EVP computation", test_evp_problem) &
                  ]
      return
   end subroutine collect_evp_testsuite

   subroutine test_spd_evp_problem(error)
      ! Error type to be returned.
      type(error_type), allocatable, intent(out) :: error
      ! Test matrix.
      class(spd_matrix), allocatable :: A
      ! Krylov subspace.
      class(rvector), allocatable :: X(:)
      ! Eigenvalues and eigenvectors.
      real(kind=wp) :: evals(test_size)
      ! Residuals.
      real(kind=wp) :: residuals(test_size)
      ! Information flag.
      integer :: info
      ! Toeplitz matrix.
      real(kind=wp) :: T(test_size, test_size), a_, b_
      ! Miscellaneous.
      integer :: i
      real(kind=wp) :: alpha
      real(kind=wp) :: true_evals(test_size), pi
      class(rvector), allocatable :: X0(1)

      ! Create the sym. pos. def. Toeplitz matrix.
      call random_number(a_); call random_number(b_); b_ = -abs(b_)
      do i = 1, test_size
         ! Diagonal entry.
         T(i, i) = a_
         ! Upper diagonal entry.
         if (i < test_size) T(i, i + 1) = b_
         ! Lower diagonal entry.
         if (i < test_size) T(i + 1, i) = b_
      end do

      ! Test matrix.
      A = spd_matrix(T)
      ! Initialize Krylov subspace.
      allocate (X(1:test_size + 1)); allocate (X0(1))
      call init_rand(X0)
      call initialize_krylov_subspace(X, X0)

      ! Initialize internal variables.
      evals = 0.0_wp; residuals = 0.0_wp

      ! Compute spectral decomposition.
      call eighs(A, X, evals, residuals, info)

      ! Analytical eigenvalues.
      true_evals = 0.0_wp; pi = 4.0_wp*atan(1.0_wp)
      do i = 1, test_size
         true_evals(i) = a_ + 2*abs(b_)*cos(i*pi/(test_size + 1))
      end do

      ! Check correctness.
      call check(error, all_close(evals, true_evals, rtol, atol))

      return
   end subroutine test_spd_evp_problem

   subroutine test_evp_problem(error)
      ! Error type to be returned.
      type(error_type), allocatable, intent(out) :: error
      ! Test matrix.
      class(rmatrix), allocatable :: A
      ! Krylov subspace.
      class(rvector), allocatable :: X(:)
      ! Eigenvalues and eigenvectors.
      complex(kind=wp) :: evals(test_size)
      ! Residuals.
      real(kind=wp) :: residuals(test_size)
      ! Information flag.
      integer :: info
      ! Toeplitz matrix.
      real(kind=wp) :: T(test_size, test_size), a_, b_
      ! Miscellaneous.
      integer :: i, k
      real(kind=wp) :: alpha
      complex(kind=wp) :: true_evals(test_size), pi
      class(rvector), allocatable :: X0(1)

      ! Create the sym. pos. def. Toeplitz matrix.
      call random_number(a_); call random_number(b_); b_ = abs(b_)
      do i = 1, test_size
         ! Diagonal entry.
         T(i, i) = a_
         ! Upper diagonal entry.
         if (i < test_size) T(i, i + 1) = b_
         ! Lower diagonal entry.
         if (i < test_size) T(i + 1, i) = -b_
      end do

      ! Test matrix.
      A = rmatrix(T)

      ! Initialize Krylov subspace.
      allocate (X(1:test_size + 1)); allocate (X0(1))
      call init_rand(X0)
      call initialize_krylov_subspace(X, X0)

      ! Initialize internal variables.
      evals = cmplx(0.0_wp, 0.0_wp, kind=wp); residuals = 0.0_wp

      ! Compute spectral decomposition.
      call eigs(A, X, evals, residuals, info)

      ! Analytical eigenvalues.
      true_evals = cmplx(0.0_wp, 0.0_wp, kind=wp); pi = 4.0_wp*atan(1.0_wp)
      k = 1
      do i = 1, test_size, 2
         true_evals(i) = a_*cmplx(1.0_wp, 0.0_wp, kind=wp) + (2_wp*b_*cos(k*pi/(test_size + 1)))*cmplx(0.0_wp, 1.0_wp, kind=wp)
         true_evals(i + 1) = a_*cmplx(1.0_wp, 0.0_wp, kind=wp) - (2_wp*b_*cos(k*pi/(test_size + 1)))*cmplx(0.0_wp, 1.0_wp, kind=wp)
         k = k + 1
      end do

      ! Check correctness.
      call check(error, norm2(abs(evals - true_evals)) < rtol)

      return
   end subroutine test_evp_problem

   !------------------------------------
   !-----                          -----
   !-----     GMRES TEST SUITE     -----
   !-----                          -----
   !------------------------------------

   subroutine collect_gmres_testsuite(testsuite)
      ! Collection of tests.
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("GMRES full computation", test_gmres_full_computation), &
                  new_unittest("GMRES full computation w. sym. pos. def. matrix", test_gmres_full_computation_spd_matrix), &
                  new_unittest("GMRES + Jacobi precond.", test_precond_gmres_full_computation) &
                  ]
      return
   end subroutine collect_gmres_testsuite

   subroutine test_gmres_full_computation(error)
      ! Error type to be returned.
      type(error_type), allocatable, intent(out) :: error
      ! Linear Problem.
      class(rmatrix), allocatable :: A ! Linear Operator.
      class(rvector), allocatable :: b ! Right-hand side vector.
      class(rvector), allocatable :: x ! Solution vector.
      ! GMRES options.
      type(gmres_opts) :: opts
      ! Information flag.
      integer :: info

      ! Initialize linear problem.
      A = rmatrix(); call init_rand(A)
      b = rvector(); call init_rand(b)
      x = rvector(); call x%zero()
      ! GMRES solver.
      opts = gmres_opts(kdim=test_size, verbose=.false.)
      call gmres(A, b, x, info, options=opts)
      ! Check convergence.
      call check(error, norm2(matmul(A%data, x%data) - b%data)**2 < rtol)

      return
   end subroutine test_gmres_full_computation

   subroutine test_gmres_full_computation_spd_matrix(error)
      ! Error type to be returned.
      type(error_type), allocatable, intent(out) :: error
      ! Linear Problem.
      class(spd_matrix), allocatable :: A ! Linear Operator.
      class(rvector), allocatable :: b ! Right-hand side vector.
      class(rvector), allocatable :: x ! Solution vector.
      ! GMRES options.
      type(gmres_opts) :: opts
      ! Information flag.
      integer :: info

      ! Initialize linear problem.
      A = spd_matrix(); call init_rand(A)
      b = rvector(); call init_rand(b)
      x = rvector(); call x%zero()
      ! GMRES solver.
      opts = gmres_opts(kdim=test_size, verbose=.false.)
      call gmres(A, b, x, info, options=opts)
      ! Check convergence.
      call check(error, norm2(matmul(A%data, x%data) - b%data)**2 < rtol)

      return
   end subroutine test_gmres_full_computation_spd_matrix

   subroutine test_precond_gmres_full_computation(error)
      ! Error type to be returned.
      type(error_type), allocatable, intent(out) :: error
      ! Linear Problem.
      class(rmatrix), allocatable :: A ! Linear Operator.
      class(rvector), allocatable :: b ! Right-hand side vector.
      class(rvector), allocatable :: x ! Solution vector.
      ! GMRES options.
      type(gmres_opts) :: opts
      ! Preconditioner.
      type(jacobi_preconditioner), allocatable :: D
      real(kind=wp) :: diag(test_size)
      ! Information flag.
      integer :: info
      ! Miscellaneous.
      integer :: i

      ! Initialize linear problem.
      A = rmatrix(); call init_rand(A)
      b = rvector(); call init_rand(b)
      x = rvector(); call x%zero()
      ! Preconditioner.
      diag = 0.0_wp
      do i = 1, test_size
         diag(i) = A%data(i, i)
      end do
      D = jacobi_preconditioner(diag)
      ! GMRES solver.
      opts = gmres_opts(kdim=test_size, verbose=.false.)
      write (*, *) "-- Running Unpreconditioned GMRES :"
      call gmres(A, b, x, info, options=opts)
      write (*, *)
      write (*, *) "-- Running GMRES with Jacobi preconditioning."
      call x%zero()
      call gmres(A, b, x, info, options=opts, preconditioner=D)
      ! Check convergence.
      call check(error, norm2(matmul(A%data, x%data) - b%data)**2 < rtol)

      return
   end subroutine test_precond_gmres_full_computation

   !----------------------------------
   !-----                        -----
   !-----     SVD TEST SUITE     -----
   !-----                        -----
   !----------------------------------

   subroutine collect_svd_testsuite(testsuite)
      ! Collection of tests.
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("SVD computation", test_svd_problem) &
                  ]
      return
   end subroutine collect_svd_testsuite

   subroutine test_svd_problem(error)
      ! Error type to be returned.
      type(error_type), allocatable, intent(out) :: error
      ! Test matrix.
      class(rmatrix), allocatable :: A
      ! Krylov subspaces.
      class(rvector), allocatable :: U(:), V(:)
      integer                     :: kdim = test_size
      ! Singular values.
      real(kind=wp) :: svdvals(1:test_size)
      ! Residuals.
      real(kind=wp) :: residuals(1:test_size)
      ! Information flag.
      integer :: info
      ! Miscellaneous.
      integer :: i
      real(kind=wp) :: alpha
      real(kind=wp) :: true_svdvals(1:test_size), pi = 4.0_wp*atan(1.0_wp)
      class(rvector), allocatable :: U0(1)

      ! Initialize matrix.
      A = rmatrix()
      do i = 1, test_size
         ! Diagonal entry.
         A%data(i, i) = 2.0_wp
         ! Upper diagonal entry.
         if (i < test_size) A%data(i, i + 1) = -1.0_wp
         ! Lower diagonal entry.
         if (i < test_size) A%data(i + 1, i) = -1.0_wp
      end do

      ! Initialize Krylov subspaces.
      allocate (U(1:kdim + 1)); allocate (V(1:kdim + 1)); allocate (U0(1))
      call init_rand(U0)
      call initialize_krylov_subspace(U, U0)
      call mat_zero(V)

      ! Initialize internal variables.
      svdvals = 0.0_wp

      ! Compute singular value decomposition.
      call svds(A, U, V, svdvals, residuals, info)

      ! Analytical singular values.
      do i = 1, test_size
         true_svdvals(i) = 2.0_wp*(1.0_wp + cos(i*pi/(test_size + 1)))
      end do

      ! Check convergence.
      call check(error, all_close(svdvals, true_svdvals, rtol, atol))

      return
   end subroutine test_svd_problem

   !------------------------------------
   !-----                          -----
   !-----     CG TEST SUITE        -----
   !-----                          -----
   !------------------------------------

   subroutine collect_cg_testsuite(testsuite)
      ! Collection of tests.
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("CG full computation w. s.p.d. matrix", test_cg_full_computation_spd_matrix) &
                  ]
      return
   end subroutine collect_cg_testsuite

   subroutine test_cg_full_computation_spd_matrix(error)
      ! Error type to be returned.
      type(error_type), allocatable, intent(out) :: error
      ! Linear Problem.
      real(kind=wp) :: Q(test_size, 10_wp*test_size)
      class(spd_matrix), allocatable :: A ! Linear Operator.
      class(rvector), allocatable :: b ! Right-hand side vector.
      class(rvector), allocatable :: x ! Solution vector.
      type(cg_opts)               :: opts
      ! Information flag.
      integer :: info

      ! Initialize linear problem.
      A = spd_matrix(); 
      call random_number(Q); 
      A%data = matmul(Q, transpose(Q)); A%data = 0.5*(A%data + transpose(A%data))
      b = rvector(); call init_rand(b)
      x = rvector(); call x%zero()
      ! CG solver.
      opts = cg_opts(verbose=.false., atol=1e-12_wp, rtol=0.0_wp)
      call cg(A, b, x, info, options=opts)
      ! Check convergence.
      call check(error, norm2(matmul(A%data, x%data) - b%data)**2 < rtol)

      return
   end subroutine test_cg_full_computation_spd_matrix

end module TestIterativeSolvers
