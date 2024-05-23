module LightKrylov_expmlib
   !! This module provides the routines for computation of dense matrix exponentiation and default routines for 
   !! the approximation of the action of the exponential propagator on a given input vector using Krylov 
   !!
   !! The module defines the procedures:
   !!
   !! - `expm` : Compute the dense matrix exponential of the square matrix \(\mathbf{A}\).
   !! - `kexpm` : Compute an approximation of the action of the exponential propagator \( e^{ \tau \mathbf{A} } \)
   !!             on a given input vector \( \mathbf{b} \).
   !!
   !! It also provides an abstract interface for the exponential propagator to which user-defined alternatives 
   !! to conform.
   Use LightKrylov
   Use lightkrylov_BaseKrylov
   Use LightKrylov_Utils

   ! Fortran standard library.
   use stdlib_optval, only: optval

   implicit none

   private
   ! Matrix operations for abstract vector types
   public :: kexpm, expm
   ! Action of the exponential propagator on a vector
   public :: abstract_exptA, k_exptA

   type, extends(abstract_opts), public :: kexpm_opts
      !! Extended `abstract_opts` type to pass options to `kexpm`
      real(kind=wp) :: tol = atol
      !! Absolute tolerance (default: `epsilon(1.0_wp)`).
      logical :: verbose = .false.
      !! Verbosity control (default: `.false.`).
      integer :: kdim = 30
      !! Dimension of the Krylov subspace (default: `kdim = 30`).
   end type kexpm_opts

   abstract interface
      subroutine abstract_exptA(vec_out, A, vec_in, tau, info, trans)
         import abstract_linop, abstract_vector, abstract_opts, wp
         !! Abstract interface to use a user-defined exponential propagator  in `LightKrylov`
         class(abstract_vector),         intent(out)   :: vec_out
         !! Solution vector
         class(abstract_linop),          intent(inout) :: A
         !! Linear operator
         class(abstract_vector),         intent(in)    :: vec_in
         !! Input vector
         real(kind=wp),                  intent(in)    :: tau
         !! Time horizon for integration
         integer,                        intent(out)   :: info
         !! Information flag
         logical, optional,              intent(in)    :: trans
         !! Use transpose?
      end subroutine abstract_exptA
   end interface

   interface kexpm
      module procedure kexpm_vec
      module procedure kexpm_mat
   end interface

contains

   subroutine k_exptA(vec_out, A, vec_in, tau, info, trans)
      
      class(abstract_vector),         intent(out)   :: vec_out
      !! Solution vector
      class(abstract_linop),          intent(inout) :: A
      !! Linear problem
      class(abstract_vector),         intent(in)    :: vec_in
      !! Input vector
      real(kind=wp),                  intent(in)    :: tau
      !! Time horizon for integration
      integer,                        intent(out)   :: info
      !! Information flag
      logical, optional,              intent(in)    :: trans
      !! Use transpose?

      ! Internals
      real(kind=wp)                                 :: tol
      integer                                       :: kdim
      logical                                       :: verbose
      logical                                       :: transpose

      tol     = atol
      verbose = .false.
      kdim    = 200
      transpose = optval(trans, .false.)

      call kexpm(vec_out, A, vec_in, tau, tol, info, trans=transpose, verbosity = verbose, kdim =kdim)

   end subroutine k_exptA

   subroutine kexpm_mat(C, A, B, tau, tol, info, trans, verbosity, kdim)
      !! Approximates the action of the exponential propagator (matrix exponential) of a linear 
      !! operator \( \mathbf{A} \) on a given matrix \( \mathbf{B} \) by computing the action of the exponential propagator
      !! on the projection of the operator on a (small) Krylov subspace.
      !! 
      !! Matrix version of the subroutine `kexpm_vec`.
      class(abstract_vector), intent(out)   :: C(:)
      !! Best approximation of the value of \( e^{\tau \mathbf{A} } \mathbf{B} \) in the computed Krylov subspace
      class(abstract_linop),  intent(inout) :: A
      !! Linear operator to be exponentiated
      class(abstract_vector), intent(in)    :: B(:)
      !! Input matrix on which to apply \( e^{\tau \mathbf{A} } \)
      real(kind=wp),          intent(in)    :: tau
      !! Time horizon for exponentiation
      real(kind=wp),          intent(in)    :: tol
      !! Solution tolerance based on error estimates
      integer,                intent(out)   :: info
      !! Information flag
      logical, optional,      intent(in)    :: trans
      !! Use transpose?
      logical, optional,      intent(in)    :: verbosity
      !! Verbosity control (default `.false.`)
      integer, optional,      intent(in)    :: kdim
      !! Maximum size of the Krylov subspace (default `kdim = 100`)
      
      ! Internal variables
      integer, parameter :: kmax = 100
      integer :: i, k, p, kpm, kp, kpp, nk
      ! Arnoldi factorisation
      class(abstract_vector), allocatable   :: X(:)
      real(kind=wp), allocatable            :: H(:,:)
      ! Normalisation & temp arrays
      real(kind=wp), allocatable            :: R(:,:), E(:,:), em(:,:)
      integer, allocatable                  :: perm(:), ptrans(:)
      class(abstract_vector), allocatable   :: Xwrk(:), Cwrk(:)
      real(kind=wp) :: err_est
      ! Optional arguments
      logical                               :: transpose
      logical                               :: verbose
      integer                               :: nstep

      ! Determine block size
      p = size(B)

      ! Deal with optional arguments
      transpose = optval(trans, .false.)
      verbose   = optval(verbosity, .false.)
      nstep     = optval(kdim, kmax)
      nk        = nstep*p
      info = 0

      ! Allocate memory
      allocate(R(1:p,1:p))
      allocate(perm(1:p))
      allocate(ptrans(1:p))
      allocate(X(1:p*(nk+1)), source=B(1)); 
      allocate(H(1:p*(nk+1),1:p*(nk+1)))
      allocate(E(1:p*(nk+1),1:nk))
      allocate(em(1:p,1:p))

      ! Scratch arrays
      allocate(Xwrk(1:p), source=B(1)); allocate(Cwrk(1:p), source=B(1))

      ! Normalize input vector and set initialise Krylov subspace
      R = 0.0_wp
      call mat_zero(Xwrk)
      call mat_copy(Xwrk(1:p), B(1:p))
      call qr_factorization(Xwrk(1:p), R(1:p,1:p), perm(1:p), info, ifpivot=.true.)
      call apply_permutation(R, perm, trans = .true.)
      !R = matmul(R, transpose(perm))
      if (norm2(R(1:p,1:p)) .eq. 0.0_wp) then
         ! input is zero => output is zero
         call mat_zero(C)
         err_est = 0.0_wp
         k = 0
         kpp = p
      else
         call initialize_krylov_subspace(X,Xwrk(1:p))
         H = 0.0_wp
         
         expm_arnoldi: do k = 1, nk 
            kpm = (k-1)*p
            kp  = kpm + p
            kpp = kp  + p

            ! reset wrk arrays
            E = 0.0_wp; call mat_zero(Xwrk)

            ! compute kth stop of the Arnoldi factorization
            call arnoldi_factorization(A, X(1:kpp), H(1:kpp,1:kp), info, kstart=k, kend=k, transpose=transpose, block_size=p)
            
            ! compute approximation
            if (info .eq. kp) then
               ! Arnoldi breakdown, do not consider extended matrix
               kpp = kp
            endif
            ! compute the (dense) matrix exponential of the extended Hessenberg matrix
            call expm(E(1:kpp,1:kpp), tau*H(1:kpp,1:kpp))
            ! project back into original space
            call mat_zero(C)
            call mat_mult(Xwrk(1:p), X(1:kpp), E(1:kpp,1:p))
            call mat_mult(C(1:p), Xwrk(1:p), R(1:p,1:p))

            ! cheap error estimate (this is actually the magnitude of the included correction thus too conservative)
            if (info .eq. kp) then
               ! approximation is exact
               err_est = 0.0_wp
            else
               em = matmul(E(kp+1:kpp,1:p), R(1:p,1:p))
               err_est = norm2(em)
            endif
            if (err_est .lt. tol) then
               ! exit the Arnoldi iteration.
               exit expm_arnoldi
            endif
         end do expm_arnoldi
      endif
   
      if (err_est .le. tol) then
         if (verbose) then
            if (p.eq.1) then
               write(*, *) 'Arnoldi approxmation of the action of the exp. propagator converged'
            else
               write(*, *) 'Block Arnoldi approxmation of the action of the exp. propagator converged'
            endif 
            write(*, *) '    n° of vectors:', k+1, 'per input vector, total:', kpp
            write(*, *) '    estimated error:   ||err_est||_2 = ', err_est
            write(*, *) '    desired tolerance:           tol = ', tol
            write(*, *)
         endif
      else
         info = -1
         if (verbose) then
            if (p.eq.1) then
               write(*, *) 'Arnoldi approxmation of the action of the exp. propagator did not converge'
            else
               write(*, *) 'Block Arnoldi approxmation of the action of the exp. propagator did not converge'
            endif
            write(*, *) '    maximum n° of vectors reached: ', nstep+1,&
                        & 'per input vector, total:', kpp
            write(*, *) '    estimated error:   ||err_est||_2 = ', err_est
            write(*, *) '    desired tolerance:           tol = ', tol
            write(*, *)
         endif
      endif
      
      return
   end subroutine kexpm_mat

   subroutine kexpm_vec(c, A, b, tau, tol, info, trans, verbosity, kdim)
      !! Approximates the action of the exponential propagator (matrix exponential) of a linear 
      !! operator \( \mathbf{A} \) on a given vector \( \mathbf{b} \) by computing the action of 
      !! the exponential propagator on the projection of the operator on a (small) Krylov subspace.
      !!
      !! **Mathematical Formulation**
      !! 
      !! Given a linear operator \( \mathbf{A} \in \mathbb{R}^{n \times n} \) and an input matrix 
      !! \( \mathbf{b}  \in \mathbb{R}^{n} \) we compute:
      !! \[
      !!     e^{\tau \mathbf{A}} \mathbf{b} \approx \mathbf{X}_{m+1} e^{\tau \mathbf{\bar{H}}_m} \mathbf{e}_1
      !! \]
      !!
      !! where
      !!
      !! - \( \mathbf{X}_{m+1} = [ \mathbf{x}_{1}, \ldots, \mathbf{x}_{m+1} ] \in \mathbb{R}^{n \times m+1 }\)
      !!       is the extended basis generated by the m-step Arnoldi factorization 
      !! - \( \mathbf{\bar{H}}_m \in \mathbb{R}^{m+1 \times m+1} \) is an augmented Hessenberg matrix with
      !! \[
      !!     \mathbf{\bar{H}}_m = \begin{bmatrix}
      !!                               \mathbf{H}_m         & \mathbf{0} \\
      !!                               c \phi(\mathbf{H}_m) & 1
      !!                          \end{bmatrix}
      !! \]
      !! where \( \mathbf{H}_m \in \mathbb{R}^{m \times m} \) is the regular Hessenberg matrix generated 
      !! by the Arnoldi factorization, \( c \) is a row vector defined as
      !! \[
      !!    c = h_{m+1,m} \textbf{e}^T_m \in \mathbb{R}^{1 \times m}
      !! \]
      !! involving the norm of the residual of the Arnoldi factorization \( h_{m+1,m} \) and the m-th 
      !! unit vector \( \textbf{e}^T_m \), the \( \varphi \)-function is defined as
      !! \[
      !!    \varphi(\mathbf{z}) = \frac{e^{\mathbf{z}} - 1}{z}
      !! \]
      !! and \( \mathbf{0} \in \mathbb{R}^{m} \) is a zero vector.
      !! 
      !! **Algorithmic Features**
      !! 
      !! - Very good approximations of the action of the exponential propagator are obtained 
      !!   with very few arnoldi steps. \( \mathcal{O}(10) \) steps are sufficient for working precision for 
      !!   "small" matrices \( \mathbf{A} \) (or \( \tau \ll 1 \) )
      !! - The small dense problem of the exponential of the Hessenberg matrix is solved using
      !!   the scaling and squaring approach combined with a rational Pade approximation.
      !! - The correction to the naive Krylov approach proposed by Gallpoulos & Saad making
      !!   use of the additional information in the last row of the Hessenberg matrix implemented
      !!
      !! **Advantages**
      !!
      !! - Very accurate for "small" matrices (in terms of their norm), i.e. for small \( \tau \).
      !! - A fairly accuracte error estimate is computed based on the Hessenberg matrix to
      !!   terminate iteration when needed (see Er3, p 223, Y. Saad, 1992)
      !! - Block arnoldi method allows for the right hand side to contain more than 1 vector
      !!
      !! **Limitations**
      !! 
      !! - No general computation of \( e^{ \tau \mathbf{A} } \) but action of the propagator on a specific 
      !!   \( \mathbf{b} \). The process must be restarted for each new \( \mathbf{b} \).
      !! - The error estimation is quite crude and conservative. In practise, the approximation
      !!   is better than the estimator suggests
      !!
      !! **References**
      !!
      !!  - Y. Saad, "Analysis of Some Krylov Subspace Approximations to the Matrix Exponent", 
      !!    SIAM Journal on Numerical Analysis, Volume 29, Number 1, Feb. 1992, pages 209-228.
      !!  - E. Gallopoulos & Y. Saad, "Efficient soltuion of parabolic equations by Krylov
      !!    approximation methods", SIAM Journal on Scientific and Statistical Computing, 
      !!    Volume 13, Number 5, September 1992, pages 1236-1264.
      !!  - C. Moler & C. VanLoan, "Nineteen Dubious Ways to Compute the Exponential of a 
      !!    Matrix, Twenty-Five Years Later", SIAM Review, Volume 45, Number 1, March 2003, 
      !!    pages 3-49.
      class(abstract_vector), intent(out)   :: c
      !! Best approximation of the value of \( e^{\tau \mathbf{A} } \mathbf{b} \) in the computed Krylov subspace
      class(abstract_linop),  intent(inout) :: A
      !! Linear operator to be exponentiated
      class(abstract_vector), intent(in)    :: b
      !! Input vector on which to apply \( \mbox{exp}(\tau \mathbf{A}) \)
      real(kind=wp),          intent(in)    :: tau
      !! Time horizon for exponentiation
      real(kind=wp),          intent(in)    :: tol
      !! Solution tolerance based on error estimates
      integer,                intent(out)   :: info
      !! Information flag
      logical, optional,      intent(in)    :: trans
      !! Use transpose?
      logical, optional,      intent(in)    :: verbosity
      !! Verbosity control (default `.false.`)
      integer, optional,      intent(in)    :: kdim
      !! Maximum size of the Krylov subspace (default `kdim = 100`)
         
      ! Internal variables
      integer, parameter :: kmax = 100
      integer :: i, k, p, km, kp, nk
      ! Arnoldi factorisation
      class(abstract_vector), allocatable   :: X(:)
      real(kind=wp), allocatable            :: H(:,:)
      ! Normalisation & temp arrays
      real(kind=wp), allocatable            :: E(:,:)
      class(abstract_vector), allocatable   :: xwrk
      real(kind=wp)                         :: err_est, beta
      ! Optinal arguments
      logical                               :: transpose
      logical                               :: verbose
      integer                               :: nstep
      
      ! Deal with optional arguemnts
      transpose = optval(trans, .false.)
      verbose   = optval(verbosity, .false.)
      nstep     = optval(kdim, kmax)
      nk        = nstep
      
      info = 0
      
      ! allocate memory
      allocate(X(1:nk+1), source=b); allocate(H(1:nk+1,1:nk+1))
      allocate(E(1:nk+1,1:nk+1))
      ! scratch arrays
      allocate(xwrk, source=b)
      
      ! normalize input vector and set initialise Krylov subspace
      beta = b%norm()
      if ( beta .eq. 0.0_wp ) then
         ! input is zero => output is zero
         call c%zero()
         err_est = 0.0_wp
         kp = 1
      else
         call mat_zero(X)
         call X(1)%axpby(0.0_wp, b, 1.0_wp)
         call X(1)%scal(1.0/beta)
         H = 0.0_wp

         expm_arnoldi: do k = 1, nk 
            km = k - 1
            kp = k + 1

            ! reset wrk arrays
            E = 0.0_wp

            ! compute kth stop of the Arnoldi factorization
            call arnoldi_factorization(A, X(1:kp), H(1:kp,1:k), info, kstart=k, kend=k, transpose=transpose)

            ! compute approximation
            if (info .eq. k) then 
               ! Arnoldi breakdown, do not consider extended matrix
               kp = k
            endif
            ! compute the (dense) matrix exponential of the extended Hessenberg matrix
            call expm(E(1:kp,1:kp), tau*H(1:kp,1:kp))
            ! project back into original space
            call get_vec(xwrk, X(1:kp), E(1:kp,1))
            call c%axpby(0.0_wp, xwrk, beta)

            ! cheap error estimate (this is actually the magnitude of the included correction thus too conservative)
            if (info .eq. k) then 
               ! approximation is exact
               err_est = 0.0_wp
            else
               err_est = abs(E(kp,1) * beta)
            endif
            if (err_est .lt. tol) then
               ! exit the Arnoldi iteration.
               exit expm_arnoldi
            endif
         end do expm_arnoldi
      endif

      if (err_est .le. tol) then
         if (verbose) then
            write(*, *) 'Arnoldi-based approxmation of the exp. propagator converged'
            write(*, *) '    n° of vectors:', kp
            write(*, *) '    estimated error:   ||err_est||_2 = ', err_est
            write(*, *) '    desired tolerance:           tol = ', tol
            write(*, *)
         endif
      else
         info = -1
         if (verbose) then
            write(*, *) 'Arnoldi-based approxmation of the exp. propagator did not converge'
            write(*, *) '    maximum n° of vectors reached: ', nk + 1
            write(*, *) '    estimated error:   ||err_est||_2 = ', err_est
            write(*, *) '    desired tolerance:           tol = ', tol
            write(*, *)
         endif
      endif
      
      return
   end subroutine kexpm_vec
   
   subroutine expm(E,A, order)
      !! Computes the matrix exponential \( \mathbf{E} = e^{\mathbf{A}} \) for a real-valued dense square matrix of 
      !! order \( n \) using the scaling and squaring approach, where the scaled problem is computed
      !! using a rational matrix Pade approximation of the form
      !! \[ 
      !!    \mathbf{R}_{pq}(\mathbf{X}) = [ \mathbf{Q}_q(\mathbf{X}) ]^{-1} \mathbf{P}_p(\mathbf{X}) 
      !! \]
      !!
      !! where \( p \) and \( q \) are the respective orders of the matrix polynomials \( \mathbf{P}(\mathbf{X}) \) and \( \mathbf{Q}(\mathbf{X}) \).
      !!
      !! **Algorithmic Features**
      !! 
      !! - Only diagonal approximations (\( p = q \) ) are used as these are more efficient to compute
      !!   than off-diagonal approximations for the same approximation order
      !! - 12th order diagonal Pade approximation by default
      !!
      !! **Advantages**
      !!
      !! - Very accurate for "small" matrices (in terms of their norm). When considering the
      !!   temporal evolution of an LTI system, i.e. \( e^{ \tau \mathbf{A} } \), the approximation is very good
      !!   for small \( t \).
      !! - Avoids numerical problems linked to the "hump" appearing for non-normal matrices \( \mathbf{A} \)
      !!
      !! **Limitations**
      !!
      !! - No error estimation/check based on the approximation order and the matrix norm
      !!
      !!  `expm` is an adaptation of the function `R8mat_expm1` by C. Moler & C. van Loan.
      !!
      !!  Original licensing:
      !!
      !!    This code is distributed under the MIT license.
      !!
      !! **References**
      !!
      !!  - C. Moler & C. Van Loan, "Nineteen Dubious Ways to Compute the Exponential of a 
      !!    Matrix, Twenty-Five Years Later", SIAM Review, Volume 45, Number 1, March 2003, 
      !!    pages 3-49.
      real(kind=wp), intent(in)      :: A(:,:)
      !! Square matrix to be exponentiated
      real(kind=wp), intent(out)     :: E(:,:)
      !! Output
      integer, optional, intent(in)  :: order
      !! Order of the Pade approximation
      
      ! Internal variables
      real(kind=wp), allocatable :: A2(:,:), Q(:,:), X(:,:), invQ(:,:), wrk(:)
      real(kind=wp)              :: a_norm, c
      integer                    :: n, ee, k, s
      logical                    :: p
      integer                    :: p_order

      ! Deal with optional arguemnts
      p_order = optval(order, 10) 

      n = size(A,1)
      ! allocate memory
      allocate(A2(1:n,1:n)); allocate(X(1:n,1:n))
      allocate(Q(1:n,1:n));  allocate(invQ(1:n,1:n))
      allocate(wrk(1:n)) 

      ! Compute the L-infinity norm.
      a_norm = norml(A)
      ! Determine scaling factor for the matrix.
      ee = int(log2(a_norm)) + 1
      s  = max(0, ee + 1)
      ! Scale the input matrix & initialize polynomial
      A2 = A / 2.0_wp**s
      X = A2 

      ! Initialize P & Q and add first step
      c = 0.5_wp
      E = 0.0_wp; forall (k=1:n) E(k, k) = 1.0_wp
      E = E + c*A2
      
      Q = 0.0_wp; forall (k=1:n) Q(k, k) = 1.0_wp
      Q = Q - c*A2
      
      ! Iteratively compute Pade approximation
      p = .true.
      do k = 2, p_order
         c = c*( p_order - k + 1 )/( k*( 2*p_order - k + 1 ) )
         X = matmul( A2, X )
        E = E + c*X
        if ( p ) then
          Q = Q + c*X
        else
          Q = Q - c*X
        end if
        p = .not. p
      end do

      !  E -> Q^(-1) * P
      invQ = Q
      call inv(invQ)
      E = matmul(invQ,E)

      !  E -> E^(2*S)
      do k = 1, s
        E = matmul ( E, E )
      end do
     
      return
   end subroutine expm

end module LightKrylov_expmlib