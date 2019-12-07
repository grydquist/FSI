MODULE SOLVEMOD
    IMPLICIT NONE

    CONTAINS

    subroutine inverse(a,c,n)
        !============================================================
        ! Inverse matrix
        ! Method: Based on Doolittle LU factorization for Ax=b
        ! Alex G. December 2009
        !-----------------------------------------------------------
        ! input ...
        ! a(n,n) - array of coefficients for matrix A
        ! n      - dimension
        ! output ...
        ! c(n,n) - inverse matrix of A
        ! comments ...
        ! the original matrix a(n,n) will be destroyed 
        ! during the calculation
        !===========================================================
        implicit none 
        integer n
        double precision a(n,n), c(n,n)
        double precision L(n,n), U(n,n), b(n), d(n), x(n)
        double precision coeff
        integer i, j, k
        
        ! step 0: initialization for matrices L and U and b
        ! Fortran 90/95 aloows such operations on matrices
        L=0.0
        U=0.0
        b=0.0
        
        ! step 1: forward elimination
        do k=1, n-1
           do i=k+1,n
              coeff=a(i,k)/a(k,k)
              L(i,k) = coeff
              do j=k+1,n
                 a(i,j) = a(i,j)-coeff*a(k,j)
              end do
           end do
        end do
        
        ! Step 2: prepare L and U matrices 
        ! L matrix is a matrix of the elimination coefficient
        ! + the diagonal elements are 1.0
        do i=1,n
          L(i,i) = 1.0
        end do
        ! U matrix is the upper triangular part of A
        do j=1,n
          do i=1,j
            U(i,j) = a(i,j)
          end do
        end do
        
        ! Step 3: compute columns of the inverse matrix C
        do k=1,n
          b(k)=1.0
          d(1) = b(1)
        ! Step 3a: Solve Ld=b using the forward substitution
          do i=2,n
            d(i)=b(i)
            do j=1,i-1
              d(i) = d(i) - L(i,j)*d(j)
            end do
          end do
        ! Step 3b: Solve Ux=d using the back substitution
          x(n)=d(n)/U(n,n)
          do i = n-1,1,-1
            x(i) = d(i)
            do j=n,i+1,-1
              x(i)=x(i)-U(i,j)*x(j)
            end do
            x(i) = x(i)/u(i,i)
          end do
        ! Step 3c: fill the solutions x(n) into column k of C
          do i=1,n
            c(i,k) = x(i)
          end do
          b(k)=0.0
        end do
        end subroutine inverse

        !*******************************************************
!*    LU decomposition routines used by test_lu.f90    *
!*                                                     *
!*                 F90 version by J-P Moreau, Paris    *
!*                        (www.jpmoreau.fr)            *
!* --------------------------------------------------- *
!* Reference:                                          *
!*                                                     *
!* "Numerical Recipes By W.H. Press, B. P. Flannery,   *
!*  S.A. Teukolsky and W.T. Vetterling, Cambridge      *
!*  University Press, 1986" [BIBLI 08].                *
!*                                                     * 
!*******************************************************
  
  !  ***************************************************************
  !  * Given an N x N matrix A, this routine replaces it by the LU *
  !  * decomposition of a rowwise permutation of itself. A and N   *
  !  * are input. INDX is an output vector which records the row   *
  !  * permutation effected by the partial pivoting; D is output   *
  !  * as -1 or 1, depending on whether the number of row inter-   *
  !  * changes was even or odd, respectively. This routine is used *
  !  * in combination with LUBKSB to solve linear equations or to  *
  !  * invert a matrix. Return code is 1, if matrix is singular.   *
  !  ***************************************************************
   Subroutine LUDCMP(A,N,INDX,D,CODE)
    INTEGER, INTENT(IN) :: N
   REAL(KIND=8),INTENT(INOUT)::  A(N,N)
   INTEGER, INTENT(OUT) :: INDX(N),D, CODE
   REAL(KIND=8) AMAX,DUM, SUM, VV(100),TINY ! 100 max iter
   INTEGER  NMAX, i,  j, imax,k
   NMAX=100
   TINY=1.5D-16
   D=1; CODE=0
  
   DO I=1,N
     AMAX=0.d0
     DO J=1,N
       IF (DABS(A(I,J)).GT.AMAX) AMAX=DABS(A(I,J))
     END DO ! j loop
     IF(AMAX.LT.TINY) THEN
       CODE = 1
       RETURN
     END IF
     VV(I) = 1.d0 / AMAX
   END DO ! i loop
  
   DO J=1,N
     DO I=1,J-1
       SUM = A(I,J)
       DO K=1,I-1
         SUM = SUM - A(I,K)*A(K,J) 
       END DO ! k loop
       A(I,J) = SUM
     END DO ! i loop
     AMAX = 0.d0
     DO I=J,N
       SUM = A(I,J)
       DO K=1,J-1
         SUM = SUM - A(I,K)*A(K,J) 
       END DO ! k loop
       A(I,J) = SUM
       DUM = VV(I)*DABS(SUM)
       IF(DUM.GE.AMAX) THEN
         IMAX = I
         AMAX = DUM
       END IF
     END DO ! i loop  
     
     IF(J.NE.IMAX) THEN
       DO K=1,N
         DUM = A(IMAX,K)
         A(IMAX,K) = A(J,K)
         A(J,K) = DUM
       END DO ! k loop
       D = -D
       VV(IMAX) = VV(J)
     END IF
  
     INDX(J) = IMAX
     IF(DABS(A(J,J)) < TINY) A(J,J) = TINY
  
     IF(J.NE.N) THEN
       DUM = 1.d0 / A(J,J)
       DO I=J+1,N
         A(I,J) = A(I,J)*DUM
       END DO ! i loop
     END IF 
   END DO ! j loop
  
   RETURN
   END subroutine LUDCMP
  
  
  !  ******************************************************************
  !  * Solves the set of N linear equations A . X = B.  Here A is     *
  !  * input, not as the matrix A but rather as its LU decomposition, *
  !  * determined by the routine LUDCMP. INDX is input as the permuta-*
  !  * tion vector returned by LUDCMP. B is input as the right-hand   *
  !  * side vector B, and returns with the solution vector X. A, N and*
  !  * INDX are not modified by this routine and can be used for suc- *
  !  * cessive calls with different right-hand sides. This routine is *
  !  * also efficient for plain matrix inversion.                     *
  !  ******************************************************************
   Subroutine LUBKSB(A,N,INDX,B)
    INTEGER, INTENT(IN):: N, INDX(N)
   REAL(KIND=8), INTENT(IN) :: A(N,N) 
   REAL(KIND=8), INTENT(INOUT) :: B(N)
   REAL(KIND=8)  SUM, TINY
   INTEGER i,  j, ii, ll
  
   TINY=1.5D-16
   II = 0
  
   DO I=1,N
     LL = INDX(I)
     SUM = B(LL)
     B(LL) = B(I)
     IF(II.NE.0) THEN
       DO J=II,I-1
         SUM = SUM - A(I,J)*B(J)
       END DO ! j loop
     ELSE IF(abs(SUM).GT.TINY) THEN
       II = I
     END IF
     B(I) = SUM
   END DO ! i loop
  
   DO I=N,1,-1
     SUM = B(I)
     IF(I < N) THEN
       DO J=I+1,N
         SUM = SUM - A(I,J)*B(J)
       END DO ! j loop
     END IF
     B(I) = SUM / A(I,I)
   END DO ! i loop
  
   RETURN
   END subroutine LUBKSB

   subroutine seidel(crit,n,mat,b,omega,x,residu,iter,rc)
    INTEGER ITERMAX            ! Maximal number of iterations
    REAL(KIND=8) :: ONE,TWO,ZERO, TINY
      integer,INTENT(IN) :: crit, n
      INTEGER, INTENT(OUT) :: iter, rc
      REAL(KIND=8), INTENT(IN) :: mat(n,n),b(n),omega
      REAL(KIND=8), INTENT(OUT) :: residu(n)
      REAL(KIND=8), INTENT(INOUT):: x(n)
      INTEGER :: i,j
    !*====================================================================*
    !*                                                                    *
    !*  seidel solves the linear system  mat * x = b  iteratively.        *
    !*  Here  mat  is a nonsingular  n x n  matrix, b is the right hand   *
    !*  side for the linear system and x is the solution.                 *
    !*                                                                    *
    !*  seidel uses the Gauss Seidel Method with relaxation for a given   *
    !*  relaxation coefficient 0 < omega < 2.                             *
    !*  If  omega = 1, the standard Gauss Seidel method (without          *
    !*  relaxation) is performed.                                         *
    !*                                                                    *
    !*====================================================================*
    !*                                                                    *
    !*   Applications:                                                    *
    !*   =============                                                    *
    !*      Solve linear systems with nonsingular system matrices that    *
    !*      satisfy one of the following criteria: row sum criterion,     *
    !*      column sum criterion or the criterion of Schmidt and v. Mises.*
    !*      Only if at least one of these criteria is satisfied for mat,  *
    !*      convergence of the scheme is guaranteed [See BIBLI 11].       *
    !*                                                                    *
    !*====================================================================*
    !*                                                                    *
    !*   Input parameters:                                                *
    !*   ================                                                 *
    !*      crit     integer crit                                         *
    !*               select criterion                                     *
    !*               =1 : row sum criterion                               *
    !*               =2 : column sum criterion                            *
    !*               =3 : criterion of Schmidt-v.Mises                    *
    !*               other : no check                                     *
    !*      n        integer n ( n > 0 )                                  *
    !*               size of mat, b and x                                 *
    !*      mat      REAL*8   mat(n,n)                                    *
    !*               Matrix of the liear system                           *
    !*      b        REAL*8 b(n)                                          *
    !*               Right hand side                                      *
    !*      omega    REAL*8 omega; (0 < omega < 2)                        *
    !*               Relaxation coefficient.                              *
    !*      x        REAL*8  x(n)                                         *
    !*               Starting vector for iteration                        *
    !*                                                                    *
    !*   Output parameters:                                               *
    !*   ==================                                               *
    !*      x        REAL*8  x(n)                                         *
    !*               solution vector                                      *
    !*      residu   REAL*8   residu(n)                                   *
    !*               residual vector  b - mat * x; close to zero vector   *
    !*      iter     integer iter                                         *
    !*               Number of iterations performed                       *
    !*      rc       integer return code                                  *
    !*               =  0     solution has been found                     *
    !*               =  1     n < 1  or omega <= 0 or omega >= 2          *
    !*               =  2     improper mat or b or x (not used here)      *
    !*               =  3     one diagonal element of mat vanishes        *
    !*               =  4     Iteration number exceeded                   *
    !*               = 11     column sum criterion violated               *
    !*               = 12     row sum criterion violated                  *
    !*               = 13     Schmidt-v.Mises criterion violated          *
    !*                                                                    *
    !*====================================================================*
      REAL(KIND=8) tmp, eps, matt(n,n), bt(n);
      TINY = 1.d-10
      ONE=1.d0
      TWO=2.d0
      ZERO=0.d0
      ITERMAX=500
    
      iter = 0                        !Initialize iteration counter
      rc = 0
      matt = mat
      bt = b
    
      if (n<1.or.omega<=ZERO.or.omega>=TWO) then
        rc=1
        return
      end if
    
      eps = 1.d-10
    
      do i=1, n                       

        if (abs(matt(i,i)) .lt. TINY) then
          rc=3
          return
        end if
        tmp = ONE / matt(i,i)
        do j=1, n
          matt(i,j)= matt(i,j)*tmp
        end do
        bt(i) = bt(i)*tmp               !adjust right hand side bt
      
      end do
    
    
      !check convergence criteria
      if (crit==1) then
         do i = 1, n                  !row sum criterion
           tmp=ZERO
           do j=1,n
             tmp = tmp + dabs(matt(i,j))
           end do
           if (tmp >= TWO) then
             rc=11
             return
           end if 
         end do
      else if (crit==2) then  
         do j=1, n                    !column sum criterion
         tmp=ZERO
           do i=1,n
             tmp = tmp + dabs(matt(i,j))
           end do
           if (tmp >= TWO) then
             rc=12
             return
           end if
         end do
      else if (crit==3) then
         tmp=ZERO
         do i=1, n
           do j=1, n                  !criterion of Schmidt
             tmp = tmp + matt(i,j)**2  !von Mises
           end do
         end do
         tmp = DSQRT(tmp - ONE)
         if (tmp >= ONE) then
           rc=13
           return
         end if
      end if
    
      do i=1, n 
        residu(i) = x(i)              !store x in residu
      end do
    
      do while (iter <= ITERMAX)      !Begin iteration
      
        iter=iter+1
    
        do i=1, n
          tmp=bt(i)
          do j=1, n
            tmp =  tmp - matt(i,j) * residu(j)
          end do 
          residu(i) = residu(i) + omega * tmp
        end do
    
        do i=1, n                     !check break-off criterion
          tmp = x(i) - residu(i)
          if (DABS (tmp) <= eps) then
            x(i) = residu(i)          !If rc = 0 at end of loop
            rc = 0                    !  -> stop iteration
          else
            do j=1, n 
              x(j) = residu(j)
            end do
            rc = 4
            goto 10
          end if
        end do
        if (rc == 0) goto 20          !solution found
    10 end do                         !End iteration
    
    20 do i=1, n                      !find residual vector
         tmp=bt(i)
         do j=1, n
           tmp = tmp - matt(i,j) * x(j)
         end do
         residu(i) = tmp
       end do
    
      return
    
    end


  !===========================================
    SUBROUTINE SOR(n,a,b,x)
      IMPLICIT NONE
  
      INTEGER,INTENT(IN)::n
      INTEGER::i,j,k,maxn
  
      REAL(KIND=8):: w,s,tol,norm
      REAL(KIND=8),DIMENSION(n,n), INTENT(IN) :: a
      REAL(KIND=8),DIMENSION(n), INTENT(IN) :: b
      REAL(KIND=8),DIMENSION(n), INTENT(INOUT) :: x
      REAL(KIND=8),DIMENSION(n) :: xn

      x=0
      maxn=100
      tol=0.00001
      w=1.4

      DO k=1,maxn
          DO i=1,n
              s=0
              DO j=1,n
                  IF (i<j) THEN
                      s=s+a(i,j)*x(j)
                  ELSE IF (i>j) THEN
                      s=s+a(i,j)*xn(j)
                  END IF
              END DO
              xn(i)=(1-w)*x(i)+w*(b(i)-s)/a(i,i)
          END DO
  
          norm=MAXVAL(ABS(x-xn))
          
          IF(norm<tol)EXIT
          x=xn
      END DO
      if (k.eq.maxn+1) print *, "oh no"
  
    END SUBROUTINE SOR

    function dnrm22 ( n, x, incx )

      !*****************************************************************************80
      !
      !! DNRM2 returns the euclidean norm of a vector.
      !
      !  Discussion:
      !
      !    This routine uses double precision real arithmetic.
      !
      !     DNRM2 ( X ) = sqrt ( X' * X )
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license. 
      !
      !  Modified:
      !
      !    16 May 2005
      !
      !  Author:
      !
      !    Original FORTRAN77 version by Charles Lawson, Richard Hanson, 
      !    David Kincaid, Fred Krogh.
      !    FORTRAN90 version by John Burkardt.
      !
      !  Reference:
      !
      !    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
      !    LINPACK User's Guide,
      !    SIAM, 1979,
      !    ISBN13: 978-0-898711-72-1,
      !    LC: QA214.L56.
      !
      !    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
      !    Algorithm 539, 
      !    Basic Linear Algebra Subprograms for Fortran Usage,
      !    ACM Transactions on Mathematical Software,
      !    Volume 5, Number 3, September 1979, pages 308-323.
      !
      !  Parameters:
      !
      !    Input, integer ( kind = 4 ) N, the number of entries in the vector.
      !
      !    Input, real ( kind = 8 ) X(*), the vector whose norm is to be computed.
      !
      !    Input, integer ( kind = 4 ) INCX, the increment between successive 
      !    entries of X.
      !
      !    Output, real ( kind = 8 ) DNRM2, the Euclidean norm of X.
      !
        implicit none
      
        real ( kind = 8 ) absxi
        real ( kind = 8 ) dnrm22
        integer ( kind = 4 ) incx
        integer ( kind = 4 ) ix
        integer ( kind = 4 ) n
        real ( kind = 8 ) norm
        real ( kind = 8 ) scale
        real ( kind = 8 ) ssq
        real ( kind = 8 ) x(*)
      
        if ( n < 1 .or. incx < 1 ) then
      
          norm  = 0.0D+00
      
        else if ( n == 1 ) then
      
          norm  = abs ( x(1) )
      
        else
      
          scale = 0.0D+00
          ssq = 1.0D+00
      
          do ix = 1, 1 + ( n - 1 ) * incx, incx
            if ( abs(x(ix)) .gt. 1d-15 ) then
              absxi = abs ( x(ix) )
              if ( scale < absxi ) then
                ssq = 1.0D+00 + ssq * ( scale / absxi )**2
                scale = absxi
              else
                ssq = ssq + ( absxi / scale )**2
              end if
            end if
          end do
          norm  = scale * sqrt ( ssq )
        end if

        norm = 0
        ssq = 0
        DO ix = 1,n
          ssq = ssq + x(ix)**2
        ENDDO
        norm =sqrt(ssq)
      
        dnrm22 = norm
      
        return
      end

END MODULE SOLVEMOD