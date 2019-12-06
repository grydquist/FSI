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

END MODULE SOLVEMOD