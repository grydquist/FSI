MODULE ELASTMOD
    USE ELEMOD
IMPLICIT NONE

CONTAINS

SUBROUTINE elastmat(kab,f,fun,el,lam,mu)

CLASS(eltype), INTENT(IN) :: el
REAL(KIND = 8), INTENT(IN) :: lam, mu
REAL(KIND = 8), INTENT(OUT) :: kab(el%eNoN*el%dof,el%eNoN*el%dof) &
& ,f(el%eNoN*el%dof)
REAL(KIND=8), ALLOCATABLE, INTENT(IN) :: fun(:,:)
INTEGER :: a,b,g,i,j

kab=0
f = 0
! Loop through element nodes
DO a = 1,el%eNoN

! If this isn't a boundary node, go ahead and calculate all the terms
IF (el%bnd(a,1).eq.0) THEN
!   Loop through the number of dof (this is looping through dof for A)
    DO i = 1,el%dof
!       Loop through neighboring nodes
        DO b = 1,el%eNoN
!           We don't add b to matrix if it is on boundary
            IF(el%bnd(b,1).eq.0) THEN
!               Loop through dof for B nodes
                DO j = 1,el%dof
!                   Calculate the integrals
                    DO g = 1,el%gp
                        kab(2*a+i-2,2*b+j-2) = kab(2*a+i-2,2*b+j-2)&
                        & + 1!Function evaluated at gauss points
                    ENDDO
                ENDDO
            ELSEIF(el%bnd(b,1).eq.1) THEN
!               If b is a dirichlet, it affects the F vector
                DO j = 1,el%dof
                    DO g = 1,el%gp
                        f(2*a+i-2) = f(2*a+i-2) + 1! Function at gauss pts
                    ENDDO
                ENDDO
            ENDIF
        ENDDO
!       Add into the F matrix
        DO g = 1,el%gp
            f(2*a+i-2) = f(2*a+i-2) + 1! Function at gauss pts fun(i,g)
        ENDDO

    ENDDO


ENDIF














END MODULE ELASTMOD