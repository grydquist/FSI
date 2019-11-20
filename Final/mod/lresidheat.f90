MODULE LRESIDHEATMOD
    USE SOLMOD
    USE ELEMOD
    USE TINTMOD
IMPLICIT NONE

CONTAINS

SUBROUTINE lresidh(G,Gt,fun,el,sol)
TYPE(eltype), INTENT(IN) :: el
REAL(KIND=8), ALLOCATABLE, INTENT(IN) :: fun(:,:)
TYPE(soltype), INTENT(IN) :: sol
REAL(KIND=8), INTENT(OUT) :: G(el%eNoN),Gt(el%eNoN)
INTEGER :: a,b,gp,i,j,ai,bi

Gt = fun(1,1)
Gt = 0
G = 0
ai = 0

! Loop through element nodes
DO a = 1,el%eNoN
! Loop through dof (for index a)
DO i = 1,el%dof
    bi = 0
    ai = ai+1
!   If this isn't a boundary node, go ahead and calculate the resid at a
    IF (el%bnd(a,1,i).eq.0) THEN
!       Go through surrounding nodes and incorporate into resid for a 
        DO b = 1,el%eNoN
!           Again through dof, now for b
            DO j =1,el%dof
                bi = bi+1
!               Loop through Gauss points to calc integral
                DO gp = 1,el%gp
                    G(a) = G(a) &
                    &    + el%Ng(a,gp)*el%Ng(b,gp)*sol%ddalphm(j,el%nds(b)) &
                    &    - sol%dalphf(j,el%nds(b))*(el%Nxg(a,1,gp)*el%Nxg(b,1,gp) &
                    &    + el%Nxg(a,2,gp)*el%Nxg(b,2,gp))*el%Wg(gp)
                ENDDO
            ENDDO
        ENDDO
!       Calculate tangent matrix
        DO gp = 1,el%gp
            Gt(a)= Gt(a) &
            &    + alphm*(el%Ng(a,gp))*el%Wg(gp)
        ENDDO
    ENDIF
ENDDO
ENDDO





END SUBROUTINE lresidh

END MODULE LRESIDHEATMOD