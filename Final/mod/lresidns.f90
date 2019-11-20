MODULE LRESIDNSMOD
    USE SOLMOD
    USE ELEMOD
    USE TINTMOD
IMPLICIT NONE

CONTAINS

SUBROUTINE lresidns(G,Gt,fun,el,u,p)
TYPE(eltype), INTENT(IN) :: el
REAL(KIND=8), ALLOCATABLE, INTENT(IN) :: fun(:,:)
TYPE(soltype), INTENT(IN) :: u,p
REAL(KIND=8), INTENT(OUT) :: G(el%eNoN*el%dof),Gt(el%eNoN*el%dof,el%dof)
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

!               If i = 1, we're dealing with x-mom, 2->y-mom, 3->cont
                IF((i.eq.1).or.(i.eq.2)) THEN
!                   Calculate Rm in direction i for node a (G(1:2))
                    DO gp = 1,el%gp
                        G(ai) = G(ai) &
                        & !   + el%Ng(a,gp)*el%Ng(b,gp)*u%ddalphm(j,el%nds(b)) &
                    !    &    - u%dalphf(j,el%nds(b))*(el%Nxg(a,1,gp)*el%Nxg(b,1,gp) &
                    !    &    + el%Nxg(a,2,gp)*el%Nxg(b,2,gp))*el%Wg(gp)
                    ENDDO
                ELSE
!                   Calulate Rc for node a (G(3))
                    DO gp = 1,el%gp

                    ENDDO
                ENDIF
            ENDDO
        ENDDO
!       Calculate tangent matrix
        DO j = 1,el%dof
            DO gp = 1,el%gp
                Gt(ai,j)= Gt(ai,j) &
                &    + alphm*(el%Ng(a,gp))*el%Wg(gp)
            ENDDO
        ENDDO
    ENDIF
ENDDO
ENDDO


END SUBROUTINE lresidns

END MODULE LRESIDNSMOD