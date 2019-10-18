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
INTEGER :: a,b,g,i,j,ai,bi

kab=0
f = 0
ai = 0
! Loop through element nodes
DO a = 1,el%eNoN
!   Loop through the number of dof (this is looping through dof for A)
DO i = 1,el%dof
    bi = 0
    ai = ai + 1
!   If this isn't a boundary node, go ahead and calculate all the terms
    IF (el%bnd(a,1,i).eq.0) THEN
!       Loop through neighboring nodes
        DO b = 1,el%eNoN
!           Loop through dof for B nodes
            DO j = 1,el%dof
                bi = bi+1
!               We don't add b to matrix if it is on boundary
                IF(el%bnd(b,1,j).eq.0) THEN
!                   Calculate the integrals
                    DO g = 1,el%gp
                        kab(ai,bi) = kab(ai,bi)&
                        & + lam*el%Nxg(a,i,g)*el%Nxg(b,j,g) &
                        & + mu *el%Nxg(a,j,g)*el%Nxg(b,i,g)
!                       Kronecker delta term
                        IF (i.eq.j) THEN
                            kab(ai,bi) = kab(ai,bi)&
                            & + mu*(el%Nxg(a,1,g)*el%Nxg(b,1,g)&
                            & + el%Nxg(a,2,g)*el%Nxg(b,2,g))
                        ENDIF
                    ENDDO
                ELSEIF(el%bnd(b,1,j).eq.1) THEN
!               If b is a dirichlet, it affects the F vector
                    DO g = 1,el%gp
                        f(ai) = f(ai) &
                        & +(lam*el%Nxg(a,i,g)*el%Nxg(b,j,g) &
                        & + mu *el%Nxg(a,j,g)*el%Nxg(b,i,g))&
                        & * el%bnd(b,2,j)
                    ENDDO
                ENDIF
            ENDDO
        ENDDO

!       Add into the F matrix
        DO g = 1,el%gp
            f(ai) = f(ai) + el%Ng(a,g)*fun(i,g)
        ENDDO
    ENDIF
ENDDO
ENDDO

END SUBROUTINE elastmat

END MODULE ELASTMOD