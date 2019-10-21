MODULE GELASTMOD
    USE MSHMOD
IMPLICIT NONE

CONTAINS

!=================================================

! Puts small kab into big KAB given element info
SUBROUTINE gelastmat(msh,kab,f,KABg,Fg,iel)
CLASS(mshtype), INTENT(IN),TARGET :: msh
INTEGER, INTENT(IN) :: iel
REAL(KIND = 8), INTENT(IN) :: kab &
& (msh%el(iel)%eNoN*msh%el(iel)%dof,msh%el(iel)%eNoN*msh%el(iel)%dof), &
& f(msh%el(iel)%eNoN*msh%el(iel)%dof)
REAL(KIND = 8), INTENT(OUT) :: KABg &
& (msh%np*msh%el(iel)%dof,msh%np*msh%el(iel)%dof), &
& Fg(msh%np*msh%el(iel)%dof)
CLASS(eltype), POINTER :: el
INTEGER :: a,i,c, locglo(6)

el => msh%el(iel)

DO a = 1,el%eNoN
    DO i = 1,el%dof
        locglo((a-1)*el%dof+i) = (el%nds(a)-1)*el%dof+i
    ENDDO
ENDDO


c = 0
DO a = 1,el%eNoN
    DO i = 1,el%dof
        c = c+1
        IF (el%bnd(a,1,i).eq.0) THEN
            KABg(locglo(c),locglo(:)) = &
            & KABg(locglo(c),locglo(:)) + kab(c,:)

            Fg(locglo(c)) = Fg(locglo(c)) +f(c)

!       This is just how I deal with Dirichlet BCs
        ELSEIF(el%bnd(a,1,i).eq.1) THEN
            KABg(locglo(c),locglo(c)) = 1
            Fg(locglo(c)) = el%bnd(a,2,i)
        ENDIF
    ENDDO

ENDDO

END SUBROUTINE gelastmat

END MODULE GELASTMOD