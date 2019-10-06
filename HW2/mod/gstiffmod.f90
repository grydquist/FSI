MODULE GSTIFFMOD
    USE MSHMOD
IMPLICIT NONE

CONTAINS

! Puts small kab into big KAB given element info
SUBROUTINE globalstiff(msh,kab,f,KABg,Fg,iel)
CLASS(mshtype), INTENT(IN),TARGET :: msh
INTEGER, INTENT(IN) :: iel
REAL(KIND = 8), INTENT(IN) :: kab(msh%el(iel)%eNoN,msh%el(iel)%eNoN),f(msh%el(iel)%eNoN)
REAL(KIND = 8), INTENT(OUT) :: KABg(msh%np,msh%np),Fg(msh%np)
CLASS(eltype), POINTER :: el
INTEGER a

el => msh%el(iel)

DO a = 1,el%eNoN
    IF (el%bnd(a,1).eq.0) THEN
        KABg(el%nds(a),el%nds(:)) = KABg(el%nds(a),el%nds(:)) + kab(a,:)
        Fg(el%nds(a)) = Fg(el%nds(a)) +f(a)

!   This is just how I deal with Dirichlet BCs
    ELSEIF(el%bnd(a,1).eq.1) THEN
        KABg(el%nds(a),el%nds(a)) = 1
        Fg(el%nds(a)) = el%bnd(a,2)
    ENDIF


ENDDO

END SUBROUTINE globalstiff

END MODULE GSTIFFMOD