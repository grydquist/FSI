MODULE LSTIFFMOD
    USE ELEMOD
IMPLICIT NONE

CONTAINS

! Calculates the local stiffness matrix for a giveen element
SUBROUTINE localstiff(kab,f,el,fun)

CLASS(eltype), INTENT(IN) :: el
REAL(KIND = 8), INTENT(OUT) :: kab(el%eNoN,el%eNoN),f(el%eNoN)
! fun is the function in evaluated at the gauss points
REAL(KIND=8), ALLOCATABLE, INTENT(IN) :: fun(:)
INTEGER :: a,b,g


kab=0
f = 0
DO a = 1,el%eNoN

if (el%bnd(a,1).eq.0) THEN

    DO b = 1,el%eNoN
        IF(el%bnd(b,1).eq.0) THEN

!           Calculate the integrals for interior points
            DO g = 1,el%gp
                kab(a,b) = kab(a,b) + el%Ng(a,g)*el%Ng(b,g)*el%J*el%Wg(g)
            ENDDO

        ELSEIF (el%bnd(b,1).eq.1) THEN
!           Deal with points on same element with boundary nodes
            DO g = 1,el%gp
                f(a) = f(a) - el%Ng(a,g)*el%Ng(b,g)*el%J*el%Wg(g)*el%bnd(b,2)
            ENDDO
        ENDIF

    ENDDO

    DO g = 1,el%gp
        f(a) = f(a)+el%Ng(a,g)*fun(g)*el%J*el%Wg(g)
    ENDDO

ENDIF

ENDDO





END SUBROUTINE localstiff


END MODULE LSTIFFMOD