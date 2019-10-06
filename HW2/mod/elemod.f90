MODULE ELEMOD
IMPLICIT NONE

TYPE :: eltype
    REAL(KIND = 8),ALLOCATABLE :: dxdxi(:,:), xig(:,:),Wg(:),Ng(:,:)
    REAL(KIND = 8) :: J
    INTEGER:: gp,shp,eNoN
    INTEGER, ALLOCATABLE :: bnd(:,:), nds(:)
END TYPE

INTERFACE eltype
    PROCEDURE :: newel
END INTERFACE eltype

CONTAINS

!=================================================

FUNCTION newel(elshp,x) RESULT(el)
    INTEGER,INTENT(IN) :: elshp
    TYPE(eltype) :: el
    REAL(KIND=8), ALLOCATABLE, INTENT(IN) :: x(:,:)
    INTEGER i

!   To hopefullly be general for different types of elements
    el%shp = elshp

!   elshp == 1 is triangular elements
    IF (el%shp .eq. 1) THEN
        el%gp = 3
        el%eNoN = 3
        ALLOCATE(el%nds(3))
        ALLOCATE(el%dxdxi(2,2),el%xig(3,2),el%Wg(3),el%Ng(3,3))

        el%xig = 1D0/6D0
        el%xig(2,2) = 2D0/3D0
        el%xig(3,1) = 2D0/3D0
        el%Wg = 1D0/6D0

        el%Ng = 1D0/6D0
        el%Ng(2,2) = 2D0/3D0
        el%Ng(1,3) = 2D0/3D0
        el%Ng(3,1) = 2D0/3D0

!      Calculate derivatives for element
        el%dxdxi(1,1) = x(1,1) - x(3,1)
        el%dxdxi(1,2) = x(2,1) - x(3,1)
        el%dxdxi(2,1) = x(1,2) - x(3,2)
        el%dxdxi(2,2) = x(2,2) - x(3,2)

!       Calculate Jacobian for element
        el%J = el%dxdxi(1,1)*el%dxdxi(2,2)-el%dxdxi(2,1)*el%dxdxi(1,2)

!       Evaluate at Gauss points for element
        DO i =1,3
            el%Ng(1,i) = el%xig(i,1)
            el%Ng(2,i) = el%xig(i,2)
            el%Ng(3,i) = 1 - el%xig(i,1) - el%xig(i,2)
        ENDDO

        ALLOCATE(el%bnd(el%eNoN,2))

        el%bnd = 0

    END IF
END FUNCTION newel

END MODULE ELEMOD
