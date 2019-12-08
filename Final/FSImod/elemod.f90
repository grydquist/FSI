MODULE ELEMOD
IMPLICIT NONE

TYPE :: eltype
    ! xig(gp,nsd), Wg(gp), Ng(eNoN,gp), Ngx(eNoN,nsd,gp), eG(eNoN,eNoN)
    REAL(KIND = 8),ALLOCATABLE :: dxdxi(:,:), xig(:,:),Wg(:),Ng(:,:),Nxg(:,:,:),eG(:,:)
    REAL(KIND = 8) :: J, dt
    INTEGER:: gp,shp,eNoN,dof,id
    INTEGER, ALLOCATABLE :: bnd(:,:,:), nds(:)
!   Is this element in the fluid or solid domain? 1 == solid, 0 == fluid
    INTEGER :: dom
    CONTAINS
    PROCEDURE :: upgeom => elupgeom
END TYPE

INTERFACE eltype
    PROCEDURE :: newel
END INTERFACE eltype

CONTAINS

!=================================================

FUNCTION newel(elshp,dof,x,id) RESULT(el)
    INTEGER,INTENT(IN) :: elshp, dof
    TYPE(eltype) :: el
    REAL(KIND=8), ALLOCATABLE, INTENT(IN) :: x(:,:)
    INTEGER, INTENT(IN) :: id
    INTEGER i,g,j
    REAL(KIND = 8), ALLOCATABLE ::Nxig(:,:,:),dxidx(:,:)

!   To hopefullly be general for different types of elements
    el%shp = elshp
    el%dof = dof
    el%id = id

!   elshp == 1 is triangular elements
    IF (el%shp .eq. 1) THEN
        el%gp = 3
        el%eNoN = 3
        ALLOCATE(el%nds(3))
        ALLOCATE(el%dxdxi(2,2),el%xig(3,2),el%Wg(3),el%Ng(3,3),dxidx(2,2))
        ALLOCATE(Nxig(el%eNoN,2,el%gp),el%Nxg(el%eNoN,2,el%gp),el%eG(2,2))

        el%xig = 1D0/6D0
        el%xig(2,2) = 2D0/3D0
        el%xig(3,1) = 2D0/3D0
        el%Wg = 1D0/6D0

!       Shape functions evaluated at gp
        el%Ng = 1D0/6D0
        el%Ng(2,2) = 2D0/3D0
        el%Ng(1,3) = 2D0/3D0
        el%Ng(3,1) = 2D0/3D0

!       Calculate derivatives for element
        el%dxdxi(1,1) = x(1,1) - x(3,1)
        el%dxdxi(1,2) = x(2,1) - x(3,1)
        el%dxdxi(2,1) = x(1,2) - x(3,2)
        el%dxdxi(2,2) = x(2,2) - x(3,2)

!       Calculate Jacobian for element
        el%J = el%dxdxi(1,1)*el%dxdxi(2,2)-el%dxdxi(2,1)*el%dxdxi(1,2)

!       Get inverse of these derivatives
        dxidx(1,1) = el%dxdxi(2,2)
        dxidx(2,2) = el%dxdxi(1,1)
        dxidx(1,2) = -el%dxdxi(1,2)
        dxidx(2,1) = -el%dxdxi(2,1)
        dxidx = dxidx/el%J

!       Evaluate at Gauss points for element
        DO i =1,3
            el%Ng(1,i) = el%xig(i,1)
            el%Ng(2,i) = el%xig(i,2)
            el%Ng(3,i) = 1 - el%xig(i,1) - el%xig(i,2)
        ENDDO

!       Get shape function derivatives at Gauss Points in parent domain
        Nxig = 0
        Nxig(1,1,:) = 1
        Nxig(2,2,:) = 1
        Nxig(3,:,:) = -1

!       Now convert that to physical domain
        DO g =1,el%gp
            DO i =1,el%eNoN
                DO j = 1,2
                    el%Nxg(i,j,g) = Nxig(i,1,g)*dxidx(1,j) + Nxig(i,2,g)*dxidx(2,j)
                ENDDO
            ENDDO
        ENDDO

!       Calculate element metric tensor
        el%eG = 0

        DO i = 1,2
            DO j = 1,2
                el%eG(i,j) = dxidx(1,j)*dxidx(1,i) + dxidx(2,j)*dxidx(2,i)
            ENDDO
        ENDDO

        ALLOCATE(el%bnd(el%eNoN,2,el%dof))

        el%bnd = 0

    END IF
END FUNCTION newel

!=====================================
SUBROUTINE elupgeom(el,x)
    CLASS(eltype), INTENT(INOUT) :: el
    REAL(KIND=8), INTENT(IN) :: x(el%eNoN,el%dof)
    INTEGER i,g,j
    REAL(KIND = 8), ALLOCATABLE ::Nxig(:,:,:),dxidx(:,:)


!   elshp == 1 is triangular elements
    IF (el%shp .eq. 1) THEN
        ALLOCATE(Nxig(el%eNoN,2,el%gp),dxidx(2,2))

!       Calculate derivatives for element
        el%dxdxi(1,1) = x(1,1) - x(3,1)
        el%dxdxi(1,2) = x(2,1) - x(3,1)
        el%dxdxi(2,1) = x(1,2) - x(3,2)
        el%dxdxi(2,2) = x(2,2) - x(3,2)

!       Calculate Jacobian for element
        el%J = el%dxdxi(1,1)*el%dxdxi(2,2)-el%dxdxi(2,1)*el%dxdxi(1,2)

!       Get inverse of these derivatives
        dxidx(1,1) = el%dxdxi(2,2)
        dxidx(2,2) = el%dxdxi(1,1)
        dxidx(1,2) = -el%dxdxi(1,2)
        dxidx(2,1) = -el%dxdxi(2,1)
        dxidx = dxidx/el%J

!       Evaluate at Gauss points for element
        DO i =1,3
            el%Ng(1,i) = el%xig(i,1)
            el%Ng(2,i) = el%xig(i,2)
            el%Ng(3,i) = 1 - el%xig(i,1) - el%xig(i,2)
        ENDDO

!       Get shape function derivatives at Gauss Points in parent domain
        Nxig = 0
        Nxig(1,1,:) = 1
        Nxig(2,2,:) = 1
        Nxig(3,:,:) = -1

!       Now convert that to physical domain
        DO g =1,el%gp
            DO i =1,el%eNoN
                DO j = 1,2
                    el%Nxg(i,j,g) = Nxig(i,1,g)*dxidx(1,j) + Nxig(i,2,g)*dxidx(2,j)
                ENDDO
            ENDDO
        ENDDO

!       Calculate element metric tensor
        el%eG = 0

        DO i = 1,2
            DO j = 1,2
                el%eG(i,j) = dxidx(1,j)*dxidx(1,i) + dxidx(2,j)*dxidx(2,i)
            ENDDO
        ENDDO

    END IF


END SUBROUTINE elupgeom

END MODULE ELEMOD
