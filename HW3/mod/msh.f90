MODULE MSHMOD
    USE ELEMOD
IMPLICIT NONE

TYPE :: mshtype
!   Connectivity (Nel, eNoN)
    INTEGER, ALLOCATABLE :: IEN(:,:)
!   Node locations (np,nsd)
    REAL(KIND = 8), ALLOCATABLE :: x(:,:)
!   Number of elements, number of nodes per element, number of nodes, number of dimensions
    INTEGER :: nEL,eNoN, np, nsd
!   Element type with info in it
    TYPE(eltype), ALLOCATABLE :: el(:)
!   List of if the node is a boundary node, and if so what type and value
!   (np,2), 1 == Dirichlet, 2 == Neumann, 0 == no boundary
    INTEGER, ALLOCATABLE :: bnd(:,:)

    CONTAINS
    PROCEDURE :: bound=>mshbound
END TYPE

INTERFACE mshType
    PROCEDURE :: newmsh
END INTERFACE mshType

CONTAINS

!=================================================

! Makes a new mesh
FUNCTION newmsh() RESULT(msh)
    TYPE(mshtype):: msh
    INTEGER np, nsd, i, eNoN, Nel,elshp, dof
    REAL(KIND =8), ALLOCATABLE :: xt(:,:)

!   Reads x and IEN
    open(10,file = 'x.txt')
    read(10,*) np,nsd
    ALLOCATE(msh%x(np,nsd))
    DO i = 1,np
            read(10,*)  msh%x(i,1),msh%x(i,2)
    ENDDO
    msh%np = np
    msh%nsd = nsd
    close(10)

    open(11,file = 'IEN.txt')
    read(11,*) Nel, EnoN
    ALLOCATE(msh%IEN(Nel,eNoN))
    DO i = 1,Nel
            read(11,*) msh%IEN(i,1),msh%IEN(i,2),msh%IEN(i,3)
    ENDDO
    msh%Nel = Nel
    msh%eNoN = eNoN
    close(11)

    ALLOCATE(msh%el(Nel))

!   Allocate every element in msh with additional info
    elshp = 1
    dof = 1!nsd
    ALLOCATE(xt(3,2))
    DO i =1,nEl
        xt = msh%x(msh%IEN(i,:),:)
        msh%el(i) = newel(elshp,dof,xt)
        msh%el(i)%nds = msh%IEN(i,:)
    ENDDO

END FUNCTION newmsh

! Deals with BCs
SUBROUTINE mshbound(msh,bnd)
    CLASS(mshtype),INTENT(INOUT) :: msh
    INTEGER, ALLOCATABLE, INTENT(IN) :: bnd(:,:)
    INTEGER i,j

    msh%bnd = bnd

    DO i=1,msh%nEL
        DO j = 1,msh%eNoN
            IF (bnd(msh%IEN(i,j),1) .ne. 0) THEN
                msh%el(i)%bnd(j,1) = bnd(msh%IEN(i,j),1)
                msh%el(i)%bnd(j,2) = bnd(msh%IEN(i,j),2)
            ENDIF
        ENDDO
    ENDDO

END SUBROUTINE mshbound

END MODULE MSHMOD