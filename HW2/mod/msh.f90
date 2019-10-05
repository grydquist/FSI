MODULE MSHMOD
    USE ELEMOD
IMPLICIT NONE

TYPE :: mshtype
    INTEGER, ALLOCATABLE :: IEN(:,:)
    REAL(KIND = 8), ALLOCATABLE :: x(:,:)
    INTEGER :: nEL,eNoN, np, nsd
    REAL(KIND = 8), ALLOCATABLE :: xig(:,:),Wg(:)
    TYPE(eltype), ALLOCATABLE :: el(:)
    !REAL(KIND = 8)
    !CONTAINS

    !PROCEDURE :: new=>newmsh

END TYPE

INTERFACE mshType
    PROCEDURE :: newmsh
END INTERFACE mshType

CONTAINS

FUNCTION newmsh() RESULT(msh)
    TYPE(mshtype):: msh
    INTEGER np, nsd, i, eNoN, Nel,elshp
    REAL(KIND =8), ALLOCATABLE :: xt(:,:)

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

    elshp = 1
    ALLOCATE(xt(3,2))
    DO i =1,nEl
        xt = msh%x(msh%IEN(i,:),:)
        msh%el(i) = newel(elshp,xt)
        print *, msh%el(i)%j
    ENDDO


END FUNCTION newmsh


END MODULE MSHMOD