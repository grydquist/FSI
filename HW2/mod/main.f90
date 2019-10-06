PROGRAM MAIN
USE MSHMOD
USE LSTIFFMOD
USE GSTIFFMOD
USE SOLVEMOD
IMPLICIT NONE
TYPE(mshtype) msh
INTEGER, ALLOCATABLE :: bnd(:,:)
INTEGER i,g,j
REAL(KIND=8) xt(2),lam,lx,ly,pi
REAL(KIND=8), ALLOCATABLE :: KABg(:,:),Fg(:), fun(:), kab(:,:),f(:),d(:),Ki(:,:)

pi = 3.14159265358979323846264337

! makes the msh
msh = mshtype()

! get some geometry and whatnot
lx = maxval(msh%x(:,1))
ly = maxval(msh%x(:,2))
lam = 4

! Allocate variables and initializa big stiffness/forcing matrices
ALLOCATE(bnd(msh%np,2),kab(msh%eNoN,msh%eNoN),f(msh%eNoN))
ALLOCATE(KABg(msh%np,msh%np),Fg(msh%np),d(msh%np),Ki(msh%np,msh%np))
KABg = 0
Fg = 0
bnd = 0

! Put in my boundary conditions
DO i = 1,msh%np

    IF(msh%x(i,1).lt.1D-8) THEN
        bnd(i,1) = 1
        bnd(i,2) = 1
    ENDIF

    IF ((msh%x(i,2) .gt. 7.999999D0).or.(msh%x(i,2) .lt. 1D-8)) THEN
        bnd(i,1) = 1
        bnd(i,2) = 0
    ENDIF
ENDDO

! Add in the BC info
CALL msh%bound(bnd)

! Big loop through elements
DO i=1,msh%nEL

!   Evaluate forcing function at location of gauss points
    ALLOCATE(fun(msh%el(i)%gp))
    fun =0D0
    DO g = 1,msh%el(i)%gp
        xt = 0D0
        DO j = 1,msh%el(i)%eNoN
            xt = xt + msh%x(msh%IEN(i,j),:)*msh%el(i)%Ng(j,g)
        ENDDO
        fun(g) = cos(lam*PI*(xt(1)/lx+xt(2)/ly))
    ENDDO

!   Get the local stiffness matrix and f with given forcing function values
    CALL localstiff(kab,f,msh%el(i),fun)
!   Put into big matrix
    CALL globalstiff(msh,kab,f,KABg,Fg,i)

    DEALLOCATE(fun)
ENDDO

! Solve!
CALL INVERSE(KABg,Ki,msh%np)
d = matmul(Ki,Fg)

open(88,file = 'd.txt',position = 'append')
DO i = 1,msh%np
write(88,*) d(i)
ENDDO
close(88)

END PROGRAM MAIN
