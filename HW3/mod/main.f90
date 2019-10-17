PROGRAM MAIN
USE MSHMOD
USE LSTIFFMOD
USE GSTIFFMOD
USE SOLVEMOD
IMPLICIT NONE
TYPE(mshtype) msh
INTEGER, ALLOCATABLE :: bnd(:,:)
INTEGER i
REAL(KIND=8) lam,lx,ly,pi,mu, E, nu
REAL(KIND=8), ALLOCATABLE :: KABg(:,:),Fg(:), fun(:,:), kab(:,:),f(:),d(:),Ki(:,:),ftmp(:)

pi = 3.14159265358979323846264337

! makes the msh
msh = mshtype()

! get some geometry and whatnot
lx = maxval(msh%x(:,1))
ly = maxval(msh%x(:,2))
nu = 0.3D0
E = 5e4
lam = nu*E/(1+nu)/(1-2d0*nu)
mu = E/(2d0+2d0*nu)

! Allocate variables and initializa big stiffness/forcing matrices
ALLOCATE(bnd(msh%np,2),kab(msh%eNoN,msh%eNoN),f(msh%eNoN))
ALLOCATE(KABg(msh%np,msh%np),Fg(msh%np),d(msh%np),Ki(msh%np,msh%np))
KABg = 0
Fg = 0
bnd = 0

! Put in my boundary conditions
DO i = 1,msh%np

    IF (ANY((msh%x(i,:) .gt. 7.999999D0)).or.ANY((msh%x(i,:) .lt. 1D-8))) THEN
        bnd(i,1) = 1
        bnd(i,2) = 0
    ENDIF
ENDDO

! Add in the BC info
CALL msh%bound(bnd)

ALLOCATE(fun(msh%el(1)%dof,msh%el(1)%gp),ftmp(msh%el(1)%gp))
fun(1,:) = 0
!fun(2,:) = -1
ftmp = fun(1,:)

! Big loop through elements
DO i=1,msh%nEL

!   Get the local stiffness matrix and f with given forcing function values
    CALL localstiff(kab,f,msh%el(i),ftmp)
!   Put into big matrix
    CALL globalstiff(msh,kab,f,KABg,Fg,i)

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
