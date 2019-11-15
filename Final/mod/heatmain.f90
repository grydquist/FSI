PROGRAM MAIN
USE TINT
USE MSHMOD
USE SOLMOD
USE LSTIFFMOD
USE GSTIFFMOD
USE SOLVEMOD
IMPLICIT NONE
TYPE(mshtype) msh
TYPE(soltype) sol
INTEGER, ALLOCATABLE :: bnd(:,:,:), bnd2(:,:,:)
INTEGER i,j, ti,dof, ts
REAL(KIND=8) lx,ly,pi,mu, E, nu, cv, k,tol, dt
REAL(KIND=8), ALLOCATABLE :: KABg(:,:),Fg(:), fun(:,:), kab(:,:),f(:),d(:),Ki(:,:)

! Constants and such
pi = 3.14159265358979323846264337
k = 1
cv = 1
dof = 1
mu=1
e=1
nu=.3
dt = 1e-3

! makes the msh
msh = mshtype(dof)


! Initialize time integrator with damping
CALL tintinit(0.5D0)

! get some geometry and whatnot
lx = maxval(msh%x(:,1))
ly = maxval(msh%x(:,2))

! Allocate variables and initializa big stiffness/forcing matrices
ALLOCATE(bnd(msh%np,2,dof),kab(msh%eNoN*dof,msh%eNoN*dof),f(msh%eNoN*dof))
ALLOCATE(KABg(msh%np*dof,msh%np*dof),Fg(msh%np*dof),d(msh%np*dof),Ki(msh%np*dof,msh%np*dof))
KABg = 0
Fg   = 0
bnd  = 0
bnd2 = bnd


! Put in my boundary conditions
DO i = 1,msh%np
    DO j =1,dof
        IF (ANY((msh%x(i,:) .gt. 7.999999D0)).or.ANY((msh%x(i,:) .lt. 1D-8))) THEN
            bnd(i,1,j) = 1
            bnd(i,2,j) = 0
        ENDIF
        IF((msh%x(i,1) .lt. 1D-8)) THEN
            bnd(i,1,j) = 1
            bnd(i,2,j) = 1
        ENDIF
    ENDDO
ENDDO

! Make solution structure
sol = soltype(dof, msh%np, bnd, bnd2)

! Add in the BC info
CALL msh%bound(bnd)

ALLOCATE(fun(msh%el(1)%dof,msh%el(1)%gp))
fun(1,:) = 0

DO ts = 1,100
!   Reset iteration counter
    ti = 0

!   Make first interation guess
    sol%d    = sol%do
    sol%ddot = sol%ddoto*(gam - 1D0)/gam

!   Iteration loop
    DO WHILE ((tol .lt. 1e-3) .or. (ti<16))
!       Get solutions at alphas
        CALL toalph(sol)

!       Big loop through elements
        DO i=1,msh%nEL
!           Get residual matrix/vector

!           Get the local stiffness matrix and f with given forcing function values
            CALL lstiff(kab,f,fun,msh%el(i),k)
!           Put into big matrix
            CALL gstiffmat(msh,kab,f,KABg,Fg,i)

        ENDDO
        ti = ti+1
    ENDDO
ENDDO

! Solve!
CALL INVERSE(KABg,Ki,msh%np*dof)
d = matmul(Ki,Fg)

open(88,file = 'd.txt',position = 'append')
DO i = 1,msh%np*dof
write(88,*) d(i)
ENDDO
close(88)

END PROGRAM MAIN
