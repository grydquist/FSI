PROGRAM MAIN
USE TINTMOD
USE MSHMOD
USE SOLMOD
USE LRESIDHEATMOD
USE SOLVEMOD
IMPLICIT NONE
TYPE(mshtype) msh
TYPE(soltype) sol
INTEGER, ALLOCATABLE :: bnd(:,:,:), bnd2(:,:,:)
INTEGER i,j, ti,dof, ts
REAL(KIND=8) lx,ly,tol, dt
!REAL(KIND=8), PARAMETER :: k=1,cv=1,mu=1,e=1,nu=.3,pi = 3.14159265358979323846264337
REAL(KIND=8), ALLOCATABLE :: KABg(:,:),Fg(:), fun(:,:), kab(:,:),&
& f(:),d(:),Ki(:,:), G(:), Gt(:),Gg(:), Ggt(:)

dof = 1
! makes the msh
msh = mshtype(dof)
dt = 1D-2

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
        IF((msh%x(i,2) .gt. 7.999999D0)) THEN
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
ALLOCATE(G(msh%el(1)%eNoN), Gg(msh%np))
Ggt = Gg
Gt = G
fun(1,:) = 0
tol = 1e-2
sol%d    = sol%do

open(88,file = 'd.txt',position = 'append')

DO ts = 1,200
!   Reset iteration counter
    ti = 0

!   Make first interation guess
    sol%ddot = sol%ddoto*(gam - 1D0)/gam

!   Iteration loop
    DO WHILE ((tol .lt. 1e-3) .or. (ti .lt. 16))
        Gg = 0
        Ggt= 0
!       Get solutions at alphas
        CALL toalph(sol)

!       Big loop through elements
        DO i=1,msh%nEL
!           Get residual matrix/vector
            CALL lresidh(G,Gt,fun,msh%el(i),sol)
            Gg(msh%IEN(i,:)) = Gg(msh%IEN(i,:)) + G
            Ggt(msh%IEN(i,:))= Ggt(msh%IEN(i,:)) + Gt
        ENDDO

        DO i=1,msh%np*dof
            IF (abs(Ggt(i)) .gt. 1D-8) THEN
                sol%ddot(1,i)  = sol%ddot(1,i) + Gg(i)/Ggt(i)
                sol%d(1,i) = sol%d(1,i) + sol%ddot(1,i)*gam*dt
            ENDIF
        ENDDO
        
        ti = ti+1

    ENDDO


!   Update Loops
    DO i=1,msh%np*dof
        sol%do(1,i) = sol%d(1,i)
        sol%d(1,i)  = sol%d(1,i) + gam*dt*sol%ddot(1,i)
    ENDDO

    DO i = 1,msh%np*dof
        write(88,*) sol%d(1,i)
    ENDDO
ENDDO

! Solve!
CALL INVERSE(KABg,Ki,msh%np*dof)
d = matmul(Ki,Fg)

close(88)

END PROGRAM MAIN
