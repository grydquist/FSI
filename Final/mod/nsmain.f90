PROGRAM MAIN
USE TINTMOD
USE MSHMOD
USE SOLMOD
USE LRESIDNSMOD
USE SOLVEMOD
IMPLICIT NONE
TYPE(mshtype) msh
TYPE(soltype) u,p
INTEGER, ALLOCATABLE :: bnd(:,:,:), bnd2(:,:,:), bndt(:,:,:),bnd2t(:,:,:)
INTEGER i,j, ti,dofu, dofp, ts, dof
REAL(KIND=8) lx,ly,tol, dt
!REAL(KIND=8), PARAMETER :: k=1,cv=1,mu=1,e=1,nu=.3,pi = 3.14159265358979323846264337
REAL(KIND=8), ALLOCATABLE :: fun(:,:), G(:), Gt(:,:),Gg(:), Ggt(:,:)

dofu = 2
dofp = 1
dof  = dofp + dofu
! makes the msh
msh = mshtype(dof)
dt = 1D-2

! Initialize time integrator with damping
CALL tintinit(0.5D0)

! get some geometry and whatnot
lx = maxval(msh%x(:,1))
ly = maxval(msh%x(:,2))

! Allocate variables and initializa big stiffness/forcing matrices
ALLOCATE(bnd(msh%np,2,dof))
bnd  = 0
bnd2 = bnd

! Put in my boundary conditions
DO i = 1,msh%np
    DO j =1,dof
        IF (ANY((msh%x(i,:) .gt. 7.999999D0)).or.ANY((msh%x(i,:) .lt. 1D-8)).and.(j.lt.3)) THEN
            bnd(i,1,j) = 1
            bnd(i,2,j) = 0
        ENDIF
!           Add in cavity-driven BC
        IF((msh%x(i,2) .gt. 7.999999D0).and.(j.eq.1)) THEN
            bnd(i,1,j) = 1
            bnd(i,2,j) = 1
        ENDIF
        IF ((msh%x(i,1).lt.1D-8).and.(msh%x(i,2).lt.1D-8).and.(j.eq.3)) THEN
            bnd(i,1,j) = 1
            bnd(i,2,j) = 0
        ENDIF
    ENDDO
ENDDO

! Make solution structure
ALLOCATE(bndt(msh%np,2,dofu),bnd2t(msh%np,2,dofu))
bndt = bnd(:,:,1:2)
bnd2t= bnd2(:,:,1:2)
u = soltype(dofu, msh%np, bndt, bnd2t)
DEALLOCATE(bndt,bnd2t)
ALLOCATE(bndt(msh%np,2,dofp),bnd2t(msh%np,2,dofp))
bndt(:,:,1) = bnd(:,:,3)
bnd2t(:,:,1)= bnd2(:,:,3)
p = soltype(dofp, msh%np, bndt, bnd2t)


! Add in the BC info
CALL msh%bound(bnd)

ALLOCATE(fun(dof,msh%el(1)%gp))
ALLOCATE(G(msh%el(1)%eNoN*dof), Gg(msh%np*dof), &
& Gt(msh%el(1)%eNoN*dof,dof), Ggt(msh%np*dof,dof))

fun(1,:) = 0
tol = 1e-2
u%d    = u%do
p%d    = p%do

open(88,file = 'd.txt',position = 'append')

DO ts = 1,200
!   Reset iteration counter
    ti = 0

!   Make first interation guess
    u%ddot = u%ddoto*(gam - 1D0)/gam
    p%ddot = p%ddoto*(gam - 1D0)/gam

!   Iteration loop
    DO WHILE ((tol .lt. 1e-3) .or. (ti .lt. 16))
        Gg = 0
        Ggt= 0
!       Get solutions at alphas
        CALL toalph(u)

!       Big loop through elements
        DO i=1,msh%nEL
!           Get residual matrix/vector
            CALL lresidns(G,Gt,fun,msh%el(i),u,p)
!           Put back into global
            DO j=1,dof
                Gg (dof*(msh%IEN(i,:)-1)+j)  = Gg(dof*(msh%IEN(i,:)-1)+j) + G
                Ggt(dof*(msh%IEN(i,:)-1)+j,:)= Ggt(dof*(msh%IEN(i,:)-1)+j,:) + Gt
            ENDDO
        ENDDO

        DO i=1,msh%np*dof
            DO j=1,dof
                IF (j.lt.3) THEN
                   ! u%ddot(j,i)  = u%ddot(j,i) + Gg(i)/Ggt(i)
                   ! u%d(j,i) = u%d(1,i) + p%ddot(1,i)*gam*dt
                ELSE

                ENDIF
            ENDDO
        ENDDO
        
        ti = ti+1

    ENDDO


!   Update Loops
    DO i=1,msh%np
        DO j = 1,dof
            !p%do(1,i) = p%d(1,i)
            !p%d(1,i)  = p%d(1,i) + gam*dt*p%ddot(1,i)
        ENDDO
    ENDDO

    DO i = 1,msh%np*dof
        !write(88,*) p%d(1,i)
    ENDDO
ENDDO

close(88)

END PROGRAM MAIN
