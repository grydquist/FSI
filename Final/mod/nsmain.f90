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
INTEGER i,j,k,l,m, ti,dofu, dofp, ts, dof
REAL(KIND=8) lx,ly,tol, dt
!REAL(KIND=8), PARAMETER :: k=1,cv=1,mu=1,e=1,nu=.3,pi = 3.14159265358979323846264337
REAL(KIND=8), ALLOCATABLE :: fun(:,:), G(:), Gt(:,:),Gg(:), Ggt(:,:), dY(:),Ggti(:,:) &
& ,Gti(:,:)

dofu = 2
dofp = 1
dof  = dofp + dofu
dt = 1D-2
! makes the msh
msh = mshtype(dof,dt)

! Initialize time integrator with damping
CALL tintinit(0.2D0)

! get some geometry and whatnot
lx = maxval(msh%x(:,1))
ly = maxval(msh%x(:,2))

! Allocate variables and initializa big stiffness/forcing matrices
ALLOCATE(bnd(msh%np,2,dof),bnd2(msh%np,2,dof))
bnd  = 0
bnd2 = bnd

! Put in my boundary conditions
DO i = 1,msh%np
    DO j =1,dof
        bnd(i,1,j) = 0
        IF (ANY((msh%x(i,:) .gt. maxval(msh%x(:,1))-1d-8)).or.ANY((msh%x(i,:) .lt. 1D-8))) THEN
            IF (j.lt.3) THEN
            bnd(i,1,j) = 1
            bnd(i,2,j) = 0
            ENDIF
        ENDIF
!           Add in cavity-driven BC
        IF((msh%x(i,2) .gt. maxval(msh%x(:,1))-1d-8).and.(j.eq.1)) THEN
            bnd(i,1,j) = 1
            bnd(i,2,j) = 1
        ENDIF
        IF ((abs(msh%x(i,1)-maxval(msh%x(:,1))/2).lt.1D-8).and.(abs(msh%x(i,2)-maxval(msh%x(:,1))/2).lt.1D-8).and.(j.eq.3)) THEN
            bnd(i,1,j) = 1
            bnd(i,2,j) = 0
        ENDIF
    ENDDO
ENDDO
bnd(:,1,3) = 1
bnd(:,2,3) = 0
bnd(5,1,3) = 0

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
& Gt(msh%el(1)%eNoN*dof,msh%el(1)%eNoN*dof), Ggt(msh%np*dof,msh%np*dof),&
& Ggti(msh%np*dof,msh%np*dof), dY(msh%np*dof), Gti(msh%el(1)%eNoN*dof,msh%el(1)%eNoN*dof))

fun(1,:) = 0
tol = 0
u%d    = u%do
p%d    = p%do

open(88,file = 'u.txt',position = 'append')

DO ts = 1,5
!   Reset iteration counter
    ti = 0

!   Make first interation guess
    u%ddot = u%ddoto*(gam - 1D0)/gam

!   Iteration loop
    DO ti = 1,15!WHILE ((ti .lt. 16).or.(tol .lt. 1e-3))
        print *, ti
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
                DO k = 1,msh%eNoN
                    Gg (dof*(msh%IEN(i,k)-1)+j)  = Gg(dof*(msh%IEN(i,k)-1)+j) + G((k-1)*dof+j)
                    DO l=1,dof
                        DO m = 1,msh%eNoN
                            Ggt(dof*(msh%IEN(i,k)-1)+j,dof*(msh%IEN(i,m)-1)+l) = &
                        &   Ggt(dof*(msh%IEN(i,k)-1)+j,dof*(msh%IEN(i,m)-1)+l) + &
                        &   Gt((k-1)*dof+j,(m-1)*dof+l)
                        ENDDO
                    ENDDO
                ENDDO
            ENDDO
        ENDDO

        DO i = 1,msh%np
            DO j = 1,dof
                IF (bnd(i,1,j) .eq. 1) THEN
                    Ggt(dof*(i-1)+j,:) = 0
                    Ggt(dof*(i-1)+j,dof*(i-1)+j) = 1
                    Gg(dof*(i-1)+j) = 0
                ENDIF
            ENDDO
        ENDDO
        DO i = 1,msh%np
            IF(ALL(GGt(:,i).lt.1d-8)) print *, i
        ENDDO
        
        CALL INVERSE(Ggt,Ggti,msh%np*dof)
        DO i = 1,msh%np
            !print *, i ,GG(i*dof-2),GG(i*dof-1),GG(i*dof)
            DO j = 1,dof
            !write(88,*) Gg((i-1)*dof+j)
                !print *, j
                !print *, Ggt((i-1)*dof+j,:)
                DO k = 1,msh%np
                    DO l = 1,dof
                        !write(88,*) Ggt((i-1)*dof+j,(k-1)*dof+l)
                    ENDDO
                ENDDO
            ENDDO
        ENDDO

        dY = matmul(Ggti,-Gg)

        DO i=1,msh%np
            !print *, dY(i*dof-2), dY(i*dof-1), dY(i*dof), i
            DO j=1,dof
                IF (j.lt.3) THEN
                    u%ddot(j,i)  = u%ddot(j,i) + dY((i-1)*dof + j)
                    u%d(j,i)     = u%d(j,i) + u%ddot(j,i)*gam*dt
                    IF(bnd(i,1,j).eq.1)  u%d(j,i) = bnd(i,2,j)
                ELSE
                    p%d(1,i)       = p%d(1,i) + dY((i-1)*dof + j)
                    IF(bnd(i,1,j).eq.1)  p%d(1,i) = bnd(i,2,j)
                ENDIF
            ENDDO
        ENDDO
        stop!if(ts.eq.2) stop
    ENDDO

!   Update Loops
    DO i=1,msh%np
        p%do(1,i) = p%d(1,i)
        DO j = 1,u%dof
            u%do(j,i) = u%d(j,i)
        ENDDO
        write(88,*) u%d(1,i)
        write(88,*) u%d(2,i)
    ENDDO
ENDDO

close(88)

END PROGRAM MAIN
