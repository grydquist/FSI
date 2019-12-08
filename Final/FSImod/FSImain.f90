PROGRAM MAIN
USE TINTMOD
USE MSHMOD
USE SOLMOD
USE LRESIDNSMOD
USE LRESIDELASTMOD
USE SOLVEMOD
IMPLICIT NONE
TYPE(mshtype) msh
TYPE(soltype) u,p,y
INTEGER, ALLOCATABLE :: bnd(:,:,:), bnd2(:,:,:), bndt(:,:,:),bnd2t(:,:,:)
INTEGER i,j,k,l,m, ti,dofu, dofp,dofy, ts, dof,cnt
REAL(KIND=8) tol, dt, maxG, rtol, t
!REAL(KIND=8), PARAMETER :: k=1,cv=1,mu=1,e=1,nu=.3,pi = 3.14159265358979323846264337
REAL(KIND=8), ALLOCATABLE :: fun(:,:), G(:), Gt(:,:),Gg(:), Ggt(:,:), dY(:),Ggti(:,:) &
& ,Gti(:,:), xog(:,:)

cnt = 0
dofu = 2
dofy = 2
dofp = 1
dof  = dofp + dofu + dofy
dt = 1D-3
! makes the msh
msh = mshtype(dof,dt)

! Initialize time integrator with damping
CALL tintinit(0.2D0)

! Allocate variables and initializa big stiffness/forcing matrices
ALLOCATE(bnd(msh%np,2,dof),bnd2(msh%np,2,dof))
bnd  = 0
bnd2 = bnd

! Put in my boundary conditions
DO i = 1,msh%np
    DO j =1,dof
        bnd(i,1,j) = 0
!       Add in BCs for proj
!       Inlet
        IF (msh%x(i,1) .lt. 1D-8) THEN
            bnd(i,1,1) = 1
            bnd(i,2,1) = 1
            bnd(i,1,2) = 1
            bnd(i,2,2) = 0
        ENDIF
            
!       Top Walls
        IF((msh%x(i,2) .gt. maxval(msh%x(:,2))-1d-8).or.(msh%x(i,2) .lt. 1d-8)) THEN
            bnd(i,1,2) = 1
            bnd(i,2,2) = 0
        ENDIF
!       Square
        IF((msh%x(i,1) .gt. 4.5 - 1d-8) .and. (msh%x(i,1) .lt. 5.5 + 1d-8) .and. &
        &  (msh%x(i,2) .gt. 5.5 - 1d-8) .and. (msh%x(i,2) .lt. 6.5 + 1d-8)) THEN
            bnd(i,1,1) = 1
            bnd(i,2,1) = 0
            bnd(i,1,2) = 1
            bnd(i,2,2) = 0
        ENDIF
!       Fluid solid boundary
        IF((msh%x(i,2) .gt. 5.97 - 1d-8) .and. (msh%x(i,2) .lt. 5.97 + 1d-8)) THEN
            bnd(i,1,4:5) = 4
        ENDIF
        IF((msh%x(i,2) .gt. 6.03 - 1d-8) .and. (msh%x(i,2) .lt. 6.03 + 1d-8)) THEN
            bnd(i,1,4:5) = 4
        ENDIF
!       Solid Bar (fixed at left edge)
        IF((msh%x(i,1) .gt. 5.5 - 1d-8)  .and. (msh%x(i,1) .lt. 5.5  + 1d-8) .and. &
        &  (msh%x(i,2) .gt. 5.97 - 1d-8) .and. (msh%x(i,2) .lt. 6.03 + 1d-8)) THEN
            bnd(i,1,4:5) = 1
            bnd(i,2,4:5) = 0
        ENDIF
    ENDDO
ENDDO
! Reference pressure
bnd(5,1,3) = 1
bnd(5,2,3) = 0

! Test BCs !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
bnd  = 0
bnd2 = bnd
DO i = 1,msh%np
    DO j =1,dof
        bnd(i,1,j) = 0
        IF (ANY((msh%x(i,:) .gt. maxval(msh%x(:,1))-1d-8)).or.ANY((msh%x(i,:) .lt. 1D-8))) THEN
            IF (j.lt.3) THEN
            bnd(i,1,j) = 1
            bnd(i,2,j) = 0
            ENDIF
        ENDIF
!       Add in cavity-driven BC
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
! Fix one point on right side
bnd(3,1,4:5) = 1
bnd(3,2,4:5) = 0
bnd(6,1,4:5) = 1
bnd(6,2,4:5) = 0
bnd(5,1,4:5) = 4

!!!!!!!!!!!!!!!!!!!!!!!!!!!


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
DEALLOCATE(bndt,bnd2t)

ALLOCATE(bndt(msh%np,2,dofy),bnd2t(msh%np,2,dofy))
bndt(:,:,1:2) = bnd(:,:,4:5)
bnd2t(:,:,1:2)= bnd2(:,:,4:5)
y = soltype(dofy, msh%np, bndt, bnd2t)

! Add in the BC info
CALL msh%bound(bnd)

ALLOCATE(fun(dof,msh%el(1)%gp))
ALLOCATE(G(msh%el(1)%eNoN*dof), Gg(msh%np*dof), &
& Gt(msh%el(1)%eNoN*dof,msh%el(1)%eNoN*dof), Ggt(msh%np*dof,msh%np*dof),&
& Ggti(msh%np*dof,msh%np*dof), dY(msh%np*dof), Gti(msh%el(1)%eNoN*dof,msh%el(1)%eNoN*dof), &
& xog(msh%np,y%dof))

! Original mesh coordinates
xog = msh%x

fun(1,:) = 0
u%d    = u%do
p%d    = p%do
y%d    = y%do

! What elements are in what domain?
DO i = 1,msh%nEl
    IF (ANY(msh%x(msh%IEN(i,:),1) .gt. 5.5  + 1D-8) .and. &
    &   ANY(msh%x(msh%IEN(i,:),1) .lt. 9.5  - 1D-8) .and. &
    &   ANY(msh%x(msh%IEN(i,:),2) .gt. 5.97 + 1D-8) .and. &
    &   ANY(msh%x(msh%IEN(i,:),2) .lt. 6.03 - 1D-8)) THEN
!       In solid
        msh%el(i)%dom = 1
    ELSE
!       In fluid
        msh%el(i)%dom = 0
    ENDIF
ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!! Test
msh%el(8)%dom = 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!


t = -dt
DO ts = 1,1000
!   Reset iteration counter
    t = t + dt
    ti = 0
    tol = 100

!   Make first interation guess
    u%ddot = u%ddoto*(gam - 1D0)/gam
    y%ddd  = y%dddo *(gam - 1D0)/gam

!   Iteration loop
    DO WHILE ((ti .lt. 16).and.(tol .gt. 1e-5))
        ti = ti+1
        Gg = 0
        Ggt= 0
        dY = 0
!       Get solutions at alphas
        CALL toalph(u)
        CALL toalph(y)

!       Big loop through elements
        DO i=1,msh%nEL
!           Get residual matrix/vector, appropriate for what domain we're in
            IF (msh%el(i)%dom .eq. 1) THEN ! solid
                !CALL lresidelast(G,Gt,fun,msh%el(i),y)
            ELSE ! liquid
                CALL lresidns(G,Gt,fun,msh%el(i),u,p)
            ENDIF
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

!       Boundary conditions: just make rows equal to 0 and 1 on diag
        DO i = 1,msh%np
            DO j = 1,dof
                IF (bnd(i,1,j) .ne. 0) THEN
                    Ggt(dof*(i-1)+j,:) = 0
                    Ggt(dof*(i-1)+j,dof*(i-1)+j) = 1
                    Gg(dof*(i-1)+j) = 0
                ENDIF
            ENDDO
        ENDDO
        
        DO i = 1,msh%np
            !print *, i ,GG(i*dof-4),GG(i*dof-3),GG(i*dof-2),GG(i*dof-1),GG(i*dof)
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

!       Get largest value of residual for relative tolerance
        IF(ts.eq.1.and.ti.eq.1) maxG = maxval(abs(Gg))

        CALL qr_solve(msh%np*dof,msh%np*dof, Ggt,-Gg,dY)

        DO i=1,msh%np
            !print *, dY(i*dof-4), dY(i*dof-3), dY(i*dof-2), i
            DO j=1,dof
                IF (j.lt.3) THEN
                    y%ddd(j,i) = y%ddd(j,i) +           dY((i-1)*dof + j)
                    y%ddot(j,i)= y%ddot(j,i)+ gam*dt*   dY((i-1)*dof + j)
                    y%d(j,i)   = y%d(j,i)   + bet*dt*dt*dY((i-1)*dof + j)
                    msh%x(i,j) = xog(i,j) + y%d(j,i)
                    u%ddot(j,i)  = u%ddot(j,i) + dY((i-1)*dof + j)
                    u%d(j,i)     = u%d(j,i) + gam*dt*dY((i-1)*dof + j)
                ELSE
                    p%d(1,i)     = p%d(1,i) + dY((i-1)*dof + j)
                ENDIF
            ENDDO
        ENDDO

!       Get updated element geometric properties
        DO i =1,msh%nEl
            CALL msh%el(i)%upgeom(msh%x(msh%IEN(i,:),:))
        ENDDO

        tol = dnrm22(msh%np*dof,GG,1)
        rtol = dnrm22(msh%np*dof,GG/maxG,1)
        print *, ti, ts, minval((dy)), tol, rtol
        !if(ti.eq.10) stop
    ENDDO

!   Update Loops

    open(99,file = 'u.txt',position = 'append')
    DO i=1,msh%np
        p%do(1,i) = p%d(1,i)
        DO j = 1,u%dof
            u%ddoto(j,i) = u%ddot(j,i)
            u%do(j,i) = u%d(j,i)
        ENDDO
        write(99,*) u%d(1,i)
        write(99,*) u%d(2,i)
    ENDDO
    close(99)
ENDDO


END PROGRAM MAIN
