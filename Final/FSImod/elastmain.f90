PROGRAM MAIN
USE TINTMOD
USE MSHMOD
USE SOLMOD
USE LRESIDELASTMOD
USE SOLVEMOD
IMPLICIT NONE
TYPE(mshtype) msh
TYPE(soltype) yd, y
INTEGER, ALLOCATABLE :: bnd(:,:,:), bnd2(:,:,:), bndt(:,:,:),bnd2t(:,:,:)
INTEGER i,j,k,l,m, ti, dofy, ts, dof, gp,a, cnt
REAL(KIND=8) lx,ly,tol, dt
!REAL(KIND=8), PARAMETER :: k=1,cv=1,mu=1,e=1,nu=.3,pi = 3.14159265358979323846264337
REAL(KIND=8), ALLOCATABLE :: fun(:,:),fung(:,:), G(:), Gt(:,:),Gg(:), Ggt(:,:), dY(:),Ggti(:,:) &
& ,Gti(:,:), yorig(:,:), xog(:,:)

cnt = 0
dofy = 2
dof  = dofy
dt = 5D-3
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
!       Left edge fixed
        IF ((msh%x(i,1) .lt. minval(msh%x(:,1))+1d-8)) THEN
            IF (j.lt.3) THEN
            bnd(i,1,j) = 1
            bnd(i,2,j) = 0
            ENDIF
        ENDIF
    ENDDO
ENDDO

! Make solution structure
ALLOCATE(bndt(msh%np,2,dofy),bnd2t(msh%np,2,dofy))
bndt = bnd(:,:,1:2)
bnd2t= bnd2(:,:,1:2)
y = soltype(dofy, msh%np, bndt, bnd2t)
yd= soltype(dofy, msh%np, bndt, bnd2t)

! Add in the BC info
CALL msh%bound(bnd)

ALLOCATE(fun(dof,msh%el(1)%gp))
ALLOCATE(G(msh%el(1)%eNoN*dof), Gg(msh%np*dof), &
& Gt(msh%el(1)%eNoN*dof,msh%el(1)%eNoN*dof), Ggt(msh%np*dof,msh%np*dof),&
& Ggti(msh%np*dof,msh%np*dof), dY(msh%np*dof), Gti(msh%el(1)%eNoN*dof,msh%el(1)%eNoN*dof), &
& yorig(msh%np,y%dof), fung(dof,msh%np), xog(msh%np,y%dof))

! Original mesh coordinates
xog = msh%x

! Put in forcing at nodes
DO i = 1, msh%np
    DO j = 1,dof
!       pulling down on bottom right
        IF ((msh%x(i,1) .gt. maxval(msh%x(:,1))-1d-8).and.((msh%x(i,2) .lt. 1D-8))) THEN
            fung(2,i) = -50D0
        ENDIF
    ENDDO
ENDDO

tol = 0
yd%d    = yd%do
y%d     = y%do
yorig  = y%d

open(88,file = 'y.txt',position = 'append')

DO ts = 1,50
!   Reset iteration counter
    ti = 0

!   Stop pulling down
    IF(ts.eq.25) fung = 0

!   Make first interation guess
    yd%ddot = yd%ddoto*(gam - 1D0)/gam

!   Iteration loop
    DO ti = 1,15!WHILE ((ti .lt. 16).and.(tol .lt. 1e-3))
        Gg = 0
        Ggt= 0
        dY = 0
!       Get solutions at alphas
        CALL toalph(y)

!       Big loop through elements
        DO i=1,msh%nEL

            fun = 0
!           Calculate forcing function at gauss points of current element
            DO j = 1,y%dof
                DO a = 1,msh%el(i)%eNoN
                    DO gp = 1,msh%el(i)%gp
                        IF (j.eq.2) THEN
                            fun(j,gp) = fun(j,gp) + fung(j,msh%el(i)%nds(a))*msh%el(i)%Ng(a,gp)
                        ENDIF
                    ENDDO
                ENDDO
            ENDDO

!           Get residual matrix/vector
            CALL lresidelast(G,Gt,fun,msh%el(i),y)
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
                IF (bnd(i,1,j) .eq. 1) THEN
                    Ggt(dof*(i-1)+j,:) = 0
                    Ggt(dof*(i-1)+j,dof*(i-1)+j) = 1
                    Gg(dof*(i-1)+j) = 0
                ENDIF
            ENDDO
        ENDDO
        
        DO i = 1,msh%np
            !print *, i, GG(i*dof-1),GG(i*dof)
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
        CALL INVERSE(Ggt,Ggti,msh%np*dof)
        dY = matmul(Ggti,-Gg)

        DO i=1,msh%np
            !print *, dY(i*dof-1), dY(i*dof), i
            DO j=1,dof
                y%ddd(j,i) = y%ddd(j,i) +           dY((i-1)*dof + j)
                y%ddot(j,i)= y%ddot(j,i)+ gam*dt*   dY((i-1)*dof + j)
                y%d(j,i)   = y%d(j,i)   + bet*dt*dt*dY((i-1)*dof + j)
                msh%x(i,j) = xog(i,j) + y%d(j,i)
            ENDDO
        ENDDO

!       Get updated element geometric properties
        DO i =1,msh%nEl
            CALL msh%el(i)%upgeom(msh%x(msh%IEN(i,:),:))
        ENDDO


        print *, ti, ts, maxval(abs(Gg)), minval((dy))
        !if(ti.eq.2) stop
    ENDDO

!   Update Loops
    DO i=1,msh%np
        DO j = 1,y%dof
            y%dddo(j,i) = y%ddd(j,i)
            y%ddoto(j,i) = y%ddot(j,i)
            y%do(j,i) = y%d(j,i)
        ENDDO
        write(88,*) msh%x(i,1)
        write(88,*) msh%x(i,2)
    ENDDO
ENDDO
close(88)

END PROGRAM MAIN
