MODULE LRESIDNSMOD
    USE SOLMOD
    USE ELEMOD
    USE TINTMOD
IMPLICIT NONE

CONTAINS

SUBROUTINE lresidns(G,Gt,fun,el,u,p)
TYPE(eltype), INTENT(IN) :: el
REAL(KIND=8), ALLOCATABLE, INTENT(IN) :: fun(:,:)
TYPE(soltype), INTENT(IN) :: u,p
REAL(KIND=8), INTENT(OUT) :: G(el%eNoN*el%dof),Gt(el%eNoN*el%dof,el%dof)
INTEGER :: a,b,gp,i,j,ai,bi, ag,bg
! Subgrid variables
REAL(KIND=8) up(u%dof,el%eNoN), pp(p%dof,el%eNoN)
! Temporary material properties !!!
REAL(KIND=8) :: rho,mu, taum(el%eNoN), nuc(el%eNoN)

rho = 1
mu = 1

Gt = 0
G = 0
ai = 0

! Loop through element nodes to calc resid without subgrid
DO a = 1,el%eNoN
ag = el%nds(a)
! Loop through dof (for index a)
DO i = 1,el%dof
    ai = ai+1
!   If this isn't a boundary node, go ahead and calculate the resid at a
    IF (el%bnd(a,1,i).eq.0) THEN
!       Go through surrounding nodes and incorporate into resid for a 
        DO b = 1,el%eNoN
            bg = el%nds(b)
!           Loop through Gauss points to calc integral
!           If i = 1, we're dealing with x-mom, 2->y-mom, 3->cont
            IF((i.eq.1).or.(i.eq.2)) THEN
!               Calculate Rm in direction i for node a (G(1:2))
                DO gp = 1,el%gp
                    G(ai) = G(ai) &
!                   Advection/forcing part
                    &     + el%Ng(a,gp)*rho*(el%Ng(b,gp)*u%ddalphm(i,bg) &
                    &     + el%Ng(b,gp)*el%Nxg(b,1,gp)*u%dalphf(i,bg)*u%dalphf(1,bg) &
                    &     + el%Ng(b,gp)*el%Nxg(b,2,gp)*u%dalphf(i,bg)*u%dalphf(2,bg) &
                    &     - fun(i,gp)) &
!                   Viscous part                   
                    &     + mu/2D0*(u%dalphf(i,bg)*(el%Nxg(a,1,gp)*el%Nxg(b,1,gp)    &
                    &     + el%Nxg(a,2,gp)*el%Nxg(b,2,gp)) &
                    &     + el%Nxg(b,i,gp)*(u%dalphf(1,bg)*el%Nxg(a,1,gp) + u%dalphf(2,bg)*el%Nxg(a,2,gp))) &
!                   Pressure part
                    &     - el%Nxg(a,i,gp)*el%Ng(b,gp)*p%d(1,bg)
!                   Traction part (0 for now)
                ENDDO
            ELSE
!               Calulate Rc for node a (G(3))
                DO gp = 1,el%gp
                    G(ai) = G(ai) &
                    &     + el%Ng(a,gp)*(el%Nxg(b,1,gp)*u%dalphf(1,bg) + el%Nxg(b,2,gp)*u%dalphf(2,bg))
                ENDDO
            ENDIF
        ENDDO
    ENDIF
ENDDO
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Weights

ai = 0

! calculate subgrid variables
DO a = 1,el%eNoN
    ag = el%nds(a)

    taum(a) =(u%dalphf(1,ag)*el%eG(1,1)*u%dalphf(1,ag) &
    &       + u%dalphf(1,ag)*el%eG(1,2)*u%dalphf(2,ag) &
    &       + u%dalphf(2,ag)*el%eG(2,1)*u%dalphf(1,ag) &
    &       + u%dalphf(2,ag)*el%eG(2,2)*u%dalphf(2,ag) &
    &       + 1D0*(mu/rho)**2 &
    &       * (el%eG(1,1)*el%eG(1,1) + el%eG(1,2)*el%eG(1,2) &
    &       +  el%eG(2,1)*el%eG(2,1) + el%eG(2,2)*el%eG(2,2))&
    &       + 4D0/el%dt/el%dt)**(-0.5D0)

    nuc(a) = 1/(el%eG(1,1)+el%eG(2,2))/taum(a)

    DO i =1,u%dof
        ai = ai + 1
        up(i,a) = -taum(a)/rho*G(ai)
    ENDDO
    ai = ai + 1
    pp(1,a) = -rho*nuc(a)*G(ai)
ENDDO

ai = 0

! Loop through element nodes to calc resid with subgrid
DO a = 1,el%eNoN
ag = el%nds(a)
! Loop through dof (for index a)
DO i = 1,el%dof
    ai = ai+1
!   If this isn't a boundary node, go ahead and calculate the resid at a
    IF (el%bnd(a,1,i).eq.0) THEN
!       Go through surrounding nodes and incorporate into resid for a 
        DO b = 1,el%eNoN
            bg = el%nds(b)
!           Loop through Gauss points to calc integral
!           If i = 1, we're dealing with x-mom, 2->y-mom, 3->cont
            IF((i.eq.1).or.(i.eq.2)) THEN
!               Calculate Rm in direction i for node a (G(1:2))
                DO gp = 1,el%gp
                    G(ai) = G(ai) &
!                   First term
                    &     + el%Ng(a,gp)*rho*(up(1,b)*el%Ng(b,gp)*u%dalphf(i,bg)*el%Nxg(b,1,gp) &
                    &     + up(2,b)*el%Ng(b,gp)*u%dalphf(i,bg)*el%Nxg(b,2,gp)) &
!                   Skip second u'*u' term, third term
                    &     - rho*up(i,b)*el%Ng(b,gp)*(el%Nxg(a,1,gp)*el%Ng(b,gp)*u%dalphf(1,bg) &
                    &     + el%Nxg(a,2,gp)*el%Ng(b,gp)*u%dalphf(2,bg)) &
                    &     - el%Nxg(a,i,gp)*pp(1,b)
                ENDDO
            ELSE
!               Calulate Rc for node a (G(3))
                DO gp = 1,el%gp
                    G(ai) = G(ai) &
                    &     - (el%Nxg(a,1,gp)*up(1,b) + el%Nxg(a,2,gp)*up(2,b))
                ENDDO
            ENDIF
        ENDDO
    ENDIF
ENDDO
ENDDO


! Now for the tangent matrix
DO a = 1,el%eNoN
ag = el%nds(a)
! Loop through dof (for index a)
DO i = 1,el%dof
    bi = 0
    ai = ai+1
!   If this isn't a boundary node, go ahead and calculate the resid at a
    IF (el%bnd(a,1,i).eq.0) THEN



!       Calculate tangent matrix
        DO j = 1,el%dof
            DO gp = 1,el%gp
                !Gt(ai,j)= Gt(ai,j) !&
                !&    + alphm*(el%Ng(a,gp))*el%Wg(gp)
            ENDDO
        ENDDO
    ENDIF
ENDDO
ENDDO


END SUBROUTINE lresidns

END MODULE LRESIDNSMOD