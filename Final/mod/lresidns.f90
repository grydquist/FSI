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
REAL(KIND=8) up(u%dof,el%eNoN), pp(p%dof,el%eNoN), ugp(u%dof,el%gp), &
&            duxgp(u%dof,u%dof,el%gp), dutgp(u%dof,el%gp), pgp(el%gp), &
&            upgp(u%dof,el%gp), ppgp(el%gp)
! Temporary material properties !!!
REAL(KIND=8) :: rho,mu, taum(el%eNoN), nuc(el%eNoN)

rho = 1
mu = 1

Gt = 0
G = 0
ai = 0
ugp = 0
duxgp = 0
dutgp = 0
pgp = 0

! First let's find u and it's derivatives at the Gauss points
DO i = 1,u%dof
    DO gp = 1,el%gp
        DO a = 1,el%eNoN
            ag = el%nds(a)
            ugp(i,gp)   = ugp(i,gp)  + el%Ng(a,gp)*u%dalphf(i,ag)
            dutgp(i,gp) = dutgp(i,gp)+ el%Ng(a,gp)*u%ddalphm(i,ag)
!           This term is deriv of ith velocity in jth direction at gp
            DO j = 1,u%dof
                duxgp(i,j,gp) = duxgp(i,j,gp) + el%Nxg(i,j,gp)*u%dalphf(i,ag)
            ENDDO
        ENDDO
    ENDDO
ENDDO

! Now p as well
DO gp = 1,el%gp
    DO a = 1,el%eNoN
        ag = el%nds(a)
        pgp(gp) = pgp(gp) + el%Ng(a,gp)*p%d(1,ag)
    ENDDO
ENDDO

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
                    &     +(el%Ng(a,gp)*rho*(dutgp(i,gp) &
                    &     + ugp(1,gp)*duxgp(i,1,gp) &
                    &     + ugp(2,gp)*duxgp(i,2,gp) &
                    &     - fun(i,gp)) &
!                   Viscous part                   
                    &     + mu/2D0 &
                    &     *(duxgp(i,1,gp)*el%Nxg(a,1,gp) &
                    &     + duxgp(i,2,gp)*el%Nxg(a,2,gp) &
                    &     + duxgp(1,i,gp)*el%Nxg(a,1,gp) &
                    &     + duxgp(2,1,gp)*el%Nxg(a,2,gp))&
!                   Pressure part
                    &     - el%Nxg(a,i,gp)*pgp(gp))*el%Wg(gp)
!                   Traction part (0 for now)
                ENDDO
            ELSE
!               Calulate Rc for node a (G(3))
                DO gp = 1,el%gp
                   G(ai) = G(ai) &
                   &     + (el%Ng(a,gp) &
                   &     * (duxgp(1,1,gp) + duxgp(2,2,gp)))*el%Wg(gp)
                ENDDO
            ENDIF
        ENDDO
    ENDIF
ENDDO
ENDDO

ai = 0
upgp = 0
ppgp = 0

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


! Calc subgrid variables at gauss points
DO i = 1,u%dof
    DO gp = 1,el%gp
        DO a = 1,el%eNoN
            upgp(i,gp) = upgp(i,gp)  + el%Ng(a,gp)*up(i,a)
        ENDDO
    ENDDO
ENDDO

! Now pp at gp as well
DO gp = 1,el%gp
    DO a = 1,el%eNoN
        ppgp(gp) = ppgp(gp) + el%Ng(a,gp)*pp(1,a)
    ENDDO
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
!           Loop through Gauss points to calc integral
!           If i = 1, we're dealing with x-mom, 2->y-mom, 3->cont
            IF((i.eq.1).or.(i.eq.2)) THEN
!               Calculate Rm in direction i for node a (G(1:2))
                DO gp = 1,el%gp
                    G(ai) = G(ai) &
!                   First term
                    &     +(el%Ng(a,gp)*rho*(upgp(1,gp)*duxgp(i,1,gp) &
                    &     + upgp(2,gp)*duxgp(i,2,gp)) &
!                   Skip second u'*u' term, third term
                    &     +(ugp(1,gp)*el%Nxg(a,1,gp)  &
                    &     + ugp(2,gp)*el%Nxg(a,2,gp)) &
                    &     * upgp(i,gp)*rho &
!                   Pressure term
                    &     - el%Nxg(a,i,gp)*ppgp(gp))*el%Wg(gp)
                ENDDO
            ELSE
!               Calulate Rc for node a (G(3))
                DO gp = 1,el%gp
                    G(ai) = G(ai) &
                    &     -(el%Nxg(a,1,gp)*upgp(1,gp) &
                    &     + el%Nxg(a,2,gp)*upgp(2,gp))*el%Wg(gp)
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
!       Go through surrounding nodes and incorporate into resid for a 
        DO b = 1,el%eNoN
!           Again through dof, now for b
            DO j =1,el%dof
                bi = bi+1
!               Loop through Gauss points to calc integral
                DO gp = 1,el%gp

                ENDDO
            ENDDO
        ENDDO
    ENDIF
ENDDO
ENDDO


END SUBROUTINE lresidns

END MODULE LRESIDNSMOD