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
REAL(KIND=8), INTENT(OUT) :: G(el%eNoN*el%dof),Gt(el%eNoN*el%dof,el%eNoN*el%dof)
INTEGER :: a,b,gp,i,j,ai,bi, ag
! Subgrid variables
REAL(KIND=8) up(u%dof,el%eNoN), pp(p%dof,el%eNoN), ugp(u%dof,el%gp), &
&            duxgp(u%dof,u%dof,el%gp), dutgp(u%dof,el%gp), pgp(el%gp), &
&            upgp(u%dof,el%gp), ppgp(el%gp), nucgp(el%gp), taumgp(el%gp), &
&            Nxg(el%eNoN,u%dof,el%gp), Ng(el%eNoN,el%gp)   
! Temporary material properties !!!
REAL(KIND=8) :: rho,mu, taum(el%eNoN), nuc(el%eNoN)


!REAL(KIND=8) :: gtemp1,gtemp3,gtemp2
!gtemp1 = 0
!gtemp2 = 0
!gtemp3 = 0

! Put shape functions and their derivs into their own variable
Nxg = el%Nxg
Ng = el%Ng

rho = 1D0
mu = 1D0

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
            ugp(i,gp)   = ugp(i,gp)  + Ng(a,gp)*u%dalphf(i,ag)
            dutgp(i,gp) = dutgp(i,gp)+ Ng(a,gp)*u%ddalphm(i,ag)
!           This term is deriv of ith velocity in jth direction at gp
            DO j = 1,u%dof
                duxgp(i,j,gp) = duxgp(i,j,gp) + Nxg(a,j,gp)*u%dalphf(i,ag)
            ENDDO
        ENDDO
    ENDDO
ENDDO

! Now p as well
DO gp = 1,el%gp
    DO a = 1,el%eNoN
        ag = el%nds(a)
        pgp(gp) = pgp(gp) + Ng(a,gp)*p%d(1,ag)
    ENDDO
ENDDO

! Loop through element nodes to calc resid without subgrid
DO a = 1,el%eNoN
! Loop through dof (for index a)
DO i = 1,el%dof
    ai = ai+1
!   If this isn't a boundary node, go ahead and calculate the resid at a
    IF (el%bnd(a,1,i).eq.0) THEN
!           Loop through Gauss points to calc integral
!           If i = 1, we're dealing with x-mom, 2->y-mom, 3->cont
            IF((i.eq.1).or.(i.eq.2)) THEN
!               Calculate Rm in direction i for node a (G(1:2))
                DO gp = 1,el%gp
                    G(ai) = G(ai) &
!                   Advection/forcing part
                    &     + (Ng(a,gp)*(dutgp(i,gp)+(ugp(1,gp)*duxgp(i,1,gp) + ugp(2,gp)*duxgp(i,2,gp)) - fun(i,gp)) &

                    !&     +(Ng(a,gp)*rho*(dutgp(i,gp) &
                    !&     + ugp(1,gp)*duxgp(i,1,gp) &
                    !&     + ugp(2,gp)*duxgp(i,2,gp) &
                    !&     - fun(i,gp)) &
!                   Viscous part
                    &     + mu/2D0 &
                    &     *(duxgp(i,1,gp)*Nxg(a,1,gp) &
                    &     + duxgp(i,2,gp)*Nxg(a,2,gp) &
                    &     + duxgp(1,i,gp)*Nxg(a,1,gp) &
                    &     + duxgp(2,i,gp)*Nxg(a,2,gp))&
!                   Pressure part
                    &     - Nxg(a,i,gp)*pgp(gp))*el%Wg(gp)
!                   Traction part (0 for now)
                ENDDO
            ELSE
!               Calulate Rc for node a (G(3))
                DO gp = 1,el%gp
                   G(ai) = G(ai) &
                   &     + (Ng(a,gp) &
                   &     * (duxgp(1,1,gp) + duxgp(2,2,gp)))*el%Wg(gp)
                ENDDO
            ENDIF
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
    taum(a) = 0
    taum(a) = 1D-7

    nuc(a) = 1/(el%eG(1,1)+el%eG(2,2))/taum(a)
    nuc(a) = 1D-7

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
            upgp(i,gp) = upgp(i,gp)  + Ng(a,gp)*up(i,a)
        ENDDO
    ENDDO
ENDDO

! Now pp at gp as well
DO gp = 1,el%gp
    DO a = 1,el%eNoN
        ppgp(gp)   = ppgp  (gp) + Ng(a,gp)*pp(1,a)
        taumgp(gp) = taumgp(gp) + Ng(a,gp)*taum(a)
        nucgp(gp)  = nucgp (gp) + Ng(a,gp)*nuc(a)
    ENDDO
ENDDO

ai = 0

! Loop through element nodes to calc resid with subgrid
DO a = 1,el%eNoN
! Loop through dof (for index a)
DO i = 1,el%dof
    ai = ai+1
!   If this isn't a boundary node, go ahead and calculate the resid at a
    IF (el%bnd(a,1,i).eq.0) THEN
!           Loop through Gauss points to calc integral
!           If i = 1, we're dealing with x-mom, 2->y-mom, 3->cont
            IF((i.eq.1).or.(i.eq.2)) THEN
!               Calculate Rm in direction i for node a (G(1:2))
                DO gp = 1,el%gp
                    G(ai) = G(ai) &
!                   First term
                    &     +(Ng(a,gp)*rho*(upgp(1,gp)*duxgp(i,1,gp) &
                    &     + upgp(2,gp)*duxgp(i,2,gp)) &
!                   Skip second u'*u' term, third term
                    &     +(ugp(1,gp)*Nxg(a,1,gp)  &
                    &     + ugp(2,gp)*Nxg(a,2,gp)) &
                    &     * upgp(i,gp)*rho &
!                   Pressure term
                    &     - Nxg(a,i,gp)*ppgp(gp))*el%Wg(gp)
                ENDDO
            ELSE
!               Calulate Rc for node a (G(3))
                DO gp = 1,el%gp
                    G(ai) = G(ai) &
                    &     -(Nxg(a,1,gp)*upgp(1,gp) &
                    &     + Nxg(a,2,gp)*upgp(2,gp))*el%Wg(gp)
                ENDDO
            ENDIF
    ENDIF
ENDDO
ENDDO
ai = 0

! Now for the tangent matrix
DO a = 1,el%eNoN
! Loop through dof of residual (i = 3) = continuity
DO i = 1,el%dof
    bi = 0
    ai = ai+1
!   Check derivatives with same or surrounding nodes (b)
    DO b = 1,el%eNoN
!       Again through dof for parameters (j = 3) = pressure
        DO j =1,el%dof
            bi = bi+1
            DO gp = 1,el%gp
!               4 derivative cases:
!               Momentum residual with respect to velocity
                IF(((i.eq.1).or.(i.eq.2)).and.((j.eq.1).or.(j.eq.2))) THEN
                    Gt(ai,bi) = Gt(ai,bi) &
!                   Second viscous term
                    &         +(mu/2D0*Nxg(a,j,gp)*Nxg(b,i,gp)*alphf*gam*el%dt &
!                   Last RBVMS term with nuc
                    &         + alphf*gam*el%dt*rho*nucgp(gp)*Nxg(a,i,gp)*Nxg(b,j,gp))*el%Wg(gp)
!                   Kronecker delta terms
                    IF(i.eq.j) THEN
                        Gt(ai,bi) = Gt(ai,bi) &
!                       Nonlinear term
                        &         +(alphm*taumgp(gp)*(ugp(1,gp)*Nxg(a,1,gp) + ugp(2,gp)*Nxg(a,2,gp))*rho*Ng(b,gp) &
                        &         + alphm*Ng(a,gp)*Ng(b,gp) &
                        &         + alphf*gam*el%dt*(Ng(a,gp)*rho*(ugp(1,gp)*Nxg(b,1,gp) + ugp(2,gp)*Nxg(b,2,gp))) &
                        &         + alphf*gam*el%dt*((Nxg(a,1,gp)*Nxg(b,1,gp) + Nxg(a,2,gp)*Nxg(b,2,gp))*mu/2D0) &
                        &         + alphf*gam*el%dt*(taumgp(gp)*(ugp(1,gp)*Nxg(a,1,gp) + ugp(2,gp)*Nxg(a,2,gp))*rho*(ugp(1,gp)*Nxg(b,1,gp) + ugp(2,gp)*Nxg(b,2,gp))) &
                        &         )*el%Wg(gp)
 
                        !&         +(rho*Ng(a,gp)*(alphm*Ng(b,gp) &
                        !&         + (ugp(1,gp)*Nxg(b,1,gp) + ugp(2,gp)*Nxg(b,2,gp)) &
                        !&         * alphf*gam*el%dt) &
!                       !First viscous term 
                        !&         + mu/2D0*(Nxg(a,1,gp)*Nxg(b,1,gp) &
                        !&         + Nxg(a,2,gp)*Nxg(b,2,gp))*alphf*gam*el%dt &
!                       !RBVMS term with taum 
                        !&         - (ugp(1,gp)*Nxg(a,1,gp) + ugp(2,gp)*Nxg(a,2,gp)) &
                        !&         * taumgp(gp)*(rho*alphm*Ng(b,gp) &
                        !&         + rho*(ugp(1,gp)*Nxg(b,1,gp) + ugp(2,gp)*Nxg(b,2,gp))&
                        !&         * alphf*gam*el%dt))*el%Wg(gp)
                    ENDIF

!               Momentum residual with respect to pressure    
                ELSEIF (((i.eq.1).or.(i.eq.2)).and.(j.eq.3)) THEN
                    Gt(ai,bi) = Gt(ai,bi) &
                    &         -(Nxg(a,i,gp)*Ng(b,gp) & 
                    &         + taumgp(gp)*(ugp(1,gp)*Nxg(a,1,gp) + ugp(2,gp)*Nxg(a,2,gp))*Nxg(b,i,gp) &
                    &         )*el%Wg(gp)
                    !&         -(Nxg(a,i,gp)*Ng(b,gp) &
                    !&         + taumgp(gp)*(Nxg(a,1,gp)*ugp(1,gp) + Nxg(a,2,gp)*ugp(2,gp)) &
                    !&         * Nxg(b,i,gp))*el%Wg(gp)

!               Continuity residual with respect to velocity    
                ELSEIF ((i.eq.3).and.((j.eq.1).or.(j.eq.2))) THEN
                    Gt(ai,bi) = Gt(ai,bi) &
                               +(alphf*gam*el%dt*(Ng(a,gp)*Nxg(b,j,gp)) &
                    &          + alphf*gam*el%dt*taumgp(gp)*Nxg(a,j,gp)*(ugp(1,gp)*Nxg(b,1,gp) + ugp(2,gp)*Nxg(b,2,gp)) &
                    &          + alphm*taumgp(gp)*Nxg(a,j,gp)*Ng(b,gp) &
                    &          )*el%Wg(gp) 
                    !&         +(Ng(a,gp)*Nxg(b,j,gp)*alphf*gam*el%dt &
                    !&         + Nxg(a,j,gp)*taumgp(gp)*(alphm*Ng(b,gp) &
                    !&         + alphf*gam*el%dt*(ugp(1,gp)*Nxg(b,1,gp) &
                    !&         + ugp(2,gp)*Nxg(b,2,gp))))*el%Wg(gp)

!               Continuity residual with respect to pressure    
                ELSEIF ((i.eq.3).and.(j.eq.3)) THEN
                    Gt(ai,bi) = Gt(ai,bi) &
                    &         +(taumgp(gp)/rho &
                    &         * (Nxg(a,1,gp)*Nxg(b,1,gp) + Nxg(a,2,gp)*Nxg(b,2,gp)))*el%Wg(gp)
                ENDIF
            ENDDO
        ENDDO
    ENDDO
ENDDO
ENDDO


END SUBROUTINE lresidns

END MODULE LRESIDNSMOD