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
REAL(KIND=8) ugp(u%dof,el%gp), duxgp(u%dof,u%dof,el%gp), dutgp(u%dof,el%gp), pgp(el%gp), &
&            upgp(u%dof,el%gp), ppgp(el%gp), nucgp(el%gp), taumgp(el%gp), &
&            Nxg(el%eNoN,u%dof,el%gp), Ng(el%eNoN,el%gp), Rm(u%dof,el%gp),Rc(el%gp),&
&            dpxgp(u%dof,el%gp)
! Temporary material properties !!!
REAL(KIND=8) :: rho,mu


!REAL(KIND=8) :: gtemp1,gtemp3,gtemp2
!gtemp1 = 0
!gtemp2 = 0
!gtemp3 = 0

! Put shape functions and their derivs into their own variable
Nxg = el%Nxg
Ng = el%Ng

rho = 1.18D-3!1D0
mu = 1.82D-4!1D0

Gt = 0
G = 0
ai = 0
ugp = 0
duxgp = 0
dutgp = 0
dpxgp = 0
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
        DO j = 1,u%dof
            dpxgp(j,gp) = dpxgp(j,gp) + Nxg(a,j,gp)*p%d(1,ag)
        ENDDO
    ENDDO
ENDDO

upgp = 0
ppgp = 0

! Calculate Rm at gauss points for subgrid stuff
DO i = 1,u%dof
    DO gp  = 1,el%gp
        Rm(i,gp) = rho*(dutgp(i,gp) + ugp(1,gp)*duxgp(i,1,gp) + ugp(2,gp)*duxgp(i,2,gp) - fun(i,gp)) & 
        &        - dpxgp(i,gp)
    ENDDO
ENDDO

! Calculate Rc at gauss points
DO gp = 1,el%gp
    Rc(gp) = duxgp(1,1,gp) + duxgp(2,2,gp)
ENDDO

! Calculate subgrid variables
DO gp = 1,el%gp
    taumgp(gp) =(ugp(1,gp)*el%eG(1,1)*ugp(1,gp) &
    &          + ugp(1,gp)*el%eG(1,2)*ugp(2,gp) &
    &          + ugp(2,gp)*el%eG(2,1)*ugp(1,gp) &
    &          + ugp(2,gp)*el%eG(2,2)*ugp(2,gp) &
    &          + 4D0*(mu/rho)**2 &
    &          * (el%eG(1,1)*el%eG(1,1) + el%eG(1,2)*el%eG(1,2) &
    &          +  el%eG(2,1)*el%eG(2,1) + el%eG(2,2)*el%eG(2,2))&
    &          + 4D0/el%dt/el%dt)**(-0.5D0)

    nucgp(gp) = 1/(el%eG(1,1)+el%eG(2,2))/taumgp(gp)

!   Add subgrid values back in
    DO i = 1,u%dof
        upgp(i,gp) = -taumgp(gp)*Rm(i,gp)/rho
        ugp(i,gp) = ugp(i,gp) + upgp(i,gp)
    ENDDO
    ppgp(gp) = -rho*nucgp(gp)*Rc(gp)
    pgp(gp) = pgp(gp) + ppgp(gp)
ENDDO

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
!                   Advection/forcing part
                    &     +(Ng(a,gp)*rho*(dutgp(i,gp) &
                    &     + ugp(1,gp)*duxgp(i,1,gp) &
                    &     + ugp(2,gp)*duxgp(i,2,gp) &
                    &     - fun(i,gp)) &
!                   Viscous part
                    &     + mu/2D0 &
                    &     *(duxgp(i,1,gp)*Nxg(a,1,gp) &
                    &     + duxgp(i,2,gp)*Nxg(a,2,gp) &
                    &     + duxgp(1,i,gp)*Nxg(a,1,gp) &
                    &     + duxgp(2,i,gp)*Nxg(a,2,gp))&
!                   Pressure part
                    &     - Nxg(a,i,gp)*pgp(gp))*el%Wg(gp)*el%J
!                   Traction part (0 for now)
                ENDDO
            ELSE
!               Calulate Rc for node a (G(3))
                DO gp = 1,el%gp
                   G(ai) = G(ai) &
                   &     + (Ng(a,gp) &
                   &     * (duxgp(1,1,gp) + duxgp(2,2,gp)))*el%Wg(gp)*el%J
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
                    Gt(ai,bi) = Gt(ai,bi) + (0&
!                   Second viscous term
                    &         +  mu/2D0*Nxg(a,j,gp)*Nxg(b,i,gp)*alphf*gam*el%dt &
!                   Last RBVMS term with nuc
                    &         + alphf*gam*el%dt*rho*nucgp(gp)*Nxg(a,i,gp)*Nxg(b,j,gp)&
                    &          )*el%Wg(gp)*el%J
!                   Kronecker delta terms
                    IF(i.eq.j) THEN
                        Gt(ai,bi) = Gt(ai,bi)+ (0&
                        &         + alphm*taumgp(gp)*(ugp(1,gp)*Nxg(a,1,gp) + ugp(2,gp)*Nxg(a,2,gp))*rho*Ng(b,gp) &
                        &         + alphm*rho*Ng(a,gp)*Ng(b,gp) &
                        &         + (Ng(a,gp)*rho*(ugp(1,gp)*Nxg(b,1,gp) + ugp(2,gp)*Nxg(b,2,gp)))*alphf*gam*el%dt &
                        &         + alphf*gam*el%dt*((Nxg(a,1,gp)*Nxg(b,1,gp) + Nxg(a,2,gp)*Nxg(b,2,gp))*mu*2D0) &
                        &         + alphf*gam*el%dt*(taumgp(gp)*(ugp(1,gp)*Nxg(a,1,gp) + ugp(2,gp)*Nxg(a,2,gp))*rho*(ugp(1,gp)*Nxg(b,1,gp) + ugp(2,gp)*Nxg(b,2,gp))) &
                        &         )*el%Wg(gp)*el%j
                    ENDIF

!               Momentum residual with respect to pressure    
                ELSEIF (((i.eq.1).or.(i.eq.2)).and.(j.eq.3)) THEN
                    Gt(ai,bi) = Gt(ai,bi) + (0 &
                    &         - Nxg(a,i,gp)*Ng(b,gp) & 
                    &         + taumgp(gp)*(ugp(1,gp)*Nxg(a,1,gp) + ugp(2,gp)*Nxg(a,2,gp))*Nxg(b,i,gp) &
                    &         )*el%Wg(gp)*el%J

!               Continuity residual with respect to velocity    
                ELSEIF ((i.eq.3).and.((j.eq.1).or.(j.eq.2))) THEN
                    Gt(ai,bi) = Gt(ai,bi) + (0 &
                               + alphf*gam*el%dt*(Ng(a,gp)*Nxg(b,j,gp)) &
                    &          + alphf*gam*el%dt*taumgp(gp)*Nxg(a,j,gp)*(ugp(1,gp)*Nxg(b,1,gp) + ugp(2,gp)*Nxg(b,2,gp)) &
                    &          + alphm*taumgp(gp)*Nxg(a,j,gp)*Ng(b,gp) &
                    &          )*el%Wg(gp)*el%J

!               Continuity residual with respect to pressure    
                ELSEIF ((i.eq.3).and.(j.eq.3)) THEN
                    Gt(ai,bi) = Gt(ai,bi) + (0 &
                    &         + taumgp(gp)/rho &
                    &         * (Nxg(a,1,gp)*Nxg(b,1,gp) + Nxg(a,2,gp)*Nxg(b,2,gp))&
                    &         )*el%Wg(gp)*el%J
                ENDIF
            ENDDO
        ENDDO
    ENDDO
ENDDO
ENDDO


END SUBROUTINE lresidns

END MODULE LRESIDNSMOD