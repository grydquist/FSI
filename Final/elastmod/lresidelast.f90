MODULE LRESIDELASTMOD
    USE SOLMOD
    USE ELEMOD
    USE TINTMOD
IMPLICIT NONE

CONTAINS

SUBROUTINE lresidelast(G,Gt,fun,el,y)
TYPE(eltype), INTENT(IN) :: el
REAL(KIND=8), ALLOCATABLE, INTENT(IN) :: fun(:,:)
TYPE(soltype), INTENT(IN) :: y
REAL(KIND=8), INTENT(OUT) :: G(el%eNoN*y%dof),Gt(el%eNoN*y%dof,el%eNoN*y%dof)
INTEGER :: a,b,gp,i,j,ai,bi, ag
! Subgrid variables
REAL(KIND=8) ygp(y%dof,el%gp), dyxgp(y%dof,y%dof,el%gp), ddytgp(y%dof,el%gp), &
&            Nxg(el%eNoN,y%dof,el%gp), Ng(el%eNoN,el%gp), dt
! Temporary material properties !!!
REAL(KIND=8) :: rho,mu, lam, mu2, Em, nu, Ba(3,2),Bb(3,2),D(3,3), tmp1(2,3), tmp2(2,2), tmp3(2)


!REAL(KIND=8) :: gtemp1,gtemp3,gtemp2
!gtemp1 = 0
!gtemp2 = 0
!gtemp3 = 0

! Put shape functions and their derivs into their own variable
Nxg = el%Nxg
Ng = el%Ng
dt = el%dt

Em = 1D8
nu = 0.3
rho = 1D0
lam = nu*Em/(1+nu)/(1-2d0*nu)
mu = Em/(2d0+2d0*nu)
mu2 = 1D0

! Make D matrix
D = 0
D(1,1) = lam+2D0*mu
D(2,2) = D(1,1)
D(1,2) = lam
D(2,1) = lam
D(3,3) = mu
Ba = 0
Bb = 0

Gt = 0
G = 0
ai = 0
ygp = 0
dyxgp = 0
ddytgp = 0

! First let's find y and it's derivatives at the Gauss points
DO i = 1,y%dof
    DO gp = 1,el%gp
        DO a = 1,el%eNoN
            ag = el%nds(a)
            ygp(i,gp)   = ygp(i,gp)   + Ng(a,gp)*y%dalphf(i,ag)
            ddytgp(i,gp)= ddytgp(1,gp)+ Ng(a,gp)*y%dddalphm(i,ag)
!           This term is deriv of ith velocity in jth direction at gp
            DO j = 1,y%dof
                dyxgp(i,j,gp) = dyxgp(i,j,gp) + Nxg(a,j,gp)*y%dalphf(i,ag)
            ENDDO
        ENDDO
    ENDDO
ENDDO

! Loop through element nodes to calc resid
DO a = 1,el%eNoN
! Loop through dof (for index a)
DO i = 1,y%dof
    ai = ai+1
!   If this isn't a boundary node, go ahead and calculate the resid at a
    IF (el%bnd(a,1,i).eq.0) THEN
!       Loop through Gauss points to calc integral. If i = 1, we're dealing with x-disp, 2->y-disp
        DO gp = 1,el%gp
            G(ai) = G(ai) + (0 &
!           Acceleration/ forcing
            &     + Ng(a,gp)*rho*(ddytgp(i,gp) - fun(i,gp))  &
!           Strain part
            &     + lam*Nxg(a,i,gp)*(dyxgp(i,1,gp)*dyxgp(i,1,gp) + dyxgp(i,2,gp)*dyxgp(i,2,gp)) &
            &     + mu*((Nxg(a,1,gp)*dyxgp(i,1,gp) + Nxg(a,2,gp)*dyxgp(i,2,gp))  &
            &     +     (Nxg(a,1,gp)*dyxgp(1,i,gp) + Nxg(a,2,gp)*dyxgp(2,i,gp))) &
!           Traction
!           Integrate
            &     )*el%Wg(gp)*el%J
        ENDDO
    ENDIF
ENDDO
ENDDO

ai = 0

! Now for the tangent matrix
DO a = 1,el%eNoN
! Loop through dof of residual
DO i = 1,y%dof
    bi = 0
    ai = ai+1
!   Check derivatives with same or surrounding nodes (b)
    DO b = 1,el%eNoN
!       Again through dof for parameters
        DO j =1,y%dof
            bi = bi+1
            DO gp = 1,el%gp
!               Strain energy residual with respect to acceleration
                Gt(ai,bi) = Gt(ai,bi) + (0&
!               Non Kronecker axial strain
                &         + alphf*bet*dt*dt*(lam*Nxg(a,i,gp)*Nxg(b,j,gp)) &
!               and shear
                &         + alphf*bet*dt*dt*(mu *Nxg(a,j,gp)*Nxg(b,i,gp)) &
                &         )*el%Wg(gp)*el%J
!               Kronecker delta terms
                IF(i.eq.j) THEN
                    Gt(ai,bi) = Gt(ai,bi)+ (0&
!                   Accel deriv
                    &         + alphm*rho*Ng(a,gp)*Ng(b,gp) &
!                   Strain dot
                    &         + alphf*bet*dt*dt*(mu*(Nxg(a,1,gp)*Nxg(b,1,gp) + Nxg(a,2,gp)*Nxg(b,2,gp))) &
                    &         )*el%Wg(gp)*el%J
                ENDIF
            ENDDO
        ENDDO
    ENDDO
ENDDO
ENDDO
 
! Different way to do the same thing
DO gp = 1,el%gp
    DO a = 1,el%eNoN
        DO b = 1,el%eNoN
!           Make B matrices
            Ba(1,1) = Nxg(a,1,gp)
            Ba(2,2) = Nxg(a,2,gp)
            Ba(3,1) = Nxg(a,2,gp)
            Ba(3,2) = Nxg(a,1,gp)

            Bb(1,1) = Nxg(b,1,gp)
            Bb(2,2) = Nxg(b,2,gp)
            Bb(3,1) = Nxg(b,2,gp)
            Bb(3,2) = Nxg(b,1,gp)

            tmp1 = matmul(transpose(Ba), D)
            tmp2 = matmul(tmp1,Bb)
            tmp2 = tmp2*alphf*bet*dt*dt*el%Wg(gp)*el%J

            DO i = 1,y%dof
                DO j = 1,y%dof
                    !Gt((a-1)*y%dof+i,(b-1)*y%dof+j) = Gt((a-1)*y%dof+i,(b-1)*y%dof+j) + tmp2(i,j)
                ENDDO
            ENDDO            
        ENDDO
    ENDDO
ENDDO



DO gp = 1,el%gp
    DO a = 1,el%eNoN
        tmp3 = 0D0
        DO b = 1,el%eNoN
!           Make B matrices
            Ba(1,1) = Nxg(a,1,gp)
            Ba(2,2) = Nxg(a,2,gp)
            Ba(3,1) = Nxg(a,2,gp)
            Ba(3,2) = Nxg(a,1,gp)

            Bb(1,1) = Nxg(b,1,gp)
            Bb(2,2) = Nxg(b,2,gp)
            Bb(3,1) = Nxg(b,2,gp)
            Bb(3,2) = Nxg(b,1,gp)

            tmp1 = matmul(transpose(Ba), D)
            tmp2 = matmul(tmp1,Bb)
            tmp2 = tmp2*alphf*bet*dt*dt*el%Wg(gp)*el%J            
        ENDDO
        tmp3 = matmul(tmp2,ygp(:,gp))
        !G((a-1)*y%dof+1) = G((a-1)*y%dof+1) + tmp3(1)
        !G((a-1)*y%dof+2) = G((a-1)*y%dof+2) + tmp3(2)
    ENDDO
ENDDO

END SUBROUTINE lresidelast

END MODULE LRESIDELASTMOD