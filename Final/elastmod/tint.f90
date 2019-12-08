MODULE TINTMOD
    USE SOLMOD
IMPLICIT NONE

    REAL(KIND=8) :: alphm, alphf, rhoin, gam, bet

CONTAINS

!=======================================================================
! Initializes the timestepping variables
SUBROUTINE tintinit(rhoi)
    REAL(KIND=8), INTENT(IN) :: rhoi
    rhoin = rhoi
    alphm = (3D0-rhoin)/2D0/(1D0+rhoin)
    alphf = 1D0/(1D0+rhoin)
    gam   = alphf
    bet = 1D0/4D0*(1+alphm-alphf)**2D0

END SUBROUTINE tintinit

!=======================================================================
! Evaluates the current solution variables at alphas
SUBROUTINE toalph(sol)
    TYPE(soltype), INTENT(INOUT) :: sol
    INTEGER i,j

!   Loop over all dof's
    DO i = 1,sol%dof
!       Loop over all nodes
        DO j = 1,sol%np
!           Get d at alphf
            sol%dalphf(i,j) = sol%do(i,j) + &
            &          alphf*(sol%d(i,j)  - sol%do(i,j))
!           Get ddot at alphm 
            sol%ddalphm(i,j) = sol%ddoto(i,j) + &
            &           alphm*(sol%ddot(i,j)  - sol%ddoto(i,j))
!           Get ddot at alphf
            sol%ddalphf(i,j) = sol%ddoto(i,j) + &
            &          alphf*(sol%ddot(i,j)  - sol%ddoto(i,j))
!           Get dddot at alphm 
            sol%dddalphm(i,j) = sol%dddo(i,j) + &
            &           alphm*(sol%ddd(i,j)  - sol%dddo(i,j))
        ENDDO
    ENDDO


END SUBROUTINE toalph

!=======================================================================
SUBROUTINE resid(sol,G)
    TYPE(solType), INTENT(IN) :: sol
    REAL(KIND=8), INTENT(OUT) :: G
    G=sol%d(1,1)
    

END SUBROUTINE resid





END MODULE TINTMOD








