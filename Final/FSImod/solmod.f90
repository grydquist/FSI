MODULE SOLMOD
IMPLICIT NONE

TYPE solType
    INTEGER :: dof,np
!   d(dof,np) (current solution), do (last solution), derivatives
    REAL(KIND=8), ALLOCATABLE :: d(:,:), do(:,:), ddot(:,:),ddoto(:,:), ddd(:,:), dddo(:,:)
!   solution/derivs at alpha
    REAL(KIND=8), ALLOCATABLE :: dalphf(:,:), ddalphm(:,:), dddalphm(:,:), ddalphf(:,:)
!   Gives the boundary conditions for both d and ddot
    INTEGER, ALLOCATABLE :: bnd(:,:,:), bnddot(:,:,:), bnddd(:,:,:)
END TYPE solType

INTERFACE solType
    PROCEDURE :: newsol
END INTERFACE solType

CONTAINS

FUNCTION newsol(dof,np,bnd1,bnd2) RESULT(sol)
    INTEGER,INTENT(IN) :: dof,np
    INTEGER, ALLOCATABLE, INTENT(IN) :: bnd1(:,:,:), bnd2(:,:,:)
    TYPE(solType) :: sol
    INTEGER i,j,shp(3)

    sol%dof = dof
    sol%np  = np
    shp =  shape(bnd1)
!   Initialize stuff
    ALLOCATE(sol%bnd(shp(1),shp(2),shp(3)),sol%bnddot(shp(1),shp(2),shp(3)),sol%bnddd(shp(1),shp(2),shp(3)))
    sol%bnd = bnd1
    sol%bnddot = bnd2
    sol%bnddd  = bnd2
    ALLOCATE(sol%d(dof,np),sol%do(dof,np),&
    & sol%ddoto(dof,np),sol%ddot(dof,np), &
    & sol%dalphf(dof,np),sol%ddalphm(dof,np),&
    & sol%ddalphf(dof,np),sol%dddalphm(dof,np),&
    & sol%ddd(dof,np), sol%dddo(dof,np))

    DO i =1,np
        DO j = 1,dof
!           Set Dirichlet BC's
            IF (sol%bnd(i,1,j) .eq. 1) THEN
                sol%d(j,i) = sol%bnd(i,2,j)
                sol%do(j,i)= sol%bnd(i,2,j)
                !IF (sol%bnd(i,2,j).eq.1) THEN!!!!!!!!!!!!!!!!!!1
                !    sol%d(j,i) = 0.001D0
                !    sol%do(j,i)= 0.001D0
                !ENDIF!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ELSE
                sol%d(j,i) = 0
                sol%do(j,i)= 0
            ENDIF

            IF (sol%bnddot(i,1,j) .eq. 1) THEN
                sol%ddot(j,i) = sol%bnddot(i,2,j)
                sol%ddoto(j,i)= sol%bnddot(i,2,j)
            ELSE
                sol%ddot(j,i) = 0
                sol%ddoto(j,i)= 0
            ENDIF

            IF (sol%bnddd(i,1,j) .eq. 1) THEN
                sol%ddd(j,i) = sol%bnddd(i,2,j)
                sol%dddo(j,i)= sol%bnddd(i,2,j)
            ELSE
                sol%ddd(j,i) = 0
                sol%dddo(j,i)= 0
            ENDIF
        ENDDO
    ENDDO

END FUNCTION newsol


END MODULE SOLMOD