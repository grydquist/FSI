MODULE SOLMOD
IMPLICIT NONE

TYPE solType
    INTEGER :: dof,np
!   d(dof,np) (current solution), do (last solution), derivatives
    REAL(KIND=8), ALLOCATABLE :: d(:,:), do(:,:), ddot(:,:),ddoto(:,:)
!   solution/derivs at alpha
    REAL(KIND=8), ALLOCATABLE :: dalphf(:,:), ddalphm(:,:)
!   Gives the boundary conditions for both d and ddot
    REAL(KIND=8), ALLOCATABLE :: bnd(:,:,:), bnddot(:,:,:)
END TYPE solType

INTERFACE solType
    PROCEDURE :: newsol
END INTERFACE solType

CONTAINS

FUNCTION newsol(dof,np,bnd1,bnd2) RESULT(sol)
    INTEGER,INTENT(IN) :: dof,np
    REAL(KIND=8), INTENT(IN) :: bnd1(:,:,:), bnd2(:,:,:)
    TYPE(solType) :: sol
    INTEGER i

    sol%dof = dof
    sol%np  = np
    sol%bnd = bnd1
    sol%bnddot = bnd2
    ALLOCATE(sol%d(dof,np),sol%do(dof,np),&
    & sol%ddoto(dof,np),sol%ddot(dof,np), &
    & sol%dalphf(dof,np),sol%ddalphm(dof,np))
    sol%do = 0
    sol%ddoto = 0
    sol%ddot = 0
    i=1

END FUNCTION newsol


END MODULE SOLMOD