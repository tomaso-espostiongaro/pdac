!----------------------------------------------------------------------
      MODULE gas_solid_density
!----------------------------------------------------------------------
      IMPLICIT NONE
!
      REAL*8, DIMENSION(:), ALLOCATABLE :: rgp, rgpn, rog
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: rlk, rlkn
!
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE allocate_density
      USE dimensions
      USE domain_mapping, ONLY: ncint, ncdom
!
      ALLOCATE(rgp(ncdom), rgpn(ncdom), rog(ncint))
      ALLOCATE(rlk(ncdom,nsolid), rlkn(ncdom,nsolid))
      rgp = 0.D0
      rgpn = 0.D0
      rog = 0.D0
      rlk = 0.D0
      rlkn = 0.D0
      
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      END MODULE
!----------------------------------------------------------------------
