!----------------------------------------------------------------------
      MODULE gas_constants
!----------------------------------------------------------------------
! ... This module contains the definition of the physical properties
! ... of each gas phase, and some constant
!
        IMPLICIT NONE
        SAVE
!
        REAL*8, DIMENSION(:), ALLOCATABLE :: ckg, mmug, mmugs, mmugek, gmw
        REAL*8, DIMENSION(:,:), ALLOCATABLE :: phij
!
        REAL*8 :: gammaair, gamn, c_joule, c_erg, rgas
        REAL*8 :: tzero, hzerog, hzeros
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE bounds_gas_constants
        USE dimensions
        IMPLICIT NONE
!
      ALLOCATE(phij(ngas,ngas))
      ALLOCATE(ckg(ngas), mmug(ngas), mmugs(ngas), mmugek(ngas), gmw(ngas))
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      END  MODULE
