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
        LOGICAL, DIMENSION(:), ALLOCATABLE :: present_gas
!
        REAL*8 :: gammaair, gamn, c_joule, c_erg, rgas
        REAL*8 :: tzero, hzerog, hzeros
        INTEGER :: default_gas
        INTEGER :: mcomp
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE allocate_gas_constants
        USE dimensions
        IMPLICIT NONE
!
      ALLOCATE(phij(ngas,ngas))
      ALLOCATE(ckg(ngas), mmug(ngas), mmugs(ngas), mmugek(ngas), gmw(ngas))
      ALLOCATE(present_gas(ngas))
!
      phij   = 0.D0
      ckg    = 0.D0
      mmug   = 0.D0
      mmugs  = 0.D0
      mmugek = 0.D0
      gmw    = 0.D0
      RETURN
      END SUBROUTINE allocate_gas_constants
!----------------------------------------------------------------------
      END  MODULE gas_constants
!----------------------------------------------------------------------
