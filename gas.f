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
        INTEGER, ALLOCATABLE :: gas_type(:)
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
      INTEGER :: ig
!
      ALLOCATE(phij(max_ngas,max_ngas))
      ALLOCATE(ckg(max_ngas), gmw(max_ngas)) 
      ALLOCATE(mmugek(max_ngas), mmug(max_ngas), mmugs(max_ngas))
      ALLOCATE(present_gas(max_ngas))
      ALLOCATE(gas_type(ngas))
!
      phij   = 0.D0
      ckg    = 0.D0
      mmug   = 0.D0
      mmugs  = 0.D0
      mmugek = 0.D0
      gmw    = 0.D0
      DO ig = 1, max_ngas
        present_gas(ig) = (ig == default_gas)
      END DO
      gas_type = 0

      RETURN
      END SUBROUTINE allocate_gas_constants
!----------------------------------------------------------------------
      END  MODULE gas_constants
!----------------------------------------------------------------------
