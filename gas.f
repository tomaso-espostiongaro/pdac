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
!  In the pdac code there are "max_ngas" predefined gas types,
!  whose constants and parameters are statically compiled
!  with pdac.
!  present_gas logical array is a musk of size "max_ngas"
!           and the i-th element is TRUE the i-th predefined
!           gas type is included in the simulation.
!  gas_type array maps the i-th simulated gas into the array 
!           of gas types statically defined at compile time.
!           i.e.:  we are simulating three gas then:
!           i1 = gas_type(1)
!           i2 = gas_type(2)
!           i3 = gas_type(3)
!           where i1, i2 and i3 are three gas out of the "max_ngas"
!           predefined in pdac.
!  
!
        REAL*8 :: gammaair, gamn, c_joule, c_erg, rgas
        REAL*8 :: tzero, hzerog, hzeros
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE allocate_gas_constants
      USE dimensions
      IMPLICIT NONE
      INTEGER :: ig
!
      ALLOCATE(phij(max_ngas,max_ngas))
      ALLOCATE(ckg(max_ngas)) 
      ALLOCATE(gmw(max_ngas)) 
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
      present_gas = .FALSE.
      gas_type = 0

      RETURN
      END SUBROUTINE allocate_gas_constants
!----------------------------------------------------------------------
      END  MODULE gas_constants
!----------------------------------------------------------------------
