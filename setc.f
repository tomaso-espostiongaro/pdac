!----------------------------------------------------------------------
      SUBROUTINE setc
!----------------------------------------------------------------------
!
      USE dimensions
      USE gas_constants, ONLY: gmw, mmugs, mmugek,     &
     &    gammaair, gamn, c_joule, c_erg, tzero, hzerog, hzeros, rgas
      USE particles_constants, ONLY: particles_constants_set, cmus
      USE reactions, ONLY: h1, h2, h3, h4, h5
      USE roughness_module, ONLY: zrough
      USE turbulence, ONLY: turbulence_setup, iturb
      USE gas_solid_viscosity, ONLY: particle_viscosity
      IMPLICIT NONE
      INTEGER :: k
!
! molecular weight of chemical components (O2,N2,CO2,H2,H2O,Air,SO2)
!
      gmw(1) = 32.D0              ! O2
      gmw(2) = 28.D0              ! N2
      gmw(3) = 44.D0              ! CO2
      gmw(4) = 2.D0               ! H2
      gmw(5) = 18.D0              ! H2O
      gmw(6) = 28.96442D0         ! Air
      gmw(7) = 64.D0              ! SO2
!
! Lennard-Jones potentials (sigma and epsilon/k)
!
      mmugs(1) = 3.467D0          ! O2
      mmugek(1) = 106.7D0
      mmugs(2) = 3.798D0          ! N2
      mmugek(2) = 71.4D0
      mmugs(3) = 3.941D0          ! CO2
      mmugek(3) = 195.2D0
      mmugs(4) = 2.827D0          ! H2
      mmugek(4) = 59.7D0
      mmugs(5) = 2.641D0          ! H2O
      mmugek(5) = 809.1D0
      mmugs(6) = 3.711D0          ! Air
      mmugek(6) = 78.6D0
      mmugs(7) = 4.112D0          ! SO2
      mmugek(7) = 335.4D0
!
! Syamlal's particle-particle interaction coefficients
!
      CALL particles_constants_set
!
! Initialize particle's viscosity
! 
      DO k = 1, nsolid
        particle_viscosity(k,:) = cmus(k)
      END DO
!
! set useful constants 
!
      gammaair = 1.33D0                          ! air adiabatic constant
      gamn     = (gammaair - 1.D0) / gammaair    ! useful constant ! 
      c_joule = 4.186D0                          ! cal per joule
      c_erg   = 4.186D7                          ! cal per erg
      tzero  = 0.D0                              ! costant reference temperature
      hzerog = 0.D0                              ! gas enthalpy at tzero
      hzeros = 0.D0                              ! particles enthalpy at tzero
      rgas  = 8.31432D7                          ! perfect gas constant
      h1 = 0.D0                                  !
      h2 = 0.D0                                  !
      h3 = 0.D0                                  ! reaction enthalpies
      h4 = 0.D0                                  !
      h5 = 0.D0                                  !
!
      IF (iturb .GT. 0) THEN
        CALL turbulence_setup( zrough )
      END IF
!
      RETURN
!----------------------------------------------------------------------
      END SUBROUTINE
!----------------------------------------------------------------------
