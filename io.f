!----------------------------------------------------------------------
      MODULE io_restart
!----------------------------------------------------------------------
      USE eos_gas, ONLY: gc_bulk_density, gc_molar_fraction,            &
                        gc_mass_fraction, gas_specific_heat
      USE gas_solid_density, ONLY: gas_density, gas_bulk_density,       &
                                 solid_bulk_density
      USE gas_solid_velocity, ONLY: gas_velocity_r, gas_velocity_x,     &
                                    gas_velocity_y, gas_velocity_z
      USE gas_solid_velocity, ONLY: solid_velocity_r, solid_velocity_x, &
                                    solid_velocity_y, solid_velocity_z
      USE gas_solid_temperature, ONLY: gas_temperature, gas_enthalpy,   &
                                     solid_temperature, solid_enthalpy
      USE parallel, ONLY: mpime, root
      USE pressure_epsilon, ONLY: gas_pressure , void_fraction
      USE specific_heat, ONLY: gc_specific_heat, solid_specific_heat
      USE time_parameters, ONLY: time
      USE control_flags, ONLY: job_type
!
        IMPLICIT NONE
        PRIVATE

        PUBLIC :: taperd, tapewr, tapebc

!----------------------------------------------------------------------
        CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE tapewr
!
      USE dimensions
      IMPLICIT NONE
!
      INTEGER :: i, j, k, ijk
      INTEGER :: ig
!
      CHARACTER*13 restart_file
      restart_file = 'pdac.res'
!
      IF( mpime .EQ. root ) THEN
!
        OPEN(UNIT=9,form='unformatted', FILE = restart_file)
!
        WRITE(9) time, nr, nx, ny, nz, nsolid

      END IF
!
! ... store the final values of the main physical variables 
! 
      IF( mpime .EQ. root ) THEN

        WRITE(9) (gas_pressure(ijk),ijk=1,ntot)

      END IF

      IF( mpime .EQ. root ) THEN

        WRITE(9) ((solid_bulk_density(ijk,k),ijk=1,ntot),k=1,nsolid)
        WRITE(9) (gas_enthalpy(ijk),ijk=1,ntot)
        IF (job_type == '2D') THEN
          WRITE(9) (gas_velocity_r(ijk),gas_velocity_z(ijk),ijk=1,ntot)
        ELSE IF (job_type == '3D') THEN
          WRITE(9) (gas_velocity_x(ijk),gas_velocity_y(ijk), &
                  gas_velocity_z(ijk),ijk=1,ntot)
        ELSE
          CALL error('tapewr','Unknown job type',1)
        END IF
        WRITE(9) ((solid_enthalpy(ijk,k),ijk=1,ntot),k=1,nsolid)
        IF (job_type == '2D') THEN
          WRITE(9) ((solid_velocity_r(ijk,k),solid_velocity_z(ijk,k), &
                                                 ijk=1,ntot),k=1,nsolid)
        ELSE IF (job_type == '3D') THEN
          WRITE(9) ((solid_velocity_x(ijk,k),solid_velocity_y(ijk,k), &
                     solid_velocity_z(ijk,k),ijk=1,ntot),k=1,nsolid)
        END IF
        WRITE(9) ((gc_mass_fraction(ig,ijk),ig=1,ngas),ijk=1,ntot)
!
        WRITE(9) (gas_bulk_density(ijk),ijk=1,ntot)
        WRITE(9) (gas_density(ijk),ijk=1,ntot)
        WRITE(9) (void_fraction(ijk),ijk=1,ntot)
        WRITE(9) (gas_temperature(ijk),ijk=1,ntot)
        WRITE(9) ((solid_temperature(ijk,k),ijk=1,ntot),k=1,nsolid)
        WRITE(9) ((gc_bulk_density(ijk,ig),ijk=1,ntot),ig=1,ngas)
        WRITE(9) ((gc_molar_fraction(ig,ijk),ig=1,ngas),ijk=1,ntot)
!
! ... store the final values of the constitutive parameters to be set up
!
        WRITE(9) (gas_specific_heat(ijk),ijk=1,ntot)
        WRITE(9) ((solid_specific_heat(k,ijk),k=1,nsolid),ijk=1,ntot)
        WRITE(9) ((gc_specific_heat(ig,ijk),ig=1,ngas),ijk=1,ntot)
!
        CLOSE(9)
! 
      END IF
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE taperd
!
      USE dimensions
      IMPLICIT NONE
!
      INTEGER :: i, j, k, ijk
      INTEGER :: nr_, nx_, ny_, nz_, nsolid_
      INTEGER :: ig

      IF( mpime .EQ. root ) THEN
!
        OPEN(UNIT=9,form='unformatted',FILE='pdac.res')
        READ(9) time, nr_, nx_, ny_, nz_, nsolid_

      END IF

      CALL bcast_integer(nr_, 1, root)
      CALL bcast_integer(nx_, 1, root)
      CALL bcast_integer(ny_, 1, root)
      CALL bcast_integer(nz_, 1, root)
      CALL bcast_integer(nsolid_, 1, root)

      WRITE(*,*) nr, nx, ny, nz, nr_, nx_, ny_, nz_
      IF( nr_ /= nr ) &
        CALL error(' taperd ',' inconsistent dimension nr ', nr_ )
      IF( nx_ /= nx ) &
        CALL error(' taperd ',' inconsistent dimension nx ', nx_ )
      IF( ny_ /= ny ) &
        CALL error(' taperd ',' inconsistent dimension ny ', ny_ )
      IF( nz_ /= nz ) &
        CALL error(' taperd ',' inconsistent dimension nz ', nz_ )
      IF( nsolid_ /= nsolid ) &
        CALL error(' taperd ',' inconsistent dimension nsolid ', nsolid_ )
      
      
      IF( mpime .EQ. root ) THEN
!
! ... read the values of the main physical variables 
! 
        READ(9) (gas_pressure(ijk),ijk=1,ntot)
        READ(9) ((solid_bulk_density(ijk,k),ijk=1,ntot),k=1,nsolid) 
        READ(9) (gas_enthalpy(ijk),ijk=1,ntot)
        IF (job_type == '2D') THEN
          READ(9) (gas_velocity_r(ijk),gas_velocity_z(ijk),ijk=1,ntot)
        ELSE IF (job_type == '3D') THEN
          READ(9) (gas_velocity_x(ijk),gas_velocity_y(ijk), &
                  gas_velocity_z(ijk),ijk=1,ntot)
        ELSE
          CALL error('taperd','Unknown job type',1)
        END IF
        READ(9) ((solid_enthalpy(ijk,k),ijk=1,ntot),k=1,nsolid)
        IF (job_type == '2D') THEN
          READ(9) ((solid_velocity_r(ijk,k),solid_velocity_z(ijk,k), &
                                                 ijk=1,ntot),k=1,nsolid)
        ELSE IF (job_type == '3D') THEN
          READ(9) ((solid_velocity_x(ijk,k),solid_velocity_y(ijk,k), &
                     solid_velocity_z(ijk,k),ijk=1,ntot), k=1,nsolid)
        END IF
        READ(9) ((gc_mass_fraction(ig,ijk),ig=1,ngas),ijk=1,ntot)
!
        READ(9) (gas_bulk_density(ijk),ijk=1,ntot)
        READ(9) (gas_density(ijk),ijk=1,ntot)
        READ(9) (void_fraction(ijk),ijk=1,ntot)
        READ(9) (gas_temperature(ijk),ijk=1,ntot)
        READ(9) ((solid_temperature(ijk,k),ijk=1,ntot),k=1,nsolid)
        READ(9) ((gc_bulk_density(ijk,ig),ijk=1,ntot),ig=1,ngas)
        READ(9) ((gc_molar_fraction(ig,ijk),ig=1,ngas),ijk=1,ntot)
!
! ... read the values of the constitutive parameters
!
        READ(9) (gas_specific_heat(ijk),ijk=1,ntot)
        READ(9) ((solid_specific_heat(k,ijk),k=1,nsolid),ijk=1,ntot)
        READ(9) ((gc_specific_heat(ig,ijk),ig=1,ngas),ijk=1,ntot)
!
        CLOSE (9)
      END IF
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE tapebc
      USE parallel, ONLY: root
      IMPLICIT NONE
!
       CALL bcast_real(gas_pressure,SIZE(gas_pressure),root)
       CALL bcast_real(solid_bulk_density,SIZE(solid_bulk_density),root)
       CALL bcast_real(gas_enthalpy,SIZE(gas_enthalpy),root)
       IF (job_type == '2D') THEN
         CALL bcast_real(gas_velocity_r,SIZE(gas_velocity_r),root)
         CALL bcast_real(gas_velocity_z,SIZE(gas_velocity_z),root)
       ELSE IF (job_type == '3D') THEN
         CALL bcast_real(gas_velocity_x,SIZE(gas_velocity_x),root)
         CALL bcast_real(gas_velocity_y,SIZE(gas_velocity_y),root)
         CALL bcast_real(gas_velocity_z,SIZE(gas_velocity_z),root)
       END IF
       CALL bcast_real(solid_enthalpy,SIZE(solid_enthalpy),root)
       IF (job_type == '2D') THEN
         CALL bcast_real(solid_velocity_r,SIZE(solid_velocity_r),root)
         CALL bcast_real(solid_velocity_z,SIZE(solid_velocity_z),root)
       ELSE IF (job_type == '3D') THEN
         CALL bcast_real(solid_velocity_x,SIZE(solid_velocity_x),root)
         CALL bcast_real(solid_velocity_y,SIZE(solid_velocity_y),root)
         CALL bcast_real(solid_velocity_z,SIZE(solid_velocity_z),root)
       END IF
       CALL bcast_real(gc_mass_fraction,SIZE(gc_mass_fraction),root)
!
       CALL bcast_real(gas_bulk_density,SIZE(gas_bulk_density),root)
       CALL bcast_real(gas_density,SIZE(gas_density),root)
       CALL bcast_real(void_fraction,SIZE(void_fraction),root)
       CALL bcast_real(gas_temperature,SIZE(gas_temperature),root)
       CALL bcast_real(solid_temperature,SIZE(solid_temperature),root)
       CALL bcast_real(gc_bulk_density,SIZE(gc_bulk_density),root)
       CALL bcast_real(gc_molar_fraction,SIZE(gc_molar_fraction),root)
!
       CALL bcast_real(gas_specific_heat,SIZE(gas_specific_heat),root)
       CALL bcast_real(gc_specific_heat,SIZE(gc_specific_heat),root)
       CALL bcast_real(solid_specific_heat,SIZE(solid_specific_heat),root)
!
      END SUBROUTINE
!----------------------------------------------------------------------
      END MODULE io_restart
!----------------------------------------------------------------------
