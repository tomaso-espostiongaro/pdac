!----------------------------------------------------------------------
      MODULE io_restart
!----------------------------------------------------------------------
      USE eos_gas, ONLY: gc_bulk_density, gc_molar_fraction,           &
                        gc_mass_fraction, gas_heat_capacity
      USE gas_solid_density, ONLY: gas_density, gas_bulk_density,      &
                                 solid_bulk_density
      USE gas_solid_velocity, ONLY: gas_velocity_r, gas_velocity_z,    &
                                  solid_velocity_r, solid_velocity_z
      USE gas_solid_temperature, ONLY: gas_temperature, gas_enthalpy,  &
                                     solid_temperature, solid_enthalpy
      USE parallel, ONLY: mpime, root
      USE pressure_epsilon, ONLY: gas_pressure, void_fraction
      USE heat_capacity, ONLY: gc_heat_capacity, solid_heat_capacity
      USE time_parameters, ONLY: time

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
      INTEGER :: i, j, k, nrnz, ij
      INTEGER :: kg
!
      CHARACTER*27 sysdump
      CHARACTER*22 sysdumpc
!
      IF( mpime .EQ. root ) THEN
!
! ... non standard fortran CALL (can be omitted)
!*************
          WRITE (sysdump,4445) char(0)
! ... IBM 
!          CALL system(sysdump)
 4445     FORMAT ('mv pdac.res.Z pdac.bak.Z',a1)
!*************
!
        OPEN(UNIT=9,form='unformatted',FILE='pdac.res')
!
        WRITE(9) time,nr,nz,nsolid
        nrnz=nr*nz
!
! ... store the final values of the main physical variables 
! 
        WRITE(9) (gas_pressure(ij),ij=1,nrnz)
        WRITE(9) (void_fraction(ij),ij=1,nrnz)
        WRITE(9) (gas_density(ij),ij=1,nrnz)
        WRITE(9) (gas_bulk_density(ij),ij=1,nrnz)
        WRITE(9) ((solid_bulk_density(k,ij),k=1,nsolid),ij=1,nrnz) 
        WRITE(9) (gas_velocity_r(ij),gas_velocity_z(ij),ij=1,nrnz)
        WRITE(9) ((solid_velocity_r(k,ij),solid_velocity_z(k,ij),k=1,nsolid),ij=1,nrnz)
        WRITE(9) (gas_enthalpy(ij),gas_temperature(ij),ij=1,nrnz)
        WRITE(9) ((solid_enthalpy(k,ij),solid_temperature(k,ij),k=1,nsolid),ij=1,nrnz)
        WRITE(9) ((gc_mass_fraction(kg,ij),gc_molar_fraction(kg,ij),kg=1,ngas),ij=1,nrnz)
        WRITE(9) ((gc_bulk_density(kg,ij),kg=1,ngas),ij=1,nrnz)
!
! ... store the final values of the constitutive parameters to be set up
!
        WRITE(9) (gas_heat_capacity(ij),ij=1,nrnz)
        WRITE(9) ((gc_heat_capacity(k,ij),k=1,nsolid),ij=1,nrnz)
        WRITE(9) ((solid_heat_capacity(kg,ij),kg=1,ngas),ij=1,nrnz)
!
        CLOSE(9)
! 
! ... non standard fortran CALL (can be omitted)
!*************
        WRITE (sysdumpc,4446) char(0)
! ... IBM 
!        CALL system(sysdumpc)
 4446   FORMAT ('compress -f pdac.res',a1)
!*************
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
      CHARACTER*23 sysdumpu
      INTEGER :: i, j, k, nrnz, ij
      INTEGER :: kg
      IF( mpime .EQ. root ) THEN
!
! non standard fortran CALL (can be omitted)
!*************
         WRITE (sysdumpu,4447) char(0)
 4447    FORMAT ('uncompress pdac.res.Z',a1)
! ... IBM 
!         CALL system(sysdumpu)
!*************
!
        OPEN(UNIT=9,form='unformatted',FILE='pdac.res')
!  
        READ(9) time,nr,nz,nsolid
        nrnz=nr*nz
!
! ... read the final values of the main physical variables 
! 
        READ(9) (gas_pressure(ij),ij=1,nrnz)
        READ(9) (void_fraction(ij),ij=1,nrnz)
        READ(9) (gas_density(ij),ij=1,nrnz)
        READ(9) (gas_bulk_density(ij),ij=1,nrnz)
        READ(9) ((solid_bulk_density(k,ij),k=1,nsolid),ij=1,nrnz) 
        READ(9) (gas_velocity_r(ij),gas_velocity_z(ij),ij=1,nrnz)
        READ(9) ((solid_velocity_r(k,ij),solid_velocity_z(k,ij),k=1,nsolid),ij=1,nrnz)
        READ(9) (gas_enthalpy(ij),gas_temperature(ij),ij=1,nrnz)
        READ(9) ((solid_enthalpy(k,ij),solid_temperature(k,ij),k=1,nsolid),ij=1,nrnz)
        READ(9) ((gc_mass_fraction(kg,ij),gc_molar_fraction(kg,ij),kg=1,ngas),ij=1,nrnz)
        READ(9) ((gc_bulk_density(kg,ij),kg=1,ngas),ij=1,nrnz)
!
! ... read the final values of the constitutive parameters set up
!
        READ(9) (gas_heat_capacity(ij),ij=1,nrnz)
        READ(9) ((gc_heat_capacity(k,ij),k=1,nsolid),ij=1,nrnz)
        READ(9) ((solid_heat_capacity(kg,ij),kg=1,ngas),ij=1,nrnz)
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
       CALL bcast_real(void_fraction,SIZE(void_fraction),root)
       CALL bcast_real(gas_density,SIZE(gas_density),root)
       CALL bcast_real(gas_bulk_density,SIZE(gas_bulk_density),root)
       CALL bcast_real(solid_bulk_density,SIZE(solid_bulk_density),root)
       CALL bcast_real(gas_velocity_r,SIZE(gas_velocity_r),root)
       CALL bcast_real(gas_velocity_z,SIZE(gas_velocity_z),root)
       CALL bcast_real(solid_velocity_r,SIZE(solid_velocity_r),root)
       CALL bcast_real(solid_velocity_z,SIZE(solid_velocity_z),root)
       CALL bcast_real(gas_enthalpy,SIZE(gas_enthalpy),root)
       CALL bcast_real(gas_temperature,SIZE(gas_temperature),root)
       CALL bcast_real(solid_enthalpy,SIZE(solid_enthalpy),root)
       CALL bcast_real(solid_temperature,SIZE(solid_temperature),root)
       CALL bcast_real(gc_mass_fraction,SIZE(gc_mass_fraction),root)
       CALL bcast_real(gc_molar_fraction,SIZE(gc_molar_fraction),root)
       CALL bcast_real(gc_bulk_density,SIZE(gc_bulk_density),root)
!
       CALL bcast_real(gas_heat_capacity,SIZE(gas_heat_capacity),root)
       CALL bcast_real(gc_heat_capacity,SIZE(gc_heat_capacity),root)
       CALL bcast_real(solid_heat_capacity,SIZE(solid_heat_capacity),root)
!
      END SUBROUTINE
!----------------------------------------------------------------------
      END MODULE io_restart
