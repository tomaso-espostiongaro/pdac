!----------------------------------------------------------------------
      MODULE output_dump
!----------------------------------------------------------------------
      USE kinds
      USE dimensions
      USE control_flags, ONLY: nfil
      IMPLICIT NONE
      SAVE
      LOGICAL :: formatted_output
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE shock_tube_out
!
      USE dimensions, ONLY: nx, nz, nsolid
      USE domain_decomposition, ONLY: fl_l, ncint
      USE gas_constants, ONLY: present_gas, default_gas, gammaair
      USE gas_solid_density, ONLY: rog
      USE gas_solid_velocity, ONLY: ug, wg
      USE gas_solid_temperature, ONLY: sieg
      USE parallel, ONLY: nproc, mpime, root, group
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: p
      USE time_parameters, ONLY: time
      USE control_flags, ONLY: job_type
      USE domain_decomposition, ONLY: data_collect, data_distribute
      USE dimensions

!
      IMPLICIT NONE
!
      CHARACTER :: filnam*15
      CHARACTER*4 :: letter
!
      INTEGER :: i,j,ijk
      INTEGER :: ig
      REAL*8 :: energy
!
      filnam='shtube.'//letter(nfil)

      IF( mpime .EQ. root ) THEN

      OPEN(UNIT=11,FILE=filnam)
!
      DO ijk = 1, ncint
      IF (fl_l(ijk) == 1) THEN
        energy = p(ijk)/(rog(ijk)*(gammaair - 1.D0))
        WRITE(11,550)rog(ijk),wg(ijk),p(ijk),energy
      END IF
      END DO
!
      CLOSE (11)
      END IF
!
      RETURN

 550  FORMAT(1x,10(1x,g12.6))
!
      END SUBROUTINE shock_tube_out
!----------------------------------------------------------------------

      SUBROUTINE outp
!
      USE dimensions, ONLY: nsolid
      USE eos_gas, ONLY: xgc
      USE gas_constants, ONLY: present_gas, default_gas
      USE gas_solid_density, ONLY: rlk
      USE gas_solid_velocity, ONLY: ug, vg, wg
      USE gas_solid_velocity, ONLY: us, vs, ws
      USE gas_solid_temperature, ONLY: tg, ts
      USE io_restart, ONLY: write_array
      USE parallel, ONLY: nproc, mpime, root, group
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: p
      USE time_parameters, ONLY: time
      USE turbulence_model, ONLY: modturbo
      USE control_flags, ONLY: job_type
!
      IMPLICIT NONE
!
      CHARACTER :: filnam*11
      CHARACTER*4 :: letter
      LOGICAL :: lform
!
      INTEGER :: ig,is
      REAL*8, ALLOCATABLE :: otmp(:)
!
    
      nfil=nfil+1
      filnam='output.'//letter(nfil)
      lform = formatted_output

      IF( mpime == root ) THEN

        IF (lform) THEN
          OPEN(UNIT=12,FILE=filnam)
          WRITE(12,*) time
        ELSE 
          OPEN(UNIT=12,FORM='UNFORMATTED',FILE=filnam)
          WRITE(12) REAL(time,4)
        END IF

        WRITE(6,fmt="('  from outp: writing file ',A20)") filnam
 
      END IF
!
      CALL write_array( 12, p, sgl, lform )  ! gas_pressure

      IF (job_type == '2D') THEN
        CALL write_array( 12, ug, sgl, lform ) ! gas_velocity_r
        CALL write_array( 12, wg, sgl, lform ) ! gas_velocity_z
      ELSE IF (job_type == '3D') THEN
        CALL write_array( 12, ug, sgl, lform ) ! gas_velocity_x
        CALL write_array( 12, vg, sgl, lform ) ! gas_velocity_y
        CALL write_array( 12, wg, sgl, lform ) ! gas_velocity_z
      ELSE
        CALL error('outp_','Unknown job type',1)
      END IF

      CALL write_array( 12, tg, sgl, lform )  ! gas_temperature
!
      ALLOCATE( otmp( SIZE( xgc, 2 ) ) )
      DO ig=1,ngas
        IF( present_gas(ig) .AND. (ig /= default_gas) ) THEN
          otmp = xgc(ig,:)
          CALL write_array( 12, otmp, sgl, lform )  ! gc_molar_fraction
        END IF
      END DO
      DEALLOCATE( otmp )
!
      ALLOCATE( otmp( SIZE( rlk, 1 ) ) )

      DO is = 1, nsolid
        otmp = rlk(:,is)*inrl(is)
        CALL write_array( 12, otmp, sgl, lform )  ! solid_bulk_density
        IF (job_type == '2D') THEN
          CALL write_array( 12, us(:,is), sgl, lform )  ! solid_velocity_r
          CALL write_array( 12, ws(:,is), sgl, lform )  ! solid_velocity_z
        ELSE IF (job_type == '3D') THEN
          CALL write_array( 12, us(:,is), sgl, lform )  ! solid_velocity_x
          CALL write_array( 12, vs(:,is), sgl, lform )  ! solid_velocity_y
          CALL write_array( 12, ws(:,is), sgl, lform )  ! solid_velocity_z
        END IF
        CALL write_array( 12, ts(:,is), sgl, lform )  ! solid_temperature
      END DO

      DEALLOCATE( otmp )

      IF( mpime == root ) THEN
        CLOSE (12)
      END IF
!
      RETURN
      END SUBROUTINE outp
!-----------------------------------------------------------------------
      END MODULE output_dump
!----------------------------------------------------------------------
