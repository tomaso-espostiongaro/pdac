!----------------------------------------------------------------------
      MODULE output_dump
!----------------------------------------------------------------------
      USE dimensions
      IMPLICIT NONE
      INTEGER nfil
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------

      SUBROUTINE outp
!
      USE dimensions, ONLY: nr, nz, nsolid
      USE eos_gas, ONLY: gc_molar_fraction
      USE gas_solid_density, ONLY: solid_bulk_density
      USE gas_solid_velocity, ONLY: gas_velocity_r, gas_velocity_z
      USE gas_solid_velocity, ONLY: solid_velocity_r, solid_velocity_z
      USE gas_solid_temperature, ONLY: gas_temperature, solid_temperature
      USE parallel, ONLY: nproc, mpime, root, group
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: gas_pressure
      USE reactions, ONLY: irex
      USE time_parameters, ONLY: time
!
      IMPLICIT NONE
!
      CHARACTER :: filnam*11
      CHARACTER*4 :: letter
      LOGICAL :: cmprs
!
! this variable is used by the non standard routine system
!*************
      CHARACTER*24 syscom
!*************
!
      INTEGER :: ij2, ij1, k, j
      INTEGER :: ijl
      INTEGER :: kg
!
!      cmprs = .TRUE.
      cmprs = .FALSE.
!
      nfil=nfil+1
      filnam='OUTPUT.'//letter(nfil)

      IF( mpime .EQ. root ) THEN

      OPEN(UNIT=3,FILE=filnam)
!
      WRITE(3,547) time
!
      WRITE(3,548)
      DO j=1,nz
        ij1=1+(nz-j)*nr
        ij2=nr+(nz-j)*nr
        WRITE(3,550)(gas_pressure(ijl),ijl=ij1,ij2)
      END DO
!
      DO k=1,nsolid
        WRITE(3,549) k
        DO j=1,nz
          ij1=1+(nz-j)*nr
          ij2=nr+(nz-j)*nr
          WRITE(3,550)(solid_bulk_density(k,ijl)*inrl(k),ijl=ij1,ij2)
        END DO
      END DO
!
      WRITE(3,552)
      DO j=1,nz
        ij1=1+(nz-j)*nr
        ij2=nr+(nz-j)*nr
        WRITE(3,550)(gas_velocity_r(ijl),ijl=ij1,ij2)
      END DO
!
      WRITE(3,553)
      DO j=1,nz
        ij1=1+(nz-j)*nr
        ij2=nr+(nz-j)*nr
        WRITE(3,550)(gas_velocity_z(ijl),ijl=ij1,ij2)
      END DO
!
      IF(irex.GE.0) THEN
        WRITE(3,560)
        DO j=1,nz
          ij1=1+(nz-j)*nr
          ij2=nr+(nz-j)*nr
          WRITE(3,550)(gas_temperature(ijl),ijl=ij1,ij2)
        END DO
!
!pdac------------
! define loop parameters
        DO 343 kg=5,5
!pdac------------
          WRITE(3,562) kg
          DO j=1,nz
            ij1=1+(nz-j)*nr
            ij2=nr+(nz-j)*nr
            WRITE(3,550)(gc_molar_fraction(kg,ijl),ijl=ij1,ij2)
           END DO
 343    CONTINUE
      ENDIF
!
      DO k=1,nsolid
!
        WRITE(3,556) k
        DO j=1,nz
          ij1=1+(nz-j)*nr 
          ij2=nr+(nz-j)*nr
          WRITE(3,550)(solid_velocity_r(k,ijl),ijl=ij1,ij2)
        END DO
!
        WRITE(3,557) k
        DO  j=1,nz
          ij1=1+(nz-j)*nr
          ij2=nr+(nz-j)*nr
          WRITE(3,550)(solid_velocity_z(k,ijl),ijl=ij1,ij2)
        END DO
!
        IF(irex.GE.0) THEN
          WRITE(3,561) k
          DO j=1,nz
            ij1=1+(nz-j)*nr
            ij2=nr+(nz-j)*nr
            WRITE(3,550)(solid_temperature(k,ijl),ijl=ij1,ij2)
          END DO
        ENDIF
!
      END DO
      CLOSE (3)
      END IF
!
      IF (cmprs) THEN
      IF( mpime .EQ. root ) THEN
        WRITE (syscom,4444) nfil,char(0)
! ....   IBM 
!        CALL system(syscom)
 4444   FORMAT ('compress -f OUTPUT.',i4.4,a1)
      END IF
      END IF
!
      RETURN

 547  FORMAT(1x,///,1x,'@@@ TIME = ',G11.4)
 548  FORMAT(1x,//,1x,'P',/)
 549  FORMAT(1x,//,1x,'EPS',i1,/)
 550  FORMAT(1x,10(1x,g12.6))
 552  FORMAT(1x,//,1x,'UG',/)
 553  FORMAT(1x,//,1x,'VG',/)
 556  FORMAT(1x,//,1x,'UP',I1,/)
 557  FORMAT(1x,//,1x,'VP',I1,/)
 560  FORMAT(1x,//,1x,'TG',/)
 561  FORMAT(1x,//,1x,'TP',I1,/)
 562  FORMAT(1x,//,1x,'XGC',I1,/)
!
      END SUBROUTINE

!----------------------------------------------------------------------
      END MODULE output_dump
!----------------------------------------------------------------------
