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
      USE gas_constants, ONLY: present_gas, default_gas
      USE gas_solid_density, ONLY: solid_bulk_density
      USE gas_solid_velocity, ONLY: gas_velocity_r, gas_velocity_z
      USE gas_solid_velocity, ONLY: solid_velocity_r, solid_velocity_z
      USE gas_solid_temperature, ONLY: gas_temperature, solid_temperature
      USE parallel, ONLY: nproc, mpime, root, group
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: gas_pressure
      USE time_parameters, ONLY: time
      USE turbulence, ONLY: smag_coeff, modturbo
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
      INTEGER :: ij2, ij1, is, j
      INTEGER :: ijl
      INTEGER :: ig
!
!      cmprs = .TRUE.
      cmprs = .FALSE.
!
! ... MODIFICARE_X3D
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
      DO is=1,nsolid
        WRITE(3,549)is 
        DO j=1,nz
          ij1=1+(nz-j)*nr
          ij2=nr+(nz-j)*nr
          WRITE(3,550)(solid_bulk_density(is,ijl)*inrl(is),ijl=ij1,ij2)
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
      WRITE(3,560)
      DO j=1,nz
        ij1=1+(nz-j)*nr
        ij2=nr+(nz-j)*nr
        WRITE(3,550)(gas_temperature(ijl),ijl=ij1,ij2)
      END DO
!
      DO ig=1,ngas
        IF( present_gas(ig) .AND. (ig /= default_gas) ) THEN
          WRITE(3,562) ig
          DO j=1,nz
            ij1=1+(nz-j)*nr
            ij2=nr+(nz-j)*nr
            WRITE(3,550)(gc_molar_fraction(ig,ijl),ijl=ij1,ij2)
          END DO
        END IF
      END DO
!
      DO is=1,nsolid
!
        WRITE(3,556)is 
        DO j=1,nz
          ij1=1+(nz-j)*nr 
          ij2=nr+(nz-j)*nr
          WRITE(3,550)(solid_velocity_r(is,ijl),ijl=ij1,ij2)
        END DO
!
        WRITE(3,557)is 
        DO  j=1,nz
          ij1=1+(nz-j)*nr
          ij2=nr+(nz-j)*nr
          WRITE(3,550)(solid_velocity_z(is,ijl),ijl=ij1,ij2)
        END DO
!
        WRITE(3,561)is 
        DO j=1,nz
          ij1=1+(nz-j)*nr
          ij2=nr+(nz-j)*nr
          WRITE(3,550)(solid_temperature(is,ijl),ijl=ij1,ij2)
        END DO
!
      END DO

      IF (modturbo > 1) THEN
        WRITE(3,555)
        DO j=1,nz
          ij1=1+(nz-j)*nr
          ij2=nr+(nz-j)*nr
          WRITE(3,550)(smag_coeff(ijl),ijl=ij1,ij2)
        END DO
      END IF
!
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
! 550  FORMAT(1x,10(1x,g14.6e3))
 550  FORMAT(1x,10(1x,g12.6))
 552  FORMAT(1x,//,1x,'UG',/)
 553  FORMAT(1x,//,1x,'VG',/)
 556  FORMAT(1x,//,1x,'UP',I1,/)
 557  FORMAT(1x,//,1x,'VP',I1,/)
 560  FORMAT(1x,//,1x,'TG',/)
 561  FORMAT(1x,//,1x,'TP',I1,/)
 562  FORMAT(1x,//,1x,'XGC',I1,/)

 555  FORMAT(1x,//,1x,'CDYN',/)
!
      END SUBROUTINE outp

!----------------------------------------------------------------------
     SUBROUTINE recover
!
      USE dimensions, ONLY: nr, nz, nsolid
      USE eos_gas, ONLY: gc_molar_fraction
      USE gas_constants, ONLY: present_gas, default_gas
      USE gas_solid_density, ONLY: solid_bulk_density
      USE gas_solid_velocity, ONLY: gas_velocity_r, gas_velocity_z
      USE gas_solid_velocity, ONLY: solid_velocity_r, solid_velocity_z
      USE gas_solid_temperature, ONLY: gas_temperature, solid_temperature
      USE parallel, ONLY: nproc, mpime, root, group
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: gas_pressure
      USE time_parameters, ONLY: time
      USE turbulence, ONLY: smag_coeff, modturbo
!
      IMPLICIT NONE
!
      CHARACTER :: filnam*11
      CHARACTER*4 :: letter
!
      INTEGER :: ij2, ij1, is, j
      INTEGER :: ijl
      INTEGER :: ig
!
! ... MODIFICARE_X3D
!
      filnam='OUTPUT.'//letter(nfil)

      IF( mpime .EQ. root ) THEN

        OPEN(UNIT=3,FILE=filnam)
!
        READ(3,647) time
!
        READ(3,640)
        DO j=1,nz
          ij1=1+(nz-j)*nr
          ij2=nr+(nz-j)*nr
          READ(3,650)(gas_pressure(ijl),ijl=ij1,ij2)
        END DO
!
        DO is=1,nsolid
          READ(3,640)
          DO j=1,nz
            ij1=1+(nz-j)*nr
            ij2=nr+(nz-j)*nr
            READ(3,650)(solid_bulk_density(is,ijl),ijl=ij1,ij2)
          END DO
          solid_bulk_density(is,:) = solid_bulk_density(is,:)*rl(is)
        END DO
!
        READ(3,640)
        DO j=1,nz
          ij1=1+(nz-j)*nr
          ij2=nr+(nz-j)*nr
          READ(3,650)(gas_velocity_z(ijl),ijl=ij1,ij2)
        END DO
!
        READ(3,640)
        DO j=1,nz
          ij1=1+(nz-j)*nr
          ij2=nr+(nz-j)*nr
          READ(3,650)(gas_velocity_r(ijl),ijl=ij1,ij2)
        END DO
!
          READ(3,640)
        DO j=1,nz
          ij1=1+(nz-j)*nr
          ij2=nr+(nz-j)*nr
          READ(3,650)(gas_temperature(ijl),ijl=ij1,ij2)
        END DO
!
        DO ig=1,ngas
          IF( present_gas(ig) .AND. (ig /= default_gas) ) THEN
            READ(3,640)
            DO j=1,nz
              ij1=1+(nz-j)*nr
              ij2=nr+(nz-j)*nr
              READ(3,650)(gc_molar_fraction(ig,ijl),ijl=ij1,ij2)
            END DO
          END IF
        END DO
!
        DO is=1,nsolid
!
          READ(3,640)
          DO  j=1,nz
            ij1=1+(nz-j)*nr
            ij2=nr+(nz-j)*nr
            READ(3,650)(solid_velocity_z(is,ijl),ijl=ij1,ij2)
          END DO
!
            READ(3,640)
          DO j=1,nz
            ij1=1+(nz-j)*nr 
            ij2=nr+(nz-j)*nr
            READ(3,650)(solid_velocity_r(is,ijl),ijl=ij1,ij2)
          END DO
!
          READ(3,640)
          DO j=1,nz
            ij1=1+(nz-j)*nr
            ij2=nr+(nz-j)*nr
            READ(3,650)(solid_temperature(is,ijl),ijl=ij1,ij2)
          END DO
!
        END DO

        IF (modturbo > 1) THEN
          READ(3,640)
          DO j=1,nz
            ij1=1+(nz-j)*nr
            ij2=nr+(nz-j)*nr
            READ(3,650)(smag_coeff(ijl),ijl=ij1,ij2)
          END DO
        END IF
!
      CLOSE (3)
      END IF

      CALL bcast_real(gas_pressure,SIZE(gas_pressure),root)
      CALL bcast_real(solid_bulk_density,SIZE(solid_bulk_density),root)
      CALL bcast_real(gas_velocity_r,SIZE(gas_velocity_r),root)
      CALL bcast_real(gas_velocity_z,SIZE(gas_velocity_z),root)
      CALL bcast_real(gas_temperature,SIZE(gas_temperature),root)
      CALL bcast_real(gc_molar_fraction,SIZE(gc_molar_fraction),root)
      CALL bcast_real(solid_velocity_r,SIZE(solid_velocity_r),root)
      CALL bcast_real(solid_velocity_z,SIZE(solid_velocity_z),root)
      CALL bcast_real(solid_temperature,SIZE(solid_temperature),root)
!
      RETURN

 647  FORMAT(1x,///,1x,'@@@ TIME = ',G11.4)
 640  FORMAT(///)
 650  FORMAT(1x,10(1x,g12.6))
! 650  FORMAT(1x,10(1x,g14.6e3))
!
     END SUBROUTINE recover

!-----------------------------------------------------------------------------------
      SUBROUTINE outp_test_fluxes(rgfr, rgft)
!
      USE dimensions, ONLY: nr, nz, nsolid
      USE tilde_momentum, ONLY: rug, rwg, rus, rws
      USE tilde_energy, ONLY: rhg, rhk
      USE parallel, ONLY: nproc, mpime, root, group
      USE particles_constants, ONLY: rl, inrl
      USE time_parameters, ONLY: time
!
      IMPLICIT NONE
!
      REAL*8, INTENT(IN) :: rgfr(:), rgft(:)
!
      INTEGER :: ij2, ij1, is, j
      INTEGER :: ijl
      INTEGER :: ig
!
! ... MODIFICARE_X3D
!
      IF( mpime .EQ. root ) THEN

      OPEN(UNIT=3,FILE='OUTPUT_TEST_FLUXES')
!
      WRITE(3,547) time
!
      WRITE(3,100)
      DO j=1,nz
        ij1=1+(nz-j)*nr
        ij2=nr+(nz-j)*nr
        WRITE(3,550)(rug(ijl),ijl=ij1,ij2)
      END DO
!
      WRITE(3,101)
      DO j=1,nz
        ij1=1+(nz-j)*nr
        ij2=nr+(nz-j)*nr
        WRITE(3,550)(rwg(ijl),ijl=ij1,ij2)
      END DO
!
      DO is=1,nsolid
        WRITE(3,102)is 
        DO j=1,nz
          ij1=1+(nz-j)*nr
          ij2=nr+(nz-j)*nr
          WRITE(3,550)(rus(is,ijl)*inrl(is),ijl=ij1,ij2)
        END DO
        WRITE(3,103)is 
        DO j=1,nz
          ij1=1+(nz-j)*nr
          ij2=nr+(nz-j)*nr
          WRITE(3,550)(rws(is,ijl)*inrl(is),ijl=ij1,ij2)
        END DO
      END DO
!
      WRITE(3,106)
      DO j=1,nz
        ij1=1+(nz-j)*nr
        ij2=nr+(nz-j)*nr
        WRITE(3,550)((rgfr(ijl)-rgfr(ijl-1)),ijl=ij1,ij2)
      END DO
!
      WRITE(3,107)
      DO j=1,nz
        ij1=1+(nz-j)*nr
        ij2=nr+(nz-j)*nr
        WRITE(3,550)((rgft(ijl)-rgft(ijl-nr)),ijl=ij1,ij2)
      END DO
!
      CLOSE (3)
      END IF
!
      RETURN

 547  FORMAT(1x,///,1x,'@@@ TIME = ',G11.4)
 550  FORMAT(1x,10(1x,g14.6e3))
 100  FORMAT(1x,//,1x,'RUG',/)
 101  FORMAT(1x,//,1x,'RWG',/)
 102  FORMAT(1x,//,1x,'RUS',I1,/)
 103  FORMAT(1x,//,1x,'RWS',I1,/)
 106  FORMAT(1x,//,1x,'MASFX',/)
 107  FORMAT(1x,//,1x,'MASFZ',/)
!
      END SUBROUTINE outp_test_fluxes
!----------------------------------------------------------------------
      END MODULE output_dump
!----------------------------------------------------------------------
