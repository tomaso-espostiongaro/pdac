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
      SUBROUTINE outp_2d
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
      USE turbulence_model, ONLY: smag_coeff, modturbo
      USE control_flags, ONLY: job_type
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
          WRITE(3,550)(solid_bulk_density(ijl,is)*inrl(is),ijl=ij1,ij2)
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
          WRITE(3,550)(solid_velocity_r(ijl,is),ijl=ij1,ij2)
        END DO
!
        WRITE(3,557)is 
        DO  j=1,nz
          ij1=1+(nz-j)*nr
          ij2=nr+(nz-j)*nr
          WRITE(3,550)(solid_velocity_z(ijl,is),ijl=ij1,ij2)
        END DO
!
        WRITE(3,561)is 
        DO j=1,nz
          ij1=1+(nz-j)*nr
          ij2=nr+(nz-j)*nr
          WRITE(3,550)(solid_temperature(ijl,is),ijl=ij1,ij2)
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
      END SUBROUTINE outp_2d
!----------------------------------------------------------------------
     SUBROUTINE recover_2d
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
      USE turbulence_model, ONLY: smag_coeff, modturbo
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
            READ(3,650)(solid_bulk_density(ijl,is),ijl=ij1,ij2)
          END DO
          solid_bulk_density(:,is) = solid_bulk_density(:,is)*rl(is)
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
            READ(3,650)(solid_velocity_z(ijl,is),ijl=ij1,ij2)
          END DO
!
            READ(3,640)
          DO j=1,nz
            ij1=1+(nz-j)*nr 
            ij2=nr+(nz-j)*nr
            READ(3,650)(solid_velocity_r(ijl,is),ijl=ij1,ij2)
          END DO
!
          READ(3,640)
          DO j=1,nz
            ij1=1+(nz-j)*nr
            ij2=nr+(nz-j)*nr
            READ(3,650)(solid_temperature(ijl,is),ijl=ij1,ij2)
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
     END SUBROUTINE recover_2d

!-----------------------------------------------------------------------------------
      SUBROUTINE outp_test_fluxes(rgfr, rgft)
!
      USE dimensions, ONLY: nr, nz, nsolid
      USE tilde_momentum, ONLY: rug, rwg, rus, rws
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
          WRITE(3,550)(rus(ijl,is)*inrl(is),ijl=ij1,ij2)
        END DO
        WRITE(3,103)is 
        DO j=1,nz
          ij1=1+(nz-j)*nr
          ij2=nr+(nz-j)*nr
          WRITE(3,550)(rws(ijl,is)*inrl(is),ijl=ij1,ij2)
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
      SUBROUTINE outp_3d
!
      USE dimensions, ONLY: nx, ny, nz, nsolid
      USE eos_gas, ONLY: gc_molar_fraction
      USE gas_constants, ONLY: present_gas, default_gas
      USE gas_solid_density, ONLY: solid_bulk_density
      USE gas_solid_velocity, ONLY: gas_velocity_x, gas_velocity_y, &
      &                             gas_velocity_z
      USE gas_solid_velocity, ONLY: solid_velocity_x, solid_velocity_y, &
      &                             solid_velocity_z
      USE gas_solid_temperature, ONLY: gas_temperature, solid_temperature
      USE parallel, ONLY: nproc, mpime, root, group
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: gas_pressure
      USE time_parameters, ONLY: time
      USE turbulence_model, ONLY: smag_coeff, modturbo
      USE control_flags, ONLY: job_type
!
      IMPLICIT NONE
!
      CHARACTER :: filnam*11
      CHARACTER*4 :: letter
!
      INTEGER :: k, ijk
      INTEGER :: ijk2, ijk1, is
      INTEGER :: ig
!
      nfil=nfil+1
      filnam='OUTPUT.'//letter(nfil)

      IF( mpime == root ) THEN

      OPEN(UNIT=3,FILE=filnam)
!
      WRITE(3,547) time
!
      WRITE(3,548)
      DO k=1,nz
        ijk1 = (k-1)*nx*ny +1
        ijk2 = k*nx*ny
        WRITE(3,550)(gas_pressure(ijk),ijk = ijk1, ijk2 )
      END DO
!
      DO is=1,nsolid
        WRITE(3,549)is 
        DO k=1,nz
          ijk1 = (k-1)*nx*ny +1
          ijk2 = k*nx*ny
          WRITE(3,550)(solid_bulk_density(ijk,is)*inrl(is),ijk=ijk1,ijk2)
        END DO
      END DO
!
      WRITE(3,552)
      DO k=1,nz
          ijk1 = (k-1)*nx*ny +1
          ijk2 = k*nx*ny
        WRITE(3,550)(gas_velocity_x(ijk),ijk=ijk1,ijk2)
      END DO
!
      WRITE(3,552)
      DO k=1,nz
          ijk1 = (k-1)*nx*ny +1
          ijk2 = k*nx*ny
        WRITE(3,550)(gas_velocity_y(ijk),ijk=ijk1,ijk2)
      END DO
!
      WRITE(3,553)
      DO k=1,nz
          ijk1 = (k-1)*nx*ny +1
          ijk2 = k*nx*ny
        WRITE(3,550)(gas_velocity_z(ijk),ijk=ijk1,ijk2)
      END DO
!
      WRITE(3,560)
      DO k=1,nz
          ijk1 = (k-1)*nx*ny +1
          ijk2 = k*nx*ny
        WRITE(3,550)(gas_temperature(ijk),ijk=ijk1,ijk2)
      END DO
!
      DO ig=1,ngas
        IF( present_gas(ig) .AND. (ig /= default_gas) ) THEN
          WRITE(3,562) ig
          DO k=1,nz
          ijk1 = (k-1)*nx*ny +1
          ijk2 = k*nx*ny
            WRITE(3,550)(gc_molar_fraction(ig,ijk),ijk=ijk1,ijk2)
          END DO
        END IF
      END DO
!
      DO is=1,nsolid
!
        WRITE(3,556)is 
        DO k=1,nz
          ijk1 = (k-1)*nx*ny +1
          ijk2 = k*nx*ny
          WRITE(3,550)(solid_velocity_x(ijk,is),ijk=ijk1,ijk2)
        END DO
!
        WRITE(3,556)is 
        DO k=1,nz
          ijk1 = (k-1)*nx*ny +1
          ijk2 = k*nx*ny
          WRITE(3,550)(solid_velocity_y(ijk,is),ijk=ijk1,ijk2)
        END DO
!
        WRITE(3,557)is 
        DO  k=1,nz
          ijk1 = (k-1)*nx*ny +1
          ijk2 = k*nx*ny
          WRITE(3,550)(solid_velocity_z(ijk,is),ijk=ijk1,ijk2)
        END DO
!
        WRITE(3,561)is 
        DO k=1,nz
          ijk1 = (k-1)*nx*ny +1
          ijk2 = k*nx*ny
          WRITE(3,550)(solid_temperature(ijk,is),ijk=ijk1,ijk2)
        END DO
!
      END DO

      IF (modturbo > 1) THEN
        WRITE(3,555)
        DO k=1,nz
          ijk1 = (k-1)*nx*ny +1
          ijk2 = k*nx*ny
          WRITE(3,550)(smag_coeff(ijk),ijk=ijk1,ijk2)
        END DO
      END IF
!
      CLOSE (3)
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
      END SUBROUTINE outp_3d
!----------------------------------------------------------------------
      SUBROUTINE outp_bin
!
      USE dimensions, ONLY: nsolid
      USE eos_gas, ONLY: gc_molar_fraction
      USE gas_constants, ONLY: present_gas, default_gas
      USE gas_solid_density, ONLY: solid_bulk_density
      USE gas_solid_velocity, ONLY: gas_velocity_r, gas_velocity_x, &
                                    gas_velocity_y, gas_velocity_z
      USE gas_solid_velocity, ONLY: solid_velocity_r, solid_velocity_x, &
                                    solid_velocity_y, solid_velocity_z
      USE gas_solid_temperature, ONLY: gas_temperature, solid_temperature
      USE parallel, ONLY: nproc, mpime, root, group
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: gas_pressure
      USE time_parameters, ONLY: time
      USE turbulence_model, ONLY: smag_coeff, modturbo
      USE control_flags, ONLY: job_type
!
      IMPLICIT NONE
!
      CHARACTER :: filnam*11
      CHARACTER*4 :: letter
!
      INTEGER :: ig,is
!
      nfil=nfil+1
      filnam='output.'//letter(nfil)

      IF( mpime == root ) THEN

        OPEN(UNIT=12,FORM='UNFORMATTED',FILE=filnam)
!
        WRITE(12) REAL(time,4)
!
        WRITE(12) REAL(gas_pressure,4)
        IF (job_type == '2D') THEN
          WRITE(12) REAL(gas_velocity_r,4)
          WRITE(12) REAL(gas_velocity_z,4)
        ELSE IF (job_type == '3D') THEN
          WRITE(12) REAL(gas_velocity_x,4)
          WRITE(12) REAL(gas_velocity_y,4)
          WRITE(12) REAL(gas_velocity_z,4)
        ELSE
          CALL error('outp_bin','Unknown job type',1)
        END IF
        WRITE(12) REAL(gas_temperature,4)
!
        DO ig=1,ngas
          IF( present_gas(ig) .AND. (ig /= default_gas) ) THEN
            WRITE(12) REAL(gc_molar_fraction(ig,:),4)
          END IF
        END DO
!
        DO is=1,nsolid
          WRITE(12) REAL(solid_bulk_density(:,is)*inrl(is),4)
          IF (job_type == '2D') THEN
            WRITE(12) REAL(solid_velocity_r(:,is),4)
            WRITE(12) REAL(solid_velocity_z(:,is),4)
          ELSE IF (job_type == '3D') THEN
            WRITE(12) REAL(solid_velocity_x(:,is),4)
            WRITE(12) REAL(solid_velocity_y(:,is),4)
            WRITE(12) REAL(solid_velocity_z(:,is),4)
          END IF
          WRITE(12) REAL(solid_temperature(:,is),4)
        END DO

        CLOSE (12)
      END IF
!
      RETURN
      END SUBROUTINE outp_bin
!-----------------------------------------------------------------------
      END MODULE output_dump
!----------------------------------------------------------------------
