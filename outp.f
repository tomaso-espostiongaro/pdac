!----------------------------------------------------------------------
      MODULE output_dump
!----------------------------------------------------------------------
      USE kinds
      USE dimensions
      IMPLICIT NONE
      INTEGER nfil
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE outp_2d
!
      USE dimensions, ONLY: nx, nz, nsolid
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
      USE domain_decomposition, ONLY: data_collect, data_distribute
      USE dimensions

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
        ij1=1+(nz-j)*nx
        ij2=nx+(nz-j)*nx
        WRITE(3,550)(gas_pressure(ijl),ijl=ij1,ij2)
      END DO
!
      DO is=1,nsolid
        WRITE(3,549)is 
        DO j=1,nz
          ij1=1+(nz-j)*nx
          ij2=nx+(nz-j)*nx
          WRITE(3,550)(solid_bulk_density(ijl,is)*inrl(is),ijl=ij1,ij2)
        END DO
      END DO
!
      WRITE(3,552)
      DO j=1,nz
        ij1=1+(nz-j)*nx
        ij2=nx+(nz-j)*nx
        WRITE(3,550)(gas_velocity_r(ijl),ijl=ij1,ij2)
      END DO
!
      WRITE(3,553)
      DO j=1,nz
        ij1=1+(nz-j)*nx
        ij2=nx+(nz-j)*nx
        WRITE(3,550)(gas_velocity_z(ijl),ijl=ij1,ij2)
      END DO
!
      WRITE(3,560)
      DO j=1,nz
        ij1=1+(nz-j)*nx
        ij2=nx+(nz-j)*nx
        WRITE(3,550)(gas_temperature(ijl),ijl=ij1,ij2)
      END DO
!
      DO ig=1,ngas
        IF( present_gas(ig) .AND. (ig /= default_gas) ) THEN
          WRITE(3,562) ig
          DO j=1,nz
            ij1=1+(nz-j)*nx
            ij2=nx+(nz-j)*nx
            WRITE(3,550)(gc_molar_fraction(ig,ijl),ijl=ij1,ij2)
          END DO
        END IF
      END DO
!
      DO is=1,nsolid
!
        WRITE(3,556)is 
        DO j=1,nz
          ij1=1+(nz-j)*nx 
          ij2=nx+(nz-j)*nx
          WRITE(3,550)(solid_velocity_r(ijl,is),ijl=ij1,ij2)
        END DO
!
        WRITE(3,557)is 
        DO  j=1,nz
          ij1=1+(nz-j)*nx
          ij2=nx+(nz-j)*nx
          WRITE(3,550)(solid_velocity_z(ijl,is),ijl=ij1,ij2)
        END DO
!
        WRITE(3,561)is 
        DO j=1,nz
          ij1=1+(nz-j)*nx
          ij2=nx+(nz-j)*nx
          WRITE(3,550)(solid_temperature(ijl,is),ijl=ij1,ij2)
        END DO
!
      END DO

      IF (modturbo > 1) THEN
        WRITE(3,555)
        DO j=1,nz
          ij1=1+(nz-j)*nx
          ij2=nx+(nz-j)*nx
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
      USE dimensions, ONLY: nx, nz, nsolid
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
          ij1=1+(nz-j)*nx
          ij2=nx+(nz-j)*nx
          READ(3,650)(gas_pressure(ijl),ijl=ij1,ij2)
        END DO
!
        DO is=1,nsolid
          READ(3,640)
          DO j=1,nz
            ij1=1+(nz-j)*nx
            ij2=nx+(nz-j)*nx
            READ(3,650)(solid_bulk_density(ijl,is),ijl=ij1,ij2)
          END DO
          solid_bulk_density(:,is) = solid_bulk_density(:,is)*rl(is)
        END DO
!
        READ(3,640)
        DO j=1,nz
          ij1=1+(nz-j)*nx
          ij2=nx+(nz-j)*nx
          READ(3,650)(gas_velocity_z(ijl),ijl=ij1,ij2)
        END DO
!
        READ(3,640)
        DO j=1,nz
          ij1=1+(nz-j)*nx
          ij2=nx+(nz-j)*nx
          READ(3,650)(gas_velocity_r(ijl),ijl=ij1,ij2)
        END DO
!
          READ(3,640)
        DO j=1,nz
          ij1=1+(nz-j)*nx
          ij2=nx+(nz-j)*nx
          READ(3,650)(gas_temperature(ijl),ijl=ij1,ij2)
        END DO
!
        DO ig=1,ngas
          IF( present_gas(ig) .AND. (ig /= default_gas) ) THEN
            READ(3,640)
            DO j=1,nz
              ij1=1+(nz-j)*nx
              ij2=nx+(nz-j)*nx
              READ(3,650)(gc_molar_fraction(ig,ijl),ijl=ij1,ij2)
            END DO
          END IF
        END DO
!
        DO is=1,nsolid
!
          READ(3,640)
          DO  j=1,nz
            ij1=1+(nz-j)*nx
            ij2=nx+(nz-j)*nx
            READ(3,650)(solid_velocity_z(ijl,is),ijl=ij1,ij2)
          END DO
!
            READ(3,640)
          DO j=1,nz
            ij1=1+(nz-j)*nx 
            ij2=nx+(nz-j)*nx
            READ(3,650)(solid_velocity_r(ijl,is),ijl=ij1,ij2)
          END DO
!
          READ(3,640)
          DO j=1,nz
            ij1=1+(nz-j)*nx
            ij2=nx+(nz-j)*nx
            READ(3,650)(solid_temperature(ijl,is),ijl=ij1,ij2)
          END DO
!
        END DO

        IF (modturbo > 1) THEN
          READ(3,640)
          DO j=1,nz
            ij1=1+(nz-j)*nx
            ij2=nx+(nz-j)*nx
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
      WRITE(3,553)
      DO k=1,nz
          ijk1 = (k-1)*nx*ny +1
          ijk2 = k*nx*ny
        WRITE(3,550)(gas_velocity_y(ijk),ijk=ijk1,ijk2)
      END DO
!
      WRITE(3,554)
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
        WRITE(3,557)is 
        DO k=1,nz
          ijk1 = (k-1)*nx*ny +1
          ijk2 = k*nx*ny
          WRITE(3,550)(solid_velocity_y(ijk,is),ijk=ijk1,ijk2)
        END DO
!
        WRITE(3,558)is 
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
 554  FORMAT(1x,//,1x,'WG',/)
 556  FORMAT(1x,//,1x,'UP',I1,/)
 557  FORMAT(1x,//,1x,'VP',I1,/)
 558  FORMAT(1x,//,1x,'WP',I1,/)
 560  FORMAT(1x,//,1x,'TG',/)
 561  FORMAT(1x,//,1x,'TP',I1,/)
 562  FORMAT(1x,//,1x,'XGC',I1,/)

 555  FORMAT(1x,//,1x,'CDYN',/)
!
      END SUBROUTINE outp_3d

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
      USE turbulence_model, ONLY: smag_coeff, modturbo
      USE control_flags, ONLY: job_type
!
      IMPLICIT NONE
!
      CHARACTER :: filnam*11
      CHARACTER*4 :: letter
!
      INTEGER :: ig,is
      LOGICAL :: lform = .TRUE.
      REAL*8, ALLOCATABLE :: otmp(:)
!
      nfil=nfil+1
      filnam='output.'//letter(nfil)

      IF( mpime == root ) THEN

        IF (lform) THEN
          OPEN(UNIT=12,FILE=filnam)
          WRITE(12,*) time
        ELSE 
          OPEN(UNIT=12,FORM='UNFORMATTED',FILE=filnam)
          WRITE(12) REAL(time,4)
        END IF
 
      END IF
!
      CALL write_array( 12, p, sgl, lform )  ! gas_pressure

      IF (job_type == '2D') THEN
        CALL write_array( 12, ug, sgl, lform ) !  REAL(gas_velocity_r,4)
        CALL write_array( 12, wg, sgl, lform ) !  REAL(gas_velocity_z,4)
      ELSE IF (job_type == '3D') THEN
        CALL write_array( 12, ug, sgl, lform ) !  REAL(gas_velocity_x,4)
        CALL write_array( 12, vg, sgl, lform ) !  REAL(gas_velocity_y,4)
        CALL write_array( 12, wg, sgl, lform ) !  REAL(gas_velocity_z,4)
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
      DO is=1,nsolid
        otmp = rlk(:,is)*inrl(is)
        CALL write_array( 12, otmp, sgl, lform )  ! solid_bulk_density
      END DO
      DEALLOCATE( otmp )

      DO is = 1, nsolid
        IF (job_type == '2D') THEN
          CALL write_array( 12, us(:,is), sgl, lform )  ! solid_velocity_r
          CALL write_array( 12, ws(:,is), sgl, lform )  ! solid_velocity_z
        ELSE IF (job_type == '3D') THEN
          CALL write_array( 12, us(:,is), sgl, lform )  ! solid_velocity_x
          CALL write_array( 12, vs(:,is), sgl, lform )  ! solid_velocity_y
          CALL write_array( 12, ws(:,is), sgl, lform )  ! solid_velocity_z
        END IF
      END DO

      DO is=1,nsolid
        CALL write_array( 12, ts(:,is), sgl, lform )  ! solid_temperature
      END DO

      IF( mpime == root ) THEN
        CLOSE (12)
      END IF
!
      RETURN
      END SUBROUTINE outp
!-----------------------------------------------------------------------
      END MODULE output_dump
!----------------------------------------------------------------------
