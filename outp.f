!----------------------------------------------------------------------
      MODULE output_dump
!----------------------------------------------------------------------
      USE kinds
      USE io_restart, ONLY: nfil
      USE io_files, ONLY: logunit

      IMPLICIT NONE
      SAVE

      LOGICAL :: formatted_output

      LOGICAL :: interpolate = .TRUE.
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE jet_out
!
      USE dimensions, ONLY: nsolid
      USE gas_constants, ONLY: gas_type
      USE gas_solid_density, ONLY: rlk, rog
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
      USE io_files, ONLY: filnam, outpunit
!
      IMPLICIT NONE
!
      CHARACTER*4 :: lettera
      LOGICAL :: lform
!
      INTEGER :: ig,is
      REAL*8, ALLOCATABLE :: otmp(:)
!
    
      nfil=nfil+1
      filnam='output.'//lettera(nfil)
      lform = formatted_output

      IF( mpime == root ) THEN

        IF (lform) THEN
          OPEN(UNIT=outpunit,FILE=filnam)
          WRITE(outpunit,*) time
        ELSE 
          OPEN(UNIT=outpunit,FORM='UNFORMATTED',FILE=filnam)
          WRITE(outpunit) REAL(time,4)
        END IF
 
      END IF
!
      CALL write_array( outpunit, p, sgl, lform )  ! gas_pressure
      CALL write_array( outpunit, tg, sgl, lform )   ! gas_temperature
      CALL write_array( outpunit, rog, sgl, lform )  ! gas_density
!
      IF (job_type == '2D') THEN
        CALL write_array( outpunit, ug, sgl, lform ) ! gas_velocity_r
        CALL write_array( outpunit, wg, sgl, lform ) ! gas_velocity_z
      ELSE IF (job_type == '3D') THEN
        CALL write_array( outpunit, ug, sgl, lform ) ! gas_velocity_x
        CALL write_array( outpunit, vg, sgl, lform ) ! gas_velocity_y
        CALL write_array( outpunit, wg, sgl, lform ) ! gas_velocity_z
      ELSE
        CALL error('outp_','Unknown job type',1)
      END IF

      IF( mpime == root ) THEN
        CLOSE (outpunit)
      END IF
!
      RETURN
      END SUBROUTINE jet_out
!----------------------------------------------------------------------
      SUBROUTINE shock_tube_out
!
      USE dimensions, ONLY: nx, nz, nsolid
      USE domain_decomposition, ONLY: flag, ncint
      USE gas_constants, ONLY: gas_type, gammaair
      USE gas_solid_density, ONLY: rog
      USE gas_solid_velocity, ONLY: ug, wg
      USE gas_solid_temperature, ONLY: sieg, tg
      USE parallel, ONLY: nproc, mpime, root, group
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: p
      USE time_parameters, ONLY: time
      USE control_flags, ONLY: job_type
      USE domain_decomposition, ONLY: data_collect, data_distribute
      USE io_files, ONLY: tempunit

!
      IMPLICIT NONE
!
      CHARACTER( LEN = 15 ) :: filnam
      CHARACTER( LEN = 4 ) :: lettera
!
      INTEGER :: i, j, ijk, ig
      REAL*8 :: energy
!
      nfil = nfil + 1
      filnam = 'shtube.' // lettera( nfil )

      IF( mpime == root ) THEN

        OPEN(UNIT=tempunit,FILE=filnam)
        DO ijk = 1, ncint
          IF (flag(ijk) == 1) THEN
            energy = p(ijk)/(rog(ijk)*(gammaair - 1.D0))
            WRITE(tempunit,550)rog(ijk),wg(ijk),p(ijk),energy,tg(ijk)
          END IF
        END DO

        CLOSE (tempunit)

      END IF
!
      RETURN

 550  FORMAT(1x,10(1x,g12.6))
!
      END SUBROUTINE shock_tube_out
!----------------------------------------------------------------------
      SUBROUTINE write_radial_profile_2d
!
      USE dimensions, ONLY: nsolid, ngas, nx, nz
      USE eos_gas, ONLY: xgc, ygc
      USE gas_constants, ONLY: gas_type
      USE gas_solid_density, ONLY: rlk
      USE gas_solid_velocity, ONLY: ug, vg, wg
      USE gas_solid_velocity, ONLY: us, vs, ws
      USE gas_solid_temperature, ONLY: tg, ts
      USE grid, ONLY: zb, r
      USE io_restart, ONLY: write_array
      USE parallel, ONLY: nproc, mpime, root, group
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: p
      USE time_parameters, ONLY: time
      USE turbulence_model, ONLY: modturbo
      USE control_flags, ONLY: job_type
      USE volcano_topography, ONLY: rim_quota
      USE io_files, ONLY: tempunit
      IMPLICIT NONE
      INTEGER :: n, is, ig, kq, k, nq

      DO k = 1, nz
        IF (zb(k) <= rim_quota) kq = k
      END DO
      
      OPEN(tempunit,FILE='radial_profile.dat')
      WRITE(tempunit,*) nx
      WRITE(tempunit,*) (r(n), n=1, nx)
      WRITE(tempunit,*) (ug(n + (kq-1)*nx), n=1, nx)
      WRITE(tempunit,*) (wg(n + (kq-1)*nx), n=1, nx)
      WRITE(tempunit,*) (tg(n + (kq-1)*nx), n=1, nx)
      WRITE(tempunit,*) (p(n + (kq-1)*nx), n=1, nx)
      DO ig = 1, ngas
        WRITE(tempunit,*) (ygc(n,ig), n=1, nx)
      END DO
      DO is = 1, nsolid
        WRITE(tempunit,*) (us(n + (kq-1)*nx,is), n=1, nx)
        WRITE(tempunit,*) (ws(n + (kq-1)*nx,is), n=1, nx)
        WRITE(tempunit,*) (ts(n + (kq-1)*nx,is), n=1, nx)
        WRITE(tempunit,*) (rlk(n + (kq-1)*nx,is)*inrl(is), n=1, nx)
      END DO
      CLOSE(tempunit)

      OPEN(tempunit,FILE='radial_profile_gas.dat')
      DO n = 2, nx-1
        nq = n + (kq-1)*nx
        WRITE(tempunit,188) r(n), ug(nq), wg(nq), tg(nq), p(nq), (ygc(nq,ig), ig=1, ngas)
      END DO
      CLOSE(tempunit)

      OPEN(tempunit,FILE='radial_profile_part.dat')
      DO n = 2, nx-1
        nq = n + (kq-1)*nx
        WRITE(tempunit,199) r(n), (us(nq,is), ws(nq,is), ts(nq,is), rlk(nq,is)*inrl(is), is=1,nsolid)
      END DO
      CLOSE(tempunit)

 188  FORMAT(7(F15.5))       
 199  FORMAT(9(F15.5))       
        
      RETURN
      END SUBROUTINE write_radial_profile_2d
!----------------------------------------------------------------------
      SUBROUTINE outp
!
      USE dimensions, ONLY: nsolid, ngas
      USE eos_gas, ONLY: xgc
      USE gas_constants, ONLY: gas_type
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
      USE io_files, ONLY: filnam, outpunit
!
      IMPLICIT NONE
!
      CHARACTER( LEN =  4 ) :: lettera
      LOGICAL :: lform
!
      INTEGER :: ig,is
      REAL*8, ALLOCATABLE :: otmp(:)
!
    
      filnam='output.'//lettera(nfil)
      lform = formatted_output

      IF( mpime == root ) THEN

        IF (lform) THEN
          OPEN(UNIT=outpunit,FILE=filnam)
          WRITE(outpunit,'(1x,///,1x,"@@@ TIME = ",g11.4)') time
        ELSE 
          OPEN(UNIT=outpunit,FORM='UNFORMATTED',FILE=filnam)
          WRITE(outpunit) REAL(time,4)
        END IF

        WRITE(logunit,fmt="('  from outp: writing file ',A20)") filnam
 
      END IF
!
      IF( lform .AND. mpime == root ) WRITE(outpunit,'(1x,//,1x,"P   ",/)')
      CALL write_array( outpunit, p, sgl, lform )  ! gas_pressure

      IF (job_type == '2D') THEN

        IF( lform .AND. mpime == root ) WRITE(outpunit,'(1x,//,1x,"UG  ",/)')
        CALL write_array( outpunit, ug, sgl, lform ) ! gas_velocity_r

        IF( lform .AND. mpime == root ) WRITE(outpunit,'(1x,//,1x,"WG  ",/)')
        CALL write_array( outpunit, wg, sgl, lform ) ! gas_velocity_z

      ELSE IF (job_type == '3D') THEN

        IF( lform .AND. mpime == root ) WRITE(outpunit,'(1x,//,1x,"UG  ",/)')
        CALL write_array( outpunit, ug, sgl, lform ) ! gas_velocity_x

        IF( lform .AND. mpime == root ) WRITE(outpunit,'(1x,//,1x,"VG  ",/)')
        CALL write_array( outpunit, vg, sgl, lform ) ! gas_velocity_y

        IF( lform .AND. mpime == root ) WRITE(outpunit,'(1x,//,1x,"WG  ",/)')
        CALL write_array( outpunit, wg, sgl, lform ) ! gas_velocity_z

      ELSE
        CALL error('outp_','Unknown job type',1)
      END IF

      IF( lform .AND. mpime == root ) WRITE(outpunit,'(1x,//,1x,"TG  ",/)')
      CALL write_array( outpunit, tg, sgl, lform )  ! gas_temperature
!
      ALLOCATE( otmp( SIZE( xgc, 1 ) ) )
      DO ig=1,ngas
          otmp = xgc(:,ig)
          IF( lform .AND. mpime == root ) WRITE(outpunit,'(1x,//,1x,"XGC",I1,/)') ig
          CALL write_array( outpunit, otmp, sgl, lform )  ! gc_molar_fraction
      END DO
      DEALLOCATE( otmp )
!
      ALLOCATE( otmp( SIZE( rlk, 1 ) ) )

      DO is = 1, nsolid

        otmp = rlk(:,is)*inrl(is)
        IF( lform .AND. mpime == root ) WRITE(outpunit,'(1x,//,1x,"EPS",I1,/)') is
        CALL write_array( outpunit, otmp, sgl, lform )  ! solid_bulk_density

        IF (job_type == '2D') THEN

          IF( lform .AND. mpime == root ) WRITE(outpunit,'(1x,//,1x," US",I1,/)') is
          CALL write_array( outpunit, us(:,is), sgl, lform )  ! solid_velocity_r

          IF( lform .AND. mpime == root ) WRITE(outpunit,'(1x,//,1x," WS",I1,/)') is
          CALL write_array( outpunit, ws(:,is), sgl, lform )  ! solid_velocity_z

        ELSE IF (job_type == '3D') THEN

          IF( lform .AND. mpime == root ) WRITE(outpunit,'(1x,//,1x," US",I1,/)') is
          CALL write_array( outpunit, us(:,is), sgl, lform )  ! solid_velocity_x

          IF( lform .AND. mpime == root ) WRITE(outpunit,'(1x,//,1x," VS",I1,/)') is
          CALL write_array( outpunit, vs(:,is), sgl, lform )  ! solid_velocity_y

          IF( lform .AND. mpime == root ) WRITE(outpunit,'(1x,//,1x," WS",I1,/)') is
          CALL write_array( outpunit, ws(:,is), sgl, lform )  ! solid_velocity_z

        END IF

        IF( lform .AND. mpime == root ) WRITE(outpunit,'(1x,//,1x," TS",I1,/)') is
        CALL write_array( outpunit, ts(:,is), sgl, lform )  ! solid_temperature

      END DO

      DEALLOCATE( otmp )

      IF( mpime == root ) THEN
        CLOSE (outpunit)
      END IF

      nfil=nfil+1
!
      RETURN
      END SUBROUTINE outp
!----------------------------------------------------------------------
      SUBROUTINE outp_recover
!
      USE dimensions, ONLY: nsolid, ngas
      USE eos_gas, ONLY: xgc
      USE gas_constants, ONLY: gas_type
      USE gas_solid_density, ONLY: rlk
      USE gas_solid_velocity, ONLY: ug, vg, wg
      USE gas_solid_velocity, ONLY: us, vs, ws
      USE gas_solid_temperature, ONLY: tg, ts
      USE io_restart, ONLY: read_array
      USE parallel, ONLY: nproc, mpime, root, group
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: p
      USE time_parameters, ONLY: time
      USE turbulence_model, ONLY: modturbo
      USE control_flags, ONLY: job_type
      USE io_files, ONLY: filnam, outpunit
!
      IMPLICIT NONE
!
      CHARACTER( LEN =  4 ) :: lettera
      LOGICAL :: lform
!
      INTEGER :: ig,is
      REAL*8, ALLOCATABLE :: otmp(:)
!
      filnam='output.'//lettera(nfil)
      lform = formatted_output

      IF( mpime == root ) THEN

        IF (lform) THEN
          OPEN(UNIT=outpunit,FILE=filnam)
          READ(outpunit,'(1x,///,1x,"@@@ TIME = ",g11.4)') time
        ELSE 
          OPEN(UNIT=outpunit,FORM='UNFORMATTED',FILE=filnam)
          READ(outpunit) time
        END IF

        WRITE(logunit,fmt="('  from outp: recovering file ',A20)") filnam
        WRITE(logunit,*) 'time = ', time
 
      END IF
!
      CALL bcast_real(time, 1, root)
!
      IF( lform .AND. mpime == root ) READ(outpunit,122)
      p = 0.D0
      CALL read_array( outpunit, p, sgl, lform )  ! gas_pressure

      IF (job_type == '2D') THEN

        IF( lform .AND. mpime == root ) READ(outpunit,122)
        ug = 0.D0
        CALL read_array( outpunit, ug, sgl, lform ) ! gas_velocity_r

        IF( lform .AND. mpime == root ) READ(outpunit,122)
        wg = 0.D0
        CALL read_array( outpunit, wg, sgl, lform ) ! gas_velocity_z

      ELSE IF (job_type == '3D') THEN

        IF( lform .AND. mpime == root ) READ(outpunit,122)
        ug = 0.D0
        CALL read_array( outpunit, ug, sgl, lform ) ! gas_velocity_x

        IF( lform .AND. mpime == root ) READ(outpunit,122)
        vg = 0.D0
        CALL read_array( outpunit, vg, sgl, lform ) ! gas_velocity_y

        IF( lform .AND. mpime == root ) READ(outpunit,122)
        wg = 0.D0
        CALL read_array( outpunit, wg, sgl, lform ) ! gas_velocity_z

      ELSE
        CALL error('outp_','Unknown job type',1)
      END IF

      IF( lform .AND. mpime == root ) READ(outpunit,122)
      tg = 0.D0
      CALL read_array( outpunit, tg, sgl, lform )  ! gas_temperature
!
      ALLOCATE( otmp( SIZE( xgc, 1 ) ) )
      DO ig=1,ngas
          IF( lform .AND. mpime == root ) READ(outpunit,122)
          otmp = 0.D0
          CALL read_array( outpunit, otmp, sgl, lform )  ! gc_molar_fraction
          xgc(:,ig) = otmp
      END DO
      DEALLOCATE( otmp )
!
      ALLOCATE( otmp( SIZE( rlk, 1 ) ) )

      DO is = 1, nsolid

        IF( lform .AND. mpime == root ) READ(outpunit,122)
        otmp = 0.D0
        CALL read_array( outpunit, otmp, sgl, lform )  ! solid_bulk_density
        rlk(:,is) = otmp * rl(is)

        IF (job_type == '2D') THEN

          IF( lform .AND. mpime == root ) READ(outpunit,122)
          us = 0.D0
          CALL read_array( outpunit, us(:,is), sgl, lform )  ! solid_velocity_r

          IF( lform .AND. mpime == root ) READ(outpunit,122)
          ws = 0.D0
          CALL read_array( outpunit, ws(:,is), sgl, lform )  ! solid_velocity_z

        ELSE IF (job_type == '3D') THEN

          IF( lform .AND. mpime == root ) READ(outpunit,122)
          us = 0.D0
          CALL read_array( outpunit, us(:,is), sgl, lform )  ! solid_velocity_x

          IF( lform .AND. mpime == root ) READ(outpunit,122)
          vs = 0.D0
          CALL read_array( outpunit, vs(:,is), sgl, lform )  ! solid_velocity_y

          IF( lform .AND. mpime == root ) READ(outpunit,122)
          ws = 0.D0
          CALL read_array( outpunit, ws(:,is), sgl, lform )  ! solid_velocity_z

        END IF

        IF( lform .AND. mpime == root ) READ(outpunit,122)
        ts = 0.D0
        CALL read_array( outpunit, ts(:,is), sgl, lform )  ! solid_temperature

      END DO

      DEALLOCATE( otmp )

      IF( mpime == root ) THEN
        CLOSE (outpunit)
      END IF

      nfil=nfil+1
!
 122  FORMAT(1x,//,6x,/)

      RETURN
      END SUBROUTINE outp_recover
!-----------------------------------------------------------------------
      END MODULE output_dump
!----------------------------------------------------------------------
