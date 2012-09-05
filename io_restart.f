!----------------------------------------------------------------------
      MODULE io_restart
!----------------------------------------------------------------------
      USE control_flags, ONLY: job_type, nfil, formatted_input
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions
      USE domain_mapping, ONLY: data_collect, data_distribute, ncint
      USE eos_gas, ONLY: xgc, ygc
      USE gas_solid_density, ONLY: rog, rgp, rlk
      USE gas_solid_velocity, ONLY: ug, vg, wg
      USE gas_solid_velocity, ONLY: us, vs, ws
      USE gas_solid_temperature, ONLY: sieg, ts, sies, tg 
      USE gas_constants, ONLY: gas_type
      USE kinds
      USE parallel, ONLY: mpime, root
      USE pressure_epsilon, ONLY: ep, p
      USE time_parameters, ONLY: time
      USE io_files, ONLY: errorunit, logunit
      USE io_parallel, ONLY: read_array, write_array
!
      IMPLICIT NONE
      PRIVATE

      REAL*8 :: max_seconds

      PUBLIC :: taperd, tapewr, outp_recover, outp_remap
      PUBLIC :: max_seconds
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE tapewr

      USE io_files, ONLY: resfile, resunit
!
      IMPLICIT NONE
!
      INTEGER :: i, j, k, ijk, is
      INTEGER :: ig
      LOGICAL :: lform = .FALSE.
!
      IF( mpime == root ) THEN
!
        WRITE(logunit,*) ' writing to restart file '
!
        OPEN(UNIT=resunit,form='unformatted', FILE = resfile)
!
        WRITE(resunit) time, nx, ny, nz, nsolid, ngas, nfil
        WRITE(resunit) gas_type

      END IF
!
! ... store the final values of the main physical variables 
! 

      CALL write_array( resunit, p, dbl, lform )
      
      IF (job_type == JOB_TYPE_2D) THEN

        CALL write_array( resunit, ug, dbl, lform )
        CALL write_array( resunit, wg, dbl, lform )

      ELSE IF (job_type == JOB_TYPE_3D) THEN

        CALL write_array( resunit, ug, dbl, lform )
        CALL write_array( resunit, vg, dbl, lform )
        CALL write_array( resunit, wg, dbl, lform )

      ELSE

          CALL error('tapewr','Unknown job type',1)

      END IF

      CALL write_array( resunit, sieg, dbl, lform )

      DO is = 1, nsolid
        CALL write_array( resunit, rlk(:,is), dbl, lform )
      END DO

      IF (job_type == JOB_TYPE_2D) THEN

        DO is = 1, nsolid
          CALL write_array( resunit, us(:,is), dbl, lform )
          CALL write_array( resunit, ws(:,is), dbl, lform )
        END DO

      ELSE IF (job_type == JOB_TYPE_3D) THEN

        DO is = 1, nsolid
          CALL write_array( resunit, us(:,is), dbl, lform )
          CALL write_array( resunit, vs(:,is), dbl, lform )
          CALL write_array( resunit, ws(:,is), dbl, lform )
        END DO

      ELSE

          CALL error('tapewr','Unknown job type',1)

      END IF

      DO is = 1, nsolid
        CALL write_array( resunit, sies(:,is), dbl, lform )
      END DO

      DO ig = 1, ngas
        CALL write_array( resunit, ygc(:,ig), dbl, lform )
      END DO
!
      CALL write_array( resunit, rgp, dbl, lform )

      CALL write_array( resunit, rog, dbl, lform )

      CALL write_array( resunit, tg, dbl, lform )
      
      DO is = 1, nsolid
        CALL write_array( resunit, ts(:,is), dbl, lform )
      END DO
!
      IF( mpime .EQ. root ) THEN
        CLOSE(resunit)
      END IF

      RETURN
      END SUBROUTINE tapewr
!----------------------------------------------------------------------
      SUBROUTINE taperd
!
      USE dimensions
      USE domain_mapping, ONLY: ncint, meshinds
      USE atmospheric_conditions, ONLY: p_atm, t_atm
      USE grid, ONLY: zb, dz
      USE specific_heat_module, ONLY: hcapg
      USE gas_constants, ONLY: tzero, hzerog
      USE io_files, ONLY: resfile, resunit 


      IMPLICIT NONE
!
      INTEGER :: i, j, k, ijk, is, imesh
      INTEGER :: nx_, ny_, nz_, nsolid_, ngas_
      INTEGER :: ig, info
      LOGICAL :: lform = .FALSE.
      REAL*8 :: zrif, trif, prif, hc
      REAL*8, ALLOCATABLE :: rtmp( :, : )
      

      IF( mpime == root ) THEN

        WRITE(logunit,*) ' reading from restart file '
!
        OPEN(UNIT=resunit,form='unformatted',FILE=resfile)

        READ(resunit,ERR=199) time, nx_, ny_, nz_, nsolid_, ngas_, nfil
        READ(resunit,ERR=199) gas_type

        WRITE(logunit,*) ' time =  ', time
        WRITE(logunit,*) ' nx   =  ', nx
        WRITE(logunit,*) ' ny   =  ', ny
        WRITE(logunit,*) ' nz   =  ', nz
        WRITE(logunit,*) ' nsolid   =  ', nsolid
        WRITE(logunit,*) ' ngas     =  ', ngas
        WRITE(logunit,*) ' nfil     =  ', nfil
        WRITE(logunit,*) ' gas_type =  ', gas_type

      END IF

      CALL bcast_real(time, 1, root)
      CALL bcast_integer(nx_, 1, root)
      CALL bcast_integer(ny_, 1, root)
      CALL bcast_integer(nz_, 1, root)
      CALL bcast_integer(nsolid_, 1, root)
      CALL bcast_integer(ngas_, 1, root)
      CALL bcast_integer(nfil, 1, root)
      CALL bcast_integer(gas_type, ngas, root)

      IF( nx_ /= nx ) &
        CALL error(' taperd ',' inconsistent dimension nx ', nx_ )
      IF( ny_ /= ny ) &
        CALL error(' taperd ',' inconsistent dimension ny ', ny_ )
      IF( nz_ /= nz ) &
        CALL error(' taperd ',' inconsistent dimension nz ', nz_ )
      IF( nsolid_ /= nsolid ) &
        CALL error(' taperd ',' inconsistent dimension nsolid ', nsolid_ )
      IF( ngas_ /= ngas ) &
        CALL error(' taperd ',' inconsistent dimension ngas ', ngas_ )

      !
      !  read pressure

      p = 0.0d0
      CALL read_array( resunit, p, dbl, lform )
      IF( ANY( p < 0 ) ) THEN
         IF( mpime == root ) WRITE(errorunit,*) 'WARNING reading restart, p < 0'
      END IF

      !
      !  read gas velocity

      IF (job_type == JOB_TYPE_2D) THEN

        ug = 0.0d0
        wg = 0.0d0
        CALL read_array( resunit, ug, dbl, lform )
        CALL read_array( resunit, wg, dbl, lform )

      ELSE IF (job_type == JOB_TYPE_3D) THEN

        ug = 0.0d0
        vg = 0.0d0
        wg = 0.0d0
        CALL read_array( resunit, ug, dbl, lform )
        CALL read_array( resunit, vg, dbl, lform )
        CALL read_array( resunit, wg, dbl, lform )

      ELSE

          CALL error('taperd','Unknown job type',1)

      END IF

      !
      !  read gas enthalpy

      sieg = 0.0d0
      CALL read_array( resunit, sieg, dbl, lform )
      IF( ANY( sieg < 0 ) ) THEN
         IF( mpime == root ) WRITE(errorunit,*) 'WARNING reading restart, sieg < 0'
      END IF

      !
      !  read particle density

      rlk = 0.0d0
      DO is = 1, nsolid
        CALL read_array( resunit, rlk(:,is), dbl, lform )
      END DO
      IF( ANY( rlk < 0 ) ) THEN
         IF( mpime == root ) WRITE(errorunit,*) 'WARNING reading restart, rlk < 0'
      END IF

      !
      !  read particles velocity

      IF (job_type == JOB_TYPE_2D) THEN

        us = 0.0d0
        ws = 0.0d0
        DO is = 1, nsolid
          CALL read_array( resunit, us(:,is), dbl, lform )
          CALL read_array( resunit, ws(:,is), dbl, lform )
        END DO

      ELSE IF (job_type == JOB_TYPE_3D) THEN

        us = 0.0d0
        vs = 0.0d0
        ws = 0.0d0
        DO is = 1, nsolid
          CALL read_array( resunit, us(:,is), dbl, lform )
          CALL read_array( resunit, vs(:,is), dbl, lform )
          CALL read_array( resunit, ws(:,is), dbl, lform )
        END DO

      ELSE

          CALL error('tapewr','Unknown job type',1)

      END IF

      !
      !  read particle enthalpy

      sies = 0.0d0
      DO is = 1, nsolid
        CALL read_array( resunit, sies(:,is), dbl, lform )
      END DO
      IF( ANY( sies < 0.0d0 ) ) THEN
         IF( mpime == root ) WRITE(errorunit,*) 'WARNING reading restart, sies < 0'
      END IF

      !
      !  read gas components

      ygc = 0.0d0
      DO ig = 1, ngas
        CALL read_array( resunit, ygc(:,ig), dbl, lform )
      END DO
      IF( ANY( ygc < 0 ) ) THEN
         IF( mpime == root ) WRITE(errorunit,*) 'WARNING reading restart, ygc < 0'
      END IF

      !
      ! read gas density

      rgp = 0.0d0
      CALL read_array( resunit, rgp, dbl, lform )
      IF( ANY( rgp < 0 ) ) THEN
         IF( mpime == root ) WRITE(errorunit,*) 'WARNING reading restart, rgp < 0'
      END IF

      rog = 0.0d0
      CALL read_array( resunit, rog, dbl, lform )
      IF( ANY( rog < 0 ) ) THEN
         IF( mpime == root ) WRITE(errorunit,*) 'WARNING reading restart, rog < 0'
      END IF

      !
      ! read gas and particle temperature

      tg = 0.0d0
      CALL read_array( resunit, tg, dbl, lform ) 
      IF( ANY( tg < 0 ) ) THEN
         IF( mpime == root ) WRITE(errorunit,*) 'WARNING reading restart, tg < 0'
      END IF

      ts = 0.0d0
      DO is = 1, nsolid
        CALL read_array( resunit, ts(:,is), dbl, lform )
      END DO
      IF( ANY( ts < 0 ) ) THEN
         IF( mpime == root ) &
           WRITE(errorunit,*) 'WARNING reading restart, ts < 0'
      END IF

      IF( mpime == root ) THEN
        CLOSE (resunit)
        WRITE(logunit,*) ' restart file has been red '
      END IF
!
      RETURN
!
 199  CALL error('io_restart.f', 'error in reading resunit', resunit)
!
      END SUBROUTINE taperd
!----------------------------------------------------------------------
      SUBROUTINE outp_recover(nf)
!
      USE dimensions, ONLY: nsolid, ngas
      USE eos_gas, ONLY: xgc
      USE gas_constants, ONLY: gas_type
      USE gas_solid_density, ONLY: rlk
      USE gas_solid_velocity, ONLY: ug, vg, wg
      USE gas_solid_velocity, ONLY: us, vs, ws
      USE gas_solid_temperature, ONLY: tg, ts
      USE io_parallel, ONLY: read_array
      USE parallel, ONLY: nproc, mpime, root, group
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: p
      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE io_files, ONLY: filnam, outpunit
!
      IMPLICIT NONE
!
      INTEGER, INTENT(INOUT) :: nf
      CHARACTER( LEN =  4 ) :: lettera
      LOGICAL :: lform
      INTEGER :: ig,is
      REAL*4 :: time4
      REAL*8, ALLOCATABLE :: otmp(:)
!
      filnam='output.'//lettera(nf)
      lform = formatted_input

      IF( mpime == root ) THEN

        IF (lform) THEN
          OPEN(UNIT=outpunit,FILE=filnam)
          !READ(outpunit,222) time
          READ(outpunit,*,ERR=199)
          READ(outpunit,*,ERR=199)
          READ(outpunit,*,ERR=199)
          READ(outpunit,'(12x,g11.4)',ERR=199) time
        ELSE 
          OPEN(UNIT=outpunit,FORM='UNFORMATTED',FILE=filnam)
          READ(outpunit,ERR=199) time4
          time = REAL(time4,dbl)
        END IF
!
        WRITE(logunit,fmt="('  from outp: recovering file ',A20)") filnam
        WRITE(logunit,*) 'time = ', time
 
      END IF
 222  FORMAT(1x,///,1x,"@@@ TIME = ",g11.4)
!
      CALL bcast_real(time, 1, root)
!
      IF( lform .AND. mpime == root ) READ(outpunit,122,ERR=199)
      p = 0.D0
      CALL read_array( outpunit, p, sgl, lform )  ! gas_pressure

      IF (job_type == JOB_TYPE_2D) THEN

        IF( lform .AND. mpime == root ) READ(outpunit,122,ERR=199)
        ug = 0.D0
        CALL read_array( outpunit, ug, sgl, lform ) ! gas_velocity_r

        IF( lform .AND. mpime == root ) READ(outpunit,122,ERR=199)
        wg = 0.D0
        CALL read_array( outpunit, wg, sgl, lform ) ! gas_velocity_z

      ELSE IF (job_type == JOB_TYPE_3D) THEN

        IF( lform .AND. mpime == root ) READ(outpunit,122,ERR=199)
        ug = 0.D0
        CALL read_array( outpunit, ug, sgl, lform ) ! gas_velocity_x

        IF( lform .AND. mpime == root ) READ(outpunit,122,ERR=199)
        vg = 0.D0
        CALL read_array( outpunit, vg, sgl, lform ) ! gas_velocity_y

        IF( lform .AND. mpime == root ) READ(outpunit,122,ERR=199)
        wg = 0.D0
        CALL read_array( outpunit, wg, sgl, lform ) ! gas_velocity_z

      ELSE
        CALL error('outp_','Unknown job type',1)
      END IF

      IF( lform .AND. mpime == root ) READ(outpunit,122,ERR=199)
      tg = 0.D0
      CALL read_array( outpunit, tg, sgl, lform )  ! gas_temperature
!
      ALLOCATE( otmp( SIZE( xgc, 1 ) ) )
      DO ig=1,ngas
          IF( lform .AND. mpime == root ) READ(outpunit,122,ERR=199)
          otmp = 0.D0
          CALL read_array( outpunit, otmp, sgl, lform )  ! gc_molar_fraction
          xgc(:,ig) = otmp
      END DO
      DEALLOCATE( otmp )
!
      ALLOCATE( otmp( SIZE( rlk, 1 ) ) )

      DO is = 1, nsolid

        IF( lform .AND. mpime == root ) READ(outpunit,122,ERR=199)
        otmp = 0.D0
        CALL read_array( outpunit, otmp, sgl, lform )  ! solid_bulk_density
        rlk(:,is) = otmp * rl(is)

        IF (job_type == JOB_TYPE_2D) THEN

          IF( lform .AND. mpime == root ) READ(outpunit,122,ERR=199)
          us(:,is) = 0.D0
          CALL read_array( outpunit, us(:,is), sgl, lform )  ! solid_velocity_r

          IF( lform .AND. mpime == root ) READ(outpunit,122,ERR=199)
          ws(:,is) = 0.D0
          CALL read_array( outpunit, ws(:,is), sgl, lform )  ! solid_velocity_z

        ELSE IF (job_type == JOB_TYPE_3D) THEN

          IF( lform .AND. mpime == root ) READ(outpunit,122,ERR=199)
          us(:,is) = 0.D0
          CALL read_array( outpunit, us(:,is), sgl, lform )  ! solid_velocity_x

          IF( lform .AND. mpime == root ) READ(outpunit,122,ERR=199)
          vs(:,is) = 0.D0
          CALL read_array( outpunit, vs(:,is), sgl, lform )  ! solid_velocity_y

          IF( lform .AND. mpime == root ) READ(outpunit,122,ERR=199)
          ws(:,is) = 0.D0
          CALL read_array( outpunit, ws(:,is), sgl, lform )  ! solid_velocity_z

        END IF

        IF( lform .AND. mpime == root ) READ(outpunit,122,ERR=199)
        ts(:,is) = 0.D0
        CALL read_array( outpunit, ts(:,is), sgl, lform )  ! solid_temperature

      END DO

      DEALLOCATE( otmp )

      IF( mpime == root ) THEN
        CLOSE (outpunit)
      END IF
!
      nf = nf + 1
!
 122  FORMAT(1x,//,6x,/)

      RETURN
!
 199    CALL error('io_restart.f', 'error reading outputunit',outpunit)
!
      END SUBROUTINE outp_recover
!----------------------------------------------------------------------
      SUBROUTINE outp_remap(nf)
!
      USE dimensions, ONLY: nsolid, ngas
      USE eos_gas, ONLY: xgc
      USE gas_constants, ONLY: gas_type
      USE gas_solid_density, ONLY: rlk
      USE gas_solid_velocity, ONLY: ug, vg, wg
      USE gas_solid_velocity, ONLY: us, vs, ws
      USE gas_solid_temperature, ONLY: tg, ts
      USE io_parallel, ONLY: read_inner_array
      USE parallel, ONLY: nproc, mpime, root, group
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: p
      USE time_parameters, ONLY: time
      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE io_files, ONLY: filnam, outpunit
!
      IMPLICIT NONE
!
      INTEGER, INTENT(INOUT) :: nf
      CHARACTER( LEN =  4 ) :: lettera
      LOGICAL :: lform
      INTEGER :: ig,is
      REAL*4 :: time4
      REAL*8, ALLOCATABLE :: otmp(:)
!
      filnam='output.'//lettera(nf)
      lform = formatted_input

      IF( mpime == root ) THEN

        IF (lform) THEN
          OPEN(UNIT=outpunit,FILE=filnam)
          READ(outpunit,'(1x,///,1x,"@@@ TIME = ",g11.4)',ERR=199) time
        ELSE 
          OPEN(UNIT=outpunit,FORM='UNFORMATTED',FILE=filnam)
          READ(outpunit,ERR=199) time4
          time = REAL(time4,dbl)
        END IF

        WRITE(logunit,fmt="('  from outp: recovering file ',A20)") filnam
        WRITE(logunit,*)  'time = ', time

      END IF
!
      CALL bcast_real(time, 1, root)
!
      IF( lform .AND. mpime == root ) READ(outpunit,122)
      CALL read_inner_array( outpunit, p, sgl, lform )  ! gas_pressure

      CALL read_inner_array( outpunit, p, sgl, lform )  ! gas_pressure

      IF (job_type == JOB_TYPE_2D) THEN

        IF( lform .AND. mpime == root ) READ(outpunit,122,ERR=199)
        CALL read_inner_array( outpunit, ug, sgl, lform ) ! gas_velocity_r

        IF( lform .AND. mpime == root ) READ(outpunit,122,ERR=199)
        CALL read_inner_array( outpunit, wg, sgl, lform ) ! gas_velocity_z

      ELSE IF (job_type == JOB_TYPE_3D) THEN

        IF( lform .AND. mpime == root ) READ(outpunit,122,ERR=199)
        CALL read_inner_array( outpunit, ug, sgl, lform ) ! gas_velocity_x

        IF( lform .AND. mpime == root ) READ(outpunit,122,ERR=199)
        CALL read_inner_array( outpunit, vg, sgl, lform ) ! gas_velocity_y

        IF( lform .AND. mpime == root ) READ(outpunit,122,ERR=199)
        CALL read_inner_array( outpunit, wg, sgl, lform ) ! gas_velocity_z

      ELSE
        CALL error('outp_','Unknown job type',1)
      END IF

      IF( lform .AND. mpime == root ) READ(outpunit,122,ERR=199)
      CALL read_inner_array( outpunit, tg, sgl, lform )  ! gas_temperature
!
      ALLOCATE( otmp( SIZE( xgc, 1 ) ) )
      DO ig=1,ngas
          IF( lform .AND. mpime == root ) READ(outpunit,122,ERR=199)
          otmp = xgc(:,ig)
          CALL read_inner_array( outpunit, otmp, sgl, lform )  ! gc_molar_fraction
          xgc(:,ig) = otmp
      END DO
      DEALLOCATE( otmp )
!
      ALLOCATE( otmp( SIZE( rlk, 1 ) ) )

      DO is = 1, nsolid

        IF( lform .AND. mpime == root ) READ(outpunit,122,ERR=199)
        otmp = rlk(:,is) * inrl(is)
        CALL read_inner_array( outpunit, otmp, sgl, lform )  ! solid_bulk_density
        rlk(:,is) = otmp * rl(is)

        IF (job_type == JOB_TYPE_2D) THEN

          IF( lform .AND. mpime == root ) READ(outpunit,122,ERR=199)
          CALL read_inner_array( outpunit, us(:,is), sgl, lform )  ! solid_velocity_r

          IF( lform .AND. mpime == root ) READ(outpunit,122,ERR=199)
          CALL read_inner_array( outpunit, ws(:,is), sgl, lform )  ! solid_velocity_z

        ELSE IF (job_type == JOB_TYPE_3D) THEN

          IF( lform .AND. mpime == root ) READ(outpunit,122,ERR=199)
          CALL read_inner_array( outpunit, us(:,is), sgl, lform )  ! solid_velocity_x

          IF( lform .AND. mpime == root ) READ(outpunit,122,ERR=199)
          CALL read_inner_array( outpunit, vs(:,is), sgl, lform )  ! solid_velocity_y

          IF( lform .AND. mpime == root ) READ(outpunit,122,ERR=199)
          CALL read_inner_array( outpunit, ws(:,is), sgl, lform )  ! solid_velocity_z

        END IF

        IF( lform .AND. mpime == root ) READ(outpunit,122,ERR=199)
        CALL read_inner_array( outpunit, ts(:,is), sgl, lform )  ! solid_temperature

      END DO

      DEALLOCATE( otmp )

      IF( mpime == root ) THEN
        CLOSE (outpunit)
      END IF
!
      nf = nf + 1
!
 122  FORMAT(1x,//,6x,/)

      RETURN
!
 199  CALL error('io_restart.f', 'error in reading outputunit', outpunit)
!
      END SUBROUTINE outp_remap
!----------------------------------------------------------------------
      END MODULE io_restart
!----------------------------------------------------------------------
