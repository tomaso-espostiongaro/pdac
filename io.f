!----------------------------------------------------------------------
      MODULE io_restart
!----------------------------------------------------------------------

      USE kinds
      USE eos_gas, ONLY: rgpgc, xgc, ygc, cg
      USE gas_solid_density, ONLY: rog, rgp, rlk
      USE gas_solid_velocity, ONLY: ug, vg, wg
      USE gas_solid_velocity, ONLY: us, vs, ws
      USE gas_solid_temperature, ONLY: sieg, ts, sies, tg 
      USE parallel, ONLY: mpime, root
      USE pressure_epsilon, ONLY: ep, p
      USE specific_heat_module, ONLY: cp, ck
      USE time_parameters, ONLY: time
      USE control_flags, ONLY: job_type, nfil
      USE domain_decomposition, ONLY: data_collect, data_distribute, ncint
      USE dimensions
!
        IMPLICIT NONE
        PRIVATE

        LOGICAL :: old_restart = .FALSE.

        REAL*8 :: max_seconds

        PUBLIC :: old_restart
        PUBLIC :: taperd, tapewr
        PUBLIC :: write_array, read_array
        PUBLIC :: max_seconds

!----------------------------------------------------------------------
        CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE tapewr
!
      IMPLICIT NONE
!
      INTEGER :: i, j, k, ijk, is
      INTEGER :: ig
      LOGICAL :: lform = .FALSE.
!
      CHARACTER*13 restart_file
      restart_file = 'pdac.res'
!
      IF( mpime .EQ. root ) THEN
!
        WRITE(6,*) ' writing to restart file '
!
        OPEN(UNIT=9,form='unformatted', FILE = restart_file)
!
        WRITE(9) time, nx, ny, nz, nsolid, nfil

      END IF
!
! ... store the final values of the main physical variables 
! 

      CALL write_array( 9, p, dbl, lform )
      
      DO is = 1, nsolid
        CALL write_array( 9, rlk(:,is), dbl, lform )
      END DO

      CALL write_array( 9, sieg, dbl, lform )

      IF (job_type == '2D') THEN

        CALL write_array( 9, ug, dbl, lform )
        CALL write_array( 9, wg, dbl, lform )

      ELSE IF (job_type == '3D') THEN

        CALL write_array( 9, ug, dbl, lform )
        CALL write_array( 9, vg, dbl, lform )
        CALL write_array( 9, wg, dbl, lform )

      ELSE

          CALL error('tapewr','Unknown job type',1)

      END IF

      DO is = 1, nsolid
        CALL write_array( 9, sies(:,is), dbl, lform )
      END DO

      IF (job_type == '2D') THEN

        DO is = 1, nsolid
          CALL write_array( 9, us(:,is), dbl, lform )
          CALL write_array( 9, ws(:,is), dbl, lform )
        END DO

      ELSE IF (job_type == '3D') THEN

        DO is = 1, nsolid
          CALL write_array( 9, us(:,is), dbl, lform )
          CALL write_array( 9, vs(:,is), dbl, lform )
          CALL write_array( 9, ws(:,is), dbl, lform )
        END DO

      ELSE

          CALL error('tapewr','Unknown job type',1)

      END IF

      DO ig = 1, ngas
        CALL write_array( 9, ygc(ig,:), dbl, lform )
      END DO
!
      CALL write_array( 9, rgp, dbl, lform )

      CALL write_array( 9, rog, dbl, lform )

      CALL write_array( 9, ep, dbl, lform )

      CALL write_array( 9, tg, dbl, lform ) 

      DO is = 1, nsolid
        CALL write_array( 9, ts(:,is), dbl, lform )
      END DO

      DO ig = 1, ngas
        CALL write_array( 9, rgpgc(:,ig), dbl, lform )
      END DO

      DO ig = 1, ngas
        CALL write_array( 9, xgc(ig,:), dbl, lform )
      END DO

      CALL write_array( 9, cg, dbl, lform ) 

      DO is = 1, nsolid
        CALL write_array( 9, ck(is,:), dbl, lform )
      END DO

      DO ig = 1, ngas
        CALL write_array( 9, cp(ig,:), dbl, lform )
      END DO
!
      IF( mpime .EQ. root ) THEN
        CLOSE(9)
      END IF

      RETURN
      END SUBROUTINE

!----------------------------------------------------------------------

      SUBROUTINE taperd
!
      USE dimensions
      IMPLICIT NONE
!
      INTEGER :: i, j, k, ijk, is
      INTEGER :: nx_, ny_, nz_, nsolid_
      INTEGER :: ig
      LOGICAL :: lform = .FALSE.

      IF( mpime .EQ. root ) THEN

        WRITE(6,*) ' reading from restart file '
!
        OPEN(UNIT=9,form='unformatted',FILE='pdac.res')

        IF( old_restart ) THEN
          READ(9) time, nx_, ny_, nz_, nsolid_
        ELSE
          READ(9) time, nx_, ny_, nz_, nsolid_, nfil
        END IF

        WRITE(6,*) ' time =  ', time
        WRITE(6,*) ' nx   =  ', nx
        WRITE(6,*) ' ny   =  ', ny
        WRITE(6,*) ' nz   =  ', nz
        WRITE(6,*) ' nsolid =  ', nsolid
        WRITE(6,*) ' nfil   =  ', nfil

      END IF

      CALL bcast_real(time, 1, root)
      CALL bcast_integer(nx_, 1, root)
      CALL bcast_integer(ny_, 1, root)
      CALL bcast_integer(nz_, 1, root)
      CALL bcast_integer(nsolid_, 1, root)

      IF( nx_ /= nx ) &
        CALL error(' taperd ',' inconsistent dimension nx ', nx_ )
      IF( ny_ /= ny ) &
        CALL error(' taperd ',' inconsistent dimension ny ', ny_ )
      IF( nz_ /= nz ) &
        CALL error(' taperd ',' inconsistent dimension nz ', nz_ )
      IF( nsolid_ /= nsolid ) &
        CALL error(' taperd ',' inconsistent dimension nsolid ', nsolid_ )

      p = 0.0d0
      CALL read_array( 9, p, dbl, lform )
      
      rlk = 0.0d0
      DO is = 1, nsolid
        CALL read_array( 9, rlk(:,is), dbl, lform )
      END DO

      sieg = 0.0d0
      CALL read_array( 9, sieg, dbl, lform )

      IF (job_type == '2D') THEN

        ug = 0.0d0
        wg = 0.0d0
        CALL read_array( 9, ug, dbl, lform )
        CALL read_array( 9, wg, dbl, lform )

      ELSE IF (job_type == '3D') THEN

        ug = 0.0d0
        vg = 0.0d0
        wg = 0.0d0
        CALL read_array( 9, ug, dbl, lform )
        CALL read_array( 9, vg, dbl, lform )
        CALL read_array( 9, wg, dbl, lform )

      ELSE

          CALL error('taperd','Unknown job type',1)

      END IF

      sies = 0.0d0
      DO is = 1, nsolid
        CALL read_array( 9, sies(:,is), dbl, lform )
      END DO

      IF (job_type == '2D') THEN

        us = 0.0d0
        ws = 0.0d0
        DO is = 1, nsolid
          CALL read_array( 9, us(:,is), dbl, lform )
          CALL read_array( 9, ws(:,is), dbl, lform )
        END DO

      ELSE IF (job_type == '3D') THEN

        us = 0.0d0
        vs = 0.0d0
        ws = 0.0d0
        DO is = 1, nsolid
          CALL read_array( 9, us(:,is), dbl, lform )
          CALL read_array( 9, vs(:,is), dbl, lform )
          CALL read_array( 9, ws(:,is), dbl, lform )
        END DO

      ELSE

          CALL error('tapewr','Unknown job type',1)

      END IF

      ygc = 0.0d0
      DO ig = 1, ngas
        CALL read_array( 9, ygc(ig,:), dbl, lform )
      END DO

      rgp = 0.0d0
      CALL read_array( 9, rgp, dbl, lform )

      rog = 0.0d0
      CALL read_array( 9, rog, dbl, lform )

      ep = 0.0d0
      CALL read_array( 9, ep, dbl, lform )

      tg = 0.0d0
      CALL read_array( 9, tg, dbl, lform ) 

      ts = 0.0d0
      DO is = 1, nsolid
        CALL read_array( 9, ts(:,is), dbl, lform )
      END DO

      rgpgc = 0.0d0
      DO ig = 1, ngas
        CALL read_array( 9, rgpgc(:,ig), dbl, lform )
      END DO

      xgc = 0.0d0
      DO ig = 1, ngas
        CALL read_array( 9, xgc(ig,:), dbl, lform )
      END DO

      cg = 0.0d0
      CALL read_array( 9, cg, dbl, lform ) 

      ck = 0.0d0
      DO is = 1, nsolid
        CALL read_array( 9, ck(is,:), dbl, lform )
      END DO

      cp = 0.0d0
      DO ig = 1, ngas
        CALL read_array( 9, cp(ig,:), dbl, lform )
      END DO

      IF( mpime .EQ. root ) THEN
        CLOSE (9)
        WRITE(6,*) ' restart file has been read '
      END IF
!
      RETURN
      END SUBROUTINE

!----------------------------------------------------------------------

      SUBROUTINE write_array( iunit, array, prec, lform )

      INTEGER, INTENT(IN) :: iunit, prec
      LOGICAL, INTENT(IN) :: lform
      REAL*8 :: array(:)

      INTEGER :: ijk, ierr
      REAL*8, ALLOCATABLE :: io_buf(:)

      IF( ntot < 1 ) &
        CALL error(' write_array ', ' ntot too small ', ntot )

      IF( prec /= sgl .AND. prec /= dbl ) &
        CALL error(' write_array ', ' unknown precision ', prec )

      ALLOCATE( io_buf( ntot ), STAT=ierr )
      IF( ierr /= 0 ) &
        CALL error(' write_array ', ' cannot allocate io_buf ', ntot )

      CALL data_collect( io_buf, array, 1, ntot )
      IF( mpime == root ) THEN
         IF( lform ) THEN
           WRITE(iunit,10) ( io_buf(ijk), ijk = 1, ntot )
         ELSE
           IF( prec == sgl ) THEN
             WRITE(iunit) REAL( io_buf( 1 : ntot ), sgl )
           ELSE 
             WRITE(iunit) ( io_buf(ijk), ijk = 1, ntot )
           END IF
         END IF
      END IF
      DEALLOCATE( io_buf )
10    FORMAT( 5(G14.6,1X) )
      RETURN
      END SUBROUTINE

!----------------------------------------------------------------------

      SUBROUTINE read_array( iunit, array, prec, lform )

      INTEGER, INTENT(IN) :: iunit, prec
      LOGICAL, INTENT(IN) :: lform
      REAL*8 :: array(:)

      INTEGER :: ijk, ierr
      REAL*8, ALLOCATABLE :: io_buf(:)
      REAL(sgl), ALLOCATABLE :: io_bufs(:)

      IF( ntot < 1 ) &
        CALL error(' read_array ', ' ntot too small ', ntot )

      IF( prec /= sgl .AND. prec /= dbl ) &
        CALL error(' read_array ', ' unknown precision ', prec )

      IF( prec /= sgl .AND. lform ) &
        CALL error(' read_array ', ' single precision with formatted file ', prec )

      IF( prec == sgl ) THEN
        ALLOCATE( io_bufs( ntot ), STAT=ierr )
      ELSE
        ALLOCATE( io_buf( ntot ), STAT=ierr )
      END IF

      IF( ierr /= 0 ) &
        CALL error(' read_array ', ' cannot allocate io_buf ', ntot )

      IF( mpime == root ) THEN
         IF( lform ) THEN
           READ(iunit,*) ( io_buf(ijk), ijk = 1, ntot )
         ELSE
           IF( prec == sgl ) THEN
             READ(iunit) ( io_bufs(ijk), ijk = 1, ntot )
           ELSE
             READ(iunit) ( io_buf(ijk), ijk = 1, ntot )
           END IF
         END IF
      END IF

      IF( prec == sgl ) THEN
        CALL data_distribute( io_bufs, array, 1, ntot )
      ELSE
        CALL data_distribute( io_buf, array, 1, ntot )
      END IF

      IF( prec == sgl ) THEN
        DEALLOCATE( io_bufs )
      ELSE
        DEALLOCATE( io_buf )
      END IF

      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      END MODULE io_restart
!----------------------------------------------------------------------
