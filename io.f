!----------------------------------------------------------------------
      MODULE io_restart
!----------------------------------------------------------------------

      USE kinds
      USE eos_gas, ONLY: rgpgc, xgc, ygc
      USE gas_solid_density, ONLY: rog, rgp, rlk
      USE gas_solid_velocity, ONLY: ug, vg, wg
      USE gas_solid_velocity, ONLY: us, vs, ws
      USE gas_solid_temperature, ONLY: sieg, ts, sies, tg 
      USE gas_constants, ONLY: gas_type
      USE parallel, ONLY: mpime, root
      USE pressure_epsilon, ONLY: ep, p
      USE time_parameters, ONLY: time
      USE control_flags, ONLY: job_type, nfil
      USE domain_decomposition, ONLY: data_collect, data_distribute, ncint
      USE dimensions
!
        IMPLICIT NONE
        PRIVATE

        LOGICAL :: old_restart = .FALSE.
        LOGICAL :: dump_all = .FALSE.

        REAL*8 :: max_seconds

        INTERFACE read_array
           MODULE PROCEDURE read_array_dpara, read_array_sserial
        END INTERFACE

        PUBLIC :: old_restart
        PUBLIC :: taperd, tapewr
        PUBLIC :: write_array, read_array
        PUBLIC :: max_seconds
        PUBLIC :: dump_all

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
        WRITE(9) time, nx, ny, nz, nsolid, ngas, nfil

      END IF
!
! ... store the final values of the main physical variables 
! 

      CALL write_array( 9, p, dbl, lform )
      
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

      CALL write_array( 9, sieg, dbl, lform )

      DO is = 1, nsolid
        CALL write_array( 9, rlk(:,is), dbl, lform )
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

      DO is = 1, nsolid
        CALL write_array( 9, sies(:,is), dbl, lform )
      END DO

      DO ig = 1, ngas
        CALL write_array( 9, ygc(ig,:), dbl, lform )
      END DO
!
      CALL write_array( 9, rgp, dbl, lform )

      CALL write_array( 9, rog, dbl, lform )


      IF( dump_all ) THEN

        CALL write_array( 9, tg, dbl, lform )
      
        DO is = 1, nsolid
          CALL write_array( 9, ts(:,is), dbl, lform )
        END DO

        CALL write_array( 9, ep, dbl, lform )

        DO ig = 1, ngas
          CALL write_array( 9, rgpgc(:,ig), dbl, lform )
        END DO

        DO ig = 1, ngas
          CALL write_array( 9, xgc(ig,:), dbl, lform )
        END DO

      END IF
!
      IF( mpime .EQ. root ) THEN
        CLOSE(9)
      END IF

      RETURN
      END SUBROUTINE tapewr
!----------------------------------------------------------------------

      SUBROUTINE taperd
!
      USE dimensions
      USE domain_decomposition, ONLY: ncint, meshinds
      USE atmospheric_conditions, ONLY: p_atm, t_atm
      USE grid, ONLY: zb, dz
      USE specific_heat_module, ONLY: hcapg
      USE gas_constants, ONLY: tzero, hzerog, gas_type


      IMPLICIT NONE
!
      INTEGER :: i, j, k, ijk, is, imesh
      INTEGER :: nx_, ny_, nz_, nsolid_, ngas_
      INTEGER :: ig, info
      LOGICAL :: lform = .FALSE.
      REAL*8 :: zrif, trif, prif, hc
      REAL*8, ALLOCATABLE :: rtmp( :, : )
      

      IF( mpime .EQ. root ) THEN

        WRITE(6,*) ' reading from restart file '
!
        OPEN(UNIT=9,form='unformatted',FILE='pdac.res')

        READ(9) time, nx_, ny_, nz_, nsolid_, ngas_, nfil

        WRITE(6,*) ' time =  ', time
        WRITE(6,*) ' nx   =  ', nx
        WRITE(6,*) ' ny   =  ', ny
        WRITE(6,*) ' nz   =  ', nz
        WRITE(6,*) ' nsolid   =  ', nsolid
        WRITE(6,*) ' ngas     =  ', ngas
        WRITE(6,*) ' nfil     =  ', nfil

      END IF

      CALL bcast_real(time, 1, root)
      CALL bcast_integer(nx_, 1, root)
      CALL bcast_integer(ny_, 1, root)
      CALL bcast_integer(nz_, 1, root)
      CALL bcast_integer(nsolid_, 1, root)
      CALL bcast_integer(ngas_, 1, root)
      CALL bcast_integer(nfil, 1, root)

      IF( nx_ /= nx ) &
        CALL error(' taperd ',' inconsistent dimension nx ', nx_ )
      IF( ny_ /= ny ) &
        CALL error(' taperd ',' inconsistent dimension ny ', ny_ )
      IF( nz_ /= nz ) &
        CALL error(' taperd ',' inconsistent dimension nz ', nz_ )
      IF( nsolid_ /= nsolid ) &
        CALL error(' taperd ',' inconsistent dimension nsolid ', nsolid_ )
      IF( ngas_ /= ngas .AND. .NOT. old_restart ) &
        CALL error(' taperd ',' inconsistent dimension ngas ', ngas_ )

      !
      !  read pressure

      p = 0.0d0
      CALL read_array( 9, p, dbl, lform )
      IF( ANY( p < 0 ) ) THEN
         WRITE(6,*) 'WARNING reading restart, p < 0'
      END IF

      !
      !  read gas velocity

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

      !
      !  read gas enthalpy

      sieg = 0.0d0
      CALL read_array( 9, sieg, dbl, lform )
      IF( ANY( sieg < 0 ) ) THEN
         WRITE(6,*) 'WARNING reading restart, sieg < 0'
      END IF

      !
      !  read particle density

      rlk = 0.0d0
      DO is = 1, nsolid
        CALL read_array( 9, rlk(:,is), dbl, lform )
      END DO
      IF( ANY( rlk < 0 ) ) THEN
         WRITE(6,*) 'WARNING reading restart, rlk < 0'
      END IF


      !
      !  read particles velocity

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

      !
      !  read particle enthalpy

      sies = 0.0d0
      DO is = 1, nsolid
        CALL read_array( 9, sies(:,is), dbl, lform )
      END DO
      IF( ANY( sies < 0.0d0 ) ) THEN
         WRITE(6,*) 'WARNING reading restart, sies < 0'
      END IF

      !
      !  read gas components

      ygc = 0.0d0
      DO ig = 1, ngas
        CALL read_array( 9, ygc(ig,:), dbl, lform )
      END DO
      IF( ANY( ygc < 0 ) ) THEN
         WRITE(6,*) 'WARNING reading restart, ygc < 0'
      END IF

      rgp = 0.0d0
      CALL read_array( 9, rgp, dbl, lform )
      IF( ANY( rgp < 0 ) ) THEN
         WRITE(6,*) 'WARNING reading restart, rgp < 0'
      END IF

      rog = 0.0d0
      CALL read_array( 9, rog, dbl, lform )
      IF( ANY( rog < 0 ) ) THEN
         WRITE(6,*) 'WARNING reading restart, rog < 0'
      END IF

      IF( dump_all ) THEN

        tg = 0.0d0
        CALL read_array( 9, tg, dbl, lform ) 
        IF( ANY( tg < 0 ) ) THEN
           WRITE(6,*) 'WARNING reading restart, tg < 0'
        END IF

        ts = 0.0d0
        DO is = 1, nsolid
          CALL read_array( 9, ts(:,is), dbl, lform )
        END DO
        IF( ANY( ts < 0 ) ) THEN
           WRITE(6,*) 'WARNING reading restart, ts < 0'
        END IF

        ep = 0.0d0
        CALL read_array( 9, ep, dbl, lform )
        IF( ANY( ep < 0 ) ) THEN
           WRITE(6,*) 'WARNING reading restart, ep < 0'
        END IF

        rgpgc = 0.0d0
        DO ig = 1, ngas
          CALL read_array( 9, rgpgc(:,ig), dbl, lform )
        END DO
        IF( ANY( rgpgc < 0 ) ) THEN
           WRITE(6,*) 'WARNING reading restart, rgpgc < 0'
        END IF

        xgc = 0.0d0
        DO ig = 1, ngas
          CALL read_array( 9, xgc(ig,:), dbl, lform )
        END DO
        IF( ANY( xgc < 0 ) ) THEN
           WRITE(6,*) 'WARNING reading restart, xgc < 0'
        END IF

      END IF

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
10    FORMAT( 5(G14.6E3,1X) )
      RETURN
      END SUBROUTINE

!----------------------------------------------------------------------

      SUBROUTINE read_array_dpara( iunit, array, prec, lform )

      !  This subroutine reads a REAL*8 array ( array ) of
      !  ntot elements from file iunit.
      !  Only root processor read the data, then the array
      !  is distributed across other processors
      !  iunit  (input)  file to be read
      !  array  (output) data read from file and distributed to processors
      !  prec   (input)  precision of the data to be read
      !                  sgl  = single precision data
      !                  dbl  = double precision data
      !  lform  (input)  format of the file 
      !                  .TRUE.  = formatted
      !                  .FALSE. = unformatted

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

      IF( lform ) THEN
         ALLOCATE( io_buf( ntot ), STAT=ierr )
         IF( ierr /= 0 ) &
            CALL error(' read_array ', ' cannot allocate io_buf ', ntot )
         IF( mpime == root ) THEN
            READ(iunit,*) ( io_buf(ijk), ijk = 1, ntot )
         END IF
         CALL data_distribute( io_buf, array, 1, ntot )
         DEALLOCATE( io_buf )
      ELSE
         IF( prec == sgl ) THEN
            ALLOCATE( io_bufs( ntot ), STAT=ierr )
            IF( ierr /= 0 ) &
               CALL error(' read_array ', ' cannot allocate io_buf ', ntot )
            IF( mpime == root ) THEN
               READ(iunit) ( io_bufs(ijk), ijk = 1, ntot )
            END IF
            CALL data_distribute( io_bufs, array, 1, ntot )
            DEALLOCATE( io_bufs )
         ELSE
            ALLOCATE( io_buf( ntot ), STAT=ierr )
            IF( ierr /= 0 ) &
               CALL error(' read_array ', ' cannot allocate io_buf ', ntot )
            IF( mpime == root ) THEN
               READ(iunit) ( io_buf(ijk), ijk = 1, ntot )
            END IF
            CALL data_distribute( io_buf, array, 1, ntot )
            DEALLOCATE( io_buf )
         END IF
      END IF

      RETURN
      END SUBROUTINE

!----------------------------------------------------------------------

      SUBROUTINE read_array_sserial( iunit, sarray, lform )

      !  This subroutine reads a REAL array ( sarray ) of
      !  ntot elements from file iunit
      !  NOTE only root processor read the data
      !  iunit  (input)  file to be read
      !  array  (output) data read from file and distributed to processors
      !  lform  (input)  format of the file 
      !                  .TRUE.  = formatted
      !                  .FALSE. = unformatted

      INTEGER, INTENT(IN) :: iunit
      LOGICAL, INTENT(IN) :: lform
      REAL :: sarray(:)

      INTEGER :: ijk, ierr

      IF( ntot < 1 ) &
        CALL error(' read_array ', ' ntot too small ', ntot )

      IF( lform ) THEN
         IF( mpime == root ) THEN
            READ(iunit,*) ( sarray(ijk), ijk = 1, ntot )
         END IF
      ELSE
         IF( mpime == root ) THEN
            READ(iunit) ( sarray(ijk), ijk = 1, ntot )
         END IF
      END IF

      RETURN
      END SUBROUTINE


!----------------------------------------------------------------------
      END MODULE io_restart
!----------------------------------------------------------------------
