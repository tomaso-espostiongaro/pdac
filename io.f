!----------------------------------------------------------------------
      MODULE io_restart
!----------------------------------------------------------------------

      USE control_flags, ONLY: job_type
      USE dimensions
      USE domain_decomposition, ONLY: data_collect, data_distribute, ncint
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
!
      IMPLICIT NONE
      PRIVATE

      INTEGER :: nfil

      REAL*8 :: max_seconds

      INTERFACE read_array
         MODULE PROCEDURE read_array_dpara, read_array_sserial
      END INTERFACE

      PUBLIC :: taperd, tapewr
      PUBLIC :: write_array, read_array
      PUBLIC :: max_seconds
      PUBLIC :: nfil

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
      
      IF (job_type == '2D') THEN

        CALL write_array( resunit, ug, dbl, lform )
        CALL write_array( resunit, wg, dbl, lform )

      ELSE IF (job_type == '3D') THEN

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

      IF (job_type == '2D') THEN

        DO is = 1, nsolid
          CALL write_array( resunit, us(:,is), dbl, lform )
          CALL write_array( resunit, ws(:,is), dbl, lform )
        END DO

      ELSE IF (job_type == '3D') THEN

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
      USE domain_decomposition, ONLY: ncint, meshinds
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

        READ(resunit) time, nx_, ny_, nz_, nsolid_, ngas_, nfil
        READ(resunit) gas_type

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

      IF (job_type == '2D') THEN

        ug = 0.0d0
        wg = 0.0d0
        CALL read_array( resunit, ug, dbl, lform )
        CALL read_array( resunit, wg, dbl, lform )

      ELSE IF (job_type == '3D') THEN

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

      IF (job_type == '2D') THEN

        us = 0.0d0
        ws = 0.0d0
        DO is = 1, nsolid
          CALL read_array( resunit, us(:,is), dbl, lform )
          CALL read_array( resunit, ws(:,is), dbl, lform )
        END DO

      ELSE IF (job_type == '3D') THEN

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
         IF( mpime == root ) WRITE(errorunit,*) 'WARNING reading restart, ts < 0'
      END IF

      IF( mpime == root ) THEN
        CLOSE (resunit)
        WRITE(logunit,*) ' restart file has been red '
      END IF
!
      RETURN
      END SUBROUTINE

!----------------------------------------------------------------------

      SUBROUTINE write_array( iunit, array, prec, lform )

      INTEGER, INTENT(IN) :: iunit, prec
      LOGICAL, INTENT(IN) :: lform
      REAL*8 :: array(:)

      INTEGER :: ijk, ij, k, ierr, first
      REAL*8, ALLOCATABLE :: io_buf(:)
      REAL, ALLOCATABLE :: io_buf_sgl(:)

      IF( ntot < 1 ) &
        CALL error(' write_array ', ' ntot too small ', ntot )

      IF( prec /= sgl .AND. prec /= dbl ) &
        CALL error(' write_array ', ' unknown precision ', prec )

      ALLOCATE( io_buf( ntot ), STAT=ierr )
      ALLOCATE( io_buf_sgl( ntot ), STAT=ierr )
      IF( ierr /= 0 ) &
        CALL error(' write_array ', ' cannot allocate io_buf ', ntot )

      CALL data_collect( io_buf, array, 1, ntot )

      IF( mpime == root ) THEN
         IF( lform ) THEN
           !WRITE(iunit,10) ( io_buf(ijk), ijk = 1, ntot )
           !
           ! ... This format is consistent with the old PDAC2D
           ! ... OUTPUT format
           DO k = 1, nz
             first = (k-1) * ntr + 1
             WRITE(iunit,10) ( io_buf(ij), ij = first, first + ntr - 1 )
           END DO
         ELSE
           IF( prec == sgl ) THEN
             io_buf_sgl(:) = REAL( io_buf(:), sgl)
             WRITE(iunit) ( io_buf_sgl(ijk), ijk = 1, ntot )
           ELSE 
             WRITE(iunit) ( io_buf(ijk), ijk = 1, ntot )
           END IF
         END IF
      END IF

      DEALLOCATE( io_buf )
10    FORMAT(1x,10(1x,G14.6E3))
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
