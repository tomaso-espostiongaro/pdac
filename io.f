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
      USE control_flags, ONLY: job_type
      USE domain_decomposition, ONLY: data_collect, data_distribute, ncint
      USE dimensions
!
        IMPLICIT NONE
        PRIVATE

        PUBLIC :: taperd, tapewr
        PUBLIC :: write_array, read_array

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
      !WRITE(6, fmt = "( /,'Subroutine TAPEWR')" )
!
      IF( mpime .EQ. root ) THEN
!
        OPEN(UNIT=9,form='unformatted', FILE = restart_file)
!
        WRITE(9) time, nx, ny, nz, nsolid

      END IF
!
! ... store the final values of the main physical variables 
! 

      !WRITE(6, fmt = "( I3,':',' gas_pressure = ',D18.8)" ) mpime, SUM( p(1:ncint) )
      CALL write_array( 9, p, dbl, lform )
      
      !WRITE(6, fmt = "( I3,':',' solid_bulk_density = ',D18.8)" ) mpime, SUM( rlk(1:ncint,1:nsolid) )
      DO is = 1, nsolid
        CALL write_array( 9, rlk(:,is), dbl, lform )
      END DO

      !WRITE(6, fmt = "( I3,':',' gas_enthalpy = ',D18.8)" ) mpime, SUM( sieg(1:ncint) )
      CALL write_array( 9, sieg, dbl, lform )

      IF (job_type == '2D') THEN

        !WRITE(6, fmt = "( I3,':',' gas_velocity_r = ',D18.8)" ) mpime, SUM( ug(1:ncint) )
        !WRITE(6, fmt = "( I3,':',' gas_velocity_z = ',D18.8)" ) mpime, SUM( wg(1:ncint) )
        CALL write_array( 9, ug, dbl, lform )
        CALL write_array( 9, wg, dbl, lform )

      ELSE IF (job_type == '3D') THEN

        !WRITE(6, fmt = "( I3,':',' gas_velocity_x = ',D18.8)" ) mpime, SUM( ug(1:ncint) )
        !WRITE(6, fmt = "( I3,':',' gas_velocity_y = ',D18.8)" ) mpime, SUM( vg(1:ncint) )
        !WRITE(6, fmt = "( I3,':',' gas_velocity_z = ',D18.8)" ) mpime, SUM( wg(1:ncint) )
        CALL write_array( 9, ug, dbl, lform )
        CALL write_array( 9, vg, dbl, lform )
        CALL write_array( 9, wg, dbl, lform )

      ELSE

          CALL error('tapewr','Unknown job type',1)

      END IF

      !WRITE(6, fmt = "( I3,':',' solid_enthalpy = ',D18.8)" ) mpime, SUM( sies(1:ncint,1:nsolid) )
      DO is = 1, nsolid
        CALL write_array( 9, sies(:,is), dbl, lform )
      END DO

      IF (job_type == '2D') THEN

        !WRITE(6, fmt = "( I3,':',' solid_velocity_r = ',D18.8)" ) mpime, SUM( us(1:ncint,1:nsolid) )
        !WRITE(6, fmt = "( I3,':',' solid_velocity_z = ',D18.8)" ) mpime, SUM( ws(1:ncint,1:nsolid) )
        DO is = 1, nsolid
          CALL write_array( 9, us(:,is), dbl, lform )
          CALL write_array( 9, ws(:,is), dbl, lform )
        END DO

      ELSE IF (job_type == '3D') THEN

        !WRITE(6, fmt = "( I3,':',' solid_velocity_x = ',D18.8)" ) mpime, SUM( us(1:ncint,1:nsolid) )
        !WRITE(6, fmt = "( I3,':',' solid_velocity_y = ',D18.8)" ) mpime, SUM( vs(1:ncint,1:nsolid) )
        !WRITE(6, fmt = "( I3,':',' solid_velocity_z = ',D18.8)" ) mpime, SUM( ws(1:ncint,1:nsolid) )
        DO is = 1, nsolid
          CALL write_array( 9, us(:,is), dbl, lform )
          CALL write_array( 9, vs(:,is), dbl, lform )
          CALL write_array( 9, ws(:,is), dbl, lform )
        END DO

      ELSE

          CALL error('tapewr','Unknown job type',1)

      END IF

      !WRITE(6, fmt = "( I3,':',' gc_mass_fraction = ',D18.8)" ) mpime, SUM( ygc(1:ngas,1:ncint) )
      DO ig = 1, ngas
        CALL write_array( 9, ygc(ig,:), dbl, lform )
      END DO
!
      !WRITE(6, fmt = "( I3,':',' gas_bulk_density = ',D18.8)" ) mpime, SUM( rgp(1:ncint) )
      CALL write_array( 9, rgp, dbl, lform )

      !WRITE(6, fmt = "( I3,':',' gas_density = ',D18.8)" ) mpime, SUM( rog(1:ncint) )
      CALL write_array( 9, rog, dbl, lform )

      !WRITE(6, fmt = "( I3,':',' void_fraction = ',D18.8)" ) mpime, SUM( ep(1:ncint) )
      CALL write_array( 9, ep, dbl, lform )

      !WRITE(6, fmt = "( I3,':',' gas_temperature = ',D18.8)" ) mpime, SUM( tg(1:ncint) )
      CALL write_array( 9, tg, dbl, lform ) 

      !WRITE(6, fmt = "( I3,':',' solid_temperature = ',D18.8)" ) mpime, SUM( ts(1:ncint,1:nsolid) )
      DO is = 1, nsolid
        CALL write_array( 9, ts(:,is), dbl, lform )
      END DO

      !WRITE(6, fmt = "( I3,':',' gc_bulk_density = ',D18.8)" ) mpime, SUM( rgpgc(1:ncint,1:ngas) )
      DO ig = 1, ngas
        CALL write_array( 9, rgpgc(:,ig), dbl, lform )
      END DO

      !WRITE(6, fmt = "( I3,':',' gc_molar_fraction = ',D18.8)" ) mpime, SUM( xgc(1:ngas,1:ncint) )
      DO ig = 1, ngas
        CALL write_array( 9, xgc(ig,:), dbl, lform )
      END DO

      !WRITE(6, fmt = "( I3,':',' gas_specific_heat = ',D18.8)" ) mpime, SUM( cg(1:ncint) )
      CALL write_array( 9, cg, dbl, lform ) 

      !WRITE(6, fmt = "( I3,':',' solid_specific_heat = ',D18.8)" ) mpime, SUM( ck(1:nsolid,1:ncint) )
      DO is = 1, nsolid
        CALL write_array( 9, ck(is,:), dbl, lform )
      END DO

      !WRITE(6, fmt = "( I3,':',' gc_specific_heat = ',D18.8)" ) mpime, SUM( cp(1:ngas,1:ncint) )
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

      !WRITE(6, fmt = "( /,'Subroutine TAPERD')" )

      IF( mpime .EQ. root ) THEN
!
        OPEN(UNIT=9,form='unformatted',FILE='pdac.res')
        READ(9) time, nx_, ny_, nz_, nsolid_

      END IF

      CALL bcast_integer(nx_, 1, root)
      CALL bcast_integer(ny_, 1, root)
      CALL bcast_integer(nz_, 1, root)
      CALL bcast_integer(nsolid_, 1, root)

      !WRITE(6, fmt = "(I3,':','dimensions = ',8I6)" ) mpime, nx, ny, nz, nx_, ny_, nz_
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
      !WRITE(6, fmt = "(I3,':',' gas_pressure = ',D18.8)" ) mpime, SUM( p(1:ncint) )
      
      rlk = 0.0d0
      DO is = 1, nsolid
        CALL read_array( 9, rlk(:,is), dbl, lform )
      END DO
      !WRITE(6, fmt = "( I3,':',' solid_bulk_density = ',D18.8)" ) mpime, SUM( rlk(1:ncint,1:nsolid) )

      sieg = 0.0d0
      CALL read_array( 9, sieg, dbl, lform )
      !WRITE(6, fmt = "( I3,':',' gas_enthalpy = ',D18.8)" ) mpime, SUM( sieg(1:ncint) )

      IF (job_type == '2D') THEN

        ug = 0.0d0
        wg = 0.0d0
        CALL read_array( 9, ug, dbl, lform )
        CALL read_array( 9, wg, dbl, lform )
        !WRITE(6, fmt = "( I3,':',' gas_velocity_r = ',D18.8)" ) mpime, SUM( ug(1:ncint) )
        !WRITE(6, fmt = "( I3,':',' gas_velocity_z = ',D18.8)" ) mpime, SUM( wg(1:ncint) )

      ELSE IF (job_type == '3D') THEN

        ug = 0.0d0
        vg = 0.0d0
        wg = 0.0d0
        CALL read_array( 9, ug, dbl, lform )
        CALL read_array( 9, vg, dbl, lform )
        CALL read_array( 9, wg, dbl, lform )
        !WRITE(6, fmt = "( I3,':',' gas_velocity_x = ',D18.8)" ) mpime, SUM( ug(1:ncint) )
        !WRITE(6, fmt = "( I3,':',' gas_velocity_y = ',D18.8)" ) mpime, SUM( vg(1:ncint) )
        !WRITE(6, fmt = "( I3,':',' gas_velocity_z = ',D18.8)" ) mpime, SUM( wg(1:ncint) )

      ELSE

          CALL error('taperd','Unknown job type',1)

      END IF

      sies = 0.0d0
      DO is = 1, nsolid
        CALL read_array( 9, sies(:,is), dbl, lform )
      END DO
      !WRITE(6, fmt = "( I3,':',' solid_enthalpy = ',D18.8)" ) mpime, SUM( sies(1:ncint,1:nsolid) )


      IF (job_type == '2D') THEN

        us = 0.0d0
        ws = 0.0d0
        DO is = 1, nsolid
          CALL read_array( 9, us(:,is), dbl, lform )
          CALL read_array( 9, ws(:,is), dbl, lform )
        END DO
        !WRITE(6, fmt = "( I3,':',' solid_velocity_r = ',D18.8)" ) mpime, SUM( us(1:ncint,1:nsolid) )
        !WRITE(6, fmt = "( I3,':',' solid_velocity_z = ',D18.8)" ) mpime, SUM( ws(1:ncint,1:nsolid) )

      ELSE IF (job_type == '3D') THEN

        us = 0.0d0
        vs = 0.0d0
        ws = 0.0d0
        DO is = 1, nsolid
          CALL read_array( 9, us(:,is), dbl, lform )
          CALL read_array( 9, vs(:,is), dbl, lform )
          CALL read_array( 9, ws(:,is), dbl, lform )
        END DO
        !WRITE(6, fmt = "( I3,':',' solid_velocity_x = ',D18.8)" ) mpime, SUM( us(1:ncint,1:nsolid) )
        !WRITE(6, fmt = "( I3,':',' solid_velocity_y = ',D18.8)" ) mpime, SUM( vs(1:ncint,1:nsolid) )
        !WRITE(6, fmt = "( I3,':',' solid_velocity_z = ',D18.8)" ) mpime, SUM( ws(1:ncint,1:nsolid) )

      ELSE

          CALL error('tapewr','Unknown job type',1)

      END IF

      ygc = 0.0d0
      DO ig = 1, ngas
        CALL read_array( 9, ygc(ig,:), dbl, lform )
      END DO
      !WRITE(6, fmt = "( I3,':',' gc_mass_fraction = ',D18.8)" ) mpime, SUM( ygc(1:ngas,1:ncint) )

      rgp = 0.0d0
      CALL read_array( 9, rgp, dbl, lform )
      !WRITE(6, fmt = "( I3,':',' gas_bulk_density = ',D18.8)" ) mpime, SUM( rgp(1:ncint) )

      rog = 0.0d0
      CALL read_array( 9, rog, dbl, lform )
      !WRITE(6, fmt = "( I3,':',' gas_density = ',D18.8)" ) mpime, SUM( rog(1:ncint) )

      ep = 0.0d0
      CALL read_array( 9, ep, dbl, lform )
      !WRITE(6, fmt = "( I3,':',' void_fraction = ',D18.8)" ) mpime, SUM( ep(1:ncint) )

      tg = 0.0d0
      CALL read_array( 9, tg, dbl, lform ) !  READ(9) (gas_temperature(ijk),ijk=1,ntot)
      !WRITE(6, fmt = "( I3,':',' gas_temperature = ',D18.8)" ) mpime, SUM( tg(1:ncint) )

      ts = 0.0d0
      DO is = 1, nsolid
        CALL read_array( 9, ts(:,is), dbl, lform )
      END DO
      !WRITE(6, fmt = "( I3,':',' solid_temperature = ',D18.8)" ) mpime, SUM( ts(1:ncint,1:nsolid) )

      rgpgc = 0.0d0
      DO ig = 1, ngas
        CALL read_array( 9, rgpgc(:,ig), dbl, lform )
      END DO
      !WRITE(6, fmt = "( I3,':',' gc_bulk_density = ',D18.8)" ) mpime, SUM( rgpgc(1:ncint,1:ngas) )

      xgc = 0.0d0
      DO ig = 1, ngas
        CALL read_array( 9, xgc(ig,:), dbl, lform )
      END DO
      !WRITE(6, fmt = "( I3,':',' gc_molar_fraction = ',D18.8)" ) mpime, SUM( xgc(1:ngas,1:ncint) )

      cg = 0.0d0
      CALL read_array( 9, cg, dbl, lform ) 
      !WRITE(6, fmt = "( I3,':',' gas_specific_heat = ',D18.8)" ) mpime, SUM( cg(1:ncint) )

      ck = 0.0d0
      DO is = 1, nsolid
        CALL read_array( 9, ck(is,:), dbl, lform )
      END DO
      !WRITE(6, fmt = "( I3,':',' solid_specific_heat = ',D18.8)" ) mpime, SUM( ck(1:nsolid,1:ncint) )

      cp = 0.0d0
      DO ig = 1, ngas
        CALL read_array( 9, cp(ig,:), dbl, lform )
      END DO
      !WRITE(6, fmt = "( I3,':',' gc_specific_heat = ',D18.8)" ) mpime, SUM( cp(1:ngas,1:ncint) )

      IF( mpime .EQ. root ) THEN
        CLOSE (9)
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
