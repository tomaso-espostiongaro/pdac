!----------------------------------------------------------------------
      MODULE output_dump
!----------------------------------------------------------------------
      USE kinds
      USE io_restart, ONLY: nfil

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
      USE eos_gas, ONLY: xgc
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
!
      IMPLICIT NONE
!
      CHARACTER :: filnam*11
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
          OPEN(UNIT=12,FILE=filnam)
          WRITE(12,*) time
        ELSE 
          OPEN(UNIT=12,FORM='UNFORMATTED',FILE=filnam)
          WRITE(12) REAL(time,4)
        END IF
 
      END IF
!
      CALL write_array( 12, p, sgl, lform )  ! gas_pressure
      CALL write_array( 12, tg, sgl, lform )   ! gas_temperature
      CALL write_array( 12, rog, sgl, lform )  ! gas_density
!
      IF (job_type == '2D') THEN
        CALL write_array( 12, ug, sgl, lform ) ! gas_velocity_r
        CALL write_array( 12, wg, sgl, lform ) ! gas_velocity_z
      ELSE IF (job_type == '3D') THEN
        CALL write_array( 12, ug, sgl, lform ) ! gas_velocity_x
        CALL write_array( 12, vg, sgl, lform ) ! gas_velocity_y
        CALL write_array( 12, wg, sgl, lform ) ! gas_velocity_z
      ELSE
        CALL error('outp_','Unknown job type',1)
      END IF

      IF( mpime == root ) THEN
        CLOSE (12)
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

!
      IMPLICIT NONE
!
      CHARACTER( LEN = 15 ) :: filnam
      CHARACTER( LEN = 4 ) :: lettera
!
      INTEGER :: i,j,ijk
      INTEGER :: ig
      REAL*8 :: energy
!
      nfil=nfil+1
      filnam='shtube.'//lettera(nfil)

      IF( mpime .EQ. root ) THEN

      OPEN(UNIT=11,FILE=filnam)
!
      DO ijk = 1, ncint
      IF (flag(ijk) == 1) THEN
        energy = p(ijk)/(rog(ijk)*(gammaair - 1.D0))
        WRITE(11,550)rog(ijk),wg(ijk),p(ijk),energy,tg(ijk)
      END IF
      END DO
!
      CLOSE (11)
      END IF
!
      RETURN

 550  FORMAT(1x,10(1x,g12.6))
!
      END SUBROUTINE shock_tube_out
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
!
      IMPLICIT NONE
!
      CHARACTER( LEN = 11 ) :: filnam
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
          OPEN(UNIT=12,FILE=filnam)
          WRITE(12,'(1x,///,1x,"@@@ TIME = ",g11.4)') time
        ELSE 
          OPEN(UNIT=12,FORM='UNFORMATTED',FILE=filnam)
          WRITE(12) REAL(time,4)
        END IF

        WRITE(6,fmt="('  from outp: writing file ',A20)") filnam
 
      END IF
!
      IF( lform .AND. mpime == root ) WRITE(12,'(1x,//,1x,"P",/)')
      CALL write_array( 12, p, sgl, lform )  ! gas_pressure

      IF (job_type == '2D') THEN

        IF( lform .AND. mpime == root ) WRITE(12,'(1x,//,1x,"UG",/)')
        CALL write_array( 12, ug, sgl, lform ) ! gas_velocity_r

        IF( lform .AND. mpime == root ) WRITE(12,'(1x,//,1x,"WG",/)')
        CALL write_array( 12, wg, sgl, lform ) ! gas_velocity_z

      ELSE IF (job_type == '3D') THEN

        IF( lform .AND. mpime == root ) WRITE(12,'(1x,//,1x,"UG",/)')
        CALL write_array( 12, ug, sgl, lform ) ! gas_velocity_x

        IF( lform .AND. mpime == root ) WRITE(12,'(1x,//,1x,"VG",/)')
        CALL write_array( 12, vg, sgl, lform ) ! gas_velocity_y

        IF( lform .AND. mpime == root ) WRITE(12,'(1x,//,1x,"WG",/)')
        CALL write_array( 12, wg, sgl, lform ) ! gas_velocity_z

      ELSE
        CALL error('outp_','Unknown job type',1)
      END IF

      IF( lform .AND. mpime == root ) WRITE(12,'(1x,//,1x,"TG",/)')
      CALL write_array( 12, tg, sgl, lform )  ! gas_temperature
!
      ALLOCATE( otmp( SIZE( xgc, 2 ) ) )
      DO ig=1,ngas
          otmp = xgc(ig,:)
          IF( lform .AND. mpime == root ) WRITE(12,'(1x,//,1x,"XGC",I1,/)') ig
          CALL write_array( 12, otmp, sgl, lform )  ! gc_molar_fraction
      END DO
      DEALLOCATE( otmp )
!
      ALLOCATE( otmp( SIZE( rlk, 1 ) ) )

      DO is = 1, nsolid

        otmp = rlk(:,is)*inrl(is)
        IF( lform .AND. mpime == root ) WRITE(12,'(1x,//,1x,"EPS",I1,/)') is
        CALL write_array( 12, otmp, sgl, lform )  ! solid_bulk_density

        IF (job_type == '2D') THEN

          IF( lform .AND. mpime == root ) WRITE(12,'(1x,//,1x,"US",I1,/)') is
          CALL write_array( 12, us(:,is), sgl, lform )  ! solid_velocity_r

          IF( lform .AND. mpime == root ) WRITE(12,'(1x,//,1x,"WS",I1,/)') is
          CALL write_array( 12, ws(:,is), sgl, lform )  ! solid_velocity_z

        ELSE IF (job_type == '3D') THEN

          IF( lform .AND. mpime == root ) WRITE(12,'(1x,//,1x,"US",I1,/)') is
          CALL write_array( 12, us(:,is), sgl, lform )  ! solid_velocity_x

          IF( lform .AND. mpime == root ) WRITE(12,'(1x,//,1x,"VS",I1,/)') is
          CALL write_array( 12, vs(:,is), sgl, lform )  ! solid_velocity_y

          IF( lform .AND. mpime == root ) WRITE(12,'(1x,//,1x,"WS",I1,/)') is
          CALL write_array( 12, ws(:,is), sgl, lform )  ! solid_velocity_z

        END IF

        IF( lform .AND. mpime == root ) WRITE(12,'(1x,//,1x,"TS",I1,/)') is
        CALL write_array( 12, ts(:,is), sgl, lform )  ! solid_temperature

      END DO

      DEALLOCATE( otmp )

      IF( mpime == root ) THEN
        CLOSE (12)
      END IF

      nfil=nfil+1
!
      RETURN
      END SUBROUTINE outp
!----------------------------------------------------------------------
      SUBROUTINE filter_outp( irest )
!
      USE dimensions, ONLY: nsolid, ngas, nx, ny, nz, ntot
      USE io_restart, ONLY: read_array
      USE parallel, ONLY: nproc, mpime, root, group
      USE control_flags, ONLY: job_type
      USE domain_decomposition, ONLY: ncint, meshinds
      USE gas_constants, ONLY: gas_type, gammaair
      USE grid, ONLY: dx, dy, dz
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: irest
      CHARACTER(LEN = 11) :: filnam
      CHARACTER(LEN = 4 ) :: lettera
      CHARACTER(LEN = 2 ) :: lettera2
      CHARACTER(LEN = 4 ) :: var
      LOGICAL :: lform
!
      INTEGER :: ig, is
      INTEGER :: iunit
      REAL*8 :: time
      REAL   :: stime

      REAL, ALLOCATABLE :: array(:)
!
!
      filnam='output.'//lettera(irest)
      WRITE(6,fmt="('  filter: reading file ',A11)") filnam

      lform = formatted_output

      iunit = 50

      IF( mpime == root ) THEN

        IF (lform) THEN
          OPEN( UNIT=12, FILE=filnam, STATUS='OLD')
          READ (12,*) time
        ELSE 
          OPEN(UNIT=12,FORM='UNFORMATTED',FILE=filnam, STATUS='OLD')
          READ (12) stime
        END IF

      END IF

      ALLOCATE( array( ntot ) )

      WRITE(6,fmt="('  filtering gas pressure ')")
!
      CALL read_array( 12, array, lform )  ! gas_pressure

      CALL crop_array( 'pgas' )  ! gas_pressure

      WRITE(6,fmt="('  filtering reading gas velocities ')")

      IF (job_type == '2D') THEN

        CALL read_array( 12, array, lform ) ! gas_velocity_r
        CALL inte_array_x( 'ugas' )  

        CALL read_array( 12, array, lform ) ! gas_velocity_z
        CALL inte_array_z( 'wgas' ) 

      ELSE IF (job_type == '3D') THEN

        CALL read_array( 12, array, lform ) ! gas_velocity_x
        CALL inte_array_x( 'ugas' )  

        CALL read_array( 12, array, lform ) ! gas_velocity_y
        CALL inte_array_y( 'vgas' )  

        CALL read_array( 12, array, lform ) ! gas_velocity_z
        CALL inte_array_z( 'wgas' )  

      ELSE
        CALL error('outp_','Unknown job type',1)
      END IF

      WRITE(6,fmt="('  filtering gas temperature ')")

      CALL read_array( 12, array, lform )  ! gas_temperature
      CALL crop_array( 'tgas' )  

      WRITE(6,fmt="('  filtering molarfraction ')")
!
      DO ig=1,ngas
          var = 'xg'//lettera2( ig )
          CALL read_array( 12, array, lform )  ! gc_molar_fraction
          CALL crop_array( var )  
      END DO

      WRITE( 6, fmt="('  filtering solid density, velocities and temperature')")
!
      DO is = 1, nsolid
        CALL read_array( 12, array, lform )  ! solid_bulk_density
        var = 'ep'//lettera2( is )
        CALL crop_array( var )  
        IF (job_type == '2D') THEN
          CALL read_array( 12, array, lform )  ! solid_velocity_r
          var = 'us'//lettera2( is )
          CALL inte_array_x( var )  
          CALL read_array( 12, array, lform )  ! solid_velocity_z
          var = 'ws'//lettera2( is )
          CALL inte_array_z( var )  
        ELSE IF (job_type == '3D') THEN
          CALL read_array( 12, array, lform )  ! solid_velocity_x
          var = 'us'//lettera2( is )
          CALL inte_array_x( var )  
          CALL read_array( 12, array, lform )  ! solid_velocity_y
          var = 'vs'//lettera2( is )
          CALL inte_array_y( var )  
          CALL read_array( 12, array, lform )  ! solid_velocity_z
          var = 'ws'//lettera2( is )
          CALL inte_array_z( var )  
        END IF
        CALL read_array( 12, array, lform )  ! solid_temperature
        var = 'ts'//lettera2( is )
        CALL crop_array( var )  
      END DO

      DEALLOCATE( array )

      IF( mpime == root ) THEN
        CLOSE (12)
      END IF
!
      RETURN
!-----------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------
      SUBROUTINE crop_array( var )
        CHARACTER(LEN = 16) :: filwri
        CHARACTER(LEN = 4) :: var
        INTEGER :: ii, i, j, k, ijk, imesh
        REAL, ALLOCATABLE :: sarray(:)
        filwri = 'filter.' // var // '.' // lettera(irest)
        IF (lform) THEN
          OPEN( UNIT=iunit, FILE=filwri, STATUS='UNKNOWN' )
        ELSE 
          OPEN( UNIT=iunit, FORM='UNFORMATTED', FILE=filwri, STATUS='UNKNOWN' )
        END IF
        WRITE(6,fmt="('  crop_array: writing file ',A16)") filwri

        ALLOCATE( sarray( ntot ) )
        sarray = 0.0
        ii = 1

        IF( job_type == '2D' ) THEN
          DO k = 2, nz-1
             DO i = 2, nx-1
                imesh = i + ( k-1 ) * nx 
                sarray( ii ) = array( imesh )
                ii = ii + 1
             END DO
          END DO
        ELSE 
          DO k = 2, nz-1
             DO j = 2, ny-1
                DO i = 2, nx-1
                   imesh = i + ( j-1 ) * nx + ( k-1 ) * nx * ny
                   sarray( ii ) = array( imesh )
                   ii = ii + 1
                END DO
             END DO
          END DO
        END IF

        IF( mpime == root ) THEN
           IF( lform ) THEN
             WRITE(iunit,10) ( sarray(ijk), ijk = 1, (ii-1) )
           ELSE
             WRITE(iunit) sarray( 1 : (ii-1) )
           END IF
        END IF

        DEALLOCATE( sarray )

        CLOSE( iunit )

10      FORMAT( 5(G14.6,1X) )
      END SUBROUTINE crop_array
!-----------------------------------------------------------------------
      SUBROUTINE inte_array_x( var )
        CHARACTER(LEN = 16) :: filwri
        CHARACTER(LEN = 4) :: var
        INTEGER :: ii, i, j, k, ijk, imesh
        REAL*8 :: u1, u2, u3, u4, uu1, uu2, uu3, uu4, s, r
        REAL, ALLOCATABLE :: sarray(:)
        filwri='filter.' // var // '.' // lettera(irest)
        IF (lform) THEN
          OPEN( UNIT=iunit, FILE=filwri, STATUS='UNKNOWN' )
        ELSE 
          OPEN( UNIT=iunit, FORM='UNFORMATTED', FILE=filwri, STATUS='UNKNOWN' )
        END IF
        WRITE(6,fmt="('  inte_array_x: writing file ',A16)") filwri
        ALLOCATE( sarray( ntot ) )
        imesh = 0
        sarray = 0.0
        ii = 1
        DO k = 1, nz-1
        DO j = 1, ny-1
        DO i = 1, nx-1
           imesh = i + ( j-1 ) * nx + ( k-1 ) * nx * ny
           u1 = array(imesh)
           imesh = i + ( j-1+1 ) * nx + ( k-1 ) * nx * ny
           u2 = array(imesh)
           imesh = i + ( j-1 ) * nx + ( k-1+1 ) * nx * ny
           u3 = array(imesh)
           imesh = i + ( j-1+1 ) * nx + ( k-1+1 ) * nx * ny
           u4 = array(imesh)

           r = dy( j ) / ( dy( j ) + dy( j + 1 ) )
           s = dz( k ) / ( dz( k ) + dz( k + 1 ) )

	   uu1 = (1-r)*(1-s);
	   uu2 = r*(1-s);
	   uu3 = (1-r)*s;
	   uu4 = r*s;
	 
	   sarray( ii ) = u1 * uu1 + u2 * uu2 + u3 * uu3 + u4 * uu4 
           ii = ii + 1
           
        END DO
        END DO
        END DO
        IF( mpime == root ) THEN
           IF( lform ) THEN
             WRITE(iunit,10) ( sarray(ijk), ijk = 1, (ii-1) )
           ELSE
             WRITE(iunit) REAL( sarray( 1 : (ii-1) ), sgl )
           END IF
        END IF
        DEALLOCATE( sarray )

        CLOSE( iunit )

10      FORMAT( 5(G14.6,1X) )
      END SUBROUTINE inte_array_x
!-----------------------------------------------------------------------
      SUBROUTINE inte_array_y( var )
        CHARACTER(LEN = 16) :: filwri
        CHARACTER(LEN = 4) :: var
        INTEGER :: ii, i, j, k, ijk, imesh
        REAL*8 :: v1, v2, v3, v4, vv1, vv2, vv3, vv4, s, r
        REAL, ALLOCATABLE :: sarray(:)
        filwri='filter.' // var // '.' // lettera(irest)
        IF (lform) THEN
          OPEN( UNIT=iunit, FILE=filwri, STATUS='UNKNOWN' )
        ELSE 
          OPEN( UNIT=iunit, FORM='UNFORMATTED', FILE=filwri, STATUS='UNKNOWN' )
        END IF
        WRITE(6,fmt="('  inte_array_y: writing file ',A16)") filwri
        ALLOCATE( sarray( ntot ) )
        imesh = 0
        sarray = 0.0
        ii = 1
        DO k = 1, nz-1
        DO j = 1, ny-1
        DO i = 1, nx-1
           imesh = i + ( j-1 ) * nx + ( k-1 ) * nx * ny
           v1 = array(imesh)
           imesh = i + ( j-1 ) * nx + ( k-1+1 ) * nx * ny
           v2 = array(imesh)
           imesh = i+1 + ( j-1 ) * nx + ( k-1 ) * nx * ny
           v3 = array(imesh)
           imesh = i+1 + ( j-1 ) * nx + ( k-1+1 ) * nx * ny
           v4 = array(imesh)

           r = dz( k ) / ( dz( k ) + dz( k + 1 ) )
           s = dx( i ) / ( dx( i ) + dx( i + 1 ) )

	   vv1 = (1-r)*(1-s);
	   vv2 = r*(1-s);
	   vv3 = (1-r)*s;
	   vv4 = r*s;
	 
	   sarray( ii ) = v1 * vv1 + v2 * vv2 + v3 * vv3 + v4 * vv4 
           ii = ii + 1
           
        END DO
        END DO
        END DO
        IF( mpime == root ) THEN
           IF( lform ) THEN
             WRITE(iunit,10) ( sarray(ijk), ijk = 1, (ii-1) )
           ELSE
             WRITE(iunit) REAL( sarray( 1 : (ii-1) ), sgl )
           END IF
        END IF
        DEALLOCATE( sarray )

        CLOSE( iunit )

10      FORMAT( 5(G14.6,1X) )
      END SUBROUTINE inte_array_y
!-----------------------------------------------------------------------
      SUBROUTINE inte_array_z( var )
        CHARACTER(LEN = 16) :: filwri
        CHARACTER(LEN = 4) :: var
        INTEGER :: ii, i, j, k, ijk, imesh
        REAL*8 :: w1, w2, w3, w4, ww1, ww2, ww3, ww4, s, r
        REAL, ALLOCATABLE :: sarray(:)
        filwri='filter.' // var // '.' // lettera(irest)
        IF (lform) THEN
          OPEN( UNIT=iunit, FILE=filwri, STATUS='UNKNOWN' )
        ELSE 
          OPEN( UNIT=iunit, FORM='UNFORMATTED', FILE=filwri, STATUS='UNKNOWN' )
        END IF
        WRITE(6,fmt="('  inte_array_z: writing file ',A16)") filwri
        ALLOCATE( sarray( ntot ) )
        imesh = 0
        sarray = 0.0
        ii = 1
        DO k = 1, nz-1
        DO j = 1, ny-1
        DO i = 1, nx-1
           imesh = i + ( j-1 ) * nx + ( k-1 ) * nx * ny
           w1 = array(imesh)
           imesh = i+1 + ( j-1 ) * nx + ( k-1 ) * nx * ny
           w2 = array(imesh)
           imesh = i + ( j-1+1 ) * nx + ( k-1 ) * nx * ny
           w3 = array(imesh)
           imesh = i+1 + ( j-1+1 ) * nx + ( k-1 ) * nx * ny
           w4 = array(imesh)

           r = dx( i ) / ( dx( i ) + dx( i + 1 ) )
           s = dy( j ) / ( dy( j ) + dy( j + 1 ) )

	   ww1 = (1-r)*(1-s);
	   ww2 = r*(1-s);
	   ww3 = (1-r)*s;
	   ww4 = r*s;
	 
	   sarray( ii ) = w1 * ww1 + w2 * ww2 + w3 * ww3 + w4 * ww4 
           ii = ii + 1
           
        END DO
        END DO
        END DO
        IF( mpime == root ) THEN
           IF( lform ) THEN
             WRITE(iunit,10) ( sarray(ijk), ijk = 1, (ii-1) )
           ELSE
             WRITE(iunit) REAL( sarray( 1 : (ii-1) ), sgl )
           END IF
        END IF
        DEALLOCATE( sarray )

        CLOSE( iunit )

10      FORMAT( 5(G14.6,1X) )
      END SUBROUTINE inte_array_z
!-----------------------------------------------------------------------
      END SUBROUTINE filter_outp
!-----------------------------------------------------------------------
      END MODULE output_dump
!----------------------------------------------------------------------
