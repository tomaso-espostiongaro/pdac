!----------------------------------------------------------------------
      MODULE filter_outp
!----------------------------------------------------------------------
      USE dimensions
      USE kinds
      USE control_flags, ONLY: job_type
      USE io_files, ONLY: logunit
      USE io_serial, ONLY: read_array
      USE output_dump, ONLY: formatted_output
!
      IMPLICIT NONE
!
! ... parameters for filters
!
      INTEGER :: downsize_x, downsize_y, downsize_z
      INTEGER :: variable_n, field_n

      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE read_points( iunit, lform, ind, time, variables )

      USE dimensions, ONLY: nx, ny, nz, ntot

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iunit
      LOGICAL, INTENT(IN) :: lform
      INTEGER, INTENT(IN) :: ind(:)
      REAL*4, INTENT(OUT) :: time
      REAL*4, INTENT(OUT) :: variables(:,:)
      REAL*4, ALLOCATABLE   :: temp_array(:)
      INTEGER :: ijk, i,j,k
      INTEGER :: nvars, nump
      INTEGER :: nv, np

      ALLOCATE(temp_array(ntot))
      temp_array = 0.D0
      nvars = SIZE(variables,DIM=2)
      nump  = SIZE(variables,DIM=1)
      
      IF (lform) THEN
        READ(iunit,'(1x,///,1x,"@@@ TIME = ",g11.4)') time
      ELSE
        READ(iunit) time
      END IF
      
      DO nv = 1, nvars
        IF( lform ) READ(iunit,'(///)')
        CALL read_array( iunit, temp_array, lform )
        DO np = 1, nump
          variables(np,nv) = temp_array(ind(np))
        END DO
      END DO
      DEALLOCATE(temp_array)

      RETURN
      END SUBROUTINE read_points
!----------------------------------------------------------------------
      SUBROUTINE filter

      USE grid, ONLY: dx, dy, dz
      USE io_files, ONLY: outpunit
      USE io_serial, ONLY: first_out, last_out, incr_out
!
      IMPLICIT NONE
!
      INTEGER :: irest
      CHARACTER(LEN = 11) :: filnam
      CHARACTER(LEN = 4 ) :: lettera
      CHARACTER(LEN = 2 ) :: lettera2
      CHARACTER(LEN = 4 ) :: var
      LOGICAL :: lform
!
      INTEGER :: ig, is
      INTEGER :: iunit
      REAL :: time
      REAL   :: stime

      REAL, ALLOCATABLE :: array(:)
!
      ALLOCATE( array( ntot ) )

      DO irest = first_out, last_out, incr_out
!
        filnam='output.'//lettera(irest)
        WRITE(logunit,fmt="('  filter: reading file ',A11)") filnam
!
        lform = formatted_output
        iunit = 50
!
        IF (lform) THEN
          OPEN( UNIT=outpunit, FILE=filnam, STATUS='OLD')
          READ (outpunit,*) time
        ELSE 
          OPEN(UNIT=outpunit,FORM='UNFORMATTED',FILE=filnam,STATUS='OLD')
          READ (outpunit) stime
        END IF
!
        CALL read_array( outpunit, array, lform )  ! gas_pressure
        CALL crop_array( 'pgas' )  ! gas_pressure
!
        IF (job_type == '2D') THEN
          CALL read_array( outpunit, array, lform ) ! gas_velocity_r
          CALL inte_array_x( 'ugas' )  
          !
          CALL read_array( outpunit, array, lform ) ! gas_velocity_z
          CALL inte_array_z( 'wgas' ) 
        ELSE IF (job_type == '3D') THEN
          CALL read_array( outpunit, array, lform ) ! gas_velocity_x
          CALL inte_array_x( 'ugas' )  
          !
          CALL read_array( outpunit, array, lform ) ! gas_velocity_y
          CALL inte_array_y( 'vgas' )  
          !
          CALL read_array( outpunit, array, lform ) ! gas_velocity_z
          CALL inte_array_z( 'wgas' )  
        ELSE
          CALL error('outp_','Unknown job type',1)
        END IF
!
        CALL read_array( outpunit, array, lform )  ! gas_temperature
        CALL crop_array( 'tgas' )  
!
        DO ig=1,ngas
            var = 'xg'//lettera2( ig )
            CALL read_array( outpunit, array, lform )  ! gc_molar_fraction
            CALL crop_array( var )  
        END DO
!
        DO is = 1, nsolid
          CALL read_array( outpunit, array, lform )  ! solid_bulk_density
          var = 'ep'//lettera2( is )
          CALL crop_array( var )  
          IF (job_type == '2D') THEN
            CALL read_array( outpunit, array, lform )  ! solid_velocity_r
            var = 'us'//lettera2( is )
            CALL inte_array_x( var )  
            CALL read_array( outpunit, array, lform )  ! solid_velocity_z
            var = 'ws'//lettera2( is )
            CALL inte_array_z( var )  
          ELSE IF (job_type == '3D') THEN
            CALL read_array( outpunit, array, lform )  ! solid_velocity_x
            var = 'us'//lettera2( is )
            CALL inte_array_x( var )  
            CALL read_array( outpunit, array, lform )  ! solid_velocity_y
            var = 'vs'//lettera2( is )
            CALL inte_array_y( var )  
            CALL read_array( outpunit, array, lform )  ! solid_velocity_z
            var = 'ws'//lettera2( is )
            CALL inte_array_z( var )  
          END IF
          CALL read_array( outpunit, array, lform )  ! solid_temperature
          var = 'ts'//lettera2( is )
          CALL crop_array( var )  
        END DO
!
        CLOSE (outpunit)
!
      END DO
      DEALLOCATE( array )
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
          OPEN(UNIT=iunit,FORM='UNFORMATTED',FILE=filwri)
        END IF
        WRITE(logunit,fmt="('  crop_array: writing file ',A16)") filwri

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

        IF( lform ) THEN
          WRITE(iunit,10) ( sarray(ijk), ijk = 1, (ii-1) )
        ELSE
          WRITE(iunit) sarray( 1 : (ii-1) )
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
         OPEN(UNIT=iunit,FORM='UNFORMATTED',FILE=filwri,STATUS='UNKNOWN')
        END IF
        WRITE(logunit,fmt="('  inte_array_x: writing file ',A16)") filwri
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

        IF( lform ) THEN
          WRITE(iunit,10) ( sarray(ijk), ijk = 1, (ii-1) )
        ELSE
          WRITE(iunit) REAL( sarray( 1 : (ii-1) ), sgl )
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
         OPEN(UNIT=iunit,FORM='UNFORMATTED',FILE=filwri,STATUS='UNKNOWN')
        END IF
        WRITE(logunit,fmt="('  inte_array_y: writing file ',A16)") filwri
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

        IF( lform ) THEN
             WRITE(iunit,10) ( sarray(ijk), ijk = 1, (ii-1) )
           ELSE
             WRITE(iunit) REAL( sarray( 1 : (ii-1) ), sgl )
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
         OPEN(UNIT=iunit,FORM='UNFORMATTED',FILE=filwri,STATUS='UNKNOWN')
        END IF
        WRITE(logunit,fmt="('  inte_array_z: writing file ',A16)") filwri
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
        
        IF( lform ) THEN
          WRITE(iunit,10) ( sarray(ijk), ijk = 1, (ii-1) )
        ELSE
          WRITE(iunit) REAL( sarray( 1 : (ii-1) ), sgl )
        END IF
       
        DEALLOCATE( sarray )

        CLOSE( iunit )

10      FORMAT( 5(G14.6,1X) )
      END SUBROUTINE inte_array_z
!-----------------------------------------------------------------------
      END SUBROUTINE filter
!----------------------------------------------------------------------
      END MODULE filter_outp
!----------------------------------------------------------------------
