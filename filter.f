!----------------------------------------------------------------------
      MODULE filter_outp
!----------------------------------------------------------------------
      USE dimensions
      USE kinds
      USE control_flags, ONLY: job_type
      USE io_files, ONLY: logunit
      USE output_dump, ONLY: formatted_output
!
      IMPLICIT NONE
!
! ... parameters for filters
!
      INTEGER :: act
      INTEGER :: first_out, last_out, incr_out
      INTEGER :: downsize_x, downsize_y, downsize_z
!
! ... parameters for sampling
!
      INTEGER :: number_of_probes
      LOGICAL :: assign_index
      INTEGER :: variable_n, field_n
      CHARACTER(LEN=80) :: probe_file
!
      REAL*8, ALLOCATABLE :: improfile(:,:,:)

      TYPE probe_point
           INTEGER :: nop
           INTEGER :: i
           INTEGER :: j
           INTEGER :: k
           REAL*8  :: x
           REAL*8  :: y
           REAL*8  :: z
      END TYPE probe_point
      !
      TYPE(probe_point), ALLOCATABLE :: probe(:)
!
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE sample
      USE control_flags, ONLY: job_type
      USE dimensions, ONLY: nsolid, ngas, nx, ny, nz
      USE grid, ONLY: x, y, z
      USE io_files, ONLY: outpunit, logunit, tempunit
!
      IMPLICIT NONE
      INTEGER :: ijk, i, j, k, ig, is, n
      CHARACTER(LEN = 11) :: filnam
      CHARACTER(LEN = 4 ) :: lettera
      CHARACTER(LEN = 20) :: probenam
      REAL*4, ALLOCATABLE :: vars(:,:)
      INTEGER, ALLOCATABLE :: ijk_probe(:), indx(:)
      INTEGER :: nop, tn, nfil, nvars, nv
      REAL*4   :: time
      LOGICAL :: lform
!
! ... Allocate probes.
! ... Read the file containing the indexes or coordinates of probe points
!
      ALLOCATE(probe(number_of_probes))
      ALLOCATE(ijk_probe(number_of_probes))
      ALLOCATE(indx(number_of_probes))
!
      OPEN(tempunit, FILE=probe_file, STATUS='OLD')
      DO nop = 1, number_of_probes
        probe(nop)%nop = nop
        IF (assign_index) THEN
                IF (job_type == '3D') THEN
                        READ(tempunit,*) probe(nop)%i, probe(nop)%j, probe(nop)%k
                        i = probe(nop)%i
                        j = probe(nop)%j
                        k = probe(nop)%k
                        ijk = i + (j-1)*nx + (k-1)*nx*ny
                        ijk_probe(nop) = ijk
                        !
                        probe(nop)%x = x(i)
                        probe(nop)%y = y(j)
                        probe(nop)%z = z(k)
                ELSE IF (job_type == '2D') THEN
                        READ(tempunit,*) probe(nop)%i,probe(nop)%k
                        i = probe(nop)%i
                        k = probe(nop)%k
                        ijk = i + (k-1)*nx
                        ijk_probe(nop) = ijk
                        !
                        probe(nop)%x = x(i)
                        probe(nop)%z = z(k)
                END IF
        ELSE
                IF (job_type == '3D') THEN
                        READ(tempunit,*) probe(nop)%x, probe(nop)%y, probe(nop)%z
                        probe(nop)%i = 1
                        probe(nop)%j = 1
                        probe(nop)%k = 1
                ELSE IF (job_type == '2D') THEN
                        READ(tempunit,*) probe(nop)%x, probe(nop)%z
                        probe(nop)%i = 1
                        probe(nop)%k = 1
                END IF
        END IF
      END DO
      CLOSE(tempunit)
!
! ... Sort the probes with progressively increasing index
! ... After the sorting, 'ijk_probe' contains the progressively
! ... increasing probe indexes, whereas 'indx' contains the
! ... probe indexes as read from the 'probe_file'
!
      CALL ikb07ad(ijk_probe(1), number_of_probes, indx(1))
!
      lform = formatted_output
!
! ... Define the total number of basic PDAC output variables
! ... and allocate arrays
!
      IF (job_type == '2D') THEN
              nvars = 4 * (nsolid + 1) + ngas
      ELSE IF (job_type == '3D') THEN
              nvars = 5 * (nsolid + 1) + ngas
      END IF
      !
      ALLOCATE(vars(nvars,number_of_probes))
      vars = 0.0
!
! ... Loop over time-steps
!
      DO tn = first_out, last_out, incr_out

        filnam='output.'//lettera(tn)
        WRITE(logunit,fmt="(/,'* Starting sampling ',I5,' * ')" ) tn

        ! ... Open PDAC output file
        !
        IF (lform) THEN
          OPEN(UNIT=outpunit, FILE=filnam)
        ELSE 
          OPEN(UNIT=outpunit,FORM='UNFORMATTED',FILE=filnam)
        END IF
        !
        ! ... Read sampling points in the progressive order
        !
        CALL read_points(outpunit,lform,ijk_probe,time,vars)
        !
        ! ... Loop over samplig points and write the corresponding files
        !
        DO nop = 1, number_of_probes
          n = indx(nop)
          i = probe(n)%i
          j = probe(n)%j
          k = probe(n)%k
          probenam ='S'//lettera(n)//'_'//lettera(i)//'_'//lettera(j)//'_'//lettera(k)
          OPEN(UNIT=tempunit, FILE=probenam, POSITION='APPEND')
            WRITE(tempunit,100) time, (vars(nv,nop), nv=1, nvars)
          CLOSE(tempunit)
        END DO

        CLOSE(outpunit)
      END DO
 100  FORMAT( F8.2, 100(G14.6E3,1X) )
!
      DEALLOCATE(vars)
      DEALLOCATE(probe)
      DEALLOCATE(ijk_probe) 
      DEALLOCATE(indx)

      RETURN
      END SUBROUTINE sample
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
      nvars = SIZE(variables,DIM=1)
      nump  = SIZE(variables,DIM=2)
      
      IF (lform) THEN
        READ(iunit,'(1x,///,1x,"@@@ TIME = ",g11.4)') time
      ELSE
        READ(iunit) time
      END IF
      
      DO nv = 1, nvars
        IF( lform ) READ(iunit,'(///)')
        CALL read_array( iunit, temp_array, lform )
        DO np = 1, nump
          variables(nv,np) = temp_array(ind(np))
        END DO
      END DO
      DEALLOCATE(temp_array)

      RETURN
      END SUBROUTINE read_points
!----------------------------------------------------------------------
      SUBROUTINE filter
      USE grid, ONLY: dx, dy, dz

      USE io_files, ONLY: outpunit
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
          OPEN(UNIT=outpunit,FORM='UNFORMATTED',FILE=filnam, STATUS='OLD')
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
          OPEN( UNIT=iunit, FORM='UNFORMATTED', FILE=filwri, STATUS='UNKNOWN' )
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
          OPEN( UNIT=iunit, FORM='UNFORMATTED', FILE=filwri, STATUS='UNKNOWN' )
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
          OPEN( UNIT=iunit, FORM='UNFORMATTED', FILE=filwri, STATUS='UNKNOWN' )
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
      SUBROUTINE read_array( iunit, sarray, lform )

      !  This subroutine reads a REAL array ( sarray ) of
      !  ntot elements from file iunit
      !  iunit  (input)  file to be read
      !  array  (output) data read from file and distributed to processors
      !  lform  (input)  format of the file 
      !                  .TRUE.  = formatted
      !                  .FALSE. = unformatted

      INTEGER, INTENT(IN) :: iunit
      LOGICAL, INTENT(IN) :: lform
      REAL :: sarray(:)

      INTEGER :: ijk, ierr

      IF( lform ) THEN
            READ(iunit,*) ( sarray(ijk), ijk = 1, ntot )
      ELSE
            READ(iunit) ( sarray(ijk), ijk = 1, ntot )
      END IF

      RETURN
      END SUBROUTINE read_array
!----------------------------------------------------------------------
      SUBROUTINE write_array( iunit, sarray, lform )

      !  This subroutine writes a REAL array ( sarray ) of
      !  ntot elements from file iunit
      !  iunit  (input)  file to be read
      !  array  (output) data read from file and distributed to processors
      !  lform  (input)  format of the file 
      !                  .TRUE.  = formatted
      !                  .FALSE. = unformatted

      INTEGER, INTENT(IN) :: iunit
      LOGICAL, INTENT(IN) :: lform
      REAL :: sarray(:)

      INTEGER :: ijk, ierr

      IF( lform ) THEN
        WRITE(iunit,100) ( sarray(ijk), ijk = 1, ntot )
      ELSE
        WRITE(iunit) ( sarray(ijk), ijk = 1, ntot )
      END IF

 100  FORMAT( 5(G14.6E3,1X) )
      RETURN
      END SUBROUTINE write_array
!----------------------------------------------------------------------
      SUBROUTINE read_implicit_profile
      USE grid, ONLY: x, y, z, xb, yb, zb
      USE io_files, ONLY: tempunit

      !  This subroutine reads the implicit profile
      !  and the mesh file computed by PDAC when a 3D
      !  volcano topography is imported from a DEM file

      INTEGER :: i,j,k
      ALLOCATE(improfile(nx,ny,nz))
!
! ... Read the georeferenced mesh
!
      OPEN(tempunit,FILE='mesh.dat',STATUS='OLD')
      READ(tempunit,*)
      READ(tempunit,*)
      READ(tempunit,*) (x(i), i=1,nx)
      READ(tempunit,*)
      READ(tempunit,*) (xb(i), i=1,nx)
      READ(tempunit,*)
      READ(tempunit,*) (y(j), j=1,ny)
      READ(tempunit,*)
      READ(tempunit,*) (yb(j), j=1,ny)
      READ(tempunit,*)
      READ(tempunit,*) (z(k), k=1,nz)
      READ(tempunit,*)
      READ(tempunit,*) (zb(k), k=1,nz)
      CLOSE(tempunit)
!      
! ... Read the Implicit Profile
!
      OPEN(tempunit,FILE='improfile.dat',STATUS='OLD')
      DO k=1,nz
        DO j=1,ny
          DO i=1,nx
            READ(tempunit,*) improfile(i,j,k)
          END DO
        END DO
      END DO
      CLOSE(tempunit)

 100  FORMAT(5(F20.6))
      RETURN
      END SUBROUTINE read_implicit_profile
!----------------------------------------------------------------------
      END MODULE filter_outp
!----------------------------------------------------------------------
