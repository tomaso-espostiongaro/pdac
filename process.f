!----------------------------------------------------------------------
      MODULE process_outp
!----------------------------------------------------------------------
      USE dimensions
      USE kinds
      USE control_flags, ONLY: job_type
      USE io_files, ONLY: logunit
!
      IMPLICIT NONE
      SAVE
!
! ... parameters for filters
!
      INTEGER :: act
      INTEGER :: first_out, last_out, incr_out
      INTEGER :: downsize_x, downsize_y, downsize_z
!
! ... parameters for sampling
!
      INTEGER :: number_of_points
      LOGICAL :: assign_index
      INTEGER :: index_i, index_j, index_k
      REAL*8  :: coordinate_x, coordinate_y, coordinate_z
      INTEGER :: variable_n, field_n
!
      LOGICAL :: interpolate = .TRUE.
      LOGICAL :: formatted_output
!
! ... main fields
!
      REAL, ALLOCATABLE, DIMENSION(:)   :: p, ug, vg, wg, tg
      REAL, ALLOCATABLE, DIMENSION(:,:) :: xgc
      REAL, ALLOCATABLE, DIMENSION(:,:) :: eps, us, vs, ws, ts
!
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE allocate_main_fields

      IMPLICIT NONE

      ALLOCATE(p(ntot), ug(ntot), vg(ntot), wg(ntot), tg(ntot))
      ALLOCATE(eps(ntot,nsolid), us(ntot,nsolid), vs(ntot,nsolid), &
                ws(ntot,nsolid), ts(ntot,nsolid))
      ALLOCATE(xgc(ntot,ngas))

      RETURN
      END SUBROUTINE allocate_main_fields
!----------------------------------------------------------------------
      SUBROUTINE read_output( tn )

      USE io_files, ONLY: outpunit
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: tn
      CHARACTER(LEN = 11) :: filnam
      CHARACTER(LEN = 4 ) :: lettera
      CHARACTER(LEN = 2 ) :: lettera2
      LOGICAL :: lform
!
      INTEGER :: ig, is, i
      INTEGER :: iunit
      REAL   :: time
      REAL   :: stime
!
      filnam='output.'//lettera(tn)

      lform = formatted_output

      iunit = 50

      IF (lform) THEN
        OPEN( UNIT=outpunit, FILE=filnam, STATUS='OLD')
        READ(outpunit,'(1x,///,1x,"@@@ TIME = ",g11.4)') time
      ELSE 
        OPEN(UNIT=outpunit,FORM='UNFORMATTED',FILE=filnam, STATUS='OLD')
        READ (outpunit) stime
      END IF

      IF( lform ) READ(outpunit,'(///)')
      CALL read_array( outpunit, p, lform )  ! gas_pressure

      IF (job_type == '2D') THEN

        IF( lform ) READ(outpunit,'(///)')
        CALL read_array( outpunit, ug, lform ) ! gas_velocity_r
        IF( lform ) READ(outpunit,'(///)')
        CALL read_array( outpunit, wg, lform ) ! gas_velocity_z

      ELSE IF (job_type == '3D') THEN

        IF( lform ) READ(outpunit,'(///)')
        CALL read_array( outpunit, ug, lform ) ! gas_velocity_x
        IF( lform ) READ(outpunit,'(///)')
        CALL read_array( outpunit, vg, lform ) ! gas_velocity_y
        IF( lform ) READ(outpunit,'(///)')
        CALL read_array( outpunit, wg, lform ) ! gas_velocity_z

      ELSE
        CALL error('outp_','Unknown job type',1)
      END IF

      IF( lform ) READ(outpunit,'(///)')
      CALL read_array( outpunit, tg, lform )  ! gas_temperature

      DO ig=1,ngas
        IF( lform ) READ(outpunit,'(///)')
        CALL read_array( outpunit, xgc(:,ig), lform )  ! gc_molar_fraction
      END DO

      DO is = 1, nsolid

        IF( lform ) READ(outpunit,'(///)')
        CALL read_array( outpunit, eps(:,is), lform )  ! solid_bulk_density

        IF (job_type == '2D') THEN

        IF( lform ) READ(outpunit,'(///)')
          CALL read_array( outpunit, us(:,is), lform )  ! solid_velocity_r
        IF( lform ) READ(outpunit,'(///)')
          CALL read_array( outpunit, ws(:,is), lform )  ! solid_velocity_z

        ELSE IF (job_type == '3D') THEN

        IF( lform ) READ(outpunit,'(///)')
          CALL read_array( outpunit, us(:,is), lform )  ! solid_velocity_x
        IF( lform ) READ(outpunit,'(///)')
          CALL read_array( outpunit, vs(:,is), lform )  ! solid_velocity_y
        IF( lform ) READ(outpunit,'(///)')
          CALL read_array( outpunit, ws(:,is), lform )  ! solid_velocity_z

        END IF

        IF( lform ) READ(outpunit,'(///)')
        CALL read_array( outpunit, ts(:,is), lform )  ! solid_temperature

      END DO

      CLOSE (outpunit)
!
      RETURN
      END SUBROUTINE read_output
!-----------------------------------------------------------------------
      SUBROUTINE process
      USE derived_fields
      USE io_files, ONLY: tempunit
      USE grid, ONLY: z

      IMPLICIT NONE
!
      INTEGER :: tn, ijk, i, k
      LOGICAL :: lform
      REAL, ALLOCATABLE, DIMENSION(:) :: rm, rg, bd, m, um, wm, mvm, c, mc 
      REAL, ALLOCATABLE, DIMENSION(:) :: epstot, lepstot
      CHARACTER(LEN = 14) :: filnam
      CHARACTER(LEN = 4 ) :: lettera
!
      lform = formatted_output

      ALLOCATE(rm(ntot))  ! Mixture Density
      ALLOCATE(rg(ntot))  ! Gas Density
      ALLOCATE(bd(ntot))  ! Bulk Density
      ALLOCATE(m(ntot))   ! Gas Component Mass Fraction
      ALLOCATE(um(ntot))  ! Mixture Velocity X
      ALLOCATE(wm(ntot))  ! Mixture Velocity Z
      ALLOCATE(mvm(ntot)) ! Mixture Velocity Modulus
      ALLOCATE(c(ntot))  ! Inverse of the Sound Speed
      ALLOCATE(mc(ntot))  ! Mach Number
      ALLOCATE(epstot(ntot))  ! Total particle fraction
      ALLOCATE(lepstot(ntot))  ! Log10 of the total part. frac.

      CALL allocate_main_fields
!
      DO tn = first_out, last_out, incr_out

        filnam = 'log10epst.'//lettera(tn)
        OPEN(tempunit,FILE=filnam, STATUS='NEW', FORM='UNFORMATTED')
        WRITE(logunit,fmt="(/,'* Starting post-processing ',I5,' * ')" ) tn
        CALL read_output ( tn )

!        rm = rhom(eps,p,tg,xgc)
!        rg = rhog(p,tg,xgc)
!        bd = rgp(eps,p,tg,xgc)
!        m  = mg(xgc)
!
!        um = velm(ug,us,eps,p,tg,xgc)
!        wm = velm(wg,ws,eps,p,tg,xgc)
!        mvm = vel(um,wm)
!        c  = cm(bd,rg,rm,m,tg)
        epstot = epst(eps)
        lepstot = leps(epstot)

!        mc = mach(mvm,c)

        CALL write_array( tempunit, lepstot, lform )
        CLOSE(tempunit)

        filnam = 'tg.'//lettera(tn)
        OPEN(tempunit,FILE=filnam, STATUS='NEW', FORM='UNFORMATTED')
        WRITE(logunit,fmt="(/,'* Starting post-processing ',I5,' * ')" ) tn

        CALL write_array( tempunit, tg, lform )
        CLOSE(tempunit)

      END DO
!
      DEALLOCATE(rm)
      DEALLOCATE(rg)
      DEALLOCATE(bd)
      DEALLOCATE(m)
      DEALLOCATE(um)
      DEALLOCATE(wm)
      DEALLOCATE(mvm)
      DEALLOCATE(c)
      DEALLOCATE(mc)
      DEALLOCATE(epstot)
      DEALLOCATE(lepstot)
! 
      RETURN
      END SUBROUTINE process
!-----------------------------------------------------------------------
      SUBROUTINE filter

      IMPLICIT NONE
!
      INTEGER :: irest
!
      DO irest = first_out, last_out, incr_out

        WRITE(logunit,fmt="(/,'* Starting filtering ',I5,' * ')" ) irest
        CALL filter_outp ( irest )

      END DO
! 
      RETURN
      END SUBROUTINE filter
!----------------------------------------------------------------------
      SUBROUTINE filter_outp( irest )
      USE grid, ONLY: dx, dy, dz

      USE io_files, ONLY: outpunit
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
      REAL :: time
      REAL   :: stime

      REAL, ALLOCATABLE :: array(:)
!
      filnam='output.'//lettera(irest)
      WRITE(logunit,fmt="('  filter: reading file ',A11)") filnam

      lform = formatted_output

      iunit = 50

      IF (lform) THEN
        OPEN( UNIT=outpunit, FILE=filnam, STATUS='OLD')
        READ (outpunit,*) time
      ELSE 
        OPEN(UNIT=outpunit,FORM='UNFORMATTED',FILE=filnam, STATUS='OLD')
        READ (outpunit) stime
      END IF

      ALLOCATE( array( ntot ) )

      WRITE(logunit,fmt="('  filtering gas pressure ')")
!
      CALL read_array( outpunit, array, lform )  ! gas_pressure

      CALL crop_array( 'pgas' )  ! gas_pressure

      WRITE(logunit,fmt="('  filtering reading gas velocities ')")

      IF (job_type == '2D') THEN

        CALL read_array( outpunit, array, lform ) ! gas_velocity_r
        CALL inte_array_x( 'ugas' )  

        CALL read_array( outpunit, array, lform ) ! gas_velocity_z
        CALL inte_array_z( 'wgas' ) 

      ELSE IF (job_type == '3D') THEN

        CALL read_array( outpunit, array, lform ) ! gas_velocity_x
        CALL inte_array_x( 'ugas' )  

        CALL read_array( outpunit, array, lform ) ! gas_velocity_y
        CALL inte_array_y( 'vgas' )  

        CALL read_array( outpunit, array, lform ) ! gas_velocity_z
        CALL inte_array_z( 'wgas' )  

      ELSE
        CALL error('outp_','Unknown job type',1)
      END IF

      WRITE(logunit,fmt="('  filtering gas temperature ')")

      CALL read_array( outpunit, array, lform )  ! gas_temperature
      CALL crop_array( 'tgas' )  

      WRITE(logunit,fmt="('  filtering molarfraction ')")
!
      DO ig=1,ngas
          var = 'xg'//lettera2( ig )
          CALL read_array( outpunit, array, lform )  ! gc_molar_fraction
          CALL crop_array( var )  
      END DO

      WRITE( 6, fmt="('  filtering solid density, velocities and  &
                        & temperature')")
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

      DEALLOCATE( array )

      CLOSE (outpunit)
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
      END SUBROUTINE filter_outp
!----------------------------------------------------------------------
      SUBROUTINE read_array( iunit, sarray, lform )

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
        READ(iunit,*) ( sarray(ijk), ijk = 1, ntot )
      ELSE
        READ(iunit) ( sarray(ijk), ijk = 1, ntot )
      END IF

      RETURN
      END SUBROUTINE read_array
!-----------------------------------------------------------------------
      SUBROUTINE write_array( iunit, sarray, lform )

      !  This subroutine writes a REAL array ( sarray ) of
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
        WRITE(iunit,100) ( sarray(ijk), ijk = 1, ntot )
      ELSE
        WRITE(iunit) ( sarray(ijk), ijk = 1, ntot )
      END IF

 100  FORMAT( 5(G14.6E3,1X) )
      RETURN
      END SUBROUTINE write_array
!-----------------------------------------------------------------------
      END MODULE process_outp
!----------------------------------------------------------------------
