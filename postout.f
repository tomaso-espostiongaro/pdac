!-----------------------------------------------------------------------
      MODULE output_module
      SAVE
!-----------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------
      SUBROUTINE write_AVS_files
      USE iotk_module
      USE grid
!
! :::::::::::::::::::::   write post-processing files   :::::::::::::::::::::::
!
      REAL*8, ALLOCATABLE, DIMENSION(:) :: grx, gry, grz

! ... I/O units
!
      INTEGER :: iuni_fld = 21
      INTEGER :: iuni_grx = 22
      INTEGER :: iuni_gry = 23
      INTEGER :: iuni_grz = 24
!
      OPEN( UNIT=iuni_grx, FILE='pdac.grx', STATUS='UNKNOWN')
      OPEN( UNIT=iuni_gry, FILE='pdac.gry', STATUS='UNKNOWN')
      OPEN( UNIT=iuni_grz, FILE='pdac.grz', STATUS='UNKNOWN')
      OPEN( UNIT=iuni_fld, FILE='pdac.fld', STATUS='UNKNOWN')

      ALLOCATE( grx( nx - 1 ) )
      ALLOCATE( gry( ny - 1 ) )
      ALLOCATE( grz( nz - 1 ) )

      grx( 1 ) = origin_x
      WRITE( iuni_grx, * ) grx( 1 )
      DO i = 2, nx - 1
        grx( i ) = grx( i - 1 ) + dx( i )
        WRITE( iuni_grx, * ) grx( i )
      END DO

      gry( 1 ) = origin_y
      WRITE( iuni_gry, * ) gry( 1 )
      DO i = 2, ny - 1
        gry( i ) = gry( i - 1 ) + dy( i )
        WRITE( iuni_gry, * ) gry( i )
      END DO

      grz( 1 ) = origin_z
      WRITE( iuni_grz, * ) grz( 1 )
      DO i = 2, nz - 1
        grz( i ) = grz( i - 1 ) + dz( i )
        WRITE( iuni_grz, * ) grz( i )
      END DO

      WRITE( iuni_fld, fmt = "('# AVS field file')" )
      WRITE( iuni_fld, fmt = "('ndim=',I3)" ) 3
      WRITE( iuni_fld, fmt = "('dim1=',I3)" ) nx+1
      WRITE( iuni_fld, fmt = "('dim2=',I3)" ) ny+1
      WRITE( iuni_fld, fmt = "('dim3=',I3)" ) nz+1
      WRITE( iuni_fld, fmt = "('nspace=',I3)" ) 3
      WRITE( iuni_fld, fmt = "('veclen=',I3)" ) 3
      WRITE( iuni_fld, fmt = "('data=double')" )
      WRITE( iuni_fld, fmt = "('field=irregular')" )
      WRITE( iuni_fld, fmt = "('coord 1 file=pdac.grx, filetype=ascii')" )
      WRITE( iuni_fld, fmt = "('coord 2 file=pdac.gry, filetype=ascii')" )
      WRITE( iuni_fld, fmt = "('coord 3 file=pdac.grz, filetype=ascii')" )

      DEALLOCATE( grx )
      DEALLOCATE( gry )
      DEALLOCATE( grz )

      CLOSE( UNIT=iuni_fld )
      CLOSE( UNIT=iuni_grx )
      CLOSE( UNIT=iuni_gry )
      CLOSE( UNIT=iuni_grz )

      RETURN
      END SUBROUTINE write_AVS_files
!-----------------------------------------------------------------------
      END MODULE output_module
!----------------------------------------------------------------------
