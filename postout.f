!-----------------------------------------------------------------------
      MODULE postp_output
      SAVE
!-----------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------
      SUBROUTINE write_AVS_files
      USE control_flags, ONLY: job_type
      USE dimensions
      USE iotk_module
      USE grid
      USE output_dump, ONLY: formatted_output
      IMPLICIT NONE
 
      CHARACTER(LEN=15) :: filetype
      INTEGER :: skip_m(3)
      INTEGER :: skip, skip_time, skip_phase, skip_block
      INTEGER :: lbl, ig, is
      INTEGER :: nfields, ndim, veclen, nlx, nly, nlz

!
! ... I/O units
!
      INTEGER :: iuni_scalar = 21
      INTEGER :: iuni_u = 22
      INTEGER :: iuni_v = 23
      INTEGER :: iuni_w = 24
!
      OPEN( UNIT=iuni_scalar, FILE='scalar.fld', STATUS='UNKNOWN')
      OPEN( UNIT=iuni_u, FILE='u.fld', STATUS='UNKNOWN')
      OPEN( UNIT=iuni_v, FILE='v.fld', STATUS='UNKNOWN')
      OPEN( UNIT=iuni_w, FILE='w.fld', STATUS='UNKNOWN')
!
! ... control parameters
!
      IF (job_type == '3D') THEN
        ndim = 3
      ELSE IF(job_type == '2D') THEN
        ndim = 2
      ELSE
        WRITE(*,*) 'Unknown job_type'
        STOP
      END IF
!
      IF (formatted_output) THEN
        filetype = 'ascii'
      ELSE
        filetype = 'unformatted'
      END IF

! ... Total number of scalar fields
!
      nfields = nphase*2 + ngas
!
      IF (formatted_output) THEN 
        ! ... number of lines for each block
        lbl = ntot / 10 + MIN(1,MOD(ntot,10)) + 4
      ELSE
        ! ... number of bytes for each block
        lbl = ntot * 4 + 8
      END IF 

! ... Number of lines in mesh file
!
        nlx = ( nx/5 + MIN(1,MOD(nx,5)) + 1 )
        nly = ( ny/5 + MIN(1,MOD(ny,5)) + 1 )
        nlz = ( nz/5 + MIN(1,MOD(nz,5)) + 1 )

      IF (formatted_output) THEN
        skip_time = 4 + 4
        skip_phase = (lbl * (ndim+2))
        skip_block = lbl
      ELSE
        skip_time = 12
        skip_phase = (lbl * (ndim+2))
        skip_block = lbl
      END IF

!***** S C A L A R S **************************************************
! ... Common Header
!
      WRITE( iuni_scalar, fmt = 1 )
      WRITE( iuni_scalar, fmt = 2 ) ndim
      IF (ndim == 3) THEN
        WRITE( iuni_scalar, fmt = 3 ) nx
        WRITE( iuni_scalar, fmt = 4 ) ny
        WRITE( iuni_scalar, fmt = 5 ) nz
      ELSE IF (ndim == 2) THEN
        WRITE( iuni_scalar, fmt = 3 ) nx
        WRITE( iuni_scalar, fmt = 4 ) nz
      END IF
      WRITE( iuni_scalar, fmt = 6 ) ndim
      WRITE( iuni_scalar, fmt = 7 )
      WRITE( iuni_scalar, fmt = 8 )
!
! ... specify variables and meshes
!
      veclen = nfields
      WRITE( iuni_scalar, fmt = 9 ) veclen
      skip_m(1) = 2
      WRITE( iuni_scalar, fmt = 10 ) skip_m(1)
      IF (ndim == 3) THEN
        skip_m(2) = 2 + 2*nlx
        WRITE( iuni_scalar, fmt = 11 ) skip_m(2)
        skip_m(3) = 2 + 2*nlx + 2*nly
        WRITE( iuni_scalar, fmt = 12 ) skip_m(3)
      ELSE IF (ndim == 2) THEN
        skip_m(3) = 2 + 2*nlx + 2*nly
        WRITE( iuni_scalar, fmt = 11 ) skip_m(3)
      END IF
!
! ... Scalar variables of GAS phase
!
      skip = skip_time
      WRITE( iuni_scalar, fmt = 13 ) 1, filetype, skip
      skip = skip + skip_block * (ndim+1)
      WRITE( iuni_scalar, fmt = 13 ) 2, filetype, skip

! ... Scalar variables of GAS COMPONENTS
!
      DO ig = 1, ngas
        skip = skip_time + skip_phase + (ig-1) * skip_block
        WRITE( iuni_scalar, fmt = 13 ) 2+ig, filetype, skip
      END DO

! ... Scalar variables of SOLID phases
!
      DO is = 1, nsolid
        skip = skip_time + skip_phase * is + ngas * skip_block
        WRITE( iuni_scalar, fmt = 13 ) 2+ngas+(is-1)*2+1, filetype, skip
        skip = skip + skip_block * (ndim+1)
        WRITE( iuni_scalar, fmt = 13 ) 2+ngas+(is)*2, filetype, skip
      END DO

      WRITE( iuni_scalar, fmt = 14 )
!
!***** V E L O C I T Y   X ********************************************
! ... Common Header
!
      WRITE( iuni_u, fmt = 1 )
      WRITE( iuni_u, fmt = 2 ) ndim
      IF (ndim == 3) THEN
        WRITE( iuni_u, fmt = 3 ) nx
        WRITE( iuni_u, fmt = 4 ) ny
        WRITE( iuni_u, fmt = 5 ) nz
      ELSE IF (ndim == 2) THEN
        WRITE( iuni_u, fmt = 3 ) nx
        WRITE( iuni_u, fmt = 4 ) nz
      END IF
      WRITE( iuni_u, fmt = 6 ) ndim
      WRITE( iuni_u, fmt = 7 )
      WRITE( iuni_u, fmt = 8 )
!
! ... specify variables and meshes
!
      veclen = nphase
      WRITE( iuni_u, fmt = 9 ) veclen
      skip_m(1) = 2 + nlx
      WRITE( iuni_u, fmt = 10 ) skip_m(1)
      IF (ndim == 3) THEN
        skip_m(2) = 2 + 2*nlx
        WRITE( iuni_u, fmt = 11 ) skip_m(2)
        skip_m(3) = 2 + 2*nlx + 2*nly
        WRITE( iuni_u, fmt = 12 ) skip_m(3)
      ELSE IF (ndim == 2) THEN
        skip_m(3) = 2 + 2*nlx + 2*nly
        WRITE( iuni_u, fmt = 11 ) skip_m(3)
      END IF
!
! ... GAS Velocity component X
!
      skip = skip_time + skip_block
      WRITE( iuni_u, fmt = 13 ) 1, filetype, skip
!
! ... SOLID Velocity component X
!
      DO is = 1, nsolid
        skip = skip_time + skip_phase * is + (ngas + 1) * skip_block
        WRITE( iuni_u, fmt = 13 ) is + 1, filetype, skip
      END DO
!
! ... EOF
!
      WRITE( iuni_u, fmt = 14 )
!
!***** V E L O C I T Y   Y ********************************************
! ... Common Header
!
      WRITE( iuni_v, fmt = 1 )
      WRITE( iuni_v, fmt = 2 ) ndim
      IF (ndim == 3) THEN
        WRITE( iuni_v, fmt = 3 ) nx
        WRITE( iuni_v, fmt = 4 ) ny
        WRITE( iuni_v, fmt = 5 ) nz
      ELSE IF (ndim == 2) THEN
        WRITE( iuni_v, fmt = 3 ) nx
        WRITE( iuni_v, fmt = 4 ) nz
      END IF
      WRITE( iuni_v, fmt = 6 ) ndim
      WRITE( iuni_v, fmt = 7 )
      WRITE( iuni_v, fmt = 8 )
!
! ... specify variables and meshes
!
      veclen = nphase
      WRITE( iuni_v, fmt = 9 ) veclen
      skip_m(1) = 2
      WRITE( iuni_v, fmt = 10 ) skip_m(1)
      IF (ndim == 3) THEN
        skip_m(2) = 2 + 2*nlx + nly
        WRITE( iuni_v, fmt = 11 ) skip_m(2)
        skip_m(3) = 2 + 2*nlx + 2*nly
        WRITE( iuni_v, fmt = 12 ) skip_m(3)
      ELSE IF (ndim == 2) THEN
        skip_m(3) = 2 + 2*nlx + 2*nly
        WRITE( iuni_v, fmt = 11 ) skip_m(3)
      END IF
!
! ... GAS Velocity component Y
!
      skip = skip_time + 2 * skip_block
      WRITE( iuni_v, fmt = 13 ) 1, filetype, skip
!
! ... SOLID Velocity component Y
!
      DO is = 1, nsolid
        skip = skip_time + skip_phase * is + (ngas + 2) * skip_block
        WRITE( iuni_v, fmt = 13 ) is + 1, filetype, skip
      END DO
!
! ... EOF
!
      WRITE( iuni_v, fmt = 14 )
!
!***** V E L O C I T Y   Z ********************************************
! ... Common Header
!
      WRITE( iuni_w, fmt = 1 )
      WRITE( iuni_w, fmt = 2 ) ndim
      IF (ndim == 3) THEN
        WRITE( iuni_w, fmt = 3 ) nx
        WRITE( iuni_w, fmt = 4 ) ny
        WRITE( iuni_w, fmt = 5 ) nz
      ELSE IF (ndim == 2) THEN
        WRITE( iuni_w, fmt = 3 ) nx
        WRITE( iuni_w, fmt = 4 ) nz
      END IF
      WRITE( iuni_w, fmt = 6 ) ndim
      WRITE( iuni_w, fmt = 7 )
      WRITE( iuni_w, fmt = 8 )
!
! ... specify variables and meshes
!
      veclen = nphase
      WRITE( iuni_w, fmt = 9 ) veclen
      skip_m(1) = 2 
      WRITE( iuni_w, fmt = 10 ) skip_m(1)
      IF (ndim == 3) THEN
        skip_m(2) = 2 + 2*nlx
        WRITE( iuni_w, fmt = 11 ) skip_m(2)
        skip_m(3) = 2 + 2*nlx + 2*nly + nlz
        WRITE( iuni_w, fmt = 12 ) skip_m(3)
      ELSE IF (ndim == 2) THEN
        skip_m(3) = 2 + 2*nlx + 2*nly + nlz
        WRITE( iuni_w, fmt = 11 ) skip_m(3)
      END IF
!
! ... GAS Velocity component Z
!
      skip = skip_time + ndim * skip_block
      WRITE( iuni_w, fmt = 13 ) 1, filetype, skip
!
! ... SOLID Velocity component Z
!
      DO is = 1, nsolid
        skip = skip_time + skip_phase * is + (ngas + ndim) * skip_block
        WRITE( iuni_w, fmt = 13 ) is + 1, filetype, skip
      END DO
!
! ... EOF
!
      WRITE( iuni_w, fmt = 14 )
!
!**********************************************************************
! ... Common Header
!
  1   FORMAT('# AVS field file')
  2   FORMAT('ndim=',I3)
  3   FORMAT('dim1=',I3)
  4   FORMAT('dim2=',I3)
  5   FORMAT('dim3=',I3)
  6   FORMAT('nspace=',I3)
  7   FORMAT('data=float')
  8   FORMAT('field=rectilinear')

! ... Specify variables and meshes
!
  9   FORMAT('veclen=',I3)
 10   FORMAT('coord 1 file=mesh.dat, filetype=ascii, skip= ',I6)
 11   FORMAT('coord 2 file=mesh.dat, filetype=ascii, skip= ',I6)
 12   FORMAT('coord 3 file=mesh.dat, filetype=ascii, skip= ',I6)
 13   FORMAT('variable',I3,' file=./output, filetype = ',A15,' skip= ',I12)

! ... Common EOF
!
 14   FORMAT('')

      CLOSE( UNIT=iuni_scalar )
      CLOSE( UNIT=iuni_u )
      CLOSE( UNIT=iuni_v )
      CLOSE( UNIT=iuni_w )

      RETURN
      END SUBROUTINE write_AVS_files
!-----------------------------------------------------------------------
      END MODULE postp_output
!----------------------------------------------------------------------
