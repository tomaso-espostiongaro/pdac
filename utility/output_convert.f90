!----------------------------------------------------------------------
      PROGRAM convert_output
!----------------------------------------------------------------------
! Using a computer with a little Endian processor, this
! program reads output data as single 32 bit / 4 byte real data element from
! a sequential access file written on a computer with a big Endian
! processor, and converts this to its intended value using
! the subroutine native_4byte_real.
!
      IMPLICIT NONE
!
      INTEGER :: is, ig, ierr, n, nf
      INTEGER :: nsolid, ngas
      INTEGER :: nx, ny, nz
      INTEGER :: first_out, last_out, incr_out
      INTEGER :: stat1, stat2
      CHARACTER(LEN=2) :: job_type
      CHARACTER( LEN =  4 ) :: lettera
      CHARACTER(LEN=80) :: filename_in, filename_out
      REAL, ALLOCATABLE :: tmp1d( : ), tmp1dOut( : )
      REAL :: time, timeOut
!
      nx = 190
      ny = 170
      nz = 180
      nsolid = 2
      ngas = 2
      job_type = '3D'
      first_out = 10
      last_out = 10
      incr_out = 1
!
      ALLOCATE( tmp1d( nx*ny*nz ), STAT=stat1 )
      ALLOCATE( tmp1dOut(nx*ny*nz) , STAT=stat2)
      WRITE(*,*) 'Allocate status: ', stat1, stat2
      !
      DO nf = first_out, last_out, incr_out
      filename_in='output.'//lettera(nf)
      filename_out='output_little_endian.'//lettera(nf)
      OPEN( UNIT=10, FORM='UNFORMATTED', STATUS='OLD', FILE=filename_in)
      OPEN( UNIT=11, FORM='UNFORMATTED', STATUS='UNKNOWN', FILE=filename_out)
!
      WRITE(*,*) ' Converting output file: ', filename_in

      READ( UNIT=10 ) time
      WRITE(*,*) ' time =  ', time
      CALL native_4byte_real( time, timeOut )
      WRITE(11) timeOut

      WRITE(*,*) ' timeOut =  ', timeOut
      !
      READ(10) tmp1d(:)  ! P
STOP
      WRITE(*,*) ' converting p... '
!
      DO n = 1, SIZE(tmp1d)
        CALL native_4byte_real( tmp1d(n), tmp1dOut(n) )
      END DO
      WRITE(11) tmp1dOut
      !
      READ(10) tmp1d  ! ug
      WRITE(*,*) ' converting ug... '
!
      DO n = 1, SIZE(tmp1d)
        CALL native_4byte_real( tmp1d(n), tmp1dOut(n) )
      END DO
      WRITE(11) tmp1dOut
      !
      IF (job_type == '3D') THEN
        READ(10) tmp1d  ! vg
        WRITE(*,*) ' converting vg... '
        DO n = 1, SIZE(tmp1d)
          CALL native_4byte_real( tmp1d(n), tmp1dOut(n) )
        END DO
        WRITE(11) tmp1dOut
      END IF
      !
      READ(10) tmp1d  ! wg
      WRITE(*,*) ' converting wg... '
      DO n = 1, SIZE(tmp1d)
        CALL native_4byte_real( tmp1d(n), tmp1dOut(n) )
      END DO
      WRITE(11) tmp1dOut
      !
      READ(10) tmp1d  ! tg
      WRITE(*,*) ' converting tg... '
      DO n = 1, SIZE(tmp1d)
        CALL native_4byte_real( tmp1d(n), tmp1dOut(n) )
      END DO
      WRITE(11) tmp1dOut
      !
      DO ig = 1, ngas
        READ(10) tmp1d  ! ygc
        WRITE(*,*) ' converting ygc',ig,' ...'
        DO n = 1, SIZE(tmp1d)
          CALL native_4byte_real( tmp1d(n), tmp1dOut(n) )
        END DO
        WRITE(11) tmp1dOut
      END DO
      !
      DO is = 1, nsolid
      !
        READ(10) tmp1d  ! eps
        WRITE(*,*) ' converting eps', is, '... '
        DO n = 1, SIZE(tmp1d)
          CALL native_4byte_real( tmp1d(n), tmp1dOut(n) )
        END DO
        WRITE(11) tmp1dOut
      !
        READ(10) tmp1d  ! us
        WRITE(*,*) ' converting us', is, '... '
        DO n = 1, SIZE(tmp1d)
          CALL native_4byte_real( tmp1d(n), tmp1dOut(n) )
        END DO
        WRITE(11) tmp1dOut
      !
        IF (job_type == '3D') THEN
          READ(10) tmp1d  ! vs
          WRITE(*,*) ' converting vs', is, '... '
          DO n = 1, SIZE(tmp1d)
            CALL native_4byte_real( tmp1d(n), tmp1dOut(n) )
          END DO
          WRITE(11) tmp1dOut
        END IF
      !
        READ(10) tmp1d  ! ws
        WRITE(*,*) ' converting ws', is, '... '
        DO n = 1, SIZE(tmp1d)
          CALL native_4byte_real( tmp1d(n), tmp1dOut(n) )
        END DO
        WRITE(11) tmp1dOut
      !
        READ(10) tmp1d  ! ts
        WRITE(*,*) ' converting ts', is, '... '
        DO n = 1, SIZE(tmp1d)
          CALL native_4byte_real( tmp1d(n), tmp1dOut(n) )
        END DO
        WRITE(11) tmp1dOut
      !
      END DO
      !
      CLOSE (10)
      CLOSE (11)

      END DO
      DEALLOCATE( tmp1d, tmp1dOut )
!
      RETURN
      END PROGRAM convert_output
!--------------------------------------------------------------------------
!     SUBPROGRAM: native_4byte_real
!
!         AUTHOR: David Stepaniak, NCAR/CGD/CAS
! DATE INITIATED: 29 April 2003 
!  LAST MODIFIED: 29 April 2003
!
!       SYNOPSIS: Converts a 32 bit, 4 byte, REAL from big Endian to
!                 little Endian, or conversely from little Endian to big
!                 Endian.
!
!    DESCRIPTION: This subprogram allows one to convert a 32 bit, 4 byte,
!                 REAL data element that was generated with, say, a big
!                 Endian processor (e.g. Sun/sparc, SGI/R10000, etc.) to its
!                 equivalent little Endian representation for use on little
!                 Endian processors (e.g. PC/Pentium running Linux). The
!                 converse, little Endian to big Endian, also holds.
!                 This conversion is accomplished by writing the 32 bits of
!                 the REAL data element into a generic 32 bit INTEGER space
!                 with the TRANSFER intrinsic, reordering the 4 bytes with
!                 the MVBITS intrinsic, and writing the reordered bytes into
!                 a new 32 bit REAL data element, again with the TRANSFER
!                 intrinsic. The following schematic illustrates the
!                 reordering process
!
!
!                  --------    --------    --------    --------
!                 |    D   |  |    C   |  |    B   |  |    A   |  4 Bytes
!                  --------    --------    --------    --------
!                                                             |
!                                                              -> 1 bit
!                                       ||
!                                     MVBITS
!                                       ||
!                                       \/
!
!                  --------    --------    --------    --------
!                 |    A   |  |    B   |  |    C   |  |    D   |  4 Bytes
!                  --------    --------    --------    --------
!                         |           |           |           |
!                         24          16          8           0   <- bit
!                                                                 position
!
!          INPUT: realIn,  a single 32 bit, 4 byte REAL data element.
!         OUTPUT: realOut, a single 32 bit, 4 byte REAL data element, with
!                 reverse byte order to that of realIn.
!    RESTRICTION: It is assumed that the default REAL data element is
!                 32 bits / 4 bytes.
!
!
      SUBROUTINE native_4byte_real( realIn, realOut )
      IMPLICIT NONE
      REAL, INTENT(IN)                              :: realIn
                                                   ! a single 32 bit, 4 byte
                                                   ! REAL data element
      REAL, INTENT(OUT)                             :: realOut
                                                   ! a single 32 bit, 4 byte
                                                   ! REAL data element, with
                                                   ! reverse byte order to
                                                   ! that of realIn
!
!... Local variables (generic 32 bit INTEGER spaces):

      INTEGER                                       :: i_element
      INTEGER                                       :: i_element_br
!
!... Transfer 32 bit of realIn to generic 32 bit INTEGER space:
      i_element = TRANSFER( realIn, 0 )
!
!... Reverse order of 4 bytes in 32 bit INTEGER space:
      CALL MVBITS( i_element, 24, 8, i_element_br, 0  )
      CALL MVBITS( i_element, 16, 8, i_element_br, 8  )
      CALL MVBITS( i_element,  8, 8, i_element_br, 16 )
      CALL MVBITS( i_element,  0, 8, i_element_br, 24 )
!
!... Transfer reversed order bytes to 32 bit REAL space (realOut):
      realOut = TRANSFER( i_element_br, 0.0 )
!
      END SUBROUTINE
!--------------------------------------------------------------------------
      CHARACTER*4 FUNCTION lettera(k)
      IMPLICIT NONE
      CHARACTER ones,tens,hund,thou
!
      INTEGER :: k
!
      INTEGER :: iten, ione, ihund, ithou
!
      ithou=INT(k/1000)
      ihund=INT((k-(ithou*1000))/100)
      iten=INT((k-(ithou*1000)-(ihund*100))/10)
      ione=k-ithou*1000-ihund*100-iten*10
      ones=CHAR(ione+48)
      tens=CHAR(iten+48)
      hund=CHAR(ihunD+48)
      thou=CHAR(ithou+48)
      lettera=thou//hund//tens//ones
!
      RETURN
      END FUNCTION
!----------------------------------------------------------------------
