!----------------------------------------------------------------------
      MODULE io_files
!----------------------------------------------------------------------
        IMPLICIT NONE
        SAVE
!
! ... I/O files
!
        INTEGER :: inputunit = 5
        INTEGER :: logunit   = 6
        INTEGER :: testunit  = 7
        INTEGER :: errorunit = 8
        INTEGER :: checkunit = 13
        INTEGER :: tempunit  = 9
        INTEGER :: blunit    = 14
        INTEGER :: iuni_nml  = 15
        INTEGER :: resunit   = 10
        INTEGER :: outpunit  = 12
        INTEGER :: postunit  = 11
        INTEGER :: ppunit    = 16

        INTEGER :: iuni_scalar = 21
        INTEGER :: iuni_u = 22
        INTEGER :: iuni_v = 23
        INTEGER :: iuni_w = 24

        INTEGER :: iunxml = 25

        CHARACTER(LEN=8)  :: inputfile = 'pdac.dat'
        CHARACTER(LEN=8)  :: logfile   = 'pdac.log'
        CHARACTER(LEN=8)  :: testfile  = 'pdac.tst'
        CHARACTER(LEN=8)  :: errorfile = 'pdac.err'
        CHARACTER(LEN=8)  :: checkfile = 'pdac.chm'
        CHARACTER(LEN=8)  :: blfile    = 'body.dat'
        CHARACTER(LEN=8)  :: nmlfile   = 'pdac.nml'
        CHARACTER(LEN=8)  :: resfile   = 'pdac.res'
        CHARACTER(LEN=10) :: xmlfile   = 'output.xml'
        CHARACTER(LEN=11) :: testnb
        CHARACTER(LEN=11) :: filnam
!----------------------------------------------------------------------
!      CONTAINS
!----------------------------------------------------------------------
      END MODULE io_files
!----------------------------------------------------------------------