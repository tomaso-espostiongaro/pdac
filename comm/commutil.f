      SUBROUTINE comm_getsiz( n )
        IMPLICIT NONE
        INTEGER n
        INCLUDE 'comm.h'
        n = nsiz
        RETURN
      END SUBROUTINE
