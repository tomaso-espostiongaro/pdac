!----------------------------------------------------------------------
      SUBROUTINE filter( first_out, last_out, incr_out )

      USE output_dump, ONLY: outp, filter_outp

      IMPLICIT NONE
!
      INTEGER :: first_out, last_out, incr_out 
      INTEGER :: irest
!
!---------------------------------------
      DO irest = first_out, last_out, incr_out
!---------------------------------------
!
       WRITE(6,fmt="(/,'* Starting filtering ',I5,' * ')" ) irest

       CALL filter_outp ( irest )
!
       CALL myflush( 6 )
!
!------------------------------------------
      END DO
!------------------------------------------
! 
      RETURN
      END SUBROUTINE filter
!----------------------------------------------------------------------
