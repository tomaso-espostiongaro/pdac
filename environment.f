!---------------------------------------------------------------------
      MODULE environment
!---------------------------------------------------------------------
        IMPLICIT NONE
        SAVE

        LOGICAL            :: tscra
        CHARACTER(LEN=256) :: cdir, inputdir, scradir, infile, outfile
        REAL*8, PRIVATE    :: max_seconds = 1.d+7
        REAL*8             :: clocks_per_second
        LOGICAL            :: timing = .TRUE.
!---------------------------------------------------------------------
        CONTAINS
!---------------------------------------------------------------------
! ...   This function return the user elapsed time in milliseconds
        FUNCTION cpclock()
          REAL*8 :: cpclock
!          INTEGER   :: mclock
          INTEGER   iclk, nclk
! ... T3E ...
          CALL SYSTEM_CLOCK(iclk,count_rate=nclk)
          clocks_per_second = DBLE(nclk)
          CALL SYSTEM_CLOCK(iclk)
          cpclock = DBLE(iclk) / clocks_per_second * 1000.d0
! ... IBM/SGI ... 
!          cpclock = DBLE(mclock()) * 10.d0
          RETURN
        END FUNCTION

!---------------------------------------------------------------------
        FUNCTION DIFFCLOCK(old_clock)
          REAL*8, INTENT(INOUT) :: old_clock
          REAL*8 :: diffclock, tmp, cclock
!          INTEGER mclock
          INTEGER   iclk
          CALL SYSTEM_CLOCK(iclk)
! ... IBM/SGI ... 
!          tmp = DBLE(mclock()) * 10.d0
! ... T3E ...
          tmp = DBLE(iclk) / clocks_per_second * 1000.d0
          diffclock = tmp - old_clock
          old_clock = tmp
          RETURN
        END FUNCTION
!---------------------------------------------------------------------
      END MODULE environment
!---------------------------------------------------------------------
