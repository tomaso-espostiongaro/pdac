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
          CALL SYSTEM_CLOCK( iclk, count_rate = nclk )
          clocks_per_second = DBLE(nclk)
          cpclock = DBLE(iclk) / clocks_per_second * 1000.d0
! ... IBM/SGI ... 
!          cpclock = DBLE(mclock()) * 10.d0
          RETURN
        END FUNCTION

! ...   This function return the walltime elapsed since the first call
! ...   to the function itself
        FUNCTION elapsed_seconds()
          REAL*8 :: elapsed_seconds
          INTEGER :: iclk, nclk, xclk
          INTEGER, SAVE :: nover = 0
          INTEGER, SAVE :: iclk_old
          INTEGER, SAVE :: iclk_first = 0
          LOGICAL, SAVE :: first = .TRUE.
          CALL SYSTEM_CLOCK( iclk, count_rate = nclk, count_max = xclk )
          IF( first ) THEN
            iclk_first = iclk
            iclk_old   = iclk
            first = .FALSE.
          END IF
          clocks_per_second = DBLE(nclk)
          IF( iclk < iclk_old ) THEN
            nover = nover + 1
          END IF
          elapsed_seconds = DBLE( iclk + nover * xclk - iclk_first ) &
                          / clocks_per_second
          iclk_old = iclk
          RETURN
        END FUNCTION

!---------------------------------------------------------------------
      END MODULE environment
!---------------------------------------------------------------------
