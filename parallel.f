!----------------------------------------------------------------------
      MODULE parallel
!----------------------------------------------------------------------
        IMPLICIT NONE
        SAVE
        INTEGER :: nproc, mpime, group, root
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
        SUBROUTINE parallel_startup
          CALL startup(nproc,mpime,root,group)
        END SUBROUTINE

        SUBROUTINE parallel_hangup
          CALL hangup
        END SUBROUTINE

!----------------------------------------------------------------------
      END MODULE
