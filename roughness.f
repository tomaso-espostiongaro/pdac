!----------------------------------------------
      MODULE roughness_module
!----------------------------------------------

        USE dimensions 
        IMPLICIT NONE
        SAVE
        TYPE roughness
          REAL*8 :: roucha
          REAL*8, POINTER :: r(:)
          INTEGER :: ir
        END TYPE
        TYPE (roughness) :: zrough

!----------------------------------------------
      CONTAINS
!----------------------------------------------

        SUBROUTINE allocate_roughness( r, nr )
          TYPE (roughness) :: r
          INTEGER, INTENT(IN) :: nr
            ALLOCATE( r%r(nr) )
          RETURN
        END SUBROUTINE

        SUBROUTINE deallocate_roughness( r )
          TYPE (roughness) :: r
            DEALLOCATE( r%r )
          RETURN
        END SUBROUTINE

!----------------------------------------------
      END MODULE
!----------------------------------------------
