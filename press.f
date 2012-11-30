!----------------------------------------------------------------------
      MODULE pressure_epsilon
!----------------------------------------------------------------------
      IMPLICIT NONE
      SAVE
!
      REAL*8, DIMENSION(:), ALLOCATABLE :: ep, epn ! (void fraction)
      REAL*8, DIMENSION(:), ALLOCATABLE :: p, pn   ! (pressure)
      INTEGER :: pmodel
!
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE allocate_press_eps
      USE dimensions
      USE domain_mapping, ONLY: ncint, ncdom
      IMPLICIT NONE
!
      ALLOCATE( ep( ncdom ), p( ncdom ), pn( ncdom ), epn( ncdom ))
      ep  = 0.0d0
      p   = 0.0d0
      pn  = 0.0d0
      epn = 0.0d0
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      END MODULE
!----------------------------------------------------------------------
