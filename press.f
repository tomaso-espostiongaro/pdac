!----------------------------------------------------------------------
      MODULE pressure_epsilon
!----------------------------------------------------------------------
      IMPLICIT NONE
      SAVE
!
      REAL*8, DIMENSION(:), ALLOCATABLE :: ep
      REAL*8, DIMENSION(:), ALLOCATABLE :: p, pn
!
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE allocate_press_eps
      USE dimensions
      USE domain_decomposition, ONLY: ncdom
      IMPLICIT NONE
!
      ALLOCATE( ep( ncdom ), p( ncdom ), pn( ncdom ))
      ep = 0.0d0
      p  = 0.0d0
      pn = 0.0d0
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      END MODULE
!----------------------------------------------------------------------
