!----------------------------------------------------------------------
      MODULE pressure_epsilon
!----------------------------------------------------------------------
      IMPLICIT NONE
      SAVE
!
      REAL*8, DIMENSION(:), ALLOCATABLE :: gas_pressure
      REAL*8, DIMENSION(:), ALLOCATABLE :: void_fraction
!
      REAL*8, DIMENSION(:), ALLOCATABLE :: ep
      REAL*8, DIMENSION(:), ALLOCATABLE :: p, pn
!
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE bounds_press_eps
      USE dimensions
      IMPLICIT NONE
!
      ALLOCATE(gas_pressure(ntot), void_fraction(ntot))
      gas_pressure  = 0.D0
      void_fraction = 0.D0
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE local_bounds_press_eps
      USE dimensions
      USE grid, ONLY: nijx_l
      IMPLICIT NONE
!
      ALLOCATE( ep( nijx_l ), p( nijx_l ), pn( nijx_l ))
      ep = 0.0d0
      p  = 0.0d0
      pn = 0.0d0
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      END MODULE
