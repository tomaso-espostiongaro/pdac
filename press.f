!----------------------------------------------------------------------
      MODULE pressure_epsilon
!----------------------------------------------------------------------
      IMPLICIT NONE
      SAVE
!
      REAL*8, DIMENSION(:), ALLOCATABLE :: p_g, pn_g
      REAL*8, DIMENSION(:), ALLOCATABLE :: ep_g
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
      ALLOCATE(p_g(ndi*ndj), pn_g(ndi*ndj),ep_g(ndi*ndj))
      p_g =0.0d0 
      pn_g = 0.0d0
      ep_g = 0.0d0

      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE local_bounds_press_eps
      USE dimensions
      USE grid, ONLY: nijx_l
      IMPLICIT NONE
!
      ALLOCATE( ep( nijx_l ), p( nijx_l ), pn ( nijx_l ))
      ep = 0.0d0
      p  = 0.0d0
      pn = 0.0d0
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      END MODULE
