!----------------------------------------------------------------------
      MODULE heat_capacity
!----------------------------------------------------------------------
      IMPLICIT NONE
      SAVE
!
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: gc_heat_capacity  
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: solid_heat_capacity
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: cp ! heat capacity of gas comp.
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: ck ! heat capacity of particles
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE bounds_hcapgs
      USE dimensions
!
      ALLOCATE(gc_heat_capacity(ngas,ndi*ndj))
      ALLOCATE(solid_heat_capacity(ncl,ndi*ndj))
      gc_heat_capacity = 0.0d0
      solid_heat_capacity = 0.0d0

      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE local_bounds_hcapgs
      USE dimensions
      USE grid, ONLY: nij_l
!
      ALLOCATE(cp(ngas,nij_l))
      ALLOCATE(ck(ncl,nij_l))
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE hcapg(cp, tg)
!----------------------------------------------------------------------
! ... This routine computes the Temperature-dependent heat capacity 
! ... for each gas phase 
!
      USE dimensions
      USE gas_constants, ONLY: gmw, c_joule
      IMPLICIT NONE
!
      REAL*8, INTENT(IN) :: tg
      REAL*8, INTENT(OUT) :: cp(:)
      REAL*8 :: c1
!
        c1=tg

!        cp(1)=(2.811D1-3.68D-6*c1+1.746D-5*c1**2-1.065D-8*c1**3)/c_joule/gmw(1)
!
!        cp(2)=(3.115D1-1.357D-2*c1+2.680D-5*c1**2-1.168D-8*c1**3)/c_joule/gmw(2)
!
!        cp(3)=(1.980D1+7.344D-2*c1-5.602D-5*c1**2+1.715D-8*c1**3)/c_joule/gmw(3)
!
!        cp(4)=(2.714D1+9.274D-3*c1-1.381D-5*c1**2+7.645D-9*c1**3)/c_joule/gmw(4)
!
        cp(5)=(3.224D1+1.924D-3*c1+1.055D-5*c1**2-3.596D-9*c1**3)/c_joule/gmw(5)
!
        cp(6)=(3.0413D1-1.059D-2*c1+2.458D-5*c1**2-1.135D-8*c1**3)/c_joule/gmw(6)
!
!        cp(7)=(2.385D1+6.699D-2*c1-4.961D-5*c1**2+1.328D-8*c1**3)/c_joule/gmw(7)
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE hcaps(ck, cps, tk)
!----------------------------------------------------------------------
! ... This routine computes the heat capacity for particles
!
      USE dimensions
      IMPLICIT NONE
!
      REAL*8, INTENT(IN) :: cps, tk
      REAL*8, INTENT(OUT) :: ck
!
! if ck depends on temperature tk, change this dependence
!
        ck = cps
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      END MODULE heat_capacity
