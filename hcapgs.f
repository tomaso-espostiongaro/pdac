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
      ALLOCATE(gc_heat_capacity(ngas,ntot))
      ALLOCATE(solid_heat_capacity(nsolid,ntot))
      gc_heat_capacity = 0.0d0
      solid_heat_capacity = 0.0d0

      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE local_bounds_hcapgs
      USE dimensions
      USE grid, ONLY: ncint
!
      ALLOCATE(cp(ngas,ncint))
      ALLOCATE(ck(nsolid,ncint))
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE hcapg(cpc, tg)
!----------------------------------------------------------------------
! ... This routine computes the Temperature-dependent heat capacity 
! ... for each gas phase 
!
      USE dimensions
      USE gas_constants, ONLY: gmw, c_joule
      IMPLICIT NONE
!
      REAL*8, INTENT(IN) :: tg
      REAL*8, INTENT(OUT) :: cpc(:)
      REAL*8 :: t1,t2,t3
!
        cpc = 0.D0

        t1=tg
        t2=tg**2
        t3=tg**3

! ... o2
        cpc(1)=(2.811D1-3.68D-6*t1+1.746D-5*t2-1.065D-8*t3)/c_joule/gmw(1)
! ... n2  
        cpc(2)=(3.115D1-1.357D-2*t1+2.680D-5*t2-1.168D-8*t3)/c_joule/gmw(2)
! ... co2
        cpc(3)=(1.980D1+7.344D-2*t1-5.602D-5*t2+1.715D-8*t3)/c_joule/gmw(3)
! ... h2
        cpc(4)=(2.714D1+9.274D-3*t1-1.381D-5*t2+7.645D-9*t3)/c_joule/gmw(4)
! ... h2o
        cpc(5)=(3.224D1+1.924D-3*t1+1.055D-5*t2-3.596D-9*t3)/c_joule/gmw(5)
! ... Air
        cpc(6)=(3.0413D1-1.059D-2*t1+2.458D-5*t2-1.135D-8*t3)/c_joule/gmw(6)
! ... so2
        cpc(7)=(2.385D1+6.699D-2*t1-4.961D-5*t2+1.328D-8*t3)/c_joule/gmw(7)

      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE hcaps(cs, cps, ts)
!----------------------------------------------------------------------
! ... This routine computes the heat capacity for particles
!
      USE dimensions
      IMPLICIT NONE
!
      REAL*8, INTENT(IN) :: cps, ts
      REAL*8, INTENT(OUT) :: cs
!
! if cs depends on temperature ts, change this dependence
!
        cs = cps
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      END MODULE heat_capacity
