!----------------------------------------------------------------------
      MODULE specific_heat_module
!----------------------------------------------------------------------
      IMPLICIT NONE
      SAVE
!
! ... Specific heat of gas components
!
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: cp 
!
! ... Specific heat of particles
!
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: ck 
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE allocate_hcapgs
      USE dimensions
      USE domain_decomposition, ONLY: ncint
!
      ALLOCATE(cp(max_ngas,ncint))
      ALLOCATE(ck(nsolid,ncint))
!
      RETURN
      END SUBROUTINE allocate_hcapgs
!----------------------------------------------------------------------
      SUBROUTINE hcapg(cpc, tg)
!----------------------------------------------------------------------
! ... This routine computes the Temperature-dependent specific heat
! ... at constant pressure per kilograms of each gas phase 
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

      !  Carlo consistency position
      !
      t1 = MAX( 0.0d0, tg )
      t2 = t1*t1
      t3 = t1*t2
!
! ... molar capacity ( joule/(kelvin mole) )
!
! ... o2
        cpc(1)=(2.811D1-3.68D-6*t1+1.746D-5*t2-1.065D-8*t3)
! ... n2  
        cpc(2)=(3.115D1-1.357D-2*t1+2.680D-5*t2-1.168D-8*t3)
! ... co2
        cpc(3)=(1.980D1+7.344D-2*t1-5.602D-5*t2+1.715D-8*t3)
! ... h2
        cpc(4)=(2.714D1+9.274D-3*t1-1.381D-5*t2+7.645D-9*t3)
! ... h2o
        cpc(5)=(3.224D1+1.924D-3*t1+1.055D-5*t2-3.596D-9*t3)
! ... Air
        cpc(6)=(3.0413D1-1.059D-2*t1+2.458D-5*t2-1.135D-8*t3)
! ... so2
        cpc(7)=(2.385D1+6.699D-2*t1-4.961D-5*t2+1.328D-8*t3)
!
! ... specific heat ( joule/(kelvin kilogram) )
        cpc = cpc / gmw

      RETURN
      END SUBROUTINE hcapg
!----------------------------------------------------------------------
      SUBROUTINE hcaps(cs, cps, ts)
!----------------------------------------------------------------------
! ... This routine computes the specific heat (MKS) of particles
!
      USE dimensions
      IMPLICIT NONE
!
      REAL*8, INTENT(IN) :: cps, ts
      REAL*8, INTENT(OUT) :: cs
!
! ... specific heat ( joule/(kelvin kilogram) )
! ... (if cs depends on temperature (ts), change this dependence)
!
        cs = cps
!
      RETURN
      END SUBROUTINE hcaps
!----------------------------------------------------------------------
      END MODULE specific_heat_module
!----------------------------------------------------------------------
