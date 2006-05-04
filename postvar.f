!-----------------------------------------------------------------------
      MODULE postp_variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
!
! ... main fields
!
      REAL*8, ALLOCATABLE, DIMENSION(:)   :: p, ug, vg, wg, tg
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: xgc
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: eps, us, vs, ws, ts
!
! ... derived fields
!
      REAL*8, ALLOCATABLE, DIMENSION(:) :: epst, vf, lepst, rhog, rgp
      REAL*8, ALLOCATABLE, DIMENSION(:) :: rhom, um, vm, wm, pd, mvm
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: rlk, ygc
!
      REAL*8 :: time
!
      SAVE
!-----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE allocate_main_fields(dime)
      USE dimensions, ONLY: nsolid, ngas
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: dime

      ALLOCATE(p(dime), ug(dime), vg(dime), wg(dime), tg(dime))
      ALLOCATE(eps(dime,nsolid), us(dime,nsolid), vs(dime,nsolid), &
                ws(dime,nsolid), ts(dime,nsolid))
      ALLOCATE(xgc(dime,ngas))

      RETURN
      END SUBROUTINE allocate_main_fields
!----------------------------------------------------------------------
      SUBROUTINE allocate_derived_fields(dime)
      USE dimensions, ONLY: nsolid, ngas
!
! ... These are the "more interesting" fields that can be derived
! ... from the primary OUTPUT fields produced by PDAC.
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: dime

      ALLOCATE(epst(dime))  ! Total particle fraction
      ALLOCATE(vf(dime))    ! Void fraction
      ALLOCATE(lepst(dime)) ! Log10 of the total part. frac.
      ALLOCATE(rhog(dime))  ! Gas Density
      ALLOCATE(rgp(dime))   ! Gas Density
      ALLOCATE(rhom(dime))  ! Mixture Density
      ALLOCATE(um(dime))  ! Mixture Velocity X
      ALLOCATE(vm(dime))  ! Mixture Velocity Y
      ALLOCATE(wm(dime))  ! Mixture Velocity Z
      ALLOCATE(mvm(dime)) ! Mixture Velocity Modulus
      ALLOCATE(pd(dime))  ! Dynamic Pressure
      ALLOCATE(rlk(dime,nsolid))  ! Solid Bulk density
      ALLOCATE(ygc(dime,ngas))  ! Gas mass fractions

      RETURN
      END SUBROUTINE allocate_derived_fields
!-----------------------------------------------------------------------
      SUBROUTINE compute_derived_fields
!
      USE control_flags, ONLY: job_type
      USE derived_fields, ONLY: total_particle_fraction, void_fraction, &
                                log10_epstot, gas_density, gas_mass_fractions, &
                                gas_bulk_density, solid_bulk_density, &
                                mixture_density, mixture_velocity, &
                                velocity_module_2D, velocity_module_3D, &
                                dynamic_pressure
      USE io_files, ONLY: logunit
!
      IMPLICIT NONE
      REAL*8, ALLOCATABLE, DIMENSION(:) :: pd1, pd2
!
! ... Compute the derived fields
!
        !
        ! ... Derived fields are computed as a function of
        ! ... primary fields and other derived fields
        !
        CALL total_particle_fraction(epst,eps)
        CALL void_fraction(vf,epst)
        CALL log10_epstot(lepst,epst)
        CALL gas_density(rhog,p,tg,xgc)
        CALL gas_mass_fractions(ygc,xgc)
        CALL gas_bulk_density(rgp,vf,rhog)
        CALL solid_bulk_density(rlk,eps)
        CALL mixture_density(rhom,rlk,rgp)

        CALL mixture_velocity(um,ug,us,rlk,rgp,rhom)
        IF (job_type == '3D') CALL mixture_velocity(vm,vg,vs,rlk,rgp,rhom)
        CALL mixture_velocity(wm,wg,ws,rlk,rgp,rhom)
        CALL velocity_module_2D(mvm,um,vm)
        !CALL velocity_module_3D(mvm,um,vm,wm)
        CALL dynamic_pressure(pd,rhom,mvm)
        !
        ! ... Dynamic pressure can be computed as the sum of 
        ! ... PARTICLE dynamic pressures
        !
        !ALLOCATE(pd1(SIZE(pd)),pd2(SIZE(pd)))
        !CALL velocity_module_2D(mvm,us(:,1),vs(:,1))
        !CALL dynamic_pressure(pd1,rlk(:,1),mvm)
        !CALL velocity_module_2D(mvm,us(:,2),vs(:,2))
        !CALL dynamic_pressure(pd2,rlk(:,2),mvm)
        !pd = pd1 + pd2
        !DEALLOCATE(pd1,pd2)
!
      RETURN
      END SUBROUTINE compute_derived_fields
!----------------------------------------------------------------------
      END MODULE postp_variables
!-----------------------------------------------------------------------
