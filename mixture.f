!----------------------------------------------------------------------
      MODULE mixture_fields
!----------------------------------------------------------------------
      IMPLICIT NONE
!
      INTEGER :: compute_mixture
      REAL*8, DIMENSION(:), ALLOCATABLE :: rhom, um, vm, wm, tm
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: yd
!
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE allocate_mixture_fields
      USE dimensions
      USE domain_mapping, ONLY: ncint, ncdom
!
      IF (compute_mixture) THEN
        ALLOCATE(rhom(ncdom), um(ncdom), vm(ncdom), wm(ncdom), tm(ncdom))
        ALLOCATE(yd(ncdom,ngas+nsolid))
        rhom = 0.D0
        um = 0.D0
        vm = 0.D0
        wm = 0.D0
        tm = 0.D0
        yd = 0.D0
      END IF
!      
      RETURN
      END SUBROUTINE allocate_mixture_fields
!----------------------------------------------------------------------
      SUBROUTINE compute_mixture_fields
!
      USE control_flags, ONLY: job_type, JOB_TYPE_3D
      USE derived_fields, ONLY: mixture_density, mixture_temperature, &
                              & mixture_velocity, mixture_mass_fractions
      USE dimensions
      USE domain_mapping, ONLY: ncint, ncdom
      USE eos_gas, ONLY: ygc
      USE gas_solid_density, ONLY: rgp, rlk
      USE gas_solid_temperature, ONLY: tg, ts
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
!
      CALL mixture_density(rhom,rlk,rgp)
      CALL mixture_temperature(tm,ts,tg,rlk,rgp,ygc)
      CALL mixture_velocity(um,ug,us,rlk,rgp,rhom)
      IF (job_type == JOB_TYPE_3D) CALL mixture_velocity(vm,vg,vs,rlk,rgp,rhom)
      CALL mixture_velocity(wm,wg,ws,rlk,rgp,rhom)
      CALL mixture_mass_fractions(yd,rgp,rlk,rhom,ygc)
!      
      RETURN
      END SUBROUTINE compute_mixture_fields
!----------------------------------------------------------------------
      END MODULE mixture_fields
!----------------------------------------------------------------------
