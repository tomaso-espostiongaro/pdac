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
      REAL*8, ALLOCATABLE, DIMENSION(:) :: tm, cm, mn, mnn, gpx, gpy, gpz
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: rlk, ygc
!
! ... derived averaged fields
!
      REAL*8, ALLOCATABLE, DIMENSION(:) :: rhom_gav, um_gav, vm_gav, wm_gav
      REAL*8, ALLOCATABLE, DIMENSION(:) :: tm_gav
!
! ... column arrays
!
      REAL*8, ALLOCATABLE :: p_axis(:,:),ug_axis(:,:),vg_axis(:,:),wg_axis(:,:),tg_axis(:,:)
      REAL*8, ALLOCATABLE :: xgc_axis(:,:,:)
      REAL*8, ALLOCATABLE :: eps_axis(:,:,:),us_axis(:,:,:),vs_axis(:,:,:),ws_axis(:,:,:), ts_axis(:,:,:)
      REAL*8, ALLOCATABLE :: rhom_axis(:,:), wm_axis(:,:), pd_axis(:,:)
! 
! ... time-average of column arrays
!
      REAL*8, ALLOCATABLE :: xgc_av(:,:)
      REAL*8, ALLOCATABLE :: ts_av(:,:)
      REAL*8, ALLOCATABLE :: p_av(:),ug_av(:),vg_av(:),wg_av(:),tg_av(:)
      REAL*8, ALLOCATABLE :: eps_av(:,:),us_av(:,:),vs_av(:,:),ws_av(:,:)
      REAL*8, ALLOCATABLE :: rhom_av(:), wm_av(:), pd_av(:)
!
! ... standard deviation of column arrays
!
      REAL*8, ALLOCATABLE :: p_sd(:),ug_sd(:),vg_sd(:),wg_sd(:),tg_sd(:)
!
      REAL*8, ALLOCATABLE :: xgc_sd(:,:)
      REAL*8, ALLOCATABLE :: eps_sd(:,:),us_sd(:,:),vs_sd(:,:),ws_sd(:,:)
      REAL*8, ALLOCATABLE :: ts_sd(:,:)
      REAL*8, ALLOCATABLE :: rhom_sd(:), wm_sd(:), pd_sd(:)
!
! ... plume variables
!
      REAL*8, ALLOCATABLE, DIMENSION(:) :: rhom_z, um_z, vm_z, wm_z
      REAL*8, ALLOCATABLE, DIMENSION(:) :: tm_z, surface
!
      REAL*8 :: time
!
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!---------------------------------------------------------------------
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
      INTEGER :: stat(20)

      ALLOCATE(epst(dime),STAT=stat(1))  ! Total particle fraction
      ALLOCATE(vf(dime),STAT=stat(2))    ! Void fraction
      ALLOCATE(lepst(dime),STAT=stat(3)) ! Log10 of the total part. frac.
      ALLOCATE(rhog(dime),STAT=stat(4))  ! Gas Density
      ALLOCATE(rgp(dime),STAT=stat(5))   ! Gas Density
      ALLOCATE(rhom(dime),STAT=stat(6))  ! Mixture Density
      ALLOCATE(tm(dime),STAT=stat(7))  ! Mixture Temperature
      ALLOCATE(um(dime),STAT=stat(8))  ! Mixture Velocity X
      ALLOCATE(vm(dime),STAT=stat(9))  ! Mixture Velocity Y
      ALLOCATE(wm(dime),STAT=stat(10))  ! Mixture Velocity Z
      ALLOCATE(mvm(dime),STAT=stat(11)) ! Mixture Velocity Modulus
      ALLOCATE(cm(dime),STAT=stat(12)) ! Mixture Sound Speed
      ALLOCATE(mn(dime),STAT=stat(13)) ! Mixture Mach Number
      ALLOCATE(mnn(dime),STAT=stat(14)) ! Mixture Normal Mach Number
      ALLOCATE(gpx(dime),STAT=stat(15)) ! Mixture Normal Mach Number
      ALLOCATE(gpy(dime),STAT=stat(16)) ! Mixture Normal Mach Number
      ALLOCATE(gpz(dime),STAT=stat(17)) ! Mixture Normal Mach Number
      ALLOCATE(pd(dime),STAT=stat(18))  ! Dynamic Pressure
      ALLOCATE(rlk(dime,nsolid),STAT=stat(19))  ! Solid Bulk density
      ALLOCATE(ygc(dime,ngas),STAT=stat(20))  ! Gas mass fractions
!
      IF (ANY(stat /= 0)) CALL error('postvar','problem allocating vars',1)
!
      epst  = 0.D0
      vf    = 1.D0
      lepst = 0.D0
      rhog  = 0.D0
      rgp   = 0.D0
      rhom  = 0.D0
      tm    = 0.D0
      um    = 0.D0
      vm    = 0.D0
      wm    = 0.D0
      mvm   = 0.D0
      cm    = 0.D0
      mn    = 0.D0
      mnn   = 0.D0
      gpx   = 0.D0
      gpy   = 0.D0
      gpz   = 0.D0
      pd    = 0.D0
      rlk   = 0.D0
      ygc   = 0.D0
!
      RETURN
      END SUBROUTINE allocate_derived_fields
!----------------------------------------------------------------------
      SUBROUTINE allocate_column_probes(dime,tdime)
! ... Time-averages and Standard deviations along vent axis
!
      USE dimensions, ONLY: ngas, nsolid
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: dime, tdime
!
! ... Allocation of the average arrays.
!
      ALLOCATE(p_axis(dime,tdime))
      ALLOCATE(ug_axis(dime,tdime))
      ALLOCATE(vg_axis(dime,tdime))
      ALLOCATE(wg_axis(dime,tdime))
      ALLOCATE(tg_axis(dime,tdime))
      ALLOCATE(xgc_axis(dime,ngas,tdime))
      ALLOCATE(eps_axis(dime,nsolid,tdime))
      ALLOCATE(us_axis(dime,nsolid,tdime))
      ALLOCATE(vs_axis(dime,nsolid,tdime))
      ALLOCATE(ws_axis(dime,nsolid,tdime))
      ALLOCATE(ts_axis(dime,nsolid,tdime))
      ALLOCATE(rhom_axis(dime,tdime))
      ALLOCATE(wm_axis(dime,tdime))
      ALLOCATE(pd_axis(dime,tdime)) 
!
! ... Allocation of the average arrays.
!
      ALLOCATE(p_av(dime))
      ALLOCATE(ug_av(dime))
      ALLOCATE(vg_av(dime))
      ALLOCATE(wg_av(dime))
      ALLOCATE(tg_av(dime))
      ALLOCATE(xgc_av(dime,ngas))
      ALLOCATE(eps_av(dime,nsolid))
      ALLOCATE(us_av(dime,nsolid))
      ALLOCATE(vs_av(dime,nsolid))
      ALLOCATE(ws_av(dime,nsolid))
      ALLOCATE(ts_av(dime,nsolid))
      ALLOCATE(rhom_av(dime))
      ALLOCATE(wm_av(dime))
      ALLOCATE(pd_av(dime)) 
! 
! ... Allocatation of the standard deviation arrays.
!
      ALLOCATE(p_sd(dime))
      ALLOCATE(ug_sd(dime))
      ALLOCATE(vg_sd(dime))
      ALLOCATE(wg_sd(dime))
      ALLOCATE(tg_sd(dime))
      ALLOCATE(xgc_sd(dime,ngas))
      ALLOCATE(eps_sd(dime,nsolid))
      ALLOCATE(us_sd(dime,nsolid))
      ALLOCATE(vs_sd(dime,nsolid))
      ALLOCATE(ws_sd(dime,nsolid))
      ALLOCATE(ts_sd(dime,nsolid))
      ALLOCATE(rhom_sd(dime))
      ALLOCATE(wm_sd(dime))
      ALLOCATE(pd_sd(dime)) 
!
! ... Initialization of the average arrays. 
!
      p_axis    = 0.D0
      ug_axis   = 0.D0
      vg_axis   = 0.D0
      wg_axis   = 0.D0
      tg_axis   = 0.D0
      xgc_axis  = 0.D0
      eps_axis  = 0.D0
      us_axis   = 0.D0
      vs_axis   = 0.D0
      ws_axis   = 0.D0
      ts_axis   = 0.D0
      rhom_axis = 0.D0
      wm_axis   = 0.D0
      pd_axis   = 0.D0
!
! ... Initialization of the average arrays. 
!
      p_av    = 0.D0
      ug_av   = 0.D0
      vg_av   = 0.D0
      wg_av   = 0.D0
      tg_av   = 0.D0
      xgc_av  = 0.D0
      eps_av  = 0.D0
      us_av   = 0.D0
      vs_av   = 0.D0
      ws_av   = 0.D0
      ts_av   = 0.D0
      rhom_av = 0.D0
      wm_av   = 0.D0
      pd_av   = 0.D0
!
! ... Initialization of the standard deviation arrays.
!
      p_sd    = 0.D0
      ug_sd   = 0.D0
      vg_sd   = 0.D0
      wg_sd   = 0.D0
      tg_sd   = 0.D0
      xgc_sd  = 0.D0
      eps_sd  = 0.D0
      us_sd   = 0.D0
      vs_sd   = 0.D0
      ws_sd   = 0.D0
      ts_sd   = 0.D0
      rhom_sd = 0.D0
      wm_sd   = 0.D0
      pd_sd   = 0.D0
!
      RETURN
      END SUBROUTINE allocate_column_probes
!----------------------------------------------------------------------
      SUBROUTINE allocate_plume_variables(dime)
      ! ... Horizontal averages on plume section
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: dime
!
      ALLOCATE(surface(dime)) ! jet plane surface
      ALLOCATE(rhom_z(dime))  ! Mixture Density along z
      ALLOCATE(tm_z(dime))    ! Mixture Density along z
      ALLOCATE(um_z(dime))    ! Mixture Velocity X along z
      ALLOCATE(vm_z(dime))    ! Mixture Velocity Y along z
      ALLOCATE(wm_z(dime))    ! Mixture Velocity Z along z
      surface = 0.0D0
      rhom_z = 0.0D0
      tm_z   = 0.0D0
      um_z   = 0.0D0
      vm_z   = 0.0D0
      wm_z   = 0.0D0
!
      RETURN
      END SUBROUTINE allocate_plume_variables
!----------------------------------------------------------------------
      SUBROUTINE allocate_derived_fields_av (dime)
      ! ... Time-averaged output fields
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: dime
!
      ALLOCATE(rhom_gav(dime))  ! Mixture Density
      ALLOCATE(tm_gav(dime))  ! Mixture Density
      ALLOCATE(um_gav(dime))    ! Mixture Velocity X
      ALLOCATE(vm_gav(dime))    ! Mixture Velocity Y
      ALLOCATE(wm_gav(dime))    ! Mixture Velocity Z
      rhom_gav = 0.0D0
      tm_gav   = 0.0D0
      um_gav   = 0.0D0
      vm_gav   = 0.0D0
      wm_gav   = 0.0D0
!
      RETURN
      END SUBROUTINE allocate_derived_fields_av
!-----------------------------------------------------------------------
      SUBROUTINE compute_derived_fields
!
      USE control_flags, ONLY: job_type
      USE derived_fields, ONLY: total_particle_fraction, void_fraction, &
                                log10_epstot, gas_density, gas_mass_fractions, &
                                gas_bulk_density, solid_bulk_density, &
                                mixture_density, mixture_velocity, &
                                velocity_module_2D, velocity_module_3D, &
                                dynamic_pressure, mixture_temperature,  &
                                mixture_sound_speed_1, mixture_sound_speed_2, &
                                mixture_sound_speed_3, &
                                mach_number,  gradient, normal_mach_number
      USE domain_mapping, ONLY: data_exchange
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
        CALL mixture_temperature(tm,ts,tg,rlk,rgp,ygc)
        !
        CALL mixture_velocity(um,ug,us,rlk,rgp,rhom)
        IF (job_type == '3D') CALL mixture_velocity(vm,vg,vs,rlk,rgp,rhom)
        CALL mixture_velocity(wm,wg,ws,rlk,rgp,rhom)
        IF (job_type == '3D') THEN
          CALL velocity_module_2D(mvm,um,vm)
          !CALL velocity_module_3D(mvm,um,vm,wm)
        ELSE IF (job_type == '2D') THEN
          CALL velocity_module_2D(mvm,um,wm)
        END IF
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
        CALL mixture_sound_speed_2(cm,xgc,rgp,rlk,rhom,rhog,epst,tg)
        CALL mach_number(mn,mvm,cm)
        !
        ! ... non-local operations require the data_exchange of ghost-cells!
        CALL data_exchange(p)
        CALL gradient(p,gpx,gpy,gpz)
        !
        CALL normal_mach_number(mnn,um,vm,wm,gpx,gpy,gpz,cm)
!
      RETURN
      END SUBROUTINE compute_derived_fields
!----------------------------------------------------------------------
      END MODULE postp_variables
!-----------------------------------------------------------------------
