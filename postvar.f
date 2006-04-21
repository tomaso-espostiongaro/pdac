!-----------------------------------------------------------------------
      MODULE postp_variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
!
! ... main fields
!
      REAL*8, ALLOCATABLE, DIMENSION(:)   :: p, ug, vg, wg, tg
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: xgc
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: eps, us, vs, ws, ts, rlk
!
! ... derived fields
!
      REAL*8, ALLOCATABLE, DIMENSION(:) :: rm, rg, bd, m, um, vm, wm, mvm, c, mc
      REAL*8, ALLOCATABLE, DIMENSION(:) :: epstot, lepstot, pd, pdp, mvmp
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: sbd
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
                ws(dime,nsolid), ts(dime,nsolid), rlk(dime,nsolid))
      ALLOCATE(xgc(dime,ngas))

      RETURN
      END SUBROUTINE allocate_main_fields
!----------------------------------------------------------------------
      SUBROUTINE allocate_derived_fields(dime)
      USE dimensions, ONLY: nsolid
!
! ... These are the "more interesting" fields that can be derived
! ... from the primary OUTPUT fields produced by PDAC.
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: dime

      ALLOCATE(rm(dime))  ! Mixture Density
      ALLOCATE(rg(dime))  ! Gas Density
      ALLOCATE(bd(dime))  ! Bulk Density
      ALLOCATE(m(dime))   ! Gas Component Mass Fraction
      ALLOCATE(um(dime))  ! Mixture Velocity X
      ALLOCATE(vm(dime))  ! Mixture Velocity Y
      ALLOCATE(wm(dime))  ! Mixture Velocity Z
      ALLOCATE(mvm(dime)) ! Mixture Velocity Modulus
      ALLOCATE(mvmp(dime)) ! Mixture Velocity Modulus for probe
      ALLOCATE(c(dime))  ! Inverse of the Sound Speed
      ALLOCATE(mc(dime))  ! Mach Number
      ALLOCATE(epstot(dime))  ! Total particle fraction
      ALLOCATE(lepstot(dime))  ! Log10 of the total part. frac.
      ALLOCATE(pd(dime))  ! Dynamic Pressure
      ALLOCATE(pdp(dime))  ! Dynamic Pressure for probe
      ALLOCATE(sbd(dime,nsolid))  ! Solid Bulk density

      RETURN
      END SUBROUTINE allocate_derived_fields
!-----------------------------------------------------------------------
      SUBROUTINE compute_derived_fields
!
      USE control_flags, ONLY: job_type
      USE derived_fields
      USE io_files, ONLY: logunit
!
      IMPLICIT NONE
!
! ... Compute the derived fields ( ...also for maps)
!
        !
        ! ... Derived fields are computed as a function of
        ! ... primary fields and other derived fields
        !
!        rg = rhog(p,tg,xgc)
!        bd = rgp(eps,p,tg,xgc)
!        m  = mg(xgc)
!        c  = cm(bd,rg,rm,m,tg)
!        mc = mach(mvm,c)
        epstot = epst(eps)
        lepstot = leps(epstot)
        rm = rhom(eps,p,tg,xgc)
        um = velm(ug,us,eps,p,tg,xgc)
        IF (job_type == '3D') vm = velm(vg,vs,eps,p,tg,xgc)
        wm = velm(wg,ws,eps,p,tg,xgc)
        mvm = vel2(um,vm)
        pd = pdyn(rm,mvm)
        !
        ! ... Compute the dynamic pressure as the sum of PARTICLE dynamic pressures
        !
        !sbd  = rlk(eps)
        !mvm = vel2(us(:,1),vs(:,1))
        !pd = pdyn(sbd(:,1),mvm)
        !mvm = vel2(us(:,2),vs(:,2))
        !pd = pd + pdyn(sbd(:,2),mvm)
!
! ... Compute the derived fields for probing (... not for maps)
!
        !
        ! ... Derived fields are computed as a function of
        ! ... primary fields and other derived fields
        !

        IF (job_type == '3D') THEN
                mvmp = vel3(um,vm,wm)
        ELSE IF (job_type == '2D') THEN
                mvmp = vel2(um,wm)
        END IF
        pdp = pdyn(rm,mvmp)
!
      RETURN
      END SUBROUTINE compute_derived_fields
!----------------------------------------------------------------------
      END MODULE postp_variables
!-----------------------------------------------------------------------
