!----------------------------------------------------------------------
      MODULE convective_fluxes_w
!
! ... This module computes convective fluxes of _____ z-momentum ______
! ... by using an upwind scheme.
! ... Standard compass notation is adopted to locate cell neighbours:
! ... %e=East, %w=West, %n=North, %s=South, %t=Top, %b=Bottom, etc.
! ... The computational stencil is defined in the `set_indexes' module
!
!----------------------------------------------------------------------
      IMPLICIT NONE
!
      REAL*8, PRIVATE :: cs                      ! convective stream   !
      REAL*8, PRIVATE :: cn                      ! Courant number      !
      REAL*8, PRIVATE :: upwnd                   ! upwinded variable   !
      REAL*8, PRIVATE :: lim                     ! limiter             !
      REAL*8, PRIVATE :: erre                    ! gradients ratio     !
      REAL*8, PRIVATE :: incr                    ! increment           !
!
      INTERFACE flw
        MODULE PROCEDURE flw_2d, flw_3d
      END INTERFACE
      INTERFACE muscl_flw
        MODULE PROCEDURE muscl_flw_2d, muscl_flw_3d
      END INTERFACE
      INTERFACE ctu1_flw
        MODULE PROCEDURE ctu1_flw_2d, ctu1_flw_3d
      END INTERFACE
      INTERFACE ctu2_flw
        MODULE PROCEDURE ctu2_flw_2d, ctu2_flw_3d
      END INTERFACE
      INTERFACE ctu3_flw
        MODULE PROCEDURE ctu3_flw_2d, ctu3_flw_3d
      END INTERFACE
     
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE flw_3d(fe, fn, ft, fw, fs, fb, dens, u, v, w, k)
!
! ... Compute the convective fluxes on East, North, and Top sides of the cell
! ... for the momentum density along z.
!
      USE grid, ONLY: dz, flag
      USE set_indexes, ONLY: imjk, ijmk, ijkm
      USE set_indexes, ONLY: stencil
      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: fe, fn, ft, fw, fs, fb
      TYPE(stencil), INTENT(IN) :: dens, u, v, w
      INTEGER, INTENT(IN) :: k

      REAL*8 :: dzp, indzp
!
      dzp=dz(k)+dz(k+1)
      indzp=1.D0/dzp
!
! ... on West volume bondary
!
      IF ( .NOT.BTEST(flag(imjk),0) ) THEN
        cs = (dz(k+1)*u%w + dz(k)*u%wt) * indzp
        fw = 0.5D0*dens%w*w%w*(cs+ABS(cs)) + 0.5D0*dens%c*w%c*(cs-ABS(cs))
      END IF
!
! ... on South volume bondary
!
      IF ( .NOT.BTEST(flag(ijmk),0) ) THEN
        cs = (dz(k+1)*v%s + dz(k)*v%st) * indzp
        fs = 0.5D0*dens%s*w%s*(cs+ABS(cs)) + 0.5D0*dens%c*w%c*(cs-ABS(cs))
      END IF
!
! ... on Bottom volume bondary
!
      IF ( .NOT.BTEST(flag(ijkm),0) ) THEN
        cs=0.5D0*(w%b+w%c)
        fb = 0.5D0*dens%b*w%b*(cs+ABS(cs)) + 0.5D0*dens%c*w%c*(cs-ABS(cs))
      END IF
!
! ... on East volume boundary
!
      cs = (dz(k+1)*u%c+dz(k)*u%t)*indzp
      fe = 0.5D0*dens%c*w%c*(cs+ABS(cs)) + 0.5D0*dens%e*w%e*(cs-ABS(cs))
!
! ... on North volume boundary
!
      cs = (dz(k+1)*v%c+dz(k)*v%t)*indzp
      fn = 0.5D0*dens%c*w%c*(cs+ABS(cs)) + 0.5D0*dens%n*w%n*(cs-ABS(cs))
!
! ... on Top volume boundary
!
      cs = 0.5D0*(w%c+w%t)
      ft = 0.5D0*dens%c*w%c*(cs+ABS(cs)) + 0.5D0*dens%t*w%t*(cs-ABS(cs))
!
      RETURN
      END SUBROUTINE flw_3d
!------------------------------------------------------    
      SUBROUTINE muscl_flw_3d(fe, fn, ft, dens, u, v, w, i, j, k)
!
! ... Compute the convective fluxes on East, North, and Top sides of the cell
! ... for the momentum density along z.
!
      USE dimensions, ONLY: nx, ny, nz
      USE flux_limiters, ONLY: limiters, lv
      USE grid, ONLY: dx, dy, dz, indz, flag
      USE set_indexes, ONLY: stencil
      USE time_parameters, ONLY: dt
      IMPLICIT NONE
!
      REAL*8, INTENT(INOUT) :: fe, fn, ft
      TYPE(stencil), INTENT(IN) :: dens, u, v, w
      INTEGER, INTENT(IN) :: i, j, k
!
      REAL*8 :: dxm, dxp, dxpp, indxpp, indxp, indxm
      REAL*8 :: dym, dyp, dypp, indypp, indyp, indym
      REAL*8 :: dzp, indzp, dzm, indzm
      REAL*8 :: gradc, grade, gradw, gradn, grads, gradt, gradb
!
      INTEGER :: ip2, jp2, kp2
!
      ip2 = MIN( nx, i+2 )
      jp2 = MIN( ny, j+2 )
      kp2 = MIN( nz, k+2 )
!
      dxm  = dx(i)  + dx(i-1)
      dxp  = dx(i)  + dx(i+1)
      dxpp = dx(i+1)+ dx(ip2)
      dym  = dy(j)  + dy(j-1)
      dyp  = dy(j)  + dy(j+1)
      dypp = dy(j+1)+ dy(jp2)
      dzp  = dz(k)  + dz(k+1)

      indxm  = 1.D0/dxm
      indxp  = 1.D0/dxp
      indxpp = 1.D0/dxpp
      indym  = 1.D0/dym
      indyp  = 1.D0/dyp
      indypp = 1.D0/dypp
      indzp  = 1.D0/dzp
!
! ... MUSCL reconstruction of momentum
!
! ... on East volume boundary
!
      gradw = 2.D0 * indxm  * (w%c *dens%c  - w%w*dens%w)
      gradc = 2.D0 * indxp  * (w%e *dens%e  - w%c*dens%c)
      grade = 2.D0 * indxpp * (w%ee*dens%ee - w%e*dens%e)
!
      lim  = 0.D0
      erre = 0.D0
!
      cs = ( dz(k+1) * u%c + dz(k) * u%t ) * indzp
      !cn = cs * dt * 2.D0 * indxp
      IF (cs >= 0.D0) THEN
        IF (gradc /= 0) erre = gradw / gradc
	incr = 0.5D0 * dx(i)
      ELSE IF (cs < 0.D0) THEN
        IF (gradc /= 0) erre = grade / gradc
	incr = - 0.5D0 * dx(i+1)
      END IF
!
      CALL limiters(lv,lim,erre)
!
      upwnd = lim * gradc * incr
!
      fe = fe + upwnd * cs
!
! ... on North volume boundary
!
      grads = 2.D0 * indym  * (w%c *dens%c  - w%s*dens%s)
      gradc = 2.D0 * indyp  * (w%n *dens%n  - w%c*dens%c)
      gradn = 2.D0 * indypp * (w%nn*dens%nn - w%n*dens%n)
!
      lim  = 0.D0
      erre = 0.D0
!
      cs = ( dz(k+1) * v%c + dz(k) * v%t ) * indzp
      IF (cs >= 0.D0) THEN
        IF (gradc /= 0) erre = grads / gradc
        incr = 0.5D0 * dy(j)
      ELSE IF (cs < 0.D0) THEN
        IF (gradc /= 0) erre = gradn / gradc
        incr = - 0.5D0 * dy(j+1)
      END IF
!
      CALL limiters(lv,lim,erre)
!
      upwnd = lim * gradc * incr
!
      fn = fn + upwnd * cs
!
! ... on Top volume boundary
!
      gradb = indz(k)   * (w%c *dens%c  - w%b*dens%b)
      gradc = indz(k+1) * (w%t *dens%t  - w%c*dens%c)
      gradt = indz(kp2) * (w%tt*dens%tt - w%t*dens%t)
!
      lim  = 0.D0
      erre = 0.D0
!
      cs = 0.5D0 * ( w%c + w%t )
      !cn = cs * dt * indz(k+1)
      IF (cs >= 0.D0) THEN
        IF (gradc /= 0) erre = gradb / gradc
        incr = 0.5D0 * dz(k+1)
      ELSE IF (cs < 0.D0) THEN
        IF (gradc /= 0) erre = gradt / gradc
        incr = - 0.5D0 * dz(k+1)
      END IF 
!
      CALL limiters(lv,lim,erre)
!
      upwnd = lim * gradc * incr
!
      ft = ft + upwnd * cs
!
      RETURN
      END SUBROUTINE muscl_flw_3d
!------------------------------------------------------    
      SUBROUTINE flw_2d(fe, ft, fw, fb, dens, u, w, i, k)
!
! ... Compute the convective fluxes on East, Top, sides of the cell
! ... for the momentum density along z.
!
      USE grid, ONLY: dz, rb, flag
      USE set_indexes, ONLY: imjk, ijkm
      USE set_indexes, ONLY: stencil

      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: fe, ft, fw, fb
      TYPE(stencil), INTENT(IN) :: dens, u, w
      INTEGER, INTENT(IN) :: i, k

      REAL*8 :: dzp, indzp
!
      dzp=dz(k)+dz(k+1)
      indzp=1.D0/dzp
!
! ... on West volume bondary
!
      IF(  .NOT.BTEST(flag(imjk),0)  ) THEN
        cs = (u%wt * dz(k) + u%w * dz(k+1)) * indzp
        fw = 0.5D0*dens%w*w%w*(cs+ABS(cs))*rb(i-1) + 0.5D0*dens%c*w%c*(cs-ABS(cs))*rb(i-1)
      END IF
!
! ... on Bottom volume bondary
!
      IF(  .NOT.BTEST(flag(ijkm),0)  ) THEN
        cs = 0.5D0 * ( w%c + w%b ) 
        fb = 0.5D0*dens%b*w%b*(cs+ABS(cs)) + 0.5D0*dens%c*w%c*(cs-ABS(cs))
      END IF
!
! ... on East volume boundary
!
      cs = indzp * (u%c * dz(k+1) + u%t * dz(k))
      fe = 0.5D0*dens%c*w%c*(cs+ABS(cs))*rb(i) + 0.5D0*dens%e*w%e*(cs-ABS(cs))*rb(i)
!
! ... on Top volume boundary
!
      cs = 0.5D0 * ( w%t + w%c )
      ft = 0.5D0*dens%c*w%c*(cs+ABS(cs)) + 0.5D0*dens%t*w%t*(cs-ABS(cs))
!
      RETURN
      END SUBROUTINE flw_2d
!------------------------------------------------------
      SUBROUTINE muscl_flw_2d(fe, ft, dens, u, w, i, k)
!
! ... Compute the convective fluxes on East, Top, sides of the cell
! ... for the momentum density along z.
!
      USE dimensions, ONLY: nx, ny, nz
      USE flux_limiters, ONLY: limiters, lv
      USE grid, ONLY: dx, rb, dz, indz
      USE set_indexes, ONLY: stencil
      USE time_parameters, ONLY: dt

      IMPLICIT NONE
!
      REAL*8, INTENT(INOUT) :: fe, ft
      TYPE(stencil), INTENT(IN) :: dens, u, w
      INTEGER, INTENT(IN) :: i, k

      REAL*8 :: dxm, dxp, dxpp, indxpp, indxp, indxm
      REAL*8 :: dzp, indzp
      REAL*8 :: gradc, grade, gradw, gradt, gradb

      INTEGER :: ip2, kp2
!
      ip2 = MIN( nx, i+2 )
      kp2 = MIN( nz, k+2 )
!
      dxm  = dx(i)  + dx(i-1)
      dxp  = dx(i)  + dx(i+1)
      dxpp = dx(i+1)+ dx(ip2)
      dzp  = dz(k)  + dz(k+1)

      indxm  = 1.D0/dxm
      indxp  = 1.D0/dxp
      indxpp = 1.D0/dxpp
      indzp  = 1.D0/dzp
!
! ... MUSCL reconstruction of momentum
!
! ... on East volume boundary
!
      gradc = 2.D0 * indxp  * (dens%e*w%e   - dens%c*w%c)
      gradw = 2.D0 * indxm  * (dens%c*w%c   - dens%w*w%w)
      grade = 2.D0 * indxpp * (dens%ee*w%ee - dens%e*w%e)
!
      lim  = 0.D0
      erre = 0.D0
!
      cs = indzp * (u%c * dz(k+1) + u%t * dz(k))
      !cn = cs * dt * 2.D0 * indzp
      IF ( cs >= 0.D0 ) THEN
        IF (gradc /= 0.D0) erre = gradw / gradc
        incr = 0.5D0 * dx(i)
      ELSE IF ( cs < 0.D0 ) THEN
        IF (gradc /= 0.D0) erre = grade / gradc
        incr = - 0.5D0 * dx(i+1)
      END IF
!
      CALL limiters(lv,lim,erre)
!
      upwnd = lim * gradc * incr
!
      fe = fe + upwnd * cs * rb(i)
!
! ... on Top volume boundary
!
      gradc = (dens%t*w%t   - dens%c*w%c) * indz(k+1)
      gradb = (dens%c*w%c   - dens%b*w%b) * indz(k)
      gradt = (dens%tt*w%tt - dens%t*w%t) * indz(kp2)
!
      lim  = 0.D0
      erre = 0.D0
!
      cs = 0.5D0 * ( w%t + w%c )
      !cn = cs * dt * indz(k+1)
      IF (cs >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb / gradc
        incr = 0.5D0 * dz(k+1)
      ELSE IF (cs < 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradt / gradc
        incr = - 0.5D0 * dz(k+1)
      END IF 
!
      CALL limiters(lv,lim,erre)
!
      upwnd = lim * gradc * incr
!
      ft = ft + upwnd * cs
!
      RETURN
      END SUBROUTINE muscl_flw_2d
!-----------------------------------------------------------------------
      SUBROUTINE ctu1_flw_2d(fe, ft, dens, u, w, i, k)
!
! ... Compute the first order Corner Transport Upwind correction (step 2 in LeVeque algorithm)
!
      USE dimensions
      USE domain_mapping, ONLY: myijk
      USE grid, ONLY: dx, dz, indz, rb
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: stencil
      USE time_parameters, ONLY: dt
      USE set_indexes, ONLY: imjk, ijkm

      REAL*8, INTENT(INOUT) :: fe, ft
      TYPE(stencil), INTENT(IN) :: dens, u, w
      INTEGER, INTENT(IN) :: i, k

      REAL*8 :: wc, uc, deltar
      REAL*8 :: dzp, dzpp, dxm, dxp, indxp, indxm, indzp, indzpp
      REAL*8 :: incrxp, incrxm, incrzp, incrzm
      INTEGER kp2
!
      kp2 = MIN( nz, k+2)
!
      dzp  = dz(k)  + dz(k+1)
      dzpp = dz(k+1)+ dz(kp2)
      dxm  = dx(i)  + dx(i-1)
      dxp  = dx(i)  + dx(i+1)
!
      indzp  = 2.D0/dzp
      indzpp = 2.D0/dzpp
      indxm  = 2.D0/dxm
      indxp  = 2.D0/dxp
!
      incrxp = 0.125D0*dt*indxp
      incrxm = 0.125D0*dt*indxm
      incrzp = 0.125D0*dt*indz(k+1)
      incrzm = 0.125D0*dt*indz(k)
!
! c-w
      uc = 0.5D0*indzp*(dz(k)*u%wt+dz(k+1)*u%w)
      wc = 0.5D0*indxm*(dx(i-1)*w%c+dx(i)*w%w)
      deltar  = dens%c*w%c - dens%w*w%w 
      ft = ft - incrxm*(uc+ABS(uc))*(wc+ABS(wc))*deltar
!
! e-c
      uc = 0.5D0*indzp*(dz(k)*u%t+dz(k+1)*u%c)
      wc = 0.5D0*indxp*(dx(i+1)*w%c+dx(i)*w%e)
      deltar  = dens%e*w%e - dens%c*w%c 
      ft = ft - incrxp*(uc-ABS(uc))*(wc+ABS(wc))*deltar
!
! t-wt
      uc = 0.5D0*indzpp*(dz(kp2)*u%wt+dz(k+1)*u%wtt)
      wc = 0.5D0*indxm *(dx(i-1)*w%t+dx(i)*w%wt)
      deltar  = dens%t*w%t - dens%wt*w%wt 
      IF (k >= nz-1) deltar = 0.D0
      ft = ft - incrxm*(uc+ABS(uc))*(wc-ABS(wc))*deltar
!
! et-t
      uc = 0.5D0*indzpp*(dz(kp2)*u%t+dz(k+1)*u%tt)
      wc = 0.5D0*indxp *(dx(i+1)*w%t+dx(i)*w%et)
      deltar  = dens%et*w%et - dens%t*w%t 
      IF (k >= nz-1) deltar = 0.D0
      ft = ft - incrxp*(uc-ABS(uc))*(wc-ABS(wc))*deltar
!
! c-b
      uc = 0.5D0*(u%c+u%w)
      wc = 0.5D0*(w%c+w%b)
      deltar  = dens%c*w%c - dens%b*w%b
      fe = fe - incrzm*(uc+ABS(uc))*(wc+ABS(wc))*deltar*rb(i)
!
! t-c
      uc = 0.5D0*(u%t+u%wt)
      wc = 0.5D0*(w%c+w%t)
      deltar  = dens%t*w%t - dens%c*w%c
      fe = fe - incrzp*(uc+ABS(uc))*(wc-ABS(wc))*deltar*rb(i)
!
! e-b
      uc = 0.5D0*(u%c+u%e)
      wc = 0.5D0*(w%e+w%eb)
      deltar  = dens%e*w%e - dens%eb*w%eb
      IF (i >= nx-1) deltar = 0.D0
      fe = fe - incrzm*(uc-ABS(uc))*(wc+ABS(wc))*deltar*rb(i)
!
! et-e
      uc = 0.5D0*(u%t+u%et)
      wc = 0.5D0*(w%e+w%et)
      deltar  = dens%et*w%et - dens%e*w%e
      IF (i >= nx-1) deltar = 0.D0
      fe = fe - incrzp*(uc-ABS(uc))*(wc-ABS(wc))*deltar*rb(i)
!
      RETURN
      END SUBROUTINE ctu1_flw_2d
!------------------------------------------------------------------------
    SUBROUTINE ctu2_flw_2d(fe, ft, dens, u, w, i, k)
!
! ... Compute the second order Corner Transport Upwind correction (step 3 in LeVeque algorithm)
!
      USE dimensions
      USE domain_mapping, ONLY: myijk
      USE grid, ONLY: dx, dz, indx, indz, flag, fluid, rb
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: stencil
      USE time_parameters, ONLY: dt
      USE set_indexes, ONLY: imjk, ijkm, ipjk, ijkp
      USE flux_limiters, ONLY: limiters, lv
!
      REAL*8, INTENT(INOUT) :: fe, ft
      TYPE(stencil), INTENT(IN) :: dens, u, w
      INTEGER, INTENT(IN) :: i, k
!
      REAL*8 :: delta, indzp, indxp, indxm, indxpp
      REAL*8 :: dxp, dxm, dxpp, dzp
      INTEGER :: ip2, kp2
      REAL*8 :: gradc, grade, gradw, gradt, gradb, deltar, uref, wref
!
      ip2 = MIN( nx, i+2 )
      kp2 = MIN( nz, k+2 )
!
      dzp  = dz(k)  + dz(k+1)
      dxm  = dx(i)  + dx(i-1)
      dxp  = dx(i)  + dx(i+1)
      dxpp = dx(i+1)+ dx(ip2)
!
      indzp  = 2.D0/dzp 
      indxm  = 2.D0/dxm
      indxp  = 2.D0/dxp
      indxpp = 2.D0/dxpp
!
! e-c
      gradc = indxp  * (dens%e *w%e  - dens%c*w%c)
      grade = indxpp * (dens%ee*w%ee - dens%e*w%e)
      gradw = indxm  * (dens%c *w%c  - dens%w*w%w)
      lim  = 0.D0
      erre = 0.D0
      uref = 0.5D0*indzp*(dz(k+1)*u%c + dz(k)*u%t)
      IF (uref >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=nx-1 .AND. flag(ipjk)==fluid) CALL limiters(lv,lim,erre)
      uref = ABS(0.5D0*indzp*(dz(k+1)*u%c + dz(k)*u%t))
      deltar = dens%e*w%e - dens%c*w%c
      fe = fe + 0.5D0*uref*(1-dt*indxp*uref)*deltar*lim*rb(i)
!
! t-c
      gradc = indz(k+1) * (dens%t *w%t  - dens%c*w%c)
      gradt = indz(kp2) * (dens%tt*w%tt - dens%t*w%t)
      gradb = indz(k)   * (dens%c *w%c  - dens%b*w%b)
      lim  = 0.D0
      erre = 0.D0
      wref = 0.5D0*(w%t + w%c)
      IF (wref >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(ijkp)==fluid) CALL limiters(lv,lim,erre)
      wref = ABS(0.5D0*(w%c + w%t))
      deltar = dens%t*w%t - dens%c*w%c
      ft = ft + 0.5D0*wref*(1-dt*indz(k+1)*wref)*deltar*lim
!
      RETURN
      END SUBROUTINE ctu2_flw_2d
!-----------------------------------------------------------------------
      SUBROUTINE ctu3_flw_2d(fe, ft, dens, u, w, i, k)
!
! ... Compute the second order Corner Transport Upwind correction (step 4 in LeVeque algorithm)
!
      USE dimensions
      USE domain_mapping, ONLY: myijk
      USE grid, ONLY: dx, dz, indx, indz, fluid, flag, rb
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: stencil, ipjk, imjk, ijkp, ijkm
      USE set_indexes, ONLY: ipjkp, ipjkm, imjkp, imjkm
      USE time_parameters, ONLY: dt
      USE flux_limiters, ONLY: limiters, lv
!
      IMPLICIT NONE
!
      REAL*8, INTENT(INOUT) :: fe, ft
      TYPE(stencil), INTENT(IN) :: u, w, dens
      INTEGER, INTENT(IN) :: i, k

      REAL*8 :: uc, wc, deltar, lim, erre
      REAL*8 :: gradc, grade, gradw, gradt, gradb
      REAL*8 :: dxp, dxm, dzp, dzm, dxpp, dxmm, dzpp, dzmm
      REAL*8 :: indxp, indxm, indzp, indzm, indxpp, indxmm, indzpp, indzmm
      REAL*8 :: incrxp, incrxpp, incrzp, incrzpp
      INTEGER :: im2, ip2, km2, kp2

      im2 = MAX(  1, i-2 )
      ip2 = MIN( nx, i+2 )
      km2 = MAX(  1, k-2 )
      kp2 = MIN( nz, k+2 )
!
      dzm  = dz(k) + dz(k-1)
      dzp  = dz(k) + dz(k+1)
      dxm  = dx(i) + dx(i-1)
      dxp  = dx(i) + dx(i+1)
      dxpp = dx(i+1)+dx(ip2) 
      dxmm = dx(i-1)+dx(im2)
      dzpp = dz(k+1)+dz(kp2)
      dzmm = dz(k-1)+dz(km2)
!
      indzm  = 2.D0/dzm
      indzp  = 2.D0/dzp
      indxm  = 2.D0/dxm
      indxp  = 2.D0/dxp
      indxpp = 2.D0/dxpp
      indxmm = 2.D0/dxmm
      indzpp = 2.D0/dzpp
      indzmm = 2.D0/dzmm
!
      incrxp  = 0.5D0*dt*indx(i)
      incrxpp = 0.5D0*dt*indx(i+1)
      incrzp  = 0.5D0*dt*indzp
      incrzpp = 0.5D0*dt*indzpp
!
! c-w
      uc = 0.5D0*indzp*(dz(k)*u%wt+dz(k+1)*u%w)
      wc = 0.5D0*indxm*(dx(i-1)*w%c+dx(i)*w%w)
      deltar = 0.5D0*ABS(uc)*(1 - dt*indxm*ABS(uc))*(w%c*dens%c - w%w*dens%w)
      gradc = indxm  * (w%c*dens%c - w%w*dens%w  )
      grade = indxp  * (w%e*dens%e - w%c*dens%c  )
      gradw = indxmm * (w%w*dens%w - w%ww*dens%ww)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=2 .AND. flag(imjk)==fluid) CALL limiters(lv,lim,erre)
      ft = ft + incrzp*(wc+ABS(wc))*deltar*lim
!
! e-c
      uc = 0.5D0*indzp*(dz(k)*u%t+dz(k+1)*u%c)
      wc = 0.5D0*indxp*(dx(i+1)*w%c+dx(i)*w%e)
      deltar = 0.5D0*ABS(uc)*(1 - dt*indxp*ABS(uc))*(w%e*dens%e - w%c*dens%c)
      gradc = indxp  * (w%e *dens%e  - w%c*dens%c)
      grade = indxpp * (w%ee*dens%ee - w%e*dens%e)
      gradw = indxm  * (w%c *dens%c  - w%w*dens%w)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=nx-1 .AND. flag(ipjk)==fluid) CALL limiters(lv,lim,erre)
      ft = ft - incrzp*(wc+ABS(wc))*deltar*lim
!
! t-wt
      uc = 0.5D0*indzpp*(dz(kp2)*u%wt+dz(k+1)*u%wtt)
      wc = 0.5D0*indxm*(dx(i-1)*w%t+dx(i)*w%wt)
      deltar = 0.5D0*ABS(uc)*(1 - dt*indxm*ABS(uc))*(w%t*dens%t - w%wt*dens%wt)
      IF (k >= nz-1) deltar = 0.D0
      gradc = indxm  * (w%t *dens%t  - w%wt *dens%wt )
      grade = indxp  * (w%et*dens%et - w%t  *dens%t  )
      gradw = indxmm * (w%wt*dens%wt - w%wwt*dens%wwt)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.0D0) erre = grade/gradc
      END IF
      IF (i/=2 .AND. flag(imjkp)==fluid) CALL limiters(lv,lim,erre)
      ft = ft + incrzpp*(wc-ABS(wc))*deltar*lim
!
! et-t
      uc = 0.5D0*indzpp*(dz(kp2)*u%t+dz(k+1)*u%tt)
      wc = 0.5D0*indxp*(dx(i+1)*w%t+dx(i)*w%et)
      deltar = 0.5D0*ABS(uc)*(1 - dt*indxp*ABS(uc))*(w%et*dens%et - w%t*dens%t)
      IF (k >= nz-1) deltar = 0.D0
      gradc = indxp  * (w%et *dens%et  - w%t *dens%t )
      grade = indxpp * (w%eet*dens%eet - w%et*dens%et)
      gradw = indxm  * (w%t  *dens%t   - w%wt*dens%wt)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=nx-1 .AND. flag(ipjkp)==fluid) CALL limiters(lv,lim,erre)
      ft = ft - incrzpp*(wc-ABS(wc))*deltar*lim
!
! c-b
      uc = 0.5D0*(u%c+u%w)
      wc = 0.5D0*(w%c+w%b)
      deltar = 0.5D0*ABS(wc)*(1 - dt*indz(k)*ABS(wc))*(w%c*dens%c - w%b*dens%b)
      gradc = indz(k)   * (w%c*dens%c - w%b *dens%b )
      gradt = indz(k+1) * (w%t*dens%t - w%c *dens%c )
      gradb = indz(k-1) * (w%b*dens%b - w%bb*dens%bb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=2 .AND. flag(ijkm)==fluid) CALL limiters(lv,lim,erre)
      fe = fe + incrxp*(uc+ABS(uc))*deltar*lim*rb(i)
!
! t-c
      uc = 0.5D0*(u%t+u%wt)
      wc = 0.5D0*(w%c+w%t)
      deltar = 0.5D0*ABS(wc)*(1 - dt*indz(k+1)*ABS(wc))*(w%t*dens%t - w%c*dens%c)
      gradc = indz(k+1) * (w%t *dens%t  - w%c*dens%c)
      gradt = indz(kp2) * (w%tt*dens%tt - w%t*dens%t)
      gradb = indz(k)   * (w%c *dens%c  - w%b*dens%b)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(ijkp)==fluid) CALL limiters(lv,lim,erre)
      fe = fe - incrxp*(uc+ABS(uc))*deltar*lim*rb(i)
!
! et-e
      uc = 0.5D0*(u%t+u%et)
      wc = 0.5D0*(w%e+w%et)
      deltar = 0.5D0*ABS(wc)*(1 - dt*indz(k+1)*ABS(wc))*(w%et*dens%et - w%e*dens%e)
      IF (i >= nx-1) deltar = 0.D0
      gradc = indz(k+1) * (w%et *dens%et  - w%e *dens%e )
      gradt = indz(kp2) * (w%ett*dens%ett - w%et*dens%et)
      gradb = indz(k-1) * (w%e  *dens%e   - w%eb*dens%eb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(ijkp)==fluid) CALL limiters(lv,lim,erre)
      fe = fe - incrxpp*(uc-ABS(uc))*deltar*lim*rb(i)
!
! e-eb
      uc = 0.5D0*(u%c+u%e)
      wc = 0.5D0*(w%e+w%eb)
      deltar = 0.5D0*ABS(wc)*(1 - dt*indz(k)*ABS(wc))*(w%e*dens%e - w%eb*dens%eb)
      IF (i >= nx-1) deltar = 0.D0
      gradc = indz(k)   * (w%e *dens%e  - w%eb *dens%eb )
      gradt = indz(k+1) * (w%et*dens%et - w%e  *dens%e  )
      gradb = indz(k-1) * (w%eb*dens%eb - w%ebb*dens%ebb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=2 .AND. flag(ijkm)==fluid) CALL limiters(lv,lim,erre)
      fe = fe + incrxpp*(uc-ABS(uc))*deltar*lim*rb(i)
!
      RETURN
      END SUBROUTINE ctu3_flw_2d
!----------------------------------------------------------------------
      SUBROUTINE ctu1_flw_3d(fe, fn, ft, dens, u, v, w, i, j, k)
!
! ... Compute the first order Corner Transport Upwind correction (step 2 in LeVeque algorithm)
!
      USE dimensions
      USE domain_mapping, ONLY: myijk, meshinds
      USE grid, ONLY:  dx, dy, dz, indx, indy, indz
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: stencil
      USE time_parameters, ONLY: dt

      IMPLICIT NONE

      REAL*8, INTENT(INOUT) :: fe, fn, ft
      TYPE(stencil), INTENT(IN) :: u, v, w, dens
      INTEGER, INTENT(IN) :: i, j, k

      REAL*8 :: uc, vc, wc, deltar, onethird
      REAL*8 :: dxp, dxm, dyp, dym, dzp, dzm, dzpp
      REAL*8 :: dxi, dxip, dxim, dyj, dyjp, dyjm, dzk, dzkp, dzkpp, dzkm
      REAL*8 :: indxp, indxm, indyp, indym, indzp, indzm, indzpp
      REAL*8 :: incrxp, incrxm, incryp, incrym, incrzp, incrzm, incrzpp
      INTEGER :: kp2
!
      onethird = 1.D0/3.D0
      kp2 = MIN( nz, k+2 )
!
      dxi   = dx(i)
      dxip  = dx(i+1)
      dxim  = dx(i-1)
      dyj   = dy(j)
      dyjp  = dy(j+1)
      dyjm  = dy(j-1)
      dzk   = dz(k)
      dzkp  = dz(k+1)
      dzkpp = dz(kp2)
      dzkm  = dz(k-1)
!
      dxm  = dxi  + dxim
      dxp  = dxi  + dxip
      dym  = dyj  + dyjm
      dyp  = dyj  + dyjp
      dzm  = dzk  + dzkm
      dzp  = dzk  + dzkp
      dzpp = dzkp + dzkpp
!
      indxp  = 2.D0/dxp
      indxm  = 2.D0/dxm
      indyp  = 2.D0/dyp
      indym  = 2.D0/dym
      indzp  = 2.D0/dzp
      indzm  = 2.D0/dzm
      indzpp = 2.D0/dzpp
!
      incrxp  = 0.5D0*dt*indxp
      incrxm  = 0.5D0*dt*indxm
      incryp  = 0.5D0*dt*indyp
      incrym  = 0.5D0*dt*indym
      incrzp  = 0.125D0*dt*indzp
      incrzm  = 0.125D0*dt*indzm
      incrzpp = 0.125D0*dt*indzpp
!
! t-c
      deltar = dens%t*w%t - dens%c*w%c
      uc = 0.5D0*(u%t+u%wt)*dt*indx(i)
      vc = 0.5D0*(v%t+v%st)*dt*indy(j)
      wc = 0.5D0*(w%t+w%c)
      fe = fe - 0.125D0*(uc+ABS(uc))*(wc-ABS(wc))*deltar
      fn = fn - 0.125D0*(vc+ABS(vc))*(wc-ABS(wc))*deltar &
              + 0.25D0*onethird*(vc+ABS(vc))*(wc-ABS(wc))*ABS(uc)*deltar 
!
! c-b
      deltar = dens%c*w%c - dens%b*w%b
      uc = 0.5D0*(u%c+u%w)*dt*indx(i)
      vc = 0.5D0*(v%c+v%s)*dt*indy(j)
      wc = 0.5D0*(w%c+w%b)
      fe = fe - 0.125D0*(uc+ABS(uc))*(wc+ABS(wc))*deltar
      fn = fn - 0.125D0*(vc+ABS(vc))*(wc+ABS(wc))*deltar &
              + 0.25D0*onethird*(vc+ABS(vc))*(wc+ABS(wc))*ABS(uc)*deltar 
!
! n-c
      deltar = dens%n*w%n - dens%c*w%c
      uc = indyp*incrzp*(dyjp*dzkp*(u%c  + u%w  ) + &
                         dyj *dzkp*(u%n  + u%wn ) + & 
                         dyjp*dzk *(u%t  + u%wt ) + & 
                         dyj *dzk *(u%nt + u%wnt))*indx(i)
      vc = 0.5D0*indzp*(dzkp*v%c+dzk*v%t)
      wc = incryp*(dyjp*w%c+dyj*w%n)*indzp
      ft = ft - 0.125D0*(wc+ABS(wc))*(vc-ABS(vc))*deltar
      fe = fe - 0.125D0*(uc+ABS(uc))*(vc-ABS(vc))*deltar &
              + 0.25D0*onethird*(uc+ABS(uc))*(vc-ABS(vc))*ABS(wc)*deltar 
!
! c-s
      deltar = dens%c*w%c - dens%s*w%s
      uc = indym*incrzp*(dyjm*dzkp*(u%c  + u%w  ) + &
                         dyj *dzkp*(u%s  + u%ws ) + & 
                         dyjm*dzk *(u%t  + u%wt ) + & 
                         dyj *dzk *(u%st + u%wst))*indx(i)
      vc = 0.5D0*indzp*(dzkp*v%s+dzk*v%st)
      wc = incrym*(dyjm*w%c+dyj*w%s)*indzp
      ft = ft - 0.125D0*(wc+ABS(wc))*(vc+ABS(vc))*deltar
      fe = fe - 0.125D0*(uc+ABS(uc))*(vc+ABS(vc))*deltar &
              + 0.25D0*onethird*(uc+ABS(uc))*(vc+ABS(vc))*ABS(wc)*deltar 
!
! e-c
      deltar = dens%e*w%e - dens%c*w%c
      uc = 0.5D0*indzp*(dzkp*u%c+dzk*u%t)
      vc = indxp*incrzp*(dxip*dzkp*(v%c + v%s  ) + &
                         dxi *dzkp*(v%e + v%es ) + &
                         dxip*dzk *(v%t + v%st ) + &
                         dxi *dzk *(v%et+ v%est))*indy(j)
      wc = incrxp*(dxip*w%c+dxi*w%e)*indzp
      fn = fn - 0.125D0*(vc+ABS(vc))*(uc-ABS(uc))*deltar
      ft = ft - 0.125D0*(wc+ABS(wc))*(uc-ABS(uc))*deltar &
              + 0.25D0*onethird*(wc+ABS(wc))*(uc-ABS(uc))*ABS(vc)*deltar 
!
! c-w
      deltar = dens%c*w%c - dens%w*w%w
      uc = 0.5D0*indzp*(dzkp*u%w+dzk*u%wt)
      vc = indxm*incrzp*(dxim*dzkp*(v%c + v%s  ) + &
                         dxi *dzkp*(v%w + v%ws ) + &
                         dxim*dzk *(v%t + v%st ) + &
                         dxi *dzk *(v%wt+ v%wst))*indy(j)
      wc = incrxm*(dxim*w%c+dxi*w%w)*indzp
      fn = fn - 0.125D0*(vc+ABS(vc))*(uc+ABS(uc))*deltar
      ft = ft - 0.125D0*(wc+ABS(wc))*(uc+ABS(uc))*deltar &
              + 0.25D0*onethird*(wc+ABS(wc))*(uc+ABS(uc))*ABS(vc)*deltar 
!
! wt-w
      deltar = dens%wt*w%wt - dens%w*w%w
      IF (i <= 2) deltar = 0.D0
      uc = 0.5D0*(u%wt+u%wwt)*dt*indx(i-1)
      vc = 0.5D0*(v%wt+v%wst)*dt*indy(j)
      wc = 0.5D0*(w%wt+w%w)
      fn = fn - 0.125D0*onethird*(wc-ABS(wc))*(vc+ABS(vc))*ABS(uc+ABS(uc))*deltar
!
! w-wb
      deltar = dens%w*w%w - dens%wb*w%wb
      IF (i <= 2) deltar = 0.D0
      uc = 0.5D0*(u%w+u%ww)*dt*indx(i-1)
      vc = 0.5D0*(v%w+v%ws)*dt*indy(j)
      wc = 0.5D0*(w%w+w%wb)
      fn = fn - 0.125D0*onethird*(wc+ABS(wc))*(vc+ABS(vc))*ABS(uc+ABS(uc))*deltar
!
! et-e
      deltar = dens%et*w%et - dens%e*w%e
      IF (i >= nx-1) deltar = 0.D0
      uc = 0.5D0*(u%et+u%t)*dt*indx(i+1)
      vc = 0.5D0*(v%et+v%est)*dt*indy(j)
      wc = 0.5D0*(w%et+w%e)
      fe = fe - 0.125D0*(wc-ABS(wc))*(uc-ABS(uc))*deltar
      fn = fn - 0.125D0*onethird*(wc-ABS(wc))*(vc+ABS(vc))*ABS(uc-ABS(uc))*deltar
!
! e-eb
      deltar = dens%e*w%e - dens%eb*w%eb
      IF (i >= nx-1) deltar = 0.D0
      uc = 0.5D0*(u%e+v%c)*dt*indx(i+1)
      vc = 0.5D0*(v%e+v%es)*dt*indy(j)
      wc = 0.5D0*(w%e+w%eb)
      fe = fe - 0.125D0*(wc+ABS(wc))*(uc-ABS(uc))*deltar
      fn = fn - 0.125D0*onethird*(wc+ABS(wc))*(vc+ABS(vc))*ABS(uc-ABS(uc))*deltar
!
! en-e
      deltar = dens%en*w%en - dens%e*w%e
      IF (i >= nx-1) deltar = 0.D0
      uc = indyp*incrzp*(dyjp*dzkp*(u%c  + u%e  ) + &
                         dyj *dzkp*(u%n  + u%en ) + & 
                         dyjp*dzk *(u%t  + u%et ) + & 
                         dyj *dzk *(u%nt + u%ent))*indx(i+1)
      vc = 0.5D0*indzp*(dzkp*v%e+dzk*v%et)
      wc = incryp*(dyjp*w%e+dyj*w%en)*indzp
      fe = fe - 0.125D0*(vc-ABS(vc))*(uc-ABS(uc))*deltar &
              + 0.25D0*onethird*(vc-ABS(vc))*(uc-ABS(uc))*ABS(wc)*deltar
!
! e-es
      deltar = dens%e*w%e - dens%es*w%es
      IF (i >= nx-1) deltar = 0.D0
      uc = indym*incrzp*(dyjm*dzkp*(u%c  + u%e  ) + &
                         dyj *dzkp*(u%s  + u%es ) + & 
                         dyjm*dzk *(u%t  + u%et ) + & 
                         dyj *dzk *(u%st + u%est))*indx(i+1)
      vc = 0.5D0*indzp*(dzkp*v%es+dzk*v%est)
      wc = incrym*(dyjm*w%e+dyj*w%es)*indzp
      fe = fe - 0.125D0*(vc+ABS(vc))*(uc-ABS(uc))*deltar &
              + 0.25D0*onethird*(vc+ABS(vc))*(uc-ABS(uc))*ABS(wc)*deltar
!
! nt-n
      deltar  = dens%nt*w%nt - dens%n*w%n
      IF (j >= ny-1) deltar = 0.D0
      uc = 0.5D0*(u%nt+u%wnt)*dt*indx(i)
      vc = 0.5D0*(v%nt+v%t)*dt*indy(j+1)
      wc = 0.5D0*(w%nt+w%n)
      fn = fn - 0.125D0*(vc-ABS(vc))*(wc-ABS(wc))*deltar &
              + 0.25D0*onethird*(vc-ABS(vc))*(wc-ABS(wc))*ABS(uc)*deltar
!
! n-nb
      deltar  = dens%n*w%n - dens%nb*w%nb
      IF (j >= ny-1) deltar = 0.D0
      uc = 0.5D0*(u%n+u%wn)*dt*indx(i)
      vc = 0.5D0*(v%n+v%c)*dt*indy(j+1)
      wc = 0.5D0*(w%n+w%nb)
      fn = fn - 0.125D0*(vc-ABS(vc))*(wc+ABS(wc))*deltar &
              + 0.25D0*onethird*(vc-ABS(vc))*(wc+ABS(wc))*ABS(uc)*deltar
!
! nt-t
      deltar  = dens%nt*w%nt - dens%t*w%t
      IF (k >= nz-1) deltar = 0.D0
      uc = indyp*incrzpp*(dyjp*dzkpp*(u%t  + u%wt  ) + &
                          dyj *dzkpp*(u%nt + u%wnt ) + & 
                          dyjp*dzkp *(u%tt + u%wtt ) + & 
                          dyj *dzkp *(u%ntt+ u%wntt))*indx(i)
      vc = 0.5D0*indzpp*(dzkpp*v%t+dzkp*v%tt)
      wc = incryp*(dyjp*w%t+dyj*w%nt)*indzpp
      ft = ft - 0.125D0*(vc-ABS(vc))*(wc-ABS(wc))*deltar
      fe = fe - 0.125D0*onethird*(vc-ABS(vc))*(uc+ABS(uc))*ABS(wc-ABS(wc))*deltar
!     
! t-st
      deltar  = dens%t*w%t - dens%st*w%st
      IF (k >= nz-1) deltar = 0.D0
      uc = indym*incrzpp*(dyjm*dzkpp*(u%t  + u%wt  ) + &
                          dyj *dzkpp*(u%st + u%wst ) + & 
                          dyjm*dzkp *(u%tt + u%wtt ) + & 
                          dyj *dzkp *(u%stt+ u%wstt))*indx(i)
      vc = 0.5D0*indzpp*(dzkpp*v%st+dzkp*v%stt)
      wc = incrym*(dyjm*w%t+dyj*w%st)*indzpp
      ft = ft - 0.125D0*(vc+ABS(vc))*(wc-ABS(wc))*deltar
      fe = fe - 0.125D0*onethird*(vc+ABS(vc))*(uc+ABS(uc))*ABS(wc-ABS(wc))*deltar
!     
! nb-b
      deltar  = dens%nb*w%nb - dens%b*w%b
      IF (k <= 2) deltar = 0.D0
      uc = indyp*incrzm*(dyjp*dzkm*(u%c  + u%w  ) + &
                         dyj *dzkm*(u%n  + u%wn ) + & 
                         dyjp*dzk *(u%b  + u%wb ) + & 
                         dyj *dzk *(u%nb + u%wnb))*indx(i)
      vc = 0.5D0*indzm*(dzkm*v%c+dzk*v%b)
      wc = incryp*(dyjp*w%b+dyj*w%nb)*indzm
      fe = fe - 0.125D0*onethird*(vc-ABS(vc))*(uc+ABS(uc))*ABS(wc+ABS(wc))*deltar
!     
! b-sb
      deltar  = dens%b*w%b - dens%sb*w%sb
      IF (k <= 2) deltar = 0.D0
      uc = indym*incrzm*(dyjm*dzkm*(u%c  + u%w  ) + &
                         dyj *dzkm*(u%s  + u%ws ) + & 
                         dyjm*dzk *(u%b  + u%wb ) + & 
                         dyj *dzk *(u%sb + u%wsb))*indx(i)
      vc = 0.5D0*indzm*(dzkm*v%s+dzk*v%sb)
      wc = incrym*(dyjm*w%b+dyj*w%sb)*indzm
      fe = fe - 0.125D0*onethird*(vc+ABS(vc))*(uc+ABS(uc))*ABS(wc+ABS(wc))*deltar
!     
! wn-n
      deltar  = dens%en*w%en - dens%n*w%n
      IF (j >= ny-1) deltar = 0.D0
      uc = 0.5D0*indzp*(dzkp*u%n+dzk*u%nt)
      vc = indxp*incrzp*(dxip*dzkp*(v%c + v%n  ) + &
                         dxi *dzkp*(v%e + v%en ) + &
                         dxip*dzk *(v%t + v%nt ) + &
                         dxi *dzk *(v%et+ v%ent))*indy(j+1)
      wc = incrxp*(dxip*w%n+dxi*w%en)*indzp
      fn = fn - 0.125D0*(vc-ABS(vc))*(uc-ABS(uc))*deltar
      ft = ft - 0.125D0*onethird*(uc-ABS(uc))*(wc+ABS(wc))*ABS(vc-ABS(vc))*deltar
!     
! n-wn
      deltar  = dens%n*w%n - dens%wn*w%wn
      IF (j >= ny-1) deltar = 0.D0
      uc = 0.5D0*indzp*(dzkp*u%wn+dzk*u%wnt)
      vc = indxm*incrzp*(dxim*dzkp*(v%c + v%n  ) + &
                         dxi *dzkp*(v%w + v%wn ) + &
                         dxim*dzk *(v%t + v%nt ) + &
                         dxi *dzk *(v%wt+ v%wnt))*indy(j+1)
      wc = incrxm*(dxim*w%n+dxi*w%wn)*indzp
      fn = fn - 0.125D0*(vc-ABS(vc))*(uc+ABS(uc))*deltar
      ft = ft - 0.125D0*onethird*(uc+ABS(uc))*(wc+ABS(wc))*ABS(vc-ABS(vc))*deltar
!     
! es-s
      deltar  = dens%es*w%es - dens%s*w%s
      IF (j <= 2) deltar = 0.D0
      uc = 0.5D0*indzp*(dzkp*u%s+dzk*u%st)
      vc = indxp*incrzp*(dxip*dzkp*(v%s  + v%ss  ) + &
                         dxi *dzkp*(v%es + v%ess ) + &
                         dxip*dzk *(v%st + v%sst ) + &
                         dxi *dzk *(v%est+ v%esst))*indy(j-1)
      wc = incrxp*(dxip*w%s+dxi*w%es)*indzp
      ft = ft - 0.125D0*onethird*(uc-ABS(uc))*(wc+ABS(wc))*ABS(vc+ABS(vc))*deltar
!     
! s-ws
      deltar  = dens%s*w%s - dens%ws*w%ws
      IF (j <= 2) deltar = 0.D0
      uc = 0.5D0*indzp*(dzkp*u%ws+dzk*u%wst)
      vc = indxm*incrzp*(dxim*dzkp*(v%s + v%ss  ) + &
                         dxi *dzkp*(v%ws + v%wss ) + &
                         dxim*dzk *(v%st + v%sst ) + &
                         dxi *dzk *(v%wst+ v%wsst))*indy(j-1)
      wc = incrxm*(dxim*w%s+dxi*w%ws)*indzp
      ft = ft - 0.125D0*onethird*(uc+ABS(uc))*(wc+ABS(wc))*ABS(vc+ABS(vc))*deltar
!     
! et-t
      deltar = dens%et*w%et - dens%t*w%t
      IF (k >= nz-1) deltar = 0.D0
      uc = 0.5D0*indzpp*(dzkpp*u%t+dzkp*u%tt)
      vc = indxp*incrzpp*(dxip*dzkpp*(v%t  + v%st  ) + &
                          dxi *dzkpp*(v%et + v%est ) + &
                          dxip*dzkp *(v%tt + v%stt ) + &
                          dxi *dzkp *(v%ett+ v%estt))*indy(j)
      wc = incrxp*(dxip*w%t+dxi*w%et)*indzpp
      ft = ft - 0.125D0*(wc-ABS(wc))*(uc-ABS(uc))*deltar &
              + 0.25D0*onethird*(wc-ABS(wc))*(uc-ABS(uc))*ABS(vc)*deltar
!     
! t-wt
      deltar  = dens%t*w%t - dens%wt*w%wt
      IF (k >= nz-1) deltar = 0.D0
      uc = 0.5D0*indzpp*(dzkpp*u%wt+dzkp*u%wtt)
      vc = indxm*incrzpp*(dxim*dzkpp*(v%t  + v%st  ) + &
                          dxi *dzkpp*(v%wt + v%wst ) + &
                          dxim*dzkp *(v%tt + v%stt ) + &
                          dxi *dzkp *(v%wtt+ v%wstt))*indy(j)
      wc = incrxm*(dxim*w%t+dxi*w%wt)*indzpp
      ft = ft - 0.125D0*(wc-ABS(wc))*(uc+ABS(uc))*deltar &
              + 0.25D0*onethird*(wc-ABS(wc))*(uc+ABS(uc))*ABS(vc)*deltar
!     
! wnt-wn
      deltar  = dens%wnt*w%wnt - dens%wn*w%wn
      IF (i <= 2) deltar = 0.D0
      IF (j >= ny-1) deltar = 0.D0
      uc = 0.5D0*(u%wnt+u%wwnt)*dt*indx(i-1)
      vc = 0.5D0*(v%wnt+w%wt)*dt*indy(j+1)
      wc = 0.5D0*(w%wnt+w%wn)
      fn = fn - 0.125D0*onethird*(wc-ABS(wc))*(vc-abs(vc))*ABS(uc+ABS(uc))*deltar
!     
! wn-wnb
      deltar  = dens%wn*w%wn - dens%wnb*w%wnb
      IF (i <= 2) deltar = 0.D0
      IF (j >= ny-1) deltar = 0.D0
      uc = 0.5D0*(w%wn+w%wwn)*dt*indx(i-1)
      vc = 0.5D0*(w%wn+w%w)*dt*indy(j+1)
      wc = 0.5D0*(w%wn+w%wnb)
      fn = fn - 0.125D0*onethird*(wc+ABS(wc))*(vc-abs(vc))*ABS(uc+ABS(uc))*deltar
!     
! ent-et
      deltar  = dens%ent*w%ent - dens%et*w%et
      IF (i >= nx-1) deltar = 0.D0
      IF (k >= nz-1) deltar = 0.D0
      uc = indyp*incrzpp*(dyjp*dzkpp*(u%t  + u%et  ) + &
                          dyj *dzkpp*(u%nt + u%ent ) + & 
                          dyjp*dzkp *(u%tt + u%ett ) + & 
                          dyj *dzkp *(u%ntt+ u%entt))*indx(i+1)
      vc = 0.5D0*indzpp*(dzkpp*v%et+dzkp*v%ett)
      wc = incryp*(dyjp*w%et+dyj*w%ent)*indzpp
      fe = fe - 0.125D0*onethird*(vc-ABS(vc))*(uc-ABS(uc))*ABS(wc-ABS(wc))*deltar
!     
! et-est
      deltar  = dens%et*w%et - dens%est*w%est
      IF (i >= nx-1) deltar = 0.D0
      IF (k >= nz-1) deltar = 0.D0
      uc = indym*incrzpp*(dyjm*dzkpp*(u%t  + u%et  ) + &
                          dyj *dzkpp*(u%st + u%est ) + & 
                          dyjm*dzkp*(u%tt + u%ett ) + & 
                          dyj *dzkp*(u%stt+ u%estt))*indx(i+1)
      vc = 0.5D0*indzpp*(dzkpp*v%est+dzkp*v%estt)
      wc = incrym*(dyjm*w%et+dyj*w%est)*indzpp
      fe = fe - 0.125D0*onethird*(vc+ABS(vc))*(uc-ABS(uc))*ABS(wc-ABS(wc))*deltar
!     
! enb-eb
      deltar  = dens%enb*w%enb - dens%eb*w%eb
      IF (i >= nx-1) deltar = 0.D0
      IF (k <= 2) deltar = 0.D0
      uc = indyp*incrzm*(dyjp*dzkm*(u%c  + u%e  ) + &
                         dyj *dzkm*(u%n  + u%en ) + & 
                         dyjp*dzk *(u%b  + u%eb ) + & 
                         dyj *dzk *(u%nb + u%enb))*indx(i+1)
      vc = 0.5D0*indzm*(dzkm*v%e+dzk*v%eb)
      wc = incryp*(dyjp*w%eb+dyj*w%enb)*indzm
      fe = fe - 0.125D0*onethird*(vc-ABS(vc))*(uc-ABS(uc))*ABS(wc+ABS(wc))*deltar
!     
! eb-esb
      deltar  = dens%eb*w%eb - dens%esb*w%esb
      IF (i >= nx-1) deltar = 0.D0
      IF (k <= 2) deltar = 0.D0
      uc = indym*incrzm*(dyjm*dzkm*(u%c  + u%e  ) + &
                         dyj *dzkm*(u%s  + u%es ) + & 
                         dyjm*dzk *(u%b  + u%eb ) + & 
                         dyj *dzk *(u%sb + u%esb))*indx(i+1)
      vc = 0.5D0*indzm*(dzkm*v%es+dzk*v%esb)
      wc = incrym*(dyjm*w%eb+dyj*w%esb)*indzm
      fe = fe - 0.125D0*onethird*(vc+ABS(vc))*(uc-ABS(uc))*ABS(wc+ABS(wc))*deltar
!     
! ent-en
      deltar  = dens%ent*w%ent - dens%en*w%en
      IF (i >= nx-1) deltar = 0.D0
      IF (j >= ny-1) deltar = 0.D0
      uc = 0.5D0*(u%ent+u%nt)*dt*indx(i+1)
      vc = 0.5D0*(v%ent+v%et)*dt*indy(j+1)
      wc = 0.5D0*(w%ent+w%en)
      fn = fn - 0.125D0*onethird*(wc-ABS(wc))*(vc-ABS(vc))*ABS(uc-ABS(uc))*deltar
!     
! en-enb
      deltar  = dens%en*w%en - dens%enb*w%enb
      IF (i >= nx-1) deltar = 0.D0
      IF (j >= ny-1) deltar = 0.D0
      uc = 0.5D0*(u%en+u%n)*dt*indx(i+1)
      vc = 0.5D0*(v%en+v%e)*dt*indy(j+1)
      wc = 0.5D0*(w%en+w%enb)
      fn = fn - 0.125D0*onethird*(wc+ABS(wc))*(vc-ABS(vc))*ABS(uc-ABS(uc))*deltar
!     
! ent-nt
      deltar  = dens%ent*w%ent - dens%nt*w%nt
      IF (k >= nz-1) deltar = 0.D0
      IF (j >= ny-1) deltar = 0.D0
      uc = 0.5D0*indzpp*(dzkpp*u%nt+dzkp*u%ntt)
      vc = indxp*incrzpp*(dxip*dzkpp*(v%t  + v%nt  ) + &
                          dxi *dzkpp*(v%et + v%ent ) + &
                          dxip*dzkp *(v%tt + v%ntt ) + &
                          dxi *dzkp *(v%ett+ v%entt))*indy(j+1)
      wc = incrxp*(dxip*w%nt+dxi*w%ent)*indzpp
      ft = ft - 0.125D0*onethird*(wc-ABS(wc))*(uc-ABS(uc))*ABS((vc-ABS(vc)))*deltar
!     
! nt-wnt
      deltar  = dens%nt*w%nt - dens%wnt*w%wnt
      IF (k >= nz-1) deltar = 0.D0
      IF (j >= ny-1) deltar = 0.D0
      uc = 0.5D0*indzpp*(dzkpp*u%wnt+dzkp*u%wntt)
      vc = indxm*incrzpp*(dxim*dzkp*(v%t + v%nt  ) + &
                          dxi *dzkp*(v%wt + v%wnt ) + &
                          dxim*dzk *(v%tt + v%ntt ) + &
                          dxi *dzk *(v%wtt+ v%wntt))*indy(j+1)
      wc = incrxm*(dxim*w%nt+dxi*w%wnt)*indzpp
      ft = ft - 0.125D0*onethird*(wc-ABS(wc))*(uc+ABS(uc))*ABS((vc-ABS(vc)))*deltar
!     
! est-st
      deltar  = dens%est*w%est - dens%st*w%st
      IF (k >= nz-1) deltar = 0.D0
      IF (j <= 2) deltar = 0.D0
      uc = 0.5D0*indzpp*(dzkpp*u%st+dzkp*u%stt)
      vc = indxp*incrzpp*(dxip*dzkpp*(v%st  + v%sst  ) + &
                          dxi *dzkpp*(v%est + v%esst ) + &
                          dxip*dzkp *(v%stt + v%sstt ) + &
                          dxi *dzkp *(v%estt+ v%esstt))*indy(j-1)
      wc = incrxp*(dxip*w%st+dxi*w%est)*indzpp
      ft = ft - 0.125D0*onethird*(wc-ABS(wc))*(uc-ABS(uc))*ABS((vc+ABS(vc)))*deltar
!     
! st-wst
      deltar  = dens%st*w%st - dens%wst*w%wst
      IF (k >= nz-1) deltar = 0.D0
      IF (j <= 2) deltar = 0.D0
      uc = 0.5D0*indzpp*(dzkpp*u%wst+dzkp*u%wstt)
      vc = indxm*incrzpp*(dxim*dzkpp*(v%st  + v%sst  ) + &
                          dxi *dzkpp*(v%wst + v%wsst ) + &
                          dxim*dzkp *(v%stt + v%sstt ) + &
                          dxi *dzkp *(v%wstt+ v%wsstt))*indy(j-1)
      wc = incrxm*(dxim*w%st+dxi*w%wst)*indzpp
      ft = ft - 0.125D0*onethird*(wc-ABS(wc))*(uc+ABS(uc))*ABS((vc+ABS(vc)))*deltar
!     
      RETURN
      END SUBROUTINE ctu1_flw_3d
!------------------------------------------------------------------------
      SUBROUTINE ctu2_flw_3d(fe, fn, ft, dens, u, v, w, i, j, k)
!
! ... Compute the second order Corner Transport Upwind correction (step 3 in LeVeque algorithm)
!
      USE dimensions
      USE domain_mapping, ONLY: myijk
      USE grid, ONLY: dx, dy, dz, indx, indy, indz, flag, fluid
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: stencil
      USE time_parameters, ONLY: dt
      USE set_indexes, ONLY: imjk, ijkm, ipjk, ijkp, ijpk, ijmk
      USE flux_limiters, ONLY: limiters, lv
!
      REAL*8, INTENT(INOUT) :: fe, fn, ft
      TYPE(stencil), INTENT(IN) :: dens, u, v, w
      INTEGER, INTENT(IN) :: i, j, k
!
      REAL*8 :: delta, indzp, indxp, indxm, indxpp, indyp, indym, indypp
      REAL*8 :: dzp, dxm, dxp, dxpp, dym, dyp, dypp
      INTEGER :: ip2, jp2, kp2
      REAL*8 :: gradc, grade, gradw, gradt, gradb, gradn, grads, deltar, uref, vref, wref
!
      ip2 = MIN( nx, i+2 )
      jp2 = MIN( ny, j+2 )
      kp2 = MIN( nz, k+2 )
!
      dzp  = dz(k)  + dz(k+1)
      dxm  = dx(i)  + dx(i-1)
      dxp  = dx(i)  + dx(i+1)
      dxpp = dx(i+1)+ dx(ip2)
      dym  = dy(j)  + dy(j-1)
      dyp  = dy(j)  + dy(j+1)
      dypp = dy(j+1)+ dy(jp2)
!
      indzp  = 2.D0/dzp
      indxm  = 2.D0/dxm
      indxp  = 2.D0/dxp
      indxpp = 2.D0/dxpp
      indym  = 2.D0/dym
      indyp  = 2.D0/dyp
      indypp = 2.D0/dypp
!
! e-c
      gradc = indxp  * (w%e *dens%e  - w%c*dens%c)
      grade = indxpp * (w%ee*dens%ee - w%e*dens%e)
      gradw = indxm  * (w%c *dens%c  - w%w*dens%w)
      lim  = 0.D0
      erre = 0.D0
      uref = 0.5D0*indzp*(dz(k+1)*u%c + dz(k)*u%t)
      IF (uref >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=nx-1 .AND. flag(ipjk)==fluid) CALL limiters(lv,lim,erre)
      uref = ABS(0.5D0*indzp*(dz(k+1)*u%c + dz(k)*u%t))
      deltar = w%e*dens%e - w%c*dens%c
      fe = fe + 0.5D0*uref*(1-dt*indxp*uref)*deltar*lim
!
! n-c
      gradc = indyp  * (w%n *dens%n  - w%c*dens%c)
      gradn = indypp * (w%nn*dens%nn - w%n*dens%n)
      grads = indym  * (w%c *dens%c  - w%s*dens%s)
      lim  = 0.D0
      erre = 0.D0
      vref = 0.5D0*indzp*(dz(k+1)*v%c + dz(k)*v%t)
      IF (vref >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = grads/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradn/gradc
      END IF
      IF (j/=ny-1 .AND. flag(ijpk)==fluid) CALL limiters(lv,lim,erre)
      vref = ABS(0.5D0*indzp*(dz(k+1)*v%c + dz(k)*v%t))
      deltar = w%n*dens%n - w%c*dens%c
      fn = fn + 0.5D0*vref*(1-dt*indyp*vref)*deltar*lim
!
! t-c
      gradc = indz(k+1) * (w%t *dens%t  - w%c*dens%c)
      gradt = indz(kp2) * (w%tt*dens%tt - w%t*dens%t)
      gradb = indz(k)   * (w%c *dens%c  - w%b*dens%b)
      lim  = 0.D0
      erre = 0.D0
      wref = 0.5D0*(w%t + w%c)
      IF (wref >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(ijkp)==fluid) CALL limiters(lv,lim,erre)
      wref = ABS(0.5D0*(w%c + w%t))
      deltar = w%t*dens%t - w%c*dens%c
      ft = ft + 0.5D0*wref*(1-dt*indz(k+1)*wref)*deltar*lim
!
      RETURN
      END SUBROUTINE ctu2_flw_3d
!----------------------------------------------------------------------
      SUBROUTINE ctu3_flw_3d(fe, fn, ft, dens, u, v, w, i, j, k)
!
! ... Compute the second order Corner Transport Upwind correction (step 4 in LeVeque algorithm)
!
      USE dimensions
      USE domain_mapping, ONLY: myijk
      USE grid, ONLY: dx, dy, dz, indx, indy, indz, fluid, flag
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: stencil
      USE set_indexes, ONLY: imjk, ijkm, ipjk, ijkp, ijmk, ijpk
      USE set_indexes, ONLY: ipjpk, ipjmk, imjpk, imjmk, ipjkp, ipjkm, imjkp, imjkm
      USE set_indexes, ONLY: ijpkp, ijpkm, ijmkp, ijmkm, ipjpkp
      USE set_indexes, ONLY: ipjpkm, ipjmkp, imjpkp, ipjmkm, imjpkm, imjmkp
      USE time_parameters, ONLY: dt
      USE flux_limiters, ONLY: limiters, lv

      IMPLICIT NONE

      REAL*8, INTENT(INOUT) :: fe, ft, fn
      TYPE(stencil), INTENT(IN) :: dens, u, v, w 
      INTEGER, INTENT(IN) :: i, j, k 

      INTEGER :: im2, ip2, jm2, jp2, km2, kp2
      REAL*8 :: uc, vc, wc, deltar, erre, lim, s

      REAL*8 :: dxm, dxp, dxmm, dxpp, indxpp, indxp, indxm, indxmm
      REAL*8 :: dym, dyp, dymm, dypp, indypp, indyp, indym, indymm
      REAL*8 :: dzm, dzp, dzmm, dzpp, indzpp, indzp, indzm, indzmm
      REAL*8 :: incrxp, incrxm, incryp, incrym, incrzp, incrzm, incrzpp
      REAL*8 :: gradc, grade, gradw, gradn, grads, gradt, gradb

      im2 = MAX(  1, i-2 )
      ip2 = MIN( nx, i+2 )
      jm2 = MAX(  1, j-2 )
      jp2 = MIN( ny, j+2 )
      km2 = MAX(  1, k-2 )
      kp2 = MIN( nz, k+2 )
!
      dxm  = dx(i)  + dx(i-1)
      dxp  = dx(i)  + dx(i+1)
      dxpp = dx(i+1)+ dx(ip2)
      dxmm = dx(i-1)+ dx(im2)
      dym  = dy(j)  + dy(j-1)
      dyp  = dy(j)  + dy(j+1)
      dypp = dy(j+1)+ dy(jp2)
      dymm = dy(j-1)+ dy(jm2)
      dzm  = dz(k)  + dz(k-1)
      dzp  = dz(k)  + dz(k+1)
      dzpp = dz(k+1)+ dz(kp2)
      dzmm = dz(k-1)+ dz(km2)
!
      indxm  = 2.D0/dxm
      indxp  = 2.D0/dxp
      indxpp = 2.D0/dxpp
      indxmm = 2.D0/dxmm
      indym  = 2.D0/dym
      indyp  = 2.D0/dyp
      indypp = 2.D0/dypp
      indymm = 2.D0/dymm
      indzm  = 2.D0/dzm
      indzp  = 2.D0/dzp
      indzpp = 2.D0/dzpp
      indzmm = 2.D0/dzmm
!
      incrxp  = 0.5D0*dt*indxp
      incrxm  = 0.5D0*dt*indxm
      incryp  = 0.5D0*dt*indyp
      incrym  = 0.5D0*dt*indym
      incrzp  = 0.125D0*dt*indzp
      incrzm  = 0.125D0*dt*indzm
      incrzpp = 0.125D0*dt*indxpp
!
! c-w
      deltar = dens%c*w%c - dens%w*w%w
      uc = 0.5D0*indzp*(dz(k+1)*u%w+dz(k)*u%wt)
      vc = indxm*incrzp*(dx(i-1)*dz(k+1)*(v%c + v%s  ) + &
                                dx(i)  *dz(k+1)*(v%w + v%ws ) + &
                                dx(i-1)*dz(k)  *(v%t + v%st ) + &
                                dx(i)  *dz(k)  *(v%wt+ v%wst))*indy(j)
      wc = incrxm*(dx(i-1)*w%c+dx(i)*w%w)*indzp
      gradc = indxm  * (w%c*dens%c - w%w *dens%w )
      grade = indxp  * (w%e*dens%e - w%c *dens%c )
      gradw = indxmm * (w%w*dens%w - w%ww*dens%ww)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=2 .AND. flag(imjk)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxm*ABS(uc))*deltar*lim
      fn = fn + 0.5D0*(vc+ABS(vc))*s
      ft = ft + 0.5D0*(1-ABS(vc))*(wc+ABS(wc))*s
!
! e-c
      deltar = dens%e*w%e - dens%c*w%c
      uc = 0.5D0*indzp*(dz(k+1)*u%c+dz(k)*u%t)
      vc = indxp*incrzp*(dx(i+1)*dz(k+1)*(v%c + v%s  ) + &
                                dx(i)  *dz(k+1)*(v%e + v%es ) + &
                                dx(i+1)*dz(k)  *(v%t + v%st ) + &
                                dx(i)  *dz(k)  *(v%et+ v%est))*indy(j)
      wc = incrxp*(dx(i+1)*w%c+dx(i)*w%e)*indzp
      gradc = indxp  * (w%e *dens%e  - w%c*dens%c)
      grade = indxpp * (w%ee*dens%ee - w%e*dens%e)
      gradw = indxm  * (w%c *dens%c  - w%w*dens%w)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=nx-1 .AND. flag(ipjk)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxp*ABS(uc))*deltar*lim
      fn = fn - 0.5D0*(vc+ABS(vc))*s
      ft = ft - 0.5D0*(1-ABS(vc))*(wc+ABS(wc))*s
!
! en-n
      deltar = dens%en*w%en - dens%n*w%n
      uc = 0.5D0*indzp*(dz(k+1)*u%n+dz(k)*u%nt)
      vc = indxp*incrzp*(dx(i+1)*dz(k+1)*(v%c + v%n  ) + &
                                dx(i)  *dz(k+1)*(v%e + v%en ) + &
                                dx(i+1)*dz(k)  *(v%t + v%nt ) + &
                                dx(i)  *dz(k)  *(v%et+ v%ent))*indy(j+1)
      wc = incrxp*(dx(i+1)*w%n+dx(i)*w%en)*indzp
      gradc = indxp  * (w%en *dens%en  - w%n *dens%n )
      grade = indxpp * (w%een*dens%een - w%en*dens%en)
      gradw = indxm  * (w%n  *dens%n   - w%wn*dens%wn)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=nx-1 .AND. flag(ipjpk)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxp*ABS(uc))*deltar*lim
      fn = fn - 0.5D0*(vc-ABS(vc))*s
      ft = ft - 0.25D0*ABS(vc-ABS(vc))*(wc+ABS(wc))*s
!
! n-wn
      deltar = dens%n*w%n - dens%wn*w%wn
      uc = 0.5D0*indzp*(dz(k+1)*u%wn+dz(k)*u%wnt)
      vc = indxm*incrzp*(dx(i-1)*dz(k+1)*(v%c + v%n  ) + &
                                dx(i)  *dz(k+1)*(v%w + v%wn ) + &
                                dx(i-1)*dz(k)  *(v%t + v%nt ) + &
                                dx(i)  *dz(k)  *(v%wt+ v%wnt))*indy(j+1)
      wc = incrxm*(dx(i-1)*w%n+dx(i)*w%wn)*indzp
      gradc = indxm  * (w%n *dens%n  - w%wn *dens%wn )
      grade = indxp  * (w%en*dens%en - w%n  *dens%n  )
      gradw = indxmm * (w%wn*dens%wn - w%wwn*dens%wwn)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=2 .AND. flag(imjpk)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxm*ABS(uc))*deltar*lim
      fn = fn + 0.5D0*(vc-ABS(vc))*s
      ft = ft + 0.25D0*ABS(vc-ABS(vc))*(wc+ABS(wc))*s
!
! es-s
      deltar = dens%es*w%es - dens%s*w%s
      uc = 0.5D0*indzp*(dz(k+1)*u%s+dz(k)*u%st)
      vc = indxp*incrzp*(dx(i+1)*dz(k+1)*(v%s  + v%ss  ) + &
                                dx(i)  *dz(k+1)*(v%es + v%ess ) + &
                                dx(i+1)*dz(k)  *(v%st + v%sst ) + &
                                dx(i)  *dz(k)  *(v%est+ v%esst))*indy(j-1)
      wc = incrxp*(dx(i+1)*w%s+dx(i)*w%es)*indzp
      gradc = indxp  * (w%es *dens%es  - w%s *dens%s )
      grade = indxpp * (w%ees*dens%ees - w%es*dens%es)
      gradw = indxm  * (w%s  *dens%s   - w%ws*dens%ws)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=nx-1 .AND. flag(ipjmk)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxp*ABS(uc))*deltar*lim
      ft = ft - 0.25D0*ABS(vc+ABS(vc))*(wc+ABS(wc))*s
!
! s-ws
      deltar = dens%s*w%s - dens%ws*w%ws
      uc = 0.5D0*indzp*(dz(k+1)*u%ws+dz(k)*u%wst)
      vc = indxm*incrzp*(dx(i-1)*dz(k+1)*(v%s + v%ss  ) + &
                                dx(i)  *dz(k+1)*(v%ws + v%wss ) + &
                                dx(i-1)*dz(k)  *(v%st + v%sst ) + &
                                dx(i)  *dz(k)  *(v%wst+ v%wsst))*indy(j-1)
      wc = incrxm*(dx(i-1)*w%s+dx(i)*w%ws)*indzp
      gradc = indxm  * (w%s *dens%s  - w%ws *dens%ws )
      grade = indxp  * (w%es*dens%es - w%s  *dens%s  )
      gradw = indxmm * (w%ws*dens%ws - w%wws*dens%wws)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=2 .AND. flag(imjmk)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxm*ABS(uc))*deltar*lim
      ft = ft + 0.25D0*ABS(vc+ABS(vc))*(wc+ABS(wc))*s
!
! et-t
      deltar = dens%et*w%et - dens%t*w%t
      uc = 0.5D0*indzpp*(dz(kp2)*u%t+dz(k+1)*u%tt)
      vc = indxp*incrzpp*(dx(i+1)*dz(kp2)*(v%t  + v%st  ) + &
                                 dx(i)  *dz(kp2)*(v%et + v%est ) + &
                                 dx(i+1)*dz(k+1)*(v%tt + v%stt ) + &
                                 dx(i)  *dz(k+1)*(v%ett+ v%estt))*indy(j)
      wc = incrxp*(dx(i+1)*w%t+dx(i)*w%et)*indzpp
      gradc = indxp  * (w%et *dens%et  - w%t *dens%t )
      grade = indxpp * (w%eet*dens%eet - w%et*dens%et)
      gradw = indxm  * (w%t  *dens%t   - w%wt*dens%wt)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=nx-1 .AND. flag(ipjkp)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxp*ABS(uc))*deltar*lim
      ft = ft - 0.5D0*(1-ABS(vc))*(wc-ABS(wc))*s
!
! t-wt
      deltar = dens%t*w%t - dens%wt*w%wt
      uc = 0.5D0*indzpp*(dz(kp2)*u%wt+dz(k+1)*u%wtt)
      vc = indxm*incrzpp*(dx(i-1)*dz(kp2)*(v%t  + v%st  ) + &
                                 dx(i)  *dz(kp2)*(v%wt + v%wst ) + &
                                 dx(i-1)*dz(k+1)*(v%tt + v%stt ) + &
                                 dx(i)  *dz(k+1)*(v%wtt+ v%wstt))*indy(j)
      wc = incrxm*(dx(i-1)*w%t+dx(i)*w%wt)*indzpp
      gradc = indxm  * (w%t *dens%t  - w%wt *dens%wt )
      grade = indxp  * (w%et*dens%et - w%t  *dens%t  )
      gradw = indxmm * (w%wt*dens%wt - w%wwt*dens%wwt)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=2 .AND. flag(imjkp)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxm*ABS(uc))*deltar*lim
      ft = ft + 0.5D0*(1-ABS(vc))*(wc-ABS(wc))*s
!
! ent-nt
      deltar = dens%ent*w%ent - dens%nt*w%nt
      uc = 0.5D0*indzpp*(dz(kp2)*u%nt+dz(k+1)*u%ntt)
      vc = indxp*incrzpp*(dx(i+1)*dz(kp2)*(v%t  + v%nt  ) + &
                                 dx(i)  *dz(kp2)*(v%et + v%ent ) + &
                                 dx(i+1)*dz(k+1)*(v%tt + v%ntt ) + &
                                 dx(i)  *dz(k+1)*(v%ett+ v%entt))*indy(j+1)
      wc = incrxp*(dx(i+1)*w%nt+dx(i)*w%ent)*indzpp
      gradc = indxp  * (w%ent *dens%ent  - w%nt *dens%nt )
      grade = indxpp * (w%eent*dens%eent - w%ent*dens%ent)
      gradw = indxm  * (w%nt  *dens%nt   - w%wnt*dens%wnt)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=nx-1 .AND. flag(ipjpkp)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxp*ABS(uc))*deltar*lim
      ft = ft - 0.25D0*ABS(vc-ABS(vc))*(wc-ABS(wc))*s
!
! nt-wnt
      deltar = dens%nt*w%nt - dens%wnt*w%wnt
      uc = 0.5D0*indzpp*(dz(kp2)*u%wnt+dz(k+1)*u%wntt)
      vc = indxm*incrzpp*(dx(i-1)*dz(k+1)*(v%t + v%nt  ) + &
                                 dx(i)  *dz(k+1)*(v%wt + v%wnt ) + &
                                 dx(i-1)*dz(k)  *(v%tt + v%ntt ) + &
                                 dx(i)  *dz(k)  *(v%wtt+ v%wntt))*indy(j+1)
      wc = incrxm*(dx(i-1)*w%nt+dx(i)*w%wnt)*indzpp
      gradc = indxm  * (w%nt *dens%nt  - w%wnt *dens%wnt )
      grade = indxp  * (w%ent*dens%ent - w%nt  *dens%nt  )
      gradw = indxmm * (w%wnt*dens%wnt - w%wwnt*dens%wwnt)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=2 .AND. flag(imjpkp)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxm*ABS(uc))*deltar*lim
      ft = ft + 0.25D0*ABS(vc-ABS(vc))*(wc-ABS(wc))*s
!
! est-st
      deltar = dens%est*w%est - dens%st*w%st
      uc = 0.5D0*indzpp*(dz(kp2)*u%st+dz(k+1)*u%stt)
      vc = indxp*incrzpp*(dx(i+1)*dz(kp2)*(v%st  + v%sst  ) + &
                                 dx(i)  *dz(kp2)*(v%est + v%esst ) + &
                                 dx(i+1)*dz(k+1)*(v%stt + v%sstt ) + &
                                 dx(i)  *dz(k+1)*(v%estt+ v%esstt))*indy(j-1)
      wc = incrxp*(dx(i+1)*w%st+dx(i)*w%est)*indzpp
      gradc = indxp  * (w%est *dens%est  - w%st *dens%st )
      grade = indxpp * (w%eest*dens%eest - w%est*dens%est)
      gradw = indxm  * (w%st  *dens%st   - w%wst*dens%wst)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=nx-1 .AND. flag(ipjmkp)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxp*ABS(uc))*deltar*lim
      ft = ft - 0.25D0*ABS(vc+ABS(vc))*(wc-ABS(wc))*s
!
! st-wst
      deltar = dens%st*w%st - dens%wst*w%wst
      uc = 0.5D0*indzpp*(dz(kp2)*u%wst+dz(k+1)*u%wstt)
      vc = indxm*incrzpp*(dx(i-1)*dz(kp2)*(v%st  + v%sst  ) + &
                                 dx(i)  *dz(kp2)*(v%wst + v%wsst ) + &
                                 dx(i-1)*dz(k+1)*(v%stt + v%sstt ) + &
                                 dx(i)  *dz(k+1)*(v%wstt+ v%wsstt))*indy(j-1)
      wc = incrxm*(dx(i-1)*w%st+dx(i)*w%wst)*indzpp
      gradc = indxm  * (w%st *dens%st  - w%wst *dens%wst )
      grade = indxp  * (w%est*dens%est - w%st  *dens%st  )
      gradw = indxmm * (w%wst*dens%wst - w%wwst*dens%wwst)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=2 .AND. flag(imjmkp)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxm*ABS(uc))*deltar*lim
      ft = ft + 0.25D0*ABS(vc+ABS(vc))*(wc-ABS(wc))*s
!
! c-s
      deltar = dens%c*w%c - dens%s*w%s
      uc = indym*incrzp*(dy(j-1)*dz(k+1)*(u%c  + u%w  ) + &
                                dy(j)  *dz(k+1)*(u%s  + u%ws ) + & 
                                dy(j-1)*dz(k)  *(u%t  + u%wt ) + & 
                                dy(j)  *dz(k)  *(u%st + u%wst))*indx(i)
      vc = 0.5D0*indzp*(dz(k+1)*v%s+dz(k)*v%st)
      wc = incrym*(dy(j-1)*w%c+dy(j)*w%s)*indzp
      gradc = indym  * (dens%c*w%c - dens%s *w%s )
      grads = indymm * (dens%s*w%s - dens%ss*w%ss)
      gradn = indyp  * (dens%n*w%n - dens%c *w%c )
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = grads / gradc
      ELSE
	IF (gradc /= 0.D0) erre = gradn / gradc
      ENDIF
      IF (j/=2 .AND. flag(ijmk)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indym*ABS(vc))*deltar*lim
      ft = ft + 0.5D0*(wc+ABS(wc))*s
      fe = fe + 0.5D0*(1-ABS(wc))*(uc+ABS(uc))*s
!
! n-c
      deltar = dens%n*w%n - dens%c*w%c
      uc = indyp*incrzp*(dy(j+1)*dz(k+1)*(u%c  + u%w  ) + &
                                dy(j)  *dz(k+1)*(u%n  + u%wn ) + & 
                                dy(j+1)*dz(k)  *(u%t  + u%wt ) + & 
                                dy(j)  *dz(k)  *(u%nt + u%wnt))*indx(i)
      vc = 0.5D0*indzp*(dz(k+1)*v%c+dz(k)*v%t)
      wc = incryp*(dy(j+1)*w%c+dy(j)*w%n)*indzp
      gradc = indyp  * (dens%n *w%n  - dens%c*w%c)
      grads = indym  * (dens%c *w%c  - dens%s*w%s)
      gradn = indypp * (dens%nn*w%nn - dens%n*w%n)
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = grads / gradc
      ELSE
	IF (gradc /= 0.D0) erre = gradn / gradc
      ENDIF
      IF (j/=ny-1 .AND. flag(ijpk)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indyp*ABS(vc))*deltar*lim
      ft = ft - 0.5D0*(wc+ABS(wc))*s
      fe = fe - 0.5D0*(1-ABS(wc))*(uc+ABS(uc))*s
!
! en-e
      deltar = dens%en*w%en - dens%e*w%e
      uc = indyp*incrzp*(dy(j+1)*dz(k+1)*(u%c  + u%e  ) + &
                                dy(j)  *dz(k+1)*(u%n  + u%en ) + & 
                                dy(j+1)*dz(k)  *(u%t  + u%et ) + & 
                                dy(j)  *dz(k)  *(u%nt + u%ent))*indx(i+1)
      vc = 0.5D0*indzp*(dz(k+1)*v%e+dz(k)*v%et)
      wc = incryp*(dy(j+1)*w%e+dy(j)*w%en)*indzp
      gradc = indyp  * (dens%en *w%en  - dens%e *w%e )
      grads = indym  * (dens%e  *w%e   - dens%es*w%es)
      gradn = indypp * (dens%enn*w%enn - dens%en*w%en)
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = grads / gradc
      ELSE
	IF (gradc /= 0.D0) erre = gradn / gradc
      ENDIF
      IF (j/=ny-1 .AND. flag(ipjpk)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indyp*ABS(vc))*deltar*lim
      fe = fe - 0.5D0*(1-ABS(wc))*(uc-ABS(uc))*s
!
! e-es
      deltar = dens%e*w%e - dens%es*w%es
      uc = indym*incrzp*(dy(j-1)*dz(k+1)*(u%c  + u%e  ) + &
                                dy(j)  *dz(k+1)*(u%s  + u%es ) + & 
                                dy(j-1)*dz(k)  *(u%t  + u%et ) + & 
                                dy(j)  *dz(k)  *(u%st + u%est))*indx(i+1)
      vc = 0.5D0*indzp*(dz(k+1)*v%es+dz(k)*v%est)
      wc = incrym*(dy(j-1)*w%e+dy(j)*w%es)*indzp
      gradc = indym  * (dens%e *w%e  - dens%es *w%es )
      grads = indymm * (dens%es*w%es - dens%ess*w%ess)
      gradn = indyp  * (dens%en*w%en - dens%e  *w%e  )
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = grads / gradc
      ELSE
	IF (gradc /= 0.D0) erre = gradn / gradc
      ENDIF
      IF (j/=2 .AND. flag(ipjmk)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indym*ABS(vc))*deltar*lim
      fe = fe + 0.5D0*(1-ABS(wc))*(uc-ABS(uc))*s
!
! nt-t
      deltar = dens%nt*w%nt - dens%t*w%t
      uc = indyp*incrzpp*(dy(j+1)*dz(kp2)*(u%t  + u%wt  ) + &
                                 dy(j)  *dz(kp2)*(u%nt + u%wnt ) + & 
                                 dy(j+1)*dz(k+1)*(u%tt + u%wtt ) + & 
                                 dy(j)  *dz(k+1)*(u%ntt+ u%wntt))*indx(i)
      vc = 0.5D0*indzpp*(dz(kp2)*v%t+dz(k+1)*v%tt)
      wc = incryp*(dy(j+1)*w%t+dy(j)*w%nt)*indzpp
      gradc = indyp  * (dens%nt *w%nt  - dens%t *w%t )
      grads = indym  * (dens%t  *w%t   - dens%st*w%st)
      gradn = indypp * (dens%nnt*w%nnt - dens%nt*w%nt)
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = grads / gradc
      ELSE
	IF (gradc /= 0.D0) erre = gradn / gradc
      ENDIF
      IF (j/=ny-1 .AND. flag(ijpkp)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indyp*ABS(vc))*deltar*lim
      ft = ft - 0.5D0*(wc-ABS(wc))*s
      fe = fe - 0.25D0*ABS(wc-ABS(wc))*(uc+ABS(uc))*s
!
! t-st
      deltar = dens%t*w%t - dens%st*w%st
      uc = indym*incrzpp*(dy(j-1)*dz(kp2)*(u%t  + u%wt  ) + &
                                 dy(j)  *dz(kp2)*(u%st + u%wst ) + & 
                                 dy(j-1)*dz(k+1)*(u%tt + u%wtt ) + & 
                                 dy(j)  *dz(k+1)*(u%stt+ u%wstt))*indx(i)
      vc = 0.5D0*indzpp*(dz(kp2)*v%st+dz(k+1)*v%stt)
      wc = incrym*(dy(j-1)*w%t+dy(j)*w%st)*indzpp
      gradc = indym  * (dens%t *w%t  - dens%st *w%st )
      grads = indymm * (dens%st*w%st - dens%sst*w%sst)
      gradn = indyp  * (dens%nt*w%nt - dens%t  *w%t  )
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = grads / gradc
      ELSE
	IF (gradc /= 0.D0) erre = gradn / gradc
      ENDIF
      IF (j/=2 .AND. flag(ijmkp)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indym*ABS(vc))*deltar*lim
      ft = ft + 0.5D0*(wc-ABS(wc))*s
      fe = fe + 0.25D0*ABS(wc-ABS(wc))*(uc+ABS(uc))*s
!
! nb-b
      deltar = dens%nb*w%nb - dens%b*w%b
      IF (k <= 2) deltar = 0.D0
      uc = indyp*incrzm*(dy(j+1)*dz(k-1)*(u%c  + u%w  ) + &
                                dy(j)  *dz(k-1)*(u%n  + u%wn ) + & 
                                dy(j+1)*dz(k)  *(u%b  + u%wb ) + & 
                                dy(j)  *dz(k)  *(u%nb + u%wnb))*indx(i)
      vc = 0.5D0*indzm*(dz(k-1)*v%c+dz(k)*v%b)
      wc = incryp*(dy(j+1)*w%b+dy(j)*w%nb)*indzm
      gradc = indyp  * (dens%nb *w%nb  - dens%b *w%b )
      grads = indym  * (dens%b  *w%b   - dens%sb*w%sb)
      gradn = indypp * (dens%nnb*w%nnb - dens%nb*w%nb)
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = grads / gradc
      ELSE
	IF (gradc /= 0.D0) erre = gradn / gradc
      ENDIF
      IF (j/=ny-1 .AND. flag(ijpkm)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indyp*ABS(vc))*deltar*lim
      fe = fe - 0.25D0*ABS(wc+ABS(wc))*(uc+ABS(uc))*s
!
! b-sb
      deltar = dens%b*w%b - dens%sb*w%sb
      IF (k <= 2) deltar = 0.D0
      uc = indym*incrzm*(dy(j-1)*dz(k-1)*(u%c  + u%w  ) + &
                                dy(j)  *dz(k-1)*(u%s  + u%ws ) + & 
                                dy(j-1)*dz(k)  *(u%b  + u%wb ) + & 
                                dy(j)  *dz(k)  *(u%sb + u%wsb))*indx(i)
      vc = 0.5D0*indzm*(dz(k-1)*v%s+dz(k)*v%sb)
      wc = incrym*(dy(j-1)*w%b+dy(j)*w%sb)*indzm
      gradc = indym  * (dens%b *w%b  - dens%sb *w%sb )
      grads = indymm * (dens%sb*w%sb - dens%ssb*w%ssb)
      gradn = indyp  * (dens%nb*w%nb - dens%b  *w%b  )
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = grads / gradc
      ELSE
	IF (gradc /= 0.D0) erre = gradn / gradc
      ENDIF
      IF (j/=2 .AND. flag(ijmkm)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indym*ABS(vc))*deltar*lim
      fe = fe + 0.25D0*ABS(wc+ABS(wc))*(uc+ABS(uc))*s
!
! ent-et
      deltar = dens%ent*w%ent - dens%et*w%et
      uc = indyp*incrzpp*(dy(j+1)*dz(kp2)*(u%t  + u%et  ) + &
                                 dy(j)  *dz(kp2)*(u%nt + u%ent ) + & 
                                 dy(j+1)*dz(k+1)*(u%tt + u%ett ) + & 
                                 dy(j)  *dz(k+1)*(u%ntt+ u%entt))*indx(i+1)
      vc = 0.5D0*indzpp*(dz(kp2)*v%et+dz(k+1)*v%ett)
      wc = incryp*(dy(j+1)*w%et+dy(j)*w%ent)*indzpp
      gradc = indyp  * (dens%ent *w%ent  - dens%et *w%et )
      grads = indym  * (dens%et  *w%et   - dens%est*w%est)
      gradn = indypp * (dens%ennt*w%ennt - dens%ent*w%ent)
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = grads / gradc
      ELSE
	IF (gradc /= 0.D0) erre = gradn / gradc
      ENDIF
      IF (j/=ny-1 .AND. flag(ipjpkp)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indyp*ABS(vc))*deltar*lim
      fe = fe - 0.25D0*ABS(wc-ABS(wc))*(uc-ABS(uc))*s
!
! et-est
      deltar = dens%et*w%et - dens%est*w%est
      uc = indym*incrzpp*(dy(j-1)*dz(kp2)*(u%t  + u%et  ) + &
                                 dy(j)  *dz(kp2)*(u%st + u%est ) + & 
                                 dy(j-1)*dz(k+1)*(u%tt + u%ett ) + & 
                                 dy(j)  *dz(k+1)*(u%stt+ u%estt))*indx(i+1)
      vc = 0.5D0*indzpp*(dz(kp2)*v%est+dz(k+1)*v%estt)
      wc = incrym*(dy(j-1)*w%et+dy(j)*w%est)*indzpp
      gradc = indym  * (dens%et *w%et  - dens%est *w%est )
      grads = indymm * (dens%est*w%est - dens%esst*w%esst)
      gradn = indyp  * (dens%ent*w%ent - dens%et  *w%et  )
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = grads / gradc
      ELSE
	IF (gradc /= 0.D0) erre = gradn / gradc
      ENDIF
      IF (j/=2 .AND. flag(ipjmkp)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indym*ABS(vc))*deltar*lim
      fe = fe + 0.25D0*ABS(wc-ABS(wc))*(uc-ABS(uc))*s
!
! enb-eb
      deltar = dens%enb*w%enb - dens%eb*w%eb
      IF (k <= 2) deltar = 0.D0
      uc = indyp*incrzm*(dy(j+1)*dz(k-1)*(u%c  + u%e  ) + &
                                dy(j)  *dz(k-1)*(u%n  + u%en ) + & 
                                dy(j+1)*dz(k)  *(u%b  + u%eb ) + & 
                                dy(j)  *dz(k)  *(u%nb + u%enb))*indx(i+1)
      vc = 0.5D0*indzm*(dz(k-1)*v%e+dz(k)*v%eb)
      wc = incryp*(dy(j+1)*w%eb+dy(j)*w%enb)*indzm
      gradc = indyp  * (dens%enb *w%enb  - dens%eb *w%eb )
      grads = indym  * (dens%eb  *w%eb   - dens%esb*w%esb)
      gradn = indypp * (dens%ennb*w%ennb - dens%enb*w%enb)
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = grads / gradc
      ELSE
	IF (gradc /= 0.D0) erre = gradn / gradc
      ENDIF
      IF (j/=ny-1 .AND. flag(ipjpkm)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indyp*ABS(vc))*deltar*lim
      fe = fe - 0.25D0*ABS(wc+ABS(wc))*(uc-ABS(uc))*s
!
! eb-esb
      deltar = dens%eb*w%eb - dens%esb*w%esb
      IF (k <= 2) deltar = 0.D0
      uc = indym*incrzm*(dy(j-1)*dz(k-1)*(u%c  + u%e  ) + &
                                dy(j)  *dz(k-1)*(u%s  + u%es ) + & 
                                dy(j-1)*dz(k)  *(u%b  + u%eb ) + & 
                                dy(j)  *dz(k)  *(u%sb + u%esb))*indx(i+1)
      vc = 0.5D0*indzm*(dz(k-1)*v%es+dz(k)*v%esb)
      wc = incrym*(dy(j-1)*w%eb+dy(j)*w%esb)*indzm
      gradc = indym  * (dens%eb *w%eb  - dens%esb *w%esb )
      grads = indymm * (dens%esb*w%esb - dens%essb*w%essb)
      gradn = indyp  * (dens%enb*w%enb - dens%eb  *w%eb  )
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = grads / gradc
      ELSE
	IF (gradc /= 0.D0) erre = gradn / gradc
      ENDIF
      IF (j/=2 .AND. flag(ipjmkm)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indym*ABS(vc))*deltar*lim
      fe = fe + 0.25D0*ABS(wc+ABS(wc))*(uc-ABS(uc))*s
!
! c-b
      deltar = dens%c*w%c - dens%b*w%b
      uc = 0.5D0*(u%c+u%w)*dt*indx(i)
      vc = 0.5D0*(v%c+v%s)*dt*indy(j)
      wc = 0.5D0*(w%c+w%b)
      gradc = indz(k)   * (dens%c *w%c - dens%b *w%b )
      gradt = indz(k+1) * (dens%t *w%t - dens%c *w%c )
      gradb = indz(k-1) * (dens%b *w%b - dens%bb*w%bb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=2 .AND. flag(ijkm)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indz(k)*ABS(wc))*deltar*lim
      fe = fe + 0.5D0*(uc+ABS(uc))*s
      fn = fn + 0.5D0*(1-ABS(uc))*(vc+ABS(vc))*s
!
! t-c
      deltar = dens%t*w%t - dens%c*w%c
      uc = 0.5D0*(u%t+u%wt)*dt*indx(i)
      vc = 0.5D0*(v%t+v%st)*dt*indy(j)
      wc = 0.5D0*(w%t+w%c)
      gradc = indz(k+1) * (dens%t *w%t  - dens%c*w%c)
      gradt = indz(kp2) * (dens%tt*w%tt - dens%t*w%t)
      gradb = indz(k)   * (dens%c *w%c  - dens%b*w%b)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(ijkp)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indz(k+1)*ABS(wc))*deltar*lim
      fe = fe - 0.5D0*(uc+ABS(uc))*s
      fn = fn - 0.5D0*(1-ABS(uc))*(vc+ABS(vc))*s
!
! wt-w
      deltar = dens%wt*w%wt - dens%w*w%w
      uc = 0.5D0*(u%wt+u%wwt)*dt*indx(i-1)
      vc = 0.5D0*(v%wt+v%wst)*dt*indy(j)
      wc = 0.5D0*(w%wt+w%w  )
      gradc = indz(k+1) * (dens%wt *w%t  - dens%w *w%w )
      gradt = indz(kp2) * (dens%wtt*w%tt - dens%wt*w%wt)
      gradb = indz(k)   * (dens%w  *w%w  - dens%wb*w%wb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(imjkp)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indz(k+1)*ABS(wc))*deltar*lim
      fn = fn - 0.25D0*ABS(uc+ABS(uc))*(vc+ABS(vc))*s
!
! w-wb
      deltar = dens%w*w%w - dens%wb*w%wb
      uc = 0.5D0*(u%w+u%ww)*dt*indx(i-1)
      vc = 0.5D0*(v%w+v%ws)*dt*indy(j)
      wc = 0.5D0*(w%w+w%wb)
      gradc = indz(k)   * (dens%w *w%w  - dens%wb *w%wb )
      gradt = indz(k+1) * (dens%wt*w%wt - dens%w  *w%w  )
      gradb = indz(k-1) * (dens%wb*w%wb - dens%wbb*w%wbb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=2 .AND. flag(imjkm)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indz(k)*ABS(wc))*deltar*lim
      fn = fn + 0.25D0*ABS(uc+ABS(uc))*(vc+ABS(vc))*s
!
! et-e
      deltar = dens%et*w%et - dens%e*w%e
      uc = 0.5D0*(u%et+u%t  )*dt*indx(i+1)
      vc = 0.5D0*(v%et+v%est)*dt*indy(j)
      wc = 0.5D0*(w%et+w%e  )
      gradc = indz(k+1) * (dens%et *w%et  - dens%e *w%e )
      gradt = indz(kp2) * (dens%ett*w%ett - dens%et*w%et)
      gradb = indz(k)   * (dens%e  *w%e   - dens%eb*w%eb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(ipjkp)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indz(k+1)*ABS(wc))*deltar*lim
      fe = fe - 0.5D0*(uc-ABS(uc))*s
      fn = fn - 0.25D0*ABS(uc-ABS(uc))*(vc+ABS(vc))*s
!
! e-eb
      deltar = dens%e*w%e - dens%eb*w%eb
      uc = 0.5D0*(u%e+v%c )*dt*indx(i+1)
      vc = 0.5D0*(v%e+v%es)*dt*indy(j)
      wc = 0.5D0*(w%e+w%eb)
      gradc = indz(k)   * (dens%e *w%e  - dens%eb *w%eb )
      gradt = indz(k+1) * (dens%et*w%et - dens%e  *w%e  )
      gradb = indz(k-1) * (dens%eb*w%eb - dens%ebb*w%ebb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=2 .AND. flag(ipjkm)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indz(k)*ABS(wc))*deltar*lim
      fe = fe + 0.5D0*(uc-ABS(uc))*s
      fn = fn + 0.25D0*ABS(uc-ABS(uc))*(vc+ABS(vc))*s
!
! nt-n
      deltar = dens%nt*w%nt - dens%n*w%n
      uc = 0.5D0*(u%nt+u%wnt)*dt*indx(i)
      vc = 0.5D0*(v%nt+v%t  )*dt*indy(j+1)
      wc = 0.5D0*(w%nt+w%n  )
      gradc = indz(k+1) * (dens%nt *w%nt  - dens%n *w%n )
      gradt = indz(kp2) * (dens%ntt*w%ntt - dens%nt*w%nt)
      gradb = indz(k)   * (dens%n  *w%n   - dens%nb*w%nb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(ijpkp)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indz(k+1)*ABS(wc))*deltar*lim
      fn = fn - 0.5D0*(1-ABS(uc))*(vc-ABS(vc))*s
!
! n-nb
      deltar = dens%n*w%n - dens%nb*w%nb
      uc = 0.5D0*(u%n+u%wn)*dt*indx(i)
      vc = 0.5D0*(v%n+v%c )*dt*indy(j+1)
      wc = 0.5D0*(w%n+w%nb)
      gradc = indz(k)   * (dens%n *w%n  - dens%nb *w%nb )
      gradt = indz(k+1) * (dens%nt*w%nt - dens%n  *w%n  )
      gradb = indz(k-1) * (dens%nb*w%nb - dens%nbb*w%nbb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=2 .AND. flag(ijpkm)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indz(k)*ABS(wc))*deltar*lim
      fn = fn + 0.5D0*(1-ABS(uc))*(vc-ABS(vc))*s
!
! wnt-wn
      deltar = dens%wnt*w%wnt - dens%wn*w%wn
      uc = 0.5D0*(u%wnt+u%wwnt)*dt*indx(i-1)
      vc = 0.5D0*(v%wnt+w%wt  )*dt*indy(j+1)
      wc = 0.5D0*(w%wnt+w%wn  )
      gradc = indz(k+1) * (dens%wnt *w%wnt  - dens%wn *w%wn )
      gradt = indz(kp2) * (dens%wntt*w%wntt - dens%wnt*w%wnt)
      gradb = indz(k)   * (dens%wn  *w%wn   - dens%wnb*w%wnb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(imjpkp)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indz(k+1)*ABS(wc))*deltar*lim
      fn = fn - 0.25D0*ABS(uc+ABS(uc))*(vc-ABS(vc))*s
!
! wn-wnb
      deltar = dens%wn*w%wn - dens%wnb*w%wnb
      uc = 0.5D0*(w%wn+w%wwn)*dt*indx(i-1)
      vc = 0.5D0*(w%wn+w%w  )*dt*indy(j+1)
      wc = 0.5D0*(w%wn+w%wnb)
      gradc = indz(k)   * (dens%wn *w%wn  - dens%wnb *w%wnb )
      gradt = indz(k+1) * (dens%wnt*w%wnt - dens%wn  *w%wn  )
      gradb = indz(k-1) * (dens%wnb*w%wnb - dens%wnbb*w%wnbb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=2 .AND. flag(imjpkm)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indz(k)*ABS(wc))*deltar*lim
      fn = fn + 0.25D0*ABS(uc+ABS(uc))*(vc-ABS(vc))*s
!
! ent-en
      deltar = dens%ent*w%ent - dens%en*w%en
      uc = 0.5D0*(u%ent+u%nt)*dt*indx(i+1)
      vc = 0.5D0*(v%ent+v%et)*dt*indy(j+1)
      wc = 0.5D0*(w%ent+w%en)
      gradc = indz(k+1) * (dens%ent *w%ent  - dens%en *w%en )
      gradt = indz(kp2) * (dens%entt*w%entt - dens%ent*w%ent)
      gradb = indz(k)   * (dens%en  *w%en   - dens%enb*w%enb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(ipjpkp)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indz(k+1)*ABS(wc))*deltar*lim
      fn = fn - 0.25D0*ABS(uc-ABS(uc))*(vc-ABS(vc))*s
!
! en-enb 
      deltar = dens%en*w%en - dens%enb*w%enb
      uc = 0.5D0*(u%en+u%n  )*dt*indx(i+1)
      vc = 0.5D0*(v%en+v%e  )*dt*indy(j+1)
      wc = 0.5D0*(w%en+w%enb)
      gradc = indz(k)   * (dens%en *w%en  - dens%enb *w%enb )
      gradt = indz(k+1) * (dens%ent*w%ent - dens%en  *w%en  )
      gradb = indz(k-1) * (dens%enb*w%enb - dens%enbb*w%enbb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=2 .AND. flag(ipjpkm)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indz(k)*ABS(wc))*deltar*lim
      fn = fn + 0.25D0*ABS(uc-ABS(uc))*(vc-ABS(vc))*s
!
      RETURN
      END SUBROUTINE ctu3_flw_3d
!----------------------------------------------------------------------
      END MODULE convective_fluxes_w
!-----------------------------------------------------------------------
