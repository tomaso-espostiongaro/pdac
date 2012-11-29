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
        MODULE PROCEDURE ctu1_flw_2d
      END INTERFACE
      INTERFACE ctu2_flw
        MODULE PROCEDURE ctu2_flw_2d
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
      END MODULE convective_fluxes_w
!-----------------------------------------------------------------------
