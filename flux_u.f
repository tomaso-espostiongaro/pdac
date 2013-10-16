!----------------------------------------------------------------------
      MODULE convective_fluxes_u
!
! ... This module computes convective fluxes of _____x(r)-momentum ____
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
      REAL*8, PRIVATE :: fou                     ! first order upwind  !
      REAL*8, PRIVATE :: upwnd                   ! upwinded variable   !
      REAL*8, PRIVATE :: lim                     ! limiter             !
      REAL*8, PRIVATE :: erre                    ! gradients ratio     !
      REAL*8, PRIVATE :: incr                    ! increment           !
!
      INTERFACE flu
        MODULE PROCEDURE flu_2d, flu_3d
      END INTERFACE
      INTERFACE muscl_flu
        MODULE PROCEDURE muscl_flu_2d, muscl_flu_3d
      END INTERFACE
      INTERFACE ctu1_flu
        MODULE PROCEDURE ctu1_flu_2d, ctu1_flu_3d
      END INTERFACE
      INTERFACE ctu2_flu
        MODULE PROCEDURE ctu2_flu_2d, ctu2_flu_3d
      END INTERFACE
      INTERFACE ctu3_flu
        MODULE PROCEDURE ctu3_flu_2d, ctu3_flu_3d
      END INTERFACE
      
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE flu_3d(fe, fn, ft, fw, fs, fb, dens, u, v, w, i )
!
! ... Compute the convective fluxes on East, North, and Top sides of the cell
! ... for the momentum density along x.
! ... The stencils must be staggered
! ... Uses First Order Donor-cell tecnique
!
      USE dimensions
      USE domain_mapping, ONLY: myijk, meshinds
      USE grid, ONLY: dx, flag
      USE set_indexes, ONLY: imjk, ijmk, ijkm
      USE set_indexes, ONLY: stencil

      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: fe, fn, ft, fw, fs, fb
      TYPE(stencil), INTENT(IN) :: dens, u, v, w
      INTEGER, INTENT(IN) :: i
!
      REAL*8 :: dxp, indxp
!
      dxp  = dx(i) + dx(i+1)
      indxp  = 1.D0 / dxp
!
! ... on West volume bondary
!
      IF( .NOT.BTEST(flag(imjk),0) ) THEN
        cs = 0.5D0 * ( u%c + u%w )
        fw = 0.5D0*(cs-ABS(cs))*dens%c*u%c + 0.5D0*(cs+ABS(cs))*dens%w*u%w
      END IF
!
! ... on South volume bondary
!
      IF( .NOT.BTEST(flag(ijmk),0) ) THEN
        cs = ( dx(i+1) * v%s + dx(i) * v%es ) * indxp
        fs = 0.5D0*(cs-ABS(cs))*dens%c*u%c + 0.5D0*(cs+ABS(cs))*dens%s*u%s
      END IF
!
! ... on Bottom volume bondary
!
      IF( .NOT.BTEST(flag(ijkm),0) ) THEN
        cs = ( dx(i+1) * w%b + dx(i) * w%eb ) * indxp
        fb = 0.5D0*(cs-ABS(cs))*dens%c*u%c + 0.5D0*(cs+ABS(cs))*dens%b*u%b
      END IF
!
! ... on East volume boundary
!
      cs = 0.5D0 * (u%c + u%e)    
      fe = 0.5D0*(cs-ABS(cs))*dens%e*u%e + 0.5D0*(cs+ABS(cs))*dens%c*u%c
!
! ... on North volume boundary
!
      cs = (dx(i+1) * v%c + dx(i) * v%e) * indxp
      fn = 0.5D0*(cs-ABS(cs))*dens%n*u%n + 0.5D0*(cs+ABS(cs))*dens%c*u%c
     
!
! ... on Top volume boundary
!
      cs = (dx(i+1) * w%c + dx(i) * w%e) * indxp
      ft = 0.5D0*(cs-ABS(cs))*dens%t*u%t + 0.5D0*(cs+ABS(cs))*dens%c*u%c
!
      RETURN
      END SUBROUTINE flu_3d
!----------------------------------------------------------------------
      SUBROUTINE muscl_flu_3d(fe, fn, ft, dens, u, v, w, i, j, k)
!
! ... Compute the convective fluxes on East, North, and Top sides of the cell
! ... for the momentum density along x.
! ... The stencils must be staggered
!
      USE dimensions
      USE flux_limiters, ONLY: limiters, lv
      USE grid, ONLY: dx, dy, dz, indx, flag
      USE set_indexes, ONLY: stencil
      USE time_parameters, ONLY: dt

      IMPLICIT NONE
!
      REAL*8, INTENT(INOUT) :: fe, fn, ft
      TYPE(stencil), INTENT(IN) :: dens, u, v, w
      INTEGER, INTENT(IN) :: i, j, k
!
      REAL*8 :: dxp, indxp
      REAL*8 :: dym, dyp, dypp, indypp, indyp, indym
      REAL*8 :: dzp, indzp, dzm, indzm, dzpp, indzpp
      REAL*8 :: gradc, grade, gradw, gradn, grads, gradt, gradb

      INTEGER :: ip2, jp2, kp2
!
      ip2 = MIN( nx, i+2 )
      jp2 = MIN( ny, j+2 )
      kp2 = MIN( nz, k+2 )
!
      dxp  = dx(i)  + dx(i+1)
      dym  = dy(j)  + dy(j-1)
      dyp  = dy(j)  + dy(j+1)
      dypp = dy(j+1)+ dy(jp2)
      dzm  = dz(k)  + dz(k-1)
      dzp  = dz(k)  + dz(k+1)
      dzpp = dz(k+1)+ dz(kp2)

      indxp  = 1.D0/dxp
      indym  = 1.D0/dym
      indyp  = 1.D0/dyp
      indypp = 1.D0/dypp
      indzm  = 1.D0/dzm
      indzp  = 1.D0/dzp
      indzpp = 1.D0/dzpp
!
! ... MUSCL reconstruction of momentum
!
! ... on East volume boundary
!
      gradw = indx(i)   * (u%c *dens%c  - u%w*dens%w)
      gradc = indx(i+1) * (u%e *dens%e  - u%c*dens%c)
      grade = indx(ip2) * (u%ee*dens%ee - u%e*dens%e)
!
      lim = 0.D0
      erre = 0.D0
!
      cs = 0.5D0 * (u%c + u%e)
      !cn = cs * dt * indx(i+1)
      IF ( cs >= 0.D0 ) THEN
        IF (gradc /= 0.0) erre = gradw / gradc    
        incr = 0.5D0 * dx(i+1)
      ELSE IF ( cs < 0.D0 ) THEN
        IF (gradc /= 0.0) erre = grade / gradc    
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
      grads = (u%c *dens%c  - u%s*dens%s) * 2.D0 * indym
      gradc = (u%n *dens%n  - u%c*dens%c) * 2.D0 * indyp
      gradn = (u%nn*dens%nn - u%n*dens%n) * 2.D0 * indypp
!
      lim = 0.D0
      erre = 0.D0
!
      cs = (dx(i+1) * v%c + dx(i) * v%e) * indxp
      !cn = cs * dt * 2.D0 * indyp
      IF (cs >= 0.D0) THEN
        IF (gradc /= 0.0) erre = grads / gradc    
        incr = 0.5D0 * dy(j)
      ELSE IF (cs < 0.D0) THEN
        IF (gradc /= 0.0) erre = gradn / gradc    
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
      gradb = (u%c *dens%c  - u%b*dens%b) * 2.D0 * indzm
      gradc = (u%t *dens%t  - u%c*dens%c) * 2.D0 * indzp
      gradt = (u%tt*dens%tt - u%t*dens%t) * 2.D0 * indzpp
!
      lim = 0.D0
      erre = 0.D0
!
      cs = (dx(i+1) * w%c + dx(i) * w%e) * indxp
      !cn = cs * dt * 2.D0 * indzp
      IF (cs >= 0.D0) THEN
        IF (gradc /= 0.0) erre = gradb / gradc    
        incr = 0.5D0 * dz(k)
      ELSE IF (cs < 0.D0) THEN
        IF (gradc /= 0.0) erre = gradt / gradc    
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
      END SUBROUTINE muscl_flu_3d
!----------------------------------------------------------------------
      SUBROUTINE flu_2d(fe, ft, fw, fb, dens, u, w, i)
!
! ... Compute the convective fluxes on East, Top sides of the cell
! ... for the momentum density along r.
! ... Uses First Order Donor-cell tecnique
!
      USE dimensions
      USE grid, ONLY: flag
      USE grid, ONLY: dx, r
      USE set_indexes, ONLY: imjk, ijkm
      USE set_indexes, ONLY: stencil

      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: fe, ft, fw, fb
      TYPE(stencil), INTENT(IN) :: dens, u, w
      INTEGER, INTENT(IN) :: i
!
      REAL*8 :: dxp, indxp 
!
      dxp = dx(i) + dx(i+1)
      indxp = 1.D0/dxp
!
! ... on West volume bondary
!
      IF( .NOT.BTEST(flag(imjk),0) ) THEN
        cs = 0.5D0*(u%c + u%w)
        fw = 0.5D0*(cs-ABS(cs))*dens%c*u%c*r(i) + 0.5D0*(cs+ABS(cs))*dens%w*u%w*r(i)
      END IF
!
! ... on Bottom volume bondary
!
      IF( .NOT.BTEST(flag(ijkm),0) ) THEN
        cs = (dx(i+1) * w%b + dx(i) * w%eb) * indxp
        fb = 0.5D0*(cs-ABS(cs))*dens%c*u%c + 0.5D0*(cs+ABS(cs))*dens%b*u%b
      END IF
!
! ... on East volume boundary
!
      cs = 0.5D0 * (u%c + u%e)
      fe = 0.5D0*(cs-ABS(cs))*dens%e*u%e*r(i+1) + 0.5D0*(cs+ABS(cs))*dens%c*u%c*r(i+1)
!
! ... on Top volume boundary
!
      cs = (dx(i+1) * w%c + dx(i) * w%e) * indxp
      ft = 0.5D0*(cs-ABS(cs))*dens%t*u%t + 0.5D0*(cs+ABS(cs))*dens%c*u%c
!
      RETURN
      END SUBROUTINE flu_2d
!----------------------------------------------------------------------
      SUBROUTINE muscl_flu_2d(fe, ft, dens, u, w, i, k)
!
! ... Compute the convective fluxes on East, Top sides of the cell
! ... for the momentum density along r.
! ... Stencils must be staggered
!
      USE dimensions
      USE flux_limiters, ONLY: limiters, lv
      USE grid, ONLY: flag
      USE grid, ONLY: dx, indx, r, dz
      USE set_indexes, ONLY: stencil
      USE time_parameters, ONLY: dt

      IMPLICIT NONE
!
      REAL*8, INTENT(INOUT) :: fe, ft
      TYPE(stencil), INTENT(IN) :: dens, u, w
      INTEGER, INTENT(IN) :: i, k

      REAL*8 :: dxp, indxp
      REAL*8 :: dzp, indzp, dzm, indzm, dzpp, indzpp
      REAL*8 :: gradc, grade, gradw, gradb, gradt

      INTEGER :: ip2, kp2
!
      ip2 = MIN( nx, i+2 )
      kp2 = MIN( nz, k+2 )
!
      dxp  = dx(i)  + dx(i+1)
      dzm  = dz(k)  + dz(k-1)
      dzp  = dz(k)  + dz(k+1)
      dzpp = dz(k+1)+ dz(kp2)

      indxp  = 1.D0/dxp
      indzm  = 1.D0/dzm
      indzp  = 1.D0/dzp
      indzpp = 1.D0/dzpp
!
! ... MUSCL reconstruction of momentum
!
! ... on East volume boundary
!
      gradc = indx(i+1) * (dens%e *u%e  - dens%c*u%c)
      gradw = indx(i)   * (dens%c *u%c  - dens%w*u%w)
      grade = indx(ip2) * (dens%ee*u%ee - dens%e*u%e)
!
      lim = 0.D0
      erre = 0.D0
!
      cs = 0.5D0 * (u%c + u%e)
      !cn = cs * dt * indx(i+1)
      IF ( cs >= 0.D0 ) THEN
        IF( gradc /= 0.D0) erre = gradw / gradc
        incr = 0.5D0 * dx(i+1)
      ELSE IF ( cs < 0.D0 ) THEN
        IF( gradc /= 0.D0) erre = grade / gradc
        incr = - 0.5D0 * dx(i+1)
      END IF
!
      CALL limiters(lv,lim,erre)
!
      upwnd = lim * gradc * incr
!
      fe = fe + upwnd * cs * r(i+1)
!
! ... on Top volume boundary
!
      gradc = (dens%t *u%t  - dens%c*u%c) * 2.D0 * indzp
      gradb = (dens%c *u%c  - dens%b*u%b) * 2.D0 * indzm
      gradt = (dens%tt*u%tt - dens%t*u%t) * 2.D0 * indzpp
!
      lim = 0.D0
      erre = 0.D0
!
      cs = (dx(i+1) * w%c + dx(i) * w%e) * indxp
      !cn = cs * dt * 2.D0 * indzp
      IF (cs >= 0.D0) THEN
        IF( gradc /= 0.D0) erre = gradb / gradc 
	incr = 0.5D0 * dz(k)
      ELSE IF (cs < 0.D0) THEN
        IF( gradc /= 0.D0) erre = gradt / gradc 
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
      END SUBROUTINE muscl_flu_2d
!-----------------------------------------------------------------------
      SUBROUTINE ctu1_flu_2d(fe, ft, dens, u, w, i, k)
!
! ... Compute the first order Corner Transport Upwind correction (step 2 in LeVeque algorithm)
!
      USE dimensions
      USE domain_mapping, ONLY: myijk
      USE grid, ONLY: dx, dz, indx, r
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: stencil
      USE time_parameters, ONLY: dt, time
      USE set_indexes, ONLY: imjk, ijkm
!
      REAL*8, INTENT(INOUT) :: fe, ft
      TYPE(stencil), INTENT(IN) :: dens, u, w
      INTEGER, INTENT(IN) :: i, k

      REAL*8 :: uc, wc, deltar
      REAL*8 :: dzm, dzp, dxm, dxp, dxpp, indzm, indzp, indxp, indxpp
      REAL*8 :: incrxp, incrxpp, incrzp, incrzm
      REAL*8 :: dxi, dxip, dxipp, dxim, dzk, dzkp, dzkm
      INTEGER ip2

      ip2 = MIN(nx, i+2)
!
      dxi   = dx(i)
      dxip  = dx(i+1)
      dxipp = dx(ip2)
      dxim  = dx(i-1)
      dzk   = dz(k)
      dzkp  = dz(k+1)
      dzkm  = dz(k-1)
!
      dxm  = dxi + dxim
      dxp  = dxi + dxip
      dxpp = dxip + dxipp
      dzm  = dzk + dzkm
      dzp  = dzk + dzkp
!
      indxp  = 2.D0/dxp 
      indxpp = 2.D0/dxpp 
      indzm  = 2.D0/dzm
      indzp  = 2.D0/dzp
!
      incrxp  = 0.125D0*dt*indx(i) 
      incrxpp = 0.125D0*dt*indx(i+1)
      incrzp  = 0.125D0*dt*indzp
      incrzm  = 0.125D0*dt*indzm
!
! c-w
      uc = 0.5D0*(u%c+u%w)
      wc = 0.5D0*(w%c+w%b)
      deltar = dens%c*u%c - dens%w*u%w
      ft = ft - incrxp*(uc+ABS(uc))*(wc+ABS(wc))*deltar
!
! e-c
      uc = 0.5D0*(u%c+u%e)
      wc = 0.5D0*(w%e+w%eb)
      deltar = dens%e*u%e - dens%c*u%c
      ft = ft - incrxpp*(uc-ABS(uc))*(wc+ABS(wc))*deltar
!
! t-wt
      uc = 0.5D0*(u%t+u%wt)
      wc = 0.5D0*(w%t+w%c)
      deltar = dens%t*u%t - dens%wt*u%wt
      IF (k >= nz-1) deltar = 0.D0
      ft = ft - incrxp*(uc+ABS(uc))*(wc-ABS(wc))*deltar
!
! et-t
      uc = 0.5D0*(u%t+u%et)
      wc = 0.5D0*(w%et+w%e)
      deltar = dens%et*u%et - dens%t*u%t
      IF (k >= nz-1) deltar = 0.D0
      ft = ft - incrxpp*(uc-ABS(uc))*(wc-ABS(wc))*deltar
!
! c-b
      uc = 0.5D0*indzm*(dzk*u%b +dzkm*u%c)
      wc = 0.5D0*indxp*(dxip*w%b+dxi*w%eb)
      deltar = dens%c*u%c - dens%b*u%b
      fe = fe - incrzm*(uc+ABS(uc))*(wc+ABS(wc))*deltar*r(i+1)
!
! t-c
      uc = 0.5D0*indzp*(dzkp*u%c+dzk*u%t)
      wc = 0.5D0*indxp*(dxip*w%c+dxi*w%e)
      deltar = dens%t*u%t - dens%c*u%c
      fe = fe - incrzp*(uc+ABS(uc))*(wc-ABS(wc))*deltar*r(i+1)
!
! e-eb
      uc = 0.5D0*indzm *(dzk*u%eb  +dzkm*u%e)
      wc = 0.5D0*indxpp*(dxipp*w%eb+dxip*w%eeb)
      deltar = dens%e*u%e-dens%eb*u%eb
      IF (i >= nx-1) deltar = 0.D0
      fe = fe - incrzm*(uc-ABS(uc))*(wc+ABS(wc))*deltar*r(i+1)
!
! et-e
      uc = 0.5D0*indzp*(dzk*u%et  +dzkp*u%e)
      wc = 0.5D0*indxpp*(dxipp*w%e+dxip*w%ee)
      deltar = dens%et*u%et-dens%e*u%e
      IF (i >= nx-1) deltar = 0.D0
      fe = fe - incrzp*(uc-ABS(uc))*(wc-ABS(wc))*deltar*r(i+1)
!
      RETURN
      END SUBROUTINE ctu1_flu_2d
!----------------------------------------------------------------------
      SUBROUTINE ctu2_flu_2d(fe, ft, dens, u, w, i, k)
!
! ... Compute the second order Corner Transport Upwind correction (step 3 in LeVeque algorithm)
!
      USE dimensions
      USE domain_mapping, ONLY: myijk
      USE grid, ONLY: dx, dz, indx, indz, flag, fluid, r
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: stencil
      USE time_parameters, ONLY: dt
      USE set_indexes, ONLY: imjk, ijkm, ijkp, ipjk
      USE flux_limiters, ONLY: limiters, lv

      REAL*8, INTENT(INOUT) :: fe, ft
      TYPE(stencil), INTENT(IN) :: dens, u, w
      INTEGER, INTENT(IN) :: i, k

      REAL*8 :: deltar, uref, wref, lim, erre
      REAL*8 :: dxp, dzm, dzp, dzpp, indzm, indzp, indzpp, indxp
      INTEGER :: ip2, kp2
      REAL*8 :: gradc, grade, gradw, gradt, gradb
!
      ip2 = MIN( nx, i+2 )
      kp2 = MIN( nz, k+2 )
!
      dxp  = dx(i)  + dx(i+1)
      dzm  = dz(k)  + dz(k-1)
      dzp  = dz(k)  + dz(k+1)
      dzpp = dz(k+1)+ dz(kp2)
!
      indxp  = 2.D0/dxp
      indzm  = 2.D0/dzm
      indzp  = 2.D0/dzp
      indzpp = 2.D0/dzpp
!
! e-c
      gradc = indx(i+1)*(dens%e *u%e  - dens%c*u%c)
      grade = indx(ip2)*(dens%ee*u%ee - dens%e*u%e)
      gradw = indx(i)  *(dens%c *u%c  - dens%w*u%w)
      lim  = 0.D0
      erre = 0.D0
      uref = 0.5D0*(u%c + u%e)
      IF (uref >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=nx-1 .AND. flag(ipjk)==fluid) CALL limiters(lv,lim,erre)
      uref = ABS(0.5D0*(u%c + u%e))
      deltar = dens%e*u%e - dens%c*u%c
      fe = fe + 0.5D0*uref*(1-dt*indx(i+1)*uref)*deltar*lim*r(i+1)
!
! t-c
      gradc = indzp *(dens%t *u%t   - dens%c*u%c)
      gradt = indzpp*(dens%tt*u%tt  - dens%t*u%t)
      gradb = indzm *(dens%c *u%c   - dens%b*u%b)
      lim  = 0.D0
      erre = 0.D0
      wref = 0.5D0*indxp*(dx(i)*w%e + dx(i+1)*w%c)
      IF (wref >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(ijkp)==fluid) CALL limiters(lv,lim,erre)
      wref = ABS(0.5D0*indxp*(dx(i)*w%e + dx(i+1)*w%c))
      deltar = dens%t*u%t - dens%c*u%c
      ft = ft + 0.5D0*wref*(1-dt*indzp*wref)*deltar*lim
!
      RETURN
      END SUBROUTINE ctu2_flu_2d
!-----------------------------------------------------------------------
      SUBROUTINE ctu3_flu_2d(fe, ft, dens, u, w, i, k)
!
! ... Compute the second order Corner Transport Upwind correction (step 4 in LeVeque algorithm)
!
      USE dimensions
      USE domain_mapping, ONLY: myijk
      USE grid, ONLY: dx, dz, indx, indz, fluid, flag, r
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: stencil, ipjk, imjk, ijkp, ijkm
      USE set_indexes, ONLY: ipjkm, ipjkp, imjkp, imjkm, ipjpk, imjpk, ipjmk, imjmk
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
      REAL*8 :: indxp, indzp, indzm, indxpp, indzpp, indzmm
      REAL*8 :: incrzp, incrzpp, incrxp, incrxpp
      REAL*8 :: dxp, dxpp, dzm, dzp, dzmm, dzpp
      INTEGER :: ip2, km2, kp2

      ip2 = MIN( nx, i+2 )
      km2 = MAX(  1, k-2 )
      kp2 = MIN( nz, k+2 )
!
      dxp  = dx(i)  + dx(i+1)
      dxpp = dx(i+1)+ dx(ip2) 
      dzm  = dz(k)  + dz(k-1)
      dzp  = dz(k)  + dz(k+1)
      dzpp = dz(k+1)+ dz(kp2)
      dzmm = dz(k-1)+ dz(km2)
!
      indxp  = 2.D0/dxp
      indxpp = 2.D0/dxpp
      indzm  = 2.D0/dzm
      indzp  = 2.D0/dzp
      indzpp = 2.D0/dzpp
      indzmm = 2.D0/dzmm
!
      incrzp  = 0.5D0*dt*indz(k)
      incrzpp = 0.5D0*dt*indz(k+1)
      incrxp  = 0.5D0*dt*indxp
      incrxpp = 0.5D0*dt*indxpp
!
! c-w
      uc = 0.5D0*(u%c + u%w)
      wc = 0.5D0*(w%c + w%b)
      deltar = 0.5D0*ABS(uc)*(1 - dt*indx(i)*ABS(uc))*(u%c*dens%c - u%w*dens%w)
      gradc = indx(i)  *(u%c*dens%c - u%w *dens%w )
      grade = indx(i+1)*(u%e*dens%e - u%c *dens%c )
      gradw = indx(i-1)*(u%w*dens%w - u%ww*dens%ww)
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
      uc = 0.5D0*(u%c + u%e)
      wc = 0.5D0*(w%e + w%eb)
      deltar = 0.5D0*ABS(uc)*(1 - dt*indx(i+1)*ABS(uc))*(u%e*dens%e - u%c*dens%c)
      gradc = indx(i+1)*(u%e *dens%e  - u%c*dens%c)
      grade = indx(ip2)*(u%ee*dens%ee - u%e*dens%e)
      gradw = indx(i)  *(u%c *dens%c  - u%w*dens%w)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=nx-1 .AND. flag(ipjk)==fluid) CALL limiters(lv,lim,erre)
      !IF (wc >= 0.D0) THEN
      !  ft = ft - dt*indz(k)*wc*deltar*lim
      !!ELSE
      !!  fb = fb - dt*indz(k)*wc*deltar*lim
      !END IF
      ft = ft - incrzp*(wc+ABS(wc))*deltar*lim
!
! t-wt
      uc = 0.5D0*(u%t + u%wt)
      wc = 0.5D0*(w%t + w%c)
      deltar = 0.5D0*ABS(uc)*(1 - dt*indx(i)*ABS(uc))*(u%t*dens%t - u%wt*dens%wt)
      IF (k >= nz-1) deltar = 0.D0
      gradc = indx(i)   * (u%t *dens%t  - u%wt *dens%wt )
      grade = indx(i+1) * (u%et*dens%et - u%t  *dens%t  )
      gradw = indx(i-1) * (u%wt*dens%wt - u%wwt*dens%wwt)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=2 .AND. flag(imjkp)==fluid) CALL limiters(lv,lim,erre)
      ft = ft + incrzpp*(wc-ABS(wc))*deltar*lim
!
! et-t
      uc = 0.5D0*(u%t + u%et)
      wc = 0.5D0*(w%et + w%e)
      deltar = 0.5D0*ABS(uc)*(1 - dt*indx(i+1)*ABS(uc))*(u%et*dens%et - u%t*dens%t)
      IF (k >= nz-1) deltar = 0.D0
      gradc = indx(i+1) * (u%et *dens%et  - u%t *dens%t )
      grade = indx(ip2) * (u%eet*dens%eet - u%et*dens%et)
      gradw = indx(i-1) * (u%t  *dens%t   - u%wt*dens%wt)
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
      uc = 0.5D0*indzm*(dz(k)*u%b+dz(k-1)*u%c)
      wc = 0.5D0*indxp*(dx(i+1)*w%b+dx(i)*w%eb)
      deltar = 0.5D0*ABS(wc)*(1 - dt*indzm*ABS(wc))*(u%c*dens%c - u%b*dens%b)
      gradc = indzm  * (u%c*dens%c - u%b *dens%b )
      gradt = indzp  * (u%t*dens%t - u%c *dens%c )
      gradb = indzmm * (u%b*dens%b - u%bb*dens%bb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=2 .AND. flag(ijkm)==fluid) CALL limiters(lv,lim,erre)
      fe = fe + incrxp*(uc+ABS(uc))*deltar*lim*r(i+1)
!
! t-c
      uc = 0.5D0*indzp*(dz(k+1)*u%c+dz(k)*u%t)
      wc = 0.5D0*indxp*(dx(i+1)*w%c+dx(i)*w%e)
      deltar = 0.5D0*ABS(wc)*(1 - dt*indzp*ABS(wc))*(u%t*dens%t - u%c*dens%c)
      gradc = indzp  * (u%t *dens%t  - u%c*dens%c)
      gradt = indzpp * (u%tt*dens%tt - u%t*dens%t)
      gradb = indzm  * (u%c *dens%c  - u%b*dens%b)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(ijkp)==fluid) CALL limiters(lv,lim,erre)
      fe = fe - incrxp*(uc+ABS(uc))*deltar*lim*r(i+1)
!
! et-e
      uc = 0.5D0*indzp*(dz(k)*u%et+dz(k+1)*u%e)
      wc = 0.5D0*indxpp*(dx(ip2)*w%e+dx(i+1)*w%ee)
      deltar = 0.5D0*ABS(wc)*(1 - dt*indzp*ABS(wc))*(u%et*dens%et - u%e*dens%e)
      IF (i >= nx-1) deltar = 0.D0
      gradc = indzp  * (u%et *dens%et  - u%e*dens%e  )
      gradt = indzpp * (u%ett*dens%ett - u%et*dens%et)
      gradb = indzm  * (u%e  *dens%e   - u%eb*dens%eb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(ipjkp)==fluid) CALL limiters(lv,lim,erre)
      fe = fe - incrxpp*(uc-ABS(uc))*deltar*lim*r(i+1)
!
! e-eb
      uc = 0.5D0*indzm*(dz(k)*u%eb+dz(k-1)*u%e)
      wc = 0.5D0*indxpp*(dx(ip2)*w%eb+dx(i+1)*w%eeb)
      deltar = 0.5D0*ABS(wc)*(1 - dt*indzm*ABS(wc))*(u%e*dens%e - u%eb*dens%eb)
      IF (i >= nx-1) deltar = 0.D0
      gradc = indzm  * (u%e *dens%e  - u%eb *dens%eb )
      gradt = indzp  * (u%et*dens%et - u%e  *dens%e  )
      gradb = indzmm * (u%eb*dens%eb - u%ebb*dens%ebb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=2 .AND. flag(ipjkm)==fluid) CALL limiters(lv,lim,erre)
      fe = fe + incrxpp*(uc-ABS(uc))*deltar*lim*r(i+1)
!
      RETURN
      END SUBROUTINE ctu3_flu_2d
!-----------------------------------------------------------------------
      SUBROUTINE ctu1_flu_3d(fe, fn, ft, dens, u, v, w, i, j, k)
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
      REAL*8 :: indxp, indxm, indyp, indym, indzp, indzm, indxpp
      REAL*8 :: incrxm, incrxp, incrxpp, incryp, incrym, incrzp, incrzm
      REAL*8 :: dxp, dxm, dyp, dym, dzp, dzm, dxpp
      REAL*8 :: dxi, dxip, dxipp, dxim, dyj, dyjp, dyjm, dzk, dzkp, dzkm
      INTEGER :: ip2
!
      onethird = 1.D0/3.D0
      ip2 = MIN( nx, i+2 )
!
      dxi   = dx(i)
      dxip  = dx(i+1)
      dxipp = dx(ip2)
      dxim  = dx(i-1)
      dyj   = dy(j)
      dyjp  = dy(j+1)
      dyjm  = dy(j-1)
      dzk   = dz(k)
      dzkp  = dz(k+1)
      dzkm  = dz(k-1)
!
      dxm  = dxi  + dxim
      dxp  = dxi  + dxip
      dxpp = dxip + dxipp
      dym  = dyj  + dyjm
      dyp  = dyj  + dyjp
      dzm  = dzk  + dzkm
      dzp  = dzk  + dzkp
!
      indxp  = 2.D0/dxp
      indxm  = 2.D0/dxm
      indxpp = 2.D0/dxpp
      indyp  = 2.D0/dyp
      indym  = 2.D0/dym
      indzp  = 2.D0/dzp
      indzm  = 2.D0/dzm
!
      incryp  = 0.5D0*dt*indyp
      incrym  = 0.5D0*dt*indyp
      incrzp  = 0.5D0*dt*indzp
      incrzm  = 0.5D0*dt*indzp
      incrxp  = 0.125D0*dt*indxp
      incrxm  = 0.125D0*dt*indxm
      incrxpp = 0.125D0*dt*indxpp
!
! t-c
      deltar  = dens%t*u%t - dens%c*u%c
      uc = incrzp*(dzkp*u%c+dzk*u%t)*indxp
      vc = indzp*incrxp*(dxip*dzkp*(v%c +v%s  ) + &
                         dxip*dzk *(v%t +v%st ) + &
                         dxi *dzkp*(v%e +v%es ) + &
                         dxi *dzk *(v%et+v%est))*indy(j)
      wc = 0.5D0*indxp*(dxip*w%c+dxi*w%e)
      fe = fe - 0.125D0*(uc+ABS(uc))*(wc-ABS(wc))*deltar
      fn = fn - 0.125D0*(vc+ABS(vc))*(wc-ABS(wc))*deltar &
              + 0.25D0*onethird*(vc+ABS(vc))*(wc-ABS(wc))*ABS(uc)*deltar 
!
! c-b
      deltar  = dens%c*u%c - dens%b*u%b
      uc = incrzm*(dzkm*u%c+dzk*u%b)*indxp
      vc = indzm*incrxp*(dxip*dzkm*(v%c +v%s  ) + &
                         dxip*dzk *(v%b +v%sb ) + &
                         dxi *dzkm*(v%e +v%es ) + &
                         dxi *dzk *(v%eb+v%esb))*indy(j)
      wc = 0.5D0*indxp*(dxip*w%b+dxi*w%eb)
      fe = fe - 0.125D0*(uc+ABS(uc))*(wc+ABS(wc))*deltar
      fn = fn - 0.125D0*(vc+ABS(vc))*(wc+ABS(wc))*deltar &
              + 0.25D0*onethird*(vc+ABS(vc))*(wc+ABS(wc))*ABS(uc)*deltar 
!
! n-c
      deltar  = dens%n*u%n - dens%c*u%c
      uc = incryp*(dyjp*u%c+dyj*u%n)*indxp
      vc = 0.5D0*indxp*(dxip*v%c+dxi*v%e)
      wc = indyp*incrxp*(dxip*dyjp*(w%c +w%b  ) + &
                         dxip*dyj *(w%n +w%nb ) + &
                         dxi *dyjp*(w%e +w%eb ) + &
                         dxi *dyj *(w%en+w%enb))*indz(k)
      ft = ft - 0.125D0*(wc+ABS(wc))*(vc-ABS(vc))*deltar
      fe = fe - 0.125D0*(uc+ABS(uc))*(vc-ABS(vc))*deltar &
              + 0.25D0*onethird*(uc+ABS(uc))*(vc-ABS(vc))*ABS(wc)*deltar 
!
! c-s
      deltar  = dens%c*u%c - dens%s*u%s
      uc = incrym*(dyjm*u%c+dyj*u%s)*indxp
      vc = 0.5D0*indxp*(dxip*v%s+dxi*v%es)
      wc = indym*incrxp*(dxip*dyjm*(w%c +w%b  ) + &
                         dxip*dyj *(w%s +w%sb ) + &
                         dxi *dyjm*(w%e +w%eb ) + &
                         dxi *dyj *(w%es+w%esb))*indz(k)
      ft = ft - 0.125D0*(wc+ABS(wc))*(vc+ABS(vc))*deltar
      fe = fe - 0.125D0*(uc+ABS(uc))*(vc+ABS(vc))*deltar &
              + 0.25D0*onethird*(uc+ABS(uc))*(vc+ABS(vc))*ABS(wc)*deltar 
!
! e-c
      deltar  = dens%e*u%e - dens%c*u%c
      uc = 0.5D0*(u%e+u%c)
      vc = 0.5D0*(v%e+v%es)*dt*indy(j)
      wc = 0.5D0*(w%e+w%eb)*dt*indz(k)
      fn = fn - 0.125D0*(vc+ABS(vc))*(uc-ABS(uc))*deltar
      ft = ft - 0.125D0*(wc+ABS(wc))*(uc-ABS(uc))*deltar &
              + 0.25D0*onethird*(wc+ABS(wc))*(uc-ABS(uc))*ABS(vc)*deltar 
!
! c-w
      deltar  = dens%c*u%c - dens%w*u%w
      uc = 0.5D0*(u%c+u%w)
      vc = 0.5D0*(v%c+v%s)*dt*indy(j)
      wc = 0.5D0*(w%c+w%b)*dt*indz(k)
      fn = fn - 0.125D0*(vc+ABS(vc))*(uc+ABS(uc))*deltar
      ft = ft - 0.125D0*(wc+ABS(wc))*(uc+ABS(uc))*deltar &
              + 0.25D0*onethird*(wc+ABS(wc))*(uc+ABS(uc))*ABS(vc)*deltar 
!
! wt-w
      deltar  = dens%wt*u%wt - dens%w*u%w
      IF (i <= 2) deltar = 0.D0
      uc = incrzp*(dzkp*u%w+dzk*u%wt)*indxm
      vc = indzp*incrxm*(dxim*dzkp*(v%c +v%s  ) + &
                         dxim*dzk *(v%t +v%st ) + &
                         dxi *dzkp*(v%w +v%ws ) + &
                         dxi *dzk *(v%wt+v%wst))*indy(j)
      wc = 0.5D0*indxm*(dxim*w%c+dxi*w%w)
      fn = fn - 0.125D0*onethird*(wc-ABS(wc))*(vc+ABS(vc))*ABS(uc+ABS(uc))*deltar
!
! w-wb
      deltar  = dens%w*u%w - dens%wb*u%wb
      IF (i <= 2) deltar = 0.D0
      uc = incrzm*(dzkm*u%w+dzk*u%wb)*indxm
      vc = indzm*incrxm*(dxim*dzkm*(v%c +v%s  ) + &
                         dxim*dzk *(v%b +v%sb ) + &
                         dxi *dzkm*(v%w +v%ws ) + &
                         dxi *dzk  *(v%wb+v%wsb))*indy(j)
      wc = 0.5D0*indxm*(dxim*w%b+dxi*w%wb)
      fn = fn - 0.125D0*onethird*(wc+ABS(wc))*(vc+ABS(vc))*ABS(uc+ABS(uc))*deltar
!
! et-e
      deltar  = dens%et*u%et - dens%e*u%e
      IF (i >= nx-1) deltar = 0.D0
      uc = incrzp*(dzkp*u%e+dzk*u%et)*indxpp
      vc = indzp*incrxpp*(dxipp*dzkp*(v%e  +v%es  ) + &
                          dxipp*dzk *(v%et +v%est ) + &
                          dxip *dzkp*(v%ee +v%ees ) + &
                          dxip *dzk *(v%eet+v%eest))*indy(j)
      wc = 0.5D0*indxpp*(dxipp*w%e+dxip*w%ee)
      fe = fe - 0.125D0*(wc-ABS(wc))*(uc-ABS(uc))*deltar
      fn = fn - 0.125D0*onethird*(wc-ABS(wc))*(vc+ABS(vc))*ABS(uc-ABS(uc))*deltar
!
! e-eb
      deltar  = dens%e*u%e - dens%eb*u%eb
      IF (i >= nx-1) deltar = 0.D0
      uc = incrzm*(dzkm*u%e+dzk*u%eb)*indxpp
      vc = indzm*incrxpp*(dxipp*dzkm*(v%e  +v%es  ) + &
                          dxipp*dzk *(v%eb +v%esb ) + &
                          dxip *dzkm*(v%ee +v%ees ) + &
                          dxip *dzk *(v%eeb+v%eesb))*indy(j)
      wc = 0.5D0*indxpp*(dxipp*w%eb+dxip*w%eeb)
      fe = fe - 0.125D0*(wc+ABS(wc))*(uc-ABS(uc))*deltar
      fn = fn - 0.125D0*onethird*(wc+ABS(wc))*(vc+ABS(vc))*ABS(uc-ABS(uc))*deltar
!
! en-e
      deltar  = dens%en*u%en - dens%e*u%e
      IF (i >= nx-1) deltar = 0.D0
      uc = incryp*(dyjp*u%e+dyj*u%en)*indxpp
      vc = 0.5D0*indxpp*(dxipp*v%e+dxip*v%ee)
      wc = indyp*incrxpp*(dxipp*dyjp*(w%e  +w%eb  ) + &
                          dxipp*dyj *(w%en +w%enb ) + &
                          dxip *dyjp*(w%ee +w%eeb ) + &
                          dxip *dyj *(w%een+w%eenb))*indz(k)
      fe = fe - 0.125D0*(vc-ABS(vc))*(uc-ABS(uc))*deltar &
              + 0.25D0*onethird*(vc-ABS(vc))*(uc-ABS(uc))*ABS(wc)*deltar
!
! e-es
      deltar  = dens%e*u%e - dens%es*u%es
      IF (i >= nx-1) deltar = 0.D0
      uc = incrym*(dyjm*u%e+dyj*u%es)*indxpp
      vc = 0.5D0*indxpp*(dxipp*v%es+dx(i+1)*v%ees)
      wc = indym*incrxpp*(dxipp*dyjm*(w%e  +w%eb  ) + &
                          dxipp*dyj *(w%es +w%esb ) + &
                          dxip *dyjm*(w%ee +w%eeb ) + &
                          dxip *dyj *(w%ees+w%eesb))*indz(k)
      fe = fe - 0.125D0*(vc+ABS(vc))*(uc-ABS(uc))*deltar &
              + 0.25D0*onethird*(vc+ABS(vc))*(uc-ABS(uc))*ABS(wc)*deltar
!
! nt-n
      deltar  = dens%nt*u%nt - dens%n*u%n
      IF (j >= ny-1) deltar = 0.D0
      uc = incrzp*(dzkp*u%n+dzk*u%nt)*indxp
      vc = indzp*incrxp*(dxip*dzkp*(v%n  +v%c ) + &
                         dxip*dzk *(v%nt +v%t ) + &
                         dxi *dzkp*(v%en +v%e ) + &
                         dxi *dzk *(v%ent+v%et))*indy(j+1)
      wc = 0.5D0*indxp*(dxip*w%n+dxi*w%en)
      fn = fn - 0.125D0*(vc-ABS(vc))*(wc-ABS(wc))*deltar &
              + 0.25D0*onethird*(vc-ABS(vc))*(wc-ABS(wc))*ABS(uc)*deltar
!
! n-nb
      deltar  = dens%n*u%n - dens%nb*u%nb
      IF (j >= ny-1) deltar = 0.D0
      uc = incrzm*(dzkm*u%n+dzk*u%nb)*indxp
      vc = indzm*incrxp*(dxip*dzkm*(v%n  +v%c ) + &
                         dxip*dzk *(v%nb +v%b ) + &
                         dxi *dzkm*(v%en +v%e ) + &
                         dxi *dzk *(v%enb+v%eb))*indy(j+1)
      wc = 0.5D0*indxp*(dxip*w%nb+dxi*w%enb)
      fn = fn - 0.125D0*(vc-ABS(vc))*(wc+ABS(wc))*deltar &
              + 0.25D0*onethird*(vc-ABS(vc))*(wc+ABS(wc))*ABS(uc)*deltar
!
! nt-t
      deltar  = dens%nt*u%nt - dens%t*u%t
      IF (k >= nz-1) deltar = 0.D0
      uc = incryp*(dyjp*u%t+dyj*u%nt)*indxp
      vc = 0.5D0*indxp*(dxip*v%t+dxi*v%et)
      wc = indyp*incrxp*(dxip*dyjp*(w%t  +w%c ) + &
                         dxip*dyj *(w%nt +w%n ) + &
                         dxi *dyjp*(w%et +w%e ) + &
                         dxi *dyj *(w%ent+w%en))*indz(k+1)
      ft = ft - 0.125D0*(vc-ABS(vc))*(wc-ABS(wc))*deltar
      fe = fe - 0.125D0*onethird*(vc-ABS(vc))*(uc+ABS(uc))*ABS(wc-ABS(wc))*deltar
!     
! t-st
      deltar  = dens%t*u%t - dens%st*u%st
      IF (k >= nz-1) deltar = 0.D0
      uc = incrym*(dyjm*u%t+dyj*u%st)*indxp
      vc = 0.5D0*indxp*(dxip*v%st+dxi*v%est)
      wc = indym*incrxp*(dxip*dyjm*(w%t  +w%c ) + &
                         dxip*dyj *(w%st +w%s ) + &
                         dxi *dyjm*(w%et +w%e ) + &
                         dxi *dyj *(w%est+w%es))*indz(k+1)
      ft = ft - 0.125D0*(vc+ABS(vc))*(wc-ABS(wc))*deltar
      fe = fe - 0.125D0*onethird*(vc+ABS(vc))*(uc+ABS(uc))*ABS(wc-ABS(wc))*deltar
!     
! nb-b
      deltar  = dens%nb*u%nb - dens%b*u%b
      IF (k <= 2) deltar = 0.D0
      uc = incryp*(dyjp*u%b+dyj*u%nb)*indxp
      vc = 0.5D0*indxp*(dxip*v%b+dxi*v%eb)
      wc = indyp*incrxp*(dxip*dyjp*(w%b  +w%c ) + &
                         dxip*dyj *(w%nb +w%n ) + &
                         dxi *dyjp*(w%eb +w%e ) + &
                         dxi *dyj *(w%enb+w%en))*indz(k-1)
      fe = fe - 0.125D0*onethird*(vc-ABS(vc))*(uc+ABS(uc))*ABS(wc+ABS(wc))*deltar
!     
! b-sb
      deltar  = dens%b*u%b - dens%sb*u%sb
      IF (k <= 2) deltar = 0.D0
      uc = incrym*(dyjm*u%b+dyj*u%sb)*indxp
      vc = 0.5D0*indxp*(dxip*v%sb+dxi*v%esb)
      wc = indym*incrxp*(dxip*dyjm*(w%b  +w%c ) + &
                         dxip*dyj *(w%sb +w%s ) + &
                         dxi  *dyjm*(w%eb +w%e ) + &
                         dxi  *dyj *(w%esb+w%es))*indz(k-1)
      fe = fe - 0.125D0*onethird*(vc+ABS(vc))*(uc+ABS(uc))*ABS(wc+ABS(wc))*deltar
!     
! en-n
      deltar  = dens%en*u%en - dens%n*u%n
      IF (j >= ny-1) deltar = 0.D0
      uc = 0.5D0*(u%en+u%n  )
      vc = 0.5D0*(v%en+v%e  )*dt*indy(j+1)
      wc = 0.5D0*(w%en+w%enb)*dt*indz(k)
      fn = fn - 0.125D0*(vc-ABS(vc))*(uc-ABS(uc))*deltar
      ft = ft - 0.125D0*onethird*(uc-ABS(uc))*(wc+ABS(wc))*ABS(vc-ABS(vc))*deltar
!     
! n-wn
      deltar  = dens%n*u%n - dens%wn*u%wn
      IF (j >= ny-1) deltar = 0.D0
      uc = 0.5D0*(u%n+u%wn)
      vc = 0.5D0*(v%n+v%c )*dt*indy(j+1)
      wc = 0.5D0*(w%n+w%nb)*dt*indz(k)
      fn = fn - 0.125D0*(vc-ABS(vc))*(uc+ABS(uc))*deltar
      ft = ft - 0.125D0*onethird*(uc+ABS(uc))*(wc+ABS(wc))*ABS(vc-ABS(vc))*deltar
!     
! es-s
      deltar  = dens%es*u%es - dens%s*u%s
      IF (j <= 2) deltar = 0.D0
      uc = 0.5D0*(u%es+u%s  )
      vc = 0.5D0*(v%es+v%ess)*dt*indy(j-1)
      wc = 0.5D0*(w%es+w%esb)*dt*indz(k)
      ft = ft - 0.125D0*onethird*(uc-ABS(uc))*(wc+ABS(wc))*ABS(vc+ABS(vc))*deltar
!     
! s-ws
      deltar  = dens%s*u%s - dens%ws*u%ws
      IF (j <= 2) deltar = 0.D0
      uc = 0.5D0*(u%s+u%ws)
      vc = 0.5D0*(v%s+v%ss)*dt*indy(j-1)
      wc = 0.5D0*(w%s+w%sb)*dt*indz(k)
      ft = ft - 0.125D0*onethird*(uc+ABS(uc))*(wc+ABS(wc))*ABS(vc+ABS(vc))*deltar
!     
! et-t
      deltar  = dens%et*u%et - dens%t*u%t
      IF (k >= nz-1) deltar = 0.D0
      uc = 0.5D0*(u%et+u%t  )
      vc = 0.5D0*(v%et+v%est)*dt*indy(j)
      wc = 0.5D0*(w%et+w%e  )*dt*indz(k+1)
      ft = ft - 0.125D0*(wc-ABS(wc))*(uc-ABS(uc))*deltar &
              + 0.25D0*onethird*(wc-ABS(wc))*(uc-ABS(uc))*ABS(vc)*deltar
!     
! t-wt
      deltar  = dens%t*u%t - dens%wt*u%wt
      IF (k >= nz-1) deltar = 0.D0
      uc = 0.5D0*(u%t+u%wt)
      vc = 0.5D0*(v%t+v%st)*dt*indy(j)
      wc = 0.5D0*(w%t+w%c )*dt*indz(k+1)
      ft = ft - 0.125D0*(wc-ABS(wc))*(uc+ABS(uc))*deltar &
              + 0.25D0*onethird*(wc-ABS(wc))*(uc+ABS(uc))*ABS(vc)*deltar
!     
! wnt-wn
      deltar  = dens%wnt*u%wnt - dens%wn*u%wn
      IF (i <= 2) deltar = 0.D0
      IF (j >= ny-1) deltar = 0.D0
      uc = incrzp*(dzkp*u%wn+dzk*u%wnt)*indxm
      vc = indzp*incrxm*(dxi *dzkp*(v%wn +v%w ) + &
                         dxi *dzk *(v%wnt+v%wt) + &
                         dxim*dzkp*(v%n  +v%c ) + &
                         dxim*dzk *(v%nt +v%t ))*indy(j+1)
      wc = 0.5D0*indxm*(dxi*w%wn+dxim*w%n)
      fn = fn - 0.125D0*onethird*(wc-ABS(wc))*(vc-ABS(vc))*ABS(uc+ABS(uc))*deltar
!     
! wn-wnb
      deltar  = dens%wn*u%wn - dens%wnb*u%wnb
      IF (i <= 2) deltar = 0.D0
      IF (j >= ny-1) deltar = 0.D0
      uc = incrzm*(dzkm*u%wn+dzk*u%wnb)*indxm
      vc = indzm*incrxm*(dxi *dzkm*(v%wn +v%w ) + &
                         dxi *dzk *(v%wnb+v%wb) + &
                         dxim*dzkm*(v%n  +v%c ) + &
                         dxim*dzk *(v%nb +v%b ))*indy(j+1)
      wc = 0.5D0*indxm*(dxi*w%wnb+dxim*w%nb)
      fn = fn - 0.125D0*onethird*(wc+ABS(wc))*(vc-ABS(vc))*ABS(uc+ABS(uc))*deltar
!     
! ent-et 
      deltar  = dens%ent*u%ent - dens%et*u%et
      IF (i >= nx-1) deltar = 0.D0
      IF (k >= nz-1) deltar = 0.D0
      uc = incryp*(dyjp*u%et+dyj*u%ent)*indxpp
      vc = 0.5D0*indxpp*(dxipp*v%et+dxip*v%eet)
      wc = indyp*incrxpp*(dxipp*dyjp*(w%et  +w%e  ) + &
                          dxipp*dyj *(w%ent +w%en ) + &
                          dxip *dyjp*(w%eet +w%ee ) + &
                          dxip *dyj *(w%eent+w%een))*indz(k+1)
      fe = fe - 0.125D0*onethird*(vc-ABS(vc))*(uc-ABS(uc))*ABS(wc-ABS(wc))*deltar
!     
! et-est
      deltar  = dens%et*u%et - dens%est*u%est
      IF (i >= nx-1) deltar = 0.D0
      IF (k >= nz-1) deltar = 0.D0
      uc = incrym*(dyjm*u%et+dyj*u%est)*indxpp
      vc = 0.5D0*indxpp*(dxipp*v%est+dxip*v%eest)
      wc = indym*incrxpp*(dxipp*dyjm*(w%et  +w%e  ) + &
                          dxipp*dyj *(w%est +w%es ) + &
                          dxip *dyjm*(w%eet +w%ee ) + &
                          dxip *dyj *(w%eest+w%ees))*indz(k+1)
      fe = fe - 0.125D0*onethird*(vc+ABS(vc))*(uc-ABS(uc))*ABS(wc-ABS(wc))*deltar
!     
! enb-eb
      deltar  = dens%enb*u%enb - dens%eb*u%eb
      IF (i >= nx-1) deltar = 0.D0
      IF (k <= 2) deltar = 0.D0
      uc = incryp*(dyjp*u%eb+dyj*u%enb)*indxpp
      vc = 0.5D0*indxpp*(dxipp*v%eb+dxip*v%eeb)
      wc = indyp*incrxpp*(dxipp*dyjp*(w%eb  +w%ebb  ) + &
                          dxipp*dyj *(w%enb +w%enbb ) + &
                          dxip *dyjp*(w%eeb +w%eebb ) + &
                          dxip *dyj *(w%eenb+w%eenbb))*indz(k-1)
      fe = fe - 0.125D0*onethird*(vc-ABS(vc))*(uc-ABS(uc))*ABS(wc+ABS(wc))*deltar
!     
! eb-esb
      deltar  = dens%eb*u%eb - dens%esb*u%esb
      IF (i >= nx-1) deltar = 0.D0
      IF (k <= 2) deltar = 0.D0
      uc = incrym*(dyjm*u%eb+dyj*u%esb)*indxpp
      vc = 0.5D0*indxpp*(dxipp*v%esb+dxip*v%eesb)
      wc = indym*incrxpp*(dxipp*dyjm*(w%eb  +w%ebb  ) + &
                          dxipp*dyj *(w%esb +w%esbb ) + &
                          dxip *dyjm*(w%eeb +w%eebb ) + &
                          dxip *dyj *(w%eesb+w%eesbb))*indz(k-1)
      fe = fe - 0.125D0*onethird*(vc+ABS(vc))*(uc-ABS(uc))*ABS(wc+ABS(wc))*deltar
!     
! ent-en
      deltar  = dens%ent*u%ent - dens%en*u%en
      IF (i >= nx-1) deltar = 0.D0
      IF (j >= ny-1) deltar = 0.D0
      uc = incrzp*(dzkp*u%en+dzk*u%ent)*indxpp
      vc = indzp*incrxpp*(dxipp*dzkp*(v%en  +v%e  ) + &
                          dxipp*dzk *(v%ent +v%et ) + &
                          dxip *dzkp*(v%een +v%ee ) + &
                          dxip *dzk *(v%eent+v%eet))*indy(j+1)
      wc = 0.5D0*indxpp*(dxipp*w%en+dxip*w%een)
      fn = fn - 0.125D0*onethird*(wc-ABS(wc))*(vc-ABS(vc))*ABS(uc-ABS(uc))*deltar
!     
! en-enb
      deltar  = dens%en*u%en - dens%enb*u%enb
      IF (i >= nx-1) deltar = 0.D0
      IF (j >= ny-1) deltar = 0.D0
      uc = incrzm*(dzkm*u%en+dzk*u%enb)*indxpp
      vc = indzm*incrxpp*(dxipp*dzkm*(v%en  +v%e  ) + &
                          dxipp*dzk *(v%enb +v%eb ) + &
                          dxip *dzkm*(v%een +v%ee ) + &
                          dxip *dzk *(v%eenb+v%eeb))*indy(j+1)
      wc = 0.5D0*indxpp*(dxipp*w%enb+dxip*w%eenb)
      fn = fn - 0.125D0*onethird*(wc+ABS(wc))*(vc-ABS(vc))*ABS(uc-ABS(uc))*deltar
!     
! ent-nt
      deltar  = dens%ent*u%ent - dens%nt*u%nt
      IF (k >= nz-1) deltar = 0.D0
      IF (j >= ny-1) deltar = 0.D0
      uc = 0.5D0*(u%ent+u%nt)
      vc = 0.5D0*(v%ent+v%et)*dt*indy(j+1)
      wc = 0.5D0*(w%ent+w%en)*dt*indz(k+1)
      ft = ft - 0.125D0*onethird*(wc-ABS(wc))*(uc-ABS(uc))*ABS(vc-ABS(vc))*deltar
!     
! nt-wnt
      deltar  = dens%nt*u%nt - dens%wnt*u%wnt
      IF (k >= nz-1) deltar = 0.D0
      IF (j >= ny-1) deltar = 0.D0
      uc = 0.5D0*(u%nt+u%wnt)
      vc = 0.5D0*(v%nt+v%t  )*dt*indy(j+1)
      wc = 0.5D0*(w%nt+w%n  )*dt*indz(k+1)
      ft = ft - 0.125D0*onethird*(wc-ABS(wc))*(uc+ABS(uc))*ABS(vc-ABS(vc))*deltar
!     
! est-st
      deltar  = dens%est*u%est - dens%st*u%st
      IF (k >= nz-1) deltar = 0.D0
      IF (j <= 2) deltar = 0.D0
      uc = 0.5D0*(u%est+u%st  )
      vc = 0.5D0*(v%est+v%esst)*dt*indy(j-1)
      wc = 0.5D0*(w%est+w%es  )*dt*indz(k+1)
      ft = ft - 0.125D0*onethird*(wc-ABS(wc))*(uc-ABS(uc))*ABS(vc+ABS(vc))*deltar
!     
! st-wst
      deltar  = dens%st*u%st - dens%wst*u%wst
      IF (k >= nz-1) deltar = 0.D0
      IF (j <= 2) deltar = 0.D0
      uc = 0.5D0*(u%st+u%wst)
      vc = 0.5D0*(v%st+v%sst)*dt*indy(j-1)
      wc = 0.5D0*(w%st+w%s  )*dt*indz(k+1)
      ft = ft - 0.125D0*onethird*(wc-ABS(wc))*(uc+ABS(uc))*ABS(vc+ABS(vc))*deltar
!
      RETURN
      END SUBROUTINE ctu1_flu_3d
!----------------------------------------------------------------------
      SUBROUTINE ctu2_flu_3d(fe, fn, ft, dens, u, v, w, i, j, k)
!
! ... Compute the second order Corner Transport Upwind correction (step 3 in LeVeque algorithm)
!
      USE dimensions
      USE domain_mapping, ONLY: myijk, meshinds
      USE grid, ONLY: dx, dy, dz, indx, indy, indz, fluid, flag
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: stencil
      USE set_indexes, ONLY: imjk, ijkm, ipjk, ijkp, ijmk, ijpk
      USE time_parameters, ONLY: dt
      USE flux_limiters, ONLY: lv, limiters
!
      IMPLICIT NONE

      REAL*8, INTENT(INOUT) :: fe, ft, fn
      TYPE(stencil), INTENT(IN) :: dens, u, v, w 
      INTEGER, INTENT(IN) :: i, j, k
!
      INTEGER :: ip2, jp2, kp2
      REAL*8 :: uref, vref, wref, deltar, erre, lim
      REAL*8 :: indxpp, indxp, indxm
      REAL*8 :: indypp, indyp, indym
      REAL*8 :: indzpp, indzp, indzm
      REAL*8 :: dxpp, dxp, dxm
      REAL*8 :: dypp, dyp, dym
      REAL*8 :: dzpp, dzp, dzm
      REAL*8 :: gradc, grade, gradw, gradn, grads, gradt, gradb
!
      ip2 = MIN( nx, i+2 )
      jp2 = MIN( ny, j+2 )
      kp2 = MIN( nz, k+2 )
!
      dxp  = dx(i)  + dx(i+1)
      dxpp = dx(i+1)+ dx(ip2)
      dym  = dy(j)  + dy(j-1)
      dyp  = dy(j)  + dy(j+1)
      dypp = dy(j+1)+ dy(jp2)
      dzm  = dz(k)  + dz(k-1)
      dzp  = dz(k)  + dz(k+1)
      dzpp = dz(k+1)+ dz(kp2)
!
      indxp  = 2.D0/dxp
      indxpp = 2.D0/dxpp
      indym  = 2.D0/dym
      indyp  = 2.D0/dyp
      indypp = 2.D0/dypp
      indzm  = 2.D0/dzm
      indzp  = 2.D0/dzp
      indzpp = 2.D0/dzpp
!
! e-c
      gradc = indx(i+1) * (u%e *dens%e  - u%c*dens%c)
      grade = indx(ip2) * (u%ee*dens%ee - u%e*dens%e)
      gradw = indx(i)   * (u%c *dens%c  - u%w*dens%w)
      lim  = 0.D0
      erre = 0.D0
      uref = 0.5D0*(u%c + u%e)
      IF (uref >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=nx-1 .AND. flag(ipjk)==fluid) CALL limiters(lv,lim,erre)
      uref = ABS(0.5D0*(u%c + u%e))
      deltar = u%e*dens%e - u%c*dens%c
      fe = fe + 0.5D0*uref*(1-dt*indx(i+1)*uref)*deltar*lim
!
! n-c
      gradc = indyp  * (u%n *dens%n  - u%c*dens%c)
      gradn = indypp * (u%nn*dens%nn - u%n*dens%n)
      grads = indym  * (u%c *dens%c  - u%s*dens%s)
      lim  = 0.D0
      erre = 0.D0
      vref = 0.5D0*indxp*(dx(i)*v%e + dx(i+1)*v%c)
      IF (vref >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = grads/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradn/gradc
      END IF
      IF (j/=ny-1 .AND. flag(ijpk)==fluid) CALL limiters(lv,lim,erre)
      vref = ABS(0.5D0*indxp*(dx(i)*v%e + dx(i+1)*v%c))
      deltar = u%n*dens%n - u%c*dens%c
      fn = fn + 0.5D0*vref*(1-dt*indyp*vref)*deltar*lim
!
! t-c
      gradc = indzp  * (u%t *dens%t  - u%c*dens%c)
      gradt = indzpp * (u%tt*dens%tt - u%t*dens%t)
      gradb = indzm  * (u%c *dens%c  - u%b*dens%b)
      lim  = 0.D0
      erre = 0.D0
      wref = 0.5D0*indxp*(dx(i)*w%e + dx(i+1)*w%c)
      IF (wref >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(ijkp)==fluid) CALL limiters(lv,lim,erre)
      wref = ABS(0.5D0*indxp*(dx(i)*w%e + dx(i+1)*w%c))
      deltar = u%t*dens%t - u%c*dens%c
      ft = ft + 0.5D0*wref*(1-dt*indzp*wref)*deltar*lim
!
      RETURN
      END SUBROUTINE ctu2_flu_3d
!----------------------------------------------------------------------
      SUBROUTINE ctu3_flu_3d(fe, fn, ft, dens, u, v, w, i, j, k)
!
! ... Compute the second order Corner Transport Upwind correction (step 4 in LeVeque algorithm)
!
      USE dimensions
      USE domain_mapping, ONLY: myijk
      USE grid, ONLY: dx, dy, dz, indx, indy, indz, fluid, flag
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: stencil
      USE set_indexes, ONLY: imjk, ijkm, ipjk, ijkp, ijmk, ijpk, ijpk, ijmk
      USE set_indexes, ONLY: ipjpk, ipjmk, imjpk, imjmk, ijpkp, ijmkm, ijpkm, ijmkp
      USE set_indexes, ONLY: ipjkp, ipjkm, imjkp, imjkm
      USE set_indexes, ONLY: ipjpkp, ipjpkm, ipjpkm, imjpkp, ipjmkp, imjmkp, imjpkm, ipjmkm
      USE time_parameters, ONLY: dt
      USE flux_limiters, ONLY: limiters, lv

      IMPLICIT NONE

      REAL*8, INTENT(INOUT) :: fe, ft, fn
      TYPE(stencil), INTENT(IN) :: dens, u, v, w 
      INTEGER, INTENT(IN) :: i, j, k 

      INTEGER :: im2, ip2, jm2, jp2, km2, kp2
      REAL*8 :: uc, vc, wc, deltar, erre, lim, s
      REAL*8 :: indxpp, indxp, indxm, indxmm
      REAL*8 :: indypp, indyp, indym, indymm
      REAL*8 :: indzpp, indzp, indzm, indzmm
      REAL*8 :: dxpp, dxp, dxm, dxmm
      REAL*8 :: dypp, dyp, dym, dymm
      REAL*8 :: dzpp, dzp, dzm, dzmm
      REAL*8 :: incrxm, incrxp, incrxpp, incrym, incryp, incrzm, incrzp
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
      incrxm  = 0.125D0*dt*indxm
      incrxp  = 0.125D0*dt*indxp
      incrxpp = 0.125D0*dt*indxpp
      incrym  = 0.5D0*dt*indym
      incryp  = 0.5D0*dt*indyp
      incrzm  = 0.5D0*dt*indzm
      incrzp  = 0.5D0*dt*indzp
!
! c-w
      deltar  = dens%c*u%c - dens%w*u%w
      uc = 0.5D0*(u%c+u%w)
      vc = 0.5D0*(v%c+v%s)*dt*indy(j)
      wc = 0.5D0*(w%c+w%b)*dt*indz(k)
      gradc = indx(i)   * (u%c*dens%c - u%w *dens%w )
      grade = indx(i+1) * (u%e*dens%e - u%c *dens%c )
      gradw = indx(i-1) * (u%w*dens%w - u%ww*dens%ww)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=2 .AND. flag(imjk)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indx(i)*ABS(uc))*deltar*lim
      fn = fn + 0.5D0*(vc+ABS(vc))*s
      ft = ft + 0.5D0*(1-ABS(vc))*(wc+ABS(wc))*s
!
! e-c
      deltar  = dens%e*u%e - dens%c*u%c
      uc = 0.5D0*(u%e+u%c)
      vc = 0.5D0*(v%e+v%es)*dt*indy(j)
      wc = 0.5D0*(w%e+w%eb)*dt*indz(k)
      gradc = indx(i+1) * (u%e *dens%e  - u%c*dens%c)
      grade = indx(ip2) * (u%ee*dens%ee - u%e*dens%e)
      gradw = indx(i)   * (u%c *dens%c  - u%w*dens%w)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=nx-1 .AND. flag(ipjk)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indx(i+1)*ABS(uc))*deltar*lim
      fn = fn - 0.5D0*(vc+ABS(vc))*s
      ft = ft - 0.5D0*(1-ABS(vc))*(wc+ABS(wc))*s
!
! en-n
      deltar  = dens%en*u%en - dens%n*u%n
      uc = 0.5D0*(u%en+u%n)
      vc = 0.5D0*(v%en+v%e)*dt*indy(j+1)
      wc = 0.5D0*(w%en+w%enb)*dt*indz(k)
      gradc = indx(i+1) * (u%en *dens%en  - u%n *dens%n )
      grade = indx(ip2) * (u%een*dens%een - u%en*dens%en)
      gradw = indx(i)   * (u%n  *dens%n   - u%wn*dens%wn)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=nx-1 .AND. flag(ipjpk)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indx(i+1)*ABS(uc))*deltar*lim
      fn = fn - 0.5D0*(vc-ABS(vc))*s
      ft = ft - 0.25D0*ABS(vc-ABS(vc))*(wc+ABS(wc))*s
!
! n-wn
      deltar  = dens%n*u%n - dens%wn*u%wn
      uc = 0.5D0*(u%n+u%wn)
      vc = 0.5D0*(v%n+v%c)*dt*indy(j+1)
      wc = 0.5D0*(w%n+w%nb)*dt*indz(k)
      gradc = indx(i)   * (u%n *dens%n  - u%wn *dens%wn )
      grade = indx(i+1) * (u%en*dens%en - u%n  *dens%n  )
      gradw = indx(i-1) * (u%wn*dens%wn - u%wwn*dens%wwn)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=2 .AND. flag(imjpk)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indx(i)*ABS(uc))*deltar*lim
      fn = fn + 0.5D0*(vc-ABS(vc))*s
      ft = ft + 0.25D0*ABS(vc-ABS(vc))*(wc+ABS(wc))*s
!
! es-s
      deltar  = dens%es*u%es - dens%s*u%s
      uc = 0.5D0*(u%es+u%s)
      vc = 0.5D0*(v%es+v%ess)*dt*indy(j-1)
      wc = 0.5D0*(w%es+w%esb)*dt*indz(k)
      gradc = indx(i+1) * (u%es *dens%es  - u%s *dens%s )
      grade = indx(ip2) * (u%ees*dens%ees - u%es*dens%es)
      gradw = indx(i)   * (u%s  *dens%s   - u%ws*dens%ws)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=nx-1 .AND. flag(ipjmk)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indx(i+1)*ABS(uc))*deltar*lim
      ft = ft - 0.25D0*ABS(vc+ABS(vc))*(wc+ABS(wc))*s
!
! s-ws
      deltar  = dens%s*u%s - dens%ws*u%ws
      uc = 0.5D0*(u%s+u%ws)
      vc = 0.5D0*(v%s+v%ss)*dt*indy(j-1)
      wc = 0.5D0*(w%s+w%sb)*dt*indz(k)
      gradc = indx(i)   * (u%s*dens%s   - u%ws *dens%ws )
      grade = indx(i+1) * (u%es*dens%es - u%s  *dens%s  )
      gradw = indx(i-1) * (u%ws*dens%ws - u%wws*dens%wws)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=2 .AND. flag(imjmk)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indx(i)*ABS(uc))*deltar*lim
      ft = ft + 0.25D0*ABS(vc+ABS(vc))*(wc+ABS(wc))*s
!
! et-t
      deltar  = dens%et*u%et - dens%t*u%t
      uc = 0.5D0*(u%et+u%t)
      vc = 0.5D0*(v%et+v%est)*dt*indy(j)
      wc = 0.5D0*(w%et+w%e)*dt*indz(k+1)
      gradc = indx(i+1) * (u%et *dens%et  - u%t *dens%t )
      grade = indx(ip2) * (u%eet*dens%eet - u%et*dens%et)
      gradw = indx(i)   * (u%t  *dens%t   - u%wt*dens%wt)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=nx-1 .AND. flag(ipjkp)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indx(i+1)*ABS(uc))*deltar*lim
      ft = ft - 0.5D0*(1-ABS(vc))*(wc-ABS(wc))*s
!
! t-wt
      deltar  = dens%t*u%t - dens%wt*u%wt
      uc = 0.5D0*(u%t+u%wt)
      vc = 0.5D0*(v%t+v%st)*dt*indy(j)
      wc = 0.5D0*(w%t+w%c)*dt*indz(k+1)
      gradc = indx(i)   * (u%t *dens%t  - u%wt *dens%wt )
      grade = indx(i+1) * (u%et*dens%et - u%t  *dens%t  )
      gradw = indx(i-1) * (u%wt*dens%wt - u%wwt*dens%wwt)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=2 .AND. flag(imjkp)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indx(i)*ABS(uc))*deltar*lim
      ft = ft + 0.5D0*(1-ABS(vc))*(wc-ABS(wc))*s
!
! ent-nt
      deltar  = dens%ent*u%ent - dens%nt*u%nt
      uc = 0.5D0*(u%ent+u%nt)
      vc = 0.5D0*(v%ent+v%et)*dt*indy(j+1)
      wc = 0.5D0*(w%ent+w%en)*dt*indz(k+1)
      gradc = indx(i+1) * (u%ent *dens%ent  - u%nt *dens%nt )
      grade = indx(ip2) * (u%eent*dens%eent - u%ent*dens%ent)
      gradw = indx(i)   * (u%nt  *dens%nt   - u%wnt*dens%wnt)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=nx-1 .AND. flag(ipjpk)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indx(i+1)*ABS(uc))*deltar*lim
      ft = ft - 0.25D0*ABS(vc-ABS(vc))*(wc-ABS(wc))*s
!
! nt-wnt
      deltar  = dens%nt*u%nt - dens%wnt*u%wnt
      uc = 0.5D0*(u%nt+u%wnt)
      vc = 0.5D0*(v%nt+v%t)*dt*indy(j+1)
      wc = 0.5D0*(w%nt+w%n)*dt*indz(k+1)
      gradc = indx(i)   * (u%nt *dens%nt  - u%wnt *dens%wnt )
      grade = indx(i+1) * (u%ent*dens%ent - u%nt  *dens%nt  )
      gradw = indx(i-1) * (u%wnt*dens%wnt - u%wwnt*dens%wwnt)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=2 .AND. flag(imjpkp)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indx(i)*ABS(uc))*deltar*lim
      ft = ft + 0.25D0*ABS(vc-ABS(vc))*(wc-ABS(wc))*s
!
! est-st
      deltar  = dens%est*u%est - dens%st*u%st
      uc = 0.5D0*(u%est+u%st)
      vc = 0.5D0*(v%est+v%esst)*dt*indy(j-1)
      wc = 0.5D0*(w%est+w%es)*dt*indz(k+1)
      gradc = indx(i+1) * (u%est *dens%est  - u%st *dens%st )
      grade = indx(ip2) * (u%eest*dens%eest - u%est*dens%est)
      gradw = indx(i)   * (u%st  *dens%st   - u%wst*dens%wst)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=nx-1 .AND. flag(ipjmkp)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indx(i+1)*ABS(uc))*deltar*lim
      ft = ft - 0.25D0*ABS(vc+ABS(vc))*(wc-ABS(wc))*s
!
! st-wst
      deltar  = dens%st*u%st - dens%wst*u%wst
      uc = 0.5D0*(u%st+u%wst)
      vc = 0.5D0*(v%st+v%sst)*dt*indy(j-1)
      wc = 0.5D0*(w%st+w%s)*dt*indz(k+1)
      gradc = indx(i)   * (u%st *dens%st  - u%wst *dens%wst )
      grade = indx(i+1) * (u%est*dens%est - u%st  *dens%st  )
      gradw = indx(i-1) * (u%wst*dens%wst - u%wwst*dens%wwst)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=2 .AND. flag(imjmkp)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indx(i)*ABS(uc))*deltar*lim
      ft = ft + 0.25D0*ABS(vc+ABS(vc))*(wc-ABS(wc))*s
!
! c-s
      deltar  = dens%c*u%c - dens%s*u%s
      uc = incrym*(dy(j-1)*u%c+dy(j)*u%s)*indxp
      vc = 0.5D0*indxp*(dx(i+1)*v%s+dx(i)*v%es)
      wc = indym*incrxp*(dx(i+1)*dy(j-1)*(w%c +w%b  ) + &
                                dx(i+1)*dy(j)  *(w%s +w%sb ) + &
                                dx(i)  *dy(j-1)*(w%e +w%eb ) + &
                                dx(i)  *dy(j)  *(w%es+w%esb))*indz(k)
      gradc = indym  * (dens%c*u%c - dens%s *u%s )
      grads = indymm * (dens%s*u%s - dens%ss*u%ss)
      gradn = indyp  * (dens%n*u%n - dens%c *u%c )
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
      deltar  = dens%n*u%n - dens%c*u%c
      uc = incryp*(dy(j+1)*u%c+dy(j)*u%n)*indxp
      vc = 0.5D0*indxp*(dx(i+1)*v%c+dx(i)*v%e)
      wc = indyp*incrxp*(dx(i+1)*dy(j+1)*(w%c +w%b  ) + &
                                dx(i+1)*dy(j)  *(w%n +w%nb ) + &
                                dx(i)  *dy(j+1)*(w%e +w%eb ) + &
                                dx(i)  *dy(j)  *(w%en+w%enb))*indz(k)
      gradc = indyp  * (dens%n *u%n  - dens%c*u%c)
      grads = indym  * (dens%c *u%c  - dens%s*u%s)
      gradn = indypp * (dens%nn*u%nn - dens%n*u%n)
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
      deltar  = dens%en*u%en - dens%e*u%e
      uc = incryp*(dy(j+1)*u%e+dy(j)*u%en)*indxpp
      vc = 0.5D0*indxpp*(dx(ip2)*v%e+dx(i+1)*v%ee)
      wc = indyp*incrxpp*(dx(ip2)*dy(j+1)*(w%e  +w%eb  ) + &
                                 dx(ip2)*dy(j)  *(w%en +w%enb ) + &
                                 dx(i+1)*dy(j+1)*(w%ee +w%eeb ) + &
                                 dx(i+1)*dy(j)  *(w%een+w%eenb))*indz(k)
      gradc = indyp  * (dens%en *u%en  - dens%e*u%e)
      grads = indym  * (dens%e  *u%e   - dens%es*u%es)
      gradn = indypp * (dens%enn*u%enn - dens%en*u%en)
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
      deltar  = dens%e*u%e - dens%es*u%es
      uc = incrym*(dy(j-1)*u%e+dy(j)*u%es)*indxpp
      vc = 0.5D0*indxpp*(dx(ip2)*v%es+dx(i+1)*v%ees)
      wc = indym*incrxpp*(dx(ip2)*dy(j-1)*(w%e  +w%eb  ) + &
                                 dx(ip2)*dy(j)  *(w%es +w%esb ) + &
                                 dx(i+1)*dy(j-1)*(w%ee +w%eeb ) + &
                                 dx(i+1)*dy(j)  *(w%ees+w%eesb))*indz(k)
      gradc = indym  * (dens%e *u%e  - dens%es *u%es )
      grads = indymm * (dens%es*u%es - dens%ess*u%ess)
      gradn = indyp  * (dens%en*u%en - dens%e  *u%e  )
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
      deltar  = dens%nt*u%nt - dens%t*u%t
      uc = incryp*(dy(j+1)*u%t+dy(j)*u%nt)*indxp
      vc = 0.5D0*indxp*(dx(i+1)*v%t+dx(i)*v%et)
      wc = indyp*incrxp*(dx(i+1)*dy(j+1)*(w%t  +w%c ) + &
                                dx(i+1)*dy(j)  *(w%nt +w%n ) + &
                                dx(i)  *dy(j+1)*(w%et +w%e ) + &
                                dx(i)  *dy(j)  *(w%ent+w%en))*indz(k+1)
      gradc = indyp  * (dens%nt *u%nt  - dens%t *u%t )
      grads = indym  * (dens%t  *u%t   - dens%st*u%st)
      gradn = indypp * (dens%nnt*u%nnt - dens%nt*u%nt)
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
      deltar  = dens%t*u%t - dens%st*u%st
      uc = incrym*(dy(j-1)*u%t+dy(j)*u%st)*indxp
      vc = 0.5D0*indxp*(dx(i+1)*v%st+dx(i)*v%est)
      wc = indym*incrxp*(dx(i+1)*dy(j-1)*(w%t  +w%c ) + &
                                dx(i+1)*dy(j)  *(w%st +w%s ) + &
                                dx(i)  *dy(j-1)*(w%et +w%e ) + &
                                dx(i)  *dy(j)  *(w%est+w%es))*indz(k+1)
      gradc = indym  * (dens%t *u%t  - dens%st *u%st )
      grads = indymm * (dens%st*u%st - dens%sst*u%sst)
      gradn = indyp  * (dens%nt*u%nt - dens%t  *u%t  )
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
      deltar  = dens%nb*u%nb - dens%b*u%b
      IF (k <= 2) deltar = 0.D0
      uc = incryp*(dy(j+1)*u%b+dy(j)*u%nb)*indxp
      vc = 0.5D0*indxp*(dx(i+1)*v%b+dx(i)*v%eb)
      wc = indyp*incrxp*(dx(i+1)*dy(j+1)*(w%b  +w%c ) + &
                                dx(i+1)*dy(j)  *(w%nb +w%n ) + &
                                dx(i)  *dy(j+1)*(w%eb +w%e ) + &
                                dx(i)  *dy(j)  *(w%enb+w%en))*indz(k-1)
      gradc = indyp  * (dens%nb *u%nb  - dens%b *u%b )
      grads = indym  * (dens%b  *u%b   - dens%sb*u%sb)
      gradn = indypp * (dens%nnb*u%nnb - dens%nb*u%nb)
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
      deltar  = dens%b*u%b - dens%sb*u%sb
      IF (k <= 2) deltar = 0.D0
      uc = incrym*(dy(j-1)*u%b+dy(j)*u%sb)*indxp
      vc = 0.5D0*indxp*(dx(i+1)*v%sb+dx(i)*v%esb)
      wc = indym*incrxp*(dx(i+1)*dy(j-1)*(w%b  +w%c ) + &
                                dx(i+1)*dy(j)  *(w%sb +w%s ) + &
                                dx(i)  *dy(j-1)*(w%eb +w%e ) + &
                                dx(i)  *dy(j)  *(w%esb+w%es))*indz(k-1)
      gradc = indym  * (dens%b *u%b  - dens%sb *u%sb )
      grads = indymm * (dens%sb*u%sb - dens%ssb*u%ssb)
      gradn = indyp  * (dens%nb*u%nb - dens%b  *u%b  )
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
      deltar  = dens%ent*u%ent - dens%et*u%et
      uc = incryp*(dy(j+1)*u%et+dy(j)*u%ent)*indxpp
      vc = 0.5D0*indxpp*(dx(ip2)*v%et+dx(i+1)*v%eet)
      wc = indyp*incrxpp*(dx(ip2)*dy(j+1)*(w%et  +w%e  ) + &
                                 dx(ip2)*dy(j)  *(w%ent +w%en ) + &
                                 dx(i+1)*dy(j+1)*(w%eet +w%ee ) + &
                                 dx(i+1)*dy(j)  *(w%eent+w%een))*indz(k+1)
      gradc = indyp  * (dens%ent *u%ent  - dens%et *u%et )
      grads = indym  * (dens%et  *u%et   - dens%est*u%est)
      gradn = indypp * (dens%ennt*u%ennt - dens%ent*u%ent)
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
      deltar  = dens%et*u%et - dens%est*u%est
      uc = incrym*(dy(j-1)*u%et+dy(j)*u%est)*indxpp
      vc = 0.5D0*indxpp*(dx(ip2)*v%est+dx(i+1)*v%eest)
      wc = indym*incrxpp*(dx(ip2)*dy(j-1)*(w%et  +w%e  ) + &
                                 dx(ip2)*dy(j)  *(w%est +w%es ) + &
                                 dx(i+1)*dy(j-1)*(w%eet +w%ee ) + &
                                 dx(i+1)*dy(j)  *(w%eest+w%ees))*indz(k+1)
      gradc = indym  * (dens%et *u%et  - dens%est *u%est )
      grads = indymm * (dens%est*u%est - dens%esst*u%esst)
      gradn = indyp  * (dens%ent*u%ent - dens%et  *u%et  )
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
      deltar  = dens%enb*u%enb - dens%eb*u%eb
      IF (k <= 2) deltar = 0.D0
      uc = incryp*(dy(j+1)*u%eb+dy(j)*u%enb)*indxpp
      vc = 0.5D0*indxpp*(dx(ip2)*v%eb+dx(i+1)*v%eeb)
      wc = indyp*incrxpp*(dx(ip2)*dy(j+1)*(w%eb  +w%ebb  ) + &
                                 dx(ip2)*dy(j)  *(w%enb +w%enbb ) + &
                                 dx(i+1)*dy(j+1)*(w%eeb +w%eebb ) + &
                                 dx(i+1)*dy(j)  *(w%eenb+w%eenbb))*indz(k-1)
      gradc = indyp  * (dens%enb *u%enb  - dens%eb *u%eb )
      grads = indym  * (dens%eb  *u%eb   - dens%esb*u%esb)
      gradn = indypp * (dens%ennb*u%ennb - dens%enb*u%enb)
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
      deltar  = dens%eb*u%eb - dens%esb*u%esb
      IF (k <= 2) deltar = 0.D0
      uc = indym*(dy(j-1)*u%eb+dy(j)*u%esb)*indxpp
      vc = 0.5D0*indxpp*(dx(ip2)*v%esb+dx(i+1)*v%eesb)
      wc = indym*incrxpp*(dx(ip2)*dy(j-1)*(w%eb  +w%ebb  ) + &
                                 dx(ip2)*dy(j)  *(w%esb +w%esbb ) + &
                                 dx(i+1)*dy(j-1)*(w%eeb +w%eebb ) + &
                                 dx(i+1)*dy(j)  *(w%eesb+w%eesbb))*indz(k-1)
      gradc = indym  * (dens%eb *u%eb  - dens%esb *u%esb )
      grads = indymm * (dens%esb*u%esb - dens%essb*u%essb)
      gradn = indyp  * (dens%enb*u%enb - dens%eb  *u%eb  )
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
      deltar  = dens%c*u%c - dens%b*u%b
      uc = indzm*(dz(k-1)*u%c+dz(k)*u%b)*indxp
      vc = indzm*incrxp*(dx(i+1)*dz(k-1)*(v%c +v%s  ) + &
                                dx(i+1)*dz(k)  *(v%b +v%sb ) + &
                                dx(i)  *dz(k-1)*(v%e +v%es ) + &
                                dx(i)  *dz(k)  *(v%eb+v%esb))*indy(j)
      wc = 0.5D0*indxp*(dx(i+1)*w%b+dx(i)*w%eb)
      gradc = indzm  * (u%c*dens%c - u%b *dens%b )
      gradt = indzp  * (u%t*dens%t - u%c *dens%c )
      gradb = indzmm * (u%b*dens%b - u%bb*dens%bb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=2 .AND. flag(ijkm)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indzm*ABS(wc))*deltar*lim
      fe = fe + 0.5D0*(uc+ABS(uc))*s
      fn = fn + 0.5D0*(1-ABS(uc))*(vc+ABS(vc))*s
!
! t-c
      deltar  = dens%t*u%t - dens%c*u%c
      uc = incrzp*(dz(k+1)*u%c+dz(k)*u%t)*indxp
      vc = indzp*incrxp*(dx(i+1)*dz(k+1)*(v%c +v%s  ) + &
                                dx(i+1)*dz(k)  *(v%t +v%st ) + &
                                dx(i)  *dz(k+1)*(v%e +v%es ) + &
                                dx(i)  *dz(k)  *(v%et+v%est))*indy(j)
      wc = 0.5D0*indxp*(dx(i+1)*w%c+dx(i)*w%e)
      gradc = indzp  * (u%t *dens%t  - u%c*dens%c)
      gradt = indzpp * (u%tt*dens%tt - u%t*dens%t)
      gradb = indzm  * (u%c *dens%c  - u%b*dens%b)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(ijkp)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indzp*ABS(wc))*deltar*lim
      fe = fe - 0.5D0*(uc+ABS(uc))*s
      fn = fn - 0.5D0*(1-ABS(uc))*(vc+ABS(vc))*s
!
! wt-w
      deltar  = dens%wt*u%wt - dens%w*u%w
      uc = incrzp*(dz(k+1)*u%w+dz(k)*u%wt)*indxm
      vc = indzp*incrxm*(dx(i-1)*dz(k+1)*(v%c +v%s  ) + &
                                dx(i-1)*dz(k)  *(v%t +v%st ) + &
                                dx(i)  *dz(k+1)*(v%w +v%ws ) + &
                                dx(i)  *dz(k)  *(v%wt+v%wst))*indy(j)
      wc = 0.5D0*indxm*(dx(i-1)*w%c+dx(i)*w%w)
      gradc = indzp  * (u%wt *dens%wt  - u%w *dens%w )
      gradt = indzpp * (u%wtt*dens%wtt - u%wt*dens%wt)
      gradb = indzm  * (u%w  *dens%w   - u%wb*dens%wb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(imjkp)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indzp*ABS(wc))*deltar*lim
      fn = fn - 0.25D0*ABS(uc+ABS(uc))*(vc+ABS(vc))*s
!
! w-wb
      deltar  = dens%w*u%w - dens%wb*u%wb
      uc = incrzm*(dz(k-1)*u%w+dz(k)*u%wb)*indxm
      vc = indzm*incrxm*(dx(i-1)*dz(k-1)*(v%c +v%s  ) + &
                                dx(i-1)*dz(k)  *(v%b +v%sb ) + &
                                dx(i)  *dz(k-1)*(v%w +v%ws ) + &
                                dx(i)  *dz(k)  *(v%wb+v%wsb))*indy(j)
      wc = 0.5D0*indxm*(dx(i-1)*w%b+dx(i)*w%wb)
      gradc = indzm  * (u%w *dens%w  - u%wb *dens%wb )
      gradt = indzp  * (u%wt*dens%wt - u%w  *dens%w  )
      gradb = indzmm * (u%wb*dens%wb - u%wbb*dens%wbb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=2 .AND. flag(imjkm)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indzm*ABS(wc))*deltar*lim
      fn = fn + 0.25D0*ABS(uc+ABS(uc))*(vc+ABS(vc))*s
!
! et-e
      deltar  = dens%et*u%et - dens%e*u%e
      uc = incrzp*(dz(k+1)*u%e+dz(k)*u%et)*indxpp
      vc = indzp*incrxpp*(dx(ip2)*dz(k+1)*(v%e  +v%es  ) + &
                                 dx(ip2)*dz(k)  *(v%et +v%est ) + &
                                 dx(i+1)*dz(k+1)*(v%ee +v%ees ) + &
                                 dx(i+1)*dz(k)  *(v%eet+v%eest))*indy(j)
      wc = 0.5D0*indxpp*(dx(ip2)*w%e+dx(i+1)*w%ee)
      gradc = indzp  * (u%et *dens%et  - u%e *dens%e )
      gradt = indzpp * (u%ett*dens%ett - u%et*dens%et)
      gradb = indzm  * (u%e  *dens%e   - u%eb*dens%eb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(ipjkp)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indzp*ABS(wc))*deltar*lim
      fe = fe - 0.5D0*(uc-ABS(uc))*s
      fn = fn - 0.25D0*ABS(uc-ABS(uc))*(vc+ABS(vc))*s
!
! e-eb
      deltar  = dens%e*u%e - dens%eb*u%eb
      uc = incrzm*(dz(k-1)*u%e+dz(k)*u%eb)*indxpp
      vc = indzm*incrxpp*(dx(ip2)*dz(k-1)*(v%e  +v%es  ) + &
                                 dx(ip2)*dz(k)  *(v%eb +v%esb ) + &
                                 dx(i+1)*dz(k-1)*(v%ee +v%ees ) + &
                                 dx(i+1)*dz(k)  *(v%eeb+v%eesb))*indy(j)
      wc = 0.5D0*indxpp*(dx(ip2)*w%eb+dx(i+1)*w%eeb)
      gradc = indzm  * (u%e *dens%e  - u%eb *dens%eb )
      gradt = indzp  * (u%et*dens%et - u%e  *dens%e  )
      gradb = indzmm * (u%eb*dens%eb - u%ebb*dens%ebb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=2 .AND. flag(ipjkm)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indzm*ABS(wc))*deltar*lim
      fe = fe + 0.5D0*(uc-ABS(uc))*s
      fn = fn + 0.25D0*ABS(uc-ABS(uc))*(vc+ABS(vc))*s
!
! nt-n
      deltar  = dens%nt*u%nt - dens%n*u%n
      uc = incrzp*(dz(k+1)*u%n+dz(k)*u%nt)*indxp
      vc = indzp*incrxp*(dx(i+1)*dz(k+1)*(v%n  +v%c ) + &
                                dx(i+1)*dz(k)  *(v%nt +v%t ) + &
                                dx(i)  *dz(k+1)*(v%en +v%e ) + &
                                dx(i)  *dz(k)  *(v%ent+v%et))*indy(j+1)
      wc = 0.5D0*indxp*(dx(i+1)*w%n+dx(i)*w%en)
      gradc = indzp  * (u%nt *dens%nt  - u%n *dens%n )
      gradt = indzpp * (u%ntt*dens%ntt - u%nt*dens%nt)
      gradb = indzm  * (u%n  *dens%n   - u%nb*dens%nb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(ijpkp)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indzp*ABS(wc))*deltar*lim
      fn = fn - 0.5D0*(1-ABS(uc))*(vc-ABS(vc))*s
!
! n-nb
      deltar  = dens%n*u%n - dens%nb*u%nb
      uc = incrzm*(dz(k-1)*u%n+dz(k)*u%nb)*indxp
      vc = indzm*incrxp*(dx(i+1)*dz(k-1)*(v%n  +v%c ) + &
                                dx(i+1)*dz(k)  *(v%nb +v%b ) + &
                                dx(i)  *dz(k-1)*(v%en +v%e ) + &
                                dx(i)  *dz(k)  *(v%enb+v%eb))*indy(j+1)
      wc = 0.5D0*indxp*(dx(i+1)*w%nb+dx(i)*w%enb)
      gradc = indzm  * (u%n *dens%n  - u%nb *dens%nb )
      gradt = indzp  * (u%nt*dens%nt - u%n  *dens%n  )
      gradb = indzmm * (u%nb*dens%nb - u%nbb*dens%nbb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=2 .AND. flag(ijpkm)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indzm*ABS(wc))*deltar*lim
      fn = fn + 0.5D0*(1-ABS(uc))*(vc-ABS(vc))*s
!
! wnt-wn
      deltar  = dens%wnt*u%wnt - dens%wn*u%wn
      uc = incrzp*(dz(k+1)*u%wn+dz(k)*u%wnt)*indxm
      vc = indzp*incrxm*(dx(i)  *dz(k+1)*(v%wn  +v%w ) + &
                                dx(i)  *dz(k)  *(v%wnt +v%wt ) + &
                                dx(i-1)*dz(k+1)*(v%n +v%c ) + &
                                dx(i-1)*dz(k)  *(v%nt+v%t))*indy(j+1)
      wc = 0.5D0*indxm*(dx(i)*w%wn+dx(i-1)*w%n)
      gradc = indzp  * (u%wnt *dens%wnt  - u%wn *dens%wn )
      gradt = indzpp * (u%wntt*dens%wntt - u%wnt*dens%wnt)
      gradb = indzm  * (u%wn  *dens%wn   - u%wnb*dens%wnb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(imjpkp)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indzp*ABS(wc))*deltar*lim
      fn = fn - 0.25D0*ABS(uc+ABS(uc))*(vc-ABS(vc))*s
!
! wn-wnb
      deltar  = dens%wn*u%wn - dens%wnb*u%wnb
      uc = incrzm*(dz(k-1)*u%wn+dz(k)*u%wnb)*indxm
      vc = indzm*incrxm*(dx(i)  *dz(k-1)*(v%wn  +v%w ) + &
                                dx(i)  *dz(k)  *(v%wnb +v%wb ) + &
                                dx(i-1)*dz(k-1)*(v%n +v%c ) + &
                                dx(i-1)*dz(k)  *(v%nb+v%b))*indy(j+1)
      wc = 0.5D0*indxm*(dx(i)*w%wnb+dx(i-1)*w%nb)
      gradc = indzm  * (u%wn *dens%wn  - u%wnb *dens%wnb )
      gradt = indzp  * (u%wnt*dens%wnt - u%wn  *dens%wn  )
      gradb = indzmm * (u%wnb*dens%wnb - u%wnbb*dens%wnbb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=2 .AND. flag(imjpkm)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indzm*ABS(wc))*deltar*lim
      fn = fn + 0.25D0*ABS(uc+ABS(uc))*(vc-ABS(vc))*s
!
! ent-en
      deltar  = dens%ent*u%ent - dens%en*u%en
      uc = incrzp*(dz(k+1)*u%en+dz(k)*u%ent)*indxpp
      vc = indzp*incrxpp*(dx(ip2)*dz(k+1)*(v%en  +v%e  ) + &
                                 dx(ip2)*dz(k)  *(v%ent +v%et ) + &
                                 dx(i+1)*dz(k+1)*(v%een +v%ee ) + &
                                 dx(i+1)*dz(k)  *(v%eent+v%eet))*indy(j+1)
      wc = 0.5D0*indxpp*(dx(ip2)*w%en+dx(i+1)*w%een)
      gradc = indzp  * (u%ent *dens%ent  - u%en *dens%en )
      gradt = indzpp * (u%entt*dens%entt - u%ent*dens%ent)
      gradb = indzm  * (u%en  *dens%en   - u%enb*dens%enb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(ipjpkp)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indzp*ABS(wc))*deltar*lim
      fn = fn - 0.25D0*ABS(uc-ABS(uc))*(vc-ABS(vc))*s
!
! en-enb
      deltar  = dens%en*u%en - dens%enb*u%enb
      uc = incrzm*(dz(k-1)*u%en+dz(k)*u%enb)*indxpp
      wc = 0.5D0*indxpp*(dx(ip2)*w%enb+dx(i+1)*w%eenb)
      vc = indzm*incrxpp*(dx(ip2)*dz(k-1)*(v%en  +v%e  ) + &
                                 dx(ip2)*dz(k)  *(v%enb +v%eb ) + &
                                 dx(i+1)*dz(k-1)*(v%een +v%ee ) + &
                                 dx(i+1)*dz(k)  *(v%eenb+v%eeb))*indy(j+1)
      gradc = indzm  * (u%en *dens%en  - u%enb *dens%enb )
      gradt = indzp  * (u%ent*dens%ent - u%en  *dens%en  )
      gradb = indzmm * (u%enb*dens%enb - u%enbb*dens%enbb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=2 .AND. flag(ipjpkm)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indzm*ABS(wc))*deltar*lim
      fn = fn + 0.25D0*ABS(uc-ABS(uc))*(vc-ABS(vc))*s
!
      RETURN
      END SUBROUTINE ctu3_flu_3d
!-----------------------------------------------------------------------
      END MODULE convective_fluxes_u
!-----------------------------------------------------------------------
