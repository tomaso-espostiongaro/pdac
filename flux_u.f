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
        MODULE PROCEDURE ctu1_flu_2d
      END INTERFACE
      INTERFACE ctu2_flu
        MODULE PROCEDURE ctu2_flu_2d
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
      END MODULE convective_fluxes_u
!-----------------------------------------------------------------------
