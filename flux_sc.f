!----------------------------------------------------------------------
      MODULE convective_fluxes_sc
!
! ... This module computes convective fluxes of ____ scalar fields ____
! ... by using an upwind scheme.
! ... Standard compass notation is adopted to locate cell neighbours:
! ... %e=East, %w=West, %n=North, %s=South, %t=Top, %b=Bottom, etc.
! ... The computational stencil is defined in the `set_indexes' module
!
!----------------------------------------------------------------------
!
      USE flux_limiters, ONLY: muscl, limiters, lm

      IMPLICIT NONE
      SAVE
!
      REAL*8, PRIVATE :: cs                      ! convective stream   !
      REAL*8, PRIVATE :: cn                      ! Courant number      !
      REAL*8, PRIVATE :: fou                     ! first order upwind  !
      REAL*8, PRIVATE :: upwnd                   ! upwinded variable   !
      REAL*8, PRIVATE :: lim                     ! limiter             !
      REAL*8, PRIVATE :: erre                    ! gradients ratio     !
      REAL*8, PRIVATE :: incr                    ! increment           !
!
      INTERFACE fsc
        MODULE PROCEDURE fsc_2d, fsc_3d
      END INTERFACE
      INTERFACE muscl_fsc
        MODULE PROCEDURE muscl_fsc_2d, muscl_fsc_3d
      END INTERFACE
      INTERFACE ctu1_fsc
        MODULE PROCEDURE ctu1_fsc_2d, ctu1_fsc_3d
      END INTERFACE
      INTERFACE ctu2_fsc
        MODULE PROCEDURE ctu2_fsc_2d, ctu2_fsc_3d
      END INTERFACE
      INTERFACE ctu3_fsc
        MODULE PROCEDURE ctu3_fsc_2d, ctu3_fsc_3d
      END INTERFACE
      
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE fsc_3d(fe, fn, ft, fw, fs, fb, dens, field, u, v, w, ijk)
!
      USE dimensions
      USE domain_mapping, ONLY: myijk
      USE grid, ONLY: flag
      USE grid, ONLY: dx, dy, dz
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: imjk, ijmk, ijkm
      USE set_indexes, ONLY: stencil
      USE time_parameters, ONLY: dt
      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: fe, fn, ft, fw, fs, fb
      TYPE(stencil), INTENT(IN) :: dens, field, u, v, w
      INTEGER, INTENT(IN) :: ijk
!
! ... on West volume boundary
!
      IF ( .NOT.BTEST(flag(imjk),0) ) THEN
        cs = u%w
        fw = 0.5D0*(cs-ABS(cs))*dens%c*field%c + 0.5D0*(cs+ABS(cs))*dens%w*field%w
      END IF
!
! ... on South volume boundary
!
      IF ( .NOT.BTEST(flag(ijmk),0) ) THEN
        cs = v%s
        fs = 0.5D0*(cs-ABS(cs))*dens%c*field%c + 0.5D0*(cs+ABS(cs))*dens%s*field%s
      END IF
!
! ... on Bottom volume boundary
!
      IF ( .NOT.BTEST(flag(ijkm),0) ) THEN
        cs = w%b
        fb = 0.5D0*(cs-ABS(cs))*dens%c*field%c + 0.5D0*(cs+ABS(cs))*dens%b*field%b
      END IF
!
! ... on East volume boundary
!
      cs = u%c
      fe = 0.5D0*(cs-ABS(cs))*dens%e*field%e + 0.5D0*(cs+ABS(cs))*dens%c*field%c
!
! ... on North volume boundary
!
      cs = v%c
      fn = 0.5D0*(cs-ABS(cs))*dens%n*field%n + 0.5D0*(cs+ABS(cs))*dens%c*field%c
!
! ... on Top volume boundary
!
      cs = w%c
      ft = 0.5D0*(cs-ABS(cs))*dens%t*field%t + 0.5D0*(cs+ABS(cs))*dens%c*field%c
!
      RETURN
      END SUBROUTINE fsc_3d
!----------------------------------------------------------------------
      SUBROUTINE muscl_fsc_3d(fe, fn, ft, dens, field, u, v, w, ijk)
!
      USE dimensions
      USE domain_mapping, ONLY: myijk
      USE grid, ONLY: flag
      USE grid, ONLY: dx, dy, dz
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: stencil
      USE time_parameters, ONLY: dt
      IMPLICIT NONE
!
      REAL*8, INTENT(INOUT) :: fe, fn, ft
      TYPE(stencil), INTENT(IN) :: dens, field, u, v, w
      INTEGER, INTENT(IN) :: ijk
!
      INTEGER :: i,j,k,imesh
      REAL*8 :: dxm, dxp, dxpp, indxpp, indxp, indxm
      REAL*8 :: dym, dyp, dypp, indypp, indyp, indym
      REAL*8 :: dzp, indzp, dzm, indzm, dzpp, indzpp
      REAL*8 :: gradc, grade, gradw, gradn, grads, gradt, gradb

      INTEGER :: im2, ip2, jm2, jp2, km2, kp2
!
      imesh = myijk( ip0_jp0_kp0_, ijk)
      i = MOD( MOD( imesh - 1, nx*ny ), nx ) + 1
      j = MOD( imesh - 1, nx*ny ) / nx + 1
      k = ( imesh - 1 ) / ( nx*ny ) + 1 
!
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
      dym  = dy(j)  + dy(j-1)
      dyp  = dy(j)  + dy(j+1)
      dypp = dy(j+1)+ dy(jp2)
      dzm  = dz(k)  + dz(k-1)
      dzp  = dz(k)  + dz(k+1)
      dzpp = dz(k+1)+ dz(kp2)

      indxm  = 1.D0/dxm
      indxp  = 1.D0/dxp
      indxpp = 1.D0/dxpp
      indym  = 1.D0/dym
      indyp  = 1.D0/dyp
      indypp = 1.D0/dypp
      indzm  = 1.D0/dzm
      indzp  = 1.D0/dzp
      indzpp = 1.D0/dzpp
!
! ... On boundaries mantain first order accuracy
! ... Otherwise, MUSCL reconstruction of fields
!
! ... on East volume boundary
!
      gradc = 2.D0 * indxp  * (dens%e *field%e  - dens%c*field%c)
      gradw = 2.D0 * indxm  * (dens%c *field%c  - dens%w*field%w)
      grade = 2.D0 * indxpp * (dens%ee*field%ee - dens%e*field%e)
!
      lim = 0.D0
      erre = 0.D0
!
      cs = u%c
      IF (cs >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = gradw / gradc
	incr = 0.5D0 * dx(i)
      ELSE IF (cs < 0.D0) THEN
	IF (gradc /= 0.D0) erre = grade / gradc
	incr = - 0.5D0 * dx(i+1)
      ENDIF
!
      IF (i /= nx-1) CALL limiters(lm,lim,erre)
!
      upwnd = lim * gradc * incr
!
      ! ... MUSCL correction
      fe = fe + upwnd * cs
!
! ... on North volume boundary
!
      gradc = 2.D0 * indyp  * (dens%n *field%n  - dens%c*field%c)
      grads = 2.D0 * indym  * (dens%c *field%c  - dens%s*field%s)
      gradn = 2.D0 * indypp * (dens%nn*field%nn - dens%n*field%n)
!
      lim = 0.D0
      erre = 0.D0
!
      cs = v%c
      IF (cs >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = grads / gradc
	incr = 0.5D0 * dy(j)
      ELSE IF (cs < 0.D0) THEN
	IF (gradc /= 0.D0) erre = gradn / gradc
	incr = - 0.5D0 * dy(j+1)
      ENDIF
!
      IF (j /= ny-1) CALL limiters(lm,lim,erre)
!
      upwnd = lim * gradc * incr
!
      ! ... MUSCL correction
      fn = fn + upwnd * cs
!
! ... on Top volume boundary
!
      gradc = 2.D0 * indzp  * (dens%t *field%t  - dens%c*field%c)
      gradb = 2.D0 * indzm  * (dens%c *field%c  - dens%b*field%b)
      gradt = 2.D0 * indzpp * (dens%tt*field%tt - dens%t*field%t)
!
      lim = 0.D0
      erre = 0.D0
!
      cs = w%c
      !cn = cs * dt * 2.D0 * indzp
      IF (cs >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = gradb / gradc
	incr = 0.5D0 * dz(k)
      ELSE IF (cs < 0.D0) THEN
	IF (gradc /= 0.D0) erre = gradt / gradc
	incr = - 0.5D0 * dz(k+1)
      ENDIF
!
      IF (k /= nz-1) CALL limiters(lm,lim,erre)
!
      upwnd = lim * gradc * incr
!
      ! ... MUSCL correction
      ft = ft + upwnd * cs
!
      RETURN
      END SUBROUTINE muscl_fsc_3d
!----------------------------------------------------------------------
      SUBROUTINE fsc_2d(fe, ft, fw, fb, dens, field, u, w, ijk)
!
      USE dimensions
      USE domain_mapping, ONLY: myijk
      USE grid, ONLY: flag
      USE grid, ONLY: dz
      USE grid, ONLY: dx, rb
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: imjk, ijkm
      USE set_indexes, ONLY: stencil
      USE time_parameters, ONLY: dt

      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: fe, ft, fw, fb
      TYPE(stencil), INTENT(IN) :: dens, field, u, w
      INTEGER, INTENT(IN) :: ijk

      INTEGER :: i,k,imesh
!
      imesh = myijk( ip0_jp0_kp0_, ijk)
      k = ( imesh - 1 ) / nx + 1
      i = MOD( ( imesh - 1 ), nx) + 1
!
! ... on West volume boundary
!
      IF ( .NOT.BTEST(flag(imjk),0) ) THEN
        cs = u%w
        fw = 0.5D0*(cs-ABS(cs))*dens%c*field%c*rb(i-1) + 0.5D0*(cs+ABS(cs))*dens%w*field%w*rb(i-1)
      END IF
!
! ... on Bottom volume boundary
!
      IF ( .NOT.BTEST(flag(ijkm),0) ) THEN
        cs = w%b
        fb = 0.5D0*(cs-ABS(cs))*dens%c*field%c + 0.5D0*(cs+ABS(cs))*dens%b*field%b
      END IF
!
! ... on East volume boundary
!
      cs = u%c
      fe = 0.5D0*(cs-ABS(cs))*dens%e*field%e*rb(i) + 0.5D0*(cs+ABS(cs))*dens%c*field%c*rb(i)
!
! ... on Top volume boundary
!
      cs = w%c
      ft = 0.5D0*(cs-ABS(cs))*dens%t*field%t + 0.5D0*(cs+ABS(cs))*dens%c*field%c
!
      RETURN
      END SUBROUTINE fsc_2d
!-----------------------------------------------------------------------
      SUBROUTINE muscl_fsc_2d(fe, ft, dens, field, u, w, ijk)
!
      USE dimensions
      USE domain_mapping, ONLY: myijk
      USE grid, ONLY: flag
      USE grid, ONLY: dz
      USE grid, ONLY: dx, rb
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: stencil
      USE time_parameters, ONLY: dt

      IMPLICIT NONE
!
      REAL*8, INTENT(INOUT) :: fe, ft
      TYPE(stencil), INTENT(IN) :: dens, field, u, w
      INTEGER, INTENT(IN) :: ijk
!
      INTEGER :: im2, ip2, km2, kp2

      INTEGER :: i,k,imesh
      REAL*8 :: dxm, dxp, dxpp, indxpp, indxp, indxm
      REAL*8 :: dzp, indzp, dzm, indzm, dzpp, indzpp
      REAL*8 :: gradc, grade, gradw, gradt, gradb
!
      imesh = myijk( ip0_jp0_kp0_, ijk)
      k = ( imesh - 1 ) / nx + 1
      i = MOD( ( imesh - 1 ), nx) + 1
!
      im2 = MAX(  1, i-2 )
      ip2 = MIN( nx, i+2 )
      km2 = MAX(  1, k-2 )
      kp2 = MIN( nz, k+2 )

!
      dxm  = dx(i)  + dx(i-1)
      dxp  = dx(i)  + dx(i+1)
      dxpp = dx(i+1)+ dx(ip2)
      dzm  = dz(k)  + dz(k-1)
      dzp  = dz(k)  + dz(k+1)
      dzpp = dz(k+1)+ dz(kp2)

      indxm  = 1.D0/dxm
      indxp  = 1.D0/dxp
      indxpp = 1.D0/dxpp
      indzm  = 1.D0/dzm
      indzp  = 1.D0/dzp
      indzpp = 1.D0/dzpp
!
! ... On boundaries mantain first order accuracy (no MUSCL correction)
! ... Otherwise, correct fluxes by using a linear interpolation of
! ... gradients
!
! ... on East volume boundary
!
      gradc = 2.D0 * indxp  * (dens%e *field%e  - dens%c*field%c)
      gradw = 2.D0 * indxm  * (dens%c *field%c  - dens%w*field%w)
      grade = 2.D0 * indxpp * (dens%ee*field%ee - dens%e*field%e)
!
      lim = 0.D0
      erre = 0.D0
!
      cs = u%c
      IF (cs >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = gradw / gradc
	incr = 0.5D0 * dx(i)
      ELSE IF (cs < 0.D0) THEN
	IF (gradc /= 0.D0) erre = grade / gradc
	incr = - 0.5D0 * dx(i+1)
      ENDIF
!
      IF (i /= nx-1) CALL limiters(lm,lim,erre)
!
      upwnd = lim * gradc * incr
!
      ! ... MUSCL correction
      fe = fe + upwnd * cs * rb(i)
!
! ... on Top volume boundary
!
      gradc = 2.D0 * indzp  * (dens%t *field%t  - dens%c*field%c)
      gradb = 2.D0 * indzm  * (dens%c *field%c  - dens%b*field%b)
      gradt = 2.D0 * indzpp * (dens%tt*field%tt - dens%t*field%t)
!
      lim = 0.D0
      erre = 0.D0
!
      cs = w%c
      IF (cs >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = gradb / gradc
	incr = 0.5D0 * dz(k)
      ELSE IF (cs < 0.D0) THEN
	IF (gradc /= 0.D0) erre = gradt / gradc
	incr = - 0.5D0 * dz(k+1)
      ENDIF
!
      IF (k /= nz-1) CALL limiters(lm,lim,erre)
!
      upwnd = lim * gradc * incr
!
      ! ... MUSCL correction
      ft = ft + upwnd * cs
!
      RETURN
      END SUBROUTINE muscl_fsc_2d
!-----------------------------------------------------------------------
    SUBROUTINE ctu1_fsc_2d(fe, ft, dens, field,  u, w, ij)
!
! ... Compute the first order Corner Transport Upwind correction (step 2 in LeVeque algorithm)
!
      USE dimensions
      USE domain_mapping, ONLY: myijk
      USE grid, ONLY: dx, dz, rb
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: imjk, ijkm
      USE set_indexes, ONLY: stencil
      USE time_parameters, ONLY: dt

      IMPLICIT NONE
!
      REAL*8, INTENT(INOUT) :: fe, ft
      TYPE(stencil), INTENT(IN) :: dens, field, u, w
      INTEGER, INTENT(IN) :: ij
!
      INTEGER :: i,k,imesh
      REAL*8 :: uc, wc, deltar
      REAL*8 :: dxp, dxm, dzp, dzm, indzp, indzm, indxp, indxm
      REAL*8 :: incrxp, incrxm, incrzp, incrzm
!
      imesh = myijk( ip0_jp0_kp0_, ij)
      k = ( imesh - 1 ) / nx + 1
      i = MOD( ( imesh - 1 ), nx) + 1
!
      dxm = dx(i) + dx(i-1)
      dxp = dx(i) + dx(i+1)
      dzm = dz(k) + dz(k-1)
      dzp = dz(k) + dz(k+1)
!
      indxm = 2.D0/dxm
      indxp = 2.D0/dxp
      indzm = 2.D0/dzm
      indzp = 2.D0/dzp
!
      incrxp = 0.125D0*dt*indxp
      incrxm = 0.125D0*dt*indxm
      incrzp = 0.125D0*dt*indzp
      incrzm = 0.125D0*dt*indzm
!
! c-w
      uc = u%w
      wc = 0.25D0*indxm*(dx(i-1)*w%c+dx(i)*w%w+dx(i-1)*w%b+dx(i)*w%wb)
      deltar  = dens%c*field%c - dens%w*field%w
      ft = ft - incrxm*(uc+ABS(uc))*(wc+ABS(wc))*deltar
!
! e-c
      uc = u%c
      wc = 0.25D0*indxp*(dx(i+1)*w%c+dx(i)*w%e+dx(i+1)*w%b+dx(i)*w%eb)
      deltar  = dens%e*field%e - dens%c*field%c
      ft = ft - incrxp*(uc-ABS(uc))*(wc+ABS(wc))*deltar
!
! t-wt
      uc = u%wt
      wc = 0.25D0*indxm*(dx(i-1)*w%c+dx(i)*w%w+dx(i-1)*w%t+dx(i)*w%wt)
      deltar  = dens%t*field%t - dens%wt*field%wt
      IF (k >= nz-1) deltar = 0.D0
      ft = ft - incrxm*(uc+ABS(uc))*(wc-ABS(wc))*deltar
!
! et-t
      uc = u%t
      wc = 0.25D0*indxp*(dx(i+1)*w%c+dx(i)*w%e+dx(i+1)*w%t+dx(i)*w%et)
      deltar  = dens%et*field%et - dens%t*field%t
      IF (k >= nz-1) deltar = 0.D0
      ft = ft - incrxp*(uc-ABS(uc))*(wc-ABS(wc))*deltar
!
! c-b
      uc = 0.25D0*indzm*(dz(k-1)*u%c+dz(k)*u%b+dz(k-1)*u%w+dz(k)*u%wb)
      wc = w%b
      deltar  = dens%c*field%c - dens%b*field%b
      fe = fe - incrzm*(uc+ABS(uc))*(wc+ABS(wc))*deltar*rb(i)
!
! t-c
      uc = 0.25D0*indzp*(dz(k+1)*u%c+dz(k)*u%t+dz(k+1)*u%w+dz(k)*u%wt)
      wc = w%c
      deltar  = dens%t*field%t - dens%c*field%c
      fe = fe - incrzp*(uc+ABS(uc))*(wc-ABS(wc))*deltar*rb(i)
!
! et-e
      uc = 0.25D0*indzp*(dz(k+1)*u%c+dz(k)*u%t+dz(k+1)*u%e+dz(k)*u%et)
      wc = w%e
      deltar  = dens%et*field%et - dens%e*field%e
      IF (i >= nx-1) deltar = 0.D0
      fe = fe - incrzp*(uc-ABS(uc))*(wc-ABS(wc))*deltar*rb(i)
!
! e-eb
      uc = 0.25D0*indzm*(dz(k-1)*u%c+dz(k)*u%b+dz(k-1)*u%e+dz(k)*u%eb)
      wc = w%eb
      deltar  = dens%e*field%e - dens%eb*field%eb
      IF (i >= nx-1) deltar = 0.D0
      fe = fe - incrzm*(uc-ABS(uc))*(wc+ABS(wc))*deltar*rb(i)
!
      RETURN
      END SUBROUTINE ctu1_fsc_2d
!-----------------------------------------------------------------------
      SUBROUTINE ctu2_fsc_2d(fe, ft, dens, field,  u, w, ij)
!
! ... Compute the second order Corner Transport Upwind correction (step 3 in LeVeque algorithm)
!
      USE dimensions
      USE domain_mapping, ONLY: myijk
      USE grid, ONLY: flag, rb, dx, dz, indx, indz, fluid
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: imjk, ijkm
      USE set_indexes, ONLY: stencil, ipjk, ijkp
      USE time_parameters, ONLY: dt
!
      IMPLICIT NONE
!
      REAL*8, INTENT(INOUT) :: fe, ft
      TYPE(stencil), INTENT(IN) :: dens, field, u, w
      INTEGER, INTENT(IN) :: ij
!
      INTEGER :: i,k,imesh, kp2, ip2, im2, km2
      REAL*8 :: uref, wref, deltar, lim, erre
      REAL*8 :: gradc, gradb, gradt, gradw, grade
      REAL*8 :: dxm, dxp, dxpp, dzm, dzp, dzpp, dxmm, dzmm
      REAL*8 :: indxm, indxp, indxpp, indzm, indzp, indzpp, indxmm, indzmm
!
      imesh = myijk( ip0_jp0_kp0_, ij)
      k = ( imesh - 1 ) / nx + 1
      i = MOD( ( imesh - 1 ), nx) + 1
!
      ip2 = MIN( nx, i+2 )
      kp2 = MIN( nz, k+2 )
      im2 = MAX(  1, i-2 )
      km2 = MAX(  1, k-2 )
!
      dxm  = dx(i)  + dx(i-1) 
      dxp  = dx(i)  + dx(i+1) 
      dxpp = dx(i+1)+ dx(ip2) 
      dxmm = dx(i-1)+ dx(im2)
      dzm  = dz(k)  + dz(k-1) 
      dzp  = dz(k)  + dz(k+1) 
      dzpp = dz(k+1)+ dz(kp2) 
      dzmm = dz(k-1)+ dz(km2) 
!
      indxm  = 2.D0/dxm 
      indxp  = 2.D0/dxp
      indxpp = 2.D0/dxpp
      indxmm = 2.D0/dxmm
      indzm  = 2.D0/dzm
      indzp  = 2.D0/dzp
      indzpp = 2.D0/dzpp
      indzmm = 2.D0/dzmm
!
! e-c
      gradc = indxp  * (dens%e *field%e  - dens%c*field%c)
      grade = indxpp * (dens%ee*field%ee - dens%e*field%e)
      gradw = indxm  * (dens%c *field%c  - dens%w*field%w)
      lim  = 0.D0
      erre = 0.D0
      uref = u%c
      IF (uref >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=nx-1 .AND. flag(ipjk)==fluid) CALL limiters(lm,lim,erre)
      uref = ABS(u%c)
      deltar = dens%e*field%e - dens%c*field%c
      fe = fe + 0.5D0*uref*(1-dt*indxp*uref)*deltar*lim*rb(i)
!
! t-c
      gradc = indzp  * (dens%t *field%t  - dens%c*field%c)
      gradt = indzpp * (dens%tt*field%tt - dens%t*field%t)
      gradb = indzm  * (dens%c *field%c  - dens%b*field%b)
      lim  = 0.D0
      erre = 0.D0
      wref = w%c
      IF (wref >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(ijkp)==fluid) CALL limiters(lm,lim,erre)
      wref = ABS(w%c)
      deltar = dens%t*field%t - dens%c*field%c
      ft = ft + 0.5D0*wref*(1-dt*indzp*wref)*deltar*lim
!
    RETURN
    END SUBROUTINE ctu2_fsc_2d
!-----------------------------------------------------------------------
    SUBROUTINE ctu3_fsc_2d(fe, ft, dens, field, u, w, ijk)
!
! ... Compute the second order Corner Transport Upwind correction (step 4 in LeVeque algorithm)
!
      USE dimensions
      USE domain_mapping, ONLY: myijk
      USE grid, ONLY: dx, dz, indx, indz, fluid, flag, rb
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: stencil, ipjk, imjk, ijkp, ijkm, ipjkm, ipjkp, imjkp
      USE time_parameters, ONLY: dt
      USE flux_limiters, ONLY: limiters, lm
!
      IMPLICIT NONE
!
      REAL*8, INTENT(INOUT) :: fe, ft
      TYPE(stencil), INTENT(IN) :: u, w, dens, field
      INTEGER, INTENT(IN) :: ijk
!
      INTEGER :: i, k, imesh
      REAL*8 :: uc, wc, deltar, lim, erre
      REAL*8 :: gradc, grade, gradw, gradt, gradb
      REAL*8 :: dxp, dxm, dzp, dzm, dxpp, dxmm, dzpp, dzmm
      REAL*8 :: indxp, indxm, indzp, indzm, indxpp, indxmm, indzpp, indzmm
      REAL*8 :: incrxp, incrxpp, incrzp, incrzpp
      INTEGER :: im2, ip2, km2, kp2
!
      imesh = myijk( ip0_jp0_kp0_, ijk)
      k = ( imesh - 1 ) / nx + 1
      i = MOD( ( imesh - 1 ), nx) + 1
!
      im2 = MAX(  1, i-2 )
      ip2 = MIN( nx, i+2 )
      km2 = MAX(  1, k-2 )
      kp2 = MIN( nz, k+2 )
!
      dxm  = dx(i)  + dx(i-1)
      dxp  = dx(i)  + dx(i+1)
      dxpp = dx(i+1)+ dx(ip2) 
      dxmm = dx(i-1)+ dx(im2)
      dzpp = dz(k+1)+ dz(kp2)
      dzmm = dz(k-1)+ dz(km2)
      dzm  = dz(k)  + dz(k-1)
      dzp  = dz(k)  + dz(k+1)
!
      indxm  = 2.D0/dxm
      indxp  = 2.D0/dxp
      indxpp = 2.D0/dxpp
      indxmm = 2.D0/dxmm
      indzp  = 2.D0/dzp
      indzm  = 2.D0/dzm
      indzpp = 2.D0/dzpp
      indzmm = 2.D0/dzmm
!
      incrzp  = 0.5D0*dt*indz(k)
      incrzpp = 0.5D0*dt*indz(k+1)
      incrxp  = 0.5D0*dt*indx(i)
      incrxpp = 0.5D0*dt*indx(i+1)
!
! c-w
      uc = u%w
      wc = 0.25D0*indxm*(dx(i-1)*w%c+dx(i)*w%w+dx(i-1)*w%b+dx(i)*w%wb)
      deltar = 0.5D0*ABS(uc)*(1 - dt*indxm*ABS(uc))*(dens%c*field%c - dens%w*field%w)
      gradc = indxm  * (dens%c*field%c - dens%w *field%w )
      grade = indxp  * (dens%e*field%e - dens%c *field%c )
      gradw = indxmm * (dens%w*field%w - dens%ww*field%ww)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=2 .AND. flag(imjk)==fluid) CALL limiters(lm,lim,erre)
      ft = ft + incrzp*(wc+ABS(wc))*deltar*lim
!
! e-c
      uc = u%c
      wc = 0.25D0*indxp*(dx(i+1)*w%c+dx(i)*w%e+dx(i+1)*w%b+dx(i)*w%eb)
      deltar = 0.5D0*ABS(uc)*(1 - dt*indxp*ABS(uc))*(dens%e*field%e - dens%c*field%c)
      gradc = indxp  * (dens%e *field%e  - dens%c*field%c)
      grade = indxpp * (dens%ee*field%ee - dens%e*field%e)
      gradw = indxm  * (dens%c *field%c  - dens%w*field%w)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=nx-1 .AND. flag(ipjk)==fluid) CALL limiters(lm,lim,erre)
      ft = ft - incrzp*(wc+ABS(wc))*deltar*lim
!
! t-wt
      uc = u%wt
      wc = 0.25D0*indxm*(dx(i-1)*w%c+dx(i)*w%w+dx(i-1)*w%t+dx(i)*w%wt)
      deltar = 0.5D0*ABS(uc)*(1 - dt*indxm*ABS(uc))*(dens%t*field%t - dens%wt*field%wt)
      IF (k >= nz-1) deltar = 0.D0
      gradc = indxm  * (dens%t *field%t  - dens%wt *field%wt )
      grade = indxp  * (dens%et*field%et - dens%t  *field%t  )
      gradw = indxmm * (dens%wt*field%wt - dens%wwt*field%wwt)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=2 .AND. flag(imjkp)==fluid) CALL limiters(lm,lim,erre)
      ft = ft + incrzpp*(wc-ABS(wc))*deltar*lim
!
! et-t
      uc = u%t
      wc = 0.25D0*indxp*(dx(i+1)*w%c+dx(i)*w%e+dx(i+1)*w%t+dx(i)*w%et)
      deltar = 0.5D0*ABS(uc)*(1 - dt*indxp*ABS(uc))*(dens%et*field%et - dens%t*field%t)
      IF (k >= nz-1) deltar = 0.D0
      gradc = indxp  * (dens%et *field%et  - dens%t *field%t )
      grade = indxpp * (dens%eet*field%eet - dens%et*field%et)
      gradw = indxm  * (dens%t  *field%t   - dens%wt*field%wt)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=nx-1 .AND. flag(ipjkp)==fluid) CALL limiters(lm,lim,erre)
      ft = ft - incrzpp*(wc-ABS(wc))*deltar*lim
!
! c-b
      uc = 0.25D0*indzm*(dz(k-1)*u%c+dz(k)*u%b+dz(k-1)*u%w+dz(k)*u%wb)
      wc = w%b
      deltar = 0.5D0*ABS(wc)*(1 - dt*indzm*ABS(wc))*(dens%c*field%c - dens%b*field%b)
      gradc = indzm  * (dens%c*field%c - dens%b *field%b )
      gradt = indzp  * (dens%t*field%t - dens%c *field%c )
      gradb = indzmm * (dens%b*field%b - dens%bb*field%bb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=2 .AND. flag(ijkm)==fluid) CALL limiters(lm,lim,erre)
      fe = fe + incrxp*(uc+ABS(uc))*deltar*lim*rb(i)
!
! t-c
      uc = 0.25D0*indzp*(dz(k+1)*u%c+dz(k)*u%t+dz(k+1)*u%w+dz(k)*u%wt)
      wc = w%c
      deltar = 0.5D0*ABS(wc)*(1 - dt*indzp*ABS(wc))*(dens%t*field%t - dens%c*field%c)
      gradc = indzp  * (dens%t *field%t  - dens%c*field%c)
      gradt = indzpp * (dens%tt*field%tt - dens%t*field%t)
      gradb = indzm  * (dens%c *field%c  - dens%b*field%b)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(ijkp)==fluid) CALL limiters(lm,lim,erre)
      fe = fe - incrxp*(uc+ABS(uc))*deltar*lim*rb(i)
!
! et-e 
      uc = 0.25D0*indzp*(dz(k+1)*u%c+dz(k)*u%t+dz(k+1)*u%e+dz(k)*u%et)
      wc = w%e
      deltar = 0.5D0*ABS(wc)*(1 - dt*indzp*ABS(wc))*(dens%et*field%et - dens%e*field%e)
      IF (i >= nx-1) deltar = 0.D0
      gradc = indzp  * (dens%et *field%et  - dens%e *field%e )
      gradt = indzpp * (dens%ett*field%ett - dens%et*field%et)
      gradb = indzm  * (dens%e  *field%e   - dens%eb*field%eb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(ipjkp)==fluid) CALL limiters(lm,lim,erre)
      fe = fe - incrxpp*(uc-ABS(uc))*deltar*lim*rb(i)
!
! e-eb
      uc = 0.25D0*indzm*(dz(k-1)*u%c+dz(k)*u%b+dz(k-1)*u%e+dz(k)*u%eb)
      wc = w%eb
      deltar = 0.5D0*ABS(wc)*(1 - dt*indzm*ABS(wc))*(dens%e*field%e - dens%eb*field%eb)
      IF (i >= nx-1) deltar = 0.D0
      gradc = indzm  * (dens%e *field%e  - dens%eb *field%eb )
      gradt = indzp  * (dens%et*field%et - dens%e  *field%e  )
      gradb = indzmm * (dens%eb*field%eb - dens%ebb*field%ebb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=2 .AND. flag(ipjkm)==fluid) CALL limiters(lm,lim,erre)
      fe = fe + incrxpp*(uc-ABS(uc))*deltar*lim*rb(i)
!
      RETURN
      END SUBROUTINE ctu3_fsc_2d
!-----------------------------------------------------------------------
    SUBROUTINE ctu1_fsc_3d( fe, fn, ft, dens, field, u, v, w, ijk )
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
      TYPE(stencil), INTENT(IN) :: u, v, w, dens, field
      INTEGER, INTENT(IN) :: ijk
!
      INTEGER :: i, j, k, imesh
      REAL*8 :: uc, vc, wc, deltar, onethird
      REAL*8 :: dxp, dxm, dyp, dym, dzp, dzm
      REAL*8 :: dxi, dxip, dxim, dyj, dyjp, dyjm, dzk, dzkp, dzkm
      REAL*8 :: indxp, indxm, indyp, indym, indzp, indzm
      REAL*8 :: incrxp, incrxm, incryp, incrym, incrzp, incrzm
!
      imesh = myijk( ip0_jp0_kp0_, ijk)
      CALL meshinds(ijk, imesh, i, j, k)
      onethird = 1.D0/3.D0
!
      dxi  = dx(i)
      dxip = dx(i+1)
      dxim = dx(i-1)
      dyj  = dy(j)
      dyjp = dy(j+1)
      dyjm = dy(j-1)
      dzk  = dz(k)
      dzkp = dz(k+1)
      dzkm = dz(k-1)
!
      dxm = dxi + dxim
      dxp = dxi + dxip
      dym = dyj + dyjm
      dyp = dyj + dyjp
      dzm = dzk + dzkm
      dzp = dzk + dzkp
!
      indxp = 2.D0/dxp
      indxm = 2.D0/dxm
      indyp = 2.D0/dyp
      indym = 2.D0/dym
      indzp = 2.D0/dzp
      indzm = 2.D0/dzm
!
      incrxp = 0.25D0*dt*indxp
      incrxm = 0.25D0*dt*indxm
      incryp = 0.25D0*dt*indyp
      incrym = 0.25D0*dt*indym
      incrzp = 0.25D0*dt*indzp
      incrzm = 0.25D0*dt*indzm
!
! t-c
      deltar  = dens%t*field%t - dens%c*field%c
      uc = (dzkp*u%c+dzk*u%t+dzkp*u%w+dzk*u%wt)*indx(i)*incrzp
      vc = (dzkp*v%c+dzk*v%t+dzkp*v%s+dzk*v%st)*indy(j)*incrzp
      wc = w%c
      fe = fe - 0.125D0*(uc+ABS(uc))*(wc-ABS(wc))*deltar
      fn = fn - 0.125D0*(vc+ABS(vc))*(wc-ABS(wc))*deltar &
              + 0.25D0*onethird*(vc+ABS(vc))*(wc-ABS(wc))*ABS(uc)*deltar 
!
! c-b
      deltar  = dens%c*field%c - dens%b*field%b
      uc = (dzkm*u%c+dzk*u%b+dzkm*u%w+dzk*u%wb)*indx(i)*incrzm
      vc = (dzkm*v%c+dzk*v%b+dzkm*v%s+dzk*v%sb)*indy(j)*incrzm
      wc = w%b
      fe = fe - 0.125D0*(uc+ABS(uc))*(wc+ABS(wc))*deltar
      fn = fn - 0.125D0*(vc+ABS(vc))*(wc+ABS(wc))*deltar &
              + 0.25D0*onethird*(vc+ABS(vc))*(wc+ABS(wc))*ABS(uc)*deltar 
!
! n-c
      deltar  = dens%n*field%n - dens%c*field%c
      uc = (dyjp*u%c+dyj*u%n+dyjp*u%w+dyj*u%wn)*indx(i)*incryp
      vc = v%c
      wc = (dyjp*w%c+dyj*w%n+dyjp*w%b+dyj*w%nb)*indz(k)*incryp
      ft = ft - 0.125D0*(wc+ABS(wc))*(vc-ABS(vc))*deltar
      fe = fe - 0.125D0*(uc+ABS(uc))*(vc-ABS(vc))*deltar &
              + 0.25D0*onethird*(uc+ABS(uc))*(vc-ABS(vc))*ABS(wc)*deltar 
!
! c-s
      deltar  = dens%c*field%c - dens%s*field%s
      uc = (dyjm*u%c+dyj*u%s+dyjm*u%w+dyj*u%ws)*indx(i)*incrym
      vc = v%s
      wc = (dyjm*w%c+dyj*w%s+dyjm*w%b+dyj*w%sb)*indz(k)*incrym
      ft = ft - 0.125D0*(wc+ABS(wc))*(vc+ABS(vc))*deltar
      fe = fe - 0.125D0*(uc+ABS(uc))*(vc+ABS(vc))*deltar &
              + 0.25D0*onethird*(uc+ABS(uc))*(vc+ABS(vc))*ABS(wc)*deltar 
!
! e-c
      deltar  = dens%e*field%e - dens%c*field%c
      uc = u%c
      vc = (dxip*v%c+dxi*v%e+dxip*v%s+dxi*v%es)*indy(j)*incrxp
      wc = (dxip*w%c+dxi*w%e+dxip*w%b+dxi*w%eb)*indz(k)*incrxp
      fn = fn - 0.125D0*(vc+ABS(vc))*(uc-ABS(uc))*deltar
      ft = ft - 0.125D0*(wc+ABS(wc))*(uc-ABS(uc))*deltar &
              + 0.25D0*onethird*(wc+ABS(wc))*(uc-ABS(uc))*ABS(vc)*deltar 
!
! c-w
      deltar  = dens%c*field%c - dens%w*field%w
      uc = u%w
      vc = (dxim*v%c+dxi*v%w+dxim*v%s+dxi*v%ws)*indy(j)*incrxm
      wc = (dxim*w%c+dxi*w%w+dxim*w%b+dxi*w%wb)*indz(k)*incrxm
      fn = fn - 0.125D0*(vc+ABS(vc))*(uc+ABS(uc))*deltar
      ft = ft - 0.125D0*(wc+ABS(wc))*(uc+ABS(uc))*deltar &
              + 0.25D0*onethird*(wc+ABS(wc))*(uc+ABS(uc))*ABS(vc)*deltar 
!
! wt-w
      deltar  = dens%wt*field%wt - dens%w*field%w
      IF (i <= 2) deltar = 0.D0
      uc = (dzkp*u%w+dzk*u%wt+dzkp*u%ww+dzk*u%wwt)*indx(i-1)*incrzp
      vc = (dzkp*v%w+dzk*v%wt+dzkp*v%ws+dzk*v%wst)*indy(j)  *incrzp
      wc = w%w
      fn = fn - 0.125D0*onethird*(wc-ABS(wc))*(vc+ABS(vc))*ABS(uc+ABS(uc))*deltar
!
! w-wb
      deltar  = dens%w*field%w - dens%wb*field%wb
      IF (i <= 2) deltar = 0.D0
      uc = (dzkm*u%w+dzk*u%wb+dzkm*u%ww+dzk*u%wwb)*indx(i-1)*incrzm
      vc = (dzkm*v%w+dzk*v%wb+dzkm*v%ws+dzk*v%wsb)*indy(j)  *incrzm
      wc = w%wb
      fn = fn - 0.125D0*onethird*(wc+ABS(wc))*(vc+ABS(vc))*ABS(uc+ABS(uc))*deltar
!
! et-e
      deltar  = dens%et*field%et - dens%e*field%e
      IF (i >= nx-1) deltar = 0.D0
      uc = (dzkp*u%e+dzk*u%et+dzkp*u%c +dzk*u%t  )*indx(i+1)*incrzp
      vc = (dzkp*v%e+dzk*v%et+dzkp*v%es+dzk*v%est)*indy(j)  *incrzp
      wc = w%e
      fe = fe - 0.125D0*(wc-ABS(wc))*(uc-ABS(uc))*deltar
      fn = fn - 0.125D0*onethird*(wc-ABS(wc))*(vc+ABS(vc))*ABS(uc-ABS(uc))*deltar
!
! e-eb
      deltar  = dens%e*field%e - dens%eb*field%eb
      IF (i >= nx-1) deltar = 0.D0
      uc = (dzkm*u%e+dzk*u%eb+dzkm*u%c +dzk*u%b  )*indx(i+1)*incrzm
      vc = (dzkm*v%e+dzk*v%eb+dzkm*v%es+dzk*v%esb)*indy(j)  *incrzm
      wc = w%eb
      fe = fe - 0.125D0*(wc+ABS(wc))*(uc-ABS(uc))*deltar
      fn = fn - 0.125D0*onethird*(wc+ABS(wc))*(vc+ABS(vc))*ABS(uc-ABS(uc))*deltar
!
! en-e
      deltar  = dens%en*field%en - dens%e*field%e
      IF (i >= nx-1) deltar = 0.D0
      uc = (dyjp*u%e+dyj*u%en+dyjp*u%c +dyj*u%n  )*indx(i+1)*incryp
      vc = v%e
      wc = (dyjp*w%e+dyj*w%en+dyjp*w%eb+dyj*w%enb)*indz(k)  *incryp
      fe = fe - 0.125D0*(vc-ABS(vc))*(uc-ABS(uc))*deltar &
              + 0.25D0*onethird*(vc-ABS(vc))*(uc-ABS(uc))*ABS(wc)*deltar
!
! e-es
      deltar  = dens%e*field%e - dens%es*field%es
      IF (i >= nx-1) deltar = 0.D0
      uc = (dyjm*u%e+dyj*u%es+dyjm*u%c +dyj*u%s  )*indx(i+1)*incrym
      vc = v%es
      wc = (dyjm*w%e+dyj*w%es+dyjm*w%eb+dyj*w%esb)*indz(k)  *incrym
      fe = fe - 0.125D0*(vc+ABS(vc))*(uc-ABS(uc))*deltar &
              + 0.25D0*onethird*(vc+ABS(vc))*(uc-ABS(uc))*ABS(wc)*deltar
!
! nt-n
      deltar  = dens%nt*field%nt - dens%n*field%n
      IF (j >= ny-1) deltar = 0.D0
      uc = (dzkp*u%n+dzk*u%nt+dzkp*u%wn+dzk*u%wnt)*indx(i)  *incrzp
      vc = (dzkp*v%n+dzk*v%nt+dzkp*v%c +dzk*v%t  )*indy(j+1)*incrzp
      wc = w%n
      fn = fn - 0.125D0*(vc-ABS(vc))*(wc-ABS(wc))*deltar &
              + 0.25D0*onethird*(vc-ABS(vc))*(wc-ABS(wc))*ABS(uc)*deltar
!
! n-nb
      deltar  = dens%n*field%n - dens%nb*field%nb
      IF (j >= ny-1) deltar = 0.D0
      uc = (dzkm*u%n+dzk*u%nb+dzkm*u%wn+dzk*u%wnb)*indx(i)  *incrzm
      vc = (dzkm*v%n+dzk*v%nb+dzkm*v%c +dzk*v%b  )*indy(j+1)*incrzm
      wc = w%nb
      fn = fn - 0.125D0*(vc-ABS(vc))*(wc+ABS(wc))*deltar &
              + 0.25D0*onethird*(vc-ABS(vc))*(wc+ABS(wc))*ABS(uc)*deltar
!
! nt-t 
      deltar  = dens%nt*field%nt - dens%t*field%t
      IF (k >= nz-1) deltar = 0.D0
      uc = (dyjp*u%t+dyj*u%nt+dyjp*u%wt+dyj*u%wnt)*indx(i)  *incryp
      vc = v%t
      wc = (dyjp*w%t+dyj*w%nt+dyjp*w%c +dyj*w%n  )*indz(k+1)*incryp
      ft = ft - 0.125D0*(vc-ABS(vc))*(wc-ABS(wc))*deltar
      fe = fe - 0.125D0*onethird*(vc-ABS(vc))*(uc+ABS(uc))*ABS(wc-ABS(wc))*deltar
!     
! t-st
      deltar  = dens%t*field%t - dens%st*field%st
      IF (k >= nz-1) deltar = 0.D0
      uc = (dyjm*u%t+dyj*u%st+dyjm*u%wt+dyj*u%wst)*indx(i)  *incrym
      vc = v%st
      wc = (dyjm*w%t+dyj*w%st+dyjm*w%c +dyj*w%s  )*indz(k+1)*incrym
      ft = ft - 0.125D0*(vc+ABS(vc))*(wc-ABS(wc))*deltar
      fe = fe - 0.125D0*onethird*(vc+ABS(vc))*(uc+ABS(uc))*ABS(wc-ABS(wc))*deltar
!     
! nb-b
      deltar  = dens%nb*field%nb - dens%b*field%b
      IF (k <= 2) deltar = 0.D0
      uc = (dyjp*u%b+dyj*u%nb+dyjp*u%wb+dyj*u%wnb)*indx(i)  *incryp
      vc = v%b
      wc = (dyjp*w%b+dyj*w%nb+dyjp*w%bb+dyj*w%nbb)*indz(k-1)*incryp
      fe = fe - 0.125D0*onethird*(vc-ABS(vc))*(uc+ABS(uc))*ABS(wc+ABS(wc))*deltar
!     
! b-sb
      deltar  = dens%b*field%b - dens%sb*field%sb
      IF (k <= 2) deltar = 0.D0
      uc = (dyjm*u%b+dyj*u%sb+dyjm*u%wb+dyj*u%wsb)*indx(i)  *incrym
      vc = v%sb
      wc = (dyjm*w%b+dyj*w%sb+dyjm*w%bb+dyj*w%sbb)*indz(k-1)*incrym
      fe = fe - 0.125D0*onethird*(vc+ABS(vc))*(uc+ABS(uc))*ABS(wc+ABS(wc))*deltar
!     
! en-n
      deltar  = dens%en*field%en - dens%n*field%n
      IF (j >= ny-1) deltar = 0.D0
      uc = u%n
      vc = (dxip*v%n+dxi*v%en+dxip*v%c +dxi*v%e  )*indy(j+1)*incrxp
      wc = (dxip*w%n+dxi*w%en+dxip*w%nb+dxi*w%enb)*indz(k)  *incrxp
      fn = fn - 0.125D0*(vc-ABS(vc))*(uc-ABS(uc))*deltar
      ft = ft - 0.125D0*onethird*(uc-ABS(uc))*(wc+ABS(wc))*ABS(vc-ABS(vc))*deltar
!     
! n-wn
      deltar  = dens%n*field%n - dens%wn*field%wn
      IF (j >= ny-1) deltar = 0.D0
      uc = u%wn
      vc = (dxim*v%n+dxi*v%wn+dxim*v%c +dxi*v%w  )*indy(j+1)*incrxm
      wc = (dxim*w%n+dxi*w%wn+dxim*w%nb+dxi*w%wnb)*indz(k)  *incrxm
      fn = fn - 0.125D0*(vc-ABS(vc))*(uc+ABS(uc))*deltar
      ft = ft - 0.125D0*onethird*(uc+ABS(uc))*(wc+ABS(wc))*ABS(vc-ABS(vc))*deltar
!     
! es-s
      deltar  = dens%es*field%es - dens%s*field%s
      IF (j <= 2) deltar = 0.D0
      uc = u%s
      vc = (dxip*v%s+dxi*v%es+dxip*v%ss+dxi*v%ess)*indy(j-1)*incrxp
      wc = (dxip*w%s+dxi*w%es+dxip*w%sb+dxi*w%esb)*indz(k)  *incrxp
      ft = ft - 0.125D0*onethird*(uc-ABS(uc))*(wc+ABS(wc))*ABS(vc+ABS(vc))*deltar
!     
! s-ws
      deltar  = dens%s*field%s - dens%ws*field%ws
      IF (j <= 2) deltar = 0.D0
      uc = u%ws
      vc = (dxim*v%s+dxi*v%ws+dxim*v%ss+dxi*v%wss)*indy(j-1)*incrxm
      wc = (dxim*w%s+dxi*w%ws+dxim*w%sb+dxi*w%wsb)*indz(k)  *incrxm
      ft = ft - 0.125D0*onethird*(uc+ABS(uc))*(wc+ABS(wc))*ABS(vc+ABS(vc))*deltar
!     
! et-t
      deltar  = dens%et*field%et - dens%t*field%t
      IF (k >= nz-1) deltar = 0.D0
      uc = u%t
      vc = (dxip*v%t+dxi*v%et+dxip*v%st+dxi*v%est)*indy(j)  *incrxp
      wc = (dxip*w%t+dxi*w%et+dxip*w%c +dxi*w%e  )*indz(k+1)*incrxp
      ft = ft - 0.125D0*(wc-ABS(wc))*(uc-ABS(uc))*deltar &
              + 0.25D0*onethird*(wc-ABS(wc))*(uc-ABS(uc))*ABS(vc)*deltar
!     
! t-wt
      deltar  = dens%t*field%t - dens%wt*field%wt
      IF (k >= nz-1) deltar = 0.D0
      uc = u%wt
      vc = (dxim*v%t+dxi*v%wt+dxim*v%st+dxi*v%wst)*indy(j)  *incrxm
      wc = (dxim*w%t+dxi*w%wt+dxim*w%c +dxi*w%w  )*indz(k+1)*incrxm
      ft = ft - 0.125D0*(wc-ABS(wc))*(uc+ABS(uc))*deltar &
              + 0.25D0*onethird*(wc-ABS(wc))*(uc+ABS(uc))*ABS(vc)*deltar
!     
! wnt-wn
      deltar  = dens%wnt*field%wnt - dens%wn*field%wn
      IF (i <= 2) deltar = 0.D0
      IF (j >= ny-1) deltar = 0.D0
      wc = w%wn
      uc = (dzkp*u%wn+dzk*u%wnt+dzkp*u%wwn+dzk*u%wwnt)*indx(i-1)*incrzp
      vc = (dzkp*v%wn+dzk*v%wnt+dzkp*v%w  +dzk*v%wt  )*indy(j+1)*incrzp
      fn = fn - 0.125D0*onethird*(wc-ABS(wc))*(vc-ABS(vc))*ABS(uc+ABS(uc))*deltar
!     
! wn-wnb
      deltar  = dens%wn*field%wn - dens%wnb*field%wnb
      IF (i <= 2) deltar = 0.D0
      IF (j >= ny-1) deltar = 0.D0
      uc = (dzkm*u%wn+dzk*u%wnb+dzkm*u%wwn+dzk*u%wwnb)*indx(i-1)*incrzm
      vc = (dzkm*v%wn+dzk*v%wnb+dzkm*v%w  +dzk*v%wb  )*indy(j+1)*incrzm
      wc = w%wnb
      fn = fn - 0.125D0*onethird*(wc+ABS(wc))*(vc-ABS(vc))*ABS(uc+ABS(uc))*deltar
!     
! ent-et
      deltar  = dens%ent*field%ent - dens%et*field%et
      uc = (dyjp*u%et+dyj*u%ent+dyjp*u%t+dyj*u%nt)*indx(i+1)*incryp
      vc = v%et
      wc = (dyjp*w%et+dyj*w%ent+dyjp*w%e+dyj*w%en)*indz(k+1)*incryp
      fe = fe - 0.125D0*onethird*(vc-ABS(vc))*(uc-ABS(uc))*ABS(wc-ABS(wc))*deltar
!     
! et-est
      deltar  = dens%et*field%et - dens%est*field%est
      IF (i >= nx-1) deltar = 0.D0
      IF (k >= nz-1) deltar = 0.D0
      uc = (dyjm*u%et+dyj*u%est+dyjm*u%t+dyj*u%st)*indx(i+1)*incrym
      vc = v%est
      wc = (dyjm*w%et+dyj*w%est+dyjm*w%e+dyj*w%es)*indz(k+1)*incrym
      fe = fe - 0.125D0*onethird*(vc+ABS(vc))*(uc-ABS(uc))*ABS(wc-ABS(wc))*deltar
!     
! enb-eb
      deltar  = dens%enb*field%enb - dens%eb*field%eb
      IF (i >= nx-1) deltar = 0.D0
      IF (k <= 2) deltar = 0.D0
      uc = (dyjp*u%eb+dyj*u%enb+dyjp*u%b  +dyj*u%nb  )*indx(i+1)*incryp
      vc = v%eb
      wc = (dyjp*w%eb+dyj*w%enb+dyjp*w%ebb+dyj*w%enbb)*indz(k-1)*incryp
      fe = fe - 0.125D0*onethird*(vc-ABS(vc))*(uc-ABS(uc))*ABS(wc+ABS(wc))*deltar
!     
! eb-esb
      deltar  = dens%eb*field%eb - dens%esb*field%esb
      IF (i >= nx-1) deltar = 0.D0
      IF (k <= 2) deltar = 0.D0
      uc = (dyjm*u%eb+dyj*u%esb+dyjm*u%b  +dyj*u%sb  )*indx(i+1)*incrym
      vc = v%esb
      wc = (dyjm*w%eb+dyj*w%esb+dyjm*w%ebb+dyj*w%esbb)*indz(k-1)*incrym
      fe = fe - 0.125D0*onethird*(vc+ABS(vc))*(uc-ABS(uc))*ABS(wc+ABS(wc))*deltar
!     
! ent-en
      deltar  = dens%ent*field%ent - dens%en*field%en
      IF (i >= nx-1) deltar = 0.D0
      IF (j >= ny-1) deltar = 0.D0
      uc = (dzkp*u%en+dzk*u%ent+dzkp*u%n+dzk*u%nt)*indx(i+1)*incrzp
      vc = (dzkp*v%en+dzk*v%ent+dzkp*v%e+dzk*v%et)*indy(j+1)*incrzp
      wc = w%en
      fn = fn - 0.125D0*onethird*(wc-ABS(wc))*(vc-ABS(vc))*ABS(uc-ABS(uc))*deltar
!     
! en-enb
      deltar  = dens%en*field%en - dens%enb*field%enb
      IF (i >= nx-1) deltar = 0.D0
      IF (j >= ny-1) deltar = 0.D0
      uc = (dzkm*u%en+dzk*u%enb+dzkm*u%n+dzk*u%nb)*indx(i+1)*incrzm
      vc = (dzkm*v%en+dzk*v%enb+dzkm*v%e+dzk*v%eb)*indy(j+1)*incrzm
      wc = w%enb
      fn = fn - 0.125D0*onethird*(wc+ABS(wc))*(vc-ABS(vc))*ABS(uc-ABS(uc))*deltar
!     
! ent-nt
      deltar  = dens%ent*field%ent - dens%nt*field%nt
      IF (k >= nz-1) deltar = 0.D0
      IF (j >= ny-1) deltar = 0.D0
      uc = u%nt
      vc = (dxip*v%nt+dxi*v%ent+dxip*v%t+dxi*v%et)*indy(j+1)*incrxp
      wc = (dxip*w%nt+dxi*w%ent+dxip*w%n+dxi*w%en)*indz(k+1)*incrxp
      ft = ft - 0.125D0*onethird*(wc-ABS(wc))*(uc-ABS(uc))*ABS(vc-ABS(vc))*deltar
!     
! nt-wnt
      deltar  = dens%nt*field%nt - dens%wnt*field%wnt
      IF (k >= nz-1) deltar = 0.D0
      IF (j >= ny-1) deltar = 0.D0
      uc = u%wnt
      vc = (dxim*v%nt+dxi*v%wnt+dxim*v%t+dxi*v%wt)*indy(j+1)*incrxm
      wc = (dxim*w%nt+dxi*w%wnt+dxim*w%n+dxi*w%wn)*indz(k+1)*incrxm
      ft = ft - 0.125D0*onethird*(wc-ABS(wc))*(uc+ABS(uc))*ABS(vc-ABS(vc))*deltar
!     
! est-st
      deltar  = dens%est*field%est - dens%st*field%st
      IF (k >= nz-1) deltar = 0.D0
      IF (j <= 2) deltar = 0.D0
      uc = u%st
      vc = (dxip*v%st+dxi*v%est+dxip*v%sst+dxi*v%esst)*indy(j-1)*incrxp
      wc = (dxip*w%st+dxi*w%est+dxip*w%s  +dxi*w%es  )*indz(k+1)*incrxp
      ft = ft - 0.125D0*onethird*(wc-ABS(wc))*(uc-ABS(uc))*ABS(vc+ABS(vc))*deltar
!     
! st-wst
      deltar  = dens%st*field%st - dens%wst*field%wst
      IF (k >= nz-1) deltar = 0.D0
      IF (j <= 2) deltar = 0.D0
      uc = u%wst
      vc = (dxim*v%st+dxi*v%wst+dxim*v%sst+dxi*v%wsst)*indy(j-1)*incrxm
      wc = (dxim*w%st+dxi*w%wst+dxim*w%s  +dxi*w%ws  )*indz(k+1)*incrxm
      ft = ft - 0.125D0*onethird*(wc-ABS(wc))*(uc+ABS(uc))*ABS(vc+ABS(vc))*deltar
!
      RETURN
      END SUBROUTINE ctu1_fsc_3d
!----------------------------------------------------------------------
      SUBROUTINE ctu2_fsc_3d(fe, fn, ft, dens, field, u, v, w, ijk)
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
!
      IMPLICIT NONE

      REAL*8, INTENT(INOUT) :: fe, ft, fn
      TYPE(stencil), INTENT(IN) :: dens, field, u, v, w 
      INTEGER, INTENT(IN) :: ijk 
!
      INTEGER :: im2, ip2, jm2, jp2, km2, kp2
      REAL*8 :: uref, vref, wref, deltar, erre, lim
      INTEGER :: i, j, k, imesh
      REAL*8 :: dxp, dxm, dxpp, dxmm, indxpp, indxp, indxm, indxmm
      REAL*8 :: dyp, dym, dypp, dymm, indypp, indyp, indym, indymm
      REAL*8 :: dzp, dzm, dzpp, dzmm, indzpp, indzp, indzm, indzmm
      REAL*8 :: dxi, dxip, dxipp, dxim, dximm
      REAL*8 :: dyj, dyjp, dyjpp, dyjm, dyjmm
      REAL*8 :: dzk, dzkp, dzkpp, dzkm, dzkmm
      REAL*8 :: gradc, grade, gradw, gradn, grads, gradt, gradb
!
      imesh = myijk( ip0_jp0_kp0_, ijk)
      i = MOD( MOD( imesh - 1, nx*ny ), nx ) + 1
      j = MOD( imesh - 1, nx*ny ) / nx + 1
      k = ( imesh - 1 ) / ( nx*ny ) + 1 
!
      im2 = MAX(  1, i-2 )
      ip2 = MIN( nx, i+2 )
      jm2 = MAX(  1, j-2 )
      jp2 = MIN( ny, j+2 )
      km2 = MAX(  1, k-2 )
      kp2 = MIN( nz, k+2 )
!
      dxi   = dx(i) 
      dxip  = dx(i+1)
      dxipp = dx(ip2)
      dxim  = dx(i-1)
      dximm = dx(im2)
      dyj   = dy(j) 
      dyjp  = dy(j+1)
      dyjpp = dy(jp2)
      dyjm  = dy(j-1)
      dyjmm = dy(jm2)
      dzk   = dz(k) 
      dzkp  = dz(k+1)
      dzkpp = dz(kp2)
      dzkm  = dz(k-1)
      dzkmm = dz(km2)
!
      dxm  = dxi  + dxim
      dxp  = dxi  + dxip
      dxpp = dxip + dxipp
      dxmm = dxim + dximm
      dym  = dyj  + dyjm
      dyp  = dyj  + dyjp
      dypp = dyjp + dyjpp
      dymm = dyjm + dyjmm
      dzm  = dzk  + dzkm
      dzp  = dzk  + dzkp
      dzpp = dzkp + dzkpp
      dzmm = dzkm + dzkmm

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
! e-c
      gradc = indxp  * (dens%e *field%e  - dens%c*field%c)
      grade = indxpp * (dens%ee*field%ee - dens%e*field%e)
      gradw = indxm  * (dens%c *field%c  - dens%w*field%w)
      lim  = 0.D0
      erre = 0.D0
      uref = u%c
      IF (uref >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=nx-1 .AND. flag(ipjk)==fluid) CALL limiters(lm,lim,erre)
      uref = ABS(u%c)
      deltar = dens%e*field%e - dens%c*field%c
      fe = fe + 0.5D0*uref*(1-dt*indxp*uref)*deltar*lim
!
! n-c
      gradc = indyp  * (dens%n *field%n  - dens%c*field%c)
      gradn = indypp * (dens%nn*field%nn - dens%n*field%n)
      grads = indym  * (dens%c *field%c  - dens%s*field%s)
      lim  = 0.D0
      erre = 0.D0
      vref = v%c
      IF (vref >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = grads/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradn/gradc
      END IF
      IF (j/=ny-1 .AND. flag(ijpk)==fluid) CALL limiters(lm,lim,erre)
      vref = ABS(v%c)
      deltar = dens%n*field%n - dens%c*field%c
      fn = fn + 0.5D0*vref*(1-dt*indyp*vref)*deltar*lim
!
! t-c
      gradc = indzp  * (dens%t *field%t  - dens%c*field%c)
      gradt = indzpp * (dens%tt*field%tt - dens%t*field%t)
      gradb = indzm  * (dens%c *field%c  - dens%b*field%b)
      lim  = 0.D0
      erre = 0.D0
      wref = w%c
      IF (wref >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(ijkp)==fluid) CALL limiters(lm,lim,erre)
      wref = ABS(w%c)
      deltar = dens%t*field%t - dens%c*field%c
      ft = ft + 0.5D0*wref*(1-dt*indzp*wref)*deltar*lim
!
      RETURN
      END SUBROUTINE ctu2_fsc_3d
!----------------------------------------------------------------------
      SUBROUTINE ctu3_fsc_3d(fe, fn, ft, dens, field, u, v, w, ijk)
!
! ... Compute the second order Corner Transport Upwind correction (step 4 in LeVeque algorithm)
!
      USE dimensions
      USE domain_mapping, ONLY: myijk, meshinds
      USE grid, ONLY: dx, dy, dz, indx, indy, indz, fluid, flag
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: stencil
      USE set_indexes, ONLY: imjk, ijkm, ipjk, ijkp, ijmk, ijpk
      USE set_indexes, ONLY: ipjpkp, ipjpkm, ipjmk, ipjmkm, ipjpk, ipjkp, imjpkm
      USE set_indexes, ONLY: ipjmkp, ipjkm, ipjmkp, imjpk, imjpkp, imjmkp, imjkm
      USE set_indexes, ONLY: imjmk, imjkp, ijpkp, ijpkm, ijmkm, ijmkp
      USE time_parameters, ONLY: dt

      IMPLICIT NONE

      REAL*8, INTENT(INOUT) :: fe, ft, fn
      TYPE(stencil), INTENT(IN) :: dens, field, u, v, w 
      INTEGER, INTENT(IN) :: ijk 
!
      INTEGER :: im2, ip2, jm2, jp2, km2, kp2
      REAL*8 :: uc, vc, wc, deltar, erre, lim, s
      INTEGER :: i, j, k, imesh
      REAL*8 :: dxp, dxm, dxpp, dxmm, indxpp, indxp, indxm, indxmm
      REAL*8 :: dyp, dym, dypp, dymm, indypp, indyp, indym, indymm
      REAL*8 :: dzp, dzm, dzpp, dzmm, indzpp, indzp, indzm, indzmm
      REAL*8 :: incrxp, incrxm, incryp, incrym, incrzp, incrzm
      REAL*8 :: gradc, grade, gradw, gradn, grads, gradt, gradb
!
      CALL meshinds(ijk, imesh, i, j, k)
!
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
      indxm  = 2.0D0/dxm
      indxp  = 2.0D0/dxp
      indxpp = 2.0D0/dxpp
      indxmm = 2.0D0/dxmm
      indym  = 2.0D0/dym
      indyp  = 2.0D0/dyp
      indypp = 2.0D0/dypp
      indymm = 2.0D0/dymm
      indzm  = 2.0D0/dzm
      indzp  = 2.0D0/dzp
      indzpp = 2.0D0/dzpp
      indzmm = 2.0D0/dzmm
!
      incrxp = 0.25D0*dt*indxp
      incrxm = 0.25D0*dt*indxm
      incryp = 0.25D0*dt*indyp
      incrym = 0.25D0*dt*indym
      incrzp = 0.25D0*dt*indzp
      incrzm = 0.25D0*dt*indzm
!
! c-w
      deltar  = dens%c*field%c - dens%w*field%w
      uc = u%w
      vc = (dx(i-1)*v%c+dx(i)*v%w+dx(i-1)*v%s+dx(i)*v%ws)*indy(j)*incrxm
      wc = (dx(i-1)*w%c+dx(i)*w%w+dx(i-1)*w%b+dx(i)*w%wb)*indz(k)*incrxm
      gradc = indxm  * (dens%c*field%c - dens%w *field%w )
      grade = indxp  * (dens%e*field%e - dens%c *field%c )
      gradw = indxmm * (dens%w*field%w - dens%ww*field%ww)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=2 .AND. flag(imjk)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxm*ABS(uc))*deltar*lim
      fn = fn + 0.5D0*(vc+ABS(vc))*s
      ft = ft + 0.5D0*(1-ABS(vc))*(wc+ABS(wc))*s
!
! e-c
      deltar  = dens%e*field%e - dens%c*field%c
      uc = u%c
      vc = (dx(i+1)*v%c+dx(i)*v%e+dx(i+1)*v%s+dx(i)*v%es)*indy(j)*incrxp
      wc = (dx(i+1)*w%c+dx(i)*w%e+dx(i+1)*w%b+dx(i)*w%eb)*indz(k)*incrxp
      gradc = indxp  * (dens%e *field%e  - dens%c*field%c)
      grade = indxpp * (dens%ee*field%ee - dens%e*field%e)
      gradw = indxm  * (dens%c *field%c  - dens%w*field%w)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=nx-1 .AND. flag(ipjk)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxp*ABS(uc))*deltar*lim
      fn = fn - 0.5D0*(vc+ABS(vc))*s
      ft = ft - 0.5D0*(1-ABS(vc))*(wc+ABS(wc))*s
!
! en-n
      deltar  = dens%en*field%en - dens%n*field%n 
      uc = u%n
      vc = (dx(i+1)*v%n+dx(i)*v%en+dx(i+1)*v%c +dx(i)*v%e  )*indy(j+1)*incrxp
      wc = (dx(i+1)*w%n+dx(i)*w%en+dx(i+1)*w%nb+dx(i)*w%enb)*indz(k)  *incrxp
      gradc = indxp  * (dens%en *field%en  - dens%n *field%n )
      grade = indxpp * (dens%een*field%een - dens%en*field%en)
      gradw = indxm  * (dens%n  *field%n   - dens%wn*field%wn)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=nx-1 .AND. flag(ipjpk)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxp*ABS(uc))*deltar*lim
      fn = fn - 0.5D0*(vc-ABS(vc))*s
      ft = ft - 0.25D0*ABS(vc-ABS(vc))*(wc+ABS(wc))*s
!
! n-wn
      deltar  = dens%n*field%n - dens%wn*field%wn
      uc = u%wn
      vc = (dx(i-1)*v%n+dx(i)*v%wn+dx(i-1)*v%c +dx(i)*v%w  )*indy(j+1)*incrxm
      wc = (dx(i-1)*w%n+dx(i)*w%wn+dx(i-1)*w%nb+dx(i)*w%wnb)*indz(k)  *incrxm
      gradc = indxm  * (dens%n *field%n  - dens%wn *field%wn )
      grade = indxp  * (dens%en*field%en - dens%n  *field%n  )
      gradw = indxmm * (dens%wn*field%wn - dens%wwn*field%wwn)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=2 .AND. flag(imjpk)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxm*ABS(uc))*deltar*lim
      fn = fn + 0.5D0*(vc-ABS(vc))*s
      ft = ft + 0.25D0*ABS(vc-ABS(vc))*(wc+ABS(wc))*s
!
! es-s
      deltar  = dens%es*field%es - dens%s*field%s
      uc = u%s
      vc = (dx(i+1)*v%s+dx(i)*v%es+dx(i+1)*v%ss+dx(i)*v%ess)*indy(j-1)*incrxp
      wc = (dx(i+1)*w%s+dx(i)*w%es+dx(i+1)*w%sb+dx(i)*w%esb)*indz(k)  *incrxp
      gradc = indxp  * (dens%es *field%es  - dens%s *field%s )
      grade = indxpp * (dens%ees*field%ees - dens%es*field%es)
      gradw = indxm  * (dens%s  *field%s   - dens%ws*field%ws)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=nx-1 .AND. flag(ipjmk)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxp*ABS(uc))*deltar*lim
      ft = ft - 0.25D0*ABS(vc+ABS(vc))*(wc+ABS(wc))*s
!
! s-ws
      deltar  = dens%s*field%s - dens%ws*field%ws
      uc = u%ws
      vc = (dx(i-1)*v%s+dx(i)*v%ws+dx(i-1)*v%ss+dx(i)*v%wss)*indy(j-1)*incrxm
      wc = (dx(i-1)*w%s+dx(i)*w%ws+dx(i-1)*w%sb+dx(i)*w%wsb)*indz(k)  *incrxm
      gradc = indxm  * (dens%s *field%s  - dens%ws *field%ws )
      grade = indxp  * (dens%es*field%es - dens%s  *field%s  )
      gradw = indxmm * (dens%ws*field%ws - dens%wws*field%wws)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=2 .AND. flag(imjmk)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxm*ABS(uc))*deltar*lim
      ft = ft + 0.25D0*ABS(vc+ABS(vc))*(wc+ABS(wc))*s
!
! et-t
      deltar  = dens%et*field%et - dens%t*field%t
      uc = u%t
      vc = (dx(i+1)*v%t+dx(i)*v%et+dx(i+1)*v%st+dx(i)*v%est)*indy(j)  *incrxp
      wc = (dx(i+1)*w%t+dx(i)*w%et+dx(i+1)*w%c +dx(i)*w%e  )*indz(k+1)*incrxp
      gradc = indxp  * (dens%et *field%et  - dens%t *field%t )
      grade = indxpp * (dens%eet*field%eet - dens%et*field%et)
      gradw = indxm  * (dens%t  *field%t   - dens%wt*field%wt)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=nx-1 .AND. flag(ipjkp)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxp*ABS(uc))*deltar*lim
      ft = ft - 0.5D0*(1-ABS(vc))*(wc-ABS(wc))*s
!
! t-wt
      deltar  = dens%t*field%t - dens%wt*field%wt
      uc = u%wt
      vc = (dx(i-1)*v%t+dx(i)*v%wt+dx(i-1)*v%st+dx(i)*v%wst)*indy(j)  *incrxm
      wc = (dx(i-1)*w%t+dx(i)*w%wt+dx(i-1)*w%c +dx(i)*w%w  )*indz(k+1)*incrxm
      gradc = indxm  * (dens%t *field%t  - dens%wt *field%wt )
      grade = indxp  * (dens%et*field%et - dens%t  *field%t  )
      gradw = indxmm * (dens%wt*field%wt - dens%wwt*field%wwt)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=2 .AND. flag(imjkp)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxm*ABS(uc))*deltar*lim
      ft = ft + 0.5D0*(1-ABS(vc))*(wc-ABS(wc))*s
!
! ent-nt
      deltar  = dens%ent*field%ent - dens%nt*field%nt
      uc = u%nt
      vc = (dx(i+1)*v%nt+dx(i)*v%ent+dx(i+1)*v%t+dx(i)*v%et)*indy(j+1)*incrxp
      wc = (dx(i+1)*w%nt+dx(i)*w%ent+dx(i+1)*w%n+dx(i)*w%en)*indz(k+1)*incrxp
      gradc = indxp  * (dens%ent *field%ent  - dens%nt *field%nt )
      grade = indxpp * (dens%eent*field%eent - dens%ent*field%ent)
      gradw = indxm  * (dens%nt  *field%nt   - dens%wnt*field%wnt)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=nx-1 .AND. flag(ipjpkp)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxp*ABS(uc))*deltar*lim
      ft = ft - 0.25D0*ABS(vc-ABS(vc))*(wc-ABS(wc))*s
!
! nt-wnt
      deltar  = dens%nt*field%nt - dens%wnt*field%wnt
      uc = u%wnt
      vc = (dx(i-1)*v%nt+dx(i)*v%wnt+dx(i-1)*v%t+dx(i)*v%wt)*indy(j+1)*incrxm
      wc = (dx(i-1)*w%nt+dx(i)*w%wnt+dx(i-1)*w%n+dx(i)*w%wn)*indz(k+1)*incrxm
      gradc = indxm  * (dens%nt*field%nt   - dens%wnt *field%wnt )
      grade = indxp  * (dens%ent*field%ent - dens%nt  *field%nt  )
      gradw = indxmm * (dens%wnt*field%wnt - dens%wwnt*field%wwnt)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=2 .AND. flag(imjpkp)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxm*ABS(uc))*deltar*lim
      ft = ft + 0.25D0*ABS(vc-ABS(vc))*(wc-ABS(wc))*s
!
! est-st
      deltar  = dens%est*field%est - dens%st*field%st
      uc = u%st
      vc = (dx(i+1)*v%st+dx(i)*v%est+dx(i+1)*v%sst+dx(i)*v%esst)*indy(j-1)*incrxp
      wc = (dx(i+1)*w%st+dx(i)*w%est+dx(i+1)*w%s  +dx(i)*w%es  )*indz(k+1)*incrxp
      gradc = indxp  * (dens%est *field%est  - dens%st *field%st )
      grade = indxpp * (dens%eest*field%eest - dens%est*field%est)
      gradw = indxm  * (dens%st  *field%st   - dens%wst*field%wst)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=nx-1 .AND. flag(ipjmkp)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxp*ABS(uc))*deltar*lim
      ft = ft - 0.25D0*ABS(vc+ABS(vc))*(wc-ABS(wc))*s
!
! st-wst
      deltar  = dens%st*field%st - dens%wst*field%wst
      uc = u%wst
      vc = (dx(i-1)*v%st+dx(i)*v%wst+dx(i-1)*v%sst+dx(i)*v%wsst)*indy(j-1)*incrxm
      wc = (dx(i-1)*w%st+dx(i)*w%wst+dx(i-1)*w%s  +dx(i)*w%ws  )*indz(k+1)*incrxm
      gradc = indxm  * (dens%st *field%st  - dens%wst *field%wst )
      grade = indxp  * (dens%est*field%est - dens%st  *field%st  )
      gradw = indxmm * (dens%wst*field%wst - dens%wwst*field%wwst)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=2 .AND. flag(imjmkp)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxm*ABS(uc))*deltar*lim
      ft = ft + 0.25D0*ABS(vc+ABS(vc))*(wc-ABS(wc))*s
!
! c-s
      deltar  = dens%c*field%c - dens%s*field%s
      uc = (dy(j-1)*u%c+dy(j)*u%s+dy(j-1)*u%w+dy(j)*u%ws)*indx(i)*incrym
      vc = v%s
      wc = (dy(j-1)*w%c+dy(j)*w%s+dy(j-1)*w%b+dy(j)*w%sb)*indz(k)*incrym
      gradc = indym  * (dens%c*field%c - dens%s *field%s )
      grads = indymm * (dens%s*field%s - dens%ss*field%ss)
      gradn = indyp  * (dens%n*field%n - dens%c *field%c )
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = grads / gradc
      ELSE
	IF (gradc /= 0.D0) erre = gradn / gradc
      ENDIF
      IF (j/=2 .AND. flag(ijmk)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indym*ABS(vc))*deltar*lim
      ft = ft + 0.5D0*(wc+ABS(wc))*s
      fe = fe + 0.5D0*(1-ABS(wc))*(uc+ABS(uc))*s
!
! n-c
      deltar  = dens%n*field%n - dens%c*field%c
      uc = (dy(j+1)*u%c+dy(j)*u%n+dy(j+1)*u%w+dy(j)*u%wn)*indx(i)*incryp
      vc = v%c
      wc = (dy(j+1)*w%c+dy(j)*w%n+dy(j+1)*w%b+dy(j)*w%nb)*indz(k)*incryp
      gradc = indyp  * (dens%n *field%n  - dens%c*field%c)
      grads = indym  * (dens%c *field%c  - dens%s*field%s)
      gradn = indypp * (dens%nn*field%nn - dens%n*field%n)
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = grads / gradc
      ELSE
	IF (gradc /= 0.D0) erre = gradn / gradc
      ENDIF
      IF (j/=ny-1 .AND. flag(ijpk)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indyp*ABS(vc))*deltar*lim
      ft = ft - 0.5D0*(wc+ABS(wc))*s
      fe = fe - 0.5D0*(1-ABS(wc))*(uc+ABS(uc))*s
!
! en-e
      deltar  = dens%en*field%en - dens%e*field%e
      uc = (dy(j+1)*u%e+dy(j)*u%en+dy(j+1)*u%c +dy(j)*u%n  )*indx(i+1)*incryp
      vc = v%e
      wc = (dy(j+1)*w%e+dy(j)*w%en+dy(j+1)*w%eb+dy(j)*w%enb)*indz(k)  *incryp
      gradc = indyp  * (dens%en *field%en  - dens%e *field%e )
      grads = indym  * (dens%e  *field%e   - dens%es*field%es)
      gradn = indypp * (dens%enn*field%enn - dens%en*field%en)
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = grads / gradc
      ELSE
	IF (gradc /= 0.D0) erre = gradn / gradc
      ENDIF
      IF (j/=ny-1 .AND. flag(ipjpk)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indyp*ABS(vc))*deltar*lim
      fe = fe - 0.5D0*(1-ABS(wc))*(uc-ABS(uc))*s
!
! e-es
      deltar  = dens%e*field%e - dens%es*field%es
      uc = (dy(j-1)*u%e+dy(j)*u%es+dy(j-1)*u%c +dy(j)*u%s  )*indx(i+1)*incrym
      vc = v%es
      wc = (dy(j-1)*w%e+dy(j)*w%es+dy(j-1)*w%eb+dy(j)*w%esb)*indz(k)  *incrym
      gradc = indym  * (dens%e *field%e  - dens%es *field%es )
      grads = indymm * (dens%es*field%es - dens%ess*field%ess)
      gradn = indyp  * (dens%en*field%en - dens%e  *field%e  )
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = grads / gradc
      ELSE
	IF (gradc /= 0.D0) erre = gradn / gradc
      ENDIF
      IF (j/=2 .AND. flag(ipjmk)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indym*ABS(vc))*deltar*lim
      fe = fe + 0.5D0*(1-ABS(wc))*(uc-ABS(uc))*s
!
! nt-t
      deltar  = dens%nt*field%nt - dens%t*field%t
      uc = (dy(j+1)*u%t+dy(j)*u%nt+dy(j+1)*u%wt+dy(j)*u%wnt)*indx(i)  *incryp
      vc = v%t
      wc = (dy(j+1)*w%t+dy(j)*w%nt+dy(j+1)*w%c +dy(j)*w%n  )*indz(k+1)*incryp
      gradc = indyp  * (dens%nt *field%nt  - dens%t *field%t )
      grads = indym  * (dens%t  *field%t   - dens%st*field%st)
      gradn = indypp * (dens%nnt*field%nnt - dens%nt*field%nt)
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = grads / gradc
      ELSE
	IF (gradc /= 0.D0) erre = gradn / gradc
      ENDIF
      IF (j/=ny-1 .AND. flag(ijpkp)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indyp*ABS(vc))*deltar*lim
      ft = ft - 0.5D0*(wc-ABS(wc))*s
      fe = fe - 0.25D0*ABS(wc-ABS(wc))*(uc+ABS(uc))*s
!
! t-st
      deltar  = dens%t*field%t - dens%st*field%st
      uc = (dy(j-1)*u%t+dy(j)*u%st+dy(j-1)*u%wt+dy(j)*u%wst)*indx(i)  *incrym
      vc = v%st
      wc = (dy(j-1)*w%t+dy(j)*w%st+dy(j-1)*w%c +dy(j)*w%s  )*indz(k+1)*incrym
      gradc = indym  * (dens%t *field%t  - dens%st *field%st )
      grads = indymm * (dens%st*field%st - dens%sst*field%sst)
      gradn = indyp  * (dens%nt*field%nt - dens%t  *field%t  )
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = grads / gradc
      ELSE
	IF (gradc /= 0.D0) erre = gradn / gradc
      ENDIF
      IF (j/=2 .AND. flag(ijmkp)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indym*ABS(vc))*deltar*lim
      ft = ft + 0.5D0*(wc-ABS(wc))*s
      fe = fe + 0.25D0*ABS(wc-ABS(wc))*(uc+ABS(uc))*s
!
! nb-b
      deltar  = dens%nb*field%nb - dens%b*field%b
      IF (k <= 2) deltar = 0.D0
      uc = (dy(j+1)*u%b+dy(j)*u%nb+dy(j+1)*u%wb+dy(j)*u%wnb)*indx(i)  *incryp
      vc = v%b
      wc = (dy(j+1)*w%b+dy(j)*w%nb+dy(j+1)*w%bb+dy(j)*w%nbb)*indz(k-1)*incryp
      gradc = indyp  * (dens%nb *field%nb  - dens%b *field%b )
      grads = indym  * (dens%b  *field%b   - dens%sb*field%sb)
      gradn = indypp * (dens%nnb*field%nnb - dens%nb*field%nb)
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = grads / gradc
      ELSE
	IF (gradc /= 0.D0) erre = gradn / gradc
      ENDIF
      IF (j/=ny-1 .AND. flag(ijpkm)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indyp*ABS(vc))*deltar*lim
      fe = fe - 0.25D0*ABS(wc+ABS(wc))*(uc+ABS(uc))*s
!
! b-sb
      deltar  = dens%b*field%b - dens%sb*field%sb
      IF (k <= 2) deltar = 0.D0
      uc = (dy(j-1)*u%b+dy(j)*u%sb+dy(j-1)*u%wb+dy(j)*u%wsb)*indx(i)  *incrym
      vc = v%sb
      wc = (dy(j-1)*w%b+dy(j)*w%sb+dy(j-1)*w%bb+dy(j)*w%sbb)*indz(k-1)*incrym
      gradc = indym  * (dens%b *field%b  - dens%sb *field%sb )
      grads = indymm * (dens%sb*field%sb - dens%ssb*field%ssb)
      gradn = indyp  * (dens%nb*field%nb - dens%b  *field%b  )
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = grads / gradc
      ELSE
	IF (gradc /= 0.D0) erre = gradn / gradc
      ENDIF
      IF (j/=2 .AND. flag(ijmkm)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indym*ABS(vc))*deltar*lim
      fe = fe + 0.25D0*ABS(wc+ABS(wc))*(uc+ABS(uc))*s
!
! ent-et
      deltar  = dens%ent*field%ent - dens%et*field%et
      uc = (dy(j+1)*u%et+dy(j)*u%ent+dy(j+1)*u%t+dy(j)*u%nt)*indx(i+1)*incryp
      vc = v%et
      wc = (dy(j+1)*w%et+dy(j)*w%ent+dy(j+1)*w%e+dy(j)*w%en)*indz(k+1)*incryp
      gradc = indyp  * (dens%ent *field%ent  - dens%et *field%et )
      grads = indym  * (dens%et  *field%et   - dens%est*field%est)
      gradn = indypp * (dens%ennt*field%ennt - dens%ent*field%ent)
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = grads / gradc
      ELSE
	IF (gradc /= 0.D0) erre = gradn / gradc
      ENDIF
      IF (j/=ny-1 .AND. flag(ipjpkp)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indyp*ABS(vc))*deltar*lim
      fe = fe - 0.25D0*ABS(wc-ABS(wc))*(uc-ABS(uc))*s
!
! et-est
      deltar  = dens%et*field%et - dens%est*field%est
      uc = (dy(j-1)*u%et+dy(j)*u%est+dy(j-1)*u%t+dy(j)*u%st)*indx(i+1)*incrym
      vc = v%est
      wc = (dy(j-1)*w%et+dy(j)*w%est+dy(j-1)*w%e+dy(j)*w%es)*indz(k+1)*incrym
      gradc = indym  * (dens%et *field%et  - dens%est *field%est )
      grads = indymm * (dens%est*field%est - dens%esst*field%esst)
      gradn = indyp  * (dens%ent*field%ent - dens%et  *field%et  )
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = grads / gradc
      ELSE
	IF (gradc /= 0.D0) erre = gradn / gradc
      ENDIF
      IF (j/=2 .AND. flag(ipjmkp)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indym*ABS(vc))*deltar*lim
      fe = fe + 0.25D0*ABS(wc-ABS(wc))*(uc-ABS(uc))*s
!
! enb-eb
      deltar  = dens%enb*field%enb - dens%eb*field%eb
      IF (k <= 2) deltar = 0.D0
      uc = (dy(j+1)*u%eb+dy(j)*u%enb+dy(j+1)*u%b  +dy(j)*u%nb  )*indx(i+1)*incryp
      vc = v%eb
      wc = (dy(j+1)*w%eb+dy(j)*w%enb+dy(j+1)*w%ebb+dy(j)*w%enbb)*indz(k-1)*incryp
      gradc = indyp  * (dens%enb *field%enb  - dens%eb *field%eb )
      grads = indym  * (dens%eb  *field%eb   - dens%esb*field%esb)
      gradn = indypp * (dens%ennb*field%ennb - dens%enb*field%enb)
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = grads / gradc
      ELSE
	IF (gradc /= 0.D0) erre = gradn / gradc
      ENDIF
      IF (j/=ny-1 .AND. flag(ipjpkm)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indyp*ABS(vc))*deltar*lim
      fe = fe - 0.25D0*ABS(wc+ABS(wc))*(uc-ABS(uc))*s
!
! eb-esb
      deltar  = dens%eb*field%eb - dens%esb*field%esb
      IF (k <= 2) deltar = 0.D0
      uc = (dy(j-1)*u%eb+dy(j)*u%esb+dy(j-1)*u%b  +dy(j)*u%sb  )*indx(i+1)*incrym
      vc = v%esb
      wc = (dy(j-1)*w%eb+dy(j)*w%esb+dy(j-1)*w%ebb+dy(j)*w%esbb)*indz(k-1)*incrym
      gradc = indym  * (dens%eb *field%eb  - dens%esb *field%esb )
      grads = indymm * (dens%esb*field%esb - dens%essb*field%essb)
      gradn = indyp  * (dens%enb*field%enb - dens%eb  *field%eb  )
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = grads / gradc
      ELSE
	IF (gradc /= 0.D0) erre = gradn / gradc
      ENDIF
      IF (j/=2 .AND. flag(ipjmkm)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indym*ABS(vc))*deltar*lim
      fe = fe + 0.25D0*ABS(wc+ABS(wc))*(uc-ABS(uc))*s
!
! c-b
      deltar  = dens%c*field%c - dens%b*field%b
      wc = w%b
      uc = (dz(k-1)*u%c+dz(k)*u%b+dz(k-1)*u%w+dz(k)*u%wb)*indx(i)*incrzm
      vc = (dz(k-1)*v%c+dz(k)*v%b+dz(k-1)*v%s+dz(k)*v%sb)*indy(j)*incrzm
      gradc = indzm  * (dens%c*field%c - dens%b *field%b )
      gradt = indzp  * (dens%t*field%t - dens%c *field%c )
      gradb = indzmm * (dens%b*field%b - dens%bb*field%bb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=2 .AND. flag(ijkm)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indzm*ABS(wc))*deltar*lim
      fe = fe + 0.5D0*(uc+ABS(uc))*s
      fn = fn + 0.5D0*(1-ABS(uc))*(vc+ABS(vc))*s
!
! t-c
      deltar  = dens%t*field%t - dens%c*field%c
      uc = (dz(k+1)*u%c+dz(k)*u%t+dz(k+1)*u%w+dz(k)*u%wt)*indx(i)*incrzp
      vc = (dz(k+1)*v%c+dz(k)*v%t+dz(k+1)*v%s+dz(k)*v%st)*indy(j)*incrzp
      wc = w%c
      gradc = indzp  * (dens%t *field%t  - dens%c*field%c)
      gradt = indzpp * (dens%tt*field%tt - dens%t*field%t)
      gradb = indzm  * (dens%c *field%c  - dens%b*field%b)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(ijkp)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indzp*ABS(wc))*deltar*lim
      fe = fe - 0.5D0*(uc+ABS(uc))*s
      fn = fn - 0.5D0*(1-ABS(uc))*(vc+ABS(vc))*s
!
! wt-w
      deltar  = dens%wt*field%wt - dens%w*field%w
      wc = w%w
      uc = (dz(k+1)*u%w+dz(k)*u%wt+dz(k+1)*u%ww+dz(k)*u%wwt)*indx(i-1)*incrzp
      vc = (dz(k+1)*v%w+dz(k)*v%wt+dz(k+1)*v%ws+dz(k)*v%wst)*indy(j)  *incrzp
      gradc = indzp  * (dens%wt *field%wt  - dens%w *field%w )
      gradt = indzpp * (dens%wtt*field%wtt - dens%wt*field%wt)
      gradb = indzm  * (dens%w  *field%w   - dens%wb*field%wb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(imjkp)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indzp*ABS(wc))*deltar*lim
      fn = fn - 0.25D0*ABS(uc+ABS(uc))*(vc+ABS(vc))*s
!
! w-wb
      deltar  = dens%w*field%w - dens%wb*field%wb
      wc = w%wb
      uc = (dz(k-1)*u%w+dz(k)*u%wb+dz(k-1)*u%ww+dz(k)*u%wwb)*indx(i-1)*incrzm
      vc = (dz(k-1)*v%w+dz(k)*v%wb+dz(k-1)*v%ws+dz(k)*v%wsb)*indy(j)  *incrzm
      gradc = indzm  * (dens%w *field%w  - dens%wb *field%wb )
      gradt = indzp  * (dens%wt*field%wt - dens%w  *field%w  )
      gradb = indzmm * (dens%wb*field%wb - dens%wbb*field%wbb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=2 .AND. flag(imjkm)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indzm*ABS(wc))*deltar*lim
      fn = fn + 0.25D0*ABS(uc+ABS(uc))*(vc+ABS(vc))*s
!
! et-e
      deltar  = dens%et*field%et - dens%e*field%e
      wc = w%e
      uc = (dz(k+1)*u%e+dz(k)*u%et+dz(k+1)*u%c +dz(k)*u%t  )*indx(i+1)*incrzp
      vc = (dz(k+1)*v%e+dz(k)*v%et+dz(k+1)*v%es+dz(k)*v%est)*indy(j)  *incrzp
      gradc = indzp  * (dens%et *field%et  - dens%e *field%e )
      gradt = indzpp * (dens%ett*field%ett - dens%et*field%et)
      gradb = indzm  * (dens%e  *field%e   - dens%eb*field%eb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(ipjkp)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indzp*ABS(wc))*deltar*lim
      fe = fe - 0.5D0*(uc-ABS(uc))*s
      fn = fn - 0.25D0*ABS(uc-ABS(uc))*(vc+ABS(vc))*s
!
! e-eb
      deltar  = dens%e*field%e - dens%eb*field%eb
      wc = w%eb
      uc = (dz(k-1)*u%e+dz(k)*u%eb+dz(k-1)*u%c +dz(k)*u%b  )*indx(i+1)*incrzm
      vc = (dz(k-1)*v%e+dz(k)*v%eb+dz(k-1)*v%es+dz(k)*v%esb)*indy(j)  *incrzm
      gradc = indzm  * (dens%e *field%e  - dens%eb *field%eb )
      gradt = indzp  * (dens%et*field%et - dens%e  *field%e  )
      gradb = indzmm * (dens%eb*field%eb - dens%ebb*field%ebb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=2 .AND. flag(ipjkm)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indzm*ABS(wc))*deltar*lim
      fe = fe + 0.5D0*(uc-ABS(uc))*s
      fn = fn + 0.25D0*ABS(uc-ABS(uc))*(vc+ABS(vc))*s
!
! nt-n
      deltar  = dens%nt*field%nt - dens%n*field%n
      wc = w%n
      uc = (dz(k+1)*u%n+dz(k)*u%nt+dz(k+1)*u%wn+dz(k)*u%wnt)*indx(i)  *incrzp
      vc = (dz(k+1)*v%n+dz(k)*v%nt+dz(k+1)*v%c +dz(k)*v%t  )*indy(j+1)*incrzp
      gradc = indzp  * (dens%nt *field%nt  - dens%n *field%n )
      gradt = indzpp * (dens%ntt*field%ntt - dens%nt*field%nt)
      gradb = indzm  * (dens%n  *field%n   - dens%nb*field%nb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(ijpkp)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indzp*ABS(wc))*deltar*lim
      fn = fn - 0.5D0*(1-ABS(uc))*(vc-ABS(vc))*s
!
! n-nb
      deltar  = dens%n*field%n - dens%nb*field%nb
      wc = w%nb
      uc = (dz(k-1)*u%n+dz(k)*u%nb+dz(k-1)*u%wn+dz(k)*u%wnb)*indx(i)  *incrzm
      vc = (dz(k-1)*v%n+dz(k)*v%nb+dz(k-1)*v%c +dz(k)*v%b  )*indy(j+1)*incrzm
      gradc = indzm  * (dens%n *field%n  - dens%nb *field%nb )
      gradt = indzp  * (dens%nt*field%nt - dens%n  *field%n  )
      gradb = indzmm * (dens%nb*field%nb - dens%nbb*field%nbb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=2 .AND. flag(ijpkm)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indzm*ABS(wc))*deltar*lim
      fn = fn + 0.5D0*(1-ABS(uc))*(vc-ABS(vc))*s
!
! wnt-wn
      deltar  = dens%wnt*field%wnt - dens%wn*field%wn
      wc = w%wn
      uc = (dz(k+1)*u%wn+dz(k)*u%wnt+dz(k+1)*u%wwn+dz(k)*u%wwnt)*indx(i-1)*incrzp
      vc = (dz(k+1)*v%wn+dz(k)*v%wnt+dz(k+1)*v%w  +dz(k)*v%wt  )*indy(j+1)*incrzp
      gradc = indzp  * (dens%wnt *field%wnt  - dens%wn *field%wn )
      gradt = indzpp * (dens%wntt*field%wntt - dens%wnt*field%wnt)
      gradb = indzm  * (dens%wn  *field%wn   - dens%wnb*field%wnb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(imjpkp)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indzp*ABS(wc))*deltar*lim
      fn = fn - 0.25D0*ABS(uc+ABS(uc))*(vc-ABS(vc))*s
!
! wn-wnb
      deltar  = dens%wn*field%wn - dens%wnb*field%wnb
      wc = w%wnb
      uc = (dz(k-1)*u%wn+dz(k)*u%wnb+dz(k-1)*u%wwn+dz(k)*u%wwnb)*indx(i-1)*incrzm
      vc = (dz(k-1)*v%wn+dz(k)*v%wnb+dz(k-1)*v%w  +dz(k)*v%wb  )*indy(j+1)*incrzm
      gradc = indzm  * (dens%wn *field%wn  - dens%wnb *field%wnb )
      gradt = indzp  * (dens%wnt*field%wnt - dens%wn  *field%wn  )
      gradb = indzmm * (dens%wnb*field%wnb - dens%wnbb*field%wnbb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=2 .AND. flag(imjpkm)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indzm*ABS(wc))*deltar*lim
      fn = fn + 0.25D0*ABS(uc+ABS(uc))*(vc-ABS(vc))*s
!
! ent-en
      deltar  = dens%ent*field%ent - dens%en*field%en
      wc = w%en
      uc = (dz(k+1)*u%en+dz(k)*u%ent+dz(k+1)*u%n+dz(k)*u%nt)*indx(i+1)*incrzp
      vc = (dz(k+1)*v%en+dz(k)*v%ent+dz(k+1)*v%e+dz(k)*v%et)*indy(j+1)*incrzp
      gradc = indzp  * (dens%ent *field%ent  - dens%en *field%en )
      gradt = indzpp * (dens%entt*field%entt - dens%ent*field%ent)
      gradb = indzm  * (dens%en  *field%en   - dens%enb*field%enb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(ipjpkp)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indzp*ABS(wc))*deltar*lim
      fn = fn - 0.25D0*ABS(uc-ABS(uc))*(vc-ABS(vc))*s
!
! en-enb
      deltar  = dens%en*field%en - dens%enb*field%enb
      wc = w%enb
      uc = (dz(k-1)*u%en+dz(k)*u%enb+dz(k-1)*u%n+dz(k)*u%nb)*indx(i+1)*incrzm
      vc = (dz(k-1)*v%en+dz(k)*v%enb+dz(k-1)*v%e+dz(k)*v%eb)*indy(j+1)*incrzm
      gradc = indzm  * (dens%en *field%en  - dens%enb *field%enb )
      gradt = indzp  * (dens%ent*field%ent - dens%en  *field%en  )
      gradb = indzmm * (dens%enb*field%enb - dens%enbb*field%enbb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=2 .AND. flag(ipjpkm)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indzm*ABS(wc))*deltar*lim
      fn = fn + 0.25D0*ABS(uc-ABS(uc))*(vc-ABS(vc))*s
!
      RETURN
      END SUBROUTINE ctu3_fsc_3d
!-----------------------------------------------------------------------
      END MODULE convective_fluxes_sc
!-----------------------------------------------------------------------
