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
        MODULE PROCEDURE ctu1_fsc_2d
      END INTERFACE
      INTERFACE ctu2_fsc
        MODULE PROCEDURE ctu2_fsc_2d
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
      END MODULE convective_fluxes_sc
!-----------------------------------------------------------------------
