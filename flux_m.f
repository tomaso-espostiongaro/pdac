!----------------------------------------------------------------------
      MODULE convective_mass_fluxes
!
! ... This module computes convective fluxes of _______ density _______
! ... by using an upwind scheme.
! ... Standard compass notation is adopted to locate cell neighbours:
! ... %e=East, %w=West, %n=North, %s=South, %t=Top, %b=Bottom, etc.
! ... The computational stencil is defined in the `set_indexes' module
!
!----------------------------------------------------------------------
!
      USE flux_limiters, ONLY: limiters, lm

      IMPLICIT NONE
      SAVE
!
      REAL*8, PRIVATE :: cs                      ! convective stream   !
      REAL*8, PRIVATE :: cn                      ! Courant number      !
      REAL*8, PRIVATE :: fou                     ! first order upwind  !
      REAL*8, PRIVATE :: upwnd                   ! upwinded variable   !
      REAL*8, PRIVATE :: centrd                  ! centered variable   !
      REAL*8, PRIVATE :: lim                     ! limiter             !
      REAL*8, PRIVATE :: erre                    ! gradients ratio     !
      REAL*8, PRIVATE :: incr                    ! increment           !
      REAL*8, PUBLIC  :: upc_e, upc_n, upc_t     ! upwinded/centered   !
      REAL*8, PUBLIC  :: upc_w, upc_s, upc_b     ! fields              !

!
      INTERFACE fmas
        MODULE PROCEDURE fmas_2d, fmas_3d
      END INTERFACE
      INTERFACE masf
        MODULE PROCEDURE masf_2d, masf_3d
      END INTERFACE
      INTERFACE ctu1_fmas
        MODULE PROCEDURE ctu1_fmas_2d
      END INTERFACE
      INTERFACE ctu2_fmas
        MODULE PROCEDURE ctu2_fmas_2d
      END INTERFACE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE masf_3d( fe, fn, ft, fw, fs, fb, dens, u, v, w, ijk )
!
! ... This routine computes convective mass fluxes by using
! ... Donor Cell technique for first order upwind
!
      USE set_indexes, ONLY: stencil
!
      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: fe, fn, ft, fw, fs, fb
      TYPE(stencil), INTENT(IN) :: u, v, w, dens
      INTEGER, INTENT(IN) :: ijk
!
! ... Correction factors used in the computation of
! ... the derivative of the mass residual with pressure
! ... (routine 'betas')
!
      upc_e = 1.D0
      upc_w = 1.D0
      upc_t = 1.D0
      upc_b = 1.D0
      upc_n = 1.D0
      upc_s = 1.D0
!
! ... West, South and Bottom fluxes
!
      cs = u%w
      fw = 0.5D0*(cs-ABS(cs))*dens%c + 0.5D0*(cs+ABS(cs))*dens%w

      cs = v%s
      fs = 0.5D0*(cs-ABS(cs))*dens%c + 0.5D0*(cs+ABS(cs))*dens%s

      cs = w%b
      fb = 0.5D0*(cs-ABS(cs))*dens%c + 0.5D0*(cs+ABS(cs))*dens%b
!
! ... East, North and Top fluxes
!
      cs = u%c
      fe = 0.5D0*(cs-ABS(cs))*dens%e + 0.5D0*(cs+ABS(cs))*dens%c

      cs = v%c
      fn = 0.5D0*(cs-ABS(cs))*dens%n + 0.5D0*(cs+ABS(cs))*dens%c

      cs = w%c
      ft = 0.5D0*(cs-ABS(cs))*dens%t + 0.5D0*(cs+ABS(cs))*dens%c
!
      RETURN
      END SUBROUTINE masf_3d
!----------------------------------------------------------------------
      SUBROUTINE fmas_3d(fe, fn, ft, fw, fs, fb, dens, u, v, w, ijk)
!
      USE dimensions
      USE domain_mapping, ONLY: myijk, meshinds
      USE grid, ONLY: dx, dy, dz, flag, immb_cell, fluid
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: stencil, imjk, ijmk, ijkm, ipjk, ijpk, ijkp
      USE time_parameters, ONLY: dt
      IMPLICIT NONE
!
      REAL*8, INTENT(INOUT) :: fe, fn, ft, fw, fs, fb
      TYPE(stencil), INTENT(IN) :: dens, u, v, w
      INTEGER, INTENT(IN) :: ijk

      INTEGER :: im2, ip2, jm2, jp2, km2, kp2
      INTEGER :: i,j,k,imesh
!
      INTEGER, SAVE :: ijk_old = -1
      REAL*8, SAVE :: dxmm, dxm, dxp, dxpp, indxpp, indxp, indxm, indxmm
      REAL*8, SAVE :: dymm, dym, dyp, dypp, indypp, indyp, indym, indymm
      REAL*8, SAVE :: dzmm, dzm, dzp, dzpp, indzpp, indzp, indzm, indzmm
      REAL*8, SAVE :: gradc, grade, gradw, gradn, grads, gradt, gradb
!
      CALL meshinds(ijk,imesh,i,j,k)

      IF( ijk /= ijk_old ) THEN
!
         im2 = MAX(  1, i-2 )
         ip2 = MIN( nx, i+2 )
         jm2 = MAX(  1, j-2 )
         jp2 = MIN( ny, j+2 )
         km2 = MAX(  1, k-2 )
         kp2 = MIN( nz, k+2 )

         dxmm = dx(i-1)+ dx(im2)
         dxm  = dx(i)  + dx(i-1)
         dxp  = dx(i)  + dx(i+1)
         dxpp = dx(i+1)+ dx(ip2)
         dymm = dy(j-1)+ dy(jm2)
         dym  = dy(j)  + dy(j-1)
         dyp  = dy(j)  + dy(j+1)
         dypp = dy(j+1)+ dy(jp2)
         dzmm = dz(k-1)+ dz(km2)
         dzm  = dz(k)  + dz(k-1)
         dzp  = dz(k)  + dz(k+1)
         dzpp = dz(k+1)+ dz(kp2)

         indxmm = 1.D0/dxmm
         indxm  = 1.D0/dxm
         indxp  = 1.D0/dxp
         indxpp = 1.D0/dxpp
         indymm = 1.D0/dymm
         indym  = 1.D0/dym
         indyp  = 1.D0/dyp
         indypp = 1.D0/dypp
         indzmm = 1.D0/dzmm
         indzm  = 1.D0/dzm
         indzp  = 1.D0/dzp
         indzpp = 1.D0/dzpp
     
         ijk_old = ijk
      END IF
!
! ... MUSCL reconstruction of fields
!
! ... on West volume boundary
!
      gradc = 2.D0 * indxm  * (dens%c - dens%w )
      gradw = 2.D0 * indxmm * (dens%w - dens%ww)
      grade = 2.D0 * indxp  * (dens%e - dens%c )
!
      lim = 0.D0
      erre = 0.D0
!
      cs = u%w
      IF (cs >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw / gradc
        incr = 0.5D0 * dx(i-1)
      ELSE IF (cs < 0.D0) THEN
        IF (gradc /= 0.D0) erre = grade / gradc
        incr = - 0.5D0 * dx(i)
      ENDIF
!
      IF (i/=2 .AND. flag(imjk)==fluid) CALL limiters(lm,lim,erre)
!
      upwnd = lim * gradc * incr
      centrd = ( dx(i) * dens%w + dx(i-1) * dens%c ) * indxm

      fw = fw + upwnd * cs

      upc_w = fw / cs / centrd
!
! ... on East volume boundary
!
      gradw = gradc
      gradc = grade
      grade = 2.D0 * indxpp * (dens%ee - dens%e)
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
      IF (i/=nx-1 .AND. flag(ipjk)==fluid) CALL limiters(lm,lim,erre)
!
      upwnd = lim * gradc * incr
      centrd = ( dx(i) * dens%e + dx(i+1) * dens%c ) * indxp

      fe = fe + upwnd * cs

      upc_e = fe / cs / centrd
!
! ... on South volume boundary
!
      gradc = 2.D0 * indym  * (dens%c - dens%s )
      grads = 2.D0 * indymm * (dens%s - dens%ss)
      gradn = 2.D0 * indyp  * (dens%n - dens%c )
!
      lim = 0.D0
      erre = 0.D0
!
      cs = v%s
      IF (cs >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = grads / gradc
	incr = 0.5D0 * dy(j-1)
      ELSE IF (cs < 0.D0) THEN
	IF (gradc /= 0.D0) erre = gradn / gradc
	incr = - 0.5D0 * dy(j)
      ENDIF
!
      IF (j/=2 .AND. flag(ijmk)==fluid) CALL limiters(lm,lim,erre)
!
      upwnd = lim * gradc * incr
      centrd = ( dy(j) * dens%s + dy(j-1) * dens%c ) * indym

      fs = fs + upwnd * cs

      upc_s = fs / cs / centrd
!
! ... on North volume boundary
!
      grads = gradc
      gradc = gradn
      gradn = 2.D0 * indypp * (dens%nn - dens%n)
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
      IF (j/=ny-1 .AND. flag(ijpk)==fluid) CALL limiters(lm,lim,erre)
!
      upwnd = lim * gradc * incr
      centrd = ( dy(j) * dens%n + dy(j+1) * dens%c ) * indyp

      fn = fn + upwnd * cs

      upc_n = fn / cs / centrd
!
! ... on Bottom volume boundary
!
      gradc = 2.D0 * indzm  * (dens%c - dens%b )
      gradb = 2.D0 * indzmm * (dens%b - dens%bb)
      gradt = 2.D0 * indzp  * (dens%t - dens%c )
!
      lim = 0.D0
      erre = 0.D0
!
      cs = w%b
      IF (cs >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = gradb / gradc
	incr = 0.5D0 * dz(k-1)
      ELSE IF (cs < 0.D0) THEN
	IF (gradc /= 0.D0) erre = gradt / gradc
	incr = - 0.5D0 * dz(k)
      ENDIF
!
      IF (k/=2 .AND. flag(ijkm)==fluid) CALL limiters(lm,lim,erre)
!
      upwnd = lim * gradc * incr
      centrd = ( dz(k) * dens%b + dz(k-1) * dens%c ) * indzm

      fb = fb + upwnd * cs

      upc_b = fb / cs / centrd
!
! ... on Top volume boundary
!
      gradb = gradc
      gradc = gradt
      gradt = 2.D0 * indzpp * (dens%tt - dens%t)
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
      IF (k/=nz-1 .AND. flag(ijkp)==fluid) CALL limiters(lm,lim,erre)
!
      upwnd = lim * gradc * incr
      centrd = ( dz(k) * dens%t + dz(k+1) * dens%c ) * indzp

      ft = ft + upwnd * cs

      upc_t = ft / cs / centrd
!
      RETURN
      END SUBROUTINE fmas_3d
!-----------------------------------------------------------------------
      SUBROUTINE masf_2d(fe, ft, fw, fb, dens, u, w, ij)
!
! ... This routine computes convective mass fluxes by using
! ... Donor Cell technique for first order upwind
!
      USE dimensions
      USE domain_mapping, ONLY: myijk
      USE grid, ONLY: rb
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: stencil
!
      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: fe, ft, fw, fb
      TYPE(stencil), INTENT(IN) :: u, w, dens
      INTEGER, INTENT(IN) :: ij
      INTEGER :: i,j,imesh
!
      imesh = myijk( ip0_jp0_kp0_, ij)
      j = ( imesh - 1 ) / nx + 1
      i = MOD( ( imesh - 1 ), nx) + 1
!
      upc_e = 1.D0
      upc_w = 1.D0
      upc_t = 1.D0
      upc_b = 1.D0
!
! ... West, Bottom fluxes
!
      cs = u%w
      fw = 0.5D0*(cs-ABS(cs))*dens%c*rb(i-1) + 0.5D0*(cs+ABS(cs))*dens%w*rb(i-1)

      cs = w%b
      fb = 0.5D0*(cs-ABS(cs))*dens%c + 0.5D0*(cs+ABS(cs))*dens%b
!
! ... East, Top fluxes
!
      cs = u%c
      fe = 0.5D0*(cs-ABS(cs))*dens%e*rb(i) + 0.5D0*(cs+ABS(cs))*dens%c*rb(i)

      cs = w%c
      ft = 0.5D0*(cs-ABS(cs))*dens%t + 0.5D0*(cs+ABS(cs))*dens%c
!
      RETURN
      END SUBROUTINE masf_2d
!----------------------------------------------------------------------
      SUBROUTINE fmas_2d(fe, ft, fw, fb, dens, u, w, ij)
!
      USE dimensions
      USE domain_mapping, ONLY: myijk
      USE grid, ONLY: dx, rb, dz, fluid, flag
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: stencil, imjk, ijkm, ipjk, ijkp
      USE time_parameters, ONLY: dt
      IMPLICIT NONE
!
      REAL*8, INTENT(INOUT) :: fe, ft, fw, fb
      TYPE(stencil), INTENT(IN) :: dens, u, w
      INTEGER, INTENT(IN) :: ij
!
      INTEGER :: i, k, imesh
      REAL*8 :: dxmm, dxm, dxp, dxpp, indxpp, indxp, indxm, indxmm
      REAL*8 :: dzmm, dzm, dzp, dzpp, indzpp, indzp, indzm, indzmm
      REAL*8 :: gradc, grade, gradw, gradt, gradb

      INTEGER :: im2, ip2, km2, kp2
!
      imesh = myijk( ip0_jp0_kp0_, ij)
      k = ( imesh - 1 ) / nx + 1
      i = MOD( ( imesh - 1 ), nx) + 1

      im2 = MAX(  1, i-2 )
      ip2 = MIN( nx, i+2 )
      km2 = MAX(  1, k-2 )
      kp2 = MIN( nz, k+2 )
!
      dxmm = dx(i-1)+ dx(im2)
      dxm  = dx(i)  + dx(i-1)
      dxp  = dx(i)  + dx(i+1)
      dxpp = dx(i+1)+ dx(ip2)
      dzmm = dz(k-1)+ dz(km2)
      dzm  = dz(k)  + dz(k-1)
      dzp  = dz(k)  + dz(k+1)
      dzpp = dz(k+1)+ dz(kp2)

      indxmm = 1.D0/dxmm
      indxm  = 1.D0/dxm
      indxp  = 1.D0/dxp
      indxpp = 1.D0/dxpp
      indzmm = 1.D0/dzmm
      indzm  = 1.D0/dzm
      indzp  = 1.D0/dzp
      indzpp = 1.D0/dzpp
!
! ... MUSCL reconstruction of fields
!
! ... on West volume boundary
!
      gradc = 2.D0 * indxm  * (dens%c - dens%w )
      gradw = 2.D0 * indxmm * (dens%w - dens%ww)
      grade = 2.D0 * indxp  * (dens%e - dens%c )
!
      lim = 0.D0
      erre = 0.D0
!
      cs = u%w
      IF (cs >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = gradw / gradc
	incr = 0.5D0 * dx(i-1)
      ELSE IF (cs < 0.D0) THEN
	IF (gradc /= 0.D0) erre = grade / gradc
	incr = - 0.5D0 * dx(i)
      ENDIF
!
      IF (i/=2 .AND. flag(imjk)==fluid) CALL limiters(lm,lim,erre)
!
      upwnd = lim * gradc * incr
      centrd = (dx(i)*dens%w+dx(i-1)*dens%c)*indxm

      fw = fw + upwnd * cs * rb(i-1)

      upc_w = fw / (cs * rb(i-1)) / centrd
!
! ... on East volume boundary
!
      gradw = gradc
      gradc = grade
      grade = 2.D0 * indxpp * (dens%ee - dens%e)
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
      IF (i/=nx-1 .AND. flag(ipjk)==fluid) CALL limiters(lm,lim,erre)
!
      upwnd = lim * gradc * incr
      centrd = (dx(i)*dens%e+dx(i+1)*dens%c)*indxp
      
      fe = fe + upwnd * cs * rb(i)

      upc_e = fe / (cs * rb(i)) / centrd
!
! ... on Bottom volume boundary
!
      gradc = 2.D0 * indzm  * (dens%c - dens%b )
      gradb = 2.D0 * indzmm * (dens%b - dens%bb)
      gradt = 2.D0 * indzp  * (dens%t - dens%c )
!
      lim = 0.D0
      erre = 0.D0
!
      cs = w%b
      IF (cs >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = gradb / gradc
	incr = 0.5D0 * dz(k-1)
      ELSE IF (cs < 0.D0) THEN
	IF (gradc /= 0.D0) erre = gradt / gradc
	incr = - 0.5D0 * dz(k)
      ENDIF
!
      IF (k/=2 .AND. flag(ijkm)==fluid) CALL limiters(lm,lim,erre)
!
      upwnd = lim * gradc * incr
      centrd = (dz(k)*dens%b+dz(k-1)*dens%c)*indzm

      fb = fb + upwnd * cs

      upc_b = fb / cs  / centrd
!
! ... on Top volume boundary
!
      gradb = gradc
      gradc = gradt
      gradt = 2.D0 * indzpp * (dens%tt - dens%t)
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
      IF (k/=nz-1 .AND. flag(ijkp)==fluid) CALL limiters(lm,lim,erre)
!
      upwnd = lim * gradc * incr
      centrd = (dz(k)*dens%t+dz(k+1)*dens%c)*indzp

      ft = ft + upwnd * cs

      upc_t = ft / cs / centrd
!
      RETURN
      END SUBROUTINE fmas_2d
!-----------------------------------------------------------------------
      SUBROUTINE ctu1_fmas_2d(fe, ft, fw, fb, dens, u, w, ij)
!
! ... Compute the first order Corner Transport Upwind correction (step 2 in LeVeque algorithm)
!
      USE dimensions
      USE domain_mapping, ONLY: myijk
      USE grid, ONLY: rb, dx, dz
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: stencil
      USE time_parameters, ONLY: dt
!
      IMPLICIT NONE
!
      REAL*8, INTENT(INOUT) :: fe, ft, fw, fb

      TYPE(stencil), INTENT(IN) :: u, w, dens
      INTEGER, INTENT(IN) :: ij
      INTEGER :: i, k, imesh

      REAL*8 :: uc, wc, deltar
      REAL*8 :: dxp, dxm, dzp, dzm, indxp, indxm, indzp, indzm
      REAL*8 :: dxi, dxip, dxim, dzk, dzkp, dzkm
      REAL*8 :: incrxp, incrxm, incrzp, incrzm
!
      imesh = myijk( ip0_jp0_kp0_, ij)
      k = ( imesh - 1 ) / nx + 1
      i = MOD( ( imesh - 1 ), nx) + 1
!
      dxi  = dx(i)
      dxip = dx(i+1)
      dxim = dx(i-1)
      dzk  = dz(k)
      dzkp = dz(k+1)
      dzkm = dz(k-1)
!
      dxm  = dxi + dxim
      dxp  = dxi + dxip
      dzm  = dzk + dzkm
      dzp  = dzk + dzkp
!
      indxm  = 2.D0/dxm
      indxp  = 2.D0/dxp
      indzm  = 2.D0/dzm
      indzp  = 2.D0/dzp
!
      incrxp = 0.125D0*dt*indxp
      incrxm = 0.125D0*dt*indxm
      incrzp = 0.125D0*dt*indzp
      incrzm = 0.125D0*dt*indzm
!
! c-w
      uc = u%w
      wc = 0.25D0*indxm*(dxim*w%c+dxi*w%w+dxim*w%b+dxi*w%wb)
      deltar  = dens%c - dens%w
      ft = ft - incrxm*(uc+ABS(uc))*(wc+ABS(wc))*deltar
      fb = fb - incrxm*(uc+ABS(uc))*(wc-ABS(wc))*deltar
!
! e-c
      uc = u%c
      wc = 0.25D0*indxp*(dxip*w%c+dxi*w%e+dxip*w%b+dxi*w%eb)
      deltar  = dens%e - dens%c
      ft = ft - incrxp*(uc-ABS(uc))*(wc+ABS(wc))*deltar
      fb = fb - incrxp*(uc-ABS(uc))*(wc-ABS(wc))*deltar
!
! t-wt
      uc = u%wt
      wc = 0.25D0*indxm*(dxim*w%c+dxi*w%w+dxim*w%t+dxi*w%wt)
      deltar  = dens%t - dens%wt
      IF (k >= nz-1) deltar = 0.D0
      ft = ft - incrxm*(uc+ABS(uc))*(wc-ABS(wc))*deltar
!
! et-t
      uc = u%t
      wc = 0.25D0*indxp*(dxip*w%c+dxi*w%e+dxip*w%t+dxi*w%et)
      deltar  = dens%et - dens%t
      IF (k >= nz-1) deltar = 0.D0
      ft = ft - incrxp*(uc-ABS(uc))*(wc-ABS(wc))*deltar
!
! b-wb
      uc = u%wb
      wc = 0.25D0*indxm*(dxim*w%b+dxi*w%wb+dxim*w%bb+dxi*w%wbb)
      deltar  = dens%b - dens%wb
      IF (k <= 2) deltar = 0.D0
      fb = fb - incrxm*(uc+ABS(uc))*(wc+ABS(wc))*deltar
! eb-b
      uc = u%b
      wc = 0.25D0*indxp*(dxip*w%b+dxi*w%eb+dxip*w%bb+dxi*w%ebb)
      deltar = dens%eb - dens%b
      IF (k <= 2) deltar = 0.D0
      fb = fb - incrxp*(uc-ABS(uc))*(wc+ABS(wc))*deltar
!
! c-b
      uc = 0.25D0*indzm*(dzkm*u%c+dzk*u%b+dzkm*u%w+dzk*u%wb)
      wc = w%b
      deltar  = dens%c - dens%b
      fe = fe - incrzm*(uc+ABS(uc))*(wc+ABS(wc))*deltar*rb(i)
      fw = fw - incrzm*(uc-ABS(uc))*(wc+ABS(wc))*deltar*rb(i-1)
!
! t-c
      uc = 0.25D0*indzp*(dzkp*u%c+dzk*u%t+dzkp*u%w+dzk*u%wt)
      wc = w%c
      deltar  = dens%t - dens%c
      fe = fe - incrzp*(uc+ABS(uc))*(wc-ABS(wc))*deltar*rb(i)
      fw = fw - incrzp*(uc-ABS(uc))*(wc-ABS(wc))*deltar*rb(i-1)
!
! wt-w
      uc = 0.25D0*indzp*(dzkp*u%w+dzk*u%wt+dzkp*u%ww+dzk*u%wwt)
      wc = w%w
      deltar  = dens%wt - dens%w
      IF (i <= 2) deltar = 0.D0
      fw = fw - incrzp*(uc+ABS(uc))*(wc-ABS(wc))*deltar*rb(i-1)
!
! w-wb
      uc = 0.25D0*indzm*(dzkm*u%w+dzk*u%wb+dzkm*u%ww+dzk*u%wwb)
      wc = w%wb
      deltar  = dens%w - dens%wb
      IF (i <= 2) deltar = 0.D0
      fw = fw - incrzm*(uc+ABS(uc))*(wc+ABS(wc))*deltar*rb(i-1)
!
! et-e
      uc = 0.25D0*indzp*(dzkp*u%c+dzk*u%t+dzkp*u%e+dzk*u%et)
      wc = w%e
      deltar = dens%et - dens%e
      IF (i >= nx-1) deltar = 0.D0
      fe = fe - incrzp*(uc-ABS(uc))*(wc-ABS(wc))*deltar*rb(i)
!
! e-eb
      uc = 0.25D0*indzm*(dzkm*u%c+dzk*u%b+dzkm*u%e+dzk*u%eb)
      wc = w%eb
      deltar  = dens%e - dens%eb
      IF (i >= nx-1) deltar = 0.D0
      fe = fe - incrzm*(uc-ABS(uc))*(wc+ABS(wc))*deltar*rb(i)
!
      RETURN
      END SUBROUTINE ctu1_fmas_2d
!----------------------------------------------------------------
      SUBROUTINE ctu2_fmas_2d(fe, ft, fw, fb, dens, u, w, ijk)
!
! ... Compute the second order Corner Transport Upwind correction (step 3 in LeVeque algorithm)
!
      USE dimensions
      USE domain_mapping, ONLY: myijk
      USE grid, ONLY: rb, dx, dz, indx, indz, fluid, flag
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: stencil
      USE set_indexes, ONLY: imjk, ijkm, ipjk, ijkp
      USE time_parameters, ONLY: dt
!
      IMPLICIT NONE
!
      REAL*8, INTENT(INOUT) :: fe, ft, fw, fb
      TYPE(stencil), INTENT(IN) :: dens, u, w 
      INTEGER, INTENT(IN) :: ijk 
!
      INTEGER :: i, k, imesh, ip2, im2, kp2, km2
      REAL*8 :: uref, wref, deltar, erre, lim
      REAL*8 :: dxm, dxp, dxpp, dxmm, dzp, dzm, dzpp, dzmm
      REAL*8 :: indxm, indxp, indxpp, indxmm, indzp, indzm, indzpp, indzmm
      REAL*8 :: gradc, gradw, grade, gradt, gradb
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
      dxmm = dx(i-1)+ dx(im2)
      dxm  = dx(i)  + dx(i-1)
      dxp  = dx(i)  + dx(i+1)
      dxpp = dx(i+1)+ dx(ip2)
      dzmm = dz(k-1)+ dz(km2)
      dzm  = dz(k)  + dz(k-1)
      dzp  = dz(k)  + dz(k+1)
      dzpp = dz(k+1)+ dz(kp2)
!
      indxm  = 2.D0/dxm
      indxmm = 2.D0/dxmm
      indxp  = 2.D0/dxp
      indxpp = 2.D0/dxpp
      indzm  = 2.D0/dzm
      indzmm = 2.D0/dzmm
      indzp  = 2.D0/dzp
      indzpp = 2.D0/dzpp
!
! e-c
      gradc = indxp  * (dens%e  - dens%c)
      grade = indxpp * (dens%ee - dens%e)
      gradw = indxm  * (dens%c  - dens%w)
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
      deltar = dens%e - dens%c
      fe = fe + 0.5D0*uref*(1-dt*indxp*uref)*deltar*lim*rb(i)
!
! c-w
      gradc = indxm  * (dens%c - dens%w )
      grade = indxp  * (dens%e - dens%c )
      gradw = indxmm * (dens%w - dens%ww)
      lim  = 0.D0
      erre = 0.D0
      uref = u%w
      IF (uref >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=2 .AND. flag(imjk)==fluid) CALL limiters(lm,lim,erre)
      uref = ABS(u%w)
      deltar = dens%c - dens%w
      fw = fw + 0.5D0*uref*(1-dt*indxm*uref)*deltar*lim*rb(i-1)
!
! t-c
      gradc = indzp  * (dens%t  - dens%c)
      gradt = indzpp * (dens%tt - dens%t)
      gradb = indzm  * (dens%c  - dens%b)
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
      deltar = dens%t - dens%c
      ft = ft + 0.5D0*wref*(1-dt*indzp*wref)*deltar*lim
!
! c-b
      gradc = indzm  * (dens%c - dens%b )
      gradt = indzp  * (dens%t - dens%c )
      gradb = indzmm * (dens%b - dens%bb)
      lim  = 0.D0
      erre = 0.D0
      wref = w%b
      IF (wref >= 0.D0) THEN
          IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
          IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=2 .AND. flag(ijkm)==fluid) CALL limiters(lm,lim,erre)
      wref = ABS(w%b)
      deltar = dens%c - dens%b
      fb = fb + 0.5D0*wref*(1-dt*indzm*wref)*deltar*lim
!
      RETURN
      END SUBROUTINE ctu2_fmas_2d
!-----------------------------------------------------------------------
      END MODULE convective_mass_fluxes
!-----------------------------------------------------------------------
