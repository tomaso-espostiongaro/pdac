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
        MODULE PROCEDURE ctu1_fmas_2d, ctu1_fmas_3d
      END INTERFACE
      INTERFACE ctu2_fmas
        MODULE PROCEDURE ctu2_fmas_2d, ctu2_fmas_3d
      END INTERFACE
      INTERFACE ctu3_fmas
        MODULE PROCEDURE ctu3_fmas_2d, ctu3_fmas_3d
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
    SUBROUTINE ctu3_fmas_2d(fe, ft, fw, fb, dens, u, w, ijk)
!
! ... Compute the second order Corner Transport Upwind correction (step 4 in LeVeque algorithm)
!
      USE dimensions
      USE domain_mapping, ONLY: myijk
      USE grid, ONLY: rb, dx, dz, indx, indz, fluid, flag
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: stencil, ipjk, imjk, ijkp, ijkm
      USE set_indexes, ONLY: ipjkp, ipjkm, imjkp, imjkm
      USE time_parameters, ONLY: dt
!
      IMPLICIT NONE
!
      REAL*8, INTENT(INOUT) :: fe, ft, fw, fb

      TYPE(stencil), INTENT(IN) :: u, w, dens
      INTEGER, INTENT(IN) :: ijk
      INTEGER :: i, k, imesh

      REAL*8 :: uc, wc, deltar, lim, erre
      REAL*8 :: gradc, grade, gradw, gradt, gradb
      REAL*8 :: dxp, dxm, dzp, dzm, dxpp, dxmm, dzpp, dzmm
      REAL*8 :: indxp, indxm, indzp, indzm, indxpp, indxmm, indzpp, indzmm
      REAL*8 :: incrxm, incrxp, incrxpp, incrzm, incrzp, incrzpp
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
      incrxm  = 0.5D0*dt*indx(i-1)
      incrxp  = 0.5D0*dt*indx(i)
      incrxpp = 0.5D0*dt*indx(i+1)
      incrzm  = 0.5D0*dt*indz(k-1)
      incrzp  = 0.5D0*dt*indz(k)
      incrzpp = 0.5D0*dt*indz(k+1)
!
! c-w
      uc = u%w
      wc = 0.25D0*indxm*(dx(i-1)*w%c+dx(i)*w%w+dx(i-1)*w%b+dx(i)*w%wb)
      deltar = 0.5D0*ABS(uc)*(1 - dt*indxm*ABS(uc))*(dens%c - dens%w)
      gradc = indxm  * (dens%c - dens%w )
      grade = indxp  * (dens%e - dens%c )
      gradw = indxmm * (dens%w - dens%ww)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=2 .AND. flag(imjk)==fluid) CALL limiters(lm,lim,erre)
      ft = ft + incrzp*(wc+ABS(wc))*deltar*lim
      fb = fb + incrzp*(wc-ABS(wc))*deltar*lim
!
! e-c 
      uc = u%c
      wc = 0.25D0*indxp*(dx(i+1)*w%c+dx(i)*w%e+dx(i+1)*w%b+dx(i)*w%eb)
      deltar = 0.5D0*ABS(uc)*(1 - dt*indxp*ABS(uc))*(dens%e - dens%c)
      gradc = indxp  * (dens%e  - dens%c)
      grade = indxpp * (dens%ee - dens%e)
      gradw = indxm  * (dens%c  - dens%w)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=nx-1 .AND. flag(ipjk)==fluid) CALL limiters(lm,lim,erre)
      ft = ft - incrzp*(wc+ABS(wc))*deltar*lim
      fb = fb - incrzp*(wc-ABS(wc))*deltar*lim
!
! t-wt
      uc = u%wt
      wc = 0.25D0*indxm*(dx(i-1)*w%c+dx(i)*w%w+dx(i-1)*w%t+dx(i)*w%wt)
      deltar = 0.5D0*ABS(uc)*(1 - dt*indxm*ABS(uc))*(dens%t - dens%wt)
      IF (k >= nz-1) deltar = 0.D0
      gradc = indxm  * (dens%t  - dens%wt )
      grade = indxp  * (dens%et - dens%t  )
      gradw = indxmm * (dens%wt - dens%wwt)
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
      deltar = 0.5D0*ABS(uc)*(1 - dt*indxp*ABS(uc))*(dens%et - dens%t)
      IF (k >= nz-1) deltar = 0.D0
      gradc = indxp  * (dens%et  - dens%t )
      grade = indxpp * (dens%eet - dens%et)
      gradw = indxm  * (dens%t   - dens%wt)
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
! b-wb
      uc = u%wb
      wc = 0.25D0*indxm*(dx(i-1)*w%b+dx(i)*w%wb+dx(i-1)*w%bb+dx(i)*w%wbb)
      deltar = 0.5D0*ABS(uc)*(1 - dt*indxm*ABS(uc))*(dens%b - dens%wb)
      IF (k <= 2) deltar = 0.D0
      gradc = indxm  * (dens%b - dens%wb  )
      grade = indxp  * (dens%eb - dens%b  )
      gradw = indxmm * (dens%wb - dens%wwb)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=2 .AND. flag(imjkm)==fluid) CALL limiters(lm,lim,erre)
      fb = fb + incrzm*(wc+ABS(wc))*deltar*lim
!
! eb-b
      uc = u%b
      wc = 0.25D0*indxp*(dx(i+1)*w%b+dx(i)*w%eb+dx(i+1)*w%bb+dx(i)*w%ebb)
      deltar = 0.5D0*ABS(uc)*(1 - dt*indxp*ABS(uc))*(dens%eb - dens%b)
      IF (k <= 2) deltar = 0.D0
      gradc = indxp  * (dens%eb - dens%b  )
      grade = indxpp * (dens%eeb - dens%eb)
      gradw = indxm  * (dens%b - dens%wb  )
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=nx-1 .AND. flag(ipjkm)==fluid) CALL limiters(lm,lim,erre)
      fb = fb - incrzm*(wc+ABS(wc))*deltar*lim
!
! c-b
      uc = 0.25D0*indzm*(dz(k-1)*u%c+dz(k)*u%b+dz(k-1)*u%w+dz(k)*u%wb)
      wc = w%b
      deltar = 0.5D0*ABS(wc)*(1 - dt*indzm*ABS(wc))*(dens%c - dens%b)
      gradc = indzm  * (dens%c - dens%b )
      gradt = indzp  * (dens%t - dens%c )
      gradb = indzmm * (dens%b - dens%bb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=2 .AND. flag(ijkm)==fluid) CALL limiters(lm,lim,erre)
      fe = fe + incrxp*(uc+ABS(uc))*deltar*lim*rb(i)
      fw = fw + incrxp*(uc-ABS(uc))*deltar*lim*rb(i-1)
!
! t-c
      uc = 0.25D0*indzp*(dz(k+1)*u%c+dz(k)*u%t+dz(k+1)*u%w+dz(k)*u%wt)
      wc = w%c
      deltar = 0.5D0*ABS(wc)*(1 - dt*indzp*ABS(wc))*(dens%t - dens%c)
      gradc = indzp  * (dens%t  - dens%c)
      gradt = indzpp * (dens%tt - dens%t)
      gradb = indzm  * (dens%c  - dens%b)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(ijkp)==fluid) CALL limiters(lm,lim,erre)
      fe = fe - incrxp*(uc+ABS(uc))*deltar*lim*rb(i)
      fw = fw - incrxp*(uc-ABS(uc))*deltar*lim*rb(i-1)
!
! wt-w 
      uc = 0.25D0*indzp*(dz(k+1)*u%w+dz(k)*u%wt+dz(k+1)*u%ww+dz(k)*u%wwt)
      wc = w%w
      deltar = 0.5D0*ABS(wc)*(1 - dt*indzp*ABS(wc))*(dens%wt - dens%w)
      IF (i <= 2) deltar = 0.D0
      gradc = indzp  * (dens%wt  - dens%w )
      gradt = indzpp * (dens%wtt - dens%wt)
      gradb = indzm  * (dens%w   - dens%wb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(imjkp)==fluid) CALL limiters(lm,lim,erre)
      fw = fw - incrxm*(uc+ABS(uc))*deltar*lim*rb(i-1)
!
! w-wb
      uc = 0.25D0*indzm*(dz(k-1)*u%w+dz(k)*u%wb+dz(k-1)*u%ww+dz(k)*u%wwb)
      wc = w%wb
      deltar = 0.5D0*ABS(wc)*(1 - dt*indzm*ABS(wc))*(dens%w - dens%wb)
      IF (i <= 2) deltar = 0.D0
      gradc = indzm  * (dens%w  - dens%wb )
      gradt = indzp  * (dens%wt - dens%w  )
      gradb = indzmm * (dens%wb - dens%wbb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=2 .AND. flag(imjkm)==fluid) CALL limiters(lm,lim,erre)
      fw = fw + incrxm*(uc+ABS(uc))*deltar*lim*rb(i-1)
!
! et-e
      uc = 0.25D0*indzp*(dz(k+1)*u%c+dz(k)*u%t+dz(k+1)*u%e+dz(k)*u%et)
      wc = w%e
      deltar = 0.5D0*ABS(wc)*(1 - dt*indzp*ABS(wc))*(dens%et - dens%e)
      IF (i >= nx-1) deltar = 0.D0
      gradc = indzp  * (dens%et  - dens%e )
      gradt = indzpp * (dens%ett - dens%et)
      gradb = indzm  * (dens%e   - dens%eb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(ipjkp)==fluid) CALL limiters(lm,lim,erre)
      fe = fe - incrxpp*(uc-ABS(uc))*deltar*lim*rb(i)

! e-eb
      uc = 0.25D0*indzm*(dz(k-1)*u%c+dz(k)*u%b+dz(k-1)*u%e+dz(k)*u%eb)
      wc = w%eb
      deltar = 0.5D0*ABS(wc)*(1 - dt*indzm*ABS(wc))*(dens%e - dens%eb)
      IF (i >= nx-1) deltar = 0.D0
      gradc = indzm  * (dens%e  - dens%eb )
      gradt = indzp  * (dens%et - dens%e  )
      gradb = indzmm * (dens%eb - dens%ebb)
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
      END SUBROUTINE ctu3_fmas_2d
!-----------------------------------------------------------------------
      SUBROUTINE ctu1_fmas_3d( fe, fn, ft, fw, fs, fb, dens, u, v, w, ijk )
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

      REAL*8, INTENT(INOUT) :: fe, fn, ft, fw, fs, fb

      TYPE(stencil), INTENT(IN) :: u, v, w, dens
      INTEGER, INTENT(IN) :: ijk
      INTEGER :: i, j, k, imesh

      REAL*8 :: uc, vc, wc, deltar, onethird
      REAL*8 :: dxp, dxm, dyp, dym, dzp, dzm
      REAL*8 :: dxi, dxip, dxim, dyj, dyjp, dyjm, dzk, dzkp, dzkm
      REAL*8 :: incrxp, incrxm, incryp, incrym, incrzp, incrzm
      REAL*8 :: indxp, indxm, indyp, indym, indzp, indzm
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
      indxm = 2.D0/dxm
      indxp = 2.D0/dxp
      indym = 2.D0/dym
      indyp = 2.D0/dyp
      indzm = 2.D0/dzm
      indzp = 2.D0/dzp
!
      incrxp = 0.25D0*dt*indxp
      incrxm = 0.25D0*dt*indxm
      incryp = 0.25D0*dt*indyp
      incrym = 0.25D0*dt*indym
      incrzp = 0.25D0*dt*indzp
      incrzm = 0.25D0*dt*indzm
!
! t-c 
      deltar  = dens%t - dens%c
      uc = (dzkp*u%c+dzk*u%t+dzkp*u%w+dzk*u%wt)*indx(i)*incrzp
      vc = (dzkp*v%c+dzk*v%t+dzkp*v%s+dzk*v%st)*indy(j)*incrzp
      wc = w%c
      fw = fw - 0.125D0*(uc-ABS(uc))*(wc-ABS(wc))*deltar 
      fe = fe - 0.125D0*(uc+ABS(uc))*(wc-ABS(wc))*deltar
      fn = fn - 0.125D0*(vc+ABS(vc))*(wc-ABS(wc))*deltar &
              + 0.25D0*onethird*(vc+ABS(vc))*(wc-ABS(wc))*ABS(uc)*deltar
      fs = fs - 0.125D0*(vc-ABS(vc))*(wc-ABS(wc))*deltar &
              + 0.25D0*onethird*(vc-ABS(vc))*(wc-ABS(wc))*ABS(uc)*deltar
!
! c-b
      deltar  = dens%c - dens%b
      uc = (dzkm*u%c+dzk*u%b+dzkm*u%w+dzk*u%wb)*indx(i)*incrzm
      vc = (dzkm*v%c+dzk*v%b+dzkm*v%s+dzk*v%sb)*indy(j)*incrzm
      wc = w%b
      fw = fw - 0.125D0*(uc-ABS(uc))*(wc+ABS(wc))*deltar
      fe = fe - 0.125D0*(uc+ABS(uc))*(wc+ABS(wc))*deltar
      fn = fn - 0.125D0*(vc+ABS(vc))*(wc+ABS(wc))*deltar &
              + 0.25D0*onethird*(vc+ABS(vc))*(wc+ABS(wc))*ABS(uc)*deltar 
      fs = fs - 0.125D0*(vc-ABS(vc))*(wc+ABS(wc))*deltar &
              + 0.25D0*onethird*(vc-ABS(vc))*(wc+ABS(wc))*ABS(uc)*deltar
!
! n-c
      deltar  = dens%n - dens%c
      uc = (dyjp*u%c+dyj*u%n+dyjp*u%w+dyj*u%wn)*indx(i)*incryp
      vc = v%c
      wc = (dyjp*w%c+dyj*w%n+dyjp*w%b+dyj*w%nb)*indz(k)*incryp
      ft = ft - 0.125D0*(wc+ABS(wc))*(vc-ABS(vc))*deltar
      fb = fb - 0.125D0*(wc-ABS(wc))*(vc-ABS(vc))*deltar
      fe = fe - 0.125D0*(uc+ABS(uc))*(vc-ABS(vc))*deltar &
              + 0.25D0*onethird*(uc+ABS(uc))*(vc-ABS(vc))*ABS(wc)*deltar 
      fw = fw - 0.125D0*(uc-ABS(uc))*(vc-ABS(vc))*deltar &
              + 0.25D0*onethird*(uc-ABS(uc))*(vc-ABS(vc))*ABS(wc)*deltar
!
! c-s
      deltar  = dens%c - dens%s
      uc = (dyjm*u%c+dyj*u%s+dyjm*u%w+dyj*u%ws)*indx(i)*incrym
      vc = v%s
      wc = (dyjm*w%c+dyj*w%s+dyjm*w%b+dyj*w%sb)*indz(k)*incrym
      ft = ft - 0.125D0*(wc+ABS(wc))*(vc+ABS(vc))*deltar
      fb = fb - 0.125D0*(wc-ABS(wc))*(vc+ABS(vc))*deltar
      fe = fe - 0.125D0*(uc+ABS(uc))*(vc+ABS(vc))*deltar &
              + 0.25D0*onethird*(uc+ABS(uc))*(vc+ABS(vc))*ABS(wc)*deltar 
      fw = fw - 0.125D0*(uc-ABS(uc))*(vc+ABS(vc))*deltar &
              + 0.25D0*onethird*(uc-ABS(uc))*(vc+ABS(vc))*ABS(wc)*deltar
!
! e-c
      deltar  = dens%e - dens%c
      uc = u%c
      vc = (dxip*v%c+dxi*v%e+dxip*v%s+dxi*v%es)*indy(j)*incrxp
      wc = (dxip*w%c+dxi*w%e+dxip*w%b+dxi*w%eb)*indz(k)*incrxp
      fn = fn - 0.125D0*(vc+ABS(vc))*(uc-ABS(uc))*deltar
      fs = fs - 0.125D0*(vc-ABS(vc))*(uc-ABS(uc))*deltar
      ft = ft - 0.125D0*(wc+ABS(wc))*(uc-ABS(uc))*deltar &
              + 0.25D0*onethird*(wc+ABS(wc))*(uc-ABS(uc))*ABS(vc)*deltar 
      fb = fb - 0.125D0*(wc-ABS(wc))*(uc-ABS(uc))*deltar &
              + 0.25D0*onethird*(wc-ABS(wc))*(uc-ABS(uc))*ABS(vc)*deltar
!
! c-w
      deltar  = dens%c - dens%w
      uc = u%w
      vc = (dxim*v%c+dxi*v%w+dxim*v%s+dxi*v%ws)*indy(j)*incrxm
      wc = (dxim*w%c+dxi*w%w+dxim*w%b+dxi*w%wb)*indz(k)*incrxm
      fn = fn - 0.125D0*(vc+ABS(vc))*(uc+ABS(uc))*deltar
      fs = fs - 0.125D0*(vc-ABS(vc))*(uc+ABS(uc))*deltar
      ft = ft - 0.125D0*(wc+ABS(wc))*(uc+ABS(uc))*deltar &
              + 0.25D0*onethird*(wc+ABS(wc))*(uc+ABS(uc))*ABS(vc)*deltar 
      fb = fb - 0.125D0*(wc-ABS(wc))*(uc+ABS(uc))*deltar &
              + 0.25D0*onethird*(wc-ABS(wc))*(uc+ABS(uc))*ABS(vc)*deltar
!
! wt-w 
      deltar  = dens%wt - dens%w
      IF (i <= 2) deltar = 0.D0
      uc = (dzkp*u%w+dzk*u%wt+dzkp*u%ww+dzk*u%wwt)*indx(i-1)*incrzp
      vc = (dzkp*v%w+dzk*v%wt+dzkp*v%ws+dzk*v%wst)*indy(j)  *incrzp
      wc = w%w
      fw = fw - 0.125D0*(wc-ABS(wc))*(uc+ABS(uc))*deltar
      fn = fn - 0.125D0*onethird*(wc-ABS(wc))*(vc+ABS(vc))*ABS(uc+ABS(uc))*deltar
      fs = fs - 0.125D0*onethird*(wc-ABS(wc))*(vc-ABS(vc))*ABS(uc+ABS(uc))*deltar
!
! w-wb 
      deltar  = dens%w - dens%wb
      IF (i <= 2) deltar = 0.D0
      uc = (dzkm*u%w+dzk*u%wb+dzkm*u%ww+dzk*u%wwb)*indx(i-1)*incrzm
      vc = (dzkm*v%w+dzk*v%wb+dzkm*v%ws+dzk*v%wsb)*indy(j)  *incrzm
      wc = w%wb
      fw = fw - 0.125D0*(wc+ABS(wc))*(uc+ABS(uc))*deltar
      fn = fn - 0.125D0*onethird*(wc+ABS(wc))*(vc+ABS(vc))*ABS(uc+ABS(uc))*deltar
      fs = fs - 0.125D0*onethird*(wc+ABS(wc))*(vc-ABS(vc))*ABS(uc+ABS(uc))*deltar
!
! et-e 
      deltar  = dens%et - dens%e
      IF (i >= nx-1) deltar = 0.D0
      uc = (dzkp*u%e+dzk*u%et+dzkp*u%c +dzk*u%t  )*indx(i+1)*incrzp
      vc = (dzkp*v%e+dzk*v%et+dzkp*v%es+dzk*v%est)*indy(j)  *incrzp
      wc = w%e
      fe = fe - 0.125D0*(wc-ABS(wc))*(uc-ABS(uc))*deltar
      fn = fn - 0.125D0*onethird*(wc-ABS(wc))*(vc+ABS(vc))*ABS(uc-ABS(uc))*deltar
      fs = fs - 0.125D0*onethird*(wc-ABS(wc))*(vc-ABS(vc))*ABS(uc-ABS(uc))*deltar
!
! e-eb
      deltar  = dens%e - dens%eb
      IF (i >= nx-1) deltar = 0.D0
      uc = (dzkm*u%e+dzk*u%eb+dzkm*u%c +dzk*u%b  )*indx(i+1)*incrzm
      vc = (dzkm*v%e+dzk*v%eb+dzkm*v%es+dzk*v%esb)*indy(j)  *incrzm
      wc = w%eb
      fe = fe - 0.125D0*(wc+ABS(wc))*(uc-ABS(uc))*deltar
      fn = fn - 0.125D0*onethird*(wc+ABS(wc))*(vc+ABS(vc))*ABS(uc-ABS(uc))*deltar
      fs = fs - 0.125D0*onethird*(wc+ABS(wc))*(vc-ABS(vc))*ABS(uc-ABS(uc))*deltar
!
! wn-w
      deltar  = dens%wn - dens%w
      IF (i <= 2) deltar = 0.D0
      uc = (dyjp*u%w+dyj*u%wn+dyjp*u%ww+dyj*u%wwn)*indx(i-1)*incryp
      vc = v%w
      wc = (dyjp*w%w+dyj*w%wn+dyjp*w%wb+dyj*w%wnb)*indz(k)  *incryp
      fw = fw - 0.125D0*(vc-ABS(vc))*(uc+ABS(uc))*deltar &
              + 0.25D0*onethird*(vc-ABS(vc))*(uc+ABS(uc))*ABS(wc)*deltar
!
! w-ws
      deltar  = dens%w - dens%ws
      IF (i <= 2) deltar = 0.D0
      uc = (dyjm*u%w+dyj*u%ws+dyjm*u%ww+dyj*u%wws)*indx(i-1)*incrym
      vc = v%ws
      wc = (dyjm*w%w+dyj*w%ws+dyjm*w%wb+dyj*w%wsb)*indz(k)  *incrym
      fw = fw - 0.125D0*(vc+ABS(vc))*(uc+ABS(uc))*deltar &
              + 0.25D0*onethird*(vc+ABS(vc))*(uc+ABS(uc))*ABS(wc)*deltar
!
! en-e 
      deltar  = dens%en - dens%e
      IF (i >= nx-1) deltar = 0.D0
      uc = (dyjp*u%e+dyj*u%en+dyjp*u%c +dyj*u%n  )*indx(i+1)*incryp
      vc = v%e
      wc = (dyjp*w%e+dyj*w%en+dyjp*w%eb+dyj*w%enb)*indz(k)  *incryp
      fe = fe - 0.125D0*(vc-ABS(vc))*(uc-ABS(uc))*deltar &
              + 0.25D0*onethird*(vc-ABS(vc))*(uc-ABS(uc))*ABS(wc)*deltar
!
! e-es 
      deltar  = dens%e - dens%es
      IF (i >= nx-1) deltar = 0.D0
      uc = (dyjm*u%e+dyj*u%es+dyjm*u%c +dyj*u%s  )*indx(i+1)*incrym
      vc = v%es
      wc = (dyjm*w%e+dyj*w%es+dyjm*w%eb+dyj*w%esb)*indz(k)  *incrym
      fe = fe - 0.125D0*(vc+ABS(vc))*(uc-ABS(uc))*deltar &
              + 0.25D0*onethird*(vc+ABS(vc))*(uc-ABS(uc))*ABS(wc)*deltar
!
! st-s
      deltar  = dens%st - dens%s
      IF (j <= 2) deltar = 0.D0
      uc = (dzkp*u%s+dzk*u%st+dzkp*u%ws+dzk*u%wst)*indx(i)  *incrzp
      vc = (dzkp*v%s+dzk*v%st+dzkp*v%ss+dzk*v%sst)*indy(j-1)*incrzp
      wc = w%s
      fs = fs - 0.125D0*(vc+ABS(vc))*(wc-ABS(wc))*deltar &
              + 0.25D0*onethird*(vc+ABS(vc))*(wc-ABS(wc))*ABS(uc)*deltar
!
! s-sb
      deltar  = dens%s - dens%sb
      IF (j <= 2) deltar = 0.D0
      uc = (dzkm*u%s+dzk*u%sb+dzkm*u%ws+dzk*u%wsb)*indx(i)  *incrzm
      vc = (dzkm*v%s+dzk*v%sb+dzkm*v%ss+dzk*v%ssb)*indy(j-1)*incrzm
      wc = w%sb
      fs = fs - 0.125D0*(vc+ABS(vc))*(wc+ABS(wc))*deltar &
              + 0.25D0*onethird*(vc+ABS(vc))*(wc+ABS(wc))*ABS(uc)*deltar
!
! nt-n
      deltar  = dens%nt - dens%n
      IF (j >= ny-1) deltar = 0.D0
      uc = (dzkp*u%n+dzk*u%nt+dzkp*u%wn+dzk*u%wnt)*indx(i)  *incrzp
      vc = (dzkp*v%n+dzk*v%nt+dzkp*v%c +dzk*v%t  )*indy(j+1)*incrzp
      wc = w%n
      fn = fn - 0.125D0*(vc-ABS(vc))*(wc-ABS(wc))*deltar &
              + 0.25D0*onethird*(vc-ABS(vc))*(wc-ABS(wc))*ABS(uc)*deltar
!
! n-nb 
      deltar  = dens%n - dens%nb
      IF (j >= ny-1) deltar = 0.D0
      uc = (dzkm*u%n+dzk*u%nb+dzkm*u%wn+dzk*u%wnb)*indx(i)  *incrzm
      vc = (dzkm*v%n+dzk*v%nb+dzkm*v%c +dzk*v%b  )*indy(j+1)*incrzm
      wc = w%nb
      fn = fn - 0.125D0*(vc-ABS(vc))*(wc+ABS(wc))*deltar &
              + 0.25D0*onethird*(vc-ABS(vc))*(wc+ABS(wc))*ABS(uc)*deltar
!
! nt-t
      deltar  = dens%nt - dens%t
      IF (k >= nz-1) deltar = 0.D0
      uc = (dyjp*u%t+dyj*u%nt+dyjp*u%wt+dyj*u%wnt)*indx(i)  *incryp
      vc = v%t
      wc = (dyjp*w%t+dyj*w%nt+dyjp*w%c +dyj*w%n  )*indz(k+1)*incryp
      ft = ft - 0.125D0*(vc-ABS(vc))*(wc-ABS(wc))*deltar
      fe = fe - 0.125D0*onethird*(vc-ABS(vc))*(uc+ABS(uc))*ABS(wc-ABS(wc))*deltar
      fw = fw - 0.125D0*onethird*(vc-ABS(vc))*(uc-ABS(uc))*ABS(wc-ABS(wc))*deltar
!
! t-st 
      deltar  = dens%t - dens%st
      IF (k >= nz-1) deltar = 0.D0
      uc = (dyjm*u%t+dyj*u%st+dyjm*u%wt+dyj*u%wst)*indx(i)  *incrym
      vc = v%st
      wc = (dyjm*w%t+dyj*w%st+dyjm*w%c +dyj*w%s  )*indz(k+1)*incrym
      ft = ft - 0.125D0*(vc+ABS(vc))*(wc-ABS(wc))*deltar
      fe = fe - 0.125D0*onethird*(vc+ABS(vc))*(uc+ABS(uc))*ABS(wc-ABS(wc))*deltar
      fw = fw - 0.125D0*onethird*(vc+ABS(vc))*(uc-ABS(uc))*ABS(wc-ABS(wc))*deltar
!
! nb-b
      deltar  = dens%nb - dens%b
      IF (k <= 2) deltar = 0.D0
      uc = (dyjp*u%b+dyj*u%nb+dyjp*u%wb+dyj*u%wnb)*indx(i)  *incryp
      vc = v%b
      wc = (dyjp*w%b+dyj*w%nb+dyjp*w%bb+dyj*w%nbb)*indz(k-1)*incryp
      fb = fb - 0.125D0*(vc-ABS(vc))*(wc+ABS(wc))*deltar
      fe = fe - 0.125D0*onethird*(vc-ABS(vc))*(uc+ABS(uc))*ABS(wc+ABS(wc))*deltar
      fw = fw - 0.125D0*onethird*(vc-ABS(vc))*(uc-ABS(uc))*ABS(wc+ABS(wc))*deltar
!
! b-sb
      deltar  = dens%b - dens%sb
      IF (k <= 2) deltar = 0.D0
      uc = (dyjm*u%b+dyj*u%sb+dyjm*u%wb+dyj*u%wsb)*indx(i)  *incrym
      vc = v%sb
      wc = (dyjm*w%b+dyj*w%sb+dyjm*w%bb+dyj*w%sbb)*indz(k-1)*incrym
      fb = fb - 0.125D0*(vc+ABS(vc))*(wc+ABS(wc))*deltar
      fe = fe - 0.125D0*onethird*(vc+ABS(vc))*(uc+ABS(uc))*ABS(wc+ABS(wc))*deltar
      fw = fw - 0.125D0*onethird*(vc+ABS(vc))*(uc-ABS(uc))*ABS(wc+ABS(wc))*deltar
!
! en-n 
      deltar  = dens%en - dens%n
      IF (j >= ny-1) deltar = 0.D0
      uc = u%n
      vc = (dxip*v%n+dxi*v%en+dxip*v%c +dxi*v%e  )*indy(j+1)*incrxp
      wc = (dxip*w%n+dxi*w%en+dxip*w%nb+dxi*w%enb)*indz(k)  *incrxp
      fn = fn - 0.125D0*(vc-ABS(vc))*(uc-ABS(uc))*deltar
      ft = ft - 0.125D0*onethird*(uc-ABS(uc))*(wc+ABS(wc))*ABS(vc-ABS(vc))*deltar
      fb = fb - 0.125D0*onethird*(uc-ABS(uc))*(wc-ABS(wc))*ABS(vc-ABS(vc))*deltar
!
! n-wn
      deltar  = dens%n - dens%wn
      IF (j >= ny-1) deltar = 0.D0
      uc = u%wn
      vc = (dxim*v%n+dxi*v%wn+dxim*v%c +dxi*v%w  )*indy(j+1)*incrxm
      wc = (dxim*w%n+dxi*w%wn+dxim*w%nb+dxi*w%wnb)*indz(k)  *incrxm
      fn = fn - 0.125D0*(vc-ABS(vc))*(uc+ABS(uc))*deltar
      ft = ft - 0.125D0*onethird*(uc+ABS(uc))*(wc+ABS(wc))*ABS(vc-ABS(vc))*deltar
      fb = fb - 0.125D0*onethird*(uc+ABS(uc))*(wc-ABS(wc))*ABS(vc-ABS(vc))*deltar
!
! es-s
      deltar  = dens%es - dens%s
      IF (j <= 2) deltar = 0.D0
      uc = u%s
      vc = (dxip*v%s+dxi*v%es+dxip*v%ss+dxi*v%ess)*indy(j-1)*incrxp
      wc = (dxip*w%s+dxi*w%es+dxip*w%sb+dxi*w%esb)*indz(k)  *incrxp
      fs = fs - 0.125D0*(vc+ABS(vc))*(uc-ABS(uc))*deltar
      ft = ft - 0.125D0*onethird*(uc-ABS(uc))*(wc+ABS(wc))*ABS(vc+ABS(vc))*deltar
      fb = fb - 0.125D0*onethird*(uc-ABS(uc))*(wc-ABS(wc))*ABS(vc+ABS(vc))*deltar
!
! s-ws
      deltar  = dens%s - dens%ws
      IF (j <= 2) deltar = 0.D0
      uc = u%ws
      vc = (dxim*v%s+dxi*v%ws+dxim*v%ss+dxi*v%wss)*indy(j-1)*incrxm
      wc = (dxim*w%s+dxi*w%ws+dxim*w%sb+dxi*w%wsb)*indz(k)  *incrxm
      fs = fs - 0.125D0*(vc+ABS(vc))*(uc+ABS(uc))*deltar
      ft = ft - 0.125D0*onethird*(uc+ABS(uc))*(wc+ABS(wc))*ABS(vc+ABS(vc))*deltar
      fb = fb - 0.125D0*onethird*(uc+ABS(uc))*(wc-ABS(wc))*ABS(vc+ABS(vc))*deltar
!
! et-t
      deltar  = dens%et - dens%t
      IF (k >= nz-1) deltar = 0.D0
      uc = u%t
      vc = (dxip*v%t+dxi*v%et+dxip*v%st+dxi*v%est)*indy(j)  *incrxp
      wc = (dxip*w%t+dxi*w%et+dxip*w%c +dxi*w%e  )*indz(k+1)*incrxp
      ft = ft - 0.125D0*(wc-ABS(wc))*(uc-ABS(uc))*deltar &
              + 0.25D0*onethird*(wc-ABS(wc))*(uc-ABS(uc))*ABS(vc)*deltar
!
! t-wt
      deltar  = dens%t - dens%wt
      IF (k >= nz-1) deltar = 0.D0
      uc = u%wt
      vc = (dxim*v%t+dxi*v%wt+dxim*v%st+dxi*v%wst)*indy(j)  *incrxm
      wc = (dxim*w%t+dxi*w%wt+dxim*w%c +dxi*w%w  )*indz(k+1)*incrxm
      ft = ft - 0.125D0*(wc-ABS(wc))*(uc+ABS(uc))*deltar &
              + 0.25D0*onethird*(wc-ABS(wc))*(uc+ABS(uc))*ABS(vc)*deltar
!
! eb-b
      deltar  = dens%eb - dens%b
      IF (k <= 2) deltar = 0.D0
      uc = u%b
      vc = (dxip*v%b+dxi*v%eb+dxip*v%sb+dxi*v%esb)*indy(j)  *incrxp
      wc = (dxip*w%b+dxi*w%eb+dxip*w%bb+dxi*w%ebb)*indz(k-1)*incrxp
      fb = fb - 0.125D0*(wc+ABS(wc))*(uc-ABS(uc))*deltar &
              + 0.25D0*onethird*(wc+ABS(wc))*(uc-ABS(uc))*ABS(vc)*deltar
!
! b-wb
      deltar  = dens%b - dens%wb
      IF (k <= 2) deltar = 0.D0
      uc = u%wb
      vc = (dxim*v%b+dxi*v%wb+dxim*v%sb+dxi*v%wsb)*indy(j)  *incrxm
      wc = (dxim*w%b+dxi*w%wb+dxim*w%bb+dxi*w%wbb)*indz(k-1)*incrxm
      fb = fb - 0.125D0*(wc+ABS(wc))*(uc+ABS(uc))*deltar &
              + 0.25D0*onethird*(wc+ABS(wc))*(uc+ABS(uc))*ABS(vc)*deltar
!
! wt-wst
      deltar  = dens%wt - dens%wst
      IF (k >= nz-1) deltar = 0.D0
      IF (i <= 2) deltar = 0.D0
      uc = (dyjm*u%wt+dyj*u%wst+dyjm*u%wwt+dyj*u%wwst)*indx(i-1)*incrym
      vc = v%wst
      wc = (dyjm*w%wt+dyj*w%wst+dyjm*w%w  +dyj*w%ws  )*indz(k+1)*incrym
      fw = fw - 0.125D0*onethird*(vc+ABS(vc))*(uc+ABS(uc))*ABS(wc-ABS(wc))*deltar
!
! wnt-wt
      deltar  = dens%wnt - dens%wt
      IF (k >= nz-1) deltar = 0.D0
      IF (i <= 2) deltar = 0.D0
      uc = (dyjp*u%wt+dyj*u%wnt+dyjp*u%wwt+dyj*u%wwnt)*indx(i-1)*incryp
      vc = v%wt
      wc = (dyjp*w%wt+dyj*w%wnt+dyjp*w%w  +dyj*w%wn  )*indz(k+1)*incryp
      fw = fw - 0.125D0*onethird*(vc-ABS(vc))*(uc+ABS(uc))*ABS(wc-ABS(wc))*deltar
!
! wnb-wb
      deltar  = dens%wnb - dens%wb
      IF (i <= 2) deltar = 0.D0
      IF (k <= 2) deltar = 0.D0
      uc = (dyjp*u%wb+dyj*u%wnb+dyjp*u%wwb+dyj*u%wwnb)*indx(i-1)*incryp
      vc = v%wb
      wc = (dyjp*w%wb+dyj*w%wnb+dyjp*w%wbb+dyj*w%wnbb)*indz(k-1)*incryp
      fw = fw - 0.125D0*onethird*(vc-ABS(vc))*(uc+ABS(uc))*ABS(wc+ABS(wc))*deltar
!
! wb-wsb 
      deltar  = dens%wb - dens%wsb
      IF (i <= 2) deltar = 0.D0
      IF (k <= 2) deltar = 0.D0
      uc = (dyjm*u%wb+dyj*u%wsb+dyjm*u%wwb+dyj*u%wwsb)*indx(i-1)*incrym
      vc = v%wsb
      wc = (dyjm*w%wb+dyj*w%wsb+dyjm*w%wbb+dyj*w%wsbb)*indz(k-1)*incrym
      fw = fw - 0.125D0*onethird*(vc+ABS(vc))*(uc+ABS(uc))*ABS(wc+ABS(wc))*deltar
!
! wsb-ws 
      deltar  = dens%wst - dens%ws
      IF (i <= 2) deltar = 0.D0
      IF (j <= 2) deltar = 0.D0
      uc = (dzkp*u%ws+dzk*u%wst+dzkp*u%wws+dzk*u%wwst)*indx(i-1)*incrzp
      vc = (dzkp*v%ws+dzk*v%wst+dzkp*v%wss+dzk*v%wsst)*indy(j-1)*incrzp
      wc = w%ws
      fs = fs - 0.125D0*onethird*(wc-ABS(wc))*(vc+ABS(vc))*ABS(uc+ABS(uc))*deltar
!
! wnt-wn
      deltar  = dens%wnt - dens%wn
      IF (i <= 2) deltar = 0.D0
      IF (j >= ny-1) deltar = 0.D0
      uc = (dzkp*u%wn+dzk*u%wnt+dzkp*u%wwn+dzk*u%wwnt)*indx(i-1)*incrzp
      vc = (dzkp*v%wn+dzk*v%wnt+dzkp*v%w  +dzk*v%wt  )*indy(j+1)*incrzp
      wc = w%wn
      fn = fn - 0.125D0*onethird*(wc-ABS(wc))*(vc-ABS(vc))*ABS(uc+ABS(uc))*deltar
!
! ws-wsb
      deltar  = dens%ws - dens%wsb
      IF (i <= 2) deltar = 0.D0
      IF (j <= 2) deltar = 0.D0
      uc = (dzkm*u%ws+dzk*u%wsb+dzkm*u%wws+dzk*u%wwsb)*indx(i-1)*incrzm
      vc = (dzkm*v%ws+dzk*v%wsb+dzkm*v%wss+dzk*v%wssb)*indy(j-1)*incrzm
      wc = w%wsb
      fs = fs - 0.125D0*onethird*(wc+ABS(wc))*(vc+ABS(vc))*ABS(uc+ABS(uc))*deltar
!
! wn-wnb
      deltar  = dens%wn - dens%wnb
      IF (i <= 2) deltar = 0.D0
      IF (j >= ny-1) deltar = 0.D0
      uc = (dzkm*u%wn+dzk*u%wnb+dzkm*u%wwn+dzk*u%wwnb)*indx(i-1)*incrzm
      vc = (dzkm*v%wn+dzk*v%wnb+dzkm*v%w  +dzk*v%wb  )*indy(j+1)*incrzm
      wc = w%wnb
      fn = fn - 0.125D0*onethird*(wc+ABS(wc))*(vc-ABS(vc))*ABS(uc+ABS(uc))*deltar
!
! ent-et
      deltar  = dens%ent - dens%et
      IF (i >= nx-1) deltar = 0.D0
      IF (k >= nz-1) deltar = 0.D0
      uc = (dyjp*u%et+dyj*u%ent+dyjp*u%t+dyj*u%nt)*indx(i+1)*incryp
      vc = v%et
      wc = (dyjp*w%et+dyj*w%ent+dyjp*w%e+dyj*w%en)*indz(k+1)*incryp
      fe = fe - 0.125D0*onethird*(vc-ABS(vc))*(uc-ABS(uc))*ABS(wc-ABS(wc))*deltar
!
! et-est
      deltar  = dens%et - dens%est
      IF (i >= nx-1) deltar = 0.D0
      IF (k >= nz-1) deltar = 0.D0
      uc = (dyjm*u%et+dyj*u%est+dyjm*u%t+dyj*u%st)*indx(i+1)*incrym
      vc = v%est
      wc = (dyjm*w%et+dyj*w%est+dyjm*w%e+dyj*w%es)*indz(k+1)*incrym
      fe = fe - 0.125D0*onethird*(vc+ABS(vc))*(uc-ABS(uc))*ABS(wc-ABS(wc))*deltar
!
! enb-eb
      deltar  = dens%enb - dens%eb
      IF (i >= nx-1) deltar = 0.D0
      IF (k <= 2) deltar = 0.D0
      uc = (dyjp*u%eb+dyj*u%enb+dyjp*u%b  +dyj*u%nb  )*indx(i+1)*incryp
      vc = v%eb
      wc = (dyjp*w%eb+dyj*w%enb+dyjp*w%ebb+dyj*w%enbb)*indz(k-1)*incryp
      fe = fe - 0.125D0*onethird*(vc-ABS(vc))*(uc-ABS(uc))*ABS(wc+ABS(wc))*deltar
!
! eb-esb
      deltar  = dens%eb - dens%esb
      IF (i >= nx-1) deltar = 0.D0
      IF (k <= 2) deltar = 0.D0
      uc = (dyjm*u%eb+dyj*u%esb+dyjm*u%b  +dyj*u%sb  )*indx(i+1)*incrym
      vc = v%esb
      wc = (dyjm*w%eb+dyj*w%esb+dyjm*w%ebb+dyj*w%esbb)*indz(k-1)*incrym
      fe = fe - 0.125D0*onethird*(vc+ABS(vc))*(uc-ABS(uc))*ABS(wc+ABS(wc))*deltar
!
! est-es
      deltar  = dens%est - dens%es
      IF (i >= nx-1) deltar = 0.D0
      IF (j <= 2) deltar = 0.D0
      uc = (dzkp*u%es+dzk*u%est+dzkp*u%s  +dzk*u%st  )*indx(i+1)*incrzp
      vc = (dzkp*v%es+dzk*v%est+dzkp*v%ess+dzk*v%esst)*indy(j-1)*incrzp
      wc = w%es
      fs = fs - 0.125D0*onethird*(wc-ABS(wc))*(vc+ABS(vc))*ABS(uc-ABS(uc))*deltar
!
! es-esb
      deltar  = dens%es - dens%esb
      IF (i >= nx-1) deltar = 0.D0
      IF (j <= 2) deltar = 0.D0
      uc = (dzkm*u%es+dzk*u%esb+dzkm*u%s  +dzk*u%sb  )*indx(i+1)*incrzm
      vc = (dzkm*v%es+dzk*v%esb+dzkm*v%ess+dzk*v%essb)*indy(j-1)*incrzm
      wc = w%esb
      fs = fs - 0.125D0*onethird*(wc+ABS(wc))*(vc+ABS(vc))*ABS(uc-ABS(uc))*deltar
!
! ent-en
      deltar  = dens%ent - dens%en
      IF (i >= nx-1) deltar = 0.D0
      IF (j >= ny-1) deltar = 0.D0
      uc = (dzkp*u%en+dzk*u%ent+dzkp*u%n+dzk*u%nt)*indx(i+1)*incrzp
      vc = (dzkp*v%en+dzk*v%ent+dzkp*v%e+dzk*v%et)*indy(j+1)*incrzp
      wc = w%en
      fn = fn - 0.125D0*onethird*(wc-ABS(wc))*(vc-ABS(vc))*ABS(uc-ABS(uc))*deltar
!
! en-enb
      deltar  = dens%en - dens%enb
      IF (i >= nx-1) deltar = 0.D0
      IF (j >= ny-1) deltar = 0.D0
      uc = (dzkm*u%en+dzk*u%enb+dzkm*u%n+dzk*u%nb)*indx(i+1)*incrzm
      vc = (dzkm*v%en+dzk*v%enb+dzkm*v%e+dzk*v%eb)*indy(j+1)*incrzm
      wc = w%enb
      fn = fn - 0.125D0*onethird*(wc+ABS(wc))*(vc-ABS(vc))*ABS(uc-ABS(uc))*deltar
!
! ent-nt
      deltar  = dens%ent - dens%nt
      IF (k >= nz-1) deltar = 0.D0
      IF (j >= ny-1) deltar = 0.D0
      uc = u%nt
      vc = (dxip*v%nt+dxi*v%ent+dxip*v%t+dxi*v%et)*indy(j+1)*incrxp
      wc = (dxip*w%nt+dxi*w%ent+dxip*w%n+dxi*w%en)*indz(k+1)*incrxp
      ft = ft - 0.125D0*onethird*(wc-ABS(wc))*(uc-ABS(uc))*ABS(vc-ABS(vc))*deltar
!
! nt-wnt
      deltar  = dens%nt - dens%wnt
      IF (k >= nz-1) deltar = 0.D0
      IF (j >= ny-1) deltar = 0.D0
      uc = u%wnt
      vc = (dxim*v%nt+dxi*v%wnt+dxim*v%t+dxi*v%wt)*indy(j+1)*incrxm
      wc = (dxim*w%nt+dxi*w%wnt+dxim*w%n+dxi*w%wn)*indz(k+1)*incrxm
      ft = ft - 0.125D0*onethird*(wc-ABS(wc))*(uc+ABS(uc))*ABS(vc-ABS(vc))*deltar
!
! est-st
      deltar  = dens%est - dens%st
      IF (k >= nz-1) deltar = 0.D0
      IF (j <= 2) deltar = 0.D0
      uc = u%st
      vc = (dxip*v%st+dxi*v%est+dxip*v%sst+dxi*v%esst)*indy(j-1)*incrxp
      wc = (dxip*w%st+dxi*w%est+dxip*w%s  +dxi*w%es  )*indz(k+1)*incrxp
      ft = ft - 0.125D0*onethird*(wc-ABS(wc))*(uc-ABS(uc))*ABS(vc+ABS(vc))*deltar
!
! st-wst
      deltar  = dens%st - dens%wst
      IF (k >= nz-1) deltar = 0.D0
      IF (j <= 2) deltar = 0.D0
      uc = u%wst
      vc = (dxim*v%st+dxi*v%wst+dxim*v%sst+dxi*v%wsst)*indy(j-1)*incrxm
      wc = (dxim*w%st+dxi*w%wst+dxim*w%s  +dxi*w%ws  )*indz(k+1)*incrxm
      ft = ft - 0.125D0*onethird*(wc-ABS(wc))*(uc+ABS(uc))*ABS(vc+ABS(vc))*deltar
!
! esb-sb
      deltar  = dens%esb - dens%sb
      IF (j <= 2) deltar = 0.D0
      IF (k <= 2) deltar = 0.D0
      uc = u%sb
      vc = (dxip*v%sb+dxi*v%esb+dxip*v%ssb+dxi*v%essb)*indy(j-1)*incrxp
      wc = (dxip*w%sb+dxi*w%esb+dxip*w%sbb+dxi*w%esbb)*indz(k-1)*incrxp
      fb = fb - 0.125D0*onethird*(wc+ABS(wc))*(uc-ABS(uc))*ABS(vc+ABS(vc))*deltar
!
! sb-wsb
      deltar  = dens%sb - dens%wsb
      IF (j <= 2) deltar = 0.D0
      IF (k <= 2) deltar = 0.D0
      uc = u%wsb
      vc = (dxim*v%sb+dxi*v%wsb+dxim*v%ssb+dxi*v%wssb)*indy(j-1)*incrxm
      wc = (dxim*w%sb+dxi*w%wsb+dxim*w%sbb+dxi*w%wsbb)*indz(k-1)*incrxm
      fb = fb - 0.125D0*onethird*(wc+ABS(wc))*(uc+ABS(uc))*ABS(vc+ABS(vc))*deltar
!
! enb-nb
      deltar  = dens%enb - dens%nb
      IF (j >= ny-1) deltar = 0.D0
      IF (k <= 2) deltar = 0.D0
      uc = u%nb
      vc = (dxip*v%nb+dxi*v%enb+dxip*v%b  +dxi*v%eb  )*indy(j+1)*incrxp
      wc = (dxip*w%nb+dxi*w%enb+dxip*w%nbb+dxi*w%enbb)*indz(k-1)*incrxp
      fb = fb - 0.125D0*onethird*(wc+ABS(wc))*(uc-ABS(uc))*ABS(vc-ABS(vc))*deltar
!
! nb-wnb
      deltar  = dens%nb - dens%wnb
      IF (j >= ny-1) deltar = 0.D0
      IF (k <= 2) deltar = 0.D0
      uc = u%wnb
      vc = (dxim*v%nb+dxi*v%wnb+dxim*v%b  +dxi*v%wb  )*indy(j+1)*incrxm
      wc = (dxim*w%nb+dxi*w%wnb+dxim*w%nbb+dxi*w%wnbb)*indz(k-1)*incrxm
      fb = fb - 0.125D0*onethird*(wc+ABS(wc))*(uc+ABS(uc))*ABS(vc-ABS(vc))*deltar
!
      RETURN
      END SUBROUTINE ctu1_fmas_3d
!----------------------------------------------------------------------
      SUBROUTINE ctu2_fmas_3d(fe, fn, ft, fw, fs, fb, dens, u, v, w, ijk)
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

      REAL*8, INTENT(INOUT) :: fe, ft, fw, fb, fn, fs
      TYPE(stencil), INTENT(IN) :: dens, u, v, w 
      INTEGER, INTENT(IN) :: ijk 
!
      INTEGER :: im2, ip2, jm2, jp2, km2, kp2
      REAL*8 :: uref, vref, wref, deltar, erre, lim
      INTEGER :: i, j, k, imesh
!
      REAL*8 :: dxp, dxm, dxpp, dxmm, indxpp, indxp, indxm, indxmm
      REAL*8 :: dyp, dym, dypp, dymm, indypp, indyp, indym, indymm
      REAL*8 :: dzp, dzm, dzpp, dzmm, indzpp, indzp, indzm, indzmm
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
!
      indxm  = 2.D0/dxm
      indxp  = 2.D0/dxp
      indxmm = 2.D0/dxmm
      indxpp = 2.D0/dxpp
      indym  = 2.D0/dym
      indyp  = 2.D0/dyp
      indymm = 2.D0/dymm
      indypp = 2.D0/dypp
      indzm  = 2.D0/dzm
      indzp  = 2.D0/dzp
      indzmm = 2.D0/dzmm
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
      fe = fe + 0.5D0*uref*(1-dt*indxp*uref)*deltar*lim
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
      fw = fw + 0.5D0*uref*(1-dt*indxm*uref)*deltar*lim
!
! n-c
      gradc = indyp  * (dens%n  - dens%c)
      gradn = indypp * (dens%nn - dens%n)
      grads = indym  * (dens%c  - dens%s)
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
      deltar = dens%n - dens%c
      fn = fn + 0.5D0*vref*(1-dt*indyp*vref)*deltar*lim
!
! c-s
      gradc = indym  * (dens%c - dens%s )
      gradn = indyp  * (dens%n - dens%c )
      grads = indymm * (dens%s - dens%ss)
      lim  = 0.D0
      erre = 0.D0
      vref = v%s
      IF (vref >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = grads/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradn/gradc
      END IF
      IF (j/=2 .AND. flag(ijmk)==fluid) CALL limiters(lm,lim,erre)
      vref = ABS(v%s)
      deltar = dens%c - dens%s
      fs = fs + 0.5D0*vref*(1-dt*indym*vref)*deltar*lim
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
      END SUBROUTINE ctu2_fmas_3d
!----------------------------------------------------------------------
      SUBROUTINE ctu3_fmas_3d(fe, fn, ft, fw, fs, fb, dens, u, v, w, ijk)
!
! ... Compute the second order Corner Transport Upwind correction (step 3 in LeVeque algorithm)
!
      USE dimensions
      USE domain_mapping, ONLY: myijk, meshinds
      USE grid, ONLY: dx, dy, dz, indx, indy, indz, fluid, flag
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: stencil
      USE set_indexes, ONLY: imjk, ijkm, ipjk, ijkp, ijmk, ijpk
      USE set_indexes, ONLY: ipjpk, ipjmk, imjpk, imjmk, ipjkp, ipjkm, imjkp, imjkm
      USE set_indexes, ONLY: ijpkp, ijpkm, ijmkp, ijmkm, ipjpkp, imjmkm
      USE set_indexes, ONLY: ipjpkm, ipjmkp, imjpkp, ipjmkm, imjpkm, imjmkp
      USE time_parameters, ONLY: dt

      IMPLICIT NONE

      REAL*8, INTENT(INOUT) :: fe, ft, fw, fb, fn, fs
      TYPE(stencil), INTENT(IN) :: dens, u, v, w 
      INTEGER, INTENT(IN) :: ijk 

      INTEGER :: im2, ip2, jm2, jp2, km2, kp2
      REAL*8 :: uc, vc, wc, deltar, erre, lim, s
      INTEGER :: i, j, k, imesh

      REAL*8 :: dxp, dxm, dxpp, dxmm, indxpp, indxp, indxm, indxmm
      REAL*8 :: dyp, dym, dypp, dymm, indypp, indyp, indym, indymm
      REAL*8 :: dzp, dzm, dzpp, dzmm, indzpp, indzp, indzm, indzmm
      REAL*8 :: incrxp, incrxm, incryp, incrym, incrzp, incrzm
      REAL*8 :: gradc, grade, gradw, gradn, grads, gradt, gradb

      CALL meshinds(ijk, imesh, i, j, k)

      im2 = MAX(  1, i-2 )
      ip2 = MIN( nx, i+2 )
      jm2 = MAX(  1, j-2 )
      jp2 = MIN( ny, j+2 )
      km2 = MAX(  1, k-2 )
      kp2 = MIN( nz, k+2 )
!
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
!
      indxm  = 2.D0/dxm
      indxp  = 2.D0/dxp
      indxmm = 2.D0/dxmm
      indxpp = 2.D0/dxpp
      indym  = 2.D0/dym
      indyp  = 2.D0/dyp
      indymm = 2.D0/dymm
      indypp = 2.D0/dypp
      indzm  = 2.D0/dzm
      indzp  = 2.D0/dzp
      indzmm = 2.D0/dzmm
      indzpp = 2.D0/dzpp
!
      incrxm  = 0.25D0*dt*indxm
      incrxp  = 0.25D0*dt*indxp
      incrym  = 0.25D0*dt*indym
      incryp  = 0.25D0*dt*indyp
      incrzm  = 0.25D0*dt*indzm
      incrzp  = 0.25D0*dt*indzp
!
! c-w
      deltar  = dens%c - dens%w
      uc = u%w
      vc = (dx(i-1)*v%c+dx(i)*v%w+dx(i-1)*v%s+dx(i)*v%ws)*indy(j)*incrxm
      wc = (dx(i-1)*w%c+dx(i)*w%w+dx(i-1)*w%b+dx(i)*w%wb)*indz(k)*incrxm
      gradc = indxm  * (dens%c - dens%w )
      grade = indxp  * (dens%e - dens%c )
      gradw = indxmm * (dens%w - dens%ww)
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
      fs = fs + 0.5D0*(vc-ABS(vc))*s
      ft = ft + 0.5D0*(1-ABS(vc))*(wc+ABS(wc))*s
      fb = fb + 0.5D0*(1-ABS(vc))*(wc-ABS(wc))*s
!
! e-c
      deltar  = dens%e - dens%c
      uc = u%c
      vc = (dx(i+1)*v%c+dx(i)*v%e+dx(i+1)*v%s+dx(i)*v%es)*indy(j)*incrxp
      wc = (dx(i+1)*w%c+dx(i)*w%e+dx(i+1)*w%b+dx(i)*w%eb)*indz(k)*incrxp
      gradc = indxp  * (dens%e - dens%c )
      grade = indxpp * (dens%ee - dens%e)
      gradw = indxm  * (dens%c - dens%w )
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
      fs = fs - 0.5D0*(vc-ABS(vc))*s
      ft = ft - 0.5D0*(1-ABS(vc))*(wc+ABS(wc))*s
      fb = fb - 0.5D0*(1-ABS(vc))*(wc-ABS(wc))*s
!
! en-n
      deltar  = dens%en - dens%n
      uc = u%n
      vc = (dx(i+1)*v%n+dx(i)*v%en+dx(i+1)*v%c +dx(i)*v%e  )*indy(j+1)*incrxp
      wc = (dx(i+1)*w%n+dx(i)*w%en+dx(i+1)*w%nb+dx(i)*w%enb)*indz(k)  *incrxp
      gradc = indxp  * (dens%en  - dens%n )
      grade = indxpp * (dens%een - dens%en)
      gradw = indxm  * (dens%n   - dens%wn)
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
      fb = fb - 0.25D0*ABS(vc-ABS(vc))*(wc-ABS(wc))*s
!
! n-wn
      deltar  = dens%n - dens%wn
      uc = u%wn
      vc = (dx(i-1)*v%n+dx(i)*v%wn+dx(i-1)*v%c +dx(i)*v%w  )*indy(j+1)*incrxm
      wc = (dx(i-1)*w%n+dx(i)*w%wn+dx(i-1)*w%nb+dx(i)*w%wnb)*indz(k)  *incrxm
      gradc = indxm  * (dens%n  - dens%wn )
      grade = indxp  * (dens%en - dens%n  )
      gradw = indxmm * (dens%wn - dens%wwn)
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
      fb = fb + 0.25D0*ABS(vc-ABS(vc))*(wc-ABS(wc))*s
!
! es-s
      deltar  = dens%es - dens%s
      uc = u%s
      vc = (dx(i+1)*v%s+dx(i)*v%es+dx(i+1)*v%ss+dx(i)*v%ess)*indy(j-1)*incrxp
      wc = (dx(i+1)*w%s+dx(i)*w%es+dx(i+1)*w%sb+dx(i)*w%esb)*indz(k)  *incrxp
      gradc = indxp  * (dens%es  - dens%s )
      grade = indxpp * (dens%ees - dens%es)
      gradw = indxm  * (dens%s   - dens%ws)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=nx-1 .AND. flag(ipjmk)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxp*ABS(uc))*deltar*lim
      fs = fs - 0.5D0*(vc+ABS(vc))*s
      ft = ft - 0.25D0*ABS(vc+ABS(vc))*(wc+ABS(wc))*s
      fb = fb - 0.25D0*ABS(vc+ABS(vc))*(wc-ABS(wc))*s
!
! s-ws
      deltar  = dens%s - dens%ws
      uc = u%ws
      vc = (dx(i-1)*v%s+dx(i)*v%ws+dx(i-1)*v%ss+dx(i)*v%wss)*indy(j-1)*incrxm
      wc = (dx(i-1)*w%s+dx(i)*w%ws+dx(i-1)*w%sb+dx(i)*w%wsb)*indz(k)  *incrxm
      gradc = indxm  * (dens%s  - dens%ws )
      grade = indxp  * (dens%es - dens%s  )
      gradw = indxmm * (dens%ws - dens%wws)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=2 .AND. flag(imjmk)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxm*ABS(uc))*deltar*lim
      fs = fs + 0.5D0*(vc+ABS(vc))*s
      ft = ft + 0.25D0*ABS(vc+ABS(vc))*(wc+ABS(wc))*s
      fb = fb + 0.25D0*ABS(vc+ABS(vc))*(wc-ABS(wc))*s
!
! et-e
      deltar  = dens%et - dens%t
      uc = u%t
      vc = (dx(i+1)*v%t+dx(i)*v%et+dx(i+1)*v%st+dx(i)*v%est)*indy(j)  *incrxp
      wc = (dx(i+1)*w%t+dx(i)*w%et+dx(i+1)*w%c +dx(i)*w%e  )*indz(k+1)*incrxp
      gradc = indxp  * (dens%et  - dens%t )
      grade = indxpp * (dens%eet - dens%et)
      gradw = indxm  * (dens%t   - dens%wt)
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
      deltar  = dens%t - dens%wt
      uc = u%wt
      vc = (dx(i-1)*v%t+dx(i)*v%wt+dx(i-1)*v%st+dx(i)*v%wst)*indy(j)  *incrxm
      wc = (dx(i-1)*w%t+dx(i)*w%wt+dx(i-1)*w%c +dx(i)*w%w  )*indz(k+1)*incrxm
      gradc = indxm  * (dens%t  - dens%wt )
      grade = indxp  * (dens%et - dens%t  )
      gradw = indxmm * (dens%wt - dens%wwt)
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
! eb-b
      deltar  = dens%eb - dens%b
      IF (k <= 2) deltar = 0.D0
      uc = u%b
      vc = (dx(i+1)*v%b+dx(i)*v%eb+dx(i+1)*v%sb+dx(i)*v%esb)*indy(j)  *incrxp
      wc = (dx(i+1)*w%b+dx(i)*w%eb+dx(i+1)*w%bb+dx(i)*w%ebb)*indz(k-1)*incrxp
      gradc = indxp  * (dens%eb  - dens%b )
      grade = indxpp * (dens%eeb - dens%eb)
      gradw = indxm  * (dens%b   - dens%wb)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=nx-1 .AND. flag(ipjkm)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxp*ABS(uc))*deltar*lim
      fb = fb - 0.5D0*(1-ABS(vc))*(wc+ABS(wc))*s
!
! b-wb
      deltar  = dens%b - dens%wb
      IF (k <= 2) deltar = 0.D0
      uc = u%wb
      vc = (dx(i-1)*v%b+dx(i)*v%wb+dx(i-1)*v%sb+dx(i)*v%wsb)*indy(j)  *incrxm
      wc = (dx(i-1)*w%b+dx(i)*w%wb+dx(i-1)*w%bb+dx(i)*w%wbb)*indz(k-1)*incrxm
      gradc = indxm  * (dens%b  - dens%wb )
      grade = indxp  * (dens%eb - dens%b  )
      gradw = indxmm * (dens%wb - dens%wwb)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=2 .AND. flag(imjkm)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxm*ABS(uc))*deltar*lim
      fb = fb + 0.5D0*(1-ABS(vc))*(wc+ABS(wc))*s
!
! ent-nt
      deltar  = dens%ent - dens%nt
      uc = u%nt
      vc = (dx(i+1)*v%nt+dx(i)*v%ent+dx(i+1)*v%t+dx(i)*v%et)*indy(j+1)*incrxp
      wc = (dx(i+1)*w%nt+dx(i)*w%ent+dx(i+1)*w%n+dx(i)*w%en)*indz(k+1)*incrxp
      gradc = indxp  * (dens%ent  - dens%nt )
      grade = indxpp * (dens%eent - dens%ent)
      gradw = indxm  * (dens%nt   - dens%wnt)
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
      deltar  = dens%nt - dens%wnt
      uc = u%wnt
      vc = (dx(i-1)*v%nt+dx(i)*v%wnt+dx(i-1)*v%t+dx(i)*v%wt)*indy(j+1)*incrxm
      wc = (dx(i-1)*w%nt+dx(i)*w%wnt+dx(i-1)*w%n+dx(i)*w%wn)*indz(k+1)*incrxm
      gradc = indxm  * (dens%nt  - dens%wnt )
      grade = indxp  * (dens%ent - dens%nt  )
      gradw = indxmm * (dens%wnt - dens%wwnt)
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
      deltar  = dens%est - dens%st
      uc = u%st
      vc = (dx(i+1)*v%st+dx(i)*v%est+dx(i+1)*v%sst+dx(i)*v%esst)*indy(j-1)*incrxp
      wc = (dx(i+1)*w%st+dx(i)*w%est+dx(i+1)*w%s  +dx(i)*w%es  )*indz(k+1)*incrxp
      gradc = indxp   * (dens%est  - dens%st )
      grade = indxpp  * (dens%eest - dens%est)
      gradw = indxm   * (dens%st   - dens%wst)
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
      deltar  = dens%st - dens%wst
      uc = u%wst
      vc = (dx(i-1)*v%st+dx(i)*v%wst+dx(i-1)*v%sst+dx(i)*v%wsst)*indy(j-1)*incrxm
      wc = (dx(i-1)*w%st+dx(i)*w%wst+dx(i-1)*w%s  +dx(i)*w%ws  )*indz(k+1)*incrxm
      gradc = indxm  * (dens%st  - dens%wst )
      grade = indxp  * (dens%est - dens%st  )
      gradw = indxmm * (dens%wst - dens%wwst)
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
! esb-sb
      deltar  = dens%esb - dens%sb
      IF (k <= 2) deltar = 0.D0
      uc = u%sb
      vc = (dx(i+1)*v%sb+dx(i)*v%esb+dx(i+1)*v%ssb+dx(i)*v%essb)*indy(j-1)*incrxp
      wc = (dx(i+1)*w%sb+dx(i)*w%esb+dx(i+1)*w%sbb+dx(i)*w%esbb)*indz(k-1)*incrxp
      gradc = indxp  * (dens%esb  - dens%sb )
      grade = indxpp * (dens%eesb - dens%esb)
      gradw = indxm  * (dens%sb   - dens%wsb)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=nx-1 .AND. flag(ipjmkm)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxp*ABS(uc))*deltar*lim
      fb = fb - 0.25D0*ABS(vc+ABS(vc))*(wc+ABS(wc))*s
!
! sb-wsb
      deltar  = dens%sb - dens%wsb
      IF (k <= 2) deltar = 0.D0
      uc = u%wsb
      vc = (dx(i-1)*v%sb+dx(i)*v%wsb+dx(i-1)*v%ssb+dx(i)*v%wssb)*indy(j-1)*incrxm
      wc = (dx(i-1)*w%sb+dx(i)*w%wsb+dx(i-1)*w%sbb+dx(i)*w%wsbb)*indz(k-1)*incrxm
      gradc = indxm  * (dens%sb  - dens%wsb )
      grade = indxp  * (dens%esb - dens%sb  )
      gradw = indxmm * (dens%wsb - dens%wwsb)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=2 .AND. flag(imjmkm)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxm*ABS(uc))*deltar*lim
      fb = fb + 0.25D0*ABS(vc+ABS(vc))*(wc+ABS(wc))*s
!
! enb-nb
      deltar = dens%enb - dens%nb
      IF (k <= 2) deltar = 0.D0
      uc = u%nb
      vc = (dx(i+1)*v%nb+dx(i)*v%enb+dx(i+1)*v%b  +dx(i)*v%eb  )*indy(j+1)*incrxp
      wc = (dx(i+1)*w%nb+dx(i)*w%enb+dx(i+1)*w%nbb+dx(i)*w%enbb)*indz(k-1)*incrxp
      gradc = indxp  * (dens%enb  - dens%nb )
      grade = indxpp * (dens%eenb - dens%enb)
      gradw = indxm  * (dens%nb   - dens%wnb)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=nx-1 .AND. flag(ipjpkm)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxp*ABS(uc))*deltar*lim
      fb = fb - 0.25D0*ABS(vc-ABS(vc))*(wc+ABS(wc))*s
!
! nb-wnb
      deltar  = dens%nb - dens%wnb
      IF (k <= 2) deltar = 0.D0
      uc = u%wnb
      vc = (dx(i-1)*v%nb+dx(i)*v%wnb+dx(i-1)*v%b  +dx(i)*v%wb  )*indy(j+1)*incrxm
      wc = (dx(i-1)*w%nb+dx(i)*w%wnb+dx(i-1)*w%nbb+dx(i)*w%wnbb)*indz(k-1)*incrxm
      gradc = indxm  * (dens%nb - dens%wnb  )
      grade = indxp  * (dens%enb - dens%nb  )
      gradw = indxmm * (dens%wnb - dens%wwnb)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=2 .AND. flag(imjpkm)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxm*ABS(uc))*deltar*lim
      fb = fb + 0.25D0*ABS(vc-ABS(vc))*(wc+ABS(wc))*s
!
! c-s
      deltar  = dens%c - dens%s
      uc = (dy(j-1)*u%c+dy(j)*u%s+dy(j-1)*u%w+dy(j)*u%ws)*indx(i)*incrym
      vc = v%s
      wc = (dy(j-1)*w%c+dy(j)*w%s+dy(j-1)*w%b+dy(j)*w%sb)*indz(k)*incrym
      gradc = indym  * (dens%c - dens%s )
      grads = indymm * (dens%s - dens%ss)
      gradn = indyp  * (dens%n - dens%c )
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
      fb = fb + 0.5D0*(wc-ABS(wc))*s
      fe = fe + 0.5D0*(1-ABS(wc))*(uc+ABS(uc))*s
      fw = fw + 0.5D0*(1-ABS(wc))*(uc-ABS(uc))*s
!
! n-c
      deltar  = dens%n - dens%c
      uc = (dy(j+1)*u%c+dy(j)*u%n+dy(j+1)*u%w+dy(j)*u%wn)*indx(i)*incryp
      vc = v%c
      wc = (dy(j+1)*w%c+dy(j)*w%n+dy(j+1)*w%b+dy(j)*w%nb)*indz(k)*incryp
      gradc = indyp  * (dens%n - dens%c )
      grads = indym  * (dens%c - dens%s )
      gradn = indypp * (dens%nn - dens%n)
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
      fb = fb - 0.5D0*(wc-ABS(wc))*s
      fe = fe - 0.5D0*(1-ABS(wc))*(uc+ABS(uc))*s
      fw = fw - 0.5D0*(1-ABS(wc))*(uc-ABS(uc))*s
!
! wn-w
      deltar  = dens%wn - dens%w
      uc = (dy(j+1)*u%w+dy(j)*u%wn+dy(j+1)*u%ww+dy(j)*u%wwn)*indx(i-1)*incryp
      vc = v%w
      wc = (dy(j+1)*w%w+dy(j)*w%wn+dy(j+1)*w%wb+dy(j)*w%wnb)*indz(k)  *incryp
      gradc = indyp  * (dens%wn  - dens%w )
      grads = indym  * (dens%w   - dens%ws)
      gradn = indypp * (dens%wnn - dens%wn)
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = grads / gradc
      ELSE
	IF (gradc /= 0.D0) erre = gradn / gradc
      ENDIF
      IF (j/=ny-1 .AND. flag(imjpk)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indyp*ABS(vc))*deltar*lim
      fw = fw - 0.5D0*(1-ABS(wc))*(uc+ABS(uc))*s
!
! w-ws
      deltar  = dens%w - dens%ws
      uc = (dy(j-1)*u%w+dy(j)*u%ws+dy(j-1)*u%ww+dy(j)*u%wws)*indx(i-1)*incrym
      vc = v%ws
      wc = (dy(j-1)*w%w+dy(j)*w%ws+dy(j-1)*w%wb+dy(j)*w%wsb)*indz(k)  *incrym
      gradc = indym  * (dens%w  - dens%ws )
      grads = indymm * (dens%ws - dens%wss)
      gradn = indyp  * (dens%wn - dens%w  )
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = grads / gradc
      ELSE
	IF (gradc /= 0.D0) erre = gradn / gradc
      ENDIF
      IF (j/=2 .AND. flag(imjmk)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indym*ABS(vc))*deltar*lim
      fw = fw + 0.5D0*(1-ABS(wc))*(uc+ABS(uc))*s
!
! en-e
      deltar  = dens%en - dens%e
      uc = (dy(j+1)*u%e+dy(j)*u%en+dy(j+1)*u%c +dy(j)*u%n  )*indx(i+1)*incryp
      vc = v%e
      wc = (dy(j+1)*w%e+dy(j)*w%en+dy(j+1)*w%eb+dy(j)*w%enb)*indz(k)  *incryp
      gradc = indyp  * (dens%en  - dens%e )
      grads = indym  * (dens%e   - dens%es)
      gradn = indypp * (dens%enn - dens%en)
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
      deltar  = dens%e - dens%es
      uc = (dy(j-1)*u%e+dy(j)*u%es+dy(j-1)*u%c +dy(j)*u%s  )*indx(i+1)*incrym
      vc = v%es
      wc = (dy(j-1)*w%e+dy(j)*w%es+dy(j-1)*w%eb+dy(j)*w%esb)*indz(k)  *incrym
      gradc = indym  * (dens%e  - dens%es )
      grads = indymm * (dens%es - dens%ess)
      gradn = indyp  * (dens%en - dens%e  )
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
      deltar  = dens%nt - dens%t
      uc = (dy(j+1)*u%t+dy(j)*u%nt+dy(j+1)*u%wt+dy(j)*u%wnt)*indx(i)  *incryp
      vc = v%t
      wc = (dy(j+1)*w%t+dy(j)*w%nt+dy(j+1)*w%c +dy(j)*w%n  )*indz(k+1)*incryp
      gradc = indyp  * (dens%nt  - dens%t )
      grads = indym  * (dens%t   - dens%st)
      gradn = indypp * (dens%nnt - dens%nt)
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
      fw = fw - 0.25D0*ABS(wc-ABS(wc))*(uc-ABS(uc))*s
!
! t-st
      deltar  = dens%t - dens%st
      uc = (dy(j-1)*u%t+dy(j)*u%st+dy(j-1)*u%wt+dy(j)*u%wst)*indx(i)  *incrym
      vc = v%st
      wc = (dy(j-1)*w%t+dy(j)*w%st+dy(j-1)*w%c +dy(j)*w%s  )*indz(k+1)*incrym
      gradc = indym  * (dens%t  - dens%st )
      grads = indymm * (dens%st - dens%sst)
      gradn = indyp  * (dens%nt - dens%t  )
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
      fw = fw + 0.25D0*ABS(wc-ABS(wc))*(uc-ABS(uc))*s
!
! nb-b 
      deltar  = dens%nb - dens%b
      IF (k <= 2) deltar = 0.D0
      uc = (dy(j+1)*u%b+dy(j)*u%nb+dy(j+1)*u%wb+dy(j)*u%wnb)*indx(i)  *incryp
      vc = v%b
      wc = (dy(j+1)*w%b+dy(j)*w%nb+dy(j+1)*w%bb+dy(j)*w%nbb)*indz(k-1)*incryp
      gradc = indyp  * (dens%nb  - dens%b )
      grads = indym  * (dens%b   - dens%sb)
      gradn = indypp * (dens%nnb - dens%nb)
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = grads / gradc
      ELSE
	IF (gradc /= 0.D0) erre = gradn / gradc
      ENDIF
      IF (j/=ny-1 .AND. flag(ijpkm)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indyp*ABS(vc))*deltar*lim
      fb = fb - 0.5D0*(wc+ABS(wc))*s
      fe = fe - 0.25D0*ABS(wc+ABS(wc))*(uc+ABS(uc))*s
      fw = fw - 0.25D0*ABS(wc+ABS(wc))*(uc-ABS(uc))*s
!
! b-sb
      deltar  = dens%b - dens%sb
      IF (k <= 2) deltar = 0.D0
      uc = (dy(j-1)*u%b+dy(j)*u%sb+dy(j-1)*u%wb+dy(j)*u%wsb)*indx(i)  *incrym
      vc = v%sb
      wc = (dy(j-1)*w%b+dy(j)*w%sb+dy(j-1)*w%bb+dy(j)*w%sbb)*indz(k-1)*incrym
      gradc = indym  * (dens%b  - dens%sb )
      grads = indymm * (dens%sb - dens%ssb)
      gradn = indyp  * (dens%nb - dens%b  )
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = grads / gradc
      ELSE
	IF (gradc /= 0.D0) erre = gradn / gradc
      ENDIF
      IF (j/=2 .AND. flag(ijmkm)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indym*ABS(vc))*deltar*lim
      fb = fb + 0.5D0*(wc+ABS(wc))*s
      fe = fe + 0.25D0*ABS(wc+ABS(wc))*(uc+ABS(uc))*s
      fw = fw + 0.25D0*ABS(wc+ABS(wc))*(uc-ABS(uc))*s
!
! wt-wst
      deltar  = dens%wt - dens%wst
      uc = (dy(j-1)*u%wt+dy(j)*u%wst+dy(j-1)*u%wwt+dy(j)*u%wwst)*indx(i-1)*incrym
      vc = v%wst
      wc = (dy(j-1)*w%wt+dy(j)*w%wst+dy(j-1)*w%w  +dy(j)*w%ws  )*indz(k+1)*incrym
      gradc = indym  * (dens%wt  - dens%wst )
      grads = indymm * (dens%wst - dens%wsst)
      gradn = indyp  * (dens%wnt - dens%wt  )
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = grads / gradc
      ELSE
	IF (gradc /= 0.D0) erre = gradn / gradc
      ENDIF
      IF (j/=2 .AND. flag(imjmkp)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indym*ABS(vc))*deltar*lim
      fw = fw + 0.25D0*ABS(wc-ABS(wc))*(uc+ABS(uc))*s
!
! wnt-wt
      deltar  = dens%wnt - dens%wt
      uc = (dy(j+1)*u%wt+dy(j)*u%wnt+dy(j+1)*u%wwt+dy(j)*u%wwnt)*indx(i-1)*incryp
      vc = v%wt
      wc = (dy(j+1)*w%wt+dy(j)*w%wnt+dy(j+1)*w%w  +dy(j)*w%wn  )*indz(k+1)*incryp
      gradc = indyp  * (dens%wnt  - dens%wt )
      grads = indym  * (dens%wt   - dens%wst)
      gradn = indypp * (dens%wnnt - dens%wnt)
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = grads / gradc
      ELSE
	IF (gradc /= 0.D0) erre = gradn / gradc
      ENDIF
      IF (j/=ny-1 .AND. flag(imjpkp)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indyp*ABS(vc))*deltar*lim
      fw = fw - 0.25D0*ABS(wc-ABS(wc))*(uc+ABS(uc))*s
!
! wnb-wb
      deltar  = dens%wnb - dens%wb
      IF (k <= 2) deltar = 0.D0
      uc = (dy(j+1)*u%wb+dy(j)*u%wnb+dy(j+1)*u%wwb+dy(j)*u%wwnb)*indx(i-1)*incryp
      vc = v%wb
      wc = (dy(j+1)*w%wb+dy(j)*w%wnb+dy(j+1)*w%wbb+dy(j)*w%wnbb)*indz(k-1)*incryp
      gradc = indyp  * (dens%wnb  - dens%wb )
      grads = indym  * (dens%wb   - dens%wsb)
      gradn = indypp * (dens%wnnb - dens%wnb)
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = grads / gradc
      ELSE
	IF (gradc /= 0.D0) erre = gradn / gradc
      ENDIF
      IF (j/=ny-1 .AND. flag(imjpkm)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indyp*ABS(vc))*deltar*lim
      fw = fw - 0.25D0*ABS(wc+ABS(wc))*(uc+ABS(uc))*s
!
! wb-wsb
      deltar  = dens%wb - dens%wsb
      IF (k <= 2) deltar = 0.D0
      uc = (dy(j-1)*u%wb+dy(j)*u%wsb+dy(j-1)*u%wwb+dy(j)*u%wwsb)*indx(i-1)*incrym
      vc = v%wsb
      wc = (dy(j-1)*w%wb+dy(j)*w%wsb+dy(j-1)*w%wbb+dy(j)*w%wsbb)*indz(k-1)*incrym
      gradc = indym  * (dens%wb  - dens%wsb )
      grads = indymm * (dens%wsb - dens%wssb)
      gradn = indyp  * (dens%wnb - dens%wb  )
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = grads / gradc
      ELSE
	IF (gradc /= 0.D0) erre = gradn / gradc
      ENDIF
      IF (j/=2 .AND. flag(imjmkm)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indym*ABS(vc))*deltar*lim
      fw = fw + 0.25D0*ABS(wc+ABS(wc))*(uc+ABS(uc))*s
!
! ent-et
      deltar  = dens%ent - dens%et
      uc = (dy(j+1)*u%et+dy(j)*u%ent+dy(j+1)*u%t+dy(j)*u%nt)*indx(i+1)*incryp
      vc = v%et
      wc = (dy(j+1)*w%et+dy(j)*w%ent+dy(j+1)*w%e+dy(j)*w%en)*indz(k+1)*incryp
      gradc = indyp  * (dens%ent  - dens%et )
      grads = indym  * (dens%et   - dens%est)
      gradn = indypp * (dens%ennt - dens%ent)
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
      deltar  = dens%et - dens%est
      uc = (dy(j-1)*u%et+dy(j)*u%est+dy(j-1)*u%t+dy(j)*u%st)*indx(i+1)*incrym
      vc = v%est
      wc = (dy(j-1)*w%et+dy(j)*w%est+dy(j-1)*w%e+dy(j)*w%es)*indz(k+1)*incrym
      gradc = indym  * (dens%et  - dens%est )
      grads = indymm * (dens%est - dens%esst)
      gradn = indyp  * (dens%ent - dens%et  )
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
      deltar  = dens%enb - dens%eb
      IF (k <= 2) deltar = 0.D0
      uc = (dy(j+1)*u%eb+dy(j)*u%enb+dy(j+1)*u%b  +dy(j)*u%nb  )*indx(i+1)*incryp
      vc = v%eb
      wc = (dy(j+1)*w%eb+dy(j)*w%enb+dy(j+1)*w%ebb+dy(j)*w%enbb)*indz(k-1)*incryp
      gradc = indyp  * (dens%enb  - dens%eb )
      grads = indym  * (dens%eb   - dens%esb)
      gradn = indypp * (dens%ennb - dens%enb)
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
      deltar  = dens%eb - dens%esb
      IF (k <= 2) deltar = 0.D0
      uc = (dy(j-1)*u%eb+dy(j)*u%esb+dy(j-1)*u%b  +dy(j)*u%sb  )*indx(i+1)*incrym
      vc = v%esb
      wc = (dy(j-1)*w%eb+dy(j)*w%esb+dy(j-1)*w%ebb+dy(j)*w%esbb)*indz(k-1)*incrym
      gradc = indym  * (dens%eb  - dens%esb )
      grads = indymm * (dens%esb - dens%essb)
      gradn = indyp  * (dens%enb - dens%eb  )
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
      deltar  = dens%c - dens%b
      uc = (dz(k-1)*u%c+dz(k)*u%b+dz(k-1)*u%w+dz(k)*u%wb)*indx(i)*incrzm
      vc = (dz(k-1)*v%c+dz(k)*v%b+dz(k-1)*v%s+dz(k)*v%sb)*indy(j)*incrzm
      wc = w%b
      gradc = indzm  * (dens%c - dens%b )
      gradt = indzp  * (dens%t - dens%c )
      gradb = indzmm * (dens%b - dens%bb)
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
      fw = fe + 0.5D0*(uc-ABS(uc))*s
      fn = fn + 0.5D0*(1-ABS(uc))*(vc+ABS(vc))*s
      fs = fs + 0.5D0*(1-ABS(uc))*(vc-ABS(vc))*s
!
! t-c
      deltar  = dens%t - dens%c
      uc = (dz(k+1)*u%c+dz(k)*u%t+dz(k+1)*u%w+dz(k)*u%wt)*indx(i)*incrzp
      vc = (dz(k+1)*v%c+dz(k)*v%t+dz(k+1)*v%s+dz(k)*v%st)*indy(j)*incrzp
      wc = w%c
      gradc = indzp  * (dens%t  - dens%c)
      gradt = indzpp * (dens%tt - dens%t)
      gradb = indzm  * (dens%c  - dens%b)
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
      fw = fe - 0.5D0*(uc-ABS(uc))*s
      fn = fn - 0.5D0*(1-ABS(uc))*(vc+ABS(vc))*s
      fs = fs - 0.5D0*(1-ABS(uc))*(vc-ABS(vc))*s
!
! wt-w
      deltar  = dens%wt - dens%w
      uc = (dz(k+1)*u%w+dz(k)*u%wt+dz(k+1)*u%ww+dz(k)*u%wwt)*indx(i-1)*incrzp
      vc = (dz(k+1)*v%w+dz(k)*v%wt+dz(k+1)*v%ws+dz(k)*v%wst)*indy(j)  *incrzp
      wc = w%w
      gradc = indzp  * (dens%wt  - dens%w )
      gradt = indzpp * (dens%wtt - dens%wt)
      gradb = indzm  * (dens%w   - dens%wb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(imjkp)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indzp*ABS(wc))*deltar*lim
      fw = fw - 0.5D0*(uc+ABS(uc))*s
      fn = fn - 0.25D0*ABS(uc+ABS(uc))*(vc+ABS(vc))*s
      fs = fs - 0.25D0*ABS(uc+ABS(uc))*(vc-ABS(vc))*s
!
! w-wb
      deltar  = dens%w - dens%wb
      uc = (dz(k-1)*u%w+dz(k)*u%wb+dz(k-1)*u%ww+dz(k)*u%wwb)*indx(i-1)*incrzm
      vc = (dz(k-1)*v%w+dz(k)*v%wb+dz(k-1)*v%ws+dz(k)*v%wsb)*indy(j)  *incrzm
      wc = w%wb
      gradc = indzm  * (dens%w  - dens%wb )
      gradt = indzp  * (dens%wt - dens%w  )
      gradb = indzmm * (dens%wb - dens%wbb)
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
      fs = fs + 0.25D0*ABS(uc+ABS(uc))*(vc-ABS(vc))*s
!
! et-e
      deltar  = dens%et - dens%e
      uc = (dz(k+1)*u%e+dz(k)*u%et+dz(k+1)*u%c +dz(k)*u%t  )*indx(i+1)*incrzp
      vc = (dz(k+1)*v%e+dz(k)*v%et+dz(k+1)*v%es+dz(k)*v%est)*indy(j)  *incrzp
      wc = w%e
      gradc = indzp  * (dens%et  - dens%e )
      gradt = indzpp * (dens%ett - dens%et)
      gradb = indzm  * (dens%e   - dens%eb)
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
      fs = fs - 0.25D0*ABS(uc-ABS(uc))*(vc-ABS(vc))*s
!
! e-eb
      deltar  = dens%e - dens%eb
      uc = (dz(k-1)*u%e+dz(k)*u%eb+dz(k-1)*u%c +dz(k)*u%b  )*indx(i+1)*incrzm
      vc = (dz(k-1)*v%e+dz(k)*v%eb+dz(k-1)*v%es+dz(k)*v%esb)*indy(j)  *incrzm
      wc = w%eb
      gradc = indzm  * (dens%e  - dens%eb )
      gradt = indzp  * (dens%et - dens%e  )
      gradb = indzmm * (dens%eb - dens%ebb)
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
      fs = fs + 0.25D0*ABS(uc-ABS(uc))*(vc-ABS(vc))*s
!
! st-s
      deltar  = dens%st - dens%s
      uc = (dz(k+1)*u%s+dz(k)*u%st+dz(k+1)*u%ws+dz(k)*u%wst)*indx(i)  *incrzp
      vc = (dz(k+1)*v%s+dz(k)*v%st+dz(k+1)*v%ss+dz(k)*v%sst)*indy(j-1)*incrzp
      wc = w%s
      gradc = indzp  * (dens%st  - dens%s )
      gradt = indzpp * (dens%stt - dens%st)
      gradb = indzm  * (dens%s   - dens%sb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(ijmkp)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indzp*ABS(wc))*deltar*lim
      fs = fs - 0.5D0*(1-ABS(uc))*(vc+ABS(vc))*s
!
! s-sb
      deltar  = dens%s - dens%sb
      uc = (dz(k-1)*u%s+dz(k)*u%sb+dz(k-1)*u%ws+dz(k)*u%wsb)*indx(i)  *incrzm
      vc = (dz(k-1)*v%s+dz(k)*v%sb+dz(k-1)*v%ss+dz(k)*v%ssb)*indy(j-1)*incrzm
      wc = w%sb
      gradc = indzm  * (dens%s  - dens%sb )
      gradt = indzp  * (dens%st - dens%s  )
      gradb = indzmm * (dens%sb - dens%sbb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=2 .AND. flag(ijmkm)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indzm*ABS(wc))*deltar*lim
      fs = fs + 0.5D0*(1-ABS(uc))*(vc+ABS(vc))*s
!
! nt-n
      deltar  = dens%nt - dens%n
      uc = (dz(k+1)*u%n+dz(k)*u%nt+dz(k+1)*u%wn+dz(k)*u%wnt)*indx(i)  *incrzp
      vc = (dz(k+1)*v%n+dz(k)*v%nt+dz(k+1)*v%c +dz(k)*v%t  )*indy(j+1)*incrzp
      wc = w%n
      gradc = indzp  * (dens%nt  - dens%n )
      gradt = indzpp * (dens%ntt - dens%nt)
      gradb = indzm  * (dens%n   - dens%nb)
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
      deltar  = dens%n - dens%nb
      uc = (dz(k-1)*u%n+dz(k)*u%nb+dz(k-1)*u%wn+dz(k)*u%wnb)*indx(i)  *incrzm
      vc = (dz(k-1)*v%n+dz(k)*v%nb+dz(k-1)*v%c +dz(k)*v%b  )*indy(j+1)*incrzm
      wc = w%nb
      gradc = indzm  * (dens%n  - dens%nb )
      gradt = indzp  * (dens%nt - dens%n  )
      gradb = indzmm * (dens%nb - dens%nbb)
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
! wst-ws
      deltar  = dens%wst - dens%ws
      uc = (dz(k+1)*u%ws+dz(k)*u%wst+dz(k+1)*u%wws+dz(k)*u%wwst)*indx(i-1)*incrzp
      vc = (dz(k+1)*v%ws+dz(k)*v%wst+dz(k+1)*v%wss+dz(k)*v%wsst)*indy(j-1)*incrzp
      wc = w%ws
      gradc = indzp  * (dens%wst  - dens%ws )
      gradt = indzpp * (dens%wstt - dens%wst)
      gradb = indzm  * (dens%ws   - dens%wsb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(imjmkp)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indzp*ABS(wc))*deltar*lim
      fs = fs - 0.25D0*ABS(uc+ABS(uc))*(vc+ABS(uc))*s
!
! wnt-wn
      deltar  = dens%wnt - dens%wn
      uc = (dz(k+1)*u%wn+dz(k)*u%wnt+dz(k+1)*u%wwn+dz(k)*u%wwnt)*indx(i-1)*incrzp
      vc = (dz(k+1)*v%wn+dz(k)*v%wnt+dz(k+1)*v%w  +dz(k)*v%wt  )*indy(j+1)*incrzp
      wc = w%wn
      gradc = indzp  * (dens%wnt  - dens%wn )
      gradt = indzpp * (dens%wntt - dens%wnt)
      gradb = indzm  * (dens%wn   - dens%wnb)
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
! ws-wsb
      deltar  = dens%ws - dens%wsb
      uc = (dz(k-1)*u%ws+dz(k)*u%wsb+dz(k-1)*u%wws+dz(k)*u%wwsb)*indx(i-1)*incrzm
      vc = (dz(k-1)*v%ws+dz(k)*v%wsb+dz(k-1)*v%wss+dz(k)*v%wssb)*indy(j-1)*incrzm
      wc = w%wsb
      gradc = indzm  * (dens%ws  - dens%wsb )
      gradt = indzp  * (dens%wst - dens%ws  )
      gradb = indzmm * (dens%wsb - dens%wsbb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=2 .AND. flag(imjmkm)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indzm*ABS(wc))*deltar*lim
      fs = fs + 0.25D0*ABS(uc+ABS(uc))*(vc+ABS(vc))*s
!
! wn-wnb
      deltar  = dens%wn - dens%wnb
      uc = (dz(k-1)*u%wn+dz(k)*u%wnb+dz(k-1)*u%wwn+dz(k)*u%wwnb)*indx(i-1)*incrzm
      vc = (dz(k-1)*v%wn+dz(k)*v%wnb+dz(k-1)*v%w  +dz(k)*v%wb  )*indy(j+1)*incrzm
      wc = w%wnb
      gradc = indzm  * (dens%wn  - dens%wnb )
      gradt = indzp  * (dens%wnt - dens%wn  )
      gradb = indzmm * (dens%wnb - dens%wnbb)
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
! est-es
      deltar  = dens%est - dens%es
      uc = (dz(k+1)*u%es+dz(k)*u%est+dz(k+1)*u%s+dz(k)*u%st)*indx(i+1)*incrzp
      vc = (dz(k+1)*v%es+dz(k)*v%est+dz(k+1)*v%ess+dz(k)*v%esst)*indy(j-1)*incrzp
      wc = w%es
      gradc = indzp  * (dens%est  - dens%es )
      gradt = indzpp * (dens%estt - dens%est)
      gradb = indzm  * (dens%es   - dens%esb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(ipjmkp)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indzp*ABS(wc))*deltar*lim
      fs = fs - 0.25D0*ABS(uc-ABS(uc))*(vc+ABS(vc))*s
!
! es-esb
      deltar  = dens%es - dens%esb
      uc = (dz(k-1)*u%es+dz(k)*u%esb+dz(k-1)*u%s  +dz(k)*u%sb  )*indx(i+1)*incrzm
      vc = (dz(k-1)*v%es+dz(k)*v%esb+dz(k-1)*v%ess+dz(k)*v%essb)*indy(j-1)*incrzm
      wc = w%esb
      gradc = indzm  * (dens%es  - dens%esb )
      gradt = indzp  * (dens%est - dens%es  )
      gradb = indzmm * (dens%esb - dens%esbb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=2 .AND. flag(ipjmkm)==fluid) CALL limiters(lm,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indzm*ABS(wc))*deltar*lim
      fs = fs + 0.25D0*ABS(uc-ABS(uc))*(vc+ABS(vc))*s
!
! ent-en
      deltar  = dens%ent - dens%en
      uc = (dz(k+1)*u%en+dz(k)*u%ent+dz(k+1)*u%n+dz(k)*u%nt)*indx(i+1)*incrzp
      vc = (dz(k+1)*v%en+dz(k)*v%ent+dz(k+1)*v%e+dz(k)*v%et)*indy(j+1)*incrzp
      wc = w%en
      gradc = indzp  * (dens%ent  - dens%en )
      gradt = indzpp * (dens%entt - dens%ent)
      gradb = indzm  * (dens%en   - dens%enb)
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
      deltar  = dens%en - dens%enb
      uc = (dz(k-1)*u%en+dz(k)*u%enb+dz(k-1)*u%n+dz(k)*u%nb)*indx(i+1)*incrzm
      vc = (dz(k-1)*v%en+dz(k)*v%enb+dz(k-1)*v%e+dz(k)*v%eb)*indy(j+1)*incrzm
      wc = w%enb
      gradc = indzm  * (dens%en  - dens%enb )
      gradt = indzp  * (dens%ent - dens%en  )
      gradb = indzmm * (dens%enb - dens%enbb)
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

      END SUBROUTINE ctu3_fmas_3d
!-----------------------------------------------------------------------
      END MODULE convective_mass_fluxes
!-----------------------------------------------------------------------
