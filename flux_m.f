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
      USE flux_limiters, ONLY: muscl, limiters

      IMPLICIT NONE
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
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE masf_3d(fe, fn, ft, fw, fs, fb, dens, u, v, w, ijk)
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
! ... West, South and Bottom fluxes
!
      IF (u%w >= 0.D0) THEN
        fw = u%w * dens%w
      ELSE
        fw = u%w * dens%c
      ENDIF

      IF (v%s >= 0.D0) THEN
        fs = v%s * dens%s
      ELSE
        fb = v%s * dens%c
      ENDIF

      IF (w%b >= 0.D0) THEN
        fb = w%b * dens%b
      ELSE
        fb = w%b * dens%c
      ENDIF
!
! ... East, North and Top fluxes
!
      IF (u%c >= 0.D0) THEN
        fe = u%c * dens%c
      ELSE
        fe = u%c * dens%e
      ENDIF

      IF (v%c >= 0.D0) THEN
        fn = v%c * dens%c
      ELSE
        fn = v%c * dens%n
      ENDIF

      IF (w%c >= 0.D0) THEN
        ft = w%c * dens%c
      ELSE
        ft = w%c * dens%t
      ENDIF
!
      RETURN
      END SUBROUTINE masf_3d
!----------------------------------------------------------------------
      SUBROUTINE fmas_3d(fe, fn, ft, fw, fs, fb, dens, u, v, w, ijk)
!
      USE dimensions
      USE grid, ONLY: myijk, fl_l
      USE grid, ONLY: dx, dy, dz
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: stencil
      USE time_parameters, ONLY: dt
      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: fe, fn, ft, fw, fs, fb
      TYPE(stencil), INTENT(IN) :: dens, u, v, w
      INTEGER, INTENT(IN) :: ijk
!
      INTEGER :: i,j,k,imesh
      REAL*8 :: dxmm, dxm, dxp, dxpp, indxpp, indxp, indxm, indxmm
      REAL*8 :: dymm, dym, dyp, dypp, indypp, indyp, indym, indymm
      REAL*8 :: dzmm, dzm, dzp, dzpp, indzpp, indzp, indzm, indzmm
      REAL*8 :: gradc, grade, gradw, gradn, grads, gradt, gradb
!
      imesh = myijk( ip0_jp0_kp0_, ijk)
      i = MOD( MOD( imesh - 1, nx*ny ), nx ) + 1
      j = MOD( imesh - 1, nx*ny ) / nx + 1
      k = ( imesh - 1 ) / ( nx*ny ) + 1
!
      dxmm=dx(i-1)+dx(i-2)
      dxm=dx(i)+dx(i-1)
      dxp=dx(i)+dx(i+1)
      dxpp=dx(i+1)+dx(i+2)
      dymm=dy(j-1)+dy(j-2)
      dym=dy(j)+dy(j-1)
      dyp=dy(j)+dy(j+1)
      dypp=dy(j+1)+dy(j+2)
      dzmm=dz(k-1)+dz(k-2)
      dzm=dz(k)+dz(k-1)
      dzp=dz(k)+dz(k+1)
      dzpp=dz(k+1)+dz(k+2)

      indxmm=1.D0/dxmm
      indxm=1.D0/dxm
      indxp=1.D0/dxp
      indxpp=1.D0/dxpp
      indymm=1.D0/dymm
      indym=1.D0/dym
      indyp=1.D0/dyp
      indypp=1.D0/dypp
      indzmm=1.D0/dzmm
      indzm=1.D0/dzm
      indzp=1.D0/dzp
      indzpp=1.D0/dzpp
!
! ... MUSCL reconstruction of fields
!
! ... on West volume boundary
!
      gradc = 2.0 * indxm * (dens%c - dens%w)
      gradw = 2.0 * indxmm * (dens%w - dens%ww)
      grade = 2.0 * indxp * (dens%e - dens%c)
!
      lim = 0.D0
      erre = 0.D0
!
      cs = u%w
      cn = cs * dt * 2.0 * indxm
      IF (cs >= 0.D0) THEN
	erre = gradw / gradc
        fou  = dens%w 
	incr = 0.5D0 * dx(i-1)
      ELSE IF (cs < 0.D0) THEN
	erre = grade / gradc
        fou  = dens%c 
	incr = 0.5D0 * dx(i)
      ENDIF
!
      IF ((muscl /= 0) .AND. (gradc /= 0.D0)) THEN
        CALL limiters(lim,erre)
      END IF
!
      upwnd = fou + lim * gradc * incr
!
      centrd = (dx(i)*dens%w+dx(i-1)*dens%c)*indxm
      upc_w = upwnd / centrd
!
      fw = upwnd * cs
!
! ... on East volume boundary
!
      gradw = gradc
      gradc = grade
      grade = 2.0 * indxpp * (dens%ee - dens%e)
!
      lim = 0.D0
      erre = 0.D0
!
      cs = u%c
      cn = cs * dt * 2.0 * indxp
      IF (cs >= 0.D0) THEN
	erre = gradw / gradc
        fou  = dens%c 
	incr = 0.5D0 * dx(i)
      ELSE IF (cs < 0.D0) THEN
	erre = grade / gradc
        fou  = dens%e 
	incr = 0.5D0 * dx(i+1)
      ENDIF
!
      IF ((muscl /= 0) .AND. (gradc /= 0.D0) .AND. (i /= nx-1)) THEN
        CALL limiters(lim,erre)
      END IF
!
      upwnd = fou + lim * gradc * incr
!
      centrd = (dx(i)*dens%e+dx(i+1)*dens%c)*indxp
      upc_e = upwnd / centrd
!
      fe = upwnd * cs
!
! ... on South volume boundary
!
      gradc = 2.0 * indym * (dens%c - dens%s)
      grads = 2.0 * indymm * (dens%s - dens%ss)
      gradn = 2.0 * indyp * (dens%n - dens%c)
!
      lim = 0.D0
      erre = 0.D0
!
      cs = v%s
      cn = cs * dt * 2.0 * indym
      IF (cs >= 0.D0) THEN
	erre = grads / gradc
        fou  = dens%s 
	incr = 0.5D0 * dy(j-1)
      ELSE IF (cs < 0.D0) THEN
	erre = gradn / gradc
        fou  = dens%c 
	incr = 0.5D0 * dy(j)
      ENDIF
!
      IF ((muscl /= 0) .AND. (gradc /= 0.D0)) THEN
        CALL limiters(lim,erre)
      END IF
!
      upwnd = fou + lim * gradc * incr
!
      centrd = (dy(j)*dens%s+dy(j-1)*dens%c)*indym
      upc_s = upwnd / centrd
!
      fs = upwnd * cs
!
! ... on North volume boundary
!
      grads = gradc
      gradc = gradn
      gradn = 2.0 * indypp * (dens%nn - dens%n)
!
      lim = 0.D0
      erre = 0.D0
!
      cs = v%c
      cn = cs * dt * 2.0 * indyp
      IF (cs >= 0.D0) THEN
	erre = grads / gradc
        fou  = dens%c 
	incr = 0.5D0 * dy(j)
      ELSE IF (cs < 0.D0) THEN
	erre = gradn / gradc
        fou  = dens%n 
	incr = 0.5D0 * dy(j+1)
      ENDIF
!
      IF ((muscl /= 0) .AND. (gradc /= 0.D0) .AND. (j /= ny-1)) THEN
        CALL limiters(lim,erre)
      END IF
!
      upwnd = fou + lim * gradc * incr
!
      centrd = (dy(j)*dens%n+dy(j+1)*dens%c)*indyp
      upc_n = upwnd / centrd
!
      fn = upwnd * cs
!
! ... on Bottom volume boundary
!
      gradc = 2.0 * indzm * (dens%c - dens%b)
      gradb = 2.0 * indzmm * (dens%b - dens%bb)
      gradt = 2.0 * indzp * (dens%t - dens%c)
!
      lim = 0.D0
      erre = 0.D0
!
      cs = w%b
      cn = cs * dt * 2.0 * indzm
      IF (cs >= 0.D0) THEN
	erre = gradb / gradc
        fou  = dens%b 
	incr = 0.5D0 * dz(k-1)
      ELSE IF (cs < 0.D0) THEN
	erre = gradt / gradc
        fou  = dens%c 
	incr = 0.5D0 * dz(k)
      ENDIF
!
      IF ((muscl /= 0) .AND. (gradc /= 0.D0)) THEN
        CALL limiters(lim,erre)
      END IF
!
      upwnd = fou + lim * gradc * incr
!
      centrd = (dz(k)*dens%b+dz(k-1)*dens%c)*indzm
      upc_b = upwnd / centrd
!
      fb = upwnd * cs
!
! ... on Top volume boundary
!
      gradb = gradc
      gradc = gradt
      gradt = 2.0 * indzpp * (dens%tt - dens%t)
!
      lim = 0.D0
      erre = 0.D0
!
      cs = w%c
      cn = cs * dt * 2.0 * indzp
      IF (cs >= 0.D0) THEN
	erre = gradb / gradc
        fou  = dens%c 
	incr = 0.5D0 * dz(k)
      ELSE IF (cs < 0.D0) THEN
	erre = gradt / gradc
        fou  = dens%t 
	incr = 0.5D0 * dz(k+1)
      ENDIF
!
      IF ((muscl /= 0) .AND. (gradc /= 0.D0) .AND. (k /= nz-1)) THEN
        CALL limiters(lim,erre)
      END IF
!
      upwnd = fou + lim * gradc * incr
!
      centrd = (dz(k)*dens%t+dz(k+1)*dens%c)*indzp
      upc_t = upwnd / centrd
!
      ft = upwnd * cs
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
      USE grid, ONLY: myijk
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
      j = ( imesh - 1 ) / nr + 1
      i = MOD( ( imesh - 1 ), nr) + 1
!
!
! ... West, Bottom fluxes
!
      IF (u%w >= 0.D0) THEN
        fw = u%w * dens%w * rb(i-1)
      ELSE
        fw = u%w * dens%c * rb(i-1)
      ENDIF

      IF (w%b >= 0.D0) THEN
        fb = w%b * dens%b
      ELSE
        fb = w%b * dens%c
      ENDIF
!
! ... East, Top fluxes
!
      IF (u%c >= 0.D0) THEN
        fe = u%c * dens%c * rb(i)
      ELSE
        fe = u%c * dens%e * rb(i)
      ENDIF

      IF (w%c >= 0.D0) THEN
        ft = w%c * dens%c
      ELSE
        ft = w%c * dens%t
      ENDIF
!
      RETURN
      END SUBROUTINE masf_2d
!----------------------------------------------------------------------
      SUBROUTINE fmas_2d(fe, ft, fw, fb, dens, u, w, ij)
!
      USE dimensions
      USE grid, ONLY: myijk, fl_l
      USE grid, ONLY: dr, rb, dz
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: stencil
      USE time_parameters, ONLY: dt
      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: fe, ft, fw, fb
      TYPE(stencil), INTENT(IN) :: dens, u, w
      INTEGER, INTENT(IN) :: ij
!
      INTEGER :: i,j,imesh
      REAL*8 :: drmm, drm, drp, drpp, indrpp, indrp, indrm, indrmm
      REAL*8 :: dzmm, dzm, dzp, dzpp, indzpp, indzp, indzm, indzmm
      REAL*8 :: gradc, grade, gradw, gradt, gradb
!
      imesh = myijk( ip0_jp0_kp0_, ij)
      j = ( imesh - 1 ) / nr + 1
      i = MOD( ( imesh - 1 ), nr) + 1
!
      drmm=dr(i-1)+dr(i-2)
      drm=dr(i)+dr(i-1)
      drp=dr(i)+dr(i+1)
      drpp=dr(i+1)+dr(i+2)
      dzmm=dz(j-1)+dz(j-2)
      dzm=dz(j)+dz(j-1)
      dzp=dz(j)+dz(j+1)
      dzpp=dz(j+1)+dz(j+2)

      indrmm=1.D0/drmm
      indrm=1.D0/drm
      indrp=1.D0/drp
      indrpp=1.D0/drpp
      indzmm=1.D0/dzmm
      indzm=1.D0/dzm
      indzp=1.D0/dzp
      indzpp=1.D0/dzpp
!
! ... MUSCL reconstruction of fields
!
! ... on West volume boundary
!
      gradc = 2.D0 * indrm * (dens%c - dens%w)
      gradw = 2.D0 * indrmm * (dens%w - dens%ww)
      grade = 2.D0 * indrp * (dens%e - dens%c)
!
      lim = 0.D0
      erre = 0.D0
!
      cs = u%w
      cn = cs * dt * 2.D0 * indrm
      IF (cs >= 0.D0) THEN
	erre = gradw / gradc
        fou  = dens%w 
	incr = 0.5D0 * dr(i-1)
      ELSE IF (cs < 0.D0) THEN
	erre = grade / gradc
        fou  = dens%c 
	incr = 0.5D0 * dr(i)
      ENDIF
!
      IF ((muscl /= 0) .AND. (gradc /= 0.D0)) THEN
        CALL limiters(lim,erre)
      END IF
!
      upwnd = fou + lim * gradc * incr
!
      centrd = (dr(i)*dens%w+dr(i-1)*dens%c)*indrm
      upc_w = upwnd / centrd
!
      fw = upwnd * cs * rb(i-1)
!
! ... on East volume boundary
!
      gradw = gradc
      gradc = grade
      grade = 2.D0 * indrpp * (dens%ee - dens%e)
!
      lim = 0.D0
      erre = 0.D0
!
      cs = u%c
      cn = cs * dt * 2.0 * indrp
      IF (cs >= 0.D0) THEN
	erre = gradw / gradc
        fou  = dens%c 
	incr = 0.5D0 * dr(i)
      ELSE IF (cs < 0.D0) THEN
	erre = grade / gradc
        fou  = dens%e 
	incr = 0.5D0 * dr(i+1)
      ENDIF
!
      IF ((muscl /= 0) .AND. (gradc /= 0.D0) .AND. (i /= nr-1)) THEN
        CALL limiters(lim,erre)
      END IF
!
      upwnd = fou + lim * gradc * incr
!
      centrd = (dr(i)*dens%e+dr(i+1)*dens%c)*indrp
      upc_e = upwnd / centrd
!
      fe = upwnd * cs * rb(i)
!
! ... on Bottom volume boundary
!
      gradc = 2.D0 * indzm * (dens%c - dens%b)
      gradb = 2.D0 * indzmm * (dens%b - dens%bb)
      gradt = 2.D0 * indzp * (dens%t - dens%c)
!
      lim = 0.D0
      erre = 0.D0
!
      cs = w%b
      cn = cs * dt * 2.D0 * indzm
      IF (cs >= 0.D0) THEN
	erre = gradb / gradc
        fou  = dens%b 
	incr = 0.5D0 * dz(j-1)
      ELSE IF (cs < 0.D0) THEN
	erre = gradt / gradc
        fou  = dens%c 
	incr = 0.5D0 * dz(j)
      ENDIF
!
      IF ((muscl /= 0) .AND. (gradc /= 0.D0)) THEN
        CALL limiters(lim,erre)
      END IF
!
      upwnd = fou + lim * gradc * incr
!
      centrd = (dz(j)*dens%b+dz(j-1)*dens%c)*indzm
      upc_s = upwnd / centrd
!
      fb = upwnd * cs
!
! ... on Top volume boundary
!
      gradt = gradc
      gradc = gradt
      gradt = 2.D0 * indzpp * (dens%tt - dens%t)
!
      lim = 0.D0
      erre = 0.D0
!
      cs = w%c
      cn = cs * dt * 2.D0 * indzp
      IF (cs >= 0.D0) THEN
	erre = gradb / gradc
        fou  = dens%c 
	incr = 0.5D0 * dz(j)
      ELSE IF (cs < 0.D0) THEN
	erre = gradt / gradc
        fou  = dens%t 
	incr = 0.5D0 * dz(j+1)
      ENDIF
!
      IF ((muscl /= 0) .AND. (gradc /= 0.D0) .AND. (j /= nz-1)) THEN
        CALL limiters(lim,erre)
      END IF
!
      upwnd = fou + lim * gradc * incr
!
      centrd = (dz(j)*dens%t+dz(j+1)*dens%c)*indzp
      upc_n = upwnd / centrd
!
      ft = upwnd * cs
!
      RETURN
      END SUBROUTINE fmas_2d
!-----------------------------------------------------------------------
      END MODULE convective_mass_fluxes
!-----------------------------------------------------------------------
