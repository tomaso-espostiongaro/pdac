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
!
      USE flux_limiters, ONLY: muscl, limiters

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
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE flu_3d(fe, fn, ft, fw, fs, fb, dens, u, v, w, ijk)
!
! ... Compute the convective fluxes on East, North, and Top sides of the cell
! ... for the momentum density along x.
!
      USE dimensions
      USE grid, ONLY: myijk, fl_l
      USE grid, ONLY: dx, dy, dz, indx
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: imjk, ijmk, ijkm
      USE set_indexes, ONLY: stencil
      USE time_parameters, ONLY: dt

      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: fe, fn, ft, fw, fs, fb
      TYPE(stencil), INTENT(IN) :: dens, u, v, w
      INTEGER, INTENT(IN) :: ijk
      INTEGER :: i,j,k,imesh
!
      REAL*8 :: dens_c, dens_e, dens_n, dens_t
      REAL*8 :: dens_w, dens_s, dens_b
      REAL*8 :: dens_ee, dens_nn, dens_tt
      REAL*8 :: dxm, dxp, dxpp, indxpp, indxp, indxm
      REAL*8 :: dym, dyp, dypp, indypp, indyp, indym
      REAL*8 :: dzp, indzp, dzm, indzm, dzpp, indzpp
      REAL*8 :: gradc, grade, gradw, gradn, grads, gradt, gradb
!
      imesh = myijk( ip0_jp0_kp0_, ijk)
      i = MOD( MOD( imesh - 1, nx*ny ), nx ) + 1
      j = MOD( imesh - 1, nx*ny ) / nx + 1
      k = ( imesh - 1 ) / ( nx*ny ) + 1
!
      dxm=dx(i)+dx(i-1)
      dxp=dx(i)+dx(i+1)
      dxpp=dx(i+1)+dx(i+2)
      dym=dy(j)+dy(j-1)
      dyp=dy(j)+dy(j+1)
      dypp=dy(j+1)+dy(j+2)
      dzm=dz(k)+dz(k-1)
      dzp=dz(k)+dz(k+1)
      dzpp=dz(k+1)+dz(k+2)

      indxm=1.D0/dxm
      indxp=1.D0/dxp
      indxpp=1.D0/dxpp
      indym=1.D0/dym
      indyp=1.D0/dyp
      indypp=1.D0/dypp
      indzm=1.D0/dzm
      indzp=1.D0/dzp
      indzpp=1.D0/dzpp
!       
! ... Compute linearly interpolated values of density on the staggered grid
!
      dens_c = (dx(i+1) * dens%c + dx(i) * dens%e) * indxp
      dens_e = (dx(i+2) * dens%e + dx(i+1) * dens%ee) * indxpp
      dens_n = (dx(i+1) * dens%n + dx(i) * dens%en) * indxp
      dens_t = (dx(i+1) * dens%t + dx(i) * dens%et) * indxp
      dens_w = (dx(i)   * dens%w + dx(i-1) * dens%c) * indxm
      dens_s = (dx(i+1) * dens%s + dx(i) * dens%es) * indxp
      dens_b = (dx(i+1) * dens%b + dx(i) * dens%eb) * indxp
!
! ... an arbitrary choice !
!
      dens_ee = dens_e
      dens_nn = dens_n
      dens_tt = dens_t
!
! ... On boundary mantain first order accuracy (1st order Upwind).
!
! ... on West volume bondary
!
      IF( fl_l(imjk) /= 1 ) THEN
        cs = 0.5D0*(u%c + u%w)
        IF ( cs >= 0.D0 ) fw = dens_w * u%w * cs
        IF ( cs <  0.D0 ) fw = dens_c * u%c * cs
      END IF
!
! ... on South volume bondary
!
      IF( fl_l(ijmk) /= 1 ) THEN
        cs = (dx(i+1) * v%s + dx(i) * v%es) * indxp
        IF ( cs >= 0.D0 ) fs = dens_s * u%s * cs
        IF ( cs <  0.D0 ) fs = dens_c * u%c * cs
      END IF
!
! ... on Bottom volume bondary
!
      IF( fl_l(ijkm) /= 1 ) THEN
        cs = (dx(i+1) * w%b + dx(i) * w%eb) * indxp
        IF ( cs >= 0.D0 ) fb = dens_b * u%b * cs
        IF ( cs <  0.D0 ) fb = dens_c * u%c * cs
      END IF
!
! ... MUSCL reconstruction of momentum
!
! ... on East volume boundary
!
      gradc = (indx(i+1) * (dens_e * u%e   - dens_c * u%c))
      gradw = (indx(i)   * (dens_c * u%c   - dens_w * u%w))
      grade = (indx(i+2) * (dens_ee * u%ee - dens_e * u%e))
!
      lim = 0.D0
      erre = 0.D0
!
      cs = 0.5D0 * (u%c + u%e)
      cn = cs * dt * indx(i+1)
      IF ( cs >= 0.D0 ) THEN
        erre = gradw / gradc
        fou  = dens_c * u%c 
        incr = 0.5D0 * dx(i+1)
      ELSE IF ( cs < 0.D0 ) THEN
        erre = grade / gradc
        fou  = dens_e * u%e
        incr = 0.5D0 * dx(i+1)
      END IF
!
      IF ((muscl /= 0) .AND. (gradc /= 0.D0) .AND. (i /= nx-1)) THEN
        CALL limiters(lim,erre)
      END IF
!
      upwnd = fou + lim * gradc * incr
!
      fe = upwnd * cs
!
! ... on North volume boundary
!
      gradc = (dens_n * u%n   - dens_c * u%c) * 2.0 * indyp
      grads = (dens_c * u%c   - dens_s * u%s) * 2.0 * indym
      gradn = (dens_nn * u%nn - dens_n * u%n) * 2.0 * indypp
!
      lim = 0.D0
      erre = 0.D0
!
      cs = (dx(i+1) * v%c + dx(i) * v%e) * indxp
      cn = cs * dt * 2.0 * indyp
      IF (cs >= 0.D0) THEN
        erre = grads / gradc 
        fou  = dens_c * u%c
	incr = 0.5D0 * dy(j)
      ELSE IF (cs < 0.D0) THEN
        erre = gradn / gradc 
        fou  = dens_n * u%n
	incr = 0.5D0 * dy(j+1)
      END IF 
!
      IF ((muscl /= 0) .AND. (gradc /= 0.D0) .AND. (j /= ny-1)) THEN
        CALL limiters(lim,erre)
      END IF
!
      upwnd = fou + lim * gradc * incr
!
      fn = upwnd * cs
!
! ... on Top volume boundary
!
      gradc = (dens_t * u%t   - dens_c * u%c) * 2.0 * indzp
      gradb = (dens_c * u%c   - dens_b * u%b) * 2.0 * indzm
      gradt = (dens_tt * u%tt - dens_t * u%t) * 2.0 * indzpp
!
      lim = 0.D0
      erre = 0.D0
!
      cs = (dx(i+1) * w%c + dx(i) * w%e) * indxp
      cn = cs * dt * 2.0 * indzp
      IF (cs >= 0.D0) THEN
	erre = gradb / gradc
        fou  = dens_c * u%c
	incr = 0.5D0 * dz(k)
      ELSE IF (cs < 0.D0) THEN
	erre = gradt / gradc
        fou  = dens_t * u%t
	incr = 0.5D0 * dz(k+1)
      END IF 
!
      IF ((muscl /= 0) .AND. (gradc /= 0.D0) .AND. (k /= nz-1)) THEN
        CALL limiters(lim,erre)
      END IF
!
      upwnd = fou + lim * gradc * incr
!
      ft = upwnd * cs
!
      RETURN
      END SUBROUTINE flu_3d
!----------------------------------------------------------------------
      SUBROUTINE flu_2d(fe, fn, fw, fs, dens, u, w, ij)
!
! ... Compute the convective fluxes on East, North sides of the cell
! ... for the momentum density along r.
!
      USE dimensions
      USE grid, ONLY: myijk, fl_l
      USE grid, ONLY: dz
      USE grid, ONLY: dr, indr, r
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: imj, ijm
      USE set_indexes, ONLY: stencil
      USE time_parameters, ONLY: dt

      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: fe, fn, fw, fs
      TYPE(stencil), INTENT(IN) :: dens, u, w
      INTEGER, INTENT(IN) :: ij
      INTEGER :: i,j,imesh
!
      REAL*8 :: dens_c, dens_e, dens_n
      REAL*8 :: dens_w, dens_s
      REAL*8 :: dens_ee, dens_nn
      REAL*8 :: drm, drp, drpp, indrpp, indrp, indrm
      REAL*8 :: dzp, indzp, dzm, indzm, dzpp, indzpp
      REAL*8 :: gradc, grade, gradw, grads, gradn
!
      imesh = myijk( ip0_jp0_kp0_, ij)
      j = ( imesh - 1 ) / nr + 1
      i = MOD( ( imesh - 1 ), nr) + 1
!
      drm=dr(i)+dr(i-1)
      drp=dr(i)+dr(i+1)
      drpp=dr(i+1)+dr(i+2)
      dzm=dz(j)+dz(j-1)
      dzp=dz(j)+dz(j+1)
      dzpp=dz(j+1)+dz(j+2)

      indrm=1.D0/drm
      indrp=1.D0/drp
      indrpp=1.D0/drpp
      indzm=1.D0/dzm
      indzp=1.D0/dzp
      indzpp=1.D0/dzpp
!       
! ... Compute linearly interpolated values of density on the staggered grid
!
      dens_c = (dr(i+1) * dens%c + dr(i) * dens%e) * indrp
      dens_e = (dr(i+2) * dens%e + dr(i+1) * dens%ee) * indrpp
      dens_n = (dr(i+1) * dens%n + dr(i) * dens%en) * indrp
      dens_w = (dr(i)   * dens%w + dr(i-1) * dens%c) * indrm
      dens_s = (dr(i+1) * dens%s + dr(i) * dens%es) * indrp
!
! ... an arbitrary choice !
!
      dens_ee = dens_e
      dens_nn = dens_n
!
! ... On boundary mantain first order accuracy (1st order Upwind).
!
! ... on West volume bondary
!
      IF( fl_l(imj) /= 1 ) THEN
        cs = 0.5D0*(u%c + u%w)
        IF ( cs >= 0.D0 ) fw = dens_w * u%w * cs * r(i)
        IF ( cs <  0.D0 ) fw = dens_c * u%c * cs * r(i)
      END IF
!
! ... on South volume bondary
!
      IF( fl_l(ijm) /= 1 ) THEN
        cs = (dr(i+1) * w%s + dr(i) * w%es) * indrp
        IF ( cs >= 0.D0 ) fs = dens_s * u%s * cs
        IF ( cs <  0.D0 ) fs = dens_c * u%c * cs
      END IF
!
! ... MUSCL reconstruction of momentum
!
! ... on East volume boundary
!
      gradc = (indr(i+1) * (dens_e * u%e   - dens_c * u%c))
      gradw = (indr(i)   * (dens_c * u%c   - dens_w * u%w))
      grade = (indr(i+2) * (dens_ee * u%ee - dens_e * u%e))
!
      lim = 0.D0
      erre = 0.D0
!
      cs = 0.5D0 * (u%c + u%e)
      cn = cs * dt * indr(i+1)
      IF ( cs >= 0.D0 ) THEN
        erre = gradw / gradc
        fou  = dens_c * u%c 
        incr = 0.5D0 * dr(i+1)
      ELSE IF ( cs < 0.D0 ) THEN
        erre = grade / gradc
        fou  = dens_e * u%e
        incr = 0.5D0 * dr(i+1)
      END IF
!
      IF ((muscl /= 0) .AND. (gradc /= 0.D0) .AND. (i /= nr-1)) THEN
        CALL limiters(lim,erre)
      END IF
!
      upwnd = fou + lim * gradc * incr
!
      fe = upwnd * cs * r(i+1)
!
! ... on North volume boundary
!
      gradc = (dens_n * u%n   - dens_c * u%c) * 2.D0 * indzp
      grads = (dens_c * u%c   - dens_s * u%s) * 2.D0 * indzm
      gradn = (dens_nn * u%nn - dens_n * u%n) * 2.D0 * indzpp
!
      lim = 0.D0
      erre = 0.D0
!
      cs = (dr(i+1) * w%c + dr(i) * w%e) * indrp
      cn = cs * dt * 2.0 * indzp
      IF (cs >= 0.D0) THEN
        erre = grads / gradc 
        fou  = dens_c * u%c
	incr = 0.5D0 * dz(j)
      ELSE IF (cs < 0.D0) THEN
        erre = gradn / gradc 
        fou  = dens_n * u%n
	incr = 0.5D0 * dz(j+1)
      END IF 
!
      IF ((muscl /= 0) .AND. (gradc /= 0.D0) .AND. (j /= nz-1)) THEN
        CALL limiters(lim,erre)
      END IF
!
      upwnd = fou + lim * gradc * incr
!
      fn = upwnd * cs
!
      RETURN
      END SUBROUTINE flu_2d
!-----------------------------------------------------------------------
      END MODULE convective_fluxes_u
!-----------------------------------------------------------------------
