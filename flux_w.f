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
      INTERFACE flw
        MODULE PROCEDURE flw_2d, flw_3d
      END INTERFACE
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE flw_3d(fe, fn, ft, fw, fs, fb, dens, u, v, w, ijk)
!
! ... Compute the convective fluxes on East, North, and Top sides of the cell
! ... for the momentum density along z.
!
      USE dimensions
      USE domain_decomposition, ONLY: myijk
      USE grid, ONLY: dx, dy, dz, indz, fl_l
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: imjk, ijmk, ijkm
      USE set_indexes, ONLY: stencil
      USE time_parameters, ONLY: dt
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: ijk
      REAL*8, INTENT(OUT) :: fe, fn, ft, fw, fs, fb
      TYPE(stencil), INTENT(IN) :: dens, u, v, w
!
      INTEGER :: i,j,k, imesh
      REAL*8 :: dxm, dxp, dxpp, indxpp, indxp, indxm
      REAL*8 :: dym, dyp, dypp, indypp, indyp, indym
      REAL*8 :: dzp, indzp, dzm, indzm, dzpp, indzpp
      REAL*8 :: dens_c, dens_e, dens_n, dens_t
      REAL*8 :: dens_w, dens_s, dens_b
      REAL*8 :: dens_ee, dens_nn, dens_tt
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
! ... values of density interpolated linearly on the staggered grid
!
      dens_c = (dz(k+1) * dens%c + dz(k) * dens%t) * indzp
      dens_e = (dz(k+1) * dens%e + dz(k) * dens%et) * indzp
      dens_n = (dz(k+1) * dens%n + dz(k) * dens%nt) * indzp
      dens_t = (dz(k+2) * dens%t + dz(k+1) * dens%tt) * indzpp
      dens_w = (dz(k+1) * dens%w  + dz(k) * dens%wt) * indzp
      dens_s = (dz(k+1) * dens%s  + dz(k) * dens%st) * indzp
      dens_b = (dz(k-1) * dens%c  + dz(k) * dens%b) * indzm
!
      dens_ee = dens_e
      dens_nn = dens_n
      dens_tt = dens_t
!
! ... On boundary mantain first order accuracy (1st order Upwind).
!
! ... on West volume bondary
!
      IF (fl_l(imjk) /= 1) THEN
        cs = (dz(k+1)*u%w + dz(k)*u%wt) * indzp
        IF ( cs >= 0.D0 ) fw = dens_w * w%w * cs
        IF ( cs <  0.D0 ) fw = dens_c * w%c * cs
      END IF
!
! ... on South volume bondary
!
      IF (fl_l(ijmk) /= 1) THEN
        cs = (dz(k+1)*v%s + dz(k)*v%st) * indzp
        IF ( cs >= 0.D0 ) fs = dens_s * w%s * cs
        IF ( cs <  0.D0 ) fs = dens_c * w%c * cs
      END IF
!
! ... on Bottom volume bondary
!
      IF (fl_l(ijkm) /= 1) THEN
        cs=0.5D0*(w%b+w%c)
        IF ( cs >= 0.D0 ) fb = dens_b * w%b * cs
        IF ( cs <  0.D0 ) fb = dens_c * w%c * cs
      END IF
!
! ... MUSCL reconstruction of momentum
!
! ... on East volume boundary
!
      gradc = (2.0 * indxp * (dens_e * w%e - dens_c * w%c))
      gradw = (2.0 * indxm * (dens_c * w%c - dens_w * w%w))
      grade = (2.0 * indxpp * (dens_ee * w%ee - dens_e * w%e))
!
      lim = 0.D0
      erre = 0.D0
!
      cs = (dz(k+1)*u%c+dz(k)*u%t)*indzp
      cn = cs * dt * 2.0 * indxp
      IF (cs >= 0.D0) THEN
	erre = gradw / gradc
        fou  = dens_c * w%c
	incr = 0.5D0 * dx(i)
      ELSE IF (cs < 0.D0) THEN
	erre = grade / gradc
        fou  = dens_e * w%e
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
      gradc = (2.0 * indyp * (dens_n * w%n - dens_c * w%c))
      grads = (2.0 * indym * (dens_c * w%c - dens_s * w%s))
      gradn = (2.0 * indypp * (dens_nn * w%nn - dens_n * w%n))
!
      lim = 0.D0
      erre = 0.D0
!
      cs = (dz(k+1)*v%c+dz(k)*v%t)*indzp
      cn = cs * dt * 2.0 * indyp
      IF (cs >= 0.D0) THEN
	erre = grads / gradc
        fou  = dens_c * w%c
	incr = 0.5D0 * dy(j)
      ELSE IF (cs < 0.D0) THEN
	erre = gradn / gradc
        fou  = dens_n * w%n
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
      gradc = (indz(k+1) * (dens_t * w%t - dens_c * w%c))
      gradt = (indz(k+2) * (dens_tt * w%tt - dens_t * w%t))
      gradb = (indz(k) * (dens_c * w%c - dens_b * w%b))
!
      lim = 0.D0
      erre = 0.D0
!
      cs = 0.5D0*(w%c+w%t)
      cn = cs * dt * indz(k+1)
      IF (cs >= 0.D0) THEN
	erre = gradb / gradc
        fou  = dens_c * w%c
	incr = 0.5D0 * dz(k+1)
      ELSE IF (cs < 0.D0) THEN
	erre = gradt / gradc
        fou  = dens_t * w%t
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
      END SUBROUTINE flw_3d
!------------------------------------------------------
      SUBROUTINE flw_2d(fe, ft, fw, fb, dens, u, w, ij)
!
! ... Compute the convective fluxes on East, Top, sides of the cell
! ... for the momentum density along z.
!
      USE dimensions
      USE domain_decomposition, ONLY: myijk
      USE grid, ONLY: dx, xb, dz, indz, fl_l
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: imjk, ijkm
      USE set_indexes, ONLY: stencil
      USE time_parameters, ONLY: dt

      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: fe, ft, fw, fb
      TYPE(stencil), INTENT(IN) :: dens, u, w
      INTEGER, INTENT(IN) :: ij
      INTEGER :: i,j,imesh
!
      REAL*8 :: dens_c, dens_e, dens_t
      REAL*8 :: dens_w, dens_b
      REAL*8 :: dens_ee, dens_tt
      REAL*8 :: dxm, dxp, dxpp, indxpp, indxp, indxm
      REAL*8 :: dzp, indzp, dzm, indzm, dzpp, indzpp
      REAL*8 :: gradc, grade, gradw, gradt, gradb
!
      imesh = myijk( ip0_jp0_kp0_, ij)
      j = ( imesh - 1 ) / nx + 1
      i = MOD( ( imesh - 1 ), nx) + 1
!
      dxm=dx(i)+dx(i-1)
      dxp=dx(i)+dx(i+1)
      dxpp=dx(i+1)+dx(i+2)
      dzm=dz(j)+dz(j-1)
      dzp=dz(j)+dz(j+1)
      dzpp=dz(j+1)+dz(j+2)

      indxm=1.D0/dxm
      indxp=1.D0/dxp
      indxpp=1.D0/dxpp
      indzm=1.D0/dzm
      indzp=1.D0/dzp
      indzpp=1.D0/dzpp
!       
! ... Compute linearly interpolated values of density on the staggered grid
!
      dens_c = (dz(j+1) * dens%c + dz(j) * dens%t) * indzp
      dens_t = (dz(j+2) * dens%t + dz(j+1) * dens%tt) * indzpp
      dens_e = (dz(j+1) * dens%e + dz(j) * dens%en) * indzp
      dens_b = (dz(j)   * dens%b + dz(j-1) * dens%c) * indzm
      dens_w = (dz(j+1) * dens%w + dz(j) * dens%wn) * indzp
!
! ... an arbitrary choice !
!
      dens_ee = dens_e
      dens_tt = dens_t
!
! ... On boundary mantain first order accuracy (1st order Upwind).
!
! ... on West volume bondary
!
      IF( fl_l(imjk) /= 1 ) THEN
        cs = (u%wn * dz(j) + u%w * dz(j+1)) * indzp
        IF ( cs >= 0.D0 ) fw = dens_w * w%w * cs * xb(i-1)
        IF ( cs <  0.D0 ) fw = dens_c * w%c * cs * xb(i-1)
      END IF
!
! ... on Bottom volume bondary
!
      IF( fl_l(ijkm) /= 1 ) THEN
        cs = 0.5D0 * ( w%c + w%b ) 
        IF ( cs >= 0.D0 ) fb = dens_b * w%b * cs
        IF ( cs <  0.D0 ) fb = dens_c * w%c * cs
      END IF
!
! ... MUSCL reconstruction of momentum
!
! ... on East volume boundary
!
      gradc = 2.D0 * indxp  * (dens_e * w%e   - dens_c * w%c)
      gradw = 2.D0 * indxm  * (dens_c * w%c   - dens_w * w%w)
      grade = 2.D0 * indxpp * (dens_ee * w%ee - dens_e * w%e)
!
      lim = 0.D0
      erre = 0.D0
!
      cs = indzp * (u%c * dz(j+1) + u%t * dz(j))
      cn = cs * dt * 2.D0 * indzp
      IF ( cs >= 0.D0 ) THEN
	erre = gradw / gradc
        fou  = dens_c * w%c
	incr = 0.5D0 * dx(i)
      ELSE IF ( cs < 0.D0 ) THEN
	erre = grade / gradc
        fou  = dens_e * w%e
	incr = 0.5D0 * dx(i+1)
      END IF
!
      IF ((muscl /= 0) .AND. (gradc /= 0.D0) .AND. (i /= nx-1)) THEN
        CALL limiters(lim,erre)
      END IF
!
      upwnd = fou + lim * gradc * incr
!
      fe = upwnd * cs * xb(i)
!
! ... on Top volume boundary
!
      gradc = (dens_t * w%t   - dens_c * w%c) * indz(j+1)
      gradb = (dens_c * w%c   - dens_b * w%b) * indz(j)
      gradt = (dens_tt * w%tt - dens_t * w%t) * indz(j+2)
!
      lim = 0.D0
      erre = 0.D0
!
      cs = 0.5D0 * ( w%t + w%c )
      cn = cs * dt * indz(j+1)
      IF (cs >= 0.D0) THEN
	erre = gradb / gradc
        fou  = dens_c * w%c
	incr = 0.5D0 * dz(j+1)
      ELSE IF (cs < 0.D0) THEN
	erre = gradt / gradc
        fou  = dens_t * w%t
	incr = 0.5D0 * dz(j+1)
      END IF 
!
      IF ((muscl /= 0) .AND. (gradc /= 0.D0) .AND. (j /= nz-1)) THEN
        CALL limiters(lim,erre)
      END IF
!
      upwnd = fou + lim * gradc * incr
!
      ft = upwnd * cs
!
      RETURN
      END SUBROUTINE flw_2d
!----------------------------------------------------------------------
      END MODULE convective_fluxes_w
!-----------------------------------------------------------------------
