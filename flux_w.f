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
      IF (flag(imjk) /= 1) THEN
        cs = (dz(k+1)*u%w + dz(k)*u%wt) * indzp
        IF ( cs >= 0.D0 ) fw = dens%w * w%w * cs
        IF ( cs <  0.D0 ) fw = dens%c * w%c * cs
      END IF
!
! ... on South volume bondary
!
      IF (flag(ijmk) /= 1) THEN
        cs = (dz(k+1)*v%s + dz(k)*v%st) * indzp
        IF ( cs >= 0.D0 ) fs = dens%s * w%s * cs
        IF ( cs <  0.D0 ) fs = dens%c * w%c * cs
      END IF
!
! ... on Bottom volume bondary
!
      IF (flag(ijkm) /= 1) THEN
        cs=0.5D0*(w%b+w%c)
        IF ( cs >= 0.D0 ) fb = dens%b * w%b * cs
        IF ( cs <  0.D0 ) fb = dens%c * w%c * cs
      END IF
!
! ... on East volume boundary
!
      cs = (dz(k+1)*u%c+dz(k)*u%t)*indzp
      IF (cs >= 0.D0) fe  = dens%c * w%c * cs
      IF (cs <  0.D0) fe  = dens%e * w%e * cs
!
! ... on North volume boundary
!
      cs = (dz(k+1)*v%c+dz(k)*v%t)*indzp
      IF (cs >= 0.D0) fn  = dens%c * w%c * cs
      IF (cs <  0.D0) fn  = dens%n * w%n * cs
!
! ... on Top volume boundary
!
      cs = 0.5D0*(w%c+w%t)
      IF (cs >= 0.D0) ft  = dens%c * w%c * cs
      IF (cs <  0.D0) ft  = dens%t * w%t * cs
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
      USE flux_limiters, ONLY: limiters
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
      dxm=dx(i)+dx(i-1)
      dxp=dx(i)+dx(i+1)
      dxpp=dx(i+1)+dx(ip2)
      dym=dy(j)+dy(j-1)
      dyp=dy(j)+dy(j+1)
      dypp=dy(j+1)+dy(jp2)
      dzp=dz(k)+dz(k+1)

      indxm=1.D0/dxm
      indxp=1.D0/dxp
      indxpp=1.D0/dxpp
      indym=1.D0/dym
      indyp=1.D0/dyp
      indypp=1.D0/dypp
      indzp=1.D0/dzp
!
! ... MUSCL reconstruction of momentum
!
! ... on East volume boundary
!
      gradc = (2.D0 * indxp * (dens%e * w%e - dens%c * w%c))
      gradw = (2.D0 * indxm * (dens%c * w%c - dens%w * w%w))
      grade = (2.D0 * indxpp * (dens%ee * w%ee - dens%e * w%e))
!
      lim = 0.D0
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
      CALL limiters(lim,erre)
!
      upwnd = lim * gradc * incr
!
      fe = fe + upwnd * cs
!
! ... on North volume boundary
!
      gradc = (2.D0 * indyp * (dens%n * w%n - dens%c * w%c))
      grads = (2.D0 * indym * (dens%c * w%c - dens%s * w%s))
      gradn = (2.D0 * indypp * (dens%nn * w%nn - dens%n * w%n))
!
      lim = 0.D0
      erre = 0.D0
!
      cs = ( dz(k+1) * v%c + dz(k) * v%t ) * indzp
      !cn = cs * dt * 2.D0 * indyp
      IF (cs >= 0.D0) THEN
        IF (gradc /= 0) erre = grads / gradc
        incr = 0.5D0 * dy(j)
      ELSE IF (cs < 0.D0) THEN
        IF (gradc /= 0) erre = gradn / gradc
        incr = - 0.5D0 * dy(j+1)
      END IF
!
      CALL limiters(lim,erre)
!
      upwnd = lim * gradc * incr
!
      fn = fn + upwnd * cs
!
! ... on Top volume boundary
!
      gradc = (indz(k+1) * (dens%t * w%t - dens%c * w%c))
      gradb = (indz(k) * (dens%c * w%c - dens%b * w%b))
      gradt = (indz(kp2) * (dens%tt * w%tt - dens%t * w%t))
!
      lim = 0.D0
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
      CALL limiters(lim,erre)
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
      IF( flag(imjk) /= 1 ) THEN
        cs = (u%wt * dz(k) + u%w * dz(k+1)) * indzp
        IF ( cs >= 0.D0 ) fw = dens%w * w%w * cs * rb(i-1)
        IF ( cs <  0.D0 ) fw = dens%c * w%c * cs * rb(i-1)
      END IF
!
! ... on Bottom volume bondary
!
      IF( flag(ijkm) /= 1 ) THEN
        cs = 0.5D0 * ( w%c + w%b ) 
        IF ( cs >= 0.D0 ) fb = dens%b * w%b * cs
        IF ( cs <  0.D0 ) fb = dens%c * w%c * cs
      END IF
!
! ... on East volume boundary
!
      cs = indzp * (u%c * dz(k+1) + u%t * dz(k))
      IF ( cs >= 0.D0 ) fe  = dens%c * w%c * cs * rb(i)
      IF ( cs <  0.D0 ) fe  = dens%e * w%e * cs * rb(i)
!
! ... on Top volume boundary
!
      cs = 0.5D0 * ( w%t + w%c )
      IF (cs >= 0.D0) ft  = dens%c * w%c * cs
      IF (cs <  0.D0) ft  = dens%t * w%t * cs
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
      USE flux_limiters, ONLY: limiters
      USE grid, ONLY: dx, rb, dz, indz, flag
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
      dxm=dx(i)+dx(i-1)
      dxp=dx(i)+dx(i+1)
      dxpp=dx(i+1)+dx(ip2)
      dzp=dz(k)+dz(k+1)

      indxm=1.D0/dxm
      indxp=1.D0/dxp
      indxpp=1.D0/dxpp
      indzp=1.D0/dzp
!
! ... MUSCL reconstruction of momentum
!
! ... on East volume boundary
!
      gradc = 2.D0 * indxp  * (dens%e * w%e   - dens%c * w%c)
      gradw = 2.D0 * indxm  * (dens%c * w%c   - dens%w * w%w)
      grade = 2.D0 * indxpp * (dens%ee * w%ee - dens%e * w%e)
!
      lim = 0.D0
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
      CALL limiters(lim,erre)
!
      upwnd = lim * gradc * incr
!
      fe = fe + upwnd * cs * rb(i)
!
! ... on Top volume boundary
!
      gradc = (dens%t * w%t   - dens%c * w%c) * indz(k+1)
      gradb = (dens%c * w%c   - dens%b * w%b) * indz(k)
      gradt = (dens%tt * w%tt - dens%t * w%t) * indz(kp2)
!
      lim = 0.D0
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
      CALL limiters(lim,erre)
!
      upwnd = lim * gradc * incr
!
      ft = ft + upwnd * cs
!
      RETURN
      END SUBROUTINE muscl_flw_2d
!----------------------------------------------------------------------
      END MODULE convective_fluxes_w
!-----------------------------------------------------------------------
