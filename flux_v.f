!----------------------------------------------------------------------
      MODULE convective_fluxes_v
!
! ... This module computes convective fluxes of _____ y-momentum ______
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
      INTERFACE muscl_flv
        MODULE PROCEDURE muscl_flv_3d
      END INTERFACE
      INTERFACE flv
        MODULE PROCEDURE flv_3d
      END INTERFACE
      
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE flv_3d(fe, fn, ft, fw, fs, fb, dens, u, v, w, j)
!
! ... Compute the convective fluxes on East, North, and Top sides of the cell
! ... for the momentum density along y.
!
      USE grid, ONLY: flag
      USE grid, ONLY: dy
      USE set_indexes, ONLY: imjk, ijmk, ijkm
      USE set_indexes, ONLY: stencil

      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: fe, fn, ft, fw, fs, fb
      TYPE(stencil), INTENT(IN) :: dens, u, v, w
      INTEGER, INTENT(IN) :: j
!
      REAL*8 :: dyp, indyp
      REAL*8 :: gradc, grade, gradw, gradn, grads, gradt, gradb

      dyp=dy(j)+dy(j+1)
      indyp=1.D0/dyp
!
! ... on West volume bondary
!
      IF( flag(imjk) /= 1 ) THEN
        cs = (u%wn * dy(j) + u%w * dy(j+1)) * indyp
        IF ( cs >= 0.D0 ) fw = dens%w * v%w * cs
        IF ( cs <  0.D0 ) fw = dens%c * v%c * cs
      END IF
!
! ... on South volume bondary
!
      IF( flag(ijmk) /= 1 ) THEN
        cs = 0.5D0 * ( v%c + v%s ) 
        IF ( cs >= 0.D0 ) fs = dens%s * v%s * cs
        IF ( cs <  0.D0 ) fs = dens%c * v%c * cs
      END IF
!
! ... on Bottom volume bondary
!
      IF( flag(ijkm) /= 1 ) THEN
        cs = (w%nb * dy(j) + w%b * dy(j+1)) * indyp
        IF ( cs >= 0.D0 ) fb = dens%b * v%b * cs
        IF ( cs <  0.D0 ) fb = dens%c * v%c * cs
      END IF
!
! ... on East volume boundary
!
      cs = indyp * (u%c * dy(j+1) + u%n * dy(j))
      IF ( cs >= 0.D0 ) fe  = dens%c * v%c * cs
      IF ( cs <  0.D0 ) fe  = dens%e * v%e * cs
!
! ... on North volume boundary
!
      cs = 0.5D0 * ( v%n + v%c )
      IF (cs >= 0.D0) fn  = dens%c * v%c * cs
      IF (cs <  0.D0) fn  = dens%n * v%n * cs
!
! ... on Top volume boundary
!
      cs = indyp * (w%c * dy(j+1) + w%n * dy(j))
      IF (cs >= 0.D0) ft  = dens%c * v%c * cs
      IF (cs <  0.D0) ft  = dens%t * v%t * cs
!
      RETURN
      END SUBROUTINE flv_3d
!----------------------------------------------------------------------
      SUBROUTINE muscl_flv_3d(fe, fn, ft, dens, u, v, w, i, j, k)
!
! ... Compute the convective fluxes on East, North, and Top sides of the cell
! ... for the momentum density along y.
!
      USE dimensions, ONLY: nx, ny, nz
      USE flux_limiters, ONLY: limiters
      USE grid, ONLY: flag
      USE grid, ONLY: dx, dy, dz, indy
      USE set_indexes, ONLY: stencil
      USE time_parameters, ONLY: dt

      IMPLICIT NONE
!
      REAL*8, INTENT(INOUT) :: fe, fn, ft
      TYPE(stencil), INTENT(IN) :: dens, u, v, w
      INTEGER, INTENT(IN) :: i, j, k

      REAL*8 :: dyp, indyp
      REAL*8 :: dxm, dxp, dxpp, indxpp, indxp, indxm
      REAL*8 :: dzp, indzp, dzm, indzm, dzpp, indzpp
      REAL*8 :: gradc, grade, gradw, gradn, grads, gradt, gradb
      REAL*8 :: dens0, dens1, dens2

      INTEGER :: ip2, jp2, kp2
!
      ip2 = MIN( nx, i+2 )
      jp2 = MIN( ny, j+2 )
      kp2 = MIN( nz, k+2 )

      dyp=dy(j)+dy(j+1)
      dxm=dx(i)+dx(i-1)
      dxp=dx(i)+dx(i+1)
      dxpp=dx(i+1)+dx(ip2)
      dzm=dz(k)+dz(k-1)
      dzp=dz(k)+dz(k+1)
      dzpp=dz(k+1)+dz(kp2)

      indyp=1.D0/dyp
      indxm=1.D0/dxm
      indxp=1.D0/dxp
      indxpp=1.D0/dxpp
      indzm=1.D0/dzm
      indzp=1.D0/dzp
      indzpp=1.D0/dzpp
!
! ... MUSCL reconstruction of momentum
!
! ... on East volume boundary
!
      dens0 = indxp * (dx(i)*dens%e + dx(i+1)*dens%c)

      gradw = 2.D0 * indxm  * dens%c * (v%c  - v%w)
      gradc = 2.D0 * indxp  * dens0 * (v%e  - v%c)
      grade = 2.D0 * indxpp * dens%e * (v%ee - v%e)
!
      lim = 0.D0
      erre = 0.D0
!
      cs = indyp * (u%c * dy(j+1) + u%n * dy(j))
      !cn = cs * dt * 2.D0 * indxp
      IF ( cs >= 0.D0 ) THEN
        IF (gradc /= 0) erre = gradw / gradc
	incr = 0.5D0 * dx(i)
      ELSE IF ( cs < 0.D0 ) THEN
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
      dens1 = 0.5D0 * (dens%n + dens%c)

      grads = dens%c * (v%c  - v%s) * indy(j)
      gradc = dens1 * (v%n  - v%c) * indy(j+1)
      gradn = dens%n * (v%nn - v%n) * indy(jp2)
!
      lim = 0.D0
      erre = 0.D0
!
      cs = 0.5D0 * ( v%n + v%c )
      !cn = cs * dt * indy(j+1)
      IF (cs >= 0.D0) THEN
        IF (gradc /= 0) erre = grads / gradc
	incr = 0.5D0 * dy(j+1)
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
      dens2 = indzp * (dz(k)*dens%t + dz(k+1)*dens%c)

      gradb = dens%c * (v%c  - v%b) * 2.D0 * indzm
      gradc = dens2 * (v%t  - v%c) * 2.D0 * indzp
      gradt = dens%t * (v%tt - v%t) * 2.D0 * indzpp
!
      lim = 0.D0
      erre = 0.D0
!
      cs = (dy(j+1) * w%c + dy(j) * w%n) * indyp
      !cn = cs * dt * 2.D0 * indzp
      IF (cs >= 0.D0) THEN
        IF (gradc /= 0) erre = gradb / gradc
	incr = 0.5D0 * dz(k)
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
      END SUBROUTINE muscl_flv_3d
!----------------------------------------------------------------------
      END MODULE convective_fluxes_v
!-----------------------------------------------------------------------
