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
      USE domain_decomposition, ONLY: myijk, meshinds
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
      IF( flag(imjk) /= 1 ) THEN
        cs = 0.5D0 * ( u%c + u%w )
        IF ( cs >= 0.D0 ) fw = dens%w * u%w * cs
        IF ( cs <  0.D0 ) fw = dens%c * u%c * cs
      END IF
!
! ... on South volume bondary
!
      IF( flag(ijmk) /= 1 ) THEN
        cs = ( dx(i+1) * v%s + dx(i) * v%es ) * indxp
        IF ( cs >= 0.D0 ) fs = dens%s * u%s * cs
        IF ( cs <  0.D0 ) fs = dens%c * u%c * cs
      END IF
!
! ... on Bottom volume bondary
!
      IF( flag(ijkm) /= 1 ) THEN
        cs = ( dx(i+1) * w%b + dx(i) * w%eb ) * indxp
        IF ( cs >= 0.D0 ) fb = dens%b * u%b * cs
        IF ( cs <  0.D0 ) fb = dens%c * u%c * cs
      END IF
!
! ... on East volume boundary
!
      cs = 0.5D0 * (u%c + u%e)    
      IF ( cs >= 0.D0 ) fe  = dens%c * u%c * cs
      IF ( cs <  0.D0 ) fe  = dens%e * u%e * cs
!
! ... on North volume boundary
!
      cs = (dx(i+1) * v%c + dx(i) * v%e) * indxp
      IF ( cs >= 0.D0 ) fn  = dens%c * u%c * cs
      IF ( cs <  0.D0 ) fn  = dens%n * u%n * cs
!
! ... on Top volume boundary
!
      cs = (dx(i+1) * w%c + dx(i) * w%e) * indxp
      IF ( cs >= 0.D0 ) ft  = dens%c * u%c * cs
      IF ( cs <  0.D0 ) ft  = dens%t * u%t * cs
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
      USE flux_limiters, ONLY: limiters
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
      dxp=dx(i)+dx(i+1)
      dym=dy(j)+dy(j-1)
      dyp=dy(j)+dy(j+1)
      dypp=dy(j+1)+dy(jp2)
      dzm=dz(k)+dz(k-1)
      dzp=dz(k)+dz(k+1)
      dzpp=dz(k+1)+dz(kp2)

      indxp=1.D0/dxp
      indym=1.D0/dym
      indyp=1.D0/dyp
      indypp=1.D0/dypp
      indzm=1.D0/dzm
      indzp=1.D0/dzp
      indzpp=1.D0/dzpp
!
! ... MUSCL reconstruction of momentum
!
! ... on East volume boundary
!
      gradc = (indx(i+1) * (dens%e * u%e   - dens%c * u%c))
      gradw = (indx(i)   * (dens%c * u%c   - dens%w * u%w))
      grade = (indx(ip2) * (dens%ee * u%ee - dens%e * u%e))
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
      CALL limiters(lim,erre)
!
      upwnd = lim * gradc * incr
!
      fe = fe + upwnd * cs
!
! ... on North volume boundary
!
      gradc = (dens%n * u%n   - dens%c * u%c) * 2.D0 * indyp
      grads = (dens%c * u%c   - dens%s * u%s) * 2.D0 * indym
      gradn = (dens%nn * u%nn - dens%n * u%n) * 2.D0 * indypp
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
      CALL limiters(lim,erre)
!
      upwnd = lim * gradc * incr
!
      fn = fn + upwnd * cs
!
! ... on Top volume boundary
!
      gradc = (dens%t * u%t   - dens%c * u%c) * 2.D0 * indzp
      gradb = (dens%c * u%c   - dens%b * u%b) * 2.D0 * indzm
      gradt = (dens%tt * u%tt - dens%t * u%t) * 2.D0 * indzpp
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
      CALL limiters(lim,erre)
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
      dxp=dx(i)+dx(i+1)
      indxp=1.D0/dxp
!
! ... on West volume bondary
!
      IF( flag(imjk) /= 1 ) THEN
        cs = 0.5D0*(u%c + u%w)
        IF ( cs >= 0.D0 ) fw = dens%w * u%w * cs * r(i)
        IF ( cs <  0.D0 ) fw = dens%c * u%c * cs * r(i)
      END IF
!
! ... on Bottom volume bondary
!
      IF( flag(ijkm) /= 1 ) THEN
        cs = (dx(i+1) * w%b + dx(i) * w%eb) * indxp
        IF ( cs >= 0.D0 ) fb = dens%b * u%b * cs
        IF ( cs <  0.D0 ) fb = dens%c * u%c * cs
      END IF
!
! ... on East volume boundary
!
      cs = 0.5D0 * (u%c + u%e)
      IF ( cs >= 0.D0 ) fe = dens%c * u%c * cs * r(i+1)
      IF ( cs <  0.D0 ) fe = dens%e * u%e * cs * r(i+1)
!
! ... on Top volume boundary
!
      cs = (dx(i+1) * w%c + dx(i) * w%e) * indxp
      IF (cs >= 0.D0) ft  = dens%c * u%c * cs
      IF (cs <  0.D0) ft  = dens%t * u%t * cs
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
      USE flux_limiters, ONLY: limiters
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
      dxp=dx(i)+dx(i+1)
      dzm=dz(k)+dz(k-1)
      dzp=dz(k)+dz(k+1)
      dzpp=dz(k+1)+dz(kp2)

      indxp=1.D0/dxp
      indzm=1.D0/dzm
      indzp=1.D0/dzp
      indzpp=1.D0/dzpp
!
! ... MUSCL reconstruction of momentum
!
! ... on East volume boundary
!
      gradc = (indx(i+1) * (dens%e * u%e   - dens%c * u%c))
      gradw = (indx(i)   * (dens%c * u%c   - dens%w * u%w))
      grade = (indx(ip2) * (dens%ee * u%ee - dens%e * u%e))
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
      CALL limiters(lim,erre)
!
      upwnd = lim * gradc * incr
!
      fe = fe + upwnd * cs * r(i+1)
!
! ... on Top volume boundary
!
      gradc = (dens%t * u%t   - dens%c * u%c) * 2.D0 * indzp
      gradb = (dens%c * u%c   - dens%b * u%b) * 2.D0 * indzm
      gradt = (dens%tt * u%tt - dens%t * u%t) * 2.D0 * indzpp
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
      CALL limiters(lim,erre)
!
      upwnd = lim * gradc * incr
!
      ft = ft + upwnd * cs
!
      RETURN
      END SUBROUTINE muscl_flu_2d
!-----------------------------------------------------------------------
      END MODULE convective_fluxes_u
!-----------------------------------------------------------------------
