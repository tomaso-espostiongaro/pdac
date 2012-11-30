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
      USE io_files, ONLY: testunit
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

      dyp = dy(j)+dy(j+1)
      indyp = 1.D0/dyp
!
! ... on West volume bondary
!
      IF( .NOT.BTEST(flag(imjk),0) ) THEN
        cs = (u%wn * dy(j) + u%w * dy(j+1)) * indyp
        fw = 0.5D0*(cs-ABS(cs))*dens%c*v%c + 0.5D0*(cs+ABS(cs))*dens%w*v%w 
      END IF
!
! ... on South volume bondary
!
      IF( .NOT.BTEST(flag(ijmk),0) ) THEN
        cs = 0.5D0 * ( v%c + v%s ) 
        fs = 0.5D0*(cs-ABS(cs))*dens%c*v%c + 0.5D0*(cs+ABS(cs))*dens%s*v%s 
      END IF
!
! ... on Bottom volume bondary
!
      IF( .NOT.BTEST(flag(ijkm),0) ) THEN
        cs = (w%nb * dy(j) + w%b * dy(j+1)) * indyp
        fb = 0.5D0*(cs-ABS(cs))*dens%c*v%c + 0.5D0*(cs+ABS(cs))*dens%b*v%b 
      END IF
!
! ... on East volume boundary
!
      cs = indyp * (u%c * dy(j+1) + u%n * dy(j))
      fe = 0.5D0*(cs-ABS(cs))*dens%e*v%e + 0.5D0*(cs+ABS(cs))*dens%c*v%c 
!
! ... on North volume boundary
!
      cs = 0.5D0 * ( v%n + v%c )
      fn = 0.5D0*(cs-ABS(cs))*dens%n*v%n + 0.5D0*(cs+ABS(cs))*dens%c*v%c 
!
! ... on Top volume boundary
!
      cs = indyp * (w%c * dy(j+1) + w%n * dy(j))
      ft = 0.5D0*(cs-ABS(cs))*dens%t*v%t + 0.5D0*(cs+ABS(cs))*dens%c*v%c 
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
      USE flux_limiters, ONLY: limiters, lv
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

      INTEGER :: ip2, jp2, kp2
!
      ip2 = MIN( nx, i+2 )
      jp2 = MIN( ny, j+2 )
      kp2 = MIN( nz, k+2 )

      dyp  = dy(j)  + dy(j+1)
      dxm  = dx(i)  + dx(i-1)
      dxp  = dx(i)  + dx(i+1)
      dxpp = dx(i+1)+ dx(ip2)
      dzm  = dz(k)  + dz(k-1)
      dzp  = dz(k)  + dz(k+1)
      dzpp = dz(k+1)+ dz(kp2)

      indyp  = 1.D0/dyp
      indxm  = 1.D0/dxm
      indxp  = 1.D0/dxp
      indxpp = 1.D0/dxpp
      indzm  = 1.D0/dzm
      indzp  = 1.D0/dzp
      indzpp = 1.D0/dzpp
!
! ... MUSCL reconstruction of momentum
!
! ... on East volume boundary
!
      gradw = 2.D0 * indxm  * (v%c *dens%c  - v%w*dens%w)
      gradc = 2.D0 * indxp  * (v%e *dens%e  - v%c*dens%c)
      grade = 2.D0 * indxpp * (v%ee*dens%ee - v%e*dens%e)
!
      lim  = 0.D0
      erre = 0.D0
!
      cs = indyp * (u%c * dy(j+1) + u%n * dy(j))
      IF ( cs >= 0.D0 ) THEN
        IF (gradc /= 0) erre = gradw / gradc
	incr = 0.5D0 * dx(i)
      ELSE IF ( cs < 0.D0 ) THEN
        IF (gradc /= 0) erre = grade / gradc
	incr = - 0.5D0 * dx(i+1)
      END IF
!
      CALL limiters(lv,lim,erre)
!
      upwnd = lim * gradc * incr
!
      fe = fe + upwnd * cs
!
! ... on North volume boundary
!
      grads = (v%c *dens%c  - v%s*dens%s) * indy(j)
      gradc = (v%n *dens%n  - v%c*dens%c) * indy(j+1)
      gradn = (v%nn*dens%nn - v%n*dens%n) * indy(jp2)
!
      lim  = 0.D0
      erre = 0.D0
!
      cs = 0.5D0 * ( v%n + v%c )
      IF (cs >= 0.D0) THEN
        IF (gradc /= 0) erre = grads / gradc
	incr = 0.5D0 * dy(j+1)
      ELSE IF (cs < 0.D0) THEN
        IF (gradc /= 0) erre = gradn / gradc
	incr = - 0.5D0 * dy(j+1)
      END IF 
!
      CALL limiters(lv,lim,erre)
!
      upwnd = lim * gradc * incr
!
      fn = fn + upwnd * cs
!
! ... on Top volume boundary
!
      gradb = (v%c *dens%c  - v%b*dens%b) * 2.D0 * indzm
      gradc = (v%t *dens%t  - v%c*dens%c) * 2.D0 * indzp
      gradt = (v%tt*dens%tt - v%t*dens%t) * 2.D0 * indzpp
!
      lim  = 0.D0
      erre = 0.D0
!
      cs = (dy(j+1) * w%c + dy(j) * w%n) * indyp
      IF (cs >= 0.D0) THEN
        IF (gradc /= 0) erre = gradb / gradc
	incr = 0.5D0 * dz(k)
      ELSE IF (cs < 0.D0) THEN
        IF (gradc /= 0) erre = gradt / gradc
	incr = - 0.5D0 * dz(k+1)
      END IF 
!
      CALL limiters(lv,lim,erre)
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
