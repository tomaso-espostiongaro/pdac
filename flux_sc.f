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
      USE flux_limiters, ONLY: muscl, limiters

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
      
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE fsc_3d(fe, fn, ft, fw, fs, fb, dens, field, u, v, w, ijk)
!
      USE dimensions
      USE domain_decomposition, ONLY: myijk
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
      IF ((flag(imjk) /= 1)) THEN
        cs = u%w
        IF (cs >= 0.D0) THEN
          upwnd = dens%w * field%w
        ELSE IF (cs < 0.D0) THEN
          upwnd = dens%c * field%c
        ENDIF
        fw = upwnd * cs
      END IF
!
! ... on South volume boundary
!
      IF ((flag(ijmk) /= 1)) THEN
        cs = v%s
        IF (cs >= 0.D0) THEN
          upwnd = dens%s * field%s
        ELSE IF (cs < 0.D0) THEN
          upwnd = dens%c * field%c
        ENDIF
        fs = upwnd * cs
      END IF
!
! ... on Bottom volume boundary
!
      IF ((flag(ijkm) /= 1)) THEN
        cs = w%b
        IF (cs >= 0.D0) THEN
          upwnd = dens%b * field%b
        ELSE IF (cs < 0.D0) THEN
          upwnd = dens%c * field%c
        ENDIF
        fb = upwnd * cs
      END IF
!
! ... on East volume boundary
!
      cs = u%c
      !cn = cs * dt * 2.0 * indxp
      IF (cs >= 0.D0) THEN
        upwnd  = dens%c*field%c 
      ELSE IF (cs < 0.D0) THEN
        upwnd  = dens%e*field%e 
      END IF
      fe = upwnd * cs
!
! ... on North volume boundary
!
      cs = v%c
      !cn = cs * dt * 2.0 * indyp
      IF (cs >= 0.D0) THEN
        upwnd  = dens%c*field%c 
      ELSE IF (cs < 0.D0) THEN
        upwnd  = dens%n*field%n 
      ENDIF
      fn = upwnd * cs
!
! ... on Top volume boundary
!
      cs = w%c
      !cn = cs * dt * 2.0 * indzp
      IF (cs >= 0.D0) THEN
        upwnd  = dens%c*field%c 
      ELSE IF (cs < 0.D0) THEN
        upwnd  = dens%t*field%t 
      ENDIF
      ft = upwnd * cs
!
      RETURN
      END SUBROUTINE fsc_3d
!----------------------------------------------------------------------
      SUBROUTINE muscl_fsc_3d(fe, fn, ft, dens, field, u, v, w, ijk)
!
      USE dimensions
      USE domain_decomposition, ONLY: myijk
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
      im2 = MAX( 1, i-2 )
      ip2 = MIN( nx, i+2 )
      jm2 = MAX( 1, j-2 )
      jp2 = MIN( ny, j+2 )
      km2 = MAX( 1, k-2 )
      kp2 = MIN( nz, k+2 )
!
      dxm=dx(i)+dx(i-1)
      dxp=dx(i)+dx(i+1)
      dxpp=dx(i+1)+dx(ip2)
      dym=dy(j)+dy(j-1)
      dyp=dy(j)+dy(j+1)
      dypp=dy(j+1)+dy(jp2)
      dzm=dz(k)+dz(k-1)
      dzp=dz(k)+dz(k+1)
      dzpp=dz(k+1)+dz(kp2)

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
! ... On boundaries mantain first order accuracy
! ... Otherwise, MUSCL reconstruction of fields
!
! ... on East volume boundary
!
      gradc = 2.0 * indxp * (dens%e*field%e - dens%c*field%c)
      gradw = 2.0 * indxm * (dens%c*field%c - dens%w*field%w)
      grade = 2.0 * indxpp * (dens%ee*field%ee - dens%e*field%e)
!
      lim = 0.D0
      erre = 0.D0
!
      cs = u%c
      !cn = cs * dt * 2.0 * indxp
      IF (cs >= 0.D0) THEN
	erre = gradw / gradc
	incr = 0.5D0 * dx(i)
      ELSE IF (cs < 0.D0) THEN
	erre = grade / gradc
	incr = 0.5D0 * dx(i+1)
      ENDIF
!
      IF ((gradc /= 0.D0) .AND. (i /= nx-1)) THEN
        CALL limiters(lim,erre)
      END IF
!
      upwnd = lim * gradc * incr
!
      ! ... MUSCL correction
      fe = fe + upwnd * cs
!
! ... on North volume boundary
!
      gradc = 2.0 * indyp * (dens%n*field%n - dens%c*field%c)
      grads = 2.0 * indym * (dens%c*field%c - dens%s*field%s)
      gradn = 2.0 * indypp * (dens%nn*field%nn - dens%n*field%n)
!
      lim = 0.D0
      erre = 0.D0
!
      cs = v%c
      !cn = cs * dt * 2.0 * indyp
      IF (cs >= 0.D0) THEN
	erre = grads / gradc
	incr = 0.5D0 * dy(j)
      ELSE IF (cs < 0.D0) THEN
	erre = gradn / gradc
	incr = 0.5D0 * dy(j+1)
      ENDIF
!
      IF ((gradc /= 0.D0) .AND. (j /= ny-1)) THEN
        CALL limiters(lim,erre)
      END IF
!
      upwnd = lim * gradc * incr
!
      ! ... MUSCL correction
      fn = fn + upwnd * cs
!
! ... on Top volume boundary
!
      gradc = 2.0 * indzp * (dens%t*field%t - dens%c*field%c)
      gradb = 2.0 * indzm * (dens%c*field%c - dens%b*field%b)
      gradt = 2.0 * indzpp * (dens%tt*field%tt - dens%t*field%t)
!
      lim = 0.D0
      erre = 0.D0
!
      cs = w%c
      !cn = cs * dt * 2.0 * indzp
      IF (cs >= 0.D0) THEN
	erre = gradb / gradc
	incr = 0.5D0 * dz(k)
      ELSE IF (cs < 0.D0) THEN
	erre = gradt / gradc
	incr = 0.5D0 * dz(k+1)
      ENDIF
!
      IF ((gradc /= 0.D0) .AND. (k /= nz-1)) THEN
        CALL limiters(lim,erre)
      END IF
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
      USE domain_decomposition, ONLY: myijk
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
      IF ((flag(imjk) /= 1)) THEN
        cs = u%w
        IF (cs >= 0.D0) THEN
          upwnd = dens%w * field%w
        ELSE IF (cs < 0.D0) THEN
          upwnd = dens%c * field%c
        ENDIF
        fw = upwnd * cs * rb(i-1)
      END IF
!
! ... on Bottom volume boundary
!
      IF ((flag(ijkm) /= 1)) THEN
        cs = w%b
        IF (cs >= 0.D0) THEN
          upwnd = dens%b * field%b
        ELSE IF (cs < 0.D0) THEN
          upwnd = dens%c * field%c
        ENDIF
        fb = upwnd * cs
      END IF
!
! ... on East volume boundary
!
      cs = u%c
      !cn = cs * dt * 2.D0 * indxp
      IF (cs >= 0.D0) THEN
        upwnd  = dens%c*field%c 
      ELSE IF (cs < 0.D0) THEN
        upwnd  = dens%e*field%e 
      ENDIF
      fe = upwnd * cs * rb(i)
!
! ... on Top volume boundary
!
      cs = w%c
      !cn = cs * dt * 2.0 * indzp
      IF (cs >= 0.D0) THEN
        upwnd  = dens%c*field%c 
      ELSE IF (cs < 0.D0) THEN
        upwnd  = dens%t*field%t 
      ENDIF
      ft = upwnd * cs
!
      RETURN
      END SUBROUTINE fsc_2d
!-----------------------------------------------------------------------
      SUBROUTINE muscl_fsc_2d(fe, ft, dens, field, u, w, ijk)
!
      USE dimensions
      USE domain_decomposition, ONLY: myijk
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
      im2 = MAX( 1, i-2 )
      ip2 = MIN( nx, i+2 )
      km2 = MAX( 1, k-2 )
      kp2 = MIN( nz, k+2 )

!
      dxm=dx(i)+dx(i-1)
      dxp=dx(i)+dx(i+1)
      dxpp=dx(i+1)+dx(ip2)
      dzm=dz(k)+dz(k-1)
      dzp=dz(k)+dz(k+1)
      dzpp=dz(k+1)+dz(kp2)

      indxm=1.D0/dxm
      indxp=1.D0/dxp
      indxpp=1.D0/dxpp
      indzm=1.D0/dzm
      indzp=1.D0/dzp
      indzpp=1.D0/dzpp
!
! ... On boundaries mantain first order accuracy (no MUSCL correction)
! ... Otherwise, correct fluxes by using a linear interpolation of
! ... gradients
!
! ... on East volume boundary
!
      gradc = 2.D0 * indxp * (dens%e*field%e - dens%c*field%c)
      gradw = 2.D0 * indxm * (dens%c*field%c - dens%w*field%w)
      grade = 2.D0 * indxpp * (dens%ee*field%ee - dens%e*field%e)
!
      lim = 0.D0
      erre = 0.D0
!
      cs = u%c
      !cn = cs * dt * 2.D0 * indxp
      IF (cs >= 0.D0) THEN
	erre = gradw / gradc
	incr = 0.5D0 * dx(i)
      ELSE IF (cs < 0.D0) THEN
	erre = grade / gradc
	incr = 0.5D0 * dx(i+1)
      ENDIF
!
      IF ((gradc /= 0.D0) .AND. (i /= nx-1)) THEN
        CALL limiters(lim,erre)
      END IF
!
      upwnd = lim * gradc * incr
!
      ! ... MUSCL correction
      fe = fe + upwnd * cs * rb(i)
!
! ... on Top volume boundary
!
      gradc = 2.D0 * indzp * (dens%t*field%t - dens%c*field%c)
      gradb = 2.D0 * indzm * (dens%c*field%c - dens%b*field%b)
      gradt = 2.D0 * indzpp * (dens%tt*field%tt - dens%t*field%t)
!
      lim = 0.D0
      erre = 0.D0
!
      cs = w%c
      !cn = cs * dt * 2.0 * indzp
      IF (cs >= 0.D0) THEN
	erre = gradb / gradc
	incr = 0.5D0 * dz(k)
      ELSE IF (cs < 0.D0) THEN
	erre = gradt / gradc
	incr = 0.5D0 * dz(k+1)
      ENDIF
!
      IF ((gradc /= 0.D0) .AND. (k /= nz-1)) THEN
        CALL limiters(lim,erre)
      END IF
!
      upwnd = lim * gradc * incr
!
      ! ... MUSCL correction
      ft = ft + upwnd * cs
!
      RETURN
      END SUBROUTINE muscl_fsc_2d
!-----------------------------------------------------------------------
      END MODULE convective_fluxes_sc
!-----------------------------------------------------------------------
