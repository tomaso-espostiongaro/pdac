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
      INTERFACE flv
        MODULE PROCEDURE flv_3d
      END INTERFACE
      INTERFACE flv_1st
        MODULE PROCEDURE flv_3d_1st
      END INTERFACE
      
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE flv_3d(fe, fn, ft, fw, fs, fb, dens, u, v, w, ijk)
!
! ... Compute the convective fluxes on East, North, and Top sides of the cell
! ... for the momentum density along y.
!
      USE dimensions
      USE domain_decomposition, ONLY: myijk
      USE grid, ONLY: fl_l
      USE grid, ONLY: dx, dy, dz, indy
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
      INTEGER :: ip2, jp2, kp2
!
      imesh = myijk( ip0_jp0_kp0_, ijk)
      i = MOD( MOD( imesh - 1, nx*ny ), nx ) + 1
      j = MOD( imesh - 1, nx*ny ) / nx + 1
      k = ( imesh - 1 ) / ( nx*ny ) + 1
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
! ... Compute linearly interpolated values of density on the staggered grid
!
      dens_c = (dy(j+1) * dens%c + dy(j) * dens%n) * indyp
      dens_n = (dy(jp2) * dens%n + dy(j+1) * dens%nn) * indypp
      dens_e = (dy(j+1) * dens%e + dy(j) * dens%en) * indyp
      dens_t = (dy(j+1) * dens%t + dy(j) * dens%nt) * indyp
      dens_s = (dy(j)   * dens%s + dy(j-1) * dens%c) * indym
      dens_w = (dy(j+1) * dens%w + dy(j) * dens%wn) * indyp
      dens_b = (dy(j+1) * dens%b + dy(j) * dens%nb) * indyp
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
        cs = (u%wn * dy(j) + u%w * dy(j+1)) * indyp
        IF ( cs >= 0.D0 ) fw = dens_w * v%w * cs
        IF ( cs <  0.D0 ) fw = dens_c * v%c * cs
      END IF
!
! ... on South volume bondary
!
      IF( fl_l(ijmk) /= 1 ) THEN
        cs = 0.5D0 * ( v%c + v%s ) 
        IF ( cs >= 0.D0 ) fs = dens_s * v%s * cs
        IF ( cs <  0.D0 ) fs = dens_c * v%c * cs
      END IF
!
! ... on Bottom volume bondary
!
      IF( fl_l(ijkm) /= 1 ) THEN
        cs = ( w%b * dy(j+1) + w%nb * dy(j) ) * indyp
        IF ( cs >= 0.D0 ) fb = dens_b * v%b * cs
        IF ( cs <  0.D0 ) fb = dens_c * v%c * cs
      END IF
!
! ... MUSCL reconstruction of momentum
!
! ... on East volume boundary
!
      gradc = 2.D0 * indxp  * (dens_e * v%e   - dens_c * v%c)
      gradw = 2.D0 * indxm  * (dens_c * v%c   - dens_w * v%w)
      grade = 2.D0 * indxpp * (dens_ee * v%ee - dens_e * v%e)
!
      lim = 0.D0
      erre = 0.D0
!
      cs = indyp * (u%c * dy(j+1) + u%n * dy(j))
      cn = cs * dt * 2.D0 * indxp
      IF ( cs >= 0.D0 ) THEN
	erre = gradw / gradc
        fou  = dens_c * v%c
	incr = 0.5D0 * dx(i)
      ELSE IF ( cs < 0.D0 ) THEN
	erre = grade / gradc
        fou  = dens_e * v%e
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
      gradc = (dens_n * v%n   - dens_c * v%c) * indy(j+1)
      grads = (dens_c * v%c   - dens_s * v%s) * indy(j)
      gradn = (dens_nn * v%nn - dens_n * v%n) * indy(jp2)
!
      lim = 0.D0
      erre = 0.D0
!
      cs = 0.5D0 * ( v%n + v%c )
      cn = cs * dt * indy(j+1)
      IF (cs >= 0.D0) THEN
	erre = grads / gradc
        fou  = dens_c * v%c
	incr = 0.5D0 * dy(j+1)
      ELSE IF (cs < 0.D0) THEN
	erre = gradn / gradc
        fou  = dens_n * v%n
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
      gradc = (dens_t * v%t   - dens_c * v%c) * 2.0 * indzp
      gradb = (dens_c * v%c   - dens_b * v%b) * 2.0 * indzm
      gradt = (dens_tt * v%tt - dens_t * v%t) * 2.0 * indzpp
!
      lim = 0.D0
      erre = 0.D0
!
      cs = (dy(j+1) * w%c + dy(j) * w%n) * indyp
      cn = cs * dt * 2.0 * indzp
      IF (cs >= 0.D0) THEN
	erre = gradb / gradc
        fou  = dens_c * v%c
	incr = 0.5D0 * dz(k)
      ELSE IF (cs < 0.D0) THEN
	erre = gradt / gradc
        fou  = dens_t * v%t
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
      END SUBROUTINE flv_3d
!xxx
!----------------------------------------------------------------------
      SUBROUTINE flv_3d_1st(fe, fn, ft, fw, fs, fb, dens, u, v, w, j)
!
! ... Compute the convective fluxes on East, North, and Top sides of the cell
! ... for the momentum density along y.
!
      USE dimensions, ONLY: nx, ny, nz
      USE grid, ONLY: fl_l
      USE grid, ONLY: dx, dy, dz, indy
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: imjk, ijmk, ijkm
      USE set_indexes, ONLY: stencil
      USE time_parameters, ONLY: dt

      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: fe, fn, ft, fw, fs, fb
      TYPE(stencil), INTENT(IN) :: dens, u, v, w
      INTEGER, INTENT(IN) :: j
      
!
      REAL*8 :: dens_c, dens_e, dens_n, dens_t
      REAL*8 :: dens_w, dens_s, dens_b    
      REAL*8 :: dym, dyp, dypp, indypp, indyp, indym
      REAL*8 :: gradc, grade, gradw, gradn, grads, gradt, gradb

      INTEGER :: jp2

      jp2 = MIN( ny, j+2 )
    
      dym=dy(j)+dy(j-1)
      dyp=dy(j)+dy(j+1)
      dypp=dy(j+1)+dy(jp2)
     
      indym=1.D0/dym
      indyp=1.D0/dyp
      indypp=1.D0/dypp
      
!       
! ... Compute linearly interpolated values of density on the staggered grid
!
      dens_c = (dy(j+1) * dens%c + dy(j) * dens%n) * indyp
      dens_n = (dy(jp2) * dens%n + dy(j+1) * dens%nn) * indypp
      dens_e = (dy(j+1) * dens%e + dy(j) * dens%en) * indyp
      dens_t = (dy(j+1) * dens%t + dy(j) * dens%nt) * indyp
      dens_s = (dy(j)   * dens%s + dy(j-1) * dens%c) * indym
      dens_w = (dy(j+1) * dens%w + dy(j) * dens%wn) * indyp
      dens_b = (dy(j+1) * dens%b + dy(j) * dens%nb) * indyp
!
! ... On boundary mantain first order accuracy (1st order Upwind).
!
! ... on West volume bondary
!
      IF( fl_l(imjk) /= 1 ) THEN
        cs = (u%wn * dy(j) + u%w * dy(j+1)) * indyp
        IF ( cs >= 0.D0 ) fw = dens_w * v%w * cs
        IF ( cs <  0.D0 ) fw = dens_c * v%c * cs
      END IF
!
! ... on South volume bondary
!
      IF( fl_l(ijmk) /= 1 ) THEN
        cs = 0.5D0 * ( v%c + v%s ) 
        IF ( cs >= 0.D0 ) fs = dens_s * v%s * cs
        IF ( cs <  0.D0 ) fs = dens_c * v%c * cs
      END IF
!
! ... on Bottom volume bondary
!
      IF( fl_l(ijkm) /= 1 ) THEN
        cs = ( w%b * dy(j+1) + w%nb * dy(j) ) * indyp
        IF ( cs >= 0.D0 ) fb = dens_b * v%b * cs
        IF ( cs <  0.D0 ) fb = dens_c * v%c * cs
      END IF
!
! ... MUSCL reconstruction of momentum
!
! ... on East volume boundary
!
      lim = 0.D0
      erre = 0.D0
!
      cs = indyp * (u%c * dy(j+1) + u%n * dy(j))
      !cn = cs * dt * 2.D0 * indxp
      IF ( cs >= 0.D0 ) THEN
        fou  = dens_c * v%c
      ELSE IF ( cs < 0.D0 ) THEN
        fou  = dens_e * v%e
      END IF
!
      fe = fou * cs
!
! ... on North volume boundary
!
      cs = 0.5D0 * ( v%n + v%c )
      !cn = cs * dt * indy(j+1)
      IF (cs >= 0.D0) THEN
        fou  = dens_c * v%c
      ELSE IF (cs < 0.D0) THEN
        fou  = dens_n * v%n
      END IF
!
      fn = fou * cs
!
! ... on Top volume boundary
!
      cs = (dy(j+1) * w%c + dy(j) * w%n) * indyp
      !cn = cs * dt * 2.0 * indzp
      IF (cs >= 0.D0) THEN
        fou  = dens_c * v%c
      ELSE IF (cs < 0.D0) THEN
        fou  = dens_t * v%t
      END IF 
!
      ft = fou * cs
!
      RETURN
      END SUBROUTINE flv_3d_1st
!----------------------------------------------------------------------
      END MODULE convective_fluxes_v
!-----------------------------------------------------------------------
