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
      IF (u%w >= 0.D0) THEN
        fw = u%w * dens%w
      ELSE
        fw = u%w * dens%c
      ENDIF

      IF (v%s >= 0.D0) THEN
        fs = v%s * dens%s
      ELSE
        fs = v%s * dens%c
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
      USE domain_decomposition, ONLY: myijk, meshinds
      USE grid, ONLY: dx, dy, dz, flag
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: stencil
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
         im2 = MAX( 1, i-2 )
         ip2 = MIN( nx, i+2 )
         jm2 = MAX( 1, j-2 )
         jp2 = MIN( ny, j+2 )
         km2 = MAX( 1, k-2 )
         kp2 = MIN( nz, k+2 )

         dxmm=dx(i-1)+dx(im2)
         dxm=dx(i)+dx(i-1)
         dxp=dx(i)+dx(i+1)
         dxpp=dx(i+1)+dx(ip2)
         dymm=dy(j-1)+dy(jm2)
         dym=dy(j)+dy(j-1)
         dyp=dy(j)+dy(j+1)
         dypp=dy(j+1)+dy(jp2)
         dzmm=dz(k-1)+dz(km2)
         dzm=dz(k)+dz(k-1)
         dzp=dz(k)+dz(k+1)
         dzpp=dz(k+1)+dz(kp2)

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
     
         ijk_old = ijk
      END IF
!
! ... MUSCL reconstruction of fields
!
! ... on West volume boundary
!
      gradc = 2.D0 * indxm * (dens%c - dens%w)
      gradw = 2.D0 * indxmm * (dens%w - dens%ww)
      grade = 2.D0 * indxp * (dens%e - dens%c)
!
      lim = 0.D0
      erre = 0.D0
!
      cs = u%w
      !cn = cs * dt * 2.D0 * indxm
      IF (cs >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = gradw / gradc
	incr = 0.5D0 * dx(i-1)
      ELSE IF (cs < 0.D0) THEN
	IF (gradc /= 0.D0) erre = grade / gradc
	incr = 0.5D0 * dx(i)
      ENDIF
!
      IF (i/=2) CALL limiters(lim,erre)
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
      !cn = cs * dt * 2.D0 * indxp
      IF (cs >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = gradw / gradc
	incr = 0.5D0 * dx(i)
      ELSE IF (cs < 0.D0) THEN
	IF (gradc /= 0.D0) erre = grade / gradc
	incr = 0.5D0 * dx(i+1)
      ENDIF
!
      IF (i/=nx-1) CALL limiters(lim,erre)
!
      upwnd = lim * gradc * incr
      centrd = ( dx(i) * dens%e + dx(i+1) * dens%c ) * indxp

      fe = fe + upwnd * cs

      upc_e = fe / cs / centrd
!
! ... on South volume boundary
!
      gradc = 2.D0 * indym * (dens%c - dens%s)
      grads = 2.D0 * indymm * (dens%s - dens%ss)
      gradn = 2.D0 * indyp * (dens%n - dens%c)
!
      lim = 0.D0
      erre = 0.D0
!
      cs = v%s
      !cn = cs * dt * 2.D0 * indym
      IF (cs >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = grads / gradc
	incr = 0.5D0 * dy(j-1)
      ELSE IF (cs < 0.D0) THEN
	IF (gradc /= 0.D0) erre = gradn / gradc
	incr = 0.5D0 * dy(j)
      ENDIF
!
      IF (j/=2) CALL limiters(lim,erre)
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
      !cn = cs * dt * 2.D0 * indyp
      IF (cs >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = grads / gradc
	incr = 0.5D0 * dy(j)
      ELSE IF (cs < 0.D0) THEN
	IF (gradc /= 0.D0) erre = gradn / gradc
	incr = 0.5D0 * dy(j+1)
      ENDIF
!
      IF (j/=ny-1) CALL limiters(lim,erre)
!
      upwnd = lim * gradc * incr
      centrd = ( dy(j) * dens%n + dy(j+1) * dens%c ) * indyp

      fn = fn + upwnd * cs

      upc_n = fn / cs / centrd
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
      !cn = cs * dt * 2.D0 * indzm
      IF (cs >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = gradb / gradc
	incr = 0.5D0 * dz(k-1)
      ELSE IF (cs < 0.D0) THEN
	IF (gradc /= 0.D0) erre = gradt / gradc
	incr = 0.5D0 * dz(k)
      ENDIF
!
      IF (k/=2) CALL limiters(lim,erre)
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
      !cn = cs * dt * 2.D0 * indzp
      IF (cs >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = gradb / gradc
	incr = 0.5D0 * dz(k)
      ELSE IF (cs < 0.D0) THEN
	IF (gradc /= 0.D0) erre = gradt / gradc
	incr = 0.5D0 * dz(k+1)
      ENDIF
!
      IF (k/=nz-1) CALL limiters(lim,erre)
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
      USE domain_decomposition, ONLY: myijk
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
      USE domain_decomposition, ONLY: myijk
      USE grid, ONLY: dx, rb, dz, flag
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: stencil
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

      im2 = MAX( 1, i-2 )
      ip2 = MIN( nx, i+2 )
      km2 = MAX( 1, k-2 )
      kp2 = MIN( nz, k+2 )
!
      dxmm=dx(i-1)+dx(im2)
      dxm=dx(i)+dx(i-1)
      dxp=dx(i)+dx(i+1)
      dxpp=dx(i+1)+dx(ip2)
      dzmm=dz(k-1)+dz(km2)
      dzm=dz(k)+dz(k-1)
      dzp=dz(k)+dz(k+1)
      dzpp=dz(k+1)+dz(kp2)

      indxmm=1.D0/dxmm
      indxm=1.D0/dxm
      indxp=1.D0/dxp
      indxpp=1.D0/dxpp
      indzmm=1.D0/dzmm
      indzm=1.D0/dzm
      indzp=1.D0/dzp
      indzpp=1.D0/dzpp
!
! ... MUSCL reconstruction of fields
!
! ... on West volume boundary
!
      gradc = 2.D0 * indxm * (dens%c - dens%w)
      gradw = 2.D0 * indxmm * (dens%w - dens%ww)
      grade = 2.D0 * indxp * (dens%e - dens%c)
!
      lim = 0.D0
      erre = 0.D0
!
      cs = u%w
      !cn = cs * dt * 2.D0 * indxm
      IF (cs >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = gradw / gradc
	incr = 0.5D0 * dx(i-1)
      ELSE IF (cs < 0.D0) THEN
	IF (gradc /= 0.D0) erre = grade / gradc
	incr = 0.5D0 * dx(i)
      ENDIF
!
      IF (i/=2) CALL limiters(lim,erre)
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
      !cn = cs * dt * 2.D0 * indxp
      IF (cs >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = gradw / gradc
	incr = 0.5D0 * dx(i)
      ELSE IF (cs < 0.D0) THEN
	IF (gradc /= 0.D0) erre = grade / gradc
	incr = 0.5D0 * dx(i+1)
      ENDIF
!
      IF (i/=nx-1) CALL limiters(lim,erre)
!
      upwnd = lim * gradc * incr
      centrd = (dx(i)*dens%e+dx(i+1)*dens%c)*indxp
      
      fe = fe + upwnd * cs * rb(i)

      upc_e = fe / (cs * rb(i)) / centrd
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
      !cn = cs * dt * 2.D0 * indzm
      IF (cs >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = gradb / gradc
	incr = 0.5D0 * dz(k-1)
      ELSE IF (cs < 0.D0) THEN
	IF (gradc /= 0.D0) erre = gradt / gradc
	incr = 0.5D0 * dz(k)
      ENDIF
!
      IF (k/=2) CALL limiters(lim,erre)
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
      !cn = cs * dt * 2.D0 * indzp
      IF (cs >= 0.D0) THEN
	IF (gradc /= 0.D0) erre = gradb / gradc
	incr = 0.5D0 * dz(k)
      ELSE IF (cs < 0.D0) THEN
	IF (gradc /= 0.D0) erre = gradt / gradc
	incr = 0.5D0 * dz(k+1)
      ENDIF
!
      IF (k/=nz-1) CALL limiters(lim,erre)
!
      upwnd = lim * gradc * incr
      centrd = (dz(k)*dens%t+dz(k+1)*dens%c)*indzp

      ft = ft + upwnd * cs

      upc_t = ft / cs / centrd
!
      RETURN
      END SUBROUTINE fmas_2d
!-----------------------------------------------------------------------
      END MODULE convective_mass_fluxes
!-----------------------------------------------------------------------
