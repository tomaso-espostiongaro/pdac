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
        MODULE PROCEDURE masf_2d, masf_3d, masf_3d_new
      END INTERFACE

!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE masf_3d( fe, fn, ft, fw, fs, fb, dens, u, v, w )
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

      SUBROUTINE masf_3d_new(fe, fn, ft, fw, fs, fb, dens, u, v, w, i, j, k, ijk)
!
! ... This routine computes convective mass fluxes by using
! ... Donor Cell technique for first order upwind
! ... This subroutine only if muscl == 0

      USE dimensions
      USE domain_decomposition, ONLY: myijk
      USE grid, ONLY: dx, dy, dz, fl_l
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: stencil
      USE time_parameters, ONLY: dt
      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: fe, fn, ft, fw, fs, fb
      TYPE(stencil), INTENT(IN) :: dens, u, v, w
      INTEGER, INTENT(IN) :: i, j, k, ijk
!
      INTEGER, SAVE :: ijk_old = -1
      REAL*8, SAVE :: indxp, indxm
      REAL*8, SAVE :: indyp, indym
      REAL*8, SAVE :: indzp, indzm
!
      IF( ijk /= ijk_old ) THEN
!
         indxm=1.D0/(dx(i)+dx(i-1))
         indxp=1.D0/(dx(i)+dx(i+1))
         indym=1.D0/(dy(j)+dy(j-1))
         indyp=1.D0/(dy(j)+dy(j+1))
         indzm=1.D0/(dz(k)+dz(k-1))
         indzp=1.D0/(dz(k)+dz(k+1))
     
         ijk_old = ijk
      END IF
!
! ... MUSCL reconstruction of fields
!
! ... on West volume boundary
!
      cs = u%w
      IF (cs >= 0.D0) THEN
        fou  = dens%w 
      ELSE IF (cs < 0.D0) THEN
        fou  = dens%c 
      ENDIF
!
      upwnd = fou 
!
      centrd = ( dx(i) * dens%w + dx(i-1) * dens%c ) * indxm
      upc_w = upwnd / centrd
!
      fw = upwnd * cs
!
! ... on East volume boundary
!
      cs = u%c
      IF (cs >= 0.D0) THEN
        fou  = dens%c 
      ELSE IF (cs < 0.D0) THEN
        fou  = dens%e 
      ENDIF
!
      upwnd = fou 
!
      centrd = ( dx(i) * dens%e + dx(i+1) * dens%c ) * indxp
      upc_e = upwnd / centrd
!
      fe = upwnd * cs
!
! ... on South volume boundary
!
!
      cs = v%s
      IF (cs >= 0.D0) THEN
        fou  = dens%s 
      ELSE IF (cs < 0.D0) THEN
        fou  = dens%c 
      ENDIF
!
      upwnd = fou
!
      centrd = ( dy(j) * dens%s + dy(j-1) * dens%c ) * indym
      upc_s = upwnd / centrd
!
      fs = upwnd * cs
!
! ... on North volume boundary
!
      cs = v%c
      IF (cs >= 0.D0) THEN
        fou  = dens%c 
      ELSE IF (cs < 0.D0) THEN
        fou  = dens%n 
      ENDIF
!
      upwnd = fou 
!
      centrd = ( dy(j) * dens%n + dy(j+1) * dens%c ) * indyp
      upc_n = upwnd / centrd
!
      fn = upwnd * cs
!
! ... on Bottom volume boundary
!
      cs = w%b
      IF (cs >= 0.D0) THEN
        fou  = dens%b 
      ELSE IF (cs < 0.D0) THEN
        fou  = dens%c 
      ENDIF
!
      upwnd = fou
!
      centrd = ( dz(k) * dens%b + dz(k-1) * dens%c ) * indzm
      upc_b = upwnd / centrd
!
      fb = upwnd * cs
!
! ... on Top volume boundary
!
      cs = w%c
      IF (cs >= 0.D0) THEN
        fou  = dens%c 
      ELSE IF (cs < 0.D0) THEN
        fou  = dens%t 
      ENDIF
!
      upwnd = fou
!
      centrd = ( dz(k) * dens%t + dz(k+1) * dens%c ) * indzp
      upc_t = upwnd / centrd
!
      ft = upwnd * cs
!
      RETURN
      END SUBROUTINE masf_3d_new


      SUBROUTINE fmas_3d(fe, fn, ft, fw, fs, fb, dens, u, v, w, i, j, k, ijk)
!
      USE dimensions
      USE domain_decomposition, ONLY: myijk
      USE grid, ONLY: dx, dy, dz, fl_l
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: stencil
      USE time_parameters, ONLY: dt
      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: fe, fn, ft, fw, fs, fb
      TYPE(stencil), INTENT(IN) :: dens, u, v, w
      INTEGER, INTENT(IN) :: i, j, k, ijk

      INTEGER :: im2, ip2, jm2, jp2, km2, kp2
!
      INTEGER, SAVE :: ijk_old = -1
      REAL*8, SAVE :: dxmm, dxm, dxp, dxpp, indxpp, indxp, indxm, indxmm
      REAL*8, SAVE :: dymm, dym, dyp, dypp, indypp, indyp, indym, indymm
      REAL*8, SAVE :: dzmm, dzm, dzp, dzpp, indzpp, indzp, indzm, indzmm
      REAL*8, SAVE :: gradc, grade, gradw, gradn, grads, gradt, gradb
!
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
      IF ((muscl /= 0) .AND. (gradc /= 0.D0) .AND. (i /= 2)) THEN
        CALL limiters(lim,erre)
      END IF
!
      upwnd = fou + lim * gradc * incr
!
      centrd = ( dx(i) * dens%w + dx(i-1) * dens%c ) * indxm
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
      centrd = ( dx(i) * dens%e + dx(i+1) * dens%c ) * indxp
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
      IF ((muscl /= 0) .AND. (gradc /= 0.D0) .AND. (j /= 2)) THEN
        CALL limiters(lim,erre)
      END IF
!
      upwnd = fou + lim * gradc * incr
!
      centrd = ( dy(j) * dens%s + dy(j-1) * dens%c ) * indym
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
      centrd = ( dy(j) * dens%n + dy(j+1) * dens%c ) * indyp
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
      IF ((muscl /= 0) .AND. (gradc /= 0.D0) .AND. (k /= 2)) THEN
        CALL limiters(lim,erre)
      END IF
!
      upwnd = fou + lim * gradc * incr
!
      centrd = ( dz(k) * dens%b + dz(k-1) * dens%c ) * indzm
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
      centrd = ( dz(k) * dens%t + dz(k+1) * dens%c ) * indzp
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
      USE domain_decomposition, ONLY: myijk
      USE grid, ONLY: xb
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
!
! ... West, Bottom fluxes
!
      IF (u%w >= 0.D0) THEN
        fw = u%w * dens%w * xb(i-1)
      ELSE
        fw = u%w * dens%c * xb(i-1)
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
        fe = u%c * dens%c * xb(i)
      ELSE
        fe = u%c * dens%e * xb(i)
      ENDIF

      IF (w%c >= 0.D0) THEN
        ft = w%c * dens%c
      ELSE
        ft = w%c * dens%t
      ENDIF
!
      upc_e = 1.D0
      upc_w = 1.D0
      upc_t = 1.D0
      upc_b = 1.D0

      RETURN
      END SUBROUTINE masf_2d
!----------------------------------------------------------------------
      SUBROUTINE fmas_2d(fe, ft, fw, fb, dens, u, w, ij)
!
      USE dimensions
      USE domain_decomposition, ONLY: myijk
      USE grid, ONLY: dx, xb, dz, fl_l
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: stencil
      USE time_parameters, ONLY: dt
      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: fe, ft, fw, fb
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
      upc_e = 1.D0
      upc_w = 1.D0
      upc_t = 1.D0
      upc_b = 1.D0
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
      cn = cs * dt * 2.D0 * indxm
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
      IF ((muscl /= 0) .AND. (gradc /= 0.D0) .AND. (i /= 2)) THEN
        CALL limiters(lim,erre)
      END IF
!
      upwnd = fou + lim * gradc * incr
!
      centrd = (dx(i)*dens%w+dx(i-1)*dens%c)*indxm
      upc_w = upwnd / centrd
!
      fw = upwnd * cs * xb(i-1)
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
      fe = upwnd * cs * xb(i)
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
	incr = 0.5D0 * dz(k-1)
      ELSE IF (cs < 0.D0) THEN
	erre = gradt / gradc
        fou  = dens%c 
	incr = 0.5D0 * dz(k)
      ENDIF
!
      IF ((muscl /= 0) .AND. (gradc /= 0.D0) .AND. (k /= 2)) THEN
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
      END SUBROUTINE fmas_2d
!-----------------------------------------------------------------------
      END MODULE convective_mass_fluxes
!-----------------------------------------------------------------------
