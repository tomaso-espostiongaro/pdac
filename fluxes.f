!----------------------------------------------------------------------
      MODULE eulerian_flux
!
! ...  this module computes eulerian fluxes by using 
! ...  1st order upwind, 
! ...  Standard compass notation is adopted to locate cell neighbours:
! ...  n=north, s=south, e=east, w=west
! ...  ne=north-east, nw=north-west, sw=south-west, se=south-east, etc.
! ...  The computational stencil is defined in the `set_indexes' module
!
!----------------------------------------------------------------------
      USE grid, ONLY: dz, dr, rb, r, inr, inrb, indr, indz
      USE set_indexes
      USE time_parameters, ONLY: dt
      IMPLICIT NONE
!
      REAL*8, PRIVATE :: cs                        ! convective stream   !
      REAL*8, PRIVATE :: cn                        ! Courant number      !
      REAL*8, PRIVATE :: upwnd                     ! upwinded variable   !
!
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE fu_lb(fl, fb, dens, u, v, ij)
!
! ... Compute the convective fluxes on left and bottom element sides
! ... for the r(x)-momentum density.
!
      USE dimensions
      USE grid, ONLY: fl_l
      USE grid, ONLY: myij

      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: fl, fb
!
      TYPE(stencil), INTENT(IN) :: dens, u, v
      INTEGER, INTENT(IN) :: ij
      INTEGER :: i, ij_g
      INTEGER :: imj, ijm
!
      REAL*8 :: dens_e, dens_w, dens_se
      REAL*8 :: drm, drp, indrp, indrm
!
      ij_g = myij( 0, 0, ij)
      i  = MOD( ( ij_g - 1 ), ndi) + 1
!
      drp=dr(i)+dr(i+1)
      drm=dr(i)+dr(i-1)
      indrp=1.D0/drp
      indrm=1.D0/drm
!
! ... Compute linearly interpolated values of density 
! ... where not available.
!
      dens_e = (dr(i+1) * dens%c + dr(i) * dens%e) * indrp
      dens_w = (dr(i-1) * dens%c + dr(i) * dens%w) * indrm
      dens_se =(dr(i+1) * dens%s + dr(i) * dens%se) * indrp
!
! ... On boundary mantain first order accuracy (1st order Upwind).
!
      imj = myij(-1, 0, ij)
      IF( fl_l(imj) .NE. 1 ) THEN
        cs=0.5D0*(u%c + u%w)
        IF(cs.GE.0.D0) fl = dens_w * u%w * cs * r(i)
        IF(cs.LT.0.D0) fl = dens_e * u%c * cs * r(i)
      END IF

      ijm = myij( 0,-1, ij)
      IF( fl_l(ijm) .NE. 1 ) THEN
        cs=(dr(i+1) * v%s + dr(i) * v%se) * indrp
        IF(cs.GE.0.D0) fb = dens_se * u%s * cs
        IF(cs.LT.0.D0) fb = dens_e * u%c * cs
      END IF
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE fu_rt(fr, ft, dens, u, v, ij)
!
! ... Compute the convective fluxes on right and top element sides
! ... for the r(x)-momentum density.
!
      USE dimensions
      USE grid, ONLY: myij

      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: fr, ft
      TYPE(stencil), INTENT(IN) :: dens, u, v
      INTEGER, INTENT(IN) :: ij
      INTEGER :: i,j,ij_g
!
      REAL*8 :: dens_e, dens_ee, dens_ne
      REAL*8 :: drp, drpp, indrpp, indrp
      REAL*8 :: dzp, indzp
!
      ij_g = myij( 0, 0, ij)
      i  = MOD( ( ij_g - 1 ), ndi) + 1
      j  = ( ij_g - 1 ) / ndi + 1
!
      drp=dr(i)+dr(i+1)
      drpp=dr(i+1)+dr(i+2)
      dzp=dz(j)+dz(j+1)
      indrp=1.D0/drp
      indrpp=1.D0/drpp
      indzp=1.D0/dzp
!       
! ... Values of density interpolated at cell boundaries
!
      dens_e = (dr(i+1) * dens%c + dr(i) * dens%e) * indrp
      dens_ne = (dr(i+1) * dens%n + dr(i) * dens%ne) * indrp
      dens_ee = (dr(i+2) * dens%e + dr(i+1) * dens%ee) * indrpp
!
! ... on right volume boundary
!
      cs = 0.5D0 * (u%c + u%e)
      cn = cs * dt * indr(i+1)
      IF (cs.GE.0.D0) THEN
        upwnd = dens_e * u%c
      ELSE IF (cs.LT.0.D0) THEN
        upwnd = dens_ee * u%e
      END IF
!
      fr = upwnd * cs * r(i+1)
!
! ... on top volume boundary
!
      cs=(dr(i+1) * v%c + dr(i) * v%e) * indrp
      cn=cs * dt * 2.0 * indzp
      IF (cs.GE.0.D0) THEN
        upwnd = dens_e * u%c
      ELSE IF (cs.LT.0.D0) THEN
        upwnd = dens_ne * u%n
      END IF 
!
      ft = upwnd * cs
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE fv_lb(fl, fb, dens, u, v, ij)
!
      USE dimensions
      USE grid, ONLY: fl_l
      USE grid, ONLY: myij
      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: fl, fb
      TYPE(stencil), INTENT(IN) :: dens, u, v
      INTEGER, INTENT(IN) :: ij
!
      INTEGER :: i,j, ij_g
      INTEGER :: ijm, imj
      REAL*8 :: dzm, dzp, indzm, indzp
      REAL*8 :: dens_n, dens_s, dens_wn
!
      ij_g = myij( 0, 0, ij)
      j  = ( ij_g - 1 ) / ndi + 1
      i  = MOD( ( ij_g - 1 ), ndi) + 1
!
! ... Compute linearly interpolated values of density
! ... where not available.
!
      dens_n = (dz(j+1) * dens%c  + dz(j) * dens%n) * indzp
      dens_s = (dz(j-1) * dens%c  + dz(j) * dens%s) * indzm
      dens_wn = (dz(j+1) * dens%w  + dz(j) * dens%nw) * indzp
!
      dzp=dz(j)+dz(j+1)
      dzm=dz(j)+dz(j-1)
      indzp=1.D0/dzp
      indzm=1.D0/dzm
!
! ... On boundary mantain first order accuracy (1st order Upwind).
!
      imj = myij(-1, 0, ij)
      IF(fl_l(imj).NE.1) THEN
        cs=(dz(j+1)*u%w+dz(j)*u%nw)*indzp
        IF(cs.GE.0.D0) fl = dens_wn * v%w * cs * rb(i-1)
        IF(cs.LT.0.D0) fl = dens_n * v%c * cs * rb(i-1)
      END IF
!
      ijm = myij( 0,-1, ij)
      IF(fl_l(ijm).NE.1) THEN
        cs=0.5D0*(v%s+v%c)
        IF(cs.GE.0.D0) fb = dens_s * v%s * cs
        IF(cs.LT.0.D0) fb = dens_n * v%c * cs
      END IF
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE fv_rt(fr, ft, dens, u, v, ij)
!
      USE dimensions
      USE grid, ONLY: fl_l
      USE grid, ONLY: myij
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: ij
      REAL*8, INTENT(OUT) :: fr, ft
      TYPE(stencil), INTENT(IN) :: dens, u, v
!
      INTEGER :: i,j, ij_g
      REAL*8 :: dzp, dzpp, indzpp, indzp
      REAL*8 :: dens_n, dens_nn, dens_en
!
      ij_g = myij( 0, 0, ij)
      j  = ( ij_g - 1 ) / ndi + 1
      i  = MOD( ( ij_g - 1 ), ndi) + 1
!
! ... values of density interpolated at cell boundaries
        dens_n = (dz(j+1) * dens%c + dz(j) * dens%n) * indzp
        dens_en = (dz(j+1) * dens%e + dz(j) * dens%ne) * indzp
        dens_nn = (dz(j+2) * dens%n + dz(j+1) * dens%nn) * indzpp
!
        dzp=dz(j)+dz(j+1)
        dzpp=dz(j+1)+dz(j+2)
        indzp=1.D0/dzp
        indzpp=1.D0/dzpp
!
! ... on right volume boundary
!
      cs = (dz(j+1)*u%c+dz(j)*u%n)*indzp
      cn = cs * dt * 2.0 * indzp
      IF (cs.GE.0.D0) THEN
        upwnd = dens_n * v%c
      ELSE IF (cs.LT.0.D0) THEN
        upwnd = dens_en * v%e
      END IF
!
      fr = upwnd * cs * rb(i)
!
! ... on top volume boundary
!
      cs = 0.5D0*(v%c+v%n)
      cn = cs * dt * indz(j+1)
      IF (cs.GE.0.D0) THEN
        upwnd = dens_n * v%c
      ELSE IF (cs.LT.0.D0 .AND. j.NE.(ndj-1)) THEN
        upwnd = dens_nn * v%n
      END IF 
!
      ft = upwnd * cs
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE fsc_lb(fl, fb, dens, field, u, v, ij)
!
      USE dimensions
      USE grid, ONLY: fl_l
      USE grid, ONLY: myij
      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: fl, fb
      TYPE(stencil), INTENT(IN) :: dens, field, u, v
      INTEGER, INTENT(IN) :: ij
!
      INTEGER :: i,j, ij_g
      INTEGER :: imj, ijm
      REAL*8  :: drp, drm, drmm, dzp, dzm, dzmm
      REAL*8  :: indrp, indrm, indrmm, indzp, indzm, indzmm
!
      ij_g = myij( 0, 0, ij)
      imj = myij(-1, 0, ij)
      ijm = myij( 0,-1, ij)
      i  = MOD( ( ij_g - 1 ), ndi) + 1
      j  = ( ij_g - 1 ) / ndi + 1
!
! ... on left volume boundary
!
      IF ((fl_l(imj) .NE. 1)) THEN
        cs = u%w
        IF (cs .GE. 0.D0) THEN
          upwnd = dens%w*field%w
        ELSE IF (cs .LT. 0.D0) THEN
          upwnd = dens%c*field%c
        ENDIF
        fl = upwnd * cs * rb(i-1)
      END IF
!
! ... on bottom volume boundary
!
      IF ((fl_l(ijm) .NE. 1)) THEN
        cs = v%s
        IF (cs .GE. 0.D0) THEN
          upwnd = dens%s*field%s
        ELSE IF (cs .LT. 0.D0) THEN
          upwnd = dens%c*field%c
        ENDIF
        fb = upwnd * cs
      END IF
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE fsc_rt(fr, ft, dens, field, u, v, ij)
!
      USE dimensions
      USE grid, ONLY: myij
      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: fr, ft
      TYPE(stencil), INTENT(IN) :: dens, field, u, v
      INTEGER, INTENT(IN) :: ij
!
      INTEGER :: i,j,ij_g
      REAL*8  :: drp, dzp
      REAL*8  :: indrp, indzp 
!
      ij_g = myij( 0, 0, ij)
      i  = MOD( ( ij_g - 1 ), ndi) + 1
      j  = ( ij_g - 1 ) / ndi + 1
!
        drp = (dr(i) + dr(i+1))
        dzp = (dz(j) + dz(j+1))
        indrp = 1.D0/drp
        indzp = 1.D0/dzp
!
! ... on right volume boundary
!
        cs = u%c
        cn = cs * dt * 2.0 * indrp
        IF (cs .GE. 0.D0) THEN
          upwnd  = dens%c*field%c
        ELSE IF (cs .LT. 0.D0) THEN
          upwnd  = dens%e*field%e
        ENDIF
!
      fr = upwnd * cs * rb(i)
!
! ... on top volume boundary
!
        cs = v%c
        cn = cs * dt * 2.0 * indzp
        IF (cs .GE. 0.D0) THEN
          upwnd = dens%c*field%c
        ELSE IF (cs .LT. 0.D0) THEN
          upwnd = dens%n*field%n
        ENDIF
!
      ft = upwnd * cs
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE masfg(rgfrm, rgftm, rgfr, rgft, ug, vg, rgp, i)
!
      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: rgfr, rgft, rgfrm, rgftm
      TYPE(stencil), INTENT(IN) :: ug, vg, rgp
      INTEGER, INTENT(IN) :: i
!
      IF (ug%w.GE.0.D0) THEN
        rgfrm = ug%w * rgp%w * rb(i-1)
      ELSE
        rgfrm = ug%w * rgp%c * rb(i-1)
      ENDIF
      IF (vg%s.GE.0.D0) THEN
        rgftm = vg%s * rgp%s
      ELSE
        rgftm = vg%s * rgp%c
      ENDIF
!
      IF (ug%c.GE.0.D0) THEN
        rgfr = ug%c * rgp%c * rb(i)
      ELSE
        rgfr = ug%c * rgp%e * rb(i)
      ENDIF
      IF (vg%c.GE.0.D0) THEN
        rgft = vg%c * rgp%c
      ELSE
        rgft = vg%c * rgp%n
      ENDIF
!
      RETURN
      END SUBROUTINE
!
!----------------------------------------------------------------------
      SUBROUTINE masfk(rlfrm, rlftm, rlfr, rlft, uk, vk, rlk, i)
!
      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: rlfrm, rlftm, rlfr, rlft
      TYPE(stencil), INTENT(IN) :: uk, vk, rlk
      INTEGER, INTENT(IN) :: i
!
        IF (uk%w.GE.0.D0) THEN
          rlfrm = uk%w * rlk%w * rb(i-1)
        ELSE
          rlfrm = uk%w * rlk%c * rb(i-1)
        END IF
        IF (vk%s.GE.0.D0) THEN
          rlftm = vk%s * rlk%s
        ELSE
          rlftm = vk%s * rlk%c
        END IF
! 
        IF (uk%c.GE.0.D0) THEN
          rlfr = uk%c * rlk%c * rb(i)
        ELSE
          rlfr = uk%c * rlk%e * rb(i)
        END IF
        IF (vk%c.GE.0.D0) THEN
          rlft = vk%c * rlk%c
        ELSE
          rlft = vk%c * rlk%n
        END IF
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      END MODULE eulerian_flux
!-------------------------------------------------------------
