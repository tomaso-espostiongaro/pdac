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

      REAL*8 :: beta, muscl
!
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE fu_lb(fl, fb, dens, u, w, ij)
!
! ... Compute the convective fluxes on left and bottom sides of the cell
! ... for the r(x)-momentum density.
!
      USE dimensions
      USE grid, ONLY: fl_l
      USE grid, ONLY: myij

      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: fl, fb
!
      TYPE(stencil), INTENT(IN) :: dens, u, w
      INTEGER, INTENT(IN) :: ij
      INTEGER :: i, imesh
      INTEGER :: imj, ijm
!
      REAL*8 :: dens_e, dens_w, dens_se
      REAL*8 :: drm, drp, indrp, indrm
!
      imesh = myij( 0, 0, ij)
      i  = MOD( ( imesh - 1 ), nr) + 1
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
        cs=(dr(i+1) * w%s + dr(i) * w%se) * indrp
        IF(cs.GE.0.D0) fb = dens_se * u%s * cs
        IF(cs.LT.0.D0) fb = dens_e * u%c * cs
      END IF
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE fu_rt(fr, ft, dens, u, w, ij)
!
! ... Compute the convective fluxes on right and top  sides of the cell
! ... for the r(x)-momentum density.
!
      USE dimensions
      USE grid, ONLY: myij

      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: fr, ft
      TYPE(stencil), INTENT(IN) :: dens, u, w
      INTEGER, INTENT(IN) :: ij
      INTEGER :: i,j,imesh
!
      REAL*8 :: dens_e, dens_ee, dens_ne
      REAL*8 :: drm, drp, drpp, indrpp, indrp, indrm
      REAL*8 :: dzp, indzp, dzm, indzm, dzpp, indzpp
      REAL*8 :: gradc, grade, gradn, gradw, grads, gradgt, gradlt
!
      imesh = myij( 0, 0, ij)
      i  = MOD( ( imesh - 1 ), nr) + 1
      j  = ( imesh - 1 ) / nr + 1
!
      drm=dr(i)+dr(i-1)
      drp=dr(i)+dr(i+1)
      drpp=dr(i+1)+dr(i+2)
      dzm=dz(j)+dz(j-1)
      dzp=dz(j)+dz(j+1)
      dzpp=dz(j+1)+dz(j+2)

      indrm=1.D0/drm
      indrp=1.D0/drp
      indrpp=1.D0/drpp
      indzm=1.D0/dzm
      indzp=1.D0/dzp
      indzpp=1.D0/dzpp
!       
! ... Values of density interpolated at cell boundaries
!
      dens_e = (dr(i+1) * dens%c + dr(i) * dens%e) * indrp
      dens_ne = (dr(i+1) * dens%n + dr(i) * dens%ne) * indrp
      dens_ee = (dr(i+2) * dens%e + dr(i+1) * dens%ee) * indrpp
!
! ... on right volume boundary
!
! ... MUSCL reconstruction of velocity (density unvaried)
!
      gradc = (indr(i+1) * (u%e - u%c))
      gradw = (indr(i) * (u%c - u%w))
      grade = (indr(i+2) * (u%ee - u%e))
      
      gradgt = (1.D0 - beta) * gradc + beta * gradw
      gradlt= (1.D0 - beta) * gradc + beta * grade

      cs = 0.5D0 * (u%c + u%e)
      cn = cs * dt * indr(i+1)
      IF (cs.GE.0.D0) THEN
        upwnd = dens_e * (u%c + muscl*(gradgt)*0.5*dr(i+1))
      ELSE IF (cs.LT.0.D0) THEN
        upwnd = dens_ee * (u%e + muscl*(gradlt)*0.5*dr(i+1))
      END IF
!
      fr = upwnd * cs * r(i+1)
!
! ... on top volume boundary
!
      gradc = (u%n - u%c) * 2.0 * indzp
      gradn = (u%nn - u%n) * 2.0 * indzpp
      grads = (u%c - u%s) * 2.0 * indzm

      gradgt = (1.D0 - beta) * gradc + beta * grads
      gradlt = (1.D0 - beta) * gradc + beta * gradn

      cs=(dr(i+1) * w%c + dr(i) * w%e) * indrp
      cn=cs * dt * 2.0 * indzp
      IF (cs.GE.0.D0) THEN
        upwnd = dens_e * (u%c + muscl*(gradgt)*0.5*dz(j))
      ELSE IF (cs.LT.0.D0) THEN
        upwnd = dens_ne * (u%n + muscl*(gradlt)*0.5*dz(j+1))
      END IF 
!
      ft = upwnd * cs
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE fw_lb(fl, fb, dens, u, w, ij)
!
      USE dimensions
      USE grid, ONLY: fl_l
      USE grid, ONLY: myij
      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: fl, fb
      TYPE(stencil), INTENT(IN) :: dens, u, w
      INTEGER, INTENT(IN) :: ij
!
      INTEGER :: i,j, imesh
      INTEGER :: ijm, imj
      REAL*8 :: dzm, dzp, indzm, indzp
      REAL*8 :: dens_n, dens_s, dens_wn
!
      imesh = myij( 0, 0, ij)
      j  = ( imesh - 1 ) / nr + 1
      i  = MOD( ( imesh - 1 ), nr) + 1
!
      dzp=dz(j)+dz(j+1)
      dzm=dz(j)+dz(j-1)
      indzp=1.D0/dzp
      indzm=1.D0/dzm
!
! ... Compute linearly interpolated values of density
! ... where not available.
!
      dens_n = (dz(j+1) * dens%c  + dz(j) * dens%n) * indzp
      dens_s = (dz(j-1) * dens%c  + dz(j) * dens%s) * indzm
      dens_wn = (dz(j+1) * dens%w  + dz(j) * dens%nw) * indzp
!
! ... On boundary mantain first order accuracy (1st order Upwind).
!
      imj = myij(-1, 0, ij)
      IF(fl_l(imj).NE.1) THEN
        cs=(dz(j+1)*u%w+dz(j)*u%nw)*indzp
        IF(cs.GE.0.D0) fl = dens_wn * w%w * cs * rb(i-1)
        IF(cs.LT.0.D0) fl = dens_n * w%c * cs * rb(i-1)
      END IF
!
      ijm = myij( 0,-1, ij)
      IF(fl_l(ijm).NE.1) THEN
        cs=0.5D0*(w%s+w%c)
        IF(cs.GE.0.D0) fb = dens_s * w%s * cs
        IF(cs.LT.0.D0) fb = dens_n * w%c * cs
      END IF
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE fw_rt(fr, ft, dens, u, w, ij)
!
      USE dimensions
      USE grid, ONLY: fl_l
      USE grid, ONLY: myij
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: ij
      REAL*8, INTENT(OUT) :: fr, ft
      TYPE(stencil), INTENT(IN) :: dens, u, w
!
      INTEGER :: i,j, imesh
      REAL*8 :: dzp, dzpp, indzpp, indzp
      REAL*8 :: drm, drp, drpp, indrm, indrp, indrpp
      REAL*8 :: dens_n, dens_nn, dens_en
      REAL*8 :: gradc, grade, gradn, gradw, grads, gradgt, gradlt
!
      imesh = myij( 0, 0, ij)
      j  = ( imesh - 1 ) / nr + 1
      i  = MOD( ( imesh - 1 ), nr) + 1
!
      dzp=dz(j)+dz(j+1)
      dzpp=dz(j+1)+dz(j+2)
      drm=dr(i)+dr(i-1)
      drp=dr(i)+dr(i+1)
      drpp=dr(i+1)+dr(i+2)

      indzp=1.D0/dzp
      indzpp=1.D0/dzpp
      indrm=1.D0/drm
      indrp=1.D0/drp
      indrpp=1.D0/drpp
!
! ... values of density interpolated at cell boundaries
      dens_n = (dz(j+1) * dens%c + dz(j) * dens%n) * indzp
      dens_en = (dz(j+1) * dens%e + dz(j) * dens%ne) * indzp
      dens_nn = (dz(j+2) * dens%n + dz(j+1) * dens%nn) * indzpp
!
! ... on right volume boundary
!
! ... MUSCL reconstruction of velocity (density unvaried)
!
      gradc = (2.0 * indrp * (w%e - w%c))
      gradw = (2.0 * indrm * (w%c - w%w))
      grade = (2.0 * indrpp * (w%ee - w%e))

      gradgt = (1.0 - beta) * gradc + beta * gradw
      gradlt = (1.0 - beta) * gradc + beta * grade

      cs = (dz(j+1)*u%c+dz(j)*u%n)*indzp
      cn = cs * dt * 2.0 * indzp
      IF (cs.GE.0.D0) THEN
        upwnd = dens_n * (w%c + muscl*(gradgt)*0.5*dr(i))
      ELSE IF (cs.LT.0.D0) THEN
        upwnd = dens_en * (w%e + muscl*(gradlt)*0.5*dr(i+1))
      END IF
!
      fr = upwnd * cs * rb(i)
!
! ... on top volume boundary
!
      gradc = (indz(j+1) * (w%n - w%c))
      gradn = (indz(j+2) * (w%nn - w%n))
      grads = (indz(j) * (w%c - w%s))

      gradgt = (1.0 - beta) * gradc + beta * grads
      gradlt = (1.0 - beta) * gradc + beta * gradn

      cs = 0.5D0*(w%c+w%n)
      cn = cs * dt * indz(j+1)
      IF (cs.GE.0.D0) THEN
        upwnd = dens_n * (w%c + muscl*(gradgt)*0.5*dz(j+1))
      ELSE IF (cs.LT.0.D0 .AND. j.NE.(nz-1)) THEN
        upwnd = dens_nn * (w%n + muscl*(gradlt)*0.5*dz(j+1))
      END IF 
!
      ft = upwnd * cs
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE fsc_lb(fl, fb, dens, field, u, w, ij)
!
      USE dimensions
      USE grid, ONLY: fl_l
      USE grid, ONLY: myij
      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: fl, fb
      TYPE(stencil), INTENT(IN) :: dens, field, u, w
      INTEGER, INTENT(IN) :: ij
!
      INTEGER :: i,j, imesh
      INTEGER :: imj, ijm
      REAL*8  :: drp, drm, drmm, dzp, dzm, dzmm
      REAL*8  :: indrp, indrm, indrmm, indzp, indzm, indzmm
!
      imesh = myij( 0, 0, ij)
      imj = myij(-1, 0, ij)
      ijm = myij( 0,-1, ij)
      i  = MOD( ( imesh - 1 ), nr) + 1
      j  = ( imesh - 1 ) / nr + 1
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
        cs = w%s
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
      SUBROUTINE fsc_rt(fr, ft, dens, field, u, w, ij)
!
      USE dimensions
      USE grid, ONLY: myij
      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: fr, ft
      TYPE(stencil), INTENT(IN) :: dens, field, u, w
      INTEGER, INTENT(IN) :: ij
!
      INTEGER :: i,j,imesh
      REAL*8  :: drm, drp, dzm, dzp, drpp, dzpp
      REAL*8  :: indrm, indrp, indzm, indzp , indrpp, indzpp
      REAL*8 :: gradc, grade, gradn, gradw, grads, gradgt, gradlt
!
      imesh = myij( 0, 0, ij)
      i  = MOD( ( imesh - 1 ), nr) + 1
      j  = ( imesh - 1 ) / nr + 1
!
        drm = (dr(i) + dr(i-1))
        drp = (dr(i) + dr(i+1))
        drpp = (dr(i+1) + dr(i+2))
        dzm = (dz(j) + dz(j-1))
        dzp = (dz(j) + dz(j+1))
        dzpp = (dz(j+1) + dz(j+2))

        indrm = 1.D0/drm
        indrp = 1.D0/drp
        indrpp = 1.D0/drpp
        indzm = 1.D0/dzm
        indzp = 1.D0/dzp
        indzpp = 1.D0/dzpp
!
! ... on right volume boundary
!
! ... MUSCL reconstruction of velocity (density unvaried)
!
        gradc = 2.0 * indrp * (dens%e*field%e - dens%c*field%c)
        gradw = 2.0 * indrm * (dens%c*field%c - dens%w*field%w)
        grade = 2.0 * indrpp * (dens%ee*field%ee - dens%e*field%e)

        gradgt = (1.D0 - beta) * gradc + beta * gradw
        gradlt = (1.D0 - beta) * gradc + beta * grade

        cs = u%c
        cn = cs * dt * 2.0 * indrp
        IF (cs .GE. 0.D0) THEN
          upwnd  = dens%c*field%c + muscl*(gradgt)*0.5*dr(i)
        ELSE IF (cs .LT. 0.D0) THEN
          upwnd  = dens%e*field%e + muscl*(gradlt)*0.5*dr(i+1)
        ENDIF
!
      fr = upwnd * cs * rb(i)
!
! ... on top volume boundary
!
        gradc = 2.0 * indzp * (dens%n*field%n - dens%c*field%c)
        grads = 2.0 * indzm * (dens%c*field%c - dens%s*field%s)
        gradn = 2.0 * indzpp * (dens%nn*field%nn - dens%n*field%n)

        gradgt = (1.D0 - beta) * gradc + beta * grads
        gradlt = (1.D0 - beta) * gradc + beta * gradn

        cs = w%c
        cn = cs * dt * 2.0 * indzp
        IF (cs .GE. 0.D0) THEN
          upwnd = dens%c*field%c + muscl*(gradgt)*0.5*dz(j)
        ELSE IF (cs .LT. 0.D0) THEN
          upwnd = dens%n*field%n + muscl*(gradlt)*0.5*dz(j+1)
        ENDIF
!
      ft = upwnd * cs
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE masfg(rgfrm, rgftm, rgfr, rgft, ug, wg, rgp, i)
!
      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: rgfr, rgft, rgfrm, rgftm
      TYPE(stencil), INTENT(IN) :: ug, wg, rgp
      INTEGER, INTENT(IN) :: i
!
      IF (ug%w.GE.0.D0) THEN
        rgfrm = ug%w * rgp%w * rb(i-1)
      ELSE
        rgfrm = ug%w * rgp%c * rb(i-1)
      ENDIF
      IF (wg%s.GE.0.D0) THEN
        rgftm = wg%s * rgp%s
      ELSE
        rgftm = wg%s * rgp%c
      ENDIF
!
      IF (ug%c.GE.0.D0) THEN
        rgfr = ug%c * rgp%c * rb(i)
      ELSE
        rgfr = ug%c * rgp%e * rb(i)
      ENDIF
      IF (wg%c.GE.0.D0) THEN
        rgft = wg%c * rgp%c
      ELSE
        rgft = wg%c * rgp%n
      ENDIF
!
      RETURN
      END SUBROUTINE
!
!----------------------------------------------------------------------
      SUBROUTINE masfs(rlfrm, rlftm, rlfr, rlft, us, ws, rlk, i)
!
      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: rlfrm, rlftm, rlfr, rlft
      TYPE(stencil), INTENT(IN) :: us, ws, rlk
      INTEGER, INTENT(IN) :: i
!
        IF (us%w.GE.0.D0) THEN
          rlfrm = us%w * rlk%w * rb(i-1)
        ELSE
          rlfrm = us%w * rlk%c * rb(i-1)
        END IF
        IF (ws%s.GE.0.D0) THEN
          rlftm = ws%s * rlk%s
        ELSE
          rlftm = ws%s * rlk%c
        END IF
! 
        IF (us%c.GE.0.D0) THEN
          rlfr = us%c * rlk%c * rb(i)
        ELSE
          rlfr = us%c * rlk%e * rb(i)
        END IF
        IF (ws%c.GE.0.D0) THEN
          rlft = ws%c * rlk%c
        ELSE
          rlft = ws%c * rlk%n
        END IF
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      END MODULE eulerian_flux
!-------------------------------------------------------------
