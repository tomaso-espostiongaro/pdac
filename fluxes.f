!----------------------------------------------------------------------
      MODULE convective_fluxes
!
! ... This module computes convective fluxes by using an upwind scheme. 
! ... Standard compass notation is adopted to locate cell neighbours:
! ... %e=East, %w=West, %n=North, %s=South, %t=Top, %b=Bottom, etc.
! ... The computational stencil is defined in the `set_indexes' module
!
!----------------------------------------------------------------------
!
      USE grid, ONLY: dx, dy, dz, indx, indy, indz
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE time_parameters, ONLY: dt

      IMPLICIT NONE
!
      REAL*8, PRIVATE :: cs                        ! convective stream   !
      REAL*8, PRIVATE :: cn                        ! Courant number      !
      REAL*8, PRIVATE :: upwnd                     ! upwinded variable   !

      REAL*8 :: beta, muscl, lmtr
!
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE flu(fe, fn, ft, fw, fs, fb, dens, u, v, w, ijk)
!
! ... Compute the convective fluxes on East, North, and Top sides of the cell
! ... for the momentum density along x.
!
      USE dimensions
      USE grid, ONLY: myijk, fl_l
      USE set_indexes, ONLY: imjk, ijmk, ijkm
      USE set_indexes, ONLY: stencil

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
      REAL*8 :: gradp, gradm
!
      imesh = myijk( ip0_jp0_kp0_, ijk)
      i = MOD( MOD( imesh - 1, nx*ny ), nx ) + 1
      j = MOD( imesh - 1, nx*ny ) / nx + 1
      k = ( imesh - 1 ) / ( nx*ny ) + 1
!
      dxm=dx(i)+dx(i-1)
      dxp=dx(i)+dx(i+1)
      dxpp=dx(i+1)+dx(i+2)
      dym=dy(j)+dy(j-1)
      dyp=dy(j)+dy(j+1)
      dypp=dy(j+1)+dy(j+2)
      dzm=dz(k)+dz(k-1)
      dzp=dz(k)+dz(k+1)
      dzpp=dz(k+1)+dz(k+2)

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
      dens_c = (dx(i+1) * dens%c + dx(i) * dens%e) * indxp
      dens_e = (dx(i+2) * dens%e + dx(i+1) * dens%ee) * indxpp
      dens_n = (dx(i+1) * dens%n + dx(i) * dens%en) * indxp
      dens_t = (dx(i+1) * dens%t + dx(i) * dens%et) * indxp
      dens_w = (dx(i)   * dens%w + dx(i-1) * dens%c) * indxm
      dens_s = (dx(i+1) * dens%s + dx(i) * dens%es) * indxp
      dens_b = (dx(i+1) * dens%b + dx(i) * dens%eb) * indxp
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
        cs = 0.5D0*(u%c + u%w)
        IF ( cs >= 0.D0 ) fw = dens_w * u%w * cs
        IF ( cs <  0.D0 ) fw = dens_c * u%c * cs
      END IF
!
! ... on South volume bondary
!
      IF( fl_l(ijmk) /= 1 ) THEN
        cs = (dx(i+1) * v%s + dx(i) * v%es) * indxp
        IF ( cs >= 0.D0 ) fs = dens_s * u%s * cs
        IF ( cs <  0.D0 ) fs = dens_c * u%c * cs
      END IF
!
! ... on Bottom volume bondary
!
      IF( fl_l(ijkm) /= 1 ) THEN
        cs = (dx(i+1) * w%b + dx(i) * w%eb) * indxp
        IF ( cs >= 0.D0 ) fb = dens_b * u%b * cs
        IF ( cs <  0.D0 ) fb = dens_c * u%c * cs
      END IF
!
! ... MUSCL reconstruction of momentum
!
! ... on East volume boundary
!
      gradc = (indx(i+1) * (dens_e * u%e   - dens_c * u%c))
      gradw = (indx(i)   * (dens_c * u%c   - dens_w * u%w))
      grade = (indx(i+2) * (dens_ee * u%ee - dens_e * u%e))
!
      gradp = (1.D0 - beta) * gradc + beta * gradw
      gradm= (1.D0 - beta) * gradc + beta * grade

      cs = 0.5D0 * (u%c + u%e)
      cn = cs * dt * indx(i+1)
      IF ( cs >= 0.D0 ) THEN
        upwnd = dens_c * u%c + muscl*(gradp)*0.5*dx(i)
      ELSE IF ( cs < 0.D0 ) THEN
        upwnd = dens_e * u%e + muscl*(gradm)*0.5*dx(i+1)
      END IF
!
      fe = upwnd * cs
!
! ... on North volume boundary
!
      gradc = (dens_n * u%n   - dens_c * u%c) * 2.0 * indyp
      grads = (dens_c * u%c   - dens_s * u%s) * 2.0 * indym
      gradn = (dens_nn * u%nn - dens_n * u%n) * 2.0 * indypp
!
      gradp = (1.D0 - beta) * gradc + beta * grads
      gradm = (1.D0 - beta) * gradc + beta * gradn

      cs = (dx(i+1) * v%c + dx(i) * v%e) * indxp
      cn = cs * dt * 2.0 * indyp
      IF (cs >= 0.D0) THEN
        upwnd = dens_c * u%c + muscl*(gradp)*0.5*dy(j)
      ELSE IF (cs < 0.D0) THEN
        upwnd = dens_n * u%n + muscl*(gradm)*0.5*dy(j+1)
      END IF 
!
      fn = upwnd * cs
!
! ... on Top volume boundary
!
      gradc = (dens_t * u%t   - dens_c * u%c) * 2.0 * indzp
      gradb = (dens_c * u%c   - dens_b * u%b) * 2.0 * indzm
      gradt = (dens_tt * u%tt - dens_t * u%t) * 2.0 * indzpp
!
      gradp = (1.D0 - beta) * gradc + beta * gradb
      gradm = (1.D0 - beta) * gradc + beta * gradt

      cs = (dx(i+1) * w%c + dx(i) * w%e) * indxp
      cn = cs * dt * 2.0 * indzp
      IF (cs >= 0.D0) THEN
        upwnd = dens_c * u%c + muscl*(gradp)*0.5*dz(k)
      ELSE IF (cs < 0.D0) THEN
        upwnd = dens_t * u%t + muscl*(gradm)*0.5*dz(k+1)
      END IF 
!
      ft = upwnd * cs
!
      RETURN
      END SUBROUTINE flu
!----------------------------------------------------------------------
      SUBROUTINE flv(fe, fn, ft, fw, fs, fb, dens, u, v, w, ijk)
!
! ... Compute the convective fluxes on East, North, and Top sides of the cell
! ... for the momentum density along y.
!
      USE dimensions
      USE grid, ONLY: myijk, fl_l
      USE set_indexes, ONLY: imjk, ijmk, ijkm
      USE set_indexes, ONLY: stencil

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
      REAL*8 :: gradp, gradm
!
      imesh = myijk( ip0_jp0_kp0_, ijk)
      i = MOD( MOD( imesh - 1, nx*ny ), nx ) + 1
      j = MOD( imesh - 1, nx*ny ) / nx + 1
      k = ( imesh - 1 ) / ( nx*ny ) + 1
!
      dxm=dx(i)+dx(i-1)
      dxp=dx(i)+dx(i+1)
      dxpp=dx(i+1)+dx(i+2)
      dym=dy(j)+dy(j-1)
      dyp=dy(j)+dy(j+1)
      dypp=dy(j+1)+dy(j+2)
      dzm=dz(k)+dz(k-1)
      dzp=dz(k)+dz(k+1)
      dzpp=dz(k+1)+dz(k+2)

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
      dens_n = (dy(j+2) * dens%n + dy(j+1) * dens%nn) * indypp
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
      gradp = (1.D0 - beta) * gradc + beta * gradw
      gradm= (1.D0 - beta) * gradc + beta * grade

      cs = indyp * (u%c * dy(j+1) + u%n * dy(j))
      cn = cs * dt * 2.D0 * indxp
      IF ( cs >= 0.D0 ) THEN
        upwnd = dens_c * v%c + muscl*(gradp)*0.5*dx(i)
      ELSE IF ( cs < 0.D0 ) THEN
        upwnd = dens_e * v%e + muscl*(gradm)*0.5*dx(i+1)
      END IF
!
      fe = upwnd * cs
!
! ... on North volume boundary
!
      gradc = (dens_n * v%n   - dens_c * v%c) * indy(j+1)
      grads = (dens_c * v%c   - dens_s * v%s) * indy(j)
      gradn = (dens_nn * v%nn - dens_n * v%n) * indy(j+2)
!
      gradp = (1.D0 - beta) * gradc + beta * grads
      gradm = (1.D0 - beta) * gradc + beta * gradn

      cs = 0.5D0 * ( v%n + v%c )
      cn = cs * dt * indy(j+1)
      IF (cs >= 0.D0) THEN
        upwnd = dens_c * v%c + muscl*(gradp)*0.5*dy(j+1)
      ELSE IF (cs < 0.D0) THEN
        upwnd = dens_n * v%n + muscl*(gradm)*0.5*dy(j+1)
      END IF 
!
      fn = upwnd * cs
!
! ... on Top volume boundary
!
      gradc = (dens_t * v%t   - dens_c * v%c) * 2.0 * indzp
      gradb = (dens_c * v%c   - dens_b * v%b) * 2.0 * indzm
      gradt = (dens_tt * v%tt - dens_t * v%t) * 2.0 * indzpp
!
      gradp = (1.D0 - beta) * gradc + beta * gradb
      gradm = (1.D0 - beta) * gradc + beta * gradt

      cs = (dy(j+1) * w%c + dy(j) * w%n) * indyp
      cn = cs * dt * 2.0 * indzp
      IF (cs >= 0.D0) THEN
        upwnd = dens_c * v%c + muscl*(gradp)*0.5*dz(k)
      ELSE IF (cs < 0.D0) THEN
        upwnd = dens_t * v%t + muscl*(gradm)*0.5*dz(k+1)
      END IF 
!
      ft = upwnd * cs
!
      RETURN
      END SUBROUTINE flv
!----------------------------------------------------------------------
      SUBROUTINE flw(fe, fn, ft, fw, fs, fb, dens, u, v, w, ijk)
!
! ... Compute the convective fluxes on East, North, and Top sides of the cell
! ... for the momentum density along z.
!
      USE dimensions
      USE grid, ONLY: myijk, fl_l
      USE set_indexes, ONLY: imjk, ijmk, ijkm
      USE set_indexes, ONLY: stencil
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: ijk
      REAL*8, INTENT(OUT) :: fe, fn, ft, fw, fs, fb
      TYPE(stencil), INTENT(IN) :: dens, u, v, w
!
      INTEGER :: i,j,k, imesh
      REAL*8 :: dxm, dxp, dxpp, indxpp, indxp, indxm
      REAL*8 :: dym, dyp, dypp, indypp, indyp, indym
      REAL*8 :: dzp, indzp, dzm, indzm, dzpp, indzpp
      REAL*8 :: dens_c, dens_e, dens_n, dens_t
      REAL*8 :: dens_w, dens_s, dens_b
      REAL*8 :: dens_ee, dens_nn, dens_tt
      REAL*8 :: gradc, grade, gradw, gradn, grads, gradt, gradb
      REAL*8 :: gradp, gradm
!
      imesh = myijk( ip0_jp0_kp0_, ijk)
      i = MOD( MOD( imesh - 1, nx*ny ), nx ) + 1
      j = MOD( imesh - 1, nx*ny ) / nx + 1
      k = ( imesh - 1 ) / ( nx*ny ) + 1
!
      dxm=dx(i)+dx(i-1)
      dxp=dx(i)+dx(i+1)
      dxpp=dx(i+1)+dx(i+2)
      dym=dy(j)+dy(j-1)
      dyp=dy(j)+dy(j+1)
      dypp=dy(j+1)+dy(j+2)
      dzm=dz(k)+dz(k-1)
      dzp=dz(k)+dz(k+1)
      dzpp=dz(k+1)+dz(k+2)

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
! ... values of density interpolated linearly on the staggered grid
!
      dens_c = (dz(k+1) * dens%c + dz(k) * dens%t) * indzp
      dens_e = (dz(k+1) * dens%e + dz(k) * dens%et) * indzp
      dens_n = (dz(k+1) * dens%n + dz(k) * dens%nt) * indzp
      dens_t = (dz(k+2) * dens%t + dz(k+1) * dens%tt) * indzpp
      dens_w = (dz(k+1) * dens%w  + dz(k) * dens%wt) * indzp
      dens_s = (dz(k+1) * dens%s  + dz(k) * dens%st) * indzp
      dens_b = (dz(k-1) * dens%c  + dz(k) * dens%b) * indzm
!
      dens_ee = dens_e
      dens_nn = dens_n
      dens_tt = dens_t
!
! ... On boundary mantain first order accuracy (1st order Upwind).
!
! ... on West volume bondary
!
      IF (fl_l(imjk) /= 1) THEN
        cs = (dz(k+1)*u%w + dz(k)*u%wt) * indzp
        IF ( cs >= 0.D0 ) fw = dens_w * w%w * cs
        IF ( cs <  0.D0 ) fw = dens_c * w%c * cs
      END IF
!
! ... on South volume bondary
!
      IF (fl_l(ijmk) /= 1) THEN
        cs = (dz(k+1)*v%s + dz(k)*v%st) * indzp
        IF ( cs >= 0.D0 ) fs = dens_s * w%s * cs
        IF ( cs <  0.D0 ) fs = dens_c * w%c * cs
      END IF
!
! ... on Bottom volume bondary
!
      IF (fl_l(ijkm) /= 1) THEN
        cs=0.5D0*(w%b+w%c)
        IF ( cs >= 0.D0 ) fb = dens_b * w%b * cs
        IF ( cs <  0.D0 ) fb = dens_c * w%c * cs
      END IF
!
! ... MUSCL reconstruction of momentum
!
! ... on East volume boundary
!
      gradc = (2.0 * indxp * (dens_e * w%e - dens_c * w%c))
      gradw = (2.0 * indxm * (dens_c * w%c - dens_w * w%w))
      grade = (2.0 * indxpp * (dens_ee * w%ee - dens_e * w%e))
!
      gradp = (1.0 - beta) * gradc + beta * gradw
      gradm = (1.0 - beta) * gradc + beta * grade

      cs = (dz(k+1)*u%c+dz(k)*u%t)*indzp
      cn = cs * dt * 2.0 * indxp
      IF (cs >= 0.D0) THEN
        upwnd = dens_c * w%c + muscl*(gradp)*0.5*dx(i)
      ELSE IF (cs < 0.D0) THEN
        upwnd = dens_e * w%e + muscl*(gradm)*0.5*dx(i+1)
      END IF
!
      fe = upwnd * cs
!
! ... on North volume boundary
!
      gradc = (2.0 * indyp * (dens_n * w%n - dens_c * w%c))
      grads = (2.0 * indym * (dens_c * w%c - dens_s * w%s))
      gradn = (2.0 * indypp * (dens_nn * w%nn - dens_n * w%n))
!
      gradp = (1.0 - beta) * gradc + beta * grads
      gradm = (1.0 - beta) * gradc + beta * gradn

      cs = (dz(k+1)*v%c+dz(k)*v%t)*indzp
      cn = cs * dt * 2.0 * indyp
      IF (cs >= 0.D0) THEN
        upwnd = dens_c * w%c + muscl*(gradp)*0.5*dy(j)
      ELSE IF (cs < 0.D0) THEN
        upwnd = dens_n * w%n + muscl*(gradm)*0.5*dy(j+1)
      END IF
!
      fn = upwnd * cs
!
! ... on Top volume boundary
!
      gradc = (indz(k+1) * (dens_t * w%t - dens_c * w%c))
      gradt = (indz(k+2) * (dens_tt * w%tt - dens_t * w%t))
      gradb = (indz(k) * (dens_c * w%c - dens_b * w%b))
!
      gradp = (1.0 - beta) * gradc + beta * gradb
      gradm = (1.0 - beta) * gradc + beta * gradt

      cs = 0.5D0*(w%c+w%t)
      cn = cs * dt * indz(k+1)
      IF (cs >= 0.D0) THEN
        upwnd = dens_c * w%c + muscl*(gradp)*0.5*dz(k+1)
      ELSE IF (cs < 0.D0) THEN
        upwnd = dens_t * w%t + muscl*(gradm)*0.5*dz(k+1)
      END IF 
!
      ft = upwnd * cs
!
      RETURN
      END SUBROUTINE flw
!----------------------------------------------------------------------
      SUBROUTINE fsc(fe, fn, ft, fw, fs, fb, dens, field, u, v, w, ijk)
!
      USE dimensions
      USE grid, ONLY: myijk, fl_l
      USE set_indexes, ONLY: imjk, ijmk, ijkm
      USE set_indexes, ONLY: stencil
      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: fe, fn, ft, fw, fs, fb
      TYPE(stencil), INTENT(IN) :: dens, field, u, v, w
      INTEGER, INTENT(IN) :: ijk
!
      INTEGER :: i,j,k,imesh
      REAL*8 :: dxm, dxp, dxpp, indxpp, indxp, indxm
      REAL*8 :: dym, dyp, dypp, indypp, indyp, indym
      REAL*8 :: dzp, indzp, dzm, indzm, dzpp, indzpp
      REAL*8 :: gradc, grade, gradw, gradn, grads, gradt, gradb
      REAL*8 :: gradp, gradm
!
      imesh = myijk( ip0_jp0_kp0_, ijk)
      i = MOD( MOD( imesh - 1, nx*ny ), nx ) + 1
      j = MOD( imesh - 1, nx*ny ) / nx + 1
      k = ( imesh - 1 ) / ( nx*ny ) + 1
!
      dxm=dx(i)+dx(i-1)
      dxp=dx(i)+dx(i+1)
      dxpp=dx(i+1)+dx(i+2)
      dym=dy(j)+dy(j-1)
      dyp=dy(j)+dy(j+1)
      dypp=dy(j+1)+dy(j+2)
      dzm=dz(k)+dz(k-1)
      dzp=dz(k)+dz(k+1)
      dzpp=dz(k+1)+dz(k+2)

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
!
! ... on West volume boundary
!
      IF ((fl_l(imjk) /= 1)) THEN
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
      IF ((fl_l(ijmk) /= 1)) THEN
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
      IF ((fl_l(ijkm) /= 1)) THEN
        cs = w%b
        IF (cs >= 0.D0) THEN
          upwnd = dens%b * field%b
        ELSE IF (cs < 0.D0) THEN
          upwnd = dens%c * field%c
        ENDIF
        fb = upwnd * cs
      END IF
!
! ... MUSCL reconstruction of fields
!
! ... on East volume boundary
!
      gradc = 2.0 * indxp * (dens%e*field%e - dens%c*field%c)
      gradw = 2.0 * indxm * (dens%c*field%c - dens%w*field%w)
      grade = 2.0 * indxpp * (dens%ee*field%ee - dens%e*field%e)

      gradp = (1.D0 - beta) * gradc + beta * gradw
      gradm = (1.D0 - beta) * gradc + beta * grade

      cs = u%c
      cn = cs * dt * 2.0 * indxp
      IF (cs >= 0.D0) THEN
        upwnd  = dens%c*field%c + muscl*(gradp)*0.5*dx(i)
      ELSE IF (cs < 0.D0) THEN
        upwnd  = dens%e*field%e + muscl*(gradm)*0.5*dx(i+1)
      ENDIF
!
      fe = upwnd * cs
!
! ... on North volume boundary
!
      gradc = 2.0 * indyp * (dens%n*field%n - dens%c*field%c)
      grads = 2.0 * indym * (dens%c*field%c - dens%s*field%s)
      gradn = 2.0 * indypp * (dens%nn*field%nn - dens%n*field%n)

      gradp = (1.D0 - beta) * gradc + beta * grads
      gradm = (1.D0 - beta) * gradc + beta * gradn

      cs = v%c
      cn = cs * dt * 2.0 * indyp
      IF (cs >= 0.D0) THEN
        upwnd = dens%c*field%c + muscl*(gradp)*0.5*dy(j)
      ELSE IF (cs < 0.D0) THEN
        upwnd = dens%n*field%n + muscl*(gradm)*0.5*dy(j+1)
      ENDIF
!
      fn = upwnd * cs
!
! ... on Top volume boundary
!
      gradc = 2.0 * indzp * (dens%t*field%t - dens%c*field%c)
      gradb = 2.0 * indzm * (dens%c*field%c - dens%b*field%b)
      gradt = 2.0 * indzpp * (dens%tt*field%tt - dens%t*field%t)

      gradp = (1.D0 - beta) * gradc + beta * gradb
      gradm = (1.D0 - beta) * gradc + beta * gradt

      cs = w%c
      cn = cs * dt * 2.0 * indzp
      IF (cs >= 0.D0) THEN
        upwnd = dens%c*field%c + muscl*(gradp)*0.5*dz(k)
      ELSE IF (cs < 0.D0) THEN
        upwnd = dens%t*field%t + muscl*(gradm)*0.5*dz(k+1)
      ENDIF
!
      ft = upwnd * cs
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE masf(fe, fn, ft, fw, fs, fb, dens, u, v, w, ijk)
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
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE fmas(fe, fn, ft, fw, fs, fb, dens, u, v, w, ijk)
!
      USE dimensions
      USE grid, ONLY: myijk, fl_l
      USE set_indexes, ONLY: imjk, ijmk, ijkm
      USE set_indexes, ONLY: stencil
      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: fe, fn, ft, fw, fs, fb
      TYPE(stencil), INTENT(IN) :: dens, u, v, w
      INTEGER, INTENT(IN) :: ijk
!
      INTEGER :: i,j,k,imesh
      REAL*8 :: dxmm, dxm, dxp, dxpp, indxpp, indxp, indxm, indxmm
      REAL*8 :: dymm, dym, dyp, dypp, indypp, indyp, indym, indymm
      REAL*8 :: dzmm, dzm, dzp, dzpp, indzpp, indzp, indzm, indzmm
      REAL*8 :: gradc, grade, gradw, gradn, grads, gradt, gradb
      REAL*8 :: gradp, gradm
!
      imesh = myijk( ip0_jp0_kp0_, ijk)
      i = MOD( MOD( imesh - 1, nx*ny ), nx ) + 1
      j = MOD( imesh - 1, nx*ny ) / nx + 1
      k = ( imesh - 1 ) / ( nx*ny ) + 1
!
      dxmm=dx(i-1)+dx(i-2)
      dxm=dx(i)+dx(i-1)
      dxp=dx(i)+dx(i+1)
      dxpp=dx(i+1)+dx(i+2)
      dymm=dy(j-1)+dy(j-2)
      dym=dy(j)+dy(j-1)
      dyp=dy(j)+dy(j+1)
      dypp=dy(j+1)+dy(j+2)
      dzmm=dz(k-1)+dz(k-2)
      dzm=dz(k)+dz(k-1)
      dzp=dz(k)+dz(k+1)
      dzpp=dz(k+1)+dz(k+2)

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
!
! ... MUSCL reconstruction of fields
!
! ... on West volume boundary
!
      gradc = 2.0 * indxm * (dens%c - dens%w)
      gradw = 2.0 * indxmm * (dens%w - dens%ww)
      grade = 2.0 * indxp * (dens%e - dens%c)

      gradp = (1.D0 - beta) * gradc + beta * gradw
      gradm = (1.D0 - beta) * gradc + beta * grade

      cs = u%w
      cn = cs * dt * 2.0 * indxm
      IF (cs >= 0.D0) THEN
        upwnd  = dens%w + lmtr*(gradp)*0.5*dx(i-1)
      ELSE IF (cs < 0.D0) THEN
        upwnd  = dens%c + lmtr*(gradm)*0.5*dx(i)
      ENDIF
!
      fw = upwnd * cs
!
! ... on East volume boundary
!
      gradw = gradc
      gradc = grade
      grade = 2.0 * indxpp * (dens%ee - dens%e)

      gradp = (1.D0 - beta) * gradc + beta * gradw
      gradm = (1.D0 - beta) * gradc + beta * grade

      cs = u%c
      cn = cs * dt * 2.0 * indxp
      IF (cs >= 0.D0) THEN
        upwnd  = dens%c + lmtr*(gradp)*0.5*dx(i)
      ELSE IF (cs < 0.D0) THEN
        upwnd  = dens%e + lmtr*(gradm)*0.5*dx(i+1)
      ENDIF
!
      fe = upwnd * cs
!
! ... on South volume boundary
!
      gradc = 2.0 * indym * (dens%c - dens%s)
      grads = 2.0 * indymm * (dens%s - dens%ss)
      gradn = 2.0 * indyp * (dens%n - dens%c)

      gradp = (1.D0 - beta) * gradc + beta * grads
      gradm = (1.D0 - beta) * gradc + beta * gradn

      cs = v%s
      cn = cs * dt * 2.0 * indym
      IF (cs >= 0.D0) THEN
        upwnd = dens%s + lmtr*(gradp)*0.5*dy(j-1)
      ELSE IF (cs < 0.D0) THEN
        upwnd = dens%c + lmtr*(gradm)*0.5*dy(j)
      ENDIF
!
      fs = upwnd * cs
!
! ... on North volume boundary
!
      grads = gradc
      gradc = gradn
      gradn = 2.0 * indypp * (dens%nn - dens%n)

      gradp = (1.D0 - beta) * gradc + beta * grads
      gradm = (1.D0 - beta) * gradc + beta * gradn

      cs = v%c
      cn = cs * dt * 2.0 * indyp
      IF (cs >= 0.D0) THEN
        upwnd = dens%c + lmtr*(gradp)*0.5*dy(j)
      ELSE IF (cs < 0.D0) THEN
        upwnd = dens%n + lmtr*(gradm)*0.5*dy(j+1)
      ENDIF
!
      fn = upwnd * cs
!
! ... on Bottom volume boundary
!
      gradc = 2.0 * indzm * (dens%c - dens%b)
      gradb = 2.0 * indzmm * (dens%b - dens%bb)
      gradt = 2.0 * indzp * (dens%t - dens%c)

      gradp = (1.D0 - beta) * gradc + beta * gradb
      gradm = (1.D0 - beta) * gradc + beta * gradt

      cs = w%b
      cn = cs * dt * 2.0 * indzm
      IF (cs >= 0.D0) THEN
        upwnd = dens%b + lmtr*(gradp)*0.5*dz(k-1)
      ELSE IF (cs < 0.D0) THEN
        upwnd = dens%c + lmtr*(gradm)*0.5*dz(k)
      ENDIF
!
      fb = upwnd * cs
!
! ... on Top volume boundary
!
      gradb = gradc
      gradc = gradt
      gradt = 2.0 * indzpp * (dens%tt - dens%t)

      gradp = (1.D0 - beta) * gradc + beta * gradb
      gradm = (1.D0 - beta) * gradc + beta * gradt

      cs = w%c
      cn = cs * dt * 2.0 * indzp
      IF (cs >= 0.D0) THEN
        upwnd = dens%c + lmtr*(gradp)*0.5*dz(k)
      ELSE IF (cs < 0.D0) THEN
        upwnd = dens%t + lmtr*(gradm)*0.5*dz(k+1)
      ENDIF
!
      ft = upwnd * cs
!
      RETURN
      END SUBROUTINE fmas
!-----------------------------------------------------------------------
      END MODULE convective_fluxes
!-----------------------------------------------------------------------
