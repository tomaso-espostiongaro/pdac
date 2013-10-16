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
      INTERFACE ctu1_flv
        MODULE PROCEDURE ctu1_flv_3d
      END INTERFACE
      INTERFACE ctu2_flv
        MODULE PROCEDURE ctu2_flv_3d
      END INTERFACE
      INTERFACE ctu3_flv
        MODULE PROCEDURE ctu3_flv_3d
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
      SUBROUTINE ctu1_flv_3d(fe, fn, ft, dens, u, v, w, i, j, k)
!
! ... Compute the first order Corner Transport Upwind correction (step 2 in LeVeque algorithm)
!
      USE dimensions
      USE domain_mapping, ONLY: myijk, meshinds
      USE grid, ONLY:  dx, dy, dz, indx, indy, indz
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: stencil
      USE time_parameters, ONLY: dt

      IMPLICIT NONE

      REAL*8, INTENT(INOUT) :: fe, fn, ft
      TYPE(stencil), INTENT(IN) :: u, v, w, dens
      INTEGER, INTENT(IN) :: i, j, k

      REAL*8 :: uc, vc, wc, deltar, onethird
      REAL*8 :: dxp, dxm, dyp, dym, dzp, dzm, dypp
      REAL*8 :: dxi, dxip, dxim, dyj, dyjp, dyjpp, dyjm, dzk, dzkp, dzkm
      REAL*8 :: indxp, indxm, indyp, indym, indzp, indzm, indypp
      REAL*8 :: incrxp, incrxm, incryp, incrym, incrzp, incrzm, incrypp
      INTEGER :: jp2
!
      onethird = 1.D0/3.D0
      jp2 = MIN( ny, j+2 )
!
      dxi   = dx(i)
      dxip  = dx(i+1)
      dxim  = dx(i-1)
      dyj   = dy(j)
      dyjp  = dy(j+1)
      dyjpp = dy(jp2)
      dyjm  = dy(j-1)
      dzk   = dz(k)
      dzkp  = dz(k+1)
      dzkm  = dz(k-1)
!
      dxm  = dxi  + dxim
      dxp  = dxi  + dxip
      dym  = dyj  + dyjm
      dyp  = dyj  + dyjp
      dypp = dyjp + dyjpp
      dzm  = dzk  + dzkm
      dzp  = dzk  + dzkp
!
      indxp  = 2.D0/dxp
      indxm  = 2.D0/dxm
      indyp  = 2.D0/dyp
      indym  = 2.D0/dym
      indypp = 2.D0/dypp
      indzp  = 2.D0/dzp
      indzm  = 2.D0/dzm
!
      incrxp  = 0.5D0*dt*indxp
      incrxm  = 0.5D0*dt*indxm
      incryp  = 0.125D0*dt*indyp
      incrym  = 0.125D0*dt*indym
      incrypp = 0.125D0*dt*indypp
      incrzp  = 0.5D0*dt*indzp
      incrzm  = 0.5D0*dt*indzm

!
! t-c
      deltar = dens%t*v%t - dens%c*v%c
      uc = indzp*incryp*(dzkp*dyjp*(u%c + u%w   ) + &
                         dzk *dyjp*(u%t + u%wt  ) + &
                         dzkp*dyj *(u%n + u%wn  ) + &
                         dzk *dyj *(u%nt+ u%wnt))*indx(i)
      vc = (dzkp*v%c+dzk*v%t )*indyp  *incrzp
      wc = 0.5D0*indyp*(dyjp*w%c+dyj*w%n)
      fe = fe - 0.125D0*(uc+ABS(uc))*(wc-ABS(wc))*deltar
      fn = fn - 0.125D0*(vc+ABS(vc))*(wc-ABS(wc))*deltar &
              + 0.25D0*onethird*(vc+ABS(vc))*(wc-ABS(wc))*ABS(uc)*deltar 
!
! c-b
      deltar = dens%c*v%c - dens%b*v%b
      uc = indzm*incryp*(dzkm*dyjp*(u%c + u%w   ) + &
                         dzk *dyjp*(u%b + u%wb  ) + &
                         dzkm*dyj *(u%n + u%wn  ) + &
                         dzk *dyj *(u%nb+ u%wnb))*indx(i)
      vc = (dzkm*v%c+dzk*v%b )*indyp  *incrzm
      wc = 0.5D0*indyp*(dyjp*w%b+dyj*w%nb)
      fe = fe - 0.125D0*(uc+ABS(uc))*(wc+ABS(wc))*deltar
      fn = fn - 0.125D0*(vc+ABS(vc))*(wc+ABS(wc))*deltar &
              + 0.25D0*onethird*(vc+ABS(vc))*(wc+ABS(wc))*ABS(uc)*deltar 
!
! n-c
      deltar = dens%n*v%n - dens%c*v%c
      uc = 0.5D0*(u%n+u%wn)*dt*indx(i)
      vc = 0.5D0*(v%c+v%n )
      wc = 0.5D0*(w%n+w%nb)*dt*indz(k)
      ft = ft - 0.125D0*(wc+ABS(wc))*(vc-ABS(vc))*deltar
      fe = fe - 0.125D0*(uc+ABS(uc))*(vc-ABS(vc))*deltar &
              + 0.25D0*onethird*(uc+ABS(uc))*(vc-ABS(vc))*ABS(wc)*deltar 
!
! c-s
      deltar = dens%c*v%c - dens%s*v%s
      uc = 0.5D0*(u%c+u%w)*dt*indx(i)
      vc = 0.5D0*(v%c+v%s)
      wc = 0.5D0*(w%c+w%b)*dt*indz(k)
      ft = ft - 0.125D0*(wc+ABS(wc))*(vc+ABS(vc))*deltar
      fe = fe - 0.125D0*(uc+ABS(uc))*(vc+ABS(vc))*deltar &
              + 0.25D0*onethird*(uc+ABS(uc))*(vc+ABS(vc))*ABS(wc)*deltar 
!
! e-c
      deltar = dens%e*v%e - dens%c*v%c
      uc = 0.5D0*indyp*(dyjp*u%c+dyj*u%n)
      vc = (dxip*v%c+dxi*v%e )*indyp  *incrxp
      wc = indxp*incryp*(dxip*dyjp*(w%c + w%b  ) + &
                         dxi *dyjp*(w%e + w%eb ) + &
                         dxip*dyj *(w%n + w%nb ) + &
                         dxi *dyj *(w%en+ w%enb))*indz(k)
      fn = fn - 0.125D0*(vc+ABS(vc))*(uc-ABS(uc))*deltar
      ft = ft - 0.125D0*(wc+ABS(wc))*(uc-ABS(uc))*deltar &
              + 0.25D0*onethird*(wc+ABS(wc))*(uc-ABS(uc))*ABS(vc)*deltar 
!
! c-w
      deltar = dens%c*v%c - dens%w*v%w
      uc = 0.5D0*indyp*(dyjp*u%w+dyj*u%wn)
      vc = (dxim*v%c+dxi*v%w )*indyp  *incrxm
      wc = indxm*incryp*(dxim*dyjp*(w%c + w%b  ) + &
                         dxi *dyjp*(w%w + w%wb ) + &
                         dxim*dyj *(w%n + w%nb ) + &
                         dxi *dyj *(w%wn+ w%wnb))*indz(k)
      fn = fn - 0.125D0*(vc+ABS(vc))*(uc+ABS(uc))*deltar
      ft = ft - 0.125D0*(wc+ABS(wc))*(uc+ABS(uc))*deltar &
              + 0.25D0*onethird*(wc+ABS(wc))*(uc+ABS(uc))*ABS(vc)*deltar 
!
! wt-w
      deltar = dens%wt*v%wt - dens%w*v%w
      IF (i <= 2) deltar = 0.D0
      uc = indzp*incryp*(dzkp*dyjp*(u%w  + u%ww   ) + &
                                dzk *dyjp*(u%wt + u%wwt  ) + &
                                dzkp*dyj *(u%wn + u%wwn  ) + &
                                dzk *dyj *(u%wnt+ u%wwnt))*indx(i-1)
      vc = (dzkp*v%w+dzk*v%wt  )*indyp    *incrzp
      wc = 0.5D0*indyp*(dyjp*w%w+dyj*w%wn)
      fn = fn - 0.125D0*onethird*(wc-ABS(wc))*(vc+ABS(vc))*ABS(uc+ABS(uc))*deltar
!
! w-wb
      deltar = dens%w*v%w - dens%wb*v%wb
      IF (i <= 2) deltar = 0.D0
      uc = indzm*incryp*(dzkm*dyjp*(u%w  + u%ww   ) + &
                         dzk *dyjp*(u%wb + u%wwb  ) + &
                         dzkm*dyj *(u%wn + u%wwn  ) + &
                         dzk *dyj *(u%wnb+ u%wwnb))*indx(i-1)
      vc = (dzkm*v%w +dzk*v%wb )*indyp    *incrzm
      wc = 0.5D0*indyp*(dyjp*w%wb+dyj*w%wnb)
      fn = fn - 0.125D0*onethird*(wc+ABS(wc))*(vc+ABS(vc))*ABS(uc+ABS(uc))*deltar
!
! et-e
      deltar = dens%et*v%et - dens%e*v%e
      IF (i >= nx-1) deltar = 0.D0
      uc = indzp*incryp*(dzkp*dyjp*(u%c + u%e   ) + &
                         dzk *dyjp*(u%t + u%et  ) + &
                         dzkp*dyj *(u%n + u%en  ) + &
                         dzk *dyj *(u%nt+ u%ent))*indx(i+1)
      vc = (dzkp*v%e +dzk*v%et )*indyp    *incrzp
      wc = 0.5D0*indyp*(dyjp*w%e+dyj*w%en)
      fe = fe - 0.125D0*(wc-ABS(wc))*(uc-ABS(uc))*deltar
      fn = fn - 0.125D0*onethird*(wc-ABS(wc))*(vc+ABS(vc))*ABS(uc-ABS(uc))*deltar
!
! e-eb
      deltar = dens%e*v%e - dens%eb*v%eb
      IF (i >= nx-1) deltar = 0.D0
      uc = indzm*incryp*(dzkm*dyjp*(u%c + u%e   ) + &
                         dzk *dyjp*(u%b + u%eb  ) + &
                         dzkm*dyj *(u%n + u%en  ) + &
                         dzk *dyj *(u%nb+ u%enb))*indx(i+1)
      vc = (dzkm*v%e+dzk*v%eb)  *indyp    *incrzm
      wc = 0.5D0*indyp*(dyjp*w%eb+dyj*w%enb)
      fe = fe - 0.125D0*(wc+ABS(wc))*(uc-ABS(uc))*deltar
      fn = fn - 0.125D0*onethird*(wc+ABS(wc))*(vc+ABS(vc))*ABS(uc-ABS(uc))*deltar
!
! en-e
      deltar = dens%en*v%en - dens%e*v%e
      IF (i >= nx-1) deltar = 0.D0
      uc = 0.5D0*(u%n + u%en )*dt*indx(i+1)
      vc = 0.5D0*(v%e + v%en )
      wc = 0.5D0*(w%en+ w%enb)*dt*indz(k)
      fe = fe - 0.125D0*(vc-ABS(vc))*(uc-ABS(uc))*deltar &
              + 0.25D0*onethird*(vc-ABS(vc))*(uc-ABS(uc))*ABS(wc)*deltar
!
! e-es
      deltar = dens%e*v%e - dens%es*v%es
      IF (i >= nx-1) deltar = 0.D0
      uc = 0.5D0*(u%c + u%e )*dt*indx(i+1)
      vc = 0.5D0*(v%e + v%es)
      wc = 0.5D0*(w%e + w%eb)*dt*indz(k)
      fe = fe - 0.125D0*(vc+ABS(vc))*(uc-ABS(uc))*deltar &
              + 0.25D0*onethird*(vc+ABS(vc))*(uc-ABS(uc))*ABS(wc)*deltar
!
! nt-n
      deltar = dens%nt*v%nt - dens%n*v%n
      IF (j >= ny-1) deltar = 0.D0
      uc = indzp*incrypp*(dzkp*dyjpp*(u%n  + u%wn  ) + &
                          dzk *dyjpp*(u%nt + u%wnt ) + &
                          dzkp*dyjp *(u%nn + u%wnn ) + &
                          dzk *dyjp *(u%nnt+ u%wnnt))*indx(i)
      vc = (dzkp*v%n +dzk *v%nt )*indypp *incrzp
      wc = 0.5D0*indypp*(dyjpp*w%n+dyjp*w%nn)
      fn = fn - 0.125D0*(vc-ABS(vc))*(wc-ABS(wc))*deltar &
              + 0.25D0*onethird*(vc-ABS(vc))*(wc-ABS(wc))*ABS(uc)*deltar
!
! n-nb
      deltar = dens%n*v%n - dens%nb*v%nb
      IF (j >= ny-1) deltar = 0.D0
      uc = indzm*incrypp*(dzkm*dyjpp*(u%n  + u%wn  ) + &
                          dzk *dyjpp*(u%nb + u%wnb ) + &
                          dzkm*dyjp *(u%nn + u%wnn ) + &
                          dzk *dyjp *(u%nnb+ u%wnnb))*indx(i)
      vc = (dzkm*v%n +dzk *v%nb )*indypp *incrzm
      wc = 0.5D0*indypp*(dyjpp*w%nb+dyjp*w%nnb)
      fn = fn - 0.125D0*(vc-ABS(vc))*(wc+ABS(wc))*deltar &
              + 0.25D0*onethird*(vc-ABS(vc))*(wc+ABS(wc))*ABS(uc)*deltar
!
! nt-t
      deltar  = dens%nt*v%nt - dens%t*v%t
      IF (k >= nz-1) deltar = 0.D0
      uc = 0.5D0*(u%nt+u%wnt)*dt*indx(i)
      vc = 0.5D0*(v%t +v%nt )
      wc = 0.5D0*(w%nt+w%n  )*dt*indz(k+1)
      ft = ft - 0.125D0*(vc-ABS(vc))*(wc-ABS(wc))*deltar
      fe = fe - 0.125D0*onethird*(vc-ABS(vc))*(uc+ABS(uc))*ABS(wc-ABS(wc))*deltar
!     
! t-st
      deltar  = dens%t*v%t - dens%st*v%st
      IF (k >= nz-1) deltar = 0.D0
      uc = 0.5D0*(u%t+u%wt)*dt*indx(i)
      vc = 0.5D0*(v%t+v%st)
      wc = 0.5D0*(w%c+w%t )*dt*indz(k+1)
      ft = ft - 0.125D0*(vc+ABS(vc))*(wc-ABS(wc))*deltar
      fe = fe - 0.125D0*onethird*(vc+ABS(vc))*(uc+ABS(uc))*ABS(wc-ABS(wc))*deltar
!     
! nb-b
      deltar = dens%nb*v%nb - dens%b*v%b
      IF (k <= 2) deltar = 0.D0
      uc = 0.5D0*(u%nb+u%wnb)*dt*indx(i)
      vc = 0.5D0*(v%nb+v%b)
      wc = 0.5D0*(w%nb+w%nbb)*dt*indz(k-1)
      fe = fe - 0.125D0*onethird*(vc-ABS(vc))*(uc+ABS(uc))*ABS(wc+ABS(wc))*deltar
!     
! b-sb
      deltar = dens%b*v%b - dens%sb*v%sb
      IF (k <= 2) deltar = 0.D0
      uc = 0.5D0*(u%b+u%wb)*dt*indx(i)
      vc = 0.5D0*(v%b+v%sb)
      wc = 0.5D0*(w%b+w%bb)*dt*indz(k-1)
      fe = fe - 0.125D0*onethird*(vc+ABS(vc))*(uc+ABS(uc))*ABS(wc+ABS(wc))*deltar
!     
! en-n
      deltar = dens%en*v%en - dens%n*v%n
      IF (j >= ny-1) deltar = 0.D0
      uc = 0.5D0*indypp*(dyjpp*u%n+dyjp*u%nn)
      vc = (dxip*v%n +dxi *v%en )*indypp *incrxp
      wc = indxp*incrypp*(dxip*dyjpp*(w%n  + w%nb  ) + &
                          dxi *dyjpp*(w%en + w%enb ) + &
                          dxip*dyjp *(w%nn + w%nnb ) + &
                          dxi *dyjp *(w%enn+ w%ennb))*indz(k)
      fn = fn - 0.125D0*(vc-ABS(vc))*(uc-ABS(uc))*deltar
      ft = ft - 0.125D0*onethird*(uc-ABS(uc))*(wc+ABS(wc))*ABS(vc-ABS(vc))*deltar
!     
! n-wn
      deltar = dens%n*v%n - dens%wn*v%wn
      IF (j >= ny-1) deltar = 0.D0
      uc = 0.5D0*indypp*(dyjpp*u%wn+dyjp*u%wnn)
      vc = (dxim*v%n +dxi *v%wn )*indypp *incrxm
      wc = indxm*incrypp*(dxim*dyjp*(w%n  + w%nb  ) + &
                          dxi *dyjp*(w%wn + w%wnb ) + &
                          dxim*dyj *(w%nn + w%nnb ) + &
                          dxi *dyj *(w%wnn+ w%wnnb))*indz(k)
      fn = fn - 0.125D0*(vc-ABS(vc))*(uc+ABS(uc))*deltar
      ft = ft - 0.125D0*onethird*(uc+ABS(uc))*(wc+ABS(wc))*ABS(vc-ABS(vc))*deltar
!     
! es-s
      deltar = dens%es*v%es - dens%s*v%s
      IF (j <= 2) deltar = 0.D0
      uc = 0.5D0*indym*(dyjm*u%c+dyj*u%s)
      vc = (dxip*v%s+dxi*v%es)*indym  *incrxp
      wc = indxp*incrym*(dxip*dyjm*(w%c + w%b  ) + &
                         dxi *dyjm*(w%e + w%eb ) + &
                         dxip*dyj *(w%s + w%sb ) + &
                         dxi *dyj *(w%es+ w%esb))*indz(k)
      ft = ft - 0.125D0*onethird*(uc-ABS(uc))*(wc+ABS(wc))*ABS(vc+ABS(vc))*deltar
!     
! s-ws
      deltar = dens%s*v%s - dens%ws*v%ws
      IF (j <= 2) deltar = 0.D0
      uc = 0.5D0*indym*(dyjm*u%w+dyj*u%ws)
      vc = (dxim*v%s+dxi*v%ws)*indym  *incrxm
      wc = indxm*incrym*(dxim*dyjm*(w%c + w%b  ) + &
                         dxi *dyjm*(w%w + w%wb ) + &
                         dxim*dyj *(w%s + w%sb ) + &
                         dxi *dyj *(w%ws+ w%wsb))*indz(k)
      ft = ft - 0.125D0*onethird*(uc+ABS(uc))*(wc+ABS(wc))*ABS(vc+ABS(vc))*deltar
!     
! et-t
      deltar = dens%et*v%et - dens%t*v%t
      IF (k >= nz-1) deltar = 0.D0
      uc = 0.5D0*indyp*(dyjp*u%t+dyj*u%nt)
      vc = (dxip*v%t+dxi*v%et)*indyp*incrxp
      wc = indxp*incryp*(dxip*dyjp*(w%c + w%t  ) + &
                         dxi *dyjp*(w%e + w%et ) + &
                         dxip*dyj *(w%n + w%nt ) + &
                         dxi *dyj *(w%en+ w%ent))*indz(k+1)
      ft = ft - 0.125D0*(wc-ABS(wc))*(uc-ABS(uc))*deltar &
              + 0.25D0*onethird*(wc-ABS(wc))*(uc-ABS(uc))*ABS(vc)*deltar
!     
! t-wt
      deltar = dens%t*v%t - dens%wt*v%wt
      IF (k >= nz-1) deltar = 0.D0
      uc = 0.5D0*indyp*(dyjp*u%wt+dyj*u%wnt)
      vc = (dxim*v%t +dxi*v%wt )*indyp    *incrxm
      wc = indxm*incryp*(dxim*dyjp*(w%c + w%t  ) + &
                         dxi *dyjp*(w%w + w%wt ) + &
                         dxim*dyj *(w%n + w%nt ) + &
                         dxi *dyj *(w%wn+ w%wnt))*indz(k+1)
      ft = ft - 0.125D0*(wc-ABS(wc))*(uc+ABS(uc))*deltar &
              + 0.25D0*onethird*(wc-ABS(wc))*(uc+ABS(uc))*ABS(vc)*deltar
!     
! wnt-wn
      deltar = dens%wnt*v%wnt - dens%wn*v%wn
      IF (i <= 2) deltar = 0.D0
      IF (j >= ny-1) deltar = 0.D0
      uc = indzp*incrypp*(dzkp*dyjpp*(u%wn  + u%wwn  ) + &
                          dzk *dyjpp*(u%wnt + u%wwnt ) + &
                          dzkp*dyjp *(u%wnn + u%wwnn ) + &
                          dzk *dyjp *(u%wnnt+ u%wwnnt))*indx(i-1)
      vc = (dzkp*v%wn +dzk  *v%wnt )*indypp   *incrzp
      wc = 0.5D0*indypp*(dyjpp*w%wn+dyjp*w%wnn)
      fn = fn - 0.125D0*onethird*(wc-ABS(wc))*(vc-ABS(vc))*ABS(uc+ABS(uc))*deltar
!     
! wn-wnb
      deltar = dens%wn*v%wn - dens%wnb*v%wnb
      IF (i <= 2) deltar = 0.D0
      IF (j >= ny-1) deltar = 0.D0
      uc = indzm*incrypp*(dzkm*dyjpp*(u%wn  + u%wwn  ) + &
                          dzk *dyjpp*(u%wnb + u%wwnb ) + &
                          dzkm*dyjp *(u%wnn + u%wwnn ) + &
                          dzk *dyjp *(u%wnnb+ u%wwnnb))*indx(i-1)
      vc = (dzkm*v%wn +dzk  *v%wnb )*indypp   *incrzm
      wc = 0.5D0*indypp*(dyjpp*w%wnb+dyjp*w%wnnb)
      fn = fn - 0.125D0*onethird*(wc+ABS(wc))*(vc-ABS(vc))*ABS(uc+ABS(uc))*deltar
!     
! ent-et
      deltar = dens%ent*v%ent - dens%et*v%et
      IF (i >= nx-1) deltar = 0.D0
      IF (k >= nz-1) deltar = 0.D0
      uc = 0.5D0*(u%ent+u%nt)*dt*indx(i+1)
      vc = 0.5D0*(v%ent+v%et)
      wc = 0.5D0*(w%ent+w%en)*dt*indz(k+1)
      fe = fe - 0.125D0*onethird*(vc-ABS(vc))*(uc-ABS(uc))*ABS(wc-ABS(wc))*deltar
!     
! et-est
      deltar = dens%et*v%et - dens%est*v%est
      IF (i >= nx-1) deltar = 0.D0
      IF (k >= nz-1) deltar = 0.D0
      uc = 0.5D0*(u%et+u%t)*dt*indx(i+1)
      vc = 0.5D0*(v%et+v%est)
      wc = 0.5D0*(w%et+w%e)*dt*indz(k+1)
      fe = fe - 0.125D0*onethird*(vc+ABS(vc))*(uc-ABS(uc))*ABS(wc-ABS(wc))*deltar
!     
! enb-eb
      deltar = dens%enb*v%enb - dens%eb*v%eb
      IF (i >= nx-1) deltar = 0.D0
      IF (k <= 2) deltar = 0.D0
      uc = 0.5D0*(u%enb+u%nb)*dt*indx(i+1)
      vc = 0.5D0*(v%enb+v%eb)
      wc = 0.5D0*(w%enb+w%enbb)*dt*indz(k-1)
      fe = fe - 0.125D0*onethird*(vc-ABS(vc))*(uc-ABS(uc))*ABS(wc+ABS(wc))*deltar
!     
! eb-esb
      deltar = dens%eb*v%eb - dens%esb*v%esb
      IF (i >= nx-1) deltar = 0.D0
      IF (k <= 2) deltar = 0.D0
      uc = 0.5D0*(u%eb+u%b)*dt*indx(i+1)
      vc = 0.5D0*(v%eb+v%esb)
      wc = 0.5D0*(w%eb+w%ebb)*dt*indz(k-1)
      fe = fe - 0.125D0*onethird*(vc+ABS(vc))*(uc-ABS(uc))*ABS(wc+ABS(wc))*deltar
!     
! ent-en
      deltar = dens%ent*v%ent - dens%en*v%en
      IF (i >= nx-1) deltar = 0.D0
      IF (j >= ny-1) deltar = 0.D0
      uc = indzp*incrypp*(dzkp*dyjpp*(u%n  + u%en  ) + &
                          dzk *dyjpp*(u%nt + u%ent ) + &
                          dzkp*dyjp *(u%nn + u%enn ) + &
                          dzk *dyjp *(u%nnt+ u%ennt))*indx(i+1)
      vc = (dzkp*v%en+dzk*v%ent)*indypp*incrzp
      fn = fn - 0.125D0*onethird*(wc-ABS(wc))*(vc-ABS(vc))*ABS(uc-ABS(uc))*deltar
!     
! en-enb
      deltar = dens%en*v%en - dens%enb*v%enb
      IF (i >= nx-1) deltar = 0.D0
      IF (j >= ny-1) deltar = 0.D0
      uc = indzm*incrypp*(dzkm*dyjpp*(u%n  + u%en  ) + &
                          dzk *dyjpp*(u%nb + u%enb ) + &
                          dzkm*dyjp *(u%nn + u%enn ) + &
                          dzk *dyjp *(u%nnb+ u%ennb))*indx(i+1)
      vc = (dzkm*v%en +dzk *v%enb )*indypp   *incrzm
      wc = 0.5D0*indypp*(dyjpp*w%enb+dyjp*w%ennb)
      fn = fn - 0.125D0*onethird*(wc+ABS(wc))*(vc-ABS(vc))*ABS(uc-ABS(uc))*deltar
!     
! ent-nt
      deltar = dens%ent*v%ent - dens%nt*v%nt
      IF (k >= nz-1) deltar = 0.D0
      IF (j >= ny-1) deltar = 0.D0
      uc = 0.5D0*indypp*(dyjpp*u%nt+dyjp*u%nnt)
      vc = (dxip*v%nt+dxi*v%ent)*indypp*incrxp
      wc = indxp*incrypp*(dxip*dyjpp*(w%n  + w%nt  ) + &
                          dxi *dyjpp*(w%en + w%ent ) + &
                          dxip*dyjp *(w%nn + w%nnt ) + &
                          dxi *dyjp *(w%enn+ w%ennt))*indz(k+1)
      ft = ft - 0.125D0*onethird*(wc-ABS(wc))*(uc-ABS(uc))*ABS(vc-ABS(vc))*deltar
!     
! nt-wnt
      deltar = dens%nt*v%nt - dens%wnt*v%wnt
      IF (k >= nz-1) deltar = 0.D0
      IF (j >= ny-1) deltar = 0.D0
      uc = 0.5D0*indypp*(dyjpp*u%wnt+dyjp*u%wnnt)
      vc = (dxim*v%nt+dxi*v%wnt)*indypp*incrxm
      wc = indxm*incrypp*(dxim*dyjpp*(w%n  + w%nt  ) + &
                          dxi *dyjpp*(w%wn + w%wnt ) + &
                          dxim*dyjp *(w%nn + w%nnt ) + &
                          dxi *dyjp *(w%wnn+ w%wnnt))*indz(k+1)
      ft = ft - 0.125D0*onethird*(wc-ABS(wc))*(uc+ABS(uc))*ABS(vc-ABS(vc))*deltar
!     
! est-st
      deltar = dens%est*v%est - dens%st*v%st
      IF (k >= nz-1) deltar = 0.D0
      IF (j <= 2) deltar = 0.D0
      uc = 0.5D0*indym*(dyjm*u%t+dyj*u%st)
      vc = (dxip*v%st+dxi*v%est)*indym*incrxp
      wc = indxp*incrym*(dxip*dyjm*(w%c + w%t  ) + &
                         dxi *dyjm*(w%e + w%et ) + &
                         dxip*dyj *(w%s + w%st ) + &
                         dxi *dyj *(w%es+ w%est))*indz(k+1)
      ft = ft - 0.125D0*onethird*(wc-ABS(wc))*(uc-ABS(uc))*ABS(vc+ABS(vc))*deltar
!     
! st-wst
      deltar = dens%st*v%st - dens%wst*v%wst
      IF (k >= nz-1) deltar = 0.D0
      IF (j <= 2) deltar = 0.D0
      uc = 0.5D0*indym*(dyjm*u%wt+dyj*u%wst)
      vc = (dxim*v%st+dxi*v%wst)*indym*incrxm
      wc = indxm*incrym*(dxim*dyjm*(w%c + w%t  ) + &
                         dxi *dyjm*(w%w + w%wt ) + &
                         dxim*dyj *(w%s + w%st ) + &
                         dxi *dyj *(w%ws+ w%wst))*indz(k+1)
      ft = ft - 0.125D0*onethird*(wc-ABS(wc))*(uc+ABS(uc))*ABS(vc+ABS(vc))*deltar
!
      RETURN
      END SUBROUTINE ctu1_flv_3d
!----------------------------------------------------------------------
      SUBROUTINE ctu2_flv_3d(fe, fn, ft, dens, u, v, w, i, j, k)
!
! ... Compute the second order Corner Transport Upwind correction (step 3 in LeVeque algorithm)
!
      USE dimensions
      USE domain_mapping, ONLY: myijk, meshinds
      USE grid, ONLY: dx, dy, dz, indx, indy, indz, fluid, flag
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: stencil
      USE set_indexes, ONLY: imjk, ijkm, ipjk, ijkp, ijmk, ijpk
      USE time_parameters, ONLY: dt
      USE flux_limiters, ONLY: lv, limiters
!
      IMPLICIT NONE

      REAL*8, INTENT(INOUT) :: fe, ft, fn
      TYPE(stencil), INTENT(IN) :: dens, u, v, w 
      INTEGER, INTENT(IN) :: i, j, k 
!
      INTEGER :: im2, ip2, jm2, jp2, km2, kp2
      REAL*8 :: uref, vref, wref, deltar, erre, lim
!
      REAL*8 :: dxpp, dxp, dxmm, dxm, indxpp, indxp, indxm, indxmm
      REAL*8 :: dypp, dyp, dymm, dym, indypp, indyp, indym, indymm
      REAL*8 :: dzpp, dzp, dzmm, dzm, indzpp, indzp, indzm, indzmm
      REAL*8 :: gradc, grade, gradw, gradn, grads, gradt, gradb
!
      im2 = MAX(  1, i-2 )
      ip2 = MIN( nx, i+2 )
      jm2 = MAX(  1, j-2 )
      jp2 = MIN( ny, j+2 )
      km2 = MAX(  1, k-2 )
      kp2 = MIN( nz, k+2 )
!
      dxm  = dx(i)  + dx(i-1)
      dxp  = dx(i)  + dx(i+1)
      dxpp = dx(i+1)+ dx(ip2)
      dxmm = dx(i-1)+ dx(im2)
      dym  = dy(j)  + dy(j-1)
      dyp  = dy(j)  + dy(j+1)
      dypp = dy(j+1)+ dy(jp2)
      dymm = dy(j-1)+ dy(jm2)
      dzm  = dz(k)  + dz(k-1)
      dzp  = dz(k)  + dz(k+1)
      dzpp = dz(k+1)+ dz(kp2)
      dzmm = dz(k-1)+ dz(km2)
!
      indxm  = 2.D0/dxm
      indxp  = 2.D0/dxp
      indxpp = 2.D0/dxpp
      indxmm = 2.D0/dxmm
      indym  = 2.D0/dym
      indyp  = 2.D0/dyp
      indypp = 2.D0/dypp
      indymm = 2.D0/dymm
      indzm  = 2.D0/dzm
      indzp  = 2.D0/dzp
      indzpp = 2.D0/dzpp
      indzmm = 2.D0/dzmm
!
! e-c
      gradw = indxm  * (v%c *dens%c  - v%w*dens%w)
      gradc = indxp  * (v%e *dens%e  - v%c*dens%c)
      grade = indxpp * (v%ee*dens%ee - v%e*dens%e)
      lim  = 0.D0
      erre = 0.D0
      uref = 0.5D0*indyp*(dy(j+1)*u%c + dy(j)*u%n)
      IF (uref >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=nx-1 .AND. flag(ipjk)==fluid) CALL limiters(lv,lim,erre)
      uref = ABS(0.5D0*indyp*(dy(j+1)*u%c + dy(j)*u%n))
      deltar = v%e*dens%e - v%c*dens%c
      fe = fe + 0.5D0*uref*(1-dt*indxp*uref)*deltar*lim
!
! n-c
      gradc = indy(j+1)*(v%n *dens%n  - v%c*dens%c)
      gradn = indy(jp2)*(v%nn*dens%nn - v%n*dens%n)
      grads = indy(j)  *(v%c *dens%c  - v%s*dens%s)
      lim  = 0.D0
      erre = 0.D0
      vref = 0.5D0*(v%c + v%n)
      IF (vref >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = grads/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradn/gradc
      END IF
      IF (j/=ny-1 .AND. flag(ijpk)==fluid) CALL limiters(lv,lim,erre)
      vref = ABS(0.5D0*(v%c + v%n))
      deltar = v%n*dens%n - v%c*dens%c
      fn = fn + 0.5D0*vref*(1-dt*indy(j+1)*vref)*deltar*lim
!
! t-c
      gradc = indzp *(v%t *dens%t  - v%c*dens%c)
      gradt = indzpp*(v%tt*dens%tt - v%t*dens%t)
      gradb = indzm *(v%c *dens%c  - v%b*dens%b)
      lim  = 0.D0
      erre = 0.D0
      wref = 0.5D0*indyp*(dy(j)*w%n + dy(j+1)*w%c)
      IF (wref >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(ijkp)==fluid) CALL limiters(lv,lim,erre)
      wref = ABS(0.5D0*indyp*(dy(j)*w%n + dy(j+1)*w%c))
      deltar = v%t*dens%t - v%c*dens%c
      ft = ft + 0.5D0*wref*(1-dt*indzp*wref)*deltar*lim
!
      RETURN
      END SUBROUTINE ctu2_flv_3d
!----------------------------------------------------------------------
      SUBROUTINE ctu3_flv_3d(fe, fn, ft, dens, u, v, w, i, j, k)
!
! ... Compute the second order Corner Transport Upwind correction (step 4 in LeVeque algorithm)
!
      USE dimensions
      USE domain_mapping, ONLY: myijk
      USE grid, ONLY: dx, dy, dz, indx, indy, indz, fluid, flag
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE set_indexes, ONLY: stencil
      USE set_indexes, ONLY: imjk, ijkm, ipjk, ijkp, ijmk, ijpk
      USE set_indexes, ONLY: ipjpk, ipjmk, imjpk, imjmk, ijpkp, ijmkp, ijpkm, ijmkm
      USE set_indexes, ONLY: ipjkp, ipjkm, imjkp, imjkm, ipjpkp
      USE set_indexes, ONLY: ipjpkm, ipjmkp, imjpkp, ipjmkm, imjpkm, imjmkp
      USE time_parameters, ONLY: dt
      USE flux_limiters, ONLY: limiters, lv

      IMPLICIT NONE

      REAL*8, INTENT(INOUT) :: fe, ft, fn
      TYPE(stencil), INTENT(IN) :: dens, u, v, w 
      INTEGER, INTENT(IN) :: i, j, k 

      INTEGER :: im2, ip2, jm2, jp2, km2, kp2
      REAL*8 :: uc, vc, wc, deltar, erre, lim, s
      REAL*8 :: dxpp, dxp, dxm, dxmm, indxpp, indxp, indxm, indxmm
      REAL*8 :: dypp, dyp, dym, indypp, indyp, indym
      REAL*8 :: dzpp, dzp, dzmm, dzm, indzpp, indzp, indzm, indzmm
      REAL*8 :: incrxp, incrxm, incryp, incrym, incrypp, incrzp, incrzm
      REAL*8 :: gradc, grade, gradw, gradn, grads, gradt, gradb
!
      im2 = MAX(  1, i-2 )
      ip2 = MIN( nx, i+2 )
      jm2 = MAX(  1, j-2 )
      jp2 = MIN( ny, j+2 )
      km2 = MAX(  1, k-2 )
      kp2 = MIN( nz, k+2 )
!
      dxm  = dx(i)  +dx(i-1)
      dxp  = dx(i)  +dx(i+1)
      dxpp = dx(i+1)+dx(ip2)
      dxmm = dx(i-1)+dx(im2)
      dym  = dy(j)  +dy(j-1)
      dyp  = dy(j)  +dy(j+1)
      dypp = dy(j+1)+dy(jp2)
      dzm  = dz(k)  +dz(k-1)
      dzp  = dz(k)  +dz(k+1)
      dzpp = dz(k+1)+dz(kp2)
      dzmm = dz(k-1)+dz(km2)
!
      indxm  = 2.D0/dxm
      indxp  = 2.D0/dxp
      indxpp = 2.D0/dxpp
      indxmm = 2.D0/dxmm
      indym  = 2.D0/dym
      indyp  = 2.D0/dyp
      indypp = 2.D0/dypp
      indzm  = 2.D0/dzm
      indzp  = 2.D0/dzp
      indzpp = 2.D0/dzpp
      indzmm = 2.D0/dzmm
!
      incrxm  = 0.5D0*dt*indxm
      incrxp  = 0.5D0*dt*indxp
      incrym  = 0.125D0*dt*indym
      incryp  = 0.125D0*dt*indyp
      incrypp = 0.125D0*dt*indypp
      incrzm  = 0.5D0*dt*indzm
      incrzp  = 0.5D0*dt*indzp
!
! c-w
      deltar  = dens%c*v%c - dens%w*v%w
      uc = 0.5D0*indyp*(dy(j+1)*u%w+dy(j)*u%wn)
      vc = (dx(i-1)*v%c+dx(i)*v%w )*indyp  *incrxm
      wc = indxm*incryp*(dx(i-1)*dy(j+1)*(w%c + w%b  ) + &
                                dx(i)  *dy(j+1)*(w%w + w%wb ) + &
                                dx(i-1)*dy(j)  *(w%n + w%nb ) + &
                                dx(i)  *dy(j)  *(w%wn+ w%wnb))*indz(k)
      gradc = indxm  * (dens%c*v%c - dens%w *v%w )
      grade = indxp  * (dens%e*v%e - dens%c *v%c )
      gradw = indxmm * (dens%w*v%w - dens%ww*v%ww)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=2 .AND. flag(imjk)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxm*ABS(uc))*deltar*lim
      fn = fn + 0.5D0*(vc+ABS(vc))*s
      ft = ft + 0.5D0*(1-ABS(vc))*(wc+ABS(wc))*s
!
! e-c
      deltar  = dens%e*v%e - dens%c*v%c
      uc = 0.5D0*indyp*(dy(j+1)*u%c+dy(j)*u%n)
      vc = (dx(i+1)*v%c+dx(i)*v%e )*indyp  *incrxp
      wc = indxp*incryp*(dx(i+1)*dy(j+1)*(w%c + w%b  ) + &
                                dx(i)  *dy(j+1)*(w%e + w%eb ) + &
                                dx(i+1)*dy(j)  *(w%n + w%nb ) + &
                                dx(i)  *dy(j)  *(w%en+ w%enb))*indz(k)
      gradc = indxp  * (dens%e *v%e  - dens%c*v%c)
      grade = indxpp * (dens%ee*v%ee - dens%e*v%e)
      gradw = indxm  * (dens%c *v%c  - dens%w*v%w)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=nx-1 .AND. flag(ipjk)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxp*ABS(uc))*deltar*lim
      fn = fn - 0.5D0*(vc+ABS(vc))*s
      ft = ft - 0.5D0*(1-ABS(vc))*(wc+ABS(wc))*s
!
! en-n
      deltar  = dens%en*v%en - dens%n*v%n
      uc = 0.5D0*indypp*(dy(jp2)*u%n+dy(j+1)*u%nn)
      vc = (dx(i+1)*v%n +dx(i)  *v%en )*indypp *incrxp
      wc = indxp*incrypp*(dx(i+1)*dy(jp2)*(w%n  + w%nb  ) + &
                                 dx(i)  *dy(jp2)*(w%en + w%enb ) + &
                                 dx(i+1)*dy(j+1)*(w%nn + w%nnb ) + &
                                 dx(i)  *dy(j+1)*(w%enn+ w%ennb))*indz(k)
      gradc = indxp  * (dens%en *v%en  - dens%n *v%n )
      grade = indxpp * (dens%een*v%een - dens%en*v%en)
      gradw = indxm  * (dens%n  *v%n   - dens%wn*v%wn)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=nx-1 .AND. flag(ipjpk)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxp*ABS(uc))*deltar*lim
      fn = fn - 0.5D0*(vc-ABS(vc))*s
      ft = ft - 0.25D0*ABS(vc-ABS(vc))*(wc+ABS(wc))*s
!
! n-wn
      deltar  = dens%n*v%n - dens%wn*v%wn
      uc = 0.5D0*indypp*(dy(jp2)*u%wn+dy(j+1)*u%wnn)
      vc = (dx(i-1)*v%n +dx(i)  *v%wn )*indypp *incrxm
      wc = indxm*incrypp*(dx(i-1)*dy(j+1)*(w%n  + w%nb  ) + &
                                dx(i)  *dy(j+1)*(w%wn + w%wnb ) + &
                                dx(i-1)*dy(j)  *(w%nn + w%nnb ) + &
                                dx(i)  *dy(j)  *(w%wnn+ w%wnnb))*indz(k)
      gradc = indxm  * (dens%n *v%n  - dens%wn *v%wn )
      grade = indxp  * (dens%en*v%en - dens%n  *v%n  )
      gradw = indxmm * (dens%wn*v%wn - dens%wwn*v%wwn)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=2 .AND. flag(imjpk)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxm*ABS(uc))*deltar*lim
      fn = fn + 0.5D0*(vc-ABS(vc))*s
      ft = ft + 0.25D0*ABS(vc-ABS(vc))*(wc+ABS(wc))*s
!
! es-s
      deltar  = dens%es*v%es - dens%s*v%s
      uc = 0.5D0*indym*(dy(j-1)*u%c+dy(j)*u%s)
      vc = (dx(i+1)*v%s+dx(i)*v%es)*indym  *incrxp
      wc = indxp*incrym*(dx(i+1)*dy(j-1)*(w%c + w%b  ) + &
                                dx(i)  *dy(j-1)*(w%e + w%eb ) + &
                                dx(i+1)*dy(j)  *(w%s + w%sb ) + &
                                dx(i)  *dy(j)  *(w%es+ w%esb))*indz(k)
      gradc = indxp  * (dens%es *v%es  - dens%s *v%s )
      grade = indxpp * (dens%ees*v%ees - dens%es*v%es)
      gradw = indxm  * (dens%s  *v%s   - dens%ws*v%ws)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=nx-1 .AND. flag(ipjmk)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxp*ABS(uc))*deltar*lim
      ft = ft - 0.25D0*ABS(vc+ABS(vc))*(wc+ABS(wc))*s
!
! s-ws
      deltar  = dens%s*v%s - dens%ws*v%ws
      uc = 0.5D0*indym*(dy(j-1)*u%w+dy(j)*u%ws)
      vc = (dx(i-1)*v%s+dx(i)*v%ws)*indym  *incrxm
      wc = indxm*incrym*(dx(i-1)*dy(j-1)*(w%c + w%b  ) + &
                                dx(i)  *dy(j-1)*(w%w + w%wb ) + &
                                dx(i-1)*dy(j)  *(w%s + w%sb ) + &
                                dx(i)  *dy(j)  *(w%ws+ w%wsb))*indz(k)
      gradc = indxm  * (dens%s *v%s  - dens%ws *v%ws )
      grade = indxp  * (dens%es*v%es - dens%s  *v%s  )
      gradw = indxmm * (dens%ws*v%ws - dens%wws*v%wws)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=2 .AND. flag(imjmk)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxm*ABS(uc))*deltar*lim
      ft = ft + 0.25D0*ABS(vc+ABS(vc))*(wc+ABS(wc))*s
!
! et-t
      deltar  = dens%et*v%et - dens%t*v%t
      uc = 0.5D0*indyp*(dy(j+1)*u%t+dy(j)*u%nt)
      vc = (dx(i+1)*v%t +dx(i)*v%et )*indyp    *incrxp
      wc = indxp*incryp*(dx(i+1)*dy(j+1)*(w%c + w%t  ) + &
                                dx(i)  *dy(j+1)*(w%e + w%et ) + &
                                dx(i+1)*dy(j)  *(w%n + w%nt ) + &
                                dx(i)  *dy(j)  *(w%en+ w%ent))*indz(k+1)
      gradc = indxp  * (dens%et *v%et  - dens%t *v%t )
      grade = indxpp * (dens%eet*v%eet - dens%et*v%et)
      gradw = indxm  * (dens%t  *v%t   - dens%wt*v%wt)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=nx-1 .AND. flag(ipjkp)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxp*ABS(uc))*deltar*lim
      ft = ft - 0.5D0*(1-ABS(vc))*(wc-ABS(wc))*s
!
! t-wt
      deltar  = dens%t*v%t - dens%wt*v%wt
      uc = 0.5D0*indyp*(dy(j+1)*u%wt+dy(j)*u%wnt)
      vc = (dx(i-1)*v%t +dx(i)*v%wt )*indyp    *incrxm
      wc = indxm*incryp*(dx(i-1)*dy(j+1)*(w%c + w%t  ) + &
                                dx(i)  *dy(j+1)*(w%w + w%wt ) + &
                                dx(i-1)*dy(j)  *(w%n + w%nt ) + &
                                dx(i)  *dy(j)  *(w%wn+ w%wnt))*indz(k+1)
      gradc = indxm  * (dens%t *v%t  - dens%wt *v%wt )
      grade = indxp  * (dens%et*v%et - dens%t  *v%t  )
      gradw = indxmm * (dens%wt*v%wt - dens%wwt*v%wwt)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=2 .AND. flag(imjkp)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxm*ABS(uc))*deltar*lim
      ft = ft + 0.5D0*(1-ABS(vc))*(wc-ABS(wc))*s
!
! ent-nt
      deltar  = dens%ent*v%ent - dens%nt*v%nt
      uc = 0.5D0*indypp*(dy(jp2)*u%nt+dy(j+1)*u%nnt)
      vc = (dx(i+1)*v%nt +dx(i)*v%ent )*indypp   *incrxp
      wc = indxp*incrypp*(dx(i+1)*dy(jp2)*(w%n  + w%nt  ) + &
                                 dx(i)  *dy(jp2)*(w%en + w%ent ) + &
                                 dx(i+1)*dy(j+1)*(w%nn + w%nnt ) + &
                                 dx(i)  *dy(j+1)*(w%enn+ w%ennt))*indz(k+1)
      gradc = indxp  * (dens%ent *v%ent  - dens%nt *v%nt )
      grade = indxpp * (dens%eent*v%eent - dens%ent*v%ent)
      gradw = indxm  * (dens%nt  *v%nt   - dens%wnt*v%wnt)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=nx-1 .AND. flag(ipjpkp)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxp*ABS(uc))*deltar*lim
      ft = ft - 0.25D0*ABS(vc-ABS(vc))*(wc-ABS(wc))*s
!
! nt-wnt
      deltar  = dens%nt*v%nt - dens%wnt*v%wnt
      uc = 0.5D0*indypp*(dy(jp2)*u%wnt+dy(j+1)*u%wnnt)
      vc = (dx(i-1)*v%nt +dx(i)  *v%wnt )*indypp   *incrxm
      wc = indxm*incrypp*(dx(i-1)*dy(jp2)*(w%n  + w%nt  ) + &
                                 dx(i)  *dy(jp2)*(w%wn + w%wnt ) + &
                                 dx(i-1)*dy(j+1)*(w%nn + w%nnt ) + &
                                 dx(i)  *dy(j+1)*(w%wnn+ w%wnnt))*indz(k+1)
      gradc = indxm  * (dens%nt *v%nt  - dens%wnt *v%wnt )
      grade = indxp  * (dens%ent*v%ent - dens%nt  *v%nt  )
      gradw = indxmm * (dens%wnt*v%wnt - dens%wwnt*v%wwnt)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=2 .AND. flag(imjpkp)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxm*ABS(uc))*deltar*lim
      ft = ft + 0.25D0*ABS(vc-ABS(vc))*(wc-ABS(wc))*s
!
! est-st
      deltar  = dens%est*v%est - dens%st*v%st
      uc = 0.5D0*indym*(dy(j-1)*u%t+dy(j)*u%st)
      vc = (dx(i+1)*v%st+dx(i)*v%est)*indym    *incrxp
      wc = indxp*incrym*(dx(i+1)*dy(j-1)*(w%c + w%t  ) + &
                                dx(i)  *dy(j-1)*(w%e + w%et ) + &
                                dx(i+1)*dy(j)  *(w%s + w%st ) + &
                                dx(i)  *dy(j)  *(w%es+ w%est))*indz(k+1)
      gradc = indxp  * (dens%est *v%est  - dens%st *v%st )
      grade = indxpp * (dens%eest*v%eest - dens%est*v%est)
      gradw = indxm  * (dens%st  *v%st   - dens%wst*v%wst)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=nx-1 .AND. flag(ipjmkp)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxp*ABS(uc))*deltar*lim
      ft = ft - 0.25D0*ABS(vc+ABS(vc))*(wc-ABS(wc))*s
!
! st-wst
      deltar  = dens%st*v%st - dens%wst*v%wst
      uc = 0.5D0*indym*(dy(j-1)*u%wt+dy(j)*u%wst)
      vc = (dx(i-1)*v%st+dx(i)*v%wst)*indym    *incrxm
      wc = indxm*incrym*(dx(i-1)*dy(j-1)*(w%c + w%t  ) + &
                                dx(i)  *dy(j-1)*(w%w + w%wt ) + &
                                dx(i-1)*dy(j)  *(w%s + w%st ) + &
                                dx(i)  *dy(j)  *(w%ws+ w%wst))*indz(k+1)
      gradc = indxm  * (dens%st *v%st  - dens%wst *v%wst )
      grade = indxp  * (dens%est*v%est - dens%st  *v%st  )
      gradw = indxmm * (dens%wst*v%wst - dens%wwst*v%wwst)
      lim  = 0.D0
      erre = 0.D0
      IF (uc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradw/gradc
      ELSE
        IF (gradc /= 0.D0) erre = grade/gradc
      END IF
      IF (i/=2 .AND. flag(imjmkp)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(uc)*(1-dt*indxm*ABS(uc))*deltar*lim
      ft = ft + 0.25D0*ABS(vc+ABS(vc))*(wc-ABS(wc))*s
!
! c-s
      deltar  = dens%c*v%c - dens%s*v%s
      uc = 0.5D0*(u%c+u%w)*dt*indx(i)
      vc = 0.5D0*(v%c+v%s)
      wc = 0.5D0*(w%c+w%b)*dt*indz(k)
      grads = (v%s *dens%s  - v%ss*dens%ss) * indy(j-1)
      gradc = (v%c *dens%c  - v%s *dens%s ) * indy(j)
      gradn = (v%n *dens%n  - v%c *dens%c ) * indy(j+1)
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = grads / gradc
      ELSE 
        IF (gradc /= 0.D0) erre = gradn / gradc
      END IF 
      IF (j/=2 .AND. flag(ijmk)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indy(j)*ABS(vc))*deltar*lim
      ft = ft + 0.5D0*(wc+ABS(wc))*s
      fe = fe + 0.5D0*(1-ABS(wc))*(uc+ABS(uc))*s
!
! n-c
      deltar  = dens%n*v%n - dens%c*v%c
      uc = 0.5D0*(u%n+u%wn)*dt*indx(i)
      vc = 0.5D0*(v%c+v%n )
      wc = 0.5D0*(w%n+w%nb)*dt*indz(k)
      grads = (v%c *dens%c  - v%s*dens%s) * indy(j)
      gradc = (v%n *dens%n  - v%c*dens%c) * indy(j+1)
      gradn = (v%nn*dens%nn - v%n*dens%n) * indy(jp2)
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = grads / gradc
      ELSE 
        IF (gradc /= 0.D0) erre = gradn / gradc
      END IF 
      IF (j/=ny-1 .AND. flag(ijpk)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indy(j+1)*ABS(vc))*deltar*lim
      ft = ft - 0.5D0*(wc+ABS(wc))*s
      fe = fe - 0.5D0*(1-ABS(wc))*(uc+ABS(uc))*s
!
! en-e
      deltar  = dens%en*v%en - dens%e*v%e
      uc = 0.5D0*(u%n +u%en )*dt*indx(i+1)
      vc = 0.5D0*(v%e +v%en )
      wc = 0.5D0*(w%en+w%enb)*dt*indz(k)
      grads = (v%e  *dens%e   - v%es*dens%es) * indy(j)
      gradc = (v%en *dens%en  - v%e *dens%e ) * indy(j+1)
      gradn = (v%enn*dens%enn - v%en*dens%en) * indy(jp2)
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = grads / gradc
      ELSE 
        IF (gradc /= 0.D0) erre = gradn / gradc
      END IF 
      IF (j/=ny-1 .AND. flag(ipjpk)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indy(j+1)*ABS(vc))*deltar*lim
      fe = fe - 0.5D0*(1-ABS(wc))*(uc-ABS(uc))*s
!
! e-es
      deltar  = dens%e*v%e - dens%es*v%es
      uc = 0.5D0*(u%c+u%e )*dt*indx(i+1)
      vc = 0.5D0*(v%e+v%es)
      wc = 0.5D0*(w%e+w%eb)*dt*indz(k)
      grads = (v%es *dens%es  - v%ess*dens%ess) * indy(j-1)
      gradc = (v%e  *dens%e   - v%es *dens%es ) * indy(j)
      gradn = (v%en *dens%en  - v%e  *dens%e  ) * indy(j+1)
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = grads / gradc
      ELSE 
        IF (gradc /= 0.D0) erre = gradn / gradc
      END IF 
      IF (j/=2 .AND. flag(ipjmk)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indy(j)*ABS(vc))*deltar*lim
      fe = fe + 0.5D0*(1-ABS(wc))*(uc-ABS(uc))*s
!
! nt-t
      deltar  = dens%nt*v%nt - dens%t*v%t
      uc = 0.5D0*(u%nt+u%wnt)*dt*indx(i)
      vc = 0.5D0*(v%t+v%nt  )
      wc = 0.5D0*(w%nt+w%n  )*dt*indz(k+1)
      grads = (v%t  *dens%t   - v%st*dens%st) * indy(j)
      gradc = (v%nt *dens%nt  - v%t *dens%t ) * indy(j+1)
      gradn = (v%nnt*dens%nnt - v%nt*dens%nt) * indy(jp2)
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = grads / gradc
      ELSE 
        IF (gradc /= 0.D0) erre = gradn / gradc
      END IF 
      IF (j/=ny-1 .AND. flag(ijpkp)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indy(j+1)*ABS(vc))*deltar*lim
      ft = ft - 0.5D0*(wc-ABS(wc))*s
      fe = fe - 0.25D0*ABS(wc-ABS(wc))*(uc+ABS(uc))*s
!
! t-st
      deltar  = dens%t*v%t - dens%st*v%st
      uc = 0.5D0*(u%t+u%wt)*dt*indx(i)
      vc = 0.5D0*(v%t+v%st)
      wc = 0.5D0*(w%c+w%t )*dt*indz(k+1)
      grads = (v%st *dens%st  - v%sst*dens%sst) * indy(j-1)
      gradc = (v%t  *dens%t   - v%st *dens%st ) * indy(j)
      gradn = (v%nt *dens%nt  - v%t  *dens%t  ) * indy(j+1)
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = grads / gradc
      ELSE 
        IF (gradc /= 0.D0) erre = gradn / gradc
      END IF 
      IF (j/=2 .AND. flag(ijmkp)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indy(j)*ABS(vc))*deltar*lim
      ft = ft + 0.5D0*(wc-ABS(wc))*s
      fe = fe + 0.25D0*ABS(wc-ABS(wc))*(uc+ABS(uc))*s
!
! nb-b
      deltar  = dens%nb*v%nb - dens%b*v%b
      IF (k <= 2) deltar = 0.D0
      uc = 0.5D0*(u%nb+u%wnb)*dt*indx(i)
      vc = 0.5D0*(v%nb+v%b  )
      wc = 0.5D0*(w%nb+w%nbb)*dt*indz(k-1)
      grads = (v%b  *dens%b   - v%sb*dens%sb) * indy(j)
      gradc = (v%nb *dens%nb  - v%b *dens%b ) * indy(j+1)
      gradn = (v%nnb*dens%nnb - v%nb*dens%nb) * indy(jp2)
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = grads / gradc
      ELSE 
        IF (gradc /= 0.D0) erre = gradn / gradc
      END IF 
      IF (j/=ny-1 .AND. flag(ijpkm)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indy(j+1)*ABS(vc))*deltar*lim
      fe = fe - 0.25D0*ABS(wc+ABS(wc))*(uc+ABS(uc))*s
!
! n-sb
      deltar  = dens%b*v%b - dens%sb*v%sb
      IF (k <= 2) deltar = 0.D0
      uc = 0.5D0*(u%b+u%wb)*dt*indx(i)
      vc = 0.5D0*(v%b+v%sb)
      wc = 0.5D0*(w%b+w%bb)*dt*indz(k-1)
      grads = (v%sb *dens%sb  - v%ssb*dens%ssb) * indy(j-1)
      gradc = (v%b  *dens%b   - v%sb *dens%sb ) * indy(j)
      gradn = (v%nb *dens%nb  - v%b  *dens%b  ) * indy(j+1)
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = grads / gradc
      ELSE 
        IF (gradc /= 0.D0) erre = gradn / gradc
      END IF 
      IF (j/=2 .AND. flag(ijmkm)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indy(j)*ABS(vc))*deltar*lim
      fe = fe + 0.25D0*ABS(wc+ABS(wc))*(uc+ABS(uc))*s
!
! ent-et
      deltar  = dens%ent*v%ent - dens%et*v%et
      uc = 0.5D0*(u%ent+u%nt)*dt*indx(i+1)
      vc = 0.5D0*(v%ent+v%et)
      wc = 0.5D0*(w%ent+w%en)*dt*indz(k+1)
      grads = (v%et  *dens%et   - v%est*dens%est) * indy(j)
      gradc = (v%ent *dens%ent  - v%et *dens%et ) * indy(j+1)
      gradn = (v%ennt*dens%ennt - v%ent*dens%ent) * indy(jp2)
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = grads / gradc
      ELSE 
        IF (gradc /= 0.D0) erre = gradn / gradc
      END IF 
      IF (j/=ny-1 .AND. flag(ipjpkp)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indy(j+1)*ABS(vc))*deltar*lim
      fe = fe - 0.25D0*ABS(wc-ABS(wc))*(uc-ABS(uc))*s
!
! et-est
      deltar  = dens%et*v%et - dens%est*v%est
      uc = 0.5D0*(u%et+u%t  )*dt*indx(i+1)
      vc = 0.5D0*(v%et+v%est)
      wc = 0.5D0*(w%et+w%e  )*dt*indz(k+1)
      grads = (v%est *dens%est  - v%esst*dens%esst) * indy(j-1)
      gradc = (v%et  *dens%et   - v%est *dens%est ) * indy(j)
      gradn = (v%ent *dens%ent  - v%et  *dens%et  ) * indy(j+1)
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = grads / gradc
      ELSE 
        IF (gradc /= 0.D0) erre = gradn / gradc
      END IF 
      IF (j/=2 .AND. flag(ipjmkp)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indy(j)*ABS(vc))*deltar*lim
      fe = fe + 0.25D0*ABS(wc-ABS(wc))*(uc-ABS(uc))*s
!
! enb-eb
      deltar  = dens%enb*v%enb - dens%eb*v%eb
      IF (k <= 2) deltar = 0.D0
      uc = 0.5D0*(u%enb+u%nb  )*dt*indx(i+1)
      vc = 0.5D0*(v%enb+v%eb  )
      wc = 0.5D0*(w%enb+w%enbb)*dt*indz(k-1)
      grads = (v%eb  *dens%eb   - v%esb*dens%esb) * indy(j)
      gradc = (v%enb *dens%enb  - v%eb *dens%eb ) * indy(j+1)
      gradn = (v%ennb*dens%ennb - v%enb*dens%enb) * indy(jp2)
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = grads / gradc
      ELSE 
        IF (gradc /= 0.D0) erre = gradn / gradc
      END IF 
      IF (j/=ny-1 .AND. flag(ipjpkm)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indy(j+1)*ABS(vc))*deltar*lim
      fe = fe - 0.25D0*ABS(wc+ABS(wc))*(uc-ABS(uc))*s
!
! eb-esb
      deltar  = dens%eb*v%eb - dens%esb*v%esb
      IF (k <= 2) deltar = 0.D0
      uc = 0.5D0*(u%eb+u%b  )*dt*indx(i+1)
      vc = 0.5D0*(v%eb+v%esb)
      wc = 0.5D0*(w%eb+w%ebb)*dt*indz(k-1)
      grads = (v%esb *dens%esb  - v%essb*dens%essb) * indy(j-1)
      gradc = (v%eb  *dens%eb   - v%esb *dens%esb ) * indy(j)
      gradn = (v%enb *dens%enb  - v%eb  *dens%eb  ) * indy(j+1)
      lim  = 0.D0
      erre = 0.D0
      IF (vc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = grads / gradc
      ELSE 
        IF (gradc /= 0.D0) erre = gradn / gradc
      END IF 
      IF (j/=2 .AND. flag(ipjmkm)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(vc)*(1-dt*indy(j)*ABS(vc))*deltar*lim
      fe = fe + 0.25D0*ABS(wc+ABS(wc))*(uc-ABS(uc))*s
!
! c-b
      deltar  = dens%c*v%c - dens%b*v%b
      uc = indzm*incryp*(dz(k-1)*dy(j+1)*(u%c + u%w   ) + &
                                dz(k)  *dy(j+1)*(u%b + u%wb  ) + &
                                dz(k-1)*dy(j)  *(u%n + u%wn  ) + &
                                dz(k)  *dy(j)  *(u%nb+ u%wnb))*indx(i)
      vc = (dz(k-1)*v%c+dz(k)*v%b )*indyp  *incrzm
      wc = 0.5D0*indyp*(dy(j+1)*w%b+dy(j)*w%nb)
      gradc = indzm  * (dens%c*v%c - dens%b *v%b )
      gradt = indzp  * (dens%t*v%t - dens%c *v%c )
      gradb = indzmm * (dens%b*v%b - dens%bb*v%bb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=2 .AND. flag(ijkm)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indzm*ABS(wc))*deltar*lim
      fe = fe + 0.5D0*(uc+ABS(uc))*s
      fn = fn + 0.5D0*(1-ABS(uc))*(vc+ABS(vc))*s
!
! t-c
      deltar  = dens%t*v%t - dens%c*v%c
      uc = indzp*incryp*(dz(k+1)*dy(j+1)*(u%c + u%w   ) + &
                                dz(k)  *dy(j+1)*(u%t + u%wt  ) + &
                                dz(k+1)*dy(j)  *(u%n + u%wn  ) + &
                                dz(k)  *dy(j)  *(u%nt+ u%wnt))*indx(i)
      vc = (dz(k+1)*v%c+dz(k)*v%t )*indyp  *incrzp
      wc = 0.5D0*indyp*(dy(j+1)*w%c+dy(j)*w%n)
      gradc = indzp  * (dens%t *v%t  - dens%c*v%c)
      gradt = indzpp * (dens%tt*v%tt - dens%t*v%t)
      gradb = indzm  * (dens%c *v%c  - dens%b*v%b)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(ijkp)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indzp*ABS(wc))*deltar*lim
      fe = fe - 0.5D0*(uc+ABS(uc))*s
      fn = fn - 0.5D0*(1-ABS(uc))*(vc+ABS(vc))*s
!
! wt-w
      deltar  = dens%wt*v%wt - dens%w*v%w
      uc = indzp*incryp*(dz(k+1)*dy(j+1)*(u%w  + u%ww   ) + &
                                dz(k)  *dy(j+1)*(u%wt + u%wwt  ) + &
                                dz(k+1)*dy(j)  *(u%wn + u%wwn  ) + &
                                dz(k)  *dy(j)  *(u%wnt+ u%wwnt))*indx(i-1)
      vc = (dz(k+1)*v%w +dz(k)*v%wt )*indyp    *incrzp
      wc = 0.5D0*indyp*(dy(j+1)*w%w+dy(j)*w%wn)
      gradc = indzp  * (dens%wt *v%wt  - dens%w *v%w )
      gradt = indzpp * (dens%wtt*v%wtt - dens%wt*v%wt)
      gradb = indzm  * (dens%w  *v%w   - dens%wb*v%wb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(imjkp)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indzp*ABS(wc))*deltar*lim
      fn = fn - 0.25D0*ABS(uc+ABS(uc))*(vc+ABS(vc))*s
!
! w-wb
      deltar  = dens%w*v%w - dens%wb*v%wb
      uc = indzm*incryp*(dz(k-1)*dy(j+1)*(u%w  + u%ww   ) + &
                                dz(k)  *dy(j+1)*(u%wb + u%wwb  ) + &
                                dz(k-1)*dy(j)  *(u%wn + u%wwn  ) + &
                                dz(k)  *dy(j)  *(u%wnb+ u%wwnb))*indx(i-1)
      vc = (dz(k-1)*v%w +dz(k)*v%wb )*indyp    *incrzm
      wc = 0.5D0*indyp*(dy(j+1)*w%wb+dy(j)*w%wnb)
      gradc = indzm  * (dens%w *v%w  - dens%wb *v%wb )
      gradt = indzp  * (dens%wt*v%wt - dens%w  *v%w  )
      gradb = indzmm * (dens%wb*v%wb - dens%wbb*v%wbb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=2 .AND. flag(imjkm)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indzm*ABS(wc))*deltar*lim
      fn = fn + 0.25D0*ABS(uc+ABS(uc))*(vc+ABS(vc))*s
!
! et-e
      deltar  = dens%et*v%et - dens%e*v%e
      uc = indzp*incryp*(dz(k+1)*dy(j+1)*(u%c + u%e   ) + &
                                dz(k)  *dy(j+1)*(u%t + u%et  ) + &
                                dz(k+1)*dy(j)  *(u%n + u%en  ) + &
                                dz(k)  *dy(j)  *(u%nt+ u%ent))*indx(i+1)
      vc = (dz(k+1)*v%e+dz(k)*v%et)*indyp*incrzp
      wc = 0.5D0*indyp*(dy(j+1)*w%e+dy(j)*w%en)
      gradc = indzp  * (dens%et *v%et  - dens%e *v%e )
      gradt = indzpp * (dens%ett*v%ett - dens%et*v%et)
      gradb = indzm  * (dens%e  *v%e   - dens%eb*v%eb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(ipjkp)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indzp*ABS(wc))*deltar*lim
      fe = fe - 0.5D0*(uc-ABS(uc))*s
      fn = fn - 0.25D0*ABS(uc-ABS(uc))*(vc+ABS(vc))*s
!
! e-eb
      deltar  = dens%e*v%e - dens%eb*v%eb
      uc = indzm*incryp*(dz(k-1)*dy(j+1)*(u%c + u%e   ) + &
                                dz(k)  *dy(j+1)*(u%b + u%eb  ) + &
                                dz(k-1)*dy(j)  *(u%n + u%en  ) + &
                                dz(k)  *dy(j)  *(u%nb+ u%enb))*indx(i+1)
      vc = (dz(k-1)*v%e +dz(k)*v%eb) *indyp    *incrzm
      wc = 0.5D0*indyp*(dy(j+1)*w%eb+dy(j)*w%enb)
      gradc = indzm  * (dens%e *v%e  - dens%eb *v%eb )
      gradt = indzp  * (dens%et*v%et - dens%e  *v%e  )
      gradb = indzmm * (dens%eb*v%eb - dens%ebb*v%ebb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=2 .AND. flag(ipjkm)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indzm*ABS(wc))*deltar*lim
      fe = fe + 0.5D0*(uc-ABS(uc))*s
      fn = fn + 0.25D0*ABS(uc-ABS(uc))*(vc+ABS(vc))*s
!
! nt-n
      deltar  = dens%nt*v%nt - dens%n*v%n
      uc = indzp*incrypp*(dz(k+1)*dy(jp2)*(u%n  + u%wn  ) + &
                                 dz(k)  *dy(jp2)*(u%nt + u%wnt ) + &
                                 dz(k+1)*dy(j+1)*(u%nn + u%wnn ) + &
                                 dz(k)  *dy(j+1)*(u%nnt+ u%wnnt))*indx(i)
      vc = (dz(k+1)*v%n +dz(k)  *v%nt )*indypp *incrzp
      wc = 0.5D0*indypp*(dy(jp2)*w%n+dy(j+1)*w%nn)
      gradc = indzp  * (dens%nt *v%nt  - dens%n *v%n )
      gradt = indzpp * (dens%ntt*v%ntt - dens%nt*v%nt)
      gradb = indzm  * (dens%n  *v%n   - dens%nb*v%nb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(ijpkp)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indzp*ABS(wc))*deltar*lim
      fn = fn - 0.5D0*(1-ABS(uc))*(vc-ABS(vc))*s
!
! n-nb
      deltar  = dens%n*v%n - dens%nb*v%nb
      uc = indzm*incrypp*(dz(k-1)*dy(jp2)*(u%n  + u%wn  ) + &
                                 dz(k)  *dy(jp2)*(u%nb + u%wnb ) + &
                                 dz(k-1)*dy(j+1)*(u%nn + u%wnn ) + &
                                 dz(k)  *dy(j+1)*(u%nnb+ u%wnnb))*indx(i)
      vc = (dz(k-1)*v%n +dz(k)  *v%nb )*indypp *incrzm
      wc = 0.5D0*indypp*(dy(jp2)*w%nb+dy(j+1)*w%nnb)
      gradc = indzm  * (dens%n *v%n  - dens%nb *v%nb )
      gradt = indzp  * (dens%nt*v%nt - dens%n  *v%n  )
      gradb = indzmm * (dens%nb*v%nb - dens%nbb*v%nbb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=2 .AND. flag(ijpkm)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indzm*ABS(wc))*deltar*lim
      fn = fn + 0.5D0*(1-ABS(uc))*(vc-ABS(vc))*s
!
! wnt-wn
      deltar  = dens%wnt*v%wnt - dens%wn*v%wn
      uc = indzp*incrypp*(dz(k+1)*dy(jp2)*(u%wn  + u%wwn  ) + &
                                 dz(k)  *dy(jp2)*(u%wnt + u%wwnt ) + &
                                 dz(k+1)*dy(j+1)*(u%wnn + u%wwnn ) + &
                                 dz(k)  *dy(j+1)*(u%wnnt+ u%wwnnt))*indx(i-1)
      vc = (dz(k+1)*v%wn +dz(k)  *v%wnt )*indypp   *incrzp
      wc = 0.5D0*indypp*(dy(jp2)*w%wn+dy(j+1)*w%wnn)
      gradc = indzp  * (dens%wnt *v%wnt  - dens%wn *v%wn )
      gradt = indzpp * (dens%wntt*v%wntt - dens%wnt*v%wnt)
      gradb = indzm  * (dens%wn  *v%wn   - dens%wnb*v%wnb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(imjpkp)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indzp*ABS(wc))*deltar*lim
      fn = fn - 0.25D0*ABS(uc+ABS(uc))*(vc-ABS(vc))*s
!
! wn-wnb
      deltar  = dens%wn*v%wn - dens%wnb*v%wnb
      uc = indzm*incrypp*(dz(k-1)*dy(jp2)*(u%wn  + u%wwn  ) + &
                                 dz(k)  *dy(jp2)*(u%wnb + u%wwnb ) + &
                                 dz(k-1)*dy(j+1)*(u%wnn + u%wwnn ) + &
                                 dz(k)  *dy(j+1)*(u%wnnb+ u%wwnnb))*indx(i-1)
      vc = (dz(k-1)*v%wn +dz(k)  *v%wnb )*indypp   *incrzm
      wc = 0.5D0*indypp*(dy(jp2)*w%wnb+dy(j+1)*w%wnnb)
      gradc = indzm  * (dens%wn *v%wn  - dens%wnb *v%wnb )
      gradt = indzp  * (dens%wnt*v%wnt - dens%wn  *v%wn  )
      gradb = indzmm * (dens%wnb*v%wnb - dens%wnbb*v%wnbb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=2 .AND. flag(imjpkm)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indzm*ABS(wc))*deltar*lim
      fn = fn + 0.25D0*ABS(uc+ABS(uc))*(vc-ABS(vc))*s
!
! ent-en
      deltar  = dens%ent*v%ent - dens%en*v%en
      uc = indzp*incrypp*(dz(k+1)*dy(jp2)*(u%n  + u%en  ) + &
                                 dz(k)  *dy(jp2)*(u%nt + u%ent ) + &
                                 dz(k+1)*dy(j+1)*(u%nn + u%enn ) + &
                                 dz(k)  *dy(j+1)*(u%nnt+ u%ennt))*indx(i+1)
      vc = (dz(k+1)*v%en +dz(k)  *v%ent )*indypp   *incrzp
      wc = 0.5D0*indypp*(dy(jp2)*w%en+dy(j+1)*w%enn)
      gradc = indzp  * (dens%ent *v%ent  - dens%en *v%en )
      gradt = indzpp * (dens%entt*v%entt - dens%ent*v%ent)
      gradb = indzm  * (dens%en  *v%en   - dens%enb*v%enb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=nz-1 .AND. flag(ipjpkp)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indzp*ABS(wc))*deltar*lim
      fn = fn - 0.25D0*ABS(uc-ABS(uc))*(vc-ABS(vc))*s
!
! en-enb
      deltar  = dens%en*v%en - dens%enb*v%enb
      uc = indzm*incrypp*(dz(k-1)*dy(jp2)*(u%n  + u%en  ) + &
                                 dz(k)  *dy(jp2)*(u%nb + u%enb ) + &
                                 dz(k-1)*dy(j+1)*(u%nn + u%enn ) + &
                                 dz(k)  *dy(j+1)*(u%nnb+ u%ennb))*indx(i+1)
      vc = (dz(k-1)*v%en +dz(k)  *v%enb) *indypp   *incrzm
      wc = 0.5D0*indypp*(dy(jp2)*w%enb+dy(j+1)*w%ennb)
      gradc = indzm  * (dens%en *v%en  - dens%enb *v%enb )
      gradt = indzp  * (dens%ent*v%ent - dens%en  *v%en  )
      gradb = indzmm * (dens%enb*v%enb - dens%enbb*v%enbb)
      lim  = 0.D0
      erre = 0.D0
      IF (wc >= 0.D0) THEN
        IF (gradc /= 0.D0) erre = gradb/gradc
      ELSE
        IF (gradc /= 0.D0) erre = gradt/gradc
      END IF
      IF (k/=2 .AND. flag(ipjpkm)==fluid) CALL limiters(lv,lim,erre)
      s = 0.5D0*ABS(wc)*(1-dt*indzm*ABS(wc))*deltar*lim
      fn = fn + 0.25D0*ABS(uc-ABS(uc))*(vc-ABS(vc))*s
!
    RETURN
    END SUBROUTINE ctu3_flv_3d
!----------------------------------------------------------------------
      END MODULE convective_fluxes_v
!-----------------------------------------------------------------------
