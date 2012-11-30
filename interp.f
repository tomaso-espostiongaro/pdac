!----------------------------------------------------------------------
      MODULE interpolate_fields
!----------------------------------------------------------------------
!
      USE dimensions
      USE grid, ONLY: dx, dy, dz

      IMPLICIT NONE
!
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE interpolate_x(field, staggered_field, i)
      USE set_indexes, ONLY: stencil
!
! ... interpolate density on the mesh staggered along x
!
      IMPLICIT NONE

      TYPE(stencil), INTENT(IN) :: field
      TYPE(stencil), INTENT(OUT) :: staggered_field
      INTEGER, INTENT(IN) :: i

      INTEGER :: ip2, ip3, im2

      REAL*8 :: dxm, dxp, dxpp, dxmm, dxppp, indxpp, indxp, indxm, indxmm, indxppp 
!
      ip2 = MIN( nx, i+2 )
      ip3 = MIN( nx, i+3 )
      im2 = MAX(  1, i-2 )

      dxm   = dx(i)   + dx(i-1)
      dxp   = dx(i)   + dx(i+1)
      dxpp  = dx(i+1) + dx(ip2)
      dxppp = dx(ip2) + dx(ip3)
      dxmm  = dx(i-1) + dx(im2)

      indxm   = 1.D0/dxm
      indxp   = 1.D0/dxp
      indxpp  = 1.D0/dxpp
      indxppp = 1.D0/dxppp
      indxmm  = 1.D0/dxmm
!       
! ... Compute linearly interpolated values of density on the staggered grid
!
      staggered_field%c    = (dx(i+1) * field%c     + dx(i)   * field%e    ) * indxp
      staggered_field%e    = (dx(ip2) * field%e     + dx(i+1) * field%ee   ) * indxpp
      staggered_field%n    = (dx(i+1) * field%n     + dx(i)   * field%en   ) * indxp
      staggered_field%t    = (dx(i+1) * field%t     + dx(i)   * field%et   ) * indxp
      staggered_field%w    = (dx(i)   * field%w     + dx(i-1) * field%c    ) * indxm
      staggered_field%s    = (dx(i+1) * field%s     + dx(i)   * field%es   ) * indxp
      staggered_field%b    = (dx(i+1) * field%b     + dx(i)   * field%eb   ) * indxp
      staggered_field%eb   = (dx(i+1) * field%eeb   + dx(ip2) * field%eb   ) * indxpp
      staggered_field%et   = (dx(i+1) * field%eet   + dx(ip2) * field%et   ) * indxpp
      staggered_field%wt   = (dx(i)   * field%wt    + dx(i-1) * field%t    ) * indxm
      staggered_field%wb   = (dx(i)   * field%wb    + dx(i-1) * field%b    ) * indxm
      staggered_field%nt   = (dx(i)   * field%ent   + dx(i+1) * field%nt   ) * indxp
      staggered_field%nb   = (dx(i)   * field%enb   + dx(i+1) * field%nb   ) * indxp
      staggered_field%st   = (dx(i)   * field%est   + dx(i+1) * field%st   ) * indxp
      staggered_field%sb   = (dx(i)   * field%esb   + dx(i+1) * field%sb   ) * indxp
      staggered_field%en   = (dx(i+1) * field%een   + dx(ip2) * field%en   ) * indxpp
      staggered_field%es   = (dx(i+1) * field%ees   + dx(ip2) * field%es   ) * indxpp
      staggered_field%wn   = (dx(i-1) * field%n     + dx(i)   * field%wn   ) * indxm
      staggered_field%ws   = (dx(i-1) * field%s     + dx(i)   * field%ws   ) * indxm
      staggered_field%wst  = (dx(i-1) * field%st    + dx(i)   * field%wst  ) * indxm
      staggered_field%wnt  = (dx(i-1) * field%nt    + dx(i)   * field%wnt  ) * indxm
      staggered_field%wsb  = (dx(i-1) * field%sb    + dx(i)   * field%wsb  ) * indxm
      staggered_field%wnb  = (dx(i-1) * field%nb    + dx(i)   * field%wnb  ) * indxm
      staggered_field%est  = (dx(i+1) * field%eest  + dx(ip2) * field%est  ) * indxpp
      staggered_field%ent  = (dx(i+1) * field%eent  + dx(ip2) * field%ent  ) * indxpp
      staggered_field%esb  = (dx(i+1) * field%eesb  + dx(ip2) * field%esb  ) * indxpp
      staggered_field%enb  = (dx(i+1) * field%eenb  + dx(ip2) * field%enb  ) * indxpp
      staggered_field%ww   = (dx(im2) * field%w     + dx(i-1) * field%ww   ) * indxmm
      staggered_field%ss   = (dx(i)   * field%ess   + dx(i+1) * field%ss   ) * indxp
      staggered_field%bb   = (dx(i)   * field%ebb   + dx(i+1) * field%bb   ) * indxp
      staggered_field%ee   = (dx(ip2) * field%eee   + dx(ip3) * field%ee   ) * indxppp
      staggered_field%nn   = (dx(i)   * field%enn   + dx(i+1) * field%nn   ) * indxp
      staggered_field%tt   = (dx(i)   * field%ett   + dx(i+1) * field%tt   ) * indxp
      staggered_field%een  = (dx(ip2) * field%eeen  + dx(ip3) * field%een  ) * indxppp
      staggered_field%ees  = (dx(ip2) * field%eees  + dx(ip3) * field%ees  ) * indxppp
      staggered_field%eet  = (dx(ip2) * field%eeet  + dx(ip3) * field%eet  ) * indxppp
      staggered_field%eent = (dx(ip2) * field%eeent + dx(ip3) * field%eent ) * indxppp
      staggered_field%eest = (dx(ip2) * field%eeest + dx(ip3) * field%eest ) * indxppp
      staggered_field%wwn  = (dx(im2) * field%wn    + dx(i-1) * field%wwn  ) * indxmm
      staggered_field%wws  = (dx(im2) * field%ws    + dx(i-1) * field%wws  ) * indxmm
      staggered_field%wwt  = (dx(im2) * field%wt    + dx(i-1) * field%wwt  ) * indxmm
      staggered_field%wwnt = (dx(im2) * field%wnt   + dx(i-1) * field%wwnt ) * indxmm
      staggered_field%wwst = (dx(im2) * field%wst   + dx(i-1) * field%wwst ) * indxmm
      staggered_field%enn  = (dx(i+1) * field%eenn  + dx(ip2) * field%enn  ) * indxpp
      staggered_field%nnt  = (dx(i)   * field%ennt  + dx(i+1) * field%nnt  ) * indxp
      staggered_field%nnb  = (dx(i)   * field%ennb  + dx(i+1) * field%nnb  ) * indxp
      staggered_field%ennt = (dx(i+1) * field%eennt + dx(ip2) * field%ennt ) * indxpp
      staggered_field%ennb = (dx(i+1) * field%eennb + dx(ip2) * field%ennb ) * indxpp
      staggered_field%ess  = (dx(i+1) * field%eess  + dx(ip2) * field%ess  ) * indxpp
      staggered_field%sst  = (dx(i)   * field%esst  + dx(i+1) * field%sst  ) * indxp
      staggered_field%ssb  = (dx(i)   * field%essb  + dx(i+1) * field%ssb  ) * indxp
      staggered_field%esst = (dx(i+1) * field%eesst + dx(ip2) * field%esst ) * indxpp
      staggered_field%essb = (dx(i+1) * field%eessb + dx(ip2) * field%essb ) * indxpp
      staggered_field%wtt  = (dx(i-1) * field%tt    + dx(i)   * field%wtt  ) * indxm
      staggered_field%ett  = (dx(i+1) * field%eett  + dx(ip2) * field%ett  ) * indxpp
      staggered_field%ntt  = (dx(i)   * field%entt  + dx(i+1) * field%ntt  ) * indxp
      staggered_field%wntt = (dx(i-1) * field%ntt   + dx(i)   * field%wntt ) * indxm
      staggered_field%entt = (dx(i+1) * field%eentt + dx(ip2) * field%entt ) * indxpp
      staggered_field%wbb  = (dx(i-1) * field%bb    + dx(i)   * field%wbb  ) * indxm
      staggered_field%ebb  = (dx(i+1) * field%eebb  + dx(ip2) * field%ebb  ) * indxpp
      staggered_field%nbb  = (dx(i)   * field%enbb  + dx(i+1) * field%nbb  ) * indxp
      staggered_field%wnbb = (dx(i-1) * field%nbb   + dx(i)   * field%wnbb ) * indxm
      staggered_field%enbb = (dx(i+1) * field%eenbb + dx(ip2) * field%enbb ) * indxpp

!
      RETURN
      END SUBROUTINE interpolate_x
!----------------------------------------------------------------------
      SUBROUTINE interpolate_y(field, staggered_field, j)
      USE set_indexes, ONLY: stencil
!
! ... interpolate density on the mesh staggered along x
!
      IMPLICIT NONE
!
      TYPE(stencil), INTENT(IN) :: field
      TYPE(stencil), INTENT(OUT) :: staggered_field
      INTEGER, INTENT(IN) :: j

      REAL*8 :: dym, dyp, dypp, dymm, dyppp, indypp, indyp, indym, indymm, indyppp

      INTEGER :: jp2, jp3, jm2
!
      jp2 = MIN( ny, j+2 )
      jp3 = MIN( ny, j+3 )
      jm2 = MAX(  1, j-2 )

      dym   = dy(j)   + dy(j-1)
      dyp   = dy(j)   + dy(j+1)
      dypp  = dy(j+1) + dy(jp2)
      dyppp = dy(jp2) + dy(jp3)
      dymm  = dy(j-1) + dy(jm2)

      indym   = 1.D0/dym
      indyp   = 1.D0/dyp
      indypp  = 1.D0/dypp
      indyppp = 1.D0/dyppp
      indymm  = 1.D0/dymm
!       
! ... Compute linearly interpolated values of density on the staggered grid
!
      staggered_field%c    = (dy(j+1) * field%c     + dy(j)   * field%n ) * indyp
      staggered_field%e    = (dy(j+1) * field%e     + dy(j)   * field%en) * indyp
      staggered_field%n    = (dy(jp2) * field%n     + dy(j+1) * field%nn) * indypp
      staggered_field%t    = (dy(j+1) * field%t     + dy(j)   * field%nt) * indyp
      staggered_field%w    = (dy(j+1) * field%w     + dy(j)   * field%wn) * indyp
      staggered_field%s    = (dy(j)   * field%s     + dy(j-1) * field%c    ) * indym
      staggered_field%b    = (dy(j+1) * field%b     + dy(j)   * field%nb   ) * indyp
      staggered_field%eb   = (dy(j)   * field%enb   + dy(j+1) * field%eb   ) * indyp
      staggered_field%et   = (dy(j)   * field%ent   + dy(j+1) * field%et   ) * indyp
      staggered_field%wb   = (dy(j)   * field%wnb   + dy(j+1) * field%wb   ) * indyp
      staggered_field%wt   = (dy(j)   * field%wnt   + dy(j+1) * field%wt   ) * indyp
      staggered_field%nt   = (dy(j+1) * field%nnt   + dy(jp2) * field%nt   ) * indypp
      staggered_field%nb   = (dy(j+1) * field%nnb   + dy(jp2) * field%nb   ) * indypp
      staggered_field%st   = (dy(j-1) * field%t     + dy(j)   * field%st   ) * indym
      staggered_field%sb   = (dy(j-1) * field%b     + dy(j)   * field%sb   ) * indym
      staggered_field%en   = (dy(j+1) * field%enn   + dy(jp2) * field%en   ) * indypp
      staggered_field%es   = (dy(j-1) * field%e     + dy(j)   * field%es   ) * indym
      staggered_field%wn   = (dy(j+1) * field%wnn   + dy(jp2) * field%wn   ) * indypp
      staggered_field%ws   = (dy(j-1) * field%w     + dy(j)   * field%ws   ) * indym
      staggered_field%wst  = (dy(j-1) * field%wt    + dy(j)   * field%wst  ) * indym
      staggered_field%wsb  = (dy(j-1) * field%wb    + dy(j)   * field%wsb  ) * indym
      staggered_field%wnt  = (dy(j+1) * field%wnnt  + dy(jp2) * field%wnt  ) * indypp
      staggered_field%wnb  = (dy(j+1) * field%wnnb  + dy(jp2) * field%wnb  ) * indypp
      staggered_field%est  = (dy(j-1) * field%et    + dy(j)   * field%est  ) * indym
      staggered_field%esb  = (dy(j-1) * field%eb    + dy(j)   * field%esb  ) * indym
      staggered_field%ent  = (dy(j+1) * field%ennt  + dy(jp2) * field%ent  ) * indypp
      staggered_field%enb  = (dy(j+1) * field%ennb  + dy(jp2) * field%enb  ) * indypp
      staggered_field%ww   = (dy(j)   * field%wwn   + dy(j+1) * field%ww   ) * indyp
      staggered_field%ss   = (dy(jm2) * field%s     + dy(j-1) * field%ss   ) * indymm
      staggered_field%bb   = (dy(j)   * field%nbb   + dy(j+1) * field%bb   ) * indyp
      staggered_field%ee   = (dy(j)   * field%een   + dy(j+1) * field%ee   ) * indyp
      staggered_field%nn   = (dy(jp2) * field%nnn   + dy(jp3) * field%nn   ) * indyppp
      staggered_field%tt   = (dy(j)   * field%ntt   + dy(j+1) * field%tt   ) * indyp
      staggered_field%een  = (dy(j+1) * field%eenn  + dy(jp2) * field%een  ) * indypp
      staggered_field%ees  = (dy(j-1) * field%ee    + dy(j)   * field%ees  ) * indym
      staggered_field%eet  = (dy(j)   * field%eent  + dy(j+1) * field%eet  ) * indyp
      staggered_field%eent = (dy(j+1) * field%eennt + dy(jp2) * field%eent ) * indypp
      staggered_field%eest = (dy(j-1) * field%eet   + dy(j)   * field%eest ) * indym
      staggered_field%wwn  = (dy(j+1) * field%wwnn  + dy(jp2) * field%wwn  ) * indypp
      staggered_field%wws  = (dy(j-1) * field%ww    + dy(j)   * field%wws  ) * indym
      staggered_field%wwt  = (dy(j)   * field%wwnt  + dy(j+1) * field%wwt  ) * indyp
      staggered_field%wwnt = (dy(j+1) * field%wwnnt + dy(jp2) * field%wwnt ) * indypp
      staggered_field%wwst = (dy(j-1) * field%wwt   + dy(j)   * field%wwst ) * indym
      staggered_field%enn  = (dy(jp2) * field%ennn  + dy(jp3) * field%enn  ) * indyppp
      staggered_field%nnt  = (dy(jp2) * field%nnnt  + dy(jp3) * field%nnt  ) * indyppp
      staggered_field%nnb  = (dy(jp2) * field%nnnb  + dy(jp3) * field%nnb  ) * indyppp
      staggered_field%ennt = (dy(jp2) * field%ennnt + dy(jp3) * field%ennt ) * indyppp
      staggered_field%ennb = (dy(jp2) * field%ennnb + dy(jp3) * field%ennb ) * indyppp
      staggered_field%ess  = (dy(jm2) * field%es    + dy(j-1) * field%ess  ) * indymm
      staggered_field%sst  = (dy(jm2) * field%st    + dy(j-1) * field%sst  ) * indymm
      staggered_field%ssb  = (dy(jm2) * field%sb    + dy(j-1) * field%ssb  ) * indymm
      staggered_field%esst = (dy(jm2) * field%est   + dy(j-1) * field%esst ) * indymm
      staggered_field%essb = (dy(jm2) * field%esb   + dy(j-1) * field%essb ) * indymm
      staggered_field%wtt  = (dy(j)   * field%wntt  + dy(j+1) * field%wtt  ) * indyp
      staggered_field%ett  = (dy(j)   * field%entt  + dy(j+1) * field%ett  ) * indyp
      staggered_field%ntt  = (dy(j+1) * field%nntt  + dy(jp2) * field%ntt  ) * indypp
      staggered_field%wntt = (dy(j+1) * field%wnntt + dy(jp2) * field%wntt ) * indypp
      staggered_field%entt = (dy(j+1) * field%enntt + dy(jp2) * field%entt ) * indypp
      staggered_field%wbb  = (dy(j)   * field%wnbb  + dy(j+1) * field%wbb  ) * indyp
      staggered_field%ebb  = (dy(j)   * field%enbb  + dy(j+1) * field%ebb  ) * indyp
      staggered_field%nbb  = (dy(j+1) * field%nnbb  + dy(jp2) * field%nbb  ) * indypp
      staggered_field%wnbb = (dy(j+1) * field%wnnbb + dy(jp2) * field%wnbb ) * indypp
      staggered_field%enbb = (dy(j+1) * field%ennbb + dy(jp2) * field%enbb ) * indypp
!
      RETURN
      END SUBROUTINE interpolate_y
!----------------------------------------------------------------------
      SUBROUTINE interpolate_z(field, staggered_field, k)
      USE set_indexes, ONLY: stencil
!
! ... interpolate density on the mesh staggered along x
!
      IMPLICIT NONE
!
      TYPE(stencil), INTENT(IN) :: field
      TYPE(stencil), INTENT(OUT) :: staggered_field
      INTEGER, INTENT(IN) :: k

      REAL*8 :: dzp, indzp, dzm, indzm, dzpp, indzpp, dzppp, indzppp, dzmm, indzmm

      INTEGER :: kp2, kp3, km2
!
      kp2 = MIN( nz, k+2 )
      kp3 = MIN( nz, k+3 )
      km2 = MAX(  1, k-2 )

      dzm   = dz(k)   + dz(k-1)
      dzp   = dz(k)   + dz(k+1)
      dzpp  = dz(k+1) + dz(kp2)
      dzppp = dz(kp2) + dz(kp3)
      dzmm  = dz(k-1) + dz(km2)

      indzm   = 1.D0/dzm
      indzp   = 1.D0/dzp
      indzpp  = 1.D0/dzpp
      indzppp = 1.D0/dzppp
      indzmm  = 1.D0/dzmm
!       
! ... values of density interpolated linearly on the staggered grid
!
      staggered_field%c    = (dz(k+1) * field%c     + dz(k)   * field%t    ) * indzp
      staggered_field%e    = (dz(k+1) * field%e     + dz(k)   * field%et   ) * indzp
      staggered_field%n    = (dz(k+1) * field%n     + dz(k)   * field%nt   ) * indzp
      staggered_field%t    = (dz(kp2) * field%t     + dz(k+1) * field%tt   ) * indzpp
      staggered_field%w    = (dz(k+1) * field%w     + dz(k)   * field%wt   ) * indzp
      staggered_field%s    = (dz(k+1) * field%s     + dz(k)   * field%st   ) * indzp
      staggered_field%b    = (dz(k)   * field%b     + dz(k-1) * field%c    ) * indzm
      staggered_field%et   = (dz(kp2) * field%et    + dz(k+1) * field%ett  ) * indzpp
      staggered_field%wt   = (dz(kp2) * field%wt    + dz(k+1) * field%wtt  ) * indzpp
      staggered_field%eb   = (dz(k)   * field%eb    + dz(k-1) * field%e    ) * indzm
      staggered_field%wb   = (dz(k)   * field%wb    + dz(k-1) * field%w    ) * indzm
      staggered_field%nt   = (dz(k+1) * field%ntt   + dz(kp2) * field%nt   ) * indzpp
      staggered_field%nb   = (dz(k-1) * field%n     + dz(k)   * field%nb   ) * indzm
      staggered_field%st   = (dz(k+1) * field%stt   + dz(kp2) * field%st   ) * indzpp
      staggered_field%sb   = (dz(k-1) * field%s     + dz(k)   * field%sb   ) * indzm
      staggered_field%en   = (dz(k)   * field%ent   + dz(k+1) * field%en   ) * indzp
      staggered_field%es   = (dz(k)   * field%est   + dz(k+1) * field%es   ) * indzp
      staggered_field%wn   = (dz(k)   * field%wnt   + dz(k+1) * field%wn   ) * indzp
      staggered_field%ws   = (dz(k)   * field%wst   + dz(k+1) * field%ws   ) * indzp
      staggered_field%wst  = (dz(k+1) * field%wstt  + dz(kp2) * field%wst  ) * indzpp
      staggered_field%wsb  = (dz(k-1) * field%ws    + dz(k)   * field%wsb  ) * indzm
      staggered_field%wnt  = (dz(k+1) * field%wntt  + dz(kp2) * field%wnt  ) * indzpp
      staggered_field%wnb  = (dz(k-1) * field%wn    + dz(k)   * field%wnb  ) * indzm
      staggered_field%ent  = (dz(k+1) * field%estt  + dz(kp2) * field%est  ) * indzpp
      staggered_field%enb  = (dz(k-1) * field%es    + dz(k)   * field%esb  ) * indzm
      staggered_field%est  = (dz(k+1) * field%entt  + dz(kp2) * field%ent  ) * indzpp
      staggered_field%esb  = (dz(k-1) * field%en    + dz(k)   * field%enb  ) * indzm
      staggered_field%ww   = (dz(k)   * field%wwt   + dz(k+1) * field%ww   ) * indzp
      staggered_field%ss   = (dz(k)   * field%sst   + dz(k+1) * field%ss   ) * indzp
      staggered_field%bb   = (dz(km2) * field%b     + dz(k-1) * field%bb   ) * indzmm
      staggered_field%ee   = (dz(k)   * field%eet   + dz(k+1) * field%ee   ) * indzp
      staggered_field%nn   = (dz(k)   * field%nnt   + dz(k+1) * field%nn   ) * indzp
      staggered_field%tt   = (dz(kp2) * field%ttt   + dz(kp3) * field%tt   ) * indzppp
      staggered_field%een  = (dz(k)   * field%eent  + dz(k+1) * field%een  ) * indzp
      staggered_field%ees  = (dz(k)   * field%eest  + dz(k+1) * field%ees  ) * indzp
      staggered_field%eet  = (dz(k+1) * field%eett  + dz(kp2) * field%eet  ) * indzpp
      staggered_field%eent = (dz(k+1) * field%eentt + dz(kp2) * field%eent ) * indzpp
      staggered_field%eest = (dz(k+1) * field%eestt + dz(kp2) * field%eest ) * indzpp
      staggered_field%wwn  = (dz(k)   * field%wwnt  + dz(k+1) * field%wwn  ) * indzp
      staggered_field%wws  = (dz(k)   * field%wwst  + dz(k+1) * field%wws  ) * indzp
      staggered_field%wwt  = (dz(k+1) * field%wwtt  + dz(kp2) * field%wwt  ) * indzpp
      staggered_field%wwnt = (dz(k+1) * field%wwntt + dz(kp2) * field%wwnt ) * indzpp
      staggered_field%wwst = (dz(k+1) * field%wwstt + dz(kp2) * field%wwst ) * indzpp
      staggered_field%enn  = (dz(k)   * field%ennt  + dz(k+1) * field%enn  ) * indzp
      staggered_field%nnt  = (dz(k+1) * field%nntt  + dz(kp2) * field%nnt  ) * indzpp
      staggered_field%nnb  = (dz(k-1) * field%nn    + dz(k)   * field%nnb  ) * indzm
      staggered_field%ennt = (dz(k+1) * field%enntt + dz(kp2) * field%ennt ) * indzpp
      staggered_field%ennb = (dz(k-1) * field%enn   + dz(k)   * field%ennb ) * indzm
      staggered_field%ess  = (dz(k)   * field%esst  + dz(k+1) * field%ess  ) * indzp
      staggered_field%sst  = (dz(k+1) * field%sstt  + dz(kp2) * field%sst  ) * indzpp
      staggered_field%ssb  = (dz(k-1) * field%ss    + dz(k)   * field%ssb  ) * indzm
      staggered_field%esst = (dz(k+1) * field%esstt + dz(kp2) * field%esst ) * indzpp
      staggered_field%essb = (dz(k-1) * field%ess   + dz(k)   * field%essb ) * indzm
      staggered_field%wtt  = (dz(kp2) * field%wttt  + dz(kp3) * field%wtt  ) * indzppp
      staggered_field%ett  = (dz(kp2) * field%ettt  + dz(kp3) * field%ett  ) * indzppp
      staggered_field%ntt  = (dz(kp2) * field%nttt  + dz(kp3) * field%ntt  ) * indzppp
      staggered_field%wntt = (dz(kp2) * field%wnttt + dz(kp3) * field%wntt ) * indzppp
      staggered_field%entt = (dz(kp2) * field%enttt + dz(kp3) * field%entt ) * indzppp
      staggered_field%wbb  = (dz(km2) * field%wb    + dz(k-1) * field%wbb  ) * indzmm
      staggered_field%ebb  = (dz(km2) * field%eb    + dz(k-1) * field%ebb  ) * indzmm
      staggered_field%nbb  = (dz(km2) * field%nb    + dz(k-1) * field%nbb  ) * indzmm
      staggered_field%wnbb = (dz(km2) * field%wnb   + dz(k-1) * field%wnbb ) * indzmm
      staggered_field%enbb = (dz(km2) * field%enb   + dz(k-1) * field%enbb ) * indzmm
!
      RETURN
      END SUBROUTINE interpolate_z
!-----------------------------------------------------------------------
!
!     I M M E R S E D   B O U N D A R I E S   
!
!     I N T E R P O L A T I O N     P R O C E D U R E S 
!
!-----------------------------------------------------------------------
      FUNCTION velint(fpt, velg, vels, ijk, cx, cy, cz) 
!
! ... Interpolate velocities on a forcing point to get no-slip
! ... conditions on a solid immersed boundary 

      USE dimensions, ONLY: nsolid, max_nsolid
      USE set_indexes, ONLY: ipjk, imjk, ippjk, immjk, ijpk, ipjpk,    &
        imjpk, ijmk, ipjmk, imjmk, ijppk, ijmmk, ijkp, ipjkp, imjkp,   &
        ijpkp, ijmkp, ijkm, ipjkm, imjkm, ijpkm, ijmkm, ijkpp, ijkmm
      USE immersed_boundaries, ONLY: forcing_point
      IMPLICIT NONE

      REAL*8, DIMENSION(max_nsolid+1) :: velint
      TYPE(forcing_point), INTENT(IN) :: fpt
      REAL*8, DIMENSION(:), INTENT(IN) :: velg
      REAL*8, DIMENSION(:,:), INTENT(IN) :: vels
      REAL*8, DIMENSION(:), INTENT(IN) :: cx, cy, cz
      INTEGER, INTENT(IN) :: ijk

      INTEGER :: i, j, k, is
      REAL*8 :: nsx, nsy, nsz
      INTEGER :: immjkpp, ippjkpp
      INTEGER :: interp

      REAL*8 :: h          !distance between the (i,j,k)-node and boundary
      REAL*8 :: zA         !distance between the boundary and the first 
                           !external node
      REAL*8 :: zB         !distance between the boundary and the second 
                           !external node
      REAL*8 :: alpha,beta !coefficients for bilinear interpolation
!
      velint = 0.D0
!
      immjkpp = imjkp
      ippjkpp = ipjkp
      
      interp = fpt%int
      i      = fpt%i
      j      = fpt%j
      k      = fpt%k
      nsx    = fpt%nsl%x
      nsy    = fpt%nsl%y
      nsz    = fpt%nsl%z

      SELECT CASE (interp)


      CASE (-3)

!===========================================
!====   interpolazione lineare con i	====
!====   valori nei nodi TW-TTWW		====
!===========================================	
   
         h=SQRT((cx(i)-nsx)**2+(cz(k)-nsz)**2)
         zA=SQRT((cx(i-1)-nsx)**2+(cz(k+1)-nsz)**2)
         IF (h <= zA) THEN
            velint(1)=-h/zA*velg(imjkp)
            DO is = 1, nsolid
              velint(is+1)=-h/zA*vels(imjkp,is)
            END DO
         ELSE
            zB=SQRT((cx(i-2)-nsx)**2+(cz(k+2)-nsz)**2)
            velint(1)=-((zB-h)*velg(imjkp)+(h-zA)*velg(immjkpp))/(zB-zA)
            DO is = 1, nsolid
              velint(is+1)=-((zB-h)*vels(imjkp,is)+(h-zA)*vels(immjkpp,is))/(zB-zA)
            END DO
         ENDIF

      CASE (-2)
   
!===========================================
!====   interpolazione lineare sin. con ====
!====   i valori nei nodi W-WW		====
!===========================================

         h=SQRT((cx(i)-nsx)**2+(cz(k)-nsz)**2)
         zA=SQRT((cx(i-1)-nsx)**2+(cz(k)-nsz)**2)
         IF (h <= zA) THEN
            velint(1)=-h/zA*velg(imjk)
            DO is= 1, nsolid
              velint(1+is)=-h/zA*vels(imjk,is)
            END DO
         ELSE
            zB=SQRT((cx(i-2)-nsx)**2+(cz(k)-nsz)**2)
            velint(1)=-((zB-h)*velg(imjk)+(h-zA)*velg(immjk))/(zB-zA)
            DO is= 1, nsolid
              velint(1+is)=-((zB-h)*vels(imjk,is)+(h-zA)*vels(immjk,is))/(zB-zA)
            END DO
         ENDIF

      CASE (-1)
   
!===========================================
!====   interpolazione bilin. con i	====
!====   valori nei nodi W-TW-T		====
!===========================================

         alpha=(cx(i-1)-nsx)/(cx(i-1)-cx(i))
         beta=(cz(k+1)-nsz)/(cz(k+1)-cz(k))
         velint(1)=-(alpha*(1-beta)*velg(ijkp)+(1-alpha)*(1-beta)*velg(imjkp)+ &
                    (1-alpha)*beta*velg(imjk))/(alpha*beta)
         DO is = 1, nsolid
           velint(1+is)=-(alpha*(1-beta)*vels(ijkp,is)+(1-alpha)*(1-beta)*vels(imjkp,is)+ &
                      (1-alpha)*beta*vels(imjk,is))/(alpha*beta)
         END DO
!
      CASE (0)
   
!===========================================
!====   interpolazione lineare con i	====
!====   valori nei nodi T-TT		====
!===========================================
   
         h=SQRT((cx(i)-nsx)**2+(cz(k)-nsz)**2)
         zA=SQRT((cx(i)-nsx)**2+(cz(k+1)-nsz)**2)
         IF (h <= zA) THEN
            velint(1)=-h/zA*velg(ijkp)
            DO is = 1, nsolid
              velint(1+is)=-h/zA*vels(ijkp,is)
            END DO
         ELSE
            zB=SQRT((cx(i)-nsx)**2+(cz(k+2)-nsz)**2)
            velint(1)=-((zB-h)*velg(ijkp)+(h-zA)*velg(ijkpp))/(zB-zA)
            DO is = 1, nsolid
              velint(1+is)=-((zB-h)*vels(ijkp,is)+(h-zA)*vels(ijkpp,is))/(zB-zA)
            END DO

! ... Quadratic interpolation
!            velint=- ((vel(ijkp)*zB-vel(ijkpp)*zA)/(zA*zB*(zA-zB))*h**2   &
!                    +(vel(ijkpp)*zA**2-vel(ijkp)*zB**2)/(zA*zB*(zA-zB))*h)
         ENDIF

      CASE (1)
   
!===========================================
!====   interpolazione bilin. con i	====
!====   valori nei nodi E-TE-T		====
!===========================================
   
         alpha=(cx(i+1)-nsx)/(cx(i+1)-cx(i))
         beta=(cz(k+1)-nsz)/(cz(k+1)-cz(k))
         velint(1)=-(alpha*(1-beta)*velg(ijkp)+(1-alpha)*(1-beta)*velg(ipjkp)+ &
                 (1-alpha)*beta*velg(ipjk))/(alpha*beta)
         DO is = 1, nsolid
           velint(1+is)=-(alpha*(1-beta)*vels(ijkp,is)+(1-alpha)*(1-beta)*vels(ipjkp,is)+ &
                   (1-alpha)*beta*vels(ipjk,is))/(alpha*beta)
         END DO
         
      CASE (2)
   
!===========================================
!====   interpolazione lineare con i	====
!====   valori nei nodi E-EE		====
!===========================================
   
         h=SQRT((cx(i)-nsx)**2+(cz(k)-nsz)**2)
         zA=SQRT((cx(i+1)-nsx)**2+(cz(k)-nsz)**2)
         IF (h <= zA) THEN
            velint(1)=-h/zA*velg(ipjk)
            DO is = 1, nsolid
              velint(1+is)=-h/zA*vels(ipjk,is)
            END DO
         ELSE
            zB=SQRT((cx(i+2)-nsx)**2+(cz(k)-nsz)**2)
            velint(1)=-((zB-h)*velg(ipjk)+(h-zA)*velg(ippjk))/(zB-zA)
            DO is = 1, nsolid
              velint(1+is)=-((zB-h)*vels(ipjk,is)+(h-zA)*vels(ippjk,is))/(zB-zA)
            END DO
         ENDIF
         
      CASE (3)
   
!===========================================
!====   interpolazione lineare con i	====
!====   valori nei nodi TE-TTEE		====
!===========================================
   
         h=SQRT((cx(i)-nsx)**2+(cz(k)-nsz)**2)
         zA=SQRT((cx(i+1)-nsx)**2+(cz(k+1)-nsz)**2)
         IF (h <= zA) THEN
            velint(1)=-h/zA*velg(ipjkp)
            DO is = 1, nsolid
              velint(1+is)=-h/zA*vels(ipjkp,is)
            END DO
         ELSE
            zB=SQRT((cx(i+2)-nsx)**2+(cz(k+2)-nsz)**2)
            velint(1)=-((zB-h)*velg(ipjkp)+(h-zA)*velg(ippjkpp))/(zB-zA)
            DO is = 1, nsolid
              velint(1+is)=-((zB-h)*vels(ipjkp,is)+(h-zA)*vels(ippjkpp,is))/(zB-zA)
            END DO
         ENDIF
   
      CASE (22)
   
!===========================================
!====   interpolazione lineare      	====
!====   esterna E           		====
!===========================================
   
         IF (nsx < cx(i)) THEN
           h=SQRT((cx(i)-nsx)**2+(cz(k)-nsz)**2)
           zA=SQRT((cx(i+1)-nsx)**2+(cz(k)-nsz)**2)
           velint(1)= + h/zA*velg(ipjk)
           DO is = 1, nsolid
             velint(1+is)= + h/zA*vels(ipjk,is)
           END DO
         ELSE
                 CALL error('bdry','control fpoint',ijk)
         END IF

      CASE (21)
   
!===========================================
!====   interpolazione lineare      	====
!====   esterna             		====
!===========================================
   
         IF (nsx > cx(i)) THEN
           h=SQRT((cx(i)-nsx)**2+(cz(k)-nsz)**2)
           zA=SQRT((cx(i-1)-nsx)**2+(cz(k)-nsz)**2)
           velint(1)= + h/zA*velg(imjk)
           DO is = 1, nsolid
             velint(1+is)= + h/zA*vels(imjk,is)
           END DO
         ELSE
                 CALL error('bdry','control fpoint',ijk)
         END IF

      CASE (20)
   
!===========================================
!====   interpolazione lineare      	====
!====   esterna TOP         		====
!===========================================
   
         h=SQRT((cx(i)-nsx)**2+(cz(k)-nsz)**2)
         zA=SQRT((cx(i)-nsx)**2+(cz(k+1)-nsz)**2)
         velint(1)= + h/zA*velg(ijkp)
         DO is = 1, nsolid
           velint(1+is)= + h/zA*vels(ijkp,is)
         END DO

! ... Quadratic interpolation
!         zB=SQRT((cx(i)-nsx)**2+(cz(k+2)-nsz)**2)
!         velint= (vel(ijkp)*zB-vel(ijkpp)*zA)/(zA*zB*(zA-zB))*h**2 &
!                +(vel(ijkpp)*zA**2-vel(ijkp)*zB**2)/(zA*zB*(zA-zB))*h

      CASE DEFAULT
   
!===========================================
!====   nessuna interpolazione  ============
!===========================================
   
         velint(1)=velg(ijk)
         DO is = 1, nsolid
           velint(1+is)=vels(ijk,is)
         END DO

      END SELECT
        
      END FUNCTION velint
!----------------------------------------------------------------------
      SUBROUTINE init_delta
! ... This routine is no longer used ...
! ... Keep it as a reminder ... !!!!!
      IMPLICIT NONE
!
      INTEGER :: delta_i_(0:30)
      INTEGER :: delta_j_(0:30)
      INTEGER :: delta_k_(0:30)
      INTEGER :: index_q_(0:30)
      INTEGER :: index_qq_(0:30)
!
! ... initialize
!
              delta_i_(:) = 0
              delta_j_(:) = 0
              delta_k_(:) = 0
!
! ... LINEAR interpolation criteria: 
! ... 0: top
! ... 1: west
! ... 2: south-west
! ... 3: south
! ... 4: south-east
! ... 5: east
! ... 6: north-east
! ... 7: north
! ... 8: north-west

              delta_i_(0) = 0
              delta_j_(0) = 0
              delta_k_(0) = 1
!
              delta_i_(1) = -1
              delta_j_(1) = 0
              delta_k_(1) = 0
!
              delta_i_(2) = -1
              delta_j_(2) = -1
              delta_k_(2) = 0
!
              delta_i_(3) = 0
              delta_j_(3) = -1
              delta_k_(3) = 0
!
              delta_i_(4) = 1
              delta_j_(4) = -1
              delta_k_(4) = 0
!
              delta_i_(5) = 1
              delta_j_(5) = 0
              delta_k_(5) = 0
!
              delta_i_(6) = 1
              delta_j_(6) = 1
              delta_k_(6) = 0
!
              delta_i_(7) = 0
              delta_j_(7) = 1
              delta_k_(7) = 0
!
              delta_i_(8) = -1
              delta_j_(8) = 1
              delta_k_(8) = 0
! ... external NE
              delta_i_(26) = 1
              delta_j_(26) = 1
              delta_k_(26) = 0
! ... external E
              delta_i_(25) = 1
              delta_j_(25) = 0
              delta_k_(25) = 0
! ... external SE
              delta_i_(24) = 1
              delta_j_(24) = -1
              delta_k_(24) = 0
! ... external N
              delta_i_(27) = 0
              delta_j_(27) = 1
              delta_k_(27) = 0
! ... external TOP
              delta_i_(20) = 0
              delta_j_(20) = 0
              delta_k_(20) = 1
! ... external S
              delta_i_(23) = 0
              delta_j_(23) = -1
              delta_k_(23) = 0
! ... external NW
              delta_i_(28) = -1
              delta_j_(28) = 1
              delta_k_(28) = 0
! ... external W
              delta_i_(21) = -1
              delta_j_(21) = 0
              delta_k_(21) = 0
! ... external SW
              delta_i_(22) = -1
              delta_j_(22) = -1
              delta_k_(22) = 0
!
      RETURN
      END SUBROUTINE init_delta
!----------------------------------------------------------------------
      FUNCTION velint3d(fpt, velg, vels, ijk, cx, cy, cz, index_q)
!
! ... Interpolate velocities on a forcing point to get no-slip
! ... conditions on a solid immersed boundary 

      USE dimensions, ONLY: nsolid, max_nsolid
      USE set_indexes
      USE immersed_boundaries, ONLY: forcing_point
      USE parallel, ONLY: mpime
      IMPLICIT NONE

      REAL*8, DIMENSION(max_nsolid+1) :: velint3d
      TYPE(forcing_point), INTENT(IN) :: fpt
      REAL*8, DIMENSION(:), INTENT(IN) :: velg
      REAL*8, DIMENSION(:,:), INTENT(IN) :: vels
      REAL*8, DIMENSION(:), INTENT(IN) :: cx, cy, cz
      INTEGER, INTENT(IN) :: ijk
      INTEGER, INTENT(OUT) :: index_q

      INTEGER :: i, j, k, is
      REAL*8 :: nsx, nsy, nsz
      INTEGER :: interp, delta_i, delta_j, delta_k
      INTEGER :: index_qq
      LOGICAL :: diagonal
      REAL*8 :: velint_1, velint_2

      REAL*8 :: h          !distance between the (i,j,k)-node and boundary
      REAL*8 :: zA         !distance between the boundary and the first 
                           !external node
      REAL*8 :: zB         !distance between the boundary and the second 
                           !external node
      REAL*8 :: hza        !ratio between the distances
      REAL*8 :: alpha,beta !coefficients for bilinear interpolation
!
      velint3d = 0.D0

      interp = fpt%int
      i      = fpt%i
      j      = fpt%j
      k      = fpt%k
      nsx    = fpt%nsl%x
      nsy    = fpt%nsl%y
      nsz    = fpt%nsl%z
!
      delta_i = fpt%delta_i
      delta_j = fpt%delta_j
      delta_k = fpt%delta_k
      index_q = fpt%index_q
      index_qq = fpt%index_qq
!
      diagonal = ( (MOD(interp,2) == 0) .AND. (interp/=0) )

      h  = SQRT( (cx(i)-nsx)**2 + (cy(j)-nsy)**2 + (cz(k)-nsz)**2 )
      
      zA = SQRT( (cx(i+delta_i)-nsx)**2 + (cy(j+delta_j)-nsy)**2 + &
                 (cz(k+delta_k)-nsz)**2 )
      hzA = h/zA

! if( index_q < 1 ) write(*,*) 'mpime,q=',mpime,index_q,i,j,k

      IF (interp >= 20) THEN
         velint3d(1) = hzA*velg(index_q)         
         DO is = 1, nsolid
           velint3d(1+is) = hzA*vels(index_q,is)         
         END DO
      ELSEIF (h <= zA .OR. diagonal) THEN
         velint3d(1) = -hzA*velg(index_q)
         DO is = 1, nsolid
           velint3d(1+is) = -hzA*vels(index_q,is)
         END DO
      ELSE 
         zB=SQRT( (cx(i+2*delta_i)-nsx)**2 + (cy(j+2*delta_j)-nsy)**2 +  &
                  (cz(k+2*delta_k)-nsz)**2 )
         !
         velint_2 = -((zB-h)*velg(index_q)+(h-zA)*velg(index_qq))/(zB-zA)
         velint_1 = -hzA*velg(index_q)
         IF (DABS(velint_1) < DABS(velint_2)) THEN
           velint3d(1) = velint_1
         ELSE
           velint3d(1) = velint_2
         END IF
         !
         DO is = 1, nsolid
           velint_2 = -((zB-h)*vels(index_q,is)+(h-zA)*vels(index_qq,is))/(zB-zA)
           velint_1 = -hzA*vels(index_q,is)
           IF (DABS(velint_1) < DABS(velint_2)) THEN
             velint3d(1+is) = velint_1
           ELSE
             velint3d(1+is) = velint_2
           END IF
         END DO
      ENDIF
        
      RETURN
      END FUNCTION velint3d
!----------------------------------------------------------------------
      FUNCTION wsink(m1, is)
      USE atmospheric_conditions, ONLY: gravz
      USE control_flags, ONLY: job_type, JOB_TYPE_2D, JOB_TYPE_3D
      USE domain_mapping, ONLY: meshinds
      USE gas_solid_velocity, ONLY: ws, wg
      USE grid, ONLY: dz
      USE particles_constants, ONLY: rl, dk
      USE gas_solid_density, ONLY: rog
      USE time_parameters, ONLY: nprint, sweep
!
      REAL*8 :: wsink
      INTEGER, INTENT(IN) :: m1, is
      INTEGER :: i, j, k, imesh, ij
      REAL*8 :: cdin
!
      CALL meshinds(m1, imesh, i, j, k)
! ... Linear extrapolation between ws(m0) (ijkp) and ws(m1) (ijk).
!      wsink = (ws(m0,is)-ws(m1,is))*(-dz(k))/dz(k+1)+ws(m1,is)
!
! ... Linear Interpolation between ws(m1) and 0.
!      wsink = dz(k-1)/(dz(k)+dz(k-1))*ws(m1,is)
!
! ... The sink velocity is set equal to the RELATIVE vertical velocity of the cell above.
!       wsink = ws(m1,is) - wg(m1)
!
! ... No interpolation/extrapolation
! ... and no resuspension
!
      wsink = MIN(0.D0, ws(m1,is))
!
! ... Constant settling velocity for single particle at high Reynolds number
!
!      cdin = 1.D0 / 0.44
!      wsink = - DSQRT(4.D0/3.D0 * cdin * 9.81 * dk(is) * rl(is) / rog(m1))
!
      RETURN
      END FUNCTION wsink
!----------------------------------------------------------------------
      REAL*8 FUNCTION extrapolate(f2,f1,d12,dx,order)
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: f1, f2, d12, dx
      INTEGER, INTENT(IN) :: order
      REAL*8 :: gradient
!
      gradient = (f2-f1)/d12 * order
      extrapolate = f2 + gradient * dx
!   
      END FUNCTION extrapolate
!----------------------------------------------------------------------
      SUBROUTINE hn(m,n)
      USE grid, ONLY: flag, filled_cell_1, filled_cell_2, fc2
      USE pressure_epsilon, ONLY: p, ep
      USE eos_gas, ONLY: ygc
      USE gas_solid_density, ONLY: rlk, rgp
      USE gas_solid_temperature, ONLY: tg, ts, sieg, sies
      USE io_files, ONLY: testunit
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: m,n
      LOGICAL :: fc2_
!
! ... m is the local cell index; 
! ... n is the index of the first neighbour used for the external forcing
!
      IF (flag(m) == filled_cell_1) THEN
        p(m) = p(n)
        ep(m) = ep(n)
        rlk(m,:) = rlk(n,:)
        rgp(m) = rgp(n)
        tg(m) = tg(n)
        ts(m,:) = ts(n,:)
        sieg(m) = sieg(n)
        sies(m,:) = sies(n,:)
        ygc(m,:) = ygc(n,:)
      ELSE IF (flag(n) == filled_cell_1) THEN
        p(n) = p(m)
        ep(n) = ep(m)
        rlk(n,:) = rlk(m,:)
        rgp(n) = rgp(m)
        tg(n) = tg(m)
        ts(n,:) = ts(m,:)
        sieg(n) = sieg(m)
        sies(n,:) = sies(m,:)
        ygc(n,:) = ygc(m,:)
      END IF
!      
      IF (flag(m) == filled_cell_2) THEN
        fc2_ = fc2(m)
        IF (.NOT.fc2_) THEN
          p(m) = 0.5D0 * p(n)
          ep(m) = 0.5D0 * ep(n)
          rlk(m,:) = 0.5D0 * rlk(n,:)
          rgp(m) = 0.5D0 * rgp(n)
          tg(m) = 0.5D0 * tg(n)
          ts(m,:) = 0.5D0 * ts(n,:)
          sieg(m) = 0.5D0 * sieg(n)
          sies(m,:) = 0.5D0 * sies(n,:)
          ygc(m,:) = 0.5D0 * ygc(n,:)
        ELSE
          p(m) = p(m) + 0.5D0 * p(n)
          ep(m) = ep(m) + 0.5D0 * ep(n)
          rlk(m,:) = rlk(m,:) + 0.5D0 * rlk(n,:)
          rgp(m) = rgp(m) + 0.5D0 * rgp(n)
          tg(m) = tg(m) + 0.5D0 * tg(n)
          ts(m,:) = ts(m,:) + 0.5D0 * ts(n,:)
          sieg(m) = sieg(m) + 0.5D0 * sieg(n)
          sies(m,:) = sies(m,:) + 0.5D0 * sies(n,:)
          ygc(m,:) = ygc(m,:) + 0.5D0 * ygc(n,:)
        END IF
        fc2(m) = .NOT.fc2_
      ELSE IF (flag(n) == filled_cell_2) THEN
        fc2_ = fc2(n)
        IF (.NOT.fc2_) THEN
          p(n) = 0.5D0 * p(m)
          ep(n) = 0.5D0 * ep(m)
          rlk(n,:) = 0.5D0 * rlk(m,:)
          rgp(n) = 0.5D0 * rgp(m)
          tg(n) = 0.5D0 * tg(m)
          ts(n,:) = 0.5D0 * ts(m,:)
          sieg(n) = 0.5D0 * sieg(m)
          sies(n,:) = 0.5D0 * sies(m,:)
          ygc(n,:) = 0.5D0 * ygc(m,:)
        ELSE
          p(n) = p(n) + 0.5D0 * p(m)
          ep(n) = ep(n) + 0.5D0 * ep(m)
          rlk(n,:) = rlk(n,:) + 0.5D0 * rlk(m,:)
          rgp(n) = rgp(n) + 0.5D0 * rgp(m)
          tg(n) = tg(n) + 0.5D0 * tg(m)
          ts(n,:) = ts(n,:) + 0.5D0 * ts(m,:)
          sieg(n) = sieg(n) + 0.5D0 * sieg(m)
          sies(n,:) = sies(n,:) + 0.5D0 * sies(m,:)
          ygc(n,:) = ygc(n,:) + 0.5D0 * ygc(m,:)
        END IF
        fc2(n) = .NOT.fc2_
      END IF 
!
      RETURN
      END SUBROUTINE hn
!----------------------------------------------------------------------
      END MODULE interpolate_fields
!----------------------------------------------------------------------
