!----------------------------------------------------------------------
      MODULE interpolate_fields
!----------------------------------------------------------------------
!
      USE dimensions
      USE grid, ONLY: dx, dy, dz
      USE set_indexes, ONLY: stencil

      IMPLICIT NONE
!
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE interpolate_x(field, staggered_field, i)
!
! ... interpolate density on the mesh staggered along x
!
      IMPLICIT NONE

      TYPE(stencil), INTENT(IN) :: field
      TYPE(stencil), INTENT(OUT) :: staggered_field
      INTEGER, INTENT(IN) :: i

      INTEGER :: ip2

      REAL*8 :: dxm, dxp, dxpp, indxpp, indxp, indxm
!
      ip2 = MIN( nx, i+2 )

      dxm  = dx(i)   + dx(i-1)
      dxp  = dx(i)   + dx(i+1)
      dxpp = dx(i+1) + dx(ip2)

      indxm  = 1.D0 / dxm
      indxp  = 1.D0 / dxp
      indxpp = 1.D0 / dxpp
!       
! ... Compute linearly interpolated values of density on the staggered grid
!
      staggered_field%c = (dx(i+1) * field%c + dx(i)   * field%e ) * indxp
      staggered_field%e = (dx(ip2) * field%e + dx(i+1) * field%ee) * indxpp
      staggered_field%n = (dx(i+1) * field%n + dx(i)   * field%en) * indxp
      staggered_field%t = (dx(i+1) * field%t + dx(i)   * field%et) * indxp
      staggered_field%w = (dx(i)   * field%w + dx(i-1) * field%c ) * indxm
      staggered_field%s = (dx(i+1) * field%s + dx(i)   * field%es) * indxp
      staggered_field%b = (dx(i+1) * field%b + dx(i)   * field%eb) * indxp
!
! ... An arbitrary choice !
!
      staggered_field%ee = staggered_field%e
      staggered_field%nn = staggered_field%n
      staggered_field%tt = staggered_field%t

      RETURN
      END SUBROUTINE interpolate_x
!----------------------------------------------------------------------
      SUBROUTINE interpolate_y(field, staggered_field, j)
!
! ... interpolate density on the mesh staggered along x
!
      IMPLICIT NONE
!
      TYPE(stencil), INTENT(IN) :: field
      TYPE(stencil), INTENT(OUT) :: staggered_field
      INTEGER, INTENT(IN) :: j

      REAL*8 :: dym, dyp, dypp, indypp, indyp, indym

      INTEGER :: jp2
!
      jp2 = MIN( ny, j+2 )

      dym = dy(j) + dy(j-1)
      dyp = dy(j) + dy(j+1)
      dypp = dy(j+1) + dy(jp2)

      indym=1.D0 / dym
      indyp=1.D0 / dyp
      indypp=1.D0 / dypp
!       
! ... Compute linearly interpolated values of density on the staggered grid
!
      staggered_field%c = (dy(j+1) * field%c + dy(j)   * field%n ) * indyp
      staggered_field%e = (dy(j+1) * field%e + dy(j)   * field%en) * indyp
      staggered_field%n = (dy(jp2) * field%n + dy(j+1) * field%nn) * indypp
      staggered_field%t = (dy(j+1) * field%t + dy(j)   * field%nt) * indyp
      staggered_field%w = (dy(j+1) * field%w + dy(j)   * field%wn) * indyp
      staggered_field%s = (dy(j)   * field%s + dy(j-1) * field%c ) * indym
      staggered_field%b = (dy(j+1) * field%b + dy(j)   * field%nb) * indyp
!
! ... an arbitrary choice !
!
      staggered_field%ee = staggered_field%e
      staggered_field%nn = staggered_field%n
      staggered_field%tt = staggered_field%t
!
      RETURN
      END SUBROUTINE interpolate_y
!----------------------------------------------------------------------
      SUBROUTINE interpolate_z(field, staggered_field, k)
!
! ... interpolate density on the mesh staggered along x
!
      IMPLICIT NONE
!
      TYPE(stencil), INTENT(IN) :: field
      TYPE(stencil), INTENT(OUT) :: staggered_field
      INTEGER, INTENT(IN) :: k

      REAL*8 :: dzp, indzp, dzm, indzm, dzpp, indzpp

      INTEGER :: kp2
!
      kp2 = MIN( nz, k+2 )

      dzm=dz(k)+dz(k-1)
      dzp=dz(k)+dz(k+1)
      dzpp=dz(k+1)+dz(kp2)

      indzm=1.D0/dzm
      indzp=1.D0/dzp
      indzpp=1.D0/dzpp
!       
! ... values of density interpolated linearly on the staggered grid
!
      staggered_field%c = (dz(k+1) * field%c + dz(k)   * field%t ) * indzp
      staggered_field%e = (dz(k+1) * field%e + dz(k)   * field%et) * indzp
      staggered_field%n = (dz(k+1) * field%n + dz(k)   * field%nt) * indzp
      staggered_field%t = (dz(kp2) * field%t + dz(k+1) * field%tt) * indzpp
      staggered_field%w = (dz(k+1) * field%w + dz(k)   * field%wt) * indzp
      staggered_field%s = (dz(k+1) * field%s + dz(k)   * field%st) * indzp
      staggered_field%b = (dz(k)   * field%b + dz(k-1) * field%c ) * indzm
!
! ... an arbitrary choice !
!
      staggered_field%ee = staggered_field%e
      staggered_field%nn = staggered_field%n
      staggered_field%tt = staggered_field%t

      RETURN
      END SUBROUTINE interpolate_z
!----------------------------------------------------------------------
      END MODULE interpolate_fields
!----------------------------------------------------------------------
