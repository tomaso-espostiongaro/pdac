!----------------------------------------------------------------------
      MODULE diffusive_fluxes
!----------------------------------------------------------------------
      USE set_indexes, ONLY: stencil
      USE indijk_module

      IMPLICIT NONE
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE hotc(hgfe, hgfn, hgft, hgfw, hgfs, hgfb, ep, tg, kapgt, ijk)
!
      USE grid, ONLY: fl_l
      USE grid, ONLY: myijk
      USE grid, ONLY: dx, dy, dz, xb
      USE dimensions
      USE set_indexes, ONLY: imjk, ijmk, ijkm
      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: hgfe, hgfn, hgft, hgfw, hgfs, hgfb
      TYPE(stencil), INTENT(IN) :: ep, tg, kapgt
      INTEGER, INTENT(IN) :: ijk
!
      REAL*8 :: dxp, dyp, dzp, dxm, dym, dzm
      REAL*8 :: indxp, indyp, indzp, indxm, indym, indzm
      REAL*8 :: kme, kmn, kmt, kmw, kms, kmb
      REAL*8 :: epme, epmn, epmt, epmw, epms, epmb
      INTEGER :: i,j,k,imesh
!
      imesh = myijk( ip0_jp0_kp0_, ijk)
      i = MOD( MOD( imesh - 1, nx*ny ), nx ) + 1
      j = MOD( imesh - 1, nx*ny ) / nx + 1
      k = ( imesh - 1 ) / ( nx*ny ) + 1
!
      IF (fl_l(imjk) /= 1) THEN
        dxm = (dx(i)+dx(i-1))
        indxm = 1.D0/dxm
        kmw = ( dx(i-1) * kapgt%c + dx(i) * kapgt%w ) * indxm
        epmw = ( dx(i-1) *ep%c + dx(i) * ep%w ) * indxm
!
        hgfw = kmw * epmw * (tg%c-tg%w) * indxm * 2.D0 * xb(i-1)
      END IF
!
      IF (fl_l(ijmk) /= 1) THEN
        dym = (dy(j)+dy(j-1))
        indym = 1.D0/dym
        kms = ( dy(j-1) * kapgt%c + dy(j) * kapgt%s ) * indym
        epms = ( dy(j-1) * ep%c + dy(j) * ep%s ) * indym
!
        hgfs = kms * epms * (tg%c-tg%s) * indym * 2.D0
      END IF
!
      IF (fl_l(ijkm) /= 1) THEN
        dzm = (dz(k)+dz(k-1))
        indzm = 1.D0/dzm
        kmb =  ( dz(k-1) * kapgt%c + dz(k) * kapgt%b ) * indzm
        epmb = ( dz(k-1) * ep%c + dz(k) * ep%b ) * indzm
!
        hgfb = kmb * epmb * (tg%c-tg%b) * indzm * 2.D0
      END IF
!
      dxp = (dx(i)+dx(i+1))
      indxp = 1.D0/dxp
      kme =  ( dx(i+1) * kapgt%c + dx(i) * kapgt%e ) * indxp
      epme = ( dx(i+1) * ep%c + dx(i) * ep%e ) * indxp
!
      hgfe = kme * epme * (tg%e-tg%c) * indxp * 2.D0 * xb(i)
!
      dyp = (dy(j)+dy(j+1))
      indyp = 1.D0/dyp
      kmn = ( dy(j+1) * kapgt%c + dy(j) * kapgt%n ) * indyp
      epmn = ( dy(j+1) * ep%c + dy(j) * ep%n ) * indyp
!
      hgfn = kmn * epmn * (tg%n-tg%c) * indyp * 2.D0
!
      dzp = (dz(k)+dz(k+1))
      indzp = 1.D0/dzp
      kmt = ( dz(k+1) * kapgt%c + dz(k) * kapgt%t ) * indzp
      epmt = ( dz(k+1) * ep%c + dz(k) * ep%t ) * indzp
!
      hgft = kmt * epmt * (tg%t-tg%c) * indzp * 2.D0
!
      RETURN
      END SUBROUTINE hotc
!----------------------------------------------------------------------
      END MODULE diffusive_fluxes
!---------------------------------------------------------------------------
