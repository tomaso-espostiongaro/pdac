!----------------------------------------------------------------------
      MODULE heat_diffusion
!----------------------------------------------------------------------
      USE set_indexes, ONLY: stencil
      IMPLICIT NONE
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE hotcg(hfgr, hfgt, hfgl, hfgb, ep, tg, kapgt, ij)
!
      USE grid, ONLY: fl_l
      USE grid, ONLY: myij
      USE grid, ONLY: dz, dr, rb
      USE dimensions
      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: hfgr, hfgt, hfgl, hfgb
      TYPE(stencil), INTENT(IN) :: ep, tg, kapgt
      INTEGER, INTENT(IN) :: ij
!
      REAL*8 :: drp3, dzp3, dzm3, drm3
      REAL*8 :: indrp3, indzp3, indzm3, indrm3
      REAL*8 :: kml, kmb, kmr, kmt
      REAL*8 :: epml, epmb, epmr, epmt
      INTEGER :: i,j,ij_g
      INTEGER :: imj, ijm, ipj, ijp
!
      ij_g = myij( 0, 0, ij)
      j  = ( ij_g - 1 ) / nr + 1
      i  = MOD( ( ij_g - 1 ), nr) + 1
!
       imj = myij(-1, 0, ij)
       IF (fl_l(imj) .NE. 1) THEN
         drm3 = (dr(i)+dr(i-1))**3
         indrm3 = 1.D0/drm3
         kml = dr(i-1) * kapgt%c + dr(i) * kapgt%w
         epml = dr(i-1) *ep%c +dr(i) * ep%w
!
         hfgl = rb(i-1) * kml * epml * (tg%c-tg%w) * indrm3 * 2.D0
       END IF
!
       ijm = myij( 0,-1, ij)
       IF (fl_l(ijm) .NE. 1) THEN
         dzm3 = (dz(j)+dz(j-1))**3
         indzm3 = 1.D0/dzm3
         kmb = dz(j-1) * kapgt%c + dz(j) * kapgt%s
         epmb = dz(j-1) * ep%c + dz(j) * ep%s
!
         hfgb = kmb * epmb * (tg%c-tg%s) * indzm3 * 2.D0
       END IF
!
       drp3 = (dr(i)+dr(i+1))**3
       indrp3 = 1.D0/drp3
       kmr = dr(i+1) * kapgt%c + dr(i) * kapgt%e
       epmr = dr(i+1) * ep%c + dr(i) * ep%e
!
       hfgr = rb(i) * kmr * epmr * (tg%e-tg%c) * indrp3 * 2.D0
       ipj = myij(+1, 0, ij)
       IF (fl_l(ipj) .EQ. 3) hfgr = 0.D0
!
       dzp3 = (dz(j)+dz(j+1))**3
       indzp3 = 1.D0/dzp3
       kmt = dz(j+1) * kapgt%c + dz(j) * kapgt%n
       epmt = dz(j+1) * ep%c + dz(j) * ep%n
!
       hfgt = kmt * epmt * (tg%n-tg%c) * indzp3 * 2.D0
       ijp = myij( 0, +1, ij)
       IF (fl_l(ijp) .EQ. 3) hfgt = 0.D0
!
      RETURN
      END SUBROUTINE hotcg
!----------------------------------------------------------------------
      SUBROUTINE hotck(hflr, hflt, hfll, hflb, rlk, tk, ij, k)
!
      USE dimensions
      USE grid, ONLY: fl_l
      USE grid, ONLY: myij
      USE grid, ONLY: dz, dr, rb
      USE particles_constants, ONLY: rl, inrl, kap
      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: hflr, hflt, hfll, hflb
      TYPE(stencil), INTENT(IN) :: rlk, tk
      INTEGER, INTENT(IN) :: ij, k
!
      REAL*8 :: drp2, dzp2, dzm2, drm2
      REAL*8 :: indrp2, indzp2, indzm2, indrm2
      REAL*8 :: rlkml, rlkmb, rlkmr, rlkmt
      INTEGER :: i, j, ij_g
      INTEGER :: imj, ijm, ipj, ijp
!
      ij_g = myij( 0, 0, ij)
      j  = ( ij_g - 1 ) / nr + 1
      i  = MOD( ( ij_g - 1 ), nr) + 1
!
      imj = myij(-1, 0, ij)
      IF (fl_l(imj) .NE. 1) THEN
        drm2 = (dr(i)+dr(i-1))**2
        indrm2 = 1.D0/drm2
        rlkml = dr(i-1) * rlk%c + dr(i) * rlk%w
!
        hfll = rb(i-1) * kap(k) * rlkml * inrl(k) * (tk%c-tk%w) * indrm2 * 2.D0
      END IF
!
      ijm = myij( 0,-1, ij)
      IF (fl_l(ijm) .NE. 1) THEN
        dzm2 = (dz(j)+dz(j-1))**2
        indzm2 = 1.D0/dzm2 
        rlkmb = dz(j-1) * rlk%c + dz(j) * rlk%s
!
        hflb = kap(k) * rlkmb * inrl(k) * (tk%c-tk%s) * indzm2 * 2.D0
      END IF
!
      drp2 = (dr(i)+dr(i+1))**2
      indrp2 = 1.D0/drp2
      rlkmr = dr(i+1) * rlk%c + dr(i) * rlk%e
!
      hflr = rb(i) * kap(k) * rlkmr * inrl(k) * (tk%e-tk%c) * indrp2 * 2.D0
      ipj = myij(+1, 0, ij)
      IF (fl_l(ipj) .EQ. 3) hflr = 0.D0
!
      dzp2 = (dz(j)+dz(j+1))**2
      indzp2 = 1.D0/dzp2
      rlkmt = dz(j+1) * rlk%c + dz(j) * rlk%n
!
      hflt = kap(k) * rlkmt * inrl(k) * (tk%n-tk%c) * indzp2 * 2.D0
      ijp = myij( 0,+1, ij)
      IF (fl_l(ijp) .EQ. 3) hflt = 0.D0
!
      RETURN
!---------------------------------------------------------------------------
      END SUBROUTINE hotck
!---------------------------------------------------------------------------
      END MODULE heat_diffusion
!---------------------------------------------------------------------------
