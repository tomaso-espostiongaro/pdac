!-----------------------------------------------------------------------
      MODULE set_indexes
!-----------------------------------------------------------------------

      USE grid, ONLY: myinds, myijk
      USE grid, ONLY: r_, t_, l_, b_, tr_, tl_, br_, bl_
      USE grid, ONLY: rr_, tt_, ll_, bb_
      USE indijk_module

      IMPLICIT NONE
!
      INTEGER :: ijr, ijt, ijl, ijb
      INTEGER :: ijtr, ijtl, ijbr, ijbl
      INTEGER :: ijtt, ijrr, ijbb, ijll
      INTEGER :: imj, ijm, ipj, ijp
      INTEGER :: ipjm, ipjp, imjm, imjp 
      INTEGER :: ijpp, ippj, ijmm, immj 
!
! ... change scheme => change stencil!
! ... the compass notation (N,W,S,E) is adopted for both myijk and myinds
!
      TYPE stencil
        REAL*8 :: c
        REAL*8 :: e
        REAL*8 :: n
        REAL*8 :: w
        REAL*8 :: s
        REAL*8 :: ne
        REAL*8 :: nw
        REAL*8 :: se
        REAL*8 :: sw
        REAL*8 :: ee
        REAL*8 :: nn
        REAL*8 :: ww
        REAL*8 :: ss
      END TYPE stencil
!
      SAVE
!-----------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------
      SUBROUTINE subscr(ij)
!
      INTEGER, INTENT(IN) :: ij
!
      ijm  = myijk( ip0_jm1_kp0_, ij)
      imj  = myijk( im1_jp0_kp0_, ij)
      ipj  = myijk( ip1_jp0_kp0_, ij)
      ijp  = myijk( ip0_jp1_kp0_, ij)
      ipjm = myijk( ip1_jm1_kp0_, ij)
      ipjp = myijk( ip1_jp1_kp0_, ij)
      imjm = myijk( im1_jm1_kp0_, ij)
      imjp = myijk( im1_jp1_kp0_, ij)
      ijpp = myijk( ip0_jp2_kp0_, ij)
      ippj = myijk( ip2_jp0_kp0_, ij)
      immj = myijk( im2_jp0_kp0_, ij)
      ijmm = myijk( ip0_jm2_kp0_, ij)
!
      ijr  = myinds(r_, ij)
      ijt  = myinds(t_, ij)
      ijl  = myinds(l_, ij)
      ijb  = myinds(b_, ij)
      ijtr = myinds(tr_, ij)
      ijtl = myinds(tl_, ij)
      ijbl = myinds(bl_, ij)
      ijbr = myinds(br_, ij)
      ijrr = myinds(rr_, ij)
      ijtt = myinds(tt_, ij) 
      ijll = myinds(ll_, ij)
      ijbb = myinds(bb_, ij) 
!
      END SUBROUTINE
!-----------------------------------------------------------------------
      FUNCTION nb(array,ij)
      USE dimensions
      IMPLICIT NONE 
!
      TYPE(stencil) :: nb
      REAL*8, INTENT(IN) :: array(:)
      INTEGER, INTENT(IN) :: ij
!
      nb%c  = array(ij)
      nb%e  = array(ijr)
      nb%n  = array(ijt)
      nb%w  = array(ijl)
      nb%s  = array(ijb)
      nb%ne  = array(ijtr)
      nb%nw  = array(ijtl)
      nb%se  = array(ijbr)
      nb%sw  = array(ijbl)
      nb%ee  = array(ijrr)
      nb%nn  = array(ijtt)
      nb%ww  = array(ijll)
      nb%ss  = array(ijbb)

      RETURN
      END FUNCTION nb
!-----------------------------------------------------------------------
      FUNCTION rnb(array,ij)
      USE dimensions
      IMPLICIT NONE 
!
      TYPE(stencil) :: rnb
      REAL*8, INTENT(IN) :: array(:)
      INTEGER, INTENT(IN) :: ij
      INTEGER :: i,j,imesh
!
      imesh = myijk( ip0_jp0_kp0_, ij)
      i  = MOD( ( imesh - 1 ), nr) + 1
      j  = ( imesh - 1 ) / nr + 1
!
      rnb%c  = array(ij)
      rnb%e  = array(ipj)
      rnb%n  = array(ijp)
      rnb%w  = array(imj)
      rnb%s  = array(ijm)
      rnb%ne  = array(ipjp)
      rnb%nw  = array(imjp)
      rnb%se  = array(ipjm)
      rnb%sw  = array(imjm)
      rnb%ee  = array(ippj)
      IF (i == (nr-1))  rnb%ee  = array(ipj)
      rnb%nn  = array(ijpp)
      IF (j == (nz-1))  rnb%nn  = array(ijp)
      rnb%ww  = array(immj)
      IF (i == 2)  rnb%ww  = array(imj)
      rnb%ss  = array(ijmm)
      IF (j == (2))  rnb%ss  = array(ijm)

      RETURN
      END FUNCTION rnb
!-----------------------------------------------------------------------
      FUNCTION const(c)
      IMPLICIT NONE 
!
      REAL*8, INTENT(IN) :: c
      TYPE(stencil) :: const
!
      const%c  = c
      const%e  = c
      const%n  = c
      const%w  = c
      const%s  = c
      const%ne  = c
      const%nw  = c
      const%se  = c
      const%sw  = c
      const%ee  = c
      const%nn  = c
      const%ww  = c
      const%ss  = c

      RETURN
      END FUNCTION const
!-----------------------------------------------------------------------
      END MODULE set_indexes
