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

      INTEGER ::  ipjk,  imjk,  ippjk,  immjk,  ijpk,  ipjpk,  imjpk,  ijmk,  &
                  ipjmk,  imjmk,  ijppk,  ijmmk,  ijkp,  ipjkp,  imjkp,  ijpkp,  &
                  ijmkp,  ijkm,  ipjkm,  imjkm,  ijpkm,  ijmkm,  ijkpp,  ijkmm

      INTEGER ::  ijke,  ijkw,  ijkee,  ijkww,  ijkn,  ijken,  ijkwn,  ijks,  ijkes,  &
                  ijkws,  ijknn,  ijkss,  ijkt,  ijket,  ijkwt,  ijknt,  ijkst,  ijkb, &
                  ijkeb,  ijkwb,  ijknb,  ijksb,  ijktt,  ijkbb


      INTEGER, PRIVATE :: job_type_flag
!
! ... change scheme => change stencil!
! ... the compass notation (N,W,S,E) is adopted for both myijk and myinds
!
      TYPE stencil
        REAL*8 :: c
        REAL*8 :: e
        REAL*8 :: w
        REAL*8 :: ee
        REAL*8 :: ww
        REAL*8 :: n
        REAL*8 :: en
        REAL*8 :: wn
        REAL*8 :: s
        REAL*8 :: es
        REAL*8 :: ws
        REAL*8 :: nn
        REAL*8 :: ss
        REAL*8 :: t
        REAL*8 :: et
        REAL*8 :: wt
        REAL*8 :: nt
        REAL*8 :: st
        REAL*8 :: b
        REAL*8 :: eb
        REAL*8 :: wb
        REAL*8 :: nb
        REAL*8 :: sb
        REAL*8 :: tt
        REAL*8 :: bb
      END TYPE stencil
!
      SAVE

!-----------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------

      SUBROUTINE subsc_setup( job_type )
        CHARACTER(LEN=80) :: job_type
        IF( job_type == '2D' ) THEN
          job_type_flag = 2
        ELSE IF( job_type == '3D' ) THEN
          job_type_flag = 3
        ELSE
          CALL error(' subsc_setup ', ' job_type not defined ', 1 )
        END IF
        RETURN
      END SUBROUTINE subsc_setup


    SUBROUTINE subscr( ijk )
!
      INTEGER, INTENT(IN) :: ijk
!
      IF( job_type_flag == 2 ) THEN

        ijm  = myijk( ip0_jm1_kp0_, ijk )
        imj  = myijk( im1_jp0_kp0_, ijk )
        ipj  = myijk( ip1_jp0_kp0_, ijk )
        ijp  = myijk( ip0_jp1_kp0_, ijk )
        ipjm = myijk( ip1_jm1_kp0_, ijk )
        ipjp = myijk( ip1_jp1_kp0_, ijk )
        imjm = myijk( im1_jm1_kp0_, ijk )
        imjp = myijk( im1_jp1_kp0_, ijk )
        ijpp = myijk( ip0_jp2_kp0_, ijk )
        ippj = myijk( ip2_jp0_kp0_, ijk )
        immj = myijk( im2_jp0_kp0_, ijk )
        ijmm = myijk( ip0_jm2_kp0_, ijk )

        ijr  = myinds(ip1_jp0_kp0_, ijk )
        ijt  = myinds(ip0_jp1_kp0_, ijk )
        ijl  = myinds(im1_jp0_kp0_, ijk )
        ijb  = myinds(ip0_jm1_kp0_, ijk )
        ijtr = myinds(ip1_jp1_kp0_, ijk )
        ijtl = myinds(im1_jp1_kp0_, ijk )
        ijbl = myinds(im1_jm1_kp0_, ijk )
        ijbr = myinds(ip1_jm1_kp0_, ijk )
        ijrr = myinds(ip2_jp0_kp0_, ijk )
        ijtt = myinds(ip0_jp2_kp0_, ijk )
        ijll = myinds(im2_jp0_kp0_, ijk )
        ijbb = myinds(ip0_jm2_kp0_, ijk )

      ELSE IF( job_type_flag == 3 ) THEN

        ipjk   = myijk( ip1_jp0_kp0_ , ijk )
        imjk   = myijk( im1_jp0_kp0_ , ijk )
        ippjk  = myijk( ip2_jp0_kp0_ , ijk )
        immjk  = myijk( im2_jp0_kp0_ , ijk )
        ijpk   = myijk( ip0_jp1_kp0_ , ijk )
        ipjpk  = myijk( ip1_jp1_kp0_ , ijk )
        imjpk  = myijk( im1_jp1_kp0_ , ijk )
        ijmk   = myijk( ip0_jm1_kp0_ , ijk )
        ipjmk  = myijk( ip1_jm1_kp0_ , ijk )
        imjmk  = myijk( im1_jm1_kp0_ , ijk )
        ijppk  = myijk( ip0_jp2_kp0_ , ijk )
        ijmmk  = myijk( ip0_jm2_kp0_ , ijk )
        ijkp   = myijk( ip0_jp0_kp1_ , ijk )
        ipjkp  = myijk( ip1_jp0_kp1_ , ijk )
        imjkp  = myijk( im1_jp0_kp1_ , ijk )
        ijpkp  = myijk( ip0_jp1_kp1_ , ijk )
        ijmkp  = myijk( ip0_jm1_kp1_ , ijk )
        ijkm   = myijk( ip0_jp0_km1_ , ijk )
        ipjkm  = myijk( ip1_jp0_km1_ , ijk )
        imjkm  = myijk( im1_jp0_km1_ , ijk )
        ijpkm  = myijk( ip0_jp1_km1_ , ijk )
        ijmkm  = myijk( ip0_jm1_km1_ , ijk )
        ijkpp  = myijk( ip0_jp0_kp2_ , ijk )
        ijkmm  = myijk( ip0_jp0_km2_ , ijk )

        ijke = myinds( ip1_jp0_kp0_ , ijk )
        ijkw = myinds( im1_jp0_kp0_ , ijk )
        ijkee = myinds( ip2_jp0_kp0_ , ijk )
        ijkww = myinds( im2_jp0_kp0_ , ijk )
        ijkn = myinds( ip0_jp1_kp0_ , ijk )
        ijken = myinds( ip1_jp1_kp0_ , ijk )
        ijkwn = myinds( im1_jp1_kp0_ , ijk )
        ijks = myinds( ip0_jm1_kp0_ , ijk )
        ijkes = myinds( ip1_jm1_kp0_ , ijk )
        ijkws = myinds( im1_jm1_kp0_ , ijk )
        ijknn = myinds( ip0_jp2_kp0_ , ijk )
        ijkss = myinds( ip0_jm2_kp0_ , ijk )
        ijkt = myinds( ip0_jp0_kp1_ , ijk )
        ijket = myinds( ip1_jp0_kp1_ , ijk )
        ijkwt = myinds( im1_jp0_kp1_ , ijk )
        ijknt = myinds( ip0_jp1_kp1_ , ijk )
        ijkst = myinds( ip0_jm1_kp1_ , ijk )
        ijkb = myinds( ip0_jp0_km1_ , ijk )
        ijkeb = myinds( ip1_jp0_km1_ , ijk )
        ijkwb = myinds( im1_jp0_km1_ , ijk )
        ijknb = myinds( ip0_jp1_km1_ , ijk )
        ijksb = myinds( ip0_jm1_km1_ , ijk )
        ijktt = myinds( ip0_jp0_kp2_ , ijk )
        ijkbb = myinds( ip0_jp0_km2_ , ijk )

      END IF

    END SUBROUTINE

!-----------------------------------------------------------------------

      FUNCTION nb( array, ijk )

      IMPLICIT NONE 
!
      TYPE(stencil) :: nb
      REAL*8, INTENT(IN) :: array(:)
      INTEGER, INTENT(IN) :: ijk
!
      IF( job_type_flag == 2 ) THEN

        nb%c  = array(ijk)
        nb%e  = array(ijr)
        nb%n  = array(ijt)
        nb%w  = array(ijl)
        nb%s  = array(ijb)
        nb%en  = array(ijtr)
        nb%wn  = array(ijtl)
        nb%es  = array(ijbr)
        nb%ws  = array(ijbl)
        nb%ee  = array(ijrr)
        nb%nn  = array(ijtt)
        nb%ww  = array(ijll)
        nb%ss  = array(ijbb)

      ELSE IF( job_type_flag == 3 ) THEN

        nb%c = array( ijk )
        nb%e = array( ijke )
        nb%w = array( ijkw )
        nb%ee = array( ijkee )
        nb%ww = array( ijkww )
        nb%n = array( ijkn )
        nb%en = array( ijken )
        nb%wn = array( ijkwn )
        nb%s = array( ijks )
        nb%es = array( ijkes )
        nb%ws = array( ijkws )
        nb%nn = array( ijknn )
        nb%ss = array( ijkss )
        nb%t = array( ijkt )
        nb%et = array( ijket )
        nb%wt = array( ijkwt )
        nb%nt = array( ijknt )
        nb%st = array( ijkst )
        nb%b = array( ijkb )
        nb%eb = array( ijkeb )
        nb%wb = array( ijkwb )
        nb%nb = array( ijknb )
        nb%sb = array( ijksb )
        nb%tt = array( ijktt )
        nb%bb = array( ijkbb )

      END IF

      RETURN
      END FUNCTION nb

!-----------------------------------------------------------------------

      FUNCTION rnb(array,ijk)
      USE dimensions
      IMPLICIT NONE 
!
      TYPE(stencil) :: rnb
      REAL*8, INTENT(IN) :: array(:)
      INTEGER, INTENT(IN) :: ijk
      INTEGER :: i,j,k,imesh
!
      imesh = myijk( ip0_jp0_kp0_, ijk)

      IF( job_type_flag == 2 ) THEN

        i  = MOD( ( imesh - 1 ), nr) + 1
        j  = ( imesh - 1 ) / nr + 1
!
        rnb%c  = array(ijk)
        rnb%e  = array(ipj)
        rnb%n  = array(ijp)
        rnb%w  = array(imj)
        rnb%s  = array(ijm)
        rnb%en  = array(ipjp)
        rnb%wn  = array(imjp)
        rnb%es  = array(ipjm)
        rnb%ws  = array(imjm)
        rnb%ee  = array(ippj)
        IF (i == (nr-1))  rnb%ee  = array(ipj)
        rnb%nn  = array(ijpp)
        IF (j == (nz-1))  rnb%nn  = array(ijp)
        rnb%ww  = array(immj)
        IF (i == 2)  rnb%ww  = array(imj)
        rnb%ss  = array(ijmm)
        IF (j == (2))  rnb%ss  = array(ijm)

      ELSE IF( job_type_flag == 3 ) THEN

         i = MOD( MOD( ijk - 1, nx*ny ), nx ) + 1
         j = MOD( ijk - 1, nx*ny ) / nx + 1
         k = ( ijk - 1 ) / ( nx*ny ) + 1

         rnb%c = array( ijk )
         rnb%e = array( ipjk )
         rnb%w = array( imjk )
         rnb%ee = array( ippjk )
         rnb%ww = array( immjk )
         rnb%n = array( ijpk )
         rnb%en = array( ipjpk )
         rnb%wn = array( imjpk )
         rnb%s = array( ijmk )
         rnb%es = array( ipjmk )
         rnb%ws = array( imjmk )
         rnb%nn = array( ijppk )
         rnb%ss = array( ijmmk )
         rnb%t = array( ijkp )
         rnb%et = array( ipjkp )
         rnb%wt = array( imjkp )
         rnb%nt = array( ijpkp )
         rnb%st = array( ijmkp )
         rnb%b = array( ijkm )
         rnb%eb = array( ipjkm )
         rnb%wb = array( imjkm )
         rnb%nb = array( ijpkm )
         rnb%sb = array( ijmkm )
         rnb%tt = array( ijkpp )
         rnb%bb = array( ijkmm )

         IF (i == (nx-1))  rnb%ee  = array(ipjk)
         IF (i == 2)       rnb%ww  = array(imjk)
         IF (j == (ny-1))  rnb%nn  = array(ijpk)
         IF (j == 2)       rnb%ss  = array(ijmk)
         IF (k == (nz-1))  rnb%tt  = array(ijkp)
         IF (k == 2)       rnb%bb  = array(ijkm)

      END IF

      RETURN
      END FUNCTION rnb

!-----------------------------------------------------------------------

      END MODULE set_indexes
