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
      INTERFACE nb
        MODULE PROCEDURE nb_rank1, nb_rank2
      END INTERFACE
      INTERFACE rnb
        MODULE PROCEDURE rnb_rank1, rnb_rank2
      END INTERFACE
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
      FUNCTION nb_rank1( array, ijk )

      IMPLICIT NONE 
!
      TYPE(stencil) :: nb_rank1
      REAL*8, INTENT(IN) :: array(:)
      INTEGER, INTENT(IN) :: ijk
!
      IF( job_type_flag == 2 ) THEN

        nb_rank1%c  = array(ijk)
        nb_rank1%e  = array(ijr)
        nb_rank1%n  = array(ijt)
        nb_rank1%w  = array(ijl)
        nb_rank1%s  = array(ijb)
        nb_rank1%en  = array(ijtr)
        nb_rank1%wn  = array(ijtl)
        nb_rank1%es  = array(ijbr)
        nb_rank1%ws  = array(ijbl)
        nb_rank1%ee  = array(ijrr)
        nb_rank1%nn  = array(ijtt)
        nb_rank1%ww  = array(ijll)
        nb_rank1%ss  = array(ijbb)

      ELSE IF( job_type_flag == 3 ) THEN

        nb_rank1%c = array( ijk )
        nb_rank1%e = array( ijke )
        nb_rank1%w = array( ijkw )
        nb_rank1%ee = array( ijkee )
        nb_rank1%ww = array( ijkww )
        nb_rank1%n = array( ijkn )
        nb_rank1%en = array( ijken )
        nb_rank1%wn = array( ijkwn )
        nb_rank1%s = array( ijks )
        nb_rank1%es = array( ijkes )
        nb_rank1%ws = array( ijkws )
        nb_rank1%nn = array( ijknn )
        nb_rank1%ss = array( ijkss )
        nb_rank1%t = array( ijkt )
        nb_rank1%et = array( ijket )
        nb_rank1%wt = array( ijkwt )
        nb_rank1%nt = array( ijknt )
        nb_rank1%st = array( ijkst )
        nb_rank1%b = array( ijkb )
        nb_rank1%eb = array( ijkeb )
        nb_rank1%wb = array( ijkwb )
        nb_rank1%nb = array( ijknb )
        nb_rank1%sb = array( ijksb )
        nb_rank1%tt = array( ijktt )
        nb_rank1%bb = array( ijkbb )

      END IF

      RETURN
      END FUNCTION nb_rank1
!-----------------------------------------------------------------------
      FUNCTION rnb_rank1(array,ijk)
      USE dimensions
      IMPLICIT NONE 
!
      TYPE(stencil) :: rnb_rank1
      REAL*8, INTENT(IN) :: array(:)
      INTEGER, INTENT(IN) :: ijk
      INTEGER :: i,j,k,imesh
!
      imesh = myijk( ip0_jp0_kp0_, ijk)

      IF( job_type_flag == 2 ) THEN

        i  = MOD( ( imesh - 1 ), nr) + 1
        j  = ( imesh - 1 ) / nr + 1
!
        rnb_rank1%c  = array(ijk)
        rnb_rank1%e  = array(ipj)
        rnb_rank1%n  = array(ijp)
        rnb_rank1%w  = array(imj)
        rnb_rank1%s  = array(ijm)
        rnb_rank1%en  = array(ipjp)
        rnb_rank1%wn  = array(imjp)
        rnb_rank1%es  = array(ipjm)
        rnb_rank1%ws  = array(imjm)
        rnb_rank1%ee  = array(ippj)
        IF (i == (nr-1))  rnb_rank1%ee  = array(ipj)
        rnb_rank1%nn  = array(ijpp)
        IF (j == (nz-1))  rnb_rank1%nn  = array(ijp)
        rnb_rank1%ww  = array(immj)
        IF (i == 2)  rnb_rank1%ww  = array(imj)
        rnb_rank1%ss  = array(ijmm)
        IF (j == (2))  rnb_rank1%ss  = array(ijm)

      ELSE IF( job_type_flag == 3 ) THEN

         i = MOD( MOD( ijk - 1, nx*ny ), nx ) + 1
         j = MOD( ijk - 1, nx*ny ) / nx + 1
         k = ( ijk - 1 ) / ( nx*ny ) + 1

         rnb_rank1%c = array( ijk )
         rnb_rank1%e = array( ipjk )
         rnb_rank1%w = array( imjk )
         rnb_rank1%ee = array( ippjk )
         rnb_rank1%ww = array( immjk )
         rnb_rank1%n = array( ijpk )
         rnb_rank1%en = array( ipjpk )
         rnb_rank1%wn = array( imjpk )
         rnb_rank1%s = array( ijmk )
         rnb_rank1%es = array( ipjmk )
         rnb_rank1%ws = array( imjmk )
         rnb_rank1%nn = array( ijppk )
         rnb_rank1%ss = array( ijmmk )
         rnb_rank1%t = array( ijkp )
         rnb_rank1%et = array( ipjkp )
         rnb_rank1%wt = array( imjkp )
         rnb_rank1%nt = array( ijpkp )
         rnb_rank1%st = array( ijmkp )
         rnb_rank1%b = array( ijkm )
         rnb_rank1%eb = array( ipjkm )
         rnb_rank1%wb = array( imjkm )
         rnb_rank1%nb = array( ijpkm )
         rnb_rank1%sb = array( ijmkm )
         rnb_rank1%tt = array( ijkpp )
         rnb_rank1%bb = array( ijkmm )

         IF (i == (nx-1))  rnb_rank1%ee  = array(ipjk)
         IF (i == 2)       rnb_rank1%ww  = array(imjk)
         IF (j == (ny-1))  rnb_rank1%nn  = array(ijpk)
         IF (j == 2)       rnb_rank1%ss  = array(ijmk)
         IF (k == (nz-1))  rnb_rank1%tt  = array(ijkp)
         IF (k == 2)       rnb_rank1%bb  = array(ijkm)

      END IF

      RETURN
      END FUNCTION rnb_rank1
!-----------------------------------------------------------------------
      FUNCTION nb_rank2( array, is, ijk )

      IMPLICIT NONE 
!
      TYPE(stencil) :: nb_rank2
      REAL*8, INTENT(IN) :: array(:,:)
      INTEGER, INTENT(IN) :: is, ijk
!
      IF( job_type_flag == 2 ) THEN

        nb_rank2%c  = array(is, ijk)
        nb_rank2%e  = array(is, ijr)
        nb_rank2%n  = array(is, ijt)
        nb_rank2%w  = array(is, ijl)
        nb_rank2%s  = array(is, ijb)
        nb_rank2%en  = array(is, ijtr)
        nb_rank2%wn  = array(is, ijtl)
        nb_rank2%es  = array(is, ijbr)
        nb_rank2%ws  = array(is, ijbl)
        nb_rank2%ee  = array(is, ijrr)
        nb_rank2%nn  = array(is, ijtt)
        nb_rank2%ww  = array(is, ijll)
        nb_rank2%ss  = array(is, ijbb)

      ELSE IF( job_type_flag == 3 ) THEN

        nb_rank2%c = array( is, ijk )
        nb_rank2%e = array( is, ijke )
        nb_rank2%w = array( is, ijkw )
        nb_rank2%ee = array( is, ijkee )
        nb_rank2%ww = array( is, ijkww )
        nb_rank2%n = array( is, ijkn )
        nb_rank2%en = array( is, ijken )
        nb_rank2%wn = array( is, ijkwn )
        nb_rank2%s = array( is, ijks )
        nb_rank2%es = array( is, ijkes )
        nb_rank2%ws = array( is, ijkws )
        nb_rank2%nn = array( is, ijknn )
        nb_rank2%ss = array( is, ijkss )
        nb_rank2%t = array( is, ijkt )
        nb_rank2%et = array( is, ijket )
        nb_rank2%wt = array( is, ijkwt )
        nb_rank2%nt = array( is, ijknt )
        nb_rank2%st = array( is, ijkst )
        nb_rank2%b = array( is, ijkb )
        nb_rank2%eb = array( is, ijkeb )
        nb_rank2%wb = array( is, ijkwb )
        nb_rank2%nb = array( is, ijknb )
        nb_rank2%sb = array( is, ijksb )
        nb_rank2%tt = array( is, ijktt )
        nb_rank2%bb = array( is, ijkbb )

      END IF

      RETURN
      END FUNCTION nb_rank2
!-----------------------------------------------------------------------
      FUNCTION rnb_rank2(array, is, ijk)
      USE dimensions
      IMPLICIT NONE 
!
      TYPE(stencil) :: rnb_rank2
      REAL*8, INTENT(IN) :: array(:, :)
      INTEGER, INTENT(IN) :: is, ijk
      INTEGER :: i,j,k,imesh
!
      imesh = myijk( ip0_jp0_kp0_, ijk)

      IF( job_type_flag == 2 ) THEN

        i  = MOD( ( imesh - 1 ), nr) + 1
        j  = ( imesh - 1 ) / nr + 1
!
        rnb_rank2%c  = array(is, ijk)
        rnb_rank2%e  = array(is, ipj)
        rnb_rank2%n  = array(is, ijp)
        rnb_rank2%w  = array(is, imj)
        rnb_rank2%s  = array(is, ijm)
        rnb_rank2%en  = array(is, ipjp)
        rnb_rank2%wn  = array(is, imjp)
        rnb_rank2%es  = array(is, ipjm)
        rnb_rank2%ws  = array(is, imjm)
        rnb_rank2%ee  = array(is, ippj)
        IF (i == (nr-1))  rnb_rank2%ee  = array(is, ipj)
        rnb_rank2%nn  = array(is, ijpp)
        IF (j == (nz-1))  rnb_rank2%nn  = array(is, ijp)
        rnb_rank2%ww  = array(is, immj)
        IF (i == 2)  rnb_rank2%ww  = array(is, imj)
        rnb_rank2%ss  = array(is, ijmm)
        IF (j == (2))  rnb_rank2%ss  = array(is, ijm)

      ELSE IF( job_type_flag == 3 ) THEN

         i = MOD( MOD( ijk - 1, nx*ny ), nx ) + 1
         j = MOD( ijk - 1, nx*ny ) / nx + 1
         k = ( ijk - 1 ) / ( nx*ny ) + 1

         rnb_rank2%c = array( is, ijk )
         rnb_rank2%e = array( is, ipjk )
         rnb_rank2%w = array( is, imjk )
         rnb_rank2%ee = array( is, ippjk )
         rnb_rank2%ww = array( is, immjk )
         rnb_rank2%n = array( is, ijpk )
         rnb_rank2%en = array( is, ipjpk )
         rnb_rank2%wn = array( is, imjpk )
         rnb_rank2%s = array( is, ijmk )
         rnb_rank2%es = array( is, ipjmk )
         rnb_rank2%ws = array( is, imjmk )
         rnb_rank2%nn = array( is, ijppk )
         rnb_rank2%ss = array( is, ijmmk )
         rnb_rank2%t = array( is, ijkp )
         rnb_rank2%et = array( is, ipjkp )
         rnb_rank2%wt = array( is, imjkp )
         rnb_rank2%nt = array( is, ijpkp )
         rnb_rank2%st = array( is, ijmkp )
         rnb_rank2%b = array( is, ijkm )
         rnb_rank2%eb = array( is, ipjkm )
         rnb_rank2%wb = array( is, imjkm )
         rnb_rank2%nb = array( is, ijpkm )
         rnb_rank2%sb = array( is, ijmkm )
         rnb_rank2%tt = array( is, ijkpp )
         rnb_rank2%bb = array( is, ijkmm )

         IF (i == (nx-1))  rnb_rank2%ee  = array(is, ipjk)
         IF (i == 2)       rnb_rank2%ww  = array(is, imjk)
         IF (j == (ny-1))  rnb_rank2%nn  = array(is, ijpk)
         IF (j == 2)       rnb_rank2%ss  = array(is, ijmk)
         IF (k == (nz-1))  rnb_rank2%tt  = array(is, ijkp)
         IF (k == 2)       rnb_rank2%bb  = array(is, ijkm)

      END IF

      RETURN
      END FUNCTION rnb_rank2
!-----------------------------------------------------------------------
      END MODULE set_indexes
!-----------------------------------------------------------------------
