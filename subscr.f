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

      INTEGER :: ipjk, imjk, ippjk, immjk, ijpk, ipjpk, imjpk, ijmk,  &
                 ipjmk, imjmk, ijppk, ijmmk, ijkp, ipjkp, imjkp, ijpkp, &
                 ijmkp, ijkm, ipjkm, imjkm, ijpkm, ijmkm, ijkpp, ijkmm

      INTEGER ::  ijke, ijkw, ijkee, ijkww, ijkn, ijken, ijkwn, ijks, ijkes, &
                  ijkws, ijknn, ijkss, ijkt, ijket, ijkwt, ijknt, ijkst, ijkb, &
                  ijkeb, ijkwb, ijknb, ijksb, ijktt, ijkbb


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

      INTERFACE OPERATOR(+)
        MODULE PROCEDURE sumstencil
      END INTERFACE
      INTERFACE OPERATOR(-)
        MODULE PROCEDURE difstencil
      END INTERFACE
      INTERFACE OPERATOR(*)
        MODULE PROCEDURE prodstencil, dotstencil
      END INTERFACE
      INTERFACE OPERATOR(/)
        MODULE PROCEDURE fracstencil
      END INTERFACE
!
      SAVE
!-----------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------

      SUBROUTINE subsc_setup( job_type )
      IMPLICIT NONE
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
!-----------------------------------------------------------------------
      SUBROUTINE subscr( ijk )
      IMPLICIT NONE
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
      FUNCTION cte(c)
      IMPLICIT NONE
        REAL*8, INTENT(IN) :: c
        TYPE(stencil) :: cte

        cte%c  = c
        cte%e  = c
        cte%w  = c
        cte%ee = c
        cte%ww = c
        cte%n  = c
        cte%en = c
        cte%wn = c
        cte%s  = c
        cte%es = c
        cte%ws = c
        cte%nn = c
        cte%ss = c
        cte%t  = c
        cte%et = c
        cte%wt = c
        cte%nt = c
        cte%st = c
        cte%b  = c
        cte%eb = c
        cte%wb = c
        cte%nb = c
        cte%sb = c
        cte%tt = c
        cte%bb = c
      
      END FUNCTION cte
!-----------------------------------------------------------------------
      FUNCTION sumstencil(st1, st2)
      IMPLICIT NONE
        TYPE(stencil), INTENT(in) :: st1, st2
        TYPE(stencil) :: sumstencil

        sumstencil%c  = st1%c  + st2%c
        sumstencil%e  = st1%e  + st2%e
        sumstencil%w  = st1%w  + st2%w
        sumstencil%ee = st1%ee + st2%ee
        sumstencil%ww = st1%ww + st2%ww
        sumstencil%n  = st1%n  + st2%n
        sumstencil%en = st1%en + st2%en
        sumstencil%wn = st1%wn + st2%wn
        sumstencil%s  = st1%s  + st2%s
        sumstencil%es = st1%es + st2%es
        sumstencil%ws = st1%ws + st2%ws
        sumstencil%nn = st1%nn + st2%nn
        sumstencil%ss = st1%ss + st2%ss
        sumstencil%t  = st1%t  + st2%t
        sumstencil%et = st1%et + st2%et
        sumstencil%wt = st1%wt + st2%wt
        sumstencil%nt = st1%nt + st2%nt
        sumstencil%st = st1%st + st2%st
        sumstencil%b  = st1%b  + st2%b
        sumstencil%eb = st1%eb + st2%eb
        sumstencil%wb = st1%wb + st2%wb
        sumstencil%nb = st1%nb + st2%nb
        sumstencil%sb = st1%sb + st2%sb
        sumstencil%tt = st1%tt + st2%tt
        sumstencil%bb = st1%bb + st2%bb
      
      END FUNCTION sumstencil
!-----------------------------------------------------------------------
      FUNCTION difstencil(st1, st2)
      IMPLICIT NONE
        TYPE(stencil), INTENT(in) :: st1, st2
        TYPE(stencil) :: difstencil

        difstencil%c  = st1%c  - st2%c
        difstencil%e  = st1%e  - st2%e
        difstencil%w  = st1%w  - st2%w
        difstencil%ee = st1%ee - st2%ee
        difstencil%ww = st1%ww - st2%ww
        difstencil%n  = st1%n  - st2%n
        difstencil%en = st1%en - st2%en
        difstencil%wn = st1%wn - st2%wn
        difstencil%s  = st1%s  - st2%s
        difstencil%es = st1%es - st2%es
        difstencil%ws = st1%ws - st2%ws
        difstencil%nn = st1%nn - st2%nn
        difstencil%ss = st1%ss - st2%ss
        difstencil%t  = st1%t  - st2%t
        difstencil%et = st1%et - st2%et
        difstencil%wt = st1%wt - st2%wt
        difstencil%nt = st1%nt - st2%nt
        difstencil%st = st1%st - st2%st
        difstencil%b  = st1%b  - st2%b
        difstencil%eb = st1%eb - st2%eb
        difstencil%wb = st1%wb - st2%wb
        difstencil%nb = st1%nb - st2%nb
        difstencil%sb = st1%sb - st2%sb
        difstencil%tt = st1%tt - st2%tt
        difstencil%bb = st1%bb - st2%bb
      
      END FUNCTION difstencil
!-----------------------------------------------------------------------
      FUNCTION prodstencil(st1, st2)
      IMPLICIT NONE
        TYPE(stencil), INTENT(in) :: st1, st2
        TYPE(stencil) :: prodstencil

        prodstencil%c  = st1%c  * st2%c
        prodstencil%e  = st1%e  * st2%e
        prodstencil%w  = st1%w  * st2%w
        prodstencil%ee = st1%ee * st2%ee
        prodstencil%ww = st1%ww * st2%ww
        prodstencil%n  = st1%n  * st2%n
        prodstencil%en = st1%en * st2%en
        prodstencil%wn = st1%wn * st2%wn
        prodstencil%s  = st1%s  * st2%s
        prodstencil%es = st1%es * st2%es
        prodstencil%ws = st1%ws * st2%ws
        prodstencil%nn = st1%nn * st2%nn
        prodstencil%ss = st1%ss * st2%ss
        prodstencil%t  = st1%t  * st2%t
        prodstencil%et = st1%et * st2%et
        prodstencil%wt = st1%wt * st2%wt
        prodstencil%nt = st1%nt * st2%nt
        prodstencil%st = st1%st * st2%st
        prodstencil%b  = st1%b  * st2%b
        prodstencil%eb = st1%eb * st2%eb
        prodstencil%wb = st1%wb * st2%wb
        prodstencil%nb = st1%nb * st2%nb
        prodstencil%sb = st1%sb * st2%sb
        prodstencil%tt = st1%tt * st2%tt
        prodstencil%bb = st1%bb * st2%bb
      
      END FUNCTION prodstencil
!-----------------------------------------------------------------------
      FUNCTION dotstencil(a, st2)
      IMPLICIT NONE
        REAL*8, INTENT(IN) :: a
        TYPE(stencil), INTENT(IN) :: st2
        TYPE(stencil) :: dotstencil

        dotstencil%c  = a * st2%c
        dotstencil%e  = a * st2%e
        dotstencil%w  = a * st2%w
        dotstencil%ee = a * st2%ee
        dotstencil%ww = a * st2%ww
        dotstencil%n  = a * st2%n
        dotstencil%en = a * st2%en
        dotstencil%wn = a * st2%wn
        dotstencil%s  = a * st2%s
        dotstencil%es = a * st2%es
        dotstencil%ws = a * st2%ws
        dotstencil%nn = a * st2%nn
        dotstencil%ss = a * st2%ss
        dotstencil%t  = a * st2%t
        dotstencil%et = a * st2%et
        dotstencil%wt = a * st2%wt
        dotstencil%nt = a * st2%nt
        dotstencil%st = a * st2%st
        dotstencil%b  = a * st2%b
        dotstencil%eb = a * st2%eb
        dotstencil%wb = a * st2%wb
        dotstencil%nb = a * st2%nb
        dotstencil%sb = a * st2%sb
        dotstencil%tt = a * st2%tt
        dotstencil%bb = a * st2%bb
      
      END FUNCTION dotstencil
!-----------------------------------------------------------------------
      FUNCTION fracstencil(st1, st2)
      IMPLICIT NONE
        TYPE(stencil), INTENT(in) :: st1, st2
        TYPE(stencil) :: fracstencil

        fracstencil%c  = st1%c  / st2%c
        fracstencil%e  = st1%e  / st2%e
        fracstencil%w  = st1%w  / st2%w
        fracstencil%ee = st1%ee / st2%ee
        fracstencil%ww = st1%ww / st2%ww
        fracstencil%n  = st1%n  / st2%n
        fracstencil%en = st1%en / st2%en
        fracstencil%wn = st1%wn / st2%wn
        fracstencil%s  = st1%s  / st2%s
        fracstencil%es = st1%es / st2%es
        fracstencil%ws = st1%ws / st2%ws
        fracstencil%nn = st1%nn / st2%nn
        fracstencil%ss = st1%ss / st2%ss
        fracstencil%t  = st1%t  / st2%t
        fracstencil%et = st1%et / st2%et
        fracstencil%wt = st1%wt / st2%wt
        fracstencil%nt = st1%nt / st2%nt
        fracstencil%st = st1%st / st2%st
        fracstencil%b  = st1%b  / st2%b
        fracstencil%eb = st1%eb / st2%eb
        fracstencil%wb = st1%wb / st2%wb
        fracstencil%nb = st1%nb / st2%nb
        fracstencil%sb = st1%sb / st2%sb
        fracstencil%tt = st1%tt / st2%tt
        fracstencil%bb = st1%bb / st2%bb
      
      END FUNCTION fracstencil
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
      FUNCTION nb_rank2( array, l, ijk )

      IMPLICIT NONE 
!
      TYPE(stencil) :: nb_rank2
      REAL*8, INTENT(IN) :: array(:,:)
      INTEGER, INTENT(IN) :: l, ijk
!
      IF( job_type_flag == 2 ) THEN

        nb_rank2%c  = array(l, ijk)
        nb_rank2%e  = array(l, ijr)
        nb_rank2%n  = array(l, ijt)
        nb_rank2%w  = array(l, ijl)
        nb_rank2%s  = array(l, ijb)
        nb_rank2%en  = array(l, ijtr)
        nb_rank2%wn  = array(l, ijtl)
        nb_rank2%es  = array(l, ijbr)
        nb_rank2%ws  = array(l, ijbl)
        nb_rank2%ee  = array(l, ijrr)
        nb_rank2%nn  = array(l, ijtt)
        nb_rank2%ww  = array(l, ijll)
        nb_rank2%ss  = array(l, ijbb)

      ELSE IF( job_type_flag == 3 ) THEN

        nb_rank2%c = array( l, ijk )
        nb_rank2%e = array( l, ijke )
        nb_rank2%w = array( l, ijkw )
        nb_rank2%ee = array( l, ijkee )
        nb_rank2%ww = array( l, ijkww )
        nb_rank2%n = array( l, ijkn )
        nb_rank2%en = array( l, ijken )
        nb_rank2%wn = array( l, ijkwn )
        nb_rank2%s = array( l, ijks )
        nb_rank2%es = array( l, ijkes )
        nb_rank2%ws = array( l, ijkws )
        nb_rank2%nn = array( l, ijknn )
        nb_rank2%ss = array( l, ijkss )
        nb_rank2%t = array( l, ijkt )
        nb_rank2%et = array( l, ijket )
        nb_rank2%wt = array( l, ijkwt )
        nb_rank2%nt = array( l, ijknt )
        nb_rank2%st = array( l, ijkst )
        nb_rank2%b = array( l, ijkb )
        nb_rank2%eb = array( l, ijkeb )
        nb_rank2%wb = array( l, ijkwb )
        nb_rank2%nb = array( l, ijknb )
        nb_rank2%sb = array( l, ijksb )
        nb_rank2%tt = array( l, ijktt )
        nb_rank2%bb = array( l, ijkbb )

      END IF

      RETURN
      END FUNCTION nb_rank2
!-----------------------------------------------------------------------
      FUNCTION rnb_rank2(array, l, ijk)
      USE dimensions
      IMPLICIT NONE 
!
      TYPE(stencil) :: rnb_rank2
      REAL*8, INTENT(IN) :: array(:, :)
      INTEGER, INTENT(IN) :: l, ijk
      INTEGER :: i,j,k,imesh
!
      imesh = myijk( ip0_jp0_kp0_, ijk)

      IF( job_type_flag == 2 ) THEN

        i  = MOD( ( imesh - 1 ), nr) + 1
        j  = ( imesh - 1 ) / nr + 1
!
        rnb_rank2%c  = array(l, ijk)
        rnb_rank2%e  = array(l, ipj)
        rnb_rank2%n  = array(l, ijp)
        rnb_rank2%w  = array(l, imj)
        rnb_rank2%s  = array(l, ijm)
        rnb_rank2%en  = array(l, ipjp)
        rnb_rank2%wn  = array(l, imjp)
        rnb_rank2%es  = array(l, ipjm)
        rnb_rank2%ws  = array(l, imjm)
        rnb_rank2%ee  = array(l, ippj)
        IF (i == (nr-1))  rnb_rank2%ee  = array(l, ipj)
        rnb_rank2%nn  = array(l, ijpp)
        IF (j == (nz-1))  rnb_rank2%nn  = array(l, ijp)
        rnb_rank2%ww  = array(l, immj)
        IF (i == 2)  rnb_rank2%ww  = array(l, imj)
        rnb_rank2%ss  = array(l, ijmm)
        IF (j == (2))  rnb_rank2%ss  = array(l, ijm)

      ELSE IF( job_type_flag == 3 ) THEN

         i = MOD( MOD( ijk - 1, nx*ny ), nx ) + 1
         j = MOD( ijk - 1, nx*ny ) / nx + 1
         k = ( ijk - 1 ) / ( nx*ny ) + 1

         rnb_rank2%c = array( l, ijk )
         rnb_rank2%e = array( l, ipjk )
         rnb_rank2%w = array( l, imjk )
         rnb_rank2%ee = array( l, ippjk )
         rnb_rank2%ww = array( l, immjk )
rnb_rank2%n = array( l, ijpk )
         rnb_rank2%en = array( l, ipjpk )
         rnb_rank2%wn = array( l, imjpk )
         rnb_rank2%s = array( l, ijmk )
         rnb_rank2%es = array( l, ipjmk )
         rnb_rank2%ws = array( l, imjmk )
         rnb_rank2%nn = array( l, ijppk )
         rnb_rank2%ss = array( l, ijmmk )
         rnb_rank2%t = array( l, ijkp )
         rnb_rank2%et = array( l, ipjkp )
         rnb_rank2%wt = array( l, imjkp )
         rnb_rank2%nt = array( l, ijpkp )
         rnb_rank2%st = array( l, ijmkp )
         rnb_rank2%b = array( l, ijkm )
         rnb_rank2%eb = array( l, ipjkm )
         rnb_rank2%wb = array( l, imjkm )
         rnb_rank2%nb = array( l, ijpkm )
         rnb_rank2%sb = array( l, ijmkm )
         rnb_rank2%tt = array( l, ijkpp )
         rnb_rank2%bb = array( l, ijkmm )

         IF (i == (nx-1))  rnb_rank2%ee  = array(l, ipjk)
         IF (i == 2)       rnb_rank2%ww  = array(l, imjk)
         IF (j == (ny-1))  rnb_rank2%nn  = array(l, ijpk)
         IF (j == 2)       rnb_rank2%ss  = array(l, ijmk)
         IF (k == (nz-1))  rnb_rank2%tt  = array(l, ijkp)
         IF (k == 2)       rnb_rank2%bb  = array(l, ijkm)

      END IF

      RETURN
      END FUNCTION rnb_rank2
!-----------------------------------------------------------------------
      END MODULE set_indexes
!-----------------------------------------------------------------------
