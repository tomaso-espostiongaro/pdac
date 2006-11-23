!-----------------------------------------------------------------------
      MODULE set_indexes
!-----------------------------------------------------------------------
      USE domain_mapping, ONLY: myinds, myijk
      USE indijk_module

      IMPLICIT NONE
!
! ... 2D/3D-indexes
!
      INTEGER :: ipjk, imjk, ippjk, immjk, ijpk, ipjpk, imjpk, ijmk,  &
                 ipjmk, imjmk, ijppk, ijmmk, ijkp, ipjkp, imjkp, ijpkp, &
                 ijmkp, ijkm, ipjkm, imjkm, ijpkm, ijmkm, ijkpp, ijkmm

      INTEGER ::  ijke, ijkw, ijkee, ijkww, ijkn, ijken, ijkwn, ijks, ijkes, &
                  ijkws, ijknn, ijkss, ijkt, ijket, ijkwt, ijknt, ijkst, ijkb, &
                  ijkeb, ijkwb, ijknb, ijksb, ijktt, ijkbb

      INTEGER, PRIVATE :: job_type_flag
!
! ... Define the computational stencil for Finite Volume schemes.
! ... The standard compass notation (east,west,north,south,top,bottom,etc.) 
! ... is adopted 
!
      TYPE stencil
        REAL*8 :: c
        REAL*8 :: e
        REAL*8 :: w
        REAL*8 :: n
        REAL*8 :: s
        REAL*8 :: b
        REAL*8 :: t
        REAL*8 :: ee
        REAL*8 :: ww
        REAL*8 :: en
        REAL*8 :: wn
        REAL*8 :: es
        REAL*8 :: ws
        REAL*8 :: nn
        REAL*8 :: ss
        REAL*8 :: et
        REAL*8 :: wt
        REAL*8 :: nt
        REAL*8 :: st
        REAL*8 :: eb
        REAL*8 :: wb
        REAL*8 :: nb
        REAL*8 :: sb
        REAL*8 :: tt
        REAL*8 :: bb
      END TYPE stencil
!
! ... Overloading of the arithmetic binary operators
! ... between stencils
!
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

        ijkm  = myijk( ip0_jp0_km1_, ijk )
        imjk  = myijk( im1_jp0_kp0_, ijk )
        ipjk  = myijk( ip1_jp0_kp0_, ijk )
        ijkp  = myijk( ip0_jp0_kp1_, ijk )
        ipjkm = myijk( ip1_jp0_km1_, ijk )
        ipjkp = myijk( ip1_jp0_kp1_, ijk )
        imjkm = myijk( im1_jp0_km1_, ijk )
        imjkp = myijk( im1_jp0_kp1_, ijk )
        ijkpp = myijk( ip0_jp0_kp2_, ijk )
        ippjk = myijk( ip2_jp0_kp0_, ijk )
        immjk = myijk( im2_jp0_kp0_, ijk )
        ijkmm = myijk( ip0_jp0_km2_, ijk )

        ijke  = myinds(ip1_jp0_kp0_, ijk )
        ijkt  = myinds(ip0_jp0_kp1_, ijk )
        ijkw  = myinds(im1_jp0_kp0_, ijk )
        ijkb  = myinds(ip0_jp0_km1_, ijk )
        ijket = myinds(ip1_jp0_kp1_, ijk )
        ijkwt = myinds(im1_jp0_kp1_, ijk )
        ijkeb = myinds(ip1_jp0_km1_, ijk )
        ijkwb = myinds(im1_jp0_km1_, ijk )
        ijkee = myinds(ip2_jp0_kp0_, ijk )
        ijktt = myinds(ip0_jp0_kp2_, ijk )
        ijkww = myinds(im2_jp0_kp0_, ijk )
        ijkbb = myinds(ip0_jp0_km2_, ijk )

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

    END SUBROUTINE subscr
!-----------------------------------------------------------------------
      SUBROUTINE first_subscr( ijk )
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: ijk
!
      IF( job_type_flag == 2 ) THEN

        ijkm  = myijk( ip0_jp0_km1_, ijk )
        imjk  = myijk( im1_jp0_kp0_, ijk )
        ipjk  = myijk( ip1_jp0_kp0_, ijk )
        ijkp  = myijk( ip0_jp0_kp1_, ijk )

        ijke  = myinds(ip1_jp0_kp0_, ijk )
        ijkt  = myinds(ip0_jp0_kp1_, ijk )
        ijkw  = myinds(im1_jp0_kp0_, ijk )
        ijkb  = myinds(ip0_jp0_km1_, ijk )

      ELSE IF( job_type_flag == 3 ) THEN

        ipjk   = myijk( ip1_jp0_kp0_ , ijk )
        imjk   = myijk( im1_jp0_kp0_ , ijk )
        ijpk   = myijk( ip0_jp1_kp0_ , ijk )
        ijmk   = myijk( ip0_jm1_kp0_ , ijk )
        ijkp   = myijk( ip0_jp0_kp1_ , ijk )
        ijkm   = myijk( ip0_jp0_km1_ , ijk )

        ijke = myinds( ip1_jp0_kp0_ , ijk )
        ijkw = myinds( im1_jp0_kp0_ , ijk )
        ijkn = myinds( ip0_jp1_kp0_ , ijk )
        ijks = myinds( ip0_jm1_kp0_ , ijk )
        ijkt = myinds( ip0_jp0_kp1_ , ijk )
        ijkb = myinds( ip0_jp0_km1_ , ijk )

      END IF

    END SUBROUTINE first_subscr
!-----------------------------------------------------------------------
      SUBROUTINE third_subscr( ijk )
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: ijk
!
      IF( job_type_flag == 2 ) THEN

        ippjk  = myijk( ip2_jp0_kp0_ , ijk )
        immjk  = myijk( im2_jp0_kp0_ , ijk )
        ijkpp  = myijk( ip0_jp0_kp2_ , ijk )
        ijkmm  = myijk( ip0_jp0_km2_ , ijk )

        ijkee = myinds( ip2_jp0_kp0_ , ijk )
        ijkww = myinds( im2_jp0_kp0_ , ijk )
        ijktt = myinds( ip0_jp0_kp2_ , ijk )
        ijkbb = myinds( ip0_jp0_km2_ , ijk )

      ELSE IF( job_type_flag == 3 ) THEN

        ippjk  = myijk( ip2_jp0_kp0_ , ijk )
        immjk  = myijk( im2_jp0_kp0_ , ijk )
        ijppk  = myijk( ip0_jp2_kp0_ , ijk )
        ijmmk  = myijk( ip0_jm2_kp0_ , ijk )
        ijkpp  = myijk( ip0_jp0_kp2_ , ijk )
        ijkmm  = myijk( ip0_jp0_km2_ , ijk )

        ijkee = myinds( ip2_jp0_kp0_ , ijk )
        ijkww = myinds( im2_jp0_kp0_ , ijk )
        ijknn = myinds( ip0_jp2_kp0_ , ijk )
        ijkss = myinds( ip0_jm2_kp0_ , ijk )
        ijktt = myinds( ip0_jp0_kp2_ , ijk )
        ijkbb = myinds( ip0_jp0_km2_ , ijk )

      END IF

    END SUBROUTINE third_subscr
!-----------------------------------------------------------------------

      FUNCTION cte(c)
! ... build a stencil with constant values c at each location
!
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
! ... sum two stencils
!
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
! ... subtract two stencils
!
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
! ... multiply two stencils, location by location
!
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
! ... multiply each stencil location by a constant value a

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
! ... compute the ratio between two stencils, location by location
!
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
      SUBROUTINE full_nb( stncl, array, ijk )
! ... This routine compute the complete stencil around a grid point
! ... considering the boundary conditions
!
      IMPLICIT NONE 
!
      TYPE(stencil) :: stncl
      REAL*8, INTENT(IN) :: array(:)
      INTEGER, INTENT(IN) :: ijk
!
      IF( job_type_flag == 2 ) THEN

        stncl%c  = array(ijk)

        stncl%e  = array(ijke)
        stncl%t  = array(ijkt)
        stncl%w  = array(ijkw)
        stncl%b  = array(ijkb)
        stncl%et  = array(ijket)
        stncl%wt  = array(ijkwt)
        stncl%eb  = array(ijkeb)
        stncl%wb  = array(ijkwb)
        stncl%ee  = array(ijkee)
        stncl%tt  = array(ijktt)
        stncl%ww  = array(ijkww)
        stncl%bb  = array(ijkbb)

      ELSE IF( job_type_flag == 3 ) THEN

        stncl%c = array( ijk )
        stncl%e = array( ijke )
        stncl%w = array( ijkw )
        stncl%ee = array( ijkee )
        stncl%ww = array( ijkww )
        stncl%n = array( ijkn )
        stncl%en = array( ijken )
        stncl%wn = array( ijkwn )
        stncl%s = array( ijks )
        stncl%es = array( ijkes )
        stncl%ws = array( ijkws )
        stncl%nn = array( ijknn )
        stncl%ss = array( ijkss )
        stncl%t = array( ijkt )
        stncl%et = array( ijket )
        stncl%wt = array( ijkwt )
        stncl%nt = array( ijknt )
        stncl%st = array( ijkst )
        stncl%b = array( ijkb )
        stncl%eb = array( ijkeb )
        stncl%wb = array( ijkwb )
        stncl%nb = array( ijknb )
        stncl%sb = array( ijksb )
        stncl%tt = array( ijktt )
        stncl%bb = array( ijkbb )

      END IF

      RETURN
      END SUBROUTINE full_nb
!-----------------------------------------------------------------------
      SUBROUTINE full_rnb(stncl,array,ijk)
! ... This routine compute the complete stencil around a grid point
! ... without considering the boundary conditions
!
      USE dimensions
      IMPLICIT NONE 
!
      TYPE(stencil) :: stncl
      REAL*8, INTENT(IN) :: array(:)
      INTEGER, INTENT(IN) :: ijk
!
      IF( job_type_flag == 2 ) THEN

        stncl%c  = array(ijk)
        stncl%e  = array(ipjk)
        stncl%t  = array(ijkp)
        stncl%w  = array(imjk)
        stncl%b  = array(ijkm)
        stncl%et  = array(ipjkp)
        stncl%wt  = array(imjkp)
        stncl%eb  = array(ipjkm)
        stncl%wb  = array(imjkm)
        stncl%ee  = array(ippjk)
        stncl%tt  = array(ijkpp)
        stncl%ww  = array(immjk)
        stncl%bb  = array(ijkmm)

      ELSE IF( job_type_flag == 3 ) THEN

         stncl%c = array( ijk )
         stncl%e = array( ipjk )
         stncl%w = array( imjk )
         stncl%ee = array( ippjk )
         stncl%ww = array( immjk )
         stncl%n = array( ijpk )
         stncl%en = array( ipjpk )
         stncl%wn = array( imjpk )
         stncl%s = array( ijmk )
         stncl%es = array( ipjmk )
         stncl%ws = array( imjmk )
         stncl%nn = array( ijppk )
         stncl%ss = array( ijmmk )
         stncl%t = array( ijkp )
         stncl%et = array( ipjkp )
         stncl%wt = array( imjkp )
         stncl%nt = array( ijpkp )
         stncl%st = array( ijmkp )
         stncl%b = array( ijkm )
         stncl%eb = array( ipjkm )
         stncl%wb = array( imjkm )
         stncl%nb = array( ijpkm )
         stncl%sb = array( ijmkm )
         stncl%tt = array( ijkpp )
         stncl%bb = array( ijkmm )

      END IF

      RETURN
      END SUBROUTINE full_rnb
!-----------------------------------------------------------------------
      SUBROUTINE first_nb( stncl, array, ijk )
! ... only first neighbours

      IMPLICIT NONE 
!
      TYPE(stencil) :: stncl
      REAL*8, INTENT(IN) :: array(:)
      INTEGER, INTENT(IN) :: ijk
!
      IF( job_type_flag == 2 ) THEN

        stncl%c  = array(ijk)
        stncl%e  = array(ijke)
        stncl%t  = array(ijkt)
        stncl%w  = array(ijkw)
        stncl%b  = array(ijkb)

      ELSE IF( job_type_flag == 3 ) THEN

        stncl%c = array( ijk )
        stncl%e = array( ijke )
        stncl%w = array( ijkw )
        stncl%n = array( ijkn )
        stncl%s = array( ijks )
        stncl%t = array( ijkt )
        stncl%b = array( ijkb )

      END IF

      RETURN
      END SUBROUTINE first_nb
!-----------------------------------------------------------------------
      SUBROUTINE third_nb( stncl, array, ijk )
! ... This routine compute the complete stencil around a grid point
! ... considering the boundary conditions
!
      IMPLICIT NONE 
!
      TYPE(stencil) :: stncl
      REAL*8, INTENT(IN) :: array(:)
      INTEGER, INTENT(IN) :: ijk
!
      IF( job_type_flag == 2 ) THEN

        stncl%ee  = array(ijkee)
        stncl%tt  = array(ijktt)
        stncl%ww  = array(ijkww)
        stncl%bb  = array(ijkbb)

      ELSE IF( job_type_flag == 3 ) THEN

        stncl%ee = array( ijkee )
        stncl%ww = array( ijkww )
        stncl%nn = array( ijknn )
        stncl%ss = array( ijkss )
        stncl%tt = array( ijktt )
        stncl%bb = array( ijkbb )

      END IF

      RETURN
      END SUBROUTINE third_nb
!-----------------------------------------------------------------------
      SUBROUTINE first_rnb(stncl,array,ijk)
! ... only first neighbours

      IMPLICIT NONE 
!
      TYPE(stencil) :: stncl
      REAL*8, INTENT(IN) :: array(:)
      INTEGER, INTENT(IN) :: ijk
!
      IF( job_type_flag == 2 ) THEN

        stncl%c  = array(ijk)
        stncl%e  = array(ipjk)
        stncl%t  = array(ijkp)
        stncl%w  = array(imjk)
        stncl%b  = array(ijkm)

      ELSE IF( job_type_flag == 3 ) THEN

        stncl%c = array( ijk )
        stncl%e = array( ipjk )
        stncl%w = array( imjk )
        stncl%n = array( ijpk )
        stncl%s = array( ijmk )
        stncl%t = array( ijkp )
        stncl%b = array( ijkm )

      END IF

      RETURN
      END SUBROUTINE first_rnb
!-----------------------------------------------------------------------
      SUBROUTINE third_rnb(stncl,array,ijk)
! ... This routine compute only the second order stencil elements
! ... without considering the boundary conditions
!
      USE dimensions
      IMPLICIT NONE 
!
      TYPE(stencil) :: stncl
      REAL*8, INTENT(IN) :: array(:)
      INTEGER, INTENT(IN) :: ijk
!
      IF( job_type_flag == 2 ) THEN

        stncl%ee  = array(ippjk)
        stncl%tt  = array(ijkpp)
        stncl%ww  = array(immjk)
        stncl%bb  = array(ijkmm)

      ELSE IF( job_type_flag == 3 ) THEN

         stncl%ee = array( ippjk )
         stncl%ww = array( immjk )
         stncl%nn = array( ijppk )
         stncl%ss = array( ijmmk )
         stncl%tt = array( ijkpp )
         stncl%bb = array( ijkmm )

      END IF

      RETURN
      END SUBROUTINE third_rnb
!-----------------------------------------------------------------------
      SUBROUTINE nb( stncl, array, ijk )
! ... This routine compute the complete stencil around a grid point
! ... considering the boundary conditions
!
      IMPLICIT NONE 
!
      TYPE(stencil) :: stncl
      REAL*8, INTENT(IN) :: array(:)
      INTEGER, INTENT(IN) :: ijk
!
      IF( job_type_flag == 2 ) THEN

        stncl%c  = array(ijk)
        stncl%e  = array(ijke)
        stncl%t  = array(ijkt)
        stncl%w  = array(ijkw)
        stncl%b  = array(ijkb)
        stncl%et  = array(ijket)
        stncl%wt  = array(ijkwt)
        stncl%eb  = array(ijkeb)
        stncl%ee  = array(ijkee)
        stncl%tt  = array(ijktt)

      ELSE IF( job_type_flag == 3 ) THEN

        stncl%c = array( ijk )
        stncl%e = array( ijke )
        stncl%w = array( ijkw )
        stncl%ee = array( ijkee )
        stncl%n = array( ijkn )
        stncl%en = array( ijken )
        stncl%wn = array( ijkwn )
        stncl%s = array( ijks )
        stncl%es = array( ijkes )
        stncl%nn = array( ijknn )
        stncl%t = array( ijkt )
        stncl%et = array( ijket )
        stncl%wt = array( ijkwt )
        stncl%nt = array( ijknt )
        stncl%st = array( ijkst )
        stncl%b = array( ijkb )
        stncl%eb = array( ijkeb )
        stncl%nb = array( ijknb )
        stncl%tt = array( ijktt )

      END IF

      RETURN
      END SUBROUTINE nb
!-----------------------------------------------------------------------
      SUBROUTINE rnb(stncl,array,ijk)
! ... This routine compute the complete stencil around a grid point
! ... without considering the boundary conditions
!
      USE dimensions
      IMPLICIT NONE 
!
      TYPE(stencil) :: stncl
      REAL*8, INTENT(IN) :: array(:)
      INTEGER, INTENT(IN) :: ijk
!
      IF( job_type_flag == 2 ) THEN

        stncl%c  = array(ijk)
        stncl%e  = array(ipjk)
        stncl%ee = array( ippjk )
        stncl%t  = array(ijkp)
        stncl%w  = array(imjk)
        stncl%b  = array(ijkm)
        stncl%wt  = array(imjkp)
        stncl%eb  = array(ipjkm)
        stncl%tt  = array( ijkpp )

      ELSE IF( job_type_flag == 3 ) THEN

        stncl%c = array( ijk )
        stncl%e = array( ipjk )
        stncl%w = array( imjk )
        stncl%ee = array( ippjk )
        stncl%n = array( ijpk )
        stncl%wn = array( imjpk )
        stncl%nn = array( ijppk )
        stncl%s = array( ijmk )
        stncl%es = array( ipjmk )
        stncl%t = array( ijkp )
        stncl%wt = array( imjkp )
        stncl%st = array( ijmkp )
        stncl%b = array( ijkm )
        stncl%eb = array( ipjkm )
        stncl%nb = array( ijpkm )
        stncl%tt = array( ijkpp )

      END IF

      RETURN
      END SUBROUTINE rnb
!-----------------------------------------------------------------------
      SUBROUTINE check_stencil(ijk,stx,sty,stz)
      USE immersed_boundaries, ONLY: numx, numy, numz
      IMPLICIT NONE
      TYPE(stencil), INTENT(INOUT) :: stx, sty, stz
      INTEGER, INTENT(IN) :: ijk
!
       IF (numx( ijk )>0) stx%c = 0.D0
       IF (numx( ipjk )>0) stx%e = 0.D0
       IF (numx( imjk )>0) stx%w = 0.D0
       IF (numx( ijpk )>0) stx%n = 0.D0
       IF (numx( imjpk )>0) stx%wn = 0.D0
       IF (numx( ijmk )>0) stx%s = 0.D0
       IF (numx( ipjmk )>0) stx%es = 0.D0
       IF (numx( ijkp )>0) stx%t = 0.D0
       IF (numx( imjkp )>0) stx%wt = 0.D0
       IF (numx( ijmkp )>0) stx%st = 0.D0
       IF (numx( ijkm )>0) stx%b = 0.D0
       IF (numx( ipjkm )>0) stx%eb = 0.D0
       IF (numx( ijpkm )>0) stx%nb = 0.D0
!     
       IF (numy( ijk )>0) sty%c = 0.D0
       IF (numy( ipjk )>0) sty%e = 0.D0
       IF (numy( imjk )>0) sty%w = 0.D0
       IF (numy( ijpk )>0) sty%n = 0.D0
       IF (numy( imjpk )>0) sty%wn = 0.D0
       IF (numy( ijmk )>0) sty%s = 0.D0
       IF (numy( ipjmk )>0) sty%es = 0.D0
       IF (numy( ijkp )>0) sty%t = 0.D0
       IF (numy( imjkp )>0) sty%wt = 0.D0
       IF (numy( ijmkp )>0) sty%st = 0.D0
       IF (numy( ijkm )>0) sty%b = 0.D0
       IF (numy( ipjkm )>0) sty%eb = 0.D0
       IF (numy( ijpkm )>0) sty%nb = 0.D0
!
       IF (numz( ijk )>0) stz%c = 0.D0
       IF (numz( ipjk )>0) stz%e = 0.D0
       IF (numz( imjk )>0) stz%w = 0.D0
       IF (numz( ijpk )>0) stz%n = 0.D0
       IF (numz( imjpk )>0) stz%wn = 0.D0
       IF (numz( ijmk )>0) stz%s = 0.D0
       IF (numz( ipjmk )>0) stz%es = 0.D0
       IF (numz( ijkp )>0) stz%t = 0.D0
       IF (numz( imjkp )>0) stz%wt = 0.D0
       IF (numz( ijmkp )>0) stz%st = 0.D0
       IF (numz( ijkm )>0) stz%b = 0.D0
       IF (numz( ipjkm )>0) stz%eb = 0.D0
       IF (numz( ijpkm )>0) stz%nb = 0.D0
!
      RETURN
      END SUBROUTINE check_stencil
!-----------------------------------------------------------------------
      END MODULE set_indexes
!-----------------------------------------------------------------------
