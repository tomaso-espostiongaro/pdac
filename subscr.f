
!-----------------------------------------------------------------------
      MODULE set_indexes
!-----------------------------------------------------------------------
      USE domain_mapping, ONLY: myinds, myijk
      USE indijk_module

      IMPLICIT NONE
!
! ... 2D/3D-indexes
!
      INTEGER :: ipjk, imjk, ijpk, ijmk, ijkp, ijkm, &
                 ippjk, immjk, ijppk, ijmmk, ijkpp, ijkmm, &
                 ipjpk, ipjmk, ipjkp, ipjkm, & 
                 imjpk, imjmk, imjkp, imjkm, &
                 ijpkp, ijpkm, ijmkp, ijmkm, &
                 ipjpkp, ipjpkm, ipjmkp, ipjmkm, imjpkp, imjpkm, imjmkp, imjmkm, &
                 ippjkp, ippjkm, ippjpk, ippjmk, immjkp, immjkm, immjpk, immjmk, &
                 ipjppk, imjppk, ijppkp, ijppkm, ipjmmk, imjmmk, ijmmkp, ijmmkm, &
                 ipjkpp, imjkpp, ijpkpp, ijmkpp, ipjkmm, imjkmm, ijpkmm, ijmkmm, &
                 ippjpkp, ippjpkm, ippjmkp, ippjmkm, immjpkp, immjpkm, immjmkp, immjmkm, &
                 ipjppkp, ipjppkm, imjppkp, imjppkm, ipjmmkp, ipjmmkm, imjmmkp, imjmmkm, &
                 ipjpkpp, ipjmkpp, imjpkpp, imjmkpp, ipjpkmm, ipjmkmm, imjpkmm, imjmkmm, &
                 ippjppk, ippjmmk, ippjkpp, &
                 ijppkpp, ijppkmm, immjppk, &
                 immjkpp, ijmmkpp, ippjkmm, &
                 ipppjk, ipppjpk, ipppjmk, ipppjkp, ipppjpkp, ipppjmkp, &
                 ijpppk, ipjpppk, ijpppkp, ijpppkm, ipjpppkp, ipjpppkm, &
                 ijkppp, ipjkppp, imjkppp, ijpkppp, ipjpkppp, imjpkppp, &
                 ippjppkp, ippjppkm, ippjmmkp, ippjmmkm, ippjpkmm, ippjmkmm, &
                 ipjppkpp, imjppkpp, ipjppkmm, imjppkmm, immjppkp, immjppkm, &
                 ippjpkpp, ippjmkpp, ipjmmkpp, immjpkpp, immjmkpp, imjmmkpp
                 

      INTEGER :: ijke, ijkw, ijkn, ijks, ijkt, ijkb, &
                 ijkee, ijkww, ijknn, ijkss, ijktt, ijkbb, &
                 ijken, ijkes, ijket, ijkeb, &
                 ijkwn, ijkws, ijkwt, ijkwb, &
                 ijknt, ijknb, ijkst, ijksb, &
                 ijkent, ijkenb, ijkest, ijkesb, ijkwnt, ijkwnb, ijkwst, ijkwsb, &
                 ijkeet, ijkeeb, ijkeen, ijkees, ijkwwt, ijkwwb, ijkwwn, ijkwws, &
                 ijkenn, ijkwnn, ijknnt, ijknnb, ijkess, ijkwss, ijksst, ijkssb, &
                 ijkett, ijkwtt, ijkntt, ijkstt, ijkebb, ijkwbb, ijknbb, ijksbb, &
                 ijkeent, ijkeenb, ijkeest, ijkeesb, ijkwwnt, ijkwwnb, ijkwwst, ijkwwsb, &
                 ijkennt, ijkennb, ijkwnnt, ijkwnnb, ijkesst, ijkessb, ijkwsst, ijkwssb, &
                 ijkentt, ijkestt, ijkwntt, ijkwstt, ijkenbb, ijkesbb, ijkwnbb, ijkwsbb, &
                 ijkeenn, ijkeess, ijkeett, &
                 ijknntt, ijknnbb, ijkwwnn, &
                 ijkwwtt, ijksstt, ijkeebb, &
                 ijkeee, ijkeeen, ijkeees, ijkeeet, ijkeeent, ijkeeest, &
                 ijknnn, ijkennn, ijknnnt, ijknnnb, ijkennnt, ijkennnb, &
                 ijkttt, ijkettt, ijkwttt, ijknttt, ijkenttt, ijkwnttt, &
                 ijkeennt, ijkeennb, ijkeesst, ijkeessb, ijkeenbb, ijkeesbb, &
                 ijkenntt, ijkwnntt, ijkennbb, ijkwnnbb, ijkwwnnt, ijkwwnnb,  &
                 ijkeentt, ijkeestt, ijkesstt, ijkwwntt, ijkwwstt, ijkwsstt

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
        REAL*8 :: nn
        REAL*8 :: ss
        REAL*8 :: bb
        REAL*8 :: tt
        REAL*8 :: et
        REAL*8 :: eb
        REAL*8 :: en
        REAL*8 :: es
        REAL*8 :: wt
        REAL*8 :: wb
        REAL*8 :: wn
        REAL*8 :: ws
        REAL*8 :: nt
        REAL*8 :: st
        REAL*8 :: nb
        REAL*8 :: sb
        REAL*8 :: wst 
        REAL*8 :: wsb
        REAL*8 :: wnt
        REAL*8 :: wnb
        REAL*8 :: est
        REAL*8 :: esb
        REAL*8 :: ent
        REAL*8 :: enb
        REAL*8 :: eet
        REAL*8 :: eeb
        REAL*8 :: een
        REAL*8 :: ees
        REAL*8 :: wwt
        REAL*8 :: wwn
        REAL*8 :: sst
        REAL*8 :: ess
        REAL*8 :: nnt
        REAL*8 :: nnb
        REAL*8 :: enn
        REAL*8 :: wnn
        REAL*8 :: ett
        REAL*8 :: wtt
        REAL*8 :: stt
        REAL*8 :: ntt
        REAL*8 :: ebb
        REAL*8 :: nbb
        REAL*8 :: eent 
        REAL*8 :: eenb
        REAL*8 :: eest
        REAL*8 :: eesb
        REAL*8 :: wwnt
        REAL*8 :: ennt
        REAL*8 :: ennb
        REAL*8 :: wnnt
        REAL*8 :: wnnb
        REAL*8 :: esst
        REAL*8 :: entt
        REAL*8 :: estt
        REAL*8 :: wntt
        REAL*8 :: wstt
        REAL*8 :: enbb
        REAL*8 :: wbb
        REAL*8 :: ssb
        REAL*8 :: wws
        REAL*8 :: wnbb
        REAL*8 :: essb
        REAL*8 :: wwst
        REAL*8 :: esbb
        REAL*8 :: wsbb
        REAL*8 :: sbb
        REAL*8 :: wssb
        REAL*8 :: wsst
        REAL*8 :: wss
        REAL*8 :: wwnb
        REAL*8 :: wwb
        REAL*8 :: wwsb
        REAL*8 :: enttt
        REAL*8 :: wnttt
        REAL*8 :: nttt
        REAL*8 :: ettt
        REAL*8 :: wttt
        REAL*8 :: ttt
        REAL*8 :: esstt
        REAL*8 :: sstt
        REAL*8 :: enntt
        REAL*8 :: nntt
        REAL*8 :: wwstt
        REAL*8 :: wwntt
        REAL*8 :: wwtt
        REAL*8 :: eestt
        REAL*8 :: eentt
        REAL*8 :: eett
        REAL*8 :: ennbb
        REAL*8 :: wnnbb
        REAL*8 :: nnbb
        REAL*8 :: wnntt
        REAL*8 :: ennnb
        REAL*8 :: ennnt
        REAL*8 :: nnnb
        REAL*8 :: nnnt
        REAL*8 :: ennn
        REAL*8 :: wwnnt
        REAL*8 :: wwnn
        REAL*8 :: eennt
        REAL*8 :: eenn
        REAL*8 :: nnn
        REAL*8 :: eenbb
        REAL*8 :: eebb
        REAL*8 :: eessb
        REAL*8 :: eesst
        REAL*8 :: eess
        REAL*8 :: eennb
        REAL*8 :: eeest
        REAL*8 :: eeent
        REAL*8 :: eeet
        REAL*8 :: eees
        REAL*8 :: eeen
        REAL*8 :: eee
        REAL*8 :: wwbb
        REAL*8 :: wwss
        REAL*8 :: ssbb
        REAL*8 :: eeeb
        REAL*8 :: eeenb
        REAL*8 :: eeesb
        REAL*8 :: wnnn
        REAL*8 :: wnnnt
        REAL*8 :: wnnnb
        REAL*8 :: sttt
        REAL*8 :: esttt
        REAL*8 :: wsttt
        REAL*8 :: wwnnb
        REAL*8 :: wwsst
        REAL*8 :: wwssb
        REAL*8 :: eesbb
        REAL*8 :: wwnbb
        REAL*8 :: wwsbb
        REAL*8 :: wsstt
        REAL*8 :: essbb
        REAL*8 :: wssbb

      END TYPE stencil
!
      TYPE masks
        INTEGER :: c
        INTEGER :: e
        INTEGER :: w
        INTEGER :: n
        INTEGER :: s
        INTEGER :: b
        INTEGER :: t
        INTEGER :: ee
        INTEGER :: ww
        INTEGER :: nn
        INTEGER :: ss
        INTEGER :: bb
        INTEGER :: tt
        INTEGER :: et
        INTEGER :: eb
        INTEGER :: en
        INTEGER :: es
        INTEGER :: wt
        INTEGER :: wb
        INTEGER :: wn
        INTEGER :: ws
        INTEGER :: nt
        INTEGER :: st
        INTEGER :: nb
        INTEGER :: sb
        INTEGER :: wst 
        INTEGER :: wsb
        INTEGER :: wnt
        INTEGER :: wnb
        INTEGER :: est
        INTEGER :: esb
        INTEGER :: ent
        INTEGER :: enb
        INTEGER :: eet
        INTEGER :: eeb
        INTEGER :: een
        INTEGER :: ees
        INTEGER :: wwt
        INTEGER :: wwn
        INTEGER :: sst
        INTEGER :: ess
        INTEGER :: nnt
        INTEGER :: nnb
        INTEGER :: enn
        INTEGER :: wnn
        INTEGER :: ett
        INTEGER :: wtt
        INTEGER :: stt
        INTEGER :: ntt
        INTEGER :: ebb
        INTEGER :: nbb
        INTEGER :: eent 
        INTEGER :: eenb
        INTEGER :: eest
        INTEGER :: eesb
        INTEGER :: wwnt
        INTEGER :: ennt
        INTEGER :: ennb
        INTEGER :: wnnt
        INTEGER :: wnnb
        INTEGER :: esst
        INTEGER :: entt
        INTEGER :: estt
        INTEGER :: wntt
        INTEGER :: wstt
        INTEGER :: enbb
        INTEGER :: wbb
        INTEGER :: ssb
        INTEGER :: wws
        INTEGER :: wnbb
        INTEGER :: essb
        INTEGER :: wwst
        INTEGER :: esbb
        INTEGER :: wsbb
        INTEGER :: sbb
        INTEGER :: wssb
        INTEGER :: wsst
        INTEGER :: wss
        INTEGER :: wwnb
        INTEGER :: wwb
        INTEGER :: wwsb
        INTEGER :: enttt
        INTEGER :: wnttt
        INTEGER :: nttt
        INTEGER :: ettt
        INTEGER :: wttt
        INTEGER :: ttt
        INTEGER :: esstt
        INTEGER :: sstt
        INTEGER :: enntt
        INTEGER :: nntt
        INTEGER :: wwstt
        INTEGER :: wwntt
        INTEGER :: wwtt
        INTEGER :: eestt
        INTEGER :: eentt
        INTEGER :: eett
        INTEGER :: ennbb
        INTEGER :: wnnbb
        INTEGER :: nnbb
        INTEGER :: wnntt
        INTEGER :: ennnb
        INTEGER :: ennnt
        INTEGER :: nnnb
        INTEGER :: nnnt
        INTEGER :: ennn
        INTEGER :: wwnnt
        INTEGER :: wwnn
        INTEGER :: eennt
        INTEGER :: eenn
        INTEGER :: nnn
        INTEGER :: eenbb
        INTEGER :: eebb
        INTEGER :: eessb
        INTEGER :: eesst
        INTEGER :: eess
        INTEGER :: eennb
        INTEGER :: eeest
        INTEGER :: eeent
        INTEGER :: eeet
        INTEGER :: eees
        INTEGER :: eeen
        INTEGER :: eee
        INTEGER :: wwbb
        INTEGER :: wwss
        INTEGER :: ssbb
        INTEGER :: eeeb
        INTEGER :: eeenb
        INTEGER :: eeesb
        INTEGER :: wnnn
        INTEGER :: wnnnt
        INTEGER :: wnnnb
        INTEGER :: sttt
        INTEGER :: esttt
        INTEGER :: wsttt
        INTEGER :: wwnnb
        INTEGER :: wwsst
        INTEGER :: wwssb
        INTEGER :: eesbb
        INTEGER :: wwnbb
        INTEGER :: wwsbb
        INTEGER :: wsstt
        INTEGER :: essbb
        INTEGER :: wssbb
      END TYPE masks
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
        MODULE PROCEDURE prodstencil, dotstencil, maskstencil
      END INTERFACE
      INTERFACE OPERATOR(/)
        MODULE PROCEDURE fracstencil
      END INTERFACE
      INTERFACE rnb
          MODULE PROCEDURE rnb_r, rnb_i
      END INTERFACE
!
      SAVE
!-----------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------
      SUBROUTINE subsc_setup( job_type )
        USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: job_type
        IF( job_type == JOB_TYPE_2D  ) THEN
          job_type_flag = 2
        ELSE IF( job_type == JOB_TYPE_3D ) THEN
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

        ijkm   = myijk( ip0_jp0_km1_, ijk )
        imjk   = myijk( im1_jp0_kp0_, ijk )
        ipjk   = myijk( ip1_jp0_kp0_, ijk )
        ijkp   = myijk( ip0_jp0_kp1_, ijk )
        ipjkm  = myijk( ip1_jp0_km1_, ijk )
        ipjkp  = myijk( ip1_jp0_kp1_, ijk )
        imjkm  = myijk( im1_jp0_km1_, ijk )
        imjkp  = myijk( im1_jp0_kp1_, ijk )
        ijkpp  = myijk( ip0_jp0_kp2_, ijk )
        ippjk  = myijk( ip2_jp0_kp0_, ijk )
        immjk  = myijk( im2_jp0_kp0_, ijk )
        ijkmm  = myijk( ip0_jp0_km2_, ijk )
        ippjkp = myijk( ip2_jp0_kp1_, ijk )
        ipjkpp = myijk( ip1_jp0_kp2_, ijk )
        ipppjk = myijk( ip3_jp0_kp0_, ijk )
        ijkppp = myijk( ip0_jp0_kp3_, ijk )

        ijke   = myinds( ip1_jp0_kp0_, ijk )
        ijkt   = myinds( ip0_jp0_kp1_, ijk )
        ijkw   = myinds( im1_jp0_kp0_, ijk )
        ijkb   = myinds( ip0_jp0_km1_, ijk )
        ijket  = myinds( ip1_jp0_kp1_, ijk )
        ijkwt  = myinds( im1_jp0_kp1_, ijk )
        ijkeb  = myinds( ip1_jp0_km1_, ijk )
        ijkwb  = myinds( im1_jp0_km1_, ijk )
        ijkee  = myinds( ip2_jp0_kp0_, ijk )
        ijktt  = myinds( ip0_jp0_kp2_, ijk )
        ijkww  = myinds( im2_jp0_kp0_, ijk )
        ijkbb  = myinds( ip0_jp0_km2_, ijk )
        ijkeet = myinds( ip2_jp0_kp1_, ijk )
        ijkett = myinds( ip1_jp0_kp2_, ijk )
        ijkeee = myinds( ip3_jp0_kp0_, ijk )
        ijkttt = myinds( ip0_jp0_kp3_, ijk )

      ELSE IF( job_type_flag == 3 ) THEN

        ipjk   = myijk( ip1_jp0_kp0_ , ijk )
        imjk   = myijk( im1_jp0_kp0_ , ijk )
        ijpk   = myijk( ip0_jp1_kp0_ , ijk )
        ijmk   = myijk( ip0_jm1_kp0_ , ijk )
        ijkp   = myijk( ip0_jp0_kp1_ , ijk )
        ijkm   = myijk( ip0_jp0_km1_ , ijk )
        ippjk  = myijk( ip2_jp0_kp0_ , ijk )
        immjk  = myijk( im2_jp0_kp0_ , ijk )
        ijppk  = myijk( ip0_jp2_kp0_ , ijk )
        ijmmk  = myijk( ip0_jm2_kp0_ , ijk )
        ijkpp  = myijk( ip0_jp0_kp2_ , ijk )
        ijkmm  = myijk( ip0_jp0_km2_ , ijk )
        ipjpk  = myijk( ip1_jp1_kp0_ , ijk )
        imjpk  = myijk( im1_jp1_kp0_ , ijk )
        ipjmk  = myijk( ip1_jm1_kp0_ , ijk )
        imjmk  = myijk( im1_jm1_kp0_ , ijk )
        ipjkp  = myijk( ip1_jp0_kp1_ , ijk )
        imjkp  = myijk( im1_jp0_kp1_ , ijk )
        ijpkp  = myijk( ip0_jp1_kp1_ , ijk )
        ijmkp  = myijk( ip0_jm1_kp1_ , ijk )
        ipjkm  = myijk( ip1_jp0_km1_ , ijk )
        imjkm  = myijk( im1_jp0_km1_ , ijk )
        ijpkm  = myijk( ip0_jp1_km1_ , ijk )
        ijmkm  = myijk( ip0_jm1_km1_ , ijk )
        ipppjk = myijk( ip3_jp0_kp0_ , ijk )
        ijpppk = myijk( ip0_jp3_kp0_ , ijk )
        ijkppp = myijk( ip0_jp0_kp3_ , ijk )
        ippjpk = myijk( ip2_jp1_kp0_ , ijk )
        ippjkp = myijk( ip2_jp0_kp1_ , ijk )
        ipjppk = myijk( ip1_jp2_kp0_ , ijk )
        ijppkp = myijk( ip0_jp2_kp1_ , ijk )
        ipjkpp = myijk( ip1_jp0_kp2_ , ijk )
        ijpkpp = myijk( ip0_jp1_kp2_ , ijk )

        ijke   = myinds( ip1_jp0_kp0_ , ijk )
        ijkw   = myinds( im1_jp0_kp0_ , ijk )
        ijkn   = myinds( ip0_jp1_kp0_ , ijk )
        ijks   = myinds( ip0_jm1_kp0_ , ijk )
        ijkt   = myinds( ip0_jp0_kp1_ , ijk )
        ijkb   = myinds( ip0_jp0_km1_ , ijk )
        ijkee  = myinds( ip2_jp0_kp0_ , ijk )
        ijkww  = myinds( im2_jp0_kp0_ , ijk )
        ijknn  = myinds( ip0_jp2_kp0_ , ijk )
        ijkss  = myinds( ip0_jm2_kp0_ , ijk )
        ijktt  = myinds( ip0_jp0_kp2_ , ijk )
        ijkbb  = myinds( ip0_jp0_km2_ , ijk )
        ijken  = myinds( ip1_jp1_kp0_ , ijk )
        ijkwn  = myinds( im1_jp1_kp0_ , ijk )
        ijkes  = myinds( ip1_jm1_kp0_ , ijk )
        ijkws  = myinds( im1_jm1_kp0_ , ijk )
        ijket  = myinds( ip1_jp0_kp1_ , ijk )
        ijkwt  = myinds( im1_jp0_kp1_ , ijk )
        ijknt  = myinds( ip0_jp1_kp1_ , ijk )
        ijkst  = myinds( ip0_jm1_kp1_ , ijk )
        ijkeb  = myinds( ip1_jp0_km1_ , ijk )
        ijkwb  = myinds( im1_jp0_km1_ , ijk )
        ijknb  = myinds( ip0_jp1_km1_ , ijk )
        ijksb  = myinds( ip0_jm1_km1_ , ijk )
        ijkeee = myinds( ip3_jp0_kp0_ , ijk )
        ijknnn = myinds( ip0_jp3_kp0_ , ijk )
        ijkttt = myinds( ip0_jp0_kp3_ , ijk )
        ijkeen = myinds( ip2_jp1_kp0_ , ijk )
        ijkeet = myinds( ip2_jp0_kp1_ , ijk )
        ijkenn = myinds( ip1_jp2_kp0_ , ijk )
        ijknnt = myinds( ip0_jp2_kp1_ , ijk )
        ijkett = myinds( ip1_jp0_kp2_ , ijk )
        ijkntt = myinds( ip0_jp1_kp2_ , ijk )

      END IF

      END SUBROUTINE subscr
!-----------------------------------------------------------------------
      SUBROUTINE first_subscr( ijk )
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: ijk
!
      IF( job_type_flag == 2 ) THEN

        ijkm = myijk( ip0_jp0_km1_, ijk )
        imjk = myijk( im1_jp0_kp0_, ijk )
        ipjk = myijk( ip1_jp0_kp0_, ijk )
        ijkp = myijk( ip0_jp0_kp1_, ijk )

        ijke = myinds( ip1_jp0_kp0_, ijk )
        ijkt = myinds( ip0_jp0_kp1_, ijk )
        ijkw = myinds( im1_jp0_kp0_, ijk )
        ijkb = myinds( ip0_jp0_km1_, ijk )

      ELSE IF( job_type_flag == 3 ) THEN

        ipjk = myijk( ip1_jp0_kp0_ , ijk )
        imjk = myijk( im1_jp0_kp0_ , ijk )
        ijpk = myijk( ip0_jp1_kp0_ , ijk )
        ijmk = myijk( ip0_jm1_kp0_ , ijk )
        ijkp = myijk( ip0_jp0_kp1_ , ijk )
        ijkm = myijk( ip0_jp0_km1_ , ijk )

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

        ippjk = myijk( ip2_jp0_kp0_ , ijk )
        immjk = myijk( im2_jp0_kp0_ , ijk )
        ijkpp = myijk( ip0_jp0_kp2_ , ijk )
        ijkmm = myijk( ip0_jp0_km2_ , ijk )

        ijkee = myinds( ip2_jp0_kp0_ , ijk )
        ijkww = myinds( im2_jp0_kp0_ , ijk )
        ijktt = myinds( ip0_jp0_kp2_ , ijk )
        ijkbb = myinds( ip0_jp0_km2_ , ijk )

      ELSE IF( job_type_flag == 3 ) THEN

        ippjk = myijk( ip2_jp0_kp0_ , ijk )
        immjk = myijk( im2_jp0_kp0_ , ijk )
        ijppk = myijk( ip0_jp2_kp0_ , ijk )
        ijmmk = myijk( ip0_jm2_kp0_ , ijk )
        ijkpp = myijk( ip0_jp0_kp2_ , ijk )
        ijkmm = myijk( ip0_jp0_km2_ , ijk )

        ijkee = myinds( ip2_jp0_kp0_ , ijk )
        ijkww = myinds( im2_jp0_kp0_ , ijk )
        ijknn = myinds( ip0_jp2_kp0_ , ijk )
        ijkss = myinds( ip0_jm2_kp0_ , ijk )
        ijktt = myinds( ip0_jp0_kp2_ , ijk )
        ijkbb = myinds( ip0_jp0_km2_ , ijk )

      END IF

      END SUBROUTINE third_subscr
!-----------------------------------------------------------------------
      SUBROUTINE ctu1_subscr( ijk )
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ijk

      IF( job_type_flag == 2 ) THEN

        ijkm   = myijk( ip0_jp0_km1_, ijk )
        imjk   = myijk( im1_jp0_kp0_, ijk )
        ipjk   = myijk( ip1_jp0_kp0_, ijk )
        ijkp   = myijk( ip0_jp0_kp1_, ijk )
        ipjkm  = myijk( ip1_jp0_km1_, ijk )
        ipjkp  = myijk( ip1_jp0_kp1_, ijk )
        imjkm  = myijk( im1_jp0_km1_, ijk )
        imjkp  = myijk( im1_jp0_kp1_, ijk )
        ijkpp  = myijk( ip0_jp0_kp2_, ijk )
        ijkmm  = myijk( ip0_jp0_km2_, ijk )
        ippjk  = myijk( ip2_jp0_kp0_, ijk )
        immjk  = myijk( im2_jp0_kp0_, ijk )
        ippjkp = myijk( ip2_jp0_kp1_, ijk )
        ippjkm = myijk( ip2_jp0_km1_, ijk )
        ipjkpp = myijk( ip1_jp0_kp2_, ijk )
        imjkpp = myijk( im1_jp0_kp2_, ijk )
        immjkm = myijk( im2_jp0_km1_, ijk )
        immjkp = myijk( im2_jp0_kp1_, ijk )
        ipjkmm = myijk( ip1_jp0_km2_, ijk )
        imjkmm = myijk( im1_jp0_km2_, ijk )

        ijke   = myinds( ip1_jp0_kp0_, ijk )
        ijkt   = myinds( ip0_jp0_kp1_, ijk )
        ijkw   = myinds( im1_jp0_kp0_, ijk )
        ijkb   = myinds( ip0_jp0_km1_, ijk )
        ijket  = myinds( ip1_jp0_kp1_, ijk )
        ijkwt  = myinds( im1_jp0_kp1_, ijk )
        ijkeb  = myinds( ip1_jp0_km1_, ijk )
        ijkwb  = myinds( im1_jp0_km1_, ijk )
        ijkee  = myinds( ip2_jp0_kp0_, ijk )
        ijkww  = myinds( im2_jp0_kp0_, ijk )
        ijktt  = myinds( ip0_jp0_kp2_, ijk )
        ijkbb  = myinds( ip0_jp0_km2_, ijk )
        ijkeeb = myinds( ip2_jp0_km1_, ijk )
        ijkeet = myinds( ip2_jp0_kp1_, ijk )
        ijkett = myinds( ip1_jp0_kp2_, ijk )
        ijkwtt = myinds( im1_jp0_kp2_, ijk )
        ijkwwb = myinds( im2_jp0_km1_, ijk )
        ijkwwt = myinds( im2_jp0_kp1_, ijk )
        ijkebb = myinds( ip1_jp0_km2_, ijk )
        ijkwbb = myinds( im1_jp0_km2_, ijk )

      ELSE

        ipjk     = myijk( ip1_jp0_kp0_ , ijk )
        imjk     = myijk( im1_jp0_kp0_ , ijk )
        ijpk     = myijk( ip0_jp1_kp0_ , ijk )
        ijmk     = myijk( ip0_jm1_kp0_ , ijk )
        ijkp     = myijk( ip0_jp0_kp1_ , ijk )
        ijkm     = myijk( ip0_jp0_km1_ , ijk )
        ipjpk    = myijk( ip1_jp1_kp0_ , ijk )
        imjpk    = myijk( im1_jp1_kp0_ , ijk )
        ipjmk    = myijk( ip1_jm1_kp0_ , ijk )
        imjmk    = myijk( im1_jm1_kp0_ , ijk )
        ipjkp    = myijk( ip1_jp0_kp1_ , ijk )
        imjkp    = myijk( im1_jp0_kp1_ , ijk )
        ijpkp    = myijk( ip0_jp1_kp1_ , ijk )
        ijmkp    = myijk( ip0_jm1_kp1_ , ijk )
        ipjkm    = myijk( ip1_jp0_km1_ , ijk )
        imjkm    = myijk( im1_jp0_km1_ , ijk )
        ijpkm    = myijk( ip0_jp1_km1_ , ijk )
        ijmkm    = myijk( ip0_jm1_km1_ , ijk )
        imjmkp   = myijk( im1_jm1_kp1_ , ijk )
        ipjmkp   = myijk( ip1_jm1_kp1_ , ijk )
        imjpkp   = myijk( im1_jp1_kp1_ , ijk )
        ipjpkp   = myijk( ip1_jp1_kp1_ , ijk )
        ipjpkm   = myijk( ip1_jp1_km1_ , ijk )
        ipjmkm   = myijk( ip1_jm1_km1_ , ijk )
        imjpkm   = myijk( im1_jp1_km1_ , ijk )
        imjmkm   = myijk( im1_jm1_km1_ , ijk )
        ijkpp    = myijk( ip0_jp0_kp2_ , ijk )
        ippjk    = myijk( ip2_jp0_kp0_ , ijk )
        immjk    = myijk( im2_jp0_kp0_ , ijk )
        ijkmm    = myijk( ip0_jp0_km2_ , ijk )
        ijppk    = myijk( ip0_jp2_kp0_ , ijk )
        ijmmk    = myijk( ip0_jm2_kp0_ , ijk )
        ijmmkp   = myijk( ip0_jm2_kp1_ , ijk )
        ippjpk   = myijk( ip2_jp1_kp0_ , ijk )
        ippjmk   = myijk( ip2_jm1_kp0_ , ijk )
        ipjmmk   = myijk( ip1_jm2_kp0_ , ijk )
        ippjkp   = myijk( ip2_jp0_kp1_ , ijk )
        ippjkm   = myijk( ip2_jp0_km1_ , ijk )
        ipjkpp   = myijk( ip1_jp0_kp2_ , ijk )
        imjkpp   = myijk( im1_jp0_kp2_ , ijk )
        ipjppk   = myijk( ip1_jp2_kp0_ , ijk )
        ijppkp   = myijk( ip0_jp2_kp1_ , ijk )
        ijpkmm   = myijk( ip0_jp1_km2_ , ijk )
        ipjkmm   = myijk( ip1_jp0_km2_ , ijk )
        ijmkpp   = myijk( ip0_jm1_kp2_ , ijk )
        ijpkpp   = myijk( ip0_jp1_kp2_ , ijk )
        ijppkm   = myijk( ip0_jp2_km1_ , ijk )
        immjpk   = myijk( im2_jp1_kp0_ , ijk )
        immjkp   = myijk( im2_jp0_kp1_ , ijk )
        imjppk   = myijk( im1_jp2_kp0_ , ijk )
        immjmk   = myijk( im2_jm1_kp0_ , ijk )
        immjkm   = myijk( im2_jp0_km1_ , ijk )
        imjmmk   = myijk( im1_jm2_kp0_ , ijk )
        ijmmkm   = myijk( ip0_jm2_km1_ , ijk )
        imjkmm   = myijk( im1_jp0_km2_ , ijk )
        ijmkmm   = myijk( ip0_jm1_km2_ , ijk )
        ipjmmkp  = myijk( ip1_jm2_kp1_ , ijk )
        ippjpkm  = myijk( ip2_jp1_km1_ , ijk )
        ippjmkm  = myijk( ip2_jm1_km1_ , ijk )
        ippjmkp  = myijk( ip2_jm1_kp1_ , ijk )
        imjppkp  = myijk( im1_jp2_kp1_ , ijk )
        ipjppkm  = myijk( ip1_jp2_km1_ , ijk )
        ipjpkmm  = myijk( ip1_jp1_km2_ , ijk )
        imjppkm  = myijk( im1_jp2_km1_ , ijk )
        imjmkpp  = myijk( im1_jm1_kp2_ , ijk )
        imjpkpp  = myijk( im1_jp1_kp2_ , ijk )
        ipjmkpp  = myijk( ip1_jm1_kp2_ , ijk )
        immjpkp  = myijk( im2_jp1_kp1_ , ijk )
        ippjpkp  = myijk( ip2_jp1_kp1_ , ijk )
        ipjppkp  = myijk( ip1_jp2_kp1_ , ijk )
        ipjpkpp  = myijk( ip1_jp1_kp2_ , ijk )
        imjpkmm  = myijk( im1_jp1_km2_ , ijk )
        ipjmkmm  = myijk( ip1_jm1_km2_ , ijk )
        imjmkmm  = myijk( im1_jm1_km2_ , ijk )
        ipjmmkm  = myijk( ip1_jm2_km1_ , ijk )
        imjmmkp  = myijk( im1_jm2_kp1_ , ijk )
        imjmmkm  = myijk( im1_jm2_km1_ , ijk )
        immjmkp  = myijk( im2_jm1_kp1_ , ijk )
        immjpkm  = myijk( im2_jp1_km1_ , ijk )
        immjmkm  = myijk( im2_jm1_km1_ , ijk )
        ippjkmm  = myijk( ip2_jp0_km2_ , ijk )
        immjppk  = myijk( im2_jp2_kp0_ , ijk )
        ijmmkpp  = myijk( ip0_jm2_kp2_ , ijk )
        ippjpkmm = myijk( ip2_jp1_km2_ , ijk )
        ippjmkmm = myijk( ip2_jm1_km2_ , ijk )
        immjppkp = myijk( im2_jp2_kp1_ , ijk )
        immjppkm = myijk( im2_jp2_km1_ , ijk )
        ipjmmkpp = myijk( ip1_jm2_kp2_ , ijk )
        imjmmkpp = myijk( im1_jm2_kp2_ , ijk )

        ijke     = myinds( ip1_jp0_kp0_ , ijk )
        ijkw     = myinds( im1_jp0_kp0_ , ijk )
        ijkn     = myinds( ip0_jp1_kp0_ , ijk )
        ijks     = myinds( ip0_jm1_kp0_ , ijk )
        ijkt     = myinds( ip0_jp0_kp1_ , ijk )
        ijkb     = myinds( ip0_jp0_km1_ , ijk )
        ijken    = myinds( ip1_jp1_kp0_ , ijk )
        ijkwn    = myinds( im1_jp1_kp0_ , ijk )
        ijkes    = myinds( ip1_jm1_kp0_ , ijk )
        ijkws    = myinds( im1_jm1_kp0_ , ijk )
        ijket    = myinds( ip1_jp0_kp1_ , ijk )
        ijkwt    = myinds( im1_jp0_kp1_ , ijk )
        ijknt    = myinds( ip0_jp1_kp1_ , ijk )
        ijkst    = myinds( ip0_jm1_kp1_ , ijk )
        ijkeb    = myinds( ip1_jp0_km1_ , ijk )
        ijkwb    = myinds( im1_jp0_km1_ , ijk )
        ijknb    = myinds( ip0_jp1_km1_ , ijk )
        ijksb    = myinds( ip0_jm1_km1_ , ijk )
        ijkwst   = myinds( im1_jm1_kp1_ , ijk )
        ijkest   = myinds( ip1_jm1_kp1_ , ijk )
        ijkwnt   = myinds( im1_jp1_kp1_ , ijk )
        ijkent   = myinds( ip1_jp1_kp1_ , ijk )
        ijkenb   = myinds( ip1_jp1_km1_ , ijk )
        ijkesb   = myinds( ip1_jm1_km1_ , ijk )
        ijkwnb   = myinds( im1_jp1_km1_ , ijk )
        ijkwsb   = myinds( im1_jm1_km1_ , ijk )
        ijkee    = myinds( ip2_jp0_kp0_ , ijk )
        ijktt    = myinds( ip0_jp0_kp2_ , ijk )
        ijkww    = myinds( im2_jp0_kp0_ , ijk )
        ijkbb    = myinds( ip0_jp0_km2_ , ijk )
        ijknn    = myinds( ip0_jp2_kp0_ , ijk )
        ijkss    = myinds( ip0_jm2_kp0_ , ijk )
        ijksst   = myinds( ip0_jm2_kp1_ , ijk )
        ijkeen   = myinds( ip2_jp1_kp0_ , ijk )
        ijkees   = myinds( ip2_jm1_kp0_ , ijk )
        ijkess   = myinds( ip1_jm2_kp0_ , ijk )
        ijkeet   = myinds( ip2_jp0_kp1_ , ijk ) 
        ijkeeb   = myinds( ip2_jp0_km1_ , ijk )
        ijkett   = myinds( ip1_jp0_kp2_ , ijk )
        ijkwtt   = myinds( im1_jp0_kp2_ , ijk )
        ijkenn   = myinds( ip1_jp2_kp0_ , ijk )
        ijknnt   = myinds( ip0_jp2_kp1_ , ijk )
        ijknbb   = myinds( ip0_jp1_km2_ , ijk )
        ijkebb   = myinds( ip1_jp0_km2_ , ijk )
        ijkstt   = myinds( ip0_jm1_kp2_ , ijk )
        ijkntt   = myinds( ip0_jp1_kp2_ , ijk )
        ijknnb   = myinds( ip0_jp2_km1_ , ijk )
        ijkwwn   = myinds( im2_jp1_kp0_ , ijk )
        ijkwwt   = myinds( im2_jp0_kp1_ , ijk )
        ijkwnn   = myinds( im1_jp2_kp0_ , ijk )
        ijkwws   = myinds( im2_jm1_kp0_ , ijk )
        ijkwwb   = myinds( im2_jp0_km1_ , ijk )
        ijkwss   = myinds( im1_jm2_kp0_ , ijk )
        ijkssb   = myinds( ip0_jm2_km1_ , ijk )
        ijkwbb   = myinds( im1_jp0_km2_ , ijk )
        ijksbb   = myinds( ip0_jm1_km2_ , ijk )
        ijkeent  = myinds( ip2_jp1_kp1_ , ijk )
        ijkeenb  = myinds( ip2_jp1_km1_ , ijk )
        ijkeesb  = myinds( ip2_jm1_km1_ , ijk )
        ijkeest  = myinds( ip2_jm1_kp1_ , ijk )
        ijkwwnt  = myinds( im2_jp1_kp1_ , ijk )
        ijkwwst  = myinds( im2_jm1_kp1_ , ijk )
        ijkwwnb  = myinds( im2_jp1_km1_ , ijk )
        ijkwwsb  = myinds( im2_jm1_km1_ , ijk )
        ijkwnnt  = myinds( im1_jp2_kp1_ , ijk )
        ijkennb  = myinds( ip1_jp2_km1_ , ijk )
        ijkennt  = myinds( ip1_jp2_kp1_ , ijk )
        ijkwnnb  = myinds( im1_jp2_km1_ , ijk )
        ijkesst  = myinds( ip1_jm2_kp1_ , ijk )
        ijkessb  = myinds( ip1_jm2_km1_ , ijk )
        ijkwsst  = myinds( im1_jm2_kp1_ , ijk )
        ijkwssb  = myinds( im1_jm2_km1_ , ijk )
        ijkwstt  = myinds( im1_jm1_kp2_ , ijk )
        ijkwntt  = myinds( im1_jp1_kp2_ , ijk )
        ijkestt  = myinds( ip1_jm1_kp2_ , ijk )
        ijkentt  = myinds( ip1_jp1_kp2_ , ijk )
        ijkenbb  = myinds( ip1_jp1_km2_ , ijk )
        ijkwnbb  = myinds( im1_jp1_km2_ , ijk )
        ijkesbb  = myinds( ip1_jm1_km2_ , ijk )
        ijkwsbb  = myinds( im1_jm1_km2_ , ijk )
        ijkeebb  = myinds( ip2_jp0_km2_ , ijk )
        ijkwwnn  = myinds( im2_jp2_kp0_ , ijk )
        ijksstt  = myinds( ip0_jm2_kp2_ , ijk )
        ijkeenbb = myinds( ip2_jp1_km2_ , ijk )
        ijkeesbb = myinds( ip2_jm1_km2_ , ijk )
        ijkwwnnt = myinds( im2_jp2_kp1_ , ijk )
        ijkwwnnb = myinds( im2_jp2_km1_ , ijk )
        ijkesstt = myinds( ip1_jm2_kp2_ , ijk )
        ijkwsstt = myinds( im1_jm2_kp2_ , ijk )

      END IF

      END SUBROUTINE ctu1_subscr
!-----------------------------------------------------------------------
      SUBROUTINE ctu2_subscr( ijk )
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ijk

      IF( job_type_flag == 2 ) THEN

        ipppjk = myijk( ip3_jp0_kp0_, ijk )
        ijkppp = myijk( ip0_jp0_kp3_, ijk )

        ijkeee = myinds( ip3_jp0_kp0_, ijk )
        ijkttt = myinds( ip0_jp0_kp3_, ijk )

      ELSE

        ipppjk = myijk( ip3_jp0_kp0_, ijk )
        ijpppk = myijk( ip0_jp3_kp0_, ijk )
        ijkppp = myijk( ip0_jp0_kp3_, ijk )

        ijkeee = myinds( ip3_jp0_kp0_, ijk )
        ijknnn = myinds( ip0_jp3_kp0_, ijk )
        ijkttt = myinds( ip0_jp0_kp3_, ijk )

      END IF

      END SUBROUTINE ctu2_subscr
!-----------------------------------------------------------------------
      SUBROUTINE ctu3_subscr( ijk )
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ijk

      IF( job_type_flag == 2 ) THEN

        ipppjkp = myijk( ip3_jp0_kp1_, ijk )
        ipjkppp = myijk( ip1_jp0_kp3_, ijk )
        ippjkpp = myijk( ip2_jp0_kp2_, ijk )
        ippjkmm = myijk( ip2_jp0_km2_, ijk )
        immjkpp = myijk( im2_jp0_kp2_, ijk )

        ijkeeet = myinds( ip3_jp0_kp1_, ijk )
        ijkettt = myinds( ip1_jp0_kp3_, ijk )
        ijkeett = myinds( ip2_jp0_kp2_, ijk )
        ijkeebb = myinds( ip2_jp0_km2_, ijk )
        ijkwwtt = myinds( im2_jp0_kp2_, ijk )

      ELSE

        immjkm   = myijk( im2_jp0_km1_ , ijk )
        imjkmm   = myijk( im1_jp0_km2_ , ijk )
        ijmmkm   = myijk( ip0_jm2_km1_ , ijk )
        immjmk   = myijk( im2_jm1_kp0_ , ijk )
        ijmkmm   = myijk( ip0_jm1_km2_ , ijk )
        imjmmk   = myijk( im1_jm2_kp0_ , ijk )
        imjpkmm  = myijk( im1_jp1_km2_ , ijk )
        ipjmmkm  = myijk( ip1_jm2_km1_ , ijk )
        immjmkp  = myijk( im2_jm1_kp1_ , ijk )
        ipjmkmm  = myijk( ip1_jm1_km2_ , ijk )
        imjmkmm  = myijk( im1_jm1_km2_ , ijk )
        imjmmkm  = myijk( im1_jm2_km1_ , ijk )
        imjmmkp  = myijk( im1_jm2_kp1_ , ijk )
        immjpkm  = myijk( im2_jp1_km1_ , ijk )
        immjmkm  = myijk( im2_jm1_km1_ , ijk )
        ijpkppp  = myijk( ip0_jp1_kp3_ , ijk )
        ipjkppp  = myijk( ip1_jp0_kp3_ , ijk )
        imjkppp  = myijk( im1_jp0_kp3_ , ijk )
        ijpppkm  = myijk( ip0_jp3_km1_ , ijk )
        ijpppkp  = myijk( ip0_jp3_kp1_ , ijk )
        ipjpppk  = myijk( ip1_jp3_kp0_ , ijk )
        ipppjkp  = myijk( ip3_jp0_kp1_ , ijk )
        ipppjmk  = myijk( ip3_jm1_kp0_ , ijk )
        ipppjpk  = myijk( ip3_jp1_kp0_ , ijk )
        ijmmkpp  = myijk( ip0_jm2_kp2_ , ijk )
        ijppkpp  = myijk( ip0_jp2_kp2_ , ijk )
        immjkpp  = myijk( im2_jp0_kp2_ , ijk )
        ippjkpp  = myijk( ip2_jp0_kp2_ , ijk )
        ijppkmm  = myijk( ip0_jp2_km2_ , ijk )
        immjppk  = myijk( im2_jp2_kp0_ , ijk )
        ippjppk  = myijk( ip2_jp2_kp0_ , ijk )
        ippjkmm  = myijk( ip2_jp0_km2_ , ijk )
        ippjmmk  = myijk( ip2_jm2_kp0_ , ijk )
        ipjpkppp = myijk( ip1_jp1_kp3_ , ijk )
        imjpkppp = myijk( im1_jp1_kp3_ , ijk )
        ipjmmkpp = myijk( ip1_jm2_kp2_ , ijk )
        ipjppkpp = myijk( ip1_jp2_kp2_ , ijk )
        immjmkpp = myijk( im2_jm1_kp2_ , ijk )
        immjpkpp = myijk( im2_jp1_kp2_ , ijk )
        ippjmkpp = myijk( ip2_jm1_kp2_ , ijk )
        ippjpkpp = myijk( ip2_jp1_kp2_ , ijk )
        ipjpppkm = myijk( ip1_jp3_km1_ , ijk )
        ipjpppkp = myijk( ip1_jp3_kp1_ , ijk )
        ipjppkmm = myijk( ip1_jp2_km2_ , ijk )
        imjppkmm = myijk( im1_jp2_km2_ , ijk )
        imjppkpp = myijk( im1_jp2_kp2_ , ijk )
        immjppkp = myijk( im2_jp2_kp1_ , ijk )
        ippjppkp = myijk( ip2_jp2_kp1_ , ijk )
        ippjpkmm = myijk( ip2_jp1_km2_ , ijk )
        ippjmmkm = myijk( ip2_jm2_km1_ , ijk )
        ippjmmkp = myijk( ip2_jm2_kp1_ , ijk )
        ippjppkm = myijk( ip2_jp2_km1_ , ijk )
        ipppjmkp = myijk( ip3_jm1_kp1_ , ijk )
        ipppjpkp = myijk( ip3_jp1_kp1_ , ijk )

        ijkwwb   = myinds( im2_jp0_km1_ , ijk )
        ijkwbb   = myinds( im1_jp0_km2_ , ijk )
        ijkssb   = myinds( ip0_jm2_km1_ , ijk )
        ijkwws   = myinds( im2_jm1_kp0_ , ijk )
        ijksbb   = myinds( ip0_jm1_km2_ , ijk )
        ijkwss   = myinds( im1_jm2_kp0_ , ijk )
        ijkwnbb  = myinds( im1_jp1_km2_ , ijk )
        ijkessb  = myinds( ip1_jm2_km1_ , ijk )
        ijkwwst  = myinds( im2_jm1_kp1_ , ijk )
        ijkesbb  = myinds( ip1_jm1_km2_ , ijk )
        ijkwsbb  = myinds( im1_jm1_km2_ , ijk )
        ijkwssb  = myinds( im1_jm2_km1_ , ijk )
        ijkwsst  = myinds( im1_jm2_kp1_ , ijk )
        ijkwwnb  = myinds( im2_jp1_km1_ , ijk )
        ijkwwsb  = myinds( im2_jm1_km1_ , ijk )
        ijkenttt = myinds( ip1_jp1_kp3_ , ijk )
        ijkwnttt = myinds( im1_jp1_kp3_ , ijk )
        ijknttt  = myinds( ip0_jp1_kp3_ , ijk )
        ijkettt  = myinds( ip1_jp0_kp3_ , ijk )
        ijkwttt  = myinds( im1_jp0_kp3_ , ijk )
        ijkesstt = myinds( ip1_jm2_kp2_ , ijk )
        ijksstt  = myinds( ip0_jm2_kp2_ , ijk )
        ijkenntt = myinds( ip1_jp2_kp2_ , ijk )
        ijknntt  = myinds( ip0_jp2_kp2_ , ijk )
        ijkwwstt = myinds( im2_jm1_kp2_ , ijk )
        ijkwwntt = myinds( im2_jp1_kp2_ , ijk )
        ijkwwtt  = myinds( im2_jp0_kp2_ , ijk )
        ijkeestt = myinds( ip2_jm1_kp2_ , ijk )
        ijkeentt = myinds( ip2_jp1_kp2_ , ijk )
        ijkeett  = myinds( ip2_jp0_kp2_ , ijk ) 
        ijkennnb = myinds( ip1_jp3_km1_ , ijk ) 
        ijkennnt = myinds( ip1_jp3_kp1_ , ijk ) 
        ijknnnb  = myinds( ip0_jp3_km1_ , ijk ) 
        ijknnnt  = myinds( ip0_jp3_kp1_ , ijk ) 
        ijkennn  = myinds( ip1_jp3_kp0_ , ijk ) 
        ijkennbb = myinds( ip1_jp2_km2_ , ijk ) 
        ijkwnnbb = myinds( im1_jp2_km2_ , ijk ) 
        ijknnbb  = myinds( ip0_jp2_km2_ , ijk ) 
        ijkwnntt = myinds( im1_jp2_kp2_ , ijk ) 
        ijkwwnnt = myinds( im2_jp2_kp1_ , ijk ) 
        ijkwwnn  = myinds( im2_jp2_kp0_ , ijk ) 
        ijkeennt = myinds( ip2_jp2_kp1_ , ijk ) 
        ijkeenn  = myinds( ip2_jp2_kp0_ , ijk ) 
        ijkeenbb = myinds( ip2_jp1_km2_ , ijk ) 
        ijkeebb  = myinds( ip2_jp0_km2_ , ijk ) 
        ijkeessb = myinds( ip2_jm2_km1_ , ijk ) 
        ijkeesst = myinds( ip2_jm2_kp1_ , ijk ) 
        ijkeess  = myinds( ip2_jm2_kp0_ , ijk ) 
        ijkeennb = myinds( ip2_jp2_km1_ , ijk ) 
        ijkeeest = myinds( ip3_jm1_kp1_ , ijk ) 
        ijkeeent = myinds( ip3_jp1_kp1_ , ijk ) 
        ijkeeet  = myinds( ip3_jp0_kp1_ , ijk ) 
        ijkeees  = myinds( ip3_jm1_kp0_ , ijk ) 
        ijkeeen  = myinds( ip3_jp1_kp0_ , ijk ) 

      END IF

      END SUBROUTINE ctu3_subscr
!-----------------------------------------------------------------------
      FUNCTION cte(c)
! ... build a stencil with constant values c at each location
!
      IMPLICIT NONE
        REAL*8, INTENT(IN) :: c
        TYPE(stencil) :: cte

        cte%c     = c
        cte%e     = c
        cte%w     = c
        cte%ee    = c
        cte%ww    = c
        cte%n     = c
        cte%en    = c
        cte%wn    = c
        cte%s     = c
        cte%es    = c
        cte%ws    = c
        cte%nn    = c
        cte%ss    = c
        cte%t     = c
        cte%et    = c
        cte%wt    = c
        cte%nt    = c
        cte%st    = c
        cte%b     = c
        cte%eb    = c
        cte%wb    = c
        cte%nb    = c
        cte%sb    = c
        cte%tt    = c
        cte%bb    = c
        cte%eet   = c
        cte%eeb   = c
        cte%wtt   = c
        cte%ett   = c
        cte%wst   = c
        cte%est   = c
        cte%wnt   = c
        cte%ent   = c
        cte%enb   = c
        cte%esb   = c
        cte%wnb   = c
        cte%sst   = c
        cte%een   = c
        cte%ees   = c
        cte%ess   = c
        cte%esst  = c
        cte%eenb  = c
        cte%eesb  = c
        cte%eest  = c
        cte%wnnt  = c
        cte%nnt   = c
        cte%ennb  = c
        cte%enn   = c
        cte%ebb   = c
        cte%enbb  = c
        cte%wnnb  = c
        cte%wnn   = c
        cte%nbb   = c
        cte%nnb   = c
        cte%wstt  = c
        cte%stt   = c
        cte%wntt  = c
        cte%ntt   = c
        cte%estt  = c
        cte%wwn   = c
        cte%wwnt  = c
        cte%wwt   = c
        cte%wsb   = c
        cte%eent  = c
        cte%ennt  = c
        cte%entt  = c
        cte%wbb   = c
        cte%ssb   = c
        cte%wws   = c
        cte%sbb   = c
        cte%wss   = c
        cte%wwb   = c
        cte%wnbb  = c
        cte%essb  = c
        cte%wwst  = c
        cte%esbb  = c
        cte%wsbb  = c
        cte%wssb  = c
        cte%wsst  = c
        cte%wwnb  = c
        cte%wwsb  = c
        cte%entt  = c
        cte%wntt  = c
        cte%nttt  = c
        cte%ettt  = c
        cte%wttt  = c
        cte%ttt   = c
        cte%esstt = c
        cte%sstt  = c
        cte%enntt = c
        cte%nntt  = c
        cte%wwstt = c
        cte%wwntt = c
        cte%wwtt  = c
        cte%eestt = c
        cte%eentt = c
        cte%eett  = c
        cte%ennbb = c
        cte%wnnbb = c
        cte%nnbb  = c
        cte%wnntt = c
        cte%ennnb = c
        cte%ennnt = c
        cte%nnnb  = c
        cte%nnnt  = c
        cte%ennn  = c
        cte%wwnnt = c
        cte%wwnn  = c
        cte%eennt = c
        cte%eenn  = c
        cte%nnn   = c
        cte%eenbb = c
        cte%eebb  = c
        cte%eessb = c
        cte%eesst = c
        cte%eess  = c
        cte%eennb = c
        cte%eeest = c
        cte%eeent = c
        cte%eeet  = c
        cte%eees  = c
        cte%eeen  = c
        cte%eee   = c

      END FUNCTION cte
!-----------------------------------------------------------------------
      FUNCTION maskval(c)
! ... build a stencil with constant values c at each location
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: c
      TYPE(masks) :: maskval

        maskval%c  = c
        maskval%e  = c
        maskval%w  = c
        maskval%ee = c
        maskval%ww = c
        maskval%n  = c
        maskval%en = c
        maskval%wn = c
        maskval%s  = c
        maskval%es = c
        maskval%ws = c
        maskval%nn = c
        maskval%ss = c
        maskval%t  = c
        maskval%et = c
        maskval%wt = c
        maskval%nt = c
        maskval%st = c
        maskval%b  = c
        maskval%eb = c
        maskval%wb = c
        maskval%nb = c
        maskval%sb = c
        maskval%tt = c
        maskval%bb = c

      END FUNCTION maskval
!-----------------------------------------------------------------------
      FUNCTION sumstencil(st1, st2)
! ... sum two stencils
!
      IMPLICIT NONE
        TYPE(stencil), INTENT(in) :: st1, st2
        TYPE(stencil) :: sumstencil

        sumstencil%c     = st1%c     + st2%c
        sumstencil%e     = st1%e     + st2%e
        sumstencil%w     = st1%w     + st2%w
        sumstencil%ee    = st1%ee    + st2%ee
        sumstencil%ww    = st1%ww    + st2%ww
        sumstencil%n     = st1%n     + st2%n
        sumstencil%en    = st1%en    + st2%en
        sumstencil%wn    = st1%wn    + st2%wn
        sumstencil%s     = st1%s     + st2%s
        sumstencil%es    = st1%es    + st2%es
        sumstencil%ws    = st1%ws    + st2%ws
        sumstencil%nn    = st1%nn    + st2%nn
        sumstencil%ss    = st1%ss    + st2%ss
        sumstencil%t     = st1%t     + st2%t
        sumstencil%et    = st1%et    + st2%et
        sumstencil%wt    = st1%wt    + st2%wt
        sumstencil%nt    = st1%nt    + st2%nt
        sumstencil%st    = st1%st    + st2%st
        sumstencil%b     = st1%b     + st2%b
        sumstencil%eb    = st1%eb    + st2%eb
        sumstencil%wb    = st1%wb    + st2%wb
        sumstencil%nb    = st1%nb    + st2%nb
        sumstencil%sb    = st1%sb    + st2%sb
        sumstencil%tt    = st1%tt    + st2%tt
        sumstencil%bb    = st1%bb    + st2%bb
        sumstencil%eet   = st1%eet   + st2%eet
        sumstencil%eeb   = st1%eeb   + st2%eeb
        sumstencil%ett   = st1%ett   + st2%ett
        sumstencil%wtt   = st1%wtt   + st2%wtt
        sumstencil%wst   = st1%wst   + st2%wst 
        sumstencil%est   = st1%est   + st2%est 
        sumstencil%wnt   = st1%wnt   + st2%wnt 
        sumstencil%ent   = st1%ent   + st2%ent 
        sumstencil%enb   = st1%enb   + st2%enb 
        sumstencil%esb   = st1%esb   + st2%esb 
        sumstencil%wnb   = st1%wnb   + st2%wnb 
        sumstencil%sst   = st1%sst   + st2%sst 
        sumstencil%een   = st1%een   + st2%een 
        sumstencil%ees   = st1%ees   + st2%ees 
        sumstencil%ess   = st1%ess   + st2%ess 
        sumstencil%esst  = st1%esst  + st2%esst 
        sumstencil%eenb  = st1%eenb  + st2%eenb 
        sumstencil%eesb  = st1%eesb  + st2%eesb 
        sumstencil%eest  = st1%eest  + st2%eest 
        sumstencil%wnnt  = st1%wnnt  + st2%wnnt 
        sumstencil%nnt   = st1%nnt   + st2%nnt 
        sumstencil%ennb  = st1%ennb  + st2%ennb 
        sumstencil%enn   = st1%enn   + st2%enn 
        sumstencil%ebb   = st1%ebb   + st2%ebb 
        sumstencil%enbb  = st1%enbb  + st2%enbb 
        sumstencil%wnnb  = st1%wnnb  + st2%wnnb 
        sumstencil%wnn   = st1%wnn   + st2%wnn 
        sumstencil%nbb   = st1%nbb   + st2%nbb 
        sumstencil%nnb   = st1%nnb   + st2%nnb 
        sumstencil%wstt  = st1%wstt  + st2%wstt 
        sumstencil%stt   = st1%stt   + st2%stt 
        sumstencil%wntt  = st1%wntt  + st2%wntt 
        sumstencil%ntt   = st1%ntt   + st2%ntt 
        sumstencil%estt  = st1%estt  + st2%estt 
        sumstencil%wwn   = st1%wwn   + st2%wwn 
        sumstencil%wwnt  = st1%wwnt  + st2%wwnt 
        sumstencil%wwt   = st1%wwt   + st2%wwt 
        sumstencil%wsb   = st1%wsb   + st2%wsb 
        sumstencil%eent  = st1%eent  + st2%eent 
        sumstencil%ennt  = st1%ennt  + st2%ennt 
        sumstencil%entt  = st1%entt  + st2%entt 
        sumstencil%wbb   = st1%wbb   + st2%wbb  
        sumstencil%ssb   = st1%ssb   + st2%ssb  
        sumstencil%wws   = st1%wws   + st2%wws  
        sumstencil%sbb   = st1%sbb   + st2%sbb  
        sumstencil%wss   = st1%wss   + st2%wss  
        sumstencil%wwb   = st1%wwb   + st2%wwb  
        sumstencil%wnbb  = st1%wnbb  + st2%wnbb  
        sumstencil%essb  = st1%essb  + st2%essb  
        sumstencil%wwst  = st1%wwst  + st2%wwst  
        sumstencil%esbb  = st1%esbb  + st2%esbb  
        sumstencil%wsbb  = st1%wsbb  + st2%wsbb  
        sumstencil%wssb  = st1%wssb  + st2%wssb  
        sumstencil%wsst  = st1%wsst  + st2%wsst  
        sumstencil%wwnb  = st1%wwnb  + st2%wwnb  
        sumstencil%wwsb  = st1%wwsb  + st2%wwsb  
        sumstencil%enttt = st1%enttt + st2%enttt  
        sumstencil%wnttt = st1%wnttt + st2%wnttt  
        sumstencil%nttt  = st1%nttt  + st2%nttt  
        sumstencil%ettt  = st1%ettt  + st2%ettt  
        sumstencil%wttt  = st1%wttt  + st2%wttt  
        sumstencil%ttt   = st1%ttt   + st2%ttt  
        sumstencil%esstt = st1%esstt + st2%esstt  
        sumstencil%sstt  = st1%sstt  + st2%sstt  
        sumstencil%enntt = st1%enntt + st2%enntt  
        sumstencil%nntt  = st1%nntt  + st2%nntt  
        sumstencil%wwstt = st1%wwstt + st2%wwstt  
        sumstencil%wwntt = st1%wwntt + st2%wwntt  
        sumstencil%wwtt  = st1%wwtt  + st2%wwtt  
        sumstencil%eestt = st1%eestt + st2%eestt  
        sumstencil%eentt = st1%eentt + st2%eentt  
        sumstencil%eett  = st1%eett  + st2%eett  
        sumstencil%ennbb = st1%ennbb + st2%ennbb  
        sumstencil%wnnbb = st1%wnnbb + st2%wnnbb  
        sumstencil%nnbb  = st1%nnbb  + st2%nnbb  
        sumstencil%wnntt = st1%wnntt + st2%wnntt  
        sumstencil%ennnb = st1%ennnb + st2%ennnb  
        sumstencil%ennnt = st1%ennnt + st2%ennnt  
        sumstencil%nnnb  = st1%nnnb  + st2%nnnb  
        sumstencil%nnnt  = st1%nnnt  + st2%nnnt  
        sumstencil%ennn  = st1%ennn  + st2%ennn  
        sumstencil%wwnnt = st1%wwnnt + st2%wwnnt  
        sumstencil%wwnn  = st1%wwnn  + st2%wwnn  
        sumstencil%eennt = st1%eennt + st2%eennt  
        sumstencil%eenn  = st1%eenn  + st2%eenn  
        sumstencil%nnn   = st1%nnn   + st2%nnn  
        sumstencil%eenbb = st1%eenbb + st2%eenbb  
        sumstencil%eebb  = st1%eebb  + st2%eebb  
        sumstencil%eessb = st1%eessb + st2%eessb  
        sumstencil%eesst = st1%eesst + st2%eesst 
        sumstencil%eess  = st1%eess  + st2%eess  
        sumstencil%eennb = st1%eennb + st2%eennb  
        sumstencil%eeest = st1%eeest + st2%eeest  
        sumstencil%eeent = st1%eeent + st2%eeent  
        sumstencil%eeet  = st1%eeet  + st2%eeet  
        sumstencil%eees  = st1%eees  + st2%eees  
        sumstencil%eeen  = st1%eeen  + st2%eeen  
        sumstencil%eee   = st1%eee   + st2%eee  

      END FUNCTION sumstencil
!-----------------------------------------------------------------------
      FUNCTION difstencil(st1, st2)
! ... subtract two stencils
!
      IMPLICIT NONE
        TYPE(stencil), INTENT(in) :: st1, st2
        TYPE(stencil) :: difstencil

        difstencil%c     = st1%c     - st2%c
        difstencil%e     = st1%e     - st2%e
        difstencil%w     = st1%w     - st2%w
        difstencil%ee    = st1%ee    - st2%ee
        difstencil%ww    = st1%ww    - st2%ww
        difstencil%n     = st1%n     - st2%n
        difstencil%en    = st1%en    - st2%en
        difstencil%wn    = st1%wn    - st2%wn
        difstencil%s     = st1%s     - st2%s
        difstencil%es    = st1%es    - st2%es
        difstencil%ws    = st1%ws    - st2%ws
        difstencil%nn    = st1%nn    - st2%nn
        difstencil%ss    = st1%ss    - st2%ss
        difstencil%t     = st1%t     - st2%t
        difstencil%et    = st1%et    - st2%et
        difstencil%wt    = st1%wt    - st2%wt
        difstencil%nt    = st1%nt    - st2%nt
        difstencil%st    = st1%st    - st2%st
        difstencil%b     = st1%b     - st2%b
        difstencil%eb    = st1%eb    - st2%eb
        difstencil%wb    = st1%wb    - st2%wb
        difstencil%nb    = st1%nb    - st2%nb
        difstencil%sb    = st1%sb    - st2%sb
        difstencil%tt    = st1%tt    - st2%tt
        difstencil%bb    = st1%bb    - st2%bb
        difstencil%eet   = st1%eet   - st2%eet
        difstencil%eeb   = st1%eeb   - st2%eeb
        difstencil%ett   = st1%ett   - st2%ett
        difstencil%wtt   = st1%wtt   - st2%wtt
        difstencil%wst   = st1%wst   - st2%wst 
        difstencil%est   = st1%est   - st2%est 
        difstencil%wnt   = st1%wnt   - st2%wnt 
        difstencil%ent   = st1%ent   - st2%ent 
        difstencil%enb   = st1%enb   - st2%enb 
        difstencil%esb   = st1%esb   - st2%esb 
        difstencil%wnb   = st1%wnb   - st2%wnb 
        difstencil%sst   = st1%sst   - st2%sst 
        difstencil%een   = st1%een   - st2%een 
        difstencil%ees   = st1%ees   - st2%ees 
        difstencil%ess   = st1%ess   - st2%ess 
        difstencil%esst  = st1%esst  - st2%esst 
        difstencil%eenb  = st1%eenb  - st2%eenb 
        difstencil%eesb  = st1%eesb  - st2%eesb 
        difstencil%eest  = st1%eest  - st2%eest 
        difstencil%wnnt  = st1%wnnt  - st2%wnnt 
        difstencil%nnt   = st1%nnt   - st2%nnt 
        difstencil%ennb  = st1%ennb  - st2%ennb 
        difstencil%enn   = st1%enn   - st2%enn 
        difstencil%ebb   = st1%ebb   - st2%ebb 
        difstencil%enbb  = st1%enbb  - st2%enbb 
        difstencil%wnnb  = st1%wnnb  - st2%wnnb 
        difstencil%wnn   = st1%wnn   - st2%wnn 
        difstencil%nbb   = st1%nbb   - st2%nbb 
        difstencil%nnb   = st1%nnb   - st2%nnb 
        difstencil%wstt  = st1%wstt  - st2%wstt 
        difstencil%stt   = st1%stt   - st2%stt 
        difstencil%wntt  = st1%wntt  - st2%wntt 
        difstencil%ntt   = st1%ntt   - st2%ntt 
        difstencil%estt  = st1%estt  - st2%estt 
        difstencil%wwn   = st1%wwn   - st2%wwn 
        difstencil%wwnt  = st1%wwnt  - st2%wwnt 
        difstencil%wwt   = st1%wwt   - st2%wwt 
        difstencil%wsb   = st1%wsb   - st2%wsb 
        difstencil%eent  = st1%eent  - st2%eent 
        difstencil%ennt  = st1%ennt  - st2%ennt 
        difstencil%entt  = st1%entt  - st2%entt 
        difstencil%wbb   = st1%wbb   - st2%wbb  
        difstencil%ssb   = st1%ssb   - st2%ssb  
        difstencil%wws   = st1%wws   - st2%wws  
        difstencil%sbb   = st1%sbb   - st2%sbb  
        difstencil%wss   = st1%wss   - st2%wss  
        difstencil%wwb   = st1%wwb   - st2%wwb  
        difstencil%wnbb  = st1%wnbb  - st2%wnbb  
        difstencil%essb  = st1%essb  - st2%essb  
        difstencil%wwst  = st1%wwst  - st2%wwst  
        difstencil%esbb  = st1%esbb  - st2%esbb  
        difstencil%wsbb  = st1%wsbb  - st2%wsbb  
        difstencil%wssb  = st1%wssb  - st2%wssb  
        difstencil%wsst  = st1%wsst  - st2%wsst  
        difstencil%wwnb  = st1%wwnb  - st2%wwnb  
        difstencil%wwsb  = st1%wwsb  - st2%wwsb  
        difstencil%enttt = st1%enttt - st2%enttt  
        difstencil%wnttt = st1%wnttt - st2%wnttt  
        difstencil%nttt  = st1%nttt  - st2%nttt  
        difstencil%ettt  = st1%ettt  - st2%ettt  
        difstencil%wttt  = st1%wttt  - st2%wttt  
        difstencil%ttt   = st1%ttt   - st2%ttt  
        difstencil%esstt = st1%esstt - st2%esstt  
        difstencil%sstt  = st1%sstt  - st2%sstt  
        difstencil%enntt = st1%enntt - st2%enntt  
        difstencil%nntt  = st1%nntt  - st2%nntt  
        difstencil%wwstt = st1%wwstt - st2%wwstt  
        difstencil%wwntt = st1%wwntt - st2%wwntt  
        difstencil%wwtt  = st1%wwtt  - st2%wwtt  
        difstencil%eestt = st1%eestt - st2%eestt  
        difstencil%eentt = st1%eentt - st2%eentt  
        difstencil%eett  = st1%eett  - st2%eett  
        difstencil%ennbb = st1%ennbb - st2%ennbb  
        difstencil%wnnbb = st1%wnnbb - st2%wnnbb  
        difstencil%nnbb  = st1%nnbb  - st2%nnbb  
        difstencil%wnntt = st1%wnntt - st2%wnntt  
        difstencil%ennnb = st1%ennnb - st2%ennnb  
        difstencil%ennnt = st1%ennnt - st2%ennnt  
        difstencil%nnnb  = st1%nnnb  - st2%nnnb  
        difstencil%nnnt  = st1%nnnt  - st2%nnnt  
        difstencil%ennn  = st1%ennn  - st2%ennn  
        difstencil%wwnnt = st1%wwnnt - st2%wwnnt  
        difstencil%wwnn  = st1%wwnn  - st2%wwnn  
        difstencil%eennt = st1%eennt - st2%eennt  
        difstencil%eenn  = st1%eenn  - st2%eenn  
        difstencil%nnn   = st1%nnn   - st2%nnn  
        difstencil%eenbb = st1%eenbb - st2%eenbb  
        difstencil%eebb  = st1%eebb  - st2%eebb  
        difstencil%eessb = st1%eessb - st2%eessb  
        difstencil%eesst = st1%eesst - st2%eesst 
        difstencil%eess  = st1%eess  - st2%eess  
        difstencil%eennb = st1%eennb - st2%eennb  
        difstencil%eeest = st1%eeest - st2%eeest  
        difstencil%eeent = st1%eeent - st2%eeent  
        difstencil%eeet  = st1%eeet  - st2%eeet  
        difstencil%eees  = st1%eees  - st2%eees  
        difstencil%eeen  = st1%eeen  - st2%eeen  
        difstencil%eee   = st1%eee   - st2%eee  

      END FUNCTION difstencil
!-----------------------------------------------------------------------
      FUNCTION prodstencil(st1, st2)
! ... multiply two stencils, location by location
!
      IMPLICIT NONE
        TYPE(stencil), INTENT(in) :: st1, st2
        TYPE(stencil) :: prodstencil

        prodstencil%c     = st1%c     * st2%c
        prodstencil%e     = st1%e     * st2%e
        prodstencil%w     = st1%w     * st2%w
        prodstencil%ee    = st1%ee    * st2%ee
        prodstencil%ww    = st1%ww    * st2%ww
        prodstencil%n     = st1%n     * st2%n
        prodstencil%en    = st1%en    * st2%en
        prodstencil%wn    = st1%wn    * st2%wn
        prodstencil%s     = st1%s     * st2%s
        prodstencil%es    = st1%es    * st2%es
        prodstencil%ws    = st1%ws    * st2%ws
        prodstencil%nn    = st1%nn    * st2%nn
        prodstencil%ss    = st1%ss    * st2%ss
        prodstencil%t     = st1%t     * st2%t
        prodstencil%et    = st1%et    * st2%et
        prodstencil%wt    = st1%wt    * st2%wt
        prodstencil%nt    = st1%nt    * st2%nt
        prodstencil%st    = st1%st    * st2%st
        prodstencil%b     = st1%b     * st2%b
        prodstencil%eb    = st1%eb    * st2%eb
        prodstencil%wb    = st1%wb    * st2%wb
        prodstencil%nb    = st1%nb    * st2%nb
        prodstencil%sb    = st1%sb    * st2%sb
        prodstencil%tt    = st1%tt    * st2%tt
        prodstencil%bb    = st1%bb    * st2%bb
        prodstencil%eet   = st1%eet   * st2%eet
        prodstencil%eeb   = st1%eeb   * st2%eeb
        prodstencil%ett   = st1%ett   * st2%ett
        prodstencil%wtt   = st1%wtt   * st2%wtt
        prodstencil%wst   = st1%wst   * st2%wst 
        prodstencil%est   = st1%est   * st2%est 
        prodstencil%wnt   = st1%wnt   * st2%wnt 
        prodstencil%ent   = st1%ent   * st2%ent 
        prodstencil%enb   = st1%enb   * st2%enb 
        prodstencil%esb   = st1%esb   * st2%esb 
        prodstencil%wnb   = st1%wnb   * st2%wnb 
        prodstencil%sst   = st1%sst   * st2%sst 
        prodstencil%een   = st1%een   * st2%een 
        prodstencil%ees   = st1%ees   * st2%ees 
        prodstencil%ess   = st1%ess   * st2%ess 
        prodstencil%esst  = st1%esst  * st2%esst 
        prodstencil%eenb  = st1%eenb  * st2%eenb 
        prodstencil%eesb  = st1%eesb  * st2%eesb 
        prodstencil%eest  = st1%eest  * st2%eest 
        prodstencil%wnnt  = st1%wnnt  * st2%wnnt 
        prodstencil%nnt   = st1%nnt   * st2%nnt 
        prodstencil%ennb  = st1%ennb  * st2%ennb 
        prodstencil%enn   = st1%enn   * st2%enn 
        prodstencil%ebb   = st1%ebb   * st2%ebb 
        prodstencil%enbb  = st1%enbb  * st2%enbb 
        prodstencil%wnnb  = st1%wnnb  * st2%wnnb 
        prodstencil%wnn   = st1%wnn   * st2%wnn 
        prodstencil%nbb   = st1%nbb   * st2%nbb 
        prodstencil%nnb   = st1%nnb   * st2%nnb 
        prodstencil%wstt  = st1%wstt  * st2%wstt 
        prodstencil%stt   = st1%stt   * st2%stt 
        prodstencil%wntt  = st1%wntt  * st2%wntt 
        prodstencil%ntt   = st1%ntt   * st2%ntt 
        prodstencil%estt  = st1%estt  * st2%estt 
        prodstencil%wwn   = st1%wwn   * st2%wwn 
        prodstencil%wwnt  = st1%wwnt  * st2%wwnt 
        prodstencil%wwt   = st1%wwt   * st2%wwt 
        prodstencil%wsb   = st1%wsb   * st2%wsb 
        prodstencil%eent  = st1%eent  * st2%eent 
        prodstencil%ennt  = st1%ennt  * st2%ennt 
        prodstencil%entt  = st1%entt  * st2%entt 
        prodstencil%wbb   = st1%wbb   * st2%wbb  
        prodstencil%ssb   = st1%ssb   * st2%ssb  
        prodstencil%wws   = st1%wws   * st2%wws  
        prodstencil%sbb   = st1%sbb   * st2%sbb  
        prodstencil%wss   = st1%wss   * st2%wss  
        prodstencil%wwb   = st1%wwb   * st2%wwb  
        prodstencil%wnbb  = st1%wnbb  * st2%wnbb  
        prodstencil%essb  = st1%essb  * st2%essb  
        prodstencil%wwst  = st1%wwst  * st2%wwst  
        prodstencil%esbb  = st1%esbb  * st2%esbb  
        prodstencil%wsbb  = st1%wsbb  * st2%wsbb  
        prodstencil%wssb  = st1%wssb  * st2%wssb  
        prodstencil%wsst  = st1%wsst  * st2%wsst  
        prodstencil%wwnb  = st1%wwnb  * st2%wwnb  
        prodstencil%wwsb  = st1%wwsb  * st2%wwsb  
        prodstencil%enttt = st1%enttt * st2%enttt  
        prodstencil%wnttt = st1%wnttt * st2%wnttt  
        prodstencil%nttt  = st1%nttt  * st2%nttt  
        prodstencil%ettt  = st1%ettt  * st2%ettt  
        prodstencil%wttt  = st1%wttt  * st2%wttt  
        prodstencil%ttt   = st1%ttt   * st2%ttt  
        prodstencil%esstt = st1%esstt * st2%esstt  
        prodstencil%sstt  = st1%sstt  * st2%sstt  
        prodstencil%enntt = st1%enntt * st2%enntt  
        prodstencil%nntt  = st1%nntt  * st2%nntt  
        prodstencil%wwstt = st1%wwstt * st2%wwstt  
        prodstencil%wwntt = st1%wwntt * st2%wwntt  
        prodstencil%wwtt  = st1%wwtt  * st2%wwtt  
        prodstencil%eestt = st1%eestt * st2%eestt  
        prodstencil%eentt = st1%eentt * st2%eentt  
        prodstencil%eett  = st1%eett  * st2%eett  
        prodstencil%ennbb = st1%ennbb * st2%ennbb  
        prodstencil%wnnbb = st1%wnnbb * st2%wnnbb  
        prodstencil%nnbb  = st1%nnbb  * st2%nnbb  
        prodstencil%wnntt = st1%wnntt * st2%wnntt  
        prodstencil%ennnb = st1%ennnb * st2%ennnb  
        prodstencil%ennnt = st1%ennnt * st2%ennnt  
        prodstencil%nnnb  = st1%nnnb  * st2%nnnb  
        prodstencil%nnnt  = st1%nnnt  * st2%nnnt  
        prodstencil%ennn  = st1%ennn  * st2%ennn  
        prodstencil%wwnnt = st1%wwnnt * st2%wwnnt  
        prodstencil%wwnn  = st1%wwnn  * st2%wwnn  
        prodstencil%eennt = st1%eennt * st2%eennt  
        prodstencil%eenn  = st1%eenn  * st2%eenn  
        prodstencil%nnn   = st1%nnn   * st2%nnn  
        prodstencil%eenbb = st1%eenbb * st2%eenbb  
        prodstencil%eebb  = st1%eebb  * st2%eebb  
        prodstencil%eessb = st1%eessb * st2%eessb  
        prodstencil%eesst = st1%eesst * st2%eesst 
        prodstencil%eess  = st1%eess  * st2%eess  
        prodstencil%eennb = st1%eennb * st2%eennb  
        prodstencil%eeest = st1%eeest * st2%eeest  
        prodstencil%eeent = st1%eeent * st2%eeent  
        prodstencil%eeet  = st1%eeet  * st2%eeet  
        prodstencil%eees  = st1%eees  * st2%eees  
        prodstencil%eeen  = st1%eeen  * st2%eeen  
        prodstencil%eee   = st1%eee   * st2%eee  

      END FUNCTION prodstencil
!-----------------------------------------------------------------------
      FUNCTION maskstencil(st1, st2)
! ... multiply a stencil and a mask, location by location
!
      IMPLICIT NONE
        TYPE(stencil), INTENT(in) :: st1
        TYPE(masks), INTENT(in) :: st2
        TYPE(stencil) :: maskstencil

        maskstencil%c  = st1%c  * st2%c
        maskstencil%e  = st1%e  * st2%e
        maskstencil%w  = st1%w  * st2%w
        maskstencil%ee = st1%ee * st2%ee
        maskstencil%ww = st1%ww * st2%ww
        maskstencil%n  = st1%n  * st2%n
        maskstencil%en = st1%en * st2%en
        maskstencil%wn = st1%wn * st2%wn
        maskstencil%s  = st1%s  * st2%s
        maskstencil%es = st1%es * st2%es
        maskstencil%ws = st1%ws * st2%ws
        maskstencil%nn = st1%nn * st2%nn
        maskstencil%ss = st1%ss * st2%ss
        maskstencil%t  = st1%t  * st2%t
        maskstencil%et = st1%et * st2%et
        maskstencil%wt = st1%wt * st2%wt
        maskstencil%nt = st1%nt * st2%nt
        maskstencil%st = st1%st * st2%st
        maskstencil%b  = st1%b  * st2%b
        maskstencil%eb = st1%eb * st2%eb
        maskstencil%wb = st1%wb * st2%wb
        maskstencil%nb = st1%nb * st2%nb
        maskstencil%sb = st1%sb * st2%sb
        maskstencil%tt = st1%tt * st2%tt
        maskstencil%bb = st1%bb * st2%bb

      END FUNCTION maskstencil
!-----------------------------------------------------------------------
      FUNCTION dotstencil(a, st2)
! ... multiply each stencil location by a constant value a

      IMPLICIT NONE
        REAL*8, INTENT(IN) :: a
        TYPE(stencil), INTENT(IN) :: st2
        TYPE(stencil) :: dotstencil

        dotstencil%c     = a * st2%c
        dotstencil%e     = a * st2%e
        dotstencil%w     = a * st2%w
        dotstencil%ee    = a * st2%ee
        dotstencil%ww    = a * st2%ww
        dotstencil%n     = a * st2%n
        dotstencil%en    = a * st2%en
        dotstencil%wn    = a * st2%wn
        dotstencil%s     = a * st2%s
        dotstencil%es    = a * st2%es
        dotstencil%ws    = a * st2%ws
        dotstencil%nn    = a * st2%nn
        dotstencil%ss    = a * st2%ss
        dotstencil%t     = a * st2%t
        dotstencil%et    = a * st2%et
        dotstencil%wt    = a * st2%wt
        dotstencil%nt    = a * st2%nt
        dotstencil%st    = a * st2%st
        dotstencil%b     = a * st2%b
        dotstencil%eb    = a * st2%eb
        dotstencil%wb    = a * st2%wb
        dotstencil%nb    = a * st2%nb
        dotstencil%sb    = a * st2%sb
        dotstencil%tt    = a * st2%tt
        dotstencil%bb    = a * st2%bb
        dotstencil%eet   = a * st2%eet
        dotstencil%eeb   = a * st2%eeb
        dotstencil%ett   = a * st2%ett
        dotstencil%wtt   = a * st2%wtt
        dotstencil%wst   = a * st2%wst 
        dotstencil%est   = a * st2%est 
        dotstencil%wnt   = a * st2%wnt 
        dotstencil%ent   = a * st2%ent 
        dotstencil%enb   = a * st2%enb 
        dotstencil%esb   = a * st2%esb 
        dotstencil%wnb   = a * st2%wnb 
        dotstencil%sst   = a * st2%sst 
        dotstencil%een   = a * st2%een 
        dotstencil%ees   = a * st2%ees 
        dotstencil%ess   = a * st2%ess 
        dotstencil%esst  = a * st2%esst 
        dotstencil%eenb  = a * st2%eenb 
        dotstencil%eesb  = a * st2%eesb 
        dotstencil%eest  = a * st2%eest 
        dotstencil%wnnt  = a * st2%wnnt  
        dotstencil%nnt   = a * st2%nnt 
        dotstencil%ennb  = a * st2%ennb 
        dotstencil%enn   = a * st2%enn 
        dotstencil%ebb   = a * st2%ebb 
        dotstencil%enbb  = a * st2%enbb 
        dotstencil%wnnb  = a * st2%wnnb 
        dotstencil%wnn   = a * st2%wnn 
        dotstencil%nbb   = a * st2%nbb 
        dotstencil%nnb   = a * st2%nnb 
        dotstencil%wstt  = a * st2%wstt 
        dotstencil%stt   = a * st2%stt 
        dotstencil%wntt  = a * st2%wntt 
        dotstencil%ntt   = a * st2%ntt 
        dotstencil%estt  = a * st2%estt 
        dotstencil%wwn   = a * st2%wwn 
        dotstencil%wwnt  = a * st2%wwnt 
        dotstencil%wwt   = a * st2%wwt 
        dotstencil%wsb   = a * st2%wsb 
        dotstencil%eent  = a * st2%eent 
        dotstencil%ennt  = a * st2%ennt 
        dotstencil%entt  = a * st2%entt 
        dotstencil%wbb   = a * st2%wbb  
        dotstencil%ssb   = a * st2%ssb  
        dotstencil%wws   = a * st2%wws  
        dotstencil%sbb   = a * st2%sbb  
        dotstencil%wss   = a * st2%wss  
        dotstencil%wwb   = a * st2%wwb  
        dotstencil%wnbb  = a * st2%wnbb  
        dotstencil%essb  = a * st2%essb  
        dotstencil%wwst  = a * st2%wwst  
        dotstencil%esbb  = a * st2%esbb  
        dotstencil%wsbb  = a * st2%wsbb  
        dotstencil%wssb  = a * st2%wssb  
        dotstencil%wsst  = a * st2%wsst  
        dotstencil%wwnb  = a * st2%wwnb  
        dotstencil%wwsb  = a * st2%wwsb  
        dotstencil%enttt = a * st2%enttt  
        dotstencil%wnttt = a * st2%wnttt  
        dotstencil%nttt  = a * st2%nttt  
        dotstencil%ettt  = a * st2%ettt  
        dotstencil%wttt  = a * st2%wttt  
        dotstencil%ttt   = a * st2%ttt  
        dotstencil%esstt = a * st2%esstt  
        dotstencil%sstt  = a * st2%sstt  
        dotstencil%enntt = a * st2%enntt  
        dotstencil%nntt  = a * st2%nntt  
        dotstencil%wwstt = a * st2%wwstt  
        dotstencil%wwntt = a * st2%wwntt  
        dotstencil%wwtt  = a * st2%wwtt  
        dotstencil%eestt = a * st2%eestt  
        dotstencil%eentt = a * st2%eentt  
        dotstencil%eett  = a * st2%eett  
        dotstencil%ennbb = a * st2%ennbb  
        dotstencil%wnnbb = a * st2%wnnbb  
        dotstencil%nnbb  = a * st2%nnbb  
        dotstencil%wnntt = a * st2%wnntt  
        dotstencil%ennnb = a * st2%ennnb  
        dotstencil%ennnt = a * st2%ennnt  
        dotstencil%nnnb  = a * st2%nnnb  
        dotstencil%nnnt  = a * st2%nnnt  
        dotstencil%ennn  = a * st2%ennn  
        dotstencil%wwnnt = a * st2%wwnnt  
        dotstencil%wwnn  = a * st2%wwnn  
        dotstencil%eennt = a * st2%eennt  
        dotstencil%eenn  = a * st2%eenn  
        dotstencil%nnn   = a * st2%nnn  
        dotstencil%eenbb = a * st2%eenbb  
        dotstencil%eebb  = a * st2%eebb  
        dotstencil%eessb = a * st2%eessb  
        dotstencil%eesst = a * st2%eesst 
        dotstencil%eess  = a * st2%eess  
        dotstencil%eennb = a * st2%eennb  
        dotstencil%eeest = a * st2%eeest  
        dotstencil%eeent = a * st2%eeent  
        dotstencil%eeet  = a * st2%eeet  
        dotstencil%eees  = a * st2%eees  
        dotstencil%eeen  = a * st2%eeen  
        dotstencil%eee   = a * st2%eee  

      END FUNCTION dotstencil
!-----------------------------------------------------------------------
      FUNCTION fracstencil(st1, st2)
! ... compute the ratio between two stencils, location by location
!
      IMPLICIT NONE
        TYPE(stencil), INTENT(in) :: st1, st2
        TYPE(stencil) :: fracstencil

        fracstencil%c     = st1%c     / st2%c
        fracstencil%e     = st1%e     / st2%e
        fracstencil%w     = st1%w     / st2%w
        fracstencil%ee    = st1%ee    / st2%ee
        fracstencil%ww    = st1%ww    / st2%ww
        fracstencil%n     = st1%n     / st2%n
        fracstencil%en    = st1%en    / st2%en
        fracstencil%wn    = st1%wn    / st2%wn
        fracstencil%s     = st1%s     / st2%s
        fracstencil%es    = st1%es    / st2%es
        fracstencil%ws    = st1%ws    / st2%ws
        fracstencil%nn    = st1%nn    / st2%nn
        fracstencil%ss    = st1%ss    / st2%ss
        fracstencil%t     = st1%t     / st2%t
        fracstencil%et    = st1%et    / st2%et
        fracstencil%wt    = st1%wt    / st2%wt
        fracstencil%nt    = st1%nt    / st2%nt
        fracstencil%st    = st1%st    / st2%st
        fracstencil%b     = st1%b     / st2%b
        fracstencil%eb    = st1%eb    / st2%eb
        fracstencil%wb    = st1%wb    / st2%wb
        fracstencil%nb    = st1%nb    / st2%nb
        fracstencil%sb    = st1%sb    / st2%sb
        fracstencil%tt    = st1%tt    / st2%tt
        fracstencil%bb    = st1%bb    / st2%bb
        fracstencil%eet   = st1%eet   / st2%eet
        fracstencil%eeb   = st1%eeb   / st2%eeb
        fracstencil%ett   = st1%ett   / st2%ett
        fracstencil%wtt   = st1%wtt   / st2%wtt
        fracstencil%wst   = st1%wst   / st2%wst 
        fracstencil%est   = st1%est   / st2%est 
        fracstencil%wnt   = st1%wnt   / st2%wnt 
        fracstencil%ent   = st1%ent   / st2%ent 
        fracstencil%enb   = st1%enb   / st2%enb 
        fracstencil%esb   = st1%esb   / st2%esb 
        fracstencil%wnb   = st1%wnb   / st2%wnb 
        fracstencil%sst   = st1%sst   / st2%sst 
        fracstencil%een   = st1%een   / st2%een 
        fracstencil%ees   = st1%ees   / st2%ees 
        fracstencil%ess   = st1%ess   / st2%ess 
        fracstencil%esst  = st1%esst  / st2%esst 
        fracstencil%eenb  = st1%eenb  / st2%eenb 
        fracstencil%eesb  = st1%eesb  / st2%eesb 
        fracstencil%eest  = st1%eest  / st2%eest 
        fracstencil%wnnt  = st1%wnnt  / st2%wnnt 
        fracstencil%nnt   = st1%nnt   / st2%nnt 
        fracstencil%ennb  = st1%ennb  / st2%ennb 
        fracstencil%enn   = st1%enn   / st2%enn 
        fracstencil%ebb   = st1%ebb   / st2%ebb 
        fracstencil%enbb  = st1%enbb  / st2%enbb 
        fracstencil%wnnb  = st1%wnnb  / st2%wnnb 
        fracstencil%wnn   = st1%wnn   / st2%wnn 
        fracstencil%nbb   = st1%nbb   / st2%nbb 
        fracstencil%nnb   = st1%nnb   / st2%nnb 
        fracstencil%wstt  = st1%wstt  / st2%wstt 
        fracstencil%stt   = st1%stt   / st2%stt 
        fracstencil%wntt  = st1%wntt  / st2%wntt 
        fracstencil%ntt   = st1%ntt   / st2%ntt 
        fracstencil%estt  = st1%estt  / st2%estt 
        fracstencil%wwn   = st1%wwn   / st2%wwn 
        fracstencil%wwnt  = st1%wwnt  / st2%wwnt 
        fracstencil%wwt   = st1%wwt   / st2%wwt 
        fracstencil%wsb   = st1%wsb   / st2%wsb 
        fracstencil%eent  = st1%eent  / st2%eent 
        fracstencil%ennt  = st1%ennt  / st2%ennt 
        fracstencil%entt  = st1%entt  / st2%entt 
        fracstencil%wbb   = st1%wbb   / st2%wbb  
        fracstencil%ssb   = st1%ssb   / st2%ssb  
        fracstencil%wws   = st1%wws   / st2%wws  
        fracstencil%sbb   = st1%sbb   / st2%sbb  
        fracstencil%wss   = st1%wss   / st2%wss  
        fracstencil%wwb   = st1%wwb   / st2%wwb  
        fracstencil%wnbb  = st1%wnbb  / st2%wnbb  
        fracstencil%essb  = st1%essb  / st2%essb  
        fracstencil%wwst  = st1%wwst  / st2%wwst  
        fracstencil%esbb  = st1%esbb  / st2%esbb  
        fracstencil%wsbb  = st1%wsbb  / st2%wsbb  
        fracstencil%wssb  = st1%wssb  / st2%wssb  
        fracstencil%wsst  = st1%wsst  / st2%wsst  
        fracstencil%wwnb  = st1%wwnb  / st2%wwnb  
        fracstencil%wwsb  = st1%wwsb  / st2%wwsb  
        fracstencil%enttt = st1%enttt / st2%enttt  
        fracstencil%wnttt = st1%wnttt / st2%wnttt  
        fracstencil%nttt  = st1%nttt  / st2%nttt  
        fracstencil%ettt  = st1%ettt  / st2%ettt  
        fracstencil%wttt  = st1%wttt  / st2%wttt  
        fracstencil%ttt   = st1%ttt   / st2%ttt  
        fracstencil%esstt = st1%esstt / st2%esstt  
        fracstencil%sstt  = st1%sstt  / st2%sstt  
        fracstencil%enntt = st1%enntt / st2%enntt  
        fracstencil%nntt  = st1%nntt  / st2%nntt  
        fracstencil%wwstt = st1%wwstt / st2%wwstt  
        fracstencil%wwntt = st1%wwntt / st2%wwntt  
        fracstencil%wwtt  = st1%wwtt  / st2%wwtt  
        fracstencil%eestt = st1%eestt / st2%eestt  
        fracstencil%eentt = st1%eentt / st2%eentt  
        fracstencil%eett  = st1%eett  / st2%eett  
        fracstencil%ennbb = st1%ennbb / st2%ennbb  
        fracstencil%wnnbb = st1%wnnbb / st2%wnnbb  
        fracstencil%nnbb  = st1%nnbb  / st2%nnbb  
        fracstencil%wnntt = st1%wnntt / st2%wnntt  
        fracstencil%ennnb = st1%ennnb / st2%ennnb  
        fracstencil%ennnt = st1%ennnt / st2%ennnt  
        fracstencil%nnnb  = st1%nnnb  / st2%nnnb  
        fracstencil%nnnt  = st1%nnnt  / st2%nnnt  
        fracstencil%ennn  = st1%ennn  / st2%ennn  
        fracstencil%wwnnt = st1%wwnnt / st2%wwnnt  
        fracstencil%wwnn  = st1%wwnn  / st2%wwnn  
        fracstencil%eennt = st1%eennt / st2%eennt  
        fracstencil%eenn  = st1%eenn  / st2%eenn  
        fracstencil%nnn   = st1%nnn   / st2%nnn  
        fracstencil%eenbb = st1%eenbb / st2%eenbb  
        fracstencil%eebb  = st1%eebb  / st2%eebb  
        fracstencil%eessb = st1%eessb / st2%eessb  
        fracstencil%eesst = st1%eesst / st2%eesst 
        fracstencil%eess  = st1%eess  / st2%eess  
        fracstencil%eennb = st1%eennb / st2%eennb  
        fracstencil%eeest = st1%eeest / st2%eeest  
        fracstencil%eeent = st1%eeent / st2%eeent  
        fracstencil%eeet  = st1%eeet  / st2%eeet  
        fracstencil%eees  = st1%eees  / st2%eees  
        fracstencil%eeen  = st1%eeen  / st2%eeen  
        fracstencil%eee   = st1%eee   / st2%eee  

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

        stncl%c    = array( ijk     )
        stncl%e    = array( ijke    )
        stncl%t    = array( ijkt    )
        stncl%w    = array( ijkw    )
        stncl%b    = array( ijkb    )
        stncl%et   = array( ijket   )
        stncl%wt   = array( ijkwt   )
        stncl%eb   = array( ijkeb   )
        stncl%wb   = array( ijkwb   )
        stncl%ee   = array( ijkee   )
        stncl%tt   = array( ijktt   )
        stncl%ww   = array( ijkww   )
        stncl%bb   = array( ijkbb   )
        stncl%eet  = array( ijkeet  )
        stncl%eeb  = array( ijkeeb  )
        stncl%ett  = array( ijkett  )
        stncl%ebb  = array( ijkebb  )
        stncl%wtt  = array( ijkwtt  )
        stncl%wbb  = array( ijkwbb  )
        stncl%wwt  = array( ijkwwt  )
        stncl%wwb  = array( ijkwwb  )
        stncl%eee  = array( ijkeee  )
        stncl%ttt  = array( ijkttt  )
        stncl%eeet = array( ijkeeet )
        stncl%ettt = array( ijkettt )
        stncl%eett = array( ijkeett )
        stncl%eebb = array( ijkeebb )
        stncl%wwtt = array( ijkwwtt )

      ELSE IF( job_type_flag == 3 ) THEN

        stncl%c     = array( ijk      )
        stncl%e     = array( ijke     )
        stncl%w     = array( ijkw     )
        stncl%n     = array( ijkn     )
        stncl%s     = array( ijks     )
        stncl%t     = array( ijkt     )
        stncl%b     = array( ijkb     )
        stncl%ee    = array( ijkee    )
        stncl%ww    = array( ijkww    )
        stncl%nn    = array( ijknn    )
        stncl%ss    = array( ijkss    )
        stncl%tt    = array( ijktt    )
        stncl%bb    = array( ijkbb    )
        stncl%en    = array( ijken    )
        stncl%wn    = array( ijkwn    )
        stncl%es    = array( ijkes    )
        stncl%ws    = array( ijkws    )
        stncl%et    = array( ijket    )
        stncl%wt    = array( ijkwt    )
        stncl%nt    = array( ijknt    )
        stncl%st    = array( ijkst    )
        stncl%eb    = array( ijkeb    )
        stncl%wb    = array( ijkwb    )
        stncl%nb    = array( ijknb    )
        stncl%sb    = array( ijksb    )
        stncl%eee   = array( ijkeee   )
        stncl%nnn   = array( ijknnn   )
        stncl%ttt   = array( ijkttt   )
        stncl%ent   = array( ijkent   )
        stncl%enb   = array( ijkenb   )
        stncl%est   = array( ijkest   )
        stncl%esb   = array( ijkesb   )
        stncl%wnt   = array( ijkwnt   )
        stncl%wnb   = array( ijkwnb   )
        stncl%wst   = array( ijkwst   )
        stncl%wsb   = array( ijkwsb   )
        stncl%eet   = array( ijkeet   )
        stncl%eeb   = array( ijkeeb   )
        stncl%een   = array( ijkeen   )
        stncl%ees   = array( ijkees   )
        stncl%enn   = array( ijkenn   )
        stncl%wnn   = array( ijkwnn   )
        stncl%nnt   = array( ijknnt   )
        stncl%nnb   = array( ijknnb   )
        stncl%ett   = array( ijkett   )
        stncl%wtt   = array( ijkwtt   )
        stncl%ntt   = array( ijkntt   )
        stncl%stt   = array( ijkstt   )
        stncl%ess   = array( ijkess   )
        stncl%sst   = array( ijksst   )
        stncl%wwt   = array( ijkwwt   )
        stncl%wwn   = array( ijkwwn   )
        stncl%ebb   = array( ijkebb   )
        stncl%nbb   = array( ijknbb   )
        stncl%wbb   = array( ijkwbb   )
        stncl%ssb   = array( ijkssb   )
        stncl%wws   = array( ijkwws   )
        stncl%sbb   = array( ijksbb   )
        stncl%wss   = array( ijkwss   )
        stncl%wwb   = array( ijkwwb   )
        stncl%eent  = array( ijkeent  )
        stncl%eenb  = array( ijkeenb  )
        stncl%eest  = array( ijkeest  )
        stncl%eesb  = array( ijkeesb  )
        stncl%wnnt  = array( ijkwnnt  )
        stncl%wnnb  = array( ijkwnnb  )
        stncl%ennt  = array( ijkennt  )
        stncl%ennb  = array( ijkennb  )
        stncl%wstt  = array( ijkwstt  )
        stncl%wntt  = array( ijkwntt  )
        stncl%estt  = array( ijkestt  )
        stncl%entt  = array( ijkentt  )
        stncl%wwnt  = array( ijkwwnt  )
        stncl%esst  = array( ijkesst  )
        stncl%enbb  = array( ijkenbb  )
        stncl%wnbb  = array( ijkwnbb  )
        stncl%essb  = array( ijkessb  )
        stncl%wwst  = array( ijkwwst  )
        stncl%esbb  = array( ijkesbb  )
        stncl%wsbb  = array( ijkwsbb  )
        stncl%wssb  = array( ijkwssb  )
        stncl%wsst  = array( ijkwsst  )
        stncl%wwnb  = array( ijkwwnb  )
        stncl%wwsb  = array( ijkwwsb  )
        stncl%nttt  = array( ijknttt  )
        stncl%ettt  = array( ijkettt  )
        stncl%wttt  = array( ijkwttt  )
        stncl%nnnb  = array( ijknnnb  )
        stncl%nnnt  = array( ijknnnt  )
        stncl%ennn  = array( ijkennn  )
        stncl%eeet  = array( ijkeeet  )
        stncl%eees  = array( ijkeees  )
        stncl%eeen  = array( ijkeeen  )
        stncl%sstt  = array( ijksstt  )
        stncl%nntt  = array( ijknntt  )
        stncl%wwtt  = array( ijkwwtt  )
        stncl%eett  = array( ijkeett  )
        stncl%nnbb  = array( ijknnbb  )
        stncl%wwnn  = array( ijkwwnn  )
        stncl%eenn  = array( ijkeenn  )
        stncl%eebb  = array( ijkeebb  )
        stncl%eess  = array( ijkeess  )
        stncl%enttt = array( ijkenttt )
        stncl%wnttt = array( ijkwnttt )
        stncl%esstt = array( ijkesstt )
        stncl%enntt = array( ijkenntt )
        stncl%wwstt = array( ijkwwstt )
        stncl%wwntt = array( ijkwwntt )
        stncl%eestt = array( ijkeestt )
        stncl%eentt = array( ijkeentt )
        stncl%ennbb = array( ijkennbb )
        stncl%wnnbb = array( ijkwnnbb )
        stncl%wnntt = array( ijkwnntt )
        stncl%ennnb = array( ijkennnb )
        stncl%ennnt = array( ijkennnt )
        stncl%wwnnt = array( ijkwwnnt )
        stncl%eennt = array( ijkeennt )
        stncl%eenbb = array( ijkeenbb )
        stncl%eessb = array( ijkeessb )
        stncl%eesst = array( ijkeesst )
        stncl%eennb = array( ijkeennb )
        stncl%eeest = array( ijkeeest )
        stncl%eeent = array( ijkeeent )

      END IF
!
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

        stncl%c    = array( ijk     )
        stncl%e    = array( ipjk    )
        stncl%t    = array( ijkp    )
        stncl%w    = array( imjk    )
        stncl%b    = array( ijk     )
        stncl%et   = array( ipjkp   )
        stncl%wt   = array( imjkp   )
        stncl%eb   = array( ipjkm   )
        stncl%wb   = array( imjkm   )
        stncl%ee   = array( ippjk   )
        stncl%tt   = array( ijkpp   )
        stncl%ww   = array( immjk   )
        stncl%bb   = array( ijkmm   )
        stncl%eet  = array( ippjkp  )
        stncl%eeb  = array( ippjkm  )
        stncl%ett  = array( ipjkpp  )
        stncl%ebb  = array( ipjkmm  )
        stncl%wtt  = array( imjkpp  )
        stncl%wbb  = array( imjkmm  )
        stncl%wwt  = array( immjkp  )
        stncl%wwb  = array( immjkm  )
        stncl%eee  = array( ipppjk  )
        stncl%ttt  = array( ijkppp  )
        stncl%eeet = array( ipppjkp )
        stncl%ettt = array( ipjkppp )
        stncl%eett = array( ippjkpp )
        stncl%eebb = array( ippjkmm )
        stncl%wwtt = array( immjkpp )

      ELSE IF( job_type_flag == 3 ) THEN

        stncl%c     = array( ijk      )
        stncl%e     = array( ipjk     )
        stncl%n     = array( ijpk     )
        stncl%t     = array( ijkp     )
        stncl%w     = array( imjk     )
        stncl%s     = array( ijmk     )
        stncl%b     = array( ijkm     )
        stncl%et    = array( ipjkp    )
        stncl%wt    = array( imjkp    )
        stncl%eb    = array( ipjkm    )
        stncl%wb    = array( imjkm    )
        stncl%nt    = array( ijpkp    )
        stncl%st    = array( ijmkp    )
        stncl%nb    = array( ijpkm    )
        stncl%sb    = array( ijmkm    )
        stncl%en    = array( ipjpk    )
        stncl%wn    = array( imjpk    )
        stncl%es    = array( ipjmk    )
        stncl%ws    = array( imjmk    )
        stncl%ent   = array( ipjpkp   )
        stncl%enb   = array( ipjpkm   )
        stncl%est   = array( ipjmkp   )
        stncl%esb   = array( ipjmkm   )
        stncl%wnt   = array( imjpkp   )
        stncl%wnb   = array( imjpkm   )
        stncl%wst   = array( imjmkp   )
        stncl%wsb   = array( imjmkm   )
        stncl%ee    = array( ippjk    )
        stncl%nn    = array( ijppk    )
        stncl%tt    = array( ijkpp    )
        stncl%ww    = array( immjk    )
        stncl%ss    = array( ijmmk    )
        stncl%bb    = array( ijkmm    )
        stncl%eet   = array( ippjkp   )
        stncl%eeb   = array( ippjkm   )
        stncl%een   = array( ippjpk   )
        stncl%ees   = array( ippjmk   )
        stncl%ett   = array( ipjkpp   )
        stncl%wtt   = array( imjkpp   )
        stncl%ntt   = array( ijpkpp   )
        stncl%stt   = array( ijmkpp   )
        stncl%enn   = array( ipjppk   )
        stncl%wnn   = array( imjppk   )
        stncl%nnt   = array( ijppkp   )
        stncl%nnb   = array( ijppkm   )
        stncl%ess   = array( ipjmmk   )
        stncl%sst   = array( ijmmkp   )
        stncl%wwt   = array( immjkp   )
        stncl%wwn   = array( immjpk   )
        stncl%ebb   = array( ipjkmm   )
        stncl%nbb   = array( ijpkmm   )
        stncl%eee   = array( ipppjk   )
        stncl%nnn   = array( ijpppk   )
        stncl%ttt   = array( ijkppp   )
        stncl%eent  = array( ippjpkp  )
        stncl%eenb  = array( ippjpkm  )
        stncl%eest  = array( ippjmkp  )
        stncl%eesb  = array( ippjmkm  )
        stncl%wnnt  = array( imjppkp  )
        stncl%wnnb  = array( imjppkm  )
        stncl%ennt  = array( ipjppkp  )
        stncl%ennb  = array( ipjppkm  )
        stncl%wstt  = array( imjmkpp  )
        stncl%wntt  = array( imjpkpp  )
        stncl%estt  = array( ipjmkpp  )
        stncl%entt  = array( ipjpkpp  )
        stncl%wwnt  = array( immjpkp  )
        stncl%esst  = array( ipjmmkp  )
        stncl%enbb  = array( ipjpkmm  )
        stncl%wbb   = array( imjkmm   )
        stncl%ssb   = array( ijmmkm   )
        stncl%wws   = array( immjmk   )
        stncl%sbb   = array( ijmkmm   )
        stncl%wss   = array( imjmmk   )
        stncl%wwb   = array( immjkm   )
        stncl%wnbb  = array( imjpkmm  )
        stncl%essb  = array( ipjmmkm  )
        stncl%wwst  = array( immjmkp  )
        stncl%esbb  = array( ipjmkmm  )
        stncl%wsbb  = array( imjmkmm  )
        stncl%wssb  = array( imjmmkm  )
        stncl%wsst  = array( imjmmkp  )
        stncl%wwnb  = array( immjpkm  )
        stncl%wwsb  = array( immjmkm  )
        stncl%enttt = array( ipjpkppp )
        stncl%wnttt = array( imjpkppp )
        stncl%nttt  = array( ijpkppp  )
        stncl%ettt  = array( ipjkppp  )
        stncl%wttt  = array( imjkppp  )
        stncl%esstt = array( ipjmmkpp )
        stncl%sstt  = array( ijmmkpp  )
        stncl%enntt = array( ipjppkpp )
        stncl%nntt  = array( ijppkpp  )
        stncl%wwstt = array( immjmkpp )
        stncl%wwntt = array( immjpkpp )
        stncl%wwtt  = array( immjkpp  )
        stncl%eestt = array( ippjmkpp )
        stncl%eentt = array( ippjpkpp )
        stncl%eett  = array( ippjkpp  )
        stncl%ennbb = array( ipjppkmm )
        stncl%wnnbb = array( imjppkmm )
        stncl%nnbb  = array( ijppkmm  )
        stncl%wnntt = array( imjppkpp )
        stncl%ennnb = array( ipjpppkm )
        stncl%ennnt = array( ipjpppkp )
        stncl%nnnb  = array( ijpppkm  )
        stncl%nnnt  = array( ijpppkp  )
        stncl%ennn  = array( ipjpppk  )
        stncl%wwnnt = array( immjppkp )
        stncl%wwnn  = array( immjppk  )
        stncl%eennt = array( ippjppkp )
        stncl%eenn  = array( ippjppk  )
        stncl%eenbb = array( ippjpkmm )
        stncl%eebb  = array( ippjkmm  )
        stncl%eessb = array( ippjmmkm )
        stncl%eesst = array( ippjmmkp )
        stncl%eess  = array( ippjmmk  )
        stncl%eennb = array( ippjppkm )
        stncl%eeest = array( ipppjmkp )
        stncl%eeent = array( ipppjpkp )
        stncl%eeet  = array( ipppjkp  )
        stncl%eees  = array( ipppjmk  )
        stncl%eeen  = array( ipppjpk  ) 

      END IF
!
      RETURN
      END SUBROUTINE full_rnb
!-----------------------------------------------------------------------
      SUBROUTINE first_nb( stncl, array, ijk )
! ... only first neighbours, considering boundary conditions

      IMPLICIT NONE 
!
      TYPE(stencil) :: stncl
      REAL*8, INTENT(IN) :: array(:)
      INTEGER, INTENT(IN) :: ijk
!
      IF( job_type_flag == 2 ) THEN

        stncl%c = array( ijk  )
        stncl%e = array( ijke )
        stncl%t = array( ijkt )
        stncl%w = array( ijkw )
        stncl%b = array( ijkb )

      ELSE IF( job_type_flag == 3 ) THEN

        stncl%c = array( ijk  )
        stncl%e = array( ijke )
        stncl%w = array( ijkw )
        stncl%n = array( ijkn )
        stncl%s = array( ijks )
        stncl%t = array( ijkt )
        stncl%b = array( ijkb )

      END IF
!
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

        stncl%ee = array( ijkee )
        stncl%tt = array( ijktt )
        stncl%ww = array( ijkww )
        stncl%bb = array( ijkbb )

      ELSE IF( job_type_flag == 3 ) THEN

        stncl%ee = array( ijkee )
        stncl%ww = array( ijkww )
        stncl%nn = array( ijknn )
        stncl%ss = array( ijkss )
        stncl%tt = array( ijktt )
        stncl%bb = array( ijkbb )

      END IF
!
      RETURN
      END SUBROUTINE third_nb
!-----------------------------------------------------------------------
      SUBROUTINE first_rnb(stncl,array,ijk)
! ... only first neighbours, without boundary conditions

      IMPLICIT NONE 
!
      TYPE(stencil) :: stncl
      REAL*8, INTENT(IN) :: array(:)
      INTEGER, INTENT(IN) :: ijk
!
      IF( job_type_flag == 2 ) THEN

        stncl%c = array( ijk  )
        stncl%e = array( ipjk )
        stncl%t = array( ijkp )
        stncl%w = array( imjk )
        stncl%b = array( ijkm )

      ELSE IF( job_type_flag == 3 ) THEN

        stncl%c = array( ijk  )
        stncl%e = array( ipjk )
        stncl%w = array( imjk )
        stncl%n = array( ijpk )
        stncl%s = array( ijmk )
        stncl%t = array( ijkp )
        stncl%b = array( ijkm )

      END IF
!
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

        stncl%ee = array( ippjk )
        stncl%tt = array( ijkpp )
        stncl%ww = array( immjk )
        stncl%bb = array( ijkmm )

      ELSE IF( job_type_flag == 3 ) THEN

        stncl%ee = array( ippjk )
        stncl%ww = array( immjk )
        stncl%nn = array( ijppk )
        stncl%ss = array( ijmmk )
        stncl%tt = array( ijkpp )
        stncl%bb = array( ijkmm )

      END IF
!
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

        stncl%c   = array( ijk    )
        stncl%e   = array( ijke   )
        stncl%t   = array( ijkt   )
        stncl%w   = array( ijkw   )
        stncl%b   = array( ijkb   )
        stncl%et  = array( ijket  )
        stncl%wt  = array( ijkwt  )
        stncl%eb  = array( ijkeb  )
        stncl%ee  = array( ijkee  )
        stncl%tt  = array( ijktt  )
        stncl%eet = array( ijkeet )
        stncl%ett = array( ijkett )
        stncl%eee = array( ijkeee )
        stncl%ttt = array( ijkttt )

      ELSE IF( job_type_flag == 3 ) THEN

        stncl%c   = array( ijk    )
        stncl%e   = array( ijke   )
        stncl%w   = array( ijkw   )
        stncl%n   = array( ijkn   )
        stncl%s   = array( ijks   )
        stncl%t   = array( ijkt   )
        stncl%b   = array( ijkb   )
        stncl%ee  = array( ijkee  )
        stncl%nn  = array( ijknn  )
        stncl%tt  = array( ijktt  )
        stncl%en  = array( ijken  )
        stncl%et  = array( ijket  )
        stncl%es  = array( ijkes  )
        stncl%eb  = array( ijkeb  )
        stncl%nt  = array( ijknt  )
        stncl%wn  = array( ijkwn  )
        stncl%nb  = array( ijknb  )
        stncl%wt  = array( ijkwt  )
        stncl%st  = array( ijkst  )
        stncl%eee = array( ijkeee )
        stncl%nnn = array( ijknnn )
        stncl%ttt = array( ijkttt )
        stncl%een = array( ijkeen )
        stncl%eet = array( ijkeet )
        stncl%enn = array( ijkenn )
        stncl%nnt = array( ijknnt )
        stncl%ett = array( ijkett )
        stncl%ntt = array( ijkntt )

      END IF
!
      RETURN
      END SUBROUTINE nb
!-----------------------------------------------------------------------
      SUBROUTINE rnb_r(stncl,array,ijk)
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

        stncl%c   = array( ijk    )
        stncl%e   = array( ipjk   )
        stncl%t   = array( ijkp   )
        stncl%w   = array( imjk   )
        stncl%b   = array( ijkm   )
        stncl%wt  = array( imjkp  )
        stncl%et  = array( ipjkp  )
        stncl%eb  = array( ipjkm  )
        stncl%ee  = array( ippjk  )
        stncl%tt  = array( ijkpp  )
        stncl%eet = array( ippjkp )
        stncl%ett = array( ipjkpp )
        stncl%eee = array( ipppjk )
        stncl%ttt = array( ijkppp )

      ELSE IF( job_type_flag == 3 ) THEN

        stncl%c   = array( ijk    )
        stncl%e   = array( ipjk   )
        stncl%n   = array( ijpk   )
        stncl%t   = array( ijkp   )
        stncl%w   = array( imjk   )
        stncl%s   = array( ijmk   )
        stncl%b   = array( ijkm   )
        stncl%ee  = array( ippjk  )
        stncl%nn  = array( ijppk  )
        stncl%tt  = array( ijkpp  )
        stncl%en  = array( ipjpk  )
        stncl%et  = array( ipjkp  )
        stncl%es  = array( ipjmk  )
        stncl%eb  = array( ipjkm  )
        stncl%nt  = array( ijpkp  )
        stncl%wn  = array( imjpk  )
        stncl%nb  = array( ijpkm  )
        stncl%wt  = array( imjkp  )
        stncl%st  = array( ijmkp  )
        stncl%eee = array( ipppjk )
        stncl%nnn = array( ijpppk )
        stncl%ttt = array( ijkppp )
        stncl%een = array( ippjpk )
        stncl%eet = array( ippjkp )
        stncl%enn = array( ipjppk )
        stncl%nnt = array( ijppkp )
        stncl%ett = array( ipjkpp )
        stncl%ntt = array( ijpkpp )

      END IF
!
      RETURN
      END SUBROUTINE rnb_r
!----------------------------------------------------------------------
      SUBROUTINE rnb_i(stncl,array,ijk)
! ... This routine compute the complete stencil around a grid point
! ... without considering the boundary conditions
!
      USE dimensions
      IMPLICIT NONE 
!
      TYPE(masks) :: stncl
      INTEGER, INTENT(IN) :: array(:)
      INTEGER, INTENT(IN) :: ijk
!
      IF( job_type_flag == 2 ) THEN

        stncl%c   = array( ijk    )
        stncl%e   = array( ipjk   )
        stncl%t   = array( ijkp   )
        stncl%w   = array( imjk   )
        stncl%b   = array( ijkm   )
        stncl%wt  = array( imjkp  )
        stncl%et  = array( ipjkp  )
        stncl%eb  = array( ipjkm  )
        stncl%ee  = array( ippjk  )
        stncl%tt  = array( ijkpp  )
        stncl%eet = array( ippjkp )
        stncl%ett = array( ipjkpp )
        stncl%eee = array( ipppjk )
        stncl%ttt = array( ijkppp )

      ELSE IF( job_type_flag == 3 ) THEN

        stncl%c   = array( ijk    )
        stncl%e   = array( ipjk   )
        stncl%n   = array( ijpk   )
        stncl%t   = array( ijkp   )
        stncl%w   = array( imjk   )
        stncl%s   = array( ijmk   )
        stncl%b   = array( ijkm   )
        stncl%ee  = array( ippjk  )
        stncl%nn  = array( ijppk  )
        stncl%tt  = array( ijkpp  )
        stncl%en  = array( ipjpk  )
        stncl%et  = array( ipjkp  )
        stncl%es  = array( ipjmk  )
        stncl%eb  = array( ipjkm  )
        stncl%nt  = array( ijpkp  )
        stncl%wn  = array( imjpk  )
        stncl%nb  = array( ijpkm  )
        stncl%wt  = array( imjkp  )
        stncl%st  = array( ijmkp  )
        stncl%eee = array( ipppjk )
        stncl%nnn = array( ijpppk )
        stncl%ttt = array( ijkppp )
        stncl%een = array( ippjpk )
        stncl%eet = array( ippjkp )
        stncl%enn = array( ipjppk )
        stncl%nnt = array( ijppkp )
        stncl%ett = array( ipjkpp )
        stncl%ntt = array( ijpkpp )

      END IF
!
      RETURN
      END SUBROUTINE rnb_i
!-----------------------------------------------------------------------
      SUBROUTINE ctu1_nb( stncl, array, ijk )
! ... This routine compute the stencil around a grid point needed by first order CTU
! ... considering the boundary conditions
!
      IMPLICIT NONE 
!
      TYPE(stencil) :: stncl
      REAL*8, INTENT(IN) :: array(:)
      INTEGER, INTENT(IN) :: ijk

      IF( job_type_flag == 2 ) THEN

        stncl%c   = array( ijk    )
        stncl%e   = array( ijke   )
        stncl%t   = array( ijkt   )
        stncl%w   = array( ijkw   )
        stncl%b   = array( ijkb   )
        stncl%et  = array( ijket  )
        stncl%wt  = array( ijkwt  )
        stncl%eb  = array( ijkeb  )
        stncl%wb  = array( ijkwb  )
        stncl%ee  = array( ijkee  )
        stncl%ww  = array( ijkww  )
        stncl%tt  = array( ijktt  )
        stncl%bb  = array( ijkbb  )
        stncl%eet = array( ijkeet )
        stncl%eeb = array( ijkeeb )
        stncl%ett = array( ijkett )
        stncl%wtt = array( ijkwtt )
        stncl%wwt = array( ijkwwt )
        stncl%wwb = array( ijkwwb )
        stncl%ebb = array( ijkebb )
        stncl%wbb = array( ijkwbb )

      ELSE

        stncl%c     = array( ijk     )
        stncl%e     = array( ijke    )
        stncl%n     = array( ijkn    )
        stncl%t     = array( ijkt    )
        stncl%w     = array( ijkw    )
        stncl%s     = array( ijks    )
        stncl%b     = array( ijkb    )
        stncl%et    = array( ijket   )
        stncl%wt    = array( ijkwt   )
        stncl%eb    = array( ijkeb   )
        stncl%wb    = array( ijkwb   )
        stncl%nt    = array( ijknt   )
        stncl%st    = array( ijkst   )
        stncl%nb    = array( ijknb   )
        stncl%sb    = array( ijksb   )
        stncl%en    = array( ijken   )
        stncl%wn    = array( ijkwn   )
        stncl%es    = array( ijkes   )
        stncl%ws    = array( ijkws   )
        stncl%ent   = array( ijkent  )
        stncl%enb   = array( ijkenb  )
        stncl%est   = array( ijkest  )
        stncl%esb   = array( ijkesb  )
        stncl%wnt   = array( ijkwnt  )
        stncl%wnb   = array( ijkwnb  )
        stncl%wst   = array( ijkwst  )
        stncl%wsb   = array( ijkwsb  )
        stncl%ee    = array( ijkee   )
        stncl%nn    = array( ijknn   )
        stncl%tt    = array( ijktt   )
        stncl%ww    = array( ijkww   )
        stncl%ss    = array( ijkss   )
        stncl%bb    = array( ijkbb   )
        stncl%eet   = array( ijkeet  )
        stncl%eeb   = array( ijkeeb  )
        stncl%een   = array( ijkeen  )
        stncl%ees   = array( ijkees  )
        stncl%ett   = array( ijkett  )
        stncl%wtt   = array( ijkwtt  )
        stncl%ntt   = array( ijkntt  )
        stncl%stt   = array( ijkstt  )
        stncl%enn   = array( ijkenn  )
        stncl%wnn   = array( ijkwnn  )
        stncl%nnt   = array( ijknnt  )
        stncl%nnb   = array( ijknnb  )
        stncl%ess   = array( ijkess  )
        stncl%sst   = array( ijksst  )
        stncl%wwt   = array( ijkwwt  )
        stncl%wwn   = array( ijkwwn  )
        stncl%ebb   = array( ijkebb  )
        stncl%nbb   = array( ijknbb  )
        stncl%eent  = array( ijkeent )
        stncl%eenb  = array( ijkeenb )
        stncl%eest  = array( ijkeest )
        stncl%eesb  = array( ijkeesb )
        stncl%wnnt  = array( ijkwnnt )
        stncl%wnnb  = array( ijkwnnb )
        stncl%ennt  = array( ijkennt )
        stncl%ennb  = array( ijkennb )
        stncl%wstt  = array( ijkwstt )
        stncl%wntt  = array( ijkwntt )
        stncl%estt  = array( ijkestt )
        stncl%entt  = array( ijkentt )
        stncl%wwnt  = array( ijkwwnt )
        stncl%esst  = array( ijkesst )
        stncl%enbb  = array( ijkenbb )
        stncl%wws   = array( ijkwws  )
        stncl%wwb   = array( ijkwwb  )
        stncl%wss   = array( ijkwss  )
        stncl%ssb   = array( ijkssb  )
        stncl%wbb   = array( ijkwbb  )
        stncl%sbb   = array( ijksbb  )
        stncl%wwst  = array( ijkwwst )
        stncl%wwsb  = array( ijkwwsb )
        stncl%wwnb  = array( ijkwwnb )
        stncl%essb  = array( ijkessb )
        stncl%wssb  = array( ijkwssb )
        stncl%wsst  = array( ijkwsst )
        stncl%esbb  = array( ijkesbb )
        stncl%wsbb  = array( ijkwsbb )
        stncl%wnbb  = array( ijkwnbb )
        stncl%eebb  = array( ijkeebb )
        stncl%wwnn  = array( ijkwwnn )
        stncl%sstt  = array( ijksstt )
        stncl%eenbb = array( ijkeenbb )
        stncl%eesbb = array( ijkeesbb )
        stncl%wwnnt = array( ijkwwnnt )
        stncl%wwnnb = array( ijkwwnnb )
        stncl%esstt = array( ijkesstt )
        stncl%wsstt = array( ijkwsstt )

      END IF
!
      RETURN
      END SUBROUTINE ctu1_nb
!-----------------------------------------------------------------------
      SUBROUTINE ctu2_nb( stncl, array, ijk )
! ... This routine compute the stencil around a grid point needed by second order CTU
! ... considering the boundary conditions
!
      IMPLICIT NONE 
!
      TYPE(stencil) :: stncl
      REAL*8, INTENT(IN) :: array(:)
      INTEGER, INTENT(IN) :: ijk

      IF( job_type_flag == 2 ) THEN

        stncl%eee = array( ijkeee )
        stncl%ttt = array( ijkttt )

      ELSE

        stncl%eee = array( ijkeee )
        stncl%nnn = array( ijknnn )
        stncl%ttt = array( ijkttt )

      END IF
!
      RETURN
      END SUBROUTINE ctu2_nb
!-----------------------------------------------------------------------
      SUBROUTINE ctu3_nb( stncl, array, ijk )
! ... This routine compute the stencil around a grid point needed by second order CTU (step 4)
! ... considering the boundary conditions
!
      IMPLICIT NONE 
!
      TYPE(stencil) :: stncl
      REAL*8, INTENT(IN) :: array(:)
      INTEGER, INTENT(IN) :: ijk

      IF( job_type_flag == 2 ) THEN

        stncl%eeet = array( ijkeeet )
        stncl%ettt = array( ijkettt )
        stncl%eett = array( ijkeett )
        stncl%eebb = array( ijkeebb )
        stncl%wwtt = array( ijkwwtt )

      ELSE

        stncl%wbb   = array( ijkwbb   )
        stncl%ssb   = array( ijkssb   )
        stncl%wws   = array( ijkwws   )
        stncl%sbb   = array( ijksbb   )
        stncl%wss   = array( ijkwss   )
        stncl%wwb   = array( ijkwwb   )
        stncl%wnbb  = array( ijkwnbb  )
        stncl%essb  = array( ijkessb  )
        stncl%wwst  = array( ijkwwst  )
        stncl%esbb  = array( ijkesbb  )
        stncl%wsbb  = array( ijkwsbb  )
        stncl%wssb  = array( ijkwssb  )
        stncl%wsst  = array( ijkwsst  )
        stncl%wwnb  = array( ijkwwnb  )
        stncl%wwsb  = array( ijkwwsb  )
        stncl%nttt  = array( ijknttt  )
        stncl%ettt  = array( ijkettt  )
        stncl%wttt  = array( ijkwttt  )
        stncl%nnnb  = array( ijknnnb  )
        stncl%nnnt  = array( ijknnnt  )
        stncl%ennn  = array( ijkennn  )
        stncl%eeet  = array( ijkeeet  )
        stncl%eees  = array( ijkeees  )
        stncl%eeen  = array( ijkeeen  )
        stncl%sstt  = array( ijksstt  )
        stncl%nntt  = array( ijknntt  )
        stncl%wwtt  = array( ijkwwtt  )
        stncl%eett  = array( ijkeett  )
        stncl%nnbb  = array( ijknnbb  )
        stncl%wwnn  = array( ijkwwnn  )
        stncl%eenn  = array( ijkeenn  )
        stncl%eebb  = array( ijkeebb  )
        stncl%eess  = array( ijkeess  )
        stncl%enttt = array( ijkenttt )
        stncl%wnttt = array( ijkwnttt )
        stncl%esstt = array( ijkesstt )
        stncl%enntt = array( ijkenntt )
        stncl%wwstt = array( ijkwwstt )
        stncl%wwntt = array( ijkwwntt )
        stncl%eestt = array( ijkeestt )
        stncl%eentt = array( ijkeentt )
        stncl%ennbb = array( ijkennbb )
        stncl%wnnbb = array( ijkwnnbb )
        stncl%wnntt = array( ijkwnntt )
        stncl%ennnb = array( ijkennnb )
        stncl%ennnt = array( ijkennnt )
        stncl%wwnnt = array( ijkwwnnt )
        stncl%eennt = array( ijkeennt )
        stncl%eenbb = array( ijkeenbb )
        stncl%eessb = array( ijkeessb )
        stncl%eesst = array( ijkeesst )
        stncl%eennb = array( ijkeennb )
        stncl%eeest = array( ijkeeest )
        stncl%eeent = array( ijkeeent )

      END IF
!
      RETURN
      END SUBROUTINE ctu3_nb
!-----------------------------------------------------------------------
      SUBROUTINE ctu1_rnb( stncl, array, ijk )
! ... This routine compute the stencil around a grid point needed by first order CTU
! ... without considering the boundary conditions
!
      USE dimensions
      USE io_files, ONLY: logunit
      IMPLICIT NONE 
!
      TYPE(stencil) :: stncl
      REAL*8, INTENT(IN) :: array(:)
      INTEGER, INTENT(IN) :: ijk

      IF( job_type_flag == 2 ) THEN

        stncl%c   = array( ijk    )
        stncl%e   = array( ipjk   )
        stncl%t   = array( ijkp   )
        stncl%w   = array( imjk   )
        stncl%b   = array( ijkm   )
        stncl%et  = array( ipjkp  )
        stncl%wt  = array( imjkp  )
        stncl%eb  = array( ipjkm  )
        stncl%wb  = array( imjkm  )
        stncl%ee  = array( ippjk  )
        stncl%ww  = array( immjk  )
        stncl%tt  = array( ijkpp  )
        stncl%bb  = array( ijkmm  )
        stncl%eet = array( ippjkp )
        stncl%eeb = array( ippjkm )
        stncl%ett = array( ipjkpp )
        stncl%wtt = array( imjkpp )
        stncl%wwt = array( immjkp )
        stncl%wwb = array( immjkm )
        stncl%ebb = array( ipjkmm )
        stncl%wbb = array( imjkmm )

      ELSE

        stncl%c     = array( ijk     )
        stncl%e     = array( ipjk    )
        stncl%n     = array( ijpk    )
        stncl%t     = array( ijkp    )
        stncl%w     = array( imjk    )
        stncl%s     = array( ijmk    )
        stncl%b     = array( ijkm    )
        stncl%et    = array( ipjkp   )
        stncl%wt    = array( imjkp   )
        stncl%eb    = array( ipjkm   )
        stncl%wb    = array( imjkm   )
        stncl%nt    = array( ijpkp   )
        stncl%st    = array( ijmkp   )
        stncl%nb    = array( ijpkm   )
        stncl%sb    = array( ijmkm   )
        stncl%en    = array( ipjpk   )
        stncl%wn    = array( imjpk   )
        stncl%es    = array( ipjmk   )
        stncl%ws    = array( imjmk   )
        stncl%ent   = array( ipjpkp  )
        stncl%enb   = array( ipjpkm  )
        stncl%est   = array( ipjmkp  )
        stncl%esb   = array( ipjmkm  )
        stncl%wnt   = array( imjpkp  )
        stncl%wnb   = array( imjpkm  )
        stncl%wst   = array( imjmkp  )
        stncl%wsb   = array( imjmkm  )
        stncl%ee    = array( ippjk   )
        stncl%nn    = array( ijppk   )
        stncl%tt    = array( ijkpp   )
        stncl%ww    = array( immjk   )
        stncl%ss    = array( ijmmk   )
        stncl%bb    = array( ijkmm   )
        stncl%eet   = array( ippjkp  )
        stncl%eeb   = array( ippjkm  )
        stncl%een   = array( ippjpk  )
        stncl%ees   = array( ippjmk  )
        stncl%ett   = array( ipjkpp  )
        stncl%wtt   = array( imjkpp  )
        stncl%ntt   = array( ijpkpp  )
        stncl%stt   = array( ijmkpp  )
        stncl%enn   = array( ipjppk  )
        stncl%wnn   = array( imjppk  )
        stncl%nnt   = array( ijppkp  )
        stncl%nnb   = array( ijppkm  )
        stncl%ess   = array( ipjmmk  )
        stncl%sst   = array( ijmmkp  )
        stncl%wwt   = array( immjkp  )
        stncl%wwn   = array( immjpk  )
        stncl%ebb   = array( ipjkmm  )
        stncl%nbb   = array( ijpkmm  )
        stncl%eent  = array( ippjpkp )
        stncl%eenb  = array( ippjpkm )
        stncl%eest  = array( ippjmkp )
        stncl%eesb  = array( ippjmkm )
        stncl%wnnt  = array( imjppkp )
        stncl%wnnb  = array( imjppkm )
        stncl%ennt  = array( ipjppkp )
        stncl%ennb  = array( ipjppkm )
        stncl%wstt  = array( imjmkpp )
        stncl%wntt  = array( imjpkpp )
        stncl%estt  = array( ipjmkpp )
        stncl%entt  = array( ipjpkpp )
        stncl%wwnt  = array( immjpkp )
        stncl%esst  = array( ipjmmkp )
        stncl%enbb  = array( ipjpkmm )
        stncl%wws   = array( immjmk  )
        stncl%wwb   = array( immjkm  )
        stncl%wss   = array( imjmmk  )
        stncl%ssb   = array( ijmmkm  )
        stncl%wbb   = array( imjkmm  )
        stncl%sbb   = array( ijmkmm  )
        stncl%wwst  = array( immjmkp )
        stncl%wwsb  = array( immjmkm )
        stncl%wwnb  = array( immjpkm )
        stncl%essb  = array( ipjmmkm )
        stncl%wssb  = array( imjmmkm )
        stncl%wsst  = array( imjmmkp )
        stncl%esbb  = array( ipjmkmm )
        stncl%wsbb  = array( imjmkmm )
        stncl%wnbb  = array( imjpkmm )
        stncl%eebb  = array( ippjkmm )
        stncl%wwnn  = array( immjppk )
        stncl%sstt  = array( ijmmkpp )
        stncl%eenbb = array( ippjpkmm )
        stncl%eesbb = array( ippjmkmm )
        stncl%wwnnt = array( immjppkp )
        stncl%wwnnb = array( immjppkm )
        stncl%esstt = array( ipjmmkpp )
        stncl%wsstt = array( imjmmkpp )

      END IF
!
      RETURN
      END SUBROUTINE ctu1_rnb
!-----------------------------------------------------------------------
      SUBROUTINE ctu2_rnb( stncl, array, ijk )
! ... This routine compute the stencil around a grid point needed by second order CTU
! ... without considering the boundary conditions
!
      IMPLICIT NONE 
!
      TYPE(stencil) :: stncl
      REAL*8, INTENT(IN) :: array(:)
      INTEGER, INTENT(IN) :: ijk

      IF( job_type_flag == 2 ) THEN

        stncl%eee = array( ipppjk )
        stncl%ttt = array( ijkppp )

      ELSE

        stncl%eee = array( ipppjk )
        stncl%nnn = array( ijpppk )
        stncl%ttt = array( ijkppp )

      END IF
!
      RETURN
      END SUBROUTINE ctu2_rnb
!-----------------------------------------------------------------------
      SUBROUTINE ctu3_rnb( stncl, array, ijk )
! ... This routine compute the stencil around a grid point needed by second order CTU (step 4)
! ... without considering the boundary conditions
!
      IMPLICIT NONE 
!
      TYPE(stencil) :: stncl
      REAL*8, INTENT(IN) :: array(:)
      INTEGER, INTENT(IN) :: ijk

      IF( job_type_flag == 2 ) THEN

        stncl%eeet = array( ipppjkp )
        stncl%ettt = array( ipjkppp )
        stncl%eett = array( ippjkpp )
        stncl%eebb = array( ippjkmm )
        stncl%wwtt = array( immjkpp )

      ELSE

        stncl%wbb   = array( imjkmm   )
        stncl%ssb   = array( ijmmkm   )
        stncl%wws   = array( immjmk   )
        stncl%sbb   = array( ijmkmm   )
        stncl%wss   = array( imjmmk   )
        stncl%wwb   = array( immjkm   )
        stncl%wnbb  = array( imjpkmm  )
        stncl%essb  = array( ipjmmkm  )
        stncl%wwst  = array( immjmkp  )
        stncl%esbb  = array( ipjmkmm  )
        stncl%wsbb  = array( imjmkmm  )
        stncl%wssb  = array( imjmmkm  )
        stncl%wsst  = array( imjmmkp  )
        stncl%wwnb  = array( immjpkm  )
        stncl%wwsb  = array( immjmkm  )
        stncl%enttt = array( ipjpkppp )
        stncl%wnttt = array( imjpkppp )
        stncl%nttt  = array( ijpkppp  )
        stncl%ettt  = array( ipjkppp  )
        stncl%wttt  = array( imjkppp  )
        stncl%esstt = array( ipjmmkpp )
        stncl%sstt  = array( ijmmkpp  )
        stncl%enntt = array( ipjppkpp )
        stncl%nntt  = array( ijppkpp  )
        stncl%wwstt = array( immjmkpp )
        stncl%wwntt = array( immjpkpp )
        stncl%wwtt  = array( immjkpp  )
        stncl%eestt = array( ippjmkpp )
        stncl%eentt = array( ippjpkpp )
        stncl%eett  = array( ippjkpp  )
        stncl%ennbb = array( ipjppkmm )
        stncl%wnnbb = array( imjppkmm )
        stncl%nnbb  = array( ijppkmm  )
        stncl%wnntt = array( imjppkpp )
        stncl%ennnb = array( ipjpppkm )
        stncl%ennnt = array( ipjpppkp )
        stncl%nnnb  = array( ijpppkm  )
        stncl%nnnt  = array( ijpppkp  )
        stncl%ennn  = array( ipjpppk  )
        stncl%wwnnt = array( immjppkp )
        stncl%wwnn  = array( immjppk  )
        stncl%eennt = array( ippjppkp )
        stncl%eenn  = array( ippjppk  )
        stncl%eenbb = array( ippjpkmm )
        stncl%eebb  = array( ippjkmm  )
        stncl%eessb = array( ippjmmkm )
        stncl%eesst = array( ippjmmkp )
        stncl%eess  = array( ippjmmk  )
        stncl%eennb = array( ippjppkm )
        stncl%eeest = array( ipppjmkp )
        stncl%eeent = array( ipppjpkp )
        stncl%eeet  = array( ipppjkp  )
        stncl%eees  = array( ipppjmk  )
        stncl%eeen  = array( ipppjpk  ) 

      END IF
!
      RETURN
      END SUBROUTINE ctu3_rnb
!-----------------------------------------------------------------------
      END MODULE set_indexes
!-----------------------------------------------------------------------
