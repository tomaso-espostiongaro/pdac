MODULE indijk_module

  INTEGER, PARAMETER :: ip0_jp0_kp0_ =   1  ! ijk 
  INTEGER, PARAMETER :: ip1_jp0_kp0_ =   2  ! ipjk 
  INTEGER, PARAMETER :: im1_jp0_kp0_ =   3  ! imjk 
  INTEGER, PARAMETER :: ip2_jp0_kp0_ =   4  ! ippjk 
  INTEGER, PARAMETER :: im2_jp0_kp0_ =   5  ! immjk 
  INTEGER, PARAMETER :: ip0_jp1_kp0_ =   6  ! ijpk 
  INTEGER, PARAMETER :: ip1_jp1_kp0_ =   7  ! ipjpk 
  INTEGER, PARAMETER :: im1_jp1_kp0_ =   8  ! imjpk 
  INTEGER, PARAMETER :: ip0_jm1_kp0_ =   9  ! ijmk 
  INTEGER, PARAMETER :: ip1_jm1_kp0_ =  10  ! ipjmk 
  INTEGER, PARAMETER :: im1_jm1_kp0_ =  11  ! imjmk 
  INTEGER, PARAMETER :: ip0_jp2_kp0_ =  12  ! ijppk 
  INTEGER, PARAMETER :: ip0_jm2_kp0_ =  13  ! ijmmk 
  INTEGER, PARAMETER :: ip0_jp0_kp1_ =  14  ! ijkp 
  INTEGER, PARAMETER :: ip1_jp0_kp1_ =  15  ! ipjkp 
  INTEGER, PARAMETER :: im1_jp0_kp1_ =  16  ! imjkp 
  INTEGER, PARAMETER :: ip0_jp1_kp1_ =  17  ! ijpkp 
  INTEGER, PARAMETER :: ip0_jm1_kp1_ =  18  ! ijmkp 
  INTEGER, PARAMETER :: ip0_jp0_km1_ =  19  ! ijkm 
  INTEGER, PARAMETER :: ip1_jp0_km1_ =  20  ! ipjkm 
  INTEGER, PARAMETER :: im1_jp0_km1_ =  21  ! imjkm 
  INTEGER, PARAMETER :: ip0_jp1_km1_ =  22  ! ijpkm 
  INTEGER, PARAMETER :: ip0_jm1_km1_ =  23  ! ijmkm 
  INTEGER, PARAMETER :: ip0_jp0_kp2_ =  24  ! ijkpp 
  INTEGER, PARAMETER :: ip0_jp0_km2_ =  25  ! ijkmm 


! Stencil Dimension
!
  INTEGER, PARAMETER :: nstdim = 25

  INTEGER :: indijk(-2:2,-2:2,-2:2)


  INTEGER, PARAMETER :: bb_ = 1
  INTEGER, PARAMETER :: bl_ = 2
  INTEGER, PARAMETER :: b_  = 3
  INTEGER, PARAMETER :: br_ = 4
  INTEGER, PARAMETER :: ll_ = 5
  INTEGER, PARAMETER :: l_  = 6
  INTEGER, PARAMETER :: r_  = 7
  INTEGER, PARAMETER :: rr_ = 8
  INTEGER, PARAMETER :: tl_ = 9
  INTEGER, PARAMETER :: t_  = 10
  INTEGER, PARAMETER :: tr_ = 11
  INTEGER, PARAMETER :: tt_ = 12


CONTAINS

  SUBROUTINE indijk_setup()
    indijk = 0
    indijk(  0,  0,  0 ) = ip0_jp0_kp0_
    indijk(  1,  0,  0 ) = ip1_jp0_kp0_
    indijk( -1,  0,  0 ) = im1_jp0_kp0_
    indijk(  2,  0,  0 ) = ip2_jp0_kp0_
    indijk( -2,  0,  0 ) = im2_jp0_kp0_
    indijk(  0,  1,  0 ) = ip0_jp1_kp0_
    indijk(  1,  1,  0 ) = ip1_jp1_kp0_
    indijk( -1,  1,  0 ) = im1_jp1_kp0_
    indijk(  0, -1,  0 ) = ip0_jm1_kp0_
    indijk(  1, -1,  0 ) = ip1_jm1_kp0_
    indijk( -1, -1,  0 ) = im1_jm1_kp0_
    indijk(  0,  2,  0 ) = ip0_jp2_kp0_
    indijk(  0, -2,  0 ) = ip0_jm2_kp0_
    indijk(  0,  0,  1 ) = ip0_jp0_kp1_
    indijk(  1,  0,  1 ) = ip1_jp0_kp1_
    indijk( -1,  0,  1 ) = im1_jp0_kp1_
    indijk(  0,  1,  1 ) = ip0_jp1_kp1_
    indijk(  0, -1,  1 ) = ip0_jm1_kp1_
    indijk(  0,  0, -1 ) = ip0_jp0_km1_
    indijk(  1,  0, -1 ) = ip1_jp0_km1_
    indijk( -1,  0, -1 ) = im1_jp0_km1_
    indijk(  0,  1, -1 ) = ip0_jp1_km1_
    indijk(  0, -1, -1 ) = ip0_jm1_km1_
    indijk(  0,  0,  2 ) = ip0_jp0_kp2_
    indijk(  0,  0, -2 ) = ip0_jp0_km2_
    RETURN
  END SUBROUTINE 

END MODULE
