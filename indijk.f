MODULE indijk_module
  INTEGER, PARAMETER :: ip0_jp0_kp0_ =   1
  INTEGER, PARAMETER :: ip1_jp0_kp0_ =   2
  INTEGER, PARAMETER :: im1_jp0_kp0_ =   3
  INTEGER, PARAMETER :: ip2_jp0_kp0_ =   4
  INTEGER, PARAMETER :: im2_jp0_kp0_ =   5
  INTEGER, PARAMETER :: ip0_jp1_kp0_ =   6
  INTEGER, PARAMETER :: ip1_jp1_kp0_ =   7
  INTEGER, PARAMETER :: im1_jp1_kp0_ =   8
  INTEGER, PARAMETER :: ip0_jm1_kp0_ =   9
  INTEGER, PARAMETER :: ip1_jm1_kp0_ =  10
  INTEGER, PARAMETER :: im1_jm1_kp0_ =  11
  INTEGER, PARAMETER :: ip0_jp2_kp0_ =  12
  INTEGER, PARAMETER :: ip0_jm2_kp0_ =  13
  INTEGER, PARAMETER :: ip0_jp0_kp1_ =  14
  INTEGER, PARAMETER :: ip1_jp0_kp1_ =  15
  INTEGER, PARAMETER :: im1_jp0_kp1_ =  16
  INTEGER, PARAMETER :: ip0_jp1_kp1_ =  17
  INTEGER, PARAMETER :: ip0_jm1_kp1_ =  18
  INTEGER, PARAMETER :: ip0_jp0_km1_ =  19
  INTEGER, PARAMETER :: ip1_jp0_km1_ =  20
  INTEGER, PARAMETER :: im1_jp0_km1_ =  21
  INTEGER, PARAMETER :: ip0_jp1_km1_ =  22
  INTEGER, PARAMETER :: ip0_jm1_km1_ =  23
  INTEGER, PARAMETER :: ip0_jp0_kp2_ =  24
  INTEGER, PARAMETER :: ip0_jp0_km2_ =  25

  INTEGER, PARAMETER :: nijk_        =  25

  INTEGER :: indijk(-2:2,-2:2,-2:2)

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
