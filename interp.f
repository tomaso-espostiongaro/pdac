!----------------------------------------------------------------------
      MODULE interpolate_fields
!----------------------------------------------------------------------
!
      USE dimensions
      USE grid, ONLY: dx, dy, dz

      IMPLICIT NONE
!
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE interpolate_x(field, staggered_field, i)
      USE set_indexes, ONLY: stencil
!
! ... interpolate density on the mesh staggered along x
!
      IMPLICIT NONE

      TYPE(stencil), INTENT(IN) :: field
      TYPE(stencil), INTENT(OUT) :: staggered_field
      INTEGER, INTENT(IN) :: i

      INTEGER :: ip2

      REAL*8 :: dxm, dxp, dxpp, indxpp, indxp, indxm
!
      ip2 = MIN( nx, i+2 )

      dxm  = dx(i)   + dx(i-1)
      dxp  = dx(i)   + dx(i+1)
      dxpp = dx(i+1) + dx(ip2)

      indxm  = 1.D0 / dxm
      indxp  = 1.D0 / dxp
      indxpp = 1.D0 / dxpp
!       
! ... Compute linearly interpolated values of density on the staggered grid
!
      staggered_field%c = (dx(i+1) * field%c + dx(i)   * field%e ) * indxp
      staggered_field%e = (dx(ip2) * field%e + dx(i+1) * field%ee) * indxpp
      staggered_field%n = (dx(i+1) * field%n + dx(i)   * field%en) * indxp
      staggered_field%t = (dx(i+1) * field%t + dx(i)   * field%et) * indxp
      staggered_field%w = (dx(i)   * field%w + dx(i-1) * field%c ) * indxm
      staggered_field%s = (dx(i+1) * field%s + dx(i)   * field%es) * indxp
      staggered_field%b = (dx(i+1) * field%b + dx(i)   * field%eb) * indxp
!
      RETURN
      END SUBROUTINE interpolate_x
!----------------------------------------------------------------------
      SUBROUTINE interpolate_y(field, staggered_field, j)
      USE set_indexes, ONLY: stencil
!
! ... interpolate density on the mesh staggered along x
!
      IMPLICIT NONE
!
      TYPE(stencil), INTENT(IN) :: field
      TYPE(stencil), INTENT(OUT) :: staggered_field
      INTEGER, INTENT(IN) :: j

      REAL*8 :: dym, dyp, dypp, indypp, indyp, indym

      INTEGER :: jp2
!
      jp2 = MIN( ny, j+2 )

      dym = dy(j) + dy(j-1)
      dyp = dy(j) + dy(j+1)
      dypp = dy(j+1) + dy(jp2)

      indym=1.D0 / dym
      indyp=1.D0 / dyp
      indypp=1.D0 / dypp
!       
! ... Compute linearly interpolated values of density on the staggered grid
!
      staggered_field%c = (dy(j+1) * field%c + dy(j)   * field%n ) * indyp
      staggered_field%e = (dy(j+1) * field%e + dy(j)   * field%en) * indyp
      staggered_field%n = (dy(jp2) * field%n + dy(j+1) * field%nn) * indypp
      staggered_field%t = (dy(j+1) * field%t + dy(j)   * field%nt) * indyp
      staggered_field%w = (dy(j+1) * field%w + dy(j)   * field%wn) * indyp
      staggered_field%s = (dy(j)   * field%s + dy(j-1) * field%c ) * indym
      staggered_field%b = (dy(j+1) * field%b + dy(j)   * field%nb) * indyp
!
      RETURN
      END SUBROUTINE interpolate_y
!----------------------------------------------------------------------
      SUBROUTINE interpolate_z(field, staggered_field, k)
      USE set_indexes, ONLY: stencil
!
! ... interpolate density on the mesh staggered along x
!
      IMPLICIT NONE
!
      TYPE(stencil), INTENT(IN) :: field
      TYPE(stencil), INTENT(OUT) :: staggered_field
      INTEGER, INTENT(IN) :: k

      REAL*8 :: dzp, indzp, dzm, indzm, dzpp, indzpp

      INTEGER :: kp2
!
      kp2 = MIN( nz, k+2 )

      dzm=dz(k)+dz(k-1)
      dzp=dz(k)+dz(k+1)
      dzpp=dz(k+1)+dz(kp2)

      indzm=1.D0/dzm
      indzp=1.D0/dzp
      indzpp=1.D0/dzpp
!       
! ... values of density interpolated linearly on the staggered grid
!
      staggered_field%c = (dz(k+1) * field%c + dz(k)   * field%t ) * indzp
      staggered_field%e = (dz(k+1) * field%e + dz(k)   * field%et) * indzp
      staggered_field%n = (dz(k+1) * field%n + dz(k)   * field%nt) * indzp
      staggered_field%t = (dz(kp2) * field%t + dz(k+1) * field%tt) * indzpp
      staggered_field%w = (dz(k+1) * field%w + dz(k)   * field%wt) * indzp
      staggered_field%s = (dz(k+1) * field%s + dz(k)   * field%st) * indzp
      staggered_field%b = (dz(k)   * field%b + dz(k-1) * field%c ) * indzm
!
      RETURN
      END SUBROUTINE interpolate_z
!-----------------------------------------------------------------------
!
!     I M M E R S E D   B O U N D A R I E S   
!
!     I N T E R P O L A T I O N     P R O C E D U R E S 
!
!-----------------------------------------------------------------------
      REAL*8 FUNCTION velint(fpt, vel, ijk, cx, cy, cz)
!
! ... Interpolate velocities on a forcing point to get no-slip
! ... conditions on a solid immersed boundary 

      USE set_indexes, ONLY: ipjk, imjk, ippjk, immjk, ijpk, ipjpk,    &
        imjpk, ijmk, ipjmk, imjmk, ijppk, ijmmk, ijkp, ipjkp, imjkp,   &
        ijpkp, ijmkp, ijkm, ipjkm, imjkm, ijpkm, ijmkm, ijkpp, ijkmm
      USE immersed_boundaries, ONLY: forcing_point
      IMPLICIT NONE

      TYPE(forcing_point), INTENT(IN) :: fpt
      REAL*8, DIMENSION(:), INTENT(IN) :: vel
      REAL*8, DIMENSION(:), INTENT(IN) :: cx, cy, cz
      INTEGER, INTENT(IN) :: ijk

      INTEGER :: i, j, k
      REAL*8 :: nsx, nsy, nsz
      INTEGER :: immjkpp, ippjkpp
      INTEGER :: interp

      REAL*8 :: h          !distance between the (i,j,k)-node and boundary
      REAL*8 :: zA         !distance between the boundary and the first 
                           !external node
      REAL*8 :: zB         !distance between the boundary and the second 
                           !external node
      REAL*8 :: alpha,beta !coefficients for bilinear interpolation

      immjkpp = imjkp
      ippjkpp = ipjkp
      
      interp = fpt%int
      i      = fpt%i
      j      = fpt%j
      k      = fpt%k
      nsx    = fpt%nsl%x
      nsy    = fpt%nsl%y
      nsz    = fpt%nsl%z

      SELECT CASE (interp)


      CASE (-3)

!===========================================
!====   interpolazione lineare con i	====
!====   valori nei nodi NW-NNWW		====
!===========================================	
   
         h=SQRT((cx(i)-nsx)**2+(cz(k)-nsz)**2)
         zA=SQRT((cx(i-1)-nsx)**2+(cz(k+1)-nsz)**2)
         IF (h <= zA) THEN
            velint=-h/zA*vel(imjkp)
         ELSE
            zB=SQRT((cx(i-2)-nsx)**2+(cz(k+2)-nsz)**2)
            velint=-((zB-h)*vel(imjkp)+(h-zA)*vel(immjkpp))/(zB-zA)
         ENDIF

      CASE (-2)
   
!===========================================
!====   interpolazione lineare sin. con ====
!====   i valori nei nodi W-WW		====
!===========================================

         h=SQRT((cx(i)-nsx)**2+(cz(k)-nsz)**2)
         zA=SQRT((cx(i-1)-nsx)**2+(cz(k)-nsz)**2)
         IF (h <= zA) THEN
            velint=-h/zA*vel(imjk)
         ELSE
            zB=SQRT((cx(i-2)-nsx)**2+(cz(k)-nsz)**2)
            velint=-((zB-h)*vel(imjk)+(h-zA)*vel(immjk))/(zB-zA)
         ENDIF

      CASE (-1)
   
!===========================================
!====   interpolazione bilin. con i	====
!====   valori nei nodi W-NW-N		====
!===========================================

         alpha=(cx(i-1)-nsx)/(cx(i-1)-cx(i))
         beta=(cz(k+1)-nsz)/(cz(k+1)-cz(k))
         velint=-(alpha*(1-beta)*vel(ijkp)+(1-alpha)*(1-beta)*vel(imjkp)+ &
                 (1-alpha)*beta*vel(imjk))/(alpha*beta)
!
      CASE (0)
   
!===========================================
!====   interpolazione lineare con i	====
!====   valori nei nodi N-NN		====
!===========================================
   
         h=SQRT((cx(i)-nsx)**2+(cz(k)-nsz)**2)
         zA=SQRT((cx(i)-nsx)**2+(cz(k+1)-nsz)**2)
         IF (h <= zA) THEN
            velint=-h/zA*vel(ijkp)
         ELSE
            zB=SQRT((cx(i)-nsx)**2+(cz(k+2)-nsz)**2)
            velint=-((zB-h)*vel(ijkp)+(h-zA)*vel(ijkpp))/(zB-zA)

! ... Quadratic interpolation
!            velint=- ((vel(ijkp)*zB-vel(ijkpp)*zA)/(zA*zB*(zA-zB))*h**2   &
!                    +(vel(ijkpp)*zA**2-vel(ijkp)*zB**2)/(zA*zB*(zA-zB))*h)
         ENDIF

      CASE (1)
   
!===========================================
!====   interpolazione bilin. con i	====
!====   valori nei nodi E-NE-N		====
!===========================================
   
         alpha=(cx(i+1)-nsx)/(cx(i+1)-cx(i))
         beta=(cz(k+1)-nsz)/(cz(k+1)-cz(k))
         velint=-(alpha*(1-beta)*vel(ijkp)+(1-alpha)*(1-beta)*vel(ipjkp)+ &
                 (1-alpha)*beta*vel(ipjk))/(alpha*beta)
         
      CASE (2)
   
!===========================================
!====   interpolazione lineare con i	====
!====   valori nei nodi E-EE		====
!===========================================
   
         h=SQRT((cx(i)-nsx)**2+(cz(k)-nsz)**2)
         zA=SQRT((cx(i+1)-nsx)**2+(cz(k)-nsz)**2)
         IF (h <= zA) THEN
            velint=-h/zA*vel(ipjk)
         ELSE
            zB=SQRT((cx(i+2)-nsx)**2+(cz(k)-nsz)**2)
            velint=-((zB-h)*vel(ipjk)+(h-zA)*vel(ippjk))/(zB-zA)
         ENDIF
         
      CASE (3)
   
!===========================================
!====   interpolazione lineare con i	====
!====   valori nei nodi NE-NNEE		====
!===========================================
   
         h=SQRT((cx(i)-nsx)**2+(cz(k)-nsz)**2)
         zA=SQRT((cx(i+1)-nsx)**2+(cz(k+1)-nsz)**2)
         IF (h <= zA) THEN
            velint=-h/zA*vel(ipjkp)
         ELSE
            zB=SQRT((cx(i+2)-nsx)**2+(cz(k+2)-nsz)**2)
            velint=-((zB-h)*vel(ipjkp)+(h-zA)*vel(ippjkpp))/(zB-zA)
         ENDIF
   
      CASE (22)
   
!===========================================
!====   interpolazione lineare      	====
!====   esterna E           		====
!===========================================
   
         IF (nsx < cx(i)) THEN
           h=SQRT((cx(i)-nsx)**2+(cz(k)-nsz)**2)
           zA=SQRT((cx(i+1)-nsx)**2+(cz(k)-nsz)**2)
           velint= + h/zA*vel(ipjk)
         ELSE
                 CALL error('bdry','control fpoint',ijk)
         END IF

      CASE (21)
   
!===========================================
!====   interpolazione lineare      	====
!====   esterna             		====
!===========================================
   
         IF (nsx > cx(i)) THEN
           h=SQRT((cx(i)-nsx)**2+(cz(k)-nsz)**2)
           zA=SQRT((cx(i-1)-nsx)**2+(cz(k)-nsz)**2)
           velint= + h/zA*vel(imjk)
         ELSE
                 CALL error('bdry','control fpoint',ijk)
         END IF

      CASE (20)
   
!===========================================
!====   interpolazione lineare      	====
!====   esterna TOP         		====
!===========================================
   
         h=SQRT((cx(i)-nsx)**2+(cz(k)-nsz)**2)
         zA=SQRT((cx(i)-nsx)**2+(cz(k+1)-nsz)**2)
         velint= + h/zA*vel(ijkp)

! ... Quadratic interpolation
!         zB=SQRT((cx(i)-nsx)**2+(cz(k+2)-nsz)**2)
!         velint= (vel(ijkp)*zB-vel(ijkpp)*zA)/(zA*zB*(zA-zB))*h**2 &
!                +(vel(ijkpp)*zA**2-vel(ijkp)*zB**2)/(zA*zB*(zA-zB))*h

      CASE DEFAULT
   
!===========================================
!====   nessuna interpolazione  ============
!===========================================
   
         velint=vel(ijk)

      END SELECT
        
      END FUNCTION velint
!----------------------------------------------------------------------
      REAL*8 FUNCTION velint3d(fpt, vel, ijk, cx, cy, cz, index_q)
!
! ... Interpolate velocities on a forcing point to get no-slip
! ... conditions on a solid immersed boundary 

      USE set_indexes
      USE immersed_boundaries, ONLY: forcing_point
      IMPLICIT NONE

      TYPE(forcing_point), INTENT(IN) :: fpt
      REAL*8, DIMENSION(:), INTENT(IN) :: vel
      REAL*8, DIMENSION(:), INTENT(IN) :: cx, cy, cz
      INTEGER, INTENT(IN) :: ijk
      INTEGER, INTENT(OUT) :: index_q

      INTEGER :: i, j, k
      REAL*8 :: nsx, nsy, nsz
      INTEGER :: interp, delta_i, delta_j, delta_k
      INTEGER :: index_qq
      LOGICAL :: diagonal

      REAL*8 :: h          !distance between the (i,j,k)-node and boundary
      REAL*8 :: zA         !distance between the boundary and the first 
                           !external node
      REAL*8 :: zB         !distance between the boundary and the second 
                           !external node
      REAL*8 :: alpha,beta !coefficients for bilinear interpolation

      velint3d = 0.D0

      interp = fpt%int
      i      = fpt%i
      j      = fpt%j
      k      = fpt%k
      nsx    = fpt%nsl%x
      nsy    = fpt%nsl%y
      nsz    = fpt%nsl%z
!
! ... LINEAR interpolation criteria: 
! ... 0: top
! ... 1: west
! ... 2: south-west
! ... 3: south
! ... 4: south-east
! ... 5: east
! ... 6: north-east
! ... 7: north
! ... 8: north-west

      SELECT CASE (interp)

        CASE (0)
              delta_i = 0
              delta_j = 0
              delta_k = 1
              index_q  = ijkp
              index_qq = ijkpp

        CASE (1)
              delta_i = -1
              delta_j = 0
              delta_k = 0
              index_q  = imjk
              index_qq = immjk

        CASE (2)
              delta_i = -1
              delta_j = -1
              delta_k = 0
              index_q  = imjmk

        CASE (3)
              delta_i = 0
              delta_j = -1
              delta_k = 0
              index_q  = ijmk
              index_qq = ijmmk

        CASE (4)
              delta_i = 1
              delta_j = -1
              delta_k = 0
              index_q  = ipjmk

        CASE (5)
              delta_i = 1
              delta_j = 0
              delta_k = 0
              index_q  = ipjk
              index_qq = ippjk

        CASE (6)
              delta_i = 1
              delta_j = 1
              delta_k = 0
              index_q  = ipjpk

        CASE (7)
              delta_i = 0
              delta_j = 1
              delta_k = 0
              index_q  = ijpk
              index_qq = ijppk

        CASE (8)
              delta_i = -1
              delta_j = 1
              delta_k = 0
              index_q  = imjpk

        ! ... external NE
        CASE(26)
              delta_i = 1
              delta_j = 1
              delta_k = 0
              index_q  = ipjpk
        ! ... external E
        CASE(25)
              delta_i = 1
              delta_j = 0
              delta_k = 0
              index_q  = ipjk
        ! ... external SE
        CASE(24)
              delta_i = 1
              delta_j = -1
              delta_k = 0
              index_q  = ipjmk
        ! ... external N
        CASE(27)
              delta_i = 0
              delta_j = 1
              delta_k = 0
              index_q  = ijpk              
        ! ... external TOP
        CASE(20)
              delta_i = 0
              delta_j = 0
              delta_k = 1
              index_q  = ijkp
        ! ... external S
        CASE(23)
              delta_i = 0
              delta_j = -1
              delta_k = 0
              index_q  = ijmk  
        ! ... external NW
        CASE(28)
              delta_i = -1
              delta_j = 1
              delta_k = 0
              index_q  = imjpk              
        ! ... external W
        CASE(21)
              delta_i = -1
              delta_j = 0
              delta_k = 0
              index_q  = imjk              
        ! ... external SW
        CASE(22)
              delta_i = -1
              delta_j = -1
              delta_k = 0
              index_q  = imjmk              

        CASE DEFAULT
              delta_i = 0
              delta_j = 0
              delta_k = 0
              index_q  = ijk
              index_qq = ijk

      END SELECT

      diagonal = ( (MOD(interp,2) == 0) .AND. (interp/=0) )

      h  = SQRT( (cx(i)-nsx)**2 + (cy(j)-nsy)**2 + (cz(k)-nsz)**2 )
      
      zA = SQRT( (cx(i+delta_i)-nsx)**2 + (cy(j+delta_j)-nsy)**2 + &
                 (cz(k+delta_k)-nsz)**2 )

      IF (interp >= 20) THEN
         velint3d = h/zA*vel(index_q)         
      ELSEIF (h <= zA .OR. diagonal) THEN
         velint3d = -h/zA*vel(index_q)
      ELSE 
         zB=SQRT( (cx(i+2*delta_i)-nsx)**2 + (cy(j+2*delta_j)-nsy)**2 +  &
                  (cz(k+2*delta_k)-nsz)**2 )
         velint3d = -((zB-h)*vel(index_q)+(h-zA)*vel(index_qq))/(zB-zA)
      ENDIF
        
      RETURN
      END FUNCTION velint3d
!----------------------------------------------------------------------
      END MODULE interpolate_fields
!----------------------------------------------------------------------
