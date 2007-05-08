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
      FUNCTION velint(fpt, velg, vels, ijk, cx, cy, cz) 
!
! ... Interpolate velocities on a forcing point to get no-slip
! ... conditions on a solid immersed boundary 

      USE dimensions, ONLY: nsolid, max_nsolid
      USE set_indexes, ONLY: ipjk, imjk, ippjk, immjk, ijpk, ipjpk,    &
        imjpk, ijmk, ipjmk, imjmk, ijppk, ijmmk, ijkp, ipjkp, imjkp,   &
        ijpkp, ijmkp, ijkm, ipjkm, imjkm, ijpkm, ijmkm, ijkpp, ijkmm
      USE immersed_boundaries, ONLY: forcing_point
      IMPLICIT NONE

      REAL*8, DIMENSION(max_nsolid+1) :: velint
      TYPE(forcing_point), INTENT(IN) :: fpt
      REAL*8, DIMENSION(:), INTENT(IN) :: velg
      REAL*8, DIMENSION(:,:), INTENT(IN) :: vels
      REAL*8, DIMENSION(:), INTENT(IN) :: cx, cy, cz
      INTEGER, INTENT(IN) :: ijk

      INTEGER :: i, j, k, is
      REAL*8 :: nsx, nsy, nsz
      INTEGER :: immjkpp, ippjkpp
      INTEGER :: interp

      REAL*8 :: h          !distance between the (i,j,k)-node and boundary
      REAL*8 :: zA         !distance between the boundary and the first 
                           !external node
      REAL*8 :: zB         !distance between the boundary and the second 
                           !external node
      REAL*8 :: alpha,beta !coefficients for bilinear interpolation
!
      velint = 0.D0
!
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
!====   valori nei nodi TW-TTWW		====
!===========================================	
   
         h=SQRT((cx(i)-nsx)**2+(cz(k)-nsz)**2)
         zA=SQRT((cx(i-1)-nsx)**2+(cz(k+1)-nsz)**2)
         IF (h <= zA) THEN
            velint(1)=-h/zA*velg(imjkp)
            DO is = 1, nsolid
              velint(is+1)=-h/zA*vels(imjkp,is)
            END DO
         ELSE
            zB=SQRT((cx(i-2)-nsx)**2+(cz(k+2)-nsz)**2)
            velint(1)=-((zB-h)*velg(imjkp)+(h-zA)*velg(immjkpp))/(zB-zA)
            DO is = 1, nsolid
              velint(is+1)=-((zB-h)*vels(imjkp,is)+(h-zA)*vels(immjkpp,is))/(zB-zA)
            END DO
         ENDIF

      CASE (-2)
   
!===========================================
!====   interpolazione lineare sin. con ====
!====   i valori nei nodi W-WW		====
!===========================================

         h=SQRT((cx(i)-nsx)**2+(cz(k)-nsz)**2)
         zA=SQRT((cx(i-1)-nsx)**2+(cz(k)-nsz)**2)
         IF (h <= zA) THEN
            velint(1)=-h/zA*velg(imjk)
            DO is= 1, nsolid
              velint(1+is)=-h/zA*vels(imjk,is)
            END DO
         ELSE
            zB=SQRT((cx(i-2)-nsx)**2+(cz(k)-nsz)**2)
            velint(1)=-((zB-h)*velg(imjk)+(h-zA)*velg(immjk))/(zB-zA)
            DO is= 1, nsolid
              velint(1+is)=-((zB-h)*vels(imjk,is)+(h-zA)*vels(immjk,is))/(zB-zA)
            END DO
         ENDIF

      CASE (-1)
   
!===========================================
!====   interpolazione bilin. con i	====
!====   valori nei nodi W-TW-T		====
!===========================================

         alpha=(cx(i-1)-nsx)/(cx(i-1)-cx(i))
         beta=(cz(k+1)-nsz)/(cz(k+1)-cz(k))
         velint(1)=-(alpha*(1-beta)*velg(ijkp)+(1-alpha)*(1-beta)*velg(imjkp)+ &
                    (1-alpha)*beta*velg(imjk))/(alpha*beta)
         DO is = 1, nsolid
           velint(1+is)=-(alpha*(1-beta)*vels(ijkp,is)+(1-alpha)*(1-beta)*vels(imjkp,is)+ &
                      (1-alpha)*beta*vels(imjk,is))/(alpha*beta)
         END DO
!
      CASE (0)
   
!===========================================
!====   interpolazione lineare con i	====
!====   valori nei nodi T-TT		====
!===========================================
   
         h=SQRT((cx(i)-nsx)**2+(cz(k)-nsz)**2)
         zA=SQRT((cx(i)-nsx)**2+(cz(k+1)-nsz)**2)
         IF (h <= zA) THEN
            velint(1)=-h/zA*velg(ijkp)
            DO is = 1, nsolid
              velint(1+is)=-h/zA*vels(ijkp,is)
            END DO
         ELSE
            zB=SQRT((cx(i)-nsx)**2+(cz(k+2)-nsz)**2)
            velint(1)=-((zB-h)*velg(ijkp)+(h-zA)*velg(ijkpp))/(zB-zA)
            DO is = 1, nsolid
              velint(1+is)=-((zB-h)*vels(ijkp,is)+(h-zA)*vels(ijkpp,is))/(zB-zA)
            END DO

! ... Quadratic interpolation
!            velint=- ((vel(ijkp)*zB-vel(ijkpp)*zA)/(zA*zB*(zA-zB))*h**2   &
!                    +(vel(ijkpp)*zA**2-vel(ijkp)*zB**2)/(zA*zB*(zA-zB))*h)
         ENDIF

      CASE (1)
   
!===========================================
!====   interpolazione bilin. con i	====
!====   valori nei nodi E-TE-T		====
!===========================================
   
         alpha=(cx(i+1)-nsx)/(cx(i+1)-cx(i))
         beta=(cz(k+1)-nsz)/(cz(k+1)-cz(k))
         velint(1)=-(alpha*(1-beta)*velg(ijkp)+(1-alpha)*(1-beta)*velg(ipjkp)+ &
                 (1-alpha)*beta*velg(ipjk))/(alpha*beta)
         DO is = 1, nsolid
           velint(1+is)=-(alpha*(1-beta)*vels(ijkp,is)+(1-alpha)*(1-beta)*vels(ipjkp,is)+ &
                   (1-alpha)*beta*vels(ipjk,is))/(alpha*beta)
         END DO
         
      CASE (2)
   
!===========================================
!====   interpolazione lineare con i	====
!====   valori nei nodi E-EE		====
!===========================================
   
         h=SQRT((cx(i)-nsx)**2+(cz(k)-nsz)**2)
         zA=SQRT((cx(i+1)-nsx)**2+(cz(k)-nsz)**2)
         IF (h <= zA) THEN
            velint(1)=-h/zA*velg(ipjk)
            DO is = 1, nsolid
              velint(1+is)=-h/zA*vels(ipjk,is)
            END DO
         ELSE
            zB=SQRT((cx(i+2)-nsx)**2+(cz(k)-nsz)**2)
            velint(1)=-((zB-h)*velg(ipjk)+(h-zA)*velg(ippjk))/(zB-zA)
            DO is = 1, nsolid
              velint(1+is)=-((zB-h)*vels(ipjk,is)+(h-zA)*vels(ippjk,is))/(zB-zA)
            END DO
         ENDIF
         
      CASE (3)
   
!===========================================
!====   interpolazione lineare con i	====
!====   valori nei nodi TE-TTEE		====
!===========================================
   
         h=SQRT((cx(i)-nsx)**2+(cz(k)-nsz)**2)
         zA=SQRT((cx(i+1)-nsx)**2+(cz(k+1)-nsz)**2)
         IF (h <= zA) THEN
            velint(1)=-h/zA*velg(ipjkp)
            DO is = 1, nsolid
              velint(1+is)=-h/zA*vels(ipjkp,is)
            END DO
         ELSE
            zB=SQRT((cx(i+2)-nsx)**2+(cz(k+2)-nsz)**2)
            velint(1)=-((zB-h)*velg(ipjkp)+(h-zA)*velg(ippjkpp))/(zB-zA)
            DO is = 1, nsolid
              velint(1+is)=-((zB-h)*vels(ipjkp,is)+(h-zA)*vels(ippjkpp,is))/(zB-zA)
            END DO
         ENDIF
   
      CASE (22)
   
!===========================================
!====   interpolazione lineare      	====
!====   esterna E           		====
!===========================================
   
         IF (nsx < cx(i)) THEN
           h=SQRT((cx(i)-nsx)**2+(cz(k)-nsz)**2)
           zA=SQRT((cx(i+1)-nsx)**2+(cz(k)-nsz)**2)
           velint(1)= + h/zA*velg(ipjk)
           DO is = 1, nsolid
             velint(1+is)= + h/zA*vels(ipjk,is)
           END DO
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
           velint(1)= + h/zA*velg(imjk)
           DO is = 1, nsolid
             velint(1+is)= + h/zA*vels(imjk,is)
           END DO
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
         velint(1)= + h/zA*velg(ijkp)
         DO is = 1, nsolid
           velint(1+is)= + h/zA*vels(ijkp,is)
         END DO

! ... Quadratic interpolation
!         zB=SQRT((cx(i)-nsx)**2+(cz(k+2)-nsz)**2)
!         velint= (vel(ijkp)*zB-vel(ijkpp)*zA)/(zA*zB*(zA-zB))*h**2 &
!                +(vel(ijkpp)*zA**2-vel(ijkp)*zB**2)/(zA*zB*(zA-zB))*h

      CASE DEFAULT
   
!===========================================
!====   nessuna interpolazione  ============
!===========================================
   
         velint(1)=velg(ijk)
         DO is = 1, nsolid
           velint(1+is)=vels(ijk,is)
         END DO

      END SELECT
        
      END FUNCTION velint
!----------------------------------------------------------------------
      SUBROUTINE init_delta
! ... This routine is no longer used ...
! ... Keep it as a reminder ... !!!!!
      IMPLICIT NONE
!
      INTEGER :: delta_i_(0:30)
      INTEGER :: delta_j_(0:30)
      INTEGER :: delta_k_(0:30)
      INTEGER :: index_q_(0:30)
      INTEGER :: index_qq_(0:30)
!
! ... initialize
!
              delta_i_(:) = 0
              delta_j_(:) = 0
              delta_k_(:) = 0
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

              delta_i_(0) = 0
              delta_j_(0) = 0
              delta_k_(0) = 1
!
              delta_i_(1) = -1
              delta_j_(1) = 0
              delta_k_(1) = 0
!
              delta_i_(2) = -1
              delta_j_(2) = -1
              delta_k_(2) = 0
!
              delta_i_(3) = 0
              delta_j_(3) = -1
              delta_k_(3) = 0
!
              delta_i_(4) = 1
              delta_j_(4) = -1
              delta_k_(4) = 0
!
              delta_i_(5) = 1
              delta_j_(5) = 0
              delta_k_(5) = 0
!
              delta_i_(6) = 1
              delta_j_(6) = 1
              delta_k_(6) = 0
!
              delta_i_(7) = 0
              delta_j_(7) = 1
              delta_k_(7) = 0
!
              delta_i_(8) = -1
              delta_j_(8) = 1
              delta_k_(8) = 0
! ... external NE
              delta_i_(26) = 1
              delta_j_(26) = 1
              delta_k_(26) = 0
! ... external E
              delta_i_(25) = 1
              delta_j_(25) = 0
              delta_k_(25) = 0
! ... external SE
              delta_i_(24) = 1
              delta_j_(24) = -1
              delta_k_(24) = 0
! ... external N
              delta_i_(27) = 0
              delta_j_(27) = 1
              delta_k_(27) = 0
! ... external TOP
              delta_i_(20) = 0
              delta_j_(20) = 0
              delta_k_(20) = 1
! ... external S
              delta_i_(23) = 0
              delta_j_(23) = -1
              delta_k_(23) = 0
! ... external NW
              delta_i_(28) = -1
              delta_j_(28) = 1
              delta_k_(28) = 0
! ... external W
              delta_i_(21) = -1
              delta_j_(21) = 0
              delta_k_(21) = 0
! ... external SW
              delta_i_(22) = -1
              delta_j_(22) = -1
              delta_k_(22) = 0
!
      RETURN
      END SUBROUTINE init_delta
!----------------------------------------------------------------------
      FUNCTION velint3d(fpt, velg, vels, ijk, cx, cy, cz, index_q)
!
! ... Interpolate velocities on a forcing point to get no-slip
! ... conditions on a solid immersed boundary 

      USE dimensions, ONLY: nsolid, max_nsolid
      USE set_indexes
      USE immersed_boundaries, ONLY: forcing_point
      IMPLICIT NONE

      REAL*8, DIMENSION(max_nsolid+1) :: velint3d
      TYPE(forcing_point), INTENT(IN) :: fpt
      REAL*8, DIMENSION(:), INTENT(IN) :: velg
      REAL*8, DIMENSION(:,:), INTENT(IN) :: vels
      REAL*8, DIMENSION(:), INTENT(IN) :: cx, cy, cz
      INTEGER, INTENT(IN) :: ijk
      INTEGER, INTENT(OUT) :: index_q

      INTEGER :: i, j, k, is
      REAL*8 :: nsx, nsy, nsz
      INTEGER :: interp, delta_i, delta_j, delta_k
      INTEGER :: index_qq
      LOGICAL :: diagonal
      REAL*8 :: velint_1, velint_2

      REAL*8 :: h          !distance between the (i,j,k)-node and boundary
      REAL*8 :: zA         !distance between the boundary and the first 
                           !external node
      REAL*8 :: zB         !distance between the boundary and the second 
                           !external node
      REAL*8 :: hza        !ratio between the distances
      REAL*8 :: alpha,beta !coefficients for bilinear interpolation
!
      velint3d = 0.D0

      interp = fpt%int
      i      = fpt%i
      j      = fpt%j
      k      = fpt%k
      nsx    = fpt%nsl%x
      nsy    = fpt%nsl%y
      nsz    = fpt%nsl%z
!
      delta_i = fpt%delta_i
      delta_j = fpt%delta_j
      delta_k = fpt%delta_k
      index_q = fpt%index_q
      index_qq = fpt%index_qq
!
      diagonal = ( (MOD(interp,2) == 0) .AND. (interp/=0) )

      h  = SQRT( (cx(i)-nsx)**2 + (cy(j)-nsy)**2 + (cz(k)-nsz)**2 )
      
      zA = SQRT( (cx(i+delta_i)-nsx)**2 + (cy(j+delta_j)-nsy)**2 + &
                 (cz(k+delta_k)-nsz)**2 )
      hzA = h/zA

      IF (interp >= 20) THEN
         velint3d(1) = hzA*velg(index_q)         
         DO is = 1, nsolid
           velint3d(1+is) = hzA*vels(index_q,is)         
         END DO
      ELSEIF (h <= zA .OR. diagonal) THEN
         velint3d(1) = -hzA*velg(index_q)
         DO is = 1, nsolid
           velint3d(1+is) = -hzA*vels(index_q,is)
         END DO
      ELSE 
         zB=SQRT( (cx(i+2*delta_i)-nsx)**2 + (cy(j+2*delta_j)-nsy)**2 +  &
                  (cz(k+2*delta_k)-nsz)**2 )
         !
         velint_2 = -((zB-h)*velg(index_q)+(h-zA)*velg(index_qq))/(zB-zA)
         velint_1 = -hzA*velg(index_q)
         IF (DABS(velint_1) < DABS(velint_2)) THEN
           velint3d(1) = velint_1
         ELSE
           velint3d(1) = velint_2
         END IF
         !
         DO is = 1, nsolid
           velint_2 = -((zB-h)*vels(index_q,is)+(h-zA)*vels(index_qq,is))/(zB-zA)
           velint_1 = -hzA*vels(index_q,is)
           IF (DABS(velint_1) < DABS(velint_2)) THEN
             velint3d(1+is) = velint_1
           ELSE
             velint3d(1+is) = velint_2
           END IF
         END DO
      ENDIF
        
      RETURN
      END FUNCTION velint3d
!----------------------------------------------------------------------
      END MODULE interpolate_fields
!----------------------------------------------------------------------
