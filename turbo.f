!----------------------------------------------------------------------
      MODULE turbulence
!----------------------------------------------------------------------
!
      USE gas_solid_velocity, ONLY: ug,wg 
      USE grid, ONLY: fl_l, myijk,  ncint, ncdom, data_exchange
      USE environment, ONLY: timing, cpclock
      USE indijk_module, ONLY: ip0_jp0_kp0_

      IMPLICIT NONE
!
      REAL*8, DIMENSION(:),   ALLOCATABLE :: mugt   ! gas turbulent viscosity
      REAL*8, DIMENSION(:),   ALLOCATABLE :: kapgt  ! gas turbulent conductivity
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: must   ! solid turbulent viscosity


      REAL*8, DIMENSION(:), ALLOCATABLE :: smag_factor, smag_coeff
      REAL*8, DIMENSION(:), ALLOCATABLE :: smag, scoeff
!
      INTEGER :: iturb, iss, modturbo
      REAL*8 :: cmut                       ! Smagorinsky constant
      REAL*8 :: pranumt                    ! Prandtl number
!
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE bounds_turbo
      USE dimensions
      IMPLICIT NONE
!
      ALLOCATE(smag_factor(ntot))
      ALLOCATE(smag_coeff(ntot))
      smag_factor = 0.0D0
      smag_coeff = 0.0D0
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE local_bounds_turbo
      USE dimensions
      IMPLICIT NONE
!
      ALLOCATE(smag(ncint))
      ALLOCATE(scoeff(ncint))
      ALLOCATE(mugt(ncdom))
      ALLOCATE(must(nsolid, ncdom))
      ALLOCATE(kapgt(ncdom))
      smag = 0.0D0
      kapgt = 0.0D0
      mugt  = 0.0D0
      must  = 0.0D0
      scoeff = 0.0D0
      RETURN
      END SUBROUTINE

! ... MODIFICARE_X3D ( fino fine file )

!----------------------------------------------------------------------

       SUBROUTINE turbulence_setup( zrough )

! ... sets the Smagorinsky turbulence length scale
!
        USE dimensions, ONLY: nr, nz, no, nx, ny
        USE grid, ONLY: iob
        USE grid, ONLY: dz, dr, zb, rb, dx, dy
        USE roughness_module, ONLY: roughness
        USE control_flags, ONLY: job_type
!
        TYPE (roughness), INTENT(IN) :: zrough

        INTEGER :: i, j, k, i1, i2, j2, n, ijk
        REAL*8  :: sl, sgsl0, zsgsl, zrou, sgslb
        REAL*8 :: delt

        INTEGER, DIMENSION(:), ALLOCATABLE :: jt, it
        REAL*8,  DIMENSION(:), ALLOCATABLE :: zbt
!
!
        ALLOCATE( jt( MAX(nr,1) ), it( MAX(nr,1) ), zbt( MAX(nr,1) ) )

        IF (iturb .EQ. 2) THEN

          IF( job_type == '2D' ) THEN

            it = 0; jt = 0
            zbt = 0.D0
            DO n=1,no
             IF( iob(n)%typ == 3 ) THEN
               it( iob(n)%rlo : iob(n)%rhi ) = 1
               jt(i)  = iob(n)%zhi
               zbt(i) = zb( iob(n)%zhi )
             ENDIF
            END DO

          END IF

        ENDIF
!
        IF( job_type == '2D' ) THEN

          DO j = 2, (nz-1)

            DO i = 2, (nr-1) 

              ijk = i+(j-1)*nr

              delt = DSQRT( dz(j)*dr(i) )

              IF ( iturb == 1 ) THEN

                sl = cmut * delt

              ELSE IF ( iturb == 2 ) THEN
!
                IF( j <= jt(i) ) THEN

                  sl = 0.D0

                ELSE

                  sgsl0 = cmut * delt

                  IF( it(i) == 1 ) THEN

                    zsgsl = zb(j)-zbt(i)-dz(j)*0.5D0
                    IF( zrough%ir == 1 ) THEN
                      zrou = zrough%r(1)
                    ELSE
                      IF( rb(i) <= zrough%roucha ) THEN
                        zrou = zrough%r(1)
                      ELSE
                        zrou = zrough%r(2)
                      ENDIF
                    ENDIF
                    sgslb=0.4D0*(zsgsl+zrou)
!
                    ! ... use the smallest between the smag length and the roughness length

                    sl = dmin1(sgsl0,sgslb)
  
                  ELSE

                    sl = sgsl0

                  ENDIF

                ENDIF

              ENDIF
!
              ! ... Squared turbulence length scale is used into Smagorinsky model

              smag_factor(ijk) = sl**2
!
            END DO

          END DO

        ELSE IF( job_type == '3D' ) THEN

          DO k = 1, nz

            DO j = 1, ny

              DO i = 1, nx

                ijk = i + (j-1)*nx + (k-1)*nx*ny
                delt = ( dz(k)*dy(j)*dx(i) )**(1.0d0/3.0d0)

                IF ( iturb == 1 ) THEN

                  sl = cmut * delt

                ELSE IF ( iturb == 2 ) THEN

                  CALL error( ' turbo_setup ' , ' undefined 3D roughness ', 1 )

                END IF

                ! ... Squared turbulence length scale is used into Smagorinsky model

                smag_factor(ijk) = sl**2

              END DO

            END DO

          END DO

        END IF

        DEALLOCATE(jt, it, zbt)

       END  SUBROUTINE 

!----------------------------------------------------------------------
      SUBROUTINE sgsg

! ... Sub Grid Stress (Gas)
! ... compute dynamic C , dynamic turbulent viscosity and 
! ... turbulent conductivity in the center of the cell.

      USE boundary_conditions, ONLY: fboundary
      USE dimensions, ONLY: nr
      USE eos_gas, ONLY: cg
      USE gas_solid_density, ONLY: rog
      USE set_indexes, ONLY: subscr
      USE grid, ONLY: dr, dz 
      IMPLICIT NONE
!
      REAL*8 ,DIMENSION(:), ALLOCATABLE :: modsr, sr1, sr2, sr12 
      REAL*8 ,DIMENSION(:), ALLOCATABLE :: p1, p2, p12 
      REAL*8, DIMENSION(:), ALLOCATABLE :: fug, fwg, fu2g, fw2g, fuwg
!
      REAL*8 :: fsr1, fsr2, fsr12, fmodsr
      REAL*8 :: fp1, fp2, fp12
      REAL*8 :: l1, l2, l12, m1, m2, m12, l1d, l2d
      REAL*8 :: cdyn 
      REAL*8 :: num, den
      REAL*8 :: delt
      INTEGER  :: i, j, ij, imesh  
!
      ALLOCATE(modsr(ncdom));     modsr = 0.0D0
      ALLOCATE(sr1(ncdom))  ;     sr1 = 0.0D0
      ALLOCATE(sr2(ncdom))  ;     sr2 = 0.0D0
      ALLOCATE(sr12(ncdom)) ;     sr12 = 0.0D0
      ALLOCATE(p1(ncdom))   ;     p1 = 0.0D0
      ALLOCATE(p2(ncdom))   ;     p2 = 0.0D0
      ALLOCATE(p12(ncdom))  ;     p12 = 0.0D0
!
      CALL data_exchange(ug)
      CALL data_exchange(wg)
!
      DO ij = 1, ncint
        IF(fl_l(ij).EQ.1) THEN
          CALL subscr(ij)
          CALL strain(ug, wg, ij, modsr(ij), sr1(ij), sr2(ij), sr12(ij), p1(ij), p2(ij), p12(ij))
        END IF
      END DO
!
      IF (modturbo .EQ. 2) THEN
        ALLOCATE(fug(ncdom))  ;     fug = 0.0D0
        ALLOCATE(fwg(ncdom))  ;     fwg = 0.0D0
        ALLOCATE(fu2g(ncdom)) ;     fu2g = 0.0D0
        ALLOCATE(fw2g(ncdom)) ;     fw2g = 0.0D0
        ALLOCATE(fuwg(ncdom)) ;     fuwg = 0.0D0
!
! ... compute filtered velocities
!
        CALL vel_hat(fug, fwg, fu2g, fw2g, fuwg)        
!
        CALL fboundary(fug, fwg)
!
        CALL data_exchange(p1)
        CALL data_exchange(p2)
        CALL data_exchange(p12)
        CALL data_exchange(fug)
        CALL data_exchange(fwg)
        CALL data_exchange(fu2g)
        CALL data_exchange(fw2g)
        CALL data_exchange(fuwg)
      END IF
!
       DO ij = 1, ncint
        IF(fl_l(ij).EQ.1) THEN
          CALL subscr(ij)
          imesh = myijk( ip0_jp0_kp0_, ij)
          j = ( imesh - 1 ) / nr + 1
          i = MOD( ( imesh - 1 ), nr) + 1
!
! ... Dynamic computation of the Smagorinsky length scale (Germano et al., 1990)
!
          IF(modturbo.EQ.2) THEN
             delt = DSQRT(dz(j)*dr(i))
             CALL strain(fug, fwg, ij, fmodsr, fsr1, fsr2, fsr12)
!
             fp1  = filter(p1,ij)
             fp2  = filter(p2,ij)
             fp12 = filter(p12,ij)   
!         
             l1  = fug(ij)**2 - fu2g(ij)
             l2  = fwg(ij)**2 - fw2g(ij)
             l12 = fug(ij)*fwg(ij) - fuwg(ij)
!
!             l1d = 0.5D0*l1 - 0.5D0*l2
!             l2d = 0.5D0*l2 - 0.5D0*l1  
!
             m1  = (4.D0*fmodsr*fsr1 - fp1)
             m2  = (4.D0*fmodsr*fsr2 - fp2)
             m12 = (4.D0*fmodsr*fsr12 - fp12)
!              
! ... min-square calculation
!
             num = l1*m1 + l2*m2 + 2*l12*m12
             den = delt**2.D0*( m1**2.D0 + m2**2.D0 + 2.D0*m12**2.D0)
!
             IF(den.EQ.0.0D0) THEN
              cdyn = 0.0D0
              ELSE 
              cdyn = 0.5D0*num/den 
             END IF
             smag(ij) = cdyn * (delt**2.D0)
          END IF
!
! ... Smagorinsky turbulent viscosity
!
          mugt(ij) = rog(ij) * smag(ij) * modsr(ij)
          IF(mugt(ij).LT.0.0D0) mugt(ij) = 0.0D0
!
          scoeff(ij) = cdyn
!-turb-heat*********************
          pranumt=0.5
          kapgt(ij)=mugt(ij)*cg(ij)/pranumt
!-turb-heat**********************

        ELSE
          mugt(ij)  = 0.D0
          kapgt(ij) = 0.D0
        END IF
       END DO
!                          
      DEALLOCATE(modsr)
      DEALLOCATE(sr1)
      DEALLOCATE(sr2)
      DEALLOCATE(sr12)
      DEALLOCATE(p1)
      DEALLOCATE(p2)
      DEALLOCATE(p12)
      IF (modturbo .EQ. 2) THEN
        DEALLOCATE(fug)
        DEALLOCATE(fwg)
        DEALLOCATE(fu2g)
        DEALLOCATE(fw2g)
        DEALLOCATE(fuwg)
      END IF

      RETURN 
      END SUBROUTINE sgsg
!----------------------------------------------------------------------
      SUBROUTINE strain(u, w, ij, modsr, sr1, sr2, sr12, p1, p2, p12 )
!
! ... here computes the components of the strain rate tensor and its module.
! ... and the components of the function pij(=modsr*srij)

      USE dimensions, ONLY: nr    
      USE grid, ONLY:dr, dz, indr, indz,itc,inr 
      USE set_indexes
      IMPLICIT NONE

      REAL*8, INTENT(IN), DIMENSION(:) :: u, w
      INTEGER,INTENT(IN)  :: ij
      REAL*8, INTENT(OUT) :: modsr, sr1, sr2, sr12
      REAL*8, OPTIONAL, INTENT(OUT) :: p1, p2, p12
!
      REAL*8 :: w1, w2, w3, u1, u2, u3, um1, um2, wm1, wm2 
      REAL*8 :: drp, dzm, drm, dzp, indrp, indzm, indrm, indzp
      REAL*8 :: d33
      INTEGER :: i, j, imesh 
!
        imesh = myijk( ip0_jp0_kp0_, ij)
        j = ( imesh - 1 ) / nr + 1
        i = MOD( ( imesh - 1 ), nr) + 1

        drp=dr(i)+dr(i+1)
        drm=dr(i)+dr(i-1)
        dzp=dz(j)+dz(j+1)
        dzm=dz(j)+dz(j-1)
        indrp=1.D0/drp
        indrm=1.D0/drm
        indzp=1.D0/dzp
        indzm=1.D0/dzm
!
! ... Compute velocity gradients in the center of the cell.
! ... Cross terms (non-diagonal) in the strain tensor are obtained
! ... by interpolating the values found on the staggered grids.
!
        sr1 = (u(ij)-u(imj))*indr(i)
        sr2 = (w(ij)-w(ijm))*indz(j)
!
        
! ... extra-term for cylindrical coordinates
!
        IF(itc.EQ.1) THEN
          d33 = (u(ij)+u(imj))*inr(i)
        ELSE
          d33 = 0.D0
        END IF
!
        u3=0.5D0*(u(ijp)+u(imjp))
        u2=0.5D0*(u(ij)+u(imj))
        u1=0.5D0*(u(ijm)+u(imjm))
        w3=0.5D0*(w(ipj)+w(ipjm))
        w2=0.5D0*(w(ij)+w(ijm))
        w1=0.5D0*(w(imj)+w(imjm))
        um2=u2+(u3-u2)*dz(j)*indzp
        um1=u1+(u2-u1)*dz(j-1)*indzm
        wm2=w2+(w3-w2)*dr(i)*indrp
        wm1=w1+(w2-w1)*dr(i-1)*indrm

        sr12 =((um2-um1)*indz(j)+(wm2-wm1)*indr(i))*0.5D0
!        
        modsr = DSQRT(2 * (sr1**2 + sr2**2 + 2*sr12**2) )

        IF (PRESENT(p1)) THEN 
          p1  =  modsr * sr1
          p2  =  modsr * sr2
          p12 =  modsr * sr12
        END IF
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE vel_hat(fug, fwg, fu2g, fw2g, fuwg) 
!
!... the filtered components of the velocity are calculated...
!... before filtering the product of the components of the velocity we 
!... have defined the function u_per_w(ij)(next subroutine)...

      USE set_indexes, ONLY: subscr
      IMPLICIT NONE
 
      REAL*8, DIMENSION(:), ALLOCATABLE :: uwg, u2g, w2g
      REAL*8, INTENT(OUT), DIMENSION(:) :: fug, fwg, fu2g, fw2g, fuwg
      INTEGER :: i, j, ij
    
      ALLOCATE(uwg(ncdom)) 
      ALLOCATE(u2g(ncdom)) 
      ALLOCATE(w2g(ncdom)) 
      uwg = 0.0D0
      u2g = 0.0D0
      w2g = 0.0D0
!
      CALL u_per_w(uwg) 
      u2g = ug * ug
      w2g = wg * wg
!
      CALL data_exchange(uwg)
!
       DO ij = 1, ncint
         IF (fl_l(ij) .EQ. 1) THEN 
          CALL subscr(ij)
          fug(ij)  = filter(ug,ij)
          fwg(ij)  = filter(wg,ij)
          fu2g(ij) = filter(u2g,ij)
          fw2g(ij) = filter(w2g,ij)
          fuwg(ij) = filter(uwg,ij)     
         END IF 
       END DO 
!
      DEALLOCATE(uwg)
      DEALLOCATE(u2g)
      DEALLOCATE(w2g)
! 
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE u_per_w(uwg)
      USE set_indexes
! 
      IMPLICIT NONE
      INTEGER :: ij
      REAL*8, INTENT(OUT) :: uwg(:)
!
      DO ij = 1, ncint
        IF(fl_l(ij) .EQ. 1) THEN
          CALL subscr(ij)
          uwg(ij)=(0.5D0*(ug(ij)+ug(imj)))*(0.5D0*(wg(ij)+wg(ijm)))
        ELSE
          uwg(ij) = ug(ij)*wg(ij)
        END IF
      END DO   
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
       FUNCTION filter(t,ij)
       USE set_indexes, ONLY: imj, ipj, ijm, ijp
!
       REAL*8,INTENT(IN),DIMENSION(:):: t
       INTEGER,INTENT(IN):: ij
       REAL*8 :: filter
!
       filter = (0.125)*(t(imj)+t(ipj)+t(ijm)+t(ijp))+0.5D0*t(ij)
!
       END FUNCTION filter
!----------------------------------------------------------------------
      SUBROUTINE sgss
!
      USE dimensions
      USE gas_solid_density, ONLY: rlk
      USE gas_solid_velocity, ONLY: us, ws
      USE grid, ONLY: itc, dz, dr, r, rb, indz, indr, inr, inrb
      USE particles_constants, ONLY: rl, inrl, dk, cmus
      USE set_indexes
      IMPLICIT NONE
!
      REAL*8 :: w1, w2, w3, u1, u2, u3, wm1, wm2, um1, um2
      REAL*8 :: d12, d11, d22, d33
      REAL*8 :: drp, dzm, drm, dzp, indrp, indzm, indrm, indzp 
      INTEGER :: ij
      INTEGER :: i, j, is, imesh
!
      CALL data_exchange(us)
      CALL data_exchange(ws)
!
      d33=0.D0
      DO ij = 1, ncint
        imesh = myijk( ip0_jp0_kp0_, ij)
        IF(fl_l(ij).EQ.1) THEN
         CALL subscr(ij)
         j = ( imesh - 1 ) / nr + 1
         i = MOD( ( imesh - 1 ), nr) + 1
!
         drp=dr(i)+dr(i+1)
         drm=dr(i)+dr(i-1)
         dzp=dz(j)+dz(j+1)
         dzm=dz(j)+dz(j-1)
!
         indrp=1.D0/drp
         indrm=1.D0/drm
         indzp=1.D0/dzp
         indzm=1.D0/dzm
!         
         DO is=1,nsolid
           d11=2.D0*(us(is,ij)-us(is,imj))*indr(i)
           d22=2.D0*(ws(is,ij)-ws(is,ijm))*indz(j)
           IF(itc.EQ.1)d33=(us(is,ij)+us(is,imj))*inr(i)
           u3=0.5D0*(us(is,ijp)+us(is,imjp))
           u2=0.5D0*(us(is,ij)+us(is,imj))
           u1=0.5D0*(us(is,ijm)+us(is,imjm))
           w3=0.5D0*(ws(is,ipj)+ws(is,ipjm))
           w2=0.5D0*(ws(is,ij)+ws(is,ijm))
           w1=0.5D0*(ws(is,imj)+ws(is,imjm))
           um2=u2+(u3-u2)*dz(j)*indzp
           um1=u1+(u2-u1)*dz(j-1)*indzm
           wm2=w2+(w3-w2)*dr(i)*indrp
           wm1=w1+(w2-w1)*dr(i-1)*indrm
           d12=(um2-um1)*indz(j)+(wm2-wm1)*indr(i)
           must(is,ij) = cmus(is)*0.7071D0*rl(is)*dk(is)**2 *            &
     &                 10.D0**(3.98D0*rlk(is,ij)*inrl(is)-1.69D0) *   &
     &                 DSQRT(d11*d11+d22*d22+d33*d33+2.0D0*d12*d12)
         END DO
        END IF
      END DO
!
      RETURN
      END SUBROUTINE
!---------------------------------------------------------------------
       END MODULE turbulence
!---------------------------------------------------------------------
