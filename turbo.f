!----------------------------------------------------------------------
      MODULE turbulence
!----------------------------------------------------------------------
!
      USE gas_solid_velocity, ONLY: ug,vg 
      USE grid, ONLY: fl_l, myij,  nij_l, nijx_l, data_exchange
      USE environment, ONLY: timing, cpclock
      IMPLICIT NONE
!
      REAL*8, DIMENSION(:),   ALLOCATABLE :: mugt   !turbulent viscosity (gas)
      REAL*8, DIMENSION(:,:),   ALLOCATABLE :: mus   !turbulent viscosity (p.)
      REAL*8, DIMENSION(:),   ALLOCATABLE :: kapgt  !turb. conductivity (gas)

      REAL*8, DIMENSION(:,:),   ALLOCATABLE :: mus_g
      REAL*8, DIMENSION(:), ALLOCATABLE :: smagl_g, zbt, scoeff_g
      REAL*8, DIMENSION(:), ALLOCATABLE :: smagl, scoeff
      REAL*8, DIMENSION(:), ALLOCATABLE :: cdy_n 
      REAL*8, DIMENSION(:), ALLOCATABLE :: delta 
      INTEGER, DIMENSION(:), ALLOCATABLE :: jt, it
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
      ALLOCATE(smagl_g(ndi*ndj))
      ALLOCATE(mus_g(ncl,ndi*ndj))
      ALLOCATE(jt(ndi), it(ndi), zbt(ndi))
      ALLOCATE(scoeff_g(ndi*ndj))
      smagl_g = 0.0D0
      mus_g = 0.0D0
      scoeff_g = 0.0D0
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE local_bounds_turbo
      USE dimensions
      IMPLICIT NONE
!
      ALLOCATE(smagl(nij_l))
      ALLOCATE(scoeff(nij_l))
      ALLOCATE(delta(nij_l))
      ALLOCATE(mugt(nijx_l))
      ALLOCATE(mus(ncl, nijx_l))
      ALLOCATE(kapgt(nijx_l))
      ALLOCATE(cdy_n(nijx_l))
      smagl = 0.0D0
      kapgt = 0.0D0
      mugt  = 0.0D0
      mus  = 0.0D0
      scoeff = 0.0D0
      cdy_n = 0.0D0
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
       SUBROUTINE turbulence_setup( ir, zrough, roucha)
!
        USE grid, ONLY: iob, nso, no
        USE grid, ONLY: ib2, ib1, jb1, dz, dr, zb, rb
!
        INTEGER, INTENT(IN) :: ir
        REAL*8, INTENT(IN) :: zrough(:), roucha
        INTEGER :: i, j, i1, i2, j2, n, ij
        REAL*8  :: sgsl0, zsgsl, zrou, sgslb, delt
        DO i=2,ib1
           it(i)=0
           DO n=1,no
             IF(nso(n).EQ.3.OR.nso(n).EQ.5) THEN
               i1=iob(1,n)
               i2=iob(2,n)
               j2=iob(4,n)
               IF(i.GE.i1.AND.i.LE.i2) THEN
                 IF(nso(n).EQ.3) it(i)=1
                 jt(i)=j2
                 zbt(i)=zb(j2)
               ENDIF
             ENDIF
           END DO
         END DO
!
         DO j=2,jb1
         DO i=2,ib1
           ij=i+(j-1)*ib2
           IF(j.LE.jt(i)) THEN
             smagl_g(ij)=0.D0
           ELSE
             delt=DSQRT(dz(j)*dr(i))
!             WRITE(*,*)'ij, delt', ij, delt
             IF(iturb.EQ.1) THEN
               smagl_g(ij)=cmut*delt
             ELSE
               sgsl0=cmut*delt
               IF(it(i).EQ.1) THEN
                 zsgsl=zb(j)-zbt(i)-dz(j)*0.5D0
                 IF(ir.EQ.1) THEN
                   zrou=zrough(1)
                 ELSE
                   IF(rb(i).LE.roucha) THEN
                     zrou=zrough(1)
                   ELSE
                     zrou=zrough(2)
                   ENDIF
                 ENDIF
                 sgslb=0.4D0*(zsgsl+zrou)
                 smagl_g(ij)=dmin1(sgsl0,sgslb)
               ELSE
                 smagl_g(ij)=sgsl0
               ENDIF
             ENDIF
           ENDIF
         END DO
         END DO

       END  SUBROUTINE 
!----------------------------------------------------------------------
      SUBROUTINE sgsgdyn 

!...compute dynamic C , dynamic turbulent viscosity and 
!...turbulent conductivity in the center of the cell.

      USE dimensions, ONLY: ndi
      USE eos_gas, ONLY: cg
      USE gas_solid_density, ONLY: rog
      USE set_indexes
      USE grid, ONLY: dr, dz 
!
      IMPLICIT NONE
      REAL*8 ,DIMENSION(:), ALLOCATABLE :: modsr, sr1, sr2, sr12 
      REAL*8 ,DIMENSION(:), ALLOCATABLE :: p1, p2, p12 
      REAL*8, DIMENSION(:), ALLOCATABLE :: fug, fvg, fu2g, fv2g, fuvg
      REAL*8 :: fsr1, fsr2, fsr12, fmodsr
      REAL*8 :: fp1, fp2, fp12
      REAL*8 :: l1, l2, l12, m1, m2, m12, l1d, l2d
      REAL*8 :: cdy 
      REAL*8 :: num, den
      INTEGER  :: i, j, ij, ij_g  
!
      ALLOCATE(modsr(nijx_l))
      ALLOCATE(sr1(nijx_l))
      ALLOCATE(sr2(nijx_l))
      ALLOCATE(sr12(nijx_l))
      ALLOCATE(p1(nijx_l))
      ALLOCATE(p2(nijx_l))
      ALLOCATE(p12(nijx_l))
      ALLOCATE(fug(nijx_l))
      ALLOCATE(fvg(nijx_l))
      ALLOCATE(fu2g(nijx_l))
      ALLOCATE(fv2g(nijx_l))
      ALLOCATE(fuvg(nijx_l))
      modsr = 0.0D0
      sr1 = 0.0D0
      sr2 = 0.0D0
      sr12 = 0.0D0
      p1 = 0.0D0
      p2 = 0.0D0
      p12 = 0.0D0
      fug = 0.0D0
      fvg = 0.0D0
      fu2g = 0.0D0
      fv2g = 0.0D0
      fuvg = 0.0D0
!
      CALL data_exchange(ug)
      CALL data_exchange(vg)
!
      CALL strain(modsr, sr1, sr2, sr12, p1, p2, p12)
!
      CALL vel_hat(fug, fvg, fu2g, fv2g, fuvg)        
!
       CALL data_exchange(p1)
       CALL data_exchange(p2)
       CALL data_exchange(p12)
       CALL data_exchange(fug)
       CALL data_exchange(fvg)
       CALL data_exchange(fu2g)
       CALL data_exchange(fv2g)
       CALL data_exchange(fuvg)
!
       DO ij = 1, nij_l
        ij_g = myij(0,0,ij)
        IF(fl_l(ij).EQ.1) THEN
        CALL subscl(ij)
         j = ( ij_g - 1 ) / ndi + 1
         i = MOD( ( ij_g - 1 ), ndi) + 1
         delta(ij) = DSQRT(dz(j)*dr(i))
!
!...here it is possible to choose between the classic or dynamic Smagorinsky...
!
         IF(modturbo.EQ.1) THEN
          cdy_n(ij)=cmut**2
         ELSE
             CALL strain_hat(fsr1, fsr2, fsr12, fmodsr,
     &                        fug, fvg, ij)
!
             fp1  = filter(p1,ij)
             fp2  = filter(p2,ij)
             fp12 = filter(p12,ij)   
!         
             l1  = fug(ij)**2 - fu2g(ij)
             l2  = fvg(ij)**2 - fv2g(ij)
             l12 = fug(ij)*fvg(ij) - fuvg(ij)
!
!             l1d = 0.5D0*l1 - 0.5D0*l2
!             l2d = 0.5D0*l2 - 0.5D0*l1  
!
             m1  = (4*fmodsr*fsr1 - fp1)
             m2  = (4*fmodsr*fsr2 - fp2)
             m12 = (4*fmodsr*fsr12 - fp12)
!              
             num = l1*m1 + l2*m2 + 2*l12*m12
             den = delta(ij)**2*( m1**2 + m2**2 + 2*m12**2)
             IF(den.EQ.0.0D0) THEN
              cdy_n(ij) = 0.0D0
              ELSE 
              cdy = 0.5D0*num/den 
              cdy_n(ij) = cdy
!   
             END IF
          END IF
          smagl(ij) =cdy_n(ij)* (delta(ij))**2
          mugt(ij)  = rog(ij)*smagl(ij)*modsr(ij)
            IF(mugt(ij).LT.0.0D0) mugt(ij) = 0.0D0
            scoeff(ij) = cdy_n(ij)
!-turb-heat*********************
            pranumt=0.5
            kapgt(ij)=mugt(ij)*cg(ij)/pranumt
!-turb-heat**********************

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
      DEALLOCATE(fug)
      DEALLOCATE(fvg)
      DEALLOCATE(fu2g)
      DEALLOCATE(fv2g)
      DEALLOCATE(fuvg)

      RETURN 
      END SUBROUTINE sgsgdyn
!----------------------------------------------------------------------
      SUBROUTINE strain(modsr, sr1, sr2, sr12, p1, p2, p12 )
!
!....here computes the components of the strain rate tensor and its module.
!...computes also the components of the function pij(=modsr*srij) that in the 
!... previous subroutine needs to be filtered...

      USE dimensions    
      USE grid, ONLY:dr, dz, indr, indz,itc,inr 
      USE set_indexes
      IMPLICIT NONE

      REAL*8 :: v1, v2, v3, u1, u2, u3, um1, um2, vm1, vm2 
      REAL*8 :: drp, dzm, drm, dzp, indrp, indzm, indrm, indzp
      REAL*8 :: d33
      REAL*8, INTENT(OUT), DIMENSION(:) :: modsr, sr1, sr2, sr12
      REAL*8, INTENT(OUT), DIMENSION(:) :: p1, p2, p12
      INTEGER :: i, j, ij, ij_g 
!
      DO ij = 1, nij_l
       ij_g = myij(0,0,ij)
       IF(fl_l(ij).EQ.1) THEN
        CALL subscl(ij)
        j = ( ij_g - 1 ) / ndi + 1
        i = MOD( ( ij_g - 1 ), ndi) + 1

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
        sr1(ij) = (ug(ij)-ug(imj))*indr(i)
        sr2(ij) = (vg(ij)-vg(ijm))*indz(j)
!
        
! ... extra-term for cylindrical coordinates
!
        IF(itc.EQ.1) d33 = (ug(ij)+ug(imj))*inr(i)
!
        u3=0.5D0*(ug(ijp)+ug(imjp))
        u2=0.5D0*(ug(ij)+ug(imj))
        u1=0.5D0*(ug(ijm)+ug(imjm))
        v3=0.5D0*(vg(ipj)+vg(ipjm))
        v2=0.5D0*(vg(ij)+vg(ijm))
        v1=0.5D0*(vg(imj)+vg(imjm))
        um2=u2+(u3-u2)*dz(j)*indzp
        um1=u1+(u2-u1)*dz(j-1)*indzm
        vm2=v2+(v3-v2)*dr(i)*indrp
        vm1=v1+(v2-v1)*dr(i-1)*indrm

        sr12(ij) =((um2-um1)*indz(j)+(vm2-vm1)*indr(i))*0.5D0
!        
        modsr(ij) = DSQRT(2*sr1(ij)**2 + 2*sr2(ij)**2 +4*sr12(ij)**2)

        p1(ij) =  modsr(ij)*sr1(ij)
        p2(ij) =  modsr(ij)*sr2(ij)
        p12(ij)=  modsr(ij)*sr12(ij)
       END IF
      END DO
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE vel_hat(fug, fvg, fu2g, fv2g, fuvg) 
!
!
!... the filtered components of the velocity are calculated...
!... before filtering the product of the components of the velocity we 
!... have defined the function u_per_v(ij)(next subroutine)...

      USE set_indexes
      IMPLICIT NONE
 
      REAL*8,DIMENSION(:),ALLOCATABLE::uvg
      REAL*8,INTENT(OUT),DIMENSION(:):: fug, fvg, fu2g, fv2g, fuvg
      INTEGER :: i, j, ij
    
      ALLOCATE(uvg(nijx_l)) 
      uvg = 0.0D0
      CALL u_per_v(uvg) 
      CALL data_exchange(uvg)
!
       DO ij = 1, nij_l
         IF (fl_l(ij) .EQ. 1) THEN 
          CALL subscl(ij)
          fug(ij) = filter(ug,ij)
          fvg(ij) = filter(vg,ij)
          fu2g(ij)= 0.125 *(ug(imj)*ug(imj) + ug(ipj)*ug(ipj) + 
     &              ug(ijm)*ug(ijm) + ug(ijp)*ug(ijp)) + 
     &              0.5D0*ug(ij)*ug(ij)
          fv2g(ij)= 0.125 *(vg(imj)*vg(imj) + vg(ipj)*vg(ipj) + 
     &              vg(ijm)*vg(ijm) + vg(ijp)*vg(ijp)) + 
     &              0.5D0*vg(ij)*vg(ij)
          fuvg(ij)= filter(uvg,ij)     

         END IF 
       END DO 
!
      DEALLOCATE(uvg)
! 
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE strain_hat(fsr1,fsr2,fsr12,fmodsr,fug,fvg,ij)
!
!...computes the filtered components of the strain rate tensor and the 
!...filtered module using not the filter but the right analitycal formula  
!...in which fug and fvg are used instead of the velocities u and v...

      USE dimensions
      USE grid, ONLY:dr, dz, indr, indz , itc, inr
      USE grid, ONLY: fl
      USE set_indexes
!
      IMPLICIT NONE
      SAVE
      INTEGER, INTENT (IN) :: ij
      REAL*8, INTENT(IN) :: fug(:), fvg(:)
      REAL*8, INTENT(OUT) :: fsr1, fsr2, fsr12, fmodsr
      REAL*8 :: fv1, fv2, fv3, fu1, fu2, fu3, fum1, fum2, fvm1, fvm2
      REAL*8 :: v1, v2, v3, u1, u2, u3, um1,um2, vm1, vm2  
      REAL*8 :: drp, dzm, drm, dzp, indrp, indzm, indrm, indzp       
      REAL*8 :: d33
      INTEGER :: ij_g, i, j
      ij_g = myij(0,0,ij)
      CALL subscl(ij)
      j = ( ij_g - 1 ) / ndi + 1
      i = MOD( ( ij_g - 1 ), ndi) + 1
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
      fsr1 = (fug(ij)-fug(imj))*indr(i)
      IF(fl(imj).NE.1) THEN
        fsr1 = (ug(ij)-ug(imj))*indr(i)
      END IF
      fsr2 = (fvg(ij)-fvg(ijm))*indz(j)
      IF(fl(ijm).NE.1) THEN
        fsr2 = (vg(ij)-vg(ijm))*indz(j)
      END IF
!
      IF(itc.EQ.1) d33 = (ug(ij)+ug(imj))*inr(i)
!
      fu3=0.5D0*(fug(ijp)+fug(imjp))
      fu2=0.5D0*(fug(ij)+fug(imj))
      fu1=0.5D0*(fug(ijm)+fug(imjm))
      fv3=0.5D0*(fvg(ipj)+fvg(ipjm))
      fv2=0.5D0*(fvg(ij)+fvg(ijm))
      fv1=0.5D0*(fvg(imj)+fvg(imjm))
      fum2=fu2+(fu3-fu2)*dz(j)*indzp
      IF(fl(ijp).NE.1) THEN
        u3=0.5D0*(ug(ijp)+ug(imjp))
        u2=0.5D0*(ug(ij)+ug(imj))
        um2=u2+(u3-u2)*dz(j)*indzp
        fum2 = um2
      END IF
      fum1=fu1+(fu2-fu1)*dz(j-1)*indzm
      IF(fl(ijm).NE.1) THEN
        u2=0.5D0*(ug(ij)+ug(imj))
        u1=0.5D0*(ug(ijm)+ug(imjm))
        um1=u1+(u2-u1)*dz(j-1)*indzm
        fum1 = um1
      END IF
      fvm2=fv2+(fv3-fv2)*dr(i)*indrp
      IF(fl(ipj).NE.1) THEN
        v3=0.5D0*(vg(ipj)+vg(ipjm))
        v2=0.5D0*(vg(ij)+vg(ijm))
        vm2=v2+(v3-v2)*dr(i)*indrp
        fvm2 = vm2
      END IF
      fvm1=fv1+(fv2-fv1)*dr(i-1)*indrm
      IF(fl(imj).NE.1) THEN
        v2=0.5D0*(vg(ij)+vg(ijm))
        v1=0.5D0*(vg(imj)+vg(imjm))
        vm1=v1+(v2-v1)*dr(i-1)*indrm
        fvm1 = vm1
      END IF
      fsr12= 0.5D0*((fum2-fum1)*indz(j)+(fvm2-fvm1)*indr(i))
      fmodsr = DSQRT(2*fsr1**2 + 2*fsr2**2 + 4*fsr12**2)
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE u_per_v(uvg)
      USE set_indexes
! 
      IMPLICIT NONE
      INTEGER :: ij
      REAL*8, INTENT(OUT) :: uvg(:)
!
      DO ij = 1, nij_l
        IF(fl_l(ij) .EQ. 1) THEN
          CALL subscl(ij)
          uvg(ij)=(0.5D0*(ug(ij)+ug(imj)))*(0.5D0*(vg(ij)+vg(ijm)))
        ELSE
          uvg(ij) = ug(ij)*vg(ij)
        END IF
      END DO   
      RETURN
      END SUBROUTINE
!
!----------------------------------------------------------------------
       FUNCTION filter(t,ij)
       USE set_indexes
!
       REAL*8,INTENT(IN),DIMENSION(:):: t
       INTEGER,INTENT(IN):: ij
       REAL*8 :: filter
!
       CALL subscl(ij)
       filter = (0.125)*(t(imj)+t(ipj)+t(ijm)+t(ijp))+0.5D0*t(ij)
!
       END FUNCTION filter
!---------------------------------------------------------------------
      SUBROUTINE sgss
!
      USE dimensions
      USE gas_solid_density, ONLY: rlk
      USE gas_solid_velocity, ONLY: uk, vk
      USE grid, ONLY: itc, dz, dr, r, rb, indz, indr, inr, inrb
      USE particles_constants, ONLY: rl, inrl, dk, cmus, nsolid
      USE set_indexes
      IMPLICIT NONE
!
      REAL*8 :: v1, v2, v3, u1, u2, u3, vm1, vm2, um1, um2
      REAL*8 :: d12, d11, d22, d33
      REAL*8 :: drp, dzm, drm, dzp, indrp, indzm, indrm, indzp 
      INTEGER :: ij
      INTEGER :: i, j, k, ij_g
!
      CALL data_exchange(uk)
      CALL data_exchange(vk)
!
      d33=0.D0
      DO ij = 1, nij_l
        ij_g = myij(0, 0, ij)
        IF(fl_l(ij).EQ.1) THEN
         CALL subscl(ij)
         j = ( ij_g - 1 ) / ndi + 1
         i = MOD( ( ij_g - 1 ), ndi) + 1
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
         DO k=1,nsolid
           d11=2.D0*(uk(k,ij)-uk(k,imj))*indr(i)
           d22=2.D0*(vk(k,ij)-vk(k,ijm))*indz(j)
           IF(itc.EQ.1)d33=(uk(k,ij)+uk(k,imj))*inr(i)
           u3=0.5D0*(uk(k,ijp)+uk(k,imjp))
           u2=0.5D0*(uk(k,ij)+uk(k,imj))
           u1=0.5D0*(uk(k,ijm)+uk(k,imjm))
           v3=0.5D0*(vk(k,ipj)+vk(k,ipjm))
           v2=0.5D0*(vk(k,ij)+vk(k,ijm))
           v1=0.5D0*(vk(k,imj)+vk(k,imjm))
           um2=u2+(u3-u2)*dz(j)*indzp
           um1=u1+(u2-u1)*dz(j-1)*indzm
           vm2=v2+(v3-v2)*dr(i)*indrp
           vm1=v1+(v2-v1)*dr(i-1)*indrm
           d12=(um2-um1)*indz(j)+(vm2-vm1)*indr(i)
           mus(k,ij) = cmus(k)*0.7071D0*rl(k)*dk(k)**2*10.D0
     &                     **(3.98D0*rlk(k,ij)*inrl(k)-1.69D0)*
     &                     DSQRT(d11*d11+d22*d22+d33*d33+2.0D0*d12*d12)
         END DO
        END IF
      END DO
!
      RETURN
      END SUBROUTINE
!---------------------------------------------------------------------
       END MODULE turbulence
!---------------------------------------------------------------------
