!----------------------------------------------------------------------
      MODULE gas_solid_viscosity
!----------------------------------------------------------------------
      IMPLICIT NONE
!
      REAL*8, DIMENSION(:),  ALLOCATABLE :: gas_viscosity
      REAL*8, DIMENSION(:),  ALLOCATABLE :: gas_thermal_conductivity
!
      REAL*8, DIMENSION(:),   ALLOCATABLE :: mug   !molecular viscosity (gas)
      REAL*8, DIMENSION(:),   ALLOCATABLE :: kapg  !thermal conductivity (gas)
!
      REAL*8, DIMENSION(:),   ALLOCATABLE :: gvisx, gvisz  !viscous stress (gas)
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: pvisx, pvisz  !viscous stress (part.)
!
      INTEGER :: icoh
!
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE bounds_viscosity
      USE dimensions
      IMPLICIT NONE
!
      ALLOCATE(gas_viscosity(nr*nz))
      ALLOCATE(gas_thermal_conductivity(nr*nz))
      gas_viscosity = 0.0d0
      gas_thermal_conductivity = 0.0d0
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE local_bounds_viscosity
      USE dimensions
      USE grid, ONLY : nijx_l
      IMPLICIT NONE
!
      ALLOCATE(mug(nijx_l))
      ALLOCATE(kapg(nijx_l))

      mug = 0.0d0
      kapg = 0.0d0

      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE viscon(mug, kapg, xgc, tg)
!----------------------------------------------------------------------
! ... This routine computes molecular viscosity and thermal conductivity
! ... of the gas mixture as a function of temperature
!
      USE dimensions
      USE gas_constants, ONLY: ckg, phij, mmugs, mmugek, mmug, gmw
      IMPLICIT NONE
!
      REAL*8, INTENT(IN) :: xgc(:), tg
      REAL*8, INTENT(OUT) :: mug, kapg
!
      REAL*8 :: bb, aa, cc, sum, c2, tst, om, tr, c1 
      INTEGER :: ik, jk
!
! ... Temperature Dependent Viscosity (Reid)
!
!pdac------------
!      DO ik=1,ngas
      DO ik=5,6
!pdac------------
        tst=tg/mmugek(ik)
        om=1.16145D0*tst**(-0.14874D0)+ 0.52487D0*DEXP(-0.77320D0*tst)+  &
     &        2.16178D0*DEXP(-2.43787D0*tst)
        mmug(ik)=26.69D0*(gmw(ik)*tg)**0.5D0/(mmugs(ik)**2*om)
      END DO
!
! ... Temperature Dependent Conductivity (Wassilijewa)
!
      c1=tg
      tr=tg/132.5D0
!pdac------------
!      ckg(1)=-3.273D-4+9.966D-5*c1-3.743D-8*c1**2+9.732D-12*c1**3
!      ckg(2)=3.919D-4+9.816D-5*c1-5.067D-8*c1**2+1.504D-11*c1**3
!      ckg(3)=-7.215D-3+8.015D-5*c1+5.477D-9*c1**2-1.053D-11*c1**3
!      ckg(4)=8.099D-3+6.689D-4*c1-4.158D-7*c1**2+1.562D-10*c1**3
!pdac------------
      ckg(5)=7.341D-3-1.013D-5*c1+1.801D-7*c1**2-9.100D-11*c1**3
      ckg(6)=25.9778D-3*(0.2395D0*tr+6.4977D-3*tr**0.5D0+1.D0      &
     &      -1.92615D0*tr**(-1.D0)+2.00383D0*tr**(-2.D0)           &
     &      -1.07553D0*tr**(-3.D0)+0.229414D0*tr**(-4.D0))
!pdac------------
!      ckg(7)=-8.086D-3+6.344D-5*c1-1.382D-8*c1**2+2.303D-12*c1**3
!pdac------------
!
!pdac------------
!      DO ik=1,ngas
!      DO jk=1,ngas
      DO ik=5,6
      DO jk=5,6
!pdac------------
        aa=gmw(ik)/gmw(jk)
        bb=mmug(ik)/mmug(jk)
        cc=1.D0+DSQRT(bb)*aa**(-0.25D0)
        phij(ik,jk)=cc**2/DSQRT(8.D0*(1.D0+aa))
      END DO
      END DO
!
      c1=0.D0
      c2=0.D0
!pdac------------
!      DO ik=1,ngas
      DO ik=5,6
!pdac------------
        sum=0.D0
!pdac------------
!        DO jk=1,ngas
        DO jk=5,6
!pdac------------
          sum=sum+xgc(jk)*phij(ik,jk)
        END DO
!
! ... Mixture Viscosity (Wilke)
!
        c1 = c1+xgc(ik) * mmug(ik) /sum     
!
! ... Mixture Conductivity (Mason and Saxema)
!
        c2 = c2+xgc(ik) * ckg(ik) /sum       
                                             !
      END DO
!
      mug = c1*1.D-6
      kapg = c2*1.D5
!
      RETURN
      END SUBROUTINE

!----------------------------------------------------------------------
      SUBROUTINE viscg
!----------------------------------------------------------------------
! ... This routine computes the components of the viscous stress tensor
! ... of the gas phase, and its spatial derivatives
!
      USE grid, ONLY: fl_l
      USE dimensions
      USE gas_solid_velocity, ONLY: ug, vg
      USE grid, ONLY: itc, dz, dr, r, rb, indz, indr, inr, inrb
      USE grid, ONLY: nij_l, myij
      USE pressure_epsilon, ONLY: ep
      USE set_indexes
      USE turbulence, ONLY: mugt
      IMPLICIT NONE
!
      INTEGER :: ij_g, i, j
      INTEGER :: ij
      REAL*8 :: ugm1, vgm2, vgm1, ugm2, d1c, gmum, t0, vmu12 
      REAL*8 :: vmu1, vdu, vmu11, vmu21, vmu22, vmu2
      REAL*8 :: txx2, tyy1, txy2, txy1, tyy2, txx1 
      REAL*8 :: tyx2, tyx1
      REAL*8 :: drm, drp, dzm, dzp, indrm, indrp, indzm, indzp
!
      mugt = mugt + mug
      
      DO ij = 1, nij_l
        ij_g = myij(0, 0, ij)
        IF(fl_l(ij).EQ.1) THEN
         CALL subscl(ij) 
         j = ( ij_g - 1 ) / nr + 1
         i = MOD( ( ij_g - 1 ), nr) + 1
!
         drp=(dr(i)+dr(i+1))
         drm=(dr(i)+dr(i-1))
         dzp=(dz(j)+dz(j+1))
         dzm=(dz(j)+dz(j-1))
!
         indrp=1.D0/drp
         indrm=1.D0/drm
         indzp=1.D0/dzp
         indzm=1.D0/dzm
!         
         txx2=2.D0*mugt(ijr)*(ug(ipj)-ug(ij))*indr(i+1)
         txx1=2.D0*mugt(ij)*(ug(ij)-ug(imj))*indr(i)
         tyx2=(ug(ijp)-ug(ij))*indzp*2.D0+(vg(ipj)-vg(ij))*indrp*2.D0
         tyx1=(ug(ij)-ug(ijm))*indzm*2.D0+(vg(ipjm)-vg(ijm))*indrp*2.D0
         tyy2=2.D0*mugt(ijt)*(vg(ijp)-vg(ij))*indz(j+1)
         tyy1=2.D0*mugt(ij)*(vg(ij)-vg(ijm))*indz(j)
         txy2=tyx2
         txy1=(ug(imjp)-ug(imj))*indzp*2.D0+(vg(ij)-vg(imj))*indrm*2.D0
         txx2=txx2-2.D0/3.D0*mug(ijr)*                          &
     &       ((rb(i+1)*ug(ipj)-rb(i)*ug(ij))*inr(i+1)*indr(i+1) &
     &       +(vg(ipj)-vg(ipjm))*indz(j))
         d1c=2.D0/3.D0*mug(ij)*((rb(i)*ug(ij)-rb(i-1)*          &
     &        ug(imj))*inr(i)*indr(i) + (vg(ij)-vg(ijm))*indz(j))
         txx1=txx1-d1c
! 
         IF(itc.EQ.1) THEN
           gmum=mugt(ij)+(mugt(ijr)-mugt(ij))*dr(i)*indrp
           t0=2.D0*gmum*ug(ij)*inrb(i)
           ugm2=0.5D0*(rb(i+1)*ug(ipj)+rb(i)*ug(ij))
           ugm1=0.5D0*(rb(i)*ug(ij)+rb(i-1)*ug(imj))
           vgm2=(dr(i)*vg(ipj)+dr(i+1)*vg(ij))*indrp
           vgm1=(dr(i)*vg(ipjm)+dr(i+1)*vg(ijm))*indrp
           gmum=mug(ij)+(mug(ijr)-mug(ij))*dr(i)*indrp
           t0=t0-2.D0/3.D0*gmum*((ugm2-ugm1)*indrp*2.D0*inrb(i)+(vgm2-vgm1)*indz(j))
         ENDIF
!
         tyy2=tyy2-2.D0/3.D0*mug(ijt)*                          &
     &       ((rb(i)*ug(ijp)-rb(i-1)*ug(imjp))*inr(i)*indr(i)+  &
     &        (vg(ijp)-vg(ij))*indz(j+1))
         tyy1=tyy1-d1c
         vmu21=(dz(j)*ep(ijt)*mugt(ijt)+dz(j+1)*ep(ij)*mugt(ij))*indzp
         vmu22=(dz(j)*ep(ijtr)*mugt(ijtr)+dz(j+1)*ep(ijr)*mugt(ijr))*indzp
         vmu2=(dr(i+1)*vmu21+dr(i)*vmu22)*indrp
         vmu11=(dz(j-1)*ep(ij)*mugt(ij)+dz(j)*ep(ijb)*mugt(ijb))*indzm
         vmu12=(dz(j-1)*ep(ijr)*mugt(ijr)+dz(j)*ep(ijbr)*mugt(ijbr))*indzm
         vmu1=(dr(i+1)*vmu11+dr(i)*vmu12)*indrp
         vdu=(dr(i)*ep(ijr)+dr(i+1)*ep(ij))*indrp
!
! gvisx and gvisz already have the gas porosity in
!
         gvisx(ij)=(txx2*ep(ijr)*r(i+1)-txx1*ep(ij)*r(i))       &
     &             *indrp*2.D0*inrb(i)+(vmu2*tyx2-vmu1*tyx1)*indz(j)-t0*vdu*inrb(i)
         vmu11=(dz(j)*ep(ijtl)*mugt(ijtl)+dz(j+1)*ep(ijl)*mugt(ijl))*indzp
         vmu1=(dr(i)*vmu11+dr(i-1)*vmu21)*indrm
         gvisz(ij)=(rb(i)*vmu2*txy2-rb(i-1)*vmu1*txy1)*inr(i)*indr(i)  &
     &             +(ep(ijt)*tyy2-ep(ij)*tyy1)*indzp*2.D0
        END IF
      END DO
!
      RETURN
      END SUBROUTINE
!
!----------------------------------------------------------------------
      SUBROUTINE viscs(k)
!----------------------------------------------------------------------
! ... This routine computes the components of the viscous stress tensor
! ... of the solid phases, and their spatial derivatives
!
      USE grid, ONLY: fl_l
      USE dimensions
      USE gas_solid_density, ONLY: rlk
      USE gas_solid_velocity, ONLY:  uk, vk
      USE grid, ONLY: itc, dz, dr, r, rb, indz, indr, inr, inrb
      USE grid, ONLY: nij_l, myij
      USE particles_constants, ONLY: rl, inrl
      USE set_indexes
      USE turbulence, ONLY: mus
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: k 
!
      INTEGER :: j, i, ij_g
      INTEGER :: ij
      REAL*8 :: ukm1, vkm2, vkm1, ukm2, d1c, gmum, t0 
      REAL*8 :: tauc1, vmu11, vmu12, vmu1, vmu2, tauc2, vmu21, vmu22 
      REAL*8 :: txx2, tyy1, txy2, txy1, tyy2 
      REAL*8 :: txx1, tyx2, tyx1 
      REAL*8 :: drp, dzm, drm, dzp, indrp, indzm, indrm, indzp 
!
      DO ij = 1, nij_l
        ij_g = myij(0, 0, ij)
        IF(fl_l(ij).EQ.1) THEN
         CALL subscl(ij)
         j = ( ij_g - 1 ) / nr + 1
         i = MOD( ( ij_g - 1 ), nr) + 1
!
         drp=(dr(i)+dr(i+1))
         drm=(dr(i)+dr(i-1))
         dzp=(dz(j)+dz(j+1))
         dzm=(dz(j)+dz(j-1))
!
         indrp=1.D0/drp
         indrm=1.D0/drm
         indzp=1.D0/dzp
         indzm=1.D0/dzm
!         
         txx2=2.D0*mus(k,ijr)*indr(i+1)*(uk(k,ipj)-uk(k,ij))
         txx1=2.D0*mus(k,ij)*indr(i)*(uk(k,ij)-uk(k,imj))
         tyx2=(uk(k,ijp)-uk(k,ij))*indzp*2.D0+(vk(k,ipj)-vk(k,ij))*indrp*2.D0 
         tyx1=(uk(k,ij)-uk(k,ijm))*indzm*2.D0+(vk(k,ipjm)-vk(k,ijm))*indrp*2.D0
         tyy2=2.D0*mus(k,ijt)*(vk(k,ijp)-vk(k,ij))*indz(j+1)
         tyy1=2.D0*mus(k,ij)*(vk(k,ij)-vk(k,ijm))*indz(j)
         txy2=tyx2
         txy1=(uk(k,imjp)-uk(k,imj))*indzp*2.D0+(vk(k,ij)-vk(k,imj))*indrm*2.D0
         txx2 = txx2 -  2.D0/3.D0*mus(k,ijr) *                         &
     &         ((rb(i+1)*uk(k,ipj)-rb(i)*uk(k,ij))*inr(i+1)*indr(i+1)  &
     &        + (vk(k,ipj)-vk(k,ipjm))*indz(j))
         d1c = 2.D0/3.D0*mus(k,ij) *                                   &
     &         ((rb(i)*uk(k,ij)-rb(i-1)*uk(k,imj))*inr(i)*indr(i)      &
     &     +(vk(k,ij)-vk(k,ijm))*indz(j))
         txx1=txx1-d1c
!
         IF(itc.EQ.1) THEN
           gmum=(dr(i)*mus(k,ijr)*rlk(k,ijr)+dr(i+1)*mus(k,ij)*rlk(k,ij))*indrp*inrl(k)
           t0=gmum*uk(k,ij)*inrb(i)
           ukm2=0.5D0*(rb(i+1)*uk(k,ipj)+rb(i)*uk(k,ij))
           ukm1=0.5D0*(rb(i)*uk(k,ij)+rb(i-1)*uk(k,imj))
           vkm2=(dr(i)*vk(k,ipj)+dr(i+1)*vk(k,ij))*indrp
           vkm1=(dr(i)*vk(k,ipjm)+dr(i+1)*vk(k,ijm))*indrp
           t0=t0-2.D0/3.D0*gmum*((ukm2-ukm1)*2.D0*indrp*inrb(i)+(vkm2-vkm1)*indz(j))
         ENDIF
!
         tyy2=tyy2-2.D0/3.D0*mus(k,ijt)                                &
     &      *((rb(i)*uk(k,ijp)-rb(i-1)*uk(k,imjp))*inr(i)*indr(i)      &
     &      +(vk(k,ijp)-vk(k,ij))*indz(j+1))
         tyy1=tyy1-d1c
!
         IF(icoh.EQ.1.AND.rl(k).LE.4.D-3) THEN
           tauc1=10.D0**(10.D0*rlk(k,ij)*inrl(k)-4.5D0)
           tauc2=10.D0**(10.D0*rlk(k,ijr)*inrl(k)-4.5D0)
           txx2=txx2+tauc2
           txx1=txx1+tauc1
           tauc2=10.D0**(10.D0*rlk(k,ijt)*inrl(k)-4.5D0)
           tyy2=tyy2+tauc2
           tyy1=tyy1+tauc1
         ENDIF
!
         vmu21=(dz(j)*mus(k,ijt)*rlk(k,ijt)+dz(j+1)*mus(k,ij)*rlk(k,ij))*indzp
         vmu22=(dz(j)*mus(k,ijtr)*rlk(k,ijtr)+dz(j+1)*mus(k,ijr)*rlk(k,ijr))*indzp
         vmu2=(dr(i+1)*vmu21+dr(i)*vmu22)*indrp*inrl(k)
         vmu11=(dz(j-1)*mus(k,ij)*rlk(k,ij)+dz(j)*mus(k,ijb)*rlk(k,ijb))*indzm
         vmu12=(dz(j-1)*mus(k,ijr)*rlk(k,ijr)+dz(j)*mus(k,ijbr)*rlk(k,ijbr))*indzm
         vmu1=(dr(i+1)*vmu11+dr(i)*vmu12)*indrp*inrl(k)
!
         pvisx(k,ij) = (rlk(k,ijr)*txx2*r(i+1)                            &
     &                 -rlk(k,ij)*txx1*r(i))*indrp*2.D0*inrb(i)*inrl(k)+  &
     &                 (vmu2*tyx2-vmu1*tyx1)*indz(j)-t0*inrb(i)
         vmu11=(dz(j)*mus(k,ijtl)*rlk(k,ijtl)+                            &
     &        dz(j+1)*mus(k,ijl)*rlk(k,ijl))*indzp
         vmu1=(dr(i)*vmu11+dr(i-1)*vmu21)*indrm*inrl(k)
!
         pvisz(k,ij) = (rb(i)*vmu2*txy2-rb(i-1)*vmu1*txy1)*inr(i)*indr(i) &
     &                 +(rlk(k,ijt)*tyy2-rlk(k,ij)*tyy1)*indzp*2.D0*inrl(k)
        END IF
      END DO
!
      RETURN
      END SUBROUTINE
!
!----------------------------------------------------------------------
      END MODULE gas_solid_viscosity
!----------------------------------------------------------------------
