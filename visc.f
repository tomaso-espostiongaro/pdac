!----------------------------------------------------------------------
      MODULE gas_solid_viscosity
!----------------------------------------------------------------------
      IMPLICIT NONE
!
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: particle_viscosity
      REAL*8, DIMENSION(:),  ALLOCATABLE :: gas_viscosity
      REAL*8, DIMENSION(:),  ALLOCATABLE :: gas_thermal_conductivity
!
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: mus   ! particle viscosity
      REAL*8, DIMENSION(:),   ALLOCATABLE :: mug   ! gas molecular viscosity
      REAL*8, DIMENSION(:),   ALLOCATABLE :: kapg  ! gas thermal conductivity
!
      REAL*8, DIMENSION(:),   ALLOCATABLE :: gvisx, gvisz  ! gas viscous stresses 
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: pvisx, pvisz  ! solid viscous stresses 
!
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE bounds_viscosity
      USE dimensions
      IMPLICIT NONE
!
      ALLOCATE(particle_viscosity(nsolid,ntot))
      ALLOCATE(gas_viscosity(ntot))
      ALLOCATE(gas_thermal_conductivity(ntot))
      particle_viscosity = 0.0d0
      gas_viscosity = 0.0d0
      gas_thermal_conductivity = 0.0d0
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE local_bounds_viscosity
      USE dimensions
      USE grid, ONLY : ncdom
      IMPLICIT NONE
!
      ALLOCATE(mus(nsolid, ncdom))
      ALLOCATE(mug(ncdom))
      ALLOCATE(kapg(ncdom))

      mus = 0.d0
      mug = 0.0d0
      kapg = 0.0d0

      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE viscon(mu, kap, xgc, tg)
!----------------------------------------------------------------------
! ... This routine computes molecular viscosity and thermal conductivity
! ... of the gas mixture as a function of temperature
!
      USE dimensions
      USE gas_constants, ONLY: ckg, phij, mmugs, mmugek, mmug, gmw
      IMPLICIT NONE
!
      REAL*8, INTENT(IN) :: xgc(:), tg
      REAL*8, INTENT(OUT) :: mu, kap
!
      REAL*8 :: bb, aa, cc, sum, c2, tst, om, tr, c1 
      INTEGER :: ig, jg
!
! ... Temperature Dependent Viscosity (Reid)
!
!pdac------------
!      DO ig=1,ngas
      DO ig=5,6
!pdac------------
        tst=tg/mmugek(ig)
        om=1.16145D0*tst**(-0.14874D0)+ 0.52487D0*DEXP(-0.77320D0*tst)+  &
     &        2.16178D0*DEXP(-2.43787D0*tst)
        mmug(ig)=26.69D0*(gmw(ig)*tg)**0.5D0/(mmugs(ig)**2*om)
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
!      DO ig=1,ngas
!      DO jg=1,ngas
      DO ig=5,6
      DO jg=5,6
!pdac------------
        aa=gmw(ig)/gmw(jg)
        bb=mmug(ig)/mmug(jg)
        cc=1.D0+DSQRT(bb)*aa**(-0.25D0)
        phij(ig,jg)=cc**2/DSQRT(8.D0*(1.D0+aa))
      END DO
      END DO
!
      c1=0.D0
      c2=0.D0
!pdac------------
!      DO ig=1,ngas
      DO ig=5,6
!pdac------------
        sum=0.D0
!pdac------------
!        DO jg=1,ngas
        DO jg=5,6
!pdac------------
          sum=sum+xgc(jg)*phij(ig,jg)
        END DO
!
! ... Mixture Viscosity (Wilke)
!
        c1 = c1+xgc(ig) * mmug(ig) /sum     
!
! ... Mixture Conductivity (Mason and Saxema)
!
        c2 = c2+xgc(ig) * ckg(ig) /sum       
                                             !
      END DO
!
      mu = c1*1.D-6
      kap = c2*1.D5
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE viscg
!----------------------------------------------------------------------
! ... This routine computes the components of the viscous diffusion terms
! ... in gas momentum transport equations
!
      USE grid, ONLY: ncint
      USE gas_solid_velocity, ONLY: ug, wg
      USE grid, ONLY: data_exchange
      USE pressure_epsilon, ONLY: ep
      USE turbulence, ONLY: mugt, iturb
      IMPLICIT NONE
      INTEGER :: ij
!
      CALL data_exchange(mug)
      CALL data_exchange(mugt)
!
      IF (iturb .GE. 1) THEN
        mugt = mug + mugt
      ELSE
        mugt = mug
      END IF
!
! ... Newtonian stress tensor
!
      CALL stress(gvisx, gvisz, mugt, mug, ep, ug, wg)
!
      RETURN
      END SUBROUTINE

! ... MODIFICARE_X3D ( fino fine file )

!----------------------------------------------------------------------
      SUBROUTINE viscs(is)
!----------------------------------------------------------------------
! ... This routine computes the components of the viscous diffusion terms
! ... in particle momentum transport equations
!
      USE dimensions
      USE grid, ONLY: fl_l, ncint, ncdom, data_exchange
      USE grid, ONLY: dz, dr
      USE gas_solid_density, ONLY: rlk
      USE gas_solid_velocity, ONLY:  us, ws
      USE particles_constants, ONLY: rl, inrl
      USE set_indexes
      USE turbulence, ONLY: must
      USE indijk_module, ONLY: ip0_jp0_kp0_

      IMPLICIT NONE
!
      REAL*8, ALLOCATABLE :: eps(:)
      INTEGER, INTENT(IN) :: is 
!
      REAL*8 :: drp, dzp, indrp, indzp
      REAL*8 :: epsx, epsz, gepx, gepz
      INTEGER :: imesh, i, j, ij
      LOGICAL :: repulsive_model
!
      ALLOCATE(eps(ncdom))
      eps(:) = rlk(is,:)*inrl(is)
!
      CALL data_exchange(mus)
!
! ... Newtonian stress tensor
!
      CALL stress(pvisx(is,:), pvisz(is,:), mus(is,:), mus(is,:), eps(:), &
                         us(is,:), ws(is,:))
!
! ... Repulsive model (Gidaspow and Ettehadieh, 1983)
!
      repulsive_model = .TRUE.
      IF ( repulsive_model ) THEN
        DO ij = 1, ncint
         imesh = myijk(ip0_jp0_kp0_, ij)
         IF(fl_l(ij).EQ.1) THEN
           CALL subscr(ij)
           j = ( imesh - 1 ) / nr + 1
           i = MOD( ( imesh - 1 ), nr) + 1
!
           drp=(dr(i)+dr(i+1))
           dzp=(dz(j)+dz(j+1))
!
           indrp=1.D0/drp
           indzp=1.D0/dzp
!
! ... Coulombic x-gradient
!
           epsx=(dr(i+1)*rlk(is,ij) + dr(i)*rlk(is,ijr)) * indrp * inrl(is)
           gepx=10.D0**(8.76D0*epsx-0.27D0)
           pvisx(is,ij) = pvisx(is,ij) -  & 
                         gepx*indrp*2.D0*(rlk(is,ijr)-rlk(is,ij))*inrl(is)
!
! ... Coulombic z-gradient
!
           epsz=(dz(j+1)*rlk(is,ij) + dz(j)*rlk(is,ijt)) * indzp * inrl(is)
           gepz=10.D0**(8.76D0*epsz-0.207D0)
           pvisz(is,ij) = pvisz(is,ij) -  & 
                         gepz*indzp*2.D0*(rlk(is,ijt)-rlk(is,ij))*inrl(is)
         END IF
        END DO
      END IF
!
      DEALLOCATE(eps)
!
      RETURN
      END SUBROUTINE

!----------------------------------------------------------------------
      SUBROUTINE stress(visx, visz, mu, lambda, eps, u, v)
!----------------------------------------------------------------------
! ... This routine computes the components of the viscous stress tensor
!
      USE dimensions
      USE grid, ONLY: itc, dz, dr, r, rb, indz, indr, inr, inrb
      USE grid, ONLY: fl_l, ncint
      USE set_indexes
      USE indijk_module, ONLY: ip0_jp0_kp0_

      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: visx(:), visz(:)
      REAL*8, INTENT(IN)  :: mu(:), lambda(:), eps(:), u(:), v(:)
      REAL*8 :: t0, mu0, lambda0, eps0
      REAL*8 :: txx2, txx1, tyy2, tyy1
      REAL*8 :: tyx2, tyx1, txy2, txy1
      REAL*8 :: vm1, vm2, dv, du
      REAL*8 :: divc, divr, divt, dive
      REAL*8 :: epsmu2, epsmu11, epsmu22, epsmu1, epsmu12, epsmu21
      REAL*8 :: drm, drp, dzm, dzp, indrm, indrp, indzm, indzp
      INTEGER :: imesh, i, j, ij
!
      visx = 0.D0
      visz = 0.D0
!
      DO ij = 1, ncint
        imesh = myijk( ip0_jp0_kp0_, ij)
        IF(fl_l(ij).EQ.1) THEN
         CALL subscr(ij)
         j = ( imesh - 1 ) / nr + 1
         i = MOD( ( imesh - 1 ), nr) + 1
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
! ... divergence of the velocity field at right, centered, top cells
!
         divr = ((rb(i+1)*u(ipj)-rb(i)*u(ij))*inr(i+1)*indr(i+1)        &
              + (v(ipj)-v(ipjm))*indz(j))
         divc = ((rb(i)*u(ij)-rb(i-1)*u(imj))*inr(i)*indr(i)            &
              + (v(ij)-v(ijm))*indz(j))
         divt = ((rb(i)*u(ijp)-rb(i-1)*u(imjp))*inr(i)*indr(i)          &
              + (v(ijp)-v(ij))*indz(j+1))
!         
! ... diagonal components of the strain tensor ...
!
         txx2 = 2.D0*(u(ipj)-u(ij))*indr(i+1)  
         txx1 = 2.D0*(u(ij)-u(imj))*indr(i)   
         tyy2 = 2.D0*(v(ijp)-v(ij))*indz(j+1) 
         tyy1 = 2.D0*(v(ij)-v(ijm))*indz(j)   
!
! ... traceless 
!
         txx2 = mu(ijr) * txx2 - 2.D0/3.D0 * lambda(ijr) * divr
         txx1 = mu(ij)  * txx1 - 2.D0/3.D0 * lambda(ij)  * divc
         tyy2 = mu(ijt) * tyy2 - 2.D0/3.D0 * lambda(ijt) * divt
         tyy1 = mu(ij)  * tyy1 - 2.D0/3.D0 * lambda(ij)  * divc
!
! ... non-diagonal component of the stress tensor
!
         tyx2 = (u(ijp)-u(ij))*indzp*2.D0 + (v(ipj)-v(ij))*indrp*2.D0
         tyx1 = (u(ij)-u(ijm))*indzm*2.D0 + (v(ipjm)-v(ijm))*indrp*2.D0
         txy2 = tyx2
         txy1 = (u(imjp)-u(imj))*indzp*2.D0+(v(ij)-v(imj))*indrm*2.D0
! 
! ... Correction for cylindrical coordinates
!
         IF(itc.EQ.1) THEN
           eps0=(eps(ij)*dr(i+1)+eps(ijr)*dr(i))*indrp
           mu0=(mu(ij)*dr(i+1)+mu(ijr)*dr(i))*indrp
           lambda0=(lambda(ij)*dr(i+1)+lambda(ijr)*dr(i))*indrp
           t0 = 2.D0 * u(ij) * inrb(i)
!
! ... divergence of the velocity field at cell boundary
!
           du = (rb(i+1)*u(ipj)-rb(i-1)*u(imj))
           vm2=(dr(i)*v(ipj)+dr(i+1)*v(ij))*indrp
           vm1=(dr(i)*v(ipjm)+dr(i+1)*v(ijm))*indrp
           dv = vm2 - vm1
           dive =  (inrb(i) * du * indrp + dv * indz(j))
!
! ... traceless
!
           t0 = mu0 * t0 - 2.D0/3.D0*lambda0 * dive
         ELSE
           t0 = 0.D0
         ENDIF
!
         epsmu21=(dz(j)*eps(ijt)*mu(ijt)+dz(j+1)*eps(ij)*mu(ij))*indzp
         epsmu22=(dz(j)*eps(ijtr)*mu(ijtr)+dz(j+1)*eps(ijr)*mu(ijr))*indzp
         epsmu2=(dr(i+1)*epsmu21+dr(i)*epsmu22)*indrp
         epsmu11=(dz(j-1)*eps(ij)*mu(ij)+dz(j)*eps(ijb)*mu(ijb))*indzm
         epsmu12=(dz(j-1)*eps(ijr)*mu(ijr)+dz(j)*eps(ijbr)*mu(ijbr))*indzm
         epsmu1=(dr(i+1)*epsmu11+dr(i)*epsmu12)*indrp
!
! ... x-gradient of the stress tensor
!
         visx(ij) = (eps(ijr)*txx2*r(i+1) - eps(ij)*txx1*r(i))*indrp*2.D0*inrb(i) + &
                   (epsmu2*tyx2-epsmu1*tyx1)*indz(j) - eps0 * t0 * inrb(i)
!
         epsmu11=(dz(j)*eps(ijtl)*mu(ijtl)+dz(j+1)*eps(ijl)*mu(ijl))*indzp
         epsmu1=(dr(i)*epsmu11+dr(i-1)*epsmu21)*indrm
!
! ... z-gradient of the stress tensor
!
         visz(ij) = (eps(ijt)*tyy2 - eps(ij)*tyy1) * indzp*2.D0   +  &
                     (rb(i)*epsmu2*txy2-rb(i-1)*epsmu1*txy1)*inr(i)*indr(i) 
        END IF
      END DO
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      END MODULE gas_solid_viscosity
!----------------------------------------------------------------------
