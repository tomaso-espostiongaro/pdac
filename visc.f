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
      REAL*8, DIMENSION(:),   ALLOCATABLE :: gvisx, gvisy, gvisz  ! gas viscous stresses 
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: pvisx, pvisy, pvisz  ! solid viscous stresses 
!
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE allocate_viscosity
      USE dimensions
      USE grid, ONLY : ncdom
      IMPLICIT NONE
!
      ALLOCATE(mus(ncdom,nsolid))
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
! ... (2D/3D-Compliant)
!
      USE dimensions
      USE gas_constants, ONLY: ckg, phij, mmugs, mmugek, mmug, gmw
      USE gas_constants, ONLY: present_gas
      IMPLICIT NONE
!
      REAL*8, INTENT(IN) :: xgc(:), tg
      REAL*8, INTENT(OUT) :: mu, kap
!
      REAL*8 :: bb, aa, cc, sum, om
      REAL*8 :: tst, tr, t1, t2, t3
      REAL*8 :: m
      INTEGER :: ig, jg
!
! ... Compute Temperature Dependent Viscosity of each gas specie (Reid)
!
      DO ig=1,ngas
        IF (present_gas(ig) ) THEN
!
! ... non-dimensional temperature 
          tst=tg/mmugek(ig)
!
! ... collision integral (omega) 
          om = 1.16145D0 * tst**(-0.14874D0) +  &
               0.52487D0*DEXP(-0.77320D0*tst)+  &
               2.16178D0*DEXP(-2.43787D0*tst)
!
! ... molecular weight (m) expressed in grams per mole
          m = 1.D3 * gmw(ig)
!
! ... semi-empirical consitutive equation for
! ... viscosity is expressed in microPoise (cgs unit system)
          mmug(ig) = 26.69D0 * (m * tg)**0.5D0 / (mmugs(ig)**2 * om)
!
! ... conversion to (Pa-s) (mks unit system)
          mmug(ig) = mmug(ig) * 1.D-7 
!
        END IF
      END DO
!
! ... Temperature Dependent Conductivity (mks unit system)
!
      t1=tg
      t2=tg**2
      t3=tg**3
      tr=tg/132.5D0

      ckg(1) = -3.273D-4 + 9.966D-5*t1 - 3.743D-8*t2 + 9.732D-12*t3
      ckg(2) = 3.919D-4 + 9.816D-5*t1 - 5.067D-8*t2 + 1.504D-11*t3
      ckg(3) = -7.215D-3 + 8.015D-5*t1 + 5.477D-9*t2 - 1.053D-11*t3
      ckg(4) = 8.099D-3 + 6.689D-4*t1 - 4.158D-7*t2 + 1.562D-10*t3
      ckg(5) = 7.341D-3 - 1.013D-5*t1 + 1.801D-7*t2 - 9.100D-11*t3
      ckg(6) = 25.9778D-3*(0.2395D0*tr+6.4977D-3*tr**0.5D0+1.D0       &
               -1.92615D0*tr**(-1.D0)+2.00383D0*tr**(-2.D0)           &
               -1.07553D0*tr**(-3.D0)+0.229414D0*tr**(-4.D0))
      ckg(7) = -8.086D-3 + 6.344D-5*t1 - 1.382D-8*t2 + 2.303D-12*t3
!
! ... Calculation of Mixture Viscosity (Wilke)
! ... and Mixture Conductivity (Mason and Saxema)
!
      DO ig=1,ngas
        IF (present_gas(ig)) THEN
          DO jg=1,ngas
            IF (present_gas(jg)) THEN
              aa=gmw(ig)/gmw(jg)
              bb=mmug(ig)/mmug(jg)
              cc=1.D0+DSQRT(bb)*aa**(-0.25D0)
              phij(ig,jg)=cc**2/DSQRT(8.D0*(1.D0+aa))
            END IF
          END DO
        END IF
      END DO
!
      mu=0.D0
      kap=0.D0
      DO ig=1,ngas
        IF (present_gas(ig)) THEN
          sum=0.D0
          DO jg=1,ngas
            sum=sum+xgc(jg)*phij(ig,jg)
          END DO
          mu = mu + xgc(ig) * mmug(ig) /sum     
          kap = kap + xgc(ig) * ckg(ig) /sum       
        END IF
      END DO
!
      RETURN
      END SUBROUTINE viscon
!----------------------------------------------------------------------
      SUBROUTINE viscg
!----------------------------------------------------------------------
! ... This routine computes the components of the viscous diffusion terms
! ... in gas momentum transport equations
!
      USE control_flags, ONLY: job_type
      USE gas_solid_velocity, ONLY: ug, vg, wg
      USE grid, ONLY: data_exchange
      USE pressure_epsilon, ONLY: ep
      USE turbulence_model, ONLY: mugt, iturb
      IMPLICIT NONE
!
      CALL data_exchange(mug)
      CALL data_exchange(mugt)
!
      IF (iturb >= 1) THEN
        mugt = mug + mugt
      ELSE IF (iturb == 0) THEN
        mugt = mug
      END IF
!
! ... Newtonian stress tensor
!
      IF (job_type == '2D') THEN
        CALL stress2D(gvisx, gvisz, mugt, mug, ep, ug, wg)
      ELSE IF (job_type == '3D') THEN
        CALL stress3D(gvisx, gvisy, gvisz, mugt, mug, ep, ug, vg, wg)
      END IF 
!
      RETURN
      END SUBROUTINE viscg
!----------------------------------------------------------------------
      SUBROUTINE viscs
!----------------------------------------------------------------------
! ... This routine computes the components of the viscous diffusion terms
! ... in particle momentum transport equations
!
      USE control_flags, ONLY: job_type
      USE dimensions
      USE grid, ONLY: fl_l, ncint, ncdom, data_exchange
      USE grid, ONLY: dx, dy, dz
      USE gas_solid_density, ONLY: rlk
      USE gas_solid_velocity, ONLY:  us, vs, ws
      USE particles_constants, ONLY: rl, inrl
      USE set_indexes
      USE turbulence_model, ONLY: must
      USE indijk_module, ONLY: ip0_jp0_kp0_

      IMPLICIT NONE
!
      INTEGER :: is 
!
      REAL*8 :: dxp, dyp, dzp, indxp, indyp, indzp
      REAL*8 :: epsx, epsy, epsz, gepx, gepy, gepz
      INTEGER :: imesh, i, j, k, ij, ijk
      LOGICAL :: repulsive_model
!
      CALL data_exchange(mus)
!
! ... Newtonian stress tensor
!
      DO is = 1, nsolid
        IF (job_type == '2D' ) THEN
          CALL stress2D(pvisx(:,is), pvisz(:,is),                 &
  	                mus(:,is), mus(:,is), rlk(:,is)*inrl(is), &
  		        us(:,is), ws(:,is))
        ELSE IF (job_type == '3D' ) THEN
          CALL stress3D(pvisx(:,is), pvisy(:,is), pvisz(:,is),    &
       	                mus(:,is), mus(:,is), rlk(:,is)*inrl(is), &
  		        us(:,is), vs(:,is), ws(:,is))
        END IF 
      END DO
!
! ... Repulsive model (Gidaspow and Ettehadieh, 1983)
!
      repulsive_model = .TRUE.
      IF ( repulsive_model ) THEN

        IF (job_type == '2D') THEN

          DO ij = 1, ncint
            imesh = myijk(ip0_jp0_kp0_, ij)
            IF(fl_l(ij) == 1) THEN
              CALL subscr(ij)
              j = ( imesh - 1 ) / nx + 1
              i = MOD( ( imesh - 1 ), nx) + 1
!
              dxp=(dx(i)+dx(i+1))
              dzp=(dz(j)+dz(j+1))
!
              indxp=1.D0/dxp
              indzp=1.D0/dzp
!
! ... Coulombic x-gradient
!
              epsx=(dx(i+1)*rlk(ij,is) + dx(i)*rlk(ijke,is)) * indxp * inrl(is)
              gepx=10.D0**(8.76D0*epsx-0.27D0)
              pvisx(ij,is) = pvisx(ij,is) -  & 
                            gepx*indxp*2.D0*(rlk(ijke,is)-rlk(ij,is))*inrl(is)
!
! ... Coulombic z-gradient
!
              epsz=(dz(j+1)*rlk(ij,is) + dz(j)*rlk(ijkt,is)) * indzp * inrl(is)
              gepz=10.D0**(8.76D0*epsz-0.207D0)
              pvisz(ij,is) = pvisz(ij,is) -  & 
                            gepz*indzp*2.D0*(rlk(ijkt,is)-rlk(ij,is))*inrl(is)
            END IF
          END DO

        ELSE IF (job_type == '3D') THEN

          DO ijk = 1, ncint
            imesh = myijk(ip0_jp0_kp0_, ijk)
            IF(fl_l(ijk) == 1) THEN
              CALL subscr(ijk)
              i = MOD( MOD( imesh - 1, nx*ny ), nx ) + 1
   	      j = MOD( imesh - 1, nx*ny ) / nx + 1
   	      k = ( imesh - 1 ) / ( nx*ny ) + 1
! 
              dxp=(dx(i)+dx(i+1))
              dyp=(dy(j)+dy(j+1))
              dzp=(dz(k)+dz(k+1))
! 
              indxp=1.D0/dxp
              indyp=1.D0/dyp
              indzp=1.D0/dzp
 
              DO is = 1, nsolid
!
! ... Coulombic x-gradient
!
                epsx=(dx(i+1)*rlk(ijk,is) + dx(i)*rlk(ijke,is)) * indxp * inrl(is)
                gepx=10.D0**(8.76D0*epsx-0.27D0)
                pvisx(ijk,is) = pvisx(ijk,is) -  & 
                            gepx*indxp*2.D0*(rlk(ijke,is)-rlk(ijk,is))*inrl(is)
!
! ... Coulombic y-gradient
!
                epsy=(dy(j+1)*rlk(ijk,is) + dy(j)*rlk(ijkn,is)) * indyp * inrl(is)
                gepy=10.D0**(8.76D0*epsy-0.207D0)
                pvisy(ijk,is) = pvisy(ijk,is) -  & 
                            gepy*indyp*2.D0*(rlk(ijkn,is)-rlk(ijk,is))*inrl(is)
!
! ... Coulombic z-gradient
!
                epsz=(dz(k+1)*rlk(ijk,is) + dz(k)*rlk(ijkt,is)) * indzp * inrl(is)
                gepz=10.D0**(8.76D0*epsz-0.207D0)
                pvisz(ijk,is) = pvisz(ijk,is) -  & 
                            gepz*indzp*2.D0*(rlk(ijkt,is)-rlk(ijk,is))*inrl(is)
             END DO

            END IF
          END DO

        END IF
      END IF
!
      RETURN
      END SUBROUTINE viscs
!----------------------------------------------------------------------
      SUBROUTINE stress3D(visx, visy, visz, mu, lambda, eps, u, v, w)
!----------------------------------------------------------------------
! ... This routine computes the components of the 3D viscous stress tensor
!
! ... 'mu' and 'lambda' are the first and second viscosity coefficients. 
! ... In general, they can be different.
!
      USE dimensions
      USE grid, ONLY: itc, dx, dy, dz, indx, indy, indz
      USE grid, ONLY: fl_l, ncint
      USE set_indexes
      USE indijk_module, ONLY: ip0_jp0_kp0_

      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: visx(:), visy(:), visz(:)
      REAL*8, INTENT(IN)  :: mu(:), lambda(:), eps(:), u(:), v(:), w(:)
      REAL*8 :: txx2, tyy2, tzz2
      REAL*8 :: txx1, tyy1, tzz1
      REAL*8 :: tyx2, txy2, txz2, tzx2, tyz2, tzy2
      REAL*8 :: tyx1, txy1, txz1, tzx1, tyz1, tzy1
      REAL*8 :: divc, dive, divn, divt
      REAL*8 :: epm1, epm2 
      REAL*8 :: epmxy2, epmxy1, epmyx2, epmyx1
      REAL*8 :: epmxz2, epmxz1, epmzx2, epmzx1
      REAL*8 :: epmyz2, epmyz1, epmzy2, epmzy1
      REAL*8 :: dxm, dxp, dym, dyp, dzm, dzp
      REAL*8 :: indxm, indxp, indym, indyp, indzm, indzp
      INTEGER :: imesh, i, j, k, ijk
!
      visx = 0.D0
      visy = 0.D0
      visz = 0.D0
!
      DO ijk = 1, ncint
        imesh = myijk( ip0_jp0_kp0_, ijk)
        IF(fl_l(ijk) == 1) THEN
         CALL subscr(ijk)
         i = MOD( MOD( imesh - 1, nx*ny ), nx ) + 1
	 j = MOD( imesh - 1, nx*ny ) / nx + 1
	 k = ( imesh - 1 ) / ( nx*ny ) + 1
!
         dxp=(dx(i)+dx(i+1))
         dxm=(dx(i)+dx(i-1))
         dyp=(dy(j)+dy(j+1))
         dym=(dy(j)+dy(j-1))
         dzp=(dz(k)+dz(k+1))
         dzm=(dz(k)+dz(k-1))
!
         indxp=1.D0/dxp
         indxm=1.D0/dxm
         indyp=1.D0/dyp
         indym=1.D0/dym
         indzp=1.D0/dzp
         indzm=1.D0/dzm
!
! ... divergence of the velocity field at centered, east, north and top cells
!
         divc = ( (u(ijk)-u(imjk)) * indx(i) + (v(ijk)-v(ijmk)) * indy(j) +        &
	          (w(ijk)-w(ijkm)) * indz(k) )
         dive = ( (u(ipjk)-u(ijk)) * indx(i+1) + (v(ipjk)-v(ipjmk)) * indy(j) +   &
	          (w(ipjk)-w(ipjkm)) * indz(k) )
         divn = ( (u(ijpk)-u(imjpk)) * indx(i) + (v(ijpk)-v(ijk)) * indy(j+1) +   &
	          (w(ijpk)-w(ijpkm)) * indz(k) )
         divt = ( (u(ijkp)-u(imjkp)) * indx(i) + (v(ijkp)-v(ijmkp)) * indy(j) +   &
	          (w(ijkp)-w(ijk)) * indz(k+1) )
!         
! ... diagonal components of the strain tensor ...
!
         txx2 = 2.D0*(u(ipjk)-u(ijk)) * indx(i+1)  
         txx1 = 2.D0*(u(ijk)-u(imjk)) * indx(i)   
         tyy2 = 2.D0*(v(ijpk)-v(ijk)) * indy(j+1) 
         tyy1 = 2.D0*(v(ijk)-v(ijmk)) * indy(j)   
         tzz2 = 2.D0*(w(ijkp)-w(ijk)) * indz(k+1) 
         tzz1 = 2.D0*(w(ijk)-w(ijkm)) * indz(k)   
!
! ... (isotropy)
!
         txx2 = mu(ijke) * txx2 - 2.D0/3.D0 * lambda(ijke) * dive
         txx1 = mu(ijk)  * txx1 - 2.D0/3.D0 * lambda(ijk)  * divc
         tyy2 = mu(ijkn) * tyy2 - 2.D0/3.D0 * lambda(ijkn) * divn
         tyy1 = mu(ijk)  * tyy1 - 2.D0/3.D0 * lambda(ijk)  * divc
         tzz2 = mu(ijkt) * tzz2 - 2.D0/3.D0 * lambda(ijkt) * divt
         tzz1 = mu(ijk)  * tzz1 - 2.D0/3.D0 * lambda(ijk)  * divc
!
! ... non-diagonal component of the stress tensor
!
         tyx2 = (u(ijpk)-u(ijk)) * indyp*2.D0 + (v(ipjk)-v(ijk)) * indxp*2.D0
         tyx1 = (u(ijk)-u(ijmk)) * indym*2.D0 + (v(ipjmk)-v(ijmk)) * indxp*2.D0
         txy2 = tyx2
         txy1 = (u(imjpk)-u(imjk)) * indyp*2.D0 + (v(ijk)-v(imjk)) * indxm*2.D0
!
         tzy2 = (v(ijkp)-v(ijk)) * indzp*2.D0 + (w(ijpk)-w(ijk)) * indyp*2.D0
         tzy1 = (v(ijk)-v(ijkm)) * indzm*2.D0 + (w(ijpkm)-w(ijkm)) * indyp*2.D0
         tyz2 = tzy2
         tyz1 = (v(ijmkp)-v(ijmk)) * indzp*2.D0 + (w(ijk)-w(ijmk)) * indym*2.D0
! 
         txz2 = (w(ipjk)-w(ijk)) * indxp*2.D0 + (u(ijkp)-u(ijk)) * indzp*2.D0
         txz1 = (w(ijk)-w(imjk)) * indxm*2.D0 + (u(imjkp)-u(imjk)) * indzp*2.D0
         tzx2 = txz2
         tzx1 = (w(ipjkm)-w(ijkm)) * indxp*2.D0 + (u(ijk)-u(ijkm)) * indzm*2.D0
!
! ... compute linearly interpolated values of viscosity on the staggered grid
!
         epm1 =   ( dy(j)*eps(ijkn)*mu(ijkn)   + dy(j+1)*eps(ijk)*mu(ijk) ) * indyp
         epm2 =   ( dy(j)*eps(ijken)*mu(ijken) + dy(j+1)*eps(ijke)*mu(ijke) ) * indyp
         epmyx2 = ( dx(i+1)*epm1 + dx(i)*epm2 ) * indxp
!
         epm1 =   ( dy(j-1)*eps(ijk)*mu(ijk)   + dy(j)*eps(ijks)*mu(ijks) ) * indym
         epm2 =   ( dy(j-1)*eps(ijke)*mu(ijke) + dy(j)*eps(ijkes)*mu(ijkes) ) * indym
         epmyx1 = ( dx(i+1)*epm1 + dx(i)*epm2 ) * indxp
!
         epm1 =   ( dz(k)*eps(ijkt)*mu(ijkt)   + dz(k+1)*eps(ijk)*mu(ijk) ) * indzp
         epm2 =   ( dz(k)*eps(ijket)*mu(ijket) + dz(k+1)*eps(ijke)*mu(ijke) ) * indzp
         epmzx2 = ( dx(i+1)*epm1 + dx(i)*epm2 ) * indxp
!
         epm1 =   ( dz(k-1)*eps(ijk)*mu(ijk)   + dz(k)*eps(ijkb)*mu(ijkb) ) * indzm
         epm2 =   ( dz(k-1)*eps(ijke)*mu(ijke) + dz(k)*eps(ijkeb)*mu(ijkeb) ) * indzm
         epmzx1 = ( dx(i+1)*epm1 + dx(i)*epm2 ) * indxp
! 
! ... X-gradient of the stress tensor
!
         visx(ijk) = (eps(ijke)*txx2 - eps(ijk)*txx1) * indxp*2.D0 +        &
                     (epmyx2*tyx2 - epmyx1*tyx1) * indy(j)         +        &
		     (epmzx2*tzx2 - epmzx1*tzx1) * indz(k)
!
! ... compute linearly interpolated values of viscosity on the staggered grid
!
         epm1 =   ( dy(j)*eps(ijkwn)*mu(ijkwn) + dy(j+1)*eps(ijkw)*mu(ijkw)) * indyp
         epm2 =   ( dy(j)*eps(ijkn)*mu(ijkn)   + dy(j+1)*eps(ijk)*mu(ijk) ) * indyp
         epmxy1 = ( dx(i)*epm1 + dx(i-1)*epm2 ) * indxm
!
         epmxy2 = epmyx2
!
         epm1 =   ( dz(k)*eps(ijkt)*mu(ijkt)   + dz(k+1)*eps(ijk)*mu(ijk) ) * indzp
         epm2 =   ( dz(k)*eps(ijknt)*mu(ijknt) + dz(k+1)*eps(ijkn)*mu(ijkn) ) * indzp
         epmzy2 = ( dy(j+1)*epm1 + dy(j)*epm2 ) * indyp
!
         epm1 =   ( dz(k-1)*eps(ijk)*mu(ijk)   + dz(k)*eps(ijkb)*mu(ijkb) ) * indzm
         epm2 =   ( dz(k-1)*eps(ijkn)*mu(ijkn) + dz(k)*eps(ijknb)*mu(ijknb) ) * indzm
         epmzy1 = ( dy(j+1)*epm1 + dy(j)*epm2 ) * indyp
!
! ... Y-gradient of the stress tensor
!
         visy(ijk) = (eps(ijkn)*tyy2 - eps(ijk)*tyy1) * indyp*2.D0    +      &
                     (epmxy2*txy2 - epmxy1*txy1) * indx(i)            +      &
                     (epmzy2*tzy2 - epmzy1*tzy1) * indz(k) 
!
! ... compute linearly interpolated values of viscosity on the staggered grid
!
         epmxz2 = epmzx2
!
         epm1 =   ( dz(k)*eps(ijkwt)*mu(ijkwt) + dz(k+1)*eps(ijkw)*mu(ijkw)) * indzp
         epm2 =   ( dz(k)*eps(ijkt)*mu(ijkt)   + dz(k+1)*eps(ijk)*mu(ijk) ) * indzp
         epmxz1 = ( dx(i)*epm1 + dx(i-1)*epm2 ) * indxm
!
         epmyz2 = epmzy2
!
         epm1 =   ( dz(k)*eps(ijkst)*mu(ijkst) + dz(k+1)*eps(ijks)*mu(ijks)) * indzp
         epm2 =   ( dz(k)*eps(ijkt)*mu(ijkt)   + dz(k+1)*eps(ijk)*mu(ijk) ) * indzp
         epmyz1 = ( dy(j)*epm1 + dy(j-1)*epm2 ) * indym
!
! ... Z-gradient of the stress tensor
!
         visz(ijk) = (eps(ijkt)*tzz2 - eps(ijk)*tzz1) * indzp*2.D0   +  &
                    (epmxz2*txz2 - epmxz1*txz1) * indx(i)            +  &
                    (epmyz2*tyz2 - epmyz1*tyz1) * indy(j) 

        END IF
      END DO
!
      RETURN
      END SUBROUTINE stress3D
!----------------------------------------------------------------------
      SUBROUTINE stress2D(visx, visz, mu, lambda, eps, u, w)
!----------------------------------------------------------------------
! ... This routine computes the components of the 
! ... 2D cylindrical and cartesian viscous stress tensor
!
! ... 'mu' and 'lambda' are the first and second viscosity coefficients. 
! ... In general, they can be different.
!
      USE dimensions
      USE grid, ONLY: itc, dz, dx, x, xb, indz, indx, inx, inxb
      USE grid, ONLY: fl_l, ncint
      USE set_indexes
      USE indijk_module, ONLY: ip0_jp0_kp0_

      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: visx(:), visz(:)
      REAL*8, INTENT(IN)  :: mu(:), lambda(:), eps(:), u(:), w(:)
      REAL*8 :: t0, mu0, lambda0, eps0
      REAL*8 :: txx2, txx1, tyy2, tyy1
      REAL*8 :: tyx2, tyx1, txy2, txy1
      REAL*8 :: wm1, wm2, dw, du
      REAL*8 :: divc, divr, divt, dive
      REAL*8 :: epsmu2, epsmu11, epsmu22, epsmu1, epsmu12, epsmu21
      REAL*8 :: dxm, dxp, dzm, dzp, indxm, indxp, indzm, indzp
      INTEGER :: imesh, i, j, ij
!
      visx = 0.D0
      visz = 0.D0
!
      DO ij = 1, ncint
        imesh = myijk( ip0_jp0_kp0_, ij)
        IF(fl_l(ij) == 1) THEN
         CALL subscr(ij)
         j = ( imesh - 1 ) / nx + 1
         i = MOD( ( imesh - 1 ), nx) + 1
!
         dxp=(dx(i)+dx(i+1))
         dxm=(dx(i)+dx(i-1))
         dzp=(dz(j)+dz(j+1))
         dzm=(dz(j)+dz(j-1))
!
         indxp=1.D0/dxp
         indxm=1.D0/dxm
         indzp=1.D0/dzp
         indzm=1.D0/dzm
!
! ... divergence of the velocity field at right, centered, top cells
!
         divr = ((xb(i+1)*u(ipjk)-xb(i)*u(ij))*inx(i+1)*indx(i+1)        &
              + (w(ipjk)-w(ipjkm))*indz(j))
         divc = ((xb(i)*u(ij)-xb(i-1)*u(imjk))*inx(i)*indx(i)            &
              + (w(ij)-w(ijkm))*indz(j))
         divt = ((xb(i)*u(ijkp)-xb(i-1)*u(imjkp))*inx(i)*indx(i)          &
              + (w(ijkp)-w(ij))*indz(j+1))
!         
! ... diagonal components of the strain tensor ...
!
         txx2 = 2.D0*(u(ipjk)-u(ij))*indx(i+1)  
         txx1 = 2.D0*(u(ij)-u(imjk))*indx(i)   
         tyy2 = 2.D0*(w(ijkp)-w(ij))*indz(j+1) 
         tyy1 = 2.D0*(w(ij)-w(ijkm))*indz(j)   
!
! ... (isotropy)
!
         txx2 = mu(ijke) * txx2 - 2.D0/3.D0 * lambda(ijke) * divr
         txx1 = mu(ij)  * txx1 - 2.D0/3.D0 * lambda(ij)  * divc
         tyy2 = mu(ijkt) * tyy2 - 2.D0/3.D0 * lambda(ijkt) * divt
         tyy1 = mu(ij)  * tyy1 - 2.D0/3.D0 * lambda(ij)  * divc
!
! ... non-diagonal component of the stress tensor
!
         tyx2 = (u(ijkp)-u(ij))*indzp*2.D0 + (w(ipjk)-w(ij))*indxp*2.D0
         tyx1 = (u(ij)-u(ijkm))*indzm*2.D0 + (w(ipjkm)-w(ijkm))*indxp*2.D0
         txy2 = tyx2
         txy1 = (u(imjkp)-u(imjk))*indzp*2.D0+(w(ij)-w(imjk))*indxm*2.D0
! 
! ... Correction for cylindrical coordinates
!
         IF(itc == 1) THEN
           eps0=(eps(ij)*dx(i+1)+eps(ijke)*dx(i))*indxp
           mu0=(mu(ij)*dx(i+1)+mu(ijke)*dx(i))*indxp
           lambda0=(lambda(ij)*dx(i+1)+lambda(ijke)*dx(i))*indxp
           t0 = 2.D0 * u(ij) * inxb(i)
!
! ... divergence of the velocity field at cell boundary
!
           du = (xb(i+1)*u(ipjk)-xb(i-1)*u(imjk))
           wm2=(dx(i)*w(ipjk)+dx(i+1)*w(ij))*indxp
           wm1=(dx(i)*w(ipjkm)+dx(i+1)*w(ijkm))*indxp
           dw = wm2 - wm1
           dive =  (inxb(i) * du * indxp + dw * indz(j))
!
! ... (isotropy)
!
           t0 = mu0 * t0 - 2.D0/3.D0*lambda0 * dive
         ELSE
           t0 = 0.D0
         ENDIF
!
         epsmu21=(dz(j)*eps(ijkt)*mu(ijkt)+dz(j+1)*eps(ij)*mu(ij))*indzp
         epsmu22=(dz(j)*eps(ijket)*mu(ijket)+dz(j+1)*eps(ijke)*mu(ijke))*indzp
         epsmu2=(dx(i+1)*epsmu21+dx(i)*epsmu22)*indxp
         epsmu11=(dz(j-1)*eps(ij)*mu(ij)+dz(j)*eps(ijkb)*mu(ijkb))*indzm
         epsmu12=(dz(j-1)*eps(ijke)*mu(ijke)+dz(j)*eps(ijkeb)*mu(ijkeb))*indzm
         epsmu1=(dx(i+1)*epsmu11+dx(i)*epsmu12)*indxp
!
! ... x-gradient of the stress tensor
!
         visx(ij) = (eps(ijke)*txx2*x(i+1) - eps(ij)*txx1*x(i))*indxp*2.D0*inxb(i) + &
                   (epsmu2*tyx2-epsmu1*tyx1)*indz(j) - eps0 * t0 * inxb(i)
!
         epsmu11=(dz(j)*eps(ijkwt)*mu(ijkwt)+dz(j+1)*eps(ijkw)*mu(ijkw))*indzp
         epsmu1=(dx(i)*epsmu11+dx(i-1)*epsmu21)*indxm
!
! ... z-gradient of the stress tensor
!
         visz(ij) = (eps(ijkt)*tyy2 - eps(ij)*tyy1) * indzp*2.D0   +  &
                     (xb(i)*epsmu2*txy2-xb(i-1)*epsmu1*txy1)*inx(i)*indx(i) 
        END IF
      END DO
!
      RETURN
      END SUBROUTINE stress2D
!----------------------------------------------------------------------
      END MODULE gas_solid_viscosity
!----------------------------------------------------------------------
