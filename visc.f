!----------------------------------------------------------------------
      MODULE gas_solid_viscosity
!----------------------------------------------------------------------
      IMPLICIT NONE
!
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: mus   ! particle viscosity
      REAL*8, DIMENSION(:),   ALLOCATABLE :: mug   ! gas molecular viscosity
      REAL*8, DIMENSION(:),   ALLOCATABLE :: kapg  ! gas thermal conductivity
!
      REAL*8, DIMENSION(:),   ALLOCATABLE :: gvisx, gvisy, gvisz  ! gas viscous stresses 
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: pvisx, pvisy, pvisz  ! solid viscous stresses 
!
      LOGICAL :: gas_viscosity
      LOGICAL :: part_viscosity
      INTEGER :: repulsive_model

      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE allocate_viscosity
      USE dimensions
      USE domain_decomposition, ONLY : ncdom
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
! ... This routine computes the molecular (dynamic) viscosity and the 
! ... thermal conductivity of the gas mixture as a function of temperature.
! ... The computed values of viscosity and conductivity are scaled by a 
! ... factor (1D-1) and (1D-5) respectively changing from CGS to MKS
! ... (International) Unit System.
! ... (2D/3D-Compliant)
!
      USE dimensions
      USE gas_constants, ONLY: ckg, phij, mmugs, mmugek, mmug, gmw
      USE gas_constants, ONLY: gas_type
      IMPLICIT NONE
!
      REAL*8, INTENT(IN) :: xgc(:), tg
      REAL*8, INTENT(OUT) :: mu, kap
!
      REAL*8 :: bb, aa, cc, sum, om
      REAL*8 :: tst, tr, t1, t2, t3
      REAL*8 :: trm1, trm2, trm3, trm4, trsrt
      REAL*8 :: m
      REAL*8 :: suminv

      INTEGER :: ig, jg, igg, jgg
!
! ... Compute Temperature Dependent Viscosity of each gas specie (Reid)
!
      DO ig = 1, ngas
        igg = gas_type(ig)
!
! ... non-dimensional temperature 
           
          tst = tg / mmugek(igg)
!
! ... collision integral (omega) 
          om = 1.16145D0 * tst ** ( -0.14874D0 ) +  &
               0.52487D0 * DEXP( -0.77320D0 * tst ) +  &
               2.16178D0 * DEXP( -2.43787D0 * tst )
!
! ... molecular weight (m) expressed in grams per mole
          m = 1.D3 * gmw(igg)
!
! ... semi-empirical consitutive equation for
! ... viscosity is expressed in microPoise (cgs unit system)
          mmug(igg) = 26.69D0 * SQRT( m * tg ) /    &
                              ( mmugs(igg)**2 * om )
!
! ... conversion from (microPoise) to (Pa-s) (mks unit system)
          mmug(igg) = mmug(igg) * 1.D-7 
!
      END DO
!
! ... Temperature Dependent Conductivity (mks unit system)
!
      t1 = tg             ! tg
      t2 = t1 * tg        ! tg**2
      t3 = t2 * tg        ! tg**3

      tr = tg / 132.5D0

      trm1 = 1.0d0 / tr   ! tr**(-1)
      trm2 = trm1 * trm1  ! tr**(-2)
      trm3 = trm2 * trm1  ! tr**(-3)
      trm4 = trm3 * trm1  ! tr**(-4)

      trsrt = SQRT( tr )

      ckg(1) = -3.273D-4 + 9.966D-5 * t1 - 3.743D-8 * t2 + 9.732D-12 * t3
      ckg(2) =  3.919D-4 + 9.816D-5 * t1 - 5.067D-8 * t2 + 1.504D-11 * t3
      ckg(3) = -7.215D-3 + 8.015D-5 * t1 + 5.477D-9 * t2 - 1.053D-11 * t3
      ckg(4) =  8.099D-3 + 6.689D-4 * t1 - 4.158D-7 * t2 + 1.562D-10 * t3
      ckg(5) =  7.341D-3 - 1.013D-5 * t1 + 1.801D-7 * t2 - 9.100D-11 * t3
      ckg(6) = 25.9778D-3 * ( &
                0.2395D0 * tr + 6.4977D-3 * trsrt + 1.D0       &
              - 1.92615D0 * trm1 + 2.00383D0  * trm2           &
              - 1.07553D0 * trm3 + 0.229414D0 * trm4 )
      ckg(7) = -8.086D-3 + 6.344D-5 * t1 - 1.382D-8 * t2 + 2.303D-12 * t3
!
! ... Calculation of Mixture Viscosity (Wilke)
! ... and Mixture Conductivity (Mason and Saxema)
!
      DO ig=1,ngas
        igg = gas_type(ig)
        DO jg=1,ngas
          jgg = gas_type(jg)

          aa = gmw(igg) / gmw(jgg)
          bb = mmug(igg) / mmug(jgg)
          cc = 1.D0 + DSQRT(bb) * aa**(-0.25D0)
          phij(igg,jgg) = cc**2 / DSQRT(8.D0 * (1.D0 + aa))

        END DO
      END DO
!
      mu=0.D0
      kap=0.D0
      DO ig=1,ngas
        igg = gas_type(ig)
        sum=0.D0
        DO jg=1,ngas
          jgg = gas_type(jg)
          sum=sum+xgc(jg)*phij(igg,jgg)
        END DO
        suminv = 1.0D0 / sum

        mu = mu + xgc(ig) * mmug(gas_type(ig)) * suminv   
        kap = kap + xgc(ig) * ckg(gas_type(ig)) * suminv       

      END DO
!
      RETURN
      END SUBROUTINE viscon
!----------------------------------------------------------------------
      SUBROUTINE viscg
! ... This routine computes the components of the viscous diffusion terms
! ... in gas momentum transport equations
!
      USE control_flags, ONLY: job_type
      USE gas_solid_velocity, ONLY: ug, vg, wg
      USE domain_decomposition, ONLY: data_exchange
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
! ... This routine computes the components of the viscous diffusion terms
! ... in particle momentum transport equations
!
      USE control_flags, ONLY: job_type
      USE dimensions
      USE domain_decomposition, ONLY: ncint, ncdom, data_exchange
      USE domain_decomposition, ONLY: myijk, meshinds
      USE grid, ONLY: dx, dy, dz, flag
      USE gas_solid_density, ONLY: rlk
      USE gas_solid_velocity, ONLY:  us, vs, ws
      USE particles_constants, ONLY: rl, inrl
      USE set_indexes, ONLY: subscr, ijke, ijkn, ijkt
      USE turbulence_model, ONLY: must, iss
      USE indijk_module, ONLY: ip0_jp0_kp0_

      IMPLICIT NONE
!
      INTEGER :: is 
!
      REAL*8 :: a, b
      REAL*8 :: dxp, dyp, dzp, indxp, indyp, indzp
      REAL*8 :: epsx, epsy, epsz, gepx, gepy, gepz
      INTEGER :: imesh, i, j, k, ijk
      REAL*8, ALLOCATABLE :: rlkinrl(:)
!
      CALL data_exchange( mus )

!
      IF (iss >= 1) THEN
        must = must + mus
      ELSE IF (iss == 0) THEN
        must = mus
      END IF

!
! ... Newtonian stress tensor
!
      ALLOCATE( rlkinrl( SIZE( rlk, 1 ) ) )
      DO is = 1, nsolid
        rlkinrl = rlk(:,is)*inrl(is)
        IF (job_type == '2D' ) THEN
          CALL stress2D(pvisx(:,is), pvisz(:,is),                 &
  	                must(:,is), mus(:,is), rlkinrl, &
  		        us(:,is), ws(:,is))
        ELSE IF (job_type == '3D' ) THEN
          CALL stress3D(pvisx(:,is), pvisy(:,is), pvisz(:,is),    &
       	                must(:,is), mus(:,is), rlkinrl, &
  		        us(:,is), vs(:,is), ws(:,is))
        END IF 
      END DO
      DEALLOCATE( rlkinrl )

!      CALL hangup
!      STOP 'qui'

!
! ... Repulsive model (Gidaspow and Ettehadieh, 1983)
!
      IF ( repulsive_model == 1 ) THEN

        a = -8.76D0
        b = -0.27D0

        DO ijk = 1, ncint
          IF(flag(ijk) == 1) THEN
            CALL meshinds(ijk,imesh,i,j,k)
            CALL subscr(ijk)

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
              gepx=10.D0**( -a * epsx + b )

              pvisx(ijk,is) = pvisx(ijk,is) -  & 
                            gepx*indxp*2.D0*(rlk(ijke,is)-rlk(ijk,is))*inrl(is)
!
! ... Coulombic z-gradient
!
              epsz=(dz(k+1)*rlk(ijk,is) + dz(k)*rlk(ijkt,is)) * indzp * inrl(is)
              gepz=10.D0**( -a * epsz + b )

              pvisz(ijk,is) = pvisz(ijk,is) -  & 
                            gepz*indzp*2.D0*(rlk(ijkt,is)-rlk(ijk,is))*inrl(is)
            
              IF (job_type == '3D') THEN
!
! ... Coulombic y-gradient
!
                epsy=(dy(j+1)*rlk(ijk,is) + dy(j)*rlk(ijkn,is))*indyp * inrl(is)
                gepy=10.D0**( -a * epsy + b )

                pvisy(ijk,is) = pvisy(ijk,is) -  & 
                            gepy*indyp*2.D0*(rlk(ijkn,is)-rlk(ijk,is))*inrl(is)

              END IF

            END DO
          END IF

        END DO
      END IF
!
      RETURN
      END SUBROUTINE viscs
!----------------------------------------------------------------------
      SUBROUTINE stress3D(visx, visy, visz, mu, lambda, eps, u, v, w)
! ... This routine computes the components of the 3D viscous stress tensor
!
! ... 'mu' and 'lambda' are the first and second viscosity coefficients. 
! ... In general, they can be different.
!
      USE dimensions
      USE domain_decomposition, ONLY: ncint
      USE grid, ONLY: itc, flag 
      USE grid, ONLY: dx, dy, dz, indx, indy, indz
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
      REAL*8 :: indxm2, indxp2, indym2, indyp2, indzm2, indzp2
      REAL*8 :: cm
      REAL*8 :: dum, dup, dvm, dvp, dwm, dwp, dvp1
      REAL*8 :: dupm1, dvpm1, dwpm1, dupm2, dvpm2, dwpm2
      REAL*8 :: epsc, epse, epsw, epsn, epss, epst, epsb
      REAL*8 :: epsen, epset, epses, epseb, epswn, epswt, epsnt, epsnb, epsst 
      REAL*8 :: muc, mue, muw, mun, mus, mut, mub
      REAL*8 :: muen, muet, mues, mueb, muwn, muwt, munt, munb, must 
       
      INTEGER :: imesh, i, j, k, ijk
!
      cm = 2.D0/3.D0

      visx = 0.D0
      visy = 0.D0
      visz = 0.D0
!
      cells: DO ijk = 1, ncint

        imesh = myijk( ip0_jp0_kp0_ , ijk )

        IF(flag(ijk) == 1) THEN

         CALL subscr(ijk)
         !CALL subscr_iter(ijk)

         i = MOD( MOD( imesh - 1, nx*ny ), nx ) + 1
	 j = MOD( imesh - 1, nx*ny ) / nx + 1
	 k = ( imesh - 1 ) / ( nx*ny ) + 1
!
         !IF( i < 2 )           WRITE(*,*) 'error 1 ', ijk
         !IF( i == SIZE( dx ) ) WRITE(*,*) 'error 2 ', ijk
         !IF( j < 2 )           WRITE(*,*) 'error 3 ', ijk
         !IF( j == SIZE( dy ) ) WRITE(*,*) 'error 4 ', ijk
         !IF( k < 2 )           WRITE(*,*) 'error 5 ', ijk
         !IF( k == SIZE( dz ) ) WRITE(*,*) 'error 6 ', ijk

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

         indxp2 = indxp * 2.0D0
         indxm2 = indxm * 2.0D0
         indyp2 = indyp * 2.0D0
         indym2 = indym * 2.0D0
         indzp2 = indzp * 2.0D0
         indzm2 = indzm * 2.0D0

         !IF( ijk < 1  .OR. ijk  > SIZE(u) ) WRITE(*,*) 'error 7 ',ijk
         !IF( imjk < 1 .OR. imjk > SIZE(u) ) WRITE(*,*) 'error 8 ',imjk
         !IF( ipjk < 1 .OR. ipjk > SIZE(u) ) WRITE(*,*) 'error 9 ',ipjk
         !IF( ijk < 1  .OR. ijk  > SIZE(v) ) WRITE(*,*) 'error 10',ijk
         !IF( ijmk < 1 .OR. ijmk > SIZE(v) ) WRITE(*,*) 'error 11',ijmk
         !IF( ijpk < 1 .OR. ijpk > SIZE(v) ) WRITE(*,*) 'error 12',ijpk
         !IF( ijk < 1  .OR. ijk  > SIZE(w) ) WRITE(*,*) 'error 13',ijk
         !IF( ijkm < 1 .OR. ijkm > SIZE(w) ) WRITE(*,*) 'error 14',ijkm
         !IF( ijkp < 1 .OR. ijkp > SIZE(w) ) WRITE(*,*) 'error 15',ijkp

         dum = u(ijk) - u(imjk)
         dup = u(ipjk) - u(ijk)
         dvm = v(ijk) - v(ijmk)
         dvp = v(ijpk) - v(ijk)    
         dwm = w(ijk) - w(ijkm)
         dwp = w(ijkp) - w(ijk)

         !IF( ijpk  < 1 .OR. ijpk  > SIZE(u) ) WRITE(*,*) 'error 16',ijpk
         !IF( imjpk < 1 .OR. imjpk > SIZE(u) ) WRITE(*,*) 'error 17',imjpk
         !IF( ijkp  < 1 .OR. ijkp  > SIZE(u) ) WRITE(*,*) 'error 18',ijkp 
         !IF( imjkp < 1 .OR. imjkp > SIZE(u) ) WRITE(*,*) 'error 19',imjkp
         !IF( ipjk  < 1 .OR. ipjk  > SIZE(v) ) WRITE(*,*) 'error 20',ipjk 
         !IF( ipjmk < 1 .OR. ipjmk > SIZE(v) ) WRITE(*,*) 'error 21',ipjmk
         !IF( ijkp  < 1 .OR. ijkp  > SIZE(v) ) WRITE(*,*) 'error 22',ijkp 
         !IF( ijmkp < 1 .OR. ijmkp > SIZE(v) ) WRITE(*,*) 'error 23',ijmkp 
         !IF( ipjk  < 1 .OR. ipjk  > SIZE(w) ) WRITE(*,*) 'error 24',ipjk  
         !IF( ipjkm < 1 .OR. ipjkm > SIZE(w) ) WRITE(*,*) 'error 25',ipjkm 
         !IF( ijpk  < 1 .OR. ijpk  > SIZE(w) ) WRITE(*,*) 'error 26',ijpk  
         !IF( ijpkm < 1 .OR. ijpkm > SIZE(w) ) WRITE(*,*) 'error 27',ijpkm  

         dupm1 = u(ijpk) - u(imjpk)
         dupm2 = u(ijkp) - u(imjkp)
         dvpm1 = v(ipjk) - v(ipjmk)
         dvpm2 = v(ijkp) - v(ijmkp)
         dwpm1 = w(ipjk) - w(ipjkm)
         dwpm2 = w(ijpk) - w(ijpkm)

!
! ... divergence of the velocity field at centered, east, north and top cells
!
         divc = dum * indx(i) + dvm * indy(j) + dwm * indz(k) 
         dive = dup * indx(i+1) + dvpm1 * indy(j) + dwpm1 * indz(k) 
         divn = dupm1 * indx(i) + dvp * indy(j+1) + dwpm2 * indz(k) 
         divt = dupm2 * indx(i) + dvpm2* indy(j) + dwp * indz(k+1) 

!         
! ... diagonal components of the strain tensor ...
!
         txx2 = 2.D0 * dup * indx(i+1)  
         txx1 = 2.D0 * dum * indx(i)   
         tyy2 = 2.D0 * dvp * indy(j+1) 
         tyy1 = 2.D0 * dvm  * indy(j) 
         tzz2 = 2.D0 * dwp * indz(k+1) 
         tzz1 = 2.D0 * dwm * indz(k) 
!
! ... (isotropy)
!
         !IF( ijke < 1 .OR. ijke > SIZE(mu) ) WRITE(*,*) 'error 28',ijke  
         !IF( ijk  < 1 .OR. ijk  > SIZE(mu) ) WRITE(*,*) 'error 29',ijk  
         !IF( ijkn < 1 .OR. ijkn > SIZE(mu) ) WRITE(*,*) 'error 30',ijkn  
         !IF( ijkt < 1 .OR. ijkt > SIZE(mu) ) WRITE(*,*) 'error 31',ijkn  
         !IF( ijke < 1 .OR. ijke > SIZE(lambda) ) WRITE(*,*) 'error 32',ijke  
         !IF( ijk  < 1 .OR. ijk  > SIZE(lambda) ) WRITE(*,*) 'error 33',ijk  
         !IF( ijkn < 1 .OR. ijkn > SIZE(lambda) ) WRITE(*,*) 'error 34',ijkn  
         !IF( ijkt < 1 .OR. ijkt > SIZE(lambda) ) WRITE(*,*) 'error 35',ijkn  

         txx2 = mu(ijke) * txx2 - cm * lambda(ijke) * dive
         txx1 = mu(ijk)  * txx1 - cm * lambda(ijk)  * divc
         tyy2 = mu(ijkn) * tyy2 - cm * lambda(ijkn) * divn
         tyy1 = mu(ijk)  * tyy1 - cm * lambda(ijk)  * divc
         tzz2 = mu(ijkt) * tzz2 - cm * lambda(ijkt) * divt
         tzz1 = mu(ijk)  * tzz1 - cm * lambda(ijk)  * divc
!
! ... non-diagonal component of the stress tensor
!
         !IF( ijmk < 1 .OR. ijmk > SIZE(u) ) WRITE(*,*) 'error 36',ijmk  
         !IF( ijpk < 1 .OR. ijpk > SIZE(u) ) WRITE(*,*) 'error 37',ijpk  
         !IF( imjk < 1 .OR. imjk > SIZE(v) ) WRITE(*,*) 'error 38',imjk  
         !IF( ipjk < 1 .OR. ipjk > SIZE(v) ) WRITE(*,*) 'error 39',ipjk  
         !IF( ijmk < 1 .OR. ijmk > SIZE(w) ) WRITE(*,*) 'error 40',ijmk  
         !IF( ijpk < 1 .OR. ijpk > SIZE(w) ) WRITE(*,*) 'error 41',ijpk  

         dum = u(ijk) - u(ijmk)
         dup = u(ijpk) - u(ijk)
         dvm = v(ijk) - v(imjk)
         dvp = v(ipjk) - v(ijk)
         dwm = w(ijk) - w(ijmk)
         dwp = w(ijpk) - w(ijk)

         !IF( imjpk < 1 .OR. imjpk > SIZE(u) ) WRITE(*,*) 'error 42',imjpk  
         !IF( imjkp < 1 .OR. imjkp > SIZE(u) ) WRITE(*,*) 'error 43',imjkp  
         !IF( ipjmk < 1 .OR. ipjmk > SIZE(v) ) WRITE(*,*) 'error 44',ipjmk  
         !IF( ijmkp < 1 .OR. ijmkp > SIZE(v) ) WRITE(*,*) 'error 45',ijmkp  
         !IF( ijpkm < 1 .OR. ijpkm > SIZE(w) ) WRITE(*,*) 'error 46',ijpkm  
         !IF( ipjkm < 1 .OR. ipjkm > SIZE(w) ) WRITE(*,*) 'error 47',ipjkm  

         dupm1 = u(imjpk) - u(imjk)
         dupm2 = u(imjkp) - u(imjk)
         dvpm1 = v(ipjmk) - v(ijmk)
         dvpm2 = v(ijmkp) - v(ijmk) 
         dwpm1 = w(ijpkm) - w(ijkm)
         dwpm2 = w(ipjkm) - w(ijkm)

         tyx2 = dup * indyp2 + dvp * indxp2
         tyx1 = dum * indym2 + dvpm1 * indxp2
         txy2 = tyx2
         txy1 = dupm1 * indyp2 + dvm * indxm2

         !IF( ijkp < 1 .OR. ijkp > SIZE(v) ) WRITE(*,*) 'error 48',ijkp  
         !IF( ijkm < 1 .OR. ijkm > SIZE(v) ) WRITE(*,*) 'error 49',ijkm  

         dvp = v(ijkp) - v(ijk)
         dvm = v(ijk) - v(ijkm)

         tzy2 = dvp * indzp2 + dwp * indyp2
         tzy1 = dvm * indzm2 + dwpm1 * indyp2
         tyz2 = tzy2
         tyz1 = dvpm2 * indzp2 + dwm * indym2

         !IF( ipjk < 1 .OR. ipjk > SIZE(w) ) WRITE(*,*) 'error 50',ipjk  
         !IF( imjk < 1 .OR. imjk > SIZE(w) ) WRITE(*,*) 'error 51',imjk  
         !IF( ijkp < 1 .OR. ijkp > SIZE(u) ) WRITE(*,*) 'error 52',ijkp  
         !IF( ijkm < 1 .OR. ijkm > SIZE(u) ) WRITE(*,*) 'error 53',ijkm  

         dwp = w(ipjk) - w(ijk)
         dwm = w(ijk) - w(imjk)
         dup = u(ijkp) - u(ijk)
         dum = u(ijk) - u(ijkm)

         txz2 = dwp * indxp2 + dup * indzp2
         txz1 = dwm * indxm2 + dupm2 * indzp2
         tzx2 = txz2
         tzx1 = dwpm2 * indxp2 + dum * indzm2

!
! ... compute linearly interpolated values of viscosity on the staggered grid
!
         !IF( ijke < 1 .OR. ijke > SIZE(eps) ) WRITE(*,*) 'error 54',ijke  
         !IF( ijkw < 1 .OR. ijkw > SIZE(eps) ) WRITE(*,*) 'error 55',ijkw  
         !IF( ijkn < 1 .OR. ijkn > SIZE(eps) ) WRITE(*,*) 'error 56',ijkn  
         !IF( ijks < 1 .OR. ijks > SIZE(eps) ) WRITE(*,*) 'error 57',ijks  
         !IF( ijkt < 1 .OR. ijkt > SIZE(eps) ) WRITE(*,*) 'error 58',ijkt  
         !IF( ijkb < 1 .OR. ijkb > SIZE(eps) ) WRITE(*,*) 'error 59',ijkb  
         !IF( ijken < 1 .OR. ijken > SIZE(eps) ) WRITE(*,*) 'error 60',ijken  
         !IF( ijket < 1 .OR. ijket > SIZE(eps) ) WRITE(*,*) 'error 61',ijket  
         !IF( ijkes < 1 .OR. ijkes > SIZE(eps) ) WRITE(*,*) 'error 62',ijkes  
         !IF( ijkeb < 1 .OR. ijkeb > SIZE(eps) ) WRITE(*,*) 'error 63',ijkeb  
         !IF( ijkwn < 1 .OR. ijkwn > SIZE(eps) ) WRITE(*,*) 'error 64',ijkwn  
         !IF( ijkwt < 1 .OR. ijkwt > SIZE(eps) ) WRITE(*,*) 'error 65',ijkwt  
         !IF( ijknt < 1 .OR. ijknt > SIZE(eps) ) WRITE(*,*) 'error 66',ijknt  
         !IF( ijknb < 1 .OR. ijknb > SIZE(eps) ) WRITE(*,*) 'error 67',ijknb  
         !IF( ijkst < 1 .OR. ijkst > SIZE(eps) ) WRITE(*,*) 'error 68',ijkst  

         epsc = eps(ijk)
         epse = eps(ijke)
         epsw = eps(ijkw)
         epsn = eps(ijkn)
         epss = eps(ijks)
         epst = eps(ijkt)
         epsb = eps(ijkb)
         epsen = eps(ijken)
         epset = eps(ijket)
         epses = eps(ijkes)
         epseb = eps(ijkeb)
         epswn = eps(ijkwn)
         epswt = eps(ijkwt)
         epsnt = eps(ijknt)
         epsnb = eps(ijknb)
         epsst = eps(ijkst)

         !IF( ijke < 1 .OR. ijke > SIZE(mu) ) WRITE(*,*) 'error 69',ijke  
         !IF( ijkw < 1 .OR. ijkw > SIZE(mu) ) WRITE(*,*) 'error 70',ijkw  
         !IF( ijkn < 1 .OR. ijkn > SIZE(mu) ) WRITE(*,*) 'error 71',ijkn  
         !IF( ijks < 1 .OR. ijks > SIZE(mu) ) WRITE(*,*) 'error 72',ijks  
         !IF( ijkt < 1 .OR. ijkt > SIZE(mu) ) WRITE(*,*) 'error 73',ijkt  
         !IF( ijkb < 1 .OR. ijkb > SIZE(mu) ) WRITE(*,*) 'error 74',ijkb  
         !IF( ijken < 1 .OR. ijken > SIZE(mu) ) WRITE(*,*) 'error 75',ijken  
         !IF( ijket < 1 .OR. ijket > SIZE(mu) ) WRITE(*,*) 'error 76',ijket  
         !IF( ijkes < 1 .OR. ijkes > SIZE(mu) ) WRITE(*,*) 'error 77',ijkes  
         !IF( ijkeb < 1 .OR. ijkeb > SIZE(mu) ) WRITE(*,*) 'error 78',ijkeb  
         !IF( ijkwn < 1 .OR. ijkwn > SIZE(mu) ) WRITE(*,*) 'error 79',ijkwn  
         !IF( ijkwt < 1 .OR. ijkwt > SIZE(mu) ) WRITE(*,*) 'error 80',ijkwt  
         !IF( ijknt < 1 .OR. ijknt > SIZE(mu) ) WRITE(*,*) 'error 81',ijknt  
         !IF( ijknb < 1 .OR. ijknb > SIZE(mu) ) WRITE(*,*) 'error 82',ijknb  
         !IF( ijkst < 1 .OR. ijkst > SIZE(mu) ) WRITE(*,*) 'error 83',ijkst  

         muc = mu(ijk)
         mue = mu(ijke)
         muw = mu(ijkw)
         mun = mu(ijkn)
         mus = mu(ijks)
         mut = mu(ijkt)
         mub = mu(ijkb)
         muen = mu(ijken)
         muet = mu(ijket)
         mues = mu(ijkes)
         mueb = mu(ijkeb)
         muwn = mu(ijkwn)
         muwt = mu(ijkwt)
         munt = mu(ijknt)
         munb = mu(ijknb)
         must = mu(ijkst)

         epm1 =   ( dy(j) * epsn * mun + dy(j+1) * epsc *muc ) * indyp
         epm2 =   ( dy(j) * epsen * muen + dy(j+1) * epse * mue ) * indyp
         epmyx2 = ( dx(i+1) * epm1 + dx(i) * epm2 ) * indxp
!
         epm1 =   ( dy(j-1) * epsc * muc + dy(j) * epss * mus ) * indym
         epm2 =   ( dy(j-1) * epse * mue + dy(j) * epses * mues ) * indym
         epmyx1 = ( dx(i+1) * epm1 + dx(i) * epm2 ) * indxp
!
         epm1 =   ( dz(k) * epst * mut + dz(k+1) * epsc * muc ) * indzp
         epm2 =   ( dz(k) * epset * muet + dz(k+1) * epse * mue ) * indzp
         epmzx2 = ( dx(i+1) * epm1 + dx(i) * epm2 ) * indxp
!
         epm1 =   ( dz(k-1) * epsc * muc + dz(k) * epsb * mub ) * indzm
         epm2 =   ( dz(k-1) * epse * mue + dz(k) * epseb * mueb ) * indzm
         epmzx1 = ( dx(i+1) * epm1 + dx(i) * epm2 ) * indxp

! 
! ... X-gradient of the stress tensor
!
         visx(ijk) = (epse * txx2 - epsc * txx1) * indxp2 
         visx(ijk) = visx(ijk) + (epmyx2 * tyx2 - epmyx1 * tyx1) * indy(j)    
         visx(ijk) = visx(ijk) + (epmzx2 * tzx2 - epmzx1 * tzx1) * indz(k)
!
! ... compute linearly interpolated values of viscosity on the staggered grid
!
         epm1 =   ( dy(j) * epswn * muwn + dy(j+1) * epsw * muw) * indyp
         epm2 =   ( dy(j) * epsn * mun + dy(j+1) * epsc * muc ) * indyp
         epmxy1 = ( dx(i) * epm1 + dx(i-1) * epm2 ) * indxm
!
         epmxy2 = epmyx2
!
         epm1 =   ( dz(k) * epst * mut + dz(k+1) * epsc * muc ) * indzp
         epm2 =   ( dz(k) * epsnt * munt + dz(k+1) * epsn * mun ) * indzp
         epmzy2 = ( dy(j+1) * epm1 + dy(j) * epm2 ) * indyp
!
         epm1 =   ( dz(k-1) * epsc * muc + dz(k) * epsb * mub ) * indzm
         epm2 =   ( dz(k-1) * epsn * mun + dz(k) * epsnb * munb ) * indzm
         epmzy1 = ( dy(j+1) * epm1 + dy(j) * epm2 ) * indyp

!
! ... Y-gradient of the stress tensor
!
         visy(ijk) = (epsn * tyy2 - epsc * tyy1) * indyp2    
         visy(ijk) = visy(ijk) + (epmxy2 * txy2 - epmxy1 * txy1) * indx(i)         
         visy(ijk) = visy(ijk) + (epmzy2 * tzy2 - epmzy1 * tzy1) * indz(k) 
!
! ... compute linearly interpolated values of viscosity on the staggered grid
!
         epmxz2 = epmzx2
!
         epm1 =   ( dz(k) * epswt * muwt + dz(k+1) * epsw * muw) * indzp
         epm2 =   ( dz(k) * epst * mut + dz(k+1) * epsc * muc ) * indzp
         epmxz1 = ( dx(i) * epm1 + dx(i-1) * epm2 ) * indxm
!
         epmyz2 = epmzy2
!
         epm1 =   ( dz(k) * epsst * must + dz(k+1) * epss * mus) * indzp
         epm2 =   ( dz(k) * epst * mut + dz(k+1) * epsc * muc ) * indzp
         epmyz1 = ( dy(j) * epm1 + dy(j-1) * epm2 ) * indym

!
! ... Z-gradient of the stress tensor
!
         visz(ijk) = (epst * tzz2 - epsc * tzz1) * indzp2   
         visz(ijk) = visz(ijk) + (epmxz2 * txz2 - epmxz1 * txz1) * indx(i)           
         visz(ijk) = visz(ijk) + (epmyz2 * tyz2 - epmyz1 * tyz1) * indy(j) 

         ! CYCLE cells

        END IF
      END DO cells


!
      RETURN
      END SUBROUTINE stress3D
!----------------------------------------------------------------------
      SUBROUTINE stress2D(visx, visz, mu, lambda, eps, u, w)
! ... This routine computes the components of the 
! ... 2D cylindrical and cartesian viscous stress tensor
!
! ... 'mu' and 'lambda' are the first and second viscosity coefficients. 
! ... In general, they can be different.
!
      USE dimensions
      USE domain_decomposition, ONLY: ncint
      USE grid, ONLY: itc, dz, dx, r, rb, indz, indx, inr, inrb
      USE grid, ONLY: flag
      USE set_indexes
      USE indijk_module, ONLY: ip0_jp0_kp0_

      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: visx(:), visz(:)
      REAL*8, INTENT(IN)  :: mu(:), lambda(:), eps(:), u(:), w(:)
      REAL*8 :: t0, mu0, lambda0, eps0
      REAL*8 :: txx2, txx1, tyy2, tyy1
      REAL*8 :: tyx2, tyx1, txy2, txy1
      REAL*8 :: gradxx, gradyy, gradxy, gradyx
      REAL*8 :: wm1, wm2, dw, du
      REAL*8 :: divc, divr, divt, dive
      REAL*8 :: epsmu2, epsmu11, epsmu22, epsmu1, epsmu12, epsmu21
      REAL*8 :: dxm, dxp, dzm, dzp, indxm, indxp, indzm, indzp
      REAL*8 :: cm
      INTEGER :: imesh, i, k, ijk
!
      cm = 2.D0/3.D0

      visx = 0.D0
      visz = 0.D0
!
      DO ijk = 1, ncint
        IF(flag(ijk) == 1) THEN
         imesh = myijk( ip0_jp0_kp0_, ijk)
         CALL subscr(ijk)
         k = ( imesh - 1 ) / nx + 1
         i = MOD( ( imesh - 1 ), nx) + 1
!
         dxp=(dx(i)+dx(i+1))
         dxm=(dx(i)+dx(i-1))
         dzp=(dz(k)+dz(k+1))
         dzm=(dz(k)+dz(k-1))
!
         indxp=1.D0/dxp
         indxm=1.D0/dxm
         indzp=1.D0/dzp
         indzm=1.D0/dzm
!
! ... divergence of the velocity field at right, centered, top cells
!
         divr = ((rb(i+1)*u(ipjk)-rb(i)*u(ijk))*inr(i+1)*indx(i+1)        &
              + (w(ipjk)-w(ipjkm))*indz(k))
         divc = ((rb(i)*u(ijk)-rb(i-1)*u(imjk))*inr(i)*indx(i)            &
              + (w(ijk)-w(ijkm))*indz(k))
         divt = ((rb(i)*u(ijkp)-rb(i-1)*u(imjkp))*inr(i)*indx(i)          &
              + (w(ijkp)-w(ijk))*indz(k+1))
!         
! ... diagonal components of the strain tensor ...
!
         txx2 = 2.D0*(u(ipjk)-u(ijk))*indx(i+1)  
         txx1 = 2.D0*(u(ijk)-u(imjk))*indx(i)   
         tyy2 = 2.D0*(w(ijkp)-w(ijk))*indz(k+1) 
         tyy1 = 2.D0*(w(ijk)-w(ijkm))*indz(k)   
!
! ... (isotropy)
!
         txx2 = mu(ijke) * txx2 - cm * lambda(ijke) * divr
         txx1 = mu(ijk)  * txx1 - cm * lambda(ijk)  * divc
         tyy2 = mu(ijkt) * tyy2 - cm * lambda(ijkt) * divt
         tyy1 = mu(ijk)  * tyy1 - cm * lambda(ijk)  * divc
!
! ... non-diagonal component of the stress tensor
!
         tyx2 = (u(ijkp)-u(ijk))*indzp*2.D0 + (w(ipjk)-w(ijk))*indxp*2.D0
         tyx1 = (u(ijk)-u(ijkm))*indzm*2.D0 + (w(ipjkm)-w(ijkm))*indxp*2.D0
         txy2 = tyx2
         txy1 = (u(imjkp)-u(imjk))*indzp*2.D0+(w(ijk)-w(imjk))*indxm*2.D0
! 
! ... Correction for cylindrical coordinates
!
         IF(itc == 1) THEN
           eps0=(eps(ijk)*dx(i+1)+eps(ijke)*dx(i))*indxp
           mu0=(mu(ijk)*dx(i+1)+mu(ijke)*dx(i))*indxp
           lambda0=(lambda(ijk)*dx(i+1)+lambda(ijke)*dx(i))*indxp
           t0 = 2.D0 * u(ijk) * inrb(i)
!
! ... divergence of the velocity field at cell boundary
!
           du = (rb(i+1)*u(ipjk)-rb(i-1)*u(imjk))
           wm2=(dx(i)*w(ipjk)+dx(i+1)*w(ijk))*indxp
           wm1=(dx(i)*w(ipjkm)+dx(i+1)*w(ijkm))*indxp
           dw = wm2 - wm1
           dive =  (inrb(i) * du * indxp + dw * indz(k))
!
! ... (isotropy)
!
           t0 = mu0 * t0 - cm *lambda0 * dive
         ELSE
           t0 = 0.D0
         ENDIF
!
         epsmu21=(dz(k)*eps(ijkt)*mu(ijkt)+dz(k+1)*eps(ijk)*mu(ijk))*indzp
         epsmu22=(dz(k)*eps(ijket)*mu(ijket)+dz(k+1)*eps(ijke)*mu(ijke))*indzp
         epsmu2=(dx(i+1)*epsmu21+dx(i)*epsmu22)*indxp
         epsmu11=(dz(k-1)*eps(ijk)*mu(ijk)+dz(k)*eps(ijkb)*mu(ijkb))*indzm
         epsmu12=(dz(k-1)*eps(ijke)*mu(ijke)+dz(k)*eps(ijkeb)*mu(ijkeb))*indzm
         epsmu1=(dx(i+1)*epsmu11+dx(i)*epsmu12)*indxp
!
! ... x-gradient of the stress tensor
!
         gradxx   = (eps(ijke)*txx2*r(i+1)-eps(ijk)*txx1*r(i))
         gradxx   = gradxx * indxp * 2.D0 * inrb(i)
         gradyx   = (epsmu2*tyx2-epsmu1*tyx1)
         gradyx   = gradyx * indz(k)

         visx(ijk) = gradxx + gradyx - eps0 * t0 * inrb(i)
!
         epsmu11=(dz(k)*eps(ijkwt)*mu(ijkwt)+dz(k+1)*eps(ijkw)*mu(ijkw))*indzp
         epsmu1=(dx(i)*epsmu11+dx(i-1)*epsmu21)*indxm
!
! ... z-gradient of the stress tensor
!
         gradyy   = (eps(ijkt)*tyy2 - eps(ijk)*tyy1)
         gradyy   = gradyy * indzp*2.D0
         gradxy   = (rb(i)*epsmu2*txy2-rb(i-1)*epsmu1*txy1)
         gradxy   = gradxy * inr(i)*indx(i)

         visz(ijk) = gradyy + gradxy
                    
        END IF
      END DO
!
      RETURN
      END SUBROUTINE stress2d
!----------------------------------------------------------------------
      END MODULE gas_solid_viscosity
!----------------------------------------------------------------------
