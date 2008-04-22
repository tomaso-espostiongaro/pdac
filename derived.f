!----------------------------------------------------------------------
      MODULE derived_fields

      USE dimensions
      USE gas_constants, ONLY: gammaair, rgas, gmw, gas_type
      USE particles_constants, ONLY: rl
      INTEGER :: ind, ig, is
      INTEGER :: arraysize
!
      PUBLIC :: total_particle_fraction, void_fraction, &
               log10_epstot, gas_density, gas_mass_fractions, &
               gas_bulk_density, solid_bulk_density, &
               mixture_density, mixture_velocity, &
               velocity_module_2D, velocity_module_3D, &
               dynamic_pressure 
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE gradient(array,grx,gry,grz)
      USE control_flags, ONLY: job_type
      USE domain_mapping, ONLY: meshinds, ncint
      USE set_indexes, ONLY: first_subscr
      USE set_indexes, ONLY: imjk,ipjk,ijmk,ijpk,ijkm,ijkp
      USE grid, ONLY: dx, dy, dz
      USE parallel, ONLY: mpime, root
      !
      ! ... computes the gradient components of a scalar array
      !
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: array(:)
      REAL*8, INTENT(OUT), DIMENSION(SIZE(array)) :: grx(:),gry(:),grz(:)
      REAL*8 :: delta
      INTEGER :: ijk,i,j,k,imesh
!       
      DO ijk = 1, ncint
        CALL first_subscr(ijk)
        CALL meshinds(ijk,imesh,i,j,k)
        IF (i/=1 .AND. i/=nx) THEN
          delta = dx(i)+0.5D0*(dx(i-1)+dx(i+1))
          grx(ijk) = (array(ipjk)-array(imjk))/delta
        ELSE IF (i==1) THEN
          delta = 0.5D0*(dx(i)+dx(i+1))
          grx(ijk) = (array(ipjk)-array(ijk))/delta
        ELSE IF (i==nx) THEN
          delta = 0.5D0*(dx(i)+dx(i-1))
          grx(ijk) = (array(ijk)-array(imjk))/delta
        END IF
        !
        IF (job_type == '3D') THEN
          IF (j/=1 .AND. j/=ny) THEN
            delta = dy(j)+0.5D0*(dy(j-1)+dy(j+1))
            gry(ijk) = (array(ijpk)-array(ijmk))/delta
          ELSE IF (j==1) THEN
            delta = 0.5D0*(dy(j)+dy(j+1))
            gry(ijk) = (array(ijpk)-array(ijk))/delta
          ELSE IF (j==ny) THEN
            delta = 0.5D0*(dy(j)+dy(j-1))
            gry(ijk) = (array(ijk)-array(ijmk))/delta
          END IF
        ELSE
          gry(ijk) = 0.D0
        END IF
        !
        IF (k/=1 .AND. k/=nz) THEN
          delta = dz(k)+0.5D0*(dz(k-1)+dz(k+1))
          grz(ijk) = (array(ijk)-array(ijkm))/delta
        ELSE IF (k==1) THEN
          delta = 0.5D0*(dz(k)+dz(k+1))
          grz(ijk) = (array(ijkp)-array(ijk))/delta
        ELSE IF (k==nz) THEN
          delta = 0.5D0*(dz(k)+dz(k-1))
          grz(ijk) = (array(ijk)-array(ijkm))/delta
        END IF
      END DO
      RETURN
      END SUBROUTINE gradient
!----------------------------------------------------------------------
      SUBROUTINE total_particle_fraction(epst,eps)
      !
      ! ... computes the total particle volumetric fraction
      
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: eps(:,:)
      REAL*8, INTENT(OUT), DIMENSION(SIZE(eps,DIM=1)) :: epst

      arraysize = SIZE(eps,DIM=1)
      DO ind = 1, arraysize
        epst(ind) = SUM(eps(ind,:))
      END DO

      RETURN
      END SUBROUTINE total_particle_fraction
!----------------------------------------------------------------------
      SUBROUTINE void_fraction(ep,epst)
      !
      ! ... computes the void fraction

      IMPLICIT NONE
      REAL*8, INTENT(IN) :: epst(:)
      REAL*8, INTENT(OUT), DIMENSION(SIZE(epst)) :: ep

      arraysize = SIZE(epst)
      DO ind = 1, arraysize
        ep(ind) = 1.D0 - epst(ind)
      END DO

      RETURN
      END SUBROUTINE void_fraction
!----------------------------------------------------------------------
      SUBROUTINE log10_epstot(lepst,epst)
      !
      ! ... computes the log10 of the particle volumetric fraction
      
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: epst(:)
      REAL*8, INTENT(OUT), DIMENSION(SIZE(epst)) :: lepst
      REAL*8 :: clamp_epst

      arraysize = SIZE(epst)
      DO ind = 1, arraysize
        clamp_epst = MAX(1.D-10,epst(ind))
        lepst(ind) = log10(clamp_epst) 
      END DO

      RETURN
      END SUBROUTINE log10_epstot
!----------------------------------------------------------------------
      SUBROUTINE gas_molecular_weight(mg,xgc)
      !
      ! ... compute the averaged gas molecular weight

      IMPLICIT NONE
      REAL*8, INTENT(IN), DIMENSION(:,:) :: xgc
      REAL*8, INTENT(OUT), DIMENSION(SIZE(xgc,DIM=1)) :: mg

      mg = 0.D0
      arraysize = SIZE(xgc,DIM=1)
      DO ind = 1, arraysize
        DO ig = 1, ngas
          mg(ind) = mg(ind) + xgc(ind,ig) * gmw(gas_type(ig))
        END DO
      END DO

      RETURN
      END SUBROUTINE gas_molecular_weight
!----------------------------------------------------------------------
      SUBROUTINE gas_density(rhog,p,tg,xgc,mgas)
      !
      ! ... computes gas thermodynamic density

      IMPLICIT NONE
      REAL*8, INTENT(IN) :: p(:), tg(:), xgc(:,:) 
      REAL*8, INTENT(OUT), DIMENSION(SIZE(p)) :: rhog
      REAL*8, DIMENSION(SIZE(p)) :: mgas

      arraysize = SIZE(p)
      CALL gas_molecular_weight(mgas,xgc)
      DO ind = 1, arraysize
        rhog(ind) = p(ind) / ( rgas * tg(ind) ) * mgas(ind)
      END DO
      
      RETURN
      END SUBROUTINE gas_density
!----------------------------------------------------------------------
      SUBROUTINE gas_mass_fractions(ygc,xgc)
      !
      ! ... computes gas components mass fractions
      
      IMPLICIT NONE
      REAL*8, INTENT(IN), DIMENSION(:,:) :: xgc
      REAL*8, INTENT(OUT), DIMENSION(SIZE(xgc,1),SIZE(xgc,2)) :: ygc
      REAL*8, DIMENSION(SIZE(xgc,DIM=1)) :: mgas
      INTEGER :: ig

      CALL gas_molecular_weight(mgas,xgc)
      arraysize = SIZE(xgc,DIM=1)
      DO ind = 1, arraysize
        DO ig=1,ngas
          ygc(ind,ig) = xgc(ind,ig) * gmw(gas_type(ig)) / mgas(ind)
        END DO
      END DO

      RETURN
      END SUBROUTINE gas_mass_fractions
!----------------------------------------------------------------------
      SUBROUTINE gas_bulk_density(rgp,ep,rhog)
      !
      ! ... computes gas bulk density

      IMPLICIT NONE
      REAL*8, INTENT(IN) :: ep(:), rhog(:)
      REAL*8, INTENT(OUT), DIMENSION(SIZE(ep)) :: rgp

      arraysize = SIZE(ep)
      DO ind = 1, arraysize
        rgp(ind) = ep(ind) * rhog(ind)
      END DO

      END SUBROUTINE gas_bulk_density
!----------------------------------------------------------------------
      SUBROUTINE solid_bulk_density(rlk,eps)
      !
      ! ... computes solid bulk densities

      IMPLICIT NONE
      REAL*8, INTENT(IN) :: eps(:,:)
      REAL*8, DIMENSION(SIZE(eps,1),SIZE(eps,2)) :: rlk

      arraysize = SIZE(eps,DIM=1)
      DO ind = 1, arraysize
        DO is = 1, nsolid
          rlk(ind,is) = eps(ind,is) * rl(is)
        END DO
      END DO

      END SUBROUTINE solid_bulk_density
!----------------------------------------------------------------------
      SUBROUTINE mixture_density(rhom,rlk,rgp)
      !
      ! ... computes gas-particle mixture density

      IMPLICIT NONE
      REAL*8, INTENT(IN) :: rlk(:,:)
      REAL*8, INTENT(IN) :: rgp(:)
      REAL*8, DIMENSION(SIZE(rgp)) :: rhom

      arraysize = SIZE(rgp)
      DO ind = 1, arraysize
        rhom(ind) = rgp(ind) + SUM(rlk(ind,:)) 
      END DO

      RETURN
      END SUBROUTINE mixture_density
!----------------------------------------------------------------------
      SUBROUTINE mixture_temperature(tm,ts,tg,rlk,rgp,ygc)
      USE gas_constants, ONLY: gas_type
      USE particles_constants, ONLY: cps
      USE specific_heat_module, ONLY: hcapg
      !
      ! ... computes gas-particle mixture temperature

      IMPLICIT NONE
      REAL*8, INTENT(IN) :: rlk(:,:), ts(:,:), ygc(:,:)
      REAL*8, INTENT(IN) :: rgp(:), tg(:)
      REAL*8, DIMENSION(SIZE(rgp)) :: tm
      REAL*8 :: den, cgas, cpgc(7)

      arraysize = SIZE(rgp)
      DO ind = 1, arraysize
        CALL hcapg( cpgc(:), tg(ind) )
        cgas = 0.D0
        DO ig = 1, ngas
          cgas = cpgc( gas_type(ig) ) * ygc(ind,ig) + cgas
        END DO
        tm(ind) = rgp(ind)*cgas*tg(ind) + &
                  SUM(rlk(ind,:)*cps(:)*ts(ind,:)) 
        den = rgp(ind)*cgas + SUM(rlk(ind,:)*cps(:))
        tm(ind) = tm(ind) / den
      END DO

      RETURN
      END SUBROUTINE mixture_temperature
!----------------------------------------------------------------------
      SUBROUTINE particle_density(rhos,rlk)
      !
      ! ... computes particle mixture density

      IMPLICIT NONE
      REAL*8, INTENT(IN) :: rlk(:,:)
      REAL*8, DIMENSION(SIZE(rlk,1)) :: rhos

      arraysize = SIZE(rlk,1)
      DO ind = 1, arraysize
        rhos(ind) = SUM(rlk(ind,:)) 
      END DO

      RETURN
      END SUBROUTINE particle_density
!----------------------------------------------------------------------
      SUBROUTINE mixture_velocity(velm,ug,us,rlk,rgp,rhom)
      !
      ! ... computes gas-particle mixture velocity (one component)

      IMPLICIT NONE
      REAL*8, INTENT(IN), DIMENSION(:,:) :: rlk, us
      REAL*8, INTENT(IN), DIMENSION(:) :: ug, rgp, rhom
      REAL*8, DIMENSION(SIZE(ug)) :: velm
      
      arraysize = SIZE(ug)
      DO ind = 1, arraysize
        velm(ind) = rgp(ind)*ug(ind) + SUM(rlk(ind,:)*us(ind,:)) 
        !velm(ind) = SUM(rlk(ind,:)*us(ind,:)) 
        velm(ind) = velm(ind) / rhom(ind)
      END DO

      RETURN
      END SUBROUTINE mixture_velocity
!----------------------------------------------------------------------
      SUBROUTINE dynamic_pressure(pdyn,rm,velom)
      !
      ! ... computes the mixture dynamic pressure
      
      IMPLICIT NONE
      REAL*8, INTENT(IN), DIMENSION(:) :: rm, velom
      REAL*8, DIMENSION(SIZE(rm)) :: pdyn

      arraysize = SIZE(rm)
      DO ind = 1, arraysize
        pdyn(ind) = 0.5D0 * rm(ind) * velom(ind)**2
      END DO

      RETURN
      END SUBROUTINE dynamic_pressure
!----------------------------------------------------------------------
      SUBROUTINE velocity_module_2D(vel2,v1,v2)
      !
      ! ... compute the modulus of a 2D velocity field
      IMPLICIT NONE
      REAL*8, INTENT(IN), DIMENSION(:) :: v1, v2
      REAL*8, DIMENSION(SIZE(v1)) :: vel2

      arraysize = SIZE(v1)
      DO ind = 1, arraysize
        vel2(ind) = SQRT( v1(ind)**2 + v2(ind)**2 )
      END DO
      
      RETURN
      END SUBROUTINE velocity_module_2D
!----------------------------------------------------------------------
      SUBROUTINE velocity_module_3D(vel3,v1,v2,v3)
      !
      ! ... compute the modulus of a velocity field
      IMPLICIT NONE
      REAL*8, INTENT(IN), DIMENSION(:) :: v1, v2, v3
      REAL*8, DIMENSION(SIZE(v1)) :: vel3

      arraysize = SIZE(v1)
      DO ind = 1, arraysize
        vel3(ind) = SQRT( v1(ind)**2 + v2(ind)**2 + v3(ind)**2 )
      END DO
      
      RETURN
      END SUBROUTINE velocity_module_3D
!----------------------------------------------------------------------
      SUBROUTINE mixture_sound_speed_1(cm,xgc,rgp,rlk,rhom,rhog,epst,tg)
      ! 
      ! ... computes the inverse of the mixture sound speed

      IMPLICIT NONE
      REAL*8, INTENT(IN), DIMENSION(:,:) :: rlk, xgc
      REAL*8, INTENT(IN), DIMENSION(:) :: rgp,rhom,rhog,epst,tg
      REAL*8, DIMENSION(SIZE(tg)) :: cm
      REAL*8, DIMENSION(SIZE(rgp)) :: mgas
      REAL*8 :: suminv, csg2, css2, epsg

      arraysize = SIZE(rgp)
!
! ... Mixture sound speed (Wallis, 1969 - Equations 2.50 and 6.110; 
! ... Gidaspow, 1994 - Equation 4.80)
!
      CALL gas_molecular_weight(mgas,xgc)
      DO ind = 1, arraysize
        IF (epst(ind) /= 0.D0) THEN
          csg2 = gammaair * rgas * tg(ind) / mgas(ind)
          css2 = 3400.D0 ! Speed of sound in rocks [m/s] 
          epsg = 1.D0 - epst(ind)
          suminv = epsg / (rhog(ind) * csg2) + SUM(rlk(ind,:)/rl(:)**2/css2)
          cm(ind) = DSQRT(1.D0 / suminv / rhom(ind))
        ELSE
          cm(ind) = SQRT(gammaair * rgas * tg(ind) / mgas(ind) )
        END IF
      END DO
!
      RETURN
      END SUBROUTINE mixture_sound_speed_1
!----------------------------------------------------------------------
      SUBROUTINE mixture_sound_speed_2(cm,xgc,rgp,rlk,rhom,rhog,epst,tg)
      ! 
      ! ... computes the inverse of the mixture sound speed

      IMPLICIT NONE
      REAL*8, INTENT(IN), DIMENSION(:,:) :: rlk, xgc
      REAL*8, INTENT(IN), DIMENSION(:) :: rgp,rhom,rhog,epst,tg
      REAL*8, DIMENSION(SIZE(tg)) :: cm
      REAL*8, DIMENSION(SIZE(rgp)) :: mgas
      REAL*8 :: fact, y, avrl

      arraysize = SIZE(rgp)
!
! ... Mixture sound speed (Wallis, 1969 - Example 6.7)
! ... (Gidaspow, 1994 - Eq. 7.41)
!
      CALL gas_molecular_weight(mgas,xgc)
      DO ind = 1, arraysize
        IF (epst(ind) /= 0.D0) THEN
          y = rgp(ind) / rhom(ind)
          avrl = SUM(rlk(ind,:))/epst(ind)
          fact = y + (1.0 - y) * rhog(ind) / avrl
          ! ... Adiabatic gas transformation
          !cm(ind) = SQRT(gammaair * rgas * tg(ind) / mgas(ind) / y )
          ! ... Isothermal gas transformation
          cm(ind) = SQRT(1.D0 * rgas * tg(ind) / mgas(ind) / y )
          cm(ind) = cm(ind) * fact
        ELSE
          cm(ind) = SQRT(gammaair * rgas * tg(ind) / mgas(ind) )
        END IF
      END DO
!
      RETURN
      END SUBROUTINE mixture_sound_speed_2
!----------------------------------------------------------------------
      SUBROUTINE mixture_sound_speed_3(cm,xgc,rgp,rlk,rhom,rhog,epst,tg)
      ! 
      ! ... computes the inverse of the mixture sound speed

      IMPLICIT NONE
      REAL*8, INTENT(IN), DIMENSION(:,:) :: rlk, xgc
      REAL*8, INTENT(IN), DIMENSION(:) :: rgp,rhom,rhog,epst,tg
      REAL*8, DIMENSION(SIZE(tg)) :: cm
      REAL*8, DIMENSION(SIZE(rgp)) :: mgas
      REAL*8 :: fact, m, epsg, csg, gammamix, cs, rtilde

      cs = 1.2D3
      arraysize = SIZE(rgp)
!
! ... Mixture sound speed (Wallis, 1969 - Example 2.2; Kieffer, 1981)
!
      CALL gas_molecular_weight(mgas,xgc)
      DO ind = 1, arraysize
        IF (epst(ind) /= 0.D0) THEN
          m = SUM(rlk(ind,:))/rgp(ind)
          rtilde = rgas / mgas(ind)
          gammamix = (3.5D0 * rgas + m * cs)/(2.5D0 * rgas + m * cs)
          cm(ind) = SQRT(gammamix * rtilde / (1.D0 + m) * tg(ind))
        ELSE
          cm(ind) = SQRT(gammaair * rtilde * tg(ind) )
        END IF
      END DO
!
      RETURN
      END SUBROUTINE mixture_sound_speed_3
!----------------------------------------------------------------------
      SUBROUTINE mach_number(mn,vel, c)
      !
      ! ... compute the Mach number for the mixture

      IMPLICIT NONE
      REAL*8, INTENT(IN), DIMENSION(:) :: vel, c
      REAL*8, DIMENSION(SIZE(vel)) :: mn

      arraysize = SIZE(vel)
      DO ind = 1, arraysize
        mn(ind) = vel(ind) / c(ind)
      END DO

      RETURN
      END SUBROUTINE mach_number
!----------------------------------------------------------------------
      SUBROUTINE normal_mach_number(mnn,v1,v2,v3,p1,p2,p3,c)
      !
      ! ... compute the Mach number for the mixture

      IMPLICIT NONE
      REAL*8, INTENT(IN), DIMENSION(:) :: v1,v2,v3,p1,p2,p3,c
      REAL*8, DIMENSION(SIZE(c)) :: mnn
      REAL*8 :: pm

      arraysize = SIZE(c)
      DO ind = 1, arraysize
        pm = DSQRT(p1(ind)**2+p2(ind)**2+p3(ind)**2)
        mnn(ind)=(v1(ind)*p1(ind)+v2(ind)*p2(ind)+v3(ind)*p3(ind))/c(ind)/pm
      END DO

      RETURN
      END SUBROUTINE normal_mach_number
!----------------------------------------------------------------------
      END MODULE derived_fields
!----------------------------------------------------------------------
