!----------------------------------------------------------------------
      MODULE derived_fields

      USE dimensions
      USE gas_constants, ONLY: gammaair, rgas, gmw, gas_type
      USE particles_constants, ONLY: rl
      INTEGER :: imesh, ig, is
      INTEGER :: meshsize
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
      SUBROUTINE total_particle_fraction(epst,eps)
      !
      ! ... computes the total particle volumetric fraction
      
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: eps(:,:)
      REAL*8, INTENT(OUT), DIMENSION(SIZE(eps,DIM=1)) :: epst

      meshsize = SIZE(eps,DIM=1)
      DO imesh = 1, meshsize
        epst(imesh) = SUM(eps(imesh,:))
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

      meshsize = SIZE(epst)
      DO imesh = 1, meshsize
        ep(imesh) = 1.D0 - epst(imesh)
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

      meshsize = SIZE(epst)
      DO imesh = 1, meshsize
        clamp_epst = MAX(1.D-10,epst(imesh))
        lepst(imesh) = log10(clamp_epst) 
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
      meshsize = SIZE(xgc,DIM=1)
      DO imesh = 1, meshsize
        DO ig = 1, ngas
          mg(imesh) = mg(imesh) + xgc(imesh,ig) * gmw(gas_type(ig))
        END DO
      END DO

      RETURN
      END SUBROUTINE gas_molecular_weight
!----------------------------------------------------------------------
      SUBROUTINE gas_density(rhog,p,tg,xgc)
      !
      ! ... computes gas thermodynamic density

      IMPLICIT NONE
      REAL*8, INTENT(IN) :: p(:), tg(:), xgc(:,:) 
      REAL*8, INTENT(OUT), DIMENSION(SIZE(p)) :: rhog
      REAL*8, DIMENSION(SIZE(p)) :: mgas

      meshsize = SIZE(p)
      CALL gas_molecular_weight(mgas,xgc)
      DO imesh = 1, meshsize
        rhog(imesh) = p(imesh) / ( rgas * tg(imesh) ) * mgas(imesh)
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
      meshsize = SIZE(xgc,DIM=1)
      DO imesh = 1, meshsize
        DO ig=1,ngas
          ygc(imesh,ig) = xgc(imesh,ig) * gmw(gas_type(ig)) / mgas(imesh)
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

      meshsize = SIZE(ep)
      DO imesh = 1, meshsize
        rgp(imesh) = ep(imesh) * rhog(imesh)
      END DO

      END SUBROUTINE gas_bulk_density
!----------------------------------------------------------------------
      SUBROUTINE solid_bulk_density(rlk,eps)
      !
      ! ... computes solid bulk densities

      IMPLICIT NONE
      REAL*8, INTENT(IN) :: eps(:,:)
      REAL*8, DIMENSION(SIZE(eps,1),SIZE(eps,2)) :: rlk

      meshsize = SIZE(eps,DIM=1)
      DO imesh = 1, meshsize
        DO is = 1, nsolid
          rlk(imesh,is) = eps(imesh,is) * rl(is)
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

      meshsize = SIZE(rgp)
      DO imesh = 1, meshsize
        rhom(imesh) = rgp(imesh) + SUM(rlk(imesh,:)) 
      END DO

      RETURN
      END SUBROUTINE mixture_density
!----------------------------------------------------------------------
      SUBROUTINE mixture_temperature(tm,ts,tg,rlk,rgp,ygc)
      USE gas_constants, ONLY: gas_type
      USE particles_constants, ONLY: cps
      USE specific_heat_module, ONLY: hcapg
      !
      ! ... computes gas-particle mixture density

      IMPLICIT NONE
      REAL*8, INTENT(IN) :: rlk(:,:), ts(:,:), ygc(:,:)
      REAL*8, INTENT(IN) :: rgp(:), tg(:)
      REAL*8, DIMENSION(SIZE(rgp)) :: tm
      REAL*8 :: den, cgas, cpgc(7)

      meshsize = SIZE(rgp)
      DO imesh = 1, meshsize
        CALL hcapg( cpgc(:), tg(imesh) )
        cgas = 0.D0
        DO ig = 1, ngas
          cgas = cpgc( gas_type(ig) ) * ygc(imesh,ig) + cgas
        END DO
        tm(imesh) = rgp(imesh)*cgas*tg(imesh) + &
                  SUM(rlk(imesh,:)*cps(:)*ts(imesh,:)) 
        den = rgp(imesh)*cgas + SUM(rlk(imesh,:)*cps(:))
        tm(imesh) = tm(imesh) / den
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

      meshsize = SIZE(rlk,1)
      DO imesh = 1, meshsize
        rhos(imesh) = SUM(rlk(imesh,:)) 
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
      
      meshsize = SIZE(ug)
      DO imesh = 1, meshsize
        velm(imesh) = rgp(imesh)*ug(imesh) + SUM(rlk(imesh,:)*us(imesh,:)) 
        velm(imesh) = velm(imesh) / rhom(imesh)
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

      meshsize = SIZE(rm)
      DO imesh = 1, meshsize
        pdyn(imesh) = 0.5D0 * rm(imesh) * velom(imesh)**2
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

      meshsize = SIZE(v1)
      DO imesh = 1, meshsize
        vel2(imesh) = SQRT( v1(imesh)**2 + v2(imesh)**2 )
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

      meshsize = SIZE(v1)
      DO imesh = 1, meshsize
        vel3(imesh) = SQRT( v1(imesh)**2 + v2(imesh)**2 + v3(imesh)**2 )
      END DO
      
      RETURN
      END SUBROUTINE velocity_module_3D
!----------------------------------------------------------------------
      FUNCTION cm(rgp,rhog,rm,mg,tg)
      ! 
      ! ... computes the inverse of the mixture sound speed

      IMPLICIT NONE
      REAL*8, INTENT(IN), DIMENSION(:) :: rgp, rhog, rm, tg, mg
      REAL*8, DIMENSION(SIZE(tg)) :: cm
      REAL*8, DIMENSION(SIZE(rgp)) :: fact, y

!
! ... Mixture sound speed (Wallis, 1969)
!
      y = rgp / rm
      fact = y + (1.0 - y) * rhog / rl(1)
      cm = SQRT(gammaair * rgas * tg / mg / y ) * fact
!
      RETURN
      END FUNCTION cm
!----------------------------------------------------------------------
      FUNCTION mach(vel, c)
      !
      ! ... compute the Mach number for the mixture

      IMPLICIT NONE
      REAL*8, INTENT(IN), DIMENSION(:) :: vel, c
      REAL*8, DIMENSION(SIZE(vel)) :: mach

      mach = vel / c

      RETURN
      END FUNCTION mach
!----------------------------------------------------------------------
      END MODULE derived_fields
!----------------------------------------------------------------------
