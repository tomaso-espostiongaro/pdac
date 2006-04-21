!----------------------------------------------------------------------
      MODULE derived_fields

      USE dimensions
      USE gas_constants, ONLY: gammaair, rgas, gmw, gas_type
      USE particles_constants, ONLY: rl
!
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      FUNCTION epst(eps)
      !
      ! ... computes the total particle volumetric fraction
      
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: eps(:,:)
      REAL*8, DIMENSION(SIZE(eps,DIM=1)) :: epst

      epst = SUM(eps,DIM=2)

      RETURN
      END FUNCTION epst
!----------------------------------------------------------------------
      FUNCTION leps(eps)
      !
      ! ... computes the log10 of the particle volumetric fraction
      
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: eps(:)
      REAL*8, DIMENSION(SIZE(eps)) :: ceps
      REAL*8, DIMENSION(SIZE(eps)) :: leps

      ceps = clamp(eps,-10)
      leps = log10(ceps)

      RETURN
      END FUNCTION leps
!----------------------------------------------------------------------
      FUNCTION clamp(array,expmin)
      !
      ! ... set minimum value

      IMPLICIT NONE
      REAL*8, INTENT(IN) :: array(:)
      INTEGER, INTENT(IN) :: expmin
      REAL*8, DIMENSION(SIZE(array)) :: clamp

      clamp = MAX(10.E0**expmin, array)

      RETURN
      END FUNCTION clamp
!----------------------------------------------------------------------
      FUNCTION mg(xgc)
      !
      ! ... compute the averaged gas molecular weight

      IMPLICIT NONE
      REAL*8, DIMENSION(:,:) :: xgc
      REAL*8, DIMENSION(SIZE(xgc,DIM=1)) :: mg
      INTEGER :: ig

      mg = 0.D0
      DO ig = 1, ngas
        mg = mg + xgc(:,ig) * gmw(gas_type(ig))
      END DO

      RETURN
      END FUNCTION mg
!----------------------------------------------------------------------
      FUNCTION rhog(p,tg,xgc)
      !
      ! ... computes gas thermodynamic density

      IMPLICIT NONE
      REAL*8, INTENT(IN) :: p(:), tg(:), xgc(:,:) 
      REAL*8, DIMENSION(SIZE(p)) :: rhog
      REAL*8, DIMENSION(SIZE(p)) :: mgas

      mgas = mg(xgc) 
      rhog = p / ( rgas * tg ) * mgas
      
      RETURN
      END FUNCTION rhog
!----------------------------------------------------------------------
      FUNCTION ygc(xgc)
      !
      ! ... computes gas components mass fractions
      
      IMPLICIT NONE
      REAL*8, INTENT(IN), DIMENSION(:,:) :: xgc
      REAL*8, DIMENSION(SIZE(xgc,1),SIZE(xgc,2)) :: ygc
      INTEGER :: ig

      DO ig=1,ngas
        ygc(ig,:) = xgc(:,ig) * gmw(gas_type(ig)) / mg(xgc)
      END DO

      RETURN
      END FUNCTION ygc
!----------------------------------------------------------------------
      FUNCTION ep(eps)
      !
      ! ... computes the void fraction

      IMPLICIT NONE
      REAL*8, INTENT(IN) :: eps(:,:)
      REAL*8, DIMENSION(SIZE(eps,DIM=1)) :: ep

      ep = 1.D0 - SUM(eps,DIM=2)

      RETURN
      END FUNCTION ep
!----------------------------------------------------------------------
      FUNCTION rgp(eps,p,tg,xgc)
      !
      ! ... computes gas bulk density

      IMPLICIT NONE
      REAL*8, INTENT(IN) :: eps(:,:)
      REAL*8, INTENT(IN) :: p(:), tg(:), xgc(:,:)
      REAL*8, DIMENSION(SIZE(p)) :: rgp

      rgp = ep(eps) * rhog(p,tg,xgc)

      END FUNCTION rgp
!----------------------------------------------------------------------
      FUNCTION rlk(eps)
      !
      ! ... computes solid bulk densities

      IMPLICIT NONE
      INTEGER :: ijk
      REAL*8, INTENT(IN) :: eps(:,:)
      REAL*8, DIMENSION(SIZE(eps,1),SIZE(eps,2)) :: rlk

      DO ijk = 1, SIZE(eps,DIM=1)
        rlk(ijk,:) = eps(ijk,:) * rl(:)
      END DO

      END FUNCTION rlk
!----------------------------------------------------------------------
      FUNCTION rhom(eps,p,tg,xgc)
      !
      ! ... computes gas-particle mixture density

      IMPLICIT NONE
      REAL*8, INTENT(IN) :: eps(:,:)
      REAL*8, INTENT(IN) :: p(:), tg(:), xgc(:,:)
      REAL*8, DIMENSION(SIZE(p)) :: rhom
      REAL*8, DIMENSION(SIZE(eps,1),SIZE(eps,2)) :: bds

      bds = rlk(eps)
      rhom = rgp(eps,p,tg,xgc) + SUM(bds,DIM=2) 

      RETURN
      END FUNCTION rhom
!----------------------------------------------------------------------
      FUNCTION rhos(eps)
      !
      ! ... computes particle mixture density

      IMPLICIT NONE
      REAL*8, INTENT(IN) :: eps(:,:)
      REAL*8, DIMENSION(SIZE(eps,1)) :: rhos
      REAL*8, DIMENSION(SIZE(eps,1),SIZE(eps,2)) :: bds

      bds = rlk(eps)
      rhos = SUM(bds,DIM=2) 

      RETURN
      END FUNCTION rhos
!----------------------------------------------------------------------
      FUNCTION velm(ug,us,eps,p,tg,xgc)
      !
      ! ... computes gas-particle mixture velocity (one component)

      IMPLICIT NONE
      REAL*8, INTENT(IN), DIMENSION(:,:) :: eps, us
      REAL*8, INTENT(IN), DIMENSION(:,:) :: xgc
      REAL*8, INTENT(IN), DIMENSION(:) :: ug, p,tg
      REAL*8, DIMENSION(SIZE(ug)) :: velm
      
      velm = rgp(eps,p,tg,xgc) * ug + SUM( rlk(eps) * us ,DIM=2) 
      velm = velm / rhom(eps,p,tg,xgc)

      RETURN
      END FUNCTION velm
!----------------------------------------------------------------------
      FUNCTION pdyn(rm,velom)
      !
      ! ... computes the mixture dynamic pressure
      
      IMPLICIT NONE
      REAL*8, INTENT(IN), DIMENSION(:) :: rm, velom
      REAL*8, DIMENSION(SIZE(rm)) :: pdyn

      pdyn = 0.5 * rm * velom**2

      RETURN
      END FUNCTION pdyn
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
      FUNCTION vel2(v1, v2)
      !
      ! ... compute the modulus of a velocity field
      IMPLICIT NONE
      REAL*8, INTENT(IN), DIMENSION(:) :: v1, v2
      REAL*8, DIMENSION(SIZE(v1)) :: vel2

      vel2 = SQRT( v1**2 + v2**2 )
      
      RETURN
      END FUNCTION vel2
!----------------------------------------------------------------------
      FUNCTION vel3(v1, v2, v3)
      !
      ! ... compute the modulus of a velocity field
      IMPLICIT NONE
      REAL*8, INTENT(IN), DIMENSION(:) :: v1, v2, v3
      REAL*8, DIMENSION(SIZE(v1)) :: vel3

      vel3 = SQRT( v1**2 + v2**2 + v3**2 )
      
      RETURN
      END FUNCTION vel3
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
