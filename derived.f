!----------------------------------------------------------------------
      MODULE derived_fields

      USE dimensions
      USE gas_constants, ONLY: gammaair, rgas, gmw, gas_type
      USE particles_constants, ONLY: rl

      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      FUNCTION epst(eps)
      !
      ! ... computes the total particle volumetric fraction
      
      IMPLICIT NONE
      REAL, INTENT(IN) :: eps(:,:)
      REAL, DIMENSION(SIZE(eps,DIM=2)) :: epst

      epst = SUM(eps,DIM=1)

      RETURN
      END FUNCTION epst
!----------------------------------------------------------------------
      FUNCTION leps(eps)
      !
      ! ... computes the log10 of the particle volumetric fraction
      
      IMPLICIT NONE
      REAL, INTENT(IN) :: eps(:)
      REAL, DIMENSION(SIZE(eps)) :: ceps
      REAL, DIMENSION(SIZE(eps)) :: leps

      ceps = clamp(eps,-10)
      leps = log10(ceps)

      RETURN
      END FUNCTION leps
!----------------------------------------------------------------------
      FUNCTION clamp(array,expmin)
      !
      ! ... set minimum value

      IMPLICIT NONE
      REAL, INTENT(IN) :: array(:)
      INTEGER, INTENT(IN) :: expmin
      REAL, DIMENSION(SIZE(array)) :: clamp

      clamp = MAX(10.E0**expmin, array)

      RETURN
      END FUNCTION clamp
!----------------------------------------------------------------------
      FUNCTION mg(xgc)
      !
      ! ... compute the averaged gas molecular weight

      IMPLICIT NONE
      REAL, DIMENSION(:,:) :: xgc
      REAL, DIMENSION(SIZE(xgc,DIM=1)) :: mg
      INTEGER :: ig

      mg = 0.D0
      DO ig = 1, ngas
        mg = mg + xgc(:,ig) * gmw(gas_type(ig))
      END DO

      WRITE(*,*) gas_type(:)
      WRITE(*,*) gmw(gas_type(:))
      WRITE(*,*) SIZE(xgc), SIZE(mg)
      RETURN

      RETURN
      END FUNCTION mg
!----------------------------------------------------------------------
      FUNCTION rhog(p,tg,xgc)
      !
      ! ... computes gas thermodynamic density

      IMPLICIT NONE
      REAL, INTENT(IN) :: p(:), tg(:), xgc(:,:) 
      REAL, DIMENSION(SIZE(p)) :: rhog
      REAL, DIMENSION(SIZE(p)) :: mgas

      mgas = mg(xgc) 
      rhog = p / ( rgas * tg ) * mgas
      
      RETURN
      END FUNCTION rhog
!----------------------------------------------------------------------
      FUNCTION ygc(xgc)
      !
      ! ... computes gas components mass fractions
      
      IMPLICIT NONE
      REAL, INTENT(IN), DIMENSION(:,:) :: xgc
      REAL, DIMENSION(SIZE(xgc,1),SIZE(xgc,2)) :: ygc
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
      REAL, INTENT(IN) :: eps(:,:)
      REAL, DIMENSION(SIZE(eps,DIM=1)) :: ep

      ep = 1.D0 - SUM(eps,DIM=2)

      RETURN
      END FUNCTION ep
!----------------------------------------------------------------------
      FUNCTION rgp(eps,p,tg,xgc)
      !
      ! ... computes gas bulk density

      IMPLICIT NONE
      REAL, INTENT(IN) :: eps(:,:)
      REAL, INTENT(IN) :: p(:), tg(:), xgc(:,:)
      REAL, DIMENSION(SIZE(p)) :: rgp

      rgp = ep(eps) * rhog(p,tg,xgc)

      END FUNCTION rgp
!----------------------------------------------------------------------
      FUNCTION rlk(eps)
      !
      ! ... computes solid bulk densities

      IMPLICIT NONE
      INTEGER :: ijk
      REAL, INTENT(IN) :: eps(:,:)
      REAL, DIMENSION(SIZE(eps,1),SIZE(eps,2)) :: rlk

      DO ijk = 1, ntot
        rlk(ijk,:) = eps(ijk,:) / rl(:)
      END DO

      END FUNCTION rlk
!----------------------------------------------------------------------
      FUNCTION rhom(eps,p,tg,xgc)
      !
      ! ... computes gas-particle mixture density

      IMPLICIT NONE
      REAL, INTENT(IN) :: eps(:,:)
      REAL, INTENT(IN) :: p(:), tg(:), xgc(:,:)
      REAL, DIMENSION(SIZE(p)) :: rhom

      rhom = rgp(eps,p,tg,xgc) + SUM(rlk(eps),DIM=2) 

      RETURN
      END FUNCTION rhom
!----------------------------------------------------------------------
      FUNCTION velm(ug,us,eps,p,tg,xgc)
      !
      ! ... computes gas-particle mixture velocity (one component)

      IMPLICIT NONE
      REAL, INTENT(IN), DIMENSION(:,:) :: eps, us
      REAL, INTENT(IN), DIMENSION(:,:) :: xgc
      REAL, INTENT(IN), DIMENSION(:) :: ug, p,tg
      REAL, DIMENSION(SIZE(ug)) :: velm
      
      velm = rgp(eps,p,tg,xgc) * ug + SUM( rlk(eps) * us ,DIM=2) 
      velm = velm / rhom(eps,p,tg,xgc)

      RETURN
      END FUNCTION velm
!----------------------------------------------------------------------
      FUNCTION pdyn(rhom,velm)
      !
      ! ... computes the mixture dynamic pressure
      
      IMPLICIT NONE
      REAL, INTENT(IN), DIMENSION(:) :: rhom, velm
      REAL, DIMENSION(SIZE(rhom)) :: pdyn

      pdyn = 0.5 * rhom * velm**2

      RETURN
      END FUNCTION pdyn
!----------------------------------------------------------------------
      FUNCTION cm(rgp,rhog,rhom,mg,tg)
      ! 
      ! ... computes the mixture sound speed

      IMPLICIT NONE
      REAL, INTENT(IN), DIMENSION(:) :: rgp, rhog, rhom, tg, mg
      REAL, DIMENSION(SIZE(tg)) :: cm
      REAL, DIMENSION(SIZE(rgp)) :: fact, y

!
! ... Mixture sound speed (Wallis, 1969)
!
!      y = rgp / rhom
!      fact = y + (1.0 - y) * rhog / rl(1)
!      cm = SQRT(gammaair * rgas * tg / mg / y ) * fact
!
! ... limit for particle density >> gas density
!
      cm = rhog / (rhom - rgp) * (gammaair * rgas * tg / mg)
      cm = cm / rgp * rhog
      
      RETURN
      END FUNCTION cm
!----------------------------------------------------------------------
      FUNCTION vel(vx, vz, vy)
      !
      ! ... compute the modulus of a velocity field
      IMPLICIT NONE
      REAL, INTENT(IN), DIMENSION(:) :: vx, vz
      REAL, INTENT(IN), DIMENSION(:), OPTIONAL :: vy
      REAL, DIMENSION(SIZE(vx)) :: vel

      IF (PRESENT(vy)) THEN
        vel = SQRT( vx**2 + vy**2 + vz**2 )
      ELSE
        vel = SQRT( vx**2 + vz**2 )
      END IF
      
      RETURN
      END FUNCTION vel
!----------------------------------------------------------------------
      FUNCTION mach(vel, c)
      !
      ! ... compute the Mach number for the mixture

      IMPLICIT NONE
      REAL, INTENT(IN), DIMENSION(:) :: vel, c
      REAL, DIMENSION(SIZE(vel)) :: mach

      mach = vel / c

      RETURN
      END FUNCTION mach
!----------------------------------------------------------------------
      END MODULE derived_fields
!----------------------------------------------------------------------
