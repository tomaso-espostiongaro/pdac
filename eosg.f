!----------------------------------------------------------------------
      MODULE eos_gas
!----------------------------------------------------------------------
      IMPLICIT NONE
      SAVE
!
      REAL*8, DIMENSION(:), ALLOCATABLE :: gas_heat_capacity
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: gc_mass_fraction
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: gc_molar_fraction
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: gc_bulk_density
      REAL*8, DIMENSION(:), ALLOCATABLE :: cg
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: ygc, xgc
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: rgpgc, rgpgcn
      REAL*8 :: rags
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE bounds_eosg
      USE dimensions
      IMPLICIT NONE
!
       ALLOCATE(gas_heat_capacity(ntot))
       ALLOCATE(gc_mass_fraction(ngas,ntot),       &
                gc_molar_fraction(ngas,ntot),      &
                gc_bulk_density(ngas,ntot))

       gas_heat_capacity     = 0.0d0
       gc_mass_fraction      = 0.0d0
       gc_molar_fraction     = 0.0d0
       gc_bulk_density       = 0.0d0

      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE local_bounds_eosg
      USE dimensions
      USE grid, ONLY: nijx_l, nij_l
      IMPLICIT NONE
!
       ALLOCATE(cg(nij_l))
       ALLOCATE(ygc(ngas,nij_l),xgc(ngas,nij_l),             &
                rgpgc(ngas,nijx_l), rgpgcn(ngas,nij_l))
!
      cg     = 0.D0
      ygc    = 0.D0
      xgc    = 0.D0
      rgpgc  = 0.D0
      rgpgcn = 0.D0
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE mole(xgc, ygc)
!
! ... computes molar fractions
!
      USE dimensions
      USE gas_constants, ONLY: gmw
      IMPLICIT NONE
!
      REAL*8, INTENT(IN) :: ygc(:)
      REAL*8, INTENT(OUT) :: xgc(:)
      INTEGER :: kg
      REAL*8 :: mol
!
      mol = 0.D0
      DO kg=1,ngas
        mol = mol + ygc(kg)/gmw(kg)
      END DO
      DO kg=1,ngas
        xgc(kg)=ygc(kg)/gmw(kg)/mol
      END DO
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE eosg(rags, rog, cp, cg, tg, ygc, xgc, sieg, p, itemp, &
                      irhog, isound, imesh)
!
! ... updates gas density with new pressure
!
      USE dimensions
      USE gas_constants, ONLY: gmw, c_erg, rgas, tzero, hzerog, gammaair
      USE heat_capacity, ONLY:  hcapg
      USE time_parameters, ONLY: time
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: imesh
      REAL*8, INTENT(IN) :: ygc(:), xgc(:), sieg, p
      REAL*8, INTENT(OUT) :: rog, rags
      REAL*8, INTENT(INOUT) :: tg
      REAL*8 :: cp(:), cg
      INTEGER :: itemp, irhog, isound
!
      REAL*8 :: tgnn, mg, hc, ratmin
      INTEGER :: ii, nlmax
      INTEGER :: kg
      PARAMETER( nlmax = 2000) 
      PARAMETER( ratmin = 1.D-8) 
!
!
      IF(itemp.GT.0) THEN
!
! ... iterative inversion of the enthalpy-temperature law
! ... (the gas thermal capacity depends on the temperature cg=cg(T) )
        DO ii = 1, nlmax
          tgnn = tg
          CALL hcapg(cp(:), tg)
          hc=0.D0
          DO kg=1,ngas
            hc=c_erg*cp(kg)*ygc(kg)+hc
          END DO
          cg = hc
          tg = tzero+(sieg-hzerog)/cg
          IF (DABS((tgnn-tg)/tgnn).LE.ratmin) GOTO 223
        END DO
        WRITE(8,*) 'max number of iteration reached in eosg'
        WRITE(8,*) 'time, imesh ',time, imesh
        CALL error( ' eosg ', ' max number of iteration reached in eosg ', 1 )
  223   CONTINUE
      ENDIF
!
      IF(irhog.GT.0) THEN
        mg=0.D0
        DO kg=1,ngas
          mg = mg + xgc(kg) * gmw(kg)
        END DO
        rog = p/(rgas*tg)*mg
      ENDIF
!
      IF(isound.GT.0) rags = rog/p/gammaair
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE cnvertg(imesh)
!
! ... computes thermodynamic mean quantities
!
      USE dimensions
      USE gas_constants, ONLY: gmw, c_erg,rgas,tzero,hzerog
      USE gas_solid_density, ONLY: gas_density, gas_bulk_density
      USE gas_solid_temperature, ONLY: gas_enthalpy, gas_temperature
      USE pressure_epsilon, ONLY: gas_pressure, void_fraction
      USE reactions, ONLY: irex
      USE heat_capacity, ONLY: gc_heat_capacity, solid_heat_capacity, hcapg
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: imesh
      INTEGER :: kg
!
      REAL*8 :: mg, hc
!
      mg = 0.D0
      DO kg = 1, ngas
        mg = gc_molar_fraction(kg,imesh) * gmw(kg) + mg
      END DO

! ... gas density (equation of state)
!
      gas_density(imesh) = gas_pressure(imesh) / (rgas*gas_temperature(imesh)) * mg
      gas_bulk_density(imesh) = gas_density(imesh) * void_fraction(imesh)

      DO kg = 1, ngas
        gc_bulk_density(kg,imesh) = gc_mass_fraction(kg,imesh) * gas_bulk_density(imesh)
      END DO

!pdac---------------
! control next statement
!      IF(irex.EQ.0) RETURN
!pdac---------------
!
! compute heat capacity (constant volume) for gas mixture
!
      CALL hcapg(gc_heat_capacity(:,imesh),gas_temperature(imesh))
      hc = 0.D0
      DO kg=1,ngas
        hc = c_erg*gc_heat_capacity(kg,imesh)*gc_mass_fraction(kg,imesh) + hc
      END DO 

      gas_heat_capacity(imesh) = hc
      gas_enthalpy(imesh)=(gas_temperature(imesh)-tzero)*gas_heat_capacity(imesh)+hzerog
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      END MODULE eos_gas
