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
                gc_bulk_density(ntot,ngas))

       gas_heat_capacity     = 0.0d0
       gc_mass_fraction      = 0.0d0
       gc_molar_fraction     = 0.0d0
       gc_bulk_density       = 0.0d0

      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE local_bounds_eosg
      USE dimensions
      USE grid, ONLY: ncdom, ncint
      IMPLICIT NONE
!
       ALLOCATE(cg(ncint))
       ALLOCATE(ygc(ngas,ncint),xgc(ngas,ncint),             &
                rgpgc(ncdom,ngas), rgpgcn(ncint,ngas))
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
      SUBROUTINE mole(xg, yg)
!
! ... computes molar fractions
!
      USE dimensions
      USE gas_constants, ONLY: gmw
      IMPLICIT NONE
!
      REAL*8, INTENT(IN) :: yg(:)
      REAL*8, INTENT(OUT) :: xg(:)
      INTEGER :: ig
      REAL*8 :: mol
!
      mol = 0.D0
      DO ig=1,ngas
        mol = mol + yg(ig)/gmw(ig)
      END DO
      DO ig=1,ngas
        xg(ig)=yg(ig)/gmw(ig)/mol
      END DO

      RETURN
      END SUBROUTINE mole
!----------------------------------------------------------------------
      SUBROUTINE mas(yg, xg)
!
! ... computes mass fractions
!
      USE dimensions
      USE gas_constants, ONLY: gmw
      IMPLICIT NONE
!
      REAL*8, INTENT(IN) :: xg(:)
      REAL*8, INTENT(OUT) :: yg(:)
      INTEGER :: ig
      REAL*8 :: mass
!
      mass = 0.D0
      DO ig=1,ngas
        mass = mass + xg(ig)*gmw(ig)
      END DO
      DO ig=1,ngas
        yg(ig)=xg(ig)*gmw(ig)/mass
        IF (mass == 0.D0) CALL error('mas','gas mass = 0',1)
      END DO

      RETURN
      END SUBROUTINE mas
!----------------------------------------------------------------------
      SUBROUTINE eosg(rags, rog, cpgc, cgas, tg, yg, xg, sieg, p, itemp, &
                      irhog, isound, imesh)
!
! ... updates gas density with new pressure
!
      USE dimensions
      USE grid, ONLY: fl
      USE gas_constants, ONLY: gmw, c_erg, rgas, tzero, hzerog, gammaair
      USE heat_capacity, ONLY:  hcapg
      USE time_parameters, ONLY: time
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: imesh
      REAL*8, INTENT(IN) :: yg(:), xg(:), sieg, p
      REAL*8, INTENT(OUT) :: rog, rags
      REAL*8, INTENT(OUT) :: cpgc(:), cgas
      REAL*8, INTENT(INOUT) :: tg
      INTEGER :: itemp, irhog, isound
!
      REAL*8 :: tgnn, mg, hc, ratmin
      REAL*8 :: tg0, sieg0, cgas0
      INTEGER :: ii, nlmax
      INTEGER :: ig
      PARAMETER( nlmax = 2000) 
      PARAMETER( ratmin = 1.D-8) 
!
      IF (fl(imesh) == 1) THEN
        IF(itemp > 0) THEN
!
! ... iterative inversion of the enthalpy-temperature law
! ... (the gas thermal capacity depends on the temperature cg=cg(T) )
          tg0=tg
          sieg0=sieg
          DO ii = 1, nlmax
           tgnn = tg
            CALL hcapg(cpgc(:), tg)
            hc=0.D0
            DO ig=1,ngas
              hc=c_erg*cpgc(ig)*yg(ig)+hc
            END DO
            IF (tg == tg0) cgas0 = hc
            cgas = hc
            tg = tzero+(sieg-hzerog)/cgas
            IF (DABS((tgnn-tg)/tgnn) <= ratmin) GOTO 223
          END DO
          WRITE(8,*) 'max number of iteration reached in eosg'
          WRITE(8,*) 'time:',time, 'cell:',imesh
          WRITE(8,*) 'temperature:',tg0, 'enthalpy:',sieg0
          WRITE(8,*) 'specific heat:',cgas0
          CALL error( ' eosg ', ' max number of iteration reached in eosg ', 1 )
  223     CONTINUE
        ENDIF
!
        IF(irhog > 0) THEN
          mg=0.D0
          DO ig=1,ngas
            mg = mg + xg(ig) * gmw(ig)
          END DO
          rog = p/(rgas*tg)*mg
        ENDIF
!
        IF(isound.GT.0) rags = rog/p/gammaair
      END IF
!
      RETURN
      END SUBROUTINE eosg
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
      USE heat_capacity, ONLY: gc_heat_capacity, solid_heat_capacity, hcapg
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: imesh
      INTEGER :: ig
!
      REAL*8 :: mg, hc
!
      mg = 0.D0
      DO ig = 1, ngas
        mg = mg + gc_molar_fraction(ig,imesh) * gmw(ig)
      END DO

! ... gas density (from equation of state)
!
      gas_density(imesh) = gas_pressure(imesh) / (rgas*gas_temperature(imesh)) * mg
      gas_bulk_density(imesh) = gas_density(imesh) * void_fraction(imesh)

      DO ig = 1, ngas
        gc_bulk_density(imesh,ig) = gc_mass_fraction(ig,imesh) * gas_bulk_density(imesh)
      END DO

!
! compute heat capacity (constant volume) for gas mixture
!
      CALL hcapg(gc_heat_capacity(:,imesh), gas_temperature(imesh))

      hc = 0.D0
      DO ig=1,ngas
        hc = hc + c_erg*gc_heat_capacity(ig,imesh)*gc_mass_fraction(ig,imesh)
      END DO 
      gas_heat_capacity(imesh) = hc

      gas_enthalpy(imesh) = (gas_temperature(imesh)-tzero) * gas_heat_capacity(imesh) + hzerog
!
      RETURN
      END SUBROUTINE cnvertg
!----------------------------------------------------------------------
      END MODULE eos_gas
!----------------------------------------------------------------------
