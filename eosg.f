!----------------------------------------------------------------------
      MODULE eos_gas
!----------------------------------------------------------------------
      IMPLICIT NONE
!
! ... gas specific heat

      REAL*8, DIMENSION(:), ALLOCATABLE :: cg 

! ... gas components mass and molar fractions

      REAL*8, DIMENSION(:,:), ALLOCATABLE :: ygc, xgc

! ... gas component bulk densities

      REAL*8, DIMENSION(:,:), ALLOCATABLE :: rgpgc, rgpgcn
!
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE allocate_eosg
      USE dimensions
      USE domain_decomposition, ONLY: ncdom, ncint
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
      USE gas_constants, ONLY: gmw, gas_type
      IMPLICIT NONE
!
      REAL*8, INTENT(IN) :: yg(:)
      REAL*8, INTENT(OUT) :: xg(:)
      INTEGER :: ig
      REAL*8 :: mol
      REAL*8 :: molinv
!
      mol = 0.D0
      DO ig = 1, ngas
        mol = mol + yg(ig) / gmw(gas_type(ig))
      END DO
      !  Carlo consistency check
      IF( ABS( mol ) <= 1.d-10 ) THEN
        WRITE(6,*) 'WARNING (mole) zero or negative mol ', mol
        molinv = 0.0d0
      ELSE
        molinv = 1.0d0 / mol
      END IF
      DO ig = 1, ngas
        xg(ig) = yg(ig) / gmw(gas_type(ig)) * molinv
      END DO

      RETURN
      END SUBROUTINE mole

!----------------------------------------------------------------------

      SUBROUTINE mas(yg, xg)
!
! ... computes mass fractions
!
      USE dimensions
      USE gas_constants, ONLY: gmw, gas_type
      IMPLICIT NONE
!
      REAL*8, INTENT(IN) :: xg(:)
      REAL*8, INTENT(OUT) :: yg(:)
      INTEGER :: ig
      REAL*8 :: mass
      REAL*8 :: massinv
!
      mass = 0.D0
      DO ig=1,ngas
        mass = mass + xg(ig) * gmw(gas_type(ig))
      END DO
      !  Carlo consistency check
      IF( ABS( mass ) < 1.0d-10 ) THEN
        massinv = 0.0d0
      ELSE
        massinv = 1.0d0 / mass
      END IF
      DO ig=1,ngas
        yg(ig) = xg(ig) * gmw(gas_type(ig)) * massinv
      END DO

      RETURN
      END SUBROUTINE mas

!----------------------------------------------------------------------
      SUBROUTINE csound(sqc, rog, p)
!
! ... Compute gas density as a function of temperature and pressure
! ... Compute squared gas sound speed
!
      USE gas_constants, ONLY: gammaair
      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: sqc
!      REAL*8, INTENT(OUT) :: sqcinv
      REAL*8, INTENT(IN) :: rog, p
!
      sqc = gammaair * p / rog
!
!      sqcinv = rog / (gammaair * p )
      
      RETURN
      END SUBROUTINE csound

!----------------------------------------------------------------------

      SUBROUTINE thermal_eosg(rog, tg, p, xg)
!
! ... Compute gas density as a function of temperature and pressure
! ... Compute sqared gas sound speed
!
      USE dimensions, ONLY: ngas
      USE gas_constants, ONLY: gmw, rgas, gas_type
      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: rog
      REAL*8, INTENT(IN) :: tg, p
      REAL*8, INTENT(IN) :: xg(:)
!
      REAL*8 :: mg
      REAL*8 :: teff
      REAL*8, PARAMETER :: tlim = 1.d-4
      INTEGER :: ig
!
      mg = 0.D0
      DO ig = 1, ngas
        mg = mg + xg(ig) * gmw(gas_type(ig))
      END DO

      !  Carlo consistency check
      teff = rgas * tg
      IF( teff > tlim ) THEN
        rog = p / ( rgas * tg ) * mg
      ELSE
        rog = p / ( tlim ) * mg
      END IF
!
      RETURN
      END SUBROUTINE thermal_eosg
!----------------------------------------------------------------------
      SUBROUTINE caloric_eosg(cpgc, cgas, tg, yg, sieg, ijk, info)
!
! ... Iterative inversion of the enthalpy-temperature law
! ... (the gas thermal capacity depends on the temperature cg=cg(T) )
!
      USE dimensions
      USE grid, ONLY: fl_l
      USE gas_constants, ONLY: gmw, c_joule, rgas, tzero, hzerog, gammaair
      USE gas_constants, ONLY: gas_type
      USE parallel, ONLY: mpime
      USE specific_heat_module, ONLY: hcapg
      USE time_parameters, ONLY: time
      IMPLICIT NONE
!
      REAL*8, INTENT(IN) :: sieg
      REAL*8, INTENT(IN) :: yg(:)
      REAL*8, INTENT(OUT) :: cpgc(:), cgas
      REAL*8, INTENT(INOUT) :: tg
      INTEGER, INTENT(IN) :: ijk
      INTEGER, INTENT(INOUT) :: info
!
      REAL*8 :: tgnn, mg, hc, ratmin
      REAL*8 :: tg0, sieg0, cgas0
      INTEGER :: ii, nlmax
      INTEGER :: ig
      PARAMETER( nlmax = 1000) 
      PARAMETER( ratmin = 1.D-8) 

      IF ( fl_l(ijk) == 1 ) THEN

          tg0   = tg
          sieg0 = sieg

          IF( sieg <= 0.0d0 ) THEN
            WRITE(6,*) 'WARNING (caloric_eosg) zero or negative enthalpy'
            WRITE(6,*) ' sieg = ', ijk, sieg
          END IF

          tgnn = tg
          CALL hcapg( cpgc(:), tg )
          cgas0 = 0.D0
          DO ig = 1, ngas
            cgas0 = cpgc( gas_type( ig ) ) * yg(ig) + cgas0
          END DO
          cgas  = cgas0
          tg = tzero + ( sieg - hzerog ) / cgas
          IF ( DABS( ( tgnn - tg ) / tgnn ) <= ratmin ) GOTO 223

          DO ii = 2, nlmax
            tgnn = tg
            CALL hcapg( cpgc(:), tg )
            hc = 0.D0
            DO ig = 1, ngas
              hc = cpgc( gas_type(ig) ) * yg(ig) + hc
            END DO
            cgas = hc
            tg = tzero + ( sieg - hzerog ) / cgas
            IF ( DABS( ( tgnn - tg ) / tgnn ) <= ratmin ) GOTO 223
          END DO

          !**********************************************************************
          !  Error report
          WRITE(8,*) 'max number of iteration reached in eosg'
          WRITE(8,*) 'time:', time, 'proc:', mpime, 'cell:', ijk 
          WRITE(8,*) 'temperature:',tg0, 'enthalpy:',sieg0
          WRITE(8,*) 'specific heat:',cgas0
          WRITE(6,*) 'max number of iteration reached in eosg'
          WRITE(6,*) 'time:', time, 'proc:', mpime, 'cell:', ijk 
          WRITE(6,*) 'temperature:',tg0, 'enthalpy:',sieg0
          WRITE(6,*) 'specific heat:',cgas0
          info = 1
          CALL error( 'eosg', 'max number of iteration reached in eosg', 1)
  223     CONTINUE

          IF( tg <= 0.0d0 ) THEN
            WRITE(6,*) 'WARNING (caloric_eosg) zero or negative temperature'
            WRITE(6,*) ' tg = ', ijk, tg
          END IF
      END IF

      RETURN
      END SUBROUTINE caloric_eosg
!----------------------------------------------------------------------
      SUBROUTINE cnvertg( ijk )
!
! ... setup thermodynamic mean quantities
!
      USE dimensions
      USE gas_constants, ONLY: gmw, c_joule, rgas, tzero, hzerog
      USE gas_constants, ONLY: gas_type
      USE gas_solid_density, ONLY: rog, rgp
      USE gas_solid_temperature, ONLY: sieg, tg
      USE pressure_epsilon, ONLY: p, ep
      USE specific_heat_module, ONLY: cp, hcapg
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: ijk
      INTEGER :: ig
!
      REAL*8 :: mg, hc
!
      mg = 0.D0
      DO ig = 1, ngas
        mg = mg + xgc( ig, ijk ) * gmw(gas_type(ig))
      END DO
!
! ... gas density (from equation of state)
!
      rog(ijk)   = p(ijk) / ( rgas * tg(ijk) ) * mg

      rgp(ijk)   = rog(ijk) * ep(ijk)

      DO ig = 1, ngas
        rgpgc(ijk,ig) = ygc(ig,ijk) * rgp(ijk)
      END DO
!
! compute heat capacity (constant volume) for gas mixture
!
      CALL hcapg( cp(:,ijk), tg(ijk))

      hc = 0.D0
      DO ig = 1, ngas
        hc = hc + cp(gas_type(ig),ijk) * ygc(ig,ijk)
      END DO 
      cg(ijk) = hc

      sieg(ijk) = (tg(ijk)-tzero) * cg(ijk) + hzerog
!
      RETURN
      END SUBROUTINE cnvertg
!----------------------------------------------------------------------
      END MODULE eos_gas
!----------------------------------------------------------------------
