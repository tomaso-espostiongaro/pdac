!----------------------------------------------------------------------
      MODULE eos_gas
!----------------------------------------------------------------------

      USE io_files, ONLY: testunit

      IMPLICIT NONE
!
! ... gas specific heat

      REAL*8, DIMENSION(:), ALLOCATABLE :: cg 

! ... gas components mass and molar fractions

      REAL*8, DIMENSION(:,:), ALLOCATABLE :: ygc, xgc
!
      LOGICAL :: update_eosg
!
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE allocate_eosg
      USE dimensions
      USE domain_mapping, ONLY: ncdom, ncint
      IMPLICIT NONE
!
       ALLOCATE(cg(ncint))
       ALLOCATE(ygc(ncdom,ngas),xgc(ncdom,ngas))
!
      cg     = 0.D0
      ygc    = 0.D0
      xgc    = 0.D0
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
      xg = 0.D0
      mol = 0.D0
      DO ig = 1, ngas
        mol = mol + yg(ig) / gmw(gas_type(ig))
      END DO
      !  Carlo consistency check
      IF( ABS( mol ) <= 1.d-10 ) THEN
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
      REAL*8 FUNCTION scsound(gamm, rog, p)
!
! ... Compute squared gas sound speed
!
      IMPLICIT NONE
!
      REAL*8, INTENT(IN) :: gamm, rog, p
!
      scsound = gamm * p / rog
!
      RETURN
      END FUNCTION scsound
!----------------------------------------------------------------------
      SUBROUTINE thermal_eosg(rog, tg, p, xg)
!
! ... Compute gas density as a function of temperature and pressure
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
      INTEGER :: ig
!
      mg = 0.D0
      DO ig = 1, ngas
        mg = mg + xg(ig) * gmw(gas_type(ig))
      END DO

      rog = p / ( rgas * tg ) * mg
!
      RETURN
      END SUBROUTINE thermal_eosg
!----------------------------------------------------------------------
      SUBROUTINE caloric_eosg(cpgc, cgas, tg, yg, sieg, ijk, info)
!
! ... Iterative inversion of the enthalpy-temperature law
! ... (the gas thermal capacity depends on the temperature cg=cg(T) )
!
      USE control_flags, ONLY: lpr
      USE dimensions
      USE domain_mapping, ONLY: meshinds
      USE grid, ONLY: z, zb
      USE gas_constants, ONLY: gmw, c_joule, rgas, tzero, hzerog, gammaair
      USE gas_constants, ONLY: gas_type
      USE parallel, ONLY: mpime
      USE pressure_epsilon, ONLY: ep, p
      USE specific_heat_module, ONLY: hcapg
      USE time_parameters, ONLY: time
      IMPLICIT NONE
!
      REAL*8, INTENT(IN) :: sieg
      REAL*8, INTENT(IN) :: yg(:)
      REAL*8, INTENT(OUT) :: cpgc(:), cgas
      REAL*8, INTENT(INOUT) :: tg
      INTEGER, INTENT(IN) :: ijk
      INTEGER, INTENT(OUT) :: info
!
      REAL*8 :: tgnn, mg, ratmin
      REAL*8 :: tg0, sieg0, cgas0
      INTEGER :: ii, nlmax
      INTEGER :: ig, i, j, k, imesh
      PARAMETER( nlmax = 2000) 
      PARAMETER( ratmin = 1.D-8) 

          tg0   = tg
          sieg0 = sieg

          IF( sieg <= 0.0D0 .AND. lpr >= 1) THEN
            CALL meshinds(ijk,imesh,i,j,k)
            WRITE(testunit,*) 'WARNING! from proc: ', mpime
            WRITE(testunit,*) 'Zero or negative enthalpy in caloric_eosg'
            WRITE(testunit,*) 'coord: ', i,j,k,' sieg= ', sieg
          END IF

          info = 1
          DO ii = 1, nlmax
            tgnn = tg
            CALL hcapg( cpgc(:), tg )
            cgas = 0.D0
            DO ig = 1, ngas
              cgas = cpgc( gas_type(ig) ) * yg(ig) + cgas
            END DO
            tg = tzero + ( sieg - hzerog ) / cgas
            IF ( DABS( ( tgnn - tg ) / tgnn ) <= ratmin ) THEN
              info = 0
              EXIT
            END IF
          END DO

          IF (info == 1 .AND. lpr >=1) THEN
!**********************************************************************
            !  Error report
            CALL meshinds(ijk,imesh,i,j,k)
            WRITE(testunit,*) 'max number of iteration reached in eosg'
            WRITE(testunit,*) 'PROC          :', mpime
            WRITE(testunit,*) 'time          :', time
            WRITE(testunit,*) 'local cell    :', ijk , i, j, k
            WRITE(testunit,*) 'temperature   :', tg0, tg
            WRITE(testunit,*) 'enthalpy      :', sieg0
            WRITE(testunit,*) 'concentrations:', yg(:)
            WRITE(testunit,*) 'void fraction :', ep(ijk)
            WRITE(testunit,*) 'pressure      :', p(ijk)
!**********************************************************************
          END IF

          IF( tg <= 0.0d0 ) THEN
            info = 1
            IF( lpr >=1 ) THEN
              CALL meshinds(ijk,imesh,i,j,k)
              WRITE(testunit,*) 'WARNING from proc: ', mpime
              WRITE(testunit,*) 'zero or negative temperature in caloric_eosg'
              WRITE(testunit,*) 'coord: ',i,j,k,' tg= ', tg
            END IF
          END IF

      RETURN
      END SUBROUTINE caloric_eosg
!----------------------------------------------------------------------
      END MODULE eos_gas
!----------------------------------------------------------------------
