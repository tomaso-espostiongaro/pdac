!----------------------------------------------------------------------
      MODULE eos_gas
!----------------------------------------------------------------------
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
      USE domain_decomposition, ONLY: ncdom, ncint
      IMPLICIT NONE
!
       ALLOCATE(cg(ncint))
       ALLOCATE(ygc(ncdom,ngas),xgc(ncint,ngas))
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
        WRITE(7,*) 'WARNING (mole) zero or negative mol ', mol
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
! ... Compute squared gas sound speed
!
      USE gas_constants, ONLY: gammaair
      IMPLICIT NONE
!
      REAL*8, INTENT(OUT) :: sqc
      REAL*8, INTENT(IN) :: rog, p
!
      sqc = gammaair * p / rog
!
      
      RETURN
      END SUBROUTINE csound
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
      USE dimensions
      USE domain_decomposition, ONLY: meshinds
      USE grid, ONLY: flag, z, zb
      USE gas_constants, ONLY: gmw, c_joule, rgas, tzero, hzerog, gammaair
      USE gas_constants, ONLY: gas_type
      USE parallel, ONLY: mpime
      USE specific_heat_module, ONLY: hcapg
      USE time_parameters, ONLY: time
      USE immersed_boundaries, ONLY: topo2d_c, topo2d_x, topo2d_y
      IMPLICIT NONE
!
      REAL*8, INTENT(IN) :: sieg
      REAL*8, INTENT(IN) :: yg(:)
      REAL*8, INTENT(OUT) :: cpgc(:), cgas
      REAL*8, INTENT(INOUT) :: tg
      INTEGER, INTENT(IN) :: ijk
      INTEGER, INTENT(INOUT) :: info
!
      REAL*8 :: tgnn, mg, ratmin
      REAL*8 :: tg0, sieg0, cgas0
      INTEGER :: ii, nlmax
      INTEGER :: ig, i, j, k, imesh
      PARAMETER( nlmax = 2000) 
      PARAMETER( ratmin = 1.D-8) 

          tg0   = tg
          sieg0 = sieg

          IF( sieg <= 0.0d0 ) THEN
            CALL meshinds(ijk,imesh,i,j,k)
            WRITE(7,*) 'WARNING (caloric_eosg) zero or negative enthalpy'
            WRITE(7,*) 'coord: ', i,j,k,' sieg= ', sieg
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

          IF (info == 1) THEN
!**********************************************************************
            !  Error report
            CALL meshinds(ijk,imesh,i,j,k)
            WRITE(7,*) 'max number of iteration reached in eosg'
            WRITE(7,*) 'time:', time
            WRITE(7,*) 'proc:', mpime
            WRITE(7,*) 'local cell:', ijk , i, j, k
            WRITE(7,*) 'temperature:', tg0, tg
            WRITE(7,*) 'enthalpy:', sieg0
            WRITE(7,*) 'concentrations:', yg(:)
!**********************************************************************
          END IF

          IF( tg <= 0.0d0 ) THEN
            CALL meshinds(ijk,imesh,i,j,k)
            WRITE(7,*) 'WARNING (caloric_eosg) zero or negative temperature'
            WRITE(7,*) 'coord: ',i,j,k,' tg= ', tg
          END IF

      RETURN
      END SUBROUTINE caloric_eosg
!----------------------------------------------------------------------
      END MODULE eos_gas
!----------------------------------------------------------------------
