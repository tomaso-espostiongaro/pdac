!-----------------------------------------------------------------------
      MODULE inflow_outflow
!-----------------------------------------------------------------------
!
      USE atmospheric_conditions, ONLY: t_atm, p_atm, atm_ygc
      USE dimensions
      USE gas_constants, ONLY: gas_type
      USE gas_solid_density, ONLY: rgp, rlk, rog
      USE gas_solid_temperature, ONLY: sieg, tg, sies, ts
      USE grid, ONLY: z
      USE eos_gas, ONLY: xgc, ygc, thermal_eosg
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: p, ep
      USE time_parameters, ONLY: dt
!
      IMPLICIT NONE
!
      INTEGER :: n0, n1, n2
      PUBLIC
!
      SAVE
!-----------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------
      SUBROUTINE inoutflow(ug, us, n2, n1, k, gnorm)
      USE eos_gas, ONLY: mole
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: us(:), ug
      INTEGER, INTENT(IN) :: n1, n2, k
      REAL*8, INTENT(IN) :: gnorm
      INTEGER :: is, ig, sgn
      REAL*8 :: rls,rhomix
!
      sgn = SIGN(1,n2-n1)
!
! ... GAS
      IF (ug*sgn > 0.D0) THEN 
          ! ... OUTFLOW
          tg(n2) = tg(n1)
          sieg(n2) = sieg(n1)
          DO ig=1, ngas
            ygc(n2,ig) = ygc(n1,ig)
            xgc(n2,ig) = xgc(n1,ig)
          END DO
          rog(n2) = rog(n1)
      ELSE IF (ug*sgn <= 0.D0) THEN
          ! ... INFLOW
          tg(n2) = t_atm(k)
          p(n2) = p_atm(k)
          DO ig=1, ngas
            ygc(n2,ig) = atm_ygc(gas_type(ig))
          END DO
          CALL mole(xgc(n2,:), ygc(n2,:))
          CALL thermal_eosg( rog(n2), tg(n2), p(n2), xgc(n2,:) )
      END IF
!
! ... PARTICLES
      rls = 0.D0
      DO is=1, nsolid
        IF (us(is)*sgn > 0.D0) THEN 
          ! ... OUTFLOW
          rlk(n2,is) = rlk(n1,is)
          ts(n2,is) = ts(n1,is)
          sies(n2,is) = sies(n1,is)
        ELSE IF (us(is)*sgn <= 0.D0) THEN
          ! ... INFLOW
          rlk(n2,is) = 0.D0
          ts(n2,is) = t_atm(k)
        END IF
        rls = rls + rlk(n2,is)*inrl(is)
      END DO
      ep(n2) = 1.D0 - rls
      rgp(n2) = rog(n2) * ep(n2)
      rhomix = rgp(n2) + SUM(rlk(n2,:))
      IF (rls >= 1.D0) &
          CALL error('bdry','control transport on boundary',1)

! ... Pressure gradient with hydrostatic contribution
!
        p(n2)=p(n1) + rhomix * gnorm 
!
      RETURN
      END SUBROUTINE inoutflow 
!-----------------------------------------------------------------------
      END MODULE inflow_outflow
!-----------------------------------------------------------------------
