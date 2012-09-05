!-----------------------------------------------------------------------
      MODULE inflow_outflow
!-----------------------------------------------------------------------
!
      USE dimensions
      USE gas_solid_density, ONLY: rgp, rlk, rog
      USE gas_solid_temperature, ONLY: sieg, tg, sies, ts
      USE grid, ONLY: z
      USE eos_gas, ONLY: xgc, ygc
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
      SUBROUTINE zero_gradient(us, n2, n1)
      USE eos_gas, ONLY: mole
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: us(:)
      INTEGER, INTENT(IN) :: n1, n2
      INTEGER :: is, ig, sgn
      REAL*8 :: rls
!
      sgn = SIGN(1,n2-n1)
!
! ... Extrapolate particle concentrations
!
      rls = 0.D0
      DO is=1, nsolid
        IF (us(is)*sgn > 0.D0) THEN 
          rlk(n2,is) = rlk(n1,is)
        ELSE
          rlk(n2,is) = 0.D0
        END IF
        rls = rls + rlk(n2,is)*inrl(is)
      END DO
      IF (rls >= 1.D0) &
          CALL error('bdry','control transport on boundary',1)

! ... Compute void fraction 
!
      ep(n2) = 1.D0 - rls

! ... Zero gradient on temperature and enthalpies
!
      tg(n2) = tg(n1)
      sieg(n2) = sieg(n1)
      !
      DO is=1,nsolid
        ts(n2,is) = ts(n1,is)
        sies(n2,is) = sies(n1,is)
      END DO

! ... Extrapolate gas components
!
      DO ig=1, ngas
        ygc(n2,ig) = ygc(n1,ig)
      END DO
      CALL mole(xgc(n2,:), ygc(n2,:))

! ... Extrapolate gas bulk density
!
      rgp(n2) = rgp(n1)
      rog(n2) = rgp(n2) / ep(n2)
!
      RETURN
      END SUBROUTINE zero_gradient 
!-----------------------------------------------------------------------
      END MODULE inflow_outflow
!-----------------------------------------------------------------------
