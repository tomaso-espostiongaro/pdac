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
      SUBROUTINE ent_inout4(umn, ucn, upn, usmn, uscn, uspn, d1, d2, grav, k)
!
! ... This routine computes the free in/outflow conditions in the boundary
! ... cell, i.e. the normal component of the velocity and the scalar fields
! ... for /// E a s t , N o r t h , T o p /// boundaries

      USE atmospheric_conditions, ONLY: p_atm, t_atm
      USE gas_constants, ONLY: gmw, gammaair, gamn, rgas
      USE gas_constants, ONLY: gas_type

      REAL*8, INTENT(IN) :: umn, ucn
      REAL*8, INTENT(INOUT) :: upn
      REAL*8, INTENT(IN) :: usmn(:), uscn(:)
      REAL*8, INTENT(INOUT) :: uspn(:)
      REAL*8, INTENT(IN) :: d1, d2, grav
      INTEGER, INTENT(IN) :: k

      REAL*8 :: prif, pnn2, p1nn
      REAL*8 :: zrif, trif, rhorif, cost, costc, rhorifm1
      REAL*8 :: rmcn, rmcnn, rmmn, rm1nn, rm1knn, rm2n, rm1n, rm0n
      REAL*8 :: u1n, u2n, ucnn, upnn
      REAL*8 :: eps, epc, epcn, ep1nn, epnn
      REAL*8 :: tcn, tmn, t1nn, t1n, t0n, t2n
      REAL*8 :: mg
      REAL*8 :: dc
      REAL*8 :: ind1,ind2,indc
      REAL*8 :: pf1, pf2, pfd

      INTEGER :: ig, is
!
! ... Definitions
!
      u1n = ( umn + ucn ) * 0.5D0
      u2n = ( upn + ucn ) * 0.5D0
      dc = ( d2  + d1 ) * 0.5D0
      ind1 = 1.0D0 / d1
      ind2 = 1.0D0 / d2
      indc = 1.0D0 / dc
!
! ... interpolation of mid-point values

      t0n  = tg( n0 )
      t1n  = tg( n1 )
      t2n  = tg( n2 )
      tmn  = ( t0n + t1n ) * 0.5D0
      tcn  = ( t1n + t2n ) * 0.5D0
      epcn = ( ep( n1 ) + ep( n2 ) ) * 0.5D0

      IF( ucn > 0.D0 ) THEN

! ...  OUTFLOW ...
!
! ...  Extrapolations
       ucnn = u1n
       upn = u2n
       !
       t1nn = tmn
!
! ... calculation of the gas volumetric fraction at time (n+1)dt
! ... from the mass balance equation of solids in cell (ijk)

        epc = 0.D0
        DO is = 1, nsolid
          eps = rlk(n1,is) * inrl(is) - ( dt * inrl(is) * ind1 )  *          &
                ( uscn(is)*rlk(n1,is) - usmn(is)*rlk(n0,is) )
          epc = epc + eps
        END DO
        ep1nn = 1.D0 - epc
!
! ... calculation of the mixture density 
! ... from a mass balance equation for the mixture.
! ... Use gas velocity for mixture.

        rm2n=0.D0
        rm1n=0.D0
        rm0n=0.D0
        DO is=1,nsolid
          rm2n=rm2n+rlk(n2,is)
          rm1n=rm1n+rlk(n1,is)
          rm0n=rm0n+rlk(n0,is)
        END DO
        rm1knn = rm1n
        rm2n = rm2n + rgp(n2)
        rm1n = rm1n + rgp(n1)
        rm0n = rm0n + rgp(n0)
        rmcn=(rm2n+rm1n)*0.5D0
        rmmn=(rm1n+rm0n)*0.5D0
!
        rmcnn = rmcn-dt*indc*(u2n*rmcn-u1n*rmmn)
!
! ... Calculation of fluid pressure from the gas equation of state and 
! ... the mixture density transport equation
!
        mg = 0.D0
        DO ig = 1, ngas
          mg = mg + xgc(n1,ig) * gmw( gas_type( ig ) )
        END DO
        rm1nn = ( rgas * t1nn ) / ( ep1nn * mg )      
        p1nn  = rm1nn * (- rm1knn + rm1n - dt * ind1 * ( rm1n * ucn - rm0n * umn ) )

        IF (epc < 1.0D-8) p1nn = p(n1)
!
! ... Calculation of the advanced-time fluid pressure from Momentum balance
! ... equation of the mixture

        p(n2) = -rmcnn*ucnn + rmcn*ucn - dt*indc*( u2n*ucn*rmcn - u1n*umn*rmmn) 
        p(n2) = p(n2) + grav * dt * rmcn
        p(n2) = p(n2) / ( dt * indc ) + p1nn
!
        ep(n2) = ep(n1)
        tg(n2) = tg(n1)
        DO ig = 1, ngas
           ygc(n2,ig) = ygc(n2,ig) - &
                          ucn * dt * ind2 * ( ygc(n2,ig) - ygc(n1,ig) )
        END DO
!
! ... Correct non-physical pressure
!
        IF ( p(n2) <= 0.0D0 ) p(n2) = p(n1)
!
      ELSE IF ( ucn <  0.D0 ) THEN

! ... INFLOW ...
!
! ... Extrapolation
        upn = ucn
!
        zrif = z(k) 
        prif = p_atm(k)
        trif = t_atm(k)

        rhorifm1 = rgas * trif / ( prif * gmw(6) )   !  rhorif ** (-1)
        cost  = prif ** ( 1.D0 / gammaair ) * rhorifm1
        costc = ( gammaair - 1.D0 ) / ( gammaair * cost )
!
! ... Adiabatic inflow
 
        pf1 = prif**gamn
        pf2 = u2n * u2n * 0.5D0 * costc
        pfd = pf1 - pf2

        IF( pfd >= 0 ) THEN
          p(n2) = pfd ** ( 1.d0 / gamn )
        ELSE
          p(n2) = prif
        END IF
!
        ep(n2) = 1.D0
        tg(n2) = trif
        DO ig = 1, ngas
          ygc(n2,ig) = ygc(n1,ig)
        END DO

! ... Correct non-physical pressure

        IF ( p(n2) < 0.0D0 ) p(n2) = prif

!
      ELSE IF ( ucn == 0.D0 ) THEN
!
        upn      = ucn
        p(n2)    = p(n1)
        ep(n2)   = ep(n1)
        tg(n2)   = tg(n1)
        DO ig=1,ngas
          ygc(n2,ig) = ygc(n1,ig)
        END DO

      ENDIF
!
! ... Set primary variables
!
      DO is = 1, nsolid
        IF( uscn(is) >= 0.D0 ) THEN
          rlk(n2,is) = rlk(n1,is)
          uspn(is)  = uscn(is)
        ELSE
          rlk(n2,is) = 0.0D0
          uspn(is)  = 0.0D0
        ENDIF
      END DO
 
      sieg(n2) = sieg(n1)
      DO is = 1, nsolid
        sies(n2,is) = sies(n1,is)
      END DO
!                
      RETURN
      END SUBROUTINE ent_inout4
!----------------------------------------------------------------------
      SUBROUTINE wsb_inout4(ucn, umn, uscn, usmn, d1, d2, k)
!
! ... This routine computes the free in/outflow conditions in the boundary
! ... cell, i.e. the normal component of the velocity and the scalar fields
! ... for /// W e s t,  S o u t h ///  boundaries

      USE atmospheric_conditions, ONLY: p_atm, t_atm
      USE gas_constants, ONLY: gmw, gammaair, gamn, rgas
      USE gas_constants, ONLY: gas_type

      REAL*8, INTENT(IN) :: ucn
      REAL*8, INTENT(INOUT) :: umn
      REAL*8, INTENT(IN) :: uscn(:)
      REAL*8, INTENT(INOUT) :: usmn(:)
      REAL*8, INTENT(IN) :: d1, d2
      INTEGER, INTENT(IN) :: k

      REAL*8 :: prif, pnn2, p1nn
      REAL*8 :: zrif, trif, rhorif, cost, costc
      REAL*8 :: rmcn, rmcnn, rmpn, rm1nn, rm1knn, rm2n, rm1n, rm0n
      REAL*8 :: u1n, u2n, umnn
      REAL*8 :: eps, epc, epcn, ep1nn
      REAL*8 :: tcn, tpn, t1nn, t1n, t0n, t2n
      REAL*8 :: mg
      REAL*8 :: dc
      REAL*8 :: ind1,ind2,indc
      REAL*8 :: pf1, pf2, pfd

      INTEGER :: ig, is
!
! ... definitions
!
      dc    = ( d1 + d2 ) * 0.5D0
      ind1 = 1.0D0 / d1
      ind2 = 1.0D0 / d2
      indc = 1.0D0 / dc 

      u1n  = ( ucn + umn ) * 0.5D0

      t2n  = tg( n2 )
      t1n  = tg( n1 )
      t0n  = tg( n0 )
      tpn  = ( t1n + t0n ) * 0.5D0
      tcn  = ( t1n + t2n ) * 0.5D0
      epcn = ( ep( n1 ) + ep( n2 ) ) * 0.5D0
!
      u2n = umn - ( ucn - umn ) * ind1 * 0.5 * d2
!
! ... OUTFLOW ...

      IF( ucn < 0.D0 ) THEN
!
! ... Extrapolations
        umnn = u1n
        t1nn = tpn
! 
! ... calculation of the gas volumetric fraction at time (n+1)
! ... from the mass balance equation of solids in cell (ijk)

        epc = 0.D0
        DO is = 1, nsolid
         eps = rlk(n1,is) * inrl(is) -                                 &
             dt*inrl(is)*ind1 * (uscn(is)*rlk(n0,is) - usmn(is)*rlk(n1,is))
         epc = epc + eps
        END DO
        ep1nn = 1.D0 - epc
!
! ... Mixture Density at time (n)dt
        rm2n=0.D0
        rm1n=0.D0
        rm0n=0.D0
        DO is=1,nsolid
          rm2n=rm2n+rlk(n2,is)
          rm1n=rm1n+rlk(n1,is)
          rm0n=rm0n+rlk(n0,is)
        END DO
        rm1knn = rm1n
        rm2n = rm2n + rgp(n2)
        rm1n = rm1n + rgp(n1)
        rm0n = rm0n + rgp(n0)
!
        rmcn=(rm2n+rm1n)*0.5D0
        rmpn=(rm1n+rm0n)*0.5D0
!
        rmcnn = rmcn - dt*indc * (u1n*rmpn-u2n*rmcn)
!
! ... calculation of fluid pressure 
! ... from a mass balance equation for the mixture.
! ... Use gas velocity for mixture.
!
        mg=0.D0
        DO ig=1,ngas
          mg = mg + xgc(n1,ig) * gmw(gas_type(ig))
        END DO
        rm1nn = ep1nn*mg/(rgas*t1nn)
        p1nn = (1.D0/rm1nn) * &
               (-rm1knn+rm1n-dt*ind1*(ucn*rm0n-rm1n*umn))
        IF (epc < 1.0D-8) p1nn = p(n1)
!
! ... Calculation of the advanced-time fluid pressure from Momentum balance
! ... equation of the mixture

        p(n2) = -rmcnn*umnn + rmcn*umn - dt*indc *     &
                 ( u1n*ucn*rmpn - u2n*umn*rmcn) + dt*indc * p1nn
        p(n2)=p(n2)/(dt*indc)
!
! ... Correct non-physical pressure
        IF (p(n2) <= 0.0D0) p(n2) = p(n1)

        ep(n2) = ep(n1)
        tg(n2) = tg(n1)
        DO ig=1,ngas
          ygc(n2,ig) = ygc(n2,ig) - &
                         ucn * dt*ind2 * (ygc(n1,ig)-ygc(n2,ig))
        END DO
!
! ... extrapolation of the temperature and solid fraction to time (n+1)dt

        DO is=1,nsolid
          ! IF(nfllt.EQ.4) ws(n2,is)=ws(n1,is)
          rlk(n2,is)=rlk(n1,is)
          usmn(is)=uscn(is)
        END DO

      ELSE IF( ucn > 0.D0 ) THEN

! ... INFLOW ...
!
        zrif = z(k)
        prif = p_atm(k)
        trif = t_atm(k)
        
        !rhorif = prif*gmw(6)/(rgas*trif)
        !cost = prif/(rhorif**gammaair)
        !costc = (gammaair*cost**(1.D0/gammaair))/(gammaair-1.D0)

        rhorif = rgas * trif/( prif*gmw(6))
        cost = prif**(1.D0/gammaair) * rhorif
        costc = (gammaair-1.D0)/ (gammaair*cost)
!
! ... Adiabatic inflow

        ! p(n2) = ( prif**gamn - ( u2n**2 ) / ( 2.D0 * costc ) ) ** ( 1.D0 / gamn )

        pf1   = prif**gamn
        !pf2   = ( u2n**2 ) / ( 2.D0 * costc )
        pf2   = u2n*u2n * 0.5D0 * costc 
        pfd   = ( pf1 - pf2 ) 
        IF( pfd >= 0 ) THEN
          p(n2) = pfd ** ( 1.d0 / gamn )
        ELSE
          p(n2) = prif
        END IF
!
! ... Correct non-physical pressure

        IF ( p(n2) < 0.0D0 ) p(n2) = prif
!
        ep(n2) = 1.D0
        tg(n2) = trif
        DO ig=1,ngas
          ygc(n2,ig) = ygc(n1,ig)
        END DO

      ENDIF
!
! .... Set primary variables
!
      DO is=1,nsolid
        IF( uscn(is) >= 0.D0 ) THEN
          rlk(n2,is)=rlk(n1,is)
          usmn(is)=uscn(is)
        ELSE
          rlk(n2,is)=0.0D0
          usmn(is) =0.0D0
        ENDIF
      END DO
!
      sieg(n2) = sieg(n1)
      DO is=1,nsolid
        sies(n2,is) = sies(n1,is)
      END DO
!
      RETURN
      END SUBROUTINE wsb_inout4
!-----------------------------------------------------------------------
      SUBROUTINE extrapolate(umn, ucn, upn, usmn, uscn, uspn, d0, d1, d2, k)
! ... This routine computes continuous inoutflow conditions by linear
! ... extrapolation
! ... This procedure is suited for low-Mach number regimes
! ... (incompressible flow)
! ... Works on East, North, Top boundaries

      USE atmospheric_conditions, ONLY: p_atm, t_atm, atm_ygc
      USE gas_constants, ONLY: gmw, gammaair, gamn, rgas
      USE gas_constants, ONLY: gas_type
      USE eos_gas, ONLY: mole, thermal_eosg

      REAL*8, INTENT(IN) :: ucn, umn
      REAL*8, INTENT(INOUT) :: upn
      REAL*8, INTENT(IN) :: uscn(:), usmn(:)
      REAL*8, INTENT(INOUT) :: uspn(:)
      REAL*8, INTENT(IN) :: d0, d1, d2
      INTEGER, INTENT(IN) :: k

      REAL*8 :: rls, mg
      REAL*8 :: delta, gradient, c1, c2
      REAL*8 :: zrif, prif, trif, rhorif
      INTEGER :: is, ig
      REAL*8 :: xgcn2(max_ngas)
!
      c1 = 0.5D0*(d0+d1)
      c2 = 0.5D0*(d1+d2)
!
! ... Extrapolate normal velocity in boundary cells
!
      gradient = (ucn-umn)/c1
      delta = gradient * c2
      upn = ucn + delta
      !
      DO is = 1, nsolid
        gradient = (uscn(is)-usmn(is))/c1
        delta = gradient * c2
        uspn(is) = uscn(is) + delta
      END DO

! ... Extrapolate particle concentrations
!
      rls = 0.D0
      DO is=1, nsolid
        gradient = (rlk(n1,is)-rlk(n0,is))/c1
        delta = gradient * c2
        rlk(n2,is) = rlk(n1,is) + delta
        rls = rls + rlk(n2,is)*inrl(is)
      END DO
      IF (rls >= 1.D0) &
          CALL error('bdry','control transport on boundary',1)

! ... Compute void fraction 
!
      ep(n2) = 1.D0 - rls

! ... Atmospheric conditions
!
      zrif = z(k)
      prif = p_atm(k)
      trif = t_atm(k)
      CALL thermal_eosg(rhorif, trif, prif, xgcn2(:))

! ... Extrapolate temperature 
!
      gradient = (tg(n1)-tg(n0))/c1
      delta = gradient * c2
      tg(n2) = tg(n1) + delta
      !
      DO is=1,nsolid
        gradient = (ts(n1,is)-ts(n0,is))/c1
        delta = gradient * c2
        ts(n2,is) = ts(n1,is) + delta
      END DO

! ... Extrapolate gas components
!
      DO ig=1, ngas
        gradient = (ygc(n1,ig)-ygc(n0,ig))/c1
        delta = gradient * c2
        ygc(n2,ig) = ygc(n1,ig) + delta
      END DO
      CALL mole(xgcn2(:), ygc(n2,:))

! ... Extrapolate gas bulk density
!
      gradient = (rgp(n1)-rgp(n0))/c1
      delta = gradient * c2
      rgp(n2) = rgp(n1) + delta
      rog(n2) = rgp(n2) / ep(n2)
!
! ... Compute new pressure by using equation of state
!
      mg = 0.D0
      DO ig = 1, ngas
        mg = mg + xgcn2(ig) * gmw(gas_type(ig))
      END DO

      p(n2) = rog(n2) * tg(n2) * rgas / mg 

      RETURN
      END SUBROUTINE extrapolate
!-----------------------------------------------------------------------
      SUBROUTINE constant_pressure
      IMPLICIT NONE
      RETURN
      END SUBROUTINE constant_pressure
!----------------------------------------------------------------------
      END MODULE inflow_outflow
!-----------------------------------------------------------------------
