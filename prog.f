!----------------------------------------------------------------------
      SUBROUTINE prog

      USE grid, ONLY: fl_l
      USE dimensions
      USE eos_gas, ONLY: mole, eosg, rags
      USE eos_gas, ONLY: ygc, rgpgc, rgpgcn, xgc, cg
      USE eos_solid, ONLY: eosl
      USE enthalpy_matrix, ONLY: ftem
      USE gas_solid_density, ONLY: rog, rgp, rgpn, rlk, rlkn
      USE gas_solid_temperature, ONLY: sieg, siegn, siek, 
     &    siekn, tk, tg
      USE grid, ONLY: nij_l, myij
      USE io_restart, ONLY: tapewr
      USE iterative_solver, ONLY: iter
      USE output_dump, ONLY: outp
      USE particles_constants, ONLY: nsolid, cps
      USE pressure_epsilon, ONLY: p, pn, ep
      USE reactions, ONLY: rexion, irex
      USE tilde_energy, ONLY: htilde
      USE tilde_momentum, ONLY: tilde, euvel
      USE time_parameters, ONLY: time, tpr, tdump, tstop, dt, itd
      USE time_parameters, ONLY: rungekut
      USE turbulence, ONLY: iturb, iss
      USE turbulence, ONLY: sgsgdyn, sgss, mugt, kapgt
      USE gas_solid_viscosity, ONLY: viscon
      USE gas_solid_viscosity, ONLY: mug, kapg
      USE gas_components, ONLY: ygas
      USE th_capacity, ONLY: ck, cp
      USE environment, ONLY: cpclock, timing
!
      IMPLICIT NONE
!
      INTEGER :: irest
      INTEGER :: i,j,k, ij_g
      INTEGER :: ij
      INTEGER :: kg, rk
      REAL*8 :: tdump1, tpri
      REAL*8 :: s0, s1, s2, s3, s4, s5, s6, s7, s8, s9
      REAL*8 :: timbdry, timout, timrestart, timtilde,
     &    timiter, timygas, timtem, timtot, dt0, tsgsgdyn
!
             IF( timing ) s0 = cpclock()
!
      timbdry = 0.0d0
      timout = 0.0d0
      timrestart= 0.0d0
      timtilde= 0.0d0
      timiter= 0.0d0
      timygas= 0.0d0
      timtem= 0.0d0
      timtot= 0.0d0

      tsgsgdyn = 0.0D0

      tdump1=time+tdump
      tpri=time+tpr
      irest = 0
!
!--------START time sweep---------------
      DO
!---------------------------------------
!
              IF( timing ) s1 = cpclock()

       irest=irest+1
       IF (itd.LE.1 .OR. irest .GT.1) THEN
!
! ... Compute Boundary Conditions
!
        CALL bdry 
!
! ... Compute explicit quantities on the mesh
! ... from primitive variables
!
        DO ij = 1, nij_l
         ij_g = myij(0, 0, ij)
         IF(fl_l(ij).EQ.1) THEN
           j = ( ij_g - 1 ) / ndi + 1
           i = MOD( ( ij_g - 1 ), ndi) + 1
! 
! ... Compute mole fractions of gas species
!
           CALL mole(xgc(:,ij), ygc(:,ij))
!
! ... Compute gas density from Equation of State
!
           CALL eosg(rags, rog(ij), cp(:,ij), cg(ij),
     &       tg(ij), ygc(:,ij), xgc(:,ij),
     &       sieg(ij), p(ij), 1, 1, 0, ij_g)
           rgp(ij)=rog(ij)*ep(ij)
! 
! ... Store explicit quantities at time n*dt
!
           rgpn(ij)=rgp(ij)
           pn(ij)=p(ij)
           DO kg=1,ngas
             rgpgc(kg,ij)=ygc(kg,ij)*rgp(ij)
             rgpgcn(kg,ij)=rgpgc(kg,ij)
           END DO
!
! ... Compute the temperature-dependent gas viscosity and conductivity
! 
           IF(irex.GE.0) THEN
             siegn(ij)=sieg(ij)
             CALL viscon(mug(ij), kapg(ij), xgc(:,ij),
     &                   tg(ij))
           ENDIF
! 
! ... Compute Thermodnamic variables for particles
! ... and store explicit quantitites
!
           DO k=1,nsolid
             rlkn(k,ij)=rlk(k,ij)
             IF (irex.GE.0) THEN
               siekn(k,ij) = siek(k,ij)
               CALL eosl(tk(k,ij),ck(k,ij),cps(k),siek(k,ij),1,1)
             ENDIF
           END DO
         END IF
        END DO
       END IF 

              IF( timing ) s2 = cpclock()
!
! ... Write OUTPUT file
!
       IF(time+0.1D0*dt.GE.tpri) THEN
         CALL collect
         CALL outp
         tpri=tpri+tpr
       ENDIF

              IF( timing ) s3 = cpclock()
!
! ... Write restart file
!
       IF ( (time+0.1D0*dt .GT. tdump1) 
!     & .OR. (time+0.1D0*dt .GT. tstop )
     &  ) THEN  
         CALL collect
         CALL tapewr
         tdump1=tdump1+tdump
       END IF

              IF( timing ) s4 = cpclock()
!
!--------Exit time sweep------------------- 
       IF (time + 0.1D0 * dt .GE. tstop) EXIT
!------------------------------------------ 
!
! ... Compute Turbulent viscosity for Smagorinsky LES model
!
         IF (iturb .GE. 1) THEN
           CALL sgsgdyn
         ELSE
           mugt = mug 
           kapgt = kapg
         END IF
         IF (iss .EQ. 1)  CALL sgss
              IF( timing ) s5 = cpclock()
! 
! ... Compute gas and particle momentum density at time n*dt
!
         CALL euvel
         dt0 = dt
!
! ... Start the Runge-Kutta iteration
!
         DO rk = 1, rungekut
           dt = dt0/(rungekut+1-rk)
!
! ... Momentum fluxes, gas-particle drag, particle-particle interaction
!
           CALL tilde

              IF( timing ) s6 = cpclock()
!
! ... Iterative solver for momentum-mass pressure coupling
! ... and explicit solver for interphase coupling
!
           CALL iter

                IF( timing ) s7 = cpclock()
! 
! ... Solve the explicit transport equation of gas species
!
           IF (irex .GE. 2) CALL rexion
           CALL ygas

                IF( timing ) s8 = cpclock()
!
! ... Solve the explicit transport equations for enthalpies
!
           IF (irex .GE. 0) THEN
             CALL htilde
             CALL ftem
           END IF
         END DO
!
! ... End the Runge-Kutta iteration
!
         dt = dt0
         time = time + dt

              IF( timing ) s9 = cpclock()
!
              timbdry = timbdry + (s2 - s1)
              timout = timout + (s3 - s2)
              timrestart = timrestart + (s4 - s3)
              tsgsgdyn = tsgsgdyn + (s5 - s4)
              timtilde = timtilde + (s6 - s5)
              timiter = timiter + (s7 - s6)
              timygas = timygas + (s8 - s7)
              timtem = timtem + (s9 - s8)

!------------------------------------------
      END DO
!--------END time sweep--------------------

              IF( timing ) s9 = cpclock()

              IF( timing ) THEN
                timtot = (s9 -s0) / 1000.0d0
                timbdry = timbdry / 1000.0d0
                timout = timout / 1000.0d0
                timrestart = timrestart / 1000.0d0
                tsgsgdyn = tsgsgdyn /1000.0D0
                timtilde = timtilde / 1000.0d0
                timiter = timiter / 1000.0d0
                timygas = timygas / 1000.0d0
                timtem = timtem / 1000.0d0
                WRITE(7,900)
                WRITE(7,999) timbdry, timout, timrestart, 
     &            tsgsgdyn, timtilde,
     &            timiter, timygas, timtem, timtot
900      FORMAT('  Bdry      Out       Restart   Dyn       Tilde     ',
     &                   'Iter      Ygas      Tem       Total')
999      FORMAT(8(1X,F9.3))

              END IF
! 
      RETURN
      END
