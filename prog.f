!----------------------------------------------------------------------
      SUBROUTINE prog

      USE boundary_conditions, ONLY: boundary
      USE dimensions
      USE eos_gas, ONLY: mole, eosg, rags
      USE eos_gas, ONLY: ygc, rgpgc, rgpgcn, xgc, cg
      USE eos_solid, ONLY: eosl
      USE enthalpy_matrix, ONLY: ftem
      USE gas_solid_density, ONLY: rog, rgp, rgpn, rlk, rlkn
      USE gas_solid_temperature, ONLY: sieg, siegn, siek, siekn, tk, tg
      USE grid, ONLY: nij_l, myij, fl_l
      USE io_restart, ONLY: tapewr
      USE iterative_solver, ONLY: iter
      USE output_dump, ONLY: outp
      USE particles_constants, ONLY: cps
      USE pressure_epsilon, ONLY: p, pn, ep
      USE reactions, ONLY: rexion, irex
      USE tilde_energy, ONLY: htilde
      USE tilde_momentum, ONLY: tilde, euvel
      USE time_parameters, ONLY: time, tpr, tdump, tstop, dt, itd
      USE time_parameters, ONLY: rungekut
      USE turbulence, ONLY: iturb, iss
      USE turbulence, ONLY: sgsg, sgss
      USE gas_solid_viscosity, ONLY: viscon, mug, kapg
      USE gas_components, ONLY: ygas
      USE heat_capacity, ONLY: ck, cp
      USE environment, ONLY: cpclock, timing
!
      IMPLICIT NONE
!
      INTEGER :: irest
      INTEGER :: is, imesh
      INTEGER :: ijk
      INTEGER :: kg, rk
      REAL*8 :: tdump1, tpri
      REAL*8 :: s0, s1, s2, s3, s4, s5, s6, s7, s8, s9
      REAL*8 :: timbdry, timout, timrestart, timtilde, timiter, &
                timygas, timtem, timtot, dt0, tsgsg
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

      tsgsg = 0.0D0

      tdump1=time+tdump
      tpri=time+tpr
      irest = 0
!
!---------------------------------------
      time_sweep: DO
!---------------------------------------
!
       IF( timing ) s1 = cpclock()

       irest=irest+1
       IF ( (itd <= 1) .OR. (irest > 1) ) THEN
!
! ... Compute Boundary Conditions
!
        CALL boundary 
!
! ... Compute derived fields from closure equations
! ... (must be dumped into restart file)
!
        DO ijk = 1, nij_l

           imesh = myij(0, 0, ijk)
! 
! ... Compute molar fractions of gas species
!
           CALL mole(xgc(:,ijk), ygc(:,ijk))
!
! ... Compute gas specific heat, density and temperature from gas Equation of State
!
           CALL eosg(rags, rog(ijk), cp(:,ijk), cg(ijk),   &
                     tg(ijk), ygc(:,ijk), xgc(:,ijk),      &
                     sieg(ijk), p(ijk), 1, 1, 0, imesh)
           rgp(ijk)=rog(ijk)*ep(ijk)
! 
           DO kg=1,ngas
             rgpgc(kg,ijk)=ygc(kg,ijk)*rgp(ijk)
           END DO
! 
! ... Compute particle specific heat and temperatures from Equation of State
!
           DO is=1,nsolid
             IF (irex.GE.0) THEN
               CALL eosl(tk(is,ijk),ck(is,ijk),cps(is),siek(is,ijk),1,1)
             ENDIF
           END DO
        END DO
!
       END IF
!
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
! ... Write RESTART file
!
       IF (time+0.1D0*dt.GT.tdump1) THEN
         CALL collect
         CALL tapewr
         tdump1=tdump1+tdump
       END IF

              IF( timing ) s4 = cpclock()
!
!------------------------------------------------------------ 
       IF (time + 0.1D0*dt .GE. tstop)     EXIT time_sweep
!------------------------------------------------------------ 
!
       DO ijk = 1, nij_l
!
! ... Store fields at time n*dt
!
           pn(ijk) = p(ijk)
           rgpn(ijk) = rgp(ijk)
           IF(irex.GE.0) siegn(ijk) = sieg(ijk)
           DO is=1,nsolid
             rlkn(is,ijk) = rlk(is,ijk)
             IF (irex.GE.0) siekn(is,ijk) = siek(is,ijk)
           END DO
           DO kg=1,ngas
             rgpgcn(kg,ijk) = rgpgc(kg,ijk)
           END DO
!
! ... Compute the temperature-dependent gas viscosity and th. conductivity
! 
           IF(irex.GE.0) THEN
             CALL viscon(mug(ijk), kapg(ijk), xgc(:,ijk), tg(ijk))
           ENDIF
       END DO
!
! ... Compute Turbulent viscosity from Smagorinsky LES model
!
         IF (iturb .GE. 1)  CALL sgsg
         IF (iss .GE. 1)    CALL sgss
!
              IF( timing ) s5 = cpclock()
! 
! ... Store gas and particle momentum densities at time n*dt
!
         CALL euvel
!
! ... Start the Runge-Kutta iteration
!
         dt0 = dt
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
              tsgsg = tsgsg + (s5 - s4)
              timtilde = timtilde + (s6 - s5)
              timiter = timiter + (s7 - s6)
              timygas = timygas + (s8 - s7)
              timtem = timtem + (s9 - s8)

!------------------------------------------
      END DO time_sweep
!------------------------------------------

              IF( timing ) s9 = cpclock()

              IF( timing ) THEN
                timtot = (s9 -s0) / 1000.0d0
                timbdry = timbdry / 1000.0d0
                timout = timout / 1000.0d0
                timrestart = timrestart / 1000.0d0
                tsgsg = tsgsg /1000.0D0
                timtilde = timtilde / 1000.0d0
                timiter = timiter / 1000.0d0
                timygas = timygas / 1000.0d0
                timtem = timtem / 1000.0d0
                WRITE(7,900)
                WRITE(7,999) timbdry, timout, timrestart, tsgsg, timtilde, &
                             timiter, timygas, timtem, timtot
900      FORMAT('  Bdry      Out       Restart   Dyn       Tilde     ',       &
                         'Iter      Ygas      Tem       Total')
999      FORMAT(8(1X,F9.3))

              END IF
! 
      RETURN
      END
