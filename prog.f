!----------------------------------------------------------------------
      SUBROUTINE prog

      USE boundary_conditions, ONLY: boundary2d, boundary3d
      USE control_flags, ONLY: job_type
      USE dimensions
      USE enthalpy_matrix, ONLY: ftem
      USE environment, ONLY: cpclock, timing
      USE eos_gas, ONLY: mole, eosg, rags
      USE eos_gas, ONLY: ygc, rgpgc, xgc, cg
      USE eos_solid, ONLY: eosl
      USE gas_components, ONLY: ygas
      USE gas_solid_density, ONLY: rog, rgp, rlk
      USE gas_solid_temperature, ONLY: sieg, sies, ts, tg
      USE glocal_arrays, ONLY: collect
      USE grid, ONLY: ncint, myijk, fl_l
      USE specific_heat, ONLY: cp, ck
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE io_restart, ONLY: tapewr
      USE iterative_solver, ONLY: iter
      USE output_dump, ONLY: outp_2d, outp_3d, outp_bin
      USE particles_constants, ONLY: cps
      USE pressure_epsilon, ONLY: p, ep
      USE reactions, ONLY: rexion, irex
      USE tilde_energy, ONLY: htilde
      USE tilde_momentum, ONLY: tilde, fieldn
      USE time_parameters, ONLY: time, tpr, tdump, tstop, dt, itd
      USE time_parameters, ONLY: rungekut
      USE turbulence_model, ONLY: iturb, iss
      USE turbulence_model, ONLY: sgsg, sgss
!
      IMPLICIT NONE
!
      INTEGER :: irest
      INTEGER :: is, imesh
      INTEGER :: ijk
      INTEGER :: ig, rk
      REAL*8 :: tdump1, tpri
      REAL*8 :: s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10
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
       IF ( (itd == 1) .OR. (irest > 1) ) THEN
!
! ... Compute Boundary Conditions
!
         IF( job_type == '2D' ) THEN
           CALL boundary2d 
         ELSE IF( job_type == '3D' ) THEN
           CALL boundary3d
         ELSE
           CALL error(' prog ',' wrong job_type ',1)
         END IF
!
! ... Compute derived fields from closure equations
! ... (all these fields must be dumped into restart file)
!
        DO ijk = 1, ncint

           imesh = myijk( ip0_jp0_kp0_, ijk)
! 
! ... Compute molar fractions of gas species
!
           CALL mole( xgc(:,ijk), ygc(:,ijk) )
!
! ... Compute gas specific heat, density and temperature from gas Equation of State
!
           CALL eosg(rags, rog(ijk), cp(:,ijk), cg(ijk),   &
                     tg(ijk), ygc(:,ijk), xgc(:,ijk),      &
                     sieg(ijk), p(ijk), 1, 1, 0, imesh)
           rgp(ijk)=rog(ijk)*ep(ijk)
! 
           DO ig=1,ngas
             rgpgc(ijk,ig)=ygc(ig,ijk)*rgp(ijk)
           END DO
! 
! ... Compute particle specific heat and temperatures from Equation of State
!
           DO is=1,nsolid
             CALL eosl(ts(ijk,is),ck(is,ijk),cps(is),sies(ijk,is),1,1)
           END DO
        END DO
!
       END IF

              IF( timing ) s2 = cpclock()
!
! ... Write OUTPUT file
!
       IF( job_type == '2D' ) THEN
         IF(time+0.1D0*dt >= tpri) THEN
           CALL collect
           CALL outp_2d
           tpri=tpri+tpr
         ENDIF
       ELSE IF( job_type == '3D' ) THEN
         IF(time+0.1D0*dt >= tpri) THEN
           CALL collect
           CALL outp_bin
!           CALL outp_3d
           tpri=tpri+tpr
         ENDIF
       ELSE
         CALL error(' prog ',' wrong job_type ',1)
       END IF

              IF( timing ) s3 = cpclock()
!
! ... Write RESTART file
!
       IF (time+0.1D0*dt > tdump1) THEN
         CALL collect
         CALL tapewr
         tdump1=tdump1+tdump
       END IF

              IF( timing ) s4 = cpclock()
!
!------------------------------------------------------------ 
       IF (time + 0.1D0*dt >= tstop)     EXIT time_sweep
!------------------------------------------------------------ 
! 
! ... Store all independent fields at time n*dt 
! ... for explicit time integration
!
       CALL fieldn
!
! ... Compute Turbulent viscosity from Smagorinsky sub-grid-stress model
!
       IF (iturb >= 1)  CALL sgsg
       IF (iss >= 1)    CALL sgss
!
              IF( timing ) s5 = cpclock()

              timbdry = timbdry + (s2 - s1)
              timout = timout + (s3 - s2)
              timrestart = timrestart + (s4 - s3)
              tsgsg = tsgsg + (s5 - s4)
!
! ... Start the explicit Runge-Kutta iteration
!
         dt0 = dt
         DO rk = 1, rungekut
!
                IF( timing ) s5 = cpclock()
!
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
!
                IF( timing ) s7 = cpclock()
! 
! ... Solve the explicit transport equation of gas species
!
           IF (irex > 2) CALL rexion
           CALL ygas

                IF( timing ) s8 = cpclock()
!
! ... Solve the explicit transport equations for enthalpies
!
           CALL htilde
           CALL ftem
!
                IF( timing ) s9 = cpclock()

                timtilde = timtilde + (s6 - s5)
                timiter = timiter + (s7 - s6)
                timygas = timygas + (s8 - s7)
                timtem = timtem + (s9 - s8)

         END DO
!
! ... End the Runge-Kutta iteration
!
         dt = dt0
         time = time + dt

!------------------------------------------
      END DO time_sweep
!------------------------------------------

              IF( timing ) s10 = cpclock()

              IF( timing ) THEN
                timtot = (s10 -s0) / 1000.0d0
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
      END SUBROUTINE prog
!----------------------------------------------------------------------
