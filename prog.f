!----------------------------------------------------------------------
      SUBROUTINE prog

      USE boundary_conditions, ONLY: boundary
      USE control_flags, ONLY: job_type
      USE dimensions
      USE domain_decomposition, ONLY: ncint, myijk
      USE enthalpy_matrix, ONLY: ftem
      USE environment, ONLY: cpclock, timing, elapsed_seconds
      USE eos_gas, ONLY: mole, caloric_eosg, thermal_eosg
      USE eos_gas, ONLY: ygc, rgpgc, xgc, cg
      USE eos_solid, ONLY: eosl
      USE gas_components, ONLY: ygas
      USE gas_solid_density, ONLY: rog, rgp, rlk
      USE gas_solid_temperature, ONLY: sieg, sies, ts, tg
      USE grid, ONLY: fl_l
      USE specific_heat_module, ONLY: cp, ck
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE io_restart, ONLY: tapewr, max_seconds
      USE iterative_solver, ONLY: iter
      USE output_dump, ONLY: outp
      USE particles_constants, ONLY: cps
      USE pressure_epsilon, ONLY: p, ep
      USE reactions, ONLY: rexion, irex
      USE tilde_energy, ONLY: htilde
      USE tilde_momentum, ONLY: allocate_fluxes, deallocate_fluxes
      USE tilde_momentum, ONLY: tilde, fieldn
      USE time_parameters, ONLY: time, tpr, tdump, tstop, dt, itd
      USE time_parameters, ONLY: rungekut
      USE turbulence_model, ONLY: iturb, iss
      USE turbulence_model, ONLY: sgsg, sgss
!
      IMPLICIT NONE
!
      INTEGER :: irest
      INTEGER :: info
      INTEGER :: is
      INTEGER :: ijk
      INTEGER :: ig, rk
      INTEGER :: myrank
      REAL*8 :: tdump1, tpri
      REAL*8 :: s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10
      REAL*8 :: w0, w1, wmax
      REAL*8 :: t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10
      REAL*8 :: p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10
      REAL*8 :: timbdry, timout, timrestart, timtilde, timiter, &
                timygas, timtem, timtot, dt0, tsgsg
      REAL*8 :: cputimbdry, cputimout, cputimrestart, cputimtilde, cputimiter, &
                cputimygas, cputimtem, cputimtot, cputsgsg
      REAL*8 :: mptimbdry, mptimout, mptimrestart, mptimtilde, mptimiter, &
                mptimygas, mptimtem, mptimtot, mptsgsg   
      
      LOGICAL :: stop_now

!
             IF( timing ) then
                s0 = cpclock()
                call cpu_time(t0)
                call MP_WALLTIME(p0,myrank)
             END IF

      w0 = elapsed_seconds()
!
      timbdry = 0.0d0
      timout = 0.0d0
      timrestart= 0.0d0 
      tsgsg = 0.0D0
      timtilde= 0.0d0
      timiter= 0.0d0
      timygas= 0.0d0
      timtem= 0.0d0
      timtot= 0.0d0
      
      cputimbdry = 0.0d0
      cputimout = 0.0d0
      cputimrestart= 0.0d0
      cputsgsg = 0.0D0
      cputimtilde= 0.0d0
      cputimiter= 0.0d0
      cputimygas= 0.0d0
      cputimtem= 0.0d0
      cputimtot= 0.0d0

      mptimbdry = 0.0d0
      mptimout = 0.0d0
      mptimrestart= 0.0d0 
      mptsgsg = 0.0D0
      mptimtilde= 0.0d0
      mptimiter= 0.0d0
      mptimygas= 0.0d0
      mptimtem= 0.0d0
      mptimtot= 0.0d0
      
      tdump1=time+tdump
      tpri=time+tpr
      irest = 0
!
!---------------------------------------
      time_sweep: DO
!---------------------------------------
!
       IF( timing ) then
           s1 = cpclock()
           call cpu_time(t1)
           call MP_WALLTIME(p1,myrank)
       END IF

       WRITE(6,fmt="(/,'* Starting iteration ',I5,' * ')" ) irest
       WRITE(6,fmt="('  Simulated time = ',F20.14)" ) time
!       WRITE(6,fmt="('  time step = ',F6.4)" ) dt
!       WRITE(6,fmt="('  tpri = ',F20.14)" ) tpri
!       WRITE(6,fmt="('  tdump1 = ',F20.14)" ) tdump1
!       WRITE(6,fmt="('  time+0.1D0*dt  = ',F20.14)" ) time+0.1D0*dt

       irest = irest + 1

       IF ( (itd == 1) .OR. (irest > 1) ) THEN
!
! ... Compute Boundary Conditions
!
         CALL boundary
!
! ... Compute derived fields from closure equations
! ... (all these fields must be dumped into restart file)
!
         info = 0
         DO ijk = 1, ncint
!
! ... Compute molar fractions of gas species
!
           CALL mole( xgc(:,ijk), ygc(:,ijk) )
!
! ... Compute gas density from thermal Equation of State
!
           CALL thermal_eosg(rog(ijk), tg(ijk), p(ijk), xgc(:,ijk) )

! ... Compute gas specific heat and gas temperature 
! ... from caloric Equation of State
!
           CALL caloric_eosg(cp(:,ijk), cg(ijk), tg(ijk), ygc(:,ijk), &
                             sieg(ijk), ijk, info)

           rgp(ijk)=rog(ijk)*ep(ijk)
! 
           DO ig=1,ngas
             rgpgc(ijk,ig)=ygc(ig,ijk)*rgp(ijk)
           END DO
! 
! ... Compute particle specific heat and temperatures from Equation of State
!
           DO is=1,nsolid
             CALL eosl(ts(ijk,is),ck(is,ijk),cps(is),sies(ijk,is))
           END DO

         END DO

         CALL parallel_sum_integer( info, 1 )
         IF( info /= 0 ) THEN
           CALL error( ' prog ',' some cells did not converge in eosg ', info )
         END IF 
!
       END IF

              IF( timing ) then
                 s2 = cpclock()
                 call cpu_time(t2)
                 call MP_WALLTIME(p2,myrank)
              END IF            
! 
! ... Write OUTPUT file
!
       IF(time+0.1D0*dt >= tpri) THEN
         CALL outp
         tpri=tpri+tpr
       ENDIF

              IF( timing ) then
                s3 = cpclock()
                call cpu_time(t3)
                call MP_WALLTIME(p3,myrank)
             END IF
!
! ... Write RESTART file
!
       wmax = elapsed_seconds() 
       call parallel_max_real( wmax, 1 )

       stop_now = ( wmax > max_seconds )

       IF( stop_now ) THEN
         WRITE(6,fmt="('  elapsed_seconds exceed max_second',/, &
                     & '  program stopping')" )
       END IF

       IF ( ( time+0.1D0*dt > tdump1 ) .OR. stop_now ) THEN
         CALL tapewr
         tdump1=tdump1+tdump
       END IF

              IF( timing ) then
                s4 = cpclock()
                call cpu_time(t4)
                call MP_WALLTIME(p4,myrank)
              END IF 

       w1 = elapsed_seconds()
       WRITE(6,fmt="('  walltime = ',F10.2,', ',F10.2)") w1, w1-w0
       w0 = w1
!
!------------------------------------------------------------ 
       IF ( ( time + 0.1D0*dt >= tstop ) .OR. stop_now )     EXIT time_sweep
!------------------------------------------------------------ 
! 
! ... Store all independent fields at time n*dt 
! ... for explicit time integration
!
       CALL fieldn
!
! ... Compute Turbulent viscosity from sub-grid-stress (sgs) model
!
       IF (iturb >= 1)  CALL sgsg
       IF (iss >= 1)    CALL sgss
!
              IF( timing ) then
                 s5 = cpclock()
                 call cpu_time(t5)
                 call MP_WALLTIME(p5,myrank)
              END IF

              timbdry = timbdry + (s2 - s1)
              timout = timout + (s3 - s2)
              timrestart = timrestart + (s4 - s3)
              tsgsg = tsgsg + (s5 - s4)
             
              cputimbdry = cputimbdry + (t2 - t1)
              cputimout = cputimout + (t3 - t2) 
              cputimrestart = cputimrestart + (t4 - t3)
              cputsgsg = cputsgsg + (t5 - t4)
          
              mptimbdry = mptimbdry + (p2 - p1)
              mptimout = mptimout + (p3 - p2)
              mptimrestart = mptimrestart + (p4 - p3)
              mptsgsg = mptsgsg + (p5 - p4)    
!
         CALL allocate_fluxes
!
! ... Start the explicit Runge-Kutta iteration
!
         dt0 = dt
         runge_kutta: DO rk = 1, rungekut
!
                IF( timing ) then 
                    s5 = cpclock()
                    call cpu_time(t5)
                    call MP_WALLTIME(p5,myrank)
                END IF
!
           dt = dt0/(rungekut+1-rk)
!
! ... Momentum fluxes, gas-particle drag, particle-particle interaction
!
           CALL tilde

                IF( timing )then
                    s6 = cpclock()
                    call cpu_time(t6)
                    call MP_WALLTIME(p6,myrank)
                END IF
!
! ... Iterative solver for momentum-mass pressure coupling
! ... and explicit solver for interphase coupling
!
           CALL iter
!
                IF( timing ) then
                    s7 = cpclock()
                    call cpu_time(t7)
                    call MP_WALLTIME(p7,myrank)
                END IF
!
! ... Solve the explicit transport equations for enthalpies
!
           CALL htilde
           CALL ftem
!
                IF( timing ) then
                    s8 = cpclock()
                    call cpu_time(t8)
                    call MP_WALLTIME(p8,myrank)
                END IF
! 
! ... Solve the explicit transport equation of gas species
!
           IF (irex > 2) CALL rexion
           CALL ygas

                IF( timing ) then
                    s9 = cpclock()
                    call cpu_time(t9)
                    call MP_WALLTIME(p9,myrank)
                END IF

                timtilde = timtilde + (s6 - s5)
                timiter = timiter + (s7 - s6)
                timtem = timtem + (s8 - s7)
                timygas = timygas + (s9 - s8)
               
                cputimtilde = cputimtilde + (t6 - t5)
                cputimiter = cputimiter + (t7 - t6)
                cputimtem = cputimtem + (t8 - t7)
                cputimygas = cputimygas + (t9 - t8)
 
                mptimtilde = mptimtilde + (p6 - p5)
                mptimiter = mptimiter + (p7 - p6)
                mptimtem = mptimtem + (p8 - p7)
                mptimygas = mptimygas + (p9 - p8)
               

         END DO runge_kutta
! ... End the Runge-Kutta iteration
         dt = dt0
!
         CALL deallocate_fluxes
!
! ... Advance time
         time = time + dt

        CALL myflush( 6 )
!
!------------------------------------------
      END DO time_sweep
!------------------------------------------

              IF( timing ) then
                  s10 = cpclock()
                  call cpu_time(t10)
                  call MP_WALLTIME(p10,myrank)
              END IF



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

                cputimtot = (t10 -t0) 
                
                mptimtot = (p10 -p0) 
                
                WRITE(7,995)'  WALL TIME computed calling SYSTEM_CLOCK (s)'
                WRITE(7,900)'Bdry','Out','Restart','Dyn','Tilde','Iter','Ygas','Tem','Total'
                WRITE(7,999) timbdry, timout, timrestart, tsgsg, timtilde, &
                             timiter, timygas, timtem, timtot
                WRITE(7,*)' WALL TIME computed calling MP_WALLTIME (s)'
                WRITE(7,900)'Bdry','Out','Restart','Dyn','Tilde','Iter','Ygas','Tem','Total'
                WRITE(7,999) mptimbdry, mptimout, mptimrestart, mptsgsg, mptimtilde, &
                             mptimiter, mptimygas, mptimtem, mptimtot
                WRITE(7,*)' CPU TIME computed calling F95 CPU_TIME (s)'
                WRITE(7,900)'Bdry','Out','Restart','Dyn','Tilde','Iter','Ygas','Tem','Total'
                WRITE(7,999) cputimbdry, cputimout, cputimrestart, cputsgsg, cputimtilde, &
                             cputimiter, cputimygas, cputimtem, cputimtot
995      FORMAT(/,A45)                   
900      FORMAT(9(1X,A10))
999      FORMAT(9(1X,F10.2),/)


              END IF
! 
      RETURN
      END SUBROUTINE prog
!----------------------------------------------------------------------
