!----------------------------------------------------------------------
      SUBROUTINE main_prog
!----------------------------------------------------------------------

      USE blunt_body, ONLY: bluntb, ibl
      USE boundary_conditions, ONLY: boundary
      USE check_residuals, ONLY: print_mass_flow_rate
      USE check_residuals, ONLY: print_mass_residuals
      USE control_flags, ONLY: job_type, lpr, imr
      USE control_flags, ONLY: implicit_enthalpy, implicit_fluxes
      USE dimensions
      USE domain_mapping, ONLY: ncint, myijk, ncdom
      USE enthalpy_matrix, ONLY: ftem
      USE environment, ONLY: cpclock, timing, elapsed_seconds
      USE eos_gas, ONLY: thermal_eosg, update_eosg
      USE eos_gas, ONLY: xgc
      USE gas_components, ONLY: ygas
      USE gas_solid_density, ONLY: rog, rgp, rlk
      USE gas_solid_temperature, ONLY: sieg, sies, ts, tg
      USE immersed_boundaries, ONLY: immb
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE io_restart, ONLY: tapewr, max_seconds
      USE iterative_solver, ONLY: iter, nit
      USE output_dump, ONLY: outp, shock_tube_out, print_volumes
      USE output_dump, ONLY: write_radial_profile_2d
      USE parallel, ONLY: mpime, root
      USE particles_constants, ONLY: cps
      USE pressure_epsilon, ONLY: p, ep
      USE reactions, ONLY: rexion, irex
      USE runtime_sampling, ONLY: set_sampling, sample_pressure, isrt
      USE tilde_energy, ONLY: htilde
      USE tilde_momentum, ONLY: allocate_fluxes, deallocate_fluxes
      USE tilde_momentum, ONLY: tilde, fieldn
      USE time_parameters, ONLY: time, tpr, tdump, tstop, dt, itd
      USE time_parameters, ONLY: rungekut, sweep, ndump, nprint
      USE turbulence_model, ONLY: iturb, iss
      USE turbulence_model, ONLY: sgsg, sgss
      USE volcano_topography, ONLY: itp
      USE io_files, ONLY: testunit, logunit
!
      IMPLICIT NONE
!
      INTEGER :: ismp, psmp
      INTEGER :: info
      INTEGER :: is
      INTEGER :: ijk,i
      INTEGER :: ig, rk, n
      INTEGER :: myrank
      REAL*8 :: dt0
      REAL*8 :: w0, w1, wmax
      REAL*8 :: s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13
      REAL*8 :: t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13
      REAL*8 :: p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13
      REAL*8 :: timbdry, timfieldn, timturbo, timtilde, timiter, timygas, &
                timtem, timout, timres, timtot
      REAL*8 :: cptimbdry, cptimfieldn, cptimturbo, cptimtilde, &
                cptimiter, cptimygas, cptimtem, cptimout, cptimres, cptimtot
      REAL*8 :: mptimbdry, mptimfieldn, mptimturbo, mptimtilde, &
                mptimiter, mptimygas, mptimtem, mptimout, mptimres, mptimtot
      REAL*8 :: xgcl(max_ngas)
      
      LOGICAL :: stop_now
!
      IF( timing ) then
         s0 = cpclock()
         !call cpu_time(t0)
         !call MP_WALLTIME(p0,myrank)
      END IF
! 
! ... Print OUTPUT every 'nprint' time-steps
! ... Print RESTART every 'ndump' time-steps
!
      nprint = NINT(tpr/dt)
      ndump  = NINT(tdump/dt)
      sweep  = NINT(time/dt)

      IF( mpime == root ) THEN
        WRITE(logunit,*) 'Print OUTPUT  every ', nprint, ' time steps'
        WRITE(logunit,*) 'Print RESTART every ', ndump, ' time steps'
      END IF
!
      !
      !   set timing variables to 0
      !
      timbdry    = 0.0d0
      timfieldn  = 0.0d0
      timturbo   = 0.0D0
      timtilde   = 0.0d0
      timiter    = 0.0d0
      timtem     = 0.0d0
      timygas    = 0.0d0
      timout     = 0.0d0
      timres     = 0.0d0 
      timtot     = 0.0d0
      
!      cptimbdry    = 0.0d0
!      cptimfieldn  = 0.0d0
!      cptimturbo   = 0.0D0
!      cptimtilde   = 0.0d0
!      cptimiter    = 0.0d0
!      cptimtem     = 0.0d0
!      cptimygas    = 0.0d0
!      cptimout     = 0.0d0
!      cptimres     = 0.0d0 
!      cptimtot     = 0.0d0
!      
!      mptimbdry    = 0.0d0
!      mptimfieldn  = 0.0d0
!      mptimturbo   = 0.0D0
!      mptimtilde   = 0.0d0
!      mptimiter    = 0.0d0
!      mptimtem     = 0.0d0
!      mptimygas    = 0.0d0
!      mptimout     = 0.0d0
!      mptimres     = 0.0d0 
!      mptimtot     = 0.0d0
!
      w0 = elapsed_seconds()
!
      IF (immb >= 1) CALL print_volumes
      IF (isrt >= 1) CALL set_sampling
!
!//////////////////////////////////////////////////////////////////////
!
      time_sweep: DO
!
!//////////////////////////////////////////////////////////////////////
!
        sweep = sweep + 1

        IF (lpr > 0) THEN
          WRITE(testunit,fmt="(/,'* Starting iteration ',I8,' * ')" ) sweep
          WRITE(testunit,fmt="('  Simulated time = ',F20.14)" ) time
        END IF
!
        IF( timing ) then
           s1 = cpclock()
           !call cpu_time(t1)
           !call MP_WALLTIME(p1,myrank)
        END IF
!
! ... Compute Boundary Conditions
!
        CALL boundary

        ! ... write initial conditions
        IF (sweep == 1) THEN
          CALL outp
          CALL print_mass_flow_rate
        END IF
!
! ... 'Check init' stops here
!
        IF (itd == 0) EXIT time_sweep
!
        IF( timing ) then
          s2 = cpclock()
          !call cpu_time(t2)
          !call MP_WALLTIME(p2,myrank)
        END IF            
!
! ... Print the total residuals of the mass conservation equation
! ... and the mass flow rate
!
        IF ( imr >= 1 ) CALL print_mass_residuals(sweep)
!
! ... If needed, update gas density (check algorithm)
! ... Notice that the update of the gas density could
! ... result in an error in the gas mass conservation ...
!
        IF (update_eosg) THEN
          DO ijk=1, ncint
            xgcl(1:ngas) = xgc(ijk,:)
            CALL thermal_eosg(rog(ijk),tg(ijk),p(ijk),xgcl(:))
            rgp(ijk) = rog(ijk) * ep(ijk)
          END DO
        END IF
!
! ... Store all independent fields at time n*dt 
! ... for explicit time integration
!
        CALL fieldn
!
        IF( timing ) then
          s3 = cpclock()
          !call cpu_time(t3)
          !call MP_WALLTIME(p3,myrank)
        END IF

! ... Compute Turbulent viscosity from sub-grid-stress (sgs) model
!
        IF ( iturb >= 1 )  CALL sgsg
        IF ( iss   == 1 )  CALL sgss

        IF( timing ) then
          s4 = cpclock()
          !call cpu_time(t4)
          !call MP_WALLTIME(p4,myrank)
        END IF
!
        timbdry    = timbdry    + (s2 - s1)
        timfieldn  = timfieldn  + (s3 - s2)
        timturbo   = timturbo   + (s4 - s3)
!
!        cptimbdry    = cptimbdry    + (t2 - t1)
!        cptimfieldn  = cptimfieldn  + (t3 - t2)
!        cptimturbo   = cptimturbo   + (t4 - t3)
!
!        mptimbdry    = mptimbdry    + (p2 - p1)
!        mptimfieldn  = mptimfieldn  + (p3 - p2)
!        mptimturbo   = mptimturbo   + (p4 - p3)
!
        CALL allocate_fluxes

! ... Start the explicit Runge-Kutta iteration
!
        dt0 = dt
        runge_kutta: DO rk = 1, rungekut

          dt = dt0 / ( rungekut + 1 - rk )
!
          IF( timing ) then 
            s5 = cpclock()
            !call cpu_time(t5)
            !call MP_WALLTIME(p5,myrank)
          END IF
!
! ... Momentum fluxes, gas-particle drag, particle-particle interaction
!
          IF ( .NOT. implicit_fluxes ) THEN
            CALL tilde
            IF ( implicit_enthalpy ) CALL htilde
          END IF
 
          IF (ibl >= 1) CALL bluntb

          IF( timing )then
            s6 = cpclock()
            !call cpu_time(t6)
            !call MP_WALLTIME(p6,myrank)
          END IF
!
! ... Iterative solver for momentum-mass pressure coupling
! ... and explicit solver for interphase coupling
!
          CALL iter
!
          IF( timing ) then
            s7 = cpclock()
            !call cpu_time(t7)
            !call MP_WALLTIME(p7,myrank)
          END IF

! ... Solve the explicit transport equation of gas species
!
          IF ( irex > 2 ) CALL rexion
          CALL ygas
!
          IF( timing ) then
            s8 = cpclock()
            !call cpu_time(t8)
            !call MP_WALLTIME(p8,myrank)
          END IF

! ... Solve explicitly transport equations for enthalpies
!
          IF (.NOT.implicit_enthalpy) THEN
            CALL htilde
            CALL ftem
          END IF
!
          IF( timing ) then
            s9 = cpclock()
            !call cpu_time(t9)
            !call MP_WALLTIME(p9,myrank)
          END IF

          timtilde = timtilde + (s6 - s5)
          timiter  = timiter  + (s7 - s6)
          timygas  = timygas  + (s8 - s7)
          timtem   = timtem   + (s9 - s8)
               
!          cptimtilde = cptimtilde + (t6 - t5)
!          cptimiter  = cptimiter  + (t7 - t6)
!          cptimygas  = cptimygas  + (t8 - t7)
!          cptimtem   = cptimtem   + (t9 - t8)
! 
!          mptimtilde = mptimtilde + (p6 - p5)
!          mptimiter  = mptimiter  + (p7 - p6)
!          mptimygas  = mptimygas  + (p8 - p7)
!          mptimtem   = mptimtem   + (p9 - p8)

        END DO runge_kutta
!
        CALL deallocate_fluxes
!
! ... End the Runge-Kutta iteration and advance time
!
        dt = dt0
 100    time = time + dt
!
! ... Force the writing on the standard output
!
        IF( mpime == root ) &
          WRITE(logunit,fmt="('Iteration: ',I6,' nit: ', I5)" ) sweep, nit
        CALL myflush( logunit )
        IF (lpr > 1) CALL myflush( testunit )
!
        IF( timing ) then
          s10 = cpclock()
          !call cpu_time(t10)
          !call MP_WALLTIME(p10,myrank)
        END IF
!
! ... Write OUTPUT file
! 
        IF(MOD(sweep,nprint) == 0) THEN
                !
                CALL outp
                !CALL shock_tube_out
                !
                ! ... Print the total residuals of the mass conservation 
                ! ... equation and the mass-flow rate
                !
                !CALL print_mass_residuals(sweep)
                !CALL write_radial_profile_2d
        END IF
        IF (isrt >= 1) CALL sample_pressure
!
        IF( timing ) then
          s11 = cpclock()
          !call cpu_time(t11)
          !call MP_WALLTIME(p11,myrank)
        END IF
!
! ... Write RESTART file
!
        wmax = elapsed_seconds() 
        call parallel_max_real( wmax, 1 )

        stop_now = ( wmax > max_seconds )

        IF( stop_now ) THEN
          IF( mpime == root ) &
            WRITE(logunit,fmt="('  elapsed_seconds exceed max_second',/,'  &
                          &  program stopping')" )
        END IF

        IF((MOD(sweep,ndump) == 0) .OR. stop_now) CALL tapewr
!
        IF( timing ) then
            s12 = cpclock()
            !call cpu_time(t12)
            !call MP_WALLTIME(p12,myrank)
        END IF
!
        timout    = timout    + (s11 - s10)
        timres    = timres    + (s12 - s11)

!        cptimout = cptimout + (t11 - t10)
!        cptimres = cptimres + (t12 - t11)
!
!        mptimout  = mptimout  + (p11 - p10)
!        mptimres  = mptimres  + (p12 - p11)
!
        w1 = elapsed_seconds()
        IF (lpr > 0) WRITE(testunit,fmt="('  walltime = ',F10.2,', ',F10.2)") w1, w1-w0
        w0 = w1
!
!
!//////////////////////////////////////////////////////////////////////
!
        IF ( ( (time+0.1D0*dt ) >= tstop ) .OR. stop_now ) &

      EXIT time_sweep
!
!//////////////////////////////////////////////////////////////////////
!
      END DO time_sweep
!
!//////////////////////////////////////////////////////////////////////
!
      IF( timing ) then
          s13 = cpclock()
          !call cpu_time(t13)
          !call MP_WALLTIME(p13,myrank)
      END IF

      IF( timing ) THEN

        timtot     = ( s13 - s0 ) / 1000.0d0
        timbdry    = timbdry / 1000.0d0
        timfieldn  = timfieldn / 1000.0d0
        timturbo   = timturbo / 1000.0D0
        timtilde   = timtilde / 1000.0d0
        timiter    = timiter / 1000.0d0
        timygas    = timygas / 1000.0d0
        timtem     = timtem / 1000.0d0
        timout     = timout / 1000.0d0
        timres     = timres / 1000.0d0

!        cptimtot  = ( t13 - t0 ) 
!        mptimtot   = ( p13 - p0 ) 
                
        IF (mpime == root) THEN
          WRITE(logunit,*) '  WALL TIME computed calling SYSTEM_CLOCK (s)'
          WRITE(logunit,900) 'Bdry','Dyn','Tilde','Iter','Ygas','Tem','Out','Restart','Total'
          WRITE(logunit,999) timbdry,  timturbo, timtilde, timiter, timygas, timtem, &
                       timout, timres, timtot
        END IF

995     FORMAT(/,A45)                   
900     FORMAT(10(1X,A10))
999     FORMAT(9(1X,F10.2),/)

      END IF
! 
      RETURN
!----------------------------------------------------------------------
      END SUBROUTINE main_prog
!----------------------------------------------------------------------
