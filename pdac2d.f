!**************************************************************
! PROGRAM: pdac2d (pyroclastic dispersion analysis code, 2d)
! description: parallel multiphase and multicomponent fluid dynamic code
!              for the simulation of pyroclastic dispersion
!              processes
! version: 3.0
! date: January 2001
! authors: A.Neri, G.Macedonio, D.Gidaspow, T. Esposti Ongaro
! parallelized by: T. Esposti Ongaro, C. Cavazzoni, A. Neri
!**************************************************************
!
      PROGRAM pdac2d

      USE atmosphere, ONLY: v0, u0, p0, temp0, uk0, vk0, &
     &                      ep0, epsmx0, gravx, gravz
      USE dimensions
      USE eos_gas, ONLY: bounds_eosg, local_bounds_eosg
      USE gas_constants, ONLY: phij, ckg, mmug, mmugs, mmugek, gmw , &
     &    bounds_gas_constants
      USE gas_solid_density, ONLY: bounds_density, &
     &    local_bounds_density
      USE gas_solid_velocity, ONLY: bounds_velocity, &
     &    local_bounds_velocity
      USE gas_solid_temperature, ONLY: bounds_temperature, &
     &    local_bounds_temperature
      USE gas_solid_viscosity, ONLY: icoh, bounds_viscosity,  &
     &    local_bounds_viscosity
      USE grid, ONLY: dz, dr, itc, ib2, ib1, ib, jb2, jb1, jb
      USE grid, ONLY: no, nso, iob, flic, partition, ghost, &
     &    bounds_blocks, bounds_grid
      USE io_restart, ONLY: taperd, tapebc
      USE iterative_solver, ONLY: inmax, maxout, omega
      USE output_dump, ONLY: nfil, outp
      USE parallel, ONLY: parallel_startup, parallel_hangup, &
     &    mpime, root
      USE particles_constants, ONLY: rl, inrl, kap, &
     &     cmus, phis, cps, dk, nsolid, bounds_part_constants
      USE phases_matrix, ONLY: rlim, bounds_matrix
      USE pressure_epsilon, ONLY: bounds_press_eps, &
     &                            local_bounds_press_eps
      USE reactions, ONLY: irex
      USE roughness_module, ONLY: zrough, deallocate_roughness
      USE initial_conditions, ONLY: setup, epsob, vpob, tpob, ygc0, &
     &    ygcob, upob, vgob, ugob, pob, tgob, epob, lpr, zzero, &
     &    bounds_setup
      USE heat_capacity, ONLY: bounds_hcapgs, local_bounds_hcapgs
      USE time_parameters, ONLY: time, tstop, dt, tpr, tdump, itd, & 
     &                            timestart, rungekut
      USE turbulence, ONLY: iturb, cmut, iss, bounds_turbo, &
     &                      local_bounds_turbo, modturbo
      USE environment, ONLY: cpclock, timing
      USE input_module
!
      IMPLICIT NONE
      CHARACTER(LEN=13) :: errnb
      CHARACTER(LEN=14) :: testnb
      CHARACTER(LEN=1)  :: gnb1, gnb2, gnb3
!
      INTEGER :: n, m, i, j, k
      INTEGER :: kg
      INTEGER :: hun, ten, one
      REAL*8 :: s0, s1, s2, s3, s4
      REAL*8 :: timtot, timprog, timdist, timsetup, timinit
!      INTEGER :: nx, ny, nz, nr
!
      IF(timing) s0 = cpclock()

!
! ... initialize parallel environment
!
      CALL parallel_startup  

! ... Initialize the IBM HW performance monitor
!             call f_hpminit( mpime, 'pdac' )
!
! ... I/O files
!
      hun = (mpime/100)
      ten = (mpime - 100*hun)/10
      one = mpime - 100*hun - 10*ten

      gnb1 = CHAR(hun + 48)      ! the ASCII position of digits
      gnb2 = CHAR(ten + 48)
      gnb3 = CHAR(one + 48)
      errnb = 'pdac2d.err'//gnb1//gnb2//gnb3
      testnb = 'pdac2d.test'//gnb1//gnb2//gnb3
      IF(mpime .EQ. root) THEN
        OPEN(UNIT=5, FILE='pdac2d.dat', STATUS='UNKNOWN')
        OPEN(UNIT=6, FILE='pdac2d.log', STATUS='UNKNOWN')
        OPEN(UNIT=7, FILE='pdac2d.test', STATUS='UNKNOWN')
        OPEN(UNIT=8, FILE='pdac2d.err', STATUS='UNKNOWN')
      ELSE
        OPEN(UNIT=7,FILE=testnb,STATUS='UNKNOWN')
        OPEN(UNIT=8,FILE=errnb,STATUS='UNKNOWN')
      END IF
!
! ... Read Input file
!
      CALL input(5)

! ... set dimensions ...
      ndi=ib2
      ndj=jb2
      timestart = time

! ... allocation of arrays ...(dependence on ndi, ndj)
      CALL bounds_grid
      CALL bounds_press_eps
!
      dr(1:ib2) = delta_r(1:ib2)
      dz(1:jb2) = delta_z(1:jb2)

      no = number_of_block
!
! ... set dimensions ...
      nnso=no
      ncl=nsolid
      nphase=ncl+1

! ... allocation of arrays ...(dependence on nnso, ngas, nsolid, ndi, ndj)
      CALL bounds_blocks
      CALL bounds_density
      CALL bounds_eosg
      CALL bounds_gas_constants
      CALL bounds_hcapgs
      CALL bounds_matrix
      CALL bounds_part_constants
      CALL bounds_setup
      CALL bounds_temperature
      CALL bounds_turbo
      CALL bounds_velocity
      CALL bounds_viscosity

      nso(1:no)   = block_type(1:no)
      iob(:,1:no) = block_bounds(:,1:no)
      ugob(1:no)  = fixed_vgas_r(1:no)
      vgob(1:no)  = fixed_vgas_z(1:no)
      pob(1:no)  = fixed_pressure(1:no)
      epob(1:no)  = fixed_gaseps(1:no)
      tgob(1:no)  = fixed_gastemp(1:no)
      upob(1:nsolid,1:no) = fixed_vpart_r(1:nsolid,1:no)
      vpob(1:nsolid,1:no) = fixed_vpart_z(1:nsolid,1:no)
      epsob(1:nsolid,1:no) = fixed_parteps(1:nsolid,1:no)
      tpob(1:nsolid,1:no) = fixed_parttemp(1:nsolid,1:no)
      ygcob(1:ngas,1:no) = fixed_gasconc(1:ngas,1:no)

      u0 = initial_vgas_r
      v0 = initial_vgas_z
      p0 = initial_pressure
      ep0 = initial_void_fraction
      epsmx0 = max_packing
      temp0 = initial_temperature
      uk0 = initial_vpart_r
      vk0 = initial_vpart_z
      ygc0(1:ngas) = initial_gasconc(1:ngas)
!
      dk(1:nsolid) = diameter(1:nsolid)
      rl(1:nsolid) = density(1:nsolid)
      phis(1:nsolid) = sphericity(1:nsolid)
      cmus(1:nsolid) = viscosity(1:nsolid)
      cps(1:nsolid) = specific_heat(1:nsolid)
      kap(1:nsolid) = thermal_conductivity(1:nsolid)
!
! ... Read restart file
!
      IF(itd.GE.2) THEN 
        CALL taperd
        CALL tapebc
      END IF
!
! ...  Writing log file
! 
      IF (mpime .EQ. root) THEN
        WRITE(6,200) run_name
        WRITE(6,220) itc,ib2,jb2,dr0,dz0,zzero,iturb,modturbo,iss,icoh,irex,ngas
        WRITE(6,611) (dr(i),i=1,ib2)
        WRITE(6,622) (dz(j),j=1,jb2)
        WRITE(6,221) cmut,inmax,maxout,omega
        IF(lpr.LT.2) THEN
         WRITE(6,250) no
          IF(no .NE. 0) THEN
           WRITE(6,252)
            DO  n=1,no
              WRITE(6,254) nso(n),(iob(m,n),m=1,4)
!pe------------------------------
!            IF(nso(n).EQ.5) THEN 
              IF(nso(n).EQ.1.OR.nso(n).eq.5) THEN 
!pe------------------------------
                WRITE(6,253) n,ugob(n),vgob(n),pob(n),epob(n)
                WRITE(6,255) (k,upob(k,n),vpob(k,n), epsob(k,n),tpob(k,n),k=1,nsolid)
                WRITE(6,256) (kg,ygcob(kg,n),kg=1,ngas)
              ENDIF
            END DO
          END IF
         WRITE(6,251) gravx,gravz
         WRITE(6,260) u0,v0,p0,ep0,epsmx0,temp0
         WRITE(6,274)
         WRITE(6,275) (dk(i),rl(i),phis(i),cmus(i),cps(i), kap(i),i=1,nsolid)
         WRITE(6,271) uk0,vk0
         WRITE(6,261) (kg,ygc0(kg),kg=1,ngas)
         WRITE(6,280) itd,time,tstop,dt,tpr,tdump
        END IF
      END IF
!
        DO k=1, nsolid
          inrl(k)=1.D0/rl(k)
        END DO
!
        ib=ib2-2
        ib1=ib2-1
        jb=jb2-2
        jb1=jb2-1
!
! ... Set boundary flags
!
      CALL flic
      IF(timing) s1 = cpclock()
!
! ... Domain decomposition for parallelization 
!
      CALL partition
      IF(timing) s2 = cpclock()
!
! ... Setting ghost cells
!
      CALL ghost
      IF(timing) s3 = cpclock()
!
      CALL local_bounds_velocity
      CALL local_bounds_density
      CALL local_bounds_press_eps
      CALL local_bounds_temperature
      CALL local_bounds_eosg
      CALL local_bounds_viscosity
      CALL local_bounds_hcapgs
      CALL local_bounds_turbo
!
! ... Set initial conditions
!
      CALL setup
!
! ... Distribute inital data among processes
!
      CALL distribute
!
! ... Time advancement loop
!
      CALL prog
!
      IF(timing ) THEN
        s4 = cpclock()
        timtot   = (s4 - s0)/1000.D0
        timprog  = (s4 - s3)/1000.D0
        timdist  = (s3 - s2)/1000.D0
        timsetup = (s2 - s1)/1000.D0
        timinit  = (s1 - s0)/1000.D0
        WRITE(7,900)
        WRITE(7,999) timinit,timsetup,timdist,timprog,timtot
999     FORMAT(6(1X,F9.3))
900     FORMAT('  Init      Part ',  '    Ghost      Prog      Total')
      END IF

! ... terminate the IBM HW performance monitor session
!      call f_hpmterminate( mpime )
      CLOSE(6)
       CLOSE(7)
      CLOSE(8)

      CALL deallocate_roughness( zrough )
! 
! ... Finalize parallel environment
!
      CALL parallel_hangup
!
      STOP
!
! ... log format
!
 100  FORMAT(a30)
 200  FORMAT(1x,'pdac_2d problem identifier -  ',a30)
 220  FORMAT(/,1x,'geometry:',/, &
      ' 1. coordinates (0=rect, 1=cylind): ',i3,/, &
      ' 2. mesh size:    x,r (ib2)=',i3,14x,'y,z (jb2)=',i3,/, &
      ' 3. unIForm cell size:    dr0=',1pe11.4,'  dz0=',1pe11.4,/, &
      '    zero = ',1pe11.4,/, &
      ' 4. turb. model (0=no, 1=sgs-s mod, 2=sgs-l0 mod) =',i2,/, &
      ' 5. model turbulence= ',i3,/, &
      ' 6. solid stress model (0=const visc, 1=sgs visc,',/, &
      '    2=visc+g(e)) =',i2,/, &
      ' 7. cohesion model (0=no cohes, 1=cohes) =',i2,/, &
      ' 8. reaction (0=hydro, 1=thermo, 2=thermo+reaction) =',i2,/, &
      ' 9. no. of gases ',i2)
 221  FORMAT(/,1x,'physical constants:',/, &
      '  cmut = ',1pe11.4,/, &
      '  inmax = ',i3,' maxout = ',i4,' omega = ',1pe11.4)
 250  FORMAT(/,1x,'obstacles no = ',i3 )
 251  FORMAT(/,1x,'gravity (x) = ',f7.2,'  gravity (z) = ',f7.2)
 252  FORMAT(/,1x,'obstacle flag',12x,'obstacle geometry')
 253  FORMAT(/,1x,'obstacle =',i3/ &
      ' uob=',1pe11.4,' vob=',1pe11.4, &
      ' pob=',1pe11.4,' epob=',1pe11.4)
 254  FORMAT(/,1x,12x,i3,6x,4i4)
 255  FORMAT(/,1x,'k-th particle =',i3/ &
      ' upob=',1pe11.4,' vpob=',1pe11.4, &
      ' epsob=',1pe11.4,'  tpob=',1pe11.4)
 256  FORMAT(/,1x,'n-th gas =',i3/ &
      ' ygasob=',1pe11.4)
 260  FORMAT(/,1x,'initial data for gas and solid'/' u0=',1pe11.4, &
      ' v0=',1pe11.4,' p0=',1pe11.4,' ep0=',1pe11.4, &
      '  max. eps=',1pe11.4,'     t0=',1pe11.4)
 261  FORMAT(/,1x,' nth-gas initial value=',i3/ ' ygasc0=',1pe11.4)
 274  FORMAT(/,1x,/,'particulate phase data',/, &
      ' diameter',3x,/, &
      ' microscopic density',3x,/, &
      ' sphericity',3x,/, &
      ' solids viscosity',3x,/, &
      ' solids specific heat',3x,/, &
      ' solids conductivity')
 275  FORMAT(/,2(6x,g10.3))
 271  FORMAT(/,1x,'particle inflow data'/' uk0=',1pe11.4,' vk0=', 1pe11.4)
 280  FORMAT(/,1x,'runtime control:',/, &
      'restart itd = ',i3,/, &
      'initial time = ',1pe11.4,/, &
      'final time = ',1pe11.4,/, &
      'time step = ',1pe11.4,/, &
      'printing = ',1pe11.4,/, &
      'restart dumping = ',1pe11.4)
 611  FORMAT('dr(i) = '/,7(2x,1p7e11.4/))
 622  FORMAT('dz(j) = '/,7(2x,1p7e11.4/))
!
      END PROGRAM pdac2d
