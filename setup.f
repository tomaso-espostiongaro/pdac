!----------------------------------------------------------------------
      MODULE initial_conditions
!----------------------------------------------------------------------
      USE io_files, ONLY: logunit
      IMPLICIT NONE
!
      REAL*8, DIMENSION(:), ALLOCATABLE   :: ugob, vgob, wgob, epob, &
                                             tgob, pob
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: upob, vpob, wpob, epsob, &
                                             tpob, ygcob
      LOGICAL :: density_specified
!
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE allocate_setup
      USE dimensions
      IMPLICIT NONE
!
! ... Specified flow blocks
!
       IF (no > 0) THEN
         ALLOCATE(ugob(no), vgob(no), wgob(no), pob(no), epob(no), tgob(no))
         ALLOCATE(upob(nsolid,no), vpob(nsolid,no), wpob(nsolid,no),        &
                  epsob(nsolid,no), tpob(nsolid,no))
         ALLOCATE(ygcob(max_ngas,no))
       END IF
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE setpar
! ... Set some computational parameters and the atmosphere
!
      USE atmospheric_conditions, ONLY: control_atmosphere, set_atmosphere
!
! ... Set constants
!
      CALL setc
!
! ... Control that atmospheric stratification is consistent
! ... Set the pressure and temperature profile
!
      CALL control_atmosphere
      CALL set_atmosphere
!
! ... Check gas species
!
      CALL gas_check
      
      RETURN
      END SUBROUTINE setpar
!----------------------------------------------------------------------
      SUBROUTINE setup
! ... Set initial conditions
! ... (2D/3D_Compliant and fully parallel)
!
      USE check_residuals, ONLY: print_mass_flow_rate
      USE compute_mean_fields, ONLY: imrt
      USE control_flags, ONLY: job_type, nfil
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions
      USE domain_mapping, ONLY: ncint, meshinds, myijk
      USE dome_conditions, ONLY: idome, set_domec
      USE gas_solid_viscosity, ONLY: mus
      USE grid, ONLY: iob
      USE io_files, ONLY: ventunit
      USE io_restart, ONLY: taperd, tapewr, outp_recover, outp_remap
      USE mass_sink, ONLY: isink, read_mass_loss
      USE output_dump, ONLY: read_mean_fields
      USE parallel, ONLY: mpime, root
      USE particles_constants, ONLY: cmus
      USE time_parameters, ONLY: itd
      USE turbulence_model, ONLY: turbulence_setup, iturb
      USE vent_conditions, ONLY: set_ventc, ivent, mass_flow_rate
      USE vent_conditions, ONLY: u_gas, v_gas, w_gas
      USE vent_conditions, ONLY: u_solid, v_solid, w_solid

      IMPLICIT NONE
!
      INTEGER :: n
      INTEGER :: ig, is
      LOGICAL :: forced = .FALSE.
      REAL*8 :: mfr, mfract
!
! ... Initialize particle's viscosity
! 
      DO is = 1, nsolid
        mus(:,is) = cmus(is)
      END DO
!
! ... Initialize Smagorinsky length scale
!
      IF (iturb >= 1) THEN
        CALL turbulence_setup
      END IF
!
! ... Set initial atmospheric conditions in all computational 
! ... domain (---> including boundary cells <---)
!
      CALL set_atmc
!
! ... Set initial conditions from prescribed input data
!
      SELECT CASE(itd)
!
      CASE (0,1)
!
        DO n = 1, no
          SELECT CASE (iob(n)%typ)
          CASE (1,5) ! ... Specified Fluid Flow in a block
            CALL specified_flow(n)
          CASE (3) ! ... Obstacles
            CONTINUE
          END SELECT
        END DO 
! 
! ... Set initial conditions in boundary cells 
! ... with specified fluid flow 
!
        IF (idome >= 1) CALL set_domec
        IF (ivent >= 1) CALL set_ventc
        CALL cnvert
!
        CALL print_mass_flow_rate(mfr)
        IF (ivent >= 1 .AND. mass_flow_rate /= 0.D0) THEN
          IF (mpime == root) WRITE(ventunit,*) &
             'Correcting velocity to fit mass flow rate'
          mfract = mass_flow_rate / mfr
          u_gas = u_gas * mfract
          v_gas = v_gas * mfract
          w_gas = w_gas * mfract
          u_solid = u_solid * mfract
          v_solid = v_solid * mfract
          w_solid = w_solid * mfract
          CALL set_ventc
          CALL cnvert
        END IF
!
      CALL print_mass_flow_rate(mfr)
!
! ... Read restart file or recover initial
! ... conditions from an output file
!
      CASE (2)
        IF (isink > 0) CALL read_mass_loss(nfil)
        IF (imrt  > 0 .AND. nfil>0) CALL read_mean_fields(nfil)
        CALL taperd
        CALL cnvert
      CASE (3)
        IF (isink > 0) CALL read_mass_loss(nfil)
        IF (imrt  > 0 .AND. nfil>0) CALL read_mean_fields(nfil)
        CALL outp_recover(nfil)
        CALL cnvert
      CASE (4)
        CALL outp_remap(nfil)
        CALL cnvert
      END SELECT
!
      RETURN
      END SUBROUTINE setup
!-----------------------------------------------------------------------
      SUBROUTINE setup_ghost
!
! ... Initialize the flow field variables in the ghost cells
!
      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions
      USE domain_mapping, ONLY: ncint, meshinds, myijk
      USE domain_mapping, ONLY: data_exchange
      USE eos_gas, ONLY: ygc, xgc
      USE gas_solid_density, ONLY: rlk, rgp, rog
      USE gas_solid_temperature, ONLY: tg, ts, sieg, sies
      USE gas_solid_velocity, ONLY: ug, wg, vg
      USE gas_solid_velocity, ONLY: us, vs, ws
      USE pressure_epsilon, ONLY: ep, p
      USE specific_heat_module, ONLY: cp, ck
      IMPLICIT NONE
!
      CALL data_exchange(p)      
      CALL data_exchange(ep)      
      CALL data_exchange(rlk)      
      CALL data_exchange(rgp)      
      CALL data_exchange(tg)      
      CALL data_exchange(ts)      
      CALL data_exchange(sieg)      
      CALL data_exchange(sies)      
      CALL data_exchange(ug)      
      CALL data_exchange(wg)      
      CALL data_exchange(us)      
      CALL data_exchange(ws)      
      IF (job_type == JOB_TYPE_3D) THEN
        CALL data_exchange(vg)      
        CALL data_exchange(vs)      
      END IF
      CALL data_exchange(ygc)      
      CALL data_exchange(xgc)      
!
      RETURN
      END SUBROUTINE setup_ghost
!-----------------------------------------------------------------------
      SUBROUTINE set_atmc
!
      USE atmospheric_conditions, ONLY: p_atm, t_atm, atm_ygc
      USE atmospheric_conditions, ONLY: wind_x, wind_y, wind_z, void_fraction
      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions
      USE domain_mapping, ONLY: ncint, meshinds, myijk
      USE eos_gas, ONLY: ygc, xgc, mole
      USE gas_constants, ONLY: gas_type
      USE gas_solid_density, ONLY: rlk
      USE gas_solid_temperature, ONLY: tg, ts
      USE gas_solid_velocity, ONLY: ug, wg, vg
      USE gas_solid_velocity, ONLY: us, vs, ws
      USE grid, ONLY: z, flag, fluid, free_io, zero_grad
      USE particles_constants, ONLY: rl
      USE pressure_epsilon, ONLY: ep, p
      IMPLICIT NONE
      REAL*8 :: zrif
      INTEGER :: ijk, imesh, i, j, k, is, ig
!
        DO ijk = 1, ncint
          CALL meshinds(ijk,imesh,i,j,k)
          
          zrif=z(k)

          p( ijk ) = p_atm(k)
          tg( ijk ) = t_atm(k)

          ! ... Set initial gas composition, particle concentrations 
          ! .... and temperature
          !
          ep(ijk) = void_fraction

          DO ig = 1, ngas 
            ygc(ijk,ig) = atm_ygc(gas_type(ig))
            ! ... Compute gas components molar fractions 
            ! ... from mass fractions
            !
            CALL mole( xgc(ijk,:), ygc(ijk,:) )
          END DO

          DO is = 1, nsolid
            rlk(ijk,is) = rl(is)*(1.D0-ep(ijk)) / nsolid
            ts(ijk,is)  = tg(ijk)
          END DO

          ! ... Set initial velocity profiles
          !
          SELECT CASE ( flag(ijk) )

          CASE (fluid, free_io, zero_grad)

            ug(ijk) = wind_x
            us(ijk,:) = ug(ijk)
            IF ( job_type == JOB_TYPE_3D) THEN
              vg(ijk) = wind_y
              vs(ijk,:) = vg(ijk)
            END IF
            wg(ijk) = wind_z
            ws(ijk,:) = wg(ijk)

          CASE DEFAULT

            CONTINUE

          END SELECT
!
       END DO
!
      RETURN
      END SUBROUTINE set_atmc
!----------------------------------------------------------------------
      SUBROUTINE cnvert
! ... Set the initial conditions in the ghost cells.
! ... Compute initial conditions depending on restart mode.
! ... Compute derived thermodynamic quantities.
!
      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions
      USE domain_mapping, ONLY: ncint, ncdom
      USE eos_gas, ONLY: mas, xgc, ygc
      USE eos_gas, ONLY: cg, thermal_eosg, mole
      USE gas_constants, ONLY: gas_type, gmw, rgas, tzero, hzerog, hzeros
      USE gas_solid_density, ONLY: rlk, rog, rgp
      USE gas_solid_temperature, ONLY: tg, ts, sieg, sies
      USE particles_constants, ONLY: inrl, cps
      USE pressure_epsilon, ONLY: ep, p
      USE time_parameters, ONLY: itd
      USE specific_heat_module, ONLY: ck, cp, hcapg, hcaps

      IMPLICIT NONE
!
      INTEGER :: ijk, info
      INTEGER :: ig, is, dfg
      REAL*8  :: xgc_def, mass, tem, rls
      REAL*8  :: hc, mg , xgcl(1:max_ngas), eps
!
! ... Compute derived thermodynamic quantities from initial conditions
! ... in the whole computational domain (--> including boundaries <--)
! ... (2D/3D_Compliant)
!
      IF (itd <= 1) THEN
!
        DO  ijk = 1, ncint
!
          ! ... compute gas components molar fractions 
          ! ... from mass fractions
          !
          CALL mole( xgc(ijk,:), ygc(ijk,:) )
          
          ! ... If density is specified tg(ijk) is a density
          ! ... that must be converted into a temperature
          !
          IF (density_specified) THEN
            mass = 0.D0
            DO ig=1,ngas
              mass = mass + xgc(ijk,ig)*gmw(gas_type(ig))
            END DO
            tem  = p(ijk) * mass / rgas / tg(ijk)
            tg(ijk) = tem
            ts(ijk,:) = tem
          END IF

          ! ... compute gas density from thermal equation of state
          ! ... compute gas and gas components bulk densities
          !
          xgcl(1:ngas) = xgc(ijk,:)
          CALL thermal_eosg( rog(ijk), tg(ijk), p(ijk), xgcl(:) )
          rgp(ijk) = rog(ijk) * ep(ijk)

          ! ... compute specific heat from temperature
          !
          CALL hcapg(cp(:,ijk), tg(ijk))
          hc = 0.D0
          DO ig = 1, ngas
            hc = hc + cp(gas_type(ig),ijk) * ygc(ijk,ig)
          END DO
          cg(ijk) = hc

          ! ... compute gas Enthalpy from temperature
          !
          sieg(ijk) = (tg(ijk)-tzero) * cg(ijk) + hzerog

          ! ... compute particle specific heat and entahlpy 
          ! ... from temperature
          !
          DO is = 1, nsolid
            CALL hcaps(ck(is,ijk), cps(is), ts(ijk,is))
            sies(ijk,is) = ( ts(ijk,is) - tzero ) * ck(is,ijk) + hzeros
          END DO

        END DO
!
      ELSE IF (itd == 2) THEN 
!
! ... Binary Restart
!
        DO ijk = 1, ncint
          ! ... compute gas components molar fractions 
          ! ... from mass fractions
          !
          CALL mole( xgc(ijk,:), ygc(ijk,:) )

          ! ... compute void fraction
          !
          rls = 0.0d0
          DO is = 1, nsolid
            rls = rls + rlk(ijk,is) * inrl(is)
          END DO
          ep(ijk) = 1.D0 - rls

          ! ... compute specific heats for gas ...
          !
          CALL hcapg(cp(:,ijk), tg(ijk))
          hc = 0.D0
          DO ig = 1, ngas
            hc = hc + cp(gas_type(ig),ijk) * ygc(ijk,ig)
          END DO
          cg(ijk) = hc

          ! ... and particles
          !
          DO is = 1, nsolid
            CALL hcaps(ck(is,ijk), cps(is), ts(ijk,is))
          END DO
!
        END DO
      ELSE IF (itd > 2) THEN 
!
! ... Output Restart
!
        DO ijk = 1, ncint
          ! ... compute gas components mass fractions 
          ! ... from molar fractions
          !
          CALL mas( ygc(ijk,:), xgc(ijk,:) )

          ! ... compute gas density from thermal equation of state
          !
          xgcl(1:ngas) = xgc(ijk,:)
          CALL thermal_eosg( rog(ijk), tg(ijk), p(ijk), xgcl(:) )

          ! ... compute void fraction
          !
          ep(ijk)  = 1.D0 
          DO is = 1, nsolid
            ep(ijk) = ep(ijk) - rlk(ijk,is)*inrl(is)
          END DO

          ! ... compute gas bulk density
          !
          rgp(ijk) = rog(ijk) * ep(ijk)

          ! ... compute specific heat from temperature
          !
          CALL hcapg(cp(:,ijk), tg(ijk))
          hc = 0.D0
          DO ig = 1, ngas
            hc = hc + cp(gas_type(ig),ijk) * ygc(ijk,ig)
          END DO
          cg(ijk) = hc

          ! ... compute gas Enthalpy from temperature
          !
          sieg(ijk) = (tg(ijk)-tzero) * cg(ijk) + hzerog

          ! ... compute particle specific heat and entahlpy 
          ! ... from temperature
          !
          DO is = 1, nsolid
            CALL hcaps(ck(is,ijk), cps(is), ts(ijk,is))
            sies(ijk,is) = ( ts(ijk,is) - tzero ) * ck(is,ijk) + hzeros
          END DO
!
        END DO
      END IF
!
! ... Set the initial conditions in the ghost cells
!
      CALL setup_ghost
!
      RETURN
      END SUBROUTINE cnvert
!----------------------------------------------------------------------
      SUBROUTINE specified_flow(n)

      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions
      USE domain_mapping, ONLY: ncint, meshinds, myijk
      USE environment, ONLY: cpclock
      USE eos_gas, ONLY: ygc
      USE gas_constants, ONLY: gmw, rgas, gammaair, gas_type
      USE gas_solid_density, ONLY: rlk
      USE gas_solid_temperature, ONLY: tg, ts
      USE gas_solid_velocity, ONLY: ug, wg, vg
      USE gas_solid_velocity, ONLY: us, vs, ws
      USE grid, ONLY: iob, x, y, z
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: ep, p
      USE indijk_module, ONLY: ip0_jp0_kp0_
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n
      REAL*8 :: ran0
      EXTERNAL :: ran0
      REAL*8 :: ygcsum
      REAL*8 :: dist, radius
      INTEGER :: ijk,i,j,k,imesh
      INTEGER :: ig, is, dfg
      INTEGER :: ic, jc, kc
      REAL*8 :: seed, rseed
!
! ... These parameters are used to set the initial conditions
! ... within a specified radius
! ... uncomment to set radius
!        dist = 0.D0
!        radius = 4.0D-1
!        ic = iob(n)%xlo
!        jc = iob(n)%ylo
!        kc = iob(n)%zlo
!
        rseed = cpclock()
        seed = INT(rseed)
!
        DO ijk = 1, ncint
          CALL meshinds(ijk,imesh,i,j,k)
          IF ( k >= iob(n)%zlo .AND. k <= iob(n)%zhi  ) THEN
            IF ( ( j >= iob(n)%ylo .AND. j <= iob(n)%yhi )    &
                                    .OR. job_type == JOB_TYPE_2D ) THEN
              IF ( i >= iob(n)%xlo .AND. i <= iob(n)%xhi  ) THEN
                !
! ... uncomment to set radius (alternative to the triple IF above)
!              dist = DSQRT((x(i)-x(ic))**2+(y(j)-y(jc))**2+(z(k)-z(kc))**2)
!              IF (dist <= radius) THEN
                ug(ijk) = ugob(n)
!
! ... Add a random noise
                !ug(ijk) = ug(ijK) + 1.D0 * (ran0(seed) - 0.5D0)
!
                IF (job_type == JOB_TYPE_3D) vg(ijk) = vgob(n)
!
! ... Add a random noise
                !IF (job_type == JOB_TYPE_3D) &
                  !vg(ijk) = vg(ijK) + 1.D0 * (ran0(seed) - 0.5D0)
!
                wg(ijk) = wgob(n)
                tg(ijk) = tgob(n)
                p(ijk)  = pob(n)
                ep(ijk) = 1.D0 - SUM(epsob(:,n))
                DO ig = 1, ngas
                  ygc(ijk,ig) = ygcob(gas_type(ig),n)
                END DO
                DO is = 1,nsolid
                  ts(ijk,is)  = tpob(is,n)
                  us(ijk,is)  = upob(is,n)
!
! ... Add a random noise
                  !us(ijk,is) = us(ijK,is) + 1.D0 * (ran0(seed) - 0.5D0)
!
                  IF (job_type == JOB_TYPE_3D) vs(ijk,is)  = vpob(is,n)
!
! ... Add a random noise
                  !IF (job_type == JOB_TYPE_3D) &
                    !vs(ijk,is) = vs(ijK,is) + 1.D0 * (ran0(seed) - 0.5D0)
!
                  ws(ijk,is)  = wpob(is,n)
                  rlk(ijk,is) = epsob(is,n)*rl(is)
                END DO
                !
                ! ... check gas components closure relation
                ygcsum = SUM(ygc(ijk,:))
                IF ( ygcsum /= 1.D0 ) THEN
                  ygc(ijk,ngas) = 1.D0 - SUM( ygc(ijk,1:ngas-1) )
              END IF
                !
!               END IF  ! uncomment to set radius
              END IF
            END IF
          END IF

        END DO 

      END SUBROUTINE specified_flow
!----------------------------------------------------------------------
      SUBROUTINE setc
!
      USE dimensions
      USE gas_constants, ONLY: gmw, mmugs, mmugek,     &
     &    gammaair, gamn, c_joule, c_erg, tzero, hzerog, hzeros, rgas
      USE particles_constants, ONLY: particles_constants_set
      USE reactions, ONLY: h1, h2, h3, h4, h5
      IMPLICIT NONE

      INTEGER :: is 
!
! ... molecular weight of chemical components (O2,N2,CO2,H2,H2O,Air,SO2)
! ... expressed in grams per mole
!
      gmw(1) = 32.D0              ! O2
      gmw(2) = 28.D0              ! N2
      gmw(3) = 44.D0              ! CO2
      gmw(4) = 2.D0               ! H2
      gmw(5) = 18.D0              ! H2O
      gmw(6) = 28.96442D0         ! Air
      gmw(7) = 64.D0              ! SO2
!
! ... international unit system (MKS)
      gmw = gmw * 1.D-3
!
! ... Lennard-Jones potential 
! ... sigma (Angstroms)
! ... epsilon/k (Kelvin/erg)
!
      mmugs(1) = 3.467D0          ! O2
      mmugek(1) = 106.7D0
      mmugs(2) = 3.798D0          ! N2
      mmugek(2) = 71.4D0
      mmugs(3) = 3.941D0          ! CO2
      mmugek(3) = 195.2D0
      mmugs(4) = 2.827D0          ! H2
      mmugek(4) = 59.7D0
      mmugs(5) = 2.641D0          ! H2O
      mmugek(5) = 809.1D0
      mmugs(6) = 3.711D0          ! Air
      mmugek(6) = 78.6D0
      mmugs(7) = 4.112D0          ! SO2
      mmugek(7) = 335.4D0
!
! Syamlal's particle-particle interaction coefficients
!
      CALL particles_constants_set
!
! set useful constants 
!
      gammaair = 1.4D0                           ! air adiabatic constant
      gamn     = (gammaair - 1.D0) / gammaair    ! useful constant ! 
      c_joule = 4.186D0                          ! cal per joule
      c_erg   = 4.186D7                          ! cal per erg
      tzero  = 0.0D0                             ! costant reference temperature
      hzerog = 0.D0                              ! gas enthalpy at tzero
      hzeros = 0.D0                              ! particles enthalpy at tzero
      rgas  = 8.31432D0                          ! perfect gas constant
      h1 = 0.D0                                  !
      h2 = 0.D0                                  !
      h3 = 0.D0                                  ! reaction enthalpies
      h4 = 0.D0                                  !
      h5 = 0.D0                                  !
!
      RETURN
      END SUBROUTINE setc
!----------------------------------------------------------------------
      SUBROUTINE gas_check
      USE atmospheric_conditions, ONLY: atm_ygc
      USE dimensions
      USE dome_conditions, ONLY: dome_ygc, idome
      USE gas_constants, ONLY: gas_type, present_gas
      USE eos_gas, ONLY: ygc
      USE time_parameters, ONLY: itd
      USE vent_conditions, ONLY: vent_ygc, ivent
      USE parallel, ONLY: mpime, root
      IMPLICIT NONE
      INTEGER :: ig, igg
      
      present_gas = .FALSE.

      IF (itd /= 2) THEN

        ig = 0
        DO igg = 1, max_ngas

          ! ... check gas species in specified-flow blocks
          !
          IF (no > 0 ) THEN
            IF( ANY(ygcob(igg,:) /= 0.0) ) THEN
              present_gas(igg) = .TRUE.
            END IF
          END IF 
          !
          ! ... check gas species in atmosphere and inlet
          !
          IF (atm_ygc(igg) /= 0.D0)  present_gas(igg) = .TRUE.
          IF (ivent > 0 .AND. vent_ygc(igg) /= 0.D0) present_gas(igg) = .TRUE.
          IF (idome > 0 .AND. dome_ygc(igg) /= 0.D0) present_gas(igg) = .TRUE.

          IF (present_gas(igg)) THEN
            ig = ig + 1
            gas_type(ig) = igg
          END IF
        END DO

        IF (ig /= ngas) CALL error('setup','wrong number of gas species',ig)

        DO ig = 1, ngas
          IF( mpime == root ) THEN
!            WRITE(logunit,*) ' Gas ', ig, ' is type ', gas_type(ig)
            IF (gas_type(ig)==0) &
              CALL error('gas_check','control gas species',1)
          END IF
        END DO

      END IF

      RETURN
      END SUBROUTINE gas_check
!----------------------------------------------------------------------
      END MODULE initial_conditions
!----------------------------------------------------------------------
