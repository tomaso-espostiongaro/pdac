!----------------------------------------------------------------------
      MODULE initial_conditions
!----------------------------------------------------------------------
      IMPLICIT NONE
      SAVE
      REAL*8, DIMENSION(:), ALLOCATABLE   :: ugob, vgob, wgob, epob, &
                                             tgob, pob
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: upob, vpob, wpob, epsob, &
                                             tpob, ygcob
      REAL*8, DIMENSION(:), ALLOCATABLE   :: ugpr, wgpr, ppr, eppr, tgpr
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: uppr,wppr,epspr,tppr, ygcpr
      LOGICAL :: density_specified

      INTEGER :: npr
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
      USE atmospheric_conditions, ONLY: wind_x, wind_y, wind_z, void_fraction
      USE atmospheric_conditions, ONLY: p_atm, t_atm, atm_ygc
      USE control_flags, ONLY: job_type
      USE dimensions
      USE domain_decomposition, ONLY: ncint, meshinds, myijk
      USE eos_gas, ONLY: mas, mole, xgc, ygc
      USE gas_constants, ONLY: gmw, rgas, gammaair
      USE gas_constants, ONLY: gas_type
      USE gas_solid_density, ONLY: rlk
      USE gas_solid_temperature, ONLY: tg, ts
      USE gas_solid_temperature, ONLY: sieg
      USE gas_solid_velocity, ONLY: ug, wg, vg
      USE gas_solid_velocity, ONLY: us, vs, ws
      USE gas_solid_viscosity, ONLY: mus
      USE grid, ONLY: zb, dz
      USE grid, ONLY: flag, iob
      USE immersed_boundaries, ONLY: numx, numy, numz, immb
      USE particles_constants, ONLY: rl, inrl, cmus
      USE pressure_epsilon, ONLY: ep, p
      USE time_parameters, ONLY: itd
      USE turbulence_model, ONLY: turbulence_setup, iturb
      USE vent_conditions, ONLY: set_ventc, ivent
      USE indijk_module, ONLY: ip0_jp0_kp0_

      IMPLICIT NONE
!
      INTEGER :: i, j, k, ijk, ikpr, kpr, n, imesh
      INTEGER :: ig, is, fx, fy, fz
      REAL*8 :: zrif
      REAL*8 :: mass, tem
      LOGICAL :: forced = .FALSE.
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
! ... Set initial conditions from prescribed input data
!
      IF (itd <= 1) THEN
!
! ... Set initial atmospheric conditions in all computational 
! ... domain (---> including boundary cells <---)
!
        DO ijk = 1, ncint
          CALL meshinds(ijk,imesh,i,j,k)
          
          zrif=zb(k)+0.5D0*(dz(1)-dz(k))

          p( ijk ) = p_atm(k)
          tg( ijk ) = t_atm(k)

          ! ... Set initial gas composition, particle concentrations 
          ! .... and temperature
          !
          ep(ijk) = void_fraction

          DO ig = 1, ngas 
            ygc(ig,ijk) = atm_ygc(gas_type(ig))
          END DO

          DO is = 1, nsolid
            rlk(ijk,is) = rl(is)*(1.D0-ep(ijk)) / nsolid
            ts(ijk,is)  = tg(ijk)
          END DO

          ! ... Set initial velocity profiles
          !
          SELECT CASE ( flag(ijk) )

          CASE (1,4,6)

            IF (immb >= 1) THEN
              fx = numx(ijk)
              IF (job_type == '2D') THEN
                fy = 0
              ELSE IF (job_type == '3D') THEN
                fy = numy(ijk)
              END IF
              fz = numz(ijk)
              forced = (fx/=0 .OR. fy/=0 .OR. fz/=0)
            END IF

            IF( .NOT.forced ) THEN
              ug(ijk) = wind_x
              wg(ijk) = wind_z
              us(ijk,:) = ug(ijk)
              ws(ijk,:) = wg(ijk)
              IF ( job_type == '3D' ) THEN
                vg(ijk) = wind_y
                vs(ijk,:) = vg(ijk)
              END IF
            END IF

          CASE DEFAULT

            CONTINUE

          END SELECT
!
        END DO
! 
! ... Set initial conditions in boundary cells 
! ... with specified fluid flow 
!
        IF (ivent >= 1) CALL set_ventc

        DO n = 1, no

          SELECT CASE (iob(n)%typ)

          CASE (1,5) ! ... Specified Fluid Flow in a block

            CALL specified_flow(n)

          CASE (3) ! ... Obstacle

            CONTINUE

          END SELECT
        
        END DO 
! 
! ... initial conditions already set from RESTART file
!
      ELSE IF (itd == 2) THEN 

        CONTINUE
!
! ... set initial conditions from OUTPUT file
!
      ELSE IF (itd >= 3) THEN 

        CONTINUE
!
      END IF
!
      RETURN
      END SUBROUTINE setup
!----------------------------------------------------------------------
      SUBROUTINE cnvert
!
! ... Compute derived thermodynamic quantities from initial conditions
! ... in the whole computational domain (--> including boundaries <--)
! ... (2D/3D_Compliant)
!
      USE control_flags, ONLY: job_type
      USE dimensions
      USE domain_decomposition, ONLY: ncint
      USE eos_gas, ONLY: mas, xgc, ygc, rgpgc
      USE eos_gas, ONLY: cg, caloric_eosg, thermal_eosg, mole
      USE eos_solid, ONLY: caloric_eosl
      USE gas_constants, ONLY: gas_type, gmw, rgas, tzero, hzerog, hzeros
      USE gas_solid_density, ONLY: rlk, rog, rgp
      USE gas_solid_temperature, ONLY: tg, ts, sieg, sies
      USE grid, ONLY: flag
      USE io_restart, ONLY: dump_all
      USE particles_constants, ONLY: inrl, cps
      USE pressure_epsilon, ONLY: ep, p
      USE time_parameters, ONLY: itd
      USE specific_heat_module, ONLY: ck, cp, hcapg, hcaps

      IMPLICIT NONE
!
      INTEGER :: ijk, info
      INTEGER :: ig, is, dfg
      REAL*8  :: xgc_def, mass, tem, rls
      REAL*8  :: hc, mg , xg(1:max_ngas), eps
!
      IF (itd <= 1) THEN
!
        DO  ijk = 1, ncint
!
          ! ... compute gas components molar fractions 
          ! ... from mass fractions
          !
          CALL mole( xgc(:,ijk), ygc(:,ijk) )
          
          ! ... If density is specified tg(ijk) is a density
          ! ... that must be converted into a temperature
          !
          IF (density_specified) THEN
            mass = 0.D0
            DO ig=1,ngas
              mass = mass + xgc(ig,ijk)*gmw(gas_type(ig))
            END DO
            tem  = p(ijk) * mass / rgas / tg(ijk)
            tg(ijk) = tem
          END IF

          ! ... compute gas density from thermal equation of state
          ! ... compute gas and gas components bulk densities
          !
          CALL thermal_eosg( rog(ijk), tg(ijk), p(ijk), xgc(:,ijk) )
          rgp(ijk) = rog(ijk) * ep(ijk)
          DO ig = 1, ngas
            rgpgc(ijk,ig) = ygc(ig,ijk) * rgp(ijk)
          END DO

          ! ... compute specific heat from temperature
          !
          CALL hcapg(cp(:,ijk), tg(ijk))
          hc = 0.D0
          DO ig = 1, ngas
            hc = hc + cp(gas_type(ig),ijk) * ygc(ig,ijk)
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
        IF( .NOT. dump_all ) THEN

          DO ijk = 1, ncint
            ! ... compute gas components molar fractions 
            ! ... from mass fractions
            !
            CALL mole( xgc(:,ijk), ygc(:,ijk) )

            ! ... compute void fraction
            !
            rls = 0.0d0
            DO is = 1, nsolid
              rls = rls + rlk(ijk,is) * inrl(is)
            END DO
            ep(ijk) = 1.D0 - rls

            ! ... gas components bulk densities
            !
            DO ig = 1, ngas
              rgpgc(ijk,ig) = ygc(ig,ijk) * rgp(ijk)
            END DO

            ! ... WARNING!: an error is introduced in gas temperature ...
            ! ... In the iterative inversion of the enthalpy equation 
            ! ... the initial value of temperature is different
            !
            CALL caloric_eosg(cp(:,ijk), cg(ijk), tg(ijk), ygc(:,ijk), &
                            sieg(ijk), ijk, info)
            DO is=1, nsolid
              CALL caloric_eosl(ts(ijk,is),cps(is),ck(is,ijk),sies(ijk,is),ijk)
            END DO
          END DO

        ELSE

          DO ijk = 1, ncint
            CALL hcapg(cp(:,ijk), tg(ijk))
            hc = 0.D0
            DO ig = 1, ngas
              hc = hc + cp(gas_type(ig),ijk) * ygc(ig,ijk)
            END DO
            cg(ijk) = hc
            DO is = 1, nsolid
              CALL hcaps(ck(is,ijk), cps(is), ts(ijk,is))
            END DO
          END DO

        END IF
!
      ELSE IF (itd >= 3) THEN 

        CONTINUE

      END IF

      RETURN
      END SUBROUTINE cnvert
!----------------------------------------------------------------------
      SUBROUTINE specified_flow(n)

      USE control_flags, ONLY: job_type
      USE dimensions
      USE domain_decomposition, ONLY: ncint, meshinds, myijk
      USE eos_gas, ONLY: ygc
      USE gas_constants, ONLY: gmw, rgas, gammaair, gas_type
      USE gas_solid_density, ONLY: rgp, rlk
      USE gas_solid_temperature, ONLY: tg, ts
      USE gas_solid_velocity, ONLY: ug, wg, vg
      USE gas_solid_velocity, ONLY: us, vs, ws
      USE grid, ONLY: iob, x, y, z
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: ep, p
      USE indijk_module, ONLY: ip0_jp0_kp0_
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n
      REAL*8 :: ygcsum, dist2
      INTEGER :: ijk,i,j,k,imesh
      INTEGER :: ig, is, dfg

        DO ijk = 1, ncint
          CALL meshinds(ijk,imesh,i,j,k)

          IF ( k >= iob(n)%zlo .AND. k <= iob(n)%zhi  ) THEN
            IF ( ( j >= iob(n)%ylo .AND. j <= iob(n)%yhi )    &
                                    .OR. job_type == '2D') THEN
              IF ( i >= iob(n)%xlo .AND. i <= iob(n)%xhi  ) THEN
                
                ug(ijk) = ugob(n)
                IF (job_type == '3D') vg(ijk) = vgob(n)
                wg(ijk) = wgob(n)
                tg(ijk) = tgob(n)
                p(ijk)  = pob(n)
                ep(ijk) = 1.D0 - SUM(epsob(:,n))
                DO ig = 1, ngas
                  ygc(ig,ijk) = ygcob(gas_type(ig),n)
                END DO
                DO is = 1,nsolid
                  ts(ijk,is)  = tpob(is,n)
                  us(ijk,is)  = upob(is,n)
                  IF (job_type == '3D') vs(ijk,is)  = vpob(is,n)
                  ws(ijk,is)  = wpob(is,n)
                  rlk(ijk,is) = epsob(is,n)*rl(is)
                END DO
                !
                ! ... check gas components closure relation
                ygcsum = SUM(ygc(:,ijk))
                IF ( ygcsum /= 1.D0 ) THEN
                  ygc(ngas,ijk) = 1.D0 - SUM( ygc(1:ngas-1,ijk) )
                END IF

              END IF
            END IF
          END IF
!
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
      USE gas_constants, ONLY: gas_type, present_gas
      USE eos_gas, ONLY: ygc
      USE time_parameters, ONLY: itd
      USE vent_conditions, ONLY: vent_ygc
      IMPLICIT NONE
      INTEGER :: ig, igg
      
      present_gas = .FALSE.

      IF (itd == 1) THEN

        ig = 0
        DO igg = 1, max_ngas

          ! ... check gas species in specified-flow blocks
          !
          IF (no > 0 ) THEN
            IF( ANY(ygcob(igg,:) /= 0.0) ) present_gas(igg) = .TRUE.
          END IF 
          !
          ! ... check gas species in atmosphere and inlet
          !
          IF ((atm_ygc(igg) /= 0.D0) .OR. (vent_ygc(igg) /= 0.D0) ) THEN
            present_gas(igg) = .TRUE.
          END IF

          IF (present_gas(igg)) THEN
            ig = ig + 1
            gas_type(ig) = igg
          END IF

        END DO
        IF (ig /= ngas) CALL error('setup','wrong number of gas species',ig)

      ELSE IF (itd == 2) THEN

        DO igg = 1, max_ngas
          IF (ANY(gas_type == igg)) present_gas(igg) = .TRUE.
        END DO

      END IF

      DO ig = 1, ngas
        WRITE(6,*) ' Gas ', ig, ' is type ', gas_type(ig)
      END DO

      RETURN
      END SUBROUTINE gas_check
!----------------------------------------------------------------------
      END MODULE initial_conditions
!----------------------------------------------------------------------
