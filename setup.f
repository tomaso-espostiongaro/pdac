!----------------------------------------------------------------------
      MODULE initial_conditions
!----------------------------------------------------------------------
      IMPLICIT NONE
      SAVE
      REAL*8, DIMENSION(:), ALLOCATABLE   :: ugob, vgob, wgob, epob, &
                                             tgob, pob, ygc0
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: upob, vpob, wpob, epsob, &
                                             tpob, ygcob
      REAL*8, DIMENSION(:), ALLOCATABLE   :: ugpr, wgpr, ppr, eppr, tgpr
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: uppr,wppr,epspr,tppr, ygcpr
      INTEGER :: npr

      INTEGER :: lpr
      REAL*8 :: zzero
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE allocate_setup
      USE dimensions
      IMPLICIT NONE
!
! ... Specified flow blocks
!
       ALLOCATE(ugob(no), vgob(no), wgob(no), pob(no), epob(no), tgob(no))
       ALLOCATE(upob(nsolid,no), vpob(nsolid,no), wpob(nsolid,no),        &
                epsob(nsolid,no), tpob(nsolid,no))
       ALLOCATE(ygc0(max_ngas))
       ALLOCATE(ygcob(max_ngas,no))
!
! ... Specified profile
!
       IF (npr > 0) THEN
         ALLOCATE(ugpr(npr), wgpr(npr), ppr(npr), eppr(npr), tgpr(npr))
         ALLOCATE(uppr(nsolid,npr), wppr(nsolid,npr))
         ALLOCATE(epspr(nsolid,npr), tppr(nsolid,npr))
         ALLOCATE(ygcpr(max_ngas,npr))
       END IF

      RETURN
      END SUBROUTINE

!----------------------------------------------------------------------

      SUBROUTINE setup
! ... Set initial conditions
! ... (2D/3D_Compliant and fully parallel)
!
      USE atmosphere, ONLY: u0, v0, w0, p0, temp0, us0, vs0, ws0, ep0
      USE atmosphere, ONLY: atm, controlatm
      USE control_flags, ONLY: job_type
      USE dimensions
      USE domain_decomposition, ONLY: ncint, meshinds, myijk
      USE eos_gas, ONLY: mas, mole, cnvertg, xgc, ygc
      USE eos_solid, ONLY: cnverts
      USE gas_constants, ONLY: gmw, rgas, gammaair
      USE gas_constants, ONLY: default_gas, gas_type
      USE gas_solid_density, ONLY: rgp, rlk
      USE gas_solid_temperature, ONLY: tg, ts
      USE gas_solid_temperature, ONLY: sieg
      USE gas_solid_velocity, ONLY: ug, wg, vg
      USE gas_solid_velocity, ONLY: us, vs, ws
      USE grid, ONLY: grid_setup, zb, dz
      USE grid, ONLY: fl_l, iob
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: ep, p
      USE time_parameters, ONLY: itd
      USE indijk_module, ONLY: ip0_jp0_kp0_

      IMPLICIT NONE
!
      INTEGER :: i, j, k, ijk, ikpr, kpr, n, imesh
      INTEGER :: ig, is
      REAL*8 :: zrif, prif, trif
      REAL*8 :: mass, tem
      LOGICAL :: density_specified
!
      CALL grid_setup( zzero )
      CALL setc
!
! ... Control that atmospheric stratification is consistent
      CALL controlatm
!
! ... Check gas species
      CALL gas_check
!
! ... set initial conditions from prescribed input data
!
      IF (itd <= 1) THEN
!
! ... Set initial atmospheric conditions in all computational 
! ... domain (including boundary cells)
!
        DO ijk = 1, ncint
          imesh = myijk( ip0_jp0_kp0_ , ijk )
          CALL meshinds(ijk,imesh,i,j,k)
          
          zrif=zb(k)+0.5D0*(dz(1)-dz(k))
            
          CALL atm( zrif, prif, trif )
            
          p( ijk ) = prif
          tg( ijk ) = trif
!
! ... Set initial gas composition and particle concentrations
!
          ep(ijk) = ep0
          DO ig = 1, ngas
            ygc(ig,ijk) = ygc0(gas_type(ig))
          END DO
          DO is = 1, nsolid
            rlk(ijk,is) = rl(is)*(1.D0-ep0) / nsolid
            ts(ijk,is)  = tg(ijk)
          END DO
!
! ... Set initial velocity profiles
!
          IF ( fl_l(ijk) == 1 .OR. fl_l(ijk) == 4 ) THEN
           ug(ijk) = u0
           wg(ijk) = w0
           DO is = 1, nsolid
             us(ijk,is) = us0
             ws(ijk,is) = ws0
           END DO
           IF ( job_type == '3D' ) THEN
             vg(ijk) = v0
             DO is = 1, nsolid
               vs(ijk,is) = vs0
             END DO
           END IF
          END IF
!
        END DO
! 
! ... Set initial conditions in boundary cells 
! ... with imposed fluid flow 
!
        DO n = 1, no

          SELECT CASE (iob(n)%typ)

          CASE (1,5) ! ... Specified Fluid Flow in a block

            CALL specified_flow(n)

          CASE (7) ! ... Assign vertical profile 
!
            IF ( job_type == '2D' ) THEN
             
              CALL specified_profile(n)

            ELSE IF ( job_type == '3D' ) THEN
            
              CALL error('setup','Flow profile can be specified only in 2D',1)
            
            END IF
            
          END SELECT
        
        END DO 
!
! ... Compute thermodynamic quantities
!
        density_specified = .FALSE.  ! specify density instead of !
                                     ! temperature                !
        DO  ijk = 1, ncint
!
          CALL mole( xgc(:,ijk), ygc(:,ijk) )

          IF (density_specified) THEN
            mass = 0.D0
            DO ig=1,ngas
              mass = mass + xgc(ig,ijk)*gmw(gas_type(ig))
            END DO
            tem  = p(ijk) * mass / rgas / tg(ijk)
            tg(ijk) = tem
          END IF

          CALL cnvertg( ijk )
          CALL cnverts( ijk )
        END DO
!
! ... initial conditions already set from RESTART file
!
      ELSE IF (itd == 2) THEN 
!
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



      SUBROUTINE resetup

! ... Set initial conditions
! ... (2D/3D_Compliant)
!
      USE control_flags, ONLY: job_type
      USE dimensions
      USE domain_decomposition, ONLY: ncint
      USE eos_gas, ONLY: mas, cnvertg, xgc, ygc
      USE eos_solid, ONLY: cnverts
      USE gas_constants, ONLY: default_gas, gas_type
      USE gas_solid_density, ONLY: rlk
      USE particles_constants, ONLY: inrl
      USE pressure_epsilon, ONLY: ep
      USE time_parameters, ONLY: itd

      IMPLICIT NONE
!
      INTEGER :: ijk
      INTEGER :: ig, is, dfg
      REAL*8  :: xgc_def
!
      IF (itd <= 1) THEN
!
        CONTINUE
!
! ... initial conditions already set from RESTART file
!
      ELSE IF (itd == 2) THEN 
!
        CONTINUE
!
! ... set initial conditions from OUTPUT file
!
      ELSE IF (itd >= 3) THEN 

        DO ijk = 1, ncint

          ep(ijk) = 1.D0
          xgc_def = 1.0d0
          DO ig = 1, ngas
            IF ( gas_type(ig) /= default_gas) THEN
              xgc_def = xgc_def - xgc(ig,ijk)
            ELSE
              dfg = ig
            END IF
          END DO
          xgc(dfg,ijk) = xgc_def

          DO is=1,nsolid
            ep(ijk) = ep(ijk) - rlk(ijk,is)*inrl(is)
          END DO

          CALL mas( ygc(:,ijk), xgc(:,ijk)) 

          CALL cnvertg(ijk)
          CALL cnverts(ijk)

        END DO

      END IF

      RETURN
      END SUBROUTINE resetup
!----------------------------------------------------------------------

!----------------------------------------------------------------------

      SUBROUTINE specified_flow(n)

      USE control_flags, ONLY: job_type
      USE dimensions
      USE domain_decomposition, ONLY: ncint, meshinds, myijk
      USE eos_gas, ONLY: ygc
      USE gas_constants, ONLY: gmw, rgas, gammaair, gas_type,&
      default_gas
      USE gas_solid_density, ONLY: rgp, rlk
      USE gas_solid_temperature, ONLY: tg, ts
      USE gas_solid_velocity, ONLY: ug, wg, vg
      USE gas_solid_velocity, ONLY: us, vs, ws
      USE grid, ONLY: iob
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: ep, p
      USE indijk_module, ONLY: ip0_jp0_kp0_
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n
      REAL*8 :: ygcsum, ygcdfg
      INTEGER :: ijk,i,j,k,imesh
      INTEGER :: ig, is, dfg

        !csnd = DSQRT(gammaair*rgas/gmw(6)*tgob(n))

        DO ijk = 1, ncint
          imesh = myijk( ip0_jp0_kp0_ , ijk )
          CALL meshinds(ijk,imesh,i,j,k)

          IF ( k >= iob(n)%zlo .AND. k <= iob(n)%zhi  ) THEN
            IF ( j >= iob(n)%ylo .AND. j <= iob(n)%yhi    &
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
!
                ygcsum = SUM(ygc(:,ijk))
                ygcdfg = 1.D0
                IF ( ygcsum /= 1.D0 ) THEN
                  DO ig = 1, ngas
                    IF (gas_type(ig) /= default_gas) THEN
                      ygcdfg = ygcdfg - ygc(ig,ijk) 
                    ELSE
                      dfg = ig
                    END IF
                  END DO
                  ygc(dfg,ijk) = ygcdfg
                END IF
!
              END IF
            END IF
          END IF

            END DO 
      END SUBROUTINE specified_flow
!----------------------------------------------------------------------
      SUBROUTINE specified_profile(n)

      USE dimensions
      USE domain_decomposition, ONLY: ncint, meshinds, myijk
      USE eos_gas, ONLY: ygc, xgc, mole, cnvertg
      USE eos_solid, ONLY: cnverts
      USE gas_constants, ONLY: rgas, gammaair
      USE gas_constants, ONLY: default_gas, present_gas, gas_type
      USE gas_solid_density, ONLY: rgp, rlk
      USE gas_solid_temperature, ONLY: tg, ts
      USE gas_solid_velocity, ONLY: ug, wg, vg
      USE gas_solid_velocity, ONLY: us, vs, ws
      USE grid, ONLY: fl_l, iob
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: ep, p
      USE indijk_module, ONLY: ip0_jp0_kp0_
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n

      REAL*8 :: ygcsum, ygcdfg
      INTEGER :: ijk,i,j,k,imesh
      INTEGER :: ig, is, np, dfg

        DO ijk = 1, ncint
          imesh = myijk( ip0_jp0_kp0_ , ijk )
          CALL meshinds(ijk,imesh,i,j,k)

          IF ( k >= iob(n)%zlo .AND. k <= iob(n)%zhi ) THEN
            IF ( i >= iob(n)%xlo .AND. i <= iob(n)%xhi  ) THEN
              np = k - iob(n)%zlo + 1 

              ug(ijk) = ugpr(np)
              wg(ijk) = wgpr(np)
              tg(ijk) = tgpr(np)+273.15
              p(ijk)  = ppr(np)
              ep(ijk) = 1.D0 - SUM(epspr(:,np))

              DO is=1,nsolid
                ts(ijk,is)=tppr(is,np)+273.15
                us(ijk,is)=uppr(is,np)
                ws(ijk,is)=wppr(is,np)
                rlk(ijk,is)=epspr(is,np)*rl(is)
              END DO

              DO ig=1,ngas
                ygc(ig,ijk) = ygcpr(gas_type(ig),np)
              END DO
!
! ... check gas components closure relation
!
              ygcsum = SUM(ygc(:,ijk))
              ygcdfg = 1.D0
              IF ( ygcsum /= 1.D0 ) THEN
                DO ig = 1, ngas
                  IF (gas_type(ig) /= default_gas) THEN
                    ygcdfg = ygcdfg - ygc(ig,ijk) 
                  ELSE
                    dfg = ig
                  END IF
                END DO
                ygc(dfg,ijk) = ygcdfg
              END IF
!
              IF (i == 1) THEN
                fl_l(ijk) = 5
              ELSE
                fl_l(ijk) = 1
              END IF
!
            END IF
          END IF
        END DO

      END SUBROUTINE specified_profile
!----------------------------------------------------------------------
      SUBROUTINE setc
!
      USE dimensions
      USE gas_constants, ONLY: gmw, mmugs, mmugek,     &
     &    gammaair, gamn, c_joule, c_erg, tzero, hzerog, hzeros, rgas
      USE particles_constants, ONLY: particles_constants_set, cmus
      USE reactions, ONLY: h1, h2, h3, h4, h5
      USE turbulence_model, ONLY: turbulence_setup, iturb
      USE gas_solid_viscosity, ONLY: mus
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
! Initialize particle's viscosity
! 
      DO is = 1, nsolid
        mus(:,is) = cmus(is)
      END DO
!
! set useful constants 
!
      gammaair = 1.33D0                          ! air adiabatic constant
      gamn     = (gammaair - 1.D0) / gammaair    ! useful constant ! 
      c_joule = 4.186D0                          ! cal per joule
      c_erg   = 4.186D7                          ! cal per erg
      tzero  = 0.0D0                             ! costant reference temperature
      hzerog = 0.D0                              ! gas enthalpy at tzero
      hzeros = 0.D0                              ! particles enthalpy at tzero
      rgas  = 8.31432                            ! perfect gas constant
      h1 = 0.D0                                  !
      h2 = 0.D0                                  !
      h3 = 0.D0                                  ! reaction enthalpies
      h4 = 0.D0                                  !
      h5 = 0.D0                                  !
!
      IF (iturb >= 1) THEN
        CALL turbulence_setup
      END IF
!
      RETURN
      END SUBROUTINE setc
!----------------------------------------------------------------------
      SUBROUTINE gas_check
      USE dimensions
      USE gas_constants, ONLY: gas_type, present_gas, default_gas
      IMPLICIT NONE
      INTEGER :: ig, igg
      
      ig = 0
      DO igg = 1, max_ngas
          IF ((ygc0(igg) /= 0.0) .OR. ANY(ygcob(igg,:) /= 0.0) ) THEN
            ig = ig + 1
            gas_type(ig) = igg
            present_gas(igg) = .TRUE.
          END IF
          IF (npr > 0) THEN
            IF( ANY(ygcpr(igg,:) /= 0.0) ) THEN
              ig = ig + 1
              gas_type(ig) = igg
              present_gas(igg) = .TRUE.
            END IF
          END IF
      END DO
      IF (ig /= ngas) CALL error('setup','wrong number of gas species',ig)
      DO ig = 1, ngas
        WRITE(6,*) ' Gas ', ig, ' is type ', gas_type(ig)
      END DO

      RETURN
      END SUBROUTINE gas_check
!----------------------------------------------------------------------
      END MODULE initial_conditions
!----------------------------------------------------------------------
