!----------------------------------------------------------------------
      MODULE initial_conditions
!----------------------------------------------------------------------
      IMPLICIT NONE
      SAVE
      REAL*8, DIMENSION(:), ALLOCATABLE   :: ugob, vgob, wgob, epob, &
                                             tgob, pob, ygc0
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: upob, vpob, wpob, epsob, &
                                             tpob, ygcob
      REAL*8 :: ugpr, wgpr, ppr, eppr, tgpr
      REAL*8, DIMENSION(:), ALLOCATABLE :: uppr,wppr,epspr,tppr, ygcpr

      INTEGER :: lpr
      REAL*8 :: zzero
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE allocate_setup
      USE dimensions
      IMPLICIT NONE
!
       ALLOCATE(ugob(no), vgob(no), wgob(no), pob(no), epob(no), tgob(no))
       ALLOCATE(upob(nsolid,no), vpob(nsolid,no), wpob(nsolid,no),        &
                epsob(nsolid,no), tpob(nsolid,no))
       ALLOCATE(ygc0(ngas))
       ALLOCATE(ygcob(ngas,no))

       ALLOCATE(ygcpr(ngas))
       ALLOCATE(uppr(nsolid), wppr(nsolid), epspr(nsolid), tppr(nsolid))
       
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE setup
! ... (2D/3D_Compliant)
!
      USE atmosphere, ONLY: u0, v0, w0, p0, temp0, us0, vs0, ws0, ep0
      USE atmosphere, ONLY: atm, controlatm
      USE control_flags, ONLY: job_type
      USE dimensions
      USE domain_decomposition, ONLY: ncint, meshinds, myijk
      USE eos_gas, ONLY: mas, mole, cnvertg, xgc, ygc
      USE eos_solid, ONLY: cnverts
      USE gas_constants, ONLY: default_gas, present_gas
      USE gas_constants, ONLY: gmw, rgas
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
      INTEGER :: x1, x2, z1, z2
      INTEGER :: ig, is
      REAL*8 :: zrif, prif, trif
      REAL*8 :: ymd, tem
!
      CALL grid_setup( zzero )
      CALL setc
!
! ... set initial conditions from prescribed input data
!
      IF (itd <= 1) THEN
!
! ... Control that atmospheric stratification is consistent
        CALL controlatm
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
            ygc(ig,ijk) = ygc0(ig)
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
          IF ( iob(n)%typ == 1 .OR. iob(n)%typ == 5 ) THEN

            tem = pob(n)/(rgas*tgob(n))*gmw(6)

            DO ijk = 1, ncint
              imesh = myijk( ip0_jp0_kp0_ , ijk )
              CALL meshinds(ijk,imesh,i,j,k)
 
              IF ( job_type == '2D' ) THEN

                IF ( k >= iob(n)%zlo .AND. k <= iob(n)%zhi  ) THEN
                  IF ( i >= iob(n)%xlo .AND. i <= iob(n)%xhi  ) THEN
                    ug(ijk) = ugob(n)
                    wg(ijk) = wgob(n)
                    tg(ijk) = tgob(n)
                    p(ijk)  = pob(n)
                    ep(ijk) = epob(n)
                    DO ig = 1, ngas
                      ygc(ig,ijk) = ygcob(ig,n)
                    END DO
                    DO is = 1,nsolid
                      ts(ijk,is)  = tpob(is,n)
                      us(ijk,is)  = upob(is,n)
                      ws(ijk,is)  = wpob(is,n)
                      rlk(ijk,is) = epsob(is,n)*rl(is)
                    END DO
                  END IF
                END IF

              ELSE IF ( job_type == '3D' ) THEN

                IF ( k >= iob(n)%zlo .AND. k <= iob(n)%zhi ) THEN
                  IF ( j >= iob(n)%ylo .AND. j <= iob(n)%yhi ) THEN
                    IF ( i >= iob(n)%xlo .AND. i <= iob(n)%xhi ) THEN
                      ug(ijk) = ugob(n)  ! 
                      vg(ijk) = vgob(n)
                      wg(ijk) = wgob(n)
                      tg(ijk) = tgob(n)
                      p(ijk)  = pob(n)
                      ep(ijk) = epob(n)
                      DO ig = 1, ngas
                        ygc(ig,ijk) = ygcob(ig,n)
                      END DO
                      DO is = 1,nsolid
                        ts(ijk,is)  = tpob(is,n)
                        us(ijk,is)  = upob(is,n)
                        vs(ijk,is)  = vpob(is,n)
                        ws(ijk,is)  = wpob(is,n)
                        rlk(ijk,is) = epsob(is,n)*rl(is)
                      END DO
                    ENDIF
                  END IF
                END IF

              END IF
            END DO 

          ELSE IF ( iob(n)%typ == 7) THEN
! ... Assign vertical profile 
!
            IF ( job_type == '2D' ) THEN

              READ(17,*) x1,x2,z1,z2
              IF ( (z2-z1) /= (iob(n)%zhi - iob(n)%zlo) .OR. (x1/=x2) )          &
                CALL error('setup','Error in input profile, block:', n)

              DO ijk = 1, ncint
                imesh = myijk( ip0_jp0_kp0_ , ijk )
                CALL meshinds(ijk,imesh,i,j,k)
 
                IF ( k >= iob(n)%zlo .AND. k <= iob(n)%zhi ) THEN
                  IF ( i >= iob(n)%xlo .AND. i <= iob(n)%xhi  ) THEN
                    IF ( i == iob(n)%xlo )                                         &
                    READ(17,*) ugpr,wgpr,ppr,eppr,tgpr,                            &
                              (uppr(is),wppr(is),epspr(is),tppr(is), is=1,nsolid), &
                              (ygcpr(ig), ig = 1,ngas)
                    eppr = 1.D0 - SUM(epspr)
                    ug(ijk)=ugpr
                    wg(ijk)=wgpr
                    tg(ijk)=tgpr+273.15
                    p(ijk)=ppr
                    ep(ijk)=eppr
                    ymd = 0.D0
                    DO ig=1,ngas
                      IF (ig /= default_gas) ymd = ymd + ygcpr(ig)
                      ygc(ig,ijk) = ygcpr(ig)
                      IF ( ygcpr(ig) /= 0.0 ) present_gas(ig) = .TRUE.
                    END DO
                    IF (.NOT.present_gas(default_gas))                      &
                      WRITE(*,*) 'default gas is not present'
                    ygc(default_gas,ijk) = 1.D0 - ymd
                    DO is=1,nsolid
                      ts(ijk,is)=tppr(is)+273.15
                      us(ijk,is)=uppr(is)
                      ws(ijk,is)=wppr(is)
                      rlk(ijk,is)=epspr(is)*rl(is)
                    END DO
                    CALL mole( xgc(:,ijk), ygc(:,ijk) )
                    CALL cnvertg(ijk)
                    CALL cnverts(ijk)
                    IF (i == 1) THEN
                      fl_l(ijk) = 5
                    ELSE
                      fl_l(ijk) = 1
                    END IF
                  END IF
                END IF
              END DO

            ELSE IF ( job_type == '3D' ) THEN
              CALL error('setup','Flow profile can be specified only in 2D',1)
            END IF
          END IF
        END DO 
!
! ... Compute thermodynamic quantities
!
        DO  ijk = 1, ncint
          CALL mole( xgc(:,ijk), ygc(:,ijk) )
          CALL cnvertg( ijk )
          CALL cnverts( ijk )
        END DO
!
! ... initial conditions already set from restart file
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
          xgc(default_gas,ijk) = 1.D0 
          DO ig = 1, ngas
            IF (ig /= default_gas) THEN
              xgc(default_gas,ijk) = xgc(default_gas,ijk) - xgc(ig,ijk)
            END IF
          END DO
          DO is=1,nsolid
            ep(ijk) = ep(ijk) - rlk(ijk,is)*inrl(is)
          END DO

          CALL mas( ygc(:,ijk), xgc(:,ijk)) 
          CALL cnvertg(ijk)
          CALL cnverts(ijk)

        END DO
!
      END IF
!
      RETURN
      END SUBROUTINE setup
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
      gammaair = 1.4D0                          ! air adiabatic constant
      gamn     = (gammaair - 1.D0) / gammaair    ! useful constant ! 
      c_joule = 4.186D0                          ! cal per joule
      c_erg   = 4.186D7                          ! cal per erg
      tzero  = 0.D0                              ! costant reference temperature
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
      END MODULE initial_conditions
!----------------------------------------------------------------------
