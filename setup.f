!----------------------------------------------------------------------
      MODULE initial_conditions
!----------------------------------------------------------------------
      IMPLICIT NONE
      SAVE
      REAL*8, DIMENSION(:), ALLOCATABLE   :: ugob, vgob, wgob, epob, &
                                             tgob, pob, ygc0
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: upob, vpob, wpob, epsob, &
                                             tpob, ygcob

      INTEGER :: lpr
      REAL*8 :: zzero
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE bounds_setup
      USE dimensions
      IMPLICIT NONE
!
       ALLOCATE(ugob(no), vgob(no), wgob(no), pob(no), epob(no), tgob(no))
       ALLOCATE(upob(nsolid,no), vpob(nsolid,no), wpob(nsolid,no), epsob(nsolid,no), tpob(nsolid,no))
       ALLOCATE(ygc0(ngas))
       ALLOCATE(ygcob(ngas,no))
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE setup
!
      USE atmosphere, ONLY: u0, w0, p0, temp0, us0, ws0, ep0, atm
      USE dimensions
      USE eos_gas, ONLY: mas, mole, cnvertg, gc_molar_fraction, gc_mass_fraction
      USE eos_solid, ONLY: cnverts
      USE gas_constants, ONLY: gmw, rgas
      USE gas_solid_density, ONLY: gas_bulk_density, solid_bulk_density
      USE gas_solid_velocity, ONLY: gas_velocity_r, gas_velocity_z
      USE gas_solid_velocity, ONLY: solid_velocity_r, solid_velocity_z
      USE gas_solid_temperature, ONLY: gas_temperature, solid_temperature
      USE gas_solid_temperature, ONLY: gas_enthalpy
      USE grid, ONLY: grid_setup, zb, dx, dy, dz, dr
      USE grid, ONLY: fl, iob
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: gas_pressure, void_fraction
      USE time_parameters, ONLY: itd
      USE gas_solid_viscosity, ONLY: gas_viscosity, gas_thermal_conductivity

      IMPLICIT NONE
!
      INTEGER :: i, j, j1, j2, i1, i2, ikpr, kpr, ij, n
      INTEGER :: ig, is, imesh
      REAL*8 :: zrif
      REAL*8 :: entemp
!
      CALL grid_setup(zzero)
      CALL setc
!
! ... MODIFICARE_X3D

      IF (itd <= 1) THEN

!
! ... Set initial ambient pressure and temperature within computational domain
! ... and boundary cells
!
        DO  j=1,nz
          zrif=zb(j)+0.5D0*(dz(1)-dz(j))
          DO i=1,nr
            ij=i+(j-1)*nr
            CALL atm(zrif,gas_pressure(ij),gas_temperature(ij))
!
! ... Set initial gas composition and particles concentration
!
            void_fraction(ij)=ep0
            DO ig=1,ngas
              gc_mass_fraction(ig,ij)=ygc0(ig)
            END DO
            DO is=1,nsolid
              solid_bulk_density(is,ij)=rl(is)*(1.D0-ep0)/DBLE(nsolid)
              solid_temperature(is,ij)=gas_temperature(ij)
            END DO
!
! ... Set velocity profiles
!
            IF(fl(ij).EQ.1 .OR. fl(ij).EQ.4) THEN
             gas_velocity_r(ij)=u0
             gas_velocity_z(ij)=w0
             DO is=1,nsolid
               solid_velocity_r(is,ij)=us0
               solid_velocity_z(is,ij)=ws0
             END DO
            END IF
!
          END DO
        END DO
! 
! ... Set initial conditions in Boundary cells with imposed fluid flow
!
        DO n=1,no
          DO j = iob(n)%zlo, iob(n)%zhi 
          DO i = iob(n)%rlo, iob(n)%rhi
            ij=i+(j-1)*nr
            IF( iob(n)%typ == 1 .OR. iob(n)%typ == 5 ) THEN
              gas_velocity_r(ij)=ugob(n)
              gas_velocity_z(ij)=wgob(n)
              gas_temperature(ij)=tgob(n)
              gas_pressure(ij)=pob(n)
              void_fraction(ij)=epob(n)
              DO ig=1,ngas
                gc_mass_fraction(ig,ij)=ygcob(ig,n)
              END DO
              DO is=1,nsolid
                solid_temperature(is,ij)=tpob(is,n)
                solid_velocity_r(is,ij)=upob(is,n)
                solid_velocity_z(is,ij)=wpob(is,n)
                solid_bulk_density(is,ij)=epsob(is,n)*rl(is)
              END DO
            ENDIF
          END DO
          END DO
        END DO
!
! ... Compute thermodynamic quantities
!
        DO  imesh = 1, ntot
          CALL mole(gc_molar_fraction(:,imesh), gc_mass_fraction(:,imesh))
          CALL cnvertg(imesh)
          CALL cnverts(imesh)
        END DO
!
      ELSE IF (itd == 2) THEN 
!
        CONTINUE
!
      ELSE IF (itd >= 3) THEN 

        void_fraction = 0.D0
        gc_molar_fraction(6,:) = 1.D0 -  gc_molar_fraction(5,:)
        DO imesh = 1, ntot
          DO is=1,nsolid
            void_fraction(imesh) = 1.D0 - solid_bulk_density(is,imesh)*inrl(is)
          END DO
          CALL mas(gc_mass_fraction(:,imesh), gc_molar_fraction(:,imesh)) 
          CALL cnvertg(imesh)
          CALL cnverts(imesh)
        END DO
!
      END IF
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      END MODULE initial_conditions
!----------------------------------------------------------------------
