!----------------------------------------------------------------------
      MODULE initial_conditions
!----------------------------------------------------------------------
      IMPLICIT NONE
      SAVE
      REAL*8, DIMENSION(:), ALLOCATABLE   :: ugob, vgob, epob,  tgob, pob, ygc0
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: upob, vpob, epsob, tpob, ygcob

      INTEGER :: lpr
      REAL*8 :: zzero
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE bounds_setup
      USE dimensions
      IMPLICIT NONE
!
       ALLOCATE(ugob(no), vgob(no), pob(no), epob(no), tgob(no))
       ALLOCATE(upob(nsolid,no), vpob(nsolid,no), epsob(nsolid,no), tpob(nsolid,no))
       ALLOCATE(ygc0(ngas))
       ALLOCATE(ygcob(ngas,no))
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE setup
!
      USE atmosphere, ONLY: u0, v0, p0, temp0, uk0, vk0, ep0, atm
      USE dimensions
      USE eos_gas, ONLY: mole, cnvertg, gc_molar_fraction, gc_mass_fraction
      USE eos_solid, ONLY: cnverts
      USE gas_constants, ONLY: gmw, rgas
      USE gas_solid_density, ONLY: gas_bulk_density, solid_bulk_density
      USE gas_solid_velocity, ONLY: gas_velocity_r, gas_velocity_z
      USE gas_solid_velocity, ONLY: solid_velocity_r, solid_velocity_z
      USE gas_solid_temperature, ONLY: gas_temperature, solid_temperature
      USE grid, ONLY: grid_setup, zb, dz, dr
      USE grid, ONLY: fl, iob, nso
      USE particles_constants, ONLY: rl
      USE pressure_epsilon, ONLY: gas_pressure, void_fraction
      USE time_parameters, ONLY: itd
      USE gas_solid_viscosity, ONLY: gas_viscosity, gas_thermal_conductivity

      IMPLICIT NONE
!
      INTEGER :: i, j, k, n, j1, j2, i1, i2, ikpr, kpr, ij
      INTEGER :: kg
      REAL*8 :: zrif
!
      CALL grid_setup(zzero)
      CALL setc
!
      IF(itd.LE.1) THEN
!
! ... Set initial ambient pressure and temperature
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
            DO kg=1,ngas
              gc_mass_fraction(kg,ij)=ygc0(kg)
            END DO
            DO k=1,nsolid
              solid_bulk_density(k,ij)=rl(k)*(1.D0-ep0)/DBLE(nsolid)
              solid_temperature(k,ij)=gas_temperature(ij)
            END DO
!
! ... Set velocity profiles
!
            IF(fl(ij).EQ.1 .OR. fl(ij).EQ.4) THEN
             gas_velocity_r(ij)=u0
             gas_velocity_z(ij)=v0
             DO k=1,nsolid
               solid_velocity_r(k,ij)=uk0
               solid_velocity_z(k,ij)=vk0
             END DO
            END IF
!
          END DO
        END DO
! 
! ... Set initial conditions in Boundary cells with imposed fluid flow
!
        DO n=1,no
          i1=iob(1,n)
          i2=iob(2,n)
          j1=iob(3,n)
          j2=iob(4,n)
          DO j=j1,j2
          DO i=i1,i2
            ij=i+(j-1)*nr
            IF(nso(n).EQ.1 .OR. nso(n).EQ.5) THEN
              gas_velocity_r(ij)=ugob(n)
              gas_velocity_z(ij)=vgob(n)
              gas_temperature(ij)=tgob(n)
              gas_pressure(ij)=pob(n)
              void_fraction(ij)=epob(n)
              DO kg=1,ngas
                gc_mass_fraction(kg,ij)=ygcob(kg,n)
              END DO
              DO k=1,nsolid
                solid_temperature(k,ij)=tpob(k,n)
                solid_velocity_r(k,ij)=upob(k,n)
                solid_velocity_z(k,ij)=vpob(k,n)
                solid_bulk_density(k,ij)=epsob(k,n)*rl(k)
              END DO
            ENDIF
          END DO
          END DO
        END DO
!
! ... Compute thermodynamic quantities
!
        DO  j=1,nz
          DO i=1,nr
           ij=i+(j-1)*nr
           CALL mole(gc_molar_fraction(:,ij), gc_mass_fraction(:,ij))
           CALL cnvertg(ij)
           CALL cnverts(ij)
          END DO
        END DO
!
      END IF
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      END MODULE initial_conditions
!----------------------------------------------------------------------
