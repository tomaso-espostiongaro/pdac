!----------------------------------------------------------------------
      MODULE initial_conditions
!----------------------------------------------------------------------
      IMPLICIT NONE
      SAVE
      REAL*8, DIMENSION(:), ALLOCATABLE   :: ugob, vgob, epob,  tgob,
     &        pob, ygc0
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: upob, vpob, epsob, tpob,
     &        ygcob

      INTEGER :: lpr
      REAL*8 :: zzero
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE bounds_setup
      USE dimensions
      IMPLICIT NONE
!
       ALLOCATE(ugob(nnso), vgob(nnso), pob(nnso), epob(nnso),
     &          tgob(nnso))
       ALLOCATE(upob(ncl,nnso), vpob(ncl,nnso), epsob(ncl,nnso),
     &          tpob(ncl,nnso))
       ALLOCATE(ygc0(ngas))
       ALLOCATE(ygcob(ngas,nnso))
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE setup
!
      USE atmosphere, ONLY: u0, v0, p0, temp0, uk0, vk0, ep0, atm
      USE dimensions
      USE eos_gas, ONLY: mole, cnvertg, xgc_g, ygc_g
      USE eos_solid, ONLY: cnverts
      USE gas_constants, ONLY: gmw, rgas
      USE gas_solid_density, ONLY: rgp_g, rlk_g
      USE gas_solid_velocity, ONLY: ug_g, vg_g, uk_g, vk_g
      USE gas_solid_temperature, ONLY: tg_g, tk_g
      USE grid, ONLY: grid_setup, zb, dz, dr, ib2, ib1, jb2, jb1
      USE grid, ONLY: fl, iob, nso, no
      USE particles_constants, ONLY: rl, nsolid
      USE pressure_epsilon, ONLY: p_g, ep_g
      USE time_parameters, ONLY: itd
      USE gas_solid_viscosity, ONLY: mug_g, kapg_g

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
! ... set initial atmospheric conditions (pressure/temperature stratification)
!
        DO  j=1,jb2
          zrif=zb(j)+0.5D0*(dz(1)-dz(j))
         DO i=1,ib2
            ij=i+(j-1)*ib2
            CALL atm(zrif,p_g(ij),tg_g(ij))
            rgp_g(ij)=p_g(ij)*ep0*gmw(6)/(rgas*tg_g(ij))
         END DO
        END DO
!
! ... set initial composition and particles concentration in atmosphere
!
        DO j=jb2,1,-1
        DO i=1,ib2
          ij=i+(j-1)*ib2
          ep_g(ij)=ep0
          DO kg=1,ngas
            ygc_g(kg,ij)=ygc0(kg)
          END DO
          DO k=1,nsolid
            rlk_g(k,ij)=rl(k)*(1.D0-ep0)/DBLE(nsolid)
            tk_g(k,ij)=tg_g(ij)
          END DO
          CALL mole(xgc_g(:,ij), ygc_g(:,ij))
          CALL cnvertg(ij)
          CALL cnverts(ij)
        END DO
        END DO
!
! ... Set the wind profile
!
        DO j=1,jb2
        DO i=1,ib2
          ij=i+(j-1)*ib2
          IF(fl(ij).EQ.1 .OR. fl(ij).EQ.4) THEN
           ug_g(ij)=u0
           vg_g(ij)=v0
           DO k=1,nsolid
             uk_g(k,ij)=uk0
             vk_g(k,ij)=vk0
           END DO
          END IF
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
            ij=i+(j-1)*ib2
            IF(nso(n).EQ.1 .OR. nso(n).EQ.5) THEN
              ug_g(ij)=ugob(n)
              vg_g(ij)=vgob(n)
              tg_g(ij)=tgob(n)
              p_g(ij)=pob(n)
              ep_g(ij)=epob(n)
              DO kg=1,ngas
                ygc_g(kg,ij)=ygcob(kg,n)
              END DO
              DO k=1,nsolid
                tk_g(k,ij)=tpob(k,n)
                uk_g(k,ij)=upob(k,n)
                vk_g(k,ij)=vpob(k,n)
                rlk_g(k,ij)=epsob(k,n)*rl(k)
              END DO
              CALL mole(xgc_g(:,ij), ygc_g(:,ij))
              CALL cnvertg(ij)
              CALL cnverts(ij)
            ENDIF
          END DO
          END DO
        END DO
!
! ... initialize molecular viscosity and thermal conductivity of gas
!
        mug_g = 0.0
        kapg_g = 0.0
!
      END IF
!
 650  FORMAT (1x,80i1)
 660  FORMAT (//)
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      END MODULE initial_conditions
!----------------------------------------------------------------------
