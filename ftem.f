!----------------------------------------------------------------------
      MODULE enthalpy_matrix
!----------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8, DIMENSION(:),   ALLOCATABLE :: bt
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: at
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: hv
!
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE ftem
!----------------------------------------------------------------------
! ... This routine computes the matrix elements for the thermal 
! ... interphase coupling, and solves for the enthalpies
!
      USE dimensions
      USE eos_gas, ONLY: cg
      USE gas_solid_density, ONLY: rgp, rgpn, rlk, rlkn, rog
      USE gas_solid_temperature, ONLY: sieg, siegn, sies, siesn
      USE gas_solid_velocity, ONLY: ug, wg, us, ws
      USE grid, ONLY: r, rb, dr, zb, dz, inr, indr, indz 
      USE grid, ONLY: fl_l
      USE grid, ONLY: ncint, data_exchange
      USE heat_transfer, ONLY: hvs
      USE tilde_energy, ONLY: rhg, rhk
      USE pressure_epsilon, ONLY: p, pn, ep
      USE reactions, ONLY: hrex, irex
      USE set_indexes
      USE heat_capacity, ONLY: ck
      USE time_parameters, ONLY: time, dt
      USE gas_solid_viscosity, ONLY: mug, kapg
      IMPLICIT NONE
!
      REAL*8 :: c3, hrexs, hrexg, c2, upxy, c1
      REAL*8 :: drc, dzc
      INTEGER :: is, m, l, is1
      INTEGER :: ij
!
      ALLOCATE(at(nphase, nphase))
      ALLOCATE(bt(nphase))
      ALLOCATE(hv(nsolid, ncint))
!
          at = 0.D0
          bt = 0.D0
          hv = 0.D0
!
!pdac------------
! Control heat reactions
      hrexg=0.D0
      hrexs=0.D0
!pdac------------
!
!
      DO ij = 1, ncint

        IF(fl_l(ij).EQ.1) THEN
          CALL subscr(ij)
!
          IF(irex.EQ.2) CALL hrex(ij,hrexg,hrexs)
!
          at(1,1) = rgp(ij)
          bt(1)   = siegn(ij) * rgpn(ij) + rhg(ij) - hrexg

          DO is=1, nsolid
            is1=is+1
!
! ... Compute gas-particle heat transfer coefficients
!
            CALL hvs(hv(is,ij), rlk(is,ij), rog(ij),       &
                 ep(ij), ug(ij), ug(imj), us(is,ij),      &
                 us(is,imj), wg(ij), wg(ijm),             &
                 ws(is,ij), ws(is,ijm), mug(ij),           &
                 kapg(ij), cg(ij), is)
!
            at(1,1)   = at(1,1)       + dt * hv(is,ij) / cg(ij)
            at(1,is1)  =               - dt * hv(is,ij) / ck(is,ij)
            at(is1,1)  =               - dt * hv(is,ij) / cg(ij)
            at(is1,is1) = rlk(is,ij) + dt * hv(is,ij) / ck(is,ij)
!
            bt(is1) = rlkn(is,ij) * siesn(is,ij) + rhk(is, ij)
          END DO
!
! ... Solve the interphase enthalpy matrix by using Gauss inversion
!
          CALL invdm(at, bt, ij)
!
          sieg(ij) = bt(1)
          DO is=1, nsolid
            sies(is,ij) = bt(is+1)
          END DO
!
        END IF

      END DO
!
      DEALLOCATE(at)
      DEALLOCATE(bt)
      DEALLOCATE(hv)
      DEALLOCATE(rhg)
      DEALLOCATE(rhk)
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE invdm(a, b, ij)
!----------------------------------------------------------------------
! ... this routine solves the N-Phases Energy-Equation Matrix
! ... (Gauss elimination)
!
      USE dimensions
      USE gas_solid_density, ONLY: rlk
      USE particles_constants, ONLY: rl, inrl
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: ij
      REAL*8 :: a(nphase,nphase),b(nphase)
!
      INTEGER ::is 
      REAL*8 :: div
!
      DO is=nphase,2,-1
!        IF(abs(a(is,is)).LE.1.D-6) THEN
        IF(rlk(is-1,ij)*inrl(is-1).LE.1.D-9) THEN
          a(1,is)=0.D0
          a(is,1)=0.D0
          b(is)=0.D0
        ELSE
!
! ... eliminate all cross elements in the gas enthalpy equation
!
          div=1.D0/a(is,is)
          a(is,1)=a(is,1)*div
          b(is)=b(is)*div
          b(1)=b(1)-a(1,is)*b(is)
          a(1,1)=a(1,1)-a(1,is)*a(is,1)
        ENDIF
      END DO
!
      b(1)=b(1)/a(1,1)
      DO is=2,nphase
        b(is)=b(is)-a(is,1)*b(1)
      END DO
!
      RETURN
      END SUBROUTINE invdm
!----------------------------------------------------------------------
      END MODULE enthalpy_matrix
