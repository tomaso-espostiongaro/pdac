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
      USE gas_solid_temperature, ONLY: sieg, siegn, siek, siekn
      USE gas_solid_velocity, ONLY: ug, vg, uk, vk
      USE grid, ONLY: r, rb, dr, zb, dz, inr, indr, indz 
      USE grid, ONLY: fl_l
      USE grid, ONLY: nij_l, myij, data_exchange
      USE heat_transfer, ONLY: hvs
      USE tilde_energy, ONLY: rhg, rhk
      USE pressure_epsilon, ONLY: p, pn, ep
      USE reactions, ONLY: hrex, irex
      USE set_indexes
      USE th_capacity, ONLY: ck
      USE time_parameters, ONLY: time, dt
      USE gas_solid_viscosity, ONLY: mug, kapg
      IMPLICIT NONE
!
      REAL*8 :: c3, hrexs, hrexg, c2, upxy, c1
      REAL*8 :: drc, dzc
      INTEGER :: k, m, l, k1, i, j, ij_g
      INTEGER :: ij
!
      ALLOCATE(at(nphase, nphase))
      ALLOCATE(bt(nphase))
      ALLOCATE(hv(ncl, nij_l))
!
!pdac------------
! Control heat reactions
      hrexg=0.D0
      hrexs=0.D0
!pdac------------
!
!
      DO ij = 1, nij_l
        ij_g = myij(0, 0, ij)
        IF(fl_l(ij).EQ.1) THEN
          CALL subscl(ij)
!
          IF(irex.EQ.2) CALL hrex(ij,hrexg,hrexs)
!
          at = 0.D0
          bt = 0.D0
! ....
          at(1,1) = rgp(ij)
          bt(1)   = siegn(ij) * rgpn(ij) + rhg(ij) - hrexg
!    &             + dt * disg
!
          DO k=1, ncl
            k1=k+1
!
            CALL hvs(hv(k,ij), rlk(k,ij), rog(ij),
     &           ep(ij), ug(ij), ug(imj), uk(k,ij),
     &           uk(k,imj), vg(ij), vg(ijm), 
     &           vk(k,ij), vk(k,ijm), mug(ij),
     &           kapg(ij), cg(ij), k)
!
            at(1,1)   = at(1,1)       + dt * hv(k,ij) / cg(ij)
            at(1,k1)  =               - dt * hv(k,ij) / ck(k,ij)
            at(k1,1)  =               - dt * hv(k,ij) / cg(ij)
            at(k1,k1) = rlk(k,ij) + dt * hv(k,ij) / ck(k,ij)
!
            bt(k1) = rlkn(k,ij) * siekn(k,ij) + rhk(k, ij)
          END DO
!
          CALL invdm(at, bt, ij)
          sieg(ij) = bt(1)
          DO k=1, ncl
            siek(k,ij) = bt(k+1)
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
      INTEGER :: k
      REAL*8 :: div
!
      DO k=nphase,2,-1
!        IF(abs(a(k,k)).LE.1.D-6) THEN
        IF(rlk(k-1,ij)*inrl(k-1).LE.1.D-9) THEN
          a(1,k)=0.D0
          a(k,1)=0.D0
          b(k)=0.D0
        ELSE
!
! ... eliminate all cross elements in the gas enthalpy equation
!
          div=1.D0/a(k,k)
          a(k,1)=a(k,1)*div
          b(k)=b(k)*div
          b(1)=b(1)-a(1,k)*b(k)
          a(1,1)=a(1,1)-a(1,k)*a(k,1)
        ENDIF
      END DO
!
      b(1)=b(1)/a(1,1)
      DO k=2,nphase
        b(k)=b(k)-a(k,1)*b(1)
      END DO
!
      RETURN
      END SUBROUTINE invdm
!----------------------------------------------------------------------
      END MODULE enthalpy_matrix
