!----------------------------------------------------------------------
      MODULE gas_components
!----------------------------------------------------------------------
      IMPLICIT NONE
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE ygas
!
      USE grid, ONLY: fl_l
      USE dimensions
      USE eos_gas, ONLY: rgpgc, rgpgcn, ygc
      USE gas_solid_velocity, ONLY: ug, vg
      USE grid, ONLY: dr, rb, dz, r, indr, inrb, indz, inr
      USE grid, ONLY: nij_l, myij, data_exchange
      USE set_indexes
      USE time_parameters, ONLY: dt
!
      IMPLICIT NONE
!
      REAL*8 :: cs(7), cs0
      INTEGER :: i, j, ij_g
      INTEGER :: ij
      INTEGER :: kg
!
      CALL data_exchange(rgpgc)
!
      DO ij = 1, nij_l
       ij_g = myij(0, 0, ij)
       IF(fl_l(ij).EQ.1) THEN
        CALL subscl(ij)
        j  = ( ij_g - 1 ) / ndi + 1
        i  = MOD( ( ij_g - 1 ), ndi) + 1
!
        DO kg=1,ngas
          cs(kg)=rgpgcn(kg,ij)
        END DO
        cs0=0.D0
        IF(vg(ij).GT.0.D0) THEN
          cs0=cs0-dt*indz(j)*vg(ij)
        ELSE
          DO kg=1,ngas
            cs(kg)=cs(kg)-dt*indz(j)*vg(ij)*rgpgc(kg,ijt)
          END DO
        ENDIF
        IF(vg(ijm).GT.0.D0) THEN
          DO kg=1,ngas
            cs(kg)=cs(kg)+dt*indz(j)*vg(ijm)*rgpgc(kg,ijb)
          END DO
        ELSE
          cs0=cs0+dt*indz(j)*vg(ijm)
        ENDIF
        IF(ug(ij).GT.0.D0) THEN
          cs0=cs0-dt*inr(i)*indr(i)*rb(i)*ug(ij)
        ELSE
          DO kg=1,ngas
            cs(kg)=cs(kg)-dt*inr(i)*indr(i)*rb(i)*ug(ij)*rgpgc(kg,ijr)
          END DO
        ENDIF
        IF(ug(imj).GT.0.D0)THEN
          DO kg=1,ngas
            cs(kg)=cs(kg) + dt * inr(i)*indr(i)*rb(i-1)*ug(imj)*rgpgc(kg,ijl)
          END DO
        ELSE
          cs0=cs0+dt*inr(i)*indr(i)*rb(i-1)*ug(imj)
        ENDIF
        cs(1)=cs(1)/(1.D0-cs0)
        IF(cs(1).LT.0.D0) cs(1)=0.D0
        cs(2)=cs(2)/(1.D0-cs0)
        IF(cs(2).LT.0.D0) cs(2)=0.D0
        cs(3)=cs(3)/(1.D0-cs0)
        IF(cs(3).LT.0.D0) cs(3)=0.D0
        cs(4)=cs(4)/(1.D0-cs0)
        IF(cs(4).LT.0.D0) cs(4)=0.D0
        cs(5)=cs(5)/(1.D0-cs0)
        IF(cs(5).LT.0.D0) cs(5)=0.D0
        cs(6)=cs(6)/(1.D0-cs0)
        cs(7)=cs(7)/(1.D0-cs0)
        IF(cs(7).LT.0.D0) cs(7)=0.D0
        cs0=cs(1)+cs(2)+cs(3)+cs(4)+cs(5)+cs(6)+cs(7)
        DO kg=1,ngas
          ygc(kg,ij)=cs(kg)/cs0
        END DO
       END IF
      END DO
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      END MODULE gas_components
!----------------------------------------------------------------------
