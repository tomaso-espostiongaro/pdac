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
      USE gas_solid_velocity, ONLY: ug, wg
      USE grid, ONLY: dr, rb, dz, r, indr, inrb, indz, inr
      USE grid, ONLY: ncint, myijk, data_exchange
      USE set_indexes
      USE time_parameters, ONLY: dt
      USE indijk_module, ONLY: ip0_jp0_kp0_
!
      IMPLICIT NONE
!
      REAL*8 :: cs(7), cs0
      INTEGER :: i, j, imesh
      INTEGER :: ijk
      INTEGER :: ig
!
      CALL data_exchange(rgpgc)
!
      DO ijk = 1, ncint

       imesh = myijk( ip0_jp0_kp0_, ijk)

       IF( fl_l(ijk) == 1 ) THEN

        CALL subscr(ijk)

        j  = ( imesh - 1 ) / nr + 1
        i  = MOD( ( imesh - 1 ), nr) + 1
!
        DO ig=1,ngas
          cs(ig)=rgpgcn(ig,ijk)
        END DO
        cs0=0.D0
        IF(wg(ijk).GT.0.D0) THEN
          cs0=cs0-dt*indz(j)*wg(ijk)
        ELSE
          DO ig=1,ngas
            cs(ig)=cs(ig)-dt*indz(j)*wg(ijk)*rgpgc(ig,ijt)
          END DO
        ENDIF
        IF(wg(ijm).GT.0.D0) THEN
          DO ig=1,ngas
            cs(ig)=cs(ig)+dt*indz(j)*wg(ijm)*rgpgc(ig,ijb)
          END DO
        ELSE
          cs0=cs0+dt*indz(j)*wg(ijm)
        ENDIF
        IF(ug(ijk).GT.0.D0) THEN
          cs0=cs0-dt*inr(i)*indr(i)*rb(i)*ug(ijk)
        ELSE
          DO ig=1,ngas
            cs(ig)=cs(ig)-dt*inr(i)*indr(i)*rb(i)*ug(ijk)*rgpgc(ig,ijr)
          END DO
        ENDIF
        IF(ug(imj).GT.0.D0)THEN
          DO ig=1,ngas
            cs(ig)=cs(ig) + dt * inr(i)*indr(i)*rb(i-1)*ug(imj)*rgpgc(ig,ijl)
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
        DO ig=1,ngas
          ygc(ig,ijk)=cs(ig)/cs0
        END DO

       END IF

      END DO
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      END MODULE gas_components
!----------------------------------------------------------------------
