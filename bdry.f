!-----------------------------------------------------------------------
      MODULE boundary_conditions
!-----------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------
      SUBROUTINE boundary
!-----------------------------------------------------------------------
! ... This routine computes boundary conditions 
!
      USE atmosphere, ONLY: gravx, gravz, atm
      USE dimensions
      USE eos_gas, ONLY: rgpgc, xgc
      USE gas_constants, ONLY: gmw, gammaair, gamn, rgas
      USE gas_solid_density, ONLY: rgp, rlk
      USE gas_solid_temperature, ONLY: sieg, tg, sies, ts
      USE gas_solid_velocity, ONLY: ug, wg, us, ws
      USE grid, ONLY: zb, dz, rb, dr, r, inr, inrb
      USE grid, ONLY: ncint, myijk
      USE grid, ONLY: fl_l
      USE parallel, ONLY: mpime
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: p, ep
      USE set_indexes
      USE time_parameters, ONLY: dt, time
      USE indijk_module, ONLY: ip0_jp0_kp0_
!
      IMPLICIT NONE
!
      REAL*8 :: prif, pnn2, p1nn
      REAL*8 :: zrif, trif, rhorif, cost, costc
      REAL*8 :: rm1n, rmcn, rmcnn, rmc1n, rm0n, rm1nn, rm1knn, rm2n
      REAL*8 :: w2n, wg2, wcnn, wcn, w1n, wc1n
      REAL*8 :: u2nn, u2, w2nn, w2
      REAL*8 :: ucn, u1n, uc1n, ug2, u2n, ucnn
      REAL*8 :: eps, epc, epcn, ep1nn, epnn
      REAL*8 :: t1nn, tc1n, tcn, t1n, t0n, t2n, tgnn
      REAL*8 :: dz1, dzc, dr1, drc, indrc
      REAL*8 :: mg

      INTEGER :: ij, i, j, imesh, ig, is
      INTEGER :: n2
!
      DO ij = 1, ncint

        imesh = myijk( ip0_jp0_kp0_, ij)

        j = ( imesh - 1 ) / nr + 1
        i = MOD( ( imesh - 1 ), nr) + 1
        IF(fl_l(ij) == 1) THEN
          CALL subscr(ij)
!
! ***** right boundary conditions
!
         n2 = ipj

         IF ( (fl_l( ipj ) /= 1) .AND. (fl_l( ipjp ) /= 1) ) THEN 
          SELECT CASE (fl_l( n2 )) 

          CASE (2) 

            wg(n2)=wg(ij)
            DO is=1,nsolid
              ws(is,n2)=ws(is,ij)
            END DO

          CASE (3)

            wg(n2)=-wg(ij)
            DO is=1,nsolid
              IF(rlk(ij,is).GT.0.D0) ws(is,n2)=-ws(is,ij)
            END DO  

          CASE (4)
!
! ... definitions
            drc=(dr(nr)+dr(nr-1))*0.5D0
            indrc = 1.0D0/drc
            dr1=dr(nr-1)
            uc1n=ug(imj)
            ucn=ug(ij)
!
! ... interpolation of mid-point values
            u1n=(ug(imj)+ug(ij))*0.5D0
            u2n=(ug(n2)+ug(ij))*0.5D0
            epcn=(ep(ij)+ep(n2))*0.5D0
            tc1n=(tg(imj)+tg(ij))*0.5D0
            tcn=(tg(ij)+tg(n2))*0.5D0
            t0n=tg(imj)
            t1n=tg(ij)
            t2n=tg(n2)
!
            IF(ucn > 0.D0) THEN
! ... OUTFLOW ...
!
! ... extrapolations
!              ucnn=u1n*r(i)*inrb(i)
!              ug(n2)=u2n*r(i+1)*inrb(i+1)

! ... transport velocity field (non-reflecting b.c.)
              ug(n2) = ug(n2) - ucn*dt/dr(nr) * (ug(n2) - ug(ij))
              ucnn = ucn * ( 1 - dt/dr(nr-1) * (ug(ij) - ug(imj)) )
!
              IF(j.EQ.2) ug(ipjm)=-ug(n2)
              t1nn=tc1n
!
! ... calculation of the gas volumetric fraction at time (n+1)dt
! ... from the mass balance equation of solids in cell (ij)
              epc=0.D0
              DO is=1,nsolid
                eps = rlk(ij,is)*inrl(is) - dt*inrl(is)*inr(i)/dr1*          &
                      (rb(i)*us(is,ij)*rlk(ij,is) - rb(i-1)*us(is,imj)*rlk(imj,is))
                epc=epc+eps
              END DO
              ep1nn=1.D0-epc
!
! ... calculation of the mixture density from a mass balance equation for the mixture.
! ... Use gas velocity for mixture.
              rm2n = 0.D0
              rm1n = 0.D0
              rm0n = 0.D0
              DO is=1,nsolid
                rm2n=rm2n+rlk(n2,is)
                rm1n=rm1n+rlk(ij,is)
                rm0n=rm0n+rlk(imj,is)
              END DO
              rm1knn = rm1n
              rm2n = rm2n + rgp(n2)
              rm1n = rm1n + rgp(ij)
              rm0n = rm0n + rgp(imj)
              rmcn=(rm2n+rm1n)*0.5D0
              rmc1n=(rm1n+rm0n)*0.5D0
!
              rmcnn = rmcn - dt*indrc*inrb(i) * ( r(i+1)*u2n*rmcn - r(i)*u1n*rmc1n )
!
! ... Calculation of fluid pressure from the gas equation of state and the mixture density
! ... transport equation 
              mg=0.D0
              DO ig=1,ngas
                mg = mg + xgc(ig,ij) * gmw(ig)
              END DO
              rm1nn=ep1nn*mg/(rgas*t1nn)
!
              p1nn = (1.D0/rm1nn) *                           &
                     ( - rm1knn + rm1n - dt/dr1*inr(i) *      &
                     ( rb(i)*rm1n*ucn - rb(i-1)*rm0n*uc1n ) )
              IF (epc .LT. 1.0D-8) p1nn = p(ij)
!
! ... Calculation of the advanced-time fluid pressure from Momentum balance
! ... equation of the mixture
              p(n2) = -rmcnn*ucnn + rmcn*ucn - dt*indrc*inrb(i) *     &
                       ( r(i+1)*u2n*ucn*rmcn - r(i)*u1n*uc1n*rmc1n )
              p(n2) = p(n2)/(dt*indrc) + p1nn
!
              ep(n2) = ep(ij)
              tg(n2) = tg(ij)
              DO ig=1,ngas
                rgpgc(ig,n2) = rgpgc(ig,n2) - ucn * dt/dr(nr) * (rgpgc(ig,n2) - rgpgc(ig,ij))
              END DO
!
! ... Correct non-physical pressure
              IF (p(n2) <= 0.0D0) p(n2) = p(ij)
!
            ELSE IF (ucn < 0.D0) THEN
! ... INFLOW ...
!
! ... extrapolations
!              ug(n2)=ug(ij)*rb(i)*inrb(i+1)

! ... transport velocity field (non-reflecting b.c.)
              ug(n2) = ug(n2) - ucn*dt/dr(nr) * (0.D0 - u2n)
              zrif=zb(j)+0.5D0*(dz(1)-dz(j))
              CALL atm(zrif,prif,trif)
              rhorif=prif*gmw(6)/(rgas*trif)
              cost=prif/(rhorif**gammaair)
              costc=(gammaair*cost**(1.D0/gammaair))/(gammaair-1.D0)
!
! ... Adiabatic inflow
              p(n2)=(prif**gamn-(ug(n2)**2)/(2.D0*costc))**(1.D0/gamn)
!
              ep(n2) = 1.D0
              tg(n2) = trif
              DO ig=1,ngas
                rgpgc(ig,n2) = rgpgc(ig,ij)
              END DO

! ... Correct non-physical pressure
              IF (p(n2).LT.0.0D0) p(n2) = prif
!
            ELSE IF (ucn == 0.D0) THEN
!
              ug(n2) = ug(ij)
              p(n2)  = p(ij)
              ep(n2) = ep(ij)
              tg(n2) = tg(ij)
              DO ig=1,ngas
                rgpgc(ig,n2) = rgpgc(ig,ij)
              END DO

            ENDIF
!
! ... Set primary variables
!
            IF( fl_l (ipjp) == 3 ) ug(ipjp) = -ug(n2)
            IF( fl_l (ipjp) == 4 ) wg(n2) = wg(ij)
            DO is=1,nsolid
              IF( fl_l (ipjp) == 4 ) ws(is,n2)=ws(is,ij)
              IF(us(is,ij).GE.0.D0) THEN
                rlk(n2,is)=rlk(ij,is)
                us(is,n2)=ug(n2)
              ELSE
                rlk(n2,is)=0.0D0
                us(is,n2) =0.0D0
              ENDIF
            END DO

            sieg(n2) = sieg(ij)
            DO is=1,nsolid
              sies(is,n2)=sies(is,ij)
            END DO
!
          CASE DEFAULT
            CONTINUE
          END SELECT
         END IF
!
! ****** left boundary conditions
!
         n2 = imj

         IF ( (fl_l( imj ) /= 1) .AND. (fl_l( imjp ) /= 1) ) THEN
          SELECT CASE ( fl_l( n2 ) )
          CASE (2)
            wg(n2)=wg(ij)
            DO is=1,nsolid
              ws(is,n2)=ws(is,ij)
            END DO 
          CASE (3)
            wg(n2)=-wg(ij)
            DO is=1,nsolid
              IF(rlk(ij,is).GT.0.D0) ws(is,n2)=-ws(is,ij)
            END DO 
          CASE (4)
!
! ... definitions
            dr1 = dr(i)
            drc=(dr(i)+dr(i-1))*0.5D0
            indrc = 1.0D0/drc
            uc1n=ug(ij)
            ucn=ug(n2)
!
! ... interpolation of mid-point values
            u1n=(uc1n+ucn)*0.5D0
            u2n = ucn - (uc1n-ucn)/dr1 * 0.5*dr(i-1)
            epcn=(ep(ij)+ep(n2))*0.5D0
            tcn=(tg(ij)+tg(n2))*0.5D0
            tc1n=(tg(ij)+tg(ipj))*0.5D0
            t1n=tg(ij)
            t0n=tg(ipj)
!
! ... OUTFLOW ...
            IF(uc1n < 0.D0) THEN
!
! ... extrapolations
            ucnn = u1n
            t1nn = tc1n

! ... transport velocity field (non-reflecting b.c.)
            ug(n2) = ug(n2) - uc1n*dt/dr1 * (ug(ij) - ug(n2))
!
            IF(j.EQ.2) ug(imjm) = - ug(n2)
!
! ... calculation of the gas volumetric fraction at time (n+1)dt
! ... from the mass balance equation of solids in cell (ij)
              epc=0.D0
              DO is=1,nsolid
               eps=rlk(ij,is)*inrl(is) -                                 &
                   dt*inrl(is)*inr(i)/dr1*(rb(i)*us(is,ij)*rlk(ipj,is) -  &
                   rb(i-1)*us(is,imj)*rlk(ij,is))
               epc=epc+eps
              END DO
              ep1nn=1.D0-epc
!
! ... Mixture Density at time (n)dt
              rm2n = 0.D0
              rm1n = 0.D0
              rm0n = 0.D0
              DO is=1,nsolid
                rm2n=rm2n+rlk(n2,is)
                rm1n=rm1n+rlk(ij,is)
                rm0n=rm0n+rlk(ipj,is)
              END DO
              rm1knn = rm1n
              rm2n=rm2n+rgp(n2)
              rm1n=rm1n+rgp(ij)
              rm0n=rm0n+rgp(ipj)
!
              rmcn=(rm2n+rm1n)*0.5D0
              rmc1n=(rm1n+rm0n)*0.5D0
!
              rmcnn=rmcn-dt*indrc*inrb(i-1)*(r(i)*u1n*rmc1n-r(i-1)*u2n*rmcn)
!
! ... Calculation of the fluid pressure from a mass balance equation for the mixture.
! ... Use gas velocity for mixture.
!
              mg=0.D0
              DO ig=1,ngas
                mg = mg + xgc(ig,ij) * gmw(ig)
              END DO
              rm1nn = ep1nn*mg/(rgas*t1nn)
              p1nn = (1.D0/rm1nn) * &
                     (-rm1knn+rm1n-dt/dr1*inr(i)*(rb(i)*uc1n*rm0n-rb(i-1)*rm1n*ucn))
              IF (epc .LT. 1.0D-8) p1nn = p(ij)
!
! ... Calculation of the advanced-time fluid pressure from Momentum balance
! ... equation of the mixture
              p(n2) = -rmcnn*ucnn + rmcn*ucn - dt*indrc*inrb(i-1) *     &
                       ( r(i)*u1n*uc1n*rmc1n - r(i-1)*u2n*ucn*rmcn) + &
                        dt*p1nn*indrc
              p(n2)=p(n2)/(dt*indrc)
!
! ... Correct non-physical pressure
              IF (p(n2) <= 0.0D0) p(n2) = p(ij)
!
              ep(n2) = ep(ij)
              tg(n2) = tg(ij)
              DO ig=1,ngas
                rgpgc(ig,n2) = rgpgc(ig,n2) - uc1n * dt/dr(1) * (rgpgc(ig,ij) - rgpgc(ig,n2))
              END DO
!
! ... extrapolation of the temperature and solid fraction to time (n+1)dt
              DO is=1,nsolid
                IF(fl_l( imjp ) == 4) ws(is,n2)=ws(is,ij)
                rlk(n2,is)=rlk(ij,is)
                us(is,n2)=us(is,ij)
              END DO
!
! ... INFLOW ...
            ELSE IF (uc1n > 0.D0) THEN
!
              ug(n2) = ug(n2) - uc1n*dt/dr1 * (u1n - 0.D0)
!
              zrif=zb(j)+0.5D0*(dz(1)-dz(j))
              CALL atm(zrif,prif,trif)
              rhorif=prif*gmw(6)/(rgas*trif)
              cost=prif/(rhorif**gammaair)
              costc=(gammaair*cost**(1.D0/gammaair))/(gammaair-1.D0)
!
! ... Adiabatic inflow
              p(n2)=(prif**gamn-(ug(n2)**2)/(2.D0*costc))**(1.D0/gamn)
!
! ... Correct non-physical pressure
              IF (p(n2).LT.0.0D0) p(n2) = prif
!
              ep(n2) = 1.D0
              tg(n2) = trif
              DO ig=1,ngas
                rgpgc(ig,n2) = rgpgc(ig,ij)
              END DO

            ELSE IF (uc1n == 0.D0) THEN
!
              ug(n2) = ug(ij)
              p(n2)  = p(ij)
              ep(n2) = ep(ij)
              tg(n2) = tg(ij)
              DO ig=1,ngas
                rgpgc(ig,n2) = rgpgc(ig,ij)
              END DO

            ENDIF
!
! ... Set primary variables
!
            IF( fl_l ( imjp ) == 4 ) wg(n2) = wg(ij)
            IF( fl_l ( imjp ) == 3 ) ug(ipjp) = -ug(n2)
            DO is=1,nsolid
              IF( fl_l ( imjp ) == 4 ) ws(is,n2)=ws(is,ij)
              IF(us(is,ij).GE.0.D0) THEN
                rlk(n2,is)=rlk(ij,is)
                us(is,n2)=ug(n2)
              ELSE
                rlk(n2,is)=0.0D0
                us(is,n2) =0.0D0
              ENDIF
              sies(is,n2)=sies(is,ij)
              ts(is,n2)=ts(is,ij)
            END DO
            sieg(n2) = sieg(ij)

          CASE DEFAULT
            CONTINUE
          END SELECT
         END IF
!
! ****** top boundary conditions
!
         n2 = ijp

         IF ( (fl_l( ijp ) /= 1) .AND. (fl_l( ipjp ) /= 1) ) THEN
          SELECT CASE ( fl_l( n2 ) )

          CASE (2)

            ug(n2)=ug(ij)
            DO is=1,nsolid
              us(is,n2)=us(is,ij)
            END DO

            IF (j .EQ. (nz-1)) THEN
              IF(i .EQ. (nr-1)) THEN
                ug(ipjp) = ug(ipj)
              ELSE IF(i .EQ. 2) THEN
                ug(imjp) = ug(imj)
              ENDIF
            END IF

          CASE (3)

            ug(n2)=-ug(ij)
            DO is=1,nsolid
              IF(rlk(ij,is).GT.0.D0) us(is,n2)=-us(is,ij)
            END DO

          CASE (4)

! ... definitions and linear interpolations
          dzc=(dz(nz)+dz(nz-1))*0.5D0
          dz1=dz(nz-1)
!
            wc1n=wg(ijm)
            wcn=wg(ij)
            w1n=(wg(ijm)+wg(ij))*0.5D0
            w2n=(wg(n2)+wg(ij))*0.5D0
            epcn=(ep(ij)+ep(n2))*0.5D0
            tc1n=(tg(ijm)+tg(ij))*0.5D0
            tcn=(tg(ij)+tg(n2))*0.5D0
            t2n=tg(n2)
            t1n=tg(ij)
            t0n=tg(ij-nz)
!
! ... OUTFLOW ...
            IF(wcn > 0.D0) THEN
! ... extrapolations
!            wg(n2) = w2n
!            wcnn = w1n
            t1nn = tc1n

! ... transport velocity field (non-reflecting b.c.)
            wg(n2) = wg(n2) - wcn*dt/dz(nz) * (wg(n2) - wg(ij))
            wcnn   = wcn * ( 1 - dt/dz(nz-1) * (wg(ij) - wg(imj)) )
!
!
! ... calculation of gas volumetric fraction at time (n+1)dt
! ... from the mass balance equation of solids, in cell (ij)
              epc=0.D0
              DO is=1,nsolid
                eps=rlk(ij,is)*inrl(is) - dt*inrl(is)/dz1*(ws(is,ij)*rlk(ij,is) -  &
                    ws(is,ijm)*rlk(ijm,is))
                epc=epc+eps
              END DO 
              ep1nn=1.D0-epc
!
! ... Mixture density at time (n)dt
              rm2n=0.D0
              rm1n=0.D0
              rm0n=0.D0
              DO is=1,nsolid
                rm2n=rm2n+rlk(n2,is)
                rm1n=rm1n+rlk(ij,is)
                rm0n=rm0n+rlk(ijm,is)
              END DO 
              rm1knn=rm1n
              rm2n=rm2n+rgp(n2)
              rm1n=rm1n+rgp(ij)
              rm0n=rm0n+rgp(ijm)
              rmcn=(rm2n+rm1n)*0.5D0
              rmc1n=(rm1n+rm0n)*0.5D0
!
              rmcnn=rmcn-dt/dzc*(w2n*rmcn-w1n*rmc1n)
!
! ... Calculation of the fluid pressure from a mass balance equation for the mixture.
! ... Use gas velocity for mixture.
              mg=0.D0
              DO ig=1,ngas
                mg = mg + xgc(ig,ij) * gmw(ig)
              END DO
              rm1nn=ep1nn*mg/(rgas*t1nn)
              p1nn=(1.D0/rm1nn)*(-rm1knn+rm1n-dt/dz1*(rm1n*wcn-rm0n*wc1n))
              IF (epc .LT. 1.0D-8) p1nn = p(ij)
!
! ... Calculation of the advanced-time fluid pressure from Momentum balance
! ... equation of the mixture
! ... 
              pnn2=p(n2)
!
              p(n2)=-rmcnn*wcnn + rmcn*wcn-dt/dzc*(w2n*wcn*rmcn-w1n*wc1n*rmc1n) + &
                     dt*p1nn/dzc+gravz*dt*rmcn
              p(n2)=p(n2)/(dt/dzc)
!
! ... (Correct non-physical pressure)
              IF(p(n2).LE.0.D0) p(n2) = pnn2
!
              ep(n2) = ep(ij)
              tg(n2) = tg(ij)
              DO ig=1,ngas
                rgpgc(ig,n2) = rgpgc(ig,n2) - wcn*dt/dz(nz) * (rgpgc(ig,n2) - rgpgc(ig,ij))
              END DO
             
            ELSEIF(wcn < 0.D0) THEN
! ... INFLOW ...
!
              wg2=wg(ij)
              wg(n2) = wg(n2) - wcn*dt/dz(nz) * (0.D0 - w2n)
!
              zrif=zb(nz)+0.5D0*(dz(1)-dz(nz))
              CALL atm(zrif,prif,trif)
              rhorif=prif*gmw(6)/(rgas*trif)
              cost=prif/(rhorif**gammaair)
              costc=(gammaair*cost**(1.D0/gammaair))/(gammaair-1.D0)
              pnn2=p(n2)
              p(n2)=(prif**gamn-(wg(n2)**2/2.D0)/costc)**(1.D0/gamn)
!
! ... (Correct non-physical pressure)
              IF(p(n2).LE.0.D0) p(n2)=pnn2
!
              ep(n2)=1.D0
              tg(n2)=trif
              DO ig=1,ngas
                rgpgc(ig,n2) = rgpgc(ig,ij)
              END DO

            ELSEIF(wcn == 0.D0) THEN
              wg(n2) = wg(ij) 
              p(n2)  = p(ij)
              ep(n2) = ep(ij)
              tg(n2) = tg(ij)
              DO ig=1,ngas
                rgpgc(ig,n2) = rgpgc(ig,ij)
              END DO
            ENDIF
!
! ... set upper corners velocities
!
            IF (j .EQ. (nz-1)) THEN
              IF(i .EQ. (nr-1)) THEN
                wg(ipjp) = wg(n2)
                ug(ipjp) = ug(n2)
              ELSE IF(i .EQ. 2) THEN
                wg(imjp) = wg(n2)
                ug(imjp) = ug(n2)
              ENDIF
            END IF
!
            IF( fl_l( ipjp ) == 4 ) ug(n2) = ug(ij)
            sieg(n2) = sieg(ij)
!
            DO is=1,nsolid
              IF( fl_l( ipjp ) == 4 ) us(is,n2)=us(is,ij)
              IF(ws(is,ij).GE.0.D0) THEN
                rlk(n2,is)=rlk(ij,is)
                ws(is,n2)=ws(is,ij)
              ELSE
                rlk(n2,is)=0.D0
                ws(is,n2)=0.D0
              ENDIF
              IF(i .EQ. (nr-1) .AND. j .EQ. (nz-1)) THEN
                ws(is,ipjp)=ws(is,n2)
                us(is,ipjp)=us(is,n2)
              ENDIF
            END DO 

            DO is=1,nsolid
              sies(is,n2)=sies(is,ij)
              ts(is,n2)=ts(is,ij)
              IF (i .EQ. (nr-1) .AND. j .EQ. (nz-1)) THEN
                ts(is,ipjp)=ts(is,n2)
              END IF
            END DO
         
          CASE DEFAULT
            CONTINUE
          END SELECT
         END IF
!
! ***** bottom boundary conditions
!
         n2 = ijm

         IF ( (fl_l( ijm ) /= 1) .AND. (fl_l( ipjm ) /= 1) ) THEN
          SELECT CASE ( fl_l( n2 ) )

          CASE (2)

            ug(n2)=ug(ij)
            DO is=1,nsolid
              us(is,n2)=us(is,ij)
            END DO 
            IF (j .EQ. (2)) THEN
              IF(i .EQ. (nr-1)) THEN
                ug(ipjm) = ug(ipj)
              ELSE IF(i .EQ. 2) THEN
                ug(imjm) = ug(imj)
              ENDIF
            END IF

          CASE (3)

            ug(n2)=-ug(ij)
            DO is=1,nsolid
              IF(rlk(ij,is).GT.0.D0) us(is,n2)=-us(is,ij)
            END DO 

          CASE DEFAULT

            CONTINUE

          END SELECT
         END IF
!
        END IF
      END DO
!
      RETURN
!-----------------------------------------------------------------------
      END SUBROUTINE boundary
!-----------------------------------------------------------------------
      SUBROUTINE boundary3d
!-----------------------------------------------------------------------
! ... This routine computes (x,y,z) boundary conditions 
!
      USE atmosphere, ONLY: gravx, gravz, atm
      USE dimensions
      USE eos_gas, ONLY: rgpgc, xgc
      USE gas_constants, ONLY: gmw, gammaair, gamn, rgas
      USE gas_solid_density, ONLY: rgp, rlk
      USE gas_solid_temperature, ONLY: sieg, tg, sies, ts
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      USE grid, ONLY: zb, dz, dx, dy
      USE grid, ONLY: ncint, myijk
      USE grid, ONLY: fl_l
      USE parallel, ONLY: mpime
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: p, ep
      USE set_indexes, ONLY: subscr, ipjk,  imjk,  ippjk,  immjk,  ijpk,  &
        ipjpk,  imjpk,  ijmk,  ipjmk,  imjmk,  ijppk,  ijmmk,  ijkp,  ipjkp,  &
        imjkp,  ijpkp,  ijmkp,  ijkm,  ipjkm,  imjkm,  ijpkm,  ijmkm,  ijkpp,  ijkmm

      USE time_parameters, ONLY: dt, time
      USE indijk_module, ONLY: ip0_jp0_kp0_
!
      IMPLICIT NONE
!
      REAL*8 :: prif, pnn2, p1nn
      REAL*8 :: zrif, trif, rhorif, cost, costc
      REAL*8 :: rm1n, rmcn, rmcnn, rmc1n, rm0n, rm1nn, rm1knn, rm2n
      REAL*8 :: w2n, wg2, wcnn, wcn, w1n, wc1n
      REAL*8 :: u2nn, u2, w2nn, w2
      REAL*8 :: ucn, u1n, uc1n, ug2, u2n, ucnn
      REAL*8 :: eps, epc, epcn, ep1nn, epnn
      REAL*8 :: t1nn, tc1n, tcn, t1n, t0n, t2n, tgnn
      REAL*8 :: dz1, dzc, dx1, dxc, dy1, dyc
      REAL*8 :: mg

      INTEGER :: ijk, i, j, k, imesh, ig, is
      INTEGER :: n2
!
      DO ijk = 1, ncint

        imesh = myijk( ip0_jp0_kp0_, ijk)

        i = MOD( MOD( ijk - 1, nx*ny ), nx ) + 1
        j = MOD( ijk - 1, nx*ny ) / nx + 1
        k = ( ijk - 1 ) / ( nx*ny ) + 1


        IF( fl_l(ijk) == 1 ) THEN

            CALL subscr(ijk)
!
! ***** East boundary conditions ***** !
!
            n2 = ipjk

            SELECT CASE ( fl_l( n2 ) ) 

            CASE (2) 

              vg( n2 ) = vg( ijk )
              wg( n2 ) = wg( ijk )
              DO is = 1, nsolid
                vs( is, n2 ) = vs( is, ijk )
                ws( is, n2 ) = ws( is, ijk )
              END DO

            CASE (3)

              vg( n2 ) = -vg( ijk )
              wg( n2 ) = -wg( ijk )
              DO is = 1, nsolid
                vs( is, n2 ) = -vs( is, ijk )
                ws( is, n2 ) = -ws( is, ijk )
              END DO  

            CASE (4)
!
! ... definitions

              dxc = ( dx( nx ) + dx( nx-1 ) ) * 0.5D0
              dx1 = dx( nx-1 )
              
              uc1n = ug( imjk )
              ucn  = ug( ijk  )
!
! ... interpolation of mid-point values

              u1n = ( ug( imjk ) + ug( ijk ) ) * 0.5D0
              u2n = ( ug( n2   ) + ug( ijk ) ) * 0.5D0
              epcn = ( ep( ijk ) + ep( n2 ) ) * 0.5D0
              tc1n = ( tg( imjk ) + tg( ijk ) ) * 0.5D0
              tcn  = ( tg( ijk ) + tg( n2 ) ) * 0.5D0
              t0n  = tg( imjk )
              t1n  = tg( ijk )
              t2n  = tg( n2 )

!
              IF( ucn > 0.D0 ) THEN

! ...  OUTFLOW ...
!
! ...  Extrapolations
!               ucnn = u1n
!               ug(n2) = u2n
!
! ...  Non-reflecting boundary conditions

                ug(n2) = ug(n2) - ucn*dt/dx(nx) * (ug(n2) - ug(ijk))
                ucnn = ucn * ( 1 - dt/dx(nx-1) * (ug(ijk) - ug(imjk)) )

! MODIFICARE X3D ....  IF(j == .EQ.2) ug(ipjm)=-ug(n2) ! 

                t1nn = tc1n
!
! ... calculation of the gas volumetric fraction at time (n+1)dt
! ... from the mass balance equation of solids in cell (ijk)

                epc = 0.D0
                DO is = 1, nsolid
                  eps = rlk(ijk,is) * inrl(is) - ( dt * inrl(is) / dx1 )  *          &
                        ( us(is,ijk)*rlk(ijk,is) - us(is,imjk)*rlk(imjk,is) )
                  epc = epc + eps
                END DO
                ep1nn = 1.D0 - epc
!
! ... calculation of the mixture density 
! ... from a mass balance equation for the mixture.
! ... Use gas velocity for mixture.

                rm2n=0.D0
                rm1n=0.D0
                rm0n=0.D0
                DO is=1,nsolid
                  rm2n=rm2n+rlk(n2,is)
                  rm1n=rm1n+rlk(ijk,is)
                  rm0n=rm0n+rlk(imjk,is)
                END DO
                rm1knn = rm1n
                rm2n = rm2n + rgp(n2)
                rm1n = rm1n + rgp(ijk)
                rm0n = rm0n + rgp(imjk)
                rmcn=(rm2n+rm1n)*0.5D0
                rmc1n=(rm1n+rm0n)*0.5D0
!
                rmcnn = rmcn-dt/dxc*(u2n*rmcn-u1n*rmc1n)
!
! ... Calculation of fluid pressure from the gas equation of state and 
! ... the mixture density transport equation
!
                mg=0.D0
                DO ig=1,ngas
                  mg = mg + xgc(ig,ijk) * gmw(ig)
                END DO
                rm1nn=ep1nn*mg/(rgas*t1nn)
                
                p1nn = (1.D0/rm1nn) * ( - rm1knn + rm1n - dt/dx1 * ( rm1n*ucn - rm0n*uc1n ) )
                IF (epc < 1.0D-8) p1nn = p(ijk)
!
! ... Calculation of the advanced-time fluid pressure from Momentum balance
! ... equation of the mixture

                p(n2) = -rmcnn*ucnn + rmcn*ucn - dt/dxc * ( u2n*ucn*rmcn - u1n*uc1n*rmc1n) 
                p(n2) = p(n2)/(dt/dxc) + p1nn
!
                ep(n2) = ep(ijk)
                tg(n2) = tg(ijk)
                DO ig=1,ngas
                   rgpgc(ig,n2) = rgpgc(ig,n2) - &
                                  ucn * dt/dx(nx) * (rgpgc(ig,n2) - rgpgc(ig,ijk))
                END DO
!
! ... Correct non-physical pressure
                IF (p(n2) <= 0.0D0) p(n2) = p(ijk)
!
              ELSE IF ( ucn <  0.D0 ) THEN

! ... INFLOW ...
!
! ... Extrapolation

!                ug(n2) = ug(ijk)

! ... Non-reflecting b.c.
                
                ug(n2) = ug(n2) - ucn*dt/dx(nx) * (0.D0 - u2n)
!
                zrif=zb(k)+0.5D0*(dz(1)-dz(k))
                CALL atm(zrif,prif,trif)
                rhorif=prif*gmw(6)/(rgas*trif)
                cost=prif/(rhorif**gammaair)
                costc=(gammaair*cost**(1.D0/gammaair))/(gammaair-1.D0)
!
! ... Adiabatic inflow
 
                p(n2)=(prif**gamn-(ug2**2)/(2.D0*costc))**(1.D0/gamn)
!
                ep(n2) = 1.D0
                tg(n2) = trif
                DO ig=1,ngas
                  rgpgc(ig,n2) = rgpgc(ig,ijk)
                END DO

! ... Correct non-physical pressure

                IF ( p(n2) < 0.0D0 ) p(n2) = prif
!
              ELSE IF ( ug(ijk) == 0.D0 ) THEN
!
                ug(n2)   = ug(ijk)
                p(n2)    = p(ijk)
                ep(n2)   = ep(ijk)
                tg(n2)   = tg(ijk)
                DO ig=1,ngas
                  rgpgc(ig,n2) = rgpgc(ig,ijk)
                END DO

              ENDIF
!
! ... Set primary variables
!
              ! ... IF(nfltr.EQ.3) ug(ipjp) = -ug(n2)
              ! ... IF(nfltr.EQ.4) wg(n2) = wg(ijk)
              DO is = 1, nsolid
                ! .. IF(nfltr.EQ.4) ws(is,n2)=ws(is,ijk)
                IF( us(is,ijk) >= 0.D0 ) THEN
                  rlk(n2,is) = rlk(ijk,is)
                  us(is,n2)  = us(is,ijk)
                ELSE
                  rlk(n2,is) = 0.0D0
                  us(is,n2)  = 0.0D0
                ENDIF
              END DO

              sieg(n2) = sieg(ijk)
              DO is = 1, nsolid
                sies(is,n2) = sies(is,ijk)
              END DO
!              
            CASE DEFAULT

              CONTINUE

            END SELECT

!
! ***** West boundary conditions ***** !
!

            n2 = imjk

            SELECT CASE ( fl_l( n2 ) )

            CASE (2)

              vg(n2)=vg(ijk)
              wg(n2)=wg(ijk)
              DO is=1,nsolid
                vs(is,n2)=vs(is,ijk)
                ws(is,n2)=ws(is,ijk)
              END DO 

            CASE (3)

              vg(n2)=-vg(ijk)
              wg(n2)=-wg(ijk)
              DO is=1,nsolid
                vs(is,n2)=-vs(is,ijk)
                ws(is,n2)=-ws(is,ijk)
              END DO 

            CASE (4)
!
! ... definitions

              dx1 = dx(i)
              dxc = (dx(i)+dx(i-1))*0.5D0
              uc1n = ug(ijk)
              ucn  = ug(n2)
!
! ... interpolation of mid-point values

              u1n  = (uc1n+ucn)*0.5D0
              epcn = (ep(ijk)+ep(n2))*0.5D0
              tcn  = (tg(ijk)+tg(n2))*0.5D0
              tc1n = (tg(ijk)+tg(ipjk))*0.5D0
              t1n  = tg(ijk)
              t0n  = tg(ipjk)
!
              u2n = ucn - (uc1n-ucn)/dx1 * 0.5*dx(i-1)
!
! ... OUTFLOW ...

              IF( uc1n < 0.D0 ) THEN
!
! ... Extrapolations
                ucnn = u1n
                t1nn = tc1n
! 
! ... Non-reflecting boundary condition
!
                ug(n2) = ug(n2) - uc1n*dt/dx1 * (ug(ijk) - ug(n2))
!
! ... MODIFICAREX3D IF(j.EQ.2) ug(imjmk) = - ug(n2)
!
! ... calculation of the gas volumetric fraction at time (n+1)dt
! ... from the mass balance equation of solids in cell (ijk)

                epc = 0.D0
                DO is = 1, nsolid
                 eps = rlk(ijk,is) * inrl(is) -                                 &
                     dt*inrl(is)/dx1*(us(is,ijk)*rlk(ipjk,is) - us(is,imjk)*rlk(ijk,is))
                 epc=epc+eps
                END DO
                ep1nn=1.D0-epc
!
! ... Mixture Density at time (n)dt
                rm2n=0.D0
                rm1n=0.D0
                rm0n=0.D0
                DO is=1,nsolid
                  rm2n=rm2n+rlk(n2,is)
                  rm1n=rm1n+rlk(ijk,is)
                  rm0n=rm0n+rlk(ipjk,is)
                END DO
                rm1knn = rm1n
                rm2n = rm2n + rgp(n2)
                rm1n = rm1n + rgp(ijk)
                rm0n = rm0n + rgp(ipjk)
!
                rmcn=(rm2n+rm1n)*0.5D0
                rmc1n=(rm1n+rm0n)*0.5D0
!
                rmcnn = rmcn - dt/dxc * (u1n*rmc1n-u2n*rmcn)
!
! ... calculation of fluid pressure 
! ... from a mass balance equation for the mixture.
! ... Use gas velocity for mixture.
!
                mg=0.D0
                DO ig=1,ngas
                  mg = mg + xgc(ig,ijk) * gmw(ig)
                END DO
                rm1nn = ep1nn*mg/(rgas*t1nn)
                p1nn = (1.D0/rm1nn) * &
                       (-rm1knn+rm1n-dt/dx1*(uc1n*rm0n-rm1n*ucn))
                IF (epc < 1.0D-8) p1nn = p(ijk)
!
! ... Calculation of the advanced-time fluid pressure from Momentum balance
! ... equation of the mixture

                p(n2) = -rmcnn*ucnn + rmcn*ucn - dt/dxc *     &
                         ( u1n*uc1n*rmc1n - u2n*ucn*rmcn) + dt/dxc * p1nn
                p(n2)=p(n2)/(dt/dxc)
!
! ... Correct non-physical pressure
                IF (p(n2) <= 0.0D0) p(n2) = p(ijk)
                ep(n2) = ep(ijk)
                tg(n2) = tg(ijk)
                DO ig=1,ngas
                  rgpgc(ig,n2) = rgpgc(ig,n2) - &
                                 uc1n * dt/dx(1) * (rgpgc(ig,ijk) - rgpgc(ig,n2))
                END DO
!
! ... extrapolation of the temperature and solid fraction to time (n+1)dt

                DO is=1,nsolid
                  ! IF(nfllt.EQ.4) ws(is,n2)=ws(is,ijk)
                  rlk(n2,is)=rlk(ijk,is)
                  us(is,n2)=us(is,ijk)
                END DO

              ELSE IF( uc1n > 0.D0 ) THEN

! ... INFLOW ...
!
                ug(n2) = ug(n2) - uc1n*dt/dx1 * (u1n - 0.D0)

                zrif = zb(k)+0.5D0*(dz(1)-dz(k))
                CALL atm(zrif,prif,trif)
                rhorif = prif*gmw(6)/(rgas*trif)
                cost = prif/(rhorif**gammaair)
                costc = (gammaair*cost**(1.D0/gammaair))/(gammaair-1.D0)
!
! ... Adiabatic inflow

                p(n2)=(prif**gamn-(ug2**2)/(2.D0*costc))**(1.D0/gamn)
!
! ... Correct non-physical pressure

                IF ( p(n2) < 0.0D0 ) p(n2) = prif
!
                ep(n2) = 1.D0
                tg(n2) = trif
                DO ig=1,ngas
                  rgpgc(ig,n2) = rgpgc(ig,ijk)
                END DO

              ENDIF
!
! .... Set primary variables
!
              ! ... IF(nfllt.EQ.4) wg(n2)=wg(ijk)
              ! ... IF(nfllt.EQ.3) ug(ipjp) = -ug(n2)
              DO is=1,nsolid
                ! .. IF(nfllt.EQ.4) ws(is,n2)=ws(is,ijk)
                IF( us(is,ijk) >= 0.D0 ) THEN
                  rlk(n2,is)=rlk(ijk,is)
                  us(is,n2)=us(is,ijk)
                ELSE
                  rlk(n2,is)=0.0D0
                  us(is,n2) =0.0D0
                ENDIF
              END DO
!
              sieg(n2) = sieg(ijk)
              DO is=1,nsolid
                sies(is,n2)=sies(is,ijk)
              END DO
!
            CASE DEFAULT

              CONTINUE

            END SELECT
!
! ***** north boundary conditions ***** !
!
!
            n2 = ijpk

            SELECT CASE (  fl_l( n2 ) )

            CASE (2)

              ug( n2 ) = ug( ijk )
              wg( n2 ) = wg( ijk )
              DO is = 1, nsolid
                us(is,n2) = us(is,ijk)
                ws(is,n2) = ws(is,ijk)
              END DO 

            CASE (3)

              ug(n2) = -ug(ijk)
              wg(n2) = -wg(ijk)
              DO is = 1, nsolid
                 us(is,n2) = -ws(is,ijk)
                 ws(is,n2) = -ws(is,ijk)
              END DO 

            CASE DEFAULT

              CONTINUE

            END SELECT
!
! ***** south boundary conditions ***** !
!
!
            n2 = ijmk

            SELECT CASE (  fl_l( n2 ) )

            CASE (2)

              ug( n2 ) = ug( ijk )
              wg( n2 ) = wg( ijk )
              DO is = 1, nsolid
                us(is,n2) = us(is,ijk)
                ws(is,n2) = ws(is,ijk)
              END DO 

            CASE (3)

              ug(n2) = -ug(ijk)
              wg(n2) = -wg(ijk)
              DO is = 1, nsolid
                 us(is,n2) = -us(is,ijk)
                 ws(is,n2) = -ws(is,ijk)
              END DO 

            CASE DEFAULT

              CONTINUE

            END SELECT
!
! ***** top boundary conditions ***** !
!

            n2 = ijkp

            SELECT CASE ( fl_l( n2 ) )

            CASE (2)

              ug(n2) = ug(ijk)
              vg(n2) = vg(ijk)
              DO is = 1, nsolid
                us(is,n2) = us(is,ijk)
                vs(is,n2) = vs(is,ijk)
              END DO

              IF ( k == (nz-1) ) THEN
                IF( ( i == (nx-1) ) .AND. ( j == (ny-1) ) ) THEN
                  ! ug(ipjp) = ug(ipj)
                ELSE IF( ( i == 2 ) .AND. ( j == (ny-1) ) ) THEN
                ELSE IF( ( i == (nx-1) ) .AND. ( j == 2 ) ) THEN
                ELSE IF( ( i == 2 ) .AND. ( j == 2 ) ) THEN
                  ! ug(imjp) = ug(imj)
                ENDIF
              END IF

            CASE (3)
  
              ug(n2) = -ug(ijk)
              vg(n2) = -vg(ijk)
              DO is = 1, nsolid
                 us(is,n2) = -us(is,ijk)
                 vs(is,n2) = -vs(is,ijk)
              END DO

            CASE (4)
!
! ... definitions and linear interpolations

              dzc=(dz(nz)+dz(nz-1))*0.5D0
              dz1=dz(nz-1)
!
              wc1n = wg( ijkm )
              wcn  = wg( ijk )
              w1n  =( wg( ijkm ) + wg( ijk ) ) * 0.5D0
              w2n  =( wg( n2 ) + wg( ijk ) ) * 0.5D0
              epcn =( ep( ijk ) + ep( n2 ) ) * 0.5D0
              tc1n =( tg( ijkm ) + tg( ijk ) ) * 0.5D0
              tcn  =( tg( ijk ) + tg( n2 ) ) * 0.5D0
              t2n  = tg( n2 )
              t1n  = tg( ijk )
              t0n  = tg( ijkm )
!
! ... OUTFLOW ...

              IF( wcn > 0.D0) THEN
!
! ... Extrapolations
!                wg(n2) = w2n
!                wcnn = w1n
                t1nn = tc1n

! ... Non-reflecting boundary condition 
!
                wg(n2) = wg(n2) - wcn*dt/dz(nz) * (wg(n2) - wg(ijk))
                wcnn   = wcn * ( 1 - dt/dz(nz-1) * (wg(ijk) - wg(ijkm)) )
!
! ... calculation of gas volumetric fraction at time (n+1)dt
! ... from the mass balance equation of solids, in cell (ijk)

                epc = 0.D0
                DO is = 1, nsolid
                  eps = rlk(ijk,is)*inrl(is) - dt*inrl(is)/dz1*(ws(is,ijk)*rlk(ijk,is) -  &
                      ws(is,ijkm)*rlk(ijkm,is))
                  epc = epc + eps
                END DO 
                ep1nn = 1.D0 - epc
!
! ... calculation of mixture density at time (n+1)dt

                rm2n=0.D0
                rm1n=0.D0
                rm0n=0.D0
                DO is=1,nsolid
                  rm2n=rm2n+rlk(n2,is)
                  rm1n=rm1n+rlk(ijk,is)
                  rm0n=rm0n+rlk(ijkm,is)
                END DO 
                rm1knn = rm1n
                rm2n=rm2n+p(n2)*ep(n2)*gmw(6)/(rgas*tg(n2))
                rm1n=rm1n+p(ijk)*ep(ijk)*gmw(6)/(rgas*tg(ijk))
                rm0n=rm0n+p(ijkm)*ep(ijkm)*gmw(6)/(rgas*tg(ijkm))
                rmcn=(rm2n+rm1n)*0.5D0
                rmc1n=(rm1n+rm0n)*0.5D0
!
                rmcnn=rmcn-dt/dzc*(w2n*rmcn-w1n*rmc1n)
!
! ... Calculation of the fluid pressure from a mass balance equation 
! ... for the mixture.
                mg=0.D0
                DO ig=1,ngas
                  mg = mg + xgc(ig,ijk) * gmw(ig)
                END DO
                rm1nn=ep1nn*mg/(rgas*t1nn)

                p1nn=(1.D0/rm1nn)*(-rm1knn+rm1n-dt/dz1*(rm1n*wcn-rm0n*wc1n))
                IF (epc < 1.0D-8) p1nn = p(ijk)
!
! ... Calculation of the advanced-time fluid pressure from Momentum balance
! ... equation of the mixture
! 
                pnn2=p(n2)
!
                p(n2) = -rmcnn*wcnn + rmcn*wcn-dt/dzc*(w2n*wcn*rmcn-w1n*wc1n*rmc1n) + &
                       dt*p1nn/dzc+gravz*dt*rmcn
                p(n2) = p(n2)/(dt/dzc)

                ep(n2) = ep(ijk)
                tg(n2) = tg(ijk)
                DO ig=1,ngas
                  rgpgc(ig,n2) = rgpgc(ig,n2) - &
                                 wcn*dt/dz(nz) * (rgpgc(ig,n2) - rgpgc(ig,ijk))
                END DO
!
! ... (Correct non-physical pressure)

                IF( p(n2) <= 0.D0 ) p(n2) = pnn2
!
              ELSEIF( wcn < 0.D0 ) THEN

! ... INFLOW ...
!
                wg2=wg(ijk)

! ... Extrapolation
!                wg(n2) = wg2
!
! ... Non-reflecting b.c.
                wg(n2) = wg(n2) - wcn*dt/dz(nz) * (0.D0 - w2n)
!
                zrif=zb(nz)+0.5D0*(dz(1)-dz(nz))
                CALL atm(zrif,prif,trif)
                rhorif=prif*gmw(6)/(rgas*trif)
                cost=prif/(rhorif**gammaair)
                costc=(gammaair*cost**(1.D0/gammaair))/(gammaair-1.D0)
!
                ep(n2)=1.D0
                tg(n2)=trif
                DO ig=1,ngas
                  rgpgc(ig,n2) = rgpgc(ig,ijk)
                END DO

                pnn2=p(n2)
                p(n2)=(prif**gamn-(wg2**2/2.D0)/costc)**(1.D0/gamn)
!
! ... (Correct non-physical pressure)

                IF( p(n2) < 0.D0 ) p(n2)=pnn2
!
              ELSEIF( wg(ijk) == 0.D0 ) THEN

                wg(n2) = wg(ijk) 
                p(n2)  = p(ijk)
                ep(n2) = ep(ijk)
                tg(n2) = tg(ijk)
                DO ig=1,ngas
                  rgpgc(ig,n2) = rgpgc(ig,ijk)
                END DO

              ENDIF
!
! ... set upper corners velocities
!
              IF ( k == (nz-1) ) THEN
                IF( ( i == (nx-1) ) .AND. ( j == (ny-1) ) ) THEN
                  ! wg(ipjp) = wg(n2)
                  ! ug(ipjp) = ug(n2)
                ELSE IF( ( i == 2 ) .AND. ( j == (ny-1) ) ) THEN
                ELSE IF( ( i == (nx-1) ) .AND. ( j == 2 ) ) THEN
                ELSE IF( ( i == 2 ) .AND. ( j == 2 ) ) THEN
                  ! wg(imjp) = wg(n2)
                  ! ug(imjp) = ug(n2)
                ENDIF
              END IF

              !  IF(nfltr.EQ.4) ug(n2) = ug(ijk)
!
              DO is=1,nsolid
                !  IF(nfltr.EQ.4) us(is,n2)=us(is,ijk)
                IF( ws(is,ijk) >= 0.D0) THEN
                  rlk(n2,is)=rlk(ijk,is)
                  ws(is,n2)=ws(is,ijk)
                ELSE
                  rlk(n2,is)=0.D0
                  ws(is,n2)=0.D0
                ENDIF
                ! IF(i .EQ. (nr-1) .AND. j .EQ. (nz-1)) THEN
                !   ws(is,ipjp)=ws(is,n2)
                !   us(is,ipjp)=us(is,n2)
                ! ENDIF
              END DO 

              sieg(n2) = sieg(ijk)
              DO is=1,nsolid
                sies(is,n2)=sies(is,ijk)
                ts(is,n2)=ts(is,ijk)
                ! IF (i .EQ. (nr-1) .AND. j .EQ. (nz-1)) THEN
                !   ts(is,ipjp)=ts(is,n2)
                ! END IF
              END DO
!
            CASE DEFAULT

              CONTINUE

            END SELECT

!
! ***** bottom boundary conditions ***** !
!
!
            n2 = ijkm

            SELECT CASE (  fl_l( n2 ) )

            CASE (2)

              ug( n2 ) = ug( ijk )
              vg( n2 ) = vg( ijk )
              DO is = 1, nsolid
                us(is,n2) = us(is,ijk)
                vs(is,n2) = vs(is,ijk)
              END DO 
              IF (k == 2) THEN
                IF( ( i == (nx-1) ) .AND. ( j == (ny-1) ) ) THEN
                  ! ug(ipjm) = ug(ipj)
                ELSE IF( ( i == 2 ) .AND. ( j == (ny-1) ) ) THEN
                ELSE IF( ( i == (nx-1) ) .AND. ( j == 2 ) ) THEN
                ELSE IF( ( i == 2 ) .AND. ( j == 2 ) ) THEN
                  ! ug(imjm) = ug(imj)
                ENDIF
              END IF

            CASE (3)

              ug(n2) = -ug(ijk)
              vg(n2) = -vg(ijk)
              DO is = 1, nsolid
                 us(is,n2) = -us(is,ijk)
                 vs(is,n2) = -vs(is,ijk)
              END DO 

            CASE DEFAULT

              CONTINUE

            END SELECT
!
        END IF
      END DO
!
      RETURN
!-----------------------------------------------------------------------
      END SUBROUTINE boundary3d
!-----------------------------------------------------------------------
      END MODULE boundary_conditions
!-----------------------------------------------------------------------
