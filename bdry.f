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
      USE eos_gas, ONLY: rgpgc
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

      INTEGER :: ij, i, j, imesh, ig, is
      INTEGER :: n2, nflr, nflt, nfll, nflb
      INTEGER :: nflbr, nfltr, nfllt, nfllb
!
      prif=0.D0; pnn2=0.D0; p1nn=0.D0
      zrif=0.D0; trif=0.D0; rhorif=0.D0; cost=0.D0; costc=0.D0
      rm1n=0.D0; rmcn=0.D0; rmcnn=0.D0; rmc1n=0.D0; rm0n=0.D0; rm1nn=0.D0; rm1knn=0.D0; rm2n=0.D0
      w2n=0.D0; wg2=0.D0; wcnn=0.D0; wcn=0.D0; w1n=0.D0; wc1n=0.D0
      u2nn=0.D0; u2=0.D0; w2nn=0.D0; w2=0.D0
      ucn=0.D0; u1n=0.D0; uc1n=0.D0; ug2=0.D0; u2n=0.D0; ucnn=0.D0
      eps=0.D0; epc=0.D0; epcn=0.D0; ep1nn=0.D0; epnn=0.D0
      t1nn=0.D0; tc1n=0.D0; tcn=0.D0; t1n=0.D0; t0n=0.D0; t2n=0.D0; tgnn=0.D0
      dz1=0.D0; dzc=0.D0; dr1=0.D0; drc=0.D0; indrc=0.D0
!

      DO ij = 1, ncint

        imesh = myijk( ip0_jp0_kp0_, ij)

        j = ( imesh - 1 ) / nr + 1
        i = MOD( ( imesh - 1 ), nr) + 1
        IF(fl_l(ij).EQ.1) THEN
          CALL subscr(ij)
!
          nflr=fl_l(ipj)
          nflt=fl_l(ijp)
          nfll=fl_l(imj)
          nflb=fl_l(ijm)
          nfltr=fl_l(ipjp)
          nfllt=fl_l(imjp)
          nflbr=fl_l(ipjm)
          nfllb=fl_l(imjm)
!
! ***** right boundary conditions
!
          IF( nflr.NE.1 .AND. nfltr.NE.1 ) THEN

            n2 = ipj
            IF (ipj .GT. ncint .AND. time .EQ. dt) THEN
             WRITE(8,*) 'right bdry'
             WRITE(8,*) 'warning: boundary cell not belonging to proc', mpime
             WRITE(8,*) 'ij, i, j, imesh', ij, i, j, imesh
            END IF
            SELECT CASE (nflr) 

            CASE (2) 

              wg(n2)=wg(ij)
              DO is=1,nsolid
                ws(is,n2)=ws(is,ij)
              END DO

            CASE (3)

              wg(n2)=-wg(ij)
              DO is=1,nsolid
                IF(rlk(is,ij).GT.0.D0) ws(is,n2)=-ws(is,ij)
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
              t0n=tg(ij-1)
              t1n=tg(ij)
              t2n=tg(n2)
!
              IF(ug(ij).GT.0.D0) THEN
! ... OUTFLOW ...
!
! ... Non-reflecting boundary conditions
                u2     = ug(n2)
                u2nn   = u2 - dt/dr(nr) * (u2*u2 - u2n*ucn)
                ug(n2) = u2nn
                ucnn   = ucn - dt/drc * (u2n*ucn - u1n*uc1n)
                IF(j.EQ.2) ug(ipjm)=-ug(n2)
!
                t1nn = t1n - dt/dr1 * (ucn*t1n - uc1n*t0n)
!
! ... calculation of the gas volumetric fraction at time (n+1)dt
! ... from the mass balance equation of solids in cell (ij)
                epc=0.D0
                DO is=1,nsolid
                  eps = rlk(is,ij)*inrl(is) - dt*inrl(is)*inr(i)/dr1*          &
                        (rb(i)*us(is,ij)*rlk(is,ij) - rb(i-1)*us(is,imj)*rlk(is,imj))
                  epc=epc+eps
                END DO
                ep1nn=1.D0-epc
!
! ... calculation of the mixture density and fluid pressure 
! ... from a mass balance equation for the mixture.
! ... Use gas velocity for mixture.
                rm1nn=ep1nn*gmw(6)/(rgas*t1nn)
                rm1knn=0.D0
                DO is=1,nsolid
                  rm1knn=rm1knn+rlk(is,ij) 
                END DO
!
                rm2n=0.D0
                rm1n=0.D0
                rm0n=0.D0
                DO is=1,nsolid
                  rm2n=rm2n+rlk(is,n2)
                  rm1n=rm1n+rlk(is,ij)
                  rm0n=rm0n+rlk(is,imj)
                END DO
                rm2n=rm2n+p(n2)*ep(n2)*gmw(6)/(rgas*tg(n2))
                rm1n=rm1n+p(ij)*ep(ij)*gmw(6)/(rgas*tg(ij))
                rm0n=rm0n+p(imj)*ep(imj)*gmw(6)/(rgas*tg(imj))
!
                rmcn=(rm2n+rm1n)*0.5D0
                rmc1n=(rm1n+rm0n)*0.5D0
!
                rmcnn = rmcn-dt*indrc*inrb(i)*(r(i+1)*u2n*rmcn-r(i)*u1n*rmc1n)
!
                p1nn = (1.D0/rm1nn) *                           &
                       ( - rm1knn + rm1n - dt/dr1*inr(i) *      &
                       ( rb(i)*rm1n*ucn - rb(i-1)*rm0n*uc1n ) )
                IF (epc .LT. 1.0D-8) p1nn = p(ij)
!
! ... Calculation of the advanced-time fluid pressure from Momentum balance
! ... equation of the mixture
                p(n2) = -rmcnn*ucnn + rmcn*ucn - dt*indrc*inrb(i) *     &
                         ( r(i+1)*u2n*ucn*rmcn - r(i)*u1n*uc1n*rmc1n) + &
                          dt*p1nn*indrc
                p(n2)=p(n2)/(dt*indrc)
!
                epnn = ep(n2) - dt/dr(nr) * (ug(n2)*ep(n2) - ug(ij)*ep(ij))
                ep(n2) = epnn
                tgnn = tg(n2) - dt/dr(nr) * (ug(n2)*tg(n2) - ug(ij)*tg(ij))
                tg(n2) = tgnn
                sieg(n2)=sieg(ij)
!
! ... Correct non-physical pressure
                IF (p(n2).LT.0.0D0) p(n2) = p(ij)
!
              ELSE IF (ug(ij) .LT. 0.D0) THEN
! ... INFLOW ...
!
                ug(n2) = ug(ij)*rb(i)*inrb(i+1)
                ug2    = ug(n2)*rb(i+1)*inr(i+1)
!
                zrif=zb(j)+0.5D0*(dz(1)-dz(j))
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
! ... Correct non-physical pressure
                IF (p(n2).LT.0.0D0) p(n2) = prif
!
              ELSE IF (ug(ij) .EQ. 0.D0) THEN
!
                ug(n2) = ug(ij)
                p(n2)  = p(ij)
                ep(n2) = ep(ij)
                tg(n2) = tg(ij)
                sieg(n2)=sieg(ij)
              ENDIF
!
! ... Set primary variables
!
              IF(nfltr.EQ.3) ug(ipjp) = -ug(n2)
              IF(nfltr.EQ.4) wg(n2) = wg(ij)
              DO is=1,nsolid
                IF(nfltr.EQ.4) ws(is,n2)=ws(is,ij)
                IF(us(is,ij).GE.0.D0) THEN
                  rlk(is,n2)=rlk(is,ij)
                  us(is,n2)=us(is,ij)
                ELSE
                  rlk(is,n2)=0.0D0
                  us(is,n2) =0.0D0
                ENDIF
              END DO

              DO ig=1,ngas
                rgpgc(ig,n2)=rgpgc(ig,ij)
              END DO

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
          IF(nfll.NE.1.AND.nfllt.NE.1) THEN

            n2 = imj
            IF (imj .GT. ncint .AND. time .EQ. dt) THEN
             WRITE(8,*) 'left bdry'
             WRITE(8,*) 'warning: boundary cell not belonging to proc', mpime
             WRITE(8,*) 'ij, i, j, imesh', ij, i, j, imesh
            END IF
            SELECT CASE (nfll)
            CASE (2)
              wg(n2)=wg(ij)
              DO is=1,nsolid
                ws(is,n2)=ws(is,ij)
              END DO 
            CASE (3)
              wg(n2)=-wg(ij)
              DO is=1,nsolid
                IF(rlk(is,ij).GT.0.D0) ws(is,n2)=-ws(is,ij)
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
              epcn=(ep(ij)+ep(n2))*0.5D0
              tcn=(tg(ij)+tg(n2))*0.5D0
              tc1n=(tg(ij)+tg(ipj))*0.5D0
              t1n=tg(ij)
              t0n=tg(ij+1)
!
              u2n = ucn - (uc1n-ucn)/dr1 * 0.5*dr(i-1)
!
! ... OUTFLOW ...
              IF(ug(ij).LE.0.D0) THEN
!
! ... Non-reflecting boundary condition
!
                ug(n2) = ucn - dt/drc * (u1n*uc1n - u2n*ucn)                
                ucnn = ug(n2)
!
                t1nn = t1n - dt/dr1 * (uc1n*t0n - ucn*t1n)
!
! ... Mixture Density at time (n)dt
                rm2n=0.D0
                rm1n=0.D0
                rm0n=0.D0
                DO is=1,nsolid
                  rm2n=rm2n+rlk(is,n2)
                  rm1n=rm1n+rlk(is,ij)
                  rm0n=rm0n+rlk(is,ipj)
                END DO
                rm1knn = rm1n
                rm2n=rm2n+p(n2)*ep(n2)*gmw(6)/(rgas*tg(n2))
                rm1n=rm1n+p(ij)*ep(ij)*gmw(6)/(rgas*tg(ij))
                rm0n=rm0n+p(ipj)*ep(ipj)*gmw(6)/(rgas*tg(ipj))
!
                rmcn=(rm2n+rm1n)*0.5D0
                rmc1n=(rm1n+rm0n)*0.5D0
!
! ... calculation of the gas volumetric fraction at time (n+1)dt
! ... from the mass balance equation of solids in cell (ij)
                epc=0.D0
                DO is=1,nsolid
                 eps=rlk(is,ij)*inrl(is) -                                 &
                     dt*inrl(is)*inr(i)/dr1*(rb(i)*us(is,ij)*rlk(is,ipj) -  &
                     rb(i-1)*us(is,imj)*rlk(is,ij))
                 epc=epc+eps
                END DO
                ep1nn=1.D0-epc
!
! ... calculation of the mixture density and fluid pressure 
! ... from a mass balance equation for the mixture.
! ... Use gas velocity for mixture.
!
! ... Mixture Density at time (n+1)dt
                rmcnn=rmcn-dt*indrc*inrb(i-1)*(r(i)*u1n*rmc1n-r(i-1)*u2n*rmcn)
!
! ... Pressure at time (n+1)dt
                rm1nn = ep1nn*gmw(6)/(rgas*t1nn)
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
                epnn = ep(n2) - u2n*dt/dr(1) * (ep(ij) -ep(n2))
                ep(n2) = epnn
                tgnn = tg(n2) - u2n*dt/dr(1) * (tg(ij) - tg(n2))
                tg(n2) = tgnn
!
                sieg(n2)=sieg(ij)
!
! ... extrapolation of the temperature and solid fraction to time (n+1)dt
                DO is=1,nsolid
                  IF(nfllt.EQ.4) ws(is,n2)=ws(is,ij)
                  rlk(is,n2)=rlk(is,ij)
                  us(is,n2)=us(is,ij)
                END DO
! ... INFLOW ...
              ELSE
!
                ug(n2) = ug(ij)
                ug2 = ug(ij)
                zrif=zb(j)+0.5D0*(dz(1)-dz(j))
                CALL atm(zrif,prif,trif)
                rhorif=prif*gmw(6)/(rgas*trif)
                cost=prif/(rhorif**gammaair)
                costc=(gammaair*cost**(1.D0/gammaair))/(gammaair-1.D0)
!
! ... Adiabatic inflow
                p(n2)=(prif**gamn-(ug2**2)/(2.D0*costc))**(1.D0/gamn)
!
! ... Correct non-physical pressure
                IF (p(n2).LT.0.0D0) p(n2) = prif
!
                ep(n2) = 1.D0
                tg(n2)=trif
              ENDIF
!
! ... Set primary variables
!
              IF(j.EQ.2) ug(imjm)=-ug(n2)
              IF(nfllt.EQ.4) wg(n2)=wg(ij)
              IF(nfllt.EQ.3) ug(ipjp) = -ug(n2)
              DO is=1,nsolid
                IF(nfllt.EQ.4) ws(is,n2)=ws(is,ij)
                IF(us(is,ij).GE.0.D0) THEN
                  rlk(is,n2)=rlk(is,ij)
                  us(is,n2)=us(is,ij)
                ELSE
                  rlk(is,n2)=0.0D0
                  us(is,n2) =0.0D0
                ENDIF
              END DO
!
              DO ig=1,ngas
                rgpgc(ig,n2)=rgpgc(ig,ij)
              END DO
!
              DO is=1,nsolid
                sies(is,n2)=sies(is,ij)
              END DO
!
            CASE DEFAULT
              CONTINUE
            END SELECT
!
        END IF
!
!
! ****** top boundary conditions
!
          IF (nflt .NE. 1 .AND. nfltr .NE. 1) THEN 

            n2 = ijp
            IF (ijp .GT. ncint .AND. time .EQ. dt) THEN
             WRITE(8,*) 'top bdry'
             WRITE(8,*) 'warning: boundary cell not belonging to proc', mpime
             WRITE(8,*) 'ij, i, j, imesh', ij, i, j, imesh
            END IF

            SELECT CASE (nflt)

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
                IF(rlk(is,ij).GT.0.D0) us(is,n2)=-us(is,ij)
              END DO

            CASE (4)
!
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
              IF(wg(ij).GT.0.D0) THEN
!
! ... Non-reflecting boundary condition 
!
                w2    = wg(n2)
                w2nn   = w2 - dt/dz(nz) * (w2*w2 - w2n*wcn)
                wg(n2) = w2nn
                wcnn   = wcn - dt/dzc * (w2n*wcn - w1n*wc1n)
!
                t1nn = t1n - dt/dz1 * (wcn*t1n - wc1n*t0n)
!
! ... calculation of gas volumetric fraction at time (n+1)dt
! ... from the mass balance equation of solids, in cell (ij)
                epc=0.D0
                DO is=1,nsolid
                  eps=rlk(is,ij)*inrl(is) - dt*inrl(is)/dz1*(ws(is,ij)*rlk(is,ij) -  &
                      ws(is,ijm)*rlk(is,ijm))
                  epc=epc+eps
                END DO 
                ep1nn=1.D0-epc
!
! ... calculation of mixture density and fluid pressure at time (n+1)dt
! ... from the mixture mass balance equation in cell (ij)
! ... Use gas velocity for mixture.
                rm1nn=ep1nn*gmw(6)/(rgas*t1nn)
                rm1knn=0.D0
                DO is=1,nsolid
                  rm1knn=rm1knn+rlk(is,ij) 
                END DO 
                rm2n=0.D0
                rm1n=0.D0
                rm0n=0.D0
                DO is=1,nsolid
                  rm2n=rm2n+rlk(is,n2)
                  rm1n=rm1n+rlk(is,ij)
                  rm0n=rm0n+rlk(is,ijm)
                END DO 
                rm2n=rm2n+p(n2)*ep(n2)*gmw(6)/(rgas*tg(n2))
                rm1n=rm1n+p(ij)*ep(ij)*gmw(6)/(rgas*tg(ij))
                rm0n=rm0n+p(ijm)*ep(ijm)*gmw(6)/(rgas*tg(ijm))
                rmcn=(rm2n+rm1n)*0.5D0
                rmc1n=(rm1n+rm0n)*0.5D0
!
                rmcnn=rmcn-dt/dzc*(w2n*rmcn-w1n*rmc1n)
!
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
                ep(n2) = ep(ij)
                tg(n2) = 2.D0*tg(ij)-tg(ijm)
                sieg(n2)=sieg(ij)
!
! ... (Correct non-physical pressure)
                IF(p(n2).LE.0.D0) p(n2) = pnn2
!
              ELSEIF(wg(ij).LT.0.D0) THEN
! ... INFLOW ...
!
                wg2=wg(ij)
                wg(n2) = wg2
!
                t1nn = t1n - w1n*dt/dz1 * (t2n - t1n)
!
                zrif=zb(nz)+0.5D0*(dz(1)-dz(nz))
                CALL atm(zrif,prif,trif)
                rhorif=prif*gmw(6)/(rgas*trif)
                cost=prif/(rhorif**gammaair)
                costc=(gammaair*cost**(1.D0/gammaair))/(gammaair-1.D0)
                ep(n2)=1.D0
                tg(n2)=trif
                pnn2=p(n2)
                p(n2)=(prif**gamn-(wg2**2/2.D0)/costc)**(1.D0/gamn)
!
! ... (Correct non-physical pressure)
                IF(p(n2).LE.0.D0) p(n2)=pnn2
!
              ELSEIF(wg(ij).EQ.0.D0) THEN
                wg(n2) = wg(ij) 
                p(n2)  = p(ij)
                ep(n2) = ep(ij)
                tg(n2) = tg(ij)
                sieg(n2)=sieg(ij)
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
              IF(nfltr.EQ.4) ug(n2) = ug(ij)
!
              DO is=1,nsolid
                IF(nfltr.EQ.4) us(is,n2)=us(is,ij)
                IF(ws(is,ij).GE.0.D0) THEN
                  rlk(is,n2)=rlk(is,ij)
                  ws(is,n2)=ws(is,ij)
                ELSE
                  rlk(is,n2)=0.D0
                  ws(is,n2)=0.D0
                ENDIF
                IF(i .EQ. (nr-1) .AND. j .EQ. (nz-1)) THEN
                  ws(is,ipjp)=ws(is,n2)
                  us(is,ipjp)=us(is,n2)
                ENDIF
              END DO 

              DO ig=1,ngas
                rgpgc(ig,n2)=rgpgc(ig,ij)
              END DO 
             
              DO is=1,nsolid
                sies(is,n2)=sies(is,ij)
                ts(is,n2)=ts(is,ij)
                IF (i .EQ. (nr-1) .AND. j .EQ. (nz-1)) THEN
                  ts(is,ipjp)=ts(is,n2)
                END IF
              END DO
!
            CASE DEFAULT
              CONTINUE
            END SELECT

          END IF
!
! ***** bottom boundary conditions
!
          IF(nflb.NE.1.AND.nflbr.NE.1) THEN
!
            n2 = ijm
            IF (ijm .GT. ncint .AND. time .EQ. dt) THEN
             WRITE(8,*) 'bottom bdry'
             WRITE(8,*) 'warning: boundary cell not belonging to proc', mpime
             WRITE(8,*) 'ij, i, j, imesh', ij, i, j, imesh
            END IF
            SELECT CASE (nflb)

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
                IF(rlk(is,ij).GT.0.D0) us(is,n2)=-us(is,ij)
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


!-----------------------------------------------------------------------
      SUBROUTINE fboundary(fug, fwg)
!-----------------------------------------------------------------------
! ... This routine computes boundary conditions for filtered velocities
!
      USE grid, ONLY: fl_l, ncint
      USE set_indexes
!
      IMPLICIT NONE
      REAL*8, INTENT(INOUT), DIMENSION(:) :: fug, fwg
      INTEGER :: ij
      INTEGER :: n2, nflr, nflt, nfll, nflb
      INTEGER :: nflbr, nfltr, nfllt, nfllb
!
      DO ij = 1, ncint
        IF (fl_l(ij) .EQ. 1) THEN
          CALL subscr(ij)
!
          nflr=fl_l(ipj)
          nflt=fl_l(ijp)
          nflb=fl_l(ijm)
          nfll=fl_l(imj)
          nfltr=fl_l(ipjp)
          nfllt=fl_l(imjp)
          nflbr=fl_l(ipjm)
          nfllb=fl_l(imjm)
!
! ***** right boundary conditions
!
          IF( nflr.NE.1 .AND. nfltr.NE.1 ) THEN
            n2 = ipj
            SELECT CASE (nflr) 
            CASE (2) 
              fwg(n2)=fwg(ij)
            CASE (3)
              fwg(n2)=-fwg(ij)
            CASE (4)
              fug(n2) = fug(ij)
              fwg(n2) = fwg(ij)
            CASE DEFAULT
              CONTINUE
            END SELECT
          END IF
!
! ***** top boundary conditions
!
          IF (nflt .NE. 1 .AND. nfltr .NE. 1) THEN 
            n2 = ijp
            SELECT CASE (nflt)
            CASE (2)
              fug(n2)=fug(ij)
            CASE (3)
              fug(n2)=-fug(ij)
            CASE (4)
              fwg(n2) = fwg(ij) 
              fug(n2) = fug(ij)
            CASE DEFAULT
              CONTINUE
            END SELECT
          END IF
!
! ***** left boundary conditions
!
          IF(nfll.NE.1.AND.nfllt.NE.1) THEN
            n2 = imj
            SELECT CASE (nfll)
            CASE (2)
              fwg(n2)=fwg(ij)
            CASE (3)
              fwg(n2)=-fwg(ij)
            CASE (4)
              fwg(n2) = fwg(ij) 
              fug(n2) = fug(ij)
            CASE DEFAULT
              CONTINUE
            END SELECT
          END IF
!
! ***** bottom boundary conditions
!
          IF(nflb.NE.1.AND.nflbr.NE.1) THEN
            n2 = ijm
            SELECT CASE (nflb)
            CASE (2)
              fug(n2)=fug(ij)
            CASE (3)
              fug(n2)=-fug(ij)
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
      END SUBROUTINE fboundary
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
      SUBROUTINE boundary3d
!-----------------------------------------------------------------------
! ... This routine computes boundary conditions 
!
      USE atmosphere, ONLY: gravx, gravz, atm
      USE dimensions
      USE eos_gas, ONLY: rgpgc
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

      INTEGER :: ijk, i, j, k, imesh, ig, is
      INTEGER :: n2
!
      prif=0.D0; pnn2=0.D0; p1nn=0.D0
      zrif=0.D0; trif=0.D0; rhorif=0.D0; cost=0.D0; costc=0.D0
      rm1n=0.D0; rmcn=0.D0; rmcnn=0.D0; rmc1n=0.D0; rm0n=0.D0; rm1nn=0.D0; rm1knn=0.D0; rm2n=0.D0
      w2n=0.D0; wg2=0.D0; wcnn=0.D0; wcn=0.D0; w1n=0.D0; wc1n=0.D0
      u2nn=0.D0; u2=0.D0; w2nn=0.D0; w2=0.D0
      ucn=0.D0; u1n=0.D0; uc1n=0.D0; ug2=0.D0; u2n=0.D0; ucnn=0.D0
      eps=0.D0; epc=0.D0; epcn=0.D0; ep1nn=0.D0; epnn=0.D0
      t1nn=0.D0; tc1n=0.D0; tcn=0.D0; t1n=0.D0; t0n=0.D0; t2n=0.D0; tgnn=0.D0
      dz1=0.D0; dzc=0.D0; 
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
              IF( ug( ijk ) > 0.D0 ) THEN

                ! ...  OUTFLOW ...
!
                ! ...  Non-reflecting boundary conditions

                u2     = ug( n2 )
                u2nn   = u2 - dt / dx(nx) * (u2*u2 - u2n*ucn)
                ug(n2) = u2nn
                ucnn   = ucn - dt / dxc * (u2n*ucn - u1n*uc1n)

                ! MODIFICARE X3D ....  IF(j == .EQ.2) ug(ipjm)=-ug(n2) ! 

                t1nn = t1n - dt/dx1 * (ucn*t1n - uc1n*t0n)
!
                ! ... calculation of the gas volumetric fraction at time (n+1)dt
                ! ... from the mass balance equation of solids in cell (ijk)

                epc = 0.D0
                DO is = 1, nsolid
                  eps = rlk(is,ijk) * inrl(is) - ( dt * inrl(is) / dx1 )  *          &
                        ( us(is,ijk)*rlk(is,ijk) - us(is,imjk)*rlk(is,imjk) )
                  epc = epc + eps
                END DO
                ep1nn = 1.D0 - epc
!
                ! ... calculation of the mixture density and fluid pressure 
                ! ... from a mass balance equation for the mixture.
                ! ... Use gas velocity for mixture.

                rm1nn = ep1nn * gmw(6) / (rgas*t1nn)
                rm1knn= 0.D0
                DO is = 1, nsolid
                  rm1knn = rm1knn + rlk(is,ijk) 
                END DO
!
                rm2n=0.D0
                rm1n=0.D0
                rm0n=0.D0
                DO is=1,nsolid
                  rm2n=rm2n+rlk(is,n2)
                  rm1n=rm1n+rlk(is,ijk)
                  rm0n=rm0n+rlk(is,imjk)
                END DO
                rm2n=rm2n+p(n2)*ep(n2)*gmw(6)/(rgas*tg(n2))
                rm1n=rm1n+p(ijk)*ep(ijk)*gmw(6)/(rgas*tg(ijk))
                rm0n=rm0n+p(imjk)*ep(imjk)*gmw(6)/(rgas*tg(imjk))
!
                rmcn=(rm2n+rm1n)*0.5D0
                rmc1n=(rm1n+rm0n)*0.5D0
!
                rmcnn = rmcn-dt/dxc*(u2n*rmcn-u1n*rmc1n)
!
                p1nn = (1.D0/rm1nn) * ( - rm1knn + rm1n - dt/dx1 * ( rm1n*ucn - rm0n*uc1n ) )
                IF (epc < 1.0D-8) p1nn = p(ijk)
!
                ! ... Calculation of the advanced-time fluid pressure from Momentum balance
                ! ... equation of the mixture

                p(n2) = -rmcnn*ucnn + rmcn*ucn - dt/dxc * ( u2n*ucn*rmcn - u1n*uc1n*rmc1n) + dt/dxc*p1nn
                p(n2) = p(n2)/(dt/dxc)
!
                epnn = ep(n2) - dt/dx(nx) * (ug(n2)*ep(n2) - ug(ijk)*ep(ijk))
                ep(n2) = epnn
                tgnn = tg(n2) - dt/dx(nx) * (ug(n2)*tg(n2) - ug(ijk)*tg(ijk))
                tg(n2) = tgnn
                sieg(n2)=sieg(ijk)
!
                ! ... Correct non-physical pressure
                IF (p(n2) < 0.0D0) p(n2) = p(ijk)
!
              ELSE IF ( ug(ijk) <  0.D0 ) THEN

                ! ... INFLOW ...
!
                ug(n2) = ug(ijk)
                ug2    = ug(n2)
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

                ! ... Correct non-physical pressure

                IF ( p(n2) < 0.0D0 ) p(n2) = prif
!
              ELSE IF ( ug(ijk) == 0.D0 ) THEN
!
                ug(n2)   = ug(ijk)
                p(n2)    = p(ijk)
                ep(n2)   = ep(ijk)
                tg(n2)   = tg(ijk)
                sieg(n2) = sieg(ijk)

              ENDIF
!
              ! ... Set primary variables
!
              ! .. IF(nfltr.EQ.3) ug(ipjp) = -ug(n2)
              ! .. IF(nfltr.EQ.4) wg(n2) = wg(ijk)

              DO is = 1, nsolid
                ! .. IF(nfltr.EQ.4) ws(is,n2)=ws(is,ijk)
                IF( us(is,ijk) >= 0.D0 ) THEN
                  rlk(is,n2) = rlk(is,ijk)
                  us(is,n2)  = us(is,ijk)
                ELSE
                  rlk(is,n2) = 0.0D0
                  us(is,n2)  = 0.0D0
                ENDIF
              END DO

              DO ig = 1, ngas
                rgpgc(ig,n2) = rgpgc(ig,ijk)
              END DO

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

              IF( ug(ijk) <= 0.D0 ) THEN
!
                ! ... Non-reflecting boundary condition
!
                ug(n2) = ucn - dt/dxc * (u1n*uc1n - u2n*ucn)                
                ucnn = ug(n2)
!
                t1nn = t1n - dt/dx1 * (uc1n*t0n - ucn*t1n)
!
                ! ... Mixture Density at time (n)dt

                rm2n=0.D0
                rm1n=0.D0
                rm0n=0.D0
                DO is=1,nsolid
                  rm2n=rm2n+rlk(is,n2)
                  rm1n=rm1n+rlk(is,ijk)
                  rm0n=rm0n+rlk(is,ipjk)
                END DO
                rm1knn = rm1n
                rm2n=rm2n+p(n2)*ep(n2)*gmw(6)/(rgas*tg(n2))
                rm1n=rm1n+p(ijk)*ep(ijk)*gmw(6)/(rgas*tg(ijk))
                rm0n=rm0n+p(ipjk)*ep(ipjk)*gmw(6)/(rgas*tg(ipjk))
!
                rmcn=(rm2n+rm1n)*0.5D0
                rmc1n=(rm1n+rm0n)*0.5D0
!
                ! ... calculation of the gas volumetric fraction at time (n+1)dt
                ! ... from the mass balance equation of solids in cell (ijk)

                epc = 0.D0
                DO is = 1, nsolid
                 eps = rlk(is,ijk) * inrl(is) -                                 &
                     dt*inrl(is)/dx1*(us(is,ijk)*rlk(is,ipjk) - us(is,imjk)*rlk(is,ijk))
                 epc=epc+eps
                END DO
                ep1nn=1.D0-epc
!
                ! ... calculation of the mixture density and fluid pressure 
                ! ... from a mass balance equation for the mixture.
                ! ... Use gas velocity for mixture.
!
                ! ... Mixture Density at time (n+1)dt

                rmcnn = rmcn - dt/dxc * (u1n*rmc1n-u2n*rmcn)
!
                ! ... Pressure at time (n+1)dt

                rm1nn = ep1nn*gmw(6)/(rgas*t1nn)
                p1nn = (1.D0/rm1nn) * &
                       (-rm1knn+rm1n-dt/dx1*(uc1n*rm0n-rm1n*ucn))
                IF (epc .LT. 1.0D-8) p1nn = p(ijk)
!
                ! ... Calculation of the advanced-time fluid pressure from Momentum balance
                ! ... equation of the mixture

                p(n2) = -rmcnn*ucnn + rmcn*ucn - dt/dxc *     &
                         ( u1n*uc1n*rmc1n - u2n*ucn*rmcn) + dt/dxc * p1nn
                p(n2)=p(n2)/(dt/dxc)
!
                epnn = ep(n2) - u2n*dt/dx(1) * (ep(ijk) -ep(n2))
                ep(n2) = epnn
                tgnn = tg(n2) - u2n*dt/dx(1) * (tg(ijk) - tg(n2))
                tg(n2) = tgnn
!
                sieg(n2)=sieg(ijk)
!
                ! ... extrapolation of the temperature and solid fraction to time (n+1)dt

                DO is=1,nsolid
                  ! IF(nfllt.EQ.4) ws(is,n2)=ws(is,ijk)
                  rlk(is,n2)=rlk(is,ijk)
                  us(is,n2)=us(is,ijk)
                END DO

              ELSE IF( ug(ijk) > 0.D0 ) THEN

                ! ... INFLOW ...
!
                ug(n2) = ug(ijk)
                ug2    = ug(ijk)
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

              ENDIF
!
              ! ... Set primary variables
!
              ! .. IF(nfllt.EQ.4) wg(n2)=wg(ijk)
              ! .. IF(nfllt.EQ.3) ug(ipjp) = -ug(n2)
              DO is=1,nsolid
                ! .. IF(nfllt.EQ.4) ws(is,n2)=ws(is,ijk)
                IF( us(is,ijk) >= 0.D0 ) THEN
                  rlk(is,n2)=rlk(is,ijk)
                  us(is,n2)=us(is,ijk)
                ELSE
                  rlk(is,n2)=0.0D0
                  us(is,n2) =0.0D0
                ENDIF
              END DO
!
              DO ig=1,ngas
                rgpgc(ig,n2)=rgpgc(ig,ijk)
              END DO
!
              DO is=1,nsolid
                sies(is,n2)=sies(is,ijk)
              END DO
!
            CASE DEFAULT

              CONTINUE

            END SELECT
!
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

              IF( wg( ijk ) > 0.D0) THEN
!
                ! ... Non-reflecting boundary condition 
!
                w2    = wg(n2)
                w2nn   = w2 - dt/dz(nz) * (w2*w2 - w2n*wcn)
                wg(n2) = w2nn
                wcnn   = wcn - dt/dzc * (w2n*wcn - w1n*wc1n)
!
                t1nn = t1n - dt/dz1 * (wcn*t1n - wc1n*t0n)
!
                ! ... calculation of gas volumetric fraction at time (n+1)dt
                ! ... from the mass balance equation of solids, in cell (ijk)

                epc = 0.D0
                DO is = 1, nsolid
                  eps = rlk(is,ijk)*inrl(is) - dt*inrl(is)/dz1*(ws(is,ijk)*rlk(is,ijk) -  &
                      ws(is,ijkm)*rlk(is,ijkm))
                  epc = epc + eps
                END DO 
                ep1nn = 1.D0 - epc
!
                ! ... calculation of mixture density and fluid pressure at time (n+1)dt
                ! ... from the mixture mass balance equation in cell (ijk)
                ! ... Use gas velocity for mixture.

                rm1nn  = ep1nn*gmw(6)/(rgas*t1nn)
                rm1knn = 0.D0
                DO is = 1, nsolid
                  rm1knn = rm1knn + rlk(is,ijk) 
                END DO 
                rm2n=0.D0
                rm1n=0.D0
                rm0n=0.D0
                DO is=1,nsolid
                  rm2n=rm2n+rlk(is,n2)
                  rm1n=rm1n+rlk(is,ijk)
                  rm0n=rm0n+rlk(is,ijkm)
                END DO 
                rm2n=rm2n+p(n2)*ep(n2)*gmw(6)/(rgas*tg(n2))
                rm1n=rm1n+p(ijk)*ep(ijk)*gmw(6)/(rgas*tg(ijk))
                rm0n=rm0n+p(ijkm)*ep(ijkm)*gmw(6)/(rgas*tg(ijkm))
                rmcn=(rm2n+rm1n)*0.5D0
                rmc1n=(rm1n+rm0n)*0.5D0
!
                rmcnn=rmcn-dt/dzc*(w2n*rmcn-w1n*rmc1n)
!
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
                tg(n2) = 2.D0*tg(ijk)-tg(ijkm)
                sieg(n2) = sieg(ijk)
!
                ! ... (Correct non-physical pressure)

                IF( p(n2) <= 0.D0 ) p(n2) = pnn2
!
              ELSEIF( wg(ijk) < 0.D0 ) THEN

                ! ... INFLOW ...
!
                wg2=wg(ijk)
                wg(n2) = wg2
!
                t1nn = t1n - w1n*dt/dz1 * (t2n - t1n)
!
                zrif=zb(nz)+0.5D0*(dz(1)-dz(nz))
                CALL atm(zrif,prif,trif)
                rhorif=prif*gmw(6)/(rgas*trif)
                cost=prif/(rhorif**gammaair)
                costc=(gammaair*cost**(1.D0/gammaair))/(gammaair-1.D0)
                ep(n2)=1.D0
                tg(n2)=trif
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
                sieg(n2)=sieg(ijk)

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
                  rlk(is,n2)=rlk(is,ijk)
                  ws(is,n2)=ws(is,ijk)
                ELSE
                  rlk(is,n2)=0.D0
                  ws(is,n2)=0.D0
                ENDIF
                ! IF(i .EQ. (nr-1) .AND. j .EQ. (nz-1)) THEN
                !   ws(is,ipjp)=ws(is,n2)
                !   us(is,ipjp)=us(is,n2)
                ! ENDIF
              END DO 

              DO ig=1,ngas
                rgpgc(ig,n2)=rgpgc(ig,ijk)
              END DO 
             
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


!-----------------------------------------------------------------------
      SUBROUTINE fboundary3d(fug, fwg)
!-----------------------------------------------------------------------
! ... This routine computes boundary conditions for filtered velocities
!
      USE grid, ONLY: fl_l, ncint
      USE set_indexes, ONLY: subscr, ipjk,  imjk,  ippjk,  immjk,  ijpk,  &
        ipjpk,  imjpk,  ijmk,  ipjmk,  imjmk,  ijppk,  ijmmk,  ijkp,  ipjkp,  &
        imjkp,  ijpkp,  ijmkp,  ijkm,  ipjkm,  imjkm,  ijpkm,  ijmkm,  ijkpp,  ijkmm
!
      IMPLICIT NONE
      REAL*8, INTENT(INOUT), DIMENSION(:) :: fug, fwg
      INTEGER :: ijk
      INTEGER :: n2, nflr, nflt, nfll, nflb
      INTEGER :: nflbr, nfltr, nfllt, nfllb
!
      DO ijk = 1, ncint

        IF ( fl_l(ijk) == 1 ) THEN

          CALL subscr(ijk)
!
          ! ***** East boundary conditions ***** !
!
          n2 = ipjk
          SELECT CASE ( fl_l( n2 ) ) 
          CASE (2) 
            fwg(n2)=fwg(ijk)
          CASE (3)
            fwg(n2)=-fwg(ijk)
          CASE (4)
            fug(n2) = fug(ijk)
            fwg(n2) = fwg(ijk)
          CASE DEFAULT
            CONTINUE
          END SELECT
!
          ! ***** top boundary conditions ***** !
!
          n2 = ijkp
          SELECT CASE ( fl_l( n2 ) )
          CASE (2)
            fug(n2)=fug(ijk)
          CASE (3)
            fug(n2)=-fug(ijk)
          CASE (4)
            fwg(n2) = fwg(ijk) 
            fug(n2) = fug(ijk)
          CASE DEFAULT
            CONTINUE
          END SELECT
!
          ! ***** West boundary conditions ***** !
!
          n2 = imjk
          SELECT CASE ( fl_l( n2 ) )
          CASE (2)
            fwg(n2)=fwg(ijk)
          CASE (3)
            fwg(n2)=-fwg(ijk)
          CASE (4)
            fwg(n2) = fwg(ijk) 
            fug(n2) = fug(ijk)
          CASE DEFAULT
            CONTINUE
          END SELECT
!
          ! ***** bottom boundary conditions ***** !
!
          n2 = ijkm
          SELECT CASE ( fl_l( n2 ) )
          CASE (2)
            fug(n2)=fug(ijk)
          CASE (3)
            fug(n2)=-fug(ijk)
          CASE DEFAULT
            CONTINUE
          END SELECT
!
        END IF

      END DO
!
      RETURN
!-----------------------------------------------------------------------
      END SUBROUTINE fboundary3d
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
      END MODULE boundary_conditions
!-----------------------------------------------------------------------
