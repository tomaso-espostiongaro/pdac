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
      USE gas_solid_temperature, ONLY: sieg, tg, siek, tk
      USE gas_solid_velocity, ONLY: ug, vg, uk, vk
      USE grid, ONLY: zb, dz, rb, dr, r, inr, inrb
      USE grid, ONLY: nij_l, myij
      USE grid, ONLY: fl_l
      USE parallel, ONLY: mpime
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: p, ep
      USE reactions, ONLY: irex
      USE set_indexes
      USE time_parameters, ONLY: dt, time
!
      IMPLICIT NONE
!
      REAL*8 :: prif, pnn2, p1nn
      REAL*8 :: zrif, trif, rhorif, cost, costc
      REAL*8 :: rm1n, rmcn, rmcnn, rmc1n, rm0n, rm1nn, rm1knn, rm2n
      REAL*8 :: v2n, vg2, vcnn, vcn, v1n, vc1n
      REAL*8 :: u2nn, u2, v2nn, v2
      REAL*8 :: ucn, u1n, uc1n, ug2, u2n, ucnn
      REAL*8 :: epk, epc, epcn, ep1nn, epnn
      REAL*8 :: t1nn, tc1n, tcn, t1n, t0n, t2n, tgnn
      REAL*8 :: dz1, dzc, dr1, drc, indrc

      INTEGER :: ij, i, j, imesh, k, kg
      INTEGER :: n2, nflr, nflt, nfll, nflb
      INTEGER :: nflbr, nfltr, nfllt, nfllb
      INTEGER :: mm
!
      DO ij = 1, nij_l
        imesh = myij(0, 0, ij)
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
            IF (ipj .GT. nij_l .AND. time .EQ. dt) THEN
             WRITE(8,*) 'right bdry'
             WRITE(8,*) 'warning: boundary cell not belonging to proc', mpime
             WRITE(8,*) 'ij, i, j, imesh', ij, i, j, imesh
            END IF
            SELECT CASE (nflr) 

            CASE (2) 

              vg(n2)=vg(ij)
              DO k=1,nsolid
                vk(k,n2)=vk(k,ij)
              END DO

            CASE (3)

              vg(n2)=-vg(ij)
              DO k=1,nsolid
                IF(rlk(k,ij).GT.0.D0) vk(k,n2)=-vk(k,ij)
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
                DO k=1,nsolid
                  epk = rlk(k,ij)*inrl(k) - dt*inrl(k)*inr(i)/dr1*          &
                        (rb(i)*uk(k,ij)*rlk(k,ij) - rb(i-1)*uk(k,imj)*rlk(k,imj))
                  epc=epc+epk
                END DO
                ep1nn=1.D0-epc
!
! ... calculation of the mixture density and fluid pressure 
! ... from a mass balance equation for the mixture.
! ... Use gas velocity for mixture.
                rm1nn=ep1nn*gmw(6)/(rgas*t1nn)
                rm1knn=0.D0
                DO k=1,nsolid
                  rm1knn=rm1knn+rlk(k,ij) 
                END DO
!
                rm2n=0.D0
                rm1n=0.D0
                rm0n=0.D0
                DO k=1,nsolid
                  rm2n=rm2n+rlk(k,n2)
                  rm1n=rm1n+rlk(k,ij)
                  rm0n=rm0n+rlk(k,imj)
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
              IF(nfltr.EQ.4) vg(n2) = vg(ij)
              DO k=1,nsolid
                IF(nfltr.EQ.4) vk(k,n2)=vk(k,ij)
                IF(uk(k,ij).GE.0.D0) THEN
                  rlk(k,n2)=rlk(k,ij)
                  uk(k,n2)=uk(k,ij)
                ELSE
                  rlk(k,n2)=0.0D0
                  uk(k,n2) =0.0D0
                ENDIF
              END DO

              DO kg=1,ngas
                rgpgc(kg,n2)=rgpgc(kg,ij)
              END DO

              DO k=1,nsolid
                siek(k,n2)=siek(k,ij)
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
            IF (imj .GT. nij_l .AND. time .EQ. dt) THEN
             WRITE(8,*) 'left bdry'
             WRITE(8,*) 'warning: boundary cell not belonging to proc', mpime
             WRITE(8,*) 'ij, i, j, imesh', ij, i, j, imesh
            END IF
            SELECT CASE (nfll)
            CASE (2)
              vg(n2)=vg(ij)
              DO k=1,nsolid
                vk(k,n2)=vk(k,ij)
              END DO 
            CASE (3)
              vg(n2)=-vg(ij)
              DO k=1,nsolid
                IF(rlk(k,ij).GT.0.D0) vk(k,n2)=-vk(k,ij)
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
                DO k=1,nsolid
                  rm2n=rm2n+rlk(k,n2)
                  rm1n=rm1n+rlk(k,ij)
                  rm0n=rm0n+rlk(k,ipj)
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
                DO k=1,nsolid
                 epk=rlk(k,ij)*inrl(k) -                                 &
                     dt*inrl(k)*inr(i)/dr1*(rb(i)*uk(k,ij)*rlk(k,ipj) -  &
                     rb(i-1)*uk(k,imj)*rlk(k,ij))
                 epc=epc+epk
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
                DO k=1,nsolid
                  IF(nfllt.EQ.4) vk(k,n2)=vk(k,ij)
                  rlk(k,n2)=rlk(k,ij)
                  uk(k,n2)=uk(k,ij)
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
              IF(nfllt.EQ.4) vg(n2)=vg(ij)
              IF(nfllt.EQ.3) ug(ipjp) = -ug(n2)
              DO k=1,nsolid
                IF(nfllt.EQ.4) vk(k,n2)=vk(k,ij)
                IF(uk(k,ij).GE.0.D0) THEN
                  rlk(k,n2)=rlk(k,ij)
                  uk(k,n2)=uk(k,ij)
                ELSE
                  rlk(k,n2)=0.0D0
                  uk(k,n2) =0.0D0
                ENDIF
              END DO
!
              DO kg=1,ngas
                rgpgc(kg,n2)=rgpgc(kg,ij)
              END DO
!
              DO k=1,nsolid
                siek(k,n2)=siek(k,ij)
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
            IF (ijp .GT. nij_l .AND. time .EQ. dt) THEN
             WRITE(8,*) 'top bdry'
             WRITE(8,*) 'warning: boundary cell not belonging to proc', mpime
             WRITE(8,*) 'ij, i, j, imesh', ij, i, j, imesh
            END IF

            SELECT CASE (nflt)

            CASE (2)

              ug(n2)=ug(ij)
              DO k=1,nsolid
                uk(k,n2)=uk(k,ij)
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
              DO k=1,nsolid
                IF(rlk(k,ij).GT.0.D0) uk(k,n2)=-uk(k,ij)
              END DO

            CASE (4)
!
! ... definitions and linear interpolations
              dzc=(dz(nz)+dz(nz-1))*0.5D0
              dz1=dz(nz-1)
!
              vc1n=vg(ijm)
              vcn=vg(ij)
              v1n=(vg(ijm)+vg(ij))*0.5D0
              v2n=(vg(n2)+vg(ij))*0.5D0
              epcn=(ep(ij)+ep(n2))*0.5D0
              tc1n=(tg(ijm)+tg(ij))*0.5D0
              tcn=(tg(ij)+tg(n2))*0.5D0
              t2n=tg(n2)
              t1n=tg(ij)
              t0n=tg(ij-nz)
!
! ... OUTFLOW ...
              IF(vg(ij).GT.0.D0) THEN
!
! ... Non-reflecting boundary condition 
!
                v2    = vg(n2)
                v2nn   = v2 - dt/dz(nz) * (v2*v2 - v2n*vcn)
                vg(n2) = v2nn
                vcnn   = vcn - dt/dzc * (v2n*vcn - v1n*vc1n)
!
                t1nn = t1n - dt/dr1 * (vcn*t1n - vc1n*t0n)
!
! ... calculation of gas volumetric fraction at time (n+1)dt
! ... from the mass balance equation of solids, in cell (ij)
                epc=0.D0
                DO k=1,nsolid
                  epk=rlk(k,ij)*inrl(k) - dt*inrl(k)/dz1*(vk(k,ij)*rlk(k,ij) -  &
                      vk(k,ijm)*rlk(k,ijm))
                  epc=epc+epk
                END DO 
                ep1nn=1.D0-epc
!
! ... calculation of mixture density and fluid pressure at time (n+1)dt
! ... from the mixture mass balance equation in cell (ij)
! ... Use gas velocity for mixture.
                rm1nn=ep1nn*gmw(6)/(rgas*t1nn)
                rm1knn=0.D0
                DO k=1,nsolid
                  rm1knn=rm1knn+rlk(k,ij) 
                END DO 
                rm2n=0.D0
                rm1n=0.D0
                rm0n=0.D0
                DO k=1,nsolid
                  rm2n=rm2n+rlk(k,n2)
                  rm1n=rm1n+rlk(k,ij)
                  rm0n=rm0n+rlk(k,ijm)
                END DO 
                rm2n=rm2n+p(n2)*ep(n2)*gmw(6)/(rgas*tg(n2))
                rm1n=rm1n+p(ij)*ep(ij)*gmw(6)/(rgas*tg(ij))
                rm0n=rm0n+p(ijm)*ep(ijm)*gmw(6)/(rgas*tg(ijm))
                rmcn=(rm2n+rm1n)*0.5D0
                rmc1n=(rm1n+rm0n)*0.5D0
!
                rmcnn=rmcn-dt/dzc*(v2n*rmcn-v1n*rmc1n)
!
                p1nn=(1.D0/rm1nn)*(-rm1knn+rm1n-dt/dz1*(rm1n*vcn-rm0n*vc1n))
                IF (epc .LT. 1.0D-8) p1nn = p(ij)
!
! ... Calculation of the advanced-time fluid pressure from Momentum balance
! ... equation of the mixture
! ... 
                pnn2=p(n2)
!
                p(n2)=-rmcnn*vcnn + rmcn*vcn-dt/dzc*(v2n*vcn*rmcn-v1n*vc1n*rmc1n) + &
                       dt*p1nn/dzc+gravz*dt*rmcn
                p(n2)=p(n2)/(dt/dzc)
                ep(n2) = ep(ij)
                tg(n2) = 2.D0*tg(ij)-tg(ijm)
                sieg(n2)=sieg(ij)
!
! ... (Correct non-physical pressure)
                IF(p(n2).LE.0.D0) p(n2) = pnn2
!
              ELSEIF(vg(ij).LT.0.D0) THEN
! ... INFLOW ...
!
                vg2=vg(ij)
                vg(n2) = vg2
!
                t1nn = t1n - v1n*dt/dr1 * (t2n - t1n)
!
                zrif=zb(nz)+0.5D0*(dz(1)-dz(nz))
                CALL atm(zrif,prif,trif)
                rhorif=prif*gmw(6)/(rgas*trif)
                cost=prif/(rhorif**gammaair)
                costc=(gammaair*cost**(1.D0/gammaair))/(gammaair-1.D0)
                ep(n2)=1.D0
                tg(n2)=trif
                pnn2=p(n2)
                p(n2)=(prif**gamn-(vg2**2/2.D0)/costc)**(1.D0/gamn)
!
! ... (Correct non-physical pressure)
                IF(p(n2).LE.0.D0) p(n2)=pnn2
!
              ELSEIF(vg(ij).EQ.0.D0) THEN
                vg(n2) = vg(ij) 
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
                  vg(ipjp) = vg(n2)
                  ug(ipjp) = ug(n2)
                ELSE IF(i .EQ. 2) THEN
                  vg(imjp) = vg(n2)
                  ug(imjp) = ug(n2)
                ENDIF
              END IF
!
              IF(nfltr.EQ.4) ug(n2) = ug(ij)
!
              DO k=1,nsolid
                IF(nfltr.EQ.4) uk(k,n2)=uk(k,ij)
                IF(vk(k,ij).GE.0.D0) THEN
                  rlk(k,n2)=rlk(k,ij)
                  vk(k,n2)=vk(k,ij)
                ELSE
                  rlk(k,n2)=0.D0
                  vk(k,n2)=0.D0
                ENDIF
                IF(i .EQ. (nr-1) .AND. j .EQ. (nz-1)) THEN
                  vk(k,ipjp)=vk(k,n2)
                  uk(k,ipjp)=uk(k,n2)
                ENDIF
              END DO 

              DO kg=1,ngas
                rgpgc(kg,n2)=rgpgc(kg,ij)
              END DO 
             
              DO k=1,nsolid
                siek(k,n2)=siek(k,ij)
                tk(k,n2)=tk(k,ij)
                IF (i .EQ. (nr-1) .AND. j .EQ. (nz-1)) THEN
                  tk(k,ipjp)=tk(k,n2)
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
            IF (ijm .GT. nij_l .AND. time .EQ. dt) THEN
             WRITE(8,*) 'bottom bdry'
             WRITE(8,*) 'warning: boundary cell not belonging to proc', mpime
             WRITE(8,*) 'ij, i, j, imesh', ij, i, j, imesh
            END IF
            SELECT CASE (nflb)

            CASE (2)

              ug(n2)=ug(ij)
              DO k=1,nsolid
                uk(k,n2)=uk(k,ij)
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
              DO k=1,nsolid
                IF(rlk(k,ij).GT.0.D0) uk(k,n2)=-uk(k,ij)
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
      SUBROUTINE fboundary(fug, fvg)
!-----------------------------------------------------------------------
! ... This routine computes boundary conditions for filtered velocities
!
      USE grid, ONLY: fl_l, nij_l
      USE set_indexes
!
      IMPLICIT NONE
      REAL*8, INTENT(INOUT), DIMENSION(:) :: fug, fvg
      INTEGER :: ij
      INTEGER :: n2, nflr, nflt, nfll, nflb
      INTEGER :: nflbr, nfltr, nfllt, nfllb
!
      DO ij = 1, nij_l
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
              fvg(n2)=fvg(ij)
            CASE (3)
              fvg(n2)=-fvg(ij)
            CASE (4)
              fug(n2) = fug(ij)
              fvg(n2) = fvg(ij)
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
              fvg(n2) = fvg(ij) 
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
              fvg(n2)=fvg(ij)
            CASE (3)
              fvg(n2)=-fvg(ij)
            CASE (4)
              fvg(n2) = fvg(ij) 
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
      END MODULE boundary_conditions
!-----------------------------------------------------------------------
