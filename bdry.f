!-----------------------------------------------------------------------
      SUBROUTINE bdry
!-----------------------------------------------------------------------
! ... This routine computes boundary conditions 
!
      USE atmosphere, ONLY: gravx, gravz, atm
      USE grid, ONLY: fl_l
      USE dimensions
      USE eos_gas, ONLY: rgpgc
      USE gas_constants, ONLY: gmw, gammaair, gamn, rgas
      USE gas_solid_density, ONLY: rgp, rlk
      USE gas_solid_temperature, ONLY: sieg, tg, siek, tk
      USE gas_solid_velocity, ONLY: ug, vg, uk, vk
      USE grid, ONLY: zb, dz, rb, dr, r, inr, inrb
      USE grid, ONLY: nij_l, myij, data_exchange
      USE parallel, ONLY: mpime
      USE particles_constants, ONLY: rl, inrl, nsolid
      USE pressure_epsilon, ONLY: p, ep
      USE reactions, ONLY: irex
      USE set_indexes
      USE time_parameters, ONLY: dt, time
!
      IMPLICIT NONE
!
      REAL*8 :: prif, pnn2, p1nn,
     & zrif, trif, rhorif, cost, costc,
     & rm1n, rmcn, rmcnn, rmc1n, rm0n, 
     & rm1nn, rm1knn, rm2n, 
     & v2n, vg2, vcnn, vcn, v1n, vc1n, 
     & ucn, u1n, uc1n, ug2, u2n, ucnn,
     & epk, epc, epcn, ep1nn,
     & t1nn, tc1n, tcn 
      REAL*8 :: dz1, dzc, dr1, drc, indrc

      INTEGER :: i, j, ij_g, k, kg
      INTEGER :: n2, nflr, nflt, nfll, nflb
      INTEGER :: nflbr, nfltr, nfllt, nfllb
      INTEGER :: ij
      INTEGER :: mm
!
      DO ij = 1, nij_l
        ij_g = myij(0, 0, ij)
        j = ( ij_g - 1 ) / ndi + 1
        i = MOD( ( ij_g - 1 ), ndi) + 1
        IF(fl_l(ij).EQ.1) THEN
          CALL subscl(ij)
!
          nflr=fl_l(ipj)
          nfltr=fl_l(ipjp)
          nflt=fl_l(ijp)
          nfll=fl_l(imj)
          nfllb=fl_l(imjm)
          nfllt=fl_l(imjp)
          nflb=fl_l(ijm)
          nflbr=fl_l(ipjm)
!
! ***** right boundary conditions
!
          IF( nflr.NE.1 .AND. nfltr.NE.1 ) THEN

            n2 = ipj
            IF (ipj .GT. nij_l .AND. time .EQ. dt) THEN
             WRITE(8,*) 'right bdry'
       WRITE(8,*) 'warning: boundary cell not belonging to proc', mpime
             WRITE(8,*) 'ij, i, j, ij_g', ij, i, j, ij_g
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

! ... OUTFLOW ...
              IF(ug(ij).GE.0.D0) THEN
!
! ... definitions
                drc=(dr(ndi)+dr(ndi-1))*0.5D0
                indrc = 1.0D0/drc
                dr1=dr(ndi-1)
                uc1n=ug(imj)
                ucn=ug(ij)
!
! ... interpolation of mid-point values
                u1n=(ug(imj)+ug(ij))*0.5D0
                u2n=(ug(n2)+ug(ij))*0.5D0
                epcn=(ep(ij)+ep(n2))*0.5D0
                tc1n=(tg(imj)+tg(ij))*0.5D0
                tcn=(tg(ij)+tg(n2))*0.5D0
!
! ... extrapolation of radial velocity to time (n+1)dt
                ucnn=u1n*r(i)*inrb(i)
                ug(n2)=u2n*r(i+1)*inrb(i+1)
                IF(j.EQ.2) ug(ipjm)=-ug(n2)
!
! ... calculation of the gas volumetric fraction at time (n+1)dt
! ... from the mass balance equation of solids in cell (ij)
                epc=0.D0
                DO k=1,nsolid
                 epk=rlk(k,ij)*inrl(k)-dt*inrl(k)*inr(i)/dr1*
     $           (rb(i)*uk(k,ij)*rlk(k,ij)-rb(i-1)*uk(k,imj)*rlk(k,imj))
                 epc=epc+epk
                END DO
                ep1nn=1.D0-epc
!
! ... extrapolation of the temperature to time (n+1)dt
                t1nn=tc1n
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
                rm0n=rm0n+p(imj)*ep(imj)*gmw(6)/
     &               (rgas*tg(imj))
                rmcn=(rm2n+rm1n)*0.5D0
                rmc1n=(rm1n+rm0n)*0.5D0
!
                rmcnn=rmcn-dt*indrc*inrb(i)*(r(i+1)*u2n*rmcn-
     $            r(i)*u1n*rmc1n)
!
                p1nn=(1.D0/rm1nn)*(-rm1knn+rm1n-dt/dr1*inr(i)*
     $            (rb(i)*rm1n*ucn-rb(i-1)*rm0n*uc1n))
!
! ... Calculation of the advanced-time fluid pressure from Momentum balance
! ... equation of the mixture
                p(n2)=-rmcnn*ucnn+
     $            rmcn*ucn-dt*indrc*inrb(i)*(r(i+1)*u2n*ucn*rmcn-
     $            r(i)*u1n*uc1n*rmc1n)+
     $            dt*p1nn*indrc
                p(n2)=p(n2)/(dt*indrc)
!
! ... Correct non-physical pressure
                IF (p(n2).LT.0.0D0) p(n2) = p(ij)
!
! ... extrapolation of the temperature and solid fraction to time (n+1)dt
                ep(n2)=epcn
                tg(n2)=tcn
! ... INFLOW ...
              ELSE
                zrif=zb(j)+0.5D0*(dz(1)-dz(j))
                CALL atm(zrif,prif,trif)
                rhorif=prif*gmw(6)/(rgas*trif)
                cost=prif/(rhorif**gammaair)
                costc=(gammaair*cost**(1.D0/gammaair))/(gammaair-1.D0)
                tg(n2)=trif
                ep(n2)=1.D0
                ug(n2)=ug(ij)*rb(i)*inrb(i+1)
                ug2=ug(n2)*rb(i+1)*inr(i+1)
!
! ... Adiabatic inflow
                p(n2)=(prif**gamn-(ug2**2)/(2.D0*costc))**(1.D0/gamn)
!
! ... Correct non-physical pressure
                IF (p(n2).LT.0.0D0) p(n2) = p(ij)
              ENDIF
!
              rgp(n2)=p(n2)*ep(n2)*gmw(6)/(rgas*tg(n2))
!
              IF(nfltr.EQ.4) vg(n2)=vg(ij)
              DO k=1,nsolid
                IF(nfltr.EQ.4) vk(k,n2)=vk(k,ij)
                IF(uk(k,ij).GT.0.D0) THEN
                  rlk(k,n2)=rlk(k,ij)
                  uk(k,n2)=uk(k,ij)*rb(i)*inrb(i+1)
                ELSE
                  rlk(k,n2)=0.0D0
                  uk(k,n2) =0.0D0
                ENDIF
              END DO

              DO kg=1,ngas
                rgpgc(kg,n2)=rgpgc(kg,ij)
              END DO

              IF(irex.GE.0) THEN
                sieg(n2)=sieg(ij)
                DO k=1,nsolid
                  siek(k,n2)=siek(k,ij)
                  tk(k,n2)=tk(k,ij)
                END DO
              ENDIF

            CASE DEFAULT
              CONTINUE
            END SELECT

          END IF
!
! ***** top boundary conditions
!
          IF (nflt .NE. 1 .AND. nfltr .NE. 1) THEN 

            n2 = ijp
            IF (ijp .GT. nij_l .AND. time .EQ. dt) THEN
             WRITE(8,*) 'top bdry'
       WRITE(8,*) 'warning: boundary cell not belonging to proc', mpime
             WRITE(8,*) 'ij, i, j, ij_g', ij, i, j, ij_g
            END IF

            SELECT CASE (nflt)

            CASE (2)

              ug(n2)=ug(ij)
              DO k=1,nsolid
                uk(k,n2)=uk(k,ij)
              END DO

            CASE (3)
  
              ug(n2)=-ug(ij)
              DO k=1,nsolid
                IF(rlk(k,ij).GT.0.D0) uk(k,n2)=-uk(k,ij)
              END DO

            CASE (4)
!
! ... OUTFLOW ...
              IF(vg(ij).GE.0.D0) THEN
!
! ... 
!                pnn2=p(n2)
                pnn2=p(ij)
!
! ... definitions and linear interpolations
                dzc=(dz(ndj)+dz(ndj-1))*0.5D0
                dz1=dz(ndj-1)
!
                vc1n=vg(ijm)
                vcn=vg(ij)
                v1n=(vg(ijm)+vg(ij))*0.5D0
                v2n=(vg(n2)+vg(ij))*0.5D0
                epcn=(ep(ij)+ep(n2))*0.5D0
                tc1n=(tg(ijm)+tg(ij))*0.5D0
                tcn=(tg(ij)+tg(n2))*0.5D0
!
! ... Extrapolation of vertical velocities to time (n+1)dt
                vcnn=v1n
                vg(n2)=v2n
!
! ... calculation of gas volumetric fraction at time (n+1)dt
! ... from the mass balance equation of solids, in cell (ij)
                epc=0.D0
                DO k=1,nsolid
                  epk=rlk(k,ij)*inrl(k)-dt*inrl(k)/dz1*
     $            (vk(k,ij)*rlk(k,ij)-vk(k,ijm)*
     &            rlk(k,ijm))
                  epc=epc+epk
                END DO 
                ep1nn=1.D0-epc
!
! ... Extrapolation of the temperature to time (n+1)dt
! ... by using a diffusion equation
                IF (ijm .LE. nij_l) THEN
                  mm = myij( 0,-1, ijm)
                ELSE
                  mm = ijm
                END IF 
                t1nn=2.D0*tg(ijm)-tg(mm)      
!                t1nn=tg(ijm)                
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
                rm0n=rm0n+p(ijm)*ep(ijm)*gmw(6)/
     &                   (rgas*tg(ijm))
                rmcn=(rm2n+rm1n)*0.5D0
                rmc1n=(rm1n+rm0n)*0.5D0
!
                rmcnn=rmcn-dt/dzc*(v2n*rmcn-v1n*rmc1n)
!
                p1nn=(1.D0/rm1nn)*
     &            (-rm1knn+rm1n-dt/dz1*(rm1n*vcn-rm0n*vc1n))
!
! ... Calculation of the advanced-time fluid pressure from Momentum balance
! ... equation of the mixture
                p(n2)=-rmcnn*vcnn+
     $            rmcn*vcn-dt/dzc*(v2n*vcn*rmcn-v1n*vc1n*rmc1n)+
     $            dt*p1nn/dzc+gravz*dt*rmcn
                p(n2)=p(n2)/(dt/dzc)
!
! ... (Correct non-physical pressure)
                IF(p(n2).LE.0.D0) THEN
                p(n2)=pnn2
                END IF
!
                ep(n2)=epcn
                tg(n2)=2.D0*tg(ij)-tg(ijm)
              ELSE
! ... INFLOW ...
!                pnn2=p(n2)
                pnn2=p(ij)
!
                zrif=zb(ndj)+0.5D0*(dz(1)-dz(ndj))
                CALL atm(zrif,prif,trif)
                rhorif=prif*gmw(6)/(rgas*trif)
                cost=prif/(rhorif**gammaair)
                costc=(gammaair*cost**(1.D0/gammaair))/(gammaair-1.D0)
                vg(n2)=vg(ij)
                vg2=vg(ij)
                ep(n2)=1.D0
                tg(n2)=trif
                p(n2)=(prif**gamn-(vg2**2/2.D0)/costc)
     $            **(1.D0/gamn)
!
! ... (Correct non-physical pressure)
                IF(p(n2).LE.0.D0) THEN
                p(n2)=pnn2
                END IF
!
              ENDIF
!
              rgp(n2)=p(n2)*ep(n2)*gmw(6)/(rgas*tg(n2))
!
              IF(nfltr.EQ.4) ug(n2)=ug(ij)
              IF(i .EQ. (ndi-1) .AND. j .EQ. (ndj-1)) THEN
                vg(ipjp)=vg(n2)
                ug(ipjp)=ug(n2)
              ENDIF
              DO k=1,nsolid
                IF(nfltr.EQ.4) uk(k,n2)=uk(k,ij)
                IF(vk(k,ij).GT.0.D0) THEN
                  rlk(k,n2)=rlk(k,ij)
                  vk(k,n2)=vk(k,ij)
                ENDIF
                IF(i .EQ. (ndi-1) .AND. j .EQ. (ndj-1)) THEN
                  vk(k,ipjp)=vk(k,n2)
                  uk(k,ipjp)=uk(k,n2)
                ENDIF
              END DO 

              DO kg=1,ngas
                rgpgc(kg,n2)=rgpgc(kg,ij)
              END DO 
             
              IF(irex.GE.0) THEN
                sieg(n2)=sieg(ij)
                IF (i .EQ. (ndi-1) .AND. j .EQ. (ndj-1)) THEN
                  tg(ipjp)=tg(n2)
                END IF
                DO k=1,nsolid
                  siek(k,n2)=siek(k,ij)
                  tk(k,n2)=tk(k,ij)
                  IF (i .EQ. (ndi-1) .AND. j .EQ. (ndj-1)) THEN
                    tk(k,ipjp)=tk(k,n2)
                  END IF
                END DO
              ENDIF
!
            CASE DEFAULT
              CONTINUE
            END SELECT

          END IF
!
! ***** left boundary conditions
!
          IF(nfll.NE.1.AND.nfllt.NE.1) THEN
            n2 = imj
            IF (imj .GT. nij_l .AND. time .EQ. dt) THEN
             WRITE(8,*) 'left bdry'
       WRITE(8,*) 'warning: boundary cell not belonging to proc', mpime
             WRITE(8,*) 'ij, i, j, ij_g', ij, i, j, ij_g
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
                drc=(dr(i)+dr(i-1))*0.5D0
                indrc = 1.0D0/drc
                dr1=dr(i)
                uc1n=ug(ij)
                ucn=ug(n2)
!
! ... interpolation of mid-point values
                u1n=(uc1n+ucn)*0.5D0
                epcn=(ep(ij)+ep(n2))*0.5D0
                tcn=(tg(ij)+tg(n2))*0.5D0
                tc1n=(tg(ij)+tg(ipj))*0.5D0
!
! ... extrapolation of radial velocity 
                u2n=ucn-0.5*(uc1n-ucn)/dr(i)*dr(1)
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
                rm2n=rm2n+p(n2)*ep(n2)*gmw(6)/(rgas*tg(n2))
                rm1n=rm1n+p(ij)*ep(ij)*gmw(6)/(rgas*tg(ij))
                rm0n=rm0n+p(ipj)*ep(ipj)*gmw(6)/
     &               (rgas*tg(ipj))
                rmcn=(rm2n+rm1n)*0.5D0
!
! ... Particle densities
                rm1knn=0.D0
                DO k=1,nsolid
                  rm1knn=rm1knn+rlk(k,ij) 
                END DO
                rmc1n=(rm1n+rm0n)*0.5D0
!
! ... OUTFLOW ...
              IF(ug(ij).LE.0.D0) THEN
                IF(j.EQ.2) ug(imjm)=-ug(n2)
!
! ... calculation of the gas volumetric fraction at time (n+1)dt
! ... from the mass balance equation of solids in cell (ij)
                epc=0.D0
                DO k=1,nsolid
                 epk=rlk(k,ij)*inrl(k)-dt*inrl(k)*inr(i)/dr1*
     &           (rb(i)*uk(k,ij)*rlk(k,ipj) - 
     &            rb(i-1)*uk(k,imj)*rlk(k,ij))
                 epc=epc+epk
                END DO
                ep1nn=1.D0-epc
!
! ... extrapolation of the temperature to time (n+1)dt
                t1nn=tc1n
!
! ... calculation of the mixture density and fluid pressure 
! ... from a mass balance equation for the mixture.
! ... Use gas velocity for mixture.
!
! ... Mixture Density at time (n+1)dt
                rmcnn=rmcn-dt*indrc*inrb(i-1)*(r(i)*u1n*rmc1n-
     &            r(i-1)*u2n*rmcn)
!
! ... Pressure at time (n+1)dt
                rm1nn=ep1nn*gmw(6)/(rgas*t1nn)
                p1nn=(1.D0/rm1nn)*(-rm1knn+rm1n-dt/dr1*inr(i)*
     &            (rb(i)*uc1n*rm0n-rb(i-1)*rm1n*ucn))
!
! ... An arbitrary choice
                p(n2)=p1nn
!
! ... extrapolation of the temperature and solid fraction to time (n+1)dt
                ep(n2)=epcn
                tg(n2)=tcn
                DO k=1,nsolid
                  IF(nfltr.EQ.4) vk(k,n2)=vk(k,ij)
                  rlk(k,n2)=rlk(k,ij)
                  tk(k,n2)=tk(k,ij)
                  uk(k,n2)=uk(k,ij)
                END DO
! ... INFLOW ...
              ELSE
                zrif=zb(j)+0.5D0*(dz(1)-dz(j))
                CALL atm(zrif,prif,trif)
                tg(n2)=trif
                ep(n2)=1.D0
!
                p(n2)=prif
                ug(n2)=ug(ij)
                DO k=1,nsolid
                  IF(nfltr.EQ.4) vk(k,n2)=vk(k,ij)
                  rlk(k,n2)=rlk(k,ij)
                  tk(k,n2)=trif
                  uk(k,n2)=uk(k,ij)
                END DO
              ENDIF
                ucnn=
     &             rmcn*ucn - dt*indrc*
     &            inrb(i-1)*(r(i)*u1n*uc1n*rmc1n -
     &            r(i-1)*u2n*ucn*rmcn)
                ug(n2) = ucnn / rmcnn
!
              rgp(n2)=p(n2)*ep(n2)*gmw(6)/(rgas*tg(n2))
              IF(nfltr.EQ.4) vg(n2)=vg(ij)
              DO kg=1,ngas
                rgpgc(kg,n2)=rgpgc(kg,ij)
              END DO

              IF(irex.GE.0) THEN
                sieg(n2)=sieg(ij)
              ENDIF
!
            CASE DEFAULT
              CONTINUE
            END SELECT
          END IF
!
! ***** bottom boundary conditions
!
          IF(nflb.NE.1.AND.nflbr.NE.1) THEN
            n2 = ijm
            IF (ijm .GT. nij_l .AND. time .EQ. dt) THEN
             WRITE(8,*) 'bottom bdry'
       WRITE(8,*) 'warning: boundary cell not belonging to proc', mpime
             WRITE(8,*) 'ij, i, j, ij_g', ij, i, j, ij_g
            END IF
            SELECT CASE (nflb)
            CASE (2)
              ug(n2)=ug(ij)
              DO k=1,nsolid
                uk(k,n2)=uk(k,ij)
              END DO 
            CASE (3)
              ug(n2)=-ug(ij)
              DO k=1,nsolid
                IF(rlk(k,ij).GT.0.D0) uk(k,n2)=-uk(k,ij)
              END DO 
            CASE DEFAULT
              CONTINUE
            END SELECT
          END IF

        END IF
!
      END DO
!
      RETURN
      END
