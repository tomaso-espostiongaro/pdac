!-----------------------------------------------------------------------
      MODULE boundary_conditions
!-----------------------------------------------------------------------
      USE atmosphere, ONLY: gravz, atm
      USE dimensions
      USE gas_constants, ONLY: gmw, gammaair, gamn, rgas
      USE gas_solid_density, ONLY: rgp, rlk
      USE gas_solid_temperature, ONLY: sieg, tg, sies, ts
      USE eos_gas, ONLY: rgpgc, xgc
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE grid, ONLY: rb, dr, r, inr, inrb
      USE grid, ONLY: zb, dz, dx, dy
      USE grid, ONLY: ncint, myijk
      USE grid, ONLY: fl_l
      USE parallel, ONLY: mpime
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: p, ep
      USE time_parameters, ONLY: dt, time

      IMPLICIT NONE

      INTEGER :: n0, n1, n2
!-----------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------
      SUBROUTINE boundary2d
!
! ... This routine computes boundary conditions 
!
      USE gas_solid_velocity, ONLY: ug, wg, us, ws
      USE set_indexes, ONLY: ipj, ipjp, imj, imjp, ijp, ijm, ipjm, imjm
      USE set_indexes, ONLY: subscr
!
      IMPLICIT NONE
!
      REAL*8:: gravx
      REAL*8 :: d1, d2
      INTEGER :: ij, i, j, imesh, ig, is
!
      DO ij = 1, ncint

        imesh = myijk( ip0_jp0_kp0_, ij)

        j = ( imesh - 1 ) / nr + 1
        i = MOD( ( imesh - 1 ), nr) + 1
        IF(fl_l(ij) == 1) THEN
          CALL subscr(ij)
!
! ***** Right boundary conditions
!
         n2 = ipj
         n1 = ij
         n0 = imj

         IF ( (fl_l( ipj ) /= 1) .AND. (fl_l( ipjp ) /= 1) ) THEN 
          SELECT CASE (fl_l( n2 )) 

          CASE (2) 

            wg(n2)=wg(ij)
            DO is=1,nsolid
              ws(n2,is)=ws(ij,is)
            END DO

          CASE (3)

            wg(n2)=-wg(ij)
            DO is=1,nsolid
              IF(rlk(ij,is).GT.0.D0) ws(n2,is)=-ws(ij,is)
            END DO  

          CASE (4)

            d1 = dr(nr-1)
            d2 = dr(nr)
            gravx = 0.D0
!
! ... Compute the normal component of the velocities and scalar fields
!
            CALL inoutflow(ug(n0),ug(n1),ug(n2),              &
                             us(n0,:),us( n1,:),us(n2,:),d1,d2,gravx,j)

! ... Compute tangential components of velocities

            IF(j == 2) ug(ipjm)=-ug(n2)
            IF( fl_l (ipjp) == 3 ) ug(ipjp) = -ug(n2)
            IF( fl_l (ipjp) == 4 ) wg(n2) = wg(ij)
            DO is=1,nsolid
              IF( fl_l (ipjp) == 4 ) ws(n2,is)=ws(ij,is)
            END DO
!
          CASE DEFAULT
            CONTINUE
          END SELECT
         END IF
!
! ****** Left boundary conditions
!
         n2 = imj

         IF ( (fl_l( imj ) /= 1) .AND. (fl_l( imjp ) /= 1) ) THEN
          SELECT CASE ( fl_l( n2 ) )
          CASE (2)
            wg(n2)=wg(ij)
            DO is=1,nsolid
              ws(n2,is)=ws(ij,is)
            END DO 
          CASE (3)
            wg(n2)=-wg(ij)
            DO is=1,nsolid
              IF(rlk(ij,is).GT.0.D0) ws(n2,is)=-ws(ij,is)
            END DO 
!
          CASE DEFAULT
            CONTINUE
          END SELECT
         END IF
!
! ****** Top boundary conditions
!
         n2 = ijp
         n1 = ij
         n0 = ijm

         IF ( (fl_l( ijp ) /= 1) .AND. (fl_l( ipjp ) /= 1) ) THEN
          SELECT CASE ( fl_l( n2 ) )

          CASE (2)

            ug(n2)=ug(ij)
            DO is=1,nsolid
              us(n2,is)=us(ij,is)
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
              IF(rlk(ij,is).GT.0.D0) us(n2,is)=-us(ij,is)
            END DO

          CASE (4)

            d1 = dx(nz-1)
            d2 = dx(nz)
!
! ... Compute the normal component of the velocities and scalar fields
!
            CALL inoutflow(wg(n0),wg(n1),wg(n2),              &
                           ws(n0,:),ws( n1,:),ws(n2,:),d1,d2,gravz,j)
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
! ... Compute tangential components of velocities
!
            IF( fl_l( ipjp ) == 4 ) ug(n2) = ug(ij)
            DO is=1,nsolid
              IF( fl_l( ipjp ) == 4 ) us(n2,is)=us(ij,is)
            END DO

          CASE DEFAULT
            CONTINUE
          END SELECT
         END IF
!
! ***** Bottom boundary conditions
!
         n2 = ijm

         IF ( (fl_l( ijm ) /= 1) .AND. (fl_l( ipjm ) /= 1) ) THEN
          SELECT CASE ( fl_l( n2 ) )

          CASE (2)

            ug(n2)=ug(ij)
            DO is=1,nsolid
              us(n2,is)=us(ij,is)
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
              IF(rlk(ij,is).GT.0.D0) us(n2,is)=-us(ij,is)
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
!
      END SUBROUTINE boundary2d
!-----------------------------------------------------------------------
      SUBROUTINE boundary3d
!
! ... This routine computes (x,y,z) boundary conditions 
!
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      USE set_indexes, ONLY: subscr
      USE set_indexes, ONLY: ipjk, imjk, ippjk, immjk, ijpk, ipjpk,    &
        imjpk, ijmk, ipjmk, imjmk, ijppk, ijmmk, ijkp, ipjkp, imjkp,   &
        ijpkp, ijmkp, ijkm, ipjkm, imjkm, ijpkm, ijmkm, ijkpp, ijkmm
!
      IMPLICIT NONE
!

      INTEGER :: ijk, i, j, k, imesh, ig, is
      REAL*8:: gravx, gravy
      REAL*8 :: d1, d2
!
      DO ijk = 1, ncint

        imesh = myijk( ip0_jp0_kp0_, ijk)

        i = MOD( MOD( imesh - 1, nx*ny ), nx ) + 1
        j = MOD( imesh - 1, nx*ny ) / nx + 1
        k = ( imesh - 1 ) / ( nx*ny ) + 1


        IF( fl_l(ijk) == 1 ) THEN

            CALL subscr(ijk)
!
! ***** East boundary conditions ***** !
!
            n2   = ipjk
	    n1   = ijk
	    n0   = imjk

            SELECT CASE ( fl_l( n2 ) ) 

            CASE (2) 

              vg( n2 ) = vg( ijk )
              wg( n2 ) = wg( ijk )
              DO is = 1, nsolid
                vs( n2, is ) = vs( ijk, is )
                ws( n2, is ) = ws( ijk, is )
              END DO

            CASE (3)

              vg( n2 ) = -vg( ijk )
              wg( n2 ) = -wg( ijk )
              DO is = 1, nsolid
                vs( n2, is ) = -vs( ijk, is )
                ws( n2, is ) = -ws( ijk, is )
              END DO  

            CASE (4)

              d1 = dx(nx-1)
	      d2 = dx(nx)
              gravx = 0.D0
!
! ... Compute the normal component of the velocities and scalar fields
!
              CALL inoutflow(ug(n0),ug(n1),ug(n2),              &
                             us(n0,:),us( n1,:),us(n2,:),d1,d2,gravx,k)

! ... Compute tangential components of velocities

! ... IF(nfltr.EQ.3) ug(ipjp) = -ug(n2)
! ... IF(nfltr.EQ.4) wg(n2) = wg(n1)
! ... IF(nfltr.EQ.4) ws(n2, is)=ws(n1, is)
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
                vs(n2,is)=vs(ijk,is)
                ws(n2,is)=ws(ijk,is)
              END DO 

            CASE (3)

              vg(n2)=-vg(ijk)
              wg(n2)=-wg(ijk)
              DO is=1,nsolid
                vs(n2,is)=-vs(ijk,is)
                ws(n2,is)=-ws(ijk,is)
              END DO 

            CASE DEFAULT

              CONTINUE

            END SELECT
!
! ***** North boundary conditions ***** !
!
!
            n2   = ijpk
	    n1   = ijk
	    n0   = ijmk

            SELECT CASE (  fl_l( n2 ) )

            CASE (2)

              ug( n2 ) = ug( ijk )
              wg( n2 ) = wg( ijk )
              DO is = 1, nsolid
                us(n2,is) = us(ijk,is)
                ws(n2,is) = ws(ijk,is)
              END DO 

            CASE (3)

              ug(n2) = -ug(ijk)
              wg(n2) = -wg(ijk)
              DO is = 1, nsolid
                 us(n2,is) = -ws(ijk,is)
                 ws(n2,is) = -ws(ijk,is)
              END DO 

            CASE (4)

              d1 = dy(ny-1)
	      d2 = dy(ny)
              gravy = 0.D0
!
! ... Compute the normal component of the velocities and scalar fields
!
              CALL inoutflow(vg(n0),vg(n1),vg(n2),              &
                             vs(n0,:),vs( n1,:),vs(n2,:),d1,d2,gravy,k)

! ... Compute tangential components of velocities

! ... IF(nfltr.EQ.3) ug(ipjp) = -ug(n2)
! ... IF(nfltr.EQ.4) wg(n2) = wg(n1)
! ... IF(nfltr.EQ.4) ws(n2, is)=ws(n1, is)
!              

            CASE DEFAULT

              CONTINUE

            END SELECT
!
! ***** South boundary conditions ***** !
!
!
            n2 = ijmk

            SELECT CASE (  fl_l( n2 ) )

            CASE (2)

              ug( n2 ) = ug( ijk )
              wg( n2 ) = wg( ijk )
              DO is = 1, nsolid
                us(n2,is) = us(ijk,is)
                ws(n2,is) = ws(ijk,is)
              END DO 

            CASE (3)

              ug(n2) = -ug(ijk)
              wg(n2) = -wg(ijk)
              DO is = 1, nsolid
                 us(n2,is) = -us(ijk,is)
                 ws(n2,is) = -ws(ijk,is)
              END DO 

            CASE DEFAULT

              CONTINUE

            END SELECT
!
! ***** Top boundary conditions ***** !
!

            n2   =  ijkp
            n1   =  ijk
            n0   =  ijkm

            SELECT CASE ( fl_l( n2 ) )

            CASE (2)

              ug(n2) = ug(ijk)
              vg(n2) = vg(ijk)
              DO is = 1, nsolid
                us(n2,is) = us(ijk,is)
                vs(n2,is) = vs(ijk,is)
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
                 us(n2,is) = -us(ijk,is)
                 vs(n2,is) = -vs(ijk,is)
              END DO

            CASE (4)

              d1 = dz(nz-1)
	      d2 = dz(nz)
!
! ... Compute the normal component of the velocities and scalar fields
!
              CALL inoutflow(wg(n0),wg(n1),wg(n2),              &
                             ws(n0,:),ws( n1,:),ws(n2,:),d1,d2,gravz,k)

! ... Compute tangential components of velocities
!
!              IF(nfltr.EQ.4) ug(n2) = ug(ijk)
!              IF(nfltr.EQ.4) us(n2,is)=us(ijk,is)
!              IF(i .EQ. (nr-1) .AND. j .EQ. (nz-1)) THEN
!                ws(ipjp,is)=ws(n2,is)
!                us(ipjp,is)=us(n2,is)
!               ENDIF
!
! ... set upper corners velocities
!
!              IF ( k == (nz-1) ) THEN
!                IF( ( i == (nx-1) ) .AND. ( j == (ny-1) ) ) THEN
!                  wg(ipjp) = wg(n2)
!                  ug(ipjp) = ug(n2)
!                ELSE IF( ( i == 2 ) .AND. ( j == (ny-1) ) ) THEN
!                ELSE IF( ( i == (nx-1) ) .AND. ( j == 2 ) ) THEN
!                ELSE IF( ( i == 2 ) .AND. ( j == 2 ) ) THEN
!                  wg(imjp) = wg(n2)
!                  ug(imjp) = ug(n2)
!                ENDIF
!              END IF
!
!
            CASE DEFAULT

              CONTINUE

            END SELECT

!
! ***** Bottom boundary conditions ***** !
!
!
            n2 = ijkm

            SELECT CASE (  fl_l( n2 ) )

            CASE (2)

              ug( n2 ) = ug( ijk )
              vg( n2 ) = vg( ijk )
              DO is = 1, nsolid
                us(n2,is) = us(ijk,is)
                vs(n2,is) = vs(ijk,is)
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
                 us(n2,is) = -us(ijk,is)
                 vs(n2,is) = -vs(ijk,is)
              END DO 

            CASE DEFAULT

              CONTINUE

            END SELECT
!
        END IF
      END DO
!
      RETURN
!
      END SUBROUTINE boundary3d
!----------------------------------------------------------------------
      SUBROUTINE inoutflow(umn, ucn, upn, usmn, uscn, uspn, d1, d2, grav, k)
!
! ... This routine computes the free in/outflow conditions in the boundary
! ... cell, i.e. the normal component of the velocity and the scalar fields

      REAL*8, INTENT(IN) :: umn, ucn
      REAL*8, INTENT(INOUT) :: upn
      REAL*8, INTENT(IN) :: usmn(:), uscn(:)
      REAL*8, INTENT(INOUT) :: uspn(:)
      REAL*8, INTENT(IN) :: d1, d2, grav
      INTEGER, INTENT(IN) :: k

      REAL*8 :: prif, pnn2, p1nn
      REAL*8 :: zrif, trif, rhorif, cost, costc
      REAL*8 :: rmcn, rmcnn, rmmn, rm1nn, rm1knn, rm2n, rm1n, rm0n
      REAL*8 :: u1n, u2n, ucnn
      REAL*8 :: eps, epc, epcn, ep1nn, epnn
      REAL*8 :: tcn, tmn, t1nn, t1n, t0n, t2n
      REAL*8 :: mg
      REAL*8 :: gravx, gravy
      REAL*8 :: dc

      INTEGER :: ig, is
!
! ... Definitions
!
      u1n = ( umn + ucn ) * 0.5D0
      u2n = ( upn + ucn ) * 0.5D0
      dc = ( d2  + d1 ) * 0.5D0
!
! ... interpolation of mid-point values

      t0n  = tg( n0 )
      t1n  = tg( n1 )
      t2n  = tg( n2 )
      tmn  = ( t0n + t1n ) * 0.5D0
      tcn  = ( t1n + t2n ) * 0.5D0
      epcn = ( ep( n1 ) + ep( n2 ) ) * 0.5D0

      IF( ucn > 0.D0 ) THEN

! ...  OUTFLOW ...
!
! ...  Extrapolations
!       ucnn = u1n
!       upn = u2n
!
! ...  Non-reflecting boundary conditions

        upn = upn - ucn*dt/d2 * (upn - ucn)
        ucnn = ucn * ( 1 - dt/d1 * (ucn - umn) )

! MODIFICARE X3D ....  IF(j == .EQ.2) ug(ipjm)=-ug(n2) ! 

        t1nn = tmn
!
! ... calculation of the gas volumetric fraction at time (n+1)dt
! ... from the mass balance equation of solids in cell (ijk)

        epc = 0.D0
        DO is = 1, nsolid
          eps = rlk(n1,is) * inrl(is) - ( dt * inrl(is) / d1 )  *          &
                ( uscn(is)*rlk(n1,is) - usmn(is)*rlk(n0,is) )
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
          rm1n=rm1n+rlk(n1,is)
          rm0n=rm0n+rlk(n0,is)
        END DO
        rm1knn = rm1n
        rm2n = rm2n + rgp(n2)
        rm1n = rm1n + rgp(n1)
        rm0n = rm0n + rgp(n0)
        rmcn=(rm2n+rm1n)*0.5D0
        rmmn=(rm1n+rm0n)*0.5D0
!
        rmcnn = rmcn-dt/dc*(u2n*rmcn-u1n*rmmn)
!
! ... Calculation of fluid pressure from the gas equation of state and 
! ... the mixture density transport equation
!
        mg=0.D0
        DO ig=1,ngas
          mg = mg + xgc(ig,n1) * gmw(ig)
        END DO
        rm1nn=ep1nn*mg/(rgas*t1nn)
                 
        p1nn = (1.D0/rm1nn) * (- rm1knn + rm1n - dt/d1*(rm1n*ucn - rm0n*umn))
        IF (epc < 1.0D-8) p1nn = p(n1)
!
! ... Calculation of the advanced-time fluid pressure from Momentum balance
! ... equation of the mixture

        p(n2) = -rmcnn*ucnn + rmcn*ucn - dt/dc*( u2n*ucn*rmcn - u1n*umn*rmmn) 
        p(n2) = p(n2) + grav * dt * rmcn
        p(n2) = p(n2)/(dt/dc) + p1nn
!
        ep(n2) = ep(n1)
        tg(n2) = tg(n1)
        DO ig=1,ngas
           rgpgc(n2,ig) = rgpgc(n2,ig) - &
                          ucn * dt/d2 * (rgpgc(n2,ig)-rgpgc(n1,ig))
        END DO
!
! ... Correct non-physical pressure
        IF (p(n2) <= 0.0D0) p(n2) = p(n1)
!
      ELSE IF ( ucn <  0.D0 ) THEN

! ... INFLOW ...
!
! ... Extrapolation

!        upn = ucn

! ... Non-reflecting b.c.
                
        upn = upn - ucn*dt/d2 * (0.D0 - u2n)
!
        zrif=zb(k)+0.5D0*(dz(1)-dz(k))
        CALL atm(zrif,prif,trif)
        rhorif=prif*gmw(6)/(rgas*trif)
        cost=prif/(rhorif**gammaair)
        costc=(gammaair*cost**(1.D0/gammaair))/(gammaair-1.D0)
!
! ... Adiabatic inflow
 
        p(n2)=(prif**gamn-(u2n**2)/(2.D0*costc))**(1.D0/gamn)
!
        ep(n2) = 1.D0
        tg(n2) = trif
        DO ig=1,ngas
          rgpgc(n2,ig) = rgpgc(n1,ig)
        END DO

! ... Correct non-physical pressure

        IF ( p(n2) < 0.0D0 ) p(n2) = prif
!
      ELSE IF ( ucn == 0.D0 ) THEN
!
        upn     = ucn
        p(n2)    = p(n1)
        ep(n2)   = ep(n1)
        tg(n2)   = tg(n1)
        DO ig=1,ngas
          rgpgc(n2,ig) = rgpgc(n1,ig)
        END DO

      ENDIF
!
! ... Set primary variables
!
      DO is = 1, nsolid
        IF( uscn(is) >= 0.D0 ) THEN
          rlk(n2,is) = rlk(n1,is)
          uspn(is)  = uscn(is)
        ELSE
          rlk(n2,is) = 0.0D0
          uspn(is)  = 0.0D0
        ENDIF
      END DO
 
      sieg(n2) = sieg(n1)
      DO is = 1, nsolid
        sies(n2,is) = sies(n1,is)
      END DO
!                
      END SUBROUTINE inoutflow
!-----------------------------------------------------------------------
      END MODULE boundary_conditions
!-----------------------------------------------------------------------
