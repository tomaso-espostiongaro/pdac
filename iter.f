
!
! ... MODIFICARE_X3D (fino fine file)
!
!----------------------------------------------------------------------
      MODULE iterative_solver
!----------------------------------------------------------------------
      USE grid, ONLY: dz, dr, r, indz, indr, inr
      IMPLICIT NONE
!
      REAL*8, DIMENSION(:),   ALLOCATABLE :: rgfr, rgft
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: rlfr, rlft
!
      REAL*8, DIMENSION(:), ALLOCATABLE :: conv, abeta
      REAL*8  :: omega, dg
      INTEGER :: inmax, maxout
!
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE iter
!
      USE dimensions
      USE eos_gas, ONLY: eosg, ygc, xgc, rags, cg
      USE gas_solid_temperature, ONLY: sieg, tg
      USE gas_solid_velocity, ONLY: ug, wg, us, ws
      USE eulerian_flux, ONLY: masfg, masfs
      USE gas_solid_density, ONLY: rog, rgp, rgpn, rlk, rlkn
      USE grid, ONLY: fl_l
      USE grid, ONLY: myijk, ncint, ncdom, data_exchange
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE tilde_momentum, ONLY: appu, appw
      USE parallel, ONLY: mpime
      USE particles_constants, ONLY: rl, inrl
      USE phases_matrix, ONLY: mats, matsa, velsk, velsk2
      USE pressure_epsilon, ONLY: p, ep
      USE set_indexes
      USE tilde_momentum, ONLY: rug, rwg, rus, rws
      USE time_parameters, ONLY: time, dt, timestart, tpr
      USE heat_capacity, ONLY: cp
!
      IMPLICIT NONE
!
! ... HW performance monitor include file
!#include "f_hpm.h" 
!
!
      INTEGER :: i, j , is, ij, imesh
      REAL*8 :: rlkx, rlx, dgorig, epx, rlkz
      REAL*8 :: omega0, dgz, dgx
      REAL*8  :: d3, p3

      INTEGER :: nit, loop, mustit, kros
      LOGICAL, ALLOCATABLE :: converge(:)
      LOGICAL :: first
!
      ALLOCATE(conv(ncint), abeta(ncint))
      conv = 0.0D0
      abeta = 0.0D0
!
      CALL betas(conv, abeta)
!
      ALLOCATE(converge(ncint))
!
! ... Allocate and initialize local mass fluxes.
!
      ALLOCATE( rgfr( ncdom ), rgft( ncdom ))
      ALLOCATE( rlfr( nsolid, ncdom ), rlft( nsolid, ncdom ))
      rgfr = 0.0D0
      rgft = 0.0D0
      rlfr = 0.0D0
      rlft = 0.0D0
!
! ... An approximate value for gas density fluxes is guessed, 
! ... using estimates of velocity and pressure.
! ... The best esitmates for new velocities and pressure
! ... are their values at previous time-step.
!
      CALL data_exchange(p)
!
! ... Solve the explicit phases matrix to get "new" velocities,
! ... biassed by the wrong pressure field.
!
      DO ij = 1, ncint
        imesh = myijk( ip0_jp0_kp0_, ij)
        IF(fl_l(ij).EQ.1) THEN
          CALL subscr(ij)
          CALL matsa(ij, ep(ij), ep(ijr), ep(ijt),            &
                     p(ij), p(ijr), p(ijt), rlk(:,ij),        &
                     rlk(:,ijr), rlk(:,ijt),                  &
                     rgp(ij),rgp(ijr), rgp(ijt) )   
          CALL velsk2(ug(ij), wg(ij), us(:,ij), ws(:,ij), ij)
        END IF
      END DO
! 
      CALL data_exchange(ug)
      CALL data_exchange(us)
      CALL data_exchange(wg)
      CALL data_exchange(ws)
!
! ... Put the new biassed velocities into the Gas Mass Balance 
! ... equation.
!
      DO ij = 1, ncint
        IF(fl_l(ij).EQ.1) THEN
          CALL subscr(ij)
          imesh = myijk( ip0_jp0_kp0_, ij)
          i = MOD( ( imesh - 1 ), nr) + 1
          CALL masfg(rgfr(imj), rgft(ijm), rgfr(ij),rgft(ij), &
                      rnb(ug,ij),rnb(wg,ij),nb(rgp,ij),i)
        END IF
      END DO
!
! ... Start external iterative sweep.
!
! ... The correction equation are iterated on the mesh
! ... to propagate the updated velocities and pressure
! ... on each cell to its neighbours.
!
      omega0 = omega
      mustit = 1
!
! ... Start the HW performance monitor
!      call f_hpmstart( 1, ' ITER ' )
!
      dgx = 0.0D0
      dgz = 0.0D0
!
      sor_loop: DO nit = 1, maxout
!
!         call f_hpmstart( 2, ' SWEEP ' )
!
         converge = .FALSE.
!
         mesh_loop: DO ij = 1, ncint

           imesh = myijk( ip0_jp0_kp0_, ij)
           j  = ( imesh - 1 ) / nr + 1
           i  = MOD( ( imesh - 1 ), nr) + 1
           
            CALL subscr(ij)

            first = .TRUE.

            IF( fl_l(ij) /= 1 ) converge(ij) = .TRUE.
            IF( fl_l(ij) == 1 ) THEN

              loop =  0
              kros = -1
!
! ... Compute the resudual of the Mass Balance equation of the gas phase
! ... using guessed velocities and pressure.
!
              dgx = (rgfr(ij) - rgfr(imj)) * inr(i) * indr(i)
              dgz = (rgft(ij) - rgft(ijm)) * indz(j)
              dg  = rgp(ij) - rgpn(ij) + dt * (dgx+dgz)
!                  - dt * (r1(ij)+r2(ij)+r3(ij)+r4(ij)+r5(ij))
              dgorig = dg
!
! ... If the residual is lower than the prescribed limit
! ... velocities and pressure have the correct values.
! ... Use the Mass Balance equation of the particulate phases
! ... to compute particle volumetric fractions.
!

              IF(DABS(dg).LE.conv(ij)) THEN
                converge(ij) = .TRUE.
                rlx = 0.D0
                DO is = 1, nsolid
                  rlkx = (rlfr(is,ij) - rlfr(is,imj)) * inr(i) * indr(i)
                  rlkz = (rlft(is,ij) - rlft(is,ijm)) * indz(j)
                  rlk(is, ij) = rlkn(is, ij) - dt * (rlkx + rlkz)
!                  - dt * (r1(ij)+r2(ij)+r3(ij)+r4(ij)+r5(ij))
                  IF( rlk(is,ij) .LT. 0.D0 ) rlk(is,ij) = 0.D0
                  rlx = rlx + rlk(is, ij) * inrl(is)
                END DO
                epx=rlx
                IF(epx.GT.1.D0) THEN
                  WRITE(8,*) 'warning: mass is not conserved'
                  WRITE(8,*) 'time,i,j,epst',time,i,j,epx
                  CALL error('iter', 'warning: mass is not conserved',1)
                ENDIF
                ep(ij) = 1.D0 - epx
                rgp(ij) = rog(ij) * ep(ij)

              ELSE
!
! ... If the residual is higher then the prescribed limit
! ... start the internal (in-cell) iterative sweep
! ... to correct pressure and velocities.
!
                d3=dg
                p3=p(ij)

 10             CONTINUE
!
! ... Correct the pressure at current location
! ... by using successively the Newton's method
! ... (or the secant method where Newton' fails)
! ... and the two-sided secant method.
! ... The first time in the loop, skip correction
! ... after calculation of particle velocities 
! ... and particle volumetric fractions ...
!
                IF( .NOT. ( first .AND. (nit.EQ.1) ) ) THEN
                  CALL padjust(p(ij), kros, d3, p3, omega, abeta(ij))
                END IF
                first = .FALSE.
!
! ... Use equation of state to calculate gas density 
! ... from new pressure and old temperature.
!
                CALL eosg(rags,rog(ij),cp(:,ij),cg(ij),      &    
                          tg(ij),ygc(:,ij), xgc(:,ij),       &
                          sieg(ij), p(ij), 0, 1, 0, imesh)

                rgp(ij) = ep(ij) * rog(ij)

                !call f_hpmstart( 3, ' ITER_CORE ' )
!
! ... Update gas and particles velocities using the corrected pressure
! ... at current location. Pressure at neighbour locations
! ... could still be wrong.
!
                CALL mats(ij, ep(ij), ep(ijr), ep(ijt),               &
                          ep(ijl), ep(ijb), p(ij), p(ijr),            &
                          p(ijt), p(ijl), p(ijb), rlk(:,ij),          &
                          rlk(:,ijr), rlk(:,ijt), rlk(:,ijl),         &
                          rlk(:,ijb), rgp(ij), rgp(ijr),              &
                          rgp(ijt), rgp(ijl), rgp(ijb) )

                !call f_hpmstop( 3 )

                CALL velsk(ug(ij), wg(ij), us(:,ij),                  & 
                           ws(:,ij), ug(imj), wg(ijm),                &
                           us(:,imj), ws(:,ijm), ij) 
!
! ... Put the new biassed velocities into the Particle Mass Balance
! ... equation
!
                DO is=1, nsolid
                CALL masfs(rlfr(is,imj), rlft(is,ijm), rlfr(is,ij),rlft(is,ij), &
                           rnb(us(is,:),ij),rnb(ws(is,:),ij),nb(rlk(is,:),ij),i)
                END DO
!
! ... and compute the corrected particle densities.
!
                rlx=0.D0
                DO is=1,nsolid
                  rlkx = (rlfr(is,ij)-rlfr(is,imj)) * inr(i) * indr(i)
                  rlkz = (rlft(is,ij)-rlft(is,ijm)) * indz(j)
                  rlk(is,ij) = rlkn(is,ij) - dt * (rlkx + rlkz)
!                             - dt * (r1(ij)+r2(ij)+r3(ij)+r4(ij)+r5(ij))
                  IF(rlk(is,ij).LT.0.D0) rlk(is,ij)=0.D0
                  rlx=rlx+rlk(is,ij)*inrl(is)
                END DO
                epx=rlx
                IF(epx.GT.1.D0) THEN
                  WRITE(8,*) 'warning: mass is not conserved'
                  WRITE(8,*) 'time,i,j,epst',time,i,j,epx
                  CALL error('iter', 'warning: mass is not conserved',1)
                ENDIF
!
! ... Update gas volumetric fraction and gas density.
!
                ep(ij)=1.D0-epx
                rgp(ij)=rog(ij)*ep(ij)
!
! ... Calculate the new value of the residual
! ... (at least one internal iteration must be done).
!
                IF(DABS(dg) .GT. conv(ij)) THEN
                   CALL masfg(rgfr(imj), rgft(ijm), rgfr(ij),rgft(ij),   &
                              rnb(ug,ij),rnb(wg,ij),nb(rgp,ij),i)
                   dgx = (rgfr(ij)-rgfr(imj))*inr(i)*indr(i)
                   dgz = (rgft(ij)-rgft(ijm))*indz(j)
                   dg = rgp(ij) - rgpn(ij) + dt * (dgx+dgz)
!                     - dt * (r1(ij)+r2(ij)+r3(ij)+r4(ij)+r5(ij))
!
! ... If the new residual is still greater then the prescribed value
! ... continue the internal sweep until the maximum number (inmax) 
! ... of inner iteration is reached.
!
                   IF ((DABS(dg).GT.conv(ij) .OR. DABS(dg).GE.DABS(dgorig))) THEN
                     IF(nit.EQ.1.AND.loop.EQ.0) dgorig=dg
                     d3=dg
                     loop=loop+1
! ... steepen the Newton's slope (accelerate)
                     IF(kros.LT.2.AND.loop.EQ.inmax) abeta(ij)=0.5D0*DBLE(inmax)*abeta(ij)
! ... new inner iteration
                     IF(loop .LT. inmax) THEN
                       GOTO 10
                     END IF
                   END IF
                END IF
!
              END IF
            END IF    

         END DO mesh_loop

      !call f_hpmstop( 2 )

! ... Exchanges all updated physical quantities on boundaries
! ... before starting a new sweep on the mesh.
!
        CALL data_exchange(rgp)
        CALL data_exchange(rlk)
        CALL data_exchange(p)
        CALL data_exchange(ep)
!
!********************************************************************
! ... check the convergence parameters and the history
!
      IF (MOD((time-timestart),tpr) <= dt) THEN
        IF (nit == 1) THEN
          WRITE(6,500) mpime, time
          WRITE(6,501)
        END IF
        WRITE(6,502) nit, ALL(converge)
 500    FORMAT('proc: ', I3, ' time: ', F8.3)
 501    FORMAT('iteration # :  convergence')
 502    FORMAT(I3,13X,L1)
      END IF
!*******************************************************************
!
! ... mustit equals zero when convergence is reached 
! ... simultaneously in each cell of the  subdomain.
!
        IF( ALL(converge) ) mustit = 0
!
! ... mustit must be zero in all subdomains.
!
        CALL parallel_sum_integer(mustit, 1)
!
        IF(mustit == 0) THEN
          omega = omega0
          EXIT sor_loop
        ENDIF

!
! ... If convergence is not reached in some cell
! ... start a new external sweep.
!
       IF(MOD(nit,1000).EQ.0) omega=omega*0.9D0
!
      END DO sor_loop

! ... stop the HW performance monitor
      !call f_hpmstop( 1 )
!
      IF( mustit /= 0) THEN
!********************************************************************
! ... report the number of cells where the procedure does not converge
!
        WRITE(6,700) time
        WRITE(6,*) 'convergence on proc ',mpime,' : ', ALL(converge)
        IF (.NOT.ALL(converge)) WRITE(6,*) 'cells not converged (ij, i, j): '
        DO ij = 1, ncint
          IF (.NOT.converge(ij)) THEN
            imesh = myijk( ip0_jp0_kp0_, ij)
            i  = MOD( ( imesh - 1 ), nr) + 1
            j  = ( imesh - 1 ) / nr + 1
            WRITE(6,*) ij, i, j
          END IF
        END DO
 700    FORMAT('max number of iterations reached at time: ', F8.3)
!*******************************************************************

        CALL error( ' iter ', 'max number of iters exceeded ', 1)
        omega=omega0
      END IF
!
      DEALLOCATE( rgfr, rgft)
      DEALLOCATE( rlfr, rlft)
      DEALLOCATE( rug, rwg )
      DEALLOCATE( rus, rws )
      DEALLOCATE( appu, appw)
      DEALLOCATE( conv, abeta)
!
      DEALLOCATE(converge)
!
      RETURN
      END SUBROUTINE
!
!----------------------------------------------------------------------
        SUBROUTINE padjust(p, kros, d3, p3, omega, abeta)
          REAL*8 :: p, d3, p3, omega, abeta
          INTEGER :: kros
          REAL*8, SAVE :: dp, d1, d2, p1, p2

          IF(kros.NE.3) THEN

            IF(d3.GT.0.D0) THEN
              d1=d3
              p1=p3
              IF(kros.EQ.-1) kros=0
              IF(kros.EQ.1) kros=2
            ELSE
              d2=d3
              p2=p3
              IF(kros.EQ.-1) kros=1
              IF(kros.EQ.0) kros=2
            END IF

            IF(kros.EQ.2) THEN
! ... Use secant method once when the sign of dg changes
              p = (d1*p2-d2*p1) / (d1-d2)
              abeta = (p1-p2) / (d1-d2)
              kros=3
            ELSE
! 
! ... Newton's method (iterated until the sign of dg changes)
! ... (abeta = -dp/dg) ...
              dp = -d3 * abeta
! ... with Under/Over-Relaxation ...
              dp =  dp * omega
! ... and a constrain on the maximum  correction.
              IF(-dp*dsign(1.D0,d3).GT.2.5D-1*p3) THEN
                dp=-2.5D-1*dsign(1.D0,d3)*p3
              END IF 
              p = p + dp

            END IF

          ELSE
! ... Use two-sided secant method
            p = newp(d1, d2, d3, p1, p2, p3)
            IF(d3.GT.0.D0) THEN
              d1=d3
              p1=p3
            ELSE
              d2=d3
              p2=p3
            END IF

          END IF
          p3=p
          RETURN
        END SUBROUTINE 
!----------------------------------------------------------------------
      REAL*8 FUNCTION newp(d1, d2, d3, p1, p2, p3)
!----------------------------------------------------------------------
! ... Bi-secant method
!
      USE dimensions
      IMPLICIT NONE
!
      REAL*8, INTENT(IN) :: d1, d2, d3, p1, p2, p3
      REAL*8 :: pa, pb
!
      IF (d1 .NE. d3) THEN
        pa = (d1*p3 - d3*p1) / (d1 - d3)
      END IF

      IF( d1*d3 .LE. 0.D0) THEN
        IF (d2 .NE. d3) THEN
          pb = (d2*p3 - d3*p2) / (d2 - d3)
        ELSE
          pb = 0.5D0 * (p1 + p3)
        END IF 
        IF(pb .LT. p3 .OR. pb .GT. p1) pb = 0.5D0 * (p1 + p3)
      ELSE
        IF(d1 .EQ. d3) pa = 0.5D0 * (p2 + p3)
        IF(pa .LT. p2 .OR. pa .GT. p3) pa = 0.5D0 * (p2 + p3)
        pb = (d2*p3 - d3*p2) / (d2-d3)
      END IF
!
      newp = 0.5D0 * (pa + pb)
!
      RETURN
      END FUNCTION
!----------------------------------------------------------------------
      SUBROUTINE betas(cnv, abt)
! 
      USE dimensions
      USE eos_gas, ONLY: rags
      USE gas_constants, ONLY: gammaair
      USE gas_solid_density, ONLY: rog, rgp
      USE grid, ONLY: fl_l
      USE grid, ONLY: r, rb, dr, dz
      USE grid, ONLY: ncint, myijk
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE pressure_epsilon, ONLY: p, ep
      USE set_indexes
      USE time_parameters, ONLY: dt, time
      IMPLICIT NONE
!
      REAL*8 :: rbtop, rbleft, rbright, rbbot, rbeta
      REAL*8 :: dt2r, dt2z, dzp, drm, drp, dzm
      REAL*8 :: indrp, indrm, indzp, indzm
      REAL*8 :: gam, csound
       
      INTEGER :: nflb, nflt, nfll, nflr
      INTEGER :: i, j, imesh
      INTEGER :: ij
      REAL*8, INTENT(OUT) :: cnv(:), abt(:)
!
      REAL*8, PARAMETER :: delg=1.D-8
!
      DO ij = 1, ncint
        IF (fl_l(ij) == 1) THEN
          imesh = myijk( ip0_jp0_kp0_, ij)
          CALL subscr(ij)
          j  = ( imesh - 1 ) / nr + 1
          i  = MOD( ( imesh - 1 ), nr) + 1

          drp=dr(i)+dr(i+1)
          drm=dr(i)+dr(i-1)
          dzp=dz(j)+dz(j+1)
          dzm=dz(j)+dz(j-1)
!
          indrp=1.D0/drp
          indrm=1.D0/drm
          indzp=1.D0/dzp
          indzm=1.D0/dzm
!
          nflr=fl_l(ipj)
          nfll=fl_l(imj)
          nflt=fl_l(ijp)
          nflb=fl_l(ijm)

          IF( (nflr .NE. 1) .AND. (nflr .ne. 4) ) THEN
            rbright = 0.D0
          ELSE
            rbright = ( dr(i+1)*ep(ij)+dr(i)*ep(ijr))*indrp*indrp*2.D0
          END IF

          IF( (nfll .NE. 1) .AND. (nfll .ne. 4) ) THEN
            rbleft = 0.D0
          ELSE
            rbleft = ( dr(i-1)*ep(ij)+dr(i)*ep(ijl) )*indrm*indrm*2.D0 
          END IF

          IF( (nflt .NE. 1) .AND. (nflt .ne. 4) ) THEN 
            rbtop = 0.0D0
          ELSE
            rbtop = ( dz(j+1)*ep(ij)+dz(j)*ep(ijt) )*indzp*indzp*2.D0 
          END IF

          IF( (nflb .NE. 1) .AND. (nflb .ne. 4) ) THEN
            rbbot=0.D0
          ELSE
            rbbot = ( dz(j-1)*ep(ij)+dz(j)*ep(ijb) )*indzm*indzm*2.D0
          END IF

          dt2z = dt * dt * indz(j)
          dt2r = dt * dt * inr(i) * indr(i)
!
! ... Inverse of the squared sound velocity
!
          gam = gammaair
          rags = rog(ij) / p(ij) / gam
          csound = DSQRT(1.0D0 /rags)
!
! ... rbeta = dD_g/dP
!
          rbeta = ep(ij) * rags + dt2z * (  rbtop + rbbot ) +    &
                  dt2r * (  rb(i) * rbright + rb(i-1) * rbleft )
!
          abt(ij) = 1.D0 / rbeta
          cnv(ij) = delg * rgp(ij)
        END IF
      END DO
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      END MODULE iterative_solver
!----------------------------------------------------------------------
