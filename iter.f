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
      USE gas_solid_velocity, ONLY: ug, vg, uk, vk
      USE eulerian_flux, ONLY: masfg, masfk
      USE gas_solid_density, ONLY: rog, rgp, rgpn, rlk, rlkn
      USE grid, ONLY: fl_l
      USE grid, ONLY: myij, nij_l, nijx_l, data_exchange
      USE tilde_momentum, ONLY: appu, appv
      USE particles_constants, ONLY: rl, inrl, nsolid
      USE phases_matrix, ONLY: mats, matsa, velsk, velsk2
      USE pressure_epsilon, ONLY: p, ep
      USE set_indexes
      USE tilde_momentum, ONLY: rug, rvg, ruk, rvk
      USE time_parameters, ONLY: time, dt, tstop, itd, timestart
      USE th_capacity, ONLY: cp
!
      IMPLICIT NONE
!
! ... HW performance monitor include file
!#include "f_hpm.h" 
!
!
      INTEGER :: i, j , k, ij, ij_g
      REAL*8 :: rlkx, rlx, dgorig, epx, rlkz
      REAL*8 :: omega0, dgz, dgx
      REAL*8  :: d3, p3
      REAL*8 :: t1, t2
      REAL*8 :: avlp

      INTEGER :: nit, loop, mustit, kros
      REAL*8, ALLOCATABLE :: avloop(:)
      LOGICAL :: first, cvg
!
      ALLOCATE(conv(nij_l), abeta(nij_l))
      conv = 0.0D0
      abeta = 0.0D0
!
      DO ij = 1, nij_l
        IF(fl_l(ij).EQ.1) THEN
          CALL betas(conv(ij), abeta(ij), ij)
        END IF
      END DO
!
! ... avloop is the per-cycle-average of inner loops in each cells
! ... avlp is the average of avloop over the proc subdomain
!
      ALLOCATE(avloop(nij_l))
      avloop = 0
      avlp = 0
!
! ... Allocate and initialize local mass fluxes.
!
      ALLOCATE( rgfr( nijx_l ), rgft( nijx_l ))
      ALLOCATE( rlfr( ncl, nijx_l ), rlft( ncl, nijx_l ))
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
      DO ij = 1, nij_l
        ij_g = myij(0, 0, ij)
        IF(fl_l(ij).EQ.1) THEN
          CALL subscl(ij)
          CALL matsa(ij, ep(ij), ep(ijr), ep(ijt),            &
                     p(ij), p(ijr), p(ijt), rlk(:,ij),        &
                     rlk(:,ijr), rlk(:,ijt),                  &
                     rgp(ij),rgp(ijr), rgp(ijt) )   
          CALL velsk2(ug(ij), vg(ij), uk(:,ij), vk(:,ij), ij)
        END IF
      END DO
! 
      CALL data_exchange(ug)
      CALL data_exchange(uk)
      CALL data_exchange(vg)
      CALL data_exchange(vk)
!
! ... Put the new biassed velocities into the Gas Mass Balance 
! ... equation.
!
      DO ij = 1, nij_l
        IF(fl_l(ij).EQ.1) THEN
          CALL subscl(ij)
          ij_g = myij(0, 0, ij)
          i = MOD( ( ij_g - 1 ), ndi) + 1
          CALL masfg(rgfr(imj), rgft(ijm), rgfr(ij),rgft(ij), &
                      rnb(ug,ij),rnb(vg,ij),nb(rgp,ij),i)
        END IF
      END DO
!
! ... Start external iterative sweep.
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
      iterative_loop: DO nit = 1, maxout
!
!         call f_hpmstart( 2, ' SWEEP ' )
!
         DO ij = 1, nij_l
           avloop(ij) = avloop(ij) + 1

           ij_g = myij(0, 0, ij)
           j  = ( ij_g - 1 ) / ndi + 1
           i  = MOD( ( ij_g - 1 ), ndi) + 1
           
            CALL subscl(ij)

            first = .TRUE.
            IF( fl_l(ij) .EQ. 1 ) THEN

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
                rlx = 0.D0
                DO k = 1, nsolid
                  rlkx = (rlfr(k,ij) - rlfr(k,imj)) * inr(i) * indr(i)
                  rlkz = (rlft(k,ij) - rlft(k,ijm)) * indz(j)
                  rlk(k, ij) = rlkn(k, ij) - dt * (rlkx + rlkz)
!                  - dt * (r1(ij)+r2(ij)+r3(ij)+r4(ij)+r5(ij))
                  IF( rlk(k,ij) .LT. 0.D0 ) rlk(k,ij) = 0.D0
                  rlx = rlx + rlk(k, ij) * inrl(k)
                END DO
                epx=rlx
                IF(epx.GT.1.D0) THEN
                  WRITE(8,*) 'warning: mass is not conserved'
                  WRITE(8,*) 'time,i,j,epst',time,i,j,epx
                ENDIF
                ep(ij) = 1.D0 - epx
                rgp(ij) = rog(ij) * ep(ij)

              ELSE
!
! ... If the residual is higher then the prescribed limit
! ... start the internal (in-cell) iterative sweep
! ... to correct pressure and velocities.
!
                mustit=1
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
                          sieg(ij), p(ij), 0, 1, 0, ij_g)

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

                CALL velsk(ug(ij), vg(ij), uk(:,ij),                  & 
                           vk(:,ij), ug(imj), vg(ijm),                &
                           uk(:,imj), vk(:,ijm), ij) 
!
! ... Put the new biassed velocities into the Particle Mass Balance
! ... equation
!
                DO k=1, nsolid
                CALL masfk(rlfr(k,imj), rlft(k,ijm), rlfr(k,ij),rlft(k,ij), &
                           rnb(uk(k,:),ij),rnb(vk(k,:),ij),nb(rlk(k,:),ij),i)
                END DO
!
! ... and compute the corrected particle densities.
!
                rlx=0.D0
                DO k=1,nsolid
                  rlkx = (rlfr(k,ij)-rlfr(k,imj)) * inr(i) * indr(i)
                  rlkz = (rlft(k,ij)-rlft(k,ijm)) * indz(j)
                  rlk(k,ij) = rlkn(k,ij) - dt * (rlkx + rlkz)
!                             - dt * (r1(ij)+r2(ij)+r3(ij)+r4(ij)+r5(ij))
                  IF(rlk(k,ij).LT.0.D0) rlk(k,ij)=0.D0
                  rlx=rlx+rlk(k,ij)*inrl(k)
                END DO
                epx=rlx
                IF(epx.GT.1.D0) THEN
                  WRITE(8,*) 'warning: mass is not conserved'
                  WRITE(8,*) 'time,i,j,epst',time,i,j,epx
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
                     rnb(ug,ij),rnb(vg,ij),nb(rgp,ij),i)
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
                       avloop(ij) = avloop(ij) + 1
                       GOTO 10
                     END IF
                   END IF
                END IF
!
              END IF
            END IF    
         END DO
         avlp = SUM(avloop)
!
      !call f_hpmstop( 2 )

! ... Exchanges all updated physical quantities on boundaries
! ... before starting a new sweep on the mesh.
!
        CALL data_exchange(rgp)
        CALL data_exchange(rlk)
        CALL data_exchange(p)
        CALL data_exchange(ep)
!  
! ... mustit equals zero when convergence is reached (cvg = .TRUE.)
! ... simultaneously in the whole computational domain.
! ... mustit must be zero for each processor.
!
        CALL parallel_sum_integer(mustit, 1)
        IF(mustit .EQ. 0) THEN
          omega = omega0
          cvg = .TRUE.                         
          EXIT iterative_loop
        ENDIF

        mustit=0
!
! ... If convergence is not reached in some cell (cvg = .FALSE.),
! ... start a new external sweep.
!
       IF(MOD(nit,1000).EQ.0) omega=omega*0.9D0
       cvg = .FALSE.                            
!
      END DO iterative_loop
                         WRITE(*,*) nit, loop, p(2933), dg

! ... stop the HW performance monitor
      !call f_hpmstop( 1 )
!********************************************************************
! ... check convergence parameters 
! ... ( number of iteration, averaged number of internal 
! ...   iterations, averaged number of external iterations ) 
!
      avloop(:) = avloop(:)/nit
      avlp = avlp/nit/nij_l
!
      t1 = tstop - 100 * dt     
!      t1 = tstop
      t2 = tstop               
      IF (time .GE. t1 .AND. time .LE. t2) THEN
          WRITE(7,500) time
          WRITE(7,*) 'nit =', nit
          WRITE(7,'(A15,F8.4)') 'avlp =', avlp
 500      FORMAT('time =', F8.3)
      END IF
      DO ij=1,nij_l
           ij_g = myij(0, 0, ij)
           j  = ( ij_g - 1 ) / ndi + 1
           i  = MOD( ( ij_g - 1 ), ndi) + 1
        IF (avloop(ij) .GT. 5*avlp) THEN
          WRITE(7,*) 'Averaged nb. of inner loop', avlp
          WRITE(7,*) 'in (i,j)=',i,j,' avloop=',avloop(ij)
        END IF
      END DO
!*******************************************************************
!
      IF(.NOT. cvg) THEN
        WRITE(8,*) 'warning: max number outer iters exceeded'
        omega=omega0
        STOP
      END IF
!
      DEALLOCATE( rgfr, rgft)
      DEALLOCATE( rlfr, rlft)
      DEALLOCATE( rug, rvg )
      DEALLOCATE( ruk, rvk )
      DEALLOCATE( appu, appv)
      DEALLOCATE( conv, abeta)
!
      DEALLOCATE(avloop)
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
      SUBROUTINE betas(cnv, abt, ij)
! 
      USE dimensions
      USE eos_gas, ONLY: rags
      USE gas_constants, ONLY: gammaair
      USE gas_solid_density, ONLY: rog, rgp
      USE grid, ONLY: fl_l
      USE grid, ONLY: r, rb, dr, dz
      USE grid, ONLY: nij_l, myij
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
      INTEGER :: i, j, ij_g
      INTEGER, INTENT(IN) :: ij
      REAL*8, INTENT(OUT) :: cnv, abt
!
      REAL*8, PARAMETER :: delg=1.D-8
!
          ij_g = myij(0, 0, ij)
          CALL subscl(ij)
          j  = ( ij_g - 1 ) / ndi + 1
          i  = MOD( ( ij_g - 1 ), ndi) + 1

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
           abt = csound
!          abt = 1.D0 / rbeta
          cnv = delg * rgp(ij)
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      END MODULE iterative_solver
!----------------------------------------------------------------------
