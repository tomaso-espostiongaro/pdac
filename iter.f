!----------------------------------------------------------------------
      MODULE iterative_solver
!----------------------------------------------------------------------
      USE grid, ONLY: dx, dy, dz, indx, indy, indz, inr
      USE set_indexes, ONLY: stencil

      IMPLICIT NONE

      SAVE
!
! ... convective mass fluxes
!
      REAL*8, DIMENSION(:),   ALLOCATABLE :: rgfe, rgfn, rgft
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: rsfe, rsfn, rsft

      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: amats

      REAL*8  :: omega, dg
      INTEGER :: inmax, maxout

      INTEGER, ALLOCATABLE, DIMENSION(:) :: b_e, b_w, b_t, b_b, b_n, b_s

      TYPE(stencil) :: u, v, w, dens         

      PRIVATE :: u, v, w, dens
!
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
!
      SUBROUTINE iter
!
!----------------------------------------------------------------------
! ... This routine is the core of the iterative procedure to solve the
! ... mass-momentum and the interphase coupling
! ... (2D-3D-Compliant)
!
      USE control_flags, ONLY: job_type, lpr
      USE control_flags, ONLY: implicit_fluxes, implicit_enthalpy
      USE dimensions
      USE domain_decomposition, ONLY: ncint, ncdom, data_exchange
      USE domain_decomposition, ONLY: myijk, meshinds
      USE enthalpy_matrix, ONLY: ftem
      USE environment, ONLY: cpclock, timing
      USE eos_gas, ONLY: xgc, cg
      USE gas_solid_temperature, ONLY: tg, ts, sieg, sies
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      USE gas_solid_density, ONLY: rog, rgp, rgpn, rlk, rlkn
      USE grid, ONLY: flag
      USE immersed_boundaries, ONLY: immb
      USE indijk_module
      USE parallel, ONLY: mpime
      USE particles_constants, ONLY: rl, inrl
      USE phases_matrix, ONLY: assemble_matrix, assemble_all_matrix
      USE phases_matrix, ONLY: solve_velocities, solve_all_velocities
      USE phases_matrix, ONLY: matspre_3phase
      USE pressure_epsilon, ONLY: p, ep
      USE set_indexes, ONLY: imjk, ijmk, ijkm, ijkn, &
                             ijke, ijkw, ijkt, ijkb
      USE set_indexes, ONLY: subscr, first_subscr
      USE tilde_energy, ONLY: htilde
      USE tilde_momentum, ONLY: tilde, appu, appv, appw
      USE time_parameters, ONLY: time, dt, timestart, tpr
!
      IMPLICIT NONE
!
      REAL*8, DIMENSION(:), ALLOCATABLE :: conv
      REAL*8, DIMENSION(:), ALLOCATABLE :: abeta
!
! ... HW performance monitor include file
!
      INTEGER :: is, imesh, ii, jj
      REAL*8 :: dgorig
      REAL*8 :: rls
      REAL*8 :: rlkx, rlky, rlkz, rlk_tmp
      REAL*8 :: omega0
      REAL*8 :: d3, p3

! IF (TIMING)
      INTEGER :: itim
      ! INTEGER :: st0,st1,ratc
      REAL*8 :: st0,st1, cclock_wall
      EXTERNAL cclock_wall
      REAL*8 ::  timiter
      REAL*8 ::  timconv(maxout)   
!
      INTEGER :: nit, nloop, mustit
      INTEGER :: n1, n2
      REAL*8  :: avloop
      INTEGER :: i, j, k, ijk

      LOGICAL, ALLOCATABLE :: converge(:)
      LOGICAL :: cvg
      LOGICAL :: optimization
!
      ALLOCATE(conv(ncint), abeta(ncint))
      conv = 0.0D0
      abeta = 0.0D0
!
      ALLOCATE( converge( ncint ) )
!
! ... Allocate and initialize local mass fluxes.
!
      ALLOCATE( rgfe( ncdom ), rgft( ncdom ))
      ALLOCATE( rsfe( ncdom, nsolid ), rsft( ncdom, nsolid ))
      rgfe = 0.0D0
      rgft = 0.0D0
      rsfe = 0.0D0
      rsft = 0.0D0
      IF (job_type == '3D') THEN
        ALLOCATE( rgfn( ncdom ))
        ALLOCATE( rsfn( ncdom, nsolid ))
        rgfn = 0.0D0
        rsfn = 0.0D0
      END IF
!
! ... Assemble and solve the explicit phase matrix 
! ... to update velocity fields. New velocities are
! ... biassed by the wrong pressure field. The best
! ... estimate for the pressure field is given by
! ... the pressure at previous time-step.
!
      DO ijk = 1, ncint
        IF ( flag( ijk ) == 1) THEN
          CALL first_subscr( ijk )
          CALL assemble_matrix( ijk )
          CALL solve_velocities( ijk)
        END IF
      END DO
! 
! ... Exchange the updated velocities
!
      CALL data_exchange(ug)
      CALL data_exchange(us)
      CALL data_exchange(wg)
      CALL data_exchange(ws)
      IF (job_type == '3D') THEN
        CALL data_exchange(vg)
        CALL data_exchange(vs)
      END IF
!
! ... In 3D use optimized routines with two particle classes
!
      optimization = ( nsolid == 2 ) .AND. ( job_type == '3D' ) .AND. &
                     ( immb < 1 )
      IF( optimization ) THEN
        ALLOCATE( amats( 6, 6, ncint ) )
      END IF
!
! ... compute the fraction of partially filled boundary cells
!
      CALL fill_cells
!
! ... Guess an approximate value for gas density fluxes
! ... using the estimates of velocity and pressure.
!
      DO ijk = 1, ncint
        IF ( flag( ijk ) == 1 ) THEN
          CALL first_subscr( ijk )
          CALL calc_gas_mass_flux( ijk )
          !
          ! ... compute the derivative of the gas mass residual
          ! ... with respect to gas pressure
          CALL betas( conv( ijk ), abeta( ijk ), ijk )

          IF(optimization) THEN
            CALL meshinds(ijk,imesh,i,j,k)
            CALL matspre_3phase(amats(:,1,ijk),amats(:,2,ijk),amats(:,3,ijk),&
                                amats(:,4,ijk),amats(:,5,ijk),amats(:,6,ijk),&
                                i, j, k, ijk )
          END IF
        END IF
      END DO
!
! ... Here Start the external iterative sweep.
!/////////////////////////////////////////////////////////////////////
!
! ... The correction equation are iterated on the mesh
! ... to propagate the updated velocities and pressure
! ... on each cell to its neighbours.
!
      omega0 = omega
      mustit = 1
!
      timconv = 0.0d0
!
      sor_loop: DO nit = 1, maxout
!
         ! ... Compute fluxes and forces in the momentum equation
         ! ... (for fully implicit solution)
         !
         IF (implicit_fluxes) THEN

           CALL data_exchange(ug)
           IF (job_type == '3D') CALL data_exchange(vg)
           CALL data_exchange(wg)
           CALL data_exchange(us)
           IF (job_type == '3D') CALL data_exchange(vs)
           CALL data_exchange(ws)

           CALL tilde

         END IF

         ! TIMING (Convergence in the ijk cell)
!
         IF( timing ) THEN
            st0 = cpclock()
            !st0 = cclock_wall()
            !call system_clock (st0,ratc)
            !CALL f_hpmstart( 1, ' sor ' )
         END IF  

         nloop  = 0
         n1 = 0
         n2 = 0

         mesh_loop: DO ijk = 1, ncint

           converge(ijk) = .FALSE.

           IF( flag(ijk) /= 1 ) THEN
             
             converge(ijk) = .TRUE.
           
           ELSE IF( flag(ijk) == 1 ) THEN

             CALL meshinds(ijk,imesh,i,j,k)
             CALL first_subscr(ijk)

             ! ... Compute locally the residual 'dg' of the
             ! ... mass balance equation of the gas phase
             !
             CALL calc_res(i,j,k,ijk,dg)

             dgorig = dg







             IF( ABS( dg ) <= conv( ijk ) ) THEN

               n1 = n1 + 1
              
               ! ... If the residual is lower then the prescribed limit
               ! ... compute the gas and particle densities and proceed 
               ! ... to next cell
               !
               converge( ijk ) = .TRUE.

               CALL calc_eps( i, j, k, ijk  )
               rgp( ijk ) = ep( ijk ) * rog( ijk )

             ELSE IF ( ABS(dg) > conv( ijk ) ) THEN

               n2 = n2 + 1

               ! ... If the residual is higher then the prescribed limit
               ! ... start the inner (in-cell) iterative loop
               ! ... to correct pressure and velocities.
               !
               d3 = dg
               p3 = p( ijk )

               IF (optimization) THEN
                 CALL opt3_inner_loop(i, j, k, ijk, nit, d3, p3, &
                                  abeta(ijk), conv(ijk), &
                                  dgorig, nloop, cvg)
               ELSE
                 CALL inner_loop(i, j, k, ijk, nit, d3, p3, &
                                 abeta(ijk),conv(ijk), &
                                 dgorig, nloop, cvg)
               END IF

             END IF

           END IF
         END DO mesh_loop

         IF( timing ) THEN
             st1 = cpclock()
             !call system_clock(st1,ratc)
             !CALL f_hpmstop( 1 )
         END IF

         timconv(nit) = timconv(nit) + real(st1-st0)/1000.D0         
!
!*******************************************************************
! ... For each iteration on the mesh, write out:
! --> n1 : the number of cells converged
! --> n2 : the number of cells not converged
! --> avloop : the averaged number of inner iterations per cells
! --> timconv: the time for the whole mesh sweep
!
         IF (lpr > 1) THEN
           IF( n2 > 0 ) THEN
             avloop = REAL(nloop) / REAL(n2)
           ELSE
             avloop = REAL(nloop)
           END IF
           WRITE(6, fmt="( I10, 2X, F4.2, 2X, I10, 2X, F10.3)" ) &
           n2, avloop, n1, timconv(nit)
         END IF
!*******************************************************************

! ... Exchange all updated physical quantities on boundaries
! ... before starting a new sweep on the mesh.
!
        CALL data_exchange(rgp)
        CALL data_exchange(rlk)
        CALL data_exchange(p)
        CALL data_exchange(ep)
!
! ... Here the enthalpy equations can be solved implicitly
!
        IF (implicit_enthalpy) THEN
          !
          ! ... enthalpy convective and diffusive fluxes
          IF (implicit_fluxes) THEN

            CALL data_exchange(sieg)
            CALL data_exchange(sies)
            CALL data_exchange(tg)
            CALL data_exchange(ts)
            
            CALL htilde

          END IF
          !
          ! ... solve pressure coupling
          CALL ftem

        END IF
!
! ... mustit equals zero when convergence is reached 
! ... simultaneously in each cell of the subdomain
!
        IF( ALL( converge ) ) mustit = 0
!
! ... mustit must be zero simultaneously in all subdomains
!
        CALL parallel_sum_integer(mustit, 1)
!
        IF( mustit == 0 ) THEN
!*******************************************************************
! ... write out the final number of iterations
!
          IF (lpr > 1) THEN
            WRITE(6,277) nit
 277        FORMAT('number of iterations: nit = ', I4)
          END IF
!*******************************************************************
          omega = omega0
          EXIT sor_loop
        ENDIF
!
! ... If convergence is not reached in some cell
! ... start a new external sweep.
!
       IF( MOD( nit, 1000 ) == 0 ) THEN
         omega = omega * 0.9D0
         IF (lpr > 1) THEN
           WRITE(6, fmt="('  reducing relaxation parameter omega')")
           WRITE(6, fmt="('  new value = ',F12.4)") omega
         END IF
       END IF
       IF (lpr > 1) call myflush( 6 )
!
      END DO sor_loop
!
! ... End the iterative sweep
!/////////////////////////////////////////////////////////////////////
!
! IF (TIMING)
      timiter = 0.0d0
      DO itim = 1, nit
        timiter = timiter + timconv(itim)
      END DO
      WRITE(6,280)'Time for iterative solver: ',timiter
 280  FORMAT(2X,A27,F8.3) 

!
!********************************************************************
! ... If the iterative sweep concluded without convergence
! ... report the number of cells where the procedure does not converge
!
      IF( mustit /= 0 ) THEN
!
        IF (lpr >= 1) THEN
          WRITE(6,700) nit, (time+dt)
          WRITE(6,*) 'convergence on proc ',mpime,' : ', ALL(converge)
          IF (.NOT.ALL(converge))   &
                          WRITE(6,*) 'cells not converged (imesh,i,j,k): '
          DO ijk = 1, ncint
            IF ( .NOT. converge( ijk ) ) THEN
              CALL meshinds( ijk , imesh, i , j , k )
              WRITE(6,*) imesh, i , j , k
            END IF
          END DO
 700      FORMAT('max number of iterations (',I5,') reached at time: ', F8.3)
 701      FORMAT(10(F14.6))
        END IF
!
        ! ... CRASH! ...
        !
        !CALL error( ' iter ', 'max number of iters exceeded ', 1)
        omega = omega0

      END IF
!*******************************************************************
!
      DEALLOCATE(rgfe)
      DEALLOCATE(rgft)
      DEALLOCATE(rsfe)
      DEALLOCATE(rsft)

      IF (job_type == '3D') THEN
        DEALLOCATE( rgfn )
        DEALLOCATE( rsfn )
      END IF

      DEALLOCATE(conv)
      DEALLOCATE(abeta)
      DEALLOCATE(converge)
      IF( ALLOCATED( amats ) )  DEALLOCATE( amats )
      DEALLOCATE(b_e, b_w, b_t, b_b)
      IF (job_type == '3D') DEALLOCATE(b_n, b_s)
!
      RETURN
      END SUBROUTINE iter
!----------------------------------------------------------------------
!
      SUBROUTINE inner_loop( i, j, k, ijk, nit, d3, p3, abeta, conv,  &
                             dgorig, nloop, cvg )
!
!----------------------------------------------------------------------
! ... iteratively correct the pressure field in a cell by minimizing
! ... the gas mass residual
!
      USE dimensions
      USE pressure_epsilon, ONLY: p, ep
      USE gas_solid_density, ONLY: rog, rgp, rgpn, rlk, rlkn
      USE gas_solid_temperature, ONLY: tg
      USE eos_gas, ONLY: thermal_eosg, xgc
      USE phases_matrix, ONLY: assemble_all_matrix
      USE phases_matrix, ONLY: solve_all_velocities
      USE set_indexes, ONLY: subscr, imjk, ijmk, ijkm
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      USE control_flags, ONLY: job_type
      USE time_parameters, ONLY: dt

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nit, ijk, i, j, k
      REAL*8, INTENT(IN)  :: conv
      REAL*8, INTENT(INOUT)  :: d3, p3, abeta, dgorig
      LOGICAL, INTENT(OUT) :: cvg

      REAL*8 :: resx, resy, resz

      INTEGER :: loop, kros, nloop

      kros = -1

      inloop: DO loop = 1, inmax
!
        ! ... Correct the pressure at current cell
        ! ... by using successively the Newton's method
        ! ... (or the secant method where Newton is not
        ! ... applicable) and the two-sided secant method.
        ! ... The first time in a cell, skip pressure correction
        ! ... to compute once the particle velocities 
        ! ... and volumetric fractions ...
!
        IF( (loop > 1) .OR. (nit > 1) ) THEN
          CALL padjust(p(ijk), kros, d3, p3, omega, abeta )
        END IF
!
        ! ... Use equation of state to calculate gas density 
        ! ... from (new) pressure and (old) temperature.
        ! ... (notice that the explicit solution of the enthalpy
        ! ... equations affects only this updating)
!
        CALL thermal_eosg(rog(ijk),tg(ijk),p(ijk),xgc(:,ijk))
!
        rgp(ijk) = ep(ijk) * rog(ijk)

        ! ... Update gas and particles velocities using the 
        ! ... corrected pressure at current location. 
        ! ... Pressure at neighbour cells could still be wrong.

        CALL assemble_all_matrix(ijk)
        CALL solve_all_velocities(ijk)
!
        ! ... update particle and gas densities and the 
        ! ... gas mass residual 'dg'
        ! ... update the particle volumetric fractions and the void fraction

        ! ... WARNING! Only the local values of the stencil have to be updated!
        ! ... This procedure can be optimized!
        !
        CALL calc_part_mass_flux( ijk )
        CALL calc_eps( i, j, k, ijk )

        rgp( ijk ) = ep( ijk ) * rog( ijk )

        ! ... update residual of the Mass Balance equation of the gas phase
        ! ... using guessed velocities and pressure.

        ! ... WARNING! Only the local values of the stencil have to be updated!
        ! ... This procedure can be optimized!
        !
        CALL calc_gas_mass_flux( ijk )
        CALL calc_res( i, j, k, ijk, dg)
!
        IF ( ( DABS(dg) > conv .OR. DABS(dg) >= DABS(dgorig) ) ) THEN

          cvg = .FALSE.
          IF( nit == 1 .AND. loop == 1 ) dgorig = dg
          d3 = dg
          ! ... steepen the Newton's slope (accelerate)
          IF( kros < 2 .AND. loop == inmax ) abeta = 0.5D0 * inmax * abeta

        ELSE IF ( DABS(dg) <= conv ) THEN

          cvg = .TRUE.
          EXIT inloop

        END IF

      END DO inloop

      nloop = nloop + loop

      RETURN
      END SUBROUTINE inner_loop
!----------------------------------------------------------------------
!
      SUBROUTINE calc_gas_mass_flux( ijk )
!
!----------------------------------------------------------------------
! ... compute the gas mass fluxes by using the stencil, as defined in the
! ... 'subscr' module
!
      USE convective_mass_fluxes, ONLY: fmas, masf
      USE dimensions
      USE flux_limiters, ONLY: muscl
      USE gas_solid_velocity, ONLY: ug, vg, wg
      USE gas_solid_density, ONLY: rgp
      USE set_indexes, ONLY: first_nb, first_rnb, third_nb, third_rnb
      USE set_indexes, ONLY: imjk, ijmk, ijkm
      USE set_indexes, ONLY: third_subscr
      USE control_flags, ONLY: job_type, implicit_fluxes, &
                               implicit_enthalpy
!
      IMPLICIT NONE

      INTEGER :: ijk, info

        IF (job_type == '2D') THEN
!
! ... Compute Fluxes by using First Order Upwind method
!
          ! ... assemble the first order computational stencils
          CALL first_nb ( dens, rgp, ijk )
          CALL first_rnb( u, ug, ijk )
          CALL first_rnb( w, wg, ijk )

          CALL masf( rgfe( ijk ), rgft( ijk ), rgfe( imjk ), rgft( ijkm ), &
                    dens, u, w, ijk )

          IF (muscl > 0) THEN
!
! ... Second order MUSCL correction
!
            ! ... add third neighbours to computational stencils
            CALL third_subscr( ijk )
            CALL third_nb( dens, rgp, ijk)
            CALL fmas( rgfe( ijk ), rgft( ijk ), rgfe( imjk ), rgft( ijkm ), &
                       dens, u, w, ijk )
          END IF


        ELSE IF (job_type == '3D') THEN
!
! ... Compute Fluxes by using First Order Upwind method
!
          ! ... assemble the first order computational stencils
          CALL first_nb ( dens, rgp, ijk )
          CALL first_rnb( u, ug, ijk )
          CALL first_rnb( v, vg, ijk )
          CALL first_rnb( w, wg, ijk )

          CALL masf( rgfe( ijk  ), rgfn( ijk  ), rgft( ijk  ),  &
                     rgfe( imjk ), rgfn( ijmk ), rgft( ijkm ),  &
                     dens, u, v, w, ijk )
            
          IF (muscl > 0) THEN
!
! ... Second order MUSCL correction
!
            ! ... add third neighbours to computational stencils
            CALL third_subscr( ijk )
            CALL third_nb( dens, rgp, ijk)
            CALL fmas( rgfe( ijk  ), rgfn( ijk  ), rgft( ijk  ),  &
                       rgfe( imjk ), rgfn( ijmk ), rgft( ijkm ),  &
                       dens, u, v, w, ijk )

          END IF
        ENDIF

      END SUBROUTINE calc_gas_mass_flux   
!----------------------------------------------------------------------
!
      SUBROUTINE calc_res( i, j, k, ijk, res )
!
!----------------------------------------------------------------------
! ... Compute the residual of the Mass Balance equation of the gas phase
!
      USE set_indexes, ONLY: subscr, imjk, ijmk, ijkm
      USE gas_solid_density, ONLY: rog, rgp, rgpn, rlk, rlkn
      USE grid, ONLY: flag
      USE control_flags, ONLY: job_type
      USE time_parameters, ONLY: dt

      IMPLICIT NONE
      REAL*8 :: resx, resy, resz
      REAL*8, INTENT(OUT) :: res
      INTEGER, INTENT(IN) :: i, j, k, ijk
      REAL*8 :: vf, ivf
!      
! ... volume fraction not occupied by the topography
!
      vf = b_e(ijk) + b_w(ijk) + b_t(ijk) + b_b(ijk)
!
      resx = ( b_e(ijk)*rgfe(ijk) - b_w(ijk)*rgfe(imjk) ) * indx(i) * inr(i)
      resz = ( b_t(ijk)*rgft(ijk) - b_b(ijk)*rgft(ijkm) ) * indz(k)
!
      IF (job_type == '2D') THEN
        resy = 0.D0
        vf = 0.25D0 * vf
      ELSE IF (job_type == '3D') THEN
        resy = (b_n(ijk)*rgfn(ijk) - b_s(ijk)*rgfn(ijmk)) * indy(j)
        vf = vf + b_n(ijk) + b_s(ijk)
        vf = vf / 6.D0
      END IF

      IF (vf > 0.D0) THEN
        ivf = 1.D0 / vf
      ELSE
        ivf = 0.D0
      END IF
!
      res  = rgp(ijk) - rgpn(ijk) + dt * ivf * (resx+resy+resz)
!          - dt * (r1(ijk)+r2(ijk)+r3(ijk)+r4(ijk)+r5(ijk))

      RETURN
      END SUBROUTINE calc_res
!-----------------------------------------------------------------------
!
      SUBROUTINE calc_part_mass_flux(ijk)
!
!-----------------------------------------------------------------------
! ... compute the particle mass fluxes by using the stencil
! ... as defined in the 'subscr' module

      USE convective_mass_fluxes, ONLY: fmas, masf
      USE dimensions
      USE flux_limiters, ONLY: muscl
      USE gas_solid_velocity, ONLY: us, vs, ws
      USE gas_solid_density, ONLY: rlk
      USE set_indexes, ONLY: first_nb, first_rnb, third_nb, third_rnb
      USE set_indexes, ONLY: imjk, ijmk, ijkm
      USE set_indexes, ONLY: third_subscr
      USE control_flags, ONLY: job_type
!
      IMPLICIT NONE

      INTEGER :: ijk
      INTEGER :: is
!
        DO is = 1, nsolid

          IF (job_type == '2D') THEN
!
! ... Compute Fluxes by using First Order Upwind method
!
            ! ... assemble the first order computational stencils
            CALL first_nb(dens,rlk(:,is),ijk)
            CALL first_rnb(u,us(:,is),ijk)
            CALL first_rnb(w,ws(:,is),ijk)

            CALL masf(rsfe(ijk,is),  rsft(ijk,is),    &
                    rsfe(imjk,is), rsft(ijkm,is),   &
                    dens, u, w, ijk)
                     
            IF (muscl > 0) THEN
!
! ... Second order MUSCL correction
!
              ! ... add third neighbours to computational stencils
              CALL third_subscr( ijk ) 
              CALL third_nb(dens,rlk(:,is),ijk)
              CALL fmas(rsfe(ijk,is), rsft(ijk,is),    &
                      rsfe(imjk,is), rsft(ijkm,is),   &
                      dens, u, w, ijk)
              
            END IF

          ELSE IF (job_type == '3D') THEN
!
! ... Compute Fluxes by using First Order Upwind method
!
            ! ... assemble the first order computational stencils
            CALL first_nb(dens,rlk(:,is),ijk)
            CALL first_rnb(u,us(:,is),ijk)
            CALL first_rnb(v,vs(:,is),ijk)
            CALL first_rnb(w,ws(:,is),ijk)

            CALL masf(rsfe(ijk,is),  rsfn(ijk,is),  rsft(ijk,is),    &
                  rsfe(imjk,is), rsfn(ijmk,is), rsft(ijkm,is),   &
                  dens, u, v, w, ijk)

            IF (muscl > 0) THEN
!
! ... Second order MUSCL correction
!
              ! ... add third neighbours to computational stencils
              CALL third_subscr( ijk ) 
              CALL third_nb(dens,rlk(:,is),ijk)
              CALL fmas(rsfe(ijk,is),  rsfn(ijk,is),  rsft(ijk,is),    &
                    rsfe(imjk,is), rsfn(ijmk,is), rsft(ijkm,is),   &
                    dens, u, v, w, ijk)

            END IF

          END IF
        END DO

      END SUBROUTINE calc_part_mass_flux
!----------------------------------------------------------------------
!
      SUBROUTINE calc_eps( i, j, k, ijk )
!
!----------------------------------------------------------------------
! ... Use the Mass Balance equation of the solids
! ... to compute particle volumetric fractions and the void fraction.
!
      USE control_flags, ONLY: job_type, lpr
      USE dimensions
      USE gas_solid_density, ONLY: rlk, rlkn
      USE immersed_boundaries, ONLY: numx, numy, numz, immb
      USE pressure_epsilon, ONLY: p, ep
      USE set_indexes, ONLY: imjk, ijmk, ijkm
      USE time_parameters, ONLY: dt, time
      USE particles_constants, ONLY: inrl

      IMPLICIT NONE
      REAL*8 :: rlkx, rlky, rlkz, rls, rlk_tmp
      INTEGER, INTENT(IN) :: i, j, k, ijk
      INTEGER :: is
      REAL*8 :: vf, vf0, ivf
      INTEGER :: fx, fy, fz
      LOGICAL :: forced = .FALSE.

      vf0 = b_e(ijk) + b_w(ijk) + b_t(ijk) + b_b(ijk)

      rls = 0.D0
      DO is = 1, nsolid

        rlkx = (b_e(ijk)*rsfe(ijk,is) - b_w(ijk)*rsfe(imjk,is)) * indx(i) * inr(i)
        rlkz = (b_t(ijk)*rsft(ijk,is) - b_b(ijk)*rsft(ijkm,is)) * indz(k)

        IF (job_type == '2D') THEN
          rlky = 0.D0
          vf = 0.25D0 * vf0
        ELSE IF (job_type == '3D') THEN
          rlky = (b_n(ijk)*rsfn(ijk,is) - b_s(ijk)*rsfn(ijmk,is)) * indy(j)
          vf = vf0 + b_n(ijk) + b_s(ijk)
          vf = vf0 / 6.D0
        END IF

        IF (vf > 0.D0) THEN
          ivf = 1.D0 / vf
        ELSE
          ivf = 0.D0
        END IF

        rlk_tmp = rlkn( ijk, is ) - dt * ivf * ( rlkx + rlky + rlkz )
              !- dt * (r1(ijk)+r2(ijk)+r3(ijk)+r4(ijk)+r5(ijk))
        rlk_tmp = MAX( 0.0d0, rlk_tmp )

        rlk( ijk, is ) = rlk_tmp
        rls = rls + rlk( ijk, is ) * inrl(is)

      END DO
      
      IF( rls > 1.D0 ) THEN
        IF (lpr >= 1) THEN
          WRITE(8,*) ' warning1: mass is not conserved'
          WRITE(8,*) ' time, i, j, k ', time, i, j, k
          WRITE(8,*) ' rls, volfrac ', rls, vf
        END IF

        IF (immb >= 1) THEN
          fx = numx(ijk)
          IF (job_type == '2D') THEN
            fy = 0
          ELSE IF (job_type == '3D') THEN
            fy = numy(ijk)
          END IF
          fz = numz(ijk)
        
          forced = (fx/=0 .OR. fy/=0 .OR. fz/=0)
        END IF

        IF (.NOT.forced) CALL error('iter','mass is not conserved',1)
      ENDIF

      ep( ijk ) = 1.D0 - rls

      END SUBROUTINE calc_eps
!----------------------------------------------------------------------
!
      SUBROUTINE padjust(p, kros, d3, p3, omega, abeta)
!
!----------------------------------------------------------------------
! ... Correct the pressure field to minimize the gas mass residual

        REAL*8 :: p, d3, p3, omega, abeta
        INTEGER :: kros
        REAL*8, SAVE :: dp, d1, d2, p1, p2, d12


        IF( kros /= 3 ) THEN

          IF( d3 > 0.D0 ) THEN
            d1 = d3
            p1 = p3
            IF( kros == -1 ) THEN
              kros = 0
            ELSE IF( kros ==  1 ) THEN
              kros = 2
            END IF
          ELSE
            d2 = d3
            p2 = p3
            IF( kros == -1 ) THEN
              kros = 1
            ELSE IF( kros ==  0 ) THEN
              kros = 2
            END IF
          END IF

          IF( kros == 2 ) THEN

            ! ... Use secant method once, when the sign of dg changes
            d12 = 1.0d0 / ( d1 - d2 )
            p = ( d1 * p2 - d2 * p1 ) * d12
            abeta = ( p1 - p2 ) * d12
            kros  = 3

          ELSE
 
            ! ... Newton's method (iterated until the sign of dg changes)
            ! ... (abeta = -dp/dg) ...
            dp = -d3 * abeta
            ! ... with Under/Over-Relaxation ...
            dp =  dp * omega
            ! ... and a constrain on the maximum  correction.
            IF( -dp * dsign(1.D0,d3) > 2.5D-1 * p3 ) THEN
              dp = - 2.5D-1 * dsign( 1.D0, d3 ) * p3
            END IF 
            p = p + dp

          END IF

        ELSE

          ! ... Use two-sided secant method
          p = newp(d1, d2, d3, p1, p2, p3)

          IF(d3 > 0.D0) THEN
            d1 = d3
            p1 = p3
          ELSE
            d2 = d3
            p2 = p3
          END IF

        END IF

        p3 = p

        RETURN

      END SUBROUTINE padjust 
!----------------------------------------------------------------------
!
      REAL*8 FUNCTION newp(d1, d2, d3, p1, p2, p3)
!
!----------------------------------------------------------------------
! ... Bi-secant method
!
      USE dimensions
      IMPLICIT NONE
!
      REAL*8, INTENT(IN) :: d1, d2, d3, p1, p2, p3
      REAL*8 :: pa, pb
    
      pa = 0.0d0
      pb = 0.0d0
!
      IF (d1 .NE. d3) THEN
        pa = (d1*p3 - d3*p1) / (d1 - d3)
      END IF

      IF( d1*d3 <= 0.D0) THEN
        IF (d2 .NE. d3) THEN
          pb = (d2*p3 - d3*p2) / (d2 - d3)
        ELSE
          pb = 0.5D0 * (p1 + p3)
        END IF 
        IF(pb < p3 .OR. pb > p1) pb = 0.5D0 * (p1 + p3)
      ELSE
        IF(d1 == d3) pa = 0.5D0 * (p2 + p3)
        IF(pa < p2 .OR. pa > p3) pa = 0.5D0 * (p2 + p3)
        pb = (d2*p3 - d3*p2) / (d2-d3)
      END IF
!
      newp = 0.5D0 * (pa + pb)
!
      RETURN
      END FUNCTION newp
!----------------------------------------------------------------------
!
      SUBROUTINE betas(cnv, abt, ijk)
!
!----------------------------------------------------------------------
! ... Compute the derivative of the gas mass residual with respect to
! ... pressure 
!
      USE dimensions
      USE domain_decomposition, ONLY: ncint, myijk, meshinds
      USE convective_mass_fluxes, ONLY: upc_e, upc_n, upc_t
      USE convective_mass_fluxes, ONLY: upc_w, upc_s, upc_b
      USE control_flags, ONLY: job_type
      USE eos_gas, ONLY: csound
      USE gas_solid_density, ONLY: rog, rgp
      USE grid, ONLY: dx, dy, dz, flag, rb
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE pressure_epsilon, ONLY: p, ep
      USE set_indexes, ONLY: ipjk, imjk, ijpk, ijmk, ijkp, ijkm
      USE set_indexes, ONLY: ijke, ijkw, ijkn, ijks, ijkt, ijkb
      USE time_parameters, ONLY: dt, time

      IMPLICIT NONE
!
      REAL*8 :: iep_e, iep_w, iep_n, iep_s, iep_t, iep_b
      REAL*8 :: iepx, iepy, iepz
      REAL*8 :: dxm, dxp, dym, dyp, dzm, dzp
      REAL*8 :: indxm, indxp, indym, indyp, indzm, indzp
      REAL*8 :: rbeta, rags
      REAL*8 :: sqc
       
      INTEGER :: nfle, nflw, nfln, nfls, nflt, nflb
      INTEGER, INTENT(IN) :: ijk
      REAL*8, INTENT(OUT) :: cnv, abt
!
      REAL*8, PARAMETER :: delg=1.D-8
      INTEGER :: i,j,k,imesh
!
          CALL meshinds(ijk,imesh,i,j,k)

          dxp=dx(i)+dx(i+1)
          dxm=dx(i)+dx(i-1)
          dzp=dz(k)+dz(k+1)
          dzm=dz(k)+dz(k-1)
!
          indxp=1.D0/dxp
          indxm=1.D0/dxm
          indzp=1.D0/dzp
          indzm=1.D0/dzm
!
          nfle=flag(ipjk)
          nflw=flag(imjk)
          nflt=flag(ijkp)
          nflb=flag(ijkm)

          SELECT CASE (nfle)
            CASE(1,4,6)
              iep_e = ( dx(i+1)*ep(ijk)+dx(i)*ep(ijke) )*indxp*indxp*2.D0
              !iep_e = iep_e * upc_e
            CASE DEFAULT
              iep_e = 0.D0
          END SELECT

          SELECT CASE (nflw)
            CASE(1,4,6)
              iep_w = ( dx(i-1)*ep(ijk)+dx(i)*ep(ijkw) )*indxm*indxm*2.D0 
              !iep_w = iep_w * upc_w
            CASE DEFAULT
              iep_w = 0.D0
          END SELECT
!
          SELECT CASE (nflt)
            CASE(1,4,6)
              iep_t = ( dz(k+1)*ep(ijk)+dz(k)*ep(ijkt) )*indzp*indzp*2.D0 
              !iep_t = iep_t * upc_t 
            CASE DEFAULT
              iep_t = 0.0D0
          END SELECT

          SELECT CASE (nflb)
            CASE (1,4,6)
              iep_b = ( dz(k-1)*ep(ijk)+dz(k)*ep(ijkb) )*indzm*indzm*2.D0
              !iep_b = iep_b * upc_b
            CASE DEFAULT
              iep_b = 0.D0
          END SELECT
!
          IF (job_type == '3D') THEN

            dyp=dy(j)+dy(j+1)
            dym=dy(j)+dy(j-1)
            indyp=1.D0/dyp
            indym=1.D0/dym
            nfln=flag(ijpk)
            nfls=flag(ijmk)

            SELECT CASE (nfln)
              CASE (1,4,6)
                iep_n = ( dy(j+1)*ep(ijk)+dy(j)*ep(ijkn) )*indyp*indyp*2.D0 
                !iep_n = iep_n * upc_n 
              CASE DEFAULT
                iep_n = 0.0D0
            END SELECT

            SELECT CASE (nfls)
              CASE (1,4,6)
                iep_s = ( dy(j-1)*ep(ijk)+dy(j)*ep(ijks) )*indym*indym*2.D0
                !iep_s = iep_s * upc_s
              CASE DEFAULT
                iep_s = 0.D0
            END SELECT

          END IF

          iepx = (  rb(i) * iep_e + rb(i-1) * iep_w ) * indx(i) * inr(i) 
          iepz = (  iep_t + iep_b ) * indz(k)
	  IF (job_type == '2D') THEN
            iepy = 0.D0
	  ELSE IF (job_type == '3D') THEN
            iepy = (  iep_n + iep_s ) * indy(j)
          END IF
!
! ... Inverse of the squared sound velocity (perfect gas)
!
          CALL csound(sqc,rog(ijk),p(ijk))
          
          rags = 1.D0 / sqc
!
! ... rbeta = dD_g/dP
!
          rbeta = ep(ijk) * rags + dt**2 * (iepx+iepy+iepz)
!
          abt = 1.D0 / rbeta
          cnv = delg * rgp(ijk)
!
      RETURN
      END SUBROUTINE betas
!----------------------------------------------------------------------
      SUBROUTINE fill_cells
!
! ... Coefficients b(1:6) are multiplied to the numerical mass
! ... fluxes to take into account that some cell face could be
! ... immersed. b=1 if the cell face is external, b=0 if the
! ... cell face is inside the topography.
!
      USE control_flags, ONLY: job_type
      USE domain_decomposition, ONLY: ncint, meshinds
      USE grid, ONLY: z, zb
      USE immersed_boundaries, ONLY: immb
      USE immersed_boundaries, ONLY: topo_c, topo_x
      USE immersed_boundaries, ONLY: topo2d_c, topo2d_x, topo2d_y
      IMPLICIT NONE

      INTEGER :: i,j,k,ijk,imesh
      LOGICAL :: fillc
!
! ... Allocate and initialize coefficients
!
      ALLOCATE(b_e(ncint)); b_e = 1
      ALLOCATE(b_w(ncint)); b_w = 1
      ALLOCATE(b_t(ncint)); b_t = 1
      ALLOCATE(b_b(ncint)); b_b = 1
      IF (job_type == '3D') THEN
        ALLOCATE(b_n(ncint)); b_n = 1
        ALLOCATE(b_s(ncint)); b_s = 1
      END IF
!
      IF (immb >= 1) THEN
        DO ijk=1, ncint
          CALL meshinds(ijk,imesh,i,j,k)

          IF (job_type == '2D') THEN
            IF (z(k) < topo_x(i))     b_e(ijk) = 0 ! East
            IF (z(k) < topo_x(i-1))   b_w(ijk) = 0 ! West
            IF (zb(k) < topo_c(i))    b_t(ijk) = 0 ! Top
            IF (zb(k-1) < topo_c(i))  b_b(ijk) = 0 ! Bottom
          ELSE IF (job_type == '3D') THEN
            IF (z(k) < topo2d_x(i,j))     b_e(ijk) = 0 ! East
            IF (z(k) < topo2d_x(i-1,j))   b_w(ijk) = 0 ! West
            IF (z(k) < topo2d_y(i,j))     b_n(ijk) = 0 ! North
            IF (z(k) < topo2d_y(i,j-1))   b_s(ijk) = 0 ! South
            IF (zb(k) < topo2d_c(i,j))    b_t(ijk) = 0 ! Top
            IF (zb(k-1) < topo2d_c(i,j))  b_b(ijk) = 0 ! Bottom
          END IF

        END DO
      END IF
      
      END SUBROUTINE fill_cells
!----------------------------------------------------------------------
!     H E R E   S T A R T   T H E  O P T I M I Z E D   R O U T I N E S 
!----------------------------------------------------------------------
!
      SUBROUTINE opt_inner_loop(i, j, k, ijk, nit, d3, p3,abeta_,conv_, &
                                dgorig, nloop, cvg )
!
!----------------------------------------------------------------------

        USE dimensions
        USE pressure_epsilon, ONLY: p, ep
        USE gas_solid_density, ONLY: rog, rgp, rgpn, rlk, rlkn
        USE gas_solid_temperature, ONLY: tg
        USE eos_gas, ONLY: thermal_eosg, xgc, cg
        USE phases_matrix, ONLY: assemble_all_matrix, solve_all_velocities
        USE set_indexes, ONLY: third_nb, third_rnb, first_rnb, first_nb
        USE set_indexes, ONLY: subscr, imjk, ijmk, ijkm
        USE set_indexes, ONLY: ijke, ijkn, ijkt, ijkw, ijks, ijkb
        USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
        USE control_flags, ONLY: job_type
        USE time_parameters, ONLY: dt
        USE flux_limiters, ONLY: muscl
        USE convective_mass_fluxes, ONLY: fmas, masf
        USE gas_constants, ONLY: gmw, rgas, gas_type
        USE particles_constants, ONLY: inrl

        IMPLICIT NONE

        INTEGER :: nit, ijk, i, j, k, nloop
        REAL*8 :: d3, p3, abeta_, conv_, dgorig
        LOGICAL, INTENT(OUT) :: cvg
        REAL*8 :: resx, resy, resz, mg
        REAL*8 :: rlkx, rlky, rlkz, rls
        REAL*8 :: fe_, fn_, ft_, fw_, fs_, fb_
        REAL*8 :: fse_(max_nsolid), fsn_(max_nsolid), fst_(max_nsolid), &
                  fsw_(max_nsolid), fss_(max_nsolid), fsb_(max_nsolid)

        TYPE(stencil) :: ug_, vg_, wg_
        TYPE(stencil) :: us_(max_nsolid), vs_(max_nsolid), ws_(max_nsolid)
        TYPE(stencil) :: rgp_, rlk_(max_nsolid) 

        INTEGER :: loop, kros, ig, is
!
        ! ... These values of the mass fluxes will be updated
        ! ... in the inner_loop.
        !
        fse_(1:nsolid) = rsfe( ijk, 1:nsolid )
        fsw_(1:nsolid) = rsfe( imjk, 1:nsolid )
        fst_(1:nsolid) = rsft( ijk, 1:nsolid )
        fsb_(1:nsolid) = rsft( ijkm, 1:nsolid )

        fe_ = rgfe( ijk )
        fw_ = rgfe( imjk )
        ft_ = rgft( ijk )
        fb_ = rgft( ijkm )
        
        IF(job_type == '3D') THEN
          fsn_(1:nsolid) = rsfn( ijk, 1:nsolid )
          fss_(1:nsolid) = rsfn( ijmk, 1:nsolid )
          fn_ = rgfn( ijk )
          fs_ = rgfn( ijmk )
        END IF

        ! ... These are the stencil used in the inner_loop
        ! ... Only some values in the stencil are updated
        ! ... in the inner_loop
        !
        CALL first_nb(rgp_,rgp,ijk)
        CALL first_rnb(ug_,ug,ijk)
        CALL first_rnb(wg_,wg,ijk)
        IF (job_type == '3D') CALL first_rnb(vg_,vg,ijk)
        IF (muscl /= 0) THEN
          CALL third_nb(rgp_,rgp,ijk)
          CALL third_rnb(ug_,ug,ijk)
          CALL third_rnb(wg_,wg,ijk)
          IF (job_type == '3D') CALL third_rnb(vg_,vg,ijk)
        END IF

        DO is = 1, nsolid
          CALL first_nb(rlk_(is),rlk(:,is),ijk)
          CALL first_rnb(us_(is),us(:,is),ijk)
          CALL first_rnb(ws_(is),ws(:,is),ijk)
          IF (job_type == '3D') CALL first_rnb(vs_(is),vs(:,is),ijk)
          IF (muscl /= 0) THEN
            CALL third_nb(rlk_(is),rlk(:,is),ijk)
            CALL third_rnb(us_(is),us(:,is),ijk)
            CALL third_rnb(ws_(is),ws(:,is),ijk)
            IF (job_type == '3D') CALL third_rnb(vs_(is),vs(:,is),ijk)
          END IF
        END DO
!
        kros = -1
        DO loop = 1, inmax

          ! ... Correct the pressure at current cell
          ! ... by using successively the Newton's method
          ! ... (or the secant method where Newton is not
          ! ... applicable) and the two-sided secant method.
          ! ... The first time in a cell, skip pressure correction
          ! ... to compute once the particle velocities 
          ! ... and volumetric fractions ...
!
          IF( (loop > 1) .OR. (nit > 1) ) THEN
            CALL padjust( p(ijk), kros, d3, p3, omega, abeta_ )
          END IF
!
          ! ... Use equation of state to calculate gas density 
          ! ... from (new) pressure and (old) temperature.
!
          CALL thermal_eosg(rog(ijk),tg(ijk),p(ijk),xgc(:,ijk))
!
          rgp(ijk) = ep(ijk) * rog(ijk)
!
          ! ... Update gas and particles velocities using the 
          ! ... corrected pressure at current location. 
          ! ... Pressure at neighbour cells could still be wrong.

          CALL assemble_all_matrix(ijk)
          CALL solve_all_velocities(ijk)
!
          ! ... Update the stencils
          !
          ug_%c = ug( ijk )
          ug_%w = ug( imjk )
          wg_%c = wg( ijk )
          wg_%b = wg( ijkm )
          IF (job_type == '3D') THEN
            vg_%c = vg( ijk )
            vg_%s = vg( ijmk )
          END IF
          DO is = 1, nsolid
            us_(is)%c = us( ijk, is )
            us_(is)%w = us( imjk, is )
            ws_(is)%c = ws( ijk, is )
            ws_(is)%b = ws( ijkm, is )
            IF (job_type == '3D') THEN
              vs_(is)%c = vs( ijk, is )
              vs_(is)%s = vs( ijmk, is )
            END IF
          END DO
!
          ! ... By using the new velocity fields
          ! ... update particle gas mass fluxes and
          ! ... solid volumetric fractions
          !
          rls = 0.D0
          DO is = 1, nsolid

            IF (job_type == '2D') THEN

              CALL masf(fse_(is), fst_(is), fsw_(is), fsb_(is), &
                        rlk_(is), us_(is), ws_(is), ijk)
                     
              IF (muscl > 0) THEN
                CALL fmas(fse_(is), fst_(is), fsw_(is), fsb_(is), &
                          rlk_(is), us_(is), ws_(is), ijk)
              END IF

            ELSE IF (job_type == '3D') THEN

              CALL masf(fse_(is), fsn_(is), fst_(is), &
                        fsw_(is), fss_(is), fsb_(is), &
                        rlk_(is), us_(is), vs_(is), ws_(is), ijk)

              IF (muscl > 0) THEN
                CALL fmas(fse_(is), fsn_(is), fst_(is), &
                          fsw_(is), fss_(is), fsb_(is), &
                          rlk_(is), us_(is), vs_(is), ws_(is), ijk)
              END IF

            END IF

            rlkx = ( fse_(is) - fsw_(is) ) * indx(i) * inr(i)
            rlkz = ( fst_(is) - fsb_(is) ) * indz(k)
            IF (job_type == '3D') THEN
              rlky = ( fsn_(is) - fss_(is) ) * indy(j)
            ELSE
              rlky = 0.D0
            END IF

            rlk( ijk, is) = rlkn( ijk, is ) - dt * ( rlkx + rlky + rlkz )
            rlk( ijk, is) = MAX( 0.0d0, rlk(ijk,is) )

            rls = rls + rlk( ijk, is ) * inrl(is)

          END DO

          IF( rls > 1.D0 ) THEN
            WRITE(8,*) ' warning1: mass is not conserved'
            WRITE(8,*) ' i, j, k, rls ', i, j, k, rls
            CALL error('iter', 'warning: mass is not conserved',1)
          ENDIF

          ep( ijk ) = 1.D0 - rls
          rgp( ijk ) = ep( ijk ) * rog( ijk )

          ! ... Update the stencil
          !
          rlk_(:)%c = rlk(ijk,:)
          rgp_%c = rgp(ijk)


          ! ... update residual of the Mass Balance equation of the gas phase
          ! ... using guessed velocities and pressure.
          !
          IF (job_type == '2D') THEN

            CALL masf( fe_, ft_, fw_, fb_, rgp_, ug_, wg_, ijk )

            IF (muscl /= 0) THEN
              CALL fmas( fe_, ft_, fw_, fb_, rgp_, ug_, wg_, ijk )
            END IF

          ELSE IF (job_type == '3D') THEN

            CALL masf( fe_, fn_, ft_, fw_, fs_, fb_, rgp_, ug_, vg_, wg_, ijk)

            IF (muscl /= 0) THEN
              CALL fmas( fe_, fn_, ft_, fw_, fs_, fb_, rgp_, ug_, vg_, wg_, ijk)
            END IF

          END IF

          resx = ( fe_ - fw_ ) * indx(i) * inr(i)
          resz = ( ft_ - fb_ ) * indz(k)
          IF (job_type == '3D') THEN
            resy = ( fn_ - fs_ ) * indy(j)
          ELSE
            resy = 0.D0
          END IF

          dg = rgp(ijk) - rgpn(ijk) + dt * ( resx + resy + resz )
!
          IF ( ( DABS(dg) > conv_ ) .OR. ( DABS(dg) >= DABS(dgorig) ) ) THEN

            cvg = .FALSE.
            IF( nit == 1 .AND. loop == 1 ) dgorig = dg
            d3 = dg
            ! ... steepen the Newton's slope (accelerate)
            IF( kros < 2 .AND. loop == inmax ) abeta_ = 0.5D0 * inmax * abeta_

          ELSE IF ( DABS( dg ) <= conv_ ) THEN

            cvg = .TRUE.
            EXIT

          END IF

        END DO 

        ! ... Update the mass fluxes
        ! ... (used in the next outer iteration)
        !
        rsfe( ijk, 1:nsolid ) = fse_(1:nsolid) 
        rsfe( imjk, 1:nsolid ) = fsw_(1:nsolid) 
        rsft( ijk, 1:nsolid ) = fst_(1:nsolid)
        rsft( ijkm, 1:nsolid ) = fsb_(1:nsolid) 

        rgfe( ijk ) = fe_
        rgfe( imjk ) = fw_ 
        rgft( ijk ) = ft_
        rgft( ijkm ) = fb_
        
        IF(job_type == '3D') THEN
          rsfn( ijk, 1:nsolid ) = fsn_(1:nsolid)
          rsfn( ijmk, 1:nsolid ) = fss_(1:nsolid)
          rgfn( ijk ) = fn_
          rgfn( ijmk ) = fs_
        END IF

        nloop = nloop + loop

        RETURN
      
      END SUBROUTINE opt_inner_loop
!----------------------------------------------------------------------
!
      SUBROUTINE opt3_inner_loop(i,j,k,ijk,nit, d3, p3, abeta_,conv_, &
                                 dgorig, nloop, cvg )
!
!----------------------------------------------------------------------

        USE dimensions
        USE pressure_epsilon, ONLY: p, ep
        USE gas_solid_density, ONLY: rog, rgp, rgpn, rlk, rlkn
        USE gas_solid_temperature, ONLY: tg
        USE eos_gas, ONLY: thermal_eosg, xgc, cg
        USE phases_matrix, ONLY: matsvels_3phase
        USE set_indexes, ONLY: third_nb, third_rnb, first_rnb, first_nb
        USE set_indexes, ONLY: subscr, imjk, ijmk, ijkm
        USE set_indexes, ONLY: ijke, ijkn, ijkt, ijkw, ijks, ijkb
        USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
        USE control_flags, ONLY: job_type
        USE time_parameters, ONLY: dt
        USE flux_limiters, ONLY: muscl
        USE convective_mass_fluxes, ONLY: fmas, masf
        USE gas_constants, ONLY: gmw, rgas, gas_type
        USE particles_constants, ONLY: inrl

        IMPLICIT NONE

        INTEGER :: nit, ijk, i, j, k, nloop
        REAL*8 :: d3, p3, abeta_, conv_, dgorig
        LOGICAL, INTENT(OUT) :: cvg
        REAL*8 :: resx, resy, resz, mg
        REAL*8 :: rlkx, rlky, rlkz, rls
        REAL*8 :: rsfe1_ , rsfe2_
        REAL*8 :: rsfem1_ , rsfem2_
        REAL*8 :: rsft1_ , rsft2_
        REAL*8 :: rsftm1_ , rsftm2_
        REAL*8 :: rsfn1_ , rsfn2_
        REAL*8 :: rsfnm1_ , rsfnm2_

        REAL*8 :: rgfe_ , rgfem_
        REAL*8 :: rgft_ , rgftm_
        REAL*8 :: rgfn_ , rgfnm_

        TYPE(stencil) :: ug_, vg_, wg_, densg_
        TYPE(stencil) :: u1, v1, w1, dens1
        TYPE(stencil) :: u2, v2, w2, dens2

        INTEGER :: loop, kros, ig

        kros = -1

        rsfe1_ = rsfe( ijk, 1 )
        rsfe2_ = rsfe( ijk, 2 )
        rsfem1_ = rsfe( imjk, 1 )
        rsfem2_ = rsfe( imjk, 2 )
        rsft1_ = rsft( ijk, 1 )
        rsft2_ = rsft( ijk, 2 )
        rsftm1_ = rsft( ijkm, 1 )
        rsftm2_ = rsft( ijkm, 2 )
        rsfn1_ = rsfn( ijk, 1 )
        rsfn2_ = rsfn( ijk, 2 )
        rsfnm1_ = rsfn( ijmk, 1 )
        rsfnm2_ = rsfn( ijmk, 2 )

        rgfe_ = rgfe( ijk )
        rgfem_ = rgfe( imjk )
        rgft_ = rgft( ijk )
        rgftm_ = rgft( ijkm )
        rgfn_ = rgfn( ijk )
        rgfnm_ = rgfn( ijmk )

        CALL first_nb(dens1,rlk(:,1),ijk)
        CALL first_nb(dens2,rlk(:,2),ijk)
        CALL first_rnb(u1,us(:,1),ijk)
        CALL first_rnb(u2,us(:,2),ijk)
        CALL first_rnb(w1,ws(:,1),ijk)
        CALL first_rnb(w2,ws(:,2),ijk)
        CALL first_rnb(v1,vs(:,1),ijk)
        CALL first_rnb(v2,vs(:,2),ijk)
        CALL first_nb( densg_, rgp, ijk )
        CALL first_rnb( ug_, ug, ijk )
        CALL first_rnb( wg_, wg, ijk )
        CALL first_rnb( vg_, vg, ijk )

        IF (muscl /= 0) THEN
          CALL third_nb(dens1,rlk(:,1),ijk)
          CALL third_nb(dens2,rlk(:,2),ijk)
          CALL third_rnb(u1,us(:,1),ijk)
          CALL third_rnb(u2,us(:,2),ijk)
          CALL third_rnb(w1,ws(:,1),ijk)
          CALL third_rnb(w2,ws(:,2),ijk)
          CALL third_rnb(v1,vs(:,1),ijk)
          CALL third_rnb(v2,vs(:,2),ijk)
          CALL third_nb( densg_, rgp, ijk )
          CALL third_rnb( ug_, ug, ijk )
          CALL third_rnb( wg_, wg, ijk )
          CALL third_rnb( vg_, vg, ijk )
        END IF

        DO loop = 1, inmax
!
          ! ... Correct the pressure at current cell
          ! ... by using successively the Newton's method
          ! ... (or the secant method where Newton is not
          ! ... applicable) and the two-sided secant method.
          ! ... The first time in a cell, skip pressure correction
          ! ... to compute once the particle velocities 
          ! ... and volumetric fractions ...
!
          IF( (loop > 1) .OR. (nit > 1) ) THEN
            CALL padjust( p(ijk), kros, d3, p3, omega, abeta_ )
          END IF
!
          ! ... Use equation of state to calculate gas density 
          ! ... from (new) pressure and (old) temperature.
!
          mg = 0.D0
          DO ig = 1, ngas
            mg = mg + xgc(ig,ijk) * gmw( gas_type(ig) )
          END DO
          rog(ijk) = p(ijk) / ( rgas * tg(ijk) ) * mg
!
          rgp(ijk) = ep(ijk) * rog(ijk)
!
          ! ... Update gas and particles velocities using the 
          ! ... corrected pressure at current location. 
          ! ... Pressure at neighbour cells could still be wrong.

          CALL matsvels_3phase( amats(:,1,ijk), amats(:,2,ijk), amats(:,3,ijk), &
            amats(:,4,ijk), amats(:,5,ijk), amats(:,6,ijk), i, j, k, ijk )
!
          ug_%c = ug( ijk )
          ug_%w = ug( imjk )
          wg_%c = wg( ijk )
          wg_%b = wg( ijkm )
          vg_%c = vg( ijk )
          vg_%s = vg( ijmk )
          u1%c = us( ijk, 1 )
          u1%w = us( imjk, 1 )
          u2%c = us( ijk, 2 )
          u2%w = us( imjk, 2 )
          w1%c = ws( ijk, 1 )
          w1%b = ws( ijkm, 1 )
          w2%c = ws( ijk, 2 )
          w2%b = ws( ijkm, 2 )
          v1%c = vs( ijk, 1 )
          v1%s = vs( ijmk, 1 )
          v2%c = vs( ijk, 2 )
          v2%s = vs( ijmk, 2 )

          ! ... update particle and gas densities and the 
          ! ... gas mass residual 'dg'

          CALL masf( rsfe1_,  rsfn1_,  rsft1_, rsfem1_, rsfnm1_, rsftm1_,   &
                 dens1, u1, v1, w1, ijk)

          IF (muscl /= 0) THEN
            CALL fmas( rsfe1_,  rsfn1_,  rsft1_, rsfem1_, rsfnm1_, rsftm1_, &
                 dens1, u1, v1, w1, ijk)
          END IF

          rlkx = ( rsfe1_ - rsfem1_ ) * indx(i) * inr(i)
          rlky = ( rsfn1_ - rsfnm1_ ) * indy(j)
          rlkz = ( rsft1_ - rsftm1_ ) * indz(k)
          rlk( ijk, 1) = rlkn( ijk, 1 ) - dt * ( rlkx + rlky + rlkz )
          rlk( ijk, 1) = MAX( 0.0d0, rlk(ijk,1) )
          rls = rlk( ijk, 1 ) * inrl(1)

          CALL masf( rsfe2_,  rsfn2_,  rsft2_, rsfem2_, rsfnm2_, rsftm2_,   &
                 dens2, u2, v2, w2, ijk)

          IF (muscl /= 0) THEN
            CALL fmas( rsfe2_,  rsfn2_,  rsft2_, rsfem2_, rsfnm2_, rsftm2_, &
                 dens2, u2, v2, w2, ijk)
          END IF

          rlkx = ( rsfe2_ - rsfem2_ ) * indx(i) * inr(i)
          rlky = ( rsfn2_ - rsfnm2_ ) * indy(j)
          rlkz = ( rsft2_ - rsftm2_ ) * indz(k)
          rlk( ijk, 2) = rlkn( ijk, 2 ) - dt * ( rlkx + rlky + rlkz )
          rlk( ijk, 2) = MAX( 0.0d0, rlk(ijk,2) )
          rls = rls + rlk( ijk, 2 ) * inrl(2)

          IF( rls > 1.D0 ) THEN
            WRITE(8,*) ' warning1: mass is not conserved'
            WRITE(8,*) ' i, j, k, rls ', i, j, k, rls
            CALL error('iter', 'warning: mass is not conserved',1)
          ENDIF

          ep( ijk ) = 1.D0 - rls
          rgp( ijk ) = ep( ijk ) * rog( ijk )
          dens1%c = rlk( ijk, 1 )
          dens2%c = rlk( ijk, 2 )
          densg_%c = rgp(ijk)

          ! ... update residual of the Mass Balance equation of the gas phase
          ! ... using guessed velocities and pressure.

          ! ... update the particle volumetric fractions and the void fraction

          CALL masf( rgfe_ ,  rgfn_ ,  rgft_ , rgfem_ , rgfnm_ , rgftm_ ,   &
                      densg_, ug_, vg_, wg_, ijk )

          IF (muscl /= 0) THEN
            CALL fmas( rgfe_ ,  rgfn_ ,  rgft_ , rgfem_ , rgfnm_ , rgftm_ ,   &
                      densg_, ug_, vg_, wg_, ijk )
          END IF

          resx = ( rgfe_ - rgfem_ ) * indx(i) * inr(i)
          resz = ( rgft_ - rgftm_ ) * indz(k)
          resy = ( rgfn_ - rgfnm_ ) * indy(j)
          dg = rgp(ijk) - rgpn(ijk) + dt * ( resx + resy + resz )
!
          IF ( ( DABS(dg) > conv_ ) .OR. ( DABS(dg) >= DABS(dgorig) ) ) THEN

            cvg = .FALSE.
            IF( nit == 1 .AND. loop == 1 ) dgorig = dg
            d3 = dg
            ! ... steepen the Newton's slope (accelerate)
            IF( kros < 2 .AND. loop == inmax ) abeta_ = 0.5D0 * inmax * abeta_

          ELSE IF ( DABS( dg ) <= conv_ ) THEN

            cvg = .TRUE.
            EXIT

          END IF

        END DO 

        rsfe( ijk, 1 ) = rsfe1_
        rsfe( ijk, 2 ) = rsfe2_
        rsfe( imjk, 1 ) = rsfem1_
        rsfe( imjk, 2 ) = rsfem2_
        rsft( ijk, 1 ) = rsft1_
        rsft( ijk, 2 ) = rsft2_
        rsft( ijkm, 1 ) = rsftm1_
        rsft( ijkm, 2 ) = rsftm2_
        rsfn( ijk, 1 ) = rsfn1_
        rsfn( ijk, 2 ) = rsfn2_
        rsfn( ijmk, 1 ) = rsfnm1_
        rsfn( ijmk, 2 ) = rsfnm2_

        rgfe( ijk ) = rgfe_
        rgfe( imjk ) = rgfem_
        rgft( ijk ) = rgft_
        rgft( ijkm ) = rgftm_
        rgfn( ijk ) = rgfn_
        rgfn( ijmk ) = rgfnm_

        nloop = nloop + loop

        RETURN
      
      END SUBROUTINE opt3_inner_loop
!----------------------------------------------------------------------
      END MODULE iterative_solver
!----------------------------------------------------------------------
