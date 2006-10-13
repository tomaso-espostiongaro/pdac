!----------------------------------------------------------------------
      MODULE iterative_solver
!----------------------------------------------------------------------
      USE grid, ONLY: dx, dy, dz
      USE grid, ONLY: indx, indy, indz, inr
      USE set_indexes, ONLY: stencil
      USE io_files, ONLY: testunit

      IMPLICIT NONE
!
! ... convective mass fluxes
!
      REAL*8, DIMENSION(:),   ALLOCATABLE :: rgfe, rgfn, rgft
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: rsfe, rsfn, rsft

      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: amats

      REAL*8  :: omega, dg, delg
      INTEGER :: inmax, maxout, nit
      INTEGER :: optimization
      INTEGER :: ierr

      INTEGER :: b_e, b_w, b_t, b_b, b_n, b_s
      REAL*8 :: ivf

      TYPE(stencil) :: u, v, w, dens         

      PRIVATE :: u, v, w, dens, b_e, b_w, b_t, b_b, b_n, b_s, ivf
      PRIVATE :: ierr

      SAVE
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
      USE domain_mapping, ONLY: ncint, ncdom, data_exchange
      USE domain_mapping, ONLY: myijk, meshinds
      USE enthalpy_matrix, ONLY: ftem
      USE environment, ONLY: cpclock, timing
      USE flux_limiters, ONLY: muscl
      USE gas_solid_temperature, ONLY: tg, ts, sieg, sies
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      USE gas_solid_density, ONLY: rog, rgp, rgpn, rlk, rlkn
      USE grid, ONLY: flag
      USE immersed_boundaries, ONLY: immb, faces
      USE indijk_module
      USE output_dump, ONLY: cell_report
      USE parallel, ONLY: mpime, root
      USE particles_constants, ONLY: rl, inrl
      USE phases_matrix, ONLY: assemble_matrix, assemble_all_matrix
      USE phases_matrix, ONLY: solve_velocities, solve_all_velocities
      USE phases_matrix, ONLY: matspre_3phase
      USE pressure_epsilon, ONLY: p, ep
      USE set_indexes, ONLY: imjk, ijmk, ijkm, ijkn, &
                             ijke, ijkw, ijkt, ijkb
      USE set_indexes, ONLY: first_subscr, third_subscr
      USE tilde_energy, ONLY: htilde
      USE tilde_momentum, ONLY: tilde, appu, appv, appw
      USE time_parameters, ONLY: time, dt, timestart, tpr, sweep
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
      INTEGER :: nloop, mustit
      INTEGER :: n1, n2
      REAL*8  :: avloop
      INTEGER :: i, j, k, ijk

      LOGICAL, ALLOCATABLE :: converge(:)
      LOGICAL :: cvg
      LOGICAL :: compute
!
! ... Initialize the cell fractions for immersed boundaries
!
      b_e = 1; b_w = 1; b_t = 1; b_b = 1; b_n = 1; b_s = 1; ivf = 1.D0
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
! ... First "predictor" step.
! ... Assemble and solve the explicit phase matrix 
! ... to update velocity fields. New velocities are
! ... biassed by the wrong pressure field. The best
! ... estimate for the pressure field is given by
! ... the pressure at previous time-step.
!
      DO ijk = 1, ncint
        compute = BTEST(flag(ijk),0)
        IF ( compute ) THEN
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
      IF( ((nsolid /= 2) .OR. (job_type == '2D')) .AND. (optimization >= 3) ) &
        optimization = 2

      IF( optimization == 3 ) THEN
        ALLOCATE( amats( 6, 6, ncint ) )
      END IF
!
! ... Guess an approximate value for gas density fluxes
! ... using the estimates of velocity and pressure.
!
      DO ijk = 1, ncint
        compute = BTEST(flag(ijk),0)
        IF ( compute  ) THEN
          CALL first_subscr( ijk )
          CALL calc_gas_mass_flux( ijk )
          !
          ! ... compute the derivative of the gas mass residual
          ! ... with respect to gas pressure
          ! ... ('abeta' is kept constant during iterative sweep)
          CALL betas( conv( ijk ), abeta( ijk ), ijk )

          IF(optimization==3) THEN
            CALL meshinds(ijk,imesh,i,j,k)
            CALL matspre_3phase(amats(:,1,ijk),amats(:,2,ijk),amats(:,3,ijk),&
                                amats(:,4,ijk),amats(:,5,ijk),amats(:,6,ijk),&
                                i, j, k, ijk )
          END IF
        END IF
      END DO

      !CALL test_fluxes
!
! ... Here Start the external iterative sweep.
!/////////////////////////////////////////////////////////////////////
!
! ... The correction equation are iterated on the mesh
! ... to propagate the updated velocities and pressure
! ... on each cell to its neighbours.
!
      omega0 = omega
!
      timconv = 0.0d0
!
      sor_loop: DO nit = 1, maxout
         mustit = 1
         ierr = 0
!
         ! ... Compute fluxes and forces in the momentum equation
         ! ... (for fully implicit solution)
         !
         IF (implicit_fluxes) THEN

           !
           ! ... Need to exchange velocities to compute
           ! ... either first or second-order momentum fluxes
           !
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
           !
           converge(ijk) = .FALSE.

           ! ... Fluid cells or internal immersed boundary cell
           compute  = BTEST(flag(ijk),0)
           IF( compute ) THEN

             CALL meshinds(ijk,imesh,i,j,k)
             CALL first_subscr(ijk)

             ! ... Compute the volumes partially filled by the
             ! ... topography
             !
             IF (immb == 1) CALL faces(ijk, b_e, b_w, b_t, b_b, b_n, b_s, ivf)
             IF (muscl > 0) CALL third_subscr(ijk)

             ! ... Compute locally the residual 'dg' of the
             ! ... mass balance equation of the gas phase and store it
             !
             CALL calc_res(i,j,k,ijk,dg)
             dgorig = dg

             IF( ABS( dg ) <= conv( ijk ) ) THEN

               ! ... If the residual is lower then the prescribed limit
               ! ... compute the gas and particle densities and proceed 
               ! ... to next cell
               !
               converge( ijk ) = .TRUE.

               ! ... The counter of converging cells is increased
               !
               n1 = n1 + 1
              
               CALL calc_eps( i, j, k, ijk  )
               rgp( ijk ) = ep( ijk ) * rog( ijk )

             ELSE IF ( ABS( dg ) > conv( ijk ) ) THEN

               ! ... The counter of non-converging cells is increased
               !
               n2 = n2 + 1

               ! ... If the residual is higher then the prescribed limit
               ! ... start the inner (in-cell) iterative loop
               ! ... to correct pressure, gas and particle velocities
               ! ... and gas density (as a function of the fixed temperature).
               !
               d3 = dg
               p3 = p( ijk )

               IF (optimization == 1) THEN
                 !
                 CALL inner_loop(i, j, k, ijk, nit, d3, p3, &
                                 abeta(ijk),conv(ijk), &
                                 dgorig, nloop, cvg)
               ELSE IF (optimization == 2) THEN
                 !
                 CALL opt_inner_loop(i, j, k, ijk, nit, d3, p3, &
                                 abeta(ijk),conv(ijk), &
                                 dgorig, nloop, cvg)
               ELSE IF (optimization == 3) THEN
                 !
                 CALL opt3_inner_loop(i, j, k, ijk, nit, d3, p3, &
                                  abeta(ijk), conv(ijk), &
                                  dgorig, nloop, cvg)
               END IF

             END IF

           ELSE 
             
             converge(ijk) = .TRUE.
           
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
! --> n1 : the n umber of cells converged
! --> n2 : the number of cells not converged
! --> avloop : the averaged number of inner iterations per cells
! --> timconv: the time for the whole mesh sweep
!
         IF (lpr > 0) THEN
           IF( n2 > 0 ) THEN
             avloop = REAL(nloop) / REAL(n2)
           ELSE
             avloop = REAL(nloop)
           END IF
           WRITE(testunit, fmt="( I10, 2X, F5.2, 2X, I10, 2X, F10.3, I10)" ) &
                            n2, avloop, n1, timconv(nit), COUNT(converge)
           CALL myflush( testunit )
         END IF
!*******************************************************************

! ... Exchange all updated physical quantities on boundaries
! ... before starting a new sweep on the mesh.
! ... Please notice that velocity MUST NOT be exchanged,
! ... since the West, South and Bottom are computed in-cell
! ... and East, North and Top are not used to compute scalar fluxes
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
          ! ... solve pressure and interphase coupling
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
        IF( mustit == 0 ) THEN
!*******************************************************************
! ... write out the final number of iterations
!
          IF (lpr > 0) THEN
            WRITE(testunit,277) nit
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
       IF( MOD( nit, 100 ) == 0 ) THEN
         omega = 0.5D0 * (omega + 1.D0)
         IF (lpr > 0) THEN
           WRITE(testunit, fmt="('  reducing relaxation parameter omega')")
           WRITE(testunit, fmt="('  new value = ',F12.4)") omega
         END IF
       END IF
!
! ... Check the closure relation for solid phases on all processors
!
        CALL parallel_sum_integer(ierr, 1)
        !IF (ierr > 1) CALL data_exchange(flag)
        !IF (ierr > 1) CALL error('iter','solid fraction exceeded 1',ierr)
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
      IF( lpr > 0 ) WRITE(testunit,280)'Time for iterative solver: ',timiter
 280  FORMAT(2X,A27,F8.3) 

!
!********************************************************************
! ... If the iterative sweep concluded without convergence
! ... report the number of cells where the procedure does not converge
!
      IF( mustit /= 0 ) THEN
!
        IF (lpr > 0) THEN
          WRITE(testunit,700) nit, (time+dt)
          WRITE(testunit,*) 'convergence on proc ',mpime,' : ', ALL(converge)
          IF(.NOT.ALL(converge)) &
            WRITE(testunit,*) 'cells not converged (imesh,i,j,k): '
          DO ijk = 1, ncint
            IF ( .NOT. converge( ijk ) ) THEN
              CALL meshinds( ijk , imesh, i , j , k )
              WRITE(testunit,*) imesh, i , j , k
              CALL cell_report(testunit, ijk, imesh, i, j, k)
              CALL correct_particles(ijk, imesh, i, j, k)
            END IF
          END DO
 700      FORMAT('max number of iterations (',I5,') reached at time: ', F8.3)
        ELSE
          DO ijk = 1, ncint
            IF ( .NOT. converge( ijk ) ) THEN
              CALL correct_particles(ijk, imesh, i, j, k)
            END IF
          END DO
        END IF
!
        ! ... CRASH! ...
        !
        !CALL error( ' iter ', 'max number of iters exceeded ', 1)
        omega = omega0
        !
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
!
      RETURN
      END SUBROUTINE iter
!----------------------------------------------------------------------
!
      SUBROUTINE inner_loop( i, j, k, ijk, nit, d3, p3, abeta_, conv_,  &
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
      USE set_indexes, ONLY: imjk, ijmk, ijkm, ipjk, ijpk, ijkp
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      USE control_flags, ONLY: job_type

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nit, ijk, i, j, k
      REAL*8, INTENT(INOUT)  :: conv_
      REAL*8, INTENT(INOUT)  :: d3, p3, abeta_, dgorig
      LOGICAL, INTENT(OUT) :: cvg

      REAL*8  :: resx, resy, resz
      INTEGER :: loop, kros, nloop

      REAL*8 :: xgcl(max_ngas)
      xgcl(1:ngas) = xgc(ijk,:)

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
          CALL padjust(p(ijk), kros, d3, p3, omega, abeta_ )
        END IF
!
        ! ... Use equation of state to calculate gas density 
        ! ... from (new) pressure and (old) temperature.
        ! ... (notice that the explicit solution of the enthalpy
        ! ... equations affects only this updating)
!
        CALL thermal_eosg(rog(ijk),tg(ijk),p(ijk),xgcl(:))
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
        ! ... This procedure is optimized (optimization = 2) !
        !
        CALL calc_gas_mass_flux( ijk )
        CALL calc_res( i, j, k, ijk, dg)
!
        IF ( ( DABS(dg) > conv_ .OR. DABS(dg) >= DABS(dgorig) ) ) THEN

          cvg = .FALSE.
          IF( nit == 1 .AND. loop == 1 ) dgorig = dg
          d3 = dg
          ! ... steepen the Newton's slope (accelerate)
          IF( kros < 2 .AND. loop == inmax ) abeta_ = 0.5D0 * inmax * abeta_

        ELSE IF ( DABS(dg) <= conv_ ) THEN

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
      USE grid, ONLY: flag
      USE set_indexes, ONLY: first_nb, first_rnb, third_nb, third_rnb
      USE set_indexes, ONLY: imjk, ijmk, ijkm
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
          ! ... (includes only the first neighbours)
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
            CALL third_nb( dens, rgp, ijk)
            CALL fmas( rgfe( ijk ), rgft( ijk ), rgfe( imjk ), rgft( ijkm ), &
                       dens, u, w, ijk )
          END IF


        ELSE IF (job_type == '3D') THEN
!
! ... Compute Fluxes by using First Order Upwind method
!
          ! ... assemble the first order computational stencils
          ! ... (includes only the first neighbours)
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
      USE set_indexes, ONLY: imjk, ijmk, ijkm
      USE gas_solid_density, ONLY: rog, rgp, rgpn, rlk, rlkn
      USE control_flags, ONLY: job_type
      USE time_parameters, ONLY: dt

      IMPLICIT NONE
      REAL*8 :: resx, resy, resz
      REAL*8, INTENT(OUT) :: res
      INTEGER, INTENT(IN) :: i, j, k, ijk
!      
      resx = ( b_e * rgfe(ijk) - b_w * rgfe(imjk) ) * indx(i) * inr(i)
      resz = ( b_t * rgft(ijk) - b_b * rgft(ijkm) ) * indz(k)
!
      IF (job_type == '2D') THEN
        resy = 0.D0
      ELSE IF (job_type == '3D') THEN
        resy = (b_n * rgfn(ijk) - b_s * rgfn(ijmk)) * indy(j)
      END IF
!
      res  = rgp(ijk) - rgpn(ijk) + dt * ivf * (resx+resy+resz)
!      
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
      USE grid, ONLY: flag
      USE set_indexes, ONLY: first_nb, first_rnb, third_nb, third_rnb
      USE set_indexes, ONLY: imjk, ijmk, ijkm
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
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      USE grid, ONLY: flag, noslip_wall
      USE pressure_epsilon, ONLY: p, ep
      USE set_indexes, ONLY: imjk, ijmk, ijkm
      USE time_parameters, ONLY: dt, time
      USE particles_constants, ONLY: inrl

      IMPLICIT NONE
      REAL*8 :: rlkx, rlky, rlkz, rls, rlk_tmp
      INTEGER, INTENT(IN) :: i, j, k, ijk
      INTEGER :: is
      INTEGER :: fx, fy, fz
      LOGICAL :: forced = .FALSE.

      rls = 0.D0
      DO is = 1, nsolid

        rlkx = (b_e * rsfe(ijk,is) - b_w * rsfe(imjk,is)) * indx(i) * inr(i)
        rlkz = (b_t * rsft(ijk,is) - b_b * rsft(ijkm,is)) * indz(k)

        IF (job_type == '2D') THEN
          rlky = 0.D0
        ELSE IF (job_type == '3D') THEN
          rlky = (b_n * rsfn(ijk,is) - b_s * rsfn(ijmk,is)) * indy(j)
        END IF

        rlk_tmp = rlkn( ijk, is ) - dt * ivf * ( rlkx + rlky + rlkz )
              !
              !- dt * (r1(ijk)+r2(ijk)+r3(ijk)+r4(ijk)+r5(ijk))

        rlk( ijk, is ) = MAX( 0.0d0, rlk_tmp )
        rls = rls + rlk( ijk, is ) * inrl(is)

      END DO

      ep( ijk ) = 1.D0 - rls
!      
!**********************************************************************
! ... Error report
!
      IF( rls > 1.D0 ) THEN
        IF (lpr > 0) THEN
          WRITE(testunit,*) ' Mass is not conserved', flag(ijk)
          WRITE(testunit,*) ' time, i, j, k ', time, i, j, k
          WRITE(testunit,*) ' rls, volfrac ', rls, 1.D0/ivf
        END IF
        ! ... cells exceeding the maximum packing are 'frozen'
        !IF (flag(ijkm) == noslip_wall) flag(ijk) = noslip_wall
        ierr = ierr + 1
      ENDIF
!
!**********************************************************************
!
      RETURN
      END SUBROUTINE calc_eps
!----------------------------------------------------------------------
!
      SUBROUTINE padjust(p, kros, d3, p3, omega, abeta)
!
!----------------------------------------------------------------------
! ... Correct the pressure field to minimize the gas mass residual

        REAL*8, INTENT(IN) :: omega
        REAL*8 :: p, d3, p3, abeta
        INTEGER :: kros
        REAL*8, SAVE :: dp, d1, d2, p1, p2, d12

        ! ... First set the parameters for non-linear zero searching
        !
        IF( kros /= 3 ) THEN

          IF( d3 > 0.D0 ) THEN
            d1 = d3
            p1 = p3
            IF( kros == -1 ) THEN
              kros = 0
            ELSE IF( kros ==  1 ) THEN
              kros = 2
            END IF
          ELSE IF ( d3 <= 0.D0 ) THEN
            d2 = d3
            p2 = p3
            IF( kros == -1 ) THEN
              kros = 1
            ELSE IF( kros ==  0 ) THEN
              kros = 2
            END IF
          END IF

          IF ( kros == 2 ) THEN

            ! ... Use secant method once, when the sign of dg changes
            d12 = 1.0d0 / ( d1 - d2 )
            p = ( d1 * p2 - d2 * p1 ) * d12
            abeta = ( p1 - p2 ) * d12
            kros  = 3

          ELSE IF ( kros < 2 ) THEN
 
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

        ELSE IF (kros == 3) THEN
          !
          ! ... Use two-sided secant method
          !
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
      USE domain_mapping, ONLY: ncint, myijk, meshinds
      USE convective_mass_fluxes, ONLY: upc_e, upc_n, upc_t
      USE convective_mass_fluxes, ONLY: upc_w, upc_s, upc_b
      USE control_flags, ONLY: job_type
      USE eos_gas, ONLY: csound
      USE gas_solid_density, ONLY: rog, rgp
      USE grid, ONLY: dx, dy, dz, flag, rb
      USE grid, ONLY: fluid, free_io, nrfree_io
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
            CASE(fluid,free_io,nrfree_io)
              iep_e = ( dx(i+1)*ep(ijk)+dx(i)*ep(ijke) )*indxp*indxp*2.D0
              !iep_e = iep_e * upc_e
            CASE DEFAULT
              iep_e = 0.D0
          END SELECT

          SELECT CASE (nflw)
            CASE(fluid,free_io,nrfree_io)
              iep_w = ( dx(i-1)*ep(ijk)+dx(i)*ep(ijkw) )*indxm*indxm*2.D0 
              !iep_w = iep_w * upc_w
            CASE DEFAULT
              iep_w = 0.D0
          END SELECT
!
          SELECT CASE (nflt)
            CASE(fluid,free_io,nrfree_io)
              iep_t = ( dz(k+1)*ep(ijk)+dz(k)*ep(ijkt) )*indzp*indzp*2.D0 
              !iep_t = iep_t * upc_t 
            CASE DEFAULT
              iep_t = 0.0D0
          END SELECT

          SELECT CASE (nflb)
            CASE (fluid,free_io,nrfree_io)
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
              CASE (fluid,free_io,nrfree_io)
                iep_n = ( dy(j+1)*ep(ijk)+dy(j)*ep(ijkn) )*indyp*indyp*2.D0 
                !iep_n = iep_n * upc_n 
              CASE DEFAULT
                iep_n = 0.0D0
            END SELECT

            SELECT CASE (nfls)
              CASE (fluid,free_io,nrfree_io)
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
          rbeta = ep(ijk) * rags 
          rbeta = rbeta + dt**2 * (iepx+iepy+iepz)
!
          abt = 1.D0 / rbeta
          cnv = delg * rgp(ijk)
!
      RETURN
      END SUBROUTINE betas
!----------------------------------------------------------------------
!
!     H E R E   S T A R T   T H E   O P T I M I Z E D   R O U T I N E S 
!
!----------------------------------------------------------------------
!
      SUBROUTINE opt_inner_loop(i, j, k, ijk, nit, d3, p3,abeta_,conv_, &
                                dgorig, nloop, cvg )
!
!----------------------------------------------------------------------
!
! ... This is an optimized version of the inner_loop, where the stencils
! ... of the various fields are COMPUTED ONCE at the beginning of the
! ... loop and only updated at the end. In fact, not every element of
! ... the stencil is updated at each inner loop.
!
        USE control_flags, ONLY: job_type, lpr
        USE dimensions
        USE pressure_epsilon, ONLY: p, ep
        USE gas_solid_density, ONLY: rog, rgp, rgpn, rlk, rlkn
        USE gas_solid_temperature, ONLY: tg
        USE grid, ONLY: flag, noslip_wall
        USE eos_gas, ONLY: thermal_eosg, xgc
        USE phases_matrix, ONLY: assemble_all_matrix, solve_all_velocities
        USE set_indexes, ONLY: third_nb, third_rnb, first_rnb, first_nb
        USE set_indexes, ONLY: imjk, ijmk, ijkm
        USE set_indexes, ONLY: ijke, ijkn, ijkt, ijkw, ijks, ijkb
        USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
        USE control_flags, ONLY: job_type
        USE time_parameters, ONLY: dt, time
        USE flux_limiters, ONLY: muscl
        USE convective_mass_fluxes, ONLY: fmas, masf
        USE gas_constants, ONLY: gmw, rgas, gas_type
        USE particles_constants, ONLY: inrl

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: ijk, i, j, k, nit
        INTEGER :: nloop
        REAL*8 :: d3, p3, abeta_, conv_, dgorig
        LOGICAL, INTENT(OUT) :: cvg
        REAL*8 :: resx, resy, resz, mg
        REAL*8 :: rlkx, rlky, rlkz, rls, rlk_tmp
        REAL*8 :: rgfe_, rgfn_, rgft_, rgfw_, rgfs_, rgfb_
        REAL*8 :: rsfe_(max_nsolid), rsfn_(max_nsolid), rsft_(max_nsolid), &
                  rsfw_(max_nsolid), rsfs_(max_nsolid), rsfb_(max_nsolid)

        TYPE(stencil) :: ug_, vg_, wg_
        TYPE(stencil) :: us_(max_nsolid), vs_(max_nsolid), ws_(max_nsolid)
        TYPE(stencil) :: rgp_, rlk_(max_nsolid) 

        INTEGER :: loop, kros, ig, is
!
        REAL*8 :: xgcl(max_ngas)
!
! ... Here prepare the stencils and the fluxes
!
        xgcl(1:ngas) = xgc(ijk,:)

        ! ... These values of the mass fluxes will be updated
        ! ... in the inner_loop. Storing them into smaller arrays
        ! ... optimizes memory access.
        !
        rsfe_(1:nsolid) = rsfe( ijk, 1:nsolid )
        rsfw_(1:nsolid) = rsfe( imjk, 1:nsolid )
        rsft_(1:nsolid) = rsft( ijk, 1:nsolid )
        rsfb_(1:nsolid) = rsft( ijkm, 1:nsolid )
        rgfe_ = rgfe( ijk )
        rgfw_ = rgfe( imjk )
        rgft_ = rgft( ijk )
        rgfb_ = rgft( ijkm )
        !
        IF(job_type == '3D') THEN
          rsfn_(1:nsolid) = rsfn( ijk, 1:nsolid )
          rsfs_(1:nsolid) = rsfn( ijmk, 1:nsolid )
          rgfn_ = rgfn( ijk )
          rgfs_ = rgfn( ijmk )
        END IF

        ! ... These are the stencil used in the inner_loop
        ! ... Only few values in the stencil are updated
        ! ... in the inner_loop
        !
        CALL first_nb(rgp_,rgp,ijk)
        CALL first_rnb(ug_,ug,ijk)
        CALL first_rnb(wg_,wg,ijk)
        IF (job_type == '3D') CALL first_rnb(vg_,vg,ijk)
        !
        IF (muscl > 0) THEN
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

          IF (muscl > 0) THEN
            CALL third_nb(rlk_(is),rlk(:,is),ijk)
            CALL third_rnb(us_(is),us(:,is),ijk)
            CALL third_rnb(ws_(is),ws(:,is),ijk)
            IF (job_type == '3D') CALL third_rnb(vs_(is),vs(:,is),ijk)
          END IF
        END DO
!
        kros = -1
!
! ... Here start the inner loop
!
        inloop: DO loop = 1, inmax

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
          CALL thermal_eosg(rog(ijk),tg(ijk),p(ijk),xgcl(:))
          rgp(ijk) = ep(ijk) * rog(ijk)
!
          ! ... Update gas and particles velocities using the 
          ! ... corrected pressure at current location. 
          ! ... Pressure at neighbour cells could still be wrong.
          !
          CALL assemble_all_matrix(ijk)
          CALL solve_all_velocities(ijk)
!
          ! ... Update the corresponding stencil elements
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
          ! ... By using the new velocity fields,
          ! ... update particle gas mass fluxes and
          ! ... solid volumetric fractions
          !
          rls = 0.D0
          DO is = 1, nsolid

            IF (job_type == '2D') THEN

              CALL masf(rsfe_(is), rsft_(is), rsfw_(is), rsfb_(is), &
                        rlk_(is), us_(is), ws_(is), ijk)
                     
              IF (muscl > 0) THEN
                CALL fmas(rsfe_(is), rsft_(is), rsfw_(is), rsfb_(is), &
                          rlk_(is), us_(is), ws_(is), ijk)
              END IF

            ELSE IF (job_type == '3D') THEN

              CALL masf(rsfe_(is), rsfn_(is), rsft_(is), &
                        rsfw_(is), rsfs_(is), rsfb_(is), &
                        rlk_(is), us_(is), vs_(is), ws_(is), ijk)

              IF (muscl > 0) THEN
                CALL fmas(rsfe_(is), rsfn_(is), rsft_(is), &
                          rsfw_(is), rsfs_(is), rsfb_(is), &
                          rlk_(is), us_(is), vs_(is), ws_(is), ijk)
              END IF

            END IF

            rlkx = ( b_e * rsfe_(is) - b_w * rsfw_(is) ) * indx(i) * inr(i)
            rlkz = ( b_t * rsft_(is) - b_b * rsfb_(is) ) * indz(k)
            IF (job_type == '3D') THEN
              rlky = ( b_n * rsfn_(is) - b_s * rsfs_(is) ) * indy(j)
            ELSE
              rlky = 0.D0
            END IF

            rlk_tmp = rlkn( ijk, is ) - dt * ivf *  ( rlkx + rlky + rlkz )
            rlk( ijk, is) = MAX( 0.0d0, rlk_tmp )

            rls = rls + rlk( ijk, is ) * inrl(is)

          END DO

          IF( rls > 1.D0 ) THEN
            IF (lpr > 0) THEN
              WRITE(testunit,*) ' Mass is not conserved', flag(ijk)
              WRITE(testunit,*) ' i, j, k, rls ', i, j, k, rls
              WRITE(testunit,*) ' rls, volfrac ', rls, 1.D0/ivf
            ENDIF
            !IF (flag(ijkm) == noslip_wall) flag(ijk) = noslip_wall
            ierr = ierr + 1
          ENDIF

          ep( ijk ) = 1.D0 - rls
          rgp( ijk ) = ep( ijk ) * rog( ijk )

          ! ... Update the corresponding stencils
          !
          rlk_(1:nsolid)%c = rlk(ijk,1:nsolid)
          rgp_%c = rgp(ijk)

          ! ... update residual of the Mass Balance equation of the gas phase
          ! ... using guessed velocities and pressure.
          !
          IF (job_type == '2D') THEN

            CALL masf( rgfe_, rgft_, rgfw_, rgfb_, rgp_, ug_, wg_, ijk )

            IF (muscl > 0) THEN
              CALL fmas( rgfe_, rgft_, rgfw_, rgfb_, rgp_, ug_, wg_, ijk )
            END IF

          ELSE IF (job_type == '3D') THEN

            CALL masf( rgfe_, rgfn_, rgft_, rgfw_, rgfs_, rgfb_, rgp_, ug_, vg_, wg_, ijk)

            IF (muscl > 0) THEN
              CALL fmas( rgfe_, rgfn_, rgft_, rgfw_, rgfs_, rgfb_, rgp_, ug_, vg_, wg_, ijk)
            END IF

          END IF

          resx = ( b_e * rgfe_ - b_w * rgfw_ ) * indx(i) * inr(i)
          resz = ( b_t * rgft_ - b_b * rgfb_ ) * indz(k)
          IF (job_type == '3D') THEN
            resy = ( b_n * rgfn_ - b_s * rgfs_ ) * indy(j)
          ELSE
            resy = 0.D0
          END IF

          dg = rgp(ijk) - rgpn(ijk) + dt * ivf * ( resx + resy + resz )
!
          IF ( ( DABS(dg) > conv_ ) .OR. ( DABS(dg) >= DABS(dgorig) ) ) THEN

            cvg = .FALSE.
            IF( nit == 1 .AND. loop == 1 ) dgorig = dg
            d3 = dg
            ! ... steepen the Newton's slope (accelerate)
            IF( kros < 2 .AND. loop == inmax ) abeta_ = 0.5D0 * inmax * abeta_

          ELSE IF ( DABS( dg ) <= conv_ ) THEN

            cvg = .TRUE.
            EXIT inloop

          END IF

        END DO inloop
!
! ... Here rebuild the arrays
!
        ! ... Update the mass fluxes
        ! ... (used in the next outer iteration)
        !
        rsfe( ijk, 1:nsolid ) = rsfe_(1:nsolid) 
        rsfe( imjk, 1:nsolid ) = rsfw_(1:nsolid) 
        rsft( ijk, 1:nsolid ) = rsft_(1:nsolid)
        rsft( ijkm, 1:nsolid ) = rsfb_(1:nsolid) 

        rgfe( ijk ) = rgfe_
        rgfe( imjk ) = rgfw_ 
        rgft( ijk ) = rgft_
        rgft( ijkm ) = rgfb_
        
        IF(job_type == '3D') THEN
          rsfn( ijk, 1:nsolid ) = rsfn_(1:nsolid)
          rsfn( ijmk, 1:nsolid ) = rsfs_(1:nsolid)
          rgfn( ijk ) = rgfn_
          rgfn( ijmk ) = rgfs_
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

        USE control_flags, ONLY: job_type, lpr
        USE dimensions
        USE pressure_epsilon, ONLY: p, ep
        USE gas_solid_density, ONLY: rog, rgp, rgpn, rlk, rlkn
        USE gas_solid_temperature, ONLY: tg
        USE grid, ONLY: flag, noslip_wall
        USE eos_gas, ONLY: thermal_eosg, xgc
        USE phases_matrix, ONLY: matsvels_3phase
        USE set_indexes, ONLY: third_nb, third_rnb, first_rnb, first_nb
        USE set_indexes, ONLY: imjk, ijmk, ijkm
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
        REAL*8 :: rsfw1_ , rsfw2_
        REAL*8 :: rsft1_ , rsft2_
        REAL*8 :: rsfb1_ , rsfb2_
        REAL*8 :: rsfn1_ , rsfn2_
        REAL*8 :: rsfs1_ , rsfs2_

        REAL*8 :: rgfe_ , rgfw_
        REAL*8 :: rgft_ , rgfb_
        REAL*8 :: rgfn_ , rgfs_

        TYPE(stencil) :: ug_, vg_, wg_, rgp_
        TYPE(stencil) :: us1_, vs1_, ws1_, rlk1_
        TYPE(stencil) :: us2_, vs2_, ws2_, rlk2_

        INTEGER :: loop, kros, ig
        REAL*8 :: xgcl(max_ngas)
        xgcl(1:ngas) = xgc(ijk,:)
!
        kros = -1

        rsfe1_ = rsfe( ijk, 1 )
        rsfe2_ = rsfe( ijk, 2 )
        rsfw1_ = rsfe( imjk, 1 )
        rsfw2_ = rsfe( imjk, 2 )
        rsft1_ = rsft( ijk, 1 )
        rsft2_ = rsft( ijk, 2 )
        rsfb1_ = rsft( ijkm, 1 )
        rsfb2_ = rsft( ijkm, 2 )
        rsfn1_ = rsfn( ijk, 1 )
        rsfn2_ = rsfn( ijk, 2 )
        rsfs1_ = rsfn( ijmk, 1 )
        rsfs2_ = rsfn( ijmk, 2 )

        rgfe_ = rgfe( ijk )
        rgfw_ = rgfe( imjk )
        rgft_ = rgft( ijk )
        rgfb_ = rgft( ijkm )
        rgfn_ = rgfn( ijk )
        rgfs_ = rgfn( ijmk )

        CALL first_nb(rlk1_,rlk(:,1),ijk)
        CALL first_nb(rlk2_,rlk(:,2),ijk)
        CALL first_rnb(us1_,us(:,1),ijk)
        CALL first_rnb(us2_,us(:,2),ijk)
        CALL first_rnb(ws1_,ws(:,1),ijk)
        CALL first_rnb(ws2_,ws(:,2),ijk)
        CALL first_rnb(vs1_,vs(:,1),ijk)
        CALL first_rnb(vs2_,vs(:,2),ijk)
        CALL first_nb( rgp_, rgp, ijk )
        CALL first_rnb( ug_, ug, ijk )
        CALL first_rnb( wg_, wg, ijk )
        CALL first_rnb( vg_, vg, ijk )

        IF (muscl > 0) THEN
          CALL third_nb(rlk1_,rlk(:,1),ijk)
          CALL third_nb(rlk2_,rlk(:,2),ijk)
          CALL third_rnb(us1_,us(:,1),ijk)
          CALL third_rnb(us2_,us(:,2),ijk)
          CALL third_rnb(ws1_,ws(:,1),ijk)
          CALL third_rnb(ws2_,ws(:,2),ijk)
          CALL third_rnb(vs1_,vs(:,1),ijk)
          CALL third_rnb(vs2_,vs(:,2),ijk)
          CALL third_nb( rgp_, rgp, ijk )
          CALL third_rnb( ug_, ug, ijk )
          CALL third_rnb( wg_, wg, ijk )
          CALL third_rnb( vg_, vg, ijk )
        END IF
!
! ... Here start the inner loop
!
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
            CALL padjust( p(ijk), kros, d3, p3, omega, abeta_ )
          END IF
!
          ! ... Use equation of state to calculate gas density 
          ! ... from (new) pressure and (old) temperature.
!
          mg = 0.D0
          DO ig = 1, ngas
            mg = mg + xgcl(ig) * gmw( gas_type(ig) )
          END DO
          rog(ijk) = p(ijk) / ( rgas * tg(ijk) ) * mg
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
          us1_%c = us( ijk, 1 )
          us1_%w = us( imjk, 1 )
          us2_%c = us( ijk, 2 )
          us2_%w = us( imjk, 2 )
          ws1_%c = ws( ijk, 1 )
          ws1_%b = ws( ijkm, 1 )
          ws2_%c = ws( ijk, 2 )
          ws2_%b = ws( ijkm, 2 )
          vs1_%c = vs( ijk, 1 )
          vs1_%s = vs( ijmk, 1 )
          vs2_%c = vs( ijk, 2 )
          vs2_%s = vs( ijmk, 2 )

          ! ... update particle and gas densities and the 
          ! ... gas mass residual 'dg'

          CALL masf( rsfe1_,  rsfn1_,  rsft1_, rsfw1_, rsfs1_, rsfb1_,   &
                 rlk1_, us1_, vs1_, ws1_, ijk)

          IF (muscl > 0) THEN
            CALL fmas( rsfe1_,  rsfn1_,  rsft1_, rsfw1_, rsfs1_, rsfb1_, &
                 rlk1_, us1_, vs1_, ws1_, ijk)
          END IF

          rlkx = ( b_e * rsfe1_ - b_w * rsfw1_ ) * indx(i) * inr(i)
          rlky = ( b_n * rsfn1_ - b_s * rsfs1_ ) * indy(j)
          rlkz = ( b_t * rsft1_ - b_b * rsfb1_ ) * indz(k)
          rlk( ijk, 1) = rlkn( ijk, 1 ) - dt * ivf *  ( rlkx + rlky + rlkz )
          rlk( ijk, 1) = MAX( 0.0d0, rlk(ijk,1) )
          rls = rlk( ijk, 1 ) * inrl(1)

          CALL masf( rsfe2_,  rsfn2_,  rsft2_, rsfw2_, rsfs2_, rsfb2_,   &
                 rlk2_, us2_, vs2_, ws2_, ijk)

          IF (muscl > 0) THEN
            CALL fmas( rsfe2_,  rsfn2_,  rsft2_, rsfw2_, rsfs2_, rsfb2_, &
                 rlk2_, us2_, vs2_, ws2_, ijk)
          END IF

          rlkx = ( b_e * rsfe2_ - b_w * rsfw2_ ) * indx(i) * inr(i)
          rlky = ( b_n * rsfn2_ - b_s * rsfs2_ ) * indy(j)
          rlkz = ( b_t * rsft2_ - b_b * rsfb2_ ) * indz(k)
          rlk( ijk, 2) = rlkn( ijk, 2 ) - dt * ivf * ( rlkx + rlky + rlkz )
          rlk( ijk, 2) = MAX( 0.0d0, rlk(ijk,2) )
          rls = rls + rlk( ijk, 2 ) * inrl(2)

          IF( rls > 1.D0 ) THEN
            IF( lpr > 0 ) THEN
              WRITE(testunit,*) ' Mass is not conserved', flag(ijk)
              WRITE(testunit,*) ' i, j, k, rls ', i, j, k, rls
              WRITE(testunit,*) ' rls, volfrac ', rls, 1.D0/ivf
            ENDIF
            !IF (flag(ijkm) == noslip_wall) flag(ijk) = noslip_wall
            ierr = ierr + 1
          ENDIF

          ep( ijk ) = 1.D0 - rls
          rgp( ijk ) = ep( ijk ) * rog( ijk )
          rlk1_%c = rlk( ijk, 1 )
          rlk2_%c = rlk( ijk, 2 )
          rgp_%c = rgp(ijk)

          ! ... update residual of the Mass Balance equation of the gas phase
          ! ... using guessed velocities and pressure.

          ! ... update the particle volumetric fractions and the void fraction

          CALL masf( rgfe_ ,  rgfn_ ,  rgft_ , rgfw_ , rgfs_ , rgfb_ ,   &
                      rgp_, ug_, vg_, wg_, ijk )

          IF (muscl > 0) THEN
            CALL fmas( rgfe_ ,  rgfn_ ,  rgft_ , rgfw_ , rgfs_ , rgfb_ ,   &
                      rgp_, ug_, vg_, wg_, ijk )
          END IF

          resx = ( b_e * rgfe_ - b_w * rgfw_ ) * indx(i) * inr(i)
          resz = ( b_t * rgft_ - b_b * rgfb_ ) * indz(k)
          resy = ( b_n * rgfn_ - b_s * rgfs_ ) * indy(j)
          dg = rgp(ijk) - rgpn(ijk) + dt * ivf * ( resx + resy + resz )
!
          IF ( ( DABS(dg) > conv_ ) .OR. ( DABS(dg) >= DABS(dgorig) ) ) THEN

            cvg = .FALSE.
            IF( nit == 1 .AND. loop == 1 ) dgorig = dg
            d3 = dg
            ! ... steepen the Newton's slope (accelerate)
            IF( kros < 2 .AND. loop == inmax ) abeta_ = 0.5D0 * inmax * abeta_

          ELSE IF ( DABS( dg ) <= conv_ ) THEN

            cvg = .TRUE.
            EXIT inloop

          END IF

        END DO inloop
!
! ... Here rebuild the arrays
!
        rsfe( ijk, 1 ) = rsfe1_
        rsfe( ijk, 2 ) = rsfe2_
        rsfe( imjk, 1 ) = rsfw1_
        rsfe( imjk, 2 ) = rsfw2_
        rsft( ijk, 1 ) = rsft1_
        rsft( ijk, 2 ) = rsft2_
        rsft( ijkm, 1 ) = rsfb1_
        rsft( ijkm, 2 ) = rsfb2_
        rsfn( ijk, 1 ) = rsfn1_
        rsfn( ijk, 2 ) = rsfn2_
        rsfn( ijmk, 1 ) = rsfs1_
        rsfn( ijmk, 2 ) = rsfs2_

        rgfe( ijk ) = rgfe_
        rgfe( imjk ) = rgfw_
        rgft( ijk ) = rgft_
        rgft( ijkm ) = rgfb_
        rgfn( ijk ) = rgfn_
        rgfn( ijmk ) = rgfs_

        nloop = nloop + loop

        RETURN
      
      END SUBROUTINE opt3_inner_loop
!----------------------------------------------------------------------
      SUBROUTINE test_fluxes
      USE domain_mapping, ONLY: ncint, meshinds
      USE set_indexes, ONLY: subscr, imjk, ijmk, ijkm
      USE io_files, ONLY: tempunit
      IMPLICIT NONE
      INTEGER :: ijk, imesh, i, j, k

      OPEN(UNIT=tempunit,FILE='pdac.fl',STATUS='UNKNOWN')
      WRITE(tempunit,*) 'Test mass fluxes ...'
      WRITE(tempunit,*)
      DO ijk = 1, ncint
        CALL subscr(ijk)
        CALL meshinds(ijk,imesh,i,j,k)

        IF (i==11 .AND. k==51) THEN
                WRITE(tempunit,101) rgfe(ijk), rgft(ijk)
                WRITE(tempunit,101) rgfe(imjk), rgft(ijkm)
        END IF
        IF (i==10 .AND. k==52) THEN
                WRITE(tempunit,101) rgfe(ijk), rgft(ijk)
                WRITE(tempunit,101) rgfe(imjk), rgft(ijkm)
        END IF
        IF (i==11 .AND. k==52) THEN
                WRITE(tempunit,101) rgfe(ijk), rgft(ijk)
                WRITE(tempunit,101) rgfe(imjk), rgft(ijkm)
        END IF

      END DO

 100  FORMAT(4(I5), 3(G20.10E2))
 101  FORMAT(3(G20.10E2))
      CLOSE(tempunit)

      RETURN
      END SUBROUTINE test_fluxes
!----------------------------------------------------------------------
      SUBROUTINE correct_particles(ijk, imesh, i, j, k)
      USE control_flags, ONLY: job_type, lpr
      USE dimensions, ONLY: nsolid, ngas
      USE gas_constants, ONLY: tzero, hzeros
      USE gas_solid_density, ONLY: rlk
      USE gas_solid_temperature, ONLY: sies, ts, tg
      USE gas_solid_velocity, ONLY: ug, vg, wg
      USE gas_solid_velocity, ONLY: us, vs, ws
      USE particles_constants, ONLY: rl, inrl, cps
      USE specific_heat_module, ONLY: ck, hcaps
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ijk, imesh, i, j, k
      INTEGER :: ig, is
!
      DO is = 1, nsolid
        IF (rlk(ijk,is)*inrl(is) <= 1.D-10) THEN
          !... Drop out particles
          rlk(ijk,is) = 0.D0
          us(ijk,is)  = 0.D0
          IF (job_type == '3D') vs(ijk,is)  = 0.D0
          ws(ijk,is) = 0.D0
        ELSE
          ! ... Assume istantaneous momentum and thermal equilibrium
          us(ijk,is)  = ug(ijk)
          IF (job_type == '3D') vs(ijk,is)  = vg(ijk)
          ws(ijk,is) = wg(ijk)
        END IF
        ts(ijk,is) = tg(ijk)
        CALL hcaps(ck(is,ijk), cps(is), ts(ijk,is))
        sies(ijk,is) = ( ts(ijk,is) - tzero ) * ck(is,ijk) + hzeros
      END DO
!
      RETURN
      END SUBROUTINE correct_particles
!----------------------------------------------------------------------
      END MODULE iterative_solver
!----------------------------------------------------------------------
