
!#if defined __HPM
!#  include "/cineca/prod/hpm/include/f_hpm.h"
!#endif

!----------------------------------------------------------------------
      MODULE iterative_solver
!----------------------------------------------------------------------
      USE grid, ONLY: dx, dy, dz, indx, indy, indz, inx
      USE set_indexes, ONLY: stencil

      IMPLICIT NONE

      SAVE
!
! ... convective mass fluxes
!
      REAL*8, DIMENSION(:),   ALLOCATABLE :: rgfe, rgfn, rgft
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: rsfe, rsfn, rsft

      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: amats
!
      REAL*8  :: omega, dg
      INTEGER :: inmax, maxout

      INTEGER :: ncnt

      ! working stencils 
      TYPE(stencil) :: u, v, w, dens         

      PRIVATE :: u, v, w, dens
      PRIVATE :: ncnt
!

!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------



!----------------------------------------------------------------------
!
!
!  This is the main iterative solver routine
!
!
!----------------------------------------------------------------------
!
   SUBROUTINE iter
!
!----------------------------------------------------------------------

! ... This routine is the core of the iterative procedure to solve the
! ... mass-momentum and the interphase coupling
! ... (2D-3D-Compliant)
!
      USE convective_mass_fluxes, ONLY: fmas, masf
      USE dimensions
      USE domain_decomposition, ONLY: ncint, ncdom, data_exchange
      USE domain_decomposition, ONLY: myijk, meshinds
      USE enthalpy_matrix, ONLY: ftem
      USE environment, ONLY: timing
      USE eos_gas, ONLY: thermal_eosg, xgc, cg
      USE gas_solid_temperature, ONLY: tg
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      USE gas_solid_density, ONLY: rog, rgp, rgpn, rlk, rlkn
      USE grid, ONLY: fl_l
      USE indijk_module
      USE output_dump, ONLY: outp
      USE parallel, ONLY: mpime
      USE particles_constants, ONLY: rl, inrl
      USE phases_matrix, ONLY: assemble_matrix, assemble_all_matrix
      USE phases_matrix, ONLY: solve_velocities, solve_all_velocities
      USE phases_matrix, ONLY: matspre_3phase
      USE pressure_epsilon, ONLY: p, ep
      USE set_indexes, ONLY: nb, rnb, first_rnb, rnb_13, nb_13
      USE set_indexes, ONLY: subscr, subscr_iter,subscr_red, imjk, ijmk, ijkm, ijkn, ijks, ijke, ijkw, ijkt, ijkb, myinds
      USE specific_heat_module, ONLY: cp
      USE tilde_energy, ONLY: htilde
      USE tilde_momentum, ONLY: tilde, appu, appv, appw
      USE time_parameters, ONLY: time, dt, timestart, tpr
      USE control_flags, ONLY: job_type
      USE flux_limiters, ONLY: muscl
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
      INTEGER :: nit, nloop, mustit, ninner
      INTEGER :: i_, j_, k_, ijk_

      LOGICAL, ALLOCATABLE :: converge(:)
      LOGICAL :: use_opt_inner
      LOGICAL :: ldg
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
! ... An approximate value for gas density fluxes is guessed, 
! ... using estimates of velocity and pressure.
! ... The best esitmates for new velocities and pressure
! ... are their values at previous time-step.
!
      CALL data_exchange(p)
!
! ... Assemble and solve the explicit phase matrix 
! ... to update velocity fields. New velocities are
! ... biassed by the wrong (old) pressure field.
!
      DO ijk_ = 1, ncint
        IF ( fl_l( ijk_ ) == 1) THEN
          !CALL meshinds( ijk_ , imesh, i_ , j_ , k_ )
          !CALL subscr( ijk_ )
          CALL subscr_iter( ijk_ )
          CALL assemble_matrix( ijk_ )
          CALL solve_velocities( ijk_)
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
! ... Compute the gas mass fluxes using the guessed velocities 
!
      use_opt_inner = ( nsolid == 2 ) .AND. ( job_type == '3D' ) 
      IF( use_opt_inner ) THEN
        ALLOCATE( amats( 6, 6, ncint ) )
      END IF

      DO ijk_ = 1, ncint

        IF ( fl_l( ijk_ ) == 1 ) THEN
 
          !CALL meshinds( ijk_ , imesh, i_ , j_ , k_ )

              imesh = myijk( ip0_jp0_kp0_, ijk_ )

              IF ( job_type == '3D' ) THEN

                 k_ = ( imesh - 1 ) / ( nx*ny ) + 1
                 j_ = MOD( imesh - 1, nx*ny) / nx + 1
                 i_ = MOD( MOD( imesh - 1, nx*ny ), nx ) + 1

              ELSE IF ( job_type == '2D' ) THEN

                 k_ = ( imesh - 1 ) / nx + 1
                 i_ = MOD( ( imesh - 1 ), nx) + 1

              END IF

          !CALL subscr( ijk_ )
          !CALL subscr_red( ijk_ )

              IF( job_type == '2D' ) THEN
                 ijkm  = myijk( ip0_jp0_km1_, ijk_ )
                 imjk  = myijk( im1_jp0_kp0_, ijk_ )

                 ijke  = myinds(ip1_jp0_kp0_, ijk_ )
                 ijkt  = myinds(ip0_jp0_kp1_, ijk_ )
                 ijkw  = myinds(im1_jp0_kp0_, ijk_ )
                 ijkb  = myinds(ip0_jp0_km1_, ijk_ )

              ELSE IF( job_type == '3D' ) THEN

                 imjk   = myijk( im1_jp0_kp0_ , ijk_ )
                 ijmk   = myijk( ip0_jm1_kp0_ , ijk_ )
                 ijkm   = myijk( ip0_jp0_km1_ , ijk_ )

                 ijke = myinds( ip1_jp0_kp0_ , ijk_ )
                 ijkw = myinds( im1_jp0_kp0_ , ijk_ )
                 ijkn = myinds( ip0_jp1_kp0_ , ijk_ )
                 ijks = myinds( ip0_jm1_kp0_ , ijk_ )
                 ijkt = myinds( ip0_jp0_kp1_ , ijk_ )
                 ijkb = myinds( ip0_jp0_km1_ , ijk_ )

              END IF


          CALL calc_gas_mass_flux( i_, j_, k_, ijk_ )
!
          ! ... compute the derivative of the gas mass residual
          ! ... with respect to gas pressure
!
          CALL betas( conv( ijk_ ), abeta( ijk_ ), i_, j_, k_, ijk_ )

          IF( use_opt_inner ) THEN
             !  prepare matrix au1 au av1 av aw1 aw from appu, appv, appw  and store them in 
             !  amats
             CALL matspre_3phase( amats(:,1,ijk_), amats(:,2,ijk_), amats(:,3,ijk_), &
               amats(:,4,ijk_), amats(:,5,ijk_), amats(:,6,ijk_), i_, j_, k_, ijk_ )
          END IF
!
        END IF
      END DO
!
!/////////////////////////////////////////////////////////////////////
! ... Here Start the external iterative sweep.
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
         !         CALL tilde

         !   *** TIMING (Convergence in the ijk cell) ***
!
         IF( timing ) THEN
            !st0 = cclock_wall()
            !CALL f_hpmstart( 1, ' sor ' )
            ! call system_clock (st0,ratc)
         END IF  

         nloop  = 0
         ninner = 0
         ncnt   = 0

         mesh_loop: DO ijk_ = 1, ncint

           converge( ijk_ ) = .FALSE.

           IF( fl_l( ijk_ ) /= 1 ) converge( ijk_ ) = .TRUE.

           IF( fl_l( ijk_ ) == 1 ) THEN

              ! ...   Calculate global index i_ , j_ , k_
              !
              !  CALL meshinds( ijk_ , imesh, i_ , j_ , k_ )

              imesh = myijk( ip0_jp0_kp0_, ijk_ )

              IF ( job_type == '3D' ) THEN

                 k_ = ( imesh - 1 ) / ( nx*ny ) + 1
                 j_ = MOD( imesh - 1, nx*ny) / nx + 1
                 i_ = MOD( MOD( imesh - 1, nx*ny ), nx ) + 1

              ELSE IF ( job_type == '2D' ) THEN

                 k_ = ( imesh - 1 ) / nx + 1
                 i_ = MOD( ( imesh - 1 ), nx) + 1

              END IF

              ! ...   Calculate indexes of the neighbouring cells

              !CALL subscr( ijk_ )
              !CALL subscr_red( ijk_ )

              IF( job_type == '2D' ) THEN
                 ijkm  = myijk( ip0_jp0_km1_, ijk_ )
                 imjk  = myijk( im1_jp0_kp0_, ijk_ )

                 ijke  = myinds(ip1_jp0_kp0_, ijk_ )
                 ijkt  = myinds(ip0_jp0_kp1_, ijk_ )
                 ijkw  = myinds(im1_jp0_kp0_, ijk_ )
                 ijkb  = myinds(ip0_jp0_km1_, ijk_ )

              ELSE IF( job_type == '3D' ) THEN

                 imjk   = myijk( im1_jp0_kp0_ , ijk_ )
                 ijmk   = myijk( ip0_jm1_kp0_ , ijk_ )
                 ijkm   = myijk( ip0_jp0_km1_ , ijk_ )

                 ijke = myinds( ip1_jp0_kp0_ , ijk_ )
                 ijkw = myinds( im1_jp0_kp0_ , ijk_ )
                 ijkn = myinds( ip0_jp1_kp0_ , ijk_ )
                 ijks = myinds( ip0_jm1_kp0_ , ijk_ )
                 ijkt = myinds( ip0_jp0_kp1_ , ijk_ )
                 ijkb = myinds( ip0_jp0_km1_ , ijk_ )

              END IF

!
              ! ...   Compute the residual of the mass balance 
              ! ...   equation of the gas phase
!
              CALL calc_res( i_ , j_ , k_ , ijk_ , dg )

              dgorig = dg
              ldg    = .TRUE.

              !IF( nit == maxout ) &
              !WRITE(6,*) 'Pre ',imesh, i_ , j_ , k_, converge( ijk_ ), conv( ijk_ ), dg

              IF( ABS( dg ) <= conv( ijk_ ) ) THEN

                ncnt = ncnt + 1
              
                ! ... If the residual is lower then the prescribed limit
                ! ... compute the gas and particle densities and proceed 
                ! ... to next cell

                converge( ijk_ ) = .TRUE.

                ! ... Compute new ep(.)

                !CALL calc_eps( i_ , j_ , k_ , ijk_  )

                rls = 0.D0 

                IF (job_type == '2D') THEN

                  DO is = 1, nsolid

                    rlkx = (rsfe(ijk_,is) - rsfe(imjk,is)) * indx(i_) * inx(i_)
                    rlkz = (rsft(ijk_,is) - rsft(ijkm,is)) * indz(k_)
                    rlky = 0.D0

                    rlk( ijk_, is) = rlkn( ijk_,is ) - dt * ( rlkx + rlky + rlkz )
!                    - dt * (r1(ijk_)+r2(ijk_)+r3(ijk_)+r4(ijk_)+r5(ijk_))
                    rlk( ijk_, is) = MAX( 0.0d0, rlk(ijk_,is) )
                    rls = rls + rlk( ijk_,is) * inrl(is)

                  END DO

                ELSE IF (job_type == '3D') THEN

                  DO is = 1, nsolid

                    rlkx = (rsfe(ijk_,is) - rsfe(imjk,is)) * indx(i_) * inx(i_)
                    rlky = (rsfn(ijk_,is) - rsfn(ijmk,is)) * indy(j_)
                    rlkz = (rsft(ijk_,is) - rsft(ijkm,is)) * indz(k_)
          
                    rlk_tmp = rlkn( ijk_, is ) - dt * ( rlkx + rlky + rlkz )
!                    - dt * (r1(ijk_)+r2(ijk_)+r3(ijk_)+r4(ijk_)+r5(ijk_))
                    rlk_tmp= MAX( 0.0d0, rlk_tmp )
                    rls = rls + rlk_tmp * inrl(is)
                    rlk( ijk_, is) = rlk_tmp

                 END DO

               END IF
      
               IF( rls > 1.D0 ) THEN
                 WRITE(8,*) ' warning1: mass is not conserved'
                 WRITE(8,*) ' time, i, j, k, rls ', time, i_, j_, k_, rls
                 CALL error('iter', 'warning: mass is not conserved',1)
               ENDIF

               ep( ijk_ ) = 1.D0 - rls


               rgp( ijk_ ) = ep( ijk_ ) * rog( ijk_ )

               !IF( nit == maxout ) &
               !WRITE(6,*) 'IF1 ',imesh, i_ , j_ , k_, converge( ijk_ ), conv( ijk_ ), dg

               ldg = .FALSE.


             ELSE IF ( ABS(dg) > conv( ijk_ ) ) THEN
!
                ! ... If the residual is higher then the prescribed limit
                ! ... start the inner (in-cell) iterative loop
                ! ... to correct pressure and velocities.
!
!
                d3 = dg
                p3 = p( ijk_ )

                !IF( nit == maxout ) &
                !WRITE(6,*) 'IF2P ',imesh, i_ , j_ , k_, converge( ijk_ ), conv( ijk_ ), d3, p3
!
                IF( use_opt_inner ) THEN
                  IF( muscl == 0 ) THEN
                    CALL inner_3phase_1st( i_ , j_ , k_ , ijk_ , nit, d3, p3, &
                                     abeta( ijk_ ), conv( ijk_ ), dgorig, nloop )
                  ELSE
                    CALL inner_3phase( i_ , j_ , k_ , ijk_ , nit, d3, p3, &
                                     abeta( ijk_ ), conv( ijk_ ), dgorig, nloop )
                  END IF
                ELSE
                  CALL inner_loop( i_ , j_ , k_ , ijk_ , nit, d3, p3, &
                                   abeta( ijk_ ), conv( ijk_ ), dgorig )
                END IF

                ninner = ninner + 1

                !IF( nit == maxout ) &
                !WRITE(6,*) 'IF2D ',imesh, i_ , j_ , k_, converge( ijk_ ), conv( ijk_ ), d3, p3

                ldg = .FALSE.

              END IF

              IF( ldg ) THEN
                WRITE(6,*) 'WARNING (iter) something wrong with dg '
                WRITE(6,*) '       ',dg,imesh, i_ , j_ , k_, converge( ijk_ ), conv( ijk_ )
              END IF
   
           END IF

         END DO mesh_loop

         IF( timing ) THEN
             !st1 = cclock_wall()
             ! call system_clock(st1,ratc)
             !CALL f_hpmstop( 1 )
         END IF


         timconv(nit) = timconv(nit) + real(st1-st0)         
         IF( ninner > 0 ) THEN
           WRITE(6, fmt="( I10, 2X, F4.2, 2X, I10, 2X, F6.3)" ) &
             ninner, REAL(nloop) / REAL(ninner), ncnt, timconv(nit)
         ELSE
           WRITE(6, fmt="( I10, 2X, F4.2, 2X, I10, 2X, F6.3)" ) &
             ninner, REAL(nloop), ncnt, timconv(nit)
         END IF
    
!
! ... Update gas and particles enthalpy
! ... (for fully implicit solution)
!
!         CALL htilde
!         CALL ftem

! ... Exchange all updated physical quantities on boundaries
! ... before starting a new sweep on the mesh.
!
        CALL data_exchange(rgp)
        CALL data_exchange(rlk)
        CALL data_exchange(p)
        CALL data_exchange(ep)
!
!********************************************************************
! ... check the convergence parameters and the convergence history
!
!      IF (MOD((time-timestart),tpr) <= dt) THEN
!        IF (nit == 1) THEN
!          WRITE(6,500) mpime, time
!          WRITE(6,501)
!        END IF
!        WRITE(6,502) nit, ALL(converge)
! 500    FORMAT('proc: ', I3, ' time: ', F8.3)
! 501    FORMAT('iteration # :  convergence')
! 502    FORMAT(I3,13X,L1)
!      END IF
!*******************************************************************
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
          omega = omega0
!*******************************************************************
! ... write out the final number of iterations
!
          WRITE(6,277) nit
 277      FORMAT('  from iter: nit = ', I4)
!*******************************************************************
          EXIT sor_loop
        ENDIF
!
! ... If convergence is not reached in some cell
! ... start a new external sweep.
!
       IF( MOD( nit, 100 ) == 0 ) THEN
         omega = omega * 0.5D0
         WRITE(6, fmt="('  reducing relaxation parameter omega')")
         WRITE(6, fmt="('  new value = ',F12.4)") omega
         call myflush( 6 )
       END IF
!
      END DO sor_loop

      IF( use_opt_inner ) THEN
        DEALLOCATE( amats )
      END IF

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
        WRITE(6,700) nit, (time+dt)
        WRITE(6,*) 'convergence on proc ',mpime,' : ', ALL(converge)
        IF (.NOT.ALL(converge)) WRITE(6,*) 'cells not converged (imesh,i,j,k): '
        DO ijk_ = 1, ncint
          IF ( .NOT. converge( ijk_ ) ) THEN
            CALL meshinds( ijk_ , imesh, i_ , j_ , k_ )
            WRITE(6,*) imesh, i_ , j_ , k_, conv( ijk_ )
          END IF
        END DO
 700    FORMAT('max number of iterations (',I5,') reached at time: ', F8.3)

        CALL error( ' iter ', 'max number of iters exceeded ', 1)
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
!
      RETURN

!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------

      SUBROUTINE test_fluxes
!
      USE control_flags, ONLY: job_type
      USE dimensions, ONLY: nsolid
      USE io_restart, ONLY: write_array
      USE kinds
      USE parallel, ONLY: nproc, mpime, root, group
      USE time_parameters, ONLY: time
!
      IMPLICIT NONE
!
      CHARACTER :: filnam*12
!
      INTEGER :: ig,is
      LOGICAL :: lform = .TRUE.
!
      filnam='output.mtest'

      IF( mpime == root ) THEN

        IF (lform) THEN
          OPEN(UNIT=12,FILE=filnam)
          WRITE(12,*) time
        ELSE 
          OPEN(UNIT=12,FORM='UNFORMATTED',FILE=filnam)
          WRITE(12) REAL(time,4)
        END IF
 
      END IF
!
      IF (job_type == '2D') THEN
        CALL write_array( 12, rgfe, sgl, lform )
        CALL write_array( 12, rgft, sgl, lform )
      ELSE IF (job_type == '3D') THEN
        CALL write_array( 12, rgfe, sgl, lform )
        CALL write_array( 12, rgfn, sgl, lform )
        CALL write_array( 12, rgft, sgl, lform )
      END IF
!
      DO is = 1, nsolid
        IF (job_type == '2D') THEN
          CALL write_array( 12, rsfe(:,is), sgl, lform )
          CALL write_array( 12, rsft(:,is), sgl, lform )
        ELSE IF (job_type == '3D') THEN
          CALL write_array( 12, rsfe(:,is), sgl, lform )
          CALL write_array( 12, rsfn(:,is), sgl, lform )
          CALL write_array( 12, rsft(:,is), sgl, lform )
        END IF
      END DO

      IF( mpime == root ) THEN
        CLOSE (12)
      END IF
!
      RETURN
      END SUBROUTINE test_fluxes

!----------------------------------------------------------------------
   END SUBROUTINE iter
!----------------------------------------------------------------------




!----------------------------------------------------------------------
!
!
!
      SUBROUTINE calc_gas_mass_flux( i, j, k, ijk )
!
!----------------------------------------------------------------------
!
! ... compute the gas mass fluxes by using the stencil, as defined in the
! ... 'subscr' module
!
      USE dimensions
      USE flux_limiters, ONLY: muscl
      USE convective_mass_fluxes, ONLY: fmas, masf
      USE gas_solid_velocity, ONLY: ug, vg, wg
      USE gas_solid_density, ONLY: rgp
      USE set_indexes, ONLY: nb, rnb, first_rnb, rnb_13, nb_13, first_nb
      USE set_indexes, ONLY: subscr, imjk, ijmk, ijkm
      USE set_indexes, ONLY: chkstencil
      USE control_flags, ONLY: job_type
!
      IMPLICIT NONE

      INTEGER :: i, j, k, ijk, info

  

        IF (job_type == '2D') THEN
        
          CALL nb_13 ( dens, rgp, ijk )
          CALL rnb_13( u, ug, ijk )
          CALL rnb_13( w, wg, ijk )

          CALL fmas( rgfe( ijk ),  rgft( ijk ), rgfe( imjk ), rgft( ijkm ),   &
                    dens, u, w, ijk )

        ELSE IF (job_type == '3D') THEN

          IF ( muscl == 0 ) THEN

            CALL first_nb ( dens, rgp, ijk )
            CALL first_rnb( u, ug, ijk )
            CALL first_rnb( w, wg, ijk )
            CALL first_rnb( v, vg, ijk )

            ! CALL chkstencil( dens, info )
            ! IF( info /= 0 ) WRITE(6,*) 'wrong dens at: ', i, j, k
            ! CALL chkstencil( u, info )
            ! IF( info /= 0 ) WRITE(6,*) 'wrong u at: ', i, j, k
            ! CALL chkstencil( v, info )
            ! IF( info /= 0 ) WRITE(6,*) 'wrong v at: ', i, j, k
            ! CALL chkstencil( w, info )
            ! IF( info /= 0 ) WRITE(6,*) 'wrong w at: ', i, j, k

            CALL masf( rgfe( ijk  ),  rgfn( ijk  ),  rgft( ijk  ),    &
                       rgfe( imjk ),  rgfn( ijmk ),  rgft( ijkm ),    &
                    dens, u, v, w, i, j, k, ijk  )

          ELSE

            CALL nb_13 ( dens, rgp, ijk )
            CALL rnb_13( u, ug, ijk )
            CALL rnb_13( w, wg, ijk )
            CALL rnb_13( v, vg, ijk )
            CALL fmas( rgfe( ijk  ),  rgfn( ijk  ),  rgft( ijk  ),    &
                       rgfe( imjk ),  rgfn( ijmk ),  rgft( ijkm ),   &
                    dens, u, v, w, i, j, k, ijk  )
          END IF

        ENDIF
      END SUBROUTINE calc_gas_mass_flux   


!----------------------------------------------------------------------
!
!
!
      SUBROUTINE calc_eps( i, j, k, ijk )
!
!----------------------------------------------------------------------

      USE dimensions
      USE gas_solid_density, ONLY: rlk, rlkn
      USE pressure_epsilon, ONLY: p, ep
      USE set_indexes, ONLY: subscr, imjk, ijmk, ijkm
      USE time_parameters, ONLY: dt, time
      USE control_flags, ONLY: job_type
      USE particles_constants, ONLY: inrl

! ... Use the Mass Balance equation of the solids
! ... to compute particle volumetric fractions and the void fraction.
!
      IMPLICIT NONE
      REAL*8 :: rlkx, rlky, rlkz, rls, rlk_tmp
      INTEGER :: i, j, k, ijk
      INTEGER :: is

      rls = 0.D0

      IF (job_type == '2D') THEN

        DO is = 1, nsolid

          rlkx = (rsfe(ijk,is) - rsfe(imjk,is)) * indx(i) * inx(i)
          rlkz = (rsft(ijk,is) - rsft(ijkm,is)) * indz(k)
          rlky = 0.D0

          rlk( ijk, is) = rlkn( ijk,is ) - dt * ( rlkx + rlky + rlkz )
!          - dt * (r1(ijk)+r2(ijk)+r3(ijk)+r4(ijk)+r5(ijk))
          rlk( ijk, is) = MAX( 0.0d0, rlk(ijk,is) )
          rls = rls + rlk( ijk,is) * inrl(is)

        END DO

      ELSE IF (job_type == '3D') THEN

        DO is = 1, nsolid

          rlkx = (rsfe(ijk,is) - rsfe(imjk,is)) * indx(i) * inx(i)
          rlky = (rsfn(ijk,is) - rsfn(ijmk,is)) * indy(j)
          rlkz = (rsft(ijk,is) - rsft(ijkm,is)) * indz(k)
           
          !rlk( ijk, is) = rlkn( ijk, is ) - dt * ( rlkx + rlky + rlkz )
!          - dt * (r1(ijk)+r2(ijk)+r3(ijk)+r4(ijk)+r5(ijk))
          !rlk( ijk, is) = MAX( 0.0d0, rlk(ijk,is) )
          !rls = rls + rlk( ijk, is ) * inrl(is)
          
          rlk_tmp = rlkn( ijk, is ) - dt * ( rlkx + rlky + rlkz )
!          - dt * (r1(ijk)+r2(ijk)+r3(ijk)+r4(ijk)+r5(ijk))
          rlk_tmp= MAX( 0.0d0, rlk_tmp )
          rls = rls + rlk_tmp * inrl(is)
          rlk( ijk, is) = rlk_tmp

        END DO

      END IF
      
      IF( rls > 1.D0 ) THEN
        WRITE(8,*) ' warning1: mass is not conserved'
        WRITE(8,*) ' time, i, j, k, rls ', time, i, j, k, rls
        CALL error('iter', 'warning: mass is not conserved',1)
      ENDIF

      ep( ijk ) = 1.D0 - rls

      END SUBROUTINE calc_eps

!----------------------------------------------------------------------
!
!
!
      SUBROUTINE calc_res( i, j, k, ijk, res )
!
!----------------------------------------------------------------------

! ... Compute the residual of the Mass Balance equation of the gas phase
!
      USE set_indexes, ONLY: subscr, imjk, ijmk, ijkm
      USE gas_solid_density, ONLY: rog, rgp, rgpn, rlk, rlkn
      USE control_flags, ONLY: job_type
      USE time_parameters, ONLY: dt

      IMPLICIT NONE
      REAL*8 :: resx, resy, resz
      REAL*8, INTENT(OUT) :: res
      INTEGER :: i, j, k, ijk

      INTEGER :: float_chk
!      
      resy = 0.D0
!
      resx = ( rgfe(ijk) - rgfe(imjk) ) * indx(i) * inx(i)
      resz = ( rgft(ijk) - rgft(ijkm) ) * indz(k)
      IF (job_type == '3D') resy = (rgfn(ijk) - rgfn(ijmk)) * indy(j)
!
      res  = rgp(ijk) - rgpn(ijk) + dt * (resx+resy+resz)
!          - dt * (r1(ijk)+r2(ijk)+r3(ijk)+r4(ijk)+r5(ijk))

      RETURN
      END SUBROUTINE calc_res

!----------------------------------------------------------------------
!
!
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
      REAL*8 FUNCTION newp(d1, d2, d3, p1, p2, p3)
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
!
!
      SUBROUTINE betas(cnv, abt, i, j, k, ijk)
!
!----------------------------------------------------------------------
! 
      USE dimensions
      USE domain_decomposition, ONLY: ncint, myijk, meshinds
      USE convective_mass_fluxes, ONLY: upc_e, upc_n, upc_t
      USE convective_mass_fluxes, ONLY: upc_w, upc_s, upc_b
      USE control_flags, ONLY: job_type
      USE eos_gas, ONLY: csound
      USE gas_solid_density, ONLY: rog, rgp
      USE grid, ONLY: dx, dy, dz, fl_l, xb
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE pressure_epsilon, ONLY: p, ep
      USE set_indexes
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
      INTEGER, INTENT(IN) :: i, j, k, ijk
      REAL*8, INTENT(OUT) :: cnv, abt
!
      REAL*8, PARAMETER :: delg=1.D-8
!
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
          nfle=fl_l(ipjk)
          nflw=fl_l(imjk)
          nflt=fl_l(ijkp)
          nflb=fl_l(ijkm)

          IF( (nfle /= 1) .AND. (nfle /= 4) ) THEN
            iep_e = 0.D0
          ELSE
            iep_e = ( dx(i+1)*ep(ijk)+dx(i)*ep(ijke))*indxp*indxp*2.D0
            iep_e = iep_e * upc_e
          END IF

          IF( (nflw /= 1) .AND. (nflw /= 4) ) THEN
            iep_w = 0.D0
          ELSE
            iep_w = ( dx(i-1)*ep(ijk)+dx(i)*ep(ijkw) )*indxm*indxm*2.D0 
            iep_w = iep_w * upc_w
          END IF
!
          IF( (nflt /= 1) .AND. (nflt /= 4) ) THEN 
            iep_t = 0.0D0
          ELSE
            iep_t = ( dz(k+1)*ep(ijk)+dz(k)*ep(ijkt) )*indzp*indzp*2.D0 
            iep_t = iep_t * upc_t 
          END IF

          IF( (nflb /= 1) .AND. (nflb /= 4) ) THEN
            iep_b = 0.D0
          ELSE
            iep_b = ( dz(k-1)*ep(ijk)+dz(k)*ep(ijkb) )*indzm*indzm*2.D0
            iep_b = iep_b * upc_b
          END IF
!
          IF (job_type == '3D') THEN

            dyp=dy(j)+dy(j+1)
            dym=dy(j)+dy(j-1)
            indyp=1.D0/dyp
            indym=1.D0/dym
            nfln=fl_l(ijpk)
            nfls=fl_l(ijmk)

            IF( (nfln /= 1) .AND. (nfln /= 4) ) THEN 
              iep_n = 0.0D0
            ELSE
              iep_n = ( dy(j+1)*ep(ijk)+dy(j)*ep(ijkn) )*indyp*indyp*2.D0 
              iep_n = iep_n * upc_n 
            END IF

            IF( (nfls /= 1) .AND. (nfls /= 4) ) THEN
              iep_s = 0.D0
            ELSE
              iep_s = ( dy(j-1)*ep(ijk)+dy(j)*ep(ijks) )*indym*indym*2.D0
              iep_s = iep_s * upc_s
            END IF

          END IF

          iepx = (  xb(i) * iep_e + xb(i-1) * iep_w ) * indx(i) * inx(i) 
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
             
!          CALL csound(rags,rog(ijk),p(ijk))
!
! ... rbeta = dD_g/dP

          rbeta = ep(ijk) * rags + dt**2 * (iepx+iepy+iepz)
!
          abt = 1.D0 / rbeta
          cnv = delg * rgp(ijk)
!
      RETURN
      END SUBROUTINE betas

!----------------------------------------------------------------------
!
!
!
      SUBROUTINE inner_3phase_1st( i, j, k, ijk, nit, d3, p3, abeta_, conv_, dgorig, nloop )
!
!----------------------------------------------------------------------

        USE dimensions
        USE pressure_epsilon, ONLY: p, ep
        USE gas_solid_density, ONLY: rog, rgp, rgpn, rlk, rlkn
        USE gas_solid_temperature, ONLY: tg
        USE eos_gas, ONLY: thermal_eosg, xgc, cg
        USE phases_matrix, ONLY: matsvels_3phase
        USE set_indexes, ONLY: imjk, ijmk, ijkm
        USE set_indexes, ONLY: ijke, ijkn, ijkt, ijkw, ijks, ijkb
        USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
        USE control_flags, ONLY: job_type
        USE time_parameters, ONLY: dt
        USE convective_mass_fluxes, ONLY: fmas, masf
        USE gas_constants, ONLY: gmw, rgas, gas_type
        USE particles_constants, ONLY: inrl

        IMPLICIT NONE

        INTEGER :: nit, ijk, i, j, k, nloop
        REAL*8 :: d3, p3, abeta_, conv_, dgorig
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

!
        mg = 0.D0
        DO ig = 1, ngas
          mg = mg + xgc( ig, ijk ) * gmw( gas_type(ig) )
        END DO

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
          ! ... update particle and gas densities and the 
          ! ... gas mass residual 'dg'

          ! ... update the particle volumetric fractions and the void fraction

          IF( loop == 1 ) THEN

               ! CALL first_nb(dens1,rlk(:,1),ijk)
               dens1%c = rlk( ijk,  1 )
               dens1%e = rlk( ijke, 1 )
               dens1%w = rlk( ijkw, 1 )
               dens1%n = rlk( ijkn, 1 )
               dens1%s = rlk( ijks, 1 )
               dens1%t = rlk( ijkt, 1 )
               dens1%b = rlk( ijkb, 1 )

               ! CALL first_nb(dens2,rlk(:,2),ijk)
               dens2%c = rlk( ijk,  2 )
               dens2%e = rlk( ijke, 2 )
               dens2%w = rlk( ijkw, 2 )
               dens2%n = rlk( ijkn, 2 )
               dens2%s = rlk( ijks, 2 )
               dens2%t = rlk( ijkt, 2 )
               dens2%b = rlk( ijkb, 2 )

               ! CALL first_rnb(u1,us(:,1),ijk)
               ! CALL first_rnb(u2,us(:,2),ijk)
               ! CALL first_rnb(w1,ws(:,1),ijk)
               ! CALL first_rnb(w2,ws(:,2),ijk)
               ! CALL first_rnb(v1,vs(:,1),ijk)
               ! CALL first_rnb(v2,vs(:,2),ijk)
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

          ELSE

             dens1%c = rlk( ijk, 1 )
             dens2%c = rlk( ijk, 2 )
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

          END IF

             ! ... West, South and Bottom fluxes
             !
             IF (u1%w >= 0.D0) THEN
                rsfem1_ = u1%w * dens1%w
             ELSE
                rsfem1_ = u1%w * dens1%c
             ENDIF

             IF (v1%s >= 0.D0) THEN
                rsfnm1_ = v1%s * dens1%s
             ELSE
                rsfnm1_ = v1%s * dens1%c  
             ENDIF

             IF (w1%b >= 0.D0) THEN
                rsftm1_ = w1%b * dens1%b
             ELSE
                rsftm1_ = w1%b * dens1%c
             ENDIF

             ! ... East, North and Top fluxes
             !
             IF (u1%c >= 0.D0) THEN
                rsfe1_ = u1%c * dens1%c
             ELSE
                rsfe1_ = u1%c * dens1%e
             ENDIF

             IF (v1%c >= 0.D0) THEN
                rsfn1_ = v1%c * dens1%c
             ELSE
                rsfn1_ = v1%c * dens1%n
             ENDIF

             IF (w1%c >= 0.D0) THEN
                rsft1_ = w1%c * dens1%c
             ELSE
                rsft1_ = w1%c * dens1%t
             ENDIF

            !CALL masf( rsfe1_ ,  rsfn1_ ,  rsft1_ , rsfem1_ , rsfnm1_ , rsftm1_ ,   &
            !     dens1, u1, v1, w1 )

          rlkx = ( rsfe1_ - rsfem1_ ) * indx(i) * inx(i)
          rlky = ( rsfn1_ - rsfnm1_ ) * indy(j)
          rlkz = ( rsft1_ - rsftm1_ ) * indz(k)
          rlk( ijk, 1) = rlkn( ijk, 1 ) - dt * ( rlkx + rlky + rlkz )
          rlk( ijk, 1) = MAX( 0.0d0, rlk(ijk,1) )
          rls = rlk( ijk, 1 ) * inrl(1)


             ! ... West, South and Bottom fluxes
             !
             IF (u2%w >= 0.D0) THEN
                rsfem2_ = u2%w * dens2%w
             ELSE
                rsfem2_ = u2%w * dens2%c
             ENDIF

             IF (v2%s >= 0.D0) THEN
                rsfnm2_ = v2%s * dens2%s
             ELSE
                rsfnm2_ = v2%s * dens2%c  
             ENDIF

             IF (w2%b >= 0.D0) THEN
                rsftm2_ = w2%b * dens2%b
             ELSE
                rsftm2_ = w2%b * dens2%c
             ENDIF

             ! ... East, North and Top fluxes
             !
             IF (u2%c >= 0.D0) THEN
                rsfe2_ = u2%c * dens2%c
             ELSE
                rsfe2_ = u2%c * dens2%e
             ENDIF

             IF (v2%c >= 0.D0) THEN
                rsfn2_ = v2%c * dens2%c
             ELSE
                rsfn2_ = v2%c * dens2%n
             ENDIF

             IF (w2%c >= 0.D0) THEN
                rsft2_ = w2%c * dens2%c
             ELSE
                rsft2_ = w2%c * dens2%t
             ENDIF

            !CALL masf( rsfe2_ ,  rsfn2_ ,  rsft2_ , rsfem2_ , rsfnm2_ , rsftm2_ ,   &
            !     dens2, u2, v2, w2 )

          rlkx = ( rsfe2_ - rsfem2_ ) * indx(i) * inx(i)
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

          ! ... update residual of the Mass Balance equation of the gas phase
          ! ... using guessed velocities and pressure.
          ! CALL update_res(dg)

          IF( loop == 1 ) THEN

               ! CALL first_nb( densg_, rgp, ijk )
               densg_%c = rgp( ijk  )
               densg_%e = rgp( ijke )
               densg_%w = rgp( ijkw )
               densg_%n = rgp( ijkn )
               densg_%s = rgp( ijks )
               densg_%t = rgp( ijkt )
               densg_%b = rgp( ijkb )

               ! CALL first_rnb( ug_, ug, ijk )
               ! CALL first_rnb( wg_, wg, ijk )
               ! CALL first_rnb( vg_, vg, ijk )

               ug_%c = ug( ijk )
               ug_%w = ug( imjk )
               wg_%c = wg( ijk )
               wg_%b = wg( ijkm )
               vg_%c = vg( ijk )
               vg_%s = vg( ijmk )

          ELSE

             densg_%c = rgp(ijk)
             ug_%c = ug( ijk )
             ug_%w = ug( imjk )
             wg_%c = wg( ijk )
             wg_%b = wg( ijkm )
             vg_%c = vg( ijk )
             vg_%s = vg( ijmk )

          END IF

             
             ! ... West, South and Bottom fluxes
             !
             IF (ug_%w >= 0.D0) THEN
                rgfem_ = ug_%w * densg_%w
             ELSE
                rgfem_ = ug_%w * densg_%c
             ENDIF

             IF (vg_%s >= 0.D0) THEN
                rgfnm_ = vg_%s * densg_%s
             ELSE
                rgfnm_ = vg_%s * densg_%c  
             ENDIF

             IF (wg_%b >= 0.D0) THEN
                rgftm_ = wg_%b * densg_%b
             ELSE
                rgftm_ = wg_%b * densg_%c
             ENDIF

             ! ... East, North and Top fluxes
             !
             IF (ug_%c >= 0.D0) THEN
               rgfe_ = ug_%c * densg_%c
             ELSE
               rgfe_ = ug_%c * densg_%e
             ENDIF

             IF (vg_%c >= 0.D0) THEN
               rgfn_ = vg_%c * densg_%c
             ELSE
               rgfn_ = vg_%c * densg_%n
             ENDIF

             IF (wg_%c >= 0.D0) THEN
               rgft_ = wg_%c * densg_%c
             ELSE
               rgft_ = wg_%c * densg_%t
             ENDIF

          resx = ( rgfe_ - rgfem_ ) * indx(i) * inx(i)
          resz = ( rgft_ - rgftm_ ) * indz(k)
          resy = ( rgfn_ - rgfnm_ ) * indy(j)
          dg = rgp(ijk) - rgpn(ijk) + dt * ( resx + resy + resz )
!
          IF ( ( DABS(dg) > conv_ ) .OR. ( DABS(dg) >= DABS(dgorig) ) ) THEN

            IF( nit == 1 .AND. loop == 1 ) dgorig = dg
            d3 = dg
            ! ... steepen the Newton's slope (accelerate)
            IF( kros < 2 .AND. loop == inmax ) abeta_ = 0.5D0 * inmax * abeta_

          ELSE IF ( DABS( dg ) <= conv_ ) THEN

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
      
      END SUBROUTINE inner_3phase_1st

!----------------------------------------------------------------------
!
!
!
      SUBROUTINE inner_3phase( i, j, k, ijk, nit, d3, p3, abeta_, conv_, dgorig, nloop )
!
!----------------------------------------------------------------------

        USE dimensions
        USE pressure_epsilon, ONLY: p, ep
        USE gas_solid_density, ONLY: rog, rgp, rgpn, rlk, rlkn
        USE gas_solid_temperature, ONLY: tg
        USE eos_gas, ONLY: thermal_eosg, xgc, cg
        USE phases_matrix, ONLY: matsvels_3phase
        USE set_indexes, ONLY: nb, rnb, first_rnb, rnb_13, nb_13, first_nb
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
          ! ... update particle and gas densities and the 
          ! ... gas mass residual 'dg'

          ! ... update the particle volumetric fractions and the void fraction

          IF( loop == 1 ) THEN

             CALL nb_13(dens1,rlk(:,1),ijk)
             CALL nb_13(dens2,rlk(:,2),ijk)
             CALL rnb_13(u1,us(:,1),ijk)
             CALL rnb_13(u2,us(:,2),ijk)
             CALL rnb_13(w1,ws(:,1),ijk)
             CALL rnb_13(w2,ws(:,2),ijk)
             CALL rnb_13(v1,vs(:,1),ijk)
             CALL rnb_13(v2,vs(:,2),ijk)
         
          ELSE

             dens1%c = rlk( ijk, 1 )
             dens2%c = rlk( ijk, 2 )
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

          END IF

          CALL fmas( rsfe1_ ,  rsfn1_ ,  rsft1_ , rsfem1_ , rsfnm1_ , rsftm1_ ,   &
                 dens1, u1, v1, w1, i, j, k, ijk)

          rlkx = ( rsfe1_ - rsfem1_ ) * indx(i) * inx(i)
          rlky = ( rsfn1_ - rsfnm1_ ) * indy(j)
          rlkz = ( rsft1_ - rsftm1_ ) * indz(k)
          rlk( ijk, 1) = rlkn( ijk, 1 ) - dt * ( rlkx + rlky + rlkz )
          rlk( ijk, 1) = MAX( 0.0d0, rlk(ijk,1) )
          rls = rlk( ijk, 1 ) * inrl(1)

          CALL fmas( rsfe2_ ,  rsfn2_ ,  rsft2_ , rsfem2_ , rsfnm2_ , rsftm2_ ,   &
                 dens2, u2, v2, w2, i, j, k, ijk)

          rlkx = ( rsfe2_ - rsfem2_ ) * indx(i) * inx(i)
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

          ! ... update residual of the Mass Balance equation of the gas phase
          ! ... using guessed velocities and pressure.
          ! CALL update_res(dg)

          IF( loop == 1 ) THEN
               CALL nb_13( densg_, rgp, ijk )
               CALL rnb_13( ug_, ug, ijk )
               CALL rnb_13( wg_, wg, ijk )
               CALL rnb_13( vg_, vg, ijk )
          ELSE
             densg_%c = rgp(ijk)
             ug_%c = ug( ijk )
             ug_%w = ug( imjk )
             wg_%c = wg( ijk )
             wg_%b = wg( ijkm )
             vg_%c = vg( ijk )
             vg_%s = vg( ijmk )
          END IF

          CALL fmas( rgfe_ ,  rgfn_ ,  rgft_ , rgfem_ , rgfnm_ , rgftm_ ,   &
                      densg_, ug_, vg_, wg_, i, j, k, ijk )

          resx = ( rgfe_ - rgfem_ ) * indx(i) * inx(i)
          resz = ( rgft_ - rgftm_ ) * indz(k)
          resy = ( rgfn_ - rgfnm_ ) * indy(j)
          dg = rgp(ijk) - rgpn(ijk) + dt * ( resx + resy + resz )
!
          IF ( ( DABS(dg) > conv_ ) .OR. ( DABS(dg) >= DABS(dgorig) ) ) THEN

            IF( nit == 1 .AND. loop == 1 ) dgorig = dg
            d3 = dg
            ! ... steepen the Newton's slope (accelerate)
            IF( kros < 2 .AND. loop == inmax ) abeta_ = 0.5D0 * inmax * abeta_

          ELSE IF ( DABS( dg ) <= conv_ ) THEN

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
      
      END SUBROUTINE inner_3phase


!----------------------------------------------------------------------
!
!
!
      SUBROUTINE inner_loop( i, j, k, ijk, nit, d3, p3, abeta_, conv_, dgorig )
!
!----------------------------------------------------------------------

        USE dimensions
        USE pressure_epsilon, ONLY: p, ep
        USE gas_solid_density, ONLY: rog, rgp, rgpn, rlk, rlkn
        USE gas_solid_temperature, ONLY: tg
        USE eos_gas, ONLY: thermal_eosg, xgc, cg
        USE phases_matrix, ONLY: assemble_matrix, assemble_all_matrix
        USE phases_matrix, ONLY: solve_velocities, solve_all_velocities
        USE set_indexes, ONLY: nb, rnb, first_rnb, rnb_13, nb_13, first_nb
        USE set_indexes, ONLY: subscr, imjk, ijmk, ijkm
        USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
        USE control_flags, ONLY: job_type
        USE time_parameters, ONLY: dt
        USE flux_limiters, ONLY: muscl
        USE convective_mass_fluxes, ONLY: fmas, masf

        IMPLICIT NONE

        INTEGER :: nit, ijk, i, j, k
        REAL*8 :: d3, p3, abeta_, conv_, dgorig
        REAL*8 :: resx, resy, resz

        INTEGER :: loop, kros

        TYPE(stencil) :: ug_, vg_, wg_, densg_

        kros = -1

        IF( muscl == 0 ) THEN
          CALL first_nb( densg_, rgp, ijk )
          CALL first_rnb( ug_, ug, ijk )
          CALL first_rnb( wg_, wg, ijk )
          IF( job_type == '3D' ) THEN
            CALL first_rnb( vg_, vg, ijk )
          END IF
        ELSE
          CALL nb_13( densg_, rgp, ijk )
          CALL rnb_13( ug_, ug, ijk )
          CALL rnb_13( wg_, wg, ijk )
          IF( job_type == '3D' ) THEN
            CALL rnb_13( vg_, vg, ijk )
          END IF
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
            CALL padjust(p(ijk), kros, d3, p3, omega, abeta_ )
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
          ! ... update particle and gas densities and the 
          ! ... gas mass residual 'dg'

          ! ... update the particle volumetric fractions and the void fraction
          ! CALL update_eps

          CALL calc_part_mass_flux
          CALL calc_eps( i, j, k, ijk )

          rgp( ijk ) = ep( ijk ) * rog( ijk )

          ! ... update residual of the Mass Balance equation of the gas phase
          ! ... using guessed velocities and pressure.
          ! CALL update_res(dg)

          densg_%c = rgp(ijk)
          ug_%c = ug( ijk )
          ug_%w = ug( imjk )
          wg_%c = wg( ijk )
          wg_%b = wg( ijkm )
          IF( job_type == '3D' ) THEN
             vg_%c = vg( ijk )
             vg_%s = vg( ijmk )
          END IF

          ! --- by hand inlining -- CALL calc_gas_mass_flux
          IF (job_type == '2D') THEN
        
             CALL fmas(rgfe(ijk),  rgft(ijk), rgfe(imjk), rgft(ijkm),   &
                      densg_, ug_, wg_, ijk)

          ELSE IF (job_type == '3D') THEN

            IF( muscl == 0 ) THEN
              CALL masf( rgfe(ijk),  rgfn(ijk),  rgft(ijk),    &
                      rgfe(imjk), rgfn(ijmk), rgft(ijkm),   &
                      densg_, ug_, vg_, wg_ )
            ELSE
              CALL fmas(rgfe(ijk),  rgfn(ijk),  rgft(ijk),    &
                      rgfe(imjk), rgfn(ijmk), rgft(ijkm),   &
                      densg_, ug_, vg_, wg_, i, j, k, ijk)
            ENDIF

          ENDIF
          ! --- by hand inlining -- CALL calc_gas_mass_flux

          ! --- by hand inlining -- CALL calc_res( i, j, k, ijk, dg)
          resx = (rgfe(ijk) - rgfe(imjk)) * indx(i) * inx(i)
          resz = (rgft(ijk) - rgft(ijkm)) * indz(k)
          resy = 0.D0
          IF (job_type == '3D') resy = (rgfn(ijk) - rgfn(ijmk)) * indy(j)
!
          dg = rgp(ijk) - rgpn(ijk) + dt * (resx+resy+resz)
!            - dt * (r1(ijk)+r2(ijk)+r3(ijk)+r4(ijk)+r5(ijk))
          ! --- by hand inlining -- CALL calc_res( i, j, k, ijk, dg)
!
          IF ( ( DABS(dg) > conv_ .OR. DABS(dg) >= DABS(dgorig) ) ) THEN

            IF( nit == 1 .AND. loop == 1 ) dgorig = dg
            d3 = dg
            ! ... steepen the Newton's slope (accelerate)
            IF( kros < 2 .AND. loop == inmax ) abeta_ = 0.5D0 * inmax * abeta_

          ELSE IF ( DABS(dg) <= conv_ ) THEN

            EXIT

          END IF


        END DO 

        RETURN
      
!-----------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------

      SUBROUTINE calc_part_mass_flux
!
! ... compute the particle mass fluxes by using the stencil
! ... as defined in the 'subscr' module

      IMPLICIT NONE

      INTEGER :: is
!
        DO is = 1, nsolid

          IF( muscl == 0 ) THEN
            CALL first_nb(dens,rlk(:,is),ijk)
            CALL first_rnb(u,us(:,is),ijk)
            CALL first_rnb(w,ws(:,is),ijk)
          ELSE
            CALL nb_13(dens,rlk(:,is),ijk)
            CALL rnb_13(u,us(:,is),ijk)
            CALL rnb_13(w,ws(:,is),ijk)
          END IF

          IF (job_type == '2D') THEN

            CALL fmas(rsfe(ijk,is),  rsft(ijk,is),    &
                    rsfe(imjk,is), rsft(ijkm,is),   &
                    dens, u, w, ijk)

          ELSE IF (job_type == '3D') THEN

           IF( muscl == 0 ) THEN
             CALL first_rnb(v,vs(:,is),ijk)
             CALL masf(rsfe(ijk,is),  rsfn(ijk,is),  rsft(ijk,is),    &
                    rsfe(imjk,is), rsfn(ijmk,is), rsft(ijkm,is),   &
                    dens, u, v, w)
           ELSE
              CALL rnb_13(v,vs(:,is),ijk)
              CALL fmas(rsfe(ijk,is),  rsfn(ijk,is),  rsft(ijk,is),    &
                    rsfe(imjk,is), rsfn(ijmk,is), rsft(ijkm,is),   &
                    dens, u, v, w, i, j, k, ijk)
           END IF

          END IF
        END DO

      END SUBROUTINE calc_part_mass_flux

      END SUBROUTINE inner_loop


!----------------------------------------------------------------------
   END MODULE iterative_solver
!----------------------------------------------------------------------
