!----------------------------------------------------------------------
      MODULE iterative_solver
!----------------------------------------------------------------------
      USE grid, ONLY: dx, dy, dz, indx, indy, indz, inx
      IMPLICIT NONE
!
! ... convective mass fluxes
!
      REAL*8, DIMENSION(:),   ALLOCATABLE :: rgfe, rgfn, rgft
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: rsfe, rsfn, rsft
!
      REAL*8  :: omega, dg
      INTEGER :: inmax, maxout
!
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE iter
! ... This routine is the core of the iterative procedure to solve the
! ... mass-momentum and the interphase coupling
! ... (2D-3D-Compliant)
!
      USE convective_mass_fluxes, ONLY: fmas, masf
      USE dimensions
      USE domain_decomposition, ONLY: ncint, ncdom, data_exchange
      USE domain_decomposition, ONLY: myijk, meshinds
      USE enthalpy_matrix, ONLY: ftem
      USE eos_gas, ONLY: thermal_eosg, xgc, cg
      USE gas_solid_temperature, ONLY: tg
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      USE gas_solid_density, ONLY: rog, rgp, rgpn, rlk, rlkn
      USE grid, ONLY: fl_l
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE output_dump, ONLY: outp
      USE parallel, ONLY: mpime
      USE particles_constants, ONLY: rl, inrl
      USE phases_matrix, ONLY: assemble_matrix, assemble_all_matrix
      USE phases_matrix, ONLY: solve_velocities, solve_all_velocities
      USE phases_matrix, ONLY: mats_3phase
      USE pressure_epsilon, ONLY: p, ep
      USE set_indexes, ONLY: stencil, nb, rnb
      USE set_indexes, ONLY: subscr, imjk, ijmk, ijkm
      USE specific_heat_module, ONLY: cp
      USE tilde_energy, ONLY: htilde
      USE tilde_momentum, ONLY: tilde
      USE time_parameters, ONLY: time, dt, timestart, tpr
      USE control_flags, ONLY: job_type
!
      IMPLICIT NONE
!
      REAL*8, DIMENSION(:), ALLOCATABLE :: conv
      REAL*8, DIMENSION(:), ALLOCATABLE :: abeta
!
! ... HW performance monitor include file

! #include "f_hpm.h" 

!
      INTEGER :: i, j , k, ijk, is, imesh
      REAL*8 :: dgorig
      REAL*8 :: rls
      REAL*8 :: omega0
      REAL*8 :: d3, p3

      INTEGER :: nit, loop, mustit, kros
      LOGICAL, ALLOCATABLE :: converge(:)
      TYPE(stencil) :: u, v, w, dens
!
      ALLOCATE(conv(ncint), abeta(ncint))
      conv = 0.0D0
      abeta = 0.0D0
!
      ALLOCATE(converge(ncint))
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
      DO ijk = 1, ncint
        IF (fl_l(ijk) == 1) THEN
          CALL subscr(ijk)
          CALL assemble_matrix(ijk)
          CALL solve_velocities(ijk)
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
      DO ijk = 1, ncint
        IF (fl_l(ijk) == 1) THEN
          CALL subscr(ijk)

          CALL calc_gas_mass_flux
!
! ... compute the derivative of the gas mass residual
! ... with respect to gas pressure
!
          CALL betas(conv(ijk), abeta(ijk), ijk)
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
! ... Start the HW performance monitor
!      call f_hpmstart( 1, ' ITER ' )
!
      sor_loop: DO nit = 1, maxout
!
!         call f_hpmstart( 2, ' SWEEP ' )
!
! ... Compute fluxes and forces in the momentum equation
! ... (for fully implicit solution)
!         CALL tilde

         mesh_loop: DO ijk = 1, ncint

           converge(ijk) = .FALSE.

           IF( fl_l(ijk) /= 1 ) converge(ijk) = .TRUE.
           IF( fl_l(ijk) == 1 ) THEN

              CALL meshinds(ijk,imesh,i,j,k)
              CALL subscr(ijk)
!
              ! ... Compute the residual of the mass balance 
              ! ... equation of the gas phase

              CALL calc_res(dg)
              dgorig = dg
!
              IF(DABS(dg) <= conv(ijk)) THEN
              
                ! ... If the residual is lower then the prescribed limit
                ! ... compute the gas and particle densities and proceed 
                ! ... to next cell

                converge(ijk) = .TRUE.

                CALL calc_eps
                rgp(ijk) = ep(ijk) * rog(ijk)

              ELSE IF (DABS(dg) > conv(ijk)) THEN
!
                ! ... If the residual is higher then the prescribed limit
                ! ... start the inner (in-cell) iterative loop
                ! ... to correct pressure and velocities.
!
                d3 = dg
                p3 = p(ijk)
                kros = -1
!
                inner_loop: DO loop = 1, inmax
!
                  ! ... Correct the pressure at current cell
                  ! ... by using successively the Newton's method
                  ! ... (or the secant method where Newton is not
                  ! ... applicable) and the two-sided secant method.
                  ! ... The first time in a cell, skip pressure correction
                  ! ... to compute once the particle velocities 
                  ! ... and volumetric fractions ...
!
!                  CALL betas(conv(ijk), abeta(ijk), ijk)
                  IF( (loop > 1) .OR. (nit > 1) ) THEN
                    CALL padjust(p(ijk), kros, d3, p3, omega, abeta(ijk))
                  END IF
!
                  ! ... Use equation of state to calculate gas density 
                  ! ... from (new) pressure and (old) temperature.
!
                  CALL thermal_eosg(rog(ijk),tg(ijk),p(ijk),xgc(:,ijk))
!
                  rgp(ijk) = ep(ijk) * rog(ijk)
!
                  !IF( MOD( ijk, 100 ) == 0 ) call f_hpmstart( 5, ' ITER_mas ')
!
                  ! ... Update gas and particles velocities using the 
                  ! ... corrected pressure at current location. 
                  ! ... Pressure at neighbour cells could still be wrong.

                  CALL assemble_all_matrix(ijk)
                  CALL solve_all_velocities(ijk)

                  ! if( MOD(ijk,100) == 0 ) call f_hpmstop( 5 )
!
                  ! ... update particle and gas densities and the 
                  ! ... gas mass residual 'dg'

                  CALL update_eps
                  rgp(ijk) = ep(ijk) * rog(ijk)
                  CALL update_res(dg)
!
                  IF ((DABS(dg) > conv(ijk) .OR. DABS(dg) >= DABS(dgorig))) THEN
                    IF(nit == 1 .AND. loop == 1) dgorig=dg
                    d3=dg
                    ! ... steepen the Newton's slope (accelerate)
                    IF(kros < 2 .AND. loop == inmax) &
                       abeta(ijk) = 0.5D0 * inmax * abeta(ijk)
                  ELSE IF (DABS(dg) <= conv(ijk)) THEN
                    EXIT inner_loop
                  END IF
                END DO inner_loop
!
              END IF

           END IF    
         END DO mesh_loop
!
! ... Update gas and particles enthalpy
! ... (for fully implicit solution)
!
!         CALL htilde
!         CALL ftem

!        call f_hpmstop( 2 )

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
        IF( ALL(converge) ) mustit = 0
!
! ... mustit must be zero simultaneously in all subdomains
!
        CALL parallel_sum_integer(mustit, 1)
!
        IF(mustit == 0) THEN
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
       IF(MOD(nit,1000) == 0) omega=omega*0.9D0
!
      END DO sor_loop
!/////////////////////////////////////////////////////////////////////

! ... stop the HW performance monitor
!      call f_hpmstop( 1 )
!
!********************************************************************
! ... If the iterative sweep concluded without convergence
! ... report the number of cells where the procedure does not converge
!
      IF( mustit /= 0) THEN
!
        WRITE(6,700) (time+dt)
        WRITE(6,*) 'convergence on proc ',mpime,' : ', ALL(converge)
        IF (.NOT.ALL(converge)) WRITE(6,*) 'cells not converged (imesh,i,j,k): '
        DO ijk = 1, ncint
          IF (.NOT.converge(ijk)) THEN
            WRITE(6,*) imesh, i, j, k
          END IF
        END DO
 700    FORMAT('max number of iterations reached at time: ', F8.3)

        CALL error( ' iter ', 'max number of iters exceeded ', 1)
        omega=omega0

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
!-------------------------------------------------------------------------
      SUBROUTINE calc_part_mass_flux
! ... compute the particle mass fluxes by using the stencil
! ... as defined in the 'subscr' module

      IMPLICIT NONE
!
      DO is = 1, nsolid

        CALL nb(dens,rlk(:,is),ijk)
        CALL rnb(u,us(:,is),ijk)
        CALL rnb(w,ws(:,is),ijk)

        IF (job_type == '2D') THEN

          CALL fmas(rsfe(ijk,is),  rsft(ijk,is),    &
                    rsfe(imjk,is), rsft(ijkm,is),   &
                    dens, u, w, ijk)

        ELSE IF (job_type == '3D') THEN

          CALL rnb(v,vs(:,is),ijk)
          CALL fmas(rsfe(ijk,is),  rsfn(ijk,is),  rsft(ijk,is),    &
                    rsfe(imjk,is), rsfn(ijmk,is), rsft(ijkm,is),   &
                    dens, u, v, w, ijk)

        END IF
      END DO

      END SUBROUTINE calc_part_mass_flux
!----------------------------------------------------------------------
      SUBROUTINE calc_eps
! ... Use the Mass Balance equation of the solids
! ... to compute particle volumetric fractions and the void fraction.
!
      IMPLICIT NONE
      REAL*8 :: rlkx, rlky, rlkz, rls

      rlkx = 0.D0
      rlky = 0.D0
      rlkz = 0.D0
      rls = 0.D0

      DO is = 1, nsolid
        rlkx = (rsfe(ijk,is) - rsfe(imjk,is)) * indx(i) * inx(i)
        rlkz = (rsft(ijk,is) - rsft(ijkm,is)) * indz(k)
        IF (job_type == '3D') rlky = (rsfn(ijk,is) - rsfn(ijmk,is)) * indy(j)

        rlk( ijk, is) = rlkn( ijk,is ) - dt * (rlkx+rlky+rlkz)
!        - dt * (r1(ijk)+r2(ijk)+r3(ijk)+r4(ijk)+r5(ijk))
        IF( rlk(ijk,is) < 0.D0 ) rlk(ijk,is) = 0.D0

        rls = rls + rlk( ijk,is) * inrl(is)

      END DO
      
      IF(rls > 1.D0) THEN
        WRITE(8,*) 'warning1: mass is not conserved'
        WRITE(8,*) 'time,i,j,k,rls',time,i,j,k,rls
        CALL error('iter', 'warning: mass is not conserved',1)
      ENDIF

      ep(ijk) = 1.D0 - rls

      END SUBROUTINE calc_eps
!----------------------------------------------------------------------
      SUBROUTINE update_eps
! ... update the particle volumetric fractions and the void fraction

      IMPLICIT NONE

      CALL calc_part_mass_flux
      CALL calc_eps

      END SUBROUTINE update_eps
!----------------------------------------------------------------------
      SUBROUTINE calc_gas_mass_flux
! ... compute the gas mass fluxes by using the stencil, as defined in the
! ... 'subscr' module
!
      IMPLICIT NONE

        CALL nb(dens,rgp,ijk)
        CALL rnb(u,ug,ijk)
        CALL rnb(w,wg,ijk)

        IF (job_type == '2D') THEN
        
          CALL fmas(rgfe(ijk),  rgft(ijk),    &
                    rgfe(imjk), rgft(ijkm),   &
                    dens, u, w, ijk)

        ELSE IF (job_type == '3D') THEN

          CALL rnb(v,vg,ijk)
          CALL fmas(rgfe(ijk),  rgfn(ijk),  rgft(ijk),    &
                    rgfe(imjk), rgfn(ijmk), rgft(ijkm),   &
                    dens, u, v, w, ijk)
        ENDIF

      END SUBROUTINE calc_gas_mass_flux
!----------------------------------------------------------------------
      SUBROUTINE calc_res(res)
! ... Compute the residual of the Mass Balance equation of the gas phase
!
      IMPLICIT NONE
      REAL*8 :: resx, resy, resz
      REAL*8, INTENT(OUT) :: res
!      
      resx = 0.D0
      resy = 0.D0
      resz = 0.D0
!
      resx = (rgfe(ijk) - rgfe(imjk)) * indx(i) * inx(i)
      resz = (rgft(ijk) - rgft(ijkm)) * indz(k)
      IF (job_type == '3D') resy = (rgfn(ijk) - rgfn(ijmk)) * indy(j)
!
      res  = rgp(ijk) - rgpn(ijk) + dt * (resx+resy+resz)
!          - dt * (r1(ijk)+r2(ijk)+r3(ijk)+r4(ijk)+r5(ijk))

      RETURN
      END SUBROUTINE calc_res
!----------------------------------------------------------------------
      SUBROUTINE update_res(residual)
! ... update residual of the Mass Balance equation of the gas phase
! ... using guessed velocities and pressure.
!
      IMPLICIT NONE
      REAL*8 :: residual
!      
      CALL calc_gas_mass_flux
      CALL calc_res(residual)
!      
      RETURN
      END SUBROUTINE update_res
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
      SUBROUTINE padjust(p, kros, d3, p3, omega, abeta)
! ... Correct the pressure field to mimize the gas mass residual

        REAL*8 :: p, d3, p3, omega, abeta
        INTEGER :: kros
        REAL*8, SAVE :: dp, d1, d2, p1, p2

        IF(kros.NE.3) THEN

          IF(d3 > 0.D0) THEN
            d1=d3
            p1=p3
            IF(kros == -1) kros=0
            IF(kros == 1) kros=2
          ELSE
            d2=d3
            p2=p3
            IF(kros == -1) kros=1
            IF(kros == 0) kros=2
          END IF

          IF(kros == 2) THEN
! ... Use secant method once, when the sign of dg changes
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
            IF(-dp*dsign(1.D0,d3) > 2.5D-1*p3) THEN
              dp=-2.5D-1*dsign(1.D0,d3)*p3
            END IF 
            p = p + dp

          END IF

        ELSE
! ... Use two-sided secant method
          p = newp(d1, d2, d3, p1, p2, p3)
          IF(d3 > 0.D0) THEN
            d1=d3
            p1=p3
          ELSE
            d2=d3
            p2=p3
          END IF

        END IF
        p3=p
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
      SUBROUTINE betas(cnv, abt, ijk)
! 
      USE dimensions
      USE domain_decomposition, ONLY: ncint, myijk, meshinds
      USE convective_mass_fluxes, ONLY: upc_e, upc_n, upc_t
      USE convective_mass_fluxes, ONLY: upc_w, upc_s, upc_b
      USE control_flags, ONLY: job_type
      USE eos_gas, ONLY: csound
      USE gas_solid_density, ONLY: rog, rgp
      USE grid, ONLY: dx, dy, dz, fl_l
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
      INTEGER :: i, j, k, imesh
      INTEGER, INTENT(IN) :: ijk
      REAL*8, INTENT(OUT) :: cnv, abt
!
      REAL*8, PARAMETER :: delg=1.D-8
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

          iepx = (  iep_e + iep_w ) * indx(i) * inx(i) 
          iepz = (  iep_t + iep_b ) * indz(i)
	  IF (job_type == '2D') THEN
            iepy = 0.D0
	  ELSE IF (job_type == '3D') THEN
            iepy = (  iep_n + iep_s ) * indy(i)
          END IF

!
! ... Inverse of the squared sound velocity (perfect gas)
!
          CALL csound(sqc,rog(ijk),p(ijk))
          rags = 1.D0 / sqc
!
! ... rbeta = dD_g/dP

          rbeta = ep(ijk) * rags + dt**2.D0 * (iepx+iepy+iepz)
!
          abt = 1.D0 / rbeta
          cnv = delg * rgp(ijk)
!
      RETURN
      END SUBROUTINE betas
!----------------------------------------------------------------------
      END MODULE iterative_solver
!----------------------------------------------------------------------
