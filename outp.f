!----------------------------------------------------------------------
      MODULE output_dump
!----------------------------------------------------------------------
      USE control_flags, ONLY: nfil, formatted_output
      USE kinds
      USE io_files, ONLY: logunit

      IMPLICIT NONE
      PRIVATE

      LOGICAL :: interpolate = .TRUE.
      PUBLIC :: outp, write_radial_profile_2D, shock_tube_out
      PUBLIC :: print_volumes, cell_report, write_mean_fields
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE jet_out
!
      USE dimensions, ONLY: nsolid
      USE gas_constants, ONLY: gas_type
      USE gas_solid_density, ONLY: rlk, rog
      USE gas_solid_velocity, ONLY: ug, vg, wg
      USE gas_solid_velocity, ONLY: us, vs, ws
      USE gas_solid_temperature, ONLY: tg, ts
      USE io_parallel, ONLY: write_array
      USE parallel, ONLY: nproc, mpime, root, group
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: p
      USE time_parameters, ONLY: time
      USE turbulence_model, ONLY: modturbo
      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE io_files, ONLY: filnam, outpunit
!
      IMPLICIT NONE
!
      CHARACTER*4 :: lettera
      LOGICAL :: lform
!
      INTEGER :: ig,is
      REAL*8, ALLOCATABLE :: otmp(:)
!
    
      nfil=nfil+1
      filnam='output.'//lettera(nfil)
      lform = formatted_output

      IF( mpime == root ) THEN

        IF (lform) THEN
          OPEN(UNIT=outpunit,FILE=filnam)
          WRITE(outpunit,*) time
        ELSE 
          OPEN(UNIT=outpunit,FORM='UNFORMATTED',FILE=filnam)
          WRITE(outpunit) REAL(time,4)
        END IF
 
      END IF
!
      CALL write_array( outpunit, p, sgl, lform )  ! gas_pressure
      CALL write_array( outpunit, tg, sgl, lform )   ! gas_temperature
      CALL write_array( outpunit, rog, sgl, lform )  ! gas_density
!
      IF (job_type == JOB_TYPE_2D) THEN
        CALL write_array( outpunit, ug, sgl, lform ) ! gas_velocity_r
        CALL write_array( outpunit, wg, sgl, lform ) ! gas_velocity_z
      ELSE IF (job_type == JOB_TYPE_3D) THEN
        CALL write_array( outpunit, ug, sgl, lform ) ! gas_velocity_x
        CALL write_array( outpunit, vg, sgl, lform ) ! gas_velocity_y
        CALL write_array( outpunit, wg, sgl, lform ) ! gas_velocity_z
      ELSE
        CALL error('outp_','Unknown job type',1)
      END IF

      IF( mpime == root ) THEN
        CLOSE (outpunit)
      END IF
!
      RETURN
      END SUBROUTINE jet_out
!----------------------------------------------------------------------
      SUBROUTINE shock_tube_out
!
      USE dimensions, ONLY: nx, nz, nsolid
      USE domain_mapping, ONLY: ncint, meshinds
      USE gas_constants, ONLY: gas_type, gammaair
      USE gas_solid_density, ONLY: rog, rlk
      USE gas_solid_velocity, ONLY: ug, wg, us, ws
      USE gas_solid_temperature, ONLY: sieg, tg, ts
      USE grid, ONLY: flag, z
      USE parallel, ONLY: nproc, mpime, root, group
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: p
      USE time_parameters, ONLY: time
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE domain_mapping, ONLY: data_collect, data_distribute
      USE io_files, ONLY: tempunit
!
      IMPLICIT NONE
!
      CHARACTER( LEN = 15 ) :: filnam
      CHARACTER( LEN = 4 ) :: lettera
!
      INTEGER :: i, j, k, ijk, ig, imesh
      REAL*8 :: energy
!
      nfil = nfil + 1
      filnam = 'shtube.' // lettera( nfil )

      IF( mpime == root ) THEN
        OPEN(UNIT=tempunit,FILE=filnam)
        DO ijk = 1, ncint
          IF ( BTEST(flag(ijk),0) ) THEN
            CALL meshinds(ijk, imesh, i, j, k)
            energy = p(ijk)/(rog(ijk)*(gammaair - 1.D0))
            WRITE(tempunit,550) z(k), rog(ijk),wg(ijk),tg(ijk),p(ijk),energy,tg(ijk),rlk(ijk,1),ws(ijk,1),ts(ijk,1)
          END IF
        END DO
        CLOSE (tempunit)
      END IF
!
      RETURN

 550  FORMAT(1x,20(1x,g12.6))
!
      END SUBROUTINE shock_tube_out
!----------------------------------------------------------------------
      SUBROUTINE write_radial_profile_2d
!
      USE dimensions, ONLY: nsolid, ngas, nx, nz
      USE eos_gas, ONLY: xgc, ygc
      USE gas_constants, ONLY: gas_type
      USE gas_solid_density, ONLY: rlk
      USE gas_solid_velocity, ONLY: ug, vg, wg
      USE gas_solid_velocity, ONLY: us, vs, ws
      USE gas_solid_temperature, ONLY: tg, ts
      USE grid, ONLY: zb, r
      USE io_parallel, ONLY: write_array
      USE parallel, ONLY: nproc, mpime, root, group
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: p
      USE time_parameters, ONLY: time
      USE turbulence_model, ONLY: modturbo
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE volcano_topography, ONLY: rim_quota
      USE io_files, ONLY: tempunit
      IMPLICIT NONE
      INTEGER :: n, is, ig, kq, k, nq

      DO k = 1, nz
        IF (zb(k) <= rim_quota) kq = k
      END DO
      
      OPEN(tempunit,FILE='radial_profile.dat')
      WRITE(tempunit,*) nx
      WRITE(tempunit,*) (r(n), n=1, nx)
      WRITE(tempunit,*) (ug(n + (kq-1)*nx), n=1, nx)
      WRITE(tempunit,*) (wg(n + (kq-1)*nx), n=1, nx)
      WRITE(tempunit,*) (tg(n + (kq-1)*nx), n=1, nx)
      WRITE(tempunit,*) (p(n + (kq-1)*nx), n=1, nx)
      DO ig = 1, ngas
        WRITE(tempunit,*) (ygc(n,ig), n=1, nx)
      END DO
      DO is = 1, nsolid
        WRITE(tempunit,*) (us(n + (kq-1)*nx,is), n=1, nx)
        WRITE(tempunit,*) (ws(n + (kq-1)*nx,is), n=1, nx)
        WRITE(tempunit,*) (ts(n + (kq-1)*nx,is), n=1, nx)
        WRITE(tempunit,*) (rlk(n + (kq-1)*nx,is)*inrl(is), n=1, nx)
      END DO
      CLOSE(tempunit)

      OPEN(tempunit,FILE='radial_profile_gas.dat')
      DO n = 2, nx-1
        nq = n + (kq-1)*nx
        WRITE(tempunit,188) r(n), ug(nq), wg(nq), tg(nq), p(nq), (ygc(nq,ig), ig=1, ngas)
      END DO
      CLOSE(tempunit)

      OPEN(tempunit,FILE='radial_profile_part.dat')
      DO n = 2, nx-1
        nq = n + (kq-1)*nx
        WRITE(tempunit,199) r(n), (us(nq,is), ws(nq,is), ts(nq,is), rlk(nq,is)*inrl(is), is=1,nsolid)
      END DO
      CLOSE(tempunit)

 188  FORMAT(7(F15.5))       
 199  FORMAT(9(F15.5))       
        
      RETURN
      END SUBROUTINE write_radial_profile_2d
!----------------------------------------------------------------------
      SUBROUTINE outp
!
      USE dimensions, ONLY: nsolid, ngas
      USE eos_gas, ONLY: xgc
      USE gas_constants, ONLY: gas_type
      USE gas_solid_density, ONLY: rlk
      USE gas_solid_velocity, ONLY: ug, vg, wg
      USE gas_solid_velocity, ONLY: us, vs, ws
      USE gas_solid_temperature, ONLY: tg, ts
      USE io_parallel, ONLY: write_array
      USE mixture_fields, ONLY: compute_mixture,rhom,um,vm,wm,tm,yd
      USE parallel, ONLY: nproc, mpime, root, group
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: p
      USE time_parameters, ONLY: time
      USE turbulence_model, ONLY: modturbo
      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE io_files, ONLY: filnam, outpunit
!
      IMPLICIT NONE
!
      CHARACTER( LEN =  4 ) :: lettera
      LOGICAL :: lform
!
      INTEGER :: ig,is
      REAL*8, ALLOCATABLE :: otmp(:)
!
    
      filnam='output.'//lettera(nfil)
      lform = formatted_output

      IF( mpime == root ) THEN

        IF (lform) THEN
          OPEN(UNIT=outpunit,FILE=filnam)
          WRITE(outpunit,'(1x,///,1x,"@@@ TIME = ",g11.4)') time
        ELSE 
          OPEN(UNIT=outpunit,FORM='UNFORMATTED',FILE=filnam)
          WRITE(outpunit) REAL(time,4)
        END IF

        WRITE(logunit,fmt="('  from outp: writing file ',A20)") filnam
 
      END IF
!
      IF( lform .AND. mpime == root ) &
                      WRITE(outpunit,'(1x,//,1x,"P   ",/)')
      CALL write_array( outpunit, p, sgl, lform )  ! gas_pressure

      IF (job_type == JOB_TYPE_2D) THEN

        IF( lform .AND. mpime == root ) &
                      WRITE(outpunit,'(1x,//,1x,"UG  ",/)')
        CALL write_array( outpunit, ug, sgl, lform ) ! gas_velocity_r

        IF( lform .AND. mpime == root ) &
                      WRITE(outpunit,'(1x,//,1x,"WG  ",/)')
        CALL write_array( outpunit, wg, sgl, lform ) ! gas_velocity_z

      ELSE IF (job_type == JOB_TYPE_3D) THEN

        IF( lform .AND. mpime == root ) &
                      WRITE(outpunit,'(1x,//,1x,"UG  ",/)')
        CALL write_array( outpunit, ug, sgl, lform ) ! gas_velocity_x

        IF( lform .AND. mpime == root ) &
                      WRITE(outpunit,'(1x,//,1x,"VG  ",/)')
        CALL write_array( outpunit, vg, sgl, lform ) ! gas_velocity_y

        IF( lform .AND. mpime == root ) &
                      WRITE(outpunit,'(1x,//,1x,"WG  ",/)')
        CALL write_array( outpunit, wg, sgl, lform ) ! gas_velocity_z

      ELSE
        CALL error('outp_','Unknown job type',1)
      END IF

      IF( lform .AND. mpime == root ) &
                      WRITE(outpunit,'(1x,//,1x,"TG  ",/)')
      CALL write_array( outpunit, tg, sgl, lform )  ! gas_temperature
!
      ALLOCATE( otmp( SIZE( xgc, 1 ) ) )
      DO ig=1,ngas
          otmp = xgc(:,ig)
          IF( lform .AND. mpime == root ) &
                      WRITE(outpunit,'(1x,//,1x,"XGC",I1,/)') ig
          CALL write_array( outpunit, otmp, sgl, lform )  ! gc_molar_fraction
      END DO
      DEALLOCATE( otmp )
!
      ALLOCATE( otmp( SIZE( rlk, 1 ) ) )

      DO is = 1, nsolid

        otmp = rlk(:,is)*inrl(is)
        IF( lform .AND. mpime == root ) &
                      WRITE(outpunit,'(1x,//,1x,"EPS",I1,/)') is
        CALL write_array( outpunit, otmp, sgl, lform )  ! solid_bulk_density

        IF (job_type == JOB_TYPE_2D) THEN

          IF( lform .AND. mpime == root ) &
                      WRITE(outpunit,'(1x,//,1x," US",I1,/)') is
          CALL write_array( outpunit, us(:,is), sgl, lform )  ! solid_velocity_r

          IF( lform .AND. mpime == root ) &
                      WRITE(outpunit,'(1x,//,1x," WS",I1,/)') is
          CALL write_array( outpunit, ws(:,is), sgl, lform )  ! solid_velocity_z

        ELSE IF (job_type == JOB_TYPE_3D) THEN

          IF( lform .AND. mpime == root ) &
                      WRITE(outpunit,'(1x,//,1x," US",I1,/)') is
          CALL write_array( outpunit, us(:,is), sgl, lform )  ! solid_velocity_x

          IF( lform .AND. mpime == root ) &
                      WRITE(outpunit,'(1x,//,1x," VS",I1,/)') is
          CALL write_array( outpunit, vs(:,is), sgl, lform )  ! solid_velocity_y

          IF( lform .AND. mpime == root ) &
                      WRITE(outpunit,'(1x,//,1x," WS",I1,/)') is
          CALL write_array( outpunit, ws(:,is), sgl, lform )  ! solid_velocity_z

        END IF

        IF( lform .AND. mpime == root ) &
                      WRITE(outpunit,'(1x,//,1x," TS",I1,/)') is
        CALL write_array( outpunit, ts(:,is), sgl, lform )  ! solid_temperature

      END DO

      DEALLOCATE( otmp )
!
! ... Write mixture fields if required
!
      IF( lform .AND. mpime == root ) &
                      WRITE(outpunit,'(1x,//,1x,"RHOM",/)')
      CALL write_array( outpunit, rhom, sgl, lform )  ! gas_pressure

      IF (job_type == JOB_TYPE_2D) THEN

        IF( lform .AND. mpime == root ) &
                      WRITE(outpunit,'(1x,//,1x,"UM  ",/)')
        CALL write_array( outpunit, um, sgl, lform ) ! gas_velocity_r

        IF( lform .AND. mpime == root ) &
                      WRITE(outpunit,'(1x,//,1x,"WM  ",/)')
        CALL write_array( outpunit, wm, sgl, lform ) ! gas_velocity_z

      ELSE IF (job_type == JOB_TYPE_3D) THEN

        IF( lform .AND. mpime == root ) &
                      WRITE(outpunit,'(1x,//,1x,"UM  ",/)')
        CALL write_array( outpunit, um, sgl, lform ) ! gas_velocity_x

        IF( lform .AND. mpime == root ) &
                      WRITE(outpunit,'(1x,//,1x,"VM  ",/)')
        CALL write_array( outpunit, vm, sgl, lform ) ! gas_velocity_y

        IF( lform .AND. mpime == root ) &
                      WRITE(outpunit,'(1x,//,1x,"WM  ",/)')
        CALL write_array( outpunit, wm, sgl, lform ) ! gas_velocity_z

      ELSE
        CALL error('outp_','Unknown job type',1)
      END IF

      IF( lform .AND. mpime == root ) &
                      WRITE(outpunit,'(1x,//,1x,"TM  ",/)')
      CALL write_array( outpunit, tm, sgl, lform )  ! gas_temperature
!
      ALLOCATE( otmp( SIZE( yd, 1 ) ) )
      DO ig=1,ngas+nsolid
          otmp = yd(:,ig)
          IF( lform .AND. mpime == root ) &
                      WRITE(outpunit,'(1x,//,1x,"YD",I1,/)') ig
          CALL write_array( outpunit, otmp, sgl, lform )  ! gc_molar_fraction
      END DO
      DEALLOCATE( otmp )
!
      IF( mpime == root ) THEN
        CLOSE (outpunit)
      END IF

      nfil=nfil+1
!
      RETURN
      END SUBROUTINE outp
!-----------------------------------------------------------------------
      SUBROUTINE write_mean_fields(filenumber)
!
! ... Write the mean field of all process fields
!
      USE compute_mean_fields, ONLY: rhom_gav, tm_gav, um_gav, vm_gav, wm_gav, yd_gav
      USE control_flags, ONLY: formatted_output,  job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE kinds
      USE dimensions, ONLY: nsolid, ngas
      USE io_files, ONLY: tempunit
      USE io_parallel, ONLY: write_array
      USE parallel, ONLY: mpime, root
      USE time_parameters, ONLY: time
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: filenumber

      CHARACTER(LEN = 14) :: filnam
      CHARACTER*4 :: lettera
      LOGICAL :: lform
!
      INTEGER :: ig,is
      REAL*8, ALLOCATABLE :: otmp(:)
!
      filnam='m_outp.'//lettera(filenumber)
      lform = formatted_output

      IF (mpime == root ) THEN
        IF (lform) THEN
          OPEN(UNIT=tempunit,FILE=filnam)
          WRITE(tempunit,*) time
        ELSE
          OPEN(UNIT=tempunit,FORM='UNFORMATTED',FILE=filnam)
          WRITE(tempunit) REAL(time,4)
        END IF
      END IF
!
      IF(lform .AND. mpime==root) WRITE(tempunit,'(1x,//,1x,"RHOM",/)')
      CALL write_array( tempunit, rhom_gav, sgl, lform )  !  mixture_density
!
      IF (job_type == JOB_TYPE_2D) THEN
      IF(lform .AND. mpime==root) WRITE(tempunit,'(1x,//,1x,"UM  ",/)')
        CALL write_array( tempunit, um_gav, sgl, lform)
      IF(lform .AND. mpime==root) WRITE(tempunit,'(1x,//,1x,"VM  ",/)')
        CALL write_array( tempunit, wm_gav, sgl, lform)
      ELSE IF (job_type == JOB_TYPE_3D) THEN
      IF(lform .AND. mpime==root) WRITE(tempunit,'(1x,//,1x,"UM  ",/)')
        CALL write_array( tempunit, um_gav, sgl, lform)
      IF(lform .AND. mpime==root) WRITE(tempunit,'(1x,//,1x,"VM  ",/)')
        CALL write_array( tempunit, vm_gav, sgl, lform)
      IF(lform .AND. mpime==root) WRITE(tempunit,'(1x,//,1x,"WM  ",/)')
        CALL write_array( tempunit, wm_gav, sgl, lform)
      ELSE
        CALL error('outp_','Unknown job type',1)
      END IF
      IF(lform .AND. mpime == root) WRITE(tempunit,'(1x,//,1x,"TM  ",/)')
      CALL write_array( tempunit, tm_gav, sgl, lform)
!
      ALLOCATE( otmp( SIZE( yd_gav, 1 ) ) )
      DO ig=1,ngas+nsolid
          otmp = yd_gav(:,ig)
          IF(lform .AND. mpime==root) &
                         WRITE(tempunit,'(1x,//,1x,"YD",I1,/)') ig
          CALL write_array( tempunit, otmp, sgl, lform )
      END DO
      DEALLOCATE( otmp )
!
      IF (mpime == root) CLOSE(tempunit)
!
      RETURN
      END SUBROUTINE write_mean_fields
!----------------------------------------------------------------------
      SUBROUTINE cell_report(iunit, ijk, imesh, i, j, k)
      USE dimensions, ONLY: nsolid, ngas
      USE eos_gas, ONLY: xgc
      USE gas_solid_density, ONLY: rlk
      USE gas_solid_velocity, ONLY: ug, vg, wg
      USE gas_solid_velocity, ONLY: us, vs, ws
      USE gas_solid_temperature, ONLY: tg, ts
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: p
      USE time_parameters, ONLY: time
      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iunit, ijk, imesh, i, j, k
      INTEGER :: ig, is

      IF (job_type == JOB_TYPE_3D) THEN
        WRITE(iunit,100) time, p(ijk), ug(ijk), vg(ijk), wg(ijk), tg(ijk),  &
        (xgc(ijk,ig), ig=1, ngas), (rlk(ijk,is)*inrl(is), us(ijk,is), &
        vs(ijk,is), ws(ijk,is), ts(ijk,is), is=1, nsolid)
      ELSE IF (job_type == JOB_TYPE_2D) THEN
        WRITE(iunit,100) time, p(ijk), ug(ijk), wg(ijk), tg(ijk),  &
        (xgc(ijk,ig), ig=1, ngas), (rlk(ijk,is)*inrl(is), us(ijk,is), &
        ws(ijk,is), ts(ijk,is), is=1, nsolid)
      END IF
 100  FORMAT( F8.2, 100(G14.6E3,1X) )

      RETURN
      END SUBROUTINE cell_report
!----------------------------------------------------------------------
      SUBROUTINE print_volumes
!
      USE immersed_boundaries, ONLY: vf
      USE io_files, ONLY: tempunit
      USE io_parallel, ONLY: write_array
      USE parallel, ONLY: mpime, root
      IMPLICIT NONE
      LOGICAL :: vform
!
      vform = .FALSE.
      IF (mpime == root) OPEN(tempunit,FILE='vf.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
!
      CALL write_array( tempunit, vf, sgl, vform )
!
      IF (mpime == root) CLOSE(tempunit)
!
      RETURN
      END SUBROUTINE print_volumes
!----------------------------------------------------------------------
      END MODULE output_dump
!----------------------------------------------------------------------
