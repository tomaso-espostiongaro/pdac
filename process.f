!----------------------------------------------------------------------
      MODULE process_outp
!----------------------------------------------------------------------
      USE dimensions
      USE filter_outp, ONLY: first_out, last_out, incr_out
      USE filter_outp, ONLY: read_array, write_array
      USE kinds
      USE control_flags, ONLY: job_type
      USE io_files, ONLY: logunit
      USE output_dump, ONLY: formatted_output
!
      IMPLICIT NONE
      SAVE
!
! ... main fields
!
      REAL, ALLOCATABLE, DIMENSION(:)   :: p, ug, vg, wg, tg
      REAL, ALLOCATABLE, DIMENSION(:,:) :: xgc
      REAL, ALLOCATABLE, DIMENSION(:,:) :: eps, us, vs, ws, ts
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: array_map
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: topo2d

      INTEGER :: imap
      REAL*8 :: deltaz
!
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE allocate_main_fields
!
! ... These are the "more interesting" fields that can be derived
! ... from the primary OUTPUT fields produced by PDAC.
!
      IMPLICIT NONE

      ALLOCATE(p(ntot), ug(ntot), vg(ntot), wg(ntot), tg(ntot))
      ALLOCATE(eps(ntot,nsolid), us(ntot,nsolid), vs(ntot,nsolid), &
                ws(ntot,nsolid), ts(ntot,nsolid))
      ALLOCATE(xgc(ntot,ngas))
      ALLOCATE(array_map(nx,ny))
      ALLOCATE(topo2d(nx,ny))
      array_map = 0.D0
      topo2d = -9999

      RETURN
      END SUBROUTINE allocate_main_fields
!----------------------------------------------------------------------
      SUBROUTINE read_output( tn )
!
! ... Read THe PDAC Output files
!
      USE io_files, ONLY: outpunit
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: tn
      CHARACTER(LEN = 11) :: filnam
      CHARACTER(LEN = 4 ) :: lettera
      CHARACTER(LEN = 2 ) :: lettera2
      LOGICAL :: lform
!
      INTEGER :: ig, is, i
      REAL*8   :: time
      REAL*4   :: stime
!
      filnam='output.'//lettera(tn)

      lform = formatted_output

      IF (lform) THEN
        OPEN(UNIT=outpunit, FILE=filnam, STATUS='OLD')
        READ(outpunit,'(1x,///,1x,"@@@ TIME = ",g11.4)') time
      ELSE 
        OPEN(UNIT=outpunit,FORM='UNFORMATTED',FILE=filnam)
        READ(outpunit) stime
      END IF

      IF( lform ) READ(outpunit,'(///)')
      CALL read_array( outpunit, p, lform )  ! gas_pressure

      IF (job_type == '2D') THEN

        IF( lform ) READ(outpunit,'(///)')
        CALL read_array( outpunit, ug, lform ) ! gas_velocity_r
        IF( lform ) READ(outpunit,'(///)')
        CALL read_array( outpunit, wg, lform ) ! gas_velocity_z

      ELSE IF (job_type == '3D') THEN

        IF( lform ) READ(outpunit,'(///)')
        CALL read_array( outpunit, ug, lform ) ! gas_velocity_x
        IF( lform ) READ(outpunit,'(///)')
        CALL read_array( outpunit, vg, lform ) ! gas_velocity_y
        IF( lform ) READ(outpunit,'(///)')
        CALL read_array( outpunit, wg, lform ) ! gas_velocity_z

      ELSE
        CALL error('outp_','Unknown job type',1)
      END IF

      IF( lform ) READ(outpunit,'(///)')
      CALL read_array( outpunit, tg, lform )  ! gas_temperature

      DO ig=1,ngas
        IF( lform ) READ(outpunit,'(///)')
        CALL read_array( outpunit, xgc(:,ig), lform )  ! gc_molar_fraction
      END DO

      DO is = 1, nsolid

        IF( lform ) READ(outpunit,'(///)')
        CALL read_array( outpunit, eps(:,is), lform )  ! solid_bulk_density

        IF (job_type == '2D') THEN

        IF( lform ) READ(outpunit,'(///)')
          CALL read_array( outpunit, us(:,is), lform )  ! solid_velocity_r
        IF( lform ) READ(outpunit,'(///)')
          CALL read_array( outpunit, ws(:,is), lform )  ! solid_velocity_z

        ELSE IF (job_type == '3D') THEN

        IF( lform ) READ(outpunit,'(///)')
          CALL read_array( outpunit, us(:,is), lform )  ! solid_velocity_x
        IF( lform ) READ(outpunit,'(///)')
          CALL read_array( outpunit, vs(:,is), lform )  ! solid_velocity_y
        IF( lform ) READ(outpunit,'(///)')
          CALL read_array( outpunit, ws(:,is), lform )  ! solid_velocity_z

        END IF

        IF( lform ) READ(outpunit,'(///)')
        CALL read_array( outpunit, ts(:,is), lform )  ! solid_temperature

      END DO

      CLOSE (outpunit)
!
      RETURN
      END SUBROUTINE read_output
!-----------------------------------------------------------------------
      SUBROUTINE process
! 
! ... Compute the derived fields
!
      USE derived_fields
      USE io_files, ONLY: tempunit
      USE grid, ONLY: z

      IMPLICIT NONE
!
      INTEGER :: tn, ijk, i, k
      LOGICAL :: lform
      REAL, ALLOCATABLE, DIMENSION(:) :: rm, rg, bd, m, um, vm, wm, mvm, c, mc 
      REAL, ALLOCATABLE, DIMENSION(:) :: epstot, lepstot, pd
      CHARACTER(LEN = 14) :: filnam
      CHARACTER(LEN = 4 ) :: lettera
!
      lform = formatted_output

      ALLOCATE(rm(ntot))  ! Mixture Density
      ALLOCATE(rg(ntot))  ! Gas Density
      ALLOCATE(bd(ntot))  ! Bulk Density
      ALLOCATE(m(ntot))   ! Gas Component Mass Fraction
      ALLOCATE(um(ntot))  ! Mixture Velocity X
      ALLOCATE(vm(ntot))  ! Mixture Velocity Y
      ALLOCATE(wm(ntot))  ! Mixture Velocity Z
      ALLOCATE(mvm(ntot)) ! Mixture Velocity Modulus
      ALLOCATE(c(ntot))  ! Inverse of the Sound Speed
      ALLOCATE(mc(ntot))  ! Mach Number
      ALLOCATE(epstot(ntot))  ! Total particle fraction
      ALLOCATE(lepstot(ntot))  ! Log10 of the total part. frac.
      ALLOCATE(pd(ntot))  ! Dynamic Pressure

      CALL allocate_main_fields
!
      DO tn = first_out, last_out, incr_out

        WRITE(logunit,fmt="(/,'* Starting post-processing ',I5,' * ')" ) tn

        ! ... Read PDAC output file
        !
        CALL read_output ( tn )

        ! ... Derived fields are computed as a function of
        ! ... primary fields and other derived fields
        !
        rm = rhom(eps,p,tg,xgc)
        rg = rhog(p,tg,xgc)
        bd = rgp(eps,p,tg,xgc)
        m  = mg(xgc)
        um = velm(ug,us,eps,p,tg,xgc)
        IF (job_type == '3D') vm = velm(vg,vs,eps,p,tg,xgc)
        wm = velm(wg,ws,eps,p,tg,xgc)
        mvm = vel(um,vm)
        c  = cm(bd,rg,rm,m,tg)
        mc = mach(mvm,c)
        pd = pdyn(rm,mvm)
        epstot = epst(eps)
        lepstot = leps(epstot)

        ! ... Write out fields of interest
        !
        filnam = 'log10epst.'//lettera(tn)
        OPEN(tempunit,FILE=filnam, STATUS='NEW', FORM='UNFORMATTED')
        CALL write_array( tempunit, lepstot, lform )
        CLOSE(tempunit)
        !
        filnam = 'tg.'//lettera(tn)
        OPEN(tempunit,FILE=filnam, STATUS='NEW', FORM='UNFORMATTED')
        CALL write_array( tempunit, tg, lform )
        CLOSE(tempunit)

        ! ... Print the map of any interesting variable above ground
        !
        IF (imap > 0) CALL write_map(tn,rm)

      END DO
!
      DEALLOCATE(rm)
      DEALLOCATE(rg)
      DEALLOCATE(bd)
      DEALLOCATE(m)
      DEALLOCATE(um)
      DEALLOCATE(vm)
      DEALLOCATE(wm)
      DEALLOCATE(mvm)
      DEALLOCATE(c)
      DEALLOCATE(mc)
      DEALLOCATE(epstot)
      DEALLOCATE(lepstot)
      DEALLOCATE(pd)
! 
      RETURN
      END SUBROUTINE process
!-----------------------------------------------------------------------
      SUBROUTINE write_map(nfil,array)

      USE control_flags, ONLY: job_type
      USE dimensions, ONLY: nx, ny, nz
      USE filter_outp, ONLY: improfile
      USE grid, ONLY: z
      USE io_files, ONLY: tempunit
      IMPLICIT NONE
      
      REAL, INTENT(IN), DIMENSION(:) :: array
      INTEGER, INTENT(IN) :: nfil
      REAL*8 :: alpha, map, quota
      INTEGER :: i, j, k, ijk, ijkm
      CHARACTER( LEN = 4 ) :: lettera
      CHARACTER( LEN = 8 ) :: filnam

      IF (job_type == '2D') RETURN

      filnam='map.'//lettera(nfil)

      DO i = 1, nx
        DO j = 1, ny
          !
          search: DO k = 1, nz
            quota = improfile(i,j,k)
            IF (quota >= 0.D0 .AND. topo2d(i,j) == -9999) topo2d(i,j) = z(k) - quota
            IF (quota >= deltaz) THEN
                ijk  = i + (j-1) * nx + (k-1) * nx * ny
                ijkm = i + (j-1) * nx + (k-2) * nx * ny
                alpha = deltaz - z(k-1)
                alpha = alpha / (z(k) - z(k-1)) 
                map = alpha * array(ijk) + (1.D0-alpha) * array(ijk)
                ! ... Map the maximum value reached at any given position in
                ! ... time
                array_map(i,j) = MAX(map, array_map(i,j))
                !array_map(i,j) = map
                EXIT search
            END IF
          END DO search
          !
        END DO
      END DO
!
! ... Print out the map and the new 2D DEM file
!
      OPEN(UNIT=tempunit,FILE=filnam)
      DO j = 1, ny
          WRITE(tempunit,122) (array_map(i,j), i=1, nx)
      END DO
      CLOSE(tempunit)
      OPEN(UNIT=tempunit,FILE='topo2d.dat')
      DO j = 1, ny
          WRITE(tempunit,122) (topo2d(i,j), i=1, nx)
      END DO
      CLOSE(tempunit)

 122  FORMAT(10(1x,G14.6E3))
 123  CONTINUE

      RETURN
      END SUBROUTINE write_map
!-----------------------------------------------------------------------
      END MODULE process_outp
!----------------------------------------------------------------------
