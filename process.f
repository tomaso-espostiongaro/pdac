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
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: array_map_max
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: topo2d
!
! ... derived fields
!
      REAL, ALLOCATABLE, DIMENSION(:) :: rm, rg, bd, m, um, vm, wm, mvm, c, mc 
      REAL, ALLOCATABLE, DIMENSION(:) :: epstot, lepstot, pd

      INTEGER :: imap
      REAL*8 :: deltaz
!
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE allocate_main_fields(dime)
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: dime

      ALLOCATE(p(dime), ug(dime), vg(dime), wg(dime), tg(dime))
      ALLOCATE(eps(dime,nsolid), us(dime,nsolid), vs(dime,nsolid), &
                ws(dime,nsolid), ts(dime,nsolid))
      ALLOCATE(xgc(dime,ngas))
      ALLOCATE(array_map(nx,ny))
      ALLOCATE(array_map_max(nx,ny))
      ALLOCATE(topo2d(nx,ny))
      array_map = 0.D0
      array_map_max = 0.D0
      topo2d = -9999

      RETURN
      END SUBROUTINE allocate_main_fields
!----------------------------------------------------------------------
      SUBROUTINE allocate_derived_fields(dime)
!
! ... These are the "more interesting" fields that can be derived
! ... from the primary OUTPUT fields produced by PDAC.
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: dime

      ALLOCATE(rm(dime))  ! Mixture Density
      ALLOCATE(rg(dime))  ! Gas Density
      ALLOCATE(bd(dime))  ! Bulk Density
      ALLOCATE(m(dime))   ! Gas Component Mass Fraction
      ALLOCATE(um(dime))  ! Mixture Velocity X
      ALLOCATE(vm(dime))  ! Mixture Velocity Y
      ALLOCATE(wm(dime))  ! Mixture Velocity Z
      ALLOCATE(mvm(dime)) ! Mixture Velocity Modulus
      ALLOCATE(c(dime))  ! Inverse of the Sound Speed
      ALLOCATE(mc(dime))  ! Mach Number
      ALLOCATE(epstot(dime))  ! Total particle fraction
      ALLOCATE(lepstot(dime))  ! Log10 of the total part. frac.
      ALLOCATE(pd(dime))  ! Dynamic Pressure

      RETURN
      END SUBROUTINE allocate_derived_fields
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
      INTEGER :: tn, ijk, i, k, nfil
      LOGICAL :: lform, ex
      CHARACTER(LEN = 14) :: filnam
      CHARACTER(LEN = 4 ) :: lettera
!
      lform = formatted_output

      CALL allocate_main_fields(ntot)
      CALL allocate_derived_fields(ntot)
!
      filnam='map_max.'//lettera(first_out-incr_out)
      INQUIRE(FILE=filnam,EXIST=ex)
      IF (ex) THEN
              OPEN(UNIT=tempunit,FILE=filnam)
              READ(tempunit,*) array_map_max(:,:)
              CLOSE(tempunit)
      END IF
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
!        rg = rhog(p,tg,xgc)
!        bd = rgp(eps,p,tg,xgc)
!        m  = mg(xgc)
        um = velm(ug,us,eps,p,tg,xgc)
        IF (job_type == '3D') vm = velm(vg,vs,eps,p,tg,xgc)
        wm = velm(wg,ws,eps,p,tg,xgc)
        mvm = vel(um,vm)
!        c  = cm(bd,rg,rm,m,tg)
!        mc = mach(mvm,c)
        pd = pdyn(rm,mvm)
        epstot = epst(eps)
        lepstot = leps(epstot)

        ! ... Write out fields of interest
        !
        filnam = 'log10epst.'//lettera(tn)
        IF (lform) THEN
          OPEN(tempunit,FILE=filnam)
        ELSE 
          OPEN(tempunit,FILE=filnam, FORM='UNFORMATTED')
        END IF
        CALL write_array( tempunit, lepstot, lform )
        CLOSE(tempunit)
        !
        filnam = 'tg.'//lettera(tn)
        IF (lform) THEN
          OPEN(tempunit,FILE=filnam)
        ELSE 
          OPEN(tempunit,FILE=filnam, FORM='UNFORMATTED')
        END IF
        CALL write_array( tempunit, tg, lform )
        CLOSE(tempunit)

        ! ... Print the map of any interesting variable above ground
        !
        IF (imap > 0) THEN
                CALL write_map_max(tn,pd)
                CALL write_map(tn,tg)
        END IF

      END DO
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
                ! ... Map the value reached at any given position at
                ! ... given time
                array_map(i,j) = map
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
      SUBROUTINE write_map_max(nfil,array)

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
      CHARACTER( LEN = 12 ) :: filnam
      LOGICAL :: ex

      IF (job_type == '2D') RETURN

      filnam='map_max.'//lettera(nfil)

      DO i = 1, nx
        DO j = 1, ny
          !
          search: DO k = 1, nz
            quota = improfile(i,j,k)
            IF (quota >= deltaz) THEN
                ijk  = i + (j-1) * nx + (k-1) * nx * ny
                ijkm = i + (j-1) * nx + (k-2) * nx * ny
                alpha = deltaz - z(k-1)
                alpha = alpha / (z(k) - z(k-1)) 
                map = alpha * array(ijk) + (1.D0-alpha) * array(ijk)
                ! ... Map the maximum value reached at any given position in
                ! ... time
                array_map_max(i,j) = MAX(map, array_map_max(i,j))
                EXIT search
            END IF
          END DO search
          !
        END DO
      END DO
!
! ... Print out the map of maxima
!
      OPEN(UNIT=tempunit,FILE=filnam)
      DO j = 1, ny
          WRITE(tempunit,122) (array_map_max(i,j), i=1, nx)
      END DO
      CLOSE(tempunit)

 122  FORMAT(10(1x,G14.6E3))
 123  CONTINUE

      RETURN
      END SUBROUTINE write_map_max
!----------------------------------------------------------------------
      SUBROUTINE sample
      USE control_flags, ONLY: job_type
      USE derived_fields
      USE dimensions, ONLY: nsolid, ngas, nx, ny, nz
      USE filter_outp, ONLY: read_points, probe_file
      USE filter_outp, ONLY: probe_point, number_of_probes, assign_index
      USE grid, ONLY: x, y, z
      USE io_files, ONLY: outpunit, logunit, tempunit
!
      IMPLICIT NONE
      INTEGER :: ijk, i, j, k, ig, is, n
      CHARACTER(LEN = 11) :: filnam
      CHARACTER(LEN = 4 ) :: lettera
      CHARACTER(LEN = 20) :: probenam
      REAL*4, ALLOCATABLE :: vars(:,:)
      INTEGER, ALLOCATABLE :: ijk_probe(:), indx(:)
      INTEGER :: nop, tn, nfil, nvars, nv, cnt
      REAL*4   :: time
      LOGICAL :: lform
      !
      TYPE(probe_point), ALLOCATABLE :: probe(:)
!
! ... Allocate probes.
! ... Read the file containing the indexes or coordinates of probe points
!
      ALLOCATE(probe(number_of_probes))
      ALLOCATE(ijk_probe(number_of_probes))
      ALLOCATE(indx(number_of_probes))
!
      CALL allocate_main_fields(number_of_probes)
      CALL allocate_derived_fields(number_of_probes)
!
      OPEN(tempunit, FILE=probe_file, STATUS='OLD')
      DO nop = 1, number_of_probes
        probe(nop)%nop = nop
        IF (assign_index) THEN
                IF (job_type == '3D') THEN
                        READ(tempunit,*) probe(nop)%i, probe(nop)%j, probe(nop)%k
                        i = probe(nop)%i
                        j = probe(nop)%j
                        k = probe(nop)%k
                        ijk = i + (j-1)*nx + (k-1)*nx*ny
                        ijk_probe(nop) = ijk
                        !
                        probe(nop)%x = x(i)
                        probe(nop)%y = y(j)
                        probe(nop)%z = z(k)
                ELSE IF (job_type == '2D') THEN
                        READ(tempunit,*) probe(nop)%i,probe(nop)%k
                        i = probe(nop)%i
                        k = probe(nop)%k
                        ijk = i + (k-1)*nx
                        ijk_probe(nop) = ijk
                        !
                        probe(nop)%x = x(i)
                        probe(nop)%z = z(k)
                END IF
        ELSE
                IF (job_type == '3D') THEN
                        READ(tempunit,*) probe(nop)%x, probe(nop)%y, probe(nop)%z
                        probe(nop)%i = 1
                        probe(nop)%j = 1
                        probe(nop)%k = 1
                ELSE IF (job_type == '2D') THEN
                        READ(tempunit,*) probe(nop)%x, probe(nop)%z
                        probe(nop)%i = 1
                        probe(nop)%k = 1
                END IF
        END IF
      END DO
      CLOSE(tempunit)
!
! ... Sort the probes with progressively increasing index
! ... After the sorting, 'ijk_probe' contains the progressively
! ... increasing probe indexes, whereas 'indx' contains the
! ... probe indexes as read from the 'probe_file'
!
      CALL ikb07ad(ijk_probe(1), number_of_probes, indx(1))
!
      lform = formatted_output
!
! ... Define the total number of basic PDAC output variables
! ... and allocate arrays
!
      IF (job_type == '2D') THEN
              nvars = 4 * (nsolid + 1) + ngas
      ELSE IF (job_type == '3D') THEN
              nvars = 5 * (nsolid + 1) + ngas
      END IF
      !
      ALLOCATE(vars(number_of_probes,nvars))
      vars = 0.0
!
! ... Loop over time-steps
!
      DO tn = first_out, last_out, incr_out

        filnam='output.'//lettera(tn)
        WRITE(logunit,fmt="(/,'* Starting sampling ',I5,' * ')" ) tn

        ! ... Open PDAC output file
        !
        IF (lform) THEN
          OPEN(UNIT=outpunit, FILE=filnam)
        ELSE 
          OPEN(UNIT=outpunit,FORM='UNFORMATTED',FILE=filnam)
        END IF
        !
        ! ... Read sampling points in the progressive order
        !
        CALL read_points(outpunit,lform,ijk_probe,time,vars)
        !
        ! ... Compute Derived Quantities
        !
        cnt = 0
        IF (job_type == '2D') THEN
          p = vars(:,1)
          ug = vars(:,2)
          wg = vars(:,3)
          tg = vars(:,4)
          DO ig = 1, ngas
            xgc(:,ig) = vars(:,4+ig)
          END DO
          cnt = 4 + ngas
          DO is = 1, nsolid
            eps(:,is) = vars(:,cnt+1)
            us(:,is) = vars(:,cnt+2)
            ws(:,is) = vars(:,cnt+3)
            ts(:,is) = vars(:,cnt+4)
            cnt = cnt + 4
          END DO
        ELSE IF (job_type == '3D') THEN
          p = vars(:,1)
          ug = vars(:,2)
          vg = vars(:,3)
          wg = vars(:,4)
          tg = vars(:,5)
          DO ig = 1, ngas
            xgc(:,ig) = vars(:,5+ig)
          END DO
          cnt = 5 + ngas
          DO is = 1, nsolid
            cnt = cnt + (is-1)*5
            eps(:,is) = vars(:,cnt+1)
            us(:,is) = vars(:,cnt+2)
            vs(:,is) = vars(:,cnt+3)
            ws(:,is) = vars(:,cnt+4)
            ts(:,is) = vars(:,cnt+5)
          END DO
        END IF
        DO is = 1, nsolid
        WRITE(*,*) eps(:,is)
        END DO
        rm = rhom(eps,p,tg,xgc)
        rg = rhog(p,tg,xgc)
        bd = rgp(eps,p,tg,xgc)
        m  = mg(xgc)
        um = velm(ug,us,eps,p,tg,xgc)
        IF (job_type == '3D') vm = velm(vg,vs,eps,p,tg,xgc)
        wm = velm(wg,ws,eps,p,tg,xgc)
        mvm = vel(um,wm)
        c  = cm(bd,rg,rm,m,tg)
        mc = mach(mvm,c)
        pd = pdyn(rm,mvm)
        epstot = epst(eps)
        lepstot = leps(epstot)
        !
        ! ... Loop over samplig points and write the corresponding files
        !
        DO nop = 1, number_of_probes
          n = indx(nop)
          i = probe(n)%i
          j = probe(n)%j
          k = probe(n)%k
          probenam ='S'//lettera(n)//'_'//lettera(i)//'_'//lettera(j)//'_'//lettera(k)
          OPEN(UNIT=tempunit, FILE=probenam, POSITION='APPEND')
            WRITE(tempunit,100) time, (vars(nop,nv), nv=1, nvars), rm(nop), pd(nop)
          CLOSE(tempunit)
        END DO

        CLOSE(outpunit)
      END DO
 100  FORMAT( F8.2, 100(G14.6E3,1X) )
!
      DEALLOCATE(vars)
      DEALLOCATE(probe)
      DEALLOCATE(ijk_probe)
      DEALLOCATE(indx)
!
      RETURN
      END SUBROUTINE sample
!----------------------------------------------------------------------
      END MODULE process_outp
!----------------------------------------------------------------------
