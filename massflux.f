!----------------------------------------------------------------------
      MODULE mass_orthoflux
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
!
! ... main fields
!
      REAL, ALLOCATABLE, DIMENSION(:)   :: p, ug, vg, wg, tg, eptemp
      REAL, ALLOCATABLE, DIMENSION(:,:) :: xgc
      REAL, ALLOCATABLE, DIMENSION(:,:) :: eps, us, vs, ws, ts
!
      TYPE orthoslice 
        INTEGER :: axis
        INTEGER :: plane
        INTEGER :: n1
        INTEGER :: n2
        INTEGER :: t1
        INTEGER :: t2
      END TYPE orthoslice 
!
      INTEGER :: number_of_planes
      CHARACTER(LEN=80) :: planes_file
      REAL*8   :: time
      REAL*4   :: stime
!
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE allocate_main_fields(dime)
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: dime

      ALLOCATE(p(dime), ug(dime), vg(dime), wg(dime), tg(dime))
      ALLOCATE(eptemp(dime))
      ALLOCATE(eps(dime,nsolid), us(dime,nsolid), vs(dime,nsolid), &
                ws(dime,nsolid), ts(dime,nsolid))
      ALLOCATE(xgc(dime,ngas))

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
      INTEGER :: ig, is, i, k, ijk
!
      filnam='output.'//lettera(tn)
!
! ... OLD OUTPUT name
!
!      filnam='OUTPUT.'//lettera(tn)

      lform = formatted_output
      IF (lform) THEN
        OPEN(UNIT=outpunit, FILE=filnam, STATUS='OLD')
        READ(outpunit,'(1x,///,1x,"@@@ TIME = ",g11.4)') time
      ELSE 
        OPEN(UNIT=outpunit,FORM='UNFORMATTED',FILE=filnam)
        READ(outpunit) stime
        time = stime
      END IF

      IF( lform ) READ(outpunit,'(///)')
      CALL read_array( outpunit, p, lform )  ! gas_pressure
!
! ... OLD OUTPUT format had volume fractions here ...
!
!      DO is = 1, nsolid
!
!        IF( lform ) READ(outpunit,'(///)')
!        CALL read_array( outpunit, eps(:,is), lform )  ! solid_bulk_density
!
!      END DO
!
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

      ! ... OLD format
      !DO ig=1,1
!      
! ... NEW output contains all gas species in
!
      DO ig=1,ngas
        IF( lform ) READ(outpunit,'(///)')
        CALL read_array( outpunit, xgc(:,ig), lform )  ! gc_molar_fraction
      END DO

      DO is = 1, nsolid
!
! ... NEW output format has volume fractions here ...
!
        IF( lform ) READ(outpunit,'(///)')
        CALL read_array( outpunit, eps(:,is), lform )  ! solid_bulk_density
!
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
! ... OLD OUTPUT were written starting from the domain top
!
!      DO is = 1, nsolid
!      eptemp = 0.D0
!      DO k = 1, nz
!        DO i = 1, nx
!          ijk = i + (k-1) * nx
!          eptemp(ijk) = eps(i + (nz-k)*nx,is)
!          END DO
!        END DO
!        eps(:,is) = eptemp(:)
!      END DO
!
      RETURN
      END SUBROUTINE read_output
!-----------------------------------------------------------------------
      SUBROUTINE fluxn
! 
      USE dimensions
      USE grid, ONLY: dx, dy, dz, r, xb, yb, zb, rb
      USE io_files, ONLY: tempunit
      USE particles_constants, ONLY: rl

      IMPLICIT NONE
!
      LOGICAL :: lform
      INTEGER :: n1, n2, t1, t2, axis, plane
      INTEGER :: i, j, k, ijk, is, n, t, tn, counter, np
      INTEGER :: ijkm, ijmk, imjk, ipjk, ijpk, ijkp
      REAL*8, ALLOCATABLE :: vel(:), epsm(:)
      REAL*8, ALLOCATABLE :: sflux_p(:,:), sflux_n(:,:)
      REAL*8 :: surface, pi, twopi, totalflux, fluxm
      REAL*8 :: cn1, cn2, ct1, ct2
      TYPE( orthoslice ), ALLOCATABLE :: slice(:)
!
      pi = 4.D0 * ATAN(1.D0)
      twopi = 2.D0 * pi
!
      lform = formatted_output

      ALLOCATE(vel(nsolid), epsm(nsolid))
      vel = 0.D0
      epsm = 0.D0

      ! ... Allocate slices.
      ! ... Read the file containing the indexes or coordinates of the slices
      !
      ALLOCATE(slice(number_of_planes))
      ALLOCATE(sflux_p(number_of_planes,nsolid))
      ALLOCATE(sflux_n(number_of_planes,nsolid))
      slice(:)%axis = 1
      slice(:)%plane = 1
      slice(:)%n1 = 1
      slice(:)%n2 = 1
      slice(:)%t1 = 1
      slice(:)%t2 = 1
!
      CALL allocate_main_fields(ntot)
!
      OPEN(tempunit, FILE=planes_file, STATUS='OLD')
!
      DO np = 1, number_of_planes
                IF (job_type == '2D') THEN
                        READ(tempunit,*) axis, plane
                        READ(tempunit,*) cn1, cn2
                        IF (axis == 1) THEN
                          slice(np)%axis = 1
                          DO i=1,nx
                            IF (xb(i) <= plane) slice(np)%plane = i
                          END DO
                          DO k=1,nz
                            IF (zb(k) <= cn1) slice(np)%n1 = k
                            IF (zb(k) <= cn2) slice(np)%n2 = k
                          END DO
                        ELSE IF (axis == 2) THEN
                          slice(np)%axis = 2
                          DO i=1,nx
                            IF (xb(i) <= cn1) slice(np)%n1 = i
                            IF (xb(i) <= cn2) slice(np)%n2 = i
                          END DO
                          DO k=1,nz
                            IF (zb(k) <= plane) slice(np)%plane= k
                          END DO
                        ELSE
                          CALL error('mass_flux','Invalid axis number',axis)
                        END IF
                ELSE IF (job_type == '3D') THEN
                        READ(tempunit,*) axis, plane
                        READ(tempunit,*) cn1, cn2, ct1, ct2
                        IF (axis == 1) THEN
                          slice(np)%axis = 1
                          DO i=1,nx
                            IF (xb(i) <= plane) slice(np)%plane = i
                          END DO
                          DO j=1,ny
                            IF (yb(j) <= cn1) slice(np)%n1 = j
                            IF (yb(j) <= cn2) slice(np)%n2 = j
                          END DO
                          DO k=1,nz
                            IF (zb(k) <= ct1) slice(np)%t1 = k
                            IF (zb(k) <= ct2) slice(np)%t2 = k
                          END DO
                        ELSE IF (axis == 2) THEN
                          slice(np)%axis = 2
                          DO j=1,ny
                            IF (yb(j) <= plane) slice(np)%plane = j
                          END DO
                          DO k=1,nz
                            IF (zb(k) <= cn1) slice(np)%n1 = k
                            IF (zb(k) <= cn2) slice(np)%n2 = k
                          END DO
                          DO i=1,nx
                            IF (xb(i) <= ct1) slice(np)%t1 = i
                            IF (xb(i) <= ct2) slice(np)%t2 = i
                          END DO
                        ELSE IF (axis == 3) THEN
                          slice(np)%axis = 3
                          DO k=1,nz
                            IF (zb(k) <= plane) slice(np)%plane= k
                          END DO
                          DO i=1,nx
                            IF (xb(i) <= cn1) slice(np)%n1 = i
                            IF (xb(i) <= cn2) slice(np)%n2 = i
                          END DO
                          DO j=1,ny
                            IF (yb(j) <= ct1) slice(np)%t1 = j
                            IF (yb(j) <= ct2) slice(np)%t2 = j
                          END DO
                        ELSE
                          CALL error('mass_flux','Invalid axis number',axis)
                        END IF
                END IF
      END DO
      CLOSE(tempunit)
      WRITE(*,*) 'Slice limits'
      WRITE(*,*) (slice(np)%axis, slice(np)%plane, &
                  slice(np)%n1, slice(np)%n2, slice(np)%t1, slice(np)%t2, &
                  np=1,number_of_planes)
!
      OPEN(tempunit, FILE='massflux.dat', STATUS='UNKNOWN', POSITION='APPEND')
      WRITE(tempunit,*) '****************************'
      DO tn = first_out, last_out, incr_out

        WRITE(logunit,fmt="(/,'* Starting post-processing ',I5,' * ')" ) tn

        ! ... Read PDAC output file
        !
        CALL read_output ( tn )

        ! ... Compute the total solid mass flux across a slice plane
        !
        sflux_p(:,:) = 0.D0
        sflux_n(:,:) = 0.D0
        totalflux = 0.D0
        DO np = 1, number_of_planes
          counter = 0
          axis  = slice(np)%axis
          plane = slice(np)%plane
          n1 = slice(np)%n1
          n2 = slice(np)%n2
          t1 = slice(np)%t1
          t2 = slice(np)%t2
          !
          DO t = t1, t2
            DO n = n1, n2
                IF (job_type == '2D') THEN
                  IF (axis == 1) THEN
                    ijk = plane + (n-1) * nx
                    ipjk = plane + 1 + (n-1) * nx
                    surface = twopi * rb(plane) * dz(n)
                    DO is = 1,nsolid
                      vel(is) = us(ijk,is)
                      epsm(is) = linterp(eps(ijk,is),eps(ipjk,is),dx(plane),dx(plane+1))
                    END DO
                  ELSE IF (axis == 2) THEN
                    ijk = n + (plane-1) * nx
                    ijkp = n + (plane) * nx
                    surface = twopi * r(n) * dx(n)
                    DO is = 1,nsolid
                      vel(is) = ws(ijk,is)
                      epsm(is) = linterp(eps(ijk,is),eps(ijkp,is),dz(plane),dz(plane+1))
                    END DO
                  END IF
                ELSE IF (job_type == '3D') THEN
                  IF (axis == 1) THEN
                    ijk = plane + (n-1) * nx + (t-1) * nx * ny
                    ipjk = plane + 1 + (n-1) * nx + (t-1) * nx * ny
                    surface = dy(n) * dz(t)
                    DO is = 1,nsolid
                      vel(is) = us(ijk,is)
                      epsm(is) = linterp(eps(ijk,is),eps(ipjk,is),dx(plane),dx(plane+1))
                    END DO
                  ELSE IF (axis == 2) THEN
                    ijk = t + (plane-1) * nx + (n-1) * nx * ny
                    ijpk = t + (plane) * nx + (n-1) * nx * ny
                    surface = dz(n) * dx(t)
                    DO is = 1,nsolid
                      vel(is) = vs(ijk,is)
                      epsm(is) = linterp(eps(ijk,is),eps(ijpk,is),dy(plane),dy(plane+1))
                    END DO
                  ELSE IF (axis == 3) THEN
                    ijk = n + (t-1) * nx + (plane-1) * nx * ny
                    ijkp = n + (t-1) * nx + (plane) * nx * ny
                    surface = dx(n) * dy(t)
                    DO is = 1,nsolid
                      vel(is) = ws(ijk,is)
                      epsm(is) = linterp(eps(ijk,is),eps(ijkp,is),dz(plane),dz(plane+1))
                    END DO
                  END IF
                END IF
                !
                counter = counter + 1
                DO is = 1, nsolid
                  fluxm = epsm(is) * rl(is) * vel(is) * surface
                  IF (vel(is) >= 0.D0) THEN
                    sflux_p(np,is)  = sflux_p(np,is) + fluxm
                  ELSE
                    sflux_n(np,is)  = sflux_n(np,is) + fluxm
                  END IF
                END DO
                !
            END DO
          END DO
          !WRITE(*,*) 'Block ', n, ' Counts: ', counter 
        END DO
        WRITE(tempunit,100) &
          time, ((sflux_p(np,is),sflux_n(np,is),is=1,nsolid),np=1,number_of_planes)
      END DO
      CLOSE(tempunit)
 100  FORMAT(20(G30.15E3))
!
      RETURN
      END SUBROUTINE fluxn
!-----------------------------------------------------------------------
      REAL FUNCTION linterp(field1, field2, dx1, dx2)
      IMPLICIT NONE
      REAL, INTENT(IN) :: field1, field2
      REAL*8, INTENT(IN) :: dx1, dx2
      REAL*8 :: dxp
!
      dxp = dx1 + dx2
      linterp = (dx1 * field2 + dx2 * field1) / dxp
! 
      RETURN
      END FUNCTION linterp
!-----------------------------------------------------------------------
      END MODULE mass_orthoflux
!----------------------------------------------------------------------
