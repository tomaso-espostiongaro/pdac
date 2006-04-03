!----------------------------------------------------------------------
      MODULE mass_orthoflux
!----------------------------------------------------------------------
      USE postp_variables, ONLY: time, eps, us, vs, ws
!
      IMPLICIT NONE
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
      INTEGER :: number_of_planes, ifluxn
      CHARACTER(LEN=80) :: planes_file
!
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE fluxn
! 
      USE control_flags, ONLY: job_type
      USE dimensions
      USE grid, ONLY: dx, dy, dz, r, xb, yb, zb, rb
      USE io_files, ONLY: tempunit
      USE particles_constants, ONLY: rl

      IMPLICIT NONE
!
      INTEGER :: n1, n2, t1, t2, axis, plane
      INTEGER :: i, j, k, ijk, is, n, t, counter, np
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
      OPEN(tempunit, FILE='massflux.dat', STATUS='UNKNOWN', &
                                          POSITION='APPEND')
!      WRITE(tempunit,*) '****************************'
!
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
!
      CLOSE(tempunit)
 100  FORMAT(150(G30.15E3))
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
