!----------------------------------------------------------------------
      MODULE mass_orthoflux
!----------------------------------------------------------------------
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
      USE domain_decomposition, ONLY: cell_owner, cell_g2l
      USE domain_mapping, ONLY: data_exchange
      USE grid, ONLY: dx, dy, dz, r, xb, yb, zb, rb
      USE io_files, ONLY: tempunit, logunit
      USE parallel, ONLY: mpime, root
      USE particles_constants, ONLY: rl
      USE postp_variables, ONLY: time, eps, us, vs, ws
      USE set_indexes, ONLY: first_subscr, ipjk, ijpk, ijkp

      IMPLICIT NONE
!
      INTEGER :: n1, n2, t1, t2, axis, plane
      INTEGER :: is, n, t, np
      INTEGER :: i, j, k, ijk, imesh
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
!
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
      IF (mpime == root) OPEN(tempunit, FILE=planes_file, STATUS='OLD')
      DO np = 1, number_of_planes
                IF (job_type == '2D') THEN
                        !
                        IF (mpime == root) READ(tempunit,*) axis, plane
                        IF (mpime == root) READ(tempunit,*) cn1, cn2
                        CALL bcast_integer(axis,1,root)
                        CALL bcast_integer(plane,1,root)
                        CALL bcast_real(cn1,1,root)
                        CALL bcast_real(cn2,1,root)
                        !
                        IF (axis == 1) THEN
                          slice(np)%axis = 1
                          DO i=1,nx
                            IF (xb(i) <= plane) slice(np)%plane = i
                          END DO
                          DO k=1,nz
                            IF (zb(k) <= cn1) slice(np)%n1 = MAX(2,k)
                            IF (zb(k) <= cn2) slice(np)%n2 = MIN(nz-1,k)
                          END DO
                        ELSE IF (axis == 2) THEN
                          slice(np)%axis = 2
                          DO i=1,nx
                            IF (xb(i) <= cn1) slice(np)%n1 = MAX(2,i)
                            IF (xb(i) <= cn2) slice(np)%n2 = MIN(nx-1,i)
                          END DO
                          DO k=1,nz
                            IF (zb(k) <= plane) slice(np)%plane= k
                          END DO
                        ELSE
                          CALL error('mass_flux','Invalid axis number',axis)
                        END IF
                ELSE IF (job_type == '3D') THEN
                        !
                        IF (mpime == root) READ(tempunit,*) axis, plane
                        IF (mpime == root) READ(tempunit,*) cn1, cn2, ct1, ct2
                        CALL bcast_integer(axis,1,root)
                        CALL bcast_integer(plane,1,root)
                        CALL bcast_real(cn1,1,root)
                        CALL bcast_real(cn2,1,root)
                        CALL bcast_real(ct1,1,root)
                        CALL bcast_real(ct2,1,root)
                        !
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
      !
      IF (mpime == root) THEN
        CLOSE(tempunit)
        WRITE(logunit,*) 'Slice limits'
        WRITE(logunit,*) (slice(np)%axis, slice(np)%plane, &
                    slice(np)%n1, slice(np)%n2, slice(np)%t1, slice(np)%t2, &
                    np=1,number_of_planes)
      END IF
!
      CALL data_exchange(eps)
!
      ! ... Compute the total solid mass flux across a slice plane
      !
      sflux_p(:,:) = 0.D0
      sflux_n(:,:) = 0.D0
      totalflux = 0.D0
      DO np = 1, number_of_planes
        axis  = slice(np)%axis
        plane = slice(np)%plane
        n1 = slice(np)%n1
        n2 = slice(np)%n2
        t1 = slice(np)%t1
        t2 = slice(np)%t2
        !
        DO t = t1, t2
          DO n = n1, n2
              epsm = 0.D0
              vel  = 0.D0
              IF (job_type == '2D') THEN
                IF (axis == 1) THEN
                  imesh = plane + (n-1) * nx
                  IF (cell_owner(imesh) == mpime) THEN
                    ijk = cell_g2l(imesh,mpime)
                    CALL first_subscr(ijk) ! Identify 'ipjk'
                    surface = twopi * rb(plane) * dz(n)
                    DO is = 1,nsolid
                      vel(is) = us(ijk,is)
                      epsm(is) = linterp(eps(ijk,is),eps(ipjk,is),dx(plane),dx(plane+1))
                    END DO
                  END IF
                ELSE IF (axis == 2) THEN
                  imesh = n + (plane-1) * nx
                  IF (cell_owner(imesh) == mpime) THEN
                    ijk = cell_g2l(imesh,mpime)
                    CALL first_subscr(ijk) ! Identify 'ijkp'
                    surface = twopi * r(n) * dx(n)
                    DO is = 1,nsolid
                      vel(is) = ws(ijk,is)
                      epsm(is) = linterp(eps(ijk,is),eps(ijkp,is),dz(plane),dz(plane+1))
                    END DO
                  END IF
                END IF
              ELSE IF (job_type == '3D') THEN
                IF (axis == 1) THEN
                  imesh = plane + (n-1) * nx + (t-1) * nx * ny
                  IF (cell_owner(imesh) == mpime) THEN
                    ijk = cell_g2l(imesh,mpime)
                    CALL first_subscr(ijk) ! Identify 'ipjk'
                    surface = dy(n) * dz(t)
                    DO is = 1,nsolid
                      vel(is) = us(ijk,is)
                      epsm(is) = linterp(eps(ijk,is),eps(ipjk,is),dx(plane),dx(plane+1))
                    END DO
                  END IF
                ELSE IF (axis == 2) THEN
                  imesh = t + (plane-1) * nx + (n-1) * nx * ny
                  IF (cell_owner(imesh) == mpime) THEN
                    ijk = cell_g2l(imesh,mpime)
                    CALL first_subscr(ijk) ! Identify 'ijpk'
                    surface = dz(n) * dx(t)
                    DO is = 1,nsolid
                      vel(is) = vs(ijk,is)
                      epsm(is) = linterp(eps(ijk,is),eps(ijpk,is),dy(plane),dy(plane+1))
                    END DO
                  END IF
                ELSE IF (axis == 3) THEN
                  imesh = n + (t-1) * nx + (plane-1) * nx * ny
                  IF (cell_owner(imesh) == mpime) THEN
                    ijk = cell_g2l(imesh,mpime)
                    CALL first_subscr(ijk) ! Identify 'ijkp'
                    surface = dx(n) * dy(t)
                    DO is = 1,nsolid
                      vel(is) = ws(ijk,is)
                      epsm(is) = linterp(eps(ijk,is),eps(ijkp,is),dz(plane),dz(plane+1))
                    END DO
                  END IF
                END IF
              END IF
              !
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
      END DO
!
      CALL parallel_sum_real(sflux_p,number_of_planes*nsolid)
      CALL parallel_sum_real(sflux_n,number_of_planes*nsolid)
!
      IF (mpime == root) THEN
        OPEN(tempunit, FILE='massflux.dat', STATUS='UNKNOWN', &
                                          POSITION='APPEND')
        WRITE(tempunit,100) &
          time, ((sflux_p(np,is),sflux_n(np,is),is=1,nsolid),np=1,number_of_planes)
        CLOSE(tempunit)
      END IF
 100  FORMAT(1500(G30.15E3))
!
      RETURN
      END SUBROUTINE fluxn
!-----------------------------------------------------------------------
      REAL FUNCTION linterp(field1, field2, dx1, dx2)
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: field1, field2
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
