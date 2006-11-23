!----------------------------------------------------------------------
      MODULE immersed_boundaries
!
! ... Identify the forcing points
! ... Set the geometric parameters characterizing each forcing point
!
!----------------------------------------------------------------------
      USE dimensions, ONLY: nx, ny, nz, ntot, ntr, max_nsolid
      USE grid, ONLY: fl, immb_cell, filled_cell_1, filled_cell_2
      USE parallel, ONLY: mpime, root
      USE volcano_topography, ONLY: xtop, ytop, ztop
      USE volcano_topography, ONLY: topo2d_c, topo2d_x, topo2d_y
      USE volcano_topography, ONLY: topo_c, topo_x
      USE io_files, ONLY: logunit, testunit

      IMPLICIT NONE
!
      TYPE point
        REAL*8 :: x
        REAL*8 :: y
        REAL*8 :: z
      END TYPE point
!
! ... Parameters identifying a forcing point: (i,j,k) are the 
! ... (x,y,z) discrete coordinates, (int) is the type of interpolation
! ... to be done:
! ... In 2D: -3=lin sx/top, -2=linsx, -1=bilsx, 
! ...        0=lintop, 1=bildx, 2=lindx, 3=lin dx/top 
! ... In 3D: always linear interpolation with the nerighbour
! ...        (from 1 to 8) that is more convenient       
! ... (nsl) are the coordinates of the noslip point
!
      TYPE forcing_point
        INTEGER :: i
        INTEGER :: j
        INTEGER :: k
        INTEGER :: int
        INTEGER :: delta_i
        INTEGER :: delta_j
        INTEGER :: delta_k
        INTEGER :: index_q
        INTEGER :: index_qq
        TYPE(point) :: nsl
        REAL*8  :: vel(max_nsolid+1)
      END TYPE forcing_point
!
      TYPE(forcing_point), ALLOCATABLE :: fptx(:), fpty(:), fptz(:)
!
! ... This function gives the increasing number of local forcing points 
      INTEGER, ALLOCATABLE :: numx(:), numy(:), numz(:)
!
! ... This array associates to each local mesh point an integer that
! ... is the decimal representation of the binary number giving the
! ... number of cell faces laying outside the topography
      INTEGER*2, ALLOCATABLE :: bd(:)
      REAL*8, ALLOCATABLE :: vf(:)

      INTEGER :: immb           
      INTEGER :: fp0
!
      PRIVATE
      PUBLIC :: immb, set_forcing
      PUBLIC :: fptx, fpty, fptz, numx, numy, numz, forcing_point, point
      PUBLIC :: bd, vf, faces
!
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE set_forcing
      USE control_flags, ONLY: job_type, lpr
      USE volcano_topography, ONLY: interpolate_profile, interpolate_dem
      USE grid, ONLY: x, xb, y, yb, z, zb
      USE grid, ONLY: inlet_cell, vent_cell
      USE parallel, ONLY: mpime, root
      USE io_files, ONLY: tempunit

      IMPLICIT NONE
      INTEGER :: p, np
      INTEGER :: i,j,k,ijk,imjk,ijmk,ijkm
      INTEGER :: nfpx, nfpy, nfpz
      INTEGER :: bds, counter
!
! ... Conditional array to be used in loops:
! ... this value is .TRUE. when force has to be computed at a given location
!
      LOGICAL, ALLOCATABLE :: forcex(:)
      LOGICAL, ALLOCATABLE :: forcey(:)
      LOGICAL, ALLOCATABLE :: forcez(:)
      LOGICAL, ALLOCATABLE :: extfx(:)
      LOGICAL, ALLOCATABLE :: extfy(:)
      LOGICAL, ALLOCATABLE :: extfz(:)
!
! ... If Immersed Boundaries are used, identify the forcing points
! ... and set interpolation parameters
!
      IF( mpime == root ) THEN
        WRITE(logunit,*)
        WRITE(logunit,*) 'Set forcing point for immersed boundaries'
      END IF

      ! ... Allocate the logical arrays that are used to 
      ! ... identify the forcing points
      !
      ALLOCATE(forcex(ntot)); forcex = .FALSE.
      ALLOCATE(forcey(ntot)); forcey = .FALSE.
      ALLOCATE(forcez(ntot)); forcez = .FALSE.
      ALLOCATE(extfx(ntot)); extfx = .FALSE.
      ALLOCATE(extfy(ntot)); extfy = .FALSE.
      ALLOCATE(extfz(ntot)); extfz = .FALSE.

      IF (job_type == '2D') THEN

        ! ... interpolate the topography on z/x-staggered mesh
        !
        topo_x = 0.D0
        topo_c = 0.D0
        CALL interpolate_profile(xb, z, topo_x, forcex)
        CALL interpolate_profile(x, zb, topo_c, forcez)

        ! ... 2D
        ! ... When all velocity components on the cell faces are
        ! ... forced except one, that component is forced externally
        ! ... and its flag is set to 'filled_cell'
        !
        DO k = 2, nz
          DO i = 2, nx-1
            ijk = (k-1) * nx + i
            imjk = (k-1) * nx + (i-1)
            ijkm = (k-2) * nx + i

            ! ... Add external Forcing in x 
            !
            IF (forcez(ijk) .AND. forcex(ijk) .AND. (z(k) > topo_x(i-1)) ) THEN
                  extfx(imjk) = .TRUE.
                  IF( fl(ijk) /= inlet_cell .AND. fl(ijk) /= vent_cell) fl(ijk) = filled_cell_1
            END IF
            !
            IF (forcez(ijk) .AND. forcex(imjk) .AND. (z(k) > topo_x(i)) ) THEN
                  extfx(ijk) = .TRUE.
                  IF( fl(ijk) /= inlet_cell .AND. fl(ijk) /= vent_cell) fl(ijk) = filled_cell_1
            END IF

            ! ... Add external Forcing in z 
            !
            IF ( forcex(imjk) .AND. forcex(ijk) .AND. forcez(ijkm) .AND. &
                 (zb(k) > topo_c(i)) ) THEN
                 extfz(ijk) = .TRUE.
                 IF( fl(ijk) /= inlet_cell .AND. fl(ijk) /= vent_cell) fl(ijk) = filled_cell_1
            END IF

          END DO
        END DO

        nfpx = COUNT(forcex) + COUNT(extfx) 
        ALLOCATE(fptx(nfpx))
        nfpz = COUNT(forcez) + COUNT(extfz)
        ALLOCATE(fptz(nfpz))

        ! ... Interpolate the topography on x-staggered mesh
        !
        CALL interpolate_profile(xb, z, topo_x, forcex)
        fp0 = 0

        ! ... Forcing along x
        !
        CALL forcing2d(xb, z, topo_x, fptx)

        ! ... External Forcing along x 
        !
        DO k = 2, nz
          DO i = 2, nx-1
            ijk = (k-1) * nx + i
            IF (extfx(ijk)) THEN
              CALL ext_forcing2d(i, k, xb, z, topo_x, fptx)
            END IF
          END DO
        END DO

        ! ... Set flags on forcing points
        !
        DO np = 1, nfpx
          i = fptx(np)%i
          k = fptx(np)%k
          ijk = i + (k-1) * nx
          IF (k>1 .AND. fl(ijk)/=inlet_cell .AND. fl(ijk)/=vent_cell .AND. &
                        fl(ijk)/=filled_cell_1) fl(ijk) = immb_cell
        END DO

        ! ... Interpolate the topography on z-staggered mesh
        !
        CALL interpolate_profile(x, zb, topo_c, forcez)
        fp0 = 0

        ! ... Forcing along z
        !
        CALL forcing2d(x, zb, topo_c, fptz)

        ! ... External Forcing along z
        !
        DO k = 2, nz
          DO i = 2, nx-1
            ijk = (k-1) * nx + i
            IF (extfz(ijk)) THEN
              CALL ext_forcing2d(i, k, x, zb, topo_c, fptz)
            END IF
          END DO
        END DO

        ! ... Set flags on forcing points
        !
        DO np = 1, nfpz
          i = fptz(np)%i
          k = fptz(np)%k
          ijk = i + (k-1) * nx
          IF (k>1 .AND. fl(ijk)/=inlet_cell .AND. fl(ijk)/=vent_cell .AND. &
                        fl(ijk)/=filled_cell_1) fl(ijk) = immb_cell
        END DO
        !
      ELSE IF (job_type == '3D') THEN

        ! ... Interpolate the topography on the y/z/x-staggered mesh
        ! ... to count the forcing points
        !
        topo2d_x = 0.D0
        topo2d_y = 0.D0
        topo2d_c = 0.D0
        CALL interpolate_dem(xb, y, z, topo2d_x, forcex)
        CALL interpolate_dem(x, yb, z, topo2d_y, forcey)
        CALL interpolate_dem(x, y, zb, topo2d_c, forcez)

        ! ... 3D
        ! ... When all velocity components on the cell faces are
        ! ... forced except ONE or TWO, those component are forced externally
        ! ... and their flag is set to 'filled_cell'
        !
        DO k = 2, nz - 1
          DO j = 2, ny - 1
            DO i = 2, nx - 1
              ijk = i + (j-1) * nx + (k-1) * nx * ny
              imjk = (i-1) + (j-1) * nx + (k-1) * nx * ny
              ijk = i + (j-1) * nx + (k-1) * nx * ny
              ijmk = i + (j-2) * nx + (k-1) * nx * ny
              ijk = i + (j-1) * nx + (k-1) * nx * ny
              ijkm = i + (j-1) * nx + (k-2) * nx * ny
 
              !
              ! ... Check the number of faces above the topography
              counter = 0
              bds = 0
              ! 
              ! East
              IF (z(k) > topo2d_x(i,j)) THEN
                    bds = bds + 1
                    counter = counter + 1
              END IF
              !
              ! West
              IF (z(k) > topo2d_x(i-1,j)) THEN
                    bds = bds + 2
                    counter = counter + 1
              END IF
              !
              ! Top
              IF (zb(k) > topo2d_c(i,j))  THEN
                    bds = bds + 4 
                    counter = counter + 1
              END IF
              !
              ! Bottom
              IF (zb(k-1) >= topo2d_c(i,j)) THEN
                   bds = bds + 8
                   counter = counter + 1
              END IF
              !
              ! North
              IF (z(k) > topo2d_y(i,j)) THEN
                   bds = bds + 16
                   counter = counter + 1
              END IF
              !
              ! South
              IF (z(k) > topo2d_y(i,j-1)) THEN
                   bds = bds + 32
                   counter = counter + 1
              END IF
!
              ! ... Add external Forcing in x
              !
              IF (bds == 1) THEN
                    extfx(ijk) = .TRUE.
                    IF( fl(ijk) /= inlet_cell .AND. fl(ijk) /= vent_cell ) fl(ijk) = filled_cell_1
              END IF

              IF (bds == 2) THEN
                    extfx(imjk) = .TRUE.
                    IF( fl(ijk) /= inlet_cell .AND. fl(ijk) /= vent_cell ) fl(ijk) = filled_cell_1
              END IF

              ! ... Add external Forcing in z
              !
              IF (bds == 4) THEN
                    extfz(ijk) = .TRUE.
                    IF( fl(ijk) /= inlet_cell .AND. fl(ijk) /= vent_cell ) fl(ijk) = filled_cell_1
              END IF

              IF (bds == 8) THEN
                    CALL error('set_forcing','inconsistent imm.b.',8)
              END IF

              ! ... Add external Forcing in y 
              !
              IF (bds == 16) THEN
                    extfy(ijk) = .TRUE.
                    IF( fl(ijk) /= inlet_cell .AND. fl(ijk) /= vent_cell ) fl(ijk) = filled_cell_1
              END IF

              IF (bds == 32) THEN
                    extfy(ijmk) = .TRUE.
                    IF( fl(ijk) /= inlet_cell .AND. fl(ijk) /= vent_cell ) fl(ijk) = filled_cell_1
              END IF

              ! ... Add external forcing in x and y
              !
              IF (bds == 17) THEN
                    extfx(ijk) = .TRUE.
                    extfy(ijk) = .TRUE.
                    IF( fl(ijk) /= inlet_cell .AND. fl(ijk) /= vent_cell ) fl(ijk) = filled_cell_2
              END IF
              !
              IF (bds == 18) THEN
                    extfx(imjk) = .TRUE.
                    extfy(ijk) = .TRUE.
                    IF( fl(ijk) /= inlet_cell .AND. fl(ijk) /= vent_cell ) fl(ijk) = filled_cell_2
              END IF
              !
              IF (bds == 34) THEN
                    extfx(imjk) = .TRUE.
                    extfy(ijmk) = .TRUE.
                    IF( fl(ijk) /= inlet_cell .AND. fl(ijk) /= vent_cell ) fl(ijk) = filled_cell_2
              END IF
              !
              IF (bds == 33) THEN
                    extfx(ijk) = .TRUE.
                    extfy(ijmk) = .TRUE.
                    IF( fl(ijk) /= inlet_cell .AND. fl(ijk) /= vent_cell ) fl(ijk) = filled_cell_2
              END IF
              !
              IF (bds == 48) THEN
                    extfy(ijk) = .TRUE.
                    extfy(ijmk) = .TRUE.
                    IF( fl(ijk) /= inlet_cell .AND. fl(ijk) /= vent_cell ) fl(ijk) = filled_cell_2
              END IF
              !
              IF (bds == 3) THEN
                    extfx(ijk) = .TRUE.
                    extfx(imjk) = .TRUE.
                    IF( fl(ijk) /= inlet_cell .AND. fl(ijk) /= vent_cell ) fl(ijk) = filled_cell_2
              END IF

            END DO
          END DO
        ENDDO

        nfpx = COUNT(forcex) + COUNT(extfx)
        ALLOCATE(fptx(nfpx))
        nfpy = COUNT(forcey) + COUNT(extfy)
        ALLOCATE(fpty(nfpy))
        nfpz = COUNT(forcez) + COUNT(extfz)
        ALLOCATE(fptz(nfpz))

        ! ... Interpolate the topography on x-staggered mesh
        !
        CALL interpolate_dem(xb, y, z, topo2d_x, forcex)
        fp0=0

        ! ... External forcing along x
        !
        DO k = 1, nz
          DO j = 1, ny
            DO i = 1, nx
              ijk = i + (j-1) * nx + (k-1) * nx * ny
              IF (extfx(ijk)) THEN
                  CALL ext_forcing3d(i, j, k, xb, y, z, topo2d_x, fptx, 1)
              END IF
            END DO
          END DO
        END DO

        ! ... Forcing along x
        !
        CALL forcing3d(xb, y, z, topo2d_x, fptx)

        ! ... Set flag on forcing points
        !
        DO np = 1, nfpx
          i = fptx(np)%i
          j = fptx(np)%j
          k = fptx(np)%k
          ijk = i + (j-1) * nx + (k-1) * nx * ny
          IF (k>1 .AND. fl(ijk)/=inlet_cell .AND. fl(ijk)/=vent_cell .AND. &
                        fl(ijk)/=filled_cell_1 .AND. fl(ijk)/=filled_cell_2) fl(ijk) = immb_cell
        END DO
        
        ! ... Interpolate the topography on y-staggered mesh.
        !
        CALL interpolate_dem(x, yb, z, topo2d_y, forcey)
        fp0=0 

        ! ... External forcing along y.
        !
        DO k = 1, nz
          DO j = 1, ny
            DO i = 1, nx
              ijk = i + (j-1) * nx + (k-1) * nx * ny
              IF (extfy(ijk)) THEN
                  CALL ext_forcing3d(i, j, k, x, yb, z, topo2d_y, fpty, 2)
              END IF
            END DO
          END DO
        END DO

        ! ... Forcing along y
        !
        CALL forcing3d(x, yb, z, topo2d_y, fpty)

        ! ... Set flag on forcing points
        !
        DO np = 1, nfpy
          i = fpty(np)%i
          j = fpty(np)%j
          k = fpty(np)%k
          ijk = i + (j-1) * nx + (k-1) * nx * ny
          IF (k>1 .AND. fl(ijk)/=inlet_cell .AND. fl(ijk)/=vent_cell .AND. &
                        fl(ijk)/=filled_cell_1 .AND. fl(ijk)/=filled_cell_2) fl(ijk) = immb_cell
        END DO
        
        ! ... Interpolate the topography on z-staggered mesh.
        !
        CALL interpolate_dem(x, y, zb, topo2d_c, forcez)
        fp0=0
        
        ! ... External forcing along z.
        !
        DO k = 1, nz
          DO j = 1, ny
            DO i = 1, nx
              ijk = i + (j-1) * nx + (k-1) * nx * ny
              IF (extfz(ijk)) CALL ext_forcing3d(i, j, k, x, y, zb, topo2d_c, fptz, 3)
            END DO
          END DO
        END DO

        ! ... Forcing along z.
        !
        CALL forcing3d(x, y, zb, topo2d_c, fptz)

        ! ... Set flag on forcing points.
        !
        DO np = 1, nfpz
          i = fptz(np)%i
          j = fptz(np)%j
          k = fptz(np)%k
          ijk = i + (j-1) * nx + (k-1) * nx * ny
          IF (k>1 .AND. fl(ijk)/=inlet_cell .AND. fl(ijk)/=vent_cell .AND. &
                        fl(ijk)/=filled_cell_1 .AND. fl(ijk)/=filled_cell_2) fl(ijk) = immb_cell
        END DO
        
      END IF
!
! ... Write out the forcing points
!
      IF (lpr > 0) THEN
        IF (mpime == root) THEN
          OPEN(UNIT=tempunit,FILE='fptx.dat',STATUS='UNKNOWN')
          DO np = 1, SIZE(fptx)
            WRITE(tempunit,33) np, fptx(np)
          END DO
          CLOSE(tempunit)
          IF (job_type == '3D') THEN
            OPEN(UNIT=tempunit,FILE='fpty.dat',STATUS='UNKNOWN')
            DO np = 1, SIZE(fpty)
              WRITE(tempunit,33) np, fpty(np)
            END DO
            CLOSE(tempunit)
          END IF
          OPEN(UNIT=tempunit,FILE='fptz.dat',STATUS='UNKNOWN')
          DO np = 1, SIZE(fptz)
            WRITE(tempunit,33) np, fptz(np)
          END DO
          CLOSE(tempunit)
        END IF
 33   FORMAT(10(I6),10(F16.3))
      END IF
!
      DEALLOCATE (forcex)
      DEALLOCATE (forcez)
      DEALLOCATE (extfx)
      DEALLOCATE (extfz)
      IF (job_type=='3D') THEN
        DEALLOCATE (forcey)
        DEALLOCATE (extfy)
      END IF

      IF( mpime == root ) THEN
        WRITE(logunit,*) 'END Set forcing'
      END IF
!
      RETURN
      END SUBROUTINE set_forcing
!----------------------------------------------------------------------
      SUBROUTINE ext_forcing2d(i, k, cx, cz, topo, fpt)
      USE volcano_topography, ONLY: ord, next

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: i,k
      REAL*8, DIMENSION(:), INTENT(IN) :: cx, cz
      REAL*8, DIMENSION(:), INTENT(IN) :: topo
      TYPE(forcing_point), DIMENSION(:), INTENT(INOUT) :: fpt
      INTEGER :: f, n
      LOGICAL :: droite, gauche
      REAL*8 :: grad
      droite = .FALSE.
      gauche = .FALSE.

      fp0 = fp0 + 1

      DO f = 1, fp0
        IF (fpt(f)%i==i+1 .AND. fpt(f)%k==k) droite = .TRUE. 
        IF (fpt(f)%i==i-1 .AND. fpt(f)%k==k) gauche = .TRUE.
      END DO
      !
      ! ... Interpolate with the external point on your left
      !
      IF (droite .AND. .NOT.gauche) THEN
        fpt(fp0)%int = 21
        fpt(fp0)%i = i
        fpt(fp0)%j = 0
        fpt(fp0)%k = k
        !
        ! ... find the no-slip point on the profile
        !
        DO n=next(i),next(i+1)
          IF ((ztop(n-1) <= cz(k)).AND.(ztop(n) >= cz(k))) THEN
            grad = (xtop(n)-xtop(n-1))/(ztop(n)-ztop(n-1))
            fpt(fp0)%nsl%x = xtop(n-1) + (cz(k)-ztop(n-1)) * grad
            fpt(fp0)%nsl%z = cz(k)
          ENDIF
        ENDDO
      !
      ! ... Interpolate with the external point on your right
      !
      ELSE IF (.NOT.droite .AND. gauche) THEN
        fpt(fp0)%int = 22
        fpt(fp0)%i = i
        fpt(fp0)%j = 0
        fpt(fp0)%k = k
        !
        ! ... find the no-slip point on the profile
        !
        DO n=next(i-1),next(i)
          IF (ztop(n-1) >= cz(k) .AND. ztop(n) <= cz(k)) THEN
            grad = (xtop(n)-xtop(n-1))/(ztop(n)-ztop(n-1))
            fpt(fp0)%nsl%x = xtop(n-1) + (cz(k)-ztop(n-1)) * grad
            fpt(fp0)%nsl%z = cz(k)
          ENDIF
        ENDDO
      !
      ! ... Interpolate with the external point above
      !
      ELSE 
        fpt(fp0)%int = 20
        fpt(fp0)%i = i
        fpt(fp0)%j = 0
        fpt(fp0)%k = k
        fpt(fp0)%nsl%x = cx(i) 
        fpt(fp0)%nsl%z = topo(i)
      END IF

      RETURN
      END SUBROUTINE ext_forcing2d
!----------------------------------------------------------------------
      SUBROUTINE forcing2d(cx, cz, topo, fpt)
      USE volcano_topography, ONLY: ord, next
!
! ... locate the forcing points and the type of interpolation
! ... locate the points where no-slip condition is applied

      IMPLICIT NONE

      REAL*8, DIMENSION(:), INTENT(IN) :: cx, cz
      REAL*8, DIMENSION(:), INTENT(IN) :: topo
      TYPE(forcing_point), DIMENSION(:), INTENT(OUT) :: fpt

      INTEGER :: i, fp
      REAL*8 :: t,s
!
! ... Skip the -- first grid point -- 
! ... Loop over grid points in x-direction
!
      fp = fp0
      DO i=2,nx-1
         
         IF ((ord(i)-ord(i-1) > 0).AND.(topo(i+1)-topo(i-1) >= 0)) THEN
            
            CALL increasing_profile(i,fp,fpt)
            
         ELSEIF ( ord(i+1)-ord(i) < 0) THEN

            CALL decreasing_profile(i,fp,fpt)

         ELSE 
            
            CALL stationnary_profile(i,fp,fpt)
            
         END IF
         
      END DO
      fp0 = fp
!     
! ... Skip the -- last grid point -- 
!
      RETURN
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE increasing_profile(i,fp,fpt)
      USE control_flags, ONLY: job_type, lpr

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: i
      INTEGER, INTENT(INOUT) :: fp
      TYPE(forcing_point), DIMENSION(:), INTENT(INOUT) :: fpt

      REAL*8 :: norm(2)
      REAL*8 :: grad, s
      REAL*8 :: dxt, dzt, dxs, dzs
      INTEGER :: n, k, err
!
! ... diagonal interpolation 
!
      fp = fp + 1

      fpt(fp)%int = -3
      fpt(fp)%i  = i
      fpt(fp)%k  = ord(i-1)
      
      norm(1) = cx(i-1) - cx(i)
      norm(2) = cz(ord(i-1)+1) - cz(ord(i-1))
!
! ... find the no-slip point on the profile 
!
      DO n=next(i-1),next(i)

        dxt = xtop(n) - xtop(n-1)
        dzt = ztop(n) - ztop(n-1)
        dxs = xtop(n) - cx(i)
        dzs = ztop(n) - cz(ord(i-1))

        s = sol2(norm(1), dxt, norm(2), dzt, dxs, dzs, err)
        IF (err > 0 .AND. lpr > 0) THEN
          WRITE(testunit,*) 'WARNING! from proc: ',mpime
          WRITE(testunit,*) 'Error in fp: ', fp
        END IF

        IF ((s >= 0).AND.(s <= 1)) THEN
          fpt(fp)%nsl%x = s*xtop(n-1) + (1-s)*xtop(n)
          fpt(fp)%nsl%z = s*ztop(n-1) + (1-s)*ztop(n)
        ENDIF
      ENDDO
!
! ... loop over grid points with same x-coordinate
! ... (left linear interpolation)
!
      DO k = ord(i-1) + 1 , ord(i) - 1
        fp = fp + 1
        fpt(fp)%int = -2
        fpt(fp)%i  = i
        fpt(fp)%k  = k
!
! ... find the no-slip point on the profile
!
        DO n=next(i-1),next(i)
          IF (ztop(n-1) <= cz(k) .AND. ztop(n) >= cz(k)) THEN
            grad = (xtop(n)-xtop(n-1))/(ztop(n)-ztop(n-1))
            fpt(fp)%nsl%x = xtop(n-1) + (cz(k)-ztop(n-1)) * grad
            fpt(fp)%nsl%z = cz(k)
          ENDIF
        ENDDO
          
      ENDDO
!!
!! ... vertical linear interpolation 
!!
!      fp = fp + 1
!
!      fpt(fp)%int=0
!
!      fpt(fp)%i  = i
!      fpt(fp)%k  = ord(i)
!
!      fpt(fp)%nsl%x = cx(i)
!      fpt(fp)%nsl%z = topo(i)
!
! ... bilinear interpolation
! ... on the last (upper) point on the same x-coord
!
      fp = fp + 1
      fpt(fp)%int = -1
      fpt(fp)%i  = i
      fpt(fp)%k  = ord(i)

      norm(1) = topo(i-1)-topo(i)
      norm(2) = cx(i)-cx(i-1)
!
! ... find the no-slip point on the profile
!
      DO n=next(i-1),next(i)

        dxt = xtop(n)-xtop(n-1)
        dzt = ztop(n)-ztop(n-1)
        dxs = xtop(n) - cx(i)
        dzs = ztop(n) - cz(ord(i))

        s = sol2( norm(1), dxt, norm(2), dzt, dxs, dzs , err)
        IF (err > 0 .AND. lpr >0) THEN
          WRITE(testunit,*) 'WARNING! from proc: ', mpime
          WRITE(testunit,*) 'Error in fp: ', fp
        END IF

        IF ((s >= 0).AND.(s <= 1))	THEN
          fpt(fp)%nsl%x = s*xtop(n-1) + (1-s)*xtop(n)
          fpt(fp)%nsl%z = s*ztop(n-1) + (1-s)*ztop(n)
        ENDIF
      ENDDO
      
      RETURN
      END SUBROUTINE increasing_profile
!----------------------------------------------------------------------
      SUBROUTINE decreasing_profile(i,fp,fpt)
      USE control_flags, ONLY: lpr

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: i
      INTEGER, INTENT(INOUT) :: fp
      TYPE(forcing_point), DIMENSION(:), INTENT(INOUT) :: fpt

      REAL*8 :: norm(2)
      REAL*8 :: grad, s
      REAL*8 :: dxt, dzt, dxs, dzs
      INTEGER :: n, k, err
!
! ... bilinear interpolation
! ... on the first (upper) point
!
      fp = fp + 1

      fpt(fp)%int=1
      fpt(fp)%i  = i
      fpt(fp)%k  = ord(i)

      norm(1) = topo(i) - topo(i+1)
      norm(2) = cx(i+1) - cx(i)
!
! ... find the no-slip point on the profile
!
      DO n=next(i),next(i+1)

        dxt = xtop(n) - xtop(n-1)
        dzt = ztop(n) - ztop(n-1)
        dxs = xtop(n) - cx(i)
        dzs = ztop(n) - cz(ord(i))

        s = sol2( norm(1), dxt, norm(2), dzt, dxs, dzs, err )
        IF (err > 0 .AND. lpr >0) THEN
          WRITE(testunit,*) 'WARNING! from proc: ', mpime
          WRITE(testunit,*) 'Error in fp: ', fp
        END IF

        IF ((s >= 0).AND.(s <= 1)) THEN
          fpt(fp)%nsl%x = s*xtop(n-1) + (1-s)*xtop(n)
          fpt(fp)%nsl%z = s*ztop(n-1) + (1-s)*ztop(n)
        ENDIF
      ENDDO
!!
!! ... vertical linear interpolation 
!!
!      fp = fp + 1
!
!      fpt(fp)%int=0
!
!      fpt(fp)%i  = i
!      fpt(fp)%k  = ord(i)
!
!      fpt(fp)%nsl%x = cx(i)
!      fpt(fp)%nsl%z = topo(i)
!
! ... loop over grid points with same x-coordinate
! ... (right linear interpolation)
!
      DO k=ord(i)-1,ord(i+1)+1,-1

       fp = fp + 1

       fpt(fp)%int = 2
       fpt(fp)%i  = i
       fpt(fp)%k  = k
!
! ... find the no-slip point on the profile
!
        DO n=next(i),next(i+1)
          IF ((ztop(n-1) >= cz(k)).AND.(ztop(n) <= cz(k))) THEN
            grad = (xtop(n)-xtop(n-1))/(ztop(n)-ztop(n-1))
            fpt(fp)%nsl%x = xtop(n-1) + (cz(k)-ztop(n-1)) * grad
            fpt(fp)%nsl%z = cz(k)
          ENDIF
        ENDDO

      ENDDO
!
! ... diagonal interpolation 
!
      fp = fp + 1

      fpt(fp)%int = 3
      fpt(fp)%i  = i
      fpt(fp)%k  = ord(i+1)

      norm(1)=cx(i+1)-cx(i)
      norm(2)=cz(ord(i+1)+1)-cz(ord(i+1))
!
! ... find the no-slip point on the profile
!
      DO n=next(i),next(i+1)
 
        dxt = xtop(n) - xtop(n-1)
        dzt = ztop(n) - ztop(n-1)
        dxs = xtop(n) - cx(i)
        dzs = ztop(n) - cz(ord(i+1))

	s = sol2( norm(1), dxt, norm(2), dzt, dxs, dzs, err )
        IF (err > 0 .AND. lpr >0) THEN
          WRITE(testunit,*) 'WARNING! from proc: ', mpime
          WRITE(testunit,*) 'Error in fp: ', fp
        END IF

 	IF ((s >= 0).AND.(s <= 1))	THEN
	  fpt(fp)%nsl%x = s*xtop(n-1) + (1-s)*xtop(n)
	  fpt(fp)%nsl%z = s*ztop(n-1) + (1-s)*ztop(n)
 	ENDIF
      ENDDO
      
      END SUBROUTINE decreasing_profile
!----------------------------------------------------------------------
      SUBROUTINE stationnary_profile(i,fp,fpt)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: i
      INTEGER, INTENT(INOUT) :: fp
      TYPE(forcing_point), DIMENSION(:), INTENT(INOUT) :: fpt

      fp = fp + 1

      fpt(fp)%int=0	

      fpt(fp)%i  = i
      fpt(fp)%k  = ord(i)

      fpt(fp)%nsl%x = cx(i)
      fpt(fp)%nsl%z = topo(i)

      END SUBROUTINE stationnary_profile
!----------------------------------------------------------------------
      FUNCTION sol1(a11,a12,a21,a22,b1,b2,info)

      REAL*8 a11,a12,a21,a22	! coefficienti della matrice
      REAL*8 b1,b2		! vettore dei termini noti
      REAL*8 sol1
      INTEGER, INTENT(OUT) :: info

      info = 0
      IF (a11*a22-a12*a21 == 0) info = 1
      sol1=(b1*a22-a12*b2)/(a11*a22-a12*a21)

      END FUNCTION sol1
!----------------------------------------------------------------------
      FUNCTION sol2(a11,a12,a21,a22,b1,b2, info)

      REAL*8 a11,a12,a21,a22	! coefficienti della matrice
      REAL*8 b1,b2		! vettore dei termini noti
      REAL*8 sol2
      INTEGER, INTENT(OUT) :: info
 
      info = 0
      IF (a11*a22-a12*a21 == 0) info = 1
      sol2=(a11*b2-b1*a21)/(a11*a22-a12*a21)

      END FUNCTION sol2
!----------------------------------------------------------------------
!
      END SUBROUTINE forcing2d
!
!----------------------------------------------------------------------
      SUBROUTINE forcing3d(cx, cy, cz, topo2d, fpt)
      USE volcano_topography, ONLY: ord2d
!
! ... locate the forcing points and the type of interpolation
! ... locate the points where no-slip condition is applied
!
      IMPLICIT NONE

      REAL*8, DIMENSION(:), INTENT(IN) :: cx, cy, cz
      REAL*8, DIMENSION(:,:), INTENT(IN) :: topo2d
      TYPE(forcing_point), DIMENSION(:), INTENT(OUT) :: fpt

      INTEGER :: i, j, k, fp
      REAL*8 :: g(8)
      INTEGER :: delta_i(8), delta_j(8)
      INTEGER :: gint
!
      fp = fp0
!
! ... LINEAR interpolation criteria ('fpt%int'):
! ... 0: top
! ... 1: west
! ... 2: south-west
! ... 3: south
! ... 4: south-east
! ... 5: east
! ... 6: north-east
! ... 7: north
! ... 8: north-west
!
      delta_i(1) = -1
      delta_j(1) = 0
      delta_i(2) = -1
      delta_j(2) = -1
      delta_i(3) = 0
      delta_j(3) = -1
      delta_i(4) = 1
      delta_j(4) = -1
      delta_i(5) = 1
      delta_j(5) = 0
      delta_i(6) = 1
      delta_j(6) = 1
      delta_i(7) = 0
      delta_j(7) = 1
      delta_i(8) = -1
      delta_j(8) = 1

! ... Skip the boundary points
!
      DO i = 2, nx - 1
        DO j = 2, ny - 1
          DO k = 2, ord2d(i,j) - 1
            g = -1
            IF( ord2d(i-1,j)   < k ) g(1) = gamma(i,j,k, 1)
            IF( ord2d(i-1,j-1) < k ) g(2) = gamma(i,j,k, 2)
            IF( ord2d(i,j-1)   < k ) g(3) = gamma(i,j,k, 3)
            IF( ord2d(i+1,j-1) < k ) g(4) = gamma(i,j,k, 4)
            IF( ord2d(i+1,j)   < k ) g(5) = gamma(i,j,k, 5)
            IF( ord2d(i+1,j+1) < k ) g(6) = gamma(i,j,k, 6)
            IF( ord2d(i,j+1)   < k ) g(7) = gamma(i,j,k, 7)
            IF( ord2d(i-1,j+1) < k ) g(8) = gamma(i,j,k, 8)

            IF (MAXVAL(g) > -1) THEN
              !
              ! ... Choose the external point for interpolation
              gint = mingam(g) 

              fp = fp + 1 
              fpt(fp)%i  = i
              fpt(fp)%j  = j
              fpt(fp)%k  = k
              fpt(fp)%int = gint
              fpt(fp)%delta_i = delta_i(gint)
              fpt(fp)%delta_j = delta_j(gint)
              fpt(fp)%delta_k = 0
              fpt(fp)%nsl%x = g(gint) * cx(i) + &
                              (1.D0 - g(gint)) * cx(i+delta_i(gint))
              fpt(fp)%nsl%y = g(gint) * cy(j) + &
                              (1.D0 - g(gint)) * cy(j+delta_j(gint))
              fpt(fp)%nsl%z = cz(k)

            END IF

          END DO
          
          k = ord2d(i,j)
          
          fp = fp + 1 
          fpt(fp)%i  = i
          fpt(fp)%j  = j
          fpt(fp)%k  = k
          fpt(fp)%int = 0
          fpt(fp)%delta_i = 0
          fpt(fp)%delta_j = 0
          fpt(fp)%delta_k = 1
          fpt(fp)%nsl%x = cx(i)
          fpt(fp)%nsl%y = cy(j)
          fpt(fp)%nsl%z = topo2d(i,j)

        END DO
      END DO
      RETURN
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      REAL*8 FUNCTION gamma(i_,j_,k_,index)
!
! ... 0 < gamma < 1
! ... accordingly to the relative horizontal position between
! ... immersed point/topographic point/external point
! ... The topography between two mesh points is interpolated
! ... linearly.
! 
      IMPLICIT NONE
      ! ... (1 < index < 8) according to the relative position of
      ! ... the neighbours
      INTEGER, INTENT(IN) :: i_, j_, k_, index

      INTEGER :: di_, dj_
      REAL*8  :: num, den

      di_ = delta_i(index)
      dj_ = delta_j(index)

      num = cz(k_) - topo2d(i_,j_)
      den = topo2d(i_+di_, j_+dj_) - topo2d(i_,j_)

      gamma = num / den
      
      RETURN
      END FUNCTION gamma
!----------------------------------------------------------------------
      INTEGER FUNCTION mingam(gam)
!
! ... Gives a criterion for the choice of the external interpolation point.
! ... Choose the external point such as the no-slip points lays as close as possible
! ... to the mid-point ('a=0.5')  between the immersed and the external point.
!
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: gam(:)
      REAL*8 :: a, ming
      INTEGER :: n

      a = 0.5D0
      ming = 1
      
      DO n = 1, 8
        IF( ABS(a-gam(n)) < ming ) THEN
          ming = ABS(a-gam(n))
          mingam = n
        END IF
      END DO
      
      RETURN
      END FUNCTION mingam
!----------------------------------------------------------------------
      END SUBROUTINE forcing3d
!----------------------------------------------------------------------
      SUBROUTINE ext_forcing3d(i, j, k, cx, cy, cz, topo2d, fpt, axis)

      USE volcano_topography, ONLY: ord2d
!
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: i, j, k
      INTEGER, INTENT(IN) :: axis
       
      REAL*8, DIMENSION(:), INTENT(IN) :: cx, cy, cz
      REAL*8, DIMENSION(:,:), INTENT(IN) :: topo2d
      TYPE(forcing_point), DIMENSION(:), INTENT(OUT) :: fpt
 
      REAL*8 :: g(8), g_EXT(4), g_min, maxg, maxg_z
      INTEGER :: delta_i(8), delta_j(8)
      INTEGER :: l,gint,fpint

      delta_i(1) = -1
      delta_j(1) = 0
      delta_i(2) = -1
      delta_j(2) = -1
      delta_i(3) = 0
      delta_j(3) = -1
      delta_i(4) = 1
      delta_j(4) = -1
      delta_i(5) = 1
      delta_j(5) = 0
      delta_i(6) = 1
      delta_j(6) = 1
      delta_i(7) = 0
      delta_j(7) = 1
      delta_i(8) = -1
      delta_j(8) = 1

      g_EXT = 0.D0
      g = -1
      maxg = 1.D0

      IF ( k > ord2d(i,j) ) THEN

         IF( ord2d(i-1,j)   >= k ) g(1) = gamma(i,j,k, 1)
         IF( ord2d(i-1,j-1) >= k ) g(2) = gamma(i,j,k, 2)
         IF( ord2d(i,j-1)   >= k ) g(3) = gamma(i,j,k, 3)
         IF( ord2d(i+1,j-1) >= k ) g(4) = gamma(i,j,k, 4)
         IF( ord2d(i+1,j)   >= k ) g(5) = gamma(i,j,k, 5)
         IF( ord2d(i+1,j+1) >= k ) g(6) = gamma(i,j,k, 6)
         IF( ord2d(i,j+1)   >= k ) g(7) = gamma(i,j,k, 7)
         IF( ord2d(i-1,j+1) >= k ) g(8) = gamma(i,j,k, 8)
!
! ... If one of this products is < 0 then the external forcing
! ... can be applied along the corresponding direction
!
         IF (axis /= 1) g_EXT(1) = g(1) * g(5)      !  E-W
         g_EXT(2) = g(2) * g(6)                     !  SE-NW
         IF (axis /= 2) g_EXT(3) = g(3) * g(7)      !  S-N
         g_EXT(4) = g(4) * g(8)                     !  SW-NE
!
! ... Select the direction which minimize the distance 
! ... between the forcing point and the no-slip point
!
         DO l = 1,4
            IF (( g_EXT(l) < 0 ) .AND. (-g_EXT(l) < maxg))  THEN
               g_min = l
               maxg = - g_EXT(l)
            ENDIF
         ENDDO
!
! ... Select the versus along the direction previously chosen
!
         IF ( MINVAL (g_EXT) < 0 ) THEN
            IF ( g_min == 1 ) THEN

               IF ( g(1) < 0 ) THEN
                  gint = 5
                  fpint = 21
               ELSE
                  gint = 1
                  fpint = 25
               ENDIF

            ELSEIF ( g_min == 2 ) THEN

               IF ( g(2) < 0 ) THEN
                  gint = 6
                  fpint = 22
               ELSE
                  gint = 2
                  fpint = 26
               ENDIF

            ELSEIF ( g_min == 3 ) THEN

               IF ( g(3) < 0 ) THEN
                  gint = 7
                  fpint = 23
               ELSE
                  gint = 3
                  fpint = 27
               ENDIF

            ELSEIF ( g_min == 4 ) THEN

               IF ( g(4) < 0 ) THEN
                  gint = 8
                  fpint = 24
               ELSE
                  gint = 4
                  fpint = 28
               ENDIF

            ENDIF
         ENDIF
      ENDIF

      IF ( k == ord2d(i,j)+1 ) THEN
         maxg_z = ( cz(k) - topo2d(i,j) ) / ( cz(k+1) - topo2d(i,j) )
      ELSE
         maxg_z = 1.D0
      END IF

      IF ( (maxg < 1) .AND. (maxg < maxg_Z) ) THEN
!
!
! ... Choose the external point for interpolation
! ... in the z-plane of the forcing point 
!         
         fp0 = fp0 + 1 
         fpt(fp0)%i  = i
         fpt(fp0)%j  = j
         fpt(fp0)%k  = k
         fpt(fp0)%int = fpint
         fpt(fp0)%delta_i = delta_i(fpint - 20)
         fpt(fp0)%delta_j = delta_j(fpint - 20)
         fpt(fp0)%delta_k = 0
         fpt(fp0)%nsl%x = (1.D0 - g(gint)) * cx(i) + &
                          g(gint) * cx(i+delta_i(gint))
         fpt(fp0)%nsl%y = (1.D0 - g(gint)) * cy(j) + &
                          g(gint) * cy(j+delta_j(gint))
         fpt(fp0)%nsl%z = cz(k)
         
      ELSE
!
! ... Choose the external forcing point for interpolation
! ... over the forcing point
!
         fp0 = fp0 + 1 
         fpt(fp0)%i  = i
         fpt(fp0)%j  = j
         fpt(fp0)%k  = k
         fpt(fp0)%int = 20
         fpt(fp0)%delta_i = 0
         fpt(fp0)%delta_j = 0
         fpt(fp0)%delta_k = 1
         fpt(fp0)%nsl%x = cx(i)
         fpt(fp0)%nsl%y = cy(j)
         fpt(fp0)%nsl%z = topo2d(i,j)
         !
      ENDIF
      
      RETURN
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      REAL*8 FUNCTION gamma(i_,j_,k_,index)
!
! ... 0 < gamma < 1
! ... accordingly to the relative horizontal position between
! ... immersed point/topographic point/external point
! ... The topography between two mesh points is interpolated
! ... linearly.
! 
      IMPLICIT NONE
      ! ... (1 < index < 8) according to the relative position of
      ! ... the neighbours
      INTEGER, INTENT(IN) :: i_, j_, k_, index

      INTEGER :: di_, dj_
      REAL*8  :: num, den

      di_ = delta_i(index)
      dj_ = delta_j(index)

      num = cz(k_) - topo2d(i_,j_)
      den = topo2d(i_+di_, j_+dj_) - topo2d(i_,j_)

      gamma = num / den
      
      RETURN
      END FUNCTION gamma
!----------------------------------------------------------------------
      END SUBROUTINE ext_forcing3d
!----------------------------------------------------------------------
      SUBROUTINE faces(ijk, b_e, b_w, b_t, b_b, b_n, b_s, ivf)
      USE control_flags, ONLY: job_type
 
      IMPLICIT NONE
 
      INTEGER, INTENT(IN) :: ijk
      INTEGER, INTENT(OUT) :: b_e, b_w, b_t, b_b, b_n, b_s
      REAL*8, INTENT(OUT) :: ivf
      INTEGER*2 :: num
 
      IF (numx(ijk)/=0 .AND. numz(ijk) /=0) THEN
        b_e = 0
        b_w = 0
        b_t = 0
        b_b = 0
        RETURN
      END IF 

      b_e = 0
      num = 1
      IF( IAND(bd(ijk),num) /= 0 )  b_e = 1 
 
      b_w = 0
      num = 2
      IF( IAND(bd(ijk),num) /= 0 )  b_w = 1
 
      b_t = 0
      num = 4
      IF( IAND(bd(ijk),num) /= 0 )  b_t = 1
 
      b_b = 0
      num = 8
      IF( IAND(bd(ijk),num) /= 0 )  b_b = 1
      
      b_n = 0
      num = 16
      IF( IAND(bd(ijk),num) /= 0 ) b_n = 1
 
      b_s = 0
      num = 32
      IF( IAND(bd(ijk),num) /= 0 ) b_s = 1
 
      IF (vf(ijk) > 0.D0) THEN
        ivf = 1.D0 / vf(ijk)
      ELSE
        ivf = 0.D0
      END IF
 
      RETURN
      END SUBROUTINE faces
!----------------------------------------------------------------------
      END MODULE immersed_boundaries
!----------------------------------------------------------------------
