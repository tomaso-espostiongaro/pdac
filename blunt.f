!----------------------------------------------------------------------
      MODULE blunt_body
!----------------------------------------------------------------------
! ... This module compute the drag and lift forces on a solid block
! ... by assuming that the main stream is flowing in the positive 
! ... x-direction
      USE dimensions
      USE grid, ONLY: iob
      USE parallel, ONLY: mpime, root
      USE pressure_epsilon, ONLY: p, ep
      USE time_parameters, ONLY: time, dt
      USE io_files, ONLY: testunit
      IMPLICIT NONE
!      
      PRIVATE
      PUBLIC :: nblu, ibl, bluntb, set_blunt

      INTEGER :: nblu(1:max_nblock) = 0
      INTEGER :: ibl
      REAL*8, ALLOCATABLE :: fdrag_xp(:), fdrag_yp(:), flift_p(:)
      REAL*8, ALLOCATABLE :: fdrag_xm(:), fdrag_ym(:), flift_m(:)
      REAL*8 :: pd
!
      TYPE surface_pressure
        INTEGER :: np
        INTEGER :: ijk
        INTEGER :: n(3)
        REAL*8  :: ds
        REAL*8  :: p
      END TYPE surface_pressure
      TYPE(surface_pressure), ALLOCATABLE :: surfp(:,:)
      INTEGER :: perim
!
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE set_blunt

      USE control_flags, ONLY: lpr
      USE control_flags, ONLY: job_type, JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions
      USE domain_decomposition, ONLY: cell_g2l, cell_owner
      USE domain_mapping, ONLY: meshinds
      USE grid, ONLY: dx, dy, dz, flag
      USE io_files, ONLY: blunit, blfile

      IMPLICIT NONE
      INTEGER :: prm, l, m, n, i, j, k, ijk, imesh
      INTEGER :: side_x, side_y, side_z, face_x, face_y, face_z
      REAL*8 :: ds
      LOGICAL :: ex

!
! ... 'perim' is the maximum dimension of the blocks surface
! ... (for allocation)
! ... Corner and edges are not included.
!
      IF (mpime == root) OPEN( blunit, FILE=blfile, POSITION='APPEND')
      perim = 0
      prm   = 0
      DO n = 1, no
        IF (nblu(n) == 1) THEN
          side_x = iob(n)%xhi - iob(n)%xlo + 1
          side_y = iob(n)%yhi - iob(n)%ylo + 1
          side_z = iob(n)%zhi - iob(n)%zlo + 1
          face_x = side_y * side_z
          face_y = side_x * side_z
          face_z = side_y * side_x
          IF (job_type == JOB_TYPE_2D) THEN
            prm = 2*( side_x + side_z )
          ELSE IF (job_type == JOB_TYPE_3D) THEN
            prm = 2*( face_x + face_y + face_z )
          END IF
          perim = MAX( prm, perim )
        END IF
      END DO
      IF (ANY( ABS(nblu(:)) > 1 )) CALL error('set_blunt','control blocks',1)
      ALLOCATE(surfp(SUM(nblu), perim))
!
      surfp(:,:)%np   = -9999
      surfp(:,:)%p   = 0.D0
      surfp(:,:)%ds  = 0.D0
      surfp(:,:)%ijk = 0

      m = 0
      blocks_loop: DO n = 1, no
      IF (nblu(n) == 1) THEN
        m = m + 1
        l = 0
        
        ! ... West wall
        !
        i = iob(n)%xlo - 1
        DO k = iob(n)%zlo, iob(n)%zhi
          DO j = iob(n)%ylo, iob(n)%yhi
            l = l+1
            IF (job_type == JOB_TYPE_2D) THEN
              imesh = i + (k-1) * nx
              ds = dz(k)
            ELSE IF (job_type == JOB_TYPE_3D) THEN
              imesh = i + (j-1) * nx + (k-1) * nx * ny
              ds = dz(k) * dy(j)
            END IF
            IF (cell_owner(imesh) == mpime) THEN
              ijk = cell_g2l(imesh,mpime)
              surfp(m,l)%np   = mpime
              surfp(m,l)%ijk  = ijk
              IF ( BTEST(flag(ijk),0) ) surfp(m,l)%ds   = ds
            END IF
            surfp(m,l)%n(1) = +1
            surfp(m,l)%n(2) = 0
            surfp(m,l)%n(3) = 0
          END DO
        END DO

        ! ... Top wall (roof)
        !
        k = iob(n)%zhi + 1
        DO i = iob(n)%xlo, iob(n)%xhi
          DO j = iob(n)%ylo, iob(n)%yhi
            l = l+1
            IF (job_type == JOB_TYPE_2D) THEN
              imesh = i + (k-1) * nx
              ds = dx(i)
            ELSE IF (job_type == JOB_TYPE_3D) THEN
              imesh = i + (j-1) * nx + (k-1) * nx * ny
              ds = dx(i) * dy(j)
            END IF
            IF (cell_owner(imesh) == mpime) THEN
              ijk = cell_g2l(imesh,mpime)
              surfp(m,l)%np   = mpime
              surfp(m,l)%ijk  = ijk
              IF ( BTEST(flag(ijk),0) ) surfp(m,l)%ds   = ds
            END IF
            surfp(m,l)%n(1) = 0
            surfp(m,l)%n(2) = 0
            surfp(m,l)%n(3) = +1
          END DO
        END DO
 
        ! ... East wall
        !
        i = iob(n)%xhi + 1
        DO k = iob(n)%zlo, iob(n)%zhi
          DO j = iob(n)%ylo, iob(n)%yhi
            l = l+1
            IF (job_type == JOB_TYPE_2D) THEN
              imesh = i + (k-1) * nx
              ds = dz(k)
            ELSE IF (job_type == JOB_TYPE_3D) THEN
              imesh = i + (j-1) * nx + (k-1) * nx * ny
              ds = dz(k) * dy(j)
            END IF
            IF (cell_owner(imesh) == mpime) THEN
              ijk = cell_g2l(imesh,mpime)
              surfp(m,l)%np   = mpime
              surfp(m,l)%ijk  = ijk
              IF ( BTEST(flag(ijk),0) ) surfp(m,l)%ds   = ds
            END IF
            surfp(m,l)%n(1) = -1
            surfp(m,l)%n(2) = 0
            surfp(m,l)%n(3) = 0
          END DO
        END DO
        
        ! ... Bottom
        !
        k = iob(n)%zlo - 1
        DO i = iob(n)%xlo, iob(n)%xhi
          DO j = iob(n)%ylo, iob(n)%yhi
            l = l+1
            IF (job_type == JOB_TYPE_2D) THEN
              imesh = i + (k-1) * nx
              ds = dx(i)
            ELSE IF (job_type == JOB_TYPE_3D) THEN
              imesh = i + (j-1) * nx + (k-1) * nx * ny
              ds = dx(i) * dy(j)
            END IF
            IF (cell_owner(imesh) == mpime) THEN
              ijk = cell_g2l(imesh,mpime)
              surfp(m,l)%np   = mpime
              surfp(m,l)%ijk  = ijk
              IF ( BTEST(flag(ijk),0) ) surfp(m,l)%ds   = ds
            END IF
            surfp(m,l)%n(1) = 0
            surfp(m,l)%n(2) = 0
            surfp(m,l)%n(3) = -1
          END DO
        END DO
        
        IF (job_type == JOB_TYPE_3D) THEN
          ! ... South
          !
          j = iob(n)%ylo - 1
          DO i = iob(n)%xlo, iob(n)%xhi 
            DO k = iob(n)%zlo, iob(n)%zhi
              l = l+1
              imesh = i + (j-1) * nx + (k-1) * nx * ny
              ds = dx(i) * dz(k)
              IF (cell_owner(imesh) == mpime) THEN
                ijk = cell_g2l(imesh,mpime)
                surfp(m,l)%np   = mpime
                surfp(m,l)%ijk  = ijk
                IF ( BTEST(flag(ijk),0) ) surfp(m,l)%ds   = ds
              END IF
              surfp(m,l)%n(1) = 0
              surfp(m,l)%n(2) = +1
              surfp(m,l)%n(3) = 0
            END DO
          END DO

          ! ... North
          !
          j = iob(n)%yhi + 1
          DO i = iob(n)%xlo, iob(n)%xhi
            DO k = iob(n)%zlo, iob(n)%zhi
              l = l+1
              imesh = i + (j-1) * nx + (k-1) * nx * ny
              ds = dx(i) * dz(k)
              IF (cell_owner(imesh) == mpime) THEN
                ijk = cell_g2l(imesh,mpime)
                surfp(m,l)%np   = mpime
                surfp(m,l)%ijk  = ijk
                IF ( BTEST(flag(ijk),0) ) surfp(m,l)%ds   = ds
              END IF
              surfp(m,l)%n(1) = 0
              surfp(m,l)%n(2) = -1
              surfp(m,l)%n(3) = 0
            END DO
          END DO
        END IF
!
      END IF
      END DO blocks_loop
!
      m = 0
      DO n = 1, no
        IF (nblu(n) == 1) THEN
          m = m + 1
          IF (lpr > 1) THEN
            WRITE(testunit,*) 
            WRITE(testunit,*) 'Computing action on block: ', n
            WRITE(testunit,*) 'Surface cells: '
            DO l = 1, perim
              IF (surfp(m,l)%np == mpime) THEN
                ijk = surfp(m,l)%ijk
                CALL meshinds(ijk,imesh,i,j,k)
                WRITE(testunit,*) l, i, j, k, ijk
              END IF
            END DO
            WRITE(testunit,*) 'END Computing action on block: ', n
          END IF
        END IF
      END DO
      IF (m /= SUM(nblu)) CALL error('set_blunt','control nblu',m)
!
      RETURN       
      END SUBROUTINE set_blunt
!----------------------------------------------------------------------
      SUBROUTINE bluntb
      USE dimensions, ONLY: no
      USE io_files, ONLY: blunit
!
      IMPLICIT NONE
      INTEGER :: n, m, pp
!
! ... Integrate the gas pressure field around
! ... the body to compute the drag and the lift 
!
      ALLOCATE(fdrag_xp(SUM(nblu)))
      ALLOCATE(fdrag_xm(SUM(nblu)))
      ALLOCATE(fdrag_yp(SUM(nblu)))
      ALLOCATE(fdrag_ym(SUM(nblu)))
      ALLOCATE(flift_p(SUM(nblu)))
      ALLOCATE(flift_m(SUM(nblu)))
!
      fdrag_xp(:) = 0.D0
      fdrag_xm(:) = 0.D0
      fdrag_yp(:) = 0.D0
      fdrag_ym(:) = 0.D0
      flift_p(:) = 0.D0
      flift_m(:) = 0.D0
!
      m = 0
      DO n = 1, no
        IF (nblu(n) == 1) THEN
          m = m + 1
          !
          CALL compute_forces(m, fdrag_xp(m), fdrag_xm(m), fdrag_yp(m), fdrag_ym(m), flift_p(m), flift_m(m))
          !
        END IF
      END DO

      ! ... 'm' is the number of blunt-bodies
      !
      IF (m /= SUM(nblu)) CALL error('set_blunt','control nblu',m)

      CALL parallel_sum_real(fdrag_xp, m)
      CALL parallel_sum_real(fdrag_xm, m)
      CALL parallel_sum_real(fdrag_yp, m)
      CALL parallel_sum_real(fdrag_ym, m)
      CALL parallel_sum_real(flift_p, m)
      CALL parallel_sum_real(flift_m, m)
!
! ... Write out the forces acting on the body
      IF (mpime == root) CALL print_action
!
      DEALLOCATE(fdrag_xp, fdrag_xm, fdrag_yp, fdrag_ym, flift_p, flift_m)

      CALL myflush( blunit )
!
      RETURN
      END SUBROUTINE bluntb
!----------------------------------------------------------------------
      SUBROUTINE compute_forces(m,fdxp,fdxm,fdyp,fdym,fdzp,fdzm)
      USE domain_decomposition, ONLY: cell_g2l, cell_owner
      USE domain_mapping, ONLY: meshinds
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: m
      REAL*8, INTENT(OUT) :: fdxp,fdxm,fdyp,fdym,fdzp,fdzm
      INTEGER :: pp, ijk, n1, n2, n3
      REAL*8 :: ds, pressure, df
!
        fdxp = 0.D0
        fdyp = 0.D0
        fdzp = 0.D0
        fdxm = 0.D0
        fdym = 0.D0
        fdzm = 0.D0
        DO pp = 1, perim
          ijk = surfp(m,pp)%ijk
          IF (mpime == surfp(m,pp)%np) surfp(m,pp)%p = p(ijk)
          ds = surfp(m,pp)%ds
          pressure = surfp(m,pp)%p
          df = ds * pressure
          n1 = surfp(m,pp)%n(1)
          n2 = surfp(m,pp)%n(2)
          n3 = surfp(m,pp)%n(3)
          IF (n1 > 0) THEN
            fdxp = fdxp + df
          ELSE IF (n1 < 0) THEN
            fdxm = fdxm - df
          END IF
          IF (n2 > 0) THEN
            fdyp = fdyp + df
          ELSE IF (n2 < 0) THEN
            fdym = fdym - df
          END IF
          IF (n3 > 0) THEN
            fdzp = fdzp + df
          ELSE IF (n3 < 0) THEN
            fdzp = fdzp - df
          END IF
        END DO
!
      END SUBROUTINE compute_forces
!----------------------------------------------------------------------
      SUBROUTINE print_action

      USE time_parameters, ONLY: time
      USE parallel, ONLY: mpime, root
      USE io_files, ONLY: blunit

      IMPLICIT NONE

      IF( mpime == root ) THEN
        WRITE(blunit,200) time, fdrag_xp(:), fdrag_xm(:), fdrag_xp(:)+fdrag_xm(:), &
                                fdrag_yp(:), fdrag_ym(:), fdrag_yp(:)+fdrag_ym(:), &
                                flift_p(:), flift_m(:), flift_p(:)+flift_m(:)
      END IF

 200  FORMAT(F14.6,100(G18.6))

      END SUBROUTINE print_action

!----------------------------------------------------------------------
      END MODULE blunt_body
!----------------------------------------------------------------------
