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
      REAL*8, ALLOCATABLE :: fdrag(:), flift(:)
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
      USE dimensions
      USE domain_decomposition, ONLY: cell_g2l, cell_owner, meshinds
      USE grid, ONLY: dx, dz, flag
      USE io_files, ONLY: blunit, blfile

      IMPLICIT NONE
      INTEGER :: prm, l, m, n, i, j, k, ijk, imesh
      LOGICAL :: ex

      IF (mpime == root) OPEN( blunit, FILE=blfile, POSITION='APPEND')
      perim = 0
      prm   = 0
      DO n = 1, no
        IF (nblu(n) == 1) THEN
          prm = 2*( iob(n)%zhi - iob(n)%zlo + iob(n)%xhi - iob(n)%xlo) + 4
          perim = MAX( prm, perim )
        END IF
      END DO

      IF (ANY( ABS(nblu(:)) > 1 )) CALL error('set_blunt','control blocks',1)
      ALLOCATE(surfp(SUM(nblu), perim))

      surfp(:,:)%p   = 0.D0
      surfp(:,:)%ds  = 0.D0
      surfp(:,:)%ijk = 0

      m = 0
      DO n = 1, no
      IF (nblu(n) == 1) THEN
        m = m + 1
        l = 0
        
        ! ... West wall
        !
        i = iob(n)%xlo - 1
        DO k = iob(n)%zlo, iob(n)%zhi
          l = l+1
          imesh = i + (k-1) * nx
          IF (cell_owner(imesh) == mpime) THEN
            ijk = cell_g2l(imesh,mpime)
            surfp(m,l)%np   = mpime
            surfp(m,l)%ijk  = ijk
            IF ( BTEST(flag(ijk),0) ) surfp(m,l)%ds   = dz(k)
          END IF
          surfp(m,l)%n(1) = +1
          surfp(m,l)%n(2) = 0
          surfp(m,l)%n(3) = 0
        END DO
        
        ! ... Roof
        !
        k = iob(n)%zhi + 1
        DO i = iob(n)%xlo, iob(n)%xhi
          l = l+1
          imesh = i + (k-1) * nx
          IF (cell_owner(imesh) == mpime) THEN
            ijk = cell_g2l(imesh,mpime)
            surfp(m,l)%np   = mpime
            surfp(m,l)%ijk  = ijk
            IF ( BTEST(flag(ijk),0) ) surfp(m,l)%ds   = dx(i)
          END IF
          surfp(m,l)%n(1) = 0
          surfp(m,l)%n(2) = 0
          surfp(m,l)%n(3) = +1
        END DO
        
        ! ... East wall
        !
        i = iob(n)%xhi + 1
        DO k = iob(n)%zhi, iob(n)%zlo, -1
          l = l+1
          imesh = i + (k-1) * nx
          IF (cell_owner(imesh) == mpime) THEN
            ijk = cell_g2l(imesh,mpime)
            surfp(m,l)%np   = mpime
            surfp(m,l)%ijk  = ijk
            IF ( BTEST(flag(ijk),0) ) surfp(m,l)%ds   = dz(k)
          END IF
          surfp(m,l)%n(1) = -1
          surfp(m,l)%n(2) = 0
          surfp(m,l)%n(3) = 0
        END DO
        
        ! ... Bottom
        !
        k = iob(n)%zlo - 1
        DO i = iob(n)%xhi, iob(n)%xlo, -1
          l = l+1
          imesh = i + (k-1) * nx
          IF (cell_owner(imesh) == mpime) THEN
            ijk = cell_g2l(imesh,mpime)
            surfp(m,l)%np   = mpime
            surfp(m,l)%ijk  = ijk
            IF ( BTEST(flag(ijk),0) ) surfp(m,l)%ds   = dx(i)
          END IF
          surfp(m,l)%n(1) = 0
          surfp(m,l)%n(2) = 0
          surfp(m,l)%n(3) = -1
        END DO

      END IF
      END DO
!
      IF (m /= SUM(nblu)) CALL error('set_blunt','control nblu',m)
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
                WRITE(testunit,*) i, j, k, ijk
              END IF
            END DO
            WRITE(testunit,*) 'END Computing action on block: ', n
          END IF
        END IF
      END DO
        
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
      ALLOCATE(fdrag(SUM(nblu)))
      ALLOCATE(flift(SUM(nblu)))

      fdrag(:) = 0.D0
      flift(:) = 0.D0

      m = 0
      DO n = 1, no
        IF (nblu(n) == 1) THEN
          m = m + 1

          CALL dragbl(m, fdrag(m))
          CALL liftbl(m, flift(m))

        END IF
      END DO

      ! ... 'm' is the number of blunt-bodies
      !
      IF (m /= SUM(nblu)) CALL error('set_blunt','control nblu',m)

      CALL parallel_sum_real(fdrag, m)
      CALL parallel_sum_real(flift, m)
!
! ... Write out the forces acting on the body
      IF (mpime == root) CALL print_action
!
      DEALLOCATE(fdrag, flift)

      CALL myflush( blunit )
!
      RETURN
      END SUBROUTINE bluntb
!----------------------------------------------------------------------
      SUBROUTINE dragbl(m,fd)
      USE domain_decomposition, ONLY: cell_g2l, cell_owner, meshinds
      IMPLICIT NONE
!
        INTEGER, INTENT(IN) :: m
	REAL*8, INTENT(OUT) :: fd
        INTEGER :: pp, ijk
!
        fd = 0.D0
        DO pp = 1, perim
          ijk = surfp(m,pp)%ijk
          IF (mpime == surfp(m,pp)%np) surfp(m,pp)%p = p(ijk)
          fd = fd + surfp(m,pp)%ds * surfp(m,pp)%p * surfp(m,pp)%n(1)
          !WRITE(testunit,*) ijk, surfp(m,pp)%p, ep(ijk)
	END DO
!
      END SUBROUTINE dragbl
!----------------------------------------------------------------------
      SUBROUTINE liftbl(m,fl)
      IMPLICIT NONE
!
        INTEGER, INTENT(IN) :: m
	REAL*8, INTENT(OUT) :: fl
        INTEGER :: pp, ijk
!
        fl = 0.D0
        DO pp = 1, perim
          ijk = surfp(m,pp)%ijk
          IF (mpime == surfp(m,pp)%np) surfp(m,pp)%p = p(ijk)
          fl = fl + surfp(m,pp)%ds * surfp(m,pp)%p * surfp(m,pp)%n(3)
	END DO
!
      END SUBROUTINE liftbl

!----------------------------------------------------------------------

      SUBROUTINE print_action

      USE time_parameters, ONLY: time
      USE parallel, ONLY: mpime, root
      USE io_files, ONLY: blunit

      IMPLICIT NONE

      IF( mpime == root ) WRITE(blunit,200) time, fdrag(:), flift(:)

 200  FORMAT(F14.6,20(G18.6))

      END SUBROUTINE print_action

!----------------------------------------------------------------------
      END MODULE blunt_body
!----------------------------------------------------------------------
