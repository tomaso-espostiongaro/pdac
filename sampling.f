!----------------------------------------------------------------------
       MODULE sample_points
       IMPLICIT NONE
       SAVE
!----------------------------------------------------------------------
       CONTAINS
!----------------------------------------------------------------------
       SUBROUTINE set_sampling(xsamp, ysamp, zsamp, ijk, proc)

       USE control_flags, ONLY: job_type
       USE dimensions, ONLY: nx, ny, nz
       USE domain_decomposition, ONLY: cell_owner, cell_g2l
       USE grid, ONLY: x,y,z

       IMPLICIT NONE

       REAL*8, INTENT(IN) :: xsamp, ysamp, zsamp
       INTEGER, INTENT(OUT) :: ijk, proc
       INTEGER :: i,j,k,imesh

       DO i = 1, nx
         IF( x(i) > xsamp ) EXIT
       END DO
       DO k = 1, nz
         IF( z(k) > zsamp ) EXIT
       END DO
       IF (job_type == '3D') THEN
         DO j = 1, ny
           IF( y(j) > ysamp ) EXIT
         END DO
         imesh = i + (j-1) * nx + (k-1) * nx * ny
       ELSE
         imesh = i + (k-1) * nx
       END IF
       proc = cell_owner(imesh)
       ijk  = cell_g2l(imesh,proc)
       
       RETURN
       END SUBROUTINE set_sampling
!----------------------------------------------------------------------
       SUBROUTINE sample_pressure(ijk)
       USE time_parameters, ONLY: time
       USE pressure_epsilon, ONLY: p
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: ijk

       WRITE(19,*) time, p(ijk)

       RETURN
       END SUBROUTINE sample_pressure
!----------------------------------------------------------------------
       END MODULE sample_points
!----------------------------------------------------------------------
