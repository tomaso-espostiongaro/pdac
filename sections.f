!----------------------------------------------------------------------
      MODULE section_outputs
!----------------------------------------------------------------------
! ... Use ALL independent and derived variables
      USE postp_variables
      IMPLICIT NONE
      PUBLIC
!
! ... parameters for sampling
!
      INTEGER :: number_of_sections
      CHARACTER(LEN=80) :: sect_file
!
      TYPE orthogonal_section
            INTEGER :: nos   ! number (ordinal) of the current section
            INTEGER :: ndim  ! dimension
            INTEGER :: n     ! number (1,2,3) of the axis parallel to the section
            INTEGER :: i     ! index of the position along the axis
            INTEGER :: j
            INTEGER :: k
       END TYPE orthogonal_section
!
      TYPE(orthogonal_section), ALLOCATABLE :: sect(:)
      INTEGER, ALLOCATABLE :: section_indices(:,:)
      REAL*8, ALLOCATABLE :: distance(:)
      INTEGER :: isect
! 
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!---------------------------------------------------------------------
      SUBROUTINE read_sections
!
      USE control_flags, ONLY: job_type, JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions, ONLY: nx, ny, nz, max_size, ngas, nsolid
      USE grid
      USE io_files, ONLY: tempunit, logunit
      USE parallel, ONLY: mpime, root
!
      IMPLICIT NONE
      INTEGER :: nos
!
! ... Allocate sections
!
      ALLOCATE(sect(number_of_sections))
      ALLOCATE(section_indices(number_of_sections,max_size))
      ALLOCATE(distance(max_size))
      distance = 0.D0
      section_indices = 0
!
! ... read sections file 
! 
      IF (mpime == root) THEN
        OPEN(tempunit, FILE=sect_file, STATUS='OLD')
        WRITE(logunit,*) 'Reading sections:'
        WRITE(logunit,*) 'number_of_sections= ', number_of_sections 
      END IF
!
      DO nos = 1, number_of_sections
        sect(nos)%nos = nos
        IF (job_type == JOB_TYPE_2D) THEN 
          IF (mpime == root) READ(tempunit,*) sect(nos)%n, sect(nos)%i, sect(nos)%k
          sect(nos)%j = 0
          CALL bcast_integer(sect(nos)%n,1,root)
          CALL bcast_integer(sect(nos)%i,1,root)
          CALL bcast_integer(sect(nos)%j,1,root)
          CALL bcast_integer(sect(nos)%k,1,root)
        ELSE IF (job_type == JOB_TYPE_3D) THEN
          IF (mpime == root) READ(tempunit,*) sect(nos)%n, sect(nos)%i, sect(nos)%j, sect(nos)%k
          IF (mpime == root) WRITE(*,*) sect(nos)%n, x(sect(nos)%i), y(sect(nos)%j), z(sect(nos)%k)
          CALL bcast_integer(sect(nos)%n,1,root)
          CALL bcast_integer(sect(nos)%i,1,root)
          CALL bcast_integer(sect(nos)%j,1,root)
          CALL bcast_integer(sect(nos)%k,1,root)
        END IF
      END DO
      IF (mpime == root) CLOSE(tempunit)
!
      RETURN
      END SUBROUTINE read_sections
!---------------------------------------------------------------------
      SUBROUTINE section(tn)
!
      USE control_flags, ONLY: job_type, JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions, ONLY: nx, ny, nz, max_size, ngas, nsolid
      USE domain_decomposition, ONLY: cell_owner, cell_g2l
      USE io_files, ONLY: tempunit, logunit
      USE parallel, ONLY: mpime, root
      USE postp_output
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: tn
      INTEGER :: imesh, i, j, k, n, ndim, ijk
      INTEGER :: nos, nv
!
! ... Extract sections from the main fields
!
      DO nos = 1, number_of_sections
        CALL compute_section_indices(nos)
        !
        p_section = 0.D0
        ug_section = 0.D0
        vg_section = 0.D0
        wg_section = 0.D0
        tg_section = 0.D0
        xgc_section = 0.D0
        eps_section = 0.D0
        us_section = 0.D0
        vs_section = 0.D0
        ws_section = 0.D0
        ts_section = 0.D0
        rhom_section = 0.D0
        mvm_section = 0.D0
        pd_section = 0.D0
        !
        DO n = 1, max_size
          IF (cell_owner(section_indices(nos,n)) == mpime) THEN
              ijk = cell_g2l(section_indices(nos,n),mpime)
              p_section(n) = p(ijk)
              ug_section(n) = ug(ijk)
              IF (job_type == JOB_TYPE_3D) vg_section(n) = vg(ijk)
              wg_section(n) = wg(ijk)
              tg_section(n) = tg(ijk)
              DO nv = 1, ngas
                xgc_section(n,nv) = xgc(ijk,nv)
              END DO
              DO nv = 1, nsolid
                eps_section(n,nv) = eps(ijk,nv)
                us_section(n,nv) = us(ijk,nv)
                IF (job_type == JOB_TYPE_3D) vs_section(n,nv) = vs(ijk,nv)
                ws_section(n,nv) = ws(ijk,nv)
                ts_section(n,nv) = ts(ijk,nv)
              END DO
              rhom_section(n) = rhom(ijk)
              mvm_section(n) = mvm(ijk)
              pd_section(n) = pd(ijk)
          END IF
        END DO
        !
        CALL parallel_sum_real(p_section,max_size)
        CALL parallel_sum_real(ug_section,max_size)
        IF (job_type == JOB_TYPE_3D) CALL parallel_sum_real(vg_section,max_size)
        CALL parallel_sum_real(wg_section,max_size)
        CALL parallel_sum_real(tg_section,max_size)
        CALL parallel_sum_real(xgc_section,max_size*ngas)
        CALL parallel_sum_real(eps_section,max_size*nsolid)
        CALL parallel_sum_real(us_section,max_size*nsolid)
        IF (job_type == JOB_TYPE_3D) CALL parallel_sum_real(vs_section,max_size*nsolid)
        CALL parallel_sum_real(ws_section,max_size*nsolid)
        CALL parallel_sum_real(ts_section,max_size*nsolid)
        CALL parallel_sum_real(rhom_section,max_size)
        CALL parallel_sum_real(mvm_section,max_size)
        CALL parallel_sum_real(pd_section,max_size)
        !
        IF (mpime == root) CALL write_section(nos,tn)
      END DO
!
      RETURN
      END SUBROUTINE section
!---------------------------------------------------------------------
      SUBROUTINE compute_section_indices(nos)
      USE control_flags, ONLY: job_type, JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions, ONLY: nx, ny, nz
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nos
      INTEGER :: ndim, n, i, j, k
      !
      IF (nos /= sect(nos)%nos) CALL error('sections','control nos:',nos)
      n = sect(nos)%n
      i =  sect(nos)%i
      j =  sect(nos)%j
      k =  sect(nos)%k
      !
      IF (job_type == JOB_TYPE_2D) THEN
          IF (n==1) THEN
              sect(nos)%ndim = nx
              DO i=1,nx
                section_indices(nos,i) = i + (k-1)*nx
              END DO
          ELSE IF (n==2 .OR. n==3)THEN
              sect(nos)%ndim = nz
              DO k=1,nz
                section_indices(nos,k) = i + (k-1)*nx
              END DO
          END IF
      ELSE IF (job_type == JOB_TYPE_3D) THEN
          IF (n==1) THEN
              sect(nos)%ndim = nx
              DO i=1,nx
                section_indices(nos,i) = i + (j-1)*nx + (k-1)*nx*ny
              END DO
          ELSE IF (n==2)THEN
              sect(nos)%ndim = ny
              DO j=1,ny
                section_indices(nos,j) =  i + (j-1)*nx + (k-1)*nx*ny
              END DO
          ELSE IF (n==3)THEN
              sect(nos)%ndim = nz
              DO k=1,nz
                section_indices(nos,k) =  i + (j-1)*nx + (k-1)*nx*ny
              END DO
          END IF
      END IF
      !
      RETURN
      END SUBROUTINE compute_section_indices
!---------------------------------------------------------------------
      SUBROUTINE write_section(nos,tn)
      USE control_flags, ONLY: job_type, JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions, ONLY: ngas, nsolid, nx, ny, nz
      USE grid, ONLY: x, y, z
      USE io_files, ONLY: tempunit
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nos,tn
      INTEGER :: nn, ndim, n, i, j, k, nv
      CHARACTER(LEN=80) :: sectnam
      CHARACTER(LEN=4) :: lettera
      !
      IF (nos /= sect(nos)%nos) CALL error('sections','control nos:',nos)
      ndim = sect(nos)%ndim
      n = sect(nos)%n
      i =  sect(nos)%i
      j =  sect(nos)%j
      k =  sect(nos)%k
      !
      IF (job_type == JOB_TYPE_2D) THEN
        IF (n==1) THEN
            sectnam = 'section_z'//lettera(k)//'.'//lettera(tn)//'.dat'
            distance(1:nx) = x(1:nx)
        ELSE IF (n==2 .OR. n==3) THEN
            sectnam = 'section_x'//lettera(i)//'.'//lettera(tn)//'.dat'
            distance(1:nz) = z(1:nz)
        END IF
        !
        OPEN(tempunit,FILE=sectnam,STATUS='UNKNOWN')
        DO nn = 1, ndim
          WRITE(tempunit,101) distance(nn), p_section(nn), ug_section(nn), wg_section(nn), &
          tg_section(nn), (xgc_section(nn,nv), nv=1,ngas), (eps_section(nn,nv), &
          us_section(nn,nv), ws_section(nn,nv), ts_section(nn,nv), nv=1,nsolid),&
          rhom_section(nn), mvm_section(nn), pd_section(nn) 
        END DO
      ELSE IF (job_type == JOB_TYPE_3D) THEN
        IF (n==1) THEN
            sectnam = 'section_y'//lettera(j)//'_z'//lettera(k)//'.'//lettera(tn)//'.dat'
            distance(1:nx) = x(1:nx)
        ELSE IF (n==2) THEN
            sectnam = 'section_z'//lettera(k)//'_x'//lettera(i)//'.'//lettera(tn)//'.dat'
            distance(1:ny) = y(1:ny)
        ELSE IF (n==3) THEN
            sectnam = 'section_x'//lettera(i)//'_y'//lettera(j)//'.'//lettera(tn)//'.dat'
            distance(1:nz) = z(1:nz)
        END IF
        !
        OPEN(tempunit,FILE=sectnam,STATUS='UNKNOWN')
        DO nn = 1, ndim
          WRITE(tempunit,101) distance(nn), p_section(nn), ug_section(nn), vg_section(nn), &
          wg_section(nn), tg_section(nn), (xgc_section(nn,nv), nv=1,ngas), &
          (eps_section(nn,nv), us_section(nn,nv), vs_section(nn,nv), ws_section(nn,nv), &
          ts_section(nn,nv), nv=1,nsolid), rhom_section(nn), mvm_section(nn), pd_section(nn)
        END DO
      END IF
 101  FORMAT(F8.2, 100(G14.6E3,1X))
      CLOSE(tempunit)
      !
      RETURN
      END SUBROUTINE write_section
!---------------------------------------------------------------------
      END MODULE section_outputs
!----------------------------------------------------------------------
