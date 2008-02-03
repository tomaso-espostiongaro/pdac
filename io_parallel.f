!----------------------------------------------------------------------
      MODULE io_parallel
!----------------------------------------------------------------------
      USE dimensions
      USE domain_mapping, ONLY: data_collect, data_distribute
      USE kinds
      USE parallel, ONLY: mpime, root
      USE io_files, ONLY: errorunit
!
      IMPLICIT NONE
!
      PRIVATE
      PUBLIC :: write_array, read_array, read_inner_array, write_crop_array
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE write_array( iunit, array, prec, lform )

      INTEGER, INTENT(IN) :: iunit, prec
      LOGICAL, INTENT(IN) :: lform
      REAL*8 :: array(:)

      INTEGER :: ijk, ij, k, ierr, first
      REAL*8, ALLOCATABLE :: io_buf(:)
      REAL, ALLOCATABLE :: io_buf_sgl(:)

      IF( ntot < 1 ) &
        CALL error(' write_array ', ' ntot too small ', ntot )

      IF( prec /= sgl .AND. prec /= dbl ) &
        CALL error(' write_array ', ' unknown precision ', prec )

      ALLOCATE( io_buf( ntot ), STAT=ierr )
      ALLOCATE( io_buf_sgl( ntot ), STAT=ierr )
      IF( ierr /= 0 ) &
        CALL error(' write_array ', ' cannot allocate io_buf ', ntot )

      CALL data_collect( io_buf, array, 1, ntot )
!
      IF( mpime == root ) THEN
                 IF( lform ) THEN
                   !WRITE(iunit,10) ( io_buf(ijk), ijk = 1, ntot )
                   !
                   ! ... This format is consistent with the old PDAC2D
                   ! ... OUTPUT format
                   DO k = 1, nz
                     first = (k-1) * ntr + 1
                     WRITE(iunit,10) ( io_buf(ij), ij = first, first + ntr - 1 )
                   END DO
                 ELSE
                   IF( prec == sgl ) THEN
                     io_buf_sgl(:) = REAL( io_buf(:), sgl)
                     WRITE(iunit) ( io_buf_sgl(ijk), ijk = 1, ntot )
                     ELSE 
                     WRITE(iunit) ( io_buf(ijk), ijk = 1, ntot )
                   END IF
                 END IF
      END IF
!
      DEALLOCATE( io_buf )
      DEALLOCATE( io_buf_sgl )
10    FORMAT(1x,10(1x,G14.6E3))
      RETURN
      END SUBROUTINE write_array
!----------------------------------------------------------------------
      SUBROUTINE write_crop_array( iunit, array, prec, lform, &
        iminc, imaxc, jminc, jmaxc, kminc, kmaxc)
      USE control_flags, ONLY: job_type
!
      INTEGER, INTENT(IN) :: iunit, prec
      INTEGER, INTENT(INOUT) :: iminc, imaxc, jminc, jmaxc, kminc, kmaxc
      LOGICAL, INTENT(IN) :: lform
      REAL*8 :: array(:)

      INTEGER :: ijk, ij, i, j, k, ierr, first
      INTEGER :: ntotc, ntrc, nxc, nyc, nzc, ijkc
      REAL*8, ALLOCATABLE :: io_buf(:)
      REAL, ALLOCATABLE :: io_buf_sgl(:)
      REAL*8, ALLOCATABLE :: io_buf_crop(:)
      REAL, ALLOCATABLE :: io_buf_sgl_crop(:)

      IF( ntot < 1 ) &
        CALL error(' write_array ', ' ntot too small ', ntot )

      IF( prec /= sgl .AND. prec /= dbl ) &
        CALL error(' write_array ', ' unknown precision ', prec )

      ALLOCATE( io_buf( ntot ), STAT=ierr )
      ALLOCATE( io_buf_sgl( ntot ), STAT=ierr )
      IF( ierr /= 0 ) &
        CALL error(' write_array ', ' cannot allocate io_buf ', ntot )

      CALL data_collect( io_buf, array, 1, ntot )
!
      IF( mpime == root ) THEN
         !
         ! ... Crop the array for post-processing
         !
         IF (iminc < 1) iminc = 1
         IF (imaxc > nx) imaxc = nx
         IF (jminc < 1) jminc = 1
         IF (jmaxc > ny) jmaxc = ny
         IF (kminc < 1) kminc = 1
         IF (kmaxc > nz) kmaxc = nz
         ! ... Bounds are included
         nxc = imaxc - iminc + 1
         nyc = jmaxc - jminc + 1
         nzc = kmaxc - kminc + 1
         !
         ntotc = nxc * nyc * nzc
         IF (job_type == '2D') THEN
                 ntrc = nxc
                 jminc = 1
                 jmaxc = 1
         ELSE IF (job_type == '3D') THEN
                 ntrc = nxc * nyc
         END IF
         ALLOCATE( io_buf_crop(ntotc) )
         ALLOCATE( io_buf_sgl_crop(ntotc) )
         ijkc = 0
         DO k = 1, nz
           DO j = 1, ny
             DO i = 1, nx
               IF (i>=iminc.AND.i<=imaxc) THEN
                 IF (j>=jminc.AND.j<=jmaxc) THEN
                   IF (k>=kminc.AND.k<=kmaxc) THEN
                           IF (job_type == '2D') THEN
                             ijk = i + (k-1) * nx
                           ELSE IF (job_type == '3D') THEN
                             ijk = i + ( (j-1) + (k-1)*ny ) * nx
                           END IF
                           ijkc = ijkc + 1
                           io_buf_crop(ijkc) = io_buf(ijk)
                   END IF
                 END IF
               END IF
             END DO
           END DO
         END DO
!
         IF( lform ) THEN
           DO k = 1, nzc
             first = (k-1) * ntrc + 1
             WRITE(iunit,10) ( io_buf_crop(ij), ij = first, first + ntrc - 1 )
           END DO
         ELSE
           IF( prec == sgl ) THEN
             io_buf_sgl_crop(:) = REAL( io_buf_crop(:), sgl)
             WRITE(iunit) ( io_buf_sgl_crop(ijk), ijk = 1, ntotc )
             ELSE 
             WRITE(iunit) ( io_buf_crop(ijk), ijk = 1, ntotc )
           END IF
         END IF
        DEALLOCATE( io_buf_crop )
        DEALLOCATE( io_buf_sgl_crop )
      END IF
!
      DEALLOCATE( io_buf )
      DEALLOCATE( io_buf_sgl )
10    FORMAT(1x,10(1x,G14.6E3))
      RETURN
      END SUBROUTINE write_crop_array
!--------------------------------------------------------------------

      SUBROUTINE read_array( iunit, array, prec, lform )

      !  This subroutine reads a REAL*8 array ( array ) of
      !  ntot elements from file iunit.
      !  Only root processor read the data, then the array
      !  is distributed across other processors
      !  iunit  (input)  file to be read
      !  array  (output) data read from file and distributed to processors
      !  prec   (input)  precision of the data to be read
      !                  sgl  = single precision data
      !                  dbl  = double precision data
      !  lform  (input)  format of the file 
      !                  .TRUE.  = formatted
      !                  .FALSE. = unformatted

      INTEGER, INTENT(IN) :: iunit, prec
      LOGICAL, INTENT(IN) :: lform
      REAL*8 :: array(:)

      INTEGER :: ijk, ierr
      REAL*8, ALLOCATABLE :: io_buf(:)
      REAL(sgl), ALLOCATABLE :: io_bufs(:)

!
      IF( ntot < 1 ) &
        CALL error(' read_array ', ' ntot too small ', ntot )

      IF( prec /= sgl .AND. prec /= dbl ) &
        CALL error(' read_array ', ' unknown precision ', prec )

      IF( lform ) THEN
         ALLOCATE( io_buf( ntot ), STAT=ierr )
         IF( ierr /= 0 ) &
            CALL error(' read_array ', ' cannot allocate io_buf ', ntot )
         IF( mpime == root ) THEN
            READ(iunit,*) ( io_buf(ijk), ijk = 1, ntot )
         END IF
         CALL data_distribute( io_buf, array, 1, ntot )
         DEALLOCATE( io_buf )
      ELSE
         IF( prec == sgl ) THEN
            ALLOCATE( io_bufs( ntot ), STAT=ierr )
            IF( ierr /= 0 ) &
               CALL error(' read_array ', ' cannot allocate io_buf ', ntot )
            IF( mpime == root ) THEN
               READ(iunit) ( io_bufs(ijk), ijk = 1, ntot )
            END IF
            CALL data_distribute( io_bufs, array, 1, ntot )
            DEALLOCATE( io_bufs )
         ELSE
            ALLOCATE( io_buf( ntot ), STAT=ierr )
            IF( ierr /= 0 ) &
               CALL error(' read_array ', ' cannot allocate io_buf ', ntot )
            IF( mpime == root ) THEN
               READ(iunit) ( io_buf(ijk), ijk = 1, ntot )
            END IF
            CALL data_distribute( io_buf, array, 1, ntot )
            DEALLOCATE( io_buf )
         END IF
      END IF
!
      RETURN
!
!
      END SUBROUTINE read_array
!----------------------------------------------------------------------
      SUBROUTINE read_inner_array( iunit, array, prec, lform )
      USE grid, ONLY: new_ijk, nx_inner, ny_inner, nz_inner

      !  This subroutine reads a REAL*8 array ( array ) of
      !  ntot elements from file iunit.
      !  The array pertains to a 'inner' mesh nested in a new, larger mesh.
      !  The array 'new_ijk' maps the old mesh onto the new one.
      !  Only root processor read the data, then the array
      !  is distributed across other processors
      !  iunit  (input)  file to be read
      !  array  (output) data read from file and distributed to processors
      !  prec   (input)  precision of the data to be read
      !                  sgl  = single precision data
      !                  dbl  = double precision data
      !  lform  (input)  format of the file 
      !                  .TRUE.  = formatted
      !                  .FALSE. = unformatted

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iunit, prec
      LOGICAL, INTENT(IN) :: lform
      REAL*8, INTENT(INOUT) :: array(:)

      INTEGER :: ijk, ierr, ntot_inner
      REAL*8, ALLOCATABLE :: io_buf(:), io_buf_inner(:)
      REAL(sgl), ALLOCATABLE :: io_bufs(:), io_bufs_inner(:)

      ntot_inner = nx_inner * ny_inner * nz_inner

      IF( ntot < 1 ) &
        CALL error(' read_array ', ' ntot too small ', ntot )

      IF( prec /= sgl .AND. prec /= dbl ) &
        CALL error(' read_array ', ' unknown precision ', prec )
!
      IF( lform ) THEN
         ALLOCATE( io_buf( ntot ), STAT=ierr )
         ALLOCATE( io_buf_inner( ntot_inner ), STAT=ierr )
         IF( ierr /= 0 ) &
            CALL error(' read_array ', ' cannot allocate io_buf ', ntot )
         IF( mpime == root ) THEN
            READ(iunit,*) ( io_buf_inner(ijk), ijk = 1, ntot_inner )
         END IF
         !
         ! Assemble global atmospheric initial conditions
         CALL data_collect( io_buf, array, 1, ntot )
         !
         ! Superimpose the initial conditions red from file
         DO ijk = 1, ntot_inner
           io_buf(new_ijk(ijk)) = io_buf_inner(ijk)
         END DO 
         !
         ! Distribute the array after merging
         CALL data_distribute( io_buf, array, 1, ntot )
         !
         DEALLOCATE( io_buf )
         DEALLOCATE( io_buf_inner )
      ELSE
         IF( prec == sgl ) THEN
           ALLOCATE( io_bufs( ntot ), STAT=ierr )
           ALLOCATE( io_bufs_inner( ntot_inner ), STAT=ierr )
           IF( ierr /= 0 ) &
              CALL error(' read_array ', ' cannot allocate io_buf ', ntot )
           IF( mpime == root ) THEN
              READ(iunit) ( io_bufs_inner(ijk), ijk = 1, ntot_inner )
           END IF
           !
           ! Assemble global atmospheric initial conditions
           CALL data_collect( io_bufs, array, 1, ntot )
           !
           ! Superimpose the initial conditions red from file
           DO ijk = 1, ntot_inner
             io_bufs(new_ijk(ijk)) = io_bufs_inner(ijk)
           END DO 
           !
           ! Distribute the array after merging
           CALL data_distribute( io_bufs, array, 1, ntot )
           !
           DEALLOCATE( io_bufs )
           DEALLOCATE( io_bufs_inner )
         ELSE
           ALLOCATE( io_buf( ntot ), STAT=ierr )
           ALLOCATE( io_buf_inner( ntot_inner ), STAT=ierr )
           IF( ierr /= 0 ) &
              CALL error(' read_array ', ' cannot allocate io_buf ', ntot )
           IF( mpime == root ) THEN
              READ(iunit) ( io_buf_inner(ijk), ijk = 1, ntot_inner )
           END IF
           !
           ! Assemble global atmospheric initial conditions
           CALL data_collect( io_buf, array, 1, ntot )
           !
           ! Superimpose the initial conditions red from file
           DO ijk = 1, ntot_inner
             io_buf(new_ijk(ijk)) = io_buf_inner(ijk)
           END DO 
           !
           ! Distribute the array after merging
           CALL data_distribute( io_buf, array, 1, ntot )
           ! 
           DEALLOCATE( io_buf )
           DEALLOCATE( io_buf_inner )
         END IF
      END IF

      RETURN
! 
      END SUBROUTINE read_inner_array
!----------------------------------------------------------------------
      END MODULE io_parallel
!----------------------------------------------------------------------
