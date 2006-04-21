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
      PRIVATE

      PUBLIC :: write_array, read_array, read_inner_array
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

      DEALLOCATE( io_buf )
      DEALLOCATE( io_buf_sgl )
10    FORMAT(1x,10(1x,G14.6E3))
      RETURN
      END SUBROUTINE write_array
!----------------------------------------------------------------------
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

      RETURN
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
      END SUBROUTINE read_inner_array
!----------------------------------------------------------------------
      END MODULE io_parallel
!----------------------------------------------------------------------
