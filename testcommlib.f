
program test_comm_lib
  IMPLICIT NONE
  INTEGER :: nproc, mpime, root, group
  INTEGER :: dim, nsiz, i
  INTEGER, ALLOCATABLE :: iarr(:)
  INTEGER :: ichk
  REAL*8, ALLOCATABLE :: darr(:), dbrr(:)
  REAL*8 :: dchk
  REAL*4, ALLOCATABLE :: sarr(:)
  REAL*4 :: schk
  COMPLEX*16, ALLOCATABLE :: zarr(:)
  COMPLEX*16 :: zchk

  CALL startup(nproc, mpime, root, group)

  CALL comm_getsiz( nsiz )
  WRITE(*,*) 'NSIZ = ', nsiz


  !============== check reduce

  dim = INT( DBLE( nsiz ) * 1.3d0 )

  CALL chk_reduce( dim, mpime, nproc )

  dim = INT( DBLE( nsiz ) * 0.3d0 )

  CALL chk_reduce( dim, mpime, nproc )

  dim = nsiz

  CALL chk_reduce( dim, mpime, nproc )


  !============== exit check 

  CALL hangup()

end program

!------------------------------------------------------------------------------!

subroutine chk_reduce( dim, mpime, nproc )
  IMPLICIT NONE
  INTEGER :: nproc, mpime
  INTEGER :: dim, i
  INTEGER, ALLOCATABLE :: iarr(:)
  INTEGER :: ichk
  REAL*8, ALLOCATABLE :: darr(:), dbrr(:)
  REAL*8 :: dchk
  REAL*4, ALLOCATABLE :: sarr(:)
  REAL*4 :: schk
  COMPLEX*16, ALLOCATABLE :: zarr(:)
  COMPLEX*16 :: zchk

  WRITE(*,*) 'DIM  = ', dim

  !================= INTEGER ===============!

  ALLOCATE( iarr( dim ) )

  WRITE(*,*) 'CHECKING parallel_sum_integer'
  iarr = mpime
  CALL parallel_sum_integer( iarr, dim )
  ichk = 0
  DO i = 0, nproc - 1
    ichk = ichk + i
  END DO
  DO i = 1, dim
    IF( iarr( i ) /= ichk ) THEN
      WRITE(*,*) ' ERROR parallel_sum_integer ( mpime, i ) ', mpime, i
    END IF
  END DO 

  WRITE(*,*) 'CHECKING parallel_max_integer'
  iarr = mpime
  CALL parallel_max_integer( iarr, dim )
  ichk = nproc - 1
  DO i = 1, dim
    IF( iarr( i ) /= ichk ) THEN
      WRITE(*,*) ' ERROR parallel_max_integer ( mpime, i ) ', mpime, i
    END IF
  END DO

  WRITE(*,*) 'CHECKING parallel_min_integer'
  iarr = mpime
  CALL parallel_min_integer( iarr, dim )
  ichk = 0
  DO i = 1, dim
    IF( iarr( i ) /= ichk ) THEN
      WRITE(*,*) ' ERROR parallel_min_integer ( mpime, i ) ', mpime, i
    END IF
  END DO

  DEALLOCATE( iarr )

  !================= REAL*8 ===============!

  ALLOCATE( darr( dim ) )
  
  WRITE(*,*) 'CHECKING parallel_sum_real'
  darr = mpime
  CALL parallel_sum_real( darr, dim )
  dchk = 0.0d0
  DO i = 0, nproc - 1
    dchk = dchk + DBLE( i )
  END DO
  DO i = 1, dim
    IF( darr( i ) /= dchk ) THEN
      WRITE(*,*) ' ERROR parallel_sum_real ( mpime, i ) ', mpime, i
    END IF
  END DO

  WRITE(*,*) 'CHECKING parallel_sum_sreal'
  ALLOCATE( dbrr( dim ) )
  darr = mpime
  CALL parallel_sum_real_to( darr, dbrr, dim )
  dchk = 0.0d0
  DO i = 0, nproc - 1
    dchk = dchk + DBLE( i )
  END DO
  DO i = 1, dim
    IF( dbrr( i ) /= dchk ) THEN
      WRITE(*,*) ' ERROR parallel_sum_real_to ( mpime, i ) ', mpime, i
    END IF
  END DO
  DEALLOCATE( dbrr )

  WRITE(*,*) 'CHECKING parallel_max_real'
  darr = mpime
  CALL parallel_max_real( darr, dim )
  dchk = DBLE( nproc - 1 )
  DO i = 1, dim
    IF( darr( i ) /= dchk ) THEN
      WRITE(*,*) ' ERROR parallel_max_real ( mpime, i ) ', mpime, i
    END IF
  END DO

  DEALLOCATE( darr )

  !================= REAL*4 ===============!

  ALLOCATE( sarr( dim ) )

  WRITE(*,*) 'CHECKING parallel_sum_sreal'
  sarr = mpime
  CALL parallel_sum_sreal( sarr, dim )
  schk = 0.0
  DO i = 0, nproc - 1
    schk = schk + REAL( i )
  END DO
  DO i = 1, dim
    IF( sarr( i ) /= schk ) THEN
      WRITE(*,*) ' ERROR parallel_sum_sreal ( mpime, i ) ', mpime, i
    END IF
  END DO

  DEALLOCATE( sarr )

  !================= COMPLEX*16 ===============!

  ALLOCATE( zarr( dim ) )

  WRITE(*,*) 'CHECKING parallel_sum_complex'
  zarr = CMPLX( mpime, mpime )
  CALL parallel_sum_complex( zarr, dim )
  zchk = 0.0
  DO i = 0, nproc - 1
    zchk = zchk + CMPLX( DBLE( i ), DBLE( i ) )
  END DO
  DO i = 1, dim
    IF( zarr( i ) /= zchk ) THEN
      WRITE(*,*) ' ERROR parallel_sum_complex ( mpime, i ) ', mpime, i
    END IF
  END DO

  DEALLOCATE( zarr )

end subroutine

!------------------------------------------------------------------------------!
