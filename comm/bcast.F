!------------------------------------------------------------------------------!

#if defined __INTEL
#define MPI_INTEGER4 MPI_INTEGER
#endif

      SUBROUTINE BCAST_REAL( ARRAY, N, root )
        IMPLICIT NONE
        REAL*8 array(*)
        INTEGER :: n, root
#if defined __MPI
        INCLUDE 'mpif.h'
        INCLUDE 'comm.h'
        INTEGER :: ib, nb, ip, nn, err
        IF( n <= nsiz ) THEN
          CALL MPI_BCAST( array, n, MPI_DOUBLE_PRECISION, root,
     &      MPI_COMM_WORLD, err)
        ELSE
          nb = n / nsiz
          DO ib = 1, nb
            ip = ( ib - 1 ) * nsiz + 1
            nn = ip + nsiz - 1
            CALL MPI_BCAST( array(ip), nsiz, MPI_DOUBLE_PRECISION, root,
     &        MPI_COMM_WORLD, err)
          END DO
          nn = MOD( n, nsiz )
          IF( nn /= 0 ) THEN
            ip = nb * nsiz + 1
            CALL MPI_BCAST( array(ip), nn, MPI_DOUBLE_PRECISION, root,
     &        MPI_COMM_WORLD, err)
          END IF
        END IF
#endif
        RETURN
      END SUBROUTINE BCAST_REAL 

      SUBROUTINE BCAST_SREAL(ARRAY,N,root)
        IMPLICIT NONE
        REAL*4  :: array(*)
        INTEGER :: n, root
#if defined __MPI
        INCLUDE 'mpif.h'
        INCLUDE 'comm.h'
        INTEGER :: ib, nb, ip, nn, err
        IF( n <= nsiz ) THEN
          CALL MPI_BCAST( array, n, MPI_REAL, root,
     &      MPI_COMM_WORLD, err)
        ELSE
          nb = n / nsiz
          DO ib = 1, nb
            ip = ( ib - 1 ) * nsiz + 1
            nn = ip + nsiz - 1
            CALL MPI_BCAST( array(ip), nsiz, MPI_REAL, root,
     &        MPI_COMM_WORLD, err)
          END DO
          nn = MOD( n, nsiz )
          IF( nn /= 0 ) THEN
            ip = nb * nsiz + 1
            CALL MPI_BCAST( array(ip), nn, MPI_REAL, root,
     &        MPI_COMM_WORLD, err)
          END IF
        END IF
#endif
        RETURN
      END SUBROUTINE BCAST_SREAL


!------------------------------------------------------------------------------!
!
      SUBROUTINE BCAST_INTEGER(ARRAY,N,root)
        IMPLICIT NONE
        INTEGER :: ARRAY(*), N
        INTEGER :: root
#if defined __MPI
        INCLUDE 'mpif.h'
        INCLUDE 'comm.h'
        INTEGER :: ib, nb, ip, nn, err
        IF( n <= nsiz ) THEN
          CALL MPI_BCAST( array, n, MPI_INTEGER, root,
     &      MPI_COMM_WORLD, err)
        ELSE
          nb = n / nsiz
          DO ib = 1, nb
            ip = ( ib - 1 ) * nsiz + 1
            nn = ip + nsiz - 1
            CALL MPI_BCAST( array(ip), nsiz, MPI_INTEGER, root,
     &        MPI_COMM_WORLD, err)
          END DO
          nn = MOD( n, nsiz )
          IF( nn /= 0 ) THEN
            ip = nb * nsiz + 1
            CALL MPI_BCAST( array(ip), nn, MPI_INTEGER, root,
     &        MPI_COMM_WORLD, err)
          END IF
        END IF

#endif
        RETURN
      END SUBROUTINE BCAST_INTEGER

!
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!

      SUBROUTINE BCAST_CHARACTER( ARRAY, N, root )
        IMPLICIT NONE
        CHARACTER(LEN=*) :: ARRAY
        INTEGER          :: N, root
#if defined __MPI
        INTEGER :: IERR, I, nn
        INTEGER, ALLOCATABLE :: IARRAY(:)
        nn = LEN( array )
        ALLOCATE( iarray( nn ) )
        DO I=1,nn
          IARRAY(I) = ICHAR( array( i:i ) )
        END DO
        CALL bcast_integer( iarray, nn, root )
        DO I=1,nn
            ARRAY(i:i) = CHAR( iarray( i ) )
        END DO
        DEALLOCATE( iarray )
#endif
        RETURN
      END  SUBROUTINE BCAST_CHARACTER  

!
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!

      SUBROUTINE BCAST_LOGICAL(ARRAY,N,root)
        IMPLICIT NONE
        LOGICAL ARRAY(*)
        INTEGER N
        INTEGER root
#if defined __MPI
        INTEGER :: IERR, I
        INTEGER, ALLOCATABLE :: IARRAY(:)
        ALLOCATE( iarray( n ) )
        DO I=1,N
          IF(ARRAY(I)) THEN
            IARRAY(I) = 1
          ELSE
            IARRAY(I) = 0
          END IF
        END DO
        CALL bcast_integer( iarray, n, root )
        DO I=1,N
          IF(IARRAY(I).EQ.1) THEN
            ARRAY(I) = .TRUE.
          ELSE
            ARRAY(I) = .FALSE.
          END IF
        END DO
        DEALLOCATE( iarray )
#endif
        RETURN
      END SUBROUTINE BCAST_LOGICAL 
