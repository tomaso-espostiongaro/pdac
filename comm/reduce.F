
! on T3e substitute MPI_DOUBLE_PRECISION with MPI_REAL
! on other architecture substitute MPI_REAL with MPI_DOUBLE_PRECISION

      SUBROUTINE PARALLEL_SUM_REAL(ARRAY,N)
        IMPLICIT NONE
        INTEGER N
        REAL*8 ARRAY(N)
#if defined __MPI 
        INCLUDE 'mpif.h'
        INCLUDE 'comm.h'
        INTEGER :: ib, nb, ip, nn, err
        REAL*8, ALLOCATABLE :: ARRAY_TMP(:)
        IF( n <= nsiz ) THEN
          ALLOCATE( ARRAY_TMP( n ) )
          CALL MPI_ALLREDUCE( ARRAY, ARRAY_TMP, N,
     &      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ERR )
          ARRAY(1:n) = ARRAY_TMP(1:n)
          DEALLOCATE( ARRAY_TMP )
        ELSE
          ALLOCATE( ARRAY_TMP( nsiz ) )
          nb = n / nsiz
          DO ib = 1, nb
            ip = ( ib - 1 ) * nsiz + 1
            nn = ip + nsiz - 1
            CALL MPI_ALLREDUCE( ARRAY(ip), ARRAY_TMP(1), nsiz, 
     &         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ERR )
            ARRAY(ip:nn) = ARRAY_TMP(1:nsiz)
          END DO
          nn = MOD( n, nsiz )
          IF( nn /= 0 ) THEN
            ip = nb * nsiz + 1
            CALL MPI_ALLREDUCE( ARRAY(ip), ARRAY_TMP(1), nn, 
     &         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ERR )
            ARRAY(ip:n) = ARRAY_TMP(1:nn)
          END IF
          DEALLOCATE( ARRAY_TMP )
        END IF
#else
#endif
        RETURN
      END SUBROUTINE PARALLEL_SUM_REAL


      SUBROUTINE PARALLEL_SUM_SREAL(ARRAY,N)
        IMPLICIT NONE
        INTEGER N
        REAL*4 ARRAY(N)
#if defined __MPI
        INCLUDE 'mpif.h'
        INCLUDE 'comm.h'
        INTEGER :: ib, nb, ip, nn, err
        REAL*4, ALLOCATABLE :: ARRAY_TMP(:)
        IF( n <= nsiz ) THEN
          ALLOCATE( ARRAY_TMP( n ) )
          CALL MPI_ALLREDUCE( ARRAY, ARRAY_TMP, N, MPI_REAL,
     &       MPI_SUM, MPI_COMM_WORLD, ERR)
          ARRAY(1:n) = ARRAY_TMP(1:n)
          DEALLOCATE( ARRAY_TMP )
        ELSE
          ALLOCATE( ARRAY_TMP( nsiz ) )
          nb = n / nsiz
          DO ib = 1, nb
            ip = ( ib - 1 ) * nsiz + 1
            nn = ip + nsiz - 1
            CALL MPI_ALLREDUCE( ARRAY(ip), ARRAY_TMP(1), nsiz,
     &         MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ERR)
            ARRAY(ip:nn) = ARRAY_TMP(1:nsiz)
          END DO
          nn = MOD( n, nsiz )
          IF( nn /= 0 ) THEN
            ip = nb * nsiz + 1
            CALL MPI_ALLREDUCE( ARRAY(ip), ARRAY_TMP(1), nn,
     &         MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ERR)
            ARRAY(ip:n) = ARRAY_TMP(1:nn)
          END IF
          DEALLOCATE( ARRAY_TMP )
        END IF
#else
#endif
        RETURN
      END SUBROUTINE PARALLEL_SUM_SREAL


      SUBROUTINE PARALLEL_SUM_COMPLEX(ARRAY,N)
        IMPLICIT NONE
        INTEGER N
        COMPLEX*16 ARRAY(N)
#if defined __MPI 
        INCLUDE 'mpif.h'
        INCLUDE 'comm.h'
        INTEGER :: ib, nb, ip, nn, err
        COMPLEX*16, ALLOCATABLE :: ARRAY_TMP(:)
        IF( n <= nsiz ) THEN
          ALLOCATE( ARRAY_TMP( n ) )
          CALL MPI_ALLREDUCE(ARRAY, ARRAY_TMP, N, 
     &      MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ERR )
          ARRAY(1:n) = ARRAY_TMP(1:n)
          DEALLOCATE( ARRAY_TMP )
        ELSE
          ALLOCATE( ARRAY_TMP( nsiz ) )
          nb = n / nsiz
          DO ib = 1, nb
            ip = ( ib - 1 ) * nsiz + 1
            nn = ip + nsiz - 1
            CALL MPI_ALLREDUCE( ARRAY(ip), ARRAY_TMP(1), nsiz,
     &         MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ERR )
            ARRAY(ip:nn) = ARRAY_TMP(1:nsiz)
          END DO
          nn = MOD( n, nsiz )
          IF( nn /= 0 ) THEN
            ip = nb * nsiz + 1
            CALL MPI_ALLREDUCE( ARRAY(ip), ARRAY_TMP(1), nn,
     &         MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ERR)
            ARRAY(ip:n) = ARRAY_TMP(1:nn)
          END IF
          DEALLOCATE( ARRAY_TMP )
        END IF
#else
#endif
        RETURN
      END SUBROUTINE PARALLEL_SUM_COMPLEX

!
      SUBROUTINE PARALLEL_SUM_INTEGER(ARRAY,N)
        IMPLICIT NONE
        INTEGER N
        INTEGER ARRAY(N)
#if defined __MPI 
        INCLUDE 'mpif.h'
        INCLUDE 'comm.h'
        INTEGER :: ib, nb, ip, nn, err
        INTEGER, ALLOCATABLE :: ARRAY_TMP(:)
        IF( n <= nsiz ) THEN
          ALLOCATE( ARRAY_TMP( n ) )
          CALL MPI_ALLREDUCE(ARRAY, ARRAY_TMP, N, MPI_INTEGER,
     &       MPI_SUM, MPI_COMM_WORLD, ERR)
          ARRAY(1:n) = ARRAY_TMP(1:n)
          DEALLOCATE( ARRAY_TMP )
        ELSE
          ALLOCATE( ARRAY_TMP( nsiz ) )
          nb = n / nsiz
          DO ib = 1, nb
            ip = ( ib - 1 ) * nsiz + 1
            nn = ip + nsiz - 1
            CALL MPI_ALLREDUCE( ARRAY(ip), ARRAY_TMP(1), nsiz,
     &         MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ERR)
            ARRAY(ip:nn) = ARRAY_TMP(1:nsiz)
          END DO
          nn = MOD( n, nsiz )
          IF( nn /= 0 ) THEN
            ip = nb * nsiz + 1
            CALL MPI_ALLREDUCE( ARRAY(ip), ARRAY_TMP(1), nn,
     &         MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ERR)
            ARRAY(ip:n) = ARRAY_TMP(1:nn)
          END IF
          DEALLOCATE( ARRAY_TMP )
        END IF
#else
#endif
      RETURN
      END SUBROUTINE PARALLEL_SUM_INTEGER
!
!
      SUBROUTINE PARALLEL_SUM_REAL_TO( ARRAY_IN, ARRAY_OUT, N )
        IMPLICIT NONE
        INTEGER N
        REAL*8 ARRAY_IN(N)
        REAL*8 ARRAY_OUT(N)
#if defined __MPI 
        INCLUDE 'mpif.h'
        INCLUDE 'comm.h'
        INTEGER :: ib, nb, ip, nn, err
        IF( n <= nsiz ) THEN
          CALL MPI_ALLREDUCE(ARRAY_IN, ARRAY_OUT, N, 
     &      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ERR)
        ELSE
          nb = n / nsiz
          DO ib = 1, nb
            ip = ( ib - 1 ) * nsiz + 1
            nn = ip + nsiz - 1
            CALL MPI_ALLREDUCE( ARRAY_IN(ip), ARRAY_OUT(ip), nsiz,
     &         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ERR)
          END DO
          nn = MOD( n, nsiz )
          IF( nn /= 0 ) THEN
            ip = nb * nsiz + 1
            CALL MPI_ALLREDUCE( ARRAY_IN(ip), ARRAY_OUT(ip), nn,
     &         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ERR)
          END IF
        END IF
#else
#endif
        RETURN
      END SUBROUTINE PARALLEL_SUM_REAL_TO


      SUBROUTINE parallel_max_integer( array, n )
        IMPLICIT NONE
        INTEGER :: n
        INTEGER :: array( n )
#if defined __MPI 
        INCLUDE 'mpif.h'
        INCLUDE 'comm.h'
        INTEGER :: ib, nb, ip, nn, err
        INTEGER, ALLOCATABLE :: ARRAY_TMP(:)
        IF( n <= nsiz ) THEN
          ALLOCATE( ARRAY_TMP( n ) )
          CALL MPI_ALLREDUCE(ARRAY, ARRAY_TMP, N, MPI_INTEGER,
     &       MPI_MAX, MPI_COMM_WORLD, ERR)
          ARRAY(1:n) = ARRAY_TMP(1:n)
          DEALLOCATE( ARRAY_TMP )
        ELSE
          ALLOCATE( ARRAY_TMP( nsiz ) )
          nb = n / nsiz
          DO ib = 1, nb
            ip = ( ib - 1 ) * nsiz + 1
            nn = ip + nsiz - 1
            CALL MPI_ALLREDUCE( ARRAY(ip), ARRAY_TMP(1), nsiz,
     &         MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ERR)
            ARRAY(ip:nn) = ARRAY_TMP(1:nsiz)
          END DO
          nn = MOD( n, nsiz )
          IF( nn /= 0 ) THEN
            ip = nb * nsiz + 1
            CALL MPI_ALLREDUCE( ARRAY(ip), ARRAY_TMP(1), nn,
     &         MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ERR)
            ARRAY(ip:n) = ARRAY_TMP(1:nn)
          END IF
          DEALLOCATE( ARRAY_TMP )
        END IF
#else
#endif
        RETURN
      END SUBROUTINE PARALLEL_MAX_INTEGER

      SUBROUTINE parallel_min_integer( array, n )
        IMPLICIT NONE
        INTEGER :: n
        INTEGER :: array( n )
#if defined __MPI 
        INCLUDE 'mpif.h'
        INCLUDE 'comm.h'
        INTEGER :: ib, nb, ip, nn, err
        INTEGER, ALLOCATABLE :: ARRAY_TMP(:)
        IF( n <= nsiz ) THEN
          ALLOCATE( ARRAY_TMP( n ) )
          CALL MPI_ALLREDUCE( ARRAY, ARRAY_TMP, N, MPI_INTEGER,
     &       MPI_MIN, MPI_COMM_WORLD, ERR)
          ARRAY(1:n) = ARRAY_TMP(1:n)
          DEALLOCATE( ARRAY_TMP )
        ELSE
          ALLOCATE( ARRAY_TMP( nsiz ) )
          nb = n / nsiz
          DO ib = 1, nb
            ip = ( ib - 1 ) * nsiz + 1
            nn = ip + nsiz - 1
            CALL MPI_ALLREDUCE( ARRAY(ip), ARRAY_TMP(1), nsiz,
     &         MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ERR)
            ARRAY(ip:nn) = ARRAY_TMP(1:nsiz)
          END DO
          nn = MOD( n, nsiz )
          IF( nn /= 0 ) THEN
            ip = nb * nsiz + 1
            CALL MPI_ALLREDUCE( ARRAY(ip), ARRAY_TMP(1), nn,
     &         MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ERR)
            ARRAY(ip:n) = ARRAY_TMP(1:nn)
          END IF
          DEALLOCATE( ARRAY_TMP )
        END IF
#else
#endif
        RETURN
      END SUBROUTINE PARALLEL_MIN_INTEGER


      SUBROUTINE parallel_max_real( array, n )
        IMPLICIT NONE
        INTEGER :: n
        REAL*8  :: array( n )
#if defined __MPI 
        INCLUDE 'mpif.h'
        INCLUDE 'comm.h'
        INTEGER :: ib, nb, ip, nn, err
        REAL*8, ALLOCATABLE :: array_tmp(:)
        IF( n <= nsiz ) THEN
          ALLOCATE( array_tmp( n ) )
          CALL MPI_ALLREDUCE(array, array_tmp, n, 
     &       MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, err )
          array( 1:n ) = array_tmp( 1:n )
          DEALLOCATE( array_tmp )
        ELSE
          ALLOCATE( array_tmp( nsiz ) )
          nb = n / nsiz
          DO ib = 1, nb
            ip = ( ib - 1 ) * nsiz + 1
            nn = ip + nsiz - 1
            CALL MPI_ALLREDUCE( array(ip), array_tmp(1), nsiz,
     &         MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, err )
            array(ip:nn) = array_tmp(1:nsiz)
          END DO
          nn = MOD( n, nsiz )
          IF( nn /= 0 ) THEN
            ip = nb * nsiz + 1
            CALL MPI_ALLREDUCE( array(ip), array_tmp(1), nn,
     &         MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, err)
            array(ip:n) = array_tmp(1:nn)
          END IF
          DEALLOCATE( array_tmp )
        END IF
#else
#endif
        RETURN
      END SUBROUTINE PARALLEL_MAX_REAL
