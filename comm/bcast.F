!------------------------------------------------------------------------------!

#if defined __INTEL
#define MPI_INTEGER4 MPI_INTEGER
#endif

      SUBROUTINE BCAST_REAL(ARRAY,N,root)
        IMPLICIT NONE
        REAL*8 array(*)
        INTEGER n, root
#if defined __MPI
        INCLUDE 'mpif.h'
        INTEGER ierr 
        CALL MPI_BCAST(array,n,MPI_DOUBLE_PRECISION,root,
     &    MPI_COMM_WORLD,ierr)
#endif
        RETURN
      END SUBROUTINE BCAST_REAL 

      SUBROUTINE BCAST_SREAL(ARRAY,N,root)
        IMPLICIT NONE
        REAL*4 array(*)
        INTEGER n, root
#if defined __MPI
        INCLUDE 'mpif.h'
        INTEGER ierr
        CALL MPI_BCAST(array,n,MPI_REAL,root,
     &    MPI_COMM_WORLD,ierr)
#endif
        RETURN
      END SUBROUTINE BCAST_SREAL


!------------------------------------------------------------------------------!
!
      SUBROUTINE BCAST_INTEGER(ARRAY,N,root)
        IMPLICIT NONE
        INTEGER ARRAY(*), N
        INTEGER root
#if defined __MPI
        INCLUDE 'mpif.h'
        INTEGER IERR
        CALL MPI_BCAST(ARRAY,N,MPI_INTEGER,root,MPI_COMM_WORLD,IERR)
#endif
        RETURN
      END SUBROUTINE BCAST_INTEGER

!
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!

      SUBROUTINE BCAST_CHARACTER(ARRAY,N,root)
        IMPLICIT NONE
        CHARACTER*(*) ARRAY
        INTEGER N, root
#if defined __MPI
        INCLUDE 'mpif.h'
        INTEGER I, IERR, NTMP
        CHARACTER*256 TMP
        INTEGER*4 TMP_CLINK(256/4)
        EQUIVALENCE (TMP,TMP_CLINK) 
        IF(N.GT.256) THEN
          WRITE(*,*)'ROOT_BCAST_CHARACTER SIZE OUT OF RANGE ',N
          CALL hangup
          stop
        END IF
        DO I=1,N
          TMP(I:I) = ARRAY(I:I)
        END DO
        NTMP = N/4
        IF(MOD(N,4).GT.0) THEN
          NTMP = NTMP + 1
        END IF
        CALL MPI_BCAST(TMP_CLINK,NTMP,MPI_INTEGER4,root, 
     &    MPI_COMM_WORLD,IERR)
        DO I=1,N
          ARRAY(I:I) = TMP(I:I)
        END DO
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
        INTEGER IERR,I
        INCLUDE 'mpif.h'
        INTEGER IARRAY(N)
        DO I=1,N
          IF(ARRAY(I)) THEN
            IARRAY(I) = 1
          ELSE
            IARRAY(I) = 0
          END IF
        END DO
        CALL MPI_BCAST(IARRAY,N,MPI_INTEGER,root,MPI_COMM_WORLD,IERR)
        DO I=1,N
          IF(IARRAY(I).EQ.1) THEN
            ARRAY(I) = .TRUE.
          ELSE
            ARRAY(I) = .FALSE.
          END IF
        END DO
#endif
        RETURN
      END SUBROUTINE BCAST_LOGICAL 
