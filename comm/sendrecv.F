
! on T3e substitute MPI_DOUBLE_PRECISION with MPI_REAL
! on other architecture substitute MPI_REAL with MPI_DOUBLE_PRECISION

      SUBROUTINE sendrecv_real(sndbuf, snd_size, idest, 
     &  rcvbuf, rcv_size, isour, ip)

        IMPLICIT none
        REAL*8  :: sndbuf(*), rcvbuf(*)
        INTEGER :: snd_size, rcv_size, idest, isour, ip
! ... 
#if defined MPI
        INCLUDE 'mpif.h'
        INTEGER ISTATUS(MPI_STATUS_SIZE), ierr
        CALL MPI_SENDRECV(sndbuf, snd_size, MPI_DOUBLE_PRECISION, 
     &    IDEST, ip, rcvbuf, rcv_size, MPI_DOUBLE_PRECISION, 
     &    ISOUR, ip, MPI_COMM_WORLD, ISTATUS, ierr)
        IF(ierr .NE. 0) THEN
          WRITE(6,*) ' ** ERROR in sendrecv ** '
          STOP
        END IF 
#else
#endif
        RETURN
      END SUBROUTINE



      SUBROUTINE sendrecv_integer(sndbuf, snd_size, idest,
     &  rcvbuf, rcv_size, isour, ip)

        IMPLICIT none
        INTEGER  :: sndbuf(*), rcvbuf(*)
        INTEGER :: snd_size, rcv_size, idest, isour, ip
! ...
#if defined MPI
        INCLUDE 'mpif.h'
        INTEGER ISTATUS(MPI_STATUS_SIZE), ierr
        CALL MPI_SENDRECV(sndbuf, snd_size, MPI_INTEGER,
     &    IDEST, ip, rcvbuf, rcv_size, MPI_INTEGER,
     &    ISOUR, ip, MPI_COMM_WORLD, ISTATUS, ierr)
        IF(ierr .NE. 0) THEN
          WRITE(6,*) ' ** ERROR in sendrecv ** '
          STOP
        END IF
#else
#endif
        RETURN
      END SUBROUTINE

