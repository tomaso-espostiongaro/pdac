       INTEGER FUNCTION  MY_PE_ID(IDUM)
         IMPLICIT NONE
#if defined __MPI
         INCLUDE 'mpif.h'
         INTEGER IDUM, RANK, ERR
         CALL MPI_COMM_RANK(MPI_COMM_WORLD,RANK,ERR)
         MY_PE_ID = RANK
#else
         MY_PE_ID = 0
#endif
         RETURN
       END 

