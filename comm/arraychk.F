       SUBROUTINE array_chk( a, n, info )
         INTEGER :: n, info
         REAL*8 :: a( n )
         INTEGER :: float_chk
         INTEGER :: i
         info = 0
         DO i = 1, n
           IF( float_chk( a( i ) ) /= 0 ) THEN
             WRITE(6,*) 'array_chk: ',i,a(i)
             a( i ) = 0.0d0
             info = info + 1
           END IF
         END DO
       END SUBROUTINE
