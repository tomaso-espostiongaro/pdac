       INTEGER FUNCTION localdim(gdim,np,me)

!   gdim = global dimension of an array to be sheared among processors
!   np   = number of processor
!   me   = index of the calling processor (starting from 0)
!  
!   this function return the number of elements of the array stored
!   in the local memory of the processor "me"

       IMPLICIT NONE
       INTEGER gdim, np, me, r, q

       q = INT(gdim / np)
       r = MOD(gdim, np)

       IF( me .LT. r ) THEN
         localdim = q+1
       ELSE
         localdim = q
       END IF
 
       RETURN
       END FUNCTION
