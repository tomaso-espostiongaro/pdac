!----------------------------------------------------------------------
      CHARACTER*4 FUNCTION lettera(k)
      IMPLICIT NONE
      CHARACTER ones,tens,hund,thou
!
      INTEGER :: k
!
      INTEGER :: iten, ione, ihund, ithou
!
      ithou=INT(k/1000)
      ihund=INT((k-(ithou*1000))/100)
      iten=INT((k-(ithou*1000)-(ihund*100))/10)
      ione=k-ithou*1000-ihund*100-iten*10
      ones=CHAR(ione+48)
      tens=CHAR(iten+48)
      hund=CHAR(ihunD+48)
      thou=CHAR(ithou+48)
      lettera=thou//hund//tens//ones
!
      RETURN
      END FUNCTION 

      CHARACTER*2 FUNCTION lettera2(k)
      IMPLICIT NONE
      CHARACTER ones,tens,hund,thou
!
      INTEGER :: k
!
      INTEGER :: iten, ione, ihund, ithou
!
      ithou=INT(k/1000)
      ihund=INT((k-(ithou*1000))/100)
      iten=INT((k-(ithou*1000)-(ihund*100))/10)
      ione=k-ithou*1000-ihund*100-iten*10
      ones=CHAR(ione+48)
      tens=CHAR(iten+48)
      hund=CHAR(ihunD+48)
      thou=CHAR(ithou+48)
      lettera2=tens//ones
!
      RETURN
      END FUNCTION

!----------------------------------------------------------------------
      CHARACTER*3 FUNCTION procnum(k)
      IMPLICIT NONE
      CHARACTER ones,tens,hund
!
      INTEGER :: k
!
      INTEGER :: iten, ione, ihund
!
      ihund=INT(k/100)
      iten=INT((k-(ihund*100))/10)
      ione=k-ihund*100-iten*10
      ones=CHAR(ione+48)
      tens=CHAR(iten+48)
      hund=CHAR(ihunD+48)
      procnum=hund//tens//ones
!
      RETURN
      END FUNCTION 
!----------------------------------------------------------------------
