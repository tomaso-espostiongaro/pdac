!----------------------------------------------------------------------
      CHARACTER*4 FUNCTION letter(k)
      IMPLICIT NONE
      CHARACTER ones,tens,hund,thou
!
      INTEGER :: k
!
      INTEGER :: iten, ione, ihund, ithou
!
      ithou=int(k/1000)
      ihund=int((k-(ithou*1000))/100)
      iten=int((k-(ithou*1000)-(ihund*100))/10)
      ione=k-ithou*1000-ihund*100-iten*10
      ones=char(ione+48)
      tens=char(iten+48)
      hund=char(ihunD+48)
      thou=char(ithou+48)
      letter=thou//hund//tens//ones
!
      RETURN
      END FUNCTION 
!----------------------------------------------------------------------
