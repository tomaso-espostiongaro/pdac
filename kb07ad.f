      LOGICAL FUNCTION cpgt(a,b)
      IMPLICIT NONE
      REAL*8 eps
      PARAMETER ( eps = 1.0D-10 )
      REAL*8 a,b,r
      r = abs(a-b)
      IF(r .LT.  eps) THEN
        cpgt = .FALSE.
      ELSE
        cpgt = (a.GT.b)
      END IF
      END

      LOGICAL FUNCTION cplt(a,b)
      IMPLICIT NONE
      REAL*8 eps
      PARAMETER ( eps = 1.0D-10 )
      REAL*8 a,b,r
      r = abs(a-b)
      IF(r .LT.  eps) THEN
        cplt = .FALSE.
      ELSE
        cplt = (a.LT.b)
      END IF
      END


!     ==================================================================
      SUBROUTINE kb07ad(count,n,index)
!     ==--------------------------------------------------------------==
!     == sorting routine for the reciprocal space vectors (g)         ==
!     == kb07ad handles DOuble precision variables                    ==
!     == standard fortran 66 (a verIFied pfort SUBROUTINE)            ==
!     == the work-space 'mark' of length 50 permits up to 2**(50/2)   ==
!     == numbers to be sorted. this is more than the ibm virtual      ==
!     == memory space will hOLD .                                     ==
!     ==--------------------------------------------------------------==
      LOGICAL cpgt,cplt
      REAL*8 count(n),av,x
      
      INTEGER index(n)
      dimension mark(50)
!     ==--------------------------------------------------------------==
!     ==  set index array to original order .                         ==
!     ==--------------------------------------------------------------==
      DO i=1,n
        index(i)=i
      ENDDO
!     ==--------------------------------------------------------------==
!     == check that a trivial CASE has not been entered.              ==
!     ==--------------------------------------------------------------==
      IF(n.EQ.1)GOTO 200
      IF(n.GE.1)GOTO 30
      GOTO 200
!     ==--------------------------------------------------------------==
!     == 'm' is the length of segment which is short enough to enter  ==
!     == the final sorting routine. it may be easily changed.         ==
!     ==--------------------------------------------------------------==
   30 m=12
!     ==--------------------------------------------------------------==
!     == set up initial values.                                       ==
!     ==--------------------------------------------------------------==
      la=2
      is=1
      IF=n
      DO 190 mloop=1,n
!     ==--------------------------------------------------------------==
!     ==  IF segment is short enough sort with final sorting routine. ==
!     ==--------------------------------------------------------------==
        IFka=if-is
        IF((ifka+1).GT.m)GOTO 70
!     ==--------------------------------------------------------------==
!     == final sorting  ( a simple bubble sort )                      ==
!     ==--------------------------------------------------------------==
        is1=is+1
        DO 60 j=is1,IF
          i=j
   40     IF(cplt(count(i-1),count(i)) )GOTO 60
          IF(cpgt(count(i-1),count(i)) )GOTO 50
          IF(index(i-1).LT.index(i))GOTO 60
   50     av=count(i-1)
          count(i-1)=count(i)
          count(i)=av
          int=index(i-1)
          index(i-1)=index(i)
          index(i)=int
          i=i-1
          IF(i.GT.is)GOTO  40
   60   CONTINUE
        la=la-2
        GOTO 170
!     ==--------------------------------------------------------------==
!     ==                *******  quicksort  ********                  ==
!     == SELECT the number in the central position in the segment as  ==
!     == the test number.replace it with the number from the segment's==
!     == highest address.                                             ==
!     ==--------------------------------------------------------------==
   70   iy=(is+IF)/2
        x=count(iy)
        intest=index(iy)
        count(iy)=count(IF)
        index(iy)=index(IF)
!     ==--------------------------------------------------------------==
!     == the markers 'i' and 'IFk' are USEd for the beginning and END ==
!     == of the section not so far tested against the PRESENT value   ==
!     == of x .                                                       ==
!     ==--------------------------------------------------------------==
        k=1
        IFk=if
!     ==--------------------------------------------------------------==
!     == we alternate between the outer loop that increases i and the ==
!     == inner loop that reduces IFk, moving numbers and indices as   ==
!     == necessary, until they meet .                                 ==
!     ==--------------------------------------------------------------==
        DO 110 i=is,IF
          IF(cpgt(x,count(i)))GOTO 110
          IF(cplt(x,count(i)))GOTO 80
          IF(intest.GT.index(i))GOTO 110
   80     IF(i.GE.ifk)GOTO 120
          count(IFk)=count(i)
          index(IFk)=index(i)
          k1=k
          DO 100 k=k1,IFka
            IFk=if-k
            IF(cpgt(count(ifk),x))GOTO 100
            IF(cplt(count(ifk),x))GOTO 90
            IF(intest.LE.index(ifk))GOTO 100
   90       IF(i.GE.ifk)GOTO 130
            count(i)=count(IFk)
            index(i)=index(IFk)
            go to 110
  100     CONTINUE
          GOTO 120
  110   CONTINUE
!     ==--------------------------------------------------------------==
!     == RETURN the test number to the position marked by the marker  ==
!     == which did not move last. it divides the initial segment into ==
!     == 2 parts. any element in the first part is less than or equal ==
!     == to any element in the second part, and they may now be sorted==
!     == indepENDently .                                              ==
!     ==--------------------------------------------------------------==
  120   count(IFk)=x
        index(IFk)=intest
        ip=IFk
        GOTO 140
  130   count(i)=x
        index(i)=intest
        ip=i
!     ==--------------------------------------------------------------==
!     ==  store the longer subdivision in workspace.                  ==
!     ==--------------------------------------------------------------==
  140   IF((ip-is).GT.(if-ip))GOTO 150
        mark(la)=IF
        mark(la-1)=ip+1
        IF=ip-1
        GOTO 160
  150   mark(la)=ip-1
        mark(la-1)=is
        is=ip+1
!     ==--------------------------------------------------------------==
!     == find the length of the shorter subdivision.                  ==
!     ==--------------------------------------------------------------==
  160   lngth=IF-is
        IF(lngth.LE.0)GOTO 180
!     ==--------------------------------------------------------------==
!     == IF it CONTAINS more than one element supply it with workspace==
!     ==--------------------------------------------------------------==
        la=la+2
        GOTO 190
  170   IF(la.LE.0)GOTO 200
!     ==--------------------------------------------------------------==
!     == obtain the address of the shortest segment awaiting quicksort==
!     ==--------------------------------------------------------------==
  180   IF=mark(la)
        is=mark(la-1)
  190 CONTINUE
!     ==--------------------------------------------------------------==
  200 RETURN
      END
!     ==================================================================


!     ==================================================================
      SUBROUTINE ikb07ad(count,n,index)
!     ==--------------------------------------------------------------==
!     == sorting routine for the reciprocal space vectors (g)         ==
!     == kb07ad handles DOuble precision variables                    ==
!     == standard fortran 66 (a verIFied pfort SUBROUTINE)            ==
!     == the work-space 'mark' of length 50 permits up to 2**(50/2)   ==
!     == numbers to be sorted. this is more than the ibm virtual      ==
!     == memory space will hOLD .                                     ==
!     ==--------------------------------------------------------------==
      INTEGER :: count(n),av,x
      
      INTEGER index(n)
      INTEGER mark(50)
!     ==--------------------------------------------------------------==
!     ==  set index array to original order .                         ==
!     ==--------------------------------------------------------------==
      DO i=1,n
        index(i)=i
      ENDDO
!     ==--------------------------------------------------------------==
!     == check that a trivial CASE has not been entered.              ==
!     ==--------------------------------------------------------------==
      IF(n.EQ.1)GOTO 200
      IF(n.GE.1)GOTO 30
      GOTO 200
!     ==--------------------------------------------------------------==
!     == 'm' is the length of segment which is short enough to enter  ==
!     == the final sorting routine. it may be easily changed.         ==
!     ==--------------------------------------------------------------==
   30 m=12
!     ==--------------------------------------------------------------==
!     == set up initial values.                                       ==
!     ==--------------------------------------------------------------==
      la=2
      is=1
      IF=n
      DO 190 mloop=1,n
!     ==--------------------------------------------------------------==
!     ==  IF segment is short enough sort with final sorting routine. ==
!     ==--------------------------------------------------------------==
        IFka=if-is
        IF((ifka+1).GT.m)GOTO 70
!     ==--------------------------------------------------------------==
!     == final sorting  ( a simple bubble sort )                      ==
!     ==--------------------------------------------------------------==
        is1=is+1
        DO 60 j=is1,IF
          i=j
   40     IF((count(i-1).LT.count(i)) )GOTO 60
          IF((count(i-1).GT.count(i)) )GOTO 50
          IF(index(i-1).LT.index(i))GOTO 60
   50     av=count(i-1)
          count(i-1)=count(i)
          count(i)=av
          int=index(i-1)
          index(i-1)=index(i)
          index(i)=int
          i=i-1
          IF(i.GT.is)GOTO  40
   60   CONTINUE
        la=la-2
        GOTO 170
!     ==--------------------------------------------------------------==
!     ==                *******  quicksort  ********                  ==
!     == SELECT the number in the central position in the segment as  ==
!     == the test number.replace it with the number from the segment's==
!     == highest address.                                             ==
!     ==--------------------------------------------------------------==
   70   iy=(is+IF)/2
        x=count(iy)
        intest=index(iy)
        count(iy)=count(IF)
        index(iy)=index(IF)
!     ==--------------------------------------------------------------==
!     == the markers 'i' and 'IFk' are USEd for the beginning and END ==
!     == of the section not so far tested against the PRESENT value   ==
!     == of x .                                                       ==
!     ==--------------------------------------------------------------==
        k=1
        IFk=if
!     ==--------------------------------------------------------------==
!     == we alternate between the outer loop that increases i and the ==
!     == inner loop that reduces IFk, moving numbers and indices as   ==
!     == necessary, until they meet .                                 ==
!     ==--------------------------------------------------------------==
        DO 110 i=is,IF
          IF((x.GT.count(i)))GOTO 110
          IF((x.LT.count(i)))GOTO 80
          IF(intest.GT.index(i))GOTO 110
   80     IF(i.GE.ifk)GOTO 120
          count(IFk)=count(i)
          index(IFk)=index(i)
          k1=k
          DO 100 k=k1,IFka
            IFk=if-k
            IF((count(ifk).GT.x))GOTO 100
            IF((count(ifk).LT.x))GOTO 90
            IF(intest.LE.index(ifk))GOTO 100
   90       IF(i.GE.ifk)GOTO 130
            count(i)=count(IFk)
            index(i)=index(IFk)
            go to 110
  100     CONTINUE
          GOTO 120
  110   CONTINUE
!     ==--------------------------------------------------------------==
!     == RETURN the test number to the position marked by the marker  ==
!     == which did not move last. it divides the initial segment into ==
!     == 2 parts. any element in the first part is less than or equal ==
!     == to any element in the second part, and they may now be sorted==
!     == indepENDently .                                              ==
!     ==--------------------------------------------------------------==
  120   count(IFk)=x
        index(IFk)=intest
        ip=IFk
        GOTO 140
  130   count(i)=x
        index(i)=intest
        ip=i
!     ==--------------------------------------------------------------==
!     ==  store the longer subdivision in workspace.                  ==
!     ==--------------------------------------------------------------==
  140   IF((ip-is).GT.(if-ip))GOTO 150
        mark(la)=IF
        mark(la-1)=ip+1
        IF=ip-1
        GOTO 160
  150   mark(la)=ip-1
        mark(la-1)=is
        is=ip+1
!     ==--------------------------------------------------------------==
!     == find the length of the shorter subdivision.                  ==
!     ==--------------------------------------------------------------==
  160   lngth=IF-is
        IF(lngth.LE.0)GOTO 180
!     ==--------------------------------------------------------------==
!     == IF it CONTAINS more than one element supply it with workspace==
!     ==--------------------------------------------------------------==
        la=la+2
        GOTO 190
  170   IF(la.LE.0)GOTO 200
!     ==--------------------------------------------------------------==
!     == obtain the address of the shortest segment awaiting quicksort==
!     ==--------------------------------------------------------------==
  180   IF=mark(la)
        is=mark(la-1)
  190 CONTINUE
!     ==--------------------------------------------------------------==
  200 RETURN
      END
!     ==================================================================
