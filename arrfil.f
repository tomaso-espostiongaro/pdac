!----------------------------------------------------------------------
      MODULE array_filters
      IMPLICIT NONE
!
      INTERFACE interp
        MODULE PROCEDURE interp_1d_array, interp_1d_scalar, interp_2d
      END INTERFACE interp
!
      INTERFACE filter
        MODULE PROCEDURE filter_1d, filter_2d
      END INTERFACE filter
!
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
! ... Filter out high frequency modes by successively subsampling,
! ... averaging and interpolating 
!
      SUBROUTINE filter_1d(x1,f1,nsm)

      IMPLICIT NONE
      REAL*8, INTENT(INOUT), DIMENSION(:) :: x1, f1
      INTEGER, INTENT(IN) :: nsm
      REAL*8, ALLOCATABLE, DIMENSION(:) :: x2, f2
      INTEGER, ALLOCATABLE :: dummy(:)
      REAL*8 :: sstep
      INTEGER :: i,j
      INTEGER :: counter, cnt

      counter = SIZE(x1)

      ! ... 'x2' is a uniform subsample of 'x1'
      ! ... 'f2' is an average of 'f1' at 'x2' locations 
      ALLOCATE( x2(nsm), f2(nsm) )
      ALLOCATE( dummy(counter) )

      x2 = 0.D0; f2 = 0.D0
      ! ... Uniformly Subsample and smooth the radial function 
      ! ... of the averaged quota. Smoothing is performed by
      ! ... averaging on a window of size 'sstep'
      !
      sstep = REAL( (x1(counter)-x1(1))/(nsm-1),8 )
      x2(1) = x1(1)
      f2(1) = f1(1)
      DO i = 2, nsm-1
        x2(i) = x2(i-1) + sstep
        cnt = 0
        jloop: DO j = 1, counter
          IF (DABS(x1(j)-x2(i)) <= sstep ) THEN
                  f2(i) = f2(i) + f1(j)
                  cnt = cnt + 1
          END IF
        END DO jloop
        f2(i) = f2(i) / cnt
      END DO
      x2(nsm) = x1(counter)
      f2(nsm) = f1(counter)
!
      ! ... Linearly interpolate quotas
      !
      CALL interp_1d_array(x2, f2, x1, f1, dummy)

      DEALLOCATE(x2, f2, dummy)
      RETURN
      END SUBROUTINE filter_1d
!----------------------------------------------------------------------
! ... Filter out high frequency modes by successively subsampling,
! ... averaging and interpolating 
!
      SUBROUTINE filter_2d(x1,y1,f1,nsmx,nsmy)

      IMPLICIT NONE
      REAL*8, INTENT(INOUT), DIMENSION(:) :: x1, y1
      REAL*8, INTENT(INOUT), DIMENSION(:,:) :: f1
      INTEGER, INTENT(IN) :: nsmx, nsmy
      REAL*8, ALLOCATABLE, DIMENSION(:) :: x2, y2
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: f2
      INTEGER, ALLOCATABLE :: dmx(:), dmy(:)
      INTEGER, ALLOCATABLE :: dcx(:,:), dcy(:,:)
      REAL*8 :: sstepx, sstepy
      INTEGER :: i,j,cx,cy,c
      INTEGER :: counterx, countery, cnt

      counterx = SIZE(x1)
      countery = SIZE(y1)

      ! ... 'x2' is a uniform subsample of 'x1'
      ! ... 'f2' is an average of 'f1' at 'x2' locations 
      ! ... dcx and dcy are the extrema of the interval 
      ! ... taken for the average around each mesh point.
      ! ... dmx and dmy are dummies!
      !
      ALLOCATE( x2(nsmx), y2(nsmy), f2(nsmx,nsmy) )
      ALLOCATE( dmx(counterx), dmy(counterx)  )
      ALLOCATE( dcx(nsmx,2), dcy(nsmy,2)  )

      ! ... Initialize the new set coordinates
      ! ... to their limit values
      !
      x2 = 0.D0; y2 = 0.D0; f2 = 0.D0
      dcx(:,1) = counterx; dcy(:,1) = countery 
      dcx(:,2) = 0; dcy(:,2) = 0

      ! ... Uniformly Subsample and smooth the x-y function 
      ! ... of the quota. 'sstep' is the size of the cells
      ! ... used for subsampling
      !
      sstepx = REAL( (x1(counterx)-x1(1))/(nsmx-1) , 8 )
      sstepy = REAL( (y1(countery)-y1(1))/(nsmy-1) , 8 )

      ! ... Coordinates of the subsampled set
      !
      x2(1) = x1(1)
      y2(1) = y1(1)
      !
      DO i = 2, nsmx-1
        x2(i) = x2(i-1) + sstepx
      END DO
      DO j = 2, nsmy-1
        y2(j) = y2(j-1) + sstepy
      END DO
      !
      x2(nsmx) = x1(counterx)
      y2(nsmy) = y1(countery)
!      
      ! ... Initialize the new elevation function
      ! ... on boundaries
      !
      f2(1,1) = f1(1,1)
      f2(nsmx,nsmy) = f1(counterx,countery)
      f2(1,nsmy) = f1(1,countery)
      f2(nsmx,1) = f1(counterx,1)

      ! ... 1D average on boundaries
      !
      DO i = 2, nsmx-1
        cnt = 0
        DO c = 1, counterx
          IF (DABS(x1(c)-x2(i)) <= sstepx ) THEN
                  f2(i,1) = f2(i,1) + f1(c,1)
                  f2(i,nsmy) = f2(i,nsmy) + f1(c,countery)
                  cnt = cnt + 1
                  dcx(i,1) = MIN(c,dcx(i,1))
                  dcx(i,2) = MAX(c,dcx(i,2))
          END IF
        END DO
        f2(i,1) = f2(i,1) / cnt
        f2(i,nsmy) = f2(i,nsmy) / cnt
      END DO
      !
      DO j = 2, nsmy-1
        cnt = 0
        DO c = 1, countery
          IF (DABS(y1(c)-y2(j)) <= sstepy ) THEN
                  f2(1,j) = f2(1,j) + f1(1,c)
                  f2(nsmx,j) = f2(nsmx,j) + f1(counterx,c)
                  cnt = cnt + 1
                  dcy(j,1) = MIN(c,dcy(j,1))
                  dcy(j,2) = MAX(c,dcy(j,2))
          END IF
        END DO
        f2(1,j) = f2(1,j) / cnt
        f2(nsmx,j) = f2(nsmx,j) / cnt
      END DO
      !
      ! ... 2D average 
      !
      DO j = 2, nsmy-1
        DO i = 2, nsmx-1
          cnt = 0
          DO cy = dcy(j,1), dcy(j,2)
            DO cx = dcx(i,1), dcx(i,2)
              IF (DABS(y1(cy)-y2(j)) <= sstepy   .AND. &
                  DABS(x1(cx)-x2(i)) <= sstepx ) THEN
                      f2(i,j) = f2(i,j) + f1(cx,cy)
                      cnt = cnt + 1
              END IF
            END DO
          END DO
          f2(i,j) = f2(i,j) / cnt
        END DO
      END DO
!
      ! ... Linearly interpolate quotas
      !
      CALL interp_2d(x2, y2, f2, x1, y1, f1, dmx, dmy)

      DEALLOCATE(y2, x2, f2, dmx, dmy)
      RETURN
      END SUBROUTINE filter_2d
!----------------------------------------------------------------------
      SUBROUTINE interp_1d_scalar(x1, f1, x2, f2)
      IMPLICIT NONE

      REAL*8, INTENT(IN), DIMENSION(:) :: x1, f1
      REAL*8, INTENT(IN) :: x2
      REAL*8, INTENT(OUT) :: f2
      INTEGER :: i, k, l, n, n1x, t
      REAL*8 :: grad

      n1x = SIZE(x1)
!
! ... locate the grid points near the topographic points
! ... and interpolate linearly the profile  
!
      t = 1
      DO n = 1, n1x
          IF (x1(n) <= x2) t = n
      END DO

      IF (t==1 .OR. t==n1x) THEN
        f2 = f1(t)
      ELSE
        grad = (f1(t+1)-f1(t))/(x1(t+1)-x1(t))
        f2 = f1(t) + (x2-x1(t)) * grad
      END IF
!
      RETURN
      END SUBROUTINE interp_1d_scalar
!----------------------------------------------------------------------
      SUBROUTINE interp_1d_array(x1, f1, x2, f2, t)
      IMPLICIT NONE

      REAL*8, INTENT(IN), DIMENSION(:) :: x1, x2, f1
      REAL*8, INTENT(OUT), DIMENSION(:) :: f2
      INTEGER, INTENT(OUT), DIMENSION(:) :: t
      INTEGER :: i, k, l, n, n1x, n2x
      REAL*8 :: grad

      n1x = SIZE(x1)
      n2x = SIZE(x2)
!
! ... locate the grid points near the topographic points
! ... and interpolate linearly the profile  
!
      l = 1
      DO n = 1, n1x
        DO i = l, n2x

          IF (x1(n) >= x2(i)) THEN

            ! ... t(i) indicates the progressive number of the
            ! ... topographic point laying on the right of a grid center 'i'
            ! ... 'l' counts the grid points
            !
            t(i) = n
 
            IF (n == 1) THEN
              f2(i) = f1(1)
            ELSE
              grad = (f1(n)-f1(n-1))/(x1(n)-x1(n-1))
              f2(i) = f1(n-1) + (x2(i)-x1(n-1)) * grad
            END IF

            l=l+1
          ENDIF

        ENDDO
      ENDDO
!
      RETURN
      END SUBROUTINE interp_1d_array
!----------------------------------------------------------------------
      SUBROUTINE interp_2d(x1, y1, f1, x2, y2, f2, tx, ty)
      IMPLICIT NONE

      REAL*8, INTENT(IN), DIMENSION(:) :: x1, y1, x2, y2
      REAL*8, INTENT(IN), DIMENSION(:,:) :: f1
      REAL*8, INTENT(OUT), DIMENSION(:,:) :: f2
      INTEGER, INTENT(OUT), DIMENSION(:) :: tx, ty

      REAL*8 :: dist1y,dist2y,dist1x,dist2x,alpha,beta
      INTEGER :: i,j,k,h,l,ii,jj
      INTEGER :: n1x, n2x, n1y, n2y
      INTEGER :: tp1, tp2

      n1x = SIZE(x1)
      n1y = SIZE(y1)
      n2x = SIZE(x2)
      n2y = SIZE(y2)
        
!C============================================
!C===    trova le posizioni dei nodi      ====
!C===    della nuova griglia rispetto     ====
!C===    alla griglia iniziale            ====
!C============================================

      l=1
      DO i = 1, n1x
        DO ii = l, n2x
      
!============================================
!===    cerca i nodi della griglia che    ===
!===    che stanno a sx. di xtop(i)       ===
!============================================
    
          IF (x1(i) >= x2(ii)) THEN
            tx(ii)=i
            l=l+1
          ENDIF
        ENDDO
      ENDDO

      l=1
      DO j = 1, n1y
        DO jj = l, n2y

!C============================================
!C===    cerca i nodi della griglia che    ===
!C===    che stanno a sotto  ytop(i)       ===
!C============================================

          IF (y1(j) >= y2(jj)) THEN
            ty(jj)=j
            l=l+1
          ENDIF
        ENDDO
      ENDDO


! il nodo della nuova griglia di indici (i,j) sara' allora
! contenuto nel rettangolo con i vertici con indici:
! P1=tx(i),ty(j)
! P2=tx(i-1),ty(j)
! P3=tx(i-1),ty(j-1)
! P4=tx(i),ty(j-1)

! sulla nuova griglia interpoliamo le quote di input ztop per 
! ottenere la quota coorZ nel punto di indici (i,j)
 

! interpolazione bilineare sui nodi interni (1<i<nodiGRIDx, 1<j<nodigGRIDy)
! utilizzando le quote nei punti P1,..,P4 definiti sopra

        DO i=1,n2x

           dist1x = x2(i) - x1(tx(i)-1)
           dist2x = x1(tx(i)) - x2(i)
           alpha  = dist1x/(dist1x+dist2x)

           DO j=1,n2y
           
              dist1y = y2(j) - y1(ty(j)-1)
              dist2y = y1(ty(j)) - y2(j)
              beta   = dist1y/(dist1y+dist2y)

              tp1    = alpha * f1(tx(i),ty(j))   + &
                       (1.D0 - alpha) * f1(tx(i)-1,ty(j))
              tp2    = alpha * f1(tx(i),ty(j)-1) + &
                       (1.D0 - alpha) * f1(tx(i)-1,ty(j)-1)
              f2(i,j) = beta * tp1 + (1.D0 - beta) * tp2
 
           ENDDO
        ENDDO
!      
      RETURN
      END SUBROUTINE interp_2d
!----------------------------------------------------------------------
      END MODULE array_filters
!----------------------------------------------------------------------
