!----------------------------------------------------------------------
      MODULE grid
!----------------------------------------------------------------------
        IMPLICIT NONE
        SAVE
!
        INTEGER, PARAMETER :: bb_ = 1
        INTEGER, PARAMETER :: bl_ = 2
        INTEGER, PARAMETER :: b_  = 3
        INTEGER, PARAMETER :: br_ = 4
        INTEGER, PARAMETER :: ll_ = 5
        INTEGER, PARAMETER :: l_  = 6
        INTEGER, PARAMETER :: r_  = 7
        INTEGER, PARAMETER :: rr_ = 8
        INTEGER, PARAMETER :: tl_ = 9
        INTEGER, PARAMETER :: t_  = 10
        INTEGER, PARAMETER :: tr_ = 11
        INTEGER, PARAMETER :: tt_ = 12
!
        TYPE rcv_map_type
          INTEGER :: nrcv               ! How Many Cells I Have to Receive
          INTEGER, POINTER :: ijrcv(:)  ! Which Cells I Have to Receive (Global Index) 
          INTEGER, POINTER :: ijl(:)    ! Where I Have to Put Them (Local Index)
        END TYPE 

        TYPE snd_map_type
          INTEGER :: nsnd               ! How Many Cells I Have to Send
          INTEGER, POINTER :: ijsnd(:)  ! Which Cells I Have to Send (Global Index)
          INTEGER, POINTER :: ijl(:)    ! Where I Have to Pick Them (Local Index)
        END TYPE 

        TYPE cells_map_type
          INTEGER :: left(2)           ! The first cell of the block (x,y)
          INTEGER :: right(2)          ! The last cell of the block (x,y)
        END TYPE

        INTEGER :: nbx,nby             ! number of blocks in x and y directions

        TYPE (rcv_map_type), ALLOCATABLE :: rcv_map(:)
        TYPE (snd_map_type), ALLOCATABLE :: snd_map(:)
        TYPE (cells_map_type), ALLOCATABLE :: block_map(:)
        INTEGER, ALLOCATABLE :: lay_map(:,:)
        INTEGER, ALLOCATABLE :: nij_lay(:), nbl_lay(:), nij1_lay(:)
        INTEGER, ALLOCATABLE :: nij1(:), n1(:), diff(:)
        INTEGER, ALLOCATABLE :: nij(:)
!
        REAL*8, DIMENSION(:), ALLOCATABLE :: r, rb, dr
        REAL*8, DIMENSION(:), ALLOCATABLE :: dx, dy, dz
        REAL*8, DIMENSION(:), ALLOCATABLE :: xb, yb, zb
        REAL*8, DIMENSION(:), ALLOCATABLE :: inr, inrb, inzb, indr, indz
!
        INTEGER :: itc, part
        INTEGER :: nij_l, nije_l, nijx_l
        INTEGER, ALLOCATABLE :: myij(:,:,:)
        INTEGER, ALLOCATABLE :: myinds(:,:)

!
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: iob
      INTEGER, DIMENSION(:), ALLOCATABLE :: nso
      INTEGER, DIMENSION(:), ALLOCATABLE :: fl
!
      INTEGER, DIMENSION(:), ALLOCATABLE :: fl_l
!
        INTEGER :: countfl(5)

        INTERFACE data_exchange
          MODULE PROCEDURE data_exchange_i, data_exchange_r, data_exchange_rm
        END INTERFACE

!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE bounds_grid
      USE dimensions
      IMPLICIT NONE
!
       ALLOCATE(fl(nr*nz))
       ALLOCATE(r(nr), rb(nr), zb(nz), dr(nr), dz(nz))
       ALLOCATE(inr(nr), inrb(nr), inzb(nz), indr(nr), indz(nz))
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE bounds_blocks
      USE dimensions
      IMPLICIT NONE
!
       ALLOCATE(iob(4,no))
       ALLOCATE(nso(no))
      RETURN
      END SUBROUTINE
!----------------------------------------------------------
        SUBROUTINE grid_setup(zzero)
        USE dimensions
          REAL*8, INTENT(IN) :: zzero
          REAL*8, PARAMETER :: VERYBIG = 1.0d+10
          INTEGER :: i, j
          IF(itc.EQ.0) THEN
            DO i=1,nr
              r(i)=1.D0
              rb(i)=1.D0
            END DO
          ELSE
            rb(1)=0.D0
            r(1)=-0.5D0*dr(1)+rb(1)
            DO i=2,nr
              rb(i)=(rb(i-1)+dr(i))
              r(i)=(rb(i)-0.5D0*dr(i))
            END DO
          END IF
!
          DO i=1,nr
            inr(i)=1.D0/r(i)
          END DO  


          IF (rb(1) .EQ. 0.0D0) THEN
            inrb(1) = VERYBIG
          ELSE
            inrb(1) = 1.0D0/rb(1)
          END IF
          IF (itc.EQ.0) inrb(1)=1.D0
          DO i=2,nr
            inrb(i)=1.D0/rb(i)
          END DO  

          DO i=1,nr
            indr(i)=1.D0/dr(i)
          END DO  
!
          zb(1) = zzero
          DO j=1,(nz-1)
            zb(j+1)=zb(j)+dz(j+1)
          END DO
!
          IF( zzero .EQ. 0.0d0 ) THEN
            inzb(1) = VERYBIG
          ELSE
            inzb(1)=1.D0/zb(1)
          END IF
          DO j=2,(nz-1)
            inzb(j)=1.D0/zb(j)
          END DO  

          DO j=1,(nz-1)
            indz(j)=1.D0/dz(j)
          END DO  
!
        END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE flic
!
      USE dimensions
      IMPLICIT NONE
!
      INTEGER :: i2, i1, j2, j1
      INTEGER :: i, j, ij, n, ijl
!
      DO j=1,nz
        DO i=1,nr
          ij=i+(j-1)*nr
          fl(ij)=1
        END DO
      END DO

!
! ... upper-right corner
!
      IF (fl((nr-1)+(nz-1)*nr) == 4 .AND. fl(nr+(nz-2)*nr) == 4) THEN
          fl(nr+(nz-1)*nr) = 4
      END IF
!
      IF(no.LE.0) RETURN
!
      DO n=1,no
        i1=iob(1,n)
        i2=iob(2,n)
        j1=iob(3,n)
        j2=iob(4,n)
        DO j=j1,j2
          DO i=i1,i2
            ij=i+(j-1)*nr
            IF(nso(n).EQ.2) fl(ij)=2
            IF(nso(n).EQ.3) fl(ij)=3
            IF(nso(n).EQ.4) fl(ij)=4
            IF(nso(n).EQ.5) fl(ij)=5
            IF(nso(n).EQ.6) fl(ij)=6
            IF(nso(n).EQ.7) fl(ij)=7
          END DO
        END DO
      END DO
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE partition
! ... a routine that subdivides a regular grid between (nproc) 
! ... processors, balancing the computational load 
! ... Gives processor maps.
!
      USE dimensions
      USE parallel, ONLY: nproc, mpime, root
      IMPLICIT NONE
      SAVE
!
      INTEGER :: i, j, ij, n, m, q, rr, k
      INTEGER :: nfl, ipe, dim
      INTEGER :: layer
      INTEGER, ALLOCATABLE :: red(:), green(:), blue(:)
      INTEGER :: localdim
!
      ALLOCATE(red(nr*nz),green(nr*nz),blue(nr*nz))
      IF (ALLOCATED(nij)) DEALLOCATE(nij)
      IF (ALLOCATED(nij1)) DEALLOCATE(nij1)
      IF (ALLOCATED(n1)) DEALLOCATE(n1)
      IF (ALLOCATED(diff)) DEALLOCATE(diff)
      ALLOCATE(nij(0:nproc-1))
      ALLOCATE(nij1(0:nproc-1))
      ALLOCATE(n1(0:nproc-1))
      ALLOCATE(diff(0:nproc-1))
      nij = 0
      nij1 = 0
      n1 = 0
      diff = 0
! ... here count the flags
      DO i = 1, SIZE(countfl)
        countfl(i) = COUNT((fl.EQ.i))
      END DO
      WRITE(7,'(a13,5i7)') ' # countfl ', countfl
!
! ... domain decomposition (build maps)
!     1=layers
!     2=blocks

      part = 2
      IF (part .EQ. 1 .OR. nproc .LE. 2) THEN
        CALL layers
      ELSE 
        CALL blocks
      END IF
      DO i= 0, nproc-1
      END DO
!
! ... n1 is the balanced number of cells to each processor:
! 
      DO ipe = 0, nproc - 1
        n1(ipe) = localdim(countfl(1), nproc, ipe)
        diff(ipe) = nij1(ipe) - n1(ipe)
      END DO
!
        WRITE(7,*) '---  partition  -----'
        WRITE(7,*) '  '
      DO ipe = 0, nproc - 1
        WRITE(7,*) ' # nij( ',ipe, ' )', nij(ipe)
        WRITE(7,*) ' # nij1( ',ipe, ' )', nij1(ipe),' -', n1(ipe),' =', diff(ipe)
      END DO
        WRITE(7,*) '  '
!
        DO ij=1,nr*nz
            red(ij) = (cell_owner(ij)+1)
            green(ij) = nproc - (cell_owner(ij)+1)
            blue(ij) = (nproc - ABS(2*(cell_owner(ij)-nproc/2))) 
          IF ((fl(ij).EQ.3 .AND. (fl(ij+1).EQ.1 .OR. fl(ij-1).EQ.1))   &
               .OR. fl(ij).EQ.5 .OR. fl(ij).EQ.2 .OR. fl(ij).EQ.4)    THEN
            red(ij) = 0
            green(ij) = 0
            blue(ij) = 0
          END IF
        END DO
!
        IF (mpime .EQ. root) THEN
          OPEN(UNIT=9,FILE='procs_map.ppm', FORM='FORMATTED', STATUS='UNKNOWN')
           WRITE(9,'(A2)') 'P3'
           WRITE(9,'(I3,A1,I3)') nr,' ',nz
           WRITE(9,'(I3)') nproc
           DO j = nz-1, 0, -1
            DO i=1,nr
              WRITE(9,100) red(i+nr*j), green(i+nr*j), blue(i+nr*j)
            END DO
           END DO
 100      FORMAT(3(I4))
          CLOSE(9)
        END IF
!
      DEALLOCATE(red,green,blue)
      RETURN
      END SUBROUTINE 
! ------------------------------------------------------------------------
      SUBROUTINE layers
! 
! ... partition N.1 (layers)
! ... decomposes the domain into N horizontal layers
!
      USE dimensions
      USE parallel, ONLY: nproc
      IMPLICIT NONE
      SAVE
!
      INTEGER :: i, j, ij
      INTEGER :: icnt, icnt_ipe, ipe
      INTEGER :: localdim
!
      IF(ALLOCATED(lay_map)) DEALLOCATE(lay_map)
      ALLOCATE(lay_map(0:nproc-1,2))
!
      nij1 = 0
      nij  = 0
!
      DO ipe = 0, nproc - 1
       nij1(ipe) = localdim(countfl(1), nproc, ipe)
      END DO
!
      lay_map(0,1) = 1
      lay_map(nproc-1,2) = nr*nz
      ipe = 0
      DO j = 1, nz
       DO i = 1, nr
          ij = i + (j-1)*nr

          IF ( fl(ij) .EQ. 1 ) THEN
            icnt_ipe = icnt_ipe + 1
            icnt     = icnt     + 1
          END IF

          nij(ipe) = nij(ipe) + 1

          IF ( icnt_ipe .EQ. nij1(ipe) .AND. i .NE. (nr-1) ) THEN
            icnt_ipe = 0
            IF( icnt .LT. countfl(1) ) THEN
              lay_map(ipe,2) = ij
              ipe = ipe + 1
              lay_map(ipe,1) = ij+1
            END IF
          END IF
       END DO 
      END DO
!
      RETURN
      END SUBROUTINE
! ------------------------------------------------------------------------
      SUBROUTINE blocks
!
! ... partition N.2 (blocks)
!
      USE dimensions
      USE parallel, ONLY: nproc
      IMPLICIT NONE
      SAVE
!
      INTEGER :: i, j, ij
      INTEGER :: i1, i2, j1, j2, ij1, ij2
      INTEGER :: nbl
      INTEGER :: icnt, ipe
      INTEGER :: bl1
      INTEGER :: size_x, rest, skipl, rrest
      INTEGER :: layer, icnt_layer, here
      INTEGER :: localdim
!
      REAL*8 side, area
      REAL*8 fact
!
        area  = DBLE(countfl(1)/nproc)
        side  = DSQRT(area)
        nbx    = NINT(nr/side)
! ... compute the number of layers nby 
        nby    = NINT(countfl(1)/nbx/side**2)     
        IF (nbx .EQ. 0) nbx = 1
        IF (nby .EQ. 0) nby = 1
        DO WHILE (nbx*nby .LT. nproc)
         IF (INT(nr/nbx) .GT. INT(nz/nby)) THEN
           nbx = nbx + 1
         ELSE
           nby = nby + 1
         END IF
        END DO
        size_x = INT(nr/nbx)
        rest = nbx*nby - nproc
!
        IF(ALLOCATED(lay_map)) DEALLOCATE(lay_map)
        ALLOCATE(lay_map(1:nby,2))
        IF(ALLOCATED(block_map)) DEALLOCATE(block_map)
        ALLOCATE(block_map(0:nproc-1))
        lay_map(:,:) = 0
        DO i = 0, nproc-1
          block_map(i)%left(:) = 0
          block_map(i)%right(:) = 0
        END DO
!
        IF (ALLOCATED(nij_lay)) DEALLOCATE(nij_lay)
        IF (ALLOCATED(nij1_lay)) DEALLOCATE(nij1_lay)
        IF (ALLOCATED(nbl_lay)) DEALLOCATE(nbl_lay)
        ALLOCATE(nij_lay(1:nby), nij1_lay(1:nby), nbl_lay(1:nby))
!
        nij_lay = 0
        nij1_lay = 0
        nbl_lay = 0
!
! ... distribute cells among layers proportionally
! ... to the number of blocks (processors) contained
! ... (compute nij1_lay(layer))
!
        IF (rest .EQ. 0) THEN
          DO layer = 1, nby
           nij1_lay(layer) = localdim(countfl(1),nby,layer)
           nbl_lay(layer) = nbx
          END DO
        ELSE
          skipl = INT(rest/nbx) 
          rrest = MOD(rest, nbx)
          IF (rest .LE. nbx) THEN
            nbl_lay(1) = nbx - rrest
            DO layer = 2, nby
              nbl_lay(layer) = nbx
            END DO
          ELSE
            DO layer = 1, skipl
              nbl_lay(layer) = 1
            END DO
              nbl_lay(skipl+1) = nbx - rrest 
            DO layer = skipl + 2, nby
              nbl_lay(layer) = nbx
            END DO
          END IF
          DO layer = 1, nby
            fact = DBLE(countfl(1)/nproc*nbl_lay(layer))
            nij1_lay(layer) = NINT(fact)
          END DO
        END IF
!
! ... build the layer maps
!
        layer = 1
        ipe = 0
        lay_map(1,1) = 1
        lay_map(nby,2) = nr*nz
        DO j = 1, nz
        DO i = 1, nr
          ij = i + (j-1)*nr
!
          IF ( fl(ij) .EQ. 1 ) THEN
            icnt_layer = icnt_layer + 1
            icnt     = icnt     + 1
          END IF

          nij_lay(layer) = nij_lay(layer) + 1

         IF ( icnt_layer .EQ. nij1_lay(layer) ) THEN
            icnt_layer = 0
            IF(layer .LT. nby) THEN
              lay_map(layer,2) = ij
              layer = layer + 1
              lay_map(layer,1) = ij+1
            END IF
          END IF
        END DO
        END DO
!
! ... cut steps
!
      ipe = -1
      DO layer = 1, nby
        ij2 = lay_map(layer,2)
        j2  = ( ij2 - 1 ) / nr + 1
        i2  = MOD( ( ij2 - 1 ), nr) + 1
        IF (i2 .LT. nr/2) j2 = j2-1
        IF (layer .LT. nby) THEN
          lay_map(layer,2) = j2 * nr
        ELSE
          lay_map(layer,2) = nz * nr
        END IF
        IF (layer .GT. 1) THEN 
          lay_map(layer,1) = lay_map(layer-1,2) + 1
        ELSE
          lay_map(layer,1) = 1
        END IF
        ij1 = lay_map(layer,1)
        ij2 = lay_map(layer,2)
        j1  = ( ij1 - 1 ) / nr + 1
        i1  = MOD( ( ij1 - 1 ), nr) + 1
        j2  = ( ij2 - 1 ) / nr + 1
        i2  = MOD( ( ij2 - 1 ), nr) + 1
        IF (i1.NE.1 .OR. i2.NE.nr) WRITE(8,*)'error in layer',layer
!
! ...   updates the number of cells with fl=1 into layers
!
        nij_lay(layer)  = 0
        nij1_lay(layer) = 0
        DO ij = ij1, ij2
          nij_lay(layer) = nij_lay(layer) + 1
          IF (fl(ij) .EQ. 1) nij1_lay(layer) = nij1_lay(layer) + 1 
        END DO
!
! ...   build block maps        
!
        bl1 = INT(nij1_lay(layer)/nbl_lay(layer))
        i = i1
        DO nbl = 1, nbl_lay(layer) 
          ipe = ipe + 1
          block_map(ipe)%left(1) = i
          block_map(ipe)%left(2) = j1
          DO WHILE (nij1(ipe) .LT. bl1)
            DO j = j1, j2
              ij = i + (j-1)*nr
              IF (fl(ij) .EQ. 1) nij1(ipe) = nij1(ipe) + 1
            END DO  
            i = i+1
          END DO
          block_map(ipe)%right(1) = i-1
          IF (nbl .EQ. nbl_lay(layer)) block_map(ipe)%right(1) = nr 
          block_map(ipe)%right(2) = j2
!
! ...     updates the number of cells with fl=1 into blocks
!
          nij(ipe)  = 0
          nij1(ipe) = 0
          DO j = block_map(ipe)%left(2), block_map(ipe)%right(2)
          DO i = block_map(ipe)%left(1), block_map(ipe)%right(1)
            ij = i + (j-1)*nr
            nij(ipe) = nij(ipe) + 1
            IF (fl(ij) .EQ. 1) nij1(ipe) = nij1(ipe) + 1
          END DO
          END DO
          WRITE(7,*)'block_map(',ipe,'):',block_map(ipe)%left(:),block_map(ipe)%right(:)
        END DO
      END DO

      RETURN
      END SUBROUTINE
! ---------------------------------------------------------------------
      SUBROUTINE ghost
!
! ... Identifies and allocates ghost cells
! 
      USE dimensions
      USE parallel, ONLY: nproc, mpime
      USE basic_types, ONLY: imatrix
      IMPLICIT NONE
!
      INTEGER :: icnt, ipe, icnt_ipe
      INTEGER :: i, j, ij, n, ijl, k
      INTEGER :: nrcv(0:nproc-1), nrcvx, jpe
      INTEGER :: nset(0:nproc-1)
      INTEGER, ALLOCATABLE :: ijrcv(:,:)
      INTEGER, ALLOCATABLE :: ijsnd(:)
      TYPE (imatrix) :: rcv_cell_set(0:nproc-1)
      INTEGER :: localdim
      INTEGER :: layer,j2,j1,i2,i1
      INTEGER :: me
!
      IF(ALLOCATED(rcv_map)) DEALLOCATE(rcv_map)
      IF(ALLOCATED(snd_map)) DEALLOCATE(snd_map)
      ALLOCATE(rcv_map( 0:(nproc-1) ) )
      ALLOCATE(snd_map( 0:(nproc-1) ) )
      

      rcv_map(:)%nrcv = 0
      snd_map(:)%nsnd = 0
      icnt     = 0
      icnt_ipe = 0
      ipe      = 0
!
! ... count the neighbouring cells on other processors
! ... nset is the number of the neighbouring cells on other processor,
! ... summed over each local cell.
!
      nset   = 0
      nije_l = 0
      IF (part .EQ. 1 .OR. nproc .LE. 2) THEN
        DO ij = lay_map(mpime,1), lay_map(mpime,2)
          IF ( fl(ij) .EQ. 1 ) THEN
            nije_l = nije_l + cell_neighbours(ij, mpime, nset)
          END IF
        END DO
      ELSE
        i1 = block_map(mpime)%left(1)
        i2 = block_map(mpime)%right(1)
        j1 = block_map(mpime)%left(2)
        j2 = block_map(mpime)%right(2)
        DO j = j1, j2
         DO i = i1, i2
          ij = i + nr*(j-1)
          IF (fl(ij).EQ.1) THEN
            nije_l = nije_l + cell_neighbours(ij, mpime, nset)
          END IF
         END DO
        END DO
      END IF
!
! ... ALLOCATE memory for the set of neighbouring cell on other proc.
      DO ipe = 0, nproc - 1
        IF(nset(ipe) .GT. 0 ) THEN
          ALLOCATE( rcv_cell_set(ipe)%i(4, nset(ipe)) )
        ELSE
          NULLIFY( rcv_cell_set(ipe)%i )
        END IF
      END DO

! ... number of local cells
      nij_l = nij(mpime)
!
! ... allocate the indexes matrix
      ALLOCATE(myij(-2:2,-2:2,nij_l))
      ALLOCATE(myinds(12,nij_l))
      myij    = 0
      myinds  = 0

! ... now fill in the receiving map and the indexes matrix
      nset  = 0
      icnt  = 0
      nrcv  = 0
      IF (part .EQ. 1 .OR. nproc .LE. 2) THEN
        DO ij = lay_map(mpime,1), lay_map(mpime,2)
          IF ( fl(ij) .EQ. 1 ) THEN
          icnt = icnt + cell_neighbours(ij, mpime, nset, rcv_cell_set, myij)
          END IF
        END DO
      ELSE
        i1 = block_map(mpime)%left(1)
        i2 = block_map(mpime)%right(1)
        j1 = block_map(mpime)%left(2)
        j2 = block_map(mpime)%right(2)
       DO j = j1, j2
        DO i = i1, i2
          ij = i + nr*(j-1)
          IF (fl(ij).EQ.1) THEN
          icnt = icnt + cell_neighbours(ij, mpime, nset, rcv_cell_set, myij)
          END IF
        END DO
       END DO
      END IF
!
      DO ijl = 1, nij_l
! ...     store the global cell index ij, of the the local cell ijl
          myij(0, 0, ijl) = cell_l2g(ijl, mpime)
      END DO

      CALL cell_set_sort(rcv_cell_set, nset, nrcv)

! ... nrcv is the number of elements to be received from other proc
! ... nrcv is equal to nset after the copies of the cells has been
! ... eliminated

      DO ipe = 0, nproc - 1
        WRITE(7,300) nset(ipe), ipe
        IF(nset(ipe) .GT. 0 ) THEN
          WRITE(7,310) rcv_cell_set(ipe)%i(1,:)
        END IF
 300    FORMAT(' # neighbours set SIZE ',i5,' from ',i3)
 310    FORMAT(10i8)
      END DO

! ... number cells required to update the physical quantities:
! ... sum of the local and their neighbouring cells
      nije_l = SUM( nrcv ) 
      nijx_l = nij_l + nije_l
      WRITE(7,* ) ' # nije_l ', nije_l
      WRITE(7,* ) ' # nijx_l ', nijx_l

! ... prepare the receive map 
      DO ipe = 0, nproc - 1
        rcv_map(ipe)%nrcv = nrcv(ipe)
        NULLIFY( rcv_map(ipe)%ijrcv )
        NULLIFY( rcv_map(ipe)%ijl )
      END DO

! ... prepare the send map 
      DO ipe = 0, nproc - 1
        NULLIFY( snd_map(ipe)%ijsnd ) 
        NULLIFY( snd_map(ipe)%ijl ) 
      END DO

! ... ALLOCATE and set the receiving map
      CALL set_rcv_map(rcv_map, myij, rcv_cell_set, nset, nrcv)

! ... print out basic information on the map
      DO ipe = 0, nproc - 1
        WRITE(7,*) ' # nrcv ', nrcv(ipe), ' from ', ipe
      END DO
! ... using the receive maps fill in the sending maps
      DO ipe = 0, nproc - 1
        IF( ASSOCIATED (snd_map(ipe)%ijsnd ) ) THEN
          DEALLOCATE( snd_map(ipe)%ijsnd)
        END IF
        IF( ASSOCIATED (snd_map(ipe)%ijl ) ) THEN
          DEALLOCATE( snd_map(ipe)%ijl)
        END IF
        NULLIFY( snd_map(ipe)%ijsnd ) 
        NULLIFY( snd_map(ipe)%ijl ) 

! ...   proc ipe send to other processor the number of elements 
! ...   he should receive from each of the other proc. (rcv_map(:)%nrcv)
        IF( ipe .EQ. mpime ) THEN
          nrcv(:) = rcv_map(:)%nrcv
          nrcvx = MAXVAL(nrcv)
        END IF
! ...   each proc store the number of element he should send to ipe
        CALL scatter_integer(nrcv, snd_map(ipe)%nsnd, 1, ipe)
        CALL bcast_integer(nrcvx, 1, ipe)

! ...   proc ipe send to other processor the array with the indexes of 
! ...   the elements he should receive from each of the other proc. 
! ...   ( rcv_map(jpe)%ijrcv )
          ALLOCATE(ijrcv(nrcvx,0:nproc-1))
          ALLOCATE(ijsnd(nrcvx))
        IF( ipe .EQ. mpime) THEN
          ijrcv = 0
          DO jpe = 0, nproc - 1
            IF( rcv_map(jpe)%nrcv .GT. 0 ) THEN
              ijrcv(1:rcv_map(jpe)%nrcv, jpe ) = rcv_map(jpe)%ijrcv(1:rcv_map(jpe)%nrcv) 
            END IF
          END DO 
        END IF 
! ...   each proc store the index of the elements he should send to ipe
        CALL scatter_integer(ijrcv, ijsnd, nrcvx, ipe)
        IF( ipe .NE. mpime) THEN
          ALLOCATE(snd_map(ipe)%ijsnd(snd_map(ipe)%nsnd))
          ALLOCATE(snd_map(ipe)%ijl(snd_map(ipe)%nsnd))
          snd_map(ipe)%ijsnd(1:snd_map(ipe)%nsnd) = ijsnd(1:snd_map(ipe)%nsnd)
          DO ij = 1, snd_map(ipe)%nsnd
            snd_map(ipe)%ijl( ij ) = cell_g2l( ijsnd( ij ), mpime )
          END DO
        END IF 
!
          DEALLOCATE(ijrcv)
          DEALLOCATE(ijsnd)

      END DO

      DO ipe = 0, nproc - 1
        WRITE(7,100) rcv_map(ipe)%nrcv, ipe
        IF(rcv_map(ipe)%nrcv .GT. 0 ) THEN
          WRITE(7,110) rcv_map(ipe)%ijrcv(:)
          WRITE(7,*) ' ---- '
          WRITE(7,110) rcv_map(ipe)%ijl(:)
        END IF
 100    FORMAT(' # receiving ',i5,' cells from ',i3)
 110    FORMAT(10i8)
      END DO

      DO ipe = 0, nproc - 1
        WRITE(7,200) snd_map(ipe)%nsnd, ipe
        IF (snd_map(ipe)%nsnd .GT. 0 ) THEN
          WRITE(7,210) snd_map(ipe)%ijsnd(:)
          WRITE(7,*) ' ---- '
          WRITE(7,210) snd_map(ipe)%ijl(:)
        END IF
 200    FORMAT(' # sending ',i5,' cells to ',i3)
 210    FORMAT(10i8)
      END DO


      CALL test_comm
      
! ... fill in the array myinds using myij
      ALLOCATE(fl_l(nijx_l))
      DO ijl = 1, nij_l
        ij = myij( 0, 0, ijl)
        fl_l(ijl) = fl(ij)
      END DO
      CALL data_exchange(fl_l)
      CALL set_myinds(myinds, myij)
!
! ... Test the indexing (following sweep direction)
!      DO ijl = 1, nij_l
!        WRITE(7,*) ijl
!        WRITE(7,*) (myinds(k,ijl), k=1,12)
!        WRITE(7,*) ((myij(i,j,ijl),i=-2,2),j=-2,2)
!      END DO
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      INTEGER FUNCTION cell_owner(ij)
! ... gives the owner of a cells from its global index (ij)
        USE dimensions
        USE parallel, ONLY: nproc
        INTEGER, INTENT(IN) :: ij
        INTEGER :: ipe, layer, cell_layer
        INTEGER :: i,j
        j = ( ij - 1 ) / nr + 1
        i = MOD( ( ij - 1 ), nr) + 1
          
        DO ipe = 0, nproc - 1
          IF (part .EQ. 1 .OR. nproc .LE. 2) THEN
            IF( ij .GE. lay_map(ipe,1) .AND. ij .LE. lay_map(ipe,2) ) THEN
              cell_owner = ipe
            END IF
          ELSE
            IF (j .GE. block_map(ipe)%left(2) .AND. j .LE. block_map(ipe)%right(2)) THEN
              IF (i .GE. block_map(ipe)%left(1) .AND.  &
                  i .LE. block_map(ipe)%right(1)) cell_owner = ipe
            END IF
          END IF
        END DO
      RETURN
      END FUNCTION     
!----------------------------------------------------------------------
      INTEGER FUNCTION cell_l2g(ijl, mpime)
      USE dimensions
      USE parallel, ONLY: nproc
      INTEGER, INTENT(IN) :: ijl, mpime
      INTEGER :: i,j,i1,i2,j1
!
      IF (part .EQ. 1 .OR. nproc .LE. 2) THEN
        cell_l2g = ijl + lay_map(mpime,1) - 1
      ELSE
        j1 = block_map(mpime)%left(2)
        i1 = block_map(mpime)%left(1)
        i2 = block_map(mpime)%right(1)
        i = MOD( ( ijl - 1 ), (i2-i1+1)) + i1
        j = ( ijl - 1 ) / (i2-i1+1) + j1
        cell_l2g = i + (j-1)*nr
      END IF
!
      RETURN
      END FUNCTION     
!----------------------------------------------------------------------
      INTEGER FUNCTION cell_g2l(ij, mpime)
      USE dimensions
      USE parallel, ONLY: nproc
      INTEGER, INTENT(IN) :: ij, mpime
      INTEGER :: i,j,i1,i2,j1
!
      IF (part .EQ. 1 .OR. nproc .LE. 2) THEN
        cell_g2l = ij - lay_map(mpime,1) + 1
      ELSE
        j = ( ij - 1 ) / nr + 1
        i = MOD( ( ij - 1 ), nr) + 1
        j1 = block_map(mpime)%left(2)
        i1 = block_map(mpime)%left(1)
        i2 = block_map(mpime)%right(1)
        cell_g2l = (i-i1+1) + (j-j1)*(i2-i1+1)
      END IF
!
      RETURN
      END FUNCTION     
!----------------------------------------------------------------------
      INTEGER FUNCTION cell_neighbours(ij, mpime, nset, rcv_cell_set, myij)
        USE dimensions
        USE basic_types, ONLY: imatrix
        INTEGER, INTENT(IN) :: ij, mpime
        INTEGER, INTENT(INOUT) :: nset(0:)
        TYPE (imatrix), OPTIONAL :: rcv_cell_set(0:)
        INTEGER, OPTIONAL :: myij(-2:,-2:,:)
        INTEGER :: ippj, ijpp, ijmm, immj
        INTEGER :: icnt, ipe, i, j, ije, ijl, ijel
        INTEGER :: im, jm
        LOGICAL :: fill

        IF ( PRESENT(rcv_cell_set) .AND. .NOT. present(myij) ) THEN
          WRITE(8,*) ' error in neighbour - 1 '
          STOP 
        END IF
        IF( .NOT. PRESENT(rcv_cell_set) .AND. present(myij) ) THEN
          WRITE(8,*) ' error in neighbour - 2 '
          STOP 
        END IF

        fill = PRESENT(rcv_cell_set)
        icnt = 0
        IF( fill ) THEN
          nij_l = SIZE(myij, 3)
          ijl   = cell_g2l(ij, mpime)
        END IF
        
! ...   loop over the first neighbouring cells
        DO im = -1, 1
          DO jm = -1, 1
            IF( (im .NE. 0) .OR. (jm .NE. 0) ) THEN
              ije = ij + im + jm*nr
              ipe =  cell_owner(ije) 
              IF( ipe .NE. mpime) THEN
! ...           the cell ije is not local, count and register its position
                nset(ipe) = nset(ipe) + 1
                icnt = icnt + 1
                IF( fill ) THEN
                  rcv_cell_set(ipe)%i(1,nset(ipe)) = ije
                  rcv_cell_set(ipe)%i(2,nset(ipe)) = ijl
                  rcv_cell_set(ipe)%i(3,nset(ipe)) = im
                  rcv_cell_set(ipe)%i(4,nset(ipe)) = jm
                END IF
              ELSE IF( fill ) THEN
! ...           the cell ije is local, set the mapping with cell ijl
                ijel = cell_g2l(ije, mpime)
                myij(im,jm,ijl) = ijel
              END IF
            ELSE IF( fill ) THEN   
! ...         store the global cell index ij, of the the local cell ijl
              myij(0, 0, ijl) = ij
            END IF
          END DO
        END DO

        j = ( ij - 1 ) / nr + 1
        i = MOD( ( ij - 1 ), nr) + 1

        immj=ij-2
        IF(i .EQ. (2)) immj = ij-1
        ipe =  cell_owner(immj) 
        IF( ipe .NE. mpime) THEN
          nset(ipe) = nset(ipe) + 1
          icnt = icnt + 1
          IF( fill ) THEN
            rcv_cell_set(ipe)%i(1,nset(ipe)) = immj
            rcv_cell_set(ipe)%i(2,nset(ipe)) = ijl
            rcv_cell_set(ipe)%i(3,nset(ipe)) = -2
            rcv_cell_set(ipe)%i(4,nset(ipe)) = 0
          END IF
        ELSE IF( fill ) THEN
          ijel = cell_g2l(immj, mpime)
          myij(-2,0,ijl) = ijel
        END IF

        ippj=ij+2
        IF(i .EQ. (nr-1)) ippj = ij+1
        ipe =  cell_owner(ippj) 
        IF( ipe .NE. mpime) THEN
          nset(ipe) = nset(ipe) + 1
          icnt = icnt + 1
          IF( fill ) THEN
            rcv_cell_set(ipe)%i(1,nset(ipe)) = ippj
            rcv_cell_set(ipe)%i(2,nset(ipe)) = ijl
            rcv_cell_set(ipe)%i(3,nset(ipe)) = 2
            rcv_cell_set(ipe)%i(4,nset(ipe)) = 0
          END IF
        ELSE IF( fill ) THEN
          ijel = cell_g2l(ippj, mpime)
          myij(2,0,ijl) = ijel
        END IF

        ijmm=ij-nr-nr
        IF(j .EQ. (2)) ijmm = ij-nr
        ipe =  cell_owner(ijmm) 
        IF( ipe .NE. mpime) THEN
          nset(ipe) = nset(ipe) + 1
          icnt = icnt + 1
          IF( fill ) THEN
            rcv_cell_set(ipe)%i(1,nset(ipe)) = ijmm
            rcv_cell_set(ipe)%i(2,nset(ipe)) = ijl
            rcv_cell_set(ipe)%i(3,nset(ipe)) = 0
            rcv_cell_set(ipe)%i(4,nset(ipe)) = -2
          END IF
        ELSE IF( fill ) THEN
          ijel = cell_g2l(ijmm, mpime)
          myij(0,-2,ijl) = ijel
        END IF

        ijpp=ij+nr+nr
        IF(j .EQ. (nz-1)) ijpp = ij+nr
        ipe =  cell_owner(ijpp) 
        IF( ipe .NE. mpime) THEN
          nset(ipe) = nset(ipe) + 1
          icnt = icnt + 1
          IF( fill ) THEN
            rcv_cell_set(ipe)%i(1,nset(ipe)) = ijpp
            rcv_cell_set(ipe)%i(2,nset(ipe)) = ijl
            rcv_cell_set(ipe)%i(3,nset(ipe)) = 0
            rcv_cell_set(ipe)%i(4,nset(ipe)) = 2
          END IF
        ELSE IF( fill ) THEN
          ijel = cell_g2l(ijpp, mpime)
          myij(0,2,ijl) = ijel
        END IF

        cell_neighbours = icnt

        RETURN
      END FUNCTION
!----------------------------------------------------------------------

      SUBROUTINE cell_set_sort(rcv_cell_set, nset, nrcv)
        USE basic_types, ONLY: imatrix
        TYPE (imatrix) :: rcv_cell_set(0:)
        INTEGER, INTENT(OUT) :: nrcv(0:)
        INTEGER, INTENT(IN) :: nset(0:)
        INTEGER :: ipe, nproc, icell(4), ic, it, icurr
        INTEGER, ALLOCATABLE :: ijsort(:), indx(:)
        nproc = SIZE(rcv_cell_set)
        DO ipe = 0, nproc - 1
          IF( nset(ipe) .GT. 0 ) THEN
            ALLOCATE( ijsort(nset(ipe)) )
            ALLOCATE( indx(nset(ipe)) )
            ijsort(:) = rcv_cell_set(ipe)%i(1,:)
            CALL ikb07ad(ijsort(1), nset(ipe), indx(1))
            DO ic = 1, nset(ipe)
              icurr = ic 
 3000         IF(indx(icurr) .NE. ic) THEN
                icell(:) = rcv_cell_set(ipe)%i(:,icurr) 
                rcv_cell_set(ipe)%i(:,icurr) = rcv_cell_set(ipe)%i(:,indx(icurr))
                rcv_cell_set(ipe)%i(:,indx(icurr)) = icell(:)
                it = icurr; icurr = indx(icurr); indx(it) = it
                IF(indx(icurr).EQ.ic) THEN
                  indx(icurr)=icurr
                ELSE
                  GOTO 3000
                END IF
              END IF
            END DO

            nrcv(ipe) = 1
            DO ic = 2, nset(ipe)
              IF( ijsort(ic) .NE. ijsort(ic-1) ) THEN
                nrcv(ipe) = nrcv(ipe) + 1
              END IF 
            END DO
            DEALLOCATE( ijsort, indx)
          END IF
        END DO

        RETURN
      END SUBROUTINE
!----------------------------------------------------------------------

      SUBROUTINE set_rcv_map(rcv_map, myij, rcv_cell_set, nset, nrcv)
        USE basic_types, ONLY: imatrix 
        TYPE (rcv_map_type) :: rcv_map(0:)
        TYPE (imatrix) :: rcv_cell_set(0:)
        INTEGER :: myij(-2:,-2:,:)
        INTEGER, INTENT(IN) :: nrcv(0:), nset(0:)
        INTEGER :: ipe, nproc, icnt, ijl, i, j, ije
        INTEGER :: ic, icnt_rcv
        INTEGER, ALLOCATABLE :: ijsort(:)

        nproc = SIZE(rcv_cell_set)
        nij_l = SIZE(myij,3)
        icnt  = 0

        DO ipe = 0, nproc - 1

          rcv_map(ipe)%nrcv = nrcv(ipe)
          IF( ASSOCIATED(rcv_map(ipe)%ijrcv) ) DEALLOCATE( rcv_map(ipe)%ijrcv )
          IF( ASSOCIATED(rcv_map(ipe)%ijl) ) DEALLOCATE( rcv_map(ipe)%ijl )

          IF( nrcv(ipe) .GT. 0 ) THEN
            ALLOCATE( rcv_map(ipe)%ijrcv(nrcv(ipe)) )
            ALLOCATE( rcv_map(ipe)%ijl(nrcv(ipe)) )

            ALLOCATE( ijsort(nset(ipe)) )
            ijsort(:) = rcv_cell_set(ipe)%i(1,:)

            ic        = 1
            icnt_rcv  = 1
            icnt      = icnt + 1
            ije = rcv_cell_set(ipe)%i(1,ic)
            ijl = rcv_cell_set(ipe)%i(2,ic)
            i   = rcv_cell_set(ipe)%i(3,ic)
            j   = rcv_cell_set(ipe)%i(4,ic)
            rcv_map(ipe)%ijrcv(icnt_rcv) = ije
            rcv_map(ipe)%ijl(icnt_rcv) = icnt + nij_l
            myij(i,j,ijl) = icnt + nij_l
            DO ic = 2, nset(ipe)
              ije = rcv_cell_set(ipe)%i(1,ic)
              ijl = rcv_cell_set(ipe)%i(2,ic)
              i   = rcv_cell_set(ipe)%i(3,ic)
              j   = rcv_cell_set(ipe)%i(4,ic)
              IF( ijsort(ic) .NE. ijsort(ic-1) ) THEN
                icnt_rcv  = icnt_rcv + 1
                icnt      = icnt + 1
                rcv_map(ipe)%ijrcv(icnt_rcv) = ije
                rcv_map(ipe)%ijl(icnt_rcv) = icnt + nij_l
              END IF 
              myij(i,j,ijl) = icnt + nij_l
            END DO

            DEALLOCATE( ijsort)

            IF( icnt_rcv .NE. nrcv(ipe) ) THEN
              WRITE(8,*) ' error in set_rcv_map - 1 '
              STOP
            END IF

          ELSE
            NULLIFY(rcv_map(ipe)%ijrcv)
            NULLIFY(rcv_map(ipe)%ijl)
          END IF

        END DO

        RETURN
      END SUBROUTINE

!---------------------------------------------------
      SUBROUTINE set_myinds(myinds, myij)
        USE dimensions
        IMPLICIT NONE
        INTEGER :: myinds(:,:)
        INTEGER :: myij(-2:,-2:,:)
!
        INTEGER :: ijl, ijr, ijb, ijt, ijbr, ijtr, ijbl, ijtl, ijrr, ijtt, ijll, ijbb
        INTEGER :: imj, ipj, ijm, ijp, ipjm, ipjp, imjm, imjp, ippj, ijpp, immj, ijmm

        INTEGER :: nflr, nflt, nfll, nflb
        INTEGER :: nfltr, nfltl, nflbr, nflbl
        INTEGER :: nflrr, nfltt, nflll, nflbb
        INTEGER :: i, j, ij, imesh
!
        DO ij = 1, nij_l 
          imesh = myij(0, 0, ij)
          j  = ( imesh - 1 ) / nr + 1
          i  = MOD( ( imesh - 1 ), nr) + 1
          IF( (i .GE. 2) .AND. (i .LE. (nr-1)) .AND.   &
              (j .GE. 2) .AND. (j .LE. (nz-1))      ) THEN
!
            ijm = myij( 0,-1, ij)
            imj = myij(-1, 0, ij)
            ipj = myij(+1, 0, ij)
            ijp = myij( 0,+1, ij)
            ipjm= myij(+1,-1, ij)
            ipjp= myij(+1,+1, ij)
            imjm= myij(-1,-1, ij)
            imjp= myij(-1,+1, ij)
            ippj= myij(+2, 0, ij)
            ijpp= myij( 0,+2, ij)
            immj= myij(-2, 0, ij)
            ijmm= myij( 0,-2, ij)
!
! ... indexes labelled Right, Top, Left, Bottom, etc., correspond to those defined 
! ...  above, except at boundaries ...
!
! ... First neighbours: near the axis or solid boundaries 
! ... impose homogeneous Neumann conditions
!
            ijr = ipj
            nflr = fl_l(ipj) 
            IF((nflr.EQ.2).OR.(nflr.EQ.3)) ijr = ij

            ijl = imj
            nfll = fl_l(imj) 
            IF((nfll.EQ.2).OR.(nfll.EQ.3)) ijl = ij

            ijt = ijp
            nflt = fl_l(ijp) 
            IF((nflt.EQ.2).OR.(nflt.EQ.3)) ijt = ij

            ijb = ijm
            nflb = fl_l(ijm) 
            IF((nflb.EQ.2).OR.(nflb.EQ.3)) ijb = ij

            ijtl = imjp
            nfltl=fl_l(imjp)
            IF(nfltl.EQ.2 .OR. nfltl.EQ.3) ijtl = ijp 

            ijbr = ipjm
            nflbr=fl_l(ipjm)
            IF(nflbr.EQ.2 .OR. nflbr.EQ.3) ijbr = ipj 

            ijtr = ipjp
            nfltr = fl_l(ipjp)
            IF (nfltr.EQ.2 .OR. nfltr.EQ.3) THEN
              nflt=fl_l(ijp)
              nflr=fl_l(ipj)
              IF (nflt.EQ.2 .OR. nflt.EQ.3) THEN
                IF (nflr.EQ.2 .OR. nflr.EQ.3) THEN
                  ijtr = ij 
                ELSE
                  ijtr = ipj
                END IF
              ELSE
                IF (nflr.EQ.2 .OR. nflr.EQ.3) THEN
                  ijtr = ijp 
                ELSE
                  ijtr = ij
                END IF
              END IF
            END IF
!
! ... Second neighbours are not available on boundaries
!
            ijrr = ippj
            IF(i .EQ. (nr-1)) ijrr = ipj
            nflrr=fl_l(ippj)
            IF( (nflrr.EQ.2) .OR. (nflrr.EQ.3) ) ijrr = ipj

            ijtt = ijpp
            IF(j.EQ.(nz-1)) ijtt = ijp
            nfltt=fl_l(ijpp)
            IF( (nfltt.EQ.2) .OR. (nfltt.EQ.3) ) ijtt = ijp

            ijll = immj
            IF(i .EQ. (2)) ijll = imj
            nflll=fl_l(immj)
            IF( (nflll.EQ.2) .OR. (nflll.EQ.3) ) ijll = imj

            ijbb = ijmm
            IF(j.EQ.(2)) ijbb = ijm
            nflbb=fl_l(ijmm)
            IF( (nflbb.EQ.2) .OR. (nflbb.EQ.3) ) ijbb = ijm
!
            myinds(bb_, ij) = ijbb
            myinds(bl_, ij) = ijbl
            myinds(b_,  ij) = ijb
            myinds(br_, ij) = ijbr
            myinds(ll_, ij) = ijll
            myinds(l_,  ij) = ijl
            myinds(r_,  ij) = ijr
            myinds(rr_, ij) = ijrr
            myinds(tl_, ij) = ijtl
            myinds(t_,  ij) = ijt
            myinds(tr_, ij) = ijtr
            myinds(tt_, ij) = ijtt

          END IF
        END DO

      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------

      SUBROUTINE data_exchange_r(array)
        USE parallel, ONLY: nproc, mpime
        IMPLICIT NONE
        REAL*8 :: array(:)
        REAL*8, ALLOCATABLE :: sndbuf(:), rcvbuf(:) 
        INTEGER :: ip, isour, idest, ib
        DO ip = 1, (nproc - 1)
          isour = MOD(mpime - ip + nproc, nproc)
          idest = MOD(mpime + ip        , nproc)
          ALLOCATE( rcvbuf( MAX(rcv_map(isour)%nrcv,1) ) )
          ALLOCATE( sndbuf( MAX(snd_map(idest)%nsnd,1) ) )
          DO ib = 1, snd_map(idest)%nsnd
            sndbuf(ib) = array( snd_map(idest)%ijl(ib) )
          END DO 
          CALL sendrecv_real(sndbuf(1), snd_map(idest)%nsnd, idest,  &
            rcvbuf, rcv_map(isour)%nrcv, isour, ip)
          DO ib = 1, rcv_map(isour)%nrcv
            array( rcv_map(isour)%ijl(ib) ) = rcvbuf(ib)
          END DO 
          DEALLOCATE( rcvbuf )
          DEALLOCATE( sndbuf )
        END DO
        RETURN
      END SUBROUTINE

      SUBROUTINE data_exchange_rm(array)
        USE parallel, ONLY: nproc, mpime
        IMPLICIT NONE
        REAL*8 :: array(:,:)
        REAL*8, ALLOCATABLE :: sndbuf(:), rcvbuf(:) 
        INTEGER :: ip, isour, idest, ib, ik, sdim, rdim
        DO ip = 1, (nproc - 1)
          isour = MOD(mpime - ip + nproc, nproc)
          idest = MOD(mpime + ip        , nproc)
          sdim = MAX( SIZE(array,1) * snd_map(idest)%nsnd, 1)
          rdim = MAX( SIZE(array,1) * rcv_map(isour)%nrcv, 1)
          ALLOCATE( rcvbuf( rdim ) )
          ALLOCATE( sndbuf( sdim ) )
!          rcvbuf = 0.0
!          sndbuf = 0.0
          DO ib = 1, snd_map(idest)%nsnd
            DO ik = 1, SIZE(array,1)
              sndbuf(ik + SIZE(array,1)*(ib-1)) = array( ik, snd_map(idest)%ijl(ib) )
            END DO
          END DO 
          CALL sendrecv_real(sndbuf(1), SIZE(sndbuf), idest,     &      
                             rcvbuf(1), SIZE(rcvbuf), isour, ip)
          DO ib = 1, rcv_map(isour)%nrcv
            DO ik = 1, SIZE(array,1)
              array( ik, rcv_map(isour)%ijl(ib) ) = rcvbuf( ik + SIZE(array,1)*(ib-1))
            END DO
          END DO 
          DEALLOCATE( rcvbuf )
          DEALLOCATE( sndbuf )
        END DO
        RETURN
      END SUBROUTINE

      SUBROUTINE data_exchange_i(array)
        USE parallel, ONLY: nproc, mpime
        IMPLICIT NONE
        INTEGER :: array(:)
        INTEGER, ALLOCATABLE :: sndbuf(:), rcvbuf(:)
        INTEGER :: ip, isour, idest, ib
        DO ip = 1, (nproc - 1)
          isour = MOD(mpime - ip + nproc, nproc)
          idest = MOD(mpime + ip        , nproc)
          ALLOCATE( rcvbuf( MAX(rcv_map(isour)%nrcv,1) ) )
          ALLOCATE( sndbuf( MAX(snd_map(idest)%nsnd,1) ) )
          DO ib = 1, snd_map(idest)%nsnd
            sndbuf(ib) = array( snd_map(idest)%ijl(ib) )
          END DO
          CALL sendrecv_integer(sndbuf, snd_map(idest)%nsnd, idest,      &
                                rcvbuf, rcv_map(isour)%nrcv, isour, ip)
          DO ib = 1, rcv_map(isour)%nrcv
            array( rcv_map(isour)%ijl(ib) ) = rcvbuf(ib)
          END DO
          DEALLOCATE( rcvbuf )
          DEALLOCATE( sndbuf )
        END DO
        RETURN
      END SUBROUTINE


      SUBROUTINE test_comm()
        IMPLICIT NONE
        INTEGER, ALLOCATABLE :: itest(:)
        INTEGER :: ijl
        
        ALLOCATE( itest(nijx_l) )

        itest = 0
        DO ijl = 1, nij_l
          itest( ijl ) = myij(0, 0, ijl)
        END DO

        CALL data_exchange(itest)
         
        WRITE(7,*)   ' index received ' 
        WRITE(7,310) itest( nij_l + 1 : nijx_l ) 
 310    FORMAT(10i8)
        
        DEALLOCATE( itest )
        
        RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      END MODULE
