!----------------------------------------------------------------------
   MODULE domain_mapping
!----------------------------------------------------------------------
        USE domain_decomposition, ONLY: nctot, proc_map
        USE domain_decomposition, ONLY: cell_g2l, cell_l2g, cell_owner
        USE domain_decomposition, ONLY: block2d_map, block3d_map, layer_map
        USE grid, ONLY: fl, flag
        USE grid, ONLY: slip_wall, noslip_wall, filled_cell_1, filled_cell_2
        USE kinds, ONLY: sgl, dbl
        USE indijk_module
        USE immersed_boundaries, ONLY: immb
        USE control_flags, ONLY: lpr
        USE io_files, ONLY: errorunit, testunit, logunit
!
        IMPLICIT NONE
!
        TYPE rcv_map_type
          INTEGER :: nrcv              ! How Many Cells I Have to Receive
          INTEGER, POINTER :: ircv(:)  ! Which Cells I Have to Receive (Global Index) 
          INTEGER, POINTER :: iloc(:)  ! Where I Have to Put Them (Local Index)
        END TYPE 
!
        TYPE snd_map_type
          INTEGER :: nsnd              ! How Many Cells I Have to Send
          INTEGER, POINTER :: isnd(:)  ! Which Cells I Have to Send (Global Index)
          INTEGER, POINTER :: iloc(:)  ! Where I Have to Pick Them (Local Index)
        END TYPE 
!
        TYPE (rcv_map_type), ALLOCATABLE :: rcv_map(:)
! ...     This array store, for each processor, the indexes of local
! ...     cell where to copy the data sent by other procesor

        TYPE (snd_map_type), ALLOCATABLE :: snd_map(:)
! ...     This array stores, for each processor, the indexes of local
! ...     cell that the present procesor has to send to the other

!       ncint (ncint) number of local cells, 
!       ncext         number of external (ghost) cells, 
!       ncdom         sum of the previous two
!
        INTEGER :: ncint, ncext, ncdom  
!
        INTEGER, ALLOCATABLE :: myijk(:,:)
        INTEGER, ALLOCATABLE :: myinds(:,:)

        INTERFACE data_exchange
          MODULE PROCEDURE data_exchange_i, data_exchange_r, data_exchange_rm, &
                           data_exchange_l
        END INTERFACE

        INTERFACE data_collect
          MODULE PROCEDURE data_collect_r, data_collect_sr
        END INTERFACE

        INTERFACE data_distribute
          MODULE PROCEDURE data_distribute_r, data_distribute_sr
        END INTERFACE
!
        PRIVATE
        PUBLIC :: ghost, meshinds
        PUBLIC :: myijk, myinds, ncint, ncdom
        PUBLIC :: data_exchange, data_collect, data_distribute
!
        SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE ghost
!
! ... Identifies and allocates ghost cells
! ... (2D/3D-Compliant)
! 
      USE dimensions
      USE basic_types, ONLY: imatrix
      USE control_flags, ONLY: job_type, prog
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE parallel, ONLY: nproc, mpime, root
!
      IMPLICIT NONE
!
      INTEGER :: icnt, ipe, icnt_ipe
      INTEGER :: i, j, n, ijkl, k, ijk
      INTEGER :: nrcv(0:nproc-1), nrcvx, jpe
      INTEGER :: nset(0:nproc-1)
      INTEGER, ALLOCATABLE :: ijkrcv(:,:)
      INTEGER, ALLOCATABLE :: ijksnd(:)
      TYPE (imatrix) :: rcv_cell_set(0:nproc-1)
      INTEGER :: localdim
      INTEGER :: numsend
      INTEGER :: layer, k2, k1, j2, j1, i2, i1, nkt
      INTEGER :: me, whose
      REAL*8  :: ratio
      LOGICAL :: lk, lj, li
!
      IF(ALLOCATED(rcv_map)) DEALLOCATE(rcv_map)
      IF(ALLOCATED(snd_map)) DEALLOCATE(snd_map)
      ALLOCATE(rcv_map( 0:(nproc-1) ) )
      ALLOCATE(snd_map( 0:(nproc-1) ) )
      
!
      rcv_map(:)%nrcv = 0
      snd_map(:)%nsnd = 0
      icnt     = 0
      icnt_ipe = 0
      ipe      = 0
!
! ... count the neighbouring cells on other processors:
! ... 'nset' is the number of the neighbouring cells on other processors,
! ... summed over each local cell.
! ... Only fluid cells, immersed_cells, and filled_cells are considered
!
      nset  = 0
      ncext = 0

      IF ( proc_map(mpime)%type == LAYER_MAP ) THEN
        DO ijk = proc_map(mpime)%lay(1), proc_map(mpime)%lay(2)
          IF ( BTEST(fl(ijk),0) .OR. BTEST(fl(ijk),9) .OR. BTEST(fl(ijk),10)) THEN
            ncext = ncext + cell_neighbours( ijk, mpime, nset)
          END IF
        END DO
      ELSE IF ( proc_map(mpime)%type == BLOCK2D_MAP ) THEN
        IF( job_type == JOB_TYPE_2D ) THEN
          i1 = proc_map(mpime)%corner1(1)
          i2 = proc_map(mpime)%corner2(1)
          k1 = proc_map(mpime)%corner1(2)
          k2 = proc_map(mpime)%corner2(2)
          DO k = k1, k2
            DO i = i1, i2
              ijk = i + (k-1) * nx
              IF ( BTEST(fl(ijk),0) .OR. BTEST(fl(ijk),9) .OR. BTEST(fl(ijk),10)) THEN
                ncext = ncext + cell_neighbours(ijk, mpime, nset)
              END IF
            END DO
          END DO
        ! ... "columns" domain-decomposition ... !
        ELSE IF( job_type == JOB_TYPE_3D ) THEN
          i1 = proc_map(mpime)%corner1(1)
          i2 = proc_map(mpime)%corner2(1)
          j1 = proc_map(mpime)%corner1(2)
          j2 = proc_map(mpime)%corner2(2)
          DO k = 1, nz
            DO j = j1, j2
              DO i = i1, i2
                ijk = i + (j-1)*nx + (k-1)*nx*ny 
                IF ( BTEST(fl(ijk),0) .OR. BTEST(fl(ijk),9) .OR. BTEST(fl(ijk),10)) THEN
                  ncext = ncext + cell_neighbours(ijk, mpime, nset)
                END IF
              END DO
            END DO
          END DO
        ELSE 
          CALL error( ' ghost ', ' wrong job_type ', job_type)
        END IF
      ELSE IF ( proc_map(mpime)%type == BLOCK3D_MAP ) THEN
        IF( job_type == JOB_TYPE_3D ) THEN
          i1 = proc_map(mpime)%blkbsw(1)
          i2 = proc_map(mpime)%blktne(1)
          j1 = proc_map(mpime)%blkbsw(2)
          j2 = proc_map(mpime)%blktne(2)
          k1 = proc_map(mpime)%blkbsw(3)
          k2 = proc_map(mpime)%blktne(3)
          DO k = k1, k2
            lk = .false.
            IF( ( k > k1 + 1 ) .AND. ( k < k2 - 1 ) ) lk = .true.
            DO j = j1, j2
              lj = .false.
              IF( ( j > j1 + 1 ) .AND. ( j < j2 - 1 ) ) lj = .true.
              DO i = i1, i2
                li = .false.
                IF( ( i > i1 + 1 ) .AND. ( i < i2 - 1 ) ) li = .true.
                ijk = i + (j-1)*nx + (k-1)*nx*ny
                IF ( BTEST(fl(ijk),0) .OR. BTEST(fl(ijk),9) .OR. BTEST(fl(ijk),10)) THEN
                  IF( li .AND. lj .AND. lk ) THEN
                    ncext = ncext + cell_neighbours(ijk, mpime, nset, is_my_cell = .true. )
                  ELSE
                    ncext = ncext + cell_neighbours(ijk, mpime, nset)
                  END IF
                END IF
              END DO
            END DO
          END DO
        ELSE 
          CALL error( ' ghost ', ' wrong job_type ', job_type)
        END IF
      ELSE
        CALL error(' ghost ', ' partition type not yet implemented ', proc_map(mpime)%type )
      END IF
!
! ... allocate memory for the set of neighbouring cell on other proc.
!
      DO ipe = 0, nproc - 1
        IF( nset(ipe) > 0 ) THEN
          ALLOCATE( rcv_cell_set(ipe)%i(5, nset(ipe)) )
        ELSE
          NULLIFY( rcv_cell_set(ipe)%i )
        END IF
      END DO

! ... number of local cells (includes boundary cells)
! ... as defined by the decomposition routine (layers or blocks)
!
      ncint = nctot( mpime )
!
! ... allocate the arrays of indexes: 'nstdim' is the stencil size 
! ... as defined in the indijk_module
!
      ALLOCATE( myijk( nstdim , ncint ) )
      ALLOCATE( myinds( nstdim, ncint ) )
!
! ... initialize arrays and set up indexes
!
      myijk   = 0
      myinds  = 0
      CALL indijk_setup()

! ... now fill in the receiving map and the indexes array 'myijk'
!
      nset  = 0
      icnt  = 0
      nrcv  = 0

      IF ( proc_map(mpime)%type == LAYER_MAP ) THEN
        DO ijk = proc_map(mpime)%lay(1), proc_map(mpime)%lay(2)
          IF ( BTEST(fl(ijk),0) .OR. BTEST(fl(ijk),9) .OR. BTEST(fl(ijk),10)) THEN
            icnt = icnt + cell_neighbours(ijk, mpime, nset, rcv_cell_set, myijk)
          END IF
        END DO
      ELSE IF ( proc_map(mpime)%type == BLOCK2D_MAP ) THEN
        IF( job_type == JOB_TYPE_2D ) THEN
          i1 = proc_map(mpime)%corner1(1)
          i2 = proc_map(mpime)%corner2(1)
          k1 = proc_map(mpime)%corner1(2)
          k2 = proc_map(mpime)%corner2(2)
          DO k = k1, k2
            DO i = i1, i2
              ijk = i + (k-1) * nx
              IF ( BTEST(fl(ijk),0) .OR. BTEST(fl(ijk),9) .OR. BTEST(fl(ijk),10)) THEN
                icnt = icnt + cell_neighbours(ijk, mpime, nset, rcv_cell_set, myijk)
              END IF
            END DO
          END DO
        ELSE IF( job_type == JOB_TYPE_3D ) THEN
          i1 = proc_map(mpime)%corner1(1)
          i2 = proc_map(mpime)%corner2(1)
          j1 = proc_map(mpime)%corner1(2)
          j2 = proc_map(mpime)%corner2(2)
          DO k = 1, nz
            DO j = j1, j2
              DO i = i1, i2
                ijk = i + (j-1)*nx + (k-1)*nx*ny
                IF ( BTEST(fl(ijk),0) .OR. BTEST(fl(ijk),9) .OR. BTEST(fl(ijk),10)) THEN
                  icnt = icnt + cell_neighbours(ijk, mpime, nset, rcv_cell_set, myijk)
                END IF
              END DO
            END DO
          END DO
        ELSE
          CALL error( ' ghost ', ' wrong job_type ', job_type)
        END IF
      ELSE IF ( proc_map(mpime)%type == BLOCK3D_MAP ) THEN
        IF( job_type == JOB_TYPE_3D ) THEN
          i1 = proc_map(mpime)%blkbsw(1)
          i2 = proc_map(mpime)%blktne(1)
          j1 = proc_map(mpime)%blkbsw(2)
          j2 = proc_map(mpime)%blktne(2)
          k1 = proc_map(mpime)%blkbsw(3)
          k2 = proc_map(mpime)%blktne(3)
          DO k = k1, k2
            lk = .false.
            IF( ( k > k1 + 1 ) .AND. ( k < k2 - 1 ) ) lk = .true.
            DO j = j1, j2
              lj = .false.
              IF( ( j > j1 + 1 ) .AND. ( j < j2 - 1 ) ) lj = .true.
              DO i = i1, i2
                li = .false.
                IF( ( i > i1 + 1 ) .AND. ( i < i2 - 1 ) ) li = .true.
                ijk = i + (j-1)*nx + (k-1)*nx*ny
                IF ( BTEST(fl(ijk),0) .OR. BTEST(fl(ijk),9) .OR. BTEST(fl(ijk),10)) THEN
                  IF( li .AND. lj .AND. lk ) THEN
                     icnt = icnt + cell_neighbours(ijk, mpime, nset, rcv_cell_set, myijk, is_my_cell = .true. )
                  ELSE
                     icnt = icnt + cell_neighbours(ijk, mpime, nset, rcv_cell_set, myijk)
                  END IF
                END IF
              END DO
            END DO
          END DO
        ELSE
          CALL error( ' ghost ', ' wrong job_type ', job_type)
        END IF
      ELSE
        CALL error(' ghost ', ' partition type not yet implemented ', &
                   proc_map(mpime)%type )
      END IF
!
! ... store the global index of the local cell ijk
!
      DO ijk = 1, ncint
        myijk( ip0_jp0_kp0_, ijk) = cell_l2g(ijk, mpime)
      END DO
!
! ... sort the neighbours cells to be received (this procedure
! ... ensures that all neighbour cells are received only once)
!
      CALL cell_set_sort(rcv_cell_set, nset, nrcv)

! ... 'nrcv' is the number of elements to be received from other procs
! ... and is equal to 'nset' once the copies of the cells has been eliminated
! ... 'rcv_cell_set' now contains, for every neighbour processor, the
! ... array of sorted global cell indexes to be sent to 'mpime' 
! .... (see test below)

      IF (lpr > 1) THEN
        DO ipe = 0, nproc - 1
          WRITE(testunit,300) nset(ipe), ipe
          IF ( nset(ipe) > 0 .AND. lpr > 2) THEN
            WRITE(testunit,310) rcv_cell_set(ipe)%i(1,:)
          END IF
 300      FORMAT(' # neighbours set SIZE ',i8,' from ',i3)
 310      FORMAT(10i8)
        END DO
      END IF

! ... the number cells required to update the physical quantities
! ... 'ncdom' is the sum of the local and neighbour cells
!
      ncext = SUM( nrcv ) 
      ncdom = ncint + ncext
      ratio = REAL(ncint,dbl)/ncext
      IF (lpr > 1) THEN
        WRITE(testunit,* ) ' # ncint ', ncint
        WRITE(testunit,* ) ' # ncext ', ncext
        WRITE(testunit,* ) ' # ncdom ', ncdom
        WRITE(testunit,* ) ' # ratio ', ratio
      END IF
!
! ... prepare the receive map 
!
      DO ipe = 0, nproc - 1
        rcv_map(ipe)%nrcv = nrcv(ipe)
        NULLIFY( rcv_map(ipe)%ircv )
        NULLIFY( rcv_map(ipe)%iloc )
      END DO
!
! ... prepare the send map 
!
      DO ipe = 0, nproc - 1
        NULLIFY( snd_map(ipe)%isnd ) 
        NULLIFY( snd_map(ipe)%iloc ) 
      END DO
!
! ... allocate and set the receiving map. This routine 
! ... assignes the ghost cells indexes.
!
      CALL set_rcv_map(rcv_map, myijk, rcv_cell_set, nset, nrcv)

! ... Each processor prints out the basic information on the map
!
      IF (lpr > 1) THEN
        DO ipe = 0, nproc - 1
          WRITE(testunit,*) ' # nrcv ', nrcv(ipe), ' from ', ipe
        END DO
      END IF
!
! ... using the receive maps fill in the sending maps
!
      DO ipe = 0, nproc - 1
!
        IF( ASSOCIATED (snd_map(ipe)%isnd ) ) THEN
          DEALLOCATE( snd_map(ipe)%isnd)
        END IF
        IF( ASSOCIATED (snd_map(ipe)%iloc ) ) THEN
          DEALLOCATE( snd_map(ipe)%iloc)
        END IF
        NULLIFY( snd_map(ipe)%isnd ) 
        NULLIFY( snd_map(ipe)%iloc ) 

! ...  proc 'ipe' send to all other processors the number of elements 
! ... 'nrcv' he must receive from each of them (rcv_map(:)%nrcv)
!
        IF( ipe .EQ. mpime ) THEN
          nrcv(:) = rcv_map(:)%nrcv
          nrcvx = MAXVAL(nrcv)
        END IF
!
! ...   each proc store the number of element he should send to ipe
!
        CALL scatter_integer(nrcv, snd_map(ipe)%nsnd, 1, ipe)
        CALL bcast_integer(nrcvx, 1, ipe)

! ...   proc ipe send to other processors the array with the indexes of 
! ...   the elements he should receive ( rcv_map(jpe)%ircv )
!
        ALLOCATE( ijkrcv(nrcvx,0:nproc-1) )
        ALLOCATE( ijksnd(nrcvx) )
!
        IF( ipe == mpime ) THEN
          ijkrcv = 0
          DO jpe = 0, nproc - 1
            IF( rcv_map(jpe)%nrcv > 0 ) THEN
              ijkrcv( 1:rcv_map(jpe)%nrcv, jpe ) = rcv_map(jpe)%ircv( 1:rcv_map(jpe)%nrcv )
            END IF
          END DO 
        END IF 
!
! ...   each proc store the index of the elements he should send to ipe
!
        CALL scatter_integer(ijkrcv, ijksnd, nrcvx, ipe)
!
        IF( ipe /= mpime ) THEN
          numsend = snd_map(ipe)%nsnd
          ALLOCATE( snd_map(ipe)%isnd( numsend ) )
          ALLOCATE( snd_map(ipe)%iloc( numsend ) )
          snd_map(ipe)%isnd( 1:numsend ) = ijksnd( 1:numsend )
          DO ijk = 1, numsend
            snd_map(ipe)%iloc( ijk ) = cell_g2l( ijksnd( ijk ), mpime )
          END DO
        END IF 
!
          DEALLOCATE(ijkrcv)
          DEALLOCATE(ijksnd)
      END DO
!
! ... Each processor prints out the information of the received cells
!
      IF (lpr > 1) THEN
        DO ipe = 0, nproc - 1
          WRITE(testunit,100) rcv_map(ipe)%nrcv, ipe
          IF( rcv_map(ipe)%nrcv > 0 .AND. lpr > 2 ) THEN
            WRITE(testunit,110) rcv_map(ipe)%ircv(:)
            WRITE(testunit,*) ' ---- '
            WRITE(testunit,110) rcv_map(ipe)%iloc(:)
          END IF
 100      FORMAT(' # receiving ',i5,' cells from ',i3)
 110      FORMAT(10i8)
        END DO
      END IF
!
! ... Each processor prints out the information of the sent cells
!
      IF (lpr > 1) THEN
        DO ipe = 0, nproc - 1
          WRITE(testunit,200) snd_map(ipe)%nsnd, ipe
          IF ( snd_map(ipe)%nsnd > 0 .AND. lpr > 2 ) THEN
            WRITE(testunit,210) snd_map(ipe)%isnd(:)
            WRITE(testunit,*) ' ---- '
            WRITE(testunit,210) snd_map(ipe)%iloc(:)
          END IF
 200      FORMAT(' # sending ',i5,' cells to ',i3)
 210      FORMAT(10i8)
        END DO
      END IF  

      IF (lpr > 2) CALL test_comm
!
! ... local flags, local arrays for forcing
!
!
      ALLOCATE( flag(ncdom) )
      DO ijkl = 1, ncint
        ijk = myijk( ip0_jp0_kp0_, ijkl)
        flag(ijkl) = fl(ijk)
      END DO
      DEALLOCATE (fl)

      CALL data_exchange(flag)
!
! ... Map the forcing points on local domains
! ... by using array numx/y/z. Scatter the array
! ... of forcing points among processors.
!
      IF (immb == 1 .AND. prog == 'PDAC') CALL local_forcing
!
! ... fill in the array myinds using myijk
!
      CALL set_myinds(myinds, myijk)
!
      RETURN
      END SUBROUTINE ghost
!----------------------------------------------------------------------
      INTEGER FUNCTION cell_neighbours(ijk, mpime, nset, rcv_cell_set, myijk, is_my_cell )
!
! ... loop over neighbour cells checking the owner. Fill in the 
! ... myijk array that identifies the ghost cells
!
        USE dimensions
        USE basic_types, ONLY: imatrix
        USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
!
        INTEGER, INTENT(IN) :: ijk
        INTEGER, INTENT(IN) :: mpime
        INTEGER, INTENT(INOUT) :: nset(0:)
        TYPE (imatrix), OPTIONAL :: rcv_cell_set(0:)
        INTEGER, OPTIONAL :: myijk(:,:)
        LOGICAL, OPTIONAL :: is_my_cell
!
        INTEGER ::  ii, jj, kk
        INTEGER :: icnt, ipe, i, j, k, ijke, ijkel, ijkl
        INTEGER :: im, jm, km
        LOGICAL :: fill, is_my_cell_

! ...   Subroutine Body

        IF ( PRESENT(rcv_cell_set) .AND. .NOT. present(myijk) ) &
          CALL error(' cell_neighbours ', ' in neighbour ', 1 )
!
        IF( .NOT. PRESENT(rcv_cell_set) .AND. present(myijk) ) &
          CALL error(' cell_neighbours ', ' in neighbour ', 2 )

        fill = PRESENT( rcv_cell_set )
        icnt = 0
        IF( fill ) THEN
          ncint = SIZE( myijk, 2 )
          ijkl  = cell_g2l( ijk, mpime )
        END IF

        is_my_cell_ = .false.
        IF( PRESENT( is_my_cell ) ) is_my_cell_ = is_my_cell

        IF( job_type == JOB_TYPE_2D ) THEN

          k = ( ijk - 1 ) / nx + 1
          i = MOD( ( ijk - 1 ), nx) + 1
        
! ...     loop over the first neighbouring cells
          DO im = -2, 2
            DO km = -2, 2
              IF( ( ABS( im ) + ABS( km ) ) <= 2 ) THEN
                IF( (im /= 0) .OR. (km /= 0) ) THEN
                  ii = im
                  kk = km 
                  IF( ( i == 2    ) .AND. ( ii == -2 ) ) ii = -1
                  IF( ( i == nx-1 ) .AND. ( ii == +2 ) ) ii = +1
                  IF( ( k == 2    ) .AND. ( kk == -2 ) ) kk = -1
                  IF( ( k == nz-1 ) .AND. ( kk == +2 ) ) kk = +1
                  ijke = ijk + ii + kk * nx
                  ipe =  cell_owner(ijke) 
                  IF( ipe /= mpime) THEN
! ...               the cell ijke is not local, count and register its position
                    nset(ipe) = nset(ipe) + 1
                    icnt = icnt + 1
                    IF( fill ) THEN
                      rcv_cell_set(ipe)%i(1,nset(ipe)) = ijke
                      rcv_cell_set(ipe)%i(2,nset(ipe)) = ijkl
                      rcv_cell_set(ipe)%i(3,nset(ipe)) = im
                      rcv_cell_set(ipe)%i(4,nset(ipe)) = 0
                      rcv_cell_set(ipe)%i(5,nset(ipe)) = km
                    END IF
                  ELSE IF( ipe==mpime .AND. fill ) THEN
! ...               the cell ijke is local, set the mapping with cell ijkl
                    ijkel = cell_g2l(ijke, mpime)
                    myijk( indijk( im, 0, km), ijkl ) = ijkel
                  END IF
                ELSE IF( fill ) THEN   
! ...             store the global cell index ijk, of the the local cell ijkl
                  myijk( ip0_jp0_kp0_, ijkl) = ijk
                END IF
              END IF
            END DO
          END DO

        ELSE IF( job_type == JOB_TYPE_3D ) THEN

          i = MOD( MOD( ijk - 1, nx*ny ), nx ) + 1
          j = MOD( ijk - 1, nx*ny ) / nx + 1
          k = ( ijk - 1 ) / ( nx*ny ) + 1

          IF( fill ) THEN
! ...        store the global cell index ij, of the the local cell ijkl
             myijk( ip0_jp0_kp0_, ijkl) = ijk
          END IF
          
! ...     loop over the neighbouring cells

          IF( is_my_cell_ ) THEN

             IF( fill ) THEN
               DO km = -2, 2
                 kk = km
                 IF( ( k == 2    ) .AND. ( kk == -2 ) ) kk = -1
                 IF( ( k == nz-1 ) .AND. ( kk == +2 ) ) kk = +1
                 DO jm = -2, 2
                   jj = jm
                   IF( ( j == 2    ) .AND. ( jj == -2 ) ) jj = -1
                   IF( ( j == ny-1 ) .AND. ( jj == +2 ) ) jj = +1
                   DO im = -2, 2
                     IF( ( ABS( im ) + ABS( jm ) + ABS( km ) ) <= 2 ) THEN
                       IF( (im /= 0) .OR. (jm /= 0) .OR. (km /= 0) ) THEN
                         ii = im
                         IF( ( i == 2    ) .AND. ( ii == -2 ) ) ii = -1
                         IF( ( i == nx-1 ) .AND. ( ii == +2 ) ) ii = +1
                         ijke = ijk + ii + jj * nx + kk * nx*ny
   ! ...                 the cell ijke is local, set the mapping with cell ijkl
                         ijkel = cell_g2l(ijke, mpime)
                         myijk( indijk( im, jm, km ), ijkl ) = ijkel
                       END IF
                     END IF
                   END DO
                 END DO
               END DO
             END IF


          ELSE


             DO km = -2, 2
               kk = km
               IF( ( k == 2    ) .AND. ( kk == -2 ) ) kk = -1
               IF( ( k == nz-1 ) .AND. ( kk == +2 ) ) kk = +1
               DO jm = -2, 2
                 jj = jm
                 IF( ( j == 2    ) .AND. ( jj == -2 ) ) jj = -1
                 IF( ( j == ny-1 ) .AND. ( jj == +2 ) ) jj = +1
                 DO im = -2, 2
                   IF( ( ABS( im ) + ABS( jm ) + ABS( km ) ) <= 2 ) THEN
                     IF( (im /= 0) .OR. (jm /= 0) .OR. (km /= 0) ) THEN
                       ii = im
                       IF( ( i == 2    ) .AND. ( ii == -2 ) ) ii = -1
                       IF( ( i == nx-1 ) .AND. ( ii == +2 ) ) ii = +1
                       ijke = ijk + ii + jj * nx + kk * nx*ny
                       ipe  =  cell_owner(ijke)
                       IF( ipe /= mpime ) THEN
   ! ...                 the cell ijke is not local, count and register its position
                         nset(ipe) = nset(ipe) + 1
                         icnt = icnt + 1
                         IF( fill ) THEN
                           rcv_cell_set(ipe)%i(1,nset(ipe)) = ijke
                           rcv_cell_set(ipe)%i(2,nset(ipe)) = ijkl
                           rcv_cell_set(ipe)%i(3,nset(ipe)) = im
                           rcv_cell_set(ipe)%i(4,nset(ipe)) = jm
                           rcv_cell_set(ipe)%i(5,nset(ipe)) = km
                         END IF
                       ELSE IF( fill ) THEN
   ! ...                 the cell ijke is local, set the mapping with cell ijkl
                         ijkel = cell_g2l(ijke, mpime)
                         myijk( indijk( im, jm, km ), ijkl ) = ijkel
                       END IF
                     END IF
                   END IF
                 END DO
               END DO
             END DO


          END IF


        ELSE 

          CALL error(' cell_neighbours ', ' wrong job_type ', 1 )

        END IF

        cell_neighbours = icnt

        RETURN
      END FUNCTION cell_neighbours
!----------------------------------------------------------------------
      SUBROUTINE cell_set_sort(rcv_cell_set, nset, nrcv)
!
        USE basic_types, ONLY: imatrix

        TYPE (imatrix) :: rcv_cell_set(0:)
        INTEGER, INTENT(OUT) :: nrcv(0:)
        INTEGER, INTENT(IN) :: nset(0:)
        INTEGER :: ipe, nproc, icell(5), ic, it, icurr
        INTEGER, ALLOCATABLE :: ijksort(:), indx(:)
!
        nproc = SIZE(rcv_cell_set)
        DO ipe = 0, nproc - 1
          IF( nset(ipe) .GT. 0 ) THEN
            ALLOCATE( ijksort(nset(ipe)) )
            ALLOCATE( indx(nset(ipe)) )
            ijksort(:) = rcv_cell_set(ipe)%i(1,:)
            CALL ikb07ad(ijksort(1), nset(ipe), indx(1))
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
              IF( ijksort(ic) .NE. ijksort(ic-1) ) THEN
                nrcv(ipe) = nrcv(ipe) + 1
              END IF 
            END DO
            DEALLOCATE( ijksort, indx)
          END IF
        END DO

        RETURN
      END SUBROUTINE cell_set_sort
!
!----------------------------------------------------------------------
!
      SUBROUTINE set_rcv_map(rcv_map, myijk, rcv_cell_set, nset, nrcv)
!
        USE basic_types, ONLY: imatrix 
!
        TYPE (rcv_map_type) :: rcv_map(0:)
        TYPE (imatrix) :: rcv_cell_set(0:)
        INTEGER :: myijk(:,:)
        INTEGER, INTENT(IN) :: nrcv(0:), nset(0:)
        INTEGER :: ipe, nproc, icnt, ijkl, i, j, k, ijke
        INTEGER :: ic, icnt_rcv
        INTEGER, ALLOCATABLE :: ijksort(:)

        nproc = SIZE( rcv_cell_set )
        ncint = SIZE( myijk, 2 )


        icnt  = 0    !  counts received cells from all processors

        DO ipe = 0, nproc - 1

          rcv_map(ipe)%nrcv = nrcv(ipe)

          IF( ASSOCIATED(rcv_map(ipe)%ircv) ) DEALLOCATE( rcv_map(ipe)%ircv )
          IF( ASSOCIATED(rcv_map(ipe)%iloc) ) DEALLOCATE( rcv_map(ipe)%iloc )

          IF( nrcv(ipe) > 0 ) THEN

            ALLOCATE( rcv_map(ipe)%ircv(nrcv(ipe)) )
            ALLOCATE( rcv_map(ipe)%iloc(nrcv(ipe)) )

            ALLOCATE( ijksort( nset(ipe) ) )

            ijksort(:) = rcv_cell_set(ipe)%i(1,:)

! ...       Manage the first cell to be received

            ic        = 1           ! index that run over neighbour cells 
                                    !   of the current processor 
            icnt_rcv  = 1           ! counts the received cells from the 
                                    !   current processor (ipe)
            icnt      = icnt + 1    ! increment the index of cell

            ijke = rcv_cell_set(ipe)%i(1,ic)
            ijkl = rcv_cell_set(ipe)%i(2,ic)
            i    = rcv_cell_set(ipe)%i(3,ic)
            j    = rcv_cell_set(ipe)%i(4,ic)
            k    = rcv_cell_set(ipe)%i(5,ic)

            rcv_map(ipe)%ircv(icnt_rcv)  = ijke
            rcv_map(ipe)%iloc(icnt_rcv)  = ncint + icnt
            myijk( indijk(i,j,k), ijkl ) = ncint + icnt

! ...       Manage the remaining cells

            DO ic = 2, nset(ipe)

              ijke = rcv_cell_set(ipe)%i(1,ic)
              ijkl = rcv_cell_set(ipe)%i(2,ic)
              i    = rcv_cell_set(ipe)%i(3,ic)
              j    = rcv_cell_set(ipe)%i(4,ic)
              k    = rcv_cell_set(ipe)%i(5,ic)

! ...         Here add to the list only new cells

              IF( ijksort( ic ) /= ijksort( ic-1 ) ) THEN
                icnt_rcv  = icnt_rcv + 1
                icnt      = icnt + 1
                rcv_map(ipe)%ircv(icnt_rcv) = ijke
                rcv_map(ipe)%iloc(icnt_rcv) = ncint + icnt
              END IF 

              myijk( indijk(i,j,k), ijkl ) = ncint + icnt

            END DO

            DEALLOCATE( ijksort )

            IF( icnt_rcv /= nrcv(ipe) ) &
              CALL error( '  set_rcv_map ', ' inconsistent number of cells ', &
                         icnt_rcv )

          ELSE

            NULLIFY( rcv_map(ipe)%ircv )
            NULLIFY( rcv_map(ipe)%iloc )

          END IF

        END DO

        RETURN
      END SUBROUTINE set_rcv_map
!----------------------------------------------------------------------
      SUBROUTINE meshinds(localindex,globalindex,i,j,k)
!
! ... computes the coordinates of a grid point in 
! ... a 2D or 3D rectilinear mesh

        USE dimensions
        USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: localindex
        INTEGER, INTENT(OUT) :: globalindex, i, j, k

        i = 0    ! index for   x(r)   coordinate
        j = 0    ! index for   y      coordinate
        k = 0    ! index for   z      coordinate

        globalindex = myijk( ip0_jp0_kp0_, localindex)

        IF (job_type == JOB_TYPE_3D) THEN

!          nxy = nx * ny
!          im1 = imesh - 1
!          mdk = im1/nxy
!          mj  = im1 - nxy * mdk
!          mdj = mj/nx
!          mi  = mj  -  nx * mdj
!
!          i   = mi  + 1      ! i = MOD( MOD( imesh - 1, nx*ny ), nx ) + 1
!          j   = mdj + 1      ! j = MOD( imesh - 1, nx*ny ) / nx + 1
!          k   = mdk + 1      ! k = ( imesh - 1 ) / ( nx*ny ) + 1
          
          k = ( globalindex - 1 ) / ( nx*ny ) + 1
          j = MOD( globalindex - 1, nx*ny) / nx + 1
          i = MOD( MOD( globalindex - 1, nx*ny ), nx ) + 1

        ELSE IF (job_type == JOB_TYPE_2D) THEN

          k = ( globalindex - 1 ) / nx + 1
          i = MOD( ( globalindex - 1 ), nx) + 1

        ELSE 

          CALL error(' meshinds ', ' unknow job_type ', job_type )

        END IF
!
      RETURN
      END SUBROUTINE meshinds
!---------------------------------------------------
      SUBROUTINE set_myinds(myinds, myijk)
!
        USE dimensions
        USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
!
        IMPLICIT NONE
        INTEGER :: myinds(:,:)
        INTEGER :: myijk(:,:)
!
        INTEGER :: i, j, k, ijk, imesh

        INTEGER ::  ipjk, imjk, ippjk, immjk, ijpk, ipjpk, imjpk, ijmk,  &
                    ipjmk, imjmk, ijppk, ijmmk, ijkp, ipjkp, imjkp, ijpkp,  &
                    ijmkp, ijkm, ipjkm, imjkm, ijpkm, ijmkm, ijkpp, ijkmm

        INTEGER ::  ijke, ijkw, ijkee, ijkww, ijkn, ijken, ijkwn, ijks,  &
                    ijkes,  ijkws, ijknn, ijkss, ijkt, ijket, ijkwt, ijknt, &
                    ijkst, ijkb, ijkeb, ijkwb, ijknb, ijksb, ijktt, ijkbb
!
        DO ijk = 1, ncint 
          CALL meshinds(ijk,imesh,i,j,k)

          IF( job_type == JOB_TYPE_2D ) THEN

            IF( (i >= 2) .AND. (i <= (nx-1)) .AND.   &
                (k >= 2) .AND. (k <= (nz-1))      ) THEN
!
              ijkm  = myijk( ip0_jp0_km1_, ijk)
              imjk  = myijk( im1_jp0_kp0_, ijk)
              ipjk  = myijk( ip1_jp0_kp0_, ijk)
              ijkp  = myijk( ip0_jp0_kp1_, ijk)
              ipjkm = myijk( ip1_jp0_km1_, ijk)
              ipjkp = myijk( ip1_jp0_kp1_, ijk)
              imjkm = myijk( im1_jp0_km1_, ijk)
              imjkp = myijk( im1_jp0_kp1_, ijk)
              ippjk = myijk( ip2_jp0_kp0_, ijk)
              ijkpp = myijk( ip0_jp0_kp2_, ijk)
              immjk = myijk( im2_jp0_kp0_, ijk)
              ijkmm = myijk( ip0_jp0_km2_, ijk)
!
! ... First neighbours: near the axis or solid boundaries 
! ... impose homogeneous Neumann conditions
!
              SELECT CASE (flag(ipjk))
                CASE (noslip_wall, slip_wall, filled_cell_1)
                        ijke = ijk
                CASE DEFAULT
                        ijke = ipjk
              END SELECT

              SELECT CASE (flag(imjk))
                CASE (noslip_wall, slip_wall, filled_cell_1)
                        ijkw = ijk
                CASE DEFAULT
                        ijkw = imjk
              END SELECT

              SELECT CASE (flag(ijkp))
                CASE (noslip_wall, slip_wall, filled_cell_1)
                        ijkt = ijk
                CASE DEFAULT
                        ijkt = ijkp
              END SELECT

              SELECT CASE (flag(ijkm))
                CASE (noslip_wall, slip_wall, filled_cell_1)
                        ijkb = ijk
                CASE DEFAULT
                        ijkb = ijkm
              END SELECT
!
! ... diagonal neighbours
!
              ijkwt = imjkp
              IF(flag(imjkp) == slip_wall .OR. flag(imjkp) == noslip_wall) ijkwt = ijkp

              ijkeb = ipjkm
              IF(flag(ipjkm) == slip_wall .OR. flag(ipjkm) == noslip_wall) ijkeb = ipjk

              ijket = ipjkp 
              IF (flag(ipjkp) == slip_wall .OR. flag(ipjkp) == noslip_wall) THEN
                IF (flag(ijkp) == slip_wall .OR. flag(ijkp) == noslip_wall) THEN
                  IF (flag(ipjk) == slip_wall .OR. flag(ipjk) == noslip_wall) THEN
                    ijket = ijk 
                  ELSE
                    ijket = ipjk
                  END IF
                ELSE
                  IF (flag(ipjk) == slip_wall .OR. flag(ipjk) == noslip_wall) THEN
                    ijket = ijkp 
                  ELSE
                    ijket = ijk
                  END IF
                END IF
              END IF
              
              ijkwb = imjkm
!
! ... Second neighbours are not available on boundaries
!
              SELECT CASE (flag(ippjk))
              CASE (slip_wall, noslip_wall, filled_cell_1)
                      ijkee = ipjk
              CASE DEFAULT
                      ijkee = ippjk
              END SELECT
              IF(i == (nx-1)) ijkee = ijke
!
              SELECT CASE (flag(ijkpp))
              CASE (slip_wall, noslip_wall, filled_cell_1)
                      ijktt = ijkp
              CASE DEFAULT
                      ijktt = ijkpp
              END SELECT
              IF(k == (nz-1)) ijktt = ijkt
!
              SELECT CASE (flag(immjk))
              CASE (slip_wall, noslip_wall, filled_cell_1)
                      ijkww = imjk
              CASE DEFAULT
                      ijkww = immjk
              END SELECT
              IF(i == 2) ijkww = ijkw
!
              SELECT CASE (flag(ijkmm))
              CASE (slip_wall, noslip_wall, filled_cell_1)
                      ijkbb = ijkm
              CASE DEFAULT
                      ijkbb = ijkmm
              END SELECT
              IF(k == 2) ijkbb = ijkb
!
              myinds(ip0_jp0_km2_, ijk) = ijkbb
              myinds(im1_jp0_km1_, ijk) = ijkwb
              myinds(ip0_jp0_km1_,  ijk) = ijkb
              myinds(ip1_jp0_km1_, ijk) = ijkeb
              myinds(im2_jp0_kp0_, ijk) = ijkww
              myinds(im1_jp0_kp0_,  ijk) = ijkw
              myinds(ip1_jp0_kp0_,  ijk) = ijke
              myinds(ip2_jp0_kp0_, ijk) = ijkee
              myinds(im1_jp0_kp1_, ijk) = ijkwt
              myinds(ip0_jp0_kp1_,  ijk) = ijkt
              myinds(ip1_jp0_kp1_, ijk) = ijket
              myinds(ip0_jp0_kp2_, ijk) = ijktt

            END IF

          ELSE IF( job_type == JOB_TYPE_3D ) THEN

            IF( (i >= 2) .AND. (i <= (nx-1)) .AND.   &
                (j >= 2) .AND. (j <= (ny-1)) .AND.   &
                (k >= 2) .AND. (k <= (nz-1))         ) THEN
!
! ... First neighbours: near the axis or solid boundaries 
! ... impose homogeneous Neumann conditions
!
              ipjk   = myijk( ip1_jp0_kp0_ , ijk )
              imjk   = myijk( im1_jp0_kp0_ , ijk )
              ippjk  = myijk( ip2_jp0_kp0_ , ijk )
              immjk  = myijk( im2_jp0_kp0_ , ijk )
              ijpk   = myijk( ip0_jp1_kp0_ , ijk )
              ipjpk  = myijk( ip1_jp1_kp0_ , ijk )
              imjpk  = myijk( im1_jp1_kp0_ , ijk )
              ijmk   = myijk( ip0_jm1_kp0_ , ijk )
              ipjmk  = myijk( ip1_jm1_kp0_ , ijk )
              imjmk  = myijk( im1_jm1_kp0_ , ijk )
              ijppk  = myijk( ip0_jp2_kp0_ , ijk )
              ijmmk  = myijk( ip0_jm2_kp0_ , ijk )
              ijkp   = myijk( ip0_jp0_kp1_ , ijk )
              ipjkp  = myijk( ip1_jp0_kp1_ , ijk )
              imjkp  = myijk( im1_jp0_kp1_ , ijk )
              ijpkp  = myijk( ip0_jp1_kp1_ , ijk )
              ijmkp  = myijk( ip0_jm1_kp1_ , ijk )
              ijkm   = myijk( ip0_jp0_km1_ , ijk )
              ipjkm  = myijk( ip1_jp0_km1_ , ijk )
              imjkm  = myijk( im1_jp0_km1_ , ijk )
              ijpkm  = myijk( ip0_jp1_km1_ , ijk )
              ijmkm  = myijk( ip0_jm1_km1_ , ijk )
              ijkpp  = myijk( ip0_jp0_kp2_ , ijk )
              ijkmm  = myijk( ip0_jp0_km2_ , ijk )
!  
              SELECT CASE (flag(ipjk))
                !CASE (noslip_wall, slip_wall, filled_cell_1, filled_cell_2)
                CASE (noslip_wall, slip_wall, filled_cell_1)
                        ijke = ijk
                CASE DEFAULT
                        ijke  =  ipjk
              END SELECT
!
              SELECT CASE (flag(imjk))
                !CASE (noslip_wall, slip_wall, filled_cell_1, filled_cell_2)
                CASE (noslip_wall, slip_wall, filled_cell_1)
                        ijkw = ijk
                CASE DEFAULT
                        ijkw  =  imjk
              END SELECT
!
              SELECT CASE (flag(ijpk))
                !CASE (noslip_wall, slip_wall, filled_cell_1, filled_cell_2)
                CASE (noslip_wall, slip_wall, filled_cell_1)
                        ijkn = ijk
                CASE DEFAULT
                        ijkn  =  ijpk
              END SELECT
!
              SELECT CASE (flag(ijmk))
                !CASE (noslip_wall, slip_wall, filled_cell_1, filled_cell_2)
                CASE (noslip_wall, slip_wall, filled_cell_1)
                        ijks = ijk
                CASE DEFAULT
                        ijks  =  ijmk
              END SELECT
!
              SELECT CASE (flag(ijkp))
                !CASE (noslip_wall, slip_wall, filled_cell_1, filled_cell_2)
                CASE (noslip_wall, slip_wall, filled_cell_1)
                        ijkt = ijk
                CASE DEFAULT
                        ijkt  =  ijkp
              END SELECT
!
              SELECT CASE (flag(ijkm))
                !CASE (noslip_wall, slip_wall, filled_cell_1, filled_cell_2)
                CASE (noslip_wall, slip_wall, filled_cell_1)
                        ijkb = ijk
                CASE DEFAULT
                        ijkb  =  ijkm
              END SELECT
!
! ... Diagonal neighbours
!
!              SELECT CASE (flag(ipjpk))
!              SELECT CASE (flag(imjpk))
!              SELECT CASE (flag(ipjmk))
!              SELECT CASE (flag(imjmk))
!              SELECT CASE (flag(ipjkp))
!              SELECT CASE (flag(imjkp))
!              SELECT CASE (flag(ijpkp))
!              SELECT CASE (flag(ijmkp))
!              SELECT CASE (flag(ipjkm))
!              SELECT CASE (flag(imjkm))
!              SELECT CASE (flag(ijpkm))
!              SELECT CASE (flag(ijmkm))
!
! ... Second neighbours are not available on boundaries
!
              SELECT CASE (flag(ippjk))
              !CASE (slip_wall, noslip_wall, filled_cell_1, filled_cell_2)
              CASE (slip_wall, noslip_wall, filled_cell_1)
                      ijkee = ipjk
              CASE DEFAULT
                      ijkee = ippjk
              END SELECT
              IF(i == (nx-1)) ijkee = ijke
!
              SELECT CASE (flag(ijkpp))
              !CASE (slip_wall, noslip_wall, filled_cell_1, filled_cell_2)
              CASE (slip_wall, noslip_wall, filled_cell_1)
                      ijktt = ijkp
              CASE DEFAULT
                      ijktt = ijkpp
              END SELECT
              IF(k == (nz-1)) ijktt = ijkt
!
              SELECT CASE (flag(ijppk))
              !CASE (slip_wall, noslip_wall, filled_cell_1, filled_cell_2)
              CASE (slip_wall, noslip_wall, filled_cell_1)
                      ijknn = ijkn
              CASE DEFAULT
                      ijknn =  ijppk
              END SELECT
              if( (j == (ny-1)) ) ijknn = ijkn

              SELECT CASE (flag(ijmmk))
              !CASE (slip_wall, noslip_wall, filled_cell_1, filled_cell_2)
              CASE (slip_wall, noslip_wall, filled_cell_1)
                      ijkss = ijks
              CASE DEFAULT
                      ijkss =  ijmmk
              END SELECT
              if( (j == 2) ) ijkss = ijks
!
              SELECT CASE (flag(immjk))
              !CASE (slip_wall, noslip_wall, filled_cell_1, filled_cell_2)
              CASE (slip_wall, noslip_wall, filled_cell_1)
                      ijkww = imjk
              CASE DEFAULT
                      ijkww = immjk
              END SELECT
              IF(i == 2) ijkww = ijkw
!
              SELECT CASE (flag(ijkmm))
              !CASE (slip_wall, noslip_wall, filled_cell_1, filled_cell_2)
              CASE (slip_wall, noslip_wall, filled_cell_1)
                      ijkbb = ijkm
              CASE DEFAULT
                      ijkbb = ijkmm
              END SELECT
              IF(k == 2) ijkbb = ijkb
!
              !!!!! check diagonals
              !
              ijken =  ipjpk
              ijkwn =  imjpk
              ijkes =  ipjmk
              ijkws =  imjmk
              ijket =  ipjkp
              ijkwt =  imjkp
              ijknt =  ijpkp
              ijkst =  ijmkp
              ijkeb =  ipjkm
              ijkwb =  imjkm
              ijknb =  ijpkm
              ijksb =  ijmkm
!
              myinds( ip1_jp0_kp0_ , ijk ) = ijke
              myinds( im1_jp0_kp0_ , ijk ) = ijkw
              myinds( ip2_jp0_kp0_ , ijk ) = ijkee
              myinds( im2_jp0_kp0_ , ijk ) = ijkww
              myinds( ip0_jp1_kp0_ , ijk ) = ijkn
              myinds( ip1_jp1_kp0_ , ijk ) = ijken
              myinds( im1_jp1_kp0_ , ijk ) = ijkwn
              myinds( ip0_jm1_kp0_ , ijk ) = ijks
              myinds( ip1_jm1_kp0_ , ijk ) = ijkes
              myinds( im1_jm1_kp0_ , ijk ) = ijkws
              myinds( ip0_jp2_kp0_ , ijk ) = ijknn
              myinds( ip0_jm2_kp0_ , ijk ) = ijkss
              myinds( ip0_jp0_kp1_ , ijk ) = ijkt
              myinds( ip1_jp0_kp1_ , ijk ) = ijket
              myinds( im1_jp0_kp1_ , ijk ) = ijkwt
              myinds( ip0_jp1_kp1_ , ijk ) = ijknt
              myinds( ip0_jm1_kp1_ , ijk ) = ijkst
              myinds( ip0_jp0_km1_ , ijk ) = ijkb
              myinds( ip1_jp0_km1_ , ijk ) = ijkeb
              myinds( im1_jp0_km1_ , ijk ) = ijkwb
              myinds( ip0_jp1_km1_ , ijk ) = ijknb
              myinds( ip0_jm1_km1_ , ijk ) = ijksb
              myinds( ip0_jp0_kp2_ , ijk ) = ijktt
              myinds( ip0_jp0_km2_ , ijk ) = ijkbb
  
            END IF

          ELSE
 
            CALL error(' set_myinds ',' unknown job_type ',1) 
            
          END IF

        END DO

      RETURN
      END SUBROUTINE set_myinds
!
!----------------------------------------------------------------------
!
      SUBROUTINE data_exchange_r(array)
!
        USE parallel, ONLY: nproc, mpime
!
        IMPLICIT NONE
        REAL*8 :: array(:)
        REAL*8, ALLOCATABLE :: sndbuf(:), rcvbuf(:) 
        INTEGER :: ip, isour, idest, ib, itag, ishand, irhand
        INTEGER :: ierr
!
        DO ip = 1, (nproc - 1)

          isour = MOD(mpime - ip + nproc, nproc)
          idest = MOD(mpime + ip        , nproc)
          ALLOCATE( rcvbuf( MAX(rcv_map(isour)%nrcv,1) ), STAT=ierr )

          IF( ierr /= 0 ) THEN
            IF (lpr > 1) WRITE(testunit,*) 'Trying to allocate ', MAX(rcv_map(isour)%nrcv,1), ' elements '
            CALL error(' data_exchange_r ', ' allocating rcvbuf ', ierr )
          END IF
          ALLOCATE( sndbuf( MAX(snd_map(idest)%nsnd,1) ), STAT=ierr )
          IF( ierr /= 0 ) THEN
            IF (lpr > 1) WRITE(testunit,*) 'Trying to allocate ', MAX(snd_map(idest)%nsnd,1), ' elements '
            CALL error(' data_exchange_r ', ' allocating sndbuf ', ierr )
          END IF

          IF( rcv_map(isour)%nrcv > 0 ) THEN
            CALL irecv_real( rcvbuf, rcv_map(isour)%nrcv, isour, ip, irhand )
          END IF

          DO ib = 1, snd_map(idest)%nsnd
            sndbuf(ib) = array( snd_map(idest)%iloc(ib) )
          END DO 

          IF( snd_map(idest)%nsnd > 0 ) THEN
            CALL isend_real( sndbuf(1), snd_map(idest)%nsnd, idest, ip, ishand )
          END IF

          !CALL sendrecv_real(sndbuf(1), snd_map(idest)%nsnd, idest,  &
          !   rcvbuf, rcv_map(isour)%nrcv, isour, ip)

          IF( rcv_map(isour)%nrcv > 0 ) THEN
            CALL mp_wait( irhand )
          END IF

          DO ib = 1, rcv_map(isour)%nrcv
            array( rcv_map(isour)%iloc(ib) ) = rcvbuf(ib)
          END DO 

          IF( snd_map(idest)%nsnd > 0 ) THEN
            CALL mp_wait( ishand )
          END IF

          DEALLOCATE( rcvbuf, STAT=ierr )
          IF( ierr /= 0 ) &
            CALL error(' data_exchange_r ', ' deallocating rcvbuf ', ierr )
          DEALLOCATE( sndbuf, STAT=ierr )
          IF( ierr /= 0 ) &
            CALL error(' data_exchange_r ', ' deallocating sndbuf ', ierr )
        END DO
!
        RETURN
      END SUBROUTINE data_exchange_r
!
!----------------------------------------------------------------------
!
      SUBROUTINE data_exchange_rm(array)
        USE parallel, ONLY: nproc, mpime
        IMPLICIT NONE
        REAL*8 :: array(:,:)
        REAL*8, ALLOCATABLE :: sndbuf(:), rcvbuf(:) 
        INTEGER :: ip, isour, idest, ib, ik, sdim, rdim, ierr
        INTEGER :: itag, ishand, irhand

        DO ip = 1, (nproc - 1)

          isour = MOD(mpime - ip + nproc, nproc)
          idest = MOD(mpime + ip        , nproc)

          sdim = MAX( SIZE(array,2) * snd_map(idest)%nsnd, 1)
          rdim = MAX( SIZE(array,2) * rcv_map(isour)%nrcv, 1)
          ALLOCATE( rcvbuf( rdim ), STAT=ierr )
          IF( ierr /= 0 ) THEN
            IF (lpr > 1) WRITE(testunit,*) 'Trying to allocate ', rdim, ' elements '
            CALL error(' data_exchange_rm ', ' allocating rcvbuf ', ierr )
          END IF
          ALLOCATE( sndbuf( sdim ), STAT=ierr )
          IF( ierr /= 0 ) THEN
            IF (lpr > 1) WRITE(testunit,*) 'Trying to allocate ', sdim, ' elements '
            CALL error(' data_exchange_rm ', ' allocating sndbuf ', ierr )
          END IF


          IF( rcv_map(isour)%nrcv > 0 ) THEN
            CALL irecv_real( rcvbuf, SIZE( rcvbuf ), isour, ip, irhand )    !!
          END IF

          DO ik = 1, SIZE(array,2)
            DO ib = 1, snd_map(idest)%nsnd
              sndbuf(ib + snd_map(idest)%nsnd*(ik-1)) = array( snd_map(idest)%iloc(ib), ik )
            END DO
          END DO 

          IF( snd_map(idest)%nsnd > 0 ) THEN
            CALL isend_real( sndbuf(1), SIZE( sndbuf ), idest, ip, ishand )  !!
          END IF
        
!          CALL sendrecv_real(sndbuf(1), SIZE(sndbuf), idest,     &      
!                             rcvbuf(1), SIZE(rcvbuf), isour, ip)
        
          IF( rcv_map(isour)%nrcv > 0 ) THEN
            CALL mp_wait( irhand )                    !!
          END IF

          DO ik = 1, SIZE(array,2)
            DO ib = 1, rcv_map(isour)%nrcv
              array( rcv_map(isour)%iloc(ib), ik ) = rcvbuf( ib + rcv_map(isour)%nrcv*(ik-1))
            END DO
          END DO 
          
          IF( snd_map(idest)%nsnd > 0 ) THEN
            CALL mp_wait( ishand )   !!
          END IF

          DEALLOCATE( rcvbuf, STAT=ierr )
          IF( ierr /= 0 ) &
            CALL error(' data_exchange_rm ', ' deallocating rcvbuf ', ierr )
          DEALLOCATE( sndbuf, STAT=ierr )
          IF( ierr /= 0 ) &
            CALL error(' data_exchange_rm ', ' deallocating sndbuf ', ierr )

        END DO
        RETURN
      END SUBROUTINE data_exchange_rm

!----------------------------------------------------------------------

      SUBROUTINE data_exchange_i(array)
        USE parallel, ONLY: nproc, mpime
        IMPLICIT NONE
        INTEGER :: array(:)
        INTEGER, ALLOCATABLE :: sndbuf(:), rcvbuf(:)
        INTEGER :: ip, isour, idest, ib, ierr
        DO ip = 1, (nproc - 1)
          isour = MOD(mpime - ip + nproc, nproc)
          idest = MOD(mpime + ip        , nproc)
          ALLOCATE( rcvbuf( MAX(rcv_map(isour)%nrcv,1) ), STAT=ierr )
          IF( ierr /= 0 ) THEN
            IF (lpr > 1) WRITE(testunit,*) 'Trying to allocate ', MAX(rcv_map(isour)%nrcv,1), ' elements '
            CALL error(' data_exchange_i ', ' allocating rcvbuf ', ierr )
          END IF
          ALLOCATE( sndbuf( MAX(snd_map(idest)%nsnd,1) ), STAT=ierr )
          IF( ierr /= 0 ) THEN
            IF (lpr > 1) WRITE(testunit,*) 'Trying to allocate ', MAX(snd_map(idest)%nsnd,1), ' elements '
            CALL error(' data_exchange_i ', ' allocating sndbuf ', ierr )
          END IF
          DO ib = 1, snd_map(idest)%nsnd
            sndbuf(ib) = array( snd_map(idest)%iloc(ib) )
          END DO
          CALL sendrecv_integer(sndbuf, snd_map(idest)%nsnd, idest,      &
                                rcvbuf, rcv_map(isour)%nrcv, isour, ip)
          DO ib = 1, rcv_map(isour)%nrcv
            array( rcv_map(isour)%iloc(ib) ) = rcvbuf(ib)
          END DO
          DEALLOCATE( rcvbuf, STAT=ierr )
          IF( ierr /= 0 ) &
            CALL error(' data_exchange_i ', ' deallocating rcvbuf ', ierr )
          DEALLOCATE( sndbuf, STAT=ierr )
          IF( ierr /= 0 ) &
            CALL error(' data_exchange_i ', ' deallocating sndbuf ', ierr )
        END DO
        RETURN
      END SUBROUTINE data_exchange_i

!----------------------------------------------------------------------

      SUBROUTINE data_exchange_l(array)
        USE parallel, ONLY: nproc, mpime
        IMPLICIT NONE
        LOGICAL :: array(:)
        INTEGER, ALLOCATABLE :: sndbuf(:), rcvbuf(:)
        INTEGER :: ip, isour, idest, ib, ierr
        DO ip = 1, (nproc - 1)
          isour = MOD(mpime - ip + nproc, nproc)
          idest = MOD(mpime + ip        , nproc)
          ALLOCATE( rcvbuf( MAX(rcv_map(isour)%nrcv,1) ), STAT=ierr )
          IF( ierr /= 0 ) THEN
            IF (lpr > 1) WRITE(testunit,*) 'Trying to allocate ', MAX(rcv_map(isour)%nrcv,1), ' elements '
            CALL error(' data_exchange_l ', ' allocating rcvbuf ', ierr )
          END IF
          ALLOCATE( sndbuf( MAX(snd_map(idest)%nsnd,1) ), STAT=ierr )
          IF( ierr /= 0 ) THEN
            IF (lpr > 1) WRITE(testunit,*) 'Trying to allocate ', MAX(snd_map(idest)%nsnd,1), ' elements '
            CALL error(' data_exchange_l ', ' allocating sndbuf ', ierr )
          END IF
          DO ib = 1, snd_map(idest)%nsnd
            IF( array( snd_map(idest)%iloc(ib) ) ) THEN
              sndbuf(ib) = 1
            ELSE
              sndbuf(ib) = 0
            END IF
          END DO
          CALL sendrecv_integer(sndbuf, snd_map(idest)%nsnd, idest,      &
                                rcvbuf, rcv_map(isour)%nrcv, isour, ip)
          DO ib = 1, rcv_map(isour)%nrcv
            IF( rcvbuf(ib) == 1 ) THEN
              array( rcv_map(isour)%iloc(ib) ) = .TRUE.
            ELSE
              array( rcv_map(isour)%iloc(ib) ) = .FALSE.
            END IF
          END DO
          DEALLOCATE( rcvbuf, STAT=ierr )
          IF( ierr /= 0 ) &
            CALL error(' data_exchange_l ', ' deallocating rcvbuf ', ierr )
          DEALLOCATE( sndbuf, STAT=ierr )
          IF( ierr /= 0 ) &
            CALL error(' data_exchange_l ', ' deallocating sndbuf ', ierr )
        END DO
        RETURN
      END SUBROUTINE data_exchange_l

!----------------------------------------------------------------------

      SUBROUTINE data_collect_r( garray, larray, imstart, imend )

        ! this subroutine collect data distributed across processors
        ! and store them in the array "garray"
        ! The subroutine is designed to allow the collection of global
        ! data using subarray of small size, useful when available memory
        ! is not enough to store a global array

        USE parallel, ONLY: nproc, mpime
        USE indijk_module, ONLY: ip0_jp0_kp0_

        IMPLICIT NONE
        REAL*8 :: garray(:) ! global array that is set with collected data
                            ! garray could be only a subarray of the whole
                            ! global array data being collected refers to.
                            ! Its first element has global index imstart
                            ! therefore its size should be = imend - imstart + 1
        REAL*8 :: larray(:) ! local array
        INTEGER :: imstart  ! global index from which we start to collect
        INTEGER :: imend    ! global index at which the collection is stopped

        INTEGER :: ijk, imesh, ncollect

        ncollect = ( imend - imstart + 1 ) 
        IF( SIZE( garray ) < ncollect ) &
          CALL error(' data_collect_r ', ' garray too small ', SIZE( garray ) )

        garray( 1 : ncollect )  = 0.0

        DO ijk = 1, ncint
          imesh = myijk( ip0_jp0_kp0_, ijk)
          IF( imesh >= imstart .AND. imesh <= imend ) THEN
            garray( imesh - imstart + 1 ) = larray( ijk )
          END IF
        END DO

        CALL parallel_sum_real( garray, ncollect )
        
        RETURN
      END SUBROUTINE data_collect_r

      SUBROUTINE data_collect_sr( garray, larray, imstart, imend )

        ! this subroutine collect data distributed across processors
        ! and store them in the array "garray"
        ! The subroutine is designed to allow the collection of global
        ! data using subarray of small size, useful when available memory
        ! is not enough to store a global array

        USE parallel, ONLY: nproc, mpime
        USE indijk_module, ONLY: ip0_jp0_kp0_

        IMPLICIT NONE
        REAL(sgl) :: garray(:) ! global array that is set with collected data
                            ! garray could be only a subarray of the whole
                            ! global array data being collected refers to.
                            ! Its first element has global index imstart
                            ! therefore its size should be = imend - imstart + 1
        REAL*8 :: larray(:) ! local array
        INTEGER :: imstart  ! global index from which we start to collect
        INTEGER :: imend    ! global index at which the collection is stopped

        INTEGER :: ijk, imesh, ncollect

        ncollect = ( imend - imstart + 1 )
        IF( SIZE( garray ) < ncollect ) &
          CALL error(' data_collect_r ', ' garray too small ', SIZE( garray ) )

        garray( 1 : ncollect )  = 0.0

        DO ijk = 1, ncint
          imesh = myijk( ip0_jp0_kp0_, ijk)
          IF( imesh >= imstart .AND. imesh <= imend ) THEN
            garray( imesh - imstart + 1 ) = larray( ijk )
          END IF
        END DO

        CALL parallel_sum_sreal( garray, ncollect )

        RETURN
      END SUBROUTINE data_collect_sr

      SUBROUTINE data_distribute_r( garray, larray, imstart, imend )

        ! this subroutine distribute a global array "garray"  across processors
        ! and store them in the local array "larray"
        ! The subroutine is designed to allow the distribution of global
        ! data using subarray of small size, useful when available memory
        ! is not enough to store a global array

        USE parallel, ONLY: nproc, mpime
        USE indijk_module, ONLY: ip0_jp0_kp0_

        IMPLICIT NONE
        REAL*8 :: garray(:) ! input global array ( on root ) with data to be distributed
                            ! garray is a subarray of the whole
                            ! global array being distributed.
                            ! Its first element has global index imstart
                            ! therefore its size should be = imend - imstart + 1
        REAL*8 :: larray(:) ! output local array, its elements are set 
                            ! with the global data belonging to the local processor
        INTEGER :: imstart  ! global index from which we start to collect
        INTEGER :: imend    ! global index at which the collection is stopped

        INTEGER :: ijk, imesh, ndistribute

        ndistribute = ( imend - imstart + 1 )
        IF( SIZE( garray ) < ndistribute ) &
          CALL error(' data_distribute_r ', ' garray too small ', SIZE( garray ) )

        CALL bcast_real( garray, ndistribute, 0 )

        DO ijk = 1, ncint
          imesh = myijk( ip0_jp0_kp0_, ijk)
          IF( imesh >= imstart .AND. imesh <= imend ) THEN
            larray( ijk ) = garray( imesh - imstart + 1 ) 
          END IF
        END DO

        RETURN
      END SUBROUTINE data_distribute_r

      SUBROUTINE data_distribute_sr( garray, larray, imstart, imend )

        ! this subroutine distribute a global array "garray"  across processors
        ! and store them in the local array "larray"
        ! The subroutine is designed to allow the distribution of global
        ! data using subarray of small size, useful when available memory
        ! is not enough to store a global array

        USE parallel, ONLY: nproc, mpime
        USE indijk_module, ONLY: ip0_jp0_kp0_

        IMPLICIT NONE
        REAL(sgl) :: garray(:) ! input global array ( on root ) with data to be distributed
                            ! garray is a subarray of the whole
                            ! global array being distributed.
                            ! Its first element has global index imstart
                            ! therefore its size should be = imend - imstart + 1
        REAL*8 :: larray(:) ! output local array, its elements are set
                            ! with the global data belonging to the local processor
        INTEGER :: imstart  ! global index from which we start to collect
        INTEGER :: imend    ! global index at which the collection is stopped

        INTEGER :: ijk, imesh, ndistribute

        ndistribute = ( imend - imstart + 1 )
        IF( SIZE( garray ) < ndistribute ) &
          CALL error(' data_distribute_r ', ' garray too small ', SIZE( garray ) )

        CALL bcast_sreal( garray, ndistribute, 0 )

        DO ijk = 1, ncint
          imesh = myijk( ip0_jp0_kp0_, ijk)
          IF( imesh >= imstart .AND. imesh <= imend ) THEN
            larray( ijk ) = garray( imesh - imstart + 1 )
          END IF
        END DO

        RETURN
      END SUBROUTINE data_distribute_sr
!----------------------------------------------------------------------
      SUBROUTINE test_comm()
        IMPLICIT NONE
        INTEGER, ALLOCATABLE :: itest(:)
        INTEGER :: ijl
        
        ALLOCATE( itest(ncdom) )

        itest = 0
        DO ijl = 1, ncint
          itest( ijl ) = myijk( ip0_jp0_kp0_, ijl)
        END DO

        CALL data_exchange(itest)
         
        WRITE(testunit,*)   ' index received ' 
        WRITE(testunit,310) itest( ncint + 1 : ncdom ) 

 310    FORMAT(10i8)
        
        DEALLOCATE( itest )
        
        RETURN
      END SUBROUTINE test_comm
!----------------------------------------------------------------------
      SUBROUTINE local_forcing
      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions, ONLY: nx, ny, nz
      USE grid, ONLY: fc2
      USE immersed_boundaries, ONLY: numx, fptx, numy, fpty, numz, fptz
      USE parallel, ONLY: mpime, root

      IMPLICIT NONE
      INTEGER :: n, i, j, k, ijk, ijkl
      INTEGER :: delta_i, delta_j, delta_k, ijk_q, ijk_qq
      INTEGER :: nfpx, nfpy, nfpz
!
! ... tag for filled_cells_2
!
        ALLOCATE( fc2(ncdom) ); fc2 = .FALSE.
!
        ALLOCATE( numx(ncdom) ); numx = 0
        nfpx = SIZE(fptx)

        IF (job_type == JOB_TYPE_3D) THEN
          ALLOCATE( numy(ncdom) ); numy = 0
          nfpy = SIZE(fpty)
        ELSE
          nfpy = 0.D0
        END IF

        ALLOCATE( numz(ncdom) ); numz = 0
        nfpz = SIZE(fptz)

        DO n = 1, nfpx
          i = fptx(n)%i
          j = fptx(n)%j
          k = fptx(n)%k
          delta_i = fptx(n)%delta_i
          delta_j = fptx(n)%delta_j
          delta_k = fptx(n)%delta_k
set_numx: IF (i/=0 .AND. k/=0) THEN
            IF (job_type == JOB_TYPE_2D) THEN
              ijk = i + (k-1) * nx
            ELSE IF (job_type == JOB_TYPE_3D) THEN
              IF (j == 0) CALL error('decomp','control numx',1)
              ijk = i + (j-1) * nx + (k-1) * nx * ny
              ijk_q = (i+delta_i) + (j+delta_j-1) * nx + (k+delta_k-1) * nx * ny
              ijk_qq = (i+2*delta_i) + (j+2*delta_j-1) * nx + (k+2*delta_k-1) * nx * ny
            END IF
            IF( cell_owner(ijk) == mpime ) THEN

              ijkl = cell_g2l(ijk,mpime)
              numx(ijkl) = n 
              !
              ! ... Identify the neighbours for imm.b. interpolations
              ! ... Second neighbours 'ijk_qq' are not required along
              ! ... diagonals
              fptx(n)%index_q  = cell_g2l(ijk_q,mpime)
!IF( fptx(n)%index_q < 0 ) write(errorunit,*) 'Eccolo x= ',fptx(n)%index_q, n, ijkl, ijk, mpime, ijk_q
              IF (MOD(fptx(n)%int,2)/=0 .OR. fptx(n)%int==0) & 
                fptx(n)%index_qq = cell_g2l(ijk_qq,mpime)

              IF (k==1 .OR. k==nz .OR. i==1 .OR. i==nx .OR. j==1 .OR. j==ny) THEN
                IF (mpime == root) THEN
                  WRITE(errorunit,*) 'WARNING! from ghost'
                  WRITE(errorunit,*) 'skipping x-forcing on boundaries', i, j, k
                END IF
              END IF
            END IF
          ELSE
            CALL error('decomp','control numx',1)
          END IF set_numx

        END DO
        !
        DO n = 1, nfpy
        
          IF (job_type == JOB_TYPE_3D) THEN
            i = fpty(n)%i
            j = fpty(n)%j
            k = fpty(n)%k
            delta_i = fpty(n)%delta_i
            delta_j = fpty(n)%delta_j
            delta_k = fpty(n)%delta_k
set_numy:   IF (i/=0 .AND. k/=0) THEN
              IF (job_type == JOB_TYPE_2D) THEN
                ijk = i + (k-1) * nx
              ELSE IF (job_type == JOB_TYPE_3D) THEN
                IF (j == 0) CALL error('decomp','control numy',1)
                ijk = i + (j-1) * nx + (k-1) * nx * ny
                ijk_q = (i+delta_i) + (j+delta_j-1) * nx + (k+delta_k-1) * nx * ny
                ijk_qq = (i+2*delta_i) + (j+2*delta_j-1) * nx + (k+2*delta_k-1) * nx * ny
              END IF
              IF( cell_owner(ijk)== mpime ) THEN

                ijkl = cell_g2l(ijk,mpime)
                numy(ijkl) = n 
                !
                ! ... Identify the neighbours for imm.b. interpolations
                ! ... Second neighbours 'ijk_qq' are not required along
                ! ... diagonals
                fpty(n)%index_q  = cell_g2l(ijk_q,mpime)
!IF( fpty(n)%index_q < 0 ) write(errorunit,*) 'Eccolo y= ',fpty(n)%index_q, n, ijkl, ijk, mpime, ijk_q
                IF (MOD(fpty(n)%int,2)/=0 .OR. fpty(n)%int==0) & 
                  fpty(n)%index_qq = cell_g2l(ijk_qq,mpime)

                IF (k==1 .OR. k==nz .OR. i==1 .OR. i==nx .OR. j==1 .OR. j==ny) THEN
                  IF (mpime == root) THEN
                    WRITE(errorunit,*) 'WARNING! from ghost'
                    WRITE(errorunit,*) 'skipping y-forcing on boundaries', i, j, k
                  END IF
                END IF
              END IF
            ELSE
              CALL error('decomp','control numy',1)
            END IF set_numy
          END IF

        END DO
        !
        DO n = 1, nfpz
        
          i = fptz(n)%i
          j = fptz(n)%j
          k = fptz(n)%k
          delta_i = fptz(n)%delta_i
          delta_j = fptz(n)%delta_j
          delta_k = fptz(n)%delta_k
set_numz: IF (i/=0 .AND. k/=0) THEN
            IF (job_type == JOB_TYPE_2D) THEN
              ijk = i + (k-1) * nx
            ELSE IF (job_type == JOB_TYPE_3D) THEN
              IF (j == 0) CALL error('decomp','control numz',n)
              ijk = i + (j-1) * nx + (k-1) * nx * ny
              ijk_q = (i+delta_i) + (j+delta_j-1) * nx + (k+delta_k-1) * nx * ny
              ijk_qq = (i+2*delta_i) + (j+2*delta_j-1) * nx + (k+2*delta_k-1) * nx * ny
            END IF
            IF( cell_owner(ijk)== mpime ) THEN

              ijkl = cell_g2l(ijk,mpime)
              numz(ijkl) = n 
              !
              ! ... Identify the neighbours for imm.b. interpolations
              ! ... Second neighbours 'ijk_qq' are not required along
              ! ... diagonals
              fptz(n)%index_q  = cell_g2l(ijk_q,mpime)
!IF( fptz(n)%index_q < 0 ) write(errorunit,*) 'Eccolo z= ',fptz(n)%index_q, n, ijkl, ijk, mpime, ijk_q
              IF (MOD(fptz(n)%int,2)/=0 .OR. fptz(n)%int==0) & 
                fptz(n)%index_qq = cell_g2l(ijk_qq,mpime)

              IF (k==1 .OR. k==nz .OR. i==1 .OR. i==nx .OR. j==1 .OR. j==ny) THEN
                IF (mpime == root) THEN
                  WRITE(errorunit,*) 'WARNING! from ghost'
                  WRITE(errorunit,*) 'skipping z-forcing on boundaries', i, j, k
                END IF
              END IF
            END IF
          ELSE
            CALL error('decomp','control numz',n)
          END IF set_numz
        
        END DO

        CALL data_exchange(numx)
        IF (job_type == JOB_TYPE_3D) CALL data_exchange(numy)
        CALL data_exchange(numz)

        CALL fill_cells

      RETURN
      END SUBROUTINE local_forcing
!----------------------------------------------------------------------
      SUBROUTINE fill_cells
!
! ... Coefficients b(1:6) are multiplied to the numerical mass
! ... fluxes to take into account that some cell face could be
! ... immersed. b=1 if the cell face is external, b=0 if the
! ... cell face is inside the topography.
!
      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE grid, ONLY: z, zb, flag, dz
      USE immersed_boundaries, ONLY: bd, vf
      USE io_files, ONLY: tempunit
      USE parallel, ONLY: mpime, root
      USE volcano_topography, ONLY: topo_c, topo_x
      USE volcano_topography, ONLY: topo2d_c, topo2d_x, topo2d_y
      USE indijk_module
      IMPLICIT NONE

      INTEGER :: i,j,k,ijk,imesh
      INTEGER :: ipjk, imjk, ijpk, ijmk, ijkp, ijkm
      INTEGER :: full, counter
!
! ... Allocate and initialize coefficients
!
      ALLOCATE(bd(ncint))
      ALLOCATE(vf(ncint))

      IF (job_type == JOB_TYPE_2D) THEN
        full = 15
      ELSE IF (job_type == JOB_TYPE_3D) THEN
        full = 63
      END IF

      bd(:) = 0
      vf(:) = 0.D0

      IF( lpr > 1 ) THEN
        WRITE( testunit, * ) 
        WRITE( testunit, * ) 'b coefficients for immersed bndaries' 
        WRITE( testunit, * ) '     ijk   i   j   k     b (int)'
      END IF
!
        DO ijk=1, ncint

          imjk  = myijk( im1_jp0_kp0_, ijk )
          ipjk  = myijk( ip1_jp0_kp0_, ijk )
          ijmk  = myijk( ip0_jm1_kp0_, ijk )
          ijpk  = myijk( ip0_jp1_kp0_, ijk )
          ijkm  = myijk( ip0_jp0_km1_, ijk )
          ijkp  = myijk( ip0_jp0_kp1_, ijk )

          IF( BTEST(flag(ijk),0) ) THEN
            CALL meshinds(ijk,imesh,i,j,k)
  
            IF (job_type == JOB_TYPE_2D) THEN
              counter = 0
              ! 
              ! East
!              IF (z(k) > topo_x(i) .AND. flag(ipjk) /= filled_cell_1) THEN
               IF (z(k) > topo_x(i)) THEN
                  IF (flag(ipjk) == filled_cell_1) THEN
                     vf(ijk) = vf(ijk) + 2.D0
                  ELSE
                     bd(ijk) = bd(ijk) + 1
                     vf(ijk) = vf(ijk) + 1.D0
                     counter = counter + 1
                  END IF
               END IF
              !
              ! West
!              IF (z(k) > topo_x(i-1) .AND. flag(imjk) /= filled_cell_1) THEN
               IF (z(k) > topo_x(i-1)) THEN
                  IF (flag(imjk) == filled_cell_1) THEN
                     vf(ijk) = vf(ijk) + 2.D0
                  ELSE
                     bd(ijk) = bd(ijk) + 2
                     vf(ijk) = vf(ijk) + 1.D0
                     counter = counter + 1
                  END IF
               END IF
              !  
              ! Top
               IF (zb(k) > topo_c(i)) THEN
                  bd(ijk) = bd(ijk) + 4 
                  vf(ijk) = vf(ijk) + 1.D0
                  counter = counter + 1
               END IF
              !
              ! Bottom
!              IF (zb(k-1) >= topo_c(i) .AND. flag(ijkm) /= filled_cell_1) THEN
               IF (zb(k-1) >= topo_c(i)) THEN
                  IF (flag(ijkm) == filled_cell_1) THEN
                     vf(ijk) = vf(ijk) + 2.D0   
                  ELSE
                     bd(ijk) = bd(ijk) + 8  
                     vf(ijk) = vf(ijk) + 1.D0
                     counter = counter + 1
                  END IF
               END IF
              !
            ELSE IF (job_type == JOB_TYPE_3D) THEN
              counter = 0
              ! 
              ! East
!              IF (z(k) > topo2d_x(i,j) .AND. flag(ipjk) /= filled_cell_1) THEN
              IF (z(k) > topo2d_x(i,j)) THEN
                 IF (flag(ipjk) == filled_cell_1) THEN
                    vf(ijk) = vf(ijk) + 2.D0
                 ELSE IF (flag(ipjk) == filled_cell_2) THEN
                    vf(ijk) = vf(ijk) + 1.D0
                 ELSE
                    bd(ijk) = bd(ijk) + 1
                    vf(ijk) = vf(ijk) + 1.D0
                    counter = counter + 1
                 END IF
              END IF
              !
              ! West
!              IF (z(k) > topo2d_x(i-1,j) .AND. flag(imjk) /= filled_cell_1) THEN
              IF (z(k) > topo2d_x(i-1,j)) THEN
                 IF (flag(imjk) == filled_cell_1) THEN
                    vf(ijk) = vf(ijk) + 2.D0  
                 ELSE IF (flag(imjk) == filled_cell_2) THEN
                    vf(ijk) = vf(ijk) + 1.D0  
                 ELSE
                    bd(ijk) = bd(ijk) + 2
                    vf(ijk) = vf(ijk) + 1.D0
                    counter = counter + 1
                 END IF
              END IF
              !
              ! Top
              IF (zb(k) > topo2d_c(i,j))  THEN
                bd(ijk) = bd(ijk) + 4 
                vf(ijk) = vf(ijk) + 1.D0
                counter = counter + 1
              END IF
              !
              ! Bottom
!              IF (zb(k-1) >= topo2d_c(i,j) .AND. flag(ijkm) /= filled_cell_1) THEN
              IF (zb(k-1) >= topo2d_c(i,j)) THEN
                IF (flag(ijkm) == filled_cell_1) THEN
                   vf(ijk) = vf(ijk) + 2.D0                   
                ELSE IF (flag(ijkm) == filled_cell_2) THEN
                   vf(ijk) = vf(ijk) + 1.D0                   
                ELSE
                   bd(ijk) = bd(ijk) + 8
                   vf(ijk) = vf(ijk) + 1.D0
                   counter = counter + 1
                END IF
             END IF
              !
              ! North
!              IF (z(k) > topo2d_y(i,j) .AND. flag(ijpk) /= filled_cell_1) THEN
              IF (z(k) > topo2d_y(i,j)) THEN
                IF (flag(ijpk) == filled_cell_1) THEN
                   vf(ijk) = vf(ijk) + 2.D0
                ELSE IF (flag(ijpk) == filled_cell_2) THEN
                   vf(ijk) = vf(ijk) + 1.D0
                ELSE
                   bd(ijk) = bd(ijk) + 16
                   vf(ijk) = vf(ijk) + 1.D0
                   counter = counter + 1
                END IF
             END IF
              !
              ! South
!              IF (z(k) > topo2d_y(i,j-1) .AND. flag(ijmk) /= filled_cell_1) THEN
              IF (z(k) > topo2d_y(i,j-1)) THEN
                 IF (flag(ijmk) == filled_cell_1) THEN
                   vf(ijk) = vf(ijk) + 2.D0
                 ELSE IF (flag(ijmk) == filled_cell_2) THEN
                   vf(ijk) = vf(ijk) + 1.D0
                 ELSE
                   bd(ijk) = bd(ijk) + 32
                   vf(ijk) = vf(ijk) + 1.D0
                   counter = counter + 1
                 END IF
              END IF
              !
            END IF

            IF (job_type == JOB_TYPE_2D) THEN
              vf(ijk) = 0.25D0 * vf(ijk)
            ELSE IF( job_type == JOB_TYPE_3D) THEN
              vf(ijk) = vf(ijk) / 6.D0
            END IF

            ! ... cells that are completely immersed
            ! ... are excluded from computation. 
            !
            IF (vf(ijk) == 0.D0) flag(ijk) = filled_cell_1
            !
            IF (lpr > 1) THEN
              IF (counter == 1 .AND. flag(ijk) /= filled_cell_1) &
                WRITE(testunit,*) 'WARNING non-filled cell', ijk, i, j, k
              IF (counter == 2 .AND. flag(ijk) /= filled_cell_2) &
                WRITE( testunit, fmt = "('NF2:')" )
            END IF

          END IF
!
          IF (lpr > 1) THEN
            IF( BTEST(flag(ijk),0) .OR. BTEST(flag(ijk),9) .OR. BTEST(flag(ijk),10)) THEN
              CALL meshinds(ijk,imesh,i,j,k)
              IF (bd(ijk) /= full) &
                WRITE( testunit, fmt = "( I8,3I4,2X,B8,I5 )" ) ijk, i, j, k, bd(ijk), flag(ijk)
            END IF
          END IF

        END DO
!
! ... Exchange modified flags
!
        CALL data_exchange(flag)
      
      RETURN
      END SUBROUTINE fill_cells
!----------------------------------------------------------------------
   END MODULE domain_mapping
!----------------------------------------------------------------------
