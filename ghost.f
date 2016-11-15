!----------------------------------------------------------------------
   MODULE domain_mapping
!----------------------------------------------------------------------
        USE domain_decomposition, ONLY: nctot, proc_map
        USE domain_decomposition, ONLY: layer_map, block2d_map, block3d_map
        USE domain_decomposition, ONLY: cell_g2l, cell_l2g, cell_owner
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
        PUBLIC :: ghost, meshinds, set_myinds
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
      USE parallel, ONLY: mpime, nproc
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
      CALL face_init
      IF (immb == 1 .AND. prog == 'PDAC') THEN
        CALL local_forcing
        CALL fill_cells
        CALL face_fractions
      END IF
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
          DO im = -2, 3
            DO km = -2, 3
              IF( ( ABS( im ) + ABS( km ) ) <= 4 ) THEN
                IF( (im /= 0) .OR. (km /= 0) ) THEN
                  ii = im
                  kk = km 
                  IF( ( i == 2    ) .AND. ( ii == -2 ) ) ii = -1
                  IF( ( i == nx-1 ) .AND. ( ii == +2 ) ) ii = +1
                  IF( ( i == nx-1 ) .AND. ( ii == +3 ) ) ii = +1
                  IF( ( i == nx-2 ) .AND. ( ii == +3 ) ) ii = +2
                  IF( ( k == 2    ) .AND. ( kk == -2 ) ) kk = -1
                  IF( ( k == nz-1 ) .AND. ( kk == +2 ) ) kk = +1
                  IF( ( k == nz-1 ) .AND. ( kk == +3 ) ) kk = +1
                  IF( ( k == nz-2 ) .AND. ( kk == +3 ) ) kk = +2
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
               DO km = -2, 3
                 kk = km
                 IF( ( k == 2    ) .AND. ( kk == -2 ) ) kk = -1
                 IF( ( k == nz-1 ) .AND. ( kk == +2 ) ) kk = +1
                 IF( ( k == nz-1 ) .AND. ( kk == +3 ) ) kk = +1
                 IF( ( k == nz-2 ) .AND. ( kk == +3 ) ) kk = +2
                 DO jm = -2, 3
                   jj = jm
                   IF( ( j == 2    ) .AND. ( jj == -2 ) ) jj = -1
                   IF( ( j == ny-1 ) .AND. ( jj == +2 ) ) jj = +1
                   IF( ( j == ny-1 ) .AND. ( jj == +3 ) ) jj = +1
                   IF( ( j == ny-2 ) .AND. ( jj == +3 ) ) jj = +2
                   DO im = -2, 3
                     IF( ( ABS( im ) + ABS( jm ) + ABS( km ) ) <= 5 ) THEN
                       IF( (im /= 0) .OR. (jm /= 0) .OR. (km /= 0) ) THEN
                         ii = im
                         IF( ( i == 2    ) .AND. ( ii == -2 ) ) ii = -1
                         IF( ( i == nx-1 ) .AND. ( ii == +2 ) ) ii = +1
                         IF( ( i == nx-1 ) .AND. ( ii == +3 ) ) ii = +1
                         IF( ( i == nx-2 ) .AND. ( ii == +3 ) ) ii = +2
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

             DO km = -2, 3
               kk = km
               IF( ( k == 2    ) .AND. ( kk == -2 ) ) kk = -1
               IF( ( k == nz-1 ) .AND. ( kk == +2 ) ) kk = +1
               IF( ( k == nz-1 ) .AND. ( kk == +3 ) ) kk = +1
               IF( ( k == nz-2 ) .AND. ( kk == +3 ) ) kk = +2
               DO jm = -2, 3
                 jj = jm
                 IF( ( j == 2    ) .AND. ( jj == -2 ) ) jj = -1
                 IF( ( j == ny-1 ) .AND. ( jj == +2 ) ) jj = +1
                 IF( ( j == ny-1 ) .AND. ( jj == +3 ) ) jj = +1
                 IF( ( j == ny-2 ) .AND. ( jj == +3 ) ) jj = +2
                 DO im = -2, 3
                   IF( ( ABS( im ) + ABS( jm ) + ABS( km ) ) <= 5 ) THEN
                     IF( (im /= 0) .OR. (jm /= 0) .OR. (km /= 0) ) THEN
                       ii = im
                       IF( ( i == 2    ) .AND. ( ii == -2 ) ) ii = -1
                       IF( ( i == nx-1 ) .AND. ( ii == +2 ) ) ii = +1
                       IF( ( i == nx-1 ) .AND. ( ii == +3 ) ) ii = +1
                       IF( ( i == nx-2 ) .AND. ( ii == +3 ) ) ii = +2
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
        USE flux_limiters, ONLY: ctu
!
        IMPLICIT NONE
        INTEGER :: myinds(:,:)
        INTEGER :: myijk(:,:)
!
        INTEGER :: i, j, k, ijk, imesh

        INTEGER :: ipjk, imjk, ijpk, ijmk, ijkp, ijkm, &
                   ippjk, immjk, ijppk, ijmmk, ijkpp, ijkmm, &
                   ipjpk, ipjmk, ipjkp, ipjkm, & 
                   imjpk, imjmk, imjkp, imjkm, &
                   ijpkp, ijpkm, ijmkp, ijmkm, &
                   ipjpkp, ipjpkm, ipjmkp, ipjmkm, imjpkp, imjpkm, imjmkp, imjmkm, &
                   ippjkp, ippjkm, ippjpk, ippjmk, immjkp, immjkm, immjpk, immjmk, &
                   ipjppk, imjppk, ijppkp, ijppkm, ipjmmk, imjmmk, ijmmkp, ijmmkm, &
                   ipjkpp, imjkpp, ijpkpp, ijmkpp, ipjkmm, imjkmm, ijpkmm, ijmkmm, &
                   ippjpkp, ippjpkm, ippjmkp, ippjmkm, immjpkp, immjpkm, immjmkp, immjmkm, &
                   ipjppkp, ipjppkm, imjppkp, imjppkm, ipjmmkp, ipjmmkm, imjmmkp, imjmmkm, &
                   ipjpkpp, ipjmkpp, imjpkpp, imjmkpp, ipjpkmm, ipjmkmm, imjpkmm, imjmkmm, &
                   ippjppk, ippjmmk, ippjkpp, &
                   ijppkpp, ijppkmm, immjppk, &
                   immjkpp, ijmmkpp, ippjkmm, &
                   ipppjk, ipppjpk, ipppjmk, ipppjkp, ipppjpkp, ipppjmkp, &
                   ijpppk, ipjpppk, ijpppkp, ijpppkm, ipjpppkp, ipjpppkm, &
                   ijkppp, ipjkppp, imjkppp, ijpkppp, ipjpkppp, imjpkppp, &
                   ippjppkp, ippjppkm, ippjmmkp, ippjmmkm, ippjpkmm, ippjmkmm, &
                   ipjppkpp, imjppkpp, ipjppkmm, imjppkmm, immjppkp, immjppkm, &
                   ippjpkpp, ippjmkpp, ipjmmkpp, immjpkpp, immjmkpp, imjmmkpp
                 

        INTEGER :: ijke, ijkw, ijkn, ijks, ijkt, ijkb, &
                   ijkee, ijkww, ijknn, ijkss, ijktt, ijkbb, &
                   ijken, ijkes, ijket, ijkeb, &
                   ijkwn, ijkws, ijkwt, ijkwb, &
                   ijknt, ijknb, ijkst, ijksb, &
                   ijkent, ijkenb, ijkest, ijkesb, ijkwnt, ijkwnb, ijkwst, ijkwsb, &
                   ijkeet, ijkeeb, ijkeen, ijkees, ijkwwt, ijkwwb, ijkwwn, ijkwws, &
                   ijkenn, ijkwnn, ijknnt, ijknnb, ijkess, ijkwss, ijksst, ijkssb, &
                   ijkett, ijkwtt, ijkntt, ijkstt, ijkebb, ijkwbb, ijknbb, ijksbb, &
                   ijkeent, ijkeenb, ijkeest, ijkeesb, ijkwwnt, ijkwwnb, ijkwwst, ijkwwsb, &
                   ijkennt, ijkennb, ijkwnnt, ijkwnnb, ijkesst, ijkessb, ijkwsst, ijkwssb, &
                   ijkentt, ijkestt, ijkwntt, ijkwstt, ijkenbb, ijkesbb, ijkwnbb, ijkwsbb, &
                   ijkeenn, ijkeess, ijkeett, &
                   ijknntt, ijknnbb, ijkwwnn, &
                   ijkwwtt, ijksstt, ijkeebb, &
                   ijkeee, ijkeeen, ijkeees, ijkeeet, ijkeeent, ijkeeest, &
                   ijknnn, ijkennn, ijknnnt, ijknnnb, ijkennnt, ijkennnb, &
                   ijkttt, ijkettt, ijkwttt, ijknttt, ijkenttt, ijkwnttt, &
                   ijkeennt, ijkeennb, ijkeesst, ijkeessb, ijkeenbb, ijkeesbb, &
                   ijkenntt, ijkwnntt, ijkennbb, ijkwnnbb, ijkwwnnt, ijkwwnnb, &
                   ijkeentt, ijkeestt, ijkesstt, ijkwwntt, ijkwwstt, ijkwsstt 

        DO ijk = 1, ncint 
          CALL meshinds(ijk,imesh,i,j,k)

          IF( job_type == JOB_TYPE_2D ) THEN

            IF( (i >= 2) .AND. (i <= (nx-1)) .AND.   &
                (k >= 2) .AND. (k <= (nz-1))      ) THEN
!
              ijkm    = myijk( ip0_jp0_km1_, ijk)
              imjk    = myijk( im1_jp0_kp0_, ijk)
              ipjk    = myijk( ip1_jp0_kp0_, ijk)
              ijkp    = myijk( ip0_jp0_kp1_, ijk)
              ipjkm   = myijk( ip1_jp0_km1_, ijk)
              ipjkp   = myijk( ip1_jp0_kp1_, ijk)
              imjkm   = myijk( im1_jp0_km1_, ijk)
              imjkp   = myijk( im1_jp0_kp1_, ijk)
              ippjk   = myijk( ip2_jp0_kp0_, ijk)
              ijkpp   = myijk( ip0_jp0_kp2_, ijk)
              immjk   = myijk( im2_jp0_kp0_, ijk)
              ijkmm   = myijk( ip0_jp0_km2_, ijk)
              ippjkp  = myijk( ip2_jp0_kp1_, ijk)
              ippjkm  = myijk( ip2_jp0_km1_, ijk)
              ipjkpp  = myijk( ip1_jp0_kp2_, ijk)
              imjkpp  = myijk( im1_jp0_kp2_, ijk)
              immjkp  = myijk( im2_jp0_kp1_, ijk)
              immjkm  = myijk( im2_jp0_km1_, ijk)
              ipjkmm  = myijk( ip1_jp0_km2_, ijk)
              imjkmm  = myijk( im1_jp0_km2_, ijk)
              ipppjk  = myijk( ip3_jp0_kp0_, ijk)
              ijkppp  = myijk( ip0_jp0_kp3_, ijk)
              ipppjkp = myijk( ip3_jp0_kp1_, ijk) 
              ipjkppp = myijk( ip1_jp0_kp3_, ijk)
              ippjkpp = myijk( ip2_jp0_kp2_, ijk) 
              ippjkmm = myijk( ip2_jp0_km2_, ijk) 
              immjkpp = myijk( im2_jp0_kp2_, ijk) 
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
              IF(flag(imjkp) == slip_wall .OR. flag(imjkp) == noslip_wall) THEN
                IF(flag(ijkp) == slip_wall .OR. flag(ijkp) == noslip_wall) THEN
                  IF(flag(imjk) == slip_wall .OR. flag(imjk) == noslip_wall) THEN
                    ijkwt = ijk
                  ELSE
                    ijkwt = imjk
                  END IF
                ELSE
                  IF(flag(imjk) == slip_wall .OR. flag(imjk) == noslip_wall) THEN
                    ijkwt = ijkp
                  ELSE
                    ijkwt = ijk
                  END IF
                END IF
              END IF
  
              ijkeb = ipjkm
              IF(flag(ipjkm) == slip_wall .OR. flag(ipjkm) == noslip_wall) THEN
                IF (flag(ipjk) == slip_wall .OR. flag(ipjk) == noslip_wall) THEN
                  IF (flag(ijkm) == slip_wall .OR. flag(ijkm) == noslip_wall) THEN
                    ijkeb = ijk
                  ELSE
                    ijkeb = ijkm
                  END IF
                ELSE
                  IF (flag(ijkm) == slip_wall .OR. flag(ijkm) == noslip_wall) THEN
                    ijkeb = ipjk
                  ELSE
                    ijkeb = ijk
                  END IF
                END IF
              END IF

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
              IF (flag(imjkm) == slip_wall .OR. flag(imjkm) == noslip_wall) THEN
                IF (flag(imjk) == slip_wall .OR. flag(imjk) == noslip_wall) THEN
                  IF (flag(ijkm) == slip_wall .OR. flag(ijkm) == noslip_wall) THEN
                    ijkwb = ijk
                  ELSE
                    ijkwb = ijkm
                  END IF
                ELSE
                  IF (flag(ijkm) == slip_wall .OR. flag(ijkm) == noslip_wall) THEN
                    ijkwb = imjk
                  ELSE
                    ijkwb = ijk
                  END IF
                END IF
              END IF
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
              SELECT CASE (flag(ippjkp))
              CASE (slip_wall, noslip_wall, filled_cell_1)
                IF (flag(ipjkp) == slip_wall .OR. flag(ipjkp) == noslip_wall) THEN
                  IF (flag(ippjk) == slip_wall .OR. flag(ippjk) == noslip_wall) THEN
                    ijkeet = ipjk
                  ELSE
                    ijkeet = ippjk
                  END IF
                ELSE
                  IF (flag(ippjk) == slip_wall .OR. flag(ippjk) == noslip_wall) THEN
                    ijkeet = ipjkp
                  ELSE
                    ijkeet = ipjk
                  END IF
                END IF
              CASE DEFAULT
                      ijkeet = ippjkp
              END SELECT
              IF(i == (nx-1)) ijkeet = ijket
!
              SELECT CASE (flag(ippjkm))
              CASE (slip_wall, noslip_wall, filled_cell_1)
                IF (flag(ipjkm) == slip_wall .OR. flag(ipjkm) == noslip_wall) THEN
                  IF (flag(ippjk) == slip_wall .OR. flag(ippjk) == noslip_wall) THEN
                    ijkeeb = ipjk
                  ELSE
                    ijkeeb = ippjk
                  END IF
                ELSE
                  IF (flag(ippjk) == slip_wall .OR. flag(ippjk) == noslip_wall) THEN
                    ijkeeb = ipjkm
                  ELSE
                    ijkeeb = ipjk
                  END IF
                END IF
              CASE DEFAULT
                      ijkeeb = ippjkm
              END SELECT
              IF(i == (nx-1)) ijkeeb = ijkeb
!
              SELECT CASE (flag(ipjkpp))
              CASE (slip_wall, noslip_wall, filled_cell_1)
                IF (flag(ipjkp) == slip_wall .OR. flag(ipjkp) == noslip_wall) THEN
                  IF (flag(ijkpp) == slip_wall .OR. flag(ijkpp) == noslip_wall) THEN
                    ijkett = ijkp
                  ELSE
                    ijkett = ijkpp
                  END IF
                ELSE
                  IF (flag(ijkpp) == slip_wall .OR. flag(ijkpp) == noslip_wall) THEN
                    ijkett = ipjkp
                  ELSE
                    ijkett = ijkp
                  END IF
                END IF
              CASE DEFAULT
                      ijkett = ipjkpp
              END SELECT
              IF(k == (nz-1)) ijkett = ijket
!
              SELECT CASE (flag(imjkpp))
              CASE (slip_wall, noslip_wall, filled_cell_1)
                IF (flag(imjkp) == slip_wall .OR. flag(imjkp) == noslip_wall) THEN
                  IF (flag(ijkpp) == slip_wall .OR. flag(ijkpp) == noslip_wall) THEN
                    ijkwtt = ijkp
                  ELSE
                    ijkwtt = ijkpp
                  END IF
                ELSE
                  IF (flag(ijkpp) == slip_wall .OR. flag(ijkpp) == noslip_wall) THEN
                    ijkwtt = imjkp
                  ELSE
                    ijkwtt = ijkp
                  END IF
                END IF
              CASE DEFAULT
                      ijkwtt = imjkpp
              END SELECT
              IF(k == (nz-1)) ijkwtt = ijkwt
!
              SELECT CASE (flag(immjkp))
              CASE (slip_wall, noslip_wall, filled_cell_1)
                IF (flag(imjkp) == slip_wall .OR. flag(imjkp) == noslip_wall) THEN
                  IF (flag(immjk) == slip_wall .OR. flag(immjk) == noslip_wall) THEN
                    ijkwwt = imjk
                  ELSE
                    ijkwwt = immjk
                  END IF
                ELSE
                  IF (flag(immjk) == slip_wall .OR. flag(immjk) == noslip_wall) THEN
                    ijkwwt = imjkp
                  ELSE
                    ijkwwt = imjk
                  END IF
                END IF
              CASE DEFAULT
                      ijkwwt = immjkp
              END SELECT
              IF(i == 2) ijkwwt = ijkwt
!
              SELECT CASE (flag(immjkm))
              CASE (slip_wall, noslip_wall, filled_cell_1)
                IF (flag(imjkm) == slip_wall .OR. flag(imjkm) == noslip_wall) THEN
                  IF (flag(immjk) == slip_wall .OR. flag(immjk) == noslip_wall) THEN
                    ijkwwb = imjk
                  ELSE
                    ijkwwb = immjk
                  END IF
                ELSE
                  IF (flag(immjk) == slip_wall .OR. flag(immjk) == noslip_wall) THEN
                    ijkwwb = imjkm
                  ELSE
                    ijkwwb = imjk
                  END IF
                END IF
              CASE DEFAULT
                      ijkwwb = immjkm
              END SELECT
              IF(i == 2) ijkwwb = ijkwb
!
              SELECT CASE (flag(ipjkmm))
              CASE (slip_wall, noslip_wall, filled_cell_1)
                IF (flag(ipjkm) == slip_wall .OR. flag(ipjkm) == noslip_wall) THEN
                  IF (flag(ijkmm) == slip_wall .OR. flag(ijkmm) == noslip_wall) THEN
                    ijkebb = ijkm
                  ELSE
                    ijkebb = ijkmm
                  END IF
                ELSE
                  IF (flag(ijkmm) == slip_wall .OR. flag(ijkmm) == noslip_wall) THEN
                    ijkebb = ipjkm
                  ELSE
                    ijkebb = ijkm
                  END IF
                END IF
              CASE DEFAULT
                      ijkebb = ipjkmm
              END SELECT
              IF(k == 2) ijkebb = ijkeb
!
              ijkwbb = imjkmm
              SELECT CASE (flag(imjkmm))
              CASE (slip_wall, noslip_wall, filled_cell_1)
                IF (flag(imjkm) == slip_wall .OR. flag(imjkm) == noslip_wall) THEN
                  IF (flag(ijkmm) == slip_wall .OR. flag(ijkmm) == noslip_wall) THEN
                    ijkwbb = ijkm
                  ELSE
                    ijkwbb = ijkmm
                  END IF
                ELSE
                  IF (flag(ijkmm) == slip_wall .OR. flag(ijkmm) == noslip_wall) THEN
                    ijkwbb = imjkm
                  ELSE
                    ijkwbb = ijkm
                  END IF
                END IF
              CASE DEFAULT
                      ijkwbb = imjkmm
              END SELECT
              IF(k == 2) ijkwbb = ijkwb
!
! ... Third neighbours?
!
              SELECT CASE (flag(ijkppp))
              CASE (slip_wall, noslip_wall, filled_cell_1)
                      ijkttt = ijkpp
              CASE DEFAULT
                      ijkttt = ijkppp
              END SELECT
              IF(k == (nz-1)) ijkttt = ijkt
              IF(k == (nz-2)) ijkttt = ijktt
!
              SELECT CASE (flag(ipppjk))
              CASE (slip_wall, noslip_wall, filled_cell_1)
                      ijkeee = ippjk
              CASE DEFAULT
                      ijkeee = ipppjk
              END SELECT
              IF(i == (nx-1)) ijkeee = ijke
              IF(i == (nx-2)) ijkeee = ijkee
!
              ijkeeet = ipppjkp
              SELECT CASE (flag(ipppjkp))
              CASE (slip_wall, noslip_wall, filled_cell_1)
                IF (flag(ippjkp) == slip_wall .OR. flag(ippjkp) == noslip_wall) THEN
                  IF (flag(ipppjk) == slip_wall .OR. flag(ipppjk) == noslip_wall) THEN
                    ijkeeet = ippjk
                  ELSE
                    ijkeeet = ipppjk
                  END IF
                ELSE
                  IF (flag(ipppjk) == slip_wall .OR. flag(ipppjk) == noslip_wall) THEN
                    ijkeeet = ippjkp
                  ELSE
                    ijkeeet = ippjk
                  END IF
                END IF
              CASE DEFAULT
                      ijkeeet = ipppjkp
              END SELECT
              IF(i == (nx-1)) ijkeeet = ijket
              IF(i == (nx-2)) ijkeeet = ijkeet
!
              SELECT CASE (flag(ipjkppp))
              CASE (slip_wall, noslip_wall, filled_cell_1)
                IF (flag(ipjkpp) == slip_wall .OR. flag(ipjkpp) == noslip_wall) THEN
                  IF (flag(ijkppp) == slip_wall .OR. flag(ijkppp) == noslip_wall) THEN
                    ijkettt = ijkpp
                  ELSE
                    ijkettt = ijkppp
                  END IF
                ELSE
                  IF (flag(ijkppp) == slip_wall .OR. flag(ijkppp) == noslip_wall) THEN
                    ijkettt = ipjkpp
                  ELSE
                    ijkettt = ijkpp
                  END IF
                END IF
              CASE DEFAULT
                      ijkettt = ipjkppp
              END SELECT
              IF(k == (nz-1)) ijkettt = ijket
              IF(k == (nz-2)) ijkettt = ijkett
!
              ijkeett = ippjkpp
              SELECT CASE (flag(ippjkpp))
              CASE (slip_wall, noslip_wall, filled_cell_1)
                IF (flag(ipjkpp) == slip_wall .OR. flag(ipjkpp) == noslip_wall) THEN
                  IF (flag(ippjkp) == slip_wall .OR. flag(ippjkp) == noslip_wall) THEN
                    ijkeett = ipjkp
                  ELSE
                    ijkeett = ippjkp
                  END IF
                ELSE
                  IF (flag(ippjkp) == slip_wall .OR. flag(ippjkp) == noslip_wall) THEN
                    ijkeett = ipjkpp
                  ELSE
                    ijkeett = ipjkp
                  END IF
                END IF
              CASE DEFAULT
                      ijkeett = ippjkpp
              END SELECT
              IF (i == nx-1) ijkeett = ijkett
              IF (k == nz-1) ijkeett = ijkeet
!
              ijkwwtt = immjkpp
              SELECT CASE (flag(immjkpp))
              CASE (slip_wall, noslip_wall, filled_cell_1)
                IF (flag(imjkpp) == slip_wall .OR. flag(imjkpp) == noslip_wall) THEN
                  IF (flag(immjkp) == slip_wall .OR. flag(immjkp) == noslip_wall) THEN
                    ijkwwtt = imjkp
                  ELSE
                    ijkwwtt = immjkp
                  END IF
                ELSE
                  IF (flag(immjkp) == slip_wall .OR. flag(immjkp) == noslip_wall) THEN
                    ijkwwtt = imjkpp
                  ELSE
                    ijkwwtt = imjkp
                  END IF
                END IF
              CASE DEFAULT
                      ijkwwtt = immjkpp
              END SELECT
              IF (i == 2)    ijkwwtt = ijkwtt
              IF (k == nz-1) ijkwwtt = ijkwwt
!
              ijkeebb = ippjkmm
              SELECT CASE (flag(ippjkmm))
              CASE (slip_wall, noslip_wall, filled_cell_1)
                IF (flag(ipjkmm) == slip_wall .OR. flag(ipjkmm) == noslip_wall) THEN
                  IF (flag(ippjkm) == slip_wall .OR. flag(ippjkm) == noslip_wall) THEN
                    ijkeebb = ipjkm
                  ELSE
                    ijkeebb = ippjkm
                  END IF
                ELSE
                  IF (flag(ippjkp) == slip_wall .OR. flag(ippjkp) == noslip_wall) THEN
                    ijkeebb = ipjkmm
                  ELSE
                    ijkeebb = ipjkm
                  END IF
                END IF
                      IF (i == nx-1) ijkeebb = ijkeb
                      IF (k == 2)    ijkeebb = ijkeeb
              CASE DEFAULT
                      ijkeebb = ippjkmm
                      IF (i == nx-1) ijkeebb = ijkebb
                      IF (k == 2) ijkeebb = ijkeeb
              END SELECT
!
              myinds( ip1_jp0_kp0_, ijk) = ijke
              myinds( im1_jp0_kp0_, ijk) = ijkw
              myinds( ip0_jp0_kp1_, ijk) = ijkt
              myinds( ip0_jp0_km1_, ijk) = ijkb
              myinds( ip2_jp0_kp0_, ijk) = ijkee
              myinds( im2_jp0_kp0_, ijk) = ijkww
              myinds( ip0_jp0_kp2_, ijk) = ijktt
              myinds( ip0_jp0_km2_, ijk) = ijkbb
              myinds( im1_jp0_km1_, ijk) = ijkwb
              myinds( ip1_jp0_km1_, ijk) = ijkeb
              myinds( im1_jp0_kp1_, ijk) = ijkwt
              myinds( ip1_jp0_kp1_, ijk) = ijket
              myinds( ip2_jp0_kp1_, ijk) = ijkeet
              myinds( ip2_jp0_km1_, ijk) = ijkeeb
              myinds( ip1_jp0_kp2_, ijk) = ijkett
              myinds( im1_jp0_kp2_, ijk) = ijkwtt
              myinds( im2_jp0_kp1_, ijk) = ijkwwt
              myinds( im2_jp0_km1_, ijk) = ijkwwb
              myinds( ip1_jp0_km2_, ijk) = ijkebb
              myinds( im1_jp0_km2_, ijk) = ijkwbb
              myinds( ip3_jp0_kp0_, ijk) = ijkeee
              myinds( ip0_jp0_kp3_, ijk) = ijkttt
              myinds( ip3_jp0_kp1_, ijk) = ijkeeet
              myinds( ip1_jp0_kp3_, ijk) = ijkettt
              myinds( ip2_jp0_kp2_, ijk) = ijkeett
              myinds( ip2_jp0_km2_, ijk) = ijkeebb
              myinds( im2_jp0_kp2_, ijk) = ijkwwtt

            END IF

          ELSE IF( job_type == JOB_TYPE_3D ) THEN

            IF( (i >= 2) .AND. (i <= (nx-1)) .AND.   &
                (j >= 2) .AND. (j <= (ny-1)) .AND.   &
                (k >= 2) .AND. (k <= (nz-1))         ) THEN
!
! ... First neighbours: near the axis or solid boundaries 
! ... impose homogeneous Neumann conditions
!
              ijkm     = myijk( ip0_jp0_km1_, ijk ) 
              imjk     = myijk( im1_jp0_kp0_, ijk ) 
              ipjk     = myijk( ip1_jp0_kp0_, ijk ) 
              ijkp     = myijk( ip0_jp0_kp1_, ijk ) 
              ipjkm    = myijk( ip1_jp0_km1_, ijk )
              ipjkp    = myijk( ip1_jp0_kp1_, ijk )
              imjkm    = myijk( im1_jp0_km1_, ijk )
              imjkp    = myijk( im1_jp0_kp1_, ijk )
              ijkpp    = myijk( ip0_jp0_kp2_, ijk )
              ippjk    = myijk( ip2_jp0_kp0_, ijk )
              immjk    = myijk( im2_jp0_kp0_, ijk )
              ijkmm    = myijk( ip0_jp0_km2_, ijk )
              ijpk     = myijk( ip0_jp1_kp0_, ijk )
              ipjpk    = myijk( ip1_jp1_kp0_, ijk )
              imjpk    = myijk( im1_jp1_kp0_, ijk )
              ijmk     = myijk( ip0_jm1_kp0_, ijk )
              ipjmk    = myijk( ip1_jm1_kp0_, ijk )
              imjmk    = myijk( im1_jm1_kp0_, ijk )
              ijppk    = myijk( ip0_jp2_kp0_, ijk )
              ijmmk    = myijk( ip0_jm2_kp0_, ijk )
              ijpkp    = myijk( ip0_jp1_kp1_, ijk )
              ijmkp    = myijk( ip0_jm1_kp1_, ijk )
              ijpkm    = myijk( ip0_jp1_km1_, ijk )
              ijmkm    = myijk( ip0_jm1_km1_, ijk )
              imjmkp   = myijk( im1_jm1_kp1_, ijk ) 
              ipjmkp   = myijk( ip1_jm1_kp1_, ijk )
              imjpkp   = myijk( im1_jp1_kp1_, ijk )
              ipjpkp   = myijk( ip1_jp1_kp1_, ijk )
              ipjpkm   = myijk( ip1_jp1_km1_, ijk )
              ipjmkm   = myijk( ip1_jm1_km1_, ijk )
              imjpkm   = myijk( im1_jp1_km1_, ijk )
              ippjkp   = myijk( ip2_jp0_kp1_, ijk )
              ippjkm   = myijk( ip2_jp0_km1_, ijk )
              ipjkpp   = myijk( ip1_jp0_kp2_, ijk )
              imjkpp   = myijk( im1_jp0_kp2_, ijk )
              ijmmkp   = myijk( ip0_jm2_kp1_, ijk )
              ippjpk   = myijk( ip2_jp1_kp0_, ijk )
              ippjmk   = myijk( ip2_jm1_kp0_, ijk )
              ipjmmk   = myijk( ip1_jm2_kp0_, ijk )
              ipjmmkp  = myijk( ip1_jm2_kp1_, ijk )
              ippjpkm  = myijk( ip2_jp1_km1_, ijk )
              ippjmkm  = myijk( ip2_jm1_km1_, ijk )
              ippjmkp  = myijk( ip2_jm1_kp1_, ijk )
              imjppkp  = myijk( im1_jp2_kp1_, ijk )
              ijppkp   = myijk( ip0_jp2_kp1_, ijk )
              ipjppkm  = myijk( ip1_jp2_km1_, ijk )
              ipjppk   = myijk( ip1_jp2_kp0_, ijk )
              ipjkmm   = myijk( ip1_jp0_km2_, ijk )
              ipjpkmm  = myijk( ip1_jp1_km2_, ijk )
              imjppkm  = myijk( im1_jp2_km1_, ijk )
              imjppk   = myijk( im1_jp2_kp0_, ijk )
              ijpkmm   = myijk( ip0_jp1_km2_, ijk )
              ijppkm   = myijk( ip0_jp2_km1_, ijk )
              imjmkpp  = myijk( im1_jm1_kp2_, ijk )
              ijmkpp   = myijk( ip0_jm1_kp2_, ijk )
              imjpkpp  = myijk( im1_jp1_kp2_, ijk )
              ijpkpp   = myijk( ip0_jp1_kp2_, ijk )
              ipjmkpp  = myijk( ip1_jm1_kp2_, ijk )
              immjpk   = myijk( im2_jp1_kp0_, ijk )
              immjpkp  = myijk( im2_jp1_kp1_, ijk )
              immjkp   = myijk( im2_jp0_kp1_, ijk )
              imjmkm   = myijk( im1_jm1_km1_, ijk )
              ippjpkp  = myijk( ip2_jp1_kp1_, ijk )
              ipjppkp  = myijk( ip1_jp2_kp1_, ijk )
              ipjpkpp  = myijk( ip1_jp1_kp2_, ijk )
              imjkmm   = myijk( im1_jp0_km2_, ijk )
              ijmmkm   = myijk( ip0_jm2_km1_, ijk )
              immjmk   = myijk( im2_jm1_kp0_, ijk )
              ijmkmm   = myijk( ip0_jm1_km2_, ijk )
              imjmmk   = myijk( im1_jm2_kp0_, ijk )
              immjkm   = myijk( im2_jp0_km1_, ijk )
              imjpkmm  = myijk( im1_jp1_km2_, ijk )
              ipjmmkm  = myijk( ip1_jm2_km1_, ijk )
              immjmkp  = myijk( im2_jm1_kp1_, ijk )
              ipjmkmm  = myijk( ip1_jm1_km2_, ijk )
              imjmkmm  = myijk( im1_jm1_km2_, ijk )
              imjmmkm  = myijk( im1_jm2_km1_, ijk )
              imjmmkp  = myijk( im1_jm2_kp1_, ijk )
              immjpkm  = myijk( im2_jp1_km1_, ijk )
              immjmkm  = myijk( im2_jm1_km1_, ijk )
              ipjpkppp = myijk( ip1_jp1_kp3_, ijk )
              imjpkppp = myijk( im1_jp1_kp3_, ijk )
              ijpkppp  = myijk( ip0_jp1_kp3_, ijk )
              ipjkppp  = myijk( ip1_jp0_kp3_, ijk )
              imjkppp  = myijk( im1_jp0_kp3_, ijk )
              ijkppp   = myijk( ip0_jp0_kp3_, ijk )
              ipjmmkpp = myijk( ip1_jm2_kp2_, ijk )
              ijmmkpp  = myijk( ip0_jm2_kp2_, ijk )
              ipjppkpp = myijk( ip1_jp2_kp2_, ijk )
              ijppkpp  = myijk( ip0_jp2_kp2_, ijk )
              immjmkpp = myijk( im2_jm1_kp2_, ijk )
              immjpkpp = myijk( im2_jp1_kp2_, ijk )
              immjkpp  = myijk( im2_jp0_kp2_, ijk )
              ippjmkpp = myijk( ip2_jm1_kp2_, ijk )
              ippjpkpp = myijk( ip2_jp1_kp2_, ijk )
              ippjkpp  = myijk( ip2_jp0_kp2_, ijk )
              ipjpppkm = myijk( ip1_jp3_km1_, ijk )
              ipjpppkp = myijk( ip1_jp3_kp1_, ijk )
              ijpppkm  = myijk( ip0_jp3_km1_, ijk )
              ijpppkp  = myijk( ip0_jp3_kp1_, ijk )
              ipjpppk  = myijk( ip1_jp3_kp0_, ijk )
              ijpppk   = myijk( ip0_jp3_kp0_, ijk )
              ipjppkmm = myijk( ip1_jp2_km2_, ijk )
              imjppkmm = myijk( im1_jp2_km2_, ijk )
              ijppkmm  = myijk( ip0_jp2_km2_, ijk )
              imjppkpp = myijk( im1_jp2_kp2_, ijk )
              immjppkp = myijk( im2_jp2_kp1_, ijk )
              immjppk  = myijk( im2_jp2_kp0_, ijk )
              ippjppkp = myijk( ip2_jp2_kp1_, ijk )
              ippjppk  = myijk( ip2_jp2_kp0_, ijk )
              ippjpkmm = myijk( ip2_jp1_km2_, ijk )
              ippjkmm  = myijk( ip2_jp0_km2_, ijk )
              ippjmmkm = myijk( ip2_jm2_km1_, ijk )
              ippjmmkp = myijk( ip2_jm2_kp1_, ijk )
              ippjmmk  = myijk( ip2_jm2_kp0_, ijk )
              ippjppkm = myijk( ip2_jp2_km1_, ijk )
              ipppjmkp = myijk( ip3_jm1_kp1_, ijk )
              ipppjpkp = myijk( ip3_jp1_kp1_, ijk )
              ipppjkp  = myijk( ip3_jp0_kp1_, ijk )
              ipppjmk  = myijk( ip3_jm1_kp0_, ijk )
              ipppjpk  = myijk( ip3_jp1_kp0_, ijk )
              ipppjk   = myijk( ip3_jp0_kp0_, ijk )
              ippjmkmm = myijk( ip2_jm1_km2_, ijk )
              immjppkm = myijk( im2_jp2_km1_, ijk )
              imjmmkpp = myijk( im1_jm2_kp2_, ijk )
!  
              SELECT CASE (flag(ipjk))
                CASE (noslip_wall, slip_wall, filled_cell_1)
                        ijke = ijk
                CASE DEFAULT
                        ijke  =  ipjk
              END SELECT
!
              SELECT CASE (flag(imjk))
                CASE (noslip_wall, slip_wall, filled_cell_1)
                        ijkw = ijk
                CASE DEFAULT
                        ijkw  =  imjk
              END SELECT
!
              SELECT CASE (flag(ijpk))
                CASE (noslip_wall, slip_wall, filled_cell_1)
                        ijkn = ijk
                CASE DEFAULT
                        ijkn  =  ijpk
              END SELECT
!
              SELECT CASE (flag(ijmk))
                CASE (noslip_wall, slip_wall, filled_cell_1)
                        ijks = ijk
                CASE DEFAULT
                        ijks  =  ijmk
              END SELECT
!
              SELECT CASE (flag(ijkp))
                CASE (noslip_wall, slip_wall, filled_cell_1)
                        ijkt = ijk
                CASE DEFAULT
                        ijkt  =  ijkp
              END SELECT
!
              SELECT CASE (flag(ijkm))
                CASE (noslip_wall, slip_wall, filled_cell_1)
                        ijkb = ijk
                CASE DEFAULT
                        ijkb  =  ijkm
              END SELECT
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
              SELECT CASE (flag(ijppk))
              CASE (slip_wall, noslip_wall, filled_cell_1)
                      ijknn = ijkn
              CASE DEFAULT
                      ijknn =  ijppk
              END SELECT
              if( (j == (ny-1)) ) ijknn = ijkn

              SELECT CASE (flag(ijmmk))
              CASE (slip_wall, noslip_wall, filled_cell_1)
                      ijkss = ijks
              CASE DEFAULT
                      ijkss =  ijmmk
              END SELECT
              if( (j == 2) ) ijkss = ijks
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
              ijken =  ipjpk
              ijkwn =  imjpk
              ijkes =  ipjmk
              ijkws =  imjmk
              ijket =  ipjkp
              ijkwt =  imjkp
              ijknt =  ijpkp
              ijkst =  ijmkp

              ijkeb =  ipjkm
              SELECT CASE (flag(ipjkm))
                CASE (noslip_wall, slip_wall, filled_cell_1)
                        ijkeb = ipjk
                CASE DEFAULT
                        ijkeb  =  ipjkm
              END SELECT

              ijkwb =  imjkm
              SELECT CASE (flag(imjkm))
                CASE (noslip_wall, slip_wall, filled_cell_1)
                        ijkwb = imjk
                CASE DEFAULT
                        ijkwb  =  imjkm
              END SELECT

              ijknb =  ijpkm
              SELECT CASE (flag(ijpkm))
                CASE (noslip_wall, slip_wall, filled_cell_1)
                        ijknb = ijpk
                CASE DEFAULT
                        ijknb  =  ijpkm
              END SELECT

              ijksb =  ijmkm
              SELECT CASE (flag(ijmkm))
                CASE (noslip_wall, slip_wall, filled_cell_1)
                        ijksb = ijmk
                CASE DEFAULT
                        ijksb  =  ijmkm
              END SELECT

              ijkeet  = ippjkp
              ijkeeb  = ippjkm
              ijkett  = ipjkpp
              ijkwtt  = imjkpp
              ijkwst  = imjmkp
              ijkest  = ipjmkp
              ijkwnt  = imjpkp
              ijkent  = ipjpkp
              ijkenb  = ipjpkm
              ijkesb  = ipjmkm
              ijkwnb  = imjpkm
              ijksst  = ijmmkp
              ijkeen  = ippjpk
              ijkees  = ippjmk
              ijkess  = ipjmmk
              ijkesst = ipjmmkp
              ijkeenb = ippjpkm
              ijkeesb = ippjmkm
              ijkeest = ippjmkp
              ijkwnnt = imjppkp
              ijknnt  = ijppkp
              ijkennb = ipjppkm
              ijkenn  = ipjppk
              ijkebb  = ipjkmm
              ijkenbb = ipjpkmm
              ijkwnnb = imjppkm
              ijkwnn  = imjppk
              ijknbb  = ijpkmm
              ijknnb  = ijppkm
              ijkwstt = imjmkpp
              ijkstt  = ijmkpp
              ijkwntt = imjpkpp
              ijkntt  = ijpkpp
              ijkestt = ipjmkpp
              ijkwwn  = immjpk
              ijkwwnt = immjpkp
              ijkwwt  = immjkp
              ijkwsb  = imjmkm
              ijkeent = ippjpkp
              ijkennt = ipjppkp
              ijkentt = ipjpkpp
!
! wsb
              SELECT CASE (flag(imjmkm))
                CASE (noslip_wall, slip_wall, filled_cell_1)
                        ijkwsb = imjmk
                CASE DEFAULT
                        ijkwsb  =  imjmkm
              END SELECT
!
! esb
              SELECT CASE (flag(ipjmkm))
                CASE (noslip_wall, slip_wall, filled_cell_1)
                        ijkesb = ipjmk
                CASE DEFAULT
                        ijkesb  =  ipjmkm
              END SELECT
!
! wnb
              SELECT CASE (flag(imjpkm))
                CASE (noslip_wall, slip_wall, filled_cell_1)
                        ijkwnb = imjpk
                CASE DEFAULT
                        ijkwnb  =  imjpkm
              END SELECT
!
! enb
              SELECT CASE (flag(ipjpkm))
                CASE (noslip_wall, slip_wall, filled_cell_1)
                        ijkenb = ipjpk
                CASE DEFAULT
                        ijkenb  =  ipjpkm
              END SELECT
!
! eet
              IF (i == nx-1) ijkeet = ijket
!
! een
              IF (i == nx-1) ijkeen = ijken
!
! ees
              IF (i == nx-1) ijkees = ijkes
!
! wwt
              IF (i == 2) ijkwwt = ijkwt
!
! wwn
              IF (i == 2) ijkwwn = ijkwn
!
! sst
              IF (j == 2) ijksst = ijkst
!
! ess
              IF (j == 2) ijkess = ijkes
!
! nnt
              IF (j == ny-1) ijknnt = ijknt
!
! enn
              IF (j == ny-1) ijkenn = ijken
!
! wnn
              IF (j == ny-1) ijkwnn = ijkwn
!
! ett
              IF (k == nz-1) ijkett = ijket
!
! wtt
              IF (k == nz-1) ijkwtt = ijkwt
!
! stt
              IF (k == nz-1) ijkstt = ijkst
!
! ntt
              IF (k == nz-1) ijkntt = ijknt
!
! eeb
              SELECT CASE (flag(ippjkm))
                CASE (noslip_wall, slip_wall, filled_cell_1)
                        ijkeeb = ippjk
                        IF (i == nx-1) ijkeeb = ijke
                CASE DEFAULT
                        ijkeeb  =  ippjkm
                        IF (i == nx-1) ijkeeb = ijkeb
              END SELECT
!
! nnb
              SELECT CASE (flag(ijppkm))
                CASE (noslip_wall, slip_wall, filled_cell_1)
                        ijknnb = ijppk
                        IF (j == ny-1) ijknnb = ijkn
                CASE DEFAULT
                        ijknnb  =  ijppkm
                        IF (j == ny-1) ijknnb = ijknb
              END SELECT
!
! ebb
              SELECT CASE (flag(ipjkmm))
              CASE (slip_wall, noslip_wall, filled_cell_1)
                      ijkebb = ipjkm
              CASE DEFAULT
                      ijkebb = ipjkmm
              END SELECT
              IF(k == 2) ijkebb = ijkeb
!
! nbb
              SELECT CASE (flag(ijpkmm))
              CASE (slip_wall, noslip_wall, filled_cell_1)
                      ijknbb = ijpkm
              CASE DEFAULT
                      ijknbb = ijpkmm
              END SELECT
              IF(k == 2) ijknbb = ijknb
!
! entt
              IF (k == nz-1) ijkentt = ijkent
!
! ennt
              IF (j == ny-1) ijkennt = ijkent
!
! eent
              IF (i == nx-1) ijkeent = ijkent
!
! wwnt
              IF (i == 2) ijkwwnt = ijkwnt
!
! estt
              IF (k == nz-1) ijkestt = ijkest
!
! wntt
              IF (k == nz-1) ijkwntt = ijkwnt
!
! wstt
              IF (k == nz-1) ijkwstt = ijkwst
!
! wnnt
              IF (j == ny-1) ijkwnnt = ijkwnt
!
! eest
              IF (i == nx-1) ijkeest = ijkest
!
! esst
              IF (j == 2) ijkesst = ijkest
!
! eenb
              SELECT CASE (flag(ippjpkm))
                CASE (noslip_wall, slip_wall, filled_cell_1)
                        ijkeenb = ippjpk
                        IF (i == nx-1) ijkeenb = ijken
                CASE DEFAULT
                        ijkeenb  =  ippjpkm
                        IF (i == nx-1) ijkeenb = ijkenb
              END SELECT
!
! eesb
              SELECT CASE (flag(ippjmkm))
                CASE (noslip_wall, slip_wall, filled_cell_1)
                        ijkeesb = ippjmk
                        IF (i == nx-1) ijkeesb = ijkes
                CASE DEFAULT
                        ijkeesb  =  ippjmkm
                        IF (i == nx-1) ijkeesb = ijkesb
              END SELECT
!
! ennb
              SELECT CASE (flag(ipjppkm))
                CASE (noslip_wall, slip_wall, filled_cell_1)
                        ijkennb = ipjppk
                        IF (j == ny-1) ijkennb = ijken
                CASE DEFAULT
                        ijkennb  =  ipjppkm
                        IF (j == ny-1) ijkennb = ijkenb
              END SELECT
!
! wnnb
              SELECT CASE (flag(imjppkm))
                CASE (noslip_wall, slip_wall, filled_cell_1)
                        ijkwnnb = imjppk
                        IF (j == ny-1) ijkwnnb = ijkwn
                CASE DEFAULT
                        ijkwnnb  =  imjppkm
                        IF (j == ny-1) ijkwnnb = ijkwnb
              END SELECT
!
! enbb
              SELECT CASE (flag(ipjpkmm))
                CASE (slip_wall, noslip_wall, filled_cell_1)
                        ijkenbb = ipjpkm
                CASE DEFAULT
                        ijkenbb = ipjpkmm
              END SELECT
              IF(k == 2) ijkenbb = ijkenb
!
! wwb
              ijkwwb = immjkm
              SELECT CASE (flag(immjkm))
                CASE (noslip_wall, slip_wall, filled_cell_1)
                        ijkwwb = immjk
                        IF (i == 2) ijkwwb = ijkw
                CASE DEFAULT
                        ijkwwb  =  immjkm
                        IF (i == 2) ijkwwb = ijkwb
              END SELECT
! wws
              ijkwws = immjmk
              IF (i == 2) ijkwws = ijkws
! wss
              ijkwss = imjmmk
              IF (j == 2) ijkwss = ijkws
! ssb
              ijkssb = ijmmkm
              SELECT CASE (flag(ijmmkm))
                CASE (noslip_wall, slip_wall, filled_cell_1)
                        ijkssb = ijmmk
                        IF (j == 2) ijkssb = ijks
                CASE DEFAULT
                        ijkssb  =  ijmmkm
                        IF (j == 2) ijkssb = ijksb
              END SELECT
! wbb
              ijkwbb = imjkmm
              SELECT CASE (flag(imjkmm))
              CASE (slip_wall, noslip_wall, filled_cell_1)
                      ijkwbb = imjkm
              CASE DEFAULT
                      ijkwbb = imjkmm
              END SELECT
              IF(k == 2) ijkwbb = ijkwb
! sbb
              ijksbb = ijmkmm
              SELECT CASE (flag(ijmkmm))
              CASE (slip_wall, noslip_wall, filled_cell_1)
                      ijksbb = ijmkm
              CASE DEFAULT
                      ijksbb = ijmkmm
              END SELECT
              IF(k == 2) ijksbb = ijksb
              
! ttt
              SELECT CASE (flag(ijkppp))
              CASE (slip_wall, noslip_wall, filled_cell_1)
                      ijkttt = ijkpp
              CASE DEFAULT
                      ijkttt = ijkppp
              END SELECT
              IF(k == (nz-1)) ijkttt = ijkt
              IF(k == (nz-2)) ijkttt = ijktt
! eee
              SELECT CASE (flag(ipppjk))
              CASE (slip_wall, noslip_wall, filled_cell_1)
                      ijkeee = ippjk
              CASE DEFAULT
                      ijkeee = ipppjk
              END SELECT
              IF(i == (nx-1)) ijkeee = ijke
              IF(i == (nx-2)) ijkeee = ijkee
! nnn
              SELECT CASE (flag(ijpppk))
              CASE (slip_wall, noslip_wall, filled_cell_1)
                      ijknnn = ijppk
              CASE DEFAULT
                      ijknnn = ijpppk
              END SELECT
              IF(j == (ny-1)) ijknnn = ijkn
              IF(j == (ny-2)) ijknnn = ijknn

              ijkwwst = immjmkp
              ijkwsst = imjmmkp
              ijkwwnb = immjpkm
              ijkwwsb = immjmkm
              ijkessb = ipjmmkm
              ijkwssb = imjmmkm
              ijkwnbb = imjpkmm
              ijkesbb = ipjmkmm
              ijkwsbb = imjmkmm

! wwst
              IF (i == 2) ijkwwst = ijkwst
!
! wsst
              IF (j == 2) ijkwsst = ijkwst
!
! wwnb
              SELECT CASE (flag(immjpkm))
                CASE (noslip_wall, slip_wall, filled_cell_1)
                        ijkwwnb = immjpk
                        IF (i == 2) ijkwwnb = ijkwn
                CASE DEFAULT
                        ijkwwnb  =  immjpkm
                        IF (i == 2) ijkwwnb = ijkwnb
              END SELECT
!
! wwsb
              SELECT CASE (flag(immjmkm))
                CASE (noslip_wall, slip_wall, filled_cell_1)
                        ijkwwsb = immjmk
                        IF (i == 2) ijkwwsb = ijkws
                CASE DEFAULT
                        ijkwwsb  =  immjmkm
                        IF (i == 2) ijkwwsb = ijkwsb
              END SELECT
!
! essb
              SELECT CASE (flag(ipjmmkm))
                CASE (noslip_wall, slip_wall, filled_cell_1)
                        ijkessb = ipjmmk
                        IF (j == 2) ijkessb = ijkes
                CASE DEFAULT
                        ijkessb  =  ipjmmkm
                        IF (j == 2) ijkessb = ijkesb
              END SELECT
!
! wssb
              SELECT CASE (flag(imjmmkm))
                CASE (noslip_wall, slip_wall, filled_cell_1)
                        ijkwssb = imjmmk
                        IF (j == 2) ijkwssb = ijkws
                CASE DEFAULT
                        ijkwssb  =  imjmmkm
                        IF (j == 2) ijkwssb = ijkwsb
              END SELECT
!
! wnbb
              SELECT CASE (flag(imjpkmm))
                CASE (slip_wall, noslip_wall, filled_cell_1)
                        ijkwnbb = imjpkm
                CASE DEFAULT
                        ijkwnbb = imjpkmm
              END SELECT
              IF(k == 2) ijkwnbb = ijkwnb
!
! esbb
              SELECT CASE (flag(ipjmkmm))
                CASE (slip_wall, noslip_wall, filled_cell_1)
                        ijkesbb = ipjmkm
                CASE DEFAULT
                        ijkesbb = ipjmkmm
              END SELECT
              IF(k == 2) ijkesbb = ijkesb
!
! wsbb
              SELECT CASE (flag(imjmkmm))
              CASE (slip_wall, noslip_wall, filled_cell_1)
                      ijkwsbb = imjmkm
              CASE DEFAULT
                      ijkwsbb = imjmkmm
              END SELECT
              IF(k == 2) ijkwsbb = ijkwsb
!
! sstt
              ijksstt = ijmmkpp
              IF (j == 2) ijksstt = ijkstt
              IF (k == nz-1) ijksstt = ijksst
!
! nntt
              ijknntt = ijppkpp
              IF (j == ny-1) ijknntt = ijkntt
              IF (k == nz-1) ijknntt = ijknnt
!
! wwtt
              ijkwwtt = immjkpp
              IF (i == 2)    ijkwwtt = ijkwtt
              IF (k == nz-1) ijkwwtt = ijkwwt
!
! eett
              ijkeett = ippjkpp
              IF (i == nx-1) ijkeett = ijkett
              IF (k == nz-1) ijkeett = ijkeet
!
! eenn
              ijkeenn = ippjppk
              IF (i == nx-1) ijkeenn = ijkenn
              IF (j == ny-1) ijkeenn = ijkeen
!
! wwnn
              ijkwwnn = immjppk
              IF (i == 2)    ijkwwnn = ijkwnn
              IF (j == ny-1) ijkwwnn = ijkwwn
!
! eess
              ijkeess = ippjmmk
              IF (i == nx-1) ijkeess = ijkess
              IF (j == 2)    ijkeess = ijkees
!
! nnbb
              ijknnbb = ijppkmm
              SELECT CASE (flag(ijppkmm))
              CASE (slip_wall, noslip_wall, filled_cell_1)
                      ijknnbb = ijppkm
                      IF (j == ny-1) ijknnbb = ijknb
                      IF (k == 2) ijknnbb = ijknnb
              CASE DEFAULT
                      ijknnbb = ijppkmm
                      IF (j == ny-1) ijknnbb = ijknbb
                      IF (k == 2) ijknnbb = ijknnb
              END SELECT
!
! eebb
              ijkeebb = ippjkmm
              SELECT CASE (flag(ippjkmm))
              CASE (slip_wall, noslip_wall, filled_cell_1)
                      ijkeebb = ippjkm
                      IF (i == nx-1) ijkeebb = ijkeb
                      IF (k == 2)    ijkeebb = ijkeeb
              CASE DEFAULT
                      ijkeebb = ippjkmm
                      IF (i == nx-1) ijkeebb = ijkebb
                      IF (k == 2) ijkeebb = ijkeeb
              END SELECT
!
! eeen
              ijkeeen = ipppjpk
              IF(i == (nx-1)) ijkeeen = ijken
              IF(i == (nx-2)) ijkeeen = ijkeen
! eees
              ijkeees = ipppjmk
              IF(i == (nx-1)) ijkeees = ijkes
              IF(i == (nx-2)) ijkeees = ijkees
! eeet
              ijkeeet = ipppjkp
              IF(i == (nx-1)) ijkeeet = ijket
              IF(i == (nx-2)) ijkeeet = ijkeet
! nnnt
              ijknnnt = ijpppkp
              IF(j == (ny-1)) ijknnnt = ijknt
              IF(j == (ny-2)) ijknnnt = ijknnt
! ennn
              ijkennn = ipjpppk
              IF(j == (ny-1)) ijkennn = ijken
              IF(j == (ny-2)) ijkennn = ijkenn
! ettt
              ijkettt = ipjkppp
              IF(k == (nz-1)) ijkettt = ijket
              IF(k == (nz-2)) ijkettt = ijkett
! wttt
              ijkwttt = imjkppp
              IF(k == (nz-1)) ijkwttt = ijkwt
              IF(k == (nz-2)) ijkwttt = ijkwtt
! nttt
              ijknttt = ijpkppp
              IF(k == (nz-1)) ijknttt = ijknt
              IF(k == (nz-2)) ijknttt = ijkntt
! nnnb
              ijknnnb = ijpppkm
              SELECT CASE (flag(ijpppkm))
                CASE (noslip_wall, slip_wall, filled_cell_1)
                        ijknnnb = ijpppk
                        IF(j == (ny-1)) ijknnnb = ijkn
                        IF(j == (ny-2)) ijknnnb = ijknn
                CASE DEFAULT
                        ijknnnb  =  ijpppkm
                        IF(j == (ny-1)) ijknnnb = ijknb
                        IF(j == (ny-2)) ijknnnb = ijknnb
              END SELECT
!
! eeent
              ijkeeent = ipppjpkp
              IF(i == (nx-1)) ijkeeent = ijkent
              IF(i == (nx-2)) ijkeeent = ijkeent
! eeest
              ijkeeest = ipppjmkp
              IF(i == (nx-1)) ijkeeest = ijkest
              IF(i == (nx-2)) ijkeeest = ijkeest
! ennnb
              ijkennnb = ipjpppkm
              SELECT CASE (flag(ipjpppkm))
                CASE (noslip_wall, slip_wall, filled_cell_1)
                        ijkennnb = ipjpppk
                        IF(j == (ny-1)) ijkennnb = ijken
                        IF(j == (ny-2)) ijkennnb = ijkenn
                CASE DEFAULT
                        ijkennnb  =  ipjpppkm
                        IF(j == (ny-1)) ijkennnb = ijkenb
                        IF(j == (ny-2)) ijkennnb = ijkennb
              END SELECT
! ennnt
              ijkennnt = ipjpppkp
              IF(j == (ny-1)) ijkennnt = ijkent
              IF(j == (ny-2)) ijkennnt = ijkennt
! enttt
              ijkenttt = ipjpkppp
              IF(k == (nz-1)) ijkenttt = ijkent
              IF(k == (nz-2)) ijkenttt = ijkentt
! wnttt
              ijkwnttt = imjpkppp
              IF(k == (nz-1)) ijkwnttt = ijkwnt
              IF(k == (nz-2)) ijkwnttt = ijkwntt
!
! eennt
              ijkeennt = ippjppkp
              IF (i == nx-1) ijkeennt = ijkennt
              IF (j == ny-1) ijkeennt = ijkeent
! eesst
              ijkeesst = ippjmmkp
              IF (i == nx-1) ijkeesst = ijkesst
              IF (j == 2)    ijkeesst = ijkeest
! wwnnt
              ijkwwnnt = immjppkp
              IF (i == 2)    ijkwwnnt = ijkwnnt
              IF (j == ny-1) ijkwwnnt = ijkwwnt
! eestt
              ijkeestt = ippjmkpp
              IF (i == nx-1) ijkeestt = ijkestt
              IF (k == nz-1) ijkeestt = ijkeest
! eentt
              ijkeentt = ippjpkpp
              IF (i == nx-1) ijkeentt = ijkentt
              IF (k == nz-1) ijkeentt = ijkeent
! wwntt
              ijkwwntt = immjpkpp
              IF (i == 2)    ijkwwntt = ijkwntt
              IF (k == nz-1) ijkwwntt = ijkwwnt
! wwstt
              ijkwwstt = immjmkpp
              IF (i == 2)    ijkwwstt = ijkwstt
              IF (k == nz-1) ijkwwstt = ijkwwst
! enntt
              ijkenntt = ipjppkpp
              IF (j == ny-1) ijkenntt = ijkentt
              IF (k == nz-1) ijkenntt = ijkennt
! wnntt
              ijkwnntt = imjppkpp
              IF (j == ny-1) ijkwnntt = ijkwntt
              IF (k == nz-1) ijkwnntt = ijkwnnt
! esstt
              ijkesstt = ipjmmkpp
              IF (j == 2)    ijkesstt = ijkestt
              IF (k == nz-1) ijkesstt = ijkesst
! wsstt
              ijkwsstt = imjmmkpp
              IF (j == 2)    ijkwsstt = ijkwstt
              IF (k == nz-1) ijkwsstt = ijkwsst
! ennbb
              ijkennbb = ipjppkmm
              SELECT CASE (flag(ipjppkmm))
              CASE (slip_wall, noslip_wall, filled_cell_1)
                      ijkennbb = ipjppkm
                      IF (j == ny-1) ijkennbb = ijkenb
                      IF (k == 2)    ijkennbb = ijkennb
              CASE DEFAULT
                      ijkennbb = ipjppkmm
                      IF (j == ny-1) ijkennbb = ijkenbb
                      IF (k == 2)    ijkennbb = ijkennb
              END SELECT
! wnnbb
              ijkwnnbb = imjppkmm
              SELECT CASE (flag(imjppkmm))
              CASE (slip_wall, noslip_wall, filled_cell_1)
                      ijkwnnbb = imjppkm
                      IF (j == ny-1) ijkwnnbb = ijkwnb
                      IF (k == 2)    ijkwnnbb = ijkwnnb
              CASE DEFAULT
                      ijkwnnbb = imjppkmm
                      IF (j == ny-1) ijkwnnbb = ijkwnbb
                      IF (k == 2)    ijkwnnbb = ijkwnnb
              END SELECT
! eenbb
              ijkeenbb = ippjpkmm
              SELECT CASE (flag(ippjpkmm))
              CASE (slip_wall, noslip_wall, filled_cell_1)
                      ijkeenbb = ippjpkm
                      IF (i == nx-1) ijkeenbb = ijkenb
                      IF (k == 2)    ijkeenbb = ijkeenb
              CASE DEFAULT
                      ijkeenbb = ippjpkmm
                      IF (i == nx-1) ijkeenbb = ijkenbb
                      IF (k == 2)    ijkeenbb = ijkeenb
              END SELECT
! eesbb
              ijkeesbb = ippjmkmm
              SELECT CASE (flag(ippjmkmm))
              CASE (slip_wall, noslip_wall, filled_cell_1)
                      ijkeesbb = ippjmkm
                      IF (i == nx-1) ijkeesbb = ijkesb
                      IF (k == 2)    ijkeesbb = ijkeesb
              CASE DEFAULT
                      ijkeesbb = ippjmkmm
                      IF (i == nx-1) ijkeesbb = ijkesbb
                      IF (k == 2)    ijkeesbb = ijkeesb
              END SELECT
! eennb
              SELECT CASE (flag(ippjppkm))
                CASE (noslip_wall, slip_wall, filled_cell_1)
                        ijkeennb = ippjppk
                        IF (j == ny-1) ijkeennb = ijkeen
                        IF (i == nx-1) ijkeennb = ijkenn
                CASE DEFAULT
                        ijkeennb  =  ippjppkm
                        IF (j == ny-1) ijkeennb = ijkeenb
                        IF (i == nx-1) ijkeennb = ijkennb
              END SELECT
! wwnnb
              SELECT CASE (flag(immjppkm))
                CASE (noslip_wall, slip_wall, filled_cell_1)
                        ijkwwnnb = immjppk
                        IF (j == ny-1) ijkwwnnb = ijkwwn
                        IF (i == 2)    ijkwwnnb = ijkwnn
                CASE DEFAULT
                        ijkwwnnb  =  immjppkm
                        IF (j == ny-1) ijkwwnnb = ijkwwnb
                        IF (i == 2)    ijkwwnnb = ijkwnnb
              END SELECT
! eessb
              SELECT CASE (flag(ippjmmkm))
                CASE (noslip_wall, slip_wall, filled_cell_1)
                        ijkeessb = ippjmmk
                        IF (j == 2)    ijkeessb = ijkees
                        IF (i == nx-1) ijkeessb = ijkess
                CASE DEFAULT
                        ijkeessb  =  ippjmmkm
                        IF (j == 2)    ijkeessb = ijkeesb
                        IF (i == nx-1) ijkeessb = ijkessb
              END SELECT
!
              myinds( ip0_jp0_km1_ , ijk ) = ijkb
              myinds( im1_jp0_kp0_ , ijk ) = ijkw
              myinds( ip1_jp0_kp0_ , ijk ) = ijke
              myinds( ip0_jp0_kp1_ , ijk ) = ijkt
              myinds( ip1_jp0_km1_ , ijk ) = ijkeb
              myinds( ip1_jp0_kp1_ , ijk ) = ijket
              myinds( im1_jp0_km1_ , ijk ) = ijkwb
              myinds( im1_jp0_kp1_ , ijk ) = ijkwt
              myinds( ip0_jp0_kp2_ , ijk ) = ijktt
              myinds( ip2_jp0_kp0_ , ijk ) = ijkee
              myinds( im2_jp0_kp0_ , ijk ) = ijkww
              myinds( ip0_jp0_km2_ , ijk ) = ijkbb
              myinds( ip0_jp1_kp0_ , ijk ) = ijkn
              myinds( ip1_jp1_kp0_ , ijk ) = ijken
              myinds( im1_jp1_kp0_ , ijk ) = ijkwn
              myinds( ip0_jm1_kp0_ , ijk ) = ijks
              myinds( ip1_jm1_kp0_ , ijk ) = ijkes
              myinds( im1_jm1_kp0_ , ijk ) = ijkws
              myinds( ip0_jp2_kp0_ , ijk ) = ijknn
              myinds( ip0_jm2_kp0_ , ijk ) = ijkss
              myinds( ip0_jp1_kp1_ , ijk ) = ijknt
              myinds( ip0_jm1_kp1_ , ijk ) = ijkst
              myinds( ip0_jp1_km1_ , ijk ) = ijknb
              myinds( ip0_jm1_km1_ , ijk ) = ijksb
              myinds( im1_jm1_kp1_ , ijk ) = ijkwst 
              myinds( ip1_jm1_kp1_ , ijk ) = ijkest
              myinds( im1_jp1_kp1_ , ijk ) = ijkwnt
              myinds( ip1_jp1_kp1_ , ijk ) = ijkent
              myinds( ip1_jp1_km1_ , ijk ) = ijkenb
              myinds( ip1_jm1_km1_ , ijk ) = ijkesb
              myinds( im1_jp1_km1_ , ijk ) = ijkwnb
              myinds( ip2_jp0_kp1_ , ijk ) = ijkeet
              myinds( ip2_jp0_km1_ , ijk ) = ijkeeb
              myinds( ip1_jp0_kp2_ , ijk ) = ijkett
              myinds( im1_jp0_kp2_ , ijk ) = ijkwtt
              myinds( ip0_jm2_kp1_ , ijk ) = ijksst
              myinds( ip2_jp1_kp0_ , ijk ) = ijkeen
              myinds( ip2_jm1_kp0_ , ijk ) = ijkees
              myinds( ip1_jm2_kp0_ , ijk ) = ijkess
              myinds( ip1_jm2_kp1_ , ijk ) = ijkesst
              myinds( ip2_jp1_km1_ , ijk ) = ijkeenb
              myinds( ip2_jm1_km1_ , ijk ) = ijkeesb
              myinds( ip2_jm1_kp1_ , ijk ) = ijkeest
              myinds( im1_jp2_kp1_ , ijk ) = ijkwnnt
              myinds( ip0_jp2_kp1_ , ijk ) = ijknnt
              myinds( ip1_jp2_km1_ , ijk ) = ijkennb
              myinds( ip1_jp2_kp0_ , ijk ) = ijkenn
              myinds( ip1_jp0_km2_ , ijk ) = ijkebb
              myinds( ip1_jp1_km2_ , ijk ) = ijkenbb
              myinds( im1_jp2_km1_ , ijk ) = ijkwnnb
              myinds( im1_jp2_kp0_ , ijk ) = ijkwnn
              myinds( ip0_jp1_km2_ , ijk ) = ijknbb
              myinds( ip0_jp2_km1_ , ijk ) = ijknnb
              myinds( im1_jm1_kp2_ , ijk ) = ijkwstt
              myinds( ip0_jm1_kp2_ , ijk ) = ijkstt
              myinds( im1_jp1_kp2_ , ijk ) = ijkwntt
              myinds( ip0_jp1_kp2_ , ijk ) = ijkntt
              myinds( ip1_jm1_kp2_ , ijk ) = ijkestt
              myinds( im2_jp1_kp0_ , ijk ) = ijkwwn
              myinds( im2_jp1_kp1_ , ijk ) = ijkwwnt
              myinds( im2_jp0_kp1_ , ijk ) = ijkwwt
              myinds( im1_jm1_km1_ , ijk ) = ijkwsb
              myinds( ip2_jp1_kp1_ , ijk ) = ijkeent
              myinds( ip1_jp2_kp1_ , ijk ) = ijkennt
              myinds( ip1_jp1_kp2_ , ijk ) = ijkentt
              myinds( im1_jp0_km2_ , ijk ) = ijkwbb
              myinds( ip0_jm2_km1_ , ijk ) = ijkssb
              myinds( im2_jm1_kp0_ , ijk ) = ijkwws
              myinds( ip0_jm1_km2_ , ijk ) = ijksbb
              myinds( im1_jm2_kp0_ , ijk ) = ijkwss
              myinds( im2_jp0_km1_ , ijk ) = ijkwwb
              myinds( ip3_jp0_kp0_ , ijk ) = ijkeee
              myinds( ip0_jp3_kp0_ , ijk ) = ijknnn
              myinds( ip0_jp0_kp3_ , ijk ) = ijkttt
              myinds( im2_jm1_kp1_ , ijk ) = ijkwwst
              myinds( im1_jm2_kp1_ , ijk ) = ijkwsst
              myinds( im2_jp1_km1_ , ijk ) = ijkwwnb
              myinds( im2_jm1_km1_ , ijk ) = ijkwwsb
              myinds( ip1_jm2_km1_ , ijk ) = ijkessb
              myinds( im1_jm2_km1_ , ijk ) = ijkwssb
              myinds( im1_jp1_km2_ , ijk ) = ijkwnbb
              myinds( ip1_jm1_km2_ , ijk ) = ijkesbb
              myinds( im1_jm1_km2_ , ijk ) = ijkwsbb
              myinds( ip0_jm2_kp2_ , ijk ) = ijksstt
              myinds( ip0_jp2_kp2_ , ijk ) = ijknntt
              myinds( im2_jp0_kp2_ , ijk ) = ijkwwtt
              myinds( ip2_jp0_kp2_ , ijk ) = ijkeett
              myinds( ip2_jp2_kp0_ , ijk ) = ijkeenn
              myinds( im2_jp2_kp0_ , ijk ) = ijkwwnn
              myinds( ip2_jm2_kp0_ , ijk ) = ijkeess
              myinds( ip0_jp2_km2_ , ijk ) = ijknnbb
              myinds( ip2_jp0_km2_ , ijk ) = ijkeebb
              myinds( ip3_jp1_kp0_ , ijk ) = ijkeeen
              myinds( ip3_jm1_kp0_ , ijk ) = ijkeees
              myinds( ip3_jp0_kp1_ , ijk ) = ijkeeet
              myinds( ip0_jp3_km1_ , ijk ) = ijknnnb
              myinds( ip0_jp3_kp1_ , ijk ) = ijknnnt
              myinds( ip1_jp3_kp0_ , ijk ) = ijkennn
              myinds( ip1_jp0_kp3_ , ijk ) = ijkettt
              myinds( im1_jp0_kp3_ , ijk ) = ijkwttt
              myinds( ip0_jp1_kp3_ , ijk ) = ijknttt
              myinds( ip3_jp1_kp1_ , ijk ) = ijkeeent
              myinds( ip3_jm1_kp1_ , ijk ) = ijkeeest
              myinds( ip1_jp3_km1_ , ijk ) = ijkennnb
              myinds( ip1_jp3_kp1_ , ijk ) = ijkennnt
              myinds( ip1_jp1_kp3_ , ijk ) = ijkenttt
              myinds( im1_jp1_kp3_ , ijk ) = ijkwnttt
              myinds( ip2_jp2_kp1_ , ijk ) = ijkeennt
              myinds( ip2_jp2_km1_ , ijk ) = ijkeennb
              myinds( ip2_jm2_kp1_ , ijk ) = ijkeesst
              myinds( ip2_jm2_km1_ , ijk ) = ijkeessb
              myinds( im2_jp2_kp1_ , ijk ) = ijkwwnnt
              myinds( ip2_jp1_kp2_ , ijk ) = ijkeentt
              myinds( ip2_jm1_kp2_ , ijk ) = ijkeestt
              myinds( ip2_jp1_km2_ , ijk ) = ijkeenbb
              myinds( im2_jp1_kp2_ , ijk ) = ijkwwntt
              myinds( im2_jm1_kp2_ , ijk ) = ijkwwstt
              myinds( ip1_jp2_kp2_ , ijk ) = ijkenntt
              myinds( im1_jp2_kp2_ , ijk ) = ijkwnntt
              myinds( ip1_jm2_kp2_ , ijk ) = ijkesstt
              myinds( ip1_jp2_km2_ , ijk ) = ijkennbb
              myinds( im1_jp2_km2_ , ijk ) = ijkwnnbb
              myinds( ip2_jm1_km2_ , ijk ) = ijkeesbb
              myinds( im2_jp2_km1_ , ijk ) = ijkwwnnb
              myinds( im1_jm2_kp2_ , ijk ) = ijkwsstt

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
              vf(ijk) = vf(ijk) / 4
            ELSE IF( job_type == JOB_TYPE_3D) THEN
              vf(ijk) = vf(ijk) / 6
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
      SUBROUTINE face_fractions
      USE grid, ONLY: flag
      USE immersed_boundaries, ONLY: numx, numy, numz, bd
      USE immersed_boundaries, ONLY: b_e_, b_w_, b_n_, b_s_, b_t_, b_b_, ivf_, vf
      INTEGER :: ijk, num
!
      DO ijk=1, ncint
      IF( BTEST(flag(ijk),0) ) THEN
  
          IF (numx(ijk)/=0 .AND. numz(ijk) /=0) THEN
            b_e_(ijk) = 0
            b_w_(ijk) = 0
            b_t_(ijk) = 0
            b_b_(ijk) = 0
            RETURN
          END IF
    
          b_e_(ijk) = 0
          num = 1
          IF( IAND(bd(ijk),num) /= 0 )  b_e_(ijk) = 1
    
          b_w_(ijk) = 0
          num = 2
          IF( IAND(bd(ijk),num) /= 0 )  b_w_(ijk) = 1
    
          b_t_(ijk) = 0
          num = 4
          IF( IAND(bd(ijk),num) /= 0 )  b_t_(ijk) = 1
    
          b_b_(ijk) = 0
          num = 8
          IF( IAND(bd(ijk),num) /= 0 )  b_b_(ijk) = 1
    
          b_n_(ijk) = 0
          num = 16
          IF( IAND(bd(ijk),num) /= 0 ) b_n_(ijk) = 1
    
          b_s_(ijk) = 0
          num = 32
          IF( IAND(bd(ijk),num) /= 0 ) b_s_(ijk) = 1
    
          IF (vf(ijk) > 0.D0) THEN
            ivf_(ijk) = 1.D0 / vf(ijk)
          ELSE
            ivf_(ijk) = 0.D0
          END IF
!
        END IF
      END DO
!
      CALL data_exchange(b_e_)
      CALL data_exchange(b_n_)
      CALL data_exchange(b_t_)
      CALL data_exchange(b_w_)
      CALL data_exchange(b_s_)
      CALL data_exchange(b_b_)
!
      RETURN
      END SUBROUTINE face_fractions
!----------------------------------------------------------------------
      SUBROUTINE face_init
      USE immersed_boundaries, ONLY: b_e_, b_w_, b_n_, b_s_, b_t_, b_b_, ivf_
      IMPLICIT NONE
!
      ALLOCATE(b_e_(ncdom))
      ALLOCATE(b_w_(ncdom))
      ALLOCATE(b_n_(ncdom))
      ALLOCATE(b_s_(ncdom))
      ALLOCATE(b_t_(ncdom))
      ALLOCATE(b_b_(ncdom))
      ALLOCATE(ivf_(ncdom))
!
      b_e_(:) = 1
      b_w_(:) = 1
      b_n_(:) = 1
      b_s_(:) = 1
      b_t_(:) = 1
      b_b_(:) = 1
      ivf_(:) = 1.D0
!
      RETURN
      END SUBROUTINE face_init
!----------------------------------------------------------------------
   END MODULE domain_mapping
!----------------------------------------------------------------------
