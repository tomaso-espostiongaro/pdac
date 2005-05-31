!----------------------------------------------------------------------
      MODULE domain_decomposition
!
! ... This module contains all the procedures used for the domain
! ... decomposition and interprocessor data exchange
!
!----------------------------------------------------------------------
        USE grid, ONLY: fl, flag
        USE grid, ONLY: slip_wall, noslip_wall, filled_cell
        USE indijk_module
        USE immersed_boundaries, ONLY: immb
        USE control_flags, ONLY: lpr
        USE io_files, ONLY: errorunit, testunit, logunit
!
        IMPLICIT NONE
        SAVE
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
        INTEGER, PARAMETER :: LAYER_MAP   = 1
        INTEGER, PARAMETER :: BLOCK2D_MAP = 2
        INTEGER, PARAMETER :: BLOCK3D_MAP = 3
!
        TYPE cells_map_type
          INTEGER :: type              ! Identify the map type (1 layer, 2 columns, 3 blocks)
          INTEGER :: lay(2)            ! lay(1) -> imesh start, lay(2) -> imesh end (layer)
          INTEGER :: corner1(2)        ! corner1(1) -> x coord., corner1(2) -> z coord. 
          INTEGER :: corner2(2)        ! corner2(1) -> x coord., corner2(2) -> z coord.
          INTEGER :: blkbsw(3)         ! blkbsw(.) -> BSW x, y, z coordinates 
          INTEGER :: blktne(3)         ! blktne(.) -> TNE x, y, z coordinates
        END TYPE
!
        TYPE (rcv_map_type), ALLOCATABLE :: rcv_map(:)
! ...     This array store, for each processor, the indexes of local
! ...     cell where to copy the data sent by other procesor

        TYPE (snd_map_type), ALLOCATABLE :: snd_map(:)
! ...     This array store, for each processor, the indexes of local
! ...     cell that the present procesor has to send to the other

        TYPE (cells_map_type), ALLOCATABLE :: proc_map(:)
! ...     This array is used to store data about which cell belong
! ...     to which processor, and is used to find owner, local and
! ...     global index of a cell

! ...   The following arrays contain information on how the cells
!       are distributed among processors.
!
        INTEGER, ALLOCATABLE :: nctot(:)  ! unbalanced number of cell
        INTEGER, ALLOCATABLE :: ncfl1(:)  ! number of cell with fl=1
        INTEGER, ALLOCATABLE :: ncell(:)  ! balanced number of cell
        INTEGER, ALLOCATABLE :: ncdif(:)  ! ncfl1 - ncell
        PRIVATE :: nctot, ncfl1, ncell, ncdif
!
        !  nctot unbalanced number of cell,
        !        this is the number of cells on each process, as if they were 
        !        equally distributed among processors, regardless of their 
        !        computational weight.
        !  ncfl1 number of cells on each process with the flag equal 1
        !  ncell bulanced number of cell, this is the number of cells on
        !        each process that it is computed to obtain the same 
        !        computational load on each process.
!
        INTEGER :: mesh_partition

!       ncint (ncint) number of local cells, 
!       ncext         number of external (ghost) cells, 
!       ncdom         sum of the previous two
!
        INTEGER :: ncint, ncext, ncdom  
!
        INTEGER, ALLOCATABLE :: myijk(:,:)
        INTEGER, ALLOCATABLE :: myinds(:,:)
!
        INTEGER :: countfl

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

!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE partition( nproc, mpime, root )
!
! ... a routine that subdivides a regular grid between (nproc) 
! ... processors, balancing the computational load.
! ... Gives out the processor maps.
! ... (2D/3D-Compliant)
!
      USE control_flags, ONLY: job_type
      USE dimensions
      IMPLICIT NONE
      SAVE

      INTEGER, INTENT(IN) :: nproc, mpime, root
!
      INTEGER :: i, j, k, ijk, n, m, q, rr
      INTEGER :: nfl, ipe, dim
      INTEGER :: layer
      INTEGER :: localdim
      LOGICAL :: ionode
!
! ... Subroutine Body
!
      ionode = ( mpime == root )

      IF ( lpr > 0 .AND. ionode ) THEN
        WRITE( 6, * ) 
        WRITE( 6, * ) 'Starting Grid Decomposition ...'
      END IF

      IF (ALLOCATED(nctot)) DEALLOCATE(nctot)
      IF (ALLOCATED(ncell)) DEALLOCATE(ncell)
      IF (ALLOCATED(ncfl1)) DEALLOCATE(ncfl1)
      IF (ALLOCATED(ncdif)) DEALLOCATE(ncdif)

      ALLOCATE(nctot(0:nproc-1))
      ALLOCATE(ncfl1(0:nproc-1))
      ALLOCATE(ncell(0:nproc-1))
      ALLOCATE(ncdif(0:nproc-1))

      nctot = 0
      ncfl1 = 0
      ncell = 0
      ncdif = 0
!
! ... here count the flags of fluid cells (including immersed boundaries and
! ... dome cells), array fl is replicated on all procs
!
      countfl = COUNT( BTEST(fl(:),0) )
!
      IF ( lpr > 0 .AND. ionode ) THEN
        WRITE( 6, fmt = '(a13)' ) ' # countfl :'
        WRITE( 6, fmt = '(10I8)') countfl
      END IF
!
! ... domain decomposition (build maps)
!
      IF( nproc <= 2 ) mesh_partition = LAYER_MAP
!
      IF( ALLOCATED( proc_map ) ) DEALLOCATE( proc_map )
      ALLOCATE( proc_map(0:nproc-1) )

      proc_map(0:nproc-1)%type = mesh_partition

      IF ( proc_map(0)%type == LAYER_MAP ) THEN
        CALL layers( ncfl1, nctot, nproc, mpime, root )
      ELSE IF ( proc_map(0)%type == BLOCK2D_MAP .AND. job_type == '2D') THEN
        CALL blocks( ncfl1, nctot, nproc, mpime, root )
      ELSE IF ( proc_map(0)%type == BLOCK2D_MAP .AND. job_type == '3D') THEN
        CALL columns( ncfl1, nctot, nproc, mpime, root )
      ELSE IF ( proc_map(0)%type == BLOCK3D_MAP .AND. job_type == '3D') THEN
        CALL blocks3d( ncfl1, nctot, nproc, mpime, root )
      ELSE
        CALL error(' partition ',' partition type not yet implemented ',proc_map(0)%type)
      END IF
!
! ... ncell is the optimal balanced number of cells for each processor
! ... which always differs from the one calculated by the subroutine that actually
! ... distribute the data, this is because there are constraints in the
! ... data distributions itself
! 
      DO ipe = 0, nproc - 1
        ncell(ipe) = localdim(countfl, nproc, ipe)
        ncdif(ipe) = ncfl1(ipe) - ncell(ipe)
      END DO
!
      IF ( lpr > 0 .AND. ionode ) THEN
        WRITE(logunit,*) 
        WRITE(logunit,*) '---  Mesh decomposition summary  ---'
        WRITE(logunit,*) 
        DO ipe = 0, nproc - 1
          WRITE(logunit,*) ' # nctot( ',ipe, ' ) = ', nctot(ipe)
          WRITE(logunit,*) ' # ncfl1( ',ipe, ' ) = ', ncfl1(ipe)
          WRITE(logunit,*) ' # ncdif( ',ipe, ' ) = ', ncfl1(ipe),' -', ncell(ipe),' =', ncdif(ipe)
          IF( proc_map(0)%type == LAYER_MAP ) THEN
            WRITE(logunit,*) '   proc(', ipe,')_map (first-last):', &
              proc_map(ipe)%lay(1), proc_map(ipe)%lay(2)
          ELSE IF(  proc_map(0)%type == BLOCK2D_MAP ) THEN
            WRITE(logunit,*) '   proc(', ipe,')_map (left-bot.):', &
              proc_map(ipe)%corner1(1), proc_map(ipe)%corner1(2)
            WRITE(logunit,*) '   proc(', ipe,')_map (right-top):', &
              proc_map(ipe)%corner2(1), proc_map(ipe)%corner2(2)
          ELSE IF(  proc_map(0)%type == BLOCK3D_MAP ) THEN
            WRITE(logunit,*) '   proc(', ipe,')_map (BSW):', &
              proc_map(ipe)%blkbsw(1), proc_map(ipe)%blkbsw(2), proc_map(ipe)%blkbsw(3)
            WRITE(logunit,*) '   proc(', ipe,')_map (TNE):', &
              proc_map(ipe)%blktne(1), proc_map(ipe)%blktne(2), proc_map(ipe)%blktne(3)
          END IF 
          WRITE(logunit,*) '-----------------------------'
        END DO
      END IF

      IF( ionode ) &
        WRITE( 6, * ) 'End of Mesh decomposition'
!
      RETURN
      END SUBROUTINE partition
!
! ------------------------------------------------------------------------
!
      SUBROUTINE layers( ncfl1, nctot, nproc, mpime, root )
! 
! ... partition N.1 (layers)
! ... decomposes the domain into N horizontal layers.
! ... The layer map is represented by the global index 
! ... of the first and the last cell belonging to the
! ... processor. The cells within the layer are ordered
! ... following the main-sweep.
! ... (2D/3D-Compliant)
!
      USE dimensions
      USE control_flags, ONLY: job_type
!
      IMPLICIT NONE
      SAVE
!
      INTEGER :: ncfl1(0:)
      INTEGER :: nctot(0:)
      INTEGER, INTENT(IN) :: nproc, mpime, root
!
      INTEGER :: i, j, k, ijk
      INTEGER :: icnt, icnt_ipe, ipe
      INTEGER :: localdim
      INTEGER, ALLOCATABLE :: lay_map(:,:)
      LOGICAL :: ionode
!
      ionode = ( mpime == root )
!
      IF( ionode ) &
        WRITE( 6, * ) 'Using Layers Mesh decomposition' 
!
      ALLOCATE(lay_map(0:nproc-1,2))
!
      ncfl1 = 0
      nctot = 0
!
! ... distribute to processors an equal number
! ... of cell of weight ONE
!
      DO ipe = 0, nproc - 1
       ncfl1(ipe) = localdim(countfl, nproc, ipe)
      END DO
!
      IF( job_type == '2D' ) THEN

        lay_map(0      ,1) = 1        !  the first cell to the first proc
        lay_map(nproc-1,2) = nx*nz    !  the last cell to the last proc
        ipe      = 0
        icnt     = 0
        icnt_ipe = 0

        DO k = 1, nz
          DO i = 1, nx

            ijk = i + (k-1) * nx
  
            IF ( BTEST(fl(ijk),0) ) THEN
              icnt_ipe = icnt_ipe + 1
              icnt     = icnt     + 1
            END IF

            nctot(ipe) = nctot(ipe) + 1

            IF ( ( icnt_ipe >= ncfl1(ipe) ) .AND. ( i /= (nx-1) ) ) THEN

! ...         The last cell has already been assigned to the last proc,
! ...         the next double condition is somehow redundant, but helps reading the code
              IF( ( icnt < countfl ) .AND. ( ipe < ( nproc - 1 ) ) ) THEN
                lay_map(ipe,2) = ijk
                ipe = ipe + 1
                lay_map(ipe,1) = ijk+1
              END IF

              icnt_ipe = 0

            END IF

          END DO 
        END DO

      ELSE IF( job_type == '3D' ) THEN

        lay_map(0,1)       = 1         ! the first cell to the first proc
        lay_map(nproc-1,2) = nx*ny*nz  ! the last cell to the last proc
        ipe      = 0
        icnt     = 0
        icnt_ipe = 0

        DO k = 1, nz
          DO j = 1, ny
            DO i = 1, nx

              ijk = i + (j-1)*nx + (k-1)*nx*ny 
  
              IF ( BTEST(fl(ijk),0) ) THEN
                icnt_ipe = icnt_ipe + 1
                icnt     = icnt     + 1
              END IF

              nctot(ipe) = nctot(ipe) + 1

              IF ( icnt_ipe >= ncfl1(ipe) ) THEN
                IF ( ( (i /= nx-1) .AND. (j <= ny-1) ) .OR. ( (i == nx) .AND. (j == ny) ) ) THEN
                  IF( icnt < countfl .AND. ( ipe < ( nproc - 1) ) ) THEN
                    lay_map(ipe,2) = ijk
                    ipe = ipe + 1
                    lay_map(ipe,1) = ijk+1
                  END IF
                  icnt_ipe = 0
                END IF
              END IF

            END DO
          END DO
        END DO

      ELSE

        CALL error(' layers ',' unknow job_type ', 1)

      END IF

      DO ipe = 0, nproc - 1
        proc_map(ipe)%lay(:) =  lay_map(ipe,:)
      END DO

      DEALLOCATE(lay_map)
      
!
      RETURN
      END SUBROUTINE layers

! ------------------------------------------------------------------------

      SUBROUTINE blocks( ncfl1, nctot, nproc, mpime, root )
!
! ... partition N.1 (layers)
! ... decomposes the domain into N rectangular blocks.
! ... The processor map consist of the coordinate pairs
! ... of the bottom-west and the top-east corners of the blocks
! ... (ONLY 2D)
!
! ... partition N.2 (blocks)
!
      USE dimensions
!
      IMPLICIT NONE
      SAVE
!
      INTEGER :: ncfl1(0:)
      INTEGER :: nctot(0:)
      INTEGER, INTENT(IN) :: nproc, mpime, root
!
      INTEGER :: i, ib, k, ijk
      INTEGER :: l
      INTEGER :: i1, i2, k1, k2, ijk1, ijk2
      INTEGER :: nbl
      INTEGER :: icnt, ipe
      INTEGER :: bl1
      INTEGER :: size_x, rest, skipl, rrest
      INTEGER :: nbx, nbz
      INTEGER :: layer, icnt_layer, here
      INTEGER :: localdim
!
      INTEGER, ALLOCATABLE :: lay_map(:,:)
      INTEGER, ALLOCATABLE :: nctot_lay(:), nbl_lay(:), ncfl1_lay(:), nbl_tot(:)
!
      REAL*8 side, area
      REAL*8 fact
      LOGICAL :: ionode

      ionode = ( mpime == root )

      IF( ionode ) &
        WRITE( 6, * ) 'Using Blocks ( 2D ) Mesh decomposition' 
!
! ... Initially subdivides the domain into 'nbz' layers. Each layer
! ... will be in turn subdivided into 'nbx' blocks.
! ... If the number of the processors 'nproc' is not factorable
! ... as 'nbx*nbz', some blocks in the lower layer will be merged together.
!
        area  = DBLE(countfl/nproc)       

! ... Here countfl is supposed to be the total area whose side
! ... are "nx" and "nz_effective" ( lower than "nz" )
! ... Remember nbx * nbz should be equal to nproc

        side  = DSQRT(area)    !  balanced area to be assigned to each proc
        nbx   = NINT(nx/side)  !  estimate number of block along x direction
        nbz   = NINT(countfl/nbx/area) ! estimate the number of block along z    
        IF (nbx == 0) nbx = 1
        IF (nbz == 0) nbz = 1
        DO WHILE (nbx*nbz <= nproc)
         IF (INT(nx/nbx) .GT. INT(nz/nbz)) THEN
           nbx = nbx + 1
         ELSE
           nbz = nbz + 1
         END IF
        END DO

! ... Guess the optimum 'size_x' of the blocks.

        size_x = INT(nx/nbx)
        rest   = nbx*nbz - nproc
!
! ... nbz is now computed, and will not be changed anymore
! ... we need now to compute the border of the layers and the
! ... number of block within each layer.
!
        ALLOCATE(lay_map(1:nbz,2))
        lay_map(:,:) = 0
!
! ... Initialize the processor maps
!
        DO i = 0, nproc-1
          proc_map(i)%corner1(:) = 0
          proc_map(i)%corner2(:) = 0
        END DO
!
        ALLOCATE(nctot_lay(1:nbz), ncfl1_lay(1:nbz), nbl_lay(1:nbz))
!
        nctot_lay = 0
        ncfl1_lay = 0
        nbl_lay   = 0
!
! ... distribute cells among layers proportionally to the 
! ... number of blocks contained. 
! ... 'ncfl1_lay' is the number of cells with fl=1 in the layer.
!
        IF (rest == 0) THEN

! ...     Every layer contain the same number of blocks ( no reminder )

          DO layer = 1, nbz
            ncfl1_lay(layer) = localdim( countfl, nbz, layer-1 )
            nbl_lay(layer) = nbx
          END DO

        ELSE

! ...     Too many blocks, some layer should remove the eceeding blocks ( "rest" )

          skipl = INT(rest/nbx) 
          rrest = MOD(rest,nbx)

          IF (rest < nbx) THEN
            nbl_lay(1) = nbx - rrest   !  the first layer (usually with topography )
            DO layer = 2, nbz          !  remove all eceeding blocks
              nbl_lay(layer) = nbx
            END DO
          ELSE
            DO layer = 1, skipl        !  rest is larger than nbx, remove all possible 
              nbl_lay(layer) = 1       !  blocks from the lower layers
            END DO
            IF( skipl < nbz ) THEN     !  this condition is somehow redundant
              nbl_lay(skipl+1) = nbx - rrest  !  remove the remaining blocks
              DO layer = skipl + 2, nbz
                nbl_lay(layer) = nbx
              END DO
            END IF
          END IF

          IF( SUM( nbl_lay ) /= nproc ) &
            CALL error(' blocks ',' non consistent number of blocks ', SUM( nbl_lay ) )

! ...     Distribute cells to each layer, the number of cells is proportional to
! ...     the number of blocks

          ALLOCATE( nbl_tot( nproc ) )

          DO ib = 1, nproc
            nbl_tot(ib) = localdim( countfl, nproc, ib-1 )
          END DO
          ib = 0
          DO layer = 1, nbz
            ncfl1_lay(layer) = 0
            DO i = 1, nbl_lay(layer)
              ib = ib + 1
              ncfl1_lay(layer) = ncfl1_lay(layer) + nbl_tot( ib )
            END DO
          END DO

          DEALLOCATE( nbl_tot )

        END IF

        IF ( SUM( ncfl1_lay(:) ) /= countfl ) &
              CALL error(' blocks ',' control blocks decomposition ',1)
!
! ... build the layer maps
!
        layer = 1
        ipe   = 0
        lay_map(1  ,1) = 1
        lay_map(nbz,2) = nx*nz
        icnt_layer = 0
        icnt       = 0

        DO k = 1, nz
          DO i = 1, nx
            ijk = i + (k-1)*nx
!
            IF ( BTEST(fl(ijk),0) ) THEN
              icnt_layer = icnt_layer + 1
              icnt       = icnt     + 1
            END IF

            nctot_lay(layer) = nctot_lay(layer) + 1

            IF ( icnt_layer >= ncfl1_lay(layer) ) THEN
              IF( layer < nbz ) THEN
                lay_map(layer,2) = ijk
                layer = layer + 1
                lay_map(layer,1) = ijk+1
              END IF
              icnt_layer = 0
            END IF
          END DO
        END DO
!
      ipe = -1
      DO layer = 1, nbz

        ijk2 = lay_map(layer,2)
        k2   = ( ijk2 - 1 ) / nx + 1
        i2   = MOD( ( ijk2 - 1 ), nx) + 1
!
! ... cut steps to layers
!
        IF ( i2 < nx/2 ) k2 = k2 - 1

        IF ( layer < nbz ) THEN
          lay_map(layer,2) = k2 * nx
        ELSE
          lay_map(layer,2) = nz * nx
        END IF
        IF (layer > 1) THEN 
          lay_map(layer,1) = lay_map(layer-1,2) + 1
        ELSE
          lay_map(layer,1) = 1
        END IF

! ... retieve layer boundaries and check them

        ijk1 = lay_map(layer,1)
        ijk2 = lay_map(layer,2)
        k1  = ( ijk1 - 1 ) / nx + 1
        i1  = MOD( ( ijk1 - 1 ), nx) + 1
        k2  = ( ijk2 - 1 ) / nx + 1
        i2  = MOD( ( ijk2 - 1 ), nx) + 1
        IF (i1 /= 1 .OR. i2 /= nx) &
          CALL error( ' blocks ', ' inconsistent layer ', layer )
!
! ...  now updates the number of cells with fl=1 
!
        nctot_lay(layer) = 0
        ncfl1_lay(layer) = 0
        DO ijk = ijk1, ijk2
          nctot_lay(layer) = nctot_lay(layer) + 1
          IF ( BTEST(fl(ijk),0) ) ncfl1_lay(layer) = ncfl1_lay(layer) + 1 
        END DO
!
! ... Estimate of the number of cells contained in each block
! ... in the layer
!
        bl1 = INT( ncfl1_lay(layer) / nbl_lay(layer) )
!
! ... Build block maps
!
        i = i1
        DO nbl = 1, nbl_lay(layer) 

          ipe = ipe + 1

          proc_map(ipe)%corner1(1) = i
          proc_map(ipe)%corner1(2) = k1

          DO WHILE ( ncfl1(ipe) < bl1 )
            DO k = k1, k2
              ijk = i + (k-1)*nx
              IF ( BTEST(fl(ijk),0) ) ncfl1(ipe) = ncfl1(ipe) + 1
            END DO  
            i = i + 1
          END DO

          proc_map(ipe)%corner2(1) = i - 1
          proc_map(ipe)%corner2(2) = k2

          IF (nbl == nbl_lay(layer)) proc_map(ipe)%corner2(1) = nx 
!
! ... Now updates the number of cells with fl=1 into blocks
!
          nctot(ipe)  = 0
          ncfl1(ipe) = 0
          DO k = proc_map(ipe)%corner1(2), proc_map(ipe)%corner2(2)
            DO i = proc_map(ipe)%corner1(1), proc_map(ipe)%corner2(1)
              ijk = i + (k-1)*nx
              nctot(ipe) = nctot(ipe) + 1
              IF ( BTEST(fl(ijk),0) ) ncfl1(ipe) = ncfl1(ipe) + 1
            END DO
          END DO
        END DO
      END DO

      DEALLOCATE(lay_map)
      DEALLOCATE( nctot_lay, ncfl1_lay, nbl_lay )

      RETURN
      END SUBROUTINE blocks

! ---------------------------------------------------------------------

      SUBROUTINE columns( ncfl1, nctot, nproc, mpime, root )
!
! ... decomposes the domain into N rectangular blocks.
! ... The processor map consist of the coordinate pairs
! ... of the bottom-west and the top-east corners of the blocks
! ... (ONLY 3D)
!
!
      USE dimensions
!
      IMPLICIT NONE
      SAVE
!
      INTEGER :: ncfl1(0:)
      INTEGER :: nctot(0:)
      INTEGER, INTENT(IN) :: nproc, mpime, root
!
      INTEGER :: i, j, k, ijk
      INTEGER :: i1, i2, j1, j2, ij1, ij2, ij
      INTEGER :: nbl, l
      INTEGER :: icnt, ipe
      INTEGER :: bl1
      INTEGER :: rest, skipl, rrest
      INTEGER :: nbx, nby
      INTEGER :: layer, icnt_layer, here
      INTEGER :: localdim
!
      INTEGER, ALLOCATABLE :: lay_map(:,:)
      INTEGER, ALLOCATABLE :: nctot_lay(:), nbl_lay(:), ncfl1_lay(:)
!
      REAL*8 side, area, volume, zavg
      REAL*8 fact
      LOGICAL :: ionode

      ionode = ( mpime == root )

      IF( ionode ) &
        WRITE( 6, * ) 'Using Colomns ( 3D ) Mesh decomposition' 

!
! ... subdivides the domain into 'nby' Slabs. Each slab
! ... will be in turn subdivided into 'nbx' column block.
! ... If the number of the processors 'nproc' is not factorable
! ... as 'nbx*nby', some blocks in the lower slab will be merged together.
!
! ... countfl represent the number of computationally demanding
! ... units of volume, that should be equally distributed across 
! ... processors.

        zavg = DBLE( countfl) / ( nx * ny )    ! average no. of counts along z dim
        volume = DBLE( countfl ) / DBLE( nproc )  ! per processor no. of counts

        area  = volume / zavg       ! expected average base area of column blocks 
        side  = DSQRT( area )       ! side of the base area
        nbx   = NINT( nx / side )   ! per processor number of "x"
        nby   = NINT( DBLE( nproc ) / nbx )
        IF (nbx == 0) nbx = 1
        IF (nby == 0) nby = 1

        DO WHILE ( nbx*nby < nproc )
         IF ( INT( nx / nbx ) .GT. INT( ny / nby ) ) THEN
           nbx = nbx + 1
         ELSE
           nby = nby + 1
         END IF
        END DO
        rest = nbx*nby - nproc

!
        ALLOCATE( lay_map( 1:nby, 2 ) )
        lay_map(:,:) = 0
!
! ... Initialize the processor maps
!
        DO i = 0, nproc-1
          proc_map(i)%corner1(:) = 0
          proc_map(i)%corner2(:) = 0
        END DO
!
        ALLOCATE( nctot_lay( 1:nby ), ncfl1_lay( 1:nby ), nbl_lay( 1:nby ) )
!
        nctot_lay = 0
        ncfl1_lay = 0
        nbl_lay   = 0
!
! ... distribute cells among layers proportionally to the 
! ... number of blocks contained. 
! ... 'ncfl1_lay' is the number of cells with fl=1 in the layer.
!
        IF ( rest == 0 ) THEN

          DO layer = 1, nby
            ncfl1_lay(layer) = localdim( countfl, nby, layer - 1 )
            nbl_lay(layer)   = nbx
          END DO

        ELSE IF (rest > 0) THEN

          IF( nbx == 1 ) THEN
            CALL error(' partition ',' (nbx == 1) AND (rest /= 0) ??? ', rest )
          END IF

          IF( rest < nby ) THEN
            DO layer = 1, rest
              nbl_lay(layer) = nbx - 1 
            END DO
            DO layer = rest + 1, nby
              nbl_lay(layer) = nbx
            END DO
          ELSE
            CALL error(' partition ',' rest >= nby ??? ', rest )
          END IF

          DO layer = 1, nby
            fact = DBLE( countfl / nproc * nbl_lay(layer) )
            ncfl1_lay(layer) = NINT(fact)   !  expected no. of fl=1 per slab (layer 2D)
          END DO

        END IF
!
! ... build the layer maps
!
        layer = 1
        ipe = 0
        lay_map(1,1)   = 1
        lay_map(nby,2) = nx*ny
        icnt_layer = 0
        icnt       = 0

        DO j = 1, ny
        DO i = 1, nx

          DO k = 1, nz
            ijk = i + (j-1)*nx + (k-1)*nx*ny
!
            IF ( BTEST(fl(ijk),0) ) THEN
              icnt_layer = icnt_layer + 1
              icnt       = icnt       + 1
            END IF

            nctot_lay(layer) = nctot_lay(layer) + 1

            IF ( icnt_layer == ncfl1_lay(layer) ) THEN
              icnt_layer = 0
              IF( layer < nby ) THEN
                lay_map(layer,2) =  i + (j-1)*nx
                layer = layer + 1
                lay_map(layer,1) =  i + (j-1)*nx + 1
              END IF
            END IF

          END DO

        END DO
        END DO
!
        ipe = -1
        DO layer = 1, nby

          ij2 = lay_map(layer,2)
          j2  = ( ij2 - 1 ) / nx + 1
          i2  = MOD( ( ij2 - 1 ), nx) + 1
!
! ... cut steps to layers
!
          IF ( i2 < nx/2 ) j2 = j2-1

          IF (layer < nby) THEN
            lay_map(layer,2) = j2 * nx
          ELSE
            lay_map(layer,2) = ny * nx
          END IF

          IF (layer > 1) THEN 
            lay_map( layer, 1 ) = lay_map( layer - 1, 2 ) + 1
          ELSE
            lay_map( layer, 1 ) = 1
          END IF

          ij1 = lay_map(layer,1)
          ij2 = lay_map(layer,2)

          j1  = ( ij1 - 1 ) / nx + 1
          i1  = MOD( ( ij1 - 1 ), nx) + 1

          j2  = ( ij2 - 1 ) / nx + 1
          i2  = MOD( ( ij2 - 1 ), nx) + 1

          IF (i1 /= 1 .OR. i2 /= nx) THEN
            CALL error( ' blocks3D ', ' error in layer ', layer )
          END IF
!
! ...  now updates the number of cells with fl=1 into layers
!
          nctot_lay(layer) = 0
          ncfl1_lay(layer) = 0
          DO ij = ij1, ij2
            DO k = 1, nz
              ijk = ij + (k-1)*nx*ny
              nctot_lay(layer) = nctot_lay(layer) + 1
              IF ( BTEST(fl(ijk),0) ) ncfl1_lay(layer) = ncfl1_lay(layer) + 1 
            END DO
          END DO

          IF( lpr > 0 .AND. ionode ) WRITE(logunit,*) ' layer = ', layer, ' no. columns = ', nbl_lay(layer)
!
! ... Estimate of the number of cells contained in each block
! ... in the layer
!
          bl1 = INT( ncfl1_lay(layer) / nbl_lay(layer) )
!
! ... Build block maps
!
          i = i1
          DO nbl = 1, nbl_lay(layer) 
            ipe = ipe + 1

            proc_map(ipe)%corner1(1) = i
            proc_map(ipe)%corner1(2) = j1

            DO WHILE ( ncfl1(ipe) < bl1 )
              DO j = j1, j2
                DO k = 1, nz
                  ijk = i + (j-1)*nx + (k-1)*nx*ny
                  IF ( BTEST(fl(ijk),0) ) ncfl1(ipe) = ncfl1(ipe) + 1
                END DO
              END DO  
              i = i + 1
            END DO

            proc_map(ipe)%corner2(1) = i - 1
            proc_map(ipe)%corner2(2) = j2

            IF ( nbl == nbl_lay(layer) ) proc_map(ipe)%corner2(1) = nx 
!
! ... Now updates the number of cells with fl=1 into blocks
!
            nctot(ipe) = 0
            ncfl1(ipe) = 0
            DO k = 1, nz
              DO j = proc_map(ipe)%corner1(2), proc_map(ipe)%corner2(2)
                DO i = proc_map(ipe)%corner1(1), proc_map(ipe)%corner2(1)
                  ijk = i + (j-1)*nx + (k-1)*nx*ny
                  nctot(ipe) = nctot(ipe) + 1
                  IF ( BTEST(fl(ijk),0) ) ncfl1(ipe) = ncfl1(ipe) + 1
                END DO
              END DO
            END DO
            IF (lpr > 1 .AND. ionode ) &
              WRITE(logunit,*) 'proc_map(',ipe,'): ', proc_map(ipe)%corner1(:), proc_map(ipe)%corner2(:)
          END DO
        END DO

        DEALLOCATE( lay_map )
        DEALLOCATE( nctot_lay, ncfl1_lay, nbl_lay )

      RETURN
      END SUBROUTINE columns
! ---------------------------------------------------------------------

      SUBROUTINE blocks3d( ncfl1, nctot, nproc, mpime, root )
!
! ... This subroutine decomposes the domain into N orthorombic blocks.
! ... The processor map consist of the cell coordinates ( i, j, k )
! ... of the bottom-south-west and the top-north-east corners of the blocks
! ... itself
!
!
      USE dimensions
!
      IMPLICIT NONE
      SAVE
!
      INTEGER :: ncfl1(0:)   !  the number of cells with fl == 1 
      INTEGER :: nctot(0:)   !  total number of local cells
      INTEGER, INTENT(IN) :: nproc, mpime, root
!
      INTEGER :: i, j, k, ijk
      INTEGER :: i1, i2, j1, j2, ij1, ij2, ij, k1, k2
      INTEGER :: nbl, l, ibz
      INTEGER :: icnt, ipe
      INTEGER :: rest, skipl, rrest
      INTEGER :: nbx, nby, nbz
      INTEGER :: layer, icnt_layer, here
      INTEGER :: icnt_column
      INTEGER :: nprocxy, nprocz

      REAL*8  :: bl1, cl1
!
      INTEGER :: localdim
!
      INTEGER, ALLOCATABLE :: lay_map(:,:)  ! layer map used as first domain partition
      INTEGER, ALLOCATABLE :: nctot_lay(:)  ! total number of cell per layer
      INTEGER, ALLOCATABLE :: nbl_lay(:)    ! number of column blocks within a layer
      INTEGER, ALLOCATABLE :: ncfl1_lay(:)  ! total number of cells in the layer with fl == 1
!
      REAL*8 :: side, area, volume, zavg, fact
      LOGICAL :: ionode

! ... Subroutine Body

      ionode = ( mpime == root )

      IF( ionode ) &
        WRITE( 6, * ) 'Using Blocks ( 3D ) Grid partition' 
!
! ... Factorize the number of processors, into plane processors nprocxy
! ... and z processors nprocz
!
      nprocz = NINT( DBLE( nproc )**(1.0d0/3.0d0) )   !  first guess

! ... ensure that nprocz is a factor of nproc
! ... we take the first factor lower than the first guess

      DO i = nprocz, 1, -1
        IF( MOD( nproc, i ) == 0 ) EXIT
      END DO 
      nprocz = i

! ... compute the number of processor in the xy plane

      nprocxy = nproc / nprocz

      IF( nprocz == 1 .AND. ionode ) THEN
        WRITE(errorunit,*) 'WARNING! from domain decomposition'
        WRITE(errorunit,*) 'Blocks3d distribution has no effect, nprocz = 1'
      END IF

! ... compute blocks in z direction 

      nbz = nprocz  

!
! ... subdivides the XY domain into 'nby' Slabs. Each slab
! ... will be in turn subdivided into 'nbx' column blocks.
! ... If the number of the processors 'nprocxy' is not factorable
! ... as 'nbx*nby', some column blocks in the lower slab will be merged together.
!
! ... countfl represent the number of computationally demanding
! ... units of volume, that should be equally distributed across 
! ... processors.
!

        zavg = DBLE( countfl) / ( nx * ny )         ! average no. of counts along z dim
        volume = DBLE( countfl ) / DBLE( nprocxy )  ! per (columns) processor no. of counts

        area  = volume / zavg       ! expected average base area of column blocks 
        side  = DSQRT( area )       ! side of the base area
        nbx   = NINT( nx / side )   ! per processor number of column blocks along "x" 
        nby   = NINT( DBLE( nprocxy ) / nbx )
        IF (nbx == 0) nbx = 1
        IF (nby == 0) nby = 1

        DO WHILE ( nbx*nby < nprocxy )
         IF ( INT( nx / nbx ) .GT. INT( ny / nby ) ) THEN
           nbx = nbx + 1
         ELSE
           nby = nby + 1
         END IF
        END DO
        rest = nbx*nby - nprocxy

        IF( lpr > 1 .AND. ionode ) THEN
          WRITE(logunit,*) 'Report on domain decomposition'
          WRITE(logunit,*) 'nbx, nby, nbz, nprocz, nprocxy = ', nbx, nby, nbz, nprocz, nprocxy
        END IF

! ... The number of layer (slab) nby is now defined, and it will not be
! ... modified further

        ALLOCATE( lay_map( 1:nby, 2 ) )
        ALLOCATE( nctot_lay( 1:nby ), ncfl1_lay( 1:nby ), nbl_lay( 1:nby ) )
        lay_map   = 0
        nctot_lay = 0
        ncfl1_lay = 0
        nbl_lay   = 0
!
! ... Initialize the processor maps
!
        DO i = 0, nproc-1
          proc_map(i)%blkbsw(:) = 0
          proc_map(i)%blktne(:) = 0
        END DO
!
! ... distribute cells among layers proportionally to the 
! ... number of blocks contained. 
! ... 'ncfl1_lay' is the number of cells with fl=1 in the layer.
!
        IF ( rest == 0 ) THEN

          DO layer = 1, nby
            ncfl1_lay(layer) = localdim( countfl, nby, layer - 1 )
            nbl_lay(layer)   = nbx
          END DO

        ELSE IF ( rest > 0 ) THEN

          IF( nbx == 1 ) THEN
            CALL error(' partition ',' (nbx == 1) AND (rest /= 0) ??? ', rest )
          END IF

! ...     since in 3D we partition the domain partitioning the XY plane,
! ...     there is no particular advantages to distinguish lower or higher
! ...     layers

          IF( rest < nby ) THEN
            DO layer = 1, rest
              nbl_lay(layer) = nbx - 1 
            END DO
            DO layer = rest + 1, nby
              nbl_lay(layer) = nbx
            END DO
          ELSE
            CALL error(' partition ',' rest >= nby ??? ', rest )
          END IF

          DO layer = 1, nby
            fact = DBLE( countfl/nprocxy * nbl_lay(layer) )
            ncfl1_lay(layer) = NINT(fact)   !  expected no. of fl=1 per slab (layer 2D)
          END DO

        END IF

        IF( lpr > 1 .AND. ionode ) THEN
          DO layer = 1, nby
            WRITE(logunit,*) ' layer ', layer, ' nbl = ', nbl_lay(layer)
          END DO
        END IF

!
! ... build the layers map
!
        layer = 1
        lay_map(1,1)   = 1
        lay_map(nby,2) = nx*ny
        icnt_layer = 0
        icnt       = 0

        DO j = 1, ny
          DO i = 1, nx

            DO k = 1, nz

              ijk = i + (j-1)*nx + (k-1)*nx*ny
!
              IF ( BTEST(fl(ijk),0) ) THEN
                icnt_layer = icnt_layer + 1
                icnt       = icnt       + 1
              END IF

              nctot_lay(layer) = nctot_lay(layer) + 1

              IF ( icnt_layer >= ncfl1_lay(layer) ) THEN
                IF( layer < nby ) THEN
                  lay_map(layer,2) =  i + (j-1)*nx
                  layer = layer + 1
                  lay_map(layer,1) =  i + (j-1)*nx + 1
                END IF
                icnt_layer = 0
              END IF

            END DO

          END DO
        END DO

        IF( lpr > 1 .AND. ionode ) THEN
          DO layer = 1, nby
            WRITE(logunit,*) ' layer_map ', layer, ' start, end = ', lay_map(layer,1), lay_map(layer,2)
          END DO
        END IF

! ... build the blocks map

        ipe = -1
        LAYERS: DO layer = 1, nby

! ...     now set y coordinate of the block corners ( j1, j2 )
! ...     using the layer map 
!
          ij2 = lay_map(layer,2)
          j2  = ( ij2 - 1 ) / nx + 1
          i2  = MOD( ( ij2 - 1 ), nx) + 1
!
! ...     cut steps to layers
!
          IF ( i2 < nx/2 ) j2 = j2-1

          IF (layer < nby) THEN
            lay_map(layer,2) = j2 * nx
          ELSE
            lay_map(layer,2) = ny * nx
          END IF

          IF (layer > 1) THEN 
            lay_map( layer, 1 ) = lay_map( layer - 1, 2 ) + 1
          ELSE
            lay_map( layer, 1 ) = 1
          END IF

          ij1 = lay_map(layer,1)
          ij2 = lay_map(layer,2)

          j1  = ( ij1 - 1 ) / nx + 1
          i1  = MOD( ( ij1 - 1 ), nx) + 1

          j2  = ( ij2 - 1 ) / nx + 1
          i2  = MOD( ( ij2 - 1 ), nx) + 1

          IF (i1 /= 1 .OR. i2 /= nx) THEN
            CALL error( ' blocks3D ', ' error in layer ', layer )
          END IF
!
! ...     now updates the number of cells with fl=1 into layers
!
          nctot_lay(layer) = 0
          ncfl1_lay(layer) = 0
          DO ij = ij1, ij2
            DO k = 1, nz
              ijk = ij + (k-1)*nx*ny
              nctot_lay(layer) = nctot_lay(layer) + 1
              IF ( BTEST(fl(ijk),0) ) ncfl1_lay(layer) = ncfl1_lay(layer) + 1 
            END DO
          END DO

          IF( lpr > 1 .AND. ionode ) &
            WRITE(logunit,*) ' layer = ', layer, ' n. column = ', nbl_lay(layer), ' n. of fl1 = ', ncfl1_lay(layer)
!
! ...     Estimate of the number of cells contained in each column block
! ...     in the layer
!
          cl1 = DBLE( ncfl1_lay(layer) ) / nbl_lay(layer)

! ...     and now estimate the number of cells in each orthorombic block

          bl1 = cl1 / nbz

          i1  = 1
!
          COLUMNS: DO nbl = 1, nbl_lay(layer) 

! ...       now set x coordinate of the block corners ( i1, i2 )
! ...       these define the column extension

            icnt_column = 0

            IF( nbl == nbl_lay(layer) ) THEN
              i2 = nx
            ELSE
              LOOPX: DO i = i1, nx
                DO j = j1, j2
                  DO k = 1, nz
                    ijk = i + (j-1)*nx + (k-1)*nx*ny
                    IF ( BTEST(fl(ijk),0) ) icnt_column = icnt_column + 1
                    IF ( icnt_column > cl1 ) THEN
                      IF( DBLE( j - j1 + 1 ) >= 0.5d0 * DBLE( j2 - j1 + 1 ) ) THEN
                        i2 = i
                      ELSE
                        i2 = i - 1
                      END IF
                      EXIT LOOPX
                    END IF
                  END DO
                END DO
              END DO LOOPX
            END IF



! ...       now set z coordinate (balancing the load) of the block corners ( k1, k2 )
! ...       These together with i1,i2,j1,j2 define the corners of the orthorombic
! ...       block

            k1 = 1

            BLOCKS: DO ibz = 1, nbz

              ipe = ipe + 1

              IF( ipe >= nproc ) &
                CALL error( ' blocks3d ', ' inconsistent proc index ', ipe )
 
              proc_map(ipe)%blkbsw(1) = i1
              proc_map(ipe)%blkbsw(2) = j1
              proc_map(ipe)%blkbsw(3) = k1

              IF( ibz < nbz ) THEN

                ncfl1(ipe) = 0

                LOOPZ: DO k = k1, nz
                  DO j = j1, j2
                    DO i = i1, i2
                      ijk = i + (j-1)*nx + (k-1)*nx*ny
                      IF ( BTEST(fl(ijk),0) ) ncfl1(ipe) = ncfl1(ipe) + 1
                      IF ( ncfl1( ipe ) > bl1 ) THEN
                        k2 = k
                        !IF( ( j - j1 + 1 ) > ( j2 - j1 + 1 ) / 2 ) THEN
                        !  k2 = k
                        !ELSE
                        !  k2 = k - 1
                        !END IF
                        EXIT LOOPZ
                      END IF
                    END DO
                  END DO  
                END DO LOOPZ


                proc_map(ipe)%blktne(1) = i2
                proc_map(ipe)%blktne(2) = j2
                proc_map(ipe)%blktne(3) = k2

              ELSE

                proc_map(ipe)%blktne(1) = i2
                proc_map(ipe)%blktne(2) = j2
                proc_map(ipe)%blktne(3) = nz 

              END IF

!
! ...         Now updates the number of cells with fl=1 into blocks
!
              nctot(ipe) = 0
              ncfl1(ipe) = 0
              DO k = proc_map(ipe)%blkbsw(3), proc_map(ipe)%blktne(3)
                DO j = proc_map(ipe)%blkbsw(2), proc_map(ipe)%blktne(2)
                  DO i = proc_map(ipe)%blkbsw(1), proc_map(ipe)%blktne(1)
                    ijk = i + (j-1)*nx + (k-1)*nx*ny
                    nctot(ipe) = nctot(ipe) + 1
                    IF ( BTEST(fl(ijk),0) ) ncfl1(ipe) = ncfl1(ipe) + 1
                  END DO
                END DO
              END DO

              IF (lpr > 1 .AND. ionode ) THEN
                 WRITE(logunit,1000) ipe, proc_map(ipe)%blkbsw(:),proc_map(ipe)%blktne(:)
  1000           FORMAT( ' proc_map(',I4,'): BSW = ', 3I4, ' TNE = ', 3I4 )
              END IF

              k1 = k2 + 1


            END DO BLOCKS

            i1 = i2 + 1

          END DO COLUMNS

        END DO LAYERS

        IF( SUM( ncfl1(0:nproc-1) ) /= countfl ) &
          CALL error(' blocks3d ', ' inconsistent number of cell ', SUM( ncfl1(0:nproc-1) ) )

        DEALLOCATE( lay_map, nctot_lay, ncfl1_lay, nbl_lay )

      RETURN
      END SUBROUTINE blocks3d
! ---------------------------------------------------------------------
      SUBROUTINE ghost
!
! ... Identifies and allocates ghost cells
! ... (2D/3D-Compliant)
! 
      USE dimensions
      USE parallel, ONLY: nproc, mpime, root
      USE basic_types, ONLY: imatrix
      USE control_flags, ONLY: job_type
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
      INTEGER :: layer, k2, k1, j2, j1, i2, i1, nkt
      INTEGER :: me, whose
!
      IF(ALLOCATED(rcv_map)) DEALLOCATE(rcv_map)
      IF(ALLOCATED(snd_map)) DEALLOCATE(snd_map)
      ALLOCATE(rcv_map( 0:(nproc-1) ) )
      ALLOCATE(snd_map( 0:(nproc-1) ) )
      
!
      WRITE( 7, * ) 'Entering Ghost ... '
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
!
      nset  = 0
      ncext = 0

      IF ( proc_map(mpime)%type == LAYER_MAP ) THEN
        DO ijk = proc_map(mpime)%lay(1), proc_map(mpime)%lay(2)
          IF ( BTEST(fl(ijk),0) .OR. BTEST(fl(ijk),9) ) THEN
            ncext = ncext + cell_neighbours( ijk, mpime, nset)
          END IF
        END DO
      ELSE IF ( proc_map(mpime)%type == BLOCK2D_MAP ) THEN
        IF( job_type == '2D' ) THEN
          i1 = proc_map(mpime)%corner1(1)
          i2 = proc_map(mpime)%corner2(1)
          k1 = proc_map(mpime)%corner1(2)
          k2 = proc_map(mpime)%corner2(2)
          DO k = k1, k2
            DO i = i1, i2
              ijk = i + (k-1) * nx
              IF ( BTEST(fl(ijk),0) .OR. BTEST(fl(ijk),9) ) THEN
                ncext = ncext + cell_neighbours(ijk, mpime, nset)
              END IF
            END DO
          END DO
        ! ... "columns" domain-decomposition ... !
        ELSE IF( job_type == '3D' ) THEN
          i1 = proc_map(mpime)%corner1(1)
          i2 = proc_map(mpime)%corner2(1)
          j1 = proc_map(mpime)%corner1(2)
          j2 = proc_map(mpime)%corner2(2)
          DO k = 1, nz
            DO j = j1, j2
              DO i = i1, i2
                ijk = i + (j-1)*nx + (k-1)*nx*ny 
                IF ( BTEST(fl(ijk),0) .OR. BTEST(fl(ijk),9) ) THEN
                  ncext = ncext + cell_neighbours(ijk, mpime, nset)
                END IF
              END DO
            END DO
          END DO
        ELSE 
          CALL error( ' ghost ', ' wrong job_type '//job_type, 1)
        END IF
      ELSE IF ( proc_map(mpime)%type == BLOCK3D_MAP ) THEN
        IF( job_type == '3D' ) THEN
          i1 = proc_map(mpime)%blkbsw(1)
          i2 = proc_map(mpime)%blktne(1)
          j1 = proc_map(mpime)%blkbsw(2)
          j2 = proc_map(mpime)%blktne(2)
          k1 = proc_map(mpime)%blkbsw(3)
          k2 = proc_map(mpime)%blktne(3)
          DO k = k1, k2
            DO j = j1, j2
              DO i = i1, i2
                ijk = i + (j-1)*nx + (k-1)*nx*ny
                IF ( BTEST(fl(ijk),0) .OR. BTEST(fl(ijk),9) ) THEN
                  ncext = ncext + cell_neighbours(ijk, mpime, nset)
                END IF
              END DO
            END DO
          END DO
        ELSE 
          CALL error( ' ghost ', ' wrong job_type '//job_type, 1)
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
          IF ( BTEST(fl(ijk),0) .OR. BTEST(fl(ijk),9) ) THEN
            icnt = icnt + cell_neighbours(ijk, mpime, nset, rcv_cell_set, myijk)
          END IF
        END DO
      ELSE IF ( proc_map(mpime)%type == BLOCK2D_MAP ) THEN
        IF( job_type == '2D' ) THEN
          i1 = proc_map(mpime)%corner1(1)
          i2 = proc_map(mpime)%corner2(1)
          k1 = proc_map(mpime)%corner1(2)
          k2 = proc_map(mpime)%corner2(2)
          DO k = k1, k2
            DO i = i1, i2
              ijk = i + (k-1) * nx
              IF ( BTEST(fl(ijk),0) .OR. BTEST(fl(ijk),9) ) THEN
                icnt = icnt + cell_neighbours(ijk, mpime, nset, rcv_cell_set, myijk)
              END IF
            END DO
          END DO
        ELSE IF( job_type == '3D' ) THEN
          i1 = proc_map(mpime)%corner1(1)
          i2 = proc_map(mpime)%corner2(1)
          j1 = proc_map(mpime)%corner1(2)
          j2 = proc_map(mpime)%corner2(2)
          DO k = 1, nz
            DO j = j1, j2
              DO i = i1, i2
                ijk = i + (j-1)*nx + (k-1)*nx*ny
                IF ( BTEST(fl(ijk),0) .OR. BTEST(fl(ijk),9) ) THEN
                  icnt = icnt + cell_neighbours(ijk, mpime, nset, rcv_cell_set, myijk)
                END IF
              END DO
            END DO
          END DO
        ELSE
          CALL error( ' ghost ', ' wrong job_type '//job_type, 1)
        END IF
      ELSE IF ( proc_map(mpime)%type == BLOCK3D_MAP ) THEN
        IF( job_type == '3D' ) THEN
          i1 = proc_map(mpime)%blkbsw(1)
          i2 = proc_map(mpime)%blktne(1)
          j1 = proc_map(mpime)%blkbsw(2)
          j2 = proc_map(mpime)%blktne(2)
          k1 = proc_map(mpime)%blkbsw(3)
          k2 = proc_map(mpime)%blktne(3)
          DO k = k1, k2
            DO j = j1, j2
              DO i = i1, i2
                ijk = i + (j-1)*nx + (k-1)*nx*ny
                IF ( BTEST(fl(ijk),0) .OR. BTEST(fl(ijk),9) ) THEN
                  icnt = icnt + cell_neighbours(ijk, mpime, nset, rcv_cell_set, myijk)
                END IF
              END DO
            END DO
          END DO
        ELSE
          CALL error( ' ghost ', ' wrong job_type '//job_type, 1)
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
 300      FORMAT(' # neighbours set SIZE ',i5,' from ',i3)
 310      FORMAT(10i8)
        END DO
      END IF

! ... the number cells required to update the physical quantities
! ... 'ncdom' is the sum of the local and neighbour cells
!
      ncext = SUM( nrcv ) 
      ncdom = ncint + ncext
      IF (lpr > 1) THEN
        WRITE(testunit,* ) ' # ncext ', ncext
        WRITE(testunit,* ) ' # ncdom ', ncdom
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

! ... print out basic information on the map
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
          ALLOCATE( snd_map(ipe)%isnd( snd_map(ipe)%nsnd ) )
          ALLOCATE( snd_map(ipe)%iloc( snd_map(ipe)%nsnd ) )
          snd_map(ipe)%isnd( 1:snd_map(ipe)%nsnd ) = ijksnd( 1:snd_map(ipe)%nsnd )
          DO ijk = 1, snd_map(ipe)%nsnd
            snd_map(ipe)%iloc( ijk ) = cell_g2l( ijksnd( ijk ), mpime )
          END DO
        END IF 
!
          DEALLOCATE(ijkrcv)
          DEALLOCATE(ijksnd)

      END DO
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
      ALLOCATE( flag(ncdom) )
      DO ijkl = 1, ncint
        ijk = myijk( ip0_jp0_kp0_, ijkl)
        flag(ijkl) = fl(ijk)
      END DO
      DEALLOCATE (fl)

      CALL data_exchange(flag)
!
! ... fill in the array myinds using myijk
!
      CALL set_myinds(myinds, myijk)
!
! ... Map the forcing points on local domains
! ... by using array numx/y/z. Scatter the array
! ... of forcing points among processors.
!
      IF (immb == 1) CALL local_forcing

      WRITE(testunit,*) 'End of Ghost'
!
      RETURN
      END SUBROUTINE ghost
!----------------------------------------------------------------------
!
      INTEGER FUNCTION cell_owner( ijk )

! ... gives the owner of a cells from its global index (ijk)

        USE dimensions
        USE parallel, ONLY: nproc
        USE control_flags, ONLY: job_type

        INTEGER, INTENT(IN) :: ijk
        INTEGER :: ipe, layer, cell_layer
        INTEGER :: i, j, k

        IF( job_type == '2D' ) THEN

          k = ( ijk - 1 ) / nx + 1
          i = MOD( ( ijk - 1 ), nx) + 1

        ELSE IF( job_type == '3D' ) THEN

          k = ( ijk - 1 ) / ( nx*ny ) + 1
          j = MOD( ijk - 1, nx*ny ) / nx + 1
          i = MOD( MOD( ijk - 1, nx*ny ), nx ) + 1

        ELSE

          CALL error(' cell_owner ', ' unknown job_type ', 1 )

        END IF
          
        DO ipe = 0, nproc - 1
          IF ( proc_map(ipe)%type == LAYER_MAP ) THEN 
            IF( ijk >= proc_map(ipe)%lay(1) .AND.  &
                ijk <= proc_map(ipe)%lay(2) )       cell_owner = ipe
          ELSE IF ( proc_map(ipe)%type == BLOCK2D_MAP ) THEN
            IF( job_type == '2D' ) THEN
              IF (k >= proc_map(ipe)%corner1(2) .AND.  &
                  k <= proc_map(ipe)%corner2(2) .AND.  &
                  i >= proc_map(ipe)%corner1(1) .AND.  &
                  i <= proc_map(ipe)%corner2(1) )   cell_owner = ipe
            ELSE IF ( job_type == '3D' ) THEN
              IF (j >= proc_map(ipe)%corner1(2) .AND.  &
                  j <= proc_map(ipe)%corner2(2) .AND.  &
                  i >= proc_map(ipe)%corner1(1) .AND.  &
                  i <= proc_map(ipe)%corner2(1) )   cell_owner = ipe
            ELSE
              CALL error(' cell_owner',' unknown job_type '//job_type, 1 )
            END IF
          ELSE IF ( proc_map(ipe)%type == BLOCK3D_MAP ) THEN
            IF( job_type == '3D' ) THEN
              IF( i >= proc_map(ipe)%blkbsw(1) .AND. &
                  i <= proc_map(ipe)%blktne(1) .AND. &
                  j >= proc_map(ipe)%blkbsw(2) .AND. &
                  j <= proc_map(ipe)%blktne(2) .AND. &
                  k >= proc_map(ipe)%blkbsw(3) .AND. &
                  k <= proc_map(ipe)%blktne(3) )  cell_owner = ipe
            ELSE
              CALL error(' cell_owner',' unknown job_type '//job_type, 1 )
            END IF
          ELSE
            CALL error(' cell_owner',' partition type not yet implemented ', &
                       proc_map(ipe)%type )
          END IF
        END DO

      RETURN
      END FUNCTION     

!----------------------------------------------------------------------
!
      INTEGER FUNCTION cell_l2g( ijkl, mpime)

        USE dimensions
        USE control_flags, ONLY: job_type
        USE parallel, ONLY: nproc

        INTEGER, INTENT(IN) :: ijkl, mpime
        INTEGER :: i,j,k, i1,i2,j1,j2,k1,k2, nxl, nyl, nzl 
!
        IF ( proc_map(mpime)%type == LAYER_MAP ) THEN

          cell_l2g = ijkl + proc_map( mpime )%lay(1) - 1

        ELSE IF ( proc_map(mpime)%type == BLOCK2D_MAP ) THEN

          IF( job_type == '2D' ) THEN

            k1 = proc_map(mpime)%corner1(2)
            i1 = proc_map(mpime)%corner1(1)
            i2 = proc_map(mpime)%corner2(1)
            i = MOD( ( ijkl - 1 ), (i2-i1+1)) + i1
            k = ( ijkl - 1 ) / (i2-i1+1) + k1
            cell_l2g = i + (k-1) * nx

          ELSE IF( job_type == '3D' ) THEN

            j1 = proc_map(mpime)%corner1(2)
            i1 = proc_map(mpime)%corner1(1)
            j2 = proc_map(mpime)%corner2(2)
            i2 = proc_map(mpime)%corner2(1)
            nxl = i2 - i1 + 1
            nyl = j2 - j1 + 1
            nzl = nz
            
            i = MOD( MOD( ijkl - 1, nxl*nyl ), nxl ) + i1
            j = MOD( ijkl - 1, nxl*nyl ) / nxl + j1
            k = ( ijkl - 1 ) / ( nxl*nyl ) + 1

            cell_l2g = i + (j-1)*nx + (k-1)*nx*ny

          ELSE

            CALL error(' cell_owner ', ' unknown job_type ', 1 )

          END IF

        ELSE IF ( proc_map(mpime)%type == BLOCK3D_MAP ) THEN

          IF( job_type == '3D' ) THEN

            k1 = proc_map(mpime)%blkbsw(3)
            j1 = proc_map(mpime)%blkbsw(2)
            i1 = proc_map(mpime)%blkbsw(1)
            k2 = proc_map(mpime)%blktne(3)
            j2 = proc_map(mpime)%blktne(2)
            i2 = proc_map(mpime)%blktne(1)
            nxl = i2 - i1 + 1
            nyl = j2 - j1 + 1
            nzl = k2 - k1 + 1
            
            i = MOD( MOD( ijkl - 1, nxl*nyl ), nxl ) + i1
            j = MOD( ijkl - 1, nxl*nyl ) / nxl + j1
            k = ( ijkl - 1 ) / ( nxl*nyl ) + k1

            cell_l2g = i + (j-1)*nx + (k-1)*nx*ny


          ELSE

            CALL error(' cell_owner ', ' unknown job_type ', 1 )

          END IF

        ELSE

          CALL error(' cell_l2g ',' partition type not yet implemented ', &
                     proc_map(mpime)%type )

        END IF
!
        RETURN
      END FUNCTION     
!
!----------------------------------------------------------------------
!
      INTEGER FUNCTION cell_g2l(ijk, mpime)

        USE dimensions
        USE parallel, ONLY: nproc
        USE control_flags, ONLY: job_type

        INTEGER, INTENT(IN) :: ijk, mpime
        INTEGER :: i, j, k, i1, i2, j1, j2, k1, k2, nxl, nyl, nzl
!
        IF ( proc_map(mpime)%type == LAYER_MAP ) THEN

          cell_g2l = ijk - proc_map(mpime)%lay(1) + 1

        ELSE IF ( proc_map(mpime)%type == BLOCK2D_MAP ) THEN

          IF( job_type == '2D' ) THEN

            k = ( ijk - 1 ) / nx + 1
            i = MOD( ( ijk - 1 ), nx) + 1

            k1 = proc_map(mpime)%corner1(2)
            i1 = proc_map(mpime)%corner1(1)
            i2 = proc_map(mpime)%corner2(1)
            cell_g2l = (i-i1+1) + (k-k1)*(i2-i1+1)

          ELSE IF( job_type == '3D' ) THEN

            i = MOD( MOD( ijk - 1, nx*ny ), nx ) + 1
            j = MOD( ijk - 1, nx*ny ) / nx + 1
            k = ( ijk - 1 ) / ( nx*ny ) + 1

            j1 = proc_map(mpime)%corner1(2)
            i1 = proc_map(mpime)%corner1(1)
            j2 = proc_map(mpime)%corner2(2)
            i2 = proc_map(mpime)%corner2(1)

            nxl = i2 - i1 + 1
            nyl = j2 - j1 + 1
            nzl = nz

            cell_g2l = 1 + (i-i1) + (j-j1)*nxl + (k-1)*nxl*nyl

          ELSE

            CALL error(' cell_g2l ', ' unknown job_type ', 1 )

          END IF

        ELSE IF ( proc_map(mpime)%type == BLOCK3D_MAP ) THEN

          IF( job_type == '3D' ) THEN

            i = MOD( MOD( ijk - 1, nx*ny ), nx ) + 1
            j = MOD( ijk - 1, nx*ny ) / nx + 1
            k = ( ijk - 1 ) / ( nx*ny ) + 1

            k1 = proc_map(mpime)%blkbsw(3)
            j1 = proc_map(mpime)%blkbsw(2)
            i1 = proc_map(mpime)%blkbsw(1)
            k2 = proc_map(mpime)%blktne(3)
            j2 = proc_map(mpime)%blktne(2)
            i2 = proc_map(mpime)%blktne(1)

            nxl = i2 - i1 + 1
            nyl = j2 - j1 + 1
            nzl = k2 - k1 + 1

            cell_g2l = 1 + (i-i1) + (j-j1)*nxl + (k-k1)*nxl*nyl

          ELSE

            CALL error(' cell_g2l ', ' unknown job_type ', 1 )

          END IF

        ELSE

          CALL error(' cell_g2l ', ' partition type not yet implemented ', proc_map(mpime)%type )

        END IF
!
        RETURN
      END FUNCTION     
!
!----------------------------------------------------------------------
      INTEGER FUNCTION cell_neighbours(ijk, mpime, nset, rcv_cell_set, myijk)
!
! ... loop over neighbour cells checking the owner. Fill in the 
! ... myijk array that identifies the ghost cells
!
        USE dimensions
        USE basic_types, ONLY: imatrix
        USE control_flags, ONLY: job_type
!
        INTEGER, INTENT(IN) :: ijk
        INTEGER, INTENT(IN) :: mpime
        INTEGER, INTENT(INOUT) :: nset(0:)
        TYPE (imatrix), OPTIONAL :: rcv_cell_set(0:)
        INTEGER, OPTIONAL :: myijk(:,:)
!
        INTEGER ::  ii, jj, kk
        INTEGER :: icnt, ipe, i, j, k, ijke, ijkel, ijkl
        INTEGER :: im, jm, km
        LOGICAL :: fill

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

        IF( job_type == '2D' ) THEN

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
                  ELSE IF( fill ) THEN
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

        ELSE IF( job_type == '3D' ) THEN

          i = MOD( MOD( ijk - 1, nx*ny ), nx ) + 1
          j = MOD( ijk - 1, nx*ny ) / nx + 1
          k = ( ijk - 1 ) / ( nx*ny ) + 1

! ...     loop over the neighbouring cells
          DO km = -2, 2
            DO jm = -2, 2
              DO im = -2, 2
                IF( ( ABS( im ) + ABS( jm ) + ABS( km ) ) <= 2 ) THEN
                  IF( (im /= 0) .OR. (jm /= 0) .OR. (km /= 0) ) THEN
                    ii = im
                    jj = jm
                    kk = km
                    IF( ( i == 2    ) .AND. ( ii == -2 ) ) ii = -1
                    IF( ( i == nx-1 ) .AND. ( ii == +2 ) ) ii = +1
                    IF( ( j == 2    ) .AND. ( jj == -2 ) ) jj = -1
                    IF( ( j == ny-1 ) .AND. ( jj == +2 ) ) jj = +1
                    IF( ( k == 2    ) .AND. ( kk == -2 ) ) kk = -1
                    IF( ( k == nz-1 ) .AND. ( kk == +2 ) ) kk = +1
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
                  ELSE IF( fill ) THEN
! ...               store the global cell index ij, of the the local cell ijkl
                    myijk( ip0_jp0_kp0_, ijkl) = ijk
                  END IF
                END IF
              END DO
            END DO
          END DO

        ELSE 

          CALL error(' cell_neighbours ', ' wrong job_type ', 1 )

        END IF

        cell_neighbours = icnt

        RETURN
      END FUNCTION
!
!----------------------------------------------------------------------
!
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

!---------------------------------------------------
      SUBROUTINE set_myinds(myinds, myijk)
!
        USE dimensions
        USE control_flags, ONLY: job_type
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

          IF( job_type == '2D' ) THEN

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
                CASE (noslip_wall, slip_wall, filled_cell)
                        ijke = ijk
                CASE DEFAULT
                        ijke = ipjk
              END SELECT

              SELECT CASE (flag(imjk))
                CASE (noslip_wall, slip_wall, filled_cell)
                        ijkw = ijk
                CASE DEFAULT
                        ijkw = imjk
              END SELECT

              SELECT CASE (flag(ijkp))
                CASE (noslip_wall, slip_wall, filled_cell)
                        ijkt = ijk
                CASE DEFAULT
                        ijkt = ijkp
              END SELECT

              SELECT CASE (flag(ijkm))
                CASE (noslip_wall, slip_wall, filled_cell)
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
!
! ... Second neighbours are not available on boundaries
!
              ijkee = ippjk
              IF( (flag(ippjk) == slip_wall) .OR. (flag(ippjk) == noslip_wall) ) ijkee = ipjk
              IF(i == (nx-1)) ijkee = ijke

              ijktt = ijkpp
              IF( (flag(ijkpp) == slip_wall) .OR. (flag(ijkpp) == noslip_wall) ) ijktt = ijkp
              IF(k == (nz-1)) ijktt = ijkt
  
              ijkww = immjk
              IF( (flag(immjk) == slip_wall) .OR. (flag(immjk) == noslip_wall) ) ijkww = imjk
              IF(i == 2) ijkww = ijkw

              ijkbb = ijkmm
              IF( (flag(ijkmm) == slip_wall) .OR. (flag(ijkmm) == noslip_wall) ) ijkbb = ijkm
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

          ELSE IF( job_type == '3D' ) THEN

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
  
              SELECT CASE (flag(ipjk))
                CASE (noslip_wall, slip_wall, filled_cell)
                        ijke = ijk
                CASE DEFAULT
                        ijke  =  ipjk
              END SELECT

              SELECT CASE (flag(imjk))
                CASE (noslip_wall, slip_wall, filled_cell)
                        ijkw = ijk
                CASE DEFAULT
                        ijkw  =  imjk
              END SELECT

              SELECT CASE (flag(ijpk))
                CASE (noslip_wall, slip_wall, filled_cell)
                        ijkn = ijk
                CASE DEFAULT
                        ijkn  =  ijpk
              END SELECT

              SELECT CASE (flag(ijmk))
                CASE (noslip_wall, slip_wall, filled_cell)
                        ijks = ijk
                CASE DEFAULT
                        ijks  =  ijmk
              END SELECT

              SELECT CASE (flag(ijkp))
                CASE (noslip_wall, slip_wall, filled_cell)
                        ijkt = ijk
                CASE DEFAULT
                        ijkt  =  ijkp
              END SELECT

              SELECT CASE (flag(ijkm))
                CASE (noslip_wall, slip_wall, filled_cell)
                        ijkb = ijk
                CASE DEFAULT
                        ijkb  =  ijkm
              END SELECT

!
! ... Second neighbours are not available on boundaries
!
              ijkee =  ippjk
              if( flag( ippjk ) == slip_wall .OR. flag( ippjk ) == noslip_wall ) ijkee = ijke
              if( i == (nx-1) ) ijkee = ijke

              ijkww =  immjk
              if( flag( immjk ) == slip_wall .OR. flag( immjk ) == noslip_wall ) ijkww = ijkw
              if( (i == 2) ) ijkww = ijkw

              ijknn =  ijppk
              if( flag( ijppk ) == slip_wall .OR. flag( ijppk ) == noslip_wall ) ijknn = ijkn
              if( (j == (ny-1)) ) ijknn = ijkn

              ijkss =  ijmmk
              if( flag( ijmmk ) == slip_wall .OR. flag( ijmmk ) == noslip_wall ) ijkss = ijks
              if( (j == 2) ) ijkss = ijks

              ijktt =  ijkpp
              if( flag( ijkpp ) == slip_wall .OR. flag( ijkpp ) == noslip_wall ) ijktt = ijkt
              if( k == (nz-1) ) ijktt = ijkt

              ijkbb =  ijkmm
              if( flag( ijkmm ) == slip_wall .OR. flag( ijkmm ) == noslip_wall ) ijkbb = ijkb
              if( k == 2 ) ijkbb = ijkb

              !!!!! check diagonals

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
            WRITE(testunit,*) 'Trying to allocate ', MAX(rcv_map(isour)%nrcv,1), ' elements '
            CALL error(' data_exchange_r ', ' allocating rcvbuf ', ierr )
          END IF
          ALLOCATE( sndbuf( MAX(snd_map(idest)%nsnd,1) ), STAT=ierr )
          IF( ierr /= 0 ) THEN
            WRITE(testunit,*) 'Trying to allocate ', MAX(snd_map(idest)%nsnd,1), ' elements '
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
        RETURN
      END SUBROUTINE data_exchange_r

!----------------------------------------------------------------------

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
            WRITE(testunit,*) 'Trying to allocate ', rdim, ' elements '
            CALL error(' data_exchange_rm ', ' allocating rcvbuf ', ierr )
          END IF
          ALLOCATE( sndbuf( sdim ), STAT=ierr )
          IF( ierr /= 0 ) THEN
            WRITE(testunit,*) 'Trying to allocate ', sdim, ' elements '
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
            WRITE(testunit,*) 'Trying to allocate ', MAX(rcv_map(isour)%nrcv,1), ' elements '
            CALL error(' data_exchange_i ', ' allocating rcvbuf ', ierr )
          END IF
          ALLOCATE( sndbuf( MAX(snd_map(idest)%nsnd,1) ), STAT=ierr )
          IF( ierr /= 0 ) THEN
            WRITE(testunit,*) 'Trying to allocate ', MAX(snd_map(idest)%nsnd,1), ' elements '
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
            WRITE(testunit,*) 'Trying to allocate ', MAX(rcv_map(isour)%nrcv,1), ' elements '
            CALL error(' data_exchange_l ', ' allocating rcvbuf ', ierr )
          END IF
          ALLOCATE( sndbuf( MAX(snd_map(idest)%nsnd,1) ), STAT=ierr )
          IF( ierr /= 0 ) THEN
            WRITE(testunit,*) 'Trying to allocate ', MAX(snd_map(idest)%nsnd,1), ' elements '
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
        USE kinds, ONLY: sgl

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
        USE kinds, ONLY: sgl

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
      SUBROUTINE meshinds(localindex,globalindex,i,j,k)
!
! ... computes the coordinates of a grid point in 
! ... a 2D or 3D rectilinear mesh

        USE dimensions
        USE control_flags, ONLY: job_type
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: localindex
        INTEGER, INTENT(OUT) :: globalindex, i, j, k

        i = 0    ! index for   x(r)   coordinate
        j = 0    ! index for   y      coordinate
        k = 0    ! index for   z      coordinate

        globalindex = myijk( ip0_jp0_kp0_, localindex)

        IF (job_type == '3D') THEN

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

        ELSE IF (job_type == '2D') THEN

          k = ( globalindex - 1 ) / nx + 1
          i = MOD( ( globalindex - 1 ), nx) + 1

        ELSE 

          CALL error(' meshinds ', ' unknow job_type '//job_type, 1 )

        END IF

      END SUBROUTINE meshinds
!----------------------------------------------------------------------
      SUBROUTINE local_forcing
      USE control_flags, ONLY: job_type
      USE dimensions, ONLY: nx, ny, nz
      USE immersed_boundaries, ONLY: numx, fptx, numy, fpty, numz, fptz
      USE parallel, ONLY: mpime

      IMPLICIT NONE
      INTEGER :: n, i, j, k, ijk, ijkl
      INTEGER :: nfpx, nfpy, nfpz
!
        ALLOCATE( numx(ncdom) ); numx = 0
        nfpx = SIZE(fptx)

        IF (job_type == '3D') THEN
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
set_numx: IF (i/=0 .AND. k/=0) THEN
            IF (job_type == '2D') THEN
              ijk = i + (k-1) * nx
            ELSE IF (job_type == '3D') THEN
              IF (j == 0) CALL error('decomp','control numx',1)
              ijk = i + (j-1) * nx + (k-1) * nx * ny
            END IF
            IF( cell_owner(ijk) == mpime ) THEN

              ijkl = cell_g2l(ijk,mpime)
              numx(ijkl) = n 

              IF (k==1 .OR. k==nz .OR. i==1 .OR. i==nx .OR. j==1 .OR. j==ny) THEN
                WRITE(errorunit,*) 'WARNING! from ghost'
                WRITE(errorunit,*) 'skipping x-forcing on boundaries', i, j, k
              END IF
            END IF
          ELSE
            CALL error('decomp','control numx',1)
          END IF set_numx

        END DO
        !
        DO n = 1, nfpy
        
          IF (job_type == '3D') THEN
            i = fpty(n)%i
            j = fpty(n)%j
            k = fpty(n)%k
set_numy:   IF (i/=0 .AND. k/=0) THEN
              IF (job_type == '2D') THEN
                ijk = i + (k-1) * nx
              ELSE IF (job_type == '3D') THEN
                IF (j == 0) CALL error('decomp','control numy',1)
                ijk = i + (j-1) * nx + (k-1) * nx * ny
              END IF
              IF( cell_owner(ijk)== mpime ) THEN

                ijkl = cell_g2l(ijk,mpime)
                numy(ijkl) = n 

                IF (k==1 .OR. k==nz .OR. i==1 .OR. i==nx .OR. j==1 .OR. j==ny) THEN
                  WRITE(errorunit,*) 'WARNING! from ghost'
                  WRITE(errorunit,*) 'skipping y-forcing on boundaries', i, j, k
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
set_numz: IF (i/=0 .AND. k/=0) THEN
            IF (job_type == '2D') THEN
              ijk = i + (k-1) * nx
            ELSE IF (job_type == '3D') THEN
              IF (j == 0) CALL error('decomp','control numz',n)
              ijk = i + (j-1) * nx + (k-1) * nx * ny
            END IF
            IF( cell_owner(ijk)== mpime ) THEN

              ijkl = cell_g2l(ijk,mpime)
              numz(ijkl) = n 

              IF (k==1 .OR. k==nz .OR. i==1 .OR. i==nx .OR. j==1 .OR. j==ny) THEN
                WRITE(errorunit,*) 'WARNING! from ghost'
                WRITE(errorunit,*) 'skipping z-forcing on boundaries', i, j, k
              END IF
            END IF
          ELSE
            CALL error('decomp','control numz',n)
          END IF set_numz
        
        END DO

        CALL data_exchange(numx)
        IF (job_type == '3D') CALL data_exchange(numy)
        CALL data_exchange(numz)

        CALL fill_cells_int

      RETURN
      END SUBROUTINE local_forcing
!----------------------------------------------------------------------
      SUBROUTINE fill_cells_int
!
! ... Coefficients b(1:6) are multiplied to the numerical mass
! ... fluxes to take into account that some cell face could be
! ... immersed. b=1 if the cell face is external, b=0 if the
! ... cell face is inside the topography.
!
      USE control_flags, ONLY: job_type
      USE grid, ONLY: z, zb, flag, dz, dome_cell
      USE immersed_boundaries, ONLY: bd, vf
      USE io_files, ONLY: tempunit
      USE volcano_topography, ONLY: topo_c, topo_x
      USE volcano_topography, ONLY: topo2d_c, topo2d_x, topo2d_y
      USE indijk_module
      IMPLICIT NONE

      INTEGER :: i,j,k,ijk,imesh
      INTEGER :: ipjk, imjk, ijpk, ijmk, ijkm
      INTEGER :: filled
!
! ... Allocate and initialize coefficients
!
      ALLOCATE(bd(ncint))
      ALLOCATE(vf(ncint))

      IF (job_type == '2D') THEN
        filled = 15
      ELSE IF (job_type == '3D') THEN
        filled = 63
      END IF

      bd(:) = 0
      vf(:) = 0.D0

      IF( lpr > 1 ) THEN
        WRITE( 7, * ) 
        WRITE( 7, * ) 'b coefficients for immersed bndaries' 
        WRITE( 7, * ) '     ijk   i   j   k     b (int)'
      END IF
!
        DO ijk=1, ncint

          imjk  = myijk( im1_jp0_kp0_, ijk )
          ipjk  = myijk( ip1_jp0_kp0_, ijk )
          ijmk  = myijk( ip0_jm1_kp0_, ijk )
          ijpk  = myijk( ip0_jp1_kp0_, ijk )
          ijkm  = myijk( ip0_jp0_km1_, ijk )

          IF( BTEST(flag(ijk),0) ) THEN
            CALL meshinds(ijk,imesh,i,j,k)
  
            IF (job_type == '2D') THEN
              ! 
              ! East
!              IF (z(k) > topo_x(i) .AND. flag(ipjk) /= filled_cell) THEN
               IF (z(k) > topo_x(i)) THEN
                  IF (flag(ipjk) == filled_cell) THEN
                     vf(ijk) = vf(ijk) + 2.D0
                  ELSE
                     bd(ijk) = bd(ijk) + 1
                     vf(ijk) = vf(ijk) + 1.D0
                  END IF
               END IF
              !
              ! West
!              IF (z(k) > topo_x(i-1) .AND. flag(imjk) /= filled_cell) THEN
               IF (z(k) > topo_x(i-1)) THEN
                  IF (flag(imjk) == filled_cell) THEN
                     vf(ijk) = vf(ijk) + 2.D0
                  ELSE
                     bd(ijk) = bd(ijk) + 2
                     vf(ijk) = vf(ijk) + 1.D0
                  END IF
               END IF
              !  
              ! Top
               IF (zb(k) > topo_c(i)) THEN
                  bd(ijk) = bd(ijk) + 4 
                  vf(ijk) = vf(ijk) + 1.D0
               END IF
              !
              ! Bottom
!              IF (zb(k-1) >= topo_c(i) .AND. flag(ijkm) /= filled_cell) THEN
               IF (zb(k-1) >= topo_c(i)) THEN
                  IF (flag(ijkm) == filled_cell) THEN
                     vf(ijk) = vf(ijk) + 2.D0   
                  ELSE
                     bd(ijk) = bd(ijk) + 8  
                     vf(ijk) = vf(ijk) + 1.D0
                  END IF
               END IF
              !
            ELSE IF (job_type == '3D') THEN
              ! 
              ! East
!              IF (z(k) > topo2d_x(i,j) .AND. flag(ipjk) /= filled_cell) THEN
              IF (z(k) > topo2d_x(i,j)) THEN
                 IF (flag(ipjk) == filled_cell) THEN
                    vf(ijk) = vf(ijk) + 2.D0
                 ELSE
                    bd(ijk) = bd(ijk) + 1
                    vf(ijk) = vf(ijk) + 1.D0
                 END IF
              END IF
              !
              ! West
!              IF (z(k) > topo2d_x(i-1,j) .AND. flag(imjk) /= filled_cell) THEN
              IF (z(k) > topo2d_x(i-1,j)) THEN
                 IF (flag(imjk) == filled_cell) THEN
                    vf(ijk) = vf(ijk) + 2.D0  
                 ELSE
                    bd(ijk) = bd(ijk) + 2
                    vf(ijk) = vf(ijk) + 1.D0
                 END IF
              END IF
              !
              ! Top
              IF (zb(k) > topo2d_c(i,j))  THEN
                bd(ijk) = bd(ijk) + 4 
                vf(ijk) = vf(ijk) + 1.D0
              END IF
              !
              ! Bottom
!              IF (zb(k-1) >= topo2d_c(i,j) .AND. flag(ijkm) /= filled_cell) THEN
              IF (zb(k-1) >= topo2d_c(i,j)) THEN
                IF (flag(ijkm) == filled_cell) THEN
                   vf(ijk) = vf(ijk) + 2.D0                   
                ELSE
                   bd(ijk) = bd(ijk) + 8
                   vf(ijk) = vf(ijk) + 1.D0
                END IF
             END IF
              !
              ! North
!              IF (z(k) > topo2d_y(i,j) .AND. flag(ijpk) /= filled_cell) THEN
              IF (z(k) > topo2d_y(i,j)) THEN
                IF (flag(ijpk) == filled_cell) THEN
                   vf(ijk) = vf(ijk) + 2.D0
                ELSE
                   bd(ijk) = bd(ijk) + 16
                   vf(ijk) = vf(ijk) + 1.D0
                END IF
             END IF
              !
              ! South
!              IF (z(k) > topo2d_y(i,j-1) .AND. flag(ijmk) /= filled_cell) THEN
              IF (z(k) > topo2d_y(i,j-1)) THEN
                 IF (flag(ijmk) == filled_cell) THEN
                    vf(ijk) = vf(ijk) + 1.D0
                 ELSE
                    bd(ijk) = bd(ijk) + 32
                    vf(ijk) = vf(ijk) + 1.D0
                 END IF
              END IF
              !
            END IF

            IF (job_type == '2D') THEN
              vf(ijk) = 0.25D0 * vf(ijk)
            ELSE IF( job_type == '3D') THEN
              vf(ijk) = vf(ijk) / 6.D0
            END IF

            ! ... cells that are completely immersed
            ! ... are excluded from computation. 
            !
            IF (vf(ijk) == 0.D0) flag(ijk) = filled_cell

            IF (bd(ijk) /= filled  .AND. lpr > 1 ) THEN
              WRITE( 7, fmt = "( I8,3I4,2X,B8,I4 )" ) ijk, i, j, k, bd(ijk), flag(ijk)
            END IF

          END IF

        END DO
        OPEN(tempunit,FILE='vf.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')  
        WRITE(tempunit) vf
        CLOSE(tempunit)
        CALL data_exchange(flag)
      
      RETURN
      END SUBROUTINE fill_cells_int
!----------------------------------------------------------------------
      SUBROUTINE fill_cells_real
!
! ... Coefficients b(1:6) are multiplied to the numerical mass
! ... fluxes to take into account that some cell face could be
! ... immersed. b=1 if the cell face is external, b=0 if the
! ... cell face is inside the topography.
!
      USE control_flags, ONLY: job_type, lpr
      USE grid, ONLY: z, zb, flag, dz
      USE immersed_boundaries, ONLY: bdr, vf
      USE volcano_topography, ONLY: topo_c, topo_x
      USE volcano_topography, ONLY: topo2d_c, topo2d_x, topo2d_y
      IMPLICIT NONE

      INTEGER :: i,j,k,ijk,imesh
      INTEGER :: filled
      REAL*8 :: alpha
!
! ... Allocate and initialize coefficients
!
      ALLOCATE(bdr(ncint,6)); bdr(:,:) = 0.D0
      ALLOCATE(vf(ncint)); vf(:) = 0.D0
!
      IF( lpr > 1 ) THEN
        WRITE( 7, * ) 
        WRITE( 7, * ) 'b coefficients for immersed bndaries' 
        WRITE( 7, * ) '     ijk   i   j   k     b (real)'
      END IF
!
        DO ijk=1, ncint
          IF( BTEST(flag(ijk),0) ) THEN
            CALL meshinds(ijk,imesh,i,j,k)
  
            alpha = 0.D0
            IF (job_type == '2D') THEN
              !
              alpha = MAX(zb(k) - topo_x(i),0.D0) ! East
              alpha = MIN(alpha / dz(k),1.D0)
              vf(ijk) = vf(ijk) + alpha
              bdr(ijk,1) = alpha**2/( 2.D0 * alpha - 1.D0 )
              !
              alpha = MAX(zb(k) - topo_x(i-1),0.D0) ! West
              alpha = MIN(alpha / dz(k),1.D0)
              vf(ijk) = vf(ijk) + alpha
              bdr(ijk,2) = alpha**2/( 2.D0 * alpha - 1.D0 )
              !
              IF (zb(k) - topo_c(i) >= 0.D0 ) bdr(ijk,3) = 1.D0 ! Top
              IF (zb(k-1) - topo_c(i) >= 0.D0) bdr(ijk,4) = 1.D0 ! Bottom
              !
            ELSE IF (job_type == '3D') THEN
              !
              alpha = MAX(zb(k) - topo2d_x(i,j),0.D0) ! East
              alpha = MIN(alpha / dz(k),1.D0)
              vf(ijk) = vf(ijk) + alpha
              bdr(ijk,1) = alpha**2/( 2.D0 * alpha - 1.D0 )
              !
              alpha = MAX(zb(k) - topo2d_x(i-1,j),0.D0) ! West
              alpha = MIN(alpha / dz(k),1.D0)
              vf(ijk) = vf(ijk) + alpha
              bdr(ijk,2) = alpha**2/( 2.D0 * alpha - 1.D0 )
              !
              IF (zb(k) > topo2d_c(i,j))      bdr(ijk,3) = 1.D0 ! Top
              IF (zb(k-1) >= topo2d_c(i,j))   bdr(ijk,4) = 1.D0 ! Bottom
              !
              alpha = MAX(zb(k) - topo2d_y(i,j),0.D0) ! North
              alpha = MIN(alpha / dz(k),1.D0)
              vf(ijk) = vf(ijk) + alpha
              bdr(ijk,5) = alpha**2/( 2.D0 * alpha - 1.D0 )
              !
              alpha = MAX(zb(k) - topo2d_y(i,j-1),0.D0) ! South
              alpha = MIN(alpha / dz(k),1.D0)
              vf(ijk) = vf(ijk) + alpha
              bdr(ijk,6) = alpha**2/( 2.D0 * alpha - 1.D0 )
              !
            END IF

            IF (job_type == '2D') THEN
              vf(ijk)  = 0.5D0 * vf(ijk)
            ELSE IF( job_type == '3D') THEN
              vf(ijk)  = 0.25D0 * vf(ijk)
            END IF

            IF( lpr > 1 ) WRITE( 7, fmt = "( I8,3I4,2X,6F8.4 )" ) ijk, i, j, k, bdr(ijk,:)

          END IF

        END DO
      
      RETURN
      END SUBROUTINE fill_cells_real
!----------------------------------------------------------------------
      END MODULE domain_decomposition
!----------------------------------------------------------------------
