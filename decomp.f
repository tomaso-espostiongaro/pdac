!----------------------------------------------------------------------
      MODULE domain_decomposition
!
! ... This module contains all the procedures used for the domain
! ... decomposition
!
!----------------------------------------------------------------------
        USE grid, ONLY: fl, flag
        USE indijk_module
        USE immersed_boundaries, ONLY: immb
        USE control_flags, ONLY: lpr
        USE io_files, ONLY: errorunit, testunit, logunit
!
        IMPLICIT NONE
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
        TYPE (cells_map_type), ALLOCATABLE :: proc_map(:)
! ...     This array is used to store data about which cell belong
! ...     to which processor, and is used to find owner, local and
! ...     global index of a cell

! ...   The following arrays contain information on how the cells
!       are distributed among processors.
!
        INTEGER, ALLOCATABLE :: nctot(:)  ! balanced number of cell
        INTEGER, ALLOCATABLE :: ncfl1(:)  ! number of cell with fl=1
        INTEGER, ALLOCATABLE :: ncell(:)  ! unbalanced number of cell
        INTEGER, ALLOCATABLE :: ncdif(:)  ! ncfl1 - ncell
!
        !  ncell unbalanced number of cell,
        !        this is the number of cells on each process, as if they were 
        !        equally distributed among processors, regardless of their 
        !        computational weight.
        !  ncfl1 number of cells on each process with the flag equal 1
        !  nctot balanced number of cell, this is the number of cells on
        !        each process that it is computed to obtain the same 
        !        computational load on each process.
!
        INTEGER :: mesh_partition
        INTEGER :: countfl
!
        PUBLIC
        PRIVATE :: ncfl1, ncell, ncdif, countfl
!
        SAVE
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
        WRITE( logunit , * ) 
        WRITE( logunit , * ) 'Starting Grid Decomposition ...'
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
      countfl = 0
      DO ijk = 1, ntot
        IF (BTEST(fl(ijk),0)) countfl = countfl + 1
      END DO
!
      IF ( lpr > 0 .AND. ionode ) THEN
        WRITE( logunit, fmt = '(a13)' ) ' # countfl :'
        WRITE( logunit, fmt = '(10I8)') countfl
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
        WRITE( logunit, * ) 'End of Mesh decomposition'
!
      RETURN
      END SUBROUTINE partition
! ------------------------------------------------------------------------
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
      USE io_files, ONLY: testunit
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
        WRITE( logunit, * ) 'Using Layers Mesh decomposition' 
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
        WRITE( logunit, * ) 'Using Blocks ( 2D ) Mesh decomposition' 
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
        WRITE( logunit, * ) 'Using Colomns ( 3D ) Mesh decomposition' 

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
      INTEGER :: divnby, modnby

! ... Subroutine Body

      ionode = ( mpime == root )

      IF( ionode ) &
        WRITE( logunit, * ) 'Using Blocks ( 3D ) Grid partition' 
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
!
        IF (nbx == 0) nbx = 1
        IF (nby == 0) nby = 1
!
        DO WHILE ( nbx*nby < nprocxy )
         IF ( INT( nx / nbx ) > INT( ny / nby ) ) THEN
           nbx = nbx + 1
         ELSE
           nby = nby + 1
         END IF
        END DO
        rest = nbx*nby - nprocxy
        !
        ! ... Check for errors: 'rest' must be strictly less than 'nby'
        !
        divnby = INT(rest/nby)
        IF (divnby >= 1) THEN
          nbx = nbx - divnby
          rest = nbx*nby - nprocxy
        END IF
        !
        modnby = MOD(rest,nby)
        IF (modnby == 0) THEN
          nbx = nbx - rest/nby
          rest = 0
        END IF
!
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
!----------------------------------------------------------------------
      INTEGER FUNCTION cell_owner( ijk )

! ... gives the owner of a cells from its ->> global <<- index (ijk)

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
      END FUNCTION cell_owner
!----------------------------------------------------------------------
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
      END FUNCTION cell_l2g
!----------------------------------------------------------------------
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
      END FUNCTION cell_g2l
!----------------------------------------------------------------------
      END MODULE domain_decomposition
!----------------------------------------------------------------------
