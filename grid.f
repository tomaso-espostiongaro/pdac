!----------------------------------------------------------------------
      MODULE grid
!----------------------------------------------------------------------

        USE indijk_module

        IMPLICIT NONE
        SAVE
!
!
!
        TYPE rcv_map_type
          INTEGER :: nrcv              ! How Many Cells I Have to Receive
          INTEGER, POINTER :: ircv(:)  ! Which Cells I Have to Receive (Global Index) 
          INTEGER, POINTER :: iloc(:)  ! Where I Have to Put Them (Local Index)
        END TYPE 
!
!
        TYPE snd_map_type
          INTEGER :: nsnd              ! How Many Cells I Have to Send
          INTEGER, POINTER :: isnd(:)  ! Which Cells I Have to Send (Global Index)
          INTEGER, POINTER :: iloc(:)  ! Where I Have to Pick Them (Local Index)
        END TYPE 
!
!
        INTEGER, PARAMETER :: LAYER_MAP   = 1
        INTEGER, PARAMETER :: BLOCK2D_MAP = 2
        INTEGER, PARAMETER :: BLOCK3D_MAP = 3
!
        TYPE cells_map_type
          INTEGER :: type              ! Identify the map type (1 layer, 2 columns, 3 blocks)
          INTEGER :: lay(2)            ! lay(1) -> imesh start, lay(2) -> imesh end (layer)
          INTEGER :: colsw(2)          ! colw(1) -> SW x coord., colw(2) -> SW y coord. 
          INTEGER :: colne(2)          ! cole(1) -> NE x coord., colw(2) -> NE y coord.
          INTEGER :: blkbsw(3)         ! blkbsw(.) -> BSW x, y, z coordinates 
          INTEGER :: blktne(3)         ! blktne(.) -> TNE x, y, z coordinates
        END TYPE

!*******MX3D da rivedere
        INTEGER :: nbx, nby            ! number of blocks in x and y directions
!
        TYPE (rcv_map_type), ALLOCATABLE :: rcv_map(:)
        TYPE (snd_map_type), ALLOCATABLE :: snd_map(:)
        TYPE (cells_map_type), ALLOCATABLE :: proc_map(:)

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
        REAL*8, DIMENSION(:), ALLOCATABLE :: r, rb, dr
        REAL*8, DIMENSION(:), ALLOCATABLE :: dx, dy, dz
        REAL*8, DIMENSION(:), ALLOCATABLE :: xb, yb, zb
        REAL*8, DIMENSION(:), ALLOCATABLE :: indx, indy, indz
        REAL*8, DIMENSION(:), ALLOCATABLE :: inr, inrb, indr
!
        INTEGER :: itc, mesh_partition

!       ncint (ncint) number of local cells, 
!       ncext         number of external (ghost) cells, 
!       ncdom         sum of the previous two
!
        INTEGER :: ncint, ncext, ncdom  
!
        INTEGER, ALLOCATABLE :: myijk(:,:)
        INTEGER, ALLOCATABLE :: myinds(:,:)
!
!

        TYPE blbody
          INTEGER :: xlo
          INTEGER :: xhi
          INTEGER :: ylo
          INTEGER :: yhi
          INTEGER :: zlo
          INTEGER :: zhi
          INTEGER :: rlo
          INTEGER :: rhi
          INTEGER :: typ
        END TYPE
!
        TYPE (blbody), ALLOCATABLE :: iob(:)
!
        LOGICAL, PRIVATE :: outppm = .FALSE.

!*******MX3D fl da eliminare in un secondo tempo
!
        INTEGER, DIMENSION(:), ALLOCATABLE :: fl

        INTEGER, DIMENSION(:), ALLOCATABLE :: fl_l
!
        INTEGER :: countfl(5)

        INTERFACE data_exchange
          MODULE PROCEDURE data_exchange_i, data_exchange_r, data_exchange_rm
        END INTERFACE

        INTERFACE data_collect
          MODULE PROCEDURE data_collect_r
        END INTERFACE

        INTERFACE data_distribute
          MODULE PROCEDURE data_distribute_r
        END INTERFACE

!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------

      SUBROUTINE bounds_grid
!
        USE dimensions
        USE control_flags, ONLY: job_type
!
        IMPLICIT NONE
!
! ...   Set the appropriate total number of cells
        IF( job_type == '2D' ) THEN
          ntot = nr*nz
        ELSE IF( job_type == '3D' ) THEN
          ntot = nx*ny*nz
        ELSE
          CALL error( ' bounds_grid ', ' wrong job_type '//job_type, 1)
        END IF

        ALLOCATE( fl(ntot) )
        ALLOCATE( r(nr), rb(nr), dr(nr) )
        ALLOCATE( inr(nr), inrb(nr), indr(nr) ) 
        ALLOCATE( xb(nx), dx(nx), yb(ny), dy(ny), zb(nz), dz(nz) )
        ALLOCATE( indx(nx), indy(ny), indz(nz) )

        RETURN
      END SUBROUTINE
!
!----------------------------------------------------------------------
!
      SUBROUTINE bounds_blbody
        USE dimensions
        IMPLICIT NONE
          ALLOCATE(iob(no))
        RETURN
      END SUBROUTINE
!
!----------------------------------------------------------
!
        SUBROUTINE grid_setup(zzero)
!
          USE dimensions
          USE control_flags, ONLY: job_type

          REAL*8, INTENT(IN) :: zzero
          REAL*8, PARAMETER :: VERYBIG = 1.0d+10
          INTEGER :: i, j

          IF( job_type == '2D' ) THEN
!
! ...       ONLY FOR CYLINDRICAL COORDINATES
!
! ...       Set "r" dimension
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
            IF (rb(1) .EQ. 0.0D0) THEN
              inrb(1) = VERYBIG
            ELSE
              inrb(1) = 1.0D0/rb(1)
            END IF
            IF (itc.EQ.0) inrb(1)=1.D0
            DO i=2,nr
              inrb(i)=1.D0/rb(i)
            END DO  
!
            DO i=1,nr
              inr(i)=1.D0/r(i)
            END DO  
            DO i=1,nr
              indr(i)=1.D0/dr(i)
            END DO  
!
! ...       Set "z" dimension
            zb(1) = zzero
            DO j=1,(nz-1)
              zb(j+1)=zb(j)+dz(j+1)
            END DO
!
            DO j=1,(nz-1)
              indz(j)=1.D0/dz(j)
            END DO  

          ELSE IF( job_type == '3D' ) THEN
!
! ...       CARTESIAN COORDINATES
! 
!                first cell
!           |<----- dx(1) ---->|<----- dx(2) ----->|<----- .... ---|
!                            xb(1)                xb(2)
!
! ...       Set "x" dimension
            xb(1) = 0.0d0
            DO j=1,(nx-1)
              xb(j+1)=xb(j)+dx(j+1)
            END DO
!
            DO j=1,(nx-1)
              indx(j)=1.D0/dx(j)
            END DO 
!
! ...       Set "y" dimension
            yb(1) = 0.0d0
            DO j=1,(ny-1)
              yb(j+1)=yb(j)+dy(j+1)
            END DO
!
            DO j=1,(ny-1)
              indy(j)=1.D0/dy(j)
            END DO
!
! ...       Set "z" dimension
            zb(1) = zzero
            DO j=1,(nz-1)
              zb(j+1)=zb(j)+dz(j+1)
            END DO
!
            DO j=1,(nz-1)
              indz(j)=1.D0/dz(j)
            END DO  
  
          ELSE
    
            CALL error(' grid_setup ', ' unknow job_type '//job_type, 1 )

          END IF
!
          RETURN 
        END SUBROUTINE

!
!----------------------------------------------------------------------
!
      SUBROUTINE flic
!
        USE dimensions
        USE control_flags, ONLY: job_type
!
        IMPLICIT NONE
!
        INTEGER :: i, j, k, ij, n, ijk
!
        fl = 1
!
        IF( no <= 0 ) RETURN
!
! ...   no is the number of cell blocks (blunt body)
!
        IF( job_type == '2D' ) THEN

          DO n = 1, no
            DO j = iob(n)%zlo, iob(n)%zhi 
              DO i = iob(n)%rlo, iob(n)%rhi
                ij=i+(j-1)*nr
                SELECT CASE ( iob(n)%typ )
                  CASE (2) 
                    fl(ij) = 2
                  CASE (3) 
                    fl(ij) = 3
                  CASE (4) 
                    fl(ij) = 4
                  CASE (5) 
                    fl(ij) = 5
                  CASE (6) 
                    fl(ij) = 6
                  CASE DEFAULT
                    fl(ij) = 1
                END SELECT
              END DO
            END DO
          END DO

        ELSE IF( job_type == '3D' ) THEN

          DO n = 1, no
            DO k = iob(n)%zlo, iob(n)%zhi
              DO j = iob(n)%ylo, iob(n)%yhi
                DO i = iob(n)%xlo, iob(n)%xhi
                  ijk = i + ( (j-1) + (k-1)*ny ) * nx
                  SELECT CASE ( iob(n)%typ )
                    CASE (2)
                      fl(ijk) = 2
                    CASE (3)
                      fl(ijk) = 3
                    CASE (4)
                      fl(ijk) = 4
                    CASE (5)
                      fl(ijk) = 5
                    CASE (6)
                      fl(ijk) = 6
                    CASE DEFAULT
                      fl(ijk) = 1
                  END SELECT
                END DO
              END DO
            END DO
          END DO

        END IF

! ... MX3D CHECK THE CORNERS
!
        IF( job_type == '2D' ) THEN

          IF ( fl(nz*nr - 1) == 4 .AND. fl((nz-1)*nr) == 4) THEN
            fl(nz*nr) = 4
          END IF
          IF ( fl(nz*nr - 1) == 6 .AND. fl((nz-1)*nr) == 6) THEN
            fl(nz*nr) = 6
          END IF

        END IF
!
        RETURN
      END SUBROUTINE
!
!----------------------------------------------------------------------
      SUBROUTINE partition
!
! ... a routine that subdivides a regular grid between (nproc) 
! ... processors, balancing the computational load 
! ... Gives processor maps.
!
      USE dimensions
      USE parallel, ONLY: nproc, mpime, root
      IMPLICIT NONE
      SAVE
!
      INTEGER :: i, j, k, ijk, n, m, q, rr
      INTEGER :: nfl, ipe, dim
      INTEGER :: layer
      INTEGER :: localdim
!
! ... Subroutine Body
!
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
! ... here count the flags
!
      DO i = 1, SIZE(countfl)
        countfl(i) = COUNT( (fl == i) )
      END DO
!
      WRITE(7,'(a13,5i7)') ' # countfl ', countfl
!
! ... domain decomposition (build maps)
!     1=layers
!     2=blocks
!
      IF( nproc <= 2 ) mesh_partition = LAYER_MAP
!
      IF( ALLOCATED( proc_map ) ) DEALLOCATE( proc_map )
      ALLOCATE( proc_map(0:nproc-1) )

      proc_map(0:nproc-1)%type = mesh_partition

      IF ( proc_map(0)%type == LAYER_MAP ) THEN
        CALL layers( ncfl1, nctot )
      ELSE IF ( proc_map(0)%type == BLOCK2D_MAP ) THEN
        CALL blocks( ncfl1, nctot )
      ELSE
        CALL error(' partition ',' partition type not yet implemented ',proc_map(0)%type)
      END IF
!
! ... ncell is the balanced number of cells to each processor:
! 
      DO ipe = 0, nproc - 1
        ncell(ipe) = localdim(countfl(1), nproc, ipe)
        ncdif(ipe) = ncfl1(ipe) - ncell(ipe)
      END DO
!
      WRITE(7,*) '---  partition  -----'
      WRITE(7,*) '  '
      DO ipe = 0, nproc - 1
        WRITE(7,*) ' # nctot( ',ipe, ' )', nctot(ipe)
        WRITE(7,*) ' # ncfl1( ',ipe, ' )', ncfl1(ipe),' -', ncell(ipe),' =', ncdif(ipe)
        WRITE(7,*) 'proc(', ipe,')_map:', proc_map(ipe)%lay(1), proc_map(ipe)%lay(2)
      END DO
!
      RETURN
      END SUBROUTINE 
!
! ------------------------------------------------------------------------
!
      SUBROUTINE layers( ncfl1, nctot )
! 
! ... partition N.1 (layers)
! ... decomposes the domain into N horizontal layers
!
      USE dimensions
      USE control_flags, ONLY: job_type
      USE parallel, ONLY: nproc
!
      IMPLICIT NONE
      SAVE
!
      INTEGER :: ncfl1(0:)
      INTEGER :: nctot(0:)
!
      INTEGER :: i, j, k, ijk
      INTEGER :: icnt, icnt_ipe, ipe
      INTEGER :: localdim
      INTEGER, ALLOCATABLE :: lay_map(:,:)
!
      ALLOCATE(lay_map(0:nproc-1,2))
!
      ncfl1 = 0
      nctot = 0
!
      DO ipe = 0, nproc - 1
       ncfl1(ipe) = localdim(countfl(1), nproc, ipe)
      END DO
!
      IF( job_type == '2D' ) THEN

        lay_map(0,1) = 1
        lay_map(nproc-1,2) = nr*nz
        ipe = 0

        DO j = 1, nz
          DO i = 1, nr

            ijk = i + (j-1)*nr
  
            IF ( fl(ijk) == 1 ) THEN
              icnt_ipe = icnt_ipe + 1
              icnt     = icnt     + 1
            END IF

            nctot(ipe) = nctot(ipe) + 1

            IF ( ( icnt_ipe == ncfl1(ipe) ) .AND. ( i /= (nr-1) ) ) THEN
              icnt_ipe = 0
              IF( icnt < countfl(1) ) THEN
                lay_map(ipe,2) = ijk
                ipe = ipe + 1
                lay_map(ipe,1) = ijk+1
              END IF
            END IF

          END DO 
        END DO

      ELSE IF( job_type == '3D' ) THEN

        lay_map(0,1)       = 1
        lay_map(nproc-1,2) = nx*ny*nz
        ipe = 0

        DO k = 1, nz
          DO j = 1, ny
            DO i = 1, nx

              ijk = i + (j-1)*nx + (k-1)*nx*ny 
  
              IF ( fl(ijk) == 1 ) THEN
                icnt_ipe = icnt_ipe + 1
                icnt     = icnt     + 1
              END IF

              nctot(ipe) = nctot(ipe) + 1

              IF ( icnt_ipe >= ncfl1(ipe) ) THEN
                IF ( (i == nx-1) .OR.  (j == ny-1) .OR. & 
                     ( (j == ny) .AND. (i /= nx) ) ) GOTO 114
                  icnt_ipe = 0
                  IF( icnt < countfl(1) ) THEN
                    lay_map(ipe,2) = ijk
                    ipe = ipe + 1
                    lay_map(ipe,1) = ijk+1
                  END IF
 114            CONTINUE
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
      END SUBROUTINE

! ------------------------------------------------------------------------

!*****MX3D da rivedere
      SUBROUTINE blocks( ncfl1, nctot )
!
! ... partition N.2 (blocks)
!
      USE dimensions
      USE parallel, ONLY: nproc
!
      IMPLICIT NONE
      SAVE
!
      INTEGER :: ncfl1(0:)
      INTEGER :: nctot(0:)
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
      INTEGER, ALLOCATABLE :: lay_map(:,:)
      INTEGER, ALLOCATABLE :: nctot_lay(:), nbl_lay(:), ncfl1_lay(:)

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
        ALLOCATE(lay_map(1:nby,2))
        lay_map(:,:) = 0
        DO i = 0, nproc-1
          proc_map(i)%colsw(:) = 0
          proc_map(i)%colne(:) = 0
        END DO
!
        IF (ALLOCATED(nctot_lay)) DEALLOCATE(nctot_lay)
        IF (ALLOCATED(ncfl1_lay)) DEALLOCATE(ncfl1_lay)
        IF (ALLOCATED(nbl_lay)) DEALLOCATE(nbl_lay)
        ALLOCATE(nctot_lay(1:nby), ncfl1_lay(1:nby), nbl_lay(1:nby))
!
        nctot_lay = 0
        ncfl1_lay = 0
        nbl_lay = 0
!
! ... distribute cells among layers proportionally
! ... to the number of blocks (processors) contained
! ... (compute ncfl1_lay(layer))
!
        IF (rest .EQ. 0) THEN
          DO layer = 1, nby
           ncfl1_lay(layer) = localdim(countfl(1),nby,layer)
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
            ncfl1_lay(layer) = NINT(fact)
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

          nctot_lay(layer) = nctot_lay(layer) + 1

         IF ( icnt_layer .EQ. ncfl1_lay(layer) ) THEN
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
        nctot_lay(layer)  = 0
        ncfl1_lay(layer) = 0
        DO ij = ij1, ij2
          nctot_lay(layer) = nctot_lay(layer) + 1
          IF (fl(ij) .EQ. 1) ncfl1_lay(layer) = ncfl1_lay(layer) + 1 
        END DO
!
! ...   build block maps        
!
        bl1 = INT(ncfl1_lay(layer)/nbl_lay(layer))
        i = i1
        DO nbl = 1, nbl_lay(layer) 
          ipe = ipe + 1
          proc_map(ipe)%colsw(1) = i
          proc_map(ipe)%colsw(2) = j1
          DO WHILE (ncfl1(ipe) .LT. bl1)
            DO j = j1, j2
              ij = i + (j-1)*nr
              IF (fl(ij) .EQ. 1) ncfl1(ipe) = ncfl1(ipe) + 1
            END DO  
            i = i+1
          END DO
          proc_map(ipe)%colne(1) = i-1
          IF (nbl .EQ. nbl_lay(layer)) proc_map(ipe)%colne(1) = nr 
          proc_map(ipe)%colne(2) = j2
!
! ...     updates the number of cells with fl=1 into blocks
!
          nctot(ipe)  = 0
          ncfl1(ipe) = 0
          DO j = proc_map(ipe)%colsw(2), proc_map(ipe)%colne(2)
          DO i = proc_map(ipe)%colsw(1), proc_map(ipe)%colne(1)
            ij = i + (j-1)*nr
            nctot(ipe) = nctot(ipe) + 1
            IF (fl(ij) .EQ. 1) ncfl1(ipe) = ncfl1(ipe) + 1
          END DO
          END DO
          WRITE(7,*)'proc_map(',ipe,'):',proc_map(ipe)%colsw(:),proc_map(ipe)%colne(:)
        END DO
      END DO

      DEALLOCATE(lay_map)

      RETURN
      END SUBROUTINE
!
! ---------------------------------------------------------------------
!
      SUBROUTINE ghost
!
! ... Identifies and allocates ghost cells
! 
      USE dimensions
      USE parallel, ONLY: nproc, mpime
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
      INTEGER :: layer, j2, j1, i2, i1, nkt
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
!
      nset  = 0
      ncext = 0

      IF ( proc_map(mpime)%type == LAYER_MAP ) THEN
        DO ijk = proc_map(mpime)%lay(1), proc_map(mpime)%lay(2)
          IF ( fl( ijk ) == 1 ) THEN
            ncext = ncext + cell_neighbours( ijk, mpime, nset)
          END IF
        END DO
      ELSE IF ( proc_map(mpime)%type == BLOCK2D_MAP ) THEN
        i1 = proc_map(mpime)%colsw(1)
        i2 = proc_map(mpime)%colne(1)
        j1 = proc_map(mpime)%colsw(2)
        j2 = proc_map(mpime)%colne(2)
        IF( job_type == '2D' ) THEN
          DO j = j1, j2
            DO i = i1, i2
              ijk = i + nr*(j-1)
              IF ( fl(ijk) == 1 ) THEN
                ncext = ncext + cell_neighbours(ijk, mpime, nset)
              END IF
            END DO
          END DO
        ELSE IF( job_type == '3D' ) THEN
          DO k = 1, nz
            DO j = j1, j2
              DO i = i1, i2
                ijk = i + (j-1)*nx + (k-1)*nx*ny 
                IF ( fl(ijk) == 1 ) THEN
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
! ... ALLOCATE memory for the set of neighbouring cell on other proc.
!
      DO ipe = 0, nproc - 1
        IF( nset(ipe) > 0 ) THEN
          ALLOCATE( rcv_cell_set(ipe)%i(5, nset(ipe)) )
        ELSE
          NULLIFY( rcv_cell_set(ipe)%i )
        END IF
      END DO

! ... number of local cells
!
      ncint = nctot( mpime )
!
! ... allocate stencil array indexes
!
      ALLOCATE( myijk( nstdim , ncint ) )
      ALLOCATE( myinds( nstdim, ncint ) )
!
! ... set up indexes
!
      myijk   = 0
      myinds  = 0
      CALL indijk_setup()

! ... now fill in the receiving map and the indexes matrix
!
      nset  = 0
      icnt  = 0
      nrcv  = 0

      IF ( proc_map(mpime)%type == LAYER_MAP ) THEN
        DO ijk = proc_map(mpime)%lay(1), proc_map(mpime)%lay(2)
          IF ( fl(ijk) == 1 ) THEN
            icnt = icnt + cell_neighbours(ijk, mpime, nset, rcv_cell_set, myijk)
          END IF
        END DO
      ELSE IF ( proc_map(mpime)%type == BLOCK2D_MAP ) THEN
        i1 = proc_map(mpime)%colsw(1)
        i2 = proc_map(mpime)%colne(1)
        j1 = proc_map(mpime)%colsw(2)
        j2 = proc_map(mpime)%colne(2)
        IF( job_type == '2D' ) THEN
          DO j = j1, j2
            DO i = i1, i2
              ijk = i + nr*(j-1)
              IF ( fl(ijk) == 1 ) THEN
                icnt = icnt + cell_neighbours(ijk, mpime, nset, rcv_cell_set, myijk)
              END IF
            END DO
          END DO
        ELSE IF( job_type == '3D' ) THEN
          DO k = 1, nz
            DO j = j1, j2
              DO i = i1, i2
                ijk = i + (j-1)*nx + (k-1)*nx*ny
                IF ( fl(ijk) == 1 ) THEN
                  icnt = icnt + cell_neighbours(ijk, mpime, nset, rcv_cell_set, myijk)
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
! ... store the global cell index ijk, of the the local cell ijk
!
      DO ijk = 1, ncint
        myijk( ip0_jp0_kp0_, ijk) = cell_l2g(ijk, mpime)
      END DO
!
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
!
      ncext = SUM( nrcv ) 
      ncdom = ncint + ncext
      WRITE(7,* ) ' # ncext ', ncext
      WRITE(7,* ) ' # ncdom ', ncdom
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
! ... ALLOCATE and set the receiving map
!
      CALL set_rcv_map(rcv_map, myijk, rcv_cell_set, nset, nrcv)

! ... print out basic information on the map
!
      DO ipe = 0, nproc - 1
        WRITE(7,*) ' # nrcv ', nrcv(ipe), ' from ', ipe
      END DO
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

! ...   proc ipe send to other processor the number of elements 
! ...   he should receive from each of the other proc. (rcv_map(:)%nrcv)
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
!*****MX3D : ok
!
      DO ipe = 0, nproc - 1
        WRITE(7,100) rcv_map(ipe)%nrcv, ipe
        IF(rcv_map(ipe)%nrcv > 0 ) THEN
          WRITE(7,110) rcv_map(ipe)%ircv(:)
          WRITE(7,*) ' ---- '
          WRITE(7,110) rcv_map(ipe)%iloc(:)
        END IF
 100    FORMAT(' # receiving ',i5,' cells from ',i3)
 110    FORMAT(10i8)
      END DO

      DO ipe = 0, nproc - 1
        WRITE(7,200) snd_map(ipe)%nsnd, ipe
        IF (snd_map(ipe)%nsnd > 0 ) THEN
          WRITE(7,210) snd_map(ipe)%isnd(:)
          WRITE(7,*) ' ---- '
          WRITE(7,210) snd_map(ipe)%iloc(:)
        END IF
 200    FORMAT(' # sending ',i5,' cells to ',i3)
 210    FORMAT(10i8)
      END DO

      CALL test_comm
      
! ... fill in the array myinds using myijk
!
      ALLOCATE( fl_l(ncdom) )
      DO ijkl = 1, ncint
        ijk = myijk( ip0_jp0_kp0_, ijkl)
        fl_l(ijkl) = fl(ijk)
      END DO
!
      CALL data_exchange(fl_l)
      CALL set_myinds(myinds, myijk)
!
      RETURN

      END SUBROUTINE

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

          j = ( ijk - 1 ) / nr + 1
          i = MOD( ( ijk - 1 ), nr) + 1

        ELSE IF( job_type == '3D' ) THEN

          i = MOD( MOD( ijk - 1, nx*ny ), nx ) + 1
          j = MOD( ijk - 1, nx*ny ) / nx + 1
          k = ( ijk - 1 ) / ( nx*ny ) + 1

        ELSE

          CALL error(' cell_owner ', ' unknown job_type ', 1 )

        END IF
          
        DO ipe = 0, nproc - 1
          IF ( proc_map(ipe)%type == LAYER_MAP ) THEN 
            IF( ijk .GE. proc_map(ipe)%lay(1) .AND. ijk .LE. proc_map(ipe)%lay(2) ) THEN
              cell_owner = ipe
            END IF
          ELSE IF ( proc_map(ipe)%type == BLOCK2D_MAP ) THEN
            IF (j .GE. proc_map(ipe)%colsw(2) .AND. j .LE. proc_map(ipe)%colne(2)) THEN
              IF (i .GE. proc_map(ipe)%colsw(1) .AND.  &
                  i .LE. proc_map(ipe)%colne(1)) cell_owner = ipe
            END IF
          ELSE
            CALL error(' cell_owner ', ' partition type not yet implemented ', proc_map(ipe)%type )
          END IF
        END DO

      RETURN
      END FUNCTION     

!----------------------------------------------------------------------
!
!
      INTEGER FUNCTION cell_l2g( ijkl, mpime)

        USE dimensions
        USE control_flags, ONLY: job_type
        USE parallel, ONLY: nproc

        INTEGER, INTENT(IN) :: ijkl, mpime
        INTEGER :: i,j,k, i1,i2,j1,j2, nxl, nyl, nzl 
!
        IF ( proc_map(mpime)%type == LAYER_MAP ) THEN

          cell_l2g = ijkl + proc_map( mpime )%lay(1) - 1

        ELSE IF ( proc_map(mpime)%type == BLOCK2D_MAP ) THEN

          IF( job_type == '2D' ) THEN

            j1 = proc_map(mpime)%colsw(2)
            i1 = proc_map(mpime)%colsw(1)
            i2 = proc_map(mpime)%colne(1)
            i = MOD( ( ijkl - 1 ), (i2-i1+1)) + i1
            j = ( ijkl - 1 ) / (i2-i1+1) + j1
            cell_l2g = i + (j-1)*nr

          ELSE IF( job_type == '3D' ) THEN

            j1 = proc_map(mpime)%colsw(2)
            i1 = proc_map(mpime)%colsw(1)
            j2 = proc_map(mpime)%colne(2)
            i2 = proc_map(mpime)%colne(1)
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

        ELSE

          CALL error(' cell_l2g ', ' partition type not yet implemented ', proc_map(mpime)%type )

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
        INTEGER :: i, j, k, i1, i2, j1, j2, nxl, nyl, nzl
!
        IF ( proc_map(mpime)%type == LAYER_MAP ) THEN

          cell_g2l = ijk - proc_map(mpime)%lay(1) + 1

        ELSE IF ( proc_map(mpime)%type == BLOCK2D_MAP ) THEN

          IF( job_type == '2D' ) THEN

            j = ( ijk - 1 ) / nr + 1
            i = MOD( ( ijk - 1 ), nr) + 1
            j1 = proc_map(mpime)%colsw(2)
            i1 = proc_map(mpime)%colsw(1)
            i2 = proc_map(mpime)%colne(1)
            cell_g2l = (i-i1+1) + (j-j1)*(i2-i1+1)

          ELSE IF( job_type == '3D' ) THEN

            i = MOD( MOD( ijk - 1, nx*ny ), nx ) + 1
            j = MOD( ijk - 1, nx*ny ) / nx + 1
            k = ( ijk - 1 ) / ( nx*ny ) + 1

            j1 = proc_map(mpime)%colsw(2)
            i1 = proc_map(mpime)%colsw(1)
            j2 = proc_map(mpime)%colne(2)
            i2 = proc_map(mpime)%colne(1)

            nxl = i2 - i1 + 1
            nyl = j2 - j1 + 1
            nzl = nz

            cell_g2l = 1 + (i-i1) + (j-j1)*nxl + (k-1)*nxl*nyl

          ELSE

            CALL error(' cell_owner ', ' unknown job_type ', 1 )

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
        INTEGER :: ippj, ijpp, ii, jj, kk
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

          j = ( ijk - 1 ) / nr + 1
          i = MOD( ( ijk - 1 ), nr) + 1
        
! ...     loop over the first neighbouring cells
          DO im = -2, 2
            DO jm = -2, 2
              IF( ( ABS( im ) + ABS( jm ) ) <= 2 ) THEN
                IF( (im /= 0) .OR. (jm /= 0) ) THEN
                  ii = im
                  jj = jm 
                  IF( ( i == 2    ) .AND. ( ii == -2 ) ) ii = -1
                  IF( ( i == nr-1 ) .AND. ( ii == +2 ) ) ii = +1
                  IF( ( j == 2    ) .AND. ( jj == -2 ) ) jj = -1
                  IF( ( j == nz-1 ) .AND. ( jj == +2 ) ) jj = +1
                  ijke = ijk + ii + jj * nr
                  ipe =  cell_owner(ijke) 
                  IF( ipe /= mpime) THEN
! ...               the cell ijke is not local, count and register its position
                    nset(ipe) = nset(ipe) + 1
                    icnt = icnt + 1
                    IF( fill ) THEN
                      rcv_cell_set(ipe)%i(1,nset(ipe)) = ijke
                      rcv_cell_set(ipe)%i(2,nset(ipe)) = ijkl
                      rcv_cell_set(ipe)%i(3,nset(ipe)) = im
                      rcv_cell_set(ipe)%i(4,nset(ipe)) = jm
                      rcv_cell_set(ipe)%i(5,nset(ipe)) = 0
                    END IF
                  ELSE IF( fill ) THEN
! ...               the cell ijke is local, set the mapping with cell ijkl
                    ijkel = cell_g2l(ijke, mpime)
                    myijk( indijk( im, jm, 0), ijkl ) = ijkel
                  END IF
                ELSE IF( fill ) THEN   
! ...             store the global cell index ij, of the the local cell ijkl
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
      END SUBROUTINE
!
!----------------------------------------------------------------------
!
      SUBROUTINE set_rcv_map(rcv_map, myijk, rcv_cell_set, nset, nrcv)
!
        USE basic_types, ONLY: imatrix 
        USE control_flags, ONLY: job_type
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

            ic        = 1           ! index that run over neighbour cells from the current processor 
            icnt_rcv  = 1           ! counts received cells from the current processor (ipe)
            icnt      = icnt + 1    ! increment the number of cell

            ijke = rcv_cell_set(ipe)%i(1,ic)
            ijkl = rcv_cell_set(ipe)%i(2,ic)
            i    = rcv_cell_set(ipe)%i(3,ic)
            j    = rcv_cell_set(ipe)%i(4,ic)
            k    = rcv_cell_set(ipe)%i(5,ic)

            rcv_map(ipe)%ircv(icnt_rcv)  = ijke
            rcv_map(ipe)%iloc(icnt_rcv)  = icnt + ncint
            myijk( indijk(i,j,k), ijkl ) = icnt + ncint

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
                rcv_map(ipe)%iloc(icnt_rcv) = icnt + ncint
              END IF 

              myijk( indijk(i,j,k), ijkl ) = icnt + ncint  

            END DO

            DEALLOCATE( ijksort )

            IF( icnt_rcv /= nrcv(ipe) ) &
              CALL error( '  set_rcv_map ', ' inconsistent numbero of cells ', icnt_rcv )

          ELSE

            NULLIFY( rcv_map(ipe)%ircv )
            NULLIFY( rcv_map(ipe)%iloc )

          END IF

        END DO

        RETURN
      END SUBROUTINE

!---------------------------------------------------
!
      SUBROUTINE set_myinds(myinds, myijk)
!
        USE dimensions
        USE control_flags, ONLY: job_type
!
        IMPLICIT NONE
        INTEGER :: myinds(:,:)
        INTEGER :: myijk(:,:)
!
        INTEGER :: ijl, ijr, ijb, ijt
        INTEGER :: ijbr, ijtr, ijbl, ijtl, ijrr, ijtt, ijll, ijbb
        INTEGER :: imj, ipj, ijm, ijp
        INTEGER :: ipjm, ipjp, imjm, imjp, ippj, ijpp, immj, ijmm

        INTEGER :: nflr, nflt, nfll, nflb
        INTEGER :: nfltr, nfltl, nflbr, nflbl
        INTEGER :: nflrr, nfltt, nflll, nflbb
        INTEGER :: i, j, k, ijk, imesh

        INTEGER ::  ipjk, imjk, ippjk, immjk, ijpk, ipjpk, imjpk, ijmk,  &
                    ipjmk, imjmk, ijppk, ijmmk, ijkp, ipjkp, imjkp, ijpkp,  &
                    ijmkp, ijkm, ipjkm, imjkm, ijpkm, ijmkm, ijkpp, ijkmm

        INTEGER ::  ijke, ijkw, ijkee, ijkww, ijkn, ijken, ijkwn, ijks,  &
                    ijkes,  ijkws, ijknn, ijkss, ijkt, ijket, ijkwt, ijknt, &
                    ijkst, ijkb, ijkeb, ijkwb, ijknb, ijksb, ijktt, ijkbb
!
        DO ijk = 1, ncint 

          imesh = myijk( ip0_jp0_kp0_, ijk)

          IF( job_type == '2D' ) THEN

            j  = ( imesh - 1 ) / nr + 1
            i  = MOD( ( imesh - 1 ), nr) + 1
            IF( (i .GE. 2) .AND. (i .LE. (nr-1)) .AND.   &
                (j .GE. 2) .AND. (j .LE. (nz-1))      ) THEN
!
              ijm = myijk( ip0_jm1_kp0_, ijk)
              imj = myijk( im1_jp0_kp0_, ijk)
              ipj = myijk( ip1_jp0_kp0_, ijk)
              ijp = myijk( ip0_jp1_kp0_, ijk)
              ipjm= myijk( ip1_jm1_kp0_, ijk)
              ipjp= myijk( ip1_jp1_kp0_, ijk)
              imjm= myijk( im1_jm1_kp0_, ijk)
              imjp= myijk( im1_jp1_kp0_, ijk)
              ippj= myijk( ip2_jp0_kp0_, ijk)
              ijpp= myijk( ip0_jp2_kp0_, ijk)
              immj= myijk( im2_jp0_kp0_, ijk)
              ijmm= myijk( ip0_jm2_kp0_, ijk)
!
! ... First neighbours: near the axis or solid boundaries 
! ... impose homogeneous Neumann conditions
!
              ijr = ipj
              nflr = fl_l(ipj) 
              IF((nflr.EQ.2).OR.(nflr.EQ.3)) ijr = ijk

              ijl = imj
              nfll = fl_l(imj) 
              IF((nfll.EQ.2).OR.(nfll.EQ.3)) ijl = ijk

              ijt = ijp
              nflt = fl_l(ijp) 
              IF((nflt.EQ.2).OR.(nflt.EQ.3)) ijt = ijk

              ijb = ijm
              nflb = fl_l(ijm) 
              IF((nflb.EQ.2).OR.(nflb.EQ.3)) ijb = ijk

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
                    ijtr = ijk 
                  ELSE
                    ijtr = ipj
                  END IF
                ELSE
                  IF (nflr.EQ.2 .OR. nflr.EQ.3) THEN
                    ijtr = ijp 
                  ELSE
                    ijtr = ijk
                  END IF
                END IF
              END IF
!
! ... Second neighbours are not available on boundaries
!
              ijrr = ippj
              IF(i .EQ. (nr-1)) ijrr = ijr
              nflrr=fl_l(ippj)
              IF( (nflrr.EQ.2) .OR. (nflrr.EQ.3) ) ijrr = ipj

              ijtt = ijpp
              IF(j.EQ.(nz-1)) ijtt = ijt
              nfltt=fl_l(ijpp)
              IF( (nfltt.EQ.2) .OR. (nfltt.EQ.3) ) ijtt = ijp
  
              ijll = immj
              IF(i .EQ. (2)) ijll = ijl
              nflll=fl_l(immj)
              IF( (nflll.EQ.2) .OR. (nflll.EQ.3) ) ijll = imj

              ijbb = ijmm
              IF(j.EQ.(2)) ijbb = ijb
              nflbb=fl_l(ijmm)
              IF( (nflbb.EQ.2) .OR. (nflbb.EQ.3) ) ijbb = ijm
!
              myinds(ip0_jm2_kp0_, ijk) = ijbb
              myinds(im1_jm1_kp0_, ijk) = ijbl
              myinds(ip0_jm1_kp0_,  ijk) = ijb
              myinds(ip1_jm1_kp0_, ijk) = ijbr
              myinds(im2_jp0_kp0_, ijk) = ijll
              myinds(im1_jp0_kp0_,  ijk) = ijl
              myinds(ip1_jp0_kp0_,  ijk) = ijr
              myinds(ip2_jp0_kp0_, ijk) = ijrr
              myinds(im1_jp1_kp0_, ijk) = ijtl
              myinds(ip0_jp1_kp0_,  ijk) = ijt
              myinds(ip1_jp1_kp0_, ijk) = ijtr
              myinds(ip0_jp2_kp0_, ijk) = ijtt

            END IF

          ELSE IF( job_type == '3D' ) THEN

            i = MOD( MOD( imesh - 1, nx*ny ), nx ) + 1
            j = MOD( imesh - 1, nx*ny ) / nx + 1
            k = ( imesh - 1 ) / ( nx*ny ) + 1

            IF( (i >= 2) .AND. (i <= (nx-1)) .AND.   &
                (j >= 2) .AND. (j <= (ny-1)) .AND.   &
                (k >= 2) .AND. (k <= (nz-1))         ) THEN
  
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
  
              ijke  =  ipjk
              ijkw  =  imjk
              ijkee =  ippjk
              ijkww =  immjk
              ijkn  =  ijpk
              ijken =  ipjpk
              ijkwn =  imjpk
              ijks  =  ijmk
              ijkes =  ipjmk
              ijkws =  imjmk
              ijknn =  ijppk
              ijkss =  ijmmk
              ijkt  =  ijkp
              ijket =  ipjkp
              ijkwt =  imjkp
              ijknt =  ijpkp
              ijkst =  ijmkp
              ijkb  =  ijkm
              ijkeb =  ipjkm
              ijkwb =  imjkm
              ijknb =  ijpkm
              ijksb =  ijmkm
              ijktt =  ijkpp
              ijkbb =  ijkmm
  
              if( fl_l( ipjk  ) == 2 .OR. fl_l( ipjk )  == 3 ) ijke = ijk
              if( fl_l( imjk  ) == 2 .OR. fl_l( imjk )  == 3 ) ijkw = ijk
              if( (i /= (nx-1)) .AND. ( fl_l( ippjk ) == 2 .OR. fl_l( ippjk ) == 3 ) ) ijkee = ijke
              if( (i == (nx-1)) ) ijkee = ijke
              if( (i /= 2) .AND. ( fl_l( immjk ) == 2 .OR. fl_l( immjk ) == 3 ) ) ijkww = ijkw
              if( (i == 2) ) ijkww = ijkw
              if( fl_l( ijpk ) == 2 .OR. fl_l( ijpk ) == 3 ) ijkn = ijk
              ! check corners   ijken ijkwn
              if( fl_l( ijmk ) == 2 .OR. fl_l( ijmk ) == 3 ) ijks = ijk
              ! check corners   ijkes ijkws
              if( (j /= (ny-1)) .AND. ( fl_l( ijppk ) == 2 .OR. fl_l( ijppk ) == 3 ) ) ijknn = ijkn
              if( (j == (ny-1)) ) ijknn = ijkn
              if( (j /= 2) .AND. ( fl_l( ijmmk ) == 2 .OR. fl_l( ijmmk ) == 3 ) ) ijkss = ijks
              if( (j == 2) ) ijkss = ijks
              if( fl_l( ijkp ) == 2 .OR. fl_l( ijkp ) == 3 ) ijkt = ijk
              ! check corners   ijket ijkwt ijknt ijkst
              if( fl_l( ijkm ) == 2 .OR. fl_l( ijkm ) == 3 ) ijkb = ijk
              ! check corners   ijkeb ijkwb ijknb ijksb
              if( (k /= (nz-1)) .AND. ( fl_l( ijkpp ) == 2 .OR. fl_l( ijkpp ) == 3 ) ) ijktt = ijkt
              if( (k == (nz-1)) ) ijktt = ijkt
              if( (k /= 2) .AND. ( fl_l( ijkmm ) == 2 .OR. fl_l( ijkmm ) == 3 ) ) ijkbb = ijkb
              if( (k == 2) ) ijkbb = ijkb

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
      END SUBROUTINE
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
        INTEGER :: ip, isour, idest, ib
!
        DO ip = 1, (nproc - 1)
          isour = MOD(mpime - ip + nproc, nproc)
          idest = MOD(mpime + ip        , nproc)
          ALLOCATE( rcvbuf( MAX(rcv_map(isour)%nrcv,1) ) )
          ALLOCATE( sndbuf( MAX(snd_map(idest)%nsnd,1) ) )
          DO ib = 1, snd_map(idest)%nsnd
            sndbuf(ib) = array( snd_map(idest)%iloc(ib) )
          END DO 
          CALL sendrecv_real(sndbuf(1), snd_map(idest)%nsnd, idest,  &
            rcvbuf, rcv_map(isour)%nrcv, isour, ip)
          DO ib = 1, rcv_map(isour)%nrcv
            array( rcv_map(isour)%iloc(ib) ) = rcvbuf(ib)
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
          sdim = MAX( SIZE(array,2) * snd_map(idest)%nsnd, 1)
          rdim = MAX( SIZE(array,2) * rcv_map(isour)%nrcv, 1)
          ALLOCATE( rcvbuf( rdim ) )
          ALLOCATE( sndbuf( sdim ) )
          DO ik = 1, SIZE(array,2)
            DO ib = 1, snd_map(idest)%nsnd
              sndbuf(ib + snd_map(idest)%nsnd*(ik-1)) = array( snd_map(idest)%iloc(ib), ik )
            END DO
          END DO 
          CALL sendrecv_real(sndbuf(1), SIZE(sndbuf), idest,     &      
                             rcvbuf(1), SIZE(rcvbuf), isour, ip)
          DO ik = 1, SIZE(array,2)
            DO ib = 1, rcv_map(isour)%nrcv
              array( rcv_map(isour)%iloc(ib), ik ) = rcvbuf( ib + rcv_map(isour)%nrcv*(ik-1))
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
            sndbuf(ib) = array( snd_map(idest)%iloc(ib) )
          END DO
          CALL sendrecv_integer(sndbuf, snd_map(idest)%nsnd, idest,      &
                                rcvbuf, rcv_map(isour)%nrcv, isour, ip)
          DO ib = 1, rcv_map(isour)%nrcv
            array( rcv_map(isour)%iloc(ib) ) = rcvbuf(ib)
          END DO
          DEALLOCATE( rcvbuf )
          DEALLOCATE( sndbuf )
        END DO
        RETURN
      END SUBROUTINE


      SUBROUTINE data_collect_r( garray, larray, imstart, imend )

        ! this subroutine collect data distributed across processors
        ! and store them in the array "garray"
        ! The subroutine is designed to allow the collection of global
        ! data using subarray of small size, useful when available memory
        ! is not enough to store a global array

        USE parallel, ONLY: nproc, mpime
        USE control_flags, ONLY: job_type
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
      END SUBROUTINE


      SUBROUTINE data_distribute_r( garray, larray, imstart, imend )

        ! this subroutine distribute a global array "garray"  across processors
        ! and store them in the local array "larray"
        ! The subroutine is designed to allow the distribution of global
        ! data using subarray of small size, useful when available memory
        ! is not enough to store a global array

        USE parallel, ONLY: nproc, mpime
        USE control_flags, ONLY: job_type
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
      END SUBROUTINE



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
         
        WRITE(7,*)   ' index received ' 
        WRITE(7,310) itest( ncint + 1 : ncdom ) 
 310    FORMAT(10i8)
        
        DEALLOCATE( itest )
        
        RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      END MODULE
!----------------------------------------------------------------------
