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
        INTEGER :: nbx,nby             ! number of blocks in x and y directions
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
        INTEGER, DIMENSION(:), ALLOCATABLE :: fl

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
        IF( job_type == '2D' ) THEN

          IF ( fl((nr-1)+(nz-1)*nr) == 4 .AND. fl(nr+(nz-2)*nr) == 4) THEN
            fl(nr+(nz-1)*nr) = 4
          END IF

        ELSE IF( job_type == '3D' ) THEN

          CALL error( ' flic ',' set the upper-right corner for X3D ', 1)

        END IF
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
                DO i = iob(n)%rlo, iob(n)%rhi
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

        RETURN
      END SUBROUTINE
!
!----------------------------------------------------------------------

! ... MODIFICARE_X3D (rivedere fino alla fine del file)

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
      INTEGER :: i, j, ij, n, m, q, rr, k
      INTEGER :: nfl, ipe, dim
      INTEGER :: layer
      INTEGER, ALLOCATABLE :: red(:), green(:), blue(:)
      INTEGER :: localdim
!
! ... Subroutine Body
!
      ALLOCATE( red(nr*nz), green(nr*nz), blue(nr*nz) )

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
      END DO
      WRITE(7,*) '  '
!

      IF( outppm ) THEN

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

      END IF

      DEALLOCATE(red,green,blue)

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
      INTEGER :: ncfl1(:)
      INTEGER :: nctot(:)
!
      INTEGER :: i, j, ij, ijk, k
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

            ij = i + (j-1)*nr
  
            IF ( fl(ij) == 1 ) THEN
              icnt_ipe = icnt_ipe + 1
              icnt     = icnt     + 1
            END IF

            nctot(ipe) = nctot(ipe) + 1

            IF ( ( icnt_ipe == ncfl1(ipe) ) .AND. ( i /= (nr-1) ) ) THEN
              icnt_ipe = 0
              IF( icnt < countfl(1) ) THEN
                lay_map(ipe,2) = ij
                ipe = ipe + 1
                lay_map(ipe,1) = ij+1
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

              IF ( ( icnt_ipe == ncfl1(ipe) ) .AND. &
                   ( i /= (nx-1) ) .AND. ( j /= (ny-1) ) ) THEN
                icnt_ipe = 0
                IF( icnt < countfl(1) ) THEN
                  lay_map(ipe,2) = ijk
                  ipe = ipe + 1
                  lay_map(ipe,1) = ijk+1
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
      INTEGER :: ncfl1(:)
      INTEGER :: nctot(:)
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
!
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
!
!*****MX3D: il loop sulle celle interne va fatto con le nuove maps
!
      nset   = 0
      ncext = 0
      IF ( proc_map(mpime)%type == LAYER_MAP ) THEN
        DO ij = proc_map(mpime)%lay(1), proc_map(mpime)%lay(2)
          IF ( fl(ij) .EQ. 1 ) THEN
            ncext = ncext + cell_neighbours(ij, mpime, nset)
          END IF
        END DO
      ELSE IF ( proc_map(mpime)%type == BLOCK2D_MAP ) THEN
        i1 = proc_map(mpime)%colsw(1)
        i2 = proc_map(mpime)%colne(1)
        j1 = proc_map(mpime)%colsw(2)
        j2 = proc_map(mpime)%colne(2)
        DO j = j1, j2
         DO i = i1, i2
          ij = i + nr*(j-1)
          IF (fl(ij).EQ.1) THEN
            ncext = ncext + cell_neighbours(ij, mpime, nset)
          END IF
         END DO
        END DO
      ELSE
        CALL error(' ghost ', ' partition type not yet implemented ', proc_map(mpime)%type )
      END IF
!
! ... ALLOCATE memory for the set of neighbouring cell on other proc.
!
!*****MX3D : cambia la prima dimensione di rcv_cell_set (4 -> 5) 
!
      DO ipe = 0, nproc - 1
        IF(nset(ipe) .GT. 0 ) THEN
          ALLOCATE( rcv_cell_set(ipe)%i(4, nset(ipe)) )
        ELSE
          NULLIFY( rcv_cell_set(ipe)%i )
        END IF
      END DO

! ... number of local cells
      ncint = nctot(mpime)
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
      nset  = 0
      icnt  = 0
      nrcv  = 0
      IF ( proc_map(mpime)%type == LAYER_MAP ) THEN
!
!*****MX3D : cambia il loop sulle celle interne (v.sopra)
!
        DO ij = proc_map(mpime)%lay(1), proc_map(mpime)%lay(2)
          IF ( fl(ij) .EQ. 1 ) THEN
          icnt = icnt + cell_neighbours(ij, mpime, nset, rcv_cell_set, myijk)
          END IF
        END DO
      ELSE IF ( proc_map(mpime)%type == BLOCK2D_MAP ) THEN
        i1 = proc_map(mpime)%colsw(1)
        i2 = proc_map(mpime)%colne(1)
        j1 = proc_map(mpime)%colsw(2)
        j2 = proc_map(mpime)%colne(2)
       DO j = j1, j2
        DO i = i1, i2
          ij = i + nr*(j-1)
          IF (fl(ij).EQ.1) THEN
          icnt = icnt + cell_neighbours(ij, mpime, nset, rcv_cell_set, myijk)
          END IF
        END DO
       END DO
      ELSE
        CALL error(' ghost ', ' partition type not yet implemented ', proc_map(mpime)%type )
      END IF
!
!
      DO ijl = 1, ncint
! ...     store the global cell index ij, of the the local cell ijl
          myijk( ip0_jp0_kp0_, ijl) = cell_l2g(ijl, mpime)
      END DO

!
!*****MX3D : ok
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
      ncext = SUM( nrcv ) 
      ncdom = ncint + ncext
      WRITE(7,* ) ' # ncext ', ncext
      WRITE(7,* ) ' # ncdom ', ncdom

!
!*****MX3D : ok
!
! ... prepare the receive map 
      DO ipe = 0, nproc - 1
        rcv_map(ipe)%nrcv = nrcv(ipe)
        NULLIFY( rcv_map(ipe)%ircv )
        NULLIFY( rcv_map(ipe)%iloc )
      END DO

!
!*****MX3D : ok
!
! ... prepare the send map 
      DO ipe = 0, nproc - 1
        NULLIFY( snd_map(ipe)%isnd ) 
        NULLIFY( snd_map(ipe)%iloc ) 
      END DO

!
!*****MX3D : piccole modifiche a set_rcv_map
!
! ... ALLOCATE and set the receiving map
      CALL set_rcv_map(rcv_map, myijk, rcv_cell_set, nset, nrcv)

! ... print out basic information on the map
      DO ipe = 0, nproc - 1
        WRITE(7,*) ' # nrcv ', nrcv(ipe), ' from ', ipe
      END DO
! ... using the receive maps fill in the sending maps
!
!*****MX3D : ok
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
        IF( ipe .EQ. mpime ) THEN
          nrcv(:) = rcv_map(:)%nrcv
          nrcvx = MAXVAL(nrcv)
        END IF
! ...   each proc store the number of element he should send to ipe
!
!*****MX3D : ok (rivedere)
!
        CALL scatter_integer(nrcv, snd_map(ipe)%nsnd, 1, ipe)
        CALL bcast_integer(nrcvx, 1, ipe)

! ...   proc ipe send to other processors the array with the indexes of 
! ...   the elements he should receive ( rcv_map(jpe)%ircv )
!
        ALLOCATE(ijrcv(nrcvx,0:nproc-1))
        ALLOCATE(ijsnd(nrcvx))
        IF( ipe .EQ. mpime) THEN
          ijrcv = 0
          DO jpe = 0, nproc - 1
            IF( rcv_map(jpe)%nrcv .GT. 0 ) THEN
              ijrcv(1:rcv_map(jpe)%nrcv, jpe ) = rcv_map(jpe)%ircv(1:rcv_map(jpe)%nrcv) 
            END IF
          END DO 
        END IF 
!
! ...   each proc store the index of the elements he should send to ipe
!
        CALL scatter_integer(ijrcv, ijsnd, nrcvx, ipe)
        IF( ipe .NE. mpime) THEN
          ALLOCATE(snd_map(ipe)%isnd(snd_map(ipe)%nsnd))
          ALLOCATE(snd_map(ipe)%iloc(snd_map(ipe)%nsnd))
          snd_map(ipe)%isnd(1:snd_map(ipe)%nsnd) = ijsnd(1:snd_map(ipe)%nsnd)
          DO ij = 1, snd_map(ipe)%nsnd
            snd_map(ipe)%iloc( ij ) = cell_g2l( ijsnd( ij ), mpime )
          END DO
        END IF 
!
          DEALLOCATE(ijrcv)
          DEALLOCATE(ijsnd)

      END DO
!
!*****MX3D : ok
!
      DO ipe = 0, nproc - 1
        WRITE(7,100) rcv_map(ipe)%nrcv, ipe
        IF(rcv_map(ipe)%nrcv .GT. 0 ) THEN
          WRITE(7,110) rcv_map(ipe)%ircv(:)
          WRITE(7,*) ' ---- '
          WRITE(7,110) rcv_map(ipe)%iloc(:)
        END IF
 100    FORMAT(' # receiving ',i5,' cells from ',i3)
 110    FORMAT(10i8)
      END DO

      DO ipe = 0, nproc - 1
        WRITE(7,200) snd_map(ipe)%nsnd, ipe
        IF (snd_map(ipe)%nsnd .GT. 0 ) THEN
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
!*****MX3D : ok
!
      ALLOCATE(fl_l(ncdom))
      DO ijl = 1, ncint
        ij = myijk( ip0_jp0_kp0_, ijl)
        fl_l(ijl) = fl(ij)
      END DO
      CALL data_exchange(fl_l)
      CALL set_myinds(myinds, myijk)
!
! ... Test the indexing (following sweep direction)
!      DO ijl = 1, ncint
!        WRITE(7,*) ijl
!        WRITE(7,*) (myinds(k,ijl), k=1,12)
! *CHANGE  (i,j,ijl) => (indijk(i,j,0),ijl) 
!        WRITE(7,*) ((myijk(i,j,ijl),i=-2,2),j=-2,2)
!      END DO
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
!
!*****MX3D : cambiare secondo le nuove maps
!
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
          IF ( proc_map(ipe)%type == LAYER_MAP ) THEN 
            IF( ij .GE. proc_map(ipe)%lay(1) .AND. ij .LE. proc_map(ipe)%lay(2) ) THEN
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
!*****MX3D : cambiare secondo le nuove maps
!
      INTEGER FUNCTION cell_l2g(ijl, mpime)
      USE dimensions
      USE parallel, ONLY: nproc
      INTEGER, INTENT(IN) :: ijl, mpime
      INTEGER :: i,j,i1,i2,j1
!
      IF ( proc_map(mpime)%type == LAYER_MAP ) THEN
        cell_l2g = ijl + proc_map(mpime)%lay(1) - 1
      ELSE IF ( proc_map(mpime)%type == BLOCK2D_MAP ) THEN
        j1 = proc_map(mpime)%colsw(2)
        i1 = proc_map(mpime)%colsw(1)
        i2 = proc_map(mpime)%colne(1)
        i = MOD( ( ijl - 1 ), (i2-i1+1)) + i1
        j = ( ijl - 1 ) / (i2-i1+1) + j1
        cell_l2g = i + (j-1)*nr
      ELSE
        CALL error(' cell_l2g ', ' partition type not yet implemented ', proc_map(mpime)%type )
      END IF
!
      RETURN
      END FUNCTION     
!----------------------------------------------------------------------
!
!*****MX3D : cambiare secondo le nuove maps
!
      INTEGER FUNCTION cell_g2l(ij, mpime)
      USE dimensions
      USE parallel, ONLY: nproc
      INTEGER, INTENT(IN) :: ij, mpime
      INTEGER :: i,j,i1,i2,j1
!
      IF ( proc_map(mpime)%type == LAYER_MAP ) THEN
        cell_g2l = ij - proc_map(mpime)%lay(1) + 1
      ELSE IF ( proc_map(mpime)%type == BLOCK2D_MAP ) THEN
        j = ( ij - 1 ) / nr + 1
        i = MOD( ( ij - 1 ), nr) + 1
        j1 = proc_map(mpime)%colsw(2)
        i1 = proc_map(mpime)%colsw(1)
        i2 = proc_map(mpime)%colne(1)
        cell_g2l = (i-i1+1) + (j-j1)*(i2-i1+1)
      ELSE
        CALL error(' cell_g2l ', ' partition type not yet implemented ', proc_map(mpime)%type )
      END IF
!
      RETURN
      END FUNCTION     
!----------------------------------------------------------------------
!
!*****MX3D : aggiungere la terza dimensione nei loop
!
      INTEGER FUNCTION cell_neighbours(ij, mpime, nset, rcv_cell_set, myijk)
        USE dimensions
        USE basic_types, ONLY: imatrix
        INTEGER, INTENT(IN) :: ij, mpime
        INTEGER, INTENT(INOUT) :: nset(0:)
        TYPE (imatrix), OPTIONAL :: rcv_cell_set(0:)
        INTEGER, OPTIONAL :: myijk(:,:)
        INTEGER :: ippj, ijpp, ijmm, immj
        INTEGER :: icnt, ipe, i, j, ije, ijl, ijel
        INTEGER :: im, jm
        LOGICAL :: fill

        IF ( PRESENT(rcv_cell_set) .AND. .NOT. present(myijk) ) THEN
          WRITE(8,*) ' error in neighbour - 1 '
          STOP 
        END IF
        IF( .NOT. PRESENT(rcv_cell_set) .AND. present(myijk) ) THEN
          WRITE(8,*) ' error in neighbour - 2 '
          STOP 
        END IF

        fill = PRESENT(rcv_cell_set)
        icnt = 0
        IF( fill ) THEN
          ncint = SIZE(myijk, 2)
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
                myijk( indijk(im,jm,0), ijl ) = ijel
              END IF
            ELSE IF( fill ) THEN   
! ...         store the global cell index ij, of the the local cell ijl
              myijk( ip0_jp0_kp0_, ijl) = ij
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
          myijk( im2_jp0_kp0_, ijl) = ijel
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
          myijk( ip2_jp0_kp0_, ijl ) = ijel
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
          myijk( ip0_jm2_kp0_, ijl ) = ijel
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
          myijk( ip0_jp2_kp0_, ijl ) = ijel
        END IF

        cell_neighbours = icnt

        RETURN
      END FUNCTION
!----------------------------------------------------------------------
!
!*****MX3D : ok
!
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
!
!*****MX3D : aggiungere la terza dimensione
!
      SUBROUTINE set_rcv_map(rcv_map, myijk, rcv_cell_set, nset, nrcv)
        USE basic_types, ONLY: imatrix 
        TYPE (rcv_map_type) :: rcv_map(0:)
        TYPE (imatrix) :: rcv_cell_set(0:)
        INTEGER :: myijk(:,:)
        INTEGER, INTENT(IN) :: nrcv(0:), nset(0:)
        INTEGER :: ipe, nproc, icnt, ijl, i, j, ije
        INTEGER :: ic, icnt_rcv
        INTEGER, ALLOCATABLE :: ijsort(:)

        nproc = SIZE(rcv_cell_set)
        ncint = SIZE(myijk,2)
        icnt  = 0                    ! counts received cells from all processors

        DO ipe = 0, nproc - 1

          rcv_map(ipe)%nrcv = nrcv(ipe)
          IF( ASSOCIATED(rcv_map(ipe)%ircv) ) DEALLOCATE( rcv_map(ipe)%ircv )
          IF( ASSOCIATED(rcv_map(ipe)%iloc) ) DEALLOCATE( rcv_map(ipe)%iloc )

          IF( nrcv(ipe) .GT. 0 ) THEN
            ALLOCATE( rcv_map(ipe)%ircv(nrcv(ipe)) )
            ALLOCATE( rcv_map(ipe)%iloc(nrcv(ipe)) )

            ALLOCATE( ijsort(nset(ipe)) )
            ijsort(:) = rcv_cell_set(ipe)%i(1,:)

            ic        = 1           ! index that run over neighbour cells from the current processor 
            icnt_rcv  = 1           ! counts received cells from the current processor (ipe)
            icnt      = icnt + 1
            ije = rcv_cell_set(ipe)%i(1,ic)
            ijl = rcv_cell_set(ipe)%i(2,ic)
            i   = rcv_cell_set(ipe)%i(3,ic)
            j   = rcv_cell_set(ipe)%i(4,ic)
            rcv_map(ipe)%ircv(icnt_rcv) = ije
            rcv_map(ipe)%iloc(icnt_rcv) = icnt + ncint
            myijk( indijk(i,j,0), ijl ) = icnt + ncint
            DO ic = 2, nset(ipe)
              ije = rcv_cell_set(ipe)%i(1,ic)
              ijl = rcv_cell_set(ipe)%i(2,ic)
              i   = rcv_cell_set(ipe)%i(3,ic)
              j   = rcv_cell_set(ipe)%i(4,ic)
              ! ... k  = rcv_cell_set(ipe)%i(5,ic)    ! aggiungere
              ! ... stind = stencil_index( i, j, k )  ! = a number from 0 and 24, aggiungere
              IF( ijsort(ic) .NE. ijsort(ic-1) ) THEN
                icnt_rcv  = icnt_rcv + 1
                icnt      = icnt + 1
                rcv_map(ipe)%ircv(icnt_rcv) = ije
                rcv_map(ipe)%iloc(icnt_rcv) = icnt + ncint
              END IF 
              myijk( indijk(i,j,0), ijl ) = icnt + ncint  ! da sostituire con
              ! ... myijk( stind, ijl ) = icnt + ncint
            END DO

            DEALLOCATE( ijsort)

            IF( icnt_rcv .NE. nrcv(ipe) ) THEN
              WRITE(8,*) ' error in set_rcv_map - 1 '
              STOP
            END IF

          ELSE
            NULLIFY(rcv_map(ipe)%ircv)
            NULLIFY(rcv_map(ipe)%iloc)
          END IF

        END DO

        RETURN
      END SUBROUTINE

!---------------------------------------------------
!
!*****MX3D : cambiare gli stencil
!
      SUBROUTINE set_myinds(myinds, myijk)
        USE dimensions
        IMPLICIT NONE
        INTEGER :: myinds(:,:)
        INTEGER :: myijk(:,:)
!
        INTEGER :: ijl, ijr, ijb, ijt, ijbr, ijtr, ijbl, ijtl, ijrr, ijtt, ijll, ijbb
        INTEGER :: imj, ipj, ijm, ijp, ipjm, ipjp, imjm, imjp, ippj, ijpp, immj, ijmm

        INTEGER :: nflr, nflt, nfll, nflb
        INTEGER :: nfltr, nfltl, nflbr, nflbl
        INTEGER :: nflrr, nfltt, nflll, nflbb
        INTEGER :: i, j, ij, imesh
!
        DO ij = 1, ncint 
          imesh = myijk( ip0_jp0_kp0_, ij)
          j  = ( imesh - 1 ) / nr + 1
          i  = MOD( ( imesh - 1 ), nr) + 1
          IF( (i .GE. 2) .AND. (i .LE. (nr-1)) .AND.   &
              (j .GE. 2) .AND. (j .LE. (nz-1))      ) THEN
!
            ijm = myijk( ip0_jm1_kp0_, ij)
            imj = myijk( im1_jp0_kp0_, ij)
            ipj = myijk( ip1_jp0_kp0_, ij)
            ijp = myijk( ip0_jp1_kp0_, ij)
            ipjm= myijk( ip1_jm1_kp0_, ij)
            ipjp= myijk( ip1_jp1_kp0_, ij)
            imjm= myijk( im1_jm1_kp0_, ij)
            imjp= myijk( im1_jp1_kp0_, ij)
            ippj= myijk( ip2_jp0_kp0_, ij)
            ijpp= myijk( ip0_jp2_kp0_, ij)
            immj= myijk( im2_jp0_kp0_, ij)
            ijmm= myijk( ip0_jm2_kp0_, ij)
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
!
!*****MX3D : ok (rivedere)
!
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
          sdim = MAX( SIZE(array,1) * snd_map(idest)%nsnd, 1)
          rdim = MAX( SIZE(array,1) * rcv_map(isour)%nrcv, 1)
          ALLOCATE( rcvbuf( rdim ) )
          ALLOCATE( sndbuf( sdim ) )
!          rcvbuf = 0.0
!          sndbuf = 0.0
          DO ib = 1, snd_map(idest)%nsnd
            DO ik = 1, SIZE(array,1)
              sndbuf(ik + SIZE(array,1)*(ib-1)) = array( ik, snd_map(idest)%iloc(ib) )
            END DO
          END DO 
          CALL sendrecv_real(sndbuf(1), SIZE(sndbuf), idest,     &      
                             rcvbuf(1), SIZE(rcvbuf), isour, ip)
          DO ib = 1, rcv_map(isour)%nrcv
            DO ik = 1, SIZE(array,1)
              array( ik, rcv_map(isour)%iloc(ib) ) = rcvbuf( ik + SIZE(array,1)*(ib-1))
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
