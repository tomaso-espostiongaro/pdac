!-----------------------------------------------------------------------
      MODULE dimensions
!-----------------------------------------------------------------------
        IMPLICIT NONE
        SAVE
        INTEGER :: nx, ny, nz
        INTEGER :: ntot ! ntot = nx*nz or nx*ny*nz
        INTEGER :: ntr  ! ntr = nx or nx*ny
        INTEGER :: no
        INTEGER :: nsolid
        INTEGER :: ngas
        INTEGER :: nphase
        INTEGER, PARAMETER :: nroughx    = 2   
        INTEGER, PARAMETER :: max_nsolid = 10  ! maximum number of solids
        INTEGER, PARAMETER :: max_size   = 512 ! maximum number of domain 
                                               ! subdivision for each direction
        INTEGER, PARAMETER :: max_ngas   = 7   ! max number of gas components
        INTEGER, PARAMETER :: max_nblock = 512 ! max number of blocks
        ! INTEGER, PARAMETER :: io_bufsiz = 2**18
!-----------------------------------------------------------------------
      END MODULE
!-----------------------------------------------------------------------
