!-----------------------------------------------------------------------
      MODULE dimensions
!-----------------------------------------------------------------------
        IMPLICIT NONE
        SAVE
        INTEGER :: nx, ny, nz
        INTEGER :: ntot ! ntot = nx*nz or nx*ny*nz
        INTEGER :: no
        INTEGER :: nsolid
        INTEGER :: ngas
        INTEGER :: nphase
        INTEGER, PARAMETER :: nroughx    = 2   
        INTEGER, PARAMETER :: max_nsolid = 10  ! maximum number of solid phase 
        INTEGER, PARAMETER :: max_size   = 512 ! maximum number of domain subdivision 
                                               ! for each direction
        INTEGER, PARAMETER :: max_ngas   = 7   ! maximum number of gas components
        INTEGER, PARAMETER :: max_nblock = 512 ! maximum number of topology blocks
        ! INTEGER, PARAMETER :: io_bufsiz = 2**18
!-----------------------------------------------------------------------
      END MODULE
!-----------------------------------------------------------------------
