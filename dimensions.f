!-----------------------------------------------------------------------
      MODULE dimensions
!-----------------------------------------------------------------------
        INTEGER :: nx, ny, nz
        INTEGER :: ntot ! ntot = nx*nz or nx*ny*nz
        INTEGER :: no
        INTEGER :: nsolid
        INTEGER :: ngas
        INTEGER :: nphase
        INTEGER, PARAMETER :: nroughx = 2
        INTEGER, PARAMETER :: max_nsolid = 10
        INTEGER, PARAMETER :: max_size = 1024
        INTEGER, PARAMETER :: max_ngas = 10
        INTEGER, PARAMETER :: max_nblock = 512
        ! INTEGER, PARAMETER :: io_bufsiz = 2**18
!-----------------------------------------------------------------------
      END MODULE
!-----------------------------------------------------------------------
