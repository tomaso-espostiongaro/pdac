!-----------------------------------------------------------------------
      MODULE dimensions
!-----------------------------------------------------------------------
! ... these variables are now completely equivalent to the input variables ...
! ... ib2, jb2, no, nsolid, ng ...
! ... this redundance will be soon eliminated ...
        INTEGER :: ndi
        INTEGER :: ndj
        INTEGER :: nnso
        INTEGER :: ncl
        INTEGER :: ngas
        INTEGER :: nphase
        INTEGER, PARAMETER :: nroughx = 2
        INTEGER, PARAMETER :: max_nsolid = 10
        INTEGER, PARAMETER :: max_size = 1024
        INTEGER, PARAMETER :: max_ngas = 10
        INTEGER, PARAMETER :: max_nblock = 512
!-----------------------------------------------------------------------
      END MODULE
