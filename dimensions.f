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
!-----------------------------------------------------------------------
      END MODULE
