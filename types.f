      MODULE basic_types
        TYPE iarray
          INTEGER, POINTER :: i(:)
        END TYPE
        TYPE imatrix
          INTEGER, POINTER :: i(:,:)
        END TYPE
      END MODULE basic_types
