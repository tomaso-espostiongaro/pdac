# 1 "sendrecv.F"
      
! on T3e substitute MPI_DOUBLE_PRECISION with MPI_REAL
! on other architecture substitute MPI_REAL with MPI_DOUBLE_PRECISION
      
      
# 5
      SUBROUTINE sendrecv_real(sndbuf, snd_size, idest, 
     &  rcvbuf, rcv_size, isour, ip)
# 8
        IMPLICIT none
        REAL*8  :: sndbuf(*), rcvbuf(*)
! ... 
! "USMID @(#) mpi_t3e/MPI.inc	11.7	01/05/2000 18:49:31"
!     
!	(C) COPYRIGHT SILICON GRAPHICS, INC.
!	UNPUBLISHED PROPRIETARY INFORMATION.
!	ALL RIGHTS RESERVED.
!     
!     
! ------------------------------
! This is mpif.h for T3E systems
! ------------------------------
!     
! Copyright Notice
!     
! (c) Copyright 1995, 1996, 1997, The University of Edinburgh.
!     
! ======================================================================
!									
!     File: MPI.inc							
!  Project: KTG-T3DMPI							
!   Author: K.Cameron & A.G. Smith & K.J. Wierenga			
!  Created: 07/11/1994							
!  Descrip: The MPI Fortran interface header file			
! 									
! ======================================================================
      
! === Follows MPI (Message-passing Interface) Standard   Annex A. ===
      
! Version of MPI standard: VERSION.SUBVERSION 
# 10
        INTEGER :: snd_size, rcv_size, idest, isour, ip
# 1 "/opt/ctl/mpt/mpt/include/mpif.h" 1
# 29
       integer    MPI_VERSION, MPI_SUBVERSION
       parameter (MPI_VERSION    = 1)
      
! MPI-1 Return codes 
# 31
       parameter (MPI_SUBVERSION = 2)
# 34
       integer MPI_SUCCESS
      
# 35
       parameter (MPI_SUCCESS			= 0)
# 37
       integer MPI_ERR_BUFFER
       integer MPI_ERR_COUNT
       integer MPI_ERR_TYPE
       integer MPI_ERR_TAG
       integer MPI_ERR_COMM
       integer MPI_ERR_RANK
       integer MPI_ERR_REQUEST
       integer MPI_ERR_ROOT
       integer MPI_ERR_GROUP
       integer MPI_ERR_OP
       integer MPI_ERR_TOPOLOGY
       integer MPI_ERR_DIMS
       integer MPI_ERR_ARG
       integer MPI_ERR_UNKNOWN
       integer MPI_ERR_TRUNCATE
       integer MPI_ERR_OTHER
       integer MPI_ERR_INTERN
       integer MPI_ERR_IN_STATUS
      
# 55
       integer MPI_ERR_PENDING
# 57
       parameter (MPI_ERR_BUFFER                = 1)
       parameter (MPI_ERR_COUNT                 = 2)
       parameter (MPI_ERR_TYPE                  = 3)
       parameter (MPI_ERR_TAG                   = 4)
       parameter (MPI_ERR_COMM                  = 5)
       parameter (MPI_ERR_RANK                  = 6)
       parameter (MPI_ERR_REQUEST               = 7)
       parameter (MPI_ERR_ROOT                  = 8)
       parameter (MPI_ERR_GROUP                 = 9)
       parameter (MPI_ERR_OP                    = 10)
       parameter (MPI_ERR_TOPOLOGY              = 11)
       parameter (MPI_ERR_DIMS                  = 12)
       parameter (MPI_ERR_ARG                   = 13)
       parameter (MPI_ERR_UNKNOWN               = 14)
       parameter (MPI_ERR_TRUNCATE              = 15)
       parameter (MPI_ERR_OTHER                 = 16)
       parameter (MPI_ERR_INTERN                = 17)
       parameter (MPI_ERR_IN_STATUS             = 18)
      
! MPI-2 Return codes
      
# 75
       parameter (MPI_ERR_PENDING               = 19)
# 79
       integer MPI_ERR_ACCESS
       integer MPI_ERR_AMODE
       integer MPI_ERR_ASSERT
       integer MPI_ERR_BAD_FILE
       integer MPI_ERR_BASE
       integer MPI_ERR_CONVERSION
       integer MPI_ERR_DISP
       integer MPI_ERR_DUP_DATAREP
       integer MPI_ERR_FILE_EXISTS
       integer MPI_ERR_FILE_IN_USE
       integer MPI_ERR_FILE
       integer MPI_ERR_INFO_KEY
       integer MPI_ERR_INFO_NOKEY
       integer MPI_ERR_INFO_VALUE
       integer MPI_ERR_INFO
       integer MPI_ERR_IO
       integer MPI_ERR_KEYVAL
       integer MPI_ERR_LOCKTYPE
       integer MPI_ERR_NAME
       integer MPI_ERR_NO_MEM
       integer MPI_ERR_NOT_SAME
       integer MPI_ERR_NO_SPACE
       integer MPI_ERR_NO_SUCH_FILE
       integer MPI_ERR_PORT
       integer MPI_ERR_QUOTA
       integer MPI_ERR_READ_ONLY
       integer MPI_ERR_RMA_CONFLICT
       integer MPI_ERR_RMA_SYNC
       integer MPI_ERR_SERVICE
       integer MPI_ERR_SIZE
       integer MPI_ERR_SPAWN
       integer MPI_ERR_UNSUPPORTED_DATAREP
       integer MPI_ERR_UNSUPPORTED_OPERATION
      
# 112
       integer MPI_ERR_WIN
# 114
       parameter (MPI_ERR_ACCESS                  = 28)
       parameter (MPI_ERR_AMODE                   = 29)
       parameter (MPI_ERR_ASSERT                  = 30)
       parameter (MPI_ERR_BAD_FILE                = 31)
       parameter (MPI_ERR_BASE                    = 32)
       parameter (MPI_ERR_CONVERSION              = 33)
       parameter (MPI_ERR_DISP                    = 34)
       parameter (MPI_ERR_DUP_DATAREP             = 35)
       parameter (MPI_ERR_FILE_EXISTS             = 36)
       parameter (MPI_ERR_FILE_IN_USE             = 37)
       parameter (MPI_ERR_FILE                    = 38)
       parameter (MPI_ERR_INFO_KEY                = 39)
       parameter (MPI_ERR_INFO_NOKEY              = 40)
       parameter (MPI_ERR_INFO_VALUE              = 41)
       parameter (MPI_ERR_INFO                    = 42)
       parameter (MPI_ERR_IO                      = 43)
       parameter (MPI_ERR_KEYVAL                  = 44)
       parameter (MPI_ERR_LOCKTYPE                = 45)
       parameter (MPI_ERR_NAME                    = 46)
       parameter (MPI_ERR_NO_MEM                  = 47)
       parameter (MPI_ERR_NOT_SAME                = 48)
       parameter (MPI_ERR_NO_SPACE                = 49)
       parameter (MPI_ERR_NO_SUCH_FILE            = 50)
       parameter (MPI_ERR_PORT                    = 51)
       parameter (MPI_ERR_QUOTA                   = 52)
       parameter (MPI_ERR_READ_ONLY               = 53)
       parameter (MPI_ERR_RMA_CONFLICT            = 54)
       parameter (MPI_ERR_RMA_SYNC                = 55)
       parameter (MPI_ERR_SERVICE                 = 56)
       parameter (MPI_ERR_SIZE                    = 57)
       parameter (MPI_ERR_SPAWN                   = 58)
       parameter (MPI_ERR_UNSUPPORTED_DATAREP     = 59)
       parameter (MPI_ERR_UNSUPPORTED_OPERATION   = 60)
      
# 147
       parameter (MPI_ERR_WIN                     = 61)
# 149
       integer MPI_ERR_LASTCODE
      
! Assorted constants 
!     [ MPI_BOTTOM must be recognised by its address:
!       it must reside in a common block, value is not used ]
# 150
       parameter (MPI_ERR_LASTCODE                = 100)
# 155
      integer MPI_BOTTOM
      common /MPIPRIVBOT/ MPI_BOTTOM
      integer    MPI_PROC_NULL
      parameter (MPI_PROC_NULL       = -1)
      integer    MPI_ANY_SOURCE
      parameter (MPI_ANY_SOURCE      = -101)
      integer    MPI_ANY_TAG
      parameter (MPI_ANY_TAG         = -102)
      integer    MPI_UNDEFINED
!     MPI_UB -- as elementary datatype below
!     MPI_LB -- as elementary datatype below
# 164
      parameter (MPI_UNDEFINED       = -1)
# 167
      integer    MPI_BSEND_OVERHEAD
      parameter (MPI_BSEND_OVERHEAD  = 0)
      integer    MPI_KEYVAL_INVALID
      
! --- Status size and reserved index values (Fortran)
# 170
      parameter (MPI_KEYVAL_INVALID  = -1)
# 173
      integer    MPI_STATUS_SIZE
      parameter (MPI_STATUS_SIZE     = 6)
      integer    MPI_SOURCE
      parameter (MPI_SOURCE          = 1)
      integer    MPI_TAG
      parameter (MPI_TAG             = 2)
      integer    MPI_ERROR
      
! --- MPI-2 STATUS_IGNORE
!     [ MPI_STATUS(ES)_IGNORE must be recognised by their addresses:
!       they must reside in a common block, values are not used.
!     ]
# 180
      parameter (MPI_ERROR           = 3)
# 186
      integer    MPI_STATUS_IGNORE(MPI_STATUS_SIZE)
      integer    MPI_STATUSES_IGNORE(MPI_STATUS_SIZE, 1)
      
! --- Error-handling specifiers
# 188
      common /MPIPRIVSTAT/ MPI_STATUS_IGNORE, MPI_STATUSES_IGNORE
# 191
      integer    MPI_ERRORS_ARE_FATAL
      parameter (MPI_ERRORS_ARE_FATAL = 0)
      integer    MPI_ERRORS_RETURN
      
! --- Maximum sizes for strings
# 194
      parameter (MPI_ERRORS_RETURN    = 1)
# 197
      integer    MPI_MAX_PROCESSOR_NAME
      parameter (MPI_MAX_PROCESSOR_NAME = 256)
      integer    MPI_MAX_ERROR_STRING
      
! --- Elementary datatypes (Fortran)
# 200
      parameter (MPI_MAX_ERROR_STRING   = 256)
# 203
      integer    MPI_INTEGER
      parameter (MPI_INTEGER            = 14)
      integer    MPI_REAL
      parameter (MPI_REAL               = 15)
      integer    MPI_DOUBLE_PRECISION
      parameter (MPI_DOUBLE_PRECISION   = 16)
      integer    MPI_COMPLEX
! ---   [ nb: DOUBLE COMPLEX is optional basic datatype (v1.2) ]
# 210
      parameter (MPI_COMPLEX            = 17)
# 212
      integer    MPI_LOGICAL
      parameter (MPI_LOGICAL            = 19)
      integer    MPI_CHARACTER
      parameter (MPI_CHARACTER          = 20)
      integer    MPI_BYTE
      parameter (MPI_BYTE               = 21)
      integer    MPI_PACKED
      parameter (MPI_PACKED             = 22)
      integer    MPI_UB
      parameter (MPI_UB                 = 29)
      integer    MPI_LB
      
! --- Optional datatypes (Fortran)
# 223
      parameter (MPI_LB                 = 30)
# 226
      integer    MPI_DOUBLE_COMPLEX
      parameter (MPI_DOUBLE_COMPLEX     = 18)
      integer    MPI_INTEGER1
      parameter (MPI_INTEGER1           = 23)
      integer    MPI_INTEGER2
      parameter (MPI_INTEGER2           = 24)
      integer    MPI_INTEGER4
      parameter (MPI_INTEGER4           = 25)
      integer    MPI_INTEGER8
      parameter (MPI_INTEGER8           = 14)
      integer    MPI_REAL4
      parameter (MPI_REAL4              = 27)
      integer    MPI_REAL8
      
# 239
      parameter (MPI_REAL8              = 28)
# 241
      integer    MPI_LOGICAL1
      parameter (MPI_LOGICAL1           = 40)
      integer    MPI_LOGICAL2
      parameter (MPI_LOGICAL2           = 41)
      integer    MPI_LOGICAL4
      parameter (MPI_LOGICAL4           = 42)
      integer    MPI_COMPLEX8
      parameter (MPI_COMPLEX8           = 43)
      integer    MPI_COMPLEX16
      
! --- Datatypes for reduction functions (Fortran)
      
# 250
      parameter (MPI_COMPLEX16          = 44)
# 254
      integer    MPI_2REAL
      parameter (MPI_2REAL               = 37)
      integer    MPI_2DOUBLE_PRECISION
      parameter (MPI_2DOUBLE_PRECISION   = 38)
      integer    MPI_2INTEGER
!       MPI_2COMPLEX not defined [MPI Errata Oct 94]
      
! Reserved communicators 
# 259
      parameter (MPI_2INTEGER            = 39)
# 263
       integer    MPI_COMM_WORLD
       parameter (MPI_COMM_WORLD      = 0)
       integer    MPI_COMM_SELF
      
! Results of communicator and group comparisons 
# 266
       parameter (MPI_COMM_SELF       = 1)
# 269
       integer    MPI_UNEQUAL
       parameter (MPI_UNEQUAL         = -1)
       integer    MPI_SIMILAR
       parameter (MPI_SIMILAR         = -2)
       integer    MPI_IDENT
       parameter (MPI_IDENT           = -3)
       integer    MPI_CONGRUENT
      
! Environmental inquiry keys 
# 276
       parameter (MPI_CONGRUENT       = -4)
# 279
       integer    MPI_TAG_UB
       parameter (MPI_TAG_UB          = -50)
       integer    MPI_IO
       parameter (MPI_IO              = -49)
       integer    MPI_HOST
       parameter (MPI_HOST            = -48)
       integer    MPI_WTIME_IS_GLOBAL
      
! Collective operations 
# 286
       parameter (MPI_WTIME_IS_GLOBAL = -47)
# 289
       integer    MPI_MAX
       parameter (MPI_MAX             = 0)
       integer    MPI_MIN
       parameter (MPI_MIN             = 1)
       integer    MPI_SUM
       parameter (MPI_SUM             = 2)
       integer    MPI_PROD
       parameter (MPI_PROD            = 3)
       integer    MPI_MAXLOC
       parameter (MPI_MAXLOC          = 4)
       integer    MPI_MINLOC
       parameter (MPI_MINLOC          = 5)
       integer    MPI_BAND
       parameter (MPI_BAND            = 6)
       integer    MPI_BOR
       parameter (MPI_BOR             = 7)
       integer    MPI_BXOR
       parameter (MPI_BXOR            = 8)
       integer    MPI_LAND
       parameter (MPI_LAND            = 9)
       integer    MPI_LOR
       parameter (MPI_LOR             = 10)
       integer    MPI_LXOR
      
! MPI-2 collective operations
# 312
       parameter (MPI_LXOR            = 11)
# 315
       integer    MPI_REPLACE
      
! Null handles 
# 316
       parameter (MPI_REPLACE         = 12)
# 319
       integer    MPI_GROUP_NULL
       parameter (MPI_GROUP_NULL      = -1)
       integer    MPI_COMM_NULL
       parameter (MPI_COMM_NULL       = -1)
       integer    MPI_DATATYPE_NULL
       parameter (MPI_DATATYPE_NULL   = -1)
       integer    MPI_REQUEST_NULL
       parameter (MPI_REQUEST_NULL    = 0)
       integer    MPI_OP_NULL
       parameter (MPI_OP_NULL         = -1)
       integer    MPI_ERRHANDLER_NULL
      
! MPI-2 null handles
# 330
       parameter (MPI_ERRHANDLER_NULL = -1)
# 333
       integer    MPI_WIN_NULL
      
!     [ Predefined attribute callbacks must be recognised by address:
!       they must reside in a common block, values are not used ]
# 334
       parameter (MPI_WIN_NULL        = 0)
# 338
      integer    MPI_NULL_COPY_FN
      integer    MPI_DUP_FN
      integer    MPI_NULL_DELETE_FN
      common /MPIPRIVATTRCOPY/ MPI_NULL_COPY_FN,  MPI_DUP_FN
      
# 342
      common /MPIPRIVATTRDEL/  MPI_NULL_DELETE_FN
# 344
      integer    MPI_COMM_NULL_COPY_FN
      integer    MPI_COMM_DUP_FN
      integer    MPI_COMM_NULL_DELETE_FN
      common /MPIPRIVCOMATRCOPY/ MPI_COMM_NULL_COPY_FN, MPI_COMM_DUP_FN
      
! --- Empty group
# 348
      common /MPIPRIVCOMATRDEL/  MPI_COMM_NULL_DELETE_FN
# 351
      integer    MPI_GROUP_EMPTY
      
! --- Topologies
# 352
      parameter (MPI_GROUP_EMPTY = 0)
# 355
      integer    MPI_GRAPH
      parameter (MPI_GRAPH = 100)
      integer    MPI_CART
      
! --- MPI-2 Assertion Constants
# 358
      parameter (MPI_CART = 101)
# 361
      integer    MPI_MODE_NOCHECK
      parameter (MPI_MODE_NOCHECK      = 1)
      integer    MPI_MODE_NOSTORE
      parameter (MPI_MODE_NOSTORE      = 2)
      integer    MPI_MODE_NOPUT
      parameter (MPI_MODE_NOPUT        = 4)
      integer    MPI_MODE_NOPRECEDE
      parameter (MPI_MODE_NOPRECEDE    = 8)
      integer    MPI_MODE_NOSUCCEED
      
! --- MPI-2 Locktype Constants
# 370
      parameter (MPI_MODE_NOSUCCEED    = 16)
# 373
      integer    MPI_LOCK_EXCLUSIVE
      parameter (MPI_LOCK_EXCLUSIVE    = 1)
      integer    MPI_LOCK_SHARED
      
! --- Timing Functions
! --- (should be DOUBLE PRECISION -- not on Cray T3D/T3E)
# 376
      parameter (MPI_LOCK_SHARED       = 2)
# 380
      external MPI_WTICK
      real*8   MPI_WTICK
      external MPI_WTIME
      real*8   MPI_WTIME
      external PMPI_WTICK
      real*8   PMPI_WTICK
      external PMPI_WTIME
      
# 387
      real*8   PMPI_WTIME
# 389
      integer MPI_OFFSET_KIND
      
# 390
      parameter (MPI_OFFSET_KIND        = 8)
# 392
      integer MPI_ADDRESS_KIND
      
      
! info parameters
      
# 393
      parameter (MPI_ADDRESS_KIND       = 8)
# 398
      INTEGER(KIND=MPI_ADDRESS_KIND), PARAMETER:: MPI_INFO_NULL=0
      INTEGER, PARAMETER :: MPI_MAX_INFO_KEY=255
      
      
! MPI-2 I/O definitions
      
! USMID @(#) mpi_t3e/mpiof.h	11.1	11/12/98 08:39:46 
!     
!	(C) COPYRIGHT SILICON GRAPHICS, INC.
!	UNPUBLISHED PROPRIETARY INFORMATION.
!	ALL RIGHTS RESERVED.
!     
!     $Id: sendrecv.i,v 1.1 2002/09/16 09:08:01 k1vpizz1 Exp $    
!     
!     Copyright (C) 1997 University of Chicago. 
!     
!     
!    user include file for Fortran MPI-IO programs 
!     
# 400
      INTEGER, PARAMETER :: MPI_MAX_INFO_VAL=1024
# 1 "/opt/ctl/mpt/mpt/include/mpiof.h" 1
# 14
      INTEGER MPI_MODE_RDONLY, MPI_MODE_RDWR, MPI_MODE_WRONLY
      INTEGER MPI_MODE_DELETE_ON_CLOSE, MPI_MODE_UNIQUE_OPEN
      INTEGER MPI_MODE_CREATE, MPI_MODE_EXCL
      INTEGER MPI_MODE_APPEND
      PARAMETER (MPI_MODE_RDONLY=2, MPI_MODE_RDWR=8, MPI_MODE_WRONLY=4)
      PARAMETER (MPI_MODE_CREATE=1, MPI_MODE_DELETE_ON_CLOSE=16)
      PARAMETER (MPI_MODE_UNIQUE_OPEN=32, MPI_MODE_EXCL=64)
!     
# 21
      PARAMETER (MPI_MODE_APPEND=128)
# 23
      INTEGER MPI_FILE_NULL
!     
# 24
      PARAMETER (MPI_FILE_NULL=0)
# 26
      INTEGER MPI_MAX_DATAREP_STRING
!     
# 27
      PARAMETER (MPI_MAX_DATAREP_STRING=128)
# 29
      INTEGER MPI_SEEK_SET, MPI_SEEK_CUR, MPI_SEEK_END
!     
# 30
      PARAMETER (MPI_SEEK_SET=600, MPI_SEEK_CUR=602, MPI_SEEK_END=604)
# 32
      INTEGER MPIO_REQUEST_NULL
!     
!      INTEGER MPI_OFFSET_KIND
!      PARAMETER (MPI_OFFSET_KIND=8)
!     
# 33
      PARAMETER (MPIO_REQUEST_NULL=0)
# 38
      INTEGER MPI_ORDER_C, MPI_ORDER_FORTRAN
      PARAMETER (MPI_ORDER_C=56, MPI_ORDER_FORTRAN=57)
      INTEGER MPI_DISTRIBUTE_BLOCK, MPI_DISTRIBUTE_CYCLIC
      INTEGER MPI_DISTRIBUTE_NONE, MPI_DISTRIBUTE_DFLT_DARG
      PARAMETER (MPI_DISTRIBUTE_BLOCK=121, MPI_DISTRIBUTE_CYCLIC=122)
      PARAMETER (MPI_DISTRIBUTE_NONE=123)
!     
!     
!     
!     
!     
# 44
      PARAMETER (MPI_DISTRIBUTE_DFLT_DARG=-49767)
# 14 "sendrecv.F" 2
# 14 "sendrecv.F"
# 14
        INTEGER ISTATUS(MPI_STATUS_SIZE), ierr
        CALL MPI_SENDRECV(sndbuf, snd_size, MPI_DOUBLE_PRECISION, 
     &    IDEST, ip, rcvbuf, rcv_size, MPI_DOUBLE_PRECISION, 
     &    ISOUR, ip, MPI_COMM_WORLD, ISTATUS, ierr)
        IF(ierr .NE. 0) THEN
          WRITE(6,*) ' ** ERROR in sendrecv ** '
          STOP
        END IF 
# 24
        RETURN
      
      
      
# 25
      END SUBROUTINE
      
# 29
      SUBROUTINE sendrecv_integer(sndbuf, snd_size, idest,
     &  rcvbuf, rcv_size, isour, ip)
# 32
        IMPLICIT none
        INTEGER  :: sndbuf(*), rcvbuf(*)
! ... 
! "USMID @(#) mpi_t3e/MPI.inc	11.7	01/05/2000 18:49:31"
!     
!	(C) COPYRIGHT SILICON GRAPHICS, INC.
!	UNPUBLISHED PROPRIETARY INFORMATION.
!	ALL RIGHTS RESERVED.
!     
!     
! ------------------------------
! This is mpif.h for T3E systems
! ------------------------------
!     
! Copyright Notice
!     
! (c) Copyright 1995, 1996, 1997, The University of Edinburgh.
!     
! ======================================================================
!									
!     File: MPI.inc							
!  Project: KTG-T3DMPI							
!   Author: K.Cameron & A.G. Smith & K.J. Wierenga			
!  Created: 07/11/1994							
!  Descrip: The MPI Fortran interface header file			
! 									
! ======================================================================
      
! === Follows MPI (Message-passing Interface) Standard   Annex A. ===
      
! Version of MPI standard: VERSION.SUBVERSION 
# 34
        INTEGER :: snd_size, rcv_size, idest, isour, ip
# 1 "/opt/ctl/mpt/mpt/include/mpif.h" 1
# 29
       integer    MPI_VERSION, MPI_SUBVERSION
       parameter (MPI_VERSION    = 1)
      
! MPI-1 Return codes 
# 31
       parameter (MPI_SUBVERSION = 2)
# 34
       integer MPI_SUCCESS
      
# 35
       parameter (MPI_SUCCESS			= 0)
# 37
       integer MPI_ERR_BUFFER
       integer MPI_ERR_COUNT
       integer MPI_ERR_TYPE
       integer MPI_ERR_TAG
       integer MPI_ERR_COMM
       integer MPI_ERR_RANK
       integer MPI_ERR_REQUEST
       integer MPI_ERR_ROOT
       integer MPI_ERR_GROUP
       integer MPI_ERR_OP
       integer MPI_ERR_TOPOLOGY
       integer MPI_ERR_DIMS
       integer MPI_ERR_ARG
       integer MPI_ERR_UNKNOWN
       integer MPI_ERR_TRUNCATE
       integer MPI_ERR_OTHER
       integer MPI_ERR_INTERN
       integer MPI_ERR_IN_STATUS
      
# 55
       integer MPI_ERR_PENDING
# 57
       parameter (MPI_ERR_BUFFER                = 1)
       parameter (MPI_ERR_COUNT                 = 2)
       parameter (MPI_ERR_TYPE                  = 3)
       parameter (MPI_ERR_TAG                   = 4)
       parameter (MPI_ERR_COMM                  = 5)
       parameter (MPI_ERR_RANK                  = 6)
       parameter (MPI_ERR_REQUEST               = 7)
       parameter (MPI_ERR_ROOT                  = 8)
       parameter (MPI_ERR_GROUP                 = 9)
       parameter (MPI_ERR_OP                    = 10)
       parameter (MPI_ERR_TOPOLOGY              = 11)
       parameter (MPI_ERR_DIMS                  = 12)
       parameter (MPI_ERR_ARG                   = 13)
       parameter (MPI_ERR_UNKNOWN               = 14)
       parameter (MPI_ERR_TRUNCATE              = 15)
       parameter (MPI_ERR_OTHER                 = 16)
       parameter (MPI_ERR_INTERN                = 17)
       parameter (MPI_ERR_IN_STATUS             = 18)
      
! MPI-2 Return codes
      
# 75
       parameter (MPI_ERR_PENDING               = 19)
# 79
       integer MPI_ERR_ACCESS
       integer MPI_ERR_AMODE
       integer MPI_ERR_ASSERT
       integer MPI_ERR_BAD_FILE
       integer MPI_ERR_BASE
       integer MPI_ERR_CONVERSION
       integer MPI_ERR_DISP
       integer MPI_ERR_DUP_DATAREP
       integer MPI_ERR_FILE_EXISTS
       integer MPI_ERR_FILE_IN_USE
       integer MPI_ERR_FILE
       integer MPI_ERR_INFO_KEY
       integer MPI_ERR_INFO_NOKEY
       integer MPI_ERR_INFO_VALUE
       integer MPI_ERR_INFO
       integer MPI_ERR_IO
       integer MPI_ERR_KEYVAL
       integer MPI_ERR_LOCKTYPE
       integer MPI_ERR_NAME
       integer MPI_ERR_NO_MEM
       integer MPI_ERR_NOT_SAME
       integer MPI_ERR_NO_SPACE
       integer MPI_ERR_NO_SUCH_FILE
       integer MPI_ERR_PORT
       integer MPI_ERR_QUOTA
       integer MPI_ERR_READ_ONLY
       integer MPI_ERR_RMA_CONFLICT
       integer MPI_ERR_RMA_SYNC
       integer MPI_ERR_SERVICE
       integer MPI_ERR_SIZE
       integer MPI_ERR_SPAWN
       integer MPI_ERR_UNSUPPORTED_DATAREP
       integer MPI_ERR_UNSUPPORTED_OPERATION
      
# 112
       integer MPI_ERR_WIN
# 114
       parameter (MPI_ERR_ACCESS                  = 28)
       parameter (MPI_ERR_AMODE                   = 29)
       parameter (MPI_ERR_ASSERT                  = 30)
       parameter (MPI_ERR_BAD_FILE                = 31)
       parameter (MPI_ERR_BASE                    = 32)
       parameter (MPI_ERR_CONVERSION              = 33)
       parameter (MPI_ERR_DISP                    = 34)
       parameter (MPI_ERR_DUP_DATAREP             = 35)
       parameter (MPI_ERR_FILE_EXISTS             = 36)
       parameter (MPI_ERR_FILE_IN_USE             = 37)
       parameter (MPI_ERR_FILE                    = 38)
       parameter (MPI_ERR_INFO_KEY                = 39)
       parameter (MPI_ERR_INFO_NOKEY              = 40)
       parameter (MPI_ERR_INFO_VALUE              = 41)
       parameter (MPI_ERR_INFO                    = 42)
       parameter (MPI_ERR_IO                      = 43)
       parameter (MPI_ERR_KEYVAL                  = 44)
       parameter (MPI_ERR_LOCKTYPE                = 45)
       parameter (MPI_ERR_NAME                    = 46)
       parameter (MPI_ERR_NO_MEM                  = 47)
       parameter (MPI_ERR_NOT_SAME                = 48)
       parameter (MPI_ERR_NO_SPACE                = 49)
       parameter (MPI_ERR_NO_SUCH_FILE            = 50)
       parameter (MPI_ERR_PORT                    = 51)
       parameter (MPI_ERR_QUOTA                   = 52)
       parameter (MPI_ERR_READ_ONLY               = 53)
       parameter (MPI_ERR_RMA_CONFLICT            = 54)
       parameter (MPI_ERR_RMA_SYNC                = 55)
       parameter (MPI_ERR_SERVICE                 = 56)
       parameter (MPI_ERR_SIZE                    = 57)
       parameter (MPI_ERR_SPAWN                   = 58)
       parameter (MPI_ERR_UNSUPPORTED_DATAREP     = 59)
       parameter (MPI_ERR_UNSUPPORTED_OPERATION   = 60)
      
# 147
       parameter (MPI_ERR_WIN                     = 61)
# 149
       integer MPI_ERR_LASTCODE
      
! Assorted constants 
!     [ MPI_BOTTOM must be recognised by its address:
!       it must reside in a common block, value is not used ]
# 150
       parameter (MPI_ERR_LASTCODE                = 100)
# 155
      integer MPI_BOTTOM
      common /MPIPRIVBOT/ MPI_BOTTOM
      integer    MPI_PROC_NULL
      parameter (MPI_PROC_NULL       = -1)
      integer    MPI_ANY_SOURCE
      parameter (MPI_ANY_SOURCE      = -101)
      integer    MPI_ANY_TAG
      parameter (MPI_ANY_TAG         = -102)
      integer    MPI_UNDEFINED
!     MPI_UB -- as elementary datatype below
!     MPI_LB -- as elementary datatype below
# 164
      parameter (MPI_UNDEFINED       = -1)
# 167
      integer    MPI_BSEND_OVERHEAD
      parameter (MPI_BSEND_OVERHEAD  = 0)
      integer    MPI_KEYVAL_INVALID
      
! --- Status size and reserved index values (Fortran)
# 170
      parameter (MPI_KEYVAL_INVALID  = -1)
# 173
      integer    MPI_STATUS_SIZE
      parameter (MPI_STATUS_SIZE     = 6)
      integer    MPI_SOURCE
      parameter (MPI_SOURCE          = 1)
      integer    MPI_TAG
      parameter (MPI_TAG             = 2)
      integer    MPI_ERROR
      
! --- MPI-2 STATUS_IGNORE
!     [ MPI_STATUS(ES)_IGNORE must be recognised by their addresses:
!       they must reside in a common block, values are not used.
!     ]
# 180
      parameter (MPI_ERROR           = 3)
# 186
      integer    MPI_STATUS_IGNORE(MPI_STATUS_SIZE)
      integer    MPI_STATUSES_IGNORE(MPI_STATUS_SIZE, 1)
      
! --- Error-handling specifiers
# 188
      common /MPIPRIVSTAT/ MPI_STATUS_IGNORE, MPI_STATUSES_IGNORE
# 191
      integer    MPI_ERRORS_ARE_FATAL
      parameter (MPI_ERRORS_ARE_FATAL = 0)
      integer    MPI_ERRORS_RETURN
      
! --- Maximum sizes for strings
# 194
      parameter (MPI_ERRORS_RETURN    = 1)
# 197
      integer    MPI_MAX_PROCESSOR_NAME
      parameter (MPI_MAX_PROCESSOR_NAME = 256)
      integer    MPI_MAX_ERROR_STRING
      
! --- Elementary datatypes (Fortran)
# 200
      parameter (MPI_MAX_ERROR_STRING   = 256)
# 203
      integer    MPI_INTEGER
      parameter (MPI_INTEGER            = 14)
      integer    MPI_REAL
      parameter (MPI_REAL               = 15)
      integer    MPI_DOUBLE_PRECISION
      parameter (MPI_DOUBLE_PRECISION   = 16)
      integer    MPI_COMPLEX
! ---   [ nb: DOUBLE COMPLEX is optional basic datatype (v1.2) ]
# 210
      parameter (MPI_COMPLEX            = 17)
# 212
      integer    MPI_LOGICAL
      parameter (MPI_LOGICAL            = 19)
      integer    MPI_CHARACTER
      parameter (MPI_CHARACTER          = 20)
      integer    MPI_BYTE
      parameter (MPI_BYTE               = 21)
      integer    MPI_PACKED
      parameter (MPI_PACKED             = 22)
      integer    MPI_UB
      parameter (MPI_UB                 = 29)
      integer    MPI_LB
      
! --- Optional datatypes (Fortran)
# 223
      parameter (MPI_LB                 = 30)
# 226
      integer    MPI_DOUBLE_COMPLEX
      parameter (MPI_DOUBLE_COMPLEX     = 18)
      integer    MPI_INTEGER1
      parameter (MPI_INTEGER1           = 23)
      integer    MPI_INTEGER2
      parameter (MPI_INTEGER2           = 24)
      integer    MPI_INTEGER4
      parameter (MPI_INTEGER4           = 25)
      integer    MPI_INTEGER8
      parameter (MPI_INTEGER8           = 14)
      integer    MPI_REAL4
      parameter (MPI_REAL4              = 27)
      integer    MPI_REAL8
      
# 239
      parameter (MPI_REAL8              = 28)
# 241
      integer    MPI_LOGICAL1
      parameter (MPI_LOGICAL1           = 40)
      integer    MPI_LOGICAL2
      parameter (MPI_LOGICAL2           = 41)
      integer    MPI_LOGICAL4
      parameter (MPI_LOGICAL4           = 42)
      integer    MPI_COMPLEX8
      parameter (MPI_COMPLEX8           = 43)
      integer    MPI_COMPLEX16
      
! --- Datatypes for reduction functions (Fortran)
      
# 250
      parameter (MPI_COMPLEX16          = 44)
# 254
      integer    MPI_2REAL
      parameter (MPI_2REAL               = 37)
      integer    MPI_2DOUBLE_PRECISION
      parameter (MPI_2DOUBLE_PRECISION   = 38)
      integer    MPI_2INTEGER
!       MPI_2COMPLEX not defined [MPI Errata Oct 94]
      
! Reserved communicators 
# 259
      parameter (MPI_2INTEGER            = 39)
# 263
       integer    MPI_COMM_WORLD
       parameter (MPI_COMM_WORLD      = 0)
       integer    MPI_COMM_SELF
      
! Results of communicator and group comparisons 
# 266
       parameter (MPI_COMM_SELF       = 1)
# 269
       integer    MPI_UNEQUAL
       parameter (MPI_UNEQUAL         = -1)
       integer    MPI_SIMILAR
       parameter (MPI_SIMILAR         = -2)
       integer    MPI_IDENT
       parameter (MPI_IDENT           = -3)
       integer    MPI_CONGRUENT
      
! Environmental inquiry keys 
# 276
       parameter (MPI_CONGRUENT       = -4)
# 279
       integer    MPI_TAG_UB
       parameter (MPI_TAG_UB          = -50)
       integer    MPI_IO
       parameter (MPI_IO              = -49)
       integer    MPI_HOST
       parameter (MPI_HOST            = -48)
       integer    MPI_WTIME_IS_GLOBAL
      
! Collective operations 
# 286
       parameter (MPI_WTIME_IS_GLOBAL = -47)
# 289
       integer    MPI_MAX
       parameter (MPI_MAX             = 0)
       integer    MPI_MIN
       parameter (MPI_MIN             = 1)
       integer    MPI_SUM
       parameter (MPI_SUM             = 2)
       integer    MPI_PROD
       parameter (MPI_PROD            = 3)
       integer    MPI_MAXLOC
       parameter (MPI_MAXLOC          = 4)
       integer    MPI_MINLOC
       parameter (MPI_MINLOC          = 5)
       integer    MPI_BAND
       parameter (MPI_BAND            = 6)
       integer    MPI_BOR
       parameter (MPI_BOR             = 7)
       integer    MPI_BXOR
       parameter (MPI_BXOR            = 8)
       integer    MPI_LAND
       parameter (MPI_LAND            = 9)
       integer    MPI_LOR
       parameter (MPI_LOR             = 10)
       integer    MPI_LXOR
      
! MPI-2 collective operations
# 312
       parameter (MPI_LXOR            = 11)
# 315
       integer    MPI_REPLACE
      
! Null handles 
# 316
       parameter (MPI_REPLACE         = 12)
# 319
       integer    MPI_GROUP_NULL
       parameter (MPI_GROUP_NULL      = -1)
       integer    MPI_COMM_NULL
       parameter (MPI_COMM_NULL       = -1)
       integer    MPI_DATATYPE_NULL
       parameter (MPI_DATATYPE_NULL   = -1)
       integer    MPI_REQUEST_NULL
       parameter (MPI_REQUEST_NULL    = 0)
       integer    MPI_OP_NULL
       parameter (MPI_OP_NULL         = -1)
       integer    MPI_ERRHANDLER_NULL
      
! MPI-2 null handles
# 330
       parameter (MPI_ERRHANDLER_NULL = -1)
# 333
       integer    MPI_WIN_NULL
      
!     [ Predefined attribute callbacks must be recognised by address:
!       they must reside in a common block, values are not used ]
# 334
       parameter (MPI_WIN_NULL        = 0)
# 338
      integer    MPI_NULL_COPY_FN
      integer    MPI_DUP_FN
      integer    MPI_NULL_DELETE_FN
      common /MPIPRIVATTRCOPY/ MPI_NULL_COPY_FN,  MPI_DUP_FN
      
# 342
      common /MPIPRIVATTRDEL/  MPI_NULL_DELETE_FN
# 344
      integer    MPI_COMM_NULL_COPY_FN
      integer    MPI_COMM_DUP_FN
      integer    MPI_COMM_NULL_DELETE_FN
      common /MPIPRIVCOMATRCOPY/ MPI_COMM_NULL_COPY_FN, MPI_COMM_DUP_FN
      
! --- Empty group
# 348
      common /MPIPRIVCOMATRDEL/  MPI_COMM_NULL_DELETE_FN
# 351
      integer    MPI_GROUP_EMPTY
      
! --- Topologies
# 352
      parameter (MPI_GROUP_EMPTY = 0)
# 355
      integer    MPI_GRAPH
      parameter (MPI_GRAPH = 100)
      integer    MPI_CART
      
! --- MPI-2 Assertion Constants
# 358
      parameter (MPI_CART = 101)
# 361
      integer    MPI_MODE_NOCHECK
      parameter (MPI_MODE_NOCHECK      = 1)
      integer    MPI_MODE_NOSTORE
      parameter (MPI_MODE_NOSTORE      = 2)
      integer    MPI_MODE_NOPUT
      parameter (MPI_MODE_NOPUT        = 4)
      integer    MPI_MODE_NOPRECEDE
      parameter (MPI_MODE_NOPRECEDE    = 8)
      integer    MPI_MODE_NOSUCCEED
      
! --- MPI-2 Locktype Constants
# 370
      parameter (MPI_MODE_NOSUCCEED    = 16)
# 373
      integer    MPI_LOCK_EXCLUSIVE
      parameter (MPI_LOCK_EXCLUSIVE    = 1)
      integer    MPI_LOCK_SHARED
      
! --- Timing Functions
! --- (should be DOUBLE PRECISION -- not on Cray T3D/T3E)
# 376
      parameter (MPI_LOCK_SHARED       = 2)
# 380
      external MPI_WTICK
      real*8   MPI_WTICK
      external MPI_WTIME
      real*8   MPI_WTIME
      external PMPI_WTICK
      real*8   PMPI_WTICK
      external PMPI_WTIME
      
# 387
      real*8   PMPI_WTIME
# 389
      integer MPI_OFFSET_KIND
      
# 390
      parameter (MPI_OFFSET_KIND        = 8)
# 392
      integer MPI_ADDRESS_KIND
      
      
! info parameters
      
# 393
      parameter (MPI_ADDRESS_KIND       = 8)
# 398
      INTEGER(KIND=MPI_ADDRESS_KIND), PARAMETER:: MPI_INFO_NULL=0
      INTEGER, PARAMETER :: MPI_MAX_INFO_KEY=255
      
      
! MPI-2 I/O definitions
      
! USMID @(#) mpi_t3e/mpiof.h	11.1	11/12/98 08:39:46 
!     
!	(C) COPYRIGHT SILICON GRAPHICS, INC.
!	UNPUBLISHED PROPRIETARY INFORMATION.
!	ALL RIGHTS RESERVED.
!     
!     $Id: sendrecv.i,v 1.1 2002/09/16 09:08:01 k1vpizz1 Exp $    
!     
!     Copyright (C) 1997 University of Chicago. 
!     
!     
!    user include file for Fortran MPI-IO programs 
!     
# 400
      INTEGER, PARAMETER :: MPI_MAX_INFO_VAL=1024
# 1 "/opt/ctl/mpt/mpt/include/mpiof.h" 1
# 14
      INTEGER MPI_MODE_RDONLY, MPI_MODE_RDWR, MPI_MODE_WRONLY
      INTEGER MPI_MODE_DELETE_ON_CLOSE, MPI_MODE_UNIQUE_OPEN
      INTEGER MPI_MODE_CREATE, MPI_MODE_EXCL
      INTEGER MPI_MODE_APPEND
      PARAMETER (MPI_MODE_RDONLY=2, MPI_MODE_RDWR=8, MPI_MODE_WRONLY=4)
      PARAMETER (MPI_MODE_CREATE=1, MPI_MODE_DELETE_ON_CLOSE=16)
      PARAMETER (MPI_MODE_UNIQUE_OPEN=32, MPI_MODE_EXCL=64)
!     
# 21
      PARAMETER (MPI_MODE_APPEND=128)
# 23
      INTEGER MPI_FILE_NULL
!     
# 24
      PARAMETER (MPI_FILE_NULL=0)
# 26
      INTEGER MPI_MAX_DATAREP_STRING
!     
# 27
      PARAMETER (MPI_MAX_DATAREP_STRING=128)
# 29
      INTEGER MPI_SEEK_SET, MPI_SEEK_CUR, MPI_SEEK_END
!     
# 30
      PARAMETER (MPI_SEEK_SET=600, MPI_SEEK_CUR=602, MPI_SEEK_END=604)
# 32
      INTEGER MPIO_REQUEST_NULL
!     
!      INTEGER MPI_OFFSET_KIND
!      PARAMETER (MPI_OFFSET_KIND=8)
!     
# 33
      PARAMETER (MPIO_REQUEST_NULL=0)
# 38
      INTEGER MPI_ORDER_C, MPI_ORDER_FORTRAN
      PARAMETER (MPI_ORDER_C=56, MPI_ORDER_FORTRAN=57)
      INTEGER MPI_DISTRIBUTE_BLOCK, MPI_DISTRIBUTE_CYCLIC
      INTEGER MPI_DISTRIBUTE_NONE, MPI_DISTRIBUTE_DFLT_DARG
      PARAMETER (MPI_DISTRIBUTE_BLOCK=121, MPI_DISTRIBUTE_CYCLIC=122)
      PARAMETER (MPI_DISTRIBUTE_NONE=123)
!     
!     
!     
!     
!     
# 44
      PARAMETER (MPI_DISTRIBUTE_DFLT_DARG=-49767)
# 38 "sendrecv.F" 2
# 38 "sendrecv.F"
# 38
        INTEGER ISTATUS(MPI_STATUS_SIZE), ierr
        CALL MPI_SENDRECV(sndbuf, snd_size, MPI_INTEGER,
     &    IDEST, ip, rcvbuf, rcv_size, MPI_INTEGER,
     &    ISOUR, ip, MPI_COMM_WORLD, ISTATUS, ierr)
        IF(ierr .NE. 0) THEN
          WRITE(6,*) ' ** ERROR in sendrecv ** '
          STOP
        END IF
# 48
        RETURN
      
# 49
      END SUBROUTINE
