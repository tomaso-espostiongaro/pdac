!----------------------------------------------------------------------
    PROGRAM restart_convert
!----------------------------------------------------------------------

      IMPLICIT NONE

      CALL bin2ascii
      ! CALL ascii2bin

    END PROGRAM


!----------------------------------------------------------------------


      SUBROUTINE bin2ascii
!
      IMPLICIT NONE
!
      INTEGER :: i, is
      INTEGER :: nx, ny, nz, nsolid, ngas, nfil
      INTEGER :: ig, gas_type( 7 )
      REAL*8, ALLOCATABLE :: tmp1d( : )
      REAL*8 :: time
      

      WRITE(*,*) ' reading from restart file '
!
      OPEN( UNIT=10, form='unformatted', status='old', FILE='pdac.res')
      OPEN( UNIT=11, form='formatted', FILE='pdac.ascii')

      READ(10) time, nx, ny, nz, nsolid, ngas, nfil
      READ(10) gas_type( 1 : ngas )

      WRITE(11,'(D20.14,6I5)') time, nx, ny, nz, nsolid, ngas, nfil
      WRITE(11,'(7I3)') gas_type( 1 : ngas )

      WRITE(*,*) ' time =  ', time
      WRITE(*,*) ' nx   =  ', nx
      WRITE(*,*) ' ny   =  ', ny
      WRITE(*,*) ' nz   =  ', nz
      WRITE(*,*) ' nsolid   =  ', nsolid
      WRITE(*,*) ' ngas     =  ', ngas
      WRITE(*,*) ' nfil     =  ', nfil
      WRITE(*,*) ' gas_type =  ', gas_type

      !
      ALLOCATE( tmp1d( nx*ny*nz ) )

      READ(10) tmp1d  ! P
      WRITE(11,'(10D22.15)') ( tmp1d( i ), i = 1, nx*ny*nz )

      READ(10) tmp1d  ! ug
      WRITE(11,'(10D22.15)') ( tmp1d( i ), i = 1, nx*ny*nz )

      READ(10) tmp1d  ! vg
      WRITE(11,'(10D22.15)') ( tmp1d( i ), i = 1, nx*ny*nz )

      READ(10) tmp1d  ! wg
      WRITE(11,'(10D22.15)') ( tmp1d( i ), i = 1, nx*ny*nz )

      READ(10) tmp1d  ! sieg
      WRITE(11,'(10D22.15)') ( tmp1d( i ), i = 1, nx*ny*nz )
      !
      DO is = 1, nsolid
        READ(10) tmp1d  ! sieg
        WRITE(11,'(10D22.15)') ( tmp1d( i ), i = 1, nx*ny*nz )
      END DO
      DO is = 1, nsolid
        READ(10) tmp1d  ! us
        WRITE(11,'(10D22.15)') ( tmp1d( i ), i = 1, nx*ny*nz )
        READ(10) tmp1d  ! vs
        WRITE(11,'(10D22.15)') ( tmp1d( i ), i = 1, nx*ny*nz )
        READ(10) tmp1d  ! ws
        WRITE(11,'(10D22.15)') ( tmp1d( i ), i = 1, nx*ny*nz )
      END DO
      DO is = 1, nsolid
        READ(10) tmp1d  ! sies
        WRITE(11,'(10D22.15)') ( tmp1d( i ), i = 1, nx*ny*nz )
      END DO
      DO ig = 1, ngas
        READ(10) tmp1d  ! ygc
        WRITE(11,'(10D22.15)') ( tmp1d( i ), i = 1, nx*ny*nz )
      END DO
      READ(10) tmp1d  ! rgp
      WRITE(11,'(10D22.15)') ( tmp1d( i ), i = 1, nx*ny*nz )
      READ(10) tmp1d  ! rog
      WRITE(11,'(10D22.15)') ( tmp1d( i ), i = 1, nx*ny*nz )
      READ(10) tmp1d  ! tg
      WRITE(11,'(10D22.15)') ( tmp1d( i ), i = 1, nx*ny*nz )
      DO is = 1, nsolid
        READ(10) tmp1d  ! ts
        WRITE(11,'(10D22.15)') ( tmp1d( i ), i = 1, nx*ny*nz )
      END DO
      CLOSE (10)
      CLOSE (11)

      DEALLOCATE( tmp1d )
!
      RETURN
      END SUBROUTINE bin2ascii


      SUBROUTINE ascii2bin
!
      IMPLICIT NONE
!
      INTEGER :: i, is
      INTEGER :: nx, ny, nz, nsolid, ngas, nfil
      INTEGER :: ig, gas_type( 100 )
      REAL*8, ALLOCATABLE :: tmp1d( : )
      REAL*8 :: time
      

      WRITE(*,*) ' writing restart file '
!
      OPEN( UNIT=10, form='unformatted', FILE='pdac.res')
      OPEN( UNIT=11, form='formatted', status='old', FILE='pdac.ascii')

      READ(11,*) time, nx, ny, nz, nsolid, ngas, nfil
      READ(11,*) gas_type( 1 : ngas )

      WRITE(10) time, nx, ny, nz, nsolid, ngas, nfil
      WRITE(10) gas_type( 1 : ngas )


      WRITE(*,*) ' time =  ', time
      WRITE(*,*) ' nx   =  ', nx
      WRITE(*,*) ' ny   =  ', ny
      WRITE(*,*) ' nz   =  ', nz
      WRITE(*,*) ' nsolid   =  ', nsolid
      WRITE(*,*) ' ngas     =  ', ngas
      WRITE(*,*) ' nfil     =  ', nfil
      WRITE(*,*) ' gas_type =  ', gas_type

      !
      ALLOCATE( tmp1d( nx*ny*nz ) )

      READ(11,'(10D22.15)') ( tmp1d( i ), i = 1, nx*ny*nz )
      WRITE(10) tmp1d  ! P

      READ(11,'(10D22.15)') ( tmp1d( i ), i = 1, nx*ny*nz )
      WRITE(10) tmp1d  ! ug

      READ(11,'(10D22.15)') ( tmp1d( i ), i = 1, nx*ny*nz )
      WRITE(10) tmp1d  ! vg

      READ(11,'(10D22.15)') ( tmp1d( i ), i = 1, nx*ny*nz )
      WRITE(10) tmp1d  ! wg

      READ(11,'(10D22.15)') ( tmp1d( i ), i = 1, nx*ny*nz )
      WRITE(10) tmp1d  ! sieg
      !
      DO is = 1, nsolid
        READ(11,'(10D22.15)') ( tmp1d( i ), i = 1, nx*ny*nz )
        WRITE(10) tmp1d  ! sieg
      END DO
      DO is = 1, nsolid
        READ(11,'(10D22.15)') ( tmp1d( i ), i = 1, nx*ny*nz )
        WRITE(10) tmp1d  ! us
        READ(11,'(10D22.15)') ( tmp1d( i ), i = 1, nx*ny*nz )
        WRITE(10) tmp1d  ! vs
        READ(11,'(10D22.15)') ( tmp1d( i ), i = 1, nx*ny*nz )
        WRITE(10) tmp1d  ! ws
      END DO
      DO is = 1, nsolid
        READ(11,'(10D22.15)') ( tmp1d( i ), i = 1, nx*ny*nz )
        WRITE(10) tmp1d  ! sies
      END DO
      DO ig = 1, ngas
        READ(11,'(10D22.15)') ( tmp1d( i ), i = 1, nx*ny*nz )
        WRITE(10) tmp1d  ! ygc
      END DO
      READ(11,'(10D22.15)') ( tmp1d( i ), i = 1, nx*ny*nz )
      WRITE(10) tmp1d  ! rgp
      READ(11,'(10D22.15)') ( tmp1d( i ), i = 1, nx*ny*nz )
      WRITE(10) tmp1d  ! rog
      READ(11,'(10D22.15)') ( tmp1d( i ), i = 1, nx*ny*nz )
      WRITE(10) tmp1d  ! tg
      DO is = 1, nsolid
        READ(11,'(10D22.15)') ( tmp1d( i ), i = 1, nx*ny*nz )
        WRITE(10) tmp1d  ! ts
      END DO
      CLOSE (10)
      CLOSE (11)

      DEALLOCATE( tmp1d )
!
      RETURN
      END SUBROUTINE ascii2bin
