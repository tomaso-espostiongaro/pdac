!----------------------------------------------------------------------
      MODULE output_dump
!----------------------------------------------------------------------
      USE dimensions
      IMPLICIT NONE
      INTEGER nfil
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------

      SUBROUTINE outp
!
      USE eos_gas, ONLY: xgc_g
      USE gas_solid_density, ONLY: rlk_g
      USE gas_solid_velocity, ONLY: ug_g, vg_g, uk_g, vk_g
      USE gas_solid_temperature, ONLY: tg_g, tk_g
      USE grid, ONLY: ib2, ib1, ib, jb2, jb1, jb
      USE parallel, ONLY: nproc, mpime, root, group
      USE particles_constants, ONLY: rl, inrl, nsolid
      USE pressure_epsilon, ONLY: p_g
      USE reactions, ONLY: irex
      USE time_parameters, ONLY: time
      USE turbulence, ONLY: scoeff_g
!
      IMPLICIT NONE
!
      CHARACTER :: filnam*11
      CHARACTER*4 :: letter
      LOGICAL :: cmprs
!
! this variable is used by the non standard routine system
!*************
      CHARACTER*24 syscom
!*************
!
      INTEGER :: ij2, ij1, k, j
      INTEGER :: ijl
      INTEGER :: kg
!
!      cmprs = .TRUE.
      cmprs = .FALSE.
!
      nfil=nfil+1
      filnam='OUTPUT.'//letter(nfil)

      IF( mpime .EQ. root ) THEN

      OPEN(UNIT=3,FILE=filnam)
!
      WRITE(3,547)time
      WRITE(3,548)
      DO 325 j=1,jb2
        ij1=1+(jb2-j)*ib2
        ij2=ib2+(jb2-j)*ib2
        WRITE(3,550)(p_g(ijl),ijl=ij1,ij2)
 325  CONTINUE
!
      DO 328 k=1,nsolid
        WRITE(3,549) k
      DO 328 j=1,jb2
        ij1=1+(jb2-j)*ib2
        ij2=ib2+(jb2-j)*ib2
        WRITE(3,550)(rlk_g(k,ijl)*inrl(k),ijl=ij1,ij2)
 328  CONTINUE
      WRITE(3,552)
      DO 332 j=1,jb2
        ij1=1+(jb2-j)*ib2
        ij2=ib2+(jb2-j)*ib2
        WRITE(3,550)(vg_g(ijl),ijl=ij1,ij2)
 332  CONTINUE
      WRITE(3,553)
      DO 333 j=1,jb2
        ij1=1+(jb2-j)*ib2
        ij2=ib2+(jb2-j)*ib2
        WRITE(3,550)(ug_g(ijl),ijl=ij1,ij2)
 333  CONTINUE
      IF(irex.GE.0) THEN
        WRITE(3,560)
        DO 340 j=1,jb2
          ij1=1+(jb2-j)*ib2
          ij2=ib2+(jb2-j)*ib2
          WRITE(3,550)(tg_g(ijl),ijl=ij1,ij2)
 340    CONTINUE
      ENDIF
      IF(irex.GE.0) THEN
!pdac------------
! define loop parameters
        DO 343 kg=5,5
!pdac------------
          WRITE(3,562) kg
          DO 342 j=1,jb2
            ij1=1+(jb2-j)*ib2
            ij2=ib2+(jb2-j)*ib2
            WRITE(3,550)(xgc_g(kg,ijl),ijl=ij1,ij2)
 342      CONTINUE
 343    CONTINUE
      ENDIF
      IF(nphase.EQ.1) RETURN
      DO 339 k=1,nsolid
        WRITE(3,556) k
        DO 336 j=1,jb2
          ij1=1+(jb2-j)*ib2 
          ij2=ib2+(jb2-j)*ib2
          WRITE(3,550)(vk_g(k,ijl),ijl=ij1,ij2)
 336    CONTINUE
        WRITE(3,557) k
        DO 337 j=1,jb2
          ij1=1+(jb2-j)*ib2
          ij2=ib2+(jb2-j)*ib2
          WRITE(3,550)(uk_g(k,ijl),ijl=ij1,ij2)
 337    CONTINUE
        IF(irex.GE.0) THEN
          WRITE(3,561) k
          DO 341 j=1,jb2
            ij1=1+(jb2-j)*ib2
            ij2=ib2+(jb2-j)*ib2
            WRITE(3,550)(tk_g(k,ijl),ijl=ij1,ij2)
 341      CONTINUE
        ENDIF
 339  CONTINUE
          WRITE(3,563) 
          DO 642 j=1,jb2
            ij1=1+(jb2-j)*ib2
            ij2=ib2+(jb2-j)*ib2
            WRITE(3,550)(scoeff_g(ijl),ijl=ij1,ij2)
 642      CONTINUE
      CLOSE (3)
      END IF
!
      IF (cmprs) THEN
      IF( mpime .EQ. root ) THEN
        WRITE (syscom,4444) nfil,char(0)
! ....   IBM 
!        CALL system(syscom)
 4444   FORMAT ('compress -f OUTPUT.',i4.4,a1)
      END IF
      END IF
!
      RETURN

 547  FORMAT(1x,///,1x,'@@@ TIME = ',G11.4)
 548  FORMAT(1x,//,1x,'P',/)
 549  FORMAT(1x,//,1x,'EPS',i1,/)
 550  FORMAT(1x,10(1x,g12.6))
 552  FORMAT(1x,//,1x,'VG',/)
 553  FORMAT(1x,//,1x,'UG',/)
 556  FORMAT(1x,//,1x,'VP',I1,/)
 557  FORMAT(1x,//,1x,'UP',I1,/)
 560  FORMAT(1x,//,1x,'TG',/)
 561  FORMAT(1x,//,1x,'TP',I1,/)
 562  FORMAT(1x,//,1x,'XGC',I1,/)
 563  FORMAT(1x,//,1x,'SMAGL',/)
!
      END SUBROUTINE

!----------------------------------------------------------------------
      END MODULE output_dump
!----------------------------------------------------------------------
