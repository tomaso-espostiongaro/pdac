!----------------------------------------------------------------------
      MODULE io_restart
!----------------------------------------------------------------------
        USE eos_gas, ONLY: rgpgc_g, rgpgcn_g, xgc_g, ygc_g, cg_g
        USE gas_solid_density, ONLY: rog_g, rgp_g, rgpn_g, rlk_g, rlkn_g
        USE gas_solid_velocity, ONLY: ug_g, vg_g, uk_g, vk_g
        USE gas_solid_temperature, ONLY: sieg_g, siegn_g, tg_g,
     &       siek_g, siekn_g, tk_g
        USE parallel, ONLY: mpime, root
        USE pressure_epsilon, ONLY: p_g, pn_g, ep_g
        USE th_capacity, ONLY: cp_g, ck_g
        USE time_parameters, ONLY: time
        USE gas_solid_viscosity, ONLY: mug_g, kapg_g

        IMPLICIT NONE
        PRIVATE

        PUBLIC :: taperd, tapewr, tapebc

!----------------------------------------------------------------------
        CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE tapewr
!
      USE dimensions
      IMPLICIT NONE
!
      INTEGER :: i, j, k, ndindj, ij
      INTEGER :: kg
!
      CHARACTER*27 sysdump
      CHARACTER*22 sysdumpc
!
      IF( mpime .EQ. root ) THEN
!
! ... non standard fortran CALL (can be omitted)
!*************
          WRITE (sysdump,4445) char(0)
! ... IBM 
!          CALL system(sysdump)
 4445     FORMAT ('mv pdac.res.Z pdac.bak.Z',a1)
!*************
!
        OPEN(UNIT=9,form='unformatted',FILE='pdac.res')
!  
        WRITE(9) time,ndi,ndj,ncl
        ndindj=ndi*ndj
!
! ... store the final values of the main physical variables 
! 
        WRITE(9) (p_g(ij),pn_g(ij),ij=1,ndindj)
        WRITE(9) (ep_g(ij),ij=1,ndindj)
        WRITE(9) (rog_g(ij),ij=1,ndindj)
        WRITE(9) (rgp_g(ij),rgpn_g(ij),ij=1,ndindj)
        WRITE(9) ((rlk_g(k,ij), rlkn_g(k,ij),k=1,ncl),ij=1,ndindj) 
        WRITE(9) (ug_g(ij),vg_g(ij),ij=1,ndindj)
        WRITE(9) ((uk_g(k,ij),vk_g(k,ij),k=1,ncl),ij=1,ndindj)
        WRITE(9) (sieg_g(ij),siegn_g(ij),tg_g(ij),ij=1,ndindj)
        WRITE(9) ((siek_g(k,ij),siekn_g(k,ij),tk_g(k,ij),k=1,ncl),
     &             ij=1,ndindj)
        WRITE(9) ((ygc_g(kg,ij),xgc_g(kg,ij),kg=1,ngas),ij=1,ndindj)
       WRITE(9) ((rgpgc_g(kg,ij),rgpgcn_g(kg,ij),kg=1,ngas),ij=1,ndindj)
!
! ... store the final values of the constitutive parameters to be set up
!
        WRITE(9) (cg_g(ij),ij=1,ndindj)
        WRITE(9) ((ck_g(k,ij),k=1,ncl),ij=1,ndindj)
        WRITE(9) ((cp_g(kg,ij),kg=1,ngas),ij=1,ndindj)
        WRITE(9) (mug_g(ij),ij=1,ndindj)
        WRITE(9) (kapg_g(ij),ij=1,ndindj)
!
        CLOSE(9)
! 
! ... non standard fortran CALL (can be omitted)
!*************
        WRITE (sysdumpc,4446) char(0)
! ... IBM 
!        CALL system(sysdumpc)
 4446   FORMAT ('compress -f pdac.res',a1)
!*************
!
      END IF
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE taperd
!
      USE dimensions
      IMPLICIT NONE
!
      CHARACTER*23 sysdumpu
      INTEGER :: i, j, k, ndindj, ij
      INTEGER :: kg
      IF( mpime .EQ. root ) THEN
!
! non standard fortran CALL (can be omitted)
!*************
         WRITE (sysdumpu,4447) char(0)
 4447    FORMAT ('uncompress pdac.res.Z',a1)
! ... IBM 
!         CALL system(sysdumpu)
!*************
!
        OPEN(UNIT=9,form='unformatted',FILE='pdac.res')
!  
        READ(9) time,ndi,ndj,ncl
        ndindj=ndi*ndj
!
! ... read the final values of the main physical variables 
! 
        READ(9) (p_g(ij),pn_g(ij),ij=1,ndindj)
        READ(9) (ep_g(ij),ij=1,ndindj)
        READ(9) (rog_g(ij),ij=1,ndindj)
        READ(9) (rgp_g(ij),rgpn_g(ij),ij=1,ndindj)
        READ(9) ((rlk_g(k,ij), rlkn_g(k,ij),k=1,ncl),ij=1,ndindj) 
        READ(9) (ug_g(ij),vg_g(ij),ij=1,ndindj)
        READ(9) ((uk_g(k,ij),vk_g(k,ij),k=1,ncl),ij=1,ndindj)
        READ(9) (sieg_g(ij),siegn_g(ij),tg_g(ij),ij=1,ndindj)
        READ(9) ((siek_g(k,ij),siekn_g(k,ij),tk_g(k,ij),k=1,ncl),
     &            ij=1,ndindj)
        READ(9) ((ygc_g(kg,ij),xgc_g(kg,ij),kg=1,ngas),ij=1,ndindj)
        READ(9) ((rgpgc_g(kg,ij),rgpgcn_g(kg,ij),kg=1,ngas),ij=1,ndindj)
!
! ... read the final values of the constitutive parameters set up
!
        READ(9) (cg_g(ij),ij=1,ndindj)
        READ(9) ((ck_g(k,ij),k=1,ncl),ij=1,ndindj)
        READ(9) ((cp_g(kg,ij),kg=1,ngas),ij=1,ndindj)
        READ(9) (mug_g(ij),ij=1,ndindj)
        READ(9) (kapg_g(ij),ij=1,ndindj)
!
        CLOSE (9)
      END IF
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE tapebc
      USE parallel, ONLY: root
      IMPLICIT NONE
!
       CALL bcast_real(p_g,SIZE(p_g),root)
       CALL bcast_real(pn_g,SIZE(pn_g),root)
       CALL bcast_real(ep_g,SIZE(ep_g),root)
       CALL bcast_real(rog_g,SIZE(rog_g),root)
       CALL bcast_real(rgp_g,SIZE(rgp_g),root)
       CALL bcast_real(rgpn_g,SIZE(rgpn_g),root)
       CALL bcast_real(rlk_g,SIZE(rlk_g),root)
       CALL bcast_real(rlkn_g,SIZE(rlkn_g),root)
       CALL bcast_real(ug_g,SIZE(ug_g),root)
       CALL bcast_real(vg_g,SIZE(vg_g),root)
       CALL bcast_real(uk_g,SIZE(uk_g),root)
       CALL bcast_real(vk_g,SIZE(vk_g),root)
       CALL bcast_real(sieg_g,SIZE(sieg_g),root)
       CALL bcast_real(siegn_g,SIZE(siegn_g),root)
       CALL bcast_real(tg_g,SIZE(tg_g),root)
       CALL bcast_real(siek_g,SIZE(siek_g),root)
       CALL bcast_real(siekn_g,SIZE(siekn_g),root)
       CALL bcast_real(tk_g,SIZE(tk_g),root)
       CALL bcast_real(ygc_g,SIZE(ygc_g),root)
       CALL bcast_real(xgc_g,SIZE(xgc_g),root)
       CALL bcast_real(rgpgc_g,SIZE(rgpgc_g),root)
       CALL bcast_real(rgpgcn_g,SIZE(rgpgcn_g),root)
!
       CALL bcast_real(cg_g,SIZE(cg_g),root)
       CALL bcast_real(ck_g,SIZE(ck_g),root)
       CALL bcast_real(cp_g,SIZE(cp_g),root)
       CALL bcast_real(mug_g,SIZE(mug_g),root)
       CALL bcast_real(kapg_g,SIZE(kapg_g),root)
!
      END SUBROUTINE
!----------------------------------------------------------------------
      END MODULE io_restart
