!----------------------------------------------------------------------
      MODULE eos_gas
!----------------------------------------------------------------------
      IMPLICIT NONE
      SAVE
!
      REAL*8, DIMENSION(:), ALLOCATABLE :: cg_g
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: ygc_g, xgc_g
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: rgpgc_g, rgpgcn_g
      REAL*8, DIMENSION(:), ALLOCATABLE :: cg
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: ygc, xgc
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: rgpgc, rgpgcn
      REAL*8 :: rags
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE bounds_eosg
      USE dimensions
      IMPLICIT NONE
!
       ALLOCATE(cg_g(ndi*ndj))
       ALLOCATE(ygc_g(ngas,ndi*ndj),xgc_g(ngas,ndi*ndj),      &
                rgpgc_g(ngas,ndi*ndj), rgpgcn_g(ngas,ndi*ndj))

       cg_g     = 0.0d0
       ygc_g    = 0.0d0
       xgc_g    = 0.0d0
       rgpgc_g  = 0.0d0
       rgpgcn_g = 0.0d0

      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE local_bounds_eosg
      USE dimensions
      USE grid, ONLY: nijx_l, nij_l
      IMPLICIT NONE
!
       ALLOCATE(cg(nij_l))
       ALLOCATE(ygc(ngas,nij_l),xgc(ngas,nij_l),             &
                rgpgc(ngas,nijx_l), rgpgcn(ngas,nij_l))
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE mole(xgc, ygc)
!
! ... computes molar fractions
!
      USE dimensions
      USE gas_constants, ONLY: gmw
      IMPLICIT NONE
!
      REAL*8, INTENT(IN) :: ygc(:)
      REAL*8, INTENT(OUT) :: xgc(:)
      INTEGER :: kg
      REAL*8 :: c1
!
      c1=0.D0
      DO kg=1,ngas
        c1=c1+ygc(kg)/gmw(kg)
      END DO
      DO kg=1,ngas
        xgc(kg)=ygc(kg)/gmw(kg)/c1
      END DO
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE eosg(rags, rog, cp, cg, tg, ygc, xgc, sieg, p, nt, nr, nc, ij)
!
! ... updates gas density with new pressure
!
      USE dimensions
      USE gas_constants, ONLY: gmw, c_erg, rgas, tzero, hzerog
      USE th_capacity, ONLY:  hcapg
      USE time_parameters, ONLY: time
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: ij
      REAL*8, INTENT(IN) :: ygc(:), xgc(:), sieg, p
      REAL*8, INTENT(OUT) :: rog, rags
      REAL*8, INTENT(INOUT) :: tg
      REAL*8 :: cp(:), cg
      INTEGER :: nr, nt, nc
!
      INTEGER :: i,j
      REAL*8 :: tgnn, c1, c2, ratmin
      INTEGER :: ii, nlmax
      INTEGER :: kg
      PARAMETER( nlmax = 2000) 
      PARAMETER( ratmin = 1.D-8) 
!
      j  = ( ij - 1 ) / ndi + 1
      i  = MOD( ( ij - 1 ), ndi) + 1
!
      IF(nt.GT.0) THEN
!
! ... iterative inversion of the enthalpy-temperature law,
! ... the gas thermal capacity depending on the temperature (cg=cg(T))
        DO ii = 1, nlmax
          tgnn = tg
          CALL hcapg(cp(:), tg)
          c2=0.D0
          DO kg=1,ngas
            c2=c_erg*cp(kg)*ygc(kg)+c2
          END DO
          cg = c2
          tg = tzero+(sieg-hzerog)/cg
          IF (DABS((tgnn-tg)/tgnn).LE.ratmin) GOTO 223
        END DO
        WRITE(8,*) 'max number of iteration reached in eosg'
        WRITE(8,*) 'time, i, j',time,i,j
  223   CONTINUE
      ENDIF
!
      IF(nr.GT.0) THEN
        c1=0.D0
        DO kg=1,ngas
          c1 = c1 + xgc(kg) * gmw(kg)
        END DO
        rog = p/(rgas*tg)*c1
      ENDIF
!
      IF(nc.GT.0) rags = rog/p
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE cnvertg(ij)
!
! ... computes thermodynamic mean quantities
!
      USE dimensions
      USE gas_constants, ONLY: gmw, c_erg,rgas,tzero,hzerog
      USE gas_solid_density, ONLY: rog_g, rgp_g, rgpn_g
      USE gas_solid_temperature, ONLY: sieg_g, siegn_g, tg_g
      USE pressure_epsilon, ONLY: p_g, ep_g
      USE reactions, ONLY: irex
      USE th_capacity, ONLY: cp_g, ck_g, hcapg
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: ij
      INTEGER :: kg
!
      REAL*8 :: c1, c2
!
      c1 = 0.D0
      DO kg = 1, ngas
        c1 = xgc_g(kg,ij) * gmw(kg) + c1
      END DO

! ... mean gas density
!
      rog_g(ij) = p_g(ij) / (rgas*tg_g(ij)) * c1
      rgp_g(ij) = rog_g(ij) * ep_g(ij)
      rgpn_g(ij) = rgp_g(ij)

      DO kg = 1, ngas
        rgpgc_g(kg,ij) = ygc_g(kg,ij) * rgp_g(ij)
        rgpgcn_g(kg,ij) = rgpgc_g(kg,ij)
      END DO

!pdac---------------
! control next statement
!      IF(irex.EQ.0) RETURN
!pdac---------------
!
! compute heat capacity (constant volume) for gas mixture
!
      CALL hcapg(cp_g(:,ij),tg_g(ij))
      c2=0.D0
      DO kg=1,ngas
        c2=c_erg*cp_g(kg,ij)*ygc_g(kg,ij)+c2
      END DO 

      cg_g(ij)=c2
      sieg_g(ij)=(tg_g(ij)-tzero)*cg_g(ij)+hzerog
      siegn_g(ij)=sieg_g(ij)
!
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      END MODULE eos_gas
