!------------------------------------------------------------------------
      MODULE atmospheric_conditions
!------------------------------------------------------------------------
      USE dimensions, ONLY: max_ngas
      USE io_files, ONLY: errorunit, logunit
      IMPLICIT NONE
      PUBLIC
!
! ... Initial atmospheric variables
!
      REAL*8 :: wind_x
      REAL*8 :: wind_y
      REAL*8 :: wind_z
      !
      ! ... Pressure and temperature at sea level
      REAL*8 :: p_ground  
      REAL*8 :: t_ground

      REAL*8 :: void_fraction
      REAL*8 :: atm_ygc(max_ngas)

      REAL*8 :: max_packing
      REAL*8 :: gravx, gravy, gravz
      REAL*8, ALLOCATABLE :: p_atm(:), t_atm(:)

      ! ... Flag for temperature stratification
      !
      LOGICAL :: stratification

      TYPE atmospheric_layer
        REAL*8 :: ztop
        REAL*8 :: ptop
        REAL*8 :: ttop
        REAL*8 :: gradt
        CHARACTER(LEN=20) :: name
      END TYPE

      INTEGER, PARAMETER :: num_layers = 7
      TYPE(atmospheric_layer) :: layer(num_layers)

      SAVE
!------------------------------------------------------------------------
      CONTAINS
!------------------------------------------------------------------------
      SUBROUTINE control_atmosphere

      USE parallel, only: mpime, root

      IMPLICIT NONE

      IF (gravz == 0.0) THEN
         stratification = .FALSE.
      ELSE IF (gravz /= -9.81D0) THEN
         IF( mpime == root ) THEN
           WRITE(errorunit,*) 'WARNING! from control atmosphere'
           WRITE(errorunit,*) 'gravz = ', gravz
         END IF
      END IF
      IF( mpime == root ) THEN
        WRITE(logunit,*) 
        WRITE(logunit,*) 'Temperature stratification: ', stratification
        WRITE(logunit,*) 'Gravity: ', gravz
      END IF
!
      END SUBROUTINE control_atmosphere
!------------------------------------------------------------------------
      SUBROUTINE set_atmosphere
!------------------------------------------------------------------------
! ... This routine computes atmospheric stratification accordingly with
! ... the description of a standard, quiet atmosphere
! ... For each atmospheric region, the upper limit (in metres) and the
! ... Temperature gradient is set from input.
!
      USE control_flags, ONLY: lpr
      USE dimensions, ONLY: nz
      USE grid, ONLY: dz
      USE parallel, only: mpime, root
!
      IMPLICIT NONE

      REAL*8 :: zbot, pbot, tbot
      REAL*8 :: ztop, ptop, ttop
      REAL*8 :: gradt

      INTEGER :: l
!
      ALLOCATE(p_atm(nz), t_atm(nz))
!
      ptop = p_ground
      ttop = t_ground
!
! ... For each layer, compute the bottom and top 
! ... pressure and temperature
!
      DO l = 1, num_layers

        IF (l>1) THEN
          zbot = layer(l-1)%ztop
        ELSE
          zbot = 0.D0
        END IF
        ztop = layer(l)%ztop
        gradt = layer(l)%gradt

        pbot = ptop
        tbot = ttop
!
! ... Temperature is assumed to vary linearly
! ... Pressure is computed by assuming hydrostatic equilibrium
!
        CALL hydrostatic(zbot, ztop, gradt, tbot, pbot, ttop, ptop)

        layer(l)%ptop = ptop
        layer(l)%ttop = ttop

        IF (lpr > 1 .AND. stratification) THEN
          IF( mpime == root ) THEN
            WRITE(logunit,*) layer(l)%name
            WRITE(logunit,*) layer(l)%gradt
            WRITE(logunit,*) layer(l)%ztop
            WRITE(logunit,*) pbot, layer(l)%ptop
            WRITE(logunit,*) tbot, layer(l)%ttop
          END IF
        END IF

      END DO
!
! ... Compute the discrete atmospheric profile
!
      CALL compute_profile

      END SUBROUTINE set_atmosphere
!------------------------------------------------------------------------
      SUBROUTINE compute_profile
!------------------------------------------------------------------------
! ... This routine computes the standard atmospheric condition in each
! ... cell of the computational domain 
!
      USE dimensions, ONLY: nz
      USE grid, ONLY: dz, zb, zzero
      USE parallel, only: mpime, root
!
      IMPLICIT NONE
!
      REAL*8 :: za, pa, ta
      REAL*8 :: zbot, pbot, tbot
      REAL*8 :: ztop, ptop, ttop
      REAL*8 :: gradt
      INTEGER :: k, l
!
! ... First layer
!
      IF( zzero <= layer(1)%ztop ) THEN
        l = 1
      ELSE
        DO l = 2, num_layers
          IF( zzero > layer(l-1)%ztop ) EXIT
        END DO
      END IF
!
! ... Loop vertically over mesh layers
!
      DO k = 1, nz
        za = zb(k) + 0.5D0*(dz(1)-dz(k))
!
        IF (za > 0.D0) THEN
          IF (.NOT.stratification) THEN
            zbot = zzero
            pbot = p_ground
            tbot = t_ground
            gradt = 0.D0
          ELSE
            IF (za > layer(l)%ztop) l=l+1
            IF (l==1) THEN
              zbot = zzero
              pbot = p_ground
              tbot = t_ground
            ELSE IF (l>1) THEN
              zbot = layer(l-1)%ztop
              pbot = layer(l-1)%ptop
              tbot = layer(l-1)%ttop
            ELSE IF (l > num_layers) THEN
              CALL error('atm','altitude out of range',1)
            END IF
            gradt = layer(l)%gradt
          END IF
          !
          ! ... Temperature is assumed to vary linearly
          ! ... Pressure is computed by assuming hydrostatic equilibrium
          !
          CALL hydrostatic(zbot, za, gradt, tbot, pbot, ta, pa)
!
        ELSE IF (za <= 0.D0) THEN
          IF ( mpime == root ) THEN
            WRITE(errorunit,*) 'WARNING! from atmospheric profile:'
            WRITE(errorunit,*) ' Row ',k, ' lays below the sea level; z = ', za
          END IF
          ta = t_ground
          pa = p_ground
        END IF
!
! ... Store the cell values of atmospheric pressure and temperature 
!
        t_atm(k) = ta
        p_atm(k) = pa

      END DO
!
      RETURN
      END SUBROUTINE compute_profile
!----------------------------------------------------------------------
      SUBROUTINE hydrostatic(zb, z, grad, tb, pb, t, p)
!
! ... integrate the hydrostatic law  dP/dz = rho(z) * gravz
! ... using the thermal equation of state to express the 
! ... density as a function of P(z) and T(z) and
! ... by assuming that T(z) = T0 + gradt * (z-z0)
! ... If gravity is ZERO, the pressure is constant.
!
      USE gas_constants, ONLY: gmw,rgas

      IMPLICIT NONE

      REAL*8, INTENT(IN) :: zb, z, grad, tb, pb 
      REAL*8, INTENT(OUT) :: t, p
      REAL*8 :: erreair, cost
      
      erreair=rgas/gmw(6)

      IF (grad /= 0.D0) THEN
        cost = gravz/(erreair*grad)
        t = tb + grad*(z-zb)
        p = pb * (t/tb)**cost
      ELSE IF (grad == 0.D0) THEN
        cost = gravz/(erreair*tb)
        t = tb  
        p = pb * DEXP(cost*(z-zb))
      END IF

      RETURN
      END SUBROUTINE hydrostatic
!----------------------------------------------------------------------
      END MODULE atmospheric_conditions
!------------------------------------------------------------------------
