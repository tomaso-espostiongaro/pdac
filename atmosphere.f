!------------------------------------------------------------------------
      MODULE atmospheric_conditions
!------------------------------------------------------------------------
      USE dimensions, ONLY: max_ngas
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

      stratification = .TRUE.
      IF (gravz == 0.0) THEN
         stratification = .FALSE.
      ELSE IF (gravz /= -9.81D0) THEN
         IF( mpime == root ) THEN
           WRITE(8,*) 'WARNING! from control atmosphere'
           WRITE(8,*) 'gravz = ', gravz
         END IF
      END IF
      IF( mpime == root ) THEN
        WRITE(6,*) 
        WRITE(6,*) 'Atmospheric stratification: ', stratification
        WRITE(6,*) 'Gravity: ', gravz
      END IF
!
      IF (stratification) THEN
        IF( mpime == root ) THEN
          IF (t_ground /= 288.15D0) WRITE(8,*) 'WARNING! control atmospheric &
                                          & temperature profile'
          IF (p_ground /= 1.01325D5) WRITE(8,*)'WARNING! control atmospheric &
                                          & pressure profile'
        END IF
      END IF	
!
      END SUBROUTINE control_atmosphere
!------------------------------------------------------------------------
      SUBROUTINE set_atmosphere
!------------------------------------------------------------------------
! ... This routine computes atmospheric stratification accordingly with
! ... the description of a standard, quiet atmosphere
! ... For each atmospheric region, the upper limit (in metres) and the
! ... Temperature gradient must be set.
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
! ... Initialize atmospheric layers
!
      layer(1)%name = 'Troposphere'   
      layer(2)%name = 'Tropopause'    
      layer(3)%name = 'Lower_Stratosphere' 
      layer(4)%name = 'Upper_Stratosphere' 
      layer(5)%name = 'Ozone_layer'   
      layer(6)%name = 'Lower_Mesosphere'   
      layer(7)%name = 'Upper_Mesosphere'   
!
! ... Top of the layer (height a.s.l.)
!
      layer(1)%ztop = 1.1D4
      layer(2)%ztop = 2.0D4
      layer(3)%ztop = 3.2D4
      layer(4)%ztop = 4.7D4
      layer(5)%ztop = 5.1D4
      layer(6)%ztop = 7.1D4
      layer(7)%ztop = 8.0D4
!
! ... Temperature gradient (T is assumed to vary linearly)
!
      layer(1)%gradt = -6.5D-3
      layer(2)%gradt = 0.D0
      layer(3)%gradt = 1.D-3
      layer(4)%gradt = 2.8D-3
      layer(5)%gradt = 0.D0
      layer(6)%gradt = -2.8D-3
      layer(7)%gradt = -2.0D-3
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
            WRITE(6,*) layer(l)%name
            WRITE(6,*) layer(l)%gradt
            WRITE(6,*) layer(l)%ztop
            WRITE(6,*) pbot, layer(l)%ptop
            WRITE(6,*) tbot, layer(l)%ttop
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

      DO k = 1, nz

        za = zb(k) + 0.5D0*(dz(1)-dz(k))

        IF ( ( za < 0.D0 ) .AND. ( mpime == root ) ) & 
          WRITE(8,*) 'WARNING! from atmospheric profile:'
          WRITE(8,*) ' Row ',k, ' lays below the sea level; z = ', za
        IF (za <= 0.D0 .OR. .NOT.stratification) THEN
          ta = t_ground
          pa = p_ground
          GOTO 100
        ELSE IF (za > layer(l)%ztop) THEN
          l = l + 1
        END IF

        IF (l > num_layers) CALL error('atm','altitude out of range',1)

        IF (l>1) THEN
          zbot = layer(l-1)%ztop
          pbot = layer(l-1)%ptop
          tbot = layer(l-1)%ttop
        ELSE
          zbot = 0.D0
          pbot = p_ground
          tbot = t_ground
        END IF
        gradt = layer(l)%gradt
!
! ... Temperature is assumed to vary linearly
! ... Pressure is computed by assuming hydrostatic equilibrium
!
        CALL hydrostatic(zbot, za, gradt, tbot, pbot, ta, pa)

 100    CONTINUE
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
