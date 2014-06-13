!------------------------------------------------------------------------
!>    Define Atmospheric conditions
!------------------------------------------------------------------------
      MODULE atmospheric_conditions
!------------------------------------------------------------------------
      USE dimensions, ONLY: max_ngas
      USE io_files, ONLY: errorunit, atmounit, atmofile
      IMPLICIT NONE
      PUBLIC
!
      REAL*8 :: wind_x  !< x-component of wind velocity
      REAL*8 :: wind_y  !< y-component of wind velocity
      REAL*8 :: wind_z  !< z-component of wind velocity
      !
      REAL*8 :: p_ground  !< Pressure at sea level
      REAL*8 :: t_ground  !< Temperature at sea level

      REAL*8 :: void_fraction !< Gas volume fraction in atmosphere
      REAL*8 :: atm_ygc(max_ngas) !< Mass fractions of atmospheric gases

      REAL*8 :: max_packing !< Maximum volumetric packing of particles
      REAL*8 :: gravx, gravy, gravz !< Components of gravity acceleration
      REAL*8, ALLOCATABLE :: p_atm(:) !< Atmospheric pressure at grid centers
      REAL*8, ALLOCATABLE :: t_atm(:) !< Atmospheric temperature at grid centers

      !> Flag for temperature stratification
      LOGICAL :: stratification

      !> Atmospheric layers
      TYPE atmospheric_layer
        REAL*8 :: ztop !< Height above sea level of layer top
        REAL*8 :: ptop !< Pressure at layer top
        REAL*8 :: ttop !< Temperature at latyer top
        REAL*8 :: gradt !< Temperature gradient [K/m]
        CHARACTER(LEN=20) :: name !< Layer name
      END TYPE

      !> Number of atmospheric layers
      INTEGER, PARAMETER :: num_layers = 7
      TYPE(atmospheric_layer) :: layer(num_layers)

      SAVE
!------------------------------------------------------------------------
      CONTAINS
!------------------------------------------------------------------------
!> Check stratification conditions.
!------------------------------------------------------------------------
      SUBROUTINE control_atmosphere
!------------------------------------------------------------------------
!
      USE control_flags, ONLY: lpr
      USE parallel, only: mpime, root
!
      IMPLICIT NONE
!
      IF (gravz == 0.0) THEN
         stratification = .FALSE.
      ELSE IF (gravz /= -9.81D0) THEN
         IF( mpime == root ) THEN
           WRITE(errorunit,*) 'WARNING! from control atmosphere'
           WRITE(errorunit,*) 'gravz = ', gravz
         END IF
      END IF
      IF(mpime == root ) THEN
        OPEN(UNIT=atmounit,FILE=atmofile,STATUS='UNKNOWN')
        WRITE(atmounit,*) 'Report of atmospheric initial conditions'
        WRITE(atmounit,*) 
        WRITE(atmounit,*) 'Temperature stratification: ', stratification
        WRITE(atmounit,*) 'Gravity: ', gravz
        WRITE(atmounit,*) 
        IF (stratification) WRITE(atmounit,*) 'Atmospheric layers:'
      END IF
!
      END SUBROUTINE control_atmosphere
!------------------------------------------------------------------------
!> Compute atmospheric stratification accordingly with
!> the description of a standard, quiet atmosphere.
!> For each atmospheric layer, the upper limit (in metres) and the
!> temperature gradient is set from input.
!------------------------------------------------------------------------
      SUBROUTINE set_atmosphere
!------------------------------------------------------------------------
!
      USE control_flags, ONLY: lpr
      USE dimensions, ONLY: nz
      USE grid, ONLY: dz, kv
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
! ... The atmospheric profile starts at sea level (z=0)
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

        IF (stratification) THEN
          IF( mpime == root ) THEN
            WRITE(atmounit,*) layer(l)%name
            WRITE(atmounit,*) layer(l)%gradt
            WRITE(atmounit,*) layer(l)%ztop
            WRITE(atmounit,*) pbot, layer(l)%ptop
            WRITE(atmounit,*) tbot, layer(l)%ttop
          END IF
        END IF
      END DO
!
! ... Compute the discrete atmospheric profile
!
      CALL compute_profile

      END SUBROUTINE set_atmosphere
!------------------------------------------------------------------------
!> Compute the standard atmospheric condition in each
!> cell of the computational domain 
!------------------------------------------------------------------------
      SUBROUTINE compute_profile
!------------------------------------------------------------------------
!
      USE control_flags, ONLY: lpr
      USE dimensions, ONLY: nz
      USE grid, ONLY: dz, z
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
! ... Identify the layer ('l') where the mesh starts
!
      IF( z(1) <= layer(1)%ztop ) THEN
        l = 1
      ELSE
        DO l = 2, num_layers
          IF( z(1) > layer(l-1)%ztop ) EXIT
        END DO
      END IF
!
! ... Loop vertically over mesh layers
!
      IF ( mpime == root ) THEN
              WRITE(atmounit,*) 
              WRITE(atmounit,*) 'Vertical stratification (k, z, p, t): '
      END IF
!
      DO k = 1, nz
        za = z(k)
!
        IF (za > 0.D0) THEN
          IF (.NOT.stratification) THEN
            zbot = 0.D0
            pbot = p_ground
            tbot = t_ground
            gradt = 0.D0
          ELSE
            IF (za > layer(l)%ztop) l=l+1
            IF (l==1) THEN
              zbot = 0.D0
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
!
        IF ( mpime == root ) THEN
          WRITE(atmounit,'(I5,3F14.2)') k, za, pa, ta
        END IF
!
      END DO
!
      RETURN
      END SUBROUTINE compute_profile
!----------------------------------------------------------------------
!> Integrate the hydrostatic law \f$ dP/dz = \rho(z) * g_z \f$
!> using the thermal equation of state to express the 
!> density \f$ \rho(z)=\rho(P(z), T(z)) \f$ and
!> by assuming that \f$ T(z) = T_0 + \nabla T \cdot (z-z0) \f$.

!> If gravity is ZERO, pressure is constant.
!----------------------------------------------------------------------
      SUBROUTINE hydrostatic(zb, z, grad, tb, pb, t, p)
!----------------------------------------------------------------------
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
