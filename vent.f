!-----------------------------------------------------------------------
      MODULE vent_conditions
!-----------------------------------------------------------------------
      USE dimensions, ONLY: max_nsolid, max_ngas
      IMPLICIT NONE

      INTEGER :: ivent
      REAL*8 :: xvent, yvent, radius
      REAL*8 :: u_gas, v_gas, w_gas, p_gas, t_gas
      REAL*8 :: u_solid(max_nsolid), v_solid(max_nsolid), w_solid(max_nsolid), &
                ep_solid(max_nsolid), t_solid(max_nsolid)
      REAL*8 :: vent_ygc(max_ngas)

      TYPE inlet_cell
        INTEGER :: imesh
        REAL*8  :: frac
      END TYPE inlet_cell
      TYPE(inlet_cell), ALLOCATABLE :: vcell(:)
      INTEGER :: nvt

      SAVE
!-----------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------
      SUBROUTINE locate_vent

      USE control_flags, ONLY: job_type
      USE dimensions, ONLY: nx, ny, nz
      USE grid, ONLY: x, y, z, fl, xb, yb, zb
      USE grid, ONLY: bottom, kv

      IMPLICIT NONE
      
      INTEGER :: i, j, k
      INTEGER :: ijk, nv
      INTEGER :: iwest, ieast, jnorth, jsouth
      REAL*8 :: dist2
      
      IF( job_type == '2D') RETURN
!
      WRITE(6,*) 'Vent conditions imposed in cells: '
!      
      DO i = 2, nx
        IF (xb(i-1) <= (xvent-radius)) iwest = i
        IF (xb(i-1) < (xvent+radius)) ieast = i
      END DO
      DO j = 2, ny
        IF (yb(j-1) <= (yvent-radius)) jsouth = j
        IF (yb(j-1) < (yvent+radius)) jnorth = j
      END DO
!
      nvt = (ieast-iwest+1)*(jnorth-jsouth+1)
      ALLOCATE(vcell(nvt))

      nv = 0
      DO j = jsouth, jnorth
        DO i = iwest, ieast
          nv = nv + 1
!
! ... topography below the vent 
!
          DO k = 1, kv-1
            ijk = i + (j-1) * nx + (k-1) * nx * ny
            fl(ijk) = bottom
          END DO
!
! ... vent cells
!
          k = kv
          ijk = i + (j-1) * nx + (k-1) * nx * ny
          !
          vcell(nv)%imesh = ijk
          vcell(nv)%frac = cell_fraction(i,j)
          IF (vcell(nv)%frac > 0.D0) THEN
            fl(ijk) = 8
          ELSE
            fl(ijk) = bottom
          END IF
          !
          WRITE(6,*) nv, ijk, i, j, k, vcell(nv)%frac, fl(ijk)
!
! ... fluid cells above the vent
!
          DO k = kv+1, nz-1
            ijk = i + (j-1) * nx + (k-1) * nx * ny
            fl(ijk) = 1
          END DO
          
        END DO
      END DO
!
      RETURN
      END SUBROUTINE locate_vent
!-----------------------------------------------------------------------
      SUBROUTINE set_ventc
! ... 
      USE atmospheric_conditions, ONLY: p_ground, t_ground
      USE control_flags, ONLY: job_type
      USE dimensions, ONLY: nsolid, ngas
      USE domain_decomposition, ONLY: ncint, meshinds
      USE eos_gas, ONLY: ygc
      USE gas_constants, ONLY: gas_type
      USE gas_solid_density, ONLY: rlk
      USE gas_solid_temperature, ONLY: tg, ts
      USE gas_solid_velocity, ONLY: ug, wg, vg
      USE gas_solid_velocity, ONLY: us, vs, ws
      USE grid, ONLY: x, y, z, flag
      USE particles_constants, ONLY: rl, inrl
      USE pressure_epsilon, ONLY: ep, p
      IMPLICIT NONE

      REAL*8 :: area, ygcsum
      REAL*8 :: mixdens, mixvel, mfr, mg
      REAL*8 :: alpha, ep0
      INTEGER :: ijk, imesh, i,j,k, is, ig, n

      area = ACOS(-1.D0) * radius**2
      
      DO ijk = 1, ncint      
        IF(flag(ijk) == 8) THEN
          CALL meshinds(ijk,imesh,i,j,k)

          DO n = 1, nvt
            IF (vcell(n)%imesh == imesh) THEN
              alpha = vcell(n)%frac
            END IF
          END DO
          
          ug(ijk) = u_gas 
          IF (job_type == '3D') vg(ijk) = v_gas 
          wg(ijk) = w_gas 
          tg(ijk) = t_gas
          ep0     = 1.D0 - SUM(ep_solid(1:nsolid))
          ep(ijk) = 1.D0 - SUM(ep_solid(1:nsolid)) * alpha
          p(ijk)  = p_gas * alpha * ep(ijk) / ep0

          DO ig = 1, ngas
            ygc(ig,ijk) = vent_ygc(gas_type(ig))
          END DO

          DO is = 1,nsolid
            us(ijk,is)  = u_solid(is) 
            IF (job_type == '3D') vs(ijk,is)  = v_solid(is) 
            ws(ijk,is)  = w_solid(is) 
            ts(ijk,is)  = t_solid(is)
            rlk(ijk,is) = ep_solid(is)*rl(is) * alpha
          END DO
          !
          ! ... check gas components closure relation
          ygcsum = SUM(ygc(:,ijk))
          IF ( ygcsum /= 1.D0 ) THEN
            ygc(ngas,ijk) = 1.D0 - SUM( ygc(1:ngas-1,ijk) )
          END IF
          
          WRITE(6,*) i,j,k, alpha
          WRITE(6,*) p(ijk), ep(ijk), rlk(ijk,:)
        END IF
      END DO

      END SUBROUTINE set_ventc
!-----------------------------------------------------------------------
      REAL*8 FUNCTION cell_fraction(i, j)
!
! ... assign a weight to partially filled vent cells:
! ... the weight is proportional to the number of corners
! ... included within the vent
!
      USE grid, ONLY: xb, yb, zb
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: i, j
      REAL*8 :: dist2, alpha
      INTEGER :: n
      TYPE corner
        REAL*8 :: x
        REAL*8 :: y
        INTEGER :: internal
      END TYPE corner
      TYPE(corner) :: cr(4) !(from west-south corner, clowckwise) 

      cr(:)%internal = 0

      cr(1)%x = xb(i-1)
      cr(1)%y = yb(j-1)

      cr(2)%x = xb(i-1)
      cr(2)%y = yb(j)

      cr(3)%x = xb(i)
      cr(3)%y = yb(j)

      cr(4)%x = xb(i)
      cr(4)%y = yb(j-1)

      alpha = 0.D0
      DO n = 1, 4
        dist2 = (cr(n)%x-xvent)**2 + (cr(n)%y-yvent)**2
        IF(dist2 <= radius**2) cr(n)%internal = 1
        alpha = alpha + cr(n)%internal
      END DO
!
      cell_fraction = alpha * 0.25D0
!
      RETURN
      END FUNCTION cell_fraction
!-----------------------------------------------------------------------
      END MODULE vent_conditions
!-----------------------------------------------------------------------
