!-----------------------------------------------------------------------
      MODULE vent_conditions
!-----------------------------------------------------------------------
      USE dimensions, ONLY: max_nsolid, max_ngas
      IMPLICIT NONE

      INTEGER :: ivent
      REAL*8 :: xvent, yvent, radius
      REAL*8 :: wrat
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
      USE grid, ONLY: bottom, iv, jv, kv, grigen
      USE grid, ONLY: center_x, center_y
      USE volcano_topography, ONLY: itp

      IMPLICIT NONE
      
      INTEGER :: i, j, k
      INTEGER :: ijk, nv
      INTEGER :: iwest, ieast, jnorth, jsouth
      REAL*8 :: quota
      
      IF( job_type == '2D') RETURN
!
      IF (grigen >= 1) THEN
        xvent = center_x
        yvent = center_y
      END IF
!
      WRITE(7,*) 'Vent conditions imposed in cells: '
!
! ... define the rectangle containing the vent
! ... 'nvt' is the number of vent cells
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
! ... define the 'quota' of the volcanic vent
! ... (considering the topography)
!
      IF( itp < 1 ) THEN
        iv = (ieast  + iwest ) / 2
        jv = (jnorth + jsouth) / 2
        DO k= 1, nz
          ijk = iv + (jv-1) * nx + (k-1) * nx * ny
          IF (fl(ijk) == 3) kv = k
        END DO
      END IF
      quota = kv

      WRITE(6,*) 'vent center: ', iv, jv, kv
      WRITE(6,*) 'center coordinates: ', x(iv), y(jv), z(kv)
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
          DO k = 1, quota - 1
            ijk = i + (j-1) * nx + (k-1) * nx * ny
            fl(ijk) = bottom
          END DO
!
! ... vent cells
!
          k = quota
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
          WRITE(7,10) nv, ijk, i, j, k, vcell(nv)%frac, fl(ijk)
 10       FORMAT(I3,I7,3(I3),F8.4,I2)
!
! ... fluid cells above the vent
!
          DO k = quota+1, nz-1
            ijk = i + (j-1) * nx + (k-1) * nx * ny
            fl(ijk) = 1
          END DO
          
        END DO
      END DO
      
      WRITE(19,'(80(I3))') (fl(ijk), ijk=1,nx*ny*nz)
!
      RETURN
      END SUBROUTINE locate_vent
!-----------------------------------------------------------------------
      SUBROUTINE set_ventc
!
! ... Compute the steady inlet conditions for a circular vent
! ... An aribtrary radial profile can be assigned, as a function
! ... of the averaged vertical velocity
!
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
      REAL*8 :: alpha, ep0, ra, beta
      INTEGER :: ijk, imesh, i,j,k, is, ig, n

      IF (job_type == '2D') RETURN
      
      DO ijk = 1, ncint      
        IF(flag(ijk) == 8) THEN
          CALL meshinds(ijk,imesh,i,j,k)

          IF (wrat > 1.D0) THEN
            beta = 1.D0 / (wrat - 1.D0)
            ra = DSQRT((x(i)-xvent)**2 + (y(j)-yvent)**2)
            ra = ra / radius
          ELSE IF (wrat <= 1.D0) THEN
            wrat = 1.D0
            beta = 1.D0
            ra = 0.D0
          END IF

          DO n = 1, nvt
            IF (vcell(n)%imesh == imesh) THEN
              alpha = vcell(n)%frac
            END IF
          END DO
          
          ug(ijk) = u_gas 
          IF (job_type == '3D') vg(ijk) = v_gas 
          
          !
          ! ... vertical velocity profile
          !
          wg(ijk) = w_gas * wrat * (1.D0 - ra ** beta)
          
          tg(ijk) = t_gas
          ep0     = 1.D0 - SUM(ep_solid(1:nsolid))
          ep(ijk) = 1.D0 - alpha * SUM(ep_solid(1:nsolid))
          p(ijk)  = p_gas * alpha * ep0 / ep(ijk)

          DO ig = 1, ngas
            ygc(ijk,ig) = vent_ygc(gas_type(ig))
          END DO

          DO is = 1,nsolid

            us(ijk,is)  = u_solid(is) 
            IF (job_type == '3D') vs(ijk,is)  = v_solid(is) 
          
            !
            ! ... vertical velocity profile
            !
            ws(ijk,is) = w_solid(is) * wrat * (1.D0 - ra ** beta)

            ws(ijk,is)  = w_solid(is) 
            ts(ijk,is)  = t_solid(is)
            rlk(ijk,is) = ep_solid(is)*rl(is) * alpha
          END DO
          !
          ! ... check gas components closure relation
          ygcsum = SUM(ygc(ijk,:))
          IF ( ygcsum /= 1.D0 ) THEN
            ygc(ijk,ngas) = 1.D0 - SUM( ygc(ijk,1:ngas-1) )
          END IF
          
        END IF
      END DO

      DEALLOCATE(vcell)

      RETURN
      END SUBROUTINE set_ventc
!-----------------------------------------------------------------------
      SUBROUTINE update_ventc(sweep,ijk)
!
! ... Update the vertical velocity profile as a function of time
! ... The time scale is 'tau'
!
      USE control_flags, ONLY: job_type
      USE dimensions, ONLY: nsolid
      USE domain_decomposition, ONLY: ncint, meshinds
      USE gas_solid_velocity, ONLY: ug, wg, vg
      USE gas_solid_velocity, ONLY: us, vs, ws
      USE grid, ONLY: flag
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: sweep, ijk
      REAL*8 :: f1, f2, fact
      INTEGER :: is

      IF (job_type == '2D') RETURN
      
      f1 = ft(sweep)
      f2 = ft(sweep-1)
      IF (sweep == 1) THEN
        fact = f1
      ELSE
        fact = f1/f2
      END IF
      
      wg(ijk) = wg(ijk) * fact
      DO is = 1,nsolid
        ws(ijk,is) = ws(ijk,is) * fact
      END DO

      RETURN
      END SUBROUTINE update_ventc
!-----------------------------------------------------------------------
      REAL*8 FUNCTION ft(n)
      USE time_parameters, ONLY: tau, dt
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n
      REAL*8 :: t

      t = n * dt

      IF (t <= tau) THEN
        ft = 3.D0 * (t / tau)**2 - 2.D0 * (t / tau)**3     
      ELSE
        ft = 1.D0
      END IF

      RETURN
      END FUNCTION ft
!-----------------------------------------------------------------------
      REAL*8 FUNCTION cell_fraction(i, j)
!
! ... assign a weight to partially filled vent cells:
! ... the weight is proportional to the area of the intersection of
! ... the cell and the vent 
!
      USE grid, ONLY: xb, yb, zb
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: i, j

	TYPE polypoint
		REAL*8 :: x
		REAL*8 :: y
		INTEGER :: side		!(0=corner, 1=W, 2=N, 3=E, 4=S)
	END TYPE polypoint
      TYPE(polypoint) :: poly_p(16) !intersection points of bound.cell and vent 
				    !(from west-south corner, clowckwise)
	INTEGER :: n,l

	REAL*8 :: delta_sol,x_sol1,x_sol2,y_sol1,y_sol2
	REAL*8 :: area_p

	n=0
	
	
	IF ( (xb(i-1)-xvent)**2 + (yb(j-1)-yvent)**2 <= radius**2 ) THEN
		n=n+1
		poly_p(n)%x=xb(i-1)
		poly_p(n)%y=yb(j-1)
		poly_p(n)%side=0
	ENDIF
	
	
	delta_sol=radius**2-(xb(i-1)-xvent)**2
	
	IF ( delta_sol > 0 ) THEN
		y_sol1= yvent-SQRT (delta_sol)
		y_sol2= yvent+SQRT (delta_sol)
		IF ( (yb(j-1) < y_sol1).AND. (y_sol1 < yb(j)) ) THEN
			n=n+1
			poly_p(n)%x=xb(i-1)
			poly_p(n)%y=y_sol1
			poly_p(n)%side=1
		ENDIF 
		
		IF ( (yb(j-1) < y_sol2) .AND. (y_sol2 < yb(j)) ) THEN
			n=n+1
			poly_p(n)%x=xb(i-1)
			poly_p(n)%y=y_sol2
			poly_p(n)%side=1
		ENDIF 
	ENDIF

	IF ( (xb(i-1)-xvent)**2 + (yb(j)-yvent)**2 <= radius**2 ) THEN
		n=n+1
		poly_p(n)%x=xb(i-1)
		poly_p(n)%y=yb(j)
		poly_p(n)%side=0
	ENDIF
	
	

	delta_sol=radius**2-(yb(j)-yvent)**2
	
	IF ( delta_sol > 0 ) THEN
		x_sol1= xvent-SQRT (delta_sol)
		x_sol2= xvent+SQRT (delta_sol)
	
		IF ( (xb(i-1) < x_sol1).AND. (x_sol1 < xb(i)) ) THEN
			n=n+1
			poly_p(n)%x=x_sol1
			poly_p(n)%y=yb(j)
			poly_p(n)%side=2
		ENDIF 
		
		IF ( (xb(i-1) < x_sol2) .AND. (x_sol2 < xb(i)) ) THEN
			n=n+1
			poly_p(n)%x=x_sol2
			poly_p(n)%y=yb(j)
			poly_p(n)%side=2
		ENDIF 
	ENDIF

	IF ( (xb(i)-xvent)**2 + (yb(j)-yvent)**2 <= radius**2 ) THEN
		n=n+1
		poly_p(n)%x=xb(i)
		poly_p(n)%y=yb(j)
		poly_p(n)%side=0
	ENDIF



	delta_sol=radius**2-(xb(i)-xvent)**2
	
	IF ( delta_sol > 0 ) THEN
		y_sol1= yvent+SQRT (delta_sol)
		y_sol2= yvent-SQRT (delta_sol)
	
		IF ( (yb(j-1) < y_sol1) .AND. (y_sol1 < yb(j)) ) THEN
			n=n+1
			poly_p(n)%x=xb(i)
			poly_p(n)%y=y_sol1
			poly_p(n)%side=3
		ENDIF 
		
		IF ( (yb(j-1) < y_sol2) .AND. ( y_sol2 < yb(j) ) ) THEN
			n=n+1
			poly_p(n)%x=xb(i)
			poly_p(n)%y=y_sol2
			poly_p(n)%side=3
		ENDIF 
	ENDIF
	
	
	IF ( (xb(i)-xvent)**2 + (yb(j-1)-yvent)**2 <= radius**2 ) THEN
		n=n+1
		poly_p(n)%x=xb(i)
		poly_p(n)%y=yb(j-1)
		poly_p(n)%side=0
	ENDIF
	
	

	delta_sol=radius**2-(yb(j-1)-yvent)**2

	IF ( delta_sol > 0 ) THEN
		x_sol1= xvent+SQRT (delta_sol)
		x_sol2= xvent-SQRT (delta_sol)
		IF ( (xb(i-1) < x_sol1).AND. (x_sol1 < xb(i)) ) THEN
			n=n+1
			poly_p(n)%x=x_sol1
			poly_p(n)%y=yb(j-1)
			poly_p(n)%side=4
		ENDIF 
		
		IF ( (xb(i-1) < x_sol2) .AND. (x_sol2 < xb(i)) ) THEN
			n=n+1
			poly_p(n)%x=x_sol2
			poly_p(n)%y=yb(j-1)
			poly_p(n)%side=4
		ENDIF 
	ENDIF
	
	DO l=1,n
		poly_p(n+l)%x=poly_p(l)%x
		poly_p(n+l)%y=poly_p(l)%y
		poly_p(n+l)%side=poly_p(l)%side
	ENDDO

	area_p=0.d0

	DO l=1,n-2
		area_p=area_p+areaT(poly_p(1)%x,poly_p(l+1)%x,poly_p(l+2)%x,   &
			            poly_p(1)%y,poly_p(l+1)%y,poly_p(l+2)%y)
	ENDDO
	
	DO l=1,n
		IF  ( ( poly_p(l)%side*poly_p(l+1)%side > 0 ) .AND. &
		      ( poly_p(l)%side /= poly_p(l+1)%side ) ) THEN
			area_p = area_p + areaSEC(poly_p(l)%x,poly_p(l+1)%x,  &
				poly_p(l)%y,poly_p(l+1)%y,xvent,yvent,radius) &
			    - areaT(poly_p(l)%x,poly_p(l+1)%x,xvent,poly_p(l)%y, &
				poly_p(l+1)%y,yvent)

		ENDIF
	ENDDO


	cell_fraction= area_p / abs( (xb(i-1)-xb(i))*(yb(j-1)-yb(j)) ) 




      RETURN
     
      END FUNCTION cell_fraction
!-----------------------------------------------------------------------
      REAL*8 FUNCTION areaT(x1,x2,x3,y1,y2,y3)

! ... evaluate the area of a triangle

	
	REAL*8 :: x1,x2,x3,y1,y2,y3
	REAL*8 :: p,a,b,c
	
	a=SQRT( (x1-x2)**2 + (y1-y2)**2 )
	b=SQRT( (x1-x3)**2 + (y1-y3)**2 )
	c=SQRT( (x3-x2)**2 + (y3-y2)**2 )
	p= (a+b+c) / 2.d0
	
	areaT= SQRT ( p * (p-a) * (p-b) * (p-c) )
	
     RETURN

	END FUNCTION areaT
	
!-----------------------------------------------------------------------
	REAL*8 FUNCTION areaSEC(x1,x2,y1,y2,xC,yC,R)

! ... evaluate the area of a sector of circle
	
	REAL*8 :: x1,x2,y1,y2,xC,yC,R
	REAL*8 :: alpha1,alpha2,alpha
        REAL*8 :: pi

        pi = 4.D0 * datan(1.D0)
	
	alpha1= datan2( (y1-yC)/R , (x1-xC)/R )
	alpha2= datan2( (y2-yC)/R , (x2-xC)/R )

	alpha= ABS (alpha1-alpha2)
	IF ( alpha > pi ) alpha = 2.D0*pi - alpha
	
	areaSEC = 0.5D0*alpha*R**2

        !WRITE(7,*) alpha1, alpha2, alpha

	RETURN
	
	END FUNCTION areaSEC


      END MODULE vent_conditions
!-----------------------------------------------------------------------
