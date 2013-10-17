!----------------------------------------------------------------------
      MODULE explicit_solver
!
! ... the module implements the explicit solver for multiphase equations
!     in particular drag and heat exchange terms are treated implicitly,
!     whereas other terms are treated explicitly with embedded third order
!     Runge-Kutta scheme that allows for time step adaptation
!
!----------------------------------------------------------------------
      USE control_flags, ONLY: job_type, JOB_TYPE_2D, JOB_TYPE_3D
      USE io_files, ONLY: testunit

      IMPLICIT NONE

      REAL*8 :: limit
!
! ... physical coefficients
!
      REAL*8, DIMENSION(:,:), ALLOCATABLE, PRIVATE :: hv
      REAL*8, DIMENSION(:,:), ALLOCATABLE, PRIVATE  :: kpgv
      REAL*8, DIMENSION(:,:), ALLOCATABLE  :: appu, appv, appw
!
! ... right hand side for RK3 method
!
      REAL*8, DIMENSION(:,:), ALLOCATABLE, PRIVATE  :: rhs3
!
! ... vectors for RK3
!
      REAL*8, DIMENSION(:,:), ALLOCATABLE, PRIVATE :: k1
      REAL*8, DIMENSION(:,:), ALLOCATABLE, PRIVATE :: k2
      REAL*8, DIMENSION(:,:), ALLOCATABLE, PRIVATE :: k3
!
! ... coefficients for Runge-Kutta method
!
      REAL*8, DIMENSION(3), PRIVATE :: rk_b, rk_c
      REAL*8, DIMENSION(3,3), PRIVATE :: rk_a
!
! ... star variables are intermediate values of variables in Runge-Kutta scheme
!
      REAL*8, DIMENSION(:),   ALLOCATABLE, PRIVATE :: rgpstar
      REAL*8, DIMENSION(:,:), ALLOCATABLE, PRIVATE :: rlkstar
      REAL*8, DIMENSION(:),   ALLOCATABLE, PRIVATE :: ugstar
      REAL*8, DIMENSION(:,:), ALLOCATABLE, PRIVATE :: usstar
      REAL*8, DIMENSION(:),   ALLOCATABLE, PRIVATE :: vgstar
      REAL*8, DIMENSION(:,:), ALLOCATABLE, PRIVATE :: vsstar
      REAL*8, DIMENSION(:),   ALLOCATABLE, PRIVATE :: wgstar
      REAL*8, DIMENSION(:,:), ALLOCATABLE, PRIVATE :: wsstar
      REAL*8, DIMENSION(:),   ALLOCATABLE, PRIVATE :: tgstar
      REAL*8, DIMENSION(:,:), ALLOCATABLE, PRIVATE :: tsstar
      REAL*8, DIMENSION(:),   ALLOCATABLE, PRIVATE :: siegstar
      REAL*8, DIMENSION(:,:), ALLOCATABLE, PRIVATE :: siesstar
      REAL*8, DIMENSION(:),   ALLOCATABLE, PRIVATE :: pstar
      REAL*8, DIMENSION(:),   ALLOCATABLE, PRIVATE :: epstar
      REAL*8, DIMENSION(:),   ALLOCATABLE, PRIVATE :: rhogstar
!
! ... variables at time step n
!
      REAL*8, DIMENSION(:),   ALLOCATABLE, PRIVATE :: ugn
      REAL*8, DIMENSION(:,:), ALLOCATABLE, PRIVATE :: usn
      REAL*8, DIMENSION(:),   ALLOCATABLE, PRIVATE :: vgn
      REAL*8, DIMENSION(:,:), ALLOCATABLE, PRIVATE :: vsn
      REAL*8, DIMENSION(:),   ALLOCATABLE, PRIVATE :: wgn
      REAL*8, DIMENSION(:,:), ALLOCATABLE, PRIVATE :: wsn
      REAL*8, DIMENSION(:),   ALLOCATABLE, PRIVATE :: tgn
      REAL*8, DIMENSION(:),   ALLOCATABLE, PRIVATE :: rhon
!
! ... matrices for momentum-x, momentum-y, momentum-z and enthalpy
!     systems of equations
!
      REAL*8, DIMENSION(:,:), ALLOCATABLE, PRIVATE :: au
      REAL*8, DIMENSION(:,:), ALLOCATABLE, PRIVATE :: av
      REAL*8, DIMENSION(:,:), ALLOCATABLE, PRIVATE :: aw
      REAL*8, DIMENSION(:,:), ALLOCATABLE, PRIVATE :: at
!
! ... Cournat number
!
      REAL*8, DIMENSION(:), ALLOCATABLE, PRIVATE :: courant
!
! ... total number of equations
!
      INTEGER :: neq

      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE allocate_fields
!
!... allocate vectors and variables needed by the explicit algorithm
!
      USE domain_mapping, ONLY: ncint, ncdom
      USE dimensions, ONLY: nsolid
      USE gas_solid_viscosity, ONLY: gvisx, gvisy, gvisz, pvisx, pvisy, pvisz
      USE iterative_solver, ONLY: rgfe, rgfn, rgft, rsfe, rsfn, rsft
      USE tilde_momentum, ONLY: ugfe, ugft, ugfn, vgfe, vgft, vgfn, wgfe, wgft, wgfn
      USE tilde_momentum, ONLY: usfe, usft, usfn, vsfe, vsft, vsfn, wsfe, wsft, wsfn
      USE tilde_energy, ONLY: egfe, egft, egfn
      USE tilde_energy, ONLY: hgfe, hgft, hgfn
      USE tilde_energy, ONLY: esfe, esft, esfn
      USE tilde_energy, ONLY: hsfe, hsft, hsfn

!... Physical coefficients
!
      ALLOCATE(kpgv(ncint,nsolid))
      kpgv = 0.D0
      ALLOCATE(hv(ncint,nsolid))
      hv = 0.D0
      ALLOCATE(appu(ncdom, ((nsolid+1)**2+(nsolid+1))/2),   &
               appw(ncdom, ((nsolid+1)**2+(nsolid+1))/2))
      appu = 0.0D0
      appw = 0.0D0
      IF (job_type == JOB_TYPE_3D) THEN
        ALLOCATE(appv(ncdom, ((nsolid+1)**2+(nsolid+1))/2) )
        appv = 0.0D0
      END IF
      ALLOCATE( gvisx(ncint), gvisz(ncint) )
      ALLOCATE( pvisx(ncint,nsolid), pvisz(ncint,nsolid) )
      gvisx = 0.D0; gvisz = 0.D0
      pvisx = 0.D0; pvisz = 0.D0

      IF ( job_type == JOB_TYPE_3D ) THEN
        ALLOCATE( gvisy(ncint) )
        ALLOCATE( pvisy(ncint,nsolid) )
        gvisy = 0.D0
        pvisy = 0.D0
      END IF
!
!... Dimensions of the problem
!
      IF (job_type == JOB_TYPE_2D) THEN
         neq = 4*(nsolid+1)
      ELSE IF (job_type == JOB_TYPE_3D) THEN
         neq = 5*(nsolid+1)
      END IF

      limit = 1.0D-12
!
!... Right hand side of the linear systems (one system on each computational cell)
!
      ALLOCATE(rhs3(ncint,neq))
      rhs3 = 0.D0
!
!... Terms k_i that are computed in the Runge-Kutta algorithm
!
      ALLOCATE(k1(ncint,neq))
      ALLOCATE(k2(ncint,neq))
      ALLOCATE(k3(ncint,neq))
      k1 = 0.D0
      k2 = 0.D0
      k3 = 0.D0
!
! ... coefficients of RK3
!
      rk_c    = 0.D0
      rk_c(1) = 0.D0
      rk_c(2) = 1.D0
      rk_c(3) = 0.5D0
!
      rk_b    = 0.D0
      rk_b(1) = 1.D0/6.D0
      rk_b(2) = 1.D0/6.D0
      rk_b(3) = 4.D0/6.D0
!
      rk_a      = 0.D0
      rk_a(2,1) = 1.D0
      rk_a(3,1) = 0.25D0
      rk_a(3,2) = 0.25D0
!
! ... Courant number
!
      ALLOCATE(courant(ncint))
      courant = 0.D0
!
!... Matrices of the linear systems (for velocity components and enthalpies)
!
      ALLOCATE(au(nsolid+1,nsolid+1))
      au = 0.D0
      IF (job_type == JOB_TYPE_3D) THEN
         ALLOCATE(av(nsolid+1,nsolid+1))
         av = 0.D0
      END IF
      ALLOCATE(aw(nsolid+1,nsolid+1))
      aw = 0.D0
      ALLOCATE(at(nsolid+1,nsolid+1))
      at = 0.D0
!
!... Star variables, that are intermediate solutions of the Runge-Kutta algorithm
!
      ALLOCATE(rgpstar(ncdom))
      rgpstar = 0.D0
      ALLOCATE(rlkstar(ncdom,nsolid))
      rlkstar = 0.D0
      ALLOCATE(ugstar(ncdom))
      ugstar = 0.D0
      ALLOCATE(usstar(ncdom,nsolid))
      usstar = 0.D0
      IF (job_type == JOB_TYPE_3D) THEN
         ALLOCATE(vgstar(ncdom))
         ALLOCATE(vsstar(ncdom,nsolid))
         vgstar = 0.D0
         vsstar = 0.D0
      END IF
      ALLOCATE(wgstar(ncdom))
      wgstar = 0.D0
      ALLOCATE(wsstar(ncdom,nsolid))
      wsstar = 0.D0
      ALLOCATE(tgstar(ncdom))
      tgstar = 0.D0
      ALLOCATE(siegstar(ncdom))
      siegstar = 0.D0
      ALLOCATE(rhogstar(ncint))
      rhogstar = 0.D0
      ALLOCATE(epstar(ncdom))
      epstar = 0.D0
      ALLOCATE(pstar(ncdom))
      pstar = 0.D0
      ALLOCATE(tsstar(ncdom,nsolid))
      tsstar = 0.D0
      ALLOCATE(siesstar(ncdom,nsolid))
      siesstar = 0.D0

      ALLOCATE(rhon(ncint))
      rhon = 0.D0
      ALLOCATE(ugn(ncint))
      ugn = 0.D0
      ALLOCATE(usn(ncint,nsolid))
      usn = 0.D0
      IF (job_type == JOB_TYPE_3D) THEN
         ALLOCATE(vgn(ncint))
         ALLOCATE(vsn(ncint,nsolid))
         vgn = 0.D0
         vsn = 0.D0
      END IF
      ALLOCATE(wgn(ncint))
      wgn = 0.D0
      ALLOCATE(wsn(ncint,nsolid))
      wsn = 0.D0
      ALLOCATE(tgn(ncint))
      tgn = 0.D0

      ALLOCATE( rgfe( ncdom ), rgft( ncdom ))
      ALLOCATE( rsfe( ncdom, nsolid ), rsft( ncdom, nsolid ))
      rgfe = 0.0D0
      rgft = 0.0D0
      rsfe = 0.0D0
      rsft = 0.0D0
      IF (job_type == JOB_TYPE_3D) THEN
        ALLOCATE( rgfn( ncdom ))
        ALLOCATE( rsfn( ncdom, nsolid ))
        rgfn = 0.0D0
        rsfn = 0.0D0
      END IF
!
! ... Allocate and initialize gas convective fluxes
!
      ALLOCATE( ugfe(ncdom), ugft(ncdom) )
      ALLOCATE( wgfe(ncdom), wgft(ncdom) )
      ugfe = 0.0D0; ugft = 0.0D0
      wgfe = 0.0D0; wgft = 0.0D0

      IF (job_type == JOB_TYPE_3D) THEN
        ALLOCATE( ugfn(ncdom) )
        ALLOCATE( vgfe(ncdom), vgfn(ncdom), vgft(ncdom) )
        ALLOCATE( wgfn(ncdom) )
        ugfn = 0.0D0
        vgfe = 0.0D0; vgfn = 0.0D0; vgft = 0.0D0
        wgfn = 0.0D0
      END IF
!
! ... Allocate and initialize particles convective fluxes
!
      ALLOCATE( usfe(ncdom,nsolid), usft(ncdom,nsolid) )
      ALLOCATE( wsfe(ncdom,nsolid), wsft(ncdom,nsolid) )

      usfe = 0.0D0; usft = 0.0D0
      wsfe = 0.0D0; wsft = 0.0D0

      IF (job_type == JOB_TYPE_3D) THEN
        ALLOCATE( usfn(ncdom,nsolid) )
        ALLOCATE( vsfe(ncdom,nsolid), vsfn(ncdom,nsolid), vsft(ncdom,nsolid) )
        ALLOCATE( wsfn(ncdom,nsolid) )
        usfn = 0.0D0
        vsfe = 0.0D0; vsfn = 0.0D0; vsft = 0.0D0
        wsfn = 0.0D0
      END IF

      ALLOCATE(egfe(ncdom), egft(ncdom))
      ALLOCATE(esfe(ncdom,nsolid), esft(ncdom,nsolid))
!
      ALLOCATE(hgfe(ncdom), hgft(ncdom))
      ALLOCATE(hsfe(ncdom,nsolid), hsft(ncdom,nsolid))
!
      egfe = 0.0D0;  hgfe = 0.0D0
      egft = 0.0D0;  hgft = 0.0D0
      esfe = 0.0D0;  hsfe = 0.0D0
      esft = 0.0D0;  hsft = 0.0D0
!
      IF (job_type == JOB_TYPE_3D) THEN
        ALLOCATE(egfn(ncdom))
        ALLOCATE(esfn(ncdom,nsolid))
!
        ALLOCATE(hgfn(ncdom))
        ALLOCATE(hsfn(ncdom,nsolid))
!
        egfn = 0.0D0;  hgfn = 0.0D0
        esfn = 0.0D0;  hsfn = 0.0D0
      END IF
!
      RETURN

      END SUBROUTINE allocate_fields
!----------------------------------------------------------------------
      SUBROUTINE deallocate_fields

      USE iterative_solver, ONLY: rgfe, rgfn, rgft, rsfe, rsfn, rsft
      USE tilde_momentum, ONLY: ugfe, ugft, ugfn, vgfe, vgft, vgfn, wgfe, wgft, wgfn
      USE tilde_momentum, ONLY: usfe, usft, usfn, vsfe, vsft, vsfn, wsfe, wsft, wsfn
      USE tilde_energy, ONLY: egfe, egft, egfn
      USE tilde_energy, ONLY: hgfe, hgft, hgfn
      USE tilde_energy, ONLY: esfe, esft, esfn
      USE tilde_energy, ONLY: hsfe, hsft, hsfn
      USE gas_solid_viscosity, ONLY: gvisx, gvisy, gvisz, pvisx, pvisy, pvisz

      IMPLICIT NONE


      DEALLOCATE(kpgv)
      DEALLOCATE(hv)
      DEALLOCATE(appu)
      IF (job_type == JOB_TYPE_3D) THEN
         DEALLOCATE(appv)
      END IF
      DEALLOCATE(appw)

      DEALLOCATE(gvisx)
      DEALLOCATE(gvisz)
      DEALLOCATE(pvisx)
      DEALLOCATE(pvisz)
      IF ( job_type == JOB_TYPE_3D ) THEN
        DEALLOCATE( gvisy )
        DEALLOCATE( pvisy )
      END IF

      DEALLOCATE(rhs3)
      DEALLOCATE(k1)
      DEALLOCATE(k2)
      DEALLOCATE(k3)
      DEALLOCATE(courant)

      DEALLOCATE(au)
      IF (job_type == JOB_TYPE_3D) THEN
         DEALLOCATE(av)
      END IF
      DEALLOCATE(aw)
      DEALLOCATE(at)

      DEALLOCATE(rgpstar)
      DEALLOCATE(rlkstar)
      DEALLOCATE(ugstar)
      DEALLOCATE(usstar)
      IF (job_type == JOB_TYPE_3D) THEN
         DEALLOCATE(vgstar)
         DEALLOCATE(vsstar)
      END IF
      DEALLOCATE(wgstar)
      DEALLOCATE(wsstar)
      DEALLOCATE(tgstar)
      DEALLOCATE(siegstar)
      DEALLOCATE(rhogstar)
      DEALLOCATE(epstar)
      DEALLOCATE(pstar)
      DEALLOCATE(tsstar)
      DEALLOCATE(siesstar)

      DEALLOCATE(rhon)
      DEALLOCATE(ugn)
      DEALLOCATE(usn)
      IF (job_type == JOB_TYPE_3D) THEN
         DEALLOCATE(vgn)
         DEALLOCATE(vsn)
      END IF
      DEALLOCATE(wgn)
      DEALLOCATE(wsn)
      DEALLOCATE(tgn)

      DEALLOCATE(rgfe, rgft)
      DEALLOCATE(rsfe, rsft)
      IF (job_type == JOB_TYPE_3D) THEN
        DEALLOCATE(rgfn)
        DEALLOCATE(rsfn)
      END IF
      DEALLOCATE(ugfe, ugft)
      DEALLOCATE(wgfe, wgft)
      IF (job_type == JOB_TYPE_3D) THEN
        DEALLOCATE(ugfn)
        DEALLOCATE(vgfe, vgfn, vgft)
        DEALLOCATE(wgfn)
      END IF
      DEALLOCATE(usfe, usft)
      DEALLOCATE(wsfe, wsft)
      IF (job_type == JOB_TYPE_3D) THEN
        DEALLOCATE(usfn)
        DEALLOCATE(vsfe, vsfn, vsft)
        DEALLOCATE(wsfn)
      END IF
      DEALLOCATE(egfe, egft)
      DEALLOCATE(esfe, esft)
      DEALLOCATE(hgfe, hgft)
      DEALLOCATE(hsfe, hsft)
      IF (job_type == JOB_TYPE_3D) THEN
        DEALLOCATE(egfn)
        DEALLOCATE(esfn)
        DEALLOCATE(hgfn)
        DEALLOCATE(hsfn)
      END IF

      RETURN

      END SUBROUTINE deallocate_fields
!----------------------------------------------------------------------
      SUBROUTINE initialize_fieldn

      USE dimensions, ONLY: ngas, nsolid
      USE domain_mapping, ONLY: data_exchange, ncint
      USE eos_gas, ONLY: ygc
      USE gas_components, ONLY: rgpgcn
      USE gas_solid_density, ONLY: rgp, rgpn, rlk, rlkn, rog
      USE gas_solid_temperature, ONLY: sieg, siegn, sies, siesn, tg, ts
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      USE gas_solid_viscosity, ONLY: viscg, viscs
      USE gas_solid_viscosity, ONLY: gas_viscosity, part_viscosity
      USE gas_solid_viscosity, ONLY: gvisx, gvisy, gvisz, pvisx, pvisy, pvisz
      USE iterative_solver, ONLY: rgfe, rgfn, rgft, rsfe, rsfn, rsft
      USE pressure_epsilon, ONLY: p, ep, pn, epn
      USE set_indexes, ONLY: first_subscr, imjk, ijmk, ijkm
      USE tilde_momentum, ONLY: ugfe, ugft, ugfn, vgfe, vgft, vgfn, wgfe, wgft, wgfn
      USE tilde_momentum, ONLY: usfe, usft, usfn, vsfe, vsft, vsfn, wsfe, wsft, wsfn
      USE tilde_energy, ONLY: egfe, egft, egfn
      USE tilde_energy, ONLY: hgfe, hgft, hgfn
      USE tilde_energy, ONLY: esfe, esft, esfn
      USE tilde_energy, ONLY: hsfe, hsft, hsfn
!
      IMPLICIT NONE

      INTEGER :: ijk, is, ig
      REAL*8 :: dugs, dvgs, dwgs
      REAL*8, DIMENSION(:), ALLOCATABLE :: nul

      ALLOCATE( nul( ( ( nsolid + 1 )**2 + ( nsolid + 1 ) ) / 2 ) )
      nul = 0.D0

      kpgv = 0.D0
      hv = 0.D0
      appu = 0.0D0
      appw = 0.0D0
      IF (job_type == JOB_TYPE_3D) THEN
        appv = 0.0D0
      END IF
      gvisx = 0.D0; gvisz = 0.D0
      pvisx = 0.D0; pvisz = 0.D0
      IF ( job_type == JOB_TYPE_3D ) THEN
        gvisy = 0.D0
        pvisy = 0.D0
      END IF

      rgfe = 0.0D0
      rgft = 0.0D0
      rsfe = 0.0D0
      rsft = 0.0D0
      IF (job_type == JOB_TYPE_3D) THEN
        rgfn = 0.0D0
        rsfn = 0.0D0
      END IF
      ugfe = 0.0D0; ugft = 0.0D0
      wgfe = 0.0D0; wgft = 0.0D0
      IF (job_type == JOB_TYPE_3D) THEN
        ugfn = 0.0D0
        vgfe = 0.0D0; vgfn = 0.0D0; vgft = 0.0D0
        wgfn = 0.0D0
      END IF
      usfe = 0.0D0; usft = 0.0D0
      wsfe = 0.0D0; wsft = 0.0D0
      IF (job_type == JOB_TYPE_3D) THEN
        usfn = 0.0D0
        vsfe = 0.0D0; vsfn = 0.0D0; vsft = 0.0D0
        wsfn = 0.0D0
      END IF
      egfe = 0.0D0;  hgfe = 0.0D0
      egft = 0.0D0;  hgft = 0.0D0
      esfe = 0.0D0;  hsfe = 0.0D0
      esft = 0.0D0;  hsft = 0.0D0
      IF (job_type == JOB_TYPE_3D) THEN
        egfn = 0.0D0;  hgfn = 0.0D0
        esfn = 0.0D0;  hsfn = 0.0D0
      END IF
!
! ... Data_exchange of primary fields
!
      CALL data_exchange(ug)   
      CALL data_exchange(us) 
      IF (job_type == JOB_TYPE_3D) THEN
        CALL data_exchange(vg) 
        CALL data_exchange(vs) 
      END IF
      CALL data_exchange(wg)  
      CALL data_exchange(ws) 

      CALL data_exchange(p) 
      CALL data_exchange(ep)
      CALL data_exchange(rgp) 
      CALL data_exchange(rlk)

      CALL data_exchange(sieg)
      CALL data_exchange(sies)
      CALL data_exchange(tg) 
      CALL data_exchange(ts)
      !
      ! ... star variables are initialized with solution at previous time step
      !
      CALL init_star

      DO ijk = 1, ncint
        tgn(ijk)   = tg(ijk)
        rhon(ijk)  = rog(ijk)
        epn(ijk)   = ep(ijk)
        pn(ijk)    = p(ijk)
        rgpn(ijk)  = rgp(ijk)
        siegn(ijk) = sieg(ijk)
        ugn(ijk) = ug(ijk)
        IF (job_type .EQ. JOB_TYPE_3D) vgn(ijk) = vg(ijk)
        wgn(ijk) = wg(ijk)
        DO is = 1, nsolid
          rlkn(ijk,is)  = rlk(ijk,is)
          siesn(ijk,is) = sies(ijk,is)
          usn(ijk,is) = us(ijk,is)
          IF (job_type .EQ. JOB_TYPE_3D) vsn(ijk,is) = vs(ijk,is)
          wsn(ijk,is) = ws(ijk,is)
        END DO

        DO ig = 1, ngas
          rgpgcn(ijk,ig) = rgpn(ijk) * ygc(ijk,ig)
        END DO

        CALL compute_physical_coeff(ijk)

      END DO
!
! ... Calculate gas viscous stress tensor
!
      IF ( gas_viscosity ) CALL viscg        
!
! ... Calculate particles viscous stress tensor
!
      IF ( part_viscosity ) CALL viscs
!
      CALL data_exchange(pn)
      CALL data_exchange(epn)
      CALL data_exchange(rgpn)
      CALL data_exchange(rlkn)
      CALL data_exchange(appu)
      IF (job_type .EQ. JOB_TYPE_3D) CALL data_exchange(appv)
      CALL data_exchange(appw)

      RETURN

      END SUBROUTINE initialize_fieldn
!----------------------------------------------------------------------
      SUBROUTINE compute_rhs(rhs,ijk)

      USE atmospheric_conditions, ONLY: gravx, gravy, gravz
      USE dimensions, ONLY: nsolid
      USE domain_mapping, ONLY: meshinds
      USE eos_gas, ONLY: cg
      USE gas_solid_viscosity, ONLY: gvisx, gvisy, gvisz, pvisx, pvisy, pvisz
      USE gas_solid_density, ONLY: rlk, rlkn
      USE grid, ONLY: indx, indy, indz, dx, dy, dz, inr, inrb
      USE iterative_solver, ONLY: rgfe, rgfn, rgft, rsfe, rsfn, rsft
      USE particles_constants, ONLY: inrl
      USE pressure_epsilon, ONLY: pn, pmodel
      USE set_indexes, ONLY: imjk, ijmk, ijkm
      USE set_indexes, ONLY: ijke, ijkn, ijkt, ijkw, ijks, ijkb
      USE set_indexes, ONLY: subscr
      USE specific_heat_module, ONLY: ck
      USE tilde_energy, ONLY: egfe, egfn, egft, esfe, esfn, esft
      USE tilde_momentum, ONLY: ugfe, ugfn, ugft, vgfe, vgfn, vgft, wgfe, wgfn, wgft
      USE tilde_momentum, ONLY: usfe, usfn, usft, vsfe, vsfn, vsft, wsfe, wsfn, wsft
      USE time_parameters, ONLY: dt

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ijk
      REAL*8, DIMENSION(:), INTENT(OUT) :: rhs
      INTEGER :: is, i, j, k, imesh, eq
      REAL*8 :: fx, fy, fz
      REAL*8 :: dxp, indxp, dyp, indyp, dzp, indzp
      REAL*8 :: epne, epnn, epnt
      REAL*8 :: epsne, epsnn, epsnt
      REAL*8 :: rhstemp, tmp, dpxyz
      REAL*8 :: ugc, vgc, wgc
      REAL*8 :: indxc, indyc, indzc

      fx = 0.D0
      fy = 0.D0
      fz = 0.D0
      rhs = 0.D0

      CALL meshinds(ijk,imesh,i,j,k)

      dxp   = dx(i) + dx(i+1)
      dzp   = dz(k) + dz(k+1)
      indxp = 1.D0 / dxp
      indzp = 1.D0 / dzp

      IF (job_type .EQ. JOB_TYPE_3D) THEN
        dyp   = dy(j) + dy(j+1)
        indyp = 1.D0 / dyp
      END IF

      CALL subscr(ijk)
!
!... RHS of the gas mass equation
!
      eq = 1
      rhstemp = 0.D0
      fx = (rgfe(ijk)-rgfe(imjk))*indx(i)*inr(i)
      fy = 0.D0
      IF (job_type == JOB_TYPE_3D) THEN
         fy = (rgfn(ijk)-rgfn(ijmk))*indy(j)
      END IF
      fz = (rgft(ijk)-rgft(ijkm))*indz(k)
      tmp = fx + fy + fz
      rhstemp = rhstemp - tmp !(fx+fy+fz)
      rhs(eq) = rhstemp
      eq = eq + 1
!
!... RHS of the particles mass equations
!
      DO is = 1, nsolid
         rhstemp = 0.D0
         fx = (rsfe(ijk,is)-rsfe(imjk,is))*indx(i)*inr(i)
         fy = 0.D0
         IF (job_type == JOB_TYPE_3D) THEN
            fy = (rsfn(ijk,is)-rsfn(ijmk,is))*indy(j)
         END IF
         fz = (rsft(ijk,is)-rsft(ijkm,is))*indz(k)
         tmp = fx + fy + fz
         rhstemp = rhstemp - tmp !(fx+fy+fz)
         rhs(eq) = rhstemp
         eq = eq + 1 
      END DO
!
      IF ( pmodel == 1 ) THEN
        epne = (dx(i)*epstar(ijke) + dx(i+1)*epstar(ijk)) * indxp
        epnn = (dy(j)*epstar(ijkn) + dy(j+1)*epstar(ijk)) *indyp
        epnt = (dz(k)*epstar(ijkt) + dz(k+1)*epstar(ijk)) * indzp
       ELSE IF ( pmodel == 2 ) THEN
        epne = 1.D0
        epnn = 1.D0
        epnt = 1.D0
       END IF
!
!... RHS of the gas momentum-x equation
!
      fx = (ugfe(ijk)-ugfe(imjk))*2.D0*indxp*inrb(i)
      fy = 0.D0
      IF (job_type == JOB_TYPE_3D) THEN
         fy = (ugfn(ijk)-ugfn(ijmk))*indy(j)
      END IF
      fz = (ugft(ijk)-ugft(ijkm))*indz(k)
      tmp = fx + fy + fz
      rhstemp = 0.D0
      rhstemp = rhstemp - tmp !(fx+fy+fz)
      rhstemp = rhstemp - epne*(pstar(ijke)-pstar(ijk))*2.D0*indxp
      rhstemp = rhstemp + (dx(i)*rgpstar(ijke)+dx(i+1)*rgpstar(ijk))*indxp*gravx
      rhstemp = rhstemp + gvisx(ijk)
      rhs(eq) = rhstemp
      eq = eq + 1
!
!... RHS of particles momentum-x equations
!
      DO is = 1, nsolid
         fx = (usfe(ijk,is)-usfe(imjk,is))*2.D0*indxp*inrb(i)
         fy = 0.D0
         IF (job_type == JOB_TYPE_3D) THEN
            fy = (usfn(ijk,is)-usfn(ijmk,is))*indy(j)
         END IF
         fz = (usft(ijk,is)-usft(ijkm,is))*indz(k)
         tmp = fx + fy + fz
         epsne = (dx(i)*rlkstar(ijke,is) + dx(i+1)*rlkstar(ijk,is))*indxp*inrl(is)
         epsne = epsne * (2-pmodel)
         rhstemp = 0.D0
         rhstemp = rhstemp - tmp !(fx+fy+fz)
         rhstemp = rhstemp - epsne*(pstar(ijke)-pstar(ijk))*2.D0*indxp
         rhstemp = rhstemp + (dx(i)*rlkstar(ijke,is)+dx(i+1)*rlkstar(ijk,is))*indxp*gravx
         rhstemp = rhstemp + pvisx(ijk,is)
         rhs(eq) = rhstemp
         eq = eq + 1
      END DO
!
!... RHS of gas momentum-y equation
!
      IF (job_type == JOB_TYPE_3D) THEN
         fx = (vgfe(ijk)-vgfe(imjk))*indx(i)
         fy = (vgfn(ijk)-vgfn(ijmk))*2.D0*indyp
         fz = (vgft(ijk)-vgft(ijkm))*indz(k)
         tmp = fx + fy + fz
         rhstemp = 0.D0
         rhstemp = rhstemp - tmp !(fx+fy+fz)
         rhstemp = rhstemp - epnn*(pstar(ijkn)-pstar(ijk))*2.D0*indyp
         rhstemp = rhstemp + (dy(j)*rgpstar(ijkn)+dy(j+1)*rgpstar(ijk))*indyp*gravy
         rhstemp = rhstemp + gvisy(ijk)
         rhs(eq) = rhstemp
         eq = eq + 1  
!
!... RHS of particles momentum-y equations
!
         DO is = 1, nsolid
            fx = (vsfe(ijk,is)-vsfe(imjk,is))*indx(i)
            fy = (vsfn(ijk,is)-vsfn(ijmk,is))*2.D0*indyp
            fz = (vsft(ijk,is)-vsft(ijkm,is))*indz(k)
            tmp = fx + fy + fz
            epsnn = (dy(j)*rlkstar(ijkn,is) + dy(j+1)*rlkstar(ijk,is))*indyp*inrl(is)
            epsnn = epsnn * (2-pmodel)
            rhstemp = 0.D0
            rhstemp = rhstemp - tmp !(fx+fy+fz)
            rhstemp = rhstemp - epsnn*(pstar(ijkn)-pstar(ijk))*2.D0*indyp
            rhstemp = rhstemp + (dy(j)*rlkstar(ijkn,is)+dy(j+1)*rlkstar(ijk,is))*indyp*gravy
            rhstemp = rhstemp + pvisy(ijk,is)
            rhs(eq) = rhstemp
            eq = eq + 1
         END DO
      END IF
!
!... RHS of gas momentum-z equation
!
      fx = (wgfe(ijk)-wgfe(imjk))*indx(i)*inr(i)
      fy = 0.D0
      IF (job_type == JOB_TYPE_3D) THEN
         fy = (wgfn(ijk)-wgfn(ijmk))*indy(j)
      END IF
      fz = (wgft(ijk)-wgft(ijkm))*2.D0*indzp
      tmp = fx + fy + fz
      rhstemp = 0.D0
      rhstemp = rhstemp - tmp !(fx+fy+fz)
      rhstemp = rhstemp - epnt*(pstar(ijkt)-pstar(ijk))*2.D0*indzp
      rhstemp = rhstemp + (dz(k)*rgpstar(ijkt)+dz(k+1)*rgpstar(ijk))*indzp*gravz
      rhstemp = rhstemp + gvisz(ijk)
      rhs(eq) = rhstemp
      eq = eq + 1
!
!... RHS of particles momentum-z equations
!
      DO is = 1, nsolid
         fx = (wsfe(ijk,is)-wsfe(imjk,is))*indx(i)*inr(i)
         fy = 0.D0
         IF (job_type == JOB_TYPE_3D) THEN
            fy = (wsfn(ijk,is)-wsfn(ijmk,is))*indy(j)
         END IF
         fz = (wsft(ijk,is)-wsft(ijkm,is))*2.D0*indzp
         tmp = fx + fy + fz
         epsnt = (dz(k)*rlkstar(ijkt,is) + dz(k+1)*rlkstar(ijk,is))*indzp*inrl(is)
         epsnt = epsnt * (2 - pmodel)
         rhstemp = 0.D0
         rhstemp = rhstemp - tmp !(fx+fy+fz)
         rhstemp = rhstemp - epsnt*(pstar(ijkt)-pstar(ijk))*2.D0*indzp
         rhstemp = rhstemp + (dz(k)*rlkstar(ijkt,is)+dz(k+1)*rlkstar(ijk,is))*indzp*gravz
         rhstemp = rhstemp + pvisz(ijk,is)
         rhs(eq) = rhstemp
         eq = eq + 1

      END DO
!
!... RHS of gas enthalpy equation
!
      indxc = 1.D0/(dx(i)+(dx(i+1)+dx(i-1))*0.5D0)
      indyc = 0.D0
      indzc = 1.D0/(dz(k)+(dz(k+1)+dz(k-1))*0.5D0)
      ugc = 0.5D0 * (ugstar(ijk) + ugstar(imjk)) 
      vgc = 0.D0
      wgc = 0.5D0 * (wgstar(ijk) + wgstar(ijkm))
      IF (job_type == JOB_TYPE_3D) THEN
         indyc = 1.D0 / (dy(j)+(dy(j+1)+dy(j-1))*0.5D0)
         vgc = 0.5D0 * (vgstar(ijk) + vgstar(ijmk))
      END IF
      fx = (egfe(ijk)-egfe(imjk))*indx(i)*inr(i)
      fy = 0.D0
      IF (job_type == JOB_TYPE_3D) THEN
         fy = (egfn(ijk)-egfn(ijmk))*indy(j)
      END IF
      fz = (egft(ijk)-egft(ijkm))*indz(k)
      tmp = fx + fy + fz
      dpxyz= 0.D0
      dpxyz = dpxyz + indxc*ugc*(pstar(ijke)-pstar(ijkw))
      dpxyz = dpxyz + indyc*vgc*(pstar(ijkn)-pstar(ijks))
      dpxyz = dpxyz + indzc*wgc*(pstar(ijkt)-pstar(ijkb))
      rhstemp = 0.D0
      rhstemp = rhstemp - tmp !(fx+fy+fz)
      tmp = 1.D0/dt
      rhstemp = rhstemp + epstar(ijk)*(pstar(ijk)-pn(ijk))*tmp
      rhstemp = rhstemp + epstar(ijk)*dpxyz
      rhs(eq) = rhstemp
      eq = eq + 1 
!
!... RHS of particles enthalpy equations
!
      DO is = 1, nsolid
         fx = (esfe(ijk,is)-esfe(imjk,is))*indx(i)*inr(i)
         fy = 0.D0
         IF (job_type == JOB_TYPE_3D) THEN
            fy = (esfn(ijk,is)-esfn(ijmk,is))*indy(j)
         END IF
         fz = (esft(ijk,is)-esft(ijkm,is))*indz(k)
         tmp = fx + fy + fz
         rhstemp = 0.D0
         rhstemp = rhstemp - tmp !(fx+fy+fz)
         rhs(eq) = rhstemp
         eq = eq + 1 
      END DO

      RETURN

      END SUBROUTINE compute_rhs
!----------------------------------------------------------------------
      SUBROUTINE compute_linear_systems(ijk)
!
!... compute the matrices au, av, aw, at of the linear systems 
!
      USE dimensions, ONLY: nsolid, ngas
      USE domain_mapping, ONLY: meshinds
      USE eos_gas, ONLY: cg
      USE grid, ONLY: dx, dy, dz
      USE set_indexes, ONLY: ijke, ijkt, ijkn, subscr
      USE specific_heat_module, ONLY: ck
      USE time_parameters, ONLY: dt, alpha

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ijk
      INTEGER :: is, is1, ig, i, j, k, imesh, nphase, l, ll, ls, ls1, eq, eqg
      REAL*8 :: dxi, dxip1, dyj, dyjp1, dzk, dzkp1
      REAL*8 :: dxp, dyp, dzp, indxp, indyp, indzp

      CALL meshinds(ijk,imesh,i,j,k)

      CALL subscr(ijk)

      dxi   = dx(i)
      dxip1 = dx(i+1)
      dyj   = dy(j)
      dyjp1 = dy(j+1)
      dzk   = dz(k)
      dzkp1 = dz(k+1)
!
      dxp = dxi + dxip1
      dyp = dyj + dyjp1
      dzp = dzk + dzkp1
      indxp = 1.D0/dxp
      indyp = 1.D0/dyp
      indzp = 1.D0/dzp

      nphase = nsolid+1
!
!... build the matrix for the momentum-x equation. Remember that only the drag terms are treated 
!    implicitly in the Runge-Kutta algorithm
!
      au(1,1) = (dxi*rgpstar(ijke)+dxip1*rgpstar(ijk))*indxp
      au(1,1) = au(1,1) + alpha*(dxi*appu(ijke,1)+dxip1*appu(ijk,1))*indxp
      IF (job_type == JOB_TYPE_3D) THEN
         av(1,1) = (dyj*rgpstar(ijkn)+dyjp1*rgpstar(ijk))*indyp
         av(1,1) = av(1,1) + alpha*(dyj*appv(ijkn,1)+dyjp1*appv(ijk,1))*indyp
      END IF
      aw(1,1) = (dzk*rgpstar(ijkt)+dzkp1*rgpstar(ijk))*indzp
      aw(1,1) = aw(1,1) + alpha*(dzk*appw(ijkt,1)+dzkp1*appw(ijk,1))*indzp

      DO l=2,nphase
         ls1=l*(l-1)/2
         DO ll=1,l
            ls=ls1+ll
            au(l,ll) = alpha*(dxi*appu(ijke,ls)+dxip1*appu(ijk,ls))*indxp
            au(ll,l) = au(l,ll)
            aw(l,ll) = alpha*(dzk*appw(ijkt,ls)+dzkp1*appw(ijk,ls))*indzp
            aw(ll,l) = aw(l,ll)
            IF (job_type == JOB_TYPE_3D) THEN
               av(l,ll) = alpha*(dyj*appv(ijkn,ls)+dyjp1*appv(ijk,ls))*indyp
               av(ll,l) = av(l,ll)
            END IF
         END DO
 
         au(l,l) = au(l,l)+(dxi*rlkstar(ijke,l-1)+dxip1*rlkstar(ijk,l-1))*indxp
         aw(l,l) = aw(l,l)+(dzk*rlkstar(ijkt,l-1)+dzkp1*rlkstar(ijk,l-1))*indzp
         IF (job_type == JOB_TYPE_3D) THEN
            av(l,l)=av(l,l)+(dyj*rlkstar(ijkn,l-1)+dyjp1*rlkstar(ijk,l-1))*indyp
         END IF
      END DO
!
!... compute the matrix related to the enthalpy system of equations
!
      at(1,1) = rgpstar(ijk)
      DO is = 1,nsolid
         is1 = is+1
         at(1,1)     = at(1,1)         + alpha*dt*hv(ijk,is)/cg(ijk)
         at(is1,1)   =                 - alpha*dt*hv(ijk,is)/cg(ijk)
         at(1,is1)   =                 - alpha*dt*hv(ijk,is)/ck(is,ijk)
         at(is1,is1) = rlkstar(ijk,is) + alpha*dt*hv(ijk,is)/ck(is,ijk)
      END DO
!
! ... if alpha < 1, then we need to add the Crank-Nicolson term (1-alpha)
!     related to the drag and the heat exchange term to the right-hans-side
!     of the linear system
!
      IF (alpha < 1.D0) THEN

        eq = nphase + 1
        eqg = eq

        ! gas momentum-x
        rhs3(ijk,eqg) = rhs3(ijk,eqg) - (1.D0-alpha)*(dxi*appu(ijke,1)+dxip1*appu(ijk,1))*indxp*ugstar(ijk) 
        eq = eq + 1       
        ! particles momentum-x
        DO l = 2,nphase
          ls1 = l*(l-1)/2
          DO ll=1,l
            ls = ls1+ll
            IF (ll .EQ. 1) THEN
              rhs3(ijk,eqg) = rhs3(ijk,eqg) - (1.D0-alpha)*(dxi*appu(ijke,ls)+dxip1*appu(ijk,ls))*indxp*usstar(ijk,l-1)
              rhs3(ijk,eqg+l-1) = rhs3(ijk,eqg+l-1)   - (1.D0-alpha)*(dxi*appu(ijke,ls)+dxip1*appu(ijk,ls))*indxp*ugstar(ijk)
            ELSE IF (ll .EQ. l) THEN
              rhs3(ijk,eqg+l-1) = rhs3(ijk,eqg+l-1) - (1.D0-alpha)*(dxi*appu(ijke,ls)+dxip1*appu(ijk,ls))*indxp*usstar(ijk,l-1)
            ELSE
              rhs3(ijk,eqg+l-1) = rhs3(ijk,eqg+l-1) - (1.D0-alpha)*(dxi*appu(ijke,ls)+dxip1*appu(ijk,ls))*indxp*usstar(ijk,ll-1)
              rhs3(ijk,eqg+ll-1) = rhs3(ijk,eqg+ll-1) - (1.D0-alpha)*(dxi*appu(ijk,ls)+dxip1*appu(ijk,ls))*indxp*usstar(ijk,l-1)
            END IF
          END DO 
          eq = eq + 1
        END DO

        IF (job_type .EQ. JOB_TYPE_3D) THEN
          eqg = eq
          ! gas momentum-y
          rhs3(ijk,eqg) = rhs3(ijk,eqg) - (1.D0-alpha)*(dyj*appv(ijkn,1)+dyjp1*appv(ijk,1))*indyp*vgstar(ijk) 
          eq = eq + 1
          ! particles momentum-y
          DO l = 2,nphase
            ls1 = l*(l-1)/2
            DO ll=1,l
              ls = ls1+ll
              IF (ll .EQ. 1) THEN
                rhs3(ijk,eqg) = rhs3(ijk,eqg) - (1.D0-alpha)*(dyj*appv(ijkn,ls)+dyjp1*appv(ijk,ls))*indyp*vsstar(ijk,l-1)
                rhs3(ijk,eqg+l-1) = rhs3(ijk,eqg+l-1)   - (1.D0-alpha)*(dyj*appv(ijkn,ls)+dyjp1*appv(ijk,ls))*indyp*vgstar(ijk)
              ELSE IF (ll .EQ. l) THEN
                rhs3(ijk,eqg+l-1) = rhs3(ijk,eqg+l-1) - (1.D0-alpha)*(dyj*appv(ijkn,ls)+dyjp1*appv(ijk,ls))*indyp*vsstar(ijk,l-1)
              ELSE
                rhs3(ijk,eqg+l-1) = rhs3(ijk,eqg+l-1) - (1.D0-alpha)*(dyj*appv(ijkn,ls)+dyjp1*appv(ijk,ls))*indyp*vsstar(ijk,ll-1)
                rhs3(ijk,eqg+ll-1) = rhs3(ijk,eqg+ll-1) - (1.D0-alpha)*(dyj*appv(ijk,ls)+dyjp1*appv(ijk,ls))*indyp*vsstar(ijk,l-1)
              END IF
            END DO
            eq = eq + 1
          END DO
        END IF

        eqg = eq
        ! gas momentum-z
        rhs3(ijk,eqg) = rhs3(ijk,eqg) - (1.D0-alpha)*(dzk*appw(ijkt,1)+dzkp1*appw(ijk,1))*indzp*wgstar(ijk) 
        eq = eq + 1
        ! particles momentum-z
        DO l = 2,nphase
          ls1 = l*(l-1)/2
          DO ll=1,l
            ls = ls1+ll
            IF (ll .EQ. 1) THEN
              rhs3(ijk,eqg) = rhs3(ijk,eqg) - (1.D0-alpha)*(dzk*appw(ijkt,ls)+dzkp1*appw(ijk,ls))*indzp*wsstar(ijk,l-1)
              rhs3(ijk,eqg+l-1) = rhs3(ijk,eqg+l-1) - (1.D0-alpha)*(dzk*appw(ijkt,ls)+dzkp1*appw(ijk,ls))*indzp*wgstar(ijk)
            ELSE IF (ll .EQ. l) THEN
              rhs3(ijk,eqg+l-1) = rhs3(ijk,eqg+l-1) - (1.D0-alpha)*(dzk*appw(ijkt,ls)+dzkp1*appw(ijk,ls))*indzp*wsstar(ijk,l-1)
            ELSE
              rhs3(ijk,eqg+l-1) = rhs3(ijk,eqg+l-1) - (1.D0-alpha)*(dzk*appw(ijkt,ls)+dzkp1*appw(ijk,ls))*indzp*wsstar(ijk,ll-1)
              rhs3(ijk,eqg+ll-1) = rhs3(ijk,eqg+ll-1) - (1.D0-alpha)*(dzk*appw(ijk,ls)+dzkp1*appw(ijk,ls))*indzp*wsstar(ijk,l-1)
            END IF
          END DO
          eq = eq + 1
        END DO

        ! gas and particles enthalpy
        eqg = eq
        eq = eq + 1
        DO is = 1,nsolid
          rhs3(ijk,eqg) = rhs3(ijk,eqg) + (1.D0-alpha)*dt*hv(ijk,is)*(tsstar(ijk,is)-tgstar(ijk))
          rhs3(ijk,eq)  = rhs3(ijk,eq)  + (1.D0-alpha)*dt*hv(ijk,is)*(tgstar(ijk)-tsstar(ijk,is))
          eq = eq + 1
        END DO

      END IF

      RETURN

      END SUBROUTINE compute_linear_systems
!----------------------------------------------------------------------
      SUBROUTINE explicit_iter
! 
! ... this subroutine compute the solution at time step n+1 from the solution
!     at time step n with Embedded 3-2 Runge-Kutta scheme.
!     Only drag terms and heat exchange terms are treated implicitly, whereas
!     all the other terms in the equations are treated explicitly
!
      USE dimensions
      USE domain_mapping, ONLY: data_exchange, ncint
      USE eos_gas, ONLY: cg, xgc, caloric_eosg, ygc, thermal_eosg
      USE eos_solid, ONLY: caloric_eosl
      USE gas_constants, ONLY: gmw, rgas, gas_type
      USE gas_solid_density, ONLY: rog, rgpn
      USE gas_solid_temperature, ONLY: siegn
      USE gas_solid_viscosity, ONLY: viscon, mug, kapg
      USE grid, ONLY: flag
      USE heat_transfer, ONLY: hvs
      USE particles_constants, ONLY: inrl, cps
      USE pressure_epsilon, ONLY: pn, epn
      USE set_indexes, ONLY: subscr, first_subscr, imjk, ijmk, ijkm
      USE specific_heat_module, ONLY: ck, cp, hcapg, hcaps
      USE time_parameters, ONLY: dt

      IMPLICIT NONE

      INTEGER :: ijk, is, eq, ig, info, num
      REAL*8 :: rtilde, epstot, epsstar, mg, tmp
      LOGICAL :: compute
      REAL*8 :: xgcl(max_ngas), rlktmp, cgas
      REAL*8 :: dugs, dvgs, dwgs
      REAL*8, DIMENSION(:), ALLOCATABLE :: nul

      ALLOCATE( nul( ( ( nsolid + 1 )**2 + ( nsolid + 1 ) ) / 2 ) )
      nul = 0.D0
      rlktmp = 0.D0
      !
      ! ... initialize rhs of ODE system
      !
      CALL init_rhs
      !
      ! ... first step of RK32 method
      !
      ! ... compute all fluxes for mass, momentum and enthalpy equations for gas and particles
      !
      CALL data_exchange(epstar)
      CALL data_exchange(pstar)
      CALL data_exchange(rgpstar)
      CALL data_exchange(rlkstar)

      CALL compute_fluxes
      !
      DO ijk = 1, ncint

        compute = BTEST(flag(ijk),0)
        IF (compute) THEN
          !
          ! ... evaluate the right hand side of the ODE using star variables
          !
          CALL compute_rhs(k1(ijk,:),ijk)
          !
          ! ... update star variables
          !
          CALL update_star(ijk,k1(ijk,:),k2(ijk,:),rk_a(2,:))
!
        END IF
      END DO
      !
      ! ... second step of RK32 scheme
      !
      ! ... compute all fluxes
      !
      CALL data_exchange(epstar)
      CALL data_exchange(pstar)
      CALL data_exchange(rgpstar)
      CALL data_exchange(rlkstar)

      CALL compute_fluxes
      !
      DO ijk = 1, ncint
        compute = BTEST(flag(ijk),0)
        IF (compute) THEN
          !
          ! ... evaluate the right hand side of the ODE using star variables
          !
          CALL compute_rhs(k2(ijk,:),ijk)
          !
          ! ... update star variables
          !
          CALL update_star(ijk,k1(ijk,:),k2(ijk,:),rk_a(3,:))

        END IF
      END DO
      !
      !
      ! ... third step of RK32 scheme
      !
      ! ... compute all fluxes
      !
      CALL data_exchange(epstar)
      CALL data_exchange(pstar)
      CALL data_exchange(rgpstar)
      CALL data_exchange(rlkstar)

      CALL compute_fluxes
      !
      DO ijk = 1, ncint
        compute = BTEST(flag(ijk),0)
        IF (compute) THEN
!
          CALL compute_rhs(k3(ijk,:),ijk)

        END IF 
      END DO
      !
      ! ... update right hand side of the final ODE system
      !
      rhs3 = rhs3 + rk_b(1)*dt*k1 &
                  + rk_b(2)*dt*k2 &
                  + rk_b(3)*dt*k3
      !
      ! ... solve ODEs
      !
      DO ijk = 1,ncint
        compute = BTEST(flag(ijk),0)
        IF (compute) THEN
          !
          ! ... mass equations can be solved directly to compute densities
          !
          eq = 1
          rgpstar(ijk) = rhs3(ijk,eq)
          eq = eq+1
          DO is = 1,nsolid
            rlktmp = rhs3(ijk,eq)

            if (rlktmp < 0.D0) then
              write(testunit,*) 'clipping in explicit iter &
                             rlk<0, ijk = ', ijk, ' rlktmp =', rlktmp
            end if
            rlkstar( ijk, is ) = MAX( 0.0d0, rlktmp )
            eq = eq+1
          END DO
        END IF
      END DO

      CALL data_exchange(rgpstar)
      CALL data_exchange(rlkstar)

      DO ijk = 1,ncint
        compute = BTEST(flag(ijk),0)

        IF (compute) CALL compute_physical_coeff(ijk)

      END DO

      CALL data_exchange(appu)
      CALL data_exchange(appw)

      num = 0

      DO ijk = 1,ncint
        compute = BTEST(flag(ijk),0)
        IF (compute) THEN

          CALL subscr(ijk)
          !
          ! ... Calculate particles viscous stress tensor
          !          
          ! ... build up matrices for momentum and enthalpy linear systems
          !
          CALL compute_linear_systems(ijk)
          !
          ! ... solve momentum and enthalpy linear systems
          !
          CALL solve_linear_systems(ijk)

          CALL caloric_eosg(cp(:,ijk), cg(ijk), tgstar(ijk), ygc(ijk,:), &
                            siegstar(ijk), ijk, info)
          IF (info > 0) then
             write(testunit,*) 'error in explicit iter'
             CALL write_star(ijk)
end if

          DO is=1, nsolid
            CALL caloric_eosl(tsstar(ijk,is),cps(is),ck(is,ijk),siesstar(ijk,is),ijk,info) 
          END DO
          num = num+info

          epstot  = 0.D0
          epsstar = 0.D0
          DO is = 1,nsolid
            epsstar = rlkstar(ijk,is)*inrl(is)
            epstot = epstot + epsstar
          END DO
          epstar(ijk) = 1.D0-epstot

          IF (epstar(ijk) > 1.D0) then
            write(testunit,*) 'ep > 1 in explicit iter in ijk =', ijk
            call write_star(ijk)
          end if
          IF (epstar(ijk) < 0.D0) then
            write(testunit,*) 'ep < 0 in explicit iter in ijk =', ijk
            call write_star(ijk)
          end if

          rhogstar(ijk) = rgpstar(ijk)/epstar(ijk)

          mg = 0.D0
          DO ig = 1, ngas
            mg = mg + xgc(ijk,ig) * gmw(gas_type(ig))
          END DO
          pstar(ijk) = rhogstar(ijk)*tgstar(ijk)*rgas/mg
          ! 
          ! ... copy star-variables in the final values
          !
          CALL save_star(ijk)
          !
          ! ... compute errors (needed only for time step adaptation)
          !
        END IF
      END DO

      CALL parallel_sum_integer( num, 1 )
      IF ( num > 1 ) &
        CALL error('explicit_iter','Error in caloric equation of state', num)

      CALL compute_courant

      END SUBROUTINE explicit_iter
!----------------------------------------------------------------------
      SUBROUTINE compute_fluxes
!
! ... the routine computes all the fluxes using star variables (that are the
!     current solution in the RK step). Fluxes are compute for both gas and particles
!     mass, momentum-x-y-z and enthalpy equations
!
      USE dimensions, ONLY: nsolid
      USE domain_mapping, ONLY: ncint, data_exchange
      USE flux_limiters, ONLY: ctu
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      USE gas_solid_temperature, ONLY: sieg, sies, tg, ts
      USE gas_solid_density, ONLY: rgp, rlk, rog
      USE grid, ONLY: flag
      USE pressure_epsilon, ONLY: p, ep
      USE set_indexes, ONLY: subscr, ctu1_subscr, ctu2_subscr, ctu3_subscr
      USE tilde_energy, ONLY: egfe, egfn, egft, esfe, esfn, esft
      USE turbulence_model, ONLY: kapgt, iturb
      USE gas_solid_viscosity, ONLY: kapg
      USE tilde_momentum, ONLY: ugfe, ugft, ugfn, vgfe, vgft, vgfn, wgfe, wgft, wgfn
      USE tilde_momentum, ONLY: usfe, usft, usfn, vsfe, vsft, vsfn, wsfe, wsft, wsfn
!
      INTEGER :: ijk, is
      LOGICAL :: compute
!
! ... these copies can be avoided if we create new routines for star variables
!
      rgp = rgpstar
      rlk = rlkstar
      ug = ugstar
      us = usstar
      IF (job_type == JOB_TYPE_3D) THEN
         vg = vgstar
         vs = vsstar
      END IF
      wg = wgstar
      ws = wsstar
      sieg = siegstar 
      sies = siesstar
      tg = tgstar
      ts = tsstar
      p = pstar
      ep = epstar
      rog = rhogstar
!
      CALL data_exchange(rgp)
      CALL data_exchange(rlk)
      CALL data_exchange(ug)
      CALL data_exchange(us)
      IF (job_type == JOB_TYPE_3D) THEN
         CALL data_exchange(vg)
         CALL data_exchange(vs)
      END IF
      CALL data_exchange(wg)
      CALL data_exchange(ws)
      CALL data_exchange(sieg)
      CALL data_exchange(sies)
      CALL data_exchange(tg)
      CALL data_exchange(ts)
      CALL data_exchange(p)
      CALL data_exchange(ep)
      CALL data_exchange(kapgt)
      IF (iturb >= 1) THEN
        kapgt = kapgt + kapg
      ELSE IF (iturb == 0) THEN
        kapgt = kapg
      END IF
!
      DO ijk = 1, ncint
        compute = BTEST(flag(ijk),0)
        IF (compute) THEN
         CALL subscr( ijk )
         IF (ctu > 0) CALL ctu1_subscr(ijk)
         IF (ctu > 1) CALL ctu2_subscr(ijk)
         IF (ctu > 2) CALL ctu3_subscr(ijk)

! soluzione ideale: una routine che calcola allo stesso tempo tutti i flussi!

         CALL calc_all_fluxes(ijk)

        END IF


      END DO

      CALL data_exchange(ugfe)
      CALL data_exchange(ugft)
      CALL data_exchange(wgfe)
      CALL data_exchange(wgft)

      CALL data_exchange(usfe)
      CALL data_exchange(usft)
      CALL data_exchange(wsfe)
      CALL data_exchange(wsft)

      IF (job_type == JOB_TYPE_3D) THEN
        CALL data_exchange(ugfn)
        CALL data_exchange(vgfe)
        CALL data_exchange(vgfn)
        CALL data_exchange(vgft)
        CALL data_exchange(wgfn)

        CALL data_exchange(usfn)
        CALL data_exchange(vsfe)
        CALL data_exchange(vsfn)
        CALL data_exchange(vsft)
        CALL data_exchange(wsfn)
      END IF

      CALL data_exchange(egfe)
      CALL data_exchange(egft)
      CALL data_exchange(esfe)
      CALL data_exchange(esft)

      IF (job_type == JOB_TYPE_3D) THEN
        CALL data_exchange(egfn)
        CALL data_exchange(esfn)
      END IF
!
      RETURN
!
      END SUBROUTINE compute_fluxes
!----------------------------------------------------------------------
      SUBROUTINE solve_linear_systems(ijk)
!
! ... Use Gauss-Jordan method for matrix inversion
!
      USE dimensions, ONLY: nsolid
      USE enthalpy_matrix, ONLY: flim
      USE grid, ONLY: flag
      USE particles_constants, ONLY: plim, dk
      USE set_indexes, ONLY: ipjk, ijpk, ijkp
      USE gas_solid_temperature, ONLY: siegn, siesn
      USE gas_solid_density, ONLY: rgpn, rlkn
      USE time_parameters, ONLY: dt
      USE eos_gas, ONLY: cg
      USE specific_heat_module, ONLY: ck

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ijk
      INTEGER :: eq, l, ll, nphase, lj, li, lp1, i0, is
      LOGICAL :: compute
      REAL*8, DIMENSION(:), ALLOCATABLE :: bu, bv, bw, bt
      REAL*8 :: amul, div, slim, num, den, coeff, tmp, rhsg, cgammag
      REAL*8, DIMENSION(nsolid) :: cgamma, rhs, alfak, alfag

      nphase = nsolid+1

      ALLOCATE(bu(nphase), bv(nphase), bw(nphase), bt(nphase))
!
!... solve linear system of momentum-x equation
!
      compute = BTEST(flag(ipjk),0)
      IF (compute) THEN
         bu = 0.D0
         DO l = 1, nphase
            eq = nphase+l
            bu(l) = rhs3(ijk,eq)
         END DO
         DO l=2,nphase
            IF(DABS(au(l,l)) < plim(l-1)) THEN
               DO ll=1,nphase
                  au(ll,ll)=au(ll,ll)+au(l,ll)
                  au(l,ll)=0.D0
                  au(ll,l)=0.D0
               END DO
               bu(l)=0.D0
            END IF
         END DO
         !
         DO l=1,nphase
            IF(au(l,l) /= 0.D0) THEN
               lp1=l+1
               div=1.D0/au(l,l)
               DO lj=lp1,nphase
                  au(l,lj)=au(l,lj)*div
               END DO
               bu(l)=bu(l)*div
               au(l,l)=0.D0
               DO li=1,nphase
                  amul=au(li,l)
                  DO lj=lp1,nphase
                     au(li,lj)=au(li,lj)-amul*au(l,lj)
                  END DO
                  bu(li)=bu(li)-amul*bu(l)
               END DO
            END IF
         END DO
!
! ... Explicitly compute East velocities
!
         ugstar(ijk)=bu(1)
         DO l=2,nphase
            usstar(ijk,l-1)=bu(l)
         END DO
         !
      END IF
!
!... solve linear system of momentum-y equations
!
      IF (job_type == JOB_TYPE_3D) THEN
         compute = BTEST(flag(ijpk),0)
         IF (compute) THEN
            bv = 0.D0
            DO l = 1,nphase
               eq = 2*nphase+l
               bv(l) = rhs3(ijk,eq)
            END DO
            DO l=2,nphase
               IF(DABS(av(l,l)) < plim(l-1)) THEN
                  DO ll=1,nphase
                     av(ll,ll)=av(ll,ll)+av(l,ll)
                     av(l,ll)=0.D0
                     av(ll,l)=0.D0
                  END DO
                  bv(l)=0.D0
               END IF
            END DO
            !
            DO l=1,nphase
               IF(av(l,l) /= 0.D0) THEN 
                  lp1=l+1
                  div=1.D0/av(l,l)
                  DO lj=lp1,nphase
                     av(l,lj)=av(l,lj)*div
                  END DO
                  bv(l)=bv(l)*div
                  av(l,l)=0.D0
                  DO li=1,nphase
                     amul=av(li,l)
                     DO lj=lp1,nphase
                        av(li,lj)=av(li,lj)-amul*av(l,lj)
                     END DO 
                     bv(li)=bv(li)-amul*bv(l)
                  END DO 
               END IF
            END DO
!
! ... Explicitly compute North velocities
!
            vgstar(ijk)=bv(1)
            DO l=2,nphase
               vsstar(ijk,l-1)=bv(l)
            END DO
            !
         END IF
      END IF
!
!... solve linear system of momentum-z equations
!
      compute = BTEST(flag(ijkp),0)
      IF (compute) THEN
         bw = 0.D0
         DO l = 1,nphase
            IF (job_type .EQ. JOB_TYPE_3D) THEN
              eq = 3*nphase+l
            ELSE
              eq = 2*nphase+l
            END IF
            bw(l) = rhs3(ijk,eq)
         END DO

         DO l=2,nphase
            IF(DABS(aw(l,l)) < plim(l-1)) THEN
               DO ll=1,nphase
                  aw(ll,ll)=aw(ll,ll)+aw(l,ll)
                  aw(l,ll)=0.D0
                  aw(ll,l)=0.D0
               END DO
               bw(l)=0.D0
            END IF
         END DO
         !
         DO l=1,nphase
            IF(aw(l,l) /= 0.D0) THEN 
               lp1=l+1
               div=1.D0/aw(l,l)
               DO lj=lp1,nphase
                  aw(l,lj)=aw(l,lj)*div
               END DO
               bw(l)=bw(l)*div
               aw(l,l)=0.D0
               DO li=1,nphase
                  amul=aw(li,l)
                  DO lj=lp1,nphase
                     aw(li,lj)=aw(li,lj)-amul*aw(l,lj)
                  END DO 
                  bw(li)=bw(li)-amul*bw(l)
               END DO 
            END IF
         END DO
!
! ... Explicitly compute Top velocities
!
         wgstar(ijk)=bw(1)
         DO l=2,nphase
            wsstar(ijk,l-1)=bw(l)
         END DO
         !
      END IF
!
!... solve linear system of enthalpy equations
!
      IF (job_type == JOB_TYPE_2D) THEN
        i0 = 3*nphase
      ELSE
        i0 = 4*nphase
      ENDIF

      DO is = 1, nphase
        eq = i0+is
        bt(is) = rhs3(ijk,eq)
      END DO

      IF (MINVAL(rlkstar(ijk,:)) > limit) THEN
!
        cgammag = rgpstar(ijk)
        rhsg = bt(1)
        DO is = 1,nsolid
          cgamma(is) = rlkstar(ijk,is)
          rhs(is) = bt(is+1)
          alfag(is) = hv(ijk,is)*dt/cg(ijk)
          alfak(is) = hv(ijk,is)*dt/ck(is,ijk)
        END DO
!
        num = rhsg
        den = cgammag
        DO is = 1,nsolid
          tmp = cgamma(is)+alfak(is)
          tmp = 1.D0/tmp 
          coeff = cgamma(is)*tmp
          num = num + rhs(is)
          num = num - coeff*rhs(is)
          den = den + coeff*alfag(is)
        END DO
        den = 1.D0/den
        siegstar(ijk) = num*den

        DO is = 1,nsolid
          num = rhs(is) + alfag(is)*siegstar(ijk)
          tmp = cgamma(is) + alfak(is)
          den = 1.0D0/tmp
          siesstar(ijk,is) = num*den
        END DO

      ELSE

        DO is=nphase,2,-1
          IF(ABS(at(is,is)) < flim) THEN
            at(1,is)=0.D0
            at(is,1)=0.D0
            bt(is)=0.D0
          ELSE
!
! ... eliminate all cross elements in the gas enthalpy equation
!
            div=1.D0/at(is,is)
            at(is,1)=at(is,1)*div
            bt(is)=bt(is)*div
            bt(1)=bt(1)-at(1,is)*bt(is)
            at(1,1)=at(1,1)-at(1,is)*at(is,1)
          ENDIF
        END DO
!
        bt(1)=bt(1)/at(1,1)
        DO is=2,nphase
          bt(is)=bt(is)-at(is,1)*bt(1)
        END DO
!
        siegstar(ijk) = bt(1)
        DO is=1, nsolid
          siesstar(ijk,is) = bt(is+1)
        END DO
      END IF
!
      DEALLOCATE(bu)
      DEALLOCATE(bv)
      DEALLOCATE(bw)
      DEALLOCATE(bt)

      RETURN

      END SUBROUTINE solve_linear_systems
!----------------------------------------------------------------------
      SUBROUTINE init_star
!
! ... star-variables are initialized with solution at the previous time step
!
      USE dimensions, ONLY: nsolid
      USE gas_solid_density, ONLY: rgp, rlk, rog
      USE gas_solid_temperature, ONLY: tg, ts, sieg, sies
      USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
      USE particles_constants, ONLY: inrl
      USE pressure_epsilon, ONLY: p, ep

      IMPLICIT NONE

      INTEGER :: is

      rgpstar = rgp
      rlkstar = rlk
      ugstar = ug
      usstar = us
      IF (job_type == JOB_TYPE_3D) THEN
         vgstar = vg
         vsstar = vs
      END IF
      wgstar = wg
      wsstar = ws
      tgstar = tg
      tsstar = ts
      siegstar = sieg
      siesstar = sies
      rhogstar = rog
      epstar = ep
      pstar = p
      RETURN

      END SUBROUTINE init_star
!----------------------------------------------------------------------
      SUBROUTINE init_rhs
!
! ... Initialize right hand side: in each computational cell we initialize the right
!     end side of the ODE system with the solution at the previous time step
!
      USE dimensions, ONLY: nsolid, ngas
      USE domain_mapping, ONLY: ncint, meshinds
      USE gas_solid_density, ONLY: rgpn, rlkn
      USE gas_solid_temperature, ONLY: siegn, siesn
      USE gas_solid_velocity, ONLY: ug, us, vg, vs, wg, ws
      USE grid, ONLY: dx, dy, dz, flag
      USE set_indexes, ONLY: ijke, ijkn, ijkt, subscr

      IMPLICIT NONE

      INTEGER :: eq, ijk, ig, is, i, j, k, imesh
      LOGICAL :: compute
      REAL*8 :: dxi, dxip1, dyj, dyjp1, dzk, dzkp1
      REAL*8 :: dxp, dyp, dzp, indxp, indyp, indzp
      REAL*8 :: rold
!
      DO ijk = 1,ncint
        compute = BTEST(flag(ijk),0)
        IF (compute) THEN

          rhs3(ijk,:) = 0.D0
          k1(ijk,:) = 0.D0
          k2(ijk,:) = 0.D0
          k3(ijk,:) = 0.D0

          CALL meshinds(ijk,imesh,i,j,k)
          CALL subscr(ijk)

          dxi   = dx(i)
          dxip1 = dx(i+1)
          dyj   = dy(j)
          dyjp1 = dy(j+1)
          dzk   = dz(k)
          dzkp1 = dz(k+1)
!
          dxp = dxi+dxip1
          dyp = dyj+dyjp1
          dzp = dzk+dzkp1
          indxp = 1.D0/dxp
          indyp = 1.D0/dyp
          indzp = 1.D0/dzp
          !
          ! ... Gas and particle mass equations
          !
          eq = 1
          rhs3(ijk,eq) = rgpn(ijk)
          eq = eq+1
          DO is = 1,nsolid
            rhs3(ijk,eq) = rlkn(ijk,is)
            eq = eq+1
          END DO
          !
          ! ... Gas and particle momentum-x equations
          !
          rold = (dxi*rgpn(ijke)+dxip1*rgpn(ijk))*indxp
          rhs3(ijk,eq) = rold*ugn(ijk)
          eq = eq+1
          DO is = 1,nsolid
            rold = (dxi*rlkn(ijke,is)+dxip1*rlkn(ijk,is))*indxp
            rhs3(ijk,eq) = rold*usn(ijk,is)
            eq = eq+1
          END DO
          !
          ! ... Gas and particle momentum-y equations
          !
          IF (job_type == JOB_TYPE_3D) THEN
            rold = (dyj*rgpn(ijkn)+dyjp1*rgpn(ijk))*indyp
            rhs3(ijk,eq) = rold*vgn(ijk)
            eq = eq+1
            DO is = 1,nsolid
              rold = (dyj*rlkn(ijkn,is)+dyjp1*rlkn(ijk,is))*indyp
              rhs3(ijk,eq) = rold*vsn(ijk,is)
              eq = eq+1
            END DO
          END IF
          !
          ! ... Gas and particle momentum-z equations
          !
          rold = (dzk*rgpn(ijkt)+dzkp1*rgpn(ijk))*indzp
          rhs3(ijk,eq) = rold*wgn(ijk)
          eq = eq+1
          DO is = 1,nsolid
            rold = (dzk*rlkn(ijkt,is)+dzkp1*rlkn(ijk,is))*indzp
            rhs3(ijk,eq) = rold*wsn(ijk,is)
            eq = eq+1
          END DO
          !
          ! ... Gas and particle enthalpy equations
          !
          rhs3(ijk,eq) = rgpn(ijk)*siegn(ijk)
          eq = eq+1
          DO is = 1,nsolid
            rhs3(ijk,eq) = rlkn(ijk,is)*siesn(ijk,is)
            eq = eq+1
          END DO     
        END IF   
      END DO

      RETURN

      END SUBROUTINE init_rhs
!----------------------------------------------------------------------
      SUBROUTINE update_star(ijk, k1, k2, rka)

      USE dimensions
      USE domain_mapping, ONLY: meshinds
      USE eos_gas, ONLY: cg, xgc, caloric_eosg, ygc, thermal_eosg
      USE eos_solid, ONLY: caloric_eosl
      USE gas_constants, ONLY: gmw, rgas, gas_type
      USE pressure_epsilon, ONLY: pn
      USE gas_solid_density, ONLY: rgpn, rlkn
      USE gas_solid_temperature, ONLY: siegn, siesn
      USE gas_solid_velocity, ONLY: ug, us, vg, vs, wg, ws
      USE grid, ONLY: dx, dy, dz
      USE particles_constants, ONLY: cps, inrl
      USE set_indexes, ONLY: ijke, ijkn, ijkt, subscr
      USE specific_heat_module, ONLY: ck, cp
      USE time_parameters, ONLY: dt
      USE pressure_epsilon, ONLY: epn, pn

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ijk
      REAL*8, DIMENSION(:), INTENT(IN) :: k1, k2, rka

      INTEGER :: eq, info, is, ig, i, j, k, imesh
      REAL*8 :: epstot, epsstar, rtilde, mg, temp, rold, rnow, tmp
      REAL*8 :: dxi, dxip1, dyj, dyjp1, dzk, dzkp1
      REAL*8 :: dxp, dyp, dzp, indxp, indyp, indzp
      REAL*8 :: xgcl(max_ngas), rlktmp, siegtmp, siestmp

      CALL meshinds(ijk,imesh,i,j,k)

      CALL subscr(ijk)

      dxi   = dx(i)
      dxip1 = dx(i+1)
      dyj   = dy(j)
      dyjp1 = dy(j+1)
      dzk   = dz(k)
      dzkp1 = dz(k+1)
!
      dxp = dxi+dxip1
      dyp = dyj+dyjp1
      dzp = dzk+dzkp1
      indxp = 1.D0/dxp
      indyp = 1.D0/dyp
      indzp = 1.D0/dzp

      rlktmp = 0.D0
      siegtmp = 0.D0
      siestmp = 0.D0
!
! ... gas mass equation
!
      eq = 1
      rgpstar(ijk) = rgpn(ijk) + rka(1)*dt*k1(eq) &
                               + rka(2)*dt*k2(eq) 
      eq = eq+1
!
! ... particles mass equation
!
      DO is = 1,nsolid
        rlktmp = rlkn(ijk,is) + rka(1)*dt*k1(eq) &
                              + rka(2)*dt*k2(eq) 
        IF (rlktmp < 0.D0) THEN
          WRITE(testunit,*) 'clipping in updatestar rlk<0, &
                             ijk = ', ijk, ' rlktmp =', rlktmp
        END IF
        rlkstar( ijk, is ) = MAX( 0.0d0, rlktmp )
        eq = eq+1
      END DO
!
! ... gas momentum-x equation
!
      rold = (dxi*rgpn(ijke) + dxip1*rgpn(ijk))*indxp
      temp = rold*ugn(ijk) + rka(1)*dt*k1(eq) &
                           + rka(2)*dt*k2(eq) 
      rnow = (dxi*rgpstar(ijke) + dxip1*rgpstar(ijk))*indxp
      tmp = 1.D0/rnow
      ugstar(ijk) = temp*tmp
      eq = eq+1
!
! ... particles momentum-x equation
!
      DO is = 1,nsolid
        rold = (dxi*rlkn(ijke,is) + dxip1*rlkn(ijk,is))*indxp
        temp = rold*usn(ijk,is) + rka(1)*dt*k1(eq) &
                                + rka(2)*dt*k2(eq) 
        rnow = (dxi*rlkstar(ijke,is) + dxip1*rlkstar(ijk,is))*indxp
        IF (rnow > limit) THEN
          tmp = 1.D0/rnow
          usstar(ijk,is) = temp*tmp
        ELSE
          usstar(ijk,is) = 0.D0
        END IF
        eq = eq+1
      END DO
!
! ... gas momentum-y equation
!
      IF (job_type == JOB_TYPE_3D) THEN
        rold = (dyj*rgpn(ijkn) + dyjp1*rgpn(ijk))*indyp
        temp = rold*vgn(ijk) + rka(1)*dt*k1(eq) &
                             + rka(2)*dt*k2(eq) 
        rnow = (dyj*rgpstar(ijkn) + dyjp1*rgpstar(ijk))*indyp
        tmp = 1.D0/rnow
        vgstar(ijk) = temp*tmp
        eq = eq+1
!
! ... particles momentum-y equation
!
        DO is = 1,nsolid
          rold = (dyj*rlkn(ijkn,is) + dyjp1*rlkn(ijk,is))*indyp
          temp = rold*vsn(ijk,is) + rka(1)*dt*k1(eq) &
                                  + rka(2)*dt*k2(eq) 
          rnow = (dyj*rlkstar(ijkn,is) + dyjp1*rlkstar(ijk,is))*indyp
          IF (rnow > limit) THEN
            tmp = 1.D0/rnow
            vsstar(ijk,is) = temp*tmp
          ELSE
            vsstar(ijk,is) = 0.D0
          END IF
          eq = eq+1
        END DO
      END IF
!
! ... gas momentum-z equation
!
      rold = (dzk*rgpn(ijkt) + dzkp1*rgpn(ijk))*indzp
      temp = rold*wgn(ijk) + rka(1)*dt*k1(eq) &
                           + rka(2)*dt*k2(eq) 
      rnow = (dzk*rgpstar(ijkt) + dzkp1*rgpstar(ijk))*indzp
      tmp = 1.D0/rnow
      wgstar(ijk) = temp*tmp
      eq = eq+1
!
! ... particles momentum-z equation
!
      DO is = 1,nsolid
        rold = (dzk*rlkn(ijkt,is) + dzkp1*rlkn(ijk,is))*indzp
        temp = rold*wsn(ijk,is) + rka(1)*dt*k1(eq) &
                                + rka(2)*dt*k2(eq) 
        rnow = (dzk*rlkstar(ijkt,is) + dzkp1*rlkstar(ijk,is))*indzp
        IF (rnow > limit) THEN
          tmp = 1.D0/rnow
          wsstar(ijk,is) = temp*tmp
        ELSE
          wsstar(ijk,is) = 0.D0
        END IF
        eq = eq+1
      END DO
!
! ... gas enthalpy equation
!
      siegtmp = rgpn(ijk)*siegn(ijk) + rka(1)*dt*k1(eq) &
                                     + rka(2)*dt*k2(eq) 
      tmp = 1.D0/rgpstar(ijk)
      siegstar(ijk) = siegtmp*tmp
      eq = eq+1
!
! ... particles enthalpy equation
!
      DO is = 1,nsolid
        siestmp = rlkn(ijk,is)*siesn(ijk,is) + rka(1)*dt*k1(eq) &
                                             + rka(2)*dt*k2(eq) 
        IF (rlkstar(ijk,is) > limit) THEN
          tmp = 1.D0/rlkstar(ijk,is)
          siesstar(ijk,is) = siestmp*tmp
        ELSE
          siesstar(ijk,is) = 0.D0
        END IF
        eq = eq+1
      END DO
!
! ... compute temperatures from EOS
!
      CALL caloric_eosg(cp(:,ijk), cg(ijk), tgstar(ijk), ygc(ijk,:), siegstar(ijk), ijk, info)

      IF (info > 0) then
             write(testunit,*) 'error in update star'
      call write_star(ijk)
      end if

      DO is=1, nsolid
        CALL caloric_eosl(tsstar(ijk,is),cps(is),ck(is,ijk),siesstar(ijk,is),ijk,info) 
      END DO
!
! ... compute volume fractions, gas density and pressure from closure relations
!
      epstot  = 0.D0
      epsstar = 0.D0
      DO is = 1,nsolid
        epsstar = rlkstar(ijk,is)*inrl(is)
        epstot = epstot + epsstar
      END DO
      epstar(ijk) = 1.0D0-epstot

      IF (epstar(ijk) > 1.D0) then
        write(testunit,*) 'ep > 1 in update star iter in ijk =', ijk
        call write_star(ijk)
      end if

      IF (epstar(ijk) < 0.D0) then
        write(testunit,*) 'ep < 0 in update star iter in ijk =', ijk
        call write_star(ijk)
      end if

      tmp = 1.D0/epstar(ijk)
      rhogstar(ijk) = rgpstar(ijk)/epstar(ijk)

      mg = 0.D0
      DO ig = 1, ngas
         mg = mg + xgc(ijk,ig) * gmw(gas_type(ig))
      END DO
      pstar(ijk) = rhogstar(ijk)*tgstar(ijk)*rgas/mg

      RETURN

      END SUBROUTINE update_star
!----------------------------------------------------------------------
      SUBROUTINE compute_errors(ijk)
!
!... compute the error between solutions obtained with Runge-Kutta order 2 and 
!    Runge-Kutta order 3. The error estimate will be used for time-step adaptation
!
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ijk

      RETURN

      END SUBROUTINE compute_errors
!----------------------------------------------------------------------
      SUBROUTINE save_star(ijk)

      USE dimensions, ONLY: nsolid
      USE gas_solid_velocity, ONLY: ug, us, vg, vs, wg, ws
      USE gas_solid_temperature, ONLY: tg, ts, sieg, sies
      USE gas_solid_density, ONLY: rgp, rlk, rog
      USE pressure_epsilon, ONLY: ep, p

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ijk

      INTEGER :: is

      ep(ijk)   = epstar(ijk)
      p(ijk)    = pstar(ijk)
      rog(ijk)  = rhogstar(ijk)
      ug(ijk)   = ugstar(ijk)
      IF (job_type .EQ. JOB_TYPE_3D) vg(ijk) = vgstar(ijk)
      wg(ijk)   = wgstar(ijk)
      rgp(ijk)  = rgpstar(ijk)
      tg(ijk)   = tgstar(ijk)
      sieg(ijk) = siegstar(ijk)
      DO is = 1,nsolid
        rlk(ijk,is)  = rlkstar(ijk,is)
        us(ijk,is)   = usstar(ijk,is)
        IF (job_type .EQ. JOB_TYPE_3D) vs(ijk,is) = vsstar(ijk,is)
        ws(ijk,is)   = wsstar(ijk,is)
        ts(ijk,is)   = tsstar(ijk,is)
        sies(ijk,is) = siesstar(ijk,is)
      END DO

      RETURN

      END SUBROUTINE save_star
!----------------------------------------------------------------------
      SUBROUTINE write_star(ijk)

      USE pressure_epsilon, ONLY: epn, pn
      USE gas_solid_temperature, ONLY: siegn, siesn
      USE gas_solid_density, ONLY: rgpn, rlkn

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ijk
      INTEGER :: cell
      cell = ijk

      WRITE(testunit,*) 'epstar = ', epstar(cell)
      WRITE(testunit,*) 'pstar = ', pstar(cell)
      WRITE(testunit,*) 'rgpstar = ', rgpstar(cell)
      WRITE(testunit,*) 'rhogstar = ', rhogstar(cell)
      WRITE(testunit,*) 'ugstar = ', ugstar(cell)
      WRITE(testunit,*) 'wgstar = ', wgstar(cell)
      WRITE(testunit,*) 'siegstar = ', siegstar(cell)
      WRITE(testunit,*) 'tgstar = ', tgstar(cell)
      WRITE(testunit,*) 'rlkstar = ', rlkstar(cell,:)
      WRITE(testunit,*) 'usstar = ', usstar(cell,:)
      WRITE(testunit,*) 'wsstar = ', wsstar(cell,:)
      WRITE(testunit,*) 'siesstar = ', siesstar(cell,:)
      WRITE(testunit,*) 'tsstar = ', tsstar(cell,:)
      WRITE(testunit,*) 'k1 = ', k1(cell,:)
      WRITE(testunit,*) 'k2 = ', k2(cell,:)
      WRITE(testunit,*) 'k3 = ', k3(cell,:)
      WRITE(testunit,*) 'rhs3 = ', rhs3(cell,:)
      WRITE(testunit,*) 'epn = ', epn(cell)
      WRITE(testunit,*) 'pn = ', pn(cell)
      WRITE(testunit,*) 'rgpn = ', rgpn(cell)
      WRITE(testunit,*) 'siegn = ', siegn(cell)
      WRITE(testunit,*) 'rlkn = ', rlkn(cell,:)
      WRITE(testunit,*) 'siesn = ', siesn(cell,:)
      WRITE(testunit,*) 'ugn = ', ugn(cell)
      WRITE(testunit,*) 'wgn = ', wgn(cell)
      WRITE(testunit,*) 'usn = ', usn(cell,:)
      WRITE(testunit,*) 'wsn = ', wsn(cell,:)
      

      END SUBROUTINE write_star
!----------------------------------------------------------------------
      SUBROUTINE compute_courant

      USE domain_mapping, ONLY: ncint
      USE gas_constants, ONLY: gammaair
      USE pressure_epsilon, ONLY: p
      USE gas_solid_density, ONLY: rog
      USE gas_solid_velocity, ONLY: ug, wg
      USE grid, ONLY: indx, indy, indz
      USE domain_mapping, ONLY: meshinds
      USE time_parameters, ONLY: dt
      USE eos_gas, ONLY: scsound
      USE parallel, ONLY: parallel_hangup

      IMPLICIT NONE

      REAL*8 :: css, cs, courantx, courantz, cfl
      INTEGER :: ijk, imesh, i, j, k

      DO ijk = 1, ncint

        CALL meshinds(ijk,imesh,i,j,k)

        css = scsound(gammaair, rog(ijk), p(ijk))
        cs  = SQRT(css)
        courantx = dt*(DABS(ug(ijk))+cs)*indx(i)
        courantz = dt*(DABS(wg(ijk))+cs)*indz(k)
        courant(ijk) = MAX(courantx,courantz)

      END DO

      cfl = MAXVAL(courant)

      IF (cfl > 1.0D0) THEN
        WRITE(testunit,*) 'Warning! CFL > 1 in explicit_iter, cfl = ', cfl
!         CALL parallel_hangup
!         STOP
      END IF

      RETURN

      END SUBROUTINE compute_courant
!----------------------------------------------------------------------
      SUBROUTINE compute_physical_coeff(ijk)

      USE dimensions, ONLY: nsolid
      USE eos_gas, ONLY: xgc, cg
      USE gas_solid_viscosity, ONLY: viscon, mug, kapg
      USE heat_transfer, ONLY: hvs
      USE momentum_transfer, ONLY: kdrags, inter
      USE set_indexes, ONLY: subscr, imjk, ijmk, ijkm

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ijk
      INTEGER :: is
      REAL*8 :: dugs, dvgs, dwgs
      REAL*8, DIMENSION(:), ALLOCATABLE :: nul

      ALLOCATE( nul( ( ( nsolid + 1 )**2 + ( nsolid + 1 ) ) / 2 ) )
      nul = 0.D0

      CALL subscr(ijk)
!
! ... compute the mixture viscosity and the thermal conductivity that are functions of
!     temperature and mixture composition
!
      CALL viscon( mug(ijk), kapg(ijk), xgc(ijk,:), tgstar(ijk) )
!
! ... compute gas-particle heat transfer coefficient and gas-particle drag coefficient
!
      DO is = 1, nsolid

        dugs = ((ugstar(ijk)-usstar(ijk,is))+(ugstar(imjk)-usstar(imjk,is)))*0.5D0
        dwgs = ((wgstar(ijk)-wsstar(ijk,is))+(wgstar(ijkm)-wsstar(ijkm,is)))*0.5D0
        IF (job_type == JOB_TYPE_2D) THEN
          dvgs = 0.D0
        ELSE IF (job_type == JOB_TYPE_3D) THEN
          dvgs = ((vgstar(ijk)-vsstar(ijk,is))+(vgstar(ijmk)-vsstar(ijmk,is)))*0.5D0
        END IF

        CALL hvs(hv(ijk,is), rlkstar(ijk,is), rhogstar(ijk), epstar(ijk), &
                 dugs, dvgs, dwgs, mug(ijk), kapg(ijk), cg(ijk), is)

        CALL kdrags(kpgv(ijk,is), dugs, dvgs, dwgs, epstar(ijk),     &
                    rgpstar(ijk), rlkstar(ijk,is), mug(ijk), is)

      END DO
!
! ... Compute the particle-particle coefficients and the interphase matrix
!
      IF (job_type == JOB_TYPE_2D ) THEN
        CALL inter(appu(ijk,:), nul(:), appw(ijk,:), kpgv(ijk,:), usstar, &
                   usstar, wsstar, rlkstar, ijk)
      ELSE IF (job_type == JOB_TYPE_3D) THEN
        CALL inter(appu(ijk,:), appv(ijk,:), appw(ijk,:), kpgv(ijk,:), &
                   usstar, vsstar, wsstar, rlkstar, ijk)
      END IF

      RETURN

      END SUBROUTINE compute_physical_coeff
!----------------------------------------------------------------------
      SUBROUTINE calc_all_fluxes(ijk)

        USE set_indexes, ONLY: nb, rnb, first_nb
        USE set_indexes, ONLY: ctu1_nb, ctu1_rnb, ctu2_nb, ctu2_rnb, ctu3_nb, ctu3_rnb
        USE set_indexes, ONLY: stencil, cte
        USE set_indexes, ONLY: maskval, maskstencil, masks
        USE set_indexes, ONLY: imjk, ijmk, ijkm
        USE gas_solid_density, ONLY: rgp, rlk
        USE gas_solid_temperature, ONLY: sieg, sies, tg, ts
        USE gas_solid_velocity, ONLY: ug, vg, wg, us, vs, ws
        USE convective_mass_fluxes, ONLY: fmas, masf, ctu1_fmas, ctu2_fmas, ctu3_fmas
        USE dimensions
        USE flux_limiters, ONLY: muscl, ctu
        USE grid, ONLY: flag, fluid, bl_cell, immb_cell
        USE control_flags, ONLY: job_type, JOB_TYPE_2D, JOB_TYPE_3D
        USE domain_mapping, ONLY: meshinds
        USE convective_fluxes_u, ONLY: flu, muscl_flu, ctu1_flu, ctu2_flu, ctu3_flu
        USE convective_fluxes_v, ONLY: flv, muscl_flv,  ctu1_flv, ctu2_flv, ctu3_flv
        USE convective_fluxes_w, ONLY: flw, muscl_flw, ctu1_flw, ctu2_flw, ctu3_flw
        USE convective_fluxes_sc, ONLY: fsc, muscl_fsc, ctu1_fsc, ctu2_fsc, ctu3_fsc
        USE diffusive_fluxes, ONLY: hotc
        USE gas_solid_viscosity, ONLY: gas_viscosity, part_viscosity
        USE iterative_solver, ONLY: rgfe, rgfn, rgft, rsfe, rsfn, rsft
        USE tilde_momentum, ONLY: ugfe, ugft, ugfn, vgfe, vgft, vgfn, wgfe, wgft, wgfn
        USE tilde_momentum, ONLY: usfe, usft, usfn, vsfe, vsft, vsfn, wsfe, wsft, wsfn
        USE tilde_energy, ONLY: egfe, egft, egfn
        USE tilde_energy, ONLY: hgfe, hgft, hgfn
        USE tilde_energy, ONLY: esfe, esft, esfn
        USE tilde_energy, ONLY: hsfe, hsft, hsfn
        USE immersed_boundaries, ONLY: immb, b_e_, b_n_, b_t_
        USE pressure_epsilon, ONLY: p, ep
        USE turbulence_model, ONLY: kapgt, iturb
        USE particles_constants, ONLY: inrl, kap
        USE set_indexes, ONLY: OPERATOR( * )
        USE interpolate_fields, ONLY: interpolate_x, interpolate_y, interpolate_z


        IMPLICIT NONE

        INTEGER, INTENT(IN) :: ijk
        INTEGER :: i, j, k, is, imesh
        TYPE(masks) :: mask_x, mask_y, mask_z
        TYPE(stencil) :: dens_stagx, dens_stagy, dens_stagz
        TYPE(stencil) :: u, v, w, dens, enth
        TYPE(stencil) :: eps, temp, kappa
        REAL*8 :: zero

        CALL meshinds(ijk,imesh,i,j,k)

        IF (job_type == JOB_TYPE_2D) THEN

          ! ... assemble computational stencils

          CALL rnb(u,ug,ijk)
          CALL rnb(w,wg,ijk)
          CALL nb(dens,rgp,ijk)
          CALL nb(enth,sieg,ijk)
          IF (ctu > 0) THEN
            CALL ctu1_nb ( dens, rgp, ijk )
            CALL ctu1_rnb( u, ug, ijk )
            CALL ctu1_rnb( w, wg, ijk )
            CALL ctu1_nb ( enth, sieg, ijk )
            IF (ctu > 1) THEN
              CALL ctu2_rnb(u,ug,ijk)
              CALL ctu2_rnb(w,wg,ijk)
              CALL ctu2_nb(dens,rgp,ijk)
              CALL ctu2_nb(enth,sieg,ijk)
              IF (ctu > 2) THEN
                CALL ctu3_rnb(u,ug,ijk)
                CALL ctu3_rnb(w,wg,ijk)
                CALL ctu3_nb(dens,rgp,ijk)
                CALL ctu3_nb(enth,sieg,ijk)
              END IF
            END IF
          END IF

! begin: First order fluxes

          ! gas mass flux

          CALL masf( rgfe( ijk ), rgft( ijk ), rgfe( imjk ), rgft( ijkm ), &
                    dens, u, w, ijk )

          ! gas momentum flux
          !
          ! ... Stencil Mask for Immersed Boundaries
          !
          IF (immb >= 1) THEN
            mask_x = maskval(1)
            mask_z = maskval(1)
            IF (flag(ijk)==bl_cell .OR. flag(ijk)==immb_cell) THEN
              !
              ! ... Stencil masks
              CALL rnb(mask_x,b_e_,ijk)
              CALL rnb(mask_z,b_t_,ijk)
            END IF
          END IF

          ! ... Mask non-physical velocities from Immersed Boundaries
          !
          IF (immb >= 1) THEN
            IF (flag(ijk)==bl_cell .OR. flag(ijk)==immb_cell) THEN
              u = maskstencil(u,mask_x)
              w = maskstencil(w,mask_z)
            END IF
          END IF
!
          ! ... Interpolate density on the staggered grid
          CALL interpolate_x(dens, dens_stagx, i)
          CALL interpolate_z(dens, dens_stagz, k)

          CALL flu(ugfe(ijk), ugft(ijk), ugfe(imjk), ugft(ijkm),  & 
                   dens_stagx, u, w, i)

          CALL flw(wgfe(ijk), wgft(ijk), wgfe(imjk), wgft(ijkm),  &
                   dens_stagz, u, w, i, k)

          ! gas energy flux

          CALL fsc(egfe(ijk), egft(ijk),    &
                   egfe(imjk), egft(ijkm),    &
                   enth, dens, u, w, ijk)

! end: First order fluxes

! begin: Muscl fluxes
!
! ... Second order MUSCL correction
!
          IF (muscl > 0 .AND. flag(ijk)==fluid) THEN

            ! gas mass flux

            CALL fmas( rgfe( ijk ), rgft( ijk ), rgfe( imjk ), rgft( ijkm ), &
                       dens, u, w, ijk )

            ! gas momentum flux

            IF ( i /= nx-1 ) CALL muscl_flu(ugfe(ijk), ugft(ijk),   &
                                            dens_stagx, u, w, i, k)
            IF ( k /= nz-1 ) CALL muscl_flw(wgfe(ijk), wgft(ijk),   &
                                            dens_stagz, u, w, i, k)

            ! gas energy flux

            CALL muscl_fsc(egfe(ijk), egft(ijk),    & 
                           enth, dens, u, w, ijk)

          END IF
!
! ... First order Corner Transport Upwind correction (step 2)
!
          IF (ctu > 0) THEN

            CALL ctu1_fmas(rgfe( ijk ), rgft( ijk ), rgfe( imjk ), rgft( ijkm ), &
                           dens, u, w, ijk)
            CALL ctu1_flu(ugfe(ijk), ugft(ijk), dens_stagx, u, w, i, k)
            CALL ctu1_flw(wgfe(ijk), wgft(ijk), dens_stagz, u, w, i, k)
            CALL ctu1_fsc(egfe(ijk), egft(ijk), dens, enth, u, w, ijk)
!
! ... Second order Corner Transport Upwind correction (step 3)
!
            IF (ctu > 1 .AND. muscl == 0 .AND. flag(ijk)==fluid) THEN

              CALL ctu2_fmas(rgfe( ijk ), rgft( ijk ), rgfe( imjk ), rgft( ijkm ), &
                      dens, u, w, ijk)
              IF (i /= nx-1) CALL ctu2_flu(ugfe(ijk), ugft(ijk), &
                                           dens_stagx, u, w, i, k)
              IF (k /= nz-1) CALL ctu2_flw(wgfe(ijk), wgft(ijk), &
                                           dens_stagz, u, w, i, k)
              CALL ctu2_fsc(egfe(ijk), egft(ijk), dens, enth, u, w, ijk)
!
! ... Second order Corner Transport Upwind transverse correction (step 4)
!
              IF (ctu > 2) THEN

                CALL ctu3_fmas(rgfe( ijk ), rgft( ijk ), rgfe( imjk ), rgft( ijkm ), &
                        dens, u, w, ijk)
                CALL ctu3_flu(ugfe(ijk), ugft(ijk), dens_stagx, u, w, i, k)
                CALL ctu3_flw(wgfe(ijk), wgft(ijk), dens_stagz, u, w, i, k)
                CALL ctu3_fsc(egfe(ijk), egft(ijk), dens, enth, u, w, ijk)
              END IF
            END IF
          END IF
!
          IF (gas_viscosity) THEN

            ! assemble first order computational stencils
            CALL first_nb(kappa,kapgt,ijk)
            CALL first_nb(eps,ep,ijk)
            CALL first_nb(temp,tg,ijk)
      
            CALL hotc(hgfe(ijk), zero, hgft(ijk),        &
                      hgfe(imjk), zero, hgft(ijkm),      &
                      eps, temp, kappa, ijk)
          END IF

! ... repeat all for particles

          DO is = 1, nsolid

            ! ... assemble computational stencils

            CALL nb(dens,rlk(:,is),ijk)
            CALL rnb(u,us(:,is),ijk)
            CALL rnb(w,ws(:,is),ijk)
            CALL nb(enth,sies(:,is),ijk)
            IF (ctu > 0) THEN
              CALL ctu1_rnb(u,us(:,is),ijk)
              CALL ctu1_rnb(w,ws(:,is),ijk)
              CALL ctu1_nb(dens,rlk(:,is),ijk)
              CALL ctu1_nb(enth,sies(:,is),ijk)
              IF (ctu > 1) THEN
                CALL ctu2_rnb(u,us(:,is),ijk)
                CALL ctu2_rnb(w,ws(:,is),ijk)
                CALL ctu2_nb(dens,rlk(:,is),ijk)
                CALL ctu2_nb(enth,sies(:,is),ijk)
                IF (ctu > 2) THEN
                  CALL ctu3_rnb(u,us(:,is),ijk)
                  CALL ctu3_rnb(w,ws(:,is),ijk)
                  CALL ctu3_nb(dens,rlk(:,is),ijk)
                  CALL ctu3_nb(enth,sies(:,is),ijk)
                END IF
              END IF
            END IF


! begin: first order fluxes

            ! particles mass flux

            CALL masf(rsfe(ijk,is),  rsft(ijk,is),    &
                    rsfe(imjk,is), rsft(ijkm,is),   &
                    dens, u, w, ijk)

            ! particles momentum flux

            ! ... Mask non-physical velocities from Immersed Boundaries
            !
            IF (immb >= 1) THEN
              IF (flag(ijk)==bl_cell .OR. flag(ijk)==immb_cell) THEN
                u = maskstencil(u,mask_x)
                w = maskstencil(w,mask_z)
              END IF
            END IF

            ! ... Interpolate density on the staggered grid
            CALL interpolate_x(dens, dens_stagx, i)
            CALL interpolate_z(dens, dens_stagz, k)

            CALL flu(usfe(ijk,is), usft(ijk,is), usfe(imjk,is), usft(ijkm,is), &
                     dens_stagx, u, w, i)

            CALL flw(wsfe(ijk,is), wsft(ijk,is), wsfe(imjk,is), wsft(ijkm,is), &
                     dens_stagz, u, w, i, k)

            ! particles energy flux

            CALL fsc(esfe(ijk, is), esft(ijk, is),  &
                     esfe(imjk, is), esft(ijkm, is),  &
                     enth, dens, u, w, ijk)
! end: first order fluxes

!
! ... Second order MUSCL correction
!
            IF (muscl > 0 .AND. flag(ijk)==fluid) THEN

              ! particles mass flux

              CALL fmas(rsfe(ijk,is), rsft(ijk,is),    &
                      rsfe(imjk,is), rsft(ijkm,is),   &
                      dens, u, w, ijk)

              ! particles momentum flux

              IF ( i /= nx-1 ) CALL muscl_flu(usfe(ijk,is), usft(ijk,is),  &
                                              dens_stagx, u, w, i, k)
              IF ( k /= nz-1 ) CALL muscl_flw(wsfe(ijk,is), wsft(ijk,is),  &
                                              dens_stagz, u, w, i, k)

              ! particles energy flux

              CALL muscl_fsc(esfe(ijk, is), esft(ijk, is),  &
                             enth, dens, u, w, ijk)
            END IF
!
! ... First order Corner Transport Upwind correction (step 2)
!
            IF (ctu > 0) THEN

              CALL ctu1_fmas(rsfe(ijk,is), rsft(ijk,is),    &
                      rsfe(imjk,is), rsft(ijkm,is), &
                      dens, u, w, ijk)
              CALL ctu1_flu(usfe(ijk,is), usft(ijk,is), dens_stagx, u, w, i, k)
              CALL ctu1_flw(wsfe(ijk,is), wsft(ijk,is), dens_stagz, u, w, i, k)
              CALL ctu1_fsc(esfe(ijk,is), esft(ijk,is), dens, enth, u, w, ijk)
!
! ... Second order Corner Transport Upwind correction (step 3)
!
              IF (ctu > 1  .AND. muscl == 0 .AND. flag(ijk) == fluid) THEN

                CALL ctu2_fmas(rsfe(ijk,is), rsft(ijk,is),    &
                        rsfe(imjk,is), rsft(ijkm,is), &
                        dens, u, w, ijk)
                IF (i/=nx-1) CALL ctu2_flu(usfe(ijk,is), usft(ijk,is), &
                                           dens_stagx, u, w, i, k)
                IF (k/=nz-1) CALL ctu2_flw(wsfe(ijk,is), wsft(ijk,is), &
                                           dens_stagz, u, w, i, k)
                CALL ctu2_fsc(esfe(ijk, is), esft(ijk, is), dens, enth, u, w, ijk)

                IF (ctu > 2) THEN
  
                  CALL ctu3_fmas(rsfe(ijk,is), rsft(ijk,is), rsfe(imjk,is), rsft(ijkm,is), &
                          dens, u, w, ijk)
                  CALL ctu3_flu(usfe(ijk,is), usft(ijk,is), dens_stagx, u, w, i, k)
                  CALL ctu3_flw(wsfe(ijk,is), wsft(ijk,is), dens_stagz, u, w, i, k)
                  CALL ctu3_fsc(esfe(ijk,is), esft(ijk,is), dens, enth, u, w, ijk)
  
                END IF
              END IF
            END IF

!
! ... Diffusive fluxes
!
            IF (part_viscosity) THEN

              ! ... assemble first order computational stencils
              eps   = inrl(is) * dens
              kappa = cte(kap(is))
              CALL first_nb(temp,ts(:,is),ijk)
            
              CALL hotc(hsfe(ijk,is), zero, hsft(ijk, is),    &
                        hsfe(imjk,is), zero, hsft(ijkm, is),    &
                        eps, temp, kappa, ijk)
            END IF
                     
          END DO

        ELSE IF (job_type == JOB_TYPE_3D) THEN

          ! ... assemble computational stencils

          CALL rnb(u,ug,ijk)
          CALL rnb(v,vg,ijk)
          CALL rnb(w,wg,ijk)
          CALL nb(dens,rgp,ijk)
          CALL nb(enth,sieg,ijk)
          IF (ctu > 0) THEN
            CALL ctu1_nb ( dens, rgp, ijk )
            CALL ctu1_rnb( u, ug, ijk )
            CALL ctu1_rnb( v, vg, ijk )
            CALL ctu1_rnb( w, wg, ijk )
            CALL ctu1_nb ( enth, sieg, ijk )
            IF (ctu > 1) THEN
              CALL ctu2_rnb(u,ug,ijk)
              CALL ctu2_rnb(v,vg,ijk)
              CALL ctu2_rnb(w,wg,ijk)
              CALL ctu2_nb(dens,rgp,ijk)
              CALL ctu2_nb(enth,sieg,ijk)
              IF (ctu > 2) THEN
                CALL ctu3_rnb(u,ug,ijk)
                CALL ctu3_rnb(v,vg,ijk)
                CALL ctu3_rnb(w,wg,ijk)
                CALL ctu3_nb(dens,rgp,ijk)
                CALL ctu3_nb(enth,sieg,ijk)
              END IF
            END IF
          END IF
!
          ! gas mass flux

          CALL masf( rgfe( ijk  ), rgfn( ijk  ), rgft( ijk  ),  &
                     rgfe( imjk ), rgfn( ijmk ), rgft( ijkm ),  &
                     dens, u, v, w, ijk )


          ! gas momentum flux
          !
          ! ... Stencil Mask for Immersed Boundaries
          !
          IF (immb >= 1) THEN
            mask_x = maskval(1)
            mask_y = maskval(1)
            mask_z = maskval(1)
            IF (flag(ijk)==bl_cell .OR. flag(ijk)==immb_cell) THEN
              !
              ! ... Stencil masks
              CALL rnb(mask_x,b_e_,ijk)
              CALL rnb(mask_y,b_n_,ijk)
              CALL rnb(mask_z,b_t_,ijk)
            END IF
          END IF

          ! ... Mask non-physical velocities from Immersed Boundaries
          !
          IF (immb >= 1) THEN
            IF (flag(ijk)==bl_cell .OR. flag(ijk)==immb_cell) THEN
              u = maskstencil(u,mask_x)
              v = maskstencil(v,mask_y)
              w = maskstencil(w,mask_z)
            END IF
          END IF

          ! ... Interpolate density on the staggered grid
          CALL interpolate_x(dens, dens_stagx, i)
          CALL interpolate_y(dens, dens_stagy, j)
          CALL interpolate_z(dens, dens_stagz, k)

          CALL flu(ugfe(ijk), ugfn(ijk), ugft(ijk),                    &
                   ugfe(imjk), ugfn(ijmk), ugft(ijkm),                 &
                   dens_stagx, u, v, w, i)

          CALL flv(vgfe(ijk), vgfn(ijk), vgft(ijk),                    &
                   vgfe(imjk), vgfn(ijmk), vgft(ijkm),                 &
                   dens_stagy, u, v, w, j)

          CALL flw(wgfe(ijk), wgfn(ijk), wgft(ijk),                    &
                   wgfe(imjk), wgfn(ijmk), wgft(ijkm),                 &
                   dens_stagz, u, v, w, k)

          ! gas energy flux

          CALL fsc(egfe(ijk), egfn(ijk), egft(ijk),    &
                   egfe(imjk), egfn(ijmk), egft(ijkm),    &
                   enth, dens, u, v, w, ijk)
!
! ... Second order MUSCL correction
!
          IF (muscl > 0 .AND. flag(ijk)==fluid) THEN

            ! gas mass flux

            CALL fmas( rgfe( ijk  ), rgfn( ijk  ), rgft( ijk  ),  &
                       rgfe( imjk ), rgfn( ijmk ), rgft( ijkm ),  &
                       dens, u, v, w, ijk )

            ! gas momentum flux

            IF ( i /= nx-1 ) &
              CALL muscl_flu(ugfe(ijk), ugfn(ijk), ugft(ijk), &
                             dens_stagx, u, v, w, i, j, k)
            IF ( j /= ny-1 ) &
              CALL muscl_flv(vgfe(ijk), vgfn(ijk), vgft(ijk), &
                             dens_stagy, u, v, w, i, j, k)
            IF ( k /= nz-1 ) &
              CALL muscl_flw(wgfe(ijk), wgfn(ijk), wgft(ijk), &
                             dens_stagz, u, v, w, i, j, k)

            ! gas energy flux

            CALL muscl_fsc(egfe(ijk), egfn(ijk), egft(ijk),    &
                           enth, dens, u, v, w, ijk)

          END IF
!
! ... First order Corner Transport Upwind correction (step 2)
!
          IF (ctu > 0) THEN

            CALL ctu1_fmas( rgfe( ijk  ), rgfn( ijk  ), rgft( ijk  ),  &
                       rgfe( imjk ), rgfn( ijmk ), rgft( ijkm ),  &
                       dens, u, v, w, ijk)

            CALL ctu1_flu(ugfe(ijk), ugfn(ijk), ugft(ijk), &
                       dens_stagx, u, v, w, i, j, k)
            CALL ctu1_flv(vgfe(ijk), vgfn(ijk), vgft(ijk), &
                       dens_stagy, u, v, w, i, j, k)
            CALL ctu1_flw(wgfe(ijk), wgfn(ijk), wgft(ijk), &
                       dens_stagz, u, v, w, i, j, k)
            CALL ctu1_fsc(egfe(ijk), egfn(ijk), egft(ijk),    &
                      enth, dens, u, v, w, ijk)
!
! ... Second order Corner Transport Upwind correction (step 3)
!
            IF (ctu > 1 .AND. muscl == 0 .AND. flag(ijk)==fluid) THEN

              CALL ctu2_fmas( rgfe( ijk  ), rgfn( ijk  ), rgft( ijk  ),  &
                         rgfe( imjk ), rgfn( ijmk ), rgft( ijkm ),  &
                         dens, u, v, w, ijk)
              IF (i /= nx-1) CALL ctu2_flu(ugfe(ijk), ugfn(ijk), ugft(ijk), &
                                         dens_stagx, u, v, w, i, j, k)
              IF (j /= ny-1) CALL ctu2_flv(vgfe(ijk), vgfn(ijk), vgft(ijk), &
                                         dens_stagy, u, v, w, i, j, k)
              IF (k /= nz-1) CALL ctu2_flw(wgfe(ijk), wgfn(ijk), wgft(ijk), &
                                         dens_stagz, u, v, w, i, j, k)
              CALL ctu2_fsc(egfe(ijk), egfn(ijk), egft(ijk),    &
                        enth, dens, u, v, w, ijk)
!
! ... Second order Corner Transport Upwind transverse correction (step 4)
!
              IF (ctu > 2) THEN
                CALL ctu3_fmas( rgfe( ijk  ), rgfn( ijk  ), rgft( ijk  ),  &
                         rgfe( imjk ), rgfn( ijmk ), rgft( ijkm ),  &
                         dens, u, v, w, ijk)
                CALL ctu3_flu(ugfe(ijk), ugfn(ijk), ugft(ijk), &
                         dens_stagx, u, v, w, i, j, k)
                CALL ctu3_flv(vgfe(ijk), vgfn(ijk), vgft(ijk), &
                         dens_stagy, u, v, w, i, j, k)
                CALL ctu3_flw(wgfe(ijk), wgfn(ijk), wgft(ijk), &
                         dens_stagz, u, v, w, i, j, k)
                CALL ctu3_fsc(egfe(ijk), egfn(ijk), egft(ijk),    &
                         enth, dens, u, v, w, ijk)
              END IF
            END IF
          END IF
!
! ... Diffusive fluxes
!
          IF (gas_viscosity) THEN

            ! assemble first order computational stencils
            CALL first_nb(kappa,kapgt,ijk)
            CALL first_nb(eps,ep,ijk)
            CALL first_nb(temp,tg,ijk)
      
            CALL hotc(hgfe(ijk), hgfn(ijk), hgft(ijk),         &
                      hgfe(imjk), hgfn(ijmk), hgft(ijkm),      &
                      eps, temp, kappa, ijk)
          END IF

          DO is = 1, nsolid

            ! ... assemble computational stencils

            CALL nb(dens,rlk(:,is),ijk)
            CALL rnb(u,us(:,is),ijk)
            CALL rnb(v,vs(:,is),ijk)
            CALL rnb(w,ws(:,is),ijk)
            CALL nb(enth,sies(:,is),ijk)
            IF (ctu > 0) THEN
              CALL ctu1_rnb(u,us(:,is),ijk)
              CALL ctu1_rnb(v,vs(:,is),ijk)
              CALL ctu1_rnb(w,ws(:,is),ijk)
              CALL ctu1_nb(dens,rlk(:,is),ijk)
              CALL ctu1_nb(enth,sies(:,is),ijk)
              IF (ctu > 1) THEN
                CALL ctu2_rnb(u,us(:,is),ijk)
                CALL ctu2_rnb(v,vs(:,is),ijk)
                CALL ctu2_rnb(w,ws(:,is),ijk)
                CALL ctu2_nb(dens,rlk(:,is),ijk)
                CALL ctu2_nb(enth,sies(:,is),ijk)
                IF (ctu > 2) THEN
                  CALL ctu3_rnb(u,us(:,is),ijk)
                  CALL ctu3_rnb(v,vs(:,is),ijk)
                  CALL ctu3_rnb(w,ws(:,is),ijk)
                  CALL ctu3_nb(dens,rlk(:,is),ijk)
                  CALL ctu3_nb(enth,sies(:,is),ijk)
                END IF
              END IF
            END IF
            ! particles mass flux

            CALL masf(rsfe(ijk,is),  rsfn(ijk,is),  rsft(ijk,is),    &
                      rsfe(imjk,is), rsfn(ijmk,is), rsft(ijkm,is),   &
                      dens, u, v, w, ijk)

            ! particles momentum flux

            ! ... Mask non-physical velocities from Immersed Boundaries
            !
            IF (immb >= 1) THEN
              IF (flag(ijk)==bl_cell .OR. flag(ijk)==immb_cell) THEN
                u = maskstencil(u,mask_x)
                v = maskstencil(v,mask_y)
                w = maskstencil(w,mask_z)
              END IF
            END IF

            ! ... Interpolate density on the staggered grid
            CALL interpolate_x(dens, dens_stagx, i)
            CALL interpolate_y(dens, dens_stagy, j)
            CALL interpolate_z(dens, dens_stagz, k)

            CALL flu(usfe(ijk,is), usfn(ijk,is), usft(ijk,is),         &
                     usfe(imjk,is), usfn(ijmk,is), usft(ijkm,is),      &
                     dens_stagx, u, v, w, i)

            CALL flv(vsfe(ijk,is), vsfn(ijk,is), vsft(ijk,is),         &
                     vsfe(imjk,is), vsfn(ijmk,is), vsft(ijkm,is),      &
                     dens_stagy, u, v, w, j)

            CALL flw(wsfe(ijk,is), wsfn(ijk,is), wsft(ijk,is),         &
                     wsfe(imjk,is), wsfn(ijmk,is), wsft(ijkm,is),      &
                     dens_stagz, u, v, w, k)

            ! particles energy flux

            CALL fsc(esfe(ijk, is), esfn(ijk, is), esft(ijk, is),  &
                     esfe(imjk, is), esfn(ijmk, is), esft(ijkm, is),  &
                     enth, dens, u, v, w, ijk)
!
! ... Second order MUSCL correction
!
            IF (muscl > 0 .AND. flag(ijk)==fluid) THEN

              ! particles mass flux

              CALL fmas(rsfe(ijk,is),  rsfn(ijk,is),  rsft(ijk,is),    &
                    rsfe(imjk,is), rsfn(ijmk,is), rsft(ijkm,is),   &
                    dens, u, v, w, ijk)

              ! particles momentum flux

              IF ( i /= nx-1 ) &
                CALL muscl_flu(usfe(ijk,is), usfn(ijk,is), usft(ijk,is), &
                               dens_stagx, u, v, w, i, j, k)
              IF ( j /= ny-1 ) &
                CALL muscl_flv(vsfe(ijk,is), vsfn(ijk,is), vsft(ijk,is), &
                               dens_stagy, u, v, w, i, j, k)
              IF ( k /= nz-1 ) &
                CALL muscl_flw(wsfe(ijk,is), wsfn(ijk,is), wsft(ijk,is), &
                               dens_stagz, u, v, w, i, j, k)

              ! particles energy flux

              CALL muscl_fsc(esfe(ijk, is), esfn(ijk, is), esft(ijk, is), &
                             enth, dens, u, v, w, ijk)

            END IF
!
! ... First order Corner Transport Upwind correction (step 2)
!
          IF (ctu > 0) THEN

            CALL ctu1_fmas(rsfe(ijk,is),  rsfn(ijk,is),  rsft(ijk,is),    &
                      rsfe(imjk,is), rsfn(ijmk,is), rsft(ijkm,is),   &
                      dens, u, v, w, ijk)
            CALL ctu1_flu(usfe(ijk,is), usfn(ijk,is), usft(ijk,is), &
                          dens_stagx, u, v, w, i, j, k)
            CALL ctu1_flv(vsfe(ijk,is), vsfn(ijk,is), vsft(ijk,is), &
                          dens_stagy, u, v, w, i, j, k)
            CALL ctu1_flw(wsfe(ijk,is), wsfn(ijk,is), wsft(ijk,is), &
                          dens_stagz, u, v, w, i, j, k)
            CALL ctu1_fsc(esfe(ijk, is), esfn(ijk, is), esft(ijk, is),  &
                        enth, dens, u, v, w, ijk)
!
! ... Second order Corner Transport Upwind correction (step 3)
!
            IF (ctu > 1 .AND. muscl == 0 .AND. flag(ijk)==fluid) THEN

              CALL ctu2_fmas(rsfe(ijk,is),  rsfn(ijk,is),  rsft(ijk,is),    &
                        rsfe(imjk,is), rsfn(ijmk,is), rsft(ijkm,is),   &
                        dens, u, v, w, ijk)
              IF (i /= nx-1) CALL ctu2_flu(usfe(ijk,is), usfn(ijk,is), usft(ijk,is), &
                                           dens_stagx, u, v, w, i, j, k)
              IF (j /= ny-1) CALL ctu2_flv(vsfe(ijk,is), vsfn(ijk,is), vsft(ijk,is), &
                                           dens_stagy, u, v, w, i, j, k)
              IF (k /= nz-1) CALL ctu2_flw(wsfe(ijk,is), wsfn(ijk,is), wsft(ijk,is), &
                                           dens_stagz, u, v, w, i, j, k)
              CALL ctu2_fsc(esfe(ijk, is), esfn(ijk, is), esft(ijk, is),  &
                          enth, dens, u, v, w, ijk)
!
! ... Second order Corner Transport Upwind transverse correction (step 4)
!
              IF (ctu > 2) THEN

                CALL ctu3_fmas(rsfe(ijk,is),  rsfn(ijk,is),  rsft(ijk,is),    &
                          rsfe(imjk,is), rsfn(ijmk,is), rsft(ijkm,is),   &
                          dens, u, v, w, ijk)
                CALL ctu3_flu(usfe(ijk,is), usfn(ijk,is), usft(ijk,is), &
                              dens_stagx, u, v, w, i, j, k)
                CALL ctu3_flv(vsfe(ijk,is), vsfn(ijk,is), vsft(ijk,is), &
                              dens_stagy, u, v, w, i, j, k)
                CALL ctu3_flw(wsfe(ijk,is), wsfn(ijk,is), wsft(ijk,is), &
                              dens_stagz, u, v, w, i, j, k)
                CALL ctu3_fsc(esfe(ijk, is), esfn(ijk, is), esft(ijk, is),  &
                            enth, dens, u, v, w, ijk)
              END IF
            END IF
          END IF

! ... Diffusive fluxes
!
            IF (part_viscosity) THEN

              ! ... assemble first order computational stencils
              eps   = inrl(is) * dens
              kappa = cte(kap(is))
              CALL first_nb(temp,ts(:,is),ijk)
            
              CALL hotc(hsfe(ijk,is), hsfn(ijk, is), hsft(ijk, is),    &
                        hsfe(imjk,is), hsfn(ijmk, is), hsft(ijkm, is),    &
                        eps, temp, kappa, ijk)
            END IF
          END DO
        END IF

        IF (gas_viscosity) THEN

          egfe(ijk) = egfe(ijk) - hgfe(ijk)
          egft(ijk) = egft(ijk) - hgft(ijk)
!
          IF ( .NOT.BTEST(flag(imjk),0) ) egfe(imjk) = egfe(imjk) - hgfe(imjk)
          IF ( .NOT.BTEST(flag(ijkm),0) ) egft(ijkm) = egft(ijkm) - hgft(ijkm)

          IF (job_type == JOB_TYPE_3D) THEN
            egfn(ijk) = egfn(ijk) - hgfn(ijk)
            IF ( .NOT.BTEST(flag(ijmk),0) ) egfn(ijmk) = egfn(ijmk) - hgfn(ijmk)
          END IF
        END IF
!
        DO is = 1, nsolid
          IF (part_viscosity) THEN
            esfe(ijk, is) = esfe(ijk, is) - hsfe(ijk,is)
            esft(ijk, is) = esft(ijk, is) - hsft(ijk, is)
!
            IF ( .NOT.BTEST(flag(imjk),0) ) esfe(imjk,is) = esfe(imjk,is) - hsfe(imjk,is)
            IF ( .NOT.BTEST(flag(ijkm),0) ) esft(ijkm,is) = esft(ijkm,is) - hsft(ijkm,is)

            IF (job_type == JOB_TYPE_3D) THEN
              esfn(ijk, is) = esfn(ijk, is) - hsfn(ijk, is)
              IF ( .NOT.BTEST(flag(ijmk),0) ) esfn(ijmk,is)=esfn(ijmk,is)-hsfn(ijmk,is)
            END IF
          END IF
        END DO

        RETURN

      END SUBROUTINE calc_all_fluxes
!----------------------------------------------------------------------
      END MODULE explicit_solver
!----------------------------------------------------------------------
