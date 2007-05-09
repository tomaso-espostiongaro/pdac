!----------------------------------------------------------------------
      MODULE turbulence_model
!----------------------------------------------------------------------
!
      USE domain_mapping, ONLY: ncint, ncdom, data_exchange
      USE environment, ONLY: timing, cpclock
      USE indijk_module, ONLY: ip0_jp0_kp0_
      USE io_files, ONLY: errorunit

      IMPLICIT NONE
!
      REAL*8, DIMENSION(:),   ALLOCATABLE :: mugt   ! gas turbulent viscosity
      REAL*8, DIMENSION(:),   ALLOCATABLE :: kapgt  ! gas turbulent conductivity
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: must   ! solid turbulent viscosity


      REAL*8, DIMENSION(:), ALLOCATABLE :: smag, scoeff
!
      INTEGER :: iturb, iss, modturbo
      REAL*8 :: cmut                       ! Smagorinsky constant
      REAL*8 :: pranumt                    ! Prandtl number
!
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE allocate_turbo
      USE dimensions
      IMPLICIT NONE
!
      ALLOCATE(smag(ncint))
      ALLOCATE(scoeff(ncint))
      ALLOCATE(mugt(ncdom))
      ALLOCATE(must(ncdom,nsolid))
      ALLOCATE(kapgt(ncdom))
      smag = 0.0D0
      kapgt = 0.0D0
      mugt  = 0.0D0
      must  = 0.0D0
      scoeff = 0.0D0
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
       SUBROUTINE turbulence_setup
!
! ... sets the Smagorinsky turbulence length scale
! ... and the roughness length
! ... (2D/3D-Compliant)
!
        USE control_flags, ONLY: job_type
        USE dimensions, ONLY: nz, no, nx, ny
        USE domain_mapping, ONLY: myijk, meshinds
        USE grid, ONLY: dz, zb, rb, dx, dy, z
        USE volcano_topography, ONLY: topo_c, topo2d_c
        USE indijk_module, ONLY: ip0_jp0_kp0_
        USE io_files, ONLY: tempunit
        USE parallel, ONLY: mpime, root
        USE roughness_module, ONLY: roughness, roughness2d,  roughness_setup
!
        INTEGER :: i, j, k, i1, i2, j2, n, ijk, imesh
        REAL*8  :: zsgs, zrou, ls, lb, lm
        REAL*8 :: delt

        INTEGER, DIMENSION(:), ALLOCATABLE :: kt, it
        REAL*8,  DIMENSION(:), ALLOCATABLE :: zbt
!
! ... Initialize the roughness matrix
!
        IF (iturb == 2) CALL roughness_setup
!
! ... Compute the turbulent mixing length by using the roughness
! ... model and the Smagorinsky model
!
        IF( job_type == '2D' ) THEN
          DO ijk = 1, ncint
            CALL meshinds(ijk,imesh,i,j,k)
            !
            ! ... Mesh length-scale
            !
            delt = DSQRT( dz(k)*dx(i) )
            ls = cmut * delt
            !
            IF ( iturb == 1 ) THEN
                    lm = ls
            ELSE IF ( iturb == 2 ) THEN
                    zsgs = z(k) - topo_c(i)
                    IF (zsgs > 0.D0) THEN
                            zrou = roughness(i)
                            lb = 0.4D0 * ( zsgs + zrou )
                            !
                            ! ... Either use the smallest between the smag length and the roughness length
                            !
                            !lm = DMIN1(ls,lb)
                            !
                            ! ... or match the Smagorinsky and the Roughness length scales
                            !
                            lm = DSQRT(1.D0 / (1.D0/ls**2 + 1.D0/lb**2))
                    END IF
            ENDIF
            !
            ! ... Squared turbulence length scale is used into Smagorinsky model
            !
            smag(ijk) = lm**2
          END DO
        ELSE IF( job_type == '3D' ) THEN
          DO ijk = 1, ncint
            CALL meshinds(ijk,imesh,i,j,k)
            !
            delt = ( dz(k)*dy(j)*dx(i) )**(1.0d0/3.0d0)
            ls = cmut * delt
            !
            IF ( iturb == 1 ) THEN
                    lm = ls
            ELSE IF ( iturb == 2 ) THEN
                    zsgs = z(k) - topo2d_c(i,j)
                    IF (zsgs > 0.D0) THEN
                            zrou = roughness2d(i,j)
                            lb = 0.4D0 * ( zsgs + zrou )
                            !
                            ! ... Either use the smallest between the smag length and the roughness length
                            !
                            !lm = DMIN1(ls,lb)
                            !
                            ! ... or match the Smagorinsky and the Roughness length scales
                            !
                            lm = DSQRT(1.D0 / (1.D0/ls**2 + 1.D0/lb**2))
                    END IF
            ENDIF
            !
            ! ... Squared turbulence length scale is used into Smagorinsky model
            !
            smag( ijk ) = lm**2
            !
            ! ... Correction for anisotropic grids (only for classical 
            ! ... Smagorinsky model)
            !
            smag(ijk) = smag(ijk) * anis(i,j,k)
          END DO
        END IF

      RETURN
      END  SUBROUTINE turbulence_setup
!----------------------------------------------------------------------
      SUBROUTINE sgsg

! ... Sub Grid Stress (Gas)
! ... compute (dynamic) Smagorinsky coefficient , turbulent viscosity and 
! ... turbulent conductivity in the center of the cell.
! ... (2D-3D-Compliant, dynamic model only computed in 3D)

      USE control_flags, ONLY: job_type
      USE dimensions, ONLY: nx, ny, nz
      USE domain_mapping, ONLY: meshinds
      USE grid, ONLY: dx, dy, dz, flag
      USE eos_gas, ONLY: cg
      USE gas_solid_density, ONLY: rog
      USE gas_solid_velocity, ONLY: ug,vg,wg 
      USE set_indexes, ONLY: subscr, stencil, first_rnb
      USE parallel, ONLY: mpime, root
      IMPLICIT NONE
!
      REAL*8, DIMENSION(:), ALLOCATABLE :: modsr, p11, p12, p22, p13, p23, p33  
      REAL*8, DIMENSION(:), ALLOCATABLE :: fug, fvg, fwg, fu2g, fv2g, fw2g,    &
                                           fuvg, fuwg, fvwg
      REAL*8 :: sr11, sr12, sr22, sr13, sr23, sr33  
      REAL*8 :: fmodsr, fsr11, fsr12, fsr22, fsr13, fsr23, fsr33
      REAL*8 :: fp11, fp12, fp22, fp13, fp23, fp33
      REAL*8 :: l11, l12, l22, l13, l23, l33, m11, m12, m22, m13, m23, m33 
      REAL*8 :: l1d, l2d, l3d
      REAL*8 :: cdyn 
      REAL*8 :: num, den
      REAL*8 :: delt
      INTEGER  :: i, j, k, ijk, imesh  

      TYPE(stencil) :: sp11, sp12, sp22, sp13, sp23, sp33
!
      ALLOCATE( modsr(ncdom) );     modsr = 0.0D0

      IF (modturbo == 2 .AND. job_type == '2D') THEN

        IF( mpime == root ) WRITE(errorunit,*) 'WARNING!: Dynamic Smagorinsky &
                                      & model cannot be applied to 2D'
        modturbo = 1

      ELSE IF (modturbo == 2 .AND. job_type == '3D' ) THEN

        ! ... Allocate and initialize the filtered velocities for the
        ! ... dynamic computation of the Smagorinsky length scale 
        ! ... (Germano et al., 1991)
        !
        ALLOCATE(p11(ncdom))  ;     p11 = 0.0D0
        ALLOCATE(p12(ncdom))  ;     p12 = 0.0D0
        ALLOCATE(p22(ncdom))  ;     p22 = 0.0D0
        ALLOCATE(p13(ncdom))  ;     p13 = 0.0D0
        ALLOCATE(p23(ncdom))  ;     p23 = 0.0D0
        ALLOCATE(p33(ncdom))  ;     p33 = 0.0D0
        ALLOCATE(fug(ncdom))  ;     fug = 0.0D0
        ALLOCATE(fvg(ncdom))  ;     fvg = 0.0D0
        ALLOCATE(fwg(ncdom))  ;     fwg = 0.0D0
        ALLOCATE(fu2g(ncdom)) ;     fu2g = 0.0D0
        ALLOCATE(fv2g(ncdom)) ;     fv2g = 0.0D0
        ALLOCATE(fw2g(ncdom)) ;     fw2g = 0.0D0
        ALLOCATE(fuvg(ncdom)) ;     fuvg = 0.0D0
        ALLOCATE(fuwg(ncdom)) ;     fuwg = 0.0D0
        ALLOCATE(fvwg(ncdom)) ;     fvwg = 0.0D0

        DO ijk = 1, ncint

          IF( BTEST(flag(ijk),0) ) THEN

            CALL subscr(ijk)
!
            ! ... compute the modulus and components of the 
            ! ... strain tensor
            !
            CALL strain3d(ug, vg, wg, modsr(ijk),  &
                          sr11,sr12,sr22,sr13,sr23,sr33, ijk)
!	  
            p11(ijk) = modsr(ijk) * sr11
            p12(ijk) = modsr(ijk) * sr12
            p22(ijk) = modsr(ijk) * sr22
            p13(ijk) = modsr(ijk) * sr13
            p23(ijk) = modsr(ijk) * sr23
            p33(ijk) = modsr(ijk) * sr33
!
            ! ... compute filtered velocities
            !
            CALL vel_hat(ijk, fug(ijk), fvg(ijk), fwg(ijk),    &
                      fu2g(ijk), fv2g(ijk), fw2g(ijk),         &
                      fuvg(ijk), fuwg(ijk), fvwg(ijk) )        
          END IF
        END DO
!
        CALL data_exchange(p11)
        CALL data_exchange(p12)
        CALL data_exchange(p22)
        CALL data_exchange(p13)
        CALL data_exchange(p23)
        CALL data_exchange(p33)
        CALL data_exchange(fug)
        CALL data_exchange(fvg)
        CALL data_exchange(fwg)
        CALL data_exchange(fu2g)
        CALL data_exchange(fv2g)
        CALL data_exchange(fw2g)
        CALL data_exchange(fuvg)
        CALL data_exchange(fuwg)
        CALL data_exchange(fvwg)

      END IF  
!
      DO ijk = 1, ncint

        IF( BTEST(flag(ijk),0) ) THEN

          CALL meshinds( ijk, imesh, i, j, k )

          CALL subscr( ijk )
!
          IF( modturbo == 1 ) THEN
             !
             ! ... 'Classical' Smagorisnky model:
             ! ... compute the modulus of strain tensor 
             !
            IF (job_type == '2D') THEN
               CALL strain2d(ug, wg, modsr(ijk), ijk)
            ELSE IF (job_type == '3D') THEN
               CALL strain3d(ug, vg, wg, modsr(ijk),   &
                             sr11,sr12,sr22,sr13,sr23,sr33, ijk)
            END IF
!
          ELSE IF(modturbo == 2) THEN
!
             CALL strain3d(fug, fvg, fwg, fmodsr, fsr11, fsr12, fsr22,  &
                           fsr13, fsr23, fsr33, ijk)
!
             ! ... compute the local neighbours of 'p' and filter
             ! ... over the stencil 'sp'
             !
             CALL first_rnb(sp11,p11,ijk)
             fp11 = filter_3d(sp11)
             
             CALL first_rnb(sp12,p12,ijk)
             fp12 = filter_3d(sp12)   
             
             CALL first_rnb(sp22,p22,ijk)
             fp22 = filter_3d(sp22)
             
             CALL first_rnb(sp13,p13,ijk)
             fp13 = filter_3d(sp13)
             
             CALL first_rnb(sp23,p23,ijk)
             fp23 = filter_3d(sp23)

             CALL first_rnb(sp33,p33,ijk)
             fp33 = filter_3d(sp33)
!         
             l11  = fug(ijk)**2 - fu2g(ijk)
             l12  = fug(ijk)*fvg(ijk) - fuvg(ijk)
             l22  = fvg(ijk)**2 - fv2g(ijk)
             l13  = fug(ijk)*fwg(ijk) - fuwg(ijk)
             l23  = fvg(ijk)*fwg(ijk) - fvwg(ijk)
             l33  = fwg(ijk)**2 - fw2g(ijk)
!
!             l1d = l11 - 0.5D0*(l11+l22+l33)
!             l2d = l22 - 0.5D0*(l11+l22+l33)  
!             l3d = l33 - 0.5D0*(l11+l22+l33)  
!
             m11 = (4.D0*fmodsr*fsr11 - fp11)
             m12 = (4.D0*fmodsr*fsr12 - fp12)
             m22 = (4.D0*fmodsr*fsr22 - fp22)
             m13 = (4.D0*fmodsr*fsr13 - fp13)
             m23 = (4.D0*fmodsr*fsr23 - fp23)
             m33 = (4.D0*fmodsr*fsr33 - fp33)

             ! ... mesh-filter length
             !
             delt = ( dz(k)*dy(j)*dx(i) )**(1.0d0/3.0d0)

             num = l11*m11 + l22*m22 + l33*m33 + &
                   2.D0*l12*m12 + 2.D0*l13*m13 + 2.D0*l23*m23

             den = delt**2 * ( m11**2 + m22**2 + m33**2 +  &
                               2.D0*m12**2 + 2.D0*m13**2 + 2.D0*m23**2 )

             IF(den == 0.0D0) THEN
               cdyn = 0.0D0
             ELSE 
               cdyn = 0.5D0 * num / den 
             END IF

             ! ... squared Smagorinsky length scale
             !
             smag(ijk) = cdyn * (delt**2)

          END IF
!
! ... Smagorinsky turbulent viscosity: for classical model
! ... smag is computed during setup
!
          mugt(ijk) = rog(ijk) * smag(ijk) * modsr(ijk)

          IF(mugt(ijk) < 0.0D0) mugt(ijk) = 0.0D0
!
          scoeff(ijk) = cdyn
!
! ... Turbulent heat conductivity
!
          pranumt=0.5
          kapgt(ijk)=mugt(ijk)*cg(ijk)/pranumt
!
        ELSE

          mugt(ijk)  = 0.D0
          kapgt(ijk) = 0.D0

        END IF

      END DO
!                          
      DEALLOCATE(modsr)
      IF (modturbo == 2) THEN
        DEALLOCATE(p11)
        DEALLOCATE(p12)
        DEALLOCATE(p22)
        DEALLOCATE(p13)
        DEALLOCATE(p23)
        DEALLOCATE(p33)
        DEALLOCATE(fug)
        DEALLOCATE(fvg)
        DEALLOCATE(fwg)
        DEALLOCATE(fu2g)
        DEALLOCATE(fv2g)
        DEALLOCATE(fw2g)
        DEALLOCATE(fuvg)
        DEALLOCATE(fuwg)
        DEALLOCATE(fvwg)
      END IF

      RETURN 
      END SUBROUTINE sgsg
!----------------------------------------------------------------------
      SUBROUTINE vel_hat(ijk, fu, fv, fw, fu2, fv2, fw2, fuv, fuw, fvw) 
!
! ... Compute the filtered components of the velocity in the cell centers
! ... Values in the cell centers are computed through the 'sumstencil' interface 
! ... Cross terms (uv, uw, vw) are computed through the 'prodstencil' interface
! ... (in set_indexes module)

      USE set_indexes
      USE gas_solid_velocity, ONLY: ug,vg,wg 
      IMPLICIT NONE
 
      INTEGER, INTENT(IN) :: ijk
      REAL*8, INTENT(OUT) :: fu, fv, fw, fu2, fv2, fw2, fuv, fuw, fvw
      TYPE(stencil) :: um, up, vm, vp, wm, wp
      TYPE(stencil) :: u, v, w, u2, v2, w2, uv, uw, vw
!
      CALL first_rnb(up,ug,ijk)
      CALL first_rnb(um,ug,imjk)
      u = 0.5D0 * (up + um) 

      CALL first_rnb(vp,vg,ijk)
      CALL first_rnb(vm,vg,ijmk)
      v = 0.5D0 * (vp + vm)

      CALL first_rnb(wp,wg,ijk)
      CALL first_rnb(wm,wg,ijkm)
      w = 0.5D0 * (wp + wm)

      u2 = u * u
      v2 = v * v
      w2 = w * w
      uv = u * v
      uw = u * w
      vw = v * w
      fu  = filter_3d(u)
      fv  = filter_3d(v)
      fw  = filter_3d(w)
      fu2 = filter_3d(u2)
      fv2 = filter_3d(v2)
      fw2 = filter_3d(w2)
      fuv = filter_3d(uv)
      fuw = filter_3d(uw)
      fvw = filter_3d(vw)
!
      RETURN
      END SUBROUTINE vel_hat
!----------------------------------------------------------------------
      REAL*8 FUNCTION filter_x(q)
      USE set_indexes
!
      TYPE(stencil),  INTENT(IN) :: q
!
      filter_x = 0.25D0 * ( q%w + 2.D0 * q%c + q%e )
!
      END FUNCTION filter_x
!----------------------------------------------------------------------
      REAL*8 FUNCTION filter_y(q)
      USE set_indexes
!
      TYPE(stencil),  INTENT(IN) :: q
!
      filter_y = 0.25D0 * ( q%s + 2.D0 * q%c + q%n )
!
      END FUNCTION filter_y
!----------------------------------------------------------------------
      REAL*8 FUNCTION filter_z(q)
      USE set_indexes
!
      TYPE(stencil),  INTENT(IN) :: q
!
      filter_z = 0.25D0 * ( q%b + 2.D0 * q%c + q%t )
!
      END FUNCTION filter_z
!----------------------------------------------------------------------
      REAL*8 FUNCTION filter_3d(field)
      USE set_indexes
!
! ... 3d filter is made by linear combination of 1d filters
!
      TYPE(stencil),  INTENT(IN) :: field
      REAL*8 :: fx, fy, fz
!
      fx = filter_x(field)
      fy = filter_y(field)
      fz = filter_z(field)
!
      filter_3d = 1.D0/3.D0 * (fx + fy + fz)
!
      END FUNCTION filter_3d
!----------------------------------------------------------------------
      SUBROUTINE sgss
! ... Computes the sub-grid-stress for solids (Gidaspow, 1994)
! ... This model is no longer used
! ... (2D-3D-Compliant)
!
      USE dimensions
      USE gas_solid_density, ONLY: rlk
      USE gas_solid_velocity, ONLY: us, vs, ws
      USE grid, ONLY: itc, flag
      USE particles_constants, ONLY: rl, inrl, dk, cmus
      USE set_indexes
      USE control_flags, ONLY: job_type
      IMPLICIT NONE
!
      INTEGER :: ijk, is
      REAL*8 :: modsr, sr11,sr12,sr22,sr13,sr23,sr33
      REAL*8 :: exp, fact, cst
!
      DO ijk = 1, ncint
        IF( BTEST(flag(ijk),0) ) THEN
         CALL subscr(ijk)
!
         DO is = 1, nsolid
           IF (job_type == '2D') THEN
             CALL strain2d(us(:,is), ws(:,is), modsr, ijk)
           ELSE IF (job_type == '3D') THEN
             CALL strain3d(us(:,is), vs(:,is), ws(:,is),                &
                           modsr,sr11,sr12,sr22,sr13,sr23,sr33, ijk)
           END IF
!
           cst  = 0.5D0*DSQRT(2.D0)
           exp  = ( 3.98D0*rlk(ijk,is)*inrl(is)-1.69D0 )
           fact = 10.D0**(exp)
           
           must(ijk,is) = cst * cmus(is) * rl(is) * dk(is)**2 * fact * modsr

         END DO
        END IF
      END DO
!
      RETURN
      END SUBROUTINE
!-----------------------------------------------------------------------
      SUBROUTINE strain3d(u, v, w, modsr, sr11,sr12,sr22,sr13,sr23,sr33, ijk)
!
! ... here computes the components of the gas strain rate tensor and its module.

      USE dimensions, ONLY: nx, ny, nz   
      USE domain_mapping, ONLY: myijk
      USE grid, ONLY: dx, dy, dz, indx, indy, indz
      USE set_indexes
      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: ijk
      REAL*8, INTENT(IN), DIMENSION(:) :: u, v, w
      REAL*8, INTENT(OUT) :: modsr
      REAL*8, INTENT(OUT) :: sr11, sr12, sr22, sr13, sr23, sr33
!
      REAL*8 :: w1, w2, w3, u1, u2, u3, v1, v2, v3
      REAL*8 :: um1, um2, vm1, vm2, wm1, wm2 
      REAL*8 :: dxp, dyp, dzp, dxm, dym, dzm
      REAL*8 :: uzp, uzm, uyp, uym, vxp, vxm, vzp, vzm, wxp, wxm, wyp, wym 
      INTEGER :: i, j, k, imesh 
!
      imesh = myijk( ip0_jp0_kp0_, ijk)
      i = MOD( MOD( imesh - 1, nx*ny ), nx ) + 1
      j = MOD( imesh - 1, nx*ny ) / nx + 1
      k = ( imesh - 1 ) / ( nx*ny ) + 1
!
      dxp=dx(i)+dx(i+1)
      dxm=dx(i)+dx(i-1)
      dyp=dy(j)+dy(j+1)
      dym=dy(j)+dy(j-1)
      dzp=dz(k)+dz(k+1)
      dzm=dz(k)+dz(k-1)
!
! ... Compute velocity gradients in the cell centers.
! ... Cross terms (non-diagonal) in the strain tensor are obtained
! ... by interpolating the values found on the staggered grid.
!
      sr11 = (u(ijk)-u(imjk))*indx(i)
      sr22 = (v(ijk)-v(ijmk))*indy(j)
      sr33 = (w(ijk)-w(ijkm))*indz(k)
!
      uzp=0.5D0*(u(ijkp)+u(imjkp))
      uzm=0.5D0*(u(ijkm)+u(imjkm))
      uyp=0.5D0*(u(ijpk)+u(imjpk))
      uym=0.5D0*(u(ijmk)+u(imjmk))
      vxp=0.5D0*(v(ipjk)+v(ipjmk))
      vxm=0.5D0*(v(imjk)+v(imjmk))
      vzp=0.5D0*(v(ijkp)+v(ijmkp))
      vzm=0.5D0*(v(ijkm)+v(ijmkm))
      wxp=0.5D0*(w(ipjk)+w(ipjkm))
      wxm=0.5D0*(w(imjk)+w(imjkm))
      wyp=0.5D0*(w(ijpk)+w(ijpkm))
      wym=0.5D0*(w(ijmk)+w(ijmkm))
!
      sr12 = (uyp - uym) / (dyp+dym) + (vxp - vxm) / (dxp+dxm)
      sr13 = (uzp - uzm) / (dzp+dzm) + (wxp - wxm) / (dxp+dxm)
      sr23 = (vzp - vzm) / (dzp+dzm) + (wyp - wym) / (dyp+dym)
!        
      modsr = DSQRT(2.D0 * (sr11**2 + sr22**2 + sr33**2 +             &
                            2.D0*sr12**2 + 2.D0*sr13**2 + 2.D0*sr23**2) )
!
      RETURN
      END SUBROUTINE strain3d
!-----------------------------------------------------------------------
      SUBROUTINE strain2d(u, w, modsr, ijk)
!
! ... here computes the components of the strain rate tensor and its module.

      USE dimensions, ONLY: nx    
      USE domain_mapping, ONLY: myijk
      USE grid, ONLY:dx, dz, indx, indz,itc,inr 
      USE set_indexes, ONLY: imjk, ijkm, ijkp, imjkp, imjkm, ipjk, ipjkm 
      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: ijk
      REAL*8, INTENT(IN), DIMENSION(:) :: u, w
      REAL*8, INTENT(OUT) :: modsr
!
      REAL*8 :: w1, w2, w3, u1, u2, u3, um1, um2, wm1, wm2 
      REAL*8 :: dxp, dzm, dxm, dzp, indxp, indzm, indxm, indzp
      REAL*8 :: sr1, sr2, sr12
      REAL*8 :: d33
      INTEGER :: i, k, imesh 
!
        imesh = myijk( ip0_jp0_kp0_, ijk)
        k = ( imesh - 1 ) / nx + 1
        i = MOD( ( imesh - 1 ), nx) + 1

        dxp=dx(i)+dx(i+1)
        dxm=dx(i)+dx(i-1)
        dzp=dz(k)+dz(k+1)
        dzm=dz(k)+dz(k-1)
        indxp=1.D0/dxp
        indxm=1.D0/dxm
        indzp=1.D0/dzp
        indzm=1.D0/dzm
!
! ... Compute velocity gradients in the center of the cell.
! ... Cross terms (non-diagonal) in the strain tensor are obtained
! ... by interpolating the values found on the staggered grids.
!
        sr1 = (u(ijk)-u(imjk))*indx(i)
        sr2 = (w(ijk)-w(ijkm))*indz(k)
        
! ... extra-term for cylindrical coordinates
!
        IF(itc == 1) THEN
          d33 = 0.5D0*(u(ijk)+u(imjk))*inr(i)
        ELSE
          d33 = 0.D0
        END IF
!
        u3=0.5D0*(u(ijkp)+u(imjkp))
        u2=0.5D0*(u(ijk)+u(imjk))
        u1=0.5D0*(u(ijkm)+u(imjkm))
        w3=0.5D0*(w(ipjk)+w(ipjkm))
        w2=0.5D0*(w(ijk)+w(ijkm))
        w1=0.5D0*(w(imjk)+w(imjkm))
        um2=u2+(u3-u2)*dz(k)*indzp
        um1=u1+(u2-u1)*dz(k-1)*indzm
        wm2=w2+(w3-w2)*dx(i)*indxp
        wm1=w1+(w2-w1)*dx(i-1)*indxm

        sr12 =((um2-um1)*indz(k)+(wm2-wm1)*indx(i))*0.5D0
!        
        modsr = DSQRT(2.D0 * (sr1**2 + sr2**2 + 2.D0 * sr12**2 +   &
                              d33**2) )
!
      RETURN
      END SUBROUTINE strain2d
!---------------------------------------------------------------------
! ... Correction of Smagorinsky length scale for anisotropic grids 
! ... [Scotti, Meneveau and Lilly, Phys. Fluids A 5(9), 1993]
!
      REAL*8 FUNCTION anis(i,j,k)
      USE grid, ONLY: dx, dy, dz
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: i,j,k
      REAL*8 :: a1, a2, a3
      REAL*8 :: la1, la2, la3, fla
      REAL*8 :: deltamax, delta1, delta2, delta3
      
      delta1 = dx(i)
      delta2 = dy(j)
      delta3 = dz(k)
      deltamax = MAX(delta1, delta2, delta3)

      a1 = delta1 / deltamax
      a2 = delta2 / deltamax
      a3 = delta3 / deltamax

      la1 = log(a1)
      la2 = log(a2)
      la3 = log(a3)

      fla = la1**2 + la2**2 + la3**2 - la1*la2 +la1*la3 - la2*la3
      fla = 4.D0 / 27.D0 * fla
      fla = DSQRT( fla )

      anis = COSH( fla )

      RETURN
      END FUNCTION anis
!---------------------------------------------------------------------
      END MODULE turbulence_model
!---------------------------------------------------------------------
