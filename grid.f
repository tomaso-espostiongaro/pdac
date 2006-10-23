!----------------------------------------------------------------------
      MODULE grid
!
! ... This module contains all the information about the rectilinear
! ... 2D or 3D mesh used for computation, and the cell types 
!
!----------------------------------------------------------------------

      USE io_files, ONLY: errorunit, testunit, logunit

      IMPLICIT NONE
!
      REAL*8, DIMENSION(:), ALLOCATABLE :: dx, dy, dz
      REAL*8, DIMENSION(:), ALLOCATABLE :: indx, indy, indz
      REAL*8, DIMENSION(:), ALLOCATABLE :: r, rb
      REAL*8, DIMENSION(:), ALLOCATABLE :: inr, inrb
      REAL*8, DIMENSION(:), ALLOCATABLE :: x, xb
      REAL*8, DIMENSION(:), ALLOCATABLE :: y, yb
      REAL*8, DIMENSION(:), ALLOCATABLE :: z, zb
!
      INTEGER :: itc

      TYPE blbody
        INTEGER :: xlo
        INTEGER :: xhi
        INTEGER :: ylo
        INTEGER :: yhi
        INTEGER :: zlo
        INTEGER :: zhi
        INTEGER :: typ
      END TYPE
!
      TYPE (blbody), ALLOCATABLE :: iob(:)
!
! ... flag defining the cell types
      INTEGER, DIMENSION(:), ALLOCATABLE :: fl   ! temporary global array
      INTEGER, DIMENSION(:), ALLOCATABLE :: flag ! local flag array
      LOGICAL, DIMENSION(:), ALLOCATABLE :: fc2
!
! ... flags for the domain boundary condition
      INTEGER :: west, east, south, north, bottom, top
!
      INTEGER, PARAMETER :: fluid          = 1            ! OLD flag = 1 NEW flag = 2**0
      INTEGER, PARAMETER :: slip_wall      = 2            ! OLD flag = 2 NEW flag = 2**1
      INTEGER, PARAMETER :: noslip_wall    = 4            ! OLD flag = 3 NEW flag = 2**2
      INTEGER, PARAMETER :: free_io        = 8            ! OLD flag = 4 NEW flag = 2**3
      INTEGER, PARAMETER :: inlet_cell     = 16           ! OLD flag = 5 NEW flag = 2**4
      INTEGER, PARAMETER :: nrfree_io      = 32           ! OLD flag = 6 NEW flag = 2**5
      INTEGER, PARAMETER :: vent_cell      = 64           ! OLD flag = 7 NEW flag = 2**6
      INTEGER, PARAMETER :: dome_cell      = 129          ! OLD flag = 8 NEW flag = 2**7 + fluid
      ! ... immersed boundaries
      INTEGER, PARAMETER :: immb_cell      = 257          !              NEW flag = 2**8 + fluid
      ! ... non-computed immersed boundaries
      INTEGER, PARAMETER :: filled_cell_1    = 512        !              NEW flag = 2**9
      INTEGER, PARAMETER :: filled_cell_2    = 1024       !              NEW flag = 2**10
      INTEGER, PARAMETER :: bl_cell          = 2049       !              NEW flag = 2**11 + fluid
!
! ... origin of atmospheric stratification
      REAL*8 :: zzero
!
! ... variables for grid generator:
! ... maximum cell size increase rate
      REAL*8 ::  maxbeta
!
! ... domain size along each axis
      REAL*8 :: domain_x, domain_y, domain_z
!
! ... minimum and maximum cell sizes
! ... number of cells with minimum size;
! ... coordinate indexes of the vent center
      REAL*8  :: dxmin, dxmax, dymin, dymax, dzmin, dzmax
      INTEGER :: n0x, n0y, n0z
      INTEGER :: iv, jv, kv
!
! ... relative ( 0.0 < c < 1.0 ) center on each axis for refinement
      REAL*8 :: alpha_x, alpha_y, alpha_z
!
! ... coordinates of the mesh 'center'
      REAL*8 :: center_x, center_y
!
! ... flag for grid input (0), automatic generation (1)
      INTEGER :: grigen
!
! ... Number of cells in the inner grid and
! ... number of cells to add in the new grid
      INTEGER :: nx_inner, ny_inner, nz_inner
      INTEGER :: npx, nmx, npy, nmy, npz, nmz
      REAL*8, DIMENSION(:), ALLOCATABLE :: dx_inner, dy_inner, dz_inner
      INTEGER, ALLOCATABLE :: new_ijk(:)
      LOGICAL :: remap_grid
!
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE allocate_grid
!
      USE dimensions, ONLY: nx, ny, nz, ntot, ntr
      USE control_flags, ONLY: job_type
      USE parallel, ONLY: mpime, root
      USE time_parameters, ONLY: itd
!
      IMPLICIT NONE
!
! ... Set the appropriate total number of cells
! ... and the coordinate system

      nx_inner = nx
      ny_inner = ny
      nz_inner = nz
      nx = nx_inner + npx + nmx
      ny = ny_inner + npy + nmy
      nz = nz_inner + npz + nmz
      remap_grid = (npx+nmx+npy+nmy+npz+nmz /= 0)
!
      IF( job_type == '2D' ) THEN
        ntot = nx*nz
        ntr  = nx
      ELSE IF( job_type == '3D' ) THEN
        ntot = nx*ny*nz
        ntr  = nx*ny
      ELSE
        CALL error( ' allocate_grid ', ' wrong job_type '//job_type, 1)
      END IF

      ALLOCATE( dx(nx), dy(ny), dz(nz) ) 
      ALLOCATE( dx_inner(nx_inner), dy_inner(ny_inner), dz_inner(nz_inner) ) 
      ALLOCATE( indx(nx), indy(ny), indz(nz) )
      ALLOCATE( r(nx), rb(nx) )
      ALLOCATE( inr(nx), inrb(nx) )
      ALLOCATE( x(nx), xb(nx) )
      ALLOCATE( y(ny), yb(ny) )
      ALLOCATE( z(nz), zb(nz) )
      ALLOCATE( fl(ntot) )

      IF( mpime == root ) THEN
        WRITE(logunit,*) 
        WRITE(logunit,*) 'Simulation Grid      : ', nx, ny, nz
        WRITE(logunit,*) 'Total number of cells: ', ntot
      END IF
!
      RETURN
      END SUBROUTINE allocate_grid
!----------------------------------------------------------------------
      SUBROUTINE grid_remap(map)
      USE dimensions, ONLY: nx, ny
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: map(:)
      INTEGER :: i, j, k, ijk
!
      DO k = 1, nz_inner
        DO j = 1, ny_inner
          DO i = 1, nx_inner
            ijk = i + (j-1) * nx_inner + (k-1) * nx_inner * ny_inner
            map(ijk) = (i+nmx) + (j+nmy-1)*nx + (k+nmz-1)*nx*ny
          END DO
        END DO
      END DO
!
      RETURN
      END SUBROUTINE grid_remap
!----------------------------------------------------------------------
      SUBROUTINE allocate_blbody
      USE dimensions
      IMPLICIT NONE
      
      IF (no > 0) ALLOCATE(iob(no))
      
      RETURN
      END SUBROUTINE allocate_blbody
!----------------------------------------------------------
      SUBROUTINE grid_setup
!
      USE control_flags, ONLY: job_type, lpr, prog
      USE dimensions, ONLY: nx, ny, nz
      USE io_files, ONLY: tempunit, testunit
      USE parallel, ONLY: mpime, root
      USE time_parameters, ONLY: itd

      REAL*8, PARAMETER :: VERYBIG = 1.0d+10
      INTEGER :: i, j, k, ijk
      REAL*8 :: zrif
!
! ... Set the cell dimensions
!
      SELECT CASE (grigen)
      CASE(0)
        !
        ! ... if the cell dimensions are prescribed in input, 
        ! ... compute the domain parameters
        !
        dxmax = MAXVAL(dx_inner)
        dymax = MAXVAL(dy_inner)
        dzmax = MAXVAL(dz_inner)
        dxmin = MINVAL(dx_inner)
        dymin = MINVAL(dy_inner)
        dzmin = MINVAL(dz_inner)
        domain_x = SUM(dx_inner)
        domain_y = SUM(dy_inner)
        domain_z = SUM(dz_inner)
        iv = nx/2+1; jv = ny/2+1; kv = nz/2+1
        !
      CASE(1)
        !
        ! ... generate the non-uniform mesh
        !
        CALL generate_grid(dx_inner,nx_inner,domain_x,alpha_x,dxmin,dxmax,n0x,iv) 
        IF (job_type == '2D') THEN
          dy = 0.D0
        ELSE IF (job_type == '3D') THEN
          CALL generate_grid(dy_inner,ny_inner,domain_y,alpha_y,dymin,dymax,n0y,jv) 
        END IF
        CALL generate_grid(dz_inner,nz_inner,domain_z,alpha_z,dzmin,dzmax,n0z,kv) 
        !
      CASE DEFAULT
        CONTINUE
      END SELECT
      !
      ! ... Set the grid
      ! ... and add new cells if prescribed
      !
      IF (nmx /= 0) dx(1:nmx) = dxmax
      dx(nmx+1:nmx+nx_inner) = dx_inner(1:nx_inner)
      IF (npx /= 0) dx(nx-npx+1:nx) = dxmax
      alpha_x = (alpha_x*domain_x + nmx*dxmax)/(domain_x + (nmx+npx)*dxmax)
      domain_x = domain_x + (nmx+npx)*dxmax
      iv = iv + nmx
      !
      IF (nmy /= 0) dy(1:nmy) = dymax
      dy(nmy+1:nmy+ny_inner) = dy_inner(1:ny_inner)
      IF (npy /= 0) dy(ny-npy+1:ny) = dymax
      alpha_y = (alpha_y*domain_y + nmy*dymax)/(domain_y + (nmy+npy)*dymax)
      domain_y = domain_y + (nmy+npy)*dymax
      jv = jv + nmy
      !
      IF (nmz /= 0) dz(1:nmz) = dzmax
      dz(nmz+1:nmz+nz_inner) = dz_inner(1:nz_inner)
      IF (npz /= 0) dz(nz-npz+1:nz) = dzmax
      alpha_z = (alpha_z*domain_z + nmz*dzmax)/(domain_z + (nmz+npz)*dzmax)
      domain_z = domain_z + (nmz+npz)*dzmax
      kv = kv + nmz
!
      IF (remap_grid) THEN
        ALLOCATE(new_ijk(nx_inner*ny_inner*nz_inner))
        CALL grid_remap(new_ijk)
      END IF
!
! ... Compute grid point (x,y,z) locations (without geographic referencing)
!
      xb(1)=0.D0
      x(1)=-0.5D0*dx(1)+xb(1)
      DO i=2,nx
        xb(i)=(xb(i-1)+dx(i))
        x(i)=(xb(i)-0.5D0*dx(i))
      END DO

      yb(1) = 0.D0
      y(1)=-0.5D0*dy(1)+yb(1)
      DO j=2,ny
        yb(j) = yb(j-1) + dy(j)
        y(j)=(yb(j)-0.5D0*dy(j))
      END DO

      zb(1) = zzero
      z(1)=-0.5D0*dz(1)+zb(1)
      DO k=2,nz
        zb(k) = zb(k-1) + dz(k)
        z(k)=(zb(k)-0.5D0*dz(k))
      END DO

! ... inverse of the cell sizes

      indx = 1.D0 / dx
      indz = 1.D0 / dz

      IF (job_type == '3D') THEN
        indy = 1.D0 / dy
      ELSE IF (job_type == '2D') THEN
        indy = 0.D0
      END IF
!
! ... Set "r" and "rb" absolute coordinate
! ... (only for 2D cylindrical coordinates)
!
      IF((itc == 0) .OR. (job_type == '3D')) THEN

        r(:) = 1.D0
        rb(:) = 1.D0
        inr(:) = 1.D0
        inrb(:) = 1.D0

      ELSE IF ((itc == 1) .AND. (job_type == '2D')) THEN

        r(:)  = x(:)
        rb(:) = xb(:)

        IF (rb(1) == 0.0D0) THEN
          inrb(1) = VERYBIG
        ELSE
          inrb(1) = 1.0D0/rb(1)
        END IF
        DO i=2,nx
          inrb(i)=1.D0/rb(i)
        END DO  
        DO i=1,nx
          inr(i)=1.D0/r(i)
        END DO  

      END IF
!
! ... Set cell flags
!
      CALL flic
!      
      IF( prog == 'PDAC' .AND. mpime == root) THEN
        OPEN( tempunit, FILE='mesh.dat')
            WRITE(tempunit,*) 'Georeferenced x-y mesh'
            WRITE(tempunit,*) 'x'
            WRITE(tempunit,17) x
            WRITE(tempunit,*) 'xb'
            WRITE(tempunit,17) xb
            WRITE(tempunit,*) 'y'
            WRITE(tempunit,17) y
            WRITE(tempunit,*) 'yb'
            WRITE(tempunit,17) yb
            WRITE(tempunit,*) 'z' 
            WRITE(tempunit,17) z
            WRITE(tempunit,*) 'zb' 
            WRITE(tempunit,17) zb
        CLOSE( tempunit )
      END IF

 17   FORMAT(5(F20.6))
      DEALLOCATE(dx_inner, dy_inner, dz_inner)

      RETURN 
      END SUBROUTINE grid_setup
!----------------------------------------------------------------------
      SUBROUTINE flic
      USE dimensions, ONLY: nx, ny, nz, no
      USE control_flags, ONLY: job_type
      USE io_files, ONLY: logunit
      USE parallel, ONLY: mpime, root
!
      IMPLICIT NONE
      INTEGER :: i, j, k, n, ijk
!
! ... Default cell type corresponds to fluid flow
!
      fl = fluid
!
      IF( job_type == '2D' ) THEN
        !
        ! ... Specify cell flags on mesh boundaries
        !
        DO i = 2, nx-1
          k = 1
          ijk = i + (k-1) * nx
          fl(ijk) = 2**(bottom-1)
          k = nz
          ijk = i + (k-1) * nx
          fl(ijk) = 2**(top-1)
        END DO
        !
        DO k = 1, nz
          i = 1
          ijk = i + (k-1) * nx
          fl(ijk) = 2**(west-1)
          i = nx
          ijk = i + (k-1) * nx
          fl(ijk) = 2**(east-1)
        END DO
        !
        ! ... Specify cell flags for blocks
        !
        DO n = 1, no
          DO k = iob(n)%zlo, iob(n)%zhi 
            DO i = iob(n)%xlo, iob(n)%xhi
              ijk = i + (k-1) * nx

              fl(ijk) = 2**(iob(n)%typ-1)
              
            END DO
          END DO
        END DO
      ELSE IF( job_type == '3D' ) THEN
        !
        ! ... Specify cell flags on mesh boundaries
        !
        ! ... Bottom and Top faces
        DO i = 2, nx-1
          DO j = 2, ny-1
            k = 1
            ijk = i + ( (j-1) + (k-1)*ny ) * nx
            fl(ijk) = 2**(bottom-1)
            k = nz
            ijk = i + ( (j-1) + (k-1)*ny ) * nx
            fl(ijk) = 2**(top-1)
          END DO
        END DO
        !
        ! ... West and East faces
        DO j = 1, ny
          DO k = 1, nz
            i = 1
            ijk = i + ( (j-1) + (k-1)*ny ) * nx
            fl(ijk) = 2**(west-1)
            i = nx
            ijk = i + ( (j-1) + (k-1)*ny ) * nx
            fl(ijk) = 2**(east-1)
          END DO
        END DO
        !
        ! ... South and North faces
        DO k = 1, nz
          DO i = 1, nx
            j = 1
            ijk = i + ( (j-1) + (k-1)*ny ) * nx
            fl(ijk) = 2**(south-1)
            j = ny
            ijk = i + ( (j-1) + (k-1)*ny ) * nx
            fl(ijk) = 2**(north-1)
          END DO
        END DO
        !
        ! ... Specify cell flags for blocks
        !
        DO n = 1, no

          DO k = iob(n)%zlo, iob(n)%zhi
            DO j = iob(n)%ylo, iob(n)%yhi
              DO i = iob(n)%xlo, iob(n)%xhi
                ijk = i + ( (j-1) + (k-1)*ny ) * nx

                fl(ijk) = 2**(iob(n)%typ-1)
              
              END DO
            END DO
          END DO
        END DO
        !
      END IF
      !IF (mpime == root) WRITE(logunit,*) fl
!
      RETURN
      END SUBROUTINE flic
!----------------------------------------------------------------------
      SUBROUTINE generate_grid(delta,nd,domain_size,alpha,demin,demax,n0,center)
      USE control_flags, ONLY: lpr
      USE parallel, ONLY: mpime, root
!
! ... this routine generates a 1D rectilinear (non)uniform mesh
! ... given the number of cells, the minimum and maximum size,
! ... and the size increase rate. Since not every set of values
! ... comes out to be consistent, priority is given to 
! ... 1) the number of cells; 2) the domain size; 3) the minimum cell
! ... size; 4) the maximum cell size; 5) the maximum increase rate 
! ... constraints.
!
      IMPLICIT NONE
!
! ... array of cell sizes
!
      REAL*8, DIMENSION(:), INTENT(OUT) :: delta
!
! ... number of cells 
!
      INTEGER, INTENT(IN)  :: nd
!
! ... size of computational domain
!
      REAL*8, INTENT(INOUT) ::  domain_size
!
! ... relative center of the mesh (for refinement)
!
      REAL*8, INTENT(INOUT) :: alpha
!
! ... minimum and maximum mesh size
!
      REAL*8, INTENT(IN) ::  demin
      REAL*8, INTENT(INOUT) ::  demax
!
! ... number of cells with minimum size beyond the center
!
      INTEGER, INTENT(IN), OPTIONAL  :: n0
!
! ... index of the central cell
!
      INTEGER, INTENT(OUT) :: center
!
      REAL*8 :: ded, der, dem, lder
!      
      REAL*8 ::  minbeta = 1.D0 + 1.D-5
      REAL*8 ::  dbeta0 = 0.99
      REAL*8  :: beta, beta1, beta2, frac, dbeta
      REAL*8  :: ntilde1, ntilde2
      REAL*8  :: l, l1, l2, lcomp
      INTEGER :: n01,n02,m1,m2,n11,n12
      INTEGER :: i, j, m, n
      LOGICAL :: print_mesh, already,maxn, ionode
      INTEGER :: idelta, ncount
!
      ionode = (mpime == root)
!
      IF (demin == demax) THEN
              delta = demin
!              alpha = 0.5D0
              center = nd/2
              domain_size = nd * demin
              RETURN
      END IF
      IF ( domain_size/demin < nd )  CALL error('grid_generator', &
          'insufficient number of cells', nd) 
      IF ( domain_size/demin == nd ) CALL error('grid_generator', &
          'number of cells is too big', nd)
      IF (domain_size*alpha < demin) alpha = 0.D0
      IF (domain_size*(1.D0-alpha) < demin) alpha = 1.D0
!
      n01 = n0 / 2
      IF( MOD(nd,2) == 0) n01 = n01 + 1
      n02 = n0 - n01
!
! ... 'l1' and 'l2' are the sizes of the left and right domain
! ... minus the uniform (fine) mesh
!
      l1 = (domain_size*alpha) - n01*demin - 0.5D0 * demin
      l2 = (domain_size-alpha*domain_size) - n02*demin - 0.5D0 * demin
      
      IF (alpha == 0.D0) THEN
  
        n01 = 0
        n02 = n0

        l1 = 0.D0
        l2 = domain_size - n02*demin - demin
        
      ELSE IF (alpha == 1.D0) THEN

        n01 = n0
        n02 = 0

        l1 = domain_size - n01*demin - demin
        l2 = 0.D0
        
      END IF
!
      der = demax / demin
      beta = 1.2D0
      dbeta = dbeta0
      IF( lpr > 2 .AND. ionode) WRITE(testunit,*) 'Initial beta = ',beta
!
! .. loop over the decreasing 'beta' to fit the domain size
!
      ncount = 1
      print_mesh = .FALSE.
      already = .FALSE.

      DO WHILE (.NOT.print_mesh)
         
        ! ... 'm' is the number of cells needed to increase the
        ! ... cell size from demin to demax at constant rate 'beta'
        !
        m = INT(LOG(der)/LOG(beta))

        ! ... 'lcomp' is the size of the non-uniform part of the mesh
        !
         frac  = ( beta**(m+1) - 1.D0 ) / ( beta - 1.D0 )
         lcomp = demin * ( frac - 1.D0 )
!
         IF (l1 < 0) THEN 
            n01 = alpha / demin
            l1  = 0.D0
            m1  = 0
            n11 = 0
         ELSE
            IF ( (l1 / demax) > nd - n01 -1 ) &
              CALL error('grid', 'number of cells is too small!', nd)
            IF (lcomp < l1) THEN
               m1  = m
               n11 = INT( (l1 - lcomp) / demax )
            ELSE
               m1  = INT( LOG( 1.D0 + (beta-1.D0)*(l1/demin+1.D0) ) / LOG(beta) - 1 )
               n11 = 0
            ENDIF
         ENDIF
!        
         IF (l2 < 0) THEN 
            n02 = (1.D0-alpha) / demin
            l2  = 0.D0
            m2  = 0
            n12 = 0
         ELSE
            IF ( l2 / demax > nd - n02 - 1 ) &
              CALL error('grid', 'number of cells is too small!', nd)
            IF (lcomp < l2) THEN
               m2  = m
               n12 = INT( ( l2 - lcomp ) / demax )
            ELSE
               m2  = INT( LOG( 1.D0 + (beta-1.D0)*(l2/demin+1.D0) ) / LOG(beta) - 1 )
               n12 = 0
            ENDIF
         ENDIF 
!
        IF ( (l2+l1) / demax > nd - n01 - n02 - 1 ) &
          CALL error('grid', 'number of cells is too small!', nd)
!
        IF ( (n01+m1+n11)+(n02+m2+n12)+1 < nd ) THEN 
           IF ( beta*dbeta < minbeta) THEN
              print_mesh = .TRUE.
           ELSE
              IF ( already ) THEN
                 dbeta = dbeta + ( 1.D0 - dbeta) / 2.D0
                 already = .FALSE.
              ENDIF
              beta = beta * dbeta 
              IF( lpr > 2 .AND. ionode) &
                 WRITE(testunit,*) 'Reducing beta = ', beta, ' n = ', &
              (n01+m1+n11)+(n02+m2+n12)+1 
           ENDIF
        ELSE IF ( (n01+m1+n11)+(n02+m2+n12)+1 > nd ) THEN
           IF ( .not.already ) THEN
              dbeta = dbeta + (1.0 - dbeta) / 2.0
              already = .TRUE.
           ENDIF
           beta = beta / dbeta
           IF( lpr > 2 .AND. ionode ) &
             WRITE(testunit,*) 'Increasing beta = ', beta, ' n = ', &
             (n01+m1+n11)+(n02+m2+n12)+1
        ELSE 
           IF ( (1.D0-dbeta) <= 1.D-15 ) THEN
              print_mesh= .TRUE.
              IF( lpr > 0 .AND. ionode) &
                WRITE(testunit,*) 'Final beta = ', beta, ' n= ', (n01+m1+n11)+(n02+m2+n12)+1
           ELSE
             already = .TRUE.
             beta = beta / dbeta
             IF( lpr > 2 .AND. ionode) &
               WRITE(testunit,*) 'Increasing beta = ', beta, ' n = ', (n01+m1+n11)+(n02+m2+n12)+1
           ENDIF
        ENDIF
      END DO
!
      IF ( n11+n12 == 0 ) THEN
        IF (ionode) THEN
          WRITE(errorunit,*) 'WARNING!!: no cells with maximum size!'
          WRITE(errorunit,*) 'Please decrease beta or dmax'
        END IF
      ENDIF 
!    
      IF ( beta >= maxbeta ) THEN
         IF (ionode) THEN
           WRITE(errorunit,*) 'WARNING!!: beta >= maximum beta!'
           WRITE(errorunit,*) 'Please increase maximum beta or number of cells'
         END IF
      ENDIF
!
      maxn = .TRUE.
      DO WHILE ( maxn )
         ntilde1 = ( l1 - n11*demax ) / demin + 1.D0
         beta1 = beta
         IF ( ntilde1 > 1.D0 ) THEN
            DO j=1,50
               beta1 = ( 1.D0 - ( m1*beta1**( m1+1 ) + 1.D0 ) / ntilde1 ) / &
               &      ( 1.D0 - ( m1+1.D0 )*( beta1**m1 ) / ntilde1 ) 
            ENDDO
         ENDIF
         IF ( demin*beta1**m1 > demax ) THEN

            IF ( (m1==0) .OR. ( beta1 > maxbeta ) ) THEN
               IF (ionode) WRITE(errorunit,*) 'WARNING!!: beta >= maximum beta!'
               maxn = .FALSE.
            ELSE
               n11 = n11 + 1
               m1 = m1 - 1
            ENDIF
         ELSE
            maxn = .FALSE.
         ENDIF
      END DO
!     
      maxn = .TRUE.
      DO WHILE ( maxn )
         ntilde2 = ( l2 - n12*demax) / demin + 1.D0
         beta2 = beta
         IF ( ntilde2 > 1.D0 ) THEN
            DO j=1,50
               beta2 = ( 1.D0 - ( m2*beta2**( m2+1.D0 ) + 1.D0 ) / ntilde2 ) / &
               ( 1.D0 - ( m2+1.D0 )*( beta2**m2 ) / ntilde2 )
            ENDDO
         ENDIF
         IF ( demin*beta2**m2 > demax ) THEN
            IF ( (m2==0) .OR. ( beta2 > maxbeta ) ) THEN
               IF(ionode) WRITE(errorunit,*) 'WARNING!!: beta >= maximum beta!'
               maxn = .FALSE.
            ELSE
               n12 = n12 + 1
               m2 = m2 - 1
            ENDIF
         ELSE
            maxn = .FALSE.
         ENDIF
      END DO          
!
      IF( lpr > 0 .AND. ionode) THEN
        WRITE(testunit,*) 'demin, demax : ', demin, demax
        WRITE(testunit,*) 'n01, m1, n11 : ', n01, m1, n11
        WRITE(testunit,*) 'n02, m2, n12 : ', n02, m2, n12
      END IF
!
      center = n11 + m1 + n01 + 1

      delta(1:n11) = demax
      DO i = 1, m1
        delta(center-n01-i) = demin * beta1**i
      END DO
      delta(center-n01:center-1) = demin
      delta(center) = demin
      delta(center+1:center+n02) = demin
      DO i = 1, m2
        delta(center+n02+i) = demin * beta2**i
      END DO
      delta(center+n02+m2+1:) = demax
!
! ... accuracy of the mesh down to centimetres (second decimal digit)
!
      DO i = 1, nd
        delta(i) = delta(i) * 10000.D0
        idelta = NINT(delta(i))
        delta(i) = idelta / 10000.D0
      END DO
      
      IF (lpr > 0 .AND. ionode) THEN
        WRITE(testunit,777) domain_size
        WRITE(testunit,888) SUM(delta)
      END IF

 777  FORMAT('domain_size = ',(F8.2)) 
 888  FORMAT('mesh_size   = ',(F8.2)) 
!
      RETURN
      END SUBROUTINE generate_grid
!----------------------------------------------------------------------
      END MODULE grid
!----------------------------------------------------------------------
