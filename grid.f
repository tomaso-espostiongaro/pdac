!----------------------------------------------------------------------
      MODULE grid
!
! ... This module contains all the information about the rectilinear
! ... 2D or 3D mesh used for computation, and the cell types 
!
!----------------------------------------------------------------------
      IMPLICIT NONE
!
      REAL*8, DIMENSION(:), ALLOCATABLE :: dx, dy, dz
      REAL*8, DIMENSION(:), ALLOCATABLE :: dxtemp
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
!
! ... flags for the domain boundary condition
      INTEGER :: west, east, south, north, bottom, top
!
! ... file name for the topography
      CHARACTER(LEN=80) :: topography
!
! ... origin of atmospheric stratification
      REAL*8 :: zzero
!
! ... variables for grid generator:
! ... maximum cell size increase rate
      REAL*8 ::  maxbeta
      REAL*8 ::  minbeta = 1.D0 + 1.D-5
      REAL*8 ::  dbeta0 = 0.99
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
! ... flag for increasing cell sizes: 1-constant rate; 2-constant slope
      INTEGER :: grigen
!
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE allocate_grid
!
      USE dimensions, ONLY: nx, ny, nz, ntot, ntr
      USE control_flags, ONLY: job_type
!
      IMPLICIT NONE
!
! ... Set the appropriate total number of cells
! ... and the coordinate system

      IF( job_type == '2D' ) THEN
        ntot = nx*nz
        ntr  = nx
      ELSE IF( job_type == '3D' ) THEN
        ntot = nx*ny*nz
        ntr  = nx*ny
      ELSE
        CALL error( ' allocate_grid ', ' wrong job_type '//job_type, 1)
      END IF

      ALLOCATE( dxtemp(nx) ) 
      ALLOCATE( dx(nx), dy(ny), dz(nz) ) 
      ALLOCATE( indx(nx), indy(ny), indz(nz) )
      ALLOCATE( r(nx), rb(nx) )
      ALLOCATE( inr(nx), inrb(nx) )
      ALLOCATE( x(nx), xb(nx) )
      ALLOCATE( y(ny), yb(ny) )
      ALLOCATE( z(nz), zb(nz) )
      ALLOCATE( fl(ntot) )

      RETURN
      END SUBROUTINE allocate_grid
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
      USE control_flags, ONLY: job_type
      USE dimensions, ONLY: nx, ny, nz

      REAL*8, PARAMETER :: VERYBIG = 1.0d+10
      INTEGER :: i, j, k, ijk
      REAL*8 :: zrif
!
! ... generate the non-uniform mesh (without referencing)
!
      IF (grigen > 0) THEN

        WRITE(7,*) 'Generate grid along x'
        CALL generate_grid(dx,nx,domain_x,alpha_x,dxmin,dxmax,n0x,iv) 

        IF (job_type == '3D') THEN

          WRITE(7,*) 'Generate grid along y'
          CALL generate_grid(dy,ny,domain_y,alpha_y,dymin,dymax,n0y,jv) 
        END IF

        WRITE(7,*) 'Generate grid along z'
        CALL generate_grid(dz,nz,domain_z,alpha_z,dzmin,dzmax,n0z,kv) 

      END IF
!
! ... Compute grid point (x,y,z) locations
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
      RETURN 
      END SUBROUTINE grid_setup
!----------------------------------------------------------------------
      SUBROUTINE flic
      USE dimensions, ONLY: nx, ny, nz, no
      USE control_flags, ONLY: job_type
!
      IMPLICIT NONE
      INTEGER :: i, j, k, n, ijk
!
! ... Default cell type corresponds to fluid flow
!
      fl = 1
!
      IF( job_type == '2D' ) THEN
        !
        ! ... Specify cell flags on mesh boundaries
        !
        DO i = 1, nx
          k = 1
          ijk = i + (k-1) * nx
          fl(ijk) = bottom
          k = nz
          ijk = i + (k-1) * nx
          fl(ijk) = top
        END DO
        !
        DO k = 1, nz
          i = 1
          ijk = i + (k-1) * nx
          fl(ijk) = west
          i = nx
          ijk = i + (k-1) * nx
          fl(ijk) = east
        END DO
        !
        ! ... Specify cell flags for blocks
        !
        DO n = 1, no
          DO k = iob(n)%zlo, iob(n)%zhi 
            DO i = iob(n)%xlo, iob(n)%xhi
              ijk = i + (k-1) * nx

              fl(ijk) = iob(n)%typ
              
            END DO
          END DO
        END DO
      ELSE IF( job_type == '3D' ) THEN
        !
        ! ... Specify cell flags on mesh boundaries
        !
        DO i = 1, nx
          DO j = 1, ny
            k = 1
            ijk = i + ( (j-1) + (k-1)*ny ) * nx
            fl(ijk) = bottom
            k = nz
            ijk = i + ( (j-1) + (k-1)*ny ) * nx
            fl(ijk) = top
          END DO
        END DO
        !
        DO j = 1, ny
          DO k = 1, nz
            i = 1
            ijk = i + ( (j-1) + (k-1)*ny ) * nx
            fl(ijk) = west
            i = nx
            ijk = i + ( (j-1) + (k-1)*ny ) * nx
            fl(ijk) = east
          END DO
        END DO
        !
        DO k = 1, nz
          DO i = 1, nx
            j = 1
            ijk = i + ( (j-1) + (k-1)*ny ) * nx
            fl(ijk) = south
            j = ny
            ijk = i + ( (j-1) + (k-1)*ny ) * nx
            fl(ijk) = north
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

                fl(ijk) = iob(n)%typ
              
              END DO
            END DO
          END DO
        END DO
        !
      END IF
!
! ... Set flags on corners
!
      IF( job_type == '2D' ) THEN
        IF ( fl(nz*nx - 1) == 4 .AND. fl((nz-1)*nx) == 4) THEN
          fl(nz*nx) = 4
        END IF
      END IF
      
      RETURN
      END SUBROUTINE flic
!----------------------------------------------------------------------
      SUBROUTINE generate_grid(delta,nd,domain_size,alpha,demin,demax,n0,center)
      USE control_flags, ONLY: lpr
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
      REAL*8, INTENT(IN) ::  domain_size
!
! ... relative center of the mesh (for refinement)
!
      REAL*8 :: alpha
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
      REAL*8 :: ded, der, dem, lder
!      
      REAL*8  :: beta, beta1, beta2, frac, dbeta
      REAL*8  :: ntilde1, ntilde2
      REAL*8  :: l, l1, l2, lcomp
      INTEGER :: n01,n02,m1,m2,n11,n12
      INTEGER :: i, j, m, n, center
      LOGICAL :: print_mesh, already,maxn
      INTEGER :: idelta, ncount
!
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

      IF( lpr > 2 ) WRITE(7,*) 'n01,n02,l1,l2',n01,n02,l1,l2
!
      der = demax / demin
      beta = 1.2D0
      dbeta = dbeta0
      IF( lpr > 2 ) WRITE(7,*) 'Initial beta = ',beta
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
            IF ( (l1 / demax) > nd - n01 -1 ) CALL error('grid', 'number of cells is too small!', nd)
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
            IF ( l2 / demax > nd - n02 - 1 ) CALL error('grid', 'number of cells is too small!', nd)
            IF (lcomp < l2) THEN
               m2  = m
               n12 = INT( ( l2 - lcomp ) / demax )
            ELSE
               m2  = INT( LOG( 1.D0 + (beta-1.D0)*(l2/demin+1.D0) ) / LOG(beta) - 1 )
               n12 = 0
            ENDIF
         ENDIF 
!
        IF ( n11+n12 == 0 ) THEN
          WRITE(7,*) 'WARNING!!: no cells with maximum size!'
          WRITE(7,*) 'Please decrease beta or dmax'
        ENDIF 
!    
        IF ( (l2+l1) / demax > nd - n01 - n02 - 1 ) CALL error('grid', 'number of cells &
        & is too small!', nd)
!
        IF ( (n01+m1+n11)+(n02+m2+n12)+1 < nd ) THEN 
           IF ( beta*dbeta < minbeta) THEN
              WRITE(7,*) 'WARNING!!: beta = minimum beta!'
              WRITE(7,*) 'Please decrease minimum beta or number of cells'
              print_mesh = .TRUE.
           ELSE
              IF ( already ) THEN
                 dbeta = dbeta + ( 1.D0 - dbeta) / 2.D0
                 already = .FALSE.
              ENDIF
              beta = beta * dbeta 
              IF( lpr > 2 ) WRITE(7,*) 'Reducing beta = ', beta, ' n = ', &
              (n01+m1+n11)+(n02+m2+n12)+1 
           ENDIF
        ELSE IF ( (n01+m1+n11)+(n02+m2+n12)+1 > nd ) THEN
           IF ( .not.already ) THEN
              dbeta = dbeta + (1.0 - dbeta) / 2.0
              already = .TRUE.
           ENDIF
           beta = beta / dbeta
           IF( lpr > 2 ) WRITE(7,*) 'Increasing beta = ', beta, ' n = ', &
             (n01+m1+n11)+(n02+m2+n12)+1
        ELSE 
           IF ( (1.D0-dbeta) <= 1.D-15 ) THEN
              print_mesh= .TRUE.
              IF( lpr > 2 ) WRITE(7,*) 'Final beta = ', beta, ' n= ', &
              (n01+m1+n11)+(n02+m2+n12)+1
           ELSE
             already = .TRUE.
             beta = beta / dbeta
             IF( lpr > 2 ) WRITE(7,*) 'Increasing beta = ', beta, ' n = ', &
              (n01+m1+n11)+(n02+m2+n12)+1
           ENDIF
        ENDIF
        
        CALL myflush( 7 )
      END DO
!
      IF ( beta >= maxbeta ) THEN
         WRITE(7,*) 'WARNING!!: beta >= maximum beta!'
         WRITE(7,*) 'Please increase maximum beta or number of cells'
      ENDIF
!

      WRITE(6,*) 'check!'

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

            WRITE(6,*) m1, demin*beta1**m1

            IF ( (m1==0) .OR. ( beta1 > maxbeta ) ) THEN
               WRITE(7,*) 'WARNING!!: beta >= maximum beta!'
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
               WRITE(7,*) 'WARNING!!: beta >= maximum beta!'
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
      IF( lpr >= 1 ) WRITE(7,*) 'beta, beta1, beta2: ', beta, beta1, beta2 
      IF( lpr >= 1 ) WRITE(7,*) 'n01,m1,n11',n01,m1,n11
      IF( lpr >= 1 ) WRITE(7,*) 'n02,m2,n12',n02,m2,n12
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
        delta(i) = delta(i) * 100.D0
        idelta = NINT(delta(i))
        delta(i) = idelta / 100.D0
        WRITE(7,'(F8.2)') delta(i)
      END DO
      
      WRITE(7,777) domain_size
      WRITE(7,888) SUM(delta)

 777  FORMAT('domain_size = ',(F8.2)) 
 888  FORMAT('mesh_size = ',(F8.2)) 
!
      RETURN
      END SUBROUTINE generate_grid
!----------------------------------------------------------------------
      END MODULE grid
!----------------------------------------------------------------------
