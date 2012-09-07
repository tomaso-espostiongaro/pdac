!--------------------------------------------------------------
      MODULE mass_sink
!--------------------------------------------------------------
! ... This module contains the sink term.
!
      USE set_indexes, ONLY: stencil
      IMPLICIT NONE
!
      REAL*8, DIMENSION(:,:), ALLOCATABLE:: sink     ! Array of instantaneous sink terms  
      REAL*8, DIMENSION(:,:), ALLOCATABLE:: massl, mass_loss ! Array of cumulative (in time) sink terms
      INTEGER :: isink
      REAL*8 :: rprox, rprox2
!
      PUBLIC :: sink_update, isink, rprox, sink, print_mass_loss, allocate_sink, read_mass_loss
      PRIVATE
!
      SAVE
! 
!--------------------------------------------------------------
      CONTAINS
!--------------------------------------------------------------
      SUBROUTINE allocate_sink
!
! ... This subroutine allocates and initialises the sink array,
! ... defined in all the computational domain, and set equal to 0.
! ... It'll be updated in bdry.f, when the flag of the bottom cell
! ... (n2=ijkm) will be considered.
!
      USE dimensions, ONLY: nsolid, ntr
      USE domain_mapping, ONLY: ncdom
      USE kinds
!
      ALLOCATE (sink(ncdom,nsolid))
      ALLOCATE (mass_loss(ntr, nsolid))
      ALLOCATE (massl(ntr, nsolid))
      sink = 0.D0
      mass_loss = 0.D0
      massl = 0.D0
!
      rprox2 = rprox**2
!
      END SUBROUTINE
!--------------------------------------------------------------
      SUBROUTINE sink_update(ijk,is)
!
      USE control_flags, ONLY: job_type, JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions, ONLY: nsolid, nx
      USE domain_mapping, ONLY: meshinds
      USE gas_solid_density, ONLY: rlk
      USE gas_solid_velocity, ONLY: ws
      USE grid, ONLY: flag, indz, x, y, dx, dy, dz
      USE grid, ONLY: center_x, center_y
      USE interpolate_fields, ONLY: wsink
      USE io_files, ONLY: testunit
      USE particles_constants, ONLY: inrl
      USE pressure_epsilon, ONLY: ep
      USE set_indexes, ONLY: subscr, ijkp, ijkm
      USE time_parameters, ONLY: sweep, nprint, dt
!
      INTEGER, INTENT(IN) :: ijk, is
      REAL*8 :: rsft_up, rsft_down, dist2
      REAL*8 :: wst, epsmax, eps, fact, nexp
      INTEGER :: i,j,k,imesh, bind
!
! ... Allocate and initialize local mass fluxes.
!
      rsft_up = 0.0D0
      rsft_down = 0.0D0
!
      CALL subscr(ijk)
      CALL meshinds(ijk, imesh, i, j, k)
!
! ......................................................................
! ... Insert here the sink term for the particle mass transport equation
!
      rsft_up = rlk(ijk,is)
      rsft_down = rlk(ijkm,is)
!
      epsmax = 0.3  ! Maximum packing
      nexp   = 4.65  ! Richardson and Zaki (1954) exponent
!
      eps = rlk(ijk,is) * inrl(is)
      fact = epsmax/eps * ((1.D0 - epsmax)/(1.D0 - eps))**nexp
      fact = MIN(fact,1.D0)
!      
      wst  = wsink(ijk,is) * fact
      wst = wsink(ijk,is)
!
! ... Sink term as in Cordoba et al., JVGR, Vol. 139 (2005)
!      sink(ijk,is) = DABS( (rsft_up*wsink(ijk,is)-rsft_down*wsink(ijkm,is))*indz(k))
!
! ... Sink term as a convective flux
      dist2 = (x(i)-center_x)**2 + (y(j)-center_y)**2
      IF (dist2 > rprox2) sink(ijk,is) = rsft_up * wst * indz(k)
! ......................................................................
!
! ... The total mass lost from the bottom is updated
!
      IF (job_type == JOB_TYPE_2D) THEN
        bind = i
        massl(bind,is) = massl(bind,is) + sink(ijk,is)*dx(i)*dz(k) * dt
      ELSE IF (job_type == JOB_TYPE_3D) THEN
        bind = i + nx*(j-1)
        massl(bind,is) = massl(bind,is) + sink(ijk,is)*dx(i)*dy(j)*dz(k) * dt
      END IF
!
! ... On restart, the massl (local) array must be incremented by the total mass-loss value
! ... Every processor adds the values on its local array positions
! ... If this is not a restart, mass_loss equals zero.
! 
      massl(bind,is) = mass_loss(bind,is) +  massl(bind,is)
      mass_loss(bind,is) = 0.D0
!
      RETURN
      END SUBROUTINE sink_update
!--------------------------------------------------------------
      SUBROUTINE print_mass_loss(nou)
      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions, ONLY: ntr, nsolid
      USE io_files, ONLY: masslunit
      USE parallel, ONLY: mpime, root
      IMPLICIT NONE
      CHARACTER*4 :: lettera
      CHARACTER(LEN=10) :: filename
      INTEGER :: is, i
      REAL*8 :: masstotal
      INTEGER, OPTIONAL, INTENT(IN) :: nou
!
      mass_loss = massl
      CALL parallel_sum_real(mass_loss,ntr*nsolid)
!
      IF (PRESENT(nou)) THEN
        filename = 'massl.'//lettera(nou)
      ELSE 
        filename = 'massl.res'
      END IF
      IF (mpime == root) THEN
        OPEN(masslunit,FILE=filename,STATUS='UNKNOWN')
        IF (job_type == JOB_TYPE_2D) THEN
          masstotal = 0.D0
          DO i = 1, ntr
            masstotal = masstotal + SUM(mass_loss(i,:))
            WRITE(masslunit,122) (mass_loss(i,is), is=1, nsolid), masstotal
          END DO
        ELSE IF (job_type == JOB_TYPE_3D) THEN
          DO is = 1, nsolid
            masstotal = SUM(mass_loss(:,is))
            WRITE(masslunit,'(1x,//,1x,"ML",I1,"=",G14.6,/)') is, masstotal
            WRITE(masslunit,122) mass_loss(:,is)
          END DO
        END IF
        CLOSE(masslunit)
      END IF
 122  FORMAT(10(1x,G14.6E3))
!
! ... Avoid accumulation!
      mass_loss = 0.D0
!
      RETURN
      END SUBROUTINE print_mass_loss
!--------------------------------------------------------------
      SUBROUTINE read_mass_loss(nou)
      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions, ONLY: ntr, nsolid
      USE io_files, ONLY: masslunit, logunit
      USE parallel, ONLY: mpime, root
      IMPLICIT NONE
      CHARACTER*4 :: lettera
      CHARACTER(LEN=10) :: filename
      INTEGER :: is, i, iss
      REAL*8 :: masstotal
      INTEGER, OPTIONAL, INTENT(IN) :: nou
!
      IF (PRESENT(nou)) THEN
        filename = 'massl.'//lettera(nou)
      ELSE 
        filename = 'massl.res'
      END IF
      IF (mpime == root) THEN
        WRITE(logunit,fmt="('  from sink: recovering file ',A20)") filename
        OPEN(masslunit,FILE=TRIM(filename),STATUS='UNKNOWN')
        IF (job_type == JOB_TYPE_2D) THEN
          DO i = 1, ntr
            READ(masslunit,122) (mass_loss(i,is), is=1, nsolid), masstotal
          END DO
        ELSE IF (job_type == JOB_TYPE_3D) THEN
          DO is = 1, nsolid
            READ(masslunit,'(1x,//,1x,"ML",I1,"=",G14.6,/)') iss, masstotal
            READ(masslunit,122) mass_loss(:,is)
          END DO
        END IF
        CLOSE(masslunit)
      END IF
 122  FORMAT(10(1x,G14.6E3))
!
      CALL bcast_real(mass_loss,ntr*nsolid,root)
!
      RETURN
      END SUBROUTINE read_mass_loss
!--------------------------------------------------------------
      END MODULE mass_sink
!--------------------------------------------------------------
