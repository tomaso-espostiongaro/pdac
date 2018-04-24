!----------------------------------------------------------------------
      MODULE compute_mean_fields
!----------------------------------------------------------------------
      IMPLICIT NONE
!
      INTEGER :: imrt
      REAL*8, DIMENSION(:), ALLOCATABLE :: rhom_gav, um_gav, vm_gav, wm_gav, tm_gav
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: yd_gav
!
      PRIVATE
      PUBLIC :: imrt, rhom_gav, um_gav, vm_gav, wm_gav, tm_gav, yd_gav, &
                average_mixture_fields, allocate_mean_fields
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE allocate_mean_fields
      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE dimensions
      USE domain_mapping, ONLY: ncdom, ncint
      ! ... Time-averaged output fields
!
      IMPLICIT NONE
!
      ALLOCATE(rhom_gav(ncdom))  ! Mixture Density
      ALLOCATE(tm_gav(ncdom))    ! Mixture Temperature
      ALLOCATE(um_gav(ncdom))    ! Mixture Velocity X
      IF (job_type == JOB_TYPE_3D) ALLOCATE(vm_gav(ncdom))    ! Mixture Velocity Y
      ALLOCATE(wm_gav(ncdom))    ! Mixture Velocity Z
      ALLOCATE(yd_gav(ncdom,ngas+nsolid)) 
      rhom_gav = 0.0D0
      tm_gav   = 0.0D0
      um_gav   = 0.0D0
      IF (job_type == JOB_TYPE_3D) vm_gav   = 0.0D0
      wm_gav   = 0.0D0
      yd_gav   = 0.0D0
!
      RETURN
      END SUBROUTINE allocate_mean_fields
!----------------------------------------------------------------------
      SUBROUTINE average_mixture_fields(time, dt, timestart)
!
      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE mixture_fields, ONLY: rhom, um, vm, wm, tm, yd
      USE io_files, ONLY: logunit
      USE parallel, ONLY: mpime, root
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: time, dt, timestart
      REAL*8 :: told, tnew, intnew
!
      told = time - timestart - dt
      tnew = time - timestart
      intnew = 1.D0 / tnew
!
      rhom_gav = intnew * (rhom_gav * told + rhom * dt) 
      um_gav = intnew * (um_gav * told + um * dt)
      IF (job_type == JOB_TYPE_3D) vm_gav = intnew * (vm_gav * told + vm * dt)
      wm_gav = intnew * (wm_gav * told + wm * dt)
      tm_gav = intnew * (tm_gav * told + tm * dt)
      yd_gav = intnew * (yd_gav * told + yd * dt)
!
      RETURN
      END SUBROUTINE average_mixture_fields
!----------------------------------------------------------------------
      END MODULE compute_mean_fields 
!----------------------------------------------------------------------
