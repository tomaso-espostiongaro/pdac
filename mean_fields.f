!----------------------------------------------------------------------
      MODULE compute_mean_fields 
!----------------------------------------------------------------------
      IMPLICIT NONE
!
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!----------------------------------------------------------------------
      SUBROUTINE compute_time_averages
!
      USE control_flags, ONLY: job_type
      USE control_flags, ONLY: JOB_TYPE_2D, JOB_TYPE_3D
      USE io_files, ONLY: logunit
      USE parallel, ONLY: mpime, root
      USE postp_variables, ONLY: rhom, um, vm, wm, tm
      USE postp_variables, ONLY: rhom_gav, um_gav, vm_gav, wm_gav, tm_gav

      IMPLICIT NONE
!
      rhom_gav = rhom_gav + rhom 
      um_gav = um_gav + um
      IF (job_type == JOB_TYPE_3D) vm_gav = vm_gav + vm
      wm_gav = wm_gav + wm
      tm_gav = tm_gav + tm
!
      RETURN
      END SUBROUTINE compute_time_averages
!----------------------------------------------------------------------
      END MODULE compute_mean_fields 
!----------------------------------------------------------------------
