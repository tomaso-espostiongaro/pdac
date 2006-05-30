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
      USE io_files, ONLY: logunit
      USE parallel, ONLY: mpime, root
      USE postp_variables, ONLY: rhom, um, vm, wm
      USE postp_variables, ONLY: rhom_av, um_av, vm_av, wm_av

      IMPLICIT NONE
!
      rhom_av = rhom_av + rhom 
      um_av = um_av + um
      IF (job_type == '3D') vm_av = vm_av + vm
      wm_av = wm_av + wm
!
      RETURN
      END SUBROUTINE compute_time_averages
!-----------------------------------------------------------------------
      END MODULE compute_mean_fields 
!----------------------------------------------------------------------
