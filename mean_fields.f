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
      USE postp_variables, ONLY: rhom, um, vm, wm, tm
      USE postp_variables, ONLY: rhom_av, um_av, vm_av, wm_av, tm_av

      IMPLICIT NONE
!
      rhom_av = rhom_av + rhom 
      um_av = um_av + um
      IF (job_type == '3D') vm_av = vm_av + vm
      wm_av = wm_av + wm
      tm_av = tm_av + tm
!
      RETURN
      END SUBROUTINE compute_time_averages
!----------------------------------------------------------------------
      SUBROUTINE compute_vertical_mean_profiles
      USE dimensions
      USE grid, ONLY: dx, dy, dz
      USE domain_mapping, ONLY: ncint, meshinds
      USE postp_variables, ONLY: rhom, um, vm, wm, tm, epst
      USE postp_variables, ONLY: rhom_z, um_z, vm_z, wm_z, tm_z, surface
      IMPLICIT NONE
!
      REAL*8 :: ds, invrhom, invsurf
      INTEGER :: ijk, imesh, i, j, k
!
      invrhom = 0.D0
      invsurf = 0.D0
      ds = 0.D0
!
      DO ijk = 1, ncint
        IF ( epst(ijk) >= 1.D-8) THEN
          CALL meshinds(ijk,imesh,i,j,k)
          ds = dx(i)*dy(j)
          rhom_z(k) = rhom_z(k) + rhom(ijk) * ds
          um_z(k) = um_z(k) + rhom(ijk) * um(ijk) * ds
          vm_z(k) = vm_z(k) + rhom(ijk) * vm(ijk) * ds
          wm_z(k) = wm_z(k) + rhom(ijk) * wm(ijk) * ds
          tm_z(k) = tm_z(k) + tm(ijk) * ds
          surface(k) = surface(k) + ds
        END IF
      END DO
!
      CALL parallel_sum_real(surface,nz)
      CALL parallel_sum_real(rhom_z,nz)
      CALL parallel_sum_real(um_z,nz)
      CALL parallel_sum_real(vm_z,nz)
      CALL parallel_sum_real(wm_z,nz)
      CALL parallel_sum_real(tm_z,nz)
!
! ... Spatial average on the plane
! ... Velocities in each cell are weighted 
! ... with the corresponding mixture density
!
      DO k=1,nz
        IF (surface(k)/=0.D0) invsurf = 1.D0 / surface(k)
        rhom_z(k) = rhom_z(k) * invsurf
        IF (rhom_z(k)/=0.D0) invrhom = 1.D0 / rhom_z(k)
        um_z(k) = um_z(k) * invrhom *  invsurf
        vm_z(k) = vm_z(k) * invrhom *  invsurf
        wm_z(k) = wm_z(k) * invrhom *  invsurf
        tm_z(k) = tm_z(k) * invsurf
      END DO
!
      RETURN
      END SUBROUTINE compute_vertical_mean_profiles
!-----------------------------------------------------------------------
      END MODULE compute_mean_fields 
!----------------------------------------------------------------------
