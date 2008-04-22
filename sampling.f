!----------------------------------------------------------------------
      MODULE sample_points
!----------------------------------------------------------------------
! ... Use ALL independent and derived variables
      USE postp_variables
!
      IMPLICIT NONE
!
! ... parameters for sampling
!
      INTEGER :: number_of_probes
      INTEGER, ALLOCATABLE :: ijk_probe(:), imesh_column_probe(:)
      LOGICAL :: assign_index
      CHARACTER(LEN=80) :: probe_file
!
      TYPE probe_point
            INTEGER :: nop
            INTEGER :: i
            INTEGER :: j
            INTEGER :: k
            REAL*8  :: x
            REAL*8  :: y
            REAL*8  :: z
       END TYPE probe_point
!
      TYPE(probe_point), ALLOCATABLE :: probe(:)
      INTEGER :: isamp
      INTEGER :: iiv, jjv
! 
      SAVE
!----------------------------------------------------------------------
      CONTAINS
!---------------------------------------------------------------------
      SUBROUTINE sample 
!
      USE control_flags, ONLY: job_type
      USE dimensions, ONLY: nsolid, ngas, nx, ny, nz
      USE domain_decomposition, ONLY: cell_owner, cell_g2l
      USE grid, ONLY: x, y, z
      USE io_files, ONLY: tempunit, logunit
      USE parallel, ONLY: mpime, root
!
      IMPLICIT NONE
      INTEGER :: ijk, imesh, i, j, k, n
      INTEGER :: nop, nprbs, nv
      INTEGER, ALLOCATABLE :: imesh_probe(:), indx(:)
      CHARACTER(LEN = 11) :: filnam
      CHARACTER(LEN = 4 ) :: lettera
      CHARACTER(LEN = 20) :: probenam
      !
      TYPE(probe_point), ALLOCATABLE :: probe(:)
!
! ... Allocate probes.
!
      ! ... If the axial points are sampled, add the number of column probes
      IF (isamp > 1) THEN
        nprbs = number_of_probes + nz
      ELSE 
        nprbs = number_of_probes
      END IF
!
      ALLOCATE(probe(nprbs))
      ALLOCATE(imesh_probe(nprbs))
      ALLOCATE(indx(nprbs))
!
! ... Assign axis probes
!
      DO nop = 1, nprbs - number_of_probes
        i = iiv
        k = nop
        probe(nop)%nop = nop
        probe(nop)%i = i
        probe(nop)%k = k
        probe(nop)%x = x(i)
        probe(nop)%z = z(k)
        IF (job_type == '2D') THEN
          imesh = i + (k-1)*nx
          imesh_probe(nop) = imesh
        ELSE IF (job_type == '3D') THEN
          j = jjv
          probe(nop)%j = j
          probe(nop)%y = y(j)
          imesh = i + (j-1)*nx + (k-1)*nx*ny
          imesh_probe(nop) = imesh
        END IF
      END DO
!
! ... If the number of "sparse" probes is not zero, read probe file
! ... and assign the remnants of the probe array
! 
      IF (number_of_probes /= 0) THEN
        IF (mpime == root) THEN
          OPEN(tempunit, FILE=probe_file, STATUS='OLD',ERR=199)
          READ(tempunit,*) 
          WRITE(logunit,*) 'Sampling probes: '
!
          WRITE(logunit,*) 'number_of_probes= ', number_of_probes 
          WRITE(logunit,*) 'nprbs= ', nprbs
!
        END IF
        DO nop = nprbs - number_of_probes + 1, nprbs
          probe(nop)%nop = nop
          IF (assign_index) THEN
                IF (job_type == '2D') THEN
                        IF (mpime == root) READ(tempunit,*) probe(nop)%i,probe(nop)%k
                        CALL bcast_integer(probe(nop)%i,1,root)
                        CALL bcast_integer(probe(nop)%k,1,root)
                        i = probe(nop)%i
                        k = probe(nop)%k
                        imesh = i + (k-1)*nx
                        imesh_probe(nop) = imesh
                        !
                        probe(nop)%x = x(i)
                        probe(nop)%z = z(k)
                ELSE IF (job_type == '3D') THEN
                        IF (mpime == root) READ(tempunit,*) probe(nop)%i, probe(nop)%j, probe(nop)%k
                        CALL bcast_integer(probe(nop)%i,1,root)
                        CALL bcast_integer(probe(nop)%j,1,root)
                        CALL bcast_integer(probe(nop)%k,1,root)
                        i = probe(nop)%i
                        j = probe(nop)%j
                        k = probe(nop)%k
                        imesh = i + (j-1)*nx + (k-1)*nx*ny
                        imesh_probe(nop) = imesh
                        !
                        probe(nop)%x = x(i)
                        probe(nop)%y = y(j)
                        probe(nop)%z = z(k)
                END IF
          ELSE
                IF (job_type == '2D') THEN
                        IF (mpime == root) READ(tempunit,*) probe(nop)%x, probe(nop)%z
                        CALL bcast_real(probe(nop)%x,1,root)
                        CALL bcast_real(probe(nop)%z,1,root)
                        DO i=1,nx
                          IF (x(i) <= probe(nop)%x) probe(nop)%i = i
                        END DO
                        DO k=1,nz
                          IF (z(k) <= probe(nop)%z) probe(nop)%k = k
                        END DO
                        i = probe(nop)%i
                        k = probe(nop)%k
                        imesh = i + (k-1)*nx
                        imesh_probe(nop) = imesh
                ELSE IF (job_type == '3D') THEN
                        IF (mpime == root) READ(tempunit,*) probe(nop)%x, probe(nop)%y, probe(nop)%z
                        CALL bcast_real(probe(nop)%x,1,root)
                        CALL bcast_real(probe(nop)%y,1,root)
                        CALL bcast_real(probe(nop)%z,1,root)
                        DO i=1,nx
                          IF (x(i) <= probe(nop)%x) probe(nop)%i = i
                        END DO
                        DO j=1,ny
                          IF (y(j) <= probe(nop)%y) probe(nop)%j = j
                        END DO
                        DO k=1,nz
                          IF (z(k) <= probe(nop)%z) probe(nop)%k = k
                        END DO
                        i = probe(nop)%i
                        j = probe(nop)%j
                        k = probe(nop)%k
                        imesh = i + (j-1)*nx + (k-1)*nx*ny
                        imesh_probe(nop) = imesh
                END IF
          END IF
        END DO
        IF (mpime == root) CLOSE(tempunit)
      END IF
!
! ... Sort the probes with progressively increasing index
! ... After the sorting, 'imesh_probe' contains the progressively
! ... increasing probe indexes, whereas 'indx' contains the
! ... probe indexes as read from the 'probe_file'
!
      CALL ikb07ad(imesh_probe(1), nprbs, indx(1))
!
      DO nop = 1, nprbs
        n = indx(nop)
        i = probe(n)%i
        j = probe(n)%j
        k = probe(n)%k
        IF (job_type == '3D') THEN
          probenam ='S'//lettera(n)//'_'//lettera(i)//'_'//&
                         lettera(j)//'_'//lettera(k)
        ELSE IF (job_type == '2D') THEN
          probenam ='S'//lettera(n)//'_'//lettera(i)//'_'//lettera(k)
        END IF
        OPEN(UNIT=tempunit, FILE=probenam, POSITION='APPEND')
        IF (cell_owner(imesh_probe(nop)) == mpime) THEN
          ijk = cell_g2l(imesh_probe(nop),mpime)
          IF (job_type == '3D') THEN
             WRITE(tempunit,100) time, p(ijk), &
             ug(ijk), vg(ijk), wg(ijk), tg(ijk), &
             (xgc(ijk,nv),nv=1, ngas), &
             (eps(ijk,nv), us(ijk,nv), & 
              vs(ijk,nv), ws(ijk,nv), &
              ts(ijk,nv), nv=1, nsolid), &
              rhom(ijk), mvm(ijk), pd(ijk)
          ELSE IF (job_type == '2D') THEN
             WRITE(tempunit,100) time, p(ijk), &
             ug(ijk), wg(ijk), tg(ijk), &
             (xgc(ijk,nv),nv=1, ngas), &
             (eps(ijk,nv), us(ijk,nv), &
              ws(ijk,nv), ts(ijk,nv),nv=1, nsolid), &
              rhom(ijk), mvm(ijk), pd(ijk)
          END IF
        END IF
        CLOSE(tempunit)
        CALL barrier
      END DO

 100  FORMAT( F8.2, 100(G14.6E3,1X) )
!
      DEALLOCATE(probe)
      DEALLOCATE(imesh_probe)
      DEALLOCATE(indx)
!
      RETURN
!
 199  CALL error ('sample','error in reading tempunit',tempunit)
!
      END SUBROUTINE sample
!---------------------------------------------------------------------
      SUBROUTINE column_axis_average(itn)
      USE control_flags, ONLY: job_type
      USE dimensions, ONLY: nsolid, ngas, nx, ny, nz
      USE domain_decomposition, ONLY: cell_owner, cell_g2l
      USE grid, ONLY: x, y, z
      USE io_files, ONLY: tempunit, logunit
      USE parallel, ONLY: mpime, root
      USE postp_output, ONLY: first_out, last_out, incr_out
!
      IMPLICIT NONE
      INTEGER :: ijk, imesh, i, j, k
      INTEGER :: nv
      INTEGER, ALLOCATABLE :: imesh_column_probe(:)
      INTEGER, INTENT(IN) :: itn
!
      ALLOCATE(imesh_column_probe(nz))
      imesh_column_probe = 0
!
! ... Set global array indexes
!
      DO k = 1, nz
        IF (job_type == '2D') THEN
          i = iiv 
          imesh = i + (k-1)*nx
          imesh_column_probe(k) = imesh
        ELSE IF (job_type == '3D') THEN
          i = iiv
          j = jjv
          imesh = i + (j-1)*nx + (k-1)*nx*ny
          imesh_column_probe(k) = imesh
        END IF
      END DO 
!
! ... Calculating the average arrays
!
      DO k = 1, nz
        IF (cell_owner(imesh_column_probe(k)) == mpime) THEN
          ijk = cell_g2l(imesh_column_probe(k),mpime)
          p_axis(k,itn) =  p(ijk)   
          p_av(k)  = p_av(k)  + p(ijk)   
          !
          ug_axis(k,itn) = ug(ijk)
          ug_av(k) = ug_av(k) + ug(ijk)
          !
          wg_axis(k,itn) = wg(ijk)
          wg_av(k) = wg_av(k) + wg(ijk)
          !
          tg_axis(k,itn) = tg(ijk)
          tg_av(k) = tg_av(k) + tg(ijk)
          !
          DO nv = 1, ngas  
            xgc_axis(k,nv,itn) = xgc(ijk,nv)
            xgc_av(k,nv) = xgc_av(k,nv) + xgc(ijk,nv)
          END DO 
          DO nv = 1, nsolid
            eps_axis(k,nv,itn) = eps(ijk,nv)
            eps_av(k,nv) = eps_av(k,nv) + eps(ijk,nv)
            !
            us_axis(k,nv,itn) = us(ijk,nv)
            us_av(k,nv)  = us_av(k,nv)  + us(ijk,nv)
            !
            ws_axis(k,nv,itn) = ws(ijk,nv)
            ws_av(k,nv)  = ws_av(k,nv)  + ws(ijk,nv)
            !
            ts_axis(k,nv,itn) = ts(ijk,nv)
            ts_av(k,nv)  = ts_av(k,nv)  + ts(ijk,nv)
          END DO
          rhom_axis(k,itn) = rhom(ijk)
          rhom_av(k) = rhom_av(k) + rhom(ijk)
          !
          wm_axis(k,itn) = wm(ijk)
          wm_av(k)  = wm_av(k)  + wm(ijk)
          !
          pd_axis(k,itn) = pd(ijk) 
          pd_av(k)  = pd_av(k)    + pd(ijk) 
          IF (job_type == '3D') THEN
            vg_axis(k,itn) = vg(ijk)
            vg_av(k) = vg_av(k)+vg(ijk)
            DO nv = 1, nsolid
              vs_axis(k,nv,itn) = vs(ijk,nv)
              vs_av(k,nv) = vs_av(k,nv)+vs(ijk,nv)
            END DO
          END IF
        END IF 
      END DO
!
      DEALLOCATE(imesh_column_probe)
!
      RETURN
      END SUBROUTINE column_axis_average
!----------------------------------------------------------------------
      SUBROUTINE column_axis_std(first_out, last_out, incr_out, cnt, md, md1)
      USE control_flags, ONLY: job_type
      USE dimensions, ONLY: nsolid, ngas, nz
!
      IMPLICIT NONE
      INTEGER :: k, itn
      INTEGER :: nv
      INTEGER, INTENT(IN) :: first_out, last_out, incr_out, cnt
      REAL*8, INTENT(IN) :: md, md1
!
! ... The average arrays are divided by the Number of outputs
!
      p_av    = p_av*md
      ug_av   = ug_av*md
      vg_av   = vg_av*md
      wg_av   = wg_av*md
      tg_av   = tg_av*md
      xgc_av  = xgc_av*md
      eps_av  = eps_av*md
      us_av   = us_av*md
      vs_av   = vs_av*md
      ws_av   = ws_av*md
      ts_av   = ts_av*md
      rhom_av = rhom_av*md
      wm_av  = wm_av*md
      pd_av   = pd_av*md

! ... Root processor gathers all data
!
      CALL parallel_sum_real(p_av,nz)
      CALL parallel_sum_real(ug_av,nz)
      CALL parallel_sum_real(vg_av,nz)
      CALL parallel_sum_real(wg_av,nz)
      CALL parallel_sum_real(tg_av,nz)
      CALL parallel_sum_real(xgc_av,nz*ngas)
      CALL parallel_sum_real(eps_av,nz*nsolid)
      CALL parallel_sum_real(us_av,nz*nsolid)
      CALL parallel_sum_real(vs_av,nz*nsolid)
      CALL parallel_sum_real(ws_av,nz*nsolid)
      CALL parallel_sum_real(ts_av,nz*nsolid)
      CALL parallel_sum_real(rhom_av,nz)
      CALL parallel_sum_real(wm_av,nz)
      CALL parallel_sum_real(pd_av,nz)
!
! ... Root processor gathers all data
!
      CALL parallel_sum_real(p_axis,nz*cnt)
      CALL parallel_sum_real(ug_axis,nz*cnt)
      CALL parallel_sum_real(vg_axis,nz*cnt)
      CALL parallel_sum_real(wg_axis,nz*cnt)
      CALL parallel_sum_real(tg_axis,nz*cnt)
      CALL parallel_sum_real(xgc_axis,nz*ngas*cnt)
      CALL parallel_sum_real(eps_axis,nz*nsolid*cnt)
      CALL parallel_sum_real(us_axis,nz*nsolid*cnt)
      CALL parallel_sum_real(vs_axis,nz*nsolid*cnt)
      CALL parallel_sum_real(ws_axis,nz*nsolid*cnt)
      CALL parallel_sum_real(ts_axis,nz*nsolid*cnt)
      CALL parallel_sum_real(rhom_axis,nz*cnt)
      CALL parallel_sum_real(wm_axis,nz*cnt)
      CALL parallel_sum_real(pd_axis,nz*cnt)
!
      DO itn = 1, cnt
      !
      ! ... Calculating the standard deviation arrays.
      !
        DO k = 1, nz
          p_sd(k)  = p_sd(k)  + (p_axis(k,itn)-p_av(k))**2
          ug_sd(k) = ug_sd(k) + (ug_axis(k,itn)-ug_av(k))**2
          wg_sd(k) = wg_sd(k) + (wg_axis(k,itn)-wg_av(k))**2
          tg_sd(k) = tg_sd(k) + (tg_axis(k,itn)-tg_av(k))**2
          DO nv = 1, ngas  
            xgc_sd(k,nv) = xgc_sd(k,nv) + (xgc_axis(k,nv,itn)-xgc_av(k,nv))**2
          END DO 
          DO nv = 1, nsolid
            eps_sd(k,nv) = eps_sd(k,nv) + (eps_axis(k,nv,itn)-eps_av(k,nv))**2
            us_sd(k,nv)  = us_sd(k,nv)  + (us_axis(k,nv,itn)-us_av(k,nv))**2
            ws_sd(k,nv)  = ws_sd(k,nv)  + (ws_axis(k,nv,itn)-ws_av(k,nv))**2
            ts_sd(k,nv)  = ts_sd(k,nv)  + (ts_axis(k,nv,itn)-ts_av(k,nv))**2
          END DO
          rhom_sd(k) = rhom_sd(k) + (rhom_axis(k,itn)-rhom_av(k))**2
          wm_sd(k)  = wm_sd(k)  + (wm_axis(k,itn)-wm_sd(k))**2
          pd_sd(k)  = pd_sd(k)    + (pd_axis(k,itn)-pd_av(k))**2
          IF (job_type == '3D') THEN
            vg_sd(k) = vg_sd(k)+(vg_axis(k,itn)-vg_av(k))**2
            DO nv = 1, nsolid
              vs_sd(k,nv) = vs_sd(k,nv)+(vs_axis(k,nv,itn)-vs_av(k,nv))**2
            END DO
          END IF
        END DO
      END DO
!
! ... The standard deviation arrays are divided by the Number of outputs
!
      p_sd    = SQRT(p_sd*md1)
      ug_sd   = SQRT(ug_sd*md1)
      vg_sd   = SQRT(vg_sd*md1)
      wg_sd   = SQRT(wg_sd*md1)
      tg_sd   = SQRT(tg_sd*md1)
      xgc_sd  = SQRT(xgc_sd*md1)
      eps_sd  = SQRT(eps_sd*md1)
      us_sd   = SQRT(us_sd*md1)
      vs_sd   = SQRT(vs_sd*md1)
      ws_sd   = SQRT(ws_sd*md1)
      ts_sd   = SQRT(ts_sd*md1)
      rhom_sd = SQRT(rhom_sd*md1)
      wm_sd   = SQRT(wm_sd*md1)
      pd_sd   = SQRT(pd_sd*md1)
!
      RETURN
      END SUBROUTINE column_axis_std
!----------------------------------------------------------------------
      SUBROUTINE compute_plume_profile
      USE control_flags, ONLY: job_type
      USE dimensions
      USE grid, ONLY: dx, dy, dz, rb, itc
      USE domain_mapping, ONLY: ncint, meshinds
      IMPLICIT NONE
!
      REAL*8 :: ds, invrhom, invsurf, pi
      INTEGER :: ijk, imesh, i, j, k
!
      pi = 4.D0*ATAN(1.D0)
!
      surface = 0.D0
      rhom_z  = 0.D0
      um_z    = 0.D0
      vm_z    = 0.D0
      wm_z    = 0.D0
      tm_z    = 0.D0
      invrhom = 0.D0
      invsurf = 0.D0
      ds = 0.D0
!
      DO ijk = 1, ncint
        IF ( epst(ijk) >= 1.D-8) THEN
          CALL meshinds(ijk,imesh,i,j,k)
          IF (job_type == '2D' .AND. itc == 1) THEN
            ds = pi*dx(i)*(2.D0*rb(i)-dx(i))
            rhom_z(k) = rhom_z(k) + rhom(ijk) * ds
            um_z(k) = um_z(k) + rhom(ijk) * um(ijk) * ds
            wm_z(k) = wm_z(k) + rhom(ijk) * wm(ijk) * ds
            ! ... Averaged temperature assumes constant mixture Cp
            tm_z(k) = tm_z(k) + tm(ijk) * ds
            surface(k) = surface(k) + ds
          ELSE IF (job_type == '3D') THEN
            ds = dx(i)*dy(j)
            rhom_z(k) = rhom_z(k) + rhom(ijk) * ds
            um_z(k) = um_z(k) + rhom(ijk) * um(ijk) * ds
            vm_z(k) = vm_z(k) + rhom(ijk) * vm(ijk) * ds
            wm_z(k) = wm_z(k) + rhom(ijk) * wm(ijk) * ds
            ! ... Averaged temperature assumes constant mixture Cp
            tm_z(k) = tm_z(k) + tm(ijk) * ds
            surface(k) = surface(k) + ds
          END IF
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
        IF (job_type == '3D') vm_z(k) = vm_z(k) * invrhom *  invsurf
        wm_z(k) = wm_z(k) * invrhom *  invsurf
        tm_z(k) = tm_z(k) * invsurf
      END DO
!
      RETURN
      END SUBROUTINE compute_plume_profile
!-----------------------------------------------------------------------
      END MODULE sample_points
!----------------------------------------------------------------------
