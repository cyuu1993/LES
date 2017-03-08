subroutine divstress_uv (divt, tx, ty, tz)
!--provides divt, jz=1:nz-1
use types,only:rprec
use param,only:ld,ny,nz, BOGUS, VERBOSE
implicit none
$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif
real(kind=rprec),dimension(ld,ny,$lbz:nz),intent(out)::divt
real (rprec), dimension (ld, ny, $lbz:nz), intent (in) :: tx, ty, tz
! sc: we should be able to save some memory here!
! do this by writing a divergence subroutine--then do not store derivs 
real(kind=rprec),dimension(ld,ny,$lbz:nz)::dtxdx,dtydy, dtzdz

logical, parameter :: DEBUG = .true.

if (VERBOSE) write (*, *) 'started divstress_uv'
 
! compute stress gradients      
!--MPI: tx 1:nz-1 => dtxdx 1:nz-1
call ddx(dtxdx, tx)  !--really should replace with ddxy (save an fft)
!$if ($MPI)
!  dtdx(:, :, 0) = BOGUS
!$endif
!dtxdx(:, :, nz) = BOGUS

!--MPI: ty 1:nz-1 => dtdy 1:nz-1
call ddy(dtydy, ty)
!$if ($MPI)
!  dtdy(:, :, 0) = BOGUS
!$endif
!dtydy(:, :, nz) = BOGUS

!--MPI: tz 1:nz => ddz_w limits dtzdz to 1:nz-1, except top process 1:nz
call ddz_w(dtzdz, tz)
!$if ($MPI)
!  dtzdz(:, :, 0) = BOGUS
!$endif
!if (USE_MPI .and. coord < nproc-1) dtzdz(:, :, nz) = BOGUS

!--MPI following comment only true at bottom process
! the following gives bad results...but it seems like i the
! issue should be taken care of somewhere
! need to correct wall level, since tz will be on uv-node there
!      dtzdz(:,:,1) = (tz(:,:,2)-tz(:,:,1))/(0.5*dz)

!--only 1:nz-1 are valid
divt(:, :, 1:nz-1) = dtxdx(:, :, 1:nz-1) + dtydy(:, :, 1:nz-1) +  &
                     dtzdz(:, :, 1:nz-1)

!--Set ld-1, ld to 0 (or could do BOGUS)
divt(ld-1:ld, :, 1:nz-1) = 0._rprec

$if ($MPI)
  divt(:, :, 0) = BOGUS
$endif
divt(:, :, nz) = BOGUS

if (VERBOSE) write (*, *) 'finished divstress_uv'

end subroutine divstress_uv
