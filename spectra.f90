subroutine spectra()
! Calls function to calculate spectra and the writes it

use sim_param,only:path,u,v,w,theta,avg_spectra_uvwT,avg_cospectra_uvwT
use param
use fft
use types,only:rprec
implicit none

real(kind=rprec),dimension(4,nx/2,nz-1)::spectra_uvwT
real(kind=rprec),dimension(4,nx/2,nz_tot-1)::spectra_uvwT_tot
real(kind=rprec),dimension(3,nx/2,nz-1)::cospectra_uvwT
real(kind=rprec),dimension(3,nx/2,nz_tot-1)::cospectra_uvwT_tot
real(kind=rprec),dimension(ld,ny)::temp_val
integer::k,jx,jy,jz,z

$if ($MPI)
  $define $lbz 0
  integer::recvcounts(nproc)
  integer::displs(nproc)
$else
  $define $lbz 1
$endif

do jz=1,nz-1
   call calc_spectra(u(:,:,jz),spectra_uvwT(1,:,jz))
   call calc_spectra(v(:,:,jz),spectra_uvwT(2,:,jz))
   call calc_spectra(w(:,:,jz),spectra_uvwT(3,:,jz))
   call calc_spectra(theta(:,:,jz),spectra_uvwT(4,:,jz))
   temp_val = u(:,:,jz)*w(:,:,jz)
   call calc_spectra(temp_val,cospectra_uvwT(1,:,jz))
   temp_val = v(:,:,jz)*w(:,:,jz)
   call calc_spectra(temp_val,cospectra_uvwT(2,:,jz))
   temp_val = theta(:,:,jz)*w(:,:,jz)
   call calc_spectra(temp_val,cospectra_uvwT(3,:,jz))
enddo


do jy=1,4
  do jx=1,nx/2
    do jz=1,nz-1
       avg_spectra_uvwT(jy,jx,jz) = avg_spectra_uvwT(jy,jx,jz) &
                                  + spectra_uvwT(jy,jx,jz)
    end do
 end do
end do
do jy=1,3
  do jx=1,nx/2
    do jz=1,nz-1
       avg_cospectra_uvwT(jy,jx,jz) = avg_cospectra_uvwT(jy,jx,jz) &
                                  + cospectra_uvwT(jy,jx,jz)
    end do
 end do
end do
end subroutine spectra
