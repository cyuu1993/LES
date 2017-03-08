subroutine rmsdiv(rms)
! actually, this is NOT the rms divergence of velocity, its like an
! l_1 norm or something.
use types,only:rprec
use param
use sim_param, only : du=>dudx, dv=>dvdy, dw=>dwdz
use io, only : jt_total
use debug_mod
implicit none
integer::jx,jy,jz
integer :: jz_max
real(kind=rprec)::rms

logical, parameter :: DEBUG = .false.
logical, parameter :: norm_magdu = .false.

logical, parameter :: write_out = .false.
integer, parameter :: nwrite_out = 2000  !--must be multiple of base (io) or
                                         !  will not get anything

real (rprec), parameter :: thresh = 0._rprec

character (64) :: file_out

integer :: n

real (rprec) :: magdu, div
$if ($MPI)
  real (rprec) :: rms_global
$endif

!---------------------------------------------------------------------

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
  jz_max = nz-1
else
  !jz_max = nz
  jz_max = nz-1
end if

rms = 0._rprec

if (norm_magdu) then

  if (write_out) then
  
    if (modulo (jt, nwrite_out) == 0) then
  
      write (file_out, '(a,i6.6,a)') 'output/div', jt_total+1, '.dat'
                                         !--+1 jt_total since jt_total not
                                         !  updated yet
      open (1, file=file_out)

      !--this assumes you are using tecplot
      write (1, *) 'variables = "jx" "jy" "jz" "div" "magdu" "div/magdu"'
      write (1, *) 'zone, f=point, i=', nx, 'j=', ny, ',k=', jz_max
      
    end if

  end if

  n = 0
  do jz=1,jz_max
  do jy=1,ny
  do jx=1,nx

     magdu = sqrt (du(jx, jy, jz)**2 + dv(jx, jy, jz)**2 + dw(jx, jy, jz)**2)
     div = du(jx, jy, jz) + dv(jx, jy, jz) + dw(jx, jy, jz) 

     if (magdu > thresh) then
       n = n + 1
       rms = rms + abs (div) / magdu
     end if

     if (write_out) then
       if (modulo (jt, nwrite_out) == 0) then
         write (1, '(3(i0,1x),3(es12.5,1x))') jx, jy, jz, div, magdu,  &
                                   div / max (magdu, epsilon (0._rprec))
       end if
     end if
     
  end do
  end do
  end do
  rms = rms / n
 
  if (write_out) then
    if (modulo (jt, nwrite_out) == 0) close (1)
  end if

else

  do jz=1,jz_max
  do jy=1,ny
  do jx=1,nx
     rms=rms+abs(du(jx,jy,jz)+dv(jx,jy,jz)+dw(jx,jy,jz))
  end do
  end do
  end do
  rms=rms/(nx*ny*(jz_max))

end if

$if ($MPI)
  call mpi_reduce (rms, rms_global, 1, MPI_RPREC, MPI_SUM, 0, comm, ierr)
  if (rank == 0) then
    rms = rms_global/nproc
    !write (*, *) 'rms_global = ', rms_global/nproc
  end if
  !if (rank == 0) rms = rms_global/nproc  !--its rank here, not coord
$endif

end subroutine rmsdiv
