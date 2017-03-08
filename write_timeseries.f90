subroutine write_timeseries()
use sim_param,only:path,u,v,w,theta
use param
use types,only:rprec
implicit none

real(kind=rprec),dimension(nz-1)::tss_u, tss_v, tss_w, tss_T
real(kind=rprec),dimension(nz_tot-1)::tss_utot, tss_vtot, tss_wtot, tss_Ttot
integer::jz

$if ($MPI)
  $define $lbz 0
  integer::recvcounts(nproc)
  integer::displs(nproc)
$else    
  $define $lbz 1
$endif

do jz=1,nz-1
     tss_u(jz) = u(64,64,jz)
     tss_v(jz) = v(64,64,jz)
     tss_w(jz) = w(64,64,jz)
     tss_T(jz) = theta(64,64,jz)
end do

$if ($MPI)
   recvcounts = size(tss_u)
   displs = coord_of_rank*recvcounts
   call mpi_gatherv (tss_u(1),size(tss_u),MPI_RPREC,   &
                     tss_utot(1),recvcounts,displs,    &
                     MPI_RPREC,rank_of_coord(0),comm,ierr)

   recvcounts = size(tss_v)
   displs = coord_of_rank*recvcounts
   call mpi_gatherv (tss_v(1),size(tss_v),MPI_RPREC,   &
                     tss_vtot(1),recvcounts,displs,    &
                     MPI_RPREC,rank_of_coord(0),comm,ierr)

   recvcounts = size(tss_w)
   displs = coord_of_rank*recvcounts
   call mpi_gatherv (tss_w(1),size(tss_w),MPI_RPREC,   &
                     tss_wtot(1),recvcounts,displs,    &
                     MPI_RPREC,rank_of_coord(0),comm,ierr)

   recvcounts = size(tss_T)
   displs = coord_of_rank*recvcounts
   call mpi_gatherv (tss_T(1),size(tss_T),MPI_RPREC,   &
                     tss_Ttot(1),recvcounts,displs,    &
                     MPI_RPREC,rank_of_coord(0),comm,ierr)
$else
   tss_utot = tss_u
   tss_vtot = tss_v
   tss_wtot = tss_w
   tss_Ttot = tss_T
$endif

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
   write(501,8003) tss_utot(3), tss_vtot(3), tss_wtot(3), tss_Ttot(3),&
                   tss_utot(32), tss_vtot(32), tss_wtot(32), tss_Ttot(32), &
                   tss_utot(96), tss_vtot(96), tss_wtot(96), tss_Ttot(96)
end if
8003 format (12(E11.4,2x))

end subroutine write_timeseries

