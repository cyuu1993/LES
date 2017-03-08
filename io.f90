module io
use messages
use types,only:rprec
use param, only : ld, nx, ny, nz, write_inflow_file, path,  &
                  use_avgslice, USE_MPI, coord, rank, nproc,      &
                  average_dim_num, nz_tot,jt_total,p_count,dt,z_i,u_star
implicit none
private
public :: jt_total
public :: openfiles,output_loop,output_final,            &
          mean_u,mean_u2,mean_v,mean_v2,mean_w,mean_w2,  &
          inflow_read, inflow_write, avg_stats
public :: average_dim_select_flag, dim1_size, dim2_size,        &
          dim1_global, dim2_global, collocate_MPI_averages     
          !compute_avg_var 

integer,parameter::num_hour_out=1
integer,parameter::base=50000,nwrite=base
! SKS
! Shouldnt be hard coded..base corresponds to the 
! time interval after which you want to write file
! So can be a factor of nsteps.
! SKS

logical,parameter:: io_spec=.false.,output_fields_3d_flag=.false.
integer,parameter::spec_write_freqz=600, fields_3d_write_freqz=p_count*6
integer,parameter::spec_write_start=1,spec_write_end=24*base
!! --------------------------------------------------------------------
!! The following block defines parameters for instantaneous slice output
!! slice_inst sets the value of the y node for the chosen x-z inst slice
!! inst_slice_freqz controls the output frequency
!! The 5 slices outputted every inst_slice_freqz (i.e. u,v,w,T,Cs in this order) ...
!! ... are saved in the 3rd dimension of inst_slice_array for inst_array_lim times ...
!! ... following which this array is outputted as a binary file and the process 
!! ... starts over
!! Therefore, each file will contain inst_array_lim x-z slices of 5 variables
!! This information is important for post-processing of these binary files
!! --------------------------------------------------------------------

logical,parameter:: inst_slice_flag=.false.
integer,parameter:: num_vars=4 ! works currently only for u,v,w,T due to the size difference in Cs
integer,parameter:: slice_inst=(nz_tot-1)/2, inst_slice_freqz=5, inst_array_lim=200

logical, parameter :: cumulative_time = .true.
character (*), parameter :: fcumulative_time = path // 'total_time.dat'

integer, parameter :: n_avg_stats = 5000000 !--interval for updates in avg_stats
character (*), parameter :: end_hdr_avg = '# end header'

!! --------------------------------------------------------------------
!! The following block defines parameters for use in avgslice and scalar_slice
!! --------------------------------------------------------------------
integer,parameter :: average_dim_select_flag=1-(average_dim_num/2) 
! The variable average_dim_select_flag generates the following values based
! on the value of average_dim_num in param.f90 :- 
! a) average_dim_num = 2 : 
! 	average_dim_select_flag = 0 ==> average the 3d array over x and y and output the z profile
! b) average_dim_num = 1 : 
! 	average_dim_select_flag = 1 ==> average the 3d array over y and output the (x,z) profile
integer, parameter :: dim1_size=average_dim_select_flag*(nx-nz+1)+nz-1
integer, parameter :: dim2_size=average_dim_select_flag*(nz-2)+1
integer, parameter :: dim1_global=average_dim_select_flag*(nx-nz_tot+1)+nz_tot-1
integer, parameter :: dim2_global=average_dim_select_flag*(nz_tot-2)+1
!! --------------------------------------------------------------------
!! --------------------------------------------------------------------

!!!!  time_spec>0 output time series spectrum (need additional calcu.)
integer,parameter::time_spec=0
integer::n_obs, jt_total_init
integer,allocatable::obs_pt(:,:)

!!!!  io_mean=.true. output small domain time-averaged velocity
logical,parameter::io_mean=.true.
integer,parameter::jx_pls=1,jx_ple=nx,width=0
integer,parameter::jy_pls=ny/2-width,jy_ple=ny/2+width
real(kind=rprec),dimension(jx_pls:jx_ple,jy_pls:jy_ple,nz)::&
     mean_u,mean_v,mean_w,mean_u2,mean_v2,mean_w2

!!!!  io_lambda2
logical,parameter::io_lambda2=.false.
real(kind=rprec),dimension(nx,ny,nz)::lam2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! file number 1-10 are used for temporary use
! 11-19 are basic output
! 20-40 are avgslice
! use >50 for debugging
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine openfiles()
use sim_param,only:path
implicit none

!--to hold file names
character (64) :: temp
character (64) :: fCS1plan, fCS2plan, fCS4plan, fVISCplan,  &
                  fDISSplan, fCS1Vplan, fCS2Vplan, fCS4Vplan

integer::i

logical :: exst

if (cumulative_time) then

  inquire (file=fcumulative_time, exist=exst)
  if (exst) then
    open (1, file=fcumulative_time)
    read (1, *) jt_total
    jt_total_init=jt_total 
    close (1)
  else  !--assume this is the first run on cumulative time
    write (*, *) 'file ', fcumulative_time, ' not found'
    write (*, *) 'assuming jt_total = 0'
    jt_total = 0
    jt_total_init=jt_total 
  end if

end if

if ((.not. USE_MPI) .or. (USE_MPI .and. rank == 0)) then
  open(13,file=path//'output/check_ke.out',status="unknown",position="append")
  ! SKS
  open(14,file=path//'output/Cs2byPr_t.out',status="unknown",position="append")
  open(84,file=path//'output/spec_u.out',status="unknown",position="append")
  open(85,file=path//'output/spec_v.out',status="unknown",position="append")
  open(86,file=path//'output/spec_w.out',status="unknown",position="append")
  open(87,file=path//'output/spec_T.out',status="unknown",position="append")
  open(88,file=path//'output/spec_uw.out',status="unknown",position="append")
  open(89,file=path//'output/spec_vw.out',status="unknown",position="append")
  open(90,file=path//'output/spec_Tw.out',status="unknown",position="append")
  open(91,file=path//'output/spec_freq.out',status="unknown",position="append")
  open(111,file=path//'output/tkeVsz.dat',status="unknown",position="append")
  ! SKS
end if

if(time_spec.gt.0)then
  open(15,file=path//'output/velspec.out',form='unformatted',position='append')
  if(jt_total.eq.0)rewind(15)
endif

if(io_mean)then
  open(51,file=path//'output/mean_u.out',form='unformatted',position='append')
  if(jt_total.eq.0)then
    rewind(51)
    write(51)jx_pls,jx_ple,jy_pls,jy_ple
  endif
endif

fCS1plan = path // 'output/CS1plan.out'
fCS2plan = path // 'output/CS2plan.out'
fCS4plan = path // 'output/CS4plan.out'
fVISCplan = path // 'output/VISCplan.out'
fDISSplan = path // 'output/DISSplan.out'
fCS1Vplan = path // 'output/CS1Vplan.out'
fCS2Vplan = path // 'output/CS2Vplan.out'
fCS4Vplan = path // 'output/CS4Vplan.out'

$if ($MPI)
  !--append coordinate identifiers
  write (temp, '(".c",i0)') coord
  fCS1plan = trim (fCS1plan) // temp
  fCS2plan = trim (fCS2plan) // temp
  fCS4plan = trim (fCS4plan) // temp
  fVISCplan = trim (fVISCplan) // temp
  fDISSplan = trim (fDISSplan) // temp
  fCS1Vplan = trim (fCS1Vplan) // temp
  fCS2Vplan = trim (fCS2Vplan) // temp
  fCS4Vplan = trim (fCS4Vplan) // temp
$endif

if(time_spec.gt.0)then
open(1,file=path//'obs.pt')
read(1,*)n_obs
allocate(obs_pt(1:2,n_obs))
do i=1,n_obs
read(1,*)obs_pt(1:2,i)
enddo
close(1)
endif

end subroutine openfiles
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_loop
use param,only:output,dt,c_count,S_FLAG,SCAL_init,jt,jan_diurnal_run
use sim_param,only:path,u,v,w,dudz,dudx,p,&
     RHSx,RHSy,RHSz,theta, txx, txy, txz, tyy, tyz, tzz
use sgsmodule,only:Cs_opt2
use scalars_module,only:sgs_t3
use scalars_module2,only:scalar_slice,budget_TKE_scalar
implicit none
real(kind=rprec),dimension(nz)::u_ndim
character(len=20)::req

! SKS
character (100) :: fname, temp  ! With 64 the code was giving an error !
! character (64) :: fname, temp ! sort of segmentation fault i guess.
! SKS

integer::jx,jy,jz
integer:: fields_3d_write_freqz_temp

$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif

real(kind=rprec),dimension(ld,$lbz:nz,num_vars*inst_array_lim),save::inst_slice_array
integer,save:: inst_slice_counter

!jt_total=jt_total+1  !--moved into main
call calculate_mean

if ((use_avgslice) .and. (mod(jt,c_count)==0)) then
       if ((S_FLAG) .and. (jt.GE.SCAL_init)) then ! Output scalar variables
         !call avgslice()
         call scalar_slice() ! Uses file unit numbers (36-47)
         call MM_budget_slice()
         ! SKS
         call budget_TKE_scalar()
         ! SKS
       elseif (.not. S_FLAG) then
         call avgslice()
       end if
end if

if (output) then
  if ((mod(jt_total,base)==0).and.(jt_total.ge.1)) then
    if (S_FLAG) then
       write (fname, '(a,i6.6,a)') path // 'output/vel_sc', jt_total, '.out'
    else
       write (fname, '(a,i6.6,a)') path // 'output/vel', jt_total, '.out'
    end if
    $if ($MPI)
       write (temp, '(".c",i0)') coord
       fname = trim (fname) // temp
    $endif

    open(1,file=fname,form='unformatted')
    call checkpoint (1)
    close(1)
    if (io_mean) call io_mean_out
  end if
end if

 if (S_FLAG) then ! If Block 1
      if ((inst_slice_flag) .AND. mod(jt_total, inst_slice_freqz)==0) then !If Block 2
        if (jt .eq. inst_slice_freqz) inst_slice_counter=1
        if (jt .eq. inst_slice_freqz) print *,'inst_slice_counter = ',inst_slice_counter
            
            if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
                print *,'inst_slice_counter = ',inst_slice_counter
            end if

       inst_slice_array(:,:,(inst_slice_counter-1)*num_vars+1) = u(:,slice_inst,:)
       inst_slice_array(:,:,(inst_slice_counter-1)*num_vars+2) = v(:,slice_inst,:)
       inst_slice_array(:,:,(inst_slice_counter-1)*num_vars+3) = w(:,slice_inst,:)
       inst_slice_array(:,:,(inst_slice_counter-1)*num_vars+4) = theta(:,slice_inst,:)

         if (mod(inst_slice_counter,inst_array_lim) .eq. 0) then !If Block 3 begins
            if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then !If Block 4
                print *,'INSIDE:inst_slice_counter = ',inst_slice_counter
            end if !If Block 4 ends
          write(fname,'(A,i6.6,A)')path//'output/fields_3d/inst_slices_uvwT_till_',jt_total,'.bin'
          $if ($MPI)
            write (temp, '(".c",i0)') coord
            fname = trim (fname) // temp
          $endif
           open(1,file=fname,form='unformatted')
           write(1) real(inst_slice_array)
           close(1)
           inst_slice_array=0._rprec; inst_slice_counter=0;
         end if ! If Block 3 ends

       inst_slice_counter = inst_slice_counter+1 !increment the slice counter by 1
            if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then !If Block 4
                print *,'inst_slice_counter = ',inst_slice_counter
            end if !If Block 4 ends
     end if ! If Block 2 ends
       
       fields_3d_write_freqz_temp=50

    if ((output_fields_3d_flag) .and. mod(jt_total,fields_3d_write_freqz_temp)==0) then !If Block 5 begins

    write(fname,'(A,i6.6,A)')path//'output/fields_3d/o_uvwT_',jt_total,'.bin'
    $if ($MPI)
      write (temp, '(".c",i0)') coord
      fname = trim (fname) // temp
    $endif
    open(1,file=fname,form='unformatted')
    write(1) real(u),real(v),real(w),real(theta); close(1)
    
    end if ! If Block 5 ends
 end if ! If Block 1 ends

  if ((io_spec) .and. (jt_total .gt. spec_write_start .and. jt_total .le. spec_write_end)) then 
   if (mod(jt_total,spec_write_freqz)==0) call post_spec(jt_total)
  end if

  if (time_spec.gt.0) call timeseries_spec

end subroutine output_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!subroutine compute_avg_var(avg_var,raw_var,output_dim)
!use param,only:nx,ny,nz,c_count,p_count
!integer :: output_dim
!real(kind=rprec),dimension(dim1_size,dim2_size)::avg_var
!real(kind=rprec),dimension(1:nx,1:ny,1:nz-1):: raw_var
!real(kind=rprec):: fr
!     
!fr=(1._rprec/real(p_count,kind=rprec))*real(c_count,kind=rprec)
!select case (output_dim)
!  case(1) !average over y and t
!      avg_var=avg_var+fr*sum(raw_var(1:nx,1:ny,1:nz-1),dim=2)/ny
!  case(2) ! average over x,y and t
!      avg_var(:,1)=avg_var(:,1)+fr*sum(sum(raw_var(1:nx,1:ny,1:nz-1),dim=1),dim=2)/(nx*ny)
!end select
!end subroutine compute_avg_var

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--The following subroutine does the collocation of the MPI arrays for
! averaging in avgslice and scalar_slice (in scalars_module2.f90)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine collocate_MPI_averages(avg_var_proc,avg_var_tot_domain,file_ind)
use param
$if ($MPI)
  integer :: recvcounts(nproc)
  integer :: displs(nproc)
$endif
integer :: ind1,ind2,file_ind
real(kind=rprec),dimension(dim1_size,dim2_size)::avg_var_proc
real(kind=rprec),dimension(dim1_global,dim2_global)::avg_var_tot_domain

  avg_var_tot_domain=0._rprec
$if ($MPI)
  recvcounts = size (avg_var_proc)
  displs = coord_of_rank * recvcounts 
  call mpi_gatherv (avg_var_proc(1,1), size (avg_var_proc), MPI_RPREC,                &
                    avg_var_tot_domain(1, 1), recvcounts, displs, MPI_RPREC,  &
                    rank_of_coord(0), comm, ierr)
$else
  avg_var_tot_domain=avg_var_proc
$endif

  if((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then  
      if (average_dim_num .eq. 1) then
         do ind2=1,nz_tot-1
          write(file_ind,5168) jt*dt,(avg_var_tot_domain(ind1,ind2),ind1=1,nx)
         end do
      else if (average_dim_num .eq. 2) then
         write(file_ind,5168) jt*dt,(avg_var_tot_domain(ind1,1),ind1=1,nz_tot-1)
      end if
         close(file_ind)
  end if

5168     format(1400(E14.5))
end subroutine collocate_MPI_averages

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--The following subroutine does the collocation of the MPI arrays for
! averaging in avgslice and scalar_slice (in scalars_module2.f90)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine collocate_MPI_averages_N(avg_var_proc,avg_var_tot_domain,file_ind,filename_str)
!subroutine collocate_MPI_averages(avg_var_proc,avg_var_tot_domain,file_ind)
use param
$if ($MPI)
  integer :: recvcounts(nproc)
  integer :: displs(nproc)
$endif
integer :: ind1,ind2,file_ind
character (*),intent(in) :: filename_str
character (len=256) :: local_filename
real(kind=rprec),dimension(dim1_size,dim2_size)::avg_var_proc
real(kind=rprec),dimension(dim1_global,dim2_global)::avg_var_tot_domain

local_filename=path//'output/aver_'//trim(filename_str)//'.out'

avg_var_tot_domain=0._rprec
$if ($MPI)
  recvcounts = size (avg_var_proc)
  displs = coord_of_rank * recvcounts 
  call mpi_gatherv (avg_var_proc(1,1), size (avg_var_proc), MPI_RPREC,       &
                    avg_var_tot_domain(1, 1), recvcounts, displs, MPI_RPREC, &
                    rank_of_coord(0), comm, ierr)
$else
  avg_var_tot_domain=avg_var_proc
$endif

  if((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then  
        open(file_ind,file=trim(local_filename),status="unknown",position="append")
        if (average_dim_num .eq. 1) then
           do ind2=1,nz_tot-1
            write(file_ind,5168) jt*dt,(avg_var_tot_domain(ind1,ind2),ind1=1,nx)
           end do
        else if (average_dim_num .eq. 2) then
           write(file_ind,5168) jt*dt,(avg_var_tot_domain(ind1,1),ind1=1,nz_tot-1)
        end if
        close(file_ind)
  end if
5168     format(1400(E14.5))
end subroutine collocate_MPI_averages_N

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
subroutine avgslice
use sim_param,only:path,u,v,w,dudz,dvdz,txx,txz,tyy,tyz,tzz,p
use param,only:dz,p_count,c_count,jt
use sgsmodule,only:Cs_opt2,Cs_Ssim,Beta_avg,Betaclip_avg
implicit none
integer::i,j,k
real(kind=rprec),dimension(nx,nz-1),save::ap,au,av,aw,p2,u2,v2,w2,auw,avw,acs
real(kind=rprec),dimension(nx,nz-1),save::adudz,advdz,aCs_Ssim,abeta_sgs,abetaclip_sgs
real(kind=rprec),dimension(nx,nz-1),save::atxx,atxz,atyy,atyz,atzz
real(kind=rprec),dimension(nx,nz-1),save::u3,v3,w3
real(kind=rprec)::tu1,tv1,tw1,ttxx,ttxz,ttyy,ttyz,ttzz,tdudz,tdvdz,&
     tu2,tv2,tw2,tp1,tp2,tuw,tvw,tCs,fr,arg1,arg2,tu3,tv3,tw3     
real(kind=rprec)::tCs_Ssim
real(kind=rprec),dimension(:,:),allocatable::avg_out

fr=(1._rprec/real(p_count,kind=rprec))*real(c_count,kind=rprec)
do k=1,Nz-1
do i=1,Nx
   tu1=0._rprec;tv1=0._rprec;tw1=0._rprec;tp1=0._rprec
   ttxx=0._rprec;ttxz=0._rprec;ttyy=0._rprec;ttyz=0._rprec
   ttzz=0._rprec;tdudz=0._rprec;tdvdz=0._rprec;tu2=0._rprec
   tv2=0._rprec;tw2=0._rprec;tp2=0._rprec;tuw=0._rprec;tvw=0._rprec
   tCs=0._rprec;tCs_Ssim=0._rprec;tu3=0._rprec;tv3=0._rprec;tw3=0._rprec

   do j=1,Ny
      tu1=tu1+u(i,j,k)
      tv1=tv1+v(i,j,k)
      tw1=tw1+w(i,j,k)
      tp1=tp1+p(i,j,k)
      ttxx=ttxx+txx(i,j,k)
      ttxz=ttxz+txz(i,j,k)
      ttyy=ttyy+tyy(i,j,k)
      ttyz=ttyz+tyz(i,j,k)
      ttzz=ttzz+tzz(i,j,k)
      tdudz=tdudz+dudz(i,j,k)
      tdvdz=tdvdz+dvdz(i,j,k)
      tu2=tu2+u(i,j,k)*u(i,j,k)
      tv2=tv2+v(i,j,k)*v(i,j,k)
      tw2=tw2+w(i,j,k)*w(i,j,k)
      tp2=tp2+p(i,j,k)*p(i,j,k)
      tCs=tCs+sqrt(Cs_opt2(i,j,k))
      tCs_Ssim=tCs_Ssim+sqrt(Cs_Ssim(i,j,k))
      if((k .eq. 1) .AND. ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0))) then  
         arg1=0._rprec
         arg2=0._rprec
      else  
         arg1=(u(i,j,k)+u(i,j,k-1))/2.
         arg2=(v(i,j,k)+v(i,j,k-1))/2.
      end if
      tuw=tuw+w(i,j,k)*arg1
      tvw=tvw+w(i,j,k)*arg2
      tu3=tu3+u(i,j,k)*u(i,j,k)*u(i,j,k)
      tv3=tv3+v(i,j,k)*v(i,j,k)*v(i,j,k)
      tw3=tw3+w(i,j,k)*w(i,j,k)*w(i,j,k)
   end do
   au(i,k)=au(i,k)+(fr)*tu1/Ny
   av(i,k)=av(i,k)+(fr)*tv1/Ny
   aw(i,k)=aw(i,k)+(fr)*tw1/Ny
   ap(i,k)=ap(i,k)+(fr)*tp1/Ny
   adudz(i,k)=adudz(i,k)+(fr)*tdudz/Ny
   advdz(i,k)=advdz(i,k)+(fr)*tdvdz/Ny
   u2(i,k)=u2(i,k)+(fr)*tu2/Ny
   v2(i,k)=v2(i,k)+(fr)*tv2/Ny
   w2(i,k)=w2(i,k)+(fr)*tw2/Ny
   atxx(i,k)=atxx(i,k)+(fr)*ttxx/Ny
   atxz(i,k)=atxz(i,k)+(fr)*ttxz/Ny
   atyy(i,k)=atyy(i,k)+(fr)*ttyy/Ny
   atyz(i,k)=atyz(i,k)+(fr)*ttyz/Ny
   atzz(i,k)=atzz(i,k)+(fr)*ttzz/Ny
   p2(i,k)=p2(i,k)+fr*tp2/Ny
   aCs(i,k)=aCs(i,k)+(fr)*tCs/Ny
   auw(i,k)=auw(i,k)+(fr)*tuw/Ny
   avw(i,k)=avw(i,k)+(fr)*tvw/Ny
   aCs_Ssim(i,k)=aCs_Ssim(i,k)+fr*tCs_Ssim/Ny
   abeta_sgs(i,k)=abeta_sgs(i,k)+fr*Beta_avg(k)
   abetaclip_sgs(i,k)=abetaclip_sgs(i,k)+fr*Betaclip_avg(k)
   u3(i,k)=u3(i,k)+(fr)*tu3/Ny
   v3(i,k)=v3(i,k)+(fr)*tv3/Ny
   w3(i,k)=w3(i,k)+(fr)*tw3/Ny
end do
end do

if (mod(jt,p_count)==0) then
        allocate(avg_out(1:nx,1:(nz_tot-1)));
        call collocate_MPI_averages_N(au,avg_out,20,'u')
        call collocate_MPI_averages_N(av,avg_out,21,'v')
        call collocate_MPI_averages_N(aw,avg_out,22,'w')
        call collocate_MPI_averages_N(ap,avg_out,23,'p')
        call collocate_MPI_averages_N(u2,avg_out,24,'u2')
        call collocate_MPI_averages_N(v2,avg_out,25,'v2')
        call collocate_MPI_averages_N(w2,avg_out,26,'w2')
        call collocate_MPI_averages_N(p2,avg_out,32,'p2')
        call collocate_MPI_averages_N(atxx,avg_out,27,'txx')
        call collocate_MPI_averages_N(atxz,avg_out,28,'txz')
        call collocate_MPI_averages_N(atyy,avg_out,29,'tyy')
        call collocate_MPI_averages_N(atyz,avg_out,30,'tyz')
        call collocate_MPI_averages_N(atzz,avg_out,31,'tzz')
        call collocate_MPI_averages_N(auw,avg_out,33,'uw')
        call collocate_MPI_averages_N(avw,avg_out,34,'vw')
        call collocate_MPI_averages_N(aCs,avg_out,35,'Cs')
        call collocate_MPI_averages_N(adudz,avg_out,36,'dudz')
        call collocate_MPI_averages_N(advdz,avg_out,37,'dvdz')
        call collocate_MPI_averages_N(aCs_Ssim,avg_out,38,'Cs_Ssim')
        call collocate_MPI_averages_N(abeta_sgs,avg_out,39,'beta_sgs')
        call collocate_MPI_averages_N(abetaclip_sgs,avg_out,40,'betaclip_sgs');
        call collocate_MPI_averages_N(u3,avg_out,41,'u3')
        call collocate_MPI_averages_N(v3,avg_out,42,'v3')
        call collocate_MPI_averages_N(w3,avg_out,43,'w3');deallocate(avg_out)

!VK Zero out the outputted averages !!
        au=0._rprec;av=0._rprec;aw=0._rprec;ap=0._rprec;u2=0._rprec;v2=0._rprec
        w2=0._rprec;atxx=0._rprec;atxz=0._rprec;atyy=0._rprec;atyz=0._rprec
        atzz=0._rprec;p2=0._rprec;auw=0._rprec;avw=0._rprec;aCs=0._rprec
        adudz=0._rprec;advdz=0._rprec;aCs_Ssim=0._rprec;abeta_sgs=0._rprec
        abetaclip_sgs=0._rprec;u3=0._rprec;v3=0._rprec;w3=0._rprec;
end if
end subroutine avgslice

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--assumes lun is open and positioned correctly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine checkpoint (lun)

use param,only:nz,S_FLAG
use sim_param,only:u,v,w,RHSx,RHSy,RHSz,theta
use sgsmodule,only:Cs_opt2,F_LM,F_MM,F_QN,F_NN,G_LM,G_MM,G_QN,G_NN,Pr_t
use scalars_module,only:RHS_T,sgs_t3,psi_m

implicit none
integer,intent(in)::lun

if (S_FLAG) then ! WITH SCALARS
   write (lun) u(:,:,1:nz),v(:,:,1:nz),w(:,:,1:nz),theta(:,:,1:nz),   &
               RHSx(:,:,1:nz),RHSy(:,:,1:nz),RHSz(:,:,1:nz),          &
               RHS_T(:,:,1:nz),sgs_t3(:,:,1),psi_m,Cs_opt2,F_LM,F_MM, &
               F_QN,F_NN,G_LM,G_MM,G_QN,G_NN,Pr_t(:,:,1:nz)
else ! No SCALARS
   write (lun) u(:,:,1:nz),v(:,:,1:nz),w(:,:,1:nz),          &
               RHSx(:,:,1:nz),RHSy(:,:,1:nz),RHSz(:,:,1:nz), &
               Cs_opt2,F_LM,F_MM,F_QN,F_NN
end if

end subroutine checkpoint

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--lun_opt gives option for writing to different unit, and is used by inflow_write
!--assumes the unit to write to (lun_default or lun_opt is already
!  open for sequential unformatted writing
!--this routine also closes the unit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_final(jt, lun_opt)
implicit none

integer,intent(in)::jt
integer, intent (in), optional :: lun_opt  !--if present, write to unit lun
integer, parameter :: lun_default = 11
integer::jx,jy,jz
integer :: lun
logical :: opn

if (present (lun_opt)) then
  lun = lun_opt
else
  lun = lun_default
end if

inquire (unit=lun, opened=opn)

if (.not. opn) then
  write (*, *) 'output_final: lun=', lun, ' is not open'
  stop
end if

rewind (lun)

call checkpoint (lun)

close (lun)

if ((cumulative_time) .and. (lun == lun_default)) then
  !--only do this for true final output, not intermediate recording
  open (1, file=fcumulative_time)
  write (1, *) jt_total
  close (1)
end if

end subroutine output_final

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine io_lambda2_out
use sim_param,only:path
implicit none
character(len=24)::fname
call lambda2()
write(fname,'(A13,i6.6,A4)')path//'output/lam-',jt_total,'.out'
open(1,file=fname,form='unformatted')
write(1)nx,ny,nz
write(1)real(lam2)
close(1)
end subroutine io_lambda2_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine lambda2()
use types,only:rprec
use sim_param,only:u,v,w,&
     dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
use param,only:dx,dy,dz
implicit none
real(kind=rprec)::S11,S22,S33,S12,O12,S13,O13,S23,O23,&
     ux,uy,uz,vx,vy,vz,wx,wy,wz
integer::jx,jy,jz
! following used for eispack call...
integer::neis,nmeis,matzeis,ierreis,iv1eis(3)
double precision,dimension(3,3)::aeis,zeis
double precision,dimension(3)::wreis,wieis,fv1eis
double precision::ave

! assignments for eispack call
neis=3
nmeis=3
matzeis=0
ierreis=0
lam2=0._rprec

! at level z=dz/2.  (Where dudz and dvdz are on UVP nodes)
jz=1
do jy=1,ny
do jx=1,nx              
   ux=dudx(jx,jy,1)  ! uvp-node
   uy=dudy(jx,jy,1)  ! uvp-node
   uz=dudz(jx,jy,1)  ! uvp-node
   vx=dvdx(jx,jy,1)  ! uvp-node
   vy=dvdy(jx,jy,1)  ! uvp-node
   vz=dvdz(jx,jy,1)  ! uvp-node 
! special case
   wx=0.5_rprec*(dwdx(jx,jy,1)+dwdx(jx,jy,2))  ! uvp-node
   wy=0.5_rprec*(dwdy(jx,jy,1)+dwdy(jx,jy,2))  ! uvp-node
   wz=dwdz(jx,jy,1)  ! uvp-node
   S11=ux          ! uvp-node
   S12=0.5_rprec*(uy+vx) ! uvp-node
! taken care of with wall stress routine
   S13=0.5_rprec*(uz+wx) ! uvp
   O12=0.5_rprec*(uy-vx) ! w-node
   O13=0.5_rprec*(uz-wx) ! w-node
   S22=vy          ! uvp-node
! taken care of with wall stress routine 
   S23=0.5_rprec*(vz+wy) ! uvp
   O23=0.5_rprec*(vz-wy) ! w-node
   S33=wz          ! uvp-node
   aeis(1,1)=s11*s11+s12*s12+s13*s13-O12*O12-O13*O13
   aeis(1,2)=s11*s12+s12*s22+s13*s23-O13*O23
   aeis(1,3)=s11*s13+s12*s23+s13*s33+O12*O23
   aeis(2,2)=s12*s12+s22*s22+s23*s23-O12*O12-O23*O23
   aeis(2,3)=s12*s13+s22*s23+s23*s33-O12*O13
   aeis(3,3)=s13*s13+s23*s23+s33*s33-O13*O13-O23*O23
   aeis(2,1)=aeis(1,2)
   aeis(3,1)=aeis(1,3)
   aeis(3,2)=aeis(2,3)
  write (*, *) 'rg temporarily removed, sorry'; stop
   if(wreis(1).ge.wreis(2).and.wreis(1).le.wreis(3)) then
      lam2(jx,jy,jz)=real(wreis(1),kind=rprec)
   elseif(wreis(2).ge.wreis(1).and.wreis(2).le.wreis(3)) then
      lam2(jx,jy,jz)=real(wreis(2),kind=rprec)
   else
      lam2(jx,jy,jz)=real(wreis(3),kind=rprec)
   endif
end do
end do
! calculate derivatives/strain on w-nodes
do jz=2,nz-1  
do jy=1,ny
do jx=1,nx              
   ux=0.5_rprec*(dudx(jx,jy,jz) + dudx(jx,jy,jz-1))  ! w-node
   uy=0.5_rprec*(dudy(jx,jy,jz) + dudy(jx,jy,jz-1))  ! w-node
   uz=dudz(jx,jy,jz)  ! w-node
   vx=0.5_rprec*(dvdx(jx,jy,jz) + dvdx(jx,jy,jz-1))  ! w-node
   vy=0.5_rprec*(dvdy(jx,jy,jz) + dvdy(jx,jy,jz-1))  ! w-node
   vz=dvdz(jx,jy,jz)  ! w-node
   wx=dwdx(jx,jy,jz)  ! w-node
   wy=dwdy(jx,jy,jz)  ! w-node
   wz=0.5_rprec*(dwdz(jx,jy,jz) + dwdz(jx,jy,jz-1))  ! w-node
   S11=ux          ! w-node
   S12=0.5_rprec*(uy+vx) ! w-node
   S13=0.5_rprec*(uz+wx) ! w-node
   O12=0.5_rprec*(uy-vx) ! w-node
   O13=0.5_rprec*(uz-wx) ! w-node
   S22=vy          ! w-node
   S23=0.5_rprec*(vz+wy) ! w-node
   O23=0.5_rprec*(vz-wy) ! w-node
   S33=wz          ! w-node
   aeis(1,1)=s11*s11+s12*s12+s13*s13-O12*O12-O13*O13
   aeis(1,2)=s11*s12+s12*s22+s13*s23-O13*O23
   aeis(1,3)=s11*s13+s12*s23+s13*s33+O12*O23
   aeis(2,2)=s12*s12+s22*s22+s23*s23-O12*O12-O23*O23
   aeis(2,3)=s12*s13+s22*s23+s23*s33-O12*O13
   aeis(3,3)=s13*s13+s23*s23+s33*s33-O13*O13-O23*O23
   aeis(2,1)=aeis(1,2)
   aeis(3,1)=aeis(1,3)
   aeis(3,2)=aeis(2,3)
   if(wreis(1).ge.wreis(2).and.wreis(1).le.wreis(3)) then
      lam2(jx,jy,jz)=real(wreis(1),kind=rprec)
   elseif(wreis(2).ge.wreis(1).and.wreis(2).le.wreis(3)) then
      lam2(jx,jy,jz)=real(wreis(2),kind=rprec)
   else
      lam2(jx,jy,jz)=real(wreis(3),kind=rprec)
   endif
end do
end do
end do

print*,'minmax',minval(lam2),maxval(lam2)
end subroutine lambda2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine io_mean_out
implicit none
write(51)real(mean_u),real(mean_u2),real(mean_v),real(mean_v2),&
     real(mean_w),real(mean_w2)
!--the nz/4*3 stuff has to go
mean_u(jx_pls:jx_ple,jy_pls:jy_ple,1:nz/4*3)=0._rprec
mean_u2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz/4*3)=0._rprec
mean_v(jx_pls:jx_ple,jy_pls:jy_ple,1:nz/4*3)=0._rprec
mean_v2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz/4*3)=0._rprec
mean_w(jx_pls:jx_ple,jy_pls:jy_ple,1:nz/4*3)=0._rprec
mean_w2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz/4*3)=0._rprec
end subroutine io_mean_out

subroutine calculate_mean
use sim_param,only:u,v,w
use sgsmodule,only:Cs_opt2,Cs_opt2_avg
implicit none
Cs_opt2_avg(:,:,:)=Cs_opt2_avg(:,:,:)+Cs_opt2(:,:,:)/nwrite
!TS
mean_u(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)=&
     mean_u(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)+&
     u(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)/nwrite
mean_u2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)=&
     mean_u2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)+&
     u(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)**2/nwrite
mean_v(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)=&
     mean_v(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)+&
     v(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)/nwrite
mean_v2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)=&
     mean_v2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)+&
     v(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)**2/nwrite
mean_w(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)=&
     mean_w(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)+&
     w(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)/nwrite
mean_w2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)=&
     mean_w2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)+&
     w(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)**2/nwrite
end subroutine calculate_mean

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine timeseries_spec
use sim_param,only:u,v,w,theta
implicit none
integer::jx,jy,jz,i
if(mod(jt_total,time_spec)==0.and.jt_total.gt.2000)then
jx=NX/8
jy=NY/2+1
jz=NZ/2
endif
end subroutine timeseries_spec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine post_spec(jt_local)
use sim_param,only:path,u,v,w,theta
use param
use fft
implicit none
real(kind=rprec),dimension(nx/2,nz)::spectra_u,spectra_v,spectra_w,&
     spectra_theta
real(kind=rprec),dimension(4,nx/2,nz-1)::spectra_uvwT
real(kind=rprec),dimension(4,nx/2,nz_tot-1)::spectra_uvwT_tot
integer,intent(in)::jt_local
integer::k,jz,z
character(len=64)::fname1,fname2,fname3,fname4
$if ($MPI)
  $define $lbz 0
  integer :: recvcounts(nproc)
  integer :: displs(nproc)
$else
  $define $lbz 1
$endif

write(fname1,'(a,a)') path//'output/spec_x','.dat'
open(82,file=fname1,form='formatted')
do jz=1,nz-1
   z=(jz-0.5_rprec)*dz*z_i
   write(82,*) (real(kx(k,1)/z_i*z),k=1,nx/2)
   call spectrum(u(:, :, jz), spectra_uvwT(1,:,jz))
   call spectrum(v(:, :, jz), spectra_uvwT(2,:,jz))
   call spectrum(w(:, :, jz), spectra_uvwT(3,:,jz))
   call spectrum(theta(:, :, jz), spectra_uvwT(4,:,jz))
enddo
   close(82)
$if ($MPI)
  recvcounts = size (spectra_uvwT)
  displs = coord_of_rank * recvcounts
  call mpi_gatherv (spectra_uvwT(1, 1,1), size (spectra_uvwT), MPI_RPREC,&
                    spectra_uvwT_tot(1, 1, 1), recvcounts, displs,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
$else
  spectra_uvwT_tot=spectra_uvwT
$endif

if((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then  
   write(fname1,'(A,i6.6,A)')path//'output/spec_uvwT_',jt_local,'.bin'
   open(83,file=fname1,form='unformatted')
   write(83) real(spectra_uvwT_tot(:,1:nx/2,:))
   close(83)
end if

end subroutine post_spec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine spectrum(u, spec)
use fft
implicit none      
real(kind=rprec),dimension(ld,ny),intent(in)::u
real(kind=rprec),dimension(nx/2),intent(out)::spec  !--assumes Nyquist is 0

integer::jy,jz,k
real(kind=rprec),dimension(nx)::vel_r,vel_c

integer*8, save :: plan
logical, save :: init = .false.

if (.not. init) then
  call rfftw_f77_create_plan(plan,nx,FFTW_REAL_TO_COMPLEX,FFTW_MEASURE)
  init = .true.
end if

! initialize
spec(:)=0._rprec
do jy=1,ny
   vel_r(:)= u(1:nx,jy)/real(nx,kind=rprec)
! check this normaliztion-part of forward; call the fft
   call rfftw_f77_one(plan,vel_r,vel_c)
! compute magnitudes the 0.5 is the 1/2, all others are taken care of! (except maybe Nyquist)
   spec(1)=spec(1)+0.5*vel_c(1)*vel_c(1)
   do k=2,nx/2
      spec(k)=spec(k)+vel_c(k)*vel_c(k)+vel_c(nx+2-k)*vel_c(nx+2-k)
   end do

   !--assume Nyquist is 0
   !spec(nx/2+1)=spec(nx/2+1)+vel_c(nx/2+1)*vel_c(nx/2+1)
end do
spec(:)=spec(:)/real(Ny,kind=rprec) ! for average over Ny
end subroutine spectrum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine avg_stats ()
use param
use sim_param, only : u, v, w, txz
use fft, only : kx
implicit none

!--choose naming convention that does not conflict with qpost
character (*), parameter :: fubar_avg = 'output/ubar-avg_stats.dat'
character (*), parameter :: fupr2bar_avg = 'output/upr2bar-avg_stats.dat'
character (*), parameter :: fstressbar_avg = 'output/stressbar-avg_stats.dat'
character (*), parameter :: fEozbar_avg = 'output/Eozbar-avg_stats.dat'

integer, parameter :: hdr_len = 256
logical, parameter :: DEBUG = .false.
character (hdr_len) :: Eozbar_hdr

$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif
$if ($MPI)
  integer :: recvcounts(nproc)
  integer :: displs(nproc)
$endif
integer, save :: n_ubar_avg
integer, save :: n_upr2bar_avg
integer, save :: n_stressbar_avg
integer, save :: n_Eozbar_avg
integer :: jz

logical, save :: init = .false.

real (rprec) :: z
real (rprec) :: zu(1, nz_tot-1)
real (rprec) :: kz_z(2, nx/2)
real (rprec), save :: ubar_avg(1, nz_tot-1)      !--1 is <u>
real (rprec), save :: upr2bar_avg(3, nz_tot-1)   !--<u'^2>, <v'^2>, <w'^2>
real (rprec), save :: stressbar_avg(3, nz_tot-1) !--1 is <u'w'>, 2 is <txz>, 3 is <u'w'> + <txz>
real (rprec), save :: Eozbar_avg(1, nx/2, nz_tot-1)  !--E11(k1,z)/z
!--tot is a temp for current stats at nz_tot size
real (rprec), save :: ubar_tot(1, nz_tot-1)      !--1 is <u>
real (rprec), save :: upr2bar_tot(3, nz_tot-1)   !--<u'^2>, <v'^2>, <w'^2>
real (rprec), save :: stressbar_tot(3, nz_tot-1) !--1 is <u'w'>, 2 is <txz>, 3 is <u'w'> + <txz>
real (rprec), save :: Eozbar_tot(1, nx/2, nz_tot-1)  !--E11(k1,z)/z
real (rprec) :: upr(nx, ny), vpr(nx, ny), wpr(nx, ny)
real (rprec) :: ubar(nz-1), vbar(nz-1), wbar(nz-1)
real (rprec) :: upr2bar(3, nz-1)
real (rprec) :: stressbar(3, nz-1)
real (rprec) :: Eozbar(nx/2, nz-1)
!---------------------------------------------------------------------

!--check whether or not to actually do anything
!--motivation for doing this way is that it cleans up interface in main
if (modulo (jt, n_avg_stats) /= 0) goto 001  !--do nothing, exit cleanly

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then  
  if (.not. init) then  !--initialization

    call init_avg (fubar_avg, 1, ubar_avg, n_ubar_avg)
    call init_avg (fupr2bar_avg, 1, upr2bar_avg, n_upr2bar_avg)
    call init_avg (fstressbar_avg, 1, stressbar_avg, n_stressbar_avg) 
    do jz = 1, nz-2
      call init_avg (fEozbar_avg, 2, Eozbar_avg(:, :, jz), n_Eozbar_avg,  &
                     leaveopn='yes')
    end do
    call init_avg (fEozbar_avg, 2, Eozbar_avg(:, :, nz-1), n_Eozbar_avg)

    init = .true.

  end if
end if

!--calculation of current stats
do jz = $lbz, nz-1

  ubar(jz) = sum (u(1:nx, 1:ny, jz)) / (nx * ny)
  vbar(jz) = sum (v(1:nx, 1:ny, jz)) / (nx * ny)
  wbar(jz) = sum (w(1:nx, 1:ny, jz)) / (nx * ny)

  if ( ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) .and.  &
       (jz == 1) ) then
    upr = 0._rprec
    vpr = 0._rprec
    wpr = 0._rprec
  else
    !--see qpost for u/w-node interpolation
    !--convention will be to put up, vp, wp on w-nodes
    upr = 0.5_rprec * (u(1:nx, 1:ny, jz) - ubar(jz) +  &
                       u(1:nx, 1:ny, jz-1) - ubar(jz-1))
    vpr = 0.5_rprec * (v(1:nx, 1:ny, jz) - vbar(jz) +  &
                       v(1:nx, 1:ny, jz-1) - vbar(jz-1))
    wpr = w(1:nx, 1:ny, jz) - wbar(jz)
  end if
 
  upr2bar(1, jz) = sum (upr**2) / (nx * ny)
  upr2bar(2, jz) = sum (vpr**2) / (nx * ny)
  upr2bar(3, jz) = sum (wpr**2) / (nx * ny)

  stressbar(1, jz) = sum (upr * wpr) / (nx * ny) 
  stressbar(2, jz) = sum (txz(1:nx, 1:ny, jz)) / (nx * ny)
  stressbar(3, jz) = sum (stressbar(1:2, jz))

  !--energy spectra
  call spectrum (u(:, :, jz), Eozbar(:, jz))  !--not /z yet
  z = (jz - 0.5_rprec) * dz
  Eozbar(:, jz) = Eozbar(:, jz) / z

end do

!--collect current stats into nz_tot sized arrays
$if ($MPI)

  if (DEBUG) then
    write (*, *) coord, ': ubar(1) = ', ubar(1)
  end if

  recvcounts = size (ubar)
  displs = coord_of_rank * recvcounts 
  call mpi_gatherv (ubar(1), size (ubar), MPI_RPREC,                &
                    ubar_tot(1, 1), recvcounts, displs, MPI_RPREC,  &
                    rank_of_coord(0), comm, ierr)

  recvcounts = size (upr2bar)
  displs = coord_of_rank * recvcounts
  call mpi_gatherv (upr2bar(1, 1), size (upr2bar), MPI_RPREC,          &
                    upr2bar_tot(1, 1), recvcounts, displs, MPI_RPREC,  &
                    rank_of_coord(0), comm, ierr)  

  recvcounts = size (stressbar)
  displs = coord_of_rank * recvcounts
  call mpi_gatherv (stressbar(1, 1), size (stressbar), MPI_RPREC,        &
                    stressbar_tot(1, 1), recvcounts, displs, MPI_RPREC,  &
                    rank_of_coord(0), comm, ierr)

  recvcounts = size (Eozbar)
  displs = coord_of_rank * recvcounts
  call mpi_gatherv (Eozbar(1, 1), size (Eozbar), MPI_RPREC,              &
                    Eozbar_tot(1, 1, 1), recvcounts, displs, MPI_RPREC,  &
                    rank_of_coord(0), comm, ierr)

$else

  ubar_tot(1, :) = ubar
  upr2bar_tot = upr2bar
  stressbar_tot = stressbar
  Eozbar_tot(1, :, :) = Eozbar

$endif

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  !--calculation of cumulative average stats
  ubar_avg = (n_ubar_avg * ubar_avg + ubar_tot) / (n_ubar_avg + 1)
  n_ubar_avg = n_ubar_avg + 1

  upr2bar_avg = (n_upr2bar_avg * upr2bar_avg + upr2bar_tot) /  &
                (n_upr2bar_avg + 1)
  n_upr2bar_avg = n_upr2bar_avg + 1

  stressbar_avg = (n_stressbar_avg * stressbar_avg + stressbar_tot) /  &
                  (n_stressbar_avg + 1)
  n_stressbar_avg = n_stressbar_avg + 1

  Eozbar_avg = (n_Eozbar_avg * Eozbar_avg + Eozbar_tot) / (n_Eozbar_avg + 1)
  n_Eozbar_avg = n_Eozbar_avg + 1

  !--prepare list of z-coordinates
  forall (jz=1:nz_tot-1) zu(1, jz) = (jz - 0.5_rprec) * dz
  !--prepare  header, optional

  !--write out to file
  call write_avg (fubar_avg, n_ubar_avg, zu, ubar_avg)
  call write_avg (fupr2bar_avg, n_upr2bar_avg, zu, upr2bar_avg)
  call write_avg (fstressbar_avg, n_stressbar_avg, zu, stressbar_avg)

  !--this is a bit awkward: maybe just have another routine to do it right
  Eozbar_hdr = 'zone' !--this is for tecplot... 
  kz_z(1, :) = kx(1:nx/2, 1) * zu(1, 1)
  kz_z(2, :) = zu(1, 1)
  call write_avg (fEozbar_avg, n_Eozbar_avg, kz_z, Eozbar_avg(:, :, 1),  &
                  hdr=Eozbar_hdr) 

  do jz = 2, nz_tot - 1
    kz_z(1, :) = kx(1:nx/2, 1) * zu(1, jz)
    kz_z(2, :) = zu(1, jz)
    call write_avg (fEozbar_avg, n_Eozbar_avg, kz_z, Eozbar_avg(:, :, jz),  &
                    hdr=Eozbar_hdr, position='append') 
  end do

end if

001 continue  !--exit cleanly

end subroutine avg_stats

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine init_avg (file_avg, n_ccol, a_avg, n_avg, leaveopn)
implicit none

character (*), intent (in) :: file_avg
integer, intent (in) :: n_ccol  !--num. coord columns: x, y, etc.
real (rprec), intent (out) :: a_avg(:, :)
integer, intent (out) :: n_avg
character (*), optional, intent (in) :: leaveopn
character (128) :: buff
logical :: exst, opn
integer :: j
real (rprec) :: z(n_ccol)

!---------------------------------------------------------------------
inquire (file=file_avg, exist=exst, opened=opn)

if (exst) then

  if (.not. opn) then
    open (1, file=file_avg)
    read (1, '(a)') buff

    if (buff(1:1) == '#') then
      read (buff(2:), *) n_avg
    else
      write (*, *) 'avg_stats: error'
      write (*, *) trim (file_avg), ' does not have expected format on line 1'
      stop  !--need to replace all stops with nice mpi exits
    end if
  end if

  !--skip data header lines here
  do
    read (1, '(a)') buff
    if (trim (buff) == trim (end_hdr_avg)) exit
  end do

  do j = 1, size (a_avg, 2)
    read (1, *) z, a_avg(:, j)  !--z is just placeholder here
  end do

  if (present (leaveopn)) then
    if (leaveopn /= 'yes') close (1)  !--case sensitive here
  else
    close (1)
  end if

else

  n_avg = 0
  a_avg = 0._rprec

end if

end subroutine init_avg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_avg (file_avg, n_avg, x, a_avg, hdr, position)
implicit none

character (*), intent (in) :: file_avg
integer, intent (in) :: n_avg
real (rprec), intent (in) :: x(:, :)  !--coord columns for x, y, etc
real (rprec), intent (in) :: a_avg(:, :)

character (*), optional, intent (in) :: hdr
character (*), optional, intent (in) :: position

character (64) :: r_fmt, fmt
character (32) :: posn

integer :: j

!---------------------------------------------------------------------

!--check sizes compatible
if (size (x, 2) /= size (a_avg, 2)) then
  write (*, *) 'write_avg: error with sizes of x, a_avg'
  stop
end if

if (present (position)) then
  posn = position
else
  posn = 'rewind'
end if

open (1, file=file_avg, position=posn)

if (trim (posn) /= 'append') then  !--case sensitive
  write (1, '(a,i0)') '# ', n_avg  !--not needed when appending
end if

if (present (hdr)) then
  !--write data header, if present
  write (1, '(a)') trim (hdr)
end if

!--write something to indicate end of header, always do this
write (1, '(a)') end_hdr_avg

!--define output format
write (r_fmt, '(2(a,i0))') 'es', precision (1._rprec) + 7,  &
                           '.', precision (1._rprec)
write (fmt, '(a,i0,3a)') '(', size (x, 1) + size (a_avg, 1),  &
                         '(1x,', trim (r_fmt), '))'

!--write to file
do j = 1, size (a_avg, 2)
  write (1, fmt) x(:, j), a_avg(:, j)
end do

close (1)

end subroutine write_avg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine inflow_write ()
use param, only : jt_total, model, jt_start_write, buff_end,  &
                  read_inflow_file, write_inflow_file
use sgsmodule, only : F_MM, F_LM, F_QN, F_NN
use sim_param, only : u, v, w
implicit none

character (*), parameter :: sub = 'inflow_write'
character (*), parameter :: inflow_file = 'output/inflow_BC.out'
character (*), parameter :: field_file = 'output/inflow.vel.out'
character (*), parameter :: MPI_suffix = '.c'

integer, parameter :: lun = 80
integer, parameter :: field_lun = 81

logical, parameter :: DEBUG = .false.

character (64) :: fname

integer, save :: rec = 0
integer :: nrec
integer :: iolen
integer :: iend, iend_w

logical, save :: initialized = .false.
logical :: opn, exst

!---------------------------------------------------------------------

!--option check
if ( read_inflow_file .and. write_inflow_file ) then
  write (*, *) sub // ': cannot have read_inflow_file and write_inflow_file'
  stop
end if

!--check consistency with inflow_cond
iend = floor (buff_end * nx + 1._rprec)
iend_w = modulo (iend - 1, nx) + 1

if (.not. initialized) then

  inquire ( unit=lun, exist=exst, opened=opn )
  if ( .not. exst ) then
    write (*, *) sub // ': lun = ', lun, ' does not exist'
    stop
  end if
  if (opn) then
    write (*, *) sub // ': lun = ', lun, ' is already open'
    stop
  end if

  if ( USE_MPI ) then
      write ( fname, '(a,a,i0)' ) trim (inflow_file), MPI_suffix, coord
  else
      write ( fname, '(a)' ) inflow_file
  end if
  
  inquire ( file=fname, exist=exst, opened=opn )
  if (exst .and. opn) then
    write (*, *) sub // ': file = ', trim (fname), ' is already open'
    stop
  end if
  
  !--figure out the record length
  if ( model.eq.4 ) then
    inquire (iolength=iolen) u(iend_w,:,:), v(iend_w,:,:), w(iend_w,:,:), F_MM(1,:,:), F_LM(1,:,:)
  else if (model.eq.5) then
    inquire (iolength=iolen) u(iend_w,:,:), v(iend_w,:,:), w(iend_w,:,:), F_MM(1,:,:), F_LM(1,:,:), F_QN(1,:,:), F_NN(1,:,:)
  else
    inquire (iolength=iolen) u(iend_w,:,:), v(iend_w,:,:), w(iend_w,:,:)
  end if

  !--always add to existing inflow_file
  !--inflow_file records always start at 1
  if ( exst ) then
      !--figure out the number of records already in file
      call len_da_file (fname, iolen, nrec)
      write (*, *) sub // ': #records in ' // trim (fname) // '= ', nrec
      rec = nrec
  else
      rec = 0
  end if
  
  !--using direct-access file to allow implementation of 'inflow recycling'
  !  more easily
  !--may want to put in some simple checks on ny, nz
  open (unit=lun, file=fname, access='direct', action='write',  &
        recl=iolen)

  initialized = .true.

end if

if (jt_total == jt_start_write) then  !--write entire flow field out
  inquire (unit=field_lun, exist=exst, opened=opn)
  if (exst .and. .not. opn) then
    open (unit=field_lun, file=field_file, form='unformatted')
    call output_final (jt_total, field_lun)
  else
    write (*, *) sub // ': problem opening ' // field_file
    stop
  end if
end if

if (jt_total >= jt_start_write) then
  rec = rec + 1
  if ( model.eq.4 ) then
    write (unit=lun, rec=rec) u(iend_w,:,:), v(iend_w,:,:), w(iend_w,:,:), F_MM(1,:,:), F_LM(1,:,:)
  else if ( model.eq.5) then 
    write (unit=lun, rec=rec) u(iend_w,:,:), v(iend_w,:,:), w(iend_w,:,:), F_MM(1,:,:), F_LM(1,:,:), F_QN(1,:,:), F_NN(1,:,:)
  else
    write (unit=lun, rec=rec) u(iend_w,:,:), v(iend_w,:,:), w(iend_w,:,:)
  end if
  if ( DEBUG ) write (*, *) sub // ': wrote record ', rec
end if

end subroutine inflow_write

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine inflow_read ()
use param, only : model, ny, nz, pi, nsteps, jt_total, buff_end
use sgsmodule, only : FMM_hold, FLM_hold, FQN_hold, FNN_hold
use sim_param, only : u, v, w
implicit none

character (*), parameter :: sub = 'inflow_read'
character (*), parameter :: inflow_file = 'output/inflow_BC.out'
character (*), parameter :: debug_file = 'inflow_read_debug.dat'
character (*), parameter :: MPI_suffix = '.c'

integer, parameter :: lun = 80  !--inflow_write for now
integer, parameter :: lun_DEBUG = 88

integer, parameter :: l_blend = 300  !--length of blending zone (recycling)
                                     !--should correspond to integral scale
                                     !--this is number of t-steps
logical, parameter :: recycle = .false.

logical, parameter :: DEBUG = .false.

character (32) :: fmt
character (64) :: fname

!--check for consistency with sim_param here
!--could define a fortran integer lbz in sim_param, and make it visible
!  here, however, this may have complications elsewhere where the name lbz
!  is used.
$if ( $MPI )
    $define $lbz 0
$else
    $define $lbz 1
$endif

integer :: jy, jz
integer :: iend, iend_w
integer :: i
integer :: iolen
integer, save :: rec
integer, save :: nrec
integer :: recp

logical, save :: init_DEBUG = .false.
logical, save :: initialized = .false.
logical :: exst, opn

real (rprec) :: wgt

real (rprec) :: u_tmp(ny, $lbz:nz), v_tmp(ny, $lbz:nz), w_tmp(ny, $lbz:nz)

!---------------------------------------------------------------------

iend = floor ( buff_end * nx + 1.0_rprec )
iend_w = modulo ( iend - 1, nx ) + 1

if ( .not. initialized ) then

    inquire ( unit=lun, exist=exst, opened=opn )
    if ( .not. exst ) then
        write (*, *) sub // ': lun = ', lun, ' does not exist'
        stop
    end if
    if ( opn ) then
        write (*, *) sub // ': lun = ', lun, ' is already open'
        stop
    end if

    if ( USE_MPI ) then
        write ( fname, '(a,a,i0)' ) trim (inflow_file), MPI_suffix, coord
    else
        write ( fname, '(a)' ) inflow_file
    end if
    
    inquire ( file=fname, exist=exst, opened=opn )
    if ( exst ) then
        if ( opn ) then
            write (*, *) sub // ': file = ', fname, ' is already open'
            stop
        end if
    else
        write (*, *) sub // ': file = ', fname, ' does not exist'
        stop
    end if

    !--can only reach this point if exst and .not. opn
  
    !--figure out the record length
    if ( model.eq.4 ) then
        inquire ( iolength=iolen ) u(iend_w,:,:), v(iend_w,:,:), w(iend_w,:,:), FMM_hold, FLM_hold
    else if ( model.eq.5 ) then
        inquire ( iolength=iolen ) u(iend_w,:,:), v(iend_w,:,:), w(iend_w,:,:), FMM_hold, FMM_hold, FQN_hold, FNN_hold 
    else
        inquire ( iolength=iolen ) u(iend_w, :, :), v(iend_w, :, :), w(iend_w, :, :)
    endif
 
    !--figure out the number of records
    call len_da_file ( fname, iolen, nrec )

    write (*, *) sub // ': number of records = ', nrec

    if ( recycle ) then
        !--check minimum length requirement
        !  checks that there are some points that will be non-blended
        
        if ( 2 * (l_blend - 1) > nrec ) then
            write (*, *) sub // ': ', fname, 'is too short to recycle'
            stop
        end if
    end if

    open ( unit=lun, file=fname, access='direct', action='read',  &
           recl=iolen )

    !--file always starts a record 1, but in continued runs, we may need to
    !  access a record that is not 1 to begin with
    !--actually, with wrap-around of records, it means the reading can start
    !  at any point in the file and be OK
    !--intended use: jt_total = 1 here at start of set of runs reading
    !  from the inflow_file, so the first record read will be record 1
    rec = jt_total - 1

    initialized = .true.

end if
rec = rec + 1
if ( recycle ) then
    rec = modulo ( rec - 1, nrec - l_blend + 1 ) + 1
else
    rec = modulo ( rec - 1, nrec ) + 1
end if

if ( model.eq.4 ) then
    read ( unit=lun, rec=rec ) u(iend_w,:,:), v(iend_w,:,:), w(iend_w,:,:), FMM_hold, FLM_hold
else if ( model.eq.5 ) then
    read ( unit=lun, rec=rec ) u(iend_w,:,:), v(iend_w,:,:), w(iend_w,:,:), FMM_hold, FLM_hold, FQN_hold, FNN_hold
else
    read ( unit=lun, rec=rec ) u(iend_w, :, :), v(iend_w, :, :), w(iend_w, :, :) 
endif
if ( DEBUG ) write (*, *) sub // ' : read record ', rec
    
if ( recycle ) then
    if ( rec < l_blend ) then
        recp = nrec - l_blend + 1 + rec
        wgt = 0.5_rprec * ( 1.0_rprec -                              &
                            cos ( pi * real (rec, rprec) / l_blend ) )
            !--wgt = 0+ when rec = 1
            !  wgt = 1- when rec = l_blend
        read ( unit=lun, rec=recp ) u_tmp, v_tmp, w_tmp
        u(iend_w, :, :) = wgt * u(iend_w, :, :) + (1.0_rprec - wgt) * u_tmp
        v(iend_w, :, :) = wgt * v(iend_w, :, :) + (1.0_rprec - wgt) * v_tmp
        w(iend_w, :, :) = wgt * w(iend_w, :, :) + (1.0_rprec - wgt) * w_tmp
    end if
end if

if ( DEBUG ) then  !--write out slices as an ascii time series
    if ( .not. init_DEBUG ) then
        inquire ( unit=lun_DEBUG, exist=exst, opened=opn )
        if ( exst .and. (.not. opn) ) then
            if ( USE_MPI ) then
                open ( unit=lun_DEBUG, file=debug_file // MPI_suffix )
            else
                open ( unit=lun_DEBUG, file=debug_file )
            end if
        
            write ( lun_DEBUG, '(a)' ) 'variables = "y" "z" "t" "u" "v" "w"'
            write ( lun_DEBUG, '(3(a,i0))' ) 'zone, f=point, i= ', ny,  &
                                             ', j= ', nz,               &
                                             ', k= ', nsteps
        else
            write (*, *) sub // ': problem opening debug file'
            stop
        end if
        init_DEBUG = .true.
    end if

    fmt = '(3(1x,i0),3(1x,es12.5))'
    do jz = 1, nz
        do jy = 1, ny
            write ( lun_DEBUG, fmt ) jy, jz, jt_total, u(iend_w, jy, jz),  &
                                     v(iend_w, jy, jz), w(iend_w, jy, jz)
        end do
    end do
end if

end subroutine inflow_read

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--finds number of records on existing direct-access unformatted file
!--taken from Clive Page's comp.lang.fortran posting (12/16/2003), 
!  under the topic counting number of records in a Fortran direct file
!--minor changes/renaming
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine len_da_file(fname, lenrec, length)
implicit none
character (*), intent(in) :: fname  ! name of existing direct-access file
integer, intent(in)       :: lenrec ! record length (O/S dependent units)
integer, intent(out) :: length      ! number of records.
!
character (1) :: cdummy
integer :: lunit, nlo, nhi, mid, kode
logical :: exists, open
!
! find a free unit on which to open the file
!
do lunit = 99, 1, -1
  !--units to skip (compiler dependent)
  select case (lunit)
    case (5:6)
      !--do nothing
    case default
      inquire(unit=lunit, exist=exists, opened=open)
      if(exists .and. .not. open) exit
  end select
end do
open(unit=lunit, file=fname, access="direct", recl=lenrec, iostat=kode)
if(kode /= 0) then
  print *, 'error in len_da_file: ', trim(fname), ' does not exist'
  return
end if
!
! expansion phase
!
mid = 1
do
  read(lunit, rec=mid, iostat=kode) cdummy
  if(kode /= 0) exit
  mid = 2 * mid
end do
!
! length is between mid/2 and mid, do binary search to refine
!
nlo = mid/2
nhi = mid
do while(nhi - nlo > 1)
  mid = (nlo + nhi) / 2
  read(lunit, rec=mid, iostat=kode) cdummy
  if(kode == 0) then
     nlo = mid
  else
     nhi = mid
  end if
end do
length = nlo
close(unit=lunit)
return
end subroutine len_da_file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine MM_budget_slice
use sim_param,only:path,u,v,w,dudz,dvdz,txx,txz,tyy,tyz,tzz,p,txy,&
         dudx,dudy,dvdx,dvdy,dwdx,dwdy,dwdz,L11t,L22t,L33t,Q11t,Q22t,Q33t
use param,only:dz,p_count,c_count,jt
use sgsmodule,only:Cs_opt2,Cs_Ssim,Beta_avg,Betaclip_avg,dissip,SS11,SS12,SS13,SS22,SS23,SS33
implicit none
integer::i,j,k,kk,kks,jx
$if ($MPI)
  !--this dimensioning adds a ghost layer for finite differences
  !--its simpler to have all arrays dimensioned the same, even though
  !--some components do not need ghost layer
  $define $lbz 0
$else
  $define $lbz 1
$endif
real(kind=rprec),dimension(nx,ny,nz-1),save::u22,v22,w22,uw,vw,u2wp,v2wp,w3p,wtau,utau,vtau,wp,t11,t12,t13,t22,t23,t33,ts11,ts12,ts13,ts22,ts23,ts33,tsgs_1,tsgs_2,tsgs_3
real(kind=rprec),dimension(nx,nz-1),save::ap,au,av,aw,p2,u2,v2,w2,auw,avw,acs,tke,SP,au2w,av2w,aw_i,w2_i,w3_i,we,TR,TR_i,adissip,awp,PC,atau,awtau,autau,avtau,Tsgs,Tpc,ap_c,awp_c,Tdissip_p,aTsgs,aTpc,aTdissip,atke,awe,aSP,aSGS_TKE,tSGS_TKE,aSGS_TKEQ,aSGS_TKEL
real(kind=rprec),dimension(nx,nz-1),save::adudz,advdz,aCs_Ssim,abeta_sgs,abetaclip_sgs
real(kind=rprec),dimension(nx,nz-1),save::atxx,atxz,atyy,atyz,atzz
real(kind=rprec),dimension(nx,nz-1),save::u3,v3,w3
real(kind=rprec)::tu1,tv1,tw1,ttxx,ttxz,ttyy,ttyz,ttzz,tdudz,tdvdz,&
     tu2,tv2,tw2,tp1,tp2,tuw,tvw,tCs,fr,arg1,arg2,tu3,tv3,tw3,tu2w,tv2w,tdissip,twp,ttau,twtau,tutau,tvtau,tp_c,twp_c,arg3,arg4,arg5,arg6,&
     arg22,arg33,arg7,L11_u,L22_u,L33_u,t11w,t12w,t22w,t33w
real(kind=rprec)::tCs_Ssim
real(kind=rprec),dimension(:,:),allocatable::avg_out
real(kind=rprec),dimension($lbz:nz)::ubar_profile,vbar_profile,pbar_profile,au_p,av_p,auw_p,avw_p,adudz_p,advdz_p,au2_p,av2_p,aw2_p,au2w_p,av2w_p,aw3_p,atxz_p,atyz_p,awp_p,awtau_p,autau_p,avtau_p,u_p,v_p,uw_p,vw_p,dudz_p,dvdz_p,u2_p,v2_p,w2_p,u2w_p,v2w_p,w3_p,txz_p,tyz_p,wp_p,wtau_p,utau_p,vtau_p,adissip_p,dissip_p,SS_p,tau_p,atau_p,aSS_p,aw_p,w_p,S11_p,S12_p,S13_p,S22_p,S23_p,S33_p,TS11_p,TS12_p,TS13_p,TS22_p,TS23_p,TS33_p,T11_p,T12_p,T13_p,T22_p,T23_p,T33_p,aS11_p,aS12_p,aS13_p,aS22_p,aS23_p,aS33_p,aTS11_p,aTS12_p,aTS13_p,aTS22_p,aTS23_p,aTS33_p,aT11_p,aT12_p,aT13_p,aT22_p,aT23_p,aT33_p,sig_1u,sig_2v,sig_3w,sig_t,sig_tQ
real(kind=rprec),dimension(nx,ny),save::SD
real(kind=rprec)::S11,S22,S33,S12,S13,S23,&
     ux,uy,uz,vx,vy,vz,wx,wy,wz


fr=(1._rprec/real(p_count,kind=rprec))*real(c_count,kind=rprec)
do k=0,nz-1
        ubar_profile(k)=sum(u(1:nx,1:ny,k))/(nx*ny)
        vbar_profile(k)=sum(v(1:nx,1:ny,k))/(nx*ny)
        pbar_profile(k)=sum(p(1:nx,1:ny,k))/(nx*ny)

end do

do k=1,Nz-1
do i=1,Nx
   tu1=0._rprec;tv1=0._rprec;tw1=0._rprec;tp1=0._rprec
   ttxx=0._rprec;ttxz=0._rprec;ttyy=0._rprec;ttyz=0._rprec
   ttzz=0._rprec;tdudz=0._rprec;tdvdz=0._rprec;tu2=0._rprec
   tv2=0._rprec;tw2=0._rprec;tp2=0._rprec;tuw=0._rprec;tvw=0._rprec
   tCs=0._rprec;tCs_Ssim=0._rprec;tu3=0._rprec;tv3=0._rprec;tw3=0._rprec
   !-- MM     
   tu2w=0._rprec;tv2w=0._rprec;tdissip=0._rprec;twp=0._rprec
   ttau=0._rprec;twtau=0._rprec;tp_c=0._rprec;twp_c=0._rprec;tutau=0._rprec;tvtau=0._rprec
   do j=1,Ny
        if (((k .eq. 1)) .AND. ((.not. USE_MPI) .or. ((USE_MPI) .and.(coord .eq. 0)))) then  
         !arg4=(p(i,j,k)-pbar_profile(k))!-(u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k)))
         arg5=(u(i,j,k)-ubar_profile(k))
         arg6=(v(i,j,k)-vbar_profile(k))
         arg4=(p(i,j,k)-pbar_profile(k)-0.5*(arg5*arg5+arg6*arg6+w(i,j,k)*w(i,j,k)))
        else
         !arg4=(p(i,j,k)-pbar_profile(k)+p(i,j,k-1)-pbar_profile(k-1))/2.!-
         !(u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k))+
         !p(i,j,k-1)-(u(i,j,k-1)*u(i,j,k-1)+v(i,j,k-1)*v(i,j,k-1)+w(i,j,k-1)*w(i,j,k-1)))/2.
         arg5=(u(i,j,k)-ubar_profile(k)+u(i,j,k-1)-ubar_profile(k-1))/2.
         arg6=(v(i,j,k)-vbar_profile(k)+v(i,j,k-1)-vbar_profile(k-1))/2.
         arg4=(p(i,j,k)-pbar_profile(k)-0.5*(arg5*arg5+arg6*arg6+w(i,j,k)*w(i,j,k)) + p(i,j,k-1)-pbar_profile(k-1)-0.5*(arg5*arg5+arg6*arg6+w(i,j,k)*w(i,j,k)) )/2
        endif
     ! end if

      tu1=tu1+u(i,j,k)
      tv1=tv1+v(i,j,k)
      tw1=tw1+w(i,j,k)
      tp1=tp1+p(i,j,k)

      ttxx=ttxx+txx(i,j,k)
      ttxz=ttxz+txz(i,j,k)
      ttyy=ttyy+tyy(i,j,k)
      ttyz=ttyz+tyz(i,j,k)
      ttzz=ttzz+tzz(i,j,k)
      tdudz=tdudz+dudz(i,j,k)
      tdvdz=tdvdz+dvdz(i,j,k)

      tu2=tu2+u(i,j,k)*u(i,j,k)!arg5*arg5
      tv2=tv2+v(i,j,k)*v(i,j,k)!arg6*arg6
      tw2=tw2+w(i,j,k)*w(i,j,k)

      tp2=tp2+p(i,j,k)*p(i,j,k)
      tCs=tCs+sqrt(Cs_opt2(i,j,k))
      tCs_Ssim=tCs_Ssim+sqrt(Cs_Ssim(i,j,k))
      if (((k .eq. Nz-1)) .AND. ((.not. USE_MPI) .or. ((USE_MPI) .and.(coord .eq. nproc-1)))) then  
         arg1=w(i,j,k)!0._rprec
         L11_u=L11t(i,j,k)! 
         L22_u=L22t(i,j,k)! 
         L33_u=L33t(i,j,k)! 
      else
         arg1=(w(i,j,k)+w(i,j,k+1))/2.0 !interpolate into uvp node
         L11_u = (L11t(i,j,k)+L11t(i,j,k+1))/2.0_rprec 
         L22_u = (L22t(i,j,k)+L22t(i,j,k+1))/2.0_rprec
         L33_u = (L33t(i,j,k)+L33t(i,j,k+1))/2.0_rprec
      end if
      if (((k .eq. 1)) .AND. ((.not. USE_MPI) .or. ((USE_MPI) .and.(coord .eq. 0)))) then  
         t11w=txx(i,j,1)
         t12w=txy(i,j,1)
         t22w=tyy(i,j,1)
         t33w=tzz(i,j,1)     
      else
         t11w=(txx(i,j,k)+txx(i,j,k-1))/2.0_rprec
         t12w=(txy(i,j,k)+txy(i,j,k-1))/2.0_rprec
         t22w=(tyy(i,j,k)+tyy(i,j,k-1))/2.0_rprec
         t33w=(tzz(i,j,k)+tzz(i,j,k-1))/2.0_rprec
      endif
      ! MM Computing Strains from velocity derivatives:
!       ux=dudx(i,j,k)
!       uy=dudy(i,j,k)
!       vx=dvdx(i,j,k)
!       vy=dvdy(i,j,k)
!
!       wz=dwdz(i,j,k)

!      if(k .eq. Nz-1) then
!         uz=dudz(i,j,k)
!         vz=dvdz(i,j,k)
!         wx=dwdx(i,j,k)
!         wy=dwdy(i,j,k)
!      else
!         uz=(dudz(i,j,k)+dudz(i,j,k+1))/2.0
!         vz=(dvdz(i,j,k)+dvdz(i,j,k+1))/2.0
!         wx=(dwdx(i,j,k)+dwdx(i,j,k+1))/2.0
!         wy=(dwdy(i,j,k)+dwdy(i,j,k+1))/2.0
!      end if
!        S11=ux
!        S12=0.5_rprec*(uy+vx)
!        S13=0.5_rprec*(uz+wx)
!        S22=vy
!        S23=0.5_rprec*(vz+wy)
!        S33=wz
      if((k .eq. nz-1) .AND. ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1))) then  
         arg2=txz(i,j,k)
         arg3=tyz(i,j,k)
      else  
         arg2=(txz(i,j,k)+txz(i,j,k+1))/2.0_rprec
         arg3=(tyz(i,j,k)+tyz(i,j,k+1))/2.0_rprec
      end if


!      SS11(i,j,k)=S11
!      SS12(i,j,k)=S12
!      SS13(i,j,k)=S13
!      SS22(i,j,k)=S22
!      SS23(i,j,k)=S23
!      SS33(i,j,k)=S33

      tuw=tuw+w(i,j,k)*arg5
      tvw=tvw+w(i,j,k)*arg6

      tu3=tu3+u(i,j,k)*u(i,j,k)*u(i,j,k)
      tv3=tv3+v(i,j,k)*v(i,j,k)*v(i,j,k)
      tw3=tw3+w(i,j,k)*w(i,j,k)*w(i,j,k)
      ! MM : added triple correlation terms required for calculating TKE budget
      tu2w=tu2w+w(i,j,k)*arg5*arg5
      tv2w=tv2w+w(i,j,k)*arg6*arg6

      tp_c=tp_c + p(i,j,k)- 0.5*(u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+arg1*arg1)
      twp_c=twp_c+arg1*(p(i,j,k)-0.5*(u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+arg1*arg1))
      tdissip=tdissip+dissip(i,j,k)
      twp=twp+arg4*w(i,j,k)

      !ttau=ttau+tau(i,j,k)
      twtau=twtau+arg1*tzz(i,j,k) !tzz is on u_node
      tutau=tutau+u(i,j,k)*arg2 !txz on u node 
      tvtau=tvtau+w(i,j,k)*arg3 !tyz on u node

      u22(i,j,k)=u(i,j,k)*u(i,j,k)
      v22(i,j,k)=v(i,j,k)*v(i,j,k)
      w22(i,j,k)=arg1*arg1 !w(i,j,k)*w(i,j,k) on uvp node
      uw(i,j,k)=u(i,j,k)*arg1 !on uvp node
      vw(i,j,k)=v(i,j,k)*arg1 !on uvp node
      u2wp(i,j,k)=u(i,j,k)*u(i,j,k)*arg1 !on uvp node
      v2wp(i,j,k)=v(i,j,k)*v(i,j,k)*arg1 !on uvp node
      w3p(i,j,k)=arg1*arg1*arg1 !on uvp node
      wtau(i,j,k)=arg1*tzz(i,j,k) !tzz is on u node
      utau(i,j,k)=u(i,j,k)*arg2 !txz is on u node as arg2
      vtau(i,j,k)=w(i,j,k)*arg3 !tyz is on u node as arg3
      wp(i,j,k)=arg1*( p(i,j,k)-0.5*(u22(i,j,k)+v22(i,j,k)+w22(i,j,k)) -(L11t(i,j,k)+L22t(i,j,k)+L33t(i,j,k)/(1.375*3.0)) )
     ! wp(i,j,k)=arg1*p(i,j,k)
      SD(i,j) = sqrt(2._rprec*(SS11(i,j,k)**2 + SS22(i,j,k)**2 +&
        SS33(i,j,k)**2 + 2._rprec*(SS12(i,j,k)**2 +&
        SS13(i,j,k)**2 + SS23(i,j,k)**2)))

      t11(i,j,k) = t11w !txx(i,j,k)!2*(SS11(i,j,k))*(SD(i,j))*CS_O(i,j,k)*(Sl(k)**2)
      t12(i,j,k) = t12w !txy(i,j,k)!2*(SS12(i,j,k))*(SD(i,j))*CS_O(i,j,k)*(Sl(k)**2)
      t13(i,j,k) = txz(i,j,k)!arg2!2*(SS13(i,j,k))*(SD(i,j))*CS_O(i,j,k)*(Sl(k)**2)
      t22(i,j,k) = t22w !tyy(i,j,k)!2*(SS22(i,j,k))*(SD(i,j))*CS_O(i,j,k)*(Sl(k)**2)
      t23(i,j,k) = tyz(i,j,k)!arg3!2*(SS23(i,j,k))*(SD(i,j))*CS_O(i,j,k)*(Sl(k)**2)
      t33(i,j,k) = t33w !tzz(i,j,k)!2*(SS33(i,j,k))*(SD(i,j))*CS_O(i,j,k)*(Sl(k)**2)

      ts11(i,j,k)= t11w*SS11(i,j,k)!2*(SS11(i,j,k)**2)*(SD(i,j))*CS_O(i,j,k)*(Sl(k)**2)
      ts12(i,j,k)= t12w*SS12(i,j,k)!2*(SS12(i,j,k)**2)*(SD(i,j))*CS_O(i,j,k)*(Sl(k)**2)
      ts13(i,j,k)= txz(i,j,k)*SS13(i,j,k)!arg2*S13!2*(SS13(i,j,k)**2)*(SD(i,j))*CS_O(i,j,k)*(Sl(k)**2)
      ts22(i,j,k)= t22w*SS22(i,j,k)!2*(SS22(i,j,k)**2)*(SD(i,j))*CS_O(i,j,k)*(Sl(k)**2)
      ts23(i,j,k)= tyz(i,j,k)*SS23(i,j,k)!arg3*S23!2*(SS23(i,j,k)**2)*(SD(i,j))*CS_O(i,j,k)*(Sl(k)**2)
      ts33(i,j,k)= t33w*SS33(i,j,k)!2*(SS33(i,j,k)**2)*(SD(i,j))*CS_O(i,j,k)*(Sl(k)**2) 

       !modify SGS TKE using correct node!!!
        tsgs_1(i,j,k)=(L11_u+L22_u+L33_u)*u(i,j,k)
        tsgs_2(i,j,k)=(L11_u+L22_u+L33_u)*v(i,j,k)
        tsgs_3(i,j,k)=(L11_u+L22_u+L33_u)*arg1 !arg1 is w in u node
        
       !SGS TKE:
!        tsgs_1(i,j,k)=(L11t(i,j,k)+L22t(i,j,k)+L33t(i,j,k))*u(i,j,k)
!        tsgs_2(i,j,k)=(L11t(i,j,k)+L22t(i,j,k)+L33t(i,j,k))*v(i,j,k)
!        tsgs_3(i,j,k)=(L11t(i,j,k)+L22t(i,j,k)+L33t(i,j,k))*w(i,j,k)

   end do

   au(i,k)=au(i,k)+(fr)*tu1/Ny
   av(i,k)=av(i,k)+(fr)*tv1/Ny
   aw(i,k)=aw(i,k)+(fr)*tw1/Ny
   awp(i,k)=awp(i,k)+(fr)*twp/Ny
   awp_c(i,k)=awp_c(i,k)+(fr)*twp_c/Ny
   !atau(i,k)=atau(i,k)+(fr)*ttau/Ny
   awtau(i,k)=awtau(i,k)+(fr)*twtau/Ny
   autau(i,k)=autau(i,k)+(fr)*tutau/Ny
   avtau(i,k)=avtau(i,k)+(fr)*tvtau/Ny
   atxx(i,k)=atxx(i,k)+(fr)*ttxx/Ny
   atxz(i,k)=atxz(i,k)+(fr)*ttxz/Ny
   atyy(i,k)=atyy(i,k)+(fr)*ttyy/Ny
   atyz(i,k)=atyz(i,k)+(fr)*ttyz/Ny
   atzz(i,k)=atzz(i,k)+(fr)*ttzz/Ny
   ap(i,k)=ap(i,k)+(fr)*tp1/Ny
   adudz(i,k)=adudz(i,k)+(fr)*tdudz/Ny
   advdz(i,k)=advdz(i,k)+(fr)*tdvdz/Ny
   u2(i,k)=u2(i,k)+(fr)*tu2/Ny
   v2(i,k)=v2(i,k)+(fr)*tv2/Ny
   w2(i,k)=w2(i,k)+(fr)*tw2/Ny
   au2w(i,j)=au2w(i,k)+(fr)*tu2w/Ny
   av2w(i,j)=av2w(i,k)+(fr)*tv2w/Ny
   p2(i,k)=p2(i,k)+fr*tp2/Ny
   aCs(i,k)=aCs(i,k)+(fr)*tCs/Ny
   auw(i,k)=auw(i,k)+(fr)*tuw/Ny
   avw(i,k)=avw(i,k)+(fr)*tvw/Ny
   aCs_Ssim(i,k)=aCs_Ssim(i,k)+fr*tCs_Ssim/Ny
   abeta_sgs(i,k)=abeta_sgs(i,k)+fr*Beta_avg(k)
   abetaclip_sgs(i,k)=abetaclip_sgs(i,k)+fr*Betaclip_avg(k)
   u3(i,k)=u3(i,k)+(fr)*tu3/Ny
   v3(i,k)=v3(i,k)+(fr)*tv3/Ny
   w3(i,k)=w3(i,k)+(fr)*tw3/Ny
  
   !--MM aw_i = interpolated aw in u-v-p nodes 
    if((k .gt. nz-2)) then !.AND. ((.not. USE_MPI) .or. (USE_MPI .and. coord ==nproc-1))) then
        aw_i(i,k)=aw(i,k)
         w2_i(i,k)=w2(i,k)
         w3_i(i,k)=w3(i,k)
      else
         aw_i(i,k)=(aw(i,k)+aw(i,k+1))/2.0
         w2_i(i,k)=(w2(i,k)+w2(i,k+1))/2.0
         w3_i(i,k)=(w3(i,k)+w3(i,k+1))/2.0
      end if

end do
end do

if((USE_MPI) .and. (coord /= nproc-1)) then
        w_p(Nz)=sum(w(1:nx,1:ny,Nz))/real(nx*ny)
        txz_p(Nz)=sum(txz(1:nx,1:ny,Nz))/real(nx*ny)
        tyz_p(Nz)=sum(tyz(1:nx,1:ny,Nz))/real(nx*ny)
end if


    do kks=1,Nz-1
        k=Nz-kks
        u_p(k)=sum(u(1:nx,1:ny,k))/real(nx*ny)
        v_p(k)=sum(v(1:nx,1:ny,k))/real(nx*ny)
        u2_p(k)=sum(u22(1:nx,1:ny,k))/real(nx*ny)
        v2_p(k)=sum(v22(1:nx,1:ny,k))/real(nx*ny)
        w2_p(k)=sum(w22(1:nx,1:ny,k))/real(nx*ny)
        uw_p(k)=sum(uw(1:nx,1:ny,k))/real(nx*ny)
        vw_p(k)=sum(vw(1:nx,1:ny,k))/real(nx*ny)
        dudz_p(k)=sum(dudz(1:nx,1:ny,k))/real(nx*ny)
        dvdz_p(k)=sum(dvdz(1:nx,1:ny,k))/real(nx*ny)
        u2w_p(k)=sum(u2wp(1:nx,1:ny,k))/real(nx*ny)
        v2w_p(k)=sum(v2wp(1:nx,1:ny,k))/real(nx*ny)
        w3_p(k)=sum(w3p(1:nx,1:ny,k))/real(nx*ny)
        wtau_p(k)=sum(wtau(1:nx,1:ny,k))/real(nx*ny) !u node
        utau_p(k)=sum(utau(1:nx,1:ny,k))/real(nx*ny) !u node
        vtau_p(k)=sum(vtau(1:nx,1:ny,k))/real(nx*ny) !u node
        txz_p(k)=sum(txz(1:nx,1:ny,k))/real(nx*ny) !w node
        tyz_p(k)=sum(tyz(1:nx,1:ny,k))/real(nx*ny) !w node
        wp_p(k)=sum(wp(1:nx,1:ny,k))/real(nx*ny)
        dissip_p(k)=sum(dissip(1:nx,1:ny,k))/real(nx*ny)
        !SS_p(k)=sum(SS(1:nx,1:ny,k))/real(nx*ny)
        !tau_p(k)=sum(tau(1:nx,1:ny,k))/real(nx*ny)
        w_p(k)=sum(w(1:nx,1:ny,k))/real(nx*ny)
        S11_p(k)=sum(SS11(1:nx,1:ny,k))/real(nx*ny)
        S12_p(k)=sum(SS12(1:nx,1:ny,k))/real(nx*ny)
        S13_p(k)=sum(SS13(1:nx,1:ny,k))/real(nx*ny)
        S22_p(k)=sum(SS22(1:nx,1:ny,k))/real(nx*ny)
        S23_p(k)=sum(SS23(1:nx,1:ny,k))/real(nx*ny)
        S33_p(k)=sum(SS33(1:nx,1:ny,k))/real(nx*ny)

        sig_1u(k)=sum(tsgs_1(1:nx,1:ny,k))/real(nx*ny) !on uvp node
        sig_2v(k)=sum(tsgs_2(1:nx,1:ny,k))/real(nx*ny) !on uvp node
        sig_3w(k)=sum(tsgs_3(1:nx,1:ny,k))/real(nx*ny) !on uvp node

        sig_t(k)=(sum(L11t(1:nx,1:ny,k))+sum(L22t(1:nx,1:ny,k))+sum(L33t(1:nx,1:ny,k)))/real(nx*ny)
        sig_tQ(k)=(sum(Q11t(1:nx,1:ny,k))+sum(Q22t(1:nx,1:ny,k))+sum(Q33t(1:nx,1:ny,k)))/real(nx*ny)

        !sig_2(k)=sum(L33(1:nx,1:ny,k))/(nz*ny)
        !sig_3(k)=sum(L33(1:nx,1:ny,k))/(nz*ny)

      if((k .eq. Nz-1) .AND. ((.not. USE_MPI) .or. ((USE_MPI) .and. (coord == nproc-1)))) then
         arg3=w_p(k)
         arg4=txz_p(k)
         arg5=tyz_p(k)
      else
         arg3=(w_p(k)+w_p(k+1))/2.0
         arg4=(txz_p(k)+txz_p(k+1))/2.0
         arg5=(tyz_p(k)+tyz_p(k+1))/2.0
      end if

      if((k .eq. Nz-1) .AND. ((.not. USE_MPI) .or. ((USE_MPI) .and. (coord == nproc-1)))) then
         arg6=dudz_p(k)
         arg7=dvdz_p(k)

      else
         arg6=(dudz_p(k)+dudz_p(k+1))/2.0
         arg7=(dvdz_p(k)+dvdz_p(k+1))/2.0
      end if


        T11_p(k)=sum(t11(1:nx,1:ny,k))/real(nx*ny)
        T12_p(k)=sum(t12(1:nx,1:ny,k))/real(nx*ny)
        T13_p(k)=sum(t13(1:nx,1:ny,k))/real(nx*ny)
        T22_p(k)=sum(t22(1:nx,1:ny,k))/real(nx*ny)
        T23_p(k)=sum(t23(1:nx,1:ny,k))/real(nx*ny)
        T33_p(k)=sum(t33(1:nx,1:ny,k))/real(nx*ny)

        TS11_p(k)=sum(ts11(1:nx,1:ny,k))/real(nx*ny)
        TS12_p(k)=sum(ts12(1:nx,1:ny,k))/real(nx*ny)
        TS13_p(k)=sum(ts13(1:nx,1:ny,k))/real(nx*ny)
        TS22_p(k)=sum(ts22(1:nx,1:ny,k))/real(nx*ny)
        TS23_p(k)=sum(ts23(1:nx,1:ny,k))/real(nx*ny)
        TS33_p(k)=sum(ts33(1:nx,1:ny,k))/real(nx*ny)

       aS11_p(k)=aS11_p(k)+(fr)*S11_p(k)
       aS12_p(k)=aS12_p(k)+(fr)*S12_p(k)
       aS13_p(k)=aS13_p(k)+(fr)*S13_p(k)
       aS22_p(k)=aS22_p(k)+(fr)*S22_p(k)
       aS23_p(k)=aS23_p(k)+(fr)*S23_p(k)
       aS33_p(k)=aS33_p(k)+(fr)*S33_p(k)

       aT11_p(k)=aT11_p(k)+(fr)*T11_p(k)
       aT12_p(k)=aT12_p(k)+(fr)*T12_p(k)
       aT13_p(k)=aT13_p(k)+(fr)*T13_p(k)
       aT22_p(k)=aT22_p(k)+(fr)*T22_p(k)
       aT23_p(k)=aT23_p(k)+(fr)*T23_p(k)
       aT33_p(k)=aT33_p(k)+(fr)*T33_p(k)

       aTS11_p(k)=aTS11_p(k)+(fr)*TS11_p(k)
       aTS12_p(k)=aTS12_p(k)+(fr)*TS12_p(k)
       aTS13_p(k)=aTS13_p(k)+(fr)*TS13_p(k)
       aTS22_p(k)=aTS22_p(k)+(fr)*TS22_p(k)
       aTS23_p(k)=aTS23_p(k)+(fr)*TS23_p(k)
       aTS33_p(k)=aTS33_p(k)+(fr)*TS33_p(k)



       au_p(k)=au_p(k)+(fr)*u_p(k)
       av_p(k)=av_p(k)+(fr)*v_p(k)
       aw_p(k)=aw_p(k)+(fr)*arg3
       auw_p(k)=auw_p(k)+(fr)*uw_p(k)
       avw_p(k)=avw_p(k)+(fr)*vw_p(k)
       adudz_p(k)=adudz_p(k)+(fr)*dudz_p(k)!arg6
       advdz_p(k)=advdz_p(k)+(fr)*dvdz_p(k)!arg7
       au2_p(k)=au2_p(k)+(fr)*u2_p(k)
       av2_p(k)=av2_p(k)+(fr)*v2_p(k)
       aw2_p(k)=aw2_p(k)+(fr)*w2_p(k)
       au2w_p(k)=au2w_p(k)+(fr)*u2w_p(k)
       av2w_p(k)=av2w_p(k)+(fr)*v2w_p(k)
       aw3_p(k)=aw3_p(k)+(fr)*w3_p(k)
       atxz_p(k)=atxz_p(k)+(fr)*arg4
       atyz_p(k)=atyz_p(k)+(fr)*arg5
       awp_p(k)=awp_p(k)+(fr)*wp_p(k)
       awtau_p(k)=awtau_p(k)+(fr)*wtau_p(k)
       autau_p(k)=autau_p(k)+(fr)*utau_p(k)
       avtau_p(k)=avtau_p(k)+(fr)*vtau_p(k)
       adissip_p(k)=adissip_p(k)+(fr)*dissip_p(k)
       !aSS_p(k)=aSS_p(k)+(fr)*SS_p(k)
       !atau_p(k)=atau_p(k)+(fr)*tau_p(k)
   end do

 do k=1,Nz-1
      if((k .eq. Nz-1) .AND. ((.not. USE_MPI) .or. ((USE_MPI) .and. (coord == nproc-1)))) then
         arg4=txz_p(k)
         arg5=tyz_p(k)
      else
         arg4=(txz_p(k)+txz_p(k+1))/2.0_rprec
         arg5=(tyz_p(k)+tyz_p(k+1))/2.0_rprec
      end if
      
 do i=1,Nx 

   tke(i,k)=( u2_p(k) - u_p(k)*u_p(k) + v2_p(k) - v_p(k)*v_p(k) + w2_p(k) )/2.0
   atke(i,k)=atke(i,k)+(fr)*tke(i,k)
      if  ((k .eq. Nz-1)) then
        SP(i,k) =-( (uw_p(k))*dudz_p(k)  + (vw_p(k))*dvdz_p(k))
      else
        SP(i,k) =-(( (uw_p(k))*dudz_p(k)  + (vw_p(k))*dvdz_p(k))/2.0+((uw_p(k))*dudz_p(k+1)  + (vw_p(k))*dvdz_p(k+1))/2.0)
      end if
   aSP(i,k)=aSP(i,k)+(fr)*SP(i,k)
   ! MM - Calculating TKE transport terms :      
    we(i,k)=-0.5*(u2w_p(k)-2*u_p(k)*uw_p(k) +v2w_p(k)-2*v_p(k)*vw_p(k) + w3_p(k))
    awe(i,k)=awe(i,k)+(fr)*we(i,k) !Tt 

    Tsgs(i,k)=-(wtau_p(k) + utau_p(k) -u_p(k)*arg4+ vtau_p(k)-v_p(k)*arg5)
    !Tsgs(i,k)=-(wtau_p(k) + utau_p(k) -u_p(k)*txz_p(k)+ vtau_p(k)-v_p(k)*tyz_p(k))
    aTsgs(i,k)=aTsgs(i,k)+(fr)*Tsgs(i,k) !second term in Tsgs

    Tpc(i,k)=-(wp_p(k))!-aw(i,k)*ap(i,k))  
    aTpc(i,k)=aTpc(i,k)+(fr)*Tpc(i,k)

    Tdissip_p(i,k)=-(  abs(TS11_p(k)-T11_p(k)*S11_p(k)) + abs(TS22_p(k)-T22_p(k)*S22_p(k)) + abs(TS33_p(k)-T33_p(k)*S33_p(k)) + 2*( abs(TS12_p(k)-T12_p(k)*S12_p(k)) + abs(TS23_p(k)-T23_p(k)*S23_p(k)) + abs(TS13_p(k)-T13_p(k)*S13_p(k)) )  )
    aTdissip(i,k)=aTdissip(i,k)+(fr)*Tdissip_p(i,k)

    tSGS_TKE(i,k)= sig_3w(k)-w_p(k)*sig_t(k)
    aSGS_TKE(i,k)=aSGS_TKE(i,k)-(fr)*tSGS_TKE(i,k)/(1.175*3) !first term in Tsgs

    aSGS_TKEL(i,k)=aSGS_TKEL(i,k)+(fr)*sig_t(k)/1.175
    aSGS_TKEQ(i,k)=aSGS_TKEQ(i,k)+(fr)*sig_tQ(k)/3.04

 end do
 end do

!  do jx = 0,nproc-1
!   call output_intmz('wtau',1,jx,200,wtau_p)
!   call output_intmz('vtau',1,jx,200,vtau_p)
!   call output_intmz('utau',1,jx,200,utau_p)
!   call output_intmz('txz',1,jx,200,txz_p)
!   call output_intmz('tyz',1,jx,200,tyz_p)
!  enddo
if (mod(jt,p_count)==0) then

        allocate(avg_out(1:nx,1:(nz_tot-1)));
        call collocate_MPI_averages_N(au,avg_out,20,'u')
        call collocate_MPI_averages_N(av,avg_out,21,'v')
        call collocate_MPI_averages_N(aw,avg_out,22,'w')
        call collocate_MPI_averages_N(ap,avg_out,23,'p')
        call collocate_MPI_averages_N(u2,avg_out,24,'u2')
        call collocate_MPI_averages_N(v2,avg_out,25,'v2')
        call collocate_MPI_averages_N(w2,avg_out,26,'w2')
        call collocate_MPI_averages_N(p2,avg_out,32,'p2')
        call collocate_MPI_averages_N(atxx,avg_out,27,'txx')
        call collocate_MPI_averages_N(atxz,avg_out,28,'txz')
        call collocate_MPI_averages_N(atyy,avg_out,29,'tyy')
        call collocate_MPI_averages_N(atyz,avg_out,30,'tyz')
        call collocate_MPI_averages_N(atzz,avg_out,31,'tzz')
        call collocate_MPI_averages_N(auw,avg_out,33,'uw')
        call collocate_MPI_averages_N(avw,avg_out,34,'vw')
        call collocate_MPI_averages_N(aCs,avg_out,35,'Cs')
        call collocate_MPI_averages_N(adudz,avg_out,36,'dudz')
        call collocate_MPI_averages_N(advdz,avg_out,37,'dvdz')
        call collocate_MPI_averages_N(aCs_Ssim,avg_out,38,'Cs_Ssim')
        call collocate_MPI_averages_N(abeta_sgs,avg_out,39,'beta_sgs')
        call collocate_MPI_averages_N(abetaclip_sgs,avg_out,40,'betaclip_sgs');
        call collocate_MPI_averages_N(u3,avg_out,41,'u3')
        call collocate_MPI_averages_N(v3,avg_out,42,'v3')
        call collocate_MPI_averages_N(w3,avg_out,43,'w3')
        call collocate_MPI_averages_N(atke,avg_out,181,'tke');
        call collocate_MPI_averages_N(aSP,avg_out,182,'SP');
        call collocate_MPI_averages_N(au2w,avg_out,183,'au2w')
        call collocate_MPI_averages_N(av2w,avg_out,184,'av2w')
        call collocate_MPI_averages_N(awe,avg_out,185,'awe')
        call collocate_MPI_averages_N(aTdissip,avg_out,186,'dissip')
        call collocate_MPI_averages_N(aTpc,avg_out,187,'awp')
        call collocate_MPI_averages_N(aw_i,avg_out,188,'aw_i');
        call collocate_MPI_averages_N(aTsgs,avg_out,200,'awtau');
        call collocate_MPI_averages_N(awtau,avg_out,201,'wtau');
        call collocate_MPI_averages_N(avtau,avg_out,202,'vtau');
        call collocate_MPI_averages_N(autau,avg_out,203,'utau');
        call collocate_MPI_averages_N(aSGS_TKE,avg_out,190,'asgs_tke');
        call collocate_MPI_averages_N(aSGS_TKEL,avg_out,191,'asgs_tkeL');
        call collocate_MPI_averages_N(aSGS_TKEQ,avg_out,192,'asgs_tkeQ');
        deallocate(avg_out)
!MM Zero out the outputted averages !!
        au=0._rprec;av=0._rprec;aw=0._rprec;ap=0._rprec;u2=0._rprec;v2=0._rprec
        w2=0._rprec;atxx=0._rprec;atxz=0._rprec;atyy=0._rprec;atyz=0._rprec
        atzz=0._rprec;p2=0._rprec;auw=0._rprec;avw=0._rprec;aCs=0._rprec
        adudz=0._rprec;advdz=0._rprec;aCs_Ssim=0._rprec;abeta_sgs=0._rprec
        abetaclip_sgs=0._rprec;u3=0._rprec;v3=0._rprec;w3=0._rprec;
        au2w=0._rprec;av2w=0._rprec;w2_i=0._rprec;w3_i=0._rprec;aw_i=0._rprec;
        adissip=0._rprec;awp=0._rprec;atau=0._rprec;awtau=0._rprec;ap_c=0._rprec;awp_c=0._rprec;autau=0._rprec;avtau=0._rprec;
        au_p=0._rprec;av_p=0._rprec;aw_p=0._rprec;auw_p=0._rprec;avw_p=0._rprec;adudz_p=0._rprec;advdz_p=0._rprec;au2_p=0._rprec;
        av2_p=0._rprec;aw2_p=0._rprec;au2w_p=0._rprec;av2w_p=0._rprec;aw3_p=0._rprec;atxz_p=0._rprec;
        atyz_p=0._rprec;awp_p=0._rprec;awtau_p=0._rprec;autau_p=0._rprec;avtau_p=0._rprec;adissip_p=0._rprec;
        atau_p=0._rprec;aSS_p=0._rprec;aw_p=0._rprec;
        aS11_p=0._rprec;aS12_p=0._rprec;aS13_p=0._rprec;aS22_p=0._rprec;aS23_p=0._rprec;aS33_p=0._rprec;
        aTS11_p=0._rprec;aTS12_p=0._rprec;aTS13_p=0._rprec;aTS22_p=0._rprec;aTS23_p=0._rprec;aTS33_p=0._rprec;
        aT11_p=0._rprec;aT12_p=0._rprec;aT13_p=0._rprec;aT22_p=0._rprec;aT23_p=0._rprec;aT33_p=0._rprec;
        aTdissip=0._rprec;aTsgs=0._rprec;aTpc=0._rprec;awe=0._rprec;atke=0._rprec;
        aSP=0._rprec;aSGS_TKE=0._rprec;aSGS_TKEQ=0._rprec;aSGS_TKEL=0._rprec;

end if
end subroutine MM_budget_slice




end module io
