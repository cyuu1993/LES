module scalars_module2
use types,only:rprec
use param !, jt_global => jt  !--rename to avoid name clashes
! could also modify all routines to access jt from param module, not argument list
use bottombc !makes obukhov functions available
use sim_param,only:u,v,w,theta,q,path,theta_diffusion
use sgsmodule,only: Nu_t
use scalars_module,only: L,wstar,dTdz,dqdz,sgs_t3,sgs_q3 
implicit none

!! --------------------------------------------------------------------
integer,parameter :: average_dim_select_flag=1-(average_dim_num/2) 
! The variable average_dim_select_flag generates the following values based
! on the value of average_dim_num in param.f90 :- 
! a) average_dim_num = 2 : 
!    average_dim_select_flag = 0 ==> average the 3d array over x and y and output the z profile
! b) average_dim_num = 1 : 
!    average_dim_select_flag = 1 ==> average the 3d array over y and output the (x,z) profile
integer, parameter :: dim1_size=average_dim_select_flag*(nx-nz+1)+nz-1
integer, parameter :: dim2_size=average_dim_select_flag*(nz-2)+1
integer, parameter :: dim1_global=average_dim_select_flag*(nx-nz_tot+1)+nz_tot-1
integer, parameter :: dim2_global=average_dim_select_flag*(nz_tot-2)+1
real(kind=rprec):: T_s_min, T_s_max,z_o_min,z_o_max
!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
! Part II of the scalar files .. also look at scalars_module.f90
!contains subroutines:
! ic_scal --- initialize velocity fields and scalars
! patch_or_remote - initialize surface boundary conditions
! scalar_in - read surface data from an external data file (remote-sensed)
! append_zeros_string - append zeros in front of a string - called by scalar_in (NOT USED !!)
! scalar_slice - outputs time-averaged x-z slices of scalar variables
! controlled by c_count & p_count from param; Uses file unit numbers from (36-47)
! obukhov_slice - outputs the obukhov variables (phi,psi, L,wstar);
! called from scalar_slice and toggled by parameter obukhov_output (specified in scalars_module.f)
! DIURNAL_FORCING - sets the surface temperature equal to the mean temperature from a data file via a 
! diurnal forcing over time changing the surface temperature at each time step
! 
! Authored by Vijayant Kumar
! Last modified - June 16, 2005
!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX

contains

!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
subroutine ic_scal()
!subroutine ic_scal(u,v,w,theta,q,sign_wt_s,theta_mean)
!c...Log profile that is modified to flatten at z=z_i
!c .. Modified by Vijayant to put scalars back
! Last modified April 14, 2004

implicit none
real(kind=rprec),dimension(nz)::ubar,vbar,wbar
real(kind=rprec)::rms, noise, arg, arg2,theta_mean
real(kind=rprec)::z,w_star,T_star,q_star,ran3
real(kind=rprec)::z_turb_limit,perturb_height_factor,z_inv
integer::jx,jy,jz,seed,jz_abs
! SKS
real(kind=rprec),dimension(nz_tot-1)::utotinit,vtotinit,wtotinit,Ttotinit
real(kind=rprec),dimension(nz-1)::uinit,vinit,winit,Tinit
$if ($MPI)
  integer::recvcounts(nproc)
  integer::displs(nproc)
$endif
! SKS

       if (lbc .eq. 0) then
           theta_mean=T_init !Set initial temperature
!          theta_mean=T_s_min-2._rprec !T_s_min is dimensional while T_s is non-dimensional
           !theta_mean=T_s_min*T_scale-0.001*T_s_min*T_scale !T_s_min is dimensional while T_s is non-dimensional
! SKS
! What ? T_s_min is dimensional while T_s is non-dimensional ? Why ? What is 0.001 ?
! SKS
           print *,'theta_mean = ',theta_mean
       else
           theta_mean=T_init
           print *,'theta_mean = ',theta_mean
       end if

       if (wt_s .lt. 0._rprec) then
         ! SKS
         perturb_height_factor=0.50_rprec
         z_inv=0.50_rprec*z_i
         ! perturb_height_factor=0.10_rprec
         ! z_inv=0.10_rprec*z_i
         ! z_inv=0.162_rprec*z_i
         ! SKS
       else
         ! SKS
         perturb_height_factor=0.60_rprec
         ! perturb_height_factor=0.3_rprec
         z_inv=0.80_rprec*z_i
         if (passive_scalar) then
         perturb_height_factor=1.00_rprec
         ! perturb_height_factor=0.3_rprec
         z_inv=1.0_rprec*z_i
         ! perturb_height_factor=0.3_rprec
         ! SKS
         endif
       end if

       z_turb_limit=perturb_height_factor*z_i
      
      if (wt_s .eq. 0.0_rprec) then
! Compute the values of w_star etc. using the default value of
! wt_s = 0.06
      w_star=(g/theta_mean*0.06_rprec*z_i)**(1._rprec/3._rprec)
! w_star is of O(1) with z_i=500 and wt_s=0.06
      T_star=0.06_rprec/w_star
      q_star=T_star
      else
      w_star=sign((g/theta_mean*abs(wt_s)*z_i)**(1._rprec/3._rprec),wt_s)
      T_star=wt_s/w_star
      q_star=T_star
      end if

         $if ($MPI)
            print *,'Modified Log Profile for IC for coord = ',coord
         $else
            print *,'Modified Log Profile for IC'
         $endif
       do jz=1,nz
         $if ($MPI)
            z = (coord*(nz-1) + jz - 0.5_rprec) * dz
         $else
            z = (real(jz)-0.5_rprec)*dz
         $endif
!c IC in equilibrium with rough surface (rough dominates in effective zo)
        arg2=z/(sum(zo)/float(nx*ny))
        arg=(1._rprec/vonk)*log(arg2)!-1./(2.*vonk*z_i*z_i)*z*z

        if (coriolis_forcing) then
        ubar(jz)=ug
        vbar(jz)=vg
        wbar(jz)=0._rprec
! Note that ug and vg have already been non-dimensionalized in param.f90
        else
        ubar(jz)=arg
        vbar(jz)=0._rprec
        wbar(jz)=0._rprec
        end if
!C sc: I changed around some parenthesis here

        if (z.gt.(z_turb_limit)) then
           ubar(jz)=ubar(jz-1)
        end if
       end do

  rms = 3._rprec
  do jz=1,nz
  $if ($MPI)
    jz_abs = coord * (nz-1) + jz
    z = (coord * (nz-1) + jz - 0.5_rprec) * dz * z_i
  $else
    jz_abs = jz
    z = (jz-.5_rprec) * dz * z_i
  $endif
  seed = -80 - jz_abs  !--trying to make consistent init for MPI
    do jy=1,ny
      do jx=1,nx
!c...Ran3 returns uniform RV between 0 and 1. (sigma_rv=0.289)
!c...Taking std dev of vel as 1 at all heights

!cVK Note that if you put wt_s = 0 symbolizing neutal conditions
!c u should also put L_z=z_i i.e. the inversion layer height should
!c be equal to the height of the domain and in that case the second
!c part of the subsequent if loop will never execute. This is
!c ensured by putting an OR statement in the if clause, which makes 
!c sure that only the first part of if block is executed and not the
!c block after else
       if (z .LE. z_turb_limit) then
         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
         u(jx,jy,jz)=noise*(1._rprec-z/z_i)*w_star/u_star+ubar(jz)
         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
         v(jx,jy,jz)=noise*(1._rprec-z/z_i)*w_star/u_star+vbar(jz) !noise
         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
         w(jx,jy,jz)=noise*(1._rprec-z/z_i)*w_star/u_star+wbar(jz)
         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
         theta(jx,jy,jz)=(theta_mean+10._rprec*noise*(1-z/z_i)*T_star)/T_scale
         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
         q(jx,jy,jz)=q_mix+50._rprec*noise*(1-z/z_i)*q_star
       else
         u(jx,jy,jz)=ubar(jz)
         v(jx,jy,jz)=vbar(jz)
         w(jx,jy,jz)=wbar(jz)
             if ((z .gt. z_turb_limit) .and. (z .le. z_inv)) then
                  theta(jx,jy,jz)=theta_mean/T_scale
             else
                  theta(jx,jy,jz)=(theta_mean+(z-z_inv)*inv_strength)/T_scale
             end if
         q(jx,jy,jz)=q_mix
        end if
       end do
     end do
  end do

  !...BC for W
  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
    w(1:nx, 1:ny, 1) = 0._rprec
  end if
  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
    w(1:nx, 1:ny, nz) = 0._rprec
  endif

  !...BC for U, V & T
  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
    u(1:nx, 1:ny, nz) = u(1:nx, 1:ny, nz-1)
    v(1:nx, 1:ny, nz) = v(1:nx, 1:ny, nz-1)
    theta(1:nx, 1:ny, nz) = theta(1:nx, 1:ny, nz-1)+dTdz_top/T_scale*z_i*dz
  end if

!VK Display the mean vertical profiles of the initialized variables on the screen
open(unit=44,file=path//'output/init_profiles.dat',status="unknown",position="append")

do jz=1,nz
     $if ($MPI)
       z = (coord*(nz-1) + jz - 0.5_rprec) * dz
     $else
       z = (jz - 0.5_rprec) * dz
     $endif
! SKS
     write(6,7781) jz,z,(sum(u(:,:,jz))/float(nx*ny)),(sum(v(:,:,jz))/&
     float(nx*ny)),(sum(w(:,:,jz))/float(nx*ny)),&
     (sum(theta(:,:,jz))/float(nx*ny))*T_scale
!    write(6,7781) jz,z,(sum(u(:,:,jz))/float(nx*ny))*u_star,(sum(v(:,:,jz))/&
!    float(nx*ny))*u_star,(sum(w(:,:,jz))/float(nx*ny))*u_star,&
!    (sum(theta(:,:,jz))/float(nx*ny))*T_scale
! SKS
     write(44,7781) jz,z,(sum(u(:,:,jz))/float(nx*ny))*u_star,(sum(v(:,:,jz))/&
     float(nx*ny))*u_star,(sum(w(:,:,jz))/float(nx*ny))*u_star,&
     (sum(theta(:,:,jz))/float(nx*ny))*T_scale
end do
close(44)
7781 format('jz, z, ubar, vbar, wbar,T_bar:',(1x,I3,1x,F11.4,1x,F11.4,1x,F11.4,1x,F11.4,1x,F11.4))

! SKS
! This part written just to get the initial profiles 

open(unit=103,file=path//'output/initial_profiles.dat',status="unknown",position="append")

do jz=1,nz-1      
   uinit(jz) = (sum(u(:,:,jz))/float(nx*ny))*u_star
   vinit(jz) = (sum(v(:,:,jz))/float(nx*ny))*u_star
   winit(jz) = (sum(w(:,:,jz))/float(nx*ny))*u_star
   Tinit(jz) = (sum(theta(:,:,jz))/float(nx*ny))*T_scale
end do
      
$if ($MPI) 
  recvcounts = size(uinit)
  displs = coord_of_rank * recvcounts 
  call mpi_gatherv (uinit(1),size(uinit),MPI_RPREC, &
                    utotinit(1),recvcounts,displs,  &
                    MPI_RPREC,rank_of_coord(0),comm,ierr)
  recvcounts = size(vinit)
  displs = coord_of_rank * recvcounts
  call mpi_gatherv (vinit(1),size(vinit),MPI_RPREC, &
                    vtotinit(1),recvcounts,displs,  &
                    MPI_RPREC,rank_of_coord(0),comm,ierr)
  recvcounts = size(winit)
  displs = coord_of_rank * recvcounts
  call mpi_gatherv (winit(1),size(winit),MPI_RPREC, &
                    wtotinit(1),recvcounts,displs,  &
                    MPI_RPREC,rank_of_coord(0),comm,ierr)
  recvcounts = size(Tinit)
  displs = coord_of_rank * recvcounts 
  call mpi_gatherv (Tinit(1),size(Tinit),MPI_RPREC, &
                    Ttotinit(1),recvcounts,displs,  &
                    MPI_RPREC,rank_of_coord(0),comm,ierr)
$else
  utotinit = uinit
  vtotinit = vinit
  wtotinit = winit
  Ttotinit = Tinit
$endif

if((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  do jz=1,nz_tot-1
     write(103,8001) utotinit(jz),vtotinit(jz),wtotinit(jz),Ttotinit(jz)
  end do
  write(103,*)
end if
8001  format(1400(E14.5))
! SKS    

end subroutine ic_scal

!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
subroutine patch_or_remote()
implicit none
real(kind=rprec):: T_s_init_mean,z_o_init_mean
integer :: i,jy
! z_os already defined in scalars_module
! zo,T_s and q_s defined in bottombc
! April 20, 2004 - so far contains both patch_or_remote
! and scalar_in subroutine

if (patch_flag .eq. 1) then
print *, 'Assigning values to patches'
call patches()

z_os(:,:)=(1._rprec/10._rprec)*zo(:,:)

! Calculate minimum and maximum values of T_s and zo for use in
! ic_scal()
   if (lbc .eq. 0) then
   T_s_min=minval(T_s); T_s_max=maxval(T_s)   
   z_o_min=minval(zo); z_o_max=maxval(zo)
   print *,'Min and Max (T_s,zo)',T_s_min,T_s_max,z_o_min,z_o_max
   end if
!c sets temperature field and roughness field at the bottom , x-y plane
!c Added by Vijayant
else if (remote_flag .eq. 1) then
print *, 'Assigning remote-sensed values to the surface'
   call scalar_in() ! Updated T_s and zo loaded from bottombc
   
! Calculate minimum and maximum values of T_s and zo for use in ic_scal()
   if (lbc .eq. 0) then
   T_s_min=minval(T_s); T_s_max=maxval(T_s)   
   z_o_min=minval(zo); z_o_max=maxval(zo)
   print *,'Min and Max (T_s,zo)',T_s_min,T_s_max,z_o_min,z_o_max
   end if

   T_s_init_mean=sum(T_s)/float(nx*ny)
!  z_o_init_mean=sum(zo)/float(nx*ny)
!  perform logarithmic average for zo
   z_o_init_mean=exp(sum(dlog(zo(1:nx,1:ny)))/float(nx*ny))

if (remote_to_patch_flag) then
   call remote_to_patch(T_s,1)
   call remote_to_patch(zo,2)
end if

if (remote_homog_flag .eq. 1) then
print *,'Homogeinizing the remote b.c.s'
   T_s=T_s_init_mean
   zo=z_o_init_mean
end if

! Non-dimensionalize T
   T_s=(T_s)/T_scale
   zo=zo/z_i

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! z_os is the scalar roughness length. Divide momentum
! roughness length by 10 to get scalar roughness length
! data (look in the scalar_in routine above)
z_os(:,:)=(1._rprec/10._rprec)*zo(:,:)

      if (GABLS_test) then
        zo(:,:)=0.1_rprec/z_i
        z_os(:,:)=zo(:,:)
        print *,'setting zo and zo_s for GABLS case = 0.1 m'
      end if
end if

! Write surface bc input to a file
open(unit=57,file=path//'output/surface_bc_input.txt',status="unknown",position="append")
do jy=1,ny
write(57,5168) (T_s(i,jy),i=1,nx)
end do
do jy=1,ny
write(57,5168) (zo(i,jy),i=1,nx)
end do
close(57)

5168     format(1400(E14.5))

end subroutine patch_or_remote

!!!xxxxxxxx--------VIJ---------XXXXXXXXXXXX-------------
!!!xxxxxxxx--------VIJ---------XXXXXXXXXXXX-------------
subroutine scalar_in()
!c This reads in the scalar input from a interpolated scalar
!c file and assigns it to the x-y plane at z=0 i.e. at the ground
!c This for the moment is used to read the momentum roughness, z0
!c and temperature from the USDA remote-sensed data set. The data 
!c can be interpolated using bilinear/cubic/spline interpolation (MATLAB)
!c for the grid size(nx,ny)
!c Authored by Vijayant Kumar
!c Last modified on April 11th, 2004

use bottombc,only:T_s,zo !Load the variables from bottombc and update in here
implicit none
integer:: ii,jj
character(len=6):: suffix2,suffix3
      
       write(suffix2,'(i6.6)') nx ! create a string from nx   
       
      if (coarse_grain_flag) then
      write(suffix3,'(i6.6)') stencil_pts    
 
open(unit=77,file='../interp_data/coarse_grained/interp_temp_cg_'&
//suffix2(4:6)//'pts_'//suffix3(4:6)//'.out',status='unknown')
open(unit=78,file='../interp_data/coarse_grained/interp_z_m_cg_'&
//suffix2(4:6)//'pts_'//suffix3(4:6)//'.out',status='unknown')
     
print *,'interp_temp_cg_'//suffix2(4:6)//'pts_'&
//suffix3(4:6)//'.out loaded from scalar_in.f'
       
      do jj=1,ny
          read(77,5169) (T_s(ii,jj),ii=1,nx)
          read(78,5169) (zo(ii,jj),ii=1,nx)
      enddo
      T_s=T_s+273.15_rprec ! Convert T to Kelvin
      close(77)
      close(78)

       else

open(unit=77,file=path//'interp_data/interp_temp_'//suffix2(4:6)&
//'X'//suffix2(4:6)//'_cubic.out',status='unknown')
open(unit=78,file=path//'interp_data/interp_z_m_'//suffix2(4:6)&
//'X'//suffix2(4:6)//'_cubic.out',status='unknown')
print *,path//'interp_data/interp_temp_'//suffix2(4:6)//'X'//suffix2(4:6)//'_cubic.out'
       
      do jj=1,ny
          read(77,5169) (T_s(ii,jj),ii=1,nx)
          read(78,5169) (zo(ii,jj),ii=1,nx)
      enddo
      print *,'Mean T_s (K) = ',sum(T_s)/float(nx*ny)
      print *,'Mean zo (m) = ',sum(zo)/float(nx*ny)
      close(77)
      close(78)
          T_s=T_s+273.15_rprec ! Convert T to Kelvin
      end if
5169     format(1400(E17.10))
end subroutine scalar_in

!!!xxxxxxxxxxxx--------VIJ-------------XXXXXXXXXXXXXXXXXXX
!!!xxxxxxxxxx-------scalar output subroutine-----XXXXXXXXXXXXXXXXXXXXX

!subroutine scalar_slice()
!!c This is exactly the same like the subroutine avgslice with the
!!c only difference being that it averages the scalar variables
!!c to find the y-averaged instantaneous x-z slices of variables
!!c t,q,sgs_t3,sgs_q3 and their variances such as t2, q2.
!!c It also outputs the average covariance between wt and wq
!!use scalars_module,only: dTdz,dqdz,sgs_t3,sgs_q3 
!!use output_slice,only: collocate_MPI_averages
!implicit none
!integer:: i,j,k
!real(kind=rprec),dimension(nx,nz-1),save:: atheta,t2,q2,asgs_t3,awt
!real(kind=rprec),dimension(nx,nz-1),save:: adTdz,anu_t,t3,var_t
!real(kind=rprec):: ttheta1,tt2,tsgst,twt,tdTdz,arg1,arg2,fr
!real(kind=rprec):: tnu_t,tt3
!real(kind=rprec),dimension(:,:),allocatable:: avg_scalar_out
!
!fr=(1._rprec/float(p_count))*float(c_count)
!
!if (jt .EQ. c_count) then
!  atheta=0._rprec;t2=0._rprec;asgs_t3=0._rprec;awt=0._rprec;adTdz=0._rprec
!  anu_t=0._rprec;t3=0._rprec;var_t=0._rprec
!end if
!
!do k=1,nz-1
!do i=1,nx
!ttheta1=0._rprec;tt2=0._rprec;tsgst=0._rprec;twt=0._rprec;tdTdz=0._rprec
!tnu_t=0._rprec;tt3=0._rprec
!
!do j=1,ny  
!    ttheta1=ttheta1+theta(i,j,k)
!    tsgst=tsgst+sgs_t3(i,j,k)
!    tdTdz=tdTdz+dTdz(i,j,k)
!    tnu_t=tnu_t+Nu_t(i,j,k)
!    tt2 = tt2+ theta(i,j,k)*theta(i,j,k)
!    tt3 = tt3 + theta(i,j,k)*theta(i,j,k)*theta(i,j,k)
!     if((k .eq. 1) .AND. ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0))) then  
!         arg1=0._rprec
!      else  
!         arg1=(theta(i,j,k)+theta(i,j,k-1))/2.
!      end if
!    twt=twt+w(i,j,k)*arg1
!end do
!
!var_t(i,k)=var_t(i,k)+fr*sum((theta(1:nx,1:ny,k)-sum(theta(1:nx,1:ny,k))/(nx*ny))**2)/(nx*ny)
!
!atheta(i,k)=atheta(i,k)+(fr)*ttheta1/ny
!t2(i,k) = t2(i,k) +(fr)*tt2/ny
!t3(i,k) = t3(i,k) +(fr)*tt3/ny
!asgs_t3(i,k)=asgs_t3(i,k)+(fr)*tsgst/ny
!awt(i,k)=awt(i,k)+(fr)*twt/ny
!adTdz(i,k)=adTdz(i,k)+(fr)*tdTdz/ny
!anu_t(i,k)=anu_t(i,k)+(fr)*tnu_t/ny
!end do
!end do
!      
!if (mod(jt,p_count)==0) then
!  allocate(avg_scalar_out(1:nx,1:nz_tot-1));
!  call collocate_MPI_averages_N(atheta,avg_scalar_out,35,'theta')
!  call collocate_MPI_averages_N(t2,avg_scalar_out,36,'t2')
!  call collocate_MPI_averages_N(asgs_t3,avg_scalar_out,37,'sgs_t3')
!  call collocate_MPI_averages_N(awt,avg_scalar_out,38,'wt')
!  call collocate_MPI_averages_N(adTdz,avg_scalar_out,39,'dTdz')
!  call collocate_MPI_averages_N(anu_t,avg_scalar_out,45,'Nu_t')
!  call collocate_MPI_averages_N(t3,avg_scalar_out,46,'t3')
!  call collocate_MPI_averages_N(var_t,avg_scalar_out,47,'var_t');deallocate(avg_scalar_out)
!
!atheta=0._rprec;t2=0._rprec;asgs_t3=0._rprec;awt=0._rprec;adTdz=0._rprec;
!anu_t=0._rprec;t3=0._rprec;
!var_t=0._rprec
!end if
!
!5168      format(1400(E19.10))
!
!end subroutine scalar_slice
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!--The following subroutine does the collocation of the MPI arrays for
!! averaging in avgslice and scalar_slice (in scalars_module2.f90)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine collocate_MPI_averages_N(avg_var_proc,avg_var_tot_domain,file_ind,filename_str)
!use param
!$if ($MPI)
!  integer :: recvcounts(nproc)
!  integer :: displs(nproc)
!$endif
!integer :: ind1,ind2,file_ind
!
!character (*),intent(in) :: filename_str
!character (len=256) :: local_filename
!real(kind=rprec),dimension(dim1_size,dim2_size)::avg_var_proc
!real(kind=rprec),dimension(dim1_global,dim2_global)::avg_var_tot_domain
!
!local_filename=path//'output/aver_'//trim(filename_str)//'.out'
!
!  avg_var_tot_domain=0._rprec
!$if ($MPI)
!  recvcounts = size (avg_var_proc)
!  displs = coord_of_rank * recvcounts 
!  call mpi_gatherv (avg_var_proc(1,1), size (avg_var_proc), MPI_RPREC,                &
!                    avg_var_tot_domain(1, 1), recvcounts, displs, MPI_RPREC,  &
!                    rank_of_coord(0), comm, ierr)
!$else
!  avg_var_tot_domain=avg_var_proc
!$endif
!
!  if((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then  
!  open(file_ind,file=trim(local_filename),status="unknown",position="append")
!      if (average_dim_num .eq. 1) then
!         do ind2=1,nz_tot-1
!          write(file_ind,5168) jt*dt,(avg_var_tot_domain(ind1,ind2),ind1=1,nx)
!         end do
!      else if (average_dim_num .eq. 2) then
!         write(file_ind,5168) jt*dt,(avg_var_tot_domain(ind1,1),ind1=1,nz_tot-1)
!      end if
!  close(file_ind)
!  end if
!
! 5168     format(1400(E14.5))
!end subroutine collocate_MPI_averages_N
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!--The following subroutine does the collocation of the MPI arrays for
!! averaging in avgslice and scalar_slice (in scalars_module2.f90)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine collocate_MPI_averages(avg_var_proc,avg_var_tot_domain,file_ind)
!use param
!$if ($MPI)
!  integer :: recvcounts(nproc)
!  integer :: displs(nproc)
!$endif
!integer :: ind1,ind2,file_ind
!real(kind=rprec),dimension(dim1_size,dim2_size)::avg_var_proc
!real(kind=rprec),dimension(dim1_global,dim2_global)::avg_var_tot_domain
!
!  avg_var_tot_domain=0._rprec
!  $if ($MPI)
!  recvcounts = size (avg_var_proc)
!  displs = coord_of_rank * recvcounts 
!  call mpi_gatherv (avg_var_proc(1,1), size (avg_var_proc), MPI_RPREC,                &
!                    avg_var_tot_domain(1, 1), recvcounts, displs, MPI_RPREC,  &
!                    rank_of_coord(0), comm, ierr)
!$else
!  avg_var_tot_domain=avg_var_proc
!$endif
!
!  if((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then  
!      if (average_dim_num .eq. 1) then
!         do ind2=1,nz_tot-1
!          write(file_ind,5168) jt*dt,(avg_var_tot_domain(ind1,ind2),ind1=1,nx)
!         end do
!      else if (average_dim_num .eq. 2) then
!         write(file_ind,5168) jt*dt,(avg_var_tot_domain(ind1,1),ind1=1,nz_tot-1)
!      end if
!         close(file_ind)
!  end if
!
!5168     format(1400(E14.5))
!end subroutine collocate_MPI_averages
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine scalar_slice()

implicit none
integer:: i,j,k
real(kind=rprec),dimension(nz-1),save::atheta,t2,q2,awt,at_diffusion
real(kind=rprec),dimension(nz-1),save::adTdz,t3,var_t,asgs_t3,aNu_t
real(kind=rprec)::ttheta1,tt2,twt,tdTdz,arg1,arg2,fr
real(kind=rprec)::tt3,tt_diffusion,tsgst3,tNu_t
real(kind=rprec),dimension(:),allocatable::avg_scalar_out

fr=(1._rprec/float(p_count))*float(c_count)

if (jt .EQ. c_count) then
   atheta=0._rprec;t2=0._rprec;awt=0._rprec;adTdz=0._rprec
   t3=0._rprec;var_t=0._rprec;at_diffusion=0._rprec;asgs_T3=0._rprec;
   aNu_t=0._rprec;
end if

do k=1,nz-1
ttheta1=0._rprec;tt2=0._rprec;twt=0._rprec;tdTdz=0._rprec;tt3=0._rprec;
tt_diffusion=0._rprec;tsgst3=0._rprec;tNu_t=0._rprec
   do i=1,nx
   do j=1,ny  
      ttheta1=ttheta1+theta(i,j,k) !*g*EkmanD*1000._rprec/(Ugeo*Ugeo)
      tdTdz=tdTdz+dTdz(i,j,k)
      tt2=tt2+theta(i,j,k)*theta(i,j,k)
      tt3=tt3+theta(i,j,k)*theta(i,j,k)*theta(i,j,k)
      tt_diffusion=tt_diffusion+theta_diffusion(i,j,k)
      tsgst3=tsgst3+sgs_t3(i,j,k)
      tNu_t=tNu_t+Nu_t(i,j,k)
      if((k .eq. 1) .AND. ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0))) then  
         arg1=0._rprec
      else  
         arg1=(theta(i,j,k)+theta(i,j,k-1))/2._rprec
      end if
      twt=twt+w(i,j,k)*arg1
   end do
   end do
   var_t(k)=var_t(k)+fr*sum((100._rprec*theta(1:nx,1:ny,k) &
           -sum(100._rprec*theta(1:nx,1:ny,k))/(real(nx*ny)))**2)/(real(nx*ny))
  ! var_t(k)=var_t(k)+fr*sum((theta(1:nx,1:ny,k) &
  !         -sum(theta(1:nx,1:ny,k))/(nx*ny))**2)/(nx*ny)
   atheta(k)=atheta(k)+(fr)*ttheta1/(nx*ny)
   t2(k) = t2(k)+(fr)*tt2/real(nx*ny)
   t3(k) = t3(k) +(fr)*tt3/real(nx*ny)
   awt(k)=awt(k)+(fr)*twt/(nx*ny)
   adTdz(k)=adTdz(k)+(fr)*tdTdz/(nx*ny)
   at_diffusion(k)=at_diffusion(k)+(fr)*tt_diffusion/(nx*ny)
   asgs_t3(k) = asgs_t3(k)+(fr)*tsgst3/(nx*ny)
   aNu_t(k) = aNu_t(k)+(fr)*tNu_t/(nx*ny)
end do
      
if (mod(jt,p_count)==0) then
  allocate(avg_scalar_out(1:nz_tot-1));
  call collocate_MPI_averages_N(atheta,avg_scalar_out,35,'theta')
  call collocate_MPI_averages_N(t2,avg_scalar_out,36,'t2')
  call collocate_MPI_averages_N(awt,avg_scalar_out,38,'wt')
  call collocate_MPI_averages_N(adTdz,avg_scalar_out,39,'dTdz')
  call collocate_MPI_averages_N(asgs_t3,avg_scalar_out,40,'sgs_t3')
  call collocate_MPI_averages_N(aNu_t,avg_scalar_out,41,'Nu_t')
  call collocate_MPI_averages_N(t3,avg_scalar_out,46,'t3')
  call collocate_MPI_averages_N(var_t,avg_scalar_out,47,'var_t');
  call collocate_MPI_averages_N(at_diffusion,avg_scalar_out,49,'theta_diffusion');
  deallocate(avg_scalar_out)

  atheta=0._rprec;t2=0._rprec;awt=0._rprec;adTdz=0._rprec;t3=0._rprec;var_t=0._rprec;
  at_diffusion=0._rprec;asgs_t3=0._rprec;aNu_t=0._rprec;
end if

5168      format(1400(E14.5))

end subroutine scalar_slice

!--------------------------------------------------------------------------------------------------
subroutine collocate_MPI_averages_N(avg_var_proc,avg_var_tot_domain,file_ind,filename_str)

use param
$if ($MPI)
  integer::recvcounts(nproc)
  integer::displs(nproc)
$endif
integer::ind1,ind2,file_ind

character(*),intent(in)::filename_str
character(len=256)::local_filename
real(kind=rprec),dimension(dim2_size)::avg_var_proc
real(kind=rprec),dimension(dim2_global)::avg_var_tot_domain

local_filename=path//'output/aver_'//trim(filename_str)//'.out'

  avg_var_tot_domain=0._rprec
$if ($MPI)
  recvcounts = size(avg_var_proc)
  displs = coord_of_rank*recvcounts 
  call mpi_gatherv (avg_var_proc(1),size(avg_var_proc),MPI_RPREC,      &
                    avg_var_tot_domain(1),recvcounts,displs,MPI_RPREC, &
                    rank_of_coord(0),comm,ierr)
$else
  avg_var_tot_domain=avg_var_proc
$endif

  if((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then  
  open(file_ind,file=trim(local_filename),status="unknown",position="append")
      if (average_dim_num .eq. 1) then
         do ind2=1,nz_tot-1
          write(file_ind,5168) jt*dt,avg_var_tot_domain(ind2)
         end do
      else if (average_dim_num .eq. 2) then
         write(file_ind,5168) jt*dt,(avg_var_tot_domain(ind1),ind1=1,nz_tot-1)
      end if
  close(file_ind)
  end if

5168     format(1400(F18.11))
end subroutine collocate_MPI_averages_N
!--------------------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------

!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
subroutine timestep_conditions(CFL,visc_stab)
! This subroutine computes CFl and viscous stability and is called every wbase timesteps from main.f90
implicit none
$if ($MPI)
  integer :: recvcounts(nproc)
  integer :: displs(nproc)
! SKS
  real(kind=rprec),dimension(4,nproc)::max_var_tot_domain
! real(kind=rprec),dimension(nproc,4)::max_var_tot_domain
! SKS
$endif

real(kind=rprec) :: delta, u_res_max, nu_max
real(kind=rprec),dimension(1,4) :: max_vels
real(kind=rprec),intent(out) :: CFL, visc_stab
 
$if ($MPI)
  recvcounts = size (max_vels)
  displs = coord_of_rank * recvcounts 
  max_vels(1,1)=maxval(u(1:nx,1:ny,1:nz-1))
  max_vels(1,2)=maxval(v(1:nx,1:ny,1:nz-1))
  max_vels(1,3)=maxval(w(1:nx,1:ny,1:nz-1))
  max_vels(1,4)=maxval(Nu_t(1:nx,1:ny,1:nz-1))
  call mpi_gatherv (max_vels(1,1), size(max_vels), MPI_RPREC,                &
                    max_var_tot_domain(1, 1), recvcounts, displs, MPI_RPREC,  &
                    rank_of_coord(0), comm, ierr)

  ! SKS
  u_res_max=sqrt(maxval(max_var_tot_domain(1,:))**2+maxval(max_var_tot_domain(2,:))**2+&
  maxval(max_var_tot_domain(3,:))**2)
  ! u_res_max=sqrt(maxval(max_var_tot_domain(:,1))**2+maxval(max_var_tot_domain(:,2))**2+&
  ! maxval(max_var_tot_domain(:,3))**2)
  ! SKS
  nu_max=maxval(max_var_tot_domain(:,4))
$else
  u_res_max = sqrt(maxval(u(1:nx,1:ny,1:nz-1)**2+v(1:nx,1:ny,1:nz-1)**2+&
  w(1:nx,1:ny,1:nz-1)**2)) ! don't bother with interp here
  nu_max=maxval(Nu_t(1:nx,1:ny,1:nz-1))
$endif  

if ((.not. USE_MPI) .or. (USE_MPI .and. rank == 0)) then
    delta = min(dx, dy, dz)
    CFL = u_res_max*dt/delta
    visc_stab = dt*nu_max/delta**2
end if

end subroutine timestep_conditions

!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
subroutine ic_scal_GABLS()
!subroutine ic_scal(u,v,w,theta,q,sign_wt_s,theta_mean)
!c...Log profile that is modified to flatten at z=z_i
!c .. Modified by Vijayant to put scalars back
! Last modified April 14, 2004

implicit none
real(kind=rprec),dimension(nz)::ubar,vbar,wbar
real(kind=rprec)::rms, noise, arg, arg2,theta_mean
real(kind=rprec)::z,w_star,T_star,q_star,ran3,perturb_height_factor
integer::jx,jy,jz,seed,jz_abs

       theta_mean=T_init
       perturb_height_factor=50._rprec/z_i
      
      if (wt_s .eq. 0.0_rprec) then
! Compute the values of w_star etc. using the default value of
! wt_s = 0.06
      w_star=(g/theta_mean*0.06_rprec*z_i)**(1._rprec/3._rprec)
! w_star is of O(1) with z_i=500 and wt_s=0.06
      T_star=0.06_rprec/w_star
      q_star=T_star
      else
      w_star=sign((g/theta_mean*abs(wt_s)*z_i)**(1._rprec/3._rprec),wt_s)
      T_star=wt_s/w_star
      q_star=T_star
      end if
        
      $if ($MPI)
         print *,'Modified Log Profile for IC for coord = ',coord
      $else
         print *,'Modified Log Profile for IC'
      $endif

      do jz=1,nz
         $if ($MPI)
            z = (coord*(nz-1) + jz - 0.5_rprec) * dz
         $else
            z=(real(jz)-0.5_rprec)*dz
         $endif
!c IC in equilibrium with rough surface (rough dominates in effective zo)
         arg2=z/(sum(zo)/float(nx*ny))
         arg=(1._rprec/vonk)*log(arg2)!-1./(2.*vonk*z_i*z_i)*z*z
         if (coriolis_forcing) then
         ubar(jz)=ug
         vbar(jz)=vg
         wbar(jz)=0._rprec
! Note that ug and vg have already been non-dimensionalized in param.f90
         else
         ubar(jz)=arg
         vbar(jz)=0._rprec
         wbar(jz)=0._rprec
         end if
!C sc: I changed around some parenthesis here
        if (z.gt.(perturb_height_factor*z_i)) then
        ubar(jz)=ubar(jz-1)
        end if
       end do

  rms = 3._rprec
  do jz=1,nz
  $if ($MPI)
    jz_abs = coord * (nz-1) + jz
    z = (coord * (nz-1) + jz - 0.5_rprec) * dz * z_i
  $else
    jz_abs = jz
    z = (jz-.5_rprec) * dz * z_i
  $endif
  seed = -80 - jz_abs  !--trying to make consistent init for MPI
    do jy=1,ny
      do jx=1,nx
!c...Ran3 returns uniform RV between 0 and 1. (sigma_rv=0.289)
!c...Taking std dev of vel as 1 at all heights

!cVK Note that if you put wt_s = 0 symbolizing neutal conditions
!c u should also put L_z=z_i i.e. the inversion layer height should
!c be equal to the height of the domain and in that case the second
!c part of the subsequent if loop will never execute. This is
!c ensured by putting an OR statement in the if clause, which makes 
!c sure that only the first part of if block is executed and not the
!c block after else

       if (z.le.perturb_height_factor*z_i) then
         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
         u(jx,jy,jz)=noise*(1._rprec-z/z_i)*w_star/u_star+ubar(jz) !noise
         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
         v(jx,jy,jz)=noise*(1._rprec-z/z_i)*w_star/u_star+vbar(jz) !noise
         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
         w(jx,jy,jz)=noise*(1._rprec-z/z_i)*w_star/u_star+wbar(jz)
         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
         theta(jx,jy,jz)=(theta_mean+10._rprec*noise*(1-z/z_i)*T_star)/T_scale
         q(jx,jy,jz)=q_mix
       else
         u(jx,jy,jz)=ubar(jz)
         v(jx,jy,jz)=vbar(jz)
         w(jx,jy,jz)=wbar(jz)
         theta(jx,jy,jz)=(theta_mean+(z-perturb_height_factor*z_i)*inv_strength)/T_scale
         q(jx,jy,jz)=q_mix
        end if
       end do
     end do
  end do

  !...BC for W
  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
    w(1:nx, 1:ny, 1) = 0._rprec
  end if
  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
    w(1:nx, 1:ny, nz) = 0._rprec
  endif

  !...BC for U, V & T
  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
    u(1:nx, 1:ny, nz) = u(1:nx, 1:ny, nz-1)
    v(1:nx, 1:ny, nz) = v(1:nx, 1:ny, nz-1)
    theta(1:nx, 1:ny, nz) = theta(1:nx, 1:ny, nz-1)+dTdz_top/T_scale*z_i*dz
  end if

!VK Display the mean vertical profiles of the initialized variables on the screen
open(unit=44,file=path//'output/init_profiles.dat',status="unknown",position="append")
do jz=1,nz
     $if ($MPI)
       z = (coord*(nz-1) + jz - 0.5_rprec) * dz
     $else
       z = (jz - 0.5_rprec) * dz
     $endif
     write(6,7781) jz,z,(sum(u(:,:,jz))/float(nx*ny))*u_star,(sum(v(:,:,jz))/&
     float(nx*ny))*u_star,(sum(w(:,:,jz))/float(nx*ny))*u_star,&
     (sum(theta(:,:,jz))/float(nx*ny))*T_scale
     write(44,7781) jz,z,(sum(u(:,:,jz))/float(nx*ny))*u_star,(sum(v(:,:,jz))/&
     float(nx*ny))*u_star,(sum(w(:,:,jz))/float(nx*ny))*u_star,&
     (sum(theta(:,:,jz))/float(nx*ny))*T_scale
end do
close(44)
7781 format('jz, z, ubar, vbar, wbar,T_bar:',(1x,I3,1x,F11.4,1x,F11.4,1x,F11.4,1x,F11.4,1x,F11.4))


end subroutine ic_scal_GABLS

!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine budget_TKE_scalar
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! To compute the budget of TKE, we need the following terms:
! <u'w'>*d/dz(<U>), <v'w'>*d/dz(<V>), (g/<T>)*<w'T'>, [0.5*d/dz(<w'e>) where
! e=<u'^2>+<v'^2>+<w'^2>], d/dz(<w'p'>), dissipation (Nu_t*S^2),d/dz(<u'\tau{13}>)+
! d/dz(v'\tau{23}), <\tau{13}>d/dz(<U>)+<\tau{23}>d/dz(<V>)
! Of the eight terms above, we can directly compute terms 1,2,3,8 from the variables
! calculated/outputted in avgslice.f90 and scalar_slice.f90
! So, the rest 4 terms will be computed/outputted here
! Similarly for temperature variance budget, we need
! <w'T'>*d/dz(<T>), d/dz(<w*T^2>) and we already have term 1 from 
! scalar_slice.f90.. so we just need to compute term 2 here 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use sim_param,only:path,u,v,w,theta,p,txz,tyz
use param,only:dz,p_count,c_count,jt,wt_s,u_star,Pr
use sgsmodule,only:dissip
implicit none
$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif

integer::i,j,k,jz_min,ubc_jz
real(kind=rprec),dimension(nx,nz-1),save::awu2,awv2,awt2,auv,awp,awtdtdz
real(kind=rprec),dimension(nx,nz-1),save::autau13,avtau23,adissip
real(kind=rprec)::twu2,twv2,twt2,tuv,twp,arg1,arg2,arg3,arg4,arg5,arg6,fr
real(kind=rprec)::arg7
real(kind=rprec)::tutau13,tvtau23,tdissip,twTdTdz,tTdq3dz
real(kind=rprec),dimension(:,:),allocatable::avg_budget_out
real(kind=rprec),dimension($lbz:nz)::ubar_profile,vbar_profile,Tbar_profile,Tz_p,dq3dz_p
real(kind=rprec),dimension(ld,ny,$lbz:nz):: dtemp,temp,wt,wT2,wTdTdz,Tdq3dz
real(kind=rprec),dimension(nx,nz-1),save::aTdq3dz,MGP,aMGP,Tttheta,aTttheta,Ttsgs,aTtsgs

fr=(1._rprec/real(p_count,kind=rprec))*real(c_count,kind=rprec)

 if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
    jz_min = 2
 else
    jz_min = 1
 end if
 
    ubc_jz = nz-1

 if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
   do j=1,Ny
     do i=1,Nx
      temp(i,j,1) = -1.*(wt_s/u_star/T_scale)
     end do
   end do
 end if

do k=jz_min,nz
  do j=1,Ny
    do i=1,Nx
       temp(i,j,k)=-1./Pr*Nu_t(i,j,k)*dTdz(i,j,k)
    enddo
  enddo
enddo

call DDZ_w(dtemp,temp)

do k=0,nz-1
ubar_profile(k)=sum(u(1:nx,1:ny,k))/(nx*ny)
vbar_profile(k)=sum(v(1:nx,1:ny,k))/(nx*ny)
Tbar_profile(k)=sum(theta(1:nx,1:ny,k))/(nx*ny)
Tz_p(k)=sum(dTdz(1:nx,1:ny,k))/(nx*ny)
dq3dz_p(k)=sum(dtemp(1:nx,1:ny,k))/(nx*ny)
end do

do k=1,Nz-1
do i=1,Nx
   twu2=0._rprec;twv2=0._rprec;twT2=0._rprec;tuv=0._rprec;twp=0._rprec;
   tutau13=0._rprec;tvtau23=0._rprec;tdissip=0._rprec;twTdTdz=0._rprec;tTdq3dz=0._rprec;

   do j=1,Ny
      if((k .eq. 1) .AND. ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0))) then  
         arg1=0._rprec;arg2=0._rprec;arg3=0._rprec;arg4=0._rprec;arg5=0._rprec;arg6=0._rprec
      else  
         arg4=(p(i,j,k)+p(i,j,k-1))/2.
         arg5=(u(i,j,k)-ubar_profile(k)+u(i,j,k-1)-ubar_profile(k-1))/2.
         arg6=(v(i,j,k)-vbar_profile(k)+v(i,j,k-1)-vbar_profile(k-1))/2.
         arg7=(theta(i,j,k)-Tbar_profile(k)+theta(i,j,k-1)-Tbar_profile(k-1))/2.
      end if
      twu2=twu2+w(i,j,k)*arg5*arg5 !directly computes <w(u')^2>
      twv2=twv2+w(i,j,k)*arg6*arg6 !directly computes <w(v')^2>
! Also note that since <w> = 0, there is no need to calculate it as
! <w^3> = <w'^3> and we are outputting <w^3> in avgslice
      twT2=twT2+w(i,j,k)*arg7*arg7
      wT2(i,j,k) = w(i,j,k)*arg7*arg7
      wt(i,j,k) = w(i,j,k)*arg7 
      twTdTdz = twTdtdz +w(i,j,k)*arg7*Tz_p(k)
      wTdTdz(i,j,k) = w(i,j,k)*arg7*Tz_p(k)
      tTdq3dz = tTdq3dz +theta(i,j,k)*dtemp(i,j,k)-Tbar_profile(k)*dq3dz_p(k)
      Tdq3dz(i,j,k) = theta(i,j,k)*dtemp(i,j,k)-Tbar_profile(k)*dq3dz_p(k)
      twp=twp+w(i,j,k)*arg4
! <u'v'> is not as simple as <u'w'> since <u> .ne. 0 whereas <w>=0
! therefore, we directly calculate <u'v'> here
      tuv=tuv+(u(i,j,k)-ubar_profile(k))*(v(i,j,k)-vbar_profile(k))
      tutau13=tutau13+arg5*txz(i,j,k) !computes SGS transport of TKE i.e. <u'\tau_{13}>
      tvtau23=tvtau23+arg6*tyz(i,j,k) !computes SGS transport of TKE i.e. <v'\tau_{13}>
      tdissip=tdissip+dissip(i,j,k) ! outputs dissip calculated in sgs stag..
   end do
   awu2(i,k)=awu2(i,k)+(fr)*twu2/Ny
   awv2(i,k)=awv2(i,k)+(fr)*twv2/Ny
   awT2(i,k)=awT2(i,k)+(fr)*twT2/Ny
   awTdTdz(i,k)=awTdTdz(i,k)+(fr)*twTdTdz/Ny
   aTdq3dz(i,k)=aTdq3dz(i,k)+(fr)*tTdq3dz/Ny
   awp(i,k)=awp(i,k)+(fr)*twp/Ny
   auv(i,k)=auv(i,k)+(fr)*tuv/Ny
   autau13(i,k)=autau13(i,k)+(fr)*tutau13/Ny
   avtau23(i,k)=avtau23(i,k)+(fr)*tvtau23/Ny
   adissip(i,k)=adissip(i,k)+(fr)*tdissip/Ny
end do
end do

     do k = 1,Nz-1
      do i=1,Nx
         MGP(i,k) = -sum(wt(:,:,k))/real(nx*ny)*Tz_p(k)
         aMGP(i,k) = aMGP(i,k) + (fr)*MGP(i,k)
         
         !turbulent transport
         Tttheta(i,k) = -sum(wT2(:,:,k))/real(nx*ny)*0.5
         aTttheta(i,k) = aTttheta(i,k)+(fr)*Tttheta(i,k)
         !resolved to sgs transport
         Ttsgs(i,k) = -sum(Tdq3dz(:,:,k))/real(nx*ny)
         aTtsgs(i,k) = aTtsgs(i,k) + (fr)*Ttsgs(i,k)
       enddo
     enddo

if (mod(jt,p_count)==0) then
  allocate(avg_budget_out(1:nx,1:(nz_tot-1)));
!  call collocate_MPI_averages_N(awu2,avg_budget_out,61,'awu2')
!  call collocate_MPI_averages_N(awv2,avg_budget_out,62,'awv2')
  call collocate_MPI_averages_N(awT2,avg_budget_out,63,'wT2')
  call collocate_MPI_averages_N(awTdTdz,avg_budget_out,64,'wTdTdz')
  call collocate_MPI_averages_N(aTdq3dz,avg_budget_out,65,'Tdq3dz')
  call collocate_MPI_averages_N(aMGP,avg_budget_out,66,'MGP')
  call collocate_MPI_averages_N(aTttheta,avg_budget_out,67,'Tttheta')
  call collocate_MPI_averages_N(aTtsgs,avg_budget_out,68,'Ttsgs')
  
! call collocate_MPI_averages_N(awp,avg_budget_out,64,'awp')
!  call collocate_MPI_averages_N(auv,avg_budget_out,65,'auv')
 ! call collocate_MPI_averages_N(autau13,avg_budget_out,66,'autau13')
!  call collocate_MPI_averages_N(avtau23,avg_budget_out,67,'avtau23')
  call collocate_MPI_averages_N(adissip,avg_budget_out,69,'dissip2');deallocate(avg_budget_out)
!VK Zero out the outputted averages !!
  awu2=0._rprec;awv2=0._rprec;awT2=0._rprec;awp=0._rprec;auv=0._rprec
  autau13=0._rprec;avtau23=0._rprec;adissip=0._rprec;awTdTdz=0._rprec;aTdq3dz=0._rprec;aMGP=0.0_rprec;aTttheta=0._rprec;aTtsgs=0._rprec
end if

end subroutine budget_TKE_scalar

!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
subroutine ic_scal_GABLS_diurnal()
!subroutine ic_scal(u,v,w,theta,q,sign_wt_s,theta_mean)
!c...Log profile that is modified to flatten at z=z_i
!c .. Modified by Vijayant to put scalars back
! Last modified April 14, 2004

implicit none
real(kind=rprec),dimension(nz)::ubar,vbar,wbar,T_bar,tke_bar,tke_sum
real(kind=rprec)::rms, noise, arg, arg2,theta_mean
real(kind=rprec)::z,w_star,T_star,q_star,ran3,perturb_height_factor
real(kind=rprec),dimension(:,:,:),allocatable::wtemp
integer::jx,jy,jz,seed,jz_abs,ii
character (64) :: fname, temp

tke_bar=0._rprec !initialize mean TKE profile as 0
tke_sum=0._rprec !initialize mean TKE profile as 0

      theta_mean=T_init
      perturb_height_factor=800._rprec/z_i
      
      if (wt_s .eq. 0.0_rprec) then
      w_star=(g/theta_mean*0.06_rprec*z_i)**(1._rprec/3._rprec)
! w_star is of O(1) with z_i=500 and wt_s=0.06
      T_star=0.06_rprec/w_star
      q_star=T_star
      else
      w_star=sign((g/theta_mean*abs(wt_s)*z_i)**(1._rprec/3._rprec),wt_s)
      T_star=wt_s/w_star
      q_star=T_star
      end if

         $if ($MPI)
            print *,'Modified Log Profile for IC for coord = ',coord
         $else
            print *,'Modified Log Profile for IC'
         $endif
       do jz=1,nz
         $if ($MPI)
            z = (coord*(nz-1) + jz - 0.5_rprec) * dz
         $else
            z=(real(jz)-0.5_rprec)*dz
         $endif
!c IC in equilibrium with rough surface (rough dominates in effective zo)
        arg2=z/(sum(zo)/float(nx*ny))
        arg=(1._rprec/vonk)*log(arg2)!-1./(2.*vonk*z_i*z_i)*z*z
        if (coriolis_forcing) then
        ubar(jz)=ug
        vbar(jz)=vg
        wbar(jz)=0._rprec
! Note that ug and vg have already been non-dimensionalized in param.f90
        else
        ubar(jz)=arg
        vbar(jz)=0._rprec
        wbar(jz)=0._rprec
        end if
!C sc: I changed around some parenthesis here
        if (z.gt.(perturb_height_factor*z_i)) then
          ubar(jz)=ubar(jz-1)
        end if
       end do
       
5169     format(1400(E17.10))

  do jz=1,nz
  $if ($MPI)
    jz_abs = coord * (nz-1) + jz
    z = (coord * (nz-1) + jz - 0.5_rprec) * dz * z_i
  $else
    jz_abs = jz
    z = (jz-.5_rprec) * dz * z_i
  $endif
  seed = -80 - jz_abs  !--trying to make consistent init for MPI
  if (z .lt. 800._rprec) tke_bar(jz)=0.5_rprec*(1-z/800._rprec)

  rms=(2._rprec*tke_bar(jz)/3._rprec)**0.5_rprec !2/3(1/2(u^2))=1/3(u^2) for each (u,v,w)
    do jy=1,ny
      do jx=1,nx
!c...Ran3 returns uniform RV between 0 and 1. (sigma_rv=0.289)
!c...Taking std dev of vel as 1 at all heights

!cVK Note that if you put wt_s = 0 symbolizing neutal conditions
!c u should also put L_z=z_i i.e. the inversion layer height should
!c be equal to the height of the domain and in that case the second
!c part of the subsequent if loop will never execute. This is
!c ensured by putting an OR statement in the if clause, which makes 
!c sure that only the first part of if block is executed and not the
!c block after else

       if (z.le.perturb_height_factor*z_i) then
         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
         u(jx,jy,jz)=noise/u_star+ubar(jz) !noise
         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
         v(jx,jy,jz)=noise/u_star+vbar(jz) !noise
         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
         w(jx,jy,jz)=noise/u_star+wbar(jz) !noise
         
         call T_pos_gabls(theta_mean,z)
         theta(jx,jy,jz)=theta_mean/T_scale
         q(jx,jy,jz)=q_mix
       else
         u(jx,jy,jz)=ubar(jz)
         v(jx,jy,jz)=vbar(jz)
         w(jx,jy,jz)=wbar(jz)
         call T_pos_gabls(theta_mean,z)
         theta(jx,jy,jz)=theta_mean/T_scale
         q(jx,jy,jz)=q_mix
        end if
       end do
     end do
  end do


  !...BC for W
  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
    w(1:nx, 1:ny, 1) = 0._rprec
  end if
  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
    w(1:nx, 1:ny, nz) = 0._rprec
  endif

  !...BC for U, V & T
  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
    u(1:nx, 1:ny, nz) = u(1:nx, 1:ny, nz-1)
    v(1:nx, 1:ny, nz) = v(1:nx, 1:ny, nz-1)
    theta(1:nx, 1:ny, nz) = theta(1:nx, 1:ny, nz-1)+dTdz_top/T_scale*z_i*dz
  end if

!calculate initial TKE profile
allocate(wtemp(ld,ny,nz)) !wtemp allocated
wtemp=0._rprec
wtemp(:,:,1:nz-1)=0.5_rprec*(w(:,:,1:nz-1)+w(:,:,2:nz));
  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
   wtemp(:,:,nz)=w(:,:,nz);
  end if
do jz=1,nz
   tke_sum(jz)=tke_sum(jz)+0.5_rprec*sum((u(1:nx,1:ny,jz)-sum(u(1:nx,1:ny,jz))/(nx*ny))**2)/(nx*ny)
   tke_sum(jz)=tke_sum(jz)+0.5_rprec*sum((v(1:nx,1:ny,jz)-sum(v(1:nx,1:ny,jz))/(nx*ny))**2)/(nx*ny)
   tke_sum(jz)=tke_sum(jz)+0.5_rprec*sum((wtemp(1:nx,1:ny,jz)-sum(wtemp(1:nx,1:ny,jz))/(nx*ny))**2)/(nx*ny)
end do
deallocate(wtemp)

!VK Display the mean vertical profiles of the initialized variables on the
!screen
    write(fname,'(A,i6.6,A)')path//'output/init_profiles.dat'
    $if ($MPI)
      write (temp, '(".c",i0)') coord
      fname = trim (fname) // temp
    $endif
open(unit=44,file=fname,status="unknown",position="append")
do jz=1,nz
     $if ($MPI)
       z = (coord*(nz-1) + jz - 0.5_rprec) * dz
     $else
       z = (jz - 0.5_rprec) * dz
     $endif
     write(6,7781) jz,z,(sum(u(:,:,jz))/float(nx*ny))*u_star,(sum(v(:,:,jz))/&
     float(nx*ny))*u_star,(sum(w(:,:,jz))/float(nx*ny))*u_star,&
     (sum(theta(:,:,jz))/float(nx*ny))*T_scale,tke_sum(jz)*u_star*u_star
     write(44,7781) jz,z,(sum(u(:,:,jz))/float(nx*ny))*u_star,(sum(v(:,:,jz))/&
     float(nx*ny))*u_star,(sum(w(:,:,jz))/float(nx*ny))*u_star,&
     (sum(theta(:,:,jz))/float(nx*ny))*T_scale,real(tke_sum(jz))*u_star*u_star
end do
close(44)

    write(fname,'(A,i6.6,A)')path//'output/init_profiles_3d.bin'
    $if ($MPI)
      write (temp, '(".c",i0)') coord
      fname = trim (fname) // temp
    $endif
open(unit=44,file=fname,form="unformatted")
write(44) real(u),real(v),real(w),real(theta);close(44)

7781 format('jz, z, ubar, vbar, wbar,T_bar,TKE:',(1x,I3,1x,F9.4,1x,F9.4,1x,F9.4,1x,F9.4,1x,F9.4,1x,F16.10))

end subroutine ic_scal_GABLS_diurnal

!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
subroutine T_pos_gabls(T1,z_loc)
implicit none
real(kind=rprec):: T1,z_loc
if (z_loc<=200._rprec) then
    T1=288._rprec-z_loc*(288._rprec-286._rprec)/200._rprec;
else if (z_loc>200._rprec .AND. z_loc<=850._rprec) then
    T1=286._rprec
else if (z_loc>850._rprec .AND. z_loc<=1000._rprec) then
    T1=286._rprec+(z_loc-850._rprec)*(292._rprec-286._rprec)/(1000._rprec-850._rprec);
else if (z_loc>1000._rprec .AND. z_loc<=2000._rprec) then
    T1=292._rprec+(z_loc-1000._rprec)*(300._rprec-292._rprec)/(2000._rprec-1000._rprec);
else if (z_loc>2000._rprec .AND. z_loc<=3500._rprec) then
    T1=300._rprec+(z_loc-2000._rprec)*(310._rprec-300._rprec)/(3500._rprec-2000._rprec);
else if (z_loc>3500._rprec .AND. z_loc<=4000._rprec) then
    T1=310._rprec+(z_loc-3500._rprec)*(312._rprec-310._rprec)/(4000._rprec-3500._rprec);
else
    print *,'z not contained in the if block range !!!'
    stop
end if
end subroutine T_pos_gabls
!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX

end module scalars_module2
