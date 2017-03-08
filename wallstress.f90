! For use with staggered grid LES
! JDA, 23 Jan 96
! zo is nondimensionalized, zo1 not!
!--provides txz, tyz, dudz, dvdz at jz=1
! calculates zo_dynamically QL
subroutine wallstress ()
use types,only:rprec
use param,only:jt,jt_total,z_i,nsteps,dz,ld,lh,nx,ny,nz,vonk,lbc_mom,WALL_HOMOG_FLAG,nu_molec,dt,path
use sim_param,only:u,v,dudz,dvdz,txz,tyz
use bottombc,only:zo,psi_m,phi_m,ustar_avg,zo_out,z_os
use test_filtermodule
implicit none
integer::jx,jy
real(kind=rprec),dimension(nx,ny)::u_avg,u_aver,denom
real(kind=rprec),dimension(ld,ny)::u1,v1
real(kind=rprec)::const
logical,parameter :: dynamic_zo = .true.
real(kind=rprec):: zo_surf,zo_mean,zo_var
if (dynamic_zo) then
  if (jt <2) then
   zo = zo_out/z_i!zo_out = 10^-5
   z_os = zo
   zo_mean=sum(zo(:,:))/real(nx*ny)
   zo_var = sum((zo(:,:)-zo_mean)**2)/real(nx*ny)
   else
   zo = (nu_molec /(9._rprec*ustar_avg))/z_i
   z_os = zo
   zo_mean=sum(zo(:,:))/real(nx*ny)
   zo_var = sum((zo(:,:)-zo_mean)**2)/real(nx*ny)
  endif
endif
!
!if (jt <10) then
!print *, 'ustaravg',ustar_avg(60,60),'jt',jt
!print *, 'zo',zo(60,60),'zos',z_os(60,60),'jt',jt
!print *, 'zo2',zo(50,60),'zos2',z_os(50,60),'jt',jt
!endif
!if (jt <10) then
!print *, 'ustaravg',ustar_avg(60,60),'jt',jt
!print *, 'zo',zo(60,60),'zos',z_os(60,60),'jt',jt
!print *, 'zo2',zo(50,60),'zos2',z_os(50,60),'jt',jt
!endif

select case (lbc_mom)

  case ('wall')

    u1=u(:,:,1)
    v1=v(:,:,1)
!   if ((patch_flag .eq. 1) .AND. (num_patch .eq. 1)) then
    if (WALL_HOMOG_FLAG) then
!     calculate u_star in the average sense !!
!     u1=sum(u1(1:nx,1:ny))/float(nx*ny);v1=sum(v1(1:nx,1:ny))/float(nx*ny)
      denom=log(0.5_rprec*dz/zo)-sum(psi_m)/float(nx*ny)
      u_avg=sum(sqrt(u1(1:nx,1:ny)**2+v1(1:nx,1:ny)**2))/float(nx*ny)
      if (jt .eq. nsteps) print *,'USED WALL HOMOG conds in wallstress'
    else
      call test_filter(u1,G_test)
      call test_filter(v1,G_test)
      denom=log(0.5_rprec*dz/zo)-psi_m
      u_avg=sqrt(u1(1:nx,1:ny)**2+v1(1:nx,1:ny)**2)
    end if
!   ustar=u_avg*vonk/denom
    ustar_avg=u_avg*vonk/denom

!if (jt <10) then
!print *, 'ustaravg',ustar_avg(60,60),'jt',jt
!print *, 'psim',psi_m(60,60),'jt',jt
!endif
if (dynamic_zo) then
  if (jt<1) then
   zo = zo_out/z_i!zo_out = 10^-5
   z_os = zo
   zo_mean=sum(zo(:,:))/real(nx*ny)
   zo_var = sum((zo(:,:)-zo_mean)**2)/real(nx*ny)
   else
   zo = (nu_molec /(9._rprec*ustar_avg))/z_i
   z_os = zo
   zo_mean=sum(zo(:,:))/real(nx*ny)
   zo_var = sum((zo(:,:)-zo_mean)**2)/real(nx*ny)
  endif
endif
!if (jt <10) then
!print *, 'ustaravg',ustar_avg(60,60),'jt',jt
!print *, 'zo',zo(60,60),'zos',z_os(60,60),'jt',jt
!print *, 'zo2',zo(50,60),'zos2',z_os(50,60),'jt',jt
!endif

    do jy=1,ny
    do jx=1,nx
       ! const=-(ustar(jx,jy)**2) /u_avg(jx,jy)
       const=-(ustar_avg(jx,jy)**2) /u_avg(jx,jy)
       txz(jx,jy,1)=const *u1(jx,jy)
       tyz(jx,jy,1)=const *v1(jx,jy)
       ! this is as in Moeng 84
       dudz(jx,jy,1)=ustar_avg(jx,jy)/(0.5_rprec*dz*vonK)*u(jx,jy,1)/u_avg(jx,jy)&
           *phi_m(jx,jy)
       dvdz(jx,jy,1)=ustar_avg(jx,jy)/(0.5_rprec*dz*vonK)*v(jx,jy,1)/u_avg(jx,jy)&
           *phi_m(jx,jy)
       dudz(jx,jy,1)=merge(0._rprec,dudz(jx,jy,1),u(jx,jy,1).eq.0._rprec)
       dvdz(jx,jy,1)=merge(0._rprec,dvdz(jx,jy,1),v(jx,jy,1).eq.0._rprec)
    end do
    end do

  case ('stress free')

    txz(:, :, 1) = 0._rprec
    tyz(:, :, 1) = 0._rprec
    dudz(:, :, 1) = 0._rprec
    dvdz(:, :, 1) = 0._rprec

  case default

    write (*, *) 'invalid lbc_mom'
    stop

end select
!it dynamic_zo then write zo 
if (dynamic_zo) then
   open (unit=57, file=path//'output/zo_dynamic.out',status="unknown",position="append")
   write(57,5168) (jt_total+1)*dt,zo_mean,zo_var
 close(57)

5168 format(E14.5,3(1x,E14.5))
endif
end subroutine wallstress
