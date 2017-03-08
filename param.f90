module param
use types,only:rprec
implicit none
$if ($MPI)
!   use mpi SHOULD NOT use 'use mpi' here!!!
  include "mpif.h"
$endif
! implicit none

save

private rprec  !--this is dumb.
public

!--mpi stuff
$if ($MPI)
  $define $MPI_LOGICAL .true.
  $define $NPROC 128
$else
  $define $MPI_LOGICAL .false.
  $define $NPROC 1
$endif

logical, parameter :: USE_MPI = $MPI_LOGICAL

$undefine $MPI_LOGICAL

$if ($MPI)
  integer :: status(MPI_STATUS_SIZE)
$endif

character(*),parameter::path='/glade/scratch/liqi1026/CBL_testmost3/'

!--this stuff must be defined, even if not using MPI
character (8) :: chcoord  !--holds character representation of coord
integer, parameter :: nproc = $NPROC  !--this must be 1 if no MPI
integer :: ierr
integer :: comm
integer :: up, down
integer :: global_rank
integer :: MPI_RPREC, MPI_CPREC
integer :: rank = -1   !--init to bogus (so its defined, even if no MPI)
integer :: coord = -1  !--same here
integer :: rank_of_coord(0:nproc-1), coord_of_rank(0:nproc-1)
!--end mpi stuff

logical,parameter::VERBOSE = .false.  !--prints small stuff to screen
                   !--use DEBUG to write lots of data/files

integer,parameter::nx=256,ny=256,nz=(257-1)/nproc + 1

integer,parameter::nz_tot=(nz-1)*nproc + 1 
integer,parameter::nx2=3*nx/2,ny2=3*ny/2
integer,parameter::lh=nx/2+1,ld=2*lh,lh_big=nx2/2+1,ld_big=2*lh_big

integer, parameter :: iBOGUS = -1234567890  !--NOT a new Apple product
real (rprec), parameter :: BOGUS = -1234567890._rprec

real(rprec),parameter::pi=3.1415926535897932384626433_rprec

real(kind=rprec),parameter::L_x=2.0_rprec*pi,L_y=2.0_rprec*pi
real(rprec),parameter::z_i=1000._rprec, L_z=1500._rprec/nproc
! L_z is not nondimensionalized by z_i yet

! set the aspect ratio of the box, already nondimensional
real(rprec),parameter::dz=L_z/z_i/(nz-1)
real(rprec),parameter::dx=L_x/nx,dy=L_y/ny

integer,parameter::nsteps=20000 ! simulation steps + SCAL_init for 32^3 with dt_dim=0.25
integer,parameter::spectraCALC=0 ! 125000

!--Coriolis stuff; coriol=non-dim coriolis parameter,
! ug=horiz geostrophic vel, vg=transverse geostrophic vel
logical,parameter::coriolis_forcing=.true.

!--fpx3 way to get around sqrt issue below:
$define $ug_dim 8.0
$define $vg_dim 0.0
$define $_rprec _rprec
real (rprec), parameter :: ug_dim = $ug_dim$_rprec
real (rprec), parameter :: vg_dim = $vg_dim$_rprec

$perl: $u_star = sqrt( $ug_dim**2 + $vg_dim**2 )
real(rprec),parameter::u_star= 1.0_rprec !$u_star$_rprec
real(rprec),parameter::Pr=.4_rprec

real(rprec),parameter::dt_dim=0.05_rprec ! this is for unstable wt_s = 0.1, 64 cube case
real(rprec)::dt=dt_dim*u_star/z_i
! real(rprec),parameter::dt=dt_dim*u_star/z_i

real(rprec),parameter::coriol=1.3942923E-04*z_i/u_star !for latitude = 73N ?
real(rprec),parameter::ug=ug_dim/u_star,vg=vg_dim/u_star

real(rprec),parameter::vonk=.4_rprec

! SKS
integer,parameter::c_count=100,p_count=10000!p_count=800 => 5 mins for 64^3;dt=0.003
! integer,parameter::c_count=10,p_count=3600 !p_count=800 => 5 mins for 64^3;dt=0.003
! SKS
integer, parameter :: cs_count = 5  !--tsteps between dynamic Cs updates
logical,parameter::output=.true.
logical, parameter :: use_avgslice = .true.
! The parameter average_dim_num controls the number of dimensions averaged
! while outputting files from avgslice .. 
! Possible values of average_dim_num: 2 => average over x,y,t ; 1 => average over y,t
integer, parameter :: average_dim_num = 1

real(rprec),parameter::nu_molec=1.14e-5_rprec

logical,parameter::use_bldg=.false.
logical,parameter::molec=.false.,sgs=.true.,dns_bc=.false.

! Model type:  1->Smagorinsky; 2->Dynamic; 3->Scale dependent
!              4->Lagrangian scale-sim     5-> Lagragian scale-dep
!              6->new Scale dependent dynamic
!              7->Lagrangian scale dependent for momentum and heat
! Models type: 1->static prandtl,          2->Dynamic
! Cs is the Smagorinsky Constant
! Co and nnn are used in the mason model for smagorisky coeff
integer,parameter::model=5,models=1,nnn=2,BETA_FLAG=1

! calc_lag_Ssim enables calculation of a Cs based on lag Ssim
! using BETA = 1 while using lag sdep for the real Cs. The Cs calc using
! calc_lag_ssim is only used for plotting purposes
! CAUTION : use only with model=5 and only when required (e.g. diurnal runs)
logical,parameter::calc_lag_Ssim=.false.

real(kind=rprec),parameter::Co=0.16_rprec

! Test filter type: 1->cut off 2->Gaussian 3->Top-hat
integer,parameter::ifilter=1

! ubc: upper boundary condition: ubc=0 stress free lid, ubc=1 sponge
! damping method = 1: use Nieuwstadt's method, = 2: use Raleigh damping method
integer,parameter::ubc=1,damping_method=2

character (*), parameter :: lbc_mom = 'wall' !--'wall', 'stress free'

! prescribed inflow: constant or read from file
logical,parameter::inflow=.false.

real (rprec), parameter :: buff_end = 1._rprec    ! position of right end of buffer region, as a fraction of L_x
real (rprec), parameter :: buff_len = 1.0_rprec/16.0_rprec ! length of buffer region as a fraction of L_x
real (rprec), parameter :: face_avg = 1.5_rprec

! SKS
logical, parameter :: read_inflow_file = .false.
! logical, parameter :: read_inflow_file = .true.
! SKS

logical, parameter :: write_inflow_file = .false. !--records at position jx_s
integer, parameter :: jt_start_write = 15000

! forcing along top and bottom bdrys if inflow is true and force_top_bot is true, 
! then the top & bottom velocities are forced to the inflow velocity
logical, parameter :: force_top_bot = .false.

logical, parameter :: use_mean_p_force = .false. ! .not.inflow
real (rprec), parameter :: mean_p_force = 1._rprec * z_i/L_z/nproc  !--usually just z_i/L_z

integer :: jt        ! global time-step counter
integer :: jt_total  ! used for cumulative time (see io module)

! time advance parameters (AB2)
real (rprec), parameter :: tadv1 = 1.5_rprec, tadv2 = 1._rprec - tadv1

!------xxxxxxxxx--SCALARS_PARAMETERS--xxxxxxxxx---------------
! SKS
! logical,parameter :: S_FLAG=.true.,spatial_flux_flag=.FALSE.,OB_homog_flag=.TRUE.,WALL_HOMOG_FLAG=.FALSE.
! logical,parameter :: S_FLAG=.TRUE.,spatial_flux_flag=.FALSE.,OB_homog_flag=.TRUE.,WALL_HOMOG_FLAG=.TRUE.
!  logical,parameter :: S_FLAG=.TRUE.,spatial_flux_flag=.FALSE.,OB_homog_flag=.FALSE.,WALL_HOMOG_FLAG=.false.
  logical,parameter :: S_FLAG=.TRUE.,spatial_flux_flag=.FALSE.,OB_homog_flag=.FALSE.,WALL_HOMOG_FLAG=.true.
! SKS

integer,parameter::DYN_init=1, SCAL_init=1

! lbc=0: prescribed surface temperature, lbc=1 prescribed surface flux
! SKS
integer,parameter :: lbc=1, patch_flag=1, remote_flag=0,remote_homog_flag=0,remote_flux_homog_flag=0
! integer,parameter :: lbc=0, patch_flag=1, remote_flag=0,remote_homog_flag=0,remote_flux_homog_flag=0
! SKS
logical,parameter :: remote_to_patch_flag=.FALSE.! create a 2 patch T_s field using the remote-sensed data
! The corresponding remote_to_patch subroutine is in bottombc.f90
integer,parameter :: diurnal_forcing_flag=0, no_days=1
logical,parameter :: jan_diurnal_run=.false.,ug_diurnal_variation=.false.
logical,parameter :: GABLS_diurnal_test=.false.
logical,parameter :: initu=.false.,initsc=.false.,inilag=.true.,interp=.FALSE.
! initu   = true to read from a file; false to create with random noise
! initlag = true to initialize cs, FLM & FMM; false to read from vel.out

! Added a new parameter - passive_scalar for passive scalars with bldngs
logical,parameter :: passive_scalar=.false.,GABLS_test=.false.
logical,parameter :: test_phase=.FALSE., vec_map=.FALSE., smag_sc=.FALSE.
logical,parameter :: check_dt=.TRUE.
integer,parameter :: stencil_pts=4
logical,parameter :: coarse_grain_flag=.FALSE.
real(kind=rprec),parameter::g=9.81_rprec, inv_strength=0.010_rprec,dTdz_top=0._rprec
! real(kind=rprec),parameter :: theta_top=300._rprec,T_scale=300._rprec,wt_s=0.1_rprec,T_init=283.15_rprec ! unstable
real(kind=rprec),parameter :: theta_top=293.15_rprec,T_scale=300._rprec,wt_s=0.05_rprec,T_init=293.15_rprec ! stable

! SKS
! Thinks these are not getting used anywhere.
real(kind=rprec),parameter::cap_thick=80._rprec, z_decay=1._rprec
real(kind=rprec),parameter::q_s1=11._rprec,q_s2=13._rprec,q_mix=12._rprec,q_top=12._rprec,wq_s=0.06
! SKS

end module param
