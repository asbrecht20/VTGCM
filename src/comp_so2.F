!
      module so2_module
      use params_module,only: nlevp1,nlonp4,nlat,dz
      use addfld_module,only: addfld
      implicit none
! -------------------------------------------------------------
! Boundary conditions, production and loss for SO2 are defined
! by comp_so2, and referenced by minor_so2. Comp_so2 is called
! from a latitude loop in dynamics. After comp_so2, dynamics calls
! minor_so2, which passes this module data to sub minor. Sub 
! minor contains 3d mpi calls and its own latitude loops.
!
!     real,dimension(nlonp4,nlat)   :: so2_ubc ! upper boundary
!     real,dimension(nlonp4,3,nlat) :: so2_lbc ! lower boundary
!     real,dimension(nlevp1,nlonp4,nlat) ::
!    |  so2_prod,     ! production of so2
!    |  so2_loss      ! loss of so2
!
! Boundary conditions and production and loss terms are allocated
! subdomains by sub alloc_so2 (called from allocdata.F).
! - last modified by SWB:  04/29/13
! 04/29/13 swb: comp_so2.F and minor_so2.F structure created
! 04/29/13 swb: P&L reactions and rates to be added
! 09/03/13 swb: phi_so2 same as for Ar
! 09/03/13 swb: mass mixing ratio for so2 at LBC added from KINETICS
! 01/27/14 swb: Adding new reactions to SO and SO2 sources/sinks
! 02/07/14 swb: Updating with new KINETICS 1-D profiles of 17 species VMRs
!               (nightside, Zhang et al, 2012 model simulation)
! 05/06/14 swb: Adding reaction specific production/loss diagnostic terms
!               for prints/plots  (2-D and 1-D)
! 05/27/14 swb: added s167 for source of SO2
! 06/16/14 swb: added s167 for source of SO2 (corrected)
! 06/16/14 swb: added s346 for source of SO2 (neglected)
! 06/18/14 swb: removed all ClCO3 reactions (zero out VMRs)
! 10/08/14 swb: s248 and s295 added; s249 dropped.
! 04/16/15 swb: detailed source and sink term diagnostics turned on
! -------------------------------------------------------------
!
      real,allocatable,dimension(:,:)   :: so2_ubc ! upper boundary (i,j)
      real,allocatable,dimension(:,:,:) :: so2_lbc ! lower boundary (i,3,j)
      real,allocatable,dimension(:,:,:) ::
     |  so2_prod,     ! production of so2 (k,i,j)
     |  so2_loss      ! loss of so2       (k,i,j)
!  Mutual diffusion coefficients (normalized) for O-SO2, CO-SO2, CO2-SO2
!  Use O-Ar, CO-Ar, CO2-Ar from mtgcm1.2/src: S. W. Bougher (9/3/2013)
!  Like Ar, set thermal diffusion coefficient = 0.0
      real :: phi_so2(3) =  (/2.398, 1.875, 0.6631/)
      real,parameter :: alfa_so2 = 0. ! thermal diffusion coefficient
!  Mass mixing ratio lower boundary condition for SO2 (at ~69 km)
!  massmix = mix*mso2/mbar = mix*64./44.
!  (1)  mix1 = 2.0E-07   ! Zhang et al (2011) Figure 8 eyeball value
!  (2)  mix2 = 3.0E-07   ! KINETICS reference simulation Zhang (2012)
!  (3)  mix3 = 1.372E-07 ! KINETICS reference simulation (Oct. 16, 2013 Tables)
!  (4)  mix4 = 2.777e-07 ! KINETICS reference simulation (Jan. 27, 2014 Tables)
!  (5)  mix5 = 2.5E-05   ! KINETICS reference simulation  Parkinson (2015) 
!     real :: psso2b =  2.91E-07   !  From Zhang et al (2011) Figure 8
!     real :: psso2b =  4.365E-07  !  From latest Zhang et al (2011) calculation tables
!     real :: psso2b =  2.000E-07  !  Updated from latest Zhang tables (2013)
!     real :: psso2b =  4.040E-07  !  Updated from latest Zhang tables (2014)
      real :: psso2b =  3.64E-05   !  KINETICS reference sim (Parkinson ea. 2015 @ 70 km)
      contains
!
!-----------------------------------------------------------------------
      subroutine alloc_so2(lon0,lon1,lat0,lat1)
!
! Allocate subdomains (without ghost cells) to module data for boundary 
! conditions and production and loss terms. This is called once per run 
! from sub allocdata (allocdata.F).
!
! Args:
      integer,intent(in) :: lon0,lon1,lat0,lat1
!
! Local:
      integer :: istat
!
! Allocate subdomains to boundary conditions:
      allocate(so2_ubc(lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_so2: error allocating',
     |  ' so2_ubc: stat=',i3)") istat
      allocate(so2_lbc(lon0:lon1,3,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_so2: error allocating',
     |  ' so2_lbc: stat=',i3)") istat
!
! Allocate subdomains to production and loss:
      allocate(so2_prod(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_so2: error allocating',
     |  ' so2_prod: stat=',i3)") istat
      allocate(so2_loss(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_so2: error allocating',
     |  ' so2_loss: stat=',i3)") istat

!     write(6,"('alloc_so2: allocated module data')")
!
      end subroutine alloc_so2
!
!-----------------------------------------------------------------------
      subroutine comp_so2(tn,o1,co,co2,n2,o2,so,so2,barm,
     |  lev0,lev1,lon0,lon1,lat)
!
!-----------------------------------------------------------------------
! Advance SO2 (minor): 
! Arguments for comp_so2 (14).
! ----------------------------------------------------------------------
! Advance so by one time step. This is called from driver at
! each subdomain latitude.
! Last modified: swb (05/29/13)
!              : add new rates from chemrates_sox
!              : finalize so2src, so2loss terms
!              : interp trace species from KINETICS to VTGCM grid
! ----------------------------------------------------------------------
!
      use qrj_module,only: qjso2a, qfjso2a, qjso2b, qfjso2b, qjso3,
     |                     mk_xnmbari
      use fields_module,only: tlbc
      use lbc,only: fb,b
      use chapman_module,only: so3vmr ! (nlevp1)
      use cons_module,only: rmassinv_o1,rmassinv_co,rmassinv_n2,p0,
     |  rmass_so, rmass_so2, rmassinv_so, rmassinv_so2,
     |  expz,expzmid_inv,boltz,rmassinv_co2,rmassinv_o2,rmass_o1,
     |  rmass_co,rmass_o2,rmass_so,rmass_so2,rmass_so3,gask,grav
      use chemrates_module,only:s114,s133,s166,s240,S244,s249,s258,
     | s261,s288,s293,s294,s300,s243,s286,s287,s167,s346,s248,s295
      use addfld_module,only: addfld
!
!     write(6,"('enter comp_so2: lat=',i2)") lat
! ----------------------------------------------------------------------
! Input args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) ::
     |  tn,      ! neutral temperature (deg K)
     |  o1,      ! atomic oxygen (mmr)
     |  co,      ! carbon monoxide (mmr)
     |  co2,     ! carbon dioxide (mmr)
     |  n2,      ! molecular nitrogen (mmr)
     |  o2,      ! molecular oxygen (mmr)
     |  so,      ! so (mmr)
     |  so2,     ! so2 (mmr)
     |  barm     ! mean molecular weight
! ----------------------------------------------------------------------
! Local:
      integer :: k,i,k0,k1,i0,i1,nk  ! for addfld (ask for check)
      real,dimension(lev0:lev1,lon0:lon1) :: 
     |  tni,            ! TN at interfaces
     |  qjso2ai,        ! Jso2a->SO+O at interfaces (qrj module): prod rates (#/cm3/sec)
     |  qfjso2ai,       ! Jso2a->SO+O at interfaces (qrj module): loss freq ((#//sec)
     |  qjso2bi,        ! Jso2b->S+O2 at interfaces (qrj module): prod rates (#/cm3/sec)
     |  qfjso2bi,       ! Jso2b->S+O2 at interfaces (qrj module): loss freq  (#/sec)
     |  qjso3i,         ! Jso3->SO2+O at interfaces (qrj module): prod rate  (#/cm3/sec)
     |  xnmbari,        ! xnmbar at interfaces (qrj module)
     |  xnmbar,         ! xnmbar at midlevels (qrj module)
     |  o1i,            ! O1 (MMR) interfaces
     |  coi,            ! CO (MMR) interfaces
     |  co2i,           ! CO2 (MMR) interfaces
     |  n2i,            ! N2 (MMR)  interfaces
     |  o2i,            ! O2 (MMR) interfaces
     |  soi,            ! SO (MMR) interfaces
     |  so2i,           ! SO2 (MMR) interfaces
     |  so3i,           ! SO3 (MMR) interfaces
     |  xtot,           ! Total number density (#/cm3)
     |  xo,             ! O1 number density (#/cm3)
     |  xco,            ! CO number density (#/cm3)
     |  xco2,           ! CO2 number density (#/cm3)
     |  xn2,            ! N2 number density (#/cm3)
     |  xo2,            ! O2 number density (#/cm3)
     |  xso,            ! SO number density (#/cm3)
     |  xso2,           ! SO2 number density (#/cm3)
     |  xso3,           ! SO3 number density (#/cm3)
     |  so2prod1,       ! SO2 production (#/cm3/s)
     |  so2prod2,       ! SO2 production (#/cm3/s)
     |  so2ps240,       ! SO2 production (#/cm3/s)
     |  so2ps133,       ! SO2 production (#/cm3/s)
     |  so2ps167,       ! SO2 production (#/cm3/s)
     |  so2ps244,       ! SO2 production (#/cm3/s)
     |  so2ps248,       ! SO2 production (#/cm3/s)
     |  so2ps293,       ! SO2 production (#/cm3/s)
     |  so2ps294,       ! SO2 production (#/cm3/s)
     |  so2ps295,       ! SO2 production (#/cm3/s)
     |  so2ps300,       ! SO2 production (#/cm3/s)
     |  so2ps286,       ! SO2 production (#/cm3/s)
     |  so2ps346,       ! SO2 production (#/cm3/s)
     |  so2loss1,       ! SO2 loss1 (#/cm3/s)
     |  so2loss2,       ! SO2 loss2 (#/cm3/s)
     |  so2ls258,       ! SO2 loss2 (#/cm3/s)
     |  so2ls114,       ! SO2 loss2 (#/cm3/s)
     |  so2ls261,       ! SO2 loss2 (#/cm3/s)
     |  xno2,           ! NO2 (#/cm3) 
     |  xcl,            ! CL (#/cm3) 
     |  xcl2,           ! CL2 (#/cm3) 
     |  xclo,           ! CLO (#/cm3) 
     |  xclco,          ! CLCO (#/cm3) 
     |  xclco3,         ! CLCO3 (#/cm3) 
     |  xclso2,         ! CLSO2 (#/cm3) 
     |  xoscl,          ! OSCl (#/cm3) 
     |  xs,             ! S (#/cm3) 
     |  xh,             ! H (#/cm3) 
     |  xh2,            ! H2 (#/cm3)
     |  xoh,            ! OH (#/cm3) 
     |  xho2,           ! HO2 (#.cm3) 
     |  xo3,            ! O3 (#/cm3) 
     |  xhcl            ! hcl (#/cm3)
!
! -----------------------------------------------------------------------
! Specificed VMR specifically for SOx, CLOx, HOx chemistry 
! From recent 1-D KINETICS calculation for nightside
! 29-MAY-2013 (private communication, C. Parkinson)
! From recent 1-D KINETICS calculation for nightside
! 28-JAN-2014 (private communication, C. Parkinson)
      real,dimension(lev0:lev1) ::
     |  no2,           ! NO2 (vmr) 
     |  cl,            ! CL (vmr) 
     |  cl2,           ! CL2 (vmr) 
     |  clo,           ! CLO (vmr) 
     |  clco,          ! CLCO (vmr) 
     |  clco3,         ! CLCO3 (vmr) 
     |  clso2,         ! CLSO2 (vmr) 
     |  oscl,          ! OSCl (vmr) 
     |  s,             ! S (vmr) 
     |  h,             ! H (vmr) 
     |  h2,            ! H2 (vmr)
     |  oh,            ! OH (vmr) 
     |  ho2,           ! HO2 (vmr) 
     |  o3,            ! O3 (vmr) 
     |  hcl            ! hcl (vmr)
!
! ----------------------------------------------------------------------
! For addfld calls: Check what is needed!
      k0=lev0
      k1=lev1
      i0=lon0
      i1=lon1
      nk = lev1-1
!
! *******************************************************************
! Low Altitude Chemistry from Chris Parkinson Venus KINETICS code
! Profiles in VMR: 22 (values) + 47 (zeros) = 69 levels
!  (07-FEB-2014 Update to VMR for 17 species)
! Profiles are identical to Zhang et al (2012) Final profiles!
! Alt range = ~69 - ~112km (zp = -16.0 to -5.5)
! Revised for latest 1-D KINETICS Model Simulations
! -- Set top 46 levels with level = 23 VMR value to avoid discontinuity
! -- recast on VTGCM zp grid
! -- SO3 = fixed (always)
! -- SO2 = from history (fixed to start)
! -- SO  = from history (fixed to start)
! -- rest = fixed always
! -- ClCO3 VMR zeroed out to remove P/L terms:  s244, s261
!********************************************************************
!** ZP Gridded structure needed for VTGCM **
!    |    -16.00,-15.50,-15.00,-14.50,-14.00,
!    |    -13.50,-13.00,-12.50,-12.00,-11.50,
!    |    -11.00,-10.50,-10.00,-9.50,-9.00,
!    |    -8.50,-8.00,-7.50,-7.00,-6.50,
!    |    -6.00,-5.50,-5.00,
!    |    -4.50,-4.00,-3.50,-3.00,-2.50,
!    |    -2.00,-1.50,-1.00,-0.50,0.00 ,
!    |    0.50 ,1.00 ,1.50 ,2.00 ,2.50 ,
!    |    3.00 ,3.50 ,4.00 ,4.50 ,5.00 ,
!    |    5.50 ,6.00 ,6.50 ,7.00 ,7.50 ,
!    |    8.00 ,8.50 ,9.00 ,9.50 ,10.00,
!    |    10.50,11.00,11.50,12.00,12.50,
!    |    13.00,13.50,14.00,14.50,15.00,
!    |    15.50,16.00,16.50,17.00,17.50,
!    |    18.00
!********************************************************************
      no2(:) = (/ 
     |    5.93e-011,8.70e-011,8.50e-011,7.64e-011,
     |    6.78e-011,8.48e-011,1.37e-010,1.35e-010,
     |    5.83e-011,8.72e-012,1.74e-012,5.07e-013,
     |    9.81e-014,4.34e-014,1.96e-014,1.08e-014,
     |    5.72e-015,3.29e-015,1.84e-015,9.91e-016,
     |    4.97e-016,2.32e-016,0.00e+000,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 /)
      cl(:) = (/ 
     |    1.62e-009,5.21e-009,9.68e-009,1.36e-008,
     |    1.41e-008,1.29e-008,1.35e-008,1.39e-008,
     |    1.73e-008,2.34e-008,3.06e-008,3.37e-008,
     |    4.20e-008,6.37e-008,9.28e-008,1.22e-007,
     |    1.49e-007,1.71e-007,1.91e-007,2.08e-007,
     |    2.24e-007,2.40e-007,0.00e+000,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 /)
      cl2(:) = (/ 
     |    2.18e-009,1.74e-009,1.36e-009,1.06e-009,
     |    7.29e-010,3.29e-010,1.71e-010,8.65e-011,
     |    6.09e-011,5.27e-011,3.66e-011,1.45e-011,
     |    8.07e-012,5.36e-012,3.53e-012,2.22e-012,
     |    1.22e-012,5.88e-013,2.77e-013,1.23e-013,
     |    5.23e-014,2.18e-014,0.00e+000,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 /)
      clco(:) = (/ 
     |    4.01e-012,1.19e-011,1.29e-011,2.25e-011,
     |    1.02e-010,1.97e-010,2.71e-010,3.53e-010,
     |    4.05e-010,3.87e-010,2.01e-010,6.45e-011,
     |    3.11e-011,1.41e-011,5.66e-012,2.20e-012,
     |    8.79e-013,3.51e-013,1.56e-013,7.38e-014,
     |    3.70e-014,1.94e-014,0.00e+000,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 /)
      clo(:) = (/ 
     |    6.21e-012,1.27e-011,1.54e-011,1.81e-011,
     |    3.04e-011,8.92e-011,2.93e-010,6.01e-010,
     |    5.87e-010,2.22e-010,5.54e-011,2.86e-011,
     |    1.80e-011,1.19e-011,7.65e-012,4.65e-012,
     |    2.58e-012,1.30e-012,6.25e-013,2.84e-013,
     |    1.26e-013,5.61e-014,0.00e+000,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 /)
      clco3(:) = (/ 
     |    1.16e-012,7.26e-012,1.12e-011,1.77e-011,
     |    5.45e-011,2.38e-010,8.43e-010,1.87e-009,
     |    2.28e-009,1.43e-009,2.11e-010,2.22e-011,
     |    3.97e-012,7.03e-013,7.92e-014,7.79e-015,
     |    9.43e-016,1.09e-016,1.64e-017,2.78e-018,
     |    5.35e-019,1.13e-019,0.00e+000,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 /)
      h(:) = (/ 
     |    2.40e-013,7.68e-013,1.61e-012,3.55e-012,
     |    1.04e-011,1.44e-011,1.13e-011,1.21e-011,
     |    2.37e-011,1.09e-010,7.44e-010,4.57e-009,
     |    2.37e-008,7.32e-008,1.51e-007,2.47e-007,
     |    3.52e-007,4.60e-007,5.66e-007,6.63e-007,
     |    7.47e-007,8.13e-007,0.00e+000,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 /)
      h2(:) = (/ 
     |    2.58e-010,2.54e-010,2.92e-010,5.83e-010,
     |    1.81e-009,3.81e-009,6.67e-009,1.07e-008,
     |    1.60e-008,2.24e-008,3.04e-008,4.01e-008,
     |    5.16e-008,6.17e-008,6.90e-008,7.36e-008,
     |    7.66e-008,7.88e-008,8.06e-008,8.24e-008,
     |    8.45e-008,8.70e-008,0.00e+000,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 /)
      hcl(:) = (/ 
     |    3.95e-007,3.92e-007,3.89e-007,3.86e-007,
     |    3.82e-007,3.78e-007,3.75e-007,3.73e-007,
     |    3.72e-007,3.70e-007,3.69e-007,3.68e-007,
     |    3.60e-007,3.38e-007,3.09e-007,2.80e-007,
     |    2.54e-007,2.31e-007,2.12e-007,1.95e-007,
     |    1.79e-007,1.64e-007,0.00e+000,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 /)
      ho2(:) = (/ 
     |    7.11e-015,3.56e-014,9.52e-014,1.59e-013,
     |    2.47e-013,6.23e-013,1.27e-012,2.27e-012,
     |    3.36e-012,4.47e-012,6.77e-012,1.42e-011,
     |    2.57e-011,2.36e-011,1.48e-011,7.48e-012,
     |    3.45e-012,1.52e-012,6.82e-013,3.09e-013,
     |    1.42e-013,6.63e-014,0.00e+000,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 /)
      o3(:) = (/ 
     |    6.01e-013,1.57e-012,4.66e-012,8.07e-012,
     |    7.31e-012,2.16e-011,6.30e-011,1.43e-010,
     |    3.98e-010,1.47e-009,4.94e-009,9.82e-009,
     |    1.21e-008,1.43e-008,1.67e-008,1.90e-008,
     |    1.85e-008,1.55e-008,1.15e-008,7.35e-009,
     |    4.19e-009,2.20e-009,0.00e+000,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 /)
      oh(:) = (/ 
     |    3.94e-016,1.94e-015,6.04e-015,1.31e-014,
     |    2.64e-014,7.05e-014,1.39e-013,2.43e-013,
     |    5.17e-013,1.67e-012,6.32e-012,2.15e-011,
     |    4.94e-011,5.60e-011,4.25e-011,2.63e-011,
     |    1.46e-011,7.65e-012,3.91e-012,1.98e-012,
     |    1.06e-012,6.90e-013,0.00e+000,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 /)
      s(:) = (/ 
     |    4.61e-012,2.34e-012,1.46e-012,1.64e-012,
     |    4.10e-012,3.16e-012,1.27e-012,7.98e-013,
     |    1.27e-012,2.94e-012,5.34e-012,7.19e-012,
     |    1.28e-011,2.25e-011,3.56e-011,5.52e-011,
     |    9.17e-011,1.63e-010,3.01e-010,5.63e-010,
     |    1.03e-009,1.81e-009,0.00e+000,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 /)
!     so(:) = (/ 
!    |    3.03e-009,2.52e-009,2.36e-009,2.43e-009,
!    |    2.90e-009,3.56e-009,3.69e-009,4.47e-009,
!    |    1.02e-008,2.58e-008,4.57e-008,5.40e-008,
!    |    7.98e-008,1.09e-007,1.32e-007,1.51e-007,
!    |    1.77e-007,2.14e-007,2.59e-007,3.08e-007,
!    |    3.50e-007,3.78e-007,0.00e+000,
!    |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
!    |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
!    |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
!    |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
!    |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
!    |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
!    |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
!    |    0.00 ,0.00 ,0.00 ,0.00 /)
!     so2(:) = (/ 
!    |    2.78e-007,1.21e-007,6.52e-008,4.18e-008,
!    |    3.03e-008,3.11e-008,4.06e-008,5.51e-008,
!    |    7.13e-008,8.43e-008,1.01e-007,1.40e-007,
!    |    1.78e-007,2.27e-007,2.61e-007,2.70e-007,
!    |    2.56e-007,2.23e-007,1.79e-007,1.30e-007,
!    |    8.61e-008,5.71e-008,0.00e+000,
!    |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
!    |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
!    |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
!    |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
!    |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
!    |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
!    |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
!    |    0.00 ,0.00 ,0.00 ,0.00 /)
!     so3(:) = (/ 
!    |    1.58e-009,7.27e-010,6.87e-010,3.91e-010,
!    |    4.92e-011,2.00e-011,2.48e-011,3.12e-011,
!    |    2.44e-011,7.19e-012,8.55e-012,5.88e-011,
!    |    6.81e-010,1.03e-008,1.27e-008,9.23e-009,
!    |    5.24e-009,2.64e-009,1.24e-009,5.98e-010,
!    |    3.21e-010,2.12e-010,0.00e+000,
!    |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
!    |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
!    |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
!    |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
!    |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
!    |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
!    |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
!    |    0.00 ,0.00 ,0.00 ,0.00 /)
      clso2(:) = (/ 
     |    5.11e-010,3.80e-010,1.41e-010,4.44e-011,
     |    2.35e-011,1.34e-011,1.08e-011,9.04e-012,
     |    5.55e-012,1.91e-012,3.35e-013,8.33e-014,
     |    3.57e-014,1.50e-014,5.39e-015,1.73e-015,
     |    5.32e-016,1.48e-016,4.30e-017,1.20e-017,
     |    3.26e-018,8.98e-019,0.00e+000,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 /)
      oscl(:) = (/ 
     |    6.72e-011,1.51e-011,3.12e-012,8.98e-013,
     |    7.02e-013,6.33e-013,4.77e-013,4.20e-013,
     |    6.16e-013,9.94e-013,7.40e-013,2.19e-013,
     |    1.21e-013,5.66e-014,2.11e-014,7.22e-015,
     |    2.62e-015,9.52e-016,4.06e-016,1.81e-016,
     |    8.28e-017,3.74e-017,0.00e+000,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,0.00 ,
     |    0.00 ,0.00 ,0.00 ,0.00 /)
!
! ----------------------------------------------------------------------
! Lower boundary:
! so2_lbc(:,1)=A, so2_lbc(:,2)=B, so2_lbc(:,3)=C define lower boundary 
!   condition where A*DPSX/DZ + B*PSX + C = 0.
!
      do i=lon0,lon1
!
! Value at bottom given by specified mass missing ratio (start)
        so2_lbc(i,1,lat) = 0.
        so2_lbc(i,2,lat) = 1.
        so2_lbc(i,3,lat) = -psso2b
!
! Upper boundary: Zero diffusive flux at top
!       so2_ubc(i,lat) = 0.
      enddo ! i=lon0,lon1
! --------------------------------------------------------------------
! Calculate tn at interfaces:
      do i=lon0,lon1
      tni(lev0,i) = tlbc(i,lat)
      	do k=lev0+1,lev1-1
         	tni(k,i) = .5*(tn(k-1,i)+tn(k,i))
      	enddo
      tni(lev1,i) = tn(lev1-1,i) ! nlevp1 <- nlev
      enddo
!
! ---------------------------------------------------------------
!
! o1,co,n2,o2,so,so2 at interface levels: mass mixing ratios
      do i=lon0,lon1
        do k=lev0,lev1-1
          o1i (k+1,i) = 0.5*(o1(k,i)+o1(k+1,i))
          coi (k+1,i) = 0.5*(co(k,i)+co(k+1,i))
          n2i (k+1,i) = 0.5*(n2(k,i)+n2(k+1,i))
          o2i (k+1,i) = 0.5*(o2(k,i)+o2(k+1,i))
          co2i (k+1,i) = 0.5*(co2(k,i)+co2(k+1,i))
          soi (k+1,i)=0.5*(so(k,i)+so(k+1,i))
          so2i(k+1,i)=0.5*(so2(k,i)+so2(k+1,i))
        enddo
      enddo
! Bottom boundary: MMR
      do i=lon0,lon1
        o1i(lev0,i)= .5*((b(i,1,1)+1.)*o1(1,i)+b(i,1,2)*co(1,i)+fb(i,1))
        coi(lev0,i)= .5*(b(i,2,1)*o1(1,i)+(b(i,2,2)+1.)*co(1,i)+fb(i,2))
        n2i(lev0,i)=  n2i(lev0+1,i)
        o2i(lev0,i)=  o2i(lev0+1,i)
        co2i(lev0,i)= co2i(lev0+1,i)
        soi(lev0,i)=  soi(lev0+1,i)
        so2i(lev0,i)= so2i(lev0+1,i)
      enddo
! --------------------------------------------------------------------
!
! Calculate p0*e(-z)*barm/kT (1/(N*MBAR))at interfaces
!     subroutine mk_xnmbari(tni,barm,xnmbari,lev0,lev1,lon0,lon1,lat)
!
      call mk_xnmbari(tni,barm(:,lon0:lon1),
     |  xnmbari,
     |  lev0,lev1,lon0,lon1,lat)
! ----------------------------------------------------------------------
! Species Number densities all at interfaces (#/cm3)
! Given Qrj dissociation rates (#/cm3/sec) at midlevels from qrj
! Calculate Qrj dissociation rates (#/sec) at interfaces: Jso2a, Jso2b
! Calculate Qrj dissociation rates (#/cm3/sec) at interfaces: Jso2a, Jso2b,Jso3
! Change k = lev0,lev-1 to k = lev0,lev1 (interfaces all fields)
!
      do i=lon0,lon1
        do k=lev0,lev1
          xo(k,i)   = xnmbari(k,i)*o1i(k,i)*rmassinv_o1
	  if (xo(k,i) .lt. 0.0) xo(k,i) = 0.0
          xco(k,i)  = xnmbari(k,i)*coi(k,i)*rmassinv_co
	  if (xco(k,i) .lt. 0.0) xco(k,i) = 0.0
          xn2(k,i)  = xnmbari(k,i)*n2i(k,i)*rmassinv_n2
	  if (xn2(k,i) .lt. 0.0) xn2(k,i) = 0.0
          xo2(k,i) = xnmbari(k,i)*o2i(k,i)*rmassinv_o2
	  if (xo2(k,i) .lt. 0.0) xo2(k,i) = 0.0
          xco2(k,i) = xnmbari(k,i)*co2i(k,i)*rmassinv_co2
	  if (xco2(k,i) .lt. 0.0) xco2(k,i) = 0.0
          xso(k,i)  = xnmbari(k,i)*soi(k,i)*rmassinv_so
	  if (xso(k,i) .lt. 0.0) xso(k,i) = 0.0
          xso2(k,i)  = xnmbari(k,i)*so2i(k,i)*rmassinv_so2
	  if (xso2(k,i) .lt. 0.0) xso2(k,i) = 0.0
          xso3(k,i)  = so3vmr(k)*xnmbari(k,i)/barm(k,i)
	  if (xso3(k,i) .lt. 0.0) xso3(k,i) = 0.0
          xtot(k,i) = xo(k,i)+xco(k,i)+xco2(k,i)+xn2(k,i)
!  2-Photolysis frequencies (#/sec) and 3-production/loss rates (#/cm3/s)
!  Already on interfaces from qrj.F
          qfjso2ai(k,i) = qfjso2a(k,i,lat)
	  if (qfjso2ai(k,i) .lt. 0.0) qfjso2ai(k,i) = 0.0
          qjso2ai(k,i) = qjso2a(k,i,lat)
	  if (qjso2ai(k,i) .lt. 0.0) qjso2ai(k,i) = 0.0
          qfjso2bi(k,i) = qfjso2b(k,i,lat)
	  if (qfjso2bi(k,i) .lt. 0.0) qfjso2bi(k,i) = 0.0
          qjso2bi(k,i) = qjso2b(k,i,lat)
	  if (qjso2bi(k,i) .lt. 0.0) qjso2bi(k,i) = 0.0
          qjso3i(k,i) = qjso3(k,i,lat)
	  if (qjso3i(k,i) .lt. 0.0) qjso3i(k,i) = 0.0
!
!        write(6,"('comp_so2: k,i,lat=',3i4,' so2(k,i)=',e12.4,
!    |      ' xnmbari(k,i)=',e12.4,' qjso2a(k)(k+1)=',2e12.4,
!    |      ' xso2(k,i)=',e12.4,' qjso2ai(k,i)=',e12.4)") k,i,lat,
!    |      so2(k,i),qjso2a(k,i,lat),qjso2a(k+1,i,lat),xso2(k,i),
!    |      qjso2ai(k,i) 
!
!   LOG10 DENSITIES FOR PLOTTING 
!         xco2l(k,i)=log10(xco2(k,i))
!         xcol(k,i)=log10(xco(k,i))
!         xn2l(k,i)=log10(xo2(k,i))
!         xol(k,i)=log10(xo(k,i))
!
!
!  Converting vmr to number density
!         xno2(k,i)   = xtot(k,i)*no2(k)
          xcl(k,i)   = xtot(k,i)*cl(k)
          xcl2(k,i)  = xtot(k,i)*cl2(k)
          xclo(k,i)  = xtot(k,i)*clo(k)
          xclco(k,i) = xtot(k,i)*clco(k)
          xclco3(k,i) = xtot(k,i)*clco3(k)
          xclso2(k,i) = xtot(k,i)*clso2(k)
          xoscl(k,i) = xtot(k,i)*oscl(k)
          xs(k,i) = xtot(k,i)*s(k)
          xh(k,i)   = xtot(k,i)*h(k)
          xh2(k,i)  = xtot(k,i)*h2(k)
          xoh(k,i)   = xtot(k,i)*oh(k)
          xho2(k,i)  = xtot(k,i)*ho2(k)
          xo3(k,i) = xtot(k,i)*o3(k)
          xhcl(k,i) = xtot(k,i)*hcl(k)
!
!        write(6,"('comp_so2: k,i,lat=',3i4,'xno2(k,i)=',e12.4)")
!    |      k,i,lat,xno2(k,i)
!
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!   Test 1: OK (all positive)
!     write(6,"('comp_so2 before calculate prod,loss: xso=',2e12.4,
!    |  ' xso2=',2e12.4,' xso3=',2e12.4)")
!    |  minval(xso),maxval(xso),
!    |  minval(xso2),maxval(xso2),
!    |  minval(xso3),maxval(xso3)
!
!   Test 2: OK, all positive 
!     write(6,"('comp_so2 before calculate prod,loss: qjso3i=',2e12.4,
!    |  ' qjso2ai=',2e12.4,' qjso2bi=',2e12.4)")
!    |  minval(qjso3i),maxval(qjso3i),
!    |  minval(qjso2ai),maxval(qjso2ai),
!    |  minval(qjso2bi),maxval(qjso2bi)
!
!     write(6,"('comp_so2 after calculate qjso2bi: lat=',i4,
!    |  ' qjso2bi mnmx=',2e12.4)") lat,minval(qjso2bi),maxval(qjso2bi)
!
! ----------------------------------------------------------------------
! (1)  Major Sources and Losses for SO2 (before NOX Chemistry Added)
!      Based upon Zhang et al (2012)  Figure 12c-d.
! ---------------------------------------------------------------------
      do i=lon0,lon1
        do k=lev0,lev1-1
! ---------------------------------------------------------------------
!     SOURCES (#/cm3.sec) AND SINKS(#/sec): Prognostics
!     SOURCES (#/cm3.sec) AND SINKS(#/cm3/sec): Diagnostics
! ---------------------------------------------------------------------
!   Add SO3 photolysis to the production of SO2! See KINETICS!
!
      so2_prod(k,i,lat) = qjso3i(k,i)               + 
     |         xo(k,i)*xso(k,i)*xco2(k,i)*s240      +
     |         xso(k,i)*xso(k,i)*s248(k,i,lat)      +
     |         xclo(k,i)*xso(k,i)*s133              +
     |         xclco3(k,i)*xso(k,i)*s244            + 
     |         xo(k,i)*xclso2(k,i)*s293             +
     |         xh(k,i)*xclso2(k,i)*s294             +
     |         xcl(k,i)*xclso2(k,i)*s295            +
     |         2.*xclso2(k,i)*xclso2(k,i)*s300      +
     |         xo(k,i)*xoscl(k,i)*s286(k,i,lat)     +      
     |         xo2(k,i)*xso(k,i)*s167(k,i,lat)      +
     |         xso(k,i)*xno2(k,i)*s346
!
! Subscript #2 of the array SO2PROD has value 39 which is greater than the upper bound of 38
          so2prod1(k,i) = qjso3i(k,i)
          so2ps240(k,i) = xo(k,i)*xso(k,i)*xco2(k,i)*s240
          so2ps248(k,i) = xso(k,i)*xso(k,i)*s248(k,i,lat)
          so2ps133(k,i) = xclo(k,i)*xso(k,i)*s133
          so2ps167(k,i) = xo2(k,i)*xso(k,i)*s167(k,i,lat)
          so2ps244(k,i) = xclco3(k,i)*xso(k,i)*s244
          so2ps293(k,i) = xo(k,i)*xclso2(k,i)*s293
          so2ps294(k,i) = xh(k,i)*xclso2(k,i)*s294
          so2ps295(k,i) = xcl(k,i)*xclso2(k,i)*s295
          so2ps300(k,i) = 2.*xclso2(k,i)*xclso2(k,i)*s300
          so2ps286(k,i) = xo(k,i)*xoscl(k,i)*s286(k,i,lat)
          so2ps346(k,i) = xso(k,i)*xno2(k,i)*s346
!
!        write(6,"('comp_so2: k,i,lat=',3i5,'so2ps346(k,i)=',e12.4)")
!    |      k,i,lat,so2ps346(k,i)
!
!
!         so2prod2(k,i) = xo(k,i)*xso(k,i)*xco2(k,i)*s240      +
!    |         xclo(k,i)*xso(k,i)*s133                         +
!    |         xclco3(k,i)*xso(k,i)*s244                       + 
!    |         xo(k,i)*xclso2(k,i)*s293                        +
!    |         xh(k,i)*xclso2(k,i)*s294                        +
!    |         2.*xclso2(k,i)*s300                             +
!    |         xo(k,i)*xoscl(k,i)*s286(k,i,lat)      
! --------------------------------------
!
          so2_loss(k,i,lat)= -(qfjso2ai(k,i) + qfjso2bi(k,i) +
     |        xo(k,i)*xco2(k,i)*s258(k,i,lat) + 
     |        xcl(k,i)*xco2(k,i)*s114(k,i,lat) +
     |        xclco3(k,i)*s261)
!         so2_loss(k,i,lat)= -(qfjso2ai(k,i) + qfjso2bi(k,i))
!
! Subscript #2 of the array SO2LOSS1 has value 1 which is less than the lower bound of 39
          so2loss1(k,i) = (qjso2ai(k,i) + qjso2bi(k,i))
          so2ls258(k,i) = xo(k,i)*xso2(k,i)*xco2(k,i)*s258(k,i,lat)
          so2ls114(k,i) = xcl(k,i)*xso2(k,i)*xco2(k,i)*s114(k,i,lat)
          so2ls261(k,i) = xclco3(k,i)*xso2(k,i)*s261
!
!         so2loss2(k,i) = xo(k,i)*xso2(k,i)*xco2(k,i)*s258(k,i,lat) +
!    |         xcl(k,i)*xso2(k,i)*xco2(k,i)*s114(k,i,lat)           +
!    |         xclco3(k,i)*xso2(k,i)*s261
!         so2loss2(k,i) = 0.0
!----------------------------------------------------------------------------------------
!  Convert to Log10 units for plotting diagnostic terms : Initial
        if (so2prod1(k,i) .GT. 1.0) THEN
	    so2prod1(k,i)= alog10(so2prod1(k,i))
	 else
	    so2prod1(k,i) = 0.
	 endif
!       if (so2prod2(k,i) .GT. 1.0) THEN
!          so2prod2(k,i)= alog10(so2prod2(k,i))
!      else
!          so2prod2(k,i) = 0.
!      endif
        if (so2loss1(k,i) .GT. 1.0) THEN
	    so2loss1(k,i)= alog10(so2loss1(k,i))
	 else
	    so2loss1(k,i) = 0.
	 endif
!       if (so2loss2(k,i) .GT. 1.0) THEN
!          so2loss2(k,i)= alog10(so2loss2(k,i))
!       else
!          so2loss2(k,i) = 0.
!       endif
!  Convert to Log10 units for plotting diagnostic terms : Expanded
!  All SO2 production terms
         if (so2ps240(k,i) .GT. 1.0) THEN
	    so2ps240(k,i)= alog10(so2ps240(k,i))
	 else
	    so2ps240(k,i) = 0.
	 endif
         if (so2ps133(k,i) .GT. 1.0) THEN
	    so2ps133(k,i)= alog10(so2ps133(k,i))
	 else
	    so2ps133(k,i) = 0.
	 endif
         if (so2ps167(k,i) .GT. 1.0) THEN
	    so2ps167(k,i)= alog10(so2ps167(k,i))
	 else
	    so2ps167(k,i) = 0.
	 endif
         if (so2ps244(k,i) .GT. 1.0) THEN
	    so2ps244(k,i)= alog10(so2ps244(k,i))
	 else
	    so2ps244(k,i) = 0.
	 endif
         if (so2ps248(k,i) .GT. 1.0) THEN
	    so2ps248(k,i)= alog10(so2ps248(k,i))
	 else
	    so2ps248(k,i) = 0.
	 endif
         if (so2ps293(k,i) .GT. 1.0) THEN
	    so2ps293(k,i)= alog10(so2ps293(k,i))
	 else
	    so2ps293(k,i) = 0.
	 endif
         if (so2ps294(k,i) .GT. 1.0) THEN
	    so2ps294(k,i)= alog10(so2ps294(k,i))
	 else
	    so2ps294(k,i) = 0.
	 endif
         if (so2ps295(k,i) .GT. 1.0) THEN
             so2ps295(k,i)= alog10(so2ps295(k,i))
         else
             so2ps295(k,i) = 0.
 	 endif
         if (so2ps300(k,i) .GT. 1.0) THEN
	    so2ps300(k,i)= alog10(so2ps300(k,i))
	 else
	    so2ps300(k,i) = 0.
	 endif
         if (so2ps286(k,i) .GT. 1.0) THEN
	    so2ps286(k,i)= alog10(so2ps286(k,i))
	 else
	    so2ps286(k,i) = 0.
	 endif
         if (so2ps346(k,i) .GT. 1.0) THEN
            so2ps346(k,i)= alog10(so2ps346(k,i))
         else
            so2ps346(k,i) = 0.
         endif
!  All SO2 loss terms
         if (so2ls258(k,i) .GT. 1.0) THEN
	    so2ls258(k,i)= alog10(so2ls258(k,i))
	 else
	    so2ls258(k,i) = 0.
	 endif
         if (so2ls114(k,i) .GT. 1.0) THEN
	    so2ls114(k,i)= alog10(so2ls114(k,i))
	 else
	    so2ls114(k,i) = 0.
	 endif
         if (so2ls261(k,i) .GT. 1.0) THEN
	    so2ls261(k,i)= alog10(so2ls261(k,i))
	 else
	    so2ls261(k,i) = 0.
	 endif
!----------------------------------------------------------------------------------------
!
!  test only: Make sure model production (+) and loss (-) terms are correct sign
!
!     if (so2_prod(k,i,lat) .LT. 0.0) so2_prod(k,i,lat) = 0.0
!     if (so2_loss(k,i,lat) .GT. 0.0) so2_loss(k,i,lat) = 0.0
!----------------------------------------------------------------------------------------
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!  Test 3:
!     write(6,"('comp_so2 after calculate prod,loss: 
!    |   so2prod1=',2e12.4, ' so2prod2=',2e12.4,
!    |  ' so2loss1=',2e12.4,' so2loss2=',2e12.4)")
!    |  minval(so2prod1),maxval(so2prod1),
!    |  minval(so2prod2),maxval(so2prod2),
!    |  minval(so2loss1),maxval(so2loss1),
!    |  minval(so2loss2),maxval(so2loss2)
!  Test 4:
!     write(6,"('comp_so2 after calculate prod,loss2: 
!    |   so2_prod=',2e12.4, ' so2_loss=',2e12.4)")
!    |  minval(so2_prod),maxval(so2_prod),
!    |  minval(so2_loss),maxval(so2_loss)
!
! ----------------------------------------------------------------------
!   * SO2 Production and Losses
      call addfld('SO2PROD1',' ',' ',so2prod1,'lev',k0,k1,'lon',i0,i1,
     |  lat)
!     call addfld('SO2PROD2',' ',' ',so2prod2,'lev',k0,k1,'lon',i0,i1,
!    |  lat)
      call addfld('SO2PS240',' ',' ',so2ps240,'lev',k0,k1,'lon',i0,i1,
     |  lat)
      call addfld('SO2PS133',' ',' ',so2ps133,'lev',k0,k1,'lon',i0,i1,
     |  lat)
      call addfld('SO2PS167',' ',' ',so2ps167,'lev',k0,k1,'lon',i0,i1,
     |  lat)
      call addfld('SO2PS244',' ',' ',so2ps244,'lev',k0,k1,'lon',i0,i1,
     |  lat)
      call addfld('SO2PS248',' ',' ',so2ps248,'lev',k0,k1,'lon',i0,i1,
     |  lat)
      call addfld('SO2PS293',' ',' ',so2ps293,'lev',k0,k1,'lon',i0,i1,
     |  lat)
      call addfld('SO2PS294',' ',' ',so2ps294,'lev',k0,k1,'lon',i0,i1,
     |  lat)
      call addfld('SO2PS295',' ',' ',so2ps295,'lev',k0,k1,'lon',i0,i1,
     |  lat)
      call addfld('SO2PS300',' ',' ',so2ps300,'lev',k0,k1,'lon',i0,i1,
     |  lat)
      call addfld('SO2PS286',' ',' ',so2ps286,'lev',k0,k1,'lon',i0,i1,
     |  lat)
      call addfld('SO2PS346',' ',' ',so2ps346,'lev',k0,k1,'lon',i0,i1,
     |  lat)
      call addfld('SO2LOSS1',' ',' ',so2loss1,'lev',k0,k1,'lon',i0,i1,
     |  lat)
!     call addfld('SO2LOSS2',' ',' ',so2loss2,'lev',k0,k1,'lon',i0,i1,
!    |  lat)
      call addfld('SO2LS258',' ',' ',so2ls258,'lev',k0,k1,'lon',i0,i1,
     |  lat)
      call addfld('SO2LS114',' ',' ',so2ls114,'lev',k0,k1,'lon',i0,i1,
     |  lat)
      call addfld('SO2LS261',' ',' ',so2ls261,'lev',k0,k1,'lon',i0,i1,
     |  lat)
      call addfld('QFJSO2AI',' ',' ',qfjso2ai,'lev',k0,k1,'lon',i0,i1,
     |  lat)
      call addfld('QFJSO2BI',' ',' ',qfjso2bi,'lev',k0,k1,'lon',i0,i1,
     |  lat)
!
!     write(6,"(' comp_so2 returning: lat=',i4)") lat

      end subroutine comp_so2
!-----------------------------------------------------------------------
      subroutine minor_so2(tn,difk,o1,co,co2,n2,so2,so2_nm,so2_out,
     |  so2_nm_out,lev0,lev1,lon0,lon1,lat0,lat1)

!     Arguments (16): o1, co, co2, n2 needed by minor.F for barm, etc.
      use cons_module,only: rmass_so, rmass_so2
!
! -----------------------------------------------------------------------
! Input args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(in) ::
     |  tn,      ! neutral temperature (deg K)
     |  difk,    ! eddy diffusion
     |  o1,      ! atomic oxygen (mmr)
     |  co,      ! carbon monoxide (mmr)
     |  co2,     ! carbon dioxide (mmr)
     |  n2,      ! molecular nitrogen (mmr)
     |  so2,     ! so2 (mmr)
     |  so2_nm   ! so2 at time n-1 
!
! -----------------------------------------------------------------------
! Output args:
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(out) ::
     |  so2_out,    ! SO2 output 
     |  so2_nm_out  ! SO2 output at time n-1
!
! -----------------------------------------------------------------------
! Local:
      integer :: lat
      integer :: i0,i1,nk ! for addfld
!
!     write(6,"('enter minor_so')")
      i0 = lon0
      i1 = lon1
      nk = lev1-1
! -----------------------------------------------------------------------
! Minor returns so_out and so_nm_out. Module data so_prod,
! so_loss, etc, were defined by comp_so.  
!
!     subroutine minor(tn,difk,o1,co,co2,n2,fcomp,fcomp_tm1,fcomp_out,
!    |  fcomp_tm1_out,sloss,sprod,flbc,fubc,rmx,phix,alfax,
!    |  lev0,lev1,lon0,lon1,lat0,lat1,idebug)
!
      call minor(tn,difk,o1,co,co2,n2,so2,so2_nm,so2_out,so2_nm_out,
     |  so2_loss,so2_prod,so2_lbc,so2_ubc,rmass_so2,phi_so2,
     |  alfa_so2,lev0,lev1,lon0,lon1,lat0,lat1,0)

!
        call addfld('SO2_OUT' ,' ',' ',so2_out(:,i0:i1,lat),
     |    'lon',lon0,lon1,'lev',lev0,lev1,lat)
        call addfld('SO2_TM1' ,' ',' ',so2_nm_out(:,i0:i1,lat),
     |    'lon',lon0,lon1,'lev',lev0,lev1,lat)
      end subroutine minor_so2
      end module so2_module
!
