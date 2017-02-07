!
      module so_module
      use params_module,only: nlevp1,nlonp4,nlat,dz
      implicit none
! -------------------------------------------------------------
! Boundary conditions, production and loss for SO are defined
! by comp_so, and referenced by minor_so. Comp_so is called
! from a latitude loop in dynamics. After comp_so, dynamics calls
! minor_so, which passes this module data to sub minor. Sub 
! minor contains 3d mpi calls and its own latitude loops.
!
!     real,dimension(nlonp4,nlat)   :: so_ubc ! upper boundary
!     real,dimension(nlonp4,3,nlat) :: so_lbc ! lower boundary
!     real,dimension(nlevp1,nlonp4,nlat) ::
!    |  so_prod,     ! production of so
!    |  so_loss      ! loss of so
!
! Boundary conditions and production and loss terms are allocated
! subdomains by sub alloc_so (called from allocdata.F).
! - last modified by SWB:  04/29/13
! 04/29/13 swb: comp_so.F and minor_so.F structure created
! 04/29/13 swb: P&L reactions and rates to be added
! 09/03/13 swb: phi_so same as for Ar with O, CO, CO2
! 09/03/13 swb: mass mixing ratio for so from 1-D KINETICS?
! 01/27/14 swb: Adding new reactions to SO and SO2 sources/sinks
! 02/07/14 swb: Updated 1-D KINETICS VMRs profiles of 17-species
!               (from best of Zhang et al. 2012)
! 05/06/14 swb: Added diagnostic plots/prints for production/loss terms
! 05/27/14 swb: Added s167 to SO losses  (on)
! 06/16/14 swb: Added s167 to SO losses  (on)
! 06/13/14 swb: Added s346 to SO losses  (off)
! 06/18/14 swb: ClCO3 VMR zeroed out to remove loss: s244
! 10/08/14 swb: s248 added.
! 04/15/15 swb: so lbc updated from Parkinson et al (2015)
! 04/15/15 swb: all source and sink diagnostics turned on
! -------------------------------------------------------------
!
      real,allocatable,dimension(:,:)   :: so_ubc ! upper boundary (i,j)
      real,allocatable,dimension(:,:,:) :: so_lbc ! lower boundary (i,3,j)
      real,allocatable,dimension(:,:,:) ::
     |  so_prod,     ! production of so (k,i,j)
     |  so_loss      ! loss of so       (k,i,j)
!  Mutual diffusion coefficients (normalized) for O-SO, CO-SO, CO2-SO
!  Use O-Ar, CO-Ar, CO2-Ar from mtgcm1.2/src: S. W. Bougher (9/3/2013)
!  Like Ar, set thermal diffusion coefficient = 0.0
      real :: phi_so(3) =  (/2.398, 1.875, 0.6631/)
      real,parameter :: alfa_so = 0. ! thermal diffusion coefficient
!  Mass mixing ratio lower boundary condition for SO (at ~69 km)
!  massmix = mix*mso/mbar = mix*48./44.
!  (1)  mix1 = 3.0E-09 ! Zhang et al (2011) Figure 8 eyeball value
!  (2)  mix2 = 2.1E-09 ! KINETICS reference simulation
!  (3)  mix3 = 2.6E-09 ! KINETICS reference simulation (Oct. 16, 2013 Tables)
!  (4)  mix4 = 3.032E-09 ! KINETICS reference simulation (Jan. 27, 2014 Tables)
!  (3)  mix5 = 3.3E-09 ! KINETICS reference simulation; Parkinson al (2015)
!     real :: pssob =  3.27E-09   !  From Zhang et al (2011) Figure 8 
!     real :: pssob =  2.29E-09   !  KINETICS reference simulation (update pls!)
!     real :: pssob =  2.83E-09   !  KINETICS reference simulation (updated 2013)
!     real :: pssob =  3.308E-09  !  KINETICS reference simulation (updated 2014)
      real :: pssob =  3.60E-09   !  KINETICS reference simulation: Parkinson (2015)
!
      contains
!-----------------------------------------------------------------------
      subroutine alloc_so(lon0,lon1,lat0,lat1)
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
      allocate(so_ubc(lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_so: error allocating',
     |  ' so_ubc: stat=',i3)") istat
      allocate(so_lbc(lon0:lon1,3,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_so: error allocating',
     |  ' so_lbc: stat=',i3)") istat
!
! Allocate subdomains to production and loss:
      allocate(so_prod(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_so: error allocating',
     |  ' so_prod: stat=',i3)") istat
      allocate(so_loss(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_so: error allocating',
     |  ' so_loss: stat=',i3)") istat

!     write(6,"('alloc_so: allocated module data')")
!
      end subroutine alloc_so
!
!-----------------------------------------------------------------------
      subroutine comp_so(tn,o1,co,co2,n2,o2,so,so2,barm,
     |  lev0,lev1,lon0,lon1,lat)
!
!-----------------------------------------------------------------------
! Advance SO (minor): 
! Arguments for comp_so (14).
! ----------------------------------------------------------------------
! Advance so by one time step. This is called from driver at
! each subdomain latitude.
! Last modified: swb (05/29/13)
!              : add new rates from chemrates_sox
!              : finalized sosrc, soloss terms
!              : next interp trace species from KINETICS to VTGCM grid
! ----------------------------------------------------------------------
!
      use qrj_module,only: qjso, qfjso, qjso2a, mk_xnmbari
      use fields_module,only: tlbc
      use cons_module,only: rmassinv_o1,rmassinv_co,rmassinv_n2,p0,
     |  rmass_so, rmass_so2, rmassinv_so, rmassinv_so2,
     |  expz,expzmid_inv,boltz,rmassinv_co2,rmassinv_o2,rmass_o1,
     |  rmass_co,rmass_o2,rmass_so,rmass_so2,gask,grav
      use chemrates_module,only:s114,s133,s166,s240,S244,s249,s258,
     |  s261,s288,s293,s294,s299,s300,s243,s286,s287,s167,s346,s248
      use addfld_module,only: addfld
!
!     write(6,"('enter comp_so: lat=',i2)") lat
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
!
! ----------------------------------------------------------------------
! Local:
      integer :: k,i,k0,k1,i0,i1,nk  ! for addfld (ask for check)
      real,dimension(lev0:lev1,lon0:lon1) :: 
     |  tni,            ! TN at interfaces
     |  qfjsoi,         ! Jso at interfaces (qrj module): loss freq (1/sec) units 
     |  qjsoi,          ! Jso at interfaces (qrj module): loss rate (#/cm3/sec) units 
     |  qjso2ai,        ! Jso2a->SO+O at interfaces (qrj module): prod rate (#/cm3/sec) units
     |  xnmbari,        ! xnmbar at interfaces (qrj module)
     |  xnmbar,         ! xnmbar at midlevels (qrj module)
     |  xtot,           ! Total number density (#/cm3)
     |  xo,             ! O1 number density (#/cm3)
     |  xco,            ! CO number density (#/cm3)
     |  xco2,           ! CO2 number density (#/cm3)
     |  xn2,            ! N2 number density (#/cm3)
     |  xo2,            ! O2 number density (#/cm3)
     |  xso,            ! SO number density (#/cm3)
     |  xso2,           ! SO2 number density (#/cm3)
     |  soprod1,        ! SO production (#/cm3/s)
     |  soprod2,        ! SO production (#/cm3/s)
     |  sops166,        ! SO production (#/cm3/s)
     |  sops287,        ! SO production (#/cm3/s)
     |  sops288,        ! SO production (#/cm3/s)
     |  soloss1,        ! SO loss1 (#/cm3/s)
     |  soloss2,        ! SO loss2 (#/cm3/s)
     |  sols240,        ! SO loss2 (#/cm3/s)
     |  sols248,        ! SO loss2 (#/cm3/s)
     |  sols249,        ! SO loss2 (#/cm3/s)
     |  sols133,        ! SO loss2 (#/cm3/s)
     |  sols167,        ! SO loss2 (#/cm3/s)
     |  sols243,        ! SO loss2 (#/cm3/s)
     |  sols299,        ! SO loss2 (#/cm3/s)
     |  sols244,        ! SO loss2 (#/cm3/s)
     |  sols346,        ! SO loss2 (#/cm3/s)
     |  xno2,           ! NO2 (#/cm3) 
     |  xcl,            ! CL (#/cm3) 
     |  xcl2,           ! CL2 (#/cm3) 
     |  xclo,           ! CLO (#/cm3) 
     |  xclco,          ! CLCO (#/cm3) 
     |  xclco3,         ! CLCO3 (#/cm3) 
     |  xclso2,         ! CLSO2 (#/cm3) 
     |  xoscl,          ! OSCl (#/cm3) 
     |  xs,             ! S (#/cm3) 
     |  xso3,           ! SO3 (#/cm3) 
     |  xh,             ! H (#/cm3) 
     |  xh2,            ! H2 (#/cm3)
     |  xoh,            ! OH (#/cm3) 
     |  xho2,           ! HO2 (#.cm3) 
     |  xo3,            ! O3 (#/cm3) 
     |  xhcl            ! hcl (#/cm3)
!
! From recent 1-D KINETICS calculation for nightside
! 29-MAY-2013 (private communication, C. Parkinson)
      real,dimension(lev0:lev1) ::
     |  no2,           ! NO2 (vmr) 
     |  cl,            ! CL (vmr) 
     |  cl2,           ! CL2 (vmr) 
     |  clo,           ! CLO (vmr) 
     |  clco,          ! CLCO (vmr) 
     |  clco3,         ! CLCO3 (vmr) 
     |  clso2,         ! CLSO2 (vmr) 
     |  oscl,          ! OSCl (vmr) 
     |  so3,           ! SO3 (vmr) 
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
! *******************************************************************
! Low Altitude Chemistry from Chris Parkinson Venus KINETICS code
! Nightside profiles in VMR: 22 (values) + 47 (zeros) = 69 levels
!  (07-FEB-2014 Update to VMR for 17 species)
! Alt rage = ~68 - ~112km (zp = -16 to -5.5)
! Revise for latest 1-D KINETICS Model Simulations
! -- Set top 47 levels with level = 22 VMR value to avoid discontinuity
! -- recast on VTGCM zp grid
! -- SO3 = fixed (always)
! -- SO2 = fixed (initial testing)
! -- SO  = fixed (initial testing)
! -- rest = fixed always
!
!********************************************************************
      no2(:) = (/
     |    5.925E-11,8.703E-11,8.508E-11,7.639E-11,
     |    6.780E-11,8.479E-11,1.368E-10,1.345E-10,
     |    5.828E-11,8.721E-12,1.743E-12,5.065E-13,
     |    9.805E-14,4.344E-14,1.959E-14,1.080E-14,
     |    5.716E-15,3.293E-15,1.841E-15,9.913E-16,
     |    4.973E-16,2.319E-16,0.00e+000,
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
! Species Number densities at interfaces (#/cm3)
! Given Qrj dissociation rates (#/cm3/sec) at midlevels from qrj
!
! ----------------------------------------------------------------------
! Lower boundary:
! so_lbc(:,1)=A, so_lbc(:,2)=B, so_lbc(:,3)=C define lower boundary 
!   condition where A*DPSX/DZ + B*PSX + C = 0.
!
      do i=lon0,lon1
!
! Value at bottom given by specified mass missing ratio (start)
        so_lbc(i,1,lat) = 0.
        so_lbc(i,2,lat) = 1.
        so_lbc(i,3,lat) = -pssob
!
! Upper boundary: Zero diffusive flux at top
!       so_ubc(i,lat) = 0.
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
! --------------------------------------------------------------------
!
! Calculate p0*e(-z)*barm/kT (1/(N*MBAR))at interfaces
!     subroutine mk_xnmbari(tni,barm,xnmbari,lev0,lev1,lon0,lon1,lat)
!
      call mk_xnmbari(tni,barm(:,lon0:lon1),
     |  xnmbari,
     |  lev0,lev1,lon0,lon1,lat)
! ----------------------------------------------------------------------
! Species Number densities at interfaces (#/cm3)
! Given Qrj dissociation rates (#/cm3/sec) at midlevels from qrj
! Calculate Qrj dissociation frequencies (#/sec) at interfaces: Jso, Jso2
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          xo(k,i)   = xnmbari(k,i)*o1(k,i)*rmassinv_o1
          xco(k,i)  = xnmbari(k,i)*co(k,i)*rmassinv_co
          xn2(k,i)  = xnmbari(k,i)*n2(k,i)*rmassinv_n2
          xo2(k,i) = xnmbari(k,i)*o2(k,i)*rmassinv_o2
          xco2(k,i) = xnmbari(k,i)*co2(k,i)*rmassinv_co2
	  if (xco2(k,i) < 0.0) xco2(k,i) = 0.0
          xso(k,i)  = xnmbari(k,i)*so(k,i)*rmassinv_so
	  if (xso(k,i) < 0.0) xso(k,i) = 0.0
          xso2(k,i)  = xnmbari(k,i)*so2(k,i)*rmassinv_so2
	  if (xso2(k,i) < 0.0) xso2(k,i) = 0.0
          qjsoi(k,i) = 0.5*(qjso(k,i,lat)+qjso(k+1,i,lat))
          qjso2ai(k,i) = 0.5*(qjso2a(k,i,lat)+qjso2a(k+1,i,lat))
          qfjsoi(k,i) = 0.5*(qfjso(k,i,lat)+qfjso(k+1,i,lat))
!
!   LOG10 DENSITIES FOR PLOTTING 
!         xco2l(k,i)=log10(xco2(k,i))
!         xcol(k,i)=log10(xco(k,i))
!         xn2l(k,i)=log10(xo2(k,i))
!         xol(k,i)=log10(xo(k,i))

!  Calcuklate Total Density
          xtot(k,i) = xo(k,i)+xco(k,i)+xco2(k,i)+xn2(k,i)
!
!  Converting vmr to number density
          xno2(k,i)   = xtot(k,i)*no2(k)
          xcl(k,i)   = xtot(k,i)*cl(k)
          xcl2(k,i)  = xtot(k,i)*cl2(k)
          xclo(k,i)  = xtot(k,i)*clo(k)
          xclco(k,i) = xtot(k,i)*clco(k)
          xclco3(k,i) = xtot(k,i)*clco3(k)
          xclso2(k,i) = xtot(k,i)*clso2(k)
          xoscl(k,i) = xtot(k,i)*oscl(k)
          xs(k,i) = xtot(k,i)*s(k)
          xso3(k,i) = xtot(k,i)*so3(k)
          xh(k,i)   = xtot(k,i)*h(k)
          xh2(k,i)  = xtot(k,i)*h2(k)
          xoh(k,i)   = xtot(k,i)*oh(k)
          xho2(k,i)  = xtot(k,i)*ho2(k)
          xo3(k,i) = xtot(k,i)*o3(k)
          xhcl(k,i) = xtot(k,i)*hcl(k)

        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!
! ----------------------------------------------------------------------
! (1)  Major Sources and Losses for SO (before NOX Chemistry Added)
!      From Zhang et al (2012)  Figures 12a-b.
! ---------------------------------------------------------------------
      do i=lon0,lon1
        do k=lev0,lev1-1
! ---------------------------------------------------------------------
!     SOURCES (#/cm3.sec) AND SINKS(#/sec): Prognostics
!     SOURCES (#/cm3.sec) AND SINKS(#/cm3/sec): Diagnostics
! ---------------------------------------------------------------------
!
          so_prod(k,i,lat) = qjso2ai(k,i) + 
     |       xs(k,i)*xo2(k,i)*s166        +              
     |       xo(k,i)*xoscl(k,i)*s287(k,i,lat) +              
     |       xcl(k,i)*xoscl(k,i)*s288

!         write(6,"('comp_so: k,i,lat=',3i4,' qjso2ai(k,i)=',e12.4,
!    |      ' xso2(k,i)=',e12.4,' xs(k,i)=',e12.4,' xo2(k,i)=',e12.4,
!    |      ' s166=',e12.4,' so_prod(k,i,lat)=',e12.4)") k,i,lat,
!    |      qjso2ai(k,i),xso2(k,i),xs(k,i),xo2(k,i),s166,
!    |      so_prod(k,i,lat)
!
          soprod1(k,i) = qjso2ai(k,i)
          sops166(k,i) = xs(k,i)*xo2(k,i)*s166
          sops287(k,i) = xo(k,i)*xoscl(k,i)*s287(k,i,lat)
          sops288(k,i) = xcl(k,i)*xoscl(k,i)*s288
!
!         soprod2(k,i) = xs(k,i)*xo2(k,i)*s166 +
!    |       xo(k,i)*xoscl(k,i)*s287(k,i,lat) +              
!    |       xcl(k,i)*xoscl(k,i)*s288
! --------------------------------------
!
	  so_loss(k,i,lat) = - (qfjsoi(k,i)     +
     |       xo(k,i)*xco2(k,i)*s240             +
     |       2.*xso(k,i)*xco2(k,i)*s249         +
     |       2.*xso(k,i)*s248(k,i,lat)          +
     |       xo2(k,i)*s167(k,i,lat)             +
     |       xclo(k,i)*s133                     +
     |       xcl(k,i)*xco2(k,i)*s243(k,i,lat)   +
     |       xclso2(k,i)*s299(k,i,lat)          +
     |       xclco3(k,i)*s244                   +
     |       xno2(k,i)*s346)
!
	  soloss1(k,i) = qjsoi(k,i)
	  sols240(k,i) = xo(k,i)*xso(k,i)*xco2(k,i)*s240
          sols249(k,i) = 2.*xso(k,i)*xso(k,i)*xco2(k,i)*s249
	  sols248(k,i) = 2.*xso(k,i)*xso(k,i)*s248(k,i,lat)
	  sols133(k,i) = xclo(k,i)*xso(k,i)*s133
	  sols167(k,i) = xo2(k,i)*xso(k,i)*s167(k,i,lat)
	  sols243(k,i) = xcl(k,i)*xso(k,i)*xco2(k,i)*s243(k,i,lat)
	  sols299(k,i) = xclso2(k,i)*xso(k,i)*s299(k,i,lat)
	  sols244(k,i) = xclco3(k,i)*xso(k,i)*s244
          sols346(k,i) = xno2(k,i)*xso(k,i)*s346
!
!        soloss2(k,i) = xo(k,i)*xso(k,i)*xco2(k,i)*s240   +
!       |       xso(k,i)*xso(k,i)*xco2(k,i)*s249              +
!       |       xclo(k,i)*xso(k,i)*s133                       +
!       |       xcl(k,i)*xso(k,i)*xco2(k,i)*s243(k,i,lat)     +
!       |       xclso2(k,i)*xso(k,i)*s299(k,i,lat)            +
!       |       xclco3(k,i)*xso(k,i)*s244
!
!  Convert to Log10 units for plotting diagnostic terms : initial
        if (soprod1(k,i) .GT. 1.0) THEN
	    soprod1(k,i)= alog10(soprod1(k,i))
	 else
	    soprod1(k,i) = 0.
	 endif
!       if (soprod2(k,i) .GT. 1.0) THEN
!           soprod2(k,i)= alog10(soprod2(k,i))
!       else
!           soprod2(k,i) = 0.
!       endif
        if (soloss1(k,i) .GT. 1.0) THEN
	    soloss1(k,i)= alog10(soloss1(k,i))
	 else
	    soloss1(k,i) = 0.
	 endif
!       if (soloss2(k,i) .GT. 1.0) THEN
!          soloss2(k,i)= alog10(soloss2(k,i))
!       else
!          soloss2(k,i) = 0.
!       endif
!  Convert to Log10 units for plotting diagnostic terms : expanded
!  Production of SO (#/cm3/s)
         if (sops166(k,i) .GT. 1.0) THEN
	    sops166(k,i)= alog10(sops166(k,i))
	 else
	    sops166(k,i) = 0.
	 endif
         if (sops287(k,i) .GT. 1.0) THEN
	    sops287(k,i)= alog10(sops287(k,i))
	 else
	    sops287(k,i) = 0.
	 endif
         if (sops288(k,i) .GT. 1.0) THEN
	    sops288(k,i)= alog10(sops288(k,i))
	 else
	    sops288(k,i) = 0.
	 endif
!  Loss of SO (#/cm3/s)
         if (sols240(k,i) .GT. 1.0) THEN
	    sols240(k,i)= alog10(sols240(k,i))
	 else
	    sols240(k,i) = 0.
	 endif
         if (sols248(k,i) .GT. 1.0) THEN
	    sols248(k,i)= alog10(sols248(k,i))
	 else
	    sols248(k,i) = 0.
	 endif
         if (sols249(k,i) .GT. 1.0) THEN
           sols249(k,i)= alog10(sols249(k,i))
         else
           sols249(k,i) = 0.
         endif
         if (sols133(k,i) .GT. 1.0) THEN
	    sols133(k,i)= alog10(sols133(k,i))
	 else
	    sols133(k,i) = 0.
	 endif
         if (sols167(k,i) .GT. 1.0) THEN
	    sols167(k,i)= alog10(sols167(k,i))
	 else
	    sols167(k,i) = 0.
	 endif
         if (sols243(k,i) .GT. 1.0) THEN
	    sols243(k,i)= alog10(sols243(k,i))
	 else
	    sols243(k,i) = 0.
	 endif
         if (sols299(k,i) .GT. 1.0) THEN
	    sols299(k,i)= alog10(sols299(k,i))
	 else
	    sols299(k,i) = 0.
	 endif
         if (sols244(k,i) .GT. 1.0) THEN
	    sols244(k,i)= alog10(sols244(k,i))
	 else
	    sols244(k,i) = 0.
	 endif
         if (sols346(k,i) .GT. 1.0) THEN
            sols346(k,i)= alog10(sols346(k,i))
         else
            sols346(k,i) = 0.
        endif
!
!  Make sure all model P(+) and L(-) terms are correct sign :  test option only! 
!       if (so_prod(k,i,lat) .LT. 0.0) so_prod(k,i,lat) = 0.0
!       if (so_loss(k,i,lat) .GT. 0.0) so_loss(k,i,lat) = 0.0
!
         enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
! ----------------------------------------------------------------------
!   * SO Production and Losses
!
!  Test 3: OK all positive
!      write(6,"('comp_so after calculate prod,loss: 
!    |	  soprod1=',2e14.4, ' soprod2=',2e14.4,
!    |  ' soloss1=',2e14.4,' soloss2=',2e14.4)")
!    |  minval(soprod1),maxval(soprod1),
!    |  minval(soprod2),maxval(soprod2),
!    |  minval(soloss1),maxval(soloss1),
!    |  minval(soloss2),maxval(soloss2)
!
!  Test 4: OK, correct signs
!     write(6,"('comp_so2 after calculate prod,loss2:
!    |   so_prod=',2e12.4, ' so_loss=',2e12.4)")
!    |  minval(so_prod),maxval(so_prod),
!    |  minval(so_loss),maxval(so_loss)
!
!
      call addfld('SOPROD1',' ',' ',soprod1,'lev',k0,k1,'lon',i0,i1,
     |  lat)
      call addfld('SOPS166',' ',' ',sops166,'lev',k0,k1,'lon',i0,i1,
     |  lat)
      call addfld('SOPS287',' ',' ',sops287,'lev',k0,k1,'lon',i0,i1,
     |  lat)
      call addfld('SOPS288',' ',' ',sops288,'lev',k0,k1,'lon',i0,i1,
     |  lat)
!     call addfld('SOPROD2',' ',' ',soprod2,'lev',k0,k1,'lon',i0,i1,
!    |  lat)
      call addfld('SOLOSS1',' ',' ',soloss1,'lev',k0,k1,'lon',i0,i1,
     |  lat)
      call addfld('SOLS240',' ',' ',sols240,'lev',k0,k1,'lon',i0,i1,
     |  lat)
      call addfld('SOLS248',' ',' ',sols248,'lev',k0,k1,'lon',i0,i1,
     |  lat)
      call addfld('SOLS249',' ',' ',sols249,'lev',k0,k1,'lon',i0,i1,
     |  lat)
      call addfld('SOLS133',' ',' ',sols133,'lev',k0,k1,'lon',i0,i1,
     |  lat)
      call addfld('SOLS167',' ',' ',sols167,'lev',k0,k1,'lon',i0,i1,
     |  lat)
      call addfld('SOLS243',' ',' ',sols243,'lev',k0,k1,'lon',i0,i1,
     |  lat)
      call addfld('SOLS299',' ',' ',sols299,'lev',k0,k1,'lon',i0,i1,
     |  lat)
      call addfld('SOLS244',' ',' ',sols244,'lev',k0,k1,'lon',i0,i1,
     |  lat)
      call addfld('SOLS346',' ',' ',sols346,'lev',k0,k1,'lon',i0,i1,
     |  lat)
!     call addfld('SOLOSS2',' ',' ',soloss2,'lev',k0,k1,'lon',i0,i1,
!    |  lat)
!
      end subroutine comp_so
!-----------------------------------------------------------------------
      subroutine minor_so(tn,difk,o1,co,co2,n2,so,so_nm,so_out,
     |  so_nm_out,lev0,lev1,lon0,lon1,lat0,lat1)

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
     |  so,      ! so (mmr)
     |  so_nm    ! so at time n-1 
!
! -----------------------------------------------------------------------
! Output args:
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(out) ::
     |  so_out,    ! SO output 
     |  so_nm_out  ! SO output at time n-1
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
!     subroutine minor(tn,difk,o2,o1,fcomp,fcomp_tm1,fcomp_out,
!    |  fcomp_tm1_out,sloss,sprod,flbc,fubc,rmx,phix,alfax,
!    |  lev0,lev1,lon0,lon1,lat0,lat1,idebug)
!
      call minor(tn,difk,o1,co,co2,n2,so,so_nm,so_out,so_nm_out,
     |  so_loss,so_prod,so_lbc,so_ubc,rmass_so,phi_so,
     |  alfa_so,lev0,lev1,lon0,lon1,lat0,lat1,0)

!
!       call addfld('SO_OUT' ,' ',' ',so_out(:,i0:i1,lat),
!    |    i0,i1,k0,k1,lat)
!       call addfld('SO_TM1' ,' ',' ',so_nm_out(:,i0:i1,lat),
!    |    i0,i1,k0,k1,lat)
      end subroutine minor_so
      end module so_module
