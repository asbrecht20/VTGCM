!
      subroutine comp_oco(tn,barm,o1,co,co2,n2,o2,op,co2p,o2p,ne,
     |  lev0,lev1,lon0,lon1,lat)
!
! Calculates fs array which gives sources and sinks for comp.F (o,co). 
! Last modified: swb (11/20/06)  
!              : add new chemrates nomenclature for rates
!              : add osrc, oloss, cosrc, coloss terms
!              : populate the fs array for comp.F (o,co)
! Last modified: swb (11/21/06)  
!              : remove IO as loss of O (minor)
!              : mid-level calculation for S and L only (plots)
! Next modified: bf (11/XX06)  add cross terms to fs  
! Last modified : swb (10/30/07) log10 sources, losses for diagnostics
! Last modified : swb (11/02/07) OLIFE and COLIFE (chemical timescales)
! Last modified : swb (12/03/07) SOURCES improved with ionden.F changes and
!                                no OSRC or COSRC below z = -10.0
! last modified; swb (12/3/07) Insert DY-NT flags for OSRC,COSRC,OLOSS,COLOSS
! last modified; swb (04/22/08) remove DY-NT flags for OSRC,COSRC,OLOSS,COLOSS
!              :     (04/22/08) ions relaxed to go to zero; electron cutoff OK
! last modified; swb (04/24/08) OPCE and COPCE secondaries put out
! last modified; swb (05/05/08) fs-array segregated terms implemented like tiegcm#2
! last modified; swb (10/31/13) cosrc (from JCO2) is doubled for testing.
! last modified; swb (11/13/13) cosrc (from JCO2) is standard for production
! last modified: swb (05/08/14) r2(k,i,lat) and r5(k,i,lat) updated (Zhang et al., 2012)
! last modified; swb (02/11/16) cosrc and coloss updated from Zhang et al (2012)
!                              (checked by C. D. Parkinson and confirmed reaction set)
!                              (trace species sources & losses updated below 110 km)
!                              (trace species sourcs & loss terms on diagnostics added)
! last modified; swb (02/11/16) coloss and oloss terms in fs array modified  (trial)
! last modified; swb (02/17/16) colossG updated from Zhang et al (2012)
!
      use params_module,only: dz
      use cons_module,only: rmassinv_o1,rmassinv_co,rmassinv_n2,p0,
     |  expz,boltz,rmassinv_co2,rmassinv_o2,rmass_o1,rmass_co,
     |  gask,grav
      use qrj_module,only: rj,qjo2,
     |  mk_xnmbari
      use fields_module,only: tlbc
      use chemrates_module,only:a13,a14,r2,r3,r5,r18,r19,r20,rk20,
     |  rk28,r104,r49,r44,r66,r108,r43,r99,r102a,r102b,r25,r144,r146,
     |  r147, r160, r111, r162, r163, fs
      use chapman_module,only: idn_lonlat
      use addfld_module,only: addfld
      use o2_module,only: integral
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
!
! All input arrays are at current time step (itp): 
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) ::
     |  tn,    ! neutral temperature (deg K)
     |  barm,  ! mean molecular weight
     |  o1,    ! o1 (mmr)
     |  co,    ! co (mmr)
     |  co2,   ! co2 (mmr)
     |  n2,    ! n2 (mmr)
     |  o2,    ! o2 (mmr)
     |  op,    ! O+ (#/cm3)
     |  co2p,  ! CO2+ (#/cm3)
     |  o2p,   ! O2+ (#/cm3)
     |  ne     ! electron density  (#/cm3)
!
! Local:
      integer :: k,i,nk,k0,k1,i0,i1
!
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  tni,            ! TN at interfaces
     |  rji,            ! JCO2 at interfaces (qrj module)
     |  rjis,            ! JCO2 at interfaces (qrj module)[/sec]
     |  qjo2i,          ! JO2 at interfaces (qrj module)
     |  qjo2is,          ! JO2 at interfaces (qrj module)[/sec]
     |  xnmbari,        ! xnmbar at interfaces (qrj module)
     |  xnmbar,         ! xnmbar at midlevels (qrj module)
     |  cosrc,          ! CO production terms (#/cm3/sec)
     |  cosrc2,         ! CO production terms (#/cm3/sec)
     |  cosrc3,         ! CO production term (#/cm3/sec)
     |  coloss,         ! CO loss terms (1/sec)
     |  coloss2,        ! CO loss terms (#/cm3/sec)
     |  cotime,         ! CO chemical loss timescale  (Earth days)
     |  cotimel,        ! CO chemical loss timescale  (log10 Earth days)
     |  osrc,           ! O1 production terms  (#/cm3/sec)
     |  osrcA,          ! O1 production terms  (#/cm3/sec)
     |  osrcB,          ! O1 production terms  (#/cm3/sec)
     |  osrcC,          ! O1 production terms  (#/cm3/sec)
     |  osrcD,          ! O1 production terms  (#/cm3/sec)
     |  osrcT,          ! Total O1 production terms  (#/cm3/sec)
     |  oloss,          ! O1 loss terms  (1/sec)
     |  oloss2,         ! O1 loss terms  (#/cm3/sec)
     |  oloss3,         ! O1 loss terms  (#/cm3/sec) 
     |  cosrcA,          ! CO production terms  (#/cm3/sec)
     |  cosrcB,          ! CO production terms  (#/cm3/sec)
     |  cosrcC,          ! CO production terms  (#/cm3/sec)
     |  cosrcD,          ! CO production terms  (#/cm3/sec)
     |  cosrcE,          ! CO production terms  (#/cm3/sec)
     |  cosrcF,          ! CO production terms  (#/cm3/sec)
     |  cosrcG,          ! CO production terms  (#/cm3/sec)
     |  cosrcH,          ! CO production terms  (#/cm3/sec)
     |  colossA,         ! CO loss terms  (#/cm3/sec) 
     |  colossB,         ! CO loss terms  (#/cm3/sec) 
     |  colossC,         ! CO loss terms  (#/cm3/sec) 
     |  colossD,         ! CO loss terms  (#/cm3/sec) 
     |  colossE,         ! CO loss terms  (#/cm3/sec) 
     |  colossF,         ! CO loss terms  (#/cm3/sec) 
     |  colossG,         ! CO loss terms  (#/cm3/sec) 
     |  olossA,         ! O1 loss terms  (#/cm3/sec) 
     |  olossB,         ! O1 loss terms  (#/cm3/sec) 
     |  olossC,         ! O1 loss terms  (#/cm3/sec) 
     |  olosDd,         ! O1 loss terms  (#/cm3/sec) 
     |  olosEd,         ! O1 loss terms  (#/cm3/sec) 
     |  olosFd,         ! O1 loss terms  (#/cm3/sec) 
     |  olosGd,         ! O1 loss terms  (#/cm3/sec) 
     |  olosHd,         ! O1 loss terms  (#/cm3/sec) 
     |  olosId,         ! O1 loss terms  (#/cm3/sec) 
     |  olosJd,         ! O1 loss terms  (#/cm3/sec) 
     |  olosKd,         ! O1 loss terms  (#/cm3/sec) 
     |  olosLd,         ! O1 loss terms  (#/cm3/sec) 
     |  olosTrdy,       ! Total O1 loss terms by trace species  (#/cm3/sec) 
     |  olosDn,         ! O1 loss terms  (#/cm3/sec) 
     |  olosEn,         ! O1 loss terms  (#/cm3/sec) 
     |  olosFn,         ! O1 loss terms  (#/cm3/sec) 
     |  olosGn,         ! O1 loss terms  (#/cm3/sec) 
     |  olosHn,         ! O1 loss terms  (#/cm3/sec) 
     |  olosIn,         ! O1 loss terms  (#/cm3/sec) 
     |  olosJn,         ! O1 loss terms  (#/cm3/sec) 
     |  olosKn,         ! O1 loss terms  (#/cm3/sec) 
     |  olosLn,         ! O1 loss terms  (#/cm3/sec) 
     |  olosTrnt,       ! Total O1 loss terms by trace species  (#/cm3/sec) 
     |  otime,          ! O1 chemical loss timescale  (Earth days)
     |  otimel,         ! O1 chemical loss timescale  (log10 Earth days)
     |  xo,             ! O1 number density (#/cm3)
     |  xopce,          ! O1 number density (#/cm3): PCE secondary
     |  xco,            ! CO number density (#/cm3)
     |  xcopce,         ! CO number density (#/cm3): PCE secondary
     |  xco2,           ! CO2 number density (#/cm3)
     |  xn2,            ! N2 number density (#/cm3)
     |  xo2,            ! O2 number density (#/cm3)
     |  xtot,           ! Total number density (#/cm3)
     |  xol,            ! LOG10-O1 number density (#/cm3)
     |  xcol,           ! LOG10-CO number density (#/cm3)
     |  xco2l,          ! LOG10-CO2 number density (#/cm3)
     |  xn2l,           ! LOG10-N2 number density (#/cm3)
     |  xo2l,           ! LOG10-O2 number density (#/cm3)
     |  slab,           ! slab width (cm)  
     |  onet,           ! O net production (osrc - oloss)
     |  oprodcol,       ! O production integrated column (#/cm2/sec)
     |  onetprodcol,    ! O production integrated column (#/cm2/sec)
     |  odestcol        ! O destruction integrated column (#/cm2/sec)
!
! -----------------------------------------------------------------------
! Low Altitude Tracers
! Specificed VMR specifically for SOx, CLOx, HOx chemistry 
! From recent 1-D KINETICS calculation for nightside
! 29-MAY-2013 (private communication, C. Parkinson)
! From recent 1-D KINETICS calculation for nightside
! 28-JAN-2014 (private communication, C. Parkinson)
      real,dimension(lev0:lev1) ::
     |  cl,            ! CL (vmr) 
     |  cl2,           ! CL2 (vmr) 
     |  clo,           ! CLO (vmr) 
     |  clco,          ! CLCO (vmr) 
     |  clco3,         ! CLCO3 (vmr) 
     |  clso2,         ! CLSO2 (vmr) 
     |  oscl,          ! OSCl (vmr) 
     |  h,             ! H (vmr) 
     |  h2,            ! H2 (vmr)
     |  oh,            ! OH (vmr) 
     |  ho2,           ! HO2 (vmr) 
     |  o3,            ! O3 (vmr) 
     |  hcl            ! hcl (vmr)
!
!   LOG10 DENSITIES FOR PLOTTING 
!
! Low Altitude Tracers  
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  xcl,            ! CL (density)
     |  xcl2,           ! CL2 (density)
     |  xclo,           ! CLO (density)
     |  xclco,          ! CLCO (density)
     |  xclco3,         ! CLCO3 (density)
     |  xclso2,         ! CLSO2 (density)
     |  xoscl,          ! OSCl (density)
     |  xh,             ! H (density)
     |  xh2,            ! H2 (density)
     |  xoh,            ! OH (density)
     |  xho2,           ! HO2 (density)
     |  xo3,            ! O3 (density)
     |  xhcl            ! hcl (density)
!
      real,dimension(lev0:lev1) ::
     | fldht,colout
#include<defs.h>
!
! For addfld calls:
      k0=lev0 
      k1=lev1
      i0=lon0
      i1=lon1
!
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
#if (NLEV==68)
!********************************************************************
! Single Res. data
!********************************************************************
! Low Altitude Chemistry from Chris Parkinson Venus KINETICS code
! Profiles in VMR: 22 (values) + 47 (zeros) = 69 levels
!  (07-FEB-2014 Update to VMR for 17 species)
!  (11-FEB-2016 confirmed OK from Parkinson for 13 species)
! Profiles are identical to Zhang et al (2012) Final profiles!
! Alt range = ~69 - ~112km (zp = -16.0 to -5.5)
! Revised for latest 1-D KINETICS Model Simulations
! -- Set top 46 levels with level = 23 VMR value to avoid discontinuity
! -- recast on VTGCM zp grid
! -- SO3 = fixed (always)
! -- SO2 = from history (fixed to start)
! -- SO  = from history (fixed to start)
! -- rest = fixed always
! -- ClCO3 VMR can be zeroed out to remove P/L terms:  s244, s261
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
#elif (NLEV==136)
!********************************************************************
! Double Res. data
!********************************************************************
!! --------------------------------------------------------------------
! Low Altitude Chemistry extracted from Krasnopolsky (2010) for the NIGHTSIDE 
! NIGHTSIDE profiles in MMR (paper has number density in Krasn.m they were converted to mmr)
!  (05-19-10)
! Alt rage = ~70km - ~112km (zp = -16 to -5)
!********************************************************************
      cl2n(:) = (/
     |  0.,0.,0.,0.,0.,0.,0.,0.,
     |  3.59E-11,8.09E-11,9.88E-11,1.10E-10,1.23E-10,1.35E-10,
     |  1.45E-10,1.55E-10,1.63E-10,1.71E-10,1.72E-10,1.71E-10,
     |  1.44E-10,1.09E-10,8.55E-11,6.44E-11,5.42E-11,4.67E-11,
     |  3.94E-11,3.23E-11,2.69E-11,2.40E-11,2.11E-11,1.89E-11,
     |  1.73E-11,1.48E-11,1.16E-11,8.82E-12,7.66E-12,6.50E-12,
     |  5.88E-12,5.39E-12,4.86E-12,4.24E-12,3.62E-12,2.11E-12,
     |  6.03E-13,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./)
!
      clcon(:) = (/
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  1.38E-13,3.12E-13,2.19E-12,4.50E-12,3.96E-12,2.71E-12,
     |  1.66E-12,6.64E-13,3.54E-13,2.17E-13,1.09E-13,1.09E-13,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./)
!
      cln(:) = (/
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  1.10e-13,2.47e-13,3.97e-13,5.39e-13,6.29e-13,7.20e-13,
     |  8.42e-11,1.88e-10,1.31e-09,2.68e-09,3.22e-09,3.55e-09,
     |  3.89e-09,4.23e-09,4.32e-09,4.35e-09,4.37e-09,4.39e-09,
     |  4.44e-09,4.55e-09,4.65e-09,4.75e-09,4.84e-09,4.81e-09,
     |  4.71e-09,4.65e-09,4.75e-09,4.85e-09,4.84e-09,4.80e-09,
     |  4.76e-09,4.78e-09,4.78e-09,5.17e-09,5.56e-09,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./)
!
      clon(:) = (/
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  1.93e-09,1.93e-09,1.93e-09,2.16e-09,2.41e-09,2.65e-09,2.85e-09,
     |  3.05e-09,3.12e-09,3.13e-09,2.02e-09,6.32e-10,2.68e-10,1.60e-10,
     |  1.08e-10,7.07e-11,5.23e-11,3.86e-11,2.75e-11,2.04e-11,1.48e-11,
     |  1.18e-11,8.72e-12,6.37e-12,4.49e-12,3.20e-12,2.31e-12,1.52e-12,
     |  1.19e-12,8.46e-13,5.19e-13,5.19e-13,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./)
!
      h2n(:) = (/
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  5.12e-11,1.15e-10,1.28e-10,1.29e-10,1.30e-10,1.31e-10,1.34e-10,
     |  1.37e-10,1.43e-10,1.48e-10,1.50e-10,1.51e-10,1.53e-10,1.55e-10,
     |  1.59e-10,1.62e-10,1.65e-10,1.66e-10,1.67e-10,1.67e-10,1.67e-10,
     |  1.68e-10,1.70e-10,1.71e-10,1.71e-10,1.72e-10,1.77e-10,1.81e-10,
     |  1.83e-10,1.84e-10,1.85e-10,1.85e-10,1.84e-10,1.94e-10,2.03e-10,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./)
!
      ho2n(:) = (/
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  2.51e-15,2.51e-15,2.51e-15,4.50e-15,7.68e-15,1.10e-14,1.38e-14,
     |  1.66e-14,2.90e-14,4.38e-14,5.04e-13,1.07e-12,9.45e-13,6.39e-13,
     |  4.44e-13,2.77e-13,2.01e-13,1.49e-13,1.04e-13,7.14e-14,4.43e-14,
     |  2.59e-14,7.39e-15,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0., 
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./)
!     
      hcln(:) = (/
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  9.47e-08,2.13e-07,2.40e-07,2.43e-07,2.44e-07,2.45e-07,2.50e-07,
     |  2.55e-07,2.62e-07,2.68e-07,2.68e-07,2.66e-07,2.70e-07,2.75e-07,
     |  2.69e-07,2.60e-07,2.61e-07,2.77e-07,2.84e-07,2.79e-07,2.75e-07,
     |  2.74e-07,2.77e-07,2.79e-07,2.82e-07,2.87e-07,2.98e-07,3.09e-07,
     |  2.99e-07,2.84e-07,2.77e-07,2.84e-07,2.90e-07,2.93e-07,2.97e-07,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./)
!
      o3n(:) = (/
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  7.39e-11,7.39e-11,7.39e-11,9.04e-11,1.19e-10,1.53e-10,1.92e-10,
     |  2.32e-10,3.02e-10,3.81e-10,5.59e-10,7.63e-10,1.07e-09,1.41e-09,
     |  1.98e-09,2.60e-09,3.38e-09,4.20e-09,5.37e-09,7.06e-09,8.88e-09,
     |  1.09e-08,1.29e-08,1.67e-08,2.16e-08,2.76e-08,3.43e-08,4.21e-08,
     |  5.45e-08,6.70e-08,8.64e-08,1.08e-07,1.35e-07,1.71e-07,2.07e-07,
     |  2.97e-07,3.87e-07,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./)
!
      ohn(:) = (/
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  1.39e-13,3.12e-13,1.77e-12,3.54e-12,4.91e-12,6.18e-12,6.42e-12,
     |  6.40e-12,6.25e-12,5.90e-12,5.54e-12,5.18e-12,4.82e-12,4.49e-12,
     |  4.18e-12,3.80e-12,3.37e-12,2.97e-12,2.69e-12,2.40e-12,2.12e-12,
     |  1.83e-12,1.58e-12,1.38e-12,1.18e-12,6.88e-13,1.97e-13,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./)
#else
 write(6,"('ERROR: unsupported NLEV = ',i6)") nlev
#endif
! --------------------------------------------------------------------
! --------------------------------------------------------------------
! Low Altitude Photochemical Sinks for O and CO from
! Yung and Demore (1982) and Massie et al. (1983)
! Possible range:  ~70 to 100 km or zp=-16 to -8

!**** Yung and Demore(1982) Model C Diurnal Average Values
!23456789012345678901234567890123456789012345678901234567890123456789012
! Data Sets
!
!     xcld(:) = (/
!    |  2.4289E+10,2.4185E+10,2.2817E+10,2.078E+10,1.7739E+10,
!    |  1.4699E+10,1.1801E+10,8.9977E+09,6.1946E+09,5.0966E+09,
!    |  3.9987E+09,3.0E+09,2.0E+09,1.5E+09,1.0E+09,7.0E+08,
!    |  5.0E+08,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
!    |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
!    |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
!    |  0.,0.,0.,0.,0.,0.,0.,0./)
!
!     xclod(:) = (/
!    |  1.2854E+08,1.3318E+08,1.1344E+08,9.4261E+07,7.592E+07,
!    |  5.758E+07,4.1573E+07,2.7122E+07,1.2671E+07,9.2612E+06,
!    |  5.8513E+06,2.4414E+06,1.6229E+06,8.0433E+05,4.0E+05,
!    |  1.5E+05,40000.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
!    |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
!    |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
!    |  0.,0.,0.,0.,0.,0.,0.,0.,0./)
!
!     xclod(:) = (/
!    |  1.2854E+08,1.3318E+08,1.1344E+08,9.4261E+07,7.592E+07,
!    |  5.758E+07,4.1573E+07,2.7122E+07,1.2671E+07,9.2612E+06,
!    |  5.8513E+06,2.4414E+06,1.6229E+06,8.0433E+05,4.0E+05,
!    |  1.5E+05,40000.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
!    |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
!    |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
!    |  0.,0.,0.,0.,0.,0.,0.,0.,0./)
!
!     xclcod(:) = (/
!    |  4.2932E+06,4.4241E+06,4.4802E+06,4.5464E+06,4.6278E+06,
!    |  4.7092E+06,5.7175E+06,7.3439E+06,8.9702E+06,7.8645E+06,
!    |  6.7587E+06,5.6529E+06,3.0E+06,1.5E+06,1.0E+06,4.0E+05,
!    |  1.0E+05,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
!    |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
!    |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
!    |  0.,0.,0.,0.,0.,0.,0.,0./)
!
!     xclcod(:) = (/
!    |  4.2932E+06,4.4241E+06,4.4802E+06,4.5464E+06,4.6278E+06,
!    |  4.7092E+06,5.7175E+06,7.3439E+06,8.9702E+06,7.8645E+06,
!    |  6.7587E+06,5.6529E+06,3.0E+06,1.5E+06,1.0E+06,4.0E+05,
!    |  1.0E+05,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
!    |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
!    |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
!    |  0.,0.,0.,0.,0.,0.,0.,0./)
!
!     xcl2d(:) = (/
!    |  4.8998E+10,2.8886E+10,1.6482E+10,8.0733E+09,5.6584E+09,
!    |  3.2436E+09,1.8366E+09,1.1017E+09,3.6671E+08,2.6566E+08,
!    |  1.646E+08,6.3551E+07,4.0762E+07,1.7974E+07,8.0E+06,
!    |  2.0E+06,5.0E+05,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
!    |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
!    |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
!    |  0.,0.,0.,0.,0.,0.,0.,0./)
!
!     xohd(:) = (/
!    |  1296.9,1591.4,1789.3,1964.8,2106.7,2248.7,2676.8,3295.7,3914.6,
!    |  5360.8,6807.,8253.2,12717.,16000.,27000.,40000.,55000.,0.,0.,0.,
!    |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
!    |  0.,0.,0.,0.,0.,0.,
!    |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
!    |  0.,0.,0.,0.,0.,0.,0.,0./)
!
!     xho2d(:) = (/
!    |  30181.,39989.,50091.,57641.,61361.,65082.,67636.,69413.,71190.,
!    |  72677.,74165.,75652.,70000.,65000.,60000.,40000.,35000.,0.,0.,
!    |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
!    |  0.,0.,0.,0.,0.,0.,0.,
!    |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
!    |  0.,0.,0.,0.,0.,0.,0.,0./)
!
!     xo3d(:) = (/
!    |  1.5752E+07,1.343E+07,1.0675E+07,8.9911E+06,8.9139E+06,8.8368E+6,
!    |  9.0825E+06,9.5435E+06,1.0004E+07,1.1116E+07,1.2227E+07,
!    |  1.3338E+07,1.4966E+07,1.6594E+07,1.687E+07,1.5121E+07,1.3371E+7,
!    |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
!    |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
!    |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
!    |  0.,0.,0.,0.,0.,0.,0.,0./)
!
!     hcld(:) = (/
!    |  7.8834E-07,7.872E-07,7.9038E-07,7.9129E-07,7.8879E-07,7.863E-07,
!    |  7.8251E-07,7.7786E-07,7.7322E-07,7.589E-07,7.4458E-07,7.3026E-7,
!    |  7.1416E-07,6.9806E-07,6.7835E-07,6.5323E-07,6.2811E-07,0.,0.,0.,
!    |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
!    |  0.,0.,0.,0.,0.,0.,
!    |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
!    |  0.,0.,0.,0.,0.,0.,0.,0./)
!
!     h2d(:) = (/
!    |  0.,0.,0.,0.,0.,0.,9.2107E-10,1.9621E-09,3.0031E-09,7.4297E-09,
!    |  1.1856E-8,1.6283E-8,2.5725E-08,3.5166E-08,4.6973E-08,6.2328E-08,
!    |  7.7682E-008,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
!    |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
!    |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
!    |  0.,0.,0.,0.,0.,0.,0.,0./)
!---------------------------------------------------------------------
!
! Calculate p0*e(-z)*barm/kT (1/(N*MBAR))at interfaces
!     subroutine mk_xnmbari(tni,barm,xnmbari,lev0,lev1,lon0,lon1,lat)
!
      call mk_xnmbari(tni,barm(:,lon0:lon1),
     |  xnmbari,
     |  lev0,lev1,lon0,lon1,lat)
!
! Species Number densities at interfaces (#/cm3)
! Qrj dissociation rates (#/cm3.sec) at interfaces: JCO2, JO2
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          xo(k,i)   = xnmbari(k,i)*o1(k,i)*rmassinv_o1
          xco(k,i)  = xnmbari(k,i)*co(k,i)*rmassinv_co
          xco2(k,i) = xnmbari(k,i)*co2(k,i)*rmassinv_co2
	  if (xco2(k,i) < 1.0) xco2(k,i) = 1.0
          xo2(k,i)  = xnmbari(k,i)*o2(k,i)*rmassinv_o2
          xn2(k,i)  = xnmbari(k,i)*n2(k,i)*rmassinv_n2
          rji  (k,i) = 0.5*(rj(k,i,lat)+rj(k+1,i,lat))
          rjis (k,i) = rji(k,i)/xco2(k,i) 
          qjo2i(k,i) = 0.5*(qjo2(k,i,lat)+qjo2(k+1,i,lat))
          qjo2is(k,i) = qjo2i(k,i)/xo2(k,i) 
!
! Creating total density
!         xtot(k,i) = xo(k,i)+xco(k,i)+xco2(k,i)+xo2(k,i)+xn2(k,i)
!         xtot(k,i) = xo(k,i)+xco(k,i)+xco2(k,i)+xn2(k,i)
!         xtot(k,i) = o1(k,i)+co(k,i)+co2(k,i)+n2(k,i)
          xtot(k,i) = xnmbari(k,i)/barm(k,i)

!  Converting vmr to number density
          xcl(k,i)   = xtot(k,i)*cl(k)
          xcl2(k,i)  = xtot(k,i)*cl2(k)
          xclo(k,i)  = xtot(k,i)*clo(k)
          xclco(k,i) = xtot(k,i)*clco(k)
          xclco3(k,i) = xtot(k,i)*clco3(k)
          xclso2(k,i) = xtot(k,i)*clso2(k)
          xoscl(k,i) = xtot(k,i)*oscl(k)
          xh(k,i)   = xtot(k,i)*h(k)
          xh2(k,i)  = xtot(k,i)*h2(k)
          xoh(k,i)   = xtot(k,i)*oh(k)
          xho2(k,i)  = xtot(k,i)*ho2(k)
          xo3(k,i) = xtot(k,i)*o3(k)
          xhcl(k,i) = xtot(k,i)*hcl(k)
!
!   LOG10 DENSITIES FOR PLOTTING 
!         xco2l(k,i)=log10(xco2(k,i))
!         xcol(k,i)=log10(xco(k,i))
!         xo2l(k,i)=log10(xo2(k,i))
!         xn2l(k,i)=log10(xn2(k,i))
!         xol(k,i)=log10(xo(k,i))
!
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

!  ** Overwriting update to test cooling/heating terms near the poles **
!     do i=lon0,lon1
!       do k=lev0,lev1
!       if (lat <= 14) then  !~56 degrees N/S
!               xo(k,i) = xo(k,i)*1.10
!       elseif (lat >= 59) then  !~56 degrees N/S
!               xo(k,i) = xo(k,i)*1.10
!       else
!               xo(k,i) = xo(k,i)
!       endif
!       enddo ! k=lev0,lev1
!     enddo ! i=lon0,lon1
!
! Fields OK: 11/21/06
!     call addfld('RJI',' ',' ',rji,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('QJO2I',' ',' ',qjo2i,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('XCO2L',' ',' ',xco2l,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('XCOL' ,' ',' ',xcol ,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('XN2L' ,' ',' ',xn2l ,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('XO2L' ,' ',' ',xo2l ,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('XOL'  ,' ',' ',xol ,'lev',k0,k1,'lon',i0,i1,lat)
!
! --------------------------------------------------------------------
      nk = lev1-lev0+1
! ---------------------------------------------------------------------
! Major sources and losses for O and CO (before NOX Chemistry Added)
! -- Checked        (10/29/07-11/02/07)
! -- Re-Checked     (11/16/07)  No OSRC or COSRC less than z = -12.0 
!                   since JCO2=JO2=0 plus ions are negligible
! -- Re-checked     (04/22/08)  No idn_latlon flags; OK from qrj and ions
! ---------------------------------------------------------------------
!
! Loop over subdomain at current latitude: for comp usage
      do i=lon0,lon1
        do k=lev0,lev1-1
!! -----------------------------------
!! Dayside and Nightside OSRC, COSRC, OLOSS, COLOSS
!! -----------------------------------
!! O sources : all [1,0] terms of fs-array
!!
!          osrc(k,i) =  rji(k,i) + 2.*qjo2i(k,i) + 2.*o2p(k,i)*ne(k,i)*
!     |      a13(k,i,lat) + co2p(k,i)*ne(k,i)*a14(k,i,lat)
!!
!! O losses : [1,1] and [2,1] terms of fs-array
!! 
!!         oloss(k,i) =  2.*xo(k,i)*xco2(k,i)*r2(k,i,lat)+
!!         oloss(k,i) =  2.*xo(k,i)*xco2(k,i)*r2+
!!    |       co2p(k,i)*(r19+r20) + xo2(k,i)*xco2(k,i)*r5(k,i,lat)
!!
!          oloss(k,i) =  2.*xo(k,i)*xco2(k,i)*r2+
!     |           co2p(k,i)*(r19+r20) + 
!     |           xo2(k,i)*xco2(k,i)*r5 +
!!       Low Altitude
!     |           xo3(k,i)*rk20(k,i,lat) +
!     |           xhcl(k,i)*r66(k,i,lat) +
!     |           xh2(k,i)*r104(k,i,lat) +
!     |           xoh(k,i)*rk28(k,i,lat) +
!     |           xho2(k,i)*r49(k,i,lat) +
!     |           xclco(k,i)*r102a +
!     |           xclco(k,i)*r102b +
!     |           xcl2(k,i)*r108(k,i,lat) +
!     |           xclo(k,i)*r44(k,i,lat)
!          oloss3(k,i) =  xo(k,i)*xco2(k,i)*r3(k,i,lat)
!          oloss2(k,i) = oloss(k,i)*xo(k,i)+oloss3(k,i)*xco(k,i)
!          otime(k,i) = 1./(oloss(k,i)+oloss3(k,i))/86400.
!          xopce(k,i) =  osrc(k,i)/oloss2(k,i)
!
! ---------------------------------------------------------------------
! Day/Night Split
         if (idn_lonlat(i,lat)==1)then
! -----------------------------------
! Dayside OSRC, OLOSS
! -----------------------------------
          osrc(k,i) =  rji(k,i) + 2.*qjo2i(k,i) + 2.*o2p(k,i)*ne(k,i)*
     |      a13(k,i,lat) + co2p(k,i)*ne(k,i)*a14(k,i,lat)
          oloss(k,i) =  2.*xo(k,i)*xco2(k,i)*r2(k,i,lat)+
     |           co2p(k,i)*(r19+r20) + 
     |           xo2(k,i)*xco2(k,i)*r5(k,i,lat) + 
!       Low Altitude
     |           xo3(k,i)*rk20(k,i,lat) +
     |           xhcl(k,i)*r66(k,i,lat) +
     |           xh2(k,i)*r104(k,i,lat) +
     |           xoh(k,i)*rk28(k,i,lat) +
     |           xho2(k,i)*r49(k,i,lat) +
     |           xclco(k,i)*r102a +
     |           xclco(k,i)*r102b +
     |           xcl2(k,i)*r108(k,i,lat) +
     |           xclo(k,i)*r44(k,i,lat)
!
! output dayside Tracers
!          olosDd(k,i) = xo(k,i)*xo3(k,i)*rk20(k,i,lat)
!          olosEd(k,i) = xo(k,i)*xhcl(k,i)*r66(k,i,lat)
!          olosFd(k,i) = xo(k,i)*xh2(k,i)*r104(k,i,lat)
!          olosGd(k,i) = xo(k,i)*xoh(k,i)*rk28(k,i,lat)
!          olosHd(k,i) = xo(k,i)*xho2(k,i)*r49(k,i,lat)
!          olosId(k,i) = xo(k,i)*xclco(k,i)*r102a
!          olosJd(k,i) = xo(k,i)*xclco(k,i)*r102b
!          olosKd(k,i) = xo(k,i)*xcl2(k,i)*r108(k,i,lat)
!          olosLd(k,i) = xo(k,i)*xclo(k,i)*r44(k,i,lat)
!          olosTrdy(k,i) = olosDd(k,i)+olosEd(k,i)+olosFd(k,i)+
!    |             olosGd(k,i)+olosHd(k,i)+olosId(k,i) +
!    |             olosJd(k,i)+olosKd(k,i)+olosLd(k,i)
         else
! -----------------------------------
! Nightside OSRC, OLOSS
! -----------------------------------
          osrc(k,i) = 0.0
          oloss(k,i) =  2.*xo(k,i)*xco2(k,i)*r2(k,i,lat)+
     |           co2p(k,i)*(r19+r20) + 
     |           xo2(k,i)*xco2(k,i)*r5(k,i,lat)+   
!       Low Altitude
     |           xo3(k,i)*rk20(k,i,lat) +
     |           xhcl(k,i)*r66(k,i,lat) +
     |           xh2(k,i)*r104(k,i,lat) +
     |           xoh(k,i)*rk28(k,i,lat) +
     |           xho2(k,i)*r49(k,i,lat) +
     |           xclco(k,i)*r102a +
     |           xclco(k,i)*r102b +
     |           xcl2(k,i)*r108(k,i,lat) +
     |           xclo(k,i)*r44(k,i,lat)
!
! output nightside Tracers
           olosDn(k,i) = xo(k,i)*xo3(k,i)*rk20(k,i,lat)
           olosEn(k,i) = xo(k,i)*xhcl(k,i)*r66(k,i,lat)
           olosFn(k,i) = xo(k,i)*xh2(k,i)*r104(k,i,lat)
           olosGn(k,i) = xo(k,i)*xoh(k,i)*rk28(k,i,lat)
           olosHn(k,i) = xo(k,i)*xho2(k,i)*r49(k,i,lat)
           olosIn(k,i) = xo(k,i)*xclco(k,i)*r102a
           olosJn(k,i) = xo(k,i)*xclco(k,i)*r102b
           olosKn(k,i) = xo(k,i)*xcl2(k,i)*r108(k,i,lat)
           olosLn(k,i) = xo(k,i)*xclo(k,i)*r44(k,i,lat)
           olosTrnt(k,i) = olosDn(k,i)+olosEn(k,i)+olosFn(k,i)+
     |             olosGn(k,i)+olosHn(k,i)+olosIn(k,i) +
     |             olosJn(k,i)+olosKn(k,i)+olosLn(k,i)
         endif
! -----------------------------------
         oloss3(k,i) =  xo(k,i)*xco2(k,i)*r3(k,i,lat)
         oloss2(k,i) = oloss(k,i)*xo(k,i)+oloss3(k,i)*xco(k,i)
         otime(k,i) = 1./(oloss(k,i)+oloss3(k,i))/86400.
         xopce(k,i) =  osrc(k,i)/oloss2(k,i)
!
!
         osrcA(k,i) = rji(k,i)
         osrcB(k,i) = 2.*qjo2i(k,i)
         osrcC(k,i) = 2.*o2p(k,i)*ne(k,i)*a13(k,i,lat)
         osrcD(k,i) = co2p(k,i)*ne(k,i)*a14(k,i,lat)
         osrcT(k,i) = osrcA(k,i)+osrcB(k,i)+osrcC(k,i)+osrcD(k,i)
!
         olossA(k,i) = xo(k,i)*2.*xo(k,i)*xco2(k,i)*r2(k,i,lat)
         olossB(k,i) = xo(k,i)*co2p(k,i)*(r19+r20)
         olossC(k,i) = xo(k,i)*xo2(k,i)*xco2(k,i)*r5(k,i,lat)
!
         onet(k,i) = osrcT(k,i) - oloss2(k,i)
! 
! ---------------------------------------------------------------------
! Day/Night Split
         if (idn_lonlat(i,lat)==1)then
! CO Day sources : [2,0] and [2,1] terms are separated for fs-array
!  ------------
!         cosrc(k,i) =  2.*rji(k,i) + co2p(k,i)*ne(k,i)*a14(k,i,lat) +
!    |     op(k,i)*xco2(k,i)*r18
!  ------------
          cosrc(k,i) =  rji(k,i) + co2p(k,i)*ne(k,i)*a14(k,i,lat) +
     |     op(k,i)*xco2(k,i)*r18 + co2p(k,i)*xo(k,i)*r20 +
!       Low Altitude
     |           xo(k,i)*xclco(k,i)*r144 +
     |           xh(k,i)*xclco(k,i)*r146 +
     |           xoh(k,i)*xclco(k,i)*r147 +
     |           xco2(k,i)*xclco(k,i)*r160(k,i,lat)
          cosrcA(k,i) =  rji(k,i)
          cosrcB(k,i) =  co2p(k,i)*ne(k,i)*a14(k,i,lat)
          cosrcC(k,i) =  op(k,i)*xco2(k,i)*r18
          cosrcD(k,i) =  co2p(k,i)*xo(k,i)*r20
          cosrcE(k,i) =  xo(k,i)*xclco(k,i)*r144 
          cosrcF(k,i) =  xh(k,i)*xclco(k,i)*r146
          cosrcG(k,i) =  xoh(k,i)*xclco(k,i)*r147
          cosrcH(k,i) =  xco2(k,i)*xclco(k,i)*r160(k,i,lat)
          cosrc2(k,i) =  co2p(k,i)*xo(k,i)*r20
          cosrc3(k,i) =  co2p(k,i)*r20
! CO Day losses :
! -------------
          coloss(k,i)  = xo(k,i)*xco2(k,i)*r3(k,i,lat) + 
!       Low Altitude
     |           xcl(k,i)*xco2(k,i)*r99(k,i,lat) +
     |           xoh(k,i)*r25 +
     |           xclo(k,i)*r43(k,i,lat) +
     |           xcl(k,i)*xn2(k,i)*r111(k,i,lat) +
     |           xo(k,i)*xo(k,i)*r163(k,i,lat) +
     |           xo(k,i)*xco(k,i)*r162(k,i,lat)
          colossA(k,i) = xco(k,i)*xo(k,i)*xco2(k,i)*r3(k,i,lat)
          colossB(k,i) = xco(k,i)*xcl(k,i)*xco2(k,i)*r99(k,i,lat)
          colossC(k,i) = xco(k,i)*xoh(k,i)*r25
          colossD(k,i) = xco(k,i)*xclo(k,i)*r43(k,i,lat)
          colossE(k,i) = xco(k,i)*xcl(k,i)*xn2(k,i)*r111(k,i,lat)
          colossF(k,i) = xco(k,i)*xo(k,i)*xo(k,i)*r163(k,i,lat)
          colossG(k,i) = xco(k,i)*xco(k,i)*xo(k,i)*r162(k,i,lat)
         else
! CO Nite sources : [2,0] and [2,1] terms are separated for fs-array
! --------------
          cosrc(k,i) =  co2p(k,i)*ne(k,i)*a14(k,i,lat) +
     |     op(k,i)*xco2(k,i)*r18 + co2p(k,i)*xo(k,i)*r20 +
!       Low Altitude
     |           xo(k,i)*xclco(k,i)*r144 +
     |           xh(k,i)*xclco(k,i)*r146 +
     |           xoh(k,i)*xclco(k,i)*r147 +
     |           xco2(k,i)*xclco(k,i)*r160(k,i,lat)
          cosrcB(k,i) =  co2p(k,i)*ne(k,i)*a14(k,i,lat)
          cosrcC(k,i) =  op(k,i)*xco2(k,i)*r18
          cosrcD(k,i) =  co2p(k,i)*xo(k,i)*r20
          cosrcE(k,i) =  xo(k,i)*xclco(k,i)*r144 
          cosrcF(k,i) =  xh(k,i)*xclco(k,i)*r146
          cosrcG(k,i) =  xoh(k,i)*xclco(k,i)*r147
          cosrcH(k,i) =  xco2(k,i)*xclco(k,i)*r160(k,i,lat)
          cosrc2(k,i) =  co2p(k,i)*xo(k,i)*r20
          cosrc3(k,i) =  co2p(k,i)*r20
! CO Nite losses :
! --------------
          coloss(k,i)  = xo(k,i)*xco2(k,i)*r3(k,i,lat) +
!       Low Altitude
     |           xcl(k,i)*xco2(k,i)*r99(k,i,lat) +
     |           xoh(k,i)*r25 +
     |           xclo(k,i)*r43(k,i,lat) +
     |           xcl(k,i)*xn2(k,i)*r111(k,i,lat) +
     |           xo(k,i)*xo(k,i)*r163(k,i,lat) +
     |           xo(k,i)*xco(k,i)*r162(k,i,lat)
          colossA(k,i) = xco(k,i)*xo(k,i)*xco2(k,i)*r3(k,i,lat)
          colossB(k,i) = xco(k,i)*xcl(k,i)*xco2(k,i)*r99(k,i,lat)
          colossC(k,i) = xco(k,i)*xoh(k,i)*r25
          colossD(k,i) = xco(k,i)*xclo(k,i)*r43(k,i,lat)
          colossE(k,i) = xco(k,i)*xcl(k,i)*xn2(k,i)*r111(k,i,lat)
          colossF(k,i) = xco(k,i)*xo(k,i)*xo(k,i)*r163(k,i,lat)
          colossG(k,i) = xco(k,i)*xco(k,i)*xo(k,i)*r162(k,i,lat)
          endif
!  Prints/plots of CO source and terms
          coloss2(k,i) = coloss(k,i)*xco(k,i)
          cotime(k,i) = (1./coloss(k,i))/86400.
          xcopce(k,i) =  (cosrc(k,i)+cosrc2(k,i))/coloss(k,i)
!
! ---------------------------------------------------------------------
! Check This Carefully Now!  January 2016!
! Matrix coefficients for O-CO-CO2 solution: from vtgcm2 code
! Cross terms are neglected here (minor). How to add back?
! Sources are translated from #/cm3.sec to mass mixing units
! Losses are already compatible mass mixing units (1/sec)
! -- coloss(k,i) could be moved to [2,1] term and scaled by 28./16.
!
! --- vtgcm2 methodology
!
!         fs(i,k,1,1,lat) = -oloss(k,i)
!         fs(i,k,1,2,lat) = -oloss3(k,i)*rmass_o1*rmassinv_co
!         fs(i,k,2,1,lat) = cosrc3(k,i)*rmass_co*rmassinv_o1 
!         fs(i,k,2,2,lat) = -coloss(k,i)
!         fs(i,k,1,0,lat) = osrc(k,i)*rmass_o1/xnmbari(k,i)
!         fs(i,k,2,0,lat) = cosrc(k,i)*rmass_co/xnmbari(k,i)
!
! --- vtgcm1.2  methodology
          fs(i,k,1,1,lat) = -oloss(k,i)
          fs(i,k,1,2,lat) = 0.0
          fs(i,k,2,1,lat) = 0.0
          fs(i,k,2,2,lat) = -coloss(k,i)
          fs(i,k,1,0,lat) = osrc(k,i)*rmass_o1/xnmbari(k,i)
          fs(i,k,2,0,lat) = cosrc(k,i)*rmass_co/xnmbari(k,i)
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!
! ---------------------------------------------------------------------
      do i = lon0,lon1
        do k = lev0,lev1-1
          slab(k,i) = gask*tni(k,i)/(barm(k,i)*grav)
        enddo
      enddo
! *O Total production INTEGRATED COLUMN RATE (#/cm2/sec)
      do i = lon0,lon1
      if (dz == 0.5) then
! Below 90 km is Zeroed out (because of O-lbc)
        do k = 1,11
           fldht(k) = 0.
        enddo
        do k = 11,lev1-1
           fldht(k) = slab(k,i)*osrcT(k,i)
        enddo
      elseif (dz == 0.25) then
! Below 90 km is Zeroed out (because of O-lbc)
        do k = 1,21
           fldht(k) = 0.
        enddo
        do k = 21,lev1-1
           fldht(k) = slab(k,i)*osrcT(k,i)
        enddo
      endif
        call integral(colout,fldht,1,dz,lev1-1)
        do k = lev0,lev1-1
           oprodcol(k,i) = colout(k)
        enddo
      enddo
! *O Net production INTEGRATED COLUMN RATE (#/cm2/sec)
      do i = lon0,lon1
      if (dz == 0.5) then
! Below 90 km is Zeroed out (because of O-lbc)
        do k = 1,11
           fldht(k) = 0.
        enddo
        do k = 11,lev1-1
           fldht(k) = slab(k,i)*onet(k,i)
        enddo
      elseif (dz == 0.25) then
! Below 90 km is Zeroed out (because of O-lbc)
        do k = 1,21
           fldht(k) = 0.
        enddo
        do k = 21,lev1-1
           fldht(k) = slab(k,i)*onet(k,i)
        enddo
      endif
        call integral(colout,fldht,1,dz,lev1-1)
        do k = lev0,lev1-1
           onetprodcol(k,i) = colout(k)
        enddo
      enddo
! *O destruction INTEGRATED COLUMN RATE (#/cm2/sec)
      do i = lon0,lon1
      if (dz == 0.5) then
! Below 90 km is Zeroed out (because of O-lbc)
        do k = 1,11
           fldht(k) = 0.
        enddo
        do k = 11,lev1-1
           fldht(k) = slab(k,i)*(oloss2(k,i))
        enddo
      elseif (dz == 0.25) then
! Below 90 km is Zeroed out (because of O-lbc)
        do k = 1,21
           fldht(k) = 0.
        enddo
        do k = 21,lev1-1
           fldht(k) = slab(k,i)*(oloss2(k,i))
        enddo
      endif
        call integral(colout,fldht,1,dz,lev1-1)
        do k = lev0,lev1-1
           odestcol(k,i) = colout(k)
        enddo
      enddo

! Loop over subdomain at current latitude: diagnostics
      do i=lon0,lon1
        do k=lev0,lev1-1
! --------------------------------------------------------
        if (osrcA(k,i) .GT. 1.0) THEN
           osrcA(k,i)= log10(osrcA(k,i))
        else
           osrcA(k,i) = 0.
        endif
        if (osrcB(k,i) .GT. 1.0) THEN
           osrcB(k,i)= log10(osrcB(k,i))
        else
           osrcB(k,i) = 0.
        endif
        if (osrcC(k,i) .GT. 1.0) THEN
           osrcC(k,i)= log10(osrcC(k,i))
        else
           osrcC(k,i) = 0.
        endif
        if (osrcD(k,i) .GT. 1.0) THEN
           osrcD(k,i)= log10(osrcD(k,i))
        else
           osrcD(k,i) = 0.
        endif
        if (osrcT(k,i) .GT. 1.0) THEN
           osrcT(k,i)= log10(osrcT(k,i))
        else
           osrcT(k,i) = 0.
        endif
! --------------------------------------------------------
        if (olossA(k,i) .GT. 1.0) THEN
           olossA(k,i)= log10(olossA(k,i))
        else
           olossA(k,i) = 0.
        endif
        if (olossB(k,i) .GT. 1.0) THEN
           olossB(k,i)= log10(olossB(k,i))
        else
           olossB(k,i) = 0.
        endif
        if (olossC(k,i) .GT. 1.0) THEN
           olossC(k,i)= log10(olossC(k,i))
        else
           olossC(k,i) = 0.
        endif
        if (olosDn(k,i) .GT. 1.0) THEN
           olosDn(k,i)= log10(olosDn(k,i))
        else
           olosDn(k,i) = 0.
        endif
        if (olosEn(k,i) .GT. 1.0) THEN
           olosEn(k,i)= log10(olosEn(k,i))
        else
           olosEn(k,i) = 0.
        endif
        if (olosFn(k,i) .GT. 1.0) THEN
           olosFn(k,i)= log10(olosFn(k,i))
        else
           olosFn(k,i) = 0.
        endif
        if (olosGn(k,i) .GT. 1.0) THEN
           olosGn(k,i)= log10(olosGn(k,i))
        else
           olosGn(k,i) = 0.
        endif
        if (olosHn(k,i) .GT. 1.0) THEN
           olosHn(k,i)= log10(olosHn(k,i))
        else
           olosHn(k,i) = 0.
        endif
        if (olosIn(k,i) .GT. 1.0) THEN
           olosIn(k,i)= log10(olosIn(k,i))
        else
           olosIn(k,i) = 0.
        endif
        if (olosJn(k,i) .GT. 1.0) THEN
           olosJn(k,i)= log10(olosJn(k,i))
        else
           olosJn(k,i) = 0.
        endif
        if (olosKn(k,i) .GT. 1.0) THEN
           olosKn(k,i)= log10(olosKn(k,i))
        else
           olosKn(k,i) = 0.
        endif
        if (olosLn(k,i) .GT. 1.0) THEN
           olosLn(k,i)= log10(olosLn(k,i))
        else
           olosLn(k,i) = 0.
        endif
        if (olosTrnt(k,i) .GT. 1.0) THEN
           olosTrnt(k,i)= log10(olosTrnt(k,i))
        else
           olosTrnt(k,i) = 0.
        endif
! --------------------------------------------------------
        if (cosrcA(k,i) .GT. 1.0) THEN
           cosrcA(k,i)= log10(cosrcA(k,i))
        else
           cosrcA(k,i) = 0.
        endif
        if (cosrcB(k,i) .GT. 1.0) THEN
           cosrcB(k,i)= log10(cosrcB(k,i))
        else
           cosrcB(k,i) = 0.
        endif
        if (cosrcC(k,i) .GT. 1.0) THEN
           cosrcC(k,i)= log10(cosrcC(k,i))
        else
           cosrcC(k,i) = 0.
        endif
        if (cosrcD(k,i) .GT. 1.0) THEN
           cosrcD(k,i)= log10(cosrcD(k,i))
        else
           cosrcD(k,i) = 0.
        endif
        if (cosrcE(k,i) .GT. 1.0) THEN
           cosrcE(k,i)= log10(cosrcE(k,i))
        else
           cosrcE(k,i) = 0.
        endif
        if (cosrcF(k,i) .GT. 1.0) THEN
           cosrcF(k,i)= log10(cosrcF(k,i))
        else
           cosrcF(k,i) = 0.
        endif
        if (cosrcG(k,i) .GT. 1.0) THEN
           cosrcG(k,i)= log10(cosrcG(k,i))
        else
           cosrcG(k,i) = 0.
        endif
        if (cosrcH(k,i) .GT. 1.0) THEN
           cosrcH(k,i)= log10(cosrcH(k,i))
        else
           cosrcH(k,i) = 0.
        endif
! --------------------------------------------------------
        if (colossA(k,i) .GT. 1.0) THEN
           colossA(k,i)= log10(colossA(k,i))
        else
           colossA(k,i) = 0.
        endif
        if (colossB(k,i) .GT. 1.0) THEN
           colossB(k,i)= log10(colossB(k,i))
        else
           colossB(k,i) = 0.
        endif
        if (colossC(k,i) .GT. 1.0) THEN
           colossC(k,i)= log10(colossC(k,i))
        else
           colossC(k,i) = 0.
        endif
        if (colossD(k,i) .GT. 1.0) THEN
           colossD(k,i)= log10(colossD(k,i))
        else
           colossD(k,i) = 0.
        endif
        if (colossE(k,i) .GT. 1.0) THEN
           colossE(k,i)= log10(colossE(k,i))
        else
           colossE(k,i) = 0.
        endif
        if (colossF(k,i) .GT. 1.0) THEN
           colossF(k,i)= log10(colossF(k,i))
        else
           colossF(k,i) = 0.
        endif
        if (colossG(k,i) .GT. 1.0) THEN
           colossG(k,i)= log10(colossG(k,i))
        else
           colossG(k,i) = 0.
        endif
! --------------------------------------------------------
!       if (osrc(k,i) .GT. 1.0) THEN
!          osrc(k,i)= log10(osrc(k,i))
!       else
!          osrc(k,i) = 0.
!       endif
!       if (oloss2(k,i) .GT. 1.0) THEN
!          oloss2(k,i)= log10(oloss2(k,i))
!       else
!          oloss2(k,i) = 0.
!       endif
!       if (otime(k,i) .GT. 0.0) THEN
!          otimel(k,i)= log10(otime(k,i))
!       else
!          otimel(k,i) = 0.
!       endif
!       if (cosrc(k,i) .GT. 1.0) THEN
!          cosrc(k,i)= log10(cosrc(k,i))
!       else
!          cosrc(k,i) = 0.
!       endif
!       if (coloss2(k,i) .GT. 1.0) THEN
!          coloss2(k,i)= log10(coloss2(k,i))
!       else
!           coloss2(k,i) = 0.
!       endif
!       if (cotime(k,i) .GT. 0.0) THEN
!          cotimel(k,i)= log10(cotime(k,i))
!       else
!          cotimel(k,i) = 0.
!       endif
!
! --------------------------------------------------------
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

!     do k=lev0,lev1-1
!       write(6,"('comp_oco: lat=',i3,' k=',i3)") lat,k
!       write(6,"('fs(:,k,1,1,lat)=',/,(6e12.4))") fs(:,k,1,1,lat)
!       write(6,"('fs(:,k,1,2,lat)=',/,(6e12.4))") fs(:,k,1,2,lat)
!       write(6,"('fs(:,k,2,1,lat)=',/,(6e12.4))") fs(:,k,2,1,lat)
!       write(6,"('fs(:,k,2,2,lat)=',/,(6e12.4))") fs(:,k,2,2,lat)
!       write(6,"('fs(:,k,1,0,lat)=',/,(6e12.4))") fs(:,k,1,0,lat)
!       write(6,"('fs(:,k,2,0,lat)=',/,(6e12.4))") fs(:,k,2,0,lat)
!
!     enddo ! k=lev0,lev1-1

!  ----------------------------------------------------------------------
!     call addfld('OSRC',' ',' ',osrc,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('OLOSS',' ',' ',oloss2,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('OLIFE',' ',' ',otimel,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('OPRODCOL',' ',' ',oprodcol,'lev',k0,k1,'lon',i0,i1,
!    |  lat)
!     call addfld('ONETPRODCOL',' ',' ',onetprodcol,'lev',k0,k1,'lon',
!    |  i0,i1,lat)
!     call addfld('ODESTCOL',' ',' ',odestcol,'lev',k0,k1,'lon',i0,i1,
!    |  lat)
!     call addfld('TAUO',' ',' ',otime,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('COSRC',' ',' ',cosrc,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('COLOSS',' ',' ',coloss2,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('COLIFE',' ',' ',cotimel,'lev',k0,k1,'lon',i0,i1,lat)
!  ----------------------------------------------------------------------
!     call addfld('XCLD',' ',' ',xcld,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('XCLN',' ',' ',xcln,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('XCL2D',' ',' ',xcl2d,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('XCL2N',' ',' ',xcl2n,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('XCLOD',' ',' ',xclod,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('XCLON',' ',' ',xclon,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('XCLCON',' ',' ',xclcon,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('XCLCOD',' ',' ',xclcod,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('XOHN',' ',' ',xohn,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('XOHD',' ',' ',xohd,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('XHO2N',' ',' ',xho2n,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('XHO2D',' ',' ',xho2d,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('XO3N',' ',' ',xo3n,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('XO3D',' ',' ',xo3d,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('XHCLN',' ',' ',xhcln,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('XHCLD',' ',' ',xhcld,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('XH2N',' ',' ',xh2n,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('XH2D',' ',' ',xh2d,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('XTOT',' ',' ',xtot,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('MeanMass',' ',' ',barm,'lev',k0,k1,'lon',i0,i1,lat)
!  ----------------------------------------------------------------------
!     call addfld('RJIS',' ',' ',rjis,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('QJO2IS',' ',' ',qjo2is,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('A13',' ',' ',a13,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('A14',' ',' ',a14,'lev',k0,k1,'lon',i0,i1,lat)
!  ----------------------------------------------------------------------
      call addfld('COSRCA',' ',' ',cosrcA,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('COSRCB',' ',' ',cosrcB,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('COSRCC',' ',' ',cosrcC,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('COSRCD',' ',' ',cosrcD,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('COSRCE',' ',' ',cosrcE,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('COSRCF',' ',' ',cosrcF,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('COSRCG',' ',' ',cosrcG,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('COSRCH',' ',' ',cosrcH,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('COLOSSA',' ',' ',colossA,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('COLOSSB',' ',' ',colossB,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('COLOSSC',' ',' ',colossC,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('COLOSSD',' ',' ',colossD,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('COLOSSE',' ',' ',colossE,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('COLOSSF',' ',' ',colossF,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('COLOSSG',' ',' ',colossG,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('OSRCA',' ',' ',osrcA,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('OSRCB',' ',' ',osrcB,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('OSRCC',' ',' ',osrcC,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('OSRCD',' ',' ',osrcD,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('OSRCT',' ',' ',osrcT,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('OLOSSA',' ',' ',olossA,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('OLOSSB',' ',' ',olossB,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('OLOSSC',' ',' ',olossC,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('OLOSDN',' ',' ',olosDn,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('OLOSEN',' ',' ',olosEn,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('OLOSFN',' ',' ',olosFn,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('OLOSGN',' ',' ',olosGn,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('OLOSHN',' ',' ',olosHn,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('OLOSIN',' ',' ',olosIn,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('OLOSJN',' ',' ',olosJn,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('OLOSKN',' ',' ',olosKn,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('OLOSLN',' ',' ',olosLn,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('OLOSTN',' ',' ',olosTrnt,'lev',
     | k0,k1,'lon',i0,i1,lat)
!  ----------------------------------------------------------------------
!     call addfld('OPCE',' ',' ',xopce,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('COPCE',' ',' ',xcopce,'lev',k0,k1,'lon',i0,i1,lat)
!  ----------------------------------------------------------------------
      end subroutine comp_oco
! ----------------------------------------------------------------------
