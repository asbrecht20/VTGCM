!
      subroutine comp_oh(tn,barm,o1,co,co2,
     |  n2,o2,op,co2p,o2p,ne,
     |  lev0,lev1,lon0,lon1,lat)
!
      use params_module,only: dz
      use cons_module,only: rmassinv_o1,rmassinv_co,rmassinv_n2,p0,
     |  expz,boltz,rmassinv_co2,rmassinv_o2,rmass_o1,rmass_co,
     |  gask,grav
      use o2_module,only: integral
      use qrj_module,only: rj,qjo2,
     |  mk_xnmbari
      use fields_module,only: tlbc
      use chemrates_module,only: r49,r41b,r66,r104,r54,
     |  r51,rk28,r53,r105,r2,rk23,r63,
     |  k50,k51,k52,r5,rk20,r53,rk23,r54,k53,k54,k55,k56,
!    ********** fractional yield  (ND)************
     |  v1,v2,v3,v4,v5,v6,v7,v8,v9,
     |  v21,v22,v23,v24,v25,v26,v27,v28,v29,
!    ********** cm3/sec ************
     |  kv98,kv97,kv96,kv95,kv94,kv93,kv92,kv91,kv9g,
     |  kv87,kv86,kv85,kv84,kv83,kv82,kv81,kv8g,
     |  kv76,kv75,kv74,kv73,kv72,kv71,kv7g,
     |  kv65,kv64,kv63,kv62,kv61,kv6g,
     |  kv54,kv53,kv52,kv51,kv5g,
     |  kv43,kv42,kv41,kv4g,
     |  kv32,kv31,kv3g,
     |  kv21,kv2g,
     |  kv1g,
!    ********** 1/sec ************
     |  sv98,sv97,sv96,sv95,sv94,sv93,sv92,sv91,sv9g,
     |  sv87,sv86,sv85,sv84,sv83,sv82,sv81,sv8g,
     |  sv76,sv75,sv74,sv73,sv72,sv71,sv7g,
     |  sv65,sv64,sv63,sv62,sv61,sv6g,
     |  sv54,sv53,sv52,sv51,sv5g,
     |  sv43,sv42,sv41,sv4g,
     |  sv32,sv31,sv3g,
     |  sv21,sv2g,
     |  sv1g
!
      use chapman_module,only:idn_lonlat
      use addfld_module,only: addfld
      implicit none

      real :: FOHIR, ksum, kssum
      real :: BETA =  1.0   !  Yield for OH IR nightglow near the peak

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

! Local:
      integer :: k,i,nk,k0,k1,i0,i1
!
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  tni,            ! TN at interfaces
     |  rji,            ! JCO2 at interfaces (qrj module)
     |  qjo2i,          ! JO2 at interfaces (qrj module)
     |  xnmbari,        ! xnmbar at interfaces (qrj module)
     |  xnmbar,         ! xnmbar at midlevels (qrj module)
     |  ohsrc,          ! OH production terms  (#/cm3/sec)
     |  ohloss,         ! OH loss terms  (1/sec)
     |  ohlossr,        ! OH loss terms  (#/cm3/sec)
     |  o3src,          ! O3 production terms  (#/cm3/sec)
     |  o3loss,         ! O3 loss terms  (1/sec)
     |  o3lossr,        ! O3 loss terms  (#/cm3/sec)
     |  ohden,          ! calculated OH number density
     |  o3den,          ! calculated O3 number density
     |  ohdenvmr,       ! calculated OH VMR
     |  o3denvmr,       ! calculated O3 VMR
     |  xovmr,          ! calculated O VMR
     |  xhvmr,          ! calculated H VMR
     |  xho2vmr,        ! calculated HO2 VMR
     |  xo,             ! O1 number density (#/cm3)
     |  xopce,          ! O1 number density (#/cm3): PCE secondary
     |  xco,            ! CO number density (#/cm3)
     |  xcopce,         ! CO number density (#/cm3): PCE secondary
     |  xco2,           ! CO2 number density (#/cm3)
     |  xn2,            ! N2 number density (#/cm3)
     |  xo2,            ! O2 number density (#/cm3)
     |  xcl,            ! CL number density (#/cm3)
     |  xclo,           ! CLO number density (#/cm3)
     |  xclco,          ! CLCO number density (#/cm3)
     |  xcl2,           ! CL2 number desity (#/cm3)
     |  xoh,            ! OH number density (#/cm3): Straight P/L for Gnd state
     |  xh,             ! H number density (#/cm3)
     |  xho2,           ! HO2 number density (#/cm3)
     |  xo3,            ! O3 number density (#/cm3); !-D KINETICS
     |  xh2,            ! h2 number density (#/cm3)
     |  xhcl,           ! hcl number density (#/cm3)
     |  xtot,	        ! Total number density (#/cm3)
     |  xol,            ! LOG10-O1 number density (#/cm3)
     |  xcol,           ! LOG10-CO number density (#/cm3)
     |  xco2l,          ! LOG10-CO2 number density (#/cm3)
     |  xn2l,           ! LOG10-N2 number density (#/cm3)
     |  xo2l,           ! LOG10-O2 number density (#/cm3)
     |  xohv9,          ! Production OH vib level (#/cm3/s): hv9
     |  xohv8,
     |  xohv7,
     |  xohv6,
     |  xohv5,
     |  xohv4,
     |  xohv3,
     |  xohv2,
     |  xohv1,
     |  xohvg,
     |  xohv9den,       ! OH number density (#/cm3): hv9
     |  xohv8den,       
     |  xohv7den,      
     |  xohv6den,     
     |  xohv5den,    
     |  xohv4den,       
     |  xohv3den,      
     |  xohv2den,     
     |  xohv1den,       
     |  xoh2,           ! OH number density (#/cm3): Cascading xohvg source + full loss for gnd state
     |  slab,           ! slab width (cm)
     |  ohirver,        ! ohir volume emission (photons/cm3/sec)
     |  ohirverl,       ! ohir volume emission (log10 photons/cm3/sec)
     |  ohirint         ! ohir integrated intensity (photons/cm2/sec in K Rayleighs)
      real,dimension(lev0:lev1) :: colout,fldht

!
! Specificed specifically for OH airglow
      real,dimension(lev0:lev1) ::
     |  cl,            ! CL (vmr) 
     |  clo,           ! CLO (vmr) 
     |  clco,          ! CLCO (vmr) 
     |  cl2,           ! CL2 (vmr) 
     |  oh,            ! OH (vmr) 
     |  h,             ! H (vmr) 
     |  ho2,           ! HO2 (vmr) 
     |  o3,            ! O3 (vmr) 
     |  hcl,           ! hcl (vmr)
     |  h2             ! h2 (vmr)
!
!
! For addfld calls:
      k0=lev0
      k1=lev1
      i0=lon0
      i1=lon1
      nk = lev1-1
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
! Starter fields for O3 adn OH (iteration)
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
! --------------------------------------------------------------------
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
          rji(k,i)  = 0.5*(rj(k,i,lat)+rj(k+1,i,lat))
          qjo2i(k,i) = 0.5*(qjo2(k,i,lat)+qjo2(k+1,i,lat))
!
!   LOG10 DENSITIES FOR PLOTTING
!         xco2l(k,i)=log10(xco2(k,i))
!         xcol(k,i)=log10(xco(k,i))
!         xo2l(k,i)=log10(xo2(k,i))
!         xn2l(k,i)=log10(xn2(k,i))
!         xol(k,i)=log10(xo(k,i))
!
! Creating total number density
!         xtot(k,i) = xo(k,i)+xco(k,i)+xco2(k,i)+xo2(k,i)+xn2(k,i)
          xtot(k,i) = xo(k,i)+xco(k,i)+xco2(k,i)+xn2(k,i)

!  Converting vmr to number density
          xcl(k,i)   = xtot(k,i)*cl(k)
          xclo(k,i)  = xtot(k,i)*clo(k)
          xclco(k,i) = xtot(k,i)*clco(k)
          xcl2(k,i)  = xtot(k,i)*cl2(k)
          xoh(k,i)   = xtot(k,i)*oh(k)
          xh(k,i)   = xtot(k,i)*h(k)
          xho2(k,i)  = xtot(k,i)*ho2(k)
          xo3(k,i) = xtot(k,i)*o3(k)
          xh2(k,i)  = xtot(k,i)*h2(k)
          xhcl(k,i) = xtot(k,i)*hcl(k)
!
!  Initializing OHvib level densities (#/cm3-sec) for each pass through LAT loop
!  Initializing OHIRVER (PH/cm3.sec) for each pass through LAT loop
!
          xohvg(k,i) = 0.0
          xohv1(k,i) = 0.0
          xohv2(k,i) = 0.0
          xohv3(k,i) = 0.0
          xohv4(k,i) = 0.0
          xohv5(k,i) = 0.0
          xohv6(k,i) = 0.0
          xohv7(k,i) = 0.0
          xohv8(k,i) = 0.0
          xohv9(k,i) = 0.0
          ohirver(k,i) = 0.0
!
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!
! -----------------------------------------------------------------------
! Loop over subdomain at current latitude: for comp usage
      do i=lon0,lon1
        do k=lev0,lev1-1
          if (idn_lonlat(i,lat)==1) then
            ohsrc(k,i) = 0.
            ohloss(k,i) = 0.
            ohden(k,i) = 0.
            xoh(k,i) = 0.0
            ohdenvmr(k,i) = 0.
            o3src(k,i) = 0. 
            o3loss(k,i) = 0. 
            o3den(k,i) = 0.
            o3denvmr(k,i) = 0.
            xohv9den(k,i) = 0.0
            xohv8den(k,i) = 0.0
            xohv7den(k,i) = 0.0
            xohv6den(k,i) = 0.0
            xohv5den(k,i) = 0.0
            xohv4den(k,i) = 0.0
            xohv3den(k,i) = 0.0
            xohv2den(k,i) = 0.0
            xohv1den(k,i) = 0.0
            xoh2(k,i) = 0.0
          else
! -----------------------------------
! Both or Nightside O3SRC, O3LOSS
! -----------------------------------
! O3 sources:
            o3src(k,i) = 
     |       xo(k,i)*xo(k,i)*xo2(k,i)*k50(k,i,lat) +
     |       xo(k,i)*xo2(k,i)*xo2(k,i)*k51(k,i,lat) +
     |       xo(k,i)*xo2(k,i)*xco(k,i)*k52(k,i,lat) +
     |       xo(k,i)*xo2(k,i)*xn2(k,i)*r63(k,i,lat) +
     |       xo(k,i)*xo2(k,i)*xco2(k,i)*r5(k,i,lat)
!
! O3 losses :
            o3loss(k,i) = 
     |       xo(k,i)*rk20(k,i,lat) +
     |       xoh(k,i)*r53(k,i,lat) +
     |       xh(k,i)*rk23(k,i,lat) +
     |       xho2(k,i)*r54(k,i,lat) +
     |       xcl(k,i)*k54(k,i,lat) +
     |       xclo(k,i)*k55 +
     |       xclo(k,i)*k56(k,i,lat)
!
! O3 Denisty (#/cm3): Simulate Self-consistently
!           o3den(k,i) = o3src(k,i)/o3loss(k,i)
! O3 Denisty (#/cm3): Simulate Self-consistently  (but scaled by 0.5 for sensitivity study)
            o3den(k,i) = 0.5*o3src(k,i)/o3loss(k,i)
! O3 Denisty (#/cm3): from 1-D KINETICS
!           o3den(k,i) = xo3(k,i)
!Converting from loss frequency to loss rate (#/cm3/sec)
!           o3lossr(k,i) =  o3loss(k,i)*o3den(k,i)
!Converting from number density to VMR 
            o3denvmr(k,i) = o3den(k,i)/xtot(k,i)
            xovmr(k,i) = xo(k,i)/xtot(k,i)
            xhvmr(k,i) = xh(k,i)/xtot(k,i)
            xho2vmr(k,i) = xho2(k,i)/xtot(k,i)
!
! -----------------------------------
! Both or Nightside OHSRC, OHLOSS
! -----------------------------------
! OH sources :
!
            ohsrc(k,i) = 
     |       xh(k,i)*o3den(k,i)*rk23(k,i,lat) +
     |       xho2(k,i)*xo(k,i)*r49(k,i,lat) +
     |       xho2(k,i)*xcl(k,i)*r41b(k,i,lat) +
     |       xhcl(k,i)*xo(k,i)*r66(k,i,lat) +
     |       xh2(k,i)*xo(k,i)*r104(k,i,lat) +
     |       xho2(k,i)*o3den(k,i)*r54(k,i,lat) 
!
! OH losses :
!
            ohloss(k,i) = 
     |       xco(k,i)*r51 +
     |       xo(k,i)*rk28(k,i,lat) +
     |       o3den(k,i)*r53(k,i,lat)
!    |       xoh(k,i)*r105(k,i,lat)
!
! OH Density(number density) :
            ohden(k,i) = ohsrc(k,i)/ohloss(k,i)
! Converting from loss frequency to loss rate (#/cm3/sec)
            ohlossr(k,i) =  ohloss(k,i)*ohden(k,i)
! Converting from number density to VMR 
            ohdenvmr(k,i) = ohden(k,i)/xtot(k,i)
!     
! ----------------------------------------------------------------------
!  ** Populating equations: Nascent Population Production (#/cm3.sec)
! ----------------------------------------------------------------------
!
!  Both channels:  HO2 + O and H + O3
!         xohvg(k,i) = ohsrc(k,i)
!         xohv1(k,i) = xho2(k,i)*xo(k,i)*v1  + xh(k,i)*o3den(k,i)*v21
!         xohv2(k,i) = xho2(k,i)*xo(k,i)*v2  + xh(k,i)*o3den(k,i)*v22
!         xohv3(k,i) = xho2(k,i)*xo(k,i)*v3(k,i,lat)  + 
!    |                 xh(k,i)*o3den(k,i)*v23
!         xohv4(k,i) = xho2(k,i)*xo(k,i)*v4(k,i,lat)  + 
!    |                 xh(k,i)*o3den(k,i)*v24
!         xohv5(k,i) = xho2(k,i)*xo(k,i)*v5(k,i,lat) +
!    |                 xh(k,i)*o3den(k,i)*v25(k,i,lat)
!         xohv6(k,i) = xho2(k,i)*xo(k,i)*v6(k,i,lat) +
!    |                 xh(k,i)*o3den(k,i)*v26(k,i,lat)
!         xohv7(k,i) = xho2(k,i)*xo(k,i)*v7 +
!    |                 xh(k,i)*o3den(k,i)*v27(k,i,lat)
!         xohv8(k,i) = xho2(k,i)*xo(k,i)*v8 +
!    |                 xh(k,i)*o3den(k,i)*v28(k,i,lat)
!         xohv9(k,i) = xho2(k,i)*xo(k,i)*v9 +
!    |                 xh(k,i)*o3den(k,i)*v29(k,i,lat)
! ----------------------------------------------------------------------
!  Only one channel:   H + O3
          xohvg(k,i) = ohsrc(k,i)
          xohv1(k,i) = xh(k,i)*o3den(k,i)*v21
          xohv2(k,i) = xh(k,i)*o3den(k,i)*v22
          xohv3(k,i) = xh(k,i)*o3den(k,i)*v23 
          xohv4(k,i) = xh(k,i)*o3den(k,i)*v24
          xohv5(k,i) = xh(k,i)*o3den(k,i)*v25(k,i,lat)
          xohv6(k,i) = xh(k,i)*o3den(k,i)*v26(k,i,lat)
          xohv7(k,i) = xh(k,i)*o3den(k,i)*v27(k,i,lat)
          xohv8(k,i) = xh(k,i)*o3den(k,i)*v28(k,i,lat)
          xohv9(k,i) = xh(k,i)*o3den(k,i)*v29(k,i,lat)
! ----------------------------------------------------------------------
!  Only one channel:   HO2 + O 
!         xohvg(k,i) = ohsrc(k,i)
!         xohv1(k,i) = xho2(k,i)*xo(k,i)*v1
!         xohv2(k,i) = xho2(k,i)*xo(k,i)*v2
!         xohv3(k,i) = xho2(k,i)*xo(k,i)*v3(k,i,lat)
!         xohv4(k,i) = xho2(k,i)*xo(k,i)*v4(k,i,lat)
!         xohv5(k,i) = xho2(k,i)*xo(k,i)*v5(k,i,lat)
!         xohv6(k,i) = xho2(k,i)*xo(k,i)*v6(k,i,lat)
!         xohv7(k,i) = xho2(k,i)*xo(k,i)*v7
!         xohv8(k,i) = xho2(k,i)*xo(k,i)*v8
!         xohv9(k,i) = xho2(k,i)*xo(k,i)*v9
!
! -----------------------------------------------------------------------
!  ** OH Vibrational states: Include both population and depopulation at each level
!                          : Cascading Production and Loss terms
!                          : Final densities at each OH-vib level
!                          : Contribution to OHIRVER from each vibrational level
! -----------------------------------------------------------------------
!         ksum = (kv98+kv97+kv96+kv95+kv94+kv93+kv92+kv91+kv9g)
!         xohv9den(k,i) = xohv9(k,i)/(xco2(k,i)*ksum)
!         xohv8(k,i) = xohv8(k,i) +  xohv9den(k,i)*xco2(k,i)*kv98
!         xohv7(k,i) = xohv7(k,i) +  xohv9den(k,i)*xco2(k,i)*kv97
!         xohv6(k,i) = xohv6(k,i) +  xohv9den(k,i)*xco2(k,i)*kv96
!         xohv5(k,i) = xohv5(k,i) +  xohv9den(k,i)*xco2(k,i)*kv95
!         xohv4(k,i) = xohv4(k,i) +  xohv9den(k,i)*xco2(k,i)*kv94
!         xohv3(k,i) = xohv3(k,i) +  xohv9den(k,i)*xco2(k,i)*kv93
!         xohv2(k,i) = xohv2(k,i) +  xohv9den(k,i)*xco2(k,i)*kv92
!         xohv1(k,i) = xohv1(k,i) +  xohv9den(k,i)*xco2(k,i)*kv91
!         xohvg(k,i) = xohvg(k,i) +  xohv9den(k,i)*xco2(k,i)*kv9g
! -----------------------------------------------------------------------
!                              EQUIVALENT
          ksum = (kv98+kv97+kv96+kv95+kv94+kv93+kv92+kv91+kv9g)
          kssum = (kv98/sv98+kv97/sv97+kv96/sv96+kv95/sv95+kv94/sv94
     |            +kv93/sv93+kv92/sv92+kv91/sv91+kv9g/sv9g)
          FOHIR = 1./(1. + kssum*xco2(k,i))
          xohv9den(k,i) = xohv9(k,i)/(xco2(k,i)*ksum)
          ohirver(k,i) =  ohirver(k,i) + xohv9(k,i)*FOHIR*BETA
          xohv8(k,i) = xohv8(k,i) +  xohv9(k,i)*kv98/ksum
          xohv7(k,i) = xohv7(k,i) +  xohv9(k,i)*kv97/ksum
          xohv6(k,i) = xohv6(k,i) +  xohv9(k,i)*kv96/ksum
          xohv5(k,i) = xohv5(k,i) +  xohv9(k,i)*kv95/ksum
          xohv4(k,i) = xohv4(k,i) +  xohv9(k,i)*kv94/ksum
          xohv3(k,i) = xohv3(k,i) +  xohv9(k,i)*kv93/ksum
          xohv2(k,i) = xohv2(k,i) +  xohv9(k,i)*kv92/ksum
          xohv1(k,i) = xohv1(k,i) +  xohv9(k,i)*kv91/ksum
          xohvg(k,i) = xohvg(k,i) +  xohv9(k,i)*kv9g/ksum
!
          ksum = (kv87+kv86+kv85+kv84+kv83+kv82+kv81+kv8g)
          kssum = (kv87/sv87+kv86/sv86+kv85/sv85+kv84/sv84+kv83/sv83
     |            +kv82/sv82+kv81/sv81+kv8g/sv8g)
          FOHIR = 1./(1. + kssum*xco2(k,i))
          xohv8den(k,i) = xohv8(k,i)/(xco2(k,i)*ksum)
          ohirver(k,i) =  ohirver(k,i) + xohv8(k,i)*FOHIR*BETA
          xohv7(k,i) = xohv7(k,i) +  xohv8den(k,i)*xco2(k,i)*kv87
          xohv6(k,i) = xohv6(k,i) +  xohv8den(k,i)*xco2(k,i)*kv86
          xohv5(k,i) = xohv5(k,i) +  xohv8den(k,i)*xco2(k,i)*kv85
          xohv4(k,i) = xohv4(k,i) +  xohv8den(k,i)*xco2(k,i)*kv84
          xohv3(k,i) = xohv3(k,i) +  xohv8den(k,i)*xco2(k,i)*kv83
          xohv2(k,i) = xohv2(k,i) +  xohv8den(k,i)*xco2(k,i)*kv82
          xohv1(k,i) = xohv1(k,i) +  xohv8den(k,i)*xco2(k,i)*kv81
          xohvg(k,i) = xohvg(k,i) +  xohv8den(k,i)*xco2(k,i)*kv8g
!
          ksum = (kv76+kv75+kv74+kv73+kv72+kv71+kv7g)
          kssum = (kv76/sv76+kv75/sv75+kv74/sv74+kv73/sv73
     |            +kv72/sv72+kv71/sv71+kv7g/sv7g)
          FOHIR = 1./(1. + kssum*xco2(k,i))
          xohv7den(k,i) = xohv7(k,i)/(xco2(k,i)*ksum)
          ohirver(k,i) =  ohirver(k,i) + xohv7(k,i)*FOHIR*BETA
          xohv6(k,i) = xohv6(k,i) +  xohv7den(k,i)*xco2(k,i)*kv76
          xohv5(k,i) = xohv5(k,i) +  xohv7den(k,i)*xco2(k,i)*kv75
          xohv4(k,i) = xohv4(k,i) +  xohv7den(k,i)*xco2(k,i)*kv74
          xohv3(k,i) = xohv3(k,i) +  xohv7den(k,i)*xco2(k,i)*kv73
          xohv2(k,i) = xohv2(k,i) +  xohv7den(k,i)*xco2(k,i)*kv72
          xohv1(k,i) = xohv1(k,i) +  xohv7den(k,i)*xco2(k,i)*kv71
          xohvg(k,i) = xohvg(k,i) +  xohv7den(k,i)*xco2(k,i)*kv7g
!
          ksum = (kv65+kv64+kv63+kv62+kv61+kv6g)
          kssum = (kv65/sv65+kv64/sv64+kv63/sv63+kv62/sv62
     |            +kv61/sv61+kv6g/sv6g)
          FOHIR = 1./(1. + kssum*xco2(k,i))
          xohv6den(k,i) = xohv6(k,i)/(xco2(k,i)*ksum)
          ohirver(k,i) =  ohirver(k,i) + xohv6(k,i)*FOHIR*BETA
          xohv5(k,i) = xohv5(k,i) + xohv6den(k,i)*xco2(k,i)*kv65
          xohv4(k,i) = xohv4(k,i) + xohv6den(k,i)*xco2(k,i)*kv64
          xohv3(k,i) = xohv3(k,i) + xohv6den(k,i)*xco2(k,i)*kv63
          xohv2(k,i) = xohv2(k,i) + xohv6den(k,i)*xco2(k,i)*kv62
          xohv1(k,i) = xohv1(k,i) + xohv6den(k,i)*xco2(k,i)*kv61
          xohvg(k,i) = xohvg(k,i) + xohv6den(k,i)*xco2(k,i)*kv6g
!
          ksum = (kv54+kv53+kv52+kv51+kv5g)
          kssum = (kv54/sv54+kv53/sv53+kv52/sv52+kv51/sv51
     |            +kv5g/sv5g)
          FOHIR = 1./(1. + kssum*xco2(k,i))
          xohv5den(k,i) = xohv5(k,i)/(xco2(k,i)*ksum)
          ohirver(k,i) =  ohirver(k,i) + xohv5(k,i)*FOHIR*BETA
          xohv4(k,i) = xohv4(k,i) + xohv5den(k,i)*xco2(k,i)*kv54
          xohv3(k,i) = xohv3(k,i) + xohv5den(k,i)*xco2(k,i)*kv53
          xohv2(k,i) = xohv2(k,i) + xohv5den(k,i)*xco2(k,i)*kv52
          xohv1(k,i) = xohv1(k,i) + xohv5den(k,i)*xco2(k,i)*kv51
          xohvg(k,i) = xohvg(k,i) + xohv5den(k,i)*xco2(k,i)*kv5g
!
          ksum = (kv43+kv42+kv41+kv4g)
          kssum = (kv43/sv43+kv42/sv42+kv41/sv41+kv4g/sv4g)
          FOHIR = 1./(1. + kssum*xco2(k,i))
          xohv4den(k,i) = xohv4(k,i)/(xco2(k,i)*ksum)
          ohirver(k,i) =  ohirver(k,i) + xohv4(k,i)*FOHIR*BETA
          xohv3(k,i) = xohv3(k,i) + xohv4den(k,i)*xco2(k,i)*kv43
          xohv2(k,i) = xohv2(k,i) + xohv4den(k,i)*xco2(k,i)*kv42
          xohv1(k,i) = xohv1(k,i) + xohv4den(k,i)*xco2(k,i)*kv41
          xohvg(k,i) = xohvg(k,i) + xohv4den(k,i)*xco2(k,i)*kv4g
!
          ksum = (kv32+kv31+kv3g)
          kssum = (kv32/sv32+kv31/sv31+kv3g/sv3g)
          FOHIR = 1./(1. + kssum*xco2(k,i))
          xohv3den(k,i) = xohv3(k,i)/(xco2(k,i)*ksum)
          ohirver(k,i) =  ohirver(k,i) + xohv3(k,i)*FOHIR*BETA
          xohv2(k,i) = xohv2(k,i) + xohv3den(k,i)*xco2(k,i)*kv32
          xohv1(k,i) = xohv1(k,i) + xohv3den(k,i)*xco2(k,i)*kv31
          xohvg(k,i) = xohvg(k,i) + xohv3den(k,i)*xco2(k,i)*kv3g
!
          ksum = (kv21+kv2g)
          kssum = (kv21/sv21+kv2g/sv2g)
          FOHIR = 1./(1. + kssum*xco2(k,i))
          xohv2den(k,i) = xohv2(k,i)/(xco2(k,i)*ksum)
          ohirver(k,i) =  ohirver(k,i) + xohv2(k,i)*FOHIR*BETA
          xohv1(k,i) = xohv1(k,i) + xohv2den(k,i)*xco2(k,i)*kv21
          xohvg(k,i) = xohvg(k,i) + xohv2den(k,i)*xco2(k,i)*kv2g
!
          kssum = (kv1g/sv1g)
          FOHIR = 1./(1. + kssum*xco2(k,i))
          xohv1den(k,i) = xohv1(k,i)/(xco2(k,i)*kv1g)
          ohirver(k,i) =  ohirver(k,i) + xohv1(k,i)*FOHIR*BETA
          xohvg(k,i) = xohvg(k,i) + xohv1(k,i)
 
          xoh2(k,i) = xohvg(k,i)/ohloss(k,i)
          endif
! --------------------------------------------------------
          if (ohsrc(k,i) .GT. 1.0) THEN
            ohsrc(k,i)= log10(ohsrc(k,i))
          else
            ohsrc(k,i) = 0.
          endif
          if (ohlossr(k,i) .GT. 1.0) THEN
            ohlossr(k,i)= log10(ohlossr(k,i))
          else
            ohlossr(k,i) = 0.
          endif
          if (ohden(k,i) .GT. 1.0) THEN
            ohden(k,i) = log10(ohden(k,i))
          else
            ohden(k,i) = 0.
          endif
!         if (o3src(k,i) .GT. 1.0) THEN
!           o3src(k,i)= log10(o3src(k,i))
!         else
!           o3src(k,i) = 0.
!         endif
!         if (o3lossr(k,i) .GT. 1.0) THEN
!           o3lossr(k,i)= log10(o3lossr(k,i))
!         else
!           o3lossr(k,i) = 0.
!         endif
          if (o3den(k,i) .GT. 1.0) THEN
            o3den(k,i) = log10(o3den(k,i))
          else
            o3den(k,i) = 0.
          endif
          if (xoh2(k,i) .GT. 1.0) THEN
            xoh2(k,i) = log10(xoh2(k,i))
          else
            xoh2(k,i) = 0.
          endif
          if (xohv1den(k,i) .GT. 1.0) THEN
            xohv1den(k,i) = log10(xohv1den(k,i))
          else
            xohv1den(k,i) = 0.
          endif
          if (xohv2den(k,i) .GT. 1.0) THEN
            xohv2den(k,i) = log10(xohv2den(k,i))
          else
            xohv2den(k,i) = 0.
          endif
          if (xohv3den(k,i) .GT. 1.0) THEN
            xohv3den(k,i) = log10(xohv3den(k,i))
          else
            xohv3den(k,i) = 0.
          endif
          if (xohv4den(k,i) .GT. 1.0) THEN
            xohv4den(k,i) = log10(xohv4den(k,i))
          else
            xohv4den(k,i) = 0.
          endif
          if (xohv5den(k,i) .GT. 1.0) THEN
            xohv5den(k,i) = log10(xohv5den(k,i))
          else
            xohv5den(k,i) = 0.
          endif
          if (xohv6den(k,i) .GT. 1.0) THEN
            xohv6den(k,i) = log10(xohv6den(k,i))
          else
            xohv6den(k,i) = 0.
          endif
          if (xohv7den(k,i) .GT. 1.0) THEN
            xohv7den(k,i) = log10(xohv7den(k,i))
          else
            xohv7den(k,i) = 0.
          endif
          if (xohv8den(k,i) .GT. 1.0) THEN
            xohv8den(k,i) = log10(xohv8den(k,i))
          else
            xohv8den(k,i) = 0.
          endif
          if (xohv9den(k,i) .GT. 1.0) THEN
            xohv9den(k,i) = log10(xohv9den(k,i))
          else
            xohv9den(k,i) = 0.
          endif
!  Background densities
          if (xho2(k,i) .GT. 1.0) THEN
            xho2(k,i) = log10(xho2(k,i))
          else
            xho2(k,i) = 0.
          endif
          if (xh(k,i) .GT. 1.0) THEN
            xh(k,i) = log10(xh(k,i))
          else
            xh(k,i) = 0.
          endif
          if (xo(k,i) .GT. 1.0) THEN
            xo(k,i) = log10(xo(k,i))
          else
            xo(k,i) = 0.
          endif
!  Calculated VMR in log10 units
          if (xovmr(k,i) .GT. 1.0E-12) THEN
            xovmr(k,i) = log10(xovmr(k,i))
          else
            xovmr(k,i) = -12.
          endif
          if (xhvmr(k,i) .GT. 1.0E-12) THEN
            xhvmr(k,i) = log10(xhvmr(k,i))
          else
            xhvmr(k,i) = -12.
          endif
          if (xho2vmr(k,i) .GT. 1.0E-12) THEN
            xho2vmr(k,i) = log10(xho2vmr(k,i))
          else
            xho2vmr(k,i) = -12.
          endif
          if (o3denvmr(k,i) .GT. 1.0E-12) THEN
            o3denvmr(k,i) = log10(o3denvmr(k,i))
          else
            o3denvmr(k,i) = -12.
          endif
          if (ohdenvmr(k,i) .GT. 1.0E-12) THEN
            ohdenvmr(k,i) = log10(ohdenvmr(k,i))
          else
            ohdenvmr(k,i) = -12.
          endif
! -----------------------------------------------------------------------
!  OH and O3 Chemical Sources and Losses
! -----------------------------------------------------------------------
!     call addfld('XOH',' ',' ',ohden,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('OHVMR',' ',' ',ohdenvmr,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('OHSRC',' ',' ',ohsrc,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('OHLOSSR',' ',' ',ohlossr,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('XO3',' ',' ',o3den,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('O3VMR',' ',' ',o3denvmr,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('O3SRC',' ',' ',o3src,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('O3LOSSR',' ',' ',o3lossr,'lev',k0,k1,'lon',i0,i1,lat)
! -----------------------------------------------------------------------
!  OH Vibrational Levels and Contributors: Cummulative end result.
! -----------------------------------------------------------------------
!     Number Densities
!     call addfld('XO3',' ',' ',o3den,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('XOH',' ',' ',ohden,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('XHO2',' ',' ',xho2,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('XH',' ',' ',xh,'lev',k0,k1,'lon',i0,i1,lat)
!     call addfld('XO',' ',' ',xo,'lev',k0,k1,'lon',i0,i1,lat)
!     Volume Mixing Ratios 
      call addfld('XO3',' ',' ',o3denvmr,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('XOH',' ',' ',ohdenvmr,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('XHO2',' ',' ',xho2vmr,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('XH',' ',' ',xhvmr,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('XO',' ',' ',xovmr,'lev',k0,k1,'lon',i0,i1,lat)
! ---------------------------------------------------------------------
      call addfld('XOHV1',' ',' ',xohv1den,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('XOHV2',' ',' ',xohv2den,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('XOHV3',' ',' ',xohv3den,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('XOHV4',' ',' ',xohv4den,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('XOHV5',' ',' ',xohv5den,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('XOHV6',' ',' ',xohv6den,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('XOHV7',' ',' ',xohv7den,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('XOHV8',' ',' ',xohv8den,'lev',k0,k1,'lon',i0,i1,lat)
      call addfld('XOHV9',' ',' ',xohv9den,'lev',k0,k1,'lon',i0,i1,lat)
! ---------------------------------------------------------------------
!  *  OH IR NIGHTGLOW (VER Processing) 
! ---------------------------------------------------------------------
!     Photons/cm3.sec
      slab(k,i)=0.5*gask*tni(k,i)/(barm(k,i)*grav)
!  *  -----------------------------------------------------------
!     Log10 (Photons/cm3.sec)
        if (ohirver(k,i) .GT. 1.0) THEN
           ohirverl(k,i)= alog10(ohirver(k,i))
        else
           ohirverl(k,i) = 0.
        endif
!  *  -----------------------------------------------------------
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
! ----------------------------------------------------------------------
!   * DISPOSE OH INTEGRATED EMISSION RATE (PH/cm2/sec)
      do i=lon0,lon1
!  Below 90 km is zeroed out
        do k=1,11
             fldht(k) = 0.
        enddo
!  Above 90 km is retained
        do k=11,lev1-1
             fldht(k) = slab(k,i)*ohirver(k,i)
        enddo
        call integral(colout,fldht,1,dz,nk)
        do k = lev0,lev1-1
!    * MEGA-RAYLEIGH UNITS FOR PLOTS (Oh IR INTEGRATED INTENSITY)
!            ohirint(k,i) = colout(k)/1.0E+12
!    * Kilo-RAYLEIGH UNITS FOR PLOTS (Oh IR INTEGRATED INTENSITY)
             ohirint(k,i) = colout(k)/1.0E+09
        enddo
      enddo
! ----------------------------------------------------------------------
!   * OH VOLUME EMISSION RATE (log10 PH/cm3/sec)
      call addfld('OHIRVER',' ',' ',ohirverl,'lev',k0,k1,'lon',i0,i1,
     | lat)
!   * OH INTEGRATED VERTICAL INTENSITY (MR = 1.0E+12*PH/cm2/sec)
      call addfld('OHIRINT',' ',' ',ohirint,'lev',k0,k1,'lon',i0,i1,
     | lat)
! ----------------------------------------------------------------------
      end subroutine comp_oh
!
