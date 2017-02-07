!
      module chemrates_module
!
! 11/2/06 btf: chemrates for venus, as per Bougher vtgcm2 notes.
! 11/15/06 swb: double-checked chemrates for venus (OK). All now
!	        compatible with Fox and Sung (2001)
! 11/20/06 btf: allocate fs array at end of tdep subroutine
! 11/30/06 swb: pass te and ti from settei for use in chemrates
! 10/09/07 swb: NOX rates updated/added for N4S, N2D, NO, NO* calculation
! 10/30/07 swb: r2(k,i,lat) updated (Slanger et al., 2006; Pejakovic et al., 2007)
! 03/17/08 swb: NOX rates updated from Gerard et al., (2008)
! 04/01/08 asb: added low altitude chem rates (Yung & DeMore, 1982)
! 04/04/08 swb: checks of low altitude chem rates (Yung & DeMore, 1982)
! 04/02/09 asb: included OH chemical rates for OH airglow
!
! btf 6/13/13: Add SOx chemical rates for SO2 and SO reactions
!              from Zhang et al (2012). Icarus, 217. Table 2.
! 01/27/14 swb: included CLx reactions for Zhang et al (2012).
! 05/08/14 swb: r2(k,i,lat) and r5(k,i,lat) updated (Zhang et al., 2012)
! 05/12/14 swb: r44, r108 updated (Zhang et al., 2012)
! 5/4/16   swb: Update r99 and  r160  reactions (T>210 K only)
!               Zhang et al (2012). Icarus, 217. Table 2.
! 05/27/14 swb: add r167 from Zhang et al., (2012)
! 06/13/14 swb: add r346 from Zhang et al., (2012)
! 10/08/14 swb: add s248 and s295 from Zhang et al (2012)
! 07/19/16 swb: checking all hox and o3 rates are in place and correct
! 07/20/16 swb: turn on OH vibrational level setup
! 07/22/16 swb: r41b updated rate
! 09/01-06/16 swb: Revised OH vibrational level setup (2-production schemess)
! 10/11/16 swb: updated OH vibrational level setup (2-production schemess)
! 11/19/16 swb: Add O2 source reaction r63 (Zhang et al 2012)
! 12/05/16 swb: adding r123 and r130 for O2 srcs
! 12/06/16 swb: adding all co_updates from 5/4/2016 baseline
! 12/06/16 swb: clean up ordering of reaction rates
!
      use params_module,only: nlevp1,nlonp4
      implicit none
      real,parameter :: preset=0.
!
! Matrices for O,CO (set by comp_oco.F, used by comp.F)
      real,allocatable :: fs(:,:,:,:,:) ! (i,k,2,0:2,j)
!
! Temperature independent reaction rate constants:
      real,parameter ::
!
! Ion-neutral chemistry (constants):
     | r18 = 1.1e-9,   ! CO2 + O+ [cm3/sec]
     | r19 = 9.6e-11,  ! O + CO2+  to O+
     | r20 = 1.64e-10, ! O + CO2+  to O2+ 
     | r23 = 1.1e-9,   ! 
     | r26 = 4.5e-10,  ! 
     | r28 = 1.0e-10,  ! N4S + O2+
!    | r28 = 0.,  ! N4S + O2+
!    | r29 = 5.0e13,   ! N2 + O+
!
! N2D Quenching reactions (s.a., rk1-rk5, minus rk2, T-dependent):
!    | rk1 = 3.6e-13,  ! N2D+CO2 [cm3/sec]
     | rk1 = 2.8e-13,  ! N2D+CO2 [cm3/sec]
     | rk3 = 1.9e-12,  ! N2D+CO
!    | rk3 = 4.0e-12,  ! N2D+CO
!    | rk4 = 6.9e-13,  ! N2D+O
!    | rk4 = 2.0e-12,  ! N2D+O
     | rk4 = 2.0e-11,  ! N2D+O
     | rk5 = 1.7e-14,  ! N2D+N2
!
! N4S losses:
!    | rg2 = 2.26e-32,   ! O + N4S + CO2 -> NO + CO2 [cm6/sec]
!    | rg2 = 3.46e-32,   ! O + N4S + CO2 -> NO + CO2 [cm6/sec]
!    | rg3 = 3.4e-11,  ! N4S + NO -> N2 + O
     | rk6 = 1.7e-16,  ! N4S+CO2
!
! Msc:
     | rk10 = 4.5e-10, ! NO+ O2+
     | rk11 = 6.7e-11, ! NO+ N2D
!
! Oxygen Loss
!    | r2 = 2.75e-32,   ! O+O+CO2 (Gerard et al., 2008 / Campbell and Grey 1973)
!    | r2 = 2.5e-32,    ! O+O+CO2 (Slanger et al., 2006)
!    | r2 = 1.8e-32,    ! O+O+CO2 (for ~200K Recommended by Huestis from lab measurments by Corey Jamieson)
!    | r2 = 1.2e-32,    ! O+O+CO2 (Recommended by Huestis from lab measurments by Corey Jamieson)
!    | r2(k,i,lat),     ! O+O+CO2 (Recommended by Zhang et al. 2012, R59) ********** Use T-dependent
!    | r5 = 1.35e-33,   ! O+O2+CO2 (Yung and DeMore 1982)
!    | r5 = 4e-33,      ! O+O2+CO2 (Recommended by Huestis JPL current evaluation for N2 x 2.5)
!    | r5(k,i,lat),     ! O+O2+CO2 (Recommended by Zhang et al. 2012, R65) ********* Use T-dependent
!
! Low altitude tracers and  OH Losses:
!    | rk28 = 3.8e-11,    ! OH+O = H+O2    (Yung and DeMore,1982)??
     | r25 = 1.5e-13,     ! OH+CO = H+CO2 [Sander et al 2006]
!    | r31 = 3.1e-11,     ! HO2+O = OH+O2  *** SEE R49 ***
     | r51 = 1.5e-13,     ! OH+CO -> H+CO2 (Liang and Yung,2009/ Saunder et al. 2006)
     | r102a = 3.0e-11,   ! CLCO +O = CO2+CL (Yung and DeMore 1982)
     | r102b = 3.0e-12,   ! CLCO+O = CO+CLO (Yung and DeMore 1982)
     | r144 = 3.0e-12,    ! O+CLCO = CLO+CO [Zhang et al 2012]
     | r146 = 1.0e-11,    ! H+CLCO = HCL+CO [Zhang et al 2012]
     | r147 = 1.5e-10,    ! OH+CLCO = CLOH+CO [Zhang et al 2012]
!
! SOx reactions:  Zhang et al. (2012). Icarus, 217. Table 2.
     | s133 = 2.8E-11,   ! ClO+SO --> SO2 + Cl  (Sander et al. 2002)
     | s166 = 2.3E-12,   ! S + O2 --> SO  + O   (Sander et al. 2006)
     | s240 = 5.1E-31,   ! O + SO + M --> SO2 + M  (Atkinson et al. 2004)
!    | s240 = 5.1E-31,   ! O + SO + M --> SO2 + M  (Singleton and Cvetanovic, 1988)
     | s244 = 1.0E-11,   ! ClCO3 + SO --> SO2 + X  (Mills, 1998)
     | s249 = 4.4E-31,   ! 2SO + M -->  X + M  (Out for curretn SOX simulations)
     | s261 = 1.0E-15,   ! ClCO3 + SO2 --> SO3 + X  (Mills, 1998)
     | s288 = 2.3E-11,   ! Cl + OSCl  --> Cl2 + SO  (Miles, 1998; Mills et al. 2007)
     | s293 = 1.0E-11,   ! O + ClSO2  --> SO2 + ClO  (Mills, 1998)
     | s294 = 1.0E-11,   ! H + ClSO2  --> SO2 + HCl  (Mills, 1998)
     | s295 = 1.0E-20,   ! Cl + ClSO2  --> SO2 + Cl2 (Moses, 2002)
     | s300 = 5.0E-13,   ! 2ClSO2 --> 2SO2 + Cl2   (Moses, 2002)
     | s346 = 1.4E-11,   ! SO + NO2 --> SO2 + NO   (Sander et al., 2002)
!
!
! O3 chem
     | k55 = 1.4e-17,  ! clo+o3->cloo+o2
!
!  OH Vibrational 
     |  v21  = 0.,
     |  v22  = 0.,
     |  v23  = 0.,
     |  v24  = 0.,
     |  v1  = 0.,
     |  v2  = 0.,
     |  v7  = 0.,
     |  v8  = 0.,
     |  v9  = 0.,
     | kv98 = 4.20e-12,
     | kv97 = 4.00e-12,
     | kv96 = 3.80e-12,
     | kv95 = 3.60e-12,
     | kv94 = 3.40e-12,
     | kv93 = 3.20e-12,
     | kv92 = 3.10e-12,
     | kv91 = 2.90e-12,
     | kv9g = 2.80e-12,
     | kv87 = 3.30e-12,
     | kv86 = 2.50e-12,
     | kv85 = 1.90e-12,
     | kv84 = 1.40e-12,
     | kv83 = 1.00e-12,
     | kv82 = 8.00e-13,
     | kv81 = 6.00e-13,
     | kv8g = 4.00e-13,
     | kv76 = 3.20e-12,
     | kv75 = 2.30e-12,
     | kv74 = 1.60e-12,
     | kv73 = 1.20e-12,
     | kv72 = 9.00e-13,
     | kv71 = 6.00e-13,
     | kv7g = 4.00e-13,
     | kv65 = 1.10e-12,
     | kv64 = 6.00e-13,
     | kv63 = 3.00e-13,
     | kv62 = 1.00e-13,
     | kv61 = 1.00e-13,
     | kv6g = 1.00e-13,
     | kv54 = 1.60e-12,
     | kv53 = 6.00e-13,
     | kv52 = 2.00e-13,
     | kv51 = 1.00e-13,
     | kv5g = 0.00,
     | kv43 = 1.00e-12,
     | kv42 = 2.00e-13,
     | kv41 = 1.00e-13,
     | kv4g = 0.00,
     | kv32 = 7.00e-13,
     | kv31 = 1.00e-13,
     | kv3g = 0.00,
     | kv21 = 4.00e-13,
     | kv2g = 0.00,
     | kv1g = 2.00e-13,
     | sv98 = 1.89E+01,
     | sv97 = 1.12E+02,
     | sv96 = 5.28E+01,
     | sv95 = 1.11E+01,
     | sv94 = 1.70E+00,
     | sv93 = 1.78E-01,
     | sv92 = 0.00E+00,
     | sv91 = 0.00E+00,
     | sv9g = 0.00E+00,
     | sv87 = 6.83E+00,
     | sv86 = 1.15E+02,
     | sv85 = 3.91E+01,
     | sv84 = 6.64E+00,
     | sv83 = 7.67E-01,
     | sv82 = 5.74E-02,
     | sv81 = 0.00E+00,
     | sv8g = 0.00E+00,
     | sv76 = 1.85E+00,
     | sv75 = 1.07E+02,
     | sv74 = 2.71E+01,
     | sv73 = 3.53E+00,
     | sv72 = 2.93E-01,
     | sv71 = 1.35E-02,
     | sv7g = 0.00E+00,
     | sv65 = 3.16E+00,
     | sv64 = 9.12E+01,
     | sv63 = 1.70E+01,
     | sv62 = 1.61E+00,
     | sv61 = 8.53E-02,
     | sv6g = 1.80E-03,
     | sv54 = 8.71E+00,
     | sv53 = 7.12E+01,
     | sv52 = 9.35E+00,
     | sv51 = 5.68E-01,
     | sv5g = 1.42E-02,
     | sv43 = 1.59E+01,
     | sv42 = 4.92E+01,
     | sv41 = 4.08E+00,
     | sv4g = 1.20E-01,
     | sv32 = 2.20E+01,
     | sv31 = 2.80E+01,
     | sv3g = 1.11E+00,
     | sv21 = 2.37E+01,
     | sv2g = 1.05E+01,
     | sv1g = 1.77E+01
!
! Temperature dependent chemical reaction rate coefficients.
! These are set by rates_tdep for every latitude at every timestep.
!
      real,dimension(:,:,:),allocatable :: ! (nlevp1,lon0:lon1,lat0:lat1)
!
!OH Chem
     |  v3,
     |  v4,
     |  v5,
     |  v6,
     |  v25,
     |  v26,
     |  v27,
     |  v28,
     |  v29,
!
! N2 quenching:
     |  rk2,   ! N2D+O2
!
! Dissociative recombination:
     |  a13,   ! O2+ +e
     |  a14,   ! CO2+ +e
     |  a15,   ! N2+ + e [cm3/sec]
     |  a17,   ! NO+ + e  (0.75 N2D, 0.25 N4S)
!
! Ti-dependent ion-netural reactions:
     |  r21,   ! CO2+ N2+ [cm3/sec]
     |  r22,   ! 
     |  r27,   ! O + N2+
     |  r29,   ! N2 + O+
!
! NOx rates from Gerard et al. (2008); Table 1
     |  rg1,   ! N4S + O -> NO* + hv
     |  rg2,   ! O + N4S + CO2 -> NO + CO2 [cm6/sec]
     |  rg3,   ! N4S + NO -> N2 + O
!
! 3-body recombination rates:
!    |  r2,    ! O+O+CO2  (Slanger et al., 2006) [cm6/sec]
!    |  r2,    ! O+O+CO2  (case 1 and 3 for Marie-Eve) [cm6/sec]
     |  r2,    ! O+O+CO2  (Recommended by Zhang et al. 2012, R59) [cm6/sec]
     |  r3,    ! O+CO+CO2 [cm6/sec]
     |  r5,    ! O+O2+CO2 (Recommended by Zhang et al. 2012, R65) [cm6/sec]
!
! OH Sources:
     |  r41b,  ! HO2+Cl->OH+ClO  (Saunder et al. 2006)
     |  r49,   ! HO2+O->OH+O2    (Liang and Yung,2009/Saunder et al. 2006)
     |  r54,   ! HO2+O3->OH+2O2  (Liang and Yung,2009/Saunder et al. 2006)
     |  r66,   ! HCl+O->OH+Cl    (Liang and Yung,2009/Saunder et al. 2006)
     |  r104,  ! H2+O->OH+H      (Liang and Yung,2009)
     |  rk23,  ! h+O3->OH+O2  (Saunder et al. 2006)
!
! OH Sinks:
     |  r53,   !OH+O3->HO2+O2    (Liang and Yung,2009/Saunder et al. 2006)
     |  r105,  !2OH->H2O+O       (Liang and Yung,2009)
     |  rk28,  !OH+O->H+O2       (Saunder et al. 2006) 
!
! O3 chem:
     |  k50,   !2o+o2->o3+o
     |  k51,   !o+2o2->o3+o2
     |  k52,   !o+o2+co->o3+co
     |  k53,   !no+o3->no2+o2  (Saunder et al. 2006)
     |  k54,   !cl+o3->clo+o2  (Atkinson et al. 2007)
     |  k56,   !clo+o3->oclo+o2
     |  r63,   !o+o2+n2->o3+n2
     |  r123,  !clo+o->cl+o2   (Zhang et al 2012)
     |  r130,  !2clo->cl2+o2   (Zhang et al 2012)
!
! NO Tracer chem:
     |  k45,   ! no+o+M (Sander et al. 2006/Krasnopolsky 2010)
     |  k46,   ! no+o3 (Sander et al. 2006/Krasnopolsky 2010)
     |  k47,   ! no+ho2  (Sander et al. 2006/Krasnopolsky 2010)
!
! Low Altitude Tracers (Yung and Demore, 1982)
     |  r15,   ! H2+O = OH+H
!    |  r16,   ! HCL+O = OH+CL
     |  r43,   ! CLO+CO = CL+CO2
     |  r44,   ! CLO+O = CL+O2 (Atkinson et al. 2007)
     |  r99,   ! CL+CO+M = CLCO+M [cm6/sec] (Niovich et al. 1990 For N2 and Corrected for w/factor)
     |  r108,  ! CL2+O = CLO+CL
     |  r111,  ! CL+N2+CO = CLCO+N2 [cm6/sec] (Zhnag et al 2012)
     |  r160,  ! CLCO+CO2 = CL+CO2+CO [cm6/sec] (Zhnag et al 2012)
     |  r162,  ! O+2CO = CO2+CO [cm6/sec] (Zhnag et al 2012)
     |  r163,  ! 2O+CO = CO2+O [cm6/sec] (Zhnag et al 2012)
     |  rk20,  ! O3+O = 2O2  (Saunder et al. 2006)
!
! SOx reactions:  Zhang et al. (2012). Icarus, 217. Table 2.
     | s114,   ! Cl + SO2 + M --> X
     | s167,   ! SO + O2  --> SO2 + O
     | s243,   ! Cl + SO  + M --> OSCl + M
     | s248,   ! 2SO  --> SO2 + S
     | s258,   ! O + SO2 + M --> SO3 + M
     | s286,   ! O + OSCl --> SO2 + Cl
     | s287,   ! O + OSCl --> SO + ClO
     | s299    ! SO + ClSO2 --> SO2 + OSCl
!
      contains
!-----------------------------------------------------------------------
      subroutine alloc_tdep
!
! Allocate temperature-dependent reaction rates for task subdomain:
! Called once per run from init_fields.
!
      use mpi_module,only: lon0,lon1,lat0,lat1
      integer :: istat
!  OH Vibrational Levels for (1-0) transition for OH NTGL
!
      allocate(v3(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' v3: stat=',i3)") istat
!
      allocate(v4(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' v4: stat=',i3)") istat
!
      allocate(v5(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' v5: stat=',i3)") istat
!
      allocate(v6(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' v6: stat=',i3)") istat
!
      allocate(v25(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' v25: stat=',i3)") istat
!
      allocate(v26(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' v26: stat=',i3)") istat
!
      allocate(v27(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' v27: stat=',i3)") istat
!
      allocate(v28(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' v28: stat=',i3)") istat
!
      allocate(v29(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' v29: stat=',i3)") istat
!
      allocate(rk2(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' rk2: stat=',i3)") istat
!
      allocate(a13(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' a13: stat=',i3)") istat
!
      allocate(a14(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' a14: stat=',i3)") istat
!
      allocate(a15(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' a15: stat=',i3)") istat
!
      allocate(a17(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' a17: stat=',i3)") istat
!
      allocate(r21(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' r21: stat=',i3)") istat
!
      allocate(r22(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' r22: stat=',i3)") istat
!
      allocate(r27(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' r27: stat=',i3)") istat
!
      allocate(r29(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' r29: stat=',i3)") istat
!
      allocate(rg1(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' rg1: stat=',i3)") istat
!
      allocate(rg2(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' rg2: stat=',i3)") istat
!
      allocate(rg3(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' rg3: stat=',i3)") istat
!
      allocate(r2(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' r2: stat=',i3)") istat
!
      allocate(r3(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' r3: stat=',i3)") istat
!
      allocate(r5(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' r5: stat=',i3)") istat
!
! LOW ALTITUDE TRACERS
      allocate(rk20(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' rk20: stat=',i3)") istat
!
      allocate(r15(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' r15: stat=',i3)") istat
!
!     allocate(r16(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
!     if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
!    |  ' r16: stat=',i3)") istat
!
      allocate(r43(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' r43: stat=',i3)") istat
!
      allocate(r44(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' r44: stat=',i3)") istat
!
      allocate(r99(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' r99: stat=',i3)") istat
!
      allocate(r108(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' r108: stat=',i3)") istat
!
      allocate(r111(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' r111: stat=',i3)") istat
!
      allocate(r160(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' r160: stat=',i3)") istat
!
      allocate(r162(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' r162: stat=',i3)") istat
!
      allocate(r163(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' r163: stat=',i3)") istat

!
! O3 CHEM
!
      allocate(k45(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' k45: stat=',i3)") istat
!
      allocate(k46(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' k46: stat=',i3)") istat
!
      allocate(k47(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' k47: stat=',i3)") istat
!
      allocate(k50(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' k50: stat=',i3)") istat
!
      allocate(k51(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' k51: stat=',i3)") istat
!
      allocate(k52(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' k52: stat=',i3)") istat
!
      allocate(k53(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' k53: stat=',i3)") istat
!
      allocate(k54(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' k54: stat=',i3)") istat
!
      allocate(k56(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' k56: stat=',i3)") istat
!
      allocate(r63(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' r63: stat=',i3)") istat
!
      allocate(r123(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' r123: stat=',i3)") istat
!
      allocate(r130(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' r130: stat=',i3)") istat
!
! OH CHEM
      allocate(r41b(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' r41b: stat=',i3)") istat
!
      allocate(r49(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' r49: stat=',i3)") istat
!
      allocate(r53(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' r53: stat=',i3)") istat
!
      allocate(r54(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' r54: stat=',i3)") istat
!
      allocate(r66(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' r66: stat=',i3)") istat
!
      allocate(r104(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' r104: stat=',i3)") istat
!
      allocate(r105(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' r105: stat=',i3)") istat
!
      allocate(rk23(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' rk23: stat=',i3)") istat
!
      allocate(rk28(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' rk28: stat=',i3)") istat
!
      allocate(fs(lon0:lon1,nlevp1,2,0:2,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' fs: stat=',i3)") istat
!
! SOX CHEM
      allocate(s114(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' s114: stat=',i3)") istat
!
      allocate(s167(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' s167: stat=',i3)") istat
!
      allocate(s243(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      allocate(s243(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' s243: stat=',i3)") istat
!
      allocate(s248(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' s248: stat=',i3)") istat
!
      allocate(s258(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' s258: stat=',i3)") istat
!
      allocate(s286(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' s286: stat=',i3)") istat
!
      allocate(s287(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' s287: stat=',i3)") istat
!
      allocate(s299(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' s299: stat=',i3)") istat
!

      end subroutine alloc_tdep
!-----------------------------------------------------------------------
! 11/2/06 btf: For now, assume tn==te==ti, so pass only tn.
!              When emperical te,ti are available, add them to
!              arg list as above.
! 11/30/06 swb: Pass te and ti from settei
!      
      subroutine chemrates_tdep(tn,te,ti,lev0,lev1,lon0,lon1,
     |  lat)
!
! Calculate temperature-dependent reaction rates (called at each latitude)
!
      use input_module,only: f107  ! 10.7 cm flux (from input and/or gpi)
      use init_module,only: sfeps  ! flux variation from orbital excentricity
      use cons_module,only: check_exp
      use addfld_module,only: addfld
!
! Args:
      integer,intent(in) ::
     |  lev0,lev1,           ! first and last level indices, this task
     |  lon0,lon1,           ! first and last longitude indices, this task
     |  lat                  ! latitude index
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) ::
     |  tn,     ! neutral temperature (deg K)
     |  te,     ! electron temperature (deg K)
     |  ti      ! ion temperature (deg K)
!
! Local:
      integer :: k,i
      real :: tncap
!     real,dimension(lev0:lev1,lon0-2:lon1+2) :: te,ti ! temporary locals
!
! expo() (util.F) is used only if check_exp is true. This will avoid
! NaNS fpe, but will degrade performance. Check_exp is in cons.F.
      real,external :: expo
!
! Set te=ti=tn until te and ti can be passed in through arg list:
!     te = tn  ! whole-array
!     ti = tn  ! whole-array
! Get empirical Te and Ti from settei.F
!
      do i=lon0,lon1
        do k=lev0,lev1-1
!         write(6,"('chemrates_tdep: k=',i2,' i=',i2,' lat=',i2)")
!    |      k,i,lat
!
!***************************************************************
!     call addfld('TECHEM',' ',' ',te(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('TICHEM',' ',' ',ti(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!***************************************************************
!    Limiter for r160 and r99 reactions (low Ts)
          tncap = tn(k,i)
!         if (tncap < 220.0) tncap = 220.0
!         if (tncap < 200.0) tncap = 200.0
          if (tncap < 210.0) tncap = 210.0
!***************************************************************
!     H + O3 reaction rates for OH(vib): v = 5-9  (revised 161011)
!     H + O3 reaction rates for OH(vib): v = 6-9  (revised 161219)
!    -----------------------------------------------------------
!         v25(k,i,lat) = 1.4e-12*exp(-470./tn(k,i))
          v25(k,i,lat) = 0.0
          v26(k,i,lat) = 4.2e-12*exp(-470./tn(k,i))
          v27(k,i,lat) = 2.1e-11*exp(-470./tn(k,i))
          v28(k,i,lat) = 4.76e-11*exp(-470./tn(k,i))
          v29(k,i,lat) = 6.58e-11*exp(-470./tn(k,i))
!
!     HO2 + O reaction rates for OH(vib): v = 3-6
!    -----------------------------------------------------------
          v3(k,i,lat) = 3.0e-13*exp(200./tn(k,i))
          v4(k,i,lat) = 4.5e-12*exp(200./tn(k,i))
          v5(k,i,lat) = 1.0e-11*exp(200./tn(k,i))
          v6(k,i,lat) = 1.4e-11*exp(200./tn(k,i))
!***************************************************************
!      DR reactions:
          a13(k,i,lat) = 1.95e-7 *(300./te(k,i))**0.7   ! O2P+e (FoxandSung,2001; MehrBiondi69
          a14(k,i,lat) = 3.5e-7  *(300./te(k,i))**0.5   ! CO2P+e
          a15(k,i,lat) = 2.196e-7*(300./te(k,i))**0.39  ! N2P+e
          a17(k,i,lat) = 4.0e-7  *(300./te(k,i))**0.5   ! NOP+e
!         a17a(k,i,lat) = 4.0e-7  *(300./te(k,i))**0.5   ! NOP+e->N(4S)
!         a17b(k,i,lat) = 4.0e-7  *(300./te(k,i))**0.5   ! NOP+e->N(2D)
!
          r21(k,i,lat) = 9.00e-10*(300./tn(k,i))**0.23
          r22(k,i,lat) = 1.50e-10*(300./tn(k,i))**1.1
          r27(k,i,lat) = 1.33e-10*(300./tn(k,i))**0.44
          r29(k,i,lat) = 1.20e-12*(300./tn(k,i))**0.45
!
!***************************************************************
!      NOX reactions:
          rg1(k,i,lat) = 1.9e-17*(300./tn(k,i))**0.5*
     |                  (1.-0.57/tn(k,i)**0.5)          !N4S+O->NO*
!         rg1(k,i,lat) = (2.8*(300./tn(k,i))**0.5)*
!    |                  (1.0*exp(-17.))                 !N4S+O->NO* (Dalgarno, 1992)
!         rg2(k,i,lat) = 2.8e-32*(300./tn(k,i))**0.5    !O+N4S+CO2 (Gerard/Campbell+Thrush)
!         rg2(k,i,lat) = 1.9e-31*(tn(k,i))**0.5         !O+N4S+M   (Yung+DeMore/Stewart)
!         rg2(k,i,lat) = 1.1e-32*(300./tn(k,i))**0.5    !O+N4S+M   (Bougher(90)/Baulch(73))
          rg2(k,i,lat) = 1.83e-32*(298./tn(k,i))**0.5   !O+N4S+CO2 (Campbell and Thrush, 66)
          rg3(k,i,lat) = 2.5e-10*(tn(k,i)/300.)**0.5*
     |                  exp(-600./tn(k,i))              ! NO+N4S->N2+O (Gerard/Fox 1994)
          rk2(k,i,lat) = 9.7e-12*exp(-185./tn(k,i))     ! N2D+O2
!
!***************************************************************
!      3-Body reactions:
!         r2(k,i,lat) = 8.6e-28/tn(k,i)**2.             ! O+O+CO2  (Yung & DeMore, 1999)
!         r2(k,i,lat) = 2.75e-32*(200./tn(k,i))**3.3    ! O+O+CO2  (Slanger,06; Pejakovic.07)
!         r2(k,i,lat) = 1.5e-27*(1./tn(k,i))**2.0       ! O+O+CO2  (Krasnopolsky,10; CampbellandGray, 1973)
!         r2(k,i,lat) = 1.2e-32*(300./tn(k,i))**2.0     ! O+O+CO2  (Krasnopolsky,95)
!         r2(k,i,lat) = 2.5*3.0e-33*(300./tn(k,i))**3.25 ! O+O+CO2 (Smith and Robertson)
          r2(k,i,lat) = 3.22e-28*(1./tn(k,i))**2.0      ! O+O+CO2  (Zhang et al. 2012, R59)
!         r3(k,i,lat) = 6.5e-33*exp(-2180./tn(k,i))     ! O+CO+CO2
          r3(k,i,lat) = 1.6e-32*exp(-2184./tn(k,i))     ! O+CO+CO2  (Yung & DeMore, 1999) CO_update
!         r5(k,i,lat) = 5.0e-28/tn(k,i)**2.3            ! O+O2+CO2
          r5(k,i,lat) = 1.4e-33*(300./tn(k,i))**2.5     ! O+O2+CO2  (Zhang et al. 2012, R65)
!
!***************************************************************
! LOW ALTITUDE TRACERS
!         rk20(k,i,lat) = 1.5e-11*exp(-2218./tn(k,i))   ! O3+O = 2O2    (Yung & DeMore, 1982)
          rk20(k,i,lat) = 8.0e-12*exp(-2060./tn(k,i))   ! O3+O = 2O2    (Saunder et al. 2006)
!         r15(k,i,lat) = 1.6e-11*exp(-4570./tn(k,i))    ! *** SEE R104 *** O+H2 = OH+H   (Yung & DeMore, 1982)
!         r16(k,i,lat) = 1.1e-11*exp(-3370./tn(k,i))    ! *** SEE R66 ***O+HCL = OH+CL (Yung & DeMore, 1982)
          r43(k,i,lat) = 1.0e-12*exp(-3700./tn(k,i))    ! CLO+CO = CL+CO2 (Yung & DeMore, 1982)
!         r44(k,i,lat) = 7.5e-11*exp(-120./tn(k,i))     ! CLO+O = CL+O2 (Yung & DeMore, 1982)
!         r44(k,i,lat) = 2.5e-11*exp(110./tn(k,i))      ! CLO+O = CL+O2 (Saunder et al. 2006)
          r44(k,i,lat) = 3.0e-11*exp(70./tn(k,i))       ! CLO+O = CL+O2 (Zhang et al. 2012)
!         r99(k,i,lat)  ----------------------------------------------------------------------------------
!         r99(k,i,lat) = 0.0                            ! CL+CO+M = CLCO+M (Zhang et al. 2012, r159)
!         r99(k,i,lat) = 1.3e-34*exp(1000./tn(k,i))     ! CL+CO+M = CLCO+M (Yung & DeMore, 1982)
!         r99(k,i,lat) = 2.6e-34*exp(810./tn(k,i))      ! CL+CO+M = CLCO+M (Niovich et al.
!         r99(k,i,lat) = 4.2e-33*(300./tn(k,i))**3.8    ! CL+CO+M = CLCO+M (Zhang et al. 2012, r159)
!         r99(k,i,lat) = 2.1e-33*(300./tn(k,i))**3.8    ! CL+CO+M = CLCO+M (Zhang et al. 2012, r159) 50%
!         r99(k,i,lat) = 1.05e-33*(300./tn(k,i))**3.8   ! CL+CO+M = CLCO+M (Zhang et al. 2012, r159) 25%
          r99(k,i,lat) = 4.2e-33*(300./tncap)**3.8      ! CL+CO+M = CLCO+M (Zhang et al. 2012, r159), TN>210K
!         r108(k,i,lat) = 4.2e-12*exp(-1370./tn(k,i))   ! CL2+O = CLO+CL (Yung & DeMore, 1982)
          r108(k,i,lat) = 7.4e-12*exp(-1650./tn(k,i))   ! CL2+O = CLO+CL (Zhang et al 2012)
          r111(k,i,lat) = 1.3e-33*(300./tn(k,i))**3.8   ! CL+N2+CO = CLCO+N2 (Zhang et al 2012)
          r123(k,i,lat) = 3.0e-11*exp(70./tn(k,i))      !clo+o->cl+o2
          r130(k,i,lat) = 1.0e-12*exp(-1590./tn(k,i))   !2clo->cl2+o2
!         r160(k,i,lat) --------------------------------------------------------------------------------
!         r160(k,i,lat) = 0.0                           ! CLCO+CO2 = CL+CO+CO2 (Zhang et al 2012)
!         r160(k,i,lat) = 1.6e-25*exp(4000./tn(k,i))    ! CLCO+CO2 = CL+CO+CO2 (Zhang et al 2012)
!         r160(k,i,lat) = 0.8e-25*exp(4000./tn(k,i))    ! CLCO+CO2 = CL+CO+CO2 (Zhang et al 2012) 50%
!         r160(k,i,lat) = 0.4e-25*exp(4000./tn(k,i))    ! CLCO+CO2 = CL+CO+CO2 (Zhang et al 2012) 25%
          r160(k,i,lat) = 1.6e-25*exp(4000./tncap)      ! CLCO+CO2 = CL+CO+CO2 (Zhang et al 2012) TN>210K
          r162(k,i,lat) = 6.5e-33*exp(-2180./tn(k,i))   ! O+2CO = CO2+CO (Zhang et al 2012)
          r163(k,i,lat) = 3.4e-33*exp(-2180./tn(k,i))   ! 2O+CO = CO2+O (Zhang et al 2012)
!
!***************************************************************
! O3 CHEM:(Chris Parkinson's chem code)
          k50(k,i,lat) = 5.9e-34*(300./tn(k,i))**2.4    !2o+o2->o3+o  (low pressure rate) 3-body
          k51(k,i,lat) = 5.9e-34*(300./tn(k,i))**2.4    !o+2O2->o3+o2a (low pressure rate) 3-body
          k52(k,i,lat) = 6.7e-34*(300./tn(k,i))**2.5    !o+o2+co->o3+co  (low pressure rate) 3-body
          k53(k,i,lat) = 3.0e-12*exp(-1500./tn(k,i))    !no+o3->no2+o2
          k54(k,i,lat) = 2.3e-11*exp(-200./tn(k,i))     !cl+o3->clo+o2  (Zhang et al, 2012 update)
          k56(k,i,lat) = 1.0e-12*exp(-4000./tn(k,i))    !clo+o3->oclo+o2
          r63(k,i,lat) = 5.95e-34*(300./tn(k,i))**2.3   !O+O2+N2->o3+n2  (Zhang et al, 2012 update)
!
!***************************************************************
! NO Loss Chem:
          k45(k,i,lat) = 1.0e-27*(1.0/tn(k,i)**1.5)     ! no+o+m
          k46(k,i,lat) = 3.0e-12*exp(-1500./tn(k,i))    ! no+o3
          k47(k,i,lat) = 3.5e-12*exp(250./tn(k,i))      ! no+ho2
!***************************************************************
! OH CHEM:	  
!         r41b(k,i,lat) = 6.3e-11*exp(-570./tn(k,i))    ! HO2+Cl->OH+ClO
          r41b(k,i,lat) = 4.1e-11*exp(-450./tn(k,i))    ! HO2+Cl->OH+ClO (Sander, 2002 and used in Zhang et al 2012)
!         r49(k,i,lat) = 2.9e-11*exp(200./tn(k,i))      ! HO2+O->OH+O2
          r49(k,i,lat) = 3.0e-11*exp(200./tn(k,i))      ! HO2+O->OH+O2
!         r53(k,i,lat) = 1.6e-12*exp(-940./tn(k,i))     ! OH+O3->HO2+O2
	  r53(k,i,lat) = 1.7e-12*exp(-940./tn(k,i))     ! OH+O3->HO2+O2
!         r54(k,i,lat) = 1.1e-14*exp(-500./tn(k,i))     ! HO2+O3->OH+2O2
	  r54(k,i,lat) = 1.0e-14*exp(-490./tn(k,i))     ! HO2+O3->OH+2O2
	  r66(k,i,lat) = 1.0e-11*exp(-3300./tn(k,i))    ! HCl+O->OH+Cl
!         r104(k,i,lat) = 9.9e-33*(tn(k,i)**6.5)*exp(-1460./tn(k,i))  ! H2+O->OH+H
	  r104(k,i,lat) = 8.5e-20*(tn(k,i)**2.7)*exp(-3160./tn(k,i))  ! H2+O->OH+H
	  r105(k,i,lat) = 4.2e-12*exp(-240./tn(k,i))    ! 2OH->H2O+O
          rk23(k,i,lat) = 1.40e-10*exp(-470./tn(k,i))   ! H+O3->OH+O2
          rk28(k,i,lat) = 2.2e-11*exp(120./tn(k,i))     ! OH+O->H+O2
!
!***************************************************************
! SOX CHEM:
!         s114(k,i,lat) = 1.30e-31*exp(940./tn(k,i))     ! Cl + SO2 + M
          s114(k,i,lat) = 1.30e-34*exp(940./tn(k,i))     ! Cl + SO2 + M   (Mills, 1998)
          s167(k,i,lat) = 1.60e-13*exp(-2280./tn(k,i))   ! SO + O2 ---> SO2 + O  (Atkinson etal 2004)
!         s167(k,i,lat) = 1.60e-13*exp(-2280./tn(k,i))   ! SO + O2 ---> SO2 + O  (Singleton and Cvetanovic, 1988)
          s243(k,i,lat) = 7.30e-21/tn(k,i)**5.0          ! Cl + SO + M --> OSCl + M (Mills et al 2007)
!         s243(k,i,lat) = 3.00e-33/(tn(k,i)/300.)**5.0   ! Cl + SO + M --> OSCl + M (Mills, 1998)
          s248(k,i,lat) = 1.00e-12*exp(-1700/tn(k,i))     ! 2SO -->  SO2 + S (Moses, 2002; Martinez and Herron, 1983)
          s258(k,i,lat) = 1.32e-31*exp(-1000./tn(k,i))   ! O + SO2 + M --> SO3 + M (Atkinson etal 2004)
          s286(k,i,lat) = 5.00e-11*exp(-600./tn(k,i))    ! O + OSCl --> SO2 + Cl  (Moses, 2002) **Not in Petrass (2013)
          s287(k,i,lat) = 2.00e-11*exp(-600./tn(k,i))    ! O + OSCl --> SO + ClO (Moses, 2002) **Not in Petrass (2013)
          s299(k,i,lat) = 5.00e-11*exp(-800./tn(k,i))    ! SO + ClSO2 --> SO2 + OSCl (Moses etal 2002) **Not in Petrass13)
!***************************************************************
		  	
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

      end subroutine chemrates_tdep
      end module chemrates_module
