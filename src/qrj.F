      module qrj_module
      use params_module,only: nlevp1,nlonp4,spval
      use addfld_module,only: addfld
      use chapman_module,only: sovmr,so2vmr,so3vmr ! (nlevp1) from 1-D KINETICS
      implicit none
!     -------------------------------------------------------------
!     Venus TGCM version : latest mods (09/27/06) S. W. BOUGHER
!     Venus TGCM version : latest mods (11/16/06) S. W. BOUGHER
!     Venus TGCM version : latest mods (03/12/08) S. W. BOUGHER
!     Venus TGCM version : latest mods (04/02/08) S. W. BOUGHER  (3.0x qjn2)
!     Venus TGCM version : latest mods (04/21/08) S. W. BOUGHER  (xuvfac 10x inc)
!     Venus TGCM version : latest mods (04/23/08) S. W. BOUGHER  (1.0x qjn2)
!     Venus TGCM version : latest mods (05/07/08) S. W. BOUGHER  (2.0x qjn2)
!     Venus TGCM version : latest mods (05/07/08) S. W. BOUGHER  (xuvfac default)
!     btf 6/13/13: Add JSO,JSO2a-b
!     swb 9/20/13: Add JSO,JSO2a-b,JSO3
!     (units of qtotal: erg/gm/sec)
!     (units of Js, Is: #/cm3//sec)
!     -- Descending wavelength order (2250-1050A):  bins = 1-31
!     -- Descending wavelength order (1050-1A)   :  bins = 32-75
!     -- lmax2 = 91 (16-MUV,23-SRC,7-SRB,1-LYA,37-EUV,7-XRAY). [3025-1A]
!     -- lmax = 75 (23-SRC,7-SRB,1-LYA,37-EUV,7-XRAY). [2250-1A]
!     -- l1-1 = 31 (23-SRC,7-SRB,1-LYA) [2250-1050A, extension #1]
!     -- l2-1 = 16 (16-MUV) [3025-2250A, extension #2]
!     -------------------------------------------------------------
!     -- X-Sections for O1(A), IO(3p)   from tables. 
!     -- X-Sections for CO(A), ICO(tot) from tables. 
!     -- X-Sections for N2(A), IN2(tot), IN2_N2P from tables. 
!     -- X-Sections for O2(A), IO2(tot), from tables; JO2 calc. 
!     -- X-Sections for CO2(A), ICO2T, ICO2I, ICO2IOP from tables. 
!                       JCO2 calc.
!     -- X-section for JSO, JSO2(a), JSO2(b), JSO3 from 1-D KINETICS code
!     -------------------------------------------------------------
!     Venus PE (secondary) ionzation rates are approx. by enhancement factors
!     -- 1.40x  ICO, IN2
!     -- 1.60x  IO
!     -- 1.50x  ICO2
!     Venus soft XRAY flux changes are tested for production of 2nd ion peak
!     -- ssflux.F tunable parameters (xuvfac option to enable)
!     -- ICO2 (above), plus SF=SF(400./wavelength) can be modeled this way
!     -------------------------------------------------------------
!     Venus TGCM version : latest mods (05/19/08) A. S. BRECHT  (PE data from Gerard)
!     Venus SRC extension for CO2, O2, SO2a, SO2b, SO, SO3: (09/20/13) Bougher
!     -- l1-1 => 1-31 (new) 2250 A - 1750A
!     Venus MUV extension for SO3 : (09/XX/13) Bougher
!     -- l2-1 => 1-16 (new) 3025-2250A
!     -------------------------------------------------------------
!     10/13/09 Added T-dependent factor to CO2 x-sections A.S.B.
!     10/16/09 put in additional ionization factor for soft XRAY (see above) A.S.B.
!     -------------------------------------------------------------
!
! Coefficients for wavelength scans in qrj.
! These are referenced in routines qrj, init_sflux, init_sigmas, and 
! set_qrj_coeff (all in qrj_module.F).
!
      integer,parameter :: lmax2=91, lmax=75, l1=32, l2=17, nsig=8
! Venus-Sun heliocentric distance scaling factor (constant)
      real,parameter :: sfactor=1.914
!     -------------------------------------------------------------
      real,save :: 
     |  euveff(nlevp1),    !
     |  srceff(nlevp1),    !
     |  sigeuv(nsig,lmax), ! o,co,n2,co2,o2,so2a,so2b,so3
     |  rlmeuv(lmax),      !
     |  feuv(lmax),        !
     |  fsrc(l1-1),        !
     |  sigsrc_co2(l1-1),  !
     |  sigsrc_o2 (l1-1),  !
     |  sigsrc_so(l1-1),   !  JSO  --->  S + O
     |  sigsrc_so2a(l1-1), !  JSO2 --->  SO + O
     |  sigsrc_so2b(l1-1), !  JSO2 --->  S + O2
     |  sigsrc_so3(l1-1),  !  JSO3 --->  SO2 + O
     |  sigmuv_so3(l2-1),  !  JSO3 --->  SO2 + O
     |  rlmsrc(l1-1),      !
     |  sigmas(nsig,lmax), !
     |  sflux(lmax),       !
     |  sfluxext(l1-16),   !
     |  ph,                !
     |  wla(lmax),         !
     |  brn2(lmax),        !
     |  brso3tso2(lmax),   !
     |  brco2(lmax),       !
     |  brco(lmax),        !
     |  bro  (lmax)        !
!
! For euvac, init_sigmas (update from EUVAC)
      integer,parameter :: neuv=37
      real,dimension(neuv),save ::
     |  wleuv1,wleuv2,
     |  sigan2, sigin2,
     |  brn2np
!
! Heating, dissociation & ioniz terms (9) set by qrj, used by other routines.
! These are allocated for task subdomains by alloc_q (called from allocdata)
!
      real,dimension(:,:,:),allocatable,save ::
     |  rj,     ! solar photodissociation of JCO2  (#/cm3/sec)
     |  qjo2,   ! solar photodissociation of JO2   (#/cm3/sec)
     |  qjn2,   ! solar EUV + PE dissociation of JN2 (#/cm3/sec)
     |  qjso,   ! solar photodissociation of JSO (#/cm3/sec)
     |  qjso2a, ! solar photodissociation of JSO2a (#/cm3/sec)
     |  qjso2b, ! solar photodissociation of JSO2b (#/cm3/sec)
     |  qjso3,  ! solar photodissociation of JSO3 (#/cm3/sec)
     |  qfjso,   ! solar photodissociation of JSO (#/sec)
     |  qfjso2a, ! solar photodissociation of JSO2a (#/sec)
     |  qfjso2b, ! solar photodissociation of JSO2b (#/sec)
     |  qfjso3,  ! solar photodissociation of JSO3 (#/sec)
     |  qjn2pd,   ! solar EUV dissociation of JN2 (#/cm3/sec)
     |  qjn2pe,   !  PE dissociation of JN2 (#/cm3/sec)
     |  qtotal, ! total euv-uv heating (erg/gm/sec)
     |  qop,    ! O+ (#/cm3/sec)
     |  qcop,   ! CO+ (#/cm3/sec)
     |  qn2p,   ! n2+ (#/cm3/sec)
     |  qco2p,  ! CO2+ (#/cm3/sec)
     |  qo2p    ! O2+ (#/cm3/sec)
!
      contains
!-----------------------------------------------------------------------
      subroutine qrj(sco1,scco,scn2,sco2,scco2,scso,scso2,scso3,
     |  tn,o1,co,n2,o2,so,so2,barm,fcp,lev0,lev1,lon0,lon1,lat)
!
!    22-arguments to qrj call-- venus
!
! Calculate heating and dissociation rates.
!
      use input_module,only: f107, f107a
      use init_module,only: sfeps,idn
      use cons_module,only: t0,avo,rmassinv_o1,rmassinv_co,rmassinv_n2,
     |  rmassinv_o2,rmassinv_co2,rmassinv_so,rmassinv_so2,rmassinv_so3,
     |  rmass_o2,rmass_co2,rmass_so,rmass_so2,rmass_so3,gask,check_exp
      use lbc,only: fb,b
      use fields_module,only: tlbc
! ----------------------------------------------------------------------
! pepdgrid is photoelectron frequencies (/sec), read from pepdgrid_init
! Call it zjn2pepd here, and stuff pepdgrid into zjn2pepd 
!
!     use pepdgrid_init_module,only: zjn2pepd=>pepdgrid, ! full domain (nlon,nlat,nlevp1)
!    |                             pepdgrid2=>pepdgrid2
      use pepdgrid_init_module,only: zjn2pepd=>pepdgrid,zjn2pe=>pegrid, 
     |                               zjn2pd=>pdgrid,pepdgrid2=>pepdgrid2
! ----------------------------------------------------------------------
!
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) ::
     |  sco1,scco,scn2,sco2,scco2, ! chapman integrals
     |  scso,scso2,scso3,          ! new chapman integrals
     |  tn,o1,co,n2,o2,            ! tn and species mass mixing ratios
     |  so,so2,                    ! new species mass mixing ratios
     |  fcp,                       ! cp from cpktkm at interfaces
     |  barm                       ! mean mass
!
! VT vampir tracing:
!
#ifdef VT
#include "VT.inc"
#endif
!
! Local:
      integer :: k,i,l,nlevs,ier
! SRB Absorption now for CO2, O2, SO3, SO2a, SO2b, SO for Venus:
! PE enhancement factors (vtgcm2)
      real,parameter ::  sfn2p = 1.40, sfop = 1.60, sfco2p = 1.50
      real,parameter ::  hc = 1.9845E-16
      real :: rlmeuvinv(lmax),rlmsrcinv(l1-1),factor,fmin,fmax,
     |        wl400r(lmax)
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  r,      ! s14
     |  rp,     ! s15
     |  taur,   ! s1
     |  tauq,   ! s2
     |  tau1,   ! s3
     |  tau2,   ! s4
     |  tau3,   ! s5
     |  etaur,  ! s6
     |  etauq   ! s7
      real,dimension(l1-1,lev0:lev1,lon0:lon1) ::
     |  cxco2   ! correction factor for CO2 x-section at src
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  o1i,    ! o1 at interfaces (s1)
     |  coi,    ! co  at interfaces (s2)
     |  n2i,    ! n2  at interfaces (s3)
     |  o2i,    ! o2  at interfaces 
     |  co2i,   ! co2  at interfaces (s4)
     |  soi,    ! so mmr  at interfaces
     |  so2i,   ! so2 mmr at interfaces
     |  so3i,   ! so3  mmr at interfaces
     |  soden,  ! so density   at interfaces
     |  so2den, ! so2 density   at interfaces
     |  so3den, ! so3 density   at interfaces
     |  sovmr2, ! so vmr   at interfaces (from either 1-D or VTGCM)
     |  so2vmr2, ! so2 vmr   at interfaces (from either 1-D or VTGCM)
     |  so3vmr2, ! so3 vmr   at interfaces (always from 1-D KINETICS)
     |  tni,    ! tn  at interfaces (s6)
     |  TotMassDen, ! xnmbari converted to Total mass density [kg/m3]
     |  xnmbari, ! p0*e(-z)*barm/kT at interfaces
     |  xnmbar  ! p0*e(-z)*barm/kT at midlevels
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  zjco2,       ! co2  dissociation 
     |  zjo2,        ! o2   dissociation 
     |  zjso,        ! JSO   dissociation
     |  zjso2a,      ! JSO2a dissociation
     |  zjso2b,      ! JSO2b dissociation
     |  fjco2,       ! co2  dissociation 
     |  fjo2,        ! o2   dissociation 
     |  fjso,        ! JSO   dissociation
     |  fjso2a,      ! JSO2a dissociation
     |  fjso2b,      ! JSO2b dissociation
     |  fjso3,       ! JSO3 dissociation
     |  zo1,         ! o1   ionization 
     |  zco,         ! co   ionization 
     |  zo2_tot,     ! o2   total ionization
     |  zn2_tot,     ! n2   total ionization
     |  zn2_brn2p,   ! n2   ionization to N2+ 
     |  zco2_tot,    ! co2  total ionization 
     |  zco2_brco2p, ! co2  ionization to CO2+ 
     |  zco2_brop,   ! co2  ionization to O+ 
     |  zco2_brcop,  ! co2  ionization to CO+ 
     |  zso3_brso2,  ! so3  dissociation to so2 
     |  sigchap,   ! 
     |  fexp,      ! 
     |  qeuv,      ! QEUV (erg/cm3/sec)
     |  qsrc,      ! QSRC (erg/cm3/sec)
     |  qtot       ! QTOT (erg/cm3/sec)
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  sum1,   ! sum(o1,co,n2,o2,co2,so2a,so2b,so3)(Fluxes)   
     |  sum2,   ! sum(o1,co,n2,o2,co2)((sigma*psi/rmass)
     |  sum3,   ! sum(co2)(sigma*Fluxes)   
     |  sum4,   ! sum(co2)(sigma*Fluxes*energy) 
     |  sum5,   ! sum(o2)(sigma*Fluxes) 
     |  sum6,   ! sum(so2a)(sigma*Fluxes)
     |  sum7,   ! sum(so2b)(sigma*Fluxes)
     |  sum8,   ! sum(so)(sigma*Fluxes)
     |  sum9    ! sum(so3)(sigma*Fluxes)

      logical,parameter :: 
     |  debug=.false.     ! insert print statements
!
! expo() (util.F) is used only if check_exp is true. This will avoid
! NaNS fpe, but will degrade performance. Check_exp is in cons.F.
!
      real,external :: expo ! used only when check_exp is set (util.F)
!
      if (debug) write(6,"('Enter qrj: lat=',i3,' lon0,1=',2i3)") 
     |  lat,lon0,lon1
!
#ifdef VT
      code = 118 ; state = 'qrj' ; activity='ModelCode'
      call vtbegin(118,ier)
#endif
!
! Exec:
      nlevs = lev1-lev0+1

!     write(6,"('qrj: lat=',i2,' lon0,1=',2i3,' idn(lon0:lon1)=',/,
!    |  (15i3))") lat,lon0,lon1,idn(lon0:lon1)
!
! --------------------------------------------------------------------
!     call addfld('SCO1',' ',' ',sco1(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('SCCO',' ',' ',scco(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('SCN2',' ',' ',scn2(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('SCO2',' ',' ',sco2(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('SCCO2',' ',' ',scCO2(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('SCSO',' ',' ',scSO(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('SCSO2',' ',' ',scSO2(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
! --------------------------------------------------------------------
!
! Depedence factor of CO2 ionization
      do i=1,lmax
         IF (rlmeuv(i) .LT. 0.4E-5) THEN
            wl400r(i) = 0.4E-5/rlmeuv(i)
!           wl400r(i) = 1.
         END IF
         IF (rlmeuv(i) .GE. 0.4E-5) THEN
            wl400r(i) = 1.
         END IF
      enddo
!
! Will multiply by inverses:
      do i=1,lmax
        rlmeuvinv(i) = 1./rlmeuv(i)
      enddo
      do i=1,l1-1 ! 1,15
        rlmsrcinv(i) = 1./rlmsrc(i)
      enddo
!
! COMPUTE S14 = RSQ, S15 = RSP
! S1 = TAU(R),  S2 = TAU(Q)
!
!     write(6,"('qrj: lat=',i2,' lon0,1=',2i3,' idn(lon0:lon1)=',/,
!    |  (12i3))") lat,lon0,lon1,idn(lon0:lon1)
!
! ---------------------------------------------------------------
! Optical Depth Calculations: Is this a PE enhancement scheme?
! Assess Carefully Before Using. Unsure for Venus!!
! ---------------------------------------------------------------
      do i=lon0,lon1
        do k=lev0,lev1
          taur(k,i) = sigeuv(1,49)*sco1(k,i)+sigeuv(2,49)*scco(k,i)+
     |                sigeuv(3,49)*scn2(k,i)
          tauq(k,i) = sigeuv(1,20)*sco1(k,i)+sigeuv(2,20)*scco(k,i)+
     |                sigeuv(3,20)*scn2(k,i)
          tau1(k,i) = 1.3*taur(k,i)
          tau2(k,i) = 2.0*taur(k,i)
          tau3(k,i) = 2.5*taur(k,i)
          if (taur(k,i) > 9.) taur(k,i) = 9.
          if (tauq(k,i) > 9.) tauq(k,i) = 9.
          if (tau1(k,i) > 9.) tau1(k,i) = 9.
          if (tau2(k,i) > 9.) tau2(k,i) = 9.
          if (tau3(k,i) > 9.) tau3(k,i) = 9.
!
! TAU(N) = EXP(-TAU(N)) FOR N = 1,3,1
!
          tau1(k,i) = exp(-tau1(k,i))
          tau2(k,i) = exp(-tau2(k,i))
          tau3(k,i) = exp(-tau3(k,i))
!
! taur = EXP(-TAU(R)), tauq = EXP(-TAU(Q))
!
          etaur(k,i) = exp(-taur(k,i))
          etauq(k,i) = exp(-tauq(k,i))
!
! r = r(n2p) :  Are these enchancements for PE yields?
!
          r(k,i) = etaur(k,i)+2.*(tau1(k,i)+tau2(k,i)+tau3(k,i))
          r(k,i) = 2.4*etaur(k,i)/r(k,i)
        enddo
      enddo
      if (debug) write(6,"(/,'qrj after tau: lat=',i3)") lat
!
!     call addfld('TAUQ',' ',' ',tauq,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('TAUR',' ',' ',taur,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('R'   ,' ',' ',r   ,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! ------------------------------------------------------------
! 9/25/06 btf: bottom boundary of TN is in tlbc:
      do i=lon0,lon1
        tni(lev0,i) = tlbc(i,lat)
        do k=lev0+1,lev1-1
          tni(k,i) = .5*(tn(k-1,i)+tn(k,i))
        enddo
        tni(lev1,i) = tn(lev1-1,i) ! nlevp1 <- nlev
      enddo

!     write(6,"('qrj: lat=',i3,' tni(lev0,:)=',/,(6e12.4))") 
!    |  lat,tni(lev0,:)

!     call addfld('TNI',' ',' ',tni,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('BARM',' ',' ',barm(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Calculate p0*e(-z)*barm/kT (1/(N*MBAR))at interfaces
      call mk_xnmbari(tni,barm(:,lon0:lon1),
     |   xnmbari,
     |   lev0,lev1,lon0,lon1,lat)

      call addfld('XNMBARI',' ',' ',xnmbari,
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! ---------------------------------------------------------------
!
! o1,co,n2,o2 at interface levels: mass mixing ratios
      do i=lon0,lon1
        do k=lev0,lev1-1
          o1i (k+1,i) = 0.5*(o1(k,i)+o1(k+1,i))
          coi (k+1,i) = 0.5*(co(k,i)+co(k+1,i))
          n2i (k+1,i) = 0.5*(n2(k,i)+n2(k+1,i))
          o2i (k+1,i) = 0.5*(o2(k,i)+o2(k+1,i))
! ----------------------------------------------------------------
! Full-up integration of so and so2 into calculation  (mmr vs vmr)
! MMR prognostics (2)
! soi and so2i must be positive!
          soi (k+1,i)=0.5*(so(k,i)+so(k+1,i))
          so2i(k+1,i)=0.5*(so2(k,i)+so2(k+1,i))
          so3i(k+1,i)=0.5*(so3vmr(k)+so3vmr(k+1))*(rmass_so3/barm(k,i))
! VMR for Diagnostics (3)
          sovmr2 (k+1,i)=soi(k+1,i)*(barm(k,i)*rmassinv_so)
          so2vmr2(k+1,i)=so2i(k+1,i)*(barm(k,i)*rmassinv_so2)
          so3vmr2(k+1,i)=0.5*(so3vmr(k)+so3vmr(k+1))
! Densities for Diagnostics
          soden(k+1,i)= soi(k+1,i)*xnmbari(k+1,i)*rmassinv_so
          so2den(k+1,i)= so2i(k+1,i)*xnmbari(k+1,i)*rmassinv_so2
          so3den(k+1,i)= so3i(k+1,i)*xnmbari(k+1,i)*rmassinv_so3
! ----------------------------------------------------------------
! Testing with 1-D KINETICS prescription of so and so2 (MMR and DEN and VMR)
!         soi (k+1,i)=0.5*(sovmr(k)+sovmr(k+1))*(rmass_so/barm(k,i))
!         so2i(k+1,i)=0.5*(so2vmr(k)+so2vmr(k+1))*(rmass_so2/barm(k,i))
!         so3i(k+1,i)=0.5*(so3vmr(k)+so3vmr(k+1))*(rmass_so3/barm(k,i))
!         sovmr2(k+1,i)=0.5*(sovmr(k)+sovmr(k+1))
!         so2vmr2(k+1,i)=0.5*(so2vmr(k)+so2vmr(k+1))
!         so3vmr2(k+1,i)=0.5*(so3vmr(k)+so3vmr(k+1))
! ----------------------------------------------------------------
        enddo
      enddo
! ----------------------------------------------------------------
!
! Bottom boundary: MMR and VMR  and DENSITIES
      do i=lon0,lon1
! MMR
        o1i(lev0,i)= .5*((b(i,1,1)+1.)*o1(1,i)+b(i,1,2)*co(1,i)+fb(i,1))
        coi(lev0,i)= .5*(b(i,2,1)*o1(1,i)+(b(i,2,2)+1.)*co(1,i)+fb(i,2))
        n2i(lev0,i)=  n2i(lev0+1,i)
        o2i(lev0,i)=  o2i(lev0+1,i)
        soi(lev0,i)=  soi(lev0+1,i)
        so2i(lev0,i)= so2i(lev0+1,i)
        so3i(lev0,i)= so3i(lev0+1,i)
! VMR
        sovmr2(lev0,i)=  sovmr2(lev0+1,i)
        so2vmr2(lev0,i)= so2vmr2(lev0+1,i)
        so3vmr2(lev0,i)= so3vmr2(lev0+1,i)
!  DENSITIES
        soden(lev0,i)=  soden(lev0+1,i)
        so2den(lev0,i)= so2den(lev0+1,i)
        so3den(lev0,i)= so3den(lev0+1,i)
      enddo
! ----------------------------------------------------------------

!     write(6,"('qrj: lat=',i4,' so2i min,max=',2e12.4)") 
!    |  lat,minval(so2i),maxval(so2i)
!
! CO2 at interfaces:  mass mixing ratios (1-PSIO-PSICO-PSIN2)
      do k=lev0,lev1
        co2i(k,:) = 1.-o1i(k,:)-coi(k,:)-n2i(k,:)
      enddo 
! ----------------------------------------------------------------
!  Checks:
!     call addfld('O1I',' ',' ',o1i,'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('COI',' ',' ',coi,'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('N2I',' ',' ',n2i,'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('O2I',' ',' ',o2i,'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('CO2I',' ',' ',co2i,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!   MMR SOx
      call addfld('SOMMR',' ',' ',soi,
     |  'ilev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('SO2MMR',' ',' ',so2i,
     |  'ilev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('SO3MMR',' ',' ',so3i,
     |  'ilev',lev0,lev1,'lon',lon0,lon1,lat)
!   VMR SOx
      call addfld('SOVMR',' ',' ',sovmr2,
     |  'ilev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('SO2VMR',' ',' ',so2vmr2,
     |  'ilev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('SO3VMR',' ',' ',so3vmr2,
     |  'ilev',lev0,lev1,'lon',lon0,lon1,lat)
!   DEN SOx
      call addfld('SODEN',' ',' ',soden,
     |  'ilev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('SO2DEN',' ',' ',so2den,
     |  'ilev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('SO3DEN',' ',' ',so3den,
     |  'ilev',lev0,lev1,'lon',lon0,lon1,lat)
!
      if (debug) write(6,"(/,'qrj after interfaces: lat=',i3)") lat
! 
! ---------------------------------------------------------------------
! Initialize heating, I, J terms on current processor domain: 8-fields
! (global module data for use by other routines)
      rj    (lev0:lev1,lon0:lon1,lat) = 0.
      qjo2  (lev0:lev1,lon0:lon1,lat) = 0.  
      qjn2  (lev0:lev1,lon0:lon1,lat) = 0.  
      qjso  (lev0:lev1,lon0:lon1,lat) = 0.
      qjso2a (lev0:lev1,lon0:lon1,lat) = 0.
      qjso2b (lev0:lev1,lon0:lon1,lat) = 0.
      qjso3 (lev0:lev1,lon0:lon1,lat) = 0.
      qjn2pe  (lev0:lev1,lon0:lon1,lat) = 0.  
      qjn2pd  (lev0:lev1,lon0:lon1,lat) = 0.  
      qtotal(lev0:lev1,lon0:lon1,lat) = 0.
      qop   (lev0:lev1,lon0:lon1,lat) = 0.
      qcop  (lev0:lev1,lon0:lon1,lat) = 0.
      qn2p  (lev0:lev1,lon0:lon1,lat) = 0.
      qco2p (lev0:lev1,lon0:lon1,lat) = 0.
      qo2p  (lev0:lev1,lon0:lon1,lat) = 0.
      if (debug) write(6,"(/,'qrj after 1st init: lat=',i3)") lat
!
! Initialize ionization, dissociation ratesm working arrays: 12-fields
! (array ops on local arrays dimensioned (lev0:lev1,lon0:lon1)):
      zjco2 = 0.        ! co2  dissociation 
      zjo2  = 0.        ! o2  dissociation 
      zjso   = 0.       ! so  dissociation
      zjso2a = 0.       ! so2a dissociation
      zjso2b = 0.       ! so2b dissociation
      zo1   = 0.        ! o   ionization 
      zco   = 0.        ! co  ionization 
      zo2_tot     = 0.  ! o2   total ionization
      zn2_tot     = 0.  ! n2   total ionization
      zn2_brn2p   = 0.  ! n2  ionization to N2+ 
      zco2_tot    = 0.  ! co2  total ionization 
      zco2_brco2p = 0.  ! co2  ionization to CO2+ 
      zco2_brop   = 0.  ! co2  ionization to O+ 
      zco2_brcop  = 0.  ! co2  ionization to CO+ 
      zso3_brso2  = 0.  ! so3  dissociation to so2 + O (major channel) 
      sigchap     = 0.  ! SRC intermediate product
      fexp        = 0.  ! SRC intermediate product

      if (debug) write(6,"(/,'qrj after initializations: lat=',i3)") lat
!
! ---------------------------------------------------------------------
!
! Wavelength loop --------------------------------------------------
! Summation over wavelength: soft-XRAY and EUV bins (44): 1050-1A
! Multiply Q by (1/r**2) mean heliocentric distance of Venus 
!
      do l=l1,lmax     ! 32,75 (44 wavelengths;37-EUV,7-XRAY) [1050A-1A]
!       if (debug) write(6,"('qrj l-loop: l=',i3,' l1=',i3,' lmax=',
!    |    i3)") l,l1,lmax
        sum1(:,:) = 0. ! sum(o1,co,n2,o2,co2,so2a,so2b,so3)(Fluxes). FLUXES.
        sum2(:,:) = 0. ! sum(o1,co,n2,o2,co2,so2,so3)(sigma*psi/rmass). DENSITIES.
        do i=lon0,lon1
          do k=lev0,lev1
            sum1(k,i) = sum1(k,i)+sigeuv(1,l)*sco1(k,i)+
     |                            sigeuv(2,l)*scco(k,i)+
     |                            sigeuv(3,l)*scn2(k,i)+
     |                            sigeuv(4,l)*scco2(k,i)+
     |                            sigeuv(5,l)*sco2(k,i)+
     |                            sigeuv(6,l)*scso2(k,i)+
     |                            sigeuv(7,l)*scso2(k,i)+
     |                            sigeuv(8,l)*scso3(k,i)

            if (.not.check_exp) then
              sum1(k,i) = feuv(l)*exp(-sum1(k,i))*sfactor
            else
              sum1(k,i) = feuv(l)*expo(-sum1(k,i),0)*sfactor
            endif

            sum2(k,i) = sum2(k,i)+sigeuv(1,l)*o1i(k,i)*rmassinv_o1+
     |                  sigeuv(2,l)*coi(k,i)*rmassinv_co+
     |                  sigeuv(3,l)*n2i(k,i)*rmassinv_n2+
     |                  sigeuv(4,l)*co2i(k,i)*rmassinv_co2+
     |                  sigeuv(5,l)*o2i(k,i)*rmassinv_o2+
     |                (sigeuv(6,l)+sigeuv(7,l))*so2i(k,i)*rmassinv_so2+
     |                  sigeuv(8,l)*so3i(k,i)*rmassinv_so3
          enddo
        enddo

!       write(6,"('qrj: lat=',i4,' l=',i4,' sum1 mnmx=',2e12.4)")
!    |    lat,l,minval(sum1),maxval(sum1)

!       if (debug) write(6,"('qrj: l=',i3,' after loop 1')") l
! -------------------------------------------------------------------
! Longitude and column domain of current process:
! Longitude-Column Loop
!    qtotal,      ! total EUV-UV heating
!    zjco2,       ! co2  dissociation 
!    zjo2,        ! o2   dissociation 
!    zjn2,        ! n2   dissociation (solar EUV) 
!    zjso         ! so  dissociation  (solar SRC/SRB)
!    zjso3        ! so3 dissociation  (solar SRC/SRB/MUV)
!    zjso2a       ! so2a dissociation  (solar EUV/SRC/SRB)
!    zjso2b       ! so2b dissociation (solar EUV/SRC/SRB)
!    zo1,         ! o1   ionization 
!    zco,         ! co   ionization  (including PE)
!    zo2_tot,     ! o2   total ionization
!    zn2_tot,     ! n2   total ionization (including PE)
!    zn2_brn2p,   ! n2   ionization to N2+ (including PE) 
!    zco2_tot,    ! co2  total ionization (including PE) 
!    zco2_brco2p, ! co2  ionization to CO2+ (including PE)
!    zco2_brop,   ! co2  ionization to O+ (including PE) 
!    zco2_brcop,  ! co2  ionization to CO+ (including PE) 
!    zso3_brso2   ! so3  dissociation to so2 + O (major channel) 
! -------------------------------------------------------------------
! ** wla (A units) goes from 2225 A - 1075 A
!
      wla = (/
     |        0.22500E+04, 0.22000E+04, 0.21500E+04, 0.21000E+04,
     |        0.20500E+04, 0.20000E+04, 0.19500E+04, 0.19200E+04,
     |        0.18900E+04, 0.18700E+04, 0.18500E+04, 0.18300E+04,
     |        0.18100E+04, 0.17900E+04, 0.17700E+04, 0.17500E+04,
     |        0.17250E+04, 0.16750E+04, 0.16250E+04, 0.15750E+04,
     |        0.15250E+04, 0.14750E+04, 0.14250E+04, 0.13750E+04,
     |        0.13250E+04, 0.12750E+04, 0.12250E+04, 0.12157E+04,
     |        0.11750E+04, 0.11250E+04, 0.10750E+04, 0.10250E+04,
     |        0.10319E+04, 0.10257E+04, 0.97500E+03, 0.97702E+03,
     |        0.92500E+03, 0.87500E+03, 0.82500E+03, 0.77500E+03,
     |        0.78936E+03, 0.77041E+03, 0.76515E+03, 0.72500E+03,
     |        0.70331E+03, 0.67500E+03, 0.62500E+03, 0.62973E+03,
     |        0.60976E+03, 0.57500E+03, 0.58433E+03, 0.55437E+03,
     |        0.52500E+03, 0.47500E+03, 0.46522E+03, 0.42500E+03,
     |        0.37500E+03, 0.36807E+03, 0.32500E+03, 0.30378E+03,
     |        0.30331E+03, 0.27500E+03, 0.28415E+03, 0.25630E+03,
     |        0.22500E+03, 0.17500E+03, 0.12500E+03, 0.75000E+02,
     |        0.41000E+02, 0.27500E+02, 0.19500E+02, 0.12000E+02,
     |        0.60000E+01, 0.30000E+01, 0.15000E+01/)
!
        do i=lon0,lon1
          do k=lev0,lev1
             qtotal(k,i,lat) = qtotal(k,i,lat)+hc*rlmeuvinv(l)*
     |         sum1(k,i)*sum2(k,i)
!
! Add ionization for o1, co, n2, o2, co2: total values (only)!
! Enhancement factors (PE) included for net impact on total ionization
             zo1(k,i) = zo1(k,i)+sigmas(1,l)*sum1(k,i)*
     |          o1i(k,i)*rmassinv_o1*sfop
             zco(k,i) = zco(k,i)+sigmas(2,l)*sum1(k,i)*
     |          coi(k,i)*rmassinv_co*sfn2p
             zn2_tot(k,i) = zn2_tot(k,i)+sigmas(3,l)*sum1(k,i)*
     |          n2i(k,i)*rmassinv_n2*sfn2p
             zo2_tot(k,i) = zo2_tot(k,i)+sigmas(5,l)*sum1(k,i)*
     |         o2i(k,i)*rmassinv_o2
             zco2_tot(k,i) = zco2_tot(k,i)+sigmas(4,l)*sum1(k,i)*
     |          co2i(k,i)*rmassinv_co2*sfco2p*wl400r(l)
!
! Photo-dissoc of o2: Must add SRC component also (using sum3 & sum4)!
            zjo2(k,i) = zjo2(k,i)+(sigeuv(5,l)-sigmas(5,l))*
     |        sum1(k,i)*o2i(k,i)*rmassinv_o2
!
! Photo-dissoc of n2: For NO* nightglow calculation
!           zjn2(k,i) = zjn2(k,i)+(sigeuv(3,l)-sigmas(3,l))*
!    |        sum1(k,i)*n2i(k,i)*rmassinv_n2
!
! Photo-dissoc of co2: Must add SRC component also (using sum3 & sum4)!
            zjco2(k,i) = zjco2(k,i)+(sigeuv(4,l)-sigmas(4,l))*
     |        sum1(k,i)*co2i(k,i)*rmassinv_co2
!
! Photo-dissoc of JSO2a: Must add SRC component also (using sum6)!

            zjso2a(k,i) = zjso2a(k,i)+sigeuv(6,l)*
     |        sum1(k,i)*so2i(k,i)*rmassinv_so2
!
! Photo-dissoc of JSO2b: Must add SRC component also (using sum7)!
            zjso2b(k,i) = zjso2b(k,i)+sigeuv(7,l)*
     |        sum1(k,i)*so2i(k,i)*rmassinv_so2
!
!           if (zjso2a(k,i) /= 0.) write(6,"('qrj: k,i,lat=',3i4,
!    |        ' zjso2a(k,i) /= 0: ',e12.4)") k,i,lat,zjso2a(k,i)
!           if (zjso2b(k,i) /= 0.) write(6,"('qrj: k,i,lat=',3i4,
!    |        ' zjso2b(k,i) /= 0: ',e12.4)") k,i,lat,zjso2b(k,i)

!
! Partitioning Ionization for  brn2p, brco2p, brop:
! Enhancement factors (PE) included for net impact on ionization paths
!  ** IN2 path to N2+ ionization ----------------------------
            zn2_brn2p(k,i) = zn2_brn2p(k,i)+
     |          sigmas(3,l)*sum1(k,i)*n2i(k,i)*rmassinv_n2*
     |          brn2(l)*sfn2p
!  ** ICO2 path to CO2+ ionization ----------------------------
            zco2_brco2p(k,i) =  zco2_brco2p(k,i)+
     |          sigmas(4,l)*sum1(k,i)*co2i(k,i)*rmassinv_co2*
     |          brco2(l)*sfco2p*wl400r(l)
!  ** ICO2 path to O+ ionization ----------------------------
            zco2_brop(k,i) =  zco2_brop(k,i)+
     |          sigmas(4,l)*sum1(k,i)*co2i(k,i)*rmassinv_co2*
     |          bro(l)*sfco2p*wl400r(l)
!  ** ICO2 path to CO+ ionization ----------------------------
            zco2_brcop(k,i) =  zco2_brcop(k,i)+
     |          sigmas(4,l)*sum1(k,i)*co2i(k,i)*rmassinv_co2*
     |          brco(l)*sfco2p*wl400r(l)
!  ** JSO3 path to SO2 + O+ dissociation ----------------------------
            zso3_brso2(k,i) =  zso3_brso2(k,i)+
     |          sigeuv(8,l)*sum1(k,i)*so3i(k,i)*rmassinv_so3*
     |          brso3tso2(l)
            enddo
          enddo  
!       if (debug) write(6,"('qrj end l-loop: l=',i3)") l
!
! End wavelength loop1
      enddo ! l=l1,lmax

!     write(6,"('qrj after l1,lmax: lat=',i4,' zjso2a=',2e12.4,
!    |  ' zjso2b=',2e12.4)") lat,
!    |  minval(zjso2a),maxval(zjso2a),minval(zjso2b),maxval(zjso2b)

      if (debug) write(6,"(/,'qrj after wavelength loop1: lat=',i3)") 
     | lat

!     call addfld('QTOTAL1',' ',' ',qtotal(:,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Multiply Q by EUV efficiency factor: Venus value from Fox (1988)
! -- Check Against vtgcm2 3-D Model Heating (qeuv) : erg/cm3/sec
! euveff(:)==0.18 (see init_qrj)
!
      do i=lon0,lon1
        do k=lev0,lev1
          qtotal(k,i,lat) = qtotal(k,i,lat)*euveff(k)*avo
          qeuv(k,i)=qtotal(k,i,lat)*xnmbari(k,i)/avo
        enddo
      enddo
!
! 1-1050A Neutral Heating at Venus (erg/cm3/sec): Scaled Incl. Efficiency
!     call addfld('QEUV',' ',' ',qeuv,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! 1-1050A Neutral Heating at Venus (erg/g/sec): Scaled Incl. Efficiency
!     call addfld('QTOTAL2',' ',' ',qtotal(:,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
!  --------------------------------------------------------------
! Final Contributions to IO+, ICO+, IN2+, ICO2+, JN2, JCO2, JO2:
!  --------------------------------------------------------------
!   **Working Array Fields to 11-slots**
!    qtotal,      ! total EUV-UV heating
!    zjco2,       ! co2  dissociation 
!    zjo2,        ! o2   dissociation 
!    zjn2,        ! n2   dissociation  (solar EUV) 
!    zjso         ! so  dissociation  (solar SRC/SRB)
!    zjso3        ! so3 dissociation  (solar SRC/SRB/MUV)
!    zjso2a       ! so2a dissociation  (solar EUV/SRC/SRB)
!    zjso2b       ! so2b dissociation (solar EUV/SRC/SRB)
!    zo1,         ! o1   ionization 
!    zco,         ! co   ionization 
!    zo2_tot,     ! o2   total ionization
!    zn2_tot,     ! n2   total ionization
!    zn2_brn2p,   ! n2   ionization to N2+ 
!    zco2_tot,    ! co2  total ionization 
!    zco2_brco2p, ! co2  ionization to CO2+ 
!    zco2_brop,   ! co2  ionization to O+ 
!    zco2_brcop,  ! co2  ionization to CO+ 
!    zso3_brso2   ! so3  dissociation to so2 + O (major channel) 
!  --------------------------------------------------------------
!   **Fields to 12-allocation slots**
!    rj    (lev0:lev1,lon0:lon1,lat) = 0.
!    qjo2  (lev0:lev1,lon0:lon1,lat) = 0.  
!    qjn2  (lev0:lev1,lon0:lon1,lat) = 0.  
!    qjso  (lev0:lev1,lon0:lon1,lat) = 0.
!    qjso2a (lev0:lev1,lon0:lon1,lat) = 0.
!    qjso2b (lev0:lev1,lon0:lon1,lat) = 0.
!    qjso3 (lev0:lev1,lon0:lon1,lat) = 0.
!    qtotal(lev0:lev1,lon0:lon1,lat) = 0.
!    qop   (lev0:lev1,lon0:lon1,lat) = 0.
!    qcop  (lev0:lev1,lon0:lon1,lat) = 0.
!    qn2p  (lev0:lev1,lon0:lon1,lat) = 0.
!    qco2p (lev0:lev1,lon0:lon1,lat) = 0.
!    qo2p  (lev0:lev1,lon0:lon1,lat) = 0.
! ---------------------------------------------------------------------
!     WRITE(*,*)'zjn2pepd', zjn2pepd
      do i=lon0,lon1
        do k=lev0,lev1
!  JCO2 Dissociation Rate (#/cm3/sec) 
          rj(k,i,lat) = zjco2(k,i)*xnmbari(k,i)
!  JO2 Dissociation Rate (#/cm3/sec) 
          qjo2(k,i,lat) = zjo2(k,i)*xnmbari(k,i)
!  JSO2a Dissociation Rate (#/cm3/sec)
          qjso2a(k,i,lat) = zjso2a(k,i)*xnmbari(k,i)
!  JSO2b Dissociation Rate (#/cm3/sec)
          qjso2b(k,i,lat) = zjso2b(k,i)*xnmbari(k,i)
!  JSO3 Dissociation Rate (#/cm3/sec)
          qjso3(k,i,lat) = zso3_brso2(k,i)*xnmbari(k,i)
!
! ------JN2 ------
!  JN2 Dissociation Rate (#/cm3/sec): both solar EUV + PE 
!         qjn2(k,i,lat) = zjn2(k,i)*xnmbari(k,i) + 
!    |                    zjn2pepd(i,lat,k)*
!    |                    xnmbari(k,i)*n2i(k,i)*rmassinv_n2
!  JN2 Dissociation Rate (#/cm3/sec): both solar EUV + PE 
!         qjn2(k,i,lat) = zjn2pepd(i,lat,k)*
!    |                    xnmbari(k,i)*n2i(k,i)*rmassinv_n2
!  JN2 Dissociation Rate (#/cm3/sec): both solar EUV + PE*3 
          qjn2(k,i,lat) = zjn2pepd(i,lat,k)*3.0*
     |                    xnmbari(k,i)*n2i(k,i)*rmassinv_n2
!  JN2 Dissociation Rate (#/cm3/sec): both solar(EUV + PE)*4 
!         qjn2(k,i,lat) = zjn2pepd(i,lat,k)*4.0*
!    |                    xnmbari(k,i)*n2i(k,i)*rmassinv_n2
!  JN2 Dissociation Rate (#/cm3/sec): both solar(EUV + PE)*2 
!         qjn2(k,i,lat) = zjn2pepd(i,lat,k)*2.0*
!    |                    xnmbari(k,i)*n2i(k,i)*rmassinv_n2
!  JN2 Dissociation Rate (#/cm3/sec): solar(EUV)*2 
          qjn2pd(k,i,lat) = zjn2pd(i,lat,k)*
     |                    xnmbari(k,i)*n2i(k,i)*rmassinv_n2
!  JN2 Dissociation Rate (#/cm3/sec): solar(PE)*2 
          qjn2pe(k,i,lat) = zjn2pe(i,lat,k)*
     |                    xnmbari(k,i)*n2i(k,i)*rmassinv_n2
! Seperate variable to see individual term
!         qjn22(k,i,lat) = zjn2(k,i)*xnmbari(k,i)  
!         qjn2pepd2(k,i,lat) = zjn2pepd(i,lat,k)*
!    |                       xnmbari(k,i)*n2i(k,i)*rmassinv_n2
!  JN2 Dissociation Rate (#/cm3/sec): solar EUV only
!         qjn2(k,i,lat) = 2.*zjn2(k,i)*xnmbari(k,i)
!
! ------JN2 ------
!  IO Ionization Rate (#/cm3/sec) 
          qop(k,i,lat) = zo1(k,i)*xnmbari(k,i)+
     |                   zco2_brop(k,i)*xnmbari(k,i)
!  ICO Ionization Rate (#/cm3/sec) 
          qcop(k,i,lat) = zco(k,i)*xnmbari(k,i)+
     |                    zco2_brcop(k,i)*xnmbari(k,i)
!  IN2 to N2+ Ionization rate
          qn2p(k,i,lat) =  zn2_brn2p(k,i)*xnmbari(k,i)
!  ICO2 to CO2+ Ionization rate
          qco2p(k,i,lat) =  zco2_brco2p(k,i)*xnmbari(k,i)
!  IO2 to O2+ Ionization rate
          qo2p(k,i,lat) =  zo2_tot(k,i)*xnmbari(k,i)
! -------------------------------------------------------------------
        enddo
      enddo

!     write(6,"('qrj after l1,lmax: lat=',i4,' qjso2a=',2e12.4,
!    |  ' qjso2b=',2e12.4)") lat,
!    |  minval(qjso2a),maxval(qjso2a),
!    |  minval(qjso2b),maxval(qjso2b)

      if (debug) write(6,"(/,'qrj after dissoc and ion rates: ',
     |  'lat=',i3)") lat
!
!
!---------------------------------------------------------------
! Wavelength loop : Extension #1:   l1-1 = 31 bins
! Summation over wave length: 2225-1050A (SRB, LYA, SRC Start)
!---------------------------------------------------------------
! Fluxes at Venus are scaled by mean Venus heliocentric distance
! Multiply Q by (1/r**2) mean heliocentric distance of Venus 
! Double Check this Scheme before proceeding to 1-D validation
! -- sigsrc_co2= co2 x-section (longward of 1050A)
! -- sigsrc_o2 = o2 x-section  (longward of 1050A)
! -- sigsrc_so2a = so2a x-section  (longward of 1050A)
! -- sigsrc_so2b = so2b x-section  (longward of 1050A)
! -- sigsrc_so3 = so3 x-section  (longward of 1050A)
!  Initializations outside lambda loop for full SRC region summation
!
      sum3(:,:) = 0. ! sum(co2)(sigma*F).    SIGMA*FLUXES.
      sum4(:,:) = 0. ! sum(co2)(sigma*F*hv). SIGMA*FLUXES*ENERGY
      sum5(:,:) = 0. ! sum(o2)(sigma*F).     SIGMA*FLUXES.
      sum6(:,:) = 0. ! sum(so2a)(sigma*F).   SIGMA*FLUXES.
      sum7(:,:) = 0. ! sum(so2b)(sigma*F).   SIGMA*FLUXES.
      sum8(:,:) = 0. ! sum(so)(sigma*F).     SIGMA*FLUXES.
      sum9(:,:) = 0. ! sum(so3)(sigma*F).    SIGMA*FLUXES.
!
!**co2 x-section longward to 1640 A have a correction factor that must be added:
!  it is wavelength and somewhat T-dependent
!
!
      do l=1,lmax
        do i=lon0,lon1
          do k=lev0,lev1
            IF(wla(l) .LT. 1640.) THEN
               ph=0.
            END IF 
            IF(wla(l) .GE. 1640. .AND. wla(l) .LT. 1900.) THEN
               ph=0.5+5.E-3*(wla(l)-1740.)
            END IF 
            IF(wla(l) .GE. 1900.) THEN
                  ph=0.5+5.E-3*(1900.-1740.)
            END IF
          enddo
        enddo
      enddo

      do l=1,l1-1 ! 1->31   2250-1050A
        do i=lon0,lon1
          do k=lev0,lev1
! Nov, 2010 btf: Move cxco2 assignment from above loop to here because 
!                its dimensioned l1-1.
            cxco2(l,k,i) = (1.-ph/100.)**(298.-tni(k,i))
            sigchap(k,i) = sigsrc_co2(l)*cxco2(l,k,i)*scco2(k,i)+
     |                     sigsrc_o2 (l)*sco2(k,i)   +
     |                     sigsrc_so2a(l)*scso2(k,i) +
     |                     sigsrc_so2b(l)*scso2(k,i) +
     |                     sigsrc_so(l)*scso(k,i)    +
     |                     sigsrc_so3(l)*scso3(k,i)
            fexp(k,i) = fsrc(l)*exp(-sigchap(k,i))*sfactor
            sum3(k,i) = sum3(k,i)+ sigsrc_co2(l)*cxco2(l,k,i)*fexp(k,i)
            sum4(k,i) = sum4(k,i)+ hc*rlmsrcinv(l)*fexp(k,i)*
     |                  sigsrc_co2(l)*cxco2(l,k,i)
            sum5(k,i) = sum5(k,i)+ sigsrc_o2 (l)*fexp(k,i)
            sum6(k,i) = sum6(k,i)+ sigsrc_so2a(l)*fexp(k,i)
            sum7(k,i) = sum7(k,i)+ sigsrc_so2b(l)*fexp(k,i)
            sum8(k,i) = sum8(k,i)+ sigsrc_so(l)*fexp(k,i)
!  ** JSO3 path to SO2 + O+ dissociation ----------------------------
            sum9(k,i) = sum9(k,i)+ sigsrc_so3(l)*fexp(k,i)*brso3tso2(l)
          enddo
        enddo
      enddo

!     write(6,"('qrj: lat=',i4,' sum6=',2e12.4,' sum7=',2e12.4)")
!    |  lat,minval(sum6),maxval(sum6),minval(sum7),maxval(sum7)

      if (debug) write(6,"(/,'qrj after wavelength loop2: lat=',i3)")
     | lat
!---------------------------------------------------------------
!  JCO2 and JO2 Dissociation Rate modification : xnmbari(k,i) factor
!  QSRC added to QTOTAL
!---------------------------------------------------------------
! Longitude and column domain of current process:
! Longitude-Column Loop
! -------------------------------------------------------------------
!
! Update qtotal : From mostly CO2 SRB-SRC Absorption Longward of 1050A:
! -- Multiply Q by SRC-SRB efficiency factor: Venus value from Fox (1988)
! -- convert to erg/gm/sec units (internal)
! -- convert to erg/cm3/sec units (comparisons)
! Update rj  : From CO2 SRB-SRC Dissociation Longward of 1050A:
! Update qjo2: From O2 SRB-SRC Dissociation Longward of 1050A:
! Update qjso2a: From SO2a SRB-SRC Dissociation Longward of 1050A:
! Update qjso2b: From SO2b SRB-SRC Dissociation Longward of 1050A:
! Update qjso: From SO SRB-SRC Dissociation Longward of 1050A:
! Update qjso3: From SO3 SRB-SRC Dissociation Longward of 1050A:
!
! --------------------------------------------------------------------
      do i=lon0,lon1
        do k=lev0,lev1
          qsrc(k,i)=sum4(k,i)*co2i(k,i)*rmassinv_co2*srceff(k)
     |       *xnmbari(k,i)
          qtotal(k,i,lat) = qtotal(k,i,lat)+ sum4(k,i)*avo
     |       *co2i(k,i)*rmassinv_co2*srceff(k)
          qtot(k,i)=  qeuv(k,i)+qsrc(k,i)
          rj(k,i,lat) = rj(k,i,lat)+sum3(k,i)*co2i(k,i)*
     |       rmassinv_co2*xnmbari(k,i)
          qjo2(k,i,lat) = qjo2(k,i,lat)+sum5(k,i)*o2i(k,i)*
     |       rmassinv_o2*xnmbari(k,i)
          qjso2a(k,i,lat) = qjso2a(k,i,lat)+sum6(k,i)*so2i(k,i)*
     |       rmassinv_so2*xnmbari(k,i)
          qjso2b(k,i,lat) = qjso2b(k,i,lat)+sum7(k,i)*so2i(k,i)*
     |       rmassinv_so2*xnmbari(k,i)
          qjso(k,i,lat) = qjso(k,i,lat)+sum8(k,i)*soi(k,i)*
     |       rmassinv_so*xnmbari(k,i)
          qjso3(k,i,lat) = qjso3(k,i,lat)+sum9(k,i)*so3i(k,i)*
     |       rmassinv_so3*xnmbari(k,i)
        enddo
      enddo
!
! 1050-1750A Neutral Heating at Venus (erg/cm3/sec): Scaled Incl. Efficiency
!     call addfld('QSRC',' ','erg/cm3/sec',qsrc(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
! TOTAL-Neutral Heating at Venus (erg/cm3/sec): Scaled Incl. Efficiency
!     call addfld('QTOT',' ','erg/cm3/sec',qtot(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
! TOTAL-Neutral Heating at Venus (erg/g/sec): Scaled Incl. Efficiency
!     call addfld('QTOTAL3',' ','erg/g/sec',qtotal(:,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
!---------------------------------------------------------------
! Contributions from MUV (3050 - 2250A) Needed Yet for JSO3
! KINETICS 1-D model XSEC and sfluxmuv used
! - Test combined EUV/SRC/SRB contributions to JSO2a, JSO2b first
! - Test combined SRC/SRB contributions to JSO first
! - Test combined SRC/SRB contributions to JSO3 first
! - Test combined SRC/SRB/MUV contributions to JSO3 last
!   (off for now: no yield for WLA > 1950A for JSO3--> SO2 + O)
!---------------------------------------------------------------
! Initialize sum10
!     sum10(:,:) = 0. ! sum(so3)(sigma*F).    SIGMA*FLUXES.
!
! Wavelength loop --------------------------------------------------
! Summation over wavelength: MUV bins (16): 3050-2250A
! Multiply Q by (1/r**2) mean heliocentric distance of Venus
!     do l=1,lmax2-1 ! 1->16   3050-2250A
!       do i=lon0,lon1
!         do k=lev0,lev1
!           sigchap(k,i) = sigmuv_so3(l)*scso3(k,i)
!          fexp(k,i) = fmuv(l)*exp(-sigchap(k,i))*sfactor
!           sum10(k,i) = sum10(k,i)+ sigmuv_so3(l)*fexp(k,i)
!         enddo
!       enddo
!     enddo
! Update Loop -----------------------------------------------------
!  Update qjso3: From SO3 MUV Dissociation Longward of 2250A
!     do i=lon0,lon1
!       do k=lev0,lev1
!         qjso3(k,i,lat) = qjso3(k,i,lat)+sum10(k,i)*so3i(k,i)*
!    |       rmassinv_so3*xnmbari(k,i)
!       enddo
!     enddo
! -------------------------------------------------------------------
!  Frequency diagnostics and allocated arrays for comp_so, comp_so2
      do i=lon0,lon1
        do k=lev0,lev1
!  Standard
        fjco2(k,i) = rj(k,i,lat)*rmass_co2/(xnmbari(k,i)*co2i(k,i))
        fjo2(k,i) = qjo2(k,i,lat)*rmass_o2/(xnmbari(k,i)*o2i(k,i))
!  SOX
        fjso2a(k,i) = qjso2a(k,i,lat)/so2den(k,i)
        qfjso2a(k,i,lat) = fjso2a(k,i)
        fjso2b(k,i) = qjso2b(k,i,lat)/so2den(k,i)
        qfjso2b(k,i,lat) = fjso2b(k,i)
        fjso(k,i) = qjso(k,i,lat)/soden(k,i)
        qfjso(k,i,lat) = fjso(k,i)
        fjso3(k,i) = qjso3(k,i,lat)/so3den(k,i)
        qfjso3(k,i,lat) = fjso3(k,i)
        enddo
      enddo
! -------------------------------------------------------------------
!   **Fields to 12-allocation slots**
!    rj    (lev0:lev1,lon0:lon1,lat) = 0.
!    qjo2  (lev0:lev1,lon0:lon1,lat) = 0.  
!    qjn2  (lev0:lev1,lon0:lon1,lat) = 0.  
!    qjso  (lev0:lev1,lon0:lon1,lat) = 0.
!    qjso2a (lev0:lev1,lon0:lon1,lat) = 0.
!    qjso2b (lev0:lev1,lon0:lon1,lat) = 0.
!    qjso3 (lev0:lev1,lon0:lon1,lat) = 0.
!    qtotal(lev0:lev1,lon0:lon1,lat) = 0.
!    qop   (lev0:lev1,lon0:lon1,lat) = 0.
!    qcop  (lev0:lev1,lon0:lon1,lat) = 0.
!    qn2p  (lev0:lev1,lon0:lon1,lat) = 0.
!    qco2p (lev0:lev1,lon0:lon1,lat) = 0.
!    qo2p  (lev0:lev1,lon0:lon1,lat) = 0.
!
!    FINAL QRJ SUBROUTINE CHECKS -------------------------------
!    -- units:  Q(erg/gm/sec), Productions (#/cm3/sec)
!     call addfld('Q'     ,' ',' ',qtotal(:,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('JCO2'    ,' ',' ',rj(:,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('JO2'   ,' ',' ',qjo2(:,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('JN2'   ,' ',' ',qjn2(:,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     write(6,"('qrj call addfld for qjso2a,b: lat=',i4,
!    |  ' qjso2a min,max=',2e12.4,' qjso2b min,max=',2e12.4)") lat,
!    |  minval(qjso2a(:,lon0:lon1,lat)),maxval(qjso2a(:,lon0:lon1,lat)),
!    |  minval(qjso2b(:,lon0:lon1,lat)),maxval(qjso2b(:,lon0:lon1,lat))
!  Photo Rates (#/cm3/s)
!     call addfld('JCO2'    ,' ',' ',rj(:,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('JO2'   ,' ',' ',qjo2(:,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('QJSO2A'   ,' ',' ',qjso2a(:,lon0:lon1,lat),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('QJSO2B'   ,' ',' ',qjso2b(:,lon0:lon1,lat),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('QJSO'   ,' ',' ',qjso(:,lon0:lon1,lat),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('QJSO3'   ,' ',' ',qjso3(:,lon0:lon1,lat),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!  Photo Frequencies (#/s)
!     call addfld('SOI',' ',' ',soi,
!    |  'ilev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('JCO2'    ,' ',' ',fjco2,
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('JO2'   ,' ',' ',fjo2,
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('JSO2A'   ,' ',' ',fjso2a,
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('JSO2B'   ,' ',' ',fjso2b,
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('JSO'   ,' ',' ',fjso,
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('JSO3'   ,' ',' ',fjso3,
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)

!     write(6,"('qrj call addfld for zjso2a,b: lat=',i4,
!    |  ' zjso2a min,max=',2e12.4,' zjso2b min,max=',2e12.4)") lat,
!    |  minval(zjso2a(:,lon0:lon1)),maxval(zjso2a(:,lon0:lon1)),
!    |  minval(zjso2b(:,lon0:lon1)),maxval(zjso2b(:,lon0:lon1))

!     call addfld('ZJSO2A'   ,' ',' ',zjso2a(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('ZJSO2B'   ,' ',' ',zjso2b(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('JN2PD'   ,' ',' ',qjn2pd(:,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('JN2PE'   ,' ',' ',qjn2pe(:,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('JN22',' ',' ',qjn22(:,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('JN2PE2',' ',' ',qjn2pe2(:,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('ZJN2PE',' ','/sec',zjn2pepd(:,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('PEPDGRID',' ','/sec',pepdgrid2(:,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('POP'   ,' ',' ',qop(:,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('PCOP'   ,' ',' ',qcop(:,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('PN2P'   ,' ',' ',qn2p(:,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('PCO2P'   ,' ',' ',qco2p(:,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('PO2P'   ,' ',' ',qo2p(:,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
      if (debug) write(6,"(/,'qrj at end of qrj sub : lat=',i3)") lat
!
!---------------------------------------------------------------
#ifdef VT
       code = 118 ; state = 'qrj' ; activity='ModelCode'
       call vtend(118,ier)
#endif
      end subroutine qrj
!---------------------------------------------------------------
      subroutine init_sflux
      use input_module,only: f107,f107a
      use init_module,only: sfeps
!
! Flux initialization once per time step, called from advance.
!
! Local:
      real :: euvflx(neuv),wave1(lmax),wave2(lmax)
      real :: flya,hlybr,fexvir,hlya,heiew,xuvfac
      integer :: iscale,n,nn
!
! Use f107,f107a from input_mod (these are either provided by
! the user, or obtained from GPI database, see input_module.F and
! gpi_module.F).
!
      hlybr = 0.
      fexvir = 0.
      hlya = 3.e+11+0.4E+10*(f107-70.)
      heiew = 0.
      xuvfac = 4.0 - (f107-68.0) / (243.0-68.0)
!     xuvfac = 20.0 - (f107-68.0) / (243.0-68.0)
      if (xuvfac .lt. 1.0) xuvfac = 1.0
!
! ISCALE   =0 for Hinteregger contrast ratio method
!          =1 for Hinteregger linear interpolation
!          =2 for Tobiska EUV91 model
!          =3 for Woods & Rottman 10 Nov. 1988 measurement
!          =4 for Woods & Rottman 20 Jun. 1989 measurement
! F107     daily 10.7 cm flux (1.E-22 W m-2 Hz-1)
! F107A    81-day centered average 10.7 cm flux
! HLYBR    ratio of H Ly-b 1026A flux to solar minimum value (optional)
! FEXVIR   ratio of Fe XVI 335A flux to solar minimum value (optional)
! HLYA     H Lyman-alpha flux (photons cm-2 s-1) (optional)
! HEIEW    He I 10830A equivalent width, (milliAngstroms) (optional)
! XUVFAC   factor for scaling flux 16-250A (optional)
! WAVE1    longwave bound of spectral intervals (Angstroms)
! WAVE2    shortwave bound of intervals (= WAVE1 for indiv. lines)
! SFLUX    scaled solar flux returned by subroutine (photons cm-2 s-1)
!
      iscale = 1
      call ssflux(iscale,f107,f107a,hlybr,fexvir,hlya,
     |            heiew, xuvfac, wave1, wave2, sflux)
      call euvac(f107,f107a,euvflx)
!
! Initializing sflux from 2250 A - 1750 A
! f10.7 = 70
! sfluxext units = photons/cm2/sec
      sfluxext = (/
     |             0.23000E+5, 0.33000E+5, 0.18200E+5, 0.14600E+5,
     |             0.47800E+4, 0.36000E+4, 0.18000E+4, 0.41400E+3,
     |             0.36000E+3, 0.31000E+3, 0.31100E+3, 0.24200E+3,
     |             0.24200E+3, 0.19200E+3, 0.14300E+3, 0.14300E+3/)
! f10.7 = 130
!     sfluxext = (/
!    |             0.25500E+5, 0.36600E+5, 0.20200E+5, 0.16200E+5,
!    |             0.53100E+4, 0.40100E+4, 0.20000E+4, 0.46000E+3,
!    |             0.40100E+3, 0.34500E+3, 0.34500E+3, 0.26800E+3,
!    |             0.26800E+3, 0.21300E+3, 0.15900E+3, 0.15900E+3/)
!
! f10.7 = 200
!     sfluxext = (/
!    |             0.27600E+5, 0.39600E+5, 0.21900E+5, 0.17500E+5,
!    |             0.57400E+4, 0.43400E+4, 0.21600E+4, 0.49700E+3,
!    |             0.43500E+3, 0.37300E+3, 0.37300E+3, 0.29000E+3,
!    |             0.29000E+3, 0.23100E+3, 0.17100E+3, 0.17100E+3/)
!
! lmax here =75, but in the ssflux routine it is stated to be 59
!  Correction factor (1.E9) applied to sfluxext b/c it is already
!   applied to sflux
!
! Transfer values of sflux to appropriate slots
!
      do n = l1,lmax
        feuv(n) = sflux(n-16)*sfeps
      enddo
      do n = 17,l1-1
        fsrc(n)   = sflux(n-16)*sfeps
      enddo
      do n = 1,16
        fsrc(n)   = sfluxext(n)*sfeps*1.E9
      enddo
      do n=1,neuv
        nn = n+31
        feuv(nn) = euvflx(n)*sfeps
      enddo
!     write(6,"('init_sflux: feuv=',/,(6e12.4))") feuv
!     write(6,"('init_sflux: fsrc=',/,(6e12.4))") fsrc
!
!     write(6,"('  n   wave1   wave2       sflux      ',
!    |  'euvflx        fsrc        feuv')")
!     do n=1,lmax                   ! 1->59
!       if (n < l1) then            ! 1->15
!         write(6,"(i3,2f8.2,e12.4,12x,e12.4,e12.4)") 
!    |      n,wave1(n),wave2(n),sflux(n),fsrc(n),feuv(n)
!       elseif (n < l1+neuv-1) then ! 16->52
!         write(6,"(i3,2f8.2,e12.4,e12.4,12x,e12.4)") 
!    |      n,wave1(n),wave2(n),sflux(n),euvflx(n-l1+1),feuv(n)
!       else                        ! 53->59
!         write(6,"(i3,2f8.2,e12.4,12x,12x,e12.4)") 
!    |      n,wave1(n),wave2(n),sflux(n),feuv(n)
!       endif 
!     enddo
      end subroutine init_sflux
!-----------------------------------------------------------------
      subroutine init_qrj
!
! Called once per run, from init_module.
!
! TOTAL ABSORPTION COEFFICIENTS
!   sigeuv(1,N) = SIGA(O1)
!   sigeuv(2,N) = SIGA(CO)
!   sigeuv(3,N) = SIGA(N2)
!   sigeuv(4,N) = SIGA(CO2)
!   sigeuv(5,N) = SIGA(O2)
!   sigeuv(6,N) = SIGA(SO2a)
!   sigeuv(7,N) = SIGA(SO2b)
! TOTAL ABSORPTION COEFFICIENTS
!   sigsrc_co2(31) = SIGA(CO2)
!   sigsrc_o2(31) = SIGA(O2)
!   sigsrc_so(31) = SIGA(SO)
!   sigsrc_so2a(31) = SIGA(SO2a)
!   sigsrc_so2b(31) = SIGA(SO2b)
!   sigsrc_so3(31) = SIGA(SO3)
!
! TOTAL IONIZATION COEFFICIENTS
!   SIGMAS(1,N) = SIGI(O1)
!   SIGMAS(2,N) = SIGI(CO)
!   SIGMAS(3,N) = SIGI(N2-TOTAL)
!   SIGMAS(4,N) = SIGI(CO2-TOTAL)
!   SIGMAS(5,N) = SIGI(O2)
!
! BRANCHING RATIOS 
!   BRN2P(L)  = IN2  to N2+
!   BRCO2P(L) = ICO2 to CO2+
!   BROP(L)   = ICO2 to O+
!
! ** replicate MTGCM - sflux bins are 20 A (1920-1750 A)
!                                      50 A (2250-1950 A)
!
! **rlmeuv units cm
!
      rlmeuv = (/
     |           0.22500E-04, 0.22000E-04, 0.21500E-04, 0.21000E-04,
     |           0.20500E-04, 0.20000E-04, 0.19500E-04, 0.19200E-04,
     |           0.18900E-04, 0.18700E-04, 0.18500E-04, 0.18300E-04,
     |           0.18100E-04, 0.17900E-04, 0.17700E-04, 0.17500E-04,
     |           0.17250E-04, 0.16750E-04, 0.16250E-04, 0.15750E-04,
     |           0.15250E-04, 0.14750E-04, 0.14250E-04, 0.13750E-04,
     |           0.13250E-04, 0.12750E-04, 0.12250E-04, 0.12157E-04,
     |           0.11750E-04, 0.11250E-04, 0.10750E-04, 0.10250E-04,
     |           0.10319E-04, 0.10257E-04, 0.97500E-05, 0.97702E-05,
     |           0.92500E-05, 0.87500E-05, 0.82500E-05, 0.77500E-05,
     |           0.78936E-05, 0.77041E-05, 0.76515E-05, 0.72500E-05,
     |           0.70331E-05, 0.67500E-05, 0.62500E-05, 0.62973E-05,
     |           0.60976E-05, 0.57500E-05, 0.58433E-05, 0.55437E-05,
     |           0.52500E-05, 0.47500E-05, 0.46522E-05, 0.42500E-05,
     |           0.37500E-05, 0.36807E-05, 0.32500E-05, 0.30378E-05,
     |           0.30331E-05, 0.27500E-05, 0.28415E-05, 0.25630E-05,
     |           0.22500E-05, 0.17500E-05, 0.12500E-05, 0.75000E-06,
     |           0.41000E-06, 0.27500E-06, 0.19500E-06, 0.12000E-06,
     |           0.60000E-07, 0.30000E-07, 0.15000E-07/)
!
! -------------------------------------------------------------------
!  Total Absorptions (O, CO, N2, CO2(195, 295), O2)
! -------------------------------------------------------------------
!
!    O Total Absorption (75 bins)S&N(2000)--------------------
!
      sigeuv(1,:) = (/
     |             0.,0.,0.,0.,0.,0.,0.,0.,
     |             0.,0.,0.,0.,0.,0.,0.,0.,
     |             0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |             0.,0.,0.,0.,0.,0.,
     |             0.00, 0.00, 1.32, 4.55, 3.50, 5.09,
     |             3.75, 3.89, 4.00,10.74,11.46,17.25,
     |            13.37,13.40,13.40,13.02,13.09,12.59,
     |            12.06,12.13,11.93,11.50, 9.69, 9.84,
     |             8.69, 7.70, 7.68, 6.46, 7.08, 6.05,
     |             5.20, 3.73, 1.84, 0.73, 0.16, 0.05,
     |             0.51, 0.07, .012, .002, .0002/)
!
!    CO Total Absorption (75-bins)S&N(2000)-----------------
!    Look for xco for xray bins (1-50 A)
!
      sigeuv(2,:) = (/
     |              0.,0.,0.,0.,0.,0.,0.,0.,
     |              0.,0.,0.,0.,0.,0.,0.,0.,
     |              0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |              0.,0.,0.,
     |              8.568,1.388,1.388,52.83,28.50,50.318,
     |              36.976,22.608,20.535,33.132,15.262,26.276,
     |              25.501,23.853,21.036,21.915,22.025,22.1,
     |              21.91,22.00,21.624,21.685,20.574,20.173,
     |              17.956,14.956,15.259,13.44,11.90,11.867,
     |              9.424,10.541,8.614,7.011,4.671,2.391,
     |              0.866,0.,0.,0.,0.,0.,0.,0./)
!
!    N2 Total Absorption (75-bins)S&N(2000)-----------------
!
      sigeuv(3,:) = (/
     |            0.,0.,0.,0.,0.,0.,0.,0.,
     |            0.,0.,0.,0.,0.,0.,0.,0.,
     |            0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |            0.,0.,0.,0.,0.,0.,
     |            50.99, 2.24, 9.68,20.25,16.99,33.58,
     |            16.49,14.18,120.49,24.66,26.54,31.76,
     |            23.34,23.37,22.79,22.79,22.40,24.13,
     |            24.50,23.47,23.16,21.68,16.40,16.91,
     |            13.86,11.70,11.67,10.49,10.90,10.21,
     |             8.39, 4.96, 2.26, 0.72, 0.24, 1.16,
     |             0.48, 0.09, .015, .003, .0003/)
!
!    CO2 Total Absorption (75-bins)S&N(2000) 
!    (Gonzalez-Galindo,06) xco2(195K, 295K) longward of 1200A
!    Check xco2 shortward of 50A
!    (vtgcm2d.v1 model interpolation is needed wrt T)
!    xco2(295k)
      sigeuv(4,:) = (/
     |             1.0E-12,1.0E-11,1.0E-9,1.0E-8,1.0E-7,
     |             1.3E-6,2.0E-5,6.0E-5,1.0E-4,2.0E-4,
     |             4.0E-4,4.5E-4,1.0E-3,2.0E-3,3.0E-3,
     |             5.0E-3,
     |             0.016,0.048,0.126,0.258,0.431,0.544,
     |             0.524,0.607,0.739,0.373,0.081,0.0496,
     |             0.882,40.80,17.40,18.24,14.20,15.10,
     |             50.31,42.87,52.08,44.67,13.98,39.83, 
     |             26.49,61.94,93.84,32.83,23.52,33.70,
     |             34.45,34.30,35.30,34.91,34.20,33.20,
     |             31.49,30.25,29.53,28.27,24.87,25.27,
     |             23.57,21.59,21.49,17.52,19.02,16.51,
     |             14.36, 9.09, 4.62, 1.55,0.358,0.358,
     |             0.358,0.358,0.358,0.358,0.358/)
!
!    xco2(195k)
!     sigeuv(4,:) = (/
!    |             0.010,0.035,0.100,0.222,0.392,0.518,
!    |             0.515,0.578,0.664,0.336,0.071,0.0433,
!    |             0.882,40.80,17.40,18.24,14.20,15.10,
!    |             50.31,42.87,52.08,44.67,13.98,39.83, 
!    |             26.49,61.94,93.84,32.83,23.52,33.70,
!    |             34.45,34.30,35.30,34.91,34.20,33.20,
!    |             31.49,30.25,29.53,28.27,24.87,25.27,
!    |             23.57,21.59,21.49,17.52,19.02,16.51,
!    |             14.36, 9.09, 4.62, 1.55,0.358 0.358,
!    |             0.358,0.358,0.358,0.358,0.358/)
!
!    O2 Total Absorption (75-bins)S&N(2000) & VTGCM2-----
      sigeuv(5,:) = (/
     |            1.0E-5,1.0E-5,1.0E-5,1.0E-5,1.0E-5,
     |            2.0E-5,5.0E-5,3.0E-4,2.0E-3,5.0E-3,
     |            1.0E-2,1.0E-2,2.0E-2,3.5E-2,5.0E-2,
     |            2.0E-1,
     |             0.50, 1.50, 3.40, 6.00,10.00,13.00,
     |            15.00,12.00, 2.20, 0.30, 3.00, 0.01,
     |             0.30, 0.10, 1.00, 1.35, 1.05, 1.63,
     |            21.11,18.73,12.82, 8.56,16.63,22.15,
     |            26.67,18.91,20.80,28.54,27.44,21.92,
     |            26.02,32.06,28.07,26.61,22.72,26.04,
     |            24.60,23.10,21.91,20.31,18.12,18.32,
     |            17.44,16.81,16.80,14.39,15.79,13.37,
     |            10.90, 7.51, 3.81, 1.32, 0.32, 0.10, 
     |             1.02, 0.14, .024, .004, .0004/)
!
!    SO2a Dissociation xsec (75-bins) KINETICS 1-D
!    Zhang et al. (2012): 629.7-2275A bins (non-zero)
      sigeuv(6,:) = (/
     |             0.000, 0.531, 2.391, 4.028, 5.875, 7.339,
     |             7.510, 6.147, 6.147, 4.139, 4.139, 4.139,
     |             2.235, 2.235, 1.145, 1.145, 0.815, 0.429,
     |             0.797, 2.350, 4.199, 4.566, 3.850, 4.260,
     |             11.04, 33.78, 36.20, 24.64, 9.540, 11.05,
     |             8.470, 15.81, 0.000, 16.28, 16.60, 18.07,
     |             10.70, 5.150, 3.950, 5.490, 0.000, 6.707,
     |             0.000, 0.800, 0.000, 0.370, 0.249, 0.000,
     |             0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     |             0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     |             0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     |             0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     |             0.000, 0.000, 0.000/)
!
!    SO2b Dissociation xsec (75-bins) KINETICS 1-D
!    Zhang et al. (2012): 629.7-2125 bins (non-zero)
      sigeuv(7,:) = (/
     |             0.000, 0.000, 0.000, 0.020, 0.059, 0.074,
     |             0.076, 0.062, 0.062, 0.042, 0.042, 0.042,
     |             0.023, 0.023, 0.017, 0.012, 0.009, 0.009,
     |             0.071, 0.120, 0.221, 0.240, 0.203, 0.224,
     |             2.410, 13.38, 15.50, 10.56, 40.90, 51.30,
     |             49.78, 9.490, 0.000, 9.770, 10.11, 10.84,
     |             6.420, 3.090, 2.370, 3.290, 0.000, 4.024,
     |             0.000, 0.480, 0.000, 0.222, 0.150, 0.000,
     |             0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     |             0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     |             0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     |             0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     |             0.000, 0.000, 0.000/)
!
!    SO3 Dissociation xsec (75-bins) KINETICS 1-D
!    Zhang et al. (2012): 1155.0-2275 bins (non-zero)
!    31-bins (SRC/LYA/SRB) non-zero
!    S. W. Bougher (10/1/2013)
!    From Moses et al (XXXX): wavelenghts below 195-nm included
!     sigeuv(8,:) = (/
!    |             10.580, 13.700, 29.570, 30.210, 28.040, 25.870,
!    |             23.700, 21.520, 21.500, 19.350, 19.350, 19.350,
!    |             17.180, 17.180, 16.100, 16.100, 13.780, 11.760,
!    |             10.050, 14.150, 24.970, 28.330, 18.520, 13.470,
!    |             38.420, 141.40, 155.00, 105.60, 40.680, 51.080,
!    |             49.790,  99.34,   0.00,   0.00, 119.80,  0.000,
!    |             0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
!    |             0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
!    |             0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
!    |             0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
!    |             0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
!    |             0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
!    |             0.000, 0.000, 0.000/)
!
!    SO3 Dissociation xsec (75-bins) KINETICS 1-D
!    Zhang et al. (2012): 1155.0-2275 bins (non-zero)
!    31-bins (SRC/LYA/SRB) non-zero
!    S. W. Bougher (10/1/2013)
!    From Zhang et al (2012): wavelenghts above 200-nm included (only)
      sigeuv(8,:) = (/
     |             0.1165, 0.1910, 0.3097, 0.4640, 0.6158, 0.7320,
     |             0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     |             0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     |             0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     |             0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     |             0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     |             0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     |             0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     |             0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     |             0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     |             0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     |             0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     |             0.000, 0.000, 0.000/)
!
!
!    SO3 Dissociation xsec (16-bins) KINETICS 1-D
!    Zhang et al. (2012): 2250-3025A bins (non-zero)
!    16-bins (MUV) non-zero
!    ** Update these with SO3 values less than 3025 A
      sigmuv_so3(:) = (/
     |             0.471, 1.863, 4.168, 7.611, 8.550, 12.24,
     |             11.13, 7.973, 7.973, 6.576, 5.178, 4.052,
     |             4.052, 2.764, 2.764, 2.603/)
!
!    SO Dissociation xsec (31-bins) KINETICS 1-D
!    Zhang et al. (2012): 1155.0-2275 bins (non-zero)
!    31-bins (SRC/LYA/SRB) non-zero
      sigsrc_so(:) = (/
     |             0.471, 1.863, 4.168, 7.611, 8.550, 12.24,
     |             11.13, 7.973, 7.973, 6.576, 5.178, 4.052,
     |             4.052, 2.764, 2.764, 2.603, 2.603, 2.603,
     |             2.603, 2.603, 2.603, 2.603, 2.603, 4.378,
     |             31.14, 93.37, 117.2, 94.38, 35.12, 0.000,
     |             0.000/)
!
!
!
! -------------------------------------------------------------------
!  Total Ionizations (O, CO, N2, CO2, O2)
! -------------------------------------------------------------------
!
!    O+ Ionization same as O absorption (75-bins)
!    S&N(2000) and VTGCM2
!
      sigmas(1,:) = (/
     |             0.,0.,0.,0.,0.,0.,0.,0.,
     |             0.,0.,0.,0.,0.,0.,0.,0.,
     |             0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |             0.,0.,0.,0.,0.,0.,
     |             0.00, 0.00, 1.32, 4.55, 3.50, 5.09,
     |             3.75, 3.89, 4.00,10.74,11.46,17.25,
     |            13.37,13.40,13.40,13.02,13.09,12.59,
     |            12.06,12.13,11.93,11.50, 9.69, 9.84,
     |             8.69, 7.70, 7.68, 6.46, 7.08, 6.05,
     |             5.20, 3.73, 1.84, 0.73, 0.16, 0.05,
     |             0.51, 0.07, .012, .002, .0002/)
!
!    CO Ionization to CO+ (75-bins)S&N (2000) ------------
!    Check for x-section in xray bins (1-50 A)
!
      sigmas(2,:) = (/
     |             0.,0.,0.,0.,0.,0.,0.,0.,
     |             0.,0.,0.,0.,0.,0.,0.,0.,
     |             0.00,  0.00, 0.00, 0.00, 0.00, 0.00,
     |             0.00,  0.00, 0.00, 0.00, 0.00, 0.00,
     |             0.00,  0.00, 0.00, 0.00, 0.00, 0.00,
     |             0.00,  0.00, 0.00,11.70,17.56,11.38,
     |             17.35, 9.88,13.03,21.08,22.12,20.10,
     |             21.85,22.03,21.92,21.90,22.00,21.57,
     |             20.45,19.67,19.09,16.88,13.61,13.98,
     |             11.53, 9.25, 9.21, 6.36, 7.71, 5.45,
     |              4.08, 2.46, 1.07, 0.29, 0.00, 0.00,
     |              0.00, 0.00, 0.00, 0.00, 0.00/)
!
!
!    N2+ Total Ionization (75-bins)S&N(2000)-------------
!
      sigmas(3,:) = (/
     |             0.,0.,0.,0.,0.,0.,0.,0.,
     |             0.,0.,0.,0.,0.,0.,0.,0.,
     |             0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |             0.,0.,0.,0.,0.,0.,
     |             0.00, 0.00, 0.00, 0.00, 0.00,14.27,
     |             8.86, 8.50,65.80,15.06,25.48,29.24,
     |            23.34,23.37,22.79,22.79,22.40,24.13,
     |            24.50,23.47,23.16,21.68,16.40,16.91,
     |            13.86,11.70,11.67,10.49,10.90,10.21,
     |             8.39, 4.96, 2.26, 0.72, 0.24, 1.16,
     |             0.48, 0.09, .015, .003, .0003/)
!
!    CO2+ Total Ionization (75-bins)S&N (2000) -----------
!    Same as CO2 absorption shortward of 200A
!    Check xco2 total ionization shortward of 50A
!--
      sigmas(4,:) = (/
     |             0.,0.,0.,0.,0.,0.,0.,0.,
     |             0.,0.,0.,0.,0.,0.,0.,0.,
     |             0.00,  0.00, 0.00, 0.00, 0.00, 0.00,
     |             0.00,  0.00, 0.00, 0.00, 0.00, 0.00,
     |             0.00,  0.00, 0.00, 0.00, 0.00, 0.00,
     |             0.00,  0.00, 0.00, 7.14,10.93,34.09,
     |             21.68,51.77,86.32,27.49,20.86,32.30,
     |             34.57,34.30,35.30,34.96,33.86,33.21,
     |             32.60,30.58,28.76,27.17,25.52,25.56,
     |             23.63,21.44,21.39,17.14,18.60,16.11,
     |             14.32, 9.09, 4.62, 1.55,0.358,0.358,
     |             0.358,0.358,0.358,0.358,0.358/)
!
!    O2+ Total Ionization (75-bins)S&N (2000)--------------
!    Same as O2 total absorption shortward of 600A
!
!     Total ionization of O2
      sigmas(5,:) = (/
     |             0.,0.,0.,0.,0.,0.,0.,0.,
     |             0.,0.,0.,0.,0.,0.,0.,0.,
     |             0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     |             0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     |             0.00, 0.00, 0.00, 0.26, 0.00, 1.05,
     |            13.94,15.54, 9.37, 5.49, 6.41,10.60,
     |            10.19, 8.47,11.72,23.81,23.75,21.31,
     |            24.94,31.10,26.07,26.61,22.72,26.04,
     |            24.60,23.10,21.91,20.31,18.12,18.32,
     |            17.44,16.81,16.80,14.39,15.79,13.37,
     |            10.90, 7.51, 3.81, 1.32, 0.32, 0.10, 
     |             1.02, 0.14, .024, .004, .0004/)
!
! -------------------------------------------------------------------
!  Branchings for Specific Ionizations (N2+, CO2+, O+)
!  Ratios WRT Total Ionizations Tablulated Above
!  All from S&N(2000)
! -------------------------------------------------------------------
!
!    N2 Ionization Branching ratio yielding N2+ (75-bins)
!    Check IN2--> N2+ Branchings for bins less than 50A.
!
      BRN2 = (/
     |             0.,0.,0.,0.,0.,0.,0.,0.,
     |             0.,0.,0.,0.,0.,0.,0.,0.,
     |             0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |             0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |             1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,
     |             1.,0.990,0.970,0.954,0.955,0.895,
     |             0.908,0.837,0.786,0.787,0.700,0.743,
     |             0.664,0.623,0.636,0.654,0.615,0.600, 
     |             0.600,0.600,0.600,0.600,0.600,0.600/)
!
!    ICO2 Ionization Branching ratio yielding CO2+ (75-bins)
!    Check ICO2--> CO2+ Branchings for bins less than 50A.
!
      BRCO2 = (/
     |             0.,0.,0.,0.,0.,0.,0.,0.,
     |             0.,0.,0.,0.,0.,0.,0.,0.,
     |             0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     |             0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     |             0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     |             0.00, 0.00, 0.00, 1.00, 1.00, 1.00,
     |             1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
     |             0.984,0.985,0.972,0.958,0.956,0.947,
     |             0.925,0.888,0.848,0.776,0.705,0.715,
     |             0.675,0.672,0.673,0.695,0.699,0.690,
     |             0.595,0.546,0.451,0.288,0.200,0.200,
     |             0.200,0.200,0.200,0.200,0.200/)
!
!    ICO2 Ionization Branching ratio yielding O+ (75-bins)
!    Check ICO2--> O+ Branchings for bins less than 50A.
!    
      BRO = (/
     |             0.,0.,0.,0.,0.,0.,0.,0.,
     |             0.,0.,0.,0.,0.,0.,0.,0.,
     |             0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     |             0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     |             0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     |             0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     |             0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     |             0.015,0.014,0.024,0.033,0.035,0.041,
     |             0.050,0.062,0.074,0.096,0.126,0.122,
     |             0.145,0.153,0.153,0.139,0.138,0.140,
     |             0.183,0.212,0.286,0.404,0.400,0.400,
     |             0.400,0.400,0.400,0.400,0.400/)
!
!    ICO2 Ionization Branching ratio yielding CO+ (75-bins)
!    Check ICO2--> CO+ Branchings for bins less than 50A.
!    
      BRCO = (/
     |             0.,0.,0.,0.,0.,0.,0.,0.,
     |             0.,0.,0.,0.,0.,0.,0.,0.,
     |             0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     |             0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     |             0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     |             0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     |             0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     |             0.001,0.001,0.004,0.009,0.009,0.012,
     |             0.025,0.049,0.078,0.124,0.117,0.119,
     |             0.103,0.092,0.092,0.090,0.088,0.091,
     |             0.113,0.116,0.110,0.105,0.100,0.100,
     |             0.100,0.100,0.100,0.100,0.100/)
!   Branching of Photolysis of SO3 to S02+O plus SO2+O(1D)
!     BRSO3TSO2 = (/
!    |             0.048, 0.241, 0.394, 0.537, 0.681, 0.825,
!    |             1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
!    |             1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
!    |             0.980, 0.980, 0.980, 0.980, 0.980, 0.980,
!    |             0.70, 0.70, 0.70, 0.70, 0.70, 0.70,
!    |             0.70, 0.70, 0.00, 0.00, 0.70, 0.00,
!    |             0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
!    |             0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
!    |             0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
!    |             0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
!    |             0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
!    |             0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
!    |             0.00, 0.00, 0.00/)
!     BRSO3TSO2 = (/
!    |             0.048, 0.241, 0.394, 0.537, 0.681, 0.825,
!    |             1.000, 0.000, 0.000, 0.000, 0.000, 0.000,
!    |             0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
!    |             0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
!    |             0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
!    |             0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
!    |             0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
!    |             0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
!    |             0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
!    |             0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
!    |             0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
!    |             0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
!    |             0.00, 0.00, 0.00/)
      BRSO3TSO2 = (/
     |             1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
     |             0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     |             0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     |             0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     |             0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     |             0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     |             0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     |             0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     |             0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     |             0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     |             0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     |             0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     |             0.00, 0.00, 0.00/)
!
!    Venus heating efficiency : 18-22%%
!    From Fox (1988) adn Keating and Bougher (1992)
!    Used in vtgcm2 code (before 1997)
!     Fox 1-D values -----------
!     euveff(:) = 0.20
!     srceff(:) = 0.22
!     Bougher-Brecht 3-D values -----------
!     euveff(:) = 0.22
      euveff(:) = 0.20
      srceff(:) = 0.22
      end subroutine init_qrj
!---------------------------------------------------------------
      subroutine init_sigmas
!
! Called once per run from init_module.
!
      integer :: m,n,nn
      do m=1,nsig
        do n=1,lmax
          sigeuv(m,n)=sigeuv(m,n)*1.E-18
          sigmas(m,n)=sigmas(m,n)*1.E-18
        enddo
      enddo
!  ------------------------------------------------------------
!  CO2 and O2 are main absorbers longward of 1050A for Venus.
!  btf: co2 and o2 for venus:
!  ------------------------------------------------------------
      do n = 1,l1-1
        rlmsrc(n) = rlmeuv(n)
! 9/25/06 btf: co2 and o2 cross-sections:
        sigsrc_co2(n) = sigeuv(4,n)
        sigsrc_o2 (n) = sigeuv(5,n)
! 6/07/13 swb: so2a and so2b cross-sections:
        sigsrc_so2a(n) = sigeuv(6,n)
        sigsrc_so2b(n) = sigeuv(7,n)
! 10/3/13 swb: so3 cross-sections:
        sigsrc_so3(n) = sigeuv(8,n)
      enddo
!  ------------------------------------------------------------
!  Add SO and SO3 as absorbers longward of 1050A for Venus.
!  swb: so and so3 for venus:
!  ------------------------------------------------------------
!  Extension #1; Needed for JSO --> S + O
        do n=1,l1-1
          sigsrc_so(n)=sigsrc_so(n)*1.E-18
        enddo
!  Extension #2 : Not needed for JSO3 --> SO2 + O
!       do n=1,l2-1
!         sigmuv_so3(n)=sigmuv_so3(n)*1.E-18
!       enddo
! -------------------------------------------------------------------
!     Updates turned off, since S&N(2000) values used in init_qrj 
! -------------------------------------------------------------------
!     Overlay EUVAC Cross-sections on sigeuv(:,:) and sigmas(:,:)
!     N2-Stuff on init-euvac only!
!     do n = 1,neuv
!       nn = n+15    ! 16:52
!     Total N2 absorption
!       sigeuv(3,nn) = sigan2(n)
!     Total N2 ionization
!       sigmas(3,nn) = sigin2(n)
!     Branching for N2 ionization yielding N+ (not used)
!       brn2 (nn) = brn2np(n)
!
! 9/25/06 btf: need data statements for brco2p and brop (like brn2np)
!              in init_euvac. When these are available, open up the
!              next 2 lines:
!       brco2(nn) = brco2p(n)
!       bro  (nn) = brop  (n)
!     enddo
! -------------------------------------------------------------------
      end subroutine init_sigmas
!-----------------------------------------------------------------------
      subroutine init_euvac
!
! Called once per run from init_module.
!
! lambdas:
      wleuv1 = 
     +    (/1000.00, 1031.91, 1025.72,  950.00,  977.02,  900.00,
     +       850.00,  800.00,  750.00,  789.36,  770.41,  765.15,
     +       700.00,  703.36,  650.00,  600.00,  629.73,  609.76,
     +       550.00,  584.33,  554.31,  500.00,  450.00,  465.22,
     +       400.00,  350.00,  368.07,  300.00,  303.78,  303.31,
     +       250.00,  284.15,  256.30,  200.00,  150.00,  100.00,
     +        50.00/)
      wleuv2 =
     +    (/1050.00, 1031.91, 1025.72, 1000.00,  977.02,  950.00,
     +       900.00,  850.00,  800.00,  789.36,  770.41,  765.15,
     +       750.00,  703.36,  700.00,  650.00,  629.73,  609.76,
     +       600.00,  584.33,  554.31,  550.00,  500.00,  465.22,
     +       450.00,  400.00,  368.07,  350.00,  303.78,  303.31,
     +       300.00,  284.15,  256.30,  250.00,  200.00,  150.00,
     +       100.00/)
c
c -------------------------------------------------------------------
c These euvac updates not used for the vtgcm. S&N(2000) values instead.
c -------------------------------------------------------------------
c sigan2 : Updated Total N2-absorption cross sections (50-1050A)
      sigan2 =
     +    (/0.00e+00, 0.00e+00, 0.00e+00, 5.10e-17, 2.24e-18, 9.68e-18,
     +      2.02e-17, 1.70e-17, 3.36e-17, 1.65e-17, 1.42e-17, 1.20e-16,
     +      2.47e-17, 2.65e-17, 3.18e-17, 2.33e-17, 2.34e-17, 2.28e-17,
     +      2.28e-17, 2.24e-17, 2.41e-17, 2.45e-17, 2.35e-17, 2.32e-17,
     +      2.17e-17, 1.64e-17, 1.69e-17, 1.39e-17, 1.17e-17, 1.17e-17,
     +      1.05e-17, 1.09e-17, 1.02e-17, 8.39e-18, 4.96e-18, 2.26e-18,
     +      7.20e-19/)
c
c sigin2 : Updated Total N2-ionization cross sections (50-1050A)
      sigin2 =
     +    (/0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     +      0.00e+00, 0.00e+00, 1.43e-17, 8.86e-18, 8.50e-18, 6.58e-17,
     +      1.51e-17, 2.55e-17, 2.92e-17, 2.33e-17, 2.34e-17, 2.28e-17,
     +      2.28e-17, 2.24e-17, 2.41e-17, 2.45e-17, 2.35e-17, 2.32e-17,
     +      2.17e-17, 1.64e-17, 1.69e-17, 1.39e-17, 1.17e-17, 1.17e-17,
     +      1.05e-17, 1.09e-17, 1.02e-17, 8.39e-18, 4.96e-18, 2.26e-18,
     +      7.20e-19/)
c
c N2 -> N+ : N2-ionization branching ratios yielding N+ (50-1050A)
      brn2np =
     +    (/0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     +      0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     +      0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     +      0.00e+00, 0.00e+00, 0.00e+00, 1.30e-03, 3.00e-02, 4.60e-02,
     +      4.50e-02, 1.05e-01, 9.00e-02, 1.63e-01, 2.13e-01, 2.13e-01,
     +      3.00e-01, 2.57e-01, 3.35e-01, 3.77e-01, 3.64e-01, 3.46e-01,
     +      3.85e-01/)
!
      end subroutine init_euvac
!-----------------------------------------------------------------------
      subroutine euvac(F107,F107A,EUVFLX)
!
! This EUV flux model uses the F74113 solar reference spectrum and
! ratios determined from Hinteregger's SERF1 model. It uses the daily
! F10.7 flux (F107) and the 81 day mean (F107A) as a proxy for scaling
! The fluxes are returned in EUVFLX and correspond to the 37 wavelength
! bins of Torr et al. [1979] Geophys. Res. Lett. p771.
! See Richards et al. [1994] J. Geophys. Res. p8981 for details.
!
! Called from init_sflux, once per timestep.
!
! F107   = input daily 10.7 cm flux index. (e.g. 74)
! F107A  = input 81 day average of daily F10.7 centered on current day
! EUVFLX = output array for EUV flux in units of photons/cm2/sec.
!
! Reduced argument list for Titan:
!     real :: wleuv1,wleuv2,sigan2,sigin2,brn2np
!
! Args:
      real,intent(in) :: f107,f107a
      real,intent(out) :: euvflx(neuv)
!
! Local:
      integer :: i
      real :: AFAC(neuv),F74113(neuv),flxfac
      real :: afac_tmp(neuv),f74113_tmp(neuv)
!
! F74113 reference spectrum (doubled below 150-250 A, tripled <150)
! Will be multiplied by 1.0E9 later
!
      F74113=
     | (/1.200,0.450,4.800,3.100,0.460,0.210,1.679,0.800,6.900,
     |   0.965,0.650,0.314,0.383,0.290,0.285,0.452,0.720,1.270,
     |   0.357,0.530,1.590,0.342,0.230,0.360,0.141,0.170,0.260,
     |   0.702,0.758,1.625,3.537,3.000,4.400,1.475,3.500,2.100,
     |   2.467/)
!
! Scaling factors(Ai) for the EUV flux
      AFAC=
     | (/1.0017E-02,7.1250E-03,1.3375E-02,1.9450E-02,2.7750E-03,
     |   1.3768E-01,2.6467E-02,2.5000E-02,3.3333E-03,2.2450E-02,
     |   6.5917E-03,3.6542E-02,7.4083E-03,7.4917E-03,2.0225E-02,
     |   8.7583E-03,3.2667E-03,5.1583E-03,3.6583E-03,1.6175E-02,
     |   3.3250E-03,1.1800E-02,4.2667E-03,3.0417E-03,4.7500E-03,
     |   3.8500E-03,1.2808E-02,3.2750E-03,4.7667E-03,4.8167E-03,
     |   5.6750E-03,4.9833E-03,3.9417E-03,4.4167E-03,5.1833E-03,
     |   5.2833E-03,4.3750E-03/)
!
! 6/20/03 bf: Reverse order of F74113 and AFAC to correspond to descending
! wavelength bins (as in wleuv1,2 and WAVEL,WAVES):
      do i=1,neuv
        afac_tmp(i) = afac(neuv-i+1)  
        f74113_tmp(i) = f74113(neuv-i+1)  
      enddo
      afac = afac_tmp      ! whole array
      f74113 = f74113_tmp
!
! Loop through the wavelengths calculating the scaling factors and
!   the resulting solar flux.
! The scaling factors are restricted to be greater than 0.8
!
      DO 50 I=1,neuv
        FLXFAC=(1.0 + AFAC(I) * (0.5*(F107+F107A) - 80.0))
        IF(FLXFAC.LT.0.8) FLXFAC=0.8
        EUVFLX(I)=F74113(I) * FLXFAC * 1.0E9
 50   CONTINUE
!     write(6,"('euvac: euvflx=',/,(6e12.4))") euvflx
      end subroutine euvac
!-----------------------------------------------------------------------
      subroutine ssflux (iscale, f107, f107a, hlybr, fexvir, hlya,
     |                   heiew, xuvfac, wave1, wave2, sflux)
      integer,parameter :: lmaxflux=59
!
! Args:
      integer,intent(in) :: iscale
      real,intent(in) :: f107,f107a,hlybr,fexvir,hlya,heiew,xuvfac
      real,intent(out) :: wave1(lmaxflux),wave2(lmaxflux),
     |                    sflux(lmaxflux) 
!
! Local:
      real ::
     |          wavel(lmaxflux), waves(lmaxflux), rflux(lmaxflux),
     |          xflux(lmaxflux),
     |          scale1(lmaxflux), scale2(lmaxflux),
     |          tchr0(lmaxflux), tchr1(lmaxflux), tchr2(lmaxflux),
     |          tcor0(lmaxflux), tcor1(lmaxflux), tcor2(lmaxflux),
     |          war1(lmaxflux), war2(lmaxflux),
     |          b1(3), b2(3)
      real :: frat,r1,r2,hlymod,heimod,xuvf
      integer :: l
!
! regression coefficients which reduce to solar min. spectrum:
      b1 = (/1.0, 0.0138, 0.005/) 
      b2 = (/1.0, 0.59425, 0.3811/)
!
! 'best fit' regression coefficients, commented out, for reference:
!     b1 = (/1.31, 0.01106, 0.00492/, B2/-6.618, 0.66159, 0.38319/)
!
      WAVEL =   (/1750.00, 1700.00, 1650.00, 1600.00, 1550.00, 1500.00,
     |            1450.00, 1400.00, 1350.00, 1300.00, 1250.00, 1215.67,
     |            1200.00, 1150.00, 1100.00, 1050.00, 1031.91, 1025.72,
     |            1000.00,  977.02,  950.00,  900.00,  850.00,  800.00,
     |             789.36,  770.41,  765.15,  750.00,  703.31,  700.00,
     |             650.00,  629.73,  609.76,  600.00,  584.33,  554.37,
     |             550.00,  500.00,  465.22,  450.00,  400.00,  368.07,
     |             350.00,  303.78,  303.31,  300.00,  284.15,  256.30,
     |             250.00,  200.00,  150.00,  100.00,   50.00,   32.00,
     |              23.00,   16.00,    8.00,    4.00,    2.00/)
      WAVES =   (/1700.00, 1650.00, 1600.00, 1550.00, 1500.00, 1450.00,
     |            1400.00, 1350.00, 1300.00, 1250.00, 1200.00, 1215.67,
     |            1150.00, 1100.00, 1050.00, 1000.00, 1031.91, 1025.72,
     |             950.00,  977.02,  900.00,  850.00,  800.00,  750.00,
     |             789.36,  770.41,  765.15,  700.00,  703.31,  650.00,
     |             600.00,  629.73,  609.76,  550.00,  584.33,  554.37,
     |             500.00,  450.00,  465.22,  400.00,  350.00,  368.07,
     |             300.00,  303.78,  303.31,  250.00,  284.15,  256.30,
     |             200.00,  150.00,  100.00,   50.00,   32.00,   23.00,
     |              16.00,    8.00,    4.00,    2.00,    1.00/)
      RFLUX =    (/322.00,  168.00,   95.00,   62.00,   44.00,   25.00,
     |              16.90,   11.80,   19.50,    4.10,   11.10,  249.00,
     |               2.78,    0.70,    3.07,    3.64,    3.18,    4.38,
     |               1.78,    5.96,    4.22,    4.43,    1.93,    0.87,
     |               0.79,    0.24,    0.20,    0.17,    0.39,    0.22,
     |               0.17,    1.50,    0.45,    0.48,    1.58,    0.80,
     |               0.51,    0.31,    0.18,    0.39,    0.21,    0.74,
     |               0.87,    6.00,    0.24,    0.84,    0.10,    0.27,
     |               0.92,    1.84,    0.13,    0.38,  0.0215,  0.0067,
     |              1.E-3,   2.E-3,   1.E-5,   5.E-8,   1.E-10/)
      XFLUX =    (/354.00,  191.00,  110.00,   76.00,   55.00,   28.00,
     |              19.60,   14.30,   25.30,    5.00,   17.20,  401.00,
     |               6.26,    1.51,    6.11,    8.66,    9.04,   13.12,
     |               4.42,   13.18,   12.03,   13.29,    5.01,    2.18,
     |               1.59,    0.67,    0.43,    0.43,    0.72,    0.46,
     |               0.48,    3.02,    1.46,    1.02,    4.86,    1.59,
     |               1.57,    1.67,    0.36,    0.99,    2.20,    1.39,
     |               5.63,   11.28,    2.50,    4.14,    3.16,    0.59,
     |               3.70,    4.85,    0.34,    1.15,    0.18,    0.08,
     |             2.5E-2,   5.E-2,   8.E-4,   3.E-5,   5.E-7/)
      SCALE1 = (/     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
     |                0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
     |            1692.09,  405.95, 1516.20, 2731.70, 3314.57, 4375.00,
     |            1316.91, 3621.91, 3908.56, 4432.54, 1541.21,  531.73,
     |             364.83,    0.00,  116.00,  129.41,  162.48,   94.07,
     |              41.29,  709.50,    0.00,  268.47, 1561.05,  367.64,
     |             290.06,  184.36,    0.00,   86.15,    7.50,    0.00,
     |               0.00, 2220.00,    0.00,   61.00,    0.00,   86.95,
     |             206.00,  135.89,   60.35,  157.12,    7.06,    0.75,
     |               0.00,    0.00,    0.00,    0.00,    0.00/)
      SCALE2 =  (/   0.00,    0.00,    0.00,    0.00,    0.00,    0.00,
     |               0.00,    0.00,    0.00,    0.00,    0.00,    0.00,
     |               0.00,    0.00,    0.00,    0.00,    0.00,    0.00,
     |               0.00,    0.00,    0.00,    0.00,    0.00,    0.00,
     |               0.00,    5.34,    0.00,    0.00,    0.00,    0.54,
     |               3.30,    0.00,   12.60,    0.00,    0.00,    0.00,
     |               5.34,   11.63,    2.28,    5.56,   24.93,    8.16,
     |              60.69,    0.00,   28.20,   45.90,   40.80,    1.27,
     |              35.47,   42.80,    1.12,    6.19,    1.26,    0.69,
     |               0.23,    0.46,  7.6E-3,  2.9E-4,  4.8E-6/)
      TCHR0 = (/
     |       0.0,       0.0,       0.0,       0.0,       0.0,       0.0,
     |       0.0,       0.0,       0.0,       0.0,       0.0,       0.0,
     |       0.0,       0.0,       0.0,-4.290E+00,-5.709E+00,-8.493E+00,
     |-1.161E+00,-3.429E+00,-5.464E+00,-6.502E+00,-1.912E+00,-4.034E-01,
     |-1.448E-01, 0.000E+00,-9.702E-02,-6.591E-02,-2.338E-02,-1.273E-01,
     |-2.406E-01,-3.351E-01, 0.000E+00,-1.465E+00,-2.405E+00,-7.975E-02,
     |-4.197E-01,-1.971E-01, 0.000E+00,-5.895E-02,-5.815E-03, 0.000E+00,
     | 0.000E+00, 2.138E-01, 0.000E+00,-7.713E-02, 0.000E+00,-3.035E-02,
     |-2.039E-01,-1.749E-01,-1.041E-01,-2.638E-01,-1.094E-02, 0.000E+00,
     |       0.0,       0.0,       0.0,       0.0,       0.0/)
      TCHR1 = (/
     |       0.0,       0.0,       0.0,       0.0,       0.0,       0.0,
     |       0.0,       0.0,       0.0,       0.0,       0.0,       0.0,
     |       0.0,       0.0,       0.0,-3.023E-13,-3.745E-13,-5.385E-13,
     |-1.211E-13,-3.868E-13,-3.646E-13,-4.125E-13,-1.527E-13,-4.753E-14,
     |-3.411E-14, 0.000E+00,-1.190E-14,-1.034E-14,-1.343E-14,-1.539E-14,
     |-5.174E-14,-6.934E-14, 0.000E+00,-1.215E-13,-1.537E-13,-2.024E-14,
     |-4.596E-14,-1.562E-14, 0.000E+00,-1.221E-14,-1.123E-15, 0.000E+00,
     | 0.000E+00,-2.263E-13, 0.000E+00,-1.508E-14, 0.000E+00,-1.744E-14,
     |-2.100E-14,-1.805E-14,-8.224E-15,-1.919E-14,-7.944E-16, 0.000E+00,
     |       0.0,       0.0,       0.0,       0.0,       0.0/)
      TCHR2 = (/
     |       0.0,       0.0,       0.0,       0.0,       0.0,       0.0,
     |       0.0,       0.0,       0.0,       0.0,       0.0,       0.0,
     |       0.0,       0.0,       0.0, 3.275E-11, 4.057E-11, 6.160E-11,
     | 1.312E-11, 4.189E-11, 4.167E-11, 4.716E-11, 1.654E-11, 5.150E-12,
     | 3.901E-12, 0.000E+00, 1.289E-12, 1.120E-12, 1.455E-12, 1.667E-12,
     | 5.604E-12, 7.931E-12, 0.000E+00, 1.317E-11, 1.757E-11, 2.194E-12,
     | 4.978E-12, 1.693E-12, 0.000E+00, 1.324E-12, 1.285E-13, 0.000E+00,
     | 0.000E+00, 2.586E-11, 0.000E+00, 1.724E-12, 0.000E+00, 1.889E-12,
     | 2.400E-12, 2.063E-12, 8.911E-13, 2.193E-12, 9.090E-14, 0.000E+00,
     |       0.0,       0.0,       0.0,       0.0,       0.0/)
      TCOR0 = (/
     |       0.0,       0.0,       0.0,       0.0,       0.0,       0.0,
     |       0.0,       0.0,       0.0,       0.0,       0.0,       0.0,
     |       0.0,       0.0,       0.0, 0.000E+00, 0.000E+00, 0.000E+00,
     | 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,-6.060E-02,
     | 0.000E+00,-3.399E-02, 0.000E+00, 0.000E+00, 0.000E+00, 4.866E-02,
     |-1.762E-01, 0.000E+00,-2.412E-01, 0.000E+00, 0.000E+00, 0.000E+00,
     |-4.743E-01,-9.713E-01, 5.891E-02,-1.263E-01,-1.246E+00, 2.870E-01,
     |-4.659E+00, 0.000E+00,-1.058E+00,-3.821E+00,-1.874E+00, 0.000E+00,
     |-1.896E+00,-8.505E-01,-2.101E-04,-2.012E-01,-6.097E-02,-2.925E-02,
     |-4.875E-03,       0.0,       0.0,       0.0,       0.0/)
      TCOR1 = (/
     |       0.0,       0.0,       0.0,       0.0,       0.0,       0.0,
     |       0.0,       0.0,       0.0,       0.0,       0.0,       0.0,
     |       0.0,       0.0,       0.0, 0.000E+00, 0.000E+00, 0.000E+00,
     | 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 2.877E-03,
     | 0.000E+00, 1.760E-03, 0.000E+00, 0.000E+00, 0.000E+00, 3.313E-04,
     | 3.643E-03, 0.000E+00, 5.225E-03, 0.000E+00, 0.000E+00, 0.000E+00,
     | 4.085E-03, 1.088E-02, 8.447E-04, 3.237E-03, 1.907E-02, 2.796E-03,
     | 4.460E-02, 0.000E+00, 1.007E-02, 3.481E-02, 1.604E-02, 0.000E+00,
     | 2.029E-02, 2.160E-02, 6.342E-04, 3.594E-03, 5.503E-04, 2.687E-04,
     | 4.479E-05,       0.0,       0.0,       0.0,       0.0/)
      TCOR2 = (/
     |       0.0,       0.0,       0.0,       0.0,       0.0,       0.0,
     |       0.0,       0.0,       0.0,       0.0,       0.0,       0.0,
     |       0.0,       0.0,       0.0, 0.000E+00, 0.000E+00, 0.000E+00,
     | 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.846E-03,
     | 0.000E+00, 1.127E-03, 0.000E+00, 0.000E+00, 0.000E+00, 1.891E-04,
     | 2.326E-03, 0.000E+00, 2.801E-03, 0.000E+00, 0.000E+00, 0.000E+00,
     | 2.446E-03, 7.121E-03, 5.204E-04, 1.983E-03, 1.204E-02, 1.721E-03,
     | 2.911E-02, 0.000E+00, 7.177E-03, 2.272E-02, 9.436E-03, 0.000E+00,
     | 1.316E-02, 1.398E-02, 4.098E-04, 2.328E-03, 3.574E-04, 1.745E-04,
     | 2.909E-05,       0.0,       0.0,       0.0,       0.0/)
      WAR1 =  (/     0.00,    0.00,    0.00,    0.00,    0.00,    0.00,
     |               0.00,    0.00,    0.00,    0.00,    0.00,    0.00,
     |               0.00,    0.00,    3.80,    6.25,    4.93,    6.06,
     |               2.70,    7.07,    8.62,    9.60,    4.54,    2.37,
     |               0.82,    0.33,    0.24,    0.67,    0.28,    0.55,
     |               1.56,    1.11,    0.77,    1.32,    1.71,    0.44,
     |               1.11,    0.95,    0.39,    0.81,    2.00,    1.49,
     |               6.81,    5.07,    1.63,    5.62,    2.08,    0.59,
     |               3.89,    5.19,    0.35,    1.18,    0.099,   0.04,
     |               0.007,   0.00,    0.00,    0.00,    0.00/)
      WAR2 =  (/     0.00,    0.00,    0.00,    0.00,    0.00,    0.00,
     |               0.00,    0.00,    0.00,    0.00,    0.00,    0.00,
     |               0.00,    0.00,   20.80,   17.90,    9.30,   14.30,
     |               6.90,   12.00,   15.60,   18.60,   10.10,    4.30,
     |              12.40,    8.00,    3.60,    1.80,    0.50,    1.40,
     |               3.90,    2.60,    1.60,    3.40,    4.10,    0.70,
     |               4.30,    4.30,    3.80,    2.60,    6.08,    1.35,
     |              12.60,    9.78,    2.96,   10.20,    4.11,    6.68,
     |               6.62,    8.07,    0.47,    1.73,    0.17,    0.075,
     |               0.012,   0.00,    0.00,    0.00,    0.00/)
C
C Linear Interpolation between SC#21REFW and F79050:
C
      FRAT = (F107-68.) / (243.-68.)
      DO 200 L=1,LMAXFLUX
        SFLUX(L) = RFLUX(L) + (XFLUX(L)-RFLUX(L)) * FRAT
  200 CONTINUE
C
C Hinteregger contrast ratio method:
C
      IF (ISCALE .EQ. 0) THEN
        IF (HLYBR .GT. 0.001) THEN
          R1 = HLYBR
        ELSE
          R1 =  B1(1) + B1(2)*(F107A-71.5) + B1(3)*(F107-F107A+3.9)
        ENDIF
        IF (FEXVIR .GT. 0.001) THEN
          R2 = FEXVIR
        ELSE
          R2 =  B2(1) + B2(2)*(F107A-71.5) + B2(3)*(F107-F107A+3.9)
        ENDIF
        DO 100 L=13,LMAXFLUX
          SFLUX(L) = (RFLUX(L) + ((R1-1.)*SCALE1(L)
     |                         +  (R2-1.)*SCALE2(L)) / 1000.)
  100   CONTINUE
      ENDIF
C
C Tobiska EUV91 Method:
C
      IF (ISCALE .EQ. 2) THEN
        IF (HLYA .GT. 0.001) THEN
          HLYMOD = HLYA
        ELSE
          IF (HEIEW .GT. 0.001) THEN
             HLYMOD = HEIEW * 3.77847E9 + 8.40317E10
          ELSE
             HLYMOD = 8.70E8 * F107 + 1.90E11
             HLYMOD = 8.70E8 * F107 + 1.90E11
          ENDIF
        ENDIF
        IF (HEIEW .GT. 0.001) THEN
          HEIMOD = HEIEW * 3.77847E9 + 8.40317E10
        ELSE
          HEIMOD = HLYMOD
        ENDIF
        DO 500 L=16,55
          SFLUX(L) = TCHR0(L) + TCHR1(L)*HLYMOD + TCHR2(L)*HEIMOD
     |             + TCOR0(L) + TCOR1(L)*F107 + TCOR2(L)*F107A
 500    CONTINUE
      ENDIF
C
C Woods and Rottman (10 Nov. 1988) spectrum:
C
      IF (ISCALE .EQ. 3) THEN
        DO 550 L=15,55
        SFLUX(L) = WAR1(L)
  550   CONTINUE
      ENDIF
C
C Woods and Rottman (20 June 1989) spectrum:
C
      IF (ISCALE .EQ. 4) THEN
        DO 560 L=15,55
        SFLUX(L) = WAR2(L)
  560   CONTINUE
      ENDIF
C
C Substitute in H Lyman-alpha and XUVFAC if provided:
C -- Update of Lyman-alpha fluxes over solar cycle
C
      IF (HLYA .GT. 0.001) SFLUX(12) = HLYA / 1.E9
C     IF (XUVFAC .GT. 0.001) THEN
      IF (XUVFAC .GT. 1.0) THEN
        XUVF = XUVFAC
      ELSE
        XUVF = 1.0
      ENDIF
C
C Convert from gigaphotons to photons, etc.:
C
      DO 600 L=1,LMAXFLUX
        WAVE1(L) = WAVEL(L)
        WAVE2(L) = WAVES(L)
        IF (SFLUX(L) .LT. 0.0) SFLUX(L) = 0.0
!       IF (WAVEL(L).LT.251.0 .AND. WAVES(L).GT.15.0)
!    |     SFLUX(L)=SFLUX(L)*XUVF
        IF (WAVEL(L).LT.251.0 .AND. WAVES(L).GT.1.0)
     |     SFLUX(L)=SFLUX(L)*XUVF
        SFLUX(L) = SFLUX(L) * 1.E9
  600 CONTINUE
C
      end subroutine ssflux
!-----------------------------------------------------------------------
      subroutine alloc_q(lon0,lon1,lat0,lat1)
!
! Args:
      integer,intent(in) :: lon0,lon1,lat0,lat1
!
! Local:
      integer :: istat
!
!  -------------------------------------------------------------- 
!  Allocate (19-fields): rj,qjo2,qjn2,qjn2pe,qjn2pd,qjso,qjso2a,qjso2b,qjso3,
!                        qfjso,qfjso2a,qfjso2b,qfjso3,
!                        qtotal,qop,qcop,qn2p,qco2p,qo2p 
!    |  rj,     ! solar photodissociation of JCO2  (#/cm3/sec)
!    |  qjo2,   ! solar photodissociation of JO2   (#/cm3/sec)
!    |  qjn2,   ! solar EUV + PE dissociation of JN2   (#/cm3/sec)
!    |  qjn2pe, ! PE dissociation of JN2 alone         (#/cm3/sec)
!    |  qjn2pd, ! PD dissociation of JN2 alone         (#/cm3/sec)
!    |  qjso,   ! solar photodissociation of SO        (#/cm3/sec)
!    |  qjso2a, ! solar photodissociation of SO2--> SO+O (#/cm3/sec)
!    |  qjso2b, ! solar photodissociation of SO2--> S+O2 (#/cm3/sec)
!    |  qjso3,  ! solar photodissociation of SO3       (#/cm3/sec)
!    |  qfjso,   ! solar photodissociation of SO        (#/sec)
!    |  qfjso2a, ! solar photodissociation of SO2--> SO+O (#/sec)
!    |  qfjso2b, ! solar photodissociation of SO2--> S+O2 (#/sec)
!    |  qfjso3,  ! solar photodissociation of SO3       (#/sec)
!    |  qtotal, ! total euv-uv heating (erg/gm/sec)
!    |  qop,    ! O+ (#/cm3/sec)
!    |  qcop,   ! CO+ (#/cm3/sec)
!    |  qn2p,   ! n2+ (#/cm3/sec)
!    |  qco2p,  ! CO2+ (#/cm3/sec)
!    |  qo2p    ! O2+ (#/cm3/sec)
!  -------------------------------------------------------------- 
      allocate(rj(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_q: error allocating',
     |  ' rj: stat=',i3)") istat
      allocate(qjo2(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_q: error allocating',
     |  ' qjo2: stat=',i3)") istat
      allocate(qjn2(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_q: error allocating',
     |  ' qjn2: stat=',i3)") istat
      allocate(qjn2pe(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_q: error allocating',
     |  ' qjn2pe: stat=',i3)") istat
      allocate(qjn2pd(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_q: error allocating',
     |  ' qjn2pd: stat=',i3)") istat
      allocate(qtotal(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_q: error allocating',
     |  ' qtotal: stat=',i3)") istat
      allocate(qop(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_q: error allocating',
     |  ' qop: stat=',i3)") istat
      allocate(qcop(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_q: error allocating',
     |  ' qcop: stat=',i3)") istat
      allocate(qn2p(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_q: error allocating',
     |  ' qn2p: stat=',i3)") istat
      allocate(qco2p(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_q: error allocating',
     |  ' qco2p: stat=',i3)") istat
      allocate(qo2p(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_q: error allocating',
     |  ' qo2p: stat=',i3)") istat
      allocate(qjso(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_q: error allocating',
     |  ' qjso: stat=',i3)") istat
      allocate(qjso2a(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_q: error allocating',
     |  ' qjso2a: stat=',i3)") istat
      allocate(qjso2b(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_q: error allocating',
     |  ' qjso2b: stat=',i3)") istat
      allocate(qjso3(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_q: error allocating',
     |  ' qjso3: stat=',i3)") istat
      allocate(qfjso(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_q: error allocating',
     |  ' qfjso: stat=',i3)") istat
      allocate(qfjso2a(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_q: error allocating',
     |  ' qfjso2a: stat=',i3)") istat
      allocate(qfjso2b(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_q: error allocating',
     |  ' qfjso2b: stat=',i3)") istat
      allocate(qfjso3(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_q: error allocating',
     |  ' qfjso3: stat=',i3)") istat
      end subroutine alloc_q
!-----------------------------------------------------------------------
      subroutine mk_xnmbari(tni,barm,xnmbari,
     |  lev0,lev1,lon0,lon1,lat)
      use cons_module,only: expz,expzmid,expzmid_inv,boltz,p0
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,dimension(lev0:lev1,lon0:lon1),intent(in) ::
     |  tni,barm
      real,dimension(lev0:lev1,lon0:lon1),intent(out) ::
     |  xnmbari
!
! Local:
      integer :: k,i
      real :: expzi(lev0:lev1,lon0:lon1) ! e(-z) at interfaces
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  TotMassDen
!
! xnmbari = p0*e(-z)*barm/kT at interfaces:
        do i=lon0,lon1
          expzi(lev0,i) = expzmid_inv*expz(lev0)
          do k=lev0+1,lev1-1
            expzi(k,i) = expzmid_inv*expz(k)
          enddo
          expzi(lev1,i) = expzmid*expz(lev1-1)
          do k=lev0,lev1
            xnmbari(k,i) = p0*expzi(k,i)*barm(k,i)/
     |        (boltz*tni(k,i))
!! Converting xnmbar to total mass density [kg/m^3]
!** (mean molecular mass @ interfaces*proton mass * unit conversion [g cm-3] to [kg m-3]) **
            TotMassDen(k,i) = xnmbari(k,i)*1.67E-24*10.**3
          enddo
        enddo
      call addfld('TOTMASDEN',' ',' ',TotMassDen,
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)

      end subroutine mk_xnmbari
!-----------------------------------------------------------------------
      end module qrj_module
