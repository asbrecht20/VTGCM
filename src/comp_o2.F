!
      module o2_module
      use params_module,only: nlevp1,nlonp4,nlat,dz
      implicit none
! -------------------------------------------------------------
! Boundary conditions, production and loss for O2 are defined
! by comp_o2, and referenced by minor_o2. Comp_o2 is called
! from a latitude loop in dynamics. After comp_o2, dynamics calls
! minor_o2, which passes this module data to sub minor. Sub 
! minor contains 3d mpi calls and its own latitude loops.
!
!     real,dimension(nlonp4,nlat)   :: o2_ubc ! upper boundary
!     real,dimension(nlonp4,3,nlat) :: o2_lbc ! lower boundary
!     real,dimension(nlevp1,nlonp4,nlat) ::
!    |  o2_prod,     ! production of o2
!    |  o2_loss      ! loss of o2
!
! Boundary conditions and production and loss terms are allocated
! subdomains by sub alloc_n4s (called from allocdata.F).
! - last modified by SWB:  02/15/07
! 2/16/07 btf: Commenting out the n4s_ubc assignment below
!              until a comp_n4s module is written.
! 5/3/07  swb: slab indices corrected from O2IRINT calculation
! 5/3/07  swb: integral limits fixed.
! 10/7/07 swb: integral limits revised to capture O2 (1-Delta) above ~90 km
!              only  (above zp = -11.0, k = 11)
! 10/18/07 swb: rates are updated accordiing to Slanger et al., (2006) and
!               Gerard et al., (2007)
! 05/06/08 swb: trsolv mod for lev1-1 in call statement
! 05/06/08 swb: qjn2/xo2 (loss freqeuncy) and calc production/loss terms
! 05/08/14 swb: r2(k,i,lat) and r5(k,i,lat) added back
! 12/05/16 swb: add r123 and r130 reactions for src of O2 (leading only to O3 production)
! -------------------------------------------------------------
!
      real,allocatable,dimension(:,:)   :: o2_ubc ! upper boundary (i,j)
      real,allocatable,dimension(:,:,:) :: o2_lbc ! lower boundary (i,3,j)
      real,allocatable,dimension(:,:,:) ::
     |  o2_prod,     ! production of o2 (k,i,j)
     |  o2_loss      ! loss of o2       (k,i,j)
      real :: phi_o2(3) =  (/1.961, 1.7316, 1.3774/)
!     real :: pso2b =  3.23E-04  !  mmr specified at bottom
      real :: pso2b =  6.0E-06   !  Y&D(1982) mmr specified at bottom
      real,parameter :: alfa_o2 = 0. ! thermal diffusion coefficient
      real :: BETA, FO2IR, FO2VIS
!
      contains
!-----------------------------------------------------------------------
      subroutine alloc_o2(lon0,lon1,lat0,lat1)
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
      allocate(o2_ubc(lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_o2: error allocating',
     |  ' o2_ubc: stat=',i3)") istat
      allocate(o2_lbc(lon0:lon1,3,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_o2: error allocating',
     |  ' o2_lbc: stat=',i3)") istat
!
! Allocate subdomains to production and loss:
      allocate(o2_prod(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_o2: error allocating',
     |  ' o2_prod: stat=',i3)") istat
      allocate(o2_loss(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_o2: error allocating',
     |  ' o2_loss: stat=',i3)") istat

!     write(6,"('alloc_o2: allocated module data')")
!
      end subroutine alloc_o2
!
!-----------------------------------------------------------------------
      subroutine comp_o2(tn,o1,co,co2,n2,o2,barm,o2p,op,co2p,ne,
     |  lev0,lev1,lon0,lon1,lat)
!
!-----------------------------------------------------------------------
! Advance o2 (minor): calculate O2-IR and O2-vis nightglows as well
! Arguments for comp_o2 (16).
! ----------------------------------------------------------------------
! Advance o2 by one time step. This is called from driver at
! each subdomain latitude.
! Last modified: swb (02/14/07)
!              : add new rates from chemrates
!              : add o2src, o2loss terms
! ----------------------------------------------------------------------
!
      use qrj_module,only: qjo2, mk_xnmbari
      use fields_module,only: tlbc
      use cons_module,only: rmassinv_o1,rmassinv_co,rmassinv_n2,p0,
     |  expz,expzmid_inv,boltz,rmassinv_co2,rmassinv_o2,rmass_o1,
     |  rmass_co,rmass_o2,gask,grav
      use chemrates_module,only:r2,r5,r123,r130
      use addfld_module,only: addfld
!
!     write(6,"('enter comp_o2: lat=',i2)") lat
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
     |  barm,    ! mean molecular weight
     |  o2p,     ! O2+ ion
     |  op,      ! O+ ion
     |  co2p,    ! CO2++ ion
     |  ne       ! electron density 
!
! ----------------------------------------------------------------------
! Local:
      integer :: k,i,k0,k1,i0,i1,nk  ! for addfld (ask for check)
      real,dimension(lev0:lev1,lon0:lon1) :: 
     |  tni,            ! TN at interfaces
     |  qjo2i,          ! JO2 at interfaces (qrj module): 1/sec units from #/cm3/sec
     |  xnmbari,        ! xnmbar at interfaces (qrj module)
     |  xnmbar,         ! xnmbar at midlevels (qrj module)
     |  xo,             ! O1 number density (#/cm3)
     |  xco,            ! CO number density (#/cm3)
     |  xclo,           ! CLO number density (#/cm3); from KINETICS (clo)
     |  xco2,           ! CO2 number density (#/cm3)
     |  xn2,            ! N2 number density (#/cm3)
     |  xo2,            ! O2 number density (#/cm3)
     |  xtot,           ! total number density (#/cm3)
     |  o2prod,         ! O2 production (#/cm3/s)
     |  o2loss1,        ! O2 loss1 (#/cm3/s)
     |  o2loss2,        ! O2 loss1 (#/cm3/s)
     |  slab,           ! slab width (cm)
     |  o2irvem,        ! o2ir volume emission (photons/cm3/sec)
     |  o2irveml,       ! o2ir volume emission (log10 photons/cm3/sec)
     |  o2irint         ! o2ir integrated intensity (photons/cm2/sec in kRayleighs)
      real,dimension(lev0:lev1) :: colout, fldht, clo 
!
! ----------------------------------------------------------------------
! For addfld calls: Check what is needed!
      k0=lev0
      k1=lev1
      i0=lon0
      i1=lon1
      nk = lev1-1
!
! ----------------------------------------------------------------------
! Lower boundary:
! o2_lbc(:,1)=A, o2_lbc(:,2)=B, o2_lbc(:,3)=C define lower boundary 
!   condition where A*DPSX/DZ + B*PSX + C = 0.
!
      do i=lon0,lon1
!
! Value at bottom given by specified mass missing ratio (start)
        o2_lbc(i,1,lat) = 0.
        o2_lbc(i,2,lat) = 1.
        o2_lbc(i,3,lat) = -pso2b
!
! Upper boundary: Zero diffusive flux at top
!       o2_ubc(i,lat) = 0.
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
!********************************************************************
! Species Number densities at interfaces (#/cm3)
! Given Qrj dissociation rates (#/cm3/sec) at midlevels from qrj
! Calculate Qrj dissociation rates (#/sec) at interfaces: JO2
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          xo(k,i)   = xnmbari(k,i)*o1(k,i)*rmassinv_o1
          xco(k,i)  = xnmbari(k,i)*co(k,i)*rmassinv_co
          xco2(k,i) = xnmbari(k,i)*co2(k,i)*rmassinv_co2
          xo2(k,i) = xnmbari(k,i)*o2(k,i)*rmassinv_o2
	  if (xco2(k,i) < 1.0) xco2(k,i) = 1.0
          xn2(k,i)  = xnmbari(k,i)*n2(k,i)*rmassinv_n2
          qjo2i(k,i) = 0.5*(qjo2(k,i,lat)+qjo2(k+1,i,lat))/xo2(k,i)
!
!   LOG10 DENSITIES FOR PLOTTING 
!         xco2l(k,i)=log10(xco2(k,i))
!         xcol(k,i)=log10(xco(k,i))
!         xn2l(k,i)=log10(xo2(k,i))
!         xol(k,i)=log10(xo(k,i))

! Creating total number density
!         xtot(k,i) = xo(k,i)+xco(k,i)+xco2(k,i)+xo2(k,i)+xn2(k,i)
          xtot(k,i) = xo(k,i)+xco(k,i)+xco2(k,i)+xn2(k,i)

!  Converting vmr to number density
          xclo(k,i)  = xtot(k,i)*clo(k)

        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!
! ----------------------------------------------------------------------
! (1)  Major Sources and Losses for O2 (before NOX Chemistry Added)

! -- Checked against mtgcm16. OK (02/12/07)
! (2) Airglow quantities calculated
! -- O2-IR (Volume Emission Rate, Integrated Vertical Intensity)
! ---------------------------------------------------------------------
      do i=lon0,lon1
        do k=lev0,lev1-1
! ---------------------------------------------------------------------
!     SOURCES (#/cm3.sec) AND SINKS(#/sec):
! ---------------------------------------------------------------------
!
          o2_prod(k,i,lat) = xo(k,i)*xo(k,i)*xco2(k,i)*r2(k,i,lat) +
     |                       xo(k,i)*xclo(k,i)*r123(k,i,lat)+ 
     |                       xclo(k,i)*xclo(k,i)*r130(k,i,lat) 
          o2_loss(k,i,lat) = -(xo(k,i)*xco2(k,i)*r5(k,i,lat) + 
     |	                      qjo2i(k,i))
	  o2prod(k,i) = o2_prod(k,i,lat)
          o2loss1(k,i) = xo(k,i)*xo(k,i)*xco2(k,i)*r5(k,i,lat)
	  o2loss2(k,i) = qjo2i(k,i)*xo2(k,i)
!
! -----------------------------------------------------------------------
!     AIRGLOW QUANTITIES (SLAB WIDTH = 0.5*SCALE HEIGHT; cm units):
!  *  O2 IR DAYGLOW FROM 3-BODY RECOMBINATION (via vtgcm2 coding)
!  *  3-BODY O-RECOMBINATION INTO VARIOUS STATES: 
!     (a) DIRECTLY to O2(a)-state and
!     (b) INDIRECTLY to O2(5pi)-state plus quenching to O2(a)-state
!     (c) INDIRECTLY to O2(c)-state plus quenching to O2(a)-state
! ---------------------------------------------------------------------
!  ** VTGCM2D **  2007-onward
!  *  3-BODY O-RECOMBINATION TO O2(c)-state plus quenching to O2(a)-state
!     FO2VIS = 0.37/(1. + 3.0E-15*xco2(k,i)*50.)
!  *  3-BODY O-RECOMBINATION TO O2(c)-state all quenched to O2(a)-state
!     FO2VIS = 0.37
!  *  3-BODY O-RECOMBINATION (a), (b), (c)  Leading to O2(a)-state 
!     (Direct, Indirect leading to O2(a) state)
!     BETA= 0.025 + 0.07 + 0.20  + FO2VIS
!  *  Quenching of O2(a)-state by CO2 and N2  (very weak)
!     FO2IR=1./(1.+ 3800.*(3.0E-20*xco2(k,i) + 2.0E-20*xn2(k,i)))
!     Photons/cm3.sec
!     o2irvem(k,i)=xo(k,i)*xo(k,i)*xco2(k,i)*r2(k,i,lat)*BETA*FO2IR
!     slab(k,i)=0.5*gask*tni(k,i)/(barm(k,i)*grav)
!  *  -----------------------------------------------------------
!  ** VTGCM2D **  Fall 2007  (Slanger et al., 2006; Gerard et al., 2007)
!  *  3-BODY O-RECOMBINATION (a), (b), (c). NET YIELD  Leading to O2(a)-state
!     (Direct, Indirect leading to O2(a) state)
      BETA= 0.75
!     BETA= 0.7
!     BETA= 1.00
!  *  Quenching of O2(a)-state by CO2 and N2  (very weak)
!  *  Quenching of O2(a)-state by CO (consider this as well?)
      FO2IR=1./(1.+ 3800.*(2.0E-20*xco2(k,i) + 2.0E-20*xn2(k,i)))
!     Photons/cm3.sec
      o2irvem(k,i)=xo(k,i)*xo(k,i)*xco2(k,i)*r2(k,i,lat)*BETA*FO2IR
      slab(k,i)=gask*tni(k,i)/(barm(k,i)*grav)
!  *  -----------------------------------------------------------
! ** Case 1 through 3 for Marie-Eve
!     BETA= 0.67
!     BETA= 0.5
!     BETA= 0.7
!     FO2IR=1./(1.+ 4470.*(5.0E-21*xco2(k,i) + 2.0E-20*xn2(k,i)))
!     FO2IR=1./(1.+ 4545.*(2.0E-20*xco2(k,i) + 2.0E-16*xo(k,i)))
!     FO2IR=1./(1.+ 4460.*(1.0E-20*xco2(k,i) + 2.0E-16*xo(k,i)))
!     Photons/cm3.sec
!     o2irvem(k,i)=xo(k,i)*xo(k,i)*xco2(k,i)*r2(k,i,lat)*BETA*FO2IR
!     o2irvem(k,i)=xo(k,i)*xo(k,i)*xco2(k,i)*r2*BETA*FO2IR
!     slab(k,i)=gask*tni(k,i)/(barm(k,i)*grav)
!  *  -----------------------------------------------------------
!  *  -----------------------------------------------------------
!     Log10 (Photons/cm3.sec)
        if (o2irvem(k,i) .GT. 1.0) THEN
	   o2irveml(k,i)= alog10(o2irvem(k,i))
        else
	   o2irveml(k,i) = 0.
        endif
!  *  -----------------------------------------------------------
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
! ----------------------------------------------------------------------
!   * DISPOSE O2(1-delta) INTEGRATED EMISSION RATE (PH/cm2/sec)
      do i=lon0,lon1
      if (dz == 0.5) then
!  Below 90 km is zeroed out
        do k=1,11
             fldht(k) = 0.
        enddo
!  Above 90 km is retained 
        do k=11,lev1-1
             fldht(k) = slab(k,i)*o2irvem(k,i)
        enddo
      elseif (dz == 0.25) then
!  Below 90 km is zeroed out
        do k=1,21
             fldht(k) = 0.
        enddo
!  Above 90 km is retained 
        do k=21,lev1-1
             fldht(k) = slab(k,i)*o2irvem(k,i)
        enddo
      endif
        call integral(colout,fldht,1,dz,nk)
        do k = lev0,lev1-1
!    * KILO-RAYLEIGH UNITS FOR PLOTS (O2 IR INTEGRATED INTENSITY)
             o2irint(k,i) = colout(k)/1.0E+9
        enddo
      enddo
! ----------------------------------------------------------------------
!   * QJOT at Interfaces
!     call addfld('QJO2I',' ',' ',qjo2i,'lev',k0,k1,'lon',i0,i1,lat)
!   * O2(1-delta) VOLUME EMISSION RATE (log10 PH/cm3/sec)
      call addfld('O2IRVEM',' ',' ',o2irveml,'lev',k0,k1,'lon',i0,i1,
     |  lat)
!   * O2(1-delta) INTEGRATED VERTICAL INTENSITY (kR = 1.0E+09*PH/cm2/sec)
      call addfld('O2IRINT',' ',' ',o2irint,'lev',k0,k1,'lon',i0,i1,lat)
! ----------------------------------------------------------------------
!   * O2 Production and Losses
      call addfld('O2PROD',' ',' ',o2prod,'lev',k0,k1,'lon',i0,i1,
     |  lat)
      call addfld('O2LOSS1',' ',' ',o2loss1,'lev',k0,k1,'lon',i0,i1,
     |  lat)
      call addfld('O2LOSS2',' ',' ',o2loss2,'lev',k0,k1,'lon',i0,i1,
     |  lat)
!
      end subroutine comp_o2
!-----------------------------------------------------------------------
      subroutine minor_o2(tn,difk,o1,co,co2,n2,o2,o2_nm,o2_out,
     |  o2_nm_out,lev0,lev1,lon0,lon1,lat0,lat1)

!     Arguments (16): o1, co, co2, n2 needed by minor.F for barm, etc.
      use cons_module,only: rmass_o2
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
     |  o2,      ! molecular oxygen (mmr)
     |  o2_nm    ! O2 at time n-1 
!
! -----------------------------------------------------------------------
! Output args:
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(out) ::
     |  o2_out,    ! O2 output 
     |  o2_nm_out  ! O2 output at time n-1
!
! -----------------------------------------------------------------------
! Local:
      integer :: lat
      integer :: i0,i1,nk ! for addfld
!
!     write(6,"('enter minor_o2')")
      i0 = lon0
      i1 = lon1
      nk = lev1-1
! -----------------------------------------------------------------------
! Minor returns o2_out and o2_nm_out. Module data o2_prod,
! o2_loss, etc, were defined by comp_o2.  
!
!     subroutine minor(tn,difk,o2,o1,fcomp,fcomp_tm1,fcomp_out,
!    |  fcomp_tm1_out,sloss,sprod,flbc,fubc,rmx,phix,alfax,
!    |  lev0,lev1,lon0,lon1,lat0,lat1,idebug)
!
      call minor(tn,difk,o1,co,co2,n2,o2,o2_nm,o2_out,o2_nm_out,
     |  o2_loss,o2_prod,o2_lbc,o2_ubc,rmass_o2,phi_o2,
     |  alfa_o2,lev0,lev1,lon0,lon1,lat0,lat1,0)

!       call addfld('O2_OUT' ,' ',' ',o2_out(:,i0:i1,lat),
!    |    i0,i1,k0,k1,lat)
!       call addfld('O2_TM1' ,' ',' ',o2_nm_out(:,i0:i1,lat),
!    |    i0,i1,k0,k1,lat)
      end subroutine minor_o2
!
! -----------------------------------------------------------------------
      subroutine integral(FN,Y,J0,delz,JEND) 
!COMPUTES INTEGRAL OF Y FROM Z(J0) TO Z(J) AND RETURNS INTEGRAL IN FN(J)
! USES 4TH ORDER ADAMS MOULTON ALGORITHM (CF HENRICI,P199)            
! Input args:
      integer,intent(in) :: j0, jend
      real,dimension(jend),intent(out) :: FN
      real,dimension(jend),intent(in) :: Y
      real,intent(in) :: delz
! Local:
      integer :: n, j, jud, k
!
      FN(J0)=0                                                       
      N=IABS(J0-JEND)+1                                             
      JUD=(JEND-J0)/IABS(JEND-J0)                                  
! JEND IS THE LARGEST J ABOVE J0 OR THE SMALLEST J BELOW J0       
      DO 1 K=1,N                                                 
      J=K+J0-1                                                  
      IF(JUD.EQ.-1) J=N+1-K                                    
      IF(K.LT.4.AND.JUD.EQ.1) GO TO 2                         
      IF(N-J.LT.3.AND.JUD.EQ.-1) GO TO 6                     
      FN(J)=FN(J-JUD)+(9.*Y(J)+19.*Y(J-JUD)-5.*Y(J-2*JUD)+Y(J-3*JUD))*  
     1delz/24                                                           
      GO TO 1                                                        
    2 IF(2-K)3,4,5                                                  
    3 FN(J)=FN(J-1)+(5.*Y(J)+8.*Y(J-1)-Y(J-2))*delz/12.              
      GO TO 1                                                     
    4 FN(J)=FN(J-1)+(Y(J)+Y(J-1))*delz/2.                          
      GO TO 1                                                   
    5 FN(J)=0                                                  
      GO TO 1                                                 
    6 IF(N-J-1) 7,8,9                                        
    7 FN(J)=0                                               
      GO TO 1                                              
    8 FN(J)=FN(J+1)+(Y(J)+Y(J+1))*delz/2.                   
      GO TO 1                                            
    9 FN(J)=FN(J+1)+(5.*Y(J)+8.*Y(J+1)-Y(J+2))*delz/12.   
    1 CONTINUE                                         
   10 CONTINUE                                        
      end subroutine integral
      end module o2_module
