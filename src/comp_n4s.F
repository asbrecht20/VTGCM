!
      module n4s_module
      use params_module,only: nlevp1,nlonp4,nlat,dz
      implicit none
! -------------------------------------------------------------
! Boundary conditions, production and loss for N4S are defined
! by comp_n4s, and referenced by minor_n4s. Comp_n4s is called
! from a latitude loop in dynamics. After comp_n4s, dynamics calls
! minor_n4s, which passes this module data to sub minor. Sub 
! minor contains 3d mpi calls and its own latitude loops.
!
!     real,dimension(nlonp4,nlat)   :: n4s_ubc ! upper boundary
!     real,dimension(nlonp4,3,nlat) :: n4s_lbc ! lower boundary
!     real,dimension(nlevp1,nlonp4,nlat) ::
!    |  n4s_prod,     ! production of n4s
!    |  n4s_loss      ! loss of n4s
!
! Boundary conditions and production and loss terms are allocated
! subdomains by sub alloc_n4s (called from allocdata.F).
! 10/11/07 swb: n4s, n2d, revised qrj.nox added
! 03/12/08 swb: n4s, n2d, revised  with latest rates (Gerard et al., 2008)
! 03/17/08 swb: n4s, n2d, revised  no separate day-night sources or losses
! 03/31/08 swb: n4s: r29(k,i,lat) usage in n4s_prod(k,i,lat)
! 04/01/08 swb: n4s: put out log10(NO) for secondary history.
! 04/16/08 swb: n4s: add N2+ reactions from Gerard et al. (1988)
! 04/17/08 swb: n4s: add diagnostic terms for NOx  (both N4S and N2D)
!                  ; simplified chemistry for major reactions
! 04/22/08 swb: n4s: N2D zeroed out below z - -8.0 (dayside)
!                  : N2D dayside and nightside prod/loss together
!                  : zero flux lbc
! 04/29/08 swb: n4s: Secondary/primary comparions (NO, N2D, N4S, ions)
!                  : Same timesteps.
! 04/29/08 swb: n4s: new tiegcm lbc tested (explicit PCE)
! 05/07/08 swb: n4s: Remove loss #4 and 5 and src #5
! 05/07/08 swb: n4s: Include O2+ loss of NO
! 05/12/08 swb: no : add C-D equation solution for NO
! 05/13/08 swb: N4S: add back O+N+CO2 reaction for N4S loss  (rg2)
! 05/14/08 swb: N4S: psn4s = 0 used for lbc
! 05/19/08 swb: NO : phi_no(3) values reasonably estimated for comp_no.F
! 06/02/08 swb: NO : phi_no(3) = large for comp_no.F  (weak NO diffusion)
! 06/02/08 swb: NO : XNOPCE for N4S loss term.
! 06/03/08 swb: NO : Add PCE Nite NO sources:  N+O-> NO*; N+O+CO2-> NO+CO2
! 06/04/08 swb: NO : Add CD Nite NO sources:  N+O-> NO*; N+O+CO2-> NO+CO2
! 06/11/08 swb: N4S : corr 2.0*(1.-gg) and 2.0*gg terms for N4S & N3D prod
! 07/11/08 swb: N4S : ff=gg=0.5
! 03/31/09 swb: N4S : bi-molecular diffusion coeff. upd (Dn-o, Dn-co2)
! 03/31/09 swb: N4S : bi-molecular diffusion coeff. upd (Dno-o, Dno-co, Dno-co2)
! -------------------------------------------------------------
!
!  Dimension both N4S and NO arrays
      real,allocatable,dimension(:,:)   :: n4s_ubc ! upper boundary (i,j)
      real,allocatable,dimension(:,:,:) :: n4s_lbc ! lower boundary (i,3,j)
      real,allocatable,dimension(:,:)   :: no_ubc ! upper boundary (i,j)
      real,allocatable,dimension(:,:,:) :: no_lbc ! lower boundary (i,3,j)
      real,allocatable,dimension(:,:,:) ::
     |  n4s_prod,     ! production of n4s (k,i,j)
     |  n4s_loss,     ! loss of n4s       (k,i,j)
     |  no_prod,      ! production of no  (k,i,j)
     |  no_loss       ! loss of no        (k,i,j)
!  N4S and NO Bi-molecular diffusion paramaters
      real :: phi_n4s(3) =  (/1.144, 1.310, 0.6200/) ! After Gerard et al (80) and E-TIGCM
!     real :: phi_n4s(3) =  (/0.315, 1.310, 0.2517/) ! Partially bad estimates
!     real :: phi_no(3) =  (/0.866, 0.814, 0.926/) ! Bad first guess
!     real :: phi_no(3) =  (/1.361, 1.737, 1.280/) ! Estmated from Earth and scaled
!     real :: phi_no(3) =  (/27.5, 15.7, 10.0/) ! Weak NO diffusion! (Dij=D0*0.1)
      real :: phi_no(3) =  (/1.144, 1.310, 0.6200/) ! Assume same coefficients as N4S
!  N4S paramaters:
      real,parameter :: psn4s =  0.0  !  mass mixing ratio specified at bottom
      real,parameter :: alfa_n4s = 0. ! thermal diffusion coefficient
      real,parameter :: ff = 0.50 !  branching n4s to n2d for JN2(eff)
!     real,parameter :: gg = 0.85 !  branching n4s to n2d for NO+ DR
      real,parameter :: gg = 0.75 !  branching n4s to n2d for NO+ DR
!     real,parameter :: gg = 0.50 !  branching n4s to n2d for NO+ DR
      real,parameter :: fnod01 = 0.2464 !  yield for NO(0,1)-d band
      real,parameter :: fnoall = 1.0    !  yield for all NO* Gamma and Delta bands
!  NO Paramaters  (tiegcm values to start; to be rechecked for Venus)
      real,parameter :: psno =  0.0  !  mass mixing ratio specified at bottom
      real,parameter :: alfa_no = 0. ! thermal diffusion coefficient
!
      contains
!-----------------------------------------------------------------------
      subroutine alloc_n4s(lon0,lon1,lat0,lat1)
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
      allocate(n4s_ubc(lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_n4s: error allocating',
     |  ' n4s_ubc: stat=',i3)") istat
      allocate(n4s_lbc(lon0:lon1,3,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_n4s: error allocating',
     |  ' n4s_lbc: stat=',i3)") istat
!
! Allocate subdomains to production and loss:
      allocate(n4s_prod(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_n4s: error allocating',
     |  ' n4s_prod: stat=',i3)") istat
      allocate(n4s_loss(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_n4s: error allocating',
     |  ' n4s_loss: stat=',i3)") istat

!     write(6,"('alloc_n4s: allocated module data')")
!
      end subroutine alloc_n4s
!-----------------------------------------------------------------------
      subroutine alloc_no(lon0,lon1,lat0,lat1)
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
      allocate(no_ubc(lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_no: error allocating',
     |  ' no_ubc: stat=',i3)") istat
      allocate(no_lbc(lon0:lon1,3,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_no: error allocating',
     |  ' no_lbc: stat=',i3)") istat
!
! Allocate subdomains to production and loss:
      allocate(no_prod(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_no: error allocating',
     |  ' no_prod: stat=',i3)") istat
      allocate(no_loss(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_no: error allocating',
     |  ' no_loss: stat=',i3)") istat

!     write(6,"('alloc_no: allocated module data')")
!
      end subroutine alloc_no
!
!-----------------------------------------------------------------------
      subroutine comp_n4s(tn,o1,co,co2,n2,barm,
     |  n4s,n2d,no,o2p,op,nop,
     |  n2p,ne,lev0,lev1,lon0,lon1,lat)
!
!-----------------------------------------------------------------------
! Advance n4s (minor): calculate NO-UV nightglows as well
! Arguments for comp_n4s (19).
! ----------------------------------------------------------------------
! Advance n4s by one time step. This is called from driver at
! each subdomain latitude.
! Last modified: swb (10/09/07)
!              : add new rates from chemrates
!              : add n4s_prod, n4s_loss terms
! ----------------------------------------------------------------------
!
      use qrj_module,only: qjn2, mk_xnmbari
      use fields_module,only: tlbc
      use cons_module,only: rmassinv_o1,rmassinv_co,rmassinv_n2,p0,
     |  expz,expzmid_inv,boltz,rmassinv_co2,rmassinv_n4s,rmass_o1,
     |  rmass_co,rmass_n2,rmass_co2,rmass_n4s,gask,grav,
     |  rmass_no,rmassinv_no
      use chemrates_module,only:rk1,rk3,rk4,rk5,rk10,r28,r29,
     | a15,a17,rg1,rg2,rg3
      use addfld_module,only: addfld
      use o2_module,only: integral
      use chapman_module,only: idn_lonlat
!
!     write(6,"('enter comp_n4s: lat=',i2)") lat
! ----------------------------------------------------------------------
! Input args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) ::
     |  tn,      ! neutral temperature (deg K)
     |  o1,      ! atomic oxygen (mmr)
     |  co,      ! carbon monoxide (mmr)
     |  co2,     ! carbon dioxide (mmr)
     |  n2,      ! molecular nitrogen (mmr)
     |  barm,    ! mean molecular weight
     |  n4s,     ! atomic nitrogen (mmr)
     |  n2d,     ! atomic nitrogen excited (mmr)
     |  no,      ! nitric oxide (mmr)
     |  o2p,     ! O2+ ion (#/cm3)
     |  op,      ! O+ ion  (#/cm3)
     |  nop,     ! NO+ ion  (#/cm3)
     |  n2p,     ! N2+ ion  (#/cm3)
     |  ne       ! electron density  (#/cm3)
!
! ----------------------------------------------------------------------
! qjn2i:  combined N2 photodissociation + PE-dissociation rate (#/cm3/sec)
! Local:
      integer :: k,i,k0,k1,i0,i1,nk  ! for addfld (ask for check)
      real,dimension(lev0:lev1,lon0:lon1) :: 
     |  tni,            ! TN at interfaces
     |  qjn2i,          ! JN2EFF at interfaces (qrj module) (#/cm3/sec)
     |  xnmbari,        ! xnmbar at interfaces (qrj module)
     |  xnmbar,         ! xnmbar at midlevel (qrj module)
     |  xo,             ! O1 number density (#/cm3)
     |  xco,            ! CO number density (#/cm3)
     |  xco2,           ! CO2 number density (#/cm3)
     |  xn2,            ! N2 number density (#/cm3)
     |  xn4s,           ! N4S number density (#/cm3)
     |  xn2d,           ! N2D number density (#/cm3)
     |  xno,            ! NO number density (#/cm3); Calc. in comp_no
     |  xnopce,         ! NO number density (#/cm3): PCE Internal/secondary
     |  xnopcel,        ! NO number density (log10 #/cm3): PCE secondary
     |  xnpce,          ! N4S number density (#/cm3): PCE secondary
     |  xnpcel,         ! N4S number density (log10 #/cm3): PCE secondary
     |  slab,           ! slab width (cm)
     |  novem,          ! no volume emission (photons/cm3/sec)
     |  noveml,         ! no volume emission (log10 photons/cm3/sec)
     |  nprodcol,       ! n4s production integrated column (#/cm2/sec)
     |  nnetprodcol,    ! net (src - loss) n4s production integrated column (#/cm2/sec)
     |  nnet,           ! 
     |  ndestcol,       ! n4s destruction integrated column (#/cm2/sec)
     |  noint           ! no integrated intensity (photons/cm2/sec in kR)
!
!
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  taun4s,            ! 1/TOTAL N(4S) loss frequncy (Earth Days)
     |  nlossT,            ! N(4S) loss term (3 major reactions)
     |  nloss1,            ! N(4S) loss term 1
     |  nloss2,            ! N(4S) loss term 2
     |  nloss3,            ! N(4S) loss term 3
     |  nloss4,            ! N(4S) loss term 4
!    |  n4s_dayloss,
!    |  n4s_nightloss,
     |  n4ssrc,
     |  n4sloss,
     |  nprod1,            ! N(4S) prod term 1
     |  nprod2,            ! N(4S) prod term 2
     |  nprod3,            ! N(4S) prod term 3
     |  nprod4,            ! N(4S) prod term 4
     |  nprod6,            ! N(4S) prod term 6
     |  nprodtot,          ! N(4S) total prod term 
     |  xne                ! electron density  (#/cm3)
!
      real,dimension(lev0:lev1) :: colout, fldht 
      real :: xnmbar_lbc
!
! ----------------------------------------------------------------------
! For addfld calls: Check what is needed!
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
! Qrj dissociation rates (#/cm3.sec) at interfaces: JO2
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          xo(k,i)   = xnmbari(k,i)*o1(k,i)*rmassinv_o1
          xco(k,i)  = xnmbari(k,i)*co(k,i)*rmassinv_co
          xco2(k,i) = xnmbari(k,i)*co2(k,i)*rmassinv_co2
	  if (xco2(k,i) < 1.0) xco2(k,i) = 1.0
          xn2(k,i)  = xnmbari(k,i)*n2(k,i)*rmassinv_n2
          xn4s(k,i)  = xnmbari(k,i)*n4s(k,i)*rmassinv_n4s
          xn2d(k,i)  = xnmbari(k,i)*n2d(k,i)*rmassinv_n4s
          xno(k,i)  = xnmbari(k,i)*no(k,i)*rmassinv_no
          qjn2i(k,i) = 0.5*(qjn2(k,i,lat)+qjn2(k+1,i,lat))
          xne(k,i)   = ne(k,i)
!
!   LOG10 DENSITIES FOR PLOTTING 
!         xco2l(k,i)=log10(xco2(k,i))
!         xcol(k,i)=log10(xco(k,i))
!         xn2l(k,i)=log10(xo2(k,i))
!         xol(k,i)=log10(xo(k,i))

        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
! ----------------------------------------------------------------------
! (1)  PCE: Major Sources and Losses for no 
! ---------------------------------------------------------------------
      do i=lon0,lon1
        do k=lev0,lev1-1
! ---------------------------------------------------------------------
!     SOURCES (#/cm3.sec) and SINKS(#/sec): for NO (#/cm3.sec
! ---------------------------------------------------------------------
! Simplified (1-production; 2-losses)
!         xnopce(k,i) = xn2d(k,i)*xco2(k,i)*rk1/
!    |                  (xn4s(k,i)*rg3(k,i,lat)+ o2p(k,i)*rk10)
          xnopce(k,i) = (xn2d(k,i)*xco2(k,i)*rk1+
     |                  xn4s(k,i)*xo(k,i)*rg1(k,i,lat)+
!    |                  xn4s(k,i)*xo(k,i)*rg1+
     |                  xn4s(k,i)*xo(k,i)*xco2(k,i)*rg2(k,i,lat)) 
!    |                  xn4s(k,i)*xo(k,i)*xco2(k,i)*rg2) 
     |                  /(xn4s(k,i)*rg3(k,i,lat)+ o2p(k,i)*rk10)
!    |                  /(xn4s(k,i)*rg3 + o2p(k,i)*rk10)
          if (xnopce(k,i) <= 1.0) then 
              xnopcel(k,i) = 0.0
          else
              xnopcel(k,i) = log10(xnopce(k,i))
          endif
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
! ----------------------------------------------------------------------
! (1)  Major Sources and Losses for n4s (NOX Chemistry Added)
! (2)  Airglow quantities calculated
! -- NO-UV (Volume Emission Rate, Integrated Vertical Intensity)
! ---------------------------------------------------------------------
      do i=lon0,lon1
        do k=lev0,lev1-1
! ---------------------------------------------------------------------
!     SOURCES (#/cm3.sec) and SINKS(#/sec): N4S (Venus day/night together)
!  ---Nightside  (No local production of n4s; simplified n4s_loss terms)
!  -- Combined, JN2-dissociation + PE dossociation scheme used 
! --- Combined jn2eff calculated in qrj.F
! ---------------------------------------------------------------------
!    Day-Night split
!         if (idn_lonlat(i,lat)==1) then
!    Day
               n4s_prod(k,i,lat) = 2.*(1-ff)*qjn2i(k,i)  +
     |	         rk3*xco(k,i)*xn2d(k,i)+rk4*xo(k,i)*xn2d(k,i) + 
     |	         rk5*xn2(k,i)*xn2d(k,i)+
     |	         a17(k,i,lat)*nop(k,i)*ne(k,i)*(1.-gg)*2.0 
!              n4ssrc(k,i) = n4s_prod(k,i,lat)
!
!!   |	         a17(k,i,lat)*nop(k,i)*ne(k,i)*(1.-gg) 
!!             n4s_loss(k,i,lat) = -o2p(k,i)*r28 -
!!   |            xno(k,i)*rg3(k,i,lat) - 
!!   |	          xo(k,i)*rg1(k,i,lat) -
!!   |            xo(k,i)*xco2(k,i)*rg2(k,i,lat) 
!              n4s_loss(k,i,lat) = -o2p(k,i)*r28 -
!    |            xnopce(k,i)*rg3(k,i,lat)
!!             n4s_loss(k,i,lat) = 2.0*(-o2p(k,i)*r28  
!!   |	         -xo(k,i)*rg1(k,i,lat) - xo(k,i)*xco2(k,i)*rg2(k,i,lat)) 
!              n4s_dayloss(k,i) = n4s_loss(k,i,lat)
!         else
!    Night
!              n4s_prod(k,i,lat) = 2.*(1-ff)*qjn2i(k,i)  +
!    |	         rk3*xco(k,i)*xn2d(k,i)+rk4*xo(k,i)*xn2d(k,i) + 
!    |	         rk5*xn2(k,i)*xn2d(k,i)+
!    |	         a17(k,i,lat)*nop(k,i)*ne(k,i)*(1.-gg)*2.0 
!
               n4s_loss(k,i,lat) = -o2p(k,i)*r28 -
     |           xnopce(k,i)*rg3(k,i,lat) - 
!    |           xnopce(k,i)*rg3 - 
     |	         xo(k,i)*rg1(k,i,lat) -
     |           xo(k,i)*xco2(k,i)*rg2(k,i,lat) 
!    |           xo(k,i)*xco2(k,i)*rg2 
               n4sloss(k,i) = n4s_loss(k,i,lat)
!         endif
!
          xnpce(k,i) = -n4s_prod(k,i,lat)/n4s_loss(k,i,lat)
          if (xnpce(k,i) <= 1.0) then 
              xnpcel(k,i) = 0.0
          else
              xnpcel(k,i) = log10(xnpce(k,i))
          endif
          nprod1(k,i) = 2.*(1-ff)*qjn2i(k,i)
          nprod2(k,i) = rk3*xco(k,i)*xn2d(k,i)
          nprod3(k,i) = rk4*xo(k,i)*xn2d(k,i)
          nprod4(k,i) = rk5*xn2(k,i)*xn2d(k,i)
          nprod6(k,i) = a17(k,i,lat)*nop(k,i)*ne(k,i)*(1.-gg)*2.0
!         nprod7(k,i) = a15(k,i,lat)*n2p(k,i)*ne(k,i)*0.5
!
          n4ssrc(k,i) = n4s_prod(k,i,lat)
          nprodtot(k,i) = nprod1(k,i)+nprod2(k,i)+nprod3(k,i)
     |         +nprod4(k,i)+nprod6(k,i)
!
          nlossT(k,i) = xn4s(k,i)*o2p(k,i)*r28
     |        +xn4s(k,i)*xnopce(k,i)*rg3(k,i,lat)
     |        +xn4s(k,i)*xo(k,i)*rg1(k,i,lat) 
     |        +xn4s(k,i)*xo(k,i)*xco2(k,i)*rg2(k,i,lat)
! 
          nloss1(k,i) = xn4s(k,i)*o2p(k,i)*r28
!         nloss2(k,i) = xn4s(k,i)*xno(k,i)*rg3(k,i,lat)
          nloss2(k,i) = xn4s(k,i)*xnopce(k,i)*rg3(k,i,lat)
!         nloss2(k,i) = xn4s(k,i)*xnopce(k,i)*rg3
          nloss3(k,i) = xn4s(k,i)*xo(k,i)*rg1(k,i,lat)
!         nloss3(k,i) = xn4s(k,i)*xo(k,i)*rg1
          nloss4(k,i) = xn4s(k,i)*xo(k,i)*xco2(k,i)*rg2(k,i,lat)
!         nloss4(k,i) = xn4s(k,i)*xo(k,i)*xco2(k,i)*rg2
          taun4s(k,i) = -(1./n4sloss(k,i))/86400.
          nnet(k,i) = n4ssrc(k,i) - (nlossT(k,i))
!
! -----------------------------------------------------------------------
!  ** VTGCM2D **  2007-onward
!     NO AIRGLOW QUANTITIES (SLAB WIDTH = 0.5*SCALE HEIGHT; cm units):
!     -- (0,1)-delta band emission yield specified (Bougher et al., 1990)
!     -- fnod01 = 0.2464  (single (0,1)-delta band)
!     -- fnoall = 1.0     (all delta and gamma  bands : 190.-260. nm)
! ---------------------------------------------------------------------
!  ** Stewart et al., (1980) chemiluminescent process (VEM = photons/cm3/sec)
!     novem(k,i) = fnod01*xo(k,i)*xn4s(k,i)*rg1(k,i,lat)
      novem(k,i) = fnoall*xo(k,i)*xn4s(k,i)*rg1(k,i,lat)
!     novem(k,i) = fnoall*xo(k,i)*xn4s(k,i)*rg1
      slab(k,i)=gask*tni(k,i)/(barm(k,i)*grav)
!  *  -----------------------------------------------------------
!     Log10 (Photons/cm3.sec)
        if (novem(k,i) .GT. 1.0) THEN
	   noveml(k,i)= alog10(novem(k,i))
        else
	   noveml(k,i) = 0.
        endif
!  *  -----------------------------------------------------------
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
! ----------------------------------------------------------------------
!
!   * N4S TOTAL production INTEGRATED COLUMN RATE (#/cm2/sec)
      do i=lon0,lon1
      if (dz == 0.5) then
!  Below 90 km is zeroed out  (because of O-lbc)
         do k=1,11
           fldht(k) = 0.
         enddo
!  Above 90 km is retained
         do k=11,lev1-1
             fldht(k) = slab(k,i)*n4ssrc(k,i)
         enddo
      elseif (dz == 0.25) then
!  Below 90 km is zeroed out  (because of O-lbc)
         do k=1,21
           fldht(k) = 0.
         enddo
!  Above 90 km is retained
         do k=21,lev1-1
             fldht(k) = slab(k,i)*n4ssrc(k,i)
         enddo
      endif
         call integral(colout,fldht,1,dz,nk)
         do k = lev0,lev1-1
             nprodcol(k,i) = colout(k)
         enddo
      enddo
!   * N4S NET production INTEGRATED COLUMN RATE (#/cm2/sec)
      do i=lon0,lon1
      if (dz == 0.5) then
!  Below 90 km is zeroed out  (because of O-lbc)
         do k=1,11
           fldht(k) = 0.
         enddo
!  Above 90 km is retained
         do k=11,lev1-1
             fldht(k) = slab(k,i)*nnet(k,i)
         enddo
      elseif (dz == 0.25) then
!  Below 90 km is zeroed out  (because of O-lbc)
         do k=1,21
           fldht(k) = 0.
         enddo
!  Above 90 km is retained
         do k=21,lev1-1
             fldht(k) = slab(k,i)*nnet(k,i)
         enddo
      endif
         call integral(colout,fldht,1,dz,nk)
         do k = lev0,lev1-1
             nnetprodcol(k,i) = colout(k)
         enddo
      enddo
!   * N4S destruction INTEGRATED COLUMN RATE (#/cm2/sec)
      do i=lon0,lon1
      if (dz == 0.5) then
!  Below 90 km is zeroed out  (because of O-lbc)
         do k=1,11
           fldht(k) = 0.
         enddo
!  Above 90 km is retained
         do k=11,lev1-1
             fldht(k) = slab(k,i)*(nlossT(k,i))
         enddo
      elseif (dz == 0.25) then
!  Below 90 km is zeroed out  (because of O-lbc)
         do k=1,21
           fldht(k) = 0.
         enddo
!  Above 90 km is retained
         do k=21,lev1-1
             fldht(k) = slab(k,i)*(nlossT(k,i))
         enddo
      endif
         call integral(colout,fldht,1,dz,nk)
         do k = lev0,lev1-1
             ndestcol(k,i) = colout(k)
         enddo
      enddo
! ----------------------------------------------------------------------
!
! Lower boundary:
! n4s_lbc(:,1)=A, n4s_lbc(:,2)=B, n4s_lbc(:,3)=C define lower boundary 
!   condition where A*DPSX/DZ + B*PSX + C = 0.
!
      do i=lon0,lon1
      xnmbar_lbc = p0*expz(1)*expzmid_inv*barm(lev0,i)/
     |             (boltz*tlbc(i,lat))
!
! Value at bottom given by specified PCE mass missing ratio (new)
! Value at bottom given by specified mass missing ratio (old)
        n4s_lbc(i,1,lat) = 0.
        n4s_lbc(i,2,lat) = 1.
!       n4s_lbc(i,3,lat) = -rmass_n4s/xnmbar_lbc*xnpce(1,i)
        n4s_lbc(i,3,lat) = -psn4s
! Value at bottom given by specified flux  (zero to start) 
!       n4s_lbc(i,1,lat) = 1.
!       n4s_lbc(i,2,lat) = 0.
!       n4s_lbc(i,3,lat) = 0.
!
! Upper boundary: Zero diffusive flux at top
        n4s_ubc(i,lat) = 0.
!
      enddo ! i=lon0,lon1
! ----------------------------------------------------------------------
!   * NO-UV (0,1)-delta band INTEGRATED EMISSION RATE (PH/cm2/sec)
!   * NO-UV all Gamma and Delta bands INTEGRATED EMISSION RATE (PH/cm2/sec)
      do i=lon0,lon1
      if (dz == 0.5) then
!  Below 90 km is zeroed out  (because of O-lbc)
         do k=1,11
	   fldht(k) = 0.
         enddo
!  Above 90 km is retained
         do k=11,lev1-1
             fldht(k) = slab(k,i)*novem(k,i)
         enddo
      elseif (dz == 0.25) then
!  Below 90 km is zeroed out  (because of O-lbc)
         do k=1,21
	   fldht(k) = 0.
         enddo
!  Above 90 km is retained
         do k=21,lev1-1
             fldht(k) = slab(k,i)*novem(k,i)
         enddo
      endif
         call integral(colout,fldht,1,dz,nk)
         do k = lev0,lev1-1
!    * KILO-RAYLEIGH UNITS FOR PLOTS (NO UV INTEGRATED INTENSITY)
             noint(k,i) = colout(k)/1.0E+9
         enddo
      enddo
! ----------------------------------------------------------------------
!   * NO(0,1)-delta band VOLUME EMISSION RATE (log10 PH/cm3/sec)
      call addfld('NOVEM',' ',' ',noveml,'lev',k0,k1,'lon',i0,i1,
     |  lat)
!   * NO(0,1)-delta band INTEGRATED VERTICAL INTENSITY (kR units))
      call addfld('NOINT',' ',' ',noint,'lev',k0,k1,'lon',i0,i1,lat)
!   * N(4S) reaction rates (terms) 
!     call addfld('N4SDAYLOSS',' ',' ',n4s_dayloss(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('N4SNIGHTLOSS',' ',' ',n4s_nightloss(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('NPROD1',' ',' ',nprod1(:,lon0:lon1),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('NPROD2',' ',' ',nprod2(:,lon0:lon1),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('NPROD3',' ',' ',nprod3(:,lon0:lon1),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('NPROD4',' ',' ',nprod4(:,lon0:lon1),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('NPROD6',' ',' ',nprod6(:,lon0:lon1),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('NPRODTOT',' ',' ',nprodtot(:,lon0:lon1),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('TAUN4S',' ',' ',taun4s(:,lon0:lon1),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('NLOSS1',' ',' ',nloss1(:,lon0:lon1),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('NLOSS2',' ',' ',nloss2(:,lon0:lon1),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('NLOSS3',' ',' ',nloss3(:,lon0:lon1),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('NLOSS4',' ',' ',nloss4(:,lon0:lon1),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('NLOSST',' ',' ',nlossT(:,lon0:lon1),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!   * NO Density (PCE): Full Moved to comp_no.F 
!   * NO Density (PCE): PCE check done here
      call addfld('NOPCE2',' ',' ',xnopcel(:,lon0:lon1),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!   * N4S Density (PCE) 
      call addfld('N4SPCE',' ',' ',xnpcel(:,lon0:lon1),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!   * Check on Ion Densities (Ne only)
!     call addfld('X1NE',' ',' ',xne(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!  *  -----------------------------------------------------------
!   * N4S Production INTEGRATED COLUMN RATE (#/cm2/sec units))
      call addfld('NPRODCOL',' ',' ',nprodcol(:,lon0:lon1),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('NNETPRODCOL',' ',' ',nnetprodcol(:,lon0:lon1),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('NDESTCOL',' ',' ',ndestcol(:,lon0:lon1),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!  Check indices for addfld?
!   * N4S Production
!     call addfld('N4S_PROD',' ',' ',n4s_prod(:,i0:i1,lat),
!    |  i0,i1,k0,k1,lat)
!   * N4S Loss Frequency
!    call addfld('N4S_LOSS',' ',' ',n4s_loss(:,i0:i1,lat),
!    |  i0,i1,k0,k1,lat)
!
      end subroutine comp_n4s
!-----------------------------------------------------------------------
      subroutine minor_n4s(tn,difk,o1,co,co2,n2,n4s,n4s_nm,n4s_out,
     |  n4s_nm_out,lev0,lev1,lon0,lon1,lat0,lat1)

!     Arguments (16): o1, co, co2, n2 needed by minor.F for barm, etc.
      use cons_module,only: rmass_n4s
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
     |  n4s,     ! atomic nitrogen (mmr)
     |  n4s_nm   ! atomic nitrogen at time n-1 
!
! -----------------------------------------------------------------------
! Output args:
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(out) ::
     |  n4s_out,   ! N4S output 
     |  n4s_nm_out ! N4S output at time n-1
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
! Minor returns n4s_out and n4s_nm_out. Module data n4s_prod,
! n4s_loss, etc, were defined by comp_n4s.  
!
!     call minor(tn,difk,o1,co,co2,n2,o2,o2_nm,o2_out,o2_nm_out,
!    |  o2_loss,o2_prod,o2_lbc,o2_ubc,rmass_o2,phi_o2,
!    |  alfa_o2,lev0,lev1,lon0,lon1,lat0,lat1,0)
!
      call minor(tn,difk,o1,co,co2,n2,n4s,n4s_nm,n4s_out,n4s_nm_out,
     |  n4s_loss,n4s_prod,n4s_lbc,n4s_ubc,rmass_n4s,phi_n4s,
     |  alfa_n4s,lev0,lev1,lon0,lon1,lat0,lat1,0)

!       call addfld('N4S_OUT' ,' ',' ',n4s_out(:,i0:i1,lat),
!    |    i0,i1,k0,k1,lat)
!       call addfld('N4S_TM1' ,' ',' ',n4s_nm_out(:,i0:i1,lat),
!    |    i0,i1,k0,k1,lat)
! -----------------------------------------------------------------------
      end subroutine minor_n4s
!-----------------------------------------------------------------------
      subroutine comp_n2d(tn,o1,co,co2,n2,barm,
     |  nop,n2p,ne,n2d,
     |  lev0,lev1,lon0,lon1,lat)
!
!-----------------------------------------------------------------------
! Advance n2d (PCE):  This is called from driver at each subdomain latitude.
! Arguments for comp_n2d (14).
! ----------------------------------------------------------------------
! Last modified: swb (04/22/08)
!              : add new rates from chemrates
!              : revised day-night and low altitude settings
! ----------------------------------------------------------------------
!
      use qrj_module,only: qjn2, mk_xnmbari
      use fields_module,only: tlbc
      use cons_module,only: rmassinv_o1,rmassinv_co,rmassinv_n2,p0,
     |  expz,expzmid_inv,boltz,rmassinv_co2,rmassinv_n4s,rmass_o1,
     |  rmass_co,rmass_n2,rmass_co2,rmass_n4s,gask,grav
      use chemrates_module,only:rk1,rk3,rk4,rk5,r27,a15,a17
      use addfld_module,only: addfld
!
!     write(6,"('enter comp_n4s: lat=',i2)") lat
! ----------------------------------------------------------------------
! Input args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) ::
     |  tn,      ! neutral temperature (deg K)
     |  o1,      ! atomic oxygen (mmr)
     |  co,      ! carbon monoxide (mmr)
     |  co2,     ! carbon dioxide (mmr)
     |  n2,      ! molecular nitrogen (mmr)
     |  barm,    ! mean molecular weight
     |  nop,     ! NO+ ion  (#/cm3)
     |  n2p,     ! N2+ ion  (#/cm3)
     |  ne       ! electron density  (#/cm3)
!
! Output args (n2d:  mass mixing mixing):
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(out) ::
     |  n2d      ! atomic nitrogen excited (mmr)
!
! ----------------------------------------------------------------------
! Local:
      integer :: k,i  ! for local indices
      real,dimension(lev0:lev1,lon0:lon1) :: 
     |  tni,            ! TN at interfaces
     |  qjn2i,          ! JN2EFF at interfaces (qrj module)
     |  xnmbari,        ! xnmbar at interfaces (qrj module)
     |  xnmbar,         ! xnmbar at midlevels (qrj module)
     |  xo,             ! O1 number density (#/cm3)
     |  xco,            ! CO number density (#/cm3)
     |  xco2,           ! CO2 number density (#/cm3)
     |  xn2,            ! N2 number density (#/cm3)
     |  xn2d,           ! N2D number density (#/cm3)
     |  xn2dl,          ! log10 N2D number density (#/cm3)
     |  n2dprod,        ! N2D production (#/cm3/sec)
     |  taun2d,         ! N2D time scale (Earth Days)
     |  n2dloss         ! N2D loss frequency (#/sec)
!
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  n2loss1,            ! N(2D) loss term 1
     |  n2loss2,            ! N(2D) loss term 2
     |  n2loss3,            ! N(2D) loss term 3
     |  n2loss4,            ! N(2D) loss term 4
     |  n2loss5,            ! N(2D) loss term 5
     |  n2prod1,            ! N(2D) prod term 1
     |  n2prod2,            ! N(2D) prod term 2
     |  n2prod3,            ! N(2D) prod term 3
     |  n2prod4,            ! N(2D) prod term 4
     |  n2prod5,            ! N(2D) prod term 5
     |  x2nop,     ! NO+ ion  (#/cm3)
     |  x2n2p      ! N2+ ion  (#/cm3)
!
!
! ----------------------------------------------------------------------
! For addfld calls: Check what is needed!
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
!
! Calculate p0*e(-z)*barm/kT (1/(N*MBAR))at interfaces
!     subroutine mk_xnmbari(tni,barm,xnmbari,lev0,lev1,lon0,lon1,lat)
!
      call mk_xnmbari(tni,barm(:,lon0:lon1),
     |  xnmbari,
     |  lev0,lev1,lon0,lon1,lat)
! ----------------------------------------------------------------------
! Species Number densities at interfaces (#/cm3)
! Qrj dissociation rates (#/cm3.sec) at interfaces: JO2
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          xo(k,i)   = xnmbari(k,i)*o1(k,i)*rmassinv_o1
          xco(k,i)  = xnmbari(k,i)*co(k,i)*rmassinv_co
          xco2(k,i) = xnmbari(k,i)*co2(k,i)*rmassinv_co2
	  if (xco2(k,i) < 1.0) xco2(k,i) = 1.0
          xn2(k,i)  = xnmbari(k,i)*n2(k,i)*rmassinv_n2
          qjn2i(k,i) = 0.5*(qjn2(k,i,lat)+qjn2(k+1,i,lat))
!         x2nop(k,i) = nop(k,i)
!         x2n2p(k,i) = n2p(k,i)
!
!   LOG10 DENSITIES FOR PLOTTING 
!         xco2l(k,i)=log10(xco2(k,i))
!         xcol(k,i)=log10(xco(k,i))
!         xn2l(k,i)=log10(xo2(k,i))
!         xol(k,i)=log10(xo(k,i))
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!
! ----------------------------------------------------------------------
! (1)  Major Sources and Losses for n2d (NOX Chemistry Added)
! ---------------------------------------------------------------------
      do i=lon0,lon1
        do k=lev0,lev1-1
! ---------------------------------------------------------------------
!     SOURCES (#/cm3.sec) and SINKS(#/sec): N2D (both Venus day & night)
! ---------------------------------------------------------------------
! Combined, JN2-dissociation + PE dossociation scheme used 
! Combined jn2eff calculated in qrj.F
          n2dprod(k,i) =2.*ff*qjn2i(k,i)  +
!    |	      a17(k,i,lat)*nop(k,i)*ne(k,i)*gg + 
     |	      a17(k,i,lat)*nop(k,i)*ne(k,i)*gg*2.0 + 
     |	      r27(k,i,lat)*n2p(k,i)*xo(k,i)
          n2dloss(k,i) = (xco2(k,i)*rk1+ 
     |	      xco(k,i)*rk3 +  xo(k,i)*rk4 + xn2(k,i)*rk5)
!  Calculating the time scale (1/N2D loss frequency)
          taun2d(k,i) = 1./n2dloss(k,i)/86400.
!
!  PCE Solution for N2D (#/cm3): Above z = -8.0 only ( ~ 100 km)
          xn2d(k,i) = n2dprod(k,i)/n2dloss(k,i)
          if (dz == 0.5) then
             if (k < 17) xn2d(k,i) = 0.
          elseif (dz == 0.25) then
             if (k < 34) xn2d(k,i) = 0.
          endif
!  Populate n2d array (mass mixing ratio)
          n2d(k,i) = xn2d(k,i)*rmass_n4s/xnmbari(k,i)
!
          n2prod1(k,i) = 2.*ff*qjn2i(k,i)
!         n2prod2(k,i) = a17(k,i,lat)*nop(k,i)*ne(k,i)*gg
          n2prod2(k,i) = a17(k,i,lat)*nop(k,i)*ne(k,i)*gg*2.0
          n2prod3(k,i) = a15(k,i,lat)*n2p(k,i)*ne(k,i)*0.5
          n2prod4(k,i) = r27(k,i,lat)*n2p(k,i)*xo(k,i)
        enddo  ! k=lev0,lev1-1
!  PCE Plotting solution for N2D:
        do k=lev0,lev1-1
	  if (xn2d(k,i) > 1.) Then
               xn2dl(k,i) = log10(xn2d(k,i))
          else
               xn2dl(k,i) = 0.0
          endif
        enddo  ! k=lev0,lev1-1
! -----------------------------------------------------------------------
	 n2d(lev1,i) = n2d(lev1-1,i)
      enddo   !  i=lon0,lon1
! ----------------------------------------------------------------------
!  addflds calls
!   * N2D Production
!     call addfld('N2DPROD',' ',' ',n2dprod(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!   * N2D Loss Frequency
!     call addfld('N2DLOSS',' ',' ',n2dloss(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!   * N2D Time Scale
      call addfld('TAUN2D',' ',' ',taun2d(:,lon0:lon1),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!   * N2D Density (PCE) 
      call addfld('N2DPCE',' ',' ',xn2dl(:,lon0:lon1),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!   * N(2D) reaction rates (terms) 
      call addfld('N2PROD1',' ',' ',n2prod1(:,lon0:lon1),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('N2PROD2',' ',' ',n2prod2(:,lon0:lon1),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('N2PROD3',' ',' ',n2prod3(:,lon0:lon1),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('N2PROD4',' ',' ',n2prod4(:,lon0:lon1),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!   * Ion densities (PCE) 
!     call addfld('X2N2P',' ',' ',x2n2p(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('X2NOP',' ',' ',x2nop(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! ----------------------------------------------------------------------
      end subroutine comp_n2d
!-----------------------------------------------------------------------
      subroutine comp_no(tn,o1,co,co2,n2,barm,
     |  n4s,n2d,no,o2p,
     |  lev0,lev1,lon0,lon1,lat)
!
!-----------------------------------------------------------------------
! Advance no (minor):  after N2D and before N4S
! Arguments for comp_no (15).
! ----------------------------------------------------------------------
! Advance no by one time step. This is called from driver at
! each subdomain latitude.
! Last modified: swb (10/09/07)
!              : add new rates from chemrates
!              : add no_prod, no_loss terms
! ----------------------------------------------------------------------
!
      use qrj_module,only: qjn2, mk_xnmbari
      use fields_module,only: tlbc
      use cons_module,only: rmassinv_o1,rmassinv_co,rmassinv_n2,p0,
     |  expz,expzmid_inv,boltz,rmassinv_co2,rmassinv_n4s,rmass_o1,
     |  rmass_co,rmass_n2,rmass_co2,rmass_n4s,rmass_no,rmassinv_no,
     |  gask,grav
      use chemrates_module,only:rk1,rk10,rg1,rg2,rg3,k45,k46,k47
      use chapman_module,only: idn_lonlat
      use addfld_module,only: addfld
!
!     write(6,"('enter comp_no: lat=',i2)") lat
! ----------------------------------------------------------------------
! Input args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) ::
     |  tn,      ! neutral temperature (deg K)
     |  o1,      ! atomic oxygen (mmr)
     |  co,      ! carbon monoxide (mmr)
     |  co2,     ! carbon dioxide (mmr)
     |  n2,      ! molecular nitrogen (mmr)
     |  barm,    ! mean molecular weight
     |  n4s,     ! atomic nitrogen (mmr)
     |  n2d,     ! atomic nitrogen excited (mmr)
     |  no,      ! nitric oxide (mmr)
     |  o2p      ! O2+ ion (#/cm3)
!
! ----------------------------------------------------------------------
! Local:
      integer :: k,i,k0,k1,i0,i1,nk  ! for addfld (ask for check)
      real,dimension(lev0:lev1,lon0:lon1) :: 
     |  tni,            ! TN at interfaces
     |  xnmbari,        ! xnmbar at interfaces (qrj module)
     |  xnmbar,         ! xnmbar at midlevels (qrj module)
     |  xo,             ! O1 number density (#/cm3)
     |  xco,            ! CO number density (#/cm3)
     |  xco2,           ! CO2 number density (#/cm3)
     |  xn2,            ! N2 number density (#/cm3)
     |  xn4s,           ! N4S number density (#/cm3)
     |  xn2d,           ! N2D number density (#/cm3)
     |  xno,            ! NO number density (#/cm3)
     |  xnol,           ! NO number density (log10 #/cm3)
     |  xho2n,          ! HO2 number density (#/cm3)
     |  xo3n            ! O3 number density (#/cm3)
!
      real,dimension(lev0:lev1)::
     |  ho2n,           ! HO2 nightside (mmr)
     |  o3n             ! O3 nightside (mmr)
!
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  noloss1,            ! NO loss term 1
     |  noloss2,            ! NO loss term 2
     |  nolossT,            ! NO loss term Total
     |  tauno,              ! NO time scale (Earth Days)
     |  noprod1             ! NO prod term 1
!
! ----------------------------------------------------------------------
! For addfld calls: Check what is needed!
      k0=lev0
      k1=lev1
      i0=lon0
      i1=lon1
      nk = lev1-1
!
#include<defs.h>
! ----------------------------------------------------------------------
! Lower boundary:
! no_lbc(:,1)=A, n4s_lbc(:,2)=B, n4s_lbc(:,3)=C define lower boundary 
!   condition where A*DPSX/DZ + B*PSX + C = 0.
!
      do i=lon0,lon1
!
! Value at bottom given by specified mass missing ratio (start)
        no_lbc(i,1,lat) = 0.
        no_lbc(i,2,lat) = 1.
        no_lbc(i,3,lat) = -psno
!
! Upper boundary: Zero diffusive flux at top
        no_ubc(i,lat) = 0.
!
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
#if (NLEV==68)
! --------------------------------------------------------------------
! Low Altitude Chemistry extracted from Krasnopolsky (2010) for the nightside
! Nightside profiles in MMR (paper has number density in Krasn.m they were converted to mmr)
!  (05-19-10)
! Alt rage = ~70km - ~112km (zp = -16 to -5)
!********************************************************************
      ho2n(:) = (/
     |  0.,0.,0.,0.,0.,2.51E-15,7.68E-15,1.38E-14,2.90E-14,5.04E-13,
     |  9.45E-13,4.44E-13,2.01E-13,1.04E-13,4.43E-14,7.39E-15,0.,0.,
     |  0.,0.,0.,0.,0.,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  0.,0.,0.,0.,0.,0./)
      o3n(:) = (/
     |  0.,0.,0.,0.,0.,7.39E-11,1.19E-10,1.92E-10,3.02E-10,5.59E-10,
     |  1.07E-09,1.98E-09,3.38E-09,5.37E-09,8.88E-09,1.29E-08,
     |  2.16E-08,3.43E-08,5.45E-08,8.64E-08,1.35E-07,2.07E-07,
     |  3.87E-07,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     |  0.,0.,0.,0.,0.,0./)
!
#elif (NLEV==136)
!********************************************************************
! Double Res. data
!********************************************************************
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
#else
 write(6,"('ERROR: unsupported NLEV = ',i6)") nlev
#endif
! ----------------------------------------------------------------------
! Species Number densities at interfaces (#/cm3)
! Qrj dissociation rates (#/cm3.sec) at interfaces: JO2
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          xo(k,i)   = xnmbari(k,i)*o1(k,i)*rmassinv_o1
          xco(k,i)  = xnmbari(k,i)*co(k,i)*rmassinv_co
          xco2(k,i) = xnmbari(k,i)*co2(k,i)*rmassinv_co2
	  if (xco2(k,i) < 1.0) xco2(k,i) = 1.0
          xn2(k,i)  = xnmbari(k,i)*n2(k,i)*rmassinv_n2
          xn4s(k,i)  = xnmbari(k,i)*n4s(k,i)*rmassinv_n4s
          xn2d(k,i)  = xnmbari(k,i)*n2d(k,i)*rmassinv_n4s
          xno(k,i)  = xnmbari(k,i)*no(k,i)*rmassinv_no
!converting mmr to number density
          xo3n(k,i)   = xnmbari(k,i)*o3n(k)*(1./48.)
          xho2n(k,i)  = xnmbari(k,i)*ho2n(k)*(1./33.)

!
!   LOG10 DENSITIES FOR PLOTTING 
!         xco2l(k,i)=log10(xco2(k,i))
!         xcol(k,i)=log10(xco(k,i))
!         xn2l(k,i)=log10(xo2(k,i))
!         xol(k,i)=log10(xo(k,i))

        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!
! --------------------------------------------------------------------
! ----------------------------------------------------------------------
! (1)  Major Sources and Losses for no (NOX Chemistry Added)
! ---------------------------------------------------------------------
      do i=lon0,lon1
        do k=lev0,lev1-1
! ---------------------------------------------------------------------
!     SOURCES (#/cm3.sec) and SINKS(#/sec): for NO (#/cm3.sec 
! ---------------------------------------------------------------------
! Simplified (1-production; 2-losses)
!         no_prod(k,i,lat) = xn2d(k,i)*xco2(k,i)*rk1
! Simplified (3-production2; 2-losses)
          no_prod(k,i,lat) = xn2d(k,i)*xco2(k,i)*rk1 +
!    |                       xn4s(k,i)*xo(k,i)*rg1+
     |                       xn4s(k,i)*xo(k,i)*rg1(k,i,lat)+
     |                       xn4s(k,i)*xo(k,i)*xco2(k,i)*rg2(k,i,lat)
!    |                       xn4s(k,i)*xo(k,i)*xco2(k,i)*rg2
! Day/Night Split
         if (idn_lonlat(i,lat)==1)then
!!        no_loss(k,i,lat) = -(xn4s(k,i)*rg3(k,i,lat)+ o2p(k,i)*rk10)
!         no_loss(k,i,lat) = -(xn4s(k,i)*rg3 + o2p(k,i)*rk10)
          no_loss(k,i,lat) = -(xn4s(k,i)*rg3(k,i,lat)+ o2p(k,i)*rk10)
         else
          no_loss(k,i,lat) = -(xn4s(k,i)*rg3(k,i,lat)+ o2p(k,i)*rk10+
     |                        xo(k,i)*xco2(k,i)*k45(k,i,lat)+
     |                        xo3n(k,i)*k46(k,i,lat)+
     |                        xho2n(k,i)*k47(k,i,lat))
         endif

          noprod1(k,i) = no_prod(k,i,lat)
          noloss1(k,i) = xno(k,i)*xn4s(k,i)*rg3(k,i,lat)
!         noloss1(k,i) = xno(k,i)*xn4s(k,i)*rg3
          noloss2(k,i) = xno(k,i)*o2p(k,i)*rk10
          nolossT(k,i) = noloss1(k,i)+noloss2(k,i)
          tauno(k,i) = 1./nolossT(k,i)/86400.
	  xno(k,i) = - no_prod(k,i,lat)/no_loss(k,i,lat)  
	  if (xno(k,i) <= 1.0) then
	      xnol(k,i) = 0.0
	  else
	      xnol(k,i) = log10(xno(k,i))
	  endif
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
! ----------------------------------------------------------------------
!   * NO reaction rates (terms) 
      call addfld('NOPROD1',' ',' ',noprod1(:,lon0:lon1),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('NOLOSS1',' ',' ',noloss1(:,lon0:lon1),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('NOLOSS2',' ',' ',noloss2(:,lon0:lon1),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('TAUNO',' ',' ',tauno(:,lon0:lon1),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!   * NO Density (PCE): 
      call addfld('NOPCE',' ',' ',xnol(:,lon0:lon1),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
! ----------------------------------------------------------------------
      end subroutine comp_no
!-----------------------------------------------------------------------
      subroutine minor_no(tn,difk,o1,co,co2,n2,no,no_nm,no_out,
     |  no_nm_out,lev0,lev1,lon0,lon1,lat0,lat1)

!     Arguments (16): o1, co, co2, n2 needed by minor.F for barm, etc.
      use cons_module,only: rmass_no
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
     |  no,      ! nitric oxide (mmr)
     |  no_nm    ! nitric oxide at time n-1 
!
! -----------------------------------------------------------------------
! Output args:
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(out) ::
     |  no_out,   ! NO output 
     |  no_nm_out ! NO output at time n-1
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
! Minor returns n4s_out and n4s_nm_out. Module data n4s_prod,
! n4s_loss, etc, were defined by comp_n4s.  
!
!     call minor(tn,difk,o1,co,co2,n2,o2,o2_nm,o2_out,o2_nm_out,
!    |  o2_loss,o2_prod,o2_lbc,o2_ubc,rmass_o2,phi_o2,
!    |  alfa_o2,lev0,lev1,lon0,lon1,lat0,lat1,0)
!
      call minor(tn,difk,o1,co,co2,n2,no,no_nm,no_out,no_nm_out,
     |  no_loss,no_prod,no_lbc,no_ubc,rmass_no,phi_no,
     |  alfa_no,lev0,lev1,lon0,lon1,lat0,lat1,0)

!       call addfld('NO_OUT' ,' ',' ',n4s_out(:,i0:i1,lat),
!    |    i0,i1,k0,k1,lat)
!       call addfld('NO_TM1' ,' ',' ',n4s_nm_out(:,i0:i1,lat),
!    |    i0,i1,k0,k1,lat)
! -----------------------------------------------------------------------
      end subroutine minor_no
      end module n4s_module
!-----------------------------------------------------------------------
