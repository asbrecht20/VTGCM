!
      module chapman_module
      use params_module,only: nlonp4,nlat,dz,spval,nlevp1
      use addfld_module,only: addfld
      use init_module,only: istep
      implicit none
      real,dimension(nlonp4,nlat) :: 
     |  chi,         ! solar zenith angle
     |  sin_chi,     ! sin(chi)
     |  cos_chi,     ! cos(chi)
     |  rt_sinchi,   ! sqrt(sin(chi))
     |  slt          ! local time
!
! Day/night index is declared/set by chapman.F: Formerly declared in init.F
      integer :: idn(nlonp4)               ! day/night index  (for line integrals)
      integer :: idn_lonlat(nlonp4,nlat)   ! day/night index  (for eddyflds)
!
! sovmr,so2vmr,so3vmr are use-associated by qrj:
      real,dimension(nlevp1) :: sovmr,so2vmr,so3vmr  ! VMR from KINETICS
!
      contains
!-----------------------------------------------------------------------
!     subroutine chapman (31)
      subroutine chapman(z,tn,o1,co,co2,n2,o2,so,so2,barm,vo1,vco,vco2,
     | vn2,vo2,vso,vso2,vso3,sco1,scco,scco2,scn2,sco2,scso,scso2,scso3,
     | lev0,lev1,lon0,lon1,lat)
!     subroutine chapman (33)
!     subroutine chapman(z,tn,o1,co,co2,n2,o2,so,so2,so_nm,so2_nm,barm,
!    | vo1,vco,vco2,
!    | vn2,vo2,vso,vso2,vso3,sco1,scco,scco2,scn2,sco2,scso,scso2,scso3,
!    | lev0,lev1,lon0,lon1,lat)
      use init_module,only: secs,sin_sundec,cos_sundec,uthr
      use cons_module,only: pi,dphi,dlamda,re,rmass_o1,rmass_co,
     |  rmass_co2,rmass_n2,rmass_o2,rmass_so,rmass_so2,rmass_so3,
     |  secperhr,rtd
      use hist_module,only: modeltime ! for diag
      use fields_module,only: tlbc
!
! -----------------------------------------------------------------
! 9/18/06 btf: adapted for vtgcm from titan2. 
! Calculate line integrals for o1,co,n2,co2,O2:
! 9/3/13 swb: adapted for vtgcm from titan2. 
! Calculate line integrals for o1,co,n2,co2,O2,SO,SO2,SO3:
! 2/7/14 swb: Updated VMR profiles for SO, SO2, SO3 from 1-D KINETICS
! 1/20/16 swb: Updated VMR profile for SO3:  2x Zhang et al (2012) Model A
!                                            over 90-110 km (trial bump)
! 1/27/16 swb: Updated VMR profile for SO3:  back to only Zhang et al (2012) Model A
! -----------------------------------------------------------------
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
!
! 6/15/13 btf: These are inout so can assign global means to so,so2
! at first timestep. This can be changed back to intent(in) after
! so,so2 are available on the source history.
!
!     real,dimension(lev0:lev1,lon0-2:lon1+2),intent(inout) ::  
!    |  so,   ! sulfur monoxide (mmr) ! 6/14/13 btf: 4d prognostic
!    |  so2,  ! sulfur dioxide (mmr)  ! 6/14/13 btf: 4d prognostic
!     real,dimension(lev0:lev1,lon0-2:lon1+2),intent(inout) ::  
!    |  so_nm,   ! sulfur monoxide (mmr) ! 11/01/13 btf: 4d prognostic
!    |  so2_nm   ! sulfur dioxide (mmr)  ! 11/01/13 btf: 4d prognostic
!
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) :: 
     |  z ,   ! geopotential height (cm)
     |  tn,   ! neutral temperature (deg K)
     |  o1,   ! methane (mmr)
     |  co,   ! molecular hydrogen (mmr)
     |  co2,  ! co2 is calculated locally as 1.-o-co-n2
     |  n2,   ! molecular nitrogen (mmr)
     |  o2,   ! molecular oxygen (mmr)
!    |  so,   ! sulfur monoxide (mmr) ! 6/14/13 btf: not on source history yet
!    |  so2,  ! sulfur dioxide (mmr)  ! 6/14/13 btf: not on source history yet
     |  so,   ! sulfur monoxide (mmr) ! 10/21/13 swb: on source history
     |  so2,  ! sulfur dioxide (mmr)  ! 10/21/13 swb: on source history
     |  barm  ! mean molecular weight
!
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(out) :: 
     |  vo1,  ! o1 vertical integration
     |  vco,  ! co vertical integration
     |  vco2, ! co2 vertical integration
     |  vn2,  ! n2 vertical integration
     |  vo2,  ! o2 vertical integration
     |  vso,  ! so vertical integration
     |  vso2, ! so2 vertical integration
     |  vso3, ! so3 vertical integration
     |  sco1, ! o1 slant column integration
     |  scco, ! co slant column integration
     |  scco2,! co2 slant column integration
     |  scn2, ! n2 slant column integration
     |  sco2, ! o2 slant column integration
     |  scso, ! so slant column integration
     |  scso2,! so2 slant column integration
     |  scso3 ! so2 slant column integration
!
! Local:
      real :: rlat,coslat,sinlat,factor,rtpi,fmin,fmax
      integer :: k,i
      real :: rp(lev0:lev1,lon0:lon1), ! z+re
     |        ti(lev0:lev1,lon0:lon1)  ! tn at interfaces
      integer :: nlevs
      real,dimension(lev0:lev1,lon0:lon1) :: ! for diagnostic plotting
     |  chi_plt, slt_plt, idn_plt, chi_deg
      real,dimension(lev0:lev1,lon0-2:lon1+2) :: sob,so2b,so3b !SO,SO2,SO3 MMR (Test)
!
! --------------------------------------------------------------------------
!     Use so, so2, and so3 for testing of qrf_sox.F intermediate products 
!     before moving onto full so and so2 photochemical calculations.
!  Old
!     so3vmr(:) = (/
!    | 6.769E-09, 3.535E-10, 2.492E-11, 6.631E-12, 4.240E-12, 3.696E-12,
!    | 4.278E-12, 4.037E-12, 2.244E-12, 5.432E-13, 1.014E-12, 1.392E-11,
!    | 2.556E-10, 2.682E-09, 2.163E-08, 2.303E-08, 1.431E-08, 8.000E-09,
!    | 4.234E-09, 2.331E-09, 1.437E-09, 1.137E-09, 1.11E-09,
!    | 1.110E-09, 1.110E-09, 1.110E-09, 1.110E-09, 1.110E-09, 1.110E-09,
!    | 1.110E-09, 1.110E-09, 1.110E-09, 1.110E-09, 1.110E-09, 1.110E-09,
!    | 1.110E-09, 1.110E-09, 1.110E-09, 1.110E-09, 1.110E-09, 1.110E-09,
!    | 1.110E-09, 1.110E-09, 1.110E-09, 1.110E-09, 1.110E-09, 1.110E-09,
!    | 1.110E-09, 1.110E-09, 1.110E-09, 1.110E-09, 1.110E-09, 1.110E-09,
!    | 1.110E-09, 1.110E-09, 1.110E-09, 1.110E-09, 1.110E-09, 1.110E-09,
!    | 1.110E-09, 1.110E-09, 1.110E-09, 1.110E-09, 1.110E-09, 1.110E-09,
!    | 1.110E-09, 1.110E-09, 1.110E-09, 1.110E-09/)
! --------------------------------------------------------------------------
!  Updated: 2/7/14:  S. W. Bougher from KINETICS
! --------------------------------------------------------------------------
      sovmr(:) = (/ 
     |    3.03e-009,2.52e-009,2.36e-009,2.43e-009,
     |    2.90e-009,3.56e-009,3.69e-009,4.47e-009,
     |    1.02e-008,2.58e-008,4.57e-008,5.40e-008,
     |    7.98e-008,1.09e-007,1.32e-007,1.51e-007,
     |    1.77e-007,2.14e-007,2.59e-007,3.08e-007,
     |    3.50e-007,3.78e-007,3.78e-007,
     |    3.78e-007,3.78e-007,3.78e-007,3.78e-007,3.78e-007,3.78e-007,
     |    3.78e-007,3.78e-007,3.78e-007,3.78e-007,3.78e-007,3.78e-007,
     |    3.78e-007,3.78e-007,3.78e-007,3.78e-007,3.78e-007,3.78e-007,
     |    3.78e-007,3.78e-007,3.78e-007,3.78e-007,3.78e-007,3.78e-007,
     |    3.78e-007,3.78e-007,3.78e-007,3.78e-007,3.78e-007,3.78e-007,
     |    3.78e-007,3.78e-007,3.78e-007,3.78e-007,3.78e-007,3.78e-007,
     |    3.78e-007,3.78e-007,3.78e-007,3.78e-007,3.78e-007,3.78e-007,
     |    3.78e-007,3.78e-007,3.78e-007,3.78e-007/)
      so2vmr(:) = (/ 
     |    2.78e-007,1.21e-007,6.52e-008,4.18e-008,
     |    3.03e-008,3.11e-008,4.06e-008,5.51e-008,
     |    7.13e-008,8.43e-008,1.01e-007,1.40e-007,
     |    1.78e-007,2.27e-007,2.61e-007,2.70e-007,
     |    2.56e-007,2.23e-007,1.79e-007,1.30e-007,
     |    8.61e-008,5.71e-008,5.71e-008,
     |    5.71e-008,5.71e-008,5.71e-008,5.71e-008,5.71e-008,5.71e-008,
     |    5.71e-008,5.71e-008,5.71e-008,5.71e-008,5.71e-008,5.71e-008,
     |    5.71e-008,5.71e-008,5.71e-008,5.71e-008,5.71e-008,5.71e-008,
     |    5.71e-008,5.71e-008,5.71e-008,5.71e-008,5.71e-008,5.71e-008,
     |    5.71e-008,5.71e-008,5.71e-008,5.71e-008,5.71e-008,5.71e-008,
     |    5.71e-008,5.71e-008,5.71e-008,5.71e-008,5.71e-008,5.71e-008,
     |    5.71e-008,5.71e-008,5.71e-008,5.71e-008,5.71e-008,5.71e-008,
     |    5.71e-008,5.71e-008,5.71e-008,5.71e-008/)
!    RESTART #1: Petrass chemistry (standard KINETICS SO3VMR)
      so3vmr(:) = (/ 
     |    1.58e-009,7.27e-010,6.87e-010,3.91e-010,
     |    4.92e-011,2.00e-011,2.48e-011,3.12e-011,
     |    2.44e-011,7.19e-012,8.55e-012,5.88e-011,
     |    6.81e-010,1.03e-008,1.27e-008,9.23e-009,
     |    5.24e-009,2.64e-009,1.24e-009,5.98e-010,
     |    3.21e-010,2.12e-010,2.12e-010,
     |    2.12e-010,2.12e-010,2.12e-010,2.12e-010,2.12e-010,2.12e-010,
     |    2.12e-010,2.12e-010,2.12e-010,2.12e-010,2.12e-010,2.12e-010,
     |    2.12e-010,2.12e-010,2.12e-010,2.12e-010,2.12e-010,2.12e-010,
     |    2.12e-010,2.12e-010,2.12e-010,2.12e-010,2.12e-010,2.12e-010,
     |    2.12e-010,2.12e-010,2.12e-010,2.12e-010,2.12e-010,2.12e-010,
     |    2.12e-010,2.12e-010,2.12e-010,2.12e-010,2.12e-010,2.12e-010,
     |    2.12e-010,2.12e-010,2.12e-010,2.12e-010,2.12e-010,2.12e-010,
     |    2.12e-010,2.12e-010,2.12e-010,2.12e-010/)
!    RESTART #2: Petrass chemistry (2x enhanced KINETICS SO3VMR over 90-110 km)
!     so3vmr(:) = (/ 
!    |    1.58e-009,7.27e-010,6.87e-010,3.91e-010,
!    |    4.92e-011,2.00e-011,2.48e-011,3.12e-011,
!    |    2.44e-011,7.19e-012,1.71e-011,1.18e-010,
!    |    1.36e-009,2.06e-008,2.54e-008,1.85e-008,
!    |    1.05e-008,5.28e-009,2.48e-009,1.20e-009,
!    |    6.40e-010,4.24e-010,4.24e-010,
!    |    4.24e-010,4.24e-010,4.24e-010,2.24e-010,2.24e-010,2.24e-010,
!    |    4.24e-010,4.24e-010,4.24e-010,2.24e-010,2.24e-010,2.24e-010,
!    |    4.24e-010,4.24e-010,4.24e-010,2.24e-010,2.24e-010,2.24e-010,
!    |    4.24e-010,4.24e-010,4.24e-010,2.24e-010,2.24e-010,2.24e-010,
!    |    4.24e-010,4.24e-010,4.24e-010,2.24e-010,2.24e-010,2.24e-010,
!    |    4.24e-010,4.24e-010,4.24e-010,2.24e-010,2.24e-010,2.24e-010,
!    |    4.24e-010,4.24e-010,4.24e-010,2.24e-010,2.24e-010,2.24e-010,
!    |    4.24e-010,4.24e-010,4.24e-010,2.24e-010/)
!
! --------------------------------------------------------------------------
!
! --------------------------------------------------------------------------
!     Start: Overlay sovmr, so2vmr, so2vmr (VMR) upon so, so2 and so3 (MMR)
      sob = 0. ; so2b = 0.
      do i=lon0,lon1
        do k=lev0,lev1-1
          sob(k,i)  = sovmr(k)*(rmass_so/barm(k,i))
          so2b(k,i)  = so2vmr(k)*(rmass_so2/barm(k,i))
          so3b(k,i)  = so3vmr(k)*(rmass_so3/barm(k,i))
        enddo ! k=lev0,lev1-1
        sob(lev1,i) = sob(lev1-1,i)
        so2b(lev1,i) = so2b(lev1-1,i)
        so3b(lev1,i) = so3b(lev1-1,i)
      enddo ! lon0,1
!
! Use global mean sox to start 
! (until so,so2 are available on source history)
!
!     if (istep==1) then
!       so = 0. ; so2 = 0.
!       do i=lon0,lon1
!        do k=lev0,lev1
!           so(k,i) = sob(k,i)
!           so2(k,i) = so2b(k,i)
!         enddo
!       enddo
!
! Use global mean sox to start !
! (until so_nm,so2_nm are available on source history)
!
!     if (istep==1) then
!      so_nm = 0. ; so2_nm = 0.
!     do i=lon0,lon1
!       do k=lev0,lev1
!          so_nm(k,i) = sob(k,i)
!          so2_nm(k,i) = so2b(k,i)
!       enddo
!     enddo
!     endif 
!
!       write(6,"('chapman: istep=',i2,' lat=',i4,
!    |    ' define so,so2 with global means: so mnmx=',2e12.4,
!    |    ' so2 mnmx=',2e12.4)") istep,lat,
!    |    minval(so (:,lon0:lon1)),maxval(so (:,lon0:lon1)),
!    |    minval(so2(:,lon0:lon1)),maxval(so2(:,lon0:lon1))
!
!     endif
!
! Set up sun related geometry according to current model date and time.
! Set day/night index idn.
!
      nlevs = lev1-lev0+1 ! for addfsech
      rtpi = sqrt(pi)
      rlat = -.5*pi+(float(lat-1)+.5)*dphi
      coslat=cos(rlat)
      sinlat=sin(rlat)
      idn(lon0:lon1) = 1
!     Titan differences (09/17/04) ----------------------------------
!     -- slt and chi formulations same as Earth
!     -- slt calculated using 57600 Earth secs for 1-uthr for Titan 
!     -- chi limit is unique for Titan extended atmosphere  (off)
      do i=lon0,lon1
!       slt(i,lat)=amod(secs/3600.+(float(i-3)*dlamda+pi)*12./pi,24.)
!       slt(i,lat)=amod(secs/57600.+(float(i-3)*dlamda+pi)*12./pi,24.)
!
! 9/20/06 btf: Use secperhr (params.F) to calculate slt:
!       slt(i,lat)=amod(secs/secperhr+(float(i-3)*dlamda+pi)*
!    |    12./pi,24.)
        slt(i,lat)=amod((float(i-3)*dlamda+pi)*
     |    12./pi,24.)
        chi(i,lat)=acos(sin_sundec*sinlat+cos_sundec*coslat*
     |    cos(pi*(slt(i,lat)-12.)/12.))
!
        sin_chi(i,lat) = sin(chi(i,lat))
        cos_chi(i,lat) = cos(chi(i,lat))
        rt_sinchi(i,lat) = sqrt(sin_chi(i,lat))
!     ---- SZA>130 Extended Atmosphere of Titan
!       if (chi(i,lat) > 2.268919) idn(i) = 0
!     ---- SZA>105  Venus and Mars nightside criteria
!       if (chi(i,lat) > 1.8326)   idn(i) = 0
!     ---------------------------------------------------------------

        if (chi(i,lat) > 1.8326) idn(i) = 0

        chi_plt(:,i) = chi(i,lat)
        chi_deg(:,i) = chi(i,lat)*rtd
        slt_plt(:,i) = slt(i,lat)
        idn_plt(:,i) = idn(i)
        idn_lonlat(i,lat) = idn(i)
      enddo                

!     write(6,"('chapman: lat=',i3,' modeltime=',4i4,' uthr=',f5.2)") 
!    |  lat,modeltime,uthr
!     write(6,"('chi=',/,(10f7.2))") chi(:,lat)
!     write(6,"('slt=',/,(10f7.2))") slt(:,lat)
!     write(6,"('idn=',/,(10i7))") idn(:)
!
!     call addfld('CHI',' ','radians',chi_plt,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('CHI_DEG',' ','degrees',chi_deg,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('SLT',' ',' ',slt_plt,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('IDN',' ',' ',idn_plt,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
!
! rp = z+re
      do k=lev0,lev1
        rp(k,:) = z(k,lon0:lon1)+re
      enddo
!     call addfld('RP',' ',' ',rp,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! ti is tn at interfaces:
!     ti(1,:) = tn(lev1,lon0:lon1) ! tn bottom boundary is stored in top slot
!
      ti(lev0,:) = tlbc(lon0:lon1,lat) ! tn bottom boundary in tlbc (fields.F)
      do k=lev0+1,lev1-1
        ti(k,:) = .5*(tn(k-1,lon0:lon1)+tn(k,lon0:lon1))
      enddo
      ti(lev1,:) = tn(lev1-1,lon0:lon1) ! nlevp1 <- nlev
!     call addfld('TNI',' ',' ',ti,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! vtgcm co2 = 1.-o-co-n2:
!     co2 = spval
!     do k=lev0,lev1
!       do i=lon0,lon1
!         co2(k,i) = (1.-o1(k,i)-co(k,i)-n2(k,i))
!         if (co2(k,i) < 0.) co2(k,i) = 0.
!       enddo
!     enddo
!
!     call addfld('CHAP_O1',' ',' ' ,o1(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('CHAP_CO',' ',' ' ,co(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('CHAP_CO2',' ',' ',co2(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('CHAP_N2',' ',' ' ,n2(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('CHAP_O2',' ',' ' ,o2(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Calculate line integrals (vo[x] and sco[x] are output):
!
      call line_integ(o1,rmass_o1,vo1,sco1,ti,barm,
     |  cos_chi(lon0:lon1,lat),sin_chi(lon0:lon1,lat),
     |  rt_sinchi(lon0:lon1,lat),rp,idn(lon0:lon1),
     |  lon0,lon1,lev0,lev1,lat,1)
      call line_integ(co,rmass_co,vco,scco,ti,barm,
     |  cos_chi(lon0:lon1,lat),sin_chi(lon0:lon1,lat),
     |  rt_sinchi(lon0:lon1,lat),rp,idn(lon0:lon1),
     |  lon0,lon1,lev0,lev1,lat,2)
      call line_integ(co2,rmass_co2,vco2,scco2,ti,barm,
     |  cos_chi(lon0:lon1,lat),sin_chi(lon0:lon1,lat),
     |  rt_sinchi(lon0:lon1,lat),rp,idn(lon0:lon1),
     |  lon0,lon1,lev0,lev1,lat,2)
      call line_integ(n2,rmass_n2,vn2,scn2,ti,barm,
     |  cos_chi(lon0:lon1,lat),sin_chi(lon0:lon1,lat),
     |  rt_sinchi(lon0:lon1,lat),rp,idn(lon0:lon1),
     |  lon0,lon1,lev0,lev1,lat,2)
      call line_integ(o2,rmass_o2,vo2,sco2,ti,barm,
     |  cos_chi(lon0:lon1,lat),sin_chi(lon0:lon1,lat),
     |  rt_sinchi(lon0:lon1,lat),rp,idn(lon0:lon1),
     |  lon0,lon1,lev0,lev1,lat,2)
      call line_integ(so,rmass_so,vso,scso,ti,barm,
     |  cos_chi(lon0:lon1,lat),sin_chi(lon0:lon1,lat),
     |  rt_sinchi(lon0:lon1,lat),rp,idn(lon0:lon1),
     |  lon0,lon1,lev0,lev1,lat,2)
      call line_integ(so2,rmass_so2,vso2,scso2,ti,barm,
     |  cos_chi(lon0:lon1,lat),sin_chi(lon0:lon1,lat),
     |  rt_sinchi(lon0:lon1,lat),rp,idn(lon0:lon1),
     |  lon0,lon1,lev0,lev1,lat,2)
      call line_integ(so3b,rmass_so3,vso3,scso3,ti,barm,
     |  cos_chi(lon0:lon1,lat),sin_chi(lon0:lon1,lat),
     |  rt_sinchi(lon0:lon1,lat),rp,idn(lon0:lon1),
     |  lon0,lon1,lev0,lev1,lat,2)
!
! Save column number densities:
!     call addfld('VO1',' ',' ' ,vo1(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('VCO',' ',' ' ,vco(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('VCO2',' ',' ',vco2(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('VN2',' ',' ' ,vn2(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('VO2',' ',' ' ,vo2(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Save slant line integrals:
!     call addfld('SCO1',' ',' ',sco1(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('SCCO',' ',' ',scco(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('SCN2',' ',' ',scn2(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('SCCO2',' ',' ',scco2(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('SCO2',' ',' ',sco2(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!      call addfld('SCSO',' ',' ',scso(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('SCSO2',' ',' ',scso2(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('SCSO3',' ',' ',scso3(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
      end subroutine chapman
!-----------------------------------------------------------------------
      subroutine line_integ(f,fmass,v,s,ti,barm,cos_chi,sin_chi,
     |  rt_sinchi,rp,idn,lon0,lon1,lev0,lev1,lat,if)
      use params_module,only: dz
      use cons_module,only: expz,grav,p0,avo,gask,pi,re
      implicit none
!
! Args:
      integer,intent(in) :: lon0,lon1,lev0,lev1,lat,if
      real,intent(in) :: fmass
      real,dimension(lev0:lev1,lon0:lon1),intent(in) :: 
     |  ti,        ! tn at interfaces
     |  rp         ! z+re         (S9)
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) ::
     |  f,         ! input species density (mmr)
     |  barm       ! mean molecular mass
      real,dimension(lon0:lon1),intent(in) ::
     |  cos_chi,   ! cosine solar zenith angle
     |  sin_chi,   ! sine solar zenith angle
     |  rt_sinchi  ! sqrt(sin_chi)
      integer,intent(in) :: 
     |  idn(lon0:lon1) ! day-night index
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(out) ::
     |  v,         ! output column density
     |  s          ! output slant-line integral
!
! Local:
      real,parameter :: big=1.e80
      integer :: k,i,nlevs
      real :: factor,rtpi,exparg
      real :: 
     |  rtrp(lev0:lev1,lon0:lon1),      ! sqrt(rp/2hp) (S6)
     |  yp  (lev0:lev1,lon0:lon1),      ! yp and ip    (S5)
     |  r2ig(lev0:lev1,lon0:lon1)       ! 2.*ig        (S4)
      real,external :: expo
!
      nlevs = lev1-lev0+1 ! for addfsech
!
! Top:
      factor = avo*p0*expz(lev1-1)*exp(-.5*dz)/(fmass**2*grav)
      do i=lon0,lon1
        v(:,i) = f(:,i)
      enddo
      where(v < 0.) v = 0.
      do i=lon0,lon1
        v(lev1,i) = factor*.5*(v(lev1-1,i)+v(lev1,i))*
     |                    barm(lev1,i)
      enddo
!
! Integrate from top down:
      do i=lon0,lon1
        do k=lev1-1,lev0,-1
          factor = avo*p0*expz(k)/(fmass*grav)*dz
          v(k,i) = v(k+1,i)+factor*v(k,i)
        enddo
      enddo
!
! Set up for slant-line integral:
      factor = fmass*grav/(2.*gask)
      rtpi = sqrt(pi)
      do i=lon0,lon1
        do k=lev0,lev1
          rtrp(k,i) = sqrt(rp(k,i)*factor/ti(k,i)) ! S6
          yp(k,i) = rtrp(k,i)*abs(cos_chi(i))      ! S5
          if (yp(k,i) >= 8.) then
            yp(k,i) = v(k,i)*rtpi*rtrp(k,i)*
     |        (0.56498823/(0.06651874+yp(k,i)))
          else
            yp(k,i) = v(k,i)*rtpi*rtrp(k,i)*
     |        ((1.0606963+0.5564383*yp(k,i))/((yp(k,i)+1.7245609)*
     |          yp(k,i)+1.0619896))
          endif
        enddo ! k=lev0,lev1
      enddo ! i=lon0,lon1
      factor = grav*fmass/gask
      do i=lon0,lon1
        do k=lev0,lev1
!
! Avoid exceeding max arg to exp on non-unicos platform:
          exparg = rp(k,i)*(1.-sin_chi(i))*factor/ti(k,i)

          if (idn(i)==1.and.exparg < 650.) then ! daytime
            r2ig(k,i) = 2.*v(k,i)*exp(exparg)*
     |        rtpi*rt_sinchi(i)*rtrp(k,i)
          else
            r2ig(k,i) = big
          endif

!         r2ig(k,i) = 2.*v(k,i)*expo(exparg,0)*
!    |      rtpi*rt_sinchi(i)*rtrp(k,i)
        enddo
      enddo
!
! Slant line integrals (0 if obscured by earth):
! Slant line integrals (0 if obscured by titan):
      do i=lon0,lon1
        do k=lev0,lev1
          if (cos_chi(i) >= 0.) then
            s(k,i) = big
          else
            s(k,i) = r2ig(k,i)-yp(k,i)
          endif
          if (rp(k,i)*sin_chi(i)-re < 0.) s(k,i) = big
          if (cos_chi(i) >= 0.) s(k,i) = yp(k,i)
        enddo
      enddo
      end subroutine line_integ
      end module chapman_module
