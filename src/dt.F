      module dt_module
      use params_module,only: nlat,nlon,nlonp2,nlonp4,dz,spval,zibot
      implicit none
      private calc_tzm
      contains
! ----------------------------------------------------------------------
      subroutine dt(tn,tn_nm, un,vn,o1,co,barm,dift,cp,kt,km,hdt,
     |  w_upd,tn_upd,tn_nm_upd,lev0,lev1,lon0,lon1,lat0,lat1)
!
! Advance neutral temperature at current latitude: 
! ----------------------------------------------------------------------
! --09/18/06  btf : adapt Titan for Venus
! --10/03/06  Bougher  (ready for diagnostics)
! --10/03/06  Bougher  (cool_exp=cool_tot on RHS)
! --10/05/06  btf: new changes for cool in dt.F
! --10/05/06  Bougher  (merge btf changes for cool in dt.F)
! --10/05/06  Bougher  (add cool_exp and cool_imp interpolation using chi)
! --10/06/06  Bougher  (add cool_imp to qpart, as in tiegcm1.8)
!                      (Units checked against tiegcm1.8 dt.F)
! --10/26/06  Foster   subroutine cool (gather, scatter to all longitudes)
!             Bougher  (addfld calls for diagnostics on: 9-fields)
!             Bougher  (cool_exp correction for rhs)
!             Bougher  (interpolation test: chiday = 95., chinite=95.)
!                      (no interoplation across the terminator)
!             Bougher  (single term test: all cool_exp on and cool_imp off)
! --10/30/06  Foster   add qnir to qtotal (units erg/g/sec from Dickinson RT)
!             Bougher  (recast QTOT, QNET, QNIR)
! --10/31/06  Bougher  Important to start with 6,0,0 histories for proper
!                      aligmnent of thermal terms and CHI and SLT.
!                      Scale QNIR down by factor of 2.0 (units problem)?
!                      Sign change for adding QIR to rhs (OK)
!                      Correction to coolexp for day and nite (OK)
!                      10% Lincoeff added to lhs and rhs of equation (sign OK)
! --11/01/06  Bougher  Refinements to match vtgcm2 temperature structure (SMIN)
!                      Test 100% QNIR 
!                      Test E(EUV)= 20% (instead of 18%)
!                      Test lincoeff impacts on total 15-micron cooling
!                      Test cool_exp only impacts on total 15-micron cooling! 
!                      Array dimension or aligment problem somewhere in cool.F
! --11/02/06  Bougher  Refinements to match vtgcm2 temperature structure (SMIN)
!                      Fix to chideg dimension in cool.F
!                      Test 100% QNIR 
!                      Test E(EUV)= 20% (instead of 18%)
!                      TEST COOL_EXP ONLY (working with DT = 60 sec only) 
! --11/03/06  Bougher  Refinements to match vtgcm2 temperature structure (SMIN)
!                      Lower eddy thermal conduction impact in Ts (kzz reduced).
! --12/04/06  Bougher  Pass dift from eddyflds_module: dift (k,i,lat)
!                      Lower eddy thermal conduction impact in Ts (kzz reduced).
! --12/22/06  Bougher  Add O+O+CO2 recombination heating to Qtotal (erg/g/sec)
! --08/03/07  Bougher  Modify qnir for K/day (Roldan et al., 2000) to erg/gm/sec (internal) 
!                      Scaling both internally and for diagnostics applied. 
!                      Retain same full domain for qnir = qirgrid
! -- 09/25/07 Bougher  Modify use statements to obtain coolref quanities from
!                      coolref_init_module. Convert to local usage units and mid-levels
! -- 02/05/08 Bougher  Modify to compare ALAPSE = -grav/cp  and DT/DZ  (for convection)
! -- 02/06/08 Bougher  Modify DT/DZ below z = -13.0 to reflect observations and mitigate
!                      lower boundary problems at z = -16.0.  See VIRA: Seiff et al. (1986).
! -- 02/12/08 Bougher  Modify DT/DZ below z = -14.0 to reflect observations and mitigate
!                      lower boundary problems at z = -16.0.  See VIRA: Seiff et al. (1986).
!                      Move location of overlay to permit all smoothing and filtering.
! -- 02/22/16 Bougher  Set lower limit T(z) below 100 km (zp < -8.0) so that T-dependent
!                      CO-reactions can be tested with more realistic T(z) profile shape.
!                      After smoothing and filtering conducted.
! -- 02/26/16 Bougher  Pseudo-aerosol heating function for testing and warming T(z) over 80-95 km
!                      Asked C. Parkinson for paramater values: aeropk, szp0, #scale heights, zp-range
!                      After smoothing and filtering conducted.
! -- 03/03/16 Bougher  Pseudo-aerosol heating function for testing and warming T(z) over 80-95 km
!                      Asked C. Parkinson for paramater values: aeropk, szp0, #scale heights, zp-range
!                      Bottom side cutoff applied to reduce low altitude heating below ~85 km.
! -- 03/09-10/16 Bougher  Pseudo-aerosol heating function for testing and warming T(z) over 80-95 km
!                      Asked C. Parkinson for paramater values: aeropk, szp0, #scale heights, zp-range
!                      Top & Bottom side cutoff applied to reduce low altitude heating below ~80 km,
!                      and high altitude heating above ~90 km. Center around 87 km, narrower layer.
! -- 03/22/16 Bougher  Diagnostic heating/cooling terms (K/day) setup clearly for future use
!                      (quev2, qnir2, cond2, qir2, totdynm, qaero2)
! -- 03/30/16 Bougher  Diagnostic heating/cooling terms (K/day) setup clearly for future use
!                      (Add qnet2=  quev2 + qnir2 + cond2 + qir2 + totdynm + qaero2)
! ----------------------------------------------------------------------
!
      use input_module,only: step,oxvgcm_ncfile,oxvgcm_zm
      use init_module,only: iter
      use cons_module,only: freq_semidi,shapiro,dtx2inv,expz,
     |  rmassinv_n2,rmassinv_o1,rmassinv_co,tsurplus,p0,boltz,avo,grav,
     |  gask,expzmid,expzmid_inv,kut_5,kutt=>kut,dtsmooth,
     |  dtsmooth_div2, rmassinv_co2, pi, dphi, aeropk
      use chemrates_module,only: r2
      use oxvgcm,only: oxvgcm_nlev,oxvgcm_nudge,t_oxvgcm,tzm_oxvgcm
      use lbc,only: t_lbc
      use qrj_module,only: qtotal ! qtotal(nlevp1,lon0:lon1,lat0:lat1)
      use fields_module,only: tlbc,ulbc,vlbc,tlbc_nm
      use addfld_module,only: addfld
!
! ----------------------------------------------------------------------
! q is near IR heating (erg/g/sec), read from DATAVOL by read_ncdata.F. 
! Call it qnir here, to avoid conflict with the local q coeffs 
!   array in this routine. Shut this down permanently!
!
!     use read_ncdata_module,only: qnir=>q  ! full domain (nlon,nlat,nlevp1)
! ----------------------------------------------------------------------
! qirgrid is near IR heating (K/day), read from qirgrid_init 
! Call it qnir here, and stuff qirgrid into qnir (like before) 
!
      use qirgrid_init_module,only: qnir=>qirgrid ! full domain (nlonp4,nlat,nlevp1)
! ----------------------------------------------------------------------
!
! Brings in the Suprathermal Electron Heating
!     use eheat_module,only: fheatfun ! Suprathermal Electrons 
! ----------------------------------------------------------------------
!
#ifdef MPI
      use mpi_module,only: mp_bndlons_f3d, mp_periodic_f3d
#endif
      implicit none
!
! Args: 15 plus 6-dimension parameters
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
!
! Input:
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),intent(in)::
     |  tn_nm, ! neutral temperature, time n-1 
     |  un,    ! neutral zonal velocity (cm/sec)
     |  vn,    ! neutral zonal velocity (cm/sec)
     |  o1,    ! atomic oxygen (mmr)
     |  co,    ! carbon monoxide (mmr)
     |  barm,  ! mean molecular weight
     |  dift,  ! eddy thermal coefficient (1/sec)  (sub eddyflds)
     |  cp,    ! specific heat (ergs/deg/gm)           (sub cpktkm)
     |  kt,    ! molecular diffusion (ergs/cm/deg/sec) (sub cpktkm)
     |  km,    ! molecular viscosity (gm/cm/sec)       (sub cpktkm)
     |  hdt,   ! horizontal diffusion of tn (from sub hdif3, hdif.F)
     |  w_upd  ! updated vertical velocity (swdot.F)
!
! tn is intent(inout) so it can be nudged at the beginning:
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(inout) :: tn     ! neutral temperature (deg K)
!
! Output:
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(out) ::
     |  tn_upd,   ! updated tn (output)
     |  tn_nm_upd ! updated tn at time n-1 (output)

! btf: Temporary for vtgcm: when qrj is ready, qtota can be use-associated 
! from qrj module.
!     real,dimension(lev0:lev1,lon0:lon1,lat0:lat1),intent(in) ::
!    |  qtotal ! will be use-associated from qrj
!
! Local: Add xlapse and dtdzp (2/06/08)
      integer :: kk,k,i,lonbeg,lonend,lat,j
      integer :: nk,nkm1,nlevs
      real :: tnlbc(lon0:lon1,lat0:lat1) ! lower boundary condition
      real :: dtdzp, xlapse
      real :: tzm(lev0:lev1,lat0:lat1)   ! zonal mean tn
      real :: tzm_pert(lon0:lon1)        ! perturations about mean
!
! diagnostics for zonal mean nudging:
      real :: tzm_vtgcm_diag(lev0:lev1,lon0:lon1,lat0:lat1)
      real :: tzm_oxvgcm_diag(lev0:lev1,lon0:lon1,lat0:lat1)
!
! Local at 2d:
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  cptn,       ! cp*exp(-s)*(V.del(T(n))-1/(2*DT))*T(n-1) (k+1/2)
     |  mbar,       ! mean molecular weight
     |  qm,         ! heating due to molecular diffusion
     |  total_heat, ! total heating  (erg/gm/sec)!
     |  dudz,       ! du/dz(k)
     |  dvdz,       ! du/dz(k)
     |  g,          ! g*KT/(p0*H*Ds**2)
     |  f,          ! g*eps/(p0*2*H*Ds)
     |  h,          ! scale height R*T/(M*g) (cm)
     |  hsc,        ! scale height R*T/(M*g) (cm)
     |  rho,        ! density
     |  tni,        ! tn at interfaces
     |  p, q, r,    ! coefficients for tridiagonal solver
     |  rhs,        ! right hand side of trsolv
     |  qpart,      ! part of q coeff (LHS)
     |  eheating,   ! The local variable for Suprathermal Electron Heating 
     |  pcp,        ! rho*cp (Converts K/s to Ergs/cm^3/s)  
     |  tkm,        ! for conduction 
     |  tkp,        ! for conduction 
     |  tk,         ! for conduction 
     |  pkn,        ! for conduction 
     |  qkn,        ! for conduction 
     |  rkn,        ! for conduction 
     |  qnir_diag,  ! for qnir to secondary histories
     |  tnlbc_diag  ! tnlbc_diag is redundant in vertical
!
! Local at 1d:
      real :: szp(lev0:lev1), glat
!
! Cooling (see sub cool):
      real,dimension(lev0:lev1,lon0:lon1,lat0:lat1) ::
     |  cool_exp,   ! explicit cooling (see sub cool) K/sec
     |  cool_imp    ! implicit cooling (see sub cool) K/sec
!
!-----------------------------------------------------
! Diagnostic 3D variables:  JMB 06-14-05
! Diagnostic 3D variables:  SWB 10-02-06
! Diagnostic 3D variables:  SWB 12-22-06
! Diagnostics are output in erg/g/sec
      real,dimension(lev0:lev1,lon0:lon1,lat0:lat1) ::
     |  xhydroh, ! Horizontal Advection Term (K/s)
     |  xhydrov, ! Vertical Advection Term (K/2)
     |  xhydrot, ! Total Advection Term (K/s)
     |  mden,    ! Mass Density (g/cm3) 
     |  adia,    ! Adiabatic Cooling/Heating 
     |  qeuv,    ! EUV, UV Heating  (erg/g/sec) 
     |  qeuvez,  ! EUV, UV Heating 
     |  qir,     ! IR Cooling  (erg/g/sec) 
     |  qirez,   ! IR Cooling 
     |  qtot,    ! qeuv + qnir (total radiative heating: erg/g/sec)
     |  qaeropk, ! peak aerosol/UV heating  (erg/gm/sec)!
     |  qaero,   ! layer aerosol heating (zp,lat)  (erg/gm/sec)!
     |  qt,      ! qtot - qir  (net heating: erg/g/sec)
     |  qrec,    ! Heating due to O+O+CO2 recombination (erg/g/sec)
     |  qmol,    ! Heating due to Molecular Diffusion 
     |  qel,     ! Heating due to Suprathermal Electrons 
     |  cond,    ! Molecular Conduction (erg/g/sec)
     |  dtdz,    ! dTdz (K/km units)
     |  alapse   ! adiabatic lapse (K/km units)
! Diagnostics are output in K/day
      real,dimension(lev0:lev1,lon0:lon1,lat0:lat1) ::
     |  qnet2,    ! Net total heating (K/day)
     |  qeuv2,    ! EUV, UV Heating  (K/day) 
     |  qnir2,    ! Near IR Heating  (K/day) 
     |  qir2,     ! IR Cooling  (K/day) 
     |  qaero2,   ! Layer aerosol heating (zp,lat)  (K/day)
     |  adia2,    ! adiabatic heating/cooling (K/day)
     |  xhydrot2, ! Total Advection Term (K/day)
     |  totdynm,  ! Total dynamical heating/cooling (adia2 + xhydrot2) (K/Day)
     |  cond2     ! Molecular Conduction (K/day)
!
! Local at 3d (tnsmooth needs lat dimension only for sub smooth):
      real,dimension(lev0:lev1,lon0:lon1,lat0:lat1) ::
     |  tnsmooth,   ! zonal and meridional smoothing of tn_nm
     |  advec_tn    ! horizontal advection (output of sub advec)
!
      nk = lev1-lev0+1
      nkm1 = nk-1
      nlevs = nk

!     rstep = float(step)
!     expt = cexp(ci*freq_semidi*rstep*iter)
!     expt2 = cexp(ci*.5*freq_semidi*rstep*iter)
!     expta = 1.

!     xlapse=3.8
      xlapse=5.0
!
! Calculate zonal means of tn for oxvgcm nudging:
      if (len_trim(oxvgcm_ncfile) > 0.and.oxvgcm_zm > 0) then
        call calc_tzm(tn(:,lon0:lon1,lat0:lat1),tzm,
     |    lev0,lev1,lon0,lon1,lat0,lat1)
!
! Plot zonal means (diags are redundant in longitude):
        do j=lat0,lat1
          do i=lon0,lon1
            tzm_vtgcm_diag (:,i,j) = tzm(:,j)        ! redundant in lon
            tzm_oxvgcm_diag(:,i,j) = tzm_oxvgcm(j,:) ! redundant in lon
          enddo
          call addfld('tzm_vtgcm',' ',' ',tzm_oxvgcm_diag(:,:,j),
     |      'lev',lev0,lev1,'lon',lon0,lon1,j)
          call addfld('tzm_oxvgcm',' ',' ',tzm_oxvgcm_diag(:,:,j),
     |      'lev',lev0,lev1,'lon',lon0,lon1,j)
        enddo
      endif
!
! First latitude scan for dt:
      do lat=lat0,lat1
!
! Nudge with Oxford VGCM (t_oxvgcm(nlonp4,nlat,oxvgcm_nlev)). See lbc.F.
! oxvgcm_nudge==1.0 at k==1, decreasing upward to 0 at k==oxvgcm_nlev.
! When nudge is 1, then we use 100% oxvgcm, and when nudge is 0, then
!   we use 100% vtgcm. Bottom boundary is already 100% oxvgcm, so
!   we start at k==2.
! Nudge with zonal means of oxvgcm if namelist parameter oxvgcm_zm > 0
! (zonal means of tn were calculated in call to calc_tzm above)
!
        if (len_trim(oxvgcm_ncfile) > 0) then

          call addfld('t_prenudge',' ',' ',tn(:,lon0:lon1,lat),
     |      'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Nudge zonal means oxvgcm (zm were calculated above, before lat loop):
          if (oxvgcm_zm > 0) then
            do k=2,oxvgcm_nlev
!
! Save perturbations about the mean:
              do i=lon0,lon1
                tzm_pert(i) = tn(k,i,lat)-tzm(k,lat) ! save perturbation about mean
              enddo
!
! Nudge vtgcm zonal mean tzm with oxvgcm zonal mean tzm_oxvgcm:
              tzm(k,lat) = tzm(k,lat)*(1.-oxvgcm_nudge(k))+
     |          tzm_oxvgcm(lat,k)*oxvgcm_nudge(k)
!
! Update tn w/ nudged zm and restore perturbations:
              do i=lon0,lon1
                tn(k,i,lat) = tzm(k,lat)+tzm_pert(i)
              enddo
            enddo ! k=2,oxvgcm_nlev

            call addfld('t_zmnudge',' ',' ',tn(:,lon0:lon1,lat),
     |        'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Normal (non-zm) nudging with t_oxvgcm:
          else 
            do k=2,oxvgcm_nlev
              do i=lon0,lon1
                tn(k,i,lat) = tn(k,i,lat)*(1.-oxvgcm_nudge(k))+
     |            t_oxvgcm(i,lat,k)*oxvgcm_nudge(k)
              enddo
            enddo
!           write(6,"('dt nudging tn with oxvgcm: lat=',i3,
!    |        ' oxvgcm_nudge=',/,(6e12.4))") lat,oxvgcm_nudge
          endif ! zm or normal nudging
        endif ! oxvgcm_ncfile
!
! Lower boundary t_lbc was calculated by sub tuvz_lbc (lbc.F)
! t_lbc is dimensioned globally, but only defined at subdomains,
! tnlbc is local.
!
        tnlbc(:,lat) = t_lbc(lon0:lon1,lat)
        do k=lev0,lev1
          tnlbc_diag(k,:) = tnlbc(:,lat) ! redundant in lev
        enddo
        call addfld('TN_LBC',' ',' ',tnlbc_diag,
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! -------------------------------------------------------------------
! Horizontal advection (pass k vs i slices at full task subdomain 
! longitudes, and the 5 latitudes centered over the current latitude).
! Vertical advection. Sub advecv adds vertical advection to advec_tn.
! -------------------------------------------------------------------
! Horizontal advection. 
        call advec(tn(:,:,lat-2:lat+2),advec_tn(:,:,lat),
     |    lev0,lev1,lon0,lon1,lat)
!
! Set Diagnostic Variable xhydroh (K/S) 
        xhydroh(:,:,lat) = advec_tn(:,:,lat)
!
! Vertical advection. Sub advecv adds vertical advection to advec_tn.
!
        call advecv(tn(:,:,lat),tnlbc(:,lat),advec_tn(:,:,lat),
     |    lev0,lev1,lon0,lon1,lat)
! Set Diagnostic Variable xhydrov, xhydrot (K/S) 
!
        xhydrov(:,:,lat) = advec_tn(:,:,lat) - xhydroh(:,:,lat)
        xhydrot(:,:,lat) = xhydrov(:,:,lat) + xhydroh(:,:,lat)
!       xhydrot(:,:,lat) = advec_tn(:,:,lat) 
!
!--------------------------------------------------------
!       call addfld('ADVEC_TN',' ',' ',advec_tn(:,:,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
        call addfld('CP'   ,' ',' ',cp(:,lon0:lon1,lat),
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! End first latitude scan:
      enddo ! lat=lat0,lat1
!
! Shapiro smoother for tn at time n-1:
      call smooth(tn_nm,tnsmooth,lev0,lev1,lon0,lon1,lat0,lat1,0)
!
! --------------------------------------------------------------
! Explicit and implicit cooling (sub cool must do gather/scatter 
!   in longitude dimension)
! Calculate cool_exp(lev0:lev1,lon0:lon1,lat:lat1):
! Calculate cool_imp(lev0:lev1,lon0:lon1,lat:lat1):
!     cool_exp = 0. ! 3d whole-array init
!     cool_imp = 0.
      call cool(
     |  tn  (:,lon0:lon1,lat0:lat1),
     |  o1  (:,lon0:lon1,lat0:lat1),
     |  barm(:,lon0:lon1,lat0:lat1),
     |  cool_exp,cool_imp,lev0,lev1,lon0,lon1,lat0,lat1)
!
      do lat=lat0,lat1
        call addfld('COOL_EXP',' ',' ',cool_exp(:,:,lat),
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
        call addfld('COOL_IMP',' ',' ',cool_imp(:,:,lat),
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
      enddo 
! --------------------------------------------------------------
!
! Begin second latitude scan:
      do lat=lat0,lat1
!
! Set cptn and mbar (k+1/2):
! (Earlier versions apparently assumed zero periodic points for
!  tnsmooth, since they were not set in smoothing. See sub smooth,
!  where the periodic points are set to zero to avoid NaNS fpe
!  in the following loop)
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          cptn(k,i) = .5*(cp(k,i,lat)+cp(k+1,i,lat))*expz(k)*
     |      (advec_tn(k,i,lat)-dtx2inv*tnsmooth(k,i,lat))
          mbar(k,i) = 1./(o1(k+1,i,lat)*rmassinv_o1 + 
     |      co(k+1,i,lat)*rmassinv_co+(1.-o1(k+1,i,lat)-co(k+1,i,lat))*
     |      rmassinv_co2)
!    |      rmassinv_n2)
        enddo ! k=lev0,lev1-1
        cptn(lev1,i) = spval
        mbar(lev1,i) = spval
      enddo ! i=lon0,lon1
!     call addfld('TNSMOOTH'   ,' ',' ',tnsmooth(:,:,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('CPTN1',' ',' ',cptn,
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('MBAR' ,' ',' ',mbar,
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('QTOTAL',' ',' ',qtotal(:,lon0:lon1,lat),
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! ---------------------------------------------------------------------
! Save diag qnir:
! -- Old Dickinson qnir (erg/g/sec). No conversion.
! -- New Roldan et al qnir (K/day). Conversion applied yielding erg/g/sec.
! -- COMPARE QNIR (Roldan et al., 2000) with vtgcm2 benchmark simulations
! ---------------------------------------------------------------------
      do i=lon0,lon1
        do k=lev0,lev1-1
!       qnir_diag(k,i) = 0.5*(qnir(i,lat,k)+ qnir(i,lat,k+1))
        qnir_diag(k,i) = 0.5*(qnir(i,lat,k)+ qnir(i,lat,k+1))*
     |     0.5*(cp(k,i,lat)+cp(k+1,i,lat))/86400.
        enddo ! k=lev0,lev1-1
	qnir_diag(lev1,i) = qnir_diag(lev1-1,i)
      enddo ! i=lon0,lon1
!
      call addfld('QNIR' ,' ',' ',qnir_diag,
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! ---------------------------------------------------------------------
! Total heat sources are in total_heat (s5).
! ---------------------------------------------------------------------
      do i=lon0,lon1
        do k=lev0,lev1-1
!
!--------------------------------------------------------
! Add radiative Solar EUV heating from qrj:  (erg/gm/sec)
! Set Diagnostic Variable QEUV : EUV Radiative Heating (Ergs/g/s)
! Set Diagnostic Variable QEUV2 : EUV Radiative Heating (K/day)
          total_heat(k,i) = .5*(qtotal(k,i,lat)+qtotal(k+1,i,lat))
          qeuv(k,i,lat) = total_heat(k,i)
          qeuv2(k,i,lat) = qeuv(k,i,lat)*86400./
     |                     (0.5*(cp(k,i,lat)+cp(k+1,i,lat)))
          qeuvez(k,i,lat) = qeuv(k,i,lat)*expz(k)
!--------------------------------------------------------
! Add near-IR heating from datavol (erg/g/sec): No scaling 
!         total_heat(k,i) = total_heat(k,i)+
!    |      0.5*(qnir(i,lat,k)+qnir(i,lat,k+1))
! Add near-IR heating from qirgrid_init: With scaling 
! to convert from K/day to erg/gm/sec
          total_heat(k,i) = total_heat(k,i)+
     |     0.5*(qnir(i,lat,k)+qnir(i,lat,k+1))*
     |     0.5*(cp(k,i,lat)+cp(k+1,i,lat))/86400.
          qnir2(k,i,lat) = 0.5*(qnir(i,lat,k)+qnir(i,lat,k+1))
          qtot(k,i,lat) = total_heat(k,i)
!--------------------------------------------------------
!
! Add heating from 4th order horizontal diffusion (hdt from sub hdif3):
          total_heat(k,i) = total_heat(k,i)+hdt(k,i,lat)
!--------------------------------------------------------
!
! Add heating due to atomic oxygen recombination: Venus (12/22/06) 
          qrec(k,i,lat) = tsurplus*r2(k,i,lat)*
     |      (p0*expz(k)*mbar(k,i)/(boltz*tn(k,i,lat))*o1(k,i,lat)*
     |      rmassinv_o1)**2*avo/mbar(k,i)
          total_heat(k,i) = total_heat(k,i)+qrec(k,i,lat)
!--------------------------------------------------------
!
! Add Parkinson Pseudo-aerosol/UV heating  : Venus (2/26/2016) 
! Updated Parkinson Pseudo-aerosol/UV heating  : Venus (3/23/2016) 
! aeropk (K/day) to qaeropk (erg/(gm.sec)
          szp(k) = zibot + (k-1)*0.5
          GLAT= -pi*0.5 + (FLOAT(lat-1)+0.5)*dphi
          qaeropk(k,i,lat) = aeropk*0.5*(cp(k,i,lat)+cp(k+1,i,lat))
     |                       /86400.
!    **case A**
!         qaero(k,i,lat) = qaeropk(k,i,lat)*exp(-abs(szp(k)+12.0)/1.75)
!    |                     *cos(GLAT)
!    **case B**
!         qaero(k,i,lat) = qaeropk(k,i,lat)*exp(-abs(szp(k)+12.5)/2.0)
!    |                     *cos(GLAT)
!    **case C**
!         qaero(k,i,lat) = qaeropk(k,i,lat)*exp(-abs(szp(k)+11.0)/2.5)
!    |                     *cos(GLAT)
!    **case D**
!         qaero(k,i,lat) = qaeropk(k,i,lat)*exp(-abs(szp(k)+11.0)/2.0)
!    |                     *cos(GLAT)
!         qaero(k,i,lat) = qaeropk(k,i,lat)*exp(-abs(szp(k)+11.0)/1.75)
!    |                     *cos(GLAT)
!    **case E**
!         qaero(k,i,lat) = qaeropk(k,i,lat)*exp(-abs(szp(k)+11.0)/1.75)
!    |                     *cos(GLAT)
!         if (szp(k) < -13.0) qaero(k,i,lat) =  qaeropk(k,i,lat)*
!    |                        exp(-abs(szp(k)+11.0)/0.8)*cos(GLAT)
!    **case F1 and F2**
          qaero(k,i,lat) = qaeropk(k,i,lat)*exp(-abs(szp(k)+12.0)/1.25)
     |                     *cos(GLAT)
          if (szp(k) < -13.0) qaero(k,i,lat) =  qaeropk(k,i,lat)*
     |                        exp(-abs(szp(k)+12.0)/0.8)*cos(GLAT)
!         qaero(k,i,lat) = qaeropk(k,i,lat)*exp(-abs(szp(k)+12.5)/1.25)
!    |                     *cos(GLAT)
!         if (szp(k) < -13.0) qaero(k,i,lat) =  qaeropk(k,i,lat)*
!    |                        exp(-abs(szp(k)+12.5)/0.8)*cos(GLAT)
!--------------------------------------------------------
          total_heat(k,i) = total_heat(k,i)+qaero(k,i,lat)
          qaero2(k,i,lat) = qaero(k,i,lat)*86400./
     |                      (0.5*(cp(k,i,lat)+cp(k+1,i,lat)))
!--------------------------------------------------------
!
        enddo ! k=lev0,lev1-1
        total_heat(lev1,i)= total_heat(lev1-1,i)
      enddo ! i=lon0,lon1

      call addfld('HEATING1',' ',' ',total_heat,
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!--------------------------------------------------------

      do i=lon0,lon1
        do k=lev0+1,lev1-1
          tni(k,i) = .5*(tn(k-1,i,lat)+tn(k,i,lat))
          h(k,i) = gask*tni(k,i)/barm(k,i,lat)
          rho(k,i) = p0*expzmid_inv*expz(k)/h(k,i)
          h(k,i) = h(k,i)/grav
        enddo ! k=lev0+1,lev1-1

        tni(lev0,i) = tlbc(i,lat)        ! bottom boundary in tlbc
        tni(lev1,i) = tn(lev1-1,i,lat)
        h(lev0,i) = gask*tni(lev0,i)/barm(lev0,i,lat)
        h(lev1,i) = gask*tni(lev1,i)/barm(lev1,i,lat)
        rho(lev0,i) = p0*expzmid_inv*expz(lev0)/h(lev0,i)
        rho(lev1,i) = p0*expzmid*expz(lev1-1)/h(lev1,i)
        h(lev0,i) = h(lev0,i)/grav
        h(lev1,i) = h(lev1,i)/grav
      enddo ! i = lon0, lon1

!     do i=lon0,lon1
!       do k=lev0+1,lev1-1
!         tni(k,i) = .5*(tn(k-1,i,lat)+tn(k,i,lat))
!         hsc(k,i) = gask*tni(k,i)/barm(k,i,lat)
!         rho(k,i) = p0*expzmid_inv*expz(k)/hsc(k,i)
!         hsc(k,i) = hsc(k,i)/grav
!       enddo ! k=lev0+1,lev1-1
!
!       tni(lev0,i) = tlbc(i,lat)        ! bottom boundary in tlbc
!       tni(lev1,i) = tn(lev1-1,i,lat)
!       hsc(lev0,i) = gask*tni(lev0,i)/barm(lev0,i,lat)
!       hsc(lev1,i) = gask*tni(lev1,i)/barm(lev1,i,lat)
!       rho(lev0,i) = p0*expzmid_inv*expz(lev0)/hsc(lev0,i)
!       rho(lev1,i) = p0*expzmid*expz(lev1-1)/hsc(lev1,i)
!       hsc(lev0,i) = hsc(lev0,i)/grav
!       hsc(lev1,i) = hsc(lev1,i)/grav
!     enddo ! i = lon0, lon1

!     call addfld('TNI',' ',' ',tni,'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('HSC',' ',' ',hsc,'lev',lev0,lev1,'lon',lon0,lon1,lat)
      call addfld('HSC',' ',' ',h,'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('DEN',' ',' ',rho,'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Set Diagnostic Variable mden (g/cm^3) 
      mden(:,:,lat) = rho(:,:)
!
! Add heating due to molecular diffusion:
! du/dz and dv/dz (s10, s11):
      do i=lon0,lon1
        do k=lev0+1,lev1-2
          dudz(k,i) = (un(k+1,i,lat)-un(k-1,i,lat))/(2.*dz)  ! s10
          dvdz(k,i) = (vn(k+1,i,lat)-vn(k-1,i,lat))/(2.*dz)  ! s11
        enddo ! k=lev0+1,lev1-2
!
! Lower boundary:
        dudz(1,i) = (un(1,i,lat)+1./3.*un(2,i,lat)-4./3.*
     |    ulbc(i,lat))/dz
        dvdz(1,i) = (vn(1,i,lat)+1./3.*vn(2,i,lat)-4./3.*
     |    vlbc(i,lat))/dz
!
! Upper boundary:
        dudz(lev1-1,i) = dudz(lev1-2,i)/3.
        dvdz(lev1-1,i) = dvdz(lev1-2,i)/3.
!
! qm = heating due to molecular diffusion:
! (km = molecular viscosity from sub cpktkm)
        do k=lev0,lev1-1
          qm(k,i) = grav**2*mbar(k,i)*.5*(km(k,i,lat)+km(k+1,i,lat))/
     |      (p0*gask*expz(k)*tn(k,i,lat))*(dudz(k,i)**2+dvdz(k,i)**2)
!--------------------------------------------------------
! Set Diagnostic Variable qmol (Ergs/g/s) 
!
          qmol(k,i,lat) = qm(k,i)
!--------------------------------------------------------
!
! Add qm to total heating:
          total_heat(k,i) = total_heat(k,i)+qm(k,i)
!
! Now Eheating can be calculated:----
!
! This next part puts eheating in Ergs/g/s
!         eheating = fheatfun(k,i,lat)/rho 
          eheating(k,i) = 0.
!--------------------------------------------------------
! Set Diagnostic Variable qel (Ergs/g/s) 
!
          qel(k,i,lat) = eheating(k,i)
!--------------------------------------------------------
!
          total_heat(k,i) = eheating(k,i) + total_heat(k,i)
!
! Complete cptn:
! -cp*exp(-s)*(T(k,n-1)/(2*Dt) - V.del(T(k,n)) +Q/cp)
!
          cptn(k,i) = cptn(k,i)-expz(k)*total_heat(k,i) ! s1
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!     call addfld('EHEAT',' ',' ',eheating,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('CPTN2',' ',' ',cptn,
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('QM'     ,' ',' ',qm,
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('HEATING2',' ',' ',total_heat,
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! H = R*T/(M*g)                            (s4)
! rho = p0*exp(-s)*M/(R*T)                 (s5)
! tni = T                                  (s6)
!
! Levels 2 through lev1-1:
      do i=lon0,lon1
!
! G = g*(kT + H**2*rho*cp*kE)/(p0*H*Ds**2) (s2)
! F = g*(kE*H**3*rho*g/T)/(p0*2*H*Ds)      (s3)
!
        do k=lev0,lev1-1 
         g(k,i) = grav*(kt(k,i,lat)+h(k,i)**2*rho(k,i)*cp(k,i,lat)*
     |      dift(k,i,lat))/(p0*h(k,i)*dz**2)
         f(k,i)=grav**2*dift(k,i,lat)*h(k,i)**2*rho(k,i)/(tni(k,i)*
     |          p0*2.*dz)
        enddo ! k=lev0,lev1-1 
      enddo ! i=lon0,lon1
!
!     call addfld('G'  ,' ',' ',g  ,'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('F'  ,' ',' ',f  ,'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Coefficients for trsolv:
! Levels 3/2 through K-3/2
      do i=lon0,lon1
        do k=lev0,lev1-2
          p(k,i) = g(k,i)-f(k,i)
          q(k,i) = -g(k,i)-g(k+1,i) - f(k,i)+f(k+1,i)
          r(k,i) = g(k+1,i) + f(k+1,i) 
          rhs(k,i) = cptn(k,i)
        enddo ! k=lev0,lev1-2
! Level k-1/2
        p(lev1-1,i) =  g(lev1-1,i)-f(lev1-1,i)
        q(lev1-1,i) = -g(lev1-1,i)-f(lev1-1,i)
        r(lev1-1,i) = 0.
        rhs(lev1-1,i) = cptn(lev1-1,i)
      enddo ! i=lon0,lon1

!--------------------------------------------------------
! Set Diagnostic Variable cond (Ergs/g/s) 
! Set Diagnostic Variable cond2 (K/day) 
! Set Diagnostics for adiabatic lapse rate comparison  (K/km)
      do i = lon0,lon1
        do k = lev0+1,lev1-1
          tkp(k,i) = tni(k+1,i)
          tkm(k,i) = tni(k-1,i)
          tk(k,i) = tni(k,i)
          cond(k,i,lat) = (tkm(k,i)*p(k,i) + tk(k,i)*q(k,i) + 
     |                     tkp(k,i)*r(k,i))/expz(k)
          cond2(k,i,lat) = cond (k,i,lat)*86400./
     |                     (0.5*(cp(k,i,lat)+cp(k+1,i,lat)))
          dtdz(k,i,lat) = (tkp(k,i)-tk(k,i))/(h(k,i)*0.5)*1.0E+05
          alapse(k,i,lat) = -2.*grav/(cp(k,i,lat)+cp(k+1,i,lat))*1.0E+05
        enddo ! k=lev0+1,lev1-1
!
        dtdz(lev0,i,lat) = (tk(lev0+1,i)-tnlbc(i,lat))/(h(lev0,i)*0.5)*
     |    1.0E+05
      enddo ! i=lon0,lon1
!
      cond(lev0,:,lat) = cond(lev0+1,:,lat)
      cond(lev1,:,lat) = 0.0D0 
      alapse(lev0,:,lat) = alapse(lev0+1,:,lat)
      alapse(lev1,:,lat) = alapse(lev1-1,:,lat)
      dtdz(lev1,:,lat) = 0.0
!--------------------------------------------------------
!
!     call addfld('P_COEF0' ,' ',' ',p,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('Q_COEF0' ,' ',' ',q,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('R_COEF0' ,' ',' ',r,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('RHS0'    ,' ',' ',rhs,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
      do i=lon0,lon1
        do k=lev0,lev1-1
!---------------------------------------------
! Set Diagnostic Variable adia, qir (Ergs/g/s)
! Set Diagnostic Variable adia2, qir2 (K/day)
          adia(k,i,lat)=.5*(w_upd(k,i,lat)+w_upd(k+1,i,lat))*
     |                   (gask/mbar(k,i))*tni(k,i) 
          adia2(k,i,lat)=-(adia(k,i,lat)/cp(k,i,lat))*86400.
          qir(k,i,lat) = 0.5*cool_exp(k,i,lat)*
     |                  (cp(k,i,lat)+cp(k+1,i,lat))
          qir2(k,i,lat) = cool_exp(k,i,lat)*86400.
          qirez(k,i,lat) = qir(k,i,lat)*expz(k)
!---------------------------------------------
!    rhs units must be     :  erg/g/sec*expz
!    q   units must be     :  erg/g/sec*expz
!    qpart units must be   :  erg/g/sec
!    cool_exp units must be:  K/sec
!    cool_imp units must be:  K/sec
!---------------------------------------------
!   LHS Correction: 10/31/06 Bougher
!      qpart(k,i) = 
!    |  .5*(cp(k,i,lat)+cp(k+1,i,lat))*(dtx2inv+cool_imp(k,i,lat))+
!    |  .5*(w_upd(k,i,lat)+w_upd(k+1,i,lat))*gask/mbar(k,i) 
       qpart(k,i) = 
     |  .5*(cp(k,i,lat)+cp(k+1,i,lat))*(dtx2inv-cool_imp(k,i,lat))+
     |  .5*(w_upd(k,i,lat)+w_upd(k+1,i,lat))*gask/mbar(k,i) 
!
!   RHS Correction: 10/26/06 Bougher
!   RHS Correction: 10/31/06 Bougher
!         rhs(k,i) = rhs(k,i) + 0.5*expz(k)*(cool_exp(k,i,lat) + 
!    |      cool_exp(k+1,i,lat))*0.5*(cp(k,i,lat)+cp(k+1,i,lat))
          rhs(k,i) = rhs(k,i) - 0.5*expz(k)*cool_exp(k,i,lat)*
     |      (cp(k,i,lat)+cp(k+1,i,lat))

!         rhs(k,i) = rhs(k,i)+cool_exp(k,i,lat) ! as in tiegcm1.8
!         (expz factor already included; positive value returned)

          q(k,i) = q(k,i)-expz(k)*qpart(k,i)
!---------------------------------------------
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

!     call addfld('QPART'  ,' ',' ',qpart,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('Q_COEF1',' ',' ',q    ,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('RHS1'   ,' ',' ',rhs  ,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     do k=lev0,lev1
!       tnlbc_diag(k,:) = tnlbc(:,lat)
!     enddo
!     call addfld('TN_LBC2',' ',' ',tnlbc_diag,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Lower boundary:
      do i=lon0,lon1
        q(lev0,i) = q(lev0,i)-p(lev0,i)
!
! Diffs in rhs lbc ??:
        rhs(lev0,i) = rhs(lev0,i)-2.*p(lev0,i)*tnlbc(i,lat)
        p(lev0,i) = 0.
      enddo ! i=lon0,lon1

!     call addfld('P_COEF2',' ',' ',p,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('Q_COEF2',' ',' ',q,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('R_COEF2',' ',' ',r,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('RHS2'   ,' ',' ',rhs,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)

!
! Solve tridiagonal system for new tn:
!
      call trsolv(p,q,r,rhs,tn_upd(:,lon0:lon1,lat),lev0,lev1,
     |  lev0,lev1-1,lon0,lon1,nlonp4,lat,0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!**** added statment to skip update
!     tn_upd(:,lon0:lon1,lat) = tn_nm(:,lon0:lon1,lat)
      
!     call addfld('TN_SOLV',' ',' ',tn_upd(:,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! DT/DZ below z = -13.0 must be -3.8K/km (ala VIRA). Stability at LBC.
! Specify Lapse rate below z = -13.0 (about 83 km), xlapse = 3.8 K/km
! DT/DZ below z = -14.0 must be -5.0K/km (ala VIRA). Stability at LBC.
! Specify Lapse rate below z = -14.0 (about 80 km). xlapse = 5.0 K/km
!
! Do this only if not nudging with oxvgcm:
      if (len_trim(oxvgcm_ncfile) == 0) then
        do i=lon0,lon1
          if (dz == 0.5) then
            do k=lev0,6
                kk = 7-k
                dtdzp = (xlapse*h(kk+1,i)*0.5/1.0E+05)
                tn_upd(kk,i,lat) = tn_upd(kk+1,i,lat)+dtdzp
            enddo
          elseif (dz == 0.25) then
            do k=lev0,12
                kk = 13-k
                dtdzp = (xlapse*h(kk+1,i)*0.5/1.0E+05)
                tn_upd(kk,i,lat) = tn_upd(kk+1,i,lat)+dtdzp
            enddo
          endif ! dz = 0.5 or 0.25
          tn_upd(lev0,i,lat)=(tn_upd(lev0+1,i,lat)+dtdzp+tnlbc(i,lat))
     |      /2.
        enddo
      endif ! if not nudging w/ oxvgcm
!
! End second latitude scan:
      enddo ! lat=lat0,lat1
!
! -----------------------------------------------------
! Begin 1st Diagnostic Latitude Scan(s):---------------
      do lat=lat0,lat1
        do i=lon0,lon1
           do k=lev0,lev1-1
!     conversions (to K/day from K/s) and sign change
        xhydrot2(k,i,lat) = -xhydrot(k,i,lat)*86400.
! Calculating the total dynamics adia2=(K/day), xhydrot2=(K/Day)
        totdynm(k,i,lat) = adia2(k,i,lat)+ xhydrot2(k,i,lat)
! Calculating qnet2 from all heating and cooling terms (K/Day)
!     qnet2(k,i,lat) = qeuv2(k,i,lat)+ qnir2(k,i,lat)  +
!    | qir2(k,i,lat) + cond2(k,i,lat) + qaero2(k,i,lat) + 
!    |   totdynm(k,i,lat)
! Calculating qnet2 from only radiative heating and cooling terms (K/Day)
      qnet2(k,i,lat) = qeuv2(k,i,lat)+ qnir2(k,i,lat)  +
     | qir2(k,i,lat) +  qaero2(k,i,lat)
!
           enddo ! k  = lev0,lev1-1
        enddo ! i  = lon0,lon1
!
! End 1st Diagnostic latitude scan:
      enddo ! lat  = lat0,lat1
! -----------------------------------------------------
!
! Begin 2nd Diagnostic Latitude Scan(s):---------------
!
      do lat=lat0,lat1
!
! Breakout total dynamics terms (K/s)
!       call addfld('XhydroH',' ',' ',xhydroh(:,:,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('XhydroV',' ',' ',xhydrov(:,:,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('XhydroT',' ',' ',xhydrot(:,:,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! ----------------------------------------------------
!  ***  Conversion (erg/g/sec to erg/cm3/sec) OFF
!       qeuv(:,:,lat) = qeuv(:,:,lat)*mden(:,:,lat)
!       qeuvez(:,:,lat) = qeuvez(:,:,lat)*mden(:,:,lat)
!       qir(:,:,lat) = qir(:,:,lat)*mden(:,:,lat)
!       qirez(:,:,lat) = qirez(:,:,lat)*mden(:,:,lat)
!       qmol(:,:,lat) = qmol(:,:,lat)*mden(:,:,lat)
!       qel(:,:,lat) = qel(:,:,lat)*mden(:,:,lat)
!       cond(:,:,lat) = cond(:,:,lat)*mden(:,:,lat)
!       adia(:,:,lat) = -1.*mden(:,:,lat)*adia(:,:,lat)
! ----------------------------------------------------
        qt(:,:,lat) = qtot(:,:,lat) + qir(:,:,lat) 
        adia(:,:,lat) = -1.*adia(:,:,lat)
! ----------------------------------------------------
! ***   Radiative Heating and Cooling (erg/g/sec)
!       call addfld('QTOT',' ',' ',qtot(:,:,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('QEUV',' ',' ',qeuv(:,:,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('QIR',' ',' ',qir(:,:,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('QNET',' ',' ',qt(:,:,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('QAERO',' ',' ',qaero(:,:,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
! ----------------------------------------------------
! ***   Radiative Heating and Cooling (K/day)
        call addfld('QEUV2',' ',' ',qeuv2(:,:,lat),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
        call addfld('QNIR2',' ',' ',qnir2(:,:,lat),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
        call addfld('QIR2',' ',' ',qir2(:,:,lat),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
        call addfld('COND2',' ',' ',cond2(:,:,lat),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
        call addfld('TOTDYNM',' ',' ',totdynm(:,:,lat),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
        call addfld('QAERO2',' ',' ',qaero2(:,:,lat),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
        call addfld('QNET2',' ',' ',qnet2(:,:,lat),
     |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
! ----------------------------------------------------
! ***   DTDZ and ALAPSE (K/k unitsm)
!       call addfld('ALAPSE',' ',' ',alapse(:,:,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('DTDZ',' ',' ',dtdz(:,:,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
! ----------------------------------------------------
! ***   Breakout Dynamical Heating and Cooling (K/day)              
!       call addfld('ADIA',' ',' ',adia(:,:,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('XHYDROT',' ',' ',xhydrot(:,:,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('TOTDYNM',' ',' ',totdynm(:,:,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
! ----------------------------------------------------
! ***   O+O+CO2 RECOMBINATION Heating (erg/g/sec)
!       call addfld('QREC',' ',' ',qrec(:,:,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
! ----------------------------------------------------
! ***   Other Products       
!       call addfld('QeuvEZ',' ',' ',qeuvez(:,:,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('QirEZ',' ',' ',qirez(:,:,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('Qmol',' ',' ',qmol(:,:,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('Qel',' ',' ',qel(:,:,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! end 2nd Diagnostic latitude scan:
      enddo ! lat  = lat0,lat1
!
! End All Diagnostics ---------------------------------
! -----------------------------------------------------
! Set kut for wave filtering according to dlat (2.5 or 5.0):
!     call set_wave_filter(36,kut_5,nlat,kutt)
!
! Filter updated tn: turn on for polar region stability
      call filter_tn(tn_upd,lev0,lev1,lon0,lon1,lat0,lat1,kutt)
!     call filter_tn(tn_upd,lev0,lev1,lon0,lon1,lat0,lat1,kut)
!     do lat=lat0,lat1
!       call addfld('TN_FILT',' ',' ',tn_upd(:,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     enddo
!
! Third latitude scan:
      do lat=lat0,lat1
!
! Smooth updated tn:
        do i=lon0,lon1 
          do k=lev0,lev1-1
            tn_nm_upd(k,i,lat) = dtsmooth_div2*(tn_nm(k,i,lat)+
     |        tn_upd(k,i,lat)) + dtsmooth*tn(k,i,lat)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1 
        tn_nm_upd(lev1,lon0:lon1,lat) = spval

!       call addfld('TN_NMOUT',' ',' ',tn_nm_upd(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! 9/19/06 btf: put spval in top slot:
        tn_upd(lev1,lon0:lon1,lat) = spval
!
! save LB for tlb c= tn(itp) in tlbc_nm
! Save lower boundary for tn(itc) in tlbc for history:
! Lower boundary is in tlbc (fields.F):
        tlbc_nm(lon0:lon1,lat) = tlbc(lon0:lon1,lat)
        tlbc(lon0:lon1,lat)    = tnlbc(:,lat)
!
#ifdef MPI
!
! Define halo longitudes in tn_upd, for use by np4
        call mp_bndlons_f3d(tn_upd,nlevs,lon0,lon1,lat0-2,lat1+2,1)
!
! Periodic points:
        call mp_periodic_f3d(tn_upd(:,lon0:lon1,lat0-1:lat1+1),
     |    lev0,lev1,lon0,lon1,lat0-1,lat1+1)
#endif
!
! Overall:
! Tn must be at least 100 deg (Everywhere):
        lonbeg = lon0-2
        if (lon0==1) lonbeg = 1
        lonend = lon1+2
        if (lon1==nlonp4) lonend = nlonp4
        do i=lonbeg,lonend
          do k=lev0,lev1
            if (tn_upd(k,i,lat) < 100.) tn_upd(k,i,lat) = 100.
          enddo
        enddo
!
        call addfld('TN_FINAL',' ',' ',tn_upd(:,lon0:lon1,lat),
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
        call addfld('TNSEC',' ',' ',tn_upd(:,lon0:lon1,lat),
     |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! End fourth lat scan:
      enddo ! lat=lat0,lat1
      end subroutine dt
!-----------------------------------------------------------------------
      subroutine filter_tn(fout,lev0,lev1,lon0,lon1,lat0,lat1,kut)
!
! Filter updated W omega:
!
      use filter_module,only: filter
#ifdef MPI
      use mpi_module,only: mp_gatherlons_f3d,mp_scatterlons_f3d,mytidi
      implicit none
#else
      implicit none
      integer :: mytidi=0
#endif
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1,kut(nlat)
      real,intent(inout) :: fout(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2)
!
! VT vampir tracing:
!
#ifdef VT
#include "VT.inc"
#endif
!
! Local:
      integer :: i,j,nlevs,nlons,nlats
      real :: fik(nlonp4,lev0:lev1),fkij(lev0:lev1,nlonp4,lat0:lat1)
      real :: fmin,fmax
!
#ifdef VT
!     code = 131 ; state = 'filter_tn' ; activity='Filtering'
      call vtbegin(131,ier)
#endif
!
      nlevs = lev1-lev0+1
      nlons = lon1-lon0+1
      nlats = lat1-lat0+1
!
! Define lons in w_ki from current task:
      fkij = 0.
      do j=lat0,lat1
        do i=lon0,lon1
          fkij(:,i,j) = fout(:,i,j)
        enddo
      enddo ! j=lat0,lat1
!
#ifdef MPI
!
! Gather longitudes into tasks in first longitude column of task table
!   (leftmost of each j-row) for global fft. (i.e., tasks with mytidi==0 
!   gather lons from other tasks in that row). This includes all latitudes.
!
      call mp_gatherlons_f3d(fkij,lev0,lev1,lon0,lon1,lat0,lat1,1)
#endif
!
! Only leftmost tasks at each j-row of tasks does the global filtering:
      if (mytidi==0) then
!
! Define 2d array with all longitudes for filter at each latitude:
        latscan: do j=lat0,lat1
          if (kut(j) >= nlon/2) cycle latscan
          do i=1,nlonp4
            fik(i,:) = fkij(:,i,j)
          enddo ! i=1,nlonp4
!
! Remove wave numbers > kut(lat):
          call filter(fik,lev0,lev1,kut(j),j)
!
! Return filtered array to fkij:
          do i=1,nlonp4
            fkij(:,i,j) = fik(i,:)
          enddo ! i=1,nlonp4
        enddo latscan ! j=lat0,lat1
      endif ! mytidi==0
#ifdef MPI
!
! Now leftmost task at each j-row must redistribute filtered data
! back to other tasks in the j-row (mytidi>0,mytidj) (includes latitude):
!
      call mp_scatterlons_f3d(fkij,lev0,lev1,lon0,lon1,lat0,lat1,1)
#endif
!
! Return filtered array to fout at current task longitudes and latitudes:
      do j=lat0,lat1
        do i=lon0,lon1
          fout(:,i,j) = fkij(:,i,j)
        enddo
      enddo
!
#ifdef VT
!     code = 131 ; state = 'filter_tn' ; activity='Filtering'
      call vtend(131,ier)
#endif
      end subroutine filter_tn
!-----------------------------------------------------------------------
      subroutine cool(tn,o1,barm,cool_exp,cool_imp,lev0,lev1,lon0,lon1,
     |  lat0,lat1)
!
! 9/20/06 btf:
! Calculate explicit cooling. This is adapted from start.f and dtuv.f
! of old vtgcm code in hao:/home/tgcm/vtgcm/src.1996.
!
      use cons_module,only: rmassinv_o1,rtd
      use addfld_module,only: addfld
!
! chi(nlonp4,nlat) is solar zenith angle, dimensioned globally, but 
! defined only at subdomains (see chapman.F):
      use chapman_module,only: chi 
!
! Arrays read from DATAVOL (see read_ncdata.F): off
!     use read_ncdata_module,only: cool1,cool2,tm1a,tm2a,coolk,
!    |  aday,anight,xnco2d,xnco2n,xnod,xnon
! Arrays read from coolref_init_module 
      use coolref_init_module,only: cool1a,cool2a,tm1a,tm2a,coolk,
     |  aday,anight,xnco2d,xnco2n,xnod,xnon
!
#ifdef MPI
      use mpi_module,only: mp_gatherlons_f3d,mp_scatterlons_f3d,mytidi,
     |                     mp_periodic_f3d
      implicit none
#else
      implicit none
      integer :: mytidi=0
#endif
!
#include<defs.h>
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
      real,dimension(lev0:lev1,lon0:lon1,lat0:lat1),intent(in) ::
     |  tn,       ! tn (K)
     |  o1,       ! O (mmr)
     |  barm      ! mbar (mean molecular wt)
      real,dimension(lev0:lev1,lon0:lon1,lat0:lat1),intent(out) :: 
     |  cool_exp, cool_imp  ! output
!
! Local:
      integer :: i,ii,k,lonbeg,lonend,lat,iday,inight
      real,dimension(lev0:lev1) :: 
     |  szp,cooldy,coolnt,fac,facdres
      real,dimension(lev0:lev1,nlonp4) ::   ! with global longitudes
     |  xno,              ! O volume mixing ratio 
     |  xnoday,           ! O volume mixing ratio with 1.3 correction factor 
     |  xnonite,         ! O volume mixing ratio w/o 1.3 correction factor
     |  tnki,o1ki,barmki, ! t,o,barm
     |  coolexp,coolimp   ! output
      real :: xint,chiday,chinight
      real :: chideg(nlonp4)
!
! External:
      integer,external :: nearest ! util.F
!
! For gather/scatter:
      real,dimension(lev0:lev1,nlonp4,lat0:lat1,4) ::
     |  fkij ! t,o,barm,chi respectively, with global longitudes
!
!     write(6,"('Enter sub cool.')")
      do k=lev0,lev1
        szp(k) = zibot + .5*dz + (k-1)*dz
      enddo
!
! Load subdomains for longitude gather:
      do lat=lat0,lat1
        do i=lon0,lon1
          fkij(:,i,lat,1) = tn(:,i,lat)
          fkij(:,i,lat,2) = o1(:,i,lat)
          fkij(:,i,lat,3) = barm(:,i,lat)
          fkij(:,i,lat,4) = chi(i,lat)    ! redundant in k
        enddo
      enddo
!
! Gather all longitudes to leftmost tasks (mytidi==0):
#ifdef MPI
      call mp_gatherlons_f3d(fkij,lev0,lev1,lon0,lon1,lat0,lat1,4)  
#endif
!
! Whole-array inits:
      cool_exp = 0.
      cool_imp = 0.
!
! Now "leftmost" tasks have t,o,barm and chideg in global longitude:
! 
      if (mytidi==0) then ! I am a "leftmost" task
        coolexp = 0.
        coolimp = 0.
!
! Latitude scan. 
        do lat=lat0,lat1
!
! Solar zenith angle in degrees (recall fkij(:,i,lat,4) is redundant 
! in vertical:
          do i=1,nlonp4
            chideg(i) = fkij(1,i,lat,4)*rtd
          enddo
! -------------------------------------------------------------------
          chiday = 95.
          chinight = 95.
! -------------------------------------------------------------------
!         iday   = nearest(chideg,nlonp4,chiday)
!         inight = nearest(chideg,nlonp4,chinight)
!
!         write(6,"('cool: lat=',i3,' chiday  =',f8.3,' iday  =',
!    |      i3)") lat,chiday,iday
!         write(6,"('cool: lat=',i3,' chinight=',f8.3,' inight=',
!    |      i3)") lat,chinight,inight
!
C     N.L. COOLING COMPLETELY IN EXPLICIT FORM OF 2-D MODEL
!
! Define local (k,i) fields (global longitude):
!** xnoday has correction factor for specific alt.
!    (30% enhancement over ~100-120 km)
!** xnonite is without correction factor
! -- xno is atomic O volume mixing ratio


          fac = (/
     |      1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,
     |      1.,1.,1.,1.,1.,1.17,1.3,1.49,1.63,1.78,1.91,2.03,
     |      2.13,2.16,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,
     |      1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,
     |      1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1./)
          facdres = (/
     |      1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,
     |      1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,
     |      1.,1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3,
     |      1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3,
     |      1.3,1.3,1.3,1.3,1.3,1.3,1.,1.,1.,1.,1.,1.,1.,
     |      1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,
     |      1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,
     |      1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,
     |      1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,
     |      1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1./)

          do i=1,nlonp4
            tnki(:,i) = fkij(:,i,lat,1)
            o1ki(:,i) = fkij(:,i,lat,2)
            barmki(:,i) = fkij(:,i,lat,3)
            do k=lev0,lev1-1
              xno(k,i) = o1ki(k,i)*0.5*(barmki(k,i)+barmki(k+1,i))*
     |          rmassinv_o1
              xnonite(k,i) = xno(k,i)
#if (NLEV==68)
! ** ZP = -8.5 to -2
                xnoday(k,i) = xno(k,i)*fac(k)
#elif(NLEV==136)
                xnoday(k,i) = xno(k,i)*facdres(k)
#else
 write(6,"('ERROR:unsupported NLEV=',i6)")nlev
#endif
            enddo
          enddo
!
! Dayside cooling
!
          do i=1,nlonp4
            if (chideg(i) <= chiday) then
              call expcool(tnki(:,i),tm1a,cool1a,xnco2d,xnod,aday,
!             call expcool(tnki(:,i),tm1a,cool1,xnco2d,xnod,aday,
     |          xno(:,i),coolexp(:,i),lev0,lev1)
!             coolimp(:,i) = aday(:)
              coolimp(:,i) = 0.
!             if (i==iday) then
!               cooldy(:) = coolexp(:,i)
!               write(6,"('cool: i=',i3,' lat=',i3,' cooldy min,max=',
!    |            2e12.4)") i,lat,minval(cooldy),maxval(cooldy)
!             endif
!             write(6,"('dayside cooling: lat=',i3,' i=',i3,
!    |          ' coolexp min,max=',2e12.4)") lat,i,
!    |          minval(coolexp(:,i)),maxval(coolexp(:,i))

            endif
          enddo ! dayside i
!
! Nightside cooling
!
          do i=1,nlonp4
!           if (chideg(i) >= chinight) then
            if (chideg(i) > chinight) then
              call expcool(tnki(:,i),tm2a,cool2a,xnco2n,xnon,anight,
!             call expcool(tnki(:,i),tm2a,cool2,xnco2n,xnon,anight,
     |          xno(:,i),coolexp(:,i),lev0,lev1)
!             coolimp(:,i) = anight(:)
              coolimp(:,i) = 0.
!             if (i==inight) then
!               coolnt(:) = coolexp(:,i)
!               write(6,"('cool: i=',i3,' lat=',i3,' coolnt min,max=',
!    |            2e12.4)") i,lat,minval(coolnt),maxval(coolnt)
!             endif
!             write(6,"('nightside cooling: lat=',i3,' i=',i3,
!    |          ' coolexp min,max=',2e12.4)") lat,i,
!    |          minval(coolexp(:,i)),maxval(coolexp(:,i))
            endif
          enddo ! nightside i
!
! Upper boundary:
          coolexp(lev1,:) = coolexp(lev1-1,:)
          coolimp(lev1,:) = coolimp(lev1-1,:)
!
! Load 1st 2 slots of fkij with cooling, for scatter to mytid /= 0:
          do i=1,nlonp4
            fkij(:,i,lat,1) = coolexp(:,i)
            fkij(:,i,lat,2) = coolimp(:,i)
          enddo
!
!         write(6,"('cool: lat=',i3,' coolexp min,max=',2e12.4)")
!    |      lat,minval(coolexp(:,:)),maxval(coolexp(:,:))

        enddo ! lat=lat0,lat1
      endif ! mytidi==0
!
! Scatter output back to subdomains:
#ifdef MPI
      call mp_scatterlons_f3d(fkij,lev0,lev1,lon0,lon1,lat0,lat1,2)
#endif
!
! Define output arrays at subdomains:
      do lat=lat0,lat1
        do i=lon0,lon1
          cool_exp(:,i,lat) = fkij(:,i,lat,1)
          cool_imp(:,i,lat) = fkij(:,i,lat,2)
        enddo
!       write(6,"('cool: lat=',i3,' cool_exp min,max=',2e12.4)")
!    |    lat,minval(cool_exp(:,:,lat)),maxval(cool_exp(:,:,lat))
#ifdef MPI
!
! Periodic points:
        call mp_periodic_f3d(cool_exp(:,:,lat),
     |    lev0,lev1,lon0,lon1,lat,lat)
        call mp_periodic_f3d(cool_imp(:,:,lat),
     |    lev0,lev1,lon0,lon1,lat,lat)
#else
!     real,dimension(lev0:lev1,lon0:lon1,lat0:lat1) ::
!    |  cool_exp,   ! explicit cooling (see sub cool) K/sec

        cool_exp(:,1:2,:) = cool_exp(:,nlon+1:nlonp2,:)
        cool_exp(:,nlon+3:nlonp4,:) = cool_exp(:,3:4,:)
        cool_imp(:,1:2,:) = cool_imp(:,nlon+1:nlonp2,:)
        cool_imp(:,nlon+3:nlonp4,:) = cool_imp(:,3:4,:)
#endif
! There are similiar calls to addfld after the dt call to cool:
!       call addfld('COOL_EXP',' ',' ',cool_exp(:,:,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('COOL_IMP',' ',' ',cool_imp(:,:,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
      enddo ! lat0,lat1
      end subroutine cool
!-----------------------------------------------------------------------
      subroutine expcool(tn,tma,cooldat,xnco2r,xnor,adn,xno,coolexp,
     |  lev0,lev1)
!     use read_ncdata_module, only: cool1,cool2,coolk
! Arrays read from coolref_init_module 
      use coolref_init_module, only: coolk
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1
      real,dimension(lev0:lev1),intent(in) :: 
     |  tn,       ! tn column at current lon,lat
     |  tma,      ! either tm1a (day), or tm2a (night)
     |  cooldat,  ! either cool1a (day), or cool2a (nite)
     |  xnco2r,   ! either xnco2d (day), or xnco2n (nite): VMR
     |  xnor,     ! either xnod (day), or xnon (nite); VMR
     |  adn,      ! either aday or anight
     |  xno       ! O volume mixing ratio
      real,dimension(lev0:lev1),intent(out) :: 
     |  coolexp   ! explicit cooling column output
!
! Local:
      integer :: k
      real :: roco2,roco2m,xxx,efd,efn,gee,szp
      real,parameter :: onethird=1./3.
!
      do k=lev0,lev1-1
        szp = zibot + .5*dz + (k-1)*dz
        roco2 = coolk*sqrt(tn(k))/(5.2e-15*(1.+2.e-3*
     |    (tn(k)-220.))*tn(k)/273.3)
        roco2m = coolk*sqrt(tma(k))/(5.2e-15*(1.+2.e-3*
     |    (tma(k)-220.))*(tma(k)/273.3))
        xxx = exp(40.6*(1./tma(k)**onethird-1./tn(k)**onethird))
!       efd = exp(szp+1.5)*((xnco2d(k)+xnod(k))/(xnco2d(k)+xnod(k)*
!     |    roco2m))
!       efn = exp(szp+1.5)*((xnco2d(k)+xno(k))/(xnco2d(k)+xno(k)*
!    |    roco2))
        efd = exp(szp+1.5)*((xnco2r(k)+xnor(k))/(xnco2r(k)+xnor(k)*
     |    roco2m))
        efn = exp(szp+1.5)*((xnco2r(k)+xno(k))/(xnco2r(k)+xno(k)*
     |    roco2))
        gee = (1.+efd)/(1.+efn*xxx)
!  ------------------------------------------------------------------
!  Test both cool_exp and cool_imp on.
!       coolexp(k) = cooldat(k)*(exp(960.*(1./tma(k)-1./tn(k)))*
!    |    gee)-adn(k)*tn(k)
!  Test all cool_exp on ; cool_imp off.
        coolexp(k) = cooldat(k)*(exp(960.*(1./tma(k)-1./tn(k)))*
     |    gee)
!  ------------------------------------------------------------------
      enddo ! k=lev0,lev1
      end subroutine expcool
!-----------------------------------------------------------------------
      subroutine calc_tzm(tn,tzm,lev0,lev1,lon0,lon1,lat0,lat1)
      use params_module,only: nlonp4
#ifdef MPI
      use mpi_module,only: mp_gatherlons_f3d,mp_scatterlons_f3d,mytidi
      implicit none
#else
      implicit none
      integer :: mytidi=0
#endif
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
      real,intent(in) :: tn(lev0:lev1,lon0:lon1,lat0:lat1)
      real,intent(out) :: tzm(lev0:lev1,lat0:lat1)
!
! Local:
      integer :: k,i,j
      real :: f_kij(lev0:lev1,nlonp4,lat0:lat1),rlon

      f_kij = 0.
      do j=lat0,lat1
        do i=lon0,lon1
          f_kij(:,i,j) = tn(:,i,j)
        enddo
      enddo
#ifdef MPI
!
! Gather global longitudes to western-most task in lat row:
      call mp_gatherlons_f3d(f_kij,lev0,lev1,lon0,lon1,lat0,lat1,1)
#endif
      tzm = 0.
      rlon = 1./real(nlonp4-4)
!
! Calculate zonal mean:
      if (mytidi==0) then
        do j=lat0,lat1
          do i=3,nlonp4-2
            tzm(:,j) = tzm(:,j)+f_kij(:,i,j)
          enddo
          tzm(:,j) = tzm(:,j)*rlon
        enddo
!
! Copy tzm redundantly to i-index of f_kij:
        do j=lat0,lat1
          do i=1,nlonp4
            f_kij(:,i,j) = tzm(:,j)
          enddo
        enddo
      endif ! mytidi==0
#ifdef MPI
!
! Scatter 3d w/ zonal means to other tasks in my lat row:
      call mp_scatterlons_f3d(f_kij,lev0,lev1,lon0,lon1,lat0,lat1,1)
#endif
!
! Finally copy zm to my subdomain tzm:
      do j=lat0,lat1
        tzm(:,j) = f_kij(:,lon0,j) ! all i's in f_kij are zm, so just use lon0
      enddo
      end subroutine calc_tzm
!-----------------------------------------------------------------------
      end module dt_module
